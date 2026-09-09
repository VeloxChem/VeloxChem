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


#include "SimdTransferHF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_hf_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t hd, const size_t id,
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

        const auto *hd_0 = buffer.data(hd + 0 * ncomps + c);
        const auto *hd_1 = buffer.data(hd + 1 * ncomps + c);
        const auto *hd_2 = buffer.data(hd + 2 * ncomps + c);
        const auto *hd_3 = buffer.data(hd + 3 * ncomps + c);
        const auto *hd_4 = buffer.data(hd + 4 * ncomps + c);
        const auto *hd_5 = buffer.data(hd + 5 * ncomps + c);
        const auto *hd_6 = buffer.data(hd + 6 * ncomps + c);
        const auto *hd_7 = buffer.data(hd + 7 * ncomps + c);
        const auto *hd_8 = buffer.data(hd + 8 * ncomps + c);
        const auto *hd_9 = buffer.data(hd + 9 * ncomps + c);
        const auto *hd_10 = buffer.data(hd + 10 * ncomps + c);
        const auto *hd_11 = buffer.data(hd + 11 * ncomps + c);
        const auto *hd_12 = buffer.data(hd + 12 * ncomps + c);
        const auto *hd_13 = buffer.data(hd + 13 * ncomps + c);
        const auto *hd_14 = buffer.data(hd + 14 * ncomps + c);
        const auto *hd_15 = buffer.data(hd + 15 * ncomps + c);
        const auto *hd_16 = buffer.data(hd + 16 * ncomps + c);
        const auto *hd_17 = buffer.data(hd + 17 * ncomps + c);
        const auto *hd_18 = buffer.data(hd + 18 * ncomps + c);
        const auto *hd_19 = buffer.data(hd + 19 * ncomps + c);
        const auto *hd_20 = buffer.data(hd + 20 * ncomps + c);
        const auto *hd_21 = buffer.data(hd + 21 * ncomps + c);
        const auto *hd_22 = buffer.data(hd + 22 * ncomps + c);
        const auto *hd_23 = buffer.data(hd + 23 * ncomps + c);
        const auto *hd_24 = buffer.data(hd + 24 * ncomps + c);
        const auto *hd_25 = buffer.data(hd + 25 * ncomps + c);
        const auto *hd_26 = buffer.data(hd + 26 * ncomps + c);
        const auto *hd_27 = buffer.data(hd + 27 * ncomps + c);
        const auto *hd_28 = buffer.data(hd + 28 * ncomps + c);
        const auto *hd_29 = buffer.data(hd + 29 * ncomps + c);
        const auto *hd_30 = buffer.data(hd + 30 * ncomps + c);
        const auto *hd_31 = buffer.data(hd + 31 * ncomps + c);
        const auto *hd_32 = buffer.data(hd + 32 * ncomps + c);
        const auto *hd_33 = buffer.data(hd + 33 * ncomps + c);
        const auto *hd_34 = buffer.data(hd + 34 * ncomps + c);
        const auto *hd_35 = buffer.data(hd + 35 * ncomps + c);
        const auto *hd_36 = buffer.data(hd + 36 * ncomps + c);
        const auto *hd_37 = buffer.data(hd + 37 * ncomps + c);
        const auto *hd_38 = buffer.data(hd + 38 * ncomps + c);
        const auto *hd_39 = buffer.data(hd + 39 * ncomps + c);
        const auto *hd_40 = buffer.data(hd + 40 * ncomps + c);
        const auto *hd_41 = buffer.data(hd + 41 * ncomps + c);
        const auto *hd_42 = buffer.data(hd + 42 * ncomps + c);
        const auto *hd_43 = buffer.data(hd + 43 * ncomps + c);
        const auto *hd_44 = buffer.data(hd + 44 * ncomps + c);
        const auto *hd_45 = buffer.data(hd + 45 * ncomps + c);
        const auto *hd_46 = buffer.data(hd + 46 * ncomps + c);
        const auto *hd_47 = buffer.data(hd + 47 * ncomps + c);
        const auto *hd_48 = buffer.data(hd + 48 * ncomps + c);
        const auto *hd_49 = buffer.data(hd + 49 * ncomps + c);
        const auto *hd_50 = buffer.data(hd + 50 * ncomps + c);
        const auto *hd_51 = buffer.data(hd + 51 * ncomps + c);
        const auto *hd_52 = buffer.data(hd + 52 * ncomps + c);
        const auto *hd_53 = buffer.data(hd + 53 * ncomps + c);
        const auto *hd_54 = buffer.data(hd + 54 * ncomps + c);
        const auto *hd_55 = buffer.data(hd + 55 * ncomps + c);
        const auto *hd_56 = buffer.data(hd + 56 * ncomps + c);
        const auto *hd_57 = buffer.data(hd + 57 * ncomps + c);
        const auto *hd_58 = buffer.data(hd + 58 * ncomps + c);
        const auto *hd_59 = buffer.data(hd + 59 * ncomps + c);
        const auto *hd_60 = buffer.data(hd + 60 * ncomps + c);
        const auto *hd_61 = buffer.data(hd + 61 * ncomps + c);
        const auto *hd_62 = buffer.data(hd + 62 * ncomps + c);
        const auto *hd_63 = buffer.data(hd + 63 * ncomps + c);
        const auto *hd_64 = buffer.data(hd + 64 * ncomps + c);
        const auto *hd_65 = buffer.data(hd + 65 * ncomps + c);
        const auto *hd_66 = buffer.data(hd + 66 * ncomps + c);
        const auto *hd_67 = buffer.data(hd + 67 * ncomps + c);
        const auto *hd_68 = buffer.data(hd + 68 * ncomps + c);
        const auto *hd_69 = buffer.data(hd + 69 * ncomps + c);
        const auto *hd_70 = buffer.data(hd + 70 * ncomps + c);
        const auto *hd_71 = buffer.data(hd + 71 * ncomps + c);
        const auto *hd_72 = buffer.data(hd + 72 * ncomps + c);
        const auto *hd_73 = buffer.data(hd + 73 * ncomps + c);
        const auto *hd_74 = buffer.data(hd + 74 * ncomps + c);
        const auto *hd_75 = buffer.data(hd + 75 * ncomps + c);
        const auto *hd_76 = buffer.data(hd + 76 * ncomps + c);
        const auto *hd_77 = buffer.data(hd + 77 * ncomps + c);
        const auto *hd_78 = buffer.data(hd + 78 * ncomps + c);
        const auto *hd_79 = buffer.data(hd + 79 * ncomps + c);
        const auto *hd_80 = buffer.data(hd + 80 * ncomps + c);
        const auto *hd_81 = buffer.data(hd + 81 * ncomps + c);
        const auto *hd_82 = buffer.data(hd + 82 * ncomps + c);
        const auto *hd_83 = buffer.data(hd + 83 * ncomps + c);
        const auto *hd_84 = buffer.data(hd + 84 * ncomps + c);
        const auto *hd_85 = buffer.data(hd + 85 * ncomps + c);
        const auto *hd_86 = buffer.data(hd + 86 * ncomps + c);
        const auto *hd_87 = buffer.data(hd + 87 * ncomps + c);
        const auto *hd_88 = buffer.data(hd + 88 * ncomps + c);

        const auto *id_0 = buffer.data(id + 0 * ncomps + c);
        const auto *id_1 = buffer.data(id + 1 * ncomps + c);
        const auto *id_2 = buffer.data(id + 2 * ncomps + c);
        const auto *id_3 = buffer.data(id + 3 * ncomps + c);
        const auto *id_4 = buffer.data(id + 4 * ncomps + c);
        const auto *id_5 = buffer.data(id + 5 * ncomps + c);
        const auto *id_6 = buffer.data(id + 6 * ncomps + c);
        const auto *id_7 = buffer.data(id + 7 * ncomps + c);
        const auto *id_8 = buffer.data(id + 8 * ncomps + c);
        const auto *id_9 = buffer.data(id + 9 * ncomps + c);
        const auto *id_10 = buffer.data(id + 10 * ncomps + c);
        const auto *id_11 = buffer.data(id + 11 * ncomps + c);
        const auto *id_12 = buffer.data(id + 12 * ncomps + c);
        const auto *id_13 = buffer.data(id + 13 * ncomps + c);
        const auto *id_14 = buffer.data(id + 14 * ncomps + c);
        const auto *id_15 = buffer.data(id + 15 * ncomps + c);
        const auto *id_16 = buffer.data(id + 16 * ncomps + c);
        const auto *id_17 = buffer.data(id + 17 * ncomps + c);
        const auto *id_18 = buffer.data(id + 18 * ncomps + c);
        const auto *id_19 = buffer.data(id + 19 * ncomps + c);
        const auto *id_20 = buffer.data(id + 20 * ncomps + c);
        const auto *id_21 = buffer.data(id + 21 * ncomps + c);
        const auto *id_22 = buffer.data(id + 22 * ncomps + c);
        const auto *id_23 = buffer.data(id + 23 * ncomps + c);
        const auto *id_24 = buffer.data(id + 24 * ncomps + c);
        const auto *id_25 = buffer.data(id + 25 * ncomps + c);
        const auto *id_26 = buffer.data(id + 26 * ncomps + c);
        const auto *id_27 = buffer.data(id + 27 * ncomps + c);
        const auto *id_28 = buffer.data(id + 28 * ncomps + c);
        const auto *id_29 = buffer.data(id + 29 * ncomps + c);
        const auto *id_30 = buffer.data(id + 30 * ncomps + c);
        const auto *id_31 = buffer.data(id + 31 * ncomps + c);
        const auto *id_32 = buffer.data(id + 32 * ncomps + c);
        const auto *id_33 = buffer.data(id + 33 * ncomps + c);
        const auto *id_34 = buffer.data(id + 34 * ncomps + c);
        const auto *id_35 = buffer.data(id + 35 * ncomps + c);
        const auto *id_36 = buffer.data(id + 36 * ncomps + c);
        const auto *id_37 = buffer.data(id + 37 * ncomps + c);
        const auto *id_38 = buffer.data(id + 38 * ncomps + c);
        const auto *id_39 = buffer.data(id + 39 * ncomps + c);
        const auto *id_40 = buffer.data(id + 40 * ncomps + c);
        const auto *id_41 = buffer.data(id + 41 * ncomps + c);
        const auto *id_42 = buffer.data(id + 42 * ncomps + c);
        const auto *id_43 = buffer.data(id + 43 * ncomps + c);
        const auto *id_44 = buffer.data(id + 44 * ncomps + c);
        const auto *id_45 = buffer.data(id + 45 * ncomps + c);
        const auto *id_46 = buffer.data(id + 46 * ncomps + c);
        const auto *id_47 = buffer.data(id + 47 * ncomps + c);
        const auto *id_48 = buffer.data(id + 48 * ncomps + c);
        const auto *id_49 = buffer.data(id + 49 * ncomps + c);
        const auto *id_50 = buffer.data(id + 50 * ncomps + c);
        const auto *id_51 = buffer.data(id + 51 * ncomps + c);
        const auto *id_52 = buffer.data(id + 52 * ncomps + c);
        const auto *id_53 = buffer.data(id + 53 * ncomps + c);
        const auto *id_54 = buffer.data(id + 54 * ncomps + c);
        const auto *id_55 = buffer.data(id + 55 * ncomps + c);
        const auto *id_56 = buffer.data(id + 56 * ncomps + c);
        const auto *id_57 = buffer.data(id + 57 * ncomps + c);
        const auto *id_58 = buffer.data(id + 58 * ncomps + c);
        const auto *id_59 = buffer.data(id + 59 * ncomps + c);
        const auto *id_60 = buffer.data(id + 60 * ncomps + c);
        const auto *id_61 = buffer.data(id + 61 * ncomps + c);
        const auto *id_62 = buffer.data(id + 62 * ncomps + c);
        const auto *id_63 = buffer.data(id + 63 * ncomps + c);
        const auto *id_64 = buffer.data(id + 64 * ncomps + c);
        const auto *id_65 = buffer.data(id + 65 * ncomps + c);
        const auto *id_66 = buffer.data(id + 66 * ncomps + c);
        const auto *id_67 = buffer.data(id + 67 * ncomps + c);
        const auto *id_68 = buffer.data(id + 68 * ncomps + c);
        const auto *id_69 = buffer.data(id + 69 * ncomps + c);
        const auto *id_70 = buffer.data(id + 70 * ncomps + c);
        const auto *id_71 = buffer.data(id + 71 * ncomps + c);
        const auto *id_72 = buffer.data(id + 72 * ncomps + c);
        const auto *id_73 = buffer.data(id + 73 * ncomps + c);
        const auto *id_74 = buffer.data(id + 74 * ncomps + c);
        const auto *id_75 = buffer.data(id + 75 * ncomps + c);
        const auto *id_76 = buffer.data(id + 76 * ncomps + c);
        const auto *id_77 = buffer.data(id + 77 * ncomps + c);
        const auto *id_78 = buffer.data(id + 78 * ncomps + c);
        const auto *id_79 = buffer.data(id + 79 * ncomps + c);
        const auto *id_80 = buffer.data(id + 80 * ncomps + c);
        const auto *id_81 = buffer.data(id + 81 * ncomps + c);
        const auto *id_82 = buffer.data(id + 82 * ncomps + c);
        const auto *id_83 = buffer.data(id + 83 * ncomps + c);
        const auto *id_84 = buffer.data(id + 84 * ncomps + c);
        const auto *id_85 = buffer.data(id + 85 * ncomps + c);
        const auto *id_86 = buffer.data(id + 86 * ncomps + c);
        const auto *id_87 = buffer.data(id + 87 * ncomps + c);
        const auto *id_88 = buffer.data(id + 88 * ncomps + c);
        const auto *id_89 = buffer.data(id + 89 * ncomps + c);
        const auto *id_93 = buffer.data(id + 93 * ncomps + c);
        const auto *id_94 = buffer.data(id + 94 * ncomps + c);
        const auto *id_95 = buffer.data(id + 95 * ncomps + c);
        const auto *id_99 = buffer.data(id + 99 * ncomps + c);
        const auto *id_100 = buffer.data(id + 100 * ncomps + c);
        const auto *id_101 = buffer.data(id + 101 * ncomps + c);
        const auto *id_105 = buffer.data(id + 105 * ncomps + c);
        const auto *id_106 = buffer.data(id + 106 * ncomps + c);
        const auto *id_107 = buffer.data(id + 107 * ncomps + c);
        const auto *id_111 = buffer.data(id + 111 * ncomps + c);
        const auto *id_112 = buffer.data(id + 112 * ncomps + c);
        const auto *id_113 = buffer.data(id + 113 * ncomps + c);
        const auto *id_119 = buffer.data(id + 119 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, hd_0, hd_1, hd_2, hd_3, hd_4, id_0, \
                         id_1, id_2, id_3, id_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * hd_0[k]
                     + id_0[k];

            t_1[k] = ab_x[k] * hd_1[k]
                     + id_1[k];

            t_2[k] = ab_x[k] * hd_2[k]
                     + id_2[k];

            t_3[k] = ab_x[k] * hd_3[k]
                     + id_3[k];

            t_4[k] = ab_x[k] * hd_4[k]
                     + id_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, hd_3, hd_4, hd_5, id_5, \
                         id_9, id_10, id_11, id_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * hd_5[k]
                     + id_5[k];

            t_6[k] = ab_y[k] * hd_3[k]
                     + id_9[k];

            t_7[k] = ab_y[k] * hd_4[k]
                     + id_10[k];

            t_8[k] = ab_y[k] * hd_5[k]
                     + id_11[k];

            t_9[k] = ab_z[k] * hd_5[k]
                     + id_17[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, hd_6, hd_7, hd_8, hd_9, hd_10, \
                         id_6, id_7, id_8, id_9, id_10 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * hd_6[k]
                      + id_6[k];

            t_11[k] = ab_x[k] * hd_7[k]
                      + id_7[k];

            t_12[k] = ab_x[k] * hd_8[k]
                      + id_8[k];

            t_13[k] = ab_x[k] * hd_9[k]
                      + id_9[k];

            t_14[k] = ab_x[k] * hd_10[k]
                      + id_10[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, hd_9, hd_10, hd_11, \
                         id_11, id_21, id_22, id_23, id_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * hd_11[k]
                      + id_11[k];

            t_16[k] = ab_y[k] * hd_9[k]
                      + id_21[k];

            t_17[k] = ab_y[k] * hd_10[k]
                      + id_22[k];

            t_18[k] = ab_y[k] * hd_11[k]
                      + id_23[k];

            t_19[k] = ab_z[k] * hd_11[k]
                      + id_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, hd_12, hd_13, hd_14, hd_15, \
                         hd_16, id_12, id_13, id_14, id_15, id_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * hd_12[k]
                      + id_12[k];

            t_21[k] = ab_x[k] * hd_13[k]
                      + id_13[k];

            t_22[k] = ab_x[k] * hd_14[k]
                      + id_14[k];

            t_23[k] = ab_x[k] * hd_15[k]
                      + id_15[k];

            t_24[k] = ab_x[k] * hd_16[k]
                      + id_16[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, hd_15, hd_16, hd_17, \
                         id_17, id_27, id_28, id_29, id_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_x[k] * hd_17[k]
                      + id_17[k];

            t_26[k] = ab_y[k] * hd_15[k]
                      + id_27[k];

            t_27[k] = ab_y[k] * hd_16[k]
                      + id_28[k];

            t_28[k] = ab_y[k] * hd_17[k]
                      + id_29[k];

            t_29[k] = ab_z[k] * hd_17[k]
                      + id_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, hd_18, hd_19, hd_20, hd_21, \
                         hd_22, id_18, id_19, id_20, id_21, id_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * hd_18[k]
                      + id_18[k];

            t_31[k] = ab_x[k] * hd_19[k]
                      + id_19[k];

            t_32[k] = ab_x[k] * hd_20[k]
                      + id_20[k];

            t_33[k] = ab_x[k] * hd_21[k]
                      + id_21[k];

            t_34[k] = ab_x[k] * hd_22[k]
                      + id_22[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ab_y, ab_z, hd_21, hd_22, hd_23, \
                         id_23, id_39, id_40, id_41, id_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * hd_23[k]
                      + id_23[k];

            t_36[k] = ab_y[k] * hd_21[k]
                      + id_39[k];

            t_37[k] = ab_y[k] * hd_22[k]
                      + id_40[k];

            t_38[k] = ab_y[k] * hd_23[k]
                      + id_41[k];

            t_39[k] = ab_z[k] * hd_23[k]
                      + id_47[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, hd_24, hd_25, hd_26, hd_27, \
                         hd_28, id_24, id_25, id_26, id_27, id_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_x[k] * hd_24[k]
                      + id_24[k];

            t_41[k] = ab_x[k] * hd_25[k]
                      + id_25[k];

            t_42[k] = ab_x[k] * hd_26[k]
                      + id_26[k];

            t_43[k] = ab_x[k] * hd_27[k]
                      + id_27[k];

            t_44[k] = ab_x[k] * hd_28[k]
                      + id_28[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, hd_27, hd_28, hd_29, \
                         id_29, id_45, id_46, id_47, id_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * hd_29[k]
                      + id_29[k];

            t_46[k] = ab_y[k] * hd_27[k]
                      + id_45[k];

            t_47[k] = ab_y[k] * hd_28[k]
                      + id_46[k];

            t_48[k] = ab_y[k] * hd_29[k]
                      + id_47[k];

            t_49[k] = ab_z[k] * hd_29[k]
                      + id_53[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, hd_30, hd_31, hd_32, hd_33, \
                         hd_34, id_30, id_31, id_32, id_33, id_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * hd_30[k]
                      + id_30[k];

            t_51[k] = ab_x[k] * hd_31[k]
                      + id_31[k];

            t_52[k] = ab_x[k] * hd_32[k]
                      + id_32[k];

            t_53[k] = ab_x[k] * hd_33[k]
                      + id_33[k];

            t_54[k] = ab_x[k] * hd_34[k]
                      + id_34[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, hd_33, hd_34, hd_35, \
                         id_35, id_51, id_52, id_53, id_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_x[k] * hd_35[k]
                      + id_35[k];

            t_56[k] = ab_y[k] * hd_33[k]
                      + id_51[k];

            t_57[k] = ab_y[k] * hd_34[k]
                      + id_52[k];

            t_58[k] = ab_y[k] * hd_35[k]
                      + id_53[k];

            t_59[k] = ab_z[k] * hd_35[k]
                      + id_59[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, hd_36, hd_37, hd_38, hd_39, \
                         hd_40, id_36, id_37, id_38, id_39, id_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * hd_36[k]
                      + id_36[k];

            t_61[k] = ab_x[k] * hd_37[k]
                      + id_37[k];

            t_62[k] = ab_x[k] * hd_38[k]
                      + id_38[k];

            t_63[k] = ab_x[k] * hd_39[k]
                      + id_39[k];

            t_64[k] = ab_x[k] * hd_40[k]
                      + id_40[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, ab_y, ab_z, hd_39, hd_40, hd_41, \
                         id_41, id_63, id_64, id_65, id_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * hd_41[k]
                      + id_41[k];

            t_66[k] = ab_y[k] * hd_39[k]
                      + id_63[k];

            t_67[k] = ab_y[k] * hd_40[k]
                      + id_64[k];

            t_68[k] = ab_y[k] * hd_41[k]
                      + id_65[k];

            t_69[k] = ab_z[k] * hd_41[k]
                      + id_71[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, hd_42, hd_43, hd_44, hd_45, \
                         hd_46, id_42, id_43, id_44, id_45, id_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_x[k] * hd_42[k]
                      + id_42[k];

            t_71[k] = ab_x[k] * hd_43[k]
                      + id_43[k];

            t_72[k] = ab_x[k] * hd_44[k]
                      + id_44[k];

            t_73[k] = ab_x[k] * hd_45[k]
                      + id_45[k];

            t_74[k] = ab_x[k] * hd_46[k]
                      + id_46[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, hd_45, hd_46, hd_47, \
                         id_47, id_69, id_70, id_71, id_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * hd_47[k]
                      + id_47[k];

            t_76[k] = ab_y[k] * hd_45[k]
                      + id_69[k];

            t_77[k] = ab_y[k] * hd_46[k]
                      + id_70[k];

            t_78[k] = ab_y[k] * hd_47[k]
                      + id_71[k];

            t_79[k] = ab_z[k] * hd_47[k]
                      + id_77[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, hd_48, hd_49, hd_50, hd_51, \
                         hd_52, id_48, id_49, id_50, id_51, id_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * hd_48[k]
                      + id_48[k];

            t_81[k] = ab_x[k] * hd_49[k]
                      + id_49[k];

            t_82[k] = ab_x[k] * hd_50[k]
                      + id_50[k];

            t_83[k] = ab_x[k] * hd_51[k]
                      + id_51[k];

            t_84[k] = ab_x[k] * hd_52[k]
                      + id_52[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, hd_51, hd_52, hd_53, \
                         id_53, id_75, id_76, id_77, id_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * hd_53[k]
                      + id_53[k];

            t_86[k] = ab_y[k] * hd_51[k]
                      + id_75[k];

            t_87[k] = ab_y[k] * hd_52[k]
                      + id_76[k];

            t_88[k] = ab_y[k] * hd_53[k]
                      + id_77[k];

            t_89[k] = ab_z[k] * hd_53[k]
                      + id_83[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, hd_54, hd_55, hd_56, hd_57, \
                         hd_58, id_54, id_55, id_56, id_57, id_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * hd_54[k]
                      + id_54[k];

            t_91[k] = ab_x[k] * hd_55[k]
                      + id_55[k];

            t_92[k] = ab_x[k] * hd_56[k]
                      + id_56[k];

            t_93[k] = ab_x[k] * hd_57[k]
                      + id_57[k];

            t_94[k] = ab_x[k] * hd_58[k]
                      + id_58[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, ab_y, ab_z, hd_57, hd_58, hd_59, \
                         id_59, id_81, id_82, id_83, id_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * hd_59[k]
                      + id_59[k];

            t_96[k] = ab_y[k] * hd_57[k]
                      + id_81[k];

            t_97[k] = ab_y[k] * hd_58[k]
                      + id_82[k];

            t_98[k] = ab_y[k] * hd_59[k]
                      + id_83[k];

            t_99[k] = ab_z[k] * hd_59[k]
                      + id_89[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, hd_60, hd_61, hd_62, hd_63, \
                         hd_64, id_60, id_61, id_62, id_63, id_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_x[k] * hd_60[k]
                       + id_60[k];

            t_101[k] = ab_x[k] * hd_61[k]
                       + id_61[k];

            t_102[k] = ab_x[k] * hd_62[k]
                       + id_62[k];

            t_103[k] = ab_x[k] * hd_63[k]
                       + id_63[k];

            t_104[k] = ab_x[k] * hd_64[k]
                       + id_64[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, hd_63, hd_64, \
                         hd_65, id_65, id_93, id_94, id_95, id_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * hd_65[k]
                       + id_65[k];

            t_106[k] = ab_y[k] * hd_63[k]
                       + id_93[k];

            t_107[k] = ab_y[k] * hd_64[k]
                       + id_94[k];

            t_108[k] = ab_y[k] * hd_65[k]
                       + id_95[k];

            t_109[k] = ab_z[k] * hd_65[k]
                       + id_101[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, hd_66, hd_67, hd_68, hd_69, \
                         hd_70, id_66, id_67, id_68, id_69, id_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_x[k] * hd_66[k]
                       + id_66[k];

            t_111[k] = ab_x[k] * hd_67[k]
                       + id_67[k];

            t_112[k] = ab_x[k] * hd_68[k]
                       + id_68[k];

            t_113[k] = ab_x[k] * hd_69[k]
                       + id_69[k];

            t_114[k] = ab_x[k] * hd_70[k]
                       + id_70[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, ab_y, ab_z, hd_69, hd_70, \
                         hd_71, id_71, id_99, id_100, id_101, id_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_x[k] * hd_71[k]
                       + id_71[k];

            t_116[k] = ab_y[k] * hd_69[k]
                       + id_99[k];

            t_117[k] = ab_y[k] * hd_70[k]
                       + id_100[k];

            t_118[k] = ab_y[k] * hd_71[k]
                       + id_101[k];

            t_119[k] = ab_z[k] * hd_71[k]
                       + id_107[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, hd_72, hd_73, hd_74, hd_75, \
                         hd_76, id_72, id_73, id_74, id_75, id_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * hd_72[k]
                       + id_72[k];

            t_121[k] = ab_x[k] * hd_73[k]
                       + id_73[k];

            t_122[k] = ab_x[k] * hd_74[k]
                       + id_74[k];

            t_123[k] = ab_x[k] * hd_75[k]
                       + id_75[k];

            t_124[k] = ab_x[k] * hd_76[k]
                       + id_76[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, ab_y, ab_z, hd_75, hd_76, \
                         hd_77, id_77, id_105, id_106, id_107, id_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * hd_77[k]
                       + id_77[k];

            t_126[k] = ab_y[k] * hd_75[k]
                       + id_105[k];

            t_127[k] = ab_y[k] * hd_76[k]
                       + id_106[k];

            t_128[k] = ab_y[k] * hd_77[k]
                       + id_107[k];

            t_129[k] = ab_z[k] * hd_77[k]
                       + id_113[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, hd_78, hd_79, hd_80, hd_81, \
                         hd_82, id_78, id_79, id_80, id_81, id_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_x[k] * hd_78[k]
                       + id_78[k];

            t_131[k] = ab_x[k] * hd_79[k]
                       + id_79[k];

            t_132[k] = ab_x[k] * hd_80[k]
                       + id_80[k];

            t_133[k] = ab_x[k] * hd_81[k]
                       + id_81[k];

            t_134[k] = ab_x[k] * hd_82[k]
                       + id_82[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, ab_y, ab_z, hd_81, hd_82, \
                         hd_83, id_83, id_111, id_112, id_113, id_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * hd_83[k]
                       + id_83[k];

            t_136[k] = ab_y[k] * hd_81[k]
                       + id_111[k];

            t_137[k] = ab_y[k] * hd_82[k]
                       + id_112[k];

            t_138[k] = ab_y[k] * hd_83[k]
                       + id_113[k];

            t_139[k] = ab_z[k] * hd_83[k]
                       + id_119[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, hd_84, hd_85, hd_86, hd_87, \
                         hd_88, id_84, id_85, id_86, id_87, id_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * hd_84[k]
                       + id_84[k];

            t_141[k] = ab_x[k] * hd_85[k]
                       + id_85[k];

            t_142[k] = ab_x[k] * hd_86[k]
                       + id_86[k];

            t_143[k] = ab_x[k] * hd_87[k]
                       + id_87[k];

            t_144[k] = ab_x[k] * hd_88[k]
                       + id_88[k];
        }
    }
}

static auto
compute_hrr_hf_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t hd, const size_t id,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *hd_87 = buffer.data(hd + 87 * ncomps + c);
        const auto *hd_88 = buffer.data(hd + 88 * ncomps + c);
        const auto *hd_89 = buffer.data(hd + 89 * ncomps + c);
        const auto *hd_90 = buffer.data(hd + 90 * ncomps + c);
        const auto *hd_91 = buffer.data(hd + 91 * ncomps + c);
        const auto *hd_92 = buffer.data(hd + 92 * ncomps + c);
        const auto *hd_93 = buffer.data(hd + 93 * ncomps + c);
        const auto *hd_94 = buffer.data(hd + 94 * ncomps + c);
        const auto *hd_95 = buffer.data(hd + 95 * ncomps + c);
        const auto *hd_96 = buffer.data(hd + 96 * ncomps + c);
        const auto *hd_97 = buffer.data(hd + 97 * ncomps + c);
        const auto *hd_98 = buffer.data(hd + 98 * ncomps + c);
        const auto *hd_99 = buffer.data(hd + 99 * ncomps + c);
        const auto *hd_100 = buffer.data(hd + 100 * ncomps + c);
        const auto *hd_101 = buffer.data(hd + 101 * ncomps + c);
        const auto *hd_102 = buffer.data(hd + 102 * ncomps + c);
        const auto *hd_103 = buffer.data(hd + 103 * ncomps + c);
        const auto *hd_104 = buffer.data(hd + 104 * ncomps + c);
        const auto *hd_105 = buffer.data(hd + 105 * ncomps + c);
        const auto *hd_106 = buffer.data(hd + 106 * ncomps + c);
        const auto *hd_107 = buffer.data(hd + 107 * ncomps + c);
        const auto *hd_108 = buffer.data(hd + 108 * ncomps + c);
        const auto *hd_109 = buffer.data(hd + 109 * ncomps + c);
        const auto *hd_110 = buffer.data(hd + 110 * ncomps + c);
        const auto *hd_111 = buffer.data(hd + 111 * ncomps + c);
        const auto *hd_112 = buffer.data(hd + 112 * ncomps + c);
        const auto *hd_113 = buffer.data(hd + 113 * ncomps + c);
        const auto *hd_114 = buffer.data(hd + 114 * ncomps + c);
        const auto *hd_115 = buffer.data(hd + 115 * ncomps + c);
        const auto *hd_116 = buffer.data(hd + 116 * ncomps + c);
        const auto *hd_117 = buffer.data(hd + 117 * ncomps + c);
        const auto *hd_118 = buffer.data(hd + 118 * ncomps + c);
        const auto *hd_119 = buffer.data(hd + 119 * ncomps + c);
        const auto *hd_120 = buffer.data(hd + 120 * ncomps + c);
        const auto *hd_121 = buffer.data(hd + 121 * ncomps + c);
        const auto *hd_122 = buffer.data(hd + 122 * ncomps + c);
        const auto *hd_123 = buffer.data(hd + 123 * ncomps + c);
        const auto *hd_124 = buffer.data(hd + 124 * ncomps + c);
        const auto *hd_125 = buffer.data(hd + 125 * ncomps + c);

        const auto *id_89 = buffer.data(id + 89 * ncomps + c);
        const auto *id_90 = buffer.data(id + 90 * ncomps + c);
        const auto *id_91 = buffer.data(id + 91 * ncomps + c);
        const auto *id_92 = buffer.data(id + 92 * ncomps + c);
        const auto *id_93 = buffer.data(id + 93 * ncomps + c);
        const auto *id_94 = buffer.data(id + 94 * ncomps + c);
        const auto *id_95 = buffer.data(id + 95 * ncomps + c);
        const auto *id_96 = buffer.data(id + 96 * ncomps + c);
        const auto *id_97 = buffer.data(id + 97 * ncomps + c);
        const auto *id_98 = buffer.data(id + 98 * ncomps + c);
        const auto *id_99 = buffer.data(id + 99 * ncomps + c);
        const auto *id_100 = buffer.data(id + 100 * ncomps + c);
        const auto *id_101 = buffer.data(id + 101 * ncomps + c);
        const auto *id_102 = buffer.data(id + 102 * ncomps + c);
        const auto *id_103 = buffer.data(id + 103 * ncomps + c);
        const auto *id_104 = buffer.data(id + 104 * ncomps + c);
        const auto *id_105 = buffer.data(id + 105 * ncomps + c);
        const auto *id_106 = buffer.data(id + 106 * ncomps + c);
        const auto *id_107 = buffer.data(id + 107 * ncomps + c);
        const auto *id_108 = buffer.data(id + 108 * ncomps + c);
        const auto *id_109 = buffer.data(id + 109 * ncomps + c);
        const auto *id_110 = buffer.data(id + 110 * ncomps + c);
        const auto *id_111 = buffer.data(id + 111 * ncomps + c);
        const auto *id_112 = buffer.data(id + 112 * ncomps + c);
        const auto *id_113 = buffer.data(id + 113 * ncomps + c);
        const auto *id_114 = buffer.data(id + 114 * ncomps + c);
        const auto *id_115 = buffer.data(id + 115 * ncomps + c);
        const auto *id_116 = buffer.data(id + 116 * ncomps + c);
        const auto *id_117 = buffer.data(id + 117 * ncomps + c);
        const auto *id_118 = buffer.data(id + 118 * ncomps + c);
        const auto *id_119 = buffer.data(id + 119 * ncomps + c);
        const auto *id_120 = buffer.data(id + 120 * ncomps + c);
        const auto *id_121 = buffer.data(id + 121 * ncomps + c);
        const auto *id_122 = buffer.data(id + 122 * ncomps + c);
        const auto *id_123 = buffer.data(id + 123 * ncomps + c);
        const auto *id_124 = buffer.data(id + 124 * ncomps + c);
        const auto *id_125 = buffer.data(id + 125 * ncomps + c);
        const auto *id_129 = buffer.data(id + 129 * ncomps + c);
        const auto *id_130 = buffer.data(id + 130 * ncomps + c);
        const auto *id_131 = buffer.data(id + 131 * ncomps + c);
        const auto *id_135 = buffer.data(id + 135 * ncomps + c);
        const auto *id_136 = buffer.data(id + 136 * ncomps + c);
        const auto *id_137 = buffer.data(id + 137 * ncomps + c);
        const auto *id_141 = buffer.data(id + 141 * ncomps + c);
        const auto *id_142 = buffer.data(id + 142 * ncomps + c);
        const auto *id_143 = buffer.data(id + 143 * ncomps + c);
        const auto *id_147 = buffer.data(id + 147 * ncomps + c);
        const auto *id_148 = buffer.data(id + 148 * ncomps + c);
        const auto *id_149 = buffer.data(id + 149 * ncomps + c);
        const auto *id_153 = buffer.data(id + 153 * ncomps + c);
        const auto *id_154 = buffer.data(id + 154 * ncomps + c);
        const auto *id_155 = buffer.data(id + 155 * ncomps + c);
        const auto *id_159 = buffer.data(id + 159 * ncomps + c);
        const auto *id_160 = buffer.data(id + 160 * ncomps + c);
        const auto *id_161 = buffer.data(id + 161 * ncomps + c);
        const auto *id_167 = buffer.data(id + 167 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, hd_87, hd_88, \
                         hd_89, id_89, id_117, id_118, id_119, id_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_x[k] * hd_89[k]
                       + id_89[k];

            t_146[k] = ab_y[k] * hd_87[k]
                       + id_117[k];

            t_147[k] = ab_y[k] * hd_88[k]
                       + id_118[k];

            t_148[k] = ab_y[k] * hd_89[k]
                       + id_119[k];

            t_149[k] = ab_z[k] * hd_89[k]
                       + id_125[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, hd_90, hd_91, hd_92, hd_93, \
                         hd_94, id_90, id_91, id_92, id_93, id_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * hd_90[k]
                       + id_90[k];

            t_151[k] = ab_x[k] * hd_91[k]
                       + id_91[k];

            t_152[k] = ab_x[k] * hd_92[k]
                       + id_92[k];

            t_153[k] = ab_x[k] * hd_93[k]
                       + id_93[k];

            t_154[k] = ab_x[k] * hd_94[k]
                       + id_94[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, ab_y, ab_z, hd_93, hd_94, \
                         hd_95, id_95, id_129, id_130, id_131, id_137 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * hd_95[k]
                       + id_95[k];

            t_156[k] = ab_y[k] * hd_93[k]
                       + id_129[k];

            t_157[k] = ab_y[k] * hd_94[k]
                       + id_130[k];

            t_158[k] = ab_y[k] * hd_95[k]
                       + id_131[k];

            t_159[k] = ab_z[k] * hd_95[k]
                       + id_137[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, hd_96, hd_97, hd_98, hd_99, \
                         hd_100, id_96, id_97, id_98, id_99, id_100 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_x[k] * hd_96[k]
                       + id_96[k];

            t_161[k] = ab_x[k] * hd_97[k]
                       + id_97[k];

            t_162[k] = ab_x[k] * hd_98[k]
                       + id_98[k];

            t_163[k] = ab_x[k] * hd_99[k]
                       + id_99[k];

            t_164[k] = ab_x[k] * hd_100[k]
                       + id_100[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, ab_y, ab_z, hd_99, hd_100, \
                         hd_101, id_101, id_135, id_136, id_137, \
                         id_143 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * hd_101[k]
                       + id_101[k];

            t_166[k] = ab_y[k] * hd_99[k]
                       + id_135[k];

            t_167[k] = ab_y[k] * hd_100[k]
                       + id_136[k];

            t_168[k] = ab_y[k] * hd_101[k]
                       + id_137[k];

            t_169[k] = ab_z[k] * hd_101[k]
                       + id_143[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, hd_102, hd_103, hd_104, \
                         hd_105, hd_106, id_102, id_103, id_104, id_105, \
                         id_106 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = ab_x[k] * hd_102[k]
                       + id_102[k];

            t_171[k] = ab_x[k] * hd_103[k]
                       + id_103[k];

            t_172[k] = ab_x[k] * hd_104[k]
                       + id_104[k];

            t_173[k] = ab_x[k] * hd_105[k]
                       + id_105[k];

            t_174[k] = ab_x[k] * hd_106[k]
                       + id_106[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, hd_105, hd_106, \
                         hd_107, id_107, id_141, id_142, id_143, \
                         id_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_x[k] * hd_107[k]
                       + id_107[k];

            t_176[k] = ab_y[k] * hd_105[k]
                       + id_141[k];

            t_177[k] = ab_y[k] * hd_106[k]
                       + id_142[k];

            t_178[k] = ab_y[k] * hd_107[k]
                       + id_143[k];

            t_179[k] = ab_z[k] * hd_107[k]
                       + id_149[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, hd_108, hd_109, hd_110, \
                         hd_111, hd_112, id_108, id_109, id_110, id_111, \
                         id_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * hd_108[k]
                       + id_108[k];

            t_181[k] = ab_x[k] * hd_109[k]
                       + id_109[k];

            t_182[k] = ab_x[k] * hd_110[k]
                       + id_110[k];

            t_183[k] = ab_x[k] * hd_111[k]
                       + id_111[k];

            t_184[k] = ab_x[k] * hd_112[k]
                       + id_112[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, ab_y, ab_z, hd_111, hd_112, \
                         hd_113, id_113, id_147, id_148, id_149, \
                         id_155 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * hd_113[k]
                       + id_113[k];

            t_186[k] = ab_y[k] * hd_111[k]
                       + id_147[k];

            t_187[k] = ab_y[k] * hd_112[k]
                       + id_148[k];

            t_188[k] = ab_y[k] * hd_113[k]
                       + id_149[k];

            t_189[k] = ab_z[k] * hd_113[k]
                       + id_155[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, hd_114, hd_115, hd_116, \
                         hd_117, hd_118, id_114, id_115, id_116, id_117, \
                         id_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_x[k] * hd_114[k]
                       + id_114[k];

            t_191[k] = ab_x[k] * hd_115[k]
                       + id_115[k];

            t_192[k] = ab_x[k] * hd_116[k]
                       + id_116[k];

            t_193[k] = ab_x[k] * hd_117[k]
                       + id_117[k];

            t_194[k] = ab_x[k] * hd_118[k]
                       + id_118[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, ab_y, ab_z, hd_117, hd_118, \
                         hd_119, id_119, id_153, id_154, id_155, \
                         id_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * hd_119[k]
                       + id_119[k];

            t_196[k] = ab_y[k] * hd_117[k]
                       + id_153[k];

            t_197[k] = ab_y[k] * hd_118[k]
                       + id_154[k];

            t_198[k] = ab_y[k] * hd_119[k]
                       + id_155[k];

            t_199[k] = ab_z[k] * hd_119[k]
                       + id_161[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, hd_120, hd_121, hd_122, \
                         hd_123, hd_124, id_120, id_121, id_122, id_123, \
                         id_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = ab_x[k] * hd_120[k]
                       + id_120[k];

            t_201[k] = ab_x[k] * hd_121[k]
                       + id_121[k];

            t_202[k] = ab_x[k] * hd_122[k]
                       + id_122[k];

            t_203[k] = ab_x[k] * hd_123[k]
                       + id_123[k];

            t_204[k] = ab_x[k] * hd_124[k]
                       + id_124[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, ab_y, ab_z, hd_123, hd_124, \
                         hd_125, id_125, id_159, id_160, id_161, \
                         id_167 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_x[k] * hd_125[k]
                       + id_125[k];

            t_206[k] = ab_y[k] * hd_123[k]
                       + id_159[k];

            t_207[k] = ab_y[k] * hd_124[k]
                       + id_160[k];

            t_208[k] = ab_y[k] * hd_125[k]
                       + id_161[k];

            t_209[k] = ab_z[k] * hd_125[k]
                       + id_167[k];
        }
    }
}

auto
compute_hrr_hf_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t hd, const size_t id,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_hf_out_of_first_piece0(buffer, coordinates, target, hd, id, ncomps, nmax);

    compute_hrr_hf_out_of_first_piece1(buffer, coordinates, target, hd, id, ncomps, nmax);
}

static auto
compute_hrr_hf_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t hd, const size_t id, const size_t ncomps,
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

        const auto *hd_0 = buffer.data(hd + 0 * ncomps + c);
        const auto *hd_1 = buffer.data(hd + 1 * ncomps + c);
        const auto *hd_2 = buffer.data(hd + 2 * ncomps + c);
        const auto *hd_3 = buffer.data(hd + 3 * ncomps + c);
        const auto *hd_4 = buffer.data(hd + 4 * ncomps + c);
        const auto *hd_5 = buffer.data(hd + 5 * ncomps + c);
        const auto *hd_6 = buffer.data(hd + 6 * ncomps + c);
        const auto *hd_7 = buffer.data(hd + 7 * ncomps + c);
        const auto *hd_8 = buffer.data(hd + 8 * ncomps + c);
        const auto *hd_9 = buffer.data(hd + 9 * ncomps + c);
        const auto *hd_10 = buffer.data(hd + 10 * ncomps + c);
        const auto *hd_11 = buffer.data(hd + 11 * ncomps + c);
        const auto *hd_12 = buffer.data(hd + 12 * ncomps + c);
        const auto *hd_13 = buffer.data(hd + 13 * ncomps + c);
        const auto *hd_14 = buffer.data(hd + 14 * ncomps + c);
        const auto *hd_15 = buffer.data(hd + 15 * ncomps + c);
        const auto *hd_16 = buffer.data(hd + 16 * ncomps + c);
        const auto *hd_17 = buffer.data(hd + 17 * ncomps + c);
        const auto *hd_18 = buffer.data(hd + 18 * ncomps + c);
        const auto *hd_19 = buffer.data(hd + 19 * ncomps + c);
        const auto *hd_20 = buffer.data(hd + 20 * ncomps + c);
        const auto *hd_21 = buffer.data(hd + 21 * ncomps + c);
        const auto *hd_22 = buffer.data(hd + 22 * ncomps + c);
        const auto *hd_23 = buffer.data(hd + 23 * ncomps + c);
        const auto *hd_24 = buffer.data(hd + 24 * ncomps + c);
        const auto *hd_25 = buffer.data(hd + 25 * ncomps + c);
        const auto *hd_26 = buffer.data(hd + 26 * ncomps + c);
        const auto *hd_27 = buffer.data(hd + 27 * ncomps + c);
        const auto *hd_28 = buffer.data(hd + 28 * ncomps + c);
        const auto *hd_29 = buffer.data(hd + 29 * ncomps + c);
        const auto *hd_30 = buffer.data(hd + 30 * ncomps + c);
        const auto *hd_31 = buffer.data(hd + 31 * ncomps + c);
        const auto *hd_32 = buffer.data(hd + 32 * ncomps + c);
        const auto *hd_33 = buffer.data(hd + 33 * ncomps + c);
        const auto *hd_34 = buffer.data(hd + 34 * ncomps + c);
        const auto *hd_35 = buffer.data(hd + 35 * ncomps + c);
        const auto *hd_36 = buffer.data(hd + 36 * ncomps + c);
        const auto *hd_37 = buffer.data(hd + 37 * ncomps + c);
        const auto *hd_38 = buffer.data(hd + 38 * ncomps + c);
        const auto *hd_39 = buffer.data(hd + 39 * ncomps + c);
        const auto *hd_40 = buffer.data(hd + 40 * ncomps + c);
        const auto *hd_41 = buffer.data(hd + 41 * ncomps + c);
        const auto *hd_42 = buffer.data(hd + 42 * ncomps + c);
        const auto *hd_43 = buffer.data(hd + 43 * ncomps + c);
        const auto *hd_44 = buffer.data(hd + 44 * ncomps + c);
        const auto *hd_45 = buffer.data(hd + 45 * ncomps + c);
        const auto *hd_46 = buffer.data(hd + 46 * ncomps + c);
        const auto *hd_47 = buffer.data(hd + 47 * ncomps + c);
        const auto *hd_48 = buffer.data(hd + 48 * ncomps + c);
        const auto *hd_49 = buffer.data(hd + 49 * ncomps + c);
        const auto *hd_50 = buffer.data(hd + 50 * ncomps + c);
        const auto *hd_51 = buffer.data(hd + 51 * ncomps + c);
        const auto *hd_52 = buffer.data(hd + 52 * ncomps + c);
        const auto *hd_53 = buffer.data(hd + 53 * ncomps + c);
        const auto *hd_54 = buffer.data(hd + 54 * ncomps + c);
        const auto *hd_55 = buffer.data(hd + 55 * ncomps + c);
        const auto *hd_56 = buffer.data(hd + 56 * ncomps + c);
        const auto *hd_57 = buffer.data(hd + 57 * ncomps + c);
        const auto *hd_58 = buffer.data(hd + 58 * ncomps + c);
        const auto *hd_59 = buffer.data(hd + 59 * ncomps + c);
        const auto *hd_60 = buffer.data(hd + 60 * ncomps + c);
        const auto *hd_61 = buffer.data(hd + 61 * ncomps + c);
        const auto *hd_62 = buffer.data(hd + 62 * ncomps + c);
        const auto *hd_63 = buffer.data(hd + 63 * ncomps + c);
        const auto *hd_64 = buffer.data(hd + 64 * ncomps + c);
        const auto *hd_65 = buffer.data(hd + 65 * ncomps + c);
        const auto *hd_66 = buffer.data(hd + 66 * ncomps + c);
        const auto *hd_67 = buffer.data(hd + 67 * ncomps + c);
        const auto *hd_68 = buffer.data(hd + 68 * ncomps + c);
        const auto *hd_69 = buffer.data(hd + 69 * ncomps + c);
        const auto *hd_70 = buffer.data(hd + 70 * ncomps + c);
        const auto *hd_71 = buffer.data(hd + 71 * ncomps + c);
        const auto *hd_72 = buffer.data(hd + 72 * ncomps + c);
        const auto *hd_73 = buffer.data(hd + 73 * ncomps + c);
        const auto *hd_74 = buffer.data(hd + 74 * ncomps + c);
        const auto *hd_75 = buffer.data(hd + 75 * ncomps + c);
        const auto *hd_76 = buffer.data(hd + 76 * ncomps + c);
        const auto *hd_77 = buffer.data(hd + 77 * ncomps + c);
        const auto *hd_78 = buffer.data(hd + 78 * ncomps + c);
        const auto *hd_79 = buffer.data(hd + 79 * ncomps + c);
        const auto *hd_80 = buffer.data(hd + 80 * ncomps + c);
        const auto *hd_81 = buffer.data(hd + 81 * ncomps + c);
        const auto *hd_82 = buffer.data(hd + 82 * ncomps + c);
        const auto *hd_83 = buffer.data(hd + 83 * ncomps + c);
        const auto *hd_84 = buffer.data(hd + 84 * ncomps + c);
        const auto *hd_85 = buffer.data(hd + 85 * ncomps + c);
        const auto *hd_86 = buffer.data(hd + 86 * ncomps + c);
        const auto *hd_87 = buffer.data(hd + 87 * ncomps + c);
        const auto *hd_88 = buffer.data(hd + 88 * ncomps + c);

        const auto *id_0 = buffer.data(id + 0 * ncomps + c);
        const auto *id_1 = buffer.data(id + 1 * ncomps + c);
        const auto *id_2 = buffer.data(id + 2 * ncomps + c);
        const auto *id_3 = buffer.data(id + 3 * ncomps + c);
        const auto *id_4 = buffer.data(id + 4 * ncomps + c);
        const auto *id_5 = buffer.data(id + 5 * ncomps + c);
        const auto *id_6 = buffer.data(id + 6 * ncomps + c);
        const auto *id_7 = buffer.data(id + 7 * ncomps + c);
        const auto *id_8 = buffer.data(id + 8 * ncomps + c);
        const auto *id_9 = buffer.data(id + 9 * ncomps + c);
        const auto *id_10 = buffer.data(id + 10 * ncomps + c);
        const auto *id_11 = buffer.data(id + 11 * ncomps + c);
        const auto *id_12 = buffer.data(id + 12 * ncomps + c);
        const auto *id_13 = buffer.data(id + 13 * ncomps + c);
        const auto *id_14 = buffer.data(id + 14 * ncomps + c);
        const auto *id_15 = buffer.data(id + 15 * ncomps + c);
        const auto *id_16 = buffer.data(id + 16 * ncomps + c);
        const auto *id_17 = buffer.data(id + 17 * ncomps + c);
        const auto *id_18 = buffer.data(id + 18 * ncomps + c);
        const auto *id_19 = buffer.data(id + 19 * ncomps + c);
        const auto *id_20 = buffer.data(id + 20 * ncomps + c);
        const auto *id_21 = buffer.data(id + 21 * ncomps + c);
        const auto *id_22 = buffer.data(id + 22 * ncomps + c);
        const auto *id_23 = buffer.data(id + 23 * ncomps + c);
        const auto *id_24 = buffer.data(id + 24 * ncomps + c);
        const auto *id_25 = buffer.data(id + 25 * ncomps + c);
        const auto *id_26 = buffer.data(id + 26 * ncomps + c);
        const auto *id_27 = buffer.data(id + 27 * ncomps + c);
        const auto *id_28 = buffer.data(id + 28 * ncomps + c);
        const auto *id_29 = buffer.data(id + 29 * ncomps + c);
        const auto *id_30 = buffer.data(id + 30 * ncomps + c);
        const auto *id_31 = buffer.data(id + 31 * ncomps + c);
        const auto *id_32 = buffer.data(id + 32 * ncomps + c);
        const auto *id_33 = buffer.data(id + 33 * ncomps + c);
        const auto *id_34 = buffer.data(id + 34 * ncomps + c);
        const auto *id_35 = buffer.data(id + 35 * ncomps + c);
        const auto *id_36 = buffer.data(id + 36 * ncomps + c);
        const auto *id_37 = buffer.data(id + 37 * ncomps + c);
        const auto *id_38 = buffer.data(id + 38 * ncomps + c);
        const auto *id_39 = buffer.data(id + 39 * ncomps + c);
        const auto *id_40 = buffer.data(id + 40 * ncomps + c);
        const auto *id_41 = buffer.data(id + 41 * ncomps + c);
        const auto *id_42 = buffer.data(id + 42 * ncomps + c);
        const auto *id_43 = buffer.data(id + 43 * ncomps + c);
        const auto *id_44 = buffer.data(id + 44 * ncomps + c);
        const auto *id_45 = buffer.data(id + 45 * ncomps + c);
        const auto *id_46 = buffer.data(id + 46 * ncomps + c);
        const auto *id_47 = buffer.data(id + 47 * ncomps + c);
        const auto *id_48 = buffer.data(id + 48 * ncomps + c);
        const auto *id_49 = buffer.data(id + 49 * ncomps + c);
        const auto *id_50 = buffer.data(id + 50 * ncomps + c);
        const auto *id_51 = buffer.data(id + 51 * ncomps + c);
        const auto *id_52 = buffer.data(id + 52 * ncomps + c);
        const auto *id_53 = buffer.data(id + 53 * ncomps + c);
        const auto *id_54 = buffer.data(id + 54 * ncomps + c);
        const auto *id_55 = buffer.data(id + 55 * ncomps + c);
        const auto *id_56 = buffer.data(id + 56 * ncomps + c);
        const auto *id_57 = buffer.data(id + 57 * ncomps + c);
        const auto *id_58 = buffer.data(id + 58 * ncomps + c);
        const auto *id_59 = buffer.data(id + 59 * ncomps + c);
        const auto *id_60 = buffer.data(id + 60 * ncomps + c);
        const auto *id_61 = buffer.data(id + 61 * ncomps + c);
        const auto *id_62 = buffer.data(id + 62 * ncomps + c);
        const auto *id_63 = buffer.data(id + 63 * ncomps + c);
        const auto *id_64 = buffer.data(id + 64 * ncomps + c);
        const auto *id_65 = buffer.data(id + 65 * ncomps + c);
        const auto *id_66 = buffer.data(id + 66 * ncomps + c);
        const auto *id_67 = buffer.data(id + 67 * ncomps + c);
        const auto *id_68 = buffer.data(id + 68 * ncomps + c);
        const auto *id_69 = buffer.data(id + 69 * ncomps + c);
        const auto *id_70 = buffer.data(id + 70 * ncomps + c);
        const auto *id_71 = buffer.data(id + 71 * ncomps + c);
        const auto *id_72 = buffer.data(id + 72 * ncomps + c);
        const auto *id_73 = buffer.data(id + 73 * ncomps + c);
        const auto *id_74 = buffer.data(id + 74 * ncomps + c);
        const auto *id_75 = buffer.data(id + 75 * ncomps + c);
        const auto *id_76 = buffer.data(id + 76 * ncomps + c);
        const auto *id_77 = buffer.data(id + 77 * ncomps + c);
        const auto *id_78 = buffer.data(id + 78 * ncomps + c);
        const auto *id_79 = buffer.data(id + 79 * ncomps + c);
        const auto *id_80 = buffer.data(id + 80 * ncomps + c);
        const auto *id_81 = buffer.data(id + 81 * ncomps + c);
        const auto *id_82 = buffer.data(id + 82 * ncomps + c);
        const auto *id_83 = buffer.data(id + 83 * ncomps + c);
        const auto *id_84 = buffer.data(id + 84 * ncomps + c);
        const auto *id_85 = buffer.data(id + 85 * ncomps + c);
        const auto *id_86 = buffer.data(id + 86 * ncomps + c);
        const auto *id_87 = buffer.data(id + 87 * ncomps + c);
        const auto *id_88 = buffer.data(id + 88 * ncomps + c);
        const auto *id_89 = buffer.data(id + 89 * ncomps + c);
        const auto *id_93 = buffer.data(id + 93 * ncomps + c);
        const auto *id_94 = buffer.data(id + 94 * ncomps + c);
        const auto *id_95 = buffer.data(id + 95 * ncomps + c);
        const auto *id_99 = buffer.data(id + 99 * ncomps + c);
        const auto *id_100 = buffer.data(id + 100 * ncomps + c);
        const auto *id_101 = buffer.data(id + 101 * ncomps + c);
        const auto *id_105 = buffer.data(id + 105 * ncomps + c);
        const auto *id_106 = buffer.data(id + 106 * ncomps + c);
        const auto *id_107 = buffer.data(id + 107 * ncomps + c);
        const auto *id_111 = buffer.data(id + 111 * ncomps + c);
        const auto *id_112 = buffer.data(id + 112 * ncomps + c);
        const auto *id_113 = buffer.data(id + 113 * ncomps + c);
        const auto *id_119 = buffer.data(id + 119 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, hd_0, hd_1, hd_2, hd_3, hd_4, id_0, \
                         id_1, id_2, id_3, id_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * hd_0[k]
                     + id_0[k];

            t_1[k] = ab_x[k] * hd_1[k]
                     + id_1[k];

            t_2[k] = ab_x[k] * hd_2[k]
                     + id_2[k];

            t_3[k] = ab_x[k] * hd_3[k]
                     + id_3[k];

            t_4[k] = ab_x[k] * hd_4[k]
                     + id_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, hd_3, hd_4, hd_5, id_5, \
                         id_9, id_10, id_11, id_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * hd_5[k]
                     + id_5[k];

            t_6[k] = ab_y[k] * hd_3[k]
                     + id_9[k];

            t_7[k] = ab_y[k] * hd_4[k]
                     + id_10[k];

            t_8[k] = ab_y[k] * hd_5[k]
                     + id_11[k];

            t_9[k] = ab_z[k] * hd_5[k]
                     + id_17[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, hd_6, hd_7, hd_8, hd_9, hd_10, \
                         id_6, id_7, id_8, id_9, id_10 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * hd_6[k]
                      + id_6[k];

            t_11[k] = ab_x[k] * hd_7[k]
                      + id_7[k];

            t_12[k] = ab_x[k] * hd_8[k]
                      + id_8[k];

            t_13[k] = ab_x[k] * hd_9[k]
                      + id_9[k];

            t_14[k] = ab_x[k] * hd_10[k]
                      + id_10[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, hd_9, hd_10, hd_11, \
                         id_11, id_21, id_22, id_23, id_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * hd_11[k]
                      + id_11[k];

            t_16[k] = ab_y[k] * hd_9[k]
                      + id_21[k];

            t_17[k] = ab_y[k] * hd_10[k]
                      + id_22[k];

            t_18[k] = ab_y[k] * hd_11[k]
                      + id_23[k];

            t_19[k] = ab_z[k] * hd_11[k]
                      + id_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, hd_12, hd_13, hd_14, hd_15, \
                         hd_16, id_12, id_13, id_14, id_15, id_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * hd_12[k]
                      + id_12[k];

            t_21[k] = ab_x[k] * hd_13[k]
                      + id_13[k];

            t_22[k] = ab_x[k] * hd_14[k]
                      + id_14[k];

            t_23[k] = ab_x[k] * hd_15[k]
                      + id_15[k];

            t_24[k] = ab_x[k] * hd_16[k]
                      + id_16[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, hd_15, hd_16, hd_17, \
                         id_17, id_27, id_28, id_29, id_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_x[k] * hd_17[k]
                      + id_17[k];

            t_26[k] = ab_y[k] * hd_15[k]
                      + id_27[k];

            t_27[k] = ab_y[k] * hd_16[k]
                      + id_28[k];

            t_28[k] = ab_y[k] * hd_17[k]
                      + id_29[k];

            t_29[k] = ab_z[k] * hd_17[k]
                      + id_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, hd_18, hd_19, hd_20, hd_21, \
                         hd_22, id_18, id_19, id_20, id_21, id_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * hd_18[k]
                      + id_18[k];

            t_31[k] = ab_x[k] * hd_19[k]
                      + id_19[k];

            t_32[k] = ab_x[k] * hd_20[k]
                      + id_20[k];

            t_33[k] = ab_x[k] * hd_21[k]
                      + id_21[k];

            t_34[k] = ab_x[k] * hd_22[k]
                      + id_22[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ab_y, ab_z, hd_21, hd_22, hd_23, \
                         id_23, id_39, id_40, id_41, id_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * hd_23[k]
                      + id_23[k];

            t_36[k] = ab_y[k] * hd_21[k]
                      + id_39[k];

            t_37[k] = ab_y[k] * hd_22[k]
                      + id_40[k];

            t_38[k] = ab_y[k] * hd_23[k]
                      + id_41[k];

            t_39[k] = ab_z[k] * hd_23[k]
                      + id_47[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, hd_24, hd_25, hd_26, hd_27, \
                         hd_28, id_24, id_25, id_26, id_27, id_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_x[k] * hd_24[k]
                      + id_24[k];

            t_41[k] = ab_x[k] * hd_25[k]
                      + id_25[k];

            t_42[k] = ab_x[k] * hd_26[k]
                      + id_26[k];

            t_43[k] = ab_x[k] * hd_27[k]
                      + id_27[k];

            t_44[k] = ab_x[k] * hd_28[k]
                      + id_28[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, hd_27, hd_28, hd_29, \
                         id_29, id_45, id_46, id_47, id_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * hd_29[k]
                      + id_29[k];

            t_46[k] = ab_y[k] * hd_27[k]
                      + id_45[k];

            t_47[k] = ab_y[k] * hd_28[k]
                      + id_46[k];

            t_48[k] = ab_y[k] * hd_29[k]
                      + id_47[k];

            t_49[k] = ab_z[k] * hd_29[k]
                      + id_53[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, hd_30, hd_31, hd_32, hd_33, \
                         hd_34, id_30, id_31, id_32, id_33, id_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * hd_30[k]
                      + id_30[k];

            t_51[k] = ab_x[k] * hd_31[k]
                      + id_31[k];

            t_52[k] = ab_x[k] * hd_32[k]
                      + id_32[k];

            t_53[k] = ab_x[k] * hd_33[k]
                      + id_33[k];

            t_54[k] = ab_x[k] * hd_34[k]
                      + id_34[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, hd_33, hd_34, hd_35, \
                         id_35, id_51, id_52, id_53, id_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_x[k] * hd_35[k]
                      + id_35[k];

            t_56[k] = ab_y[k] * hd_33[k]
                      + id_51[k];

            t_57[k] = ab_y[k] * hd_34[k]
                      + id_52[k];

            t_58[k] = ab_y[k] * hd_35[k]
                      + id_53[k];

            t_59[k] = ab_z[k] * hd_35[k]
                      + id_59[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, hd_36, hd_37, hd_38, hd_39, \
                         hd_40, id_36, id_37, id_38, id_39, id_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * hd_36[k]
                      + id_36[k];

            t_61[k] = ab_x[k] * hd_37[k]
                      + id_37[k];

            t_62[k] = ab_x[k] * hd_38[k]
                      + id_38[k];

            t_63[k] = ab_x[k] * hd_39[k]
                      + id_39[k];

            t_64[k] = ab_x[k] * hd_40[k]
                      + id_40[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, ab_y, ab_z, hd_39, hd_40, hd_41, \
                         id_41, id_63, id_64, id_65, id_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * hd_41[k]
                      + id_41[k];

            t_66[k] = ab_y[k] * hd_39[k]
                      + id_63[k];

            t_67[k] = ab_y[k] * hd_40[k]
                      + id_64[k];

            t_68[k] = ab_y[k] * hd_41[k]
                      + id_65[k];

            t_69[k] = ab_z[k] * hd_41[k]
                      + id_71[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, hd_42, hd_43, hd_44, hd_45, \
                         hd_46, id_42, id_43, id_44, id_45, id_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_x[k] * hd_42[k]
                      + id_42[k];

            t_71[k] = ab_x[k] * hd_43[k]
                      + id_43[k];

            t_72[k] = ab_x[k] * hd_44[k]
                      + id_44[k];

            t_73[k] = ab_x[k] * hd_45[k]
                      + id_45[k];

            t_74[k] = ab_x[k] * hd_46[k]
                      + id_46[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, hd_45, hd_46, hd_47, \
                         id_47, id_69, id_70, id_71, id_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * hd_47[k]
                      + id_47[k];

            t_76[k] = ab_y[k] * hd_45[k]
                      + id_69[k];

            t_77[k] = ab_y[k] * hd_46[k]
                      + id_70[k];

            t_78[k] = ab_y[k] * hd_47[k]
                      + id_71[k];

            t_79[k] = ab_z[k] * hd_47[k]
                      + id_77[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, hd_48, hd_49, hd_50, hd_51, \
                         hd_52, id_48, id_49, id_50, id_51, id_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * hd_48[k]
                      + id_48[k];

            t_81[k] = ab_x[k] * hd_49[k]
                      + id_49[k];

            t_82[k] = ab_x[k] * hd_50[k]
                      + id_50[k];

            t_83[k] = ab_x[k] * hd_51[k]
                      + id_51[k];

            t_84[k] = ab_x[k] * hd_52[k]
                      + id_52[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, hd_51, hd_52, hd_53, \
                         id_53, id_75, id_76, id_77, id_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * hd_53[k]
                      + id_53[k];

            t_86[k] = ab_y[k] * hd_51[k]
                      + id_75[k];

            t_87[k] = ab_y[k] * hd_52[k]
                      + id_76[k];

            t_88[k] = ab_y[k] * hd_53[k]
                      + id_77[k];

            t_89[k] = ab_z[k] * hd_53[k]
                      + id_83[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, hd_54, hd_55, hd_56, hd_57, \
                         hd_58, id_54, id_55, id_56, id_57, id_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * hd_54[k]
                      + id_54[k];

            t_91[k] = ab_x[k] * hd_55[k]
                      + id_55[k];

            t_92[k] = ab_x[k] * hd_56[k]
                      + id_56[k];

            t_93[k] = ab_x[k] * hd_57[k]
                      + id_57[k];

            t_94[k] = ab_x[k] * hd_58[k]
                      + id_58[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, ab_y, ab_z, hd_57, hd_58, hd_59, \
                         id_59, id_81, id_82, id_83, id_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * hd_59[k]
                      + id_59[k];

            t_96[k] = ab_y[k] * hd_57[k]
                      + id_81[k];

            t_97[k] = ab_y[k] * hd_58[k]
                      + id_82[k];

            t_98[k] = ab_y[k] * hd_59[k]
                      + id_83[k];

            t_99[k] = ab_z[k] * hd_59[k]
                      + id_89[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, hd_60, hd_61, hd_62, hd_63, \
                         hd_64, id_60, id_61, id_62, id_63, id_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_x[k] * hd_60[k]
                       + id_60[k];

            t_101[k] = ab_x[k] * hd_61[k]
                       + id_61[k];

            t_102[k] = ab_x[k] * hd_62[k]
                       + id_62[k];

            t_103[k] = ab_x[k] * hd_63[k]
                       + id_63[k];

            t_104[k] = ab_x[k] * hd_64[k]
                       + id_64[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, hd_63, hd_64, \
                         hd_65, id_65, id_93, id_94, id_95, id_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * hd_65[k]
                       + id_65[k];

            t_106[k] = ab_y[k] * hd_63[k]
                       + id_93[k];

            t_107[k] = ab_y[k] * hd_64[k]
                       + id_94[k];

            t_108[k] = ab_y[k] * hd_65[k]
                       + id_95[k];

            t_109[k] = ab_z[k] * hd_65[k]
                       + id_101[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, hd_66, hd_67, hd_68, hd_69, \
                         hd_70, id_66, id_67, id_68, id_69, id_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_x[k] * hd_66[k]
                       + id_66[k];

            t_111[k] = ab_x[k] * hd_67[k]
                       + id_67[k];

            t_112[k] = ab_x[k] * hd_68[k]
                       + id_68[k];

            t_113[k] = ab_x[k] * hd_69[k]
                       + id_69[k];

            t_114[k] = ab_x[k] * hd_70[k]
                       + id_70[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, ab_y, ab_z, hd_69, hd_70, \
                         hd_71, id_71, id_99, id_100, id_101, id_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_x[k] * hd_71[k]
                       + id_71[k];

            t_116[k] = ab_y[k] * hd_69[k]
                       + id_99[k];

            t_117[k] = ab_y[k] * hd_70[k]
                       + id_100[k];

            t_118[k] = ab_y[k] * hd_71[k]
                       + id_101[k];

            t_119[k] = ab_z[k] * hd_71[k]
                       + id_107[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, hd_72, hd_73, hd_74, hd_75, \
                         hd_76, id_72, id_73, id_74, id_75, id_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * hd_72[k]
                       + id_72[k];

            t_121[k] = ab_x[k] * hd_73[k]
                       + id_73[k];

            t_122[k] = ab_x[k] * hd_74[k]
                       + id_74[k];

            t_123[k] = ab_x[k] * hd_75[k]
                       + id_75[k];

            t_124[k] = ab_x[k] * hd_76[k]
                       + id_76[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, ab_y, ab_z, hd_75, hd_76, \
                         hd_77, id_77, id_105, id_106, id_107, id_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * hd_77[k]
                       + id_77[k];

            t_126[k] = ab_y[k] * hd_75[k]
                       + id_105[k];

            t_127[k] = ab_y[k] * hd_76[k]
                       + id_106[k];

            t_128[k] = ab_y[k] * hd_77[k]
                       + id_107[k];

            t_129[k] = ab_z[k] * hd_77[k]
                       + id_113[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, hd_78, hd_79, hd_80, hd_81, \
                         hd_82, id_78, id_79, id_80, id_81, id_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_x[k] * hd_78[k]
                       + id_78[k];

            t_131[k] = ab_x[k] * hd_79[k]
                       + id_79[k];

            t_132[k] = ab_x[k] * hd_80[k]
                       + id_80[k];

            t_133[k] = ab_x[k] * hd_81[k]
                       + id_81[k];

            t_134[k] = ab_x[k] * hd_82[k]
                       + id_82[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, ab_y, ab_z, hd_81, hd_82, \
                         hd_83, id_83, id_111, id_112, id_113, id_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * hd_83[k]
                       + id_83[k];

            t_136[k] = ab_y[k] * hd_81[k]
                       + id_111[k];

            t_137[k] = ab_y[k] * hd_82[k]
                       + id_112[k];

            t_138[k] = ab_y[k] * hd_83[k]
                       + id_113[k];

            t_139[k] = ab_z[k] * hd_83[k]
                       + id_119[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, hd_84, hd_85, hd_86, hd_87, \
                         hd_88, id_84, id_85, id_86, id_87, id_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * hd_84[k]
                       + id_84[k];

            t_141[k] = ab_x[k] * hd_85[k]
                       + id_85[k];

            t_142[k] = ab_x[k] * hd_86[k]
                       + id_86[k];

            t_143[k] = ab_x[k] * hd_87[k]
                       + id_87[k];

            t_144[k] = ab_x[k] * hd_88[k]
                       + id_88[k];
        }
    }
}

static auto
compute_hrr_hf_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t hd, const size_t id, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *hd_87 = buffer.data(hd + 87 * ncomps + c);
        const auto *hd_88 = buffer.data(hd + 88 * ncomps + c);
        const auto *hd_89 = buffer.data(hd + 89 * ncomps + c);
        const auto *hd_90 = buffer.data(hd + 90 * ncomps + c);
        const auto *hd_91 = buffer.data(hd + 91 * ncomps + c);
        const auto *hd_92 = buffer.data(hd + 92 * ncomps + c);
        const auto *hd_93 = buffer.data(hd + 93 * ncomps + c);
        const auto *hd_94 = buffer.data(hd + 94 * ncomps + c);
        const auto *hd_95 = buffer.data(hd + 95 * ncomps + c);
        const auto *hd_96 = buffer.data(hd + 96 * ncomps + c);
        const auto *hd_97 = buffer.data(hd + 97 * ncomps + c);
        const auto *hd_98 = buffer.data(hd + 98 * ncomps + c);
        const auto *hd_99 = buffer.data(hd + 99 * ncomps + c);
        const auto *hd_100 = buffer.data(hd + 100 * ncomps + c);
        const auto *hd_101 = buffer.data(hd + 101 * ncomps + c);
        const auto *hd_102 = buffer.data(hd + 102 * ncomps + c);
        const auto *hd_103 = buffer.data(hd + 103 * ncomps + c);
        const auto *hd_104 = buffer.data(hd + 104 * ncomps + c);
        const auto *hd_105 = buffer.data(hd + 105 * ncomps + c);
        const auto *hd_106 = buffer.data(hd + 106 * ncomps + c);
        const auto *hd_107 = buffer.data(hd + 107 * ncomps + c);
        const auto *hd_108 = buffer.data(hd + 108 * ncomps + c);
        const auto *hd_109 = buffer.data(hd + 109 * ncomps + c);
        const auto *hd_110 = buffer.data(hd + 110 * ncomps + c);
        const auto *hd_111 = buffer.data(hd + 111 * ncomps + c);
        const auto *hd_112 = buffer.data(hd + 112 * ncomps + c);
        const auto *hd_113 = buffer.data(hd + 113 * ncomps + c);
        const auto *hd_114 = buffer.data(hd + 114 * ncomps + c);
        const auto *hd_115 = buffer.data(hd + 115 * ncomps + c);
        const auto *hd_116 = buffer.data(hd + 116 * ncomps + c);
        const auto *hd_117 = buffer.data(hd + 117 * ncomps + c);
        const auto *hd_118 = buffer.data(hd + 118 * ncomps + c);
        const auto *hd_119 = buffer.data(hd + 119 * ncomps + c);
        const auto *hd_120 = buffer.data(hd + 120 * ncomps + c);
        const auto *hd_121 = buffer.data(hd + 121 * ncomps + c);
        const auto *hd_122 = buffer.data(hd + 122 * ncomps + c);
        const auto *hd_123 = buffer.data(hd + 123 * ncomps + c);
        const auto *hd_124 = buffer.data(hd + 124 * ncomps + c);
        const auto *hd_125 = buffer.data(hd + 125 * ncomps + c);

        const auto *id_89 = buffer.data(id + 89 * ncomps + c);
        const auto *id_90 = buffer.data(id + 90 * ncomps + c);
        const auto *id_91 = buffer.data(id + 91 * ncomps + c);
        const auto *id_92 = buffer.data(id + 92 * ncomps + c);
        const auto *id_93 = buffer.data(id + 93 * ncomps + c);
        const auto *id_94 = buffer.data(id + 94 * ncomps + c);
        const auto *id_95 = buffer.data(id + 95 * ncomps + c);
        const auto *id_96 = buffer.data(id + 96 * ncomps + c);
        const auto *id_97 = buffer.data(id + 97 * ncomps + c);
        const auto *id_98 = buffer.data(id + 98 * ncomps + c);
        const auto *id_99 = buffer.data(id + 99 * ncomps + c);
        const auto *id_100 = buffer.data(id + 100 * ncomps + c);
        const auto *id_101 = buffer.data(id + 101 * ncomps + c);
        const auto *id_102 = buffer.data(id + 102 * ncomps + c);
        const auto *id_103 = buffer.data(id + 103 * ncomps + c);
        const auto *id_104 = buffer.data(id + 104 * ncomps + c);
        const auto *id_105 = buffer.data(id + 105 * ncomps + c);
        const auto *id_106 = buffer.data(id + 106 * ncomps + c);
        const auto *id_107 = buffer.data(id + 107 * ncomps + c);
        const auto *id_108 = buffer.data(id + 108 * ncomps + c);
        const auto *id_109 = buffer.data(id + 109 * ncomps + c);
        const auto *id_110 = buffer.data(id + 110 * ncomps + c);
        const auto *id_111 = buffer.data(id + 111 * ncomps + c);
        const auto *id_112 = buffer.data(id + 112 * ncomps + c);
        const auto *id_113 = buffer.data(id + 113 * ncomps + c);
        const auto *id_114 = buffer.data(id + 114 * ncomps + c);
        const auto *id_115 = buffer.data(id + 115 * ncomps + c);
        const auto *id_116 = buffer.data(id + 116 * ncomps + c);
        const auto *id_117 = buffer.data(id + 117 * ncomps + c);
        const auto *id_118 = buffer.data(id + 118 * ncomps + c);
        const auto *id_119 = buffer.data(id + 119 * ncomps + c);
        const auto *id_120 = buffer.data(id + 120 * ncomps + c);
        const auto *id_121 = buffer.data(id + 121 * ncomps + c);
        const auto *id_122 = buffer.data(id + 122 * ncomps + c);
        const auto *id_123 = buffer.data(id + 123 * ncomps + c);
        const auto *id_124 = buffer.data(id + 124 * ncomps + c);
        const auto *id_125 = buffer.data(id + 125 * ncomps + c);
        const auto *id_129 = buffer.data(id + 129 * ncomps + c);
        const auto *id_130 = buffer.data(id + 130 * ncomps + c);
        const auto *id_131 = buffer.data(id + 131 * ncomps + c);
        const auto *id_135 = buffer.data(id + 135 * ncomps + c);
        const auto *id_136 = buffer.data(id + 136 * ncomps + c);
        const auto *id_137 = buffer.data(id + 137 * ncomps + c);
        const auto *id_141 = buffer.data(id + 141 * ncomps + c);
        const auto *id_142 = buffer.data(id + 142 * ncomps + c);
        const auto *id_143 = buffer.data(id + 143 * ncomps + c);
        const auto *id_147 = buffer.data(id + 147 * ncomps + c);
        const auto *id_148 = buffer.data(id + 148 * ncomps + c);
        const auto *id_149 = buffer.data(id + 149 * ncomps + c);
        const auto *id_153 = buffer.data(id + 153 * ncomps + c);
        const auto *id_154 = buffer.data(id + 154 * ncomps + c);
        const auto *id_155 = buffer.data(id + 155 * ncomps + c);
        const auto *id_159 = buffer.data(id + 159 * ncomps + c);
        const auto *id_160 = buffer.data(id + 160 * ncomps + c);
        const auto *id_161 = buffer.data(id + 161 * ncomps + c);
        const auto *id_167 = buffer.data(id + 167 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, hd_87, hd_88, \
                         hd_89, id_89, id_117, id_118, id_119, id_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_x[k] * hd_89[k]
                       + id_89[k];

            t_146[k] = ab_y[k] * hd_87[k]
                       + id_117[k];

            t_147[k] = ab_y[k] * hd_88[k]
                       + id_118[k];

            t_148[k] = ab_y[k] * hd_89[k]
                       + id_119[k];

            t_149[k] = ab_z[k] * hd_89[k]
                       + id_125[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, hd_90, hd_91, hd_92, hd_93, \
                         hd_94, id_90, id_91, id_92, id_93, id_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * hd_90[k]
                       + id_90[k];

            t_151[k] = ab_x[k] * hd_91[k]
                       + id_91[k];

            t_152[k] = ab_x[k] * hd_92[k]
                       + id_92[k];

            t_153[k] = ab_x[k] * hd_93[k]
                       + id_93[k];

            t_154[k] = ab_x[k] * hd_94[k]
                       + id_94[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, ab_y, ab_z, hd_93, hd_94, \
                         hd_95, id_95, id_129, id_130, id_131, id_137 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * hd_95[k]
                       + id_95[k];

            t_156[k] = ab_y[k] * hd_93[k]
                       + id_129[k];

            t_157[k] = ab_y[k] * hd_94[k]
                       + id_130[k];

            t_158[k] = ab_y[k] * hd_95[k]
                       + id_131[k];

            t_159[k] = ab_z[k] * hd_95[k]
                       + id_137[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, hd_96, hd_97, hd_98, hd_99, \
                         hd_100, id_96, id_97, id_98, id_99, id_100 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_x[k] * hd_96[k]
                       + id_96[k];

            t_161[k] = ab_x[k] * hd_97[k]
                       + id_97[k];

            t_162[k] = ab_x[k] * hd_98[k]
                       + id_98[k];

            t_163[k] = ab_x[k] * hd_99[k]
                       + id_99[k];

            t_164[k] = ab_x[k] * hd_100[k]
                       + id_100[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, ab_y, ab_z, hd_99, hd_100, \
                         hd_101, id_101, id_135, id_136, id_137, \
                         id_143 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * hd_101[k]
                       + id_101[k];

            t_166[k] = ab_y[k] * hd_99[k]
                       + id_135[k];

            t_167[k] = ab_y[k] * hd_100[k]
                       + id_136[k];

            t_168[k] = ab_y[k] * hd_101[k]
                       + id_137[k];

            t_169[k] = ab_z[k] * hd_101[k]
                       + id_143[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, hd_102, hd_103, hd_104, \
                         hd_105, hd_106, id_102, id_103, id_104, id_105, \
                         id_106 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = ab_x[k] * hd_102[k]
                       + id_102[k];

            t_171[k] = ab_x[k] * hd_103[k]
                       + id_103[k];

            t_172[k] = ab_x[k] * hd_104[k]
                       + id_104[k];

            t_173[k] = ab_x[k] * hd_105[k]
                       + id_105[k];

            t_174[k] = ab_x[k] * hd_106[k]
                       + id_106[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, hd_105, hd_106, \
                         hd_107, id_107, id_141, id_142, id_143, \
                         id_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_x[k] * hd_107[k]
                       + id_107[k];

            t_176[k] = ab_y[k] * hd_105[k]
                       + id_141[k];

            t_177[k] = ab_y[k] * hd_106[k]
                       + id_142[k];

            t_178[k] = ab_y[k] * hd_107[k]
                       + id_143[k];

            t_179[k] = ab_z[k] * hd_107[k]
                       + id_149[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, hd_108, hd_109, hd_110, \
                         hd_111, hd_112, id_108, id_109, id_110, id_111, \
                         id_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * hd_108[k]
                       + id_108[k];

            t_181[k] = ab_x[k] * hd_109[k]
                       + id_109[k];

            t_182[k] = ab_x[k] * hd_110[k]
                       + id_110[k];

            t_183[k] = ab_x[k] * hd_111[k]
                       + id_111[k];

            t_184[k] = ab_x[k] * hd_112[k]
                       + id_112[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, ab_y, ab_z, hd_111, hd_112, \
                         hd_113, id_113, id_147, id_148, id_149, \
                         id_155 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * hd_113[k]
                       + id_113[k];

            t_186[k] = ab_y[k] * hd_111[k]
                       + id_147[k];

            t_187[k] = ab_y[k] * hd_112[k]
                       + id_148[k];

            t_188[k] = ab_y[k] * hd_113[k]
                       + id_149[k];

            t_189[k] = ab_z[k] * hd_113[k]
                       + id_155[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, hd_114, hd_115, hd_116, \
                         hd_117, hd_118, id_114, id_115, id_116, id_117, \
                         id_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_x[k] * hd_114[k]
                       + id_114[k];

            t_191[k] = ab_x[k] * hd_115[k]
                       + id_115[k];

            t_192[k] = ab_x[k] * hd_116[k]
                       + id_116[k];

            t_193[k] = ab_x[k] * hd_117[k]
                       + id_117[k];

            t_194[k] = ab_x[k] * hd_118[k]
                       + id_118[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, ab_y, ab_z, hd_117, hd_118, \
                         hd_119, id_119, id_153, id_154, id_155, \
                         id_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * hd_119[k]
                       + id_119[k];

            t_196[k] = ab_y[k] * hd_117[k]
                       + id_153[k];

            t_197[k] = ab_y[k] * hd_118[k]
                       + id_154[k];

            t_198[k] = ab_y[k] * hd_119[k]
                       + id_155[k];

            t_199[k] = ab_z[k] * hd_119[k]
                       + id_161[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, hd_120, hd_121, hd_122, \
                         hd_123, hd_124, id_120, id_121, id_122, id_123, \
                         id_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = ab_x[k] * hd_120[k]
                       + id_120[k];

            t_201[k] = ab_x[k] * hd_121[k]
                       + id_121[k];

            t_202[k] = ab_x[k] * hd_122[k]
                       + id_122[k];

            t_203[k] = ab_x[k] * hd_123[k]
                       + id_123[k];

            t_204[k] = ab_x[k] * hd_124[k]
                       + id_124[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, ab_y, ab_z, hd_123, hd_124, \
                         hd_125, id_125, id_159, id_160, id_161, \
                         id_167 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_x[k] * hd_125[k]
                       + id_125[k];

            t_206[k] = ab_y[k] * hd_123[k]
                       + id_159[k];

            t_207[k] = ab_y[k] * hd_124[k]
                       + id_160[k];

            t_208[k] = ab_y[k] * hd_125[k]
                       + id_161[k];

            t_209[k] = ab_z[k] * hd_125[k]
                       + id_167[k];
        }
    }
}

auto
compute_hrr_hf(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t hd, const size_t id, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_hf_piece0(buffer, coordinates, target, hd, id, ncomps, nmax);

    compute_hrr_hf_piece1(buffer, coordinates, target, hd, id, ncomps, nmax);
}

}  // namespace simdtrf
