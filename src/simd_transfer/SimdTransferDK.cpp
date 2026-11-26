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


#include "SimdTransferDK.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_dk_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t pk, const size_t pl, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);

        const auto *pk_0 = buffer.data(pk + 0 * ncomps + c);
        const auto *pk_1 = buffer.data(pk + 1 * ncomps + c);
        const auto *pk_2 = buffer.data(pk + 2 * ncomps + c);
        const auto *pk_3 = buffer.data(pk + 3 * ncomps + c);
        const auto *pk_4 = buffer.data(pk + 4 * ncomps + c);
        const auto *pk_5 = buffer.data(pk + 5 * ncomps + c);
        const auto *pk_6 = buffer.data(pk + 6 * ncomps + c);
        const auto *pk_7 = buffer.data(pk + 7 * ncomps + c);
        const auto *pk_8 = buffer.data(pk + 8 * ncomps + c);
        const auto *pk_9 = buffer.data(pk + 9 * ncomps + c);
        const auto *pk_10 = buffer.data(pk + 10 * ncomps + c);
        const auto *pk_11 = buffer.data(pk + 11 * ncomps + c);
        const auto *pk_12 = buffer.data(pk + 12 * ncomps + c);
        const auto *pk_13 = buffer.data(pk + 13 * ncomps + c);
        const auto *pk_14 = buffer.data(pk + 14 * ncomps + c);
        const auto *pk_15 = buffer.data(pk + 15 * ncomps + c);
        const auto *pk_16 = buffer.data(pk + 16 * ncomps + c);
        const auto *pk_17 = buffer.data(pk + 17 * ncomps + c);
        const auto *pk_18 = buffer.data(pk + 18 * ncomps + c);
        const auto *pk_19 = buffer.data(pk + 19 * ncomps + c);
        const auto *pk_20 = buffer.data(pk + 20 * ncomps + c);
        const auto *pk_21 = buffer.data(pk + 21 * ncomps + c);
        const auto *pk_22 = buffer.data(pk + 22 * ncomps + c);
        const auto *pk_23 = buffer.data(pk + 23 * ncomps + c);
        const auto *pk_24 = buffer.data(pk + 24 * ncomps + c);
        const auto *pk_25 = buffer.data(pk + 25 * ncomps + c);
        const auto *pk_26 = buffer.data(pk + 26 * ncomps + c);
        const auto *pk_27 = buffer.data(pk + 27 * ncomps + c);
        const auto *pk_28 = buffer.data(pk + 28 * ncomps + c);
        const auto *pk_29 = buffer.data(pk + 29 * ncomps + c);
        const auto *pk_30 = buffer.data(pk + 30 * ncomps + c);
        const auto *pk_31 = buffer.data(pk + 31 * ncomps + c);
        const auto *pk_32 = buffer.data(pk + 32 * ncomps + c);
        const auto *pk_33 = buffer.data(pk + 33 * ncomps + c);
        const auto *pk_34 = buffer.data(pk + 34 * ncomps + c);
        const auto *pk_35 = buffer.data(pk + 35 * ncomps + c);
        const auto *pk_36 = buffer.data(pk + 36 * ncomps + c);
        const auto *pk_37 = buffer.data(pk + 37 * ncomps + c);
        const auto *pk_38 = buffer.data(pk + 38 * ncomps + c);
        const auto *pk_39 = buffer.data(pk + 39 * ncomps + c);
        const auto *pk_40 = buffer.data(pk + 40 * ncomps + c);
        const auto *pk_41 = buffer.data(pk + 41 * ncomps + c);
        const auto *pk_42 = buffer.data(pk + 42 * ncomps + c);
        const auto *pk_43 = buffer.data(pk + 43 * ncomps + c);
        const auto *pk_44 = buffer.data(pk + 44 * ncomps + c);
        const auto *pk_45 = buffer.data(pk + 45 * ncomps + c);
        const auto *pk_46 = buffer.data(pk + 46 * ncomps + c);
        const auto *pk_47 = buffer.data(pk + 47 * ncomps + c);
        const auto *pk_48 = buffer.data(pk + 48 * ncomps + c);
        const auto *pk_49 = buffer.data(pk + 49 * ncomps + c);
        const auto *pk_50 = buffer.data(pk + 50 * ncomps + c);
        const auto *pk_51 = buffer.data(pk + 51 * ncomps + c);
        const auto *pk_52 = buffer.data(pk + 52 * ncomps + c);
        const auto *pk_53 = buffer.data(pk + 53 * ncomps + c);
        const auto *pk_54 = buffer.data(pk + 54 * ncomps + c);
        const auto *pk_55 = buffer.data(pk + 55 * ncomps + c);
        const auto *pk_56 = buffer.data(pk + 56 * ncomps + c);
        const auto *pk_57 = buffer.data(pk + 57 * ncomps + c);
        const auto *pk_58 = buffer.data(pk + 58 * ncomps + c);
        const auto *pk_59 = buffer.data(pk + 59 * ncomps + c);
        const auto *pk_60 = buffer.data(pk + 60 * ncomps + c);
        const auto *pk_61 = buffer.data(pk + 61 * ncomps + c);
        const auto *pk_62 = buffer.data(pk + 62 * ncomps + c);
        const auto *pk_63 = buffer.data(pk + 63 * ncomps + c);
        const auto *pk_64 = buffer.data(pk + 64 * ncomps + c);
        const auto *pk_65 = buffer.data(pk + 65 * ncomps + c);
        const auto *pk_66 = buffer.data(pk + 66 * ncomps + c);
        const auto *pk_67 = buffer.data(pk + 67 * ncomps + c);
        const auto *pk_68 = buffer.data(pk + 68 * ncomps + c);
        const auto *pk_69 = buffer.data(pk + 69 * ncomps + c);
        const auto *pk_70 = buffer.data(pk + 70 * ncomps + c);
        const auto *pk_71 = buffer.data(pk + 71 * ncomps + c);
        const auto *pk_72 = buffer.data(pk + 72 * ncomps + c);
        const auto *pk_73 = buffer.data(pk + 73 * ncomps + c);
        const auto *pk_74 = buffer.data(pk + 74 * ncomps + c);
        const auto *pk_75 = buffer.data(pk + 75 * ncomps + c);
        const auto *pk_76 = buffer.data(pk + 76 * ncomps + c);
        const auto *pk_77 = buffer.data(pk + 77 * ncomps + c);
        const auto *pk_78 = buffer.data(pk + 78 * ncomps + c);
        const auto *pk_79 = buffer.data(pk + 79 * ncomps + c);
        const auto *pk_80 = buffer.data(pk + 80 * ncomps + c);
        const auto *pk_81 = buffer.data(pk + 81 * ncomps + c);
        const auto *pk_82 = buffer.data(pk + 82 * ncomps + c);
        const auto *pk_83 = buffer.data(pk + 83 * ncomps + c);
        const auto *pk_84 = buffer.data(pk + 84 * ncomps + c);
        const auto *pk_85 = buffer.data(pk + 85 * ncomps + c);
        const auto *pk_86 = buffer.data(pk + 86 * ncomps + c);
        const auto *pk_87 = buffer.data(pk + 87 * ncomps + c);
        const auto *pk_88 = buffer.data(pk + 88 * ncomps + c);
        const auto *pk_89 = buffer.data(pk + 89 * ncomps + c);
        const auto *pk_90 = buffer.data(pk + 90 * ncomps + c);
        const auto *pk_91 = buffer.data(pk + 91 * ncomps + c);
        const auto *pk_92 = buffer.data(pk + 92 * ncomps + c);
        const auto *pk_93 = buffer.data(pk + 93 * ncomps + c);
        const auto *pk_94 = buffer.data(pk + 94 * ncomps + c);
        const auto *pk_95 = buffer.data(pk + 95 * ncomps + c);
        const auto *pk_96 = buffer.data(pk + 96 * ncomps + c);
        const auto *pk_97 = buffer.data(pk + 97 * ncomps + c);
        const auto *pk_98 = buffer.data(pk + 98 * ncomps + c);
        const auto *pk_99 = buffer.data(pk + 99 * ncomps + c);
        const auto *pk_100 = buffer.data(pk + 100 * ncomps + c);
        const auto *pk_101 = buffer.data(pk + 101 * ncomps + c);
        const auto *pk_102 = buffer.data(pk + 102 * ncomps + c);
        const auto *pk_103 = buffer.data(pk + 103 * ncomps + c);
        const auto *pk_104 = buffer.data(pk + 104 * ncomps + c);
        const auto *pk_105 = buffer.data(pk + 105 * ncomps + c);
        const auto *pk_106 = buffer.data(pk + 106 * ncomps + c);
        const auto *pk_107 = buffer.data(pk + 107 * ncomps + c);

        const auto *pl_0 = buffer.data(pl + 0 * ncomps + c);
        const auto *pl_1 = buffer.data(pl + 1 * ncomps + c);
        const auto *pl_2 = buffer.data(pl + 2 * ncomps + c);
        const auto *pl_3 = buffer.data(pl + 3 * ncomps + c);
        const auto *pl_4 = buffer.data(pl + 4 * ncomps + c);
        const auto *pl_5 = buffer.data(pl + 5 * ncomps + c);
        const auto *pl_6 = buffer.data(pl + 6 * ncomps + c);
        const auto *pl_7 = buffer.data(pl + 7 * ncomps + c);
        const auto *pl_8 = buffer.data(pl + 8 * ncomps + c);
        const auto *pl_9 = buffer.data(pl + 9 * ncomps + c);
        const auto *pl_10 = buffer.data(pl + 10 * ncomps + c);
        const auto *pl_11 = buffer.data(pl + 11 * ncomps + c);
        const auto *pl_12 = buffer.data(pl + 12 * ncomps + c);
        const auto *pl_13 = buffer.data(pl + 13 * ncomps + c);
        const auto *pl_14 = buffer.data(pl + 14 * ncomps + c);
        const auto *pl_15 = buffer.data(pl + 15 * ncomps + c);
        const auto *pl_16 = buffer.data(pl + 16 * ncomps + c);
        const auto *pl_17 = buffer.data(pl + 17 * ncomps + c);
        const auto *pl_18 = buffer.data(pl + 18 * ncomps + c);
        const auto *pl_19 = buffer.data(pl + 19 * ncomps + c);
        const auto *pl_20 = buffer.data(pl + 20 * ncomps + c);
        const auto *pl_21 = buffer.data(pl + 21 * ncomps + c);
        const auto *pl_22 = buffer.data(pl + 22 * ncomps + c);
        const auto *pl_23 = buffer.data(pl + 23 * ncomps + c);
        const auto *pl_24 = buffer.data(pl + 24 * ncomps + c);
        const auto *pl_25 = buffer.data(pl + 25 * ncomps + c);
        const auto *pl_26 = buffer.data(pl + 26 * ncomps + c);
        const auto *pl_27 = buffer.data(pl + 27 * ncomps + c);
        const auto *pl_28 = buffer.data(pl + 28 * ncomps + c);
        const auto *pl_29 = buffer.data(pl + 29 * ncomps + c);
        const auto *pl_30 = buffer.data(pl + 30 * ncomps + c);
        const auto *pl_31 = buffer.data(pl + 31 * ncomps + c);
        const auto *pl_32 = buffer.data(pl + 32 * ncomps + c);
        const auto *pl_33 = buffer.data(pl + 33 * ncomps + c);
        const auto *pl_34 = buffer.data(pl + 34 * ncomps + c);
        const auto *pl_35 = buffer.data(pl + 35 * ncomps + c);
        const auto *pl_45 = buffer.data(pl + 45 * ncomps + c);
        const auto *pl_46 = buffer.data(pl + 46 * ncomps + c);
        const auto *pl_47 = buffer.data(pl + 47 * ncomps + c);
        const auto *pl_48 = buffer.data(pl + 48 * ncomps + c);
        const auto *pl_49 = buffer.data(pl + 49 * ncomps + c);
        const auto *pl_50 = buffer.data(pl + 50 * ncomps + c);
        const auto *pl_51 = buffer.data(pl + 51 * ncomps + c);
        const auto *pl_52 = buffer.data(pl + 52 * ncomps + c);
        const auto *pl_53 = buffer.data(pl + 53 * ncomps + c);
        const auto *pl_54 = buffer.data(pl + 54 * ncomps + c);
        const auto *pl_55 = buffer.data(pl + 55 * ncomps + c);
        const auto *pl_56 = buffer.data(pl + 56 * ncomps + c);
        const auto *pl_57 = buffer.data(pl + 57 * ncomps + c);
        const auto *pl_58 = buffer.data(pl + 58 * ncomps + c);
        const auto *pl_59 = buffer.data(pl + 59 * ncomps + c);
        const auto *pl_60 = buffer.data(pl + 60 * ncomps + c);
        const auto *pl_61 = buffer.data(pl + 61 * ncomps + c);
        const auto *pl_62 = buffer.data(pl + 62 * ncomps + c);
        const auto *pl_63 = buffer.data(pl + 63 * ncomps + c);
        const auto *pl_64 = buffer.data(pl + 64 * ncomps + c);
        const auto *pl_65 = buffer.data(pl + 65 * ncomps + c);
        const auto *pl_66 = buffer.data(pl + 66 * ncomps + c);
        const auto *pl_67 = buffer.data(pl + 67 * ncomps + c);
        const auto *pl_68 = buffer.data(pl + 68 * ncomps + c);
        const auto *pl_69 = buffer.data(pl + 69 * ncomps + c);
        const auto *pl_70 = buffer.data(pl + 70 * ncomps + c);
        const auto *pl_71 = buffer.data(pl + 71 * ncomps + c);
        const auto *pl_72 = buffer.data(pl + 72 * ncomps + c);
        const auto *pl_73 = buffer.data(pl + 73 * ncomps + c);
        const auto *pl_74 = buffer.data(pl + 74 * ncomps + c);
        const auto *pl_75 = buffer.data(pl + 75 * ncomps + c);
        const auto *pl_76 = buffer.data(pl + 76 * ncomps + c);
        const auto *pl_77 = buffer.data(pl + 77 * ncomps + c);
        const auto *pl_78 = buffer.data(pl + 78 * ncomps + c);
        const auto *pl_79 = buffer.data(pl + 79 * ncomps + c);
        const auto *pl_80 = buffer.data(pl + 80 * ncomps + c);
        const auto *pl_81 = buffer.data(pl + 81 * ncomps + c);
        const auto *pl_82 = buffer.data(pl + 82 * ncomps + c);
        const auto *pl_83 = buffer.data(pl + 83 * ncomps + c);
        const auto *pl_84 = buffer.data(pl + 84 * ncomps + c);
        const auto *pl_85 = buffer.data(pl + 85 * ncomps + c);
        const auto *pl_86 = buffer.data(pl + 86 * ncomps + c);
        const auto *pl_87 = buffer.data(pl + 87 * ncomps + c);
        const auto *pl_88 = buffer.data(pl + 88 * ncomps + c);
        const auto *pl_90 = buffer.data(pl + 90 * ncomps + c);
        const auto *pl_91 = buffer.data(pl + 91 * ncomps + c);
        const auto *pl_92 = buffer.data(pl + 92 * ncomps + c);
        const auto *pl_93 = buffer.data(pl + 93 * ncomps + c);
        const auto *pl_94 = buffer.data(pl + 94 * ncomps + c);
        const auto *pl_95 = buffer.data(pl + 95 * ncomps + c);
        const auto *pl_96 = buffer.data(pl + 96 * ncomps + c);
        const auto *pl_97 = buffer.data(pl + 97 * ncomps + c);
        const auto *pl_98 = buffer.data(pl + 98 * ncomps + c);
        const auto *pl_99 = buffer.data(pl + 99 * ncomps + c);
        const auto *pl_100 = buffer.data(pl + 100 * ncomps + c);
        const auto *pl_101 = buffer.data(pl + 101 * ncomps + c);
        const auto *pl_102 = buffer.data(pl + 102 * ncomps + c);
        const auto *pl_103 = buffer.data(pl + 103 * ncomps + c);
        const auto *pl_104 = buffer.data(pl + 104 * ncomps + c);
        const auto *pl_105 = buffer.data(pl + 105 * ncomps + c);
        const auto *pl_106 = buffer.data(pl + 106 * ncomps + c);
        const auto *pl_107 = buffer.data(pl + 107 * ncomps + c);
        const auto *pl_108 = buffer.data(pl + 108 * ncomps + c);
        const auto *pl_109 = buffer.data(pl + 109 * ncomps + c);
        const auto *pl_110 = buffer.data(pl + 110 * ncomps + c);
        const auto *pl_111 = buffer.data(pl + 111 * ncomps + c);
        const auto *pl_112 = buffer.data(pl + 112 * ncomps + c);
        const auto *pl_113 = buffer.data(pl + 113 * ncomps + c);
        const auto *pl_114 = buffer.data(pl + 114 * ncomps + c);
        const auto *pl_115 = buffer.data(pl + 115 * ncomps + c);
        const auto *pl_116 = buffer.data(pl + 116 * ncomps + c);
        const auto *pl_117 = buffer.data(pl + 117 * ncomps + c);
        const auto *pl_118 = buffer.data(pl + 118 * ncomps + c);
        const auto *pl_119 = buffer.data(pl + 119 * ncomps + c);
        const auto *pl_120 = buffer.data(pl + 120 * ncomps + c);
        const auto *pl_121 = buffer.data(pl + 121 * ncomps + c);
        const auto *pl_122 = buffer.data(pl + 122 * ncomps + c);
        const auto *pl_123 = buffer.data(pl + 123 * ncomps + c);
        const auto *pl_124 = buffer.data(pl + 124 * ncomps + c);
        const auto *pl_125 = buffer.data(pl + 125 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pk_0, pk_1, pk_2, pk_3, pk_4, pl_0, \
                         pl_1, pl_2, pl_3, pl_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * pk_0[k]
                     + pl_0[k];

            t_1[k] = -ab_x[k] * pk_1[k]
                     + pl_1[k];

            t_2[k] = -ab_x[k] * pk_2[k]
                     + pl_2[k];

            t_3[k] = -ab_x[k] * pk_3[k]
                     + pl_3[k];

            t_4[k] = -ab_x[k] * pk_4[k]
                     + pl_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pk_5, pk_6, pk_7, pk_8, pk_9, pl_5, \
                         pl_6, pl_7, pl_8, pl_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * pk_5[k]
                     + pl_5[k];

            t_6[k] = -ab_x[k] * pk_6[k]
                     + pl_6[k];

            t_7[k] = -ab_x[k] * pk_7[k]
                     + pl_7[k];

            t_8[k] = -ab_x[k] * pk_8[k]
                     + pl_8[k];

            t_9[k] = -ab_x[k] * pk_9[k]
                     + pl_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pk_10, pk_11, pk_12, pk_13, \
                         pk_14, pl_10, pl_11, pl_12, pl_13, pl_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * pk_10[k]
                      + pl_10[k];

            t_11[k] = -ab_x[k] * pk_11[k]
                      + pl_11[k];

            t_12[k] = -ab_x[k] * pk_12[k]
                      + pl_12[k];

            t_13[k] = -ab_x[k] * pk_13[k]
                      + pl_13[k];

            t_14[k] = -ab_x[k] * pk_14[k]
                      + pl_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, pk_15, pk_16, pk_17, pk_18, \
                         pk_19, pl_15, pl_16, pl_17, pl_18, pl_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * pk_15[k]
                      + pl_15[k];

            t_16[k] = -ab_x[k] * pk_16[k]
                      + pl_16[k];

            t_17[k] = -ab_x[k] * pk_17[k]
                      + pl_17[k];

            t_18[k] = -ab_x[k] * pk_18[k]
                      + pl_18[k];

            t_19[k] = -ab_x[k] * pk_19[k]
                      + pl_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, pk_20, pk_21, pk_22, pk_23, \
                         pk_24, pl_20, pl_21, pl_22, pl_23, pl_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * pk_20[k]
                      + pl_20[k];

            t_21[k] = -ab_x[k] * pk_21[k]
                      + pl_21[k];

            t_22[k] = -ab_x[k] * pk_22[k]
                      + pl_22[k];

            t_23[k] = -ab_x[k] * pk_23[k]
                      + pl_23[k];

            t_24[k] = -ab_x[k] * pk_24[k]
                      + pl_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, pk_25, pk_26, pk_27, pk_28, \
                         pk_29, pl_25, pl_26, pl_27, pl_28, pl_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * pk_25[k]
                      + pl_25[k];

            t_26[k] = -ab_x[k] * pk_26[k]
                      + pl_26[k];

            t_27[k] = -ab_x[k] * pk_27[k]
                      + pl_27[k];

            t_28[k] = -ab_x[k] * pk_28[k]
                      + pl_28[k];

            t_29[k] = -ab_x[k] * pk_29[k]
                      + pl_29[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, pk_30, pk_31, pk_32, pk_33, \
                         pk_34, pl_30, pl_31, pl_32, pl_33, pl_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * pk_30[k]
                      + pl_30[k];

            t_31[k] = -ab_x[k] * pk_31[k]
                      + pl_31[k];

            t_32[k] = -ab_x[k] * pk_32[k]
                      + pl_32[k];

            t_33[k] = -ab_x[k] * pk_33[k]
                      + pl_33[k];

            t_34[k] = -ab_x[k] * pk_34[k]
                      + pl_34[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, pk_35, pk_36, pk_37, pk_38, \
                         pk_39, pl_35, pl_45, pl_46, pl_47, pl_48 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * pk_35[k]
                      + pl_35[k];

            t_36[k] = -ab_x[k] * pk_36[k]
                      + pl_45[k];

            t_37[k] = -ab_x[k] * pk_37[k]
                      + pl_46[k];

            t_38[k] = -ab_x[k] * pk_38[k]
                      + pl_47[k];

            t_39[k] = -ab_x[k] * pk_39[k]
                      + pl_48[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, pk_40, pk_41, pk_42, pk_43, \
                         pk_44, pl_49, pl_50, pl_51, pl_52, pl_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * pk_40[k]
                      + pl_49[k];

            t_41[k] = -ab_x[k] * pk_41[k]
                      + pl_50[k];

            t_42[k] = -ab_x[k] * pk_42[k]
                      + pl_51[k];

            t_43[k] = -ab_x[k] * pk_43[k]
                      + pl_52[k];

            t_44[k] = -ab_x[k] * pk_44[k]
                      + pl_53[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, pk_45, pk_46, pk_47, pk_48, \
                         pk_49, pl_54, pl_55, pl_56, pl_57, pl_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * pk_45[k]
                      + pl_54[k];

            t_46[k] = -ab_x[k] * pk_46[k]
                      + pl_55[k];

            t_47[k] = -ab_x[k] * pk_47[k]
                      + pl_56[k];

            t_48[k] = -ab_x[k] * pk_48[k]
                      + pl_57[k];

            t_49[k] = -ab_x[k] * pk_49[k]
                      + pl_58[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, pk_50, pk_51, pk_52, pk_53, \
                         pk_54, pl_59, pl_60, pl_61, pl_62, pl_63 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * pk_50[k]
                      + pl_59[k];

            t_51[k] = -ab_x[k] * pk_51[k]
                      + pl_60[k];

            t_52[k] = -ab_x[k] * pk_52[k]
                      + pl_61[k];

            t_53[k] = -ab_x[k] * pk_53[k]
                      + pl_62[k];

            t_54[k] = -ab_x[k] * pk_54[k]
                      + pl_63[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, pk_55, pk_56, pk_57, pk_58, \
                         pk_59, pl_64, pl_65, pl_66, pl_67, pl_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * pk_55[k]
                      + pl_64[k];

            t_56[k] = -ab_x[k] * pk_56[k]
                      + pl_65[k];

            t_57[k] = -ab_x[k] * pk_57[k]
                      + pl_66[k];

            t_58[k] = -ab_x[k] * pk_58[k]
                      + pl_67[k];

            t_59[k] = -ab_x[k] * pk_59[k]
                      + pl_68[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, pk_60, pk_61, pk_62, pk_63, \
                         pk_64, pl_69, pl_70, pl_71, pl_72, pl_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * pk_60[k]
                      + pl_69[k];

            t_61[k] = -ab_x[k] * pk_61[k]
                      + pl_70[k];

            t_62[k] = -ab_x[k] * pk_62[k]
                      + pl_71[k];

            t_63[k] = -ab_x[k] * pk_63[k]
                      + pl_72[k];

            t_64[k] = -ab_x[k] * pk_64[k]
                      + pl_73[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, pk_65, pk_66, pk_67, pk_68, \
                         pk_69, pl_74, pl_75, pl_76, pl_77, pl_78 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * pk_65[k]
                      + pl_74[k];

            t_66[k] = -ab_x[k] * pk_66[k]
                      + pl_75[k];

            t_67[k] = -ab_x[k] * pk_67[k]
                      + pl_76[k];

            t_68[k] = -ab_x[k] * pk_68[k]
                      + pl_77[k];

            t_69[k] = -ab_x[k] * pk_69[k]
                      + pl_78[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, pk_70, pk_71, pk_72, pk_73, \
                         pk_74, pl_79, pl_80, pl_90, pl_91, pl_92 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * pk_70[k]
                      + pl_79[k];

            t_71[k] = -ab_x[k] * pk_71[k]
                      + pl_80[k];

            t_72[k] = -ab_x[k] * pk_72[k]
                      + pl_90[k];

            t_73[k] = -ab_x[k] * pk_73[k]
                      + pl_91[k];

            t_74[k] = -ab_x[k] * pk_74[k]
                      + pl_92[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, pk_75, pk_76, pk_77, pk_78, \
                         pk_79, pl_93, pl_94, pl_95, pl_96, pl_97 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * pk_75[k]
                      + pl_93[k];

            t_76[k] = -ab_x[k] * pk_76[k]
                      + pl_94[k];

            t_77[k] = -ab_x[k] * pk_77[k]
                      + pl_95[k];

            t_78[k] = -ab_x[k] * pk_78[k]
                      + pl_96[k];

            t_79[k] = -ab_x[k] * pk_79[k]
                      + pl_97[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, pk_80, pk_81, pk_82, pk_83, \
                         pk_84, pl_98, pl_99, pl_100, pl_101, pl_102 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * pk_80[k]
                      + pl_98[k];

            t_81[k] = -ab_x[k] * pk_81[k]
                      + pl_99[k];

            t_82[k] = -ab_x[k] * pk_82[k]
                      + pl_100[k];

            t_83[k] = -ab_x[k] * pk_83[k]
                      + pl_101[k];

            t_84[k] = -ab_x[k] * pk_84[k]
                      + pl_102[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, pk_85, pk_86, pk_87, pk_88, \
                         pk_89, pl_103, pl_104, pl_105, pl_106, \
                         pl_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * pk_85[k]
                      + pl_103[k];

            t_86[k] = -ab_x[k] * pk_86[k]
                      + pl_104[k];

            t_87[k] = -ab_x[k] * pk_87[k]
                      + pl_105[k];

            t_88[k] = -ab_x[k] * pk_88[k]
                      + pl_106[k];

            t_89[k] = -ab_x[k] * pk_89[k]
                      + pl_107[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, pk_90, pk_91, pk_92, pk_93, \
                         pk_94, pl_108, pl_109, pl_110, pl_111, \
                         pl_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * pk_90[k]
                      + pl_108[k];

            t_91[k] = -ab_x[k] * pk_91[k]
                      + pl_109[k];

            t_92[k] = -ab_x[k] * pk_92[k]
                      + pl_110[k];

            t_93[k] = -ab_x[k] * pk_93[k]
                      + pl_111[k];

            t_94[k] = -ab_x[k] * pk_94[k]
                      + pl_112[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, pk_95, pk_96, pk_97, pk_98, \
                         pk_99, pl_113, pl_114, pl_115, pl_116, \
                         pl_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * pk_95[k]
                      + pl_113[k];

            t_96[k] = -ab_x[k] * pk_96[k]
                      + pl_114[k];

            t_97[k] = -ab_x[k] * pk_97[k]
                      + pl_115[k];

            t_98[k] = -ab_x[k] * pk_98[k]
                      + pl_116[k];

            t_99[k] = -ab_x[k] * pk_99[k]
                      + pl_117[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, pk_100, pk_101, pk_102, \
                         pk_103, pk_104, pl_118, pl_119, pl_120, pl_121, \
                         pl_122 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * pk_100[k]
                       + pl_118[k];

            t_101[k] = -ab_x[k] * pk_101[k]
                       + pl_119[k];

            t_102[k] = -ab_x[k] * pk_102[k]
                       + pl_120[k];

            t_103[k] = -ab_x[k] * pk_103[k]
                       + pl_121[k];

            t_104[k] = -ab_x[k] * pk_104[k]
                       + pl_122[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, ab_x, ab_y, pk_36, pk_105, pk_106, \
                         pk_107, pl_46, pl_123, pl_124, pl_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * pk_105[k]
                       + pl_123[k];

            t_106[k] = -ab_x[k] * pk_106[k]
                       + pl_124[k];

            t_107[k] = -ab_x[k] * pk_107[k]
                       + pl_125[k];

            t_108[k] = -ab_y[k] * pk_36[k]
                       + pl_46[k];
        }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, ab_y, pk_37, pk_38, pk_39, pk_40, \
                         pk_41, pl_48, pl_49, pl_51, pl_52, pl_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_109[k] = -ab_y[k] * pk_37[k]
                       + pl_48[k];

            t_110[k] = -ab_y[k] * pk_38[k]
                       + pl_49[k];

            t_111[k] = -ab_y[k] * pk_39[k]
                       + pl_51[k];

            t_112[k] = -ab_y[k] * pk_40[k]
                       + pl_52[k];

            t_113[k] = -ab_y[k] * pk_41[k]
                       + pl_53[k];
        }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, ab_y, pk_42, pk_43, pk_44, pk_45, \
                         pk_46, pl_55, pl_56, pl_57, pl_58, pl_60 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_114[k] = -ab_y[k] * pk_42[k]
                       + pl_55[k];

            t_115[k] = -ab_y[k] * pk_43[k]
                       + pl_56[k];

            t_116[k] = -ab_y[k] * pk_44[k]
                       + pl_57[k];

            t_117[k] = -ab_y[k] * pk_45[k]
                       + pl_58[k];

            t_118[k] = -ab_y[k] * pk_46[k]
                       + pl_60[k];
        }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, ab_y, pk_47, pk_48, pk_49, pk_50, \
                         pk_51, pl_61, pl_62, pl_63, pl_64, pl_66 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_119[k] = -ab_y[k] * pk_47[k]
                       + pl_61[k];

            t_120[k] = -ab_y[k] * pk_48[k]
                       + pl_62[k];

            t_121[k] = -ab_y[k] * pk_49[k]
                       + pl_63[k];

            t_122[k] = -ab_y[k] * pk_50[k]
                       + pl_64[k];

            t_123[k] = -ab_y[k] * pk_51[k]
                       + pl_66[k];
        }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, ab_y, pk_52, pk_53, pk_54, pk_55, \
                         pk_56, pl_67, pl_68, pl_69, pl_70, pl_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_124[k] = -ab_y[k] * pk_52[k]
                       + pl_67[k];

            t_125[k] = -ab_y[k] * pk_53[k]
                       + pl_68[k];

            t_126[k] = -ab_y[k] * pk_54[k]
                       + pl_69[k];

            t_127[k] = -ab_y[k] * pk_55[k]
                       + pl_70[k];

            t_128[k] = -ab_y[k] * pk_56[k]
                       + pl_71[k];
        }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, ab_y, pk_57, pk_58, pk_59, pk_60, \
                         pk_61, pl_73, pl_74, pl_75, pl_76, pl_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_129[k] = -ab_y[k] * pk_57[k]
                       + pl_73[k];

            t_130[k] = -ab_y[k] * pk_58[k]
                       + pl_74[k];

            t_131[k] = -ab_y[k] * pk_59[k]
                       + pl_75[k];

            t_132[k] = -ab_y[k] * pk_60[k]
                       + pl_76[k];

            t_133[k] = -ab_y[k] * pk_61[k]
                       + pl_77[k];
        }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ab_y, pk_62, pk_63, pk_64, pk_65, \
                         pk_66, pl_78, pl_79, pl_81, pl_82, pl_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_134[k] = -ab_y[k] * pk_62[k]
                       + pl_78[k];

            t_135[k] = -ab_y[k] * pk_63[k]
                       + pl_79[k];

            t_136[k] = -ab_y[k] * pk_64[k]
                       + pl_81[k];

            t_137[k] = -ab_y[k] * pk_65[k]
                       + pl_82[k];

            t_138[k] = -ab_y[k] * pk_66[k]
                       + pl_83[k];
        }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_y, pk_67, pk_68, pk_69, pk_70, \
                         pk_71, pl_84, pl_85, pl_86, pl_87, pl_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_139[k] = -ab_y[k] * pk_67[k]
                       + pl_84[k];

            t_140[k] = -ab_y[k] * pk_68[k]
                       + pl_85[k];

            t_141[k] = -ab_y[k] * pk_69[k]
                       + pl_86[k];

            t_142[k] = -ab_y[k] * pk_70[k]
                       + pl_87[k];

            t_143[k] = -ab_y[k] * pk_71[k]
                       + pl_88[k];
        }
    }
}

static auto
compute_hrr_dk_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t pk, const size_t pl, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *pk_72 = buffer.data(pk + 72 * ncomps + c);
        const auto *pk_73 = buffer.data(pk + 73 * ncomps + c);
        const auto *pk_74 = buffer.data(pk + 74 * ncomps + c);
        const auto *pk_75 = buffer.data(pk + 75 * ncomps + c);
        const auto *pk_76 = buffer.data(pk + 76 * ncomps + c);
        const auto *pk_77 = buffer.data(pk + 77 * ncomps + c);
        const auto *pk_78 = buffer.data(pk + 78 * ncomps + c);
        const auto *pk_79 = buffer.data(pk + 79 * ncomps + c);
        const auto *pk_80 = buffer.data(pk + 80 * ncomps + c);
        const auto *pk_81 = buffer.data(pk + 81 * ncomps + c);
        const auto *pk_82 = buffer.data(pk + 82 * ncomps + c);
        const auto *pk_83 = buffer.data(pk + 83 * ncomps + c);
        const auto *pk_84 = buffer.data(pk + 84 * ncomps + c);
        const auto *pk_85 = buffer.data(pk + 85 * ncomps + c);
        const auto *pk_86 = buffer.data(pk + 86 * ncomps + c);
        const auto *pk_87 = buffer.data(pk + 87 * ncomps + c);
        const auto *pk_88 = buffer.data(pk + 88 * ncomps + c);
        const auto *pk_89 = buffer.data(pk + 89 * ncomps + c);
        const auto *pk_90 = buffer.data(pk + 90 * ncomps + c);
        const auto *pk_91 = buffer.data(pk + 91 * ncomps + c);
        const auto *pk_92 = buffer.data(pk + 92 * ncomps + c);
        const auto *pk_93 = buffer.data(pk + 93 * ncomps + c);
        const auto *pk_94 = buffer.data(pk + 94 * ncomps + c);
        const auto *pk_95 = buffer.data(pk + 95 * ncomps + c);
        const auto *pk_96 = buffer.data(pk + 96 * ncomps + c);
        const auto *pk_97 = buffer.data(pk + 97 * ncomps + c);
        const auto *pk_98 = buffer.data(pk + 98 * ncomps + c);
        const auto *pk_99 = buffer.data(pk + 99 * ncomps + c);
        const auto *pk_100 = buffer.data(pk + 100 * ncomps + c);
        const auto *pk_101 = buffer.data(pk + 101 * ncomps + c);
        const auto *pk_102 = buffer.data(pk + 102 * ncomps + c);
        const auto *pk_103 = buffer.data(pk + 103 * ncomps + c);
        const auto *pk_104 = buffer.data(pk + 104 * ncomps + c);
        const auto *pk_105 = buffer.data(pk + 105 * ncomps + c);
        const auto *pk_106 = buffer.data(pk + 106 * ncomps + c);
        const auto *pk_107 = buffer.data(pk + 107 * ncomps + c);

        const auto *pl_91 = buffer.data(pl + 91 * ncomps + c);
        const auto *pl_92 = buffer.data(pl + 92 * ncomps + c);
        const auto *pl_93 = buffer.data(pl + 93 * ncomps + c);
        const auto *pl_94 = buffer.data(pl + 94 * ncomps + c);
        const auto *pl_95 = buffer.data(pl + 95 * ncomps + c);
        const auto *pl_96 = buffer.data(pl + 96 * ncomps + c);
        const auto *pl_97 = buffer.data(pl + 97 * ncomps + c);
        const auto *pl_98 = buffer.data(pl + 98 * ncomps + c);
        const auto *pl_99 = buffer.data(pl + 99 * ncomps + c);
        const auto *pl_100 = buffer.data(pl + 100 * ncomps + c);
        const auto *pl_101 = buffer.data(pl + 101 * ncomps + c);
        const auto *pl_102 = buffer.data(pl + 102 * ncomps + c);
        const auto *pl_103 = buffer.data(pl + 103 * ncomps + c);
        const auto *pl_104 = buffer.data(pl + 104 * ncomps + c);
        const auto *pl_105 = buffer.data(pl + 105 * ncomps + c);
        const auto *pl_106 = buffer.data(pl + 106 * ncomps + c);
        const auto *pl_107 = buffer.data(pl + 107 * ncomps + c);
        const auto *pl_108 = buffer.data(pl + 108 * ncomps + c);
        const auto *pl_109 = buffer.data(pl + 109 * ncomps + c);
        const auto *pl_110 = buffer.data(pl + 110 * ncomps + c);
        const auto *pl_111 = buffer.data(pl + 111 * ncomps + c);
        const auto *pl_112 = buffer.data(pl + 112 * ncomps + c);
        const auto *pl_113 = buffer.data(pl + 113 * ncomps + c);
        const auto *pl_114 = buffer.data(pl + 114 * ncomps + c);
        const auto *pl_115 = buffer.data(pl + 115 * ncomps + c);
        const auto *pl_116 = buffer.data(pl + 116 * ncomps + c);
        const auto *pl_117 = buffer.data(pl + 117 * ncomps + c);
        const auto *pl_118 = buffer.data(pl + 118 * ncomps + c);
        const auto *pl_119 = buffer.data(pl + 119 * ncomps + c);
        const auto *pl_120 = buffer.data(pl + 120 * ncomps + c);
        const auto *pl_121 = buffer.data(pl + 121 * ncomps + c);
        const auto *pl_122 = buffer.data(pl + 122 * ncomps + c);
        const auto *pl_123 = buffer.data(pl + 123 * ncomps + c);
        const auto *pl_124 = buffer.data(pl + 124 * ncomps + c);
        const auto *pl_125 = buffer.data(pl + 125 * ncomps + c);
        const auto *pl_126 = buffer.data(pl + 126 * ncomps + c);
        const auto *pl_127 = buffer.data(pl + 127 * ncomps + c);
        const auto *pl_128 = buffer.data(pl + 128 * ncomps + c);
        const auto *pl_129 = buffer.data(pl + 129 * ncomps + c);
        const auto *pl_130 = buffer.data(pl + 130 * ncomps + c);
        const auto *pl_131 = buffer.data(pl + 131 * ncomps + c);
        const auto *pl_132 = buffer.data(pl + 132 * ncomps + c);
        const auto *pl_133 = buffer.data(pl + 133 * ncomps + c);
        const auto *pl_134 = buffer.data(pl + 134 * ncomps + c);

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, ab_y, pk_72, pk_73, pk_74, pk_75, \
                         pk_76, pl_91, pl_93, pl_94, pl_96, pl_97 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = -ab_y[k] * pk_72[k]
                       + pl_91[k];

            t_145[k] = -ab_y[k] * pk_73[k]
                       + pl_93[k];

            t_146[k] = -ab_y[k] * pk_74[k]
                       + pl_94[k];

            t_147[k] = -ab_y[k] * pk_75[k]
                       + pl_96[k];

            t_148[k] = -ab_y[k] * pk_76[k]
                       + pl_97[k];
        }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, ab_y, pk_77, pk_78, pk_79, pk_80, \
                         pk_81, pl_98, pl_100, pl_101, pl_102, pl_103 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_149[k] = -ab_y[k] * pk_77[k]
                       + pl_98[k];

            t_150[k] = -ab_y[k] * pk_78[k]
                       + pl_100[k];

            t_151[k] = -ab_y[k] * pk_79[k]
                       + pl_101[k];

            t_152[k] = -ab_y[k] * pk_80[k]
                       + pl_102[k];

            t_153[k] = -ab_y[k] * pk_81[k]
                       + pl_103[k];
        }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, ab_y, pk_82, pk_83, pk_84, pk_85, \
                         pk_86, pl_105, pl_106, pl_107, pl_108, \
                         pl_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_154[k] = -ab_y[k] * pk_82[k]
                       + pl_105[k];

            t_155[k] = -ab_y[k] * pk_83[k]
                       + pl_106[k];

            t_156[k] = -ab_y[k] * pk_84[k]
                       + pl_107[k];

            t_157[k] = -ab_y[k] * pk_85[k]
                       + pl_108[k];

            t_158[k] = -ab_y[k] * pk_86[k]
                       + pl_109[k];
        }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, ab_y, pk_87, pk_88, pk_89, pk_90, \
                         pk_91, pl_111, pl_112, pl_113, pl_114, \
                         pl_115 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_159[k] = -ab_y[k] * pk_87[k]
                       + pl_111[k];

            t_160[k] = -ab_y[k] * pk_88[k]
                       + pl_112[k];

            t_161[k] = -ab_y[k] * pk_89[k]
                       + pl_113[k];

            t_162[k] = -ab_y[k] * pk_90[k]
                       + pl_114[k];

            t_163[k] = -ab_y[k] * pk_91[k]
                       + pl_115[k];
        }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, ab_y, pk_92, pk_93, pk_94, pk_95, \
                         pk_96, pl_116, pl_118, pl_119, pl_120, \
                         pl_121 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_164[k] = -ab_y[k] * pk_92[k]
                       + pl_116[k];

            t_165[k] = -ab_y[k] * pk_93[k]
                       + pl_118[k];

            t_166[k] = -ab_y[k] * pk_94[k]
                       + pl_119[k];

            t_167[k] = -ab_y[k] * pk_95[k]
                       + pl_120[k];

            t_168[k] = -ab_y[k] * pk_96[k]
                       + pl_121[k];
        }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, ab_y, pk_97, pk_98, pk_99, pk_100, \
                         pk_101, pl_122, pl_123, pl_124, pl_126, \
                         pl_127 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_169[k] = -ab_y[k] * pk_97[k]
                       + pl_122[k];

            t_170[k] = -ab_y[k] * pk_98[k]
                       + pl_123[k];

            t_171[k] = -ab_y[k] * pk_99[k]
                       + pl_124[k];

            t_172[k] = -ab_y[k] * pk_100[k]
                       + pl_126[k];

            t_173[k] = -ab_y[k] * pk_101[k]
                       + pl_127[k];
        }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, ab_y, pk_102, pk_103, pk_104, \
                         pk_105, pk_106, pl_128, pl_129, pl_130, pl_131, \
                         pl_132 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_174[k] = -ab_y[k] * pk_102[k]
                       + pl_128[k];

            t_175[k] = -ab_y[k] * pk_103[k]
                       + pl_129[k];

            t_176[k] = -ab_y[k] * pk_104[k]
                       + pl_130[k];

            t_177[k] = -ab_y[k] * pk_105[k]
                       + pl_131[k];

            t_178[k] = -ab_y[k] * pk_106[k]
                       + pl_132[k];
        }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, ab_y, ab_z, pk_72, pk_73, pk_74, pk_107, \
                         pl_92, pl_94, pl_95, pl_133 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_179[k] = -ab_y[k] * pk_107[k]
                       + pl_133[k];

            t_180[k] = -ab_z[k] * pk_72[k]
                       + pl_92[k];

            t_181[k] = -ab_z[k] * pk_73[k]
                       + pl_94[k];

            t_182[k] = -ab_z[k] * pk_74[k]
                       + pl_95[k];
        }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, ab_z, pk_75, pk_76, pk_77, pk_78, \
                         pk_79, pl_97, pl_98, pl_99, pl_101, pl_102 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_183[k] = -ab_z[k] * pk_75[k]
                       + pl_97[k];

            t_184[k] = -ab_z[k] * pk_76[k]
                       + pl_98[k];

            t_185[k] = -ab_z[k] * pk_77[k]
                       + pl_99[k];

            t_186[k] = -ab_z[k] * pk_78[k]
                       + pl_101[k];

            t_187[k] = -ab_z[k] * pk_79[k]
                       + pl_102[k];
        }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, ab_z, pk_80, pk_81, pk_82, pk_83, \
                         pk_84, pl_103, pl_104, pl_106, pl_107, \
                         pl_108 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_188[k] = -ab_z[k] * pk_80[k]
                       + pl_103[k];

            t_189[k] = -ab_z[k] * pk_81[k]
                       + pl_104[k];

            t_190[k] = -ab_z[k] * pk_82[k]
                       + pl_106[k];

            t_191[k] = -ab_z[k] * pk_83[k]
                       + pl_107[k];

            t_192[k] = -ab_z[k] * pk_84[k]
                       + pl_108[k];
        }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, ab_z, pk_85, pk_86, pk_87, pk_88, \
                         pk_89, pl_109, pl_110, pl_112, pl_113, \
                         pl_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_193[k] = -ab_z[k] * pk_85[k]
                       + pl_109[k];

            t_194[k] = -ab_z[k] * pk_86[k]
                       + pl_110[k];

            t_195[k] = -ab_z[k] * pk_87[k]
                       + pl_112[k];

            t_196[k] = -ab_z[k] * pk_88[k]
                       + pl_113[k];

            t_197[k] = -ab_z[k] * pk_89[k]
                       + pl_114[k];
        }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, ab_z, pk_90, pk_91, pk_92, pk_93, \
                         pk_94, pl_115, pl_116, pl_117, pl_119, \
                         pl_120 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_198[k] = -ab_z[k] * pk_90[k]
                       + pl_115[k];

            t_199[k] = -ab_z[k] * pk_91[k]
                       + pl_116[k];

            t_200[k] = -ab_z[k] * pk_92[k]
                       + pl_117[k];

            t_201[k] = -ab_z[k] * pk_93[k]
                       + pl_119[k];

            t_202[k] = -ab_z[k] * pk_94[k]
                       + pl_120[k];
        }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, ab_z, pk_95, pk_96, pk_97, pk_98, \
                         pk_99, pl_121, pl_122, pl_123, pl_124, \
                         pl_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_203[k] = -ab_z[k] * pk_95[k]
                       + pl_121[k];

            t_204[k] = -ab_z[k] * pk_96[k]
                       + pl_122[k];

            t_205[k] = -ab_z[k] * pk_97[k]
                       + pl_123[k];

            t_206[k] = -ab_z[k] * pk_98[k]
                       + pl_124[k];

            t_207[k] = -ab_z[k] * pk_99[k]
                       + pl_125[k];
        }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, ab_z, pk_100, pk_101, pk_102, \
                         pk_103, pk_104, pl_127, pl_128, pl_129, pl_130, \
                         pl_131 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_208[k] = -ab_z[k] * pk_100[k]
                       + pl_127[k];

            t_209[k] = -ab_z[k] * pk_101[k]
                       + pl_128[k];

            t_210[k] = -ab_z[k] * pk_102[k]
                       + pl_129[k];

            t_211[k] = -ab_z[k] * pk_103[k]
                       + pl_130[k];

            t_212[k] = -ab_z[k] * pk_104[k]
                       + pl_131[k];
        }

#pragma omp simd aligned(t_213, t_214, t_215, ab_z, pk_105, pk_106, pk_107, pl_132, pl_133, \
                         pl_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_213[k] = -ab_z[k] * pk_105[k]
                       + pl_132[k];

            t_214[k] = -ab_z[k] * pk_106[k]
                       + pl_133[k];

            t_215[k] = -ab_z[k] * pk_107[k]
                       + pl_134[k];
        }
    }
}

auto
compute_hrr_dk(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t pk, const size_t pl, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_dk_piece0(buffer, coordinates, target, pk, pl, ncomps, nmax);

    compute_hrr_dk_piece1(buffer, coordinates, target, pk, pl, ncomps, nmax);
}

}  // namespace simdtrf
