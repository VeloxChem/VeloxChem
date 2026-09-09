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


#include "SimdTransferPO.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_po_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t so, const size_t sq, const size_t ncomps,
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

        const auto *so_0 = buffer.data(so + 0 * ncomps + c);
        const auto *so_1 = buffer.data(so + 1 * ncomps + c);
        const auto *so_2 = buffer.data(so + 2 * ncomps + c);
        const auto *so_3 = buffer.data(so + 3 * ncomps + c);
        const auto *so_4 = buffer.data(so + 4 * ncomps + c);
        const auto *so_5 = buffer.data(so + 5 * ncomps + c);
        const auto *so_6 = buffer.data(so + 6 * ncomps + c);
        const auto *so_7 = buffer.data(so + 7 * ncomps + c);
        const auto *so_8 = buffer.data(so + 8 * ncomps + c);
        const auto *so_9 = buffer.data(so + 9 * ncomps + c);
        const auto *so_10 = buffer.data(so + 10 * ncomps + c);
        const auto *so_11 = buffer.data(so + 11 * ncomps + c);
        const auto *so_12 = buffer.data(so + 12 * ncomps + c);
        const auto *so_13 = buffer.data(so + 13 * ncomps + c);
        const auto *so_14 = buffer.data(so + 14 * ncomps + c);
        const auto *so_15 = buffer.data(so + 15 * ncomps + c);
        const auto *so_16 = buffer.data(so + 16 * ncomps + c);
        const auto *so_17 = buffer.data(so + 17 * ncomps + c);
        const auto *so_18 = buffer.data(so + 18 * ncomps + c);
        const auto *so_19 = buffer.data(so + 19 * ncomps + c);
        const auto *so_20 = buffer.data(so + 20 * ncomps + c);
        const auto *so_21 = buffer.data(so + 21 * ncomps + c);
        const auto *so_22 = buffer.data(so + 22 * ncomps + c);
        const auto *so_23 = buffer.data(so + 23 * ncomps + c);
        const auto *so_24 = buffer.data(so + 24 * ncomps + c);
        const auto *so_25 = buffer.data(so + 25 * ncomps + c);
        const auto *so_26 = buffer.data(so + 26 * ncomps + c);
        const auto *so_27 = buffer.data(so + 27 * ncomps + c);
        const auto *so_28 = buffer.data(so + 28 * ncomps + c);
        const auto *so_29 = buffer.data(so + 29 * ncomps + c);
        const auto *so_30 = buffer.data(so + 30 * ncomps + c);
        const auto *so_31 = buffer.data(so + 31 * ncomps + c);
        const auto *so_32 = buffer.data(so + 32 * ncomps + c);
        const auto *so_33 = buffer.data(so + 33 * ncomps + c);
        const auto *so_34 = buffer.data(so + 34 * ncomps + c);
        const auto *so_35 = buffer.data(so + 35 * ncomps + c);
        const auto *so_36 = buffer.data(so + 36 * ncomps + c);
        const auto *so_37 = buffer.data(so + 37 * ncomps + c);
        const auto *so_38 = buffer.data(so + 38 * ncomps + c);
        const auto *so_39 = buffer.data(so + 39 * ncomps + c);
        const auto *so_40 = buffer.data(so + 40 * ncomps + c);
        const auto *so_41 = buffer.data(so + 41 * ncomps + c);
        const auto *so_42 = buffer.data(so + 42 * ncomps + c);
        const auto *so_43 = buffer.data(so + 43 * ncomps + c);
        const auto *so_44 = buffer.data(so + 44 * ncomps + c);
        const auto *so_45 = buffer.data(so + 45 * ncomps + c);
        const auto *so_46 = buffer.data(so + 46 * ncomps + c);
        const auto *so_47 = buffer.data(so + 47 * ncomps + c);
        const auto *so_48 = buffer.data(so + 48 * ncomps + c);
        const auto *so_49 = buffer.data(so + 49 * ncomps + c);
        const auto *so_50 = buffer.data(so + 50 * ncomps + c);
        const auto *so_51 = buffer.data(so + 51 * ncomps + c);
        const auto *so_52 = buffer.data(so + 52 * ncomps + c);
        const auto *so_53 = buffer.data(so + 53 * ncomps + c);
        const auto *so_54 = buffer.data(so + 54 * ncomps + c);
        const auto *so_55 = buffer.data(so + 55 * ncomps + c);
        const auto *so_56 = buffer.data(so + 56 * ncomps + c);
        const auto *so_57 = buffer.data(so + 57 * ncomps + c);
        const auto *so_58 = buffer.data(so + 58 * ncomps + c);
        const auto *so_59 = buffer.data(so + 59 * ncomps + c);
        const auto *so_60 = buffer.data(so + 60 * ncomps + c);
        const auto *so_61 = buffer.data(so + 61 * ncomps + c);
        const auto *so_62 = buffer.data(so + 62 * ncomps + c);
        const auto *so_63 = buffer.data(so + 63 * ncomps + c);
        const auto *so_64 = buffer.data(so + 64 * ncomps + c);
        const auto *so_65 = buffer.data(so + 65 * ncomps + c);
        const auto *so_66 = buffer.data(so + 66 * ncomps + c);
        const auto *so_67 = buffer.data(so + 67 * ncomps + c);
        const auto *so_68 = buffer.data(so + 68 * ncomps + c);
        const auto *so_69 = buffer.data(so + 69 * ncomps + c);
        const auto *so_70 = buffer.data(so + 70 * ncomps + c);
        const auto *so_71 = buffer.data(so + 71 * ncomps + c);
        const auto *so_72 = buffer.data(so + 72 * ncomps + c);
        const auto *so_73 = buffer.data(so + 73 * ncomps + c);
        const auto *so_74 = buffer.data(so + 74 * ncomps + c);
        const auto *so_75 = buffer.data(so + 75 * ncomps + c);
        const auto *so_76 = buffer.data(so + 76 * ncomps + c);
        const auto *so_77 = buffer.data(so + 77 * ncomps + c);

        const auto *sq_0 = buffer.data(sq + 0 * ncomps + c);
        const auto *sq_1 = buffer.data(sq + 1 * ncomps + c);
        const auto *sq_2 = buffer.data(sq + 2 * ncomps + c);
        const auto *sq_3 = buffer.data(sq + 3 * ncomps + c);
        const auto *sq_4 = buffer.data(sq + 4 * ncomps + c);
        const auto *sq_5 = buffer.data(sq + 5 * ncomps + c);
        const auto *sq_6 = buffer.data(sq + 6 * ncomps + c);
        const auto *sq_7 = buffer.data(sq + 7 * ncomps + c);
        const auto *sq_8 = buffer.data(sq + 8 * ncomps + c);
        const auto *sq_9 = buffer.data(sq + 9 * ncomps + c);
        const auto *sq_10 = buffer.data(sq + 10 * ncomps + c);
        const auto *sq_11 = buffer.data(sq + 11 * ncomps + c);
        const auto *sq_12 = buffer.data(sq + 12 * ncomps + c);
        const auto *sq_13 = buffer.data(sq + 13 * ncomps + c);
        const auto *sq_14 = buffer.data(sq + 14 * ncomps + c);
        const auto *sq_15 = buffer.data(sq + 15 * ncomps + c);
        const auto *sq_16 = buffer.data(sq + 16 * ncomps + c);
        const auto *sq_17 = buffer.data(sq + 17 * ncomps + c);
        const auto *sq_18 = buffer.data(sq + 18 * ncomps + c);
        const auto *sq_19 = buffer.data(sq + 19 * ncomps + c);
        const auto *sq_20 = buffer.data(sq + 20 * ncomps + c);
        const auto *sq_21 = buffer.data(sq + 21 * ncomps + c);
        const auto *sq_22 = buffer.data(sq + 22 * ncomps + c);
        const auto *sq_23 = buffer.data(sq + 23 * ncomps + c);
        const auto *sq_24 = buffer.data(sq + 24 * ncomps + c);
        const auto *sq_25 = buffer.data(sq + 25 * ncomps + c);
        const auto *sq_26 = buffer.data(sq + 26 * ncomps + c);
        const auto *sq_27 = buffer.data(sq + 27 * ncomps + c);
        const auto *sq_28 = buffer.data(sq + 28 * ncomps + c);
        const auto *sq_29 = buffer.data(sq + 29 * ncomps + c);
        const auto *sq_30 = buffer.data(sq + 30 * ncomps + c);
        const auto *sq_31 = buffer.data(sq + 31 * ncomps + c);
        const auto *sq_32 = buffer.data(sq + 32 * ncomps + c);
        const auto *sq_33 = buffer.data(sq + 33 * ncomps + c);
        const auto *sq_34 = buffer.data(sq + 34 * ncomps + c);
        const auto *sq_35 = buffer.data(sq + 35 * ncomps + c);
        const auto *sq_36 = buffer.data(sq + 36 * ncomps + c);
        const auto *sq_37 = buffer.data(sq + 37 * ncomps + c);
        const auto *sq_38 = buffer.data(sq + 38 * ncomps + c);
        const auto *sq_39 = buffer.data(sq + 39 * ncomps + c);
        const auto *sq_40 = buffer.data(sq + 40 * ncomps + c);
        const auto *sq_41 = buffer.data(sq + 41 * ncomps + c);
        const auto *sq_42 = buffer.data(sq + 42 * ncomps + c);
        const auto *sq_43 = buffer.data(sq + 43 * ncomps + c);
        const auto *sq_44 = buffer.data(sq + 44 * ncomps + c);
        const auto *sq_45 = buffer.data(sq + 45 * ncomps + c);
        const auto *sq_46 = buffer.data(sq + 46 * ncomps + c);
        const auto *sq_47 = buffer.data(sq + 47 * ncomps + c);
        const auto *sq_48 = buffer.data(sq + 48 * ncomps + c);
        const auto *sq_49 = buffer.data(sq + 49 * ncomps + c);
        const auto *sq_50 = buffer.data(sq + 50 * ncomps + c);
        const auto *sq_51 = buffer.data(sq + 51 * ncomps + c);
        const auto *sq_52 = buffer.data(sq + 52 * ncomps + c);
        const auto *sq_53 = buffer.data(sq + 53 * ncomps + c);
        const auto *sq_54 = buffer.data(sq + 54 * ncomps + c);
        const auto *sq_55 = buffer.data(sq + 55 * ncomps + c);
        const auto *sq_56 = buffer.data(sq + 56 * ncomps + c);
        const auto *sq_57 = buffer.data(sq + 57 * ncomps + c);
        const auto *sq_58 = buffer.data(sq + 58 * ncomps + c);
        const auto *sq_59 = buffer.data(sq + 59 * ncomps + c);
        const auto *sq_60 = buffer.data(sq + 60 * ncomps + c);
        const auto *sq_61 = buffer.data(sq + 61 * ncomps + c);
        const auto *sq_62 = buffer.data(sq + 62 * ncomps + c);
        const auto *sq_63 = buffer.data(sq + 63 * ncomps + c);
        const auto *sq_64 = buffer.data(sq + 64 * ncomps + c);
        const auto *sq_65 = buffer.data(sq + 65 * ncomps + c);
        const auto *sq_66 = buffer.data(sq + 66 * ncomps + c);
        const auto *sq_67 = buffer.data(sq + 67 * ncomps + c);
        const auto *sq_68 = buffer.data(sq + 68 * ncomps + c);
        const auto *sq_69 = buffer.data(sq + 69 * ncomps + c);
        const auto *sq_70 = buffer.data(sq + 70 * ncomps + c);
        const auto *sq_71 = buffer.data(sq + 71 * ncomps + c);
        const auto *sq_72 = buffer.data(sq + 72 * ncomps + c);
        const auto *sq_73 = buffer.data(sq + 73 * ncomps + c);
        const auto *sq_74 = buffer.data(sq + 74 * ncomps + c);
        const auto *sq_75 = buffer.data(sq + 75 * ncomps + c);
        const auto *sq_76 = buffer.data(sq + 76 * ncomps + c);
        const auto *sq_77 = buffer.data(sq + 77 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, so_0, so_1, so_2, so_3, so_4, sq_0, \
                         sq_1, sq_2, sq_3, sq_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * so_0[k]
                     + sq_0[k];

            t_1[k] = -ab_x[k] * so_1[k]
                     + sq_1[k];

            t_2[k] = -ab_x[k] * so_2[k]
                     + sq_2[k];

            t_3[k] = -ab_x[k] * so_3[k]
                     + sq_3[k];

            t_4[k] = -ab_x[k] * so_4[k]
                     + sq_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, so_5, so_6, so_7, so_8, so_9, sq_5, \
                         sq_6, sq_7, sq_8, sq_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * so_5[k]
                     + sq_5[k];

            t_6[k] = -ab_x[k] * so_6[k]
                     + sq_6[k];

            t_7[k] = -ab_x[k] * so_7[k]
                     + sq_7[k];

            t_8[k] = -ab_x[k] * so_8[k]
                     + sq_8[k];

            t_9[k] = -ab_x[k] * so_9[k]
                     + sq_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, so_10, so_11, so_12, so_13, \
                         so_14, sq_10, sq_11, sq_12, sq_13, sq_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * so_10[k]
                      + sq_10[k];

            t_11[k] = -ab_x[k] * so_11[k]
                      + sq_11[k];

            t_12[k] = -ab_x[k] * so_12[k]
                      + sq_12[k];

            t_13[k] = -ab_x[k] * so_13[k]
                      + sq_13[k];

            t_14[k] = -ab_x[k] * so_14[k]
                      + sq_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, so_15, so_16, so_17, so_18, \
                         so_19, sq_15, sq_16, sq_17, sq_18, sq_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * so_15[k]
                      + sq_15[k];

            t_16[k] = -ab_x[k] * so_16[k]
                      + sq_16[k];

            t_17[k] = -ab_x[k] * so_17[k]
                      + sq_17[k];

            t_18[k] = -ab_x[k] * so_18[k]
                      + sq_18[k];

            t_19[k] = -ab_x[k] * so_19[k]
                      + sq_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, so_20, so_21, so_22, so_23, \
                         so_24, sq_20, sq_21, sq_22, sq_23, sq_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * so_20[k]
                      + sq_20[k];

            t_21[k] = -ab_x[k] * so_21[k]
                      + sq_21[k];

            t_22[k] = -ab_x[k] * so_22[k]
                      + sq_22[k];

            t_23[k] = -ab_x[k] * so_23[k]
                      + sq_23[k];

            t_24[k] = -ab_x[k] * so_24[k]
                      + sq_24[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, so_25, so_26, so_27, so_28, \
                         so_29, sq_25, sq_26, sq_27, sq_28, sq_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * so_25[k]
                      + sq_25[k];

            t_26[k] = -ab_x[k] * so_26[k]
                      + sq_26[k];

            t_27[k] = -ab_x[k] * so_27[k]
                      + sq_27[k];

            t_28[k] = -ab_x[k] * so_28[k]
                      + sq_28[k];

            t_29[k] = -ab_x[k] * so_29[k]
                      + sq_29[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, so_30, so_31, so_32, so_33, \
                         so_34, sq_30, sq_31, sq_32, sq_33, sq_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * so_30[k]
                      + sq_30[k];

            t_31[k] = -ab_x[k] * so_31[k]
                      + sq_31[k];

            t_32[k] = -ab_x[k] * so_32[k]
                      + sq_32[k];

            t_33[k] = -ab_x[k] * so_33[k]
                      + sq_33[k];

            t_34[k] = -ab_x[k] * so_34[k]
                      + sq_34[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, so_35, so_36, so_37, so_38, \
                         so_39, sq_35, sq_36, sq_37, sq_38, sq_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * so_35[k]
                      + sq_35[k];

            t_36[k] = -ab_x[k] * so_36[k]
                      + sq_36[k];

            t_37[k] = -ab_x[k] * so_37[k]
                      + sq_37[k];

            t_38[k] = -ab_x[k] * so_38[k]
                      + sq_38[k];

            t_39[k] = -ab_x[k] * so_39[k]
                      + sq_39[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, so_40, so_41, so_42, so_43, \
                         so_44, sq_40, sq_41, sq_42, sq_43, sq_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * so_40[k]
                      + sq_40[k];

            t_41[k] = -ab_x[k] * so_41[k]
                      + sq_41[k];

            t_42[k] = -ab_x[k] * so_42[k]
                      + sq_42[k];

            t_43[k] = -ab_x[k] * so_43[k]
                      + sq_43[k];

            t_44[k] = -ab_x[k] * so_44[k]
                      + sq_44[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, so_45, so_46, so_47, so_48, \
                         so_49, sq_45, sq_46, sq_47, sq_48, sq_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * so_45[k]
                      + sq_45[k];

            t_46[k] = -ab_x[k] * so_46[k]
                      + sq_46[k];

            t_47[k] = -ab_x[k] * so_47[k]
                      + sq_47[k];

            t_48[k] = -ab_x[k] * so_48[k]
                      + sq_48[k];

            t_49[k] = -ab_x[k] * so_49[k]
                      + sq_49[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, so_50, so_51, so_52, so_53, \
                         so_54, sq_50, sq_51, sq_52, sq_53, sq_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * so_50[k]
                      + sq_50[k];

            t_51[k] = -ab_x[k] * so_51[k]
                      + sq_51[k];

            t_52[k] = -ab_x[k] * so_52[k]
                      + sq_52[k];

            t_53[k] = -ab_x[k] * so_53[k]
                      + sq_53[k];

            t_54[k] = -ab_x[k] * so_54[k]
                      + sq_54[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, so_55, so_56, so_57, so_58, \
                         so_59, sq_55, sq_56, sq_57, sq_58, sq_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * so_55[k]
                      + sq_55[k];

            t_56[k] = -ab_x[k] * so_56[k]
                      + sq_56[k];

            t_57[k] = -ab_x[k] * so_57[k]
                      + sq_57[k];

            t_58[k] = -ab_x[k] * so_58[k]
                      + sq_58[k];

            t_59[k] = -ab_x[k] * so_59[k]
                      + sq_59[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, so_60, so_61, so_62, so_63, \
                         so_64, sq_60, sq_61, sq_62, sq_63, sq_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * so_60[k]
                      + sq_60[k];

            t_61[k] = -ab_x[k] * so_61[k]
                      + sq_61[k];

            t_62[k] = -ab_x[k] * so_62[k]
                      + sq_62[k];

            t_63[k] = -ab_x[k] * so_63[k]
                      + sq_63[k];

            t_64[k] = -ab_x[k] * so_64[k]
                      + sq_64[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, so_65, so_66, so_67, so_68, \
                         so_69, sq_65, sq_66, sq_67, sq_68, sq_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * so_65[k]
                      + sq_65[k];

            t_66[k] = -ab_x[k] * so_66[k]
                      + sq_66[k];

            t_67[k] = -ab_x[k] * so_67[k]
                      + sq_67[k];

            t_68[k] = -ab_x[k] * so_68[k]
                      + sq_68[k];

            t_69[k] = -ab_x[k] * so_69[k]
                      + sq_69[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, so_70, so_71, so_72, so_73, \
                         so_74, sq_70, sq_71, sq_72, sq_73, sq_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * so_70[k]
                      + sq_70[k];

            t_71[k] = -ab_x[k] * so_71[k]
                      + sq_71[k];

            t_72[k] = -ab_x[k] * so_72[k]
                      + sq_72[k];

            t_73[k] = -ab_x[k] * so_73[k]
                      + sq_73[k];

            t_74[k] = -ab_x[k] * so_74[k]
                      + sq_74[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, ab_x, ab_y, so_0, so_75, so_76, so_77, sq_1, \
                         sq_75, sq_76, sq_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * so_75[k]
                      + sq_75[k];

            t_76[k] = -ab_x[k] * so_76[k]
                      + sq_76[k];

            t_77[k] = -ab_x[k] * so_77[k]
                      + sq_77[k];

            t_78[k] = -ab_y[k] * so_0[k]
                      + sq_1[k];
        }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, ab_y, so_1, so_2, so_3, so_4, so_5, \
                         sq_3, sq_4, sq_6, sq_7, sq_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_79[k] = -ab_y[k] * so_1[k]
                      + sq_3[k];

            t_80[k] = -ab_y[k] * so_2[k]
                      + sq_4[k];

            t_81[k] = -ab_y[k] * so_3[k]
                      + sq_6[k];

            t_82[k] = -ab_y[k] * so_4[k]
                      + sq_7[k];

            t_83[k] = -ab_y[k] * so_5[k]
                      + sq_8[k];
        }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ab_y, so_6, so_7, so_8, so_9, so_10, \
                         sq_10, sq_11, sq_12, sq_13, sq_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_84[k] = -ab_y[k] * so_6[k]
                      + sq_10[k];

            t_85[k] = -ab_y[k] * so_7[k]
                      + sq_11[k];

            t_86[k] = -ab_y[k] * so_8[k]
                      + sq_12[k];

            t_87[k] = -ab_y[k] * so_9[k]
                      + sq_13[k];

            t_88[k] = -ab_y[k] * so_10[k]
                      + sq_15[k];
        }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ab_y, so_11, so_12, so_13, so_14, \
                         so_15, sq_16, sq_17, sq_18, sq_19, sq_21 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_89[k] = -ab_y[k] * so_11[k]
                      + sq_16[k];

            t_90[k] = -ab_y[k] * so_12[k]
                      + sq_17[k];

            t_91[k] = -ab_y[k] * so_13[k]
                      + sq_18[k];

            t_92[k] = -ab_y[k] * so_14[k]
                      + sq_19[k];

            t_93[k] = -ab_y[k] * so_15[k]
                      + sq_21[k];
        }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ab_y, so_16, so_17, so_18, so_19, \
                         so_20, sq_22, sq_23, sq_24, sq_25, sq_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_94[k] = -ab_y[k] * so_16[k]
                      + sq_22[k];

            t_95[k] = -ab_y[k] * so_17[k]
                      + sq_23[k];

            t_96[k] = -ab_y[k] * so_18[k]
                      + sq_24[k];

            t_97[k] = -ab_y[k] * so_19[k]
                      + sq_25[k];

            t_98[k] = -ab_y[k] * so_20[k]
                      + sq_26[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ab_y, so_21, so_22, so_23, so_24, \
                         so_25, sq_28, sq_29, sq_30, sq_31, sq_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = -ab_y[k] * so_21[k]
                      + sq_28[k];

            t_100[k] = -ab_y[k] * so_22[k]
                       + sq_29[k];

            t_101[k] = -ab_y[k] * so_23[k]
                       + sq_30[k];

            t_102[k] = -ab_y[k] * so_24[k]
                       + sq_31[k];

            t_103[k] = -ab_y[k] * so_25[k]
                       + sq_32[k];
        }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, ab_y, so_26, so_27, so_28, so_29, \
                         so_30, sq_33, sq_34, sq_36, sq_37, sq_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_104[k] = -ab_y[k] * so_26[k]
                       + sq_33[k];

            t_105[k] = -ab_y[k] * so_27[k]
                       + sq_34[k];

            t_106[k] = -ab_y[k] * so_28[k]
                       + sq_36[k];

            t_107[k] = -ab_y[k] * so_29[k]
                       + sq_37[k];

            t_108[k] = -ab_y[k] * so_30[k]
                       + sq_38[k];
        }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, ab_y, so_31, so_32, so_33, so_34, \
                         so_35, sq_39, sq_40, sq_41, sq_42, sq_43 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_109[k] = -ab_y[k] * so_31[k]
                       + sq_39[k];

            t_110[k] = -ab_y[k] * so_32[k]
                       + sq_40[k];

            t_111[k] = -ab_y[k] * so_33[k]
                       + sq_41[k];

            t_112[k] = -ab_y[k] * so_34[k]
                       + sq_42[k];

            t_113[k] = -ab_y[k] * so_35[k]
                       + sq_43[k];
        }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, ab_y, so_36, so_37, so_38, so_39, \
                         so_40, sq_45, sq_46, sq_47, sq_48, sq_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_114[k] = -ab_y[k] * so_36[k]
                       + sq_45[k];

            t_115[k] = -ab_y[k] * so_37[k]
                       + sq_46[k];

            t_116[k] = -ab_y[k] * so_38[k]
                       + sq_47[k];

            t_117[k] = -ab_y[k] * so_39[k]
                       + sq_48[k];

            t_118[k] = -ab_y[k] * so_40[k]
                       + sq_49[k];
        }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, ab_y, so_41, so_42, so_43, so_44, \
                         so_45, sq_50, sq_51, sq_52, sq_53, sq_55 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_119[k] = -ab_y[k] * so_41[k]
                       + sq_50[k];

            t_120[k] = -ab_y[k] * so_42[k]
                       + sq_51[k];

            t_121[k] = -ab_y[k] * so_43[k]
                       + sq_52[k];

            t_122[k] = -ab_y[k] * so_44[k]
                       + sq_53[k];

            t_123[k] = -ab_y[k] * so_45[k]
                       + sq_55[k];
        }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, ab_y, so_46, so_47, so_48, so_49, \
                         so_50, sq_56, sq_57, sq_58, sq_59, sq_60 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_124[k] = -ab_y[k] * so_46[k]
                       + sq_56[k];

            t_125[k] = -ab_y[k] * so_47[k]
                       + sq_57[k];

            t_126[k] = -ab_y[k] * so_48[k]
                       + sq_58[k];

            t_127[k] = -ab_y[k] * so_49[k]
                       + sq_59[k];

            t_128[k] = -ab_y[k] * so_50[k]
                       + sq_60[k];
        }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, ab_y, so_51, so_52, so_53, so_54, \
                         so_55, sq_61, sq_62, sq_63, sq_64, sq_66 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_129[k] = -ab_y[k] * so_51[k]
                       + sq_61[k];

            t_130[k] = -ab_y[k] * so_52[k]
                       + sq_62[k];

            t_131[k] = -ab_y[k] * so_53[k]
                       + sq_63[k];

            t_132[k] = -ab_y[k] * so_54[k]
                       + sq_64[k];

            t_133[k] = -ab_y[k] * so_55[k]
                       + sq_66[k];
        }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ab_y, so_56, so_57, so_58, so_59, \
                         so_60, sq_67, sq_68, sq_69, sq_70, sq_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_134[k] = -ab_y[k] * so_56[k]
                       + sq_67[k];

            t_135[k] = -ab_y[k] * so_57[k]
                       + sq_68[k];

            t_136[k] = -ab_y[k] * so_58[k]
                       + sq_69[k];

            t_137[k] = -ab_y[k] * so_59[k]
                       + sq_70[k];

            t_138[k] = -ab_y[k] * so_60[k]
                       + sq_71[k];
        }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_y, so_61, so_62, so_63, so_64, \
                         so_65, sq_72, sq_73, sq_74, sq_75, sq_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_139[k] = -ab_y[k] * so_61[k]
                       + sq_72[k];

            t_140[k] = -ab_y[k] * so_62[k]
                       + sq_73[k];

            t_141[k] = -ab_y[k] * so_63[k]
                       + sq_74[k];

            t_142[k] = -ab_y[k] * so_64[k]
                       + sq_75[k];

            t_143[k] = -ab_y[k] * so_65[k]
                       + sq_76[k];
        }
    }
}

static auto
compute_hrr_po_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t so, const size_t sq, const size_t ncomps,
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

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *so_0 = buffer.data(so + 0 * ncomps + c);
        const auto *so_1 = buffer.data(so + 1 * ncomps + c);
        const auto *so_2 = buffer.data(so + 2 * ncomps + c);
        const auto *so_3 = buffer.data(so + 3 * ncomps + c);
        const auto *so_4 = buffer.data(so + 4 * ncomps + c);
        const auto *so_5 = buffer.data(so + 5 * ncomps + c);
        const auto *so_6 = buffer.data(so + 6 * ncomps + c);
        const auto *so_7 = buffer.data(so + 7 * ncomps + c);
        const auto *so_8 = buffer.data(so + 8 * ncomps + c);
        const auto *so_9 = buffer.data(so + 9 * ncomps + c);
        const auto *so_10 = buffer.data(so + 10 * ncomps + c);
        const auto *so_11 = buffer.data(so + 11 * ncomps + c);
        const auto *so_12 = buffer.data(so + 12 * ncomps + c);
        const auto *so_13 = buffer.data(so + 13 * ncomps + c);
        const auto *so_14 = buffer.data(so + 14 * ncomps + c);
        const auto *so_15 = buffer.data(so + 15 * ncomps + c);
        const auto *so_16 = buffer.data(so + 16 * ncomps + c);
        const auto *so_17 = buffer.data(so + 17 * ncomps + c);
        const auto *so_18 = buffer.data(so + 18 * ncomps + c);
        const auto *so_19 = buffer.data(so + 19 * ncomps + c);
        const auto *so_20 = buffer.data(so + 20 * ncomps + c);
        const auto *so_21 = buffer.data(so + 21 * ncomps + c);
        const auto *so_22 = buffer.data(so + 22 * ncomps + c);
        const auto *so_23 = buffer.data(so + 23 * ncomps + c);
        const auto *so_24 = buffer.data(so + 24 * ncomps + c);
        const auto *so_25 = buffer.data(so + 25 * ncomps + c);
        const auto *so_26 = buffer.data(so + 26 * ncomps + c);
        const auto *so_27 = buffer.data(so + 27 * ncomps + c);
        const auto *so_28 = buffer.data(so + 28 * ncomps + c);
        const auto *so_29 = buffer.data(so + 29 * ncomps + c);
        const auto *so_30 = buffer.data(so + 30 * ncomps + c);
        const auto *so_31 = buffer.data(so + 31 * ncomps + c);
        const auto *so_32 = buffer.data(so + 32 * ncomps + c);
        const auto *so_33 = buffer.data(so + 33 * ncomps + c);
        const auto *so_34 = buffer.data(so + 34 * ncomps + c);
        const auto *so_35 = buffer.data(so + 35 * ncomps + c);
        const auto *so_36 = buffer.data(so + 36 * ncomps + c);
        const auto *so_37 = buffer.data(so + 37 * ncomps + c);
        const auto *so_38 = buffer.data(so + 38 * ncomps + c);
        const auto *so_39 = buffer.data(so + 39 * ncomps + c);
        const auto *so_40 = buffer.data(so + 40 * ncomps + c);
        const auto *so_41 = buffer.data(so + 41 * ncomps + c);
        const auto *so_42 = buffer.data(so + 42 * ncomps + c);
        const auto *so_43 = buffer.data(so + 43 * ncomps + c);
        const auto *so_44 = buffer.data(so + 44 * ncomps + c);
        const auto *so_45 = buffer.data(so + 45 * ncomps + c);
        const auto *so_46 = buffer.data(so + 46 * ncomps + c);
        const auto *so_47 = buffer.data(so + 47 * ncomps + c);
        const auto *so_48 = buffer.data(so + 48 * ncomps + c);
        const auto *so_49 = buffer.data(so + 49 * ncomps + c);
        const auto *so_50 = buffer.data(so + 50 * ncomps + c);
        const auto *so_51 = buffer.data(so + 51 * ncomps + c);
        const auto *so_52 = buffer.data(so + 52 * ncomps + c);
        const auto *so_53 = buffer.data(so + 53 * ncomps + c);
        const auto *so_54 = buffer.data(so + 54 * ncomps + c);
        const auto *so_55 = buffer.data(so + 55 * ncomps + c);
        const auto *so_56 = buffer.data(so + 56 * ncomps + c);
        const auto *so_57 = buffer.data(so + 57 * ncomps + c);
        const auto *so_58 = buffer.data(so + 58 * ncomps + c);
        const auto *so_59 = buffer.data(so + 59 * ncomps + c);
        const auto *so_60 = buffer.data(so + 60 * ncomps + c);
        const auto *so_61 = buffer.data(so + 61 * ncomps + c);
        const auto *so_62 = buffer.data(so + 62 * ncomps + c);
        const auto *so_63 = buffer.data(so + 63 * ncomps + c);
        const auto *so_64 = buffer.data(so + 64 * ncomps + c);
        const auto *so_65 = buffer.data(so + 65 * ncomps + c);
        const auto *so_66 = buffer.data(so + 66 * ncomps + c);
        const auto *so_67 = buffer.data(so + 67 * ncomps + c);
        const auto *so_68 = buffer.data(so + 68 * ncomps + c);
        const auto *so_69 = buffer.data(so + 69 * ncomps + c);
        const auto *so_70 = buffer.data(so + 70 * ncomps + c);
        const auto *so_71 = buffer.data(so + 71 * ncomps + c);
        const auto *so_72 = buffer.data(so + 72 * ncomps + c);
        const auto *so_73 = buffer.data(so + 73 * ncomps + c);
        const auto *so_74 = buffer.data(so + 74 * ncomps + c);
        const auto *so_75 = buffer.data(so + 75 * ncomps + c);
        const auto *so_76 = buffer.data(so + 76 * ncomps + c);
        const auto *so_77 = buffer.data(so + 77 * ncomps + c);

        const auto *sq_2 = buffer.data(sq + 2 * ncomps + c);
        const auto *sq_4 = buffer.data(sq + 4 * ncomps + c);
        const auto *sq_5 = buffer.data(sq + 5 * ncomps + c);
        const auto *sq_7 = buffer.data(sq + 7 * ncomps + c);
        const auto *sq_8 = buffer.data(sq + 8 * ncomps + c);
        const auto *sq_9 = buffer.data(sq + 9 * ncomps + c);
        const auto *sq_11 = buffer.data(sq + 11 * ncomps + c);
        const auto *sq_12 = buffer.data(sq + 12 * ncomps + c);
        const auto *sq_13 = buffer.data(sq + 13 * ncomps + c);
        const auto *sq_14 = buffer.data(sq + 14 * ncomps + c);
        const auto *sq_16 = buffer.data(sq + 16 * ncomps + c);
        const auto *sq_17 = buffer.data(sq + 17 * ncomps + c);
        const auto *sq_18 = buffer.data(sq + 18 * ncomps + c);
        const auto *sq_19 = buffer.data(sq + 19 * ncomps + c);
        const auto *sq_20 = buffer.data(sq + 20 * ncomps + c);
        const auto *sq_22 = buffer.data(sq + 22 * ncomps + c);
        const auto *sq_23 = buffer.data(sq + 23 * ncomps + c);
        const auto *sq_24 = buffer.data(sq + 24 * ncomps + c);
        const auto *sq_25 = buffer.data(sq + 25 * ncomps + c);
        const auto *sq_26 = buffer.data(sq + 26 * ncomps + c);
        const auto *sq_27 = buffer.data(sq + 27 * ncomps + c);
        const auto *sq_29 = buffer.data(sq + 29 * ncomps + c);
        const auto *sq_30 = buffer.data(sq + 30 * ncomps + c);
        const auto *sq_31 = buffer.data(sq + 31 * ncomps + c);
        const auto *sq_32 = buffer.data(sq + 32 * ncomps + c);
        const auto *sq_33 = buffer.data(sq + 33 * ncomps + c);
        const auto *sq_34 = buffer.data(sq + 34 * ncomps + c);
        const auto *sq_35 = buffer.data(sq + 35 * ncomps + c);
        const auto *sq_37 = buffer.data(sq + 37 * ncomps + c);
        const auto *sq_38 = buffer.data(sq + 38 * ncomps + c);
        const auto *sq_39 = buffer.data(sq + 39 * ncomps + c);
        const auto *sq_40 = buffer.data(sq + 40 * ncomps + c);
        const auto *sq_41 = buffer.data(sq + 41 * ncomps + c);
        const auto *sq_42 = buffer.data(sq + 42 * ncomps + c);
        const auto *sq_43 = buffer.data(sq + 43 * ncomps + c);
        const auto *sq_44 = buffer.data(sq + 44 * ncomps + c);
        const auto *sq_46 = buffer.data(sq + 46 * ncomps + c);
        const auto *sq_47 = buffer.data(sq + 47 * ncomps + c);
        const auto *sq_48 = buffer.data(sq + 48 * ncomps + c);
        const auto *sq_49 = buffer.data(sq + 49 * ncomps + c);
        const auto *sq_50 = buffer.data(sq + 50 * ncomps + c);
        const auto *sq_51 = buffer.data(sq + 51 * ncomps + c);
        const auto *sq_52 = buffer.data(sq + 52 * ncomps + c);
        const auto *sq_53 = buffer.data(sq + 53 * ncomps + c);
        const auto *sq_54 = buffer.data(sq + 54 * ncomps + c);
        const auto *sq_56 = buffer.data(sq + 56 * ncomps + c);
        const auto *sq_57 = buffer.data(sq + 57 * ncomps + c);
        const auto *sq_58 = buffer.data(sq + 58 * ncomps + c);
        const auto *sq_59 = buffer.data(sq + 59 * ncomps + c);
        const auto *sq_60 = buffer.data(sq + 60 * ncomps + c);
        const auto *sq_61 = buffer.data(sq + 61 * ncomps + c);
        const auto *sq_62 = buffer.data(sq + 62 * ncomps + c);
        const auto *sq_63 = buffer.data(sq + 63 * ncomps + c);
        const auto *sq_64 = buffer.data(sq + 64 * ncomps + c);
        const auto *sq_65 = buffer.data(sq + 65 * ncomps + c);
        const auto *sq_67 = buffer.data(sq + 67 * ncomps + c);
        const auto *sq_68 = buffer.data(sq + 68 * ncomps + c);
        const auto *sq_69 = buffer.data(sq + 69 * ncomps + c);
        const auto *sq_70 = buffer.data(sq + 70 * ncomps + c);
        const auto *sq_71 = buffer.data(sq + 71 * ncomps + c);
        const auto *sq_72 = buffer.data(sq + 72 * ncomps + c);
        const auto *sq_73 = buffer.data(sq + 73 * ncomps + c);
        const auto *sq_74 = buffer.data(sq + 74 * ncomps + c);
        const auto *sq_75 = buffer.data(sq + 75 * ncomps + c);
        const auto *sq_76 = buffer.data(sq + 76 * ncomps + c);
        const auto *sq_77 = buffer.data(sq + 77 * ncomps + c);
        const auto *sq_78 = buffer.data(sq + 78 * ncomps + c);
        const auto *sq_79 = buffer.data(sq + 79 * ncomps + c);
        const auto *sq_80 = buffer.data(sq + 80 * ncomps + c);
        const auto *sq_81 = buffer.data(sq + 81 * ncomps + c);
        const auto *sq_82 = buffer.data(sq + 82 * ncomps + c);
        const auto *sq_83 = buffer.data(sq + 83 * ncomps + c);
        const auto *sq_84 = buffer.data(sq + 84 * ncomps + c);
        const auto *sq_85 = buffer.data(sq + 85 * ncomps + c);
        const auto *sq_86 = buffer.data(sq + 86 * ncomps + c);
        const auto *sq_87 = buffer.data(sq + 87 * ncomps + c);
        const auto *sq_88 = buffer.data(sq + 88 * ncomps + c);
        const auto *sq_89 = buffer.data(sq + 89 * ncomps + c);
        const auto *sq_90 = buffer.data(sq + 90 * ncomps + c);

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, ab_y, so_66, so_67, so_68, so_69, \
                         so_70, sq_78, sq_79, sq_80, sq_81, sq_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = -ab_y[k] * so_66[k]
                       + sq_78[k];

            t_145[k] = -ab_y[k] * so_67[k]
                       + sq_79[k];

            t_146[k] = -ab_y[k] * so_68[k]
                       + sq_80[k];

            t_147[k] = -ab_y[k] * so_69[k]
                       + sq_81[k];

            t_148[k] = -ab_y[k] * so_70[k]
                       + sq_82[k];
        }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, ab_y, so_71, so_72, so_73, so_74, \
                         so_75, sq_83, sq_84, sq_85, sq_86, sq_87 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_149[k] = -ab_y[k] * so_71[k]
                       + sq_83[k];

            t_150[k] = -ab_y[k] * so_72[k]
                       + sq_84[k];

            t_151[k] = -ab_y[k] * so_73[k]
                       + sq_85[k];

            t_152[k] = -ab_y[k] * so_74[k]
                       + sq_86[k];

            t_153[k] = -ab_y[k] * so_75[k]
                       + sq_87[k];
        }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, ab_y, ab_z, so_0, so_1, so_76, so_77, \
                         sq_2, sq_4, sq_88, sq_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_154[k] = -ab_y[k] * so_76[k]
                       + sq_88[k];

            t_155[k] = -ab_y[k] * so_77[k]
                       + sq_89[k];

            t_156[k] = -ab_z[k] * so_0[k]
                       + sq_2[k];

            t_157[k] = -ab_z[k] * so_1[k]
                       + sq_4[k];
        }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, ab_z, so_2, so_3, so_4, so_5, \
                         so_6, sq_5, sq_7, sq_8, sq_9, sq_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_158[k] = -ab_z[k] * so_2[k]
                       + sq_5[k];

            t_159[k] = -ab_z[k] * so_3[k]
                       + sq_7[k];

            t_160[k] = -ab_z[k] * so_4[k]
                       + sq_8[k];

            t_161[k] = -ab_z[k] * so_5[k]
                       + sq_9[k];

            t_162[k] = -ab_z[k] * so_6[k]
                       + sq_11[k];
        }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, ab_z, so_7, so_8, so_9, so_10, \
                         so_11, sq_12, sq_13, sq_14, sq_16, sq_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_163[k] = -ab_z[k] * so_7[k]
                       + sq_12[k];

            t_164[k] = -ab_z[k] * so_8[k]
                       + sq_13[k];

            t_165[k] = -ab_z[k] * so_9[k]
                       + sq_14[k];

            t_166[k] = -ab_z[k] * so_10[k]
                       + sq_16[k];

            t_167[k] = -ab_z[k] * so_11[k]
                       + sq_17[k];
        }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, ab_z, so_12, so_13, so_14, so_15, \
                         so_16, sq_18, sq_19, sq_20, sq_22, sq_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_168[k] = -ab_z[k] * so_12[k]
                       + sq_18[k];

            t_169[k] = -ab_z[k] * so_13[k]
                       + sq_19[k];

            t_170[k] = -ab_z[k] * so_14[k]
                       + sq_20[k];

            t_171[k] = -ab_z[k] * so_15[k]
                       + sq_22[k];

            t_172[k] = -ab_z[k] * so_16[k]
                       + sq_23[k];
        }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, ab_z, so_17, so_18, so_19, so_20, \
                         so_21, sq_24, sq_25, sq_26, sq_27, sq_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_173[k] = -ab_z[k] * so_17[k]
                       + sq_24[k];

            t_174[k] = -ab_z[k] * so_18[k]
                       + sq_25[k];

            t_175[k] = -ab_z[k] * so_19[k]
                       + sq_26[k];

            t_176[k] = -ab_z[k] * so_20[k]
                       + sq_27[k];

            t_177[k] = -ab_z[k] * so_21[k]
                       + sq_29[k];
        }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, ab_z, so_22, so_23, so_24, so_25, \
                         so_26, sq_30, sq_31, sq_32, sq_33, sq_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_178[k] = -ab_z[k] * so_22[k]
                       + sq_30[k];

            t_179[k] = -ab_z[k] * so_23[k]
                       + sq_31[k];

            t_180[k] = -ab_z[k] * so_24[k]
                       + sq_32[k];

            t_181[k] = -ab_z[k] * so_25[k]
                       + sq_33[k];

            t_182[k] = -ab_z[k] * so_26[k]
                       + sq_34[k];
        }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, ab_z, so_27, so_28, so_29, so_30, \
                         so_31, sq_35, sq_37, sq_38, sq_39, sq_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_183[k] = -ab_z[k] * so_27[k]
                       + sq_35[k];

            t_184[k] = -ab_z[k] * so_28[k]
                       + sq_37[k];

            t_185[k] = -ab_z[k] * so_29[k]
                       + sq_38[k];

            t_186[k] = -ab_z[k] * so_30[k]
                       + sq_39[k];

            t_187[k] = -ab_z[k] * so_31[k]
                       + sq_40[k];
        }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, ab_z, so_32, so_33, so_34, so_35, \
                         so_36, sq_41, sq_42, sq_43, sq_44, sq_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_188[k] = -ab_z[k] * so_32[k]
                       + sq_41[k];

            t_189[k] = -ab_z[k] * so_33[k]
                       + sq_42[k];

            t_190[k] = -ab_z[k] * so_34[k]
                       + sq_43[k];

            t_191[k] = -ab_z[k] * so_35[k]
                       + sq_44[k];

            t_192[k] = -ab_z[k] * so_36[k]
                       + sq_46[k];
        }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, ab_z, so_37, so_38, so_39, so_40, \
                         so_41, sq_47, sq_48, sq_49, sq_50, sq_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_193[k] = -ab_z[k] * so_37[k]
                       + sq_47[k];

            t_194[k] = -ab_z[k] * so_38[k]
                       + sq_48[k];

            t_195[k] = -ab_z[k] * so_39[k]
                       + sq_49[k];

            t_196[k] = -ab_z[k] * so_40[k]
                       + sq_50[k];

            t_197[k] = -ab_z[k] * so_41[k]
                       + sq_51[k];
        }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, ab_z, so_42, so_43, so_44, so_45, \
                         so_46, sq_52, sq_53, sq_54, sq_56, sq_57 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_198[k] = -ab_z[k] * so_42[k]
                       + sq_52[k];

            t_199[k] = -ab_z[k] * so_43[k]
                       + sq_53[k];

            t_200[k] = -ab_z[k] * so_44[k]
                       + sq_54[k];

            t_201[k] = -ab_z[k] * so_45[k]
                       + sq_56[k];

            t_202[k] = -ab_z[k] * so_46[k]
                       + sq_57[k];
        }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, ab_z, so_47, so_48, so_49, so_50, \
                         so_51, sq_58, sq_59, sq_60, sq_61, sq_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_203[k] = -ab_z[k] * so_47[k]
                       + sq_58[k];

            t_204[k] = -ab_z[k] * so_48[k]
                       + sq_59[k];

            t_205[k] = -ab_z[k] * so_49[k]
                       + sq_60[k];

            t_206[k] = -ab_z[k] * so_50[k]
                       + sq_61[k];

            t_207[k] = -ab_z[k] * so_51[k]
                       + sq_62[k];
        }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, ab_z, so_52, so_53, so_54, so_55, \
                         so_56, sq_63, sq_64, sq_65, sq_67, sq_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_208[k] = -ab_z[k] * so_52[k]
                       + sq_63[k];

            t_209[k] = -ab_z[k] * so_53[k]
                       + sq_64[k];

            t_210[k] = -ab_z[k] * so_54[k]
                       + sq_65[k];

            t_211[k] = -ab_z[k] * so_55[k]
                       + sq_67[k];

            t_212[k] = -ab_z[k] * so_56[k]
                       + sq_68[k];
        }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, ab_z, so_57, so_58, so_59, so_60, \
                         so_61, sq_69, sq_70, sq_71, sq_72, sq_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_213[k] = -ab_z[k] * so_57[k]
                       + sq_69[k];

            t_214[k] = -ab_z[k] * so_58[k]
                       + sq_70[k];

            t_215[k] = -ab_z[k] * so_59[k]
                       + sq_71[k];

            t_216[k] = -ab_z[k] * so_60[k]
                       + sq_72[k];

            t_217[k] = -ab_z[k] * so_61[k]
                       + sq_73[k];
        }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, ab_z, so_62, so_63, so_64, so_65, \
                         so_66, sq_74, sq_75, sq_76, sq_77, sq_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_218[k] = -ab_z[k] * so_62[k]
                       + sq_74[k];

            t_219[k] = -ab_z[k] * so_63[k]
                       + sq_75[k];

            t_220[k] = -ab_z[k] * so_64[k]
                       + sq_76[k];

            t_221[k] = -ab_z[k] * so_65[k]
                       + sq_77[k];

            t_222[k] = -ab_z[k] * so_66[k]
                       + sq_79[k];
        }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, ab_z, so_67, so_68, so_69, so_70, \
                         so_71, sq_80, sq_81, sq_82, sq_83, sq_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_223[k] = -ab_z[k] * so_67[k]
                       + sq_80[k];

            t_224[k] = -ab_z[k] * so_68[k]
                       + sq_81[k];

            t_225[k] = -ab_z[k] * so_69[k]
                       + sq_82[k];

            t_226[k] = -ab_z[k] * so_70[k]
                       + sq_83[k];

            t_227[k] = -ab_z[k] * so_71[k]
                       + sq_84[k];
        }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, ab_z, so_72, so_73, so_74, so_75, \
                         so_76, sq_85, sq_86, sq_87, sq_88, sq_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_228[k] = -ab_z[k] * so_72[k]
                       + sq_85[k];

            t_229[k] = -ab_z[k] * so_73[k]
                       + sq_86[k];

            t_230[k] = -ab_z[k] * so_74[k]
                       + sq_87[k];

            t_231[k] = -ab_z[k] * so_75[k]
                       + sq_88[k];

            t_232[k] = -ab_z[k] * so_76[k]
                       + sq_89[k];
        }

#pragma omp simd aligned(t_233, ab_z, so_77, sq_90 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_233[k] = -ab_z[k] * so_77[k]
                       + sq_90[k];
        }
    }
}

auto
compute_hrr_po(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t so, const size_t sq, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_po_piece0(buffer, coordinates, target, so, sq, ncomps, nmax);

    compute_hrr_po_piece1(buffer, coordinates, target, so, sq, ncomps, nmax);
}

}  // namespace simdtrf
