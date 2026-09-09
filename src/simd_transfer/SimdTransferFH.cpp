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


#include "SimdTransferFH.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_fh_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t dh, const size_t di, const size_t ncomps,
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
        const auto *dh_15 = buffer.data(dh + 15 * ncomps + c);
        const auto *dh_16 = buffer.data(dh + 16 * ncomps + c);
        const auto *dh_17 = buffer.data(dh + 17 * ncomps + c);
        const auto *dh_18 = buffer.data(dh + 18 * ncomps + c);
        const auto *dh_19 = buffer.data(dh + 19 * ncomps + c);
        const auto *dh_20 = buffer.data(dh + 20 * ncomps + c);
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
        const auto *dh_36 = buffer.data(dh + 36 * ncomps + c);
        const auto *dh_37 = buffer.data(dh + 37 * ncomps + c);
        const auto *dh_38 = buffer.data(dh + 38 * ncomps + c);
        const auto *dh_39 = buffer.data(dh + 39 * ncomps + c);
        const auto *dh_40 = buffer.data(dh + 40 * ncomps + c);
        const auto *dh_41 = buffer.data(dh + 41 * ncomps + c);
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
        const auto *dh_57 = buffer.data(dh + 57 * ncomps + c);
        const auto *dh_58 = buffer.data(dh + 58 * ncomps + c);
        const auto *dh_59 = buffer.data(dh + 59 * ncomps + c);
        const auto *dh_60 = buffer.data(dh + 60 * ncomps + c);
        const auto *dh_61 = buffer.data(dh + 61 * ncomps + c);
        const auto *dh_62 = buffer.data(dh + 62 * ncomps + c);
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
        const auto *dh_83 = buffer.data(dh + 83 * ncomps + c);
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
        const auto *dh_104 = buffer.data(dh + 104 * ncomps + c);
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
        const auto *dh_125 = buffer.data(dh + 125 * ncomps + c);

        const auto *di_0 = buffer.data(di + 0 * ncomps + c);
        const auto *di_1 = buffer.data(di + 1 * ncomps + c);
        const auto *di_2 = buffer.data(di + 2 * ncomps + c);
        const auto *di_3 = buffer.data(di + 3 * ncomps + c);
        const auto *di_4 = buffer.data(di + 4 * ncomps + c);
        const auto *di_5 = buffer.data(di + 5 * ncomps + c);
        const auto *di_6 = buffer.data(di + 6 * ncomps + c);
        const auto *di_7 = buffer.data(di + 7 * ncomps + c);
        const auto *di_8 = buffer.data(di + 8 * ncomps + c);
        const auto *di_9 = buffer.data(di + 9 * ncomps + c);
        const auto *di_10 = buffer.data(di + 10 * ncomps + c);
        const auto *di_11 = buffer.data(di + 11 * ncomps + c);
        const auto *di_12 = buffer.data(di + 12 * ncomps + c);
        const auto *di_13 = buffer.data(di + 13 * ncomps + c);
        const auto *di_14 = buffer.data(di + 14 * ncomps + c);
        const auto *di_15 = buffer.data(di + 15 * ncomps + c);
        const auto *di_16 = buffer.data(di + 16 * ncomps + c);
        const auto *di_17 = buffer.data(di + 17 * ncomps + c);
        const auto *di_18 = buffer.data(di + 18 * ncomps + c);
        const auto *di_19 = buffer.data(di + 19 * ncomps + c);
        const auto *di_20 = buffer.data(di + 20 * ncomps + c);
        const auto *di_28 = buffer.data(di + 28 * ncomps + c);
        const auto *di_29 = buffer.data(di + 29 * ncomps + c);
        const auto *di_30 = buffer.data(di + 30 * ncomps + c);
        const auto *di_31 = buffer.data(di + 31 * ncomps + c);
        const auto *di_32 = buffer.data(di + 32 * ncomps + c);
        const auto *di_33 = buffer.data(di + 33 * ncomps + c);
        const auto *di_34 = buffer.data(di + 34 * ncomps + c);
        const auto *di_35 = buffer.data(di + 35 * ncomps + c);
        const auto *di_36 = buffer.data(di + 36 * ncomps + c);
        const auto *di_37 = buffer.data(di + 37 * ncomps + c);
        const auto *di_38 = buffer.data(di + 38 * ncomps + c);
        const auto *di_39 = buffer.data(di + 39 * ncomps + c);
        const auto *di_40 = buffer.data(di + 40 * ncomps + c);
        const auto *di_41 = buffer.data(di + 41 * ncomps + c);
        const auto *di_42 = buffer.data(di + 42 * ncomps + c);
        const auto *di_43 = buffer.data(di + 43 * ncomps + c);
        const auto *di_44 = buffer.data(di + 44 * ncomps + c);
        const auto *di_45 = buffer.data(di + 45 * ncomps + c);
        const auto *di_46 = buffer.data(di + 46 * ncomps + c);
        const auto *di_47 = buffer.data(di + 47 * ncomps + c);
        const auto *di_48 = buffer.data(di + 48 * ncomps + c);
        const auto *di_56 = buffer.data(di + 56 * ncomps + c);
        const auto *di_57 = buffer.data(di + 57 * ncomps + c);
        const auto *di_58 = buffer.data(di + 58 * ncomps + c);
        const auto *di_59 = buffer.data(di + 59 * ncomps + c);
        const auto *di_60 = buffer.data(di + 60 * ncomps + c);
        const auto *di_61 = buffer.data(di + 61 * ncomps + c);
        const auto *di_62 = buffer.data(di + 62 * ncomps + c);
        const auto *di_63 = buffer.data(di + 63 * ncomps + c);
        const auto *di_64 = buffer.data(di + 64 * ncomps + c);
        const auto *di_65 = buffer.data(di + 65 * ncomps + c);
        const auto *di_66 = buffer.data(di + 66 * ncomps + c);
        const auto *di_67 = buffer.data(di + 67 * ncomps + c);
        const auto *di_68 = buffer.data(di + 68 * ncomps + c);
        const auto *di_69 = buffer.data(di + 69 * ncomps + c);
        const auto *di_70 = buffer.data(di + 70 * ncomps + c);
        const auto *di_71 = buffer.data(di + 71 * ncomps + c);
        const auto *di_72 = buffer.data(di + 72 * ncomps + c);
        const auto *di_73 = buffer.data(di + 73 * ncomps + c);
        const auto *di_74 = buffer.data(di + 74 * ncomps + c);
        const auto *di_75 = buffer.data(di + 75 * ncomps + c);
        const auto *di_76 = buffer.data(di + 76 * ncomps + c);
        const auto *di_84 = buffer.data(di + 84 * ncomps + c);
        const auto *di_85 = buffer.data(di + 85 * ncomps + c);
        const auto *di_86 = buffer.data(di + 86 * ncomps + c);
        const auto *di_87 = buffer.data(di + 87 * ncomps + c);
        const auto *di_88 = buffer.data(di + 88 * ncomps + c);
        const auto *di_89 = buffer.data(di + 89 * ncomps + c);
        const auto *di_90 = buffer.data(di + 90 * ncomps + c);
        const auto *di_91 = buffer.data(di + 91 * ncomps + c);
        const auto *di_92 = buffer.data(di + 92 * ncomps + c);
        const auto *di_93 = buffer.data(di + 93 * ncomps + c);
        const auto *di_94 = buffer.data(di + 94 * ncomps + c);
        const auto *di_95 = buffer.data(di + 95 * ncomps + c);
        const auto *di_96 = buffer.data(di + 96 * ncomps + c);
        const auto *di_97 = buffer.data(di + 97 * ncomps + c);
        const auto *di_98 = buffer.data(di + 98 * ncomps + c);
        const auto *di_99 = buffer.data(di + 99 * ncomps + c);
        const auto *di_100 = buffer.data(di + 100 * ncomps + c);
        const auto *di_101 = buffer.data(di + 101 * ncomps + c);
        const auto *di_102 = buffer.data(di + 102 * ncomps + c);
        const auto *di_103 = buffer.data(di + 103 * ncomps + c);
        const auto *di_104 = buffer.data(di + 104 * ncomps + c);
        const auto *di_105 = buffer.data(di + 105 * ncomps + c);
        const auto *di_106 = buffer.data(di + 106 * ncomps + c);
        const auto *di_107 = buffer.data(di + 107 * ncomps + c);
        const auto *di_112 = buffer.data(di + 112 * ncomps + c);
        const auto *di_113 = buffer.data(di + 113 * ncomps + c);
        const auto *di_114 = buffer.data(di + 114 * ncomps + c);
        const auto *di_115 = buffer.data(di + 115 * ncomps + c);
        const auto *di_116 = buffer.data(di + 116 * ncomps + c);
        const auto *di_117 = buffer.data(di + 117 * ncomps + c);
        const auto *di_118 = buffer.data(di + 118 * ncomps + c);
        const auto *di_119 = buffer.data(di + 119 * ncomps + c);
        const auto *di_120 = buffer.data(di + 120 * ncomps + c);
        const auto *di_121 = buffer.data(di + 121 * ncomps + c);
        const auto *di_122 = buffer.data(di + 122 * ncomps + c);
        const auto *di_123 = buffer.data(di + 123 * ncomps + c);
        const auto *di_124 = buffer.data(di + 124 * ncomps + c);
        const auto *di_125 = buffer.data(di + 125 * ncomps + c);
        const auto *di_126 = buffer.data(di + 126 * ncomps + c);
        const auto *di_127 = buffer.data(di + 127 * ncomps + c);
        const auto *di_128 = buffer.data(di + 128 * ncomps + c);
        const auto *di_129 = buffer.data(di + 129 * ncomps + c);
        const auto *di_130 = buffer.data(di + 130 * ncomps + c);
        const auto *di_131 = buffer.data(di + 131 * ncomps + c);
        const auto *di_132 = buffer.data(di + 132 * ncomps + c);
        const auto *di_140 = buffer.data(di + 140 * ncomps + c);
        const auto *di_141 = buffer.data(di + 141 * ncomps + c);
        const auto *di_142 = buffer.data(di + 142 * ncomps + c);
        const auto *di_143 = buffer.data(di + 143 * ncomps + c);
        const auto *di_144 = buffer.data(di + 144 * ncomps + c);
        const auto *di_145 = buffer.data(di + 145 * ncomps + c);
        const auto *di_146 = buffer.data(di + 146 * ncomps + c);
        const auto *di_147 = buffer.data(di + 147 * ncomps + c);
        const auto *di_148 = buffer.data(di + 148 * ncomps + c);
        const auto *di_149 = buffer.data(di + 149 * ncomps + c);
        const auto *di_150 = buffer.data(di + 150 * ncomps + c);
        const auto *di_151 = buffer.data(di + 151 * ncomps + c);
        const auto *di_152 = buffer.data(di + 152 * ncomps + c);
        const auto *di_153 = buffer.data(di + 153 * ncomps + c);
        const auto *di_154 = buffer.data(di + 154 * ncomps + c);
        const auto *di_155 = buffer.data(di + 155 * ncomps + c);
        const auto *di_156 = buffer.data(di + 156 * ncomps + c);
        const auto *di_157 = buffer.data(di + 157 * ncomps + c);
        const auto *di_158 = buffer.data(di + 158 * ncomps + c);
        const auto *di_159 = buffer.data(di + 159 * ncomps + c);
        const auto *di_160 = buffer.data(di + 160 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, dh_0, dh_1, dh_2, dh_3, dh_4, di_0, \
                         di_1, di_2, di_3, di_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * dh_0[k]
                     + di_0[k];

            t_1[k] = -ab_x[k] * dh_1[k]
                     + di_1[k];

            t_2[k] = -ab_x[k] * dh_2[k]
                     + di_2[k];

            t_3[k] = -ab_x[k] * dh_3[k]
                     + di_3[k];

            t_4[k] = -ab_x[k] * dh_4[k]
                     + di_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, dh_5, dh_6, dh_7, dh_8, dh_9, di_5, \
                         di_6, di_7, di_8, di_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * dh_5[k]
                     + di_5[k];

            t_6[k] = -ab_x[k] * dh_6[k]
                     + di_6[k];

            t_7[k] = -ab_x[k] * dh_7[k]
                     + di_7[k];

            t_8[k] = -ab_x[k] * dh_8[k]
                     + di_8[k];

            t_9[k] = -ab_x[k] * dh_9[k]
                     + di_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, dh_10, dh_11, dh_12, dh_13, \
                         dh_14, di_10, di_11, di_12, di_13, di_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * dh_10[k]
                      + di_10[k];

            t_11[k] = -ab_x[k] * dh_11[k]
                      + di_11[k];

            t_12[k] = -ab_x[k] * dh_12[k]
                      + di_12[k];

            t_13[k] = -ab_x[k] * dh_13[k]
                      + di_13[k];

            t_14[k] = -ab_x[k] * dh_14[k]
                      + di_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, dh_15, dh_16, dh_17, dh_18, \
                         dh_19, di_15, di_16, di_17, di_18, di_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * dh_15[k]
                      + di_15[k];

            t_16[k] = -ab_x[k] * dh_16[k]
                      + di_16[k];

            t_17[k] = -ab_x[k] * dh_17[k]
                      + di_17[k];

            t_18[k] = -ab_x[k] * dh_18[k]
                      + di_18[k];

            t_19[k] = -ab_x[k] * dh_19[k]
                      + di_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, dh_20, dh_21, dh_22, dh_23, \
                         dh_24, di_20, di_28, di_29, di_30, di_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * dh_20[k]
                      + di_20[k];

            t_21[k] = -ab_x[k] * dh_21[k]
                      + di_28[k];

            t_22[k] = -ab_x[k] * dh_22[k]
                      + di_29[k];

            t_23[k] = -ab_x[k] * dh_23[k]
                      + di_30[k];

            t_24[k] = -ab_x[k] * dh_24[k]
                      + di_31[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, dh_25, dh_26, dh_27, dh_28, \
                         dh_29, di_32, di_33, di_34, di_35, di_36 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * dh_25[k]
                      + di_32[k];

            t_26[k] = -ab_x[k] * dh_26[k]
                      + di_33[k];

            t_27[k] = -ab_x[k] * dh_27[k]
                      + di_34[k];

            t_28[k] = -ab_x[k] * dh_28[k]
                      + di_35[k];

            t_29[k] = -ab_x[k] * dh_29[k]
                      + di_36[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, dh_30, dh_31, dh_32, dh_33, \
                         dh_34, di_37, di_38, di_39, di_40, di_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * dh_30[k]
                      + di_37[k];

            t_31[k] = -ab_x[k] * dh_31[k]
                      + di_38[k];

            t_32[k] = -ab_x[k] * dh_32[k]
                      + di_39[k];

            t_33[k] = -ab_x[k] * dh_33[k]
                      + di_40[k];

            t_34[k] = -ab_x[k] * dh_34[k]
                      + di_41[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, dh_35, dh_36, dh_37, dh_38, \
                         dh_39, di_42, di_43, di_44, di_45, di_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * dh_35[k]
                      + di_42[k];

            t_36[k] = -ab_x[k] * dh_36[k]
                      + di_43[k];

            t_37[k] = -ab_x[k] * dh_37[k]
                      + di_44[k];

            t_38[k] = -ab_x[k] * dh_38[k]
                      + di_45[k];

            t_39[k] = -ab_x[k] * dh_39[k]
                      + di_46[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, dh_40, dh_41, dh_42, dh_43, \
                         dh_44, di_47, di_48, di_56, di_57, di_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * dh_40[k]
                      + di_47[k];

            t_41[k] = -ab_x[k] * dh_41[k]
                      + di_48[k];

            t_42[k] = -ab_x[k] * dh_42[k]
                      + di_56[k];

            t_43[k] = -ab_x[k] * dh_43[k]
                      + di_57[k];

            t_44[k] = -ab_x[k] * dh_44[k]
                      + di_58[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, dh_45, dh_46, dh_47, dh_48, \
                         dh_49, di_59, di_60, di_61, di_62, di_63 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * dh_45[k]
                      + di_59[k];

            t_46[k] = -ab_x[k] * dh_46[k]
                      + di_60[k];

            t_47[k] = -ab_x[k] * dh_47[k]
                      + di_61[k];

            t_48[k] = -ab_x[k] * dh_48[k]
                      + di_62[k];

            t_49[k] = -ab_x[k] * dh_49[k]
                      + di_63[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, dh_50, dh_51, dh_52, dh_53, \
                         dh_54, di_64, di_65, di_66, di_67, di_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * dh_50[k]
                      + di_64[k];

            t_51[k] = -ab_x[k] * dh_51[k]
                      + di_65[k];

            t_52[k] = -ab_x[k] * dh_52[k]
                      + di_66[k];

            t_53[k] = -ab_x[k] * dh_53[k]
                      + di_67[k];

            t_54[k] = -ab_x[k] * dh_54[k]
                      + di_68[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, dh_55, dh_56, dh_57, dh_58, \
                         dh_59, di_69, di_70, di_71, di_72, di_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * dh_55[k]
                      + di_69[k];

            t_56[k] = -ab_x[k] * dh_56[k]
                      + di_70[k];

            t_57[k] = -ab_x[k] * dh_57[k]
                      + di_71[k];

            t_58[k] = -ab_x[k] * dh_58[k]
                      + di_72[k];

            t_59[k] = -ab_x[k] * dh_59[k]
                      + di_73[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, dh_60, dh_61, dh_62, dh_63, \
                         dh_64, di_74, di_75, di_76, di_84, di_85 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * dh_60[k]
                      + di_74[k];

            t_61[k] = -ab_x[k] * dh_61[k]
                      + di_75[k];

            t_62[k] = -ab_x[k] * dh_62[k]
                      + di_76[k];

            t_63[k] = -ab_x[k] * dh_63[k]
                      + di_84[k];

            t_64[k] = -ab_x[k] * dh_64[k]
                      + di_85[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, dh_65, dh_66, dh_67, dh_68, \
                         dh_69, di_86, di_87, di_88, di_89, di_90 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * dh_65[k]
                      + di_86[k];

            t_66[k] = -ab_x[k] * dh_66[k]
                      + di_87[k];

            t_67[k] = -ab_x[k] * dh_67[k]
                      + di_88[k];

            t_68[k] = -ab_x[k] * dh_68[k]
                      + di_89[k];

            t_69[k] = -ab_x[k] * dh_69[k]
                      + di_90[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, dh_70, dh_71, dh_72, dh_73, \
                         dh_74, di_91, di_92, di_93, di_94, di_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * dh_70[k]
                      + di_91[k];

            t_71[k] = -ab_x[k] * dh_71[k]
                      + di_92[k];

            t_72[k] = -ab_x[k] * dh_72[k]
                      + di_93[k];

            t_73[k] = -ab_x[k] * dh_73[k]
                      + di_94[k];

            t_74[k] = -ab_x[k] * dh_74[k]
                      + di_95[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, dh_75, dh_76, dh_77, dh_78, \
                         dh_79, di_96, di_97, di_98, di_99, di_100 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * dh_75[k]
                      + di_96[k];

            t_76[k] = -ab_x[k] * dh_76[k]
                      + di_97[k];

            t_77[k] = -ab_x[k] * dh_77[k]
                      + di_98[k];

            t_78[k] = -ab_x[k] * dh_78[k]
                      + di_99[k];

            t_79[k] = -ab_x[k] * dh_79[k]
                      + di_100[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, dh_80, dh_81, dh_82, dh_83, \
                         dh_84, di_101, di_102, di_103, di_104, \
                         di_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * dh_80[k]
                      + di_101[k];

            t_81[k] = -ab_x[k] * dh_81[k]
                      + di_102[k];

            t_82[k] = -ab_x[k] * dh_82[k]
                      + di_103[k];

            t_83[k] = -ab_x[k] * dh_83[k]
                      + di_104[k];

            t_84[k] = -ab_x[k] * dh_84[k]
                      + di_112[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, dh_85, dh_86, dh_87, dh_88, \
                         dh_89, di_113, di_114, di_115, di_116, \
                         di_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * dh_85[k]
                      + di_113[k];

            t_86[k] = -ab_x[k] * dh_86[k]
                      + di_114[k];

            t_87[k] = -ab_x[k] * dh_87[k]
                      + di_115[k];

            t_88[k] = -ab_x[k] * dh_88[k]
                      + di_116[k];

            t_89[k] = -ab_x[k] * dh_89[k]
                      + di_117[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, dh_90, dh_91, dh_92, dh_93, \
                         dh_94, di_118, di_119, di_120, di_121, \
                         di_122 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * dh_90[k]
                      + di_118[k];

            t_91[k] = -ab_x[k] * dh_91[k]
                      + di_119[k];

            t_92[k] = -ab_x[k] * dh_92[k]
                      + di_120[k];

            t_93[k] = -ab_x[k] * dh_93[k]
                      + di_121[k];

            t_94[k] = -ab_x[k] * dh_94[k]
                      + di_122[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, dh_95, dh_96, dh_97, dh_98, \
                         dh_99, di_123, di_124, di_125, di_126, \
                         di_127 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * dh_95[k]
                      + di_123[k];

            t_96[k] = -ab_x[k] * dh_96[k]
                      + di_124[k];

            t_97[k] = -ab_x[k] * dh_97[k]
                      + di_125[k];

            t_98[k] = -ab_x[k] * dh_98[k]
                      + di_126[k];

            t_99[k] = -ab_x[k] * dh_99[k]
                      + di_127[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, dh_100, dh_101, dh_102, \
                         dh_103, dh_104, di_128, di_129, di_130, di_131, \
                         di_132 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * dh_100[k]
                       + di_128[k];

            t_101[k] = -ab_x[k] * dh_101[k]
                       + di_129[k];

            t_102[k] = -ab_x[k] * dh_102[k]
                       + di_130[k];

            t_103[k] = -ab_x[k] * dh_103[k]
                       + di_131[k];

            t_104[k] = -ab_x[k] * dh_104[k]
                       + di_132[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, dh_105, dh_106, dh_107, \
                         dh_108, dh_109, di_140, di_141, di_142, di_143, \
                         di_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * dh_105[k]
                       + di_140[k];

            t_106[k] = -ab_x[k] * dh_106[k]
                       + di_141[k];

            t_107[k] = -ab_x[k] * dh_107[k]
                       + di_142[k];

            t_108[k] = -ab_x[k] * dh_108[k]
                       + di_143[k];

            t_109[k] = -ab_x[k] * dh_109[k]
                       + di_144[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, dh_110, dh_111, dh_112, \
                         dh_113, dh_114, di_145, di_146, di_147, di_148, \
                         di_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * dh_110[k]
                       + di_145[k];

            t_111[k] = -ab_x[k] * dh_111[k]
                       + di_146[k];

            t_112[k] = -ab_x[k] * dh_112[k]
                       + di_147[k];

            t_113[k] = -ab_x[k] * dh_113[k]
                       + di_148[k];

            t_114[k] = -ab_x[k] * dh_114[k]
                       + di_149[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, dh_115, dh_116, dh_117, \
                         dh_118, dh_119, di_150, di_151, di_152, di_153, \
                         di_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * dh_115[k]
                       + di_150[k];

            t_116[k] = -ab_x[k] * dh_116[k]
                       + di_151[k];

            t_117[k] = -ab_x[k] * dh_117[k]
                       + di_152[k];

            t_118[k] = -ab_x[k] * dh_118[k]
                       + di_153[k];

            t_119[k] = -ab_x[k] * dh_119[k]
                       + di_154[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, dh_120, dh_121, dh_122, \
                         dh_123, dh_124, di_155, di_156, di_157, di_158, \
                         di_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * dh_120[k]
                       + di_155[k];

            t_121[k] = -ab_x[k] * dh_121[k]
                       + di_156[k];

            t_122[k] = -ab_x[k] * dh_122[k]
                       + di_157[k];

            t_123[k] = -ab_x[k] * dh_123[k]
                       + di_158[k];

            t_124[k] = -ab_x[k] * dh_124[k]
                       + di_159[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, ab_x, ab_y, dh_63, dh_64, dh_65, dh_125, \
                         di_85, di_87, di_88, di_160 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * dh_125[k]
                       + di_160[k];

            t_126[k] = -ab_y[k] * dh_63[k]
                       + di_85[k];

            t_127[k] = -ab_y[k] * dh_64[k]
                       + di_87[k];

            t_128[k] = -ab_y[k] * dh_65[k]
                       + di_88[k];
        }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, ab_y, dh_66, dh_67, dh_68, dh_69, \
                         dh_70, di_90, di_91, di_92, di_94, di_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_129[k] = -ab_y[k] * dh_66[k]
                       + di_90[k];

            t_130[k] = -ab_y[k] * dh_67[k]
                       + di_91[k];

            t_131[k] = -ab_y[k] * dh_68[k]
                       + di_92[k];

            t_132[k] = -ab_y[k] * dh_69[k]
                       + di_94[k];

            t_133[k] = -ab_y[k] * dh_70[k]
                       + di_95[k];
        }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ab_y, dh_71, dh_72, dh_73, dh_74, \
                         dh_75, di_96, di_97, di_99, di_100, di_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_134[k] = -ab_y[k] * dh_71[k]
                       + di_96[k];

            t_135[k] = -ab_y[k] * dh_72[k]
                       + di_97[k];

            t_136[k] = -ab_y[k] * dh_73[k]
                       + di_99[k];

            t_137[k] = -ab_y[k] * dh_74[k]
                       + di_100[k];

            t_138[k] = -ab_y[k] * dh_75[k]
                       + di_101[k];
        }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_y, dh_76, dh_77, dh_78, dh_79, \
                         dh_80, di_102, di_103, di_105, di_106, \
                         di_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_139[k] = -ab_y[k] * dh_76[k]
                       + di_102[k];

            t_140[k] = -ab_y[k] * dh_77[k]
                       + di_103[k];

            t_141[k] = -ab_y[k] * dh_78[k]
                       + di_105[k];

            t_142[k] = -ab_y[k] * dh_79[k]
                       + di_106[k];

            t_143[k] = -ab_y[k] * dh_80[k]
                       + di_107[k];
        }
    }
}

static auto
compute_hrr_fh_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t dh, const size_t di, const size_t ncomps,
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

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *dh_81 = buffer.data(dh + 81 * ncomps + c);
        const auto *dh_82 = buffer.data(dh + 82 * ncomps + c);
        const auto *dh_83 = buffer.data(dh + 83 * ncomps + c);
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
        const auto *dh_104 = buffer.data(dh + 104 * ncomps + c);
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
        const auto *dh_125 = buffer.data(dh + 125 * ncomps + c);

        const auto *di_108 = buffer.data(di + 108 * ncomps + c);
        const auto *di_109 = buffer.data(di + 109 * ncomps + c);
        const auto *di_110 = buffer.data(di + 110 * ncomps + c);
        const auto *di_113 = buffer.data(di + 113 * ncomps + c);
        const auto *di_115 = buffer.data(di + 115 * ncomps + c);
        const auto *di_116 = buffer.data(di + 116 * ncomps + c);
        const auto *di_118 = buffer.data(di + 118 * ncomps + c);
        const auto *di_119 = buffer.data(di + 119 * ncomps + c);
        const auto *di_120 = buffer.data(di + 120 * ncomps + c);
        const auto *di_122 = buffer.data(di + 122 * ncomps + c);
        const auto *di_123 = buffer.data(di + 123 * ncomps + c);
        const auto *di_124 = buffer.data(di + 124 * ncomps + c);
        const auto *di_125 = buffer.data(di + 125 * ncomps + c);
        const auto *di_127 = buffer.data(di + 127 * ncomps + c);
        const auto *di_128 = buffer.data(di + 128 * ncomps + c);
        const auto *di_129 = buffer.data(di + 129 * ncomps + c);
        const auto *di_130 = buffer.data(di + 130 * ncomps + c);
        const auto *di_131 = buffer.data(di + 131 * ncomps + c);
        const auto *di_133 = buffer.data(di + 133 * ncomps + c);
        const auto *di_134 = buffer.data(di + 134 * ncomps + c);
        const auto *di_135 = buffer.data(di + 135 * ncomps + c);
        const auto *di_136 = buffer.data(di + 136 * ncomps + c);
        const auto *di_137 = buffer.data(di + 137 * ncomps + c);
        const auto *di_138 = buffer.data(di + 138 * ncomps + c);
        const auto *di_141 = buffer.data(di + 141 * ncomps + c);
        const auto *di_142 = buffer.data(di + 142 * ncomps + c);
        const auto *di_143 = buffer.data(di + 143 * ncomps + c);
        const auto *di_144 = buffer.data(di + 144 * ncomps + c);
        const auto *di_145 = buffer.data(di + 145 * ncomps + c);
        const auto *di_146 = buffer.data(di + 146 * ncomps + c);
        const auto *di_147 = buffer.data(di + 147 * ncomps + c);
        const auto *di_148 = buffer.data(di + 148 * ncomps + c);
        const auto *di_149 = buffer.data(di + 149 * ncomps + c);
        const auto *di_150 = buffer.data(di + 150 * ncomps + c);
        const auto *di_151 = buffer.data(di + 151 * ncomps + c);
        const auto *di_152 = buffer.data(di + 152 * ncomps + c);
        const auto *di_153 = buffer.data(di + 153 * ncomps + c);
        const auto *di_154 = buffer.data(di + 154 * ncomps + c);
        const auto *di_155 = buffer.data(di + 155 * ncomps + c);
        const auto *di_156 = buffer.data(di + 156 * ncomps + c);
        const auto *di_157 = buffer.data(di + 157 * ncomps + c);
        const auto *di_158 = buffer.data(di + 158 * ncomps + c);
        const auto *di_159 = buffer.data(di + 159 * ncomps + c);
        const auto *di_160 = buffer.data(di + 160 * ncomps + c);
        const auto *di_161 = buffer.data(di + 161 * ncomps + c);
        const auto *di_162 = buffer.data(di + 162 * ncomps + c);
        const auto *di_163 = buffer.data(di + 163 * ncomps + c);
        const auto *di_164 = buffer.data(di + 164 * ncomps + c);
        const auto *di_165 = buffer.data(di + 165 * ncomps + c);
        const auto *di_166 = buffer.data(di + 166 * ncomps + c);
        const auto *di_167 = buffer.data(di + 167 * ncomps + c);

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, ab_y, dh_81, dh_82, dh_83, dh_84, \
                         dh_85, di_108, di_109, di_110, di_113, \
                         di_115 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = -ab_y[k] * dh_81[k]
                       + di_108[k];

            t_145[k] = -ab_y[k] * dh_82[k]
                       + di_109[k];

            t_146[k] = -ab_y[k] * dh_83[k]
                       + di_110[k];

            t_147[k] = -ab_y[k] * dh_84[k]
                       + di_113[k];

            t_148[k] = -ab_y[k] * dh_85[k]
                       + di_115[k];
        }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, ab_y, dh_86, dh_87, dh_88, dh_89, \
                         dh_90, di_116, di_118, di_119, di_120, \
                         di_122 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_149[k] = -ab_y[k] * dh_86[k]
                       + di_116[k];

            t_150[k] = -ab_y[k] * dh_87[k]
                       + di_118[k];

            t_151[k] = -ab_y[k] * dh_88[k]
                       + di_119[k];

            t_152[k] = -ab_y[k] * dh_89[k]
                       + di_120[k];

            t_153[k] = -ab_y[k] * dh_90[k]
                       + di_122[k];
        }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, ab_y, dh_91, dh_92, dh_93, dh_94, \
                         dh_95, di_123, di_124, di_125, di_127, \
                         di_128 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_154[k] = -ab_y[k] * dh_91[k]
                       + di_123[k];

            t_155[k] = -ab_y[k] * dh_92[k]
                       + di_124[k];

            t_156[k] = -ab_y[k] * dh_93[k]
                       + di_125[k];

            t_157[k] = -ab_y[k] * dh_94[k]
                       + di_127[k];

            t_158[k] = -ab_y[k] * dh_95[k]
                       + di_128[k];
        }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, ab_y, dh_96, dh_97, dh_98, dh_99, \
                         dh_100, di_129, di_130, di_131, di_133, \
                         di_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_159[k] = -ab_y[k] * dh_96[k]
                       + di_129[k];

            t_160[k] = -ab_y[k] * dh_97[k]
                       + di_130[k];

            t_161[k] = -ab_y[k] * dh_98[k]
                       + di_131[k];

            t_162[k] = -ab_y[k] * dh_99[k]
                       + di_133[k];

            t_163[k] = -ab_y[k] * dh_100[k]
                       + di_134[k];
        }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, ab_y, dh_101, dh_102, dh_103, \
                         dh_104, dh_105, di_135, di_136, di_137, di_138, \
                         di_141 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_164[k] = -ab_y[k] * dh_101[k]
                       + di_135[k];

            t_165[k] = -ab_y[k] * dh_102[k]
                       + di_136[k];

            t_166[k] = -ab_y[k] * dh_103[k]
                       + di_137[k];

            t_167[k] = -ab_y[k] * dh_104[k]
                       + di_138[k];

            t_168[k] = -ab_y[k] * dh_105[k]
                       + di_141[k];
        }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, ab_y, dh_106, dh_107, dh_108, \
                         dh_109, dh_110, di_143, di_144, di_146, di_147, \
                         di_148 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_169[k] = -ab_y[k] * dh_106[k]
                       + di_143[k];

            t_170[k] = -ab_y[k] * dh_107[k]
                       + di_144[k];

            t_171[k] = -ab_y[k] * dh_108[k]
                       + di_146[k];

            t_172[k] = -ab_y[k] * dh_109[k]
                       + di_147[k];

            t_173[k] = -ab_y[k] * dh_110[k]
                       + di_148[k];
        }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, ab_y, dh_111, dh_112, dh_113, \
                         dh_114, dh_115, di_150, di_151, di_152, di_153, \
                         di_155 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_174[k] = -ab_y[k] * dh_111[k]
                       + di_150[k];

            t_175[k] = -ab_y[k] * dh_112[k]
                       + di_151[k];

            t_176[k] = -ab_y[k] * dh_113[k]
                       + di_152[k];

            t_177[k] = -ab_y[k] * dh_114[k]
                       + di_153[k];

            t_178[k] = -ab_y[k] * dh_115[k]
                       + di_155[k];
        }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, ab_y, dh_116, dh_117, dh_118, \
                         dh_119, dh_120, di_156, di_157, di_158, di_159, \
                         di_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_179[k] = -ab_y[k] * dh_116[k]
                       + di_156[k];

            t_180[k] = -ab_y[k] * dh_117[k]
                       + di_157[k];

            t_181[k] = -ab_y[k] * dh_118[k]
                       + di_158[k];

            t_182[k] = -ab_y[k] * dh_119[k]
                       + di_159[k];

            t_183[k] = -ab_y[k] * dh_120[k]
                       + di_161[k];
        }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, ab_y, dh_121, dh_122, dh_123, \
                         dh_124, dh_125, di_162, di_163, di_164, di_165, \
                         di_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_184[k] = -ab_y[k] * dh_121[k]
                       + di_162[k];

            t_185[k] = -ab_y[k] * dh_122[k]
                       + di_163[k];

            t_186[k] = -ab_y[k] * dh_123[k]
                       + di_164[k];

            t_187[k] = -ab_y[k] * dh_124[k]
                       + di_165[k];

            t_188[k] = -ab_y[k] * dh_125[k]
                       + di_166[k];
        }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, ab_z, dh_105, dh_106, dh_107, \
                         dh_108, dh_109, di_142, di_144, di_145, di_147, \
                         di_148 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_189[k] = -ab_z[k] * dh_105[k]
                       + di_142[k];

            t_190[k] = -ab_z[k] * dh_106[k]
                       + di_144[k];

            t_191[k] = -ab_z[k] * dh_107[k]
                       + di_145[k];

            t_192[k] = -ab_z[k] * dh_108[k]
                       + di_147[k];

            t_193[k] = -ab_z[k] * dh_109[k]
                       + di_148[k];
        }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, ab_z, dh_110, dh_111, dh_112, \
                         dh_113, dh_114, di_149, di_151, di_152, di_153, \
                         di_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_194[k] = -ab_z[k] * dh_110[k]
                       + di_149[k];

            t_195[k] = -ab_z[k] * dh_111[k]
                       + di_151[k];

            t_196[k] = -ab_z[k] * dh_112[k]
                       + di_152[k];

            t_197[k] = -ab_z[k] * dh_113[k]
                       + di_153[k];

            t_198[k] = -ab_z[k] * dh_114[k]
                       + di_154[k];
        }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, ab_z, dh_115, dh_116, dh_117, \
                         dh_118, dh_119, di_156, di_157, di_158, di_159, \
                         di_160 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_199[k] = -ab_z[k] * dh_115[k]
                       + di_156[k];

            t_200[k] = -ab_z[k] * dh_116[k]
                       + di_157[k];

            t_201[k] = -ab_z[k] * dh_117[k]
                       + di_158[k];

            t_202[k] = -ab_z[k] * dh_118[k]
                       + di_159[k];

            t_203[k] = -ab_z[k] * dh_119[k]
                       + di_160[k];
        }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, ab_z, dh_120, dh_121, dh_122, \
                         dh_123, dh_124, di_162, di_163, di_164, di_165, \
                         di_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_204[k] = -ab_z[k] * dh_120[k]
                       + di_162[k];

            t_205[k] = -ab_z[k] * dh_121[k]
                       + di_163[k];

            t_206[k] = -ab_z[k] * dh_122[k]
                       + di_164[k];

            t_207[k] = -ab_z[k] * dh_123[k]
                       + di_165[k];

            t_208[k] = -ab_z[k] * dh_124[k]
                       + di_166[k];
        }

#pragma omp simd aligned(t_209, ab_z, dh_125, di_167 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_209[k] = -ab_z[k] * dh_125[k]
                       + di_167[k];
        }
    }
}

auto
compute_hrr_fh(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t dh, const size_t di, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_fh_piece0(buffer, coordinates, target, dh, di, ncomps, nmax);

    compute_hrr_fh_piece1(buffer, coordinates, target, dh, di, ncomps, nmax);
}

}  // namespace simdtrf
