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


#include "SimdTransferIG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_ig_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t if_, const size_t kf,
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

        const auto *if__0 = buffer.data(if_ + 0 * ncomps + c);
        const auto *if__1 = buffer.data(if_ + 1 * ncomps + c);
        const auto *if__2 = buffer.data(if_ + 2 * ncomps + c);
        const auto *if__3 = buffer.data(if_ + 3 * ncomps + c);
        const auto *if__4 = buffer.data(if_ + 4 * ncomps + c);
        const auto *if__5 = buffer.data(if_ + 5 * ncomps + c);
        const auto *if__6 = buffer.data(if_ + 6 * ncomps + c);
        const auto *if__7 = buffer.data(if_ + 7 * ncomps + c);
        const auto *if__8 = buffer.data(if_ + 8 * ncomps + c);
        const auto *if__9 = buffer.data(if_ + 9 * ncomps + c);
        const auto *if__10 = buffer.data(if_ + 10 * ncomps + c);
        const auto *if__11 = buffer.data(if_ + 11 * ncomps + c);
        const auto *if__12 = buffer.data(if_ + 12 * ncomps + c);
        const auto *if__13 = buffer.data(if_ + 13 * ncomps + c);
        const auto *if__14 = buffer.data(if_ + 14 * ncomps + c);
        const auto *if__15 = buffer.data(if_ + 15 * ncomps + c);
        const auto *if__16 = buffer.data(if_ + 16 * ncomps + c);
        const auto *if__17 = buffer.data(if_ + 17 * ncomps + c);
        const auto *if__18 = buffer.data(if_ + 18 * ncomps + c);
        const auto *if__19 = buffer.data(if_ + 19 * ncomps + c);
        const auto *if__20 = buffer.data(if_ + 20 * ncomps + c);
        const auto *if__21 = buffer.data(if_ + 21 * ncomps + c);
        const auto *if__22 = buffer.data(if_ + 22 * ncomps + c);
        const auto *if__23 = buffer.data(if_ + 23 * ncomps + c);
        const auto *if__24 = buffer.data(if_ + 24 * ncomps + c);
        const auto *if__25 = buffer.data(if_ + 25 * ncomps + c);
        const auto *if__26 = buffer.data(if_ + 26 * ncomps + c);
        const auto *if__27 = buffer.data(if_ + 27 * ncomps + c);
        const auto *if__28 = buffer.data(if_ + 28 * ncomps + c);
        const auto *if__29 = buffer.data(if_ + 29 * ncomps + c);
        const auto *if__30 = buffer.data(if_ + 30 * ncomps + c);
        const auto *if__31 = buffer.data(if_ + 31 * ncomps + c);
        const auto *if__32 = buffer.data(if_ + 32 * ncomps + c);
        const auto *if__33 = buffer.data(if_ + 33 * ncomps + c);
        const auto *if__34 = buffer.data(if_ + 34 * ncomps + c);
        const auto *if__35 = buffer.data(if_ + 35 * ncomps + c);
        const auto *if__36 = buffer.data(if_ + 36 * ncomps + c);
        const auto *if__37 = buffer.data(if_ + 37 * ncomps + c);
        const auto *if__38 = buffer.data(if_ + 38 * ncomps + c);
        const auto *if__39 = buffer.data(if_ + 39 * ncomps + c);
        const auto *if__40 = buffer.data(if_ + 40 * ncomps + c);
        const auto *if__41 = buffer.data(if_ + 41 * ncomps + c);
        const auto *if__42 = buffer.data(if_ + 42 * ncomps + c);
        const auto *if__43 = buffer.data(if_ + 43 * ncomps + c);
        const auto *if__44 = buffer.data(if_ + 44 * ncomps + c);
        const auto *if__45 = buffer.data(if_ + 45 * ncomps + c);
        const auto *if__46 = buffer.data(if_ + 46 * ncomps + c);
        const auto *if__47 = buffer.data(if_ + 47 * ncomps + c);
        const auto *if__48 = buffer.data(if_ + 48 * ncomps + c);
        const auto *if__49 = buffer.data(if_ + 49 * ncomps + c);
        const auto *if__50 = buffer.data(if_ + 50 * ncomps + c);
        const auto *if__51 = buffer.data(if_ + 51 * ncomps + c);
        const auto *if__52 = buffer.data(if_ + 52 * ncomps + c);
        const auto *if__53 = buffer.data(if_ + 53 * ncomps + c);
        const auto *if__54 = buffer.data(if_ + 54 * ncomps + c);
        const auto *if__55 = buffer.data(if_ + 55 * ncomps + c);
        const auto *if__56 = buffer.data(if_ + 56 * ncomps + c);
        const auto *if__57 = buffer.data(if_ + 57 * ncomps + c);
        const auto *if__58 = buffer.data(if_ + 58 * ncomps + c);
        const auto *if__59 = buffer.data(if_ + 59 * ncomps + c);
        const auto *if__60 = buffer.data(if_ + 60 * ncomps + c);
        const auto *if__61 = buffer.data(if_ + 61 * ncomps + c);
        const auto *if__62 = buffer.data(if_ + 62 * ncomps + c);
        const auto *if__63 = buffer.data(if_ + 63 * ncomps + c);
        const auto *if__64 = buffer.data(if_ + 64 * ncomps + c);
        const auto *if__65 = buffer.data(if_ + 65 * ncomps + c);
        const auto *if__66 = buffer.data(if_ + 66 * ncomps + c);
        const auto *if__67 = buffer.data(if_ + 67 * ncomps + c);
        const auto *if__68 = buffer.data(if_ + 68 * ncomps + c);
        const auto *if__69 = buffer.data(if_ + 69 * ncomps + c);
        const auto *if__70 = buffer.data(if_ + 70 * ncomps + c);
        const auto *if__71 = buffer.data(if_ + 71 * ncomps + c);
        const auto *if__72 = buffer.data(if_ + 72 * ncomps + c);
        const auto *if__73 = buffer.data(if_ + 73 * ncomps + c);
        const auto *if__74 = buffer.data(if_ + 74 * ncomps + c);
        const auto *if__75 = buffer.data(if_ + 75 * ncomps + c);
        const auto *if__76 = buffer.data(if_ + 76 * ncomps + c);
        const auto *if__77 = buffer.data(if_ + 77 * ncomps + c);
        const auto *if__78 = buffer.data(if_ + 78 * ncomps + c);
        const auto *if__79 = buffer.data(if_ + 79 * ncomps + c);
        const auto *if__80 = buffer.data(if_ + 80 * ncomps + c);
        const auto *if__81 = buffer.data(if_ + 81 * ncomps + c);
        const auto *if__82 = buffer.data(if_ + 82 * ncomps + c);
        const auto *if__83 = buffer.data(if_ + 83 * ncomps + c);
        const auto *if__84 = buffer.data(if_ + 84 * ncomps + c);
        const auto *if__85 = buffer.data(if_ + 85 * ncomps + c);
        const auto *if__86 = buffer.data(if_ + 86 * ncomps + c);
        const auto *if__87 = buffer.data(if_ + 87 * ncomps + c);
        const auto *if__88 = buffer.data(if_ + 88 * ncomps + c);
        const auto *if__89 = buffer.data(if_ + 89 * ncomps + c);
        const auto *if__90 = buffer.data(if_ + 90 * ncomps + c);
        const auto *if__91 = buffer.data(if_ + 91 * ncomps + c);
        const auto *if__92 = buffer.data(if_ + 92 * ncomps + c);
        const auto *if__93 = buffer.data(if_ + 93 * ncomps + c);
        const auto *if__94 = buffer.data(if_ + 94 * ncomps + c);
        const auto *if__95 = buffer.data(if_ + 95 * ncomps + c);
        const auto *if__96 = buffer.data(if_ + 96 * ncomps + c);
        const auto *if__97 = buffer.data(if_ + 97 * ncomps + c);
        const auto *if__98 = buffer.data(if_ + 98 * ncomps + c);
        const auto *if__99 = buffer.data(if_ + 99 * ncomps + c);

        const auto *kf_0 = buffer.data(kf + 0 * ncomps + c);
        const auto *kf_1 = buffer.data(kf + 1 * ncomps + c);
        const auto *kf_2 = buffer.data(kf + 2 * ncomps + c);
        const auto *kf_3 = buffer.data(kf + 3 * ncomps + c);
        const auto *kf_4 = buffer.data(kf + 4 * ncomps + c);
        const auto *kf_5 = buffer.data(kf + 5 * ncomps + c);
        const auto *kf_6 = buffer.data(kf + 6 * ncomps + c);
        const auto *kf_7 = buffer.data(kf + 7 * ncomps + c);
        const auto *kf_8 = buffer.data(kf + 8 * ncomps + c);
        const auto *kf_9 = buffer.data(kf + 9 * ncomps + c);
        const auto *kf_10 = buffer.data(kf + 10 * ncomps + c);
        const auto *kf_11 = buffer.data(kf + 11 * ncomps + c);
        const auto *kf_12 = buffer.data(kf + 12 * ncomps + c);
        const auto *kf_13 = buffer.data(kf + 13 * ncomps + c);
        const auto *kf_14 = buffer.data(kf + 14 * ncomps + c);
        const auto *kf_15 = buffer.data(kf + 15 * ncomps + c);
        const auto *kf_16 = buffer.data(kf + 16 * ncomps + c);
        const auto *kf_17 = buffer.data(kf + 17 * ncomps + c);
        const auto *kf_18 = buffer.data(kf + 18 * ncomps + c);
        const auto *kf_19 = buffer.data(kf + 19 * ncomps + c);
        const auto *kf_20 = buffer.data(kf + 20 * ncomps + c);
        const auto *kf_21 = buffer.data(kf + 21 * ncomps + c);
        const auto *kf_22 = buffer.data(kf + 22 * ncomps + c);
        const auto *kf_23 = buffer.data(kf + 23 * ncomps + c);
        const auto *kf_24 = buffer.data(kf + 24 * ncomps + c);
        const auto *kf_25 = buffer.data(kf + 25 * ncomps + c);
        const auto *kf_26 = buffer.data(kf + 26 * ncomps + c);
        const auto *kf_27 = buffer.data(kf + 27 * ncomps + c);
        const auto *kf_28 = buffer.data(kf + 28 * ncomps + c);
        const auto *kf_29 = buffer.data(kf + 29 * ncomps + c);
        const auto *kf_30 = buffer.data(kf + 30 * ncomps + c);
        const auto *kf_31 = buffer.data(kf + 31 * ncomps + c);
        const auto *kf_32 = buffer.data(kf + 32 * ncomps + c);
        const auto *kf_33 = buffer.data(kf + 33 * ncomps + c);
        const auto *kf_34 = buffer.data(kf + 34 * ncomps + c);
        const auto *kf_35 = buffer.data(kf + 35 * ncomps + c);
        const auto *kf_36 = buffer.data(kf + 36 * ncomps + c);
        const auto *kf_37 = buffer.data(kf + 37 * ncomps + c);
        const auto *kf_38 = buffer.data(kf + 38 * ncomps + c);
        const auto *kf_39 = buffer.data(kf + 39 * ncomps + c);
        const auto *kf_40 = buffer.data(kf + 40 * ncomps + c);
        const auto *kf_41 = buffer.data(kf + 41 * ncomps + c);
        const auto *kf_42 = buffer.data(kf + 42 * ncomps + c);
        const auto *kf_43 = buffer.data(kf + 43 * ncomps + c);
        const auto *kf_44 = buffer.data(kf + 44 * ncomps + c);
        const auto *kf_45 = buffer.data(kf + 45 * ncomps + c);
        const auto *kf_46 = buffer.data(kf + 46 * ncomps + c);
        const auto *kf_47 = buffer.data(kf + 47 * ncomps + c);
        const auto *kf_48 = buffer.data(kf + 48 * ncomps + c);
        const auto *kf_49 = buffer.data(kf + 49 * ncomps + c);
        const auto *kf_50 = buffer.data(kf + 50 * ncomps + c);
        const auto *kf_51 = buffer.data(kf + 51 * ncomps + c);
        const auto *kf_52 = buffer.data(kf + 52 * ncomps + c);
        const auto *kf_53 = buffer.data(kf + 53 * ncomps + c);
        const auto *kf_54 = buffer.data(kf + 54 * ncomps + c);
        const auto *kf_55 = buffer.data(kf + 55 * ncomps + c);
        const auto *kf_56 = buffer.data(kf + 56 * ncomps + c);
        const auto *kf_57 = buffer.data(kf + 57 * ncomps + c);
        const auto *kf_58 = buffer.data(kf + 58 * ncomps + c);
        const auto *kf_59 = buffer.data(kf + 59 * ncomps + c);
        const auto *kf_60 = buffer.data(kf + 60 * ncomps + c);
        const auto *kf_61 = buffer.data(kf + 61 * ncomps + c);
        const auto *kf_62 = buffer.data(kf + 62 * ncomps + c);
        const auto *kf_63 = buffer.data(kf + 63 * ncomps + c);
        const auto *kf_64 = buffer.data(kf + 64 * ncomps + c);
        const auto *kf_65 = buffer.data(kf + 65 * ncomps + c);
        const auto *kf_66 = buffer.data(kf + 66 * ncomps + c);
        const auto *kf_67 = buffer.data(kf + 67 * ncomps + c);
        const auto *kf_68 = buffer.data(kf + 68 * ncomps + c);
        const auto *kf_69 = buffer.data(kf + 69 * ncomps + c);
        const auto *kf_70 = buffer.data(kf + 70 * ncomps + c);
        const auto *kf_71 = buffer.data(kf + 71 * ncomps + c);
        const auto *kf_72 = buffer.data(kf + 72 * ncomps + c);
        const auto *kf_73 = buffer.data(kf + 73 * ncomps + c);
        const auto *kf_74 = buffer.data(kf + 74 * ncomps + c);
        const auto *kf_75 = buffer.data(kf + 75 * ncomps + c);
        const auto *kf_76 = buffer.data(kf + 76 * ncomps + c);
        const auto *kf_77 = buffer.data(kf + 77 * ncomps + c);
        const auto *kf_78 = buffer.data(kf + 78 * ncomps + c);
        const auto *kf_79 = buffer.data(kf + 79 * ncomps + c);
        const auto *kf_80 = buffer.data(kf + 80 * ncomps + c);
        const auto *kf_81 = buffer.data(kf + 81 * ncomps + c);
        const auto *kf_82 = buffer.data(kf + 82 * ncomps + c);
        const auto *kf_83 = buffer.data(kf + 83 * ncomps + c);
        const auto *kf_84 = buffer.data(kf + 84 * ncomps + c);
        const auto *kf_85 = buffer.data(kf + 85 * ncomps + c);
        const auto *kf_86 = buffer.data(kf + 86 * ncomps + c);
        const auto *kf_87 = buffer.data(kf + 87 * ncomps + c);
        const auto *kf_88 = buffer.data(kf + 88 * ncomps + c);
        const auto *kf_89 = buffer.data(kf + 89 * ncomps + c);
        const auto *kf_90 = buffer.data(kf + 90 * ncomps + c);
        const auto *kf_91 = buffer.data(kf + 91 * ncomps + c);
        const auto *kf_92 = buffer.data(kf + 92 * ncomps + c);
        const auto *kf_93 = buffer.data(kf + 93 * ncomps + c);
        const auto *kf_94 = buffer.data(kf + 94 * ncomps + c);
        const auto *kf_95 = buffer.data(kf + 95 * ncomps + c);
        const auto *kf_96 = buffer.data(kf + 96 * ncomps + c);
        const auto *kf_97 = buffer.data(kf + 97 * ncomps + c);
        const auto *kf_98 = buffer.data(kf + 98 * ncomps + c);
        const auto *kf_99 = buffer.data(kf + 99 * ncomps + c);
        const auto *kf_106 = buffer.data(kf + 106 * ncomps + c);
        const auto *kf_107 = buffer.data(kf + 107 * ncomps + c);
        const auto *kf_108 = buffer.data(kf + 108 * ncomps + c);
        const auto *kf_109 = buffer.data(kf + 109 * ncomps + c);
        const auto *kf_116 = buffer.data(kf + 116 * ncomps + c);
        const auto *kf_117 = buffer.data(kf + 117 * ncomps + c);
        const auto *kf_118 = buffer.data(kf + 118 * ncomps + c);
        const auto *kf_119 = buffer.data(kf + 119 * ncomps + c);
        const auto *kf_126 = buffer.data(kf + 126 * ncomps + c);
        const auto *kf_127 = buffer.data(kf + 127 * ncomps + c);
        const auto *kf_128 = buffer.data(kf + 128 * ncomps + c);
        const auto *kf_129 = buffer.data(kf + 129 * ncomps + c);
        const auto *kf_139 = buffer.data(kf + 139 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, if__0, if__1, if__2, if__3, if__4, \
                         kf_0, kf_1, kf_2, kf_3, kf_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * if__0[k]
                     + kf_0[k];

            t_1[k] = ab_x[k] * if__1[k]
                     + kf_1[k];

            t_2[k] = ab_x[k] * if__2[k]
                     + kf_2[k];

            t_3[k] = ab_x[k] * if__3[k]
                     + kf_3[k];

            t_4[k] = ab_x[k] * if__4[k]
                     + kf_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, if__5, if__6, if__7, if__8, if__9, \
                         kf_5, kf_6, kf_7, kf_8, kf_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * if__5[k]
                     + kf_5[k];

            t_6[k] = ab_x[k] * if__6[k]
                     + kf_6[k];

            t_7[k] = ab_x[k] * if__7[k]
                     + kf_7[k];

            t_8[k] = ab_x[k] * if__8[k]
                     + kf_8[k];

            t_9[k] = ab_x[k] * if__9[k]
                     + kf_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_y, ab_z, if__6, if__7, if__8, if__9, \
                         kf_16, kf_17, kf_18, kf_19, kf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_y[k] * if__6[k]
                      + kf_16[k];

            t_11[k] = ab_y[k] * if__7[k]
                      + kf_17[k];

            t_12[k] = ab_y[k] * if__8[k]
                      + kf_18[k];

            t_13[k] = ab_y[k] * if__9[k]
                      + kf_19[k];

            t_14[k] = ab_z[k] * if__9[k]
                      + kf_29[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, if__10, if__11, if__12, if__13, \
                         if__14, kf_10, kf_11, kf_12, kf_13, kf_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * if__10[k]
                      + kf_10[k];

            t_16[k] = ab_x[k] * if__11[k]
                      + kf_11[k];

            t_17[k] = ab_x[k] * if__12[k]
                      + kf_12[k];

            t_18[k] = ab_x[k] * if__13[k]
                      + kf_13[k];

            t_19[k] = ab_x[k] * if__14[k]
                      + kf_14[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, if__15, if__16, if__17, if__18, \
                         if__19, kf_15, kf_16, kf_17, kf_18, kf_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * if__15[k]
                      + kf_15[k];

            t_21[k] = ab_x[k] * if__16[k]
                      + kf_16[k];

            t_22[k] = ab_x[k] * if__17[k]
                      + kf_17[k];

            t_23[k] = ab_x[k] * if__18[k]
                      + kf_18[k];

            t_24[k] = ab_x[k] * if__19[k]
                      + kf_19[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_y, ab_z, if__16, if__17, if__18, \
                         if__19, kf_36, kf_37, kf_38, kf_39, kf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_y[k] * if__16[k]
                      + kf_36[k];

            t_26[k] = ab_y[k] * if__17[k]
                      + kf_37[k];

            t_27[k] = ab_y[k] * if__18[k]
                      + kf_38[k];

            t_28[k] = ab_y[k] * if__19[k]
                      + kf_39[k];

            t_29[k] = ab_z[k] * if__19[k]
                      + kf_49[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, if__20, if__21, if__22, if__23, \
                         if__24, kf_20, kf_21, kf_22, kf_23, kf_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * if__20[k]
                      + kf_20[k];

            t_31[k] = ab_x[k] * if__21[k]
                      + kf_21[k];

            t_32[k] = ab_x[k] * if__22[k]
                      + kf_22[k];

            t_33[k] = ab_x[k] * if__23[k]
                      + kf_23[k];

            t_34[k] = ab_x[k] * if__24[k]
                      + kf_24[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, if__25, if__26, if__27, if__28, \
                         if__29, kf_25, kf_26, kf_27, kf_28, kf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * if__25[k]
                      + kf_25[k];

            t_36[k] = ab_x[k] * if__26[k]
                      + kf_26[k];

            t_37[k] = ab_x[k] * if__27[k]
                      + kf_27[k];

            t_38[k] = ab_x[k] * if__28[k]
                      + kf_28[k];

            t_39[k] = ab_x[k] * if__29[k]
                      + kf_29[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_y, ab_z, if__26, if__27, if__28, \
                         if__29, kf_46, kf_47, kf_48, kf_49, kf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_y[k] * if__26[k]
                      + kf_46[k];

            t_41[k] = ab_y[k] * if__27[k]
                      + kf_47[k];

            t_42[k] = ab_y[k] * if__28[k]
                      + kf_48[k];

            t_43[k] = ab_y[k] * if__29[k]
                      + kf_49[k];

            t_44[k] = ab_z[k] * if__29[k]
                      + kf_59[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, if__30, if__31, if__32, if__33, \
                         if__34, kf_30, kf_31, kf_32, kf_33, kf_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * if__30[k]
                      + kf_30[k];

            t_46[k] = ab_x[k] * if__31[k]
                      + kf_31[k];

            t_47[k] = ab_x[k] * if__32[k]
                      + kf_32[k];

            t_48[k] = ab_x[k] * if__33[k]
                      + kf_33[k];

            t_49[k] = ab_x[k] * if__34[k]
                      + kf_34[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, if__35, if__36, if__37, if__38, \
                         if__39, kf_35, kf_36, kf_37, kf_38, kf_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * if__35[k]
                      + kf_35[k];

            t_51[k] = ab_x[k] * if__36[k]
                      + kf_36[k];

            t_52[k] = ab_x[k] * if__37[k]
                      + kf_37[k];

            t_53[k] = ab_x[k] * if__38[k]
                      + kf_38[k];

            t_54[k] = ab_x[k] * if__39[k]
                      + kf_39[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_y, ab_z, if__36, if__37, if__38, \
                         if__39, kf_66, kf_67, kf_68, kf_69, kf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_y[k] * if__36[k]
                      + kf_66[k];

            t_56[k] = ab_y[k] * if__37[k]
                      + kf_67[k];

            t_57[k] = ab_y[k] * if__38[k]
                      + kf_68[k];

            t_58[k] = ab_y[k] * if__39[k]
                      + kf_69[k];

            t_59[k] = ab_z[k] * if__39[k]
                      + kf_79[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, if__40, if__41, if__42, if__43, \
                         if__44, kf_40, kf_41, kf_42, kf_43, kf_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * if__40[k]
                      + kf_40[k];

            t_61[k] = ab_x[k] * if__41[k]
                      + kf_41[k];

            t_62[k] = ab_x[k] * if__42[k]
                      + kf_42[k];

            t_63[k] = ab_x[k] * if__43[k]
                      + kf_43[k];

            t_64[k] = ab_x[k] * if__44[k]
                      + kf_44[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, if__45, if__46, if__47, if__48, \
                         if__49, kf_45, kf_46, kf_47, kf_48, kf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * if__45[k]
                      + kf_45[k];

            t_66[k] = ab_x[k] * if__46[k]
                      + kf_46[k];

            t_67[k] = ab_x[k] * if__47[k]
                      + kf_47[k];

            t_68[k] = ab_x[k] * if__48[k]
                      + kf_48[k];

            t_69[k] = ab_x[k] * if__49[k]
                      + kf_49[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, ab_z, if__46, if__47, if__48, \
                         if__49, kf_76, kf_77, kf_78, kf_79, kf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_y[k] * if__46[k]
                      + kf_76[k];

            t_71[k] = ab_y[k] * if__47[k]
                      + kf_77[k];

            t_72[k] = ab_y[k] * if__48[k]
                      + kf_78[k];

            t_73[k] = ab_y[k] * if__49[k]
                      + kf_79[k];

            t_74[k] = ab_z[k] * if__49[k]
                      + kf_89[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, if__50, if__51, if__52, if__53, \
                         if__54, kf_50, kf_51, kf_52, kf_53, kf_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * if__50[k]
                      + kf_50[k];

            t_76[k] = ab_x[k] * if__51[k]
                      + kf_51[k];

            t_77[k] = ab_x[k] * if__52[k]
                      + kf_52[k];

            t_78[k] = ab_x[k] * if__53[k]
                      + kf_53[k];

            t_79[k] = ab_x[k] * if__54[k]
                      + kf_54[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, if__55, if__56, if__57, if__58, \
                         if__59, kf_55, kf_56, kf_57, kf_58, kf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * if__55[k]
                      + kf_55[k];

            t_81[k] = ab_x[k] * if__56[k]
                      + kf_56[k];

            t_82[k] = ab_x[k] * if__57[k]
                      + kf_57[k];

            t_83[k] = ab_x[k] * if__58[k]
                      + kf_58[k];

            t_84[k] = ab_x[k] * if__59[k]
                      + kf_59[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_y, ab_z, if__56, if__57, if__58, \
                         if__59, kf_86, kf_87, kf_88, kf_89, kf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_y[k] * if__56[k]
                      + kf_86[k];

            t_86[k] = ab_y[k] * if__57[k]
                      + kf_87[k];

            t_87[k] = ab_y[k] * if__58[k]
                      + kf_88[k];

            t_88[k] = ab_y[k] * if__59[k]
                      + kf_89[k];

            t_89[k] = ab_z[k] * if__59[k]
                      + kf_99[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, if__60, if__61, if__62, if__63, \
                         if__64, kf_60, kf_61, kf_62, kf_63, kf_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * if__60[k]
                      + kf_60[k];

            t_91[k] = ab_x[k] * if__61[k]
                      + kf_61[k];

            t_92[k] = ab_x[k] * if__62[k]
                      + kf_62[k];

            t_93[k] = ab_x[k] * if__63[k]
                      + kf_63[k];

            t_94[k] = ab_x[k] * if__64[k]
                      + kf_64[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, if__65, if__66, if__67, if__68, \
                         if__69, kf_65, kf_66, kf_67, kf_68, kf_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * if__65[k]
                      + kf_65[k];

            t_96[k] = ab_x[k] * if__66[k]
                      + kf_66[k];

            t_97[k] = ab_x[k] * if__67[k]
                      + kf_67[k];

            t_98[k] = ab_x[k] * if__68[k]
                      + kf_68[k];

            t_99[k] = ab_x[k] * if__69[k]
                      + kf_69[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_y, ab_z, if__66, if__67, \
                         if__68, if__69, kf_106, kf_107, kf_108, kf_109, \
                         kf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_y[k] * if__66[k]
                       + kf_106[k];

            t_101[k] = ab_y[k] * if__67[k]
                       + kf_107[k];

            t_102[k] = ab_y[k] * if__68[k]
                       + kf_108[k];

            t_103[k] = ab_y[k] * if__69[k]
                       + kf_109[k];

            t_104[k] = ab_z[k] * if__69[k]
                       + kf_119[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, if__70, if__71, if__72, \
                         if__73, if__74, kf_70, kf_71, kf_72, kf_73, \
                         kf_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * if__70[k]
                       + kf_70[k];

            t_106[k] = ab_x[k] * if__71[k]
                       + kf_71[k];

            t_107[k] = ab_x[k] * if__72[k]
                       + kf_72[k];

            t_108[k] = ab_x[k] * if__73[k]
                       + kf_73[k];

            t_109[k] = ab_x[k] * if__74[k]
                       + kf_74[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, if__75, if__76, if__77, \
                         if__78, if__79, kf_75, kf_76, kf_77, kf_78, \
                         kf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_x[k] * if__75[k]
                       + kf_75[k];

            t_111[k] = ab_x[k] * if__76[k]
                       + kf_76[k];

            t_112[k] = ab_x[k] * if__77[k]
                       + kf_77[k];

            t_113[k] = ab_x[k] * if__78[k]
                       + kf_78[k];

            t_114[k] = ab_x[k] * if__79[k]
                       + kf_79[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_y, ab_z, if__76, if__77, \
                         if__78, if__79, kf_116, kf_117, kf_118, kf_119, \
                         kf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_y[k] * if__76[k]
                       + kf_116[k];

            t_116[k] = ab_y[k] * if__77[k]
                       + kf_117[k];

            t_117[k] = ab_y[k] * if__78[k]
                       + kf_118[k];

            t_118[k] = ab_y[k] * if__79[k]
                       + kf_119[k];

            t_119[k] = ab_z[k] * if__79[k]
                       + kf_129[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, if__80, if__81, if__82, \
                         if__83, if__84, kf_80, kf_81, kf_82, kf_83, \
                         kf_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * if__80[k]
                       + kf_80[k];

            t_121[k] = ab_x[k] * if__81[k]
                       + kf_81[k];

            t_122[k] = ab_x[k] * if__82[k]
                       + kf_82[k];

            t_123[k] = ab_x[k] * if__83[k]
                       + kf_83[k];

            t_124[k] = ab_x[k] * if__84[k]
                       + kf_84[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, if__85, if__86, if__87, \
                         if__88, if__89, kf_85, kf_86, kf_87, kf_88, \
                         kf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * if__85[k]
                       + kf_85[k];

            t_126[k] = ab_x[k] * if__86[k]
                       + kf_86[k];

            t_127[k] = ab_x[k] * if__87[k]
                       + kf_87[k];

            t_128[k] = ab_x[k] * if__88[k]
                       + kf_88[k];

            t_129[k] = ab_x[k] * if__89[k]
                       + kf_89[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_y, ab_z, if__86, if__87, \
                         if__88, if__89, kf_126, kf_127, kf_128, kf_129, \
                         kf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_y[k] * if__86[k]
                       + kf_126[k];

            t_131[k] = ab_y[k] * if__87[k]
                       + kf_127[k];

            t_132[k] = ab_y[k] * if__88[k]
                       + kf_128[k];

            t_133[k] = ab_y[k] * if__89[k]
                       + kf_129[k];

            t_134[k] = ab_z[k] * if__89[k]
                       + kf_139[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, if__90, if__91, if__92, \
                         if__93, if__94, kf_90, kf_91, kf_92, kf_93, \
                         kf_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * if__90[k]
                       + kf_90[k];

            t_136[k] = ab_x[k] * if__91[k]
                       + kf_91[k];

            t_137[k] = ab_x[k] * if__92[k]
                       + kf_92[k];

            t_138[k] = ab_x[k] * if__93[k]
                       + kf_93[k];

            t_139[k] = ab_x[k] * if__94[k]
                       + kf_94[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, if__95, if__96, if__97, \
                         if__98, if__99, kf_95, kf_96, kf_97, kf_98, \
                         kf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * if__95[k]
                       + kf_95[k];

            t_141[k] = ab_x[k] * if__96[k]
                       + kf_96[k];

            t_142[k] = ab_x[k] * if__97[k]
                       + kf_97[k];

            t_143[k] = ab_x[k] * if__98[k]
                       + kf_98[k];

            t_144[k] = ab_x[k] * if__99[k]
                       + kf_99[k];
        }
    }
}

static auto
compute_hrr_ig_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t if_, const size_t kf,
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

        const auto *if__96 = buffer.data(if_ + 96 * ncomps + c);
        const auto *if__97 = buffer.data(if_ + 97 * ncomps + c);
        const auto *if__98 = buffer.data(if_ + 98 * ncomps + c);
        const auto *if__99 = buffer.data(if_ + 99 * ncomps + c);
        const auto *if__100 = buffer.data(if_ + 100 * ncomps + c);
        const auto *if__101 = buffer.data(if_ + 101 * ncomps + c);
        const auto *if__102 = buffer.data(if_ + 102 * ncomps + c);
        const auto *if__103 = buffer.data(if_ + 103 * ncomps + c);
        const auto *if__104 = buffer.data(if_ + 104 * ncomps + c);
        const auto *if__105 = buffer.data(if_ + 105 * ncomps + c);
        const auto *if__106 = buffer.data(if_ + 106 * ncomps + c);
        const auto *if__107 = buffer.data(if_ + 107 * ncomps + c);
        const auto *if__108 = buffer.data(if_ + 108 * ncomps + c);
        const auto *if__109 = buffer.data(if_ + 109 * ncomps + c);
        const auto *if__110 = buffer.data(if_ + 110 * ncomps + c);
        const auto *if__111 = buffer.data(if_ + 111 * ncomps + c);
        const auto *if__112 = buffer.data(if_ + 112 * ncomps + c);
        const auto *if__113 = buffer.data(if_ + 113 * ncomps + c);
        const auto *if__114 = buffer.data(if_ + 114 * ncomps + c);
        const auto *if__115 = buffer.data(if_ + 115 * ncomps + c);
        const auto *if__116 = buffer.data(if_ + 116 * ncomps + c);
        const auto *if__117 = buffer.data(if_ + 117 * ncomps + c);
        const auto *if__118 = buffer.data(if_ + 118 * ncomps + c);
        const auto *if__119 = buffer.data(if_ + 119 * ncomps + c);
        const auto *if__120 = buffer.data(if_ + 120 * ncomps + c);
        const auto *if__121 = buffer.data(if_ + 121 * ncomps + c);
        const auto *if__122 = buffer.data(if_ + 122 * ncomps + c);
        const auto *if__123 = buffer.data(if_ + 123 * ncomps + c);
        const auto *if__124 = buffer.data(if_ + 124 * ncomps + c);
        const auto *if__125 = buffer.data(if_ + 125 * ncomps + c);
        const auto *if__126 = buffer.data(if_ + 126 * ncomps + c);
        const auto *if__127 = buffer.data(if_ + 127 * ncomps + c);
        const auto *if__128 = buffer.data(if_ + 128 * ncomps + c);
        const auto *if__129 = buffer.data(if_ + 129 * ncomps + c);
        const auto *if__130 = buffer.data(if_ + 130 * ncomps + c);
        const auto *if__131 = buffer.data(if_ + 131 * ncomps + c);
        const auto *if__132 = buffer.data(if_ + 132 * ncomps + c);
        const auto *if__133 = buffer.data(if_ + 133 * ncomps + c);
        const auto *if__134 = buffer.data(if_ + 134 * ncomps + c);
        const auto *if__135 = buffer.data(if_ + 135 * ncomps + c);
        const auto *if__136 = buffer.data(if_ + 136 * ncomps + c);
        const auto *if__137 = buffer.data(if_ + 137 * ncomps + c);
        const auto *if__138 = buffer.data(if_ + 138 * ncomps + c);
        const auto *if__139 = buffer.data(if_ + 139 * ncomps + c);
        const auto *if__140 = buffer.data(if_ + 140 * ncomps + c);
        const auto *if__141 = buffer.data(if_ + 141 * ncomps + c);
        const auto *if__142 = buffer.data(if_ + 142 * ncomps + c);
        const auto *if__143 = buffer.data(if_ + 143 * ncomps + c);
        const auto *if__144 = buffer.data(if_ + 144 * ncomps + c);
        const auto *if__145 = buffer.data(if_ + 145 * ncomps + c);
        const auto *if__146 = buffer.data(if_ + 146 * ncomps + c);
        const auto *if__147 = buffer.data(if_ + 147 * ncomps + c);
        const auto *if__148 = buffer.data(if_ + 148 * ncomps + c);
        const auto *if__149 = buffer.data(if_ + 149 * ncomps + c);
        const auto *if__150 = buffer.data(if_ + 150 * ncomps + c);
        const auto *if__151 = buffer.data(if_ + 151 * ncomps + c);
        const auto *if__152 = buffer.data(if_ + 152 * ncomps + c);
        const auto *if__153 = buffer.data(if_ + 153 * ncomps + c);
        const auto *if__154 = buffer.data(if_ + 154 * ncomps + c);
        const auto *if__155 = buffer.data(if_ + 155 * ncomps + c);
        const auto *if__156 = buffer.data(if_ + 156 * ncomps + c);
        const auto *if__157 = buffer.data(if_ + 157 * ncomps + c);
        const auto *if__158 = buffer.data(if_ + 158 * ncomps + c);
        const auto *if__159 = buffer.data(if_ + 159 * ncomps + c);
        const auto *if__160 = buffer.data(if_ + 160 * ncomps + c);
        const auto *if__161 = buffer.data(if_ + 161 * ncomps + c);
        const auto *if__162 = buffer.data(if_ + 162 * ncomps + c);
        const auto *if__163 = buffer.data(if_ + 163 * ncomps + c);
        const auto *if__164 = buffer.data(if_ + 164 * ncomps + c);
        const auto *if__165 = buffer.data(if_ + 165 * ncomps + c);
        const auto *if__166 = buffer.data(if_ + 166 * ncomps + c);
        const auto *if__167 = buffer.data(if_ + 167 * ncomps + c);
        const auto *if__168 = buffer.data(if_ + 168 * ncomps + c);
        const auto *if__169 = buffer.data(if_ + 169 * ncomps + c);
        const auto *if__170 = buffer.data(if_ + 170 * ncomps + c);
        const auto *if__171 = buffer.data(if_ + 171 * ncomps + c);
        const auto *if__172 = buffer.data(if_ + 172 * ncomps + c);
        const auto *if__173 = buffer.data(if_ + 173 * ncomps + c);
        const auto *if__174 = buffer.data(if_ + 174 * ncomps + c);
        const auto *if__175 = buffer.data(if_ + 175 * ncomps + c);
        const auto *if__176 = buffer.data(if_ + 176 * ncomps + c);
        const auto *if__177 = buffer.data(if_ + 177 * ncomps + c);
        const auto *if__178 = buffer.data(if_ + 178 * ncomps + c);
        const auto *if__179 = buffer.data(if_ + 179 * ncomps + c);
        const auto *if__180 = buffer.data(if_ + 180 * ncomps + c);
        const auto *if__181 = buffer.data(if_ + 181 * ncomps + c);
        const auto *if__182 = buffer.data(if_ + 182 * ncomps + c);
        const auto *if__183 = buffer.data(if_ + 183 * ncomps + c);
        const auto *if__184 = buffer.data(if_ + 184 * ncomps + c);
        const auto *if__185 = buffer.data(if_ + 185 * ncomps + c);
        const auto *if__186 = buffer.data(if_ + 186 * ncomps + c);
        const auto *if__187 = buffer.data(if_ + 187 * ncomps + c);
        const auto *if__188 = buffer.data(if_ + 188 * ncomps + c);
        const auto *if__189 = buffer.data(if_ + 189 * ncomps + c);
        const auto *if__190 = buffer.data(if_ + 190 * ncomps + c);
        const auto *if__191 = buffer.data(if_ + 191 * ncomps + c);
        const auto *if__192 = buffer.data(if_ + 192 * ncomps + c);
        const auto *if__193 = buffer.data(if_ + 193 * ncomps + c);
        const auto *if__194 = buffer.data(if_ + 194 * ncomps + c);

        const auto *kf_100 = buffer.data(kf + 100 * ncomps + c);
        const auto *kf_101 = buffer.data(kf + 101 * ncomps + c);
        const auto *kf_102 = buffer.data(kf + 102 * ncomps + c);
        const auto *kf_103 = buffer.data(kf + 103 * ncomps + c);
        const auto *kf_104 = buffer.data(kf + 104 * ncomps + c);
        const auto *kf_105 = buffer.data(kf + 105 * ncomps + c);
        const auto *kf_106 = buffer.data(kf + 106 * ncomps + c);
        const auto *kf_107 = buffer.data(kf + 107 * ncomps + c);
        const auto *kf_108 = buffer.data(kf + 108 * ncomps + c);
        const auto *kf_109 = buffer.data(kf + 109 * ncomps + c);
        const auto *kf_110 = buffer.data(kf + 110 * ncomps + c);
        const auto *kf_111 = buffer.data(kf + 111 * ncomps + c);
        const auto *kf_112 = buffer.data(kf + 112 * ncomps + c);
        const auto *kf_113 = buffer.data(kf + 113 * ncomps + c);
        const auto *kf_114 = buffer.data(kf + 114 * ncomps + c);
        const auto *kf_115 = buffer.data(kf + 115 * ncomps + c);
        const auto *kf_116 = buffer.data(kf + 116 * ncomps + c);
        const auto *kf_117 = buffer.data(kf + 117 * ncomps + c);
        const auto *kf_118 = buffer.data(kf + 118 * ncomps + c);
        const auto *kf_119 = buffer.data(kf + 119 * ncomps + c);
        const auto *kf_120 = buffer.data(kf + 120 * ncomps + c);
        const auto *kf_121 = buffer.data(kf + 121 * ncomps + c);
        const auto *kf_122 = buffer.data(kf + 122 * ncomps + c);
        const auto *kf_123 = buffer.data(kf + 123 * ncomps + c);
        const auto *kf_124 = buffer.data(kf + 124 * ncomps + c);
        const auto *kf_125 = buffer.data(kf + 125 * ncomps + c);
        const auto *kf_126 = buffer.data(kf + 126 * ncomps + c);
        const auto *kf_127 = buffer.data(kf + 127 * ncomps + c);
        const auto *kf_128 = buffer.data(kf + 128 * ncomps + c);
        const auto *kf_129 = buffer.data(kf + 129 * ncomps + c);
        const auto *kf_130 = buffer.data(kf + 130 * ncomps + c);
        const auto *kf_131 = buffer.data(kf + 131 * ncomps + c);
        const auto *kf_132 = buffer.data(kf + 132 * ncomps + c);
        const auto *kf_133 = buffer.data(kf + 133 * ncomps + c);
        const auto *kf_134 = buffer.data(kf + 134 * ncomps + c);
        const auto *kf_135 = buffer.data(kf + 135 * ncomps + c);
        const auto *kf_136 = buffer.data(kf + 136 * ncomps + c);
        const auto *kf_137 = buffer.data(kf + 137 * ncomps + c);
        const auto *kf_138 = buffer.data(kf + 138 * ncomps + c);
        const auto *kf_139 = buffer.data(kf + 139 * ncomps + c);
        const auto *kf_140 = buffer.data(kf + 140 * ncomps + c);
        const auto *kf_141 = buffer.data(kf + 141 * ncomps + c);
        const auto *kf_142 = buffer.data(kf + 142 * ncomps + c);
        const auto *kf_143 = buffer.data(kf + 143 * ncomps + c);
        const auto *kf_144 = buffer.data(kf + 144 * ncomps + c);
        const auto *kf_145 = buffer.data(kf + 145 * ncomps + c);
        const auto *kf_146 = buffer.data(kf + 146 * ncomps + c);
        const auto *kf_147 = buffer.data(kf + 147 * ncomps + c);
        const auto *kf_148 = buffer.data(kf + 148 * ncomps + c);
        const auto *kf_149 = buffer.data(kf + 149 * ncomps + c);
        const auto *kf_150 = buffer.data(kf + 150 * ncomps + c);
        const auto *kf_151 = buffer.data(kf + 151 * ncomps + c);
        const auto *kf_152 = buffer.data(kf + 152 * ncomps + c);
        const auto *kf_153 = buffer.data(kf + 153 * ncomps + c);
        const auto *kf_154 = buffer.data(kf + 154 * ncomps + c);
        const auto *kf_155 = buffer.data(kf + 155 * ncomps + c);
        const auto *kf_156 = buffer.data(kf + 156 * ncomps + c);
        const auto *kf_157 = buffer.data(kf + 157 * ncomps + c);
        const auto *kf_158 = buffer.data(kf + 158 * ncomps + c);
        const auto *kf_159 = buffer.data(kf + 159 * ncomps + c);
        const auto *kf_160 = buffer.data(kf + 160 * ncomps + c);
        const auto *kf_161 = buffer.data(kf + 161 * ncomps + c);
        const auto *kf_162 = buffer.data(kf + 162 * ncomps + c);
        const auto *kf_163 = buffer.data(kf + 163 * ncomps + c);
        const auto *kf_164 = buffer.data(kf + 164 * ncomps + c);
        const auto *kf_165 = buffer.data(kf + 165 * ncomps + c);
        const auto *kf_166 = buffer.data(kf + 166 * ncomps + c);
        const auto *kf_167 = buffer.data(kf + 167 * ncomps + c);
        const auto *kf_168 = buffer.data(kf + 168 * ncomps + c);
        const auto *kf_169 = buffer.data(kf + 169 * ncomps + c);
        const auto *kf_170 = buffer.data(kf + 170 * ncomps + c);
        const auto *kf_171 = buffer.data(kf + 171 * ncomps + c);
        const auto *kf_172 = buffer.data(kf + 172 * ncomps + c);
        const auto *kf_173 = buffer.data(kf + 173 * ncomps + c);
        const auto *kf_174 = buffer.data(kf + 174 * ncomps + c);
        const auto *kf_175 = buffer.data(kf + 175 * ncomps + c);
        const auto *kf_176 = buffer.data(kf + 176 * ncomps + c);
        const auto *kf_177 = buffer.data(kf + 177 * ncomps + c);
        const auto *kf_178 = buffer.data(kf + 178 * ncomps + c);
        const auto *kf_179 = buffer.data(kf + 179 * ncomps + c);
        const auto *kf_180 = buffer.data(kf + 180 * ncomps + c);
        const auto *kf_181 = buffer.data(kf + 181 * ncomps + c);
        const auto *kf_182 = buffer.data(kf + 182 * ncomps + c);
        const auto *kf_183 = buffer.data(kf + 183 * ncomps + c);
        const auto *kf_184 = buffer.data(kf + 184 * ncomps + c);
        const auto *kf_185 = buffer.data(kf + 185 * ncomps + c);
        const auto *kf_186 = buffer.data(kf + 186 * ncomps + c);
        const auto *kf_187 = buffer.data(kf + 187 * ncomps + c);
        const auto *kf_188 = buffer.data(kf + 188 * ncomps + c);
        const auto *kf_189 = buffer.data(kf + 189 * ncomps + c);
        const auto *kf_190 = buffer.data(kf + 190 * ncomps + c);
        const auto *kf_191 = buffer.data(kf + 191 * ncomps + c);
        const auto *kf_192 = buffer.data(kf + 192 * ncomps + c);
        const auto *kf_193 = buffer.data(kf + 193 * ncomps + c);
        const auto *kf_194 = buffer.data(kf + 194 * ncomps + c);
        const auto *kf_196 = buffer.data(kf + 196 * ncomps + c);
        const auto *kf_197 = buffer.data(kf + 197 * ncomps + c);
        const auto *kf_198 = buffer.data(kf + 198 * ncomps + c);
        const auto *kf_199 = buffer.data(kf + 199 * ncomps + c);
        const auto *kf_209 = buffer.data(kf + 209 * ncomps + c);
        const auto *kf_216 = buffer.data(kf + 216 * ncomps + c);
        const auto *kf_217 = buffer.data(kf + 217 * ncomps + c);
        const auto *kf_218 = buffer.data(kf + 218 * ncomps + c);
        const auto *kf_219 = buffer.data(kf + 219 * ncomps + c);
        const auto *kf_226 = buffer.data(kf + 226 * ncomps + c);
        const auto *kf_227 = buffer.data(kf + 227 * ncomps + c);
        const auto *kf_228 = buffer.data(kf + 228 * ncomps + c);
        const auto *kf_229 = buffer.data(kf + 229 * ncomps + c);
        const auto *kf_236 = buffer.data(kf + 236 * ncomps + c);
        const auto *kf_237 = buffer.data(kf + 237 * ncomps + c);
        const auto *kf_238 = buffer.data(kf + 238 * ncomps + c);
        const auto *kf_239 = buffer.data(kf + 239 * ncomps + c);
        const auto *kf_246 = buffer.data(kf + 246 * ncomps + c);
        const auto *kf_247 = buffer.data(kf + 247 * ncomps + c);
        const auto *kf_248 = buffer.data(kf + 248 * ncomps + c);
        const auto *kf_249 = buffer.data(kf + 249 * ncomps + c);
        const auto *kf_259 = buffer.data(kf + 259 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_y, ab_z, if__96, if__97, \
                         if__98, if__99, kf_136, kf_137, kf_138, kf_139, \
                         kf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_y[k] * if__96[k]
                       + kf_136[k];

            t_146[k] = ab_y[k] * if__97[k]
                       + kf_137[k];

            t_147[k] = ab_y[k] * if__98[k]
                       + kf_138[k];

            t_148[k] = ab_y[k] * if__99[k]
                       + kf_139[k];

            t_149[k] = ab_z[k] * if__99[k]
                       + kf_149[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, if__100, if__101, if__102, \
                         if__103, if__104, kf_100, kf_101, kf_102, kf_103, \
                         kf_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * if__100[k]
                       + kf_100[k];

            t_151[k] = ab_x[k] * if__101[k]
                       + kf_101[k];

            t_152[k] = ab_x[k] * if__102[k]
                       + kf_102[k];

            t_153[k] = ab_x[k] * if__103[k]
                       + kf_103[k];

            t_154[k] = ab_x[k] * if__104[k]
                       + kf_104[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, if__105, if__106, if__107, \
                         if__108, if__109, kf_105, kf_106, kf_107, kf_108, \
                         kf_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * if__105[k]
                       + kf_105[k];

            t_156[k] = ab_x[k] * if__106[k]
                       + kf_106[k];

            t_157[k] = ab_x[k] * if__107[k]
                       + kf_107[k];

            t_158[k] = ab_x[k] * if__108[k]
                       + kf_108[k];

            t_159[k] = ab_x[k] * if__109[k]
                       + kf_109[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_y, ab_z, if__106, if__107, \
                         if__108, if__109, kf_156, kf_157, kf_158, kf_159, \
                         kf_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_y[k] * if__106[k]
                       + kf_156[k];

            t_161[k] = ab_y[k] * if__107[k]
                       + kf_157[k];

            t_162[k] = ab_y[k] * if__108[k]
                       + kf_158[k];

            t_163[k] = ab_y[k] * if__109[k]
                       + kf_159[k];

            t_164[k] = ab_z[k] * if__109[k]
                       + kf_169[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, if__110, if__111, if__112, \
                         if__113, if__114, kf_110, kf_111, kf_112, kf_113, \
                         kf_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * if__110[k]
                       + kf_110[k];

            t_166[k] = ab_x[k] * if__111[k]
                       + kf_111[k];

            t_167[k] = ab_x[k] * if__112[k]
                       + kf_112[k];

            t_168[k] = ab_x[k] * if__113[k]
                       + kf_113[k];

            t_169[k] = ab_x[k] * if__114[k]
                       + kf_114[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, if__115, if__116, if__117, \
                         if__118, if__119, kf_115, kf_116, kf_117, kf_118, \
                         kf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = ab_x[k] * if__115[k]
                       + kf_115[k];

            t_171[k] = ab_x[k] * if__116[k]
                       + kf_116[k];

            t_172[k] = ab_x[k] * if__117[k]
                       + kf_117[k];

            t_173[k] = ab_x[k] * if__118[k]
                       + kf_118[k];

            t_174[k] = ab_x[k] * if__119[k]
                       + kf_119[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_y, ab_z, if__116, if__117, \
                         if__118, if__119, kf_166, kf_167, kf_168, kf_169, \
                         kf_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_y[k] * if__116[k]
                       + kf_166[k];

            t_176[k] = ab_y[k] * if__117[k]
                       + kf_167[k];

            t_177[k] = ab_y[k] * if__118[k]
                       + kf_168[k];

            t_178[k] = ab_y[k] * if__119[k]
                       + kf_169[k];

            t_179[k] = ab_z[k] * if__119[k]
                       + kf_179[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, if__120, if__121, if__122, \
                         if__123, if__124, kf_120, kf_121, kf_122, kf_123, \
                         kf_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * if__120[k]
                       + kf_120[k];

            t_181[k] = ab_x[k] * if__121[k]
                       + kf_121[k];

            t_182[k] = ab_x[k] * if__122[k]
                       + kf_122[k];

            t_183[k] = ab_x[k] * if__123[k]
                       + kf_123[k];

            t_184[k] = ab_x[k] * if__124[k]
                       + kf_124[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, if__125, if__126, if__127, \
                         if__128, if__129, kf_125, kf_126, kf_127, kf_128, \
                         kf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * if__125[k]
                       + kf_125[k];

            t_186[k] = ab_x[k] * if__126[k]
                       + kf_126[k];

            t_187[k] = ab_x[k] * if__127[k]
                       + kf_127[k];

            t_188[k] = ab_x[k] * if__128[k]
                       + kf_128[k];

            t_189[k] = ab_x[k] * if__129[k]
                       + kf_129[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_y, ab_z, if__126, if__127, \
                         if__128, if__129, kf_176, kf_177, kf_178, kf_179, \
                         kf_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_y[k] * if__126[k]
                       + kf_176[k];

            t_191[k] = ab_y[k] * if__127[k]
                       + kf_177[k];

            t_192[k] = ab_y[k] * if__128[k]
                       + kf_178[k];

            t_193[k] = ab_y[k] * if__129[k]
                       + kf_179[k];

            t_194[k] = ab_z[k] * if__129[k]
                       + kf_189[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, if__130, if__131, if__132, \
                         if__133, if__134, kf_130, kf_131, kf_132, kf_133, \
                         kf_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * if__130[k]
                       + kf_130[k];

            t_196[k] = ab_x[k] * if__131[k]
                       + kf_131[k];

            t_197[k] = ab_x[k] * if__132[k]
                       + kf_132[k];

            t_198[k] = ab_x[k] * if__133[k]
                       + kf_133[k];

            t_199[k] = ab_x[k] * if__134[k]
                       + kf_134[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, if__135, if__136, if__137, \
                         if__138, if__139, kf_135, kf_136, kf_137, kf_138, \
                         kf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = ab_x[k] * if__135[k]
                       + kf_135[k];

            t_201[k] = ab_x[k] * if__136[k]
                       + kf_136[k];

            t_202[k] = ab_x[k] * if__137[k]
                       + kf_137[k];

            t_203[k] = ab_x[k] * if__138[k]
                       + kf_138[k];

            t_204[k] = ab_x[k] * if__139[k]
                       + kf_139[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_y, ab_z, if__136, if__137, \
                         if__138, if__139, kf_186, kf_187, kf_188, kf_189, \
                         kf_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_y[k] * if__136[k]
                       + kf_186[k];

            t_206[k] = ab_y[k] * if__137[k]
                       + kf_187[k];

            t_207[k] = ab_y[k] * if__138[k]
                       + kf_188[k];

            t_208[k] = ab_y[k] * if__139[k]
                       + kf_189[k];

            t_209[k] = ab_z[k] * if__139[k]
                       + kf_199[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, if__140, if__141, if__142, \
                         if__143, if__144, kf_140, kf_141, kf_142, kf_143, \
                         kf_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * if__140[k]
                       + kf_140[k];

            t_211[k] = ab_x[k] * if__141[k]
                       + kf_141[k];

            t_212[k] = ab_x[k] * if__142[k]
                       + kf_142[k];

            t_213[k] = ab_x[k] * if__143[k]
                       + kf_143[k];

            t_214[k] = ab_x[k] * if__144[k]
                       + kf_144[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, if__145, if__146, if__147, \
                         if__148, if__149, kf_145, kf_146, kf_147, kf_148, \
                         kf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * if__145[k]
                       + kf_145[k];

            t_216[k] = ab_x[k] * if__146[k]
                       + kf_146[k];

            t_217[k] = ab_x[k] * if__147[k]
                       + kf_147[k];

            t_218[k] = ab_x[k] * if__148[k]
                       + kf_148[k];

            t_219[k] = ab_x[k] * if__149[k]
                       + kf_149[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_y, ab_z, if__146, if__147, \
                         if__148, if__149, kf_196, kf_197, kf_198, kf_199, \
                         kf_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_y[k] * if__146[k]
                       + kf_196[k];

            t_221[k] = ab_y[k] * if__147[k]
                       + kf_197[k];

            t_222[k] = ab_y[k] * if__148[k]
                       + kf_198[k];

            t_223[k] = ab_y[k] * if__149[k]
                       + kf_199[k];

            t_224[k] = ab_z[k] * if__149[k]
                       + kf_209[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, if__150, if__151, if__152, \
                         if__153, if__154, kf_150, kf_151, kf_152, kf_153, \
                         kf_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_x[k] * if__150[k]
                       + kf_150[k];

            t_226[k] = ab_x[k] * if__151[k]
                       + kf_151[k];

            t_227[k] = ab_x[k] * if__152[k]
                       + kf_152[k];

            t_228[k] = ab_x[k] * if__153[k]
                       + kf_153[k];

            t_229[k] = ab_x[k] * if__154[k]
                       + kf_154[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, if__155, if__156, if__157, \
                         if__158, if__159, kf_155, kf_156, kf_157, kf_158, \
                         kf_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_x[k] * if__155[k]
                       + kf_155[k];

            t_231[k] = ab_x[k] * if__156[k]
                       + kf_156[k];

            t_232[k] = ab_x[k] * if__157[k]
                       + kf_157[k];

            t_233[k] = ab_x[k] * if__158[k]
                       + kf_158[k];

            t_234[k] = ab_x[k] * if__159[k]
                       + kf_159[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_y, ab_z, if__156, if__157, \
                         if__158, if__159, kf_216, kf_217, kf_218, kf_219, \
                         kf_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = ab_y[k] * if__156[k]
                       + kf_216[k];

            t_236[k] = ab_y[k] * if__157[k]
                       + kf_217[k];

            t_237[k] = ab_y[k] * if__158[k]
                       + kf_218[k];

            t_238[k] = ab_y[k] * if__159[k]
                       + kf_219[k];

            t_239[k] = ab_z[k] * if__159[k]
                       + kf_229[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, if__160, if__161, if__162, \
                         if__163, if__164, kf_160, kf_161, kf_162, kf_163, \
                         kf_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = ab_x[k] * if__160[k]
                       + kf_160[k];

            t_241[k] = ab_x[k] * if__161[k]
                       + kf_161[k];

            t_242[k] = ab_x[k] * if__162[k]
                       + kf_162[k];

            t_243[k] = ab_x[k] * if__163[k]
                       + kf_163[k];

            t_244[k] = ab_x[k] * if__164[k]
                       + kf_164[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, if__165, if__166, if__167, \
                         if__168, if__169, kf_165, kf_166, kf_167, kf_168, \
                         kf_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = ab_x[k] * if__165[k]
                       + kf_165[k];

            t_246[k] = ab_x[k] * if__166[k]
                       + kf_166[k];

            t_247[k] = ab_x[k] * if__167[k]
                       + kf_167[k];

            t_248[k] = ab_x[k] * if__168[k]
                       + kf_168[k];

            t_249[k] = ab_x[k] * if__169[k]
                       + kf_169[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_y, ab_z, if__166, if__167, \
                         if__168, if__169, kf_226, kf_227, kf_228, kf_229, \
                         kf_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = ab_y[k] * if__166[k]
                       + kf_226[k];

            t_251[k] = ab_y[k] * if__167[k]
                       + kf_227[k];

            t_252[k] = ab_y[k] * if__168[k]
                       + kf_228[k];

            t_253[k] = ab_y[k] * if__169[k]
                       + kf_229[k];

            t_254[k] = ab_z[k] * if__169[k]
                       + kf_239[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, if__170, if__171, if__172, \
                         if__173, if__174, kf_170, kf_171, kf_172, kf_173, \
                         kf_174 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = ab_x[k] * if__170[k]
                       + kf_170[k];

            t_256[k] = ab_x[k] * if__171[k]
                       + kf_171[k];

            t_257[k] = ab_x[k] * if__172[k]
                       + kf_172[k];

            t_258[k] = ab_x[k] * if__173[k]
                       + kf_173[k];

            t_259[k] = ab_x[k] * if__174[k]
                       + kf_174[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, if__175, if__176, if__177, \
                         if__178, if__179, kf_175, kf_176, kf_177, kf_178, \
                         kf_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = ab_x[k] * if__175[k]
                       + kf_175[k];

            t_261[k] = ab_x[k] * if__176[k]
                       + kf_176[k];

            t_262[k] = ab_x[k] * if__177[k]
                       + kf_177[k];

            t_263[k] = ab_x[k] * if__178[k]
                       + kf_178[k];

            t_264[k] = ab_x[k] * if__179[k]
                       + kf_179[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_y, ab_z, if__176, if__177, \
                         if__178, if__179, kf_236, kf_237, kf_238, kf_239, \
                         kf_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = ab_y[k] * if__176[k]
                       + kf_236[k];

            t_266[k] = ab_y[k] * if__177[k]
                       + kf_237[k];

            t_267[k] = ab_y[k] * if__178[k]
                       + kf_238[k];

            t_268[k] = ab_y[k] * if__179[k]
                       + kf_239[k];

            t_269[k] = ab_z[k] * if__179[k]
                       + kf_249[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, if__180, if__181, if__182, \
                         if__183, if__184, kf_180, kf_181, kf_182, kf_183, \
                         kf_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = ab_x[k] * if__180[k]
                       + kf_180[k];

            t_271[k] = ab_x[k] * if__181[k]
                       + kf_181[k];

            t_272[k] = ab_x[k] * if__182[k]
                       + kf_182[k];

            t_273[k] = ab_x[k] * if__183[k]
                       + kf_183[k];

            t_274[k] = ab_x[k] * if__184[k]
                       + kf_184[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, if__185, if__186, if__187, \
                         if__188, if__189, kf_185, kf_186, kf_187, kf_188, \
                         kf_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = ab_x[k] * if__185[k]
                       + kf_185[k];

            t_276[k] = ab_x[k] * if__186[k]
                       + kf_186[k];

            t_277[k] = ab_x[k] * if__187[k]
                       + kf_187[k];

            t_278[k] = ab_x[k] * if__188[k]
                       + kf_188[k];

            t_279[k] = ab_x[k] * if__189[k]
                       + kf_189[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_y, ab_z, if__186, if__187, \
                         if__188, if__189, kf_246, kf_247, kf_248, kf_249, \
                         kf_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = ab_y[k] * if__186[k]
                       + kf_246[k];

            t_281[k] = ab_y[k] * if__187[k]
                       + kf_247[k];

            t_282[k] = ab_y[k] * if__188[k]
                       + kf_248[k];

            t_283[k] = ab_y[k] * if__189[k]
                       + kf_249[k];

            t_284[k] = ab_z[k] * if__189[k]
                       + kf_259[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, if__190, if__191, if__192, \
                         if__193, if__194, kf_190, kf_191, kf_192, kf_193, \
                         kf_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = ab_x[k] * if__190[k]
                       + kf_190[k];

            t_286[k] = ab_x[k] * if__191[k]
                       + kf_191[k];

            t_287[k] = ab_x[k] * if__192[k]
                       + kf_192[k];

            t_288[k] = ab_x[k] * if__193[k]
                       + kf_193[k];

            t_289[k] = ab_x[k] * if__194[k]
                       + kf_194[k];
        }
    }
}

static auto
compute_hrr_ig_out_of_first_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t if_, const size_t kf,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *if__195 = buffer.data(if_ + 195 * ncomps + c);
        const auto *if__196 = buffer.data(if_ + 196 * ncomps + c);
        const auto *if__197 = buffer.data(if_ + 197 * ncomps + c);
        const auto *if__198 = buffer.data(if_ + 198 * ncomps + c);
        const auto *if__199 = buffer.data(if_ + 199 * ncomps + c);
        const auto *if__200 = buffer.data(if_ + 200 * ncomps + c);
        const auto *if__201 = buffer.data(if_ + 201 * ncomps + c);
        const auto *if__202 = buffer.data(if_ + 202 * ncomps + c);
        const auto *if__203 = buffer.data(if_ + 203 * ncomps + c);
        const auto *if__204 = buffer.data(if_ + 204 * ncomps + c);
        const auto *if__205 = buffer.data(if_ + 205 * ncomps + c);
        const auto *if__206 = buffer.data(if_ + 206 * ncomps + c);
        const auto *if__207 = buffer.data(if_ + 207 * ncomps + c);
        const auto *if__208 = buffer.data(if_ + 208 * ncomps + c);
        const auto *if__209 = buffer.data(if_ + 209 * ncomps + c);
        const auto *if__210 = buffer.data(if_ + 210 * ncomps + c);
        const auto *if__211 = buffer.data(if_ + 211 * ncomps + c);
        const auto *if__212 = buffer.data(if_ + 212 * ncomps + c);
        const auto *if__213 = buffer.data(if_ + 213 * ncomps + c);
        const auto *if__214 = buffer.data(if_ + 214 * ncomps + c);
        const auto *if__215 = buffer.data(if_ + 215 * ncomps + c);
        const auto *if__216 = buffer.data(if_ + 216 * ncomps + c);
        const auto *if__217 = buffer.data(if_ + 217 * ncomps + c);
        const auto *if__218 = buffer.data(if_ + 218 * ncomps + c);
        const auto *if__219 = buffer.data(if_ + 219 * ncomps + c);
        const auto *if__220 = buffer.data(if_ + 220 * ncomps + c);
        const auto *if__221 = buffer.data(if_ + 221 * ncomps + c);
        const auto *if__222 = buffer.data(if_ + 222 * ncomps + c);
        const auto *if__223 = buffer.data(if_ + 223 * ncomps + c);
        const auto *if__224 = buffer.data(if_ + 224 * ncomps + c);
        const auto *if__225 = buffer.data(if_ + 225 * ncomps + c);
        const auto *if__226 = buffer.data(if_ + 226 * ncomps + c);
        const auto *if__227 = buffer.data(if_ + 227 * ncomps + c);
        const auto *if__228 = buffer.data(if_ + 228 * ncomps + c);
        const auto *if__229 = buffer.data(if_ + 229 * ncomps + c);
        const auto *if__230 = buffer.data(if_ + 230 * ncomps + c);
        const auto *if__231 = buffer.data(if_ + 231 * ncomps + c);
        const auto *if__232 = buffer.data(if_ + 232 * ncomps + c);
        const auto *if__233 = buffer.data(if_ + 233 * ncomps + c);
        const auto *if__234 = buffer.data(if_ + 234 * ncomps + c);
        const auto *if__235 = buffer.data(if_ + 235 * ncomps + c);
        const auto *if__236 = buffer.data(if_ + 236 * ncomps + c);
        const auto *if__237 = buffer.data(if_ + 237 * ncomps + c);
        const auto *if__238 = buffer.data(if_ + 238 * ncomps + c);
        const auto *if__239 = buffer.data(if_ + 239 * ncomps + c);
        const auto *if__240 = buffer.data(if_ + 240 * ncomps + c);
        const auto *if__241 = buffer.data(if_ + 241 * ncomps + c);
        const auto *if__242 = buffer.data(if_ + 242 * ncomps + c);
        const auto *if__243 = buffer.data(if_ + 243 * ncomps + c);
        const auto *if__244 = buffer.data(if_ + 244 * ncomps + c);
        const auto *if__245 = buffer.data(if_ + 245 * ncomps + c);
        const auto *if__246 = buffer.data(if_ + 246 * ncomps + c);
        const auto *if__247 = buffer.data(if_ + 247 * ncomps + c);
        const auto *if__248 = buffer.data(if_ + 248 * ncomps + c);
        const auto *if__249 = buffer.data(if_ + 249 * ncomps + c);
        const auto *if__250 = buffer.data(if_ + 250 * ncomps + c);
        const auto *if__251 = buffer.data(if_ + 251 * ncomps + c);
        const auto *if__252 = buffer.data(if_ + 252 * ncomps + c);
        const auto *if__253 = buffer.data(if_ + 253 * ncomps + c);
        const auto *if__254 = buffer.data(if_ + 254 * ncomps + c);
        const auto *if__255 = buffer.data(if_ + 255 * ncomps + c);
        const auto *if__256 = buffer.data(if_ + 256 * ncomps + c);
        const auto *if__257 = buffer.data(if_ + 257 * ncomps + c);
        const auto *if__258 = buffer.data(if_ + 258 * ncomps + c);
        const auto *if__259 = buffer.data(if_ + 259 * ncomps + c);
        const auto *if__260 = buffer.data(if_ + 260 * ncomps + c);
        const auto *if__261 = buffer.data(if_ + 261 * ncomps + c);
        const auto *if__262 = buffer.data(if_ + 262 * ncomps + c);
        const auto *if__263 = buffer.data(if_ + 263 * ncomps + c);
        const auto *if__264 = buffer.data(if_ + 264 * ncomps + c);
        const auto *if__265 = buffer.data(if_ + 265 * ncomps + c);
        const auto *if__266 = buffer.data(if_ + 266 * ncomps + c);
        const auto *if__267 = buffer.data(if_ + 267 * ncomps + c);
        const auto *if__268 = buffer.data(if_ + 268 * ncomps + c);
        const auto *if__269 = buffer.data(if_ + 269 * ncomps + c);
        const auto *if__270 = buffer.data(if_ + 270 * ncomps + c);
        const auto *if__271 = buffer.data(if_ + 271 * ncomps + c);
        const auto *if__272 = buffer.data(if_ + 272 * ncomps + c);
        const auto *if__273 = buffer.data(if_ + 273 * ncomps + c);
        const auto *if__274 = buffer.data(if_ + 274 * ncomps + c);
        const auto *if__275 = buffer.data(if_ + 275 * ncomps + c);
        const auto *if__276 = buffer.data(if_ + 276 * ncomps + c);
        const auto *if__277 = buffer.data(if_ + 277 * ncomps + c);
        const auto *if__278 = buffer.data(if_ + 278 * ncomps + c);
        const auto *if__279 = buffer.data(if_ + 279 * ncomps + c);

        const auto *kf_195 = buffer.data(kf + 195 * ncomps + c);
        const auto *kf_196 = buffer.data(kf + 196 * ncomps + c);
        const auto *kf_197 = buffer.data(kf + 197 * ncomps + c);
        const auto *kf_198 = buffer.data(kf + 198 * ncomps + c);
        const auto *kf_199 = buffer.data(kf + 199 * ncomps + c);
        const auto *kf_200 = buffer.data(kf + 200 * ncomps + c);
        const auto *kf_201 = buffer.data(kf + 201 * ncomps + c);
        const auto *kf_202 = buffer.data(kf + 202 * ncomps + c);
        const auto *kf_203 = buffer.data(kf + 203 * ncomps + c);
        const auto *kf_204 = buffer.data(kf + 204 * ncomps + c);
        const auto *kf_205 = buffer.data(kf + 205 * ncomps + c);
        const auto *kf_206 = buffer.data(kf + 206 * ncomps + c);
        const auto *kf_207 = buffer.data(kf + 207 * ncomps + c);
        const auto *kf_208 = buffer.data(kf + 208 * ncomps + c);
        const auto *kf_209 = buffer.data(kf + 209 * ncomps + c);
        const auto *kf_210 = buffer.data(kf + 210 * ncomps + c);
        const auto *kf_211 = buffer.data(kf + 211 * ncomps + c);
        const auto *kf_212 = buffer.data(kf + 212 * ncomps + c);
        const auto *kf_213 = buffer.data(kf + 213 * ncomps + c);
        const auto *kf_214 = buffer.data(kf + 214 * ncomps + c);
        const auto *kf_215 = buffer.data(kf + 215 * ncomps + c);
        const auto *kf_216 = buffer.data(kf + 216 * ncomps + c);
        const auto *kf_217 = buffer.data(kf + 217 * ncomps + c);
        const auto *kf_218 = buffer.data(kf + 218 * ncomps + c);
        const auto *kf_219 = buffer.data(kf + 219 * ncomps + c);
        const auto *kf_220 = buffer.data(kf + 220 * ncomps + c);
        const auto *kf_221 = buffer.data(kf + 221 * ncomps + c);
        const auto *kf_222 = buffer.data(kf + 222 * ncomps + c);
        const auto *kf_223 = buffer.data(kf + 223 * ncomps + c);
        const auto *kf_224 = buffer.data(kf + 224 * ncomps + c);
        const auto *kf_225 = buffer.data(kf + 225 * ncomps + c);
        const auto *kf_226 = buffer.data(kf + 226 * ncomps + c);
        const auto *kf_227 = buffer.data(kf + 227 * ncomps + c);
        const auto *kf_228 = buffer.data(kf + 228 * ncomps + c);
        const auto *kf_229 = buffer.data(kf + 229 * ncomps + c);
        const auto *kf_230 = buffer.data(kf + 230 * ncomps + c);
        const auto *kf_231 = buffer.data(kf + 231 * ncomps + c);
        const auto *kf_232 = buffer.data(kf + 232 * ncomps + c);
        const auto *kf_233 = buffer.data(kf + 233 * ncomps + c);
        const auto *kf_234 = buffer.data(kf + 234 * ncomps + c);
        const auto *kf_235 = buffer.data(kf + 235 * ncomps + c);
        const auto *kf_236 = buffer.data(kf + 236 * ncomps + c);
        const auto *kf_237 = buffer.data(kf + 237 * ncomps + c);
        const auto *kf_238 = buffer.data(kf + 238 * ncomps + c);
        const auto *kf_239 = buffer.data(kf + 239 * ncomps + c);
        const auto *kf_240 = buffer.data(kf + 240 * ncomps + c);
        const auto *kf_241 = buffer.data(kf + 241 * ncomps + c);
        const auto *kf_242 = buffer.data(kf + 242 * ncomps + c);
        const auto *kf_243 = buffer.data(kf + 243 * ncomps + c);
        const auto *kf_244 = buffer.data(kf + 244 * ncomps + c);
        const auto *kf_245 = buffer.data(kf + 245 * ncomps + c);
        const auto *kf_246 = buffer.data(kf + 246 * ncomps + c);
        const auto *kf_247 = buffer.data(kf + 247 * ncomps + c);
        const auto *kf_248 = buffer.data(kf + 248 * ncomps + c);
        const auto *kf_249 = buffer.data(kf + 249 * ncomps + c);
        const auto *kf_250 = buffer.data(kf + 250 * ncomps + c);
        const auto *kf_251 = buffer.data(kf + 251 * ncomps + c);
        const auto *kf_252 = buffer.data(kf + 252 * ncomps + c);
        const auto *kf_253 = buffer.data(kf + 253 * ncomps + c);
        const auto *kf_254 = buffer.data(kf + 254 * ncomps + c);
        const auto *kf_255 = buffer.data(kf + 255 * ncomps + c);
        const auto *kf_256 = buffer.data(kf + 256 * ncomps + c);
        const auto *kf_257 = buffer.data(kf + 257 * ncomps + c);
        const auto *kf_258 = buffer.data(kf + 258 * ncomps + c);
        const auto *kf_259 = buffer.data(kf + 259 * ncomps + c);
        const auto *kf_260 = buffer.data(kf + 260 * ncomps + c);
        const auto *kf_261 = buffer.data(kf + 261 * ncomps + c);
        const auto *kf_262 = buffer.data(kf + 262 * ncomps + c);
        const auto *kf_263 = buffer.data(kf + 263 * ncomps + c);
        const auto *kf_264 = buffer.data(kf + 264 * ncomps + c);
        const auto *kf_265 = buffer.data(kf + 265 * ncomps + c);
        const auto *kf_266 = buffer.data(kf + 266 * ncomps + c);
        const auto *kf_267 = buffer.data(kf + 267 * ncomps + c);
        const auto *kf_268 = buffer.data(kf + 268 * ncomps + c);
        const auto *kf_269 = buffer.data(kf + 269 * ncomps + c);
        const auto *kf_270 = buffer.data(kf + 270 * ncomps + c);
        const auto *kf_271 = buffer.data(kf + 271 * ncomps + c);
        const auto *kf_272 = buffer.data(kf + 272 * ncomps + c);
        const auto *kf_273 = buffer.data(kf + 273 * ncomps + c);
        const auto *kf_274 = buffer.data(kf + 274 * ncomps + c);
        const auto *kf_275 = buffer.data(kf + 275 * ncomps + c);
        const auto *kf_276 = buffer.data(kf + 276 * ncomps + c);
        const auto *kf_277 = buffer.data(kf + 277 * ncomps + c);
        const auto *kf_278 = buffer.data(kf + 278 * ncomps + c);
        const auto *kf_279 = buffer.data(kf + 279 * ncomps + c);
        const auto *kf_286 = buffer.data(kf + 286 * ncomps + c);
        const auto *kf_287 = buffer.data(kf + 287 * ncomps + c);
        const auto *kf_288 = buffer.data(kf + 288 * ncomps + c);
        const auto *kf_289 = buffer.data(kf + 289 * ncomps + c);
        const auto *kf_296 = buffer.data(kf + 296 * ncomps + c);
        const auto *kf_297 = buffer.data(kf + 297 * ncomps + c);
        const auto *kf_298 = buffer.data(kf + 298 * ncomps + c);
        const auto *kf_299 = buffer.data(kf + 299 * ncomps + c);
        const auto *kf_306 = buffer.data(kf + 306 * ncomps + c);
        const auto *kf_307 = buffer.data(kf + 307 * ncomps + c);
        const auto *kf_308 = buffer.data(kf + 308 * ncomps + c);
        const auto *kf_309 = buffer.data(kf + 309 * ncomps + c);
        const auto *kf_316 = buffer.data(kf + 316 * ncomps + c);
        const auto *kf_317 = buffer.data(kf + 317 * ncomps + c);
        const auto *kf_318 = buffer.data(kf + 318 * ncomps + c);
        const auto *kf_319 = buffer.data(kf + 319 * ncomps + c);
        const auto *kf_326 = buffer.data(kf + 326 * ncomps + c);
        const auto *kf_327 = buffer.data(kf + 327 * ncomps + c);
        const auto *kf_328 = buffer.data(kf + 328 * ncomps + c);
        const auto *kf_329 = buffer.data(kf + 329 * ncomps + c);
        const auto *kf_336 = buffer.data(kf + 336 * ncomps + c);
        const auto *kf_337 = buffer.data(kf + 337 * ncomps + c);
        const auto *kf_338 = buffer.data(kf + 338 * ncomps + c);
        const auto *kf_339 = buffer.data(kf + 339 * ncomps + c);
        const auto *kf_346 = buffer.data(kf + 346 * ncomps + c);
        const auto *kf_347 = buffer.data(kf + 347 * ncomps + c);
        const auto *kf_348 = buffer.data(kf + 348 * ncomps + c);
        const auto *kf_349 = buffer.data(kf + 349 * ncomps + c);
        const auto *kf_359 = buffer.data(kf + 359 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, if__195, if__196, if__197, \
                         if__198, if__199, kf_195, kf_196, kf_197, kf_198, \
                         kf_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_x[k] * if__195[k]
                       + kf_195[k];

            t_291[k] = ab_x[k] * if__196[k]
                       + kf_196[k];

            t_292[k] = ab_x[k] * if__197[k]
                       + kf_197[k];

            t_293[k] = ab_x[k] * if__198[k]
                       + kf_198[k];

            t_294[k] = ab_x[k] * if__199[k]
                       + kf_199[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_y, ab_z, if__196, if__197, \
                         if__198, if__199, kf_256, kf_257, kf_258, kf_259, \
                         kf_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = ab_y[k] * if__196[k]
                       + kf_256[k];

            t_296[k] = ab_y[k] * if__197[k]
                       + kf_257[k];

            t_297[k] = ab_y[k] * if__198[k]
                       + kf_258[k];

            t_298[k] = ab_y[k] * if__199[k]
                       + kf_259[k];

            t_299[k] = ab_z[k] * if__199[k]
                       + kf_269[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, if__200, if__201, if__202, \
                         if__203, if__204, kf_200, kf_201, kf_202, kf_203, \
                         kf_204 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = ab_x[k] * if__200[k]
                       + kf_200[k];

            t_301[k] = ab_x[k] * if__201[k]
                       + kf_201[k];

            t_302[k] = ab_x[k] * if__202[k]
                       + kf_202[k];

            t_303[k] = ab_x[k] * if__203[k]
                       + kf_203[k];

            t_304[k] = ab_x[k] * if__204[k]
                       + kf_204[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, if__205, if__206, if__207, \
                         if__208, if__209, kf_205, kf_206, kf_207, kf_208, \
                         kf_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = ab_x[k] * if__205[k]
                       + kf_205[k];

            t_306[k] = ab_x[k] * if__206[k]
                       + kf_206[k];

            t_307[k] = ab_x[k] * if__207[k]
                       + kf_207[k];

            t_308[k] = ab_x[k] * if__208[k]
                       + kf_208[k];

            t_309[k] = ab_x[k] * if__209[k]
                       + kf_209[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_y, ab_z, if__206, if__207, \
                         if__208, if__209, kf_266, kf_267, kf_268, kf_269, \
                         kf_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = ab_y[k] * if__206[k]
                       + kf_266[k];

            t_311[k] = ab_y[k] * if__207[k]
                       + kf_267[k];

            t_312[k] = ab_y[k] * if__208[k]
                       + kf_268[k];

            t_313[k] = ab_y[k] * if__209[k]
                       + kf_269[k];

            t_314[k] = ab_z[k] * if__209[k]
                       + kf_279[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, if__210, if__211, if__212, \
                         if__213, if__214, kf_210, kf_211, kf_212, kf_213, \
                         kf_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = ab_x[k] * if__210[k]
                       + kf_210[k];

            t_316[k] = ab_x[k] * if__211[k]
                       + kf_211[k];

            t_317[k] = ab_x[k] * if__212[k]
                       + kf_212[k];

            t_318[k] = ab_x[k] * if__213[k]
                       + kf_213[k];

            t_319[k] = ab_x[k] * if__214[k]
                       + kf_214[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, if__215, if__216, if__217, \
                         if__218, if__219, kf_215, kf_216, kf_217, kf_218, \
                         kf_219 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = ab_x[k] * if__215[k]
                       + kf_215[k];

            t_321[k] = ab_x[k] * if__216[k]
                       + kf_216[k];

            t_322[k] = ab_x[k] * if__217[k]
                       + kf_217[k];

            t_323[k] = ab_x[k] * if__218[k]
                       + kf_218[k];

            t_324[k] = ab_x[k] * if__219[k]
                       + kf_219[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_y, ab_z, if__216, if__217, \
                         if__218, if__219, kf_286, kf_287, kf_288, kf_289, \
                         kf_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = ab_y[k] * if__216[k]
                       + kf_286[k];

            t_326[k] = ab_y[k] * if__217[k]
                       + kf_287[k];

            t_327[k] = ab_y[k] * if__218[k]
                       + kf_288[k];

            t_328[k] = ab_y[k] * if__219[k]
                       + kf_289[k];

            t_329[k] = ab_z[k] * if__219[k]
                       + kf_299[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, if__220, if__221, if__222, \
                         if__223, if__224, kf_220, kf_221, kf_222, kf_223, \
                         kf_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = ab_x[k] * if__220[k]
                       + kf_220[k];

            t_331[k] = ab_x[k] * if__221[k]
                       + kf_221[k];

            t_332[k] = ab_x[k] * if__222[k]
                       + kf_222[k];

            t_333[k] = ab_x[k] * if__223[k]
                       + kf_223[k];

            t_334[k] = ab_x[k] * if__224[k]
                       + kf_224[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, if__225, if__226, if__227, \
                         if__228, if__229, kf_225, kf_226, kf_227, kf_228, \
                         kf_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = ab_x[k] * if__225[k]
                       + kf_225[k];

            t_336[k] = ab_x[k] * if__226[k]
                       + kf_226[k];

            t_337[k] = ab_x[k] * if__227[k]
                       + kf_227[k];

            t_338[k] = ab_x[k] * if__228[k]
                       + kf_228[k];

            t_339[k] = ab_x[k] * if__229[k]
                       + kf_229[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_y, ab_z, if__226, if__227, \
                         if__228, if__229, kf_296, kf_297, kf_298, kf_299, \
                         kf_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = ab_y[k] * if__226[k]
                       + kf_296[k];

            t_341[k] = ab_y[k] * if__227[k]
                       + kf_297[k];

            t_342[k] = ab_y[k] * if__228[k]
                       + kf_298[k];

            t_343[k] = ab_y[k] * if__229[k]
                       + kf_299[k];

            t_344[k] = ab_z[k] * if__229[k]
                       + kf_309[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, if__230, if__231, if__232, \
                         if__233, if__234, kf_230, kf_231, kf_232, kf_233, \
                         kf_234 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = ab_x[k] * if__230[k]
                       + kf_230[k];

            t_346[k] = ab_x[k] * if__231[k]
                       + kf_231[k];

            t_347[k] = ab_x[k] * if__232[k]
                       + kf_232[k];

            t_348[k] = ab_x[k] * if__233[k]
                       + kf_233[k];

            t_349[k] = ab_x[k] * if__234[k]
                       + kf_234[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, if__235, if__236, if__237, \
                         if__238, if__239, kf_235, kf_236, kf_237, kf_238, \
                         kf_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = ab_x[k] * if__235[k]
                       + kf_235[k];

            t_351[k] = ab_x[k] * if__236[k]
                       + kf_236[k];

            t_352[k] = ab_x[k] * if__237[k]
                       + kf_237[k];

            t_353[k] = ab_x[k] * if__238[k]
                       + kf_238[k];

            t_354[k] = ab_x[k] * if__239[k]
                       + kf_239[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_y, ab_z, if__236, if__237, \
                         if__238, if__239, kf_306, kf_307, kf_308, kf_309, \
                         kf_319 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = ab_y[k] * if__236[k]
                       + kf_306[k];

            t_356[k] = ab_y[k] * if__237[k]
                       + kf_307[k];

            t_357[k] = ab_y[k] * if__238[k]
                       + kf_308[k];

            t_358[k] = ab_y[k] * if__239[k]
                       + kf_309[k];

            t_359[k] = ab_z[k] * if__239[k]
                       + kf_319[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, if__240, if__241, if__242, \
                         if__243, if__244, kf_240, kf_241, kf_242, kf_243, \
                         kf_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * if__240[k]
                       + kf_240[k];

            t_361[k] = ab_x[k] * if__241[k]
                       + kf_241[k];

            t_362[k] = ab_x[k] * if__242[k]
                       + kf_242[k];

            t_363[k] = ab_x[k] * if__243[k]
                       + kf_243[k];

            t_364[k] = ab_x[k] * if__244[k]
                       + kf_244[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, if__245, if__246, if__247, \
                         if__248, if__249, kf_245, kf_246, kf_247, kf_248, \
                         kf_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * if__245[k]
                       + kf_245[k];

            t_366[k] = ab_x[k] * if__246[k]
                       + kf_246[k];

            t_367[k] = ab_x[k] * if__247[k]
                       + kf_247[k];

            t_368[k] = ab_x[k] * if__248[k]
                       + kf_248[k];

            t_369[k] = ab_x[k] * if__249[k]
                       + kf_249[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_y, ab_z, if__246, if__247, \
                         if__248, if__249, kf_316, kf_317, kf_318, kf_319, \
                         kf_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_y[k] * if__246[k]
                       + kf_316[k];

            t_371[k] = ab_y[k] * if__247[k]
                       + kf_317[k];

            t_372[k] = ab_y[k] * if__248[k]
                       + kf_318[k];

            t_373[k] = ab_y[k] * if__249[k]
                       + kf_319[k];

            t_374[k] = ab_z[k] * if__249[k]
                       + kf_329[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, if__250, if__251, if__252, \
                         if__253, if__254, kf_250, kf_251, kf_252, kf_253, \
                         kf_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = ab_x[k] * if__250[k]
                       + kf_250[k];

            t_376[k] = ab_x[k] * if__251[k]
                       + kf_251[k];

            t_377[k] = ab_x[k] * if__252[k]
                       + kf_252[k];

            t_378[k] = ab_x[k] * if__253[k]
                       + kf_253[k];

            t_379[k] = ab_x[k] * if__254[k]
                       + kf_254[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, if__255, if__256, if__257, \
                         if__258, if__259, kf_255, kf_256, kf_257, kf_258, \
                         kf_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = ab_x[k] * if__255[k]
                       + kf_255[k];

            t_381[k] = ab_x[k] * if__256[k]
                       + kf_256[k];

            t_382[k] = ab_x[k] * if__257[k]
                       + kf_257[k];

            t_383[k] = ab_x[k] * if__258[k]
                       + kf_258[k];

            t_384[k] = ab_x[k] * if__259[k]
                       + kf_259[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_y, ab_z, if__256, if__257, \
                         if__258, if__259, kf_326, kf_327, kf_328, kf_329, \
                         kf_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = ab_y[k] * if__256[k]
                       + kf_326[k];

            t_386[k] = ab_y[k] * if__257[k]
                       + kf_327[k];

            t_387[k] = ab_y[k] * if__258[k]
                       + kf_328[k];

            t_388[k] = ab_y[k] * if__259[k]
                       + kf_329[k];

            t_389[k] = ab_z[k] * if__259[k]
                       + kf_339[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, if__260, if__261, if__262, \
                         if__263, if__264, kf_260, kf_261, kf_262, kf_263, \
                         kf_264 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = ab_x[k] * if__260[k]
                       + kf_260[k];

            t_391[k] = ab_x[k] * if__261[k]
                       + kf_261[k];

            t_392[k] = ab_x[k] * if__262[k]
                       + kf_262[k];

            t_393[k] = ab_x[k] * if__263[k]
                       + kf_263[k];

            t_394[k] = ab_x[k] * if__264[k]
                       + kf_264[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, if__265, if__266, if__267, \
                         if__268, if__269, kf_265, kf_266, kf_267, kf_268, \
                         kf_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = ab_x[k] * if__265[k]
                       + kf_265[k];

            t_396[k] = ab_x[k] * if__266[k]
                       + kf_266[k];

            t_397[k] = ab_x[k] * if__267[k]
                       + kf_267[k];

            t_398[k] = ab_x[k] * if__268[k]
                       + kf_268[k];

            t_399[k] = ab_x[k] * if__269[k]
                       + kf_269[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_y, ab_z, if__266, if__267, \
                         if__268, if__269, kf_336, kf_337, kf_338, kf_339, \
                         kf_349 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = ab_y[k] * if__266[k]
                       + kf_336[k];

            t_401[k] = ab_y[k] * if__267[k]
                       + kf_337[k];

            t_402[k] = ab_y[k] * if__268[k]
                       + kf_338[k];

            t_403[k] = ab_y[k] * if__269[k]
                       + kf_339[k];

            t_404[k] = ab_z[k] * if__269[k]
                       + kf_349[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, if__270, if__271, if__272, \
                         if__273, if__274, kf_270, kf_271, kf_272, kf_273, \
                         kf_274 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = ab_x[k] * if__270[k]
                       + kf_270[k];

            t_406[k] = ab_x[k] * if__271[k]
                       + kf_271[k];

            t_407[k] = ab_x[k] * if__272[k]
                       + kf_272[k];

            t_408[k] = ab_x[k] * if__273[k]
                       + kf_273[k];

            t_409[k] = ab_x[k] * if__274[k]
                       + kf_274[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, if__275, if__276, if__277, \
                         if__278, if__279, kf_275, kf_276, kf_277, kf_278, \
                         kf_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = ab_x[k] * if__275[k]
                       + kf_275[k];

            t_411[k] = ab_x[k] * if__276[k]
                       + kf_276[k];

            t_412[k] = ab_x[k] * if__277[k]
                       + kf_277[k];

            t_413[k] = ab_x[k] * if__278[k]
                       + kf_278[k];

            t_414[k] = ab_x[k] * if__279[k]
                       + kf_279[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_y, ab_z, if__276, if__277, \
                         if__278, if__279, kf_346, kf_347, kf_348, kf_349, \
                         kf_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = ab_y[k] * if__276[k]
                       + kf_346[k];

            t_416[k] = ab_y[k] * if__277[k]
                       + kf_347[k];

            t_417[k] = ab_y[k] * if__278[k]
                       + kf_348[k];

            t_418[k] = ab_y[k] * if__279[k]
                       + kf_349[k];

            t_419[k] = ab_z[k] * if__279[k]
                       + kf_359[k];
        }
    }
}

auto
compute_hrr_ig_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t if_, const size_t kf,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_ig_out_of_first_piece0(buffer, coordinates, target, if_, kf, ncomps, nmax);

    compute_hrr_ig_out_of_first_piece1(buffer, coordinates, target, if_, kf, ncomps, nmax);

    compute_hrr_ig_out_of_first_piece2(buffer, coordinates, target, if_, kf, ncomps, nmax);
}

static auto
compute_hrr_ig_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t if_, const size_t kf, const size_t ncomps,
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

        const auto *if__0 = buffer.data(if_ + 0 * ncomps + c);
        const auto *if__1 = buffer.data(if_ + 1 * ncomps + c);
        const auto *if__2 = buffer.data(if_ + 2 * ncomps + c);
        const auto *if__3 = buffer.data(if_ + 3 * ncomps + c);
        const auto *if__4 = buffer.data(if_ + 4 * ncomps + c);
        const auto *if__5 = buffer.data(if_ + 5 * ncomps + c);
        const auto *if__6 = buffer.data(if_ + 6 * ncomps + c);
        const auto *if__7 = buffer.data(if_ + 7 * ncomps + c);
        const auto *if__8 = buffer.data(if_ + 8 * ncomps + c);
        const auto *if__9 = buffer.data(if_ + 9 * ncomps + c);
        const auto *if__10 = buffer.data(if_ + 10 * ncomps + c);
        const auto *if__11 = buffer.data(if_ + 11 * ncomps + c);
        const auto *if__12 = buffer.data(if_ + 12 * ncomps + c);
        const auto *if__13 = buffer.data(if_ + 13 * ncomps + c);
        const auto *if__14 = buffer.data(if_ + 14 * ncomps + c);
        const auto *if__15 = buffer.data(if_ + 15 * ncomps + c);
        const auto *if__16 = buffer.data(if_ + 16 * ncomps + c);
        const auto *if__17 = buffer.data(if_ + 17 * ncomps + c);
        const auto *if__18 = buffer.data(if_ + 18 * ncomps + c);
        const auto *if__19 = buffer.data(if_ + 19 * ncomps + c);
        const auto *if__20 = buffer.data(if_ + 20 * ncomps + c);
        const auto *if__21 = buffer.data(if_ + 21 * ncomps + c);
        const auto *if__22 = buffer.data(if_ + 22 * ncomps + c);
        const auto *if__23 = buffer.data(if_ + 23 * ncomps + c);
        const auto *if__24 = buffer.data(if_ + 24 * ncomps + c);
        const auto *if__25 = buffer.data(if_ + 25 * ncomps + c);
        const auto *if__26 = buffer.data(if_ + 26 * ncomps + c);
        const auto *if__27 = buffer.data(if_ + 27 * ncomps + c);
        const auto *if__28 = buffer.data(if_ + 28 * ncomps + c);
        const auto *if__29 = buffer.data(if_ + 29 * ncomps + c);
        const auto *if__30 = buffer.data(if_ + 30 * ncomps + c);
        const auto *if__31 = buffer.data(if_ + 31 * ncomps + c);
        const auto *if__32 = buffer.data(if_ + 32 * ncomps + c);
        const auto *if__33 = buffer.data(if_ + 33 * ncomps + c);
        const auto *if__34 = buffer.data(if_ + 34 * ncomps + c);
        const auto *if__35 = buffer.data(if_ + 35 * ncomps + c);
        const auto *if__36 = buffer.data(if_ + 36 * ncomps + c);
        const auto *if__37 = buffer.data(if_ + 37 * ncomps + c);
        const auto *if__38 = buffer.data(if_ + 38 * ncomps + c);
        const auto *if__39 = buffer.data(if_ + 39 * ncomps + c);
        const auto *if__40 = buffer.data(if_ + 40 * ncomps + c);
        const auto *if__41 = buffer.data(if_ + 41 * ncomps + c);
        const auto *if__42 = buffer.data(if_ + 42 * ncomps + c);
        const auto *if__43 = buffer.data(if_ + 43 * ncomps + c);
        const auto *if__44 = buffer.data(if_ + 44 * ncomps + c);
        const auto *if__45 = buffer.data(if_ + 45 * ncomps + c);
        const auto *if__46 = buffer.data(if_ + 46 * ncomps + c);
        const auto *if__47 = buffer.data(if_ + 47 * ncomps + c);
        const auto *if__48 = buffer.data(if_ + 48 * ncomps + c);
        const auto *if__49 = buffer.data(if_ + 49 * ncomps + c);
        const auto *if__50 = buffer.data(if_ + 50 * ncomps + c);
        const auto *if__51 = buffer.data(if_ + 51 * ncomps + c);
        const auto *if__52 = buffer.data(if_ + 52 * ncomps + c);
        const auto *if__53 = buffer.data(if_ + 53 * ncomps + c);
        const auto *if__54 = buffer.data(if_ + 54 * ncomps + c);
        const auto *if__55 = buffer.data(if_ + 55 * ncomps + c);
        const auto *if__56 = buffer.data(if_ + 56 * ncomps + c);
        const auto *if__57 = buffer.data(if_ + 57 * ncomps + c);
        const auto *if__58 = buffer.data(if_ + 58 * ncomps + c);
        const auto *if__59 = buffer.data(if_ + 59 * ncomps + c);
        const auto *if__60 = buffer.data(if_ + 60 * ncomps + c);
        const auto *if__61 = buffer.data(if_ + 61 * ncomps + c);
        const auto *if__62 = buffer.data(if_ + 62 * ncomps + c);
        const auto *if__63 = buffer.data(if_ + 63 * ncomps + c);
        const auto *if__64 = buffer.data(if_ + 64 * ncomps + c);
        const auto *if__65 = buffer.data(if_ + 65 * ncomps + c);
        const auto *if__66 = buffer.data(if_ + 66 * ncomps + c);
        const auto *if__67 = buffer.data(if_ + 67 * ncomps + c);
        const auto *if__68 = buffer.data(if_ + 68 * ncomps + c);
        const auto *if__69 = buffer.data(if_ + 69 * ncomps + c);
        const auto *if__70 = buffer.data(if_ + 70 * ncomps + c);
        const auto *if__71 = buffer.data(if_ + 71 * ncomps + c);
        const auto *if__72 = buffer.data(if_ + 72 * ncomps + c);
        const auto *if__73 = buffer.data(if_ + 73 * ncomps + c);
        const auto *if__74 = buffer.data(if_ + 74 * ncomps + c);
        const auto *if__75 = buffer.data(if_ + 75 * ncomps + c);
        const auto *if__76 = buffer.data(if_ + 76 * ncomps + c);
        const auto *if__77 = buffer.data(if_ + 77 * ncomps + c);
        const auto *if__78 = buffer.data(if_ + 78 * ncomps + c);
        const auto *if__79 = buffer.data(if_ + 79 * ncomps + c);
        const auto *if__80 = buffer.data(if_ + 80 * ncomps + c);
        const auto *if__81 = buffer.data(if_ + 81 * ncomps + c);
        const auto *if__82 = buffer.data(if_ + 82 * ncomps + c);
        const auto *if__83 = buffer.data(if_ + 83 * ncomps + c);
        const auto *if__84 = buffer.data(if_ + 84 * ncomps + c);
        const auto *if__85 = buffer.data(if_ + 85 * ncomps + c);
        const auto *if__86 = buffer.data(if_ + 86 * ncomps + c);
        const auto *if__87 = buffer.data(if_ + 87 * ncomps + c);
        const auto *if__88 = buffer.data(if_ + 88 * ncomps + c);
        const auto *if__89 = buffer.data(if_ + 89 * ncomps + c);
        const auto *if__90 = buffer.data(if_ + 90 * ncomps + c);
        const auto *if__91 = buffer.data(if_ + 91 * ncomps + c);
        const auto *if__92 = buffer.data(if_ + 92 * ncomps + c);
        const auto *if__93 = buffer.data(if_ + 93 * ncomps + c);
        const auto *if__94 = buffer.data(if_ + 94 * ncomps + c);
        const auto *if__95 = buffer.data(if_ + 95 * ncomps + c);
        const auto *if__96 = buffer.data(if_ + 96 * ncomps + c);
        const auto *if__97 = buffer.data(if_ + 97 * ncomps + c);
        const auto *if__98 = buffer.data(if_ + 98 * ncomps + c);
        const auto *if__99 = buffer.data(if_ + 99 * ncomps + c);

        const auto *kf_0 = buffer.data(kf + 0 * ncomps + c);
        const auto *kf_1 = buffer.data(kf + 1 * ncomps + c);
        const auto *kf_2 = buffer.data(kf + 2 * ncomps + c);
        const auto *kf_3 = buffer.data(kf + 3 * ncomps + c);
        const auto *kf_4 = buffer.data(kf + 4 * ncomps + c);
        const auto *kf_5 = buffer.data(kf + 5 * ncomps + c);
        const auto *kf_6 = buffer.data(kf + 6 * ncomps + c);
        const auto *kf_7 = buffer.data(kf + 7 * ncomps + c);
        const auto *kf_8 = buffer.data(kf + 8 * ncomps + c);
        const auto *kf_9 = buffer.data(kf + 9 * ncomps + c);
        const auto *kf_10 = buffer.data(kf + 10 * ncomps + c);
        const auto *kf_11 = buffer.data(kf + 11 * ncomps + c);
        const auto *kf_12 = buffer.data(kf + 12 * ncomps + c);
        const auto *kf_13 = buffer.data(kf + 13 * ncomps + c);
        const auto *kf_14 = buffer.data(kf + 14 * ncomps + c);
        const auto *kf_15 = buffer.data(kf + 15 * ncomps + c);
        const auto *kf_16 = buffer.data(kf + 16 * ncomps + c);
        const auto *kf_17 = buffer.data(kf + 17 * ncomps + c);
        const auto *kf_18 = buffer.data(kf + 18 * ncomps + c);
        const auto *kf_19 = buffer.data(kf + 19 * ncomps + c);
        const auto *kf_20 = buffer.data(kf + 20 * ncomps + c);
        const auto *kf_21 = buffer.data(kf + 21 * ncomps + c);
        const auto *kf_22 = buffer.data(kf + 22 * ncomps + c);
        const auto *kf_23 = buffer.data(kf + 23 * ncomps + c);
        const auto *kf_24 = buffer.data(kf + 24 * ncomps + c);
        const auto *kf_25 = buffer.data(kf + 25 * ncomps + c);
        const auto *kf_26 = buffer.data(kf + 26 * ncomps + c);
        const auto *kf_27 = buffer.data(kf + 27 * ncomps + c);
        const auto *kf_28 = buffer.data(kf + 28 * ncomps + c);
        const auto *kf_29 = buffer.data(kf + 29 * ncomps + c);
        const auto *kf_30 = buffer.data(kf + 30 * ncomps + c);
        const auto *kf_31 = buffer.data(kf + 31 * ncomps + c);
        const auto *kf_32 = buffer.data(kf + 32 * ncomps + c);
        const auto *kf_33 = buffer.data(kf + 33 * ncomps + c);
        const auto *kf_34 = buffer.data(kf + 34 * ncomps + c);
        const auto *kf_35 = buffer.data(kf + 35 * ncomps + c);
        const auto *kf_36 = buffer.data(kf + 36 * ncomps + c);
        const auto *kf_37 = buffer.data(kf + 37 * ncomps + c);
        const auto *kf_38 = buffer.data(kf + 38 * ncomps + c);
        const auto *kf_39 = buffer.data(kf + 39 * ncomps + c);
        const auto *kf_40 = buffer.data(kf + 40 * ncomps + c);
        const auto *kf_41 = buffer.data(kf + 41 * ncomps + c);
        const auto *kf_42 = buffer.data(kf + 42 * ncomps + c);
        const auto *kf_43 = buffer.data(kf + 43 * ncomps + c);
        const auto *kf_44 = buffer.data(kf + 44 * ncomps + c);
        const auto *kf_45 = buffer.data(kf + 45 * ncomps + c);
        const auto *kf_46 = buffer.data(kf + 46 * ncomps + c);
        const auto *kf_47 = buffer.data(kf + 47 * ncomps + c);
        const auto *kf_48 = buffer.data(kf + 48 * ncomps + c);
        const auto *kf_49 = buffer.data(kf + 49 * ncomps + c);
        const auto *kf_50 = buffer.data(kf + 50 * ncomps + c);
        const auto *kf_51 = buffer.data(kf + 51 * ncomps + c);
        const auto *kf_52 = buffer.data(kf + 52 * ncomps + c);
        const auto *kf_53 = buffer.data(kf + 53 * ncomps + c);
        const auto *kf_54 = buffer.data(kf + 54 * ncomps + c);
        const auto *kf_55 = buffer.data(kf + 55 * ncomps + c);
        const auto *kf_56 = buffer.data(kf + 56 * ncomps + c);
        const auto *kf_57 = buffer.data(kf + 57 * ncomps + c);
        const auto *kf_58 = buffer.data(kf + 58 * ncomps + c);
        const auto *kf_59 = buffer.data(kf + 59 * ncomps + c);
        const auto *kf_60 = buffer.data(kf + 60 * ncomps + c);
        const auto *kf_61 = buffer.data(kf + 61 * ncomps + c);
        const auto *kf_62 = buffer.data(kf + 62 * ncomps + c);
        const auto *kf_63 = buffer.data(kf + 63 * ncomps + c);
        const auto *kf_64 = buffer.data(kf + 64 * ncomps + c);
        const auto *kf_65 = buffer.data(kf + 65 * ncomps + c);
        const auto *kf_66 = buffer.data(kf + 66 * ncomps + c);
        const auto *kf_67 = buffer.data(kf + 67 * ncomps + c);
        const auto *kf_68 = buffer.data(kf + 68 * ncomps + c);
        const auto *kf_69 = buffer.data(kf + 69 * ncomps + c);
        const auto *kf_70 = buffer.data(kf + 70 * ncomps + c);
        const auto *kf_71 = buffer.data(kf + 71 * ncomps + c);
        const auto *kf_72 = buffer.data(kf + 72 * ncomps + c);
        const auto *kf_73 = buffer.data(kf + 73 * ncomps + c);
        const auto *kf_74 = buffer.data(kf + 74 * ncomps + c);
        const auto *kf_75 = buffer.data(kf + 75 * ncomps + c);
        const auto *kf_76 = buffer.data(kf + 76 * ncomps + c);
        const auto *kf_77 = buffer.data(kf + 77 * ncomps + c);
        const auto *kf_78 = buffer.data(kf + 78 * ncomps + c);
        const auto *kf_79 = buffer.data(kf + 79 * ncomps + c);
        const auto *kf_80 = buffer.data(kf + 80 * ncomps + c);
        const auto *kf_81 = buffer.data(kf + 81 * ncomps + c);
        const auto *kf_82 = buffer.data(kf + 82 * ncomps + c);
        const auto *kf_83 = buffer.data(kf + 83 * ncomps + c);
        const auto *kf_84 = buffer.data(kf + 84 * ncomps + c);
        const auto *kf_85 = buffer.data(kf + 85 * ncomps + c);
        const auto *kf_86 = buffer.data(kf + 86 * ncomps + c);
        const auto *kf_87 = buffer.data(kf + 87 * ncomps + c);
        const auto *kf_88 = buffer.data(kf + 88 * ncomps + c);
        const auto *kf_89 = buffer.data(kf + 89 * ncomps + c);
        const auto *kf_90 = buffer.data(kf + 90 * ncomps + c);
        const auto *kf_91 = buffer.data(kf + 91 * ncomps + c);
        const auto *kf_92 = buffer.data(kf + 92 * ncomps + c);
        const auto *kf_93 = buffer.data(kf + 93 * ncomps + c);
        const auto *kf_94 = buffer.data(kf + 94 * ncomps + c);
        const auto *kf_95 = buffer.data(kf + 95 * ncomps + c);
        const auto *kf_96 = buffer.data(kf + 96 * ncomps + c);
        const auto *kf_97 = buffer.data(kf + 97 * ncomps + c);
        const auto *kf_98 = buffer.data(kf + 98 * ncomps + c);
        const auto *kf_99 = buffer.data(kf + 99 * ncomps + c);
        const auto *kf_106 = buffer.data(kf + 106 * ncomps + c);
        const auto *kf_107 = buffer.data(kf + 107 * ncomps + c);
        const auto *kf_108 = buffer.data(kf + 108 * ncomps + c);
        const auto *kf_109 = buffer.data(kf + 109 * ncomps + c);
        const auto *kf_116 = buffer.data(kf + 116 * ncomps + c);
        const auto *kf_117 = buffer.data(kf + 117 * ncomps + c);
        const auto *kf_118 = buffer.data(kf + 118 * ncomps + c);
        const auto *kf_119 = buffer.data(kf + 119 * ncomps + c);
        const auto *kf_126 = buffer.data(kf + 126 * ncomps + c);
        const auto *kf_127 = buffer.data(kf + 127 * ncomps + c);
        const auto *kf_128 = buffer.data(kf + 128 * ncomps + c);
        const auto *kf_129 = buffer.data(kf + 129 * ncomps + c);
        const auto *kf_139 = buffer.data(kf + 139 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, if__0, if__1, if__2, if__3, if__4, \
                         kf_0, kf_1, kf_2, kf_3, kf_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * if__0[k]
                     + kf_0[k];

            t_1[k] = ab_x[k] * if__1[k]
                     + kf_1[k];

            t_2[k] = ab_x[k] * if__2[k]
                     + kf_2[k];

            t_3[k] = ab_x[k] * if__3[k]
                     + kf_3[k];

            t_4[k] = ab_x[k] * if__4[k]
                     + kf_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, if__5, if__6, if__7, if__8, if__9, \
                         kf_5, kf_6, kf_7, kf_8, kf_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * if__5[k]
                     + kf_5[k];

            t_6[k] = ab_x[k] * if__6[k]
                     + kf_6[k];

            t_7[k] = ab_x[k] * if__7[k]
                     + kf_7[k];

            t_8[k] = ab_x[k] * if__8[k]
                     + kf_8[k];

            t_9[k] = ab_x[k] * if__9[k]
                     + kf_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_y, ab_z, if__6, if__7, if__8, if__9, \
                         kf_16, kf_17, kf_18, kf_19, kf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_y[k] * if__6[k]
                      + kf_16[k];

            t_11[k] = ab_y[k] * if__7[k]
                      + kf_17[k];

            t_12[k] = ab_y[k] * if__8[k]
                      + kf_18[k];

            t_13[k] = ab_y[k] * if__9[k]
                      + kf_19[k];

            t_14[k] = ab_z[k] * if__9[k]
                      + kf_29[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, if__10, if__11, if__12, if__13, \
                         if__14, kf_10, kf_11, kf_12, kf_13, kf_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * if__10[k]
                      + kf_10[k];

            t_16[k] = ab_x[k] * if__11[k]
                      + kf_11[k];

            t_17[k] = ab_x[k] * if__12[k]
                      + kf_12[k];

            t_18[k] = ab_x[k] * if__13[k]
                      + kf_13[k];

            t_19[k] = ab_x[k] * if__14[k]
                      + kf_14[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, if__15, if__16, if__17, if__18, \
                         if__19, kf_15, kf_16, kf_17, kf_18, kf_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * if__15[k]
                      + kf_15[k];

            t_21[k] = ab_x[k] * if__16[k]
                      + kf_16[k];

            t_22[k] = ab_x[k] * if__17[k]
                      + kf_17[k];

            t_23[k] = ab_x[k] * if__18[k]
                      + kf_18[k];

            t_24[k] = ab_x[k] * if__19[k]
                      + kf_19[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_y, ab_z, if__16, if__17, if__18, \
                         if__19, kf_36, kf_37, kf_38, kf_39, kf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_y[k] * if__16[k]
                      + kf_36[k];

            t_26[k] = ab_y[k] * if__17[k]
                      + kf_37[k];

            t_27[k] = ab_y[k] * if__18[k]
                      + kf_38[k];

            t_28[k] = ab_y[k] * if__19[k]
                      + kf_39[k];

            t_29[k] = ab_z[k] * if__19[k]
                      + kf_49[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, if__20, if__21, if__22, if__23, \
                         if__24, kf_20, kf_21, kf_22, kf_23, kf_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * if__20[k]
                      + kf_20[k];

            t_31[k] = ab_x[k] * if__21[k]
                      + kf_21[k];

            t_32[k] = ab_x[k] * if__22[k]
                      + kf_22[k];

            t_33[k] = ab_x[k] * if__23[k]
                      + kf_23[k];

            t_34[k] = ab_x[k] * if__24[k]
                      + kf_24[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, if__25, if__26, if__27, if__28, \
                         if__29, kf_25, kf_26, kf_27, kf_28, kf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * if__25[k]
                      + kf_25[k];

            t_36[k] = ab_x[k] * if__26[k]
                      + kf_26[k];

            t_37[k] = ab_x[k] * if__27[k]
                      + kf_27[k];

            t_38[k] = ab_x[k] * if__28[k]
                      + kf_28[k];

            t_39[k] = ab_x[k] * if__29[k]
                      + kf_29[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_y, ab_z, if__26, if__27, if__28, \
                         if__29, kf_46, kf_47, kf_48, kf_49, kf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_y[k] * if__26[k]
                      + kf_46[k];

            t_41[k] = ab_y[k] * if__27[k]
                      + kf_47[k];

            t_42[k] = ab_y[k] * if__28[k]
                      + kf_48[k];

            t_43[k] = ab_y[k] * if__29[k]
                      + kf_49[k];

            t_44[k] = ab_z[k] * if__29[k]
                      + kf_59[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, if__30, if__31, if__32, if__33, \
                         if__34, kf_30, kf_31, kf_32, kf_33, kf_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * if__30[k]
                      + kf_30[k];

            t_46[k] = ab_x[k] * if__31[k]
                      + kf_31[k];

            t_47[k] = ab_x[k] * if__32[k]
                      + kf_32[k];

            t_48[k] = ab_x[k] * if__33[k]
                      + kf_33[k];

            t_49[k] = ab_x[k] * if__34[k]
                      + kf_34[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, if__35, if__36, if__37, if__38, \
                         if__39, kf_35, kf_36, kf_37, kf_38, kf_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * if__35[k]
                      + kf_35[k];

            t_51[k] = ab_x[k] * if__36[k]
                      + kf_36[k];

            t_52[k] = ab_x[k] * if__37[k]
                      + kf_37[k];

            t_53[k] = ab_x[k] * if__38[k]
                      + kf_38[k];

            t_54[k] = ab_x[k] * if__39[k]
                      + kf_39[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_y, ab_z, if__36, if__37, if__38, \
                         if__39, kf_66, kf_67, kf_68, kf_69, kf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_y[k] * if__36[k]
                      + kf_66[k];

            t_56[k] = ab_y[k] * if__37[k]
                      + kf_67[k];

            t_57[k] = ab_y[k] * if__38[k]
                      + kf_68[k];

            t_58[k] = ab_y[k] * if__39[k]
                      + kf_69[k];

            t_59[k] = ab_z[k] * if__39[k]
                      + kf_79[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, if__40, if__41, if__42, if__43, \
                         if__44, kf_40, kf_41, kf_42, kf_43, kf_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * if__40[k]
                      + kf_40[k];

            t_61[k] = ab_x[k] * if__41[k]
                      + kf_41[k];

            t_62[k] = ab_x[k] * if__42[k]
                      + kf_42[k];

            t_63[k] = ab_x[k] * if__43[k]
                      + kf_43[k];

            t_64[k] = ab_x[k] * if__44[k]
                      + kf_44[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, if__45, if__46, if__47, if__48, \
                         if__49, kf_45, kf_46, kf_47, kf_48, kf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * if__45[k]
                      + kf_45[k];

            t_66[k] = ab_x[k] * if__46[k]
                      + kf_46[k];

            t_67[k] = ab_x[k] * if__47[k]
                      + kf_47[k];

            t_68[k] = ab_x[k] * if__48[k]
                      + kf_48[k];

            t_69[k] = ab_x[k] * if__49[k]
                      + kf_49[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, ab_z, if__46, if__47, if__48, \
                         if__49, kf_76, kf_77, kf_78, kf_79, kf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_y[k] * if__46[k]
                      + kf_76[k];

            t_71[k] = ab_y[k] * if__47[k]
                      + kf_77[k];

            t_72[k] = ab_y[k] * if__48[k]
                      + kf_78[k];

            t_73[k] = ab_y[k] * if__49[k]
                      + kf_79[k];

            t_74[k] = ab_z[k] * if__49[k]
                      + kf_89[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, if__50, if__51, if__52, if__53, \
                         if__54, kf_50, kf_51, kf_52, kf_53, kf_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * if__50[k]
                      + kf_50[k];

            t_76[k] = ab_x[k] * if__51[k]
                      + kf_51[k];

            t_77[k] = ab_x[k] * if__52[k]
                      + kf_52[k];

            t_78[k] = ab_x[k] * if__53[k]
                      + kf_53[k];

            t_79[k] = ab_x[k] * if__54[k]
                      + kf_54[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, if__55, if__56, if__57, if__58, \
                         if__59, kf_55, kf_56, kf_57, kf_58, kf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * if__55[k]
                      + kf_55[k];

            t_81[k] = ab_x[k] * if__56[k]
                      + kf_56[k];

            t_82[k] = ab_x[k] * if__57[k]
                      + kf_57[k];

            t_83[k] = ab_x[k] * if__58[k]
                      + kf_58[k];

            t_84[k] = ab_x[k] * if__59[k]
                      + kf_59[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_y, ab_z, if__56, if__57, if__58, \
                         if__59, kf_86, kf_87, kf_88, kf_89, kf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_y[k] * if__56[k]
                      + kf_86[k];

            t_86[k] = ab_y[k] * if__57[k]
                      + kf_87[k];

            t_87[k] = ab_y[k] * if__58[k]
                      + kf_88[k];

            t_88[k] = ab_y[k] * if__59[k]
                      + kf_89[k];

            t_89[k] = ab_z[k] * if__59[k]
                      + kf_99[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, if__60, if__61, if__62, if__63, \
                         if__64, kf_60, kf_61, kf_62, kf_63, kf_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * if__60[k]
                      + kf_60[k];

            t_91[k] = ab_x[k] * if__61[k]
                      + kf_61[k];

            t_92[k] = ab_x[k] * if__62[k]
                      + kf_62[k];

            t_93[k] = ab_x[k] * if__63[k]
                      + kf_63[k];

            t_94[k] = ab_x[k] * if__64[k]
                      + kf_64[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, if__65, if__66, if__67, if__68, \
                         if__69, kf_65, kf_66, kf_67, kf_68, kf_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * if__65[k]
                      + kf_65[k];

            t_96[k] = ab_x[k] * if__66[k]
                      + kf_66[k];

            t_97[k] = ab_x[k] * if__67[k]
                      + kf_67[k];

            t_98[k] = ab_x[k] * if__68[k]
                      + kf_68[k];

            t_99[k] = ab_x[k] * if__69[k]
                      + kf_69[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_y, ab_z, if__66, if__67, \
                         if__68, if__69, kf_106, kf_107, kf_108, kf_109, \
                         kf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_y[k] * if__66[k]
                       + kf_106[k];

            t_101[k] = ab_y[k] * if__67[k]
                       + kf_107[k];

            t_102[k] = ab_y[k] * if__68[k]
                       + kf_108[k];

            t_103[k] = ab_y[k] * if__69[k]
                       + kf_109[k];

            t_104[k] = ab_z[k] * if__69[k]
                       + kf_119[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, if__70, if__71, if__72, \
                         if__73, if__74, kf_70, kf_71, kf_72, kf_73, \
                         kf_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * if__70[k]
                       + kf_70[k];

            t_106[k] = ab_x[k] * if__71[k]
                       + kf_71[k];

            t_107[k] = ab_x[k] * if__72[k]
                       + kf_72[k];

            t_108[k] = ab_x[k] * if__73[k]
                       + kf_73[k];

            t_109[k] = ab_x[k] * if__74[k]
                       + kf_74[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, if__75, if__76, if__77, \
                         if__78, if__79, kf_75, kf_76, kf_77, kf_78, \
                         kf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_x[k] * if__75[k]
                       + kf_75[k];

            t_111[k] = ab_x[k] * if__76[k]
                       + kf_76[k];

            t_112[k] = ab_x[k] * if__77[k]
                       + kf_77[k];

            t_113[k] = ab_x[k] * if__78[k]
                       + kf_78[k];

            t_114[k] = ab_x[k] * if__79[k]
                       + kf_79[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_y, ab_z, if__76, if__77, \
                         if__78, if__79, kf_116, kf_117, kf_118, kf_119, \
                         kf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_y[k] * if__76[k]
                       + kf_116[k];

            t_116[k] = ab_y[k] * if__77[k]
                       + kf_117[k];

            t_117[k] = ab_y[k] * if__78[k]
                       + kf_118[k];

            t_118[k] = ab_y[k] * if__79[k]
                       + kf_119[k];

            t_119[k] = ab_z[k] * if__79[k]
                       + kf_129[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, if__80, if__81, if__82, \
                         if__83, if__84, kf_80, kf_81, kf_82, kf_83, \
                         kf_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * if__80[k]
                       + kf_80[k];

            t_121[k] = ab_x[k] * if__81[k]
                       + kf_81[k];

            t_122[k] = ab_x[k] * if__82[k]
                       + kf_82[k];

            t_123[k] = ab_x[k] * if__83[k]
                       + kf_83[k];

            t_124[k] = ab_x[k] * if__84[k]
                       + kf_84[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, if__85, if__86, if__87, \
                         if__88, if__89, kf_85, kf_86, kf_87, kf_88, \
                         kf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * if__85[k]
                       + kf_85[k];

            t_126[k] = ab_x[k] * if__86[k]
                       + kf_86[k];

            t_127[k] = ab_x[k] * if__87[k]
                       + kf_87[k];

            t_128[k] = ab_x[k] * if__88[k]
                       + kf_88[k];

            t_129[k] = ab_x[k] * if__89[k]
                       + kf_89[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_y, ab_z, if__86, if__87, \
                         if__88, if__89, kf_126, kf_127, kf_128, kf_129, \
                         kf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_y[k] * if__86[k]
                       + kf_126[k];

            t_131[k] = ab_y[k] * if__87[k]
                       + kf_127[k];

            t_132[k] = ab_y[k] * if__88[k]
                       + kf_128[k];

            t_133[k] = ab_y[k] * if__89[k]
                       + kf_129[k];

            t_134[k] = ab_z[k] * if__89[k]
                       + kf_139[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, if__90, if__91, if__92, \
                         if__93, if__94, kf_90, kf_91, kf_92, kf_93, \
                         kf_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * if__90[k]
                       + kf_90[k];

            t_136[k] = ab_x[k] * if__91[k]
                       + kf_91[k];

            t_137[k] = ab_x[k] * if__92[k]
                       + kf_92[k];

            t_138[k] = ab_x[k] * if__93[k]
                       + kf_93[k];

            t_139[k] = ab_x[k] * if__94[k]
                       + kf_94[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, if__95, if__96, if__97, \
                         if__98, if__99, kf_95, kf_96, kf_97, kf_98, \
                         kf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * if__95[k]
                       + kf_95[k];

            t_141[k] = ab_x[k] * if__96[k]
                       + kf_96[k];

            t_142[k] = ab_x[k] * if__97[k]
                       + kf_97[k];

            t_143[k] = ab_x[k] * if__98[k]
                       + kf_98[k];

            t_144[k] = ab_x[k] * if__99[k]
                       + kf_99[k];
        }
    }
}

static auto
compute_hrr_ig_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t if_, const size_t kf, const size_t ncomps,
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
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *if__96 = buffer.data(if_ + 96 * ncomps + c);
        const auto *if__97 = buffer.data(if_ + 97 * ncomps + c);
        const auto *if__98 = buffer.data(if_ + 98 * ncomps + c);
        const auto *if__99 = buffer.data(if_ + 99 * ncomps + c);
        const auto *if__100 = buffer.data(if_ + 100 * ncomps + c);
        const auto *if__101 = buffer.data(if_ + 101 * ncomps + c);
        const auto *if__102 = buffer.data(if_ + 102 * ncomps + c);
        const auto *if__103 = buffer.data(if_ + 103 * ncomps + c);
        const auto *if__104 = buffer.data(if_ + 104 * ncomps + c);
        const auto *if__105 = buffer.data(if_ + 105 * ncomps + c);
        const auto *if__106 = buffer.data(if_ + 106 * ncomps + c);
        const auto *if__107 = buffer.data(if_ + 107 * ncomps + c);
        const auto *if__108 = buffer.data(if_ + 108 * ncomps + c);
        const auto *if__109 = buffer.data(if_ + 109 * ncomps + c);
        const auto *if__110 = buffer.data(if_ + 110 * ncomps + c);
        const auto *if__111 = buffer.data(if_ + 111 * ncomps + c);
        const auto *if__112 = buffer.data(if_ + 112 * ncomps + c);
        const auto *if__113 = buffer.data(if_ + 113 * ncomps + c);
        const auto *if__114 = buffer.data(if_ + 114 * ncomps + c);
        const auto *if__115 = buffer.data(if_ + 115 * ncomps + c);
        const auto *if__116 = buffer.data(if_ + 116 * ncomps + c);
        const auto *if__117 = buffer.data(if_ + 117 * ncomps + c);
        const auto *if__118 = buffer.data(if_ + 118 * ncomps + c);
        const auto *if__119 = buffer.data(if_ + 119 * ncomps + c);
        const auto *if__120 = buffer.data(if_ + 120 * ncomps + c);
        const auto *if__121 = buffer.data(if_ + 121 * ncomps + c);
        const auto *if__122 = buffer.data(if_ + 122 * ncomps + c);
        const auto *if__123 = buffer.data(if_ + 123 * ncomps + c);
        const auto *if__124 = buffer.data(if_ + 124 * ncomps + c);
        const auto *if__125 = buffer.data(if_ + 125 * ncomps + c);
        const auto *if__126 = buffer.data(if_ + 126 * ncomps + c);
        const auto *if__127 = buffer.data(if_ + 127 * ncomps + c);
        const auto *if__128 = buffer.data(if_ + 128 * ncomps + c);
        const auto *if__129 = buffer.data(if_ + 129 * ncomps + c);
        const auto *if__130 = buffer.data(if_ + 130 * ncomps + c);
        const auto *if__131 = buffer.data(if_ + 131 * ncomps + c);
        const auto *if__132 = buffer.data(if_ + 132 * ncomps + c);
        const auto *if__133 = buffer.data(if_ + 133 * ncomps + c);
        const auto *if__134 = buffer.data(if_ + 134 * ncomps + c);
        const auto *if__135 = buffer.data(if_ + 135 * ncomps + c);
        const auto *if__136 = buffer.data(if_ + 136 * ncomps + c);
        const auto *if__137 = buffer.data(if_ + 137 * ncomps + c);
        const auto *if__138 = buffer.data(if_ + 138 * ncomps + c);
        const auto *if__139 = buffer.data(if_ + 139 * ncomps + c);
        const auto *if__140 = buffer.data(if_ + 140 * ncomps + c);
        const auto *if__141 = buffer.data(if_ + 141 * ncomps + c);
        const auto *if__142 = buffer.data(if_ + 142 * ncomps + c);
        const auto *if__143 = buffer.data(if_ + 143 * ncomps + c);
        const auto *if__144 = buffer.data(if_ + 144 * ncomps + c);
        const auto *if__145 = buffer.data(if_ + 145 * ncomps + c);
        const auto *if__146 = buffer.data(if_ + 146 * ncomps + c);
        const auto *if__147 = buffer.data(if_ + 147 * ncomps + c);
        const auto *if__148 = buffer.data(if_ + 148 * ncomps + c);
        const auto *if__149 = buffer.data(if_ + 149 * ncomps + c);
        const auto *if__150 = buffer.data(if_ + 150 * ncomps + c);
        const auto *if__151 = buffer.data(if_ + 151 * ncomps + c);
        const auto *if__152 = buffer.data(if_ + 152 * ncomps + c);
        const auto *if__153 = buffer.data(if_ + 153 * ncomps + c);
        const auto *if__154 = buffer.data(if_ + 154 * ncomps + c);
        const auto *if__155 = buffer.data(if_ + 155 * ncomps + c);
        const auto *if__156 = buffer.data(if_ + 156 * ncomps + c);
        const auto *if__157 = buffer.data(if_ + 157 * ncomps + c);
        const auto *if__158 = buffer.data(if_ + 158 * ncomps + c);
        const auto *if__159 = buffer.data(if_ + 159 * ncomps + c);
        const auto *if__160 = buffer.data(if_ + 160 * ncomps + c);
        const auto *if__161 = buffer.data(if_ + 161 * ncomps + c);
        const auto *if__162 = buffer.data(if_ + 162 * ncomps + c);
        const auto *if__163 = buffer.data(if_ + 163 * ncomps + c);
        const auto *if__164 = buffer.data(if_ + 164 * ncomps + c);
        const auto *if__165 = buffer.data(if_ + 165 * ncomps + c);
        const auto *if__166 = buffer.data(if_ + 166 * ncomps + c);
        const auto *if__167 = buffer.data(if_ + 167 * ncomps + c);
        const auto *if__168 = buffer.data(if_ + 168 * ncomps + c);
        const auto *if__169 = buffer.data(if_ + 169 * ncomps + c);
        const auto *if__170 = buffer.data(if_ + 170 * ncomps + c);
        const auto *if__171 = buffer.data(if_ + 171 * ncomps + c);
        const auto *if__172 = buffer.data(if_ + 172 * ncomps + c);
        const auto *if__173 = buffer.data(if_ + 173 * ncomps + c);
        const auto *if__174 = buffer.data(if_ + 174 * ncomps + c);
        const auto *if__175 = buffer.data(if_ + 175 * ncomps + c);
        const auto *if__176 = buffer.data(if_ + 176 * ncomps + c);
        const auto *if__177 = buffer.data(if_ + 177 * ncomps + c);
        const auto *if__178 = buffer.data(if_ + 178 * ncomps + c);
        const auto *if__179 = buffer.data(if_ + 179 * ncomps + c);
        const auto *if__180 = buffer.data(if_ + 180 * ncomps + c);
        const auto *if__181 = buffer.data(if_ + 181 * ncomps + c);
        const auto *if__182 = buffer.data(if_ + 182 * ncomps + c);
        const auto *if__183 = buffer.data(if_ + 183 * ncomps + c);
        const auto *if__184 = buffer.data(if_ + 184 * ncomps + c);
        const auto *if__185 = buffer.data(if_ + 185 * ncomps + c);
        const auto *if__186 = buffer.data(if_ + 186 * ncomps + c);
        const auto *if__187 = buffer.data(if_ + 187 * ncomps + c);
        const auto *if__188 = buffer.data(if_ + 188 * ncomps + c);
        const auto *if__189 = buffer.data(if_ + 189 * ncomps + c);
        const auto *if__190 = buffer.data(if_ + 190 * ncomps + c);
        const auto *if__191 = buffer.data(if_ + 191 * ncomps + c);
        const auto *if__192 = buffer.data(if_ + 192 * ncomps + c);
        const auto *if__193 = buffer.data(if_ + 193 * ncomps + c);
        const auto *if__194 = buffer.data(if_ + 194 * ncomps + c);

        const auto *kf_100 = buffer.data(kf + 100 * ncomps + c);
        const auto *kf_101 = buffer.data(kf + 101 * ncomps + c);
        const auto *kf_102 = buffer.data(kf + 102 * ncomps + c);
        const auto *kf_103 = buffer.data(kf + 103 * ncomps + c);
        const auto *kf_104 = buffer.data(kf + 104 * ncomps + c);
        const auto *kf_105 = buffer.data(kf + 105 * ncomps + c);
        const auto *kf_106 = buffer.data(kf + 106 * ncomps + c);
        const auto *kf_107 = buffer.data(kf + 107 * ncomps + c);
        const auto *kf_108 = buffer.data(kf + 108 * ncomps + c);
        const auto *kf_109 = buffer.data(kf + 109 * ncomps + c);
        const auto *kf_110 = buffer.data(kf + 110 * ncomps + c);
        const auto *kf_111 = buffer.data(kf + 111 * ncomps + c);
        const auto *kf_112 = buffer.data(kf + 112 * ncomps + c);
        const auto *kf_113 = buffer.data(kf + 113 * ncomps + c);
        const auto *kf_114 = buffer.data(kf + 114 * ncomps + c);
        const auto *kf_115 = buffer.data(kf + 115 * ncomps + c);
        const auto *kf_116 = buffer.data(kf + 116 * ncomps + c);
        const auto *kf_117 = buffer.data(kf + 117 * ncomps + c);
        const auto *kf_118 = buffer.data(kf + 118 * ncomps + c);
        const auto *kf_119 = buffer.data(kf + 119 * ncomps + c);
        const auto *kf_120 = buffer.data(kf + 120 * ncomps + c);
        const auto *kf_121 = buffer.data(kf + 121 * ncomps + c);
        const auto *kf_122 = buffer.data(kf + 122 * ncomps + c);
        const auto *kf_123 = buffer.data(kf + 123 * ncomps + c);
        const auto *kf_124 = buffer.data(kf + 124 * ncomps + c);
        const auto *kf_125 = buffer.data(kf + 125 * ncomps + c);
        const auto *kf_126 = buffer.data(kf + 126 * ncomps + c);
        const auto *kf_127 = buffer.data(kf + 127 * ncomps + c);
        const auto *kf_128 = buffer.data(kf + 128 * ncomps + c);
        const auto *kf_129 = buffer.data(kf + 129 * ncomps + c);
        const auto *kf_130 = buffer.data(kf + 130 * ncomps + c);
        const auto *kf_131 = buffer.data(kf + 131 * ncomps + c);
        const auto *kf_132 = buffer.data(kf + 132 * ncomps + c);
        const auto *kf_133 = buffer.data(kf + 133 * ncomps + c);
        const auto *kf_134 = buffer.data(kf + 134 * ncomps + c);
        const auto *kf_135 = buffer.data(kf + 135 * ncomps + c);
        const auto *kf_136 = buffer.data(kf + 136 * ncomps + c);
        const auto *kf_137 = buffer.data(kf + 137 * ncomps + c);
        const auto *kf_138 = buffer.data(kf + 138 * ncomps + c);
        const auto *kf_139 = buffer.data(kf + 139 * ncomps + c);
        const auto *kf_140 = buffer.data(kf + 140 * ncomps + c);
        const auto *kf_141 = buffer.data(kf + 141 * ncomps + c);
        const auto *kf_142 = buffer.data(kf + 142 * ncomps + c);
        const auto *kf_143 = buffer.data(kf + 143 * ncomps + c);
        const auto *kf_144 = buffer.data(kf + 144 * ncomps + c);
        const auto *kf_145 = buffer.data(kf + 145 * ncomps + c);
        const auto *kf_146 = buffer.data(kf + 146 * ncomps + c);
        const auto *kf_147 = buffer.data(kf + 147 * ncomps + c);
        const auto *kf_148 = buffer.data(kf + 148 * ncomps + c);
        const auto *kf_149 = buffer.data(kf + 149 * ncomps + c);
        const auto *kf_150 = buffer.data(kf + 150 * ncomps + c);
        const auto *kf_151 = buffer.data(kf + 151 * ncomps + c);
        const auto *kf_152 = buffer.data(kf + 152 * ncomps + c);
        const auto *kf_153 = buffer.data(kf + 153 * ncomps + c);
        const auto *kf_154 = buffer.data(kf + 154 * ncomps + c);
        const auto *kf_155 = buffer.data(kf + 155 * ncomps + c);
        const auto *kf_156 = buffer.data(kf + 156 * ncomps + c);
        const auto *kf_157 = buffer.data(kf + 157 * ncomps + c);
        const auto *kf_158 = buffer.data(kf + 158 * ncomps + c);
        const auto *kf_159 = buffer.data(kf + 159 * ncomps + c);
        const auto *kf_160 = buffer.data(kf + 160 * ncomps + c);
        const auto *kf_161 = buffer.data(kf + 161 * ncomps + c);
        const auto *kf_162 = buffer.data(kf + 162 * ncomps + c);
        const auto *kf_163 = buffer.data(kf + 163 * ncomps + c);
        const auto *kf_164 = buffer.data(kf + 164 * ncomps + c);
        const auto *kf_165 = buffer.data(kf + 165 * ncomps + c);
        const auto *kf_166 = buffer.data(kf + 166 * ncomps + c);
        const auto *kf_167 = buffer.data(kf + 167 * ncomps + c);
        const auto *kf_168 = buffer.data(kf + 168 * ncomps + c);
        const auto *kf_169 = buffer.data(kf + 169 * ncomps + c);
        const auto *kf_170 = buffer.data(kf + 170 * ncomps + c);
        const auto *kf_171 = buffer.data(kf + 171 * ncomps + c);
        const auto *kf_172 = buffer.data(kf + 172 * ncomps + c);
        const auto *kf_173 = buffer.data(kf + 173 * ncomps + c);
        const auto *kf_174 = buffer.data(kf + 174 * ncomps + c);
        const auto *kf_175 = buffer.data(kf + 175 * ncomps + c);
        const auto *kf_176 = buffer.data(kf + 176 * ncomps + c);
        const auto *kf_177 = buffer.data(kf + 177 * ncomps + c);
        const auto *kf_178 = buffer.data(kf + 178 * ncomps + c);
        const auto *kf_179 = buffer.data(kf + 179 * ncomps + c);
        const auto *kf_180 = buffer.data(kf + 180 * ncomps + c);
        const auto *kf_181 = buffer.data(kf + 181 * ncomps + c);
        const auto *kf_182 = buffer.data(kf + 182 * ncomps + c);
        const auto *kf_183 = buffer.data(kf + 183 * ncomps + c);
        const auto *kf_184 = buffer.data(kf + 184 * ncomps + c);
        const auto *kf_185 = buffer.data(kf + 185 * ncomps + c);
        const auto *kf_186 = buffer.data(kf + 186 * ncomps + c);
        const auto *kf_187 = buffer.data(kf + 187 * ncomps + c);
        const auto *kf_188 = buffer.data(kf + 188 * ncomps + c);
        const auto *kf_189 = buffer.data(kf + 189 * ncomps + c);
        const auto *kf_190 = buffer.data(kf + 190 * ncomps + c);
        const auto *kf_191 = buffer.data(kf + 191 * ncomps + c);
        const auto *kf_192 = buffer.data(kf + 192 * ncomps + c);
        const auto *kf_193 = buffer.data(kf + 193 * ncomps + c);
        const auto *kf_194 = buffer.data(kf + 194 * ncomps + c);
        const auto *kf_196 = buffer.data(kf + 196 * ncomps + c);
        const auto *kf_197 = buffer.data(kf + 197 * ncomps + c);
        const auto *kf_198 = buffer.data(kf + 198 * ncomps + c);
        const auto *kf_199 = buffer.data(kf + 199 * ncomps + c);
        const auto *kf_209 = buffer.data(kf + 209 * ncomps + c);
        const auto *kf_216 = buffer.data(kf + 216 * ncomps + c);
        const auto *kf_217 = buffer.data(kf + 217 * ncomps + c);
        const auto *kf_218 = buffer.data(kf + 218 * ncomps + c);
        const auto *kf_219 = buffer.data(kf + 219 * ncomps + c);
        const auto *kf_226 = buffer.data(kf + 226 * ncomps + c);
        const auto *kf_227 = buffer.data(kf + 227 * ncomps + c);
        const auto *kf_228 = buffer.data(kf + 228 * ncomps + c);
        const auto *kf_229 = buffer.data(kf + 229 * ncomps + c);
        const auto *kf_236 = buffer.data(kf + 236 * ncomps + c);
        const auto *kf_237 = buffer.data(kf + 237 * ncomps + c);
        const auto *kf_238 = buffer.data(kf + 238 * ncomps + c);
        const auto *kf_239 = buffer.data(kf + 239 * ncomps + c);
        const auto *kf_246 = buffer.data(kf + 246 * ncomps + c);
        const auto *kf_247 = buffer.data(kf + 247 * ncomps + c);
        const auto *kf_248 = buffer.data(kf + 248 * ncomps + c);
        const auto *kf_249 = buffer.data(kf + 249 * ncomps + c);
        const auto *kf_259 = buffer.data(kf + 259 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_y, ab_z, if__96, if__97, \
                         if__98, if__99, kf_136, kf_137, kf_138, kf_139, \
                         kf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_y[k] * if__96[k]
                       + kf_136[k];

            t_146[k] = ab_y[k] * if__97[k]
                       + kf_137[k];

            t_147[k] = ab_y[k] * if__98[k]
                       + kf_138[k];

            t_148[k] = ab_y[k] * if__99[k]
                       + kf_139[k];

            t_149[k] = ab_z[k] * if__99[k]
                       + kf_149[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, if__100, if__101, if__102, \
                         if__103, if__104, kf_100, kf_101, kf_102, kf_103, \
                         kf_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * if__100[k]
                       + kf_100[k];

            t_151[k] = ab_x[k] * if__101[k]
                       + kf_101[k];

            t_152[k] = ab_x[k] * if__102[k]
                       + kf_102[k];

            t_153[k] = ab_x[k] * if__103[k]
                       + kf_103[k];

            t_154[k] = ab_x[k] * if__104[k]
                       + kf_104[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, if__105, if__106, if__107, \
                         if__108, if__109, kf_105, kf_106, kf_107, kf_108, \
                         kf_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * if__105[k]
                       + kf_105[k];

            t_156[k] = ab_x[k] * if__106[k]
                       + kf_106[k];

            t_157[k] = ab_x[k] * if__107[k]
                       + kf_107[k];

            t_158[k] = ab_x[k] * if__108[k]
                       + kf_108[k];

            t_159[k] = ab_x[k] * if__109[k]
                       + kf_109[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_y, ab_z, if__106, if__107, \
                         if__108, if__109, kf_156, kf_157, kf_158, kf_159, \
                         kf_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_y[k] * if__106[k]
                       + kf_156[k];

            t_161[k] = ab_y[k] * if__107[k]
                       + kf_157[k];

            t_162[k] = ab_y[k] * if__108[k]
                       + kf_158[k];

            t_163[k] = ab_y[k] * if__109[k]
                       + kf_159[k];

            t_164[k] = ab_z[k] * if__109[k]
                       + kf_169[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, if__110, if__111, if__112, \
                         if__113, if__114, kf_110, kf_111, kf_112, kf_113, \
                         kf_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * if__110[k]
                       + kf_110[k];

            t_166[k] = ab_x[k] * if__111[k]
                       + kf_111[k];

            t_167[k] = ab_x[k] * if__112[k]
                       + kf_112[k];

            t_168[k] = ab_x[k] * if__113[k]
                       + kf_113[k];

            t_169[k] = ab_x[k] * if__114[k]
                       + kf_114[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, if__115, if__116, if__117, \
                         if__118, if__119, kf_115, kf_116, kf_117, kf_118, \
                         kf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = ab_x[k] * if__115[k]
                       + kf_115[k];

            t_171[k] = ab_x[k] * if__116[k]
                       + kf_116[k];

            t_172[k] = ab_x[k] * if__117[k]
                       + kf_117[k];

            t_173[k] = ab_x[k] * if__118[k]
                       + kf_118[k];

            t_174[k] = ab_x[k] * if__119[k]
                       + kf_119[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_y, ab_z, if__116, if__117, \
                         if__118, if__119, kf_166, kf_167, kf_168, kf_169, \
                         kf_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_y[k] * if__116[k]
                       + kf_166[k];

            t_176[k] = ab_y[k] * if__117[k]
                       + kf_167[k];

            t_177[k] = ab_y[k] * if__118[k]
                       + kf_168[k];

            t_178[k] = ab_y[k] * if__119[k]
                       + kf_169[k];

            t_179[k] = ab_z[k] * if__119[k]
                       + kf_179[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, if__120, if__121, if__122, \
                         if__123, if__124, kf_120, kf_121, kf_122, kf_123, \
                         kf_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * if__120[k]
                       + kf_120[k];

            t_181[k] = ab_x[k] * if__121[k]
                       + kf_121[k];

            t_182[k] = ab_x[k] * if__122[k]
                       + kf_122[k];

            t_183[k] = ab_x[k] * if__123[k]
                       + kf_123[k];

            t_184[k] = ab_x[k] * if__124[k]
                       + kf_124[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, if__125, if__126, if__127, \
                         if__128, if__129, kf_125, kf_126, kf_127, kf_128, \
                         kf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * if__125[k]
                       + kf_125[k];

            t_186[k] = ab_x[k] * if__126[k]
                       + kf_126[k];

            t_187[k] = ab_x[k] * if__127[k]
                       + kf_127[k];

            t_188[k] = ab_x[k] * if__128[k]
                       + kf_128[k];

            t_189[k] = ab_x[k] * if__129[k]
                       + kf_129[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_y, ab_z, if__126, if__127, \
                         if__128, if__129, kf_176, kf_177, kf_178, kf_179, \
                         kf_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_y[k] * if__126[k]
                       + kf_176[k];

            t_191[k] = ab_y[k] * if__127[k]
                       + kf_177[k];

            t_192[k] = ab_y[k] * if__128[k]
                       + kf_178[k];

            t_193[k] = ab_y[k] * if__129[k]
                       + kf_179[k];

            t_194[k] = ab_z[k] * if__129[k]
                       + kf_189[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, if__130, if__131, if__132, \
                         if__133, if__134, kf_130, kf_131, kf_132, kf_133, \
                         kf_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * if__130[k]
                       + kf_130[k];

            t_196[k] = ab_x[k] * if__131[k]
                       + kf_131[k];

            t_197[k] = ab_x[k] * if__132[k]
                       + kf_132[k];

            t_198[k] = ab_x[k] * if__133[k]
                       + kf_133[k];

            t_199[k] = ab_x[k] * if__134[k]
                       + kf_134[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, if__135, if__136, if__137, \
                         if__138, if__139, kf_135, kf_136, kf_137, kf_138, \
                         kf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = ab_x[k] * if__135[k]
                       + kf_135[k];

            t_201[k] = ab_x[k] * if__136[k]
                       + kf_136[k];

            t_202[k] = ab_x[k] * if__137[k]
                       + kf_137[k];

            t_203[k] = ab_x[k] * if__138[k]
                       + kf_138[k];

            t_204[k] = ab_x[k] * if__139[k]
                       + kf_139[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_y, ab_z, if__136, if__137, \
                         if__138, if__139, kf_186, kf_187, kf_188, kf_189, \
                         kf_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_y[k] * if__136[k]
                       + kf_186[k];

            t_206[k] = ab_y[k] * if__137[k]
                       + kf_187[k];

            t_207[k] = ab_y[k] * if__138[k]
                       + kf_188[k];

            t_208[k] = ab_y[k] * if__139[k]
                       + kf_189[k];

            t_209[k] = ab_z[k] * if__139[k]
                       + kf_199[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, if__140, if__141, if__142, \
                         if__143, if__144, kf_140, kf_141, kf_142, kf_143, \
                         kf_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * if__140[k]
                       + kf_140[k];

            t_211[k] = ab_x[k] * if__141[k]
                       + kf_141[k];

            t_212[k] = ab_x[k] * if__142[k]
                       + kf_142[k];

            t_213[k] = ab_x[k] * if__143[k]
                       + kf_143[k];

            t_214[k] = ab_x[k] * if__144[k]
                       + kf_144[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, if__145, if__146, if__147, \
                         if__148, if__149, kf_145, kf_146, kf_147, kf_148, \
                         kf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * if__145[k]
                       + kf_145[k];

            t_216[k] = ab_x[k] * if__146[k]
                       + kf_146[k];

            t_217[k] = ab_x[k] * if__147[k]
                       + kf_147[k];

            t_218[k] = ab_x[k] * if__148[k]
                       + kf_148[k];

            t_219[k] = ab_x[k] * if__149[k]
                       + kf_149[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_y, ab_z, if__146, if__147, \
                         if__148, if__149, kf_196, kf_197, kf_198, kf_199, \
                         kf_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_y[k] * if__146[k]
                       + kf_196[k];

            t_221[k] = ab_y[k] * if__147[k]
                       + kf_197[k];

            t_222[k] = ab_y[k] * if__148[k]
                       + kf_198[k];

            t_223[k] = ab_y[k] * if__149[k]
                       + kf_199[k];

            t_224[k] = ab_z[k] * if__149[k]
                       + kf_209[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, if__150, if__151, if__152, \
                         if__153, if__154, kf_150, kf_151, kf_152, kf_153, \
                         kf_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_x[k] * if__150[k]
                       + kf_150[k];

            t_226[k] = ab_x[k] * if__151[k]
                       + kf_151[k];

            t_227[k] = ab_x[k] * if__152[k]
                       + kf_152[k];

            t_228[k] = ab_x[k] * if__153[k]
                       + kf_153[k];

            t_229[k] = ab_x[k] * if__154[k]
                       + kf_154[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, if__155, if__156, if__157, \
                         if__158, if__159, kf_155, kf_156, kf_157, kf_158, \
                         kf_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_x[k] * if__155[k]
                       + kf_155[k];

            t_231[k] = ab_x[k] * if__156[k]
                       + kf_156[k];

            t_232[k] = ab_x[k] * if__157[k]
                       + kf_157[k];

            t_233[k] = ab_x[k] * if__158[k]
                       + kf_158[k];

            t_234[k] = ab_x[k] * if__159[k]
                       + kf_159[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_y, ab_z, if__156, if__157, \
                         if__158, if__159, kf_216, kf_217, kf_218, kf_219, \
                         kf_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = ab_y[k] * if__156[k]
                       + kf_216[k];

            t_236[k] = ab_y[k] * if__157[k]
                       + kf_217[k];

            t_237[k] = ab_y[k] * if__158[k]
                       + kf_218[k];

            t_238[k] = ab_y[k] * if__159[k]
                       + kf_219[k];

            t_239[k] = ab_z[k] * if__159[k]
                       + kf_229[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, if__160, if__161, if__162, \
                         if__163, if__164, kf_160, kf_161, kf_162, kf_163, \
                         kf_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = ab_x[k] * if__160[k]
                       + kf_160[k];

            t_241[k] = ab_x[k] * if__161[k]
                       + kf_161[k];

            t_242[k] = ab_x[k] * if__162[k]
                       + kf_162[k];

            t_243[k] = ab_x[k] * if__163[k]
                       + kf_163[k];

            t_244[k] = ab_x[k] * if__164[k]
                       + kf_164[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, if__165, if__166, if__167, \
                         if__168, if__169, kf_165, kf_166, kf_167, kf_168, \
                         kf_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = ab_x[k] * if__165[k]
                       + kf_165[k];

            t_246[k] = ab_x[k] * if__166[k]
                       + kf_166[k];

            t_247[k] = ab_x[k] * if__167[k]
                       + kf_167[k];

            t_248[k] = ab_x[k] * if__168[k]
                       + kf_168[k];

            t_249[k] = ab_x[k] * if__169[k]
                       + kf_169[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_y, ab_z, if__166, if__167, \
                         if__168, if__169, kf_226, kf_227, kf_228, kf_229, \
                         kf_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = ab_y[k] * if__166[k]
                       + kf_226[k];

            t_251[k] = ab_y[k] * if__167[k]
                       + kf_227[k];

            t_252[k] = ab_y[k] * if__168[k]
                       + kf_228[k];

            t_253[k] = ab_y[k] * if__169[k]
                       + kf_229[k];

            t_254[k] = ab_z[k] * if__169[k]
                       + kf_239[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, if__170, if__171, if__172, \
                         if__173, if__174, kf_170, kf_171, kf_172, kf_173, \
                         kf_174 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = ab_x[k] * if__170[k]
                       + kf_170[k];

            t_256[k] = ab_x[k] * if__171[k]
                       + kf_171[k];

            t_257[k] = ab_x[k] * if__172[k]
                       + kf_172[k];

            t_258[k] = ab_x[k] * if__173[k]
                       + kf_173[k];

            t_259[k] = ab_x[k] * if__174[k]
                       + kf_174[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, if__175, if__176, if__177, \
                         if__178, if__179, kf_175, kf_176, kf_177, kf_178, \
                         kf_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = ab_x[k] * if__175[k]
                       + kf_175[k];

            t_261[k] = ab_x[k] * if__176[k]
                       + kf_176[k];

            t_262[k] = ab_x[k] * if__177[k]
                       + kf_177[k];

            t_263[k] = ab_x[k] * if__178[k]
                       + kf_178[k];

            t_264[k] = ab_x[k] * if__179[k]
                       + kf_179[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_y, ab_z, if__176, if__177, \
                         if__178, if__179, kf_236, kf_237, kf_238, kf_239, \
                         kf_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = ab_y[k] * if__176[k]
                       + kf_236[k];

            t_266[k] = ab_y[k] * if__177[k]
                       + kf_237[k];

            t_267[k] = ab_y[k] * if__178[k]
                       + kf_238[k];

            t_268[k] = ab_y[k] * if__179[k]
                       + kf_239[k];

            t_269[k] = ab_z[k] * if__179[k]
                       + kf_249[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, if__180, if__181, if__182, \
                         if__183, if__184, kf_180, kf_181, kf_182, kf_183, \
                         kf_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = ab_x[k] * if__180[k]
                       + kf_180[k];

            t_271[k] = ab_x[k] * if__181[k]
                       + kf_181[k];

            t_272[k] = ab_x[k] * if__182[k]
                       + kf_182[k];

            t_273[k] = ab_x[k] * if__183[k]
                       + kf_183[k];

            t_274[k] = ab_x[k] * if__184[k]
                       + kf_184[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, if__185, if__186, if__187, \
                         if__188, if__189, kf_185, kf_186, kf_187, kf_188, \
                         kf_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = ab_x[k] * if__185[k]
                       + kf_185[k];

            t_276[k] = ab_x[k] * if__186[k]
                       + kf_186[k];

            t_277[k] = ab_x[k] * if__187[k]
                       + kf_187[k];

            t_278[k] = ab_x[k] * if__188[k]
                       + kf_188[k];

            t_279[k] = ab_x[k] * if__189[k]
                       + kf_189[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_y, ab_z, if__186, if__187, \
                         if__188, if__189, kf_246, kf_247, kf_248, kf_249, \
                         kf_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = ab_y[k] * if__186[k]
                       + kf_246[k];

            t_281[k] = ab_y[k] * if__187[k]
                       + kf_247[k];

            t_282[k] = ab_y[k] * if__188[k]
                       + kf_248[k];

            t_283[k] = ab_y[k] * if__189[k]
                       + kf_249[k];

            t_284[k] = ab_z[k] * if__189[k]
                       + kf_259[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, if__190, if__191, if__192, \
                         if__193, if__194, kf_190, kf_191, kf_192, kf_193, \
                         kf_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = ab_x[k] * if__190[k]
                       + kf_190[k];

            t_286[k] = ab_x[k] * if__191[k]
                       + kf_191[k];

            t_287[k] = ab_x[k] * if__192[k]
                       + kf_192[k];

            t_288[k] = ab_x[k] * if__193[k]
                       + kf_193[k];

            t_289[k] = ab_x[k] * if__194[k]
                       + kf_194[k];
        }
    }
}

static auto
compute_hrr_ig_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t if_, const size_t kf, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *if__195 = buffer.data(if_ + 195 * ncomps + c);
        const auto *if__196 = buffer.data(if_ + 196 * ncomps + c);
        const auto *if__197 = buffer.data(if_ + 197 * ncomps + c);
        const auto *if__198 = buffer.data(if_ + 198 * ncomps + c);
        const auto *if__199 = buffer.data(if_ + 199 * ncomps + c);
        const auto *if__200 = buffer.data(if_ + 200 * ncomps + c);
        const auto *if__201 = buffer.data(if_ + 201 * ncomps + c);
        const auto *if__202 = buffer.data(if_ + 202 * ncomps + c);
        const auto *if__203 = buffer.data(if_ + 203 * ncomps + c);
        const auto *if__204 = buffer.data(if_ + 204 * ncomps + c);
        const auto *if__205 = buffer.data(if_ + 205 * ncomps + c);
        const auto *if__206 = buffer.data(if_ + 206 * ncomps + c);
        const auto *if__207 = buffer.data(if_ + 207 * ncomps + c);
        const auto *if__208 = buffer.data(if_ + 208 * ncomps + c);
        const auto *if__209 = buffer.data(if_ + 209 * ncomps + c);
        const auto *if__210 = buffer.data(if_ + 210 * ncomps + c);
        const auto *if__211 = buffer.data(if_ + 211 * ncomps + c);
        const auto *if__212 = buffer.data(if_ + 212 * ncomps + c);
        const auto *if__213 = buffer.data(if_ + 213 * ncomps + c);
        const auto *if__214 = buffer.data(if_ + 214 * ncomps + c);
        const auto *if__215 = buffer.data(if_ + 215 * ncomps + c);
        const auto *if__216 = buffer.data(if_ + 216 * ncomps + c);
        const auto *if__217 = buffer.data(if_ + 217 * ncomps + c);
        const auto *if__218 = buffer.data(if_ + 218 * ncomps + c);
        const auto *if__219 = buffer.data(if_ + 219 * ncomps + c);
        const auto *if__220 = buffer.data(if_ + 220 * ncomps + c);
        const auto *if__221 = buffer.data(if_ + 221 * ncomps + c);
        const auto *if__222 = buffer.data(if_ + 222 * ncomps + c);
        const auto *if__223 = buffer.data(if_ + 223 * ncomps + c);
        const auto *if__224 = buffer.data(if_ + 224 * ncomps + c);
        const auto *if__225 = buffer.data(if_ + 225 * ncomps + c);
        const auto *if__226 = buffer.data(if_ + 226 * ncomps + c);
        const auto *if__227 = buffer.data(if_ + 227 * ncomps + c);
        const auto *if__228 = buffer.data(if_ + 228 * ncomps + c);
        const auto *if__229 = buffer.data(if_ + 229 * ncomps + c);
        const auto *if__230 = buffer.data(if_ + 230 * ncomps + c);
        const auto *if__231 = buffer.data(if_ + 231 * ncomps + c);
        const auto *if__232 = buffer.data(if_ + 232 * ncomps + c);
        const auto *if__233 = buffer.data(if_ + 233 * ncomps + c);
        const auto *if__234 = buffer.data(if_ + 234 * ncomps + c);
        const auto *if__235 = buffer.data(if_ + 235 * ncomps + c);
        const auto *if__236 = buffer.data(if_ + 236 * ncomps + c);
        const auto *if__237 = buffer.data(if_ + 237 * ncomps + c);
        const auto *if__238 = buffer.data(if_ + 238 * ncomps + c);
        const auto *if__239 = buffer.data(if_ + 239 * ncomps + c);
        const auto *if__240 = buffer.data(if_ + 240 * ncomps + c);
        const auto *if__241 = buffer.data(if_ + 241 * ncomps + c);
        const auto *if__242 = buffer.data(if_ + 242 * ncomps + c);
        const auto *if__243 = buffer.data(if_ + 243 * ncomps + c);
        const auto *if__244 = buffer.data(if_ + 244 * ncomps + c);
        const auto *if__245 = buffer.data(if_ + 245 * ncomps + c);
        const auto *if__246 = buffer.data(if_ + 246 * ncomps + c);
        const auto *if__247 = buffer.data(if_ + 247 * ncomps + c);
        const auto *if__248 = buffer.data(if_ + 248 * ncomps + c);
        const auto *if__249 = buffer.data(if_ + 249 * ncomps + c);
        const auto *if__250 = buffer.data(if_ + 250 * ncomps + c);
        const auto *if__251 = buffer.data(if_ + 251 * ncomps + c);
        const auto *if__252 = buffer.data(if_ + 252 * ncomps + c);
        const auto *if__253 = buffer.data(if_ + 253 * ncomps + c);
        const auto *if__254 = buffer.data(if_ + 254 * ncomps + c);
        const auto *if__255 = buffer.data(if_ + 255 * ncomps + c);
        const auto *if__256 = buffer.data(if_ + 256 * ncomps + c);
        const auto *if__257 = buffer.data(if_ + 257 * ncomps + c);
        const auto *if__258 = buffer.data(if_ + 258 * ncomps + c);
        const auto *if__259 = buffer.data(if_ + 259 * ncomps + c);
        const auto *if__260 = buffer.data(if_ + 260 * ncomps + c);
        const auto *if__261 = buffer.data(if_ + 261 * ncomps + c);
        const auto *if__262 = buffer.data(if_ + 262 * ncomps + c);
        const auto *if__263 = buffer.data(if_ + 263 * ncomps + c);
        const auto *if__264 = buffer.data(if_ + 264 * ncomps + c);
        const auto *if__265 = buffer.data(if_ + 265 * ncomps + c);
        const auto *if__266 = buffer.data(if_ + 266 * ncomps + c);
        const auto *if__267 = buffer.data(if_ + 267 * ncomps + c);
        const auto *if__268 = buffer.data(if_ + 268 * ncomps + c);
        const auto *if__269 = buffer.data(if_ + 269 * ncomps + c);
        const auto *if__270 = buffer.data(if_ + 270 * ncomps + c);
        const auto *if__271 = buffer.data(if_ + 271 * ncomps + c);
        const auto *if__272 = buffer.data(if_ + 272 * ncomps + c);
        const auto *if__273 = buffer.data(if_ + 273 * ncomps + c);
        const auto *if__274 = buffer.data(if_ + 274 * ncomps + c);
        const auto *if__275 = buffer.data(if_ + 275 * ncomps + c);
        const auto *if__276 = buffer.data(if_ + 276 * ncomps + c);
        const auto *if__277 = buffer.data(if_ + 277 * ncomps + c);
        const auto *if__278 = buffer.data(if_ + 278 * ncomps + c);
        const auto *if__279 = buffer.data(if_ + 279 * ncomps + c);

        const auto *kf_195 = buffer.data(kf + 195 * ncomps + c);
        const auto *kf_196 = buffer.data(kf + 196 * ncomps + c);
        const auto *kf_197 = buffer.data(kf + 197 * ncomps + c);
        const auto *kf_198 = buffer.data(kf + 198 * ncomps + c);
        const auto *kf_199 = buffer.data(kf + 199 * ncomps + c);
        const auto *kf_200 = buffer.data(kf + 200 * ncomps + c);
        const auto *kf_201 = buffer.data(kf + 201 * ncomps + c);
        const auto *kf_202 = buffer.data(kf + 202 * ncomps + c);
        const auto *kf_203 = buffer.data(kf + 203 * ncomps + c);
        const auto *kf_204 = buffer.data(kf + 204 * ncomps + c);
        const auto *kf_205 = buffer.data(kf + 205 * ncomps + c);
        const auto *kf_206 = buffer.data(kf + 206 * ncomps + c);
        const auto *kf_207 = buffer.data(kf + 207 * ncomps + c);
        const auto *kf_208 = buffer.data(kf + 208 * ncomps + c);
        const auto *kf_209 = buffer.data(kf + 209 * ncomps + c);
        const auto *kf_210 = buffer.data(kf + 210 * ncomps + c);
        const auto *kf_211 = buffer.data(kf + 211 * ncomps + c);
        const auto *kf_212 = buffer.data(kf + 212 * ncomps + c);
        const auto *kf_213 = buffer.data(kf + 213 * ncomps + c);
        const auto *kf_214 = buffer.data(kf + 214 * ncomps + c);
        const auto *kf_215 = buffer.data(kf + 215 * ncomps + c);
        const auto *kf_216 = buffer.data(kf + 216 * ncomps + c);
        const auto *kf_217 = buffer.data(kf + 217 * ncomps + c);
        const auto *kf_218 = buffer.data(kf + 218 * ncomps + c);
        const auto *kf_219 = buffer.data(kf + 219 * ncomps + c);
        const auto *kf_220 = buffer.data(kf + 220 * ncomps + c);
        const auto *kf_221 = buffer.data(kf + 221 * ncomps + c);
        const auto *kf_222 = buffer.data(kf + 222 * ncomps + c);
        const auto *kf_223 = buffer.data(kf + 223 * ncomps + c);
        const auto *kf_224 = buffer.data(kf + 224 * ncomps + c);
        const auto *kf_225 = buffer.data(kf + 225 * ncomps + c);
        const auto *kf_226 = buffer.data(kf + 226 * ncomps + c);
        const auto *kf_227 = buffer.data(kf + 227 * ncomps + c);
        const auto *kf_228 = buffer.data(kf + 228 * ncomps + c);
        const auto *kf_229 = buffer.data(kf + 229 * ncomps + c);
        const auto *kf_230 = buffer.data(kf + 230 * ncomps + c);
        const auto *kf_231 = buffer.data(kf + 231 * ncomps + c);
        const auto *kf_232 = buffer.data(kf + 232 * ncomps + c);
        const auto *kf_233 = buffer.data(kf + 233 * ncomps + c);
        const auto *kf_234 = buffer.data(kf + 234 * ncomps + c);
        const auto *kf_235 = buffer.data(kf + 235 * ncomps + c);
        const auto *kf_236 = buffer.data(kf + 236 * ncomps + c);
        const auto *kf_237 = buffer.data(kf + 237 * ncomps + c);
        const auto *kf_238 = buffer.data(kf + 238 * ncomps + c);
        const auto *kf_239 = buffer.data(kf + 239 * ncomps + c);
        const auto *kf_240 = buffer.data(kf + 240 * ncomps + c);
        const auto *kf_241 = buffer.data(kf + 241 * ncomps + c);
        const auto *kf_242 = buffer.data(kf + 242 * ncomps + c);
        const auto *kf_243 = buffer.data(kf + 243 * ncomps + c);
        const auto *kf_244 = buffer.data(kf + 244 * ncomps + c);
        const auto *kf_245 = buffer.data(kf + 245 * ncomps + c);
        const auto *kf_246 = buffer.data(kf + 246 * ncomps + c);
        const auto *kf_247 = buffer.data(kf + 247 * ncomps + c);
        const auto *kf_248 = buffer.data(kf + 248 * ncomps + c);
        const auto *kf_249 = buffer.data(kf + 249 * ncomps + c);
        const auto *kf_250 = buffer.data(kf + 250 * ncomps + c);
        const auto *kf_251 = buffer.data(kf + 251 * ncomps + c);
        const auto *kf_252 = buffer.data(kf + 252 * ncomps + c);
        const auto *kf_253 = buffer.data(kf + 253 * ncomps + c);
        const auto *kf_254 = buffer.data(kf + 254 * ncomps + c);
        const auto *kf_255 = buffer.data(kf + 255 * ncomps + c);
        const auto *kf_256 = buffer.data(kf + 256 * ncomps + c);
        const auto *kf_257 = buffer.data(kf + 257 * ncomps + c);
        const auto *kf_258 = buffer.data(kf + 258 * ncomps + c);
        const auto *kf_259 = buffer.data(kf + 259 * ncomps + c);
        const auto *kf_260 = buffer.data(kf + 260 * ncomps + c);
        const auto *kf_261 = buffer.data(kf + 261 * ncomps + c);
        const auto *kf_262 = buffer.data(kf + 262 * ncomps + c);
        const auto *kf_263 = buffer.data(kf + 263 * ncomps + c);
        const auto *kf_264 = buffer.data(kf + 264 * ncomps + c);
        const auto *kf_265 = buffer.data(kf + 265 * ncomps + c);
        const auto *kf_266 = buffer.data(kf + 266 * ncomps + c);
        const auto *kf_267 = buffer.data(kf + 267 * ncomps + c);
        const auto *kf_268 = buffer.data(kf + 268 * ncomps + c);
        const auto *kf_269 = buffer.data(kf + 269 * ncomps + c);
        const auto *kf_270 = buffer.data(kf + 270 * ncomps + c);
        const auto *kf_271 = buffer.data(kf + 271 * ncomps + c);
        const auto *kf_272 = buffer.data(kf + 272 * ncomps + c);
        const auto *kf_273 = buffer.data(kf + 273 * ncomps + c);
        const auto *kf_274 = buffer.data(kf + 274 * ncomps + c);
        const auto *kf_275 = buffer.data(kf + 275 * ncomps + c);
        const auto *kf_276 = buffer.data(kf + 276 * ncomps + c);
        const auto *kf_277 = buffer.data(kf + 277 * ncomps + c);
        const auto *kf_278 = buffer.data(kf + 278 * ncomps + c);
        const auto *kf_279 = buffer.data(kf + 279 * ncomps + c);
        const auto *kf_286 = buffer.data(kf + 286 * ncomps + c);
        const auto *kf_287 = buffer.data(kf + 287 * ncomps + c);
        const auto *kf_288 = buffer.data(kf + 288 * ncomps + c);
        const auto *kf_289 = buffer.data(kf + 289 * ncomps + c);
        const auto *kf_296 = buffer.data(kf + 296 * ncomps + c);
        const auto *kf_297 = buffer.data(kf + 297 * ncomps + c);
        const auto *kf_298 = buffer.data(kf + 298 * ncomps + c);
        const auto *kf_299 = buffer.data(kf + 299 * ncomps + c);
        const auto *kf_306 = buffer.data(kf + 306 * ncomps + c);
        const auto *kf_307 = buffer.data(kf + 307 * ncomps + c);
        const auto *kf_308 = buffer.data(kf + 308 * ncomps + c);
        const auto *kf_309 = buffer.data(kf + 309 * ncomps + c);
        const auto *kf_316 = buffer.data(kf + 316 * ncomps + c);
        const auto *kf_317 = buffer.data(kf + 317 * ncomps + c);
        const auto *kf_318 = buffer.data(kf + 318 * ncomps + c);
        const auto *kf_319 = buffer.data(kf + 319 * ncomps + c);
        const auto *kf_326 = buffer.data(kf + 326 * ncomps + c);
        const auto *kf_327 = buffer.data(kf + 327 * ncomps + c);
        const auto *kf_328 = buffer.data(kf + 328 * ncomps + c);
        const auto *kf_329 = buffer.data(kf + 329 * ncomps + c);
        const auto *kf_336 = buffer.data(kf + 336 * ncomps + c);
        const auto *kf_337 = buffer.data(kf + 337 * ncomps + c);
        const auto *kf_338 = buffer.data(kf + 338 * ncomps + c);
        const auto *kf_339 = buffer.data(kf + 339 * ncomps + c);
        const auto *kf_346 = buffer.data(kf + 346 * ncomps + c);
        const auto *kf_347 = buffer.data(kf + 347 * ncomps + c);
        const auto *kf_348 = buffer.data(kf + 348 * ncomps + c);
        const auto *kf_349 = buffer.data(kf + 349 * ncomps + c);
        const auto *kf_359 = buffer.data(kf + 359 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, if__195, if__196, if__197, \
                         if__198, if__199, kf_195, kf_196, kf_197, kf_198, \
                         kf_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_x[k] * if__195[k]
                       + kf_195[k];

            t_291[k] = ab_x[k] * if__196[k]
                       + kf_196[k];

            t_292[k] = ab_x[k] * if__197[k]
                       + kf_197[k];

            t_293[k] = ab_x[k] * if__198[k]
                       + kf_198[k];

            t_294[k] = ab_x[k] * if__199[k]
                       + kf_199[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_y, ab_z, if__196, if__197, \
                         if__198, if__199, kf_256, kf_257, kf_258, kf_259, \
                         kf_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = ab_y[k] * if__196[k]
                       + kf_256[k];

            t_296[k] = ab_y[k] * if__197[k]
                       + kf_257[k];

            t_297[k] = ab_y[k] * if__198[k]
                       + kf_258[k];

            t_298[k] = ab_y[k] * if__199[k]
                       + kf_259[k];

            t_299[k] = ab_z[k] * if__199[k]
                       + kf_269[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, if__200, if__201, if__202, \
                         if__203, if__204, kf_200, kf_201, kf_202, kf_203, \
                         kf_204 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = ab_x[k] * if__200[k]
                       + kf_200[k];

            t_301[k] = ab_x[k] * if__201[k]
                       + kf_201[k];

            t_302[k] = ab_x[k] * if__202[k]
                       + kf_202[k];

            t_303[k] = ab_x[k] * if__203[k]
                       + kf_203[k];

            t_304[k] = ab_x[k] * if__204[k]
                       + kf_204[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, if__205, if__206, if__207, \
                         if__208, if__209, kf_205, kf_206, kf_207, kf_208, \
                         kf_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = ab_x[k] * if__205[k]
                       + kf_205[k];

            t_306[k] = ab_x[k] * if__206[k]
                       + kf_206[k];

            t_307[k] = ab_x[k] * if__207[k]
                       + kf_207[k];

            t_308[k] = ab_x[k] * if__208[k]
                       + kf_208[k];

            t_309[k] = ab_x[k] * if__209[k]
                       + kf_209[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_y, ab_z, if__206, if__207, \
                         if__208, if__209, kf_266, kf_267, kf_268, kf_269, \
                         kf_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = ab_y[k] * if__206[k]
                       + kf_266[k];

            t_311[k] = ab_y[k] * if__207[k]
                       + kf_267[k];

            t_312[k] = ab_y[k] * if__208[k]
                       + kf_268[k];

            t_313[k] = ab_y[k] * if__209[k]
                       + kf_269[k];

            t_314[k] = ab_z[k] * if__209[k]
                       + kf_279[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, if__210, if__211, if__212, \
                         if__213, if__214, kf_210, kf_211, kf_212, kf_213, \
                         kf_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = ab_x[k] * if__210[k]
                       + kf_210[k];

            t_316[k] = ab_x[k] * if__211[k]
                       + kf_211[k];

            t_317[k] = ab_x[k] * if__212[k]
                       + kf_212[k];

            t_318[k] = ab_x[k] * if__213[k]
                       + kf_213[k];

            t_319[k] = ab_x[k] * if__214[k]
                       + kf_214[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, if__215, if__216, if__217, \
                         if__218, if__219, kf_215, kf_216, kf_217, kf_218, \
                         kf_219 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = ab_x[k] * if__215[k]
                       + kf_215[k];

            t_321[k] = ab_x[k] * if__216[k]
                       + kf_216[k];

            t_322[k] = ab_x[k] * if__217[k]
                       + kf_217[k];

            t_323[k] = ab_x[k] * if__218[k]
                       + kf_218[k];

            t_324[k] = ab_x[k] * if__219[k]
                       + kf_219[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_y, ab_z, if__216, if__217, \
                         if__218, if__219, kf_286, kf_287, kf_288, kf_289, \
                         kf_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = ab_y[k] * if__216[k]
                       + kf_286[k];

            t_326[k] = ab_y[k] * if__217[k]
                       + kf_287[k];

            t_327[k] = ab_y[k] * if__218[k]
                       + kf_288[k];

            t_328[k] = ab_y[k] * if__219[k]
                       + kf_289[k];

            t_329[k] = ab_z[k] * if__219[k]
                       + kf_299[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, if__220, if__221, if__222, \
                         if__223, if__224, kf_220, kf_221, kf_222, kf_223, \
                         kf_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = ab_x[k] * if__220[k]
                       + kf_220[k];

            t_331[k] = ab_x[k] * if__221[k]
                       + kf_221[k];

            t_332[k] = ab_x[k] * if__222[k]
                       + kf_222[k];

            t_333[k] = ab_x[k] * if__223[k]
                       + kf_223[k];

            t_334[k] = ab_x[k] * if__224[k]
                       + kf_224[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, if__225, if__226, if__227, \
                         if__228, if__229, kf_225, kf_226, kf_227, kf_228, \
                         kf_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = ab_x[k] * if__225[k]
                       + kf_225[k];

            t_336[k] = ab_x[k] * if__226[k]
                       + kf_226[k];

            t_337[k] = ab_x[k] * if__227[k]
                       + kf_227[k];

            t_338[k] = ab_x[k] * if__228[k]
                       + kf_228[k];

            t_339[k] = ab_x[k] * if__229[k]
                       + kf_229[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_y, ab_z, if__226, if__227, \
                         if__228, if__229, kf_296, kf_297, kf_298, kf_299, \
                         kf_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = ab_y[k] * if__226[k]
                       + kf_296[k];

            t_341[k] = ab_y[k] * if__227[k]
                       + kf_297[k];

            t_342[k] = ab_y[k] * if__228[k]
                       + kf_298[k];

            t_343[k] = ab_y[k] * if__229[k]
                       + kf_299[k];

            t_344[k] = ab_z[k] * if__229[k]
                       + kf_309[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, if__230, if__231, if__232, \
                         if__233, if__234, kf_230, kf_231, kf_232, kf_233, \
                         kf_234 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = ab_x[k] * if__230[k]
                       + kf_230[k];

            t_346[k] = ab_x[k] * if__231[k]
                       + kf_231[k];

            t_347[k] = ab_x[k] * if__232[k]
                       + kf_232[k];

            t_348[k] = ab_x[k] * if__233[k]
                       + kf_233[k];

            t_349[k] = ab_x[k] * if__234[k]
                       + kf_234[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, if__235, if__236, if__237, \
                         if__238, if__239, kf_235, kf_236, kf_237, kf_238, \
                         kf_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = ab_x[k] * if__235[k]
                       + kf_235[k];

            t_351[k] = ab_x[k] * if__236[k]
                       + kf_236[k];

            t_352[k] = ab_x[k] * if__237[k]
                       + kf_237[k];

            t_353[k] = ab_x[k] * if__238[k]
                       + kf_238[k];

            t_354[k] = ab_x[k] * if__239[k]
                       + kf_239[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_y, ab_z, if__236, if__237, \
                         if__238, if__239, kf_306, kf_307, kf_308, kf_309, \
                         kf_319 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = ab_y[k] * if__236[k]
                       + kf_306[k];

            t_356[k] = ab_y[k] * if__237[k]
                       + kf_307[k];

            t_357[k] = ab_y[k] * if__238[k]
                       + kf_308[k];

            t_358[k] = ab_y[k] * if__239[k]
                       + kf_309[k];

            t_359[k] = ab_z[k] * if__239[k]
                       + kf_319[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, if__240, if__241, if__242, \
                         if__243, if__244, kf_240, kf_241, kf_242, kf_243, \
                         kf_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * if__240[k]
                       + kf_240[k];

            t_361[k] = ab_x[k] * if__241[k]
                       + kf_241[k];

            t_362[k] = ab_x[k] * if__242[k]
                       + kf_242[k];

            t_363[k] = ab_x[k] * if__243[k]
                       + kf_243[k];

            t_364[k] = ab_x[k] * if__244[k]
                       + kf_244[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, if__245, if__246, if__247, \
                         if__248, if__249, kf_245, kf_246, kf_247, kf_248, \
                         kf_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * if__245[k]
                       + kf_245[k];

            t_366[k] = ab_x[k] * if__246[k]
                       + kf_246[k];

            t_367[k] = ab_x[k] * if__247[k]
                       + kf_247[k];

            t_368[k] = ab_x[k] * if__248[k]
                       + kf_248[k];

            t_369[k] = ab_x[k] * if__249[k]
                       + kf_249[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_y, ab_z, if__246, if__247, \
                         if__248, if__249, kf_316, kf_317, kf_318, kf_319, \
                         kf_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_y[k] * if__246[k]
                       + kf_316[k];

            t_371[k] = ab_y[k] * if__247[k]
                       + kf_317[k];

            t_372[k] = ab_y[k] * if__248[k]
                       + kf_318[k];

            t_373[k] = ab_y[k] * if__249[k]
                       + kf_319[k];

            t_374[k] = ab_z[k] * if__249[k]
                       + kf_329[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, if__250, if__251, if__252, \
                         if__253, if__254, kf_250, kf_251, kf_252, kf_253, \
                         kf_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = ab_x[k] * if__250[k]
                       + kf_250[k];

            t_376[k] = ab_x[k] * if__251[k]
                       + kf_251[k];

            t_377[k] = ab_x[k] * if__252[k]
                       + kf_252[k];

            t_378[k] = ab_x[k] * if__253[k]
                       + kf_253[k];

            t_379[k] = ab_x[k] * if__254[k]
                       + kf_254[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, if__255, if__256, if__257, \
                         if__258, if__259, kf_255, kf_256, kf_257, kf_258, \
                         kf_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = ab_x[k] * if__255[k]
                       + kf_255[k];

            t_381[k] = ab_x[k] * if__256[k]
                       + kf_256[k];

            t_382[k] = ab_x[k] * if__257[k]
                       + kf_257[k];

            t_383[k] = ab_x[k] * if__258[k]
                       + kf_258[k];

            t_384[k] = ab_x[k] * if__259[k]
                       + kf_259[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_y, ab_z, if__256, if__257, \
                         if__258, if__259, kf_326, kf_327, kf_328, kf_329, \
                         kf_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = ab_y[k] * if__256[k]
                       + kf_326[k];

            t_386[k] = ab_y[k] * if__257[k]
                       + kf_327[k];

            t_387[k] = ab_y[k] * if__258[k]
                       + kf_328[k];

            t_388[k] = ab_y[k] * if__259[k]
                       + kf_329[k];

            t_389[k] = ab_z[k] * if__259[k]
                       + kf_339[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, if__260, if__261, if__262, \
                         if__263, if__264, kf_260, kf_261, kf_262, kf_263, \
                         kf_264 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = ab_x[k] * if__260[k]
                       + kf_260[k];

            t_391[k] = ab_x[k] * if__261[k]
                       + kf_261[k];

            t_392[k] = ab_x[k] * if__262[k]
                       + kf_262[k];

            t_393[k] = ab_x[k] * if__263[k]
                       + kf_263[k];

            t_394[k] = ab_x[k] * if__264[k]
                       + kf_264[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, if__265, if__266, if__267, \
                         if__268, if__269, kf_265, kf_266, kf_267, kf_268, \
                         kf_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = ab_x[k] * if__265[k]
                       + kf_265[k];

            t_396[k] = ab_x[k] * if__266[k]
                       + kf_266[k];

            t_397[k] = ab_x[k] * if__267[k]
                       + kf_267[k];

            t_398[k] = ab_x[k] * if__268[k]
                       + kf_268[k];

            t_399[k] = ab_x[k] * if__269[k]
                       + kf_269[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_y, ab_z, if__266, if__267, \
                         if__268, if__269, kf_336, kf_337, kf_338, kf_339, \
                         kf_349 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = ab_y[k] * if__266[k]
                       + kf_336[k];

            t_401[k] = ab_y[k] * if__267[k]
                       + kf_337[k];

            t_402[k] = ab_y[k] * if__268[k]
                       + kf_338[k];

            t_403[k] = ab_y[k] * if__269[k]
                       + kf_339[k];

            t_404[k] = ab_z[k] * if__269[k]
                       + kf_349[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, if__270, if__271, if__272, \
                         if__273, if__274, kf_270, kf_271, kf_272, kf_273, \
                         kf_274 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = ab_x[k] * if__270[k]
                       + kf_270[k];

            t_406[k] = ab_x[k] * if__271[k]
                       + kf_271[k];

            t_407[k] = ab_x[k] * if__272[k]
                       + kf_272[k];

            t_408[k] = ab_x[k] * if__273[k]
                       + kf_273[k];

            t_409[k] = ab_x[k] * if__274[k]
                       + kf_274[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, if__275, if__276, if__277, \
                         if__278, if__279, kf_275, kf_276, kf_277, kf_278, \
                         kf_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = ab_x[k] * if__275[k]
                       + kf_275[k];

            t_411[k] = ab_x[k] * if__276[k]
                       + kf_276[k];

            t_412[k] = ab_x[k] * if__277[k]
                       + kf_277[k];

            t_413[k] = ab_x[k] * if__278[k]
                       + kf_278[k];

            t_414[k] = ab_x[k] * if__279[k]
                       + kf_279[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_y, ab_z, if__276, if__277, \
                         if__278, if__279, kf_346, kf_347, kf_348, kf_349, \
                         kf_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = ab_y[k] * if__276[k]
                       + kf_346[k];

            t_416[k] = ab_y[k] * if__277[k]
                       + kf_347[k];

            t_417[k] = ab_y[k] * if__278[k]
                       + kf_348[k];

            t_418[k] = ab_y[k] * if__279[k]
                       + kf_349[k];

            t_419[k] = ab_z[k] * if__279[k]
                       + kf_359[k];
        }
    }
}

auto
compute_hrr_ig(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t if_, const size_t kf, const size_t ncomps,
               const size_t nmax) -> void
{
    compute_hrr_ig_piece0(buffer, coordinates, target, if_, kf, ncomps, nmax);

    compute_hrr_ig_piece1(buffer, coordinates, target, if_, kf, ncomps, nmax);

    compute_hrr_ig_piece2(buffer, coordinates, target, if_, kf, ncomps, nmax);
}

}  // namespace simdtrf
