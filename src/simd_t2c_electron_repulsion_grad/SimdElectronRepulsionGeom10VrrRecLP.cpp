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


#include "SimdElectronRepulsionGeom10VrrRecLP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_lp_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t kp, const size_t mp,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);
    auto *t_9 = buffer.data(target + 9);
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);
    auto *t_45 = buffer.data(target + 45);
    auto *t_46 = buffer.data(target + 46);
    auto *t_47 = buffer.data(target + 47);
    auto *t_48 = buffer.data(target + 48);
    auto *t_49 = buffer.data(target + 49);
    auto *t_50 = buffer.data(target + 50);
    auto *t_51 = buffer.data(target + 51);
    auto *t_52 = buffer.data(target + 52);
    auto *t_53 = buffer.data(target + 53);
    auto *t_54 = buffer.data(target + 54);
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);
    auto *t_100 = buffer.data(target + 100);
    auto *t_101 = buffer.data(target + 101);
    auto *t_102 = buffer.data(target + 102);
    auto *t_103 = buffer.data(target + 103);
    auto *t_104 = buffer.data(target + 104);
    auto *t_105 = buffer.data(target + 105);
    auto *t_106 = buffer.data(target + 106);
    auto *t_107 = buffer.data(target + 107);
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);
    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);
    const auto *kp_49 = buffer.data(kp + 49);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_51 = buffer.data(kp + 51);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_54 = buffer.data(kp + 54);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_57 = buffer.data(kp + 57);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_60 = buffer.data(kp + 60);
    const auto *kp_61 = buffer.data(kp + 61);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_63 = buffer.data(kp + 63);
    const auto *kp_64 = buffer.data(kp + 64);
    const auto *kp_65 = buffer.data(kp + 65);
    const auto *kp_66 = buffer.data(kp + 66);
    const auto *kp_67 = buffer.data(kp + 67);
    const auto *kp_68 = buffer.data(kp + 68);
    const auto *kp_69 = buffer.data(kp + 69);
    const auto *kp_70 = buffer.data(kp + 70);
    const auto *kp_71 = buffer.data(kp + 71);
    const auto *kp_72 = buffer.data(kp + 72);
    const auto *kp_73 = buffer.data(kp + 73);
    const auto *kp_74 = buffer.data(kp + 74);
    const auto *kp_75 = buffer.data(kp + 75);
    const auto *kp_76 = buffer.data(kp + 76);
    const auto *kp_77 = buffer.data(kp + 77);
    const auto *kp_78 = buffer.data(kp + 78);
    const auto *kp_79 = buffer.data(kp + 79);
    const auto *kp_80 = buffer.data(kp + 80);
    const auto *kp_81 = buffer.data(kp + 81);
    const auto *kp_82 = buffer.data(kp + 82);
    const auto *kp_83 = buffer.data(kp + 83);
    const auto *kp_84 = buffer.data(kp + 84);
    const auto *kp_85 = buffer.data(kp + 85);
    const auto *kp_86 = buffer.data(kp + 86);
    const auto *kp_87 = buffer.data(kp + 87);
    const auto *kp_88 = buffer.data(kp + 88);
    const auto *kp_89 = buffer.data(kp + 89);
    const auto *kp_90 = buffer.data(kp + 90);
    const auto *kp_91 = buffer.data(kp + 91);
    const auto *kp_92 = buffer.data(kp + 92);
    const auto *kp_93 = buffer.data(kp + 93);
    const auto *kp_94 = buffer.data(kp + 94);
    const auto *kp_95 = buffer.data(kp + 95);
    const auto *kp_96 = buffer.data(kp + 96);
    const auto *kp_97 = buffer.data(kp + 97);
    const auto *kp_98 = buffer.data(kp + 98);
    const auto *kp_99 = buffer.data(kp + 99);
    const auto *kp_100 = buffer.data(kp + 100);
    const auto *kp_101 = buffer.data(kp + 101);
    const auto *kp_102 = buffer.data(kp + 102);
    const auto *kp_103 = buffer.data(kp + 103);
    const auto *kp_104 = buffer.data(kp + 104);
    const auto *kp_105 = buffer.data(kp + 105);
    const auto *kp_106 = buffer.data(kp + 106);
    const auto *kp_107 = buffer.data(kp + 107);

    const auto *mp_0 = buffer.data(mp + 0);
    const auto *mp_1 = buffer.data(mp + 1);
    const auto *mp_2 = buffer.data(mp + 2);
    const auto *mp_3 = buffer.data(mp + 3);
    const auto *mp_4 = buffer.data(mp + 4);
    const auto *mp_5 = buffer.data(mp + 5);
    const auto *mp_6 = buffer.data(mp + 6);
    const auto *mp_7 = buffer.data(mp + 7);
    const auto *mp_8 = buffer.data(mp + 8);
    const auto *mp_9 = buffer.data(mp + 9);
    const auto *mp_10 = buffer.data(mp + 10);
    const auto *mp_11 = buffer.data(mp + 11);
    const auto *mp_12 = buffer.data(mp + 12);
    const auto *mp_13 = buffer.data(mp + 13);
    const auto *mp_14 = buffer.data(mp + 14);
    const auto *mp_15 = buffer.data(mp + 15);
    const auto *mp_16 = buffer.data(mp + 16);
    const auto *mp_17 = buffer.data(mp + 17);
    const auto *mp_18 = buffer.data(mp + 18);
    const auto *mp_19 = buffer.data(mp + 19);
    const auto *mp_20 = buffer.data(mp + 20);
    const auto *mp_21 = buffer.data(mp + 21);
    const auto *mp_22 = buffer.data(mp + 22);
    const auto *mp_23 = buffer.data(mp + 23);
    const auto *mp_24 = buffer.data(mp + 24);
    const auto *mp_25 = buffer.data(mp + 25);
    const auto *mp_26 = buffer.data(mp + 26);
    const auto *mp_27 = buffer.data(mp + 27);
    const auto *mp_28 = buffer.data(mp + 28);
    const auto *mp_29 = buffer.data(mp + 29);
    const auto *mp_30 = buffer.data(mp + 30);
    const auto *mp_31 = buffer.data(mp + 31);
    const auto *mp_32 = buffer.data(mp + 32);
    const auto *mp_33 = buffer.data(mp + 33);
    const auto *mp_34 = buffer.data(mp + 34);
    const auto *mp_35 = buffer.data(mp + 35);
    const auto *mp_36 = buffer.data(mp + 36);
    const auto *mp_37 = buffer.data(mp + 37);
    const auto *mp_38 = buffer.data(mp + 38);
    const auto *mp_39 = buffer.data(mp + 39);
    const auto *mp_40 = buffer.data(mp + 40);
    const auto *mp_41 = buffer.data(mp + 41);
    const auto *mp_42 = buffer.data(mp + 42);
    const auto *mp_43 = buffer.data(mp + 43);
    const auto *mp_44 = buffer.data(mp + 44);
    const auto *mp_45 = buffer.data(mp + 45);
    const auto *mp_46 = buffer.data(mp + 46);
    const auto *mp_47 = buffer.data(mp + 47);
    const auto *mp_48 = buffer.data(mp + 48);
    const auto *mp_49 = buffer.data(mp + 49);
    const auto *mp_50 = buffer.data(mp + 50);
    const auto *mp_51 = buffer.data(mp + 51);
    const auto *mp_52 = buffer.data(mp + 52);
    const auto *mp_53 = buffer.data(mp + 53);
    const auto *mp_54 = buffer.data(mp + 54);
    const auto *mp_55 = buffer.data(mp + 55);
    const auto *mp_56 = buffer.data(mp + 56);
    const auto *mp_57 = buffer.data(mp + 57);
    const auto *mp_58 = buffer.data(mp + 58);
    const auto *mp_59 = buffer.data(mp + 59);
    const auto *mp_60 = buffer.data(mp + 60);
    const auto *mp_61 = buffer.data(mp + 61);
    const auto *mp_62 = buffer.data(mp + 62);
    const auto *mp_63 = buffer.data(mp + 63);
    const auto *mp_64 = buffer.data(mp + 64);
    const auto *mp_65 = buffer.data(mp + 65);
    const auto *mp_66 = buffer.data(mp + 66);
    const auto *mp_67 = buffer.data(mp + 67);
    const auto *mp_68 = buffer.data(mp + 68);
    const auto *mp_69 = buffer.data(mp + 69);
    const auto *mp_70 = buffer.data(mp + 70);
    const auto *mp_71 = buffer.data(mp + 71);
    const auto *mp_72 = buffer.data(mp + 72);
    const auto *mp_73 = buffer.data(mp + 73);
    const auto *mp_74 = buffer.data(mp + 74);
    const auto *mp_75 = buffer.data(mp + 75);
    const auto *mp_76 = buffer.data(mp + 76);
    const auto *mp_77 = buffer.data(mp + 77);
    const auto *mp_78 = buffer.data(mp + 78);
    const auto *mp_79 = buffer.data(mp + 79);
    const auto *mp_80 = buffer.data(mp + 80);
    const auto *mp_81 = buffer.data(mp + 81);
    const auto *mp_82 = buffer.data(mp + 82);
    const auto *mp_83 = buffer.data(mp + 83);
    const auto *mp_84 = buffer.data(mp + 84);
    const auto *mp_85 = buffer.data(mp + 85);
    const auto *mp_86 = buffer.data(mp + 86);
    const auto *mp_87 = buffer.data(mp + 87);
    const auto *mp_88 = buffer.data(mp + 88);
    const auto *mp_89 = buffer.data(mp + 89);
    const auto *mp_90 = buffer.data(mp + 90);
    const auto *mp_91 = buffer.data(mp + 91);
    const auto *mp_92 = buffer.data(mp + 92);
    const auto *mp_93 = buffer.data(mp + 93);
    const auto *mp_94 = buffer.data(mp + 94);
    const auto *mp_95 = buffer.data(mp + 95);
    const auto *mp_96 = buffer.data(mp + 96);
    const auto *mp_97 = buffer.data(mp + 97);
    const auto *mp_98 = buffer.data(mp + 98);
    const auto *mp_99 = buffer.data(mp + 99);
    const auto *mp_100 = buffer.data(mp + 100);
    const auto *mp_101 = buffer.data(mp + 101);
    const auto *mp_102 = buffer.data(mp + 102);
    const auto *mp_103 = buffer.data(mp + 103);
    const auto *mp_104 = buffer.data(mp + 104);
    const auto *mp_105 = buffer.data(mp + 105);
    const auto *mp_106 = buffer.data(mp + 106);
    const auto *mp_107 = buffer.data(mp + 107);
    const auto *mp_108 = buffer.data(mp + 108);
    const auto *mp_109 = buffer.data(mp + 109);
    const auto *mp_110 = buffer.data(mp + 110);
    const auto *mp_111 = buffer.data(mp + 111);
    const auto *mp_112 = buffer.data(mp + 112);
    const auto *mp_113 = buffer.data(mp + 113);
    const auto *mp_114 = buffer.data(mp + 114);
    const auto *mp_115 = buffer.data(mp + 115);
    const auto *mp_116 = buffer.data(mp + 116);
    const auto *mp_117 = buffer.data(mp + 117);
    const auto *mp_118 = buffer.data(mp + 118);
    const auto *mp_119 = buffer.data(mp + 119);
    const auto *mp_120 = buffer.data(mp + 120);
    const auto *mp_121 = buffer.data(mp + 121);
    const auto *mp_122 = buffer.data(mp + 122);
    const auto *mp_123 = buffer.data(mp + 123);
    const auto *mp_124 = buffer.data(mp + 124);
    const auto *mp_125 = buffer.data(mp + 125);
    const auto *mp_126 = buffer.data(mp + 126);
    const auto *mp_127 = buffer.data(mp + 127);
    const auto *mp_128 = buffer.data(mp + 128);
    const auto *mp_129 = buffer.data(mp + 129);
    const auto *mp_130 = buffer.data(mp + 130);
    const auto *mp_131 = buffer.data(mp + 131);
    const auto *mp_132 = buffer.data(mp + 132);
    const auto *mp_133 = buffer.data(mp + 133);
    const auto *mp_134 = buffer.data(mp + 134);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, kp_0, kp_1, kp_2, kp_3, kp_4, mp_0, mp_1, \
                         mp_2, mp_3, mp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -8.0 * kp_0[k]
                 + f_0 * mp_0[k];

        t_1[k] = -8.0 * kp_1[k]
                 + f_0 * mp_1[k];

        t_2[k] = -8.0 * kp_2[k]
                 + f_0 * mp_2[k];

        t_3[k] = -7.0 * kp_3[k]
                 + f_0 * mp_3[k];

        t_4[k] = -7.0 * kp_4[k]
                 + f_0 * mp_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, kp_5, kp_6, kp_7, kp_8, kp_9, mp_5, mp_6, \
                         mp_7, mp_8, mp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -7.0 * kp_5[k]
                 + f_0 * mp_5[k];

        t_6[k] = -7.0 * kp_6[k]
                 + f_0 * mp_6[k];

        t_7[k] = -7.0 * kp_7[k]
                 + f_0 * mp_7[k];

        t_8[k] = -7.0 * kp_8[k]
                 + f_0 * mp_8[k];

        t_9[k] = -6.0 * kp_9[k]
                 + f_0 * mp_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, kp_10, kp_11, kp_12, kp_13, kp_14, \
                         mp_10, mp_11, mp_12, mp_13, mp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -6.0 * kp_10[k]
                  + f_0 * mp_10[k];

        t_11[k] = -6.0 * kp_11[k]
                  + f_0 * mp_11[k];

        t_12[k] = -6.0 * kp_12[k]
                  + f_0 * mp_12[k];

        t_13[k] = -6.0 * kp_13[k]
                  + f_0 * mp_13[k];

        t_14[k] = -6.0 * kp_14[k]
                  + f_0 * mp_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, kp_15, kp_16, kp_17, kp_18, kp_19, \
                         mp_15, mp_16, mp_17, mp_18, mp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -6.0 * kp_15[k]
                  + f_0 * mp_15[k];

        t_16[k] = -6.0 * kp_16[k]
                  + f_0 * mp_16[k];

        t_17[k] = -6.0 * kp_17[k]
                  + f_0 * mp_17[k];

        t_18[k] = -5.0 * kp_18[k]
                  + f_0 * mp_18[k];

        t_19[k] = -5.0 * kp_19[k]
                  + f_0 * mp_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, kp_20, kp_21, kp_22, kp_23, kp_24, \
                         mp_20, mp_21, mp_22, mp_23, mp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -5.0 * kp_20[k]
                  + f_0 * mp_20[k];

        t_21[k] = -5.0 * kp_21[k]
                  + f_0 * mp_21[k];

        t_22[k] = -5.0 * kp_22[k]
                  + f_0 * mp_22[k];

        t_23[k] = -5.0 * kp_23[k]
                  + f_0 * mp_23[k];

        t_24[k] = -5.0 * kp_24[k]
                  + f_0 * mp_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, kp_25, kp_26, kp_27, kp_28, kp_29, \
                         mp_25, mp_26, mp_27, mp_28, mp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -5.0 * kp_25[k]
                  + f_0 * mp_25[k];

        t_26[k] = -5.0 * kp_26[k]
                  + f_0 * mp_26[k];

        t_27[k] = -5.0 * kp_27[k]
                  + f_0 * mp_27[k];

        t_28[k] = -5.0 * kp_28[k]
                  + f_0 * mp_28[k];

        t_29[k] = -5.0 * kp_29[k]
                  + f_0 * mp_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, kp_30, kp_31, kp_32, kp_33, kp_34, \
                         mp_30, mp_31, mp_32, mp_33, mp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -4.0 * kp_30[k]
                  + f_0 * mp_30[k];

        t_31[k] = -4.0 * kp_31[k]
                  + f_0 * mp_31[k];

        t_32[k] = -4.0 * kp_32[k]
                  + f_0 * mp_32[k];

        t_33[k] = -4.0 * kp_33[k]
                  + f_0 * mp_33[k];

        t_34[k] = -4.0 * kp_34[k]
                  + f_0 * mp_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, kp_35, kp_36, kp_37, kp_38, kp_39, \
                         mp_35, mp_36, mp_37, mp_38, mp_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -4.0 * kp_35[k]
                  + f_0 * mp_35[k];

        t_36[k] = -4.0 * kp_36[k]
                  + f_0 * mp_36[k];

        t_37[k] = -4.0 * kp_37[k]
                  + f_0 * mp_37[k];

        t_38[k] = -4.0 * kp_38[k]
                  + f_0 * mp_38[k];

        t_39[k] = -4.0 * kp_39[k]
                  + f_0 * mp_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, kp_40, kp_41, kp_42, kp_43, kp_44, \
                         mp_40, mp_41, mp_42, mp_43, mp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -4.0 * kp_40[k]
                  + f_0 * mp_40[k];

        t_41[k] = -4.0 * kp_41[k]
                  + f_0 * mp_41[k];

        t_42[k] = -4.0 * kp_42[k]
                  + f_0 * mp_42[k];

        t_43[k] = -4.0 * kp_43[k]
                  + f_0 * mp_43[k];

        t_44[k] = -4.0 * kp_44[k]
                  + f_0 * mp_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, kp_45, kp_46, kp_47, kp_48, kp_49, \
                         mp_45, mp_46, mp_47, mp_48, mp_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -3.0 * kp_45[k]
                  + f_0 * mp_45[k];

        t_46[k] = -3.0 * kp_46[k]
                  + f_0 * mp_46[k];

        t_47[k] = -3.0 * kp_47[k]
                  + f_0 * mp_47[k];

        t_48[k] = -3.0 * kp_48[k]
                  + f_0 * mp_48[k];

        t_49[k] = -3.0 * kp_49[k]
                  + f_0 * mp_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, kp_50, kp_51, kp_52, kp_53, kp_54, \
                         mp_50, mp_51, mp_52, mp_53, mp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -3.0 * kp_50[k]
                  + f_0 * mp_50[k];

        t_51[k] = -3.0 * kp_51[k]
                  + f_0 * mp_51[k];

        t_52[k] = -3.0 * kp_52[k]
                  + f_0 * mp_52[k];

        t_53[k] = -3.0 * kp_53[k]
                  + f_0 * mp_53[k];

        t_54[k] = -3.0 * kp_54[k]
                  + f_0 * mp_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, kp_55, kp_56, kp_57, kp_58, kp_59, \
                         mp_55, mp_56, mp_57, mp_58, mp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -3.0 * kp_55[k]
                  + f_0 * mp_55[k];

        t_56[k] = -3.0 * kp_56[k]
                  + f_0 * mp_56[k];

        t_57[k] = -3.0 * kp_57[k]
                  + f_0 * mp_57[k];

        t_58[k] = -3.0 * kp_58[k]
                  + f_0 * mp_58[k];

        t_59[k] = -3.0 * kp_59[k]
                  + f_0 * mp_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, kp_60, kp_61, kp_62, kp_63, kp_64, \
                         mp_60, mp_61, mp_62, mp_63, mp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * kp_60[k]
                  + f_0 * mp_60[k];

        t_61[k] = -3.0 * kp_61[k]
                  + f_0 * mp_61[k];

        t_62[k] = -3.0 * kp_62[k]
                  + f_0 * mp_62[k];

        t_63[k] = -2.0 * kp_63[k]
                  + f_0 * mp_63[k];

        t_64[k] = -2.0 * kp_64[k]
                  + f_0 * mp_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, kp_65, kp_66, kp_67, kp_68, kp_69, \
                         mp_65, mp_66, mp_67, mp_68, mp_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -2.0 * kp_65[k]
                  + f_0 * mp_65[k];

        t_66[k] = -2.0 * kp_66[k]
                  + f_0 * mp_66[k];

        t_67[k] = -2.0 * kp_67[k]
                  + f_0 * mp_67[k];

        t_68[k] = -2.0 * kp_68[k]
                  + f_0 * mp_68[k];

        t_69[k] = -2.0 * kp_69[k]
                  + f_0 * mp_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, kp_70, kp_71, kp_72, kp_73, kp_74, \
                         mp_70, mp_71, mp_72, mp_73, mp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * kp_70[k]
                  + f_0 * mp_70[k];

        t_71[k] = -2.0 * kp_71[k]
                  + f_0 * mp_71[k];

        t_72[k] = -2.0 * kp_72[k]
                  + f_0 * mp_72[k];

        t_73[k] = -2.0 * kp_73[k]
                  + f_0 * mp_73[k];

        t_74[k] = -2.0 * kp_74[k]
                  + f_0 * mp_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, kp_75, kp_76, kp_77, kp_78, kp_79, \
                         mp_75, mp_76, mp_77, mp_78, mp_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * kp_75[k]
                  + f_0 * mp_75[k];

        t_76[k] = -2.0 * kp_76[k]
                  + f_0 * mp_76[k];

        t_77[k] = -2.0 * kp_77[k]
                  + f_0 * mp_77[k];

        t_78[k] = -2.0 * kp_78[k]
                  + f_0 * mp_78[k];

        t_79[k] = -2.0 * kp_79[k]
                  + f_0 * mp_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, kp_80, kp_81, kp_82, kp_83, kp_84, \
                         mp_80, mp_81, mp_82, mp_83, mp_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -2.0 * kp_80[k]
                  + f_0 * mp_80[k];

        t_81[k] = -2.0 * kp_81[k]
                  + f_0 * mp_81[k];

        t_82[k] = -2.0 * kp_82[k]
                  + f_0 * mp_82[k];

        t_83[k] = -2.0 * kp_83[k]
                  + f_0 * mp_83[k];

        t_84[k] = -kp_84[k]
                  + f_0 * mp_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, kp_85, kp_86, kp_87, kp_88, kp_89, \
                         mp_85, mp_86, mp_87, mp_88, mp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -kp_85[k]
                  + f_0 * mp_85[k];

        t_86[k] = -kp_86[k]
                  + f_0 * mp_86[k];

        t_87[k] = -kp_87[k]
                  + f_0 * mp_87[k];

        t_88[k] = -kp_88[k]
                  + f_0 * mp_88[k];

        t_89[k] = -kp_89[k]
                  + f_0 * mp_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, kp_90, kp_91, kp_92, kp_93, kp_94, \
                         mp_90, mp_91, mp_92, mp_93, mp_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -kp_90[k]
                  + f_0 * mp_90[k];

        t_91[k] = -kp_91[k]
                  + f_0 * mp_91[k];

        t_92[k] = -kp_92[k]
                  + f_0 * mp_92[k];

        t_93[k] = -kp_93[k]
                  + f_0 * mp_93[k];

        t_94[k] = -kp_94[k]
                  + f_0 * mp_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, kp_95, kp_96, kp_97, kp_98, kp_99, \
                         mp_95, mp_96, mp_97, mp_98, mp_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -kp_95[k]
                  + f_0 * mp_95[k];

        t_96[k] = -kp_96[k]
                  + f_0 * mp_96[k];

        t_97[k] = -kp_97[k]
                  + f_0 * mp_97[k];

        t_98[k] = -kp_98[k]
                  + f_0 * mp_98[k];

        t_99[k] = -kp_99[k]
                  + f_0 * mp_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, kp_100, kp_101, kp_102, kp_103, \
                         kp_104, mp_100, mp_101, mp_102, mp_103, \
                         mp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -kp_100[k]
                   + f_0 * mp_100[k];

        t_101[k] = -kp_101[k]
                   + f_0 * mp_101[k];

        t_102[k] = -kp_102[k]
                   + f_0 * mp_102[k];

        t_103[k] = -kp_103[k]
                   + f_0 * mp_103[k];

        t_104[k] = -kp_104[k]
                   + f_0 * mp_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, t_110, kp_105, kp_106, kp_107, \
                         mp_105, mp_106, mp_107, mp_108, mp_109, \
                         mp_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -kp_105[k]
                   + f_0 * mp_105[k];

        t_106[k] = -kp_106[k]
                   + f_0 * mp_106[k];

        t_107[k] = -kp_107[k]
                   + f_0 * mp_107[k];

        t_108[k] = f_0 * mp_108[k];

        t_109[k] = f_0 * mp_109[k];

        t_110[k] = f_0 * mp_110[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, t_116, t_117, t_118, mp_111, \
                         mp_112, mp_113, mp_114, mp_115, mp_116, mp_117, \
                         mp_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_0 * mp_111[k];

        t_112[k] = f_0 * mp_112[k];

        t_113[k] = f_0 * mp_113[k];

        t_114[k] = f_0 * mp_114[k];

        t_115[k] = f_0 * mp_115[k];

        t_116[k] = f_0 * mp_116[k];

        t_117[k] = f_0 * mp_117[k];

        t_118[k] = f_0 * mp_118[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, t_124, t_125, t_126, mp_119, \
                         mp_120, mp_121, mp_122, mp_123, mp_124, mp_125, \
                         mp_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_0 * mp_119[k];

        t_120[k] = f_0 * mp_120[k];

        t_121[k] = f_0 * mp_121[k];

        t_122[k] = f_0 * mp_122[k];

        t_123[k] = f_0 * mp_123[k];

        t_124[k] = f_0 * mp_124[k];

        t_125[k] = f_0 * mp_125[k];

        t_126[k] = f_0 * mp_126[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, t_132, t_133, t_134, mp_127, \
                         mp_128, mp_129, mp_130, mp_131, mp_132, mp_133, \
                         mp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_0 * mp_127[k];

        t_128[k] = f_0 * mp_128[k];

        t_129[k] = f_0 * mp_129[k];

        t_130[k] = f_0 * mp_130[k];

        t_131[k] = f_0 * mp_131[k];

        t_132[k] = f_0 * mp_132[k];

        t_133[k] = f_0 * mp_133[k];

        t_134[k] = f_0 * mp_134[k];
    }
}

auto
compute_prim_geom_10_lp_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t kp, const size_t mp,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);
    auto *t_9 = buffer.data(target + 9);
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);
    auto *t_45 = buffer.data(target + 45);
    auto *t_46 = buffer.data(target + 46);
    auto *t_47 = buffer.data(target + 47);
    auto *t_48 = buffer.data(target + 48);
    auto *t_49 = buffer.data(target + 49);
    auto *t_50 = buffer.data(target + 50);
    auto *t_51 = buffer.data(target + 51);
    auto *t_52 = buffer.data(target + 52);
    auto *t_53 = buffer.data(target + 53);
    auto *t_54 = buffer.data(target + 54);
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);
    auto *t_100 = buffer.data(target + 100);
    auto *t_101 = buffer.data(target + 101);
    auto *t_102 = buffer.data(target + 102);
    auto *t_103 = buffer.data(target + 103);
    auto *t_104 = buffer.data(target + 104);
    auto *t_105 = buffer.data(target + 105);
    auto *t_106 = buffer.data(target + 106);
    auto *t_107 = buffer.data(target + 107);
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);
    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);
    const auto *kp_49 = buffer.data(kp + 49);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_51 = buffer.data(kp + 51);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_54 = buffer.data(kp + 54);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_57 = buffer.data(kp + 57);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_60 = buffer.data(kp + 60);
    const auto *kp_61 = buffer.data(kp + 61);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_63 = buffer.data(kp + 63);
    const auto *kp_64 = buffer.data(kp + 64);
    const auto *kp_65 = buffer.data(kp + 65);
    const auto *kp_66 = buffer.data(kp + 66);
    const auto *kp_67 = buffer.data(kp + 67);
    const auto *kp_68 = buffer.data(kp + 68);
    const auto *kp_69 = buffer.data(kp + 69);
    const auto *kp_70 = buffer.data(kp + 70);
    const auto *kp_71 = buffer.data(kp + 71);
    const auto *kp_72 = buffer.data(kp + 72);
    const auto *kp_73 = buffer.data(kp + 73);
    const auto *kp_74 = buffer.data(kp + 74);
    const auto *kp_75 = buffer.data(kp + 75);
    const auto *kp_76 = buffer.data(kp + 76);
    const auto *kp_77 = buffer.data(kp + 77);
    const auto *kp_78 = buffer.data(kp + 78);
    const auto *kp_79 = buffer.data(kp + 79);
    const auto *kp_80 = buffer.data(kp + 80);
    const auto *kp_81 = buffer.data(kp + 81);
    const auto *kp_82 = buffer.data(kp + 82);
    const auto *kp_83 = buffer.data(kp + 83);
    const auto *kp_84 = buffer.data(kp + 84);
    const auto *kp_85 = buffer.data(kp + 85);
    const auto *kp_86 = buffer.data(kp + 86);
    const auto *kp_87 = buffer.data(kp + 87);
    const auto *kp_88 = buffer.data(kp + 88);
    const auto *kp_89 = buffer.data(kp + 89);
    const auto *kp_90 = buffer.data(kp + 90);
    const auto *kp_91 = buffer.data(kp + 91);
    const auto *kp_92 = buffer.data(kp + 92);
    const auto *kp_93 = buffer.data(kp + 93);
    const auto *kp_94 = buffer.data(kp + 94);
    const auto *kp_95 = buffer.data(kp + 95);
    const auto *kp_96 = buffer.data(kp + 96);
    const auto *kp_97 = buffer.data(kp + 97);
    const auto *kp_98 = buffer.data(kp + 98);
    const auto *kp_99 = buffer.data(kp + 99);
    const auto *kp_100 = buffer.data(kp + 100);
    const auto *kp_101 = buffer.data(kp + 101);
    const auto *kp_102 = buffer.data(kp + 102);
    const auto *kp_103 = buffer.data(kp + 103);
    const auto *kp_104 = buffer.data(kp + 104);
    const auto *kp_105 = buffer.data(kp + 105);
    const auto *kp_106 = buffer.data(kp + 106);
    const auto *kp_107 = buffer.data(kp + 107);

    const auto *mp_3 = buffer.data(mp + 3);
    const auto *mp_4 = buffer.data(mp + 4);
    const auto *mp_5 = buffer.data(mp + 5);
    const auto *mp_9 = buffer.data(mp + 9);
    const auto *mp_10 = buffer.data(mp + 10);
    const auto *mp_11 = buffer.data(mp + 11);
    const auto *mp_12 = buffer.data(mp + 12);
    const auto *mp_13 = buffer.data(mp + 13);
    const auto *mp_14 = buffer.data(mp + 14);
    const auto *mp_18 = buffer.data(mp + 18);
    const auto *mp_19 = buffer.data(mp + 19);
    const auto *mp_20 = buffer.data(mp + 20);
    const auto *mp_21 = buffer.data(mp + 21);
    const auto *mp_22 = buffer.data(mp + 22);
    const auto *mp_23 = buffer.data(mp + 23);
    const auto *mp_24 = buffer.data(mp + 24);
    const auto *mp_25 = buffer.data(mp + 25);
    const auto *mp_26 = buffer.data(mp + 26);
    const auto *mp_30 = buffer.data(mp + 30);
    const auto *mp_31 = buffer.data(mp + 31);
    const auto *mp_32 = buffer.data(mp + 32);
    const auto *mp_33 = buffer.data(mp + 33);
    const auto *mp_34 = buffer.data(mp + 34);
    const auto *mp_35 = buffer.data(mp + 35);
    const auto *mp_36 = buffer.data(mp + 36);
    const auto *mp_37 = buffer.data(mp + 37);
    const auto *mp_38 = buffer.data(mp + 38);
    const auto *mp_39 = buffer.data(mp + 39);
    const auto *mp_40 = buffer.data(mp + 40);
    const auto *mp_41 = buffer.data(mp + 41);
    const auto *mp_45 = buffer.data(mp + 45);
    const auto *mp_46 = buffer.data(mp + 46);
    const auto *mp_47 = buffer.data(mp + 47);
    const auto *mp_48 = buffer.data(mp + 48);
    const auto *mp_49 = buffer.data(mp + 49);
    const auto *mp_50 = buffer.data(mp + 50);
    const auto *mp_51 = buffer.data(mp + 51);
    const auto *mp_52 = buffer.data(mp + 52);
    const auto *mp_53 = buffer.data(mp + 53);
    const auto *mp_54 = buffer.data(mp + 54);
    const auto *mp_55 = buffer.data(mp + 55);
    const auto *mp_56 = buffer.data(mp + 56);
    const auto *mp_57 = buffer.data(mp + 57);
    const auto *mp_58 = buffer.data(mp + 58);
    const auto *mp_59 = buffer.data(mp + 59);
    const auto *mp_63 = buffer.data(mp + 63);
    const auto *mp_64 = buffer.data(mp + 64);
    const auto *mp_65 = buffer.data(mp + 65);
    const auto *mp_66 = buffer.data(mp + 66);
    const auto *mp_67 = buffer.data(mp + 67);
    const auto *mp_68 = buffer.data(mp + 68);
    const auto *mp_69 = buffer.data(mp + 69);
    const auto *mp_70 = buffer.data(mp + 70);
    const auto *mp_71 = buffer.data(mp + 71);
    const auto *mp_72 = buffer.data(mp + 72);
    const auto *mp_73 = buffer.data(mp + 73);
    const auto *mp_74 = buffer.data(mp + 74);
    const auto *mp_75 = buffer.data(mp + 75);
    const auto *mp_76 = buffer.data(mp + 76);
    const auto *mp_77 = buffer.data(mp + 77);
    const auto *mp_78 = buffer.data(mp + 78);
    const auto *mp_79 = buffer.data(mp + 79);
    const auto *mp_80 = buffer.data(mp + 80);
    const auto *mp_84 = buffer.data(mp + 84);
    const auto *mp_85 = buffer.data(mp + 85);
    const auto *mp_86 = buffer.data(mp + 86);
    const auto *mp_87 = buffer.data(mp + 87);
    const auto *mp_88 = buffer.data(mp + 88);
    const auto *mp_89 = buffer.data(mp + 89);
    const auto *mp_90 = buffer.data(mp + 90);
    const auto *mp_91 = buffer.data(mp + 91);
    const auto *mp_92 = buffer.data(mp + 92);
    const auto *mp_93 = buffer.data(mp + 93);
    const auto *mp_94 = buffer.data(mp + 94);
    const auto *mp_95 = buffer.data(mp + 95);
    const auto *mp_96 = buffer.data(mp + 96);
    const auto *mp_97 = buffer.data(mp + 97);
    const auto *mp_98 = buffer.data(mp + 98);
    const auto *mp_99 = buffer.data(mp + 99);
    const auto *mp_100 = buffer.data(mp + 100);
    const auto *mp_101 = buffer.data(mp + 101);
    const auto *mp_102 = buffer.data(mp + 102);
    const auto *mp_103 = buffer.data(mp + 103);
    const auto *mp_104 = buffer.data(mp + 104);
    const auto *mp_108 = buffer.data(mp + 108);
    const auto *mp_109 = buffer.data(mp + 109);
    const auto *mp_110 = buffer.data(mp + 110);
    const auto *mp_111 = buffer.data(mp + 111);
    const auto *mp_112 = buffer.data(mp + 112);
    const auto *mp_113 = buffer.data(mp + 113);
    const auto *mp_114 = buffer.data(mp + 114);
    const auto *mp_115 = buffer.data(mp + 115);
    const auto *mp_116 = buffer.data(mp + 116);
    const auto *mp_117 = buffer.data(mp + 117);
    const auto *mp_118 = buffer.data(mp + 118);
    const auto *mp_119 = buffer.data(mp + 119);
    const auto *mp_120 = buffer.data(mp + 120);
    const auto *mp_121 = buffer.data(mp + 121);
    const auto *mp_122 = buffer.data(mp + 122);
    const auto *mp_123 = buffer.data(mp + 123);
    const auto *mp_124 = buffer.data(mp + 124);
    const auto *mp_125 = buffer.data(mp + 125);
    const auto *mp_126 = buffer.data(mp + 126);
    const auto *mp_127 = buffer.data(mp + 127);
    const auto *mp_128 = buffer.data(mp + 128);
    const auto *mp_129 = buffer.data(mp + 129);
    const auto *mp_130 = buffer.data(mp + 130);
    const auto *mp_131 = buffer.data(mp + 131);
    const auto *mp_135 = buffer.data(mp + 135);
    const auto *mp_136 = buffer.data(mp + 136);
    const auto *mp_137 = buffer.data(mp + 137);
    const auto *mp_138 = buffer.data(mp + 138);
    const auto *mp_139 = buffer.data(mp + 139);
    const auto *mp_140 = buffer.data(mp + 140);
    const auto *mp_141 = buffer.data(mp + 141);
    const auto *mp_142 = buffer.data(mp + 142);
    const auto *mp_143 = buffer.data(mp + 143);
    const auto *mp_144 = buffer.data(mp + 144);
    const auto *mp_145 = buffer.data(mp + 145);
    const auto *mp_146 = buffer.data(mp + 146);
    const auto *mp_147 = buffer.data(mp + 147);
    const auto *mp_148 = buffer.data(mp + 148);
    const auto *mp_149 = buffer.data(mp + 149);
    const auto *mp_150 = buffer.data(mp + 150);
    const auto *mp_151 = buffer.data(mp + 151);
    const auto *mp_152 = buffer.data(mp + 152);
    const auto *mp_153 = buffer.data(mp + 153);
    const auto *mp_154 = buffer.data(mp + 154);
    const auto *mp_155 = buffer.data(mp + 155);
    const auto *mp_156 = buffer.data(mp + 156);
    const auto *mp_157 = buffer.data(mp + 157);
    const auto *mp_158 = buffer.data(mp + 158);
    const auto *mp_159 = buffer.data(mp + 159);
    const auto *mp_160 = buffer.data(mp + 160);
    const auto *mp_161 = buffer.data(mp + 161);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, kp_0, kp_1, kp_2, mp_3, mp_4, mp_5, \
                         mp_9, mp_10, mp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * mp_3[k];

        t_1[k] = f_0 * mp_4[k];

        t_2[k] = f_0 * mp_5[k];

        t_3[k] = -kp_0[k]
                 + f_0 * mp_9[k];

        t_4[k] = -kp_1[k]
                 + f_0 * mp_10[k];

        t_5[k] = -kp_2[k]
                 + f_0 * mp_11[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, kp_3, kp_4, kp_5, mp_12, mp_13, \
                         mp_14, mp_18, mp_19, mp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * mp_12[k];

        t_7[k] = f_0 * mp_13[k];

        t_8[k] = f_0 * mp_14[k];

        t_9[k] = -2.0 * kp_3[k]
                 + f_0 * mp_18[k];

        t_10[k] = -2.0 * kp_4[k]
                  + f_0 * mp_19[k];

        t_11[k] = -2.0 * kp_5[k]
                  + f_0 * mp_20[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, kp_6, kp_7, kp_8, mp_21, mp_22, \
                         mp_23, mp_24, mp_25, mp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -kp_6[k]
                  + f_0 * mp_21[k];

        t_13[k] = -kp_7[k]
                  + f_0 * mp_22[k];

        t_14[k] = -kp_8[k]
                  + f_0 * mp_23[k];

        t_15[k] = f_0 * mp_24[k];

        t_16[k] = f_0 * mp_25[k];

        t_17[k] = f_0 * mp_26[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, kp_9, kp_10, kp_11, kp_12, kp_13, \
                         mp_30, mp_31, mp_32, mp_33, mp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -3.0 * kp_9[k]
                  + f_0 * mp_30[k];

        t_19[k] = -3.0 * kp_10[k]
                  + f_0 * mp_31[k];

        t_20[k] = -3.0 * kp_11[k]
                  + f_0 * mp_32[k];

        t_21[k] = -2.0 * kp_12[k]
                  + f_0 * mp_33[k];

        t_22[k] = -2.0 * kp_13[k]
                  + f_0 * mp_34[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, t_28, kp_14, kp_15, kp_16, kp_17, \
                         mp_35, mp_36, mp_37, mp_38, mp_39, mp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -2.0 * kp_14[k]
                  + f_0 * mp_35[k];

        t_24[k] = -kp_15[k]
                  + f_0 * mp_36[k];

        t_25[k] = -kp_16[k]
                  + f_0 * mp_37[k];

        t_26[k] = -kp_17[k]
                  + f_0 * mp_38[k];

        t_27[k] = f_0 * mp_39[k];

        t_28[k] = f_0 * mp_40[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, kp_18, kp_19, kp_20, kp_21, mp_41, \
                         mp_45, mp_46, mp_47, mp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * mp_41[k];

        t_30[k] = -4.0 * kp_18[k]
                  + f_0 * mp_45[k];

        t_31[k] = -4.0 * kp_19[k]
                  + f_0 * mp_46[k];

        t_32[k] = -4.0 * kp_20[k]
                  + f_0 * mp_47[k];

        t_33[k] = -3.0 * kp_21[k]
                  + f_0 * mp_48[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, kp_22, kp_23, kp_24, kp_25, kp_26, \
                         mp_49, mp_50, mp_51, mp_52, mp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = -3.0 * kp_22[k]
                  + f_0 * mp_49[k];

        t_35[k] = -3.0 * kp_23[k]
                  + f_0 * mp_50[k];

        t_36[k] = -2.0 * kp_24[k]
                  + f_0 * mp_51[k];

        t_37[k] = -2.0 * kp_25[k]
                  + f_0 * mp_52[k];

        t_38[k] = -2.0 * kp_26[k]
                  + f_0 * mp_53[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, t_44, kp_27, kp_28, kp_29, mp_54, \
                         mp_55, mp_56, mp_57, mp_58, mp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -kp_27[k]
                  + f_0 * mp_54[k];

        t_40[k] = -kp_28[k]
                  + f_0 * mp_55[k];

        t_41[k] = -kp_29[k]
                  + f_0 * mp_56[k];

        t_42[k] = f_0 * mp_57[k];

        t_43[k] = f_0 * mp_58[k];

        t_44[k] = f_0 * mp_59[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, kp_30, kp_31, kp_32, kp_33, kp_34, \
                         mp_63, mp_64, mp_65, mp_66, mp_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -5.0 * kp_30[k]
                  + f_0 * mp_63[k];

        t_46[k] = -5.0 * kp_31[k]
                  + f_0 * mp_64[k];

        t_47[k] = -5.0 * kp_32[k]
                  + f_0 * mp_65[k];

        t_48[k] = -4.0 * kp_33[k]
                  + f_0 * mp_66[k];

        t_49[k] = -4.0 * kp_34[k]
                  + f_0 * mp_67[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, kp_35, kp_36, kp_37, kp_38, kp_39, \
                         mp_68, mp_69, mp_70, mp_71, mp_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -4.0 * kp_35[k]
                  + f_0 * mp_68[k];

        t_51[k] = -3.0 * kp_36[k]
                  + f_0 * mp_69[k];

        t_52[k] = -3.0 * kp_37[k]
                  + f_0 * mp_70[k];

        t_53[k] = -3.0 * kp_38[k]
                  + f_0 * mp_71[k];

        t_54[k] = -2.0 * kp_39[k]
                  + f_0 * mp_72[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, kp_40, kp_41, kp_42, kp_43, kp_44, \
                         mp_73, mp_74, mp_75, mp_76, mp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * kp_40[k]
                  + f_0 * mp_73[k];

        t_56[k] = -2.0 * kp_41[k]
                  + f_0 * mp_74[k];

        t_57[k] = -kp_42[k]
                  + f_0 * mp_75[k];

        t_58[k] = -kp_43[k]
                  + f_0 * mp_76[k];

        t_59[k] = -kp_44[k]
                  + f_0 * mp_77[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, kp_45, kp_46, kp_47, mp_78, \
                         mp_79, mp_80, mp_84, mp_85, mp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * mp_78[k];

        t_61[k] = f_0 * mp_79[k];

        t_62[k] = f_0 * mp_80[k];

        t_63[k] = -6.0 * kp_45[k]
                  + f_0 * mp_84[k];

        t_64[k] = -6.0 * kp_46[k]
                  + f_0 * mp_85[k];

        t_65[k] = -6.0 * kp_47[k]
                  + f_0 * mp_86[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, kp_48, kp_49, kp_50, kp_51, kp_52, \
                         mp_87, mp_88, mp_89, mp_90, mp_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -5.0 * kp_48[k]
                  + f_0 * mp_87[k];

        t_67[k] = -5.0 * kp_49[k]
                  + f_0 * mp_88[k];

        t_68[k] = -5.0 * kp_50[k]
                  + f_0 * mp_89[k];

        t_69[k] = -4.0 * kp_51[k]
                  + f_0 * mp_90[k];

        t_70[k] = -4.0 * kp_52[k]
                  + f_0 * mp_91[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, kp_53, kp_54, kp_55, kp_56, kp_57, \
                         mp_92, mp_93, mp_94, mp_95, mp_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -4.0 * kp_53[k]
                  + f_0 * mp_92[k];

        t_72[k] = -3.0 * kp_54[k]
                  + f_0 * mp_93[k];

        t_73[k] = -3.0 * kp_55[k]
                  + f_0 * mp_94[k];

        t_74[k] = -3.0 * kp_56[k]
                  + f_0 * mp_95[k];

        t_75[k] = -2.0 * kp_57[k]
                  + f_0 * mp_96[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, kp_58, kp_59, kp_60, kp_61, kp_62, \
                         mp_97, mp_98, mp_99, mp_100, mp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -2.0 * kp_58[k]
                  + f_0 * mp_97[k];

        t_77[k] = -2.0 * kp_59[k]
                  + f_0 * mp_98[k];

        t_78[k] = -kp_60[k]
                  + f_0 * mp_99[k];

        t_79[k] = -kp_61[k]
                  + f_0 * mp_100[k];

        t_80[k] = -kp_62[k]
                  + f_0 * mp_101[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, kp_63, kp_64, kp_65, mp_102, \
                         mp_103, mp_104, mp_108, mp_109, mp_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_0 * mp_102[k];

        t_82[k] = f_0 * mp_103[k];

        t_83[k] = f_0 * mp_104[k];

        t_84[k] = -7.0 * kp_63[k]
                  + f_0 * mp_108[k];

        t_85[k] = -7.0 * kp_64[k]
                  + f_0 * mp_109[k];

        t_86[k] = -7.0 * kp_65[k]
                  + f_0 * mp_110[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, kp_66, kp_67, kp_68, kp_69, kp_70, \
                         mp_111, mp_112, mp_113, mp_114, mp_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -6.0 * kp_66[k]
                  + f_0 * mp_111[k];

        t_88[k] = -6.0 * kp_67[k]
                  + f_0 * mp_112[k];

        t_89[k] = -6.0 * kp_68[k]
                  + f_0 * mp_113[k];

        t_90[k] = -5.0 * kp_69[k]
                  + f_0 * mp_114[k];

        t_91[k] = -5.0 * kp_70[k]
                  + f_0 * mp_115[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, kp_71, kp_72, kp_73, kp_74, kp_75, \
                         mp_116, mp_117, mp_118, mp_119, mp_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = -5.0 * kp_71[k]
                  + f_0 * mp_116[k];

        t_93[k] = -4.0 * kp_72[k]
                  + f_0 * mp_117[k];

        t_94[k] = -4.0 * kp_73[k]
                  + f_0 * mp_118[k];

        t_95[k] = -4.0 * kp_74[k]
                  + f_0 * mp_119[k];

        t_96[k] = -3.0 * kp_75[k]
                  + f_0 * mp_120[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, kp_76, kp_77, kp_78, kp_79, kp_80, \
                         mp_121, mp_122, mp_123, mp_124, mp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -3.0 * kp_76[k]
                  + f_0 * mp_121[k];

        t_98[k] = -3.0 * kp_77[k]
                  + f_0 * mp_122[k];

        t_99[k] = -2.0 * kp_78[k]
                  + f_0 * mp_123[k];

        t_100[k] = -2.0 * kp_79[k]
                   + f_0 * mp_124[k];

        t_101[k] = -2.0 * kp_80[k]
                   + f_0 * mp_125[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, kp_81, kp_82, kp_83, \
                         mp_126, mp_127, mp_128, mp_129, mp_130, \
                         mp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -kp_81[k]
                   + f_0 * mp_126[k];

        t_103[k] = -kp_82[k]
                   + f_0 * mp_127[k];

        t_104[k] = -kp_83[k]
                   + f_0 * mp_128[k];

        t_105[k] = f_0 * mp_129[k];

        t_106[k] = f_0 * mp_130[k];

        t_107[k] = f_0 * mp_131[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, kp_84, kp_85, kp_86, kp_87, kp_88, \
                         mp_135, mp_136, mp_137, mp_138, mp_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = -8.0 * kp_84[k]
                   + f_0 * mp_135[k];

        t_109[k] = -8.0 * kp_85[k]
                   + f_0 * mp_136[k];

        t_110[k] = -8.0 * kp_86[k]
                   + f_0 * mp_137[k];

        t_111[k] = -7.0 * kp_87[k]
                   + f_0 * mp_138[k];

        t_112[k] = -7.0 * kp_88[k]
                   + f_0 * mp_139[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, kp_89, kp_90, kp_91, kp_92, kp_93, \
                         mp_140, mp_141, mp_142, mp_143, mp_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = -7.0 * kp_89[k]
                   + f_0 * mp_140[k];

        t_114[k] = -6.0 * kp_90[k]
                   + f_0 * mp_141[k];

        t_115[k] = -6.0 * kp_91[k]
                   + f_0 * mp_142[k];

        t_116[k] = -6.0 * kp_92[k]
                   + f_0 * mp_143[k];

        t_117[k] = -5.0 * kp_93[k]
                   + f_0 * mp_144[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, kp_94, kp_95, kp_96, kp_97, kp_98, \
                         mp_145, mp_146, mp_147, mp_148, mp_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = -5.0 * kp_94[k]
                   + f_0 * mp_145[k];

        t_119[k] = -5.0 * kp_95[k]
                   + f_0 * mp_146[k];

        t_120[k] = -4.0 * kp_96[k]
                   + f_0 * mp_147[k];

        t_121[k] = -4.0 * kp_97[k]
                   + f_0 * mp_148[k];

        t_122[k] = -4.0 * kp_98[k]
                   + f_0 * mp_149[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, t_127, kp_99, kp_100, kp_101, kp_102, \
                         kp_103, mp_150, mp_151, mp_152, mp_153, \
                         mp_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = -3.0 * kp_99[k]
                   + f_0 * mp_150[k];

        t_124[k] = -3.0 * kp_100[k]
                   + f_0 * mp_151[k];

        t_125[k] = -3.0 * kp_101[k]
                   + f_0 * mp_152[k];

        t_126[k] = -2.0 * kp_102[k]
                   + f_0 * mp_153[k];

        t_127[k] = -2.0 * kp_103[k]
                   + f_0 * mp_154[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, kp_104, kp_105, kp_106, \
                         kp_107, mp_155, mp_156, mp_157, mp_158, mp_159, \
                         mp_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = -2.0 * kp_104[k]
                   + f_0 * mp_155[k];

        t_129[k] = -kp_105[k]
                   + f_0 * mp_156[k];

        t_130[k] = -kp_106[k]
                   + f_0 * mp_157[k];

        t_131[k] = -kp_107[k]
                   + f_0 * mp_158[k];

        t_132[k] = f_0 * mp_159[k];

        t_133[k] = f_0 * mp_160[k];
    }

#pragma omp simd aligned(t_134, mp_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_0 * mp_161[k];
    }
}

auto
compute_prim_geom_10_lp_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t kp, const size_t mp,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);
    auto *t_9 = buffer.data(target + 9);
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);
    auto *t_45 = buffer.data(target + 45);
    auto *t_46 = buffer.data(target + 46);
    auto *t_47 = buffer.data(target + 47);
    auto *t_48 = buffer.data(target + 48);
    auto *t_49 = buffer.data(target + 49);
    auto *t_50 = buffer.data(target + 50);
    auto *t_51 = buffer.data(target + 51);
    auto *t_52 = buffer.data(target + 52);
    auto *t_53 = buffer.data(target + 53);
    auto *t_54 = buffer.data(target + 54);
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);
    auto *t_100 = buffer.data(target + 100);
    auto *t_101 = buffer.data(target + 101);
    auto *t_102 = buffer.data(target + 102);
    auto *t_103 = buffer.data(target + 103);
    auto *t_104 = buffer.data(target + 104);
    auto *t_105 = buffer.data(target + 105);
    auto *t_106 = buffer.data(target + 106);
    auto *t_107 = buffer.data(target + 107);
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);
    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);
    const auto *kp_49 = buffer.data(kp + 49);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_51 = buffer.data(kp + 51);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_54 = buffer.data(kp + 54);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_57 = buffer.data(kp + 57);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_60 = buffer.data(kp + 60);
    const auto *kp_61 = buffer.data(kp + 61);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_63 = buffer.data(kp + 63);
    const auto *kp_64 = buffer.data(kp + 64);
    const auto *kp_65 = buffer.data(kp + 65);
    const auto *kp_66 = buffer.data(kp + 66);
    const auto *kp_67 = buffer.data(kp + 67);
    const auto *kp_68 = buffer.data(kp + 68);
    const auto *kp_69 = buffer.data(kp + 69);
    const auto *kp_70 = buffer.data(kp + 70);
    const auto *kp_71 = buffer.data(kp + 71);
    const auto *kp_72 = buffer.data(kp + 72);
    const auto *kp_73 = buffer.data(kp + 73);
    const auto *kp_74 = buffer.data(kp + 74);
    const auto *kp_75 = buffer.data(kp + 75);
    const auto *kp_76 = buffer.data(kp + 76);
    const auto *kp_77 = buffer.data(kp + 77);
    const auto *kp_78 = buffer.data(kp + 78);
    const auto *kp_79 = buffer.data(kp + 79);
    const auto *kp_80 = buffer.data(kp + 80);
    const auto *kp_81 = buffer.data(kp + 81);
    const auto *kp_82 = buffer.data(kp + 82);
    const auto *kp_83 = buffer.data(kp + 83);
    const auto *kp_84 = buffer.data(kp + 84);
    const auto *kp_85 = buffer.data(kp + 85);
    const auto *kp_86 = buffer.data(kp + 86);
    const auto *kp_87 = buffer.data(kp + 87);
    const auto *kp_88 = buffer.data(kp + 88);
    const auto *kp_89 = buffer.data(kp + 89);
    const auto *kp_90 = buffer.data(kp + 90);
    const auto *kp_91 = buffer.data(kp + 91);
    const auto *kp_92 = buffer.data(kp + 92);
    const auto *kp_93 = buffer.data(kp + 93);
    const auto *kp_94 = buffer.data(kp + 94);
    const auto *kp_95 = buffer.data(kp + 95);
    const auto *kp_96 = buffer.data(kp + 96);
    const auto *kp_97 = buffer.data(kp + 97);
    const auto *kp_98 = buffer.data(kp + 98);
    const auto *kp_99 = buffer.data(kp + 99);
    const auto *kp_100 = buffer.data(kp + 100);
    const auto *kp_101 = buffer.data(kp + 101);
    const auto *kp_102 = buffer.data(kp + 102);
    const auto *kp_103 = buffer.data(kp + 103);
    const auto *kp_104 = buffer.data(kp + 104);
    const auto *kp_105 = buffer.data(kp + 105);
    const auto *kp_106 = buffer.data(kp + 106);
    const auto *kp_107 = buffer.data(kp + 107);

    const auto *mp_6 = buffer.data(mp + 6);
    const auto *mp_7 = buffer.data(mp + 7);
    const auto *mp_8 = buffer.data(mp + 8);
    const auto *mp_12 = buffer.data(mp + 12);
    const auto *mp_13 = buffer.data(mp + 13);
    const auto *mp_14 = buffer.data(mp + 14);
    const auto *mp_15 = buffer.data(mp + 15);
    const auto *mp_16 = buffer.data(mp + 16);
    const auto *mp_17 = buffer.data(mp + 17);
    const auto *mp_21 = buffer.data(mp + 21);
    const auto *mp_22 = buffer.data(mp + 22);
    const auto *mp_23 = buffer.data(mp + 23);
    const auto *mp_24 = buffer.data(mp + 24);
    const auto *mp_25 = buffer.data(mp + 25);
    const auto *mp_26 = buffer.data(mp + 26);
    const auto *mp_27 = buffer.data(mp + 27);
    const auto *mp_28 = buffer.data(mp + 28);
    const auto *mp_29 = buffer.data(mp + 29);
    const auto *mp_33 = buffer.data(mp + 33);
    const auto *mp_34 = buffer.data(mp + 34);
    const auto *mp_35 = buffer.data(mp + 35);
    const auto *mp_36 = buffer.data(mp + 36);
    const auto *mp_37 = buffer.data(mp + 37);
    const auto *mp_38 = buffer.data(mp + 38);
    const auto *mp_39 = buffer.data(mp + 39);
    const auto *mp_40 = buffer.data(mp + 40);
    const auto *mp_41 = buffer.data(mp + 41);
    const auto *mp_42 = buffer.data(mp + 42);
    const auto *mp_43 = buffer.data(mp + 43);
    const auto *mp_44 = buffer.data(mp + 44);
    const auto *mp_48 = buffer.data(mp + 48);
    const auto *mp_49 = buffer.data(mp + 49);
    const auto *mp_50 = buffer.data(mp + 50);
    const auto *mp_51 = buffer.data(mp + 51);
    const auto *mp_52 = buffer.data(mp + 52);
    const auto *mp_53 = buffer.data(mp + 53);
    const auto *mp_54 = buffer.data(mp + 54);
    const auto *mp_55 = buffer.data(mp + 55);
    const auto *mp_56 = buffer.data(mp + 56);
    const auto *mp_57 = buffer.data(mp + 57);
    const auto *mp_58 = buffer.data(mp + 58);
    const auto *mp_59 = buffer.data(mp + 59);
    const auto *mp_60 = buffer.data(mp + 60);
    const auto *mp_61 = buffer.data(mp + 61);
    const auto *mp_62 = buffer.data(mp + 62);
    const auto *mp_66 = buffer.data(mp + 66);
    const auto *mp_67 = buffer.data(mp + 67);
    const auto *mp_68 = buffer.data(mp + 68);
    const auto *mp_69 = buffer.data(mp + 69);
    const auto *mp_70 = buffer.data(mp + 70);
    const auto *mp_71 = buffer.data(mp + 71);
    const auto *mp_72 = buffer.data(mp + 72);
    const auto *mp_73 = buffer.data(mp + 73);
    const auto *mp_74 = buffer.data(mp + 74);
    const auto *mp_75 = buffer.data(mp + 75);
    const auto *mp_76 = buffer.data(mp + 76);
    const auto *mp_77 = buffer.data(mp + 77);
    const auto *mp_78 = buffer.data(mp + 78);
    const auto *mp_79 = buffer.data(mp + 79);
    const auto *mp_80 = buffer.data(mp + 80);
    const auto *mp_81 = buffer.data(mp + 81);
    const auto *mp_82 = buffer.data(mp + 82);
    const auto *mp_83 = buffer.data(mp + 83);
    const auto *mp_87 = buffer.data(mp + 87);
    const auto *mp_88 = buffer.data(mp + 88);
    const auto *mp_89 = buffer.data(mp + 89);
    const auto *mp_90 = buffer.data(mp + 90);
    const auto *mp_91 = buffer.data(mp + 91);
    const auto *mp_92 = buffer.data(mp + 92);
    const auto *mp_93 = buffer.data(mp + 93);
    const auto *mp_94 = buffer.data(mp + 94);
    const auto *mp_95 = buffer.data(mp + 95);
    const auto *mp_96 = buffer.data(mp + 96);
    const auto *mp_97 = buffer.data(mp + 97);
    const auto *mp_98 = buffer.data(mp + 98);
    const auto *mp_99 = buffer.data(mp + 99);
    const auto *mp_100 = buffer.data(mp + 100);
    const auto *mp_101 = buffer.data(mp + 101);
    const auto *mp_102 = buffer.data(mp + 102);
    const auto *mp_103 = buffer.data(mp + 103);
    const auto *mp_104 = buffer.data(mp + 104);
    const auto *mp_105 = buffer.data(mp + 105);
    const auto *mp_106 = buffer.data(mp + 106);
    const auto *mp_107 = buffer.data(mp + 107);
    const auto *mp_111 = buffer.data(mp + 111);
    const auto *mp_112 = buffer.data(mp + 112);
    const auto *mp_113 = buffer.data(mp + 113);
    const auto *mp_114 = buffer.data(mp + 114);
    const auto *mp_115 = buffer.data(mp + 115);
    const auto *mp_116 = buffer.data(mp + 116);
    const auto *mp_117 = buffer.data(mp + 117);
    const auto *mp_118 = buffer.data(mp + 118);
    const auto *mp_119 = buffer.data(mp + 119);
    const auto *mp_120 = buffer.data(mp + 120);
    const auto *mp_121 = buffer.data(mp + 121);
    const auto *mp_122 = buffer.data(mp + 122);
    const auto *mp_123 = buffer.data(mp + 123);
    const auto *mp_124 = buffer.data(mp + 124);
    const auto *mp_125 = buffer.data(mp + 125);
    const auto *mp_126 = buffer.data(mp + 126);
    const auto *mp_127 = buffer.data(mp + 127);
    const auto *mp_128 = buffer.data(mp + 128);
    const auto *mp_129 = buffer.data(mp + 129);
    const auto *mp_130 = buffer.data(mp + 130);
    const auto *mp_131 = buffer.data(mp + 131);
    const auto *mp_132 = buffer.data(mp + 132);
    const auto *mp_133 = buffer.data(mp + 133);
    const auto *mp_134 = buffer.data(mp + 134);
    const auto *mp_138 = buffer.data(mp + 138);
    const auto *mp_139 = buffer.data(mp + 139);
    const auto *mp_140 = buffer.data(mp + 140);
    const auto *mp_141 = buffer.data(mp + 141);
    const auto *mp_142 = buffer.data(mp + 142);
    const auto *mp_143 = buffer.data(mp + 143);
    const auto *mp_144 = buffer.data(mp + 144);
    const auto *mp_145 = buffer.data(mp + 145);
    const auto *mp_146 = buffer.data(mp + 146);
    const auto *mp_147 = buffer.data(mp + 147);
    const auto *mp_148 = buffer.data(mp + 148);
    const auto *mp_149 = buffer.data(mp + 149);
    const auto *mp_150 = buffer.data(mp + 150);
    const auto *mp_151 = buffer.data(mp + 151);
    const auto *mp_152 = buffer.data(mp + 152);
    const auto *mp_153 = buffer.data(mp + 153);
    const auto *mp_154 = buffer.data(mp + 154);
    const auto *mp_155 = buffer.data(mp + 155);
    const auto *mp_156 = buffer.data(mp + 156);
    const auto *mp_157 = buffer.data(mp + 157);
    const auto *mp_158 = buffer.data(mp + 158);
    const auto *mp_159 = buffer.data(mp + 159);
    const auto *mp_160 = buffer.data(mp + 160);
    const auto *mp_161 = buffer.data(mp + 161);
    const auto *mp_162 = buffer.data(mp + 162);
    const auto *mp_163 = buffer.data(mp + 163);
    const auto *mp_164 = buffer.data(mp + 164);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, kp_0, mp_6, mp_7, mp_8, mp_12, \
                         mp_13, mp_14, mp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * mp_6[k];

        t_1[k] = f_0 * mp_7[k];

        t_2[k] = f_0 * mp_8[k];

        t_3[k] = f_0 * mp_12[k];

        t_4[k] = f_0 * mp_13[k];

        t_5[k] = f_0 * mp_14[k];

        t_6[k] = -kp_0[k]
                 + f_0 * mp_15[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, kp_1, kp_2, kp_3, mp_16, mp_17, \
                         mp_21, mp_22, mp_23, mp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -kp_1[k]
                 + f_0 * mp_16[k];

        t_8[k] = -kp_2[k]
                 + f_0 * mp_17[k];

        t_9[k] = f_0 * mp_21[k];

        t_10[k] = f_0 * mp_22[k];

        t_11[k] = f_0 * mp_23[k];

        t_12[k] = -kp_3[k]
                  + f_0 * mp_24[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, kp_4, kp_5, kp_6, kp_7, kp_8, mp_25, \
                         mp_26, mp_27, mp_28, mp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -kp_4[k]
                  + f_0 * mp_25[k];

        t_14[k] = -kp_5[k]
                  + f_0 * mp_26[k];

        t_15[k] = -2.0 * kp_6[k]
                  + f_0 * mp_27[k];

        t_16[k] = -2.0 * kp_7[k]
                  + f_0 * mp_28[k];

        t_17[k] = -2.0 * kp_8[k]
                  + f_0 * mp_29[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, kp_9, kp_10, kp_11, mp_33, mp_34, \
                         mp_35, mp_36, mp_37, mp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * mp_33[k];

        t_19[k] = f_0 * mp_34[k];

        t_20[k] = f_0 * mp_35[k];

        t_21[k] = -kp_9[k]
                  + f_0 * mp_36[k];

        t_22[k] = -kp_10[k]
                  + f_0 * mp_37[k];

        t_23[k] = -kp_11[k]
                  + f_0 * mp_38[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, kp_12, kp_13, kp_14, kp_15, kp_16, \
                         mp_39, mp_40, mp_41, mp_42, mp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -2.0 * kp_12[k]
                  + f_0 * mp_39[k];

        t_25[k] = -2.0 * kp_13[k]
                  + f_0 * mp_40[k];

        t_26[k] = -2.0 * kp_14[k]
                  + f_0 * mp_41[k];

        t_27[k] = -3.0 * kp_15[k]
                  + f_0 * mp_42[k];

        t_28[k] = -3.0 * kp_16[k]
                  + f_0 * mp_43[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, kp_17, kp_18, kp_19, mp_44, \
                         mp_48, mp_49, mp_50, mp_51, mp_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -3.0 * kp_17[k]
                  + f_0 * mp_44[k];

        t_30[k] = f_0 * mp_48[k];

        t_31[k] = f_0 * mp_49[k];

        t_32[k] = f_0 * mp_50[k];

        t_33[k] = -kp_18[k]
                  + f_0 * mp_51[k];

        t_34[k] = -kp_19[k]
                  + f_0 * mp_52[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, kp_20, kp_21, kp_22, kp_23, kp_24, \
                         mp_53, mp_54, mp_55, mp_56, mp_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -kp_20[k]
                  + f_0 * mp_53[k];

        t_36[k] = -2.0 * kp_21[k]
                  + f_0 * mp_54[k];

        t_37[k] = -2.0 * kp_22[k]
                  + f_0 * mp_55[k];

        t_38[k] = -2.0 * kp_23[k]
                  + f_0 * mp_56[k];

        t_39[k] = -3.0 * kp_24[k]
                  + f_0 * mp_57[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, kp_25, kp_26, kp_27, kp_28, kp_29, \
                         mp_58, mp_59, mp_60, mp_61, mp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * kp_25[k]
                  + f_0 * mp_58[k];

        t_41[k] = -3.0 * kp_26[k]
                  + f_0 * mp_59[k];

        t_42[k] = -4.0 * kp_27[k]
                  + f_0 * mp_60[k];

        t_43[k] = -4.0 * kp_28[k]
                  + f_0 * mp_61[k];

        t_44[k] = -4.0 * kp_29[k]
                  + f_0 * mp_62[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, kp_30, kp_31, kp_32, mp_66, \
                         mp_67, mp_68, mp_69, mp_70, mp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * mp_66[k];

        t_46[k] = f_0 * mp_67[k];

        t_47[k] = f_0 * mp_68[k];

        t_48[k] = -kp_30[k]
                  + f_0 * mp_69[k];

        t_49[k] = -kp_31[k]
                  + f_0 * mp_70[k];

        t_50[k] = -kp_32[k]
                  + f_0 * mp_71[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, kp_33, kp_34, kp_35, kp_36, kp_37, \
                         mp_72, mp_73, mp_74, mp_75, mp_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -2.0 * kp_33[k]
                  + f_0 * mp_72[k];

        t_52[k] = -2.0 * kp_34[k]
                  + f_0 * mp_73[k];

        t_53[k] = -2.0 * kp_35[k]
                  + f_0 * mp_74[k];

        t_54[k] = -3.0 * kp_36[k]
                  + f_0 * mp_75[k];

        t_55[k] = -3.0 * kp_37[k]
                  + f_0 * mp_76[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, kp_38, kp_39, kp_40, kp_41, kp_42, \
                         mp_77, mp_78, mp_79, mp_80, mp_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -3.0 * kp_38[k]
                  + f_0 * mp_77[k];

        t_57[k] = -4.0 * kp_39[k]
                  + f_0 * mp_78[k];

        t_58[k] = -4.0 * kp_40[k]
                  + f_0 * mp_79[k];

        t_59[k] = -4.0 * kp_41[k]
                  + f_0 * mp_80[k];

        t_60[k] = -5.0 * kp_42[k]
                  + f_0 * mp_81[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, kp_43, kp_44, kp_45, mp_82, \
                         mp_83, mp_87, mp_88, mp_89, mp_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -5.0 * kp_43[k]
                  + f_0 * mp_82[k];

        t_62[k] = -5.0 * kp_44[k]
                  + f_0 * mp_83[k];

        t_63[k] = f_0 * mp_87[k];

        t_64[k] = f_0 * mp_88[k];

        t_65[k] = f_0 * mp_89[k];

        t_66[k] = -kp_45[k]
                  + f_0 * mp_90[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, kp_46, kp_47, kp_48, kp_49, kp_50, \
                         mp_91, mp_92, mp_93, mp_94, mp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -kp_46[k]
                  + f_0 * mp_91[k];

        t_68[k] = -kp_47[k]
                  + f_0 * mp_92[k];

        t_69[k] = -2.0 * kp_48[k]
                  + f_0 * mp_93[k];

        t_70[k] = -2.0 * kp_49[k]
                  + f_0 * mp_94[k];

        t_71[k] = -2.0 * kp_50[k]
                  + f_0 * mp_95[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, kp_51, kp_52, kp_53, kp_54, kp_55, \
                         mp_96, mp_97, mp_98, mp_99, mp_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -3.0 * kp_51[k]
                  + f_0 * mp_96[k];

        t_73[k] = -3.0 * kp_52[k]
                  + f_0 * mp_97[k];

        t_74[k] = -3.0 * kp_53[k]
                  + f_0 * mp_98[k];

        t_75[k] = -4.0 * kp_54[k]
                  + f_0 * mp_99[k];

        t_76[k] = -4.0 * kp_55[k]
                  + f_0 * mp_100[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, kp_56, kp_57, kp_58, kp_59, kp_60, \
                         mp_101, mp_102, mp_103, mp_104, mp_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -4.0 * kp_56[k]
                  + f_0 * mp_101[k];

        t_78[k] = -5.0 * kp_57[k]
                  + f_0 * mp_102[k];

        t_79[k] = -5.0 * kp_58[k]
                  + f_0 * mp_103[k];

        t_80[k] = -5.0 * kp_59[k]
                  + f_0 * mp_104[k];

        t_81[k] = -6.0 * kp_60[k]
                  + f_0 * mp_105[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, t_87, kp_61, kp_62, kp_63, mp_106, \
                         mp_107, mp_111, mp_112, mp_113, mp_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -6.0 * kp_61[k]
                  + f_0 * mp_106[k];

        t_83[k] = -6.0 * kp_62[k]
                  + f_0 * mp_107[k];

        t_84[k] = f_0 * mp_111[k];

        t_85[k] = f_0 * mp_112[k];

        t_86[k] = f_0 * mp_113[k];

        t_87[k] = -kp_63[k]
                  + f_0 * mp_114[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, kp_64, kp_65, kp_66, kp_67, kp_68, \
                         mp_115, mp_116, mp_117, mp_118, mp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = -kp_64[k]
                  + f_0 * mp_115[k];

        t_89[k] = -kp_65[k]
                  + f_0 * mp_116[k];

        t_90[k] = -2.0 * kp_66[k]
                  + f_0 * mp_117[k];

        t_91[k] = -2.0 * kp_67[k]
                  + f_0 * mp_118[k];

        t_92[k] = -2.0 * kp_68[k]
                  + f_0 * mp_119[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, kp_69, kp_70, kp_71, kp_72, kp_73, \
                         mp_120, mp_121, mp_122, mp_123, mp_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = -3.0 * kp_69[k]
                  + f_0 * mp_120[k];

        t_94[k] = -3.0 * kp_70[k]
                  + f_0 * mp_121[k];

        t_95[k] = -3.0 * kp_71[k]
                  + f_0 * mp_122[k];

        t_96[k] = -4.0 * kp_72[k]
                  + f_0 * mp_123[k];

        t_97[k] = -4.0 * kp_73[k]
                  + f_0 * mp_124[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, kp_74, kp_75, kp_76, kp_77, kp_78, \
                         mp_125, mp_126, mp_127, mp_128, mp_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = -4.0 * kp_74[k]
                  + f_0 * mp_125[k];

        t_99[k] = -5.0 * kp_75[k]
                  + f_0 * mp_126[k];

        t_100[k] = -5.0 * kp_76[k]
                   + f_0 * mp_127[k];

        t_101[k] = -5.0 * kp_77[k]
                   + f_0 * mp_128[k];

        t_102[k] = -6.0 * kp_78[k]
                   + f_0 * mp_129[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, kp_79, kp_80, kp_81, kp_82, kp_83, \
                         mp_130, mp_131, mp_132, mp_133, mp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = -6.0 * kp_79[k]
                   + f_0 * mp_130[k];

        t_104[k] = -6.0 * kp_80[k]
                   + f_0 * mp_131[k];

        t_105[k] = -7.0 * kp_81[k]
                   + f_0 * mp_132[k];

        t_106[k] = -7.0 * kp_82[k]
                   + f_0 * mp_133[k];

        t_107[k] = -7.0 * kp_83[k]
                   + f_0 * mp_134[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, t_113, kp_84, kp_85, kp_86, \
                         mp_138, mp_139, mp_140, mp_141, mp_142, \
                         mp_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_0 * mp_138[k];

        t_109[k] = f_0 * mp_139[k];

        t_110[k] = f_0 * mp_140[k];

        t_111[k] = -kp_84[k]
                   + f_0 * mp_141[k];

        t_112[k] = -kp_85[k]
                   + f_0 * mp_142[k];

        t_113[k] = -kp_86[k]
                   + f_0 * mp_143[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, kp_87, kp_88, kp_89, kp_90, kp_91, \
                         mp_144, mp_145, mp_146, mp_147, mp_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -2.0 * kp_87[k]
                   + f_0 * mp_144[k];

        t_115[k] = -2.0 * kp_88[k]
                   + f_0 * mp_145[k];

        t_116[k] = -2.0 * kp_89[k]
                   + f_0 * mp_146[k];

        t_117[k] = -3.0 * kp_90[k]
                   + f_0 * mp_147[k];

        t_118[k] = -3.0 * kp_91[k]
                   + f_0 * mp_148[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, kp_92, kp_93, kp_94, kp_95, kp_96, \
                         mp_149, mp_150, mp_151, mp_152, mp_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -3.0 * kp_92[k]
                   + f_0 * mp_149[k];

        t_120[k] = -4.0 * kp_93[k]
                   + f_0 * mp_150[k];

        t_121[k] = -4.0 * kp_94[k]
                   + f_0 * mp_151[k];

        t_122[k] = -4.0 * kp_95[k]
                   + f_0 * mp_152[k];

        t_123[k] = -5.0 * kp_96[k]
                   + f_0 * mp_153[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, kp_97, kp_98, kp_99, kp_100, \
                         kp_101, mp_154, mp_155, mp_156, mp_157, \
                         mp_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -5.0 * kp_97[k]
                   + f_0 * mp_154[k];

        t_125[k] = -5.0 * kp_98[k]
                   + f_0 * mp_155[k];

        t_126[k] = -6.0 * kp_99[k]
                   + f_0 * mp_156[k];

        t_127[k] = -6.0 * kp_100[k]
                   + f_0 * mp_157[k];

        t_128[k] = -6.0 * kp_101[k]
                   + f_0 * mp_158[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, kp_102, kp_103, kp_104, kp_105, \
                         kp_106, mp_159, mp_160, mp_161, mp_162, \
                         mp_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -7.0 * kp_102[k]
                   + f_0 * mp_159[k];

        t_130[k] = -7.0 * kp_103[k]
                   + f_0 * mp_160[k];

        t_131[k] = -7.0 * kp_104[k]
                   + f_0 * mp_161[k];

        t_132[k] = -8.0 * kp_105[k]
                   + f_0 * mp_162[k];

        t_133[k] = -8.0 * kp_106[k]
                   + f_0 * mp_163[k];
    }

#pragma omp simd aligned(t_134, kp_107, mp_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -8.0 * kp_107[k]
                   + f_0 * mp_164[k];
    }
}

}  // namespace simdt2ceri
