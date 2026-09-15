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


#include "SimdElectronRepulsionGeom10VrrRecLH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_lh_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
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
    auto *t_135 = buffer.data(target + 135);
    auto *t_136 = buffer.data(target + 136);
    auto *t_137 = buffer.data(target + 137);
    auto *t_138 = buffer.data(target + 138);
    auto *t_139 = buffer.data(target + 139);
    auto *t_140 = buffer.data(target + 140);
    auto *t_141 = buffer.data(target + 141);
    auto *t_142 = buffer.data(target + 142);
    auto *t_143 = buffer.data(target + 143);
    auto *t_144 = buffer.data(target + 144);
    auto *t_145 = buffer.data(target + 145);
    auto *t_146 = buffer.data(target + 146);
    auto *t_147 = buffer.data(target + 147);
    auto *t_148 = buffer.data(target + 148);
    auto *t_149 = buffer.data(target + 149);

    const auto *kh_0 = buffer.data(kh + 0);
    const auto *kh_1 = buffer.data(kh + 1);
    const auto *kh_2 = buffer.data(kh + 2);
    const auto *kh_3 = buffer.data(kh + 3);
    const auto *kh_4 = buffer.data(kh + 4);
    const auto *kh_5 = buffer.data(kh + 5);
    const auto *kh_6 = buffer.data(kh + 6);
    const auto *kh_7 = buffer.data(kh + 7);
    const auto *kh_8 = buffer.data(kh + 8);
    const auto *kh_9 = buffer.data(kh + 9);
    const auto *kh_10 = buffer.data(kh + 10);
    const auto *kh_11 = buffer.data(kh + 11);
    const auto *kh_12 = buffer.data(kh + 12);
    const auto *kh_13 = buffer.data(kh + 13);
    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_15 = buffer.data(kh + 15);
    const auto *kh_16 = buffer.data(kh + 16);
    const auto *kh_17 = buffer.data(kh + 17);
    const auto *kh_18 = buffer.data(kh + 18);
    const auto *kh_19 = buffer.data(kh + 19);
    const auto *kh_20 = buffer.data(kh + 20);
    const auto *kh_21 = buffer.data(kh + 21);
    const auto *kh_22 = buffer.data(kh + 22);
    const auto *kh_23 = buffer.data(kh + 23);
    const auto *kh_24 = buffer.data(kh + 24);
    const auto *kh_25 = buffer.data(kh + 25);
    const auto *kh_26 = buffer.data(kh + 26);
    const auto *kh_27 = buffer.data(kh + 27);
    const auto *kh_28 = buffer.data(kh + 28);
    const auto *kh_29 = buffer.data(kh + 29);
    const auto *kh_30 = buffer.data(kh + 30);
    const auto *kh_31 = buffer.data(kh + 31);
    const auto *kh_32 = buffer.data(kh + 32);
    const auto *kh_33 = buffer.data(kh + 33);
    const auto *kh_34 = buffer.data(kh + 34);
    const auto *kh_35 = buffer.data(kh + 35);
    const auto *kh_36 = buffer.data(kh + 36);
    const auto *kh_37 = buffer.data(kh + 37);
    const auto *kh_38 = buffer.data(kh + 38);
    const auto *kh_39 = buffer.data(kh + 39);
    const auto *kh_40 = buffer.data(kh + 40);
    const auto *kh_41 = buffer.data(kh + 41);
    const auto *kh_42 = buffer.data(kh + 42);
    const auto *kh_43 = buffer.data(kh + 43);
    const auto *kh_44 = buffer.data(kh + 44);
    const auto *kh_45 = buffer.data(kh + 45);
    const auto *kh_46 = buffer.data(kh + 46);
    const auto *kh_47 = buffer.data(kh + 47);
    const auto *kh_48 = buffer.data(kh + 48);
    const auto *kh_49 = buffer.data(kh + 49);
    const auto *kh_50 = buffer.data(kh + 50);
    const auto *kh_51 = buffer.data(kh + 51);
    const auto *kh_52 = buffer.data(kh + 52);
    const auto *kh_53 = buffer.data(kh + 53);
    const auto *kh_54 = buffer.data(kh + 54);
    const auto *kh_55 = buffer.data(kh + 55);
    const auto *kh_56 = buffer.data(kh + 56);
    const auto *kh_57 = buffer.data(kh + 57);
    const auto *kh_58 = buffer.data(kh + 58);
    const auto *kh_59 = buffer.data(kh + 59);
    const auto *kh_60 = buffer.data(kh + 60);
    const auto *kh_61 = buffer.data(kh + 61);
    const auto *kh_62 = buffer.data(kh + 62);
    const auto *kh_63 = buffer.data(kh + 63);
    const auto *kh_64 = buffer.data(kh + 64);
    const auto *kh_65 = buffer.data(kh + 65);
    const auto *kh_66 = buffer.data(kh + 66);
    const auto *kh_67 = buffer.data(kh + 67);
    const auto *kh_68 = buffer.data(kh + 68);
    const auto *kh_69 = buffer.data(kh + 69);
    const auto *kh_70 = buffer.data(kh + 70);
    const auto *kh_71 = buffer.data(kh + 71);
    const auto *kh_72 = buffer.data(kh + 72);
    const auto *kh_73 = buffer.data(kh + 73);
    const auto *kh_74 = buffer.data(kh + 74);
    const auto *kh_75 = buffer.data(kh + 75);
    const auto *kh_76 = buffer.data(kh + 76);
    const auto *kh_77 = buffer.data(kh + 77);
    const auto *kh_78 = buffer.data(kh + 78);
    const auto *kh_79 = buffer.data(kh + 79);
    const auto *kh_80 = buffer.data(kh + 80);
    const auto *kh_81 = buffer.data(kh + 81);
    const auto *kh_82 = buffer.data(kh + 82);
    const auto *kh_83 = buffer.data(kh + 83);
    const auto *kh_84 = buffer.data(kh + 84);
    const auto *kh_85 = buffer.data(kh + 85);
    const auto *kh_86 = buffer.data(kh + 86);
    const auto *kh_87 = buffer.data(kh + 87);
    const auto *kh_88 = buffer.data(kh + 88);
    const auto *kh_89 = buffer.data(kh + 89);
    const auto *kh_90 = buffer.data(kh + 90);
    const auto *kh_91 = buffer.data(kh + 91);
    const auto *kh_92 = buffer.data(kh + 92);
    const auto *kh_93 = buffer.data(kh + 93);
    const auto *kh_94 = buffer.data(kh + 94);
    const auto *kh_95 = buffer.data(kh + 95);
    const auto *kh_96 = buffer.data(kh + 96);
    const auto *kh_97 = buffer.data(kh + 97);
    const auto *kh_98 = buffer.data(kh + 98);
    const auto *kh_99 = buffer.data(kh + 99);
    const auto *kh_100 = buffer.data(kh + 100);
    const auto *kh_101 = buffer.data(kh + 101);
    const auto *kh_102 = buffer.data(kh + 102);
    const auto *kh_103 = buffer.data(kh + 103);
    const auto *kh_104 = buffer.data(kh + 104);
    const auto *kh_105 = buffer.data(kh + 105);
    const auto *kh_106 = buffer.data(kh + 106);
    const auto *kh_107 = buffer.data(kh + 107);
    const auto *kh_108 = buffer.data(kh + 108);
    const auto *kh_109 = buffer.data(kh + 109);
    const auto *kh_110 = buffer.data(kh + 110);
    const auto *kh_111 = buffer.data(kh + 111);
    const auto *kh_112 = buffer.data(kh + 112);
    const auto *kh_113 = buffer.data(kh + 113);
    const auto *kh_114 = buffer.data(kh + 114);
    const auto *kh_115 = buffer.data(kh + 115);
    const auto *kh_116 = buffer.data(kh + 116);
    const auto *kh_117 = buffer.data(kh + 117);
    const auto *kh_118 = buffer.data(kh + 118);
    const auto *kh_119 = buffer.data(kh + 119);
    const auto *kh_120 = buffer.data(kh + 120);
    const auto *kh_121 = buffer.data(kh + 121);
    const auto *kh_122 = buffer.data(kh + 122);
    const auto *kh_123 = buffer.data(kh + 123);
    const auto *kh_124 = buffer.data(kh + 124);
    const auto *kh_125 = buffer.data(kh + 125);
    const auto *kh_126 = buffer.data(kh + 126);
    const auto *kh_127 = buffer.data(kh + 127);
    const auto *kh_128 = buffer.data(kh + 128);
    const auto *kh_129 = buffer.data(kh + 129);
    const auto *kh_130 = buffer.data(kh + 130);
    const auto *kh_131 = buffer.data(kh + 131);
    const auto *kh_132 = buffer.data(kh + 132);
    const auto *kh_133 = buffer.data(kh + 133);
    const auto *kh_134 = buffer.data(kh + 134);
    const auto *kh_135 = buffer.data(kh + 135);
    const auto *kh_136 = buffer.data(kh + 136);
    const auto *kh_137 = buffer.data(kh + 137);
    const auto *kh_138 = buffer.data(kh + 138);
    const auto *kh_139 = buffer.data(kh + 139);
    const auto *kh_140 = buffer.data(kh + 140);
    const auto *kh_141 = buffer.data(kh + 141);
    const auto *kh_142 = buffer.data(kh + 142);
    const auto *kh_143 = buffer.data(kh + 143);
    const auto *kh_144 = buffer.data(kh + 144);
    const auto *kh_145 = buffer.data(kh + 145);
    const auto *kh_146 = buffer.data(kh + 146);
    const auto *kh_147 = buffer.data(kh + 147);
    const auto *kh_148 = buffer.data(kh + 148);
    const auto *kh_149 = buffer.data(kh + 149);

    const auto *mh_0 = buffer.data(mh + 0);
    const auto *mh_1 = buffer.data(mh + 1);
    const auto *mh_2 = buffer.data(mh + 2);
    const auto *mh_3 = buffer.data(mh + 3);
    const auto *mh_4 = buffer.data(mh + 4);
    const auto *mh_5 = buffer.data(mh + 5);
    const auto *mh_6 = buffer.data(mh + 6);
    const auto *mh_7 = buffer.data(mh + 7);
    const auto *mh_8 = buffer.data(mh + 8);
    const auto *mh_9 = buffer.data(mh + 9);
    const auto *mh_10 = buffer.data(mh + 10);
    const auto *mh_11 = buffer.data(mh + 11);
    const auto *mh_12 = buffer.data(mh + 12);
    const auto *mh_13 = buffer.data(mh + 13);
    const auto *mh_14 = buffer.data(mh + 14);
    const auto *mh_15 = buffer.data(mh + 15);
    const auto *mh_16 = buffer.data(mh + 16);
    const auto *mh_17 = buffer.data(mh + 17);
    const auto *mh_18 = buffer.data(mh + 18);
    const auto *mh_19 = buffer.data(mh + 19);
    const auto *mh_20 = buffer.data(mh + 20);
    const auto *mh_21 = buffer.data(mh + 21);
    const auto *mh_22 = buffer.data(mh + 22);
    const auto *mh_23 = buffer.data(mh + 23);
    const auto *mh_24 = buffer.data(mh + 24);
    const auto *mh_25 = buffer.data(mh + 25);
    const auto *mh_26 = buffer.data(mh + 26);
    const auto *mh_27 = buffer.data(mh + 27);
    const auto *mh_28 = buffer.data(mh + 28);
    const auto *mh_29 = buffer.data(mh + 29);
    const auto *mh_30 = buffer.data(mh + 30);
    const auto *mh_31 = buffer.data(mh + 31);
    const auto *mh_32 = buffer.data(mh + 32);
    const auto *mh_33 = buffer.data(mh + 33);
    const auto *mh_34 = buffer.data(mh + 34);
    const auto *mh_35 = buffer.data(mh + 35);
    const auto *mh_36 = buffer.data(mh + 36);
    const auto *mh_37 = buffer.data(mh + 37);
    const auto *mh_38 = buffer.data(mh + 38);
    const auto *mh_39 = buffer.data(mh + 39);
    const auto *mh_40 = buffer.data(mh + 40);
    const auto *mh_41 = buffer.data(mh + 41);
    const auto *mh_42 = buffer.data(mh + 42);
    const auto *mh_43 = buffer.data(mh + 43);
    const auto *mh_44 = buffer.data(mh + 44);
    const auto *mh_45 = buffer.data(mh + 45);
    const auto *mh_46 = buffer.data(mh + 46);
    const auto *mh_47 = buffer.data(mh + 47);
    const auto *mh_48 = buffer.data(mh + 48);
    const auto *mh_49 = buffer.data(mh + 49);
    const auto *mh_50 = buffer.data(mh + 50);
    const auto *mh_51 = buffer.data(mh + 51);
    const auto *mh_52 = buffer.data(mh + 52);
    const auto *mh_53 = buffer.data(mh + 53);
    const auto *mh_54 = buffer.data(mh + 54);
    const auto *mh_55 = buffer.data(mh + 55);
    const auto *mh_56 = buffer.data(mh + 56);
    const auto *mh_57 = buffer.data(mh + 57);
    const auto *mh_58 = buffer.data(mh + 58);
    const auto *mh_59 = buffer.data(mh + 59);
    const auto *mh_60 = buffer.data(mh + 60);
    const auto *mh_61 = buffer.data(mh + 61);
    const auto *mh_62 = buffer.data(mh + 62);
    const auto *mh_63 = buffer.data(mh + 63);
    const auto *mh_64 = buffer.data(mh + 64);
    const auto *mh_65 = buffer.data(mh + 65);
    const auto *mh_66 = buffer.data(mh + 66);
    const auto *mh_67 = buffer.data(mh + 67);
    const auto *mh_68 = buffer.data(mh + 68);
    const auto *mh_69 = buffer.data(mh + 69);
    const auto *mh_70 = buffer.data(mh + 70);
    const auto *mh_71 = buffer.data(mh + 71);
    const auto *mh_72 = buffer.data(mh + 72);
    const auto *mh_73 = buffer.data(mh + 73);
    const auto *mh_74 = buffer.data(mh + 74);
    const auto *mh_75 = buffer.data(mh + 75);
    const auto *mh_76 = buffer.data(mh + 76);
    const auto *mh_77 = buffer.data(mh + 77);
    const auto *mh_78 = buffer.data(mh + 78);
    const auto *mh_79 = buffer.data(mh + 79);
    const auto *mh_80 = buffer.data(mh + 80);
    const auto *mh_81 = buffer.data(mh + 81);
    const auto *mh_82 = buffer.data(mh + 82);
    const auto *mh_83 = buffer.data(mh + 83);
    const auto *mh_84 = buffer.data(mh + 84);
    const auto *mh_85 = buffer.data(mh + 85);
    const auto *mh_86 = buffer.data(mh + 86);
    const auto *mh_87 = buffer.data(mh + 87);
    const auto *mh_88 = buffer.data(mh + 88);
    const auto *mh_89 = buffer.data(mh + 89);
    const auto *mh_90 = buffer.data(mh + 90);
    const auto *mh_91 = buffer.data(mh + 91);
    const auto *mh_92 = buffer.data(mh + 92);
    const auto *mh_93 = buffer.data(mh + 93);
    const auto *mh_94 = buffer.data(mh + 94);
    const auto *mh_95 = buffer.data(mh + 95);
    const auto *mh_96 = buffer.data(mh + 96);
    const auto *mh_97 = buffer.data(mh + 97);
    const auto *mh_98 = buffer.data(mh + 98);
    const auto *mh_99 = buffer.data(mh + 99);
    const auto *mh_100 = buffer.data(mh + 100);
    const auto *mh_101 = buffer.data(mh + 101);
    const auto *mh_102 = buffer.data(mh + 102);
    const auto *mh_103 = buffer.data(mh + 103);
    const auto *mh_104 = buffer.data(mh + 104);
    const auto *mh_105 = buffer.data(mh + 105);
    const auto *mh_106 = buffer.data(mh + 106);
    const auto *mh_107 = buffer.data(mh + 107);
    const auto *mh_108 = buffer.data(mh + 108);
    const auto *mh_109 = buffer.data(mh + 109);
    const auto *mh_110 = buffer.data(mh + 110);
    const auto *mh_111 = buffer.data(mh + 111);
    const auto *mh_112 = buffer.data(mh + 112);
    const auto *mh_113 = buffer.data(mh + 113);
    const auto *mh_114 = buffer.data(mh + 114);
    const auto *mh_115 = buffer.data(mh + 115);
    const auto *mh_116 = buffer.data(mh + 116);
    const auto *mh_117 = buffer.data(mh + 117);
    const auto *mh_118 = buffer.data(mh + 118);
    const auto *mh_119 = buffer.data(mh + 119);
    const auto *mh_120 = buffer.data(mh + 120);
    const auto *mh_121 = buffer.data(mh + 121);
    const auto *mh_122 = buffer.data(mh + 122);
    const auto *mh_123 = buffer.data(mh + 123);
    const auto *mh_124 = buffer.data(mh + 124);
    const auto *mh_125 = buffer.data(mh + 125);
    const auto *mh_126 = buffer.data(mh + 126);
    const auto *mh_127 = buffer.data(mh + 127);
    const auto *mh_128 = buffer.data(mh + 128);
    const auto *mh_129 = buffer.data(mh + 129);
    const auto *mh_130 = buffer.data(mh + 130);
    const auto *mh_131 = buffer.data(mh + 131);
    const auto *mh_132 = buffer.data(mh + 132);
    const auto *mh_133 = buffer.data(mh + 133);
    const auto *mh_134 = buffer.data(mh + 134);
    const auto *mh_135 = buffer.data(mh + 135);
    const auto *mh_136 = buffer.data(mh + 136);
    const auto *mh_137 = buffer.data(mh + 137);
    const auto *mh_138 = buffer.data(mh + 138);
    const auto *mh_139 = buffer.data(mh + 139);
    const auto *mh_140 = buffer.data(mh + 140);
    const auto *mh_141 = buffer.data(mh + 141);
    const auto *mh_142 = buffer.data(mh + 142);
    const auto *mh_143 = buffer.data(mh + 143);
    const auto *mh_144 = buffer.data(mh + 144);
    const auto *mh_145 = buffer.data(mh + 145);
    const auto *mh_146 = buffer.data(mh + 146);
    const auto *mh_147 = buffer.data(mh + 147);
    const auto *mh_148 = buffer.data(mh + 148);
    const auto *mh_149 = buffer.data(mh + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, kh_0, kh_1, kh_2, kh_3, kh_4, mh_0, mh_1, \
                         mh_2, mh_3, mh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -8.0 * kh_0[k]
                 + f_0 * mh_0[k];

        t_1[k] = -8.0 * kh_1[k]
                 + f_0 * mh_1[k];

        t_2[k] = -8.0 * kh_2[k]
                 + f_0 * mh_2[k];

        t_3[k] = -8.0 * kh_3[k]
                 + f_0 * mh_3[k];

        t_4[k] = -8.0 * kh_4[k]
                 + f_0 * mh_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, kh_5, kh_6, kh_7, kh_8, kh_9, mh_5, mh_6, \
                         mh_7, mh_8, mh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -8.0 * kh_5[k]
                 + f_0 * mh_5[k];

        t_6[k] = -8.0 * kh_6[k]
                 + f_0 * mh_6[k];

        t_7[k] = -8.0 * kh_7[k]
                 + f_0 * mh_7[k];

        t_8[k] = -8.0 * kh_8[k]
                 + f_0 * mh_8[k];

        t_9[k] = -8.0 * kh_9[k]
                 + f_0 * mh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, kh_10, kh_11, kh_12, kh_13, kh_14, \
                         mh_10, mh_11, mh_12, mh_13, mh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -8.0 * kh_10[k]
                  + f_0 * mh_10[k];

        t_11[k] = -8.0 * kh_11[k]
                  + f_0 * mh_11[k];

        t_12[k] = -8.0 * kh_12[k]
                  + f_0 * mh_12[k];

        t_13[k] = -8.0 * kh_13[k]
                  + f_0 * mh_13[k];

        t_14[k] = -8.0 * kh_14[k]
                  + f_0 * mh_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, kh_15, kh_16, kh_17, kh_18, kh_19, \
                         mh_15, mh_16, mh_17, mh_18, mh_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -8.0 * kh_15[k]
                  + f_0 * mh_15[k];

        t_16[k] = -8.0 * kh_16[k]
                  + f_0 * mh_16[k];

        t_17[k] = -8.0 * kh_17[k]
                  + f_0 * mh_17[k];

        t_18[k] = -8.0 * kh_18[k]
                  + f_0 * mh_18[k];

        t_19[k] = -8.0 * kh_19[k]
                  + f_0 * mh_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, kh_20, kh_21, kh_22, kh_23, kh_24, \
                         mh_20, mh_21, mh_22, mh_23, mh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -8.0 * kh_20[k]
                  + f_0 * mh_20[k];

        t_21[k] = -7.0 * kh_21[k]
                  + f_0 * mh_21[k];

        t_22[k] = -7.0 * kh_22[k]
                  + f_0 * mh_22[k];

        t_23[k] = -7.0 * kh_23[k]
                  + f_0 * mh_23[k];

        t_24[k] = -7.0 * kh_24[k]
                  + f_0 * mh_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, kh_25, kh_26, kh_27, kh_28, kh_29, \
                         mh_25, mh_26, mh_27, mh_28, mh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -7.0 * kh_25[k]
                  + f_0 * mh_25[k];

        t_26[k] = -7.0 * kh_26[k]
                  + f_0 * mh_26[k];

        t_27[k] = -7.0 * kh_27[k]
                  + f_0 * mh_27[k];

        t_28[k] = -7.0 * kh_28[k]
                  + f_0 * mh_28[k];

        t_29[k] = -7.0 * kh_29[k]
                  + f_0 * mh_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, kh_30, kh_31, kh_32, kh_33, kh_34, \
                         mh_30, mh_31, mh_32, mh_33, mh_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -7.0 * kh_30[k]
                  + f_0 * mh_30[k];

        t_31[k] = -7.0 * kh_31[k]
                  + f_0 * mh_31[k];

        t_32[k] = -7.0 * kh_32[k]
                  + f_0 * mh_32[k];

        t_33[k] = -7.0 * kh_33[k]
                  + f_0 * mh_33[k];

        t_34[k] = -7.0 * kh_34[k]
                  + f_0 * mh_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, kh_35, kh_36, kh_37, kh_38, kh_39, \
                         mh_35, mh_36, mh_37, mh_38, mh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -7.0 * kh_35[k]
                  + f_0 * mh_35[k];

        t_36[k] = -7.0 * kh_36[k]
                  + f_0 * mh_36[k];

        t_37[k] = -7.0 * kh_37[k]
                  + f_0 * mh_37[k];

        t_38[k] = -7.0 * kh_38[k]
                  + f_0 * mh_38[k];

        t_39[k] = -7.0 * kh_39[k]
                  + f_0 * mh_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, kh_40, kh_41, kh_42, kh_43, kh_44, \
                         mh_40, mh_41, mh_42, mh_43, mh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -7.0 * kh_40[k]
                  + f_0 * mh_40[k];

        t_41[k] = -7.0 * kh_41[k]
                  + f_0 * mh_41[k];

        t_42[k] = -7.0 * kh_42[k]
                  + f_0 * mh_42[k];

        t_43[k] = -7.0 * kh_43[k]
                  + f_0 * mh_43[k];

        t_44[k] = -7.0 * kh_44[k]
                  + f_0 * mh_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, kh_45, kh_46, kh_47, kh_48, kh_49, \
                         mh_45, mh_46, mh_47, mh_48, mh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -7.0 * kh_45[k]
                  + f_0 * mh_45[k];

        t_46[k] = -7.0 * kh_46[k]
                  + f_0 * mh_46[k];

        t_47[k] = -7.0 * kh_47[k]
                  + f_0 * mh_47[k];

        t_48[k] = -7.0 * kh_48[k]
                  + f_0 * mh_48[k];

        t_49[k] = -7.0 * kh_49[k]
                  + f_0 * mh_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, kh_50, kh_51, kh_52, kh_53, kh_54, \
                         mh_50, mh_51, mh_52, mh_53, mh_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -7.0 * kh_50[k]
                  + f_0 * mh_50[k];

        t_51[k] = -7.0 * kh_51[k]
                  + f_0 * mh_51[k];

        t_52[k] = -7.0 * kh_52[k]
                  + f_0 * mh_52[k];

        t_53[k] = -7.0 * kh_53[k]
                  + f_0 * mh_53[k];

        t_54[k] = -7.0 * kh_54[k]
                  + f_0 * mh_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, kh_55, kh_56, kh_57, kh_58, kh_59, \
                         mh_55, mh_56, mh_57, mh_58, mh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -7.0 * kh_55[k]
                  + f_0 * mh_55[k];

        t_56[k] = -7.0 * kh_56[k]
                  + f_0 * mh_56[k];

        t_57[k] = -7.0 * kh_57[k]
                  + f_0 * mh_57[k];

        t_58[k] = -7.0 * kh_58[k]
                  + f_0 * mh_58[k];

        t_59[k] = -7.0 * kh_59[k]
                  + f_0 * mh_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, kh_60, kh_61, kh_62, kh_63, kh_64, \
                         mh_60, mh_61, mh_62, mh_63, mh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -7.0 * kh_60[k]
                  + f_0 * mh_60[k];

        t_61[k] = -7.0 * kh_61[k]
                  + f_0 * mh_61[k];

        t_62[k] = -7.0 * kh_62[k]
                  + f_0 * mh_62[k];

        t_63[k] = -6.0 * kh_63[k]
                  + f_0 * mh_63[k];

        t_64[k] = -6.0 * kh_64[k]
                  + f_0 * mh_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, kh_65, kh_66, kh_67, kh_68, kh_69, \
                         mh_65, mh_66, mh_67, mh_68, mh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -6.0 * kh_65[k]
                  + f_0 * mh_65[k];

        t_66[k] = -6.0 * kh_66[k]
                  + f_0 * mh_66[k];

        t_67[k] = -6.0 * kh_67[k]
                  + f_0 * mh_67[k];

        t_68[k] = -6.0 * kh_68[k]
                  + f_0 * mh_68[k];

        t_69[k] = -6.0 * kh_69[k]
                  + f_0 * mh_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, kh_70, kh_71, kh_72, kh_73, kh_74, \
                         mh_70, mh_71, mh_72, mh_73, mh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -6.0 * kh_70[k]
                  + f_0 * mh_70[k];

        t_71[k] = -6.0 * kh_71[k]
                  + f_0 * mh_71[k];

        t_72[k] = -6.0 * kh_72[k]
                  + f_0 * mh_72[k];

        t_73[k] = -6.0 * kh_73[k]
                  + f_0 * mh_73[k];

        t_74[k] = -6.0 * kh_74[k]
                  + f_0 * mh_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, kh_75, kh_76, kh_77, kh_78, kh_79, \
                         mh_75, mh_76, mh_77, mh_78, mh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -6.0 * kh_75[k]
                  + f_0 * mh_75[k];

        t_76[k] = -6.0 * kh_76[k]
                  + f_0 * mh_76[k];

        t_77[k] = -6.0 * kh_77[k]
                  + f_0 * mh_77[k];

        t_78[k] = -6.0 * kh_78[k]
                  + f_0 * mh_78[k];

        t_79[k] = -6.0 * kh_79[k]
                  + f_0 * mh_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, kh_80, kh_81, kh_82, kh_83, kh_84, \
                         mh_80, mh_81, mh_82, mh_83, mh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -6.0 * kh_80[k]
                  + f_0 * mh_80[k];

        t_81[k] = -6.0 * kh_81[k]
                  + f_0 * mh_81[k];

        t_82[k] = -6.0 * kh_82[k]
                  + f_0 * mh_82[k];

        t_83[k] = -6.0 * kh_83[k]
                  + f_0 * mh_83[k];

        t_84[k] = -6.0 * kh_84[k]
                  + f_0 * mh_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, kh_85, kh_86, kh_87, kh_88, kh_89, \
                         mh_85, mh_86, mh_87, mh_88, mh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -6.0 * kh_85[k]
                  + f_0 * mh_85[k];

        t_86[k] = -6.0 * kh_86[k]
                  + f_0 * mh_86[k];

        t_87[k] = -6.0 * kh_87[k]
                  + f_0 * mh_87[k];

        t_88[k] = -6.0 * kh_88[k]
                  + f_0 * mh_88[k];

        t_89[k] = -6.0 * kh_89[k]
                  + f_0 * mh_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, kh_90, kh_91, kh_92, kh_93, kh_94, \
                         mh_90, mh_91, mh_92, mh_93, mh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -6.0 * kh_90[k]
                  + f_0 * mh_90[k];

        t_91[k] = -6.0 * kh_91[k]
                  + f_0 * mh_91[k];

        t_92[k] = -6.0 * kh_92[k]
                  + f_0 * mh_92[k];

        t_93[k] = -6.0 * kh_93[k]
                  + f_0 * mh_93[k];

        t_94[k] = -6.0 * kh_94[k]
                  + f_0 * mh_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, kh_95, kh_96, kh_97, kh_98, kh_99, \
                         mh_95, mh_96, mh_97, mh_98, mh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -6.0 * kh_95[k]
                  + f_0 * mh_95[k];

        t_96[k] = -6.0 * kh_96[k]
                  + f_0 * mh_96[k];

        t_97[k] = -6.0 * kh_97[k]
                  + f_0 * mh_97[k];

        t_98[k] = -6.0 * kh_98[k]
                  + f_0 * mh_98[k];

        t_99[k] = -6.0 * kh_99[k]
                  + f_0 * mh_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, kh_100, kh_101, kh_102, kh_103, \
                         kh_104, mh_100, mh_101, mh_102, mh_103, \
                         mh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -6.0 * kh_100[k]
                   + f_0 * mh_100[k];

        t_101[k] = -6.0 * kh_101[k]
                   + f_0 * mh_101[k];

        t_102[k] = -6.0 * kh_102[k]
                   + f_0 * mh_102[k];

        t_103[k] = -6.0 * kh_103[k]
                   + f_0 * mh_103[k];

        t_104[k] = -6.0 * kh_104[k]
                   + f_0 * mh_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, kh_105, kh_106, kh_107, kh_108, \
                         kh_109, mh_105, mh_106, mh_107, mh_108, \
                         mh_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -6.0 * kh_105[k]
                   + f_0 * mh_105[k];

        t_106[k] = -6.0 * kh_106[k]
                   + f_0 * mh_106[k];

        t_107[k] = -6.0 * kh_107[k]
                   + f_0 * mh_107[k];

        t_108[k] = -6.0 * kh_108[k]
                   + f_0 * mh_108[k];

        t_109[k] = -6.0 * kh_109[k]
                   + f_0 * mh_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, kh_110, kh_111, kh_112, kh_113, \
                         kh_114, mh_110, mh_111, mh_112, mh_113, \
                         mh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -6.0 * kh_110[k]
                   + f_0 * mh_110[k];

        t_111[k] = -6.0 * kh_111[k]
                   + f_0 * mh_111[k];

        t_112[k] = -6.0 * kh_112[k]
                   + f_0 * mh_112[k];

        t_113[k] = -6.0 * kh_113[k]
                   + f_0 * mh_113[k];

        t_114[k] = -6.0 * kh_114[k]
                   + f_0 * mh_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, kh_115, kh_116, kh_117, kh_118, \
                         kh_119, mh_115, mh_116, mh_117, mh_118, \
                         mh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -6.0 * kh_115[k]
                   + f_0 * mh_115[k];

        t_116[k] = -6.0 * kh_116[k]
                   + f_0 * mh_116[k];

        t_117[k] = -6.0 * kh_117[k]
                   + f_0 * mh_117[k];

        t_118[k] = -6.0 * kh_118[k]
                   + f_0 * mh_118[k];

        t_119[k] = -6.0 * kh_119[k]
                   + f_0 * mh_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, kh_120, kh_121, kh_122, kh_123, \
                         kh_124, mh_120, mh_121, mh_122, mh_123, \
                         mh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -6.0 * kh_120[k]
                   + f_0 * mh_120[k];

        t_121[k] = -6.0 * kh_121[k]
                   + f_0 * mh_121[k];

        t_122[k] = -6.0 * kh_122[k]
                   + f_0 * mh_122[k];

        t_123[k] = -6.0 * kh_123[k]
                   + f_0 * mh_123[k];

        t_124[k] = -6.0 * kh_124[k]
                   + f_0 * mh_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, kh_125, kh_126, kh_127, kh_128, \
                         kh_129, mh_125, mh_126, mh_127, mh_128, \
                         mh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -6.0 * kh_125[k]
                   + f_0 * mh_125[k];

        t_126[k] = -5.0 * kh_126[k]
                   + f_0 * mh_126[k];

        t_127[k] = -5.0 * kh_127[k]
                   + f_0 * mh_127[k];

        t_128[k] = -5.0 * kh_128[k]
                   + f_0 * mh_128[k];

        t_129[k] = -5.0 * kh_129[k]
                   + f_0 * mh_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, kh_130, kh_131, kh_132, kh_133, \
                         kh_134, mh_130, mh_131, mh_132, mh_133, \
                         mh_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -5.0 * kh_130[k]
                   + f_0 * mh_130[k];

        t_131[k] = -5.0 * kh_131[k]
                   + f_0 * mh_131[k];

        t_132[k] = -5.0 * kh_132[k]
                   + f_0 * mh_132[k];

        t_133[k] = -5.0 * kh_133[k]
                   + f_0 * mh_133[k];

        t_134[k] = -5.0 * kh_134[k]
                   + f_0 * mh_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, kh_135, kh_136, kh_137, kh_138, \
                         kh_139, mh_135, mh_136, mh_137, mh_138, \
                         mh_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -5.0 * kh_135[k]
                   + f_0 * mh_135[k];

        t_136[k] = -5.0 * kh_136[k]
                   + f_0 * mh_136[k];

        t_137[k] = -5.0 * kh_137[k]
                   + f_0 * mh_137[k];

        t_138[k] = -5.0 * kh_138[k]
                   + f_0 * mh_138[k];

        t_139[k] = -5.0 * kh_139[k]
                   + f_0 * mh_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, kh_140, kh_141, kh_142, kh_143, \
                         kh_144, mh_140, mh_141, mh_142, mh_143, \
                         mh_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -5.0 * kh_140[k]
                   + f_0 * mh_140[k];

        t_141[k] = -5.0 * kh_141[k]
                   + f_0 * mh_141[k];

        t_142[k] = -5.0 * kh_142[k]
                   + f_0 * mh_142[k];

        t_143[k] = -5.0 * kh_143[k]
                   + f_0 * mh_143[k];

        t_144[k] = -5.0 * kh_144[k]
                   + f_0 * mh_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, kh_145, kh_146, kh_147, kh_148, \
                         kh_149, mh_145, mh_146, mh_147, mh_148, \
                         mh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -5.0 * kh_145[k]
                   + f_0 * mh_145[k];

        t_146[k] = -5.0 * kh_146[k]
                   + f_0 * mh_146[k];

        t_147[k] = -5.0 * kh_147[k]
                   + f_0 * mh_147[k];

        t_148[k] = -5.0 * kh_148[k]
                   + f_0 * mh_148[k];

        t_149[k] = -5.0 * kh_149[k]
                   + f_0 * mh_149[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_150 = buffer.data(target + 150);
    auto *t_151 = buffer.data(target + 151);
    auto *t_152 = buffer.data(target + 152);
    auto *t_153 = buffer.data(target + 153);
    auto *t_154 = buffer.data(target + 154);
    auto *t_155 = buffer.data(target + 155);
    auto *t_156 = buffer.data(target + 156);
    auto *t_157 = buffer.data(target + 157);
    auto *t_158 = buffer.data(target + 158);
    auto *t_159 = buffer.data(target + 159);
    auto *t_160 = buffer.data(target + 160);
    auto *t_161 = buffer.data(target + 161);
    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);
    auto *t_168 = buffer.data(target + 168);
    auto *t_169 = buffer.data(target + 169);
    auto *t_170 = buffer.data(target + 170);
    auto *t_171 = buffer.data(target + 171);
    auto *t_172 = buffer.data(target + 172);
    auto *t_173 = buffer.data(target + 173);
    auto *t_174 = buffer.data(target + 174);
    auto *t_175 = buffer.data(target + 175);
    auto *t_176 = buffer.data(target + 176);
    auto *t_177 = buffer.data(target + 177);
    auto *t_178 = buffer.data(target + 178);
    auto *t_179 = buffer.data(target + 179);
    auto *t_180 = buffer.data(target + 180);
    auto *t_181 = buffer.data(target + 181);
    auto *t_182 = buffer.data(target + 182);
    auto *t_183 = buffer.data(target + 183);
    auto *t_184 = buffer.data(target + 184);
    auto *t_185 = buffer.data(target + 185);
    auto *t_186 = buffer.data(target + 186);
    auto *t_187 = buffer.data(target + 187);
    auto *t_188 = buffer.data(target + 188);
    auto *t_189 = buffer.data(target + 189);
    auto *t_190 = buffer.data(target + 190);
    auto *t_191 = buffer.data(target + 191);
    auto *t_192 = buffer.data(target + 192);
    auto *t_193 = buffer.data(target + 193);
    auto *t_194 = buffer.data(target + 194);
    auto *t_195 = buffer.data(target + 195);
    auto *t_196 = buffer.data(target + 196);
    auto *t_197 = buffer.data(target + 197);
    auto *t_198 = buffer.data(target + 198);
    auto *t_199 = buffer.data(target + 199);
    auto *t_200 = buffer.data(target + 200);
    auto *t_201 = buffer.data(target + 201);
    auto *t_202 = buffer.data(target + 202);
    auto *t_203 = buffer.data(target + 203);
    auto *t_204 = buffer.data(target + 204);
    auto *t_205 = buffer.data(target + 205);
    auto *t_206 = buffer.data(target + 206);
    auto *t_207 = buffer.data(target + 207);
    auto *t_208 = buffer.data(target + 208);
    auto *t_209 = buffer.data(target + 209);
    auto *t_210 = buffer.data(target + 210);
    auto *t_211 = buffer.data(target + 211);
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);
    auto *t_215 = buffer.data(target + 215);
    auto *t_216 = buffer.data(target + 216);
    auto *t_217 = buffer.data(target + 217);
    auto *t_218 = buffer.data(target + 218);
    auto *t_219 = buffer.data(target + 219);
    auto *t_220 = buffer.data(target + 220);
    auto *t_221 = buffer.data(target + 221);
    auto *t_222 = buffer.data(target + 222);
    auto *t_223 = buffer.data(target + 223);
    auto *t_224 = buffer.data(target + 224);
    auto *t_225 = buffer.data(target + 225);
    auto *t_226 = buffer.data(target + 226);
    auto *t_227 = buffer.data(target + 227);
    auto *t_228 = buffer.data(target + 228);
    auto *t_229 = buffer.data(target + 229);
    auto *t_230 = buffer.data(target + 230);
    auto *t_231 = buffer.data(target + 231);
    auto *t_232 = buffer.data(target + 232);
    auto *t_233 = buffer.data(target + 233);
    auto *t_234 = buffer.data(target + 234);
    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);
    auto *t_239 = buffer.data(target + 239);
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);
    auto *t_262 = buffer.data(target + 262);
    auto *t_263 = buffer.data(target + 263);
    auto *t_264 = buffer.data(target + 264);
    auto *t_265 = buffer.data(target + 265);
    auto *t_266 = buffer.data(target + 266);
    auto *t_267 = buffer.data(target + 267);
    auto *t_268 = buffer.data(target + 268);
    auto *t_269 = buffer.data(target + 269);
    auto *t_270 = buffer.data(target + 270);
    auto *t_271 = buffer.data(target + 271);
    auto *t_272 = buffer.data(target + 272);
    auto *t_273 = buffer.data(target + 273);
    auto *t_274 = buffer.data(target + 274);
    auto *t_275 = buffer.data(target + 275);
    auto *t_276 = buffer.data(target + 276);
    auto *t_277 = buffer.data(target + 277);
    auto *t_278 = buffer.data(target + 278);
    auto *t_279 = buffer.data(target + 279);
    auto *t_280 = buffer.data(target + 280);
    auto *t_281 = buffer.data(target + 281);
    auto *t_282 = buffer.data(target + 282);
    auto *t_283 = buffer.data(target + 283);
    auto *t_284 = buffer.data(target + 284);
    auto *t_285 = buffer.data(target + 285);
    auto *t_286 = buffer.data(target + 286);
    auto *t_287 = buffer.data(target + 287);
    auto *t_288 = buffer.data(target + 288);
    auto *t_289 = buffer.data(target + 289);
    auto *t_290 = buffer.data(target + 290);
    auto *t_291 = buffer.data(target + 291);
    auto *t_292 = buffer.data(target + 292);
    auto *t_293 = buffer.data(target + 293);
    auto *t_294 = buffer.data(target + 294);
    auto *t_295 = buffer.data(target + 295);
    auto *t_296 = buffer.data(target + 296);
    auto *t_297 = buffer.data(target + 297);
    auto *t_298 = buffer.data(target + 298);
    auto *t_299 = buffer.data(target + 299);

    const auto *kh_150 = buffer.data(kh + 150);
    const auto *kh_151 = buffer.data(kh + 151);
    const auto *kh_152 = buffer.data(kh + 152);
    const auto *kh_153 = buffer.data(kh + 153);
    const auto *kh_154 = buffer.data(kh + 154);
    const auto *kh_155 = buffer.data(kh + 155);
    const auto *kh_156 = buffer.data(kh + 156);
    const auto *kh_157 = buffer.data(kh + 157);
    const auto *kh_158 = buffer.data(kh + 158);
    const auto *kh_159 = buffer.data(kh + 159);
    const auto *kh_160 = buffer.data(kh + 160);
    const auto *kh_161 = buffer.data(kh + 161);
    const auto *kh_162 = buffer.data(kh + 162);
    const auto *kh_163 = buffer.data(kh + 163);
    const auto *kh_164 = buffer.data(kh + 164);
    const auto *kh_165 = buffer.data(kh + 165);
    const auto *kh_166 = buffer.data(kh + 166);
    const auto *kh_167 = buffer.data(kh + 167);
    const auto *kh_168 = buffer.data(kh + 168);
    const auto *kh_169 = buffer.data(kh + 169);
    const auto *kh_170 = buffer.data(kh + 170);
    const auto *kh_171 = buffer.data(kh + 171);
    const auto *kh_172 = buffer.data(kh + 172);
    const auto *kh_173 = buffer.data(kh + 173);
    const auto *kh_174 = buffer.data(kh + 174);
    const auto *kh_175 = buffer.data(kh + 175);
    const auto *kh_176 = buffer.data(kh + 176);
    const auto *kh_177 = buffer.data(kh + 177);
    const auto *kh_178 = buffer.data(kh + 178);
    const auto *kh_179 = buffer.data(kh + 179);
    const auto *kh_180 = buffer.data(kh + 180);
    const auto *kh_181 = buffer.data(kh + 181);
    const auto *kh_182 = buffer.data(kh + 182);
    const auto *kh_183 = buffer.data(kh + 183);
    const auto *kh_184 = buffer.data(kh + 184);
    const auto *kh_185 = buffer.data(kh + 185);
    const auto *kh_186 = buffer.data(kh + 186);
    const auto *kh_187 = buffer.data(kh + 187);
    const auto *kh_188 = buffer.data(kh + 188);
    const auto *kh_189 = buffer.data(kh + 189);
    const auto *kh_190 = buffer.data(kh + 190);
    const auto *kh_191 = buffer.data(kh + 191);
    const auto *kh_192 = buffer.data(kh + 192);
    const auto *kh_193 = buffer.data(kh + 193);
    const auto *kh_194 = buffer.data(kh + 194);
    const auto *kh_195 = buffer.data(kh + 195);
    const auto *kh_196 = buffer.data(kh + 196);
    const auto *kh_197 = buffer.data(kh + 197);
    const auto *kh_198 = buffer.data(kh + 198);
    const auto *kh_199 = buffer.data(kh + 199);
    const auto *kh_200 = buffer.data(kh + 200);
    const auto *kh_201 = buffer.data(kh + 201);
    const auto *kh_202 = buffer.data(kh + 202);
    const auto *kh_203 = buffer.data(kh + 203);
    const auto *kh_204 = buffer.data(kh + 204);
    const auto *kh_205 = buffer.data(kh + 205);
    const auto *kh_206 = buffer.data(kh + 206);
    const auto *kh_207 = buffer.data(kh + 207);
    const auto *kh_208 = buffer.data(kh + 208);
    const auto *kh_209 = buffer.data(kh + 209);
    const auto *kh_210 = buffer.data(kh + 210);
    const auto *kh_211 = buffer.data(kh + 211);
    const auto *kh_212 = buffer.data(kh + 212);
    const auto *kh_213 = buffer.data(kh + 213);
    const auto *kh_214 = buffer.data(kh + 214);
    const auto *kh_215 = buffer.data(kh + 215);
    const auto *kh_216 = buffer.data(kh + 216);
    const auto *kh_217 = buffer.data(kh + 217);
    const auto *kh_218 = buffer.data(kh + 218);
    const auto *kh_219 = buffer.data(kh + 219);
    const auto *kh_220 = buffer.data(kh + 220);
    const auto *kh_221 = buffer.data(kh + 221);
    const auto *kh_222 = buffer.data(kh + 222);
    const auto *kh_223 = buffer.data(kh + 223);
    const auto *kh_224 = buffer.data(kh + 224);
    const auto *kh_225 = buffer.data(kh + 225);
    const auto *kh_226 = buffer.data(kh + 226);
    const auto *kh_227 = buffer.data(kh + 227);
    const auto *kh_228 = buffer.data(kh + 228);
    const auto *kh_229 = buffer.data(kh + 229);
    const auto *kh_230 = buffer.data(kh + 230);
    const auto *kh_231 = buffer.data(kh + 231);
    const auto *kh_232 = buffer.data(kh + 232);
    const auto *kh_233 = buffer.data(kh + 233);
    const auto *kh_234 = buffer.data(kh + 234);
    const auto *kh_235 = buffer.data(kh + 235);
    const auto *kh_236 = buffer.data(kh + 236);
    const auto *kh_237 = buffer.data(kh + 237);
    const auto *kh_238 = buffer.data(kh + 238);
    const auto *kh_239 = buffer.data(kh + 239);
    const auto *kh_240 = buffer.data(kh + 240);
    const auto *kh_241 = buffer.data(kh + 241);
    const auto *kh_242 = buffer.data(kh + 242);
    const auto *kh_243 = buffer.data(kh + 243);
    const auto *kh_244 = buffer.data(kh + 244);
    const auto *kh_245 = buffer.data(kh + 245);
    const auto *kh_246 = buffer.data(kh + 246);
    const auto *kh_247 = buffer.data(kh + 247);
    const auto *kh_248 = buffer.data(kh + 248);
    const auto *kh_249 = buffer.data(kh + 249);
    const auto *kh_250 = buffer.data(kh + 250);
    const auto *kh_251 = buffer.data(kh + 251);
    const auto *kh_252 = buffer.data(kh + 252);
    const auto *kh_253 = buffer.data(kh + 253);
    const auto *kh_254 = buffer.data(kh + 254);
    const auto *kh_255 = buffer.data(kh + 255);
    const auto *kh_256 = buffer.data(kh + 256);
    const auto *kh_257 = buffer.data(kh + 257);
    const auto *kh_258 = buffer.data(kh + 258);
    const auto *kh_259 = buffer.data(kh + 259);
    const auto *kh_260 = buffer.data(kh + 260);
    const auto *kh_261 = buffer.data(kh + 261);
    const auto *kh_262 = buffer.data(kh + 262);
    const auto *kh_263 = buffer.data(kh + 263);
    const auto *kh_264 = buffer.data(kh + 264);
    const auto *kh_265 = buffer.data(kh + 265);
    const auto *kh_266 = buffer.data(kh + 266);
    const auto *kh_267 = buffer.data(kh + 267);
    const auto *kh_268 = buffer.data(kh + 268);
    const auto *kh_269 = buffer.data(kh + 269);
    const auto *kh_270 = buffer.data(kh + 270);
    const auto *kh_271 = buffer.data(kh + 271);
    const auto *kh_272 = buffer.data(kh + 272);
    const auto *kh_273 = buffer.data(kh + 273);
    const auto *kh_274 = buffer.data(kh + 274);
    const auto *kh_275 = buffer.data(kh + 275);
    const auto *kh_276 = buffer.data(kh + 276);
    const auto *kh_277 = buffer.data(kh + 277);
    const auto *kh_278 = buffer.data(kh + 278);
    const auto *kh_279 = buffer.data(kh + 279);
    const auto *kh_280 = buffer.data(kh + 280);
    const auto *kh_281 = buffer.data(kh + 281);
    const auto *kh_282 = buffer.data(kh + 282);
    const auto *kh_283 = buffer.data(kh + 283);
    const auto *kh_284 = buffer.data(kh + 284);
    const auto *kh_285 = buffer.data(kh + 285);
    const auto *kh_286 = buffer.data(kh + 286);
    const auto *kh_287 = buffer.data(kh + 287);
    const auto *kh_288 = buffer.data(kh + 288);
    const auto *kh_289 = buffer.data(kh + 289);
    const auto *kh_290 = buffer.data(kh + 290);
    const auto *kh_291 = buffer.data(kh + 291);
    const auto *kh_292 = buffer.data(kh + 292);
    const auto *kh_293 = buffer.data(kh + 293);
    const auto *kh_294 = buffer.data(kh + 294);
    const auto *kh_295 = buffer.data(kh + 295);
    const auto *kh_296 = buffer.data(kh + 296);
    const auto *kh_297 = buffer.data(kh + 297);
    const auto *kh_298 = buffer.data(kh + 298);
    const auto *kh_299 = buffer.data(kh + 299);

    const auto *mh_150 = buffer.data(mh + 150);
    const auto *mh_151 = buffer.data(mh + 151);
    const auto *mh_152 = buffer.data(mh + 152);
    const auto *mh_153 = buffer.data(mh + 153);
    const auto *mh_154 = buffer.data(mh + 154);
    const auto *mh_155 = buffer.data(mh + 155);
    const auto *mh_156 = buffer.data(mh + 156);
    const auto *mh_157 = buffer.data(mh + 157);
    const auto *mh_158 = buffer.data(mh + 158);
    const auto *mh_159 = buffer.data(mh + 159);
    const auto *mh_160 = buffer.data(mh + 160);
    const auto *mh_161 = buffer.data(mh + 161);
    const auto *mh_162 = buffer.data(mh + 162);
    const auto *mh_163 = buffer.data(mh + 163);
    const auto *mh_164 = buffer.data(mh + 164);
    const auto *mh_165 = buffer.data(mh + 165);
    const auto *mh_166 = buffer.data(mh + 166);
    const auto *mh_167 = buffer.data(mh + 167);
    const auto *mh_168 = buffer.data(mh + 168);
    const auto *mh_169 = buffer.data(mh + 169);
    const auto *mh_170 = buffer.data(mh + 170);
    const auto *mh_171 = buffer.data(mh + 171);
    const auto *mh_172 = buffer.data(mh + 172);
    const auto *mh_173 = buffer.data(mh + 173);
    const auto *mh_174 = buffer.data(mh + 174);
    const auto *mh_175 = buffer.data(mh + 175);
    const auto *mh_176 = buffer.data(mh + 176);
    const auto *mh_177 = buffer.data(mh + 177);
    const auto *mh_178 = buffer.data(mh + 178);
    const auto *mh_179 = buffer.data(mh + 179);
    const auto *mh_180 = buffer.data(mh + 180);
    const auto *mh_181 = buffer.data(mh + 181);
    const auto *mh_182 = buffer.data(mh + 182);
    const auto *mh_183 = buffer.data(mh + 183);
    const auto *mh_184 = buffer.data(mh + 184);
    const auto *mh_185 = buffer.data(mh + 185);
    const auto *mh_186 = buffer.data(mh + 186);
    const auto *mh_187 = buffer.data(mh + 187);
    const auto *mh_188 = buffer.data(mh + 188);
    const auto *mh_189 = buffer.data(mh + 189);
    const auto *mh_190 = buffer.data(mh + 190);
    const auto *mh_191 = buffer.data(mh + 191);
    const auto *mh_192 = buffer.data(mh + 192);
    const auto *mh_193 = buffer.data(mh + 193);
    const auto *mh_194 = buffer.data(mh + 194);
    const auto *mh_195 = buffer.data(mh + 195);
    const auto *mh_196 = buffer.data(mh + 196);
    const auto *mh_197 = buffer.data(mh + 197);
    const auto *mh_198 = buffer.data(mh + 198);
    const auto *mh_199 = buffer.data(mh + 199);
    const auto *mh_200 = buffer.data(mh + 200);
    const auto *mh_201 = buffer.data(mh + 201);
    const auto *mh_202 = buffer.data(mh + 202);
    const auto *mh_203 = buffer.data(mh + 203);
    const auto *mh_204 = buffer.data(mh + 204);
    const auto *mh_205 = buffer.data(mh + 205);
    const auto *mh_206 = buffer.data(mh + 206);
    const auto *mh_207 = buffer.data(mh + 207);
    const auto *mh_208 = buffer.data(mh + 208);
    const auto *mh_209 = buffer.data(mh + 209);
    const auto *mh_210 = buffer.data(mh + 210);
    const auto *mh_211 = buffer.data(mh + 211);
    const auto *mh_212 = buffer.data(mh + 212);
    const auto *mh_213 = buffer.data(mh + 213);
    const auto *mh_214 = buffer.data(mh + 214);
    const auto *mh_215 = buffer.data(mh + 215);
    const auto *mh_216 = buffer.data(mh + 216);
    const auto *mh_217 = buffer.data(mh + 217);
    const auto *mh_218 = buffer.data(mh + 218);
    const auto *mh_219 = buffer.data(mh + 219);
    const auto *mh_220 = buffer.data(mh + 220);
    const auto *mh_221 = buffer.data(mh + 221);
    const auto *mh_222 = buffer.data(mh + 222);
    const auto *mh_223 = buffer.data(mh + 223);
    const auto *mh_224 = buffer.data(mh + 224);
    const auto *mh_225 = buffer.data(mh + 225);
    const auto *mh_226 = buffer.data(mh + 226);
    const auto *mh_227 = buffer.data(mh + 227);
    const auto *mh_228 = buffer.data(mh + 228);
    const auto *mh_229 = buffer.data(mh + 229);
    const auto *mh_230 = buffer.data(mh + 230);
    const auto *mh_231 = buffer.data(mh + 231);
    const auto *mh_232 = buffer.data(mh + 232);
    const auto *mh_233 = buffer.data(mh + 233);
    const auto *mh_234 = buffer.data(mh + 234);
    const auto *mh_235 = buffer.data(mh + 235);
    const auto *mh_236 = buffer.data(mh + 236);
    const auto *mh_237 = buffer.data(mh + 237);
    const auto *mh_238 = buffer.data(mh + 238);
    const auto *mh_239 = buffer.data(mh + 239);
    const auto *mh_240 = buffer.data(mh + 240);
    const auto *mh_241 = buffer.data(mh + 241);
    const auto *mh_242 = buffer.data(mh + 242);
    const auto *mh_243 = buffer.data(mh + 243);
    const auto *mh_244 = buffer.data(mh + 244);
    const auto *mh_245 = buffer.data(mh + 245);
    const auto *mh_246 = buffer.data(mh + 246);
    const auto *mh_247 = buffer.data(mh + 247);
    const auto *mh_248 = buffer.data(mh + 248);
    const auto *mh_249 = buffer.data(mh + 249);
    const auto *mh_250 = buffer.data(mh + 250);
    const auto *mh_251 = buffer.data(mh + 251);
    const auto *mh_252 = buffer.data(mh + 252);
    const auto *mh_253 = buffer.data(mh + 253);
    const auto *mh_254 = buffer.data(mh + 254);
    const auto *mh_255 = buffer.data(mh + 255);
    const auto *mh_256 = buffer.data(mh + 256);
    const auto *mh_257 = buffer.data(mh + 257);
    const auto *mh_258 = buffer.data(mh + 258);
    const auto *mh_259 = buffer.data(mh + 259);
    const auto *mh_260 = buffer.data(mh + 260);
    const auto *mh_261 = buffer.data(mh + 261);
    const auto *mh_262 = buffer.data(mh + 262);
    const auto *mh_263 = buffer.data(mh + 263);
    const auto *mh_264 = buffer.data(mh + 264);
    const auto *mh_265 = buffer.data(mh + 265);
    const auto *mh_266 = buffer.data(mh + 266);
    const auto *mh_267 = buffer.data(mh + 267);
    const auto *mh_268 = buffer.data(mh + 268);
    const auto *mh_269 = buffer.data(mh + 269);
    const auto *mh_270 = buffer.data(mh + 270);
    const auto *mh_271 = buffer.data(mh + 271);
    const auto *mh_272 = buffer.data(mh + 272);
    const auto *mh_273 = buffer.data(mh + 273);
    const auto *mh_274 = buffer.data(mh + 274);
    const auto *mh_275 = buffer.data(mh + 275);
    const auto *mh_276 = buffer.data(mh + 276);
    const auto *mh_277 = buffer.data(mh + 277);
    const auto *mh_278 = buffer.data(mh + 278);
    const auto *mh_279 = buffer.data(mh + 279);
    const auto *mh_280 = buffer.data(mh + 280);
    const auto *mh_281 = buffer.data(mh + 281);
    const auto *mh_282 = buffer.data(mh + 282);
    const auto *mh_283 = buffer.data(mh + 283);
    const auto *mh_284 = buffer.data(mh + 284);
    const auto *mh_285 = buffer.data(mh + 285);
    const auto *mh_286 = buffer.data(mh + 286);
    const auto *mh_287 = buffer.data(mh + 287);
    const auto *mh_288 = buffer.data(mh + 288);
    const auto *mh_289 = buffer.data(mh + 289);
    const auto *mh_290 = buffer.data(mh + 290);
    const auto *mh_291 = buffer.data(mh + 291);
    const auto *mh_292 = buffer.data(mh + 292);
    const auto *mh_293 = buffer.data(mh + 293);
    const auto *mh_294 = buffer.data(mh + 294);
    const auto *mh_295 = buffer.data(mh + 295);
    const auto *mh_296 = buffer.data(mh + 296);
    const auto *mh_297 = buffer.data(mh + 297);
    const auto *mh_298 = buffer.data(mh + 298);
    const auto *mh_299 = buffer.data(mh + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, kh_150, kh_151, kh_152, kh_153, \
                         kh_154, mh_150, mh_151, mh_152, mh_153, \
                         mh_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -5.0 * kh_150[k]
                   + f_0 * mh_150[k];

        t_151[k] = -5.0 * kh_151[k]
                   + f_0 * mh_151[k];

        t_152[k] = -5.0 * kh_152[k]
                   + f_0 * mh_152[k];

        t_153[k] = -5.0 * kh_153[k]
                   + f_0 * mh_153[k];

        t_154[k] = -5.0 * kh_154[k]
                   + f_0 * mh_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, kh_155, kh_156, kh_157, kh_158, \
                         kh_159, mh_155, mh_156, mh_157, mh_158, \
                         mh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -5.0 * kh_155[k]
                   + f_0 * mh_155[k];

        t_156[k] = -5.0 * kh_156[k]
                   + f_0 * mh_156[k];

        t_157[k] = -5.0 * kh_157[k]
                   + f_0 * mh_157[k];

        t_158[k] = -5.0 * kh_158[k]
                   + f_0 * mh_158[k];

        t_159[k] = -5.0 * kh_159[k]
                   + f_0 * mh_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, kh_160, kh_161, kh_162, kh_163, \
                         kh_164, mh_160, mh_161, mh_162, mh_163, \
                         mh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -5.0 * kh_160[k]
                   + f_0 * mh_160[k];

        t_161[k] = -5.0 * kh_161[k]
                   + f_0 * mh_161[k];

        t_162[k] = -5.0 * kh_162[k]
                   + f_0 * mh_162[k];

        t_163[k] = -5.0 * kh_163[k]
                   + f_0 * mh_163[k];

        t_164[k] = -5.0 * kh_164[k]
                   + f_0 * mh_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, kh_165, kh_166, kh_167, kh_168, \
                         kh_169, mh_165, mh_166, mh_167, mh_168, \
                         mh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -5.0 * kh_165[k]
                   + f_0 * mh_165[k];

        t_166[k] = -5.0 * kh_166[k]
                   + f_0 * mh_166[k];

        t_167[k] = -5.0 * kh_167[k]
                   + f_0 * mh_167[k];

        t_168[k] = -5.0 * kh_168[k]
                   + f_0 * mh_168[k];

        t_169[k] = -5.0 * kh_169[k]
                   + f_0 * mh_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, kh_170, kh_171, kh_172, kh_173, \
                         kh_174, mh_170, mh_171, mh_172, mh_173, \
                         mh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -5.0 * kh_170[k]
                   + f_0 * mh_170[k];

        t_171[k] = -5.0 * kh_171[k]
                   + f_0 * mh_171[k];

        t_172[k] = -5.0 * kh_172[k]
                   + f_0 * mh_172[k];

        t_173[k] = -5.0 * kh_173[k]
                   + f_0 * mh_173[k];

        t_174[k] = -5.0 * kh_174[k]
                   + f_0 * mh_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, kh_175, kh_176, kh_177, kh_178, \
                         kh_179, mh_175, mh_176, mh_177, mh_178, \
                         mh_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -5.0 * kh_175[k]
                   + f_0 * mh_175[k];

        t_176[k] = -5.0 * kh_176[k]
                   + f_0 * mh_176[k];

        t_177[k] = -5.0 * kh_177[k]
                   + f_0 * mh_177[k];

        t_178[k] = -5.0 * kh_178[k]
                   + f_0 * mh_178[k];

        t_179[k] = -5.0 * kh_179[k]
                   + f_0 * mh_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, kh_180, kh_181, kh_182, kh_183, \
                         kh_184, mh_180, mh_181, mh_182, mh_183, \
                         mh_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -5.0 * kh_180[k]
                   + f_0 * mh_180[k];

        t_181[k] = -5.0 * kh_181[k]
                   + f_0 * mh_181[k];

        t_182[k] = -5.0 * kh_182[k]
                   + f_0 * mh_182[k];

        t_183[k] = -5.0 * kh_183[k]
                   + f_0 * mh_183[k];

        t_184[k] = -5.0 * kh_184[k]
                   + f_0 * mh_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, kh_185, kh_186, kh_187, kh_188, \
                         kh_189, mh_185, mh_186, mh_187, mh_188, \
                         mh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -5.0 * kh_185[k]
                   + f_0 * mh_185[k];

        t_186[k] = -5.0 * kh_186[k]
                   + f_0 * mh_186[k];

        t_187[k] = -5.0 * kh_187[k]
                   + f_0 * mh_187[k];

        t_188[k] = -5.0 * kh_188[k]
                   + f_0 * mh_188[k];

        t_189[k] = -5.0 * kh_189[k]
                   + f_0 * mh_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, kh_190, kh_191, kh_192, kh_193, \
                         kh_194, mh_190, mh_191, mh_192, mh_193, \
                         mh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -5.0 * kh_190[k]
                   + f_0 * mh_190[k];

        t_191[k] = -5.0 * kh_191[k]
                   + f_0 * mh_191[k];

        t_192[k] = -5.0 * kh_192[k]
                   + f_0 * mh_192[k];

        t_193[k] = -5.0 * kh_193[k]
                   + f_0 * mh_193[k];

        t_194[k] = -5.0 * kh_194[k]
                   + f_0 * mh_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, kh_195, kh_196, kh_197, kh_198, \
                         kh_199, mh_195, mh_196, mh_197, mh_198, \
                         mh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -5.0 * kh_195[k]
                   + f_0 * mh_195[k];

        t_196[k] = -5.0 * kh_196[k]
                   + f_0 * mh_196[k];

        t_197[k] = -5.0 * kh_197[k]
                   + f_0 * mh_197[k];

        t_198[k] = -5.0 * kh_198[k]
                   + f_0 * mh_198[k];

        t_199[k] = -5.0 * kh_199[k]
                   + f_0 * mh_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, kh_200, kh_201, kh_202, kh_203, \
                         kh_204, mh_200, mh_201, mh_202, mh_203, \
                         mh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -5.0 * kh_200[k]
                   + f_0 * mh_200[k];

        t_201[k] = -5.0 * kh_201[k]
                   + f_0 * mh_201[k];

        t_202[k] = -5.0 * kh_202[k]
                   + f_0 * mh_202[k];

        t_203[k] = -5.0 * kh_203[k]
                   + f_0 * mh_203[k];

        t_204[k] = -5.0 * kh_204[k]
                   + f_0 * mh_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, kh_205, kh_206, kh_207, kh_208, \
                         kh_209, mh_205, mh_206, mh_207, mh_208, \
                         mh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -5.0 * kh_205[k]
                   + f_0 * mh_205[k];

        t_206[k] = -5.0 * kh_206[k]
                   + f_0 * mh_206[k];

        t_207[k] = -5.0 * kh_207[k]
                   + f_0 * mh_207[k];

        t_208[k] = -5.0 * kh_208[k]
                   + f_0 * mh_208[k];

        t_209[k] = -5.0 * kh_209[k]
                   + f_0 * mh_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, kh_210, kh_211, kh_212, kh_213, \
                         kh_214, mh_210, mh_211, mh_212, mh_213, \
                         mh_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -4.0 * kh_210[k]
                   + f_0 * mh_210[k];

        t_211[k] = -4.0 * kh_211[k]
                   + f_0 * mh_211[k];

        t_212[k] = -4.0 * kh_212[k]
                   + f_0 * mh_212[k];

        t_213[k] = -4.0 * kh_213[k]
                   + f_0 * mh_213[k];

        t_214[k] = -4.0 * kh_214[k]
                   + f_0 * mh_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, kh_215, kh_216, kh_217, kh_218, \
                         kh_219, mh_215, mh_216, mh_217, mh_218, \
                         mh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -4.0 * kh_215[k]
                   + f_0 * mh_215[k];

        t_216[k] = -4.0 * kh_216[k]
                   + f_0 * mh_216[k];

        t_217[k] = -4.0 * kh_217[k]
                   + f_0 * mh_217[k];

        t_218[k] = -4.0 * kh_218[k]
                   + f_0 * mh_218[k];

        t_219[k] = -4.0 * kh_219[k]
                   + f_0 * mh_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, kh_220, kh_221, kh_222, kh_223, \
                         kh_224, mh_220, mh_221, mh_222, mh_223, \
                         mh_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -4.0 * kh_220[k]
                   + f_0 * mh_220[k];

        t_221[k] = -4.0 * kh_221[k]
                   + f_0 * mh_221[k];

        t_222[k] = -4.0 * kh_222[k]
                   + f_0 * mh_222[k];

        t_223[k] = -4.0 * kh_223[k]
                   + f_0 * mh_223[k];

        t_224[k] = -4.0 * kh_224[k]
                   + f_0 * mh_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, kh_225, kh_226, kh_227, kh_228, \
                         kh_229, mh_225, mh_226, mh_227, mh_228, \
                         mh_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -4.0 * kh_225[k]
                   + f_0 * mh_225[k];

        t_226[k] = -4.0 * kh_226[k]
                   + f_0 * mh_226[k];

        t_227[k] = -4.0 * kh_227[k]
                   + f_0 * mh_227[k];

        t_228[k] = -4.0 * kh_228[k]
                   + f_0 * mh_228[k];

        t_229[k] = -4.0 * kh_229[k]
                   + f_0 * mh_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, kh_230, kh_231, kh_232, kh_233, \
                         kh_234, mh_230, mh_231, mh_232, mh_233, \
                         mh_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -4.0 * kh_230[k]
                   + f_0 * mh_230[k];

        t_231[k] = -4.0 * kh_231[k]
                   + f_0 * mh_231[k];

        t_232[k] = -4.0 * kh_232[k]
                   + f_0 * mh_232[k];

        t_233[k] = -4.0 * kh_233[k]
                   + f_0 * mh_233[k];

        t_234[k] = -4.0 * kh_234[k]
                   + f_0 * mh_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, kh_235, kh_236, kh_237, kh_238, \
                         kh_239, mh_235, mh_236, mh_237, mh_238, \
                         mh_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -4.0 * kh_235[k]
                   + f_0 * mh_235[k];

        t_236[k] = -4.0 * kh_236[k]
                   + f_0 * mh_236[k];

        t_237[k] = -4.0 * kh_237[k]
                   + f_0 * mh_237[k];

        t_238[k] = -4.0 * kh_238[k]
                   + f_0 * mh_238[k];

        t_239[k] = -4.0 * kh_239[k]
                   + f_0 * mh_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, kh_240, kh_241, kh_242, kh_243, \
                         kh_244, mh_240, mh_241, mh_242, mh_243, \
                         mh_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -4.0 * kh_240[k]
                   + f_0 * mh_240[k];

        t_241[k] = -4.0 * kh_241[k]
                   + f_0 * mh_241[k];

        t_242[k] = -4.0 * kh_242[k]
                   + f_0 * mh_242[k];

        t_243[k] = -4.0 * kh_243[k]
                   + f_0 * mh_243[k];

        t_244[k] = -4.0 * kh_244[k]
                   + f_0 * mh_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, kh_245, kh_246, kh_247, kh_248, \
                         kh_249, mh_245, mh_246, mh_247, mh_248, \
                         mh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -4.0 * kh_245[k]
                   + f_0 * mh_245[k];

        t_246[k] = -4.0 * kh_246[k]
                   + f_0 * mh_246[k];

        t_247[k] = -4.0 * kh_247[k]
                   + f_0 * mh_247[k];

        t_248[k] = -4.0 * kh_248[k]
                   + f_0 * mh_248[k];

        t_249[k] = -4.0 * kh_249[k]
                   + f_0 * mh_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, kh_250, kh_251, kh_252, kh_253, \
                         kh_254, mh_250, mh_251, mh_252, mh_253, \
                         mh_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -4.0 * kh_250[k]
                   + f_0 * mh_250[k];

        t_251[k] = -4.0 * kh_251[k]
                   + f_0 * mh_251[k];

        t_252[k] = -4.0 * kh_252[k]
                   + f_0 * mh_252[k];

        t_253[k] = -4.0 * kh_253[k]
                   + f_0 * mh_253[k];

        t_254[k] = -4.0 * kh_254[k]
                   + f_0 * mh_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, kh_255, kh_256, kh_257, kh_258, \
                         kh_259, mh_255, mh_256, mh_257, mh_258, \
                         mh_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -4.0 * kh_255[k]
                   + f_0 * mh_255[k];

        t_256[k] = -4.0 * kh_256[k]
                   + f_0 * mh_256[k];

        t_257[k] = -4.0 * kh_257[k]
                   + f_0 * mh_257[k];

        t_258[k] = -4.0 * kh_258[k]
                   + f_0 * mh_258[k];

        t_259[k] = -4.0 * kh_259[k]
                   + f_0 * mh_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, kh_260, kh_261, kh_262, kh_263, \
                         kh_264, mh_260, mh_261, mh_262, mh_263, \
                         mh_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -4.0 * kh_260[k]
                   + f_0 * mh_260[k];

        t_261[k] = -4.0 * kh_261[k]
                   + f_0 * mh_261[k];

        t_262[k] = -4.0 * kh_262[k]
                   + f_0 * mh_262[k];

        t_263[k] = -4.0 * kh_263[k]
                   + f_0 * mh_263[k];

        t_264[k] = -4.0 * kh_264[k]
                   + f_0 * mh_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, kh_265, kh_266, kh_267, kh_268, \
                         kh_269, mh_265, mh_266, mh_267, mh_268, \
                         mh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -4.0 * kh_265[k]
                   + f_0 * mh_265[k];

        t_266[k] = -4.0 * kh_266[k]
                   + f_0 * mh_266[k];

        t_267[k] = -4.0 * kh_267[k]
                   + f_0 * mh_267[k];

        t_268[k] = -4.0 * kh_268[k]
                   + f_0 * mh_268[k];

        t_269[k] = -4.0 * kh_269[k]
                   + f_0 * mh_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, kh_270, kh_271, kh_272, kh_273, \
                         kh_274, mh_270, mh_271, mh_272, mh_273, \
                         mh_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -4.0 * kh_270[k]
                   + f_0 * mh_270[k];

        t_271[k] = -4.0 * kh_271[k]
                   + f_0 * mh_271[k];

        t_272[k] = -4.0 * kh_272[k]
                   + f_0 * mh_272[k];

        t_273[k] = -4.0 * kh_273[k]
                   + f_0 * mh_273[k];

        t_274[k] = -4.0 * kh_274[k]
                   + f_0 * mh_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, kh_275, kh_276, kh_277, kh_278, \
                         kh_279, mh_275, mh_276, mh_277, mh_278, \
                         mh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -4.0 * kh_275[k]
                   + f_0 * mh_275[k];

        t_276[k] = -4.0 * kh_276[k]
                   + f_0 * mh_276[k];

        t_277[k] = -4.0 * kh_277[k]
                   + f_0 * mh_277[k];

        t_278[k] = -4.0 * kh_278[k]
                   + f_0 * mh_278[k];

        t_279[k] = -4.0 * kh_279[k]
                   + f_0 * mh_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, kh_280, kh_281, kh_282, kh_283, \
                         kh_284, mh_280, mh_281, mh_282, mh_283, \
                         mh_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -4.0 * kh_280[k]
                   + f_0 * mh_280[k];

        t_281[k] = -4.0 * kh_281[k]
                   + f_0 * mh_281[k];

        t_282[k] = -4.0 * kh_282[k]
                   + f_0 * mh_282[k];

        t_283[k] = -4.0 * kh_283[k]
                   + f_0 * mh_283[k];

        t_284[k] = -4.0 * kh_284[k]
                   + f_0 * mh_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, kh_285, kh_286, kh_287, kh_288, \
                         kh_289, mh_285, mh_286, mh_287, mh_288, \
                         mh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -4.0 * kh_285[k]
                   + f_0 * mh_285[k];

        t_286[k] = -4.0 * kh_286[k]
                   + f_0 * mh_286[k];

        t_287[k] = -4.0 * kh_287[k]
                   + f_0 * mh_287[k];

        t_288[k] = -4.0 * kh_288[k]
                   + f_0 * mh_288[k];

        t_289[k] = -4.0 * kh_289[k]
                   + f_0 * mh_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, kh_290, kh_291, kh_292, kh_293, \
                         kh_294, mh_290, mh_291, mh_292, mh_293, \
                         mh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -4.0 * kh_290[k]
                   + f_0 * mh_290[k];

        t_291[k] = -4.0 * kh_291[k]
                   + f_0 * mh_291[k];

        t_292[k] = -4.0 * kh_292[k]
                   + f_0 * mh_292[k];

        t_293[k] = -4.0 * kh_293[k]
                   + f_0 * mh_293[k];

        t_294[k] = -4.0 * kh_294[k]
                   + f_0 * mh_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, kh_295, kh_296, kh_297, kh_298, \
                         kh_299, mh_295, mh_296, mh_297, mh_298, \
                         mh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * kh_295[k]
                   + f_0 * mh_295[k];

        t_296[k] = -4.0 * kh_296[k]
                   + f_0 * mh_296[k];

        t_297[k] = -4.0 * kh_297[k]
                   + f_0 * mh_297[k];

        t_298[k] = -4.0 * kh_298[k]
                   + f_0 * mh_298[k];

        t_299[k] = -4.0 * kh_299[k]
                   + f_0 * mh_299[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_300 = buffer.data(target + 300);
    auto *t_301 = buffer.data(target + 301);
    auto *t_302 = buffer.data(target + 302);
    auto *t_303 = buffer.data(target + 303);
    auto *t_304 = buffer.data(target + 304);
    auto *t_305 = buffer.data(target + 305);
    auto *t_306 = buffer.data(target + 306);
    auto *t_307 = buffer.data(target + 307);
    auto *t_308 = buffer.data(target + 308);
    auto *t_309 = buffer.data(target + 309);
    auto *t_310 = buffer.data(target + 310);
    auto *t_311 = buffer.data(target + 311);
    auto *t_312 = buffer.data(target + 312);
    auto *t_313 = buffer.data(target + 313);
    auto *t_314 = buffer.data(target + 314);
    auto *t_315 = buffer.data(target + 315);
    auto *t_316 = buffer.data(target + 316);
    auto *t_317 = buffer.data(target + 317);
    auto *t_318 = buffer.data(target + 318);
    auto *t_319 = buffer.data(target + 319);
    auto *t_320 = buffer.data(target + 320);
    auto *t_321 = buffer.data(target + 321);
    auto *t_322 = buffer.data(target + 322);
    auto *t_323 = buffer.data(target + 323);
    auto *t_324 = buffer.data(target + 324);
    auto *t_325 = buffer.data(target + 325);
    auto *t_326 = buffer.data(target + 326);
    auto *t_327 = buffer.data(target + 327);
    auto *t_328 = buffer.data(target + 328);
    auto *t_329 = buffer.data(target + 329);
    auto *t_330 = buffer.data(target + 330);
    auto *t_331 = buffer.data(target + 331);
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);
    auto *t_334 = buffer.data(target + 334);
    auto *t_335 = buffer.data(target + 335);
    auto *t_336 = buffer.data(target + 336);
    auto *t_337 = buffer.data(target + 337);
    auto *t_338 = buffer.data(target + 338);
    auto *t_339 = buffer.data(target + 339);
    auto *t_340 = buffer.data(target + 340);
    auto *t_341 = buffer.data(target + 341);
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);
    auto *t_347 = buffer.data(target + 347);
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);
    auto *t_380 = buffer.data(target + 380);
    auto *t_381 = buffer.data(target + 381);
    auto *t_382 = buffer.data(target + 382);
    auto *t_383 = buffer.data(target + 383);
    auto *t_384 = buffer.data(target + 384);
    auto *t_385 = buffer.data(target + 385);
    auto *t_386 = buffer.data(target + 386);
    auto *t_387 = buffer.data(target + 387);
    auto *t_388 = buffer.data(target + 388);
    auto *t_389 = buffer.data(target + 389);
    auto *t_390 = buffer.data(target + 390);
    auto *t_391 = buffer.data(target + 391);
    auto *t_392 = buffer.data(target + 392);
    auto *t_393 = buffer.data(target + 393);
    auto *t_394 = buffer.data(target + 394);
    auto *t_395 = buffer.data(target + 395);
    auto *t_396 = buffer.data(target + 396);
    auto *t_397 = buffer.data(target + 397);
    auto *t_398 = buffer.data(target + 398);
    auto *t_399 = buffer.data(target + 399);
    auto *t_400 = buffer.data(target + 400);
    auto *t_401 = buffer.data(target + 401);
    auto *t_402 = buffer.data(target + 402);
    auto *t_403 = buffer.data(target + 403);
    auto *t_404 = buffer.data(target + 404);
    auto *t_405 = buffer.data(target + 405);
    auto *t_406 = buffer.data(target + 406);
    auto *t_407 = buffer.data(target + 407);
    auto *t_408 = buffer.data(target + 408);
    auto *t_409 = buffer.data(target + 409);
    auto *t_410 = buffer.data(target + 410);
    auto *t_411 = buffer.data(target + 411);
    auto *t_412 = buffer.data(target + 412);
    auto *t_413 = buffer.data(target + 413);
    auto *t_414 = buffer.data(target + 414);
    auto *t_415 = buffer.data(target + 415);
    auto *t_416 = buffer.data(target + 416);
    auto *t_417 = buffer.data(target + 417);
    auto *t_418 = buffer.data(target + 418);
    auto *t_419 = buffer.data(target + 419);
    auto *t_420 = buffer.data(target + 420);
    auto *t_421 = buffer.data(target + 421);
    auto *t_422 = buffer.data(target + 422);
    auto *t_423 = buffer.data(target + 423);
    auto *t_424 = buffer.data(target + 424);
    auto *t_425 = buffer.data(target + 425);
    auto *t_426 = buffer.data(target + 426);
    auto *t_427 = buffer.data(target + 427);
    auto *t_428 = buffer.data(target + 428);
    auto *t_429 = buffer.data(target + 429);
    auto *t_430 = buffer.data(target + 430);
    auto *t_431 = buffer.data(target + 431);
    auto *t_432 = buffer.data(target + 432);
    auto *t_433 = buffer.data(target + 433);
    auto *t_434 = buffer.data(target + 434);
    auto *t_435 = buffer.data(target + 435);
    auto *t_436 = buffer.data(target + 436);
    auto *t_437 = buffer.data(target + 437);
    auto *t_438 = buffer.data(target + 438);
    auto *t_439 = buffer.data(target + 439);
    auto *t_440 = buffer.data(target + 440);
    auto *t_441 = buffer.data(target + 441);
    auto *t_442 = buffer.data(target + 442);
    auto *t_443 = buffer.data(target + 443);
    auto *t_444 = buffer.data(target + 444);
    auto *t_445 = buffer.data(target + 445);
    auto *t_446 = buffer.data(target + 446);
    auto *t_447 = buffer.data(target + 447);
    auto *t_448 = buffer.data(target + 448);
    auto *t_449 = buffer.data(target + 449);

    const auto *kh_300 = buffer.data(kh + 300);
    const auto *kh_301 = buffer.data(kh + 301);
    const auto *kh_302 = buffer.data(kh + 302);
    const auto *kh_303 = buffer.data(kh + 303);
    const auto *kh_304 = buffer.data(kh + 304);
    const auto *kh_305 = buffer.data(kh + 305);
    const auto *kh_306 = buffer.data(kh + 306);
    const auto *kh_307 = buffer.data(kh + 307);
    const auto *kh_308 = buffer.data(kh + 308);
    const auto *kh_309 = buffer.data(kh + 309);
    const auto *kh_310 = buffer.data(kh + 310);
    const auto *kh_311 = buffer.data(kh + 311);
    const auto *kh_312 = buffer.data(kh + 312);
    const auto *kh_313 = buffer.data(kh + 313);
    const auto *kh_314 = buffer.data(kh + 314);
    const auto *kh_315 = buffer.data(kh + 315);
    const auto *kh_316 = buffer.data(kh + 316);
    const auto *kh_317 = buffer.data(kh + 317);
    const auto *kh_318 = buffer.data(kh + 318);
    const auto *kh_319 = buffer.data(kh + 319);
    const auto *kh_320 = buffer.data(kh + 320);
    const auto *kh_321 = buffer.data(kh + 321);
    const auto *kh_322 = buffer.data(kh + 322);
    const auto *kh_323 = buffer.data(kh + 323);
    const auto *kh_324 = buffer.data(kh + 324);
    const auto *kh_325 = buffer.data(kh + 325);
    const auto *kh_326 = buffer.data(kh + 326);
    const auto *kh_327 = buffer.data(kh + 327);
    const auto *kh_328 = buffer.data(kh + 328);
    const auto *kh_329 = buffer.data(kh + 329);
    const auto *kh_330 = buffer.data(kh + 330);
    const auto *kh_331 = buffer.data(kh + 331);
    const auto *kh_332 = buffer.data(kh + 332);
    const auto *kh_333 = buffer.data(kh + 333);
    const auto *kh_334 = buffer.data(kh + 334);
    const auto *kh_335 = buffer.data(kh + 335);
    const auto *kh_336 = buffer.data(kh + 336);
    const auto *kh_337 = buffer.data(kh + 337);
    const auto *kh_338 = buffer.data(kh + 338);
    const auto *kh_339 = buffer.data(kh + 339);
    const auto *kh_340 = buffer.data(kh + 340);
    const auto *kh_341 = buffer.data(kh + 341);
    const auto *kh_342 = buffer.data(kh + 342);
    const auto *kh_343 = buffer.data(kh + 343);
    const auto *kh_344 = buffer.data(kh + 344);
    const auto *kh_345 = buffer.data(kh + 345);
    const auto *kh_346 = buffer.data(kh + 346);
    const auto *kh_347 = buffer.data(kh + 347);
    const auto *kh_348 = buffer.data(kh + 348);
    const auto *kh_349 = buffer.data(kh + 349);
    const auto *kh_350 = buffer.data(kh + 350);
    const auto *kh_351 = buffer.data(kh + 351);
    const auto *kh_352 = buffer.data(kh + 352);
    const auto *kh_353 = buffer.data(kh + 353);
    const auto *kh_354 = buffer.data(kh + 354);
    const auto *kh_355 = buffer.data(kh + 355);
    const auto *kh_356 = buffer.data(kh + 356);
    const auto *kh_357 = buffer.data(kh + 357);
    const auto *kh_358 = buffer.data(kh + 358);
    const auto *kh_359 = buffer.data(kh + 359);
    const auto *kh_360 = buffer.data(kh + 360);
    const auto *kh_361 = buffer.data(kh + 361);
    const auto *kh_362 = buffer.data(kh + 362);
    const auto *kh_363 = buffer.data(kh + 363);
    const auto *kh_364 = buffer.data(kh + 364);
    const auto *kh_365 = buffer.data(kh + 365);
    const auto *kh_366 = buffer.data(kh + 366);
    const auto *kh_367 = buffer.data(kh + 367);
    const auto *kh_368 = buffer.data(kh + 368);
    const auto *kh_369 = buffer.data(kh + 369);
    const auto *kh_370 = buffer.data(kh + 370);
    const auto *kh_371 = buffer.data(kh + 371);
    const auto *kh_372 = buffer.data(kh + 372);
    const auto *kh_373 = buffer.data(kh + 373);
    const auto *kh_374 = buffer.data(kh + 374);
    const auto *kh_375 = buffer.data(kh + 375);
    const auto *kh_376 = buffer.data(kh + 376);
    const auto *kh_377 = buffer.data(kh + 377);
    const auto *kh_378 = buffer.data(kh + 378);
    const auto *kh_379 = buffer.data(kh + 379);
    const auto *kh_380 = buffer.data(kh + 380);
    const auto *kh_381 = buffer.data(kh + 381);
    const auto *kh_382 = buffer.data(kh + 382);
    const auto *kh_383 = buffer.data(kh + 383);
    const auto *kh_384 = buffer.data(kh + 384);
    const auto *kh_385 = buffer.data(kh + 385);
    const auto *kh_386 = buffer.data(kh + 386);
    const auto *kh_387 = buffer.data(kh + 387);
    const auto *kh_388 = buffer.data(kh + 388);
    const auto *kh_389 = buffer.data(kh + 389);
    const auto *kh_390 = buffer.data(kh + 390);
    const auto *kh_391 = buffer.data(kh + 391);
    const auto *kh_392 = buffer.data(kh + 392);
    const auto *kh_393 = buffer.data(kh + 393);
    const auto *kh_394 = buffer.data(kh + 394);
    const auto *kh_395 = buffer.data(kh + 395);
    const auto *kh_396 = buffer.data(kh + 396);
    const auto *kh_397 = buffer.data(kh + 397);
    const auto *kh_398 = buffer.data(kh + 398);
    const auto *kh_399 = buffer.data(kh + 399);
    const auto *kh_400 = buffer.data(kh + 400);
    const auto *kh_401 = buffer.data(kh + 401);
    const auto *kh_402 = buffer.data(kh + 402);
    const auto *kh_403 = buffer.data(kh + 403);
    const auto *kh_404 = buffer.data(kh + 404);
    const auto *kh_405 = buffer.data(kh + 405);
    const auto *kh_406 = buffer.data(kh + 406);
    const auto *kh_407 = buffer.data(kh + 407);
    const auto *kh_408 = buffer.data(kh + 408);
    const auto *kh_409 = buffer.data(kh + 409);
    const auto *kh_410 = buffer.data(kh + 410);
    const auto *kh_411 = buffer.data(kh + 411);
    const auto *kh_412 = buffer.data(kh + 412);
    const auto *kh_413 = buffer.data(kh + 413);
    const auto *kh_414 = buffer.data(kh + 414);
    const auto *kh_415 = buffer.data(kh + 415);
    const auto *kh_416 = buffer.data(kh + 416);
    const auto *kh_417 = buffer.data(kh + 417);
    const auto *kh_418 = buffer.data(kh + 418);
    const auto *kh_419 = buffer.data(kh + 419);
    const auto *kh_420 = buffer.data(kh + 420);
    const auto *kh_421 = buffer.data(kh + 421);
    const auto *kh_422 = buffer.data(kh + 422);
    const auto *kh_423 = buffer.data(kh + 423);
    const auto *kh_424 = buffer.data(kh + 424);
    const auto *kh_425 = buffer.data(kh + 425);
    const auto *kh_426 = buffer.data(kh + 426);
    const auto *kh_427 = buffer.data(kh + 427);
    const auto *kh_428 = buffer.data(kh + 428);
    const auto *kh_429 = buffer.data(kh + 429);
    const auto *kh_430 = buffer.data(kh + 430);
    const auto *kh_431 = buffer.data(kh + 431);
    const auto *kh_432 = buffer.data(kh + 432);
    const auto *kh_433 = buffer.data(kh + 433);
    const auto *kh_434 = buffer.data(kh + 434);
    const auto *kh_435 = buffer.data(kh + 435);
    const auto *kh_436 = buffer.data(kh + 436);
    const auto *kh_437 = buffer.data(kh + 437);
    const auto *kh_438 = buffer.data(kh + 438);
    const auto *kh_439 = buffer.data(kh + 439);
    const auto *kh_440 = buffer.data(kh + 440);
    const auto *kh_441 = buffer.data(kh + 441);
    const auto *kh_442 = buffer.data(kh + 442);
    const auto *kh_443 = buffer.data(kh + 443);
    const auto *kh_444 = buffer.data(kh + 444);
    const auto *kh_445 = buffer.data(kh + 445);
    const auto *kh_446 = buffer.data(kh + 446);
    const auto *kh_447 = buffer.data(kh + 447);
    const auto *kh_448 = buffer.data(kh + 448);
    const auto *kh_449 = buffer.data(kh + 449);

    const auto *mh_300 = buffer.data(mh + 300);
    const auto *mh_301 = buffer.data(mh + 301);
    const auto *mh_302 = buffer.data(mh + 302);
    const auto *mh_303 = buffer.data(mh + 303);
    const auto *mh_304 = buffer.data(mh + 304);
    const auto *mh_305 = buffer.data(mh + 305);
    const auto *mh_306 = buffer.data(mh + 306);
    const auto *mh_307 = buffer.data(mh + 307);
    const auto *mh_308 = buffer.data(mh + 308);
    const auto *mh_309 = buffer.data(mh + 309);
    const auto *mh_310 = buffer.data(mh + 310);
    const auto *mh_311 = buffer.data(mh + 311);
    const auto *mh_312 = buffer.data(mh + 312);
    const auto *mh_313 = buffer.data(mh + 313);
    const auto *mh_314 = buffer.data(mh + 314);
    const auto *mh_315 = buffer.data(mh + 315);
    const auto *mh_316 = buffer.data(mh + 316);
    const auto *mh_317 = buffer.data(mh + 317);
    const auto *mh_318 = buffer.data(mh + 318);
    const auto *mh_319 = buffer.data(mh + 319);
    const auto *mh_320 = buffer.data(mh + 320);
    const auto *mh_321 = buffer.data(mh + 321);
    const auto *mh_322 = buffer.data(mh + 322);
    const auto *mh_323 = buffer.data(mh + 323);
    const auto *mh_324 = buffer.data(mh + 324);
    const auto *mh_325 = buffer.data(mh + 325);
    const auto *mh_326 = buffer.data(mh + 326);
    const auto *mh_327 = buffer.data(mh + 327);
    const auto *mh_328 = buffer.data(mh + 328);
    const auto *mh_329 = buffer.data(mh + 329);
    const auto *mh_330 = buffer.data(mh + 330);
    const auto *mh_331 = buffer.data(mh + 331);
    const auto *mh_332 = buffer.data(mh + 332);
    const auto *mh_333 = buffer.data(mh + 333);
    const auto *mh_334 = buffer.data(mh + 334);
    const auto *mh_335 = buffer.data(mh + 335);
    const auto *mh_336 = buffer.data(mh + 336);
    const auto *mh_337 = buffer.data(mh + 337);
    const auto *mh_338 = buffer.data(mh + 338);
    const auto *mh_339 = buffer.data(mh + 339);
    const auto *mh_340 = buffer.data(mh + 340);
    const auto *mh_341 = buffer.data(mh + 341);
    const auto *mh_342 = buffer.data(mh + 342);
    const auto *mh_343 = buffer.data(mh + 343);
    const auto *mh_344 = buffer.data(mh + 344);
    const auto *mh_345 = buffer.data(mh + 345);
    const auto *mh_346 = buffer.data(mh + 346);
    const auto *mh_347 = buffer.data(mh + 347);
    const auto *mh_348 = buffer.data(mh + 348);
    const auto *mh_349 = buffer.data(mh + 349);
    const auto *mh_350 = buffer.data(mh + 350);
    const auto *mh_351 = buffer.data(mh + 351);
    const auto *mh_352 = buffer.data(mh + 352);
    const auto *mh_353 = buffer.data(mh + 353);
    const auto *mh_354 = buffer.data(mh + 354);
    const auto *mh_355 = buffer.data(mh + 355);
    const auto *mh_356 = buffer.data(mh + 356);
    const auto *mh_357 = buffer.data(mh + 357);
    const auto *mh_358 = buffer.data(mh + 358);
    const auto *mh_359 = buffer.data(mh + 359);
    const auto *mh_360 = buffer.data(mh + 360);
    const auto *mh_361 = buffer.data(mh + 361);
    const auto *mh_362 = buffer.data(mh + 362);
    const auto *mh_363 = buffer.data(mh + 363);
    const auto *mh_364 = buffer.data(mh + 364);
    const auto *mh_365 = buffer.data(mh + 365);
    const auto *mh_366 = buffer.data(mh + 366);
    const auto *mh_367 = buffer.data(mh + 367);
    const auto *mh_368 = buffer.data(mh + 368);
    const auto *mh_369 = buffer.data(mh + 369);
    const auto *mh_370 = buffer.data(mh + 370);
    const auto *mh_371 = buffer.data(mh + 371);
    const auto *mh_372 = buffer.data(mh + 372);
    const auto *mh_373 = buffer.data(mh + 373);
    const auto *mh_374 = buffer.data(mh + 374);
    const auto *mh_375 = buffer.data(mh + 375);
    const auto *mh_376 = buffer.data(mh + 376);
    const auto *mh_377 = buffer.data(mh + 377);
    const auto *mh_378 = buffer.data(mh + 378);
    const auto *mh_379 = buffer.data(mh + 379);
    const auto *mh_380 = buffer.data(mh + 380);
    const auto *mh_381 = buffer.data(mh + 381);
    const auto *mh_382 = buffer.data(mh + 382);
    const auto *mh_383 = buffer.data(mh + 383);
    const auto *mh_384 = buffer.data(mh + 384);
    const auto *mh_385 = buffer.data(mh + 385);
    const auto *mh_386 = buffer.data(mh + 386);
    const auto *mh_387 = buffer.data(mh + 387);
    const auto *mh_388 = buffer.data(mh + 388);
    const auto *mh_389 = buffer.data(mh + 389);
    const auto *mh_390 = buffer.data(mh + 390);
    const auto *mh_391 = buffer.data(mh + 391);
    const auto *mh_392 = buffer.data(mh + 392);
    const auto *mh_393 = buffer.data(mh + 393);
    const auto *mh_394 = buffer.data(mh + 394);
    const auto *mh_395 = buffer.data(mh + 395);
    const auto *mh_396 = buffer.data(mh + 396);
    const auto *mh_397 = buffer.data(mh + 397);
    const auto *mh_398 = buffer.data(mh + 398);
    const auto *mh_399 = buffer.data(mh + 399);
    const auto *mh_400 = buffer.data(mh + 400);
    const auto *mh_401 = buffer.data(mh + 401);
    const auto *mh_402 = buffer.data(mh + 402);
    const auto *mh_403 = buffer.data(mh + 403);
    const auto *mh_404 = buffer.data(mh + 404);
    const auto *mh_405 = buffer.data(mh + 405);
    const auto *mh_406 = buffer.data(mh + 406);
    const auto *mh_407 = buffer.data(mh + 407);
    const auto *mh_408 = buffer.data(mh + 408);
    const auto *mh_409 = buffer.data(mh + 409);
    const auto *mh_410 = buffer.data(mh + 410);
    const auto *mh_411 = buffer.data(mh + 411);
    const auto *mh_412 = buffer.data(mh + 412);
    const auto *mh_413 = buffer.data(mh + 413);
    const auto *mh_414 = buffer.data(mh + 414);
    const auto *mh_415 = buffer.data(mh + 415);
    const auto *mh_416 = buffer.data(mh + 416);
    const auto *mh_417 = buffer.data(mh + 417);
    const auto *mh_418 = buffer.data(mh + 418);
    const auto *mh_419 = buffer.data(mh + 419);
    const auto *mh_420 = buffer.data(mh + 420);
    const auto *mh_421 = buffer.data(mh + 421);
    const auto *mh_422 = buffer.data(mh + 422);
    const auto *mh_423 = buffer.data(mh + 423);
    const auto *mh_424 = buffer.data(mh + 424);
    const auto *mh_425 = buffer.data(mh + 425);
    const auto *mh_426 = buffer.data(mh + 426);
    const auto *mh_427 = buffer.data(mh + 427);
    const auto *mh_428 = buffer.data(mh + 428);
    const auto *mh_429 = buffer.data(mh + 429);
    const auto *mh_430 = buffer.data(mh + 430);
    const auto *mh_431 = buffer.data(mh + 431);
    const auto *mh_432 = buffer.data(mh + 432);
    const auto *mh_433 = buffer.data(mh + 433);
    const auto *mh_434 = buffer.data(mh + 434);
    const auto *mh_435 = buffer.data(mh + 435);
    const auto *mh_436 = buffer.data(mh + 436);
    const auto *mh_437 = buffer.data(mh + 437);
    const auto *mh_438 = buffer.data(mh + 438);
    const auto *mh_439 = buffer.data(mh + 439);
    const auto *mh_440 = buffer.data(mh + 440);
    const auto *mh_441 = buffer.data(mh + 441);
    const auto *mh_442 = buffer.data(mh + 442);
    const auto *mh_443 = buffer.data(mh + 443);
    const auto *mh_444 = buffer.data(mh + 444);
    const auto *mh_445 = buffer.data(mh + 445);
    const auto *mh_446 = buffer.data(mh + 446);
    const auto *mh_447 = buffer.data(mh + 447);
    const auto *mh_448 = buffer.data(mh + 448);
    const auto *mh_449 = buffer.data(mh + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, kh_300, kh_301, kh_302, kh_303, \
                         kh_304, mh_300, mh_301, mh_302, mh_303, \
                         mh_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * kh_300[k]
                   + f_0 * mh_300[k];

        t_301[k] = -4.0 * kh_301[k]
                   + f_0 * mh_301[k];

        t_302[k] = -4.0 * kh_302[k]
                   + f_0 * mh_302[k];

        t_303[k] = -4.0 * kh_303[k]
                   + f_0 * mh_303[k];

        t_304[k] = -4.0 * kh_304[k]
                   + f_0 * mh_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, kh_305, kh_306, kh_307, kh_308, \
                         kh_309, mh_305, mh_306, mh_307, mh_308, \
                         mh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * kh_305[k]
                   + f_0 * mh_305[k];

        t_306[k] = -4.0 * kh_306[k]
                   + f_0 * mh_306[k];

        t_307[k] = -4.0 * kh_307[k]
                   + f_0 * mh_307[k];

        t_308[k] = -4.0 * kh_308[k]
                   + f_0 * mh_308[k];

        t_309[k] = -4.0 * kh_309[k]
                   + f_0 * mh_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, kh_310, kh_311, kh_312, kh_313, \
                         kh_314, mh_310, mh_311, mh_312, mh_313, \
                         mh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -4.0 * kh_310[k]
                   + f_0 * mh_310[k];

        t_311[k] = -4.0 * kh_311[k]
                   + f_0 * mh_311[k];

        t_312[k] = -4.0 * kh_312[k]
                   + f_0 * mh_312[k];

        t_313[k] = -4.0 * kh_313[k]
                   + f_0 * mh_313[k];

        t_314[k] = -4.0 * kh_314[k]
                   + f_0 * mh_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, kh_315, kh_316, kh_317, kh_318, \
                         kh_319, mh_315, mh_316, mh_317, mh_318, \
                         mh_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -3.0 * kh_315[k]
                   + f_0 * mh_315[k];

        t_316[k] = -3.0 * kh_316[k]
                   + f_0 * mh_316[k];

        t_317[k] = -3.0 * kh_317[k]
                   + f_0 * mh_317[k];

        t_318[k] = -3.0 * kh_318[k]
                   + f_0 * mh_318[k];

        t_319[k] = -3.0 * kh_319[k]
                   + f_0 * mh_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, kh_320, kh_321, kh_322, kh_323, \
                         kh_324, mh_320, mh_321, mh_322, mh_323, \
                         mh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -3.0 * kh_320[k]
                   + f_0 * mh_320[k];

        t_321[k] = -3.0 * kh_321[k]
                   + f_0 * mh_321[k];

        t_322[k] = -3.0 * kh_322[k]
                   + f_0 * mh_322[k];

        t_323[k] = -3.0 * kh_323[k]
                   + f_0 * mh_323[k];

        t_324[k] = -3.0 * kh_324[k]
                   + f_0 * mh_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, kh_325, kh_326, kh_327, kh_328, \
                         kh_329, mh_325, mh_326, mh_327, mh_328, \
                         mh_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -3.0 * kh_325[k]
                   + f_0 * mh_325[k];

        t_326[k] = -3.0 * kh_326[k]
                   + f_0 * mh_326[k];

        t_327[k] = -3.0 * kh_327[k]
                   + f_0 * mh_327[k];

        t_328[k] = -3.0 * kh_328[k]
                   + f_0 * mh_328[k];

        t_329[k] = -3.0 * kh_329[k]
                   + f_0 * mh_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, kh_330, kh_331, kh_332, kh_333, \
                         kh_334, mh_330, mh_331, mh_332, mh_333, \
                         mh_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -3.0 * kh_330[k]
                   + f_0 * mh_330[k];

        t_331[k] = -3.0 * kh_331[k]
                   + f_0 * mh_331[k];

        t_332[k] = -3.0 * kh_332[k]
                   + f_0 * mh_332[k];

        t_333[k] = -3.0 * kh_333[k]
                   + f_0 * mh_333[k];

        t_334[k] = -3.0 * kh_334[k]
                   + f_0 * mh_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, kh_335, kh_336, kh_337, kh_338, \
                         kh_339, mh_335, mh_336, mh_337, mh_338, \
                         mh_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -3.0 * kh_335[k]
                   + f_0 * mh_335[k];

        t_336[k] = -3.0 * kh_336[k]
                   + f_0 * mh_336[k];

        t_337[k] = -3.0 * kh_337[k]
                   + f_0 * mh_337[k];

        t_338[k] = -3.0 * kh_338[k]
                   + f_0 * mh_338[k];

        t_339[k] = -3.0 * kh_339[k]
                   + f_0 * mh_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, kh_340, kh_341, kh_342, kh_343, \
                         kh_344, mh_340, mh_341, mh_342, mh_343, \
                         mh_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -3.0 * kh_340[k]
                   + f_0 * mh_340[k];

        t_341[k] = -3.0 * kh_341[k]
                   + f_0 * mh_341[k];

        t_342[k] = -3.0 * kh_342[k]
                   + f_0 * mh_342[k];

        t_343[k] = -3.0 * kh_343[k]
                   + f_0 * mh_343[k];

        t_344[k] = -3.0 * kh_344[k]
                   + f_0 * mh_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, kh_345, kh_346, kh_347, kh_348, \
                         kh_349, mh_345, mh_346, mh_347, mh_348, \
                         mh_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -3.0 * kh_345[k]
                   + f_0 * mh_345[k];

        t_346[k] = -3.0 * kh_346[k]
                   + f_0 * mh_346[k];

        t_347[k] = -3.0 * kh_347[k]
                   + f_0 * mh_347[k];

        t_348[k] = -3.0 * kh_348[k]
                   + f_0 * mh_348[k];

        t_349[k] = -3.0 * kh_349[k]
                   + f_0 * mh_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, kh_350, kh_351, kh_352, kh_353, \
                         kh_354, mh_350, mh_351, mh_352, mh_353, \
                         mh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -3.0 * kh_350[k]
                   + f_0 * mh_350[k];

        t_351[k] = -3.0 * kh_351[k]
                   + f_0 * mh_351[k];

        t_352[k] = -3.0 * kh_352[k]
                   + f_0 * mh_352[k];

        t_353[k] = -3.0 * kh_353[k]
                   + f_0 * mh_353[k];

        t_354[k] = -3.0 * kh_354[k]
                   + f_0 * mh_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, kh_355, kh_356, kh_357, kh_358, \
                         kh_359, mh_355, mh_356, mh_357, mh_358, \
                         mh_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -3.0 * kh_355[k]
                   + f_0 * mh_355[k];

        t_356[k] = -3.0 * kh_356[k]
                   + f_0 * mh_356[k];

        t_357[k] = -3.0 * kh_357[k]
                   + f_0 * mh_357[k];

        t_358[k] = -3.0 * kh_358[k]
                   + f_0 * mh_358[k];

        t_359[k] = -3.0 * kh_359[k]
                   + f_0 * mh_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, kh_360, kh_361, kh_362, kh_363, \
                         kh_364, mh_360, mh_361, mh_362, mh_363, \
                         mh_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -3.0 * kh_360[k]
                   + f_0 * mh_360[k];

        t_361[k] = -3.0 * kh_361[k]
                   + f_0 * mh_361[k];

        t_362[k] = -3.0 * kh_362[k]
                   + f_0 * mh_362[k];

        t_363[k] = -3.0 * kh_363[k]
                   + f_0 * mh_363[k];

        t_364[k] = -3.0 * kh_364[k]
                   + f_0 * mh_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, kh_365, kh_366, kh_367, kh_368, \
                         kh_369, mh_365, mh_366, mh_367, mh_368, \
                         mh_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -3.0 * kh_365[k]
                   + f_0 * mh_365[k];

        t_366[k] = -3.0 * kh_366[k]
                   + f_0 * mh_366[k];

        t_367[k] = -3.0 * kh_367[k]
                   + f_0 * mh_367[k];

        t_368[k] = -3.0 * kh_368[k]
                   + f_0 * mh_368[k];

        t_369[k] = -3.0 * kh_369[k]
                   + f_0 * mh_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, kh_370, kh_371, kh_372, kh_373, \
                         kh_374, mh_370, mh_371, mh_372, mh_373, \
                         mh_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -3.0 * kh_370[k]
                   + f_0 * mh_370[k];

        t_371[k] = -3.0 * kh_371[k]
                   + f_0 * mh_371[k];

        t_372[k] = -3.0 * kh_372[k]
                   + f_0 * mh_372[k];

        t_373[k] = -3.0 * kh_373[k]
                   + f_0 * mh_373[k];

        t_374[k] = -3.0 * kh_374[k]
                   + f_0 * mh_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, kh_375, kh_376, kh_377, kh_378, \
                         kh_379, mh_375, mh_376, mh_377, mh_378, \
                         mh_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -3.0 * kh_375[k]
                   + f_0 * mh_375[k];

        t_376[k] = -3.0 * kh_376[k]
                   + f_0 * mh_376[k];

        t_377[k] = -3.0 * kh_377[k]
                   + f_0 * mh_377[k];

        t_378[k] = -3.0 * kh_378[k]
                   + f_0 * mh_378[k];

        t_379[k] = -3.0 * kh_379[k]
                   + f_0 * mh_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, kh_380, kh_381, kh_382, kh_383, \
                         kh_384, mh_380, mh_381, mh_382, mh_383, \
                         mh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -3.0 * kh_380[k]
                   + f_0 * mh_380[k];

        t_381[k] = -3.0 * kh_381[k]
                   + f_0 * mh_381[k];

        t_382[k] = -3.0 * kh_382[k]
                   + f_0 * mh_382[k];

        t_383[k] = -3.0 * kh_383[k]
                   + f_0 * mh_383[k];

        t_384[k] = -3.0 * kh_384[k]
                   + f_0 * mh_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, kh_385, kh_386, kh_387, kh_388, \
                         kh_389, mh_385, mh_386, mh_387, mh_388, \
                         mh_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -3.0 * kh_385[k]
                   + f_0 * mh_385[k];

        t_386[k] = -3.0 * kh_386[k]
                   + f_0 * mh_386[k];

        t_387[k] = -3.0 * kh_387[k]
                   + f_0 * mh_387[k];

        t_388[k] = -3.0 * kh_388[k]
                   + f_0 * mh_388[k];

        t_389[k] = -3.0 * kh_389[k]
                   + f_0 * mh_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, kh_390, kh_391, kh_392, kh_393, \
                         kh_394, mh_390, mh_391, mh_392, mh_393, \
                         mh_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -3.0 * kh_390[k]
                   + f_0 * mh_390[k];

        t_391[k] = -3.0 * kh_391[k]
                   + f_0 * mh_391[k];

        t_392[k] = -3.0 * kh_392[k]
                   + f_0 * mh_392[k];

        t_393[k] = -3.0 * kh_393[k]
                   + f_0 * mh_393[k];

        t_394[k] = -3.0 * kh_394[k]
                   + f_0 * mh_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, kh_395, kh_396, kh_397, kh_398, \
                         kh_399, mh_395, mh_396, mh_397, mh_398, \
                         mh_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -3.0 * kh_395[k]
                   + f_0 * mh_395[k];

        t_396[k] = -3.0 * kh_396[k]
                   + f_0 * mh_396[k];

        t_397[k] = -3.0 * kh_397[k]
                   + f_0 * mh_397[k];

        t_398[k] = -3.0 * kh_398[k]
                   + f_0 * mh_398[k];

        t_399[k] = -3.0 * kh_399[k]
                   + f_0 * mh_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, kh_400, kh_401, kh_402, kh_403, \
                         kh_404, mh_400, mh_401, mh_402, mh_403, \
                         mh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -3.0 * kh_400[k]
                   + f_0 * mh_400[k];

        t_401[k] = -3.0 * kh_401[k]
                   + f_0 * mh_401[k];

        t_402[k] = -3.0 * kh_402[k]
                   + f_0 * mh_402[k];

        t_403[k] = -3.0 * kh_403[k]
                   + f_0 * mh_403[k];

        t_404[k] = -3.0 * kh_404[k]
                   + f_0 * mh_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, kh_405, kh_406, kh_407, kh_408, \
                         kh_409, mh_405, mh_406, mh_407, mh_408, \
                         mh_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -3.0 * kh_405[k]
                   + f_0 * mh_405[k];

        t_406[k] = -3.0 * kh_406[k]
                   + f_0 * mh_406[k];

        t_407[k] = -3.0 * kh_407[k]
                   + f_0 * mh_407[k];

        t_408[k] = -3.0 * kh_408[k]
                   + f_0 * mh_408[k];

        t_409[k] = -3.0 * kh_409[k]
                   + f_0 * mh_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, kh_410, kh_411, kh_412, kh_413, \
                         kh_414, mh_410, mh_411, mh_412, mh_413, \
                         mh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -3.0 * kh_410[k]
                   + f_0 * mh_410[k];

        t_411[k] = -3.0 * kh_411[k]
                   + f_0 * mh_411[k];

        t_412[k] = -3.0 * kh_412[k]
                   + f_0 * mh_412[k];

        t_413[k] = -3.0 * kh_413[k]
                   + f_0 * mh_413[k];

        t_414[k] = -3.0 * kh_414[k]
                   + f_0 * mh_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, kh_415, kh_416, kh_417, kh_418, \
                         kh_419, mh_415, mh_416, mh_417, mh_418, \
                         mh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -3.0 * kh_415[k]
                   + f_0 * mh_415[k];

        t_416[k] = -3.0 * kh_416[k]
                   + f_0 * mh_416[k];

        t_417[k] = -3.0 * kh_417[k]
                   + f_0 * mh_417[k];

        t_418[k] = -3.0 * kh_418[k]
                   + f_0 * mh_418[k];

        t_419[k] = -3.0 * kh_419[k]
                   + f_0 * mh_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, kh_420, kh_421, kh_422, kh_423, \
                         kh_424, mh_420, mh_421, mh_422, mh_423, \
                         mh_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -3.0 * kh_420[k]
                   + f_0 * mh_420[k];

        t_421[k] = -3.0 * kh_421[k]
                   + f_0 * mh_421[k];

        t_422[k] = -3.0 * kh_422[k]
                   + f_0 * mh_422[k];

        t_423[k] = -3.0 * kh_423[k]
                   + f_0 * mh_423[k];

        t_424[k] = -3.0 * kh_424[k]
                   + f_0 * mh_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, kh_425, kh_426, kh_427, kh_428, \
                         kh_429, mh_425, mh_426, mh_427, mh_428, \
                         mh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -3.0 * kh_425[k]
                   + f_0 * mh_425[k];

        t_426[k] = -3.0 * kh_426[k]
                   + f_0 * mh_426[k];

        t_427[k] = -3.0 * kh_427[k]
                   + f_0 * mh_427[k];

        t_428[k] = -3.0 * kh_428[k]
                   + f_0 * mh_428[k];

        t_429[k] = -3.0 * kh_429[k]
                   + f_0 * mh_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, kh_430, kh_431, kh_432, kh_433, \
                         kh_434, mh_430, mh_431, mh_432, mh_433, \
                         mh_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -3.0 * kh_430[k]
                   + f_0 * mh_430[k];

        t_431[k] = -3.0 * kh_431[k]
                   + f_0 * mh_431[k];

        t_432[k] = -3.0 * kh_432[k]
                   + f_0 * mh_432[k];

        t_433[k] = -3.0 * kh_433[k]
                   + f_0 * mh_433[k];

        t_434[k] = -3.0 * kh_434[k]
                   + f_0 * mh_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, kh_435, kh_436, kh_437, kh_438, \
                         kh_439, mh_435, mh_436, mh_437, mh_438, \
                         mh_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -3.0 * kh_435[k]
                   + f_0 * mh_435[k];

        t_436[k] = -3.0 * kh_436[k]
                   + f_0 * mh_436[k];

        t_437[k] = -3.0 * kh_437[k]
                   + f_0 * mh_437[k];

        t_438[k] = -3.0 * kh_438[k]
                   + f_0 * mh_438[k];

        t_439[k] = -3.0 * kh_439[k]
                   + f_0 * mh_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, kh_440, kh_441, kh_442, kh_443, \
                         kh_444, mh_440, mh_441, mh_442, mh_443, \
                         mh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -3.0 * kh_440[k]
                   + f_0 * mh_440[k];

        t_441[k] = -2.0 * kh_441[k]
                   + f_0 * mh_441[k];

        t_442[k] = -2.0 * kh_442[k]
                   + f_0 * mh_442[k];

        t_443[k] = -2.0 * kh_443[k]
                   + f_0 * mh_443[k];

        t_444[k] = -2.0 * kh_444[k]
                   + f_0 * mh_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, kh_445, kh_446, kh_447, kh_448, \
                         kh_449, mh_445, mh_446, mh_447, mh_448, \
                         mh_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -2.0 * kh_445[k]
                   + f_0 * mh_445[k];

        t_446[k] = -2.0 * kh_446[k]
                   + f_0 * mh_446[k];

        t_447[k] = -2.0 * kh_447[k]
                   + f_0 * mh_447[k];

        t_448[k] = -2.0 * kh_448[k]
                   + f_0 * mh_448[k];

        t_449[k] = -2.0 * kh_449[k]
                   + f_0 * mh_449[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_450 = buffer.data(target + 450);
    auto *t_451 = buffer.data(target + 451);
    auto *t_452 = buffer.data(target + 452);
    auto *t_453 = buffer.data(target + 453);
    auto *t_454 = buffer.data(target + 454);
    auto *t_455 = buffer.data(target + 455);
    auto *t_456 = buffer.data(target + 456);
    auto *t_457 = buffer.data(target + 457);
    auto *t_458 = buffer.data(target + 458);
    auto *t_459 = buffer.data(target + 459);
    auto *t_460 = buffer.data(target + 460);
    auto *t_461 = buffer.data(target + 461);
    auto *t_462 = buffer.data(target + 462);
    auto *t_463 = buffer.data(target + 463);
    auto *t_464 = buffer.data(target + 464);
    auto *t_465 = buffer.data(target + 465);
    auto *t_466 = buffer.data(target + 466);
    auto *t_467 = buffer.data(target + 467);
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);
    auto *t_474 = buffer.data(target + 474);
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);
    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);
    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);
    auto *t_500 = buffer.data(target + 500);
    auto *t_501 = buffer.data(target + 501);
    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);
    auto *t_504 = buffer.data(target + 504);
    auto *t_505 = buffer.data(target + 505);
    auto *t_506 = buffer.data(target + 506);
    auto *t_507 = buffer.data(target + 507);
    auto *t_508 = buffer.data(target + 508);
    auto *t_509 = buffer.data(target + 509);
    auto *t_510 = buffer.data(target + 510);
    auto *t_511 = buffer.data(target + 511);
    auto *t_512 = buffer.data(target + 512);
    auto *t_513 = buffer.data(target + 513);
    auto *t_514 = buffer.data(target + 514);
    auto *t_515 = buffer.data(target + 515);
    auto *t_516 = buffer.data(target + 516);
    auto *t_517 = buffer.data(target + 517);
    auto *t_518 = buffer.data(target + 518);
    auto *t_519 = buffer.data(target + 519);
    auto *t_520 = buffer.data(target + 520);
    auto *t_521 = buffer.data(target + 521);
    auto *t_522 = buffer.data(target + 522);
    auto *t_523 = buffer.data(target + 523);
    auto *t_524 = buffer.data(target + 524);
    auto *t_525 = buffer.data(target + 525);
    auto *t_526 = buffer.data(target + 526);
    auto *t_527 = buffer.data(target + 527);
    auto *t_528 = buffer.data(target + 528);
    auto *t_529 = buffer.data(target + 529);
    auto *t_530 = buffer.data(target + 530);
    auto *t_531 = buffer.data(target + 531);
    auto *t_532 = buffer.data(target + 532);
    auto *t_533 = buffer.data(target + 533);
    auto *t_534 = buffer.data(target + 534);
    auto *t_535 = buffer.data(target + 535);
    auto *t_536 = buffer.data(target + 536);
    auto *t_537 = buffer.data(target + 537);
    auto *t_538 = buffer.data(target + 538);
    auto *t_539 = buffer.data(target + 539);
    auto *t_540 = buffer.data(target + 540);
    auto *t_541 = buffer.data(target + 541);
    auto *t_542 = buffer.data(target + 542);
    auto *t_543 = buffer.data(target + 543);
    auto *t_544 = buffer.data(target + 544);
    auto *t_545 = buffer.data(target + 545);
    auto *t_546 = buffer.data(target + 546);
    auto *t_547 = buffer.data(target + 547);
    auto *t_548 = buffer.data(target + 548);
    auto *t_549 = buffer.data(target + 549);
    auto *t_550 = buffer.data(target + 550);
    auto *t_551 = buffer.data(target + 551);
    auto *t_552 = buffer.data(target + 552);
    auto *t_553 = buffer.data(target + 553);
    auto *t_554 = buffer.data(target + 554);
    auto *t_555 = buffer.data(target + 555);
    auto *t_556 = buffer.data(target + 556);
    auto *t_557 = buffer.data(target + 557);
    auto *t_558 = buffer.data(target + 558);
    auto *t_559 = buffer.data(target + 559);
    auto *t_560 = buffer.data(target + 560);
    auto *t_561 = buffer.data(target + 561);
    auto *t_562 = buffer.data(target + 562);
    auto *t_563 = buffer.data(target + 563);
    auto *t_564 = buffer.data(target + 564);
    auto *t_565 = buffer.data(target + 565);
    auto *t_566 = buffer.data(target + 566);
    auto *t_567 = buffer.data(target + 567);
    auto *t_568 = buffer.data(target + 568);
    auto *t_569 = buffer.data(target + 569);
    auto *t_570 = buffer.data(target + 570);
    auto *t_571 = buffer.data(target + 571);
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);
    auto *t_575 = buffer.data(target + 575);
    auto *t_576 = buffer.data(target + 576);
    auto *t_577 = buffer.data(target + 577);
    auto *t_578 = buffer.data(target + 578);
    auto *t_579 = buffer.data(target + 579);
    auto *t_580 = buffer.data(target + 580);
    auto *t_581 = buffer.data(target + 581);
    auto *t_582 = buffer.data(target + 582);
    auto *t_583 = buffer.data(target + 583);
    auto *t_584 = buffer.data(target + 584);
    auto *t_585 = buffer.data(target + 585);
    auto *t_586 = buffer.data(target + 586);
    auto *t_587 = buffer.data(target + 587);
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);

    const auto *kh_450 = buffer.data(kh + 450);
    const auto *kh_451 = buffer.data(kh + 451);
    const auto *kh_452 = buffer.data(kh + 452);
    const auto *kh_453 = buffer.data(kh + 453);
    const auto *kh_454 = buffer.data(kh + 454);
    const auto *kh_455 = buffer.data(kh + 455);
    const auto *kh_456 = buffer.data(kh + 456);
    const auto *kh_457 = buffer.data(kh + 457);
    const auto *kh_458 = buffer.data(kh + 458);
    const auto *kh_459 = buffer.data(kh + 459);
    const auto *kh_460 = buffer.data(kh + 460);
    const auto *kh_461 = buffer.data(kh + 461);
    const auto *kh_462 = buffer.data(kh + 462);
    const auto *kh_463 = buffer.data(kh + 463);
    const auto *kh_464 = buffer.data(kh + 464);
    const auto *kh_465 = buffer.data(kh + 465);
    const auto *kh_466 = buffer.data(kh + 466);
    const auto *kh_467 = buffer.data(kh + 467);
    const auto *kh_468 = buffer.data(kh + 468);
    const auto *kh_469 = buffer.data(kh + 469);
    const auto *kh_470 = buffer.data(kh + 470);
    const auto *kh_471 = buffer.data(kh + 471);
    const auto *kh_472 = buffer.data(kh + 472);
    const auto *kh_473 = buffer.data(kh + 473);
    const auto *kh_474 = buffer.data(kh + 474);
    const auto *kh_475 = buffer.data(kh + 475);
    const auto *kh_476 = buffer.data(kh + 476);
    const auto *kh_477 = buffer.data(kh + 477);
    const auto *kh_478 = buffer.data(kh + 478);
    const auto *kh_479 = buffer.data(kh + 479);
    const auto *kh_480 = buffer.data(kh + 480);
    const auto *kh_481 = buffer.data(kh + 481);
    const auto *kh_482 = buffer.data(kh + 482);
    const auto *kh_483 = buffer.data(kh + 483);
    const auto *kh_484 = buffer.data(kh + 484);
    const auto *kh_485 = buffer.data(kh + 485);
    const auto *kh_486 = buffer.data(kh + 486);
    const auto *kh_487 = buffer.data(kh + 487);
    const auto *kh_488 = buffer.data(kh + 488);
    const auto *kh_489 = buffer.data(kh + 489);
    const auto *kh_490 = buffer.data(kh + 490);
    const auto *kh_491 = buffer.data(kh + 491);
    const auto *kh_492 = buffer.data(kh + 492);
    const auto *kh_493 = buffer.data(kh + 493);
    const auto *kh_494 = buffer.data(kh + 494);
    const auto *kh_495 = buffer.data(kh + 495);
    const auto *kh_496 = buffer.data(kh + 496);
    const auto *kh_497 = buffer.data(kh + 497);
    const auto *kh_498 = buffer.data(kh + 498);
    const auto *kh_499 = buffer.data(kh + 499);
    const auto *kh_500 = buffer.data(kh + 500);
    const auto *kh_501 = buffer.data(kh + 501);
    const auto *kh_502 = buffer.data(kh + 502);
    const auto *kh_503 = buffer.data(kh + 503);
    const auto *kh_504 = buffer.data(kh + 504);
    const auto *kh_505 = buffer.data(kh + 505);
    const auto *kh_506 = buffer.data(kh + 506);
    const auto *kh_507 = buffer.data(kh + 507);
    const auto *kh_508 = buffer.data(kh + 508);
    const auto *kh_509 = buffer.data(kh + 509);
    const auto *kh_510 = buffer.data(kh + 510);
    const auto *kh_511 = buffer.data(kh + 511);
    const auto *kh_512 = buffer.data(kh + 512);
    const auto *kh_513 = buffer.data(kh + 513);
    const auto *kh_514 = buffer.data(kh + 514);
    const auto *kh_515 = buffer.data(kh + 515);
    const auto *kh_516 = buffer.data(kh + 516);
    const auto *kh_517 = buffer.data(kh + 517);
    const auto *kh_518 = buffer.data(kh + 518);
    const auto *kh_519 = buffer.data(kh + 519);
    const auto *kh_520 = buffer.data(kh + 520);
    const auto *kh_521 = buffer.data(kh + 521);
    const auto *kh_522 = buffer.data(kh + 522);
    const auto *kh_523 = buffer.data(kh + 523);
    const auto *kh_524 = buffer.data(kh + 524);
    const auto *kh_525 = buffer.data(kh + 525);
    const auto *kh_526 = buffer.data(kh + 526);
    const auto *kh_527 = buffer.data(kh + 527);
    const auto *kh_528 = buffer.data(kh + 528);
    const auto *kh_529 = buffer.data(kh + 529);
    const auto *kh_530 = buffer.data(kh + 530);
    const auto *kh_531 = buffer.data(kh + 531);
    const auto *kh_532 = buffer.data(kh + 532);
    const auto *kh_533 = buffer.data(kh + 533);
    const auto *kh_534 = buffer.data(kh + 534);
    const auto *kh_535 = buffer.data(kh + 535);
    const auto *kh_536 = buffer.data(kh + 536);
    const auto *kh_537 = buffer.data(kh + 537);
    const auto *kh_538 = buffer.data(kh + 538);
    const auto *kh_539 = buffer.data(kh + 539);
    const auto *kh_540 = buffer.data(kh + 540);
    const auto *kh_541 = buffer.data(kh + 541);
    const auto *kh_542 = buffer.data(kh + 542);
    const auto *kh_543 = buffer.data(kh + 543);
    const auto *kh_544 = buffer.data(kh + 544);
    const auto *kh_545 = buffer.data(kh + 545);
    const auto *kh_546 = buffer.data(kh + 546);
    const auto *kh_547 = buffer.data(kh + 547);
    const auto *kh_548 = buffer.data(kh + 548);
    const auto *kh_549 = buffer.data(kh + 549);
    const auto *kh_550 = buffer.data(kh + 550);
    const auto *kh_551 = buffer.data(kh + 551);
    const auto *kh_552 = buffer.data(kh + 552);
    const auto *kh_553 = buffer.data(kh + 553);
    const auto *kh_554 = buffer.data(kh + 554);
    const auto *kh_555 = buffer.data(kh + 555);
    const auto *kh_556 = buffer.data(kh + 556);
    const auto *kh_557 = buffer.data(kh + 557);
    const auto *kh_558 = buffer.data(kh + 558);
    const auto *kh_559 = buffer.data(kh + 559);
    const auto *kh_560 = buffer.data(kh + 560);
    const auto *kh_561 = buffer.data(kh + 561);
    const auto *kh_562 = buffer.data(kh + 562);
    const auto *kh_563 = buffer.data(kh + 563);
    const auto *kh_564 = buffer.data(kh + 564);
    const auto *kh_565 = buffer.data(kh + 565);
    const auto *kh_566 = buffer.data(kh + 566);
    const auto *kh_567 = buffer.data(kh + 567);
    const auto *kh_568 = buffer.data(kh + 568);
    const auto *kh_569 = buffer.data(kh + 569);
    const auto *kh_570 = buffer.data(kh + 570);
    const auto *kh_571 = buffer.data(kh + 571);
    const auto *kh_572 = buffer.data(kh + 572);
    const auto *kh_573 = buffer.data(kh + 573);
    const auto *kh_574 = buffer.data(kh + 574);
    const auto *kh_575 = buffer.data(kh + 575);
    const auto *kh_576 = buffer.data(kh + 576);
    const auto *kh_577 = buffer.data(kh + 577);
    const auto *kh_578 = buffer.data(kh + 578);
    const auto *kh_579 = buffer.data(kh + 579);
    const auto *kh_580 = buffer.data(kh + 580);
    const auto *kh_581 = buffer.data(kh + 581);
    const auto *kh_582 = buffer.data(kh + 582);
    const auto *kh_583 = buffer.data(kh + 583);
    const auto *kh_584 = buffer.data(kh + 584);
    const auto *kh_585 = buffer.data(kh + 585);
    const auto *kh_586 = buffer.data(kh + 586);
    const auto *kh_587 = buffer.data(kh + 587);
    const auto *kh_588 = buffer.data(kh + 588);
    const auto *kh_589 = buffer.data(kh + 589);
    const auto *kh_590 = buffer.data(kh + 590);
    const auto *kh_591 = buffer.data(kh + 591);
    const auto *kh_592 = buffer.data(kh + 592);
    const auto *kh_593 = buffer.data(kh + 593);
    const auto *kh_594 = buffer.data(kh + 594);
    const auto *kh_595 = buffer.data(kh + 595);
    const auto *kh_596 = buffer.data(kh + 596);
    const auto *kh_597 = buffer.data(kh + 597);
    const auto *kh_598 = buffer.data(kh + 598);
    const auto *kh_599 = buffer.data(kh + 599);

    const auto *mh_450 = buffer.data(mh + 450);
    const auto *mh_451 = buffer.data(mh + 451);
    const auto *mh_452 = buffer.data(mh + 452);
    const auto *mh_453 = buffer.data(mh + 453);
    const auto *mh_454 = buffer.data(mh + 454);
    const auto *mh_455 = buffer.data(mh + 455);
    const auto *mh_456 = buffer.data(mh + 456);
    const auto *mh_457 = buffer.data(mh + 457);
    const auto *mh_458 = buffer.data(mh + 458);
    const auto *mh_459 = buffer.data(mh + 459);
    const auto *mh_460 = buffer.data(mh + 460);
    const auto *mh_461 = buffer.data(mh + 461);
    const auto *mh_462 = buffer.data(mh + 462);
    const auto *mh_463 = buffer.data(mh + 463);
    const auto *mh_464 = buffer.data(mh + 464);
    const auto *mh_465 = buffer.data(mh + 465);
    const auto *mh_466 = buffer.data(mh + 466);
    const auto *mh_467 = buffer.data(mh + 467);
    const auto *mh_468 = buffer.data(mh + 468);
    const auto *mh_469 = buffer.data(mh + 469);
    const auto *mh_470 = buffer.data(mh + 470);
    const auto *mh_471 = buffer.data(mh + 471);
    const auto *mh_472 = buffer.data(mh + 472);
    const auto *mh_473 = buffer.data(mh + 473);
    const auto *mh_474 = buffer.data(mh + 474);
    const auto *mh_475 = buffer.data(mh + 475);
    const auto *mh_476 = buffer.data(mh + 476);
    const auto *mh_477 = buffer.data(mh + 477);
    const auto *mh_478 = buffer.data(mh + 478);
    const auto *mh_479 = buffer.data(mh + 479);
    const auto *mh_480 = buffer.data(mh + 480);
    const auto *mh_481 = buffer.data(mh + 481);
    const auto *mh_482 = buffer.data(mh + 482);
    const auto *mh_483 = buffer.data(mh + 483);
    const auto *mh_484 = buffer.data(mh + 484);
    const auto *mh_485 = buffer.data(mh + 485);
    const auto *mh_486 = buffer.data(mh + 486);
    const auto *mh_487 = buffer.data(mh + 487);
    const auto *mh_488 = buffer.data(mh + 488);
    const auto *mh_489 = buffer.data(mh + 489);
    const auto *mh_490 = buffer.data(mh + 490);
    const auto *mh_491 = buffer.data(mh + 491);
    const auto *mh_492 = buffer.data(mh + 492);
    const auto *mh_493 = buffer.data(mh + 493);
    const auto *mh_494 = buffer.data(mh + 494);
    const auto *mh_495 = buffer.data(mh + 495);
    const auto *mh_496 = buffer.data(mh + 496);
    const auto *mh_497 = buffer.data(mh + 497);
    const auto *mh_498 = buffer.data(mh + 498);
    const auto *mh_499 = buffer.data(mh + 499);
    const auto *mh_500 = buffer.data(mh + 500);
    const auto *mh_501 = buffer.data(mh + 501);
    const auto *mh_502 = buffer.data(mh + 502);
    const auto *mh_503 = buffer.data(mh + 503);
    const auto *mh_504 = buffer.data(mh + 504);
    const auto *mh_505 = buffer.data(mh + 505);
    const auto *mh_506 = buffer.data(mh + 506);
    const auto *mh_507 = buffer.data(mh + 507);
    const auto *mh_508 = buffer.data(mh + 508);
    const auto *mh_509 = buffer.data(mh + 509);
    const auto *mh_510 = buffer.data(mh + 510);
    const auto *mh_511 = buffer.data(mh + 511);
    const auto *mh_512 = buffer.data(mh + 512);
    const auto *mh_513 = buffer.data(mh + 513);
    const auto *mh_514 = buffer.data(mh + 514);
    const auto *mh_515 = buffer.data(mh + 515);
    const auto *mh_516 = buffer.data(mh + 516);
    const auto *mh_517 = buffer.data(mh + 517);
    const auto *mh_518 = buffer.data(mh + 518);
    const auto *mh_519 = buffer.data(mh + 519);
    const auto *mh_520 = buffer.data(mh + 520);
    const auto *mh_521 = buffer.data(mh + 521);
    const auto *mh_522 = buffer.data(mh + 522);
    const auto *mh_523 = buffer.data(mh + 523);
    const auto *mh_524 = buffer.data(mh + 524);
    const auto *mh_525 = buffer.data(mh + 525);
    const auto *mh_526 = buffer.data(mh + 526);
    const auto *mh_527 = buffer.data(mh + 527);
    const auto *mh_528 = buffer.data(mh + 528);
    const auto *mh_529 = buffer.data(mh + 529);
    const auto *mh_530 = buffer.data(mh + 530);
    const auto *mh_531 = buffer.data(mh + 531);
    const auto *mh_532 = buffer.data(mh + 532);
    const auto *mh_533 = buffer.data(mh + 533);
    const auto *mh_534 = buffer.data(mh + 534);
    const auto *mh_535 = buffer.data(mh + 535);
    const auto *mh_536 = buffer.data(mh + 536);
    const auto *mh_537 = buffer.data(mh + 537);
    const auto *mh_538 = buffer.data(mh + 538);
    const auto *mh_539 = buffer.data(mh + 539);
    const auto *mh_540 = buffer.data(mh + 540);
    const auto *mh_541 = buffer.data(mh + 541);
    const auto *mh_542 = buffer.data(mh + 542);
    const auto *mh_543 = buffer.data(mh + 543);
    const auto *mh_544 = buffer.data(mh + 544);
    const auto *mh_545 = buffer.data(mh + 545);
    const auto *mh_546 = buffer.data(mh + 546);
    const auto *mh_547 = buffer.data(mh + 547);
    const auto *mh_548 = buffer.data(mh + 548);
    const auto *mh_549 = buffer.data(mh + 549);
    const auto *mh_550 = buffer.data(mh + 550);
    const auto *mh_551 = buffer.data(mh + 551);
    const auto *mh_552 = buffer.data(mh + 552);
    const auto *mh_553 = buffer.data(mh + 553);
    const auto *mh_554 = buffer.data(mh + 554);
    const auto *mh_555 = buffer.data(mh + 555);
    const auto *mh_556 = buffer.data(mh + 556);
    const auto *mh_557 = buffer.data(mh + 557);
    const auto *mh_558 = buffer.data(mh + 558);
    const auto *mh_559 = buffer.data(mh + 559);
    const auto *mh_560 = buffer.data(mh + 560);
    const auto *mh_561 = buffer.data(mh + 561);
    const auto *mh_562 = buffer.data(mh + 562);
    const auto *mh_563 = buffer.data(mh + 563);
    const auto *mh_564 = buffer.data(mh + 564);
    const auto *mh_565 = buffer.data(mh + 565);
    const auto *mh_566 = buffer.data(mh + 566);
    const auto *mh_567 = buffer.data(mh + 567);
    const auto *mh_568 = buffer.data(mh + 568);
    const auto *mh_569 = buffer.data(mh + 569);
    const auto *mh_570 = buffer.data(mh + 570);
    const auto *mh_571 = buffer.data(mh + 571);
    const auto *mh_572 = buffer.data(mh + 572);
    const auto *mh_573 = buffer.data(mh + 573);
    const auto *mh_574 = buffer.data(mh + 574);
    const auto *mh_575 = buffer.data(mh + 575);
    const auto *mh_576 = buffer.data(mh + 576);
    const auto *mh_577 = buffer.data(mh + 577);
    const auto *mh_578 = buffer.data(mh + 578);
    const auto *mh_579 = buffer.data(mh + 579);
    const auto *mh_580 = buffer.data(mh + 580);
    const auto *mh_581 = buffer.data(mh + 581);
    const auto *mh_582 = buffer.data(mh + 582);
    const auto *mh_583 = buffer.data(mh + 583);
    const auto *mh_584 = buffer.data(mh + 584);
    const auto *mh_585 = buffer.data(mh + 585);
    const auto *mh_586 = buffer.data(mh + 586);
    const auto *mh_587 = buffer.data(mh + 587);
    const auto *mh_588 = buffer.data(mh + 588);
    const auto *mh_589 = buffer.data(mh + 589);
    const auto *mh_590 = buffer.data(mh + 590);
    const auto *mh_591 = buffer.data(mh + 591);
    const auto *mh_592 = buffer.data(mh + 592);
    const auto *mh_593 = buffer.data(mh + 593);
    const auto *mh_594 = buffer.data(mh + 594);
    const auto *mh_595 = buffer.data(mh + 595);
    const auto *mh_596 = buffer.data(mh + 596);
    const auto *mh_597 = buffer.data(mh + 597);
    const auto *mh_598 = buffer.data(mh + 598);
    const auto *mh_599 = buffer.data(mh + 599);

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, kh_450, kh_451, kh_452, kh_453, \
                         kh_454, mh_450, mh_451, mh_452, mh_453, \
                         mh_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -2.0 * kh_450[k]
                   + f_0 * mh_450[k];

        t_451[k] = -2.0 * kh_451[k]
                   + f_0 * mh_451[k];

        t_452[k] = -2.0 * kh_452[k]
                   + f_0 * mh_452[k];

        t_453[k] = -2.0 * kh_453[k]
                   + f_0 * mh_453[k];

        t_454[k] = -2.0 * kh_454[k]
                   + f_0 * mh_454[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, kh_455, kh_456, kh_457, kh_458, \
                         kh_459, mh_455, mh_456, mh_457, mh_458, \
                         mh_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -2.0 * kh_455[k]
                   + f_0 * mh_455[k];

        t_456[k] = -2.0 * kh_456[k]
                   + f_0 * mh_456[k];

        t_457[k] = -2.0 * kh_457[k]
                   + f_0 * mh_457[k];

        t_458[k] = -2.0 * kh_458[k]
                   + f_0 * mh_458[k];

        t_459[k] = -2.0 * kh_459[k]
                   + f_0 * mh_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, kh_460, kh_461, kh_462, kh_463, \
                         kh_464, mh_460, mh_461, mh_462, mh_463, \
                         mh_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -2.0 * kh_460[k]
                   + f_0 * mh_460[k];

        t_461[k] = -2.0 * kh_461[k]
                   + f_0 * mh_461[k];

        t_462[k] = -2.0 * kh_462[k]
                   + f_0 * mh_462[k];

        t_463[k] = -2.0 * kh_463[k]
                   + f_0 * mh_463[k];

        t_464[k] = -2.0 * kh_464[k]
                   + f_0 * mh_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, kh_465, kh_466, kh_467, kh_468, \
                         kh_469, mh_465, mh_466, mh_467, mh_468, \
                         mh_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -2.0 * kh_465[k]
                   + f_0 * mh_465[k];

        t_466[k] = -2.0 * kh_466[k]
                   + f_0 * mh_466[k];

        t_467[k] = -2.0 * kh_467[k]
                   + f_0 * mh_467[k];

        t_468[k] = -2.0 * kh_468[k]
                   + f_0 * mh_468[k];

        t_469[k] = -2.0 * kh_469[k]
                   + f_0 * mh_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, kh_470, kh_471, kh_472, kh_473, \
                         kh_474, mh_470, mh_471, mh_472, mh_473, \
                         mh_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -2.0 * kh_470[k]
                   + f_0 * mh_470[k];

        t_471[k] = -2.0 * kh_471[k]
                   + f_0 * mh_471[k];

        t_472[k] = -2.0 * kh_472[k]
                   + f_0 * mh_472[k];

        t_473[k] = -2.0 * kh_473[k]
                   + f_0 * mh_473[k];

        t_474[k] = -2.0 * kh_474[k]
                   + f_0 * mh_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, kh_475, kh_476, kh_477, kh_478, \
                         kh_479, mh_475, mh_476, mh_477, mh_478, \
                         mh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -2.0 * kh_475[k]
                   + f_0 * mh_475[k];

        t_476[k] = -2.0 * kh_476[k]
                   + f_0 * mh_476[k];

        t_477[k] = -2.0 * kh_477[k]
                   + f_0 * mh_477[k];

        t_478[k] = -2.0 * kh_478[k]
                   + f_0 * mh_478[k];

        t_479[k] = -2.0 * kh_479[k]
                   + f_0 * mh_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, kh_480, kh_481, kh_482, kh_483, \
                         kh_484, mh_480, mh_481, mh_482, mh_483, \
                         mh_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -2.0 * kh_480[k]
                   + f_0 * mh_480[k];

        t_481[k] = -2.0 * kh_481[k]
                   + f_0 * mh_481[k];

        t_482[k] = -2.0 * kh_482[k]
                   + f_0 * mh_482[k];

        t_483[k] = -2.0 * kh_483[k]
                   + f_0 * mh_483[k];

        t_484[k] = -2.0 * kh_484[k]
                   + f_0 * mh_484[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, kh_485, kh_486, kh_487, kh_488, \
                         kh_489, mh_485, mh_486, mh_487, mh_488, \
                         mh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -2.0 * kh_485[k]
                   + f_0 * mh_485[k];

        t_486[k] = -2.0 * kh_486[k]
                   + f_0 * mh_486[k];

        t_487[k] = -2.0 * kh_487[k]
                   + f_0 * mh_487[k];

        t_488[k] = -2.0 * kh_488[k]
                   + f_0 * mh_488[k];

        t_489[k] = -2.0 * kh_489[k]
                   + f_0 * mh_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, kh_490, kh_491, kh_492, kh_493, \
                         kh_494, mh_490, mh_491, mh_492, mh_493, \
                         mh_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -2.0 * kh_490[k]
                   + f_0 * mh_490[k];

        t_491[k] = -2.0 * kh_491[k]
                   + f_0 * mh_491[k];

        t_492[k] = -2.0 * kh_492[k]
                   + f_0 * mh_492[k];

        t_493[k] = -2.0 * kh_493[k]
                   + f_0 * mh_493[k];

        t_494[k] = -2.0 * kh_494[k]
                   + f_0 * mh_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, kh_495, kh_496, kh_497, kh_498, \
                         kh_499, mh_495, mh_496, mh_497, mh_498, \
                         mh_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -2.0 * kh_495[k]
                   + f_0 * mh_495[k];

        t_496[k] = -2.0 * kh_496[k]
                   + f_0 * mh_496[k];

        t_497[k] = -2.0 * kh_497[k]
                   + f_0 * mh_497[k];

        t_498[k] = -2.0 * kh_498[k]
                   + f_0 * mh_498[k];

        t_499[k] = -2.0 * kh_499[k]
                   + f_0 * mh_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, kh_500, kh_501, kh_502, kh_503, \
                         kh_504, mh_500, mh_501, mh_502, mh_503, \
                         mh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -2.0 * kh_500[k]
                   + f_0 * mh_500[k];

        t_501[k] = -2.0 * kh_501[k]
                   + f_0 * mh_501[k];

        t_502[k] = -2.0 * kh_502[k]
                   + f_0 * mh_502[k];

        t_503[k] = -2.0 * kh_503[k]
                   + f_0 * mh_503[k];

        t_504[k] = -2.0 * kh_504[k]
                   + f_0 * mh_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, kh_505, kh_506, kh_507, kh_508, \
                         kh_509, mh_505, mh_506, mh_507, mh_508, \
                         mh_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -2.0 * kh_505[k]
                   + f_0 * mh_505[k];

        t_506[k] = -2.0 * kh_506[k]
                   + f_0 * mh_506[k];

        t_507[k] = -2.0 * kh_507[k]
                   + f_0 * mh_507[k];

        t_508[k] = -2.0 * kh_508[k]
                   + f_0 * mh_508[k];

        t_509[k] = -2.0 * kh_509[k]
                   + f_0 * mh_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, kh_510, kh_511, kh_512, kh_513, \
                         kh_514, mh_510, mh_511, mh_512, mh_513, \
                         mh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -2.0 * kh_510[k]
                   + f_0 * mh_510[k];

        t_511[k] = -2.0 * kh_511[k]
                   + f_0 * mh_511[k];

        t_512[k] = -2.0 * kh_512[k]
                   + f_0 * mh_512[k];

        t_513[k] = -2.0 * kh_513[k]
                   + f_0 * mh_513[k];

        t_514[k] = -2.0 * kh_514[k]
                   + f_0 * mh_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, kh_515, kh_516, kh_517, kh_518, \
                         kh_519, mh_515, mh_516, mh_517, mh_518, \
                         mh_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -2.0 * kh_515[k]
                   + f_0 * mh_515[k];

        t_516[k] = -2.0 * kh_516[k]
                   + f_0 * mh_516[k];

        t_517[k] = -2.0 * kh_517[k]
                   + f_0 * mh_517[k];

        t_518[k] = -2.0 * kh_518[k]
                   + f_0 * mh_518[k];

        t_519[k] = -2.0 * kh_519[k]
                   + f_0 * mh_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, kh_520, kh_521, kh_522, kh_523, \
                         kh_524, mh_520, mh_521, mh_522, mh_523, \
                         mh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -2.0 * kh_520[k]
                   + f_0 * mh_520[k];

        t_521[k] = -2.0 * kh_521[k]
                   + f_0 * mh_521[k];

        t_522[k] = -2.0 * kh_522[k]
                   + f_0 * mh_522[k];

        t_523[k] = -2.0 * kh_523[k]
                   + f_0 * mh_523[k];

        t_524[k] = -2.0 * kh_524[k]
                   + f_0 * mh_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, kh_525, kh_526, kh_527, kh_528, \
                         kh_529, mh_525, mh_526, mh_527, mh_528, \
                         mh_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -2.0 * kh_525[k]
                   + f_0 * mh_525[k];

        t_526[k] = -2.0 * kh_526[k]
                   + f_0 * mh_526[k];

        t_527[k] = -2.0 * kh_527[k]
                   + f_0 * mh_527[k];

        t_528[k] = -2.0 * kh_528[k]
                   + f_0 * mh_528[k];

        t_529[k] = -2.0 * kh_529[k]
                   + f_0 * mh_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, kh_530, kh_531, kh_532, kh_533, \
                         kh_534, mh_530, mh_531, mh_532, mh_533, \
                         mh_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -2.0 * kh_530[k]
                   + f_0 * mh_530[k];

        t_531[k] = -2.0 * kh_531[k]
                   + f_0 * mh_531[k];

        t_532[k] = -2.0 * kh_532[k]
                   + f_0 * mh_532[k];

        t_533[k] = -2.0 * kh_533[k]
                   + f_0 * mh_533[k];

        t_534[k] = -2.0 * kh_534[k]
                   + f_0 * mh_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, kh_535, kh_536, kh_537, kh_538, \
                         kh_539, mh_535, mh_536, mh_537, mh_538, \
                         mh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -2.0 * kh_535[k]
                   + f_0 * mh_535[k];

        t_536[k] = -2.0 * kh_536[k]
                   + f_0 * mh_536[k];

        t_537[k] = -2.0 * kh_537[k]
                   + f_0 * mh_537[k];

        t_538[k] = -2.0 * kh_538[k]
                   + f_0 * mh_538[k];

        t_539[k] = -2.0 * kh_539[k]
                   + f_0 * mh_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, kh_540, kh_541, kh_542, kh_543, \
                         kh_544, mh_540, mh_541, mh_542, mh_543, \
                         mh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -2.0 * kh_540[k]
                   + f_0 * mh_540[k];

        t_541[k] = -2.0 * kh_541[k]
                   + f_0 * mh_541[k];

        t_542[k] = -2.0 * kh_542[k]
                   + f_0 * mh_542[k];

        t_543[k] = -2.0 * kh_543[k]
                   + f_0 * mh_543[k];

        t_544[k] = -2.0 * kh_544[k]
                   + f_0 * mh_544[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, kh_545, kh_546, kh_547, kh_548, \
                         kh_549, mh_545, mh_546, mh_547, mh_548, \
                         mh_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -2.0 * kh_545[k]
                   + f_0 * mh_545[k];

        t_546[k] = -2.0 * kh_546[k]
                   + f_0 * mh_546[k];

        t_547[k] = -2.0 * kh_547[k]
                   + f_0 * mh_547[k];

        t_548[k] = -2.0 * kh_548[k]
                   + f_0 * mh_548[k];

        t_549[k] = -2.0 * kh_549[k]
                   + f_0 * mh_549[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, kh_550, kh_551, kh_552, kh_553, \
                         kh_554, mh_550, mh_551, mh_552, mh_553, \
                         mh_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -2.0 * kh_550[k]
                   + f_0 * mh_550[k];

        t_551[k] = -2.0 * kh_551[k]
                   + f_0 * mh_551[k];

        t_552[k] = -2.0 * kh_552[k]
                   + f_0 * mh_552[k];

        t_553[k] = -2.0 * kh_553[k]
                   + f_0 * mh_553[k];

        t_554[k] = -2.0 * kh_554[k]
                   + f_0 * mh_554[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, kh_555, kh_556, kh_557, kh_558, \
                         kh_559, mh_555, mh_556, mh_557, mh_558, \
                         mh_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -2.0 * kh_555[k]
                   + f_0 * mh_555[k];

        t_556[k] = -2.0 * kh_556[k]
                   + f_0 * mh_556[k];

        t_557[k] = -2.0 * kh_557[k]
                   + f_0 * mh_557[k];

        t_558[k] = -2.0 * kh_558[k]
                   + f_0 * mh_558[k];

        t_559[k] = -2.0 * kh_559[k]
                   + f_0 * mh_559[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, kh_560, kh_561, kh_562, kh_563, \
                         kh_564, mh_560, mh_561, mh_562, mh_563, \
                         mh_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = -2.0 * kh_560[k]
                   + f_0 * mh_560[k];

        t_561[k] = -2.0 * kh_561[k]
                   + f_0 * mh_561[k];

        t_562[k] = -2.0 * kh_562[k]
                   + f_0 * mh_562[k];

        t_563[k] = -2.0 * kh_563[k]
                   + f_0 * mh_563[k];

        t_564[k] = -2.0 * kh_564[k]
                   + f_0 * mh_564[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, kh_565, kh_566, kh_567, kh_568, \
                         kh_569, mh_565, mh_566, mh_567, mh_568, \
                         mh_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = -2.0 * kh_565[k]
                   + f_0 * mh_565[k];

        t_566[k] = -2.0 * kh_566[k]
                   + f_0 * mh_566[k];

        t_567[k] = -2.0 * kh_567[k]
                   + f_0 * mh_567[k];

        t_568[k] = -2.0 * kh_568[k]
                   + f_0 * mh_568[k];

        t_569[k] = -2.0 * kh_569[k]
                   + f_0 * mh_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, kh_570, kh_571, kh_572, kh_573, \
                         kh_574, mh_570, mh_571, mh_572, mh_573, \
                         mh_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -2.0 * kh_570[k]
                   + f_0 * mh_570[k];

        t_571[k] = -2.0 * kh_571[k]
                   + f_0 * mh_571[k];

        t_572[k] = -2.0 * kh_572[k]
                   + f_0 * mh_572[k];

        t_573[k] = -2.0 * kh_573[k]
                   + f_0 * mh_573[k];

        t_574[k] = -2.0 * kh_574[k]
                   + f_0 * mh_574[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, kh_575, kh_576, kh_577, kh_578, \
                         kh_579, mh_575, mh_576, mh_577, mh_578, \
                         mh_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = -2.0 * kh_575[k]
                   + f_0 * mh_575[k];

        t_576[k] = -2.0 * kh_576[k]
                   + f_0 * mh_576[k];

        t_577[k] = -2.0 * kh_577[k]
                   + f_0 * mh_577[k];

        t_578[k] = -2.0 * kh_578[k]
                   + f_0 * mh_578[k];

        t_579[k] = -2.0 * kh_579[k]
                   + f_0 * mh_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, kh_580, kh_581, kh_582, kh_583, \
                         kh_584, mh_580, mh_581, mh_582, mh_583, \
                         mh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -2.0 * kh_580[k]
                   + f_0 * mh_580[k];

        t_581[k] = -2.0 * kh_581[k]
                   + f_0 * mh_581[k];

        t_582[k] = -2.0 * kh_582[k]
                   + f_0 * mh_582[k];

        t_583[k] = -2.0 * kh_583[k]
                   + f_0 * mh_583[k];

        t_584[k] = -2.0 * kh_584[k]
                   + f_0 * mh_584[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, kh_585, kh_586, kh_587, kh_588, \
                         kh_589, mh_585, mh_586, mh_587, mh_588, \
                         mh_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -2.0 * kh_585[k]
                   + f_0 * mh_585[k];

        t_586[k] = -2.0 * kh_586[k]
                   + f_0 * mh_586[k];

        t_587[k] = -2.0 * kh_587[k]
                   + f_0 * mh_587[k];

        t_588[k] = -kh_588[k]
                   + f_0 * mh_588[k];

        t_589[k] = -kh_589[k]
                   + f_0 * mh_589[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, kh_590, kh_591, kh_592, kh_593, \
                         kh_594, mh_590, mh_591, mh_592, mh_593, \
                         mh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -kh_590[k]
                   + f_0 * mh_590[k];

        t_591[k] = -kh_591[k]
                   + f_0 * mh_591[k];

        t_592[k] = -kh_592[k]
                   + f_0 * mh_592[k];

        t_593[k] = -kh_593[k]
                   + f_0 * mh_593[k];

        t_594[k] = -kh_594[k]
                   + f_0 * mh_594[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, kh_595, kh_596, kh_597, kh_598, \
                         kh_599, mh_595, mh_596, mh_597, mh_598, \
                         mh_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -kh_595[k]
                   + f_0 * mh_595[k];

        t_596[k] = -kh_596[k]
                   + f_0 * mh_596[k];

        t_597[k] = -kh_597[k]
                   + f_0 * mh_597[k];

        t_598[k] = -kh_598[k]
                   + f_0 * mh_598[k];

        t_599[k] = -kh_599[k]
                   + f_0 * mh_599[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);
    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);
    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);
    auto *t_618 = buffer.data(target + 618);
    auto *t_619 = buffer.data(target + 619);
    auto *t_620 = buffer.data(target + 620);
    auto *t_621 = buffer.data(target + 621);
    auto *t_622 = buffer.data(target + 622);
    auto *t_623 = buffer.data(target + 623);
    auto *t_624 = buffer.data(target + 624);
    auto *t_625 = buffer.data(target + 625);
    auto *t_626 = buffer.data(target + 626);
    auto *t_627 = buffer.data(target + 627);
    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);
    auto *t_637 = buffer.data(target + 637);
    auto *t_638 = buffer.data(target + 638);
    auto *t_639 = buffer.data(target + 639);
    auto *t_640 = buffer.data(target + 640);
    auto *t_641 = buffer.data(target + 641);
    auto *t_642 = buffer.data(target + 642);
    auto *t_643 = buffer.data(target + 643);
    auto *t_644 = buffer.data(target + 644);
    auto *t_645 = buffer.data(target + 645);
    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);
    auto *t_654 = buffer.data(target + 654);
    auto *t_655 = buffer.data(target + 655);
    auto *t_656 = buffer.data(target + 656);
    auto *t_657 = buffer.data(target + 657);
    auto *t_658 = buffer.data(target + 658);
    auto *t_659 = buffer.data(target + 659);
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);
    auto *t_663 = buffer.data(target + 663);
    auto *t_664 = buffer.data(target + 664);
    auto *t_665 = buffer.data(target + 665);
    auto *t_666 = buffer.data(target + 666);
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);
    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);
    auto *t_686 = buffer.data(target + 686);
    auto *t_687 = buffer.data(target + 687);
    auto *t_688 = buffer.data(target + 688);
    auto *t_689 = buffer.data(target + 689);
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);
    auto *t_693 = buffer.data(target + 693);
    auto *t_694 = buffer.data(target + 694);
    auto *t_695 = buffer.data(target + 695);
    auto *t_696 = buffer.data(target + 696);
    auto *t_697 = buffer.data(target + 697);
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);
    auto *t_708 = buffer.data(target + 708);
    auto *t_709 = buffer.data(target + 709);
    auto *t_710 = buffer.data(target + 710);
    auto *t_711 = buffer.data(target + 711);
    auto *t_712 = buffer.data(target + 712);
    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
    auto *t_722 = buffer.data(target + 722);
    auto *t_723 = buffer.data(target + 723);
    auto *t_724 = buffer.data(target + 724);
    auto *t_725 = buffer.data(target + 725);
    auto *t_726 = buffer.data(target + 726);
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);
    auto *t_732 = buffer.data(target + 732);
    auto *t_733 = buffer.data(target + 733);
    auto *t_734 = buffer.data(target + 734);
    auto *t_735 = buffer.data(target + 735);
    auto *t_736 = buffer.data(target + 736);
    auto *t_737 = buffer.data(target + 737);
    auto *t_738 = buffer.data(target + 738);
    auto *t_739 = buffer.data(target + 739);
    auto *t_740 = buffer.data(target + 740);
    auto *t_741 = buffer.data(target + 741);
    auto *t_742 = buffer.data(target + 742);
    auto *t_743 = buffer.data(target + 743);
    auto *t_744 = buffer.data(target + 744);
    auto *t_745 = buffer.data(target + 745);
    auto *t_746 = buffer.data(target + 746);
    auto *t_747 = buffer.data(target + 747);
    auto *t_748 = buffer.data(target + 748);
    auto *t_749 = buffer.data(target + 749);

    const auto *kh_600 = buffer.data(kh + 600);
    const auto *kh_601 = buffer.data(kh + 601);
    const auto *kh_602 = buffer.data(kh + 602);
    const auto *kh_603 = buffer.data(kh + 603);
    const auto *kh_604 = buffer.data(kh + 604);
    const auto *kh_605 = buffer.data(kh + 605);
    const auto *kh_606 = buffer.data(kh + 606);
    const auto *kh_607 = buffer.data(kh + 607);
    const auto *kh_608 = buffer.data(kh + 608);
    const auto *kh_609 = buffer.data(kh + 609);
    const auto *kh_610 = buffer.data(kh + 610);
    const auto *kh_611 = buffer.data(kh + 611);
    const auto *kh_612 = buffer.data(kh + 612);
    const auto *kh_613 = buffer.data(kh + 613);
    const auto *kh_614 = buffer.data(kh + 614);
    const auto *kh_615 = buffer.data(kh + 615);
    const auto *kh_616 = buffer.data(kh + 616);
    const auto *kh_617 = buffer.data(kh + 617);
    const auto *kh_618 = buffer.data(kh + 618);
    const auto *kh_619 = buffer.data(kh + 619);
    const auto *kh_620 = buffer.data(kh + 620);
    const auto *kh_621 = buffer.data(kh + 621);
    const auto *kh_622 = buffer.data(kh + 622);
    const auto *kh_623 = buffer.data(kh + 623);
    const auto *kh_624 = buffer.data(kh + 624);
    const auto *kh_625 = buffer.data(kh + 625);
    const auto *kh_626 = buffer.data(kh + 626);
    const auto *kh_627 = buffer.data(kh + 627);
    const auto *kh_628 = buffer.data(kh + 628);
    const auto *kh_629 = buffer.data(kh + 629);
    const auto *kh_630 = buffer.data(kh + 630);
    const auto *kh_631 = buffer.data(kh + 631);
    const auto *kh_632 = buffer.data(kh + 632);
    const auto *kh_633 = buffer.data(kh + 633);
    const auto *kh_634 = buffer.data(kh + 634);
    const auto *kh_635 = buffer.data(kh + 635);
    const auto *kh_636 = buffer.data(kh + 636);
    const auto *kh_637 = buffer.data(kh + 637);
    const auto *kh_638 = buffer.data(kh + 638);
    const auto *kh_639 = buffer.data(kh + 639);
    const auto *kh_640 = buffer.data(kh + 640);
    const auto *kh_641 = buffer.data(kh + 641);
    const auto *kh_642 = buffer.data(kh + 642);
    const auto *kh_643 = buffer.data(kh + 643);
    const auto *kh_644 = buffer.data(kh + 644);
    const auto *kh_645 = buffer.data(kh + 645);
    const auto *kh_646 = buffer.data(kh + 646);
    const auto *kh_647 = buffer.data(kh + 647);
    const auto *kh_648 = buffer.data(kh + 648);
    const auto *kh_649 = buffer.data(kh + 649);
    const auto *kh_650 = buffer.data(kh + 650);
    const auto *kh_651 = buffer.data(kh + 651);
    const auto *kh_652 = buffer.data(kh + 652);
    const auto *kh_653 = buffer.data(kh + 653);
    const auto *kh_654 = buffer.data(kh + 654);
    const auto *kh_655 = buffer.data(kh + 655);
    const auto *kh_656 = buffer.data(kh + 656);
    const auto *kh_657 = buffer.data(kh + 657);
    const auto *kh_658 = buffer.data(kh + 658);
    const auto *kh_659 = buffer.data(kh + 659);
    const auto *kh_660 = buffer.data(kh + 660);
    const auto *kh_661 = buffer.data(kh + 661);
    const auto *kh_662 = buffer.data(kh + 662);
    const auto *kh_663 = buffer.data(kh + 663);
    const auto *kh_664 = buffer.data(kh + 664);
    const auto *kh_665 = buffer.data(kh + 665);
    const auto *kh_666 = buffer.data(kh + 666);
    const auto *kh_667 = buffer.data(kh + 667);
    const auto *kh_668 = buffer.data(kh + 668);
    const auto *kh_669 = buffer.data(kh + 669);
    const auto *kh_670 = buffer.data(kh + 670);
    const auto *kh_671 = buffer.data(kh + 671);
    const auto *kh_672 = buffer.data(kh + 672);
    const auto *kh_673 = buffer.data(kh + 673);
    const auto *kh_674 = buffer.data(kh + 674);
    const auto *kh_675 = buffer.data(kh + 675);
    const auto *kh_676 = buffer.data(kh + 676);
    const auto *kh_677 = buffer.data(kh + 677);
    const auto *kh_678 = buffer.data(kh + 678);
    const auto *kh_679 = buffer.data(kh + 679);
    const auto *kh_680 = buffer.data(kh + 680);
    const auto *kh_681 = buffer.data(kh + 681);
    const auto *kh_682 = buffer.data(kh + 682);
    const auto *kh_683 = buffer.data(kh + 683);
    const auto *kh_684 = buffer.data(kh + 684);
    const auto *kh_685 = buffer.data(kh + 685);
    const auto *kh_686 = buffer.data(kh + 686);
    const auto *kh_687 = buffer.data(kh + 687);
    const auto *kh_688 = buffer.data(kh + 688);
    const auto *kh_689 = buffer.data(kh + 689);
    const auto *kh_690 = buffer.data(kh + 690);
    const auto *kh_691 = buffer.data(kh + 691);
    const auto *kh_692 = buffer.data(kh + 692);
    const auto *kh_693 = buffer.data(kh + 693);
    const auto *kh_694 = buffer.data(kh + 694);
    const auto *kh_695 = buffer.data(kh + 695);
    const auto *kh_696 = buffer.data(kh + 696);
    const auto *kh_697 = buffer.data(kh + 697);
    const auto *kh_698 = buffer.data(kh + 698);
    const auto *kh_699 = buffer.data(kh + 699);
    const auto *kh_700 = buffer.data(kh + 700);
    const auto *kh_701 = buffer.data(kh + 701);
    const auto *kh_702 = buffer.data(kh + 702);
    const auto *kh_703 = buffer.data(kh + 703);
    const auto *kh_704 = buffer.data(kh + 704);
    const auto *kh_705 = buffer.data(kh + 705);
    const auto *kh_706 = buffer.data(kh + 706);
    const auto *kh_707 = buffer.data(kh + 707);
    const auto *kh_708 = buffer.data(kh + 708);
    const auto *kh_709 = buffer.data(kh + 709);
    const auto *kh_710 = buffer.data(kh + 710);
    const auto *kh_711 = buffer.data(kh + 711);
    const auto *kh_712 = buffer.data(kh + 712);
    const auto *kh_713 = buffer.data(kh + 713);
    const auto *kh_714 = buffer.data(kh + 714);
    const auto *kh_715 = buffer.data(kh + 715);
    const auto *kh_716 = buffer.data(kh + 716);
    const auto *kh_717 = buffer.data(kh + 717);
    const auto *kh_718 = buffer.data(kh + 718);
    const auto *kh_719 = buffer.data(kh + 719);
    const auto *kh_720 = buffer.data(kh + 720);
    const auto *kh_721 = buffer.data(kh + 721);
    const auto *kh_722 = buffer.data(kh + 722);
    const auto *kh_723 = buffer.data(kh + 723);
    const auto *kh_724 = buffer.data(kh + 724);
    const auto *kh_725 = buffer.data(kh + 725);
    const auto *kh_726 = buffer.data(kh + 726);
    const auto *kh_727 = buffer.data(kh + 727);
    const auto *kh_728 = buffer.data(kh + 728);
    const auto *kh_729 = buffer.data(kh + 729);
    const auto *kh_730 = buffer.data(kh + 730);
    const auto *kh_731 = buffer.data(kh + 731);
    const auto *kh_732 = buffer.data(kh + 732);
    const auto *kh_733 = buffer.data(kh + 733);
    const auto *kh_734 = buffer.data(kh + 734);
    const auto *kh_735 = buffer.data(kh + 735);
    const auto *kh_736 = buffer.data(kh + 736);
    const auto *kh_737 = buffer.data(kh + 737);
    const auto *kh_738 = buffer.data(kh + 738);
    const auto *kh_739 = buffer.data(kh + 739);
    const auto *kh_740 = buffer.data(kh + 740);
    const auto *kh_741 = buffer.data(kh + 741);
    const auto *kh_742 = buffer.data(kh + 742);
    const auto *kh_743 = buffer.data(kh + 743);
    const auto *kh_744 = buffer.data(kh + 744);
    const auto *kh_745 = buffer.data(kh + 745);
    const auto *kh_746 = buffer.data(kh + 746);
    const auto *kh_747 = buffer.data(kh + 747);
    const auto *kh_748 = buffer.data(kh + 748);
    const auto *kh_749 = buffer.data(kh + 749);

    const auto *mh_600 = buffer.data(mh + 600);
    const auto *mh_601 = buffer.data(mh + 601);
    const auto *mh_602 = buffer.data(mh + 602);
    const auto *mh_603 = buffer.data(mh + 603);
    const auto *mh_604 = buffer.data(mh + 604);
    const auto *mh_605 = buffer.data(mh + 605);
    const auto *mh_606 = buffer.data(mh + 606);
    const auto *mh_607 = buffer.data(mh + 607);
    const auto *mh_608 = buffer.data(mh + 608);
    const auto *mh_609 = buffer.data(mh + 609);
    const auto *mh_610 = buffer.data(mh + 610);
    const auto *mh_611 = buffer.data(mh + 611);
    const auto *mh_612 = buffer.data(mh + 612);
    const auto *mh_613 = buffer.data(mh + 613);
    const auto *mh_614 = buffer.data(mh + 614);
    const auto *mh_615 = buffer.data(mh + 615);
    const auto *mh_616 = buffer.data(mh + 616);
    const auto *mh_617 = buffer.data(mh + 617);
    const auto *mh_618 = buffer.data(mh + 618);
    const auto *mh_619 = buffer.data(mh + 619);
    const auto *mh_620 = buffer.data(mh + 620);
    const auto *mh_621 = buffer.data(mh + 621);
    const auto *mh_622 = buffer.data(mh + 622);
    const auto *mh_623 = buffer.data(mh + 623);
    const auto *mh_624 = buffer.data(mh + 624);
    const auto *mh_625 = buffer.data(mh + 625);
    const auto *mh_626 = buffer.data(mh + 626);
    const auto *mh_627 = buffer.data(mh + 627);
    const auto *mh_628 = buffer.data(mh + 628);
    const auto *mh_629 = buffer.data(mh + 629);
    const auto *mh_630 = buffer.data(mh + 630);
    const auto *mh_631 = buffer.data(mh + 631);
    const auto *mh_632 = buffer.data(mh + 632);
    const auto *mh_633 = buffer.data(mh + 633);
    const auto *mh_634 = buffer.data(mh + 634);
    const auto *mh_635 = buffer.data(mh + 635);
    const auto *mh_636 = buffer.data(mh + 636);
    const auto *mh_637 = buffer.data(mh + 637);
    const auto *mh_638 = buffer.data(mh + 638);
    const auto *mh_639 = buffer.data(mh + 639);
    const auto *mh_640 = buffer.data(mh + 640);
    const auto *mh_641 = buffer.data(mh + 641);
    const auto *mh_642 = buffer.data(mh + 642);
    const auto *mh_643 = buffer.data(mh + 643);
    const auto *mh_644 = buffer.data(mh + 644);
    const auto *mh_645 = buffer.data(mh + 645);
    const auto *mh_646 = buffer.data(mh + 646);
    const auto *mh_647 = buffer.data(mh + 647);
    const auto *mh_648 = buffer.data(mh + 648);
    const auto *mh_649 = buffer.data(mh + 649);
    const auto *mh_650 = buffer.data(mh + 650);
    const auto *mh_651 = buffer.data(mh + 651);
    const auto *mh_652 = buffer.data(mh + 652);
    const auto *mh_653 = buffer.data(mh + 653);
    const auto *mh_654 = buffer.data(mh + 654);
    const auto *mh_655 = buffer.data(mh + 655);
    const auto *mh_656 = buffer.data(mh + 656);
    const auto *mh_657 = buffer.data(mh + 657);
    const auto *mh_658 = buffer.data(mh + 658);
    const auto *mh_659 = buffer.data(mh + 659);
    const auto *mh_660 = buffer.data(mh + 660);
    const auto *mh_661 = buffer.data(mh + 661);
    const auto *mh_662 = buffer.data(mh + 662);
    const auto *mh_663 = buffer.data(mh + 663);
    const auto *mh_664 = buffer.data(mh + 664);
    const auto *mh_665 = buffer.data(mh + 665);
    const auto *mh_666 = buffer.data(mh + 666);
    const auto *mh_667 = buffer.data(mh + 667);
    const auto *mh_668 = buffer.data(mh + 668);
    const auto *mh_669 = buffer.data(mh + 669);
    const auto *mh_670 = buffer.data(mh + 670);
    const auto *mh_671 = buffer.data(mh + 671);
    const auto *mh_672 = buffer.data(mh + 672);
    const auto *mh_673 = buffer.data(mh + 673);
    const auto *mh_674 = buffer.data(mh + 674);
    const auto *mh_675 = buffer.data(mh + 675);
    const auto *mh_676 = buffer.data(mh + 676);
    const auto *mh_677 = buffer.data(mh + 677);
    const auto *mh_678 = buffer.data(mh + 678);
    const auto *mh_679 = buffer.data(mh + 679);
    const auto *mh_680 = buffer.data(mh + 680);
    const auto *mh_681 = buffer.data(mh + 681);
    const auto *mh_682 = buffer.data(mh + 682);
    const auto *mh_683 = buffer.data(mh + 683);
    const auto *mh_684 = buffer.data(mh + 684);
    const auto *mh_685 = buffer.data(mh + 685);
    const auto *mh_686 = buffer.data(mh + 686);
    const auto *mh_687 = buffer.data(mh + 687);
    const auto *mh_688 = buffer.data(mh + 688);
    const auto *mh_689 = buffer.data(mh + 689);
    const auto *mh_690 = buffer.data(mh + 690);
    const auto *mh_691 = buffer.data(mh + 691);
    const auto *mh_692 = buffer.data(mh + 692);
    const auto *mh_693 = buffer.data(mh + 693);
    const auto *mh_694 = buffer.data(mh + 694);
    const auto *mh_695 = buffer.data(mh + 695);
    const auto *mh_696 = buffer.data(mh + 696);
    const auto *mh_697 = buffer.data(mh + 697);
    const auto *mh_698 = buffer.data(mh + 698);
    const auto *mh_699 = buffer.data(mh + 699);
    const auto *mh_700 = buffer.data(mh + 700);
    const auto *mh_701 = buffer.data(mh + 701);
    const auto *mh_702 = buffer.data(mh + 702);
    const auto *mh_703 = buffer.data(mh + 703);
    const auto *mh_704 = buffer.data(mh + 704);
    const auto *mh_705 = buffer.data(mh + 705);
    const auto *mh_706 = buffer.data(mh + 706);
    const auto *mh_707 = buffer.data(mh + 707);
    const auto *mh_708 = buffer.data(mh + 708);
    const auto *mh_709 = buffer.data(mh + 709);
    const auto *mh_710 = buffer.data(mh + 710);
    const auto *mh_711 = buffer.data(mh + 711);
    const auto *mh_712 = buffer.data(mh + 712);
    const auto *mh_713 = buffer.data(mh + 713);
    const auto *mh_714 = buffer.data(mh + 714);
    const auto *mh_715 = buffer.data(mh + 715);
    const auto *mh_716 = buffer.data(mh + 716);
    const auto *mh_717 = buffer.data(mh + 717);
    const auto *mh_718 = buffer.data(mh + 718);
    const auto *mh_719 = buffer.data(mh + 719);
    const auto *mh_720 = buffer.data(mh + 720);
    const auto *mh_721 = buffer.data(mh + 721);
    const auto *mh_722 = buffer.data(mh + 722);
    const auto *mh_723 = buffer.data(mh + 723);
    const auto *mh_724 = buffer.data(mh + 724);
    const auto *mh_725 = buffer.data(mh + 725);
    const auto *mh_726 = buffer.data(mh + 726);
    const auto *mh_727 = buffer.data(mh + 727);
    const auto *mh_728 = buffer.data(mh + 728);
    const auto *mh_729 = buffer.data(mh + 729);
    const auto *mh_730 = buffer.data(mh + 730);
    const auto *mh_731 = buffer.data(mh + 731);
    const auto *mh_732 = buffer.data(mh + 732);
    const auto *mh_733 = buffer.data(mh + 733);
    const auto *mh_734 = buffer.data(mh + 734);
    const auto *mh_735 = buffer.data(mh + 735);
    const auto *mh_736 = buffer.data(mh + 736);
    const auto *mh_737 = buffer.data(mh + 737);
    const auto *mh_738 = buffer.data(mh + 738);
    const auto *mh_739 = buffer.data(mh + 739);
    const auto *mh_740 = buffer.data(mh + 740);
    const auto *mh_741 = buffer.data(mh + 741);
    const auto *mh_742 = buffer.data(mh + 742);
    const auto *mh_743 = buffer.data(mh + 743);
    const auto *mh_744 = buffer.data(mh + 744);
    const auto *mh_745 = buffer.data(mh + 745);
    const auto *mh_746 = buffer.data(mh + 746);
    const auto *mh_747 = buffer.data(mh + 747);
    const auto *mh_748 = buffer.data(mh + 748);
    const auto *mh_749 = buffer.data(mh + 749);

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, kh_600, kh_601, kh_602, kh_603, \
                         kh_604, mh_600, mh_601, mh_602, mh_603, \
                         mh_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -kh_600[k]
                   + f_0 * mh_600[k];

        t_601[k] = -kh_601[k]
                   + f_0 * mh_601[k];

        t_602[k] = -kh_602[k]
                   + f_0 * mh_602[k];

        t_603[k] = -kh_603[k]
                   + f_0 * mh_603[k];

        t_604[k] = -kh_604[k]
                   + f_0 * mh_604[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, kh_605, kh_606, kh_607, kh_608, \
                         kh_609, mh_605, mh_606, mh_607, mh_608, \
                         mh_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -kh_605[k]
                   + f_0 * mh_605[k];

        t_606[k] = -kh_606[k]
                   + f_0 * mh_606[k];

        t_607[k] = -kh_607[k]
                   + f_0 * mh_607[k];

        t_608[k] = -kh_608[k]
                   + f_0 * mh_608[k];

        t_609[k] = -kh_609[k]
                   + f_0 * mh_609[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, kh_610, kh_611, kh_612, kh_613, \
                         kh_614, mh_610, mh_611, mh_612, mh_613, \
                         mh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -kh_610[k]
                   + f_0 * mh_610[k];

        t_611[k] = -kh_611[k]
                   + f_0 * mh_611[k];

        t_612[k] = -kh_612[k]
                   + f_0 * mh_612[k];

        t_613[k] = -kh_613[k]
                   + f_0 * mh_613[k];

        t_614[k] = -kh_614[k]
                   + f_0 * mh_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, kh_615, kh_616, kh_617, kh_618, \
                         kh_619, mh_615, mh_616, mh_617, mh_618, \
                         mh_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -kh_615[k]
                   + f_0 * mh_615[k];

        t_616[k] = -kh_616[k]
                   + f_0 * mh_616[k];

        t_617[k] = -kh_617[k]
                   + f_0 * mh_617[k];

        t_618[k] = -kh_618[k]
                   + f_0 * mh_618[k];

        t_619[k] = -kh_619[k]
                   + f_0 * mh_619[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, kh_620, kh_621, kh_622, kh_623, \
                         kh_624, mh_620, mh_621, mh_622, mh_623, \
                         mh_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -kh_620[k]
                   + f_0 * mh_620[k];

        t_621[k] = -kh_621[k]
                   + f_0 * mh_621[k];

        t_622[k] = -kh_622[k]
                   + f_0 * mh_622[k];

        t_623[k] = -kh_623[k]
                   + f_0 * mh_623[k];

        t_624[k] = -kh_624[k]
                   + f_0 * mh_624[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, kh_625, kh_626, kh_627, kh_628, \
                         kh_629, mh_625, mh_626, mh_627, mh_628, \
                         mh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -kh_625[k]
                   + f_0 * mh_625[k];

        t_626[k] = -kh_626[k]
                   + f_0 * mh_626[k];

        t_627[k] = -kh_627[k]
                   + f_0 * mh_627[k];

        t_628[k] = -kh_628[k]
                   + f_0 * mh_628[k];

        t_629[k] = -kh_629[k]
                   + f_0 * mh_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, kh_630, kh_631, kh_632, kh_633, \
                         kh_634, mh_630, mh_631, mh_632, mh_633, \
                         mh_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -kh_630[k]
                   + f_0 * mh_630[k];

        t_631[k] = -kh_631[k]
                   + f_0 * mh_631[k];

        t_632[k] = -kh_632[k]
                   + f_0 * mh_632[k];

        t_633[k] = -kh_633[k]
                   + f_0 * mh_633[k];

        t_634[k] = -kh_634[k]
                   + f_0 * mh_634[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, kh_635, kh_636, kh_637, kh_638, \
                         kh_639, mh_635, mh_636, mh_637, mh_638, \
                         mh_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -kh_635[k]
                   + f_0 * mh_635[k];

        t_636[k] = -kh_636[k]
                   + f_0 * mh_636[k];

        t_637[k] = -kh_637[k]
                   + f_0 * mh_637[k];

        t_638[k] = -kh_638[k]
                   + f_0 * mh_638[k];

        t_639[k] = -kh_639[k]
                   + f_0 * mh_639[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, kh_640, kh_641, kh_642, kh_643, \
                         kh_644, mh_640, mh_641, mh_642, mh_643, \
                         mh_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -kh_640[k]
                   + f_0 * mh_640[k];

        t_641[k] = -kh_641[k]
                   + f_0 * mh_641[k];

        t_642[k] = -kh_642[k]
                   + f_0 * mh_642[k];

        t_643[k] = -kh_643[k]
                   + f_0 * mh_643[k];

        t_644[k] = -kh_644[k]
                   + f_0 * mh_644[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, kh_645, kh_646, kh_647, kh_648, \
                         kh_649, mh_645, mh_646, mh_647, mh_648, \
                         mh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -kh_645[k]
                   + f_0 * mh_645[k];

        t_646[k] = -kh_646[k]
                   + f_0 * mh_646[k];

        t_647[k] = -kh_647[k]
                   + f_0 * mh_647[k];

        t_648[k] = -kh_648[k]
                   + f_0 * mh_648[k];

        t_649[k] = -kh_649[k]
                   + f_0 * mh_649[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, kh_650, kh_651, kh_652, kh_653, \
                         kh_654, mh_650, mh_651, mh_652, mh_653, \
                         mh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -kh_650[k]
                   + f_0 * mh_650[k];

        t_651[k] = -kh_651[k]
                   + f_0 * mh_651[k];

        t_652[k] = -kh_652[k]
                   + f_0 * mh_652[k];

        t_653[k] = -kh_653[k]
                   + f_0 * mh_653[k];

        t_654[k] = -kh_654[k]
                   + f_0 * mh_654[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, kh_655, kh_656, kh_657, kh_658, \
                         kh_659, mh_655, mh_656, mh_657, mh_658, \
                         mh_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -kh_655[k]
                   + f_0 * mh_655[k];

        t_656[k] = -kh_656[k]
                   + f_0 * mh_656[k];

        t_657[k] = -kh_657[k]
                   + f_0 * mh_657[k];

        t_658[k] = -kh_658[k]
                   + f_0 * mh_658[k];

        t_659[k] = -kh_659[k]
                   + f_0 * mh_659[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, kh_660, kh_661, kh_662, kh_663, \
                         kh_664, mh_660, mh_661, mh_662, mh_663, \
                         mh_664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -kh_660[k]
                   + f_0 * mh_660[k];

        t_661[k] = -kh_661[k]
                   + f_0 * mh_661[k];

        t_662[k] = -kh_662[k]
                   + f_0 * mh_662[k];

        t_663[k] = -kh_663[k]
                   + f_0 * mh_663[k];

        t_664[k] = -kh_664[k]
                   + f_0 * mh_664[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, kh_665, kh_666, kh_667, kh_668, \
                         kh_669, mh_665, mh_666, mh_667, mh_668, \
                         mh_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -kh_665[k]
                   + f_0 * mh_665[k];

        t_666[k] = -kh_666[k]
                   + f_0 * mh_666[k];

        t_667[k] = -kh_667[k]
                   + f_0 * mh_667[k];

        t_668[k] = -kh_668[k]
                   + f_0 * mh_668[k];

        t_669[k] = -kh_669[k]
                   + f_0 * mh_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, kh_670, kh_671, kh_672, kh_673, \
                         kh_674, mh_670, mh_671, mh_672, mh_673, \
                         mh_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -kh_670[k]
                   + f_0 * mh_670[k];

        t_671[k] = -kh_671[k]
                   + f_0 * mh_671[k];

        t_672[k] = -kh_672[k]
                   + f_0 * mh_672[k];

        t_673[k] = -kh_673[k]
                   + f_0 * mh_673[k];

        t_674[k] = -kh_674[k]
                   + f_0 * mh_674[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, kh_675, kh_676, kh_677, kh_678, \
                         kh_679, mh_675, mh_676, mh_677, mh_678, \
                         mh_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -kh_675[k]
                   + f_0 * mh_675[k];

        t_676[k] = -kh_676[k]
                   + f_0 * mh_676[k];

        t_677[k] = -kh_677[k]
                   + f_0 * mh_677[k];

        t_678[k] = -kh_678[k]
                   + f_0 * mh_678[k];

        t_679[k] = -kh_679[k]
                   + f_0 * mh_679[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, kh_680, kh_681, kh_682, kh_683, \
                         kh_684, mh_680, mh_681, mh_682, mh_683, \
                         mh_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -kh_680[k]
                   + f_0 * mh_680[k];

        t_681[k] = -kh_681[k]
                   + f_0 * mh_681[k];

        t_682[k] = -kh_682[k]
                   + f_0 * mh_682[k];

        t_683[k] = -kh_683[k]
                   + f_0 * mh_683[k];

        t_684[k] = -kh_684[k]
                   + f_0 * mh_684[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, kh_685, kh_686, kh_687, kh_688, \
                         kh_689, mh_685, mh_686, mh_687, mh_688, \
                         mh_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -kh_685[k]
                   + f_0 * mh_685[k];

        t_686[k] = -kh_686[k]
                   + f_0 * mh_686[k];

        t_687[k] = -kh_687[k]
                   + f_0 * mh_687[k];

        t_688[k] = -kh_688[k]
                   + f_0 * mh_688[k];

        t_689[k] = -kh_689[k]
                   + f_0 * mh_689[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, kh_690, kh_691, kh_692, kh_693, \
                         kh_694, mh_690, mh_691, mh_692, mh_693, \
                         mh_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -kh_690[k]
                   + f_0 * mh_690[k];

        t_691[k] = -kh_691[k]
                   + f_0 * mh_691[k];

        t_692[k] = -kh_692[k]
                   + f_0 * mh_692[k];

        t_693[k] = -kh_693[k]
                   + f_0 * mh_693[k];

        t_694[k] = -kh_694[k]
                   + f_0 * mh_694[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, kh_695, kh_696, kh_697, kh_698, \
                         kh_699, mh_695, mh_696, mh_697, mh_698, \
                         mh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -kh_695[k]
                   + f_0 * mh_695[k];

        t_696[k] = -kh_696[k]
                   + f_0 * mh_696[k];

        t_697[k] = -kh_697[k]
                   + f_0 * mh_697[k];

        t_698[k] = -kh_698[k]
                   + f_0 * mh_698[k];

        t_699[k] = -kh_699[k]
                   + f_0 * mh_699[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, kh_700, kh_701, kh_702, kh_703, \
                         kh_704, mh_700, mh_701, mh_702, mh_703, \
                         mh_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -kh_700[k]
                   + f_0 * mh_700[k];

        t_701[k] = -kh_701[k]
                   + f_0 * mh_701[k];

        t_702[k] = -kh_702[k]
                   + f_0 * mh_702[k];

        t_703[k] = -kh_703[k]
                   + f_0 * mh_703[k];

        t_704[k] = -kh_704[k]
                   + f_0 * mh_704[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, kh_705, kh_706, kh_707, kh_708, \
                         kh_709, mh_705, mh_706, mh_707, mh_708, \
                         mh_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -kh_705[k]
                   + f_0 * mh_705[k];

        t_706[k] = -kh_706[k]
                   + f_0 * mh_706[k];

        t_707[k] = -kh_707[k]
                   + f_0 * mh_707[k];

        t_708[k] = -kh_708[k]
                   + f_0 * mh_708[k];

        t_709[k] = -kh_709[k]
                   + f_0 * mh_709[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, kh_710, kh_711, kh_712, kh_713, \
                         kh_714, mh_710, mh_711, mh_712, mh_713, \
                         mh_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -kh_710[k]
                   + f_0 * mh_710[k];

        t_711[k] = -kh_711[k]
                   + f_0 * mh_711[k];

        t_712[k] = -kh_712[k]
                   + f_0 * mh_712[k];

        t_713[k] = -kh_713[k]
                   + f_0 * mh_713[k];

        t_714[k] = -kh_714[k]
                   + f_0 * mh_714[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, kh_715, kh_716, kh_717, kh_718, \
                         kh_719, mh_715, mh_716, mh_717, mh_718, \
                         mh_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -kh_715[k]
                   + f_0 * mh_715[k];

        t_716[k] = -kh_716[k]
                   + f_0 * mh_716[k];

        t_717[k] = -kh_717[k]
                   + f_0 * mh_717[k];

        t_718[k] = -kh_718[k]
                   + f_0 * mh_718[k];

        t_719[k] = -kh_719[k]
                   + f_0 * mh_719[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, kh_720, kh_721, kh_722, kh_723, \
                         kh_724, mh_720, mh_721, mh_722, mh_723, \
                         mh_724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -kh_720[k]
                   + f_0 * mh_720[k];

        t_721[k] = -kh_721[k]
                   + f_0 * mh_721[k];

        t_722[k] = -kh_722[k]
                   + f_0 * mh_722[k];

        t_723[k] = -kh_723[k]
                   + f_0 * mh_723[k];

        t_724[k] = -kh_724[k]
                   + f_0 * mh_724[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, kh_725, kh_726, kh_727, kh_728, \
                         kh_729, mh_725, mh_726, mh_727, mh_728, \
                         mh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -kh_725[k]
                   + f_0 * mh_725[k];

        t_726[k] = -kh_726[k]
                   + f_0 * mh_726[k];

        t_727[k] = -kh_727[k]
                   + f_0 * mh_727[k];

        t_728[k] = -kh_728[k]
                   + f_0 * mh_728[k];

        t_729[k] = -kh_729[k]
                   + f_0 * mh_729[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, kh_730, kh_731, kh_732, kh_733, \
                         kh_734, mh_730, mh_731, mh_732, mh_733, \
                         mh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -kh_730[k]
                   + f_0 * mh_730[k];

        t_731[k] = -kh_731[k]
                   + f_0 * mh_731[k];

        t_732[k] = -kh_732[k]
                   + f_0 * mh_732[k];

        t_733[k] = -kh_733[k]
                   + f_0 * mh_733[k];

        t_734[k] = -kh_734[k]
                   + f_0 * mh_734[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, kh_735, kh_736, kh_737, kh_738, \
                         kh_739, mh_735, mh_736, mh_737, mh_738, \
                         mh_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -kh_735[k]
                   + f_0 * mh_735[k];

        t_736[k] = -kh_736[k]
                   + f_0 * mh_736[k];

        t_737[k] = -kh_737[k]
                   + f_0 * mh_737[k];

        t_738[k] = -kh_738[k]
                   + f_0 * mh_738[k];

        t_739[k] = -kh_739[k]
                   + f_0 * mh_739[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, kh_740, kh_741, kh_742, kh_743, \
                         kh_744, mh_740, mh_741, mh_742, mh_743, \
                         mh_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -kh_740[k]
                   + f_0 * mh_740[k];

        t_741[k] = -kh_741[k]
                   + f_0 * mh_741[k];

        t_742[k] = -kh_742[k]
                   + f_0 * mh_742[k];

        t_743[k] = -kh_743[k]
                   + f_0 * mh_743[k];

        t_744[k] = -kh_744[k]
                   + f_0 * mh_744[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, kh_745, kh_746, kh_747, kh_748, \
                         kh_749, mh_745, mh_746, mh_747, mh_748, \
                         mh_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -kh_745[k]
                   + f_0 * mh_745[k];

        t_746[k] = -kh_746[k]
                   + f_0 * mh_746[k];

        t_747[k] = -kh_747[k]
                   + f_0 * mh_747[k];

        t_748[k] = -kh_748[k]
                   + f_0 * mh_748[k];

        t_749[k] = -kh_749[k]
                   + f_0 * mh_749[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_750 = buffer.data(target + 750);
    auto *t_751 = buffer.data(target + 751);
    auto *t_752 = buffer.data(target + 752);
    auto *t_753 = buffer.data(target + 753);
    auto *t_754 = buffer.data(target + 754);
    auto *t_755 = buffer.data(target + 755);
    auto *t_756 = buffer.data(target + 756);
    auto *t_757 = buffer.data(target + 757);
    auto *t_758 = buffer.data(target + 758);
    auto *t_759 = buffer.data(target + 759);
    auto *t_760 = buffer.data(target + 760);
    auto *t_761 = buffer.data(target + 761);
    auto *t_762 = buffer.data(target + 762);
    auto *t_763 = buffer.data(target + 763);
    auto *t_764 = buffer.data(target + 764);
    auto *t_765 = buffer.data(target + 765);
    auto *t_766 = buffer.data(target + 766);
    auto *t_767 = buffer.data(target + 767);
    auto *t_768 = buffer.data(target + 768);
    auto *t_769 = buffer.data(target + 769);
    auto *t_770 = buffer.data(target + 770);
    auto *t_771 = buffer.data(target + 771);
    auto *t_772 = buffer.data(target + 772);
    auto *t_773 = buffer.data(target + 773);
    auto *t_774 = buffer.data(target + 774);
    auto *t_775 = buffer.data(target + 775);
    auto *t_776 = buffer.data(target + 776);
    auto *t_777 = buffer.data(target + 777);
    auto *t_778 = buffer.data(target + 778);
    auto *t_779 = buffer.data(target + 779);
    auto *t_780 = buffer.data(target + 780);
    auto *t_781 = buffer.data(target + 781);
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);
    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);
    auto *t_801 = buffer.data(target + 801);
    auto *t_802 = buffer.data(target + 802);
    auto *t_803 = buffer.data(target + 803);
    auto *t_804 = buffer.data(target + 804);
    auto *t_805 = buffer.data(target + 805);
    auto *t_806 = buffer.data(target + 806);
    auto *t_807 = buffer.data(target + 807);
    auto *t_808 = buffer.data(target + 808);
    auto *t_809 = buffer.data(target + 809);
    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);
    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);
    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);
    auto *t_840 = buffer.data(target + 840);
    auto *t_841 = buffer.data(target + 841);
    auto *t_842 = buffer.data(target + 842);
    auto *t_843 = buffer.data(target + 843);
    auto *t_844 = buffer.data(target + 844);
    auto *t_845 = buffer.data(target + 845);
    auto *t_846 = buffer.data(target + 846);
    auto *t_847 = buffer.data(target + 847);
    auto *t_848 = buffer.data(target + 848);
    auto *t_849 = buffer.data(target + 849);
    auto *t_850 = buffer.data(target + 850);
    auto *t_851 = buffer.data(target + 851);
    auto *t_852 = buffer.data(target + 852);
    auto *t_853 = buffer.data(target + 853);
    auto *t_854 = buffer.data(target + 854);
    auto *t_855 = buffer.data(target + 855);
    auto *t_856 = buffer.data(target + 856);
    auto *t_857 = buffer.data(target + 857);
    auto *t_858 = buffer.data(target + 858);
    auto *t_859 = buffer.data(target + 859);
    auto *t_860 = buffer.data(target + 860);
    auto *t_861 = buffer.data(target + 861);
    auto *t_862 = buffer.data(target + 862);
    auto *t_863 = buffer.data(target + 863);
    auto *t_864 = buffer.data(target + 864);
    auto *t_865 = buffer.data(target + 865);
    auto *t_866 = buffer.data(target + 866);
    auto *t_867 = buffer.data(target + 867);
    auto *t_868 = buffer.data(target + 868);
    auto *t_869 = buffer.data(target + 869);
    auto *t_870 = buffer.data(target + 870);
    auto *t_871 = buffer.data(target + 871);
    auto *t_872 = buffer.data(target + 872);
    auto *t_873 = buffer.data(target + 873);
    auto *t_874 = buffer.data(target + 874);
    auto *t_875 = buffer.data(target + 875);
    auto *t_876 = buffer.data(target + 876);
    auto *t_877 = buffer.data(target + 877);
    auto *t_878 = buffer.data(target + 878);
    auto *t_879 = buffer.data(target + 879);
    auto *t_880 = buffer.data(target + 880);
    auto *t_881 = buffer.data(target + 881);
    auto *t_882 = buffer.data(target + 882);
    auto *t_883 = buffer.data(target + 883);
    auto *t_884 = buffer.data(target + 884);
    auto *t_885 = buffer.data(target + 885);
    auto *t_886 = buffer.data(target + 886);
    auto *t_887 = buffer.data(target + 887);
    auto *t_888 = buffer.data(target + 888);
    auto *t_889 = buffer.data(target + 889);
    auto *t_890 = buffer.data(target + 890);
    auto *t_891 = buffer.data(target + 891);
    auto *t_892 = buffer.data(target + 892);
    auto *t_893 = buffer.data(target + 893);
    auto *t_894 = buffer.data(target + 894);
    auto *t_895 = buffer.data(target + 895);
    auto *t_896 = buffer.data(target + 896);
    auto *t_897 = buffer.data(target + 897);
    auto *t_898 = buffer.data(target + 898);
    auto *t_899 = buffer.data(target + 899);
    auto *t_900 = buffer.data(target + 900);
    auto *t_901 = buffer.data(target + 901);
    auto *t_902 = buffer.data(target + 902);
    auto *t_903 = buffer.data(target + 903);
    auto *t_904 = buffer.data(target + 904);
    auto *t_905 = buffer.data(target + 905);
    auto *t_906 = buffer.data(target + 906);
    auto *t_907 = buffer.data(target + 907);
    auto *t_908 = buffer.data(target + 908);
    auto *t_909 = buffer.data(target + 909);
    auto *t_910 = buffer.data(target + 910);
    auto *t_911 = buffer.data(target + 911);
    auto *t_912 = buffer.data(target + 912);
    auto *t_913 = buffer.data(target + 913);
    auto *t_914 = buffer.data(target + 914);
    auto *t_915 = buffer.data(target + 915);
    auto *t_916 = buffer.data(target + 916);
    auto *t_917 = buffer.data(target + 917);
    auto *t_918 = buffer.data(target + 918);
    auto *t_919 = buffer.data(target + 919);
    auto *t_920 = buffer.data(target + 920);
    auto *t_921 = buffer.data(target + 921);
    auto *t_922 = buffer.data(target + 922);
    auto *t_923 = buffer.data(target + 923);
    auto *t_924 = buffer.data(target + 924);
    auto *t_925 = buffer.data(target + 925);
    auto *t_926 = buffer.data(target + 926);
    auto *t_927 = buffer.data(target + 927);
    auto *t_928 = buffer.data(target + 928);
    auto *t_929 = buffer.data(target + 929);
    auto *t_930 = buffer.data(target + 930);
    auto *t_931 = buffer.data(target + 931);
    auto *t_932 = buffer.data(target + 932);
    auto *t_933 = buffer.data(target + 933);
    auto *t_934 = buffer.data(target + 934);
    auto *t_935 = buffer.data(target + 935);
    auto *t_936 = buffer.data(target + 936);
    auto *t_937 = buffer.data(target + 937);
    auto *t_938 = buffer.data(target + 938);
    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);

    const auto *kh_750 = buffer.data(kh + 750);
    const auto *kh_751 = buffer.data(kh + 751);
    const auto *kh_752 = buffer.data(kh + 752);
    const auto *kh_753 = buffer.data(kh + 753);
    const auto *kh_754 = buffer.data(kh + 754);
    const auto *kh_755 = buffer.data(kh + 755);

    const auto *mh_750 = buffer.data(mh + 750);
    const auto *mh_751 = buffer.data(mh + 751);
    const auto *mh_752 = buffer.data(mh + 752);
    const auto *mh_753 = buffer.data(mh + 753);
    const auto *mh_754 = buffer.data(mh + 754);
    const auto *mh_755 = buffer.data(mh + 755);
    const auto *mh_756 = buffer.data(mh + 756);
    const auto *mh_757 = buffer.data(mh + 757);
    const auto *mh_758 = buffer.data(mh + 758);
    const auto *mh_759 = buffer.data(mh + 759);
    const auto *mh_760 = buffer.data(mh + 760);
    const auto *mh_761 = buffer.data(mh + 761);
    const auto *mh_762 = buffer.data(mh + 762);
    const auto *mh_763 = buffer.data(mh + 763);
    const auto *mh_764 = buffer.data(mh + 764);
    const auto *mh_765 = buffer.data(mh + 765);
    const auto *mh_766 = buffer.data(mh + 766);
    const auto *mh_767 = buffer.data(mh + 767);
    const auto *mh_768 = buffer.data(mh + 768);
    const auto *mh_769 = buffer.data(mh + 769);
    const auto *mh_770 = buffer.data(mh + 770);
    const auto *mh_771 = buffer.data(mh + 771);
    const auto *mh_772 = buffer.data(mh + 772);
    const auto *mh_773 = buffer.data(mh + 773);
    const auto *mh_774 = buffer.data(mh + 774);
    const auto *mh_775 = buffer.data(mh + 775);
    const auto *mh_776 = buffer.data(mh + 776);
    const auto *mh_777 = buffer.data(mh + 777);
    const auto *mh_778 = buffer.data(mh + 778);
    const auto *mh_779 = buffer.data(mh + 779);
    const auto *mh_780 = buffer.data(mh + 780);
    const auto *mh_781 = buffer.data(mh + 781);
    const auto *mh_782 = buffer.data(mh + 782);
    const auto *mh_783 = buffer.data(mh + 783);
    const auto *mh_784 = buffer.data(mh + 784);
    const auto *mh_785 = buffer.data(mh + 785);
    const auto *mh_786 = buffer.data(mh + 786);
    const auto *mh_787 = buffer.data(mh + 787);
    const auto *mh_788 = buffer.data(mh + 788);
    const auto *mh_789 = buffer.data(mh + 789);
    const auto *mh_790 = buffer.data(mh + 790);
    const auto *mh_791 = buffer.data(mh + 791);
    const auto *mh_792 = buffer.data(mh + 792);
    const auto *mh_793 = buffer.data(mh + 793);
    const auto *mh_794 = buffer.data(mh + 794);
    const auto *mh_795 = buffer.data(mh + 795);
    const auto *mh_796 = buffer.data(mh + 796);
    const auto *mh_797 = buffer.data(mh + 797);
    const auto *mh_798 = buffer.data(mh + 798);
    const auto *mh_799 = buffer.data(mh + 799);
    const auto *mh_800 = buffer.data(mh + 800);
    const auto *mh_801 = buffer.data(mh + 801);
    const auto *mh_802 = buffer.data(mh + 802);
    const auto *mh_803 = buffer.data(mh + 803);
    const auto *mh_804 = buffer.data(mh + 804);
    const auto *mh_805 = buffer.data(mh + 805);
    const auto *mh_806 = buffer.data(mh + 806);
    const auto *mh_807 = buffer.data(mh + 807);
    const auto *mh_808 = buffer.data(mh + 808);
    const auto *mh_809 = buffer.data(mh + 809);
    const auto *mh_810 = buffer.data(mh + 810);
    const auto *mh_811 = buffer.data(mh + 811);
    const auto *mh_812 = buffer.data(mh + 812);
    const auto *mh_813 = buffer.data(mh + 813);
    const auto *mh_814 = buffer.data(mh + 814);
    const auto *mh_815 = buffer.data(mh + 815);
    const auto *mh_816 = buffer.data(mh + 816);
    const auto *mh_817 = buffer.data(mh + 817);
    const auto *mh_818 = buffer.data(mh + 818);
    const auto *mh_819 = buffer.data(mh + 819);
    const auto *mh_820 = buffer.data(mh + 820);
    const auto *mh_821 = buffer.data(mh + 821);
    const auto *mh_822 = buffer.data(mh + 822);
    const auto *mh_823 = buffer.data(mh + 823);
    const auto *mh_824 = buffer.data(mh + 824);
    const auto *mh_825 = buffer.data(mh + 825);
    const auto *mh_826 = buffer.data(mh + 826);
    const auto *mh_827 = buffer.data(mh + 827);
    const auto *mh_828 = buffer.data(mh + 828);
    const auto *mh_829 = buffer.data(mh + 829);
    const auto *mh_830 = buffer.data(mh + 830);
    const auto *mh_831 = buffer.data(mh + 831);
    const auto *mh_832 = buffer.data(mh + 832);
    const auto *mh_833 = buffer.data(mh + 833);
    const auto *mh_834 = buffer.data(mh + 834);
    const auto *mh_835 = buffer.data(mh + 835);
    const auto *mh_836 = buffer.data(mh + 836);
    const auto *mh_837 = buffer.data(mh + 837);
    const auto *mh_838 = buffer.data(mh + 838);
    const auto *mh_839 = buffer.data(mh + 839);
    const auto *mh_840 = buffer.data(mh + 840);
    const auto *mh_841 = buffer.data(mh + 841);
    const auto *mh_842 = buffer.data(mh + 842);
    const auto *mh_843 = buffer.data(mh + 843);
    const auto *mh_844 = buffer.data(mh + 844);
    const auto *mh_845 = buffer.data(mh + 845);
    const auto *mh_846 = buffer.data(mh + 846);
    const auto *mh_847 = buffer.data(mh + 847);
    const auto *mh_848 = buffer.data(mh + 848);
    const auto *mh_849 = buffer.data(mh + 849);
    const auto *mh_850 = buffer.data(mh + 850);
    const auto *mh_851 = buffer.data(mh + 851);
    const auto *mh_852 = buffer.data(mh + 852);
    const auto *mh_853 = buffer.data(mh + 853);
    const auto *mh_854 = buffer.data(mh + 854);
    const auto *mh_855 = buffer.data(mh + 855);
    const auto *mh_856 = buffer.data(mh + 856);
    const auto *mh_857 = buffer.data(mh + 857);
    const auto *mh_858 = buffer.data(mh + 858);
    const auto *mh_859 = buffer.data(mh + 859);
    const auto *mh_860 = buffer.data(mh + 860);
    const auto *mh_861 = buffer.data(mh + 861);
    const auto *mh_862 = buffer.data(mh + 862);
    const auto *mh_863 = buffer.data(mh + 863);
    const auto *mh_864 = buffer.data(mh + 864);
    const auto *mh_865 = buffer.data(mh + 865);
    const auto *mh_866 = buffer.data(mh + 866);
    const auto *mh_867 = buffer.data(mh + 867);
    const auto *mh_868 = buffer.data(mh + 868);
    const auto *mh_869 = buffer.data(mh + 869);
    const auto *mh_870 = buffer.data(mh + 870);
    const auto *mh_871 = buffer.data(mh + 871);
    const auto *mh_872 = buffer.data(mh + 872);
    const auto *mh_873 = buffer.data(mh + 873);
    const auto *mh_874 = buffer.data(mh + 874);
    const auto *mh_875 = buffer.data(mh + 875);
    const auto *mh_876 = buffer.data(mh + 876);
    const auto *mh_877 = buffer.data(mh + 877);
    const auto *mh_878 = buffer.data(mh + 878);
    const auto *mh_879 = buffer.data(mh + 879);
    const auto *mh_880 = buffer.data(mh + 880);
    const auto *mh_881 = buffer.data(mh + 881);
    const auto *mh_882 = buffer.data(mh + 882);
    const auto *mh_883 = buffer.data(mh + 883);
    const auto *mh_884 = buffer.data(mh + 884);
    const auto *mh_885 = buffer.data(mh + 885);
    const auto *mh_886 = buffer.data(mh + 886);
    const auto *mh_887 = buffer.data(mh + 887);
    const auto *mh_888 = buffer.data(mh + 888);
    const auto *mh_889 = buffer.data(mh + 889);
    const auto *mh_890 = buffer.data(mh + 890);
    const auto *mh_891 = buffer.data(mh + 891);
    const auto *mh_892 = buffer.data(mh + 892);
    const auto *mh_893 = buffer.data(mh + 893);
    const auto *mh_894 = buffer.data(mh + 894);
    const auto *mh_895 = buffer.data(mh + 895);
    const auto *mh_896 = buffer.data(mh + 896);
    const auto *mh_897 = buffer.data(mh + 897);
    const auto *mh_898 = buffer.data(mh + 898);
    const auto *mh_899 = buffer.data(mh + 899);
    const auto *mh_900 = buffer.data(mh + 900);
    const auto *mh_901 = buffer.data(mh + 901);
    const auto *mh_902 = buffer.data(mh + 902);
    const auto *mh_903 = buffer.data(mh + 903);
    const auto *mh_904 = buffer.data(mh + 904);
    const auto *mh_905 = buffer.data(mh + 905);
    const auto *mh_906 = buffer.data(mh + 906);
    const auto *mh_907 = buffer.data(mh + 907);
    const auto *mh_908 = buffer.data(mh + 908);
    const auto *mh_909 = buffer.data(mh + 909);
    const auto *mh_910 = buffer.data(mh + 910);
    const auto *mh_911 = buffer.data(mh + 911);
    const auto *mh_912 = buffer.data(mh + 912);
    const auto *mh_913 = buffer.data(mh + 913);
    const auto *mh_914 = buffer.data(mh + 914);
    const auto *mh_915 = buffer.data(mh + 915);
    const auto *mh_916 = buffer.data(mh + 916);
    const auto *mh_917 = buffer.data(mh + 917);
    const auto *mh_918 = buffer.data(mh + 918);
    const auto *mh_919 = buffer.data(mh + 919);
    const auto *mh_920 = buffer.data(mh + 920);
    const auto *mh_921 = buffer.data(mh + 921);
    const auto *mh_922 = buffer.data(mh + 922);
    const auto *mh_923 = buffer.data(mh + 923);
    const auto *mh_924 = buffer.data(mh + 924);
    const auto *mh_925 = buffer.data(mh + 925);
    const auto *mh_926 = buffer.data(mh + 926);
    const auto *mh_927 = buffer.data(mh + 927);
    const auto *mh_928 = buffer.data(mh + 928);
    const auto *mh_929 = buffer.data(mh + 929);
    const auto *mh_930 = buffer.data(mh + 930);
    const auto *mh_931 = buffer.data(mh + 931);
    const auto *mh_932 = buffer.data(mh + 932);
    const auto *mh_933 = buffer.data(mh + 933);
    const auto *mh_934 = buffer.data(mh + 934);
    const auto *mh_935 = buffer.data(mh + 935);
    const auto *mh_936 = buffer.data(mh + 936);
    const auto *mh_937 = buffer.data(mh + 937);
    const auto *mh_938 = buffer.data(mh + 938);
    const auto *mh_939 = buffer.data(mh + 939);
    const auto *mh_940 = buffer.data(mh + 940);
    const auto *mh_941 = buffer.data(mh + 941);
    const auto *mh_942 = buffer.data(mh + 942);
    const auto *mh_943 = buffer.data(mh + 943);
    const auto *mh_944 = buffer.data(mh + 944);

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, kh_750, kh_751, kh_752, kh_753, \
                         kh_754, mh_750, mh_751, mh_752, mh_753, \
                         mh_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -kh_750[k]
                   + f_0 * mh_750[k];

        t_751[k] = -kh_751[k]
                   + f_0 * mh_751[k];

        t_752[k] = -kh_752[k]
                   + f_0 * mh_752[k];

        t_753[k] = -kh_753[k]
                   + f_0 * mh_753[k];

        t_754[k] = -kh_754[k]
                   + f_0 * mh_754[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, t_760, t_761, kh_755, mh_755, \
                         mh_756, mh_757, mh_758, mh_759, mh_760, \
                         mh_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -kh_755[k]
                   + f_0 * mh_755[k];

        t_756[k] = f_0 * mh_756[k];

        t_757[k] = f_0 * mh_757[k];

        t_758[k] = f_0 * mh_758[k];

        t_759[k] = f_0 * mh_759[k];

        t_760[k] = f_0 * mh_760[k];

        t_761[k] = f_0 * mh_761[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, t_766, t_767, t_768, t_769, mh_762, \
                         mh_763, mh_764, mh_765, mh_766, mh_767, mh_768, \
                         mh_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_0 * mh_762[k];

        t_763[k] = f_0 * mh_763[k];

        t_764[k] = f_0 * mh_764[k];

        t_765[k] = f_0 * mh_765[k];

        t_766[k] = f_0 * mh_766[k];

        t_767[k] = f_0 * mh_767[k];

        t_768[k] = f_0 * mh_768[k];

        t_769[k] = f_0 * mh_769[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, t_775, t_776, t_777, mh_770, \
                         mh_771, mh_772, mh_773, mh_774, mh_775, mh_776, \
                         mh_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_0 * mh_770[k];

        t_771[k] = f_0 * mh_771[k];

        t_772[k] = f_0 * mh_772[k];

        t_773[k] = f_0 * mh_773[k];

        t_774[k] = f_0 * mh_774[k];

        t_775[k] = f_0 * mh_775[k];

        t_776[k] = f_0 * mh_776[k];

        t_777[k] = f_0 * mh_777[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, t_783, t_784, t_785, mh_778, \
                         mh_779, mh_780, mh_781, mh_782, mh_783, mh_784, \
                         mh_785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_0 * mh_778[k];

        t_779[k] = f_0 * mh_779[k];

        t_780[k] = f_0 * mh_780[k];

        t_781[k] = f_0 * mh_781[k];

        t_782[k] = f_0 * mh_782[k];

        t_783[k] = f_0 * mh_783[k];

        t_784[k] = f_0 * mh_784[k];

        t_785[k] = f_0 * mh_785[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, t_790, t_791, t_792, t_793, mh_786, \
                         mh_787, mh_788, mh_789, mh_790, mh_791, mh_792, \
                         mh_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_0 * mh_786[k];

        t_787[k] = f_0 * mh_787[k];

        t_788[k] = f_0 * mh_788[k];

        t_789[k] = f_0 * mh_789[k];

        t_790[k] = f_0 * mh_790[k];

        t_791[k] = f_0 * mh_791[k];

        t_792[k] = f_0 * mh_792[k];

        t_793[k] = f_0 * mh_793[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, t_798, t_799, t_800, t_801, mh_794, \
                         mh_795, mh_796, mh_797, mh_798, mh_799, mh_800, \
                         mh_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_0 * mh_794[k];

        t_795[k] = f_0 * mh_795[k];

        t_796[k] = f_0 * mh_796[k];

        t_797[k] = f_0 * mh_797[k];

        t_798[k] = f_0 * mh_798[k];

        t_799[k] = f_0 * mh_799[k];

        t_800[k] = f_0 * mh_800[k];

        t_801[k] = f_0 * mh_801[k];
    }

#pragma omp simd aligned(t_802, t_803, t_804, t_805, t_806, t_807, t_808, t_809, mh_802, \
                         mh_803, mh_804, mh_805, mh_806, mh_807, mh_808, \
                         mh_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_802[k] = f_0 * mh_802[k];

        t_803[k] = f_0 * mh_803[k];

        t_804[k] = f_0 * mh_804[k];

        t_805[k] = f_0 * mh_805[k];

        t_806[k] = f_0 * mh_806[k];

        t_807[k] = f_0 * mh_807[k];

        t_808[k] = f_0 * mh_808[k];

        t_809[k] = f_0 * mh_809[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, t_815, t_816, t_817, mh_810, \
                         mh_811, mh_812, mh_813, mh_814, mh_815, mh_816, \
                         mh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_0 * mh_810[k];

        t_811[k] = f_0 * mh_811[k];

        t_812[k] = f_0 * mh_812[k];

        t_813[k] = f_0 * mh_813[k];

        t_814[k] = f_0 * mh_814[k];

        t_815[k] = f_0 * mh_815[k];

        t_816[k] = f_0 * mh_816[k];

        t_817[k] = f_0 * mh_817[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, t_821, t_822, t_823, t_824, t_825, mh_818, \
                         mh_819, mh_820, mh_821, mh_822, mh_823, mh_824, \
                         mh_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_0 * mh_818[k];

        t_819[k] = f_0 * mh_819[k];

        t_820[k] = f_0 * mh_820[k];

        t_821[k] = f_0 * mh_821[k];

        t_822[k] = f_0 * mh_822[k];

        t_823[k] = f_0 * mh_823[k];

        t_824[k] = f_0 * mh_824[k];

        t_825[k] = f_0 * mh_825[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, t_830, t_831, t_832, t_833, mh_826, \
                         mh_827, mh_828, mh_829, mh_830, mh_831, mh_832, \
                         mh_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_0 * mh_826[k];

        t_827[k] = f_0 * mh_827[k];

        t_828[k] = f_0 * mh_828[k];

        t_829[k] = f_0 * mh_829[k];

        t_830[k] = f_0 * mh_830[k];

        t_831[k] = f_0 * mh_831[k];

        t_832[k] = f_0 * mh_832[k];

        t_833[k] = f_0 * mh_833[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, t_837, t_838, t_839, t_840, t_841, mh_834, \
                         mh_835, mh_836, mh_837, mh_838, mh_839, mh_840, \
                         mh_841 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_0 * mh_834[k];

        t_835[k] = f_0 * mh_835[k];

        t_836[k] = f_0 * mh_836[k];

        t_837[k] = f_0 * mh_837[k];

        t_838[k] = f_0 * mh_838[k];

        t_839[k] = f_0 * mh_839[k];

        t_840[k] = f_0 * mh_840[k];

        t_841[k] = f_0 * mh_841[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, t_845, t_846, t_847, t_848, t_849, mh_842, \
                         mh_843, mh_844, mh_845, mh_846, mh_847, mh_848, \
                         mh_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_0 * mh_842[k];

        t_843[k] = f_0 * mh_843[k];

        t_844[k] = f_0 * mh_844[k];

        t_845[k] = f_0 * mh_845[k];

        t_846[k] = f_0 * mh_846[k];

        t_847[k] = f_0 * mh_847[k];

        t_848[k] = f_0 * mh_848[k];

        t_849[k] = f_0 * mh_849[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, t_855, t_856, t_857, mh_850, \
                         mh_851, mh_852, mh_853, mh_854, mh_855, mh_856, \
                         mh_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = f_0 * mh_850[k];

        t_851[k] = f_0 * mh_851[k];

        t_852[k] = f_0 * mh_852[k];

        t_853[k] = f_0 * mh_853[k];

        t_854[k] = f_0 * mh_854[k];

        t_855[k] = f_0 * mh_855[k];

        t_856[k] = f_0 * mh_856[k];

        t_857[k] = f_0 * mh_857[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, t_862, t_863, t_864, t_865, mh_858, \
                         mh_859, mh_860, mh_861, mh_862, mh_863, mh_864, \
                         mh_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_0 * mh_858[k];

        t_859[k] = f_0 * mh_859[k];

        t_860[k] = f_0 * mh_860[k];

        t_861[k] = f_0 * mh_861[k];

        t_862[k] = f_0 * mh_862[k];

        t_863[k] = f_0 * mh_863[k];

        t_864[k] = f_0 * mh_864[k];

        t_865[k] = f_0 * mh_865[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, t_870, t_871, t_872, t_873, mh_866, \
                         mh_867, mh_868, mh_869, mh_870, mh_871, mh_872, \
                         mh_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_0 * mh_866[k];

        t_867[k] = f_0 * mh_867[k];

        t_868[k] = f_0 * mh_868[k];

        t_869[k] = f_0 * mh_869[k];

        t_870[k] = f_0 * mh_870[k];

        t_871[k] = f_0 * mh_871[k];

        t_872[k] = f_0 * mh_872[k];

        t_873[k] = f_0 * mh_873[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, t_878, t_879, t_880, t_881, mh_874, \
                         mh_875, mh_876, mh_877, mh_878, mh_879, mh_880, \
                         mh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_0 * mh_874[k];

        t_875[k] = f_0 * mh_875[k];

        t_876[k] = f_0 * mh_876[k];

        t_877[k] = f_0 * mh_877[k];

        t_878[k] = f_0 * mh_878[k];

        t_879[k] = f_0 * mh_879[k];

        t_880[k] = f_0 * mh_880[k];

        t_881[k] = f_0 * mh_881[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, t_886, t_887, t_888, t_889, mh_882, \
                         mh_883, mh_884, mh_885, mh_886, mh_887, mh_888, \
                         mh_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_0 * mh_882[k];

        t_883[k] = f_0 * mh_883[k];

        t_884[k] = f_0 * mh_884[k];

        t_885[k] = f_0 * mh_885[k];

        t_886[k] = f_0 * mh_886[k];

        t_887[k] = f_0 * mh_887[k];

        t_888[k] = f_0 * mh_888[k];

        t_889[k] = f_0 * mh_889[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, t_895, t_896, t_897, mh_890, \
                         mh_891, mh_892, mh_893, mh_894, mh_895, mh_896, \
                         mh_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_0 * mh_890[k];

        t_891[k] = f_0 * mh_891[k];

        t_892[k] = f_0 * mh_892[k];

        t_893[k] = f_0 * mh_893[k];

        t_894[k] = f_0 * mh_894[k];

        t_895[k] = f_0 * mh_895[k];

        t_896[k] = f_0 * mh_896[k];

        t_897[k] = f_0 * mh_897[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, t_902, t_903, t_904, t_905, mh_898, \
                         mh_899, mh_900, mh_901, mh_902, mh_903, mh_904, \
                         mh_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_0 * mh_898[k];

        t_899[k] = f_0 * mh_899[k];

        t_900[k] = f_0 * mh_900[k];

        t_901[k] = f_0 * mh_901[k];

        t_902[k] = f_0 * mh_902[k];

        t_903[k] = f_0 * mh_903[k];

        t_904[k] = f_0 * mh_904[k];

        t_905[k] = f_0 * mh_905[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, t_909, t_910, t_911, t_912, t_913, mh_906, \
                         mh_907, mh_908, mh_909, mh_910, mh_911, mh_912, \
                         mh_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_0 * mh_906[k];

        t_907[k] = f_0 * mh_907[k];

        t_908[k] = f_0 * mh_908[k];

        t_909[k] = f_0 * mh_909[k];

        t_910[k] = f_0 * mh_910[k];

        t_911[k] = f_0 * mh_911[k];

        t_912[k] = f_0 * mh_912[k];

        t_913[k] = f_0 * mh_913[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, t_918, t_919, t_920, t_921, mh_914, \
                         mh_915, mh_916, mh_917, mh_918, mh_919, mh_920, \
                         mh_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_0 * mh_914[k];

        t_915[k] = f_0 * mh_915[k];

        t_916[k] = f_0 * mh_916[k];

        t_917[k] = f_0 * mh_917[k];

        t_918[k] = f_0 * mh_918[k];

        t_919[k] = f_0 * mh_919[k];

        t_920[k] = f_0 * mh_920[k];

        t_921[k] = f_0 * mh_921[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, t_926, t_927, t_928, t_929, mh_922, \
                         mh_923, mh_924, mh_925, mh_926, mh_927, mh_928, \
                         mh_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_0 * mh_922[k];

        t_923[k] = f_0 * mh_923[k];

        t_924[k] = f_0 * mh_924[k];

        t_925[k] = f_0 * mh_925[k];

        t_926[k] = f_0 * mh_926[k];

        t_927[k] = f_0 * mh_927[k];

        t_928[k] = f_0 * mh_928[k];

        t_929[k] = f_0 * mh_929[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, t_935, t_936, t_937, mh_930, \
                         mh_931, mh_932, mh_933, mh_934, mh_935, mh_936, \
                         mh_937 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_0 * mh_930[k];

        t_931[k] = f_0 * mh_931[k];

        t_932[k] = f_0 * mh_932[k];

        t_933[k] = f_0 * mh_933[k];

        t_934[k] = f_0 * mh_934[k];

        t_935[k] = f_0 * mh_935[k];

        t_936[k] = f_0 * mh_936[k];

        t_937[k] = f_0 * mh_937[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, t_942, t_943, t_944, mh_938, mh_939, \
                         mh_940, mh_941, mh_942, mh_943, mh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_0 * mh_938[k];

        t_939[k] = f_0 * mh_939[k];

        t_940[k] = f_0 * mh_940[k];

        t_941[k] = f_0 * mh_941[k];

        t_942[k] = f_0 * mh_942[k];

        t_943[k] = f_0 * mh_943[k];

        t_944[k] = f_0 * mh_944[k];
    }
}

auto
compute_prim_geom_10_lh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t kh, const size_t mh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_lh_electron_repulsion_0_piece0(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_0_piece1(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_0_piece2(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_0_piece3(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_0_piece4(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_0_piece5(buffer, target, kh, mh, ncols, alpha);
}

static auto
compute_prim_geom_10_lh_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
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
    auto *t_135 = buffer.data(target + 135);
    auto *t_136 = buffer.data(target + 136);
    auto *t_137 = buffer.data(target + 137);
    auto *t_138 = buffer.data(target + 138);
    auto *t_139 = buffer.data(target + 139);
    auto *t_140 = buffer.data(target + 140);
    auto *t_141 = buffer.data(target + 141);
    auto *t_142 = buffer.data(target + 142);
    auto *t_143 = buffer.data(target + 143);
    auto *t_144 = buffer.data(target + 144);
    auto *t_145 = buffer.data(target + 145);
    auto *t_146 = buffer.data(target + 146);
    auto *t_147 = buffer.data(target + 147);
    auto *t_148 = buffer.data(target + 148);
    auto *t_149 = buffer.data(target + 149);
    auto *t_150 = buffer.data(target + 150);
    auto *t_151 = buffer.data(target + 151);
    auto *t_152 = buffer.data(target + 152);
    auto *t_153 = buffer.data(target + 153);
    auto *t_154 = buffer.data(target + 154);
    auto *t_155 = buffer.data(target + 155);
    auto *t_156 = buffer.data(target + 156);
    auto *t_157 = buffer.data(target + 157);
    auto *t_158 = buffer.data(target + 158);
    auto *t_159 = buffer.data(target + 159);
    auto *t_160 = buffer.data(target + 160);
    auto *t_161 = buffer.data(target + 161);
    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);
    auto *t_168 = buffer.data(target + 168);

    const auto *kh_0 = buffer.data(kh + 0);
    const auto *kh_1 = buffer.data(kh + 1);
    const auto *kh_2 = buffer.data(kh + 2);
    const auto *kh_3 = buffer.data(kh + 3);
    const auto *kh_4 = buffer.data(kh + 4);
    const auto *kh_5 = buffer.data(kh + 5);
    const auto *kh_6 = buffer.data(kh + 6);
    const auto *kh_7 = buffer.data(kh + 7);
    const auto *kh_8 = buffer.data(kh + 8);
    const auto *kh_9 = buffer.data(kh + 9);
    const auto *kh_10 = buffer.data(kh + 10);
    const auto *kh_11 = buffer.data(kh + 11);
    const auto *kh_12 = buffer.data(kh + 12);
    const auto *kh_13 = buffer.data(kh + 13);
    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_15 = buffer.data(kh + 15);
    const auto *kh_16 = buffer.data(kh + 16);
    const auto *kh_17 = buffer.data(kh + 17);
    const auto *kh_18 = buffer.data(kh + 18);
    const auto *kh_19 = buffer.data(kh + 19);
    const auto *kh_20 = buffer.data(kh + 20);
    const auto *kh_21 = buffer.data(kh + 21);
    const auto *kh_22 = buffer.data(kh + 22);
    const auto *kh_23 = buffer.data(kh + 23);
    const auto *kh_24 = buffer.data(kh + 24);
    const auto *kh_25 = buffer.data(kh + 25);
    const auto *kh_26 = buffer.data(kh + 26);
    const auto *kh_27 = buffer.data(kh + 27);
    const auto *kh_28 = buffer.data(kh + 28);
    const auto *kh_29 = buffer.data(kh + 29);
    const auto *kh_30 = buffer.data(kh + 30);
    const auto *kh_31 = buffer.data(kh + 31);
    const auto *kh_32 = buffer.data(kh + 32);
    const auto *kh_33 = buffer.data(kh + 33);
    const auto *kh_34 = buffer.data(kh + 34);
    const auto *kh_35 = buffer.data(kh + 35);
    const auto *kh_36 = buffer.data(kh + 36);
    const auto *kh_37 = buffer.data(kh + 37);
    const auto *kh_38 = buffer.data(kh + 38);
    const auto *kh_39 = buffer.data(kh + 39);
    const auto *kh_40 = buffer.data(kh + 40);
    const auto *kh_41 = buffer.data(kh + 41);
    const auto *kh_42 = buffer.data(kh + 42);
    const auto *kh_43 = buffer.data(kh + 43);
    const auto *kh_44 = buffer.data(kh + 44);
    const auto *kh_45 = buffer.data(kh + 45);
    const auto *kh_46 = buffer.data(kh + 46);
    const auto *kh_47 = buffer.data(kh + 47);
    const auto *kh_48 = buffer.data(kh + 48);
    const auto *kh_49 = buffer.data(kh + 49);
    const auto *kh_50 = buffer.data(kh + 50);
    const auto *kh_51 = buffer.data(kh + 51);
    const auto *kh_52 = buffer.data(kh + 52);
    const auto *kh_53 = buffer.data(kh + 53);
    const auto *kh_54 = buffer.data(kh + 54);
    const auto *kh_55 = buffer.data(kh + 55);
    const auto *kh_56 = buffer.data(kh + 56);
    const auto *kh_57 = buffer.data(kh + 57);
    const auto *kh_58 = buffer.data(kh + 58);
    const auto *kh_59 = buffer.data(kh + 59);
    const auto *kh_60 = buffer.data(kh + 60);
    const auto *kh_61 = buffer.data(kh + 61);
    const auto *kh_62 = buffer.data(kh + 62);
    const auto *kh_63 = buffer.data(kh + 63);
    const auto *kh_64 = buffer.data(kh + 64);
    const auto *kh_65 = buffer.data(kh + 65);
    const auto *kh_66 = buffer.data(kh + 66);
    const auto *kh_67 = buffer.data(kh + 67);
    const auto *kh_68 = buffer.data(kh + 68);
    const auto *kh_69 = buffer.data(kh + 69);
    const auto *kh_70 = buffer.data(kh + 70);
    const auto *kh_71 = buffer.data(kh + 71);
    const auto *kh_72 = buffer.data(kh + 72);
    const auto *kh_73 = buffer.data(kh + 73);
    const auto *kh_74 = buffer.data(kh + 74);
    const auto *kh_75 = buffer.data(kh + 75);
    const auto *kh_76 = buffer.data(kh + 76);
    const auto *kh_77 = buffer.data(kh + 77);
    const auto *kh_78 = buffer.data(kh + 78);
    const auto *kh_79 = buffer.data(kh + 79);
    const auto *kh_80 = buffer.data(kh + 80);
    const auto *kh_81 = buffer.data(kh + 81);
    const auto *kh_82 = buffer.data(kh + 82);
    const auto *kh_83 = buffer.data(kh + 83);
    const auto *kh_84 = buffer.data(kh + 84);
    const auto *kh_85 = buffer.data(kh + 85);
    const auto *kh_86 = buffer.data(kh + 86);
    const auto *kh_87 = buffer.data(kh + 87);
    const auto *kh_88 = buffer.data(kh + 88);
    const auto *kh_89 = buffer.data(kh + 89);
    const auto *kh_90 = buffer.data(kh + 90);
    const auto *kh_91 = buffer.data(kh + 91);
    const auto *kh_92 = buffer.data(kh + 92);
    const auto *kh_93 = buffer.data(kh + 93);
    const auto *kh_94 = buffer.data(kh + 94);
    const auto *kh_95 = buffer.data(kh + 95);
    const auto *kh_96 = buffer.data(kh + 96);
    const auto *kh_97 = buffer.data(kh + 97);
    const auto *kh_98 = buffer.data(kh + 98);
    const auto *kh_99 = buffer.data(kh + 99);
    const auto *kh_100 = buffer.data(kh + 100);
    const auto *kh_101 = buffer.data(kh + 101);
    const auto *kh_102 = buffer.data(kh + 102);
    const auto *kh_103 = buffer.data(kh + 103);
    const auto *kh_104 = buffer.data(kh + 104);
    const auto *kh_105 = buffer.data(kh + 105);

    const auto *mh_21 = buffer.data(mh + 21);
    const auto *mh_22 = buffer.data(mh + 22);
    const auto *mh_23 = buffer.data(mh + 23);
    const auto *mh_24 = buffer.data(mh + 24);
    const auto *mh_25 = buffer.data(mh + 25);
    const auto *mh_26 = buffer.data(mh + 26);
    const auto *mh_27 = buffer.data(mh + 27);
    const auto *mh_28 = buffer.data(mh + 28);
    const auto *mh_29 = buffer.data(mh + 29);
    const auto *mh_30 = buffer.data(mh + 30);
    const auto *mh_31 = buffer.data(mh + 31);
    const auto *mh_32 = buffer.data(mh + 32);
    const auto *mh_33 = buffer.data(mh + 33);
    const auto *mh_34 = buffer.data(mh + 34);
    const auto *mh_35 = buffer.data(mh + 35);
    const auto *mh_36 = buffer.data(mh + 36);
    const auto *mh_37 = buffer.data(mh + 37);
    const auto *mh_38 = buffer.data(mh + 38);
    const auto *mh_39 = buffer.data(mh + 39);
    const auto *mh_40 = buffer.data(mh + 40);
    const auto *mh_41 = buffer.data(mh + 41);
    const auto *mh_63 = buffer.data(mh + 63);
    const auto *mh_64 = buffer.data(mh + 64);
    const auto *mh_65 = buffer.data(mh + 65);
    const auto *mh_66 = buffer.data(mh + 66);
    const auto *mh_67 = buffer.data(mh + 67);
    const auto *mh_68 = buffer.data(mh + 68);
    const auto *mh_69 = buffer.data(mh + 69);
    const auto *mh_70 = buffer.data(mh + 70);
    const auto *mh_71 = buffer.data(mh + 71);
    const auto *mh_72 = buffer.data(mh + 72);
    const auto *mh_73 = buffer.data(mh + 73);
    const auto *mh_74 = buffer.data(mh + 74);
    const auto *mh_75 = buffer.data(mh + 75);
    const auto *mh_76 = buffer.data(mh + 76);
    const auto *mh_77 = buffer.data(mh + 77);
    const auto *mh_78 = buffer.data(mh + 78);
    const auto *mh_79 = buffer.data(mh + 79);
    const auto *mh_80 = buffer.data(mh + 80);
    const auto *mh_81 = buffer.data(mh + 81);
    const auto *mh_82 = buffer.data(mh + 82);
    const auto *mh_83 = buffer.data(mh + 83);
    const auto *mh_84 = buffer.data(mh + 84);
    const auto *mh_85 = buffer.data(mh + 85);
    const auto *mh_86 = buffer.data(mh + 86);
    const auto *mh_87 = buffer.data(mh + 87);
    const auto *mh_88 = buffer.data(mh + 88);
    const auto *mh_89 = buffer.data(mh + 89);
    const auto *mh_90 = buffer.data(mh + 90);
    const auto *mh_91 = buffer.data(mh + 91);
    const auto *mh_92 = buffer.data(mh + 92);
    const auto *mh_93 = buffer.data(mh + 93);
    const auto *mh_94 = buffer.data(mh + 94);
    const auto *mh_95 = buffer.data(mh + 95);
    const auto *mh_96 = buffer.data(mh + 96);
    const auto *mh_97 = buffer.data(mh + 97);
    const auto *mh_98 = buffer.data(mh + 98);
    const auto *mh_99 = buffer.data(mh + 99);
    const auto *mh_100 = buffer.data(mh + 100);
    const auto *mh_101 = buffer.data(mh + 101);
    const auto *mh_102 = buffer.data(mh + 102);
    const auto *mh_103 = buffer.data(mh + 103);
    const auto *mh_104 = buffer.data(mh + 104);
    const auto *mh_126 = buffer.data(mh + 126);
    const auto *mh_127 = buffer.data(mh + 127);
    const auto *mh_128 = buffer.data(mh + 128);
    const auto *mh_129 = buffer.data(mh + 129);
    const auto *mh_130 = buffer.data(mh + 130);
    const auto *mh_131 = buffer.data(mh + 131);
    const auto *mh_132 = buffer.data(mh + 132);
    const auto *mh_133 = buffer.data(mh + 133);
    const auto *mh_134 = buffer.data(mh + 134);
    const auto *mh_135 = buffer.data(mh + 135);
    const auto *mh_136 = buffer.data(mh + 136);
    const auto *mh_137 = buffer.data(mh + 137);
    const auto *mh_138 = buffer.data(mh + 138);
    const auto *mh_139 = buffer.data(mh + 139);
    const auto *mh_140 = buffer.data(mh + 140);
    const auto *mh_141 = buffer.data(mh + 141);
    const auto *mh_142 = buffer.data(mh + 142);
    const auto *mh_143 = buffer.data(mh + 143);
    const auto *mh_144 = buffer.data(mh + 144);
    const auto *mh_145 = buffer.data(mh + 145);
    const auto *mh_146 = buffer.data(mh + 146);
    const auto *mh_147 = buffer.data(mh + 147);
    const auto *mh_148 = buffer.data(mh + 148);
    const auto *mh_149 = buffer.data(mh + 149);
    const auto *mh_150 = buffer.data(mh + 150);
    const auto *mh_151 = buffer.data(mh + 151);
    const auto *mh_152 = buffer.data(mh + 152);
    const auto *mh_153 = buffer.data(mh + 153);
    const auto *mh_154 = buffer.data(mh + 154);
    const auto *mh_155 = buffer.data(mh + 155);
    const auto *mh_156 = buffer.data(mh + 156);
    const auto *mh_157 = buffer.data(mh + 157);
    const auto *mh_158 = buffer.data(mh + 158);
    const auto *mh_159 = buffer.data(mh + 159);
    const auto *mh_160 = buffer.data(mh + 160);
    const auto *mh_161 = buffer.data(mh + 161);
    const auto *mh_162 = buffer.data(mh + 162);
    const auto *mh_163 = buffer.data(mh + 163);
    const auto *mh_164 = buffer.data(mh + 164);
    const auto *mh_165 = buffer.data(mh + 165);
    const auto *mh_166 = buffer.data(mh + 166);
    const auto *mh_167 = buffer.data(mh + 167);
    const auto *mh_168 = buffer.data(mh + 168);
    const auto *mh_169 = buffer.data(mh + 169);
    const auto *mh_170 = buffer.data(mh + 170);
    const auto *mh_171 = buffer.data(mh + 171);
    const auto *mh_172 = buffer.data(mh + 172);
    const auto *mh_173 = buffer.data(mh + 173);
    const auto *mh_174 = buffer.data(mh + 174);
    const auto *mh_175 = buffer.data(mh + 175);
    const auto *mh_176 = buffer.data(mh + 176);
    const auto *mh_177 = buffer.data(mh + 177);
    const auto *mh_178 = buffer.data(mh + 178);
    const auto *mh_179 = buffer.data(mh + 179);
    const auto *mh_180 = buffer.data(mh + 180);
    const auto *mh_181 = buffer.data(mh + 181);
    const auto *mh_182 = buffer.data(mh + 182);
    const auto *mh_183 = buffer.data(mh + 183);
    const auto *mh_184 = buffer.data(mh + 184);
    const auto *mh_185 = buffer.data(mh + 185);
    const auto *mh_186 = buffer.data(mh + 186);
    const auto *mh_187 = buffer.data(mh + 187);
    const auto *mh_188 = buffer.data(mh + 188);
    const auto *mh_210 = buffer.data(mh + 210);
    const auto *mh_211 = buffer.data(mh + 211);
    const auto *mh_212 = buffer.data(mh + 212);
    const auto *mh_213 = buffer.data(mh + 213);
    const auto *mh_214 = buffer.data(mh + 214);
    const auto *mh_215 = buffer.data(mh + 215);
    const auto *mh_216 = buffer.data(mh + 216);
    const auto *mh_217 = buffer.data(mh + 217);
    const auto *mh_218 = buffer.data(mh + 218);
    const auto *mh_219 = buffer.data(mh + 219);
    const auto *mh_220 = buffer.data(mh + 220);
    const auto *mh_221 = buffer.data(mh + 221);
    const auto *mh_222 = buffer.data(mh + 222);
    const auto *mh_223 = buffer.data(mh + 223);
    const auto *mh_224 = buffer.data(mh + 224);
    const auto *mh_225 = buffer.data(mh + 225);
    const auto *mh_226 = buffer.data(mh + 226);
    const auto *mh_227 = buffer.data(mh + 227);
    const auto *mh_228 = buffer.data(mh + 228);
    const auto *mh_229 = buffer.data(mh + 229);
    const auto *mh_230 = buffer.data(mh + 230);
    const auto *mh_231 = buffer.data(mh + 231);
    const auto *mh_232 = buffer.data(mh + 232);
    const auto *mh_233 = buffer.data(mh + 233);
    const auto *mh_234 = buffer.data(mh + 234);
    const auto *mh_235 = buffer.data(mh + 235);
    const auto *mh_236 = buffer.data(mh + 236);
    const auto *mh_237 = buffer.data(mh + 237);
    const auto *mh_238 = buffer.data(mh + 238);
    const auto *mh_239 = buffer.data(mh + 239);
    const auto *mh_240 = buffer.data(mh + 240);
    const auto *mh_241 = buffer.data(mh + 241);
    const auto *mh_242 = buffer.data(mh + 242);
    const auto *mh_243 = buffer.data(mh + 243);
    const auto *mh_244 = buffer.data(mh + 244);
    const auto *mh_245 = buffer.data(mh + 245);
    const auto *mh_246 = buffer.data(mh + 246);
    const auto *mh_247 = buffer.data(mh + 247);
    const auto *mh_248 = buffer.data(mh + 248);
    const auto *mh_249 = buffer.data(mh + 249);
    const auto *mh_250 = buffer.data(mh + 250);
    const auto *mh_251 = buffer.data(mh + 251);
    const auto *mh_252 = buffer.data(mh + 252);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, mh_21, mh_22, mh_23, mh_24, \
                         mh_25, mh_26, mh_27, mh_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * mh_21[k];

        t_1[k] = f_0 * mh_22[k];

        t_2[k] = f_0 * mh_23[k];

        t_3[k] = f_0 * mh_24[k];

        t_4[k] = f_0 * mh_25[k];

        t_5[k] = f_0 * mh_26[k];

        t_6[k] = f_0 * mh_27[k];

        t_7[k] = f_0 * mh_28[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, mh_29, mh_30, mh_31, \
                         mh_32, mh_33, mh_34, mh_35, mh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * mh_29[k];

        t_9[k] = f_0 * mh_30[k];

        t_10[k] = f_0 * mh_31[k];

        t_11[k] = f_0 * mh_32[k];

        t_12[k] = f_0 * mh_33[k];

        t_13[k] = f_0 * mh_34[k];

        t_14[k] = f_0 * mh_35[k];

        t_15[k] = f_0 * mh_36[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, kh_0, kh_1, mh_37, mh_38, \
                         mh_39, mh_40, mh_41, mh_63, mh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * mh_37[k];

        t_17[k] = f_0 * mh_38[k];

        t_18[k] = f_0 * mh_39[k];

        t_19[k] = f_0 * mh_40[k];

        t_20[k] = f_0 * mh_41[k];

        t_21[k] = -kh_0[k]
                  + f_0 * mh_63[k];

        t_22[k] = -kh_1[k]
                  + f_0 * mh_64[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, kh_2, kh_3, kh_4, kh_5, kh_6, mh_65, \
                         mh_66, mh_67, mh_68, mh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -kh_2[k]
                  + f_0 * mh_65[k];

        t_24[k] = -kh_3[k]
                  + f_0 * mh_66[k];

        t_25[k] = -kh_4[k]
                  + f_0 * mh_67[k];

        t_26[k] = -kh_5[k]
                  + f_0 * mh_68[k];

        t_27[k] = -kh_6[k]
                  + f_0 * mh_69[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, kh_7, kh_8, kh_9, kh_10, kh_11, mh_70, \
                         mh_71, mh_72, mh_73, mh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -kh_7[k]
                  + f_0 * mh_70[k];

        t_29[k] = -kh_8[k]
                  + f_0 * mh_71[k];

        t_30[k] = -kh_9[k]
                  + f_0 * mh_72[k];

        t_31[k] = -kh_10[k]
                  + f_0 * mh_73[k];

        t_32[k] = -kh_11[k]
                  + f_0 * mh_74[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, kh_12, kh_13, kh_14, kh_15, kh_16, \
                         mh_75, mh_76, mh_77, mh_78, mh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -kh_12[k]
                  + f_0 * mh_75[k];

        t_34[k] = -kh_13[k]
                  + f_0 * mh_76[k];

        t_35[k] = -kh_14[k]
                  + f_0 * mh_77[k];

        t_36[k] = -kh_15[k]
                  + f_0 * mh_78[k];

        t_37[k] = -kh_16[k]
                  + f_0 * mh_79[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, kh_17, kh_18, kh_19, kh_20, \
                         mh_80, mh_81, mh_82, mh_83, mh_84, mh_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -kh_17[k]
                  + f_0 * mh_80[k];

        t_39[k] = -kh_18[k]
                  + f_0 * mh_81[k];

        t_40[k] = -kh_19[k]
                  + f_0 * mh_82[k];

        t_41[k] = -kh_20[k]
                  + f_0 * mh_83[k];

        t_42[k] = f_0 * mh_84[k];

        t_43[k] = f_0 * mh_85[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, t_51, mh_86, mh_87, mh_88, \
                         mh_89, mh_90, mh_91, mh_92, mh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * mh_86[k];

        t_45[k] = f_0 * mh_87[k];

        t_46[k] = f_0 * mh_88[k];

        t_47[k] = f_0 * mh_89[k];

        t_48[k] = f_0 * mh_90[k];

        t_49[k] = f_0 * mh_91[k];

        t_50[k] = f_0 * mh_92[k];

        t_51[k] = f_0 * mh_93[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, t_58, t_59, mh_94, mh_95, mh_96, \
                         mh_97, mh_98, mh_99, mh_100, mh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * mh_94[k];

        t_53[k] = f_0 * mh_95[k];

        t_54[k] = f_0 * mh_96[k];

        t_55[k] = f_0 * mh_97[k];

        t_56[k] = f_0 * mh_98[k];

        t_57[k] = f_0 * mh_99[k];

        t_58[k] = f_0 * mh_100[k];

        t_59[k] = f_0 * mh_101[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, kh_21, kh_22, kh_23, mh_102, \
                         mh_103, mh_104, mh_126, mh_127, mh_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * mh_102[k];

        t_61[k] = f_0 * mh_103[k];

        t_62[k] = f_0 * mh_104[k];

        t_63[k] = -2.0 * kh_21[k]
                  + f_0 * mh_126[k];

        t_64[k] = -2.0 * kh_22[k]
                  + f_0 * mh_127[k];

        t_65[k] = -2.0 * kh_23[k]
                  + f_0 * mh_128[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, kh_24, kh_25, kh_26, kh_27, kh_28, \
                         mh_129, mh_130, mh_131, mh_132, mh_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -2.0 * kh_24[k]
                  + f_0 * mh_129[k];

        t_67[k] = -2.0 * kh_25[k]
                  + f_0 * mh_130[k];

        t_68[k] = -2.0 * kh_26[k]
                  + f_0 * mh_131[k];

        t_69[k] = -2.0 * kh_27[k]
                  + f_0 * mh_132[k];

        t_70[k] = -2.0 * kh_28[k]
                  + f_0 * mh_133[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, kh_29, kh_30, kh_31, kh_32, kh_33, \
                         mh_134, mh_135, mh_136, mh_137, mh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -2.0 * kh_29[k]
                  + f_0 * mh_134[k];

        t_72[k] = -2.0 * kh_30[k]
                  + f_0 * mh_135[k];

        t_73[k] = -2.0 * kh_31[k]
                  + f_0 * mh_136[k];

        t_74[k] = -2.0 * kh_32[k]
                  + f_0 * mh_137[k];

        t_75[k] = -2.0 * kh_33[k]
                  + f_0 * mh_138[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, kh_34, kh_35, kh_36, kh_37, kh_38, \
                         mh_139, mh_140, mh_141, mh_142, mh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -2.0 * kh_34[k]
                  + f_0 * mh_139[k];

        t_77[k] = -2.0 * kh_35[k]
                  + f_0 * mh_140[k];

        t_78[k] = -2.0 * kh_36[k]
                  + f_0 * mh_141[k];

        t_79[k] = -2.0 * kh_37[k]
                  + f_0 * mh_142[k];

        t_80[k] = -2.0 * kh_38[k]
                  + f_0 * mh_143[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, kh_39, kh_40, kh_41, kh_42, kh_43, \
                         mh_144, mh_145, mh_146, mh_147, mh_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -2.0 * kh_39[k]
                  + f_0 * mh_144[k];

        t_82[k] = -2.0 * kh_40[k]
                  + f_0 * mh_145[k];

        t_83[k] = -2.0 * kh_41[k]
                  + f_0 * mh_146[k];

        t_84[k] = -kh_42[k]
                  + f_0 * mh_147[k];

        t_85[k] = -kh_43[k]
                  + f_0 * mh_148[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, kh_44, kh_45, kh_46, kh_47, kh_48, \
                         mh_149, mh_150, mh_151, mh_152, mh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -kh_44[k]
                  + f_0 * mh_149[k];

        t_87[k] = -kh_45[k]
                  + f_0 * mh_150[k];

        t_88[k] = -kh_46[k]
                  + f_0 * mh_151[k];

        t_89[k] = -kh_47[k]
                  + f_0 * mh_152[k];

        t_90[k] = -kh_48[k]
                  + f_0 * mh_153[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, kh_49, kh_50, kh_51, kh_52, kh_53, \
                         mh_154, mh_155, mh_156, mh_157, mh_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = -kh_49[k]
                  + f_0 * mh_154[k];

        t_92[k] = -kh_50[k]
                  + f_0 * mh_155[k];

        t_93[k] = -kh_51[k]
                  + f_0 * mh_156[k];

        t_94[k] = -kh_52[k]
                  + f_0 * mh_157[k];

        t_95[k] = -kh_53[k]
                  + f_0 * mh_158[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, kh_54, kh_55, kh_56, kh_57, kh_58, \
                         mh_159, mh_160, mh_161, mh_162, mh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = -kh_54[k]
                  + f_0 * mh_159[k];

        t_97[k] = -kh_55[k]
                  + f_0 * mh_160[k];

        t_98[k] = -kh_56[k]
                  + f_0 * mh_161[k];

        t_99[k] = -kh_57[k]
                  + f_0 * mh_162[k];

        t_100[k] = -kh_58[k]
                   + f_0 * mh_163[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, kh_59, kh_60, kh_61, kh_62, \
                         mh_164, mh_165, mh_166, mh_167, mh_168, \
                         mh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = -kh_59[k]
                   + f_0 * mh_164[k];

        t_102[k] = -kh_60[k]
                   + f_0 * mh_165[k];

        t_103[k] = -kh_61[k]
                   + f_0 * mh_166[k];

        t_104[k] = -kh_62[k]
                   + f_0 * mh_167[k];

        t_105[k] = f_0 * mh_168[k];

        t_106[k] = f_0 * mh_169[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, t_114, mh_170, \
                         mh_171, mh_172, mh_173, mh_174, mh_175, mh_176, \
                         mh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * mh_170[k];

        t_108[k] = f_0 * mh_171[k];

        t_109[k] = f_0 * mh_172[k];

        t_110[k] = f_0 * mh_173[k];

        t_111[k] = f_0 * mh_174[k];

        t_112[k] = f_0 * mh_175[k];

        t_113[k] = f_0 * mh_176[k];

        t_114[k] = f_0 * mh_177[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, t_121, t_122, mh_178, \
                         mh_179, mh_180, mh_181, mh_182, mh_183, mh_184, \
                         mh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_0 * mh_178[k];

        t_116[k] = f_0 * mh_179[k];

        t_117[k] = f_0 * mh_180[k];

        t_118[k] = f_0 * mh_181[k];

        t_119[k] = f_0 * mh_182[k];

        t_120[k] = f_0 * mh_183[k];

        t_121[k] = f_0 * mh_184[k];

        t_122[k] = f_0 * mh_185[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, t_127, t_128, kh_63, kh_64, kh_65, \
                         mh_186, mh_187, mh_188, mh_210, mh_211, \
                         mh_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_0 * mh_186[k];

        t_124[k] = f_0 * mh_187[k];

        t_125[k] = f_0 * mh_188[k];

        t_126[k] = -3.0 * kh_63[k]
                   + f_0 * mh_210[k];

        t_127[k] = -3.0 * kh_64[k]
                   + f_0 * mh_211[k];

        t_128[k] = -3.0 * kh_65[k]
                   + f_0 * mh_212[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, kh_66, kh_67, kh_68, kh_69, kh_70, \
                         mh_213, mh_214, mh_215, mh_216, mh_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -3.0 * kh_66[k]
                   + f_0 * mh_213[k];

        t_130[k] = -3.0 * kh_67[k]
                   + f_0 * mh_214[k];

        t_131[k] = -3.0 * kh_68[k]
                   + f_0 * mh_215[k];

        t_132[k] = -3.0 * kh_69[k]
                   + f_0 * mh_216[k];

        t_133[k] = -3.0 * kh_70[k]
                   + f_0 * mh_217[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, kh_71, kh_72, kh_73, kh_74, kh_75, \
                         mh_218, mh_219, mh_220, mh_221, mh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -3.0 * kh_71[k]
                   + f_0 * mh_218[k];

        t_135[k] = -3.0 * kh_72[k]
                   + f_0 * mh_219[k];

        t_136[k] = -3.0 * kh_73[k]
                   + f_0 * mh_220[k];

        t_137[k] = -3.0 * kh_74[k]
                   + f_0 * mh_221[k];

        t_138[k] = -3.0 * kh_75[k]
                   + f_0 * mh_222[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, kh_76, kh_77, kh_78, kh_79, kh_80, \
                         mh_223, mh_224, mh_225, mh_226, mh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -3.0 * kh_76[k]
                   + f_0 * mh_223[k];

        t_140[k] = -3.0 * kh_77[k]
                   + f_0 * mh_224[k];

        t_141[k] = -3.0 * kh_78[k]
                   + f_0 * mh_225[k];

        t_142[k] = -3.0 * kh_79[k]
                   + f_0 * mh_226[k];

        t_143[k] = -3.0 * kh_80[k]
                   + f_0 * mh_227[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, kh_81, kh_82, kh_83, kh_84, kh_85, \
                         mh_228, mh_229, mh_230, mh_231, mh_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -3.0 * kh_81[k]
                   + f_0 * mh_228[k];

        t_145[k] = -3.0 * kh_82[k]
                   + f_0 * mh_229[k];

        t_146[k] = -3.0 * kh_83[k]
                   + f_0 * mh_230[k];

        t_147[k] = -2.0 * kh_84[k]
                   + f_0 * mh_231[k];

        t_148[k] = -2.0 * kh_85[k]
                   + f_0 * mh_232[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, kh_86, kh_87, kh_88, kh_89, kh_90, \
                         mh_233, mh_234, mh_235, mh_236, mh_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -2.0 * kh_86[k]
                   + f_0 * mh_233[k];

        t_150[k] = -2.0 * kh_87[k]
                   + f_0 * mh_234[k];

        t_151[k] = -2.0 * kh_88[k]
                   + f_0 * mh_235[k];

        t_152[k] = -2.0 * kh_89[k]
                   + f_0 * mh_236[k];

        t_153[k] = -2.0 * kh_90[k]
                   + f_0 * mh_237[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, kh_91, kh_92, kh_93, kh_94, kh_95, \
                         mh_238, mh_239, mh_240, mh_241, mh_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -2.0 * kh_91[k]
                   + f_0 * mh_238[k];

        t_155[k] = -2.0 * kh_92[k]
                   + f_0 * mh_239[k];

        t_156[k] = -2.0 * kh_93[k]
                   + f_0 * mh_240[k];

        t_157[k] = -2.0 * kh_94[k]
                   + f_0 * mh_241[k];

        t_158[k] = -2.0 * kh_95[k]
                   + f_0 * mh_242[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, kh_96, kh_97, kh_98, kh_99, \
                         kh_100, mh_243, mh_244, mh_245, mh_246, \
                         mh_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -2.0 * kh_96[k]
                   + f_0 * mh_243[k];

        t_160[k] = -2.0 * kh_97[k]
                   + f_0 * mh_244[k];

        t_161[k] = -2.0 * kh_98[k]
                   + f_0 * mh_245[k];

        t_162[k] = -2.0 * kh_99[k]
                   + f_0 * mh_246[k];

        t_163[k] = -2.0 * kh_100[k]
                   + f_0 * mh_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, kh_101, kh_102, kh_103, kh_104, \
                         kh_105, mh_248, mh_249, mh_250, mh_251, \
                         mh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -2.0 * kh_101[k]
                   + f_0 * mh_248[k];

        t_165[k] = -2.0 * kh_102[k]
                   + f_0 * mh_249[k];

        t_166[k] = -2.0 * kh_103[k]
                   + f_0 * mh_250[k];

        t_167[k] = -2.0 * kh_104[k]
                   + f_0 * mh_251[k];

        t_168[k] = -kh_105[k]
                   + f_0 * mh_252[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_169 = buffer.data(target + 169);
    auto *t_170 = buffer.data(target + 170);
    auto *t_171 = buffer.data(target + 171);
    auto *t_172 = buffer.data(target + 172);
    auto *t_173 = buffer.data(target + 173);
    auto *t_174 = buffer.data(target + 174);
    auto *t_175 = buffer.data(target + 175);
    auto *t_176 = buffer.data(target + 176);
    auto *t_177 = buffer.data(target + 177);
    auto *t_178 = buffer.data(target + 178);
    auto *t_179 = buffer.data(target + 179);
    auto *t_180 = buffer.data(target + 180);
    auto *t_181 = buffer.data(target + 181);
    auto *t_182 = buffer.data(target + 182);
    auto *t_183 = buffer.data(target + 183);
    auto *t_184 = buffer.data(target + 184);
    auto *t_185 = buffer.data(target + 185);
    auto *t_186 = buffer.data(target + 186);
    auto *t_187 = buffer.data(target + 187);
    auto *t_188 = buffer.data(target + 188);
    auto *t_189 = buffer.data(target + 189);
    auto *t_190 = buffer.data(target + 190);
    auto *t_191 = buffer.data(target + 191);
    auto *t_192 = buffer.data(target + 192);
    auto *t_193 = buffer.data(target + 193);
    auto *t_194 = buffer.data(target + 194);
    auto *t_195 = buffer.data(target + 195);
    auto *t_196 = buffer.data(target + 196);
    auto *t_197 = buffer.data(target + 197);
    auto *t_198 = buffer.data(target + 198);
    auto *t_199 = buffer.data(target + 199);
    auto *t_200 = buffer.data(target + 200);
    auto *t_201 = buffer.data(target + 201);
    auto *t_202 = buffer.data(target + 202);
    auto *t_203 = buffer.data(target + 203);
    auto *t_204 = buffer.data(target + 204);
    auto *t_205 = buffer.data(target + 205);
    auto *t_206 = buffer.data(target + 206);
    auto *t_207 = buffer.data(target + 207);
    auto *t_208 = buffer.data(target + 208);
    auto *t_209 = buffer.data(target + 209);
    auto *t_210 = buffer.data(target + 210);
    auto *t_211 = buffer.data(target + 211);
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);
    auto *t_215 = buffer.data(target + 215);
    auto *t_216 = buffer.data(target + 216);
    auto *t_217 = buffer.data(target + 217);
    auto *t_218 = buffer.data(target + 218);
    auto *t_219 = buffer.data(target + 219);
    auto *t_220 = buffer.data(target + 220);
    auto *t_221 = buffer.data(target + 221);
    auto *t_222 = buffer.data(target + 222);
    auto *t_223 = buffer.data(target + 223);
    auto *t_224 = buffer.data(target + 224);
    auto *t_225 = buffer.data(target + 225);
    auto *t_226 = buffer.data(target + 226);
    auto *t_227 = buffer.data(target + 227);
    auto *t_228 = buffer.data(target + 228);
    auto *t_229 = buffer.data(target + 229);
    auto *t_230 = buffer.data(target + 230);
    auto *t_231 = buffer.data(target + 231);
    auto *t_232 = buffer.data(target + 232);
    auto *t_233 = buffer.data(target + 233);
    auto *t_234 = buffer.data(target + 234);
    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);
    auto *t_239 = buffer.data(target + 239);
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);
    auto *t_262 = buffer.data(target + 262);
    auto *t_263 = buffer.data(target + 263);
    auto *t_264 = buffer.data(target + 264);
    auto *t_265 = buffer.data(target + 265);
    auto *t_266 = buffer.data(target + 266);
    auto *t_267 = buffer.data(target + 267);
    auto *t_268 = buffer.data(target + 268);
    auto *t_269 = buffer.data(target + 269);
    auto *t_270 = buffer.data(target + 270);
    auto *t_271 = buffer.data(target + 271);
    auto *t_272 = buffer.data(target + 272);
    auto *t_273 = buffer.data(target + 273);
    auto *t_274 = buffer.data(target + 274);
    auto *t_275 = buffer.data(target + 275);
    auto *t_276 = buffer.data(target + 276);
    auto *t_277 = buffer.data(target + 277);
    auto *t_278 = buffer.data(target + 278);
    auto *t_279 = buffer.data(target + 279);
    auto *t_280 = buffer.data(target + 280);
    auto *t_281 = buffer.data(target + 281);
    auto *t_282 = buffer.data(target + 282);
    auto *t_283 = buffer.data(target + 283);
    auto *t_284 = buffer.data(target + 284);
    auto *t_285 = buffer.data(target + 285);
    auto *t_286 = buffer.data(target + 286);
    auto *t_287 = buffer.data(target + 287);
    auto *t_288 = buffer.data(target + 288);
    auto *t_289 = buffer.data(target + 289);
    auto *t_290 = buffer.data(target + 290);
    auto *t_291 = buffer.data(target + 291);
    auto *t_292 = buffer.data(target + 292);
    auto *t_293 = buffer.data(target + 293);
    auto *t_294 = buffer.data(target + 294);
    auto *t_295 = buffer.data(target + 295);
    auto *t_296 = buffer.data(target + 296);
    auto *t_297 = buffer.data(target + 297);
    auto *t_298 = buffer.data(target + 298);
    auto *t_299 = buffer.data(target + 299);
    auto *t_300 = buffer.data(target + 300);
    auto *t_301 = buffer.data(target + 301);
    auto *t_302 = buffer.data(target + 302);
    auto *t_303 = buffer.data(target + 303);
    auto *t_304 = buffer.data(target + 304);
    auto *t_305 = buffer.data(target + 305);
    auto *t_306 = buffer.data(target + 306);
    auto *t_307 = buffer.data(target + 307);
    auto *t_308 = buffer.data(target + 308);
    auto *t_309 = buffer.data(target + 309);
    auto *t_310 = buffer.data(target + 310);
    auto *t_311 = buffer.data(target + 311);
    auto *t_312 = buffer.data(target + 312);
    auto *t_313 = buffer.data(target + 313);
    auto *t_314 = buffer.data(target + 314);
    auto *t_315 = buffer.data(target + 315);
    auto *t_316 = buffer.data(target + 316);
    auto *t_317 = buffer.data(target + 317);
    auto *t_318 = buffer.data(target + 318);
    auto *t_319 = buffer.data(target + 319);
    auto *t_320 = buffer.data(target + 320);
    auto *t_321 = buffer.data(target + 321);
    auto *t_322 = buffer.data(target + 322);
    auto *t_323 = buffer.data(target + 323);
    auto *t_324 = buffer.data(target + 324);
    auto *t_325 = buffer.data(target + 325);
    auto *t_326 = buffer.data(target + 326);
    auto *t_327 = buffer.data(target + 327);
    auto *t_328 = buffer.data(target + 328);
    auto *t_329 = buffer.data(target + 329);

    const auto *kh_106 = buffer.data(kh + 106);
    const auto *kh_107 = buffer.data(kh + 107);
    const auto *kh_108 = buffer.data(kh + 108);
    const auto *kh_109 = buffer.data(kh + 109);
    const auto *kh_110 = buffer.data(kh + 110);
    const auto *kh_111 = buffer.data(kh + 111);
    const auto *kh_112 = buffer.data(kh + 112);
    const auto *kh_113 = buffer.data(kh + 113);
    const auto *kh_114 = buffer.data(kh + 114);
    const auto *kh_115 = buffer.data(kh + 115);
    const auto *kh_116 = buffer.data(kh + 116);
    const auto *kh_117 = buffer.data(kh + 117);
    const auto *kh_118 = buffer.data(kh + 118);
    const auto *kh_119 = buffer.data(kh + 119);
    const auto *kh_120 = buffer.data(kh + 120);
    const auto *kh_121 = buffer.data(kh + 121);
    const auto *kh_122 = buffer.data(kh + 122);
    const auto *kh_123 = buffer.data(kh + 123);
    const auto *kh_124 = buffer.data(kh + 124);
    const auto *kh_125 = buffer.data(kh + 125);
    const auto *kh_126 = buffer.data(kh + 126);
    const auto *kh_127 = buffer.data(kh + 127);
    const auto *kh_128 = buffer.data(kh + 128);
    const auto *kh_129 = buffer.data(kh + 129);
    const auto *kh_130 = buffer.data(kh + 130);
    const auto *kh_131 = buffer.data(kh + 131);
    const auto *kh_132 = buffer.data(kh + 132);
    const auto *kh_133 = buffer.data(kh + 133);
    const auto *kh_134 = buffer.data(kh + 134);
    const auto *kh_135 = buffer.data(kh + 135);
    const auto *kh_136 = buffer.data(kh + 136);
    const auto *kh_137 = buffer.data(kh + 137);
    const auto *kh_138 = buffer.data(kh + 138);
    const auto *kh_139 = buffer.data(kh + 139);
    const auto *kh_140 = buffer.data(kh + 140);
    const auto *kh_141 = buffer.data(kh + 141);
    const auto *kh_142 = buffer.data(kh + 142);
    const auto *kh_143 = buffer.data(kh + 143);
    const auto *kh_144 = buffer.data(kh + 144);
    const auto *kh_145 = buffer.data(kh + 145);
    const auto *kh_146 = buffer.data(kh + 146);
    const auto *kh_147 = buffer.data(kh + 147);
    const auto *kh_148 = buffer.data(kh + 148);
    const auto *kh_149 = buffer.data(kh + 149);
    const auto *kh_150 = buffer.data(kh + 150);
    const auto *kh_151 = buffer.data(kh + 151);
    const auto *kh_152 = buffer.data(kh + 152);
    const auto *kh_153 = buffer.data(kh + 153);
    const auto *kh_154 = buffer.data(kh + 154);
    const auto *kh_155 = buffer.data(kh + 155);
    const auto *kh_156 = buffer.data(kh + 156);
    const auto *kh_157 = buffer.data(kh + 157);
    const auto *kh_158 = buffer.data(kh + 158);
    const auto *kh_159 = buffer.data(kh + 159);
    const auto *kh_160 = buffer.data(kh + 160);
    const auto *kh_161 = buffer.data(kh + 161);
    const auto *kh_162 = buffer.data(kh + 162);
    const auto *kh_163 = buffer.data(kh + 163);
    const auto *kh_164 = buffer.data(kh + 164);
    const auto *kh_165 = buffer.data(kh + 165);
    const auto *kh_166 = buffer.data(kh + 166);
    const auto *kh_167 = buffer.data(kh + 167);
    const auto *kh_168 = buffer.data(kh + 168);
    const auto *kh_169 = buffer.data(kh + 169);
    const auto *kh_170 = buffer.data(kh + 170);
    const auto *kh_171 = buffer.data(kh + 171);
    const auto *kh_172 = buffer.data(kh + 172);
    const auto *kh_173 = buffer.data(kh + 173);
    const auto *kh_174 = buffer.data(kh + 174);
    const auto *kh_175 = buffer.data(kh + 175);
    const auto *kh_176 = buffer.data(kh + 176);
    const auto *kh_177 = buffer.data(kh + 177);
    const auto *kh_178 = buffer.data(kh + 178);
    const auto *kh_179 = buffer.data(kh + 179);
    const auto *kh_180 = buffer.data(kh + 180);
    const auto *kh_181 = buffer.data(kh + 181);
    const auto *kh_182 = buffer.data(kh + 182);
    const auto *kh_183 = buffer.data(kh + 183);
    const auto *kh_184 = buffer.data(kh + 184);
    const auto *kh_185 = buffer.data(kh + 185);
    const auto *kh_186 = buffer.data(kh + 186);
    const auto *kh_187 = buffer.data(kh + 187);
    const auto *kh_188 = buffer.data(kh + 188);
    const auto *kh_189 = buffer.data(kh + 189);
    const auto *kh_190 = buffer.data(kh + 190);
    const auto *kh_191 = buffer.data(kh + 191);
    const auto *kh_192 = buffer.data(kh + 192);
    const auto *kh_193 = buffer.data(kh + 193);
    const auto *kh_194 = buffer.data(kh + 194);
    const auto *kh_195 = buffer.data(kh + 195);
    const auto *kh_196 = buffer.data(kh + 196);
    const auto *kh_197 = buffer.data(kh + 197);
    const auto *kh_198 = buffer.data(kh + 198);
    const auto *kh_199 = buffer.data(kh + 199);
    const auto *kh_200 = buffer.data(kh + 200);
    const auto *kh_201 = buffer.data(kh + 201);
    const auto *kh_202 = buffer.data(kh + 202);
    const auto *kh_203 = buffer.data(kh + 203);
    const auto *kh_204 = buffer.data(kh + 204);
    const auto *kh_205 = buffer.data(kh + 205);
    const auto *kh_206 = buffer.data(kh + 206);
    const auto *kh_207 = buffer.data(kh + 207);
    const auto *kh_208 = buffer.data(kh + 208);
    const auto *kh_209 = buffer.data(kh + 209);
    const auto *kh_210 = buffer.data(kh + 210);
    const auto *kh_211 = buffer.data(kh + 211);
    const auto *kh_212 = buffer.data(kh + 212);
    const auto *kh_213 = buffer.data(kh + 213);
    const auto *kh_214 = buffer.data(kh + 214);
    const auto *kh_215 = buffer.data(kh + 215);
    const auto *kh_216 = buffer.data(kh + 216);
    const auto *kh_217 = buffer.data(kh + 217);
    const auto *kh_218 = buffer.data(kh + 218);
    const auto *kh_219 = buffer.data(kh + 219);
    const auto *kh_220 = buffer.data(kh + 220);
    const auto *kh_221 = buffer.data(kh + 221);
    const auto *kh_222 = buffer.data(kh + 222);
    const auto *kh_223 = buffer.data(kh + 223);
    const auto *kh_224 = buffer.data(kh + 224);

    const auto *mh_253 = buffer.data(mh + 253);
    const auto *mh_254 = buffer.data(mh + 254);
    const auto *mh_255 = buffer.data(mh + 255);
    const auto *mh_256 = buffer.data(mh + 256);
    const auto *mh_257 = buffer.data(mh + 257);
    const auto *mh_258 = buffer.data(mh + 258);
    const auto *mh_259 = buffer.data(mh + 259);
    const auto *mh_260 = buffer.data(mh + 260);
    const auto *mh_261 = buffer.data(mh + 261);
    const auto *mh_262 = buffer.data(mh + 262);
    const auto *mh_263 = buffer.data(mh + 263);
    const auto *mh_264 = buffer.data(mh + 264);
    const auto *mh_265 = buffer.data(mh + 265);
    const auto *mh_266 = buffer.data(mh + 266);
    const auto *mh_267 = buffer.data(mh + 267);
    const auto *mh_268 = buffer.data(mh + 268);
    const auto *mh_269 = buffer.data(mh + 269);
    const auto *mh_270 = buffer.data(mh + 270);
    const auto *mh_271 = buffer.data(mh + 271);
    const auto *mh_272 = buffer.data(mh + 272);
    const auto *mh_273 = buffer.data(mh + 273);
    const auto *mh_274 = buffer.data(mh + 274);
    const auto *mh_275 = buffer.data(mh + 275);
    const auto *mh_276 = buffer.data(mh + 276);
    const auto *mh_277 = buffer.data(mh + 277);
    const auto *mh_278 = buffer.data(mh + 278);
    const auto *mh_279 = buffer.data(mh + 279);
    const auto *mh_280 = buffer.data(mh + 280);
    const auto *mh_281 = buffer.data(mh + 281);
    const auto *mh_282 = buffer.data(mh + 282);
    const auto *mh_283 = buffer.data(mh + 283);
    const auto *mh_284 = buffer.data(mh + 284);
    const auto *mh_285 = buffer.data(mh + 285);
    const auto *mh_286 = buffer.data(mh + 286);
    const auto *mh_287 = buffer.data(mh + 287);
    const auto *mh_288 = buffer.data(mh + 288);
    const auto *mh_289 = buffer.data(mh + 289);
    const auto *mh_290 = buffer.data(mh + 290);
    const auto *mh_291 = buffer.data(mh + 291);
    const auto *mh_292 = buffer.data(mh + 292);
    const auto *mh_293 = buffer.data(mh + 293);
    const auto *mh_315 = buffer.data(mh + 315);
    const auto *mh_316 = buffer.data(mh + 316);
    const auto *mh_317 = buffer.data(mh + 317);
    const auto *mh_318 = buffer.data(mh + 318);
    const auto *mh_319 = buffer.data(mh + 319);
    const auto *mh_320 = buffer.data(mh + 320);
    const auto *mh_321 = buffer.data(mh + 321);
    const auto *mh_322 = buffer.data(mh + 322);
    const auto *mh_323 = buffer.data(mh + 323);
    const auto *mh_324 = buffer.data(mh + 324);
    const auto *mh_325 = buffer.data(mh + 325);
    const auto *mh_326 = buffer.data(mh + 326);
    const auto *mh_327 = buffer.data(mh + 327);
    const auto *mh_328 = buffer.data(mh + 328);
    const auto *mh_329 = buffer.data(mh + 329);
    const auto *mh_330 = buffer.data(mh + 330);
    const auto *mh_331 = buffer.data(mh + 331);
    const auto *mh_332 = buffer.data(mh + 332);
    const auto *mh_333 = buffer.data(mh + 333);
    const auto *mh_334 = buffer.data(mh + 334);
    const auto *mh_335 = buffer.data(mh + 335);
    const auto *mh_336 = buffer.data(mh + 336);
    const auto *mh_337 = buffer.data(mh + 337);
    const auto *mh_338 = buffer.data(mh + 338);
    const auto *mh_339 = buffer.data(mh + 339);
    const auto *mh_340 = buffer.data(mh + 340);
    const auto *mh_341 = buffer.data(mh + 341);
    const auto *mh_342 = buffer.data(mh + 342);
    const auto *mh_343 = buffer.data(mh + 343);
    const auto *mh_344 = buffer.data(mh + 344);
    const auto *mh_345 = buffer.data(mh + 345);
    const auto *mh_346 = buffer.data(mh + 346);
    const auto *mh_347 = buffer.data(mh + 347);
    const auto *mh_348 = buffer.data(mh + 348);
    const auto *mh_349 = buffer.data(mh + 349);
    const auto *mh_350 = buffer.data(mh + 350);
    const auto *mh_351 = buffer.data(mh + 351);
    const auto *mh_352 = buffer.data(mh + 352);
    const auto *mh_353 = buffer.data(mh + 353);
    const auto *mh_354 = buffer.data(mh + 354);
    const auto *mh_355 = buffer.data(mh + 355);
    const auto *mh_356 = buffer.data(mh + 356);
    const auto *mh_357 = buffer.data(mh + 357);
    const auto *mh_358 = buffer.data(mh + 358);
    const auto *mh_359 = buffer.data(mh + 359);
    const auto *mh_360 = buffer.data(mh + 360);
    const auto *mh_361 = buffer.data(mh + 361);
    const auto *mh_362 = buffer.data(mh + 362);
    const auto *mh_363 = buffer.data(mh + 363);
    const auto *mh_364 = buffer.data(mh + 364);
    const auto *mh_365 = buffer.data(mh + 365);
    const auto *mh_366 = buffer.data(mh + 366);
    const auto *mh_367 = buffer.data(mh + 367);
    const auto *mh_368 = buffer.data(mh + 368);
    const auto *mh_369 = buffer.data(mh + 369);
    const auto *mh_370 = buffer.data(mh + 370);
    const auto *mh_371 = buffer.data(mh + 371);
    const auto *mh_372 = buffer.data(mh + 372);
    const auto *mh_373 = buffer.data(mh + 373);
    const auto *mh_374 = buffer.data(mh + 374);
    const auto *mh_375 = buffer.data(mh + 375);
    const auto *mh_376 = buffer.data(mh + 376);
    const auto *mh_377 = buffer.data(mh + 377);
    const auto *mh_378 = buffer.data(mh + 378);
    const auto *mh_379 = buffer.data(mh + 379);
    const auto *mh_380 = buffer.data(mh + 380);
    const auto *mh_381 = buffer.data(mh + 381);
    const auto *mh_382 = buffer.data(mh + 382);
    const auto *mh_383 = buffer.data(mh + 383);
    const auto *mh_384 = buffer.data(mh + 384);
    const auto *mh_385 = buffer.data(mh + 385);
    const auto *mh_386 = buffer.data(mh + 386);
    const auto *mh_387 = buffer.data(mh + 387);
    const auto *mh_388 = buffer.data(mh + 388);
    const auto *mh_389 = buffer.data(mh + 389);
    const auto *mh_390 = buffer.data(mh + 390);
    const auto *mh_391 = buffer.data(mh + 391);
    const auto *mh_392 = buffer.data(mh + 392);
    const auto *mh_393 = buffer.data(mh + 393);
    const auto *mh_394 = buffer.data(mh + 394);
    const auto *mh_395 = buffer.data(mh + 395);
    const auto *mh_396 = buffer.data(mh + 396);
    const auto *mh_397 = buffer.data(mh + 397);
    const auto *mh_398 = buffer.data(mh + 398);
    const auto *mh_399 = buffer.data(mh + 399);
    const auto *mh_400 = buffer.data(mh + 400);
    const auto *mh_401 = buffer.data(mh + 401);
    const auto *mh_402 = buffer.data(mh + 402);
    const auto *mh_403 = buffer.data(mh + 403);
    const auto *mh_404 = buffer.data(mh + 404);
    const auto *mh_405 = buffer.data(mh + 405);
    const auto *mh_406 = buffer.data(mh + 406);
    const auto *mh_407 = buffer.data(mh + 407);
    const auto *mh_408 = buffer.data(mh + 408);
    const auto *mh_409 = buffer.data(mh + 409);
    const auto *mh_410 = buffer.data(mh + 410);
    const auto *mh_411 = buffer.data(mh + 411);
    const auto *mh_412 = buffer.data(mh + 412);
    const auto *mh_413 = buffer.data(mh + 413);
    const auto *mh_414 = buffer.data(mh + 414);
    const auto *mh_415 = buffer.data(mh + 415);
    const auto *mh_416 = buffer.data(mh + 416);
    const auto *mh_417 = buffer.data(mh + 417);
    const auto *mh_418 = buffer.data(mh + 418);
    const auto *mh_419 = buffer.data(mh + 419);
    const auto *mh_441 = buffer.data(mh + 441);
    const auto *mh_442 = buffer.data(mh + 442);
    const auto *mh_443 = buffer.data(mh + 443);
    const auto *mh_444 = buffer.data(mh + 444);
    const auto *mh_445 = buffer.data(mh + 445);
    const auto *mh_446 = buffer.data(mh + 446);
    const auto *mh_447 = buffer.data(mh + 447);
    const auto *mh_448 = buffer.data(mh + 448);
    const auto *mh_449 = buffer.data(mh + 449);
    const auto *mh_450 = buffer.data(mh + 450);
    const auto *mh_451 = buffer.data(mh + 451);
    const auto *mh_452 = buffer.data(mh + 452);
    const auto *mh_453 = buffer.data(mh + 453);
    const auto *mh_454 = buffer.data(mh + 454);
    const auto *mh_455 = buffer.data(mh + 455);

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, kh_106, kh_107, kh_108, kh_109, \
                         kh_110, mh_253, mh_254, mh_255, mh_256, \
                         mh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -kh_106[k]
                   + f_0 * mh_253[k];

        t_170[k] = -kh_107[k]
                   + f_0 * mh_254[k];

        t_171[k] = -kh_108[k]
                   + f_0 * mh_255[k];

        t_172[k] = -kh_109[k]
                   + f_0 * mh_256[k];

        t_173[k] = -kh_110[k]
                   + f_0 * mh_257[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, kh_111, kh_112, kh_113, kh_114, \
                         kh_115, mh_258, mh_259, mh_260, mh_261, \
                         mh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -kh_111[k]
                   + f_0 * mh_258[k];

        t_175[k] = -kh_112[k]
                   + f_0 * mh_259[k];

        t_176[k] = -kh_113[k]
                   + f_0 * mh_260[k];

        t_177[k] = -kh_114[k]
                   + f_0 * mh_261[k];

        t_178[k] = -kh_115[k]
                   + f_0 * mh_262[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, kh_116, kh_117, kh_118, kh_119, \
                         kh_120, mh_263, mh_264, mh_265, mh_266, \
                         mh_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -kh_116[k]
                   + f_0 * mh_263[k];

        t_180[k] = -kh_117[k]
                   + f_0 * mh_264[k];

        t_181[k] = -kh_118[k]
                   + f_0 * mh_265[k];

        t_182[k] = -kh_119[k]
                   + f_0 * mh_266[k];

        t_183[k] = -kh_120[k]
                   + f_0 * mh_267[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, kh_121, kh_122, kh_123, kh_124, \
                         kh_125, mh_268, mh_269, mh_270, mh_271, \
                         mh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -kh_121[k]
                   + f_0 * mh_268[k];

        t_185[k] = -kh_122[k]
                   + f_0 * mh_269[k];

        t_186[k] = -kh_123[k]
                   + f_0 * mh_270[k];

        t_187[k] = -kh_124[k]
                   + f_0 * mh_271[k];

        t_188[k] = -kh_125[k]
                   + f_0 * mh_272[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, t_196, mh_273, \
                         mh_274, mh_275, mh_276, mh_277, mh_278, mh_279, \
                         mh_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * mh_273[k];

        t_190[k] = f_0 * mh_274[k];

        t_191[k] = f_0 * mh_275[k];

        t_192[k] = f_0 * mh_276[k];

        t_193[k] = f_0 * mh_277[k];

        t_194[k] = f_0 * mh_278[k];

        t_195[k] = f_0 * mh_279[k];

        t_196[k] = f_0 * mh_280[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, t_202, t_203, t_204, mh_281, \
                         mh_282, mh_283, mh_284, mh_285, mh_286, mh_287, \
                         mh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_0 * mh_281[k];

        t_198[k] = f_0 * mh_282[k];

        t_199[k] = f_0 * mh_283[k];

        t_200[k] = f_0 * mh_284[k];

        t_201[k] = f_0 * mh_285[k];

        t_202[k] = f_0 * mh_286[k];

        t_203[k] = f_0 * mh_287[k];

        t_204[k] = f_0 * mh_288[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, t_210, t_211, kh_126, kh_127, \
                         mh_289, mh_290, mh_291, mh_292, mh_293, mh_315, \
                         mh_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_0 * mh_289[k];

        t_206[k] = f_0 * mh_290[k];

        t_207[k] = f_0 * mh_291[k];

        t_208[k] = f_0 * mh_292[k];

        t_209[k] = f_0 * mh_293[k];

        t_210[k] = -4.0 * kh_126[k]
                   + f_0 * mh_315[k];

        t_211[k] = -4.0 * kh_127[k]
                   + f_0 * mh_316[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, kh_128, kh_129, kh_130, kh_131, \
                         kh_132, mh_317, mh_318, mh_319, mh_320, \
                         mh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = -4.0 * kh_128[k]
                   + f_0 * mh_317[k];

        t_213[k] = -4.0 * kh_129[k]
                   + f_0 * mh_318[k];

        t_214[k] = -4.0 * kh_130[k]
                   + f_0 * mh_319[k];

        t_215[k] = -4.0 * kh_131[k]
                   + f_0 * mh_320[k];

        t_216[k] = -4.0 * kh_132[k]
                   + f_0 * mh_321[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, kh_133, kh_134, kh_135, kh_136, \
                         kh_137, mh_322, mh_323, mh_324, mh_325, \
                         mh_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = -4.0 * kh_133[k]
                   + f_0 * mh_322[k];

        t_218[k] = -4.0 * kh_134[k]
                   + f_0 * mh_323[k];

        t_219[k] = -4.0 * kh_135[k]
                   + f_0 * mh_324[k];

        t_220[k] = -4.0 * kh_136[k]
                   + f_0 * mh_325[k];

        t_221[k] = -4.0 * kh_137[k]
                   + f_0 * mh_326[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, kh_138, kh_139, kh_140, kh_141, \
                         kh_142, mh_327, mh_328, mh_329, mh_330, \
                         mh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -4.0 * kh_138[k]
                   + f_0 * mh_327[k];

        t_223[k] = -4.0 * kh_139[k]
                   + f_0 * mh_328[k];

        t_224[k] = -4.0 * kh_140[k]
                   + f_0 * mh_329[k];

        t_225[k] = -4.0 * kh_141[k]
                   + f_0 * mh_330[k];

        t_226[k] = -4.0 * kh_142[k]
                   + f_0 * mh_331[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, t_231, kh_143, kh_144, kh_145, kh_146, \
                         kh_147, mh_332, mh_333, mh_334, mh_335, \
                         mh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = -4.0 * kh_143[k]
                   + f_0 * mh_332[k];

        t_228[k] = -4.0 * kh_144[k]
                   + f_0 * mh_333[k];

        t_229[k] = -4.0 * kh_145[k]
                   + f_0 * mh_334[k];

        t_230[k] = -4.0 * kh_146[k]
                   + f_0 * mh_335[k];

        t_231[k] = -3.0 * kh_147[k]
                   + f_0 * mh_336[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, kh_148, kh_149, kh_150, kh_151, \
                         kh_152, mh_337, mh_338, mh_339, mh_340, \
                         mh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = -3.0 * kh_148[k]
                   + f_0 * mh_337[k];

        t_233[k] = -3.0 * kh_149[k]
                   + f_0 * mh_338[k];

        t_234[k] = -3.0 * kh_150[k]
                   + f_0 * mh_339[k];

        t_235[k] = -3.0 * kh_151[k]
                   + f_0 * mh_340[k];

        t_236[k] = -3.0 * kh_152[k]
                   + f_0 * mh_341[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, kh_153, kh_154, kh_155, kh_156, \
                         kh_157, mh_342, mh_343, mh_344, mh_345, \
                         mh_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = -3.0 * kh_153[k]
                   + f_0 * mh_342[k];

        t_238[k] = -3.0 * kh_154[k]
                   + f_0 * mh_343[k];

        t_239[k] = -3.0 * kh_155[k]
                   + f_0 * mh_344[k];

        t_240[k] = -3.0 * kh_156[k]
                   + f_0 * mh_345[k];

        t_241[k] = -3.0 * kh_157[k]
                   + f_0 * mh_346[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, kh_158, kh_159, kh_160, kh_161, \
                         kh_162, mh_347, mh_348, mh_349, mh_350, \
                         mh_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = -3.0 * kh_158[k]
                   + f_0 * mh_347[k];

        t_243[k] = -3.0 * kh_159[k]
                   + f_0 * mh_348[k];

        t_244[k] = -3.0 * kh_160[k]
                   + f_0 * mh_349[k];

        t_245[k] = -3.0 * kh_161[k]
                   + f_0 * mh_350[k];

        t_246[k] = -3.0 * kh_162[k]
                   + f_0 * mh_351[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, kh_163, kh_164, kh_165, kh_166, \
                         kh_167, mh_352, mh_353, mh_354, mh_355, \
                         mh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = -3.0 * kh_163[k]
                   + f_0 * mh_352[k];

        t_248[k] = -3.0 * kh_164[k]
                   + f_0 * mh_353[k];

        t_249[k] = -3.0 * kh_165[k]
                   + f_0 * mh_354[k];

        t_250[k] = -3.0 * kh_166[k]
                   + f_0 * mh_355[k];

        t_251[k] = -3.0 * kh_167[k]
                   + f_0 * mh_356[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, kh_168, kh_169, kh_170, kh_171, \
                         kh_172, mh_357, mh_358, mh_359, mh_360, \
                         mh_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -2.0 * kh_168[k]
                   + f_0 * mh_357[k];

        t_253[k] = -2.0 * kh_169[k]
                   + f_0 * mh_358[k];

        t_254[k] = -2.0 * kh_170[k]
                   + f_0 * mh_359[k];

        t_255[k] = -2.0 * kh_171[k]
                   + f_0 * mh_360[k];

        t_256[k] = -2.0 * kh_172[k]
                   + f_0 * mh_361[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, kh_173, kh_174, kh_175, kh_176, \
                         kh_177, mh_362, mh_363, mh_364, mh_365, \
                         mh_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = -2.0 * kh_173[k]
                   + f_0 * mh_362[k];

        t_258[k] = -2.0 * kh_174[k]
                   + f_0 * mh_363[k];

        t_259[k] = -2.0 * kh_175[k]
                   + f_0 * mh_364[k];

        t_260[k] = -2.0 * kh_176[k]
                   + f_0 * mh_365[k];

        t_261[k] = -2.0 * kh_177[k]
                   + f_0 * mh_366[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, kh_178, kh_179, kh_180, kh_181, \
                         kh_182, mh_367, mh_368, mh_369, mh_370, \
                         mh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = -2.0 * kh_178[k]
                   + f_0 * mh_367[k];

        t_263[k] = -2.0 * kh_179[k]
                   + f_0 * mh_368[k];

        t_264[k] = -2.0 * kh_180[k]
                   + f_0 * mh_369[k];

        t_265[k] = -2.0 * kh_181[k]
                   + f_0 * mh_370[k];

        t_266[k] = -2.0 * kh_182[k]
                   + f_0 * mh_371[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, kh_183, kh_184, kh_185, kh_186, \
                         kh_187, mh_372, mh_373, mh_374, mh_375, \
                         mh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = -2.0 * kh_183[k]
                   + f_0 * mh_372[k];

        t_268[k] = -2.0 * kh_184[k]
                   + f_0 * mh_373[k];

        t_269[k] = -2.0 * kh_185[k]
                   + f_0 * mh_374[k];

        t_270[k] = -2.0 * kh_186[k]
                   + f_0 * mh_375[k];

        t_271[k] = -2.0 * kh_187[k]
                   + f_0 * mh_376[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, kh_188, kh_189, kh_190, kh_191, \
                         kh_192, mh_377, mh_378, mh_379, mh_380, \
                         mh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -2.0 * kh_188[k]
                   + f_0 * mh_377[k];

        t_273[k] = -kh_189[k]
                   + f_0 * mh_378[k];

        t_274[k] = -kh_190[k]
                   + f_0 * mh_379[k];

        t_275[k] = -kh_191[k]
                   + f_0 * mh_380[k];

        t_276[k] = -kh_192[k]
                   + f_0 * mh_381[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, t_281, kh_193, kh_194, kh_195, kh_196, \
                         kh_197, mh_382, mh_383, mh_384, mh_385, \
                         mh_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = -kh_193[k]
                   + f_0 * mh_382[k];

        t_278[k] = -kh_194[k]
                   + f_0 * mh_383[k];

        t_279[k] = -kh_195[k]
                   + f_0 * mh_384[k];

        t_280[k] = -kh_196[k]
                   + f_0 * mh_385[k];

        t_281[k] = -kh_197[k]
                   + f_0 * mh_386[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, kh_198, kh_199, kh_200, kh_201, \
                         kh_202, mh_387, mh_388, mh_389, mh_390, \
                         mh_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = -kh_198[k]
                   + f_0 * mh_387[k];

        t_283[k] = -kh_199[k]
                   + f_0 * mh_388[k];

        t_284[k] = -kh_200[k]
                   + f_0 * mh_389[k];

        t_285[k] = -kh_201[k]
                   + f_0 * mh_390[k];

        t_286[k] = -kh_202[k]
                   + f_0 * mh_391[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, kh_203, kh_204, kh_205, kh_206, \
                         kh_207, mh_392, mh_393, mh_394, mh_395, \
                         mh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = -kh_203[k]
                   + f_0 * mh_392[k];

        t_288[k] = -kh_204[k]
                   + f_0 * mh_393[k];

        t_289[k] = -kh_205[k]
                   + f_0 * mh_394[k];

        t_290[k] = -kh_206[k]
                   + f_0 * mh_395[k];

        t_291[k] = -kh_207[k]
                   + f_0 * mh_396[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, t_297, t_298, kh_208, kh_209, \
                         mh_397, mh_398, mh_399, mh_400, mh_401, mh_402, \
                         mh_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = -kh_208[k]
                   + f_0 * mh_397[k];

        t_293[k] = -kh_209[k]
                   + f_0 * mh_398[k];

        t_294[k] = f_0 * mh_399[k];

        t_295[k] = f_0 * mh_400[k];

        t_296[k] = f_0 * mh_401[k];

        t_297[k] = f_0 * mh_402[k];

        t_298[k] = f_0 * mh_403[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, t_304, t_305, t_306, mh_404, \
                         mh_405, mh_406, mh_407, mh_408, mh_409, mh_410, \
                         mh_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_0 * mh_404[k];

        t_300[k] = f_0 * mh_405[k];

        t_301[k] = f_0 * mh_406[k];

        t_302[k] = f_0 * mh_407[k];

        t_303[k] = f_0 * mh_408[k];

        t_304[k] = f_0 * mh_409[k];

        t_305[k] = f_0 * mh_410[k];

        t_306[k] = f_0 * mh_411[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, t_312, t_313, t_314, mh_412, \
                         mh_413, mh_414, mh_415, mh_416, mh_417, mh_418, \
                         mh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_0 * mh_412[k];

        t_308[k] = f_0 * mh_413[k];

        t_309[k] = f_0 * mh_414[k];

        t_310[k] = f_0 * mh_415[k];

        t_311[k] = f_0 * mh_416[k];

        t_312[k] = f_0 * mh_417[k];

        t_313[k] = f_0 * mh_418[k];

        t_314[k] = f_0 * mh_419[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, kh_210, kh_211, kh_212, kh_213, \
                         kh_214, mh_441, mh_442, mh_443, mh_444, \
                         mh_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -5.0 * kh_210[k]
                   + f_0 * mh_441[k];

        t_316[k] = -5.0 * kh_211[k]
                   + f_0 * mh_442[k];

        t_317[k] = -5.0 * kh_212[k]
                   + f_0 * mh_443[k];

        t_318[k] = -5.0 * kh_213[k]
                   + f_0 * mh_444[k];

        t_319[k] = -5.0 * kh_214[k]
                   + f_0 * mh_445[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, kh_215, kh_216, kh_217, kh_218, \
                         kh_219, mh_446, mh_447, mh_448, mh_449, \
                         mh_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -5.0 * kh_215[k]
                   + f_0 * mh_446[k];

        t_321[k] = -5.0 * kh_216[k]
                   + f_0 * mh_447[k];

        t_322[k] = -5.0 * kh_217[k]
                   + f_0 * mh_448[k];

        t_323[k] = -5.0 * kh_218[k]
                   + f_0 * mh_449[k];

        t_324[k] = -5.0 * kh_219[k]
                   + f_0 * mh_450[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, kh_220, kh_221, kh_222, kh_223, \
                         kh_224, mh_451, mh_452, mh_453, mh_454, \
                         mh_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -5.0 * kh_220[k]
                   + f_0 * mh_451[k];

        t_326[k] = -5.0 * kh_221[k]
                   + f_0 * mh_452[k];

        t_327[k] = -5.0 * kh_222[k]
                   + f_0 * mh_453[k];

        t_328[k] = -5.0 * kh_223[k]
                   + f_0 * mh_454[k];

        t_329[k] = -5.0 * kh_224[k]
                   + f_0 * mh_455[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_330 = buffer.data(target + 330);
    auto *t_331 = buffer.data(target + 331);
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);
    auto *t_334 = buffer.data(target + 334);
    auto *t_335 = buffer.data(target + 335);
    auto *t_336 = buffer.data(target + 336);
    auto *t_337 = buffer.data(target + 337);
    auto *t_338 = buffer.data(target + 338);
    auto *t_339 = buffer.data(target + 339);
    auto *t_340 = buffer.data(target + 340);
    auto *t_341 = buffer.data(target + 341);
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);
    auto *t_347 = buffer.data(target + 347);
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);
    auto *t_380 = buffer.data(target + 380);
    auto *t_381 = buffer.data(target + 381);
    auto *t_382 = buffer.data(target + 382);
    auto *t_383 = buffer.data(target + 383);
    auto *t_384 = buffer.data(target + 384);
    auto *t_385 = buffer.data(target + 385);
    auto *t_386 = buffer.data(target + 386);
    auto *t_387 = buffer.data(target + 387);
    auto *t_388 = buffer.data(target + 388);
    auto *t_389 = buffer.data(target + 389);
    auto *t_390 = buffer.data(target + 390);
    auto *t_391 = buffer.data(target + 391);
    auto *t_392 = buffer.data(target + 392);
    auto *t_393 = buffer.data(target + 393);
    auto *t_394 = buffer.data(target + 394);
    auto *t_395 = buffer.data(target + 395);
    auto *t_396 = buffer.data(target + 396);
    auto *t_397 = buffer.data(target + 397);
    auto *t_398 = buffer.data(target + 398);
    auto *t_399 = buffer.data(target + 399);
    auto *t_400 = buffer.data(target + 400);
    auto *t_401 = buffer.data(target + 401);
    auto *t_402 = buffer.data(target + 402);
    auto *t_403 = buffer.data(target + 403);
    auto *t_404 = buffer.data(target + 404);
    auto *t_405 = buffer.data(target + 405);
    auto *t_406 = buffer.data(target + 406);
    auto *t_407 = buffer.data(target + 407);
    auto *t_408 = buffer.data(target + 408);
    auto *t_409 = buffer.data(target + 409);
    auto *t_410 = buffer.data(target + 410);
    auto *t_411 = buffer.data(target + 411);
    auto *t_412 = buffer.data(target + 412);
    auto *t_413 = buffer.data(target + 413);
    auto *t_414 = buffer.data(target + 414);
    auto *t_415 = buffer.data(target + 415);
    auto *t_416 = buffer.data(target + 416);
    auto *t_417 = buffer.data(target + 417);
    auto *t_418 = buffer.data(target + 418);
    auto *t_419 = buffer.data(target + 419);
    auto *t_420 = buffer.data(target + 420);
    auto *t_421 = buffer.data(target + 421);
    auto *t_422 = buffer.data(target + 422);
    auto *t_423 = buffer.data(target + 423);
    auto *t_424 = buffer.data(target + 424);
    auto *t_425 = buffer.data(target + 425);
    auto *t_426 = buffer.data(target + 426);
    auto *t_427 = buffer.data(target + 427);
    auto *t_428 = buffer.data(target + 428);
    auto *t_429 = buffer.data(target + 429);
    auto *t_430 = buffer.data(target + 430);
    auto *t_431 = buffer.data(target + 431);
    auto *t_432 = buffer.data(target + 432);
    auto *t_433 = buffer.data(target + 433);
    auto *t_434 = buffer.data(target + 434);
    auto *t_435 = buffer.data(target + 435);
    auto *t_436 = buffer.data(target + 436);
    auto *t_437 = buffer.data(target + 437);
    auto *t_438 = buffer.data(target + 438);
    auto *t_439 = buffer.data(target + 439);
    auto *t_440 = buffer.data(target + 440);
    auto *t_441 = buffer.data(target + 441);
    auto *t_442 = buffer.data(target + 442);
    auto *t_443 = buffer.data(target + 443);
    auto *t_444 = buffer.data(target + 444);
    auto *t_445 = buffer.data(target + 445);
    auto *t_446 = buffer.data(target + 446);
    auto *t_447 = buffer.data(target + 447);
    auto *t_448 = buffer.data(target + 448);
    auto *t_449 = buffer.data(target + 449);
    auto *t_450 = buffer.data(target + 450);
    auto *t_451 = buffer.data(target + 451);
    auto *t_452 = buffer.data(target + 452);
    auto *t_453 = buffer.data(target + 453);
    auto *t_454 = buffer.data(target + 454);
    auto *t_455 = buffer.data(target + 455);
    auto *t_456 = buffer.data(target + 456);
    auto *t_457 = buffer.data(target + 457);
    auto *t_458 = buffer.data(target + 458);
    auto *t_459 = buffer.data(target + 459);
    auto *t_460 = buffer.data(target + 460);
    auto *t_461 = buffer.data(target + 461);
    auto *t_462 = buffer.data(target + 462);
    auto *t_463 = buffer.data(target + 463);
    auto *t_464 = buffer.data(target + 464);
    auto *t_465 = buffer.data(target + 465);
    auto *t_466 = buffer.data(target + 466);
    auto *t_467 = buffer.data(target + 467);
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);
    auto *t_474 = buffer.data(target + 474);
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);

    const auto *kh_225 = buffer.data(kh + 225);
    const auto *kh_226 = buffer.data(kh + 226);
    const auto *kh_227 = buffer.data(kh + 227);
    const auto *kh_228 = buffer.data(kh + 228);
    const auto *kh_229 = buffer.data(kh + 229);
    const auto *kh_230 = buffer.data(kh + 230);
    const auto *kh_231 = buffer.data(kh + 231);
    const auto *kh_232 = buffer.data(kh + 232);
    const auto *kh_233 = buffer.data(kh + 233);
    const auto *kh_234 = buffer.data(kh + 234);
    const auto *kh_235 = buffer.data(kh + 235);
    const auto *kh_236 = buffer.data(kh + 236);
    const auto *kh_237 = buffer.data(kh + 237);
    const auto *kh_238 = buffer.data(kh + 238);
    const auto *kh_239 = buffer.data(kh + 239);
    const auto *kh_240 = buffer.data(kh + 240);
    const auto *kh_241 = buffer.data(kh + 241);
    const auto *kh_242 = buffer.data(kh + 242);
    const auto *kh_243 = buffer.data(kh + 243);
    const auto *kh_244 = buffer.data(kh + 244);
    const auto *kh_245 = buffer.data(kh + 245);
    const auto *kh_246 = buffer.data(kh + 246);
    const auto *kh_247 = buffer.data(kh + 247);
    const auto *kh_248 = buffer.data(kh + 248);
    const auto *kh_249 = buffer.data(kh + 249);
    const auto *kh_250 = buffer.data(kh + 250);
    const auto *kh_251 = buffer.data(kh + 251);
    const auto *kh_252 = buffer.data(kh + 252);
    const auto *kh_253 = buffer.data(kh + 253);
    const auto *kh_254 = buffer.data(kh + 254);
    const auto *kh_255 = buffer.data(kh + 255);
    const auto *kh_256 = buffer.data(kh + 256);
    const auto *kh_257 = buffer.data(kh + 257);
    const auto *kh_258 = buffer.data(kh + 258);
    const auto *kh_259 = buffer.data(kh + 259);
    const auto *kh_260 = buffer.data(kh + 260);
    const auto *kh_261 = buffer.data(kh + 261);
    const auto *kh_262 = buffer.data(kh + 262);
    const auto *kh_263 = buffer.data(kh + 263);
    const auto *kh_264 = buffer.data(kh + 264);
    const auto *kh_265 = buffer.data(kh + 265);
    const auto *kh_266 = buffer.data(kh + 266);
    const auto *kh_267 = buffer.data(kh + 267);
    const auto *kh_268 = buffer.data(kh + 268);
    const auto *kh_269 = buffer.data(kh + 269);
    const auto *kh_270 = buffer.data(kh + 270);
    const auto *kh_271 = buffer.data(kh + 271);
    const auto *kh_272 = buffer.data(kh + 272);
    const auto *kh_273 = buffer.data(kh + 273);
    const auto *kh_274 = buffer.data(kh + 274);
    const auto *kh_275 = buffer.data(kh + 275);
    const auto *kh_276 = buffer.data(kh + 276);
    const auto *kh_277 = buffer.data(kh + 277);
    const auto *kh_278 = buffer.data(kh + 278);
    const auto *kh_279 = buffer.data(kh + 279);
    const auto *kh_280 = buffer.data(kh + 280);
    const auto *kh_281 = buffer.data(kh + 281);
    const auto *kh_282 = buffer.data(kh + 282);
    const auto *kh_283 = buffer.data(kh + 283);
    const auto *kh_284 = buffer.data(kh + 284);
    const auto *kh_285 = buffer.data(kh + 285);
    const auto *kh_286 = buffer.data(kh + 286);
    const auto *kh_287 = buffer.data(kh + 287);
    const auto *kh_288 = buffer.data(kh + 288);
    const auto *kh_289 = buffer.data(kh + 289);
    const auto *kh_290 = buffer.data(kh + 290);
    const auto *kh_291 = buffer.data(kh + 291);
    const auto *kh_292 = buffer.data(kh + 292);
    const auto *kh_293 = buffer.data(kh + 293);
    const auto *kh_294 = buffer.data(kh + 294);
    const auto *kh_295 = buffer.data(kh + 295);
    const auto *kh_296 = buffer.data(kh + 296);
    const auto *kh_297 = buffer.data(kh + 297);
    const auto *kh_298 = buffer.data(kh + 298);
    const auto *kh_299 = buffer.data(kh + 299);
    const auto *kh_300 = buffer.data(kh + 300);
    const auto *kh_301 = buffer.data(kh + 301);
    const auto *kh_302 = buffer.data(kh + 302);
    const auto *kh_303 = buffer.data(kh + 303);
    const auto *kh_304 = buffer.data(kh + 304);
    const auto *kh_305 = buffer.data(kh + 305);
    const auto *kh_306 = buffer.data(kh + 306);
    const auto *kh_307 = buffer.data(kh + 307);
    const auto *kh_308 = buffer.data(kh + 308);
    const auto *kh_309 = buffer.data(kh + 309);
    const auto *kh_310 = buffer.data(kh + 310);
    const auto *kh_311 = buffer.data(kh + 311);
    const auto *kh_312 = buffer.data(kh + 312);
    const auto *kh_313 = buffer.data(kh + 313);
    const auto *kh_314 = buffer.data(kh + 314);
    const auto *kh_315 = buffer.data(kh + 315);
    const auto *kh_316 = buffer.data(kh + 316);
    const auto *kh_317 = buffer.data(kh + 317);
    const auto *kh_318 = buffer.data(kh + 318);
    const auto *kh_319 = buffer.data(kh + 319);
    const auto *kh_320 = buffer.data(kh + 320);
    const auto *kh_321 = buffer.data(kh + 321);
    const auto *kh_322 = buffer.data(kh + 322);
    const auto *kh_323 = buffer.data(kh + 323);
    const auto *kh_324 = buffer.data(kh + 324);
    const auto *kh_325 = buffer.data(kh + 325);
    const auto *kh_326 = buffer.data(kh + 326);
    const auto *kh_327 = buffer.data(kh + 327);
    const auto *kh_328 = buffer.data(kh + 328);
    const auto *kh_329 = buffer.data(kh + 329);
    const auto *kh_330 = buffer.data(kh + 330);
    const auto *kh_331 = buffer.data(kh + 331);
    const auto *kh_332 = buffer.data(kh + 332);
    const auto *kh_333 = buffer.data(kh + 333);
    const auto *kh_334 = buffer.data(kh + 334);
    const auto *kh_335 = buffer.data(kh + 335);
    const auto *kh_336 = buffer.data(kh + 336);
    const auto *kh_337 = buffer.data(kh + 337);
    const auto *kh_338 = buffer.data(kh + 338);
    const auto *kh_339 = buffer.data(kh + 339);
    const auto *kh_340 = buffer.data(kh + 340);
    const auto *kh_341 = buffer.data(kh + 341);
    const auto *kh_342 = buffer.data(kh + 342);
    const auto *kh_343 = buffer.data(kh + 343);
    const auto *kh_344 = buffer.data(kh + 344);
    const auto *kh_345 = buffer.data(kh + 345);
    const auto *kh_346 = buffer.data(kh + 346);
    const auto *kh_347 = buffer.data(kh + 347);
    const auto *kh_348 = buffer.data(kh + 348);
    const auto *kh_349 = buffer.data(kh + 349);
    const auto *kh_350 = buffer.data(kh + 350);
    const auto *kh_351 = buffer.data(kh + 351);
    const auto *kh_352 = buffer.data(kh + 352);
    const auto *kh_353 = buffer.data(kh + 353);
    const auto *kh_354 = buffer.data(kh + 354);
    const auto *kh_355 = buffer.data(kh + 355);
    const auto *kh_356 = buffer.data(kh + 356);
    const auto *kh_357 = buffer.data(kh + 357);
    const auto *kh_358 = buffer.data(kh + 358);
    const auto *kh_359 = buffer.data(kh + 359);
    const auto *kh_360 = buffer.data(kh + 360);
    const auto *kh_361 = buffer.data(kh + 361);

    const auto *mh_456 = buffer.data(mh + 456);
    const auto *mh_457 = buffer.data(mh + 457);
    const auto *mh_458 = buffer.data(mh + 458);
    const auto *mh_459 = buffer.data(mh + 459);
    const auto *mh_460 = buffer.data(mh + 460);
    const auto *mh_461 = buffer.data(mh + 461);
    const auto *mh_462 = buffer.data(mh + 462);
    const auto *mh_463 = buffer.data(mh + 463);
    const auto *mh_464 = buffer.data(mh + 464);
    const auto *mh_465 = buffer.data(mh + 465);
    const auto *mh_466 = buffer.data(mh + 466);
    const auto *mh_467 = buffer.data(mh + 467);
    const auto *mh_468 = buffer.data(mh + 468);
    const auto *mh_469 = buffer.data(mh + 469);
    const auto *mh_470 = buffer.data(mh + 470);
    const auto *mh_471 = buffer.data(mh + 471);
    const auto *mh_472 = buffer.data(mh + 472);
    const auto *mh_473 = buffer.data(mh + 473);
    const auto *mh_474 = buffer.data(mh + 474);
    const auto *mh_475 = buffer.data(mh + 475);
    const auto *mh_476 = buffer.data(mh + 476);
    const auto *mh_477 = buffer.data(mh + 477);
    const auto *mh_478 = buffer.data(mh + 478);
    const auto *mh_479 = buffer.data(mh + 479);
    const auto *mh_480 = buffer.data(mh + 480);
    const auto *mh_481 = buffer.data(mh + 481);
    const auto *mh_482 = buffer.data(mh + 482);
    const auto *mh_483 = buffer.data(mh + 483);
    const auto *mh_484 = buffer.data(mh + 484);
    const auto *mh_485 = buffer.data(mh + 485);
    const auto *mh_486 = buffer.data(mh + 486);
    const auto *mh_487 = buffer.data(mh + 487);
    const auto *mh_488 = buffer.data(mh + 488);
    const auto *mh_489 = buffer.data(mh + 489);
    const auto *mh_490 = buffer.data(mh + 490);
    const auto *mh_491 = buffer.data(mh + 491);
    const auto *mh_492 = buffer.data(mh + 492);
    const auto *mh_493 = buffer.data(mh + 493);
    const auto *mh_494 = buffer.data(mh + 494);
    const auto *mh_495 = buffer.data(mh + 495);
    const auto *mh_496 = buffer.data(mh + 496);
    const auto *mh_497 = buffer.data(mh + 497);
    const auto *mh_498 = buffer.data(mh + 498);
    const auto *mh_499 = buffer.data(mh + 499);
    const auto *mh_500 = buffer.data(mh + 500);
    const auto *mh_501 = buffer.data(mh + 501);
    const auto *mh_502 = buffer.data(mh + 502);
    const auto *mh_503 = buffer.data(mh + 503);
    const auto *mh_504 = buffer.data(mh + 504);
    const auto *mh_505 = buffer.data(mh + 505);
    const auto *mh_506 = buffer.data(mh + 506);
    const auto *mh_507 = buffer.data(mh + 507);
    const auto *mh_508 = buffer.data(mh + 508);
    const auto *mh_509 = buffer.data(mh + 509);
    const auto *mh_510 = buffer.data(mh + 510);
    const auto *mh_511 = buffer.data(mh + 511);
    const auto *mh_512 = buffer.data(mh + 512);
    const auto *mh_513 = buffer.data(mh + 513);
    const auto *mh_514 = buffer.data(mh + 514);
    const auto *mh_515 = buffer.data(mh + 515);
    const auto *mh_516 = buffer.data(mh + 516);
    const auto *mh_517 = buffer.data(mh + 517);
    const auto *mh_518 = buffer.data(mh + 518);
    const auto *mh_519 = buffer.data(mh + 519);
    const auto *mh_520 = buffer.data(mh + 520);
    const auto *mh_521 = buffer.data(mh + 521);
    const auto *mh_522 = buffer.data(mh + 522);
    const auto *mh_523 = buffer.data(mh + 523);
    const auto *mh_524 = buffer.data(mh + 524);
    const auto *mh_525 = buffer.data(mh + 525);
    const auto *mh_526 = buffer.data(mh + 526);
    const auto *mh_527 = buffer.data(mh + 527);
    const auto *mh_528 = buffer.data(mh + 528);
    const auto *mh_529 = buffer.data(mh + 529);
    const auto *mh_530 = buffer.data(mh + 530);
    const auto *mh_531 = buffer.data(mh + 531);
    const auto *mh_532 = buffer.data(mh + 532);
    const auto *mh_533 = buffer.data(mh + 533);
    const auto *mh_534 = buffer.data(mh + 534);
    const auto *mh_535 = buffer.data(mh + 535);
    const auto *mh_536 = buffer.data(mh + 536);
    const auto *mh_537 = buffer.data(mh + 537);
    const auto *mh_538 = buffer.data(mh + 538);
    const auto *mh_539 = buffer.data(mh + 539);
    const auto *mh_540 = buffer.data(mh + 540);
    const auto *mh_541 = buffer.data(mh + 541);
    const auto *mh_542 = buffer.data(mh + 542);
    const auto *mh_543 = buffer.data(mh + 543);
    const auto *mh_544 = buffer.data(mh + 544);
    const auto *mh_545 = buffer.data(mh + 545);
    const auto *mh_546 = buffer.data(mh + 546);
    const auto *mh_547 = buffer.data(mh + 547);
    const auto *mh_548 = buffer.data(mh + 548);
    const auto *mh_549 = buffer.data(mh + 549);
    const auto *mh_550 = buffer.data(mh + 550);
    const auto *mh_551 = buffer.data(mh + 551);
    const auto *mh_552 = buffer.data(mh + 552);
    const auto *mh_553 = buffer.data(mh + 553);
    const auto *mh_554 = buffer.data(mh + 554);
    const auto *mh_555 = buffer.data(mh + 555);
    const auto *mh_556 = buffer.data(mh + 556);
    const auto *mh_557 = buffer.data(mh + 557);
    const auto *mh_558 = buffer.data(mh + 558);
    const auto *mh_559 = buffer.data(mh + 559);
    const auto *mh_560 = buffer.data(mh + 560);
    const auto *mh_561 = buffer.data(mh + 561);
    const auto *mh_562 = buffer.data(mh + 562);
    const auto *mh_563 = buffer.data(mh + 563);
    const auto *mh_564 = buffer.data(mh + 564);
    const auto *mh_565 = buffer.data(mh + 565);
    const auto *mh_566 = buffer.data(mh + 566);
    const auto *mh_588 = buffer.data(mh + 588);
    const auto *mh_589 = buffer.data(mh + 589);
    const auto *mh_590 = buffer.data(mh + 590);
    const auto *mh_591 = buffer.data(mh + 591);
    const auto *mh_592 = buffer.data(mh + 592);
    const auto *mh_593 = buffer.data(mh + 593);
    const auto *mh_594 = buffer.data(mh + 594);
    const auto *mh_595 = buffer.data(mh + 595);
    const auto *mh_596 = buffer.data(mh + 596);
    const auto *mh_597 = buffer.data(mh + 597);
    const auto *mh_598 = buffer.data(mh + 598);
    const auto *mh_599 = buffer.data(mh + 599);
    const auto *mh_600 = buffer.data(mh + 600);
    const auto *mh_601 = buffer.data(mh + 601);
    const auto *mh_602 = buffer.data(mh + 602);
    const auto *mh_603 = buffer.data(mh + 603);
    const auto *mh_604 = buffer.data(mh + 604);
    const auto *mh_605 = buffer.data(mh + 605);
    const auto *mh_606 = buffer.data(mh + 606);
    const auto *mh_607 = buffer.data(mh + 607);
    const auto *mh_608 = buffer.data(mh + 608);
    const auto *mh_609 = buffer.data(mh + 609);
    const auto *mh_610 = buffer.data(mh + 610);
    const auto *mh_611 = buffer.data(mh + 611);
    const auto *mh_612 = buffer.data(mh + 612);
    const auto *mh_613 = buffer.data(mh + 613);
    const auto *mh_614 = buffer.data(mh + 614);
    const auto *mh_615 = buffer.data(mh + 615);
    const auto *mh_616 = buffer.data(mh + 616);
    const auto *mh_617 = buffer.data(mh + 617);
    const auto *mh_618 = buffer.data(mh + 618);
    const auto *mh_619 = buffer.data(mh + 619);
    const auto *mh_620 = buffer.data(mh + 620);
    const auto *mh_621 = buffer.data(mh + 621);
    const auto *mh_622 = buffer.data(mh + 622);
    const auto *mh_623 = buffer.data(mh + 623);
    const auto *mh_624 = buffer.data(mh + 624);
    const auto *mh_625 = buffer.data(mh + 625);
    const auto *mh_626 = buffer.data(mh + 626);
    const auto *mh_627 = buffer.data(mh + 627);
    const auto *mh_628 = buffer.data(mh + 628);
    const auto *mh_629 = buffer.data(mh + 629);
    const auto *mh_630 = buffer.data(mh + 630);
    const auto *mh_631 = buffer.data(mh + 631);
    const auto *mh_632 = buffer.data(mh + 632);
    const auto *mh_633 = buffer.data(mh + 633);
    const auto *mh_634 = buffer.data(mh + 634);

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, kh_225, kh_226, kh_227, kh_228, \
                         kh_229, mh_456, mh_457, mh_458, mh_459, \
                         mh_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -5.0 * kh_225[k]
                   + f_0 * mh_456[k];

        t_331[k] = -5.0 * kh_226[k]
                   + f_0 * mh_457[k];

        t_332[k] = -5.0 * kh_227[k]
                   + f_0 * mh_458[k];

        t_333[k] = -5.0 * kh_228[k]
                   + f_0 * mh_459[k];

        t_334[k] = -5.0 * kh_229[k]
                   + f_0 * mh_460[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, kh_230, kh_231, kh_232, kh_233, \
                         kh_234, mh_461, mh_462, mh_463, mh_464, \
                         mh_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -5.0 * kh_230[k]
                   + f_0 * mh_461[k];

        t_336[k] = -4.0 * kh_231[k]
                   + f_0 * mh_462[k];

        t_337[k] = -4.0 * kh_232[k]
                   + f_0 * mh_463[k];

        t_338[k] = -4.0 * kh_233[k]
                   + f_0 * mh_464[k];

        t_339[k] = -4.0 * kh_234[k]
                   + f_0 * mh_465[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, kh_235, kh_236, kh_237, kh_238, \
                         kh_239, mh_466, mh_467, mh_468, mh_469, \
                         mh_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -4.0 * kh_235[k]
                   + f_0 * mh_466[k];

        t_341[k] = -4.0 * kh_236[k]
                   + f_0 * mh_467[k];

        t_342[k] = -4.0 * kh_237[k]
                   + f_0 * mh_468[k];

        t_343[k] = -4.0 * kh_238[k]
                   + f_0 * mh_469[k];

        t_344[k] = -4.0 * kh_239[k]
                   + f_0 * mh_470[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, kh_240, kh_241, kh_242, kh_243, \
                         kh_244, mh_471, mh_472, mh_473, mh_474, \
                         mh_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -4.0 * kh_240[k]
                   + f_0 * mh_471[k];

        t_346[k] = -4.0 * kh_241[k]
                   + f_0 * mh_472[k];

        t_347[k] = -4.0 * kh_242[k]
                   + f_0 * mh_473[k];

        t_348[k] = -4.0 * kh_243[k]
                   + f_0 * mh_474[k];

        t_349[k] = -4.0 * kh_244[k]
                   + f_0 * mh_475[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, kh_245, kh_246, kh_247, kh_248, \
                         kh_249, mh_476, mh_477, mh_478, mh_479, \
                         mh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -4.0 * kh_245[k]
                   + f_0 * mh_476[k];

        t_351[k] = -4.0 * kh_246[k]
                   + f_0 * mh_477[k];

        t_352[k] = -4.0 * kh_247[k]
                   + f_0 * mh_478[k];

        t_353[k] = -4.0 * kh_248[k]
                   + f_0 * mh_479[k];

        t_354[k] = -4.0 * kh_249[k]
                   + f_0 * mh_480[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, kh_250, kh_251, kh_252, kh_253, \
                         kh_254, mh_481, mh_482, mh_483, mh_484, \
                         mh_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -4.0 * kh_250[k]
                   + f_0 * mh_481[k];

        t_356[k] = -4.0 * kh_251[k]
                   + f_0 * mh_482[k];

        t_357[k] = -3.0 * kh_252[k]
                   + f_0 * mh_483[k];

        t_358[k] = -3.0 * kh_253[k]
                   + f_0 * mh_484[k];

        t_359[k] = -3.0 * kh_254[k]
                   + f_0 * mh_485[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, kh_255, kh_256, kh_257, kh_258, \
                         kh_259, mh_486, mh_487, mh_488, mh_489, \
                         mh_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -3.0 * kh_255[k]
                   + f_0 * mh_486[k];

        t_361[k] = -3.0 * kh_256[k]
                   + f_0 * mh_487[k];

        t_362[k] = -3.0 * kh_257[k]
                   + f_0 * mh_488[k];

        t_363[k] = -3.0 * kh_258[k]
                   + f_0 * mh_489[k];

        t_364[k] = -3.0 * kh_259[k]
                   + f_0 * mh_490[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, kh_260, kh_261, kh_262, kh_263, \
                         kh_264, mh_491, mh_492, mh_493, mh_494, \
                         mh_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -3.0 * kh_260[k]
                   + f_0 * mh_491[k];

        t_366[k] = -3.0 * kh_261[k]
                   + f_0 * mh_492[k];

        t_367[k] = -3.0 * kh_262[k]
                   + f_0 * mh_493[k];

        t_368[k] = -3.0 * kh_263[k]
                   + f_0 * mh_494[k];

        t_369[k] = -3.0 * kh_264[k]
                   + f_0 * mh_495[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, kh_265, kh_266, kh_267, kh_268, \
                         kh_269, mh_496, mh_497, mh_498, mh_499, \
                         mh_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -3.0 * kh_265[k]
                   + f_0 * mh_496[k];

        t_371[k] = -3.0 * kh_266[k]
                   + f_0 * mh_497[k];

        t_372[k] = -3.0 * kh_267[k]
                   + f_0 * mh_498[k];

        t_373[k] = -3.0 * kh_268[k]
                   + f_0 * mh_499[k];

        t_374[k] = -3.0 * kh_269[k]
                   + f_0 * mh_500[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, kh_270, kh_271, kh_272, kh_273, \
                         kh_274, mh_501, mh_502, mh_503, mh_504, \
                         mh_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -3.0 * kh_270[k]
                   + f_0 * mh_501[k];

        t_376[k] = -3.0 * kh_271[k]
                   + f_0 * mh_502[k];

        t_377[k] = -3.0 * kh_272[k]
                   + f_0 * mh_503[k];

        t_378[k] = -2.0 * kh_273[k]
                   + f_0 * mh_504[k];

        t_379[k] = -2.0 * kh_274[k]
                   + f_0 * mh_505[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, kh_275, kh_276, kh_277, kh_278, \
                         kh_279, mh_506, mh_507, mh_508, mh_509, \
                         mh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * kh_275[k]
                   + f_0 * mh_506[k];

        t_381[k] = -2.0 * kh_276[k]
                   + f_0 * mh_507[k];

        t_382[k] = -2.0 * kh_277[k]
                   + f_0 * mh_508[k];

        t_383[k] = -2.0 * kh_278[k]
                   + f_0 * mh_509[k];

        t_384[k] = -2.0 * kh_279[k]
                   + f_0 * mh_510[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, kh_280, kh_281, kh_282, kh_283, \
                         kh_284, mh_511, mh_512, mh_513, mh_514, \
                         mh_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * kh_280[k]
                   + f_0 * mh_511[k];

        t_386[k] = -2.0 * kh_281[k]
                   + f_0 * mh_512[k];

        t_387[k] = -2.0 * kh_282[k]
                   + f_0 * mh_513[k];

        t_388[k] = -2.0 * kh_283[k]
                   + f_0 * mh_514[k];

        t_389[k] = -2.0 * kh_284[k]
                   + f_0 * mh_515[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, kh_285, kh_286, kh_287, kh_288, \
                         kh_289, mh_516, mh_517, mh_518, mh_519, \
                         mh_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -2.0 * kh_285[k]
                   + f_0 * mh_516[k];

        t_391[k] = -2.0 * kh_286[k]
                   + f_0 * mh_517[k];

        t_392[k] = -2.0 * kh_287[k]
                   + f_0 * mh_518[k];

        t_393[k] = -2.0 * kh_288[k]
                   + f_0 * mh_519[k];

        t_394[k] = -2.0 * kh_289[k]
                   + f_0 * mh_520[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, kh_290, kh_291, kh_292, kh_293, \
                         kh_294, mh_521, mh_522, mh_523, mh_524, \
                         mh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -2.0 * kh_290[k]
                   + f_0 * mh_521[k];

        t_396[k] = -2.0 * kh_291[k]
                   + f_0 * mh_522[k];

        t_397[k] = -2.0 * kh_292[k]
                   + f_0 * mh_523[k];

        t_398[k] = -2.0 * kh_293[k]
                   + f_0 * mh_524[k];

        t_399[k] = -kh_294[k]
                   + f_0 * mh_525[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, kh_295, kh_296, kh_297, kh_298, \
                         kh_299, mh_526, mh_527, mh_528, mh_529, \
                         mh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -kh_295[k]
                   + f_0 * mh_526[k];

        t_401[k] = -kh_296[k]
                   + f_0 * mh_527[k];

        t_402[k] = -kh_297[k]
                   + f_0 * mh_528[k];

        t_403[k] = -kh_298[k]
                   + f_0 * mh_529[k];

        t_404[k] = -kh_299[k]
                   + f_0 * mh_530[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, kh_300, kh_301, kh_302, kh_303, \
                         kh_304, mh_531, mh_532, mh_533, mh_534, \
                         mh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -kh_300[k]
                   + f_0 * mh_531[k];

        t_406[k] = -kh_301[k]
                   + f_0 * mh_532[k];

        t_407[k] = -kh_302[k]
                   + f_0 * mh_533[k];

        t_408[k] = -kh_303[k]
                   + f_0 * mh_534[k];

        t_409[k] = -kh_304[k]
                   + f_0 * mh_535[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, kh_305, kh_306, kh_307, kh_308, \
                         kh_309, mh_536, mh_537, mh_538, mh_539, \
                         mh_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -kh_305[k]
                   + f_0 * mh_536[k];

        t_411[k] = -kh_306[k]
                   + f_0 * mh_537[k];

        t_412[k] = -kh_307[k]
                   + f_0 * mh_538[k];

        t_413[k] = -kh_308[k]
                   + f_0 * mh_539[k];

        t_414[k] = -kh_309[k]
                   + f_0 * mh_540[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, kh_310, kh_311, kh_312, kh_313, \
                         kh_314, mh_541, mh_542, mh_543, mh_544, \
                         mh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -kh_310[k]
                   + f_0 * mh_541[k];

        t_416[k] = -kh_311[k]
                   + f_0 * mh_542[k];

        t_417[k] = -kh_312[k]
                   + f_0 * mh_543[k];

        t_418[k] = -kh_313[k]
                   + f_0 * mh_544[k];

        t_419[k] = -kh_314[k]
                   + f_0 * mh_545[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, t_425, t_426, t_427, mh_546, \
                         mh_547, mh_548, mh_549, mh_550, mh_551, mh_552, \
                         mh_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_0 * mh_546[k];

        t_421[k] = f_0 * mh_547[k];

        t_422[k] = f_0 * mh_548[k];

        t_423[k] = f_0 * mh_549[k];

        t_424[k] = f_0 * mh_550[k];

        t_425[k] = f_0 * mh_551[k];

        t_426[k] = f_0 * mh_552[k];

        t_427[k] = f_0 * mh_553[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, t_433, t_434, t_435, mh_554, \
                         mh_555, mh_556, mh_557, mh_558, mh_559, mh_560, \
                         mh_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_0 * mh_554[k];

        t_429[k] = f_0 * mh_555[k];

        t_430[k] = f_0 * mh_556[k];

        t_431[k] = f_0 * mh_557[k];

        t_432[k] = f_0 * mh_558[k];

        t_433[k] = f_0 * mh_559[k];

        t_434[k] = f_0 * mh_560[k];

        t_435[k] = f_0 * mh_561[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, t_441, t_442, kh_315, kh_316, \
                         mh_562, mh_563, mh_564, mh_565, mh_566, mh_588, \
                         mh_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_0 * mh_562[k];

        t_437[k] = f_0 * mh_563[k];

        t_438[k] = f_0 * mh_564[k];

        t_439[k] = f_0 * mh_565[k];

        t_440[k] = f_0 * mh_566[k];

        t_441[k] = -6.0 * kh_315[k]
                   + f_0 * mh_588[k];

        t_442[k] = -6.0 * kh_316[k]
                   + f_0 * mh_589[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, t_447, kh_317, kh_318, kh_319, kh_320, \
                         kh_321, mh_590, mh_591, mh_592, mh_593, \
                         mh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = -6.0 * kh_317[k]
                   + f_0 * mh_590[k];

        t_444[k] = -6.0 * kh_318[k]
                   + f_0 * mh_591[k];

        t_445[k] = -6.0 * kh_319[k]
                   + f_0 * mh_592[k];

        t_446[k] = -6.0 * kh_320[k]
                   + f_0 * mh_593[k];

        t_447[k] = -6.0 * kh_321[k]
                   + f_0 * mh_594[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, t_452, kh_322, kh_323, kh_324, kh_325, \
                         kh_326, mh_595, mh_596, mh_597, mh_598, \
                         mh_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = -6.0 * kh_322[k]
                   + f_0 * mh_595[k];

        t_449[k] = -6.0 * kh_323[k]
                   + f_0 * mh_596[k];

        t_450[k] = -6.0 * kh_324[k]
                   + f_0 * mh_597[k];

        t_451[k] = -6.0 * kh_325[k]
                   + f_0 * mh_598[k];

        t_452[k] = -6.0 * kh_326[k]
                   + f_0 * mh_599[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, t_457, kh_327, kh_328, kh_329, kh_330, \
                         kh_331, mh_600, mh_601, mh_602, mh_603, \
                         mh_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = -6.0 * kh_327[k]
                   + f_0 * mh_600[k];

        t_454[k] = -6.0 * kh_328[k]
                   + f_0 * mh_601[k];

        t_455[k] = -6.0 * kh_329[k]
                   + f_0 * mh_602[k];

        t_456[k] = -6.0 * kh_330[k]
                   + f_0 * mh_603[k];

        t_457[k] = -6.0 * kh_331[k]
                   + f_0 * mh_604[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, t_462, kh_332, kh_333, kh_334, kh_335, \
                         kh_336, mh_605, mh_606, mh_607, mh_608, \
                         mh_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = -6.0 * kh_332[k]
                   + f_0 * mh_605[k];

        t_459[k] = -6.0 * kh_333[k]
                   + f_0 * mh_606[k];

        t_460[k] = -6.0 * kh_334[k]
                   + f_0 * mh_607[k];

        t_461[k] = -6.0 * kh_335[k]
                   + f_0 * mh_608[k];

        t_462[k] = -5.0 * kh_336[k]
                   + f_0 * mh_609[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, t_467, kh_337, kh_338, kh_339, kh_340, \
                         kh_341, mh_610, mh_611, mh_612, mh_613, \
                         mh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = -5.0 * kh_337[k]
                   + f_0 * mh_610[k];

        t_464[k] = -5.0 * kh_338[k]
                   + f_0 * mh_611[k];

        t_465[k] = -5.0 * kh_339[k]
                   + f_0 * mh_612[k];

        t_466[k] = -5.0 * kh_340[k]
                   + f_0 * mh_613[k];

        t_467[k] = -5.0 * kh_341[k]
                   + f_0 * mh_614[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, t_471, t_472, kh_342, kh_343, kh_344, kh_345, \
                         kh_346, mh_615, mh_616, mh_617, mh_618, \
                         mh_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = -5.0 * kh_342[k]
                   + f_0 * mh_615[k];

        t_469[k] = -5.0 * kh_343[k]
                   + f_0 * mh_616[k];

        t_470[k] = -5.0 * kh_344[k]
                   + f_0 * mh_617[k];

        t_471[k] = -5.0 * kh_345[k]
                   + f_0 * mh_618[k];

        t_472[k] = -5.0 * kh_346[k]
                   + f_0 * mh_619[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, t_476, t_477, kh_347, kh_348, kh_349, kh_350, \
                         kh_351, mh_620, mh_621, mh_622, mh_623, \
                         mh_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = -5.0 * kh_347[k]
                   + f_0 * mh_620[k];

        t_474[k] = -5.0 * kh_348[k]
                   + f_0 * mh_621[k];

        t_475[k] = -5.0 * kh_349[k]
                   + f_0 * mh_622[k];

        t_476[k] = -5.0 * kh_350[k]
                   + f_0 * mh_623[k];

        t_477[k] = -5.0 * kh_351[k]
                   + f_0 * mh_624[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, t_482, kh_352, kh_353, kh_354, kh_355, \
                         kh_356, mh_625, mh_626, mh_627, mh_628, \
                         mh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = -5.0 * kh_352[k]
                   + f_0 * mh_625[k];

        t_479[k] = -5.0 * kh_353[k]
                   + f_0 * mh_626[k];

        t_480[k] = -5.0 * kh_354[k]
                   + f_0 * mh_627[k];

        t_481[k] = -5.0 * kh_355[k]
                   + f_0 * mh_628[k];

        t_482[k] = -5.0 * kh_356[k]
                   + f_0 * mh_629[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, t_487, kh_357, kh_358, kh_359, kh_360, \
                         kh_361, mh_630, mh_631, mh_632, mh_633, \
                         mh_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = -4.0 * kh_357[k]
                   + f_0 * mh_630[k];

        t_484[k] = -4.0 * kh_358[k]
                   + f_0 * mh_631[k];

        t_485[k] = -4.0 * kh_359[k]
                   + f_0 * mh_632[k];

        t_486[k] = -4.0 * kh_360[k]
                   + f_0 * mh_633[k];

        t_487[k] = -4.0 * kh_361[k]
                   + f_0 * mh_634[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);
    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);
    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);
    auto *t_500 = buffer.data(target + 500);
    auto *t_501 = buffer.data(target + 501);
    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);
    auto *t_504 = buffer.data(target + 504);
    auto *t_505 = buffer.data(target + 505);
    auto *t_506 = buffer.data(target + 506);
    auto *t_507 = buffer.data(target + 507);
    auto *t_508 = buffer.data(target + 508);
    auto *t_509 = buffer.data(target + 509);
    auto *t_510 = buffer.data(target + 510);
    auto *t_511 = buffer.data(target + 511);
    auto *t_512 = buffer.data(target + 512);
    auto *t_513 = buffer.data(target + 513);
    auto *t_514 = buffer.data(target + 514);
    auto *t_515 = buffer.data(target + 515);
    auto *t_516 = buffer.data(target + 516);
    auto *t_517 = buffer.data(target + 517);
    auto *t_518 = buffer.data(target + 518);
    auto *t_519 = buffer.data(target + 519);
    auto *t_520 = buffer.data(target + 520);
    auto *t_521 = buffer.data(target + 521);
    auto *t_522 = buffer.data(target + 522);
    auto *t_523 = buffer.data(target + 523);
    auto *t_524 = buffer.data(target + 524);
    auto *t_525 = buffer.data(target + 525);
    auto *t_526 = buffer.data(target + 526);
    auto *t_527 = buffer.data(target + 527);
    auto *t_528 = buffer.data(target + 528);
    auto *t_529 = buffer.data(target + 529);
    auto *t_530 = buffer.data(target + 530);
    auto *t_531 = buffer.data(target + 531);
    auto *t_532 = buffer.data(target + 532);
    auto *t_533 = buffer.data(target + 533);
    auto *t_534 = buffer.data(target + 534);
    auto *t_535 = buffer.data(target + 535);
    auto *t_536 = buffer.data(target + 536);
    auto *t_537 = buffer.data(target + 537);
    auto *t_538 = buffer.data(target + 538);
    auto *t_539 = buffer.data(target + 539);
    auto *t_540 = buffer.data(target + 540);
    auto *t_541 = buffer.data(target + 541);
    auto *t_542 = buffer.data(target + 542);
    auto *t_543 = buffer.data(target + 543);
    auto *t_544 = buffer.data(target + 544);
    auto *t_545 = buffer.data(target + 545);
    auto *t_546 = buffer.data(target + 546);
    auto *t_547 = buffer.data(target + 547);
    auto *t_548 = buffer.data(target + 548);
    auto *t_549 = buffer.data(target + 549);
    auto *t_550 = buffer.data(target + 550);
    auto *t_551 = buffer.data(target + 551);
    auto *t_552 = buffer.data(target + 552);
    auto *t_553 = buffer.data(target + 553);
    auto *t_554 = buffer.data(target + 554);
    auto *t_555 = buffer.data(target + 555);
    auto *t_556 = buffer.data(target + 556);
    auto *t_557 = buffer.data(target + 557);
    auto *t_558 = buffer.data(target + 558);
    auto *t_559 = buffer.data(target + 559);
    auto *t_560 = buffer.data(target + 560);
    auto *t_561 = buffer.data(target + 561);
    auto *t_562 = buffer.data(target + 562);
    auto *t_563 = buffer.data(target + 563);
    auto *t_564 = buffer.data(target + 564);
    auto *t_565 = buffer.data(target + 565);
    auto *t_566 = buffer.data(target + 566);
    auto *t_567 = buffer.data(target + 567);
    auto *t_568 = buffer.data(target + 568);
    auto *t_569 = buffer.data(target + 569);
    auto *t_570 = buffer.data(target + 570);
    auto *t_571 = buffer.data(target + 571);
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);
    auto *t_575 = buffer.data(target + 575);
    auto *t_576 = buffer.data(target + 576);
    auto *t_577 = buffer.data(target + 577);
    auto *t_578 = buffer.data(target + 578);
    auto *t_579 = buffer.data(target + 579);
    auto *t_580 = buffer.data(target + 580);
    auto *t_581 = buffer.data(target + 581);
    auto *t_582 = buffer.data(target + 582);
    auto *t_583 = buffer.data(target + 583);
    auto *t_584 = buffer.data(target + 584);
    auto *t_585 = buffer.data(target + 585);
    auto *t_586 = buffer.data(target + 586);
    auto *t_587 = buffer.data(target + 587);
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);
    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);
    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);
    auto *t_618 = buffer.data(target + 618);
    auto *t_619 = buffer.data(target + 619);
    auto *t_620 = buffer.data(target + 620);
    auto *t_621 = buffer.data(target + 621);
    auto *t_622 = buffer.data(target + 622);
    auto *t_623 = buffer.data(target + 623);
    auto *t_624 = buffer.data(target + 624);
    auto *t_625 = buffer.data(target + 625);
    auto *t_626 = buffer.data(target + 626);
    auto *t_627 = buffer.data(target + 627);
    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);
    auto *t_637 = buffer.data(target + 637);
    auto *t_638 = buffer.data(target + 638);
    auto *t_639 = buffer.data(target + 639);
    auto *t_640 = buffer.data(target + 640);
    auto *t_641 = buffer.data(target + 641);
    auto *t_642 = buffer.data(target + 642);
    auto *t_643 = buffer.data(target + 643);
    auto *t_644 = buffer.data(target + 644);
    auto *t_645 = buffer.data(target + 645);

    const auto *kh_362 = buffer.data(kh + 362);
    const auto *kh_363 = buffer.data(kh + 363);
    const auto *kh_364 = buffer.data(kh + 364);
    const auto *kh_365 = buffer.data(kh + 365);
    const auto *kh_366 = buffer.data(kh + 366);
    const auto *kh_367 = buffer.data(kh + 367);
    const auto *kh_368 = buffer.data(kh + 368);
    const auto *kh_369 = buffer.data(kh + 369);
    const auto *kh_370 = buffer.data(kh + 370);
    const auto *kh_371 = buffer.data(kh + 371);
    const auto *kh_372 = buffer.data(kh + 372);
    const auto *kh_373 = buffer.data(kh + 373);
    const auto *kh_374 = buffer.data(kh + 374);
    const auto *kh_375 = buffer.data(kh + 375);
    const auto *kh_376 = buffer.data(kh + 376);
    const auto *kh_377 = buffer.data(kh + 377);
    const auto *kh_378 = buffer.data(kh + 378);
    const auto *kh_379 = buffer.data(kh + 379);
    const auto *kh_380 = buffer.data(kh + 380);
    const auto *kh_381 = buffer.data(kh + 381);
    const auto *kh_382 = buffer.data(kh + 382);
    const auto *kh_383 = buffer.data(kh + 383);
    const auto *kh_384 = buffer.data(kh + 384);
    const auto *kh_385 = buffer.data(kh + 385);
    const auto *kh_386 = buffer.data(kh + 386);
    const auto *kh_387 = buffer.data(kh + 387);
    const auto *kh_388 = buffer.data(kh + 388);
    const auto *kh_389 = buffer.data(kh + 389);
    const auto *kh_390 = buffer.data(kh + 390);
    const auto *kh_391 = buffer.data(kh + 391);
    const auto *kh_392 = buffer.data(kh + 392);
    const auto *kh_393 = buffer.data(kh + 393);
    const auto *kh_394 = buffer.data(kh + 394);
    const auto *kh_395 = buffer.data(kh + 395);
    const auto *kh_396 = buffer.data(kh + 396);
    const auto *kh_397 = buffer.data(kh + 397);
    const auto *kh_398 = buffer.data(kh + 398);
    const auto *kh_399 = buffer.data(kh + 399);
    const auto *kh_400 = buffer.data(kh + 400);
    const auto *kh_401 = buffer.data(kh + 401);
    const auto *kh_402 = buffer.data(kh + 402);
    const auto *kh_403 = buffer.data(kh + 403);
    const auto *kh_404 = buffer.data(kh + 404);
    const auto *kh_405 = buffer.data(kh + 405);
    const auto *kh_406 = buffer.data(kh + 406);
    const auto *kh_407 = buffer.data(kh + 407);
    const auto *kh_408 = buffer.data(kh + 408);
    const auto *kh_409 = buffer.data(kh + 409);
    const auto *kh_410 = buffer.data(kh + 410);
    const auto *kh_411 = buffer.data(kh + 411);
    const auto *kh_412 = buffer.data(kh + 412);
    const auto *kh_413 = buffer.data(kh + 413);
    const auto *kh_414 = buffer.data(kh + 414);
    const auto *kh_415 = buffer.data(kh + 415);
    const auto *kh_416 = buffer.data(kh + 416);
    const auto *kh_417 = buffer.data(kh + 417);
    const auto *kh_418 = buffer.data(kh + 418);
    const auto *kh_419 = buffer.data(kh + 419);
    const auto *kh_420 = buffer.data(kh + 420);
    const auto *kh_421 = buffer.data(kh + 421);
    const auto *kh_422 = buffer.data(kh + 422);
    const auto *kh_423 = buffer.data(kh + 423);
    const auto *kh_424 = buffer.data(kh + 424);
    const auto *kh_425 = buffer.data(kh + 425);
    const auto *kh_426 = buffer.data(kh + 426);
    const auto *kh_427 = buffer.data(kh + 427);
    const auto *kh_428 = buffer.data(kh + 428);
    const auto *kh_429 = buffer.data(kh + 429);
    const auto *kh_430 = buffer.data(kh + 430);
    const auto *kh_431 = buffer.data(kh + 431);
    const auto *kh_432 = buffer.data(kh + 432);
    const auto *kh_433 = buffer.data(kh + 433);
    const auto *kh_434 = buffer.data(kh + 434);
    const auto *kh_435 = buffer.data(kh + 435);
    const auto *kh_436 = buffer.data(kh + 436);
    const auto *kh_437 = buffer.data(kh + 437);
    const auto *kh_438 = buffer.data(kh + 438);
    const auto *kh_439 = buffer.data(kh + 439);
    const auto *kh_440 = buffer.data(kh + 440);
    const auto *kh_441 = buffer.data(kh + 441);
    const auto *kh_442 = buffer.data(kh + 442);
    const auto *kh_443 = buffer.data(kh + 443);
    const auto *kh_444 = buffer.data(kh + 444);
    const auto *kh_445 = buffer.data(kh + 445);
    const auto *kh_446 = buffer.data(kh + 446);
    const auto *kh_447 = buffer.data(kh + 447);
    const auto *kh_448 = buffer.data(kh + 448);
    const auto *kh_449 = buffer.data(kh + 449);
    const auto *kh_450 = buffer.data(kh + 450);
    const auto *kh_451 = buffer.data(kh + 451);
    const auto *kh_452 = buffer.data(kh + 452);
    const auto *kh_453 = buffer.data(kh + 453);
    const auto *kh_454 = buffer.data(kh + 454);
    const auto *kh_455 = buffer.data(kh + 455);
    const auto *kh_456 = buffer.data(kh + 456);
    const auto *kh_457 = buffer.data(kh + 457);
    const auto *kh_458 = buffer.data(kh + 458);
    const auto *kh_459 = buffer.data(kh + 459);
    const auto *kh_460 = buffer.data(kh + 460);
    const auto *kh_461 = buffer.data(kh + 461);
    const auto *kh_462 = buffer.data(kh + 462);
    const auto *kh_463 = buffer.data(kh + 463);
    const auto *kh_464 = buffer.data(kh + 464);
    const auto *kh_465 = buffer.data(kh + 465);
    const auto *kh_466 = buffer.data(kh + 466);
    const auto *kh_467 = buffer.data(kh + 467);
    const auto *kh_468 = buffer.data(kh + 468);
    const auto *kh_469 = buffer.data(kh + 469);
    const auto *kh_470 = buffer.data(kh + 470);
    const auto *kh_471 = buffer.data(kh + 471);
    const auto *kh_472 = buffer.data(kh + 472);
    const auto *kh_473 = buffer.data(kh + 473);
    const auto *kh_474 = buffer.data(kh + 474);
    const auto *kh_475 = buffer.data(kh + 475);
    const auto *kh_476 = buffer.data(kh + 476);
    const auto *kh_477 = buffer.data(kh + 477);
    const auto *kh_478 = buffer.data(kh + 478);
    const auto *kh_479 = buffer.data(kh + 479);
    const auto *kh_480 = buffer.data(kh + 480);
    const auto *kh_481 = buffer.data(kh + 481);
    const auto *kh_482 = buffer.data(kh + 482);
    const auto *kh_483 = buffer.data(kh + 483);
    const auto *kh_484 = buffer.data(kh + 484);
    const auto *kh_485 = buffer.data(kh + 485);
    const auto *kh_486 = buffer.data(kh + 486);
    const auto *kh_487 = buffer.data(kh + 487);
    const auto *kh_488 = buffer.data(kh + 488);
    const auto *kh_489 = buffer.data(kh + 489);
    const auto *kh_490 = buffer.data(kh + 490);
    const auto *kh_491 = buffer.data(kh + 491);
    const auto *kh_492 = buffer.data(kh + 492);
    const auto *kh_493 = buffer.data(kh + 493);
    const auto *kh_494 = buffer.data(kh + 494);
    const auto *kh_495 = buffer.data(kh + 495);
    const auto *kh_496 = buffer.data(kh + 496);
    const auto *kh_497 = buffer.data(kh + 497);
    const auto *kh_498 = buffer.data(kh + 498);

    const auto *mh_635 = buffer.data(mh + 635);
    const auto *mh_636 = buffer.data(mh + 636);
    const auto *mh_637 = buffer.data(mh + 637);
    const auto *mh_638 = buffer.data(mh + 638);
    const auto *mh_639 = buffer.data(mh + 639);
    const auto *mh_640 = buffer.data(mh + 640);
    const auto *mh_641 = buffer.data(mh + 641);
    const auto *mh_642 = buffer.data(mh + 642);
    const auto *mh_643 = buffer.data(mh + 643);
    const auto *mh_644 = buffer.data(mh + 644);
    const auto *mh_645 = buffer.data(mh + 645);
    const auto *mh_646 = buffer.data(mh + 646);
    const auto *mh_647 = buffer.data(mh + 647);
    const auto *mh_648 = buffer.data(mh + 648);
    const auto *mh_649 = buffer.data(mh + 649);
    const auto *mh_650 = buffer.data(mh + 650);
    const auto *mh_651 = buffer.data(mh + 651);
    const auto *mh_652 = buffer.data(mh + 652);
    const auto *mh_653 = buffer.data(mh + 653);
    const auto *mh_654 = buffer.data(mh + 654);
    const auto *mh_655 = buffer.data(mh + 655);
    const auto *mh_656 = buffer.data(mh + 656);
    const auto *mh_657 = buffer.data(mh + 657);
    const auto *mh_658 = buffer.data(mh + 658);
    const auto *mh_659 = buffer.data(mh + 659);
    const auto *mh_660 = buffer.data(mh + 660);
    const auto *mh_661 = buffer.data(mh + 661);
    const auto *mh_662 = buffer.data(mh + 662);
    const auto *mh_663 = buffer.data(mh + 663);
    const auto *mh_664 = buffer.data(mh + 664);
    const auto *mh_665 = buffer.data(mh + 665);
    const auto *mh_666 = buffer.data(mh + 666);
    const auto *mh_667 = buffer.data(mh + 667);
    const auto *mh_668 = buffer.data(mh + 668);
    const auto *mh_669 = buffer.data(mh + 669);
    const auto *mh_670 = buffer.data(mh + 670);
    const auto *mh_671 = buffer.data(mh + 671);
    const auto *mh_672 = buffer.data(mh + 672);
    const auto *mh_673 = buffer.data(mh + 673);
    const auto *mh_674 = buffer.data(mh + 674);
    const auto *mh_675 = buffer.data(mh + 675);
    const auto *mh_676 = buffer.data(mh + 676);
    const auto *mh_677 = buffer.data(mh + 677);
    const auto *mh_678 = buffer.data(mh + 678);
    const auto *mh_679 = buffer.data(mh + 679);
    const auto *mh_680 = buffer.data(mh + 680);
    const auto *mh_681 = buffer.data(mh + 681);
    const auto *mh_682 = buffer.data(mh + 682);
    const auto *mh_683 = buffer.data(mh + 683);
    const auto *mh_684 = buffer.data(mh + 684);
    const auto *mh_685 = buffer.data(mh + 685);
    const auto *mh_686 = buffer.data(mh + 686);
    const auto *mh_687 = buffer.data(mh + 687);
    const auto *mh_688 = buffer.data(mh + 688);
    const auto *mh_689 = buffer.data(mh + 689);
    const auto *mh_690 = buffer.data(mh + 690);
    const auto *mh_691 = buffer.data(mh + 691);
    const auto *mh_692 = buffer.data(mh + 692);
    const auto *mh_693 = buffer.data(mh + 693);
    const auto *mh_694 = buffer.data(mh + 694);
    const auto *mh_695 = buffer.data(mh + 695);
    const auto *mh_696 = buffer.data(mh + 696);
    const auto *mh_697 = buffer.data(mh + 697);
    const auto *mh_698 = buffer.data(mh + 698);
    const auto *mh_699 = buffer.data(mh + 699);
    const auto *mh_700 = buffer.data(mh + 700);
    const auto *mh_701 = buffer.data(mh + 701);
    const auto *mh_702 = buffer.data(mh + 702);
    const auto *mh_703 = buffer.data(mh + 703);
    const auto *mh_704 = buffer.data(mh + 704);
    const auto *mh_705 = buffer.data(mh + 705);
    const auto *mh_706 = buffer.data(mh + 706);
    const auto *mh_707 = buffer.data(mh + 707);
    const auto *mh_708 = buffer.data(mh + 708);
    const auto *mh_709 = buffer.data(mh + 709);
    const auto *mh_710 = buffer.data(mh + 710);
    const auto *mh_711 = buffer.data(mh + 711);
    const auto *mh_712 = buffer.data(mh + 712);
    const auto *mh_713 = buffer.data(mh + 713);
    const auto *mh_714 = buffer.data(mh + 714);
    const auto *mh_715 = buffer.data(mh + 715);
    const auto *mh_716 = buffer.data(mh + 716);
    const auto *mh_717 = buffer.data(mh + 717);
    const auto *mh_718 = buffer.data(mh + 718);
    const auto *mh_719 = buffer.data(mh + 719);
    const auto *mh_720 = buffer.data(mh + 720);
    const auto *mh_721 = buffer.data(mh + 721);
    const auto *mh_722 = buffer.data(mh + 722);
    const auto *mh_723 = buffer.data(mh + 723);
    const auto *mh_724 = buffer.data(mh + 724);
    const auto *mh_725 = buffer.data(mh + 725);
    const auto *mh_726 = buffer.data(mh + 726);
    const auto *mh_727 = buffer.data(mh + 727);
    const auto *mh_728 = buffer.data(mh + 728);
    const auto *mh_729 = buffer.data(mh + 729);
    const auto *mh_730 = buffer.data(mh + 730);
    const auto *mh_731 = buffer.data(mh + 731);
    const auto *mh_732 = buffer.data(mh + 732);
    const auto *mh_733 = buffer.data(mh + 733);
    const auto *mh_734 = buffer.data(mh + 734);
    const auto *mh_756 = buffer.data(mh + 756);
    const auto *mh_757 = buffer.data(mh + 757);
    const auto *mh_758 = buffer.data(mh + 758);
    const auto *mh_759 = buffer.data(mh + 759);
    const auto *mh_760 = buffer.data(mh + 760);
    const auto *mh_761 = buffer.data(mh + 761);
    const auto *mh_762 = buffer.data(mh + 762);
    const auto *mh_763 = buffer.data(mh + 763);
    const auto *mh_764 = buffer.data(mh + 764);
    const auto *mh_765 = buffer.data(mh + 765);
    const auto *mh_766 = buffer.data(mh + 766);
    const auto *mh_767 = buffer.data(mh + 767);
    const auto *mh_768 = buffer.data(mh + 768);
    const auto *mh_769 = buffer.data(mh + 769);
    const auto *mh_770 = buffer.data(mh + 770);
    const auto *mh_771 = buffer.data(mh + 771);
    const auto *mh_772 = buffer.data(mh + 772);
    const auto *mh_773 = buffer.data(mh + 773);
    const auto *mh_774 = buffer.data(mh + 774);
    const auto *mh_775 = buffer.data(mh + 775);
    const auto *mh_776 = buffer.data(mh + 776);
    const auto *mh_777 = buffer.data(mh + 777);
    const auto *mh_778 = buffer.data(mh + 778);
    const auto *mh_779 = buffer.data(mh + 779);
    const auto *mh_780 = buffer.data(mh + 780);
    const auto *mh_781 = buffer.data(mh + 781);
    const auto *mh_782 = buffer.data(mh + 782);
    const auto *mh_783 = buffer.data(mh + 783);
    const auto *mh_784 = buffer.data(mh + 784);
    const auto *mh_785 = buffer.data(mh + 785);
    const auto *mh_786 = buffer.data(mh + 786);
    const auto *mh_787 = buffer.data(mh + 787);
    const auto *mh_788 = buffer.data(mh + 788);
    const auto *mh_789 = buffer.data(mh + 789);
    const auto *mh_790 = buffer.data(mh + 790);
    const auto *mh_791 = buffer.data(mh + 791);
    const auto *mh_792 = buffer.data(mh + 792);
    const auto *mh_793 = buffer.data(mh + 793);
    const auto *mh_794 = buffer.data(mh + 794);
    const auto *mh_795 = buffer.data(mh + 795);
    const auto *mh_796 = buffer.data(mh + 796);
    const auto *mh_797 = buffer.data(mh + 797);
    const auto *mh_798 = buffer.data(mh + 798);
    const auto *mh_799 = buffer.data(mh + 799);
    const auto *mh_800 = buffer.data(mh + 800);
    const auto *mh_801 = buffer.data(mh + 801);
    const auto *mh_802 = buffer.data(mh + 802);
    const auto *mh_803 = buffer.data(mh + 803);
    const auto *mh_804 = buffer.data(mh + 804);
    const auto *mh_805 = buffer.data(mh + 805);
    const auto *mh_806 = buffer.data(mh + 806);
    const auto *mh_807 = buffer.data(mh + 807);
    const auto *mh_808 = buffer.data(mh + 808);
    const auto *mh_809 = buffer.data(mh + 809);
    const auto *mh_810 = buffer.data(mh + 810);
    const auto *mh_811 = buffer.data(mh + 811);
    const auto *mh_812 = buffer.data(mh + 812);
    const auto *mh_813 = buffer.data(mh + 813);

#pragma omp simd aligned(t_488, t_489, t_490, t_491, t_492, kh_362, kh_363, kh_364, kh_365, \
                         kh_366, mh_635, mh_636, mh_637, mh_638, \
                         mh_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = -4.0 * kh_362[k]
                   + f_0 * mh_635[k];

        t_489[k] = -4.0 * kh_363[k]
                   + f_0 * mh_636[k];

        t_490[k] = -4.0 * kh_364[k]
                   + f_0 * mh_637[k];

        t_491[k] = -4.0 * kh_365[k]
                   + f_0 * mh_638[k];

        t_492[k] = -4.0 * kh_366[k]
                   + f_0 * mh_639[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, t_496, t_497, kh_367, kh_368, kh_369, kh_370, \
                         kh_371, mh_640, mh_641, mh_642, mh_643, \
                         mh_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = -4.0 * kh_367[k]
                   + f_0 * mh_640[k];

        t_494[k] = -4.0 * kh_368[k]
                   + f_0 * mh_641[k];

        t_495[k] = -4.0 * kh_369[k]
                   + f_0 * mh_642[k];

        t_496[k] = -4.0 * kh_370[k]
                   + f_0 * mh_643[k];

        t_497[k] = -4.0 * kh_371[k]
                   + f_0 * mh_644[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, t_501, t_502, kh_372, kh_373, kh_374, kh_375, \
                         kh_376, mh_645, mh_646, mh_647, mh_648, \
                         mh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = -4.0 * kh_372[k]
                   + f_0 * mh_645[k];

        t_499[k] = -4.0 * kh_373[k]
                   + f_0 * mh_646[k];

        t_500[k] = -4.0 * kh_374[k]
                   + f_0 * mh_647[k];

        t_501[k] = -4.0 * kh_375[k]
                   + f_0 * mh_648[k];

        t_502[k] = -4.0 * kh_376[k]
                   + f_0 * mh_649[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, t_507, kh_377, kh_378, kh_379, kh_380, \
                         kh_381, mh_650, mh_651, mh_652, mh_653, \
                         mh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = -4.0 * kh_377[k]
                   + f_0 * mh_650[k];

        t_504[k] = -3.0 * kh_378[k]
                   + f_0 * mh_651[k];

        t_505[k] = -3.0 * kh_379[k]
                   + f_0 * mh_652[k];

        t_506[k] = -3.0 * kh_380[k]
                   + f_0 * mh_653[k];

        t_507[k] = -3.0 * kh_381[k]
                   + f_0 * mh_654[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, t_512, kh_382, kh_383, kh_384, kh_385, \
                         kh_386, mh_655, mh_656, mh_657, mh_658, \
                         mh_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = -3.0 * kh_382[k]
                   + f_0 * mh_655[k];

        t_509[k] = -3.0 * kh_383[k]
                   + f_0 * mh_656[k];

        t_510[k] = -3.0 * kh_384[k]
                   + f_0 * mh_657[k];

        t_511[k] = -3.0 * kh_385[k]
                   + f_0 * mh_658[k];

        t_512[k] = -3.0 * kh_386[k]
                   + f_0 * mh_659[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, t_516, t_517, kh_387, kh_388, kh_389, kh_390, \
                         kh_391, mh_660, mh_661, mh_662, mh_663, \
                         mh_664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = -3.0 * kh_387[k]
                   + f_0 * mh_660[k];

        t_514[k] = -3.0 * kh_388[k]
                   + f_0 * mh_661[k];

        t_515[k] = -3.0 * kh_389[k]
                   + f_0 * mh_662[k];

        t_516[k] = -3.0 * kh_390[k]
                   + f_0 * mh_663[k];

        t_517[k] = -3.0 * kh_391[k]
                   + f_0 * mh_664[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, t_522, kh_392, kh_393, kh_394, kh_395, \
                         kh_396, mh_665, mh_666, mh_667, mh_668, \
                         mh_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = -3.0 * kh_392[k]
                   + f_0 * mh_665[k];

        t_519[k] = -3.0 * kh_393[k]
                   + f_0 * mh_666[k];

        t_520[k] = -3.0 * kh_394[k]
                   + f_0 * mh_667[k];

        t_521[k] = -3.0 * kh_395[k]
                   + f_0 * mh_668[k];

        t_522[k] = -3.0 * kh_396[k]
                   + f_0 * mh_669[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, t_527, kh_397, kh_398, kh_399, kh_400, \
                         kh_401, mh_670, mh_671, mh_672, mh_673, \
                         mh_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = -3.0 * kh_397[k]
                   + f_0 * mh_670[k];

        t_524[k] = -3.0 * kh_398[k]
                   + f_0 * mh_671[k];

        t_525[k] = -2.0 * kh_399[k]
                   + f_0 * mh_672[k];

        t_526[k] = -2.0 * kh_400[k]
                   + f_0 * mh_673[k];

        t_527[k] = -2.0 * kh_401[k]
                   + f_0 * mh_674[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, kh_402, kh_403, kh_404, kh_405, \
                         kh_406, mh_675, mh_676, mh_677, mh_678, \
                         mh_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = -2.0 * kh_402[k]
                   + f_0 * mh_675[k];

        t_529[k] = -2.0 * kh_403[k]
                   + f_0 * mh_676[k];

        t_530[k] = -2.0 * kh_404[k]
                   + f_0 * mh_677[k];

        t_531[k] = -2.0 * kh_405[k]
                   + f_0 * mh_678[k];

        t_532[k] = -2.0 * kh_406[k]
                   + f_0 * mh_679[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, t_537, kh_407, kh_408, kh_409, kh_410, \
                         kh_411, mh_680, mh_681, mh_682, mh_683, \
                         mh_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = -2.0 * kh_407[k]
                   + f_0 * mh_680[k];

        t_534[k] = -2.0 * kh_408[k]
                   + f_0 * mh_681[k];

        t_535[k] = -2.0 * kh_409[k]
                   + f_0 * mh_682[k];

        t_536[k] = -2.0 * kh_410[k]
                   + f_0 * mh_683[k];

        t_537[k] = -2.0 * kh_411[k]
                   + f_0 * mh_684[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, t_542, kh_412, kh_413, kh_414, kh_415, \
                         kh_416, mh_685, mh_686, mh_687, mh_688, \
                         mh_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = -2.0 * kh_412[k]
                   + f_0 * mh_685[k];

        t_539[k] = -2.0 * kh_413[k]
                   + f_0 * mh_686[k];

        t_540[k] = -2.0 * kh_414[k]
                   + f_0 * mh_687[k];

        t_541[k] = -2.0 * kh_415[k]
                   + f_0 * mh_688[k];

        t_542[k] = -2.0 * kh_416[k]
                   + f_0 * mh_689[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, t_547, kh_417, kh_418, kh_419, kh_420, \
                         kh_421, mh_690, mh_691, mh_692, mh_693, \
                         mh_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = -2.0 * kh_417[k]
                   + f_0 * mh_690[k];

        t_544[k] = -2.0 * kh_418[k]
                   + f_0 * mh_691[k];

        t_545[k] = -2.0 * kh_419[k]
                   + f_0 * mh_692[k];

        t_546[k] = -kh_420[k]
                   + f_0 * mh_693[k];

        t_547[k] = -kh_421[k]
                   + f_0 * mh_694[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, kh_422, kh_423, kh_424, kh_425, \
                         kh_426, mh_695, mh_696, mh_697, mh_698, \
                         mh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = -kh_422[k]
                   + f_0 * mh_695[k];

        t_549[k] = -kh_423[k]
                   + f_0 * mh_696[k];

        t_550[k] = -kh_424[k]
                   + f_0 * mh_697[k];

        t_551[k] = -kh_425[k]
                   + f_0 * mh_698[k];

        t_552[k] = -kh_426[k]
                   + f_0 * mh_699[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, t_557, kh_427, kh_428, kh_429, kh_430, \
                         kh_431, mh_700, mh_701, mh_702, mh_703, \
                         mh_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = -kh_427[k]
                   + f_0 * mh_700[k];

        t_554[k] = -kh_428[k]
                   + f_0 * mh_701[k];

        t_555[k] = -kh_429[k]
                   + f_0 * mh_702[k];

        t_556[k] = -kh_430[k]
                   + f_0 * mh_703[k];

        t_557[k] = -kh_431[k]
                   + f_0 * mh_704[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, t_562, kh_432, kh_433, kh_434, kh_435, \
                         kh_436, mh_705, mh_706, mh_707, mh_708, \
                         mh_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = -kh_432[k]
                   + f_0 * mh_705[k];

        t_559[k] = -kh_433[k]
                   + f_0 * mh_706[k];

        t_560[k] = -kh_434[k]
                   + f_0 * mh_707[k];

        t_561[k] = -kh_435[k]
                   + f_0 * mh_708[k];

        t_562[k] = -kh_436[k]
                   + f_0 * mh_709[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, t_567, t_568, kh_437, kh_438, kh_439, \
                         kh_440, mh_710, mh_711, mh_712, mh_713, mh_714, \
                         mh_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = -kh_437[k]
                   + f_0 * mh_710[k];

        t_564[k] = -kh_438[k]
                   + f_0 * mh_711[k];

        t_565[k] = -kh_439[k]
                   + f_0 * mh_712[k];

        t_566[k] = -kh_440[k]
                   + f_0 * mh_713[k];

        t_567[k] = f_0 * mh_714[k];

        t_568[k] = f_0 * mh_715[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, t_573, t_574, t_575, t_576, mh_716, \
                         mh_717, mh_718, mh_719, mh_720, mh_721, mh_722, \
                         mh_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_0 * mh_716[k];

        t_570[k] = f_0 * mh_717[k];

        t_571[k] = f_0 * mh_718[k];

        t_572[k] = f_0 * mh_719[k];

        t_573[k] = f_0 * mh_720[k];

        t_574[k] = f_0 * mh_721[k];

        t_575[k] = f_0 * mh_722[k];

        t_576[k] = f_0 * mh_723[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, t_582, t_583, t_584, mh_724, \
                         mh_725, mh_726, mh_727, mh_728, mh_729, mh_730, \
                         mh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_0 * mh_724[k];

        t_578[k] = f_0 * mh_725[k];

        t_579[k] = f_0 * mh_726[k];

        t_580[k] = f_0 * mh_727[k];

        t_581[k] = f_0 * mh_728[k];

        t_582[k] = f_0 * mh_729[k];

        t_583[k] = f_0 * mh_730[k];

        t_584[k] = f_0 * mh_731[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, t_590, kh_441, kh_442, kh_443, \
                         mh_732, mh_733, mh_734, mh_756, mh_757, \
                         mh_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_0 * mh_732[k];

        t_586[k] = f_0 * mh_733[k];

        t_587[k] = f_0 * mh_734[k];

        t_588[k] = -7.0 * kh_441[k]
                   + f_0 * mh_756[k];

        t_589[k] = -7.0 * kh_442[k]
                   + f_0 * mh_757[k];

        t_590[k] = -7.0 * kh_443[k]
                   + f_0 * mh_758[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, t_595, kh_444, kh_445, kh_446, kh_447, \
                         kh_448, mh_759, mh_760, mh_761, mh_762, \
                         mh_763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = -7.0 * kh_444[k]
                   + f_0 * mh_759[k];

        t_592[k] = -7.0 * kh_445[k]
                   + f_0 * mh_760[k];

        t_593[k] = -7.0 * kh_446[k]
                   + f_0 * mh_761[k];

        t_594[k] = -7.0 * kh_447[k]
                   + f_0 * mh_762[k];

        t_595[k] = -7.0 * kh_448[k]
                   + f_0 * mh_763[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, t_600, kh_449, kh_450, kh_451, kh_452, \
                         kh_453, mh_764, mh_765, mh_766, mh_767, \
                         mh_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = -7.0 * kh_449[k]
                   + f_0 * mh_764[k];

        t_597[k] = -7.0 * kh_450[k]
                   + f_0 * mh_765[k];

        t_598[k] = -7.0 * kh_451[k]
                   + f_0 * mh_766[k];

        t_599[k] = -7.0 * kh_452[k]
                   + f_0 * mh_767[k];

        t_600[k] = -7.0 * kh_453[k]
                   + f_0 * mh_768[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, t_604, t_605, kh_454, kh_455, kh_456, kh_457, \
                         kh_458, mh_769, mh_770, mh_771, mh_772, \
                         mh_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = -7.0 * kh_454[k]
                   + f_0 * mh_769[k];

        t_602[k] = -7.0 * kh_455[k]
                   + f_0 * mh_770[k];

        t_603[k] = -7.0 * kh_456[k]
                   + f_0 * mh_771[k];

        t_604[k] = -7.0 * kh_457[k]
                   + f_0 * mh_772[k];

        t_605[k] = -7.0 * kh_458[k]
                   + f_0 * mh_773[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, t_609, t_610, kh_459, kh_460, kh_461, kh_462, \
                         kh_463, mh_774, mh_775, mh_776, mh_777, \
                         mh_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = -7.0 * kh_459[k]
                   + f_0 * mh_774[k];

        t_607[k] = -7.0 * kh_460[k]
                   + f_0 * mh_775[k];

        t_608[k] = -7.0 * kh_461[k]
                   + f_0 * mh_776[k];

        t_609[k] = -6.0 * kh_462[k]
                   + f_0 * mh_777[k];

        t_610[k] = -6.0 * kh_463[k]
                   + f_0 * mh_778[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, t_614, t_615, kh_464, kh_465, kh_466, kh_467, \
                         kh_468, mh_779, mh_780, mh_781, mh_782, \
                         mh_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = -6.0 * kh_464[k]
                   + f_0 * mh_779[k];

        t_612[k] = -6.0 * kh_465[k]
                   + f_0 * mh_780[k];

        t_613[k] = -6.0 * kh_466[k]
                   + f_0 * mh_781[k];

        t_614[k] = -6.0 * kh_467[k]
                   + f_0 * mh_782[k];

        t_615[k] = -6.0 * kh_468[k]
                   + f_0 * mh_783[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, t_620, kh_469, kh_470, kh_471, kh_472, \
                         kh_473, mh_784, mh_785, mh_786, mh_787, \
                         mh_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = -6.0 * kh_469[k]
                   + f_0 * mh_784[k];

        t_617[k] = -6.0 * kh_470[k]
                   + f_0 * mh_785[k];

        t_618[k] = -6.0 * kh_471[k]
                   + f_0 * mh_786[k];

        t_619[k] = -6.0 * kh_472[k]
                   + f_0 * mh_787[k];

        t_620[k] = -6.0 * kh_473[k]
                   + f_0 * mh_788[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, t_624, t_625, kh_474, kh_475, kh_476, kh_477, \
                         kh_478, mh_789, mh_790, mh_791, mh_792, \
                         mh_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = -6.0 * kh_474[k]
                   + f_0 * mh_789[k];

        t_622[k] = -6.0 * kh_475[k]
                   + f_0 * mh_790[k];

        t_623[k] = -6.0 * kh_476[k]
                   + f_0 * mh_791[k];

        t_624[k] = -6.0 * kh_477[k]
                   + f_0 * mh_792[k];

        t_625[k] = -6.0 * kh_478[k]
                   + f_0 * mh_793[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, t_630, kh_479, kh_480, kh_481, kh_482, \
                         kh_483, mh_794, mh_795, mh_796, mh_797, \
                         mh_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = -6.0 * kh_479[k]
                   + f_0 * mh_794[k];

        t_627[k] = -6.0 * kh_480[k]
                   + f_0 * mh_795[k];

        t_628[k] = -6.0 * kh_481[k]
                   + f_0 * mh_796[k];

        t_629[k] = -6.0 * kh_482[k]
                   + f_0 * mh_797[k];

        t_630[k] = -5.0 * kh_483[k]
                   + f_0 * mh_798[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, t_635, kh_484, kh_485, kh_486, kh_487, \
                         kh_488, mh_799, mh_800, mh_801, mh_802, \
                         mh_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = -5.0 * kh_484[k]
                   + f_0 * mh_799[k];

        t_632[k] = -5.0 * kh_485[k]
                   + f_0 * mh_800[k];

        t_633[k] = -5.0 * kh_486[k]
                   + f_0 * mh_801[k];

        t_634[k] = -5.0 * kh_487[k]
                   + f_0 * mh_802[k];

        t_635[k] = -5.0 * kh_488[k]
                   + f_0 * mh_803[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, t_640, kh_489, kh_490, kh_491, kh_492, \
                         kh_493, mh_804, mh_805, mh_806, mh_807, \
                         mh_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = -5.0 * kh_489[k]
                   + f_0 * mh_804[k];

        t_637[k] = -5.0 * kh_490[k]
                   + f_0 * mh_805[k];

        t_638[k] = -5.0 * kh_491[k]
                   + f_0 * mh_806[k];

        t_639[k] = -5.0 * kh_492[k]
                   + f_0 * mh_807[k];

        t_640[k] = -5.0 * kh_493[k]
                   + f_0 * mh_808[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, t_645, kh_494, kh_495, kh_496, kh_497, \
                         kh_498, mh_809, mh_810, mh_811, mh_812, \
                         mh_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = -5.0 * kh_494[k]
                   + f_0 * mh_809[k];

        t_642[k] = -5.0 * kh_495[k]
                   + f_0 * mh_810[k];

        t_643[k] = -5.0 * kh_496[k]
                   + f_0 * mh_811[k];

        t_644[k] = -5.0 * kh_497[k]
                   + f_0 * mh_812[k];

        t_645[k] = -5.0 * kh_498[k]
                   + f_0 * mh_813[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);
    auto *t_654 = buffer.data(target + 654);
    auto *t_655 = buffer.data(target + 655);
    auto *t_656 = buffer.data(target + 656);
    auto *t_657 = buffer.data(target + 657);
    auto *t_658 = buffer.data(target + 658);
    auto *t_659 = buffer.data(target + 659);
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);
    auto *t_663 = buffer.data(target + 663);
    auto *t_664 = buffer.data(target + 664);
    auto *t_665 = buffer.data(target + 665);
    auto *t_666 = buffer.data(target + 666);
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);
    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);
    auto *t_686 = buffer.data(target + 686);
    auto *t_687 = buffer.data(target + 687);
    auto *t_688 = buffer.data(target + 688);
    auto *t_689 = buffer.data(target + 689);
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);
    auto *t_693 = buffer.data(target + 693);
    auto *t_694 = buffer.data(target + 694);
    auto *t_695 = buffer.data(target + 695);
    auto *t_696 = buffer.data(target + 696);
    auto *t_697 = buffer.data(target + 697);
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);
    auto *t_708 = buffer.data(target + 708);
    auto *t_709 = buffer.data(target + 709);
    auto *t_710 = buffer.data(target + 710);
    auto *t_711 = buffer.data(target + 711);
    auto *t_712 = buffer.data(target + 712);
    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
    auto *t_722 = buffer.data(target + 722);
    auto *t_723 = buffer.data(target + 723);
    auto *t_724 = buffer.data(target + 724);
    auto *t_725 = buffer.data(target + 725);
    auto *t_726 = buffer.data(target + 726);
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);
    auto *t_732 = buffer.data(target + 732);
    auto *t_733 = buffer.data(target + 733);
    auto *t_734 = buffer.data(target + 734);
    auto *t_735 = buffer.data(target + 735);
    auto *t_736 = buffer.data(target + 736);
    auto *t_737 = buffer.data(target + 737);
    auto *t_738 = buffer.data(target + 738);
    auto *t_739 = buffer.data(target + 739);
    auto *t_740 = buffer.data(target + 740);
    auto *t_741 = buffer.data(target + 741);
    auto *t_742 = buffer.data(target + 742);
    auto *t_743 = buffer.data(target + 743);
    auto *t_744 = buffer.data(target + 744);
    auto *t_745 = buffer.data(target + 745);
    auto *t_746 = buffer.data(target + 746);
    auto *t_747 = buffer.data(target + 747);
    auto *t_748 = buffer.data(target + 748);
    auto *t_749 = buffer.data(target + 749);
    auto *t_750 = buffer.data(target + 750);
    auto *t_751 = buffer.data(target + 751);
    auto *t_752 = buffer.data(target + 752);
    auto *t_753 = buffer.data(target + 753);
    auto *t_754 = buffer.data(target + 754);
    auto *t_755 = buffer.data(target + 755);
    auto *t_756 = buffer.data(target + 756);
    auto *t_757 = buffer.data(target + 757);
    auto *t_758 = buffer.data(target + 758);
    auto *t_759 = buffer.data(target + 759);
    auto *t_760 = buffer.data(target + 760);
    auto *t_761 = buffer.data(target + 761);
    auto *t_762 = buffer.data(target + 762);
    auto *t_763 = buffer.data(target + 763);
    auto *t_764 = buffer.data(target + 764);
    auto *t_765 = buffer.data(target + 765);
    auto *t_766 = buffer.data(target + 766);
    auto *t_767 = buffer.data(target + 767);
    auto *t_768 = buffer.data(target + 768);
    auto *t_769 = buffer.data(target + 769);
    auto *t_770 = buffer.data(target + 770);
    auto *t_771 = buffer.data(target + 771);
    auto *t_772 = buffer.data(target + 772);
    auto *t_773 = buffer.data(target + 773);
    auto *t_774 = buffer.data(target + 774);
    auto *t_775 = buffer.data(target + 775);
    auto *t_776 = buffer.data(target + 776);
    auto *t_777 = buffer.data(target + 777);
    auto *t_778 = buffer.data(target + 778);
    auto *t_779 = buffer.data(target + 779);
    auto *t_780 = buffer.data(target + 780);
    auto *t_781 = buffer.data(target + 781);
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);
    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);
    auto *t_801 = buffer.data(target + 801);
    auto *t_802 = buffer.data(target + 802);
    auto *t_803 = buffer.data(target + 803);

    const auto *kh_499 = buffer.data(kh + 499);
    const auto *kh_500 = buffer.data(kh + 500);
    const auto *kh_501 = buffer.data(kh + 501);
    const auto *kh_502 = buffer.data(kh + 502);
    const auto *kh_503 = buffer.data(kh + 503);
    const auto *kh_504 = buffer.data(kh + 504);
    const auto *kh_505 = buffer.data(kh + 505);
    const auto *kh_506 = buffer.data(kh + 506);
    const auto *kh_507 = buffer.data(kh + 507);
    const auto *kh_508 = buffer.data(kh + 508);
    const auto *kh_509 = buffer.data(kh + 509);
    const auto *kh_510 = buffer.data(kh + 510);
    const auto *kh_511 = buffer.data(kh + 511);
    const auto *kh_512 = buffer.data(kh + 512);
    const auto *kh_513 = buffer.data(kh + 513);
    const auto *kh_514 = buffer.data(kh + 514);
    const auto *kh_515 = buffer.data(kh + 515);
    const auto *kh_516 = buffer.data(kh + 516);
    const auto *kh_517 = buffer.data(kh + 517);
    const auto *kh_518 = buffer.data(kh + 518);
    const auto *kh_519 = buffer.data(kh + 519);
    const auto *kh_520 = buffer.data(kh + 520);
    const auto *kh_521 = buffer.data(kh + 521);
    const auto *kh_522 = buffer.data(kh + 522);
    const auto *kh_523 = buffer.data(kh + 523);
    const auto *kh_524 = buffer.data(kh + 524);
    const auto *kh_525 = buffer.data(kh + 525);
    const auto *kh_526 = buffer.data(kh + 526);
    const auto *kh_527 = buffer.data(kh + 527);
    const auto *kh_528 = buffer.data(kh + 528);
    const auto *kh_529 = buffer.data(kh + 529);
    const auto *kh_530 = buffer.data(kh + 530);
    const auto *kh_531 = buffer.data(kh + 531);
    const auto *kh_532 = buffer.data(kh + 532);
    const auto *kh_533 = buffer.data(kh + 533);
    const auto *kh_534 = buffer.data(kh + 534);
    const auto *kh_535 = buffer.data(kh + 535);
    const auto *kh_536 = buffer.data(kh + 536);
    const auto *kh_537 = buffer.data(kh + 537);
    const auto *kh_538 = buffer.data(kh + 538);
    const auto *kh_539 = buffer.data(kh + 539);
    const auto *kh_540 = buffer.data(kh + 540);
    const auto *kh_541 = buffer.data(kh + 541);
    const auto *kh_542 = buffer.data(kh + 542);
    const auto *kh_543 = buffer.data(kh + 543);
    const auto *kh_544 = buffer.data(kh + 544);
    const auto *kh_545 = buffer.data(kh + 545);
    const auto *kh_546 = buffer.data(kh + 546);
    const auto *kh_547 = buffer.data(kh + 547);
    const auto *kh_548 = buffer.data(kh + 548);
    const auto *kh_549 = buffer.data(kh + 549);
    const auto *kh_550 = buffer.data(kh + 550);
    const auto *kh_551 = buffer.data(kh + 551);
    const auto *kh_552 = buffer.data(kh + 552);
    const auto *kh_553 = buffer.data(kh + 553);
    const auto *kh_554 = buffer.data(kh + 554);
    const auto *kh_555 = buffer.data(kh + 555);
    const auto *kh_556 = buffer.data(kh + 556);
    const auto *kh_557 = buffer.data(kh + 557);
    const auto *kh_558 = buffer.data(kh + 558);
    const auto *kh_559 = buffer.data(kh + 559);
    const auto *kh_560 = buffer.data(kh + 560);
    const auto *kh_561 = buffer.data(kh + 561);
    const auto *kh_562 = buffer.data(kh + 562);
    const auto *kh_563 = buffer.data(kh + 563);
    const auto *kh_564 = buffer.data(kh + 564);
    const auto *kh_565 = buffer.data(kh + 565);
    const auto *kh_566 = buffer.data(kh + 566);
    const auto *kh_567 = buffer.data(kh + 567);
    const auto *kh_568 = buffer.data(kh + 568);
    const auto *kh_569 = buffer.data(kh + 569);
    const auto *kh_570 = buffer.data(kh + 570);
    const auto *kh_571 = buffer.data(kh + 571);
    const auto *kh_572 = buffer.data(kh + 572);
    const auto *kh_573 = buffer.data(kh + 573);
    const auto *kh_574 = buffer.data(kh + 574);
    const auto *kh_575 = buffer.data(kh + 575);
    const auto *kh_576 = buffer.data(kh + 576);
    const auto *kh_577 = buffer.data(kh + 577);
    const auto *kh_578 = buffer.data(kh + 578);
    const auto *kh_579 = buffer.data(kh + 579);
    const auto *kh_580 = buffer.data(kh + 580);
    const auto *kh_581 = buffer.data(kh + 581);
    const auto *kh_582 = buffer.data(kh + 582);
    const auto *kh_583 = buffer.data(kh + 583);
    const auto *kh_584 = buffer.data(kh + 584);
    const auto *kh_585 = buffer.data(kh + 585);
    const auto *kh_586 = buffer.data(kh + 586);
    const auto *kh_587 = buffer.data(kh + 587);
    const auto *kh_588 = buffer.data(kh + 588);
    const auto *kh_589 = buffer.data(kh + 589);
    const auto *kh_590 = buffer.data(kh + 590);
    const auto *kh_591 = buffer.data(kh + 591);
    const auto *kh_592 = buffer.data(kh + 592);
    const auto *kh_593 = buffer.data(kh + 593);
    const auto *kh_594 = buffer.data(kh + 594);
    const auto *kh_595 = buffer.data(kh + 595);
    const auto *kh_596 = buffer.data(kh + 596);
    const auto *kh_597 = buffer.data(kh + 597);
    const auto *kh_598 = buffer.data(kh + 598);
    const auto *kh_599 = buffer.data(kh + 599);
    const auto *kh_600 = buffer.data(kh + 600);
    const auto *kh_601 = buffer.data(kh + 601);
    const auto *kh_602 = buffer.data(kh + 602);
    const auto *kh_603 = buffer.data(kh + 603);
    const auto *kh_604 = buffer.data(kh + 604);
    const auto *kh_605 = buffer.data(kh + 605);
    const auto *kh_606 = buffer.data(kh + 606);
    const auto *kh_607 = buffer.data(kh + 607);
    const auto *kh_608 = buffer.data(kh + 608);
    const auto *kh_609 = buffer.data(kh + 609);
    const auto *kh_610 = buffer.data(kh + 610);
    const auto *kh_611 = buffer.data(kh + 611);
    const auto *kh_612 = buffer.data(kh + 612);
    const auto *kh_613 = buffer.data(kh + 613);
    const auto *kh_614 = buffer.data(kh + 614);
    const auto *kh_615 = buffer.data(kh + 615);
    const auto *kh_616 = buffer.data(kh + 616);
    const auto *kh_617 = buffer.data(kh + 617);
    const auto *kh_618 = buffer.data(kh + 618);
    const auto *kh_619 = buffer.data(kh + 619);
    const auto *kh_620 = buffer.data(kh + 620);
    const auto *kh_621 = buffer.data(kh + 621);
    const auto *kh_622 = buffer.data(kh + 622);
    const auto *kh_623 = buffer.data(kh + 623);
    const auto *kh_624 = buffer.data(kh + 624);
    const auto *kh_625 = buffer.data(kh + 625);
    const auto *kh_626 = buffer.data(kh + 626);
    const auto *kh_627 = buffer.data(kh + 627);
    const auto *kh_628 = buffer.data(kh + 628);
    const auto *kh_629 = buffer.data(kh + 629);
    const auto *kh_630 = buffer.data(kh + 630);
    const auto *kh_631 = buffer.data(kh + 631);
    const auto *kh_632 = buffer.data(kh + 632);
    const auto *kh_633 = buffer.data(kh + 633);
    const auto *kh_634 = buffer.data(kh + 634);
    const auto *kh_635 = buffer.data(kh + 635);

    const auto *mh_814 = buffer.data(mh + 814);
    const auto *mh_815 = buffer.data(mh + 815);
    const auto *mh_816 = buffer.data(mh + 816);
    const auto *mh_817 = buffer.data(mh + 817);
    const auto *mh_818 = buffer.data(mh + 818);
    const auto *mh_819 = buffer.data(mh + 819);
    const auto *mh_820 = buffer.data(mh + 820);
    const auto *mh_821 = buffer.data(mh + 821);
    const auto *mh_822 = buffer.data(mh + 822);
    const auto *mh_823 = buffer.data(mh + 823);
    const auto *mh_824 = buffer.data(mh + 824);
    const auto *mh_825 = buffer.data(mh + 825);
    const auto *mh_826 = buffer.data(mh + 826);
    const auto *mh_827 = buffer.data(mh + 827);
    const auto *mh_828 = buffer.data(mh + 828);
    const auto *mh_829 = buffer.data(mh + 829);
    const auto *mh_830 = buffer.data(mh + 830);
    const auto *mh_831 = buffer.data(mh + 831);
    const auto *mh_832 = buffer.data(mh + 832);
    const auto *mh_833 = buffer.data(mh + 833);
    const auto *mh_834 = buffer.data(mh + 834);
    const auto *mh_835 = buffer.data(mh + 835);
    const auto *mh_836 = buffer.data(mh + 836);
    const auto *mh_837 = buffer.data(mh + 837);
    const auto *mh_838 = buffer.data(mh + 838);
    const auto *mh_839 = buffer.data(mh + 839);
    const auto *mh_840 = buffer.data(mh + 840);
    const auto *mh_841 = buffer.data(mh + 841);
    const auto *mh_842 = buffer.data(mh + 842);
    const auto *mh_843 = buffer.data(mh + 843);
    const auto *mh_844 = buffer.data(mh + 844);
    const auto *mh_845 = buffer.data(mh + 845);
    const auto *mh_846 = buffer.data(mh + 846);
    const auto *mh_847 = buffer.data(mh + 847);
    const auto *mh_848 = buffer.data(mh + 848);
    const auto *mh_849 = buffer.data(mh + 849);
    const auto *mh_850 = buffer.data(mh + 850);
    const auto *mh_851 = buffer.data(mh + 851);
    const auto *mh_852 = buffer.data(mh + 852);
    const auto *mh_853 = buffer.data(mh + 853);
    const auto *mh_854 = buffer.data(mh + 854);
    const auto *mh_855 = buffer.data(mh + 855);
    const auto *mh_856 = buffer.data(mh + 856);
    const auto *mh_857 = buffer.data(mh + 857);
    const auto *mh_858 = buffer.data(mh + 858);
    const auto *mh_859 = buffer.data(mh + 859);
    const auto *mh_860 = buffer.data(mh + 860);
    const auto *mh_861 = buffer.data(mh + 861);
    const auto *mh_862 = buffer.data(mh + 862);
    const auto *mh_863 = buffer.data(mh + 863);
    const auto *mh_864 = buffer.data(mh + 864);
    const auto *mh_865 = buffer.data(mh + 865);
    const auto *mh_866 = buffer.data(mh + 866);
    const auto *mh_867 = buffer.data(mh + 867);
    const auto *mh_868 = buffer.data(mh + 868);
    const auto *mh_869 = buffer.data(mh + 869);
    const auto *mh_870 = buffer.data(mh + 870);
    const auto *mh_871 = buffer.data(mh + 871);
    const auto *mh_872 = buffer.data(mh + 872);
    const auto *mh_873 = buffer.data(mh + 873);
    const auto *mh_874 = buffer.data(mh + 874);
    const auto *mh_875 = buffer.data(mh + 875);
    const auto *mh_876 = buffer.data(mh + 876);
    const auto *mh_877 = buffer.data(mh + 877);
    const auto *mh_878 = buffer.data(mh + 878);
    const auto *mh_879 = buffer.data(mh + 879);
    const auto *mh_880 = buffer.data(mh + 880);
    const auto *mh_881 = buffer.data(mh + 881);
    const auto *mh_882 = buffer.data(mh + 882);
    const auto *mh_883 = buffer.data(mh + 883);
    const auto *mh_884 = buffer.data(mh + 884);
    const auto *mh_885 = buffer.data(mh + 885);
    const auto *mh_886 = buffer.data(mh + 886);
    const auto *mh_887 = buffer.data(mh + 887);
    const auto *mh_888 = buffer.data(mh + 888);
    const auto *mh_889 = buffer.data(mh + 889);
    const auto *mh_890 = buffer.data(mh + 890);
    const auto *mh_891 = buffer.data(mh + 891);
    const auto *mh_892 = buffer.data(mh + 892);
    const auto *mh_893 = buffer.data(mh + 893);
    const auto *mh_894 = buffer.data(mh + 894);
    const auto *mh_895 = buffer.data(mh + 895);
    const auto *mh_896 = buffer.data(mh + 896);
    const auto *mh_897 = buffer.data(mh + 897);
    const auto *mh_898 = buffer.data(mh + 898);
    const auto *mh_899 = buffer.data(mh + 899);
    const auto *mh_900 = buffer.data(mh + 900);
    const auto *mh_901 = buffer.data(mh + 901);
    const auto *mh_902 = buffer.data(mh + 902);
    const auto *mh_903 = buffer.data(mh + 903);
    const auto *mh_904 = buffer.data(mh + 904);
    const auto *mh_905 = buffer.data(mh + 905);
    const auto *mh_906 = buffer.data(mh + 906);
    const auto *mh_907 = buffer.data(mh + 907);
    const auto *mh_908 = buffer.data(mh + 908);
    const auto *mh_909 = buffer.data(mh + 909);
    const auto *mh_910 = buffer.data(mh + 910);
    const auto *mh_911 = buffer.data(mh + 911);
    const auto *mh_912 = buffer.data(mh + 912);
    const auto *mh_913 = buffer.data(mh + 913);
    const auto *mh_914 = buffer.data(mh + 914);
    const auto *mh_915 = buffer.data(mh + 915);
    const auto *mh_916 = buffer.data(mh + 916);
    const auto *mh_917 = buffer.data(mh + 917);
    const auto *mh_918 = buffer.data(mh + 918);
    const auto *mh_919 = buffer.data(mh + 919);
    const auto *mh_920 = buffer.data(mh + 920);
    const auto *mh_921 = buffer.data(mh + 921);
    const auto *mh_922 = buffer.data(mh + 922);
    const auto *mh_923 = buffer.data(mh + 923);
    const auto *mh_945 = buffer.data(mh + 945);
    const auto *mh_946 = buffer.data(mh + 946);
    const auto *mh_947 = buffer.data(mh + 947);
    const auto *mh_948 = buffer.data(mh + 948);
    const auto *mh_949 = buffer.data(mh + 949);
    const auto *mh_950 = buffer.data(mh + 950);
    const auto *mh_951 = buffer.data(mh + 951);
    const auto *mh_952 = buffer.data(mh + 952);
    const auto *mh_953 = buffer.data(mh + 953);
    const auto *mh_954 = buffer.data(mh + 954);
    const auto *mh_955 = buffer.data(mh + 955);
    const auto *mh_956 = buffer.data(mh + 956);
    const auto *mh_957 = buffer.data(mh + 957);
    const auto *mh_958 = buffer.data(mh + 958);
    const auto *mh_959 = buffer.data(mh + 959);
    const auto *mh_960 = buffer.data(mh + 960);
    const auto *mh_961 = buffer.data(mh + 961);
    const auto *mh_962 = buffer.data(mh + 962);
    const auto *mh_963 = buffer.data(mh + 963);
    const auto *mh_964 = buffer.data(mh + 964);
    const auto *mh_965 = buffer.data(mh + 965);
    const auto *mh_966 = buffer.data(mh + 966);
    const auto *mh_967 = buffer.data(mh + 967);
    const auto *mh_968 = buffer.data(mh + 968);
    const auto *mh_969 = buffer.data(mh + 969);
    const auto *mh_970 = buffer.data(mh + 970);
    const auto *mh_971 = buffer.data(mh + 971);
    const auto *mh_972 = buffer.data(mh + 972);
    const auto *mh_973 = buffer.data(mh + 973);
    const auto *mh_974 = buffer.data(mh + 974);
    const auto *mh_975 = buffer.data(mh + 975);
    const auto *mh_976 = buffer.data(mh + 976);
    const auto *mh_977 = buffer.data(mh + 977);
    const auto *mh_978 = buffer.data(mh + 978);
    const auto *mh_979 = buffer.data(mh + 979);
    const auto *mh_980 = buffer.data(mh + 980);
    const auto *mh_981 = buffer.data(mh + 981);
    const auto *mh_982 = buffer.data(mh + 982);
    const auto *mh_983 = buffer.data(mh + 983);
    const auto *mh_984 = buffer.data(mh + 984);
    const auto *mh_985 = buffer.data(mh + 985);
    const auto *mh_986 = buffer.data(mh + 986);
    const auto *mh_987 = buffer.data(mh + 987);
    const auto *mh_988 = buffer.data(mh + 988);
    const auto *mh_989 = buffer.data(mh + 989);
    const auto *mh_990 = buffer.data(mh + 990);
    const auto *mh_991 = buffer.data(mh + 991);
    const auto *mh_992 = buffer.data(mh + 992);

#pragma omp simd aligned(t_646, t_647, t_648, t_649, t_650, kh_499, kh_500, kh_501, kh_502, \
                         kh_503, mh_814, mh_815, mh_816, mh_817, \
                         mh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = -5.0 * kh_499[k]
                   + f_0 * mh_814[k];

        t_647[k] = -5.0 * kh_500[k]
                   + f_0 * mh_815[k];

        t_648[k] = -5.0 * kh_501[k]
                   + f_0 * mh_816[k];

        t_649[k] = -5.0 * kh_502[k]
                   + f_0 * mh_817[k];

        t_650[k] = -5.0 * kh_503[k]
                   + f_0 * mh_818[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, t_654, t_655, kh_504, kh_505, kh_506, kh_507, \
                         kh_508, mh_819, mh_820, mh_821, mh_822, \
                         mh_823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = -4.0 * kh_504[k]
                   + f_0 * mh_819[k];

        t_652[k] = -4.0 * kh_505[k]
                   + f_0 * mh_820[k];

        t_653[k] = -4.0 * kh_506[k]
                   + f_0 * mh_821[k];

        t_654[k] = -4.0 * kh_507[k]
                   + f_0 * mh_822[k];

        t_655[k] = -4.0 * kh_508[k]
                   + f_0 * mh_823[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, t_659, t_660, kh_509, kh_510, kh_511, kh_512, \
                         kh_513, mh_824, mh_825, mh_826, mh_827, \
                         mh_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = -4.0 * kh_509[k]
                   + f_0 * mh_824[k];

        t_657[k] = -4.0 * kh_510[k]
                   + f_0 * mh_825[k];

        t_658[k] = -4.0 * kh_511[k]
                   + f_0 * mh_826[k];

        t_659[k] = -4.0 * kh_512[k]
                   + f_0 * mh_827[k];

        t_660[k] = -4.0 * kh_513[k]
                   + f_0 * mh_828[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, kh_514, kh_515, kh_516, kh_517, \
                         kh_518, mh_829, mh_830, mh_831, mh_832, \
                         mh_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = -4.0 * kh_514[k]
                   + f_0 * mh_829[k];

        t_662[k] = -4.0 * kh_515[k]
                   + f_0 * mh_830[k];

        t_663[k] = -4.0 * kh_516[k]
                   + f_0 * mh_831[k];

        t_664[k] = -4.0 * kh_517[k]
                   + f_0 * mh_832[k];

        t_665[k] = -4.0 * kh_518[k]
                   + f_0 * mh_833[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, t_670, kh_519, kh_520, kh_521, kh_522, \
                         kh_523, mh_834, mh_835, mh_836, mh_837, \
                         mh_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = -4.0 * kh_519[k]
                   + f_0 * mh_834[k];

        t_667[k] = -4.0 * kh_520[k]
                   + f_0 * mh_835[k];

        t_668[k] = -4.0 * kh_521[k]
                   + f_0 * mh_836[k];

        t_669[k] = -4.0 * kh_522[k]
                   + f_0 * mh_837[k];

        t_670[k] = -4.0 * kh_523[k]
                   + f_0 * mh_838[k];
    }

#pragma omp simd aligned(t_671, t_672, t_673, t_674, t_675, kh_524, kh_525, kh_526, kh_527, \
                         kh_528, mh_839, mh_840, mh_841, mh_842, \
                         mh_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_671[k] = -4.0 * kh_524[k]
                   + f_0 * mh_839[k];

        t_672[k] = -3.0 * kh_525[k]
                   + f_0 * mh_840[k];

        t_673[k] = -3.0 * kh_526[k]
                   + f_0 * mh_841[k];

        t_674[k] = -3.0 * kh_527[k]
                   + f_0 * mh_842[k];

        t_675[k] = -3.0 * kh_528[k]
                   + f_0 * mh_843[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, t_680, kh_529, kh_530, kh_531, kh_532, \
                         kh_533, mh_844, mh_845, mh_846, mh_847, \
                         mh_848 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = -3.0 * kh_529[k]
                   + f_0 * mh_844[k];

        t_677[k] = -3.0 * kh_530[k]
                   + f_0 * mh_845[k];

        t_678[k] = -3.0 * kh_531[k]
                   + f_0 * mh_846[k];

        t_679[k] = -3.0 * kh_532[k]
                   + f_0 * mh_847[k];

        t_680[k] = -3.0 * kh_533[k]
                   + f_0 * mh_848[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, t_684, t_685, kh_534, kh_535, kh_536, kh_537, \
                         kh_538, mh_849, mh_850, mh_851, mh_852, \
                         mh_853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = -3.0 * kh_534[k]
                   + f_0 * mh_849[k];

        t_682[k] = -3.0 * kh_535[k]
                   + f_0 * mh_850[k];

        t_683[k] = -3.0 * kh_536[k]
                   + f_0 * mh_851[k];

        t_684[k] = -3.0 * kh_537[k]
                   + f_0 * mh_852[k];

        t_685[k] = -3.0 * kh_538[k]
                   + f_0 * mh_853[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, t_690, kh_539, kh_540, kh_541, kh_542, \
                         kh_543, mh_854, mh_855, mh_856, mh_857, \
                         mh_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = -3.0 * kh_539[k]
                   + f_0 * mh_854[k];

        t_687[k] = -3.0 * kh_540[k]
                   + f_0 * mh_855[k];

        t_688[k] = -3.0 * kh_541[k]
                   + f_0 * mh_856[k];

        t_689[k] = -3.0 * kh_542[k]
                   + f_0 * mh_857[k];

        t_690[k] = -3.0 * kh_543[k]
                   + f_0 * mh_858[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, t_695, kh_544, kh_545, kh_546, kh_547, \
                         kh_548, mh_859, mh_860, mh_861, mh_862, \
                         mh_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = -3.0 * kh_544[k]
                   + f_0 * mh_859[k];

        t_692[k] = -3.0 * kh_545[k]
                   + f_0 * mh_860[k];

        t_693[k] = -2.0 * kh_546[k]
                   + f_0 * mh_861[k];

        t_694[k] = -2.0 * kh_547[k]
                   + f_0 * mh_862[k];

        t_695[k] = -2.0 * kh_548[k]
                   + f_0 * mh_863[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, t_699, t_700, kh_549, kh_550, kh_551, kh_552, \
                         kh_553, mh_864, mh_865, mh_866, mh_867, \
                         mh_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = -2.0 * kh_549[k]
                   + f_0 * mh_864[k];

        t_697[k] = -2.0 * kh_550[k]
                   + f_0 * mh_865[k];

        t_698[k] = -2.0 * kh_551[k]
                   + f_0 * mh_866[k];

        t_699[k] = -2.0 * kh_552[k]
                   + f_0 * mh_867[k];

        t_700[k] = -2.0 * kh_553[k]
                   + f_0 * mh_868[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, t_705, kh_554, kh_555, kh_556, kh_557, \
                         kh_558, mh_869, mh_870, mh_871, mh_872, \
                         mh_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = -2.0 * kh_554[k]
                   + f_0 * mh_869[k];

        t_702[k] = -2.0 * kh_555[k]
                   + f_0 * mh_870[k];

        t_703[k] = -2.0 * kh_556[k]
                   + f_0 * mh_871[k];

        t_704[k] = -2.0 * kh_557[k]
                   + f_0 * mh_872[k];

        t_705[k] = -2.0 * kh_558[k]
                   + f_0 * mh_873[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, t_710, kh_559, kh_560, kh_561, kh_562, \
                         kh_563, mh_874, mh_875, mh_876, mh_877, \
                         mh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = -2.0 * kh_559[k]
                   + f_0 * mh_874[k];

        t_707[k] = -2.0 * kh_560[k]
                   + f_0 * mh_875[k];

        t_708[k] = -2.0 * kh_561[k]
                   + f_0 * mh_876[k];

        t_709[k] = -2.0 * kh_562[k]
                   + f_0 * mh_877[k];

        t_710[k] = -2.0 * kh_563[k]
                   + f_0 * mh_878[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, t_714, t_715, kh_564, kh_565, kh_566, kh_567, \
                         kh_568, mh_879, mh_880, mh_881, mh_882, \
                         mh_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = -2.0 * kh_564[k]
                   + f_0 * mh_879[k];

        t_712[k] = -2.0 * kh_565[k]
                   + f_0 * mh_880[k];

        t_713[k] = -2.0 * kh_566[k]
                   + f_0 * mh_881[k];

        t_714[k] = -kh_567[k]
                   + f_0 * mh_882[k];

        t_715[k] = -kh_568[k]
                   + f_0 * mh_883[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, t_720, kh_569, kh_570, kh_571, kh_572, \
                         kh_573, mh_884, mh_885, mh_886, mh_887, \
                         mh_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = -kh_569[k]
                   + f_0 * mh_884[k];

        t_717[k] = -kh_570[k]
                   + f_0 * mh_885[k];

        t_718[k] = -kh_571[k]
                   + f_0 * mh_886[k];

        t_719[k] = -kh_572[k]
                   + f_0 * mh_887[k];

        t_720[k] = -kh_573[k]
                   + f_0 * mh_888[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, t_725, kh_574, kh_575, kh_576, kh_577, \
                         kh_578, mh_889, mh_890, mh_891, mh_892, \
                         mh_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = -kh_574[k]
                   + f_0 * mh_889[k];

        t_722[k] = -kh_575[k]
                   + f_0 * mh_890[k];

        t_723[k] = -kh_576[k]
                   + f_0 * mh_891[k];

        t_724[k] = -kh_577[k]
                   + f_0 * mh_892[k];

        t_725[k] = -kh_578[k]
                   + f_0 * mh_893[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, t_730, kh_579, kh_580, kh_581, kh_582, \
                         kh_583, mh_894, mh_895, mh_896, mh_897, \
                         mh_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = -kh_579[k]
                   + f_0 * mh_894[k];

        t_727[k] = -kh_580[k]
                   + f_0 * mh_895[k];

        t_728[k] = -kh_581[k]
                   + f_0 * mh_896[k];

        t_729[k] = -kh_582[k]
                   + f_0 * mh_897[k];

        t_730[k] = -kh_583[k]
                   + f_0 * mh_898[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, t_734, t_735, t_736, kh_584, kh_585, kh_586, \
                         kh_587, mh_899, mh_900, mh_901, mh_902, mh_903, \
                         mh_904 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = -kh_584[k]
                   + f_0 * mh_899[k];

        t_732[k] = -kh_585[k]
                   + f_0 * mh_900[k];

        t_733[k] = -kh_586[k]
                   + f_0 * mh_901[k];

        t_734[k] = -kh_587[k]
                   + f_0 * mh_902[k];

        t_735[k] = f_0 * mh_903[k];

        t_736[k] = f_0 * mh_904[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, t_740, t_741, t_742, t_743, t_744, mh_905, \
                         mh_906, mh_907, mh_908, mh_909, mh_910, mh_911, \
                         mh_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_0 * mh_905[k];

        t_738[k] = f_0 * mh_906[k];

        t_739[k] = f_0 * mh_907[k];

        t_740[k] = f_0 * mh_908[k];

        t_741[k] = f_0 * mh_909[k];

        t_742[k] = f_0 * mh_910[k];

        t_743[k] = f_0 * mh_911[k];

        t_744[k] = f_0 * mh_912[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, t_750, t_751, t_752, mh_913, \
                         mh_914, mh_915, mh_916, mh_917, mh_918, mh_919, \
                         mh_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_0 * mh_913[k];

        t_746[k] = f_0 * mh_914[k];

        t_747[k] = f_0 * mh_915[k];

        t_748[k] = f_0 * mh_916[k];

        t_749[k] = f_0 * mh_917[k];

        t_750[k] = f_0 * mh_918[k];

        t_751[k] = f_0 * mh_919[k];

        t_752[k] = f_0 * mh_920[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, t_756, t_757, t_758, kh_588, kh_589, kh_590, \
                         mh_921, mh_922, mh_923, mh_945, mh_946, \
                         mh_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_0 * mh_921[k];

        t_754[k] = f_0 * mh_922[k];

        t_755[k] = f_0 * mh_923[k];

        t_756[k] = -8.0 * kh_588[k]
                   + f_0 * mh_945[k];

        t_757[k] = -8.0 * kh_589[k]
                   + f_0 * mh_946[k];

        t_758[k] = -8.0 * kh_590[k]
                   + f_0 * mh_947[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, t_762, t_763, kh_591, kh_592, kh_593, kh_594, \
                         kh_595, mh_948, mh_949, mh_950, mh_951, \
                         mh_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = -8.0 * kh_591[k]
                   + f_0 * mh_948[k];

        t_760[k] = -8.0 * kh_592[k]
                   + f_0 * mh_949[k];

        t_761[k] = -8.0 * kh_593[k]
                   + f_0 * mh_950[k];

        t_762[k] = -8.0 * kh_594[k]
                   + f_0 * mh_951[k];

        t_763[k] = -8.0 * kh_595[k]
                   + f_0 * mh_952[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, t_768, kh_596, kh_597, kh_598, kh_599, \
                         kh_600, mh_953, mh_954, mh_955, mh_956, \
                         mh_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = -8.0 * kh_596[k]
                   + f_0 * mh_953[k];

        t_765[k] = -8.0 * kh_597[k]
                   + f_0 * mh_954[k];

        t_766[k] = -8.0 * kh_598[k]
                   + f_0 * mh_955[k];

        t_767[k] = -8.0 * kh_599[k]
                   + f_0 * mh_956[k];

        t_768[k] = -8.0 * kh_600[k]
                   + f_0 * mh_957[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, t_772, t_773, kh_601, kh_602, kh_603, kh_604, \
                         kh_605, mh_958, mh_959, mh_960, mh_961, \
                         mh_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = -8.0 * kh_601[k]
                   + f_0 * mh_958[k];

        t_770[k] = -8.0 * kh_602[k]
                   + f_0 * mh_959[k];

        t_771[k] = -8.0 * kh_603[k]
                   + f_0 * mh_960[k];

        t_772[k] = -8.0 * kh_604[k]
                   + f_0 * mh_961[k];

        t_773[k] = -8.0 * kh_605[k]
                   + f_0 * mh_962[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, t_778, kh_606, kh_607, kh_608, kh_609, \
                         kh_610, mh_963, mh_964, mh_965, mh_966, \
                         mh_967 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = -8.0 * kh_606[k]
                   + f_0 * mh_963[k];

        t_775[k] = -8.0 * kh_607[k]
                   + f_0 * mh_964[k];

        t_776[k] = -8.0 * kh_608[k]
                   + f_0 * mh_965[k];

        t_777[k] = -7.0 * kh_609[k]
                   + f_0 * mh_966[k];

        t_778[k] = -7.0 * kh_610[k]
                   + f_0 * mh_967[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, t_782, t_783, kh_611, kh_612, kh_613, kh_614, \
                         kh_615, mh_968, mh_969, mh_970, mh_971, \
                         mh_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = -7.0 * kh_611[k]
                   + f_0 * mh_968[k];

        t_780[k] = -7.0 * kh_612[k]
                   + f_0 * mh_969[k];

        t_781[k] = -7.0 * kh_613[k]
                   + f_0 * mh_970[k];

        t_782[k] = -7.0 * kh_614[k]
                   + f_0 * mh_971[k];

        t_783[k] = -7.0 * kh_615[k]
                   + f_0 * mh_972[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, t_788, kh_616, kh_617, kh_618, kh_619, \
                         kh_620, mh_973, mh_974, mh_975, mh_976, \
                         mh_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = -7.0 * kh_616[k]
                   + f_0 * mh_973[k];

        t_785[k] = -7.0 * kh_617[k]
                   + f_0 * mh_974[k];

        t_786[k] = -7.0 * kh_618[k]
                   + f_0 * mh_975[k];

        t_787[k] = -7.0 * kh_619[k]
                   + f_0 * mh_976[k];

        t_788[k] = -7.0 * kh_620[k]
                   + f_0 * mh_977[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, t_792, t_793, kh_621, kh_622, kh_623, kh_624, \
                         kh_625, mh_978, mh_979, mh_980, mh_981, \
                         mh_982 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = -7.0 * kh_621[k]
                   + f_0 * mh_978[k];

        t_790[k] = -7.0 * kh_622[k]
                   + f_0 * mh_979[k];

        t_791[k] = -7.0 * kh_623[k]
                   + f_0 * mh_980[k];

        t_792[k] = -7.0 * kh_624[k]
                   + f_0 * mh_981[k];

        t_793[k] = -7.0 * kh_625[k]
                   + f_0 * mh_982[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, t_798, kh_626, kh_627, kh_628, kh_629, \
                         kh_630, mh_983, mh_984, mh_985, mh_986, \
                         mh_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = -7.0 * kh_626[k]
                   + f_0 * mh_983[k];

        t_795[k] = -7.0 * kh_627[k]
                   + f_0 * mh_984[k];

        t_796[k] = -7.0 * kh_628[k]
                   + f_0 * mh_985[k];

        t_797[k] = -7.0 * kh_629[k]
                   + f_0 * mh_986[k];

        t_798[k] = -6.0 * kh_630[k]
                   + f_0 * mh_987[k];
    }

#pragma omp simd aligned(t_799, t_800, t_801, t_802, t_803, kh_631, kh_632, kh_633, kh_634, \
                         kh_635, mh_988, mh_989, mh_990, mh_991, \
                         mh_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = -6.0 * kh_631[k]
                   + f_0 * mh_988[k];

        t_800[k] = -6.0 * kh_632[k]
                   + f_0 * mh_989[k];

        t_801[k] = -6.0 * kh_633[k]
                   + f_0 * mh_990[k];

        t_802[k] = -6.0 * kh_634[k]
                   + f_0 * mh_991[k];

        t_803[k] = -6.0 * kh_635[k]
                   + f_0 * mh_992[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_1_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_804 = buffer.data(target + 804);
    auto *t_805 = buffer.data(target + 805);
    auto *t_806 = buffer.data(target + 806);
    auto *t_807 = buffer.data(target + 807);
    auto *t_808 = buffer.data(target + 808);
    auto *t_809 = buffer.data(target + 809);
    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);
    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);
    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);
    auto *t_840 = buffer.data(target + 840);
    auto *t_841 = buffer.data(target + 841);
    auto *t_842 = buffer.data(target + 842);
    auto *t_843 = buffer.data(target + 843);
    auto *t_844 = buffer.data(target + 844);
    auto *t_845 = buffer.data(target + 845);
    auto *t_846 = buffer.data(target + 846);
    auto *t_847 = buffer.data(target + 847);
    auto *t_848 = buffer.data(target + 848);
    auto *t_849 = buffer.data(target + 849);
    auto *t_850 = buffer.data(target + 850);
    auto *t_851 = buffer.data(target + 851);
    auto *t_852 = buffer.data(target + 852);
    auto *t_853 = buffer.data(target + 853);
    auto *t_854 = buffer.data(target + 854);
    auto *t_855 = buffer.data(target + 855);
    auto *t_856 = buffer.data(target + 856);
    auto *t_857 = buffer.data(target + 857);
    auto *t_858 = buffer.data(target + 858);
    auto *t_859 = buffer.data(target + 859);
    auto *t_860 = buffer.data(target + 860);
    auto *t_861 = buffer.data(target + 861);
    auto *t_862 = buffer.data(target + 862);
    auto *t_863 = buffer.data(target + 863);
    auto *t_864 = buffer.data(target + 864);
    auto *t_865 = buffer.data(target + 865);
    auto *t_866 = buffer.data(target + 866);
    auto *t_867 = buffer.data(target + 867);
    auto *t_868 = buffer.data(target + 868);
    auto *t_869 = buffer.data(target + 869);
    auto *t_870 = buffer.data(target + 870);
    auto *t_871 = buffer.data(target + 871);
    auto *t_872 = buffer.data(target + 872);
    auto *t_873 = buffer.data(target + 873);
    auto *t_874 = buffer.data(target + 874);
    auto *t_875 = buffer.data(target + 875);
    auto *t_876 = buffer.data(target + 876);
    auto *t_877 = buffer.data(target + 877);
    auto *t_878 = buffer.data(target + 878);
    auto *t_879 = buffer.data(target + 879);
    auto *t_880 = buffer.data(target + 880);
    auto *t_881 = buffer.data(target + 881);
    auto *t_882 = buffer.data(target + 882);
    auto *t_883 = buffer.data(target + 883);
    auto *t_884 = buffer.data(target + 884);
    auto *t_885 = buffer.data(target + 885);
    auto *t_886 = buffer.data(target + 886);
    auto *t_887 = buffer.data(target + 887);
    auto *t_888 = buffer.data(target + 888);
    auto *t_889 = buffer.data(target + 889);
    auto *t_890 = buffer.data(target + 890);
    auto *t_891 = buffer.data(target + 891);
    auto *t_892 = buffer.data(target + 892);
    auto *t_893 = buffer.data(target + 893);
    auto *t_894 = buffer.data(target + 894);
    auto *t_895 = buffer.data(target + 895);
    auto *t_896 = buffer.data(target + 896);
    auto *t_897 = buffer.data(target + 897);
    auto *t_898 = buffer.data(target + 898);
    auto *t_899 = buffer.data(target + 899);
    auto *t_900 = buffer.data(target + 900);
    auto *t_901 = buffer.data(target + 901);
    auto *t_902 = buffer.data(target + 902);
    auto *t_903 = buffer.data(target + 903);
    auto *t_904 = buffer.data(target + 904);
    auto *t_905 = buffer.data(target + 905);
    auto *t_906 = buffer.data(target + 906);
    auto *t_907 = buffer.data(target + 907);
    auto *t_908 = buffer.data(target + 908);
    auto *t_909 = buffer.data(target + 909);
    auto *t_910 = buffer.data(target + 910);
    auto *t_911 = buffer.data(target + 911);
    auto *t_912 = buffer.data(target + 912);
    auto *t_913 = buffer.data(target + 913);
    auto *t_914 = buffer.data(target + 914);
    auto *t_915 = buffer.data(target + 915);
    auto *t_916 = buffer.data(target + 916);
    auto *t_917 = buffer.data(target + 917);
    auto *t_918 = buffer.data(target + 918);
    auto *t_919 = buffer.data(target + 919);
    auto *t_920 = buffer.data(target + 920);
    auto *t_921 = buffer.data(target + 921);
    auto *t_922 = buffer.data(target + 922);
    auto *t_923 = buffer.data(target + 923);
    auto *t_924 = buffer.data(target + 924);
    auto *t_925 = buffer.data(target + 925);
    auto *t_926 = buffer.data(target + 926);
    auto *t_927 = buffer.data(target + 927);
    auto *t_928 = buffer.data(target + 928);
    auto *t_929 = buffer.data(target + 929);
    auto *t_930 = buffer.data(target + 930);
    auto *t_931 = buffer.data(target + 931);
    auto *t_932 = buffer.data(target + 932);
    auto *t_933 = buffer.data(target + 933);
    auto *t_934 = buffer.data(target + 934);
    auto *t_935 = buffer.data(target + 935);
    auto *t_936 = buffer.data(target + 936);
    auto *t_937 = buffer.data(target + 937);
    auto *t_938 = buffer.data(target + 938);
    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);

    const auto *kh_636 = buffer.data(kh + 636);
    const auto *kh_637 = buffer.data(kh + 637);
    const auto *kh_638 = buffer.data(kh + 638);
    const auto *kh_639 = buffer.data(kh + 639);
    const auto *kh_640 = buffer.data(kh + 640);
    const auto *kh_641 = buffer.data(kh + 641);
    const auto *kh_642 = buffer.data(kh + 642);
    const auto *kh_643 = buffer.data(kh + 643);
    const auto *kh_644 = buffer.data(kh + 644);
    const auto *kh_645 = buffer.data(kh + 645);
    const auto *kh_646 = buffer.data(kh + 646);
    const auto *kh_647 = buffer.data(kh + 647);
    const auto *kh_648 = buffer.data(kh + 648);
    const auto *kh_649 = buffer.data(kh + 649);
    const auto *kh_650 = buffer.data(kh + 650);
    const auto *kh_651 = buffer.data(kh + 651);
    const auto *kh_652 = buffer.data(kh + 652);
    const auto *kh_653 = buffer.data(kh + 653);
    const auto *kh_654 = buffer.data(kh + 654);
    const auto *kh_655 = buffer.data(kh + 655);
    const auto *kh_656 = buffer.data(kh + 656);
    const auto *kh_657 = buffer.data(kh + 657);
    const auto *kh_658 = buffer.data(kh + 658);
    const auto *kh_659 = buffer.data(kh + 659);
    const auto *kh_660 = buffer.data(kh + 660);
    const auto *kh_661 = buffer.data(kh + 661);
    const auto *kh_662 = buffer.data(kh + 662);
    const auto *kh_663 = buffer.data(kh + 663);
    const auto *kh_664 = buffer.data(kh + 664);
    const auto *kh_665 = buffer.data(kh + 665);
    const auto *kh_666 = buffer.data(kh + 666);
    const auto *kh_667 = buffer.data(kh + 667);
    const auto *kh_668 = buffer.data(kh + 668);
    const auto *kh_669 = buffer.data(kh + 669);
    const auto *kh_670 = buffer.data(kh + 670);
    const auto *kh_671 = buffer.data(kh + 671);
    const auto *kh_672 = buffer.data(kh + 672);
    const auto *kh_673 = buffer.data(kh + 673);
    const auto *kh_674 = buffer.data(kh + 674);
    const auto *kh_675 = buffer.data(kh + 675);
    const auto *kh_676 = buffer.data(kh + 676);
    const auto *kh_677 = buffer.data(kh + 677);
    const auto *kh_678 = buffer.data(kh + 678);
    const auto *kh_679 = buffer.data(kh + 679);
    const auto *kh_680 = buffer.data(kh + 680);
    const auto *kh_681 = buffer.data(kh + 681);
    const auto *kh_682 = buffer.data(kh + 682);
    const auto *kh_683 = buffer.data(kh + 683);
    const auto *kh_684 = buffer.data(kh + 684);
    const auto *kh_685 = buffer.data(kh + 685);
    const auto *kh_686 = buffer.data(kh + 686);
    const auto *kh_687 = buffer.data(kh + 687);
    const auto *kh_688 = buffer.data(kh + 688);
    const auto *kh_689 = buffer.data(kh + 689);
    const auto *kh_690 = buffer.data(kh + 690);
    const auto *kh_691 = buffer.data(kh + 691);
    const auto *kh_692 = buffer.data(kh + 692);
    const auto *kh_693 = buffer.data(kh + 693);
    const auto *kh_694 = buffer.data(kh + 694);
    const auto *kh_695 = buffer.data(kh + 695);
    const auto *kh_696 = buffer.data(kh + 696);
    const auto *kh_697 = buffer.data(kh + 697);
    const auto *kh_698 = buffer.data(kh + 698);
    const auto *kh_699 = buffer.data(kh + 699);
    const auto *kh_700 = buffer.data(kh + 700);
    const auto *kh_701 = buffer.data(kh + 701);
    const auto *kh_702 = buffer.data(kh + 702);
    const auto *kh_703 = buffer.data(kh + 703);
    const auto *kh_704 = buffer.data(kh + 704);
    const auto *kh_705 = buffer.data(kh + 705);
    const auto *kh_706 = buffer.data(kh + 706);
    const auto *kh_707 = buffer.data(kh + 707);
    const auto *kh_708 = buffer.data(kh + 708);
    const auto *kh_709 = buffer.data(kh + 709);
    const auto *kh_710 = buffer.data(kh + 710);
    const auto *kh_711 = buffer.data(kh + 711);
    const auto *kh_712 = buffer.data(kh + 712);
    const auto *kh_713 = buffer.data(kh + 713);
    const auto *kh_714 = buffer.data(kh + 714);
    const auto *kh_715 = buffer.data(kh + 715);
    const auto *kh_716 = buffer.data(kh + 716);
    const auto *kh_717 = buffer.data(kh + 717);
    const auto *kh_718 = buffer.data(kh + 718);
    const auto *kh_719 = buffer.data(kh + 719);
    const auto *kh_720 = buffer.data(kh + 720);
    const auto *kh_721 = buffer.data(kh + 721);
    const auto *kh_722 = buffer.data(kh + 722);
    const auto *kh_723 = buffer.data(kh + 723);
    const auto *kh_724 = buffer.data(kh + 724);
    const auto *kh_725 = buffer.data(kh + 725);
    const auto *kh_726 = buffer.data(kh + 726);
    const auto *kh_727 = buffer.data(kh + 727);
    const auto *kh_728 = buffer.data(kh + 728);
    const auto *kh_729 = buffer.data(kh + 729);
    const auto *kh_730 = buffer.data(kh + 730);
    const auto *kh_731 = buffer.data(kh + 731);
    const auto *kh_732 = buffer.data(kh + 732);
    const auto *kh_733 = buffer.data(kh + 733);
    const auto *kh_734 = buffer.data(kh + 734);
    const auto *kh_735 = buffer.data(kh + 735);
    const auto *kh_736 = buffer.data(kh + 736);
    const auto *kh_737 = buffer.data(kh + 737);
    const auto *kh_738 = buffer.data(kh + 738);
    const auto *kh_739 = buffer.data(kh + 739);
    const auto *kh_740 = buffer.data(kh + 740);
    const auto *kh_741 = buffer.data(kh + 741);
    const auto *kh_742 = buffer.data(kh + 742);
    const auto *kh_743 = buffer.data(kh + 743);
    const auto *kh_744 = buffer.data(kh + 744);
    const auto *kh_745 = buffer.data(kh + 745);
    const auto *kh_746 = buffer.data(kh + 746);
    const auto *kh_747 = buffer.data(kh + 747);
    const auto *kh_748 = buffer.data(kh + 748);
    const auto *kh_749 = buffer.data(kh + 749);
    const auto *kh_750 = buffer.data(kh + 750);
    const auto *kh_751 = buffer.data(kh + 751);
    const auto *kh_752 = buffer.data(kh + 752);
    const auto *kh_753 = buffer.data(kh + 753);
    const auto *kh_754 = buffer.data(kh + 754);
    const auto *kh_755 = buffer.data(kh + 755);

    const auto *mh_993 = buffer.data(mh + 993);
    const auto *mh_994 = buffer.data(mh + 994);
    const auto *mh_995 = buffer.data(mh + 995);
    const auto *mh_996 = buffer.data(mh + 996);
    const auto *mh_997 = buffer.data(mh + 997);
    const auto *mh_998 = buffer.data(mh + 998);
    const auto *mh_999 = buffer.data(mh + 999);
    const auto *mh_1000 = buffer.data(mh + 1000);
    const auto *mh_1001 = buffer.data(mh + 1001);
    const auto *mh_1002 = buffer.data(mh + 1002);
    const auto *mh_1003 = buffer.data(mh + 1003);
    const auto *mh_1004 = buffer.data(mh + 1004);
    const auto *mh_1005 = buffer.data(mh + 1005);
    const auto *mh_1006 = buffer.data(mh + 1006);
    const auto *mh_1007 = buffer.data(mh + 1007);
    const auto *mh_1008 = buffer.data(mh + 1008);
    const auto *mh_1009 = buffer.data(mh + 1009);
    const auto *mh_1010 = buffer.data(mh + 1010);
    const auto *mh_1011 = buffer.data(mh + 1011);
    const auto *mh_1012 = buffer.data(mh + 1012);
    const auto *mh_1013 = buffer.data(mh + 1013);
    const auto *mh_1014 = buffer.data(mh + 1014);
    const auto *mh_1015 = buffer.data(mh + 1015);
    const auto *mh_1016 = buffer.data(mh + 1016);
    const auto *mh_1017 = buffer.data(mh + 1017);
    const auto *mh_1018 = buffer.data(mh + 1018);
    const auto *mh_1019 = buffer.data(mh + 1019);
    const auto *mh_1020 = buffer.data(mh + 1020);
    const auto *mh_1021 = buffer.data(mh + 1021);
    const auto *mh_1022 = buffer.data(mh + 1022);
    const auto *mh_1023 = buffer.data(mh + 1023);
    const auto *mh_1024 = buffer.data(mh + 1024);
    const auto *mh_1025 = buffer.data(mh + 1025);
    const auto *mh_1026 = buffer.data(mh + 1026);
    const auto *mh_1027 = buffer.data(mh + 1027);
    const auto *mh_1028 = buffer.data(mh + 1028);
    const auto *mh_1029 = buffer.data(mh + 1029);
    const auto *mh_1030 = buffer.data(mh + 1030);
    const auto *mh_1031 = buffer.data(mh + 1031);
    const auto *mh_1032 = buffer.data(mh + 1032);
    const auto *mh_1033 = buffer.data(mh + 1033);
    const auto *mh_1034 = buffer.data(mh + 1034);
    const auto *mh_1035 = buffer.data(mh + 1035);
    const auto *mh_1036 = buffer.data(mh + 1036);
    const auto *mh_1037 = buffer.data(mh + 1037);
    const auto *mh_1038 = buffer.data(mh + 1038);
    const auto *mh_1039 = buffer.data(mh + 1039);
    const auto *mh_1040 = buffer.data(mh + 1040);
    const auto *mh_1041 = buffer.data(mh + 1041);
    const auto *mh_1042 = buffer.data(mh + 1042);
    const auto *mh_1043 = buffer.data(mh + 1043);
    const auto *mh_1044 = buffer.data(mh + 1044);
    const auto *mh_1045 = buffer.data(mh + 1045);
    const auto *mh_1046 = buffer.data(mh + 1046);
    const auto *mh_1047 = buffer.data(mh + 1047);
    const auto *mh_1048 = buffer.data(mh + 1048);
    const auto *mh_1049 = buffer.data(mh + 1049);
    const auto *mh_1050 = buffer.data(mh + 1050);
    const auto *mh_1051 = buffer.data(mh + 1051);
    const auto *mh_1052 = buffer.data(mh + 1052);
    const auto *mh_1053 = buffer.data(mh + 1053);
    const auto *mh_1054 = buffer.data(mh + 1054);
    const auto *mh_1055 = buffer.data(mh + 1055);
    const auto *mh_1056 = buffer.data(mh + 1056);
    const auto *mh_1057 = buffer.data(mh + 1057);
    const auto *mh_1058 = buffer.data(mh + 1058);
    const auto *mh_1059 = buffer.data(mh + 1059);
    const auto *mh_1060 = buffer.data(mh + 1060);
    const auto *mh_1061 = buffer.data(mh + 1061);
    const auto *mh_1062 = buffer.data(mh + 1062);
    const auto *mh_1063 = buffer.data(mh + 1063);
    const auto *mh_1064 = buffer.data(mh + 1064);
    const auto *mh_1065 = buffer.data(mh + 1065);
    const auto *mh_1066 = buffer.data(mh + 1066);
    const auto *mh_1067 = buffer.data(mh + 1067);
    const auto *mh_1068 = buffer.data(mh + 1068);
    const auto *mh_1069 = buffer.data(mh + 1069);
    const auto *mh_1070 = buffer.data(mh + 1070);
    const auto *mh_1071 = buffer.data(mh + 1071);
    const auto *mh_1072 = buffer.data(mh + 1072);
    const auto *mh_1073 = buffer.data(mh + 1073);
    const auto *mh_1074 = buffer.data(mh + 1074);
    const auto *mh_1075 = buffer.data(mh + 1075);
    const auto *mh_1076 = buffer.data(mh + 1076);
    const auto *mh_1077 = buffer.data(mh + 1077);
    const auto *mh_1078 = buffer.data(mh + 1078);
    const auto *mh_1079 = buffer.data(mh + 1079);
    const auto *mh_1080 = buffer.data(mh + 1080);
    const auto *mh_1081 = buffer.data(mh + 1081);
    const auto *mh_1082 = buffer.data(mh + 1082);
    const auto *mh_1083 = buffer.data(mh + 1083);
    const auto *mh_1084 = buffer.data(mh + 1084);
    const auto *mh_1085 = buffer.data(mh + 1085);
    const auto *mh_1086 = buffer.data(mh + 1086);
    const auto *mh_1087 = buffer.data(mh + 1087);
    const auto *mh_1088 = buffer.data(mh + 1088);
    const auto *mh_1089 = buffer.data(mh + 1089);
    const auto *mh_1090 = buffer.data(mh + 1090);
    const auto *mh_1091 = buffer.data(mh + 1091);
    const auto *mh_1092 = buffer.data(mh + 1092);
    const auto *mh_1093 = buffer.data(mh + 1093);
    const auto *mh_1094 = buffer.data(mh + 1094);
    const auto *mh_1095 = buffer.data(mh + 1095);
    const auto *mh_1096 = buffer.data(mh + 1096);
    const auto *mh_1097 = buffer.data(mh + 1097);
    const auto *mh_1098 = buffer.data(mh + 1098);
    const auto *mh_1099 = buffer.data(mh + 1099);
    const auto *mh_1100 = buffer.data(mh + 1100);
    const auto *mh_1101 = buffer.data(mh + 1101);
    const auto *mh_1102 = buffer.data(mh + 1102);
    const auto *mh_1103 = buffer.data(mh + 1103);
    const auto *mh_1104 = buffer.data(mh + 1104);
    const auto *mh_1105 = buffer.data(mh + 1105);
    const auto *mh_1106 = buffer.data(mh + 1106);
    const auto *mh_1107 = buffer.data(mh + 1107);
    const auto *mh_1108 = buffer.data(mh + 1108);
    const auto *mh_1109 = buffer.data(mh + 1109);
    const auto *mh_1110 = buffer.data(mh + 1110);
    const auto *mh_1111 = buffer.data(mh + 1111);
    const auto *mh_1112 = buffer.data(mh + 1112);
    const auto *mh_1113 = buffer.data(mh + 1113);
    const auto *mh_1114 = buffer.data(mh + 1114);
    const auto *mh_1115 = buffer.data(mh + 1115);
    const auto *mh_1116 = buffer.data(mh + 1116);
    const auto *mh_1117 = buffer.data(mh + 1117);
    const auto *mh_1118 = buffer.data(mh + 1118);
    const auto *mh_1119 = buffer.data(mh + 1119);
    const auto *mh_1120 = buffer.data(mh + 1120);
    const auto *mh_1121 = buffer.data(mh + 1121);
    const auto *mh_1122 = buffer.data(mh + 1122);
    const auto *mh_1123 = buffer.data(mh + 1123);
    const auto *mh_1124 = buffer.data(mh + 1124);
    const auto *mh_1125 = buffer.data(mh + 1125);
    const auto *mh_1126 = buffer.data(mh + 1126);
    const auto *mh_1127 = buffer.data(mh + 1127);
    const auto *mh_1128 = buffer.data(mh + 1128);
    const auto *mh_1129 = buffer.data(mh + 1129);
    const auto *mh_1130 = buffer.data(mh + 1130);
    const auto *mh_1131 = buffer.data(mh + 1131);
    const auto *mh_1132 = buffer.data(mh + 1132);
    const auto *mh_1133 = buffer.data(mh + 1133);

#pragma omp simd aligned(t_804, t_805, t_806, t_807, t_808, kh_636, kh_637, kh_638, kh_639, \
                         kh_640, mh_993, mh_994, mh_995, mh_996, \
                         mh_997 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = -6.0 * kh_636[k]
                   + f_0 * mh_993[k];

        t_805[k] = -6.0 * kh_637[k]
                   + f_0 * mh_994[k];

        t_806[k] = -6.0 * kh_638[k]
                   + f_0 * mh_995[k];

        t_807[k] = -6.0 * kh_639[k]
                   + f_0 * mh_996[k];

        t_808[k] = -6.0 * kh_640[k]
                   + f_0 * mh_997[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, t_812, t_813, kh_641, kh_642, kh_643, kh_644, \
                         kh_645, mh_998, mh_999, mh_1000, mh_1001, \
                         mh_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = -6.0 * kh_641[k]
                   + f_0 * mh_998[k];

        t_810[k] = -6.0 * kh_642[k]
                   + f_0 * mh_999[k];

        t_811[k] = -6.0 * kh_643[k]
                   + f_0 * mh_1000[k];

        t_812[k] = -6.0 * kh_644[k]
                   + f_0 * mh_1001[k];

        t_813[k] = -6.0 * kh_645[k]
                   + f_0 * mh_1002[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, t_817, t_818, kh_646, kh_647, kh_648, kh_649, \
                         kh_650, mh_1003, mh_1004, mh_1005, mh_1006, \
                         mh_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = -6.0 * kh_646[k]
                   + f_0 * mh_1003[k];

        t_815[k] = -6.0 * kh_647[k]
                   + f_0 * mh_1004[k];

        t_816[k] = -6.0 * kh_648[k]
                   + f_0 * mh_1005[k];

        t_817[k] = -6.0 * kh_649[k]
                   + f_0 * mh_1006[k];

        t_818[k] = -6.0 * kh_650[k]
                   + f_0 * mh_1007[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, t_822, t_823, kh_651, kh_652, kh_653, kh_654, \
                         kh_655, mh_1008, mh_1009, mh_1010, mh_1011, \
                         mh_1012 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = -5.0 * kh_651[k]
                   + f_0 * mh_1008[k];

        t_820[k] = -5.0 * kh_652[k]
                   + f_0 * mh_1009[k];

        t_821[k] = -5.0 * kh_653[k]
                   + f_0 * mh_1010[k];

        t_822[k] = -5.0 * kh_654[k]
                   + f_0 * mh_1011[k];

        t_823[k] = -5.0 * kh_655[k]
                   + f_0 * mh_1012[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, t_827, t_828, kh_656, kh_657, kh_658, kh_659, \
                         kh_660, mh_1013, mh_1014, mh_1015, mh_1016, \
                         mh_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = -5.0 * kh_656[k]
                   + f_0 * mh_1013[k];

        t_825[k] = -5.0 * kh_657[k]
                   + f_0 * mh_1014[k];

        t_826[k] = -5.0 * kh_658[k]
                   + f_0 * mh_1015[k];

        t_827[k] = -5.0 * kh_659[k]
                   + f_0 * mh_1016[k];

        t_828[k] = -5.0 * kh_660[k]
                   + f_0 * mh_1017[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, t_833, kh_661, kh_662, kh_663, kh_664, \
                         kh_665, mh_1018, mh_1019, mh_1020, mh_1021, \
                         mh_1022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = -5.0 * kh_661[k]
                   + f_0 * mh_1018[k];

        t_830[k] = -5.0 * kh_662[k]
                   + f_0 * mh_1019[k];

        t_831[k] = -5.0 * kh_663[k]
                   + f_0 * mh_1020[k];

        t_832[k] = -5.0 * kh_664[k]
                   + f_0 * mh_1021[k];

        t_833[k] = -5.0 * kh_665[k]
                   + f_0 * mh_1022[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, t_837, t_838, kh_666, kh_667, kh_668, kh_669, \
                         kh_670, mh_1023, mh_1024, mh_1025, mh_1026, \
                         mh_1027 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = -5.0 * kh_666[k]
                   + f_0 * mh_1023[k];

        t_835[k] = -5.0 * kh_667[k]
                   + f_0 * mh_1024[k];

        t_836[k] = -5.0 * kh_668[k]
                   + f_0 * mh_1025[k];

        t_837[k] = -5.0 * kh_669[k]
                   + f_0 * mh_1026[k];

        t_838[k] = -5.0 * kh_670[k]
                   + f_0 * mh_1027[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, t_842, t_843, kh_671, kh_672, kh_673, kh_674, \
                         kh_675, mh_1028, mh_1029, mh_1030, mh_1031, \
                         mh_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = -5.0 * kh_671[k]
                   + f_0 * mh_1028[k];

        t_840[k] = -4.0 * kh_672[k]
                   + f_0 * mh_1029[k];

        t_841[k] = -4.0 * kh_673[k]
                   + f_0 * mh_1030[k];

        t_842[k] = -4.0 * kh_674[k]
                   + f_0 * mh_1031[k];

        t_843[k] = -4.0 * kh_675[k]
                   + f_0 * mh_1032[k];
    }

#pragma omp simd aligned(t_844, t_845, t_846, t_847, t_848, kh_676, kh_677, kh_678, kh_679, \
                         kh_680, mh_1033, mh_1034, mh_1035, mh_1036, \
                         mh_1037 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_844[k] = -4.0 * kh_676[k]
                   + f_0 * mh_1033[k];

        t_845[k] = -4.0 * kh_677[k]
                   + f_0 * mh_1034[k];

        t_846[k] = -4.0 * kh_678[k]
                   + f_0 * mh_1035[k];

        t_847[k] = -4.0 * kh_679[k]
                   + f_0 * mh_1036[k];

        t_848[k] = -4.0 * kh_680[k]
                   + f_0 * mh_1037[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, t_853, kh_681, kh_682, kh_683, kh_684, \
                         kh_685, mh_1038, mh_1039, mh_1040, mh_1041, \
                         mh_1042 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = -4.0 * kh_681[k]
                   + f_0 * mh_1038[k];

        t_850[k] = -4.0 * kh_682[k]
                   + f_0 * mh_1039[k];

        t_851[k] = -4.0 * kh_683[k]
                   + f_0 * mh_1040[k];

        t_852[k] = -4.0 * kh_684[k]
                   + f_0 * mh_1041[k];

        t_853[k] = -4.0 * kh_685[k]
                   + f_0 * mh_1042[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, t_858, kh_686, kh_687, kh_688, kh_689, \
                         kh_690, mh_1043, mh_1044, mh_1045, mh_1046, \
                         mh_1047 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = -4.0 * kh_686[k]
                   + f_0 * mh_1043[k];

        t_855[k] = -4.0 * kh_687[k]
                   + f_0 * mh_1044[k];

        t_856[k] = -4.0 * kh_688[k]
                   + f_0 * mh_1045[k];

        t_857[k] = -4.0 * kh_689[k]
                   + f_0 * mh_1046[k];

        t_858[k] = -4.0 * kh_690[k]
                   + f_0 * mh_1047[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, t_862, t_863, kh_691, kh_692, kh_693, kh_694, \
                         kh_695, mh_1048, mh_1049, mh_1050, mh_1051, \
                         mh_1052 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = -4.0 * kh_691[k]
                   + f_0 * mh_1048[k];

        t_860[k] = -4.0 * kh_692[k]
                   + f_0 * mh_1049[k];

        t_861[k] = -3.0 * kh_693[k]
                   + f_0 * mh_1050[k];

        t_862[k] = -3.0 * kh_694[k]
                   + f_0 * mh_1051[k];

        t_863[k] = -3.0 * kh_695[k]
                   + f_0 * mh_1052[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, t_867, t_868, kh_696, kh_697, kh_698, kh_699, \
                         kh_700, mh_1053, mh_1054, mh_1055, mh_1056, \
                         mh_1057 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = -3.0 * kh_696[k]
                   + f_0 * mh_1053[k];

        t_865[k] = -3.0 * kh_697[k]
                   + f_0 * mh_1054[k];

        t_866[k] = -3.0 * kh_698[k]
                   + f_0 * mh_1055[k];

        t_867[k] = -3.0 * kh_699[k]
                   + f_0 * mh_1056[k];

        t_868[k] = -3.0 * kh_700[k]
                   + f_0 * mh_1057[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, t_873, kh_701, kh_702, kh_703, kh_704, \
                         kh_705, mh_1058, mh_1059, mh_1060, mh_1061, \
                         mh_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = -3.0 * kh_701[k]
                   + f_0 * mh_1058[k];

        t_870[k] = -3.0 * kh_702[k]
                   + f_0 * mh_1059[k];

        t_871[k] = -3.0 * kh_703[k]
                   + f_0 * mh_1060[k];

        t_872[k] = -3.0 * kh_704[k]
                   + f_0 * mh_1061[k];

        t_873[k] = -3.0 * kh_705[k]
                   + f_0 * mh_1062[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, t_878, kh_706, kh_707, kh_708, kh_709, \
                         kh_710, mh_1063, mh_1064, mh_1065, mh_1066, \
                         mh_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = -3.0 * kh_706[k]
                   + f_0 * mh_1063[k];

        t_875[k] = -3.0 * kh_707[k]
                   + f_0 * mh_1064[k];

        t_876[k] = -3.0 * kh_708[k]
                   + f_0 * mh_1065[k];

        t_877[k] = -3.0 * kh_709[k]
                   + f_0 * mh_1066[k];

        t_878[k] = -3.0 * kh_710[k]
                   + f_0 * mh_1067[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, t_882, t_883, kh_711, kh_712, kh_713, kh_714, \
                         kh_715, mh_1068, mh_1069, mh_1070, mh_1071, \
                         mh_1072 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = -3.0 * kh_711[k]
                   + f_0 * mh_1068[k];

        t_880[k] = -3.0 * kh_712[k]
                   + f_0 * mh_1069[k];

        t_881[k] = -3.0 * kh_713[k]
                   + f_0 * mh_1070[k];

        t_882[k] = -2.0 * kh_714[k]
                   + f_0 * mh_1071[k];

        t_883[k] = -2.0 * kh_715[k]
                   + f_0 * mh_1072[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, t_888, kh_716, kh_717, kh_718, kh_719, \
                         kh_720, mh_1073, mh_1074, mh_1075, mh_1076, \
                         mh_1077 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = -2.0 * kh_716[k]
                   + f_0 * mh_1073[k];

        t_885[k] = -2.0 * kh_717[k]
                   + f_0 * mh_1074[k];

        t_886[k] = -2.0 * kh_718[k]
                   + f_0 * mh_1075[k];

        t_887[k] = -2.0 * kh_719[k]
                   + f_0 * mh_1076[k];

        t_888[k] = -2.0 * kh_720[k]
                   + f_0 * mh_1077[k];
    }

#pragma omp simd aligned(t_889, t_890, t_891, t_892, t_893, kh_721, kh_722, kh_723, kh_724, \
                         kh_725, mh_1078, mh_1079, mh_1080, mh_1081, \
                         mh_1082 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = -2.0 * kh_721[k]
                   + f_0 * mh_1078[k];

        t_890[k] = -2.0 * kh_722[k]
                   + f_0 * mh_1079[k];

        t_891[k] = -2.0 * kh_723[k]
                   + f_0 * mh_1080[k];

        t_892[k] = -2.0 * kh_724[k]
                   + f_0 * mh_1081[k];

        t_893[k] = -2.0 * kh_725[k]
                   + f_0 * mh_1082[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, t_898, kh_726, kh_727, kh_728, kh_729, \
                         kh_730, mh_1083, mh_1084, mh_1085, mh_1086, \
                         mh_1087 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = -2.0 * kh_726[k]
                   + f_0 * mh_1083[k];

        t_895[k] = -2.0 * kh_727[k]
                   + f_0 * mh_1084[k];

        t_896[k] = -2.0 * kh_728[k]
                   + f_0 * mh_1085[k];

        t_897[k] = -2.0 * kh_729[k]
                   + f_0 * mh_1086[k];

        t_898[k] = -2.0 * kh_730[k]
                   + f_0 * mh_1087[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, t_903, kh_731, kh_732, kh_733, kh_734, \
                         kh_735, mh_1088, mh_1089, mh_1090, mh_1091, \
                         mh_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = -2.0 * kh_731[k]
                   + f_0 * mh_1088[k];

        t_900[k] = -2.0 * kh_732[k]
                   + f_0 * mh_1089[k];

        t_901[k] = -2.0 * kh_733[k]
                   + f_0 * mh_1090[k];

        t_902[k] = -2.0 * kh_734[k]
                   + f_0 * mh_1091[k];

        t_903[k] = -kh_735[k]
                   + f_0 * mh_1092[k];
    }

#pragma omp simd aligned(t_904, t_905, t_906, t_907, t_908, kh_736, kh_737, kh_738, kh_739, \
                         kh_740, mh_1093, mh_1094, mh_1095, mh_1096, \
                         mh_1097 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_904[k] = -kh_736[k]
                   + f_0 * mh_1093[k];

        t_905[k] = -kh_737[k]
                   + f_0 * mh_1094[k];

        t_906[k] = -kh_738[k]
                   + f_0 * mh_1095[k];

        t_907[k] = -kh_739[k]
                   + f_0 * mh_1096[k];

        t_908[k] = -kh_740[k]
                   + f_0 * mh_1097[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, t_912, t_913, kh_741, kh_742, kh_743, kh_744, \
                         kh_745, mh_1098, mh_1099, mh_1100, mh_1101, \
                         mh_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = -kh_741[k]
                   + f_0 * mh_1098[k];

        t_910[k] = -kh_742[k]
                   + f_0 * mh_1099[k];

        t_911[k] = -kh_743[k]
                   + f_0 * mh_1100[k];

        t_912[k] = -kh_744[k]
                   + f_0 * mh_1101[k];

        t_913[k] = -kh_745[k]
                   + f_0 * mh_1102[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, t_918, kh_746, kh_747, kh_748, kh_749, \
                         kh_750, mh_1103, mh_1104, mh_1105, mh_1106, \
                         mh_1107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = -kh_746[k]
                   + f_0 * mh_1103[k];

        t_915[k] = -kh_747[k]
                   + f_0 * mh_1104[k];

        t_916[k] = -kh_748[k]
                   + f_0 * mh_1105[k];

        t_917[k] = -kh_749[k]
                   + f_0 * mh_1106[k];

        t_918[k] = -kh_750[k]
                   + f_0 * mh_1107[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, t_922, t_923, kh_751, kh_752, kh_753, kh_754, \
                         kh_755, mh_1108, mh_1109, mh_1110, mh_1111, \
                         mh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = -kh_751[k]
                   + f_0 * mh_1108[k];

        t_920[k] = -kh_752[k]
                   + f_0 * mh_1109[k];

        t_921[k] = -kh_753[k]
                   + f_0 * mh_1110[k];

        t_922[k] = -kh_754[k]
                   + f_0 * mh_1111[k];

        t_923[k] = -kh_755[k]
                   + f_0 * mh_1112[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, t_928, t_929, t_930, t_931, mh_1113, \
                         mh_1114, mh_1115, mh_1116, mh_1117, mh_1118, mh_1119, \
                         mh_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_0 * mh_1113[k];

        t_925[k] = f_0 * mh_1114[k];

        t_926[k] = f_0 * mh_1115[k];

        t_927[k] = f_0 * mh_1116[k];

        t_928[k] = f_0 * mh_1117[k];

        t_929[k] = f_0 * mh_1118[k];

        t_930[k] = f_0 * mh_1119[k];

        t_931[k] = f_0 * mh_1120[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, t_935, t_936, t_937, t_938, t_939, mh_1121, \
                         mh_1122, mh_1123, mh_1124, mh_1125, mh_1126, mh_1127, \
                         mh_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_0 * mh_1121[k];

        t_933[k] = f_0 * mh_1122[k];

        t_934[k] = f_0 * mh_1123[k];

        t_935[k] = f_0 * mh_1124[k];

        t_936[k] = f_0 * mh_1125[k];

        t_937[k] = f_0 * mh_1126[k];

        t_938[k] = f_0 * mh_1127[k];

        t_939[k] = f_0 * mh_1128[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, t_944, mh_1129, mh_1130, mh_1131, \
                         mh_1132, mh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = f_0 * mh_1129[k];

        t_941[k] = f_0 * mh_1130[k];

        t_942[k] = f_0 * mh_1131[k];

        t_943[k] = f_0 * mh_1132[k];

        t_944[k] = f_0 * mh_1133[k];
    }
}

auto
compute_prim_geom_10_lh_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t kh, const size_t mh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_lh_electron_repulsion_1_piece0(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_1_piece1(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_1_piece2(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_1_piece3(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_1_piece4(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_1_piece5(buffer, target, kh, mh, ncols, alpha);
}

static auto
compute_prim_geom_10_lh_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
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
    auto *t_135 = buffer.data(target + 135);
    auto *t_136 = buffer.data(target + 136);
    auto *t_137 = buffer.data(target + 137);
    auto *t_138 = buffer.data(target + 138);
    auto *t_139 = buffer.data(target + 139);
    auto *t_140 = buffer.data(target + 140);
    auto *t_141 = buffer.data(target + 141);
    auto *t_142 = buffer.data(target + 142);
    auto *t_143 = buffer.data(target + 143);
    auto *t_144 = buffer.data(target + 144);
    auto *t_145 = buffer.data(target + 145);
    auto *t_146 = buffer.data(target + 146);
    auto *t_147 = buffer.data(target + 147);
    auto *t_148 = buffer.data(target + 148);
    auto *t_149 = buffer.data(target + 149);
    auto *t_150 = buffer.data(target + 150);
    auto *t_151 = buffer.data(target + 151);
    auto *t_152 = buffer.data(target + 152);
    auto *t_153 = buffer.data(target + 153);
    auto *t_154 = buffer.data(target + 154);
    auto *t_155 = buffer.data(target + 155);
    auto *t_156 = buffer.data(target + 156);
    auto *t_157 = buffer.data(target + 157);
    auto *t_158 = buffer.data(target + 158);
    auto *t_159 = buffer.data(target + 159);
    auto *t_160 = buffer.data(target + 160);
    auto *t_161 = buffer.data(target + 161);
    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);
    auto *t_168 = buffer.data(target + 168);
    auto *t_169 = buffer.data(target + 169);
    auto *t_170 = buffer.data(target + 170);
    auto *t_171 = buffer.data(target + 171);
    auto *t_172 = buffer.data(target + 172);
    auto *t_173 = buffer.data(target + 173);
    auto *t_174 = buffer.data(target + 174);
    auto *t_175 = buffer.data(target + 175);
    auto *t_176 = buffer.data(target + 176);

    const auto *kh_0 = buffer.data(kh + 0);
    const auto *kh_1 = buffer.data(kh + 1);
    const auto *kh_2 = buffer.data(kh + 2);
    const auto *kh_3 = buffer.data(kh + 3);
    const auto *kh_4 = buffer.data(kh + 4);
    const auto *kh_5 = buffer.data(kh + 5);
    const auto *kh_6 = buffer.data(kh + 6);
    const auto *kh_7 = buffer.data(kh + 7);
    const auto *kh_8 = buffer.data(kh + 8);
    const auto *kh_9 = buffer.data(kh + 9);
    const auto *kh_10 = buffer.data(kh + 10);
    const auto *kh_11 = buffer.data(kh + 11);
    const auto *kh_12 = buffer.data(kh + 12);
    const auto *kh_13 = buffer.data(kh + 13);
    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_15 = buffer.data(kh + 15);
    const auto *kh_16 = buffer.data(kh + 16);
    const auto *kh_17 = buffer.data(kh + 17);
    const auto *kh_18 = buffer.data(kh + 18);
    const auto *kh_19 = buffer.data(kh + 19);
    const auto *kh_20 = buffer.data(kh + 20);
    const auto *kh_21 = buffer.data(kh + 21);
    const auto *kh_22 = buffer.data(kh + 22);
    const auto *kh_23 = buffer.data(kh + 23);
    const auto *kh_24 = buffer.data(kh + 24);
    const auto *kh_25 = buffer.data(kh + 25);
    const auto *kh_26 = buffer.data(kh + 26);
    const auto *kh_27 = buffer.data(kh + 27);
    const auto *kh_28 = buffer.data(kh + 28);
    const auto *kh_29 = buffer.data(kh + 29);
    const auto *kh_30 = buffer.data(kh + 30);
    const auto *kh_31 = buffer.data(kh + 31);
    const auto *kh_32 = buffer.data(kh + 32);
    const auto *kh_33 = buffer.data(kh + 33);
    const auto *kh_34 = buffer.data(kh + 34);
    const auto *kh_35 = buffer.data(kh + 35);
    const auto *kh_36 = buffer.data(kh + 36);
    const auto *kh_37 = buffer.data(kh + 37);
    const auto *kh_38 = buffer.data(kh + 38);
    const auto *kh_39 = buffer.data(kh + 39);
    const auto *kh_40 = buffer.data(kh + 40);
    const auto *kh_41 = buffer.data(kh + 41);
    const auto *kh_42 = buffer.data(kh + 42);
    const auto *kh_43 = buffer.data(kh + 43);
    const auto *kh_44 = buffer.data(kh + 44);
    const auto *kh_45 = buffer.data(kh + 45);
    const auto *kh_46 = buffer.data(kh + 46);
    const auto *kh_47 = buffer.data(kh + 47);
    const auto *kh_48 = buffer.data(kh + 48);
    const auto *kh_49 = buffer.data(kh + 49);
    const auto *kh_50 = buffer.data(kh + 50);
    const auto *kh_51 = buffer.data(kh + 51);
    const auto *kh_52 = buffer.data(kh + 52);
    const auto *kh_53 = buffer.data(kh + 53);
    const auto *kh_54 = buffer.data(kh + 54);
    const auto *kh_55 = buffer.data(kh + 55);
    const auto *kh_56 = buffer.data(kh + 56);
    const auto *kh_57 = buffer.data(kh + 57);
    const auto *kh_58 = buffer.data(kh + 58);
    const auto *kh_59 = buffer.data(kh + 59);
    const auto *kh_60 = buffer.data(kh + 60);
    const auto *kh_61 = buffer.data(kh + 61);
    const auto *kh_62 = buffer.data(kh + 62);
    const auto *kh_63 = buffer.data(kh + 63);
    const auto *kh_64 = buffer.data(kh + 64);
    const auto *kh_65 = buffer.data(kh + 65);
    const auto *kh_66 = buffer.data(kh + 66);
    const auto *kh_67 = buffer.data(kh + 67);
    const auto *kh_68 = buffer.data(kh + 68);
    const auto *kh_69 = buffer.data(kh + 69);
    const auto *kh_70 = buffer.data(kh + 70);
    const auto *kh_71 = buffer.data(kh + 71);
    const auto *kh_72 = buffer.data(kh + 72);
    const auto *kh_73 = buffer.data(kh + 73);
    const auto *kh_74 = buffer.data(kh + 74);
    const auto *kh_75 = buffer.data(kh + 75);
    const auto *kh_76 = buffer.data(kh + 76);
    const auto *kh_77 = buffer.data(kh + 77);
    const auto *kh_78 = buffer.data(kh + 78);
    const auto *kh_79 = buffer.data(kh + 79);
    const auto *kh_80 = buffer.data(kh + 80);
    const auto *kh_81 = buffer.data(kh + 81);
    const auto *kh_82 = buffer.data(kh + 82);
    const auto *kh_83 = buffer.data(kh + 83);
    const auto *kh_84 = buffer.data(kh + 84);
    const auto *kh_85 = buffer.data(kh + 85);
    const auto *kh_86 = buffer.data(kh + 86);
    const auto *kh_87 = buffer.data(kh + 87);
    const auto *kh_88 = buffer.data(kh + 88);
    const auto *kh_89 = buffer.data(kh + 89);
    const auto *kh_90 = buffer.data(kh + 90);
    const auto *kh_91 = buffer.data(kh + 91);
    const auto *kh_92 = buffer.data(kh + 92);

    const auto *mh_42 = buffer.data(mh + 42);
    const auto *mh_43 = buffer.data(mh + 43);
    const auto *mh_44 = buffer.data(mh + 44);
    const auto *mh_45 = buffer.data(mh + 45);
    const auto *mh_46 = buffer.data(mh + 46);
    const auto *mh_47 = buffer.data(mh + 47);
    const auto *mh_48 = buffer.data(mh + 48);
    const auto *mh_49 = buffer.data(mh + 49);
    const auto *mh_50 = buffer.data(mh + 50);
    const auto *mh_51 = buffer.data(mh + 51);
    const auto *mh_52 = buffer.data(mh + 52);
    const auto *mh_53 = buffer.data(mh + 53);
    const auto *mh_54 = buffer.data(mh + 54);
    const auto *mh_55 = buffer.data(mh + 55);
    const auto *mh_56 = buffer.data(mh + 56);
    const auto *mh_57 = buffer.data(mh + 57);
    const auto *mh_58 = buffer.data(mh + 58);
    const auto *mh_59 = buffer.data(mh + 59);
    const auto *mh_60 = buffer.data(mh + 60);
    const auto *mh_61 = buffer.data(mh + 61);
    const auto *mh_62 = buffer.data(mh + 62);
    const auto *mh_84 = buffer.data(mh + 84);
    const auto *mh_85 = buffer.data(mh + 85);
    const auto *mh_86 = buffer.data(mh + 86);
    const auto *mh_87 = buffer.data(mh + 87);
    const auto *mh_88 = buffer.data(mh + 88);
    const auto *mh_89 = buffer.data(mh + 89);
    const auto *mh_90 = buffer.data(mh + 90);
    const auto *mh_91 = buffer.data(mh + 91);
    const auto *mh_92 = buffer.data(mh + 92);
    const auto *mh_93 = buffer.data(mh + 93);
    const auto *mh_94 = buffer.data(mh + 94);
    const auto *mh_95 = buffer.data(mh + 95);
    const auto *mh_96 = buffer.data(mh + 96);
    const auto *mh_97 = buffer.data(mh + 97);
    const auto *mh_98 = buffer.data(mh + 98);
    const auto *mh_99 = buffer.data(mh + 99);
    const auto *mh_100 = buffer.data(mh + 100);
    const auto *mh_101 = buffer.data(mh + 101);
    const auto *mh_102 = buffer.data(mh + 102);
    const auto *mh_103 = buffer.data(mh + 103);
    const auto *mh_104 = buffer.data(mh + 104);
    const auto *mh_105 = buffer.data(mh + 105);
    const auto *mh_106 = buffer.data(mh + 106);
    const auto *mh_107 = buffer.data(mh + 107);
    const auto *mh_108 = buffer.data(mh + 108);
    const auto *mh_109 = buffer.data(mh + 109);
    const auto *mh_110 = buffer.data(mh + 110);
    const auto *mh_111 = buffer.data(mh + 111);
    const auto *mh_112 = buffer.data(mh + 112);
    const auto *mh_113 = buffer.data(mh + 113);
    const auto *mh_114 = buffer.data(mh + 114);
    const auto *mh_115 = buffer.data(mh + 115);
    const auto *mh_116 = buffer.data(mh + 116);
    const auto *mh_117 = buffer.data(mh + 117);
    const auto *mh_118 = buffer.data(mh + 118);
    const auto *mh_119 = buffer.data(mh + 119);
    const auto *mh_120 = buffer.data(mh + 120);
    const auto *mh_121 = buffer.data(mh + 121);
    const auto *mh_122 = buffer.data(mh + 122);
    const auto *mh_123 = buffer.data(mh + 123);
    const auto *mh_124 = buffer.data(mh + 124);
    const auto *mh_125 = buffer.data(mh + 125);
    const auto *mh_147 = buffer.data(mh + 147);
    const auto *mh_148 = buffer.data(mh + 148);
    const auto *mh_149 = buffer.data(mh + 149);
    const auto *mh_150 = buffer.data(mh + 150);
    const auto *mh_151 = buffer.data(mh + 151);
    const auto *mh_152 = buffer.data(mh + 152);
    const auto *mh_153 = buffer.data(mh + 153);
    const auto *mh_154 = buffer.data(mh + 154);
    const auto *mh_155 = buffer.data(mh + 155);
    const auto *mh_156 = buffer.data(mh + 156);
    const auto *mh_157 = buffer.data(mh + 157);
    const auto *mh_158 = buffer.data(mh + 158);
    const auto *mh_159 = buffer.data(mh + 159);
    const auto *mh_160 = buffer.data(mh + 160);
    const auto *mh_161 = buffer.data(mh + 161);
    const auto *mh_162 = buffer.data(mh + 162);
    const auto *mh_163 = buffer.data(mh + 163);
    const auto *mh_164 = buffer.data(mh + 164);
    const auto *mh_165 = buffer.data(mh + 165);
    const auto *mh_166 = buffer.data(mh + 166);
    const auto *mh_167 = buffer.data(mh + 167);
    const auto *mh_168 = buffer.data(mh + 168);
    const auto *mh_169 = buffer.data(mh + 169);
    const auto *mh_170 = buffer.data(mh + 170);
    const auto *mh_171 = buffer.data(mh + 171);
    const auto *mh_172 = buffer.data(mh + 172);
    const auto *mh_173 = buffer.data(mh + 173);
    const auto *mh_174 = buffer.data(mh + 174);
    const auto *mh_175 = buffer.data(mh + 175);
    const auto *mh_176 = buffer.data(mh + 176);
    const auto *mh_177 = buffer.data(mh + 177);
    const auto *mh_178 = buffer.data(mh + 178);
    const auto *mh_179 = buffer.data(mh + 179);
    const auto *mh_180 = buffer.data(mh + 180);
    const auto *mh_181 = buffer.data(mh + 181);
    const auto *mh_182 = buffer.data(mh + 182);
    const auto *mh_183 = buffer.data(mh + 183);
    const auto *mh_184 = buffer.data(mh + 184);
    const auto *mh_185 = buffer.data(mh + 185);
    const auto *mh_186 = buffer.data(mh + 186);
    const auto *mh_187 = buffer.data(mh + 187);
    const auto *mh_188 = buffer.data(mh + 188);
    const auto *mh_189 = buffer.data(mh + 189);
    const auto *mh_190 = buffer.data(mh + 190);
    const auto *mh_191 = buffer.data(mh + 191);
    const auto *mh_192 = buffer.data(mh + 192);
    const auto *mh_193 = buffer.data(mh + 193);
    const auto *mh_194 = buffer.data(mh + 194);
    const auto *mh_195 = buffer.data(mh + 195);
    const auto *mh_196 = buffer.data(mh + 196);
    const auto *mh_197 = buffer.data(mh + 197);
    const auto *mh_198 = buffer.data(mh + 198);
    const auto *mh_199 = buffer.data(mh + 199);
    const auto *mh_200 = buffer.data(mh + 200);
    const auto *mh_201 = buffer.data(mh + 201);
    const auto *mh_202 = buffer.data(mh + 202);
    const auto *mh_203 = buffer.data(mh + 203);
    const auto *mh_204 = buffer.data(mh + 204);
    const auto *mh_205 = buffer.data(mh + 205);
    const auto *mh_206 = buffer.data(mh + 206);
    const auto *mh_207 = buffer.data(mh + 207);
    const auto *mh_208 = buffer.data(mh + 208);
    const auto *mh_209 = buffer.data(mh + 209);
    const auto *mh_231 = buffer.data(mh + 231);
    const auto *mh_232 = buffer.data(mh + 232);
    const auto *mh_233 = buffer.data(mh + 233);
    const auto *mh_234 = buffer.data(mh + 234);
    const auto *mh_235 = buffer.data(mh + 235);
    const auto *mh_236 = buffer.data(mh + 236);
    const auto *mh_237 = buffer.data(mh + 237);
    const auto *mh_238 = buffer.data(mh + 238);
    const auto *mh_239 = buffer.data(mh + 239);
    const auto *mh_240 = buffer.data(mh + 240);
    const auto *mh_241 = buffer.data(mh + 241);
    const auto *mh_242 = buffer.data(mh + 242);
    const auto *mh_243 = buffer.data(mh + 243);
    const auto *mh_244 = buffer.data(mh + 244);
    const auto *mh_245 = buffer.data(mh + 245);
    const auto *mh_246 = buffer.data(mh + 246);
    const auto *mh_247 = buffer.data(mh + 247);
    const auto *mh_248 = buffer.data(mh + 248);
    const auto *mh_249 = buffer.data(mh + 249);
    const auto *mh_250 = buffer.data(mh + 250);
    const auto *mh_251 = buffer.data(mh + 251);
    const auto *mh_252 = buffer.data(mh + 252);
    const auto *mh_253 = buffer.data(mh + 253);
    const auto *mh_254 = buffer.data(mh + 254);
    const auto *mh_255 = buffer.data(mh + 255);
    const auto *mh_256 = buffer.data(mh + 256);
    const auto *mh_257 = buffer.data(mh + 257);
    const auto *mh_258 = buffer.data(mh + 258);
    const auto *mh_259 = buffer.data(mh + 259);
    const auto *mh_260 = buffer.data(mh + 260);
    const auto *mh_261 = buffer.data(mh + 261);
    const auto *mh_262 = buffer.data(mh + 262);
    const auto *mh_263 = buffer.data(mh + 263);
    const auto *mh_264 = buffer.data(mh + 264);
    const auto *mh_265 = buffer.data(mh + 265);
    const auto *mh_266 = buffer.data(mh + 266);
    const auto *mh_267 = buffer.data(mh + 267);
    const auto *mh_268 = buffer.data(mh + 268);
    const auto *mh_269 = buffer.data(mh + 269);
    const auto *mh_270 = buffer.data(mh + 270);
    const auto *mh_271 = buffer.data(mh + 271);
    const auto *mh_272 = buffer.data(mh + 272);
    const auto *mh_273 = buffer.data(mh + 273);
    const auto *mh_274 = buffer.data(mh + 274);
    const auto *mh_275 = buffer.data(mh + 275);
    const auto *mh_276 = buffer.data(mh + 276);
    const auto *mh_277 = buffer.data(mh + 277);
    const auto *mh_278 = buffer.data(mh + 278);
    const auto *mh_279 = buffer.data(mh + 279);
    const auto *mh_280 = buffer.data(mh + 280);
    const auto *mh_281 = buffer.data(mh + 281);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, mh_42, mh_43, mh_44, mh_45, \
                         mh_46, mh_47, mh_48, mh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * mh_42[k];

        t_1[k] = f_0 * mh_43[k];

        t_2[k] = f_0 * mh_44[k];

        t_3[k] = f_0 * mh_45[k];

        t_4[k] = f_0 * mh_46[k];

        t_5[k] = f_0 * mh_47[k];

        t_6[k] = f_0 * mh_48[k];

        t_7[k] = f_0 * mh_49[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, mh_50, mh_51, mh_52, \
                         mh_53, mh_54, mh_55, mh_56, mh_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * mh_50[k];

        t_9[k] = f_0 * mh_51[k];

        t_10[k] = f_0 * mh_52[k];

        t_11[k] = f_0 * mh_53[k];

        t_12[k] = f_0 * mh_54[k];

        t_13[k] = f_0 * mh_55[k];

        t_14[k] = f_0 * mh_56[k];

        t_15[k] = f_0 * mh_57[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, mh_58, mh_59, mh_60, \
                         mh_61, mh_62, mh_84, mh_85, mh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * mh_58[k];

        t_17[k] = f_0 * mh_59[k];

        t_18[k] = f_0 * mh_60[k];

        t_19[k] = f_0 * mh_61[k];

        t_20[k] = f_0 * mh_62[k];

        t_21[k] = f_0 * mh_84[k];

        t_22[k] = f_0 * mh_85[k];

        t_23[k] = f_0 * mh_86[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, mh_87, mh_88, mh_89, \
                         mh_90, mh_91, mh_92, mh_93, mh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * mh_87[k];

        t_25[k] = f_0 * mh_88[k];

        t_26[k] = f_0 * mh_89[k];

        t_27[k] = f_0 * mh_90[k];

        t_28[k] = f_0 * mh_91[k];

        t_29[k] = f_0 * mh_92[k];

        t_30[k] = f_0 * mh_93[k];

        t_31[k] = f_0 * mh_94[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, mh_95, mh_96, mh_97, \
                         mh_98, mh_99, mh_100, mh_101, mh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * mh_95[k];

        t_33[k] = f_0 * mh_96[k];

        t_34[k] = f_0 * mh_97[k];

        t_35[k] = f_0 * mh_98[k];

        t_36[k] = f_0 * mh_99[k];

        t_37[k] = f_0 * mh_100[k];

        t_38[k] = f_0 * mh_101[k];

        t_39[k] = f_0 * mh_102[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, kh_0, kh_1, kh_2, kh_3, mh_103, \
                         mh_104, mh_105, mh_106, mh_107, mh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * mh_103[k];

        t_41[k] = f_0 * mh_104[k];

        t_42[k] = -kh_0[k]
                  + f_0 * mh_105[k];

        t_43[k] = -kh_1[k]
                  + f_0 * mh_106[k];

        t_44[k] = -kh_2[k]
                  + f_0 * mh_107[k];

        t_45[k] = -kh_3[k]
                  + f_0 * mh_108[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, kh_4, kh_5, kh_6, kh_7, kh_8, mh_109, \
                         mh_110, mh_111, mh_112, mh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -kh_4[k]
                  + f_0 * mh_109[k];

        t_47[k] = -kh_5[k]
                  + f_0 * mh_110[k];

        t_48[k] = -kh_6[k]
                  + f_0 * mh_111[k];

        t_49[k] = -kh_7[k]
                  + f_0 * mh_112[k];

        t_50[k] = -kh_8[k]
                  + f_0 * mh_113[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, kh_9, kh_10, kh_11, kh_12, kh_13, \
                         mh_114, mh_115, mh_116, mh_117, mh_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -kh_9[k]
                  + f_0 * mh_114[k];

        t_52[k] = -kh_10[k]
                  + f_0 * mh_115[k];

        t_53[k] = -kh_11[k]
                  + f_0 * mh_116[k];

        t_54[k] = -kh_12[k]
                  + f_0 * mh_117[k];

        t_55[k] = -kh_13[k]
                  + f_0 * mh_118[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, kh_14, kh_15, kh_16, kh_17, kh_18, \
                         mh_119, mh_120, mh_121, mh_122, mh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -kh_14[k]
                  + f_0 * mh_119[k];

        t_57[k] = -kh_15[k]
                  + f_0 * mh_120[k];

        t_58[k] = -kh_16[k]
                  + f_0 * mh_121[k];

        t_59[k] = -kh_17[k]
                  + f_0 * mh_122[k];

        t_60[k] = -kh_18[k]
                  + f_0 * mh_123[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, t_67, kh_19, kh_20, mh_124, \
                         mh_125, mh_147, mh_148, mh_149, mh_150, \
                         mh_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -kh_19[k]
                  + f_0 * mh_124[k];

        t_62[k] = -kh_20[k]
                  + f_0 * mh_125[k];

        t_63[k] = f_0 * mh_147[k];

        t_64[k] = f_0 * mh_148[k];

        t_65[k] = f_0 * mh_149[k];

        t_66[k] = f_0 * mh_150[k];

        t_67[k] = f_0 * mh_151[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, t_74, t_75, mh_152, mh_153, \
                         mh_154, mh_155, mh_156, mh_157, mh_158, \
                         mh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * mh_152[k];

        t_69[k] = f_0 * mh_153[k];

        t_70[k] = f_0 * mh_154[k];

        t_71[k] = f_0 * mh_155[k];

        t_72[k] = f_0 * mh_156[k];

        t_73[k] = f_0 * mh_157[k];

        t_74[k] = f_0 * mh_158[k];

        t_75[k] = f_0 * mh_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, t_82, t_83, mh_160, mh_161, \
                         mh_162, mh_163, mh_164, mh_165, mh_166, \
                         mh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_0 * mh_160[k];

        t_77[k] = f_0 * mh_161[k];

        t_78[k] = f_0 * mh_162[k];

        t_79[k] = f_0 * mh_163[k];

        t_80[k] = f_0 * mh_164[k];

        t_81[k] = f_0 * mh_165[k];

        t_82[k] = f_0 * mh_166[k];

        t_83[k] = f_0 * mh_167[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, kh_21, kh_22, kh_23, kh_24, kh_25, \
                         mh_168, mh_169, mh_170, mh_171, mh_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -kh_21[k]
                  + f_0 * mh_168[k];

        t_85[k] = -kh_22[k]
                  + f_0 * mh_169[k];

        t_86[k] = -kh_23[k]
                  + f_0 * mh_170[k];

        t_87[k] = -kh_24[k]
                  + f_0 * mh_171[k];

        t_88[k] = -kh_25[k]
                  + f_0 * mh_172[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, kh_26, kh_27, kh_28, kh_29, kh_30, \
                         mh_173, mh_174, mh_175, mh_176, mh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -kh_26[k]
                  + f_0 * mh_173[k];

        t_90[k] = -kh_27[k]
                  + f_0 * mh_174[k];

        t_91[k] = -kh_28[k]
                  + f_0 * mh_175[k];

        t_92[k] = -kh_29[k]
                  + f_0 * mh_176[k];

        t_93[k] = -kh_30[k]
                  + f_0 * mh_177[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, kh_31, kh_32, kh_33, kh_34, kh_35, \
                         mh_178, mh_179, mh_180, mh_181, mh_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -kh_31[k]
                  + f_0 * mh_178[k];

        t_95[k] = -kh_32[k]
                  + f_0 * mh_179[k];

        t_96[k] = -kh_33[k]
                  + f_0 * mh_180[k];

        t_97[k] = -kh_34[k]
                  + f_0 * mh_181[k];

        t_98[k] = -kh_35[k]
                  + f_0 * mh_182[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, kh_36, kh_37, kh_38, kh_39, kh_40, \
                         mh_183, mh_184, mh_185, mh_186, mh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -kh_36[k]
                  + f_0 * mh_183[k];

        t_100[k] = -kh_37[k]
                   + f_0 * mh_184[k];

        t_101[k] = -kh_38[k]
                   + f_0 * mh_185[k];

        t_102[k] = -kh_39[k]
                   + f_0 * mh_186[k];

        t_103[k] = -kh_40[k]
                   + f_0 * mh_187[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, kh_41, kh_42, kh_43, kh_44, kh_45, \
                         mh_188, mh_189, mh_190, mh_191, mh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -kh_41[k]
                   + f_0 * mh_188[k];

        t_105[k] = -2.0 * kh_42[k]
                   + f_0 * mh_189[k];

        t_106[k] = -2.0 * kh_43[k]
                   + f_0 * mh_190[k];

        t_107[k] = -2.0 * kh_44[k]
                   + f_0 * mh_191[k];

        t_108[k] = -2.0 * kh_45[k]
                   + f_0 * mh_192[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, kh_46, kh_47, kh_48, kh_49, kh_50, \
                         mh_193, mh_194, mh_195, mh_196, mh_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -2.0 * kh_46[k]
                   + f_0 * mh_193[k];

        t_110[k] = -2.0 * kh_47[k]
                   + f_0 * mh_194[k];

        t_111[k] = -2.0 * kh_48[k]
                   + f_0 * mh_195[k];

        t_112[k] = -2.0 * kh_49[k]
                   + f_0 * mh_196[k];

        t_113[k] = -2.0 * kh_50[k]
                   + f_0 * mh_197[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, kh_51, kh_52, kh_53, kh_54, kh_55, \
                         mh_198, mh_199, mh_200, mh_201, mh_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -2.0 * kh_51[k]
                   + f_0 * mh_198[k];

        t_115[k] = -2.0 * kh_52[k]
                   + f_0 * mh_199[k];

        t_116[k] = -2.0 * kh_53[k]
                   + f_0 * mh_200[k];

        t_117[k] = -2.0 * kh_54[k]
                   + f_0 * mh_201[k];

        t_118[k] = -2.0 * kh_55[k]
                   + f_0 * mh_202[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, kh_56, kh_57, kh_58, kh_59, kh_60, \
                         mh_203, mh_204, mh_205, mh_206, mh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -2.0 * kh_56[k]
                   + f_0 * mh_203[k];

        t_120[k] = -2.0 * kh_57[k]
                   + f_0 * mh_204[k];

        t_121[k] = -2.0 * kh_58[k]
                   + f_0 * mh_205[k];

        t_122[k] = -2.0 * kh_59[k]
                   + f_0 * mh_206[k];

        t_123[k] = -2.0 * kh_60[k]
                   + f_0 * mh_207[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, t_129, t_130, kh_61, kh_62, \
                         mh_208, mh_209, mh_231, mh_232, mh_233, mh_234, \
                         mh_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * kh_61[k]
                   + f_0 * mh_208[k];

        t_125[k] = -2.0 * kh_62[k]
                   + f_0 * mh_209[k];

        t_126[k] = f_0 * mh_231[k];

        t_127[k] = f_0 * mh_232[k];

        t_128[k] = f_0 * mh_233[k];

        t_129[k] = f_0 * mh_234[k];

        t_130[k] = f_0 * mh_235[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, t_136, t_137, t_138, mh_236, \
                         mh_237, mh_238, mh_239, mh_240, mh_241, mh_242, \
                         mh_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_0 * mh_236[k];

        t_132[k] = f_0 * mh_237[k];

        t_133[k] = f_0 * mh_238[k];

        t_134[k] = f_0 * mh_239[k];

        t_135[k] = f_0 * mh_240[k];

        t_136[k] = f_0 * mh_241[k];

        t_137[k] = f_0 * mh_242[k];

        t_138[k] = f_0 * mh_243[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, t_145, t_146, mh_244, \
                         mh_245, mh_246, mh_247, mh_248, mh_249, mh_250, \
                         mh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_0 * mh_244[k];

        t_140[k] = f_0 * mh_245[k];

        t_141[k] = f_0 * mh_246[k];

        t_142[k] = f_0 * mh_247[k];

        t_143[k] = f_0 * mh_248[k];

        t_144[k] = f_0 * mh_249[k];

        t_145[k] = f_0 * mh_250[k];

        t_146[k] = f_0 * mh_251[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, kh_63, kh_64, kh_65, kh_66, kh_67, \
                         mh_252, mh_253, mh_254, mh_255, mh_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -kh_63[k]
                   + f_0 * mh_252[k];

        t_148[k] = -kh_64[k]
                   + f_0 * mh_253[k];

        t_149[k] = -kh_65[k]
                   + f_0 * mh_254[k];

        t_150[k] = -kh_66[k]
                   + f_0 * mh_255[k];

        t_151[k] = -kh_67[k]
                   + f_0 * mh_256[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, kh_68, kh_69, kh_70, kh_71, kh_72, \
                         mh_257, mh_258, mh_259, mh_260, mh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -kh_68[k]
                   + f_0 * mh_257[k];

        t_153[k] = -kh_69[k]
                   + f_0 * mh_258[k];

        t_154[k] = -kh_70[k]
                   + f_0 * mh_259[k];

        t_155[k] = -kh_71[k]
                   + f_0 * mh_260[k];

        t_156[k] = -kh_72[k]
                   + f_0 * mh_261[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, kh_73, kh_74, kh_75, kh_76, kh_77, \
                         mh_262, mh_263, mh_264, mh_265, mh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -kh_73[k]
                   + f_0 * mh_262[k];

        t_158[k] = -kh_74[k]
                   + f_0 * mh_263[k];

        t_159[k] = -kh_75[k]
                   + f_0 * mh_264[k];

        t_160[k] = -kh_76[k]
                   + f_0 * mh_265[k];

        t_161[k] = -kh_77[k]
                   + f_0 * mh_266[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, kh_78, kh_79, kh_80, kh_81, kh_82, \
                         mh_267, mh_268, mh_269, mh_270, mh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -kh_78[k]
                   + f_0 * mh_267[k];

        t_163[k] = -kh_79[k]
                   + f_0 * mh_268[k];

        t_164[k] = -kh_80[k]
                   + f_0 * mh_269[k];

        t_165[k] = -kh_81[k]
                   + f_0 * mh_270[k];

        t_166[k] = -kh_82[k]
                   + f_0 * mh_271[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, kh_83, kh_84, kh_85, kh_86, kh_87, \
                         mh_272, mh_273, mh_274, mh_275, mh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -kh_83[k]
                   + f_0 * mh_272[k];

        t_168[k] = -2.0 * kh_84[k]
                   + f_0 * mh_273[k];

        t_169[k] = -2.0 * kh_85[k]
                   + f_0 * mh_274[k];

        t_170[k] = -2.0 * kh_86[k]
                   + f_0 * mh_275[k];

        t_171[k] = -2.0 * kh_87[k]
                   + f_0 * mh_276[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, kh_88, kh_89, kh_90, kh_91, kh_92, \
                         mh_277, mh_278, mh_279, mh_280, mh_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -2.0 * kh_88[k]
                   + f_0 * mh_277[k];

        t_173[k] = -2.0 * kh_89[k]
                   + f_0 * mh_278[k];

        t_174[k] = -2.0 * kh_90[k]
                   + f_0 * mh_279[k];

        t_175[k] = -2.0 * kh_91[k]
                   + f_0 * mh_280[k];

        t_176[k] = -2.0 * kh_92[k]
                   + f_0 * mh_281[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_177 = buffer.data(target + 177);
    auto *t_178 = buffer.data(target + 178);
    auto *t_179 = buffer.data(target + 179);
    auto *t_180 = buffer.data(target + 180);
    auto *t_181 = buffer.data(target + 181);
    auto *t_182 = buffer.data(target + 182);
    auto *t_183 = buffer.data(target + 183);
    auto *t_184 = buffer.data(target + 184);
    auto *t_185 = buffer.data(target + 185);
    auto *t_186 = buffer.data(target + 186);
    auto *t_187 = buffer.data(target + 187);
    auto *t_188 = buffer.data(target + 188);
    auto *t_189 = buffer.data(target + 189);
    auto *t_190 = buffer.data(target + 190);
    auto *t_191 = buffer.data(target + 191);
    auto *t_192 = buffer.data(target + 192);
    auto *t_193 = buffer.data(target + 193);
    auto *t_194 = buffer.data(target + 194);
    auto *t_195 = buffer.data(target + 195);
    auto *t_196 = buffer.data(target + 196);
    auto *t_197 = buffer.data(target + 197);
    auto *t_198 = buffer.data(target + 198);
    auto *t_199 = buffer.data(target + 199);
    auto *t_200 = buffer.data(target + 200);
    auto *t_201 = buffer.data(target + 201);
    auto *t_202 = buffer.data(target + 202);
    auto *t_203 = buffer.data(target + 203);
    auto *t_204 = buffer.data(target + 204);
    auto *t_205 = buffer.data(target + 205);
    auto *t_206 = buffer.data(target + 206);
    auto *t_207 = buffer.data(target + 207);
    auto *t_208 = buffer.data(target + 208);
    auto *t_209 = buffer.data(target + 209);
    auto *t_210 = buffer.data(target + 210);
    auto *t_211 = buffer.data(target + 211);
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);
    auto *t_215 = buffer.data(target + 215);
    auto *t_216 = buffer.data(target + 216);
    auto *t_217 = buffer.data(target + 217);
    auto *t_218 = buffer.data(target + 218);
    auto *t_219 = buffer.data(target + 219);
    auto *t_220 = buffer.data(target + 220);
    auto *t_221 = buffer.data(target + 221);
    auto *t_222 = buffer.data(target + 222);
    auto *t_223 = buffer.data(target + 223);
    auto *t_224 = buffer.data(target + 224);
    auto *t_225 = buffer.data(target + 225);
    auto *t_226 = buffer.data(target + 226);
    auto *t_227 = buffer.data(target + 227);
    auto *t_228 = buffer.data(target + 228);
    auto *t_229 = buffer.data(target + 229);
    auto *t_230 = buffer.data(target + 230);
    auto *t_231 = buffer.data(target + 231);
    auto *t_232 = buffer.data(target + 232);
    auto *t_233 = buffer.data(target + 233);
    auto *t_234 = buffer.data(target + 234);
    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);
    auto *t_239 = buffer.data(target + 239);
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);
    auto *t_262 = buffer.data(target + 262);
    auto *t_263 = buffer.data(target + 263);
    auto *t_264 = buffer.data(target + 264);
    auto *t_265 = buffer.data(target + 265);
    auto *t_266 = buffer.data(target + 266);
    auto *t_267 = buffer.data(target + 267);
    auto *t_268 = buffer.data(target + 268);
    auto *t_269 = buffer.data(target + 269);
    auto *t_270 = buffer.data(target + 270);
    auto *t_271 = buffer.data(target + 271);
    auto *t_272 = buffer.data(target + 272);
    auto *t_273 = buffer.data(target + 273);
    auto *t_274 = buffer.data(target + 274);
    auto *t_275 = buffer.data(target + 275);
    auto *t_276 = buffer.data(target + 276);
    auto *t_277 = buffer.data(target + 277);
    auto *t_278 = buffer.data(target + 278);
    auto *t_279 = buffer.data(target + 279);
    auto *t_280 = buffer.data(target + 280);
    auto *t_281 = buffer.data(target + 281);
    auto *t_282 = buffer.data(target + 282);
    auto *t_283 = buffer.data(target + 283);
    auto *t_284 = buffer.data(target + 284);
    auto *t_285 = buffer.data(target + 285);
    auto *t_286 = buffer.data(target + 286);
    auto *t_287 = buffer.data(target + 287);
    auto *t_288 = buffer.data(target + 288);
    auto *t_289 = buffer.data(target + 289);
    auto *t_290 = buffer.data(target + 290);
    auto *t_291 = buffer.data(target + 291);
    auto *t_292 = buffer.data(target + 292);
    auto *t_293 = buffer.data(target + 293);
    auto *t_294 = buffer.data(target + 294);
    auto *t_295 = buffer.data(target + 295);
    auto *t_296 = buffer.data(target + 296);
    auto *t_297 = buffer.data(target + 297);
    auto *t_298 = buffer.data(target + 298);
    auto *t_299 = buffer.data(target + 299);
    auto *t_300 = buffer.data(target + 300);
    auto *t_301 = buffer.data(target + 301);
    auto *t_302 = buffer.data(target + 302);
    auto *t_303 = buffer.data(target + 303);
    auto *t_304 = buffer.data(target + 304);
    auto *t_305 = buffer.data(target + 305);
    auto *t_306 = buffer.data(target + 306);
    auto *t_307 = buffer.data(target + 307);
    auto *t_308 = buffer.data(target + 308);
    auto *t_309 = buffer.data(target + 309);
    auto *t_310 = buffer.data(target + 310);
    auto *t_311 = buffer.data(target + 311);
    auto *t_312 = buffer.data(target + 312);
    auto *t_313 = buffer.data(target + 313);
    auto *t_314 = buffer.data(target + 314);
    auto *t_315 = buffer.data(target + 315);
    auto *t_316 = buffer.data(target + 316);
    auto *t_317 = buffer.data(target + 317);
    auto *t_318 = buffer.data(target + 318);
    auto *t_319 = buffer.data(target + 319);
    auto *t_320 = buffer.data(target + 320);
    auto *t_321 = buffer.data(target + 321);
    auto *t_322 = buffer.data(target + 322);
    auto *t_323 = buffer.data(target + 323);
    auto *t_324 = buffer.data(target + 324);
    auto *t_325 = buffer.data(target + 325);
    auto *t_326 = buffer.data(target + 326);
    auto *t_327 = buffer.data(target + 327);
    auto *t_328 = buffer.data(target + 328);
    auto *t_329 = buffer.data(target + 329);
    auto *t_330 = buffer.data(target + 330);
    auto *t_331 = buffer.data(target + 331);
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);
    auto *t_334 = buffer.data(target + 334);
    auto *t_335 = buffer.data(target + 335);
    auto *t_336 = buffer.data(target + 336);
    auto *t_337 = buffer.data(target + 337);

    const auto *kh_93 = buffer.data(kh + 93);
    const auto *kh_94 = buffer.data(kh + 94);
    const auto *kh_95 = buffer.data(kh + 95);
    const auto *kh_96 = buffer.data(kh + 96);
    const auto *kh_97 = buffer.data(kh + 97);
    const auto *kh_98 = buffer.data(kh + 98);
    const auto *kh_99 = buffer.data(kh + 99);
    const auto *kh_100 = buffer.data(kh + 100);
    const auto *kh_101 = buffer.data(kh + 101);
    const auto *kh_102 = buffer.data(kh + 102);
    const auto *kh_103 = buffer.data(kh + 103);
    const auto *kh_104 = buffer.data(kh + 104);
    const auto *kh_105 = buffer.data(kh + 105);
    const auto *kh_106 = buffer.data(kh + 106);
    const auto *kh_107 = buffer.data(kh + 107);
    const auto *kh_108 = buffer.data(kh + 108);
    const auto *kh_109 = buffer.data(kh + 109);
    const auto *kh_110 = buffer.data(kh + 110);
    const auto *kh_111 = buffer.data(kh + 111);
    const auto *kh_112 = buffer.data(kh + 112);
    const auto *kh_113 = buffer.data(kh + 113);
    const auto *kh_114 = buffer.data(kh + 114);
    const auto *kh_115 = buffer.data(kh + 115);
    const auto *kh_116 = buffer.data(kh + 116);
    const auto *kh_117 = buffer.data(kh + 117);
    const auto *kh_118 = buffer.data(kh + 118);
    const auto *kh_119 = buffer.data(kh + 119);
    const auto *kh_120 = buffer.data(kh + 120);
    const auto *kh_121 = buffer.data(kh + 121);
    const auto *kh_122 = buffer.data(kh + 122);
    const auto *kh_123 = buffer.data(kh + 123);
    const auto *kh_124 = buffer.data(kh + 124);
    const auto *kh_125 = buffer.data(kh + 125);
    const auto *kh_126 = buffer.data(kh + 126);
    const auto *kh_127 = buffer.data(kh + 127);
    const auto *kh_128 = buffer.data(kh + 128);
    const auto *kh_129 = buffer.data(kh + 129);
    const auto *kh_130 = buffer.data(kh + 130);
    const auto *kh_131 = buffer.data(kh + 131);
    const auto *kh_132 = buffer.data(kh + 132);
    const auto *kh_133 = buffer.data(kh + 133);
    const auto *kh_134 = buffer.data(kh + 134);
    const auto *kh_135 = buffer.data(kh + 135);
    const auto *kh_136 = buffer.data(kh + 136);
    const auto *kh_137 = buffer.data(kh + 137);
    const auto *kh_138 = buffer.data(kh + 138);
    const auto *kh_139 = buffer.data(kh + 139);
    const auto *kh_140 = buffer.data(kh + 140);
    const auto *kh_141 = buffer.data(kh + 141);
    const auto *kh_142 = buffer.data(kh + 142);
    const auto *kh_143 = buffer.data(kh + 143);
    const auto *kh_144 = buffer.data(kh + 144);
    const auto *kh_145 = buffer.data(kh + 145);
    const auto *kh_146 = buffer.data(kh + 146);
    const auto *kh_147 = buffer.data(kh + 147);
    const auto *kh_148 = buffer.data(kh + 148);
    const auto *kh_149 = buffer.data(kh + 149);
    const auto *kh_150 = buffer.data(kh + 150);
    const auto *kh_151 = buffer.data(kh + 151);
    const auto *kh_152 = buffer.data(kh + 152);
    const auto *kh_153 = buffer.data(kh + 153);
    const auto *kh_154 = buffer.data(kh + 154);
    const auto *kh_155 = buffer.data(kh + 155);
    const auto *kh_156 = buffer.data(kh + 156);
    const auto *kh_157 = buffer.data(kh + 157);
    const auto *kh_158 = buffer.data(kh + 158);
    const auto *kh_159 = buffer.data(kh + 159);
    const auto *kh_160 = buffer.data(kh + 160);
    const auto *kh_161 = buffer.data(kh + 161);
    const auto *kh_162 = buffer.data(kh + 162);
    const auto *kh_163 = buffer.data(kh + 163);
    const auto *kh_164 = buffer.data(kh + 164);
    const auto *kh_165 = buffer.data(kh + 165);
    const auto *kh_166 = buffer.data(kh + 166);
    const auto *kh_167 = buffer.data(kh + 167);
    const auto *kh_168 = buffer.data(kh + 168);
    const auto *kh_169 = buffer.data(kh + 169);
    const auto *kh_170 = buffer.data(kh + 170);
    const auto *kh_171 = buffer.data(kh + 171);
    const auto *kh_172 = buffer.data(kh + 172);
    const auto *kh_173 = buffer.data(kh + 173);
    const auto *kh_174 = buffer.data(kh + 174);
    const auto *kh_175 = buffer.data(kh + 175);
    const auto *kh_176 = buffer.data(kh + 176);
    const auto *kh_177 = buffer.data(kh + 177);
    const auto *kh_178 = buffer.data(kh + 178);
    const auto *kh_179 = buffer.data(kh + 179);
    const auto *kh_180 = buffer.data(kh + 180);
    const auto *kh_181 = buffer.data(kh + 181);
    const auto *kh_182 = buffer.data(kh + 182);
    const auto *kh_183 = buffer.data(kh + 183);
    const auto *kh_184 = buffer.data(kh + 184);
    const auto *kh_185 = buffer.data(kh + 185);
    const auto *kh_186 = buffer.data(kh + 186);
    const auto *kh_187 = buffer.data(kh + 187);
    const auto *kh_188 = buffer.data(kh + 188);
    const auto *kh_189 = buffer.data(kh + 189);
    const auto *kh_190 = buffer.data(kh + 190);
    const auto *kh_191 = buffer.data(kh + 191);
    const auto *kh_192 = buffer.data(kh + 192);
    const auto *kh_193 = buffer.data(kh + 193);
    const auto *kh_194 = buffer.data(kh + 194);
    const auto *kh_195 = buffer.data(kh + 195);
    const auto *kh_196 = buffer.data(kh + 196);
    const auto *kh_197 = buffer.data(kh + 197);
    const auto *kh_198 = buffer.data(kh + 198);
    const auto *kh_199 = buffer.data(kh + 199);
    const auto *kh_200 = buffer.data(kh + 200);
    const auto *kh_201 = buffer.data(kh + 201);
    const auto *kh_202 = buffer.data(kh + 202);
    const auto *kh_203 = buffer.data(kh + 203);
    const auto *kh_204 = buffer.data(kh + 204);
    const auto *kh_205 = buffer.data(kh + 205);
    const auto *kh_206 = buffer.data(kh + 206);
    const auto *kh_207 = buffer.data(kh + 207);
    const auto *kh_208 = buffer.data(kh + 208);
    const auto *kh_209 = buffer.data(kh + 209);
    const auto *kh_210 = buffer.data(kh + 210);
    const auto *kh_211 = buffer.data(kh + 211);

    const auto *mh_282 = buffer.data(mh + 282);
    const auto *mh_283 = buffer.data(mh + 283);
    const auto *mh_284 = buffer.data(mh + 284);
    const auto *mh_285 = buffer.data(mh + 285);
    const auto *mh_286 = buffer.data(mh + 286);
    const auto *mh_287 = buffer.data(mh + 287);
    const auto *mh_288 = buffer.data(mh + 288);
    const auto *mh_289 = buffer.data(mh + 289);
    const auto *mh_290 = buffer.data(mh + 290);
    const auto *mh_291 = buffer.data(mh + 291);
    const auto *mh_292 = buffer.data(mh + 292);
    const auto *mh_293 = buffer.data(mh + 293);
    const auto *mh_294 = buffer.data(mh + 294);
    const auto *mh_295 = buffer.data(mh + 295);
    const auto *mh_296 = buffer.data(mh + 296);
    const auto *mh_297 = buffer.data(mh + 297);
    const auto *mh_298 = buffer.data(mh + 298);
    const auto *mh_299 = buffer.data(mh + 299);
    const auto *mh_300 = buffer.data(mh + 300);
    const auto *mh_301 = buffer.data(mh + 301);
    const auto *mh_302 = buffer.data(mh + 302);
    const auto *mh_303 = buffer.data(mh + 303);
    const auto *mh_304 = buffer.data(mh + 304);
    const auto *mh_305 = buffer.data(mh + 305);
    const auto *mh_306 = buffer.data(mh + 306);
    const auto *mh_307 = buffer.data(mh + 307);
    const auto *mh_308 = buffer.data(mh + 308);
    const auto *mh_309 = buffer.data(mh + 309);
    const auto *mh_310 = buffer.data(mh + 310);
    const auto *mh_311 = buffer.data(mh + 311);
    const auto *mh_312 = buffer.data(mh + 312);
    const auto *mh_313 = buffer.data(mh + 313);
    const auto *mh_314 = buffer.data(mh + 314);
    const auto *mh_336 = buffer.data(mh + 336);
    const auto *mh_337 = buffer.data(mh + 337);
    const auto *mh_338 = buffer.data(mh + 338);
    const auto *mh_339 = buffer.data(mh + 339);
    const auto *mh_340 = buffer.data(mh + 340);
    const auto *mh_341 = buffer.data(mh + 341);
    const auto *mh_342 = buffer.data(mh + 342);
    const auto *mh_343 = buffer.data(mh + 343);
    const auto *mh_344 = buffer.data(mh + 344);
    const auto *mh_345 = buffer.data(mh + 345);
    const auto *mh_346 = buffer.data(mh + 346);
    const auto *mh_347 = buffer.data(mh + 347);
    const auto *mh_348 = buffer.data(mh + 348);
    const auto *mh_349 = buffer.data(mh + 349);
    const auto *mh_350 = buffer.data(mh + 350);
    const auto *mh_351 = buffer.data(mh + 351);
    const auto *mh_352 = buffer.data(mh + 352);
    const auto *mh_353 = buffer.data(mh + 353);
    const auto *mh_354 = buffer.data(mh + 354);
    const auto *mh_355 = buffer.data(mh + 355);
    const auto *mh_356 = buffer.data(mh + 356);
    const auto *mh_357 = buffer.data(mh + 357);
    const auto *mh_358 = buffer.data(mh + 358);
    const auto *mh_359 = buffer.data(mh + 359);
    const auto *mh_360 = buffer.data(mh + 360);
    const auto *mh_361 = buffer.data(mh + 361);
    const auto *mh_362 = buffer.data(mh + 362);
    const auto *mh_363 = buffer.data(mh + 363);
    const auto *mh_364 = buffer.data(mh + 364);
    const auto *mh_365 = buffer.data(mh + 365);
    const auto *mh_366 = buffer.data(mh + 366);
    const auto *mh_367 = buffer.data(mh + 367);
    const auto *mh_368 = buffer.data(mh + 368);
    const auto *mh_369 = buffer.data(mh + 369);
    const auto *mh_370 = buffer.data(mh + 370);
    const auto *mh_371 = buffer.data(mh + 371);
    const auto *mh_372 = buffer.data(mh + 372);
    const auto *mh_373 = buffer.data(mh + 373);
    const auto *mh_374 = buffer.data(mh + 374);
    const auto *mh_375 = buffer.data(mh + 375);
    const auto *mh_376 = buffer.data(mh + 376);
    const auto *mh_377 = buffer.data(mh + 377);
    const auto *mh_378 = buffer.data(mh + 378);
    const auto *mh_379 = buffer.data(mh + 379);
    const auto *mh_380 = buffer.data(mh + 380);
    const auto *mh_381 = buffer.data(mh + 381);
    const auto *mh_382 = buffer.data(mh + 382);
    const auto *mh_383 = buffer.data(mh + 383);
    const auto *mh_384 = buffer.data(mh + 384);
    const auto *mh_385 = buffer.data(mh + 385);
    const auto *mh_386 = buffer.data(mh + 386);
    const auto *mh_387 = buffer.data(mh + 387);
    const auto *mh_388 = buffer.data(mh + 388);
    const auto *mh_389 = buffer.data(mh + 389);
    const auto *mh_390 = buffer.data(mh + 390);
    const auto *mh_391 = buffer.data(mh + 391);
    const auto *mh_392 = buffer.data(mh + 392);
    const auto *mh_393 = buffer.data(mh + 393);
    const auto *mh_394 = buffer.data(mh + 394);
    const auto *mh_395 = buffer.data(mh + 395);
    const auto *mh_396 = buffer.data(mh + 396);
    const auto *mh_397 = buffer.data(mh + 397);
    const auto *mh_398 = buffer.data(mh + 398);
    const auto *mh_399 = buffer.data(mh + 399);
    const auto *mh_400 = buffer.data(mh + 400);
    const auto *mh_401 = buffer.data(mh + 401);
    const auto *mh_402 = buffer.data(mh + 402);
    const auto *mh_403 = buffer.data(mh + 403);
    const auto *mh_404 = buffer.data(mh + 404);
    const auto *mh_405 = buffer.data(mh + 405);
    const auto *mh_406 = buffer.data(mh + 406);
    const auto *mh_407 = buffer.data(mh + 407);
    const auto *mh_408 = buffer.data(mh + 408);
    const auto *mh_409 = buffer.data(mh + 409);
    const auto *mh_410 = buffer.data(mh + 410);
    const auto *mh_411 = buffer.data(mh + 411);
    const auto *mh_412 = buffer.data(mh + 412);
    const auto *mh_413 = buffer.data(mh + 413);
    const auto *mh_414 = buffer.data(mh + 414);
    const auto *mh_415 = buffer.data(mh + 415);
    const auto *mh_416 = buffer.data(mh + 416);
    const auto *mh_417 = buffer.data(mh + 417);
    const auto *mh_418 = buffer.data(mh + 418);
    const auto *mh_419 = buffer.data(mh + 419);
    const auto *mh_420 = buffer.data(mh + 420);
    const auto *mh_421 = buffer.data(mh + 421);
    const auto *mh_422 = buffer.data(mh + 422);
    const auto *mh_423 = buffer.data(mh + 423);
    const auto *mh_424 = buffer.data(mh + 424);
    const auto *mh_425 = buffer.data(mh + 425);
    const auto *mh_426 = buffer.data(mh + 426);
    const auto *mh_427 = buffer.data(mh + 427);
    const auto *mh_428 = buffer.data(mh + 428);
    const auto *mh_429 = buffer.data(mh + 429);
    const auto *mh_430 = buffer.data(mh + 430);
    const auto *mh_431 = buffer.data(mh + 431);
    const auto *mh_432 = buffer.data(mh + 432);
    const auto *mh_433 = buffer.data(mh + 433);
    const auto *mh_434 = buffer.data(mh + 434);
    const auto *mh_435 = buffer.data(mh + 435);
    const auto *mh_436 = buffer.data(mh + 436);
    const auto *mh_437 = buffer.data(mh + 437);
    const auto *mh_438 = buffer.data(mh + 438);
    const auto *mh_439 = buffer.data(mh + 439);
    const auto *mh_440 = buffer.data(mh + 440);
    const auto *mh_462 = buffer.data(mh + 462);
    const auto *mh_463 = buffer.data(mh + 463);
    const auto *mh_464 = buffer.data(mh + 464);
    const auto *mh_465 = buffer.data(mh + 465);
    const auto *mh_466 = buffer.data(mh + 466);
    const auto *mh_467 = buffer.data(mh + 467);
    const auto *mh_468 = buffer.data(mh + 468);
    const auto *mh_469 = buffer.data(mh + 469);
    const auto *mh_470 = buffer.data(mh + 470);
    const auto *mh_471 = buffer.data(mh + 471);
    const auto *mh_472 = buffer.data(mh + 472);
    const auto *mh_473 = buffer.data(mh + 473);
    const auto *mh_474 = buffer.data(mh + 474);
    const auto *mh_475 = buffer.data(mh + 475);
    const auto *mh_476 = buffer.data(mh + 476);
    const auto *mh_477 = buffer.data(mh + 477);
    const auto *mh_478 = buffer.data(mh + 478);
    const auto *mh_479 = buffer.data(mh + 479);
    const auto *mh_480 = buffer.data(mh + 480);
    const auto *mh_481 = buffer.data(mh + 481);
    const auto *mh_482 = buffer.data(mh + 482);
    const auto *mh_483 = buffer.data(mh + 483);
    const auto *mh_484 = buffer.data(mh + 484);

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, kh_93, kh_94, kh_95, kh_96, kh_97, \
                         mh_282, mh_283, mh_284, mh_285, mh_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -2.0 * kh_93[k]
                   + f_0 * mh_282[k];

        t_178[k] = -2.0 * kh_94[k]
                   + f_0 * mh_283[k];

        t_179[k] = -2.0 * kh_95[k]
                   + f_0 * mh_284[k];

        t_180[k] = -2.0 * kh_96[k]
                   + f_0 * mh_285[k];

        t_181[k] = -2.0 * kh_97[k]
                   + f_0 * mh_286[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, kh_98, kh_99, kh_100, kh_101, \
                         kh_102, mh_287, mh_288, mh_289, mh_290, \
                         mh_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * kh_98[k]
                   + f_0 * mh_287[k];

        t_183[k] = -2.0 * kh_99[k]
                   + f_0 * mh_288[k];

        t_184[k] = -2.0 * kh_100[k]
                   + f_0 * mh_289[k];

        t_185[k] = -2.0 * kh_101[k]
                   + f_0 * mh_290[k];

        t_186[k] = -2.0 * kh_102[k]
                   + f_0 * mh_291[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, kh_103, kh_104, kh_105, kh_106, \
                         kh_107, mh_292, mh_293, mh_294, mh_295, \
                         mh_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * kh_103[k]
                   + f_0 * mh_292[k];

        t_188[k] = -2.0 * kh_104[k]
                   + f_0 * mh_293[k];

        t_189[k] = -3.0 * kh_105[k]
                   + f_0 * mh_294[k];

        t_190[k] = -3.0 * kh_106[k]
                   + f_0 * mh_295[k];

        t_191[k] = -3.0 * kh_107[k]
                   + f_0 * mh_296[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, kh_108, kh_109, kh_110, kh_111, \
                         kh_112, mh_297, mh_298, mh_299, mh_300, \
                         mh_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -3.0 * kh_108[k]
                   + f_0 * mh_297[k];

        t_193[k] = -3.0 * kh_109[k]
                   + f_0 * mh_298[k];

        t_194[k] = -3.0 * kh_110[k]
                   + f_0 * mh_299[k];

        t_195[k] = -3.0 * kh_111[k]
                   + f_0 * mh_300[k];

        t_196[k] = -3.0 * kh_112[k]
                   + f_0 * mh_301[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, kh_113, kh_114, kh_115, kh_116, \
                         kh_117, mh_302, mh_303, mh_304, mh_305, \
                         mh_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -3.0 * kh_113[k]
                   + f_0 * mh_302[k];

        t_198[k] = -3.0 * kh_114[k]
                   + f_0 * mh_303[k];

        t_199[k] = -3.0 * kh_115[k]
                   + f_0 * mh_304[k];

        t_200[k] = -3.0 * kh_116[k]
                   + f_0 * mh_305[k];

        t_201[k] = -3.0 * kh_117[k]
                   + f_0 * mh_306[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, kh_118, kh_119, kh_120, kh_121, \
                         kh_122, mh_307, mh_308, mh_309, mh_310, \
                         mh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -3.0 * kh_118[k]
                   + f_0 * mh_307[k];

        t_203[k] = -3.0 * kh_119[k]
                   + f_0 * mh_308[k];

        t_204[k] = -3.0 * kh_120[k]
                   + f_0 * mh_309[k];

        t_205[k] = -3.0 * kh_121[k]
                   + f_0 * mh_310[k];

        t_206[k] = -3.0 * kh_122[k]
                   + f_0 * mh_311[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, kh_123, kh_124, kh_125, \
                         mh_312, mh_313, mh_314, mh_336, mh_337, \
                         mh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -3.0 * kh_123[k]
                   + f_0 * mh_312[k];

        t_208[k] = -3.0 * kh_124[k]
                   + f_0 * mh_313[k];

        t_209[k] = -3.0 * kh_125[k]
                   + f_0 * mh_314[k];

        t_210[k] = f_0 * mh_336[k];

        t_211[k] = f_0 * mh_337[k];

        t_212[k] = f_0 * mh_338[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, t_218, t_219, t_220, mh_339, \
                         mh_340, mh_341, mh_342, mh_343, mh_344, mh_345, \
                         mh_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_0 * mh_339[k];

        t_214[k] = f_0 * mh_340[k];

        t_215[k] = f_0 * mh_341[k];

        t_216[k] = f_0 * mh_342[k];

        t_217[k] = f_0 * mh_343[k];

        t_218[k] = f_0 * mh_344[k];

        t_219[k] = f_0 * mh_345[k];

        t_220[k] = f_0 * mh_346[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, t_226, t_227, t_228, mh_347, \
                         mh_348, mh_349, mh_350, mh_351, mh_352, mh_353, \
                         mh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_0 * mh_347[k];

        t_222[k] = f_0 * mh_348[k];

        t_223[k] = f_0 * mh_349[k];

        t_224[k] = f_0 * mh_350[k];

        t_225[k] = f_0 * mh_351[k];

        t_226[k] = f_0 * mh_352[k];

        t_227[k] = f_0 * mh_353[k];

        t_228[k] = f_0 * mh_354[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, t_234, kh_126, kh_127, kh_128, \
                         kh_129, mh_355, mh_356, mh_357, mh_358, mh_359, \
                         mh_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_0 * mh_355[k];

        t_230[k] = f_0 * mh_356[k];

        t_231[k] = -kh_126[k]
                   + f_0 * mh_357[k];

        t_232[k] = -kh_127[k]
                   + f_0 * mh_358[k];

        t_233[k] = -kh_128[k]
                   + f_0 * mh_359[k];

        t_234[k] = -kh_129[k]
                   + f_0 * mh_360[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, kh_130, kh_131, kh_132, kh_133, \
                         kh_134, mh_361, mh_362, mh_363, mh_364, \
                         mh_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -kh_130[k]
                   + f_0 * mh_361[k];

        t_236[k] = -kh_131[k]
                   + f_0 * mh_362[k];

        t_237[k] = -kh_132[k]
                   + f_0 * mh_363[k];

        t_238[k] = -kh_133[k]
                   + f_0 * mh_364[k];

        t_239[k] = -kh_134[k]
                   + f_0 * mh_365[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, kh_135, kh_136, kh_137, kh_138, \
                         kh_139, mh_366, mh_367, mh_368, mh_369, \
                         mh_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -kh_135[k]
                   + f_0 * mh_366[k];

        t_241[k] = -kh_136[k]
                   + f_0 * mh_367[k];

        t_242[k] = -kh_137[k]
                   + f_0 * mh_368[k];

        t_243[k] = -kh_138[k]
                   + f_0 * mh_369[k];

        t_244[k] = -kh_139[k]
                   + f_0 * mh_370[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, kh_140, kh_141, kh_142, kh_143, \
                         kh_144, mh_371, mh_372, mh_373, mh_374, \
                         mh_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -kh_140[k]
                   + f_0 * mh_371[k];

        t_246[k] = -kh_141[k]
                   + f_0 * mh_372[k];

        t_247[k] = -kh_142[k]
                   + f_0 * mh_373[k];

        t_248[k] = -kh_143[k]
                   + f_0 * mh_374[k];

        t_249[k] = -kh_144[k]
                   + f_0 * mh_375[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, kh_145, kh_146, kh_147, kh_148, \
                         kh_149, mh_376, mh_377, mh_378, mh_379, \
                         mh_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -kh_145[k]
                   + f_0 * mh_376[k];

        t_251[k] = -kh_146[k]
                   + f_0 * mh_377[k];

        t_252[k] = -2.0 * kh_147[k]
                   + f_0 * mh_378[k];

        t_253[k] = -2.0 * kh_148[k]
                   + f_0 * mh_379[k];

        t_254[k] = -2.0 * kh_149[k]
                   + f_0 * mh_380[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, kh_150, kh_151, kh_152, kh_153, \
                         kh_154, mh_381, mh_382, mh_383, mh_384, \
                         mh_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * kh_150[k]
                   + f_0 * mh_381[k];

        t_256[k] = -2.0 * kh_151[k]
                   + f_0 * mh_382[k];

        t_257[k] = -2.0 * kh_152[k]
                   + f_0 * mh_383[k];

        t_258[k] = -2.0 * kh_153[k]
                   + f_0 * mh_384[k];

        t_259[k] = -2.0 * kh_154[k]
                   + f_0 * mh_385[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, kh_155, kh_156, kh_157, kh_158, \
                         kh_159, mh_386, mh_387, mh_388, mh_389, \
                         mh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * kh_155[k]
                   + f_0 * mh_386[k];

        t_261[k] = -2.0 * kh_156[k]
                   + f_0 * mh_387[k];

        t_262[k] = -2.0 * kh_157[k]
                   + f_0 * mh_388[k];

        t_263[k] = -2.0 * kh_158[k]
                   + f_0 * mh_389[k];

        t_264[k] = -2.0 * kh_159[k]
                   + f_0 * mh_390[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, kh_160, kh_161, kh_162, kh_163, \
                         kh_164, mh_391, mh_392, mh_393, mh_394, \
                         mh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * kh_160[k]
                   + f_0 * mh_391[k];

        t_266[k] = -2.0 * kh_161[k]
                   + f_0 * mh_392[k];

        t_267[k] = -2.0 * kh_162[k]
                   + f_0 * mh_393[k];

        t_268[k] = -2.0 * kh_163[k]
                   + f_0 * mh_394[k];

        t_269[k] = -2.0 * kh_164[k]
                   + f_0 * mh_395[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, kh_165, kh_166, kh_167, kh_168, \
                         kh_169, mh_396, mh_397, mh_398, mh_399, \
                         mh_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * kh_165[k]
                   + f_0 * mh_396[k];

        t_271[k] = -2.0 * kh_166[k]
                   + f_0 * mh_397[k];

        t_272[k] = -2.0 * kh_167[k]
                   + f_0 * mh_398[k];

        t_273[k] = -3.0 * kh_168[k]
                   + f_0 * mh_399[k];

        t_274[k] = -3.0 * kh_169[k]
                   + f_0 * mh_400[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, kh_170, kh_171, kh_172, kh_173, \
                         kh_174, mh_401, mh_402, mh_403, mh_404, \
                         mh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * kh_170[k]
                   + f_0 * mh_401[k];

        t_276[k] = -3.0 * kh_171[k]
                   + f_0 * mh_402[k];

        t_277[k] = -3.0 * kh_172[k]
                   + f_0 * mh_403[k];

        t_278[k] = -3.0 * kh_173[k]
                   + f_0 * mh_404[k];

        t_279[k] = -3.0 * kh_174[k]
                   + f_0 * mh_405[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, kh_175, kh_176, kh_177, kh_178, \
                         kh_179, mh_406, mh_407, mh_408, mh_409, \
                         mh_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * kh_175[k]
                   + f_0 * mh_406[k];

        t_281[k] = -3.0 * kh_176[k]
                   + f_0 * mh_407[k];

        t_282[k] = -3.0 * kh_177[k]
                   + f_0 * mh_408[k];

        t_283[k] = -3.0 * kh_178[k]
                   + f_0 * mh_409[k];

        t_284[k] = -3.0 * kh_179[k]
                   + f_0 * mh_410[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, kh_180, kh_181, kh_182, kh_183, \
                         kh_184, mh_411, mh_412, mh_413, mh_414, \
                         mh_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * kh_180[k]
                   + f_0 * mh_411[k];

        t_286[k] = -3.0 * kh_181[k]
                   + f_0 * mh_412[k];

        t_287[k] = -3.0 * kh_182[k]
                   + f_0 * mh_413[k];

        t_288[k] = -3.0 * kh_183[k]
                   + f_0 * mh_414[k];

        t_289[k] = -3.0 * kh_184[k]
                   + f_0 * mh_415[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, kh_185, kh_186, kh_187, kh_188, \
                         kh_189, mh_416, mh_417, mh_418, mh_419, \
                         mh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * kh_185[k]
                   + f_0 * mh_416[k];

        t_291[k] = -3.0 * kh_186[k]
                   + f_0 * mh_417[k];

        t_292[k] = -3.0 * kh_187[k]
                   + f_0 * mh_418[k];

        t_293[k] = -3.0 * kh_188[k]
                   + f_0 * mh_419[k];

        t_294[k] = -4.0 * kh_189[k]
                   + f_0 * mh_420[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, kh_190, kh_191, kh_192, kh_193, \
                         kh_194, mh_421, mh_422, mh_423, mh_424, \
                         mh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * kh_190[k]
                   + f_0 * mh_421[k];

        t_296[k] = -4.0 * kh_191[k]
                   + f_0 * mh_422[k];

        t_297[k] = -4.0 * kh_192[k]
                   + f_0 * mh_423[k];

        t_298[k] = -4.0 * kh_193[k]
                   + f_0 * mh_424[k];

        t_299[k] = -4.0 * kh_194[k]
                   + f_0 * mh_425[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, kh_195, kh_196, kh_197, kh_198, \
                         kh_199, mh_426, mh_427, mh_428, mh_429, \
                         mh_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * kh_195[k]
                   + f_0 * mh_426[k];

        t_301[k] = -4.0 * kh_196[k]
                   + f_0 * mh_427[k];

        t_302[k] = -4.0 * kh_197[k]
                   + f_0 * mh_428[k];

        t_303[k] = -4.0 * kh_198[k]
                   + f_0 * mh_429[k];

        t_304[k] = -4.0 * kh_199[k]
                   + f_0 * mh_430[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, kh_200, kh_201, kh_202, kh_203, \
                         kh_204, mh_431, mh_432, mh_433, mh_434, \
                         mh_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * kh_200[k]
                   + f_0 * mh_431[k];

        t_306[k] = -4.0 * kh_201[k]
                   + f_0 * mh_432[k];

        t_307[k] = -4.0 * kh_202[k]
                   + f_0 * mh_433[k];

        t_308[k] = -4.0 * kh_203[k]
                   + f_0 * mh_434[k];

        t_309[k] = -4.0 * kh_204[k]
                   + f_0 * mh_435[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, kh_205, kh_206, kh_207, kh_208, \
                         kh_209, mh_436, mh_437, mh_438, mh_439, \
                         mh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -4.0 * kh_205[k]
                   + f_0 * mh_436[k];

        t_311[k] = -4.0 * kh_206[k]
                   + f_0 * mh_437[k];

        t_312[k] = -4.0 * kh_207[k]
                   + f_0 * mh_438[k];

        t_313[k] = -4.0 * kh_208[k]
                   + f_0 * mh_439[k];

        t_314[k] = -4.0 * kh_209[k]
                   + f_0 * mh_440[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, t_320, t_321, t_322, mh_462, \
                         mh_463, mh_464, mh_465, mh_466, mh_467, mh_468, \
                         mh_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_0 * mh_462[k];

        t_316[k] = f_0 * mh_463[k];

        t_317[k] = f_0 * mh_464[k];

        t_318[k] = f_0 * mh_465[k];

        t_319[k] = f_0 * mh_466[k];

        t_320[k] = f_0 * mh_467[k];

        t_321[k] = f_0 * mh_468[k];

        t_322[k] = f_0 * mh_469[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, t_328, t_329, t_330, mh_470, \
                         mh_471, mh_472, mh_473, mh_474, mh_475, mh_476, \
                         mh_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_0 * mh_470[k];

        t_324[k] = f_0 * mh_471[k];

        t_325[k] = f_0 * mh_472[k];

        t_326[k] = f_0 * mh_473[k];

        t_327[k] = f_0 * mh_474[k];

        t_328[k] = f_0 * mh_475[k];

        t_329[k] = f_0 * mh_476[k];

        t_330[k] = f_0 * mh_477[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, t_335, t_336, t_337, kh_210, kh_211, \
                         mh_478, mh_479, mh_480, mh_481, mh_482, mh_483, \
                         mh_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_0 * mh_478[k];

        t_332[k] = f_0 * mh_479[k];

        t_333[k] = f_0 * mh_480[k];

        t_334[k] = f_0 * mh_481[k];

        t_335[k] = f_0 * mh_482[k];

        t_336[k] = -kh_210[k]
                   + f_0 * mh_483[k];

        t_337[k] = -kh_211[k]
                   + f_0 * mh_484[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_338 = buffer.data(target + 338);
    auto *t_339 = buffer.data(target + 339);
    auto *t_340 = buffer.data(target + 340);
    auto *t_341 = buffer.data(target + 341);
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);
    auto *t_347 = buffer.data(target + 347);
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);
    auto *t_380 = buffer.data(target + 380);
    auto *t_381 = buffer.data(target + 381);
    auto *t_382 = buffer.data(target + 382);
    auto *t_383 = buffer.data(target + 383);
    auto *t_384 = buffer.data(target + 384);
    auto *t_385 = buffer.data(target + 385);
    auto *t_386 = buffer.data(target + 386);
    auto *t_387 = buffer.data(target + 387);
    auto *t_388 = buffer.data(target + 388);
    auto *t_389 = buffer.data(target + 389);
    auto *t_390 = buffer.data(target + 390);
    auto *t_391 = buffer.data(target + 391);
    auto *t_392 = buffer.data(target + 392);
    auto *t_393 = buffer.data(target + 393);
    auto *t_394 = buffer.data(target + 394);
    auto *t_395 = buffer.data(target + 395);
    auto *t_396 = buffer.data(target + 396);
    auto *t_397 = buffer.data(target + 397);
    auto *t_398 = buffer.data(target + 398);
    auto *t_399 = buffer.data(target + 399);
    auto *t_400 = buffer.data(target + 400);
    auto *t_401 = buffer.data(target + 401);
    auto *t_402 = buffer.data(target + 402);
    auto *t_403 = buffer.data(target + 403);
    auto *t_404 = buffer.data(target + 404);
    auto *t_405 = buffer.data(target + 405);
    auto *t_406 = buffer.data(target + 406);
    auto *t_407 = buffer.data(target + 407);
    auto *t_408 = buffer.data(target + 408);
    auto *t_409 = buffer.data(target + 409);
    auto *t_410 = buffer.data(target + 410);
    auto *t_411 = buffer.data(target + 411);
    auto *t_412 = buffer.data(target + 412);
    auto *t_413 = buffer.data(target + 413);
    auto *t_414 = buffer.data(target + 414);
    auto *t_415 = buffer.data(target + 415);
    auto *t_416 = buffer.data(target + 416);
    auto *t_417 = buffer.data(target + 417);
    auto *t_418 = buffer.data(target + 418);
    auto *t_419 = buffer.data(target + 419);
    auto *t_420 = buffer.data(target + 420);
    auto *t_421 = buffer.data(target + 421);
    auto *t_422 = buffer.data(target + 422);
    auto *t_423 = buffer.data(target + 423);
    auto *t_424 = buffer.data(target + 424);
    auto *t_425 = buffer.data(target + 425);
    auto *t_426 = buffer.data(target + 426);
    auto *t_427 = buffer.data(target + 427);
    auto *t_428 = buffer.data(target + 428);
    auto *t_429 = buffer.data(target + 429);
    auto *t_430 = buffer.data(target + 430);
    auto *t_431 = buffer.data(target + 431);
    auto *t_432 = buffer.data(target + 432);
    auto *t_433 = buffer.data(target + 433);
    auto *t_434 = buffer.data(target + 434);
    auto *t_435 = buffer.data(target + 435);
    auto *t_436 = buffer.data(target + 436);
    auto *t_437 = buffer.data(target + 437);
    auto *t_438 = buffer.data(target + 438);
    auto *t_439 = buffer.data(target + 439);
    auto *t_440 = buffer.data(target + 440);
    auto *t_441 = buffer.data(target + 441);
    auto *t_442 = buffer.data(target + 442);
    auto *t_443 = buffer.data(target + 443);
    auto *t_444 = buffer.data(target + 444);
    auto *t_445 = buffer.data(target + 445);
    auto *t_446 = buffer.data(target + 446);
    auto *t_447 = buffer.data(target + 447);
    auto *t_448 = buffer.data(target + 448);
    auto *t_449 = buffer.data(target + 449);
    auto *t_450 = buffer.data(target + 450);
    auto *t_451 = buffer.data(target + 451);
    auto *t_452 = buffer.data(target + 452);
    auto *t_453 = buffer.data(target + 453);
    auto *t_454 = buffer.data(target + 454);
    auto *t_455 = buffer.data(target + 455);
    auto *t_456 = buffer.data(target + 456);
    auto *t_457 = buffer.data(target + 457);
    auto *t_458 = buffer.data(target + 458);
    auto *t_459 = buffer.data(target + 459);
    auto *t_460 = buffer.data(target + 460);
    auto *t_461 = buffer.data(target + 461);
    auto *t_462 = buffer.data(target + 462);
    auto *t_463 = buffer.data(target + 463);
    auto *t_464 = buffer.data(target + 464);
    auto *t_465 = buffer.data(target + 465);
    auto *t_466 = buffer.data(target + 466);
    auto *t_467 = buffer.data(target + 467);
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);
    auto *t_474 = buffer.data(target + 474);
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);

    const auto *kh_212 = buffer.data(kh + 212);
    const auto *kh_213 = buffer.data(kh + 213);
    const auto *kh_214 = buffer.data(kh + 214);
    const auto *kh_215 = buffer.data(kh + 215);
    const auto *kh_216 = buffer.data(kh + 216);
    const auto *kh_217 = buffer.data(kh + 217);
    const auto *kh_218 = buffer.data(kh + 218);
    const auto *kh_219 = buffer.data(kh + 219);
    const auto *kh_220 = buffer.data(kh + 220);
    const auto *kh_221 = buffer.data(kh + 221);
    const auto *kh_222 = buffer.data(kh + 222);
    const auto *kh_223 = buffer.data(kh + 223);
    const auto *kh_224 = buffer.data(kh + 224);
    const auto *kh_225 = buffer.data(kh + 225);
    const auto *kh_226 = buffer.data(kh + 226);
    const auto *kh_227 = buffer.data(kh + 227);
    const auto *kh_228 = buffer.data(kh + 228);
    const auto *kh_229 = buffer.data(kh + 229);
    const auto *kh_230 = buffer.data(kh + 230);
    const auto *kh_231 = buffer.data(kh + 231);
    const auto *kh_232 = buffer.data(kh + 232);
    const auto *kh_233 = buffer.data(kh + 233);
    const auto *kh_234 = buffer.data(kh + 234);
    const auto *kh_235 = buffer.data(kh + 235);
    const auto *kh_236 = buffer.data(kh + 236);
    const auto *kh_237 = buffer.data(kh + 237);
    const auto *kh_238 = buffer.data(kh + 238);
    const auto *kh_239 = buffer.data(kh + 239);
    const auto *kh_240 = buffer.data(kh + 240);
    const auto *kh_241 = buffer.data(kh + 241);
    const auto *kh_242 = buffer.data(kh + 242);
    const auto *kh_243 = buffer.data(kh + 243);
    const auto *kh_244 = buffer.data(kh + 244);
    const auto *kh_245 = buffer.data(kh + 245);
    const auto *kh_246 = buffer.data(kh + 246);
    const auto *kh_247 = buffer.data(kh + 247);
    const auto *kh_248 = buffer.data(kh + 248);
    const auto *kh_249 = buffer.data(kh + 249);
    const auto *kh_250 = buffer.data(kh + 250);
    const auto *kh_251 = buffer.data(kh + 251);
    const auto *kh_252 = buffer.data(kh + 252);
    const auto *kh_253 = buffer.data(kh + 253);
    const auto *kh_254 = buffer.data(kh + 254);
    const auto *kh_255 = buffer.data(kh + 255);
    const auto *kh_256 = buffer.data(kh + 256);
    const auto *kh_257 = buffer.data(kh + 257);
    const auto *kh_258 = buffer.data(kh + 258);
    const auto *kh_259 = buffer.data(kh + 259);
    const auto *kh_260 = buffer.data(kh + 260);
    const auto *kh_261 = buffer.data(kh + 261);
    const auto *kh_262 = buffer.data(kh + 262);
    const auto *kh_263 = buffer.data(kh + 263);
    const auto *kh_264 = buffer.data(kh + 264);
    const auto *kh_265 = buffer.data(kh + 265);
    const auto *kh_266 = buffer.data(kh + 266);
    const auto *kh_267 = buffer.data(kh + 267);
    const auto *kh_268 = buffer.data(kh + 268);
    const auto *kh_269 = buffer.data(kh + 269);
    const auto *kh_270 = buffer.data(kh + 270);
    const auto *kh_271 = buffer.data(kh + 271);
    const auto *kh_272 = buffer.data(kh + 272);
    const auto *kh_273 = buffer.data(kh + 273);
    const auto *kh_274 = buffer.data(kh + 274);
    const auto *kh_275 = buffer.data(kh + 275);
    const auto *kh_276 = buffer.data(kh + 276);
    const auto *kh_277 = buffer.data(kh + 277);
    const auto *kh_278 = buffer.data(kh + 278);
    const auto *kh_279 = buffer.data(kh + 279);
    const auto *kh_280 = buffer.data(kh + 280);
    const auto *kh_281 = buffer.data(kh + 281);
    const auto *kh_282 = buffer.data(kh + 282);
    const auto *kh_283 = buffer.data(kh + 283);
    const auto *kh_284 = buffer.data(kh + 284);
    const auto *kh_285 = buffer.data(kh + 285);
    const auto *kh_286 = buffer.data(kh + 286);
    const auto *kh_287 = buffer.data(kh + 287);
    const auto *kh_288 = buffer.data(kh + 288);
    const auto *kh_289 = buffer.data(kh + 289);
    const auto *kh_290 = buffer.data(kh + 290);
    const auto *kh_291 = buffer.data(kh + 291);
    const auto *kh_292 = buffer.data(kh + 292);
    const auto *kh_293 = buffer.data(kh + 293);
    const auto *kh_294 = buffer.data(kh + 294);
    const auto *kh_295 = buffer.data(kh + 295);
    const auto *kh_296 = buffer.data(kh + 296);
    const auto *kh_297 = buffer.data(kh + 297);
    const auto *kh_298 = buffer.data(kh + 298);
    const auto *kh_299 = buffer.data(kh + 299);
    const auto *kh_300 = buffer.data(kh + 300);
    const auto *kh_301 = buffer.data(kh + 301);
    const auto *kh_302 = buffer.data(kh + 302);
    const auto *kh_303 = buffer.data(kh + 303);
    const auto *kh_304 = buffer.data(kh + 304);
    const auto *kh_305 = buffer.data(kh + 305);
    const auto *kh_306 = buffer.data(kh + 306);
    const auto *kh_307 = buffer.data(kh + 307);
    const auto *kh_308 = buffer.data(kh + 308);
    const auto *kh_309 = buffer.data(kh + 309);
    const auto *kh_310 = buffer.data(kh + 310);
    const auto *kh_311 = buffer.data(kh + 311);
    const auto *kh_312 = buffer.data(kh + 312);
    const auto *kh_313 = buffer.data(kh + 313);
    const auto *kh_314 = buffer.data(kh + 314);
    const auto *kh_315 = buffer.data(kh + 315);
    const auto *kh_316 = buffer.data(kh + 316);
    const auto *kh_317 = buffer.data(kh + 317);
    const auto *kh_318 = buffer.data(kh + 318);
    const auto *kh_319 = buffer.data(kh + 319);
    const auto *kh_320 = buffer.data(kh + 320);
    const auto *kh_321 = buffer.data(kh + 321);
    const auto *kh_322 = buffer.data(kh + 322);
    const auto *kh_323 = buffer.data(kh + 323);
    const auto *kh_324 = buffer.data(kh + 324);
    const auto *kh_325 = buffer.data(kh + 325);
    const auto *kh_326 = buffer.data(kh + 326);
    const auto *kh_327 = buffer.data(kh + 327);
    const auto *kh_328 = buffer.data(kh + 328);
    const auto *kh_329 = buffer.data(kh + 329);
    const auto *kh_330 = buffer.data(kh + 330);
    const auto *kh_331 = buffer.data(kh + 331);
    const auto *kh_332 = buffer.data(kh + 332);
    const auto *kh_333 = buffer.data(kh + 333);
    const auto *kh_334 = buffer.data(kh + 334);
    const auto *kh_335 = buffer.data(kh + 335);
    const auto *kh_336 = buffer.data(kh + 336);
    const auto *kh_337 = buffer.data(kh + 337);
    const auto *kh_338 = buffer.data(kh + 338);
    const auto *kh_339 = buffer.data(kh + 339);
    const auto *kh_340 = buffer.data(kh + 340);
    const auto *kh_341 = buffer.data(kh + 341);
    const auto *kh_342 = buffer.data(kh + 342);
    const auto *kh_343 = buffer.data(kh + 343);
    const auto *kh_344 = buffer.data(kh + 344);
    const auto *kh_345 = buffer.data(kh + 345);
    const auto *kh_346 = buffer.data(kh + 346);
    const auto *kh_347 = buffer.data(kh + 347);
    const auto *kh_348 = buffer.data(kh + 348);

    const auto *mh_485 = buffer.data(mh + 485);
    const auto *mh_486 = buffer.data(mh + 486);
    const auto *mh_487 = buffer.data(mh + 487);
    const auto *mh_488 = buffer.data(mh + 488);
    const auto *mh_489 = buffer.data(mh + 489);
    const auto *mh_490 = buffer.data(mh + 490);
    const auto *mh_491 = buffer.data(mh + 491);
    const auto *mh_492 = buffer.data(mh + 492);
    const auto *mh_493 = buffer.data(mh + 493);
    const auto *mh_494 = buffer.data(mh + 494);
    const auto *mh_495 = buffer.data(mh + 495);
    const auto *mh_496 = buffer.data(mh + 496);
    const auto *mh_497 = buffer.data(mh + 497);
    const auto *mh_498 = buffer.data(mh + 498);
    const auto *mh_499 = buffer.data(mh + 499);
    const auto *mh_500 = buffer.data(mh + 500);
    const auto *mh_501 = buffer.data(mh + 501);
    const auto *mh_502 = buffer.data(mh + 502);
    const auto *mh_503 = buffer.data(mh + 503);
    const auto *mh_504 = buffer.data(mh + 504);
    const auto *mh_505 = buffer.data(mh + 505);
    const auto *mh_506 = buffer.data(mh + 506);
    const auto *mh_507 = buffer.data(mh + 507);
    const auto *mh_508 = buffer.data(mh + 508);
    const auto *mh_509 = buffer.data(mh + 509);
    const auto *mh_510 = buffer.data(mh + 510);
    const auto *mh_511 = buffer.data(mh + 511);
    const auto *mh_512 = buffer.data(mh + 512);
    const auto *mh_513 = buffer.data(mh + 513);
    const auto *mh_514 = buffer.data(mh + 514);
    const auto *mh_515 = buffer.data(mh + 515);
    const auto *mh_516 = buffer.data(mh + 516);
    const auto *mh_517 = buffer.data(mh + 517);
    const auto *mh_518 = buffer.data(mh + 518);
    const auto *mh_519 = buffer.data(mh + 519);
    const auto *mh_520 = buffer.data(mh + 520);
    const auto *mh_521 = buffer.data(mh + 521);
    const auto *mh_522 = buffer.data(mh + 522);
    const auto *mh_523 = buffer.data(mh + 523);
    const auto *mh_524 = buffer.data(mh + 524);
    const auto *mh_525 = buffer.data(mh + 525);
    const auto *mh_526 = buffer.data(mh + 526);
    const auto *mh_527 = buffer.data(mh + 527);
    const auto *mh_528 = buffer.data(mh + 528);
    const auto *mh_529 = buffer.data(mh + 529);
    const auto *mh_530 = buffer.data(mh + 530);
    const auto *mh_531 = buffer.data(mh + 531);
    const auto *mh_532 = buffer.data(mh + 532);
    const auto *mh_533 = buffer.data(mh + 533);
    const auto *mh_534 = buffer.data(mh + 534);
    const auto *mh_535 = buffer.data(mh + 535);
    const auto *mh_536 = buffer.data(mh + 536);
    const auto *mh_537 = buffer.data(mh + 537);
    const auto *mh_538 = buffer.data(mh + 538);
    const auto *mh_539 = buffer.data(mh + 539);
    const auto *mh_540 = buffer.data(mh + 540);
    const auto *mh_541 = buffer.data(mh + 541);
    const auto *mh_542 = buffer.data(mh + 542);
    const auto *mh_543 = buffer.data(mh + 543);
    const auto *mh_544 = buffer.data(mh + 544);
    const auto *mh_545 = buffer.data(mh + 545);
    const auto *mh_546 = buffer.data(mh + 546);
    const auto *mh_547 = buffer.data(mh + 547);
    const auto *mh_548 = buffer.data(mh + 548);
    const auto *mh_549 = buffer.data(mh + 549);
    const auto *mh_550 = buffer.data(mh + 550);
    const auto *mh_551 = buffer.data(mh + 551);
    const auto *mh_552 = buffer.data(mh + 552);
    const auto *mh_553 = buffer.data(mh + 553);
    const auto *mh_554 = buffer.data(mh + 554);
    const auto *mh_555 = buffer.data(mh + 555);
    const auto *mh_556 = buffer.data(mh + 556);
    const auto *mh_557 = buffer.data(mh + 557);
    const auto *mh_558 = buffer.data(mh + 558);
    const auto *mh_559 = buffer.data(mh + 559);
    const auto *mh_560 = buffer.data(mh + 560);
    const auto *mh_561 = buffer.data(mh + 561);
    const auto *mh_562 = buffer.data(mh + 562);
    const auto *mh_563 = buffer.data(mh + 563);
    const auto *mh_564 = buffer.data(mh + 564);
    const auto *mh_565 = buffer.data(mh + 565);
    const auto *mh_566 = buffer.data(mh + 566);
    const auto *mh_567 = buffer.data(mh + 567);
    const auto *mh_568 = buffer.data(mh + 568);
    const auto *mh_569 = buffer.data(mh + 569);
    const auto *mh_570 = buffer.data(mh + 570);
    const auto *mh_571 = buffer.data(mh + 571);
    const auto *mh_572 = buffer.data(mh + 572);
    const auto *mh_573 = buffer.data(mh + 573);
    const auto *mh_574 = buffer.data(mh + 574);
    const auto *mh_575 = buffer.data(mh + 575);
    const auto *mh_576 = buffer.data(mh + 576);
    const auto *mh_577 = buffer.data(mh + 577);
    const auto *mh_578 = buffer.data(mh + 578);
    const auto *mh_579 = buffer.data(mh + 579);
    const auto *mh_580 = buffer.data(mh + 580);
    const auto *mh_581 = buffer.data(mh + 581);
    const auto *mh_582 = buffer.data(mh + 582);
    const auto *mh_583 = buffer.data(mh + 583);
    const auto *mh_584 = buffer.data(mh + 584);
    const auto *mh_585 = buffer.data(mh + 585);
    const auto *mh_586 = buffer.data(mh + 586);
    const auto *mh_587 = buffer.data(mh + 587);
    const auto *mh_609 = buffer.data(mh + 609);
    const auto *mh_610 = buffer.data(mh + 610);
    const auto *mh_611 = buffer.data(mh + 611);
    const auto *mh_612 = buffer.data(mh + 612);
    const auto *mh_613 = buffer.data(mh + 613);
    const auto *mh_614 = buffer.data(mh + 614);
    const auto *mh_615 = buffer.data(mh + 615);
    const auto *mh_616 = buffer.data(mh + 616);
    const auto *mh_617 = buffer.data(mh + 617);
    const auto *mh_618 = buffer.data(mh + 618);
    const auto *mh_619 = buffer.data(mh + 619);
    const auto *mh_620 = buffer.data(mh + 620);
    const auto *mh_621 = buffer.data(mh + 621);
    const auto *mh_622 = buffer.data(mh + 622);
    const auto *mh_623 = buffer.data(mh + 623);
    const auto *mh_624 = buffer.data(mh + 624);
    const auto *mh_625 = buffer.data(mh + 625);
    const auto *mh_626 = buffer.data(mh + 626);
    const auto *mh_627 = buffer.data(mh + 627);
    const auto *mh_628 = buffer.data(mh + 628);
    const auto *mh_629 = buffer.data(mh + 629);
    const auto *mh_630 = buffer.data(mh + 630);
    const auto *mh_631 = buffer.data(mh + 631);
    const auto *mh_632 = buffer.data(mh + 632);
    const auto *mh_633 = buffer.data(mh + 633);
    const auto *mh_634 = buffer.data(mh + 634);
    const auto *mh_635 = buffer.data(mh + 635);
    const auto *mh_636 = buffer.data(mh + 636);
    const auto *mh_637 = buffer.data(mh + 637);
    const auto *mh_638 = buffer.data(mh + 638);
    const auto *mh_639 = buffer.data(mh + 639);
    const auto *mh_640 = buffer.data(mh + 640);
    const auto *mh_641 = buffer.data(mh + 641);
    const auto *mh_642 = buffer.data(mh + 642);
    const auto *mh_643 = buffer.data(mh + 643);
    const auto *mh_644 = buffer.data(mh + 644);
    const auto *mh_645 = buffer.data(mh + 645);
    const auto *mh_646 = buffer.data(mh + 646);
    const auto *mh_647 = buffer.data(mh + 647);
    const auto *mh_648 = buffer.data(mh + 648);
    const auto *mh_649 = buffer.data(mh + 649);
    const auto *mh_650 = buffer.data(mh + 650);
    const auto *mh_651 = buffer.data(mh + 651);
    const auto *mh_652 = buffer.data(mh + 652);
    const auto *mh_653 = buffer.data(mh + 653);
    const auto *mh_654 = buffer.data(mh + 654);
    const auto *mh_655 = buffer.data(mh + 655);
    const auto *mh_656 = buffer.data(mh + 656);
    const auto *mh_657 = buffer.data(mh + 657);
    const auto *mh_658 = buffer.data(mh + 658);
    const auto *mh_659 = buffer.data(mh + 659);
    const auto *mh_660 = buffer.data(mh + 660);
    const auto *mh_661 = buffer.data(mh + 661);
    const auto *mh_662 = buffer.data(mh + 662);
    const auto *mh_663 = buffer.data(mh + 663);

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, kh_212, kh_213, kh_214, kh_215, \
                         kh_216, mh_485, mh_486, mh_487, mh_488, \
                         mh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = -kh_212[k]
                   + f_0 * mh_485[k];

        t_339[k] = -kh_213[k]
                   + f_0 * mh_486[k];

        t_340[k] = -kh_214[k]
                   + f_0 * mh_487[k];

        t_341[k] = -kh_215[k]
                   + f_0 * mh_488[k];

        t_342[k] = -kh_216[k]
                   + f_0 * mh_489[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, kh_217, kh_218, kh_219, kh_220, \
                         kh_221, mh_490, mh_491, mh_492, mh_493, \
                         mh_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = -kh_217[k]
                   + f_0 * mh_490[k];

        t_344[k] = -kh_218[k]
                   + f_0 * mh_491[k];

        t_345[k] = -kh_219[k]
                   + f_0 * mh_492[k];

        t_346[k] = -kh_220[k]
                   + f_0 * mh_493[k];

        t_347[k] = -kh_221[k]
                   + f_0 * mh_494[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, kh_222, kh_223, kh_224, kh_225, \
                         kh_226, mh_495, mh_496, mh_497, mh_498, \
                         mh_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = -kh_222[k]
                   + f_0 * mh_495[k];

        t_349[k] = -kh_223[k]
                   + f_0 * mh_496[k];

        t_350[k] = -kh_224[k]
                   + f_0 * mh_497[k];

        t_351[k] = -kh_225[k]
                   + f_0 * mh_498[k];

        t_352[k] = -kh_226[k]
                   + f_0 * mh_499[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, kh_227, kh_228, kh_229, kh_230, \
                         kh_231, mh_500, mh_501, mh_502, mh_503, \
                         mh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = -kh_227[k]
                   + f_0 * mh_500[k];

        t_354[k] = -kh_228[k]
                   + f_0 * mh_501[k];

        t_355[k] = -kh_229[k]
                   + f_0 * mh_502[k];

        t_356[k] = -kh_230[k]
                   + f_0 * mh_503[k];

        t_357[k] = -2.0 * kh_231[k]
                   + f_0 * mh_504[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, kh_232, kh_233, kh_234, kh_235, \
                         kh_236, mh_505, mh_506, mh_507, mh_508, \
                         mh_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = -2.0 * kh_232[k]
                   + f_0 * mh_505[k];

        t_359[k] = -2.0 * kh_233[k]
                   + f_0 * mh_506[k];

        t_360[k] = -2.0 * kh_234[k]
                   + f_0 * mh_507[k];

        t_361[k] = -2.0 * kh_235[k]
                   + f_0 * mh_508[k];

        t_362[k] = -2.0 * kh_236[k]
                   + f_0 * mh_509[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, kh_237, kh_238, kh_239, kh_240, \
                         kh_241, mh_510, mh_511, mh_512, mh_513, \
                         mh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = -2.0 * kh_237[k]
                   + f_0 * mh_510[k];

        t_364[k] = -2.0 * kh_238[k]
                   + f_0 * mh_511[k];

        t_365[k] = -2.0 * kh_239[k]
                   + f_0 * mh_512[k];

        t_366[k] = -2.0 * kh_240[k]
                   + f_0 * mh_513[k];

        t_367[k] = -2.0 * kh_241[k]
                   + f_0 * mh_514[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, kh_242, kh_243, kh_244, kh_245, \
                         kh_246, mh_515, mh_516, mh_517, mh_518, \
                         mh_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = -2.0 * kh_242[k]
                   + f_0 * mh_515[k];

        t_369[k] = -2.0 * kh_243[k]
                   + f_0 * mh_516[k];

        t_370[k] = -2.0 * kh_244[k]
                   + f_0 * mh_517[k];

        t_371[k] = -2.0 * kh_245[k]
                   + f_0 * mh_518[k];

        t_372[k] = -2.0 * kh_246[k]
                   + f_0 * mh_519[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, kh_247, kh_248, kh_249, kh_250, \
                         kh_251, mh_520, mh_521, mh_522, mh_523, \
                         mh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = -2.0 * kh_247[k]
                   + f_0 * mh_520[k];

        t_374[k] = -2.0 * kh_248[k]
                   + f_0 * mh_521[k];

        t_375[k] = -2.0 * kh_249[k]
                   + f_0 * mh_522[k];

        t_376[k] = -2.0 * kh_250[k]
                   + f_0 * mh_523[k];

        t_377[k] = -2.0 * kh_251[k]
                   + f_0 * mh_524[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, kh_252, kh_253, kh_254, kh_255, \
                         kh_256, mh_525, mh_526, mh_527, mh_528, \
                         mh_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = -3.0 * kh_252[k]
                   + f_0 * mh_525[k];

        t_379[k] = -3.0 * kh_253[k]
                   + f_0 * mh_526[k];

        t_380[k] = -3.0 * kh_254[k]
                   + f_0 * mh_527[k];

        t_381[k] = -3.0 * kh_255[k]
                   + f_0 * mh_528[k];

        t_382[k] = -3.0 * kh_256[k]
                   + f_0 * mh_529[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, kh_257, kh_258, kh_259, kh_260, \
                         kh_261, mh_530, mh_531, mh_532, mh_533, \
                         mh_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = -3.0 * kh_257[k]
                   + f_0 * mh_530[k];

        t_384[k] = -3.0 * kh_258[k]
                   + f_0 * mh_531[k];

        t_385[k] = -3.0 * kh_259[k]
                   + f_0 * mh_532[k];

        t_386[k] = -3.0 * kh_260[k]
                   + f_0 * mh_533[k];

        t_387[k] = -3.0 * kh_261[k]
                   + f_0 * mh_534[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, kh_262, kh_263, kh_264, kh_265, \
                         kh_266, mh_535, mh_536, mh_537, mh_538, \
                         mh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = -3.0 * kh_262[k]
                   + f_0 * mh_535[k];

        t_389[k] = -3.0 * kh_263[k]
                   + f_0 * mh_536[k];

        t_390[k] = -3.0 * kh_264[k]
                   + f_0 * mh_537[k];

        t_391[k] = -3.0 * kh_265[k]
                   + f_0 * mh_538[k];

        t_392[k] = -3.0 * kh_266[k]
                   + f_0 * mh_539[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, t_397, kh_267, kh_268, kh_269, kh_270, \
                         kh_271, mh_540, mh_541, mh_542, mh_543, \
                         mh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = -3.0 * kh_267[k]
                   + f_0 * mh_540[k];

        t_394[k] = -3.0 * kh_268[k]
                   + f_0 * mh_541[k];

        t_395[k] = -3.0 * kh_269[k]
                   + f_0 * mh_542[k];

        t_396[k] = -3.0 * kh_270[k]
                   + f_0 * mh_543[k];

        t_397[k] = -3.0 * kh_271[k]
                   + f_0 * mh_544[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, t_402, kh_272, kh_273, kh_274, kh_275, \
                         kh_276, mh_545, mh_546, mh_547, mh_548, \
                         mh_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = -3.0 * kh_272[k]
                   + f_0 * mh_545[k];

        t_399[k] = -4.0 * kh_273[k]
                   + f_0 * mh_546[k];

        t_400[k] = -4.0 * kh_274[k]
                   + f_0 * mh_547[k];

        t_401[k] = -4.0 * kh_275[k]
                   + f_0 * mh_548[k];

        t_402[k] = -4.0 * kh_276[k]
                   + f_0 * mh_549[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, t_407, kh_277, kh_278, kh_279, kh_280, \
                         kh_281, mh_550, mh_551, mh_552, mh_553, \
                         mh_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = -4.0 * kh_277[k]
                   + f_0 * mh_550[k];

        t_404[k] = -4.0 * kh_278[k]
                   + f_0 * mh_551[k];

        t_405[k] = -4.0 * kh_279[k]
                   + f_0 * mh_552[k];

        t_406[k] = -4.0 * kh_280[k]
                   + f_0 * mh_553[k];

        t_407[k] = -4.0 * kh_281[k]
                   + f_0 * mh_554[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, kh_282, kh_283, kh_284, kh_285, \
                         kh_286, mh_555, mh_556, mh_557, mh_558, \
                         mh_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = -4.0 * kh_282[k]
                   + f_0 * mh_555[k];

        t_409[k] = -4.0 * kh_283[k]
                   + f_0 * mh_556[k];

        t_410[k] = -4.0 * kh_284[k]
                   + f_0 * mh_557[k];

        t_411[k] = -4.0 * kh_285[k]
                   + f_0 * mh_558[k];

        t_412[k] = -4.0 * kh_286[k]
                   + f_0 * mh_559[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, kh_287, kh_288, kh_289, kh_290, \
                         kh_291, mh_560, mh_561, mh_562, mh_563, \
                         mh_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = -4.0 * kh_287[k]
                   + f_0 * mh_560[k];

        t_414[k] = -4.0 * kh_288[k]
                   + f_0 * mh_561[k];

        t_415[k] = -4.0 * kh_289[k]
                   + f_0 * mh_562[k];

        t_416[k] = -4.0 * kh_290[k]
                   + f_0 * mh_563[k];

        t_417[k] = -4.0 * kh_291[k]
                   + f_0 * mh_564[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, t_422, kh_292, kh_293, kh_294, kh_295, \
                         kh_296, mh_565, mh_566, mh_567, mh_568, \
                         mh_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = -4.0 * kh_292[k]
                   + f_0 * mh_565[k];

        t_419[k] = -4.0 * kh_293[k]
                   + f_0 * mh_566[k];

        t_420[k] = -5.0 * kh_294[k]
                   + f_0 * mh_567[k];

        t_421[k] = -5.0 * kh_295[k]
                   + f_0 * mh_568[k];

        t_422[k] = -5.0 * kh_296[k]
                   + f_0 * mh_569[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, kh_297, kh_298, kh_299, kh_300, \
                         kh_301, mh_570, mh_571, mh_572, mh_573, \
                         mh_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = -5.0 * kh_297[k]
                   + f_0 * mh_570[k];

        t_424[k] = -5.0 * kh_298[k]
                   + f_0 * mh_571[k];

        t_425[k] = -5.0 * kh_299[k]
                   + f_0 * mh_572[k];

        t_426[k] = -5.0 * kh_300[k]
                   + f_0 * mh_573[k];

        t_427[k] = -5.0 * kh_301[k]
                   + f_0 * mh_574[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, kh_302, kh_303, kh_304, kh_305, \
                         kh_306, mh_575, mh_576, mh_577, mh_578, \
                         mh_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = -5.0 * kh_302[k]
                   + f_0 * mh_575[k];

        t_429[k] = -5.0 * kh_303[k]
                   + f_0 * mh_576[k];

        t_430[k] = -5.0 * kh_304[k]
                   + f_0 * mh_577[k];

        t_431[k] = -5.0 * kh_305[k]
                   + f_0 * mh_578[k];

        t_432[k] = -5.0 * kh_306[k]
                   + f_0 * mh_579[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, kh_307, kh_308, kh_309, kh_310, \
                         kh_311, mh_580, mh_581, mh_582, mh_583, \
                         mh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = -5.0 * kh_307[k]
                   + f_0 * mh_580[k];

        t_434[k] = -5.0 * kh_308[k]
                   + f_0 * mh_581[k];

        t_435[k] = -5.0 * kh_309[k]
                   + f_0 * mh_582[k];

        t_436[k] = -5.0 * kh_310[k]
                   + f_0 * mh_583[k];

        t_437[k] = -5.0 * kh_311[k]
                   + f_0 * mh_584[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, t_442, t_443, kh_312, kh_313, kh_314, \
                         mh_585, mh_586, mh_587, mh_609, mh_610, \
                         mh_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = -5.0 * kh_312[k]
                   + f_0 * mh_585[k];

        t_439[k] = -5.0 * kh_313[k]
                   + f_0 * mh_586[k];

        t_440[k] = -5.0 * kh_314[k]
                   + f_0 * mh_587[k];

        t_441[k] = f_0 * mh_609[k];

        t_442[k] = f_0 * mh_610[k];

        t_443[k] = f_0 * mh_611[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, t_449, t_450, t_451, mh_612, \
                         mh_613, mh_614, mh_615, mh_616, mh_617, mh_618, \
                         mh_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_0 * mh_612[k];

        t_445[k] = f_0 * mh_613[k];

        t_446[k] = f_0 * mh_614[k];

        t_447[k] = f_0 * mh_615[k];

        t_448[k] = f_0 * mh_616[k];

        t_449[k] = f_0 * mh_617[k];

        t_450[k] = f_0 * mh_618[k];

        t_451[k] = f_0 * mh_619[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, t_457, t_458, t_459, mh_620, \
                         mh_621, mh_622, mh_623, mh_624, mh_625, mh_626, \
                         mh_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_0 * mh_620[k];

        t_453[k] = f_0 * mh_621[k];

        t_454[k] = f_0 * mh_622[k];

        t_455[k] = f_0 * mh_623[k];

        t_456[k] = f_0 * mh_624[k];

        t_457[k] = f_0 * mh_625[k];

        t_458[k] = f_0 * mh_626[k];

        t_459[k] = f_0 * mh_627[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, t_465, kh_315, kh_316, kh_317, \
                         kh_318, mh_628, mh_629, mh_630, mh_631, mh_632, \
                         mh_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_0 * mh_628[k];

        t_461[k] = f_0 * mh_629[k];

        t_462[k] = -kh_315[k]
                   + f_0 * mh_630[k];

        t_463[k] = -kh_316[k]
                   + f_0 * mh_631[k];

        t_464[k] = -kh_317[k]
                   + f_0 * mh_632[k];

        t_465[k] = -kh_318[k]
                   + f_0 * mh_633[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, kh_319, kh_320, kh_321, kh_322, \
                         kh_323, mh_634, mh_635, mh_636, mh_637, \
                         mh_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = -kh_319[k]
                   + f_0 * mh_634[k];

        t_467[k] = -kh_320[k]
                   + f_0 * mh_635[k];

        t_468[k] = -kh_321[k]
                   + f_0 * mh_636[k];

        t_469[k] = -kh_322[k]
                   + f_0 * mh_637[k];

        t_470[k] = -kh_323[k]
                   + f_0 * mh_638[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, t_475, kh_324, kh_325, kh_326, kh_327, \
                         kh_328, mh_639, mh_640, mh_641, mh_642, \
                         mh_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = -kh_324[k]
                   + f_0 * mh_639[k];

        t_472[k] = -kh_325[k]
                   + f_0 * mh_640[k];

        t_473[k] = -kh_326[k]
                   + f_0 * mh_641[k];

        t_474[k] = -kh_327[k]
                   + f_0 * mh_642[k];

        t_475[k] = -kh_328[k]
                   + f_0 * mh_643[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, kh_329, kh_330, kh_331, kh_332, \
                         kh_333, mh_644, mh_645, mh_646, mh_647, \
                         mh_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = -kh_329[k]
                   + f_0 * mh_644[k];

        t_477[k] = -kh_330[k]
                   + f_0 * mh_645[k];

        t_478[k] = -kh_331[k]
                   + f_0 * mh_646[k];

        t_479[k] = -kh_332[k]
                   + f_0 * mh_647[k];

        t_480[k] = -kh_333[k]
                   + f_0 * mh_648[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, kh_334, kh_335, kh_336, kh_337, \
                         kh_338, mh_649, mh_650, mh_651, mh_652, \
                         mh_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = -kh_334[k]
                   + f_0 * mh_649[k];

        t_482[k] = -kh_335[k]
                   + f_0 * mh_650[k];

        t_483[k] = -2.0 * kh_336[k]
                   + f_0 * mh_651[k];

        t_484[k] = -2.0 * kh_337[k]
                   + f_0 * mh_652[k];

        t_485[k] = -2.0 * kh_338[k]
                   + f_0 * mh_653[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, kh_339, kh_340, kh_341, kh_342, \
                         kh_343, mh_654, mh_655, mh_656, mh_657, \
                         mh_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = -2.0 * kh_339[k]
                   + f_0 * mh_654[k];

        t_487[k] = -2.0 * kh_340[k]
                   + f_0 * mh_655[k];

        t_488[k] = -2.0 * kh_341[k]
                   + f_0 * mh_656[k];

        t_489[k] = -2.0 * kh_342[k]
                   + f_0 * mh_657[k];

        t_490[k] = -2.0 * kh_343[k]
                   + f_0 * mh_658[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, kh_344, kh_345, kh_346, kh_347, \
                         kh_348, mh_659, mh_660, mh_661, mh_662, \
                         mh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = -2.0 * kh_344[k]
                   + f_0 * mh_659[k];

        t_492[k] = -2.0 * kh_345[k]
                   + f_0 * mh_660[k];

        t_493[k] = -2.0 * kh_346[k]
                   + f_0 * mh_661[k];

        t_494[k] = -2.0 * kh_347[k]
                   + f_0 * mh_662[k];

        t_495[k] = -2.0 * kh_348[k]
                   + f_0 * mh_663[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);
    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);
    auto *t_500 = buffer.data(target + 500);
    auto *t_501 = buffer.data(target + 501);
    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);
    auto *t_504 = buffer.data(target + 504);
    auto *t_505 = buffer.data(target + 505);
    auto *t_506 = buffer.data(target + 506);
    auto *t_507 = buffer.data(target + 507);
    auto *t_508 = buffer.data(target + 508);
    auto *t_509 = buffer.data(target + 509);
    auto *t_510 = buffer.data(target + 510);
    auto *t_511 = buffer.data(target + 511);
    auto *t_512 = buffer.data(target + 512);
    auto *t_513 = buffer.data(target + 513);
    auto *t_514 = buffer.data(target + 514);
    auto *t_515 = buffer.data(target + 515);
    auto *t_516 = buffer.data(target + 516);
    auto *t_517 = buffer.data(target + 517);
    auto *t_518 = buffer.data(target + 518);
    auto *t_519 = buffer.data(target + 519);
    auto *t_520 = buffer.data(target + 520);
    auto *t_521 = buffer.data(target + 521);
    auto *t_522 = buffer.data(target + 522);
    auto *t_523 = buffer.data(target + 523);
    auto *t_524 = buffer.data(target + 524);
    auto *t_525 = buffer.data(target + 525);
    auto *t_526 = buffer.data(target + 526);
    auto *t_527 = buffer.data(target + 527);
    auto *t_528 = buffer.data(target + 528);
    auto *t_529 = buffer.data(target + 529);
    auto *t_530 = buffer.data(target + 530);
    auto *t_531 = buffer.data(target + 531);
    auto *t_532 = buffer.data(target + 532);
    auto *t_533 = buffer.data(target + 533);
    auto *t_534 = buffer.data(target + 534);
    auto *t_535 = buffer.data(target + 535);
    auto *t_536 = buffer.data(target + 536);
    auto *t_537 = buffer.data(target + 537);
    auto *t_538 = buffer.data(target + 538);
    auto *t_539 = buffer.data(target + 539);
    auto *t_540 = buffer.data(target + 540);
    auto *t_541 = buffer.data(target + 541);
    auto *t_542 = buffer.data(target + 542);
    auto *t_543 = buffer.data(target + 543);
    auto *t_544 = buffer.data(target + 544);
    auto *t_545 = buffer.data(target + 545);
    auto *t_546 = buffer.data(target + 546);
    auto *t_547 = buffer.data(target + 547);
    auto *t_548 = buffer.data(target + 548);
    auto *t_549 = buffer.data(target + 549);
    auto *t_550 = buffer.data(target + 550);
    auto *t_551 = buffer.data(target + 551);
    auto *t_552 = buffer.data(target + 552);
    auto *t_553 = buffer.data(target + 553);
    auto *t_554 = buffer.data(target + 554);
    auto *t_555 = buffer.data(target + 555);
    auto *t_556 = buffer.data(target + 556);
    auto *t_557 = buffer.data(target + 557);
    auto *t_558 = buffer.data(target + 558);
    auto *t_559 = buffer.data(target + 559);
    auto *t_560 = buffer.data(target + 560);
    auto *t_561 = buffer.data(target + 561);
    auto *t_562 = buffer.data(target + 562);
    auto *t_563 = buffer.data(target + 563);
    auto *t_564 = buffer.data(target + 564);
    auto *t_565 = buffer.data(target + 565);
    auto *t_566 = buffer.data(target + 566);
    auto *t_567 = buffer.data(target + 567);
    auto *t_568 = buffer.data(target + 568);
    auto *t_569 = buffer.data(target + 569);
    auto *t_570 = buffer.data(target + 570);
    auto *t_571 = buffer.data(target + 571);
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);
    auto *t_575 = buffer.data(target + 575);
    auto *t_576 = buffer.data(target + 576);
    auto *t_577 = buffer.data(target + 577);
    auto *t_578 = buffer.data(target + 578);
    auto *t_579 = buffer.data(target + 579);
    auto *t_580 = buffer.data(target + 580);
    auto *t_581 = buffer.data(target + 581);
    auto *t_582 = buffer.data(target + 582);
    auto *t_583 = buffer.data(target + 583);
    auto *t_584 = buffer.data(target + 584);
    auto *t_585 = buffer.data(target + 585);
    auto *t_586 = buffer.data(target + 586);
    auto *t_587 = buffer.data(target + 587);
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);
    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);
    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);
    auto *t_618 = buffer.data(target + 618);
    auto *t_619 = buffer.data(target + 619);
    auto *t_620 = buffer.data(target + 620);
    auto *t_621 = buffer.data(target + 621);
    auto *t_622 = buffer.data(target + 622);
    auto *t_623 = buffer.data(target + 623);
    auto *t_624 = buffer.data(target + 624);
    auto *t_625 = buffer.data(target + 625);
    auto *t_626 = buffer.data(target + 626);
    auto *t_627 = buffer.data(target + 627);
    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);
    auto *t_637 = buffer.data(target + 637);
    auto *t_638 = buffer.data(target + 638);
    auto *t_639 = buffer.data(target + 639);
    auto *t_640 = buffer.data(target + 640);
    auto *t_641 = buffer.data(target + 641);
    auto *t_642 = buffer.data(target + 642);
    auto *t_643 = buffer.data(target + 643);
    auto *t_644 = buffer.data(target + 644);
    auto *t_645 = buffer.data(target + 645);
    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);

    const auto *kh_349 = buffer.data(kh + 349);
    const auto *kh_350 = buffer.data(kh + 350);
    const auto *kh_351 = buffer.data(kh + 351);
    const auto *kh_352 = buffer.data(kh + 352);
    const auto *kh_353 = buffer.data(kh + 353);
    const auto *kh_354 = buffer.data(kh + 354);
    const auto *kh_355 = buffer.data(kh + 355);
    const auto *kh_356 = buffer.data(kh + 356);
    const auto *kh_357 = buffer.data(kh + 357);
    const auto *kh_358 = buffer.data(kh + 358);
    const auto *kh_359 = buffer.data(kh + 359);
    const auto *kh_360 = buffer.data(kh + 360);
    const auto *kh_361 = buffer.data(kh + 361);
    const auto *kh_362 = buffer.data(kh + 362);
    const auto *kh_363 = buffer.data(kh + 363);
    const auto *kh_364 = buffer.data(kh + 364);
    const auto *kh_365 = buffer.data(kh + 365);
    const auto *kh_366 = buffer.data(kh + 366);
    const auto *kh_367 = buffer.data(kh + 367);
    const auto *kh_368 = buffer.data(kh + 368);
    const auto *kh_369 = buffer.data(kh + 369);
    const auto *kh_370 = buffer.data(kh + 370);
    const auto *kh_371 = buffer.data(kh + 371);
    const auto *kh_372 = buffer.data(kh + 372);
    const auto *kh_373 = buffer.data(kh + 373);
    const auto *kh_374 = buffer.data(kh + 374);
    const auto *kh_375 = buffer.data(kh + 375);
    const auto *kh_376 = buffer.data(kh + 376);
    const auto *kh_377 = buffer.data(kh + 377);
    const auto *kh_378 = buffer.data(kh + 378);
    const auto *kh_379 = buffer.data(kh + 379);
    const auto *kh_380 = buffer.data(kh + 380);
    const auto *kh_381 = buffer.data(kh + 381);
    const auto *kh_382 = buffer.data(kh + 382);
    const auto *kh_383 = buffer.data(kh + 383);
    const auto *kh_384 = buffer.data(kh + 384);
    const auto *kh_385 = buffer.data(kh + 385);
    const auto *kh_386 = buffer.data(kh + 386);
    const auto *kh_387 = buffer.data(kh + 387);
    const auto *kh_388 = buffer.data(kh + 388);
    const auto *kh_389 = buffer.data(kh + 389);
    const auto *kh_390 = buffer.data(kh + 390);
    const auto *kh_391 = buffer.data(kh + 391);
    const auto *kh_392 = buffer.data(kh + 392);
    const auto *kh_393 = buffer.data(kh + 393);
    const auto *kh_394 = buffer.data(kh + 394);
    const auto *kh_395 = buffer.data(kh + 395);
    const auto *kh_396 = buffer.data(kh + 396);
    const auto *kh_397 = buffer.data(kh + 397);
    const auto *kh_398 = buffer.data(kh + 398);
    const auto *kh_399 = buffer.data(kh + 399);
    const auto *kh_400 = buffer.data(kh + 400);
    const auto *kh_401 = buffer.data(kh + 401);
    const auto *kh_402 = buffer.data(kh + 402);
    const auto *kh_403 = buffer.data(kh + 403);
    const auto *kh_404 = buffer.data(kh + 404);
    const auto *kh_405 = buffer.data(kh + 405);
    const auto *kh_406 = buffer.data(kh + 406);
    const auto *kh_407 = buffer.data(kh + 407);
    const auto *kh_408 = buffer.data(kh + 408);
    const auto *kh_409 = buffer.data(kh + 409);
    const auto *kh_410 = buffer.data(kh + 410);
    const auto *kh_411 = buffer.data(kh + 411);
    const auto *kh_412 = buffer.data(kh + 412);
    const auto *kh_413 = buffer.data(kh + 413);
    const auto *kh_414 = buffer.data(kh + 414);
    const auto *kh_415 = buffer.data(kh + 415);
    const auto *kh_416 = buffer.data(kh + 416);
    const auto *kh_417 = buffer.data(kh + 417);
    const auto *kh_418 = buffer.data(kh + 418);
    const auto *kh_419 = buffer.data(kh + 419);
    const auto *kh_420 = buffer.data(kh + 420);
    const auto *kh_421 = buffer.data(kh + 421);
    const auto *kh_422 = buffer.data(kh + 422);
    const auto *kh_423 = buffer.data(kh + 423);
    const auto *kh_424 = buffer.data(kh + 424);
    const auto *kh_425 = buffer.data(kh + 425);
    const auto *kh_426 = buffer.data(kh + 426);
    const auto *kh_427 = buffer.data(kh + 427);
    const auto *kh_428 = buffer.data(kh + 428);
    const auto *kh_429 = buffer.data(kh + 429);
    const auto *kh_430 = buffer.data(kh + 430);
    const auto *kh_431 = buffer.data(kh + 431);
    const auto *kh_432 = buffer.data(kh + 432);
    const auto *kh_433 = buffer.data(kh + 433);
    const auto *kh_434 = buffer.data(kh + 434);
    const auto *kh_435 = buffer.data(kh + 435);
    const auto *kh_436 = buffer.data(kh + 436);
    const auto *kh_437 = buffer.data(kh + 437);
    const auto *kh_438 = buffer.data(kh + 438);
    const auto *kh_439 = buffer.data(kh + 439);
    const auto *kh_440 = buffer.data(kh + 440);
    const auto *kh_441 = buffer.data(kh + 441);
    const auto *kh_442 = buffer.data(kh + 442);
    const auto *kh_443 = buffer.data(kh + 443);
    const auto *kh_444 = buffer.data(kh + 444);
    const auto *kh_445 = buffer.data(kh + 445);
    const auto *kh_446 = buffer.data(kh + 446);
    const auto *kh_447 = buffer.data(kh + 447);
    const auto *kh_448 = buffer.data(kh + 448);
    const auto *kh_449 = buffer.data(kh + 449);
    const auto *kh_450 = buffer.data(kh + 450);
    const auto *kh_451 = buffer.data(kh + 451);
    const auto *kh_452 = buffer.data(kh + 452);
    const auto *kh_453 = buffer.data(kh + 453);
    const auto *kh_454 = buffer.data(kh + 454);
    const auto *kh_455 = buffer.data(kh + 455);
    const auto *kh_456 = buffer.data(kh + 456);
    const auto *kh_457 = buffer.data(kh + 457);
    const auto *kh_458 = buffer.data(kh + 458);
    const auto *kh_459 = buffer.data(kh + 459);
    const auto *kh_460 = buffer.data(kh + 460);
    const auto *kh_461 = buffer.data(kh + 461);
    const auto *kh_462 = buffer.data(kh + 462);
    const auto *kh_463 = buffer.data(kh + 463);
    const auto *kh_464 = buffer.data(kh + 464);
    const auto *kh_465 = buffer.data(kh + 465);
    const auto *kh_466 = buffer.data(kh + 466);
    const auto *kh_467 = buffer.data(kh + 467);
    const auto *kh_468 = buffer.data(kh + 468);
    const auto *kh_469 = buffer.data(kh + 469);
    const auto *kh_470 = buffer.data(kh + 470);
    const auto *kh_471 = buffer.data(kh + 471);
    const auto *kh_472 = buffer.data(kh + 472);
    const auto *kh_473 = buffer.data(kh + 473);
    const auto *kh_474 = buffer.data(kh + 474);
    const auto *kh_475 = buffer.data(kh + 475);
    const auto *kh_476 = buffer.data(kh + 476);
    const auto *kh_477 = buffer.data(kh + 477);
    const auto *kh_478 = buffer.data(kh + 478);
    const auto *kh_479 = buffer.data(kh + 479);
    const auto *kh_480 = buffer.data(kh + 480);
    const auto *kh_481 = buffer.data(kh + 481);
    const auto *kh_482 = buffer.data(kh + 482);
    const auto *kh_483 = buffer.data(kh + 483);
    const auto *kh_484 = buffer.data(kh + 484);
    const auto *kh_485 = buffer.data(kh + 485);

    const auto *mh_664 = buffer.data(mh + 664);
    const auto *mh_665 = buffer.data(mh + 665);
    const auto *mh_666 = buffer.data(mh + 666);
    const auto *mh_667 = buffer.data(mh + 667);
    const auto *mh_668 = buffer.data(mh + 668);
    const auto *mh_669 = buffer.data(mh + 669);
    const auto *mh_670 = buffer.data(mh + 670);
    const auto *mh_671 = buffer.data(mh + 671);
    const auto *mh_672 = buffer.data(mh + 672);
    const auto *mh_673 = buffer.data(mh + 673);
    const auto *mh_674 = buffer.data(mh + 674);
    const auto *mh_675 = buffer.data(mh + 675);
    const auto *mh_676 = buffer.data(mh + 676);
    const auto *mh_677 = buffer.data(mh + 677);
    const auto *mh_678 = buffer.data(mh + 678);
    const auto *mh_679 = buffer.data(mh + 679);
    const auto *mh_680 = buffer.data(mh + 680);
    const auto *mh_681 = buffer.data(mh + 681);
    const auto *mh_682 = buffer.data(mh + 682);
    const auto *mh_683 = buffer.data(mh + 683);
    const auto *mh_684 = buffer.data(mh + 684);
    const auto *mh_685 = buffer.data(mh + 685);
    const auto *mh_686 = buffer.data(mh + 686);
    const auto *mh_687 = buffer.data(mh + 687);
    const auto *mh_688 = buffer.data(mh + 688);
    const auto *mh_689 = buffer.data(mh + 689);
    const auto *mh_690 = buffer.data(mh + 690);
    const auto *mh_691 = buffer.data(mh + 691);
    const auto *mh_692 = buffer.data(mh + 692);
    const auto *mh_693 = buffer.data(mh + 693);
    const auto *mh_694 = buffer.data(mh + 694);
    const auto *mh_695 = buffer.data(mh + 695);
    const auto *mh_696 = buffer.data(mh + 696);
    const auto *mh_697 = buffer.data(mh + 697);
    const auto *mh_698 = buffer.data(mh + 698);
    const auto *mh_699 = buffer.data(mh + 699);
    const auto *mh_700 = buffer.data(mh + 700);
    const auto *mh_701 = buffer.data(mh + 701);
    const auto *mh_702 = buffer.data(mh + 702);
    const auto *mh_703 = buffer.data(mh + 703);
    const auto *mh_704 = buffer.data(mh + 704);
    const auto *mh_705 = buffer.data(mh + 705);
    const auto *mh_706 = buffer.data(mh + 706);
    const auto *mh_707 = buffer.data(mh + 707);
    const auto *mh_708 = buffer.data(mh + 708);
    const auto *mh_709 = buffer.data(mh + 709);
    const auto *mh_710 = buffer.data(mh + 710);
    const auto *mh_711 = buffer.data(mh + 711);
    const auto *mh_712 = buffer.data(mh + 712);
    const auto *mh_713 = buffer.data(mh + 713);
    const auto *mh_714 = buffer.data(mh + 714);
    const auto *mh_715 = buffer.data(mh + 715);
    const auto *mh_716 = buffer.data(mh + 716);
    const auto *mh_717 = buffer.data(mh + 717);
    const auto *mh_718 = buffer.data(mh + 718);
    const auto *mh_719 = buffer.data(mh + 719);
    const auto *mh_720 = buffer.data(mh + 720);
    const auto *mh_721 = buffer.data(mh + 721);
    const auto *mh_722 = buffer.data(mh + 722);
    const auto *mh_723 = buffer.data(mh + 723);
    const auto *mh_724 = buffer.data(mh + 724);
    const auto *mh_725 = buffer.data(mh + 725);
    const auto *mh_726 = buffer.data(mh + 726);
    const auto *mh_727 = buffer.data(mh + 727);
    const auto *mh_728 = buffer.data(mh + 728);
    const auto *mh_729 = buffer.data(mh + 729);
    const auto *mh_730 = buffer.data(mh + 730);
    const auto *mh_731 = buffer.data(mh + 731);
    const auto *mh_732 = buffer.data(mh + 732);
    const auto *mh_733 = buffer.data(mh + 733);
    const auto *mh_734 = buffer.data(mh + 734);
    const auto *mh_735 = buffer.data(mh + 735);
    const auto *mh_736 = buffer.data(mh + 736);
    const auto *mh_737 = buffer.data(mh + 737);
    const auto *mh_738 = buffer.data(mh + 738);
    const auto *mh_739 = buffer.data(mh + 739);
    const auto *mh_740 = buffer.data(mh + 740);
    const auto *mh_741 = buffer.data(mh + 741);
    const auto *mh_742 = buffer.data(mh + 742);
    const auto *mh_743 = buffer.data(mh + 743);
    const auto *mh_744 = buffer.data(mh + 744);
    const auto *mh_745 = buffer.data(mh + 745);
    const auto *mh_746 = buffer.data(mh + 746);
    const auto *mh_747 = buffer.data(mh + 747);
    const auto *mh_748 = buffer.data(mh + 748);
    const auto *mh_749 = buffer.data(mh + 749);
    const auto *mh_750 = buffer.data(mh + 750);
    const auto *mh_751 = buffer.data(mh + 751);
    const auto *mh_752 = buffer.data(mh + 752);
    const auto *mh_753 = buffer.data(mh + 753);
    const auto *mh_754 = buffer.data(mh + 754);
    const auto *mh_755 = buffer.data(mh + 755);
    const auto *mh_777 = buffer.data(mh + 777);
    const auto *mh_778 = buffer.data(mh + 778);
    const auto *mh_779 = buffer.data(mh + 779);
    const auto *mh_780 = buffer.data(mh + 780);
    const auto *mh_781 = buffer.data(mh + 781);
    const auto *mh_782 = buffer.data(mh + 782);
    const auto *mh_783 = buffer.data(mh + 783);
    const auto *mh_784 = buffer.data(mh + 784);
    const auto *mh_785 = buffer.data(mh + 785);
    const auto *mh_786 = buffer.data(mh + 786);
    const auto *mh_787 = buffer.data(mh + 787);
    const auto *mh_788 = buffer.data(mh + 788);
    const auto *mh_789 = buffer.data(mh + 789);
    const auto *mh_790 = buffer.data(mh + 790);
    const auto *mh_791 = buffer.data(mh + 791);
    const auto *mh_792 = buffer.data(mh + 792);
    const auto *mh_793 = buffer.data(mh + 793);
    const auto *mh_794 = buffer.data(mh + 794);
    const auto *mh_795 = buffer.data(mh + 795);
    const auto *mh_796 = buffer.data(mh + 796);
    const auto *mh_797 = buffer.data(mh + 797);
    const auto *mh_798 = buffer.data(mh + 798);
    const auto *mh_799 = buffer.data(mh + 799);
    const auto *mh_800 = buffer.data(mh + 800);
    const auto *mh_801 = buffer.data(mh + 801);
    const auto *mh_802 = buffer.data(mh + 802);
    const auto *mh_803 = buffer.data(mh + 803);
    const auto *mh_804 = buffer.data(mh + 804);
    const auto *mh_805 = buffer.data(mh + 805);
    const auto *mh_806 = buffer.data(mh + 806);
    const auto *mh_807 = buffer.data(mh + 807);
    const auto *mh_808 = buffer.data(mh + 808);
    const auto *mh_809 = buffer.data(mh + 809);
    const auto *mh_810 = buffer.data(mh + 810);
    const auto *mh_811 = buffer.data(mh + 811);
    const auto *mh_812 = buffer.data(mh + 812);
    const auto *mh_813 = buffer.data(mh + 813);
    const auto *mh_814 = buffer.data(mh + 814);
    const auto *mh_815 = buffer.data(mh + 815);
    const auto *mh_816 = buffer.data(mh + 816);
    const auto *mh_817 = buffer.data(mh + 817);
    const auto *mh_818 = buffer.data(mh + 818);
    const auto *mh_819 = buffer.data(mh + 819);
    const auto *mh_820 = buffer.data(mh + 820);
    const auto *mh_821 = buffer.data(mh + 821);
    const auto *mh_822 = buffer.data(mh + 822);
    const auto *mh_823 = buffer.data(mh + 823);
    const auto *mh_824 = buffer.data(mh + 824);
    const auto *mh_825 = buffer.data(mh + 825);
    const auto *mh_826 = buffer.data(mh + 826);
    const auto *mh_827 = buffer.data(mh + 827);
    const auto *mh_828 = buffer.data(mh + 828);
    const auto *mh_829 = buffer.data(mh + 829);
    const auto *mh_830 = buffer.data(mh + 830);
    const auto *mh_831 = buffer.data(mh + 831);
    const auto *mh_832 = buffer.data(mh + 832);
    const auto *mh_833 = buffer.data(mh + 833);
    const auto *mh_834 = buffer.data(mh + 834);
    const auto *mh_835 = buffer.data(mh + 835);
    const auto *mh_836 = buffer.data(mh + 836);
    const auto *mh_837 = buffer.data(mh + 837);
    const auto *mh_838 = buffer.data(mh + 838);
    const auto *mh_839 = buffer.data(mh + 839);
    const auto *mh_840 = buffer.data(mh + 840);
    const auto *mh_841 = buffer.data(mh + 841);
    const auto *mh_842 = buffer.data(mh + 842);

#pragma omp simd aligned(t_496, t_497, t_498, t_499, t_500, kh_349, kh_350, kh_351, kh_352, \
                         kh_353, mh_664, mh_665, mh_666, mh_667, \
                         mh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = -2.0 * kh_349[k]
                   + f_0 * mh_664[k];

        t_497[k] = -2.0 * kh_350[k]
                   + f_0 * mh_665[k];

        t_498[k] = -2.0 * kh_351[k]
                   + f_0 * mh_666[k];

        t_499[k] = -2.0 * kh_352[k]
                   + f_0 * mh_667[k];

        t_500[k] = -2.0 * kh_353[k]
                   + f_0 * mh_668[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, t_505, kh_354, kh_355, kh_356, kh_357, \
                         kh_358, mh_669, mh_670, mh_671, mh_672, \
                         mh_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = -2.0 * kh_354[k]
                   + f_0 * mh_669[k];

        t_502[k] = -2.0 * kh_355[k]
                   + f_0 * mh_670[k];

        t_503[k] = -2.0 * kh_356[k]
                   + f_0 * mh_671[k];

        t_504[k] = -3.0 * kh_357[k]
                   + f_0 * mh_672[k];

        t_505[k] = -3.0 * kh_358[k]
                   + f_0 * mh_673[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, t_509, t_510, kh_359, kh_360, kh_361, kh_362, \
                         kh_363, mh_674, mh_675, mh_676, mh_677, \
                         mh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = -3.0 * kh_359[k]
                   + f_0 * mh_674[k];

        t_507[k] = -3.0 * kh_360[k]
                   + f_0 * mh_675[k];

        t_508[k] = -3.0 * kh_361[k]
                   + f_0 * mh_676[k];

        t_509[k] = -3.0 * kh_362[k]
                   + f_0 * mh_677[k];

        t_510[k] = -3.0 * kh_363[k]
                   + f_0 * mh_678[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, t_515, kh_364, kh_365, kh_366, kh_367, \
                         kh_368, mh_679, mh_680, mh_681, mh_682, \
                         mh_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = -3.0 * kh_364[k]
                   + f_0 * mh_679[k];

        t_512[k] = -3.0 * kh_365[k]
                   + f_0 * mh_680[k];

        t_513[k] = -3.0 * kh_366[k]
                   + f_0 * mh_681[k];

        t_514[k] = -3.0 * kh_367[k]
                   + f_0 * mh_682[k];

        t_515[k] = -3.0 * kh_368[k]
                   + f_0 * mh_683[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, t_520, kh_369, kh_370, kh_371, kh_372, \
                         kh_373, mh_684, mh_685, mh_686, mh_687, \
                         mh_688 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = -3.0 * kh_369[k]
                   + f_0 * mh_684[k];

        t_517[k] = -3.0 * kh_370[k]
                   + f_0 * mh_685[k];

        t_518[k] = -3.0 * kh_371[k]
                   + f_0 * mh_686[k];

        t_519[k] = -3.0 * kh_372[k]
                   + f_0 * mh_687[k];

        t_520[k] = -3.0 * kh_373[k]
                   + f_0 * mh_688[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, t_525, kh_374, kh_375, kh_376, kh_377, \
                         kh_378, mh_689, mh_690, mh_691, mh_692, \
                         mh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = -3.0 * kh_374[k]
                   + f_0 * mh_689[k];

        t_522[k] = -3.0 * kh_375[k]
                   + f_0 * mh_690[k];

        t_523[k] = -3.0 * kh_376[k]
                   + f_0 * mh_691[k];

        t_524[k] = -3.0 * kh_377[k]
                   + f_0 * mh_692[k];

        t_525[k] = -4.0 * kh_378[k]
                   + f_0 * mh_693[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, kh_379, kh_380, kh_381, kh_382, \
                         kh_383, mh_694, mh_695, mh_696, mh_697, \
                         mh_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = -4.0 * kh_379[k]
                   + f_0 * mh_694[k];

        t_527[k] = -4.0 * kh_380[k]
                   + f_0 * mh_695[k];

        t_528[k] = -4.0 * kh_381[k]
                   + f_0 * mh_696[k];

        t_529[k] = -4.0 * kh_382[k]
                   + f_0 * mh_697[k];

        t_530[k] = -4.0 * kh_383[k]
                   + f_0 * mh_698[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, t_535, kh_384, kh_385, kh_386, kh_387, \
                         kh_388, mh_699, mh_700, mh_701, mh_702, \
                         mh_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = -4.0 * kh_384[k]
                   + f_0 * mh_699[k];

        t_532[k] = -4.0 * kh_385[k]
                   + f_0 * mh_700[k];

        t_533[k] = -4.0 * kh_386[k]
                   + f_0 * mh_701[k];

        t_534[k] = -4.0 * kh_387[k]
                   + f_0 * mh_702[k];

        t_535[k] = -4.0 * kh_388[k]
                   + f_0 * mh_703[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, kh_389, kh_390, kh_391, kh_392, \
                         kh_393, mh_704, mh_705, mh_706, mh_707, \
                         mh_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = -4.0 * kh_389[k]
                   + f_0 * mh_704[k];

        t_537[k] = -4.0 * kh_390[k]
                   + f_0 * mh_705[k];

        t_538[k] = -4.0 * kh_391[k]
                   + f_0 * mh_706[k];

        t_539[k] = -4.0 * kh_392[k]
                   + f_0 * mh_707[k];

        t_540[k] = -4.0 * kh_393[k]
                   + f_0 * mh_708[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, t_545, kh_394, kh_395, kh_396, kh_397, \
                         kh_398, mh_709, mh_710, mh_711, mh_712, \
                         mh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = -4.0 * kh_394[k]
                   + f_0 * mh_709[k];

        t_542[k] = -4.0 * kh_395[k]
                   + f_0 * mh_710[k];

        t_543[k] = -4.0 * kh_396[k]
                   + f_0 * mh_711[k];

        t_544[k] = -4.0 * kh_397[k]
                   + f_0 * mh_712[k];

        t_545[k] = -4.0 * kh_398[k]
                   + f_0 * mh_713[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, t_550, kh_399, kh_400, kh_401, kh_402, \
                         kh_403, mh_714, mh_715, mh_716, mh_717, \
                         mh_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = -5.0 * kh_399[k]
                   + f_0 * mh_714[k];

        t_547[k] = -5.0 * kh_400[k]
                   + f_0 * mh_715[k];

        t_548[k] = -5.0 * kh_401[k]
                   + f_0 * mh_716[k];

        t_549[k] = -5.0 * kh_402[k]
                   + f_0 * mh_717[k];

        t_550[k] = -5.0 * kh_403[k]
                   + f_0 * mh_718[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, t_554, t_555, kh_404, kh_405, kh_406, kh_407, \
                         kh_408, mh_719, mh_720, mh_721, mh_722, \
                         mh_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = -5.0 * kh_404[k]
                   + f_0 * mh_719[k];

        t_552[k] = -5.0 * kh_405[k]
                   + f_0 * mh_720[k];

        t_553[k] = -5.0 * kh_406[k]
                   + f_0 * mh_721[k];

        t_554[k] = -5.0 * kh_407[k]
                   + f_0 * mh_722[k];

        t_555[k] = -5.0 * kh_408[k]
                   + f_0 * mh_723[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, kh_409, kh_410, kh_411, kh_412, \
                         kh_413, mh_724, mh_725, mh_726, mh_727, \
                         mh_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = -5.0 * kh_409[k]
                   + f_0 * mh_724[k];

        t_557[k] = -5.0 * kh_410[k]
                   + f_0 * mh_725[k];

        t_558[k] = -5.0 * kh_411[k]
                   + f_0 * mh_726[k];

        t_559[k] = -5.0 * kh_412[k]
                   + f_0 * mh_727[k];

        t_560[k] = -5.0 * kh_413[k]
                   + f_0 * mh_728[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, t_565, kh_414, kh_415, kh_416, kh_417, \
                         kh_418, mh_729, mh_730, mh_731, mh_732, \
                         mh_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = -5.0 * kh_414[k]
                   + f_0 * mh_729[k];

        t_562[k] = -5.0 * kh_415[k]
                   + f_0 * mh_730[k];

        t_563[k] = -5.0 * kh_416[k]
                   + f_0 * mh_731[k];

        t_564[k] = -5.0 * kh_417[k]
                   + f_0 * mh_732[k];

        t_565[k] = -5.0 * kh_418[k]
                   + f_0 * mh_733[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, t_570, kh_419, kh_420, kh_421, kh_422, \
                         kh_423, mh_734, mh_735, mh_736, mh_737, \
                         mh_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = -5.0 * kh_419[k]
                   + f_0 * mh_734[k];

        t_567[k] = -6.0 * kh_420[k]
                   + f_0 * mh_735[k];

        t_568[k] = -6.0 * kh_421[k]
                   + f_0 * mh_736[k];

        t_569[k] = -6.0 * kh_422[k]
                   + f_0 * mh_737[k];

        t_570[k] = -6.0 * kh_423[k]
                   + f_0 * mh_738[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, t_574, t_575, kh_424, kh_425, kh_426, kh_427, \
                         kh_428, mh_739, mh_740, mh_741, mh_742, \
                         mh_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = -6.0 * kh_424[k]
                   + f_0 * mh_739[k];

        t_572[k] = -6.0 * kh_425[k]
                   + f_0 * mh_740[k];

        t_573[k] = -6.0 * kh_426[k]
                   + f_0 * mh_741[k];

        t_574[k] = -6.0 * kh_427[k]
                   + f_0 * mh_742[k];

        t_575[k] = -6.0 * kh_428[k]
                   + f_0 * mh_743[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, kh_429, kh_430, kh_431, kh_432, \
                         kh_433, mh_744, mh_745, mh_746, mh_747, \
                         mh_748 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = -6.0 * kh_429[k]
                   + f_0 * mh_744[k];

        t_577[k] = -6.0 * kh_430[k]
                   + f_0 * mh_745[k];

        t_578[k] = -6.0 * kh_431[k]
                   + f_0 * mh_746[k];

        t_579[k] = -6.0 * kh_432[k]
                   + f_0 * mh_747[k];

        t_580[k] = -6.0 * kh_433[k]
                   + f_0 * mh_748[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, t_584, t_585, kh_434, kh_435, kh_436, kh_437, \
                         kh_438, mh_749, mh_750, mh_751, mh_752, \
                         mh_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = -6.0 * kh_434[k]
                   + f_0 * mh_749[k];

        t_582[k] = -6.0 * kh_435[k]
                   + f_0 * mh_750[k];

        t_583[k] = -6.0 * kh_436[k]
                   + f_0 * mh_751[k];

        t_584[k] = -6.0 * kh_437[k]
                   + f_0 * mh_752[k];

        t_585[k] = -6.0 * kh_438[k]
                   + f_0 * mh_753[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, t_590, t_591, t_592, kh_439, kh_440, \
                         mh_754, mh_755, mh_777, mh_778, mh_779, mh_780, \
                         mh_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = -6.0 * kh_439[k]
                   + f_0 * mh_754[k];

        t_587[k] = -6.0 * kh_440[k]
                   + f_0 * mh_755[k];

        t_588[k] = f_0 * mh_777[k];

        t_589[k] = f_0 * mh_778[k];

        t_590[k] = f_0 * mh_779[k];

        t_591[k] = f_0 * mh_780[k];

        t_592[k] = f_0 * mh_781[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, t_596, t_597, t_598, t_599, t_600, mh_782, \
                         mh_783, mh_784, mh_785, mh_786, mh_787, mh_788, \
                         mh_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = f_0 * mh_782[k];

        t_594[k] = f_0 * mh_783[k];

        t_595[k] = f_0 * mh_784[k];

        t_596[k] = f_0 * mh_785[k];

        t_597[k] = f_0 * mh_786[k];

        t_598[k] = f_0 * mh_787[k];

        t_599[k] = f_0 * mh_788[k];

        t_600[k] = f_0 * mh_789[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, t_604, t_605, t_606, t_607, t_608, mh_790, \
                         mh_791, mh_792, mh_793, mh_794, mh_795, mh_796, \
                         mh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_0 * mh_790[k];

        t_602[k] = f_0 * mh_791[k];

        t_603[k] = f_0 * mh_792[k];

        t_604[k] = f_0 * mh_793[k];

        t_605[k] = f_0 * mh_794[k];

        t_606[k] = f_0 * mh_795[k];

        t_607[k] = f_0 * mh_796[k];

        t_608[k] = f_0 * mh_797[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, t_613, kh_441, kh_442, kh_443, kh_444, \
                         kh_445, mh_798, mh_799, mh_800, mh_801, \
                         mh_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = -kh_441[k]
                   + f_0 * mh_798[k];

        t_610[k] = -kh_442[k]
                   + f_0 * mh_799[k];

        t_611[k] = -kh_443[k]
                   + f_0 * mh_800[k];

        t_612[k] = -kh_444[k]
                   + f_0 * mh_801[k];

        t_613[k] = -kh_445[k]
                   + f_0 * mh_802[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, kh_446, kh_447, kh_448, kh_449, \
                         kh_450, mh_803, mh_804, mh_805, mh_806, \
                         mh_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = -kh_446[k]
                   + f_0 * mh_803[k];

        t_615[k] = -kh_447[k]
                   + f_0 * mh_804[k];

        t_616[k] = -kh_448[k]
                   + f_0 * mh_805[k];

        t_617[k] = -kh_449[k]
                   + f_0 * mh_806[k];

        t_618[k] = -kh_450[k]
                   + f_0 * mh_807[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, kh_451, kh_452, kh_453, kh_454, \
                         kh_455, mh_808, mh_809, mh_810, mh_811, \
                         mh_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = -kh_451[k]
                   + f_0 * mh_808[k];

        t_620[k] = -kh_452[k]
                   + f_0 * mh_809[k];

        t_621[k] = -kh_453[k]
                   + f_0 * mh_810[k];

        t_622[k] = -kh_454[k]
                   + f_0 * mh_811[k];

        t_623[k] = -kh_455[k]
                   + f_0 * mh_812[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, t_628, kh_456, kh_457, kh_458, kh_459, \
                         kh_460, mh_813, mh_814, mh_815, mh_816, \
                         mh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = -kh_456[k]
                   + f_0 * mh_813[k];

        t_625[k] = -kh_457[k]
                   + f_0 * mh_814[k];

        t_626[k] = -kh_458[k]
                   + f_0 * mh_815[k];

        t_627[k] = -kh_459[k]
                   + f_0 * mh_816[k];

        t_628[k] = -kh_460[k]
                   + f_0 * mh_817[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, t_633, kh_461, kh_462, kh_463, kh_464, \
                         kh_465, mh_818, mh_819, mh_820, mh_821, \
                         mh_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = -kh_461[k]
                   + f_0 * mh_818[k];

        t_630[k] = -2.0 * kh_462[k]
                   + f_0 * mh_819[k];

        t_631[k] = -2.0 * kh_463[k]
                   + f_0 * mh_820[k];

        t_632[k] = -2.0 * kh_464[k]
                   + f_0 * mh_821[k];

        t_633[k] = -2.0 * kh_465[k]
                   + f_0 * mh_822[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, kh_466, kh_467, kh_468, kh_469, \
                         kh_470, mh_823, mh_824, mh_825, mh_826, \
                         mh_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = -2.0 * kh_466[k]
                   + f_0 * mh_823[k];

        t_635[k] = -2.0 * kh_467[k]
                   + f_0 * mh_824[k];

        t_636[k] = -2.0 * kh_468[k]
                   + f_0 * mh_825[k];

        t_637[k] = -2.0 * kh_469[k]
                   + f_0 * mh_826[k];

        t_638[k] = -2.0 * kh_470[k]
                   + f_0 * mh_827[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, t_643, kh_471, kh_472, kh_473, kh_474, \
                         kh_475, mh_828, mh_829, mh_830, mh_831, \
                         mh_832 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = -2.0 * kh_471[k]
                   + f_0 * mh_828[k];

        t_640[k] = -2.0 * kh_472[k]
                   + f_0 * mh_829[k];

        t_641[k] = -2.0 * kh_473[k]
                   + f_0 * mh_830[k];

        t_642[k] = -2.0 * kh_474[k]
                   + f_0 * mh_831[k];

        t_643[k] = -2.0 * kh_475[k]
                   + f_0 * mh_832[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, t_648, kh_476, kh_477, kh_478, kh_479, \
                         kh_480, mh_833, mh_834, mh_835, mh_836, \
                         mh_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = -2.0 * kh_476[k]
                   + f_0 * mh_833[k];

        t_645[k] = -2.0 * kh_477[k]
                   + f_0 * mh_834[k];

        t_646[k] = -2.0 * kh_478[k]
                   + f_0 * mh_835[k];

        t_647[k] = -2.0 * kh_479[k]
                   + f_0 * mh_836[k];

        t_648[k] = -2.0 * kh_480[k]
                   + f_0 * mh_837[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, kh_481, kh_482, kh_483, kh_484, \
                         kh_485, mh_838, mh_839, mh_840, mh_841, \
                         mh_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = -2.0 * kh_481[k]
                   + f_0 * mh_838[k];

        t_650[k] = -2.0 * kh_482[k]
                   + f_0 * mh_839[k];

        t_651[k] = -3.0 * kh_483[k]
                   + f_0 * mh_840[k];

        t_652[k] = -3.0 * kh_484[k]
                   + f_0 * mh_841[k];

        t_653[k] = -3.0 * kh_485[k]
                   + f_0 * mh_842[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_654 = buffer.data(target + 654);
    auto *t_655 = buffer.data(target + 655);
    auto *t_656 = buffer.data(target + 656);
    auto *t_657 = buffer.data(target + 657);
    auto *t_658 = buffer.data(target + 658);
    auto *t_659 = buffer.data(target + 659);
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);
    auto *t_663 = buffer.data(target + 663);
    auto *t_664 = buffer.data(target + 664);
    auto *t_665 = buffer.data(target + 665);
    auto *t_666 = buffer.data(target + 666);
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);
    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);
    auto *t_686 = buffer.data(target + 686);
    auto *t_687 = buffer.data(target + 687);
    auto *t_688 = buffer.data(target + 688);
    auto *t_689 = buffer.data(target + 689);
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);
    auto *t_693 = buffer.data(target + 693);
    auto *t_694 = buffer.data(target + 694);
    auto *t_695 = buffer.data(target + 695);
    auto *t_696 = buffer.data(target + 696);
    auto *t_697 = buffer.data(target + 697);
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);
    auto *t_708 = buffer.data(target + 708);
    auto *t_709 = buffer.data(target + 709);
    auto *t_710 = buffer.data(target + 710);
    auto *t_711 = buffer.data(target + 711);
    auto *t_712 = buffer.data(target + 712);
    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
    auto *t_722 = buffer.data(target + 722);
    auto *t_723 = buffer.data(target + 723);
    auto *t_724 = buffer.data(target + 724);
    auto *t_725 = buffer.data(target + 725);
    auto *t_726 = buffer.data(target + 726);
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);
    auto *t_732 = buffer.data(target + 732);
    auto *t_733 = buffer.data(target + 733);
    auto *t_734 = buffer.data(target + 734);
    auto *t_735 = buffer.data(target + 735);
    auto *t_736 = buffer.data(target + 736);
    auto *t_737 = buffer.data(target + 737);
    auto *t_738 = buffer.data(target + 738);
    auto *t_739 = buffer.data(target + 739);
    auto *t_740 = buffer.data(target + 740);
    auto *t_741 = buffer.data(target + 741);
    auto *t_742 = buffer.data(target + 742);
    auto *t_743 = buffer.data(target + 743);
    auto *t_744 = buffer.data(target + 744);
    auto *t_745 = buffer.data(target + 745);
    auto *t_746 = buffer.data(target + 746);
    auto *t_747 = buffer.data(target + 747);
    auto *t_748 = buffer.data(target + 748);
    auto *t_749 = buffer.data(target + 749);
    auto *t_750 = buffer.data(target + 750);
    auto *t_751 = buffer.data(target + 751);
    auto *t_752 = buffer.data(target + 752);
    auto *t_753 = buffer.data(target + 753);
    auto *t_754 = buffer.data(target + 754);
    auto *t_755 = buffer.data(target + 755);
    auto *t_756 = buffer.data(target + 756);
    auto *t_757 = buffer.data(target + 757);
    auto *t_758 = buffer.data(target + 758);
    auto *t_759 = buffer.data(target + 759);
    auto *t_760 = buffer.data(target + 760);
    auto *t_761 = buffer.data(target + 761);
    auto *t_762 = buffer.data(target + 762);
    auto *t_763 = buffer.data(target + 763);
    auto *t_764 = buffer.data(target + 764);
    auto *t_765 = buffer.data(target + 765);
    auto *t_766 = buffer.data(target + 766);
    auto *t_767 = buffer.data(target + 767);
    auto *t_768 = buffer.data(target + 768);
    auto *t_769 = buffer.data(target + 769);
    auto *t_770 = buffer.data(target + 770);
    auto *t_771 = buffer.data(target + 771);
    auto *t_772 = buffer.data(target + 772);
    auto *t_773 = buffer.data(target + 773);
    auto *t_774 = buffer.data(target + 774);
    auto *t_775 = buffer.data(target + 775);
    auto *t_776 = buffer.data(target + 776);
    auto *t_777 = buffer.data(target + 777);
    auto *t_778 = buffer.data(target + 778);
    auto *t_779 = buffer.data(target + 779);
    auto *t_780 = buffer.data(target + 780);
    auto *t_781 = buffer.data(target + 781);
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);
    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);
    auto *t_801 = buffer.data(target + 801);
    auto *t_802 = buffer.data(target + 802);
    auto *t_803 = buffer.data(target + 803);
    auto *t_804 = buffer.data(target + 804);
    auto *t_805 = buffer.data(target + 805);
    auto *t_806 = buffer.data(target + 806);
    auto *t_807 = buffer.data(target + 807);
    auto *t_808 = buffer.data(target + 808);
    auto *t_809 = buffer.data(target + 809);
    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);

    const auto *kh_486 = buffer.data(kh + 486);
    const auto *kh_487 = buffer.data(kh + 487);
    const auto *kh_488 = buffer.data(kh + 488);
    const auto *kh_489 = buffer.data(kh + 489);
    const auto *kh_490 = buffer.data(kh + 490);
    const auto *kh_491 = buffer.data(kh + 491);
    const auto *kh_492 = buffer.data(kh + 492);
    const auto *kh_493 = buffer.data(kh + 493);
    const auto *kh_494 = buffer.data(kh + 494);
    const auto *kh_495 = buffer.data(kh + 495);
    const auto *kh_496 = buffer.data(kh + 496);
    const auto *kh_497 = buffer.data(kh + 497);
    const auto *kh_498 = buffer.data(kh + 498);
    const auto *kh_499 = buffer.data(kh + 499);
    const auto *kh_500 = buffer.data(kh + 500);
    const auto *kh_501 = buffer.data(kh + 501);
    const auto *kh_502 = buffer.data(kh + 502);
    const auto *kh_503 = buffer.data(kh + 503);
    const auto *kh_504 = buffer.data(kh + 504);
    const auto *kh_505 = buffer.data(kh + 505);
    const auto *kh_506 = buffer.data(kh + 506);
    const auto *kh_507 = buffer.data(kh + 507);
    const auto *kh_508 = buffer.data(kh + 508);
    const auto *kh_509 = buffer.data(kh + 509);
    const auto *kh_510 = buffer.data(kh + 510);
    const auto *kh_511 = buffer.data(kh + 511);
    const auto *kh_512 = buffer.data(kh + 512);
    const auto *kh_513 = buffer.data(kh + 513);
    const auto *kh_514 = buffer.data(kh + 514);
    const auto *kh_515 = buffer.data(kh + 515);
    const auto *kh_516 = buffer.data(kh + 516);
    const auto *kh_517 = buffer.data(kh + 517);
    const auto *kh_518 = buffer.data(kh + 518);
    const auto *kh_519 = buffer.data(kh + 519);
    const auto *kh_520 = buffer.data(kh + 520);
    const auto *kh_521 = buffer.data(kh + 521);
    const auto *kh_522 = buffer.data(kh + 522);
    const auto *kh_523 = buffer.data(kh + 523);
    const auto *kh_524 = buffer.data(kh + 524);
    const auto *kh_525 = buffer.data(kh + 525);
    const auto *kh_526 = buffer.data(kh + 526);
    const auto *kh_527 = buffer.data(kh + 527);
    const auto *kh_528 = buffer.data(kh + 528);
    const auto *kh_529 = buffer.data(kh + 529);
    const auto *kh_530 = buffer.data(kh + 530);
    const auto *kh_531 = buffer.data(kh + 531);
    const auto *kh_532 = buffer.data(kh + 532);
    const auto *kh_533 = buffer.data(kh + 533);
    const auto *kh_534 = buffer.data(kh + 534);
    const auto *kh_535 = buffer.data(kh + 535);
    const auto *kh_536 = buffer.data(kh + 536);
    const auto *kh_537 = buffer.data(kh + 537);
    const auto *kh_538 = buffer.data(kh + 538);
    const auto *kh_539 = buffer.data(kh + 539);
    const auto *kh_540 = buffer.data(kh + 540);
    const auto *kh_541 = buffer.data(kh + 541);
    const auto *kh_542 = buffer.data(kh + 542);
    const auto *kh_543 = buffer.data(kh + 543);
    const auto *kh_544 = buffer.data(kh + 544);
    const auto *kh_545 = buffer.data(kh + 545);
    const auto *kh_546 = buffer.data(kh + 546);
    const auto *kh_547 = buffer.data(kh + 547);
    const auto *kh_548 = buffer.data(kh + 548);
    const auto *kh_549 = buffer.data(kh + 549);
    const auto *kh_550 = buffer.data(kh + 550);
    const auto *kh_551 = buffer.data(kh + 551);
    const auto *kh_552 = buffer.data(kh + 552);
    const auto *kh_553 = buffer.data(kh + 553);
    const auto *kh_554 = buffer.data(kh + 554);
    const auto *kh_555 = buffer.data(kh + 555);
    const auto *kh_556 = buffer.data(kh + 556);
    const auto *kh_557 = buffer.data(kh + 557);
    const auto *kh_558 = buffer.data(kh + 558);
    const auto *kh_559 = buffer.data(kh + 559);
    const auto *kh_560 = buffer.data(kh + 560);
    const auto *kh_561 = buffer.data(kh + 561);
    const auto *kh_562 = buffer.data(kh + 562);
    const auto *kh_563 = buffer.data(kh + 563);
    const auto *kh_564 = buffer.data(kh + 564);
    const auto *kh_565 = buffer.data(kh + 565);
    const auto *kh_566 = buffer.data(kh + 566);
    const auto *kh_567 = buffer.data(kh + 567);
    const auto *kh_568 = buffer.data(kh + 568);
    const auto *kh_569 = buffer.data(kh + 569);
    const auto *kh_570 = buffer.data(kh + 570);
    const auto *kh_571 = buffer.data(kh + 571);
    const auto *kh_572 = buffer.data(kh + 572);
    const auto *kh_573 = buffer.data(kh + 573);
    const auto *kh_574 = buffer.data(kh + 574);
    const auto *kh_575 = buffer.data(kh + 575);
    const auto *kh_576 = buffer.data(kh + 576);
    const auto *kh_577 = buffer.data(kh + 577);
    const auto *kh_578 = buffer.data(kh + 578);
    const auto *kh_579 = buffer.data(kh + 579);
    const auto *kh_580 = buffer.data(kh + 580);
    const auto *kh_581 = buffer.data(kh + 581);
    const auto *kh_582 = buffer.data(kh + 582);
    const auto *kh_583 = buffer.data(kh + 583);
    const auto *kh_584 = buffer.data(kh + 584);
    const auto *kh_585 = buffer.data(kh + 585);
    const auto *kh_586 = buffer.data(kh + 586);
    const auto *kh_587 = buffer.data(kh + 587);
    const auto *kh_588 = buffer.data(kh + 588);
    const auto *kh_589 = buffer.data(kh + 589);
    const auto *kh_590 = buffer.data(kh + 590);
    const auto *kh_591 = buffer.data(kh + 591);
    const auto *kh_592 = buffer.data(kh + 592);
    const auto *kh_593 = buffer.data(kh + 593);
    const auto *kh_594 = buffer.data(kh + 594);
    const auto *kh_595 = buffer.data(kh + 595);
    const auto *kh_596 = buffer.data(kh + 596);
    const auto *kh_597 = buffer.data(kh + 597);
    const auto *kh_598 = buffer.data(kh + 598);
    const auto *kh_599 = buffer.data(kh + 599);
    const auto *kh_600 = buffer.data(kh + 600);
    const auto *kh_601 = buffer.data(kh + 601);
    const auto *kh_602 = buffer.data(kh + 602);
    const auto *kh_603 = buffer.data(kh + 603);
    const auto *kh_604 = buffer.data(kh + 604);
    const auto *kh_605 = buffer.data(kh + 605);
    const auto *kh_606 = buffer.data(kh + 606);
    const auto *kh_607 = buffer.data(kh + 607);
    const auto *kh_608 = buffer.data(kh + 608);
    const auto *kh_609 = buffer.data(kh + 609);
    const auto *kh_610 = buffer.data(kh + 610);
    const auto *kh_611 = buffer.data(kh + 611);
    const auto *kh_612 = buffer.data(kh + 612);
    const auto *kh_613 = buffer.data(kh + 613);
    const auto *kh_614 = buffer.data(kh + 614);
    const auto *kh_615 = buffer.data(kh + 615);
    const auto *kh_616 = buffer.data(kh + 616);
    const auto *kh_617 = buffer.data(kh + 617);
    const auto *kh_618 = buffer.data(kh + 618);
    const auto *kh_619 = buffer.data(kh + 619);
    const auto *kh_620 = buffer.data(kh + 620);
    const auto *kh_621 = buffer.data(kh + 621);
    const auto *kh_622 = buffer.data(kh + 622);

    const auto *mh_843 = buffer.data(mh + 843);
    const auto *mh_844 = buffer.data(mh + 844);
    const auto *mh_845 = buffer.data(mh + 845);
    const auto *mh_846 = buffer.data(mh + 846);
    const auto *mh_847 = buffer.data(mh + 847);
    const auto *mh_848 = buffer.data(mh + 848);
    const auto *mh_849 = buffer.data(mh + 849);
    const auto *mh_850 = buffer.data(mh + 850);
    const auto *mh_851 = buffer.data(mh + 851);
    const auto *mh_852 = buffer.data(mh + 852);
    const auto *mh_853 = buffer.data(mh + 853);
    const auto *mh_854 = buffer.data(mh + 854);
    const auto *mh_855 = buffer.data(mh + 855);
    const auto *mh_856 = buffer.data(mh + 856);
    const auto *mh_857 = buffer.data(mh + 857);
    const auto *mh_858 = buffer.data(mh + 858);
    const auto *mh_859 = buffer.data(mh + 859);
    const auto *mh_860 = buffer.data(mh + 860);
    const auto *mh_861 = buffer.data(mh + 861);
    const auto *mh_862 = buffer.data(mh + 862);
    const auto *mh_863 = buffer.data(mh + 863);
    const auto *mh_864 = buffer.data(mh + 864);
    const auto *mh_865 = buffer.data(mh + 865);
    const auto *mh_866 = buffer.data(mh + 866);
    const auto *mh_867 = buffer.data(mh + 867);
    const auto *mh_868 = buffer.data(mh + 868);
    const auto *mh_869 = buffer.data(mh + 869);
    const auto *mh_870 = buffer.data(mh + 870);
    const auto *mh_871 = buffer.data(mh + 871);
    const auto *mh_872 = buffer.data(mh + 872);
    const auto *mh_873 = buffer.data(mh + 873);
    const auto *mh_874 = buffer.data(mh + 874);
    const auto *mh_875 = buffer.data(mh + 875);
    const auto *mh_876 = buffer.data(mh + 876);
    const auto *mh_877 = buffer.data(mh + 877);
    const auto *mh_878 = buffer.data(mh + 878);
    const auto *mh_879 = buffer.data(mh + 879);
    const auto *mh_880 = buffer.data(mh + 880);
    const auto *mh_881 = buffer.data(mh + 881);
    const auto *mh_882 = buffer.data(mh + 882);
    const auto *mh_883 = buffer.data(mh + 883);
    const auto *mh_884 = buffer.data(mh + 884);
    const auto *mh_885 = buffer.data(mh + 885);
    const auto *mh_886 = buffer.data(mh + 886);
    const auto *mh_887 = buffer.data(mh + 887);
    const auto *mh_888 = buffer.data(mh + 888);
    const auto *mh_889 = buffer.data(mh + 889);
    const auto *mh_890 = buffer.data(mh + 890);
    const auto *mh_891 = buffer.data(mh + 891);
    const auto *mh_892 = buffer.data(mh + 892);
    const auto *mh_893 = buffer.data(mh + 893);
    const auto *mh_894 = buffer.data(mh + 894);
    const auto *mh_895 = buffer.data(mh + 895);
    const auto *mh_896 = buffer.data(mh + 896);
    const auto *mh_897 = buffer.data(mh + 897);
    const auto *mh_898 = buffer.data(mh + 898);
    const auto *mh_899 = buffer.data(mh + 899);
    const auto *mh_900 = buffer.data(mh + 900);
    const auto *mh_901 = buffer.data(mh + 901);
    const auto *mh_902 = buffer.data(mh + 902);
    const auto *mh_903 = buffer.data(mh + 903);
    const auto *mh_904 = buffer.data(mh + 904);
    const auto *mh_905 = buffer.data(mh + 905);
    const auto *mh_906 = buffer.data(mh + 906);
    const auto *mh_907 = buffer.data(mh + 907);
    const auto *mh_908 = buffer.data(mh + 908);
    const auto *mh_909 = buffer.data(mh + 909);
    const auto *mh_910 = buffer.data(mh + 910);
    const auto *mh_911 = buffer.data(mh + 911);
    const auto *mh_912 = buffer.data(mh + 912);
    const auto *mh_913 = buffer.data(mh + 913);
    const auto *mh_914 = buffer.data(mh + 914);
    const auto *mh_915 = buffer.data(mh + 915);
    const auto *mh_916 = buffer.data(mh + 916);
    const auto *mh_917 = buffer.data(mh + 917);
    const auto *mh_918 = buffer.data(mh + 918);
    const auto *mh_919 = buffer.data(mh + 919);
    const auto *mh_920 = buffer.data(mh + 920);
    const auto *mh_921 = buffer.data(mh + 921);
    const auto *mh_922 = buffer.data(mh + 922);
    const auto *mh_923 = buffer.data(mh + 923);
    const auto *mh_924 = buffer.data(mh + 924);
    const auto *mh_925 = buffer.data(mh + 925);
    const auto *mh_926 = buffer.data(mh + 926);
    const auto *mh_927 = buffer.data(mh + 927);
    const auto *mh_928 = buffer.data(mh + 928);
    const auto *mh_929 = buffer.data(mh + 929);
    const auto *mh_930 = buffer.data(mh + 930);
    const auto *mh_931 = buffer.data(mh + 931);
    const auto *mh_932 = buffer.data(mh + 932);
    const auto *mh_933 = buffer.data(mh + 933);
    const auto *mh_934 = buffer.data(mh + 934);
    const auto *mh_935 = buffer.data(mh + 935);
    const auto *mh_936 = buffer.data(mh + 936);
    const auto *mh_937 = buffer.data(mh + 937);
    const auto *mh_938 = buffer.data(mh + 938);
    const auto *mh_939 = buffer.data(mh + 939);
    const auto *mh_940 = buffer.data(mh + 940);
    const auto *mh_941 = buffer.data(mh + 941);
    const auto *mh_942 = buffer.data(mh + 942);
    const auto *mh_943 = buffer.data(mh + 943);
    const auto *mh_944 = buffer.data(mh + 944);
    const auto *mh_966 = buffer.data(mh + 966);
    const auto *mh_967 = buffer.data(mh + 967);
    const auto *mh_968 = buffer.data(mh + 968);
    const auto *mh_969 = buffer.data(mh + 969);
    const auto *mh_970 = buffer.data(mh + 970);
    const auto *mh_971 = buffer.data(mh + 971);
    const auto *mh_972 = buffer.data(mh + 972);
    const auto *mh_973 = buffer.data(mh + 973);
    const auto *mh_974 = buffer.data(mh + 974);
    const auto *mh_975 = buffer.data(mh + 975);
    const auto *mh_976 = buffer.data(mh + 976);
    const auto *mh_977 = buffer.data(mh + 977);
    const auto *mh_978 = buffer.data(mh + 978);
    const auto *mh_979 = buffer.data(mh + 979);
    const auto *mh_980 = buffer.data(mh + 980);
    const auto *mh_981 = buffer.data(mh + 981);
    const auto *mh_982 = buffer.data(mh + 982);
    const auto *mh_983 = buffer.data(mh + 983);
    const auto *mh_984 = buffer.data(mh + 984);
    const auto *mh_985 = buffer.data(mh + 985);
    const auto *mh_986 = buffer.data(mh + 986);
    const auto *mh_987 = buffer.data(mh + 987);
    const auto *mh_988 = buffer.data(mh + 988);
    const auto *mh_989 = buffer.data(mh + 989);
    const auto *mh_990 = buffer.data(mh + 990);
    const auto *mh_991 = buffer.data(mh + 991);
    const auto *mh_992 = buffer.data(mh + 992);
    const auto *mh_993 = buffer.data(mh + 993);
    const auto *mh_994 = buffer.data(mh + 994);
    const auto *mh_995 = buffer.data(mh + 995);
    const auto *mh_996 = buffer.data(mh + 996);
    const auto *mh_997 = buffer.data(mh + 997);
    const auto *mh_998 = buffer.data(mh + 998);
    const auto *mh_999 = buffer.data(mh + 999);
    const auto *mh_1000 = buffer.data(mh + 1000);
    const auto *mh_1001 = buffer.data(mh + 1001);
    const auto *mh_1002 = buffer.data(mh + 1002);
    const auto *mh_1003 = buffer.data(mh + 1003);
    const auto *mh_1004 = buffer.data(mh + 1004);
    const auto *mh_1005 = buffer.data(mh + 1005);
    const auto *mh_1006 = buffer.data(mh + 1006);
    const auto *mh_1007 = buffer.data(mh + 1007);
    const auto *mh_1008 = buffer.data(mh + 1008);
    const auto *mh_1009 = buffer.data(mh + 1009);
    const auto *mh_1010 = buffer.data(mh + 1010);
    const auto *mh_1011 = buffer.data(mh + 1011);
    const auto *mh_1012 = buffer.data(mh + 1012);
    const auto *mh_1013 = buffer.data(mh + 1013);
    const auto *mh_1014 = buffer.data(mh + 1014);
    const auto *mh_1015 = buffer.data(mh + 1015);
    const auto *mh_1016 = buffer.data(mh + 1016);
    const auto *mh_1017 = buffer.data(mh + 1017);
    const auto *mh_1018 = buffer.data(mh + 1018);
    const auto *mh_1019 = buffer.data(mh + 1019);
    const auto *mh_1020 = buffer.data(mh + 1020);
    const auto *mh_1021 = buffer.data(mh + 1021);

#pragma omp simd aligned(t_654, t_655, t_656, t_657, t_658, kh_486, kh_487, kh_488, kh_489, \
                         kh_490, mh_843, mh_844, mh_845, mh_846, \
                         mh_847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = -3.0 * kh_486[k]
                   + f_0 * mh_843[k];

        t_655[k] = -3.0 * kh_487[k]
                   + f_0 * mh_844[k];

        t_656[k] = -3.0 * kh_488[k]
                   + f_0 * mh_845[k];

        t_657[k] = -3.0 * kh_489[k]
                   + f_0 * mh_846[k];

        t_658[k] = -3.0 * kh_490[k]
                   + f_0 * mh_847[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, kh_491, kh_492, kh_493, kh_494, \
                         kh_495, mh_848, mh_849, mh_850, mh_851, \
                         mh_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = -3.0 * kh_491[k]
                   + f_0 * mh_848[k];

        t_660[k] = -3.0 * kh_492[k]
                   + f_0 * mh_849[k];

        t_661[k] = -3.0 * kh_493[k]
                   + f_0 * mh_850[k];

        t_662[k] = -3.0 * kh_494[k]
                   + f_0 * mh_851[k];

        t_663[k] = -3.0 * kh_495[k]
                   + f_0 * mh_852[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, t_668, kh_496, kh_497, kh_498, kh_499, \
                         kh_500, mh_853, mh_854, mh_855, mh_856, \
                         mh_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = -3.0 * kh_496[k]
                   + f_0 * mh_853[k];

        t_665[k] = -3.0 * kh_497[k]
                   + f_0 * mh_854[k];

        t_666[k] = -3.0 * kh_498[k]
                   + f_0 * mh_855[k];

        t_667[k] = -3.0 * kh_499[k]
                   + f_0 * mh_856[k];

        t_668[k] = -3.0 * kh_500[k]
                   + f_0 * mh_857[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, t_672, t_673, kh_501, kh_502, kh_503, kh_504, \
                         kh_505, mh_858, mh_859, mh_860, mh_861, \
                         mh_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = -3.0 * kh_501[k]
                   + f_0 * mh_858[k];

        t_670[k] = -3.0 * kh_502[k]
                   + f_0 * mh_859[k];

        t_671[k] = -3.0 * kh_503[k]
                   + f_0 * mh_860[k];

        t_672[k] = -4.0 * kh_504[k]
                   + f_0 * mh_861[k];

        t_673[k] = -4.0 * kh_505[k]
                   + f_0 * mh_862[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, t_677, t_678, kh_506, kh_507, kh_508, kh_509, \
                         kh_510, mh_863, mh_864, mh_865, mh_866, \
                         mh_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = -4.0 * kh_506[k]
                   + f_0 * mh_863[k];

        t_675[k] = -4.0 * kh_507[k]
                   + f_0 * mh_864[k];

        t_676[k] = -4.0 * kh_508[k]
                   + f_0 * mh_865[k];

        t_677[k] = -4.0 * kh_509[k]
                   + f_0 * mh_866[k];

        t_678[k] = -4.0 * kh_510[k]
                   + f_0 * mh_867[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, t_682, t_683, kh_511, kh_512, kh_513, kh_514, \
                         kh_515, mh_868, mh_869, mh_870, mh_871, \
                         mh_872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = -4.0 * kh_511[k]
                   + f_0 * mh_868[k];

        t_680[k] = -4.0 * kh_512[k]
                   + f_0 * mh_869[k];

        t_681[k] = -4.0 * kh_513[k]
                   + f_0 * mh_870[k];

        t_682[k] = -4.0 * kh_514[k]
                   + f_0 * mh_871[k];

        t_683[k] = -4.0 * kh_515[k]
                   + f_0 * mh_872[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, t_688, kh_516, kh_517, kh_518, kh_519, \
                         kh_520, mh_873, mh_874, mh_875, mh_876, \
                         mh_877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = -4.0 * kh_516[k]
                   + f_0 * mh_873[k];

        t_685[k] = -4.0 * kh_517[k]
                   + f_0 * mh_874[k];

        t_686[k] = -4.0 * kh_518[k]
                   + f_0 * mh_875[k];

        t_687[k] = -4.0 * kh_519[k]
                   + f_0 * mh_876[k];

        t_688[k] = -4.0 * kh_520[k]
                   + f_0 * mh_877[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, t_693, kh_521, kh_522, kh_523, kh_524, \
                         kh_525, mh_878, mh_879, mh_880, mh_881, \
                         mh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = -4.0 * kh_521[k]
                   + f_0 * mh_878[k];

        t_690[k] = -4.0 * kh_522[k]
                   + f_0 * mh_879[k];

        t_691[k] = -4.0 * kh_523[k]
                   + f_0 * mh_880[k];

        t_692[k] = -4.0 * kh_524[k]
                   + f_0 * mh_881[k];

        t_693[k] = -5.0 * kh_525[k]
                   + f_0 * mh_882[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, t_698, kh_526, kh_527, kh_528, kh_529, \
                         kh_530, mh_883, mh_884, mh_885, mh_886, \
                         mh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = -5.0 * kh_526[k]
                   + f_0 * mh_883[k];

        t_695[k] = -5.0 * kh_527[k]
                   + f_0 * mh_884[k];

        t_696[k] = -5.0 * kh_528[k]
                   + f_0 * mh_885[k];

        t_697[k] = -5.0 * kh_529[k]
                   + f_0 * mh_886[k];

        t_698[k] = -5.0 * kh_530[k]
                   + f_0 * mh_887[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, t_703, kh_531, kh_532, kh_533, kh_534, \
                         kh_535, mh_888, mh_889, mh_890, mh_891, \
                         mh_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = -5.0 * kh_531[k]
                   + f_0 * mh_888[k];

        t_700[k] = -5.0 * kh_532[k]
                   + f_0 * mh_889[k];

        t_701[k] = -5.0 * kh_533[k]
                   + f_0 * mh_890[k];

        t_702[k] = -5.0 * kh_534[k]
                   + f_0 * mh_891[k];

        t_703[k] = -5.0 * kh_535[k]
                   + f_0 * mh_892[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, kh_536, kh_537, kh_538, kh_539, \
                         kh_540, mh_893, mh_894, mh_895, mh_896, \
                         mh_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = -5.0 * kh_536[k]
                   + f_0 * mh_893[k];

        t_705[k] = -5.0 * kh_537[k]
                   + f_0 * mh_894[k];

        t_706[k] = -5.0 * kh_538[k]
                   + f_0 * mh_895[k];

        t_707[k] = -5.0 * kh_539[k]
                   + f_0 * mh_896[k];

        t_708[k] = -5.0 * kh_540[k]
                   + f_0 * mh_897[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, t_712, t_713, kh_541, kh_542, kh_543, kh_544, \
                         kh_545, mh_898, mh_899, mh_900, mh_901, \
                         mh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = -5.0 * kh_541[k]
                   + f_0 * mh_898[k];

        t_710[k] = -5.0 * kh_542[k]
                   + f_0 * mh_899[k];

        t_711[k] = -5.0 * kh_543[k]
                   + f_0 * mh_900[k];

        t_712[k] = -5.0 * kh_544[k]
                   + f_0 * mh_901[k];

        t_713[k] = -5.0 * kh_545[k]
                   + f_0 * mh_902[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, t_718, kh_546, kh_547, kh_548, kh_549, \
                         kh_550, mh_903, mh_904, mh_905, mh_906, \
                         mh_907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = -6.0 * kh_546[k]
                   + f_0 * mh_903[k];

        t_715[k] = -6.0 * kh_547[k]
                   + f_0 * mh_904[k];

        t_716[k] = -6.0 * kh_548[k]
                   + f_0 * mh_905[k];

        t_717[k] = -6.0 * kh_549[k]
                   + f_0 * mh_906[k];

        t_718[k] = -6.0 * kh_550[k]
                   + f_0 * mh_907[k];
    }

#pragma omp simd aligned(t_719, t_720, t_721, t_722, t_723, kh_551, kh_552, kh_553, kh_554, \
                         kh_555, mh_908, mh_909, mh_910, mh_911, \
                         mh_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_719[k] = -6.0 * kh_551[k]
                   + f_0 * mh_908[k];

        t_720[k] = -6.0 * kh_552[k]
                   + f_0 * mh_909[k];

        t_721[k] = -6.0 * kh_553[k]
                   + f_0 * mh_910[k];

        t_722[k] = -6.0 * kh_554[k]
                   + f_0 * mh_911[k];

        t_723[k] = -6.0 * kh_555[k]
                   + f_0 * mh_912[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, kh_556, kh_557, kh_558, kh_559, \
                         kh_560, mh_913, mh_914, mh_915, mh_916, \
                         mh_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = -6.0 * kh_556[k]
                   + f_0 * mh_913[k];

        t_725[k] = -6.0 * kh_557[k]
                   + f_0 * mh_914[k];

        t_726[k] = -6.0 * kh_558[k]
                   + f_0 * mh_915[k];

        t_727[k] = -6.0 * kh_559[k]
                   + f_0 * mh_916[k];

        t_728[k] = -6.0 * kh_560[k]
                   + f_0 * mh_917[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, t_733, kh_561, kh_562, kh_563, kh_564, \
                         kh_565, mh_918, mh_919, mh_920, mh_921, \
                         mh_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = -6.0 * kh_561[k]
                   + f_0 * mh_918[k];

        t_730[k] = -6.0 * kh_562[k]
                   + f_0 * mh_919[k];

        t_731[k] = -6.0 * kh_563[k]
                   + f_0 * mh_920[k];

        t_732[k] = -6.0 * kh_564[k]
                   + f_0 * mh_921[k];

        t_733[k] = -6.0 * kh_565[k]
                   + f_0 * mh_922[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, t_738, kh_566, kh_567, kh_568, kh_569, \
                         kh_570, mh_923, mh_924, mh_925, mh_926, \
                         mh_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = -6.0 * kh_566[k]
                   + f_0 * mh_923[k];

        t_735[k] = -7.0 * kh_567[k]
                   + f_0 * mh_924[k];

        t_736[k] = -7.0 * kh_568[k]
                   + f_0 * mh_925[k];

        t_737[k] = -7.0 * kh_569[k]
                   + f_0 * mh_926[k];

        t_738[k] = -7.0 * kh_570[k]
                   + f_0 * mh_927[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, t_743, kh_571, kh_572, kh_573, kh_574, \
                         kh_575, mh_928, mh_929, mh_930, mh_931, \
                         mh_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = -7.0 * kh_571[k]
                   + f_0 * mh_928[k];

        t_740[k] = -7.0 * kh_572[k]
                   + f_0 * mh_929[k];

        t_741[k] = -7.0 * kh_573[k]
                   + f_0 * mh_930[k];

        t_742[k] = -7.0 * kh_574[k]
                   + f_0 * mh_931[k];

        t_743[k] = -7.0 * kh_575[k]
                   + f_0 * mh_932[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, t_748, kh_576, kh_577, kh_578, kh_579, \
                         kh_580, mh_933, mh_934, mh_935, mh_936, \
                         mh_937 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = -7.0 * kh_576[k]
                   + f_0 * mh_933[k];

        t_745[k] = -7.0 * kh_577[k]
                   + f_0 * mh_934[k];

        t_746[k] = -7.0 * kh_578[k]
                   + f_0 * mh_935[k];

        t_747[k] = -7.0 * kh_579[k]
                   + f_0 * mh_936[k];

        t_748[k] = -7.0 * kh_580[k]
                   + f_0 * mh_937[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, t_752, t_753, kh_581, kh_582, kh_583, kh_584, \
                         kh_585, mh_938, mh_939, mh_940, mh_941, \
                         mh_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = -7.0 * kh_581[k]
                   + f_0 * mh_938[k];

        t_750[k] = -7.0 * kh_582[k]
                   + f_0 * mh_939[k];

        t_751[k] = -7.0 * kh_583[k]
                   + f_0 * mh_940[k];

        t_752[k] = -7.0 * kh_584[k]
                   + f_0 * mh_941[k];

        t_753[k] = -7.0 * kh_585[k]
                   + f_0 * mh_942[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, t_758, t_759, t_760, kh_586, kh_587, \
                         mh_943, mh_944, mh_966, mh_967, mh_968, mh_969, \
                         mh_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = -7.0 * kh_586[k]
                   + f_0 * mh_943[k];

        t_755[k] = -7.0 * kh_587[k]
                   + f_0 * mh_944[k];

        t_756[k] = f_0 * mh_966[k];

        t_757[k] = f_0 * mh_967[k];

        t_758[k] = f_0 * mh_968[k];

        t_759[k] = f_0 * mh_969[k];

        t_760[k] = f_0 * mh_970[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, t_764, t_765, t_766, t_767, t_768, mh_971, \
                         mh_972, mh_973, mh_974, mh_975, mh_976, mh_977, \
                         mh_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_0 * mh_971[k];

        t_762[k] = f_0 * mh_972[k];

        t_763[k] = f_0 * mh_973[k];

        t_764[k] = f_0 * mh_974[k];

        t_765[k] = f_0 * mh_975[k];

        t_766[k] = f_0 * mh_976[k];

        t_767[k] = f_0 * mh_977[k];

        t_768[k] = f_0 * mh_978[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, t_772, t_773, t_774, t_775, t_776, mh_979, \
                         mh_980, mh_981, mh_982, mh_983, mh_984, mh_985, \
                         mh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_0 * mh_979[k];

        t_770[k] = f_0 * mh_980[k];

        t_771[k] = f_0 * mh_981[k];

        t_772[k] = f_0 * mh_982[k];

        t_773[k] = f_0 * mh_983[k];

        t_774[k] = f_0 * mh_984[k];

        t_775[k] = f_0 * mh_985[k];

        t_776[k] = f_0 * mh_986[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, t_780, t_781, kh_588, kh_589, kh_590, kh_591, \
                         kh_592, mh_987, mh_988, mh_989, mh_990, \
                         mh_991 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = -kh_588[k]
                   + f_0 * mh_987[k];

        t_778[k] = -kh_589[k]
                   + f_0 * mh_988[k];

        t_779[k] = -kh_590[k]
                   + f_0 * mh_989[k];

        t_780[k] = -kh_591[k]
                   + f_0 * mh_990[k];

        t_781[k] = -kh_592[k]
                   + f_0 * mh_991[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, t_785, t_786, kh_593, kh_594, kh_595, kh_596, \
                         kh_597, mh_992, mh_993, mh_994, mh_995, \
                         mh_996 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = -kh_593[k]
                   + f_0 * mh_992[k];

        t_783[k] = -kh_594[k]
                   + f_0 * mh_993[k];

        t_784[k] = -kh_595[k]
                   + f_0 * mh_994[k];

        t_785[k] = -kh_596[k]
                   + f_0 * mh_995[k];

        t_786[k] = -kh_597[k]
                   + f_0 * mh_996[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, t_790, t_791, kh_598, kh_599, kh_600, kh_601, \
                         kh_602, mh_997, mh_998, mh_999, mh_1000, \
                         mh_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = -kh_598[k]
                   + f_0 * mh_997[k];

        t_788[k] = -kh_599[k]
                   + f_0 * mh_998[k];

        t_789[k] = -kh_600[k]
                   + f_0 * mh_999[k];

        t_790[k] = -kh_601[k]
                   + f_0 * mh_1000[k];

        t_791[k] = -kh_602[k]
                   + f_0 * mh_1001[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, t_796, kh_603, kh_604, kh_605, kh_606, \
                         kh_607, mh_1002, mh_1003, mh_1004, mh_1005, \
                         mh_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = -kh_603[k]
                   + f_0 * mh_1002[k];

        t_793[k] = -kh_604[k]
                   + f_0 * mh_1003[k];

        t_794[k] = -kh_605[k]
                   + f_0 * mh_1004[k];

        t_795[k] = -kh_606[k]
                   + f_0 * mh_1005[k];

        t_796[k] = -kh_607[k]
                   + f_0 * mh_1006[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, t_800, t_801, kh_608, kh_609, kh_610, kh_611, \
                         kh_612, mh_1007, mh_1008, mh_1009, mh_1010, \
                         mh_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = -kh_608[k]
                   + f_0 * mh_1007[k];

        t_798[k] = -2.0 * kh_609[k]
                   + f_0 * mh_1008[k];

        t_799[k] = -2.0 * kh_610[k]
                   + f_0 * mh_1009[k];

        t_800[k] = -2.0 * kh_611[k]
                   + f_0 * mh_1010[k];

        t_801[k] = -2.0 * kh_612[k]
                   + f_0 * mh_1011[k];
    }

#pragma omp simd aligned(t_802, t_803, t_804, t_805, t_806, kh_613, kh_614, kh_615, kh_616, \
                         kh_617, mh_1012, mh_1013, mh_1014, mh_1015, \
                         mh_1016 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_802[k] = -2.0 * kh_613[k]
                   + f_0 * mh_1012[k];

        t_803[k] = -2.0 * kh_614[k]
                   + f_0 * mh_1013[k];

        t_804[k] = -2.0 * kh_615[k]
                   + f_0 * mh_1014[k];

        t_805[k] = -2.0 * kh_616[k]
                   + f_0 * mh_1015[k];

        t_806[k] = -2.0 * kh_617[k]
                   + f_0 * mh_1016[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, t_810, t_811, kh_618, kh_619, kh_620, kh_621, \
                         kh_622, mh_1017, mh_1018, mh_1019, mh_1020, \
                         mh_1021 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = -2.0 * kh_618[k]
                   + f_0 * mh_1017[k];

        t_808[k] = -2.0 * kh_619[k]
                   + f_0 * mh_1018[k];

        t_809[k] = -2.0 * kh_620[k]
                   + f_0 * mh_1019[k];

        t_810[k] = -2.0 * kh_621[k]
                   + f_0 * mh_1020[k];

        t_811[k] = -2.0 * kh_622[k]
                   + f_0 * mh_1021[k];
    }
}

static auto
compute_prim_geom_10_lh_electron_repulsion_2_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t mh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);
    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);
    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);
    auto *t_840 = buffer.data(target + 840);
    auto *t_841 = buffer.data(target + 841);
    auto *t_842 = buffer.data(target + 842);
    auto *t_843 = buffer.data(target + 843);
    auto *t_844 = buffer.data(target + 844);
    auto *t_845 = buffer.data(target + 845);
    auto *t_846 = buffer.data(target + 846);
    auto *t_847 = buffer.data(target + 847);
    auto *t_848 = buffer.data(target + 848);
    auto *t_849 = buffer.data(target + 849);
    auto *t_850 = buffer.data(target + 850);
    auto *t_851 = buffer.data(target + 851);
    auto *t_852 = buffer.data(target + 852);
    auto *t_853 = buffer.data(target + 853);
    auto *t_854 = buffer.data(target + 854);
    auto *t_855 = buffer.data(target + 855);
    auto *t_856 = buffer.data(target + 856);
    auto *t_857 = buffer.data(target + 857);
    auto *t_858 = buffer.data(target + 858);
    auto *t_859 = buffer.data(target + 859);
    auto *t_860 = buffer.data(target + 860);
    auto *t_861 = buffer.data(target + 861);
    auto *t_862 = buffer.data(target + 862);
    auto *t_863 = buffer.data(target + 863);
    auto *t_864 = buffer.data(target + 864);
    auto *t_865 = buffer.data(target + 865);
    auto *t_866 = buffer.data(target + 866);
    auto *t_867 = buffer.data(target + 867);
    auto *t_868 = buffer.data(target + 868);
    auto *t_869 = buffer.data(target + 869);
    auto *t_870 = buffer.data(target + 870);
    auto *t_871 = buffer.data(target + 871);
    auto *t_872 = buffer.data(target + 872);
    auto *t_873 = buffer.data(target + 873);
    auto *t_874 = buffer.data(target + 874);
    auto *t_875 = buffer.data(target + 875);
    auto *t_876 = buffer.data(target + 876);
    auto *t_877 = buffer.data(target + 877);
    auto *t_878 = buffer.data(target + 878);
    auto *t_879 = buffer.data(target + 879);
    auto *t_880 = buffer.data(target + 880);
    auto *t_881 = buffer.data(target + 881);
    auto *t_882 = buffer.data(target + 882);
    auto *t_883 = buffer.data(target + 883);
    auto *t_884 = buffer.data(target + 884);
    auto *t_885 = buffer.data(target + 885);
    auto *t_886 = buffer.data(target + 886);
    auto *t_887 = buffer.data(target + 887);
    auto *t_888 = buffer.data(target + 888);
    auto *t_889 = buffer.data(target + 889);
    auto *t_890 = buffer.data(target + 890);
    auto *t_891 = buffer.data(target + 891);
    auto *t_892 = buffer.data(target + 892);
    auto *t_893 = buffer.data(target + 893);
    auto *t_894 = buffer.data(target + 894);
    auto *t_895 = buffer.data(target + 895);
    auto *t_896 = buffer.data(target + 896);
    auto *t_897 = buffer.data(target + 897);
    auto *t_898 = buffer.data(target + 898);
    auto *t_899 = buffer.data(target + 899);
    auto *t_900 = buffer.data(target + 900);
    auto *t_901 = buffer.data(target + 901);
    auto *t_902 = buffer.data(target + 902);
    auto *t_903 = buffer.data(target + 903);
    auto *t_904 = buffer.data(target + 904);
    auto *t_905 = buffer.data(target + 905);
    auto *t_906 = buffer.data(target + 906);
    auto *t_907 = buffer.data(target + 907);
    auto *t_908 = buffer.data(target + 908);
    auto *t_909 = buffer.data(target + 909);
    auto *t_910 = buffer.data(target + 910);
    auto *t_911 = buffer.data(target + 911);
    auto *t_912 = buffer.data(target + 912);
    auto *t_913 = buffer.data(target + 913);
    auto *t_914 = buffer.data(target + 914);
    auto *t_915 = buffer.data(target + 915);
    auto *t_916 = buffer.data(target + 916);
    auto *t_917 = buffer.data(target + 917);
    auto *t_918 = buffer.data(target + 918);
    auto *t_919 = buffer.data(target + 919);
    auto *t_920 = buffer.data(target + 920);
    auto *t_921 = buffer.data(target + 921);
    auto *t_922 = buffer.data(target + 922);
    auto *t_923 = buffer.data(target + 923);
    auto *t_924 = buffer.data(target + 924);
    auto *t_925 = buffer.data(target + 925);
    auto *t_926 = buffer.data(target + 926);
    auto *t_927 = buffer.data(target + 927);
    auto *t_928 = buffer.data(target + 928);
    auto *t_929 = buffer.data(target + 929);
    auto *t_930 = buffer.data(target + 930);
    auto *t_931 = buffer.data(target + 931);
    auto *t_932 = buffer.data(target + 932);
    auto *t_933 = buffer.data(target + 933);
    auto *t_934 = buffer.data(target + 934);
    auto *t_935 = buffer.data(target + 935);
    auto *t_936 = buffer.data(target + 936);
    auto *t_937 = buffer.data(target + 937);
    auto *t_938 = buffer.data(target + 938);
    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);

    const auto *kh_623 = buffer.data(kh + 623);
    const auto *kh_624 = buffer.data(kh + 624);
    const auto *kh_625 = buffer.data(kh + 625);
    const auto *kh_626 = buffer.data(kh + 626);
    const auto *kh_627 = buffer.data(kh + 627);
    const auto *kh_628 = buffer.data(kh + 628);
    const auto *kh_629 = buffer.data(kh + 629);
    const auto *kh_630 = buffer.data(kh + 630);
    const auto *kh_631 = buffer.data(kh + 631);
    const auto *kh_632 = buffer.data(kh + 632);
    const auto *kh_633 = buffer.data(kh + 633);
    const auto *kh_634 = buffer.data(kh + 634);
    const auto *kh_635 = buffer.data(kh + 635);
    const auto *kh_636 = buffer.data(kh + 636);
    const auto *kh_637 = buffer.data(kh + 637);
    const auto *kh_638 = buffer.data(kh + 638);
    const auto *kh_639 = buffer.data(kh + 639);
    const auto *kh_640 = buffer.data(kh + 640);
    const auto *kh_641 = buffer.data(kh + 641);
    const auto *kh_642 = buffer.data(kh + 642);
    const auto *kh_643 = buffer.data(kh + 643);
    const auto *kh_644 = buffer.data(kh + 644);
    const auto *kh_645 = buffer.data(kh + 645);
    const auto *kh_646 = buffer.data(kh + 646);
    const auto *kh_647 = buffer.data(kh + 647);
    const auto *kh_648 = buffer.data(kh + 648);
    const auto *kh_649 = buffer.data(kh + 649);
    const auto *kh_650 = buffer.data(kh + 650);
    const auto *kh_651 = buffer.data(kh + 651);
    const auto *kh_652 = buffer.data(kh + 652);
    const auto *kh_653 = buffer.data(kh + 653);
    const auto *kh_654 = buffer.data(kh + 654);
    const auto *kh_655 = buffer.data(kh + 655);
    const auto *kh_656 = buffer.data(kh + 656);
    const auto *kh_657 = buffer.data(kh + 657);
    const auto *kh_658 = buffer.data(kh + 658);
    const auto *kh_659 = buffer.data(kh + 659);
    const auto *kh_660 = buffer.data(kh + 660);
    const auto *kh_661 = buffer.data(kh + 661);
    const auto *kh_662 = buffer.data(kh + 662);
    const auto *kh_663 = buffer.data(kh + 663);
    const auto *kh_664 = buffer.data(kh + 664);
    const auto *kh_665 = buffer.data(kh + 665);
    const auto *kh_666 = buffer.data(kh + 666);
    const auto *kh_667 = buffer.data(kh + 667);
    const auto *kh_668 = buffer.data(kh + 668);
    const auto *kh_669 = buffer.data(kh + 669);
    const auto *kh_670 = buffer.data(kh + 670);
    const auto *kh_671 = buffer.data(kh + 671);
    const auto *kh_672 = buffer.data(kh + 672);
    const auto *kh_673 = buffer.data(kh + 673);
    const auto *kh_674 = buffer.data(kh + 674);
    const auto *kh_675 = buffer.data(kh + 675);
    const auto *kh_676 = buffer.data(kh + 676);
    const auto *kh_677 = buffer.data(kh + 677);
    const auto *kh_678 = buffer.data(kh + 678);
    const auto *kh_679 = buffer.data(kh + 679);
    const auto *kh_680 = buffer.data(kh + 680);
    const auto *kh_681 = buffer.data(kh + 681);
    const auto *kh_682 = buffer.data(kh + 682);
    const auto *kh_683 = buffer.data(kh + 683);
    const auto *kh_684 = buffer.data(kh + 684);
    const auto *kh_685 = buffer.data(kh + 685);
    const auto *kh_686 = buffer.data(kh + 686);
    const auto *kh_687 = buffer.data(kh + 687);
    const auto *kh_688 = buffer.data(kh + 688);
    const auto *kh_689 = buffer.data(kh + 689);
    const auto *kh_690 = buffer.data(kh + 690);
    const auto *kh_691 = buffer.data(kh + 691);
    const auto *kh_692 = buffer.data(kh + 692);
    const auto *kh_693 = buffer.data(kh + 693);
    const auto *kh_694 = buffer.data(kh + 694);
    const auto *kh_695 = buffer.data(kh + 695);
    const auto *kh_696 = buffer.data(kh + 696);
    const auto *kh_697 = buffer.data(kh + 697);
    const auto *kh_698 = buffer.data(kh + 698);
    const auto *kh_699 = buffer.data(kh + 699);
    const auto *kh_700 = buffer.data(kh + 700);
    const auto *kh_701 = buffer.data(kh + 701);
    const auto *kh_702 = buffer.data(kh + 702);
    const auto *kh_703 = buffer.data(kh + 703);
    const auto *kh_704 = buffer.data(kh + 704);
    const auto *kh_705 = buffer.data(kh + 705);
    const auto *kh_706 = buffer.data(kh + 706);
    const auto *kh_707 = buffer.data(kh + 707);
    const auto *kh_708 = buffer.data(kh + 708);
    const auto *kh_709 = buffer.data(kh + 709);
    const auto *kh_710 = buffer.data(kh + 710);
    const auto *kh_711 = buffer.data(kh + 711);
    const auto *kh_712 = buffer.data(kh + 712);
    const auto *kh_713 = buffer.data(kh + 713);
    const auto *kh_714 = buffer.data(kh + 714);
    const auto *kh_715 = buffer.data(kh + 715);
    const auto *kh_716 = buffer.data(kh + 716);
    const auto *kh_717 = buffer.data(kh + 717);
    const auto *kh_718 = buffer.data(kh + 718);
    const auto *kh_719 = buffer.data(kh + 719);
    const auto *kh_720 = buffer.data(kh + 720);
    const auto *kh_721 = buffer.data(kh + 721);
    const auto *kh_722 = buffer.data(kh + 722);
    const auto *kh_723 = buffer.data(kh + 723);
    const auto *kh_724 = buffer.data(kh + 724);
    const auto *kh_725 = buffer.data(kh + 725);
    const auto *kh_726 = buffer.data(kh + 726);
    const auto *kh_727 = buffer.data(kh + 727);
    const auto *kh_728 = buffer.data(kh + 728);
    const auto *kh_729 = buffer.data(kh + 729);
    const auto *kh_730 = buffer.data(kh + 730);
    const auto *kh_731 = buffer.data(kh + 731);
    const auto *kh_732 = buffer.data(kh + 732);
    const auto *kh_733 = buffer.data(kh + 733);
    const auto *kh_734 = buffer.data(kh + 734);
    const auto *kh_735 = buffer.data(kh + 735);
    const auto *kh_736 = buffer.data(kh + 736);
    const auto *kh_737 = buffer.data(kh + 737);
    const auto *kh_738 = buffer.data(kh + 738);
    const auto *kh_739 = buffer.data(kh + 739);
    const auto *kh_740 = buffer.data(kh + 740);
    const auto *kh_741 = buffer.data(kh + 741);
    const auto *kh_742 = buffer.data(kh + 742);
    const auto *kh_743 = buffer.data(kh + 743);
    const auto *kh_744 = buffer.data(kh + 744);
    const auto *kh_745 = buffer.data(kh + 745);
    const auto *kh_746 = buffer.data(kh + 746);
    const auto *kh_747 = buffer.data(kh + 747);
    const auto *kh_748 = buffer.data(kh + 748);
    const auto *kh_749 = buffer.data(kh + 749);
    const auto *kh_750 = buffer.data(kh + 750);
    const auto *kh_751 = buffer.data(kh + 751);
    const auto *kh_752 = buffer.data(kh + 752);
    const auto *kh_753 = buffer.data(kh + 753);
    const auto *kh_754 = buffer.data(kh + 754);
    const auto *kh_755 = buffer.data(kh + 755);

    const auto *mh_1022 = buffer.data(mh + 1022);
    const auto *mh_1023 = buffer.data(mh + 1023);
    const auto *mh_1024 = buffer.data(mh + 1024);
    const auto *mh_1025 = buffer.data(mh + 1025);
    const auto *mh_1026 = buffer.data(mh + 1026);
    const auto *mh_1027 = buffer.data(mh + 1027);
    const auto *mh_1028 = buffer.data(mh + 1028);
    const auto *mh_1029 = buffer.data(mh + 1029);
    const auto *mh_1030 = buffer.data(mh + 1030);
    const auto *mh_1031 = buffer.data(mh + 1031);
    const auto *mh_1032 = buffer.data(mh + 1032);
    const auto *mh_1033 = buffer.data(mh + 1033);
    const auto *mh_1034 = buffer.data(mh + 1034);
    const auto *mh_1035 = buffer.data(mh + 1035);
    const auto *mh_1036 = buffer.data(mh + 1036);
    const auto *mh_1037 = buffer.data(mh + 1037);
    const auto *mh_1038 = buffer.data(mh + 1038);
    const auto *mh_1039 = buffer.data(mh + 1039);
    const auto *mh_1040 = buffer.data(mh + 1040);
    const auto *mh_1041 = buffer.data(mh + 1041);
    const auto *mh_1042 = buffer.data(mh + 1042);
    const auto *mh_1043 = buffer.data(mh + 1043);
    const auto *mh_1044 = buffer.data(mh + 1044);
    const auto *mh_1045 = buffer.data(mh + 1045);
    const auto *mh_1046 = buffer.data(mh + 1046);
    const auto *mh_1047 = buffer.data(mh + 1047);
    const auto *mh_1048 = buffer.data(mh + 1048);
    const auto *mh_1049 = buffer.data(mh + 1049);
    const auto *mh_1050 = buffer.data(mh + 1050);
    const auto *mh_1051 = buffer.data(mh + 1051);
    const auto *mh_1052 = buffer.data(mh + 1052);
    const auto *mh_1053 = buffer.data(mh + 1053);
    const auto *mh_1054 = buffer.data(mh + 1054);
    const auto *mh_1055 = buffer.data(mh + 1055);
    const auto *mh_1056 = buffer.data(mh + 1056);
    const auto *mh_1057 = buffer.data(mh + 1057);
    const auto *mh_1058 = buffer.data(mh + 1058);
    const auto *mh_1059 = buffer.data(mh + 1059);
    const auto *mh_1060 = buffer.data(mh + 1060);
    const auto *mh_1061 = buffer.data(mh + 1061);
    const auto *mh_1062 = buffer.data(mh + 1062);
    const auto *mh_1063 = buffer.data(mh + 1063);
    const auto *mh_1064 = buffer.data(mh + 1064);
    const auto *mh_1065 = buffer.data(mh + 1065);
    const auto *mh_1066 = buffer.data(mh + 1066);
    const auto *mh_1067 = buffer.data(mh + 1067);
    const auto *mh_1068 = buffer.data(mh + 1068);
    const auto *mh_1069 = buffer.data(mh + 1069);
    const auto *mh_1070 = buffer.data(mh + 1070);
    const auto *mh_1071 = buffer.data(mh + 1071);
    const auto *mh_1072 = buffer.data(mh + 1072);
    const auto *mh_1073 = buffer.data(mh + 1073);
    const auto *mh_1074 = buffer.data(mh + 1074);
    const auto *mh_1075 = buffer.data(mh + 1075);
    const auto *mh_1076 = buffer.data(mh + 1076);
    const auto *mh_1077 = buffer.data(mh + 1077);
    const auto *mh_1078 = buffer.data(mh + 1078);
    const auto *mh_1079 = buffer.data(mh + 1079);
    const auto *mh_1080 = buffer.data(mh + 1080);
    const auto *mh_1081 = buffer.data(mh + 1081);
    const auto *mh_1082 = buffer.data(mh + 1082);
    const auto *mh_1083 = buffer.data(mh + 1083);
    const auto *mh_1084 = buffer.data(mh + 1084);
    const auto *mh_1085 = buffer.data(mh + 1085);
    const auto *mh_1086 = buffer.data(mh + 1086);
    const auto *mh_1087 = buffer.data(mh + 1087);
    const auto *mh_1088 = buffer.data(mh + 1088);
    const auto *mh_1089 = buffer.data(mh + 1089);
    const auto *mh_1090 = buffer.data(mh + 1090);
    const auto *mh_1091 = buffer.data(mh + 1091);
    const auto *mh_1092 = buffer.data(mh + 1092);
    const auto *mh_1093 = buffer.data(mh + 1093);
    const auto *mh_1094 = buffer.data(mh + 1094);
    const auto *mh_1095 = buffer.data(mh + 1095);
    const auto *mh_1096 = buffer.data(mh + 1096);
    const auto *mh_1097 = buffer.data(mh + 1097);
    const auto *mh_1098 = buffer.data(mh + 1098);
    const auto *mh_1099 = buffer.data(mh + 1099);
    const auto *mh_1100 = buffer.data(mh + 1100);
    const auto *mh_1101 = buffer.data(mh + 1101);
    const auto *mh_1102 = buffer.data(mh + 1102);
    const auto *mh_1103 = buffer.data(mh + 1103);
    const auto *mh_1104 = buffer.data(mh + 1104);
    const auto *mh_1105 = buffer.data(mh + 1105);
    const auto *mh_1106 = buffer.data(mh + 1106);
    const auto *mh_1107 = buffer.data(mh + 1107);
    const auto *mh_1108 = buffer.data(mh + 1108);
    const auto *mh_1109 = buffer.data(mh + 1109);
    const auto *mh_1110 = buffer.data(mh + 1110);
    const auto *mh_1111 = buffer.data(mh + 1111);
    const auto *mh_1112 = buffer.data(mh + 1112);
    const auto *mh_1113 = buffer.data(mh + 1113);
    const auto *mh_1114 = buffer.data(mh + 1114);
    const auto *mh_1115 = buffer.data(mh + 1115);
    const auto *mh_1116 = buffer.data(mh + 1116);
    const auto *mh_1117 = buffer.data(mh + 1117);
    const auto *mh_1118 = buffer.data(mh + 1118);
    const auto *mh_1119 = buffer.data(mh + 1119);
    const auto *mh_1120 = buffer.data(mh + 1120);
    const auto *mh_1121 = buffer.data(mh + 1121);
    const auto *mh_1122 = buffer.data(mh + 1122);
    const auto *mh_1123 = buffer.data(mh + 1123);
    const auto *mh_1124 = buffer.data(mh + 1124);
    const auto *mh_1125 = buffer.data(mh + 1125);
    const auto *mh_1126 = buffer.data(mh + 1126);
    const auto *mh_1127 = buffer.data(mh + 1127);
    const auto *mh_1128 = buffer.data(mh + 1128);
    const auto *mh_1129 = buffer.data(mh + 1129);
    const auto *mh_1130 = buffer.data(mh + 1130);
    const auto *mh_1131 = buffer.data(mh + 1131);
    const auto *mh_1132 = buffer.data(mh + 1132);
    const auto *mh_1133 = buffer.data(mh + 1133);
    const auto *mh_1134 = buffer.data(mh + 1134);
    const auto *mh_1135 = buffer.data(mh + 1135);
    const auto *mh_1136 = buffer.data(mh + 1136);
    const auto *mh_1137 = buffer.data(mh + 1137);
    const auto *mh_1138 = buffer.data(mh + 1138);
    const auto *mh_1139 = buffer.data(mh + 1139);
    const auto *mh_1140 = buffer.data(mh + 1140);
    const auto *mh_1141 = buffer.data(mh + 1141);
    const auto *mh_1142 = buffer.data(mh + 1142);
    const auto *mh_1143 = buffer.data(mh + 1143);
    const auto *mh_1144 = buffer.data(mh + 1144);
    const auto *mh_1145 = buffer.data(mh + 1145);
    const auto *mh_1146 = buffer.data(mh + 1146);
    const auto *mh_1147 = buffer.data(mh + 1147);
    const auto *mh_1148 = buffer.data(mh + 1148);
    const auto *mh_1149 = buffer.data(mh + 1149);
    const auto *mh_1150 = buffer.data(mh + 1150);
    const auto *mh_1151 = buffer.data(mh + 1151);
    const auto *mh_1152 = buffer.data(mh + 1152);
    const auto *mh_1153 = buffer.data(mh + 1153);
    const auto *mh_1154 = buffer.data(mh + 1154);

#pragma omp simd aligned(t_812, t_813, t_814, t_815, t_816, kh_623, kh_624, kh_625, kh_626, \
                         kh_627, mh_1022, mh_1023, mh_1024, mh_1025, \
                         mh_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = -2.0 * kh_623[k]
                   + f_0 * mh_1022[k];

        t_813[k] = -2.0 * kh_624[k]
                   + f_0 * mh_1023[k];

        t_814[k] = -2.0 * kh_625[k]
                   + f_0 * mh_1024[k];

        t_815[k] = -2.0 * kh_626[k]
                   + f_0 * mh_1025[k];

        t_816[k] = -2.0 * kh_627[k]
                   + f_0 * mh_1026[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, t_820, t_821, kh_628, kh_629, kh_630, kh_631, \
                         kh_632, mh_1027, mh_1028, mh_1029, mh_1030, \
                         mh_1031 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = -2.0 * kh_628[k]
                   + f_0 * mh_1027[k];

        t_818[k] = -2.0 * kh_629[k]
                   + f_0 * mh_1028[k];

        t_819[k] = -3.0 * kh_630[k]
                   + f_0 * mh_1029[k];

        t_820[k] = -3.0 * kh_631[k]
                   + f_0 * mh_1030[k];

        t_821[k] = -3.0 * kh_632[k]
                   + f_0 * mh_1031[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, t_825, t_826, kh_633, kh_634, kh_635, kh_636, \
                         kh_637, mh_1032, mh_1033, mh_1034, mh_1035, \
                         mh_1036 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = -3.0 * kh_633[k]
                   + f_0 * mh_1032[k];

        t_823[k] = -3.0 * kh_634[k]
                   + f_0 * mh_1033[k];

        t_824[k] = -3.0 * kh_635[k]
                   + f_0 * mh_1034[k];

        t_825[k] = -3.0 * kh_636[k]
                   + f_0 * mh_1035[k];

        t_826[k] = -3.0 * kh_637[k]
                   + f_0 * mh_1036[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, t_830, t_831, kh_638, kh_639, kh_640, kh_641, \
                         kh_642, mh_1037, mh_1038, mh_1039, mh_1040, \
                         mh_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = -3.0 * kh_638[k]
                   + f_0 * mh_1037[k];

        t_828[k] = -3.0 * kh_639[k]
                   + f_0 * mh_1038[k];

        t_829[k] = -3.0 * kh_640[k]
                   + f_0 * mh_1039[k];

        t_830[k] = -3.0 * kh_641[k]
                   + f_0 * mh_1040[k];

        t_831[k] = -3.0 * kh_642[k]
                   + f_0 * mh_1041[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, t_835, t_836, kh_643, kh_644, kh_645, kh_646, \
                         kh_647, mh_1042, mh_1043, mh_1044, mh_1045, \
                         mh_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = -3.0 * kh_643[k]
                   + f_0 * mh_1042[k];

        t_833[k] = -3.0 * kh_644[k]
                   + f_0 * mh_1043[k];

        t_834[k] = -3.0 * kh_645[k]
                   + f_0 * mh_1044[k];

        t_835[k] = -3.0 * kh_646[k]
                   + f_0 * mh_1045[k];

        t_836[k] = -3.0 * kh_647[k]
                   + f_0 * mh_1046[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, t_840, t_841, kh_648, kh_649, kh_650, kh_651, \
                         kh_652, mh_1047, mh_1048, mh_1049, mh_1050, \
                         mh_1051 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = -3.0 * kh_648[k]
                   + f_0 * mh_1047[k];

        t_838[k] = -3.0 * kh_649[k]
                   + f_0 * mh_1048[k];

        t_839[k] = -3.0 * kh_650[k]
                   + f_0 * mh_1049[k];

        t_840[k] = -4.0 * kh_651[k]
                   + f_0 * mh_1050[k];

        t_841[k] = -4.0 * kh_652[k]
                   + f_0 * mh_1051[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, t_845, t_846, kh_653, kh_654, kh_655, kh_656, \
                         kh_657, mh_1052, mh_1053, mh_1054, mh_1055, \
                         mh_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = -4.0 * kh_653[k]
                   + f_0 * mh_1052[k];

        t_843[k] = -4.0 * kh_654[k]
                   + f_0 * mh_1053[k];

        t_844[k] = -4.0 * kh_655[k]
                   + f_0 * mh_1054[k];

        t_845[k] = -4.0 * kh_656[k]
                   + f_0 * mh_1055[k];

        t_846[k] = -4.0 * kh_657[k]
                   + f_0 * mh_1056[k];
    }

#pragma omp simd aligned(t_847, t_848, t_849, t_850, t_851, kh_658, kh_659, kh_660, kh_661, \
                         kh_662, mh_1057, mh_1058, mh_1059, mh_1060, \
                         mh_1061 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_847[k] = -4.0 * kh_658[k]
                   + f_0 * mh_1057[k];

        t_848[k] = -4.0 * kh_659[k]
                   + f_0 * mh_1058[k];

        t_849[k] = -4.0 * kh_660[k]
                   + f_0 * mh_1059[k];

        t_850[k] = -4.0 * kh_661[k]
                   + f_0 * mh_1060[k];

        t_851[k] = -4.0 * kh_662[k]
                   + f_0 * mh_1061[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, t_856, kh_663, kh_664, kh_665, kh_666, \
                         kh_667, mh_1062, mh_1063, mh_1064, mh_1065, \
                         mh_1066 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = -4.0 * kh_663[k]
                   + f_0 * mh_1062[k];

        t_853[k] = -4.0 * kh_664[k]
                   + f_0 * mh_1063[k];

        t_854[k] = -4.0 * kh_665[k]
                   + f_0 * mh_1064[k];

        t_855[k] = -4.0 * kh_666[k]
                   + f_0 * mh_1065[k];

        t_856[k] = -4.0 * kh_667[k]
                   + f_0 * mh_1066[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, t_860, t_861, kh_668, kh_669, kh_670, kh_671, \
                         kh_672, mh_1067, mh_1068, mh_1069, mh_1070, \
                         mh_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = -4.0 * kh_668[k]
                   + f_0 * mh_1067[k];

        t_858[k] = -4.0 * kh_669[k]
                   + f_0 * mh_1068[k];

        t_859[k] = -4.0 * kh_670[k]
                   + f_0 * mh_1069[k];

        t_860[k] = -4.0 * kh_671[k]
                   + f_0 * mh_1070[k];

        t_861[k] = -5.0 * kh_672[k]
                   + f_0 * mh_1071[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, t_866, kh_673, kh_674, kh_675, kh_676, \
                         kh_677, mh_1072, mh_1073, mh_1074, mh_1075, \
                         mh_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = -5.0 * kh_673[k]
                   + f_0 * mh_1072[k];

        t_863[k] = -5.0 * kh_674[k]
                   + f_0 * mh_1073[k];

        t_864[k] = -5.0 * kh_675[k]
                   + f_0 * mh_1074[k];

        t_865[k] = -5.0 * kh_676[k]
                   + f_0 * mh_1075[k];

        t_866[k] = -5.0 * kh_677[k]
                   + f_0 * mh_1076[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, t_870, t_871, kh_678, kh_679, kh_680, kh_681, \
                         kh_682, mh_1077, mh_1078, mh_1079, mh_1080, \
                         mh_1081 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = -5.0 * kh_678[k]
                   + f_0 * mh_1077[k];

        t_868[k] = -5.0 * kh_679[k]
                   + f_0 * mh_1078[k];

        t_869[k] = -5.0 * kh_680[k]
                   + f_0 * mh_1079[k];

        t_870[k] = -5.0 * kh_681[k]
                   + f_0 * mh_1080[k];

        t_871[k] = -5.0 * kh_682[k]
                   + f_0 * mh_1081[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, t_875, t_876, kh_683, kh_684, kh_685, kh_686, \
                         kh_687, mh_1082, mh_1083, mh_1084, mh_1085, \
                         mh_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = -5.0 * kh_683[k]
                   + f_0 * mh_1082[k];

        t_873[k] = -5.0 * kh_684[k]
                   + f_0 * mh_1083[k];

        t_874[k] = -5.0 * kh_685[k]
                   + f_0 * mh_1084[k];

        t_875[k] = -5.0 * kh_686[k]
                   + f_0 * mh_1085[k];

        t_876[k] = -5.0 * kh_687[k]
                   + f_0 * mh_1086[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, t_880, t_881, kh_688, kh_689, kh_690, kh_691, \
                         kh_692, mh_1087, mh_1088, mh_1089, mh_1090, \
                         mh_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = -5.0 * kh_688[k]
                   + f_0 * mh_1087[k];

        t_878[k] = -5.0 * kh_689[k]
                   + f_0 * mh_1088[k];

        t_879[k] = -5.0 * kh_690[k]
                   + f_0 * mh_1089[k];

        t_880[k] = -5.0 * kh_691[k]
                   + f_0 * mh_1090[k];

        t_881[k] = -5.0 * kh_692[k]
                   + f_0 * mh_1091[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, t_886, kh_693, kh_694, kh_695, kh_696, \
                         kh_697, mh_1092, mh_1093, mh_1094, mh_1095, \
                         mh_1096 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = -6.0 * kh_693[k]
                   + f_0 * mh_1092[k];

        t_883[k] = -6.0 * kh_694[k]
                   + f_0 * mh_1093[k];

        t_884[k] = -6.0 * kh_695[k]
                   + f_0 * mh_1094[k];

        t_885[k] = -6.0 * kh_696[k]
                   + f_0 * mh_1095[k];

        t_886[k] = -6.0 * kh_697[k]
                   + f_0 * mh_1096[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, t_890, t_891, kh_698, kh_699, kh_700, kh_701, \
                         kh_702, mh_1097, mh_1098, mh_1099, mh_1100, \
                         mh_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = -6.0 * kh_698[k]
                   + f_0 * mh_1097[k];

        t_888[k] = -6.0 * kh_699[k]
                   + f_0 * mh_1098[k];

        t_889[k] = -6.0 * kh_700[k]
                   + f_0 * mh_1099[k];

        t_890[k] = -6.0 * kh_701[k]
                   + f_0 * mh_1100[k];

        t_891[k] = -6.0 * kh_702[k]
                   + f_0 * mh_1101[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, t_895, t_896, kh_703, kh_704, kh_705, kh_706, \
                         kh_707, mh_1102, mh_1103, mh_1104, mh_1105, \
                         mh_1106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = -6.0 * kh_703[k]
                   + f_0 * mh_1102[k];

        t_893[k] = -6.0 * kh_704[k]
                   + f_0 * mh_1103[k];

        t_894[k] = -6.0 * kh_705[k]
                   + f_0 * mh_1104[k];

        t_895[k] = -6.0 * kh_706[k]
                   + f_0 * mh_1105[k];

        t_896[k] = -6.0 * kh_707[k]
                   + f_0 * mh_1106[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, t_900, t_901, kh_708, kh_709, kh_710, kh_711, \
                         kh_712, mh_1107, mh_1108, mh_1109, mh_1110, \
                         mh_1111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = -6.0 * kh_708[k]
                   + f_0 * mh_1107[k];

        t_898[k] = -6.0 * kh_709[k]
                   + f_0 * mh_1108[k];

        t_899[k] = -6.0 * kh_710[k]
                   + f_0 * mh_1109[k];

        t_900[k] = -6.0 * kh_711[k]
                   + f_0 * mh_1110[k];

        t_901[k] = -6.0 * kh_712[k]
                   + f_0 * mh_1111[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, t_905, t_906, kh_713, kh_714, kh_715, kh_716, \
                         kh_717, mh_1112, mh_1113, mh_1114, mh_1115, \
                         mh_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = -6.0 * kh_713[k]
                   + f_0 * mh_1112[k];

        t_903[k] = -7.0 * kh_714[k]
                   + f_0 * mh_1113[k];

        t_904[k] = -7.0 * kh_715[k]
                   + f_0 * mh_1114[k];

        t_905[k] = -7.0 * kh_716[k]
                   + f_0 * mh_1115[k];

        t_906[k] = -7.0 * kh_717[k]
                   + f_0 * mh_1116[k];
    }

#pragma omp simd aligned(t_907, t_908, t_909, t_910, t_911, kh_718, kh_719, kh_720, kh_721, \
                         kh_722, mh_1117, mh_1118, mh_1119, mh_1120, \
                         mh_1121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_907[k] = -7.0 * kh_718[k]
                   + f_0 * mh_1117[k];

        t_908[k] = -7.0 * kh_719[k]
                   + f_0 * mh_1118[k];

        t_909[k] = -7.0 * kh_720[k]
                   + f_0 * mh_1119[k];

        t_910[k] = -7.0 * kh_721[k]
                   + f_0 * mh_1120[k];

        t_911[k] = -7.0 * kh_722[k]
                   + f_0 * mh_1121[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, t_915, t_916, kh_723, kh_724, kh_725, kh_726, \
                         kh_727, mh_1122, mh_1123, mh_1124, mh_1125, \
                         mh_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = -7.0 * kh_723[k]
                   + f_0 * mh_1122[k];

        t_913[k] = -7.0 * kh_724[k]
                   + f_0 * mh_1123[k];

        t_914[k] = -7.0 * kh_725[k]
                   + f_0 * mh_1124[k];

        t_915[k] = -7.0 * kh_726[k]
                   + f_0 * mh_1125[k];

        t_916[k] = -7.0 * kh_727[k]
                   + f_0 * mh_1126[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, t_920, t_921, kh_728, kh_729, kh_730, kh_731, \
                         kh_732, mh_1127, mh_1128, mh_1129, mh_1130, \
                         mh_1131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = -7.0 * kh_728[k]
                   + f_0 * mh_1127[k];

        t_918[k] = -7.0 * kh_729[k]
                   + f_0 * mh_1128[k];

        t_919[k] = -7.0 * kh_730[k]
                   + f_0 * mh_1129[k];

        t_920[k] = -7.0 * kh_731[k]
                   + f_0 * mh_1130[k];

        t_921[k] = -7.0 * kh_732[k]
                   + f_0 * mh_1131[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, t_926, kh_733, kh_734, kh_735, kh_736, \
                         kh_737, mh_1132, mh_1133, mh_1134, mh_1135, \
                         mh_1136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = -7.0 * kh_733[k]
                   + f_0 * mh_1132[k];

        t_923[k] = -7.0 * kh_734[k]
                   + f_0 * mh_1133[k];

        t_924[k] = -8.0 * kh_735[k]
                   + f_0 * mh_1134[k];

        t_925[k] = -8.0 * kh_736[k]
                   + f_0 * mh_1135[k];

        t_926[k] = -8.0 * kh_737[k]
                   + f_0 * mh_1136[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, t_931, kh_738, kh_739, kh_740, kh_741, \
                         kh_742, mh_1137, mh_1138, mh_1139, mh_1140, \
                         mh_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = -8.0 * kh_738[k]
                   + f_0 * mh_1137[k];

        t_928[k] = -8.0 * kh_739[k]
                   + f_0 * mh_1138[k];

        t_929[k] = -8.0 * kh_740[k]
                   + f_0 * mh_1139[k];

        t_930[k] = -8.0 * kh_741[k]
                   + f_0 * mh_1140[k];

        t_931[k] = -8.0 * kh_742[k]
                   + f_0 * mh_1141[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, t_935, t_936, kh_743, kh_744, kh_745, kh_746, \
                         kh_747, mh_1142, mh_1143, mh_1144, mh_1145, \
                         mh_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = -8.0 * kh_743[k]
                   + f_0 * mh_1142[k];

        t_933[k] = -8.0 * kh_744[k]
                   + f_0 * mh_1143[k];

        t_934[k] = -8.0 * kh_745[k]
                   + f_0 * mh_1144[k];

        t_935[k] = -8.0 * kh_746[k]
                   + f_0 * mh_1145[k];

        t_936[k] = -8.0 * kh_747[k]
                   + f_0 * mh_1146[k];
    }

#pragma omp simd aligned(t_937, t_938, t_939, t_940, t_941, kh_748, kh_749, kh_750, kh_751, \
                         kh_752, mh_1147, mh_1148, mh_1149, mh_1150, \
                         mh_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_937[k] = -8.0 * kh_748[k]
                   + f_0 * mh_1147[k];

        t_938[k] = -8.0 * kh_749[k]
                   + f_0 * mh_1148[k];

        t_939[k] = -8.0 * kh_750[k]
                   + f_0 * mh_1149[k];

        t_940[k] = -8.0 * kh_751[k]
                   + f_0 * mh_1150[k];

        t_941[k] = -8.0 * kh_752[k]
                   + f_0 * mh_1151[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, kh_753, kh_754, kh_755, mh_1152, mh_1153, \
                         mh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = -8.0 * kh_753[k]
                   + f_0 * mh_1152[k];

        t_943[k] = -8.0 * kh_754[k]
                   + f_0 * mh_1153[k];

        t_944[k] = -8.0 * kh_755[k]
                   + f_0 * mh_1154[k];
    }
}

auto
compute_prim_geom_10_lh_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t kh, const size_t mh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_lh_electron_repulsion_2_piece0(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_2_piece1(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_2_piece2(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_2_piece3(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_2_piece4(buffer, target, kh, mh, ncols, alpha);

    compute_prim_geom_10_lh_electron_repulsion_2_piece5(buffer, target, kh, mh, ncols, alpha);
}

}  // namespace simdt2ceri
