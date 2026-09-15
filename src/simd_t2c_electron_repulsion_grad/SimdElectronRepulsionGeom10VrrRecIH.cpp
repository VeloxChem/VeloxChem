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


#include "SimdElectronRepulsionGeom10VrrRecIH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_ih_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hh, const size_t kh,
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

    const auto *hh_0 = buffer.data(hh + 0);
    const auto *hh_1 = buffer.data(hh + 1);
    const auto *hh_2 = buffer.data(hh + 2);
    const auto *hh_3 = buffer.data(hh + 3);
    const auto *hh_4 = buffer.data(hh + 4);
    const auto *hh_5 = buffer.data(hh + 5);
    const auto *hh_6 = buffer.data(hh + 6);
    const auto *hh_7 = buffer.data(hh + 7);
    const auto *hh_8 = buffer.data(hh + 8);
    const auto *hh_9 = buffer.data(hh + 9);
    const auto *hh_10 = buffer.data(hh + 10);
    const auto *hh_11 = buffer.data(hh + 11);
    const auto *hh_12 = buffer.data(hh + 12);
    const auto *hh_13 = buffer.data(hh + 13);
    const auto *hh_14 = buffer.data(hh + 14);
    const auto *hh_15 = buffer.data(hh + 15);
    const auto *hh_16 = buffer.data(hh + 16);
    const auto *hh_17 = buffer.data(hh + 17);
    const auto *hh_18 = buffer.data(hh + 18);
    const auto *hh_19 = buffer.data(hh + 19);
    const auto *hh_20 = buffer.data(hh + 20);
    const auto *hh_21 = buffer.data(hh + 21);
    const auto *hh_22 = buffer.data(hh + 22);
    const auto *hh_23 = buffer.data(hh + 23);
    const auto *hh_24 = buffer.data(hh + 24);
    const auto *hh_25 = buffer.data(hh + 25);
    const auto *hh_26 = buffer.data(hh + 26);
    const auto *hh_27 = buffer.data(hh + 27);
    const auto *hh_28 = buffer.data(hh + 28);
    const auto *hh_29 = buffer.data(hh + 29);
    const auto *hh_30 = buffer.data(hh + 30);
    const auto *hh_31 = buffer.data(hh + 31);
    const auto *hh_32 = buffer.data(hh + 32);
    const auto *hh_33 = buffer.data(hh + 33);
    const auto *hh_34 = buffer.data(hh + 34);
    const auto *hh_35 = buffer.data(hh + 35);
    const auto *hh_36 = buffer.data(hh + 36);
    const auto *hh_37 = buffer.data(hh + 37);
    const auto *hh_38 = buffer.data(hh + 38);
    const auto *hh_39 = buffer.data(hh + 39);
    const auto *hh_40 = buffer.data(hh + 40);
    const auto *hh_41 = buffer.data(hh + 41);
    const auto *hh_42 = buffer.data(hh + 42);
    const auto *hh_43 = buffer.data(hh + 43);
    const auto *hh_44 = buffer.data(hh + 44);
    const auto *hh_45 = buffer.data(hh + 45);
    const auto *hh_46 = buffer.data(hh + 46);
    const auto *hh_47 = buffer.data(hh + 47);
    const auto *hh_48 = buffer.data(hh + 48);
    const auto *hh_49 = buffer.data(hh + 49);
    const auto *hh_50 = buffer.data(hh + 50);
    const auto *hh_51 = buffer.data(hh + 51);
    const auto *hh_52 = buffer.data(hh + 52);
    const auto *hh_53 = buffer.data(hh + 53);
    const auto *hh_54 = buffer.data(hh + 54);
    const auto *hh_55 = buffer.data(hh + 55);
    const auto *hh_56 = buffer.data(hh + 56);
    const auto *hh_57 = buffer.data(hh + 57);
    const auto *hh_58 = buffer.data(hh + 58);
    const auto *hh_59 = buffer.data(hh + 59);
    const auto *hh_60 = buffer.data(hh + 60);
    const auto *hh_61 = buffer.data(hh + 61);
    const auto *hh_62 = buffer.data(hh + 62);
    const auto *hh_63 = buffer.data(hh + 63);
    const auto *hh_64 = buffer.data(hh + 64);
    const auto *hh_65 = buffer.data(hh + 65);
    const auto *hh_66 = buffer.data(hh + 66);
    const auto *hh_67 = buffer.data(hh + 67);
    const auto *hh_68 = buffer.data(hh + 68);
    const auto *hh_69 = buffer.data(hh + 69);
    const auto *hh_70 = buffer.data(hh + 70);
    const auto *hh_71 = buffer.data(hh + 71);
    const auto *hh_72 = buffer.data(hh + 72);
    const auto *hh_73 = buffer.data(hh + 73);
    const auto *hh_74 = buffer.data(hh + 74);
    const auto *hh_75 = buffer.data(hh + 75);
    const auto *hh_76 = buffer.data(hh + 76);
    const auto *hh_77 = buffer.data(hh + 77);
    const auto *hh_78 = buffer.data(hh + 78);
    const auto *hh_79 = buffer.data(hh + 79);
    const auto *hh_80 = buffer.data(hh + 80);
    const auto *hh_81 = buffer.data(hh + 81);
    const auto *hh_82 = buffer.data(hh + 82);
    const auto *hh_83 = buffer.data(hh + 83);
    const auto *hh_84 = buffer.data(hh + 84);
    const auto *hh_85 = buffer.data(hh + 85);
    const auto *hh_86 = buffer.data(hh + 86);
    const auto *hh_87 = buffer.data(hh + 87);
    const auto *hh_88 = buffer.data(hh + 88);
    const auto *hh_89 = buffer.data(hh + 89);
    const auto *hh_90 = buffer.data(hh + 90);
    const auto *hh_91 = buffer.data(hh + 91);
    const auto *hh_92 = buffer.data(hh + 92);
    const auto *hh_93 = buffer.data(hh + 93);
    const auto *hh_94 = buffer.data(hh + 94);
    const auto *hh_95 = buffer.data(hh + 95);
    const auto *hh_96 = buffer.data(hh + 96);
    const auto *hh_97 = buffer.data(hh + 97);
    const auto *hh_98 = buffer.data(hh + 98);
    const auto *hh_99 = buffer.data(hh + 99);
    const auto *hh_100 = buffer.data(hh + 100);
    const auto *hh_101 = buffer.data(hh + 101);
    const auto *hh_102 = buffer.data(hh + 102);
    const auto *hh_103 = buffer.data(hh + 103);
    const auto *hh_104 = buffer.data(hh + 104);
    const auto *hh_105 = buffer.data(hh + 105);
    const auto *hh_106 = buffer.data(hh + 106);
    const auto *hh_107 = buffer.data(hh + 107);
    const auto *hh_108 = buffer.data(hh + 108);
    const auto *hh_109 = buffer.data(hh + 109);
    const auto *hh_110 = buffer.data(hh + 110);
    const auto *hh_111 = buffer.data(hh + 111);
    const auto *hh_112 = buffer.data(hh + 112);
    const auto *hh_113 = buffer.data(hh + 113);
    const auto *hh_114 = buffer.data(hh + 114);
    const auto *hh_115 = buffer.data(hh + 115);
    const auto *hh_116 = buffer.data(hh + 116);
    const auto *hh_117 = buffer.data(hh + 117);
    const auto *hh_118 = buffer.data(hh + 118);
    const auto *hh_119 = buffer.data(hh + 119);
    const auto *hh_120 = buffer.data(hh + 120);
    const auto *hh_121 = buffer.data(hh + 121);
    const auto *hh_122 = buffer.data(hh + 122);
    const auto *hh_123 = buffer.data(hh + 123);
    const auto *hh_124 = buffer.data(hh + 124);
    const auto *hh_125 = buffer.data(hh + 125);
    const auto *hh_126 = buffer.data(hh + 126);
    const auto *hh_127 = buffer.data(hh + 127);
    const auto *hh_128 = buffer.data(hh + 128);
    const auto *hh_129 = buffer.data(hh + 129);
    const auto *hh_130 = buffer.data(hh + 130);
    const auto *hh_131 = buffer.data(hh + 131);
    const auto *hh_132 = buffer.data(hh + 132);
    const auto *hh_133 = buffer.data(hh + 133);
    const auto *hh_134 = buffer.data(hh + 134);
    const auto *hh_135 = buffer.data(hh + 135);
    const auto *hh_136 = buffer.data(hh + 136);
    const auto *hh_137 = buffer.data(hh + 137);
    const auto *hh_138 = buffer.data(hh + 138);
    const auto *hh_139 = buffer.data(hh + 139);
    const auto *hh_140 = buffer.data(hh + 140);
    const auto *hh_141 = buffer.data(hh + 141);
    const auto *hh_142 = buffer.data(hh + 142);
    const auto *hh_143 = buffer.data(hh + 143);
    const auto *hh_144 = buffer.data(hh + 144);
    const auto *hh_145 = buffer.data(hh + 145);
    const auto *hh_146 = buffer.data(hh + 146);
    const auto *hh_147 = buffer.data(hh + 147);
    const auto *hh_148 = buffer.data(hh + 148);
    const auto *hh_149 = buffer.data(hh + 149);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, hh_0, hh_1, hh_2, hh_3, hh_4, kh_0, kh_1, \
                         kh_2, kh_3, kh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -6.0 * hh_0[k]
                 + f_0 * kh_0[k];

        t_1[k] = -6.0 * hh_1[k]
                 + f_0 * kh_1[k];

        t_2[k] = -6.0 * hh_2[k]
                 + f_0 * kh_2[k];

        t_3[k] = -6.0 * hh_3[k]
                 + f_0 * kh_3[k];

        t_4[k] = -6.0 * hh_4[k]
                 + f_0 * kh_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, hh_5, hh_6, hh_7, hh_8, hh_9, kh_5, kh_6, \
                         kh_7, kh_8, kh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -6.0 * hh_5[k]
                 + f_0 * kh_5[k];

        t_6[k] = -6.0 * hh_6[k]
                 + f_0 * kh_6[k];

        t_7[k] = -6.0 * hh_7[k]
                 + f_0 * kh_7[k];

        t_8[k] = -6.0 * hh_8[k]
                 + f_0 * kh_8[k];

        t_9[k] = -6.0 * hh_9[k]
                 + f_0 * kh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, hh_10, hh_11, hh_12, hh_13, hh_14, \
                         kh_10, kh_11, kh_12, kh_13, kh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -6.0 * hh_10[k]
                  + f_0 * kh_10[k];

        t_11[k] = -6.0 * hh_11[k]
                  + f_0 * kh_11[k];

        t_12[k] = -6.0 * hh_12[k]
                  + f_0 * kh_12[k];

        t_13[k] = -6.0 * hh_13[k]
                  + f_0 * kh_13[k];

        t_14[k] = -6.0 * hh_14[k]
                  + f_0 * kh_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, hh_15, hh_16, hh_17, hh_18, hh_19, \
                         kh_15, kh_16, kh_17, kh_18, kh_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -6.0 * hh_15[k]
                  + f_0 * kh_15[k];

        t_16[k] = -6.0 * hh_16[k]
                  + f_0 * kh_16[k];

        t_17[k] = -6.0 * hh_17[k]
                  + f_0 * kh_17[k];

        t_18[k] = -6.0 * hh_18[k]
                  + f_0 * kh_18[k];

        t_19[k] = -6.0 * hh_19[k]
                  + f_0 * kh_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, hh_20, hh_21, hh_22, hh_23, hh_24, \
                         kh_20, kh_21, kh_22, kh_23, kh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -6.0 * hh_20[k]
                  + f_0 * kh_20[k];

        t_21[k] = -5.0 * hh_21[k]
                  + f_0 * kh_21[k];

        t_22[k] = -5.0 * hh_22[k]
                  + f_0 * kh_22[k];

        t_23[k] = -5.0 * hh_23[k]
                  + f_0 * kh_23[k];

        t_24[k] = -5.0 * hh_24[k]
                  + f_0 * kh_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, hh_25, hh_26, hh_27, hh_28, hh_29, \
                         kh_25, kh_26, kh_27, kh_28, kh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -5.0 * hh_25[k]
                  + f_0 * kh_25[k];

        t_26[k] = -5.0 * hh_26[k]
                  + f_0 * kh_26[k];

        t_27[k] = -5.0 * hh_27[k]
                  + f_0 * kh_27[k];

        t_28[k] = -5.0 * hh_28[k]
                  + f_0 * kh_28[k];

        t_29[k] = -5.0 * hh_29[k]
                  + f_0 * kh_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, hh_30, hh_31, hh_32, hh_33, hh_34, \
                         kh_30, kh_31, kh_32, kh_33, kh_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -5.0 * hh_30[k]
                  + f_0 * kh_30[k];

        t_31[k] = -5.0 * hh_31[k]
                  + f_0 * kh_31[k];

        t_32[k] = -5.0 * hh_32[k]
                  + f_0 * kh_32[k];

        t_33[k] = -5.0 * hh_33[k]
                  + f_0 * kh_33[k];

        t_34[k] = -5.0 * hh_34[k]
                  + f_0 * kh_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, hh_35, hh_36, hh_37, hh_38, hh_39, \
                         kh_35, kh_36, kh_37, kh_38, kh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -5.0 * hh_35[k]
                  + f_0 * kh_35[k];

        t_36[k] = -5.0 * hh_36[k]
                  + f_0 * kh_36[k];

        t_37[k] = -5.0 * hh_37[k]
                  + f_0 * kh_37[k];

        t_38[k] = -5.0 * hh_38[k]
                  + f_0 * kh_38[k];

        t_39[k] = -5.0 * hh_39[k]
                  + f_0 * kh_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, hh_40, hh_41, hh_42, hh_43, hh_44, \
                         kh_40, kh_41, kh_42, kh_43, kh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -5.0 * hh_40[k]
                  + f_0 * kh_40[k];

        t_41[k] = -5.0 * hh_41[k]
                  + f_0 * kh_41[k];

        t_42[k] = -5.0 * hh_42[k]
                  + f_0 * kh_42[k];

        t_43[k] = -5.0 * hh_43[k]
                  + f_0 * kh_43[k];

        t_44[k] = -5.0 * hh_44[k]
                  + f_0 * kh_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hh_45, hh_46, hh_47, hh_48, hh_49, \
                         kh_45, kh_46, kh_47, kh_48, kh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -5.0 * hh_45[k]
                  + f_0 * kh_45[k];

        t_46[k] = -5.0 * hh_46[k]
                  + f_0 * kh_46[k];

        t_47[k] = -5.0 * hh_47[k]
                  + f_0 * kh_47[k];

        t_48[k] = -5.0 * hh_48[k]
                  + f_0 * kh_48[k];

        t_49[k] = -5.0 * hh_49[k]
                  + f_0 * kh_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hh_50, hh_51, hh_52, hh_53, hh_54, \
                         kh_50, kh_51, kh_52, kh_53, kh_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -5.0 * hh_50[k]
                  + f_0 * kh_50[k];

        t_51[k] = -5.0 * hh_51[k]
                  + f_0 * kh_51[k];

        t_52[k] = -5.0 * hh_52[k]
                  + f_0 * kh_52[k];

        t_53[k] = -5.0 * hh_53[k]
                  + f_0 * kh_53[k];

        t_54[k] = -5.0 * hh_54[k]
                  + f_0 * kh_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, hh_55, hh_56, hh_57, hh_58, hh_59, \
                         kh_55, kh_56, kh_57, kh_58, kh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -5.0 * hh_55[k]
                  + f_0 * kh_55[k];

        t_56[k] = -5.0 * hh_56[k]
                  + f_0 * kh_56[k];

        t_57[k] = -5.0 * hh_57[k]
                  + f_0 * kh_57[k];

        t_58[k] = -5.0 * hh_58[k]
                  + f_0 * kh_58[k];

        t_59[k] = -5.0 * hh_59[k]
                  + f_0 * kh_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, hh_60, hh_61, hh_62, hh_63, hh_64, \
                         kh_60, kh_61, kh_62, kh_63, kh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -5.0 * hh_60[k]
                  + f_0 * kh_60[k];

        t_61[k] = -5.0 * hh_61[k]
                  + f_0 * kh_61[k];

        t_62[k] = -5.0 * hh_62[k]
                  + f_0 * kh_62[k];

        t_63[k] = -4.0 * hh_63[k]
                  + f_0 * kh_63[k];

        t_64[k] = -4.0 * hh_64[k]
                  + f_0 * kh_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, hh_65, hh_66, hh_67, hh_68, hh_69, \
                         kh_65, kh_66, kh_67, kh_68, kh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -4.0 * hh_65[k]
                  + f_0 * kh_65[k];

        t_66[k] = -4.0 * hh_66[k]
                  + f_0 * kh_66[k];

        t_67[k] = -4.0 * hh_67[k]
                  + f_0 * kh_67[k];

        t_68[k] = -4.0 * hh_68[k]
                  + f_0 * kh_68[k];

        t_69[k] = -4.0 * hh_69[k]
                  + f_0 * kh_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, hh_70, hh_71, hh_72, hh_73, hh_74, \
                         kh_70, kh_71, kh_72, kh_73, kh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -4.0 * hh_70[k]
                  + f_0 * kh_70[k];

        t_71[k] = -4.0 * hh_71[k]
                  + f_0 * kh_71[k];

        t_72[k] = -4.0 * hh_72[k]
                  + f_0 * kh_72[k];

        t_73[k] = -4.0 * hh_73[k]
                  + f_0 * kh_73[k];

        t_74[k] = -4.0 * hh_74[k]
                  + f_0 * kh_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, hh_75, hh_76, hh_77, hh_78, hh_79, \
                         kh_75, kh_76, kh_77, kh_78, kh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -4.0 * hh_75[k]
                  + f_0 * kh_75[k];

        t_76[k] = -4.0 * hh_76[k]
                  + f_0 * kh_76[k];

        t_77[k] = -4.0 * hh_77[k]
                  + f_0 * kh_77[k];

        t_78[k] = -4.0 * hh_78[k]
                  + f_0 * kh_78[k];

        t_79[k] = -4.0 * hh_79[k]
                  + f_0 * kh_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, hh_80, hh_81, hh_82, hh_83, hh_84, \
                         kh_80, kh_81, kh_82, kh_83, kh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -4.0 * hh_80[k]
                  + f_0 * kh_80[k];

        t_81[k] = -4.0 * hh_81[k]
                  + f_0 * kh_81[k];

        t_82[k] = -4.0 * hh_82[k]
                  + f_0 * kh_82[k];

        t_83[k] = -4.0 * hh_83[k]
                  + f_0 * kh_83[k];

        t_84[k] = -4.0 * hh_84[k]
                  + f_0 * kh_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, hh_85, hh_86, hh_87, hh_88, hh_89, \
                         kh_85, kh_86, kh_87, kh_88, kh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * hh_85[k]
                  + f_0 * kh_85[k];

        t_86[k] = -4.0 * hh_86[k]
                  + f_0 * kh_86[k];

        t_87[k] = -4.0 * hh_87[k]
                  + f_0 * kh_87[k];

        t_88[k] = -4.0 * hh_88[k]
                  + f_0 * kh_88[k];

        t_89[k] = -4.0 * hh_89[k]
                  + f_0 * kh_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, hh_90, hh_91, hh_92, hh_93, hh_94, \
                         kh_90, kh_91, kh_92, kh_93, kh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -4.0 * hh_90[k]
                  + f_0 * kh_90[k];

        t_91[k] = -4.0 * hh_91[k]
                  + f_0 * kh_91[k];

        t_92[k] = -4.0 * hh_92[k]
                  + f_0 * kh_92[k];

        t_93[k] = -4.0 * hh_93[k]
                  + f_0 * kh_93[k];

        t_94[k] = -4.0 * hh_94[k]
                  + f_0 * kh_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, hh_95, hh_96, hh_97, hh_98, hh_99, \
                         kh_95, kh_96, kh_97, kh_98, kh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -4.0 * hh_95[k]
                  + f_0 * kh_95[k];

        t_96[k] = -4.0 * hh_96[k]
                  + f_0 * kh_96[k];

        t_97[k] = -4.0 * hh_97[k]
                  + f_0 * kh_97[k];

        t_98[k] = -4.0 * hh_98[k]
                  + f_0 * kh_98[k];

        t_99[k] = -4.0 * hh_99[k]
                  + f_0 * kh_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, hh_100, hh_101, hh_102, hh_103, \
                         hh_104, kh_100, kh_101, kh_102, kh_103, \
                         kh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -4.0 * hh_100[k]
                   + f_0 * kh_100[k];

        t_101[k] = -4.0 * hh_101[k]
                   + f_0 * kh_101[k];

        t_102[k] = -4.0 * hh_102[k]
                   + f_0 * kh_102[k];

        t_103[k] = -4.0 * hh_103[k]
                   + f_0 * kh_103[k];

        t_104[k] = -4.0 * hh_104[k]
                   + f_0 * kh_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, hh_105, hh_106, hh_107, hh_108, \
                         hh_109, kh_105, kh_106, kh_107, kh_108, \
                         kh_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -4.0 * hh_105[k]
                   + f_0 * kh_105[k];

        t_106[k] = -4.0 * hh_106[k]
                   + f_0 * kh_106[k];

        t_107[k] = -4.0 * hh_107[k]
                   + f_0 * kh_107[k];

        t_108[k] = -4.0 * hh_108[k]
                   + f_0 * kh_108[k];

        t_109[k] = -4.0 * hh_109[k]
                   + f_0 * kh_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, hh_110, hh_111, hh_112, hh_113, \
                         hh_114, kh_110, kh_111, kh_112, kh_113, \
                         kh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -4.0 * hh_110[k]
                   + f_0 * kh_110[k];

        t_111[k] = -4.0 * hh_111[k]
                   + f_0 * kh_111[k];

        t_112[k] = -4.0 * hh_112[k]
                   + f_0 * kh_112[k];

        t_113[k] = -4.0 * hh_113[k]
                   + f_0 * kh_113[k];

        t_114[k] = -4.0 * hh_114[k]
                   + f_0 * kh_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, hh_115, hh_116, hh_117, hh_118, \
                         hh_119, kh_115, kh_116, kh_117, kh_118, \
                         kh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -4.0 * hh_115[k]
                   + f_0 * kh_115[k];

        t_116[k] = -4.0 * hh_116[k]
                   + f_0 * kh_116[k];

        t_117[k] = -4.0 * hh_117[k]
                   + f_0 * kh_117[k];

        t_118[k] = -4.0 * hh_118[k]
                   + f_0 * kh_118[k];

        t_119[k] = -4.0 * hh_119[k]
                   + f_0 * kh_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, hh_120, hh_121, hh_122, hh_123, \
                         hh_124, kh_120, kh_121, kh_122, kh_123, \
                         kh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -4.0 * hh_120[k]
                   + f_0 * kh_120[k];

        t_121[k] = -4.0 * hh_121[k]
                   + f_0 * kh_121[k];

        t_122[k] = -4.0 * hh_122[k]
                   + f_0 * kh_122[k];

        t_123[k] = -4.0 * hh_123[k]
                   + f_0 * kh_123[k];

        t_124[k] = -4.0 * hh_124[k]
                   + f_0 * kh_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, hh_125, hh_126, hh_127, hh_128, \
                         hh_129, kh_125, kh_126, kh_127, kh_128, \
                         kh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -4.0 * hh_125[k]
                   + f_0 * kh_125[k];

        t_126[k] = -3.0 * hh_126[k]
                   + f_0 * kh_126[k];

        t_127[k] = -3.0 * hh_127[k]
                   + f_0 * kh_127[k];

        t_128[k] = -3.0 * hh_128[k]
                   + f_0 * kh_128[k];

        t_129[k] = -3.0 * hh_129[k]
                   + f_0 * kh_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, hh_130, hh_131, hh_132, hh_133, \
                         hh_134, kh_130, kh_131, kh_132, kh_133, \
                         kh_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -3.0 * hh_130[k]
                   + f_0 * kh_130[k];

        t_131[k] = -3.0 * hh_131[k]
                   + f_0 * kh_131[k];

        t_132[k] = -3.0 * hh_132[k]
                   + f_0 * kh_132[k];

        t_133[k] = -3.0 * hh_133[k]
                   + f_0 * kh_133[k];

        t_134[k] = -3.0 * hh_134[k]
                   + f_0 * kh_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, hh_135, hh_136, hh_137, hh_138, \
                         hh_139, kh_135, kh_136, kh_137, kh_138, \
                         kh_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -3.0 * hh_135[k]
                   + f_0 * kh_135[k];

        t_136[k] = -3.0 * hh_136[k]
                   + f_0 * kh_136[k];

        t_137[k] = -3.0 * hh_137[k]
                   + f_0 * kh_137[k];

        t_138[k] = -3.0 * hh_138[k]
                   + f_0 * kh_138[k];

        t_139[k] = -3.0 * hh_139[k]
                   + f_0 * kh_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, hh_140, hh_141, hh_142, hh_143, \
                         hh_144, kh_140, kh_141, kh_142, kh_143, \
                         kh_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -3.0 * hh_140[k]
                   + f_0 * kh_140[k];

        t_141[k] = -3.0 * hh_141[k]
                   + f_0 * kh_141[k];

        t_142[k] = -3.0 * hh_142[k]
                   + f_0 * kh_142[k];

        t_143[k] = -3.0 * hh_143[k]
                   + f_0 * kh_143[k];

        t_144[k] = -3.0 * hh_144[k]
                   + f_0 * kh_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, hh_145, hh_146, hh_147, hh_148, \
                         hh_149, kh_145, kh_146, kh_147, kh_148, \
                         kh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -3.0 * hh_145[k]
                   + f_0 * kh_145[k];

        t_146[k] = -3.0 * hh_146[k]
                   + f_0 * kh_146[k];

        t_147[k] = -3.0 * hh_147[k]
                   + f_0 * kh_147[k];

        t_148[k] = -3.0 * hh_148[k]
                   + f_0 * kh_148[k];

        t_149[k] = -3.0 * hh_149[k]
                   + f_0 * kh_149[k];
    }
}

static auto
compute_prim_geom_10_ih_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hh, const size_t kh,
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

    const auto *hh_150 = buffer.data(hh + 150);
    const auto *hh_151 = buffer.data(hh + 151);
    const auto *hh_152 = buffer.data(hh + 152);
    const auto *hh_153 = buffer.data(hh + 153);
    const auto *hh_154 = buffer.data(hh + 154);
    const auto *hh_155 = buffer.data(hh + 155);
    const auto *hh_156 = buffer.data(hh + 156);
    const auto *hh_157 = buffer.data(hh + 157);
    const auto *hh_158 = buffer.data(hh + 158);
    const auto *hh_159 = buffer.data(hh + 159);
    const auto *hh_160 = buffer.data(hh + 160);
    const auto *hh_161 = buffer.data(hh + 161);
    const auto *hh_162 = buffer.data(hh + 162);
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_164 = buffer.data(hh + 164);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_166 = buffer.data(hh + 166);
    const auto *hh_167 = buffer.data(hh + 167);
    const auto *hh_168 = buffer.data(hh + 168);
    const auto *hh_169 = buffer.data(hh + 169);
    const auto *hh_170 = buffer.data(hh + 170);
    const auto *hh_171 = buffer.data(hh + 171);
    const auto *hh_172 = buffer.data(hh + 172);
    const auto *hh_173 = buffer.data(hh + 173);
    const auto *hh_174 = buffer.data(hh + 174);
    const auto *hh_175 = buffer.data(hh + 175);
    const auto *hh_176 = buffer.data(hh + 176);
    const auto *hh_177 = buffer.data(hh + 177);
    const auto *hh_178 = buffer.data(hh + 178);
    const auto *hh_179 = buffer.data(hh + 179);
    const auto *hh_180 = buffer.data(hh + 180);
    const auto *hh_181 = buffer.data(hh + 181);
    const auto *hh_182 = buffer.data(hh + 182);
    const auto *hh_183 = buffer.data(hh + 183);
    const auto *hh_184 = buffer.data(hh + 184);
    const auto *hh_185 = buffer.data(hh + 185);
    const auto *hh_186 = buffer.data(hh + 186);
    const auto *hh_187 = buffer.data(hh + 187);
    const auto *hh_188 = buffer.data(hh + 188);
    const auto *hh_189 = buffer.data(hh + 189);
    const auto *hh_190 = buffer.data(hh + 190);
    const auto *hh_191 = buffer.data(hh + 191);
    const auto *hh_192 = buffer.data(hh + 192);
    const auto *hh_193 = buffer.data(hh + 193);
    const auto *hh_194 = buffer.data(hh + 194);
    const auto *hh_195 = buffer.data(hh + 195);
    const auto *hh_196 = buffer.data(hh + 196);
    const auto *hh_197 = buffer.data(hh + 197);
    const auto *hh_198 = buffer.data(hh + 198);
    const auto *hh_199 = buffer.data(hh + 199);
    const auto *hh_200 = buffer.data(hh + 200);
    const auto *hh_201 = buffer.data(hh + 201);
    const auto *hh_202 = buffer.data(hh + 202);
    const auto *hh_203 = buffer.data(hh + 203);
    const auto *hh_204 = buffer.data(hh + 204);
    const auto *hh_205 = buffer.data(hh + 205);
    const auto *hh_206 = buffer.data(hh + 206);
    const auto *hh_207 = buffer.data(hh + 207);
    const auto *hh_208 = buffer.data(hh + 208);
    const auto *hh_209 = buffer.data(hh + 209);
    const auto *hh_210 = buffer.data(hh + 210);
    const auto *hh_211 = buffer.data(hh + 211);
    const auto *hh_212 = buffer.data(hh + 212);
    const auto *hh_213 = buffer.data(hh + 213);
    const auto *hh_214 = buffer.data(hh + 214);
    const auto *hh_215 = buffer.data(hh + 215);
    const auto *hh_216 = buffer.data(hh + 216);
    const auto *hh_217 = buffer.data(hh + 217);
    const auto *hh_218 = buffer.data(hh + 218);
    const auto *hh_219 = buffer.data(hh + 219);
    const auto *hh_220 = buffer.data(hh + 220);
    const auto *hh_221 = buffer.data(hh + 221);
    const auto *hh_222 = buffer.data(hh + 222);
    const auto *hh_223 = buffer.data(hh + 223);
    const auto *hh_224 = buffer.data(hh + 224);
    const auto *hh_225 = buffer.data(hh + 225);
    const auto *hh_226 = buffer.data(hh + 226);
    const auto *hh_227 = buffer.data(hh + 227);
    const auto *hh_228 = buffer.data(hh + 228);
    const auto *hh_229 = buffer.data(hh + 229);
    const auto *hh_230 = buffer.data(hh + 230);
    const auto *hh_231 = buffer.data(hh + 231);
    const auto *hh_232 = buffer.data(hh + 232);
    const auto *hh_233 = buffer.data(hh + 233);
    const auto *hh_234 = buffer.data(hh + 234);
    const auto *hh_235 = buffer.data(hh + 235);
    const auto *hh_236 = buffer.data(hh + 236);
    const auto *hh_237 = buffer.data(hh + 237);
    const auto *hh_238 = buffer.data(hh + 238);
    const auto *hh_239 = buffer.data(hh + 239);
    const auto *hh_240 = buffer.data(hh + 240);
    const auto *hh_241 = buffer.data(hh + 241);
    const auto *hh_242 = buffer.data(hh + 242);
    const auto *hh_243 = buffer.data(hh + 243);
    const auto *hh_244 = buffer.data(hh + 244);
    const auto *hh_245 = buffer.data(hh + 245);
    const auto *hh_246 = buffer.data(hh + 246);
    const auto *hh_247 = buffer.data(hh + 247);
    const auto *hh_248 = buffer.data(hh + 248);
    const auto *hh_249 = buffer.data(hh + 249);
    const auto *hh_250 = buffer.data(hh + 250);
    const auto *hh_251 = buffer.data(hh + 251);
    const auto *hh_252 = buffer.data(hh + 252);
    const auto *hh_253 = buffer.data(hh + 253);
    const auto *hh_254 = buffer.data(hh + 254);
    const auto *hh_255 = buffer.data(hh + 255);
    const auto *hh_256 = buffer.data(hh + 256);
    const auto *hh_257 = buffer.data(hh + 257);
    const auto *hh_258 = buffer.data(hh + 258);
    const auto *hh_259 = buffer.data(hh + 259);
    const auto *hh_260 = buffer.data(hh + 260);
    const auto *hh_261 = buffer.data(hh + 261);
    const auto *hh_262 = buffer.data(hh + 262);
    const auto *hh_263 = buffer.data(hh + 263);
    const auto *hh_264 = buffer.data(hh + 264);
    const auto *hh_265 = buffer.data(hh + 265);
    const auto *hh_266 = buffer.data(hh + 266);
    const auto *hh_267 = buffer.data(hh + 267);
    const auto *hh_268 = buffer.data(hh + 268);
    const auto *hh_269 = buffer.data(hh + 269);
    const auto *hh_270 = buffer.data(hh + 270);
    const auto *hh_271 = buffer.data(hh + 271);
    const auto *hh_272 = buffer.data(hh + 272);
    const auto *hh_273 = buffer.data(hh + 273);
    const auto *hh_274 = buffer.data(hh + 274);
    const auto *hh_275 = buffer.data(hh + 275);
    const auto *hh_276 = buffer.data(hh + 276);
    const auto *hh_277 = buffer.data(hh + 277);
    const auto *hh_278 = buffer.data(hh + 278);
    const auto *hh_279 = buffer.data(hh + 279);
    const auto *hh_280 = buffer.data(hh + 280);
    const auto *hh_281 = buffer.data(hh + 281);
    const auto *hh_282 = buffer.data(hh + 282);
    const auto *hh_283 = buffer.data(hh + 283);
    const auto *hh_284 = buffer.data(hh + 284);
    const auto *hh_285 = buffer.data(hh + 285);
    const auto *hh_286 = buffer.data(hh + 286);
    const auto *hh_287 = buffer.data(hh + 287);
    const auto *hh_288 = buffer.data(hh + 288);
    const auto *hh_289 = buffer.data(hh + 289);
    const auto *hh_290 = buffer.data(hh + 290);
    const auto *hh_291 = buffer.data(hh + 291);
    const auto *hh_292 = buffer.data(hh + 292);
    const auto *hh_293 = buffer.data(hh + 293);
    const auto *hh_294 = buffer.data(hh + 294);
    const auto *hh_295 = buffer.data(hh + 295);
    const auto *hh_296 = buffer.data(hh + 296);
    const auto *hh_297 = buffer.data(hh + 297);
    const auto *hh_298 = buffer.data(hh + 298);
    const auto *hh_299 = buffer.data(hh + 299);

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

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, hh_150, hh_151, hh_152, hh_153, \
                         hh_154, kh_150, kh_151, kh_152, kh_153, \
                         kh_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -3.0 * hh_150[k]
                   + f_0 * kh_150[k];

        t_151[k] = -3.0 * hh_151[k]
                   + f_0 * kh_151[k];

        t_152[k] = -3.0 * hh_152[k]
                   + f_0 * kh_152[k];

        t_153[k] = -3.0 * hh_153[k]
                   + f_0 * kh_153[k];

        t_154[k] = -3.0 * hh_154[k]
                   + f_0 * kh_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, hh_155, hh_156, hh_157, hh_158, \
                         hh_159, kh_155, kh_156, kh_157, kh_158, \
                         kh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -3.0 * hh_155[k]
                   + f_0 * kh_155[k];

        t_156[k] = -3.0 * hh_156[k]
                   + f_0 * kh_156[k];

        t_157[k] = -3.0 * hh_157[k]
                   + f_0 * kh_157[k];

        t_158[k] = -3.0 * hh_158[k]
                   + f_0 * kh_158[k];

        t_159[k] = -3.0 * hh_159[k]
                   + f_0 * kh_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, hh_160, hh_161, hh_162, hh_163, \
                         hh_164, kh_160, kh_161, kh_162, kh_163, \
                         kh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -3.0 * hh_160[k]
                   + f_0 * kh_160[k];

        t_161[k] = -3.0 * hh_161[k]
                   + f_0 * kh_161[k];

        t_162[k] = -3.0 * hh_162[k]
                   + f_0 * kh_162[k];

        t_163[k] = -3.0 * hh_163[k]
                   + f_0 * kh_163[k];

        t_164[k] = -3.0 * hh_164[k]
                   + f_0 * kh_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, hh_165, hh_166, hh_167, hh_168, \
                         hh_169, kh_165, kh_166, kh_167, kh_168, \
                         kh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -3.0 * hh_165[k]
                   + f_0 * kh_165[k];

        t_166[k] = -3.0 * hh_166[k]
                   + f_0 * kh_166[k];

        t_167[k] = -3.0 * hh_167[k]
                   + f_0 * kh_167[k];

        t_168[k] = -3.0 * hh_168[k]
                   + f_0 * kh_168[k];

        t_169[k] = -3.0 * hh_169[k]
                   + f_0 * kh_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, hh_170, hh_171, hh_172, hh_173, \
                         hh_174, kh_170, kh_171, kh_172, kh_173, \
                         kh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * hh_170[k]
                   + f_0 * kh_170[k];

        t_171[k] = -3.0 * hh_171[k]
                   + f_0 * kh_171[k];

        t_172[k] = -3.0 * hh_172[k]
                   + f_0 * kh_172[k];

        t_173[k] = -3.0 * hh_173[k]
                   + f_0 * kh_173[k];

        t_174[k] = -3.0 * hh_174[k]
                   + f_0 * kh_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, hh_175, hh_176, hh_177, hh_178, \
                         hh_179, kh_175, kh_176, kh_177, kh_178, \
                         kh_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * hh_175[k]
                   + f_0 * kh_175[k];

        t_176[k] = -3.0 * hh_176[k]
                   + f_0 * kh_176[k];

        t_177[k] = -3.0 * hh_177[k]
                   + f_0 * kh_177[k];

        t_178[k] = -3.0 * hh_178[k]
                   + f_0 * kh_178[k];

        t_179[k] = -3.0 * hh_179[k]
                   + f_0 * kh_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, hh_180, hh_181, hh_182, hh_183, \
                         hh_184, kh_180, kh_181, kh_182, kh_183, \
                         kh_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * hh_180[k]
                   + f_0 * kh_180[k];

        t_181[k] = -3.0 * hh_181[k]
                   + f_0 * kh_181[k];

        t_182[k] = -3.0 * hh_182[k]
                   + f_0 * kh_182[k];

        t_183[k] = -3.0 * hh_183[k]
                   + f_0 * kh_183[k];

        t_184[k] = -3.0 * hh_184[k]
                   + f_0 * kh_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, hh_185, hh_186, hh_187, hh_188, \
                         hh_189, kh_185, kh_186, kh_187, kh_188, \
                         kh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * hh_185[k]
                   + f_0 * kh_185[k];

        t_186[k] = -3.0 * hh_186[k]
                   + f_0 * kh_186[k];

        t_187[k] = -3.0 * hh_187[k]
                   + f_0 * kh_187[k];

        t_188[k] = -3.0 * hh_188[k]
                   + f_0 * kh_188[k];

        t_189[k] = -3.0 * hh_189[k]
                   + f_0 * kh_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, hh_190, hh_191, hh_192, hh_193, \
                         hh_194, kh_190, kh_191, kh_192, kh_193, \
                         kh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * hh_190[k]
                   + f_0 * kh_190[k];

        t_191[k] = -3.0 * hh_191[k]
                   + f_0 * kh_191[k];

        t_192[k] = -3.0 * hh_192[k]
                   + f_0 * kh_192[k];

        t_193[k] = -3.0 * hh_193[k]
                   + f_0 * kh_193[k];

        t_194[k] = -3.0 * hh_194[k]
                   + f_0 * kh_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, hh_195, hh_196, hh_197, hh_198, \
                         hh_199, kh_195, kh_196, kh_197, kh_198, \
                         kh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * hh_195[k]
                   + f_0 * kh_195[k];

        t_196[k] = -3.0 * hh_196[k]
                   + f_0 * kh_196[k];

        t_197[k] = -3.0 * hh_197[k]
                   + f_0 * kh_197[k];

        t_198[k] = -3.0 * hh_198[k]
                   + f_0 * kh_198[k];

        t_199[k] = -3.0 * hh_199[k]
                   + f_0 * kh_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, hh_200, hh_201, hh_202, hh_203, \
                         hh_204, kh_200, kh_201, kh_202, kh_203, \
                         kh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -3.0 * hh_200[k]
                   + f_0 * kh_200[k];

        t_201[k] = -3.0 * hh_201[k]
                   + f_0 * kh_201[k];

        t_202[k] = -3.0 * hh_202[k]
                   + f_0 * kh_202[k];

        t_203[k] = -3.0 * hh_203[k]
                   + f_0 * kh_203[k];

        t_204[k] = -3.0 * hh_204[k]
                   + f_0 * kh_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, hh_205, hh_206, hh_207, hh_208, \
                         hh_209, kh_205, kh_206, kh_207, kh_208, \
                         kh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -3.0 * hh_205[k]
                   + f_0 * kh_205[k];

        t_206[k] = -3.0 * hh_206[k]
                   + f_0 * kh_206[k];

        t_207[k] = -3.0 * hh_207[k]
                   + f_0 * kh_207[k];

        t_208[k] = -3.0 * hh_208[k]
                   + f_0 * kh_208[k];

        t_209[k] = -3.0 * hh_209[k]
                   + f_0 * kh_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, hh_210, hh_211, hh_212, hh_213, \
                         hh_214, kh_210, kh_211, kh_212, kh_213, \
                         kh_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * hh_210[k]
                   + f_0 * kh_210[k];

        t_211[k] = -2.0 * hh_211[k]
                   + f_0 * kh_211[k];

        t_212[k] = -2.0 * hh_212[k]
                   + f_0 * kh_212[k];

        t_213[k] = -2.0 * hh_213[k]
                   + f_0 * kh_213[k];

        t_214[k] = -2.0 * hh_214[k]
                   + f_0 * kh_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, hh_215, hh_216, hh_217, hh_218, \
                         hh_219, kh_215, kh_216, kh_217, kh_218, \
                         kh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * hh_215[k]
                   + f_0 * kh_215[k];

        t_216[k] = -2.0 * hh_216[k]
                   + f_0 * kh_216[k];

        t_217[k] = -2.0 * hh_217[k]
                   + f_0 * kh_217[k];

        t_218[k] = -2.0 * hh_218[k]
                   + f_0 * kh_218[k];

        t_219[k] = -2.0 * hh_219[k]
                   + f_0 * kh_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, hh_220, hh_221, hh_222, hh_223, \
                         hh_224, kh_220, kh_221, kh_222, kh_223, \
                         kh_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * hh_220[k]
                   + f_0 * kh_220[k];

        t_221[k] = -2.0 * hh_221[k]
                   + f_0 * kh_221[k];

        t_222[k] = -2.0 * hh_222[k]
                   + f_0 * kh_222[k];

        t_223[k] = -2.0 * hh_223[k]
                   + f_0 * kh_223[k];

        t_224[k] = -2.0 * hh_224[k]
                   + f_0 * kh_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, hh_225, hh_226, hh_227, hh_228, \
                         hh_229, kh_225, kh_226, kh_227, kh_228, \
                         kh_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -2.0 * hh_225[k]
                   + f_0 * kh_225[k];

        t_226[k] = -2.0 * hh_226[k]
                   + f_0 * kh_226[k];

        t_227[k] = -2.0 * hh_227[k]
                   + f_0 * kh_227[k];

        t_228[k] = -2.0 * hh_228[k]
                   + f_0 * kh_228[k];

        t_229[k] = -2.0 * hh_229[k]
                   + f_0 * kh_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, hh_230, hh_231, hh_232, hh_233, \
                         hh_234, kh_230, kh_231, kh_232, kh_233, \
                         kh_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -2.0 * hh_230[k]
                   + f_0 * kh_230[k];

        t_231[k] = -2.0 * hh_231[k]
                   + f_0 * kh_231[k];

        t_232[k] = -2.0 * hh_232[k]
                   + f_0 * kh_232[k];

        t_233[k] = -2.0 * hh_233[k]
                   + f_0 * kh_233[k];

        t_234[k] = -2.0 * hh_234[k]
                   + f_0 * kh_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, hh_235, hh_236, hh_237, hh_238, \
                         hh_239, kh_235, kh_236, kh_237, kh_238, \
                         kh_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -2.0 * hh_235[k]
                   + f_0 * kh_235[k];

        t_236[k] = -2.0 * hh_236[k]
                   + f_0 * kh_236[k];

        t_237[k] = -2.0 * hh_237[k]
                   + f_0 * kh_237[k];

        t_238[k] = -2.0 * hh_238[k]
                   + f_0 * kh_238[k];

        t_239[k] = -2.0 * hh_239[k]
                   + f_0 * kh_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, hh_240, hh_241, hh_242, hh_243, \
                         hh_244, kh_240, kh_241, kh_242, kh_243, \
                         kh_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -2.0 * hh_240[k]
                   + f_0 * kh_240[k];

        t_241[k] = -2.0 * hh_241[k]
                   + f_0 * kh_241[k];

        t_242[k] = -2.0 * hh_242[k]
                   + f_0 * kh_242[k];

        t_243[k] = -2.0 * hh_243[k]
                   + f_0 * kh_243[k];

        t_244[k] = -2.0 * hh_244[k]
                   + f_0 * kh_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, hh_245, hh_246, hh_247, hh_248, \
                         hh_249, kh_245, kh_246, kh_247, kh_248, \
                         kh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -2.0 * hh_245[k]
                   + f_0 * kh_245[k];

        t_246[k] = -2.0 * hh_246[k]
                   + f_0 * kh_246[k];

        t_247[k] = -2.0 * hh_247[k]
                   + f_0 * kh_247[k];

        t_248[k] = -2.0 * hh_248[k]
                   + f_0 * kh_248[k];

        t_249[k] = -2.0 * hh_249[k]
                   + f_0 * kh_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, hh_250, hh_251, hh_252, hh_253, \
                         hh_254, kh_250, kh_251, kh_252, kh_253, \
                         kh_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * hh_250[k]
                   + f_0 * kh_250[k];

        t_251[k] = -2.0 * hh_251[k]
                   + f_0 * kh_251[k];

        t_252[k] = -2.0 * hh_252[k]
                   + f_0 * kh_252[k];

        t_253[k] = -2.0 * hh_253[k]
                   + f_0 * kh_253[k];

        t_254[k] = -2.0 * hh_254[k]
                   + f_0 * kh_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, hh_255, hh_256, hh_257, hh_258, \
                         hh_259, kh_255, kh_256, kh_257, kh_258, \
                         kh_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * hh_255[k]
                   + f_0 * kh_255[k];

        t_256[k] = -2.0 * hh_256[k]
                   + f_0 * kh_256[k];

        t_257[k] = -2.0 * hh_257[k]
                   + f_0 * kh_257[k];

        t_258[k] = -2.0 * hh_258[k]
                   + f_0 * kh_258[k];

        t_259[k] = -2.0 * hh_259[k]
                   + f_0 * kh_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, hh_260, hh_261, hh_262, hh_263, \
                         hh_264, kh_260, kh_261, kh_262, kh_263, \
                         kh_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * hh_260[k]
                   + f_0 * kh_260[k];

        t_261[k] = -2.0 * hh_261[k]
                   + f_0 * kh_261[k];

        t_262[k] = -2.0 * hh_262[k]
                   + f_0 * kh_262[k];

        t_263[k] = -2.0 * hh_263[k]
                   + f_0 * kh_263[k];

        t_264[k] = -2.0 * hh_264[k]
                   + f_0 * kh_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, hh_265, hh_266, hh_267, hh_268, \
                         hh_269, kh_265, kh_266, kh_267, kh_268, \
                         kh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * hh_265[k]
                   + f_0 * kh_265[k];

        t_266[k] = -2.0 * hh_266[k]
                   + f_0 * kh_266[k];

        t_267[k] = -2.0 * hh_267[k]
                   + f_0 * kh_267[k];

        t_268[k] = -2.0 * hh_268[k]
                   + f_0 * kh_268[k];

        t_269[k] = -2.0 * hh_269[k]
                   + f_0 * kh_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, hh_270, hh_271, hh_272, hh_273, \
                         hh_274, kh_270, kh_271, kh_272, kh_273, \
                         kh_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * hh_270[k]
                   + f_0 * kh_270[k];

        t_271[k] = -2.0 * hh_271[k]
                   + f_0 * kh_271[k];

        t_272[k] = -2.0 * hh_272[k]
                   + f_0 * kh_272[k];

        t_273[k] = -2.0 * hh_273[k]
                   + f_0 * kh_273[k];

        t_274[k] = -2.0 * hh_274[k]
                   + f_0 * kh_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, hh_275, hh_276, hh_277, hh_278, \
                         hh_279, kh_275, kh_276, kh_277, kh_278, \
                         kh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * hh_275[k]
                   + f_0 * kh_275[k];

        t_276[k] = -2.0 * hh_276[k]
                   + f_0 * kh_276[k];

        t_277[k] = -2.0 * hh_277[k]
                   + f_0 * kh_277[k];

        t_278[k] = -2.0 * hh_278[k]
                   + f_0 * kh_278[k];

        t_279[k] = -2.0 * hh_279[k]
                   + f_0 * kh_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, hh_280, hh_281, hh_282, hh_283, \
                         hh_284, kh_280, kh_281, kh_282, kh_283, \
                         kh_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * hh_280[k]
                   + f_0 * kh_280[k];

        t_281[k] = -2.0 * hh_281[k]
                   + f_0 * kh_281[k];

        t_282[k] = -2.0 * hh_282[k]
                   + f_0 * kh_282[k];

        t_283[k] = -2.0 * hh_283[k]
                   + f_0 * kh_283[k];

        t_284[k] = -2.0 * hh_284[k]
                   + f_0 * kh_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, hh_285, hh_286, hh_287, hh_288, \
                         hh_289, kh_285, kh_286, kh_287, kh_288, \
                         kh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -2.0 * hh_285[k]
                   + f_0 * kh_285[k];

        t_286[k] = -2.0 * hh_286[k]
                   + f_0 * kh_286[k];

        t_287[k] = -2.0 * hh_287[k]
                   + f_0 * kh_287[k];

        t_288[k] = -2.0 * hh_288[k]
                   + f_0 * kh_288[k];

        t_289[k] = -2.0 * hh_289[k]
                   + f_0 * kh_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, hh_290, hh_291, hh_292, hh_293, \
                         hh_294, kh_290, kh_291, kh_292, kh_293, \
                         kh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -2.0 * hh_290[k]
                   + f_0 * kh_290[k];

        t_291[k] = -2.0 * hh_291[k]
                   + f_0 * kh_291[k];

        t_292[k] = -2.0 * hh_292[k]
                   + f_0 * kh_292[k];

        t_293[k] = -2.0 * hh_293[k]
                   + f_0 * kh_293[k];

        t_294[k] = -2.0 * hh_294[k]
                   + f_0 * kh_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, hh_295, hh_296, hh_297, hh_298, \
                         hh_299, kh_295, kh_296, kh_297, kh_298, \
                         kh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -2.0 * hh_295[k]
                   + f_0 * kh_295[k];

        t_296[k] = -2.0 * hh_296[k]
                   + f_0 * kh_296[k];

        t_297[k] = -2.0 * hh_297[k]
                   + f_0 * kh_297[k];

        t_298[k] = -2.0 * hh_298[k]
                   + f_0 * kh_298[k];

        t_299[k] = -2.0 * hh_299[k]
                   + f_0 * kh_299[k];
    }
}

static auto
compute_prim_geom_10_ih_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hh, const size_t kh,
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
    auto *t_450 = buffer.data(target + 450);
    auto *t_451 = buffer.data(target + 451);
    auto *t_452 = buffer.data(target + 452);
    auto *t_453 = buffer.data(target + 453);
    auto *t_454 = buffer.data(target + 454);

    const auto *hh_300 = buffer.data(hh + 300);
    const auto *hh_301 = buffer.data(hh + 301);
    const auto *hh_302 = buffer.data(hh + 302);
    const auto *hh_303 = buffer.data(hh + 303);
    const auto *hh_304 = buffer.data(hh + 304);
    const auto *hh_305 = buffer.data(hh + 305);
    const auto *hh_306 = buffer.data(hh + 306);
    const auto *hh_307 = buffer.data(hh + 307);
    const auto *hh_308 = buffer.data(hh + 308);
    const auto *hh_309 = buffer.data(hh + 309);
    const auto *hh_310 = buffer.data(hh + 310);
    const auto *hh_311 = buffer.data(hh + 311);
    const auto *hh_312 = buffer.data(hh + 312);
    const auto *hh_313 = buffer.data(hh + 313);
    const auto *hh_314 = buffer.data(hh + 314);
    const auto *hh_315 = buffer.data(hh + 315);
    const auto *hh_316 = buffer.data(hh + 316);
    const auto *hh_317 = buffer.data(hh + 317);
    const auto *hh_318 = buffer.data(hh + 318);
    const auto *hh_319 = buffer.data(hh + 319);
    const auto *hh_320 = buffer.data(hh + 320);
    const auto *hh_321 = buffer.data(hh + 321);
    const auto *hh_322 = buffer.data(hh + 322);
    const auto *hh_323 = buffer.data(hh + 323);
    const auto *hh_324 = buffer.data(hh + 324);
    const auto *hh_325 = buffer.data(hh + 325);
    const auto *hh_326 = buffer.data(hh + 326);
    const auto *hh_327 = buffer.data(hh + 327);
    const auto *hh_328 = buffer.data(hh + 328);
    const auto *hh_329 = buffer.data(hh + 329);
    const auto *hh_330 = buffer.data(hh + 330);
    const auto *hh_331 = buffer.data(hh + 331);
    const auto *hh_332 = buffer.data(hh + 332);
    const auto *hh_333 = buffer.data(hh + 333);
    const auto *hh_334 = buffer.data(hh + 334);
    const auto *hh_335 = buffer.data(hh + 335);
    const auto *hh_336 = buffer.data(hh + 336);
    const auto *hh_337 = buffer.data(hh + 337);
    const auto *hh_338 = buffer.data(hh + 338);
    const auto *hh_339 = buffer.data(hh + 339);
    const auto *hh_340 = buffer.data(hh + 340);
    const auto *hh_341 = buffer.data(hh + 341);
    const auto *hh_342 = buffer.data(hh + 342);
    const auto *hh_343 = buffer.data(hh + 343);
    const auto *hh_344 = buffer.data(hh + 344);
    const auto *hh_345 = buffer.data(hh + 345);
    const auto *hh_346 = buffer.data(hh + 346);
    const auto *hh_347 = buffer.data(hh + 347);
    const auto *hh_348 = buffer.data(hh + 348);
    const auto *hh_349 = buffer.data(hh + 349);
    const auto *hh_350 = buffer.data(hh + 350);
    const auto *hh_351 = buffer.data(hh + 351);
    const auto *hh_352 = buffer.data(hh + 352);
    const auto *hh_353 = buffer.data(hh + 353);
    const auto *hh_354 = buffer.data(hh + 354);
    const auto *hh_355 = buffer.data(hh + 355);
    const auto *hh_356 = buffer.data(hh + 356);
    const auto *hh_357 = buffer.data(hh + 357);
    const auto *hh_358 = buffer.data(hh + 358);
    const auto *hh_359 = buffer.data(hh + 359);
    const auto *hh_360 = buffer.data(hh + 360);
    const auto *hh_361 = buffer.data(hh + 361);
    const auto *hh_362 = buffer.data(hh + 362);
    const auto *hh_363 = buffer.data(hh + 363);
    const auto *hh_364 = buffer.data(hh + 364);
    const auto *hh_365 = buffer.data(hh + 365);
    const auto *hh_366 = buffer.data(hh + 366);
    const auto *hh_367 = buffer.data(hh + 367);
    const auto *hh_368 = buffer.data(hh + 368);
    const auto *hh_369 = buffer.data(hh + 369);
    const auto *hh_370 = buffer.data(hh + 370);
    const auto *hh_371 = buffer.data(hh + 371);
    const auto *hh_372 = buffer.data(hh + 372);
    const auto *hh_373 = buffer.data(hh + 373);
    const auto *hh_374 = buffer.data(hh + 374);
    const auto *hh_375 = buffer.data(hh + 375);
    const auto *hh_376 = buffer.data(hh + 376);
    const auto *hh_377 = buffer.data(hh + 377);
    const auto *hh_378 = buffer.data(hh + 378);
    const auto *hh_379 = buffer.data(hh + 379);
    const auto *hh_380 = buffer.data(hh + 380);
    const auto *hh_381 = buffer.data(hh + 381);
    const auto *hh_382 = buffer.data(hh + 382);
    const auto *hh_383 = buffer.data(hh + 383);
    const auto *hh_384 = buffer.data(hh + 384);
    const auto *hh_385 = buffer.data(hh + 385);
    const auto *hh_386 = buffer.data(hh + 386);
    const auto *hh_387 = buffer.data(hh + 387);
    const auto *hh_388 = buffer.data(hh + 388);
    const auto *hh_389 = buffer.data(hh + 389);
    const auto *hh_390 = buffer.data(hh + 390);
    const auto *hh_391 = buffer.data(hh + 391);
    const auto *hh_392 = buffer.data(hh + 392);
    const auto *hh_393 = buffer.data(hh + 393);
    const auto *hh_394 = buffer.data(hh + 394);
    const auto *hh_395 = buffer.data(hh + 395);
    const auto *hh_396 = buffer.data(hh + 396);
    const auto *hh_397 = buffer.data(hh + 397);
    const auto *hh_398 = buffer.data(hh + 398);
    const auto *hh_399 = buffer.data(hh + 399);
    const auto *hh_400 = buffer.data(hh + 400);
    const auto *hh_401 = buffer.data(hh + 401);
    const auto *hh_402 = buffer.data(hh + 402);
    const auto *hh_403 = buffer.data(hh + 403);
    const auto *hh_404 = buffer.data(hh + 404);
    const auto *hh_405 = buffer.data(hh + 405);
    const auto *hh_406 = buffer.data(hh + 406);
    const auto *hh_407 = buffer.data(hh + 407);
    const auto *hh_408 = buffer.data(hh + 408);
    const auto *hh_409 = buffer.data(hh + 409);
    const auto *hh_410 = buffer.data(hh + 410);
    const auto *hh_411 = buffer.data(hh + 411);
    const auto *hh_412 = buffer.data(hh + 412);
    const auto *hh_413 = buffer.data(hh + 413);
    const auto *hh_414 = buffer.data(hh + 414);
    const auto *hh_415 = buffer.data(hh + 415);
    const auto *hh_416 = buffer.data(hh + 416);
    const auto *hh_417 = buffer.data(hh + 417);
    const auto *hh_418 = buffer.data(hh + 418);
    const auto *hh_419 = buffer.data(hh + 419);
    const auto *hh_420 = buffer.data(hh + 420);
    const auto *hh_421 = buffer.data(hh + 421);
    const auto *hh_422 = buffer.data(hh + 422);
    const auto *hh_423 = buffer.data(hh + 423);
    const auto *hh_424 = buffer.data(hh + 424);
    const auto *hh_425 = buffer.data(hh + 425);
    const auto *hh_426 = buffer.data(hh + 426);
    const auto *hh_427 = buffer.data(hh + 427);
    const auto *hh_428 = buffer.data(hh + 428);
    const auto *hh_429 = buffer.data(hh + 429);
    const auto *hh_430 = buffer.data(hh + 430);
    const auto *hh_431 = buffer.data(hh + 431);
    const auto *hh_432 = buffer.data(hh + 432);
    const auto *hh_433 = buffer.data(hh + 433);
    const auto *hh_434 = buffer.data(hh + 434);
    const auto *hh_435 = buffer.data(hh + 435);
    const auto *hh_436 = buffer.data(hh + 436);
    const auto *hh_437 = buffer.data(hh + 437);
    const auto *hh_438 = buffer.data(hh + 438);
    const auto *hh_439 = buffer.data(hh + 439);
    const auto *hh_440 = buffer.data(hh + 440);

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
    const auto *kh_450 = buffer.data(kh + 450);
    const auto *kh_451 = buffer.data(kh + 451);
    const auto *kh_452 = buffer.data(kh + 452);
    const auto *kh_453 = buffer.data(kh + 453);
    const auto *kh_454 = buffer.data(kh + 454);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, hh_300, hh_301, hh_302, hh_303, \
                         hh_304, kh_300, kh_301, kh_302, kh_303, \
                         kh_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -2.0 * hh_300[k]
                   + f_0 * kh_300[k];

        t_301[k] = -2.0 * hh_301[k]
                   + f_0 * kh_301[k];

        t_302[k] = -2.0 * hh_302[k]
                   + f_0 * kh_302[k];

        t_303[k] = -2.0 * hh_303[k]
                   + f_0 * kh_303[k];

        t_304[k] = -2.0 * hh_304[k]
                   + f_0 * kh_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, hh_305, hh_306, hh_307, hh_308, \
                         hh_309, kh_305, kh_306, kh_307, kh_308, \
                         kh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -2.0 * hh_305[k]
                   + f_0 * kh_305[k];

        t_306[k] = -2.0 * hh_306[k]
                   + f_0 * kh_306[k];

        t_307[k] = -2.0 * hh_307[k]
                   + f_0 * kh_307[k];

        t_308[k] = -2.0 * hh_308[k]
                   + f_0 * kh_308[k];

        t_309[k] = -2.0 * hh_309[k]
                   + f_0 * kh_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, hh_310, hh_311, hh_312, hh_313, \
                         hh_314, kh_310, kh_311, kh_312, kh_313, \
                         kh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -2.0 * hh_310[k]
                   + f_0 * kh_310[k];

        t_311[k] = -2.0 * hh_311[k]
                   + f_0 * kh_311[k];

        t_312[k] = -2.0 * hh_312[k]
                   + f_0 * kh_312[k];

        t_313[k] = -2.0 * hh_313[k]
                   + f_0 * kh_313[k];

        t_314[k] = -2.0 * hh_314[k]
                   + f_0 * kh_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, hh_315, hh_316, hh_317, hh_318, \
                         hh_319, kh_315, kh_316, kh_317, kh_318, \
                         kh_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -hh_315[k]
                   + f_0 * kh_315[k];

        t_316[k] = -hh_316[k]
                   + f_0 * kh_316[k];

        t_317[k] = -hh_317[k]
                   + f_0 * kh_317[k];

        t_318[k] = -hh_318[k]
                   + f_0 * kh_318[k];

        t_319[k] = -hh_319[k]
                   + f_0 * kh_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, hh_320, hh_321, hh_322, hh_323, \
                         hh_324, kh_320, kh_321, kh_322, kh_323, \
                         kh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -hh_320[k]
                   + f_0 * kh_320[k];

        t_321[k] = -hh_321[k]
                   + f_0 * kh_321[k];

        t_322[k] = -hh_322[k]
                   + f_0 * kh_322[k];

        t_323[k] = -hh_323[k]
                   + f_0 * kh_323[k];

        t_324[k] = -hh_324[k]
                   + f_0 * kh_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, hh_325, hh_326, hh_327, hh_328, \
                         hh_329, kh_325, kh_326, kh_327, kh_328, \
                         kh_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -hh_325[k]
                   + f_0 * kh_325[k];

        t_326[k] = -hh_326[k]
                   + f_0 * kh_326[k];

        t_327[k] = -hh_327[k]
                   + f_0 * kh_327[k];

        t_328[k] = -hh_328[k]
                   + f_0 * kh_328[k];

        t_329[k] = -hh_329[k]
                   + f_0 * kh_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, hh_330, hh_331, hh_332, hh_333, \
                         hh_334, kh_330, kh_331, kh_332, kh_333, \
                         kh_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -hh_330[k]
                   + f_0 * kh_330[k];

        t_331[k] = -hh_331[k]
                   + f_0 * kh_331[k];

        t_332[k] = -hh_332[k]
                   + f_0 * kh_332[k];

        t_333[k] = -hh_333[k]
                   + f_0 * kh_333[k];

        t_334[k] = -hh_334[k]
                   + f_0 * kh_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, hh_335, hh_336, hh_337, hh_338, \
                         hh_339, kh_335, kh_336, kh_337, kh_338, \
                         kh_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -hh_335[k]
                   + f_0 * kh_335[k];

        t_336[k] = -hh_336[k]
                   + f_0 * kh_336[k];

        t_337[k] = -hh_337[k]
                   + f_0 * kh_337[k];

        t_338[k] = -hh_338[k]
                   + f_0 * kh_338[k];

        t_339[k] = -hh_339[k]
                   + f_0 * kh_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, hh_340, hh_341, hh_342, hh_343, \
                         hh_344, kh_340, kh_341, kh_342, kh_343, \
                         kh_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -hh_340[k]
                   + f_0 * kh_340[k];

        t_341[k] = -hh_341[k]
                   + f_0 * kh_341[k];

        t_342[k] = -hh_342[k]
                   + f_0 * kh_342[k];

        t_343[k] = -hh_343[k]
                   + f_0 * kh_343[k];

        t_344[k] = -hh_344[k]
                   + f_0 * kh_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, hh_345, hh_346, hh_347, hh_348, \
                         hh_349, kh_345, kh_346, kh_347, kh_348, \
                         kh_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -hh_345[k]
                   + f_0 * kh_345[k];

        t_346[k] = -hh_346[k]
                   + f_0 * kh_346[k];

        t_347[k] = -hh_347[k]
                   + f_0 * kh_347[k];

        t_348[k] = -hh_348[k]
                   + f_0 * kh_348[k];

        t_349[k] = -hh_349[k]
                   + f_0 * kh_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, hh_350, hh_351, hh_352, hh_353, \
                         hh_354, kh_350, kh_351, kh_352, kh_353, \
                         kh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -hh_350[k]
                   + f_0 * kh_350[k];

        t_351[k] = -hh_351[k]
                   + f_0 * kh_351[k];

        t_352[k] = -hh_352[k]
                   + f_0 * kh_352[k];

        t_353[k] = -hh_353[k]
                   + f_0 * kh_353[k];

        t_354[k] = -hh_354[k]
                   + f_0 * kh_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, hh_355, hh_356, hh_357, hh_358, \
                         hh_359, kh_355, kh_356, kh_357, kh_358, \
                         kh_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -hh_355[k]
                   + f_0 * kh_355[k];

        t_356[k] = -hh_356[k]
                   + f_0 * kh_356[k];

        t_357[k] = -hh_357[k]
                   + f_0 * kh_357[k];

        t_358[k] = -hh_358[k]
                   + f_0 * kh_358[k];

        t_359[k] = -hh_359[k]
                   + f_0 * kh_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, hh_360, hh_361, hh_362, hh_363, \
                         hh_364, kh_360, kh_361, kh_362, kh_363, \
                         kh_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -hh_360[k]
                   + f_0 * kh_360[k];

        t_361[k] = -hh_361[k]
                   + f_0 * kh_361[k];

        t_362[k] = -hh_362[k]
                   + f_0 * kh_362[k];

        t_363[k] = -hh_363[k]
                   + f_0 * kh_363[k];

        t_364[k] = -hh_364[k]
                   + f_0 * kh_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, hh_365, hh_366, hh_367, hh_368, \
                         hh_369, kh_365, kh_366, kh_367, kh_368, \
                         kh_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -hh_365[k]
                   + f_0 * kh_365[k];

        t_366[k] = -hh_366[k]
                   + f_0 * kh_366[k];

        t_367[k] = -hh_367[k]
                   + f_0 * kh_367[k];

        t_368[k] = -hh_368[k]
                   + f_0 * kh_368[k];

        t_369[k] = -hh_369[k]
                   + f_0 * kh_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, hh_370, hh_371, hh_372, hh_373, \
                         hh_374, kh_370, kh_371, kh_372, kh_373, \
                         kh_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -hh_370[k]
                   + f_0 * kh_370[k];

        t_371[k] = -hh_371[k]
                   + f_0 * kh_371[k];

        t_372[k] = -hh_372[k]
                   + f_0 * kh_372[k];

        t_373[k] = -hh_373[k]
                   + f_0 * kh_373[k];

        t_374[k] = -hh_374[k]
                   + f_0 * kh_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, hh_375, hh_376, hh_377, hh_378, \
                         hh_379, kh_375, kh_376, kh_377, kh_378, \
                         kh_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -hh_375[k]
                   + f_0 * kh_375[k];

        t_376[k] = -hh_376[k]
                   + f_0 * kh_376[k];

        t_377[k] = -hh_377[k]
                   + f_0 * kh_377[k];

        t_378[k] = -hh_378[k]
                   + f_0 * kh_378[k];

        t_379[k] = -hh_379[k]
                   + f_0 * kh_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, hh_380, hh_381, hh_382, hh_383, \
                         hh_384, kh_380, kh_381, kh_382, kh_383, \
                         kh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -hh_380[k]
                   + f_0 * kh_380[k];

        t_381[k] = -hh_381[k]
                   + f_0 * kh_381[k];

        t_382[k] = -hh_382[k]
                   + f_0 * kh_382[k];

        t_383[k] = -hh_383[k]
                   + f_0 * kh_383[k];

        t_384[k] = -hh_384[k]
                   + f_0 * kh_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, hh_385, hh_386, hh_387, hh_388, \
                         hh_389, kh_385, kh_386, kh_387, kh_388, \
                         kh_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -hh_385[k]
                   + f_0 * kh_385[k];

        t_386[k] = -hh_386[k]
                   + f_0 * kh_386[k];

        t_387[k] = -hh_387[k]
                   + f_0 * kh_387[k];

        t_388[k] = -hh_388[k]
                   + f_0 * kh_388[k];

        t_389[k] = -hh_389[k]
                   + f_0 * kh_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, hh_390, hh_391, hh_392, hh_393, \
                         hh_394, kh_390, kh_391, kh_392, kh_393, \
                         kh_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -hh_390[k]
                   + f_0 * kh_390[k];

        t_391[k] = -hh_391[k]
                   + f_0 * kh_391[k];

        t_392[k] = -hh_392[k]
                   + f_0 * kh_392[k];

        t_393[k] = -hh_393[k]
                   + f_0 * kh_393[k];

        t_394[k] = -hh_394[k]
                   + f_0 * kh_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, hh_395, hh_396, hh_397, hh_398, \
                         hh_399, kh_395, kh_396, kh_397, kh_398, \
                         kh_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -hh_395[k]
                   + f_0 * kh_395[k];

        t_396[k] = -hh_396[k]
                   + f_0 * kh_396[k];

        t_397[k] = -hh_397[k]
                   + f_0 * kh_397[k];

        t_398[k] = -hh_398[k]
                   + f_0 * kh_398[k];

        t_399[k] = -hh_399[k]
                   + f_0 * kh_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, hh_400, hh_401, hh_402, hh_403, \
                         hh_404, kh_400, kh_401, kh_402, kh_403, \
                         kh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -hh_400[k]
                   + f_0 * kh_400[k];

        t_401[k] = -hh_401[k]
                   + f_0 * kh_401[k];

        t_402[k] = -hh_402[k]
                   + f_0 * kh_402[k];

        t_403[k] = -hh_403[k]
                   + f_0 * kh_403[k];

        t_404[k] = -hh_404[k]
                   + f_0 * kh_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, hh_405, hh_406, hh_407, hh_408, \
                         hh_409, kh_405, kh_406, kh_407, kh_408, \
                         kh_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -hh_405[k]
                   + f_0 * kh_405[k];

        t_406[k] = -hh_406[k]
                   + f_0 * kh_406[k];

        t_407[k] = -hh_407[k]
                   + f_0 * kh_407[k];

        t_408[k] = -hh_408[k]
                   + f_0 * kh_408[k];

        t_409[k] = -hh_409[k]
                   + f_0 * kh_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, hh_410, hh_411, hh_412, hh_413, \
                         hh_414, kh_410, kh_411, kh_412, kh_413, \
                         kh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -hh_410[k]
                   + f_0 * kh_410[k];

        t_411[k] = -hh_411[k]
                   + f_0 * kh_411[k];

        t_412[k] = -hh_412[k]
                   + f_0 * kh_412[k];

        t_413[k] = -hh_413[k]
                   + f_0 * kh_413[k];

        t_414[k] = -hh_414[k]
                   + f_0 * kh_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, hh_415, hh_416, hh_417, hh_418, \
                         hh_419, kh_415, kh_416, kh_417, kh_418, \
                         kh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -hh_415[k]
                   + f_0 * kh_415[k];

        t_416[k] = -hh_416[k]
                   + f_0 * kh_416[k];

        t_417[k] = -hh_417[k]
                   + f_0 * kh_417[k];

        t_418[k] = -hh_418[k]
                   + f_0 * kh_418[k];

        t_419[k] = -hh_419[k]
                   + f_0 * kh_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, hh_420, hh_421, hh_422, hh_423, \
                         hh_424, kh_420, kh_421, kh_422, kh_423, \
                         kh_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -hh_420[k]
                   + f_0 * kh_420[k];

        t_421[k] = -hh_421[k]
                   + f_0 * kh_421[k];

        t_422[k] = -hh_422[k]
                   + f_0 * kh_422[k];

        t_423[k] = -hh_423[k]
                   + f_0 * kh_423[k];

        t_424[k] = -hh_424[k]
                   + f_0 * kh_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, hh_425, hh_426, hh_427, hh_428, \
                         hh_429, kh_425, kh_426, kh_427, kh_428, \
                         kh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -hh_425[k]
                   + f_0 * kh_425[k];

        t_426[k] = -hh_426[k]
                   + f_0 * kh_426[k];

        t_427[k] = -hh_427[k]
                   + f_0 * kh_427[k];

        t_428[k] = -hh_428[k]
                   + f_0 * kh_428[k];

        t_429[k] = -hh_429[k]
                   + f_0 * kh_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, hh_430, hh_431, hh_432, hh_433, \
                         hh_434, kh_430, kh_431, kh_432, kh_433, \
                         kh_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -hh_430[k]
                   + f_0 * kh_430[k];

        t_431[k] = -hh_431[k]
                   + f_0 * kh_431[k];

        t_432[k] = -hh_432[k]
                   + f_0 * kh_432[k];

        t_433[k] = -hh_433[k]
                   + f_0 * kh_433[k];

        t_434[k] = -hh_434[k]
                   + f_0 * kh_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, hh_435, hh_436, hh_437, hh_438, \
                         hh_439, kh_435, kh_436, kh_437, kh_438, \
                         kh_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -hh_435[k]
                   + f_0 * kh_435[k];

        t_436[k] = -hh_436[k]
                   + f_0 * kh_436[k];

        t_437[k] = -hh_437[k]
                   + f_0 * kh_437[k];

        t_438[k] = -hh_438[k]
                   + f_0 * kh_438[k];

        t_439[k] = -hh_439[k]
                   + f_0 * kh_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, t_445, t_446, hh_440, kh_440, \
                         kh_441, kh_442, kh_443, kh_444, kh_445, \
                         kh_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -hh_440[k]
                   + f_0 * kh_440[k];

        t_441[k] = f_0 * kh_441[k];

        t_442[k] = f_0 * kh_442[k];

        t_443[k] = f_0 * kh_443[k];

        t_444[k] = f_0 * kh_444[k];

        t_445[k] = f_0 * kh_445[k];

        t_446[k] = f_0 * kh_446[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, t_451, t_452, t_453, t_454, kh_447, \
                         kh_448, kh_449, kh_450, kh_451, kh_452, kh_453, \
                         kh_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_0 * kh_447[k];

        t_448[k] = f_0 * kh_448[k];

        t_449[k] = f_0 * kh_449[k];

        t_450[k] = f_0 * kh_450[k];

        t_451[k] = f_0 * kh_451[k];

        t_452[k] = f_0 * kh_452[k];

        t_453[k] = f_0 * kh_453[k];

        t_454[k] = f_0 * kh_454[k];
    }
}

static auto
compute_prim_geom_10_ih_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kh, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, t_460, t_461, t_462, kh_455, \
                         kh_456, kh_457, kh_458, kh_459, kh_460, kh_461, \
                         kh_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_0 * kh_455[k];

        t_456[k] = f_0 * kh_456[k];

        t_457[k] = f_0 * kh_457[k];

        t_458[k] = f_0 * kh_458[k];

        t_459[k] = f_0 * kh_459[k];

        t_460[k] = f_0 * kh_460[k];

        t_461[k] = f_0 * kh_461[k];

        t_462[k] = f_0 * kh_462[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, t_467, t_468, t_469, t_470, kh_463, \
                         kh_464, kh_465, kh_466, kh_467, kh_468, kh_469, \
                         kh_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_0 * kh_463[k];

        t_464[k] = f_0 * kh_464[k];

        t_465[k] = f_0 * kh_465[k];

        t_466[k] = f_0 * kh_466[k];

        t_467[k] = f_0 * kh_467[k];

        t_468[k] = f_0 * kh_468[k];

        t_469[k] = f_0 * kh_469[k];

        t_470[k] = f_0 * kh_470[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, t_475, t_476, t_477, t_478, kh_471, \
                         kh_472, kh_473, kh_474, kh_475, kh_476, kh_477, \
                         kh_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_0 * kh_471[k];

        t_472[k] = f_0 * kh_472[k];

        t_473[k] = f_0 * kh_473[k];

        t_474[k] = f_0 * kh_474[k];

        t_475[k] = f_0 * kh_475[k];

        t_476[k] = f_0 * kh_476[k];

        t_477[k] = f_0 * kh_477[k];

        t_478[k] = f_0 * kh_478[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, t_483, t_484, t_485, t_486, kh_479, \
                         kh_480, kh_481, kh_482, kh_483, kh_484, kh_485, \
                         kh_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_0 * kh_479[k];

        t_480[k] = f_0 * kh_480[k];

        t_481[k] = f_0 * kh_481[k];

        t_482[k] = f_0 * kh_482[k];

        t_483[k] = f_0 * kh_483[k];

        t_484[k] = f_0 * kh_484[k];

        t_485[k] = f_0 * kh_485[k];

        t_486[k] = f_0 * kh_486[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, t_492, t_493, t_494, kh_487, \
                         kh_488, kh_489, kh_490, kh_491, kh_492, kh_493, \
                         kh_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_0 * kh_487[k];

        t_488[k] = f_0 * kh_488[k];

        t_489[k] = f_0 * kh_489[k];

        t_490[k] = f_0 * kh_490[k];

        t_491[k] = f_0 * kh_491[k];

        t_492[k] = f_0 * kh_492[k];

        t_493[k] = f_0 * kh_493[k];

        t_494[k] = f_0 * kh_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, t_500, t_501, t_502, kh_495, \
                         kh_496, kh_497, kh_498, kh_499, kh_500, kh_501, \
                         kh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_0 * kh_495[k];

        t_496[k] = f_0 * kh_496[k];

        t_497[k] = f_0 * kh_497[k];

        t_498[k] = f_0 * kh_498[k];

        t_499[k] = f_0 * kh_499[k];

        t_500[k] = f_0 * kh_500[k];

        t_501[k] = f_0 * kh_501[k];

        t_502[k] = f_0 * kh_502[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, t_507, t_508, t_509, t_510, kh_503, \
                         kh_504, kh_505, kh_506, kh_507, kh_508, kh_509, \
                         kh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_0 * kh_503[k];

        t_504[k] = f_0 * kh_504[k];

        t_505[k] = f_0 * kh_505[k];

        t_506[k] = f_0 * kh_506[k];

        t_507[k] = f_0 * kh_507[k];

        t_508[k] = f_0 * kh_508[k];

        t_509[k] = f_0 * kh_509[k];

        t_510[k] = f_0 * kh_510[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, t_515, t_516, t_517, t_518, kh_511, \
                         kh_512, kh_513, kh_514, kh_515, kh_516, kh_517, \
                         kh_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_0 * kh_511[k];

        t_512[k] = f_0 * kh_512[k];

        t_513[k] = f_0 * kh_513[k];

        t_514[k] = f_0 * kh_514[k];

        t_515[k] = f_0 * kh_515[k];

        t_516[k] = f_0 * kh_516[k];

        t_517[k] = f_0 * kh_517[k];

        t_518[k] = f_0 * kh_518[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, t_523, t_524, t_525, t_526, kh_519, \
                         kh_520, kh_521, kh_522, kh_523, kh_524, kh_525, \
                         kh_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_0 * kh_519[k];

        t_520[k] = f_0 * kh_520[k];

        t_521[k] = f_0 * kh_521[k];

        t_522[k] = f_0 * kh_522[k];

        t_523[k] = f_0 * kh_523[k];

        t_524[k] = f_0 * kh_524[k];

        t_525[k] = f_0 * kh_525[k];

        t_526[k] = f_0 * kh_526[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, t_532, t_533, t_534, kh_527, \
                         kh_528, kh_529, kh_530, kh_531, kh_532, kh_533, \
                         kh_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_0 * kh_527[k];

        t_528[k] = f_0 * kh_528[k];

        t_529[k] = f_0 * kh_529[k];

        t_530[k] = f_0 * kh_530[k];

        t_531[k] = f_0 * kh_531[k];

        t_532[k] = f_0 * kh_532[k];

        t_533[k] = f_0 * kh_533[k];

        t_534[k] = f_0 * kh_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, t_540, t_541, t_542, kh_535, \
                         kh_536, kh_537, kh_538, kh_539, kh_540, kh_541, \
                         kh_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = f_0 * kh_535[k];

        t_536[k] = f_0 * kh_536[k];

        t_537[k] = f_0 * kh_537[k];

        t_538[k] = f_0 * kh_538[k];

        t_539[k] = f_0 * kh_539[k];

        t_540[k] = f_0 * kh_540[k];

        t_541[k] = f_0 * kh_541[k];

        t_542[k] = f_0 * kh_542[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, t_547, t_548, t_549, t_550, kh_543, \
                         kh_544, kh_545, kh_546, kh_547, kh_548, kh_549, \
                         kh_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_0 * kh_543[k];

        t_544[k] = f_0 * kh_544[k];

        t_545[k] = f_0 * kh_545[k];

        t_546[k] = f_0 * kh_546[k];

        t_547[k] = f_0 * kh_547[k];

        t_548[k] = f_0 * kh_548[k];

        t_549[k] = f_0 * kh_549[k];

        t_550[k] = f_0 * kh_550[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, t_554, t_555, t_556, t_557, t_558, kh_551, \
                         kh_552, kh_553, kh_554, kh_555, kh_556, kh_557, \
                         kh_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_0 * kh_551[k];

        t_552[k] = f_0 * kh_552[k];

        t_553[k] = f_0 * kh_553[k];

        t_554[k] = f_0 * kh_554[k];

        t_555[k] = f_0 * kh_555[k];

        t_556[k] = f_0 * kh_556[k];

        t_557[k] = f_0 * kh_557[k];

        t_558[k] = f_0 * kh_558[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, t_563, t_564, t_565, t_566, kh_559, \
                         kh_560, kh_561, kh_562, kh_563, kh_564, kh_565, \
                         kh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_0 * kh_559[k];

        t_560[k] = f_0 * kh_560[k];

        t_561[k] = f_0 * kh_561[k];

        t_562[k] = f_0 * kh_562[k];

        t_563[k] = f_0 * kh_563[k];

        t_564[k] = f_0 * kh_564[k];

        t_565[k] = f_0 * kh_565[k];

        t_566[k] = f_0 * kh_566[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, t_572, t_573, t_574, kh_567, \
                         kh_568, kh_569, kh_570, kh_571, kh_572, kh_573, \
                         kh_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_0 * kh_567[k];

        t_568[k] = f_0 * kh_568[k];

        t_569[k] = f_0 * kh_569[k];

        t_570[k] = f_0 * kh_570[k];

        t_571[k] = f_0 * kh_571[k];

        t_572[k] = f_0 * kh_572[k];

        t_573[k] = f_0 * kh_573[k];

        t_574[k] = f_0 * kh_574[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, t_580, t_581, t_582, kh_575, \
                         kh_576, kh_577, kh_578, kh_579, kh_580, kh_581, \
                         kh_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_0 * kh_575[k];

        t_576[k] = f_0 * kh_576[k];

        t_577[k] = f_0 * kh_577[k];

        t_578[k] = f_0 * kh_578[k];

        t_579[k] = f_0 * kh_579[k];

        t_580[k] = f_0 * kh_580[k];

        t_581[k] = f_0 * kh_581[k];

        t_582[k] = f_0 * kh_582[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, t_586, t_587, kh_583, kh_584, kh_585, kh_586, \
                         kh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_0 * kh_583[k];

        t_584[k] = f_0 * kh_584[k];

        t_585[k] = f_0 * kh_585[k];

        t_586[k] = f_0 * kh_586[k];

        t_587[k] = f_0 * kh_587[k];
    }
}

auto
compute_prim_geom_10_ih_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t hh, const size_t kh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ih_electron_repulsion_0_piece0(buffer, target, hh, kh, ncols, alpha);

    compute_prim_geom_10_ih_electron_repulsion_0_piece1(buffer, target, hh, kh, ncols, alpha);

    compute_prim_geom_10_ih_electron_repulsion_0_piece2(buffer, target, hh, kh, ncols, alpha);

    compute_prim_geom_10_ih_electron_repulsion_0_piece3(buffer, target, kh, ncols, alpha);
}

static auto
compute_prim_geom_10_ih_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hh, const size_t kh,
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

    const auto *hh_0 = buffer.data(hh + 0);
    const auto *hh_1 = buffer.data(hh + 1);
    const auto *hh_2 = buffer.data(hh + 2);
    const auto *hh_3 = buffer.data(hh + 3);
    const auto *hh_4 = buffer.data(hh + 4);
    const auto *hh_5 = buffer.data(hh + 5);
    const auto *hh_6 = buffer.data(hh + 6);
    const auto *hh_7 = buffer.data(hh + 7);
    const auto *hh_8 = buffer.data(hh + 8);
    const auto *hh_9 = buffer.data(hh + 9);
    const auto *hh_10 = buffer.data(hh + 10);
    const auto *hh_11 = buffer.data(hh + 11);
    const auto *hh_12 = buffer.data(hh + 12);
    const auto *hh_13 = buffer.data(hh + 13);
    const auto *hh_14 = buffer.data(hh + 14);
    const auto *hh_15 = buffer.data(hh + 15);
    const auto *hh_16 = buffer.data(hh + 16);
    const auto *hh_17 = buffer.data(hh + 17);
    const auto *hh_18 = buffer.data(hh + 18);
    const auto *hh_19 = buffer.data(hh + 19);
    const auto *hh_20 = buffer.data(hh + 20);
    const auto *hh_21 = buffer.data(hh + 21);
    const auto *hh_22 = buffer.data(hh + 22);
    const auto *hh_23 = buffer.data(hh + 23);
    const auto *hh_24 = buffer.data(hh + 24);
    const auto *hh_25 = buffer.data(hh + 25);
    const auto *hh_26 = buffer.data(hh + 26);
    const auto *hh_27 = buffer.data(hh + 27);
    const auto *hh_28 = buffer.data(hh + 28);
    const auto *hh_29 = buffer.data(hh + 29);
    const auto *hh_30 = buffer.data(hh + 30);
    const auto *hh_31 = buffer.data(hh + 31);
    const auto *hh_32 = buffer.data(hh + 32);
    const auto *hh_33 = buffer.data(hh + 33);
    const auto *hh_34 = buffer.data(hh + 34);
    const auto *hh_35 = buffer.data(hh + 35);
    const auto *hh_36 = buffer.data(hh + 36);
    const auto *hh_37 = buffer.data(hh + 37);
    const auto *hh_38 = buffer.data(hh + 38);
    const auto *hh_39 = buffer.data(hh + 39);
    const auto *hh_40 = buffer.data(hh + 40);
    const auto *hh_41 = buffer.data(hh + 41);
    const auto *hh_42 = buffer.data(hh + 42);
    const auto *hh_43 = buffer.data(hh + 43);
    const auto *hh_44 = buffer.data(hh + 44);
    const auto *hh_45 = buffer.data(hh + 45);
    const auto *hh_46 = buffer.data(hh + 46);
    const auto *hh_47 = buffer.data(hh + 47);
    const auto *hh_48 = buffer.data(hh + 48);
    const auto *hh_49 = buffer.data(hh + 49);
    const auto *hh_50 = buffer.data(hh + 50);
    const auto *hh_51 = buffer.data(hh + 51);
    const auto *hh_52 = buffer.data(hh + 52);
    const auto *hh_53 = buffer.data(hh + 53);
    const auto *hh_54 = buffer.data(hh + 54);
    const auto *hh_55 = buffer.data(hh + 55);
    const auto *hh_56 = buffer.data(hh + 56);
    const auto *hh_57 = buffer.data(hh + 57);
    const auto *hh_58 = buffer.data(hh + 58);
    const auto *hh_59 = buffer.data(hh + 59);
    const auto *hh_60 = buffer.data(hh + 60);
    const auto *hh_61 = buffer.data(hh + 61);
    const auto *hh_62 = buffer.data(hh + 62);
    const auto *hh_63 = buffer.data(hh + 63);
    const auto *hh_64 = buffer.data(hh + 64);
    const auto *hh_65 = buffer.data(hh + 65);
    const auto *hh_66 = buffer.data(hh + 66);
    const auto *hh_67 = buffer.data(hh + 67);
    const auto *hh_68 = buffer.data(hh + 68);
    const auto *hh_69 = buffer.data(hh + 69);
    const auto *hh_70 = buffer.data(hh + 70);
    const auto *hh_71 = buffer.data(hh + 71);
    const auto *hh_72 = buffer.data(hh + 72);
    const auto *hh_73 = buffer.data(hh + 73);
    const auto *hh_74 = buffer.data(hh + 74);
    const auto *hh_75 = buffer.data(hh + 75);
    const auto *hh_76 = buffer.data(hh + 76);
    const auto *hh_77 = buffer.data(hh + 77);
    const auto *hh_78 = buffer.data(hh + 78);
    const auto *hh_79 = buffer.data(hh + 79);
    const auto *hh_80 = buffer.data(hh + 80);
    const auto *hh_81 = buffer.data(hh + 81);
    const auto *hh_82 = buffer.data(hh + 82);
    const auto *hh_83 = buffer.data(hh + 83);
    const auto *hh_84 = buffer.data(hh + 84);
    const auto *hh_85 = buffer.data(hh + 85);
    const auto *hh_86 = buffer.data(hh + 86);
    const auto *hh_87 = buffer.data(hh + 87);
    const auto *hh_88 = buffer.data(hh + 88);
    const auto *hh_89 = buffer.data(hh + 89);
    const auto *hh_90 = buffer.data(hh + 90);
    const auto *hh_91 = buffer.data(hh + 91);
    const auto *hh_92 = buffer.data(hh + 92);
    const auto *hh_93 = buffer.data(hh + 93);
    const auto *hh_94 = buffer.data(hh + 94);
    const auto *hh_95 = buffer.data(hh + 95);
    const auto *hh_96 = buffer.data(hh + 96);
    const auto *hh_97 = buffer.data(hh + 97);
    const auto *hh_98 = buffer.data(hh + 98);
    const auto *hh_99 = buffer.data(hh + 99);
    const auto *hh_100 = buffer.data(hh + 100);
    const auto *hh_101 = buffer.data(hh + 101);
    const auto *hh_102 = buffer.data(hh + 102);
    const auto *hh_103 = buffer.data(hh + 103);
    const auto *hh_104 = buffer.data(hh + 104);
    const auto *hh_105 = buffer.data(hh + 105);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, kh_21, kh_22, kh_23, kh_24, \
                         kh_25, kh_26, kh_27, kh_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kh_21[k];

        t_1[k] = f_0 * kh_22[k];

        t_2[k] = f_0 * kh_23[k];

        t_3[k] = f_0 * kh_24[k];

        t_4[k] = f_0 * kh_25[k];

        t_5[k] = f_0 * kh_26[k];

        t_6[k] = f_0 * kh_27[k];

        t_7[k] = f_0 * kh_28[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, kh_29, kh_30, kh_31, \
                         kh_32, kh_33, kh_34, kh_35, kh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * kh_29[k];

        t_9[k] = f_0 * kh_30[k];

        t_10[k] = f_0 * kh_31[k];

        t_11[k] = f_0 * kh_32[k];

        t_12[k] = f_0 * kh_33[k];

        t_13[k] = f_0 * kh_34[k];

        t_14[k] = f_0 * kh_35[k];

        t_15[k] = f_0 * kh_36[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, hh_0, hh_1, kh_37, kh_38, \
                         kh_39, kh_40, kh_41, kh_63, kh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * kh_37[k];

        t_17[k] = f_0 * kh_38[k];

        t_18[k] = f_0 * kh_39[k];

        t_19[k] = f_0 * kh_40[k];

        t_20[k] = f_0 * kh_41[k];

        t_21[k] = -hh_0[k]
                  + f_0 * kh_63[k];

        t_22[k] = -hh_1[k]
                  + f_0 * kh_64[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, hh_2, hh_3, hh_4, hh_5, hh_6, kh_65, \
                         kh_66, kh_67, kh_68, kh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -hh_2[k]
                  + f_0 * kh_65[k];

        t_24[k] = -hh_3[k]
                  + f_0 * kh_66[k];

        t_25[k] = -hh_4[k]
                  + f_0 * kh_67[k];

        t_26[k] = -hh_5[k]
                  + f_0 * kh_68[k];

        t_27[k] = -hh_6[k]
                  + f_0 * kh_69[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, hh_7, hh_8, hh_9, hh_10, hh_11, kh_70, \
                         kh_71, kh_72, kh_73, kh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -hh_7[k]
                  + f_0 * kh_70[k];

        t_29[k] = -hh_8[k]
                  + f_0 * kh_71[k];

        t_30[k] = -hh_9[k]
                  + f_0 * kh_72[k];

        t_31[k] = -hh_10[k]
                  + f_0 * kh_73[k];

        t_32[k] = -hh_11[k]
                  + f_0 * kh_74[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, hh_12, hh_13, hh_14, hh_15, hh_16, \
                         kh_75, kh_76, kh_77, kh_78, kh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -hh_12[k]
                  + f_0 * kh_75[k];

        t_34[k] = -hh_13[k]
                  + f_0 * kh_76[k];

        t_35[k] = -hh_14[k]
                  + f_0 * kh_77[k];

        t_36[k] = -hh_15[k]
                  + f_0 * kh_78[k];

        t_37[k] = -hh_16[k]
                  + f_0 * kh_79[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, hh_17, hh_18, hh_19, hh_20, \
                         kh_80, kh_81, kh_82, kh_83, kh_84, kh_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -hh_17[k]
                  + f_0 * kh_80[k];

        t_39[k] = -hh_18[k]
                  + f_0 * kh_81[k];

        t_40[k] = -hh_19[k]
                  + f_0 * kh_82[k];

        t_41[k] = -hh_20[k]
                  + f_0 * kh_83[k];

        t_42[k] = f_0 * kh_84[k];

        t_43[k] = f_0 * kh_85[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, t_51, kh_86, kh_87, kh_88, \
                         kh_89, kh_90, kh_91, kh_92, kh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * kh_86[k];

        t_45[k] = f_0 * kh_87[k];

        t_46[k] = f_0 * kh_88[k];

        t_47[k] = f_0 * kh_89[k];

        t_48[k] = f_0 * kh_90[k];

        t_49[k] = f_0 * kh_91[k];

        t_50[k] = f_0 * kh_92[k];

        t_51[k] = f_0 * kh_93[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, t_58, t_59, kh_94, kh_95, kh_96, \
                         kh_97, kh_98, kh_99, kh_100, kh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * kh_94[k];

        t_53[k] = f_0 * kh_95[k];

        t_54[k] = f_0 * kh_96[k];

        t_55[k] = f_0 * kh_97[k];

        t_56[k] = f_0 * kh_98[k];

        t_57[k] = f_0 * kh_99[k];

        t_58[k] = f_0 * kh_100[k];

        t_59[k] = f_0 * kh_101[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, hh_21, hh_22, hh_23, kh_102, \
                         kh_103, kh_104, kh_126, kh_127, kh_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * kh_102[k];

        t_61[k] = f_0 * kh_103[k];

        t_62[k] = f_0 * kh_104[k];

        t_63[k] = -2.0 * hh_21[k]
                  + f_0 * kh_126[k];

        t_64[k] = -2.0 * hh_22[k]
                  + f_0 * kh_127[k];

        t_65[k] = -2.0 * hh_23[k]
                  + f_0 * kh_128[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, hh_24, hh_25, hh_26, hh_27, hh_28, \
                         kh_129, kh_130, kh_131, kh_132, kh_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -2.0 * hh_24[k]
                  + f_0 * kh_129[k];

        t_67[k] = -2.0 * hh_25[k]
                  + f_0 * kh_130[k];

        t_68[k] = -2.0 * hh_26[k]
                  + f_0 * kh_131[k];

        t_69[k] = -2.0 * hh_27[k]
                  + f_0 * kh_132[k];

        t_70[k] = -2.0 * hh_28[k]
                  + f_0 * kh_133[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, hh_29, hh_30, hh_31, hh_32, hh_33, \
                         kh_134, kh_135, kh_136, kh_137, kh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -2.0 * hh_29[k]
                  + f_0 * kh_134[k];

        t_72[k] = -2.0 * hh_30[k]
                  + f_0 * kh_135[k];

        t_73[k] = -2.0 * hh_31[k]
                  + f_0 * kh_136[k];

        t_74[k] = -2.0 * hh_32[k]
                  + f_0 * kh_137[k];

        t_75[k] = -2.0 * hh_33[k]
                  + f_0 * kh_138[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, hh_34, hh_35, hh_36, hh_37, hh_38, \
                         kh_139, kh_140, kh_141, kh_142, kh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -2.0 * hh_34[k]
                  + f_0 * kh_139[k];

        t_77[k] = -2.0 * hh_35[k]
                  + f_0 * kh_140[k];

        t_78[k] = -2.0 * hh_36[k]
                  + f_0 * kh_141[k];

        t_79[k] = -2.0 * hh_37[k]
                  + f_0 * kh_142[k];

        t_80[k] = -2.0 * hh_38[k]
                  + f_0 * kh_143[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, hh_39, hh_40, hh_41, hh_42, hh_43, \
                         kh_144, kh_145, kh_146, kh_147, kh_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -2.0 * hh_39[k]
                  + f_0 * kh_144[k];

        t_82[k] = -2.0 * hh_40[k]
                  + f_0 * kh_145[k];

        t_83[k] = -2.0 * hh_41[k]
                  + f_0 * kh_146[k];

        t_84[k] = -hh_42[k]
                  + f_0 * kh_147[k];

        t_85[k] = -hh_43[k]
                  + f_0 * kh_148[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, hh_44, hh_45, hh_46, hh_47, hh_48, \
                         kh_149, kh_150, kh_151, kh_152, kh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -hh_44[k]
                  + f_0 * kh_149[k];

        t_87[k] = -hh_45[k]
                  + f_0 * kh_150[k];

        t_88[k] = -hh_46[k]
                  + f_0 * kh_151[k];

        t_89[k] = -hh_47[k]
                  + f_0 * kh_152[k];

        t_90[k] = -hh_48[k]
                  + f_0 * kh_153[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, hh_49, hh_50, hh_51, hh_52, hh_53, \
                         kh_154, kh_155, kh_156, kh_157, kh_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = -hh_49[k]
                  + f_0 * kh_154[k];

        t_92[k] = -hh_50[k]
                  + f_0 * kh_155[k];

        t_93[k] = -hh_51[k]
                  + f_0 * kh_156[k];

        t_94[k] = -hh_52[k]
                  + f_0 * kh_157[k];

        t_95[k] = -hh_53[k]
                  + f_0 * kh_158[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, hh_54, hh_55, hh_56, hh_57, hh_58, \
                         kh_159, kh_160, kh_161, kh_162, kh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = -hh_54[k]
                  + f_0 * kh_159[k];

        t_97[k] = -hh_55[k]
                  + f_0 * kh_160[k];

        t_98[k] = -hh_56[k]
                  + f_0 * kh_161[k];

        t_99[k] = -hh_57[k]
                  + f_0 * kh_162[k];

        t_100[k] = -hh_58[k]
                   + f_0 * kh_163[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, hh_59, hh_60, hh_61, hh_62, \
                         kh_164, kh_165, kh_166, kh_167, kh_168, \
                         kh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = -hh_59[k]
                   + f_0 * kh_164[k];

        t_102[k] = -hh_60[k]
                   + f_0 * kh_165[k];

        t_103[k] = -hh_61[k]
                   + f_0 * kh_166[k];

        t_104[k] = -hh_62[k]
                   + f_0 * kh_167[k];

        t_105[k] = f_0 * kh_168[k];

        t_106[k] = f_0 * kh_169[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, t_114, kh_170, \
                         kh_171, kh_172, kh_173, kh_174, kh_175, kh_176, \
                         kh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * kh_170[k];

        t_108[k] = f_0 * kh_171[k];

        t_109[k] = f_0 * kh_172[k];

        t_110[k] = f_0 * kh_173[k];

        t_111[k] = f_0 * kh_174[k];

        t_112[k] = f_0 * kh_175[k];

        t_113[k] = f_0 * kh_176[k];

        t_114[k] = f_0 * kh_177[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, t_121, t_122, kh_178, \
                         kh_179, kh_180, kh_181, kh_182, kh_183, kh_184, \
                         kh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_0 * kh_178[k];

        t_116[k] = f_0 * kh_179[k];

        t_117[k] = f_0 * kh_180[k];

        t_118[k] = f_0 * kh_181[k];

        t_119[k] = f_0 * kh_182[k];

        t_120[k] = f_0 * kh_183[k];

        t_121[k] = f_0 * kh_184[k];

        t_122[k] = f_0 * kh_185[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, t_127, t_128, hh_63, hh_64, hh_65, \
                         kh_186, kh_187, kh_188, kh_210, kh_211, \
                         kh_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_0 * kh_186[k];

        t_124[k] = f_0 * kh_187[k];

        t_125[k] = f_0 * kh_188[k];

        t_126[k] = -3.0 * hh_63[k]
                   + f_0 * kh_210[k];

        t_127[k] = -3.0 * hh_64[k]
                   + f_0 * kh_211[k];

        t_128[k] = -3.0 * hh_65[k]
                   + f_0 * kh_212[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, hh_66, hh_67, hh_68, hh_69, hh_70, \
                         kh_213, kh_214, kh_215, kh_216, kh_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -3.0 * hh_66[k]
                   + f_0 * kh_213[k];

        t_130[k] = -3.0 * hh_67[k]
                   + f_0 * kh_214[k];

        t_131[k] = -3.0 * hh_68[k]
                   + f_0 * kh_215[k];

        t_132[k] = -3.0 * hh_69[k]
                   + f_0 * kh_216[k];

        t_133[k] = -3.0 * hh_70[k]
                   + f_0 * kh_217[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, hh_71, hh_72, hh_73, hh_74, hh_75, \
                         kh_218, kh_219, kh_220, kh_221, kh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -3.0 * hh_71[k]
                   + f_0 * kh_218[k];

        t_135[k] = -3.0 * hh_72[k]
                   + f_0 * kh_219[k];

        t_136[k] = -3.0 * hh_73[k]
                   + f_0 * kh_220[k];

        t_137[k] = -3.0 * hh_74[k]
                   + f_0 * kh_221[k];

        t_138[k] = -3.0 * hh_75[k]
                   + f_0 * kh_222[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, hh_76, hh_77, hh_78, hh_79, hh_80, \
                         kh_223, kh_224, kh_225, kh_226, kh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -3.0 * hh_76[k]
                   + f_0 * kh_223[k];

        t_140[k] = -3.0 * hh_77[k]
                   + f_0 * kh_224[k];

        t_141[k] = -3.0 * hh_78[k]
                   + f_0 * kh_225[k];

        t_142[k] = -3.0 * hh_79[k]
                   + f_0 * kh_226[k];

        t_143[k] = -3.0 * hh_80[k]
                   + f_0 * kh_227[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, hh_81, hh_82, hh_83, hh_84, hh_85, \
                         kh_228, kh_229, kh_230, kh_231, kh_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -3.0 * hh_81[k]
                   + f_0 * kh_228[k];

        t_145[k] = -3.0 * hh_82[k]
                   + f_0 * kh_229[k];

        t_146[k] = -3.0 * hh_83[k]
                   + f_0 * kh_230[k];

        t_147[k] = -2.0 * hh_84[k]
                   + f_0 * kh_231[k];

        t_148[k] = -2.0 * hh_85[k]
                   + f_0 * kh_232[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, hh_86, hh_87, hh_88, hh_89, hh_90, \
                         kh_233, kh_234, kh_235, kh_236, kh_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -2.0 * hh_86[k]
                   + f_0 * kh_233[k];

        t_150[k] = -2.0 * hh_87[k]
                   + f_0 * kh_234[k];

        t_151[k] = -2.0 * hh_88[k]
                   + f_0 * kh_235[k];

        t_152[k] = -2.0 * hh_89[k]
                   + f_0 * kh_236[k];

        t_153[k] = -2.0 * hh_90[k]
                   + f_0 * kh_237[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, hh_91, hh_92, hh_93, hh_94, hh_95, \
                         kh_238, kh_239, kh_240, kh_241, kh_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -2.0 * hh_91[k]
                   + f_0 * kh_238[k];

        t_155[k] = -2.0 * hh_92[k]
                   + f_0 * kh_239[k];

        t_156[k] = -2.0 * hh_93[k]
                   + f_0 * kh_240[k];

        t_157[k] = -2.0 * hh_94[k]
                   + f_0 * kh_241[k];

        t_158[k] = -2.0 * hh_95[k]
                   + f_0 * kh_242[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, hh_96, hh_97, hh_98, hh_99, \
                         hh_100, kh_243, kh_244, kh_245, kh_246, \
                         kh_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -2.0 * hh_96[k]
                   + f_0 * kh_243[k];

        t_160[k] = -2.0 * hh_97[k]
                   + f_0 * kh_244[k];

        t_161[k] = -2.0 * hh_98[k]
                   + f_0 * kh_245[k];

        t_162[k] = -2.0 * hh_99[k]
                   + f_0 * kh_246[k];

        t_163[k] = -2.0 * hh_100[k]
                   + f_0 * kh_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, hh_101, hh_102, hh_103, hh_104, \
                         hh_105, kh_248, kh_249, kh_250, kh_251, \
                         kh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -2.0 * hh_101[k]
                   + f_0 * kh_248[k];

        t_165[k] = -2.0 * hh_102[k]
                   + f_0 * kh_249[k];

        t_166[k] = -2.0 * hh_103[k]
                   + f_0 * kh_250[k];

        t_167[k] = -2.0 * hh_104[k]
                   + f_0 * kh_251[k];

        t_168[k] = -hh_105[k]
                   + f_0 * kh_252[k];
    }
}

static auto
compute_prim_geom_10_ih_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hh, const size_t kh,
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

    const auto *hh_106 = buffer.data(hh + 106);
    const auto *hh_107 = buffer.data(hh + 107);
    const auto *hh_108 = buffer.data(hh + 108);
    const auto *hh_109 = buffer.data(hh + 109);
    const auto *hh_110 = buffer.data(hh + 110);
    const auto *hh_111 = buffer.data(hh + 111);
    const auto *hh_112 = buffer.data(hh + 112);
    const auto *hh_113 = buffer.data(hh + 113);
    const auto *hh_114 = buffer.data(hh + 114);
    const auto *hh_115 = buffer.data(hh + 115);
    const auto *hh_116 = buffer.data(hh + 116);
    const auto *hh_117 = buffer.data(hh + 117);
    const auto *hh_118 = buffer.data(hh + 118);
    const auto *hh_119 = buffer.data(hh + 119);
    const auto *hh_120 = buffer.data(hh + 120);
    const auto *hh_121 = buffer.data(hh + 121);
    const auto *hh_122 = buffer.data(hh + 122);
    const auto *hh_123 = buffer.data(hh + 123);
    const auto *hh_124 = buffer.data(hh + 124);
    const auto *hh_125 = buffer.data(hh + 125);
    const auto *hh_126 = buffer.data(hh + 126);
    const auto *hh_127 = buffer.data(hh + 127);
    const auto *hh_128 = buffer.data(hh + 128);
    const auto *hh_129 = buffer.data(hh + 129);
    const auto *hh_130 = buffer.data(hh + 130);
    const auto *hh_131 = buffer.data(hh + 131);
    const auto *hh_132 = buffer.data(hh + 132);
    const auto *hh_133 = buffer.data(hh + 133);
    const auto *hh_134 = buffer.data(hh + 134);
    const auto *hh_135 = buffer.data(hh + 135);
    const auto *hh_136 = buffer.data(hh + 136);
    const auto *hh_137 = buffer.data(hh + 137);
    const auto *hh_138 = buffer.data(hh + 138);
    const auto *hh_139 = buffer.data(hh + 139);
    const auto *hh_140 = buffer.data(hh + 140);
    const auto *hh_141 = buffer.data(hh + 141);
    const auto *hh_142 = buffer.data(hh + 142);
    const auto *hh_143 = buffer.data(hh + 143);
    const auto *hh_144 = buffer.data(hh + 144);
    const auto *hh_145 = buffer.data(hh + 145);
    const auto *hh_146 = buffer.data(hh + 146);
    const auto *hh_147 = buffer.data(hh + 147);
    const auto *hh_148 = buffer.data(hh + 148);
    const auto *hh_149 = buffer.data(hh + 149);
    const auto *hh_150 = buffer.data(hh + 150);
    const auto *hh_151 = buffer.data(hh + 151);
    const auto *hh_152 = buffer.data(hh + 152);
    const auto *hh_153 = buffer.data(hh + 153);
    const auto *hh_154 = buffer.data(hh + 154);
    const auto *hh_155 = buffer.data(hh + 155);
    const auto *hh_156 = buffer.data(hh + 156);
    const auto *hh_157 = buffer.data(hh + 157);
    const auto *hh_158 = buffer.data(hh + 158);
    const auto *hh_159 = buffer.data(hh + 159);
    const auto *hh_160 = buffer.data(hh + 160);
    const auto *hh_161 = buffer.data(hh + 161);
    const auto *hh_162 = buffer.data(hh + 162);
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_164 = buffer.data(hh + 164);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_166 = buffer.data(hh + 166);
    const auto *hh_167 = buffer.data(hh + 167);
    const auto *hh_168 = buffer.data(hh + 168);
    const auto *hh_169 = buffer.data(hh + 169);
    const auto *hh_170 = buffer.data(hh + 170);
    const auto *hh_171 = buffer.data(hh + 171);
    const auto *hh_172 = buffer.data(hh + 172);
    const auto *hh_173 = buffer.data(hh + 173);
    const auto *hh_174 = buffer.data(hh + 174);
    const auto *hh_175 = buffer.data(hh + 175);
    const auto *hh_176 = buffer.data(hh + 176);
    const auto *hh_177 = buffer.data(hh + 177);
    const auto *hh_178 = buffer.data(hh + 178);
    const auto *hh_179 = buffer.data(hh + 179);
    const auto *hh_180 = buffer.data(hh + 180);
    const auto *hh_181 = buffer.data(hh + 181);
    const auto *hh_182 = buffer.data(hh + 182);
    const auto *hh_183 = buffer.data(hh + 183);
    const auto *hh_184 = buffer.data(hh + 184);
    const auto *hh_185 = buffer.data(hh + 185);
    const auto *hh_186 = buffer.data(hh + 186);
    const auto *hh_187 = buffer.data(hh + 187);
    const auto *hh_188 = buffer.data(hh + 188);
    const auto *hh_189 = buffer.data(hh + 189);
    const auto *hh_190 = buffer.data(hh + 190);
    const auto *hh_191 = buffer.data(hh + 191);
    const auto *hh_192 = buffer.data(hh + 192);
    const auto *hh_193 = buffer.data(hh + 193);
    const auto *hh_194 = buffer.data(hh + 194);
    const auto *hh_195 = buffer.data(hh + 195);
    const auto *hh_196 = buffer.data(hh + 196);
    const auto *hh_197 = buffer.data(hh + 197);
    const auto *hh_198 = buffer.data(hh + 198);
    const auto *hh_199 = buffer.data(hh + 199);
    const auto *hh_200 = buffer.data(hh + 200);
    const auto *hh_201 = buffer.data(hh + 201);
    const auto *hh_202 = buffer.data(hh + 202);
    const auto *hh_203 = buffer.data(hh + 203);
    const auto *hh_204 = buffer.data(hh + 204);
    const auto *hh_205 = buffer.data(hh + 205);
    const auto *hh_206 = buffer.data(hh + 206);
    const auto *hh_207 = buffer.data(hh + 207);
    const auto *hh_208 = buffer.data(hh + 208);
    const auto *hh_209 = buffer.data(hh + 209);
    const auto *hh_210 = buffer.data(hh + 210);
    const auto *hh_211 = buffer.data(hh + 211);
    const auto *hh_212 = buffer.data(hh + 212);
    const auto *hh_213 = buffer.data(hh + 213);
    const auto *hh_214 = buffer.data(hh + 214);
    const auto *hh_215 = buffer.data(hh + 215);
    const auto *hh_216 = buffer.data(hh + 216);
    const auto *hh_217 = buffer.data(hh + 217);
    const auto *hh_218 = buffer.data(hh + 218);
    const auto *hh_219 = buffer.data(hh + 219);
    const auto *hh_220 = buffer.data(hh + 220);
    const auto *hh_221 = buffer.data(hh + 221);
    const auto *hh_222 = buffer.data(hh + 222);
    const auto *hh_223 = buffer.data(hh + 223);
    const auto *hh_224 = buffer.data(hh + 224);

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

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, hh_106, hh_107, hh_108, hh_109, \
                         hh_110, kh_253, kh_254, kh_255, kh_256, \
                         kh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -hh_106[k]
                   + f_0 * kh_253[k];

        t_170[k] = -hh_107[k]
                   + f_0 * kh_254[k];

        t_171[k] = -hh_108[k]
                   + f_0 * kh_255[k];

        t_172[k] = -hh_109[k]
                   + f_0 * kh_256[k];

        t_173[k] = -hh_110[k]
                   + f_0 * kh_257[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, hh_111, hh_112, hh_113, hh_114, \
                         hh_115, kh_258, kh_259, kh_260, kh_261, \
                         kh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -hh_111[k]
                   + f_0 * kh_258[k];

        t_175[k] = -hh_112[k]
                   + f_0 * kh_259[k];

        t_176[k] = -hh_113[k]
                   + f_0 * kh_260[k];

        t_177[k] = -hh_114[k]
                   + f_0 * kh_261[k];

        t_178[k] = -hh_115[k]
                   + f_0 * kh_262[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, hh_116, hh_117, hh_118, hh_119, \
                         hh_120, kh_263, kh_264, kh_265, kh_266, \
                         kh_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -hh_116[k]
                   + f_0 * kh_263[k];

        t_180[k] = -hh_117[k]
                   + f_0 * kh_264[k];

        t_181[k] = -hh_118[k]
                   + f_0 * kh_265[k];

        t_182[k] = -hh_119[k]
                   + f_0 * kh_266[k];

        t_183[k] = -hh_120[k]
                   + f_0 * kh_267[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, hh_121, hh_122, hh_123, hh_124, \
                         hh_125, kh_268, kh_269, kh_270, kh_271, \
                         kh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -hh_121[k]
                   + f_0 * kh_268[k];

        t_185[k] = -hh_122[k]
                   + f_0 * kh_269[k];

        t_186[k] = -hh_123[k]
                   + f_0 * kh_270[k];

        t_187[k] = -hh_124[k]
                   + f_0 * kh_271[k];

        t_188[k] = -hh_125[k]
                   + f_0 * kh_272[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, t_196, kh_273, \
                         kh_274, kh_275, kh_276, kh_277, kh_278, kh_279, \
                         kh_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * kh_273[k];

        t_190[k] = f_0 * kh_274[k];

        t_191[k] = f_0 * kh_275[k];

        t_192[k] = f_0 * kh_276[k];

        t_193[k] = f_0 * kh_277[k];

        t_194[k] = f_0 * kh_278[k];

        t_195[k] = f_0 * kh_279[k];

        t_196[k] = f_0 * kh_280[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, t_202, t_203, t_204, kh_281, \
                         kh_282, kh_283, kh_284, kh_285, kh_286, kh_287, \
                         kh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_0 * kh_281[k];

        t_198[k] = f_0 * kh_282[k];

        t_199[k] = f_0 * kh_283[k];

        t_200[k] = f_0 * kh_284[k];

        t_201[k] = f_0 * kh_285[k];

        t_202[k] = f_0 * kh_286[k];

        t_203[k] = f_0 * kh_287[k];

        t_204[k] = f_0 * kh_288[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, t_210, t_211, hh_126, hh_127, \
                         kh_289, kh_290, kh_291, kh_292, kh_293, kh_315, \
                         kh_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_0 * kh_289[k];

        t_206[k] = f_0 * kh_290[k];

        t_207[k] = f_0 * kh_291[k];

        t_208[k] = f_0 * kh_292[k];

        t_209[k] = f_0 * kh_293[k];

        t_210[k] = -4.0 * hh_126[k]
                   + f_0 * kh_315[k];

        t_211[k] = -4.0 * hh_127[k]
                   + f_0 * kh_316[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, hh_128, hh_129, hh_130, hh_131, \
                         hh_132, kh_317, kh_318, kh_319, kh_320, \
                         kh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = -4.0 * hh_128[k]
                   + f_0 * kh_317[k];

        t_213[k] = -4.0 * hh_129[k]
                   + f_0 * kh_318[k];

        t_214[k] = -4.0 * hh_130[k]
                   + f_0 * kh_319[k];

        t_215[k] = -4.0 * hh_131[k]
                   + f_0 * kh_320[k];

        t_216[k] = -4.0 * hh_132[k]
                   + f_0 * kh_321[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, hh_133, hh_134, hh_135, hh_136, \
                         hh_137, kh_322, kh_323, kh_324, kh_325, \
                         kh_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = -4.0 * hh_133[k]
                   + f_0 * kh_322[k];

        t_218[k] = -4.0 * hh_134[k]
                   + f_0 * kh_323[k];

        t_219[k] = -4.0 * hh_135[k]
                   + f_0 * kh_324[k];

        t_220[k] = -4.0 * hh_136[k]
                   + f_0 * kh_325[k];

        t_221[k] = -4.0 * hh_137[k]
                   + f_0 * kh_326[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, hh_138, hh_139, hh_140, hh_141, \
                         hh_142, kh_327, kh_328, kh_329, kh_330, \
                         kh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -4.0 * hh_138[k]
                   + f_0 * kh_327[k];

        t_223[k] = -4.0 * hh_139[k]
                   + f_0 * kh_328[k];

        t_224[k] = -4.0 * hh_140[k]
                   + f_0 * kh_329[k];

        t_225[k] = -4.0 * hh_141[k]
                   + f_0 * kh_330[k];

        t_226[k] = -4.0 * hh_142[k]
                   + f_0 * kh_331[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, t_231, hh_143, hh_144, hh_145, hh_146, \
                         hh_147, kh_332, kh_333, kh_334, kh_335, \
                         kh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = -4.0 * hh_143[k]
                   + f_0 * kh_332[k];

        t_228[k] = -4.0 * hh_144[k]
                   + f_0 * kh_333[k];

        t_229[k] = -4.0 * hh_145[k]
                   + f_0 * kh_334[k];

        t_230[k] = -4.0 * hh_146[k]
                   + f_0 * kh_335[k];

        t_231[k] = -3.0 * hh_147[k]
                   + f_0 * kh_336[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, hh_148, hh_149, hh_150, hh_151, \
                         hh_152, kh_337, kh_338, kh_339, kh_340, \
                         kh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = -3.0 * hh_148[k]
                   + f_0 * kh_337[k];

        t_233[k] = -3.0 * hh_149[k]
                   + f_0 * kh_338[k];

        t_234[k] = -3.0 * hh_150[k]
                   + f_0 * kh_339[k];

        t_235[k] = -3.0 * hh_151[k]
                   + f_0 * kh_340[k];

        t_236[k] = -3.0 * hh_152[k]
                   + f_0 * kh_341[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, hh_153, hh_154, hh_155, hh_156, \
                         hh_157, kh_342, kh_343, kh_344, kh_345, \
                         kh_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = -3.0 * hh_153[k]
                   + f_0 * kh_342[k];

        t_238[k] = -3.0 * hh_154[k]
                   + f_0 * kh_343[k];

        t_239[k] = -3.0 * hh_155[k]
                   + f_0 * kh_344[k];

        t_240[k] = -3.0 * hh_156[k]
                   + f_0 * kh_345[k];

        t_241[k] = -3.0 * hh_157[k]
                   + f_0 * kh_346[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, hh_158, hh_159, hh_160, hh_161, \
                         hh_162, kh_347, kh_348, kh_349, kh_350, \
                         kh_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = -3.0 * hh_158[k]
                   + f_0 * kh_347[k];

        t_243[k] = -3.0 * hh_159[k]
                   + f_0 * kh_348[k];

        t_244[k] = -3.0 * hh_160[k]
                   + f_0 * kh_349[k];

        t_245[k] = -3.0 * hh_161[k]
                   + f_0 * kh_350[k];

        t_246[k] = -3.0 * hh_162[k]
                   + f_0 * kh_351[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, hh_163, hh_164, hh_165, hh_166, \
                         hh_167, kh_352, kh_353, kh_354, kh_355, \
                         kh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = -3.0 * hh_163[k]
                   + f_0 * kh_352[k];

        t_248[k] = -3.0 * hh_164[k]
                   + f_0 * kh_353[k];

        t_249[k] = -3.0 * hh_165[k]
                   + f_0 * kh_354[k];

        t_250[k] = -3.0 * hh_166[k]
                   + f_0 * kh_355[k];

        t_251[k] = -3.0 * hh_167[k]
                   + f_0 * kh_356[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, hh_168, hh_169, hh_170, hh_171, \
                         hh_172, kh_357, kh_358, kh_359, kh_360, \
                         kh_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -2.0 * hh_168[k]
                   + f_0 * kh_357[k];

        t_253[k] = -2.0 * hh_169[k]
                   + f_0 * kh_358[k];

        t_254[k] = -2.0 * hh_170[k]
                   + f_0 * kh_359[k];

        t_255[k] = -2.0 * hh_171[k]
                   + f_0 * kh_360[k];

        t_256[k] = -2.0 * hh_172[k]
                   + f_0 * kh_361[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, hh_173, hh_174, hh_175, hh_176, \
                         hh_177, kh_362, kh_363, kh_364, kh_365, \
                         kh_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = -2.0 * hh_173[k]
                   + f_0 * kh_362[k];

        t_258[k] = -2.0 * hh_174[k]
                   + f_0 * kh_363[k];

        t_259[k] = -2.0 * hh_175[k]
                   + f_0 * kh_364[k];

        t_260[k] = -2.0 * hh_176[k]
                   + f_0 * kh_365[k];

        t_261[k] = -2.0 * hh_177[k]
                   + f_0 * kh_366[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, hh_178, hh_179, hh_180, hh_181, \
                         hh_182, kh_367, kh_368, kh_369, kh_370, \
                         kh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = -2.0 * hh_178[k]
                   + f_0 * kh_367[k];

        t_263[k] = -2.0 * hh_179[k]
                   + f_0 * kh_368[k];

        t_264[k] = -2.0 * hh_180[k]
                   + f_0 * kh_369[k];

        t_265[k] = -2.0 * hh_181[k]
                   + f_0 * kh_370[k];

        t_266[k] = -2.0 * hh_182[k]
                   + f_0 * kh_371[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, hh_183, hh_184, hh_185, hh_186, \
                         hh_187, kh_372, kh_373, kh_374, kh_375, \
                         kh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = -2.0 * hh_183[k]
                   + f_0 * kh_372[k];

        t_268[k] = -2.0 * hh_184[k]
                   + f_0 * kh_373[k];

        t_269[k] = -2.0 * hh_185[k]
                   + f_0 * kh_374[k];

        t_270[k] = -2.0 * hh_186[k]
                   + f_0 * kh_375[k];

        t_271[k] = -2.0 * hh_187[k]
                   + f_0 * kh_376[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, hh_188, hh_189, hh_190, hh_191, \
                         hh_192, kh_377, kh_378, kh_379, kh_380, \
                         kh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -2.0 * hh_188[k]
                   + f_0 * kh_377[k];

        t_273[k] = -hh_189[k]
                   + f_0 * kh_378[k];

        t_274[k] = -hh_190[k]
                   + f_0 * kh_379[k];

        t_275[k] = -hh_191[k]
                   + f_0 * kh_380[k];

        t_276[k] = -hh_192[k]
                   + f_0 * kh_381[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, t_281, hh_193, hh_194, hh_195, hh_196, \
                         hh_197, kh_382, kh_383, kh_384, kh_385, \
                         kh_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = -hh_193[k]
                   + f_0 * kh_382[k];

        t_278[k] = -hh_194[k]
                   + f_0 * kh_383[k];

        t_279[k] = -hh_195[k]
                   + f_0 * kh_384[k];

        t_280[k] = -hh_196[k]
                   + f_0 * kh_385[k];

        t_281[k] = -hh_197[k]
                   + f_0 * kh_386[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, hh_198, hh_199, hh_200, hh_201, \
                         hh_202, kh_387, kh_388, kh_389, kh_390, \
                         kh_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = -hh_198[k]
                   + f_0 * kh_387[k];

        t_283[k] = -hh_199[k]
                   + f_0 * kh_388[k];

        t_284[k] = -hh_200[k]
                   + f_0 * kh_389[k];

        t_285[k] = -hh_201[k]
                   + f_0 * kh_390[k];

        t_286[k] = -hh_202[k]
                   + f_0 * kh_391[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, hh_203, hh_204, hh_205, hh_206, \
                         hh_207, kh_392, kh_393, kh_394, kh_395, \
                         kh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = -hh_203[k]
                   + f_0 * kh_392[k];

        t_288[k] = -hh_204[k]
                   + f_0 * kh_393[k];

        t_289[k] = -hh_205[k]
                   + f_0 * kh_394[k];

        t_290[k] = -hh_206[k]
                   + f_0 * kh_395[k];

        t_291[k] = -hh_207[k]
                   + f_0 * kh_396[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, t_297, t_298, hh_208, hh_209, \
                         kh_397, kh_398, kh_399, kh_400, kh_401, kh_402, \
                         kh_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = -hh_208[k]
                   + f_0 * kh_397[k];

        t_293[k] = -hh_209[k]
                   + f_0 * kh_398[k];

        t_294[k] = f_0 * kh_399[k];

        t_295[k] = f_0 * kh_400[k];

        t_296[k] = f_0 * kh_401[k];

        t_297[k] = f_0 * kh_402[k];

        t_298[k] = f_0 * kh_403[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, t_304, t_305, t_306, kh_404, \
                         kh_405, kh_406, kh_407, kh_408, kh_409, kh_410, \
                         kh_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_0 * kh_404[k];

        t_300[k] = f_0 * kh_405[k];

        t_301[k] = f_0 * kh_406[k];

        t_302[k] = f_0 * kh_407[k];

        t_303[k] = f_0 * kh_408[k];

        t_304[k] = f_0 * kh_409[k];

        t_305[k] = f_0 * kh_410[k];

        t_306[k] = f_0 * kh_411[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, t_312, t_313, t_314, kh_412, \
                         kh_413, kh_414, kh_415, kh_416, kh_417, kh_418, \
                         kh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_0 * kh_412[k];

        t_308[k] = f_0 * kh_413[k];

        t_309[k] = f_0 * kh_414[k];

        t_310[k] = f_0 * kh_415[k];

        t_311[k] = f_0 * kh_416[k];

        t_312[k] = f_0 * kh_417[k];

        t_313[k] = f_0 * kh_418[k];

        t_314[k] = f_0 * kh_419[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, hh_210, hh_211, hh_212, hh_213, \
                         hh_214, kh_441, kh_442, kh_443, kh_444, \
                         kh_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -5.0 * hh_210[k]
                   + f_0 * kh_441[k];

        t_316[k] = -5.0 * hh_211[k]
                   + f_0 * kh_442[k];

        t_317[k] = -5.0 * hh_212[k]
                   + f_0 * kh_443[k];

        t_318[k] = -5.0 * hh_213[k]
                   + f_0 * kh_444[k];

        t_319[k] = -5.0 * hh_214[k]
                   + f_0 * kh_445[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, hh_215, hh_216, hh_217, hh_218, \
                         hh_219, kh_446, kh_447, kh_448, kh_449, \
                         kh_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -5.0 * hh_215[k]
                   + f_0 * kh_446[k];

        t_321[k] = -5.0 * hh_216[k]
                   + f_0 * kh_447[k];

        t_322[k] = -5.0 * hh_217[k]
                   + f_0 * kh_448[k];

        t_323[k] = -5.0 * hh_218[k]
                   + f_0 * kh_449[k];

        t_324[k] = -5.0 * hh_219[k]
                   + f_0 * kh_450[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, hh_220, hh_221, hh_222, hh_223, \
                         hh_224, kh_451, kh_452, kh_453, kh_454, \
                         kh_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -5.0 * hh_220[k]
                   + f_0 * kh_451[k];

        t_326[k] = -5.0 * hh_221[k]
                   + f_0 * kh_452[k];

        t_327[k] = -5.0 * hh_222[k]
                   + f_0 * kh_453[k];

        t_328[k] = -5.0 * hh_223[k]
                   + f_0 * kh_454[k];

        t_329[k] = -5.0 * hh_224[k]
                   + f_0 * kh_455[k];
    }
}

static auto
compute_prim_geom_10_ih_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hh, const size_t kh,
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

    const auto *hh_225 = buffer.data(hh + 225);
    const auto *hh_226 = buffer.data(hh + 226);
    const auto *hh_227 = buffer.data(hh + 227);
    const auto *hh_228 = buffer.data(hh + 228);
    const auto *hh_229 = buffer.data(hh + 229);
    const auto *hh_230 = buffer.data(hh + 230);
    const auto *hh_231 = buffer.data(hh + 231);
    const auto *hh_232 = buffer.data(hh + 232);
    const auto *hh_233 = buffer.data(hh + 233);
    const auto *hh_234 = buffer.data(hh + 234);
    const auto *hh_235 = buffer.data(hh + 235);
    const auto *hh_236 = buffer.data(hh + 236);
    const auto *hh_237 = buffer.data(hh + 237);
    const auto *hh_238 = buffer.data(hh + 238);
    const auto *hh_239 = buffer.data(hh + 239);
    const auto *hh_240 = buffer.data(hh + 240);
    const auto *hh_241 = buffer.data(hh + 241);
    const auto *hh_242 = buffer.data(hh + 242);
    const auto *hh_243 = buffer.data(hh + 243);
    const auto *hh_244 = buffer.data(hh + 244);
    const auto *hh_245 = buffer.data(hh + 245);
    const auto *hh_246 = buffer.data(hh + 246);
    const auto *hh_247 = buffer.data(hh + 247);
    const auto *hh_248 = buffer.data(hh + 248);
    const auto *hh_249 = buffer.data(hh + 249);
    const auto *hh_250 = buffer.data(hh + 250);
    const auto *hh_251 = buffer.data(hh + 251);
    const auto *hh_252 = buffer.data(hh + 252);
    const auto *hh_253 = buffer.data(hh + 253);
    const auto *hh_254 = buffer.data(hh + 254);
    const auto *hh_255 = buffer.data(hh + 255);
    const auto *hh_256 = buffer.data(hh + 256);
    const auto *hh_257 = buffer.data(hh + 257);
    const auto *hh_258 = buffer.data(hh + 258);
    const auto *hh_259 = buffer.data(hh + 259);
    const auto *hh_260 = buffer.data(hh + 260);
    const auto *hh_261 = buffer.data(hh + 261);
    const auto *hh_262 = buffer.data(hh + 262);
    const auto *hh_263 = buffer.data(hh + 263);
    const auto *hh_264 = buffer.data(hh + 264);
    const auto *hh_265 = buffer.data(hh + 265);
    const auto *hh_266 = buffer.data(hh + 266);
    const auto *hh_267 = buffer.data(hh + 267);
    const auto *hh_268 = buffer.data(hh + 268);
    const auto *hh_269 = buffer.data(hh + 269);
    const auto *hh_270 = buffer.data(hh + 270);
    const auto *hh_271 = buffer.data(hh + 271);
    const auto *hh_272 = buffer.data(hh + 272);
    const auto *hh_273 = buffer.data(hh + 273);
    const auto *hh_274 = buffer.data(hh + 274);
    const auto *hh_275 = buffer.data(hh + 275);
    const auto *hh_276 = buffer.data(hh + 276);
    const auto *hh_277 = buffer.data(hh + 277);
    const auto *hh_278 = buffer.data(hh + 278);
    const auto *hh_279 = buffer.data(hh + 279);
    const auto *hh_280 = buffer.data(hh + 280);
    const auto *hh_281 = buffer.data(hh + 281);
    const auto *hh_282 = buffer.data(hh + 282);
    const auto *hh_283 = buffer.data(hh + 283);
    const auto *hh_284 = buffer.data(hh + 284);
    const auto *hh_285 = buffer.data(hh + 285);
    const auto *hh_286 = buffer.data(hh + 286);
    const auto *hh_287 = buffer.data(hh + 287);
    const auto *hh_288 = buffer.data(hh + 288);
    const auto *hh_289 = buffer.data(hh + 289);
    const auto *hh_290 = buffer.data(hh + 290);
    const auto *hh_291 = buffer.data(hh + 291);
    const auto *hh_292 = buffer.data(hh + 292);
    const auto *hh_293 = buffer.data(hh + 293);
    const auto *hh_294 = buffer.data(hh + 294);
    const auto *hh_295 = buffer.data(hh + 295);
    const auto *hh_296 = buffer.data(hh + 296);
    const auto *hh_297 = buffer.data(hh + 297);
    const auto *hh_298 = buffer.data(hh + 298);
    const auto *hh_299 = buffer.data(hh + 299);
    const auto *hh_300 = buffer.data(hh + 300);
    const auto *hh_301 = buffer.data(hh + 301);
    const auto *hh_302 = buffer.data(hh + 302);
    const auto *hh_303 = buffer.data(hh + 303);
    const auto *hh_304 = buffer.data(hh + 304);
    const auto *hh_305 = buffer.data(hh + 305);
    const auto *hh_306 = buffer.data(hh + 306);
    const auto *hh_307 = buffer.data(hh + 307);
    const auto *hh_308 = buffer.data(hh + 308);
    const auto *hh_309 = buffer.data(hh + 309);
    const auto *hh_310 = buffer.data(hh + 310);
    const auto *hh_311 = buffer.data(hh + 311);
    const auto *hh_312 = buffer.data(hh + 312);
    const auto *hh_313 = buffer.data(hh + 313);
    const auto *hh_314 = buffer.data(hh + 314);
    const auto *hh_315 = buffer.data(hh + 315);
    const auto *hh_316 = buffer.data(hh + 316);
    const auto *hh_317 = buffer.data(hh + 317);
    const auto *hh_318 = buffer.data(hh + 318);
    const auto *hh_319 = buffer.data(hh + 319);
    const auto *hh_320 = buffer.data(hh + 320);
    const auto *hh_321 = buffer.data(hh + 321);
    const auto *hh_322 = buffer.data(hh + 322);
    const auto *hh_323 = buffer.data(hh + 323);
    const auto *hh_324 = buffer.data(hh + 324);
    const auto *hh_325 = buffer.data(hh + 325);
    const auto *hh_326 = buffer.data(hh + 326);
    const auto *hh_327 = buffer.data(hh + 327);
    const auto *hh_328 = buffer.data(hh + 328);
    const auto *hh_329 = buffer.data(hh + 329);
    const auto *hh_330 = buffer.data(hh + 330);
    const auto *hh_331 = buffer.data(hh + 331);
    const auto *hh_332 = buffer.data(hh + 332);
    const auto *hh_333 = buffer.data(hh + 333);
    const auto *hh_334 = buffer.data(hh + 334);
    const auto *hh_335 = buffer.data(hh + 335);
    const auto *hh_336 = buffer.data(hh + 336);
    const auto *hh_337 = buffer.data(hh + 337);
    const auto *hh_338 = buffer.data(hh + 338);
    const auto *hh_339 = buffer.data(hh + 339);
    const auto *hh_340 = buffer.data(hh + 340);
    const auto *hh_341 = buffer.data(hh + 341);
    const auto *hh_342 = buffer.data(hh + 342);
    const auto *hh_343 = buffer.data(hh + 343);
    const auto *hh_344 = buffer.data(hh + 344);
    const auto *hh_345 = buffer.data(hh + 345);
    const auto *hh_346 = buffer.data(hh + 346);
    const auto *hh_347 = buffer.data(hh + 347);
    const auto *hh_348 = buffer.data(hh + 348);
    const auto *hh_349 = buffer.data(hh + 349);
    const auto *hh_350 = buffer.data(hh + 350);
    const auto *hh_351 = buffer.data(hh + 351);
    const auto *hh_352 = buffer.data(hh + 352);
    const auto *hh_353 = buffer.data(hh + 353);
    const auto *hh_354 = buffer.data(hh + 354);
    const auto *hh_355 = buffer.data(hh + 355);
    const auto *hh_356 = buffer.data(hh + 356);
    const auto *hh_357 = buffer.data(hh + 357);
    const auto *hh_358 = buffer.data(hh + 358);
    const auto *hh_359 = buffer.data(hh + 359);
    const auto *hh_360 = buffer.data(hh + 360);
    const auto *hh_361 = buffer.data(hh + 361);

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

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, hh_225, hh_226, hh_227, hh_228, \
                         hh_229, kh_456, kh_457, kh_458, kh_459, \
                         kh_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -5.0 * hh_225[k]
                   + f_0 * kh_456[k];

        t_331[k] = -5.0 * hh_226[k]
                   + f_0 * kh_457[k];

        t_332[k] = -5.0 * hh_227[k]
                   + f_0 * kh_458[k];

        t_333[k] = -5.0 * hh_228[k]
                   + f_0 * kh_459[k];

        t_334[k] = -5.0 * hh_229[k]
                   + f_0 * kh_460[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, hh_230, hh_231, hh_232, hh_233, \
                         hh_234, kh_461, kh_462, kh_463, kh_464, \
                         kh_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -5.0 * hh_230[k]
                   + f_0 * kh_461[k];

        t_336[k] = -4.0 * hh_231[k]
                   + f_0 * kh_462[k];

        t_337[k] = -4.0 * hh_232[k]
                   + f_0 * kh_463[k];

        t_338[k] = -4.0 * hh_233[k]
                   + f_0 * kh_464[k];

        t_339[k] = -4.0 * hh_234[k]
                   + f_0 * kh_465[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, hh_235, hh_236, hh_237, hh_238, \
                         hh_239, kh_466, kh_467, kh_468, kh_469, \
                         kh_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -4.0 * hh_235[k]
                   + f_0 * kh_466[k];

        t_341[k] = -4.0 * hh_236[k]
                   + f_0 * kh_467[k];

        t_342[k] = -4.0 * hh_237[k]
                   + f_0 * kh_468[k];

        t_343[k] = -4.0 * hh_238[k]
                   + f_0 * kh_469[k];

        t_344[k] = -4.0 * hh_239[k]
                   + f_0 * kh_470[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, hh_240, hh_241, hh_242, hh_243, \
                         hh_244, kh_471, kh_472, kh_473, kh_474, \
                         kh_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -4.0 * hh_240[k]
                   + f_0 * kh_471[k];

        t_346[k] = -4.0 * hh_241[k]
                   + f_0 * kh_472[k];

        t_347[k] = -4.0 * hh_242[k]
                   + f_0 * kh_473[k];

        t_348[k] = -4.0 * hh_243[k]
                   + f_0 * kh_474[k];

        t_349[k] = -4.0 * hh_244[k]
                   + f_0 * kh_475[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, hh_245, hh_246, hh_247, hh_248, \
                         hh_249, kh_476, kh_477, kh_478, kh_479, \
                         kh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -4.0 * hh_245[k]
                   + f_0 * kh_476[k];

        t_351[k] = -4.0 * hh_246[k]
                   + f_0 * kh_477[k];

        t_352[k] = -4.0 * hh_247[k]
                   + f_0 * kh_478[k];

        t_353[k] = -4.0 * hh_248[k]
                   + f_0 * kh_479[k];

        t_354[k] = -4.0 * hh_249[k]
                   + f_0 * kh_480[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, hh_250, hh_251, hh_252, hh_253, \
                         hh_254, kh_481, kh_482, kh_483, kh_484, \
                         kh_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -4.0 * hh_250[k]
                   + f_0 * kh_481[k];

        t_356[k] = -4.0 * hh_251[k]
                   + f_0 * kh_482[k];

        t_357[k] = -3.0 * hh_252[k]
                   + f_0 * kh_483[k];

        t_358[k] = -3.0 * hh_253[k]
                   + f_0 * kh_484[k];

        t_359[k] = -3.0 * hh_254[k]
                   + f_0 * kh_485[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, hh_255, hh_256, hh_257, hh_258, \
                         hh_259, kh_486, kh_487, kh_488, kh_489, \
                         kh_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -3.0 * hh_255[k]
                   + f_0 * kh_486[k];

        t_361[k] = -3.0 * hh_256[k]
                   + f_0 * kh_487[k];

        t_362[k] = -3.0 * hh_257[k]
                   + f_0 * kh_488[k];

        t_363[k] = -3.0 * hh_258[k]
                   + f_0 * kh_489[k];

        t_364[k] = -3.0 * hh_259[k]
                   + f_0 * kh_490[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, hh_260, hh_261, hh_262, hh_263, \
                         hh_264, kh_491, kh_492, kh_493, kh_494, \
                         kh_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -3.0 * hh_260[k]
                   + f_0 * kh_491[k];

        t_366[k] = -3.0 * hh_261[k]
                   + f_0 * kh_492[k];

        t_367[k] = -3.0 * hh_262[k]
                   + f_0 * kh_493[k];

        t_368[k] = -3.0 * hh_263[k]
                   + f_0 * kh_494[k];

        t_369[k] = -3.0 * hh_264[k]
                   + f_0 * kh_495[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, hh_265, hh_266, hh_267, hh_268, \
                         hh_269, kh_496, kh_497, kh_498, kh_499, \
                         kh_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -3.0 * hh_265[k]
                   + f_0 * kh_496[k];

        t_371[k] = -3.0 * hh_266[k]
                   + f_0 * kh_497[k];

        t_372[k] = -3.0 * hh_267[k]
                   + f_0 * kh_498[k];

        t_373[k] = -3.0 * hh_268[k]
                   + f_0 * kh_499[k];

        t_374[k] = -3.0 * hh_269[k]
                   + f_0 * kh_500[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, hh_270, hh_271, hh_272, hh_273, \
                         hh_274, kh_501, kh_502, kh_503, kh_504, \
                         kh_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -3.0 * hh_270[k]
                   + f_0 * kh_501[k];

        t_376[k] = -3.0 * hh_271[k]
                   + f_0 * kh_502[k];

        t_377[k] = -3.0 * hh_272[k]
                   + f_0 * kh_503[k];

        t_378[k] = -2.0 * hh_273[k]
                   + f_0 * kh_504[k];

        t_379[k] = -2.0 * hh_274[k]
                   + f_0 * kh_505[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, hh_275, hh_276, hh_277, hh_278, \
                         hh_279, kh_506, kh_507, kh_508, kh_509, \
                         kh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * hh_275[k]
                   + f_0 * kh_506[k];

        t_381[k] = -2.0 * hh_276[k]
                   + f_0 * kh_507[k];

        t_382[k] = -2.0 * hh_277[k]
                   + f_0 * kh_508[k];

        t_383[k] = -2.0 * hh_278[k]
                   + f_0 * kh_509[k];

        t_384[k] = -2.0 * hh_279[k]
                   + f_0 * kh_510[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, hh_280, hh_281, hh_282, hh_283, \
                         hh_284, kh_511, kh_512, kh_513, kh_514, \
                         kh_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * hh_280[k]
                   + f_0 * kh_511[k];

        t_386[k] = -2.0 * hh_281[k]
                   + f_0 * kh_512[k];

        t_387[k] = -2.0 * hh_282[k]
                   + f_0 * kh_513[k];

        t_388[k] = -2.0 * hh_283[k]
                   + f_0 * kh_514[k];

        t_389[k] = -2.0 * hh_284[k]
                   + f_0 * kh_515[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, hh_285, hh_286, hh_287, hh_288, \
                         hh_289, kh_516, kh_517, kh_518, kh_519, \
                         kh_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -2.0 * hh_285[k]
                   + f_0 * kh_516[k];

        t_391[k] = -2.0 * hh_286[k]
                   + f_0 * kh_517[k];

        t_392[k] = -2.0 * hh_287[k]
                   + f_0 * kh_518[k];

        t_393[k] = -2.0 * hh_288[k]
                   + f_0 * kh_519[k];

        t_394[k] = -2.0 * hh_289[k]
                   + f_0 * kh_520[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, hh_290, hh_291, hh_292, hh_293, \
                         hh_294, kh_521, kh_522, kh_523, kh_524, \
                         kh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -2.0 * hh_290[k]
                   + f_0 * kh_521[k];

        t_396[k] = -2.0 * hh_291[k]
                   + f_0 * kh_522[k];

        t_397[k] = -2.0 * hh_292[k]
                   + f_0 * kh_523[k];

        t_398[k] = -2.0 * hh_293[k]
                   + f_0 * kh_524[k];

        t_399[k] = -hh_294[k]
                   + f_0 * kh_525[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, hh_295, hh_296, hh_297, hh_298, \
                         hh_299, kh_526, kh_527, kh_528, kh_529, \
                         kh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -hh_295[k]
                   + f_0 * kh_526[k];

        t_401[k] = -hh_296[k]
                   + f_0 * kh_527[k];

        t_402[k] = -hh_297[k]
                   + f_0 * kh_528[k];

        t_403[k] = -hh_298[k]
                   + f_0 * kh_529[k];

        t_404[k] = -hh_299[k]
                   + f_0 * kh_530[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, hh_300, hh_301, hh_302, hh_303, \
                         hh_304, kh_531, kh_532, kh_533, kh_534, \
                         kh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -hh_300[k]
                   + f_0 * kh_531[k];

        t_406[k] = -hh_301[k]
                   + f_0 * kh_532[k];

        t_407[k] = -hh_302[k]
                   + f_0 * kh_533[k];

        t_408[k] = -hh_303[k]
                   + f_0 * kh_534[k];

        t_409[k] = -hh_304[k]
                   + f_0 * kh_535[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, hh_305, hh_306, hh_307, hh_308, \
                         hh_309, kh_536, kh_537, kh_538, kh_539, \
                         kh_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -hh_305[k]
                   + f_0 * kh_536[k];

        t_411[k] = -hh_306[k]
                   + f_0 * kh_537[k];

        t_412[k] = -hh_307[k]
                   + f_0 * kh_538[k];

        t_413[k] = -hh_308[k]
                   + f_0 * kh_539[k];

        t_414[k] = -hh_309[k]
                   + f_0 * kh_540[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, hh_310, hh_311, hh_312, hh_313, \
                         hh_314, kh_541, kh_542, kh_543, kh_544, \
                         kh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -hh_310[k]
                   + f_0 * kh_541[k];

        t_416[k] = -hh_311[k]
                   + f_0 * kh_542[k];

        t_417[k] = -hh_312[k]
                   + f_0 * kh_543[k];

        t_418[k] = -hh_313[k]
                   + f_0 * kh_544[k];

        t_419[k] = -hh_314[k]
                   + f_0 * kh_545[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, t_425, t_426, t_427, kh_546, \
                         kh_547, kh_548, kh_549, kh_550, kh_551, kh_552, \
                         kh_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_0 * kh_546[k];

        t_421[k] = f_0 * kh_547[k];

        t_422[k] = f_0 * kh_548[k];

        t_423[k] = f_0 * kh_549[k];

        t_424[k] = f_0 * kh_550[k];

        t_425[k] = f_0 * kh_551[k];

        t_426[k] = f_0 * kh_552[k];

        t_427[k] = f_0 * kh_553[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, t_433, t_434, t_435, kh_554, \
                         kh_555, kh_556, kh_557, kh_558, kh_559, kh_560, \
                         kh_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_0 * kh_554[k];

        t_429[k] = f_0 * kh_555[k];

        t_430[k] = f_0 * kh_556[k];

        t_431[k] = f_0 * kh_557[k];

        t_432[k] = f_0 * kh_558[k];

        t_433[k] = f_0 * kh_559[k];

        t_434[k] = f_0 * kh_560[k];

        t_435[k] = f_0 * kh_561[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, t_441, t_442, hh_315, hh_316, \
                         kh_562, kh_563, kh_564, kh_565, kh_566, kh_588, \
                         kh_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_0 * kh_562[k];

        t_437[k] = f_0 * kh_563[k];

        t_438[k] = f_0 * kh_564[k];

        t_439[k] = f_0 * kh_565[k];

        t_440[k] = f_0 * kh_566[k];

        t_441[k] = -6.0 * hh_315[k]
                   + f_0 * kh_588[k];

        t_442[k] = -6.0 * hh_316[k]
                   + f_0 * kh_589[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, t_447, hh_317, hh_318, hh_319, hh_320, \
                         hh_321, kh_590, kh_591, kh_592, kh_593, \
                         kh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = -6.0 * hh_317[k]
                   + f_0 * kh_590[k];

        t_444[k] = -6.0 * hh_318[k]
                   + f_0 * kh_591[k];

        t_445[k] = -6.0 * hh_319[k]
                   + f_0 * kh_592[k];

        t_446[k] = -6.0 * hh_320[k]
                   + f_0 * kh_593[k];

        t_447[k] = -6.0 * hh_321[k]
                   + f_0 * kh_594[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, t_452, hh_322, hh_323, hh_324, hh_325, \
                         hh_326, kh_595, kh_596, kh_597, kh_598, \
                         kh_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = -6.0 * hh_322[k]
                   + f_0 * kh_595[k];

        t_449[k] = -6.0 * hh_323[k]
                   + f_0 * kh_596[k];

        t_450[k] = -6.0 * hh_324[k]
                   + f_0 * kh_597[k];

        t_451[k] = -6.0 * hh_325[k]
                   + f_0 * kh_598[k];

        t_452[k] = -6.0 * hh_326[k]
                   + f_0 * kh_599[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, t_457, hh_327, hh_328, hh_329, hh_330, \
                         hh_331, kh_600, kh_601, kh_602, kh_603, \
                         kh_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = -6.0 * hh_327[k]
                   + f_0 * kh_600[k];

        t_454[k] = -6.0 * hh_328[k]
                   + f_0 * kh_601[k];

        t_455[k] = -6.0 * hh_329[k]
                   + f_0 * kh_602[k];

        t_456[k] = -6.0 * hh_330[k]
                   + f_0 * kh_603[k];

        t_457[k] = -6.0 * hh_331[k]
                   + f_0 * kh_604[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, t_462, hh_332, hh_333, hh_334, hh_335, \
                         hh_336, kh_605, kh_606, kh_607, kh_608, \
                         kh_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = -6.0 * hh_332[k]
                   + f_0 * kh_605[k];

        t_459[k] = -6.0 * hh_333[k]
                   + f_0 * kh_606[k];

        t_460[k] = -6.0 * hh_334[k]
                   + f_0 * kh_607[k];

        t_461[k] = -6.0 * hh_335[k]
                   + f_0 * kh_608[k];

        t_462[k] = -5.0 * hh_336[k]
                   + f_0 * kh_609[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, t_467, hh_337, hh_338, hh_339, hh_340, \
                         hh_341, kh_610, kh_611, kh_612, kh_613, \
                         kh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = -5.0 * hh_337[k]
                   + f_0 * kh_610[k];

        t_464[k] = -5.0 * hh_338[k]
                   + f_0 * kh_611[k];

        t_465[k] = -5.0 * hh_339[k]
                   + f_0 * kh_612[k];

        t_466[k] = -5.0 * hh_340[k]
                   + f_0 * kh_613[k];

        t_467[k] = -5.0 * hh_341[k]
                   + f_0 * kh_614[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, t_471, t_472, hh_342, hh_343, hh_344, hh_345, \
                         hh_346, kh_615, kh_616, kh_617, kh_618, \
                         kh_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = -5.0 * hh_342[k]
                   + f_0 * kh_615[k];

        t_469[k] = -5.0 * hh_343[k]
                   + f_0 * kh_616[k];

        t_470[k] = -5.0 * hh_344[k]
                   + f_0 * kh_617[k];

        t_471[k] = -5.0 * hh_345[k]
                   + f_0 * kh_618[k];

        t_472[k] = -5.0 * hh_346[k]
                   + f_0 * kh_619[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, t_476, t_477, hh_347, hh_348, hh_349, hh_350, \
                         hh_351, kh_620, kh_621, kh_622, kh_623, \
                         kh_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = -5.0 * hh_347[k]
                   + f_0 * kh_620[k];

        t_474[k] = -5.0 * hh_348[k]
                   + f_0 * kh_621[k];

        t_475[k] = -5.0 * hh_349[k]
                   + f_0 * kh_622[k];

        t_476[k] = -5.0 * hh_350[k]
                   + f_0 * kh_623[k];

        t_477[k] = -5.0 * hh_351[k]
                   + f_0 * kh_624[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, t_482, hh_352, hh_353, hh_354, hh_355, \
                         hh_356, kh_625, kh_626, kh_627, kh_628, \
                         kh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = -5.0 * hh_352[k]
                   + f_0 * kh_625[k];

        t_479[k] = -5.0 * hh_353[k]
                   + f_0 * kh_626[k];

        t_480[k] = -5.0 * hh_354[k]
                   + f_0 * kh_627[k];

        t_481[k] = -5.0 * hh_355[k]
                   + f_0 * kh_628[k];

        t_482[k] = -5.0 * hh_356[k]
                   + f_0 * kh_629[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, t_487, hh_357, hh_358, hh_359, hh_360, \
                         hh_361, kh_630, kh_631, kh_632, kh_633, \
                         kh_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = -4.0 * hh_357[k]
                   + f_0 * kh_630[k];

        t_484[k] = -4.0 * hh_358[k]
                   + f_0 * kh_631[k];

        t_485[k] = -4.0 * hh_359[k]
                   + f_0 * kh_632[k];

        t_486[k] = -4.0 * hh_360[k]
                   + f_0 * kh_633[k];

        t_487[k] = -4.0 * hh_361[k]
                   + f_0 * kh_634[k];
    }
}

static auto
compute_prim_geom_10_ih_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hh, const size_t kh,
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

    const auto *hh_362 = buffer.data(hh + 362);
    const auto *hh_363 = buffer.data(hh + 363);
    const auto *hh_364 = buffer.data(hh + 364);
    const auto *hh_365 = buffer.data(hh + 365);
    const auto *hh_366 = buffer.data(hh + 366);
    const auto *hh_367 = buffer.data(hh + 367);
    const auto *hh_368 = buffer.data(hh + 368);
    const auto *hh_369 = buffer.data(hh + 369);
    const auto *hh_370 = buffer.data(hh + 370);
    const auto *hh_371 = buffer.data(hh + 371);
    const auto *hh_372 = buffer.data(hh + 372);
    const auto *hh_373 = buffer.data(hh + 373);
    const auto *hh_374 = buffer.data(hh + 374);
    const auto *hh_375 = buffer.data(hh + 375);
    const auto *hh_376 = buffer.data(hh + 376);
    const auto *hh_377 = buffer.data(hh + 377);
    const auto *hh_378 = buffer.data(hh + 378);
    const auto *hh_379 = buffer.data(hh + 379);
    const auto *hh_380 = buffer.data(hh + 380);
    const auto *hh_381 = buffer.data(hh + 381);
    const auto *hh_382 = buffer.data(hh + 382);
    const auto *hh_383 = buffer.data(hh + 383);
    const auto *hh_384 = buffer.data(hh + 384);
    const auto *hh_385 = buffer.data(hh + 385);
    const auto *hh_386 = buffer.data(hh + 386);
    const auto *hh_387 = buffer.data(hh + 387);
    const auto *hh_388 = buffer.data(hh + 388);
    const auto *hh_389 = buffer.data(hh + 389);
    const auto *hh_390 = buffer.data(hh + 390);
    const auto *hh_391 = buffer.data(hh + 391);
    const auto *hh_392 = buffer.data(hh + 392);
    const auto *hh_393 = buffer.data(hh + 393);
    const auto *hh_394 = buffer.data(hh + 394);
    const auto *hh_395 = buffer.data(hh + 395);
    const auto *hh_396 = buffer.data(hh + 396);
    const auto *hh_397 = buffer.data(hh + 397);
    const auto *hh_398 = buffer.data(hh + 398);
    const auto *hh_399 = buffer.data(hh + 399);
    const auto *hh_400 = buffer.data(hh + 400);
    const auto *hh_401 = buffer.data(hh + 401);
    const auto *hh_402 = buffer.data(hh + 402);
    const auto *hh_403 = buffer.data(hh + 403);
    const auto *hh_404 = buffer.data(hh + 404);
    const auto *hh_405 = buffer.data(hh + 405);
    const auto *hh_406 = buffer.data(hh + 406);
    const auto *hh_407 = buffer.data(hh + 407);
    const auto *hh_408 = buffer.data(hh + 408);
    const auto *hh_409 = buffer.data(hh + 409);
    const auto *hh_410 = buffer.data(hh + 410);
    const auto *hh_411 = buffer.data(hh + 411);
    const auto *hh_412 = buffer.data(hh + 412);
    const auto *hh_413 = buffer.data(hh + 413);
    const auto *hh_414 = buffer.data(hh + 414);
    const auto *hh_415 = buffer.data(hh + 415);
    const auto *hh_416 = buffer.data(hh + 416);
    const auto *hh_417 = buffer.data(hh + 417);
    const auto *hh_418 = buffer.data(hh + 418);
    const auto *hh_419 = buffer.data(hh + 419);
    const auto *hh_420 = buffer.data(hh + 420);
    const auto *hh_421 = buffer.data(hh + 421);
    const auto *hh_422 = buffer.data(hh + 422);
    const auto *hh_423 = buffer.data(hh + 423);
    const auto *hh_424 = buffer.data(hh + 424);
    const auto *hh_425 = buffer.data(hh + 425);
    const auto *hh_426 = buffer.data(hh + 426);
    const auto *hh_427 = buffer.data(hh + 427);
    const auto *hh_428 = buffer.data(hh + 428);
    const auto *hh_429 = buffer.data(hh + 429);
    const auto *hh_430 = buffer.data(hh + 430);
    const auto *hh_431 = buffer.data(hh + 431);
    const auto *hh_432 = buffer.data(hh + 432);
    const auto *hh_433 = buffer.data(hh + 433);
    const auto *hh_434 = buffer.data(hh + 434);
    const auto *hh_435 = buffer.data(hh + 435);
    const auto *hh_436 = buffer.data(hh + 436);
    const auto *hh_437 = buffer.data(hh + 437);
    const auto *hh_438 = buffer.data(hh + 438);
    const auto *hh_439 = buffer.data(hh + 439);
    const auto *hh_440 = buffer.data(hh + 440);

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

#pragma omp simd aligned(t_488, t_489, t_490, t_491, t_492, hh_362, hh_363, hh_364, hh_365, \
                         hh_366, kh_635, kh_636, kh_637, kh_638, \
                         kh_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = -4.0 * hh_362[k]
                   + f_0 * kh_635[k];

        t_489[k] = -4.0 * hh_363[k]
                   + f_0 * kh_636[k];

        t_490[k] = -4.0 * hh_364[k]
                   + f_0 * kh_637[k];

        t_491[k] = -4.0 * hh_365[k]
                   + f_0 * kh_638[k];

        t_492[k] = -4.0 * hh_366[k]
                   + f_0 * kh_639[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, t_496, t_497, hh_367, hh_368, hh_369, hh_370, \
                         hh_371, kh_640, kh_641, kh_642, kh_643, \
                         kh_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = -4.0 * hh_367[k]
                   + f_0 * kh_640[k];

        t_494[k] = -4.0 * hh_368[k]
                   + f_0 * kh_641[k];

        t_495[k] = -4.0 * hh_369[k]
                   + f_0 * kh_642[k];

        t_496[k] = -4.0 * hh_370[k]
                   + f_0 * kh_643[k];

        t_497[k] = -4.0 * hh_371[k]
                   + f_0 * kh_644[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, t_501, t_502, hh_372, hh_373, hh_374, hh_375, \
                         hh_376, kh_645, kh_646, kh_647, kh_648, \
                         kh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = -4.0 * hh_372[k]
                   + f_0 * kh_645[k];

        t_499[k] = -4.0 * hh_373[k]
                   + f_0 * kh_646[k];

        t_500[k] = -4.0 * hh_374[k]
                   + f_0 * kh_647[k];

        t_501[k] = -4.0 * hh_375[k]
                   + f_0 * kh_648[k];

        t_502[k] = -4.0 * hh_376[k]
                   + f_0 * kh_649[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, t_507, hh_377, hh_378, hh_379, hh_380, \
                         hh_381, kh_650, kh_651, kh_652, kh_653, \
                         kh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = -4.0 * hh_377[k]
                   + f_0 * kh_650[k];

        t_504[k] = -3.0 * hh_378[k]
                   + f_0 * kh_651[k];

        t_505[k] = -3.0 * hh_379[k]
                   + f_0 * kh_652[k];

        t_506[k] = -3.0 * hh_380[k]
                   + f_0 * kh_653[k];

        t_507[k] = -3.0 * hh_381[k]
                   + f_0 * kh_654[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, t_512, hh_382, hh_383, hh_384, hh_385, \
                         hh_386, kh_655, kh_656, kh_657, kh_658, \
                         kh_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = -3.0 * hh_382[k]
                   + f_0 * kh_655[k];

        t_509[k] = -3.0 * hh_383[k]
                   + f_0 * kh_656[k];

        t_510[k] = -3.0 * hh_384[k]
                   + f_0 * kh_657[k];

        t_511[k] = -3.0 * hh_385[k]
                   + f_0 * kh_658[k];

        t_512[k] = -3.0 * hh_386[k]
                   + f_0 * kh_659[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, t_516, t_517, hh_387, hh_388, hh_389, hh_390, \
                         hh_391, kh_660, kh_661, kh_662, kh_663, \
                         kh_664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = -3.0 * hh_387[k]
                   + f_0 * kh_660[k];

        t_514[k] = -3.0 * hh_388[k]
                   + f_0 * kh_661[k];

        t_515[k] = -3.0 * hh_389[k]
                   + f_0 * kh_662[k];

        t_516[k] = -3.0 * hh_390[k]
                   + f_0 * kh_663[k];

        t_517[k] = -3.0 * hh_391[k]
                   + f_0 * kh_664[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, t_522, hh_392, hh_393, hh_394, hh_395, \
                         hh_396, kh_665, kh_666, kh_667, kh_668, \
                         kh_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = -3.0 * hh_392[k]
                   + f_0 * kh_665[k];

        t_519[k] = -3.0 * hh_393[k]
                   + f_0 * kh_666[k];

        t_520[k] = -3.0 * hh_394[k]
                   + f_0 * kh_667[k];

        t_521[k] = -3.0 * hh_395[k]
                   + f_0 * kh_668[k];

        t_522[k] = -3.0 * hh_396[k]
                   + f_0 * kh_669[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, t_527, hh_397, hh_398, hh_399, hh_400, \
                         hh_401, kh_670, kh_671, kh_672, kh_673, \
                         kh_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = -3.0 * hh_397[k]
                   + f_0 * kh_670[k];

        t_524[k] = -3.0 * hh_398[k]
                   + f_0 * kh_671[k];

        t_525[k] = -2.0 * hh_399[k]
                   + f_0 * kh_672[k];

        t_526[k] = -2.0 * hh_400[k]
                   + f_0 * kh_673[k];

        t_527[k] = -2.0 * hh_401[k]
                   + f_0 * kh_674[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, hh_402, hh_403, hh_404, hh_405, \
                         hh_406, kh_675, kh_676, kh_677, kh_678, \
                         kh_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = -2.0 * hh_402[k]
                   + f_0 * kh_675[k];

        t_529[k] = -2.0 * hh_403[k]
                   + f_0 * kh_676[k];

        t_530[k] = -2.0 * hh_404[k]
                   + f_0 * kh_677[k];

        t_531[k] = -2.0 * hh_405[k]
                   + f_0 * kh_678[k];

        t_532[k] = -2.0 * hh_406[k]
                   + f_0 * kh_679[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, t_537, hh_407, hh_408, hh_409, hh_410, \
                         hh_411, kh_680, kh_681, kh_682, kh_683, \
                         kh_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = -2.0 * hh_407[k]
                   + f_0 * kh_680[k];

        t_534[k] = -2.0 * hh_408[k]
                   + f_0 * kh_681[k];

        t_535[k] = -2.0 * hh_409[k]
                   + f_0 * kh_682[k];

        t_536[k] = -2.0 * hh_410[k]
                   + f_0 * kh_683[k];

        t_537[k] = -2.0 * hh_411[k]
                   + f_0 * kh_684[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, t_542, hh_412, hh_413, hh_414, hh_415, \
                         hh_416, kh_685, kh_686, kh_687, kh_688, \
                         kh_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = -2.0 * hh_412[k]
                   + f_0 * kh_685[k];

        t_539[k] = -2.0 * hh_413[k]
                   + f_0 * kh_686[k];

        t_540[k] = -2.0 * hh_414[k]
                   + f_0 * kh_687[k];

        t_541[k] = -2.0 * hh_415[k]
                   + f_0 * kh_688[k];

        t_542[k] = -2.0 * hh_416[k]
                   + f_0 * kh_689[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, t_547, hh_417, hh_418, hh_419, hh_420, \
                         hh_421, kh_690, kh_691, kh_692, kh_693, \
                         kh_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = -2.0 * hh_417[k]
                   + f_0 * kh_690[k];

        t_544[k] = -2.0 * hh_418[k]
                   + f_0 * kh_691[k];

        t_545[k] = -2.0 * hh_419[k]
                   + f_0 * kh_692[k];

        t_546[k] = -hh_420[k]
                   + f_0 * kh_693[k];

        t_547[k] = -hh_421[k]
                   + f_0 * kh_694[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, hh_422, hh_423, hh_424, hh_425, \
                         hh_426, kh_695, kh_696, kh_697, kh_698, \
                         kh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = -hh_422[k]
                   + f_0 * kh_695[k];

        t_549[k] = -hh_423[k]
                   + f_0 * kh_696[k];

        t_550[k] = -hh_424[k]
                   + f_0 * kh_697[k];

        t_551[k] = -hh_425[k]
                   + f_0 * kh_698[k];

        t_552[k] = -hh_426[k]
                   + f_0 * kh_699[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, t_557, hh_427, hh_428, hh_429, hh_430, \
                         hh_431, kh_700, kh_701, kh_702, kh_703, \
                         kh_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = -hh_427[k]
                   + f_0 * kh_700[k];

        t_554[k] = -hh_428[k]
                   + f_0 * kh_701[k];

        t_555[k] = -hh_429[k]
                   + f_0 * kh_702[k];

        t_556[k] = -hh_430[k]
                   + f_0 * kh_703[k];

        t_557[k] = -hh_431[k]
                   + f_0 * kh_704[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, t_562, hh_432, hh_433, hh_434, hh_435, \
                         hh_436, kh_705, kh_706, kh_707, kh_708, \
                         kh_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = -hh_432[k]
                   + f_0 * kh_705[k];

        t_559[k] = -hh_433[k]
                   + f_0 * kh_706[k];

        t_560[k] = -hh_434[k]
                   + f_0 * kh_707[k];

        t_561[k] = -hh_435[k]
                   + f_0 * kh_708[k];

        t_562[k] = -hh_436[k]
                   + f_0 * kh_709[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, t_567, t_568, hh_437, hh_438, hh_439, \
                         hh_440, kh_710, kh_711, kh_712, kh_713, kh_714, \
                         kh_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = -hh_437[k]
                   + f_0 * kh_710[k];

        t_564[k] = -hh_438[k]
                   + f_0 * kh_711[k];

        t_565[k] = -hh_439[k]
                   + f_0 * kh_712[k];

        t_566[k] = -hh_440[k]
                   + f_0 * kh_713[k];

        t_567[k] = f_0 * kh_714[k];

        t_568[k] = f_0 * kh_715[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, t_573, t_574, t_575, t_576, kh_716, \
                         kh_717, kh_718, kh_719, kh_720, kh_721, kh_722, \
                         kh_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_0 * kh_716[k];

        t_570[k] = f_0 * kh_717[k];

        t_571[k] = f_0 * kh_718[k];

        t_572[k] = f_0 * kh_719[k];

        t_573[k] = f_0 * kh_720[k];

        t_574[k] = f_0 * kh_721[k];

        t_575[k] = f_0 * kh_722[k];

        t_576[k] = f_0 * kh_723[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, t_582, t_583, t_584, kh_724, \
                         kh_725, kh_726, kh_727, kh_728, kh_729, kh_730, \
                         kh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_0 * kh_724[k];

        t_578[k] = f_0 * kh_725[k];

        t_579[k] = f_0 * kh_726[k];

        t_580[k] = f_0 * kh_727[k];

        t_581[k] = f_0 * kh_728[k];

        t_582[k] = f_0 * kh_729[k];

        t_583[k] = f_0 * kh_730[k];

        t_584[k] = f_0 * kh_731[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, kh_732, kh_733, kh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_0 * kh_732[k];

        t_586[k] = f_0 * kh_733[k];

        t_587[k] = f_0 * kh_734[k];
    }
}

auto
compute_prim_geom_10_ih_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t hh, const size_t kh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ih_electron_repulsion_1_piece0(buffer, target, hh, kh, ncols, alpha);

    compute_prim_geom_10_ih_electron_repulsion_1_piece1(buffer, target, hh, kh, ncols, alpha);

    compute_prim_geom_10_ih_electron_repulsion_1_piece2(buffer, target, hh, kh, ncols, alpha);

    compute_prim_geom_10_ih_electron_repulsion_1_piece3(buffer, target, hh, kh, ncols, alpha);
}

static auto
compute_prim_geom_10_ih_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hh, const size_t kh,
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

    const auto *hh_0 = buffer.data(hh + 0);
    const auto *hh_1 = buffer.data(hh + 1);
    const auto *hh_2 = buffer.data(hh + 2);
    const auto *hh_3 = buffer.data(hh + 3);
    const auto *hh_4 = buffer.data(hh + 4);
    const auto *hh_5 = buffer.data(hh + 5);
    const auto *hh_6 = buffer.data(hh + 6);
    const auto *hh_7 = buffer.data(hh + 7);
    const auto *hh_8 = buffer.data(hh + 8);
    const auto *hh_9 = buffer.data(hh + 9);
    const auto *hh_10 = buffer.data(hh + 10);
    const auto *hh_11 = buffer.data(hh + 11);
    const auto *hh_12 = buffer.data(hh + 12);
    const auto *hh_13 = buffer.data(hh + 13);
    const auto *hh_14 = buffer.data(hh + 14);
    const auto *hh_15 = buffer.data(hh + 15);
    const auto *hh_16 = buffer.data(hh + 16);
    const auto *hh_17 = buffer.data(hh + 17);
    const auto *hh_18 = buffer.data(hh + 18);
    const auto *hh_19 = buffer.data(hh + 19);
    const auto *hh_20 = buffer.data(hh + 20);
    const auto *hh_21 = buffer.data(hh + 21);
    const auto *hh_22 = buffer.data(hh + 22);
    const auto *hh_23 = buffer.data(hh + 23);
    const auto *hh_24 = buffer.data(hh + 24);
    const auto *hh_25 = buffer.data(hh + 25);
    const auto *hh_26 = buffer.data(hh + 26);
    const auto *hh_27 = buffer.data(hh + 27);
    const auto *hh_28 = buffer.data(hh + 28);
    const auto *hh_29 = buffer.data(hh + 29);
    const auto *hh_30 = buffer.data(hh + 30);
    const auto *hh_31 = buffer.data(hh + 31);
    const auto *hh_32 = buffer.data(hh + 32);
    const auto *hh_33 = buffer.data(hh + 33);
    const auto *hh_34 = buffer.data(hh + 34);
    const auto *hh_35 = buffer.data(hh + 35);
    const auto *hh_36 = buffer.data(hh + 36);
    const auto *hh_37 = buffer.data(hh + 37);
    const auto *hh_38 = buffer.data(hh + 38);
    const auto *hh_39 = buffer.data(hh + 39);
    const auto *hh_40 = buffer.data(hh + 40);
    const auto *hh_41 = buffer.data(hh + 41);
    const auto *hh_42 = buffer.data(hh + 42);
    const auto *hh_43 = buffer.data(hh + 43);
    const auto *hh_44 = buffer.data(hh + 44);
    const auto *hh_45 = buffer.data(hh + 45);
    const auto *hh_46 = buffer.data(hh + 46);
    const auto *hh_47 = buffer.data(hh + 47);
    const auto *hh_48 = buffer.data(hh + 48);
    const auto *hh_49 = buffer.data(hh + 49);
    const auto *hh_50 = buffer.data(hh + 50);
    const auto *hh_51 = buffer.data(hh + 51);
    const auto *hh_52 = buffer.data(hh + 52);
    const auto *hh_53 = buffer.data(hh + 53);
    const auto *hh_54 = buffer.data(hh + 54);
    const auto *hh_55 = buffer.data(hh + 55);
    const auto *hh_56 = buffer.data(hh + 56);
    const auto *hh_57 = buffer.data(hh + 57);
    const auto *hh_58 = buffer.data(hh + 58);
    const auto *hh_59 = buffer.data(hh + 59);
    const auto *hh_60 = buffer.data(hh + 60);
    const auto *hh_61 = buffer.data(hh + 61);
    const auto *hh_62 = buffer.data(hh + 62);
    const auto *hh_63 = buffer.data(hh + 63);
    const auto *hh_64 = buffer.data(hh + 64);
    const auto *hh_65 = buffer.data(hh + 65);
    const auto *hh_66 = buffer.data(hh + 66);
    const auto *hh_67 = buffer.data(hh + 67);
    const auto *hh_68 = buffer.data(hh + 68);
    const auto *hh_69 = buffer.data(hh + 69);
    const auto *hh_70 = buffer.data(hh + 70);
    const auto *hh_71 = buffer.data(hh + 71);
    const auto *hh_72 = buffer.data(hh + 72);
    const auto *hh_73 = buffer.data(hh + 73);
    const auto *hh_74 = buffer.data(hh + 74);
    const auto *hh_75 = buffer.data(hh + 75);
    const auto *hh_76 = buffer.data(hh + 76);
    const auto *hh_77 = buffer.data(hh + 77);
    const auto *hh_78 = buffer.data(hh + 78);
    const auto *hh_79 = buffer.data(hh + 79);
    const auto *hh_80 = buffer.data(hh + 80);
    const auto *hh_81 = buffer.data(hh + 81);
    const auto *hh_82 = buffer.data(hh + 82);
    const auto *hh_83 = buffer.data(hh + 83);
    const auto *hh_84 = buffer.data(hh + 84);
    const auto *hh_85 = buffer.data(hh + 85);
    const auto *hh_86 = buffer.data(hh + 86);
    const auto *hh_87 = buffer.data(hh + 87);
    const auto *hh_88 = buffer.data(hh + 88);
    const auto *hh_89 = buffer.data(hh + 89);
    const auto *hh_90 = buffer.data(hh + 90);
    const auto *hh_91 = buffer.data(hh + 91);
    const auto *hh_92 = buffer.data(hh + 92);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, kh_42, kh_43, kh_44, kh_45, \
                         kh_46, kh_47, kh_48, kh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kh_42[k];

        t_1[k] = f_0 * kh_43[k];

        t_2[k] = f_0 * kh_44[k];

        t_3[k] = f_0 * kh_45[k];

        t_4[k] = f_0 * kh_46[k];

        t_5[k] = f_0 * kh_47[k];

        t_6[k] = f_0 * kh_48[k];

        t_7[k] = f_0 * kh_49[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, kh_50, kh_51, kh_52, \
                         kh_53, kh_54, kh_55, kh_56, kh_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * kh_50[k];

        t_9[k] = f_0 * kh_51[k];

        t_10[k] = f_0 * kh_52[k];

        t_11[k] = f_0 * kh_53[k];

        t_12[k] = f_0 * kh_54[k];

        t_13[k] = f_0 * kh_55[k];

        t_14[k] = f_0 * kh_56[k];

        t_15[k] = f_0 * kh_57[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, kh_58, kh_59, kh_60, \
                         kh_61, kh_62, kh_84, kh_85, kh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * kh_58[k];

        t_17[k] = f_0 * kh_59[k];

        t_18[k] = f_0 * kh_60[k];

        t_19[k] = f_0 * kh_61[k];

        t_20[k] = f_0 * kh_62[k];

        t_21[k] = f_0 * kh_84[k];

        t_22[k] = f_0 * kh_85[k];

        t_23[k] = f_0 * kh_86[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, kh_87, kh_88, kh_89, \
                         kh_90, kh_91, kh_92, kh_93, kh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * kh_87[k];

        t_25[k] = f_0 * kh_88[k];

        t_26[k] = f_0 * kh_89[k];

        t_27[k] = f_0 * kh_90[k];

        t_28[k] = f_0 * kh_91[k];

        t_29[k] = f_0 * kh_92[k];

        t_30[k] = f_0 * kh_93[k];

        t_31[k] = f_0 * kh_94[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, kh_95, kh_96, kh_97, \
                         kh_98, kh_99, kh_100, kh_101, kh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * kh_95[k];

        t_33[k] = f_0 * kh_96[k];

        t_34[k] = f_0 * kh_97[k];

        t_35[k] = f_0 * kh_98[k];

        t_36[k] = f_0 * kh_99[k];

        t_37[k] = f_0 * kh_100[k];

        t_38[k] = f_0 * kh_101[k];

        t_39[k] = f_0 * kh_102[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, hh_0, hh_1, hh_2, hh_3, kh_103, \
                         kh_104, kh_105, kh_106, kh_107, kh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * kh_103[k];

        t_41[k] = f_0 * kh_104[k];

        t_42[k] = -hh_0[k]
                  + f_0 * kh_105[k];

        t_43[k] = -hh_1[k]
                  + f_0 * kh_106[k];

        t_44[k] = -hh_2[k]
                  + f_0 * kh_107[k];

        t_45[k] = -hh_3[k]
                  + f_0 * kh_108[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, hh_4, hh_5, hh_6, hh_7, hh_8, kh_109, \
                         kh_110, kh_111, kh_112, kh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -hh_4[k]
                  + f_0 * kh_109[k];

        t_47[k] = -hh_5[k]
                  + f_0 * kh_110[k];

        t_48[k] = -hh_6[k]
                  + f_0 * kh_111[k];

        t_49[k] = -hh_7[k]
                  + f_0 * kh_112[k];

        t_50[k] = -hh_8[k]
                  + f_0 * kh_113[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, hh_9, hh_10, hh_11, hh_12, hh_13, \
                         kh_114, kh_115, kh_116, kh_117, kh_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -hh_9[k]
                  + f_0 * kh_114[k];

        t_52[k] = -hh_10[k]
                  + f_0 * kh_115[k];

        t_53[k] = -hh_11[k]
                  + f_0 * kh_116[k];

        t_54[k] = -hh_12[k]
                  + f_0 * kh_117[k];

        t_55[k] = -hh_13[k]
                  + f_0 * kh_118[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, hh_14, hh_15, hh_16, hh_17, hh_18, \
                         kh_119, kh_120, kh_121, kh_122, kh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -hh_14[k]
                  + f_0 * kh_119[k];

        t_57[k] = -hh_15[k]
                  + f_0 * kh_120[k];

        t_58[k] = -hh_16[k]
                  + f_0 * kh_121[k];

        t_59[k] = -hh_17[k]
                  + f_0 * kh_122[k];

        t_60[k] = -hh_18[k]
                  + f_0 * kh_123[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, t_67, hh_19, hh_20, kh_124, \
                         kh_125, kh_147, kh_148, kh_149, kh_150, \
                         kh_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -hh_19[k]
                  + f_0 * kh_124[k];

        t_62[k] = -hh_20[k]
                  + f_0 * kh_125[k];

        t_63[k] = f_0 * kh_147[k];

        t_64[k] = f_0 * kh_148[k];

        t_65[k] = f_0 * kh_149[k];

        t_66[k] = f_0 * kh_150[k];

        t_67[k] = f_0 * kh_151[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, t_74, t_75, kh_152, kh_153, \
                         kh_154, kh_155, kh_156, kh_157, kh_158, \
                         kh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * kh_152[k];

        t_69[k] = f_0 * kh_153[k];

        t_70[k] = f_0 * kh_154[k];

        t_71[k] = f_0 * kh_155[k];

        t_72[k] = f_0 * kh_156[k];

        t_73[k] = f_0 * kh_157[k];

        t_74[k] = f_0 * kh_158[k];

        t_75[k] = f_0 * kh_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, t_82, t_83, kh_160, kh_161, \
                         kh_162, kh_163, kh_164, kh_165, kh_166, \
                         kh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_0 * kh_160[k];

        t_77[k] = f_0 * kh_161[k];

        t_78[k] = f_0 * kh_162[k];

        t_79[k] = f_0 * kh_163[k];

        t_80[k] = f_0 * kh_164[k];

        t_81[k] = f_0 * kh_165[k];

        t_82[k] = f_0 * kh_166[k];

        t_83[k] = f_0 * kh_167[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, hh_21, hh_22, hh_23, hh_24, hh_25, \
                         kh_168, kh_169, kh_170, kh_171, kh_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -hh_21[k]
                  + f_0 * kh_168[k];

        t_85[k] = -hh_22[k]
                  + f_0 * kh_169[k];

        t_86[k] = -hh_23[k]
                  + f_0 * kh_170[k];

        t_87[k] = -hh_24[k]
                  + f_0 * kh_171[k];

        t_88[k] = -hh_25[k]
                  + f_0 * kh_172[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, hh_26, hh_27, hh_28, hh_29, hh_30, \
                         kh_173, kh_174, kh_175, kh_176, kh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -hh_26[k]
                  + f_0 * kh_173[k];

        t_90[k] = -hh_27[k]
                  + f_0 * kh_174[k];

        t_91[k] = -hh_28[k]
                  + f_0 * kh_175[k];

        t_92[k] = -hh_29[k]
                  + f_0 * kh_176[k];

        t_93[k] = -hh_30[k]
                  + f_0 * kh_177[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, hh_31, hh_32, hh_33, hh_34, hh_35, \
                         kh_178, kh_179, kh_180, kh_181, kh_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -hh_31[k]
                  + f_0 * kh_178[k];

        t_95[k] = -hh_32[k]
                  + f_0 * kh_179[k];

        t_96[k] = -hh_33[k]
                  + f_0 * kh_180[k];

        t_97[k] = -hh_34[k]
                  + f_0 * kh_181[k];

        t_98[k] = -hh_35[k]
                  + f_0 * kh_182[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, hh_36, hh_37, hh_38, hh_39, hh_40, \
                         kh_183, kh_184, kh_185, kh_186, kh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -hh_36[k]
                  + f_0 * kh_183[k];

        t_100[k] = -hh_37[k]
                   + f_0 * kh_184[k];

        t_101[k] = -hh_38[k]
                   + f_0 * kh_185[k];

        t_102[k] = -hh_39[k]
                   + f_0 * kh_186[k];

        t_103[k] = -hh_40[k]
                   + f_0 * kh_187[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, hh_41, hh_42, hh_43, hh_44, hh_45, \
                         kh_188, kh_189, kh_190, kh_191, kh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -hh_41[k]
                   + f_0 * kh_188[k];

        t_105[k] = -2.0 * hh_42[k]
                   + f_0 * kh_189[k];

        t_106[k] = -2.0 * hh_43[k]
                   + f_0 * kh_190[k];

        t_107[k] = -2.0 * hh_44[k]
                   + f_0 * kh_191[k];

        t_108[k] = -2.0 * hh_45[k]
                   + f_0 * kh_192[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, hh_46, hh_47, hh_48, hh_49, hh_50, \
                         kh_193, kh_194, kh_195, kh_196, kh_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -2.0 * hh_46[k]
                   + f_0 * kh_193[k];

        t_110[k] = -2.0 * hh_47[k]
                   + f_0 * kh_194[k];

        t_111[k] = -2.0 * hh_48[k]
                   + f_0 * kh_195[k];

        t_112[k] = -2.0 * hh_49[k]
                   + f_0 * kh_196[k];

        t_113[k] = -2.0 * hh_50[k]
                   + f_0 * kh_197[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, hh_51, hh_52, hh_53, hh_54, hh_55, \
                         kh_198, kh_199, kh_200, kh_201, kh_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -2.0 * hh_51[k]
                   + f_0 * kh_198[k];

        t_115[k] = -2.0 * hh_52[k]
                   + f_0 * kh_199[k];

        t_116[k] = -2.0 * hh_53[k]
                   + f_0 * kh_200[k];

        t_117[k] = -2.0 * hh_54[k]
                   + f_0 * kh_201[k];

        t_118[k] = -2.0 * hh_55[k]
                   + f_0 * kh_202[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, hh_56, hh_57, hh_58, hh_59, hh_60, \
                         kh_203, kh_204, kh_205, kh_206, kh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -2.0 * hh_56[k]
                   + f_0 * kh_203[k];

        t_120[k] = -2.0 * hh_57[k]
                   + f_0 * kh_204[k];

        t_121[k] = -2.0 * hh_58[k]
                   + f_0 * kh_205[k];

        t_122[k] = -2.0 * hh_59[k]
                   + f_0 * kh_206[k];

        t_123[k] = -2.0 * hh_60[k]
                   + f_0 * kh_207[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, t_129, t_130, hh_61, hh_62, \
                         kh_208, kh_209, kh_231, kh_232, kh_233, kh_234, \
                         kh_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * hh_61[k]
                   + f_0 * kh_208[k];

        t_125[k] = -2.0 * hh_62[k]
                   + f_0 * kh_209[k];

        t_126[k] = f_0 * kh_231[k];

        t_127[k] = f_0 * kh_232[k];

        t_128[k] = f_0 * kh_233[k];

        t_129[k] = f_0 * kh_234[k];

        t_130[k] = f_0 * kh_235[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, t_136, t_137, t_138, kh_236, \
                         kh_237, kh_238, kh_239, kh_240, kh_241, kh_242, \
                         kh_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_0 * kh_236[k];

        t_132[k] = f_0 * kh_237[k];

        t_133[k] = f_0 * kh_238[k];

        t_134[k] = f_0 * kh_239[k];

        t_135[k] = f_0 * kh_240[k];

        t_136[k] = f_0 * kh_241[k];

        t_137[k] = f_0 * kh_242[k];

        t_138[k] = f_0 * kh_243[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, t_145, t_146, kh_244, \
                         kh_245, kh_246, kh_247, kh_248, kh_249, kh_250, \
                         kh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_0 * kh_244[k];

        t_140[k] = f_0 * kh_245[k];

        t_141[k] = f_0 * kh_246[k];

        t_142[k] = f_0 * kh_247[k];

        t_143[k] = f_0 * kh_248[k];

        t_144[k] = f_0 * kh_249[k];

        t_145[k] = f_0 * kh_250[k];

        t_146[k] = f_0 * kh_251[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, hh_63, hh_64, hh_65, hh_66, hh_67, \
                         kh_252, kh_253, kh_254, kh_255, kh_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -hh_63[k]
                   + f_0 * kh_252[k];

        t_148[k] = -hh_64[k]
                   + f_0 * kh_253[k];

        t_149[k] = -hh_65[k]
                   + f_0 * kh_254[k];

        t_150[k] = -hh_66[k]
                   + f_0 * kh_255[k];

        t_151[k] = -hh_67[k]
                   + f_0 * kh_256[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, hh_68, hh_69, hh_70, hh_71, hh_72, \
                         kh_257, kh_258, kh_259, kh_260, kh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -hh_68[k]
                   + f_0 * kh_257[k];

        t_153[k] = -hh_69[k]
                   + f_0 * kh_258[k];

        t_154[k] = -hh_70[k]
                   + f_0 * kh_259[k];

        t_155[k] = -hh_71[k]
                   + f_0 * kh_260[k];

        t_156[k] = -hh_72[k]
                   + f_0 * kh_261[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, hh_73, hh_74, hh_75, hh_76, hh_77, \
                         kh_262, kh_263, kh_264, kh_265, kh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -hh_73[k]
                   + f_0 * kh_262[k];

        t_158[k] = -hh_74[k]
                   + f_0 * kh_263[k];

        t_159[k] = -hh_75[k]
                   + f_0 * kh_264[k];

        t_160[k] = -hh_76[k]
                   + f_0 * kh_265[k];

        t_161[k] = -hh_77[k]
                   + f_0 * kh_266[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, hh_78, hh_79, hh_80, hh_81, hh_82, \
                         kh_267, kh_268, kh_269, kh_270, kh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -hh_78[k]
                   + f_0 * kh_267[k];

        t_163[k] = -hh_79[k]
                   + f_0 * kh_268[k];

        t_164[k] = -hh_80[k]
                   + f_0 * kh_269[k];

        t_165[k] = -hh_81[k]
                   + f_0 * kh_270[k];

        t_166[k] = -hh_82[k]
                   + f_0 * kh_271[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, hh_83, hh_84, hh_85, hh_86, hh_87, \
                         kh_272, kh_273, kh_274, kh_275, kh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -hh_83[k]
                   + f_0 * kh_272[k];

        t_168[k] = -2.0 * hh_84[k]
                   + f_0 * kh_273[k];

        t_169[k] = -2.0 * hh_85[k]
                   + f_0 * kh_274[k];

        t_170[k] = -2.0 * hh_86[k]
                   + f_0 * kh_275[k];

        t_171[k] = -2.0 * hh_87[k]
                   + f_0 * kh_276[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, hh_88, hh_89, hh_90, hh_91, hh_92, \
                         kh_277, kh_278, kh_279, kh_280, kh_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -2.0 * hh_88[k]
                   + f_0 * kh_277[k];

        t_173[k] = -2.0 * hh_89[k]
                   + f_0 * kh_278[k];

        t_174[k] = -2.0 * hh_90[k]
                   + f_0 * kh_279[k];

        t_175[k] = -2.0 * hh_91[k]
                   + f_0 * kh_280[k];

        t_176[k] = -2.0 * hh_92[k]
                   + f_0 * kh_281[k];
    }
}

static auto
compute_prim_geom_10_ih_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hh, const size_t kh,
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

    const auto *hh_93 = buffer.data(hh + 93);
    const auto *hh_94 = buffer.data(hh + 94);
    const auto *hh_95 = buffer.data(hh + 95);
    const auto *hh_96 = buffer.data(hh + 96);
    const auto *hh_97 = buffer.data(hh + 97);
    const auto *hh_98 = buffer.data(hh + 98);
    const auto *hh_99 = buffer.data(hh + 99);
    const auto *hh_100 = buffer.data(hh + 100);
    const auto *hh_101 = buffer.data(hh + 101);
    const auto *hh_102 = buffer.data(hh + 102);
    const auto *hh_103 = buffer.data(hh + 103);
    const auto *hh_104 = buffer.data(hh + 104);
    const auto *hh_105 = buffer.data(hh + 105);
    const auto *hh_106 = buffer.data(hh + 106);
    const auto *hh_107 = buffer.data(hh + 107);
    const auto *hh_108 = buffer.data(hh + 108);
    const auto *hh_109 = buffer.data(hh + 109);
    const auto *hh_110 = buffer.data(hh + 110);
    const auto *hh_111 = buffer.data(hh + 111);
    const auto *hh_112 = buffer.data(hh + 112);
    const auto *hh_113 = buffer.data(hh + 113);
    const auto *hh_114 = buffer.data(hh + 114);
    const auto *hh_115 = buffer.data(hh + 115);
    const auto *hh_116 = buffer.data(hh + 116);
    const auto *hh_117 = buffer.data(hh + 117);
    const auto *hh_118 = buffer.data(hh + 118);
    const auto *hh_119 = buffer.data(hh + 119);
    const auto *hh_120 = buffer.data(hh + 120);
    const auto *hh_121 = buffer.data(hh + 121);
    const auto *hh_122 = buffer.data(hh + 122);
    const auto *hh_123 = buffer.data(hh + 123);
    const auto *hh_124 = buffer.data(hh + 124);
    const auto *hh_125 = buffer.data(hh + 125);
    const auto *hh_126 = buffer.data(hh + 126);
    const auto *hh_127 = buffer.data(hh + 127);
    const auto *hh_128 = buffer.data(hh + 128);
    const auto *hh_129 = buffer.data(hh + 129);
    const auto *hh_130 = buffer.data(hh + 130);
    const auto *hh_131 = buffer.data(hh + 131);
    const auto *hh_132 = buffer.data(hh + 132);
    const auto *hh_133 = buffer.data(hh + 133);
    const auto *hh_134 = buffer.data(hh + 134);
    const auto *hh_135 = buffer.data(hh + 135);
    const auto *hh_136 = buffer.data(hh + 136);
    const auto *hh_137 = buffer.data(hh + 137);
    const auto *hh_138 = buffer.data(hh + 138);
    const auto *hh_139 = buffer.data(hh + 139);
    const auto *hh_140 = buffer.data(hh + 140);
    const auto *hh_141 = buffer.data(hh + 141);
    const auto *hh_142 = buffer.data(hh + 142);
    const auto *hh_143 = buffer.data(hh + 143);
    const auto *hh_144 = buffer.data(hh + 144);
    const auto *hh_145 = buffer.data(hh + 145);
    const auto *hh_146 = buffer.data(hh + 146);
    const auto *hh_147 = buffer.data(hh + 147);
    const auto *hh_148 = buffer.data(hh + 148);
    const auto *hh_149 = buffer.data(hh + 149);
    const auto *hh_150 = buffer.data(hh + 150);
    const auto *hh_151 = buffer.data(hh + 151);
    const auto *hh_152 = buffer.data(hh + 152);
    const auto *hh_153 = buffer.data(hh + 153);
    const auto *hh_154 = buffer.data(hh + 154);
    const auto *hh_155 = buffer.data(hh + 155);
    const auto *hh_156 = buffer.data(hh + 156);
    const auto *hh_157 = buffer.data(hh + 157);
    const auto *hh_158 = buffer.data(hh + 158);
    const auto *hh_159 = buffer.data(hh + 159);
    const auto *hh_160 = buffer.data(hh + 160);
    const auto *hh_161 = buffer.data(hh + 161);
    const auto *hh_162 = buffer.data(hh + 162);
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_164 = buffer.data(hh + 164);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_166 = buffer.data(hh + 166);
    const auto *hh_167 = buffer.data(hh + 167);
    const auto *hh_168 = buffer.data(hh + 168);
    const auto *hh_169 = buffer.data(hh + 169);
    const auto *hh_170 = buffer.data(hh + 170);
    const auto *hh_171 = buffer.data(hh + 171);
    const auto *hh_172 = buffer.data(hh + 172);
    const auto *hh_173 = buffer.data(hh + 173);
    const auto *hh_174 = buffer.data(hh + 174);
    const auto *hh_175 = buffer.data(hh + 175);
    const auto *hh_176 = buffer.data(hh + 176);
    const auto *hh_177 = buffer.data(hh + 177);
    const auto *hh_178 = buffer.data(hh + 178);
    const auto *hh_179 = buffer.data(hh + 179);
    const auto *hh_180 = buffer.data(hh + 180);
    const auto *hh_181 = buffer.data(hh + 181);
    const auto *hh_182 = buffer.data(hh + 182);
    const auto *hh_183 = buffer.data(hh + 183);
    const auto *hh_184 = buffer.data(hh + 184);
    const auto *hh_185 = buffer.data(hh + 185);
    const auto *hh_186 = buffer.data(hh + 186);
    const auto *hh_187 = buffer.data(hh + 187);
    const auto *hh_188 = buffer.data(hh + 188);
    const auto *hh_189 = buffer.data(hh + 189);
    const auto *hh_190 = buffer.data(hh + 190);
    const auto *hh_191 = buffer.data(hh + 191);
    const auto *hh_192 = buffer.data(hh + 192);
    const auto *hh_193 = buffer.data(hh + 193);
    const auto *hh_194 = buffer.data(hh + 194);
    const auto *hh_195 = buffer.data(hh + 195);
    const auto *hh_196 = buffer.data(hh + 196);
    const auto *hh_197 = buffer.data(hh + 197);
    const auto *hh_198 = buffer.data(hh + 198);
    const auto *hh_199 = buffer.data(hh + 199);
    const auto *hh_200 = buffer.data(hh + 200);
    const auto *hh_201 = buffer.data(hh + 201);
    const auto *hh_202 = buffer.data(hh + 202);
    const auto *hh_203 = buffer.data(hh + 203);
    const auto *hh_204 = buffer.data(hh + 204);
    const auto *hh_205 = buffer.data(hh + 205);
    const auto *hh_206 = buffer.data(hh + 206);
    const auto *hh_207 = buffer.data(hh + 207);
    const auto *hh_208 = buffer.data(hh + 208);
    const auto *hh_209 = buffer.data(hh + 209);
    const auto *hh_210 = buffer.data(hh + 210);
    const auto *hh_211 = buffer.data(hh + 211);

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

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, hh_93, hh_94, hh_95, hh_96, hh_97, \
                         kh_282, kh_283, kh_284, kh_285, kh_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -2.0 * hh_93[k]
                   + f_0 * kh_282[k];

        t_178[k] = -2.0 * hh_94[k]
                   + f_0 * kh_283[k];

        t_179[k] = -2.0 * hh_95[k]
                   + f_0 * kh_284[k];

        t_180[k] = -2.0 * hh_96[k]
                   + f_0 * kh_285[k];

        t_181[k] = -2.0 * hh_97[k]
                   + f_0 * kh_286[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, hh_98, hh_99, hh_100, hh_101, \
                         hh_102, kh_287, kh_288, kh_289, kh_290, \
                         kh_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * hh_98[k]
                   + f_0 * kh_287[k];

        t_183[k] = -2.0 * hh_99[k]
                   + f_0 * kh_288[k];

        t_184[k] = -2.0 * hh_100[k]
                   + f_0 * kh_289[k];

        t_185[k] = -2.0 * hh_101[k]
                   + f_0 * kh_290[k];

        t_186[k] = -2.0 * hh_102[k]
                   + f_0 * kh_291[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, hh_103, hh_104, hh_105, hh_106, \
                         hh_107, kh_292, kh_293, kh_294, kh_295, \
                         kh_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * hh_103[k]
                   + f_0 * kh_292[k];

        t_188[k] = -2.0 * hh_104[k]
                   + f_0 * kh_293[k];

        t_189[k] = -3.0 * hh_105[k]
                   + f_0 * kh_294[k];

        t_190[k] = -3.0 * hh_106[k]
                   + f_0 * kh_295[k];

        t_191[k] = -3.0 * hh_107[k]
                   + f_0 * kh_296[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, hh_108, hh_109, hh_110, hh_111, \
                         hh_112, kh_297, kh_298, kh_299, kh_300, \
                         kh_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -3.0 * hh_108[k]
                   + f_0 * kh_297[k];

        t_193[k] = -3.0 * hh_109[k]
                   + f_0 * kh_298[k];

        t_194[k] = -3.0 * hh_110[k]
                   + f_0 * kh_299[k];

        t_195[k] = -3.0 * hh_111[k]
                   + f_0 * kh_300[k];

        t_196[k] = -3.0 * hh_112[k]
                   + f_0 * kh_301[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, hh_113, hh_114, hh_115, hh_116, \
                         hh_117, kh_302, kh_303, kh_304, kh_305, \
                         kh_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -3.0 * hh_113[k]
                   + f_0 * kh_302[k];

        t_198[k] = -3.0 * hh_114[k]
                   + f_0 * kh_303[k];

        t_199[k] = -3.0 * hh_115[k]
                   + f_0 * kh_304[k];

        t_200[k] = -3.0 * hh_116[k]
                   + f_0 * kh_305[k];

        t_201[k] = -3.0 * hh_117[k]
                   + f_0 * kh_306[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, hh_118, hh_119, hh_120, hh_121, \
                         hh_122, kh_307, kh_308, kh_309, kh_310, \
                         kh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -3.0 * hh_118[k]
                   + f_0 * kh_307[k];

        t_203[k] = -3.0 * hh_119[k]
                   + f_0 * kh_308[k];

        t_204[k] = -3.0 * hh_120[k]
                   + f_0 * kh_309[k];

        t_205[k] = -3.0 * hh_121[k]
                   + f_0 * kh_310[k];

        t_206[k] = -3.0 * hh_122[k]
                   + f_0 * kh_311[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, hh_123, hh_124, hh_125, \
                         kh_312, kh_313, kh_314, kh_336, kh_337, \
                         kh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -3.0 * hh_123[k]
                   + f_0 * kh_312[k];

        t_208[k] = -3.0 * hh_124[k]
                   + f_0 * kh_313[k];

        t_209[k] = -3.0 * hh_125[k]
                   + f_0 * kh_314[k];

        t_210[k] = f_0 * kh_336[k];

        t_211[k] = f_0 * kh_337[k];

        t_212[k] = f_0 * kh_338[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, t_218, t_219, t_220, kh_339, \
                         kh_340, kh_341, kh_342, kh_343, kh_344, kh_345, \
                         kh_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_0 * kh_339[k];

        t_214[k] = f_0 * kh_340[k];

        t_215[k] = f_0 * kh_341[k];

        t_216[k] = f_0 * kh_342[k];

        t_217[k] = f_0 * kh_343[k];

        t_218[k] = f_0 * kh_344[k];

        t_219[k] = f_0 * kh_345[k];

        t_220[k] = f_0 * kh_346[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, t_226, t_227, t_228, kh_347, \
                         kh_348, kh_349, kh_350, kh_351, kh_352, kh_353, \
                         kh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_0 * kh_347[k];

        t_222[k] = f_0 * kh_348[k];

        t_223[k] = f_0 * kh_349[k];

        t_224[k] = f_0 * kh_350[k];

        t_225[k] = f_0 * kh_351[k];

        t_226[k] = f_0 * kh_352[k];

        t_227[k] = f_0 * kh_353[k];

        t_228[k] = f_0 * kh_354[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, t_234, hh_126, hh_127, hh_128, \
                         hh_129, kh_355, kh_356, kh_357, kh_358, kh_359, \
                         kh_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_0 * kh_355[k];

        t_230[k] = f_0 * kh_356[k];

        t_231[k] = -hh_126[k]
                   + f_0 * kh_357[k];

        t_232[k] = -hh_127[k]
                   + f_0 * kh_358[k];

        t_233[k] = -hh_128[k]
                   + f_0 * kh_359[k];

        t_234[k] = -hh_129[k]
                   + f_0 * kh_360[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, hh_130, hh_131, hh_132, hh_133, \
                         hh_134, kh_361, kh_362, kh_363, kh_364, \
                         kh_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -hh_130[k]
                   + f_0 * kh_361[k];

        t_236[k] = -hh_131[k]
                   + f_0 * kh_362[k];

        t_237[k] = -hh_132[k]
                   + f_0 * kh_363[k];

        t_238[k] = -hh_133[k]
                   + f_0 * kh_364[k];

        t_239[k] = -hh_134[k]
                   + f_0 * kh_365[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, hh_135, hh_136, hh_137, hh_138, \
                         hh_139, kh_366, kh_367, kh_368, kh_369, \
                         kh_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -hh_135[k]
                   + f_0 * kh_366[k];

        t_241[k] = -hh_136[k]
                   + f_0 * kh_367[k];

        t_242[k] = -hh_137[k]
                   + f_0 * kh_368[k];

        t_243[k] = -hh_138[k]
                   + f_0 * kh_369[k];

        t_244[k] = -hh_139[k]
                   + f_0 * kh_370[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, hh_140, hh_141, hh_142, hh_143, \
                         hh_144, kh_371, kh_372, kh_373, kh_374, \
                         kh_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -hh_140[k]
                   + f_0 * kh_371[k];

        t_246[k] = -hh_141[k]
                   + f_0 * kh_372[k];

        t_247[k] = -hh_142[k]
                   + f_0 * kh_373[k];

        t_248[k] = -hh_143[k]
                   + f_0 * kh_374[k];

        t_249[k] = -hh_144[k]
                   + f_0 * kh_375[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, hh_145, hh_146, hh_147, hh_148, \
                         hh_149, kh_376, kh_377, kh_378, kh_379, \
                         kh_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -hh_145[k]
                   + f_0 * kh_376[k];

        t_251[k] = -hh_146[k]
                   + f_0 * kh_377[k];

        t_252[k] = -2.0 * hh_147[k]
                   + f_0 * kh_378[k];

        t_253[k] = -2.0 * hh_148[k]
                   + f_0 * kh_379[k];

        t_254[k] = -2.0 * hh_149[k]
                   + f_0 * kh_380[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, hh_150, hh_151, hh_152, hh_153, \
                         hh_154, kh_381, kh_382, kh_383, kh_384, \
                         kh_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * hh_150[k]
                   + f_0 * kh_381[k];

        t_256[k] = -2.0 * hh_151[k]
                   + f_0 * kh_382[k];

        t_257[k] = -2.0 * hh_152[k]
                   + f_0 * kh_383[k];

        t_258[k] = -2.0 * hh_153[k]
                   + f_0 * kh_384[k];

        t_259[k] = -2.0 * hh_154[k]
                   + f_0 * kh_385[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, hh_155, hh_156, hh_157, hh_158, \
                         hh_159, kh_386, kh_387, kh_388, kh_389, \
                         kh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * hh_155[k]
                   + f_0 * kh_386[k];

        t_261[k] = -2.0 * hh_156[k]
                   + f_0 * kh_387[k];

        t_262[k] = -2.0 * hh_157[k]
                   + f_0 * kh_388[k];

        t_263[k] = -2.0 * hh_158[k]
                   + f_0 * kh_389[k];

        t_264[k] = -2.0 * hh_159[k]
                   + f_0 * kh_390[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, hh_160, hh_161, hh_162, hh_163, \
                         hh_164, kh_391, kh_392, kh_393, kh_394, \
                         kh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * hh_160[k]
                   + f_0 * kh_391[k];

        t_266[k] = -2.0 * hh_161[k]
                   + f_0 * kh_392[k];

        t_267[k] = -2.0 * hh_162[k]
                   + f_0 * kh_393[k];

        t_268[k] = -2.0 * hh_163[k]
                   + f_0 * kh_394[k];

        t_269[k] = -2.0 * hh_164[k]
                   + f_0 * kh_395[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, hh_165, hh_166, hh_167, hh_168, \
                         hh_169, kh_396, kh_397, kh_398, kh_399, \
                         kh_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * hh_165[k]
                   + f_0 * kh_396[k];

        t_271[k] = -2.0 * hh_166[k]
                   + f_0 * kh_397[k];

        t_272[k] = -2.0 * hh_167[k]
                   + f_0 * kh_398[k];

        t_273[k] = -3.0 * hh_168[k]
                   + f_0 * kh_399[k];

        t_274[k] = -3.0 * hh_169[k]
                   + f_0 * kh_400[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, hh_170, hh_171, hh_172, hh_173, \
                         hh_174, kh_401, kh_402, kh_403, kh_404, \
                         kh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * hh_170[k]
                   + f_0 * kh_401[k];

        t_276[k] = -3.0 * hh_171[k]
                   + f_0 * kh_402[k];

        t_277[k] = -3.0 * hh_172[k]
                   + f_0 * kh_403[k];

        t_278[k] = -3.0 * hh_173[k]
                   + f_0 * kh_404[k];

        t_279[k] = -3.0 * hh_174[k]
                   + f_0 * kh_405[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, hh_175, hh_176, hh_177, hh_178, \
                         hh_179, kh_406, kh_407, kh_408, kh_409, \
                         kh_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * hh_175[k]
                   + f_0 * kh_406[k];

        t_281[k] = -3.0 * hh_176[k]
                   + f_0 * kh_407[k];

        t_282[k] = -3.0 * hh_177[k]
                   + f_0 * kh_408[k];

        t_283[k] = -3.0 * hh_178[k]
                   + f_0 * kh_409[k];

        t_284[k] = -3.0 * hh_179[k]
                   + f_0 * kh_410[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, hh_180, hh_181, hh_182, hh_183, \
                         hh_184, kh_411, kh_412, kh_413, kh_414, \
                         kh_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * hh_180[k]
                   + f_0 * kh_411[k];

        t_286[k] = -3.0 * hh_181[k]
                   + f_0 * kh_412[k];

        t_287[k] = -3.0 * hh_182[k]
                   + f_0 * kh_413[k];

        t_288[k] = -3.0 * hh_183[k]
                   + f_0 * kh_414[k];

        t_289[k] = -3.0 * hh_184[k]
                   + f_0 * kh_415[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, hh_185, hh_186, hh_187, hh_188, \
                         hh_189, kh_416, kh_417, kh_418, kh_419, \
                         kh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * hh_185[k]
                   + f_0 * kh_416[k];

        t_291[k] = -3.0 * hh_186[k]
                   + f_0 * kh_417[k];

        t_292[k] = -3.0 * hh_187[k]
                   + f_0 * kh_418[k];

        t_293[k] = -3.0 * hh_188[k]
                   + f_0 * kh_419[k];

        t_294[k] = -4.0 * hh_189[k]
                   + f_0 * kh_420[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, hh_190, hh_191, hh_192, hh_193, \
                         hh_194, kh_421, kh_422, kh_423, kh_424, \
                         kh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * hh_190[k]
                   + f_0 * kh_421[k];

        t_296[k] = -4.0 * hh_191[k]
                   + f_0 * kh_422[k];

        t_297[k] = -4.0 * hh_192[k]
                   + f_0 * kh_423[k];

        t_298[k] = -4.0 * hh_193[k]
                   + f_0 * kh_424[k];

        t_299[k] = -4.0 * hh_194[k]
                   + f_0 * kh_425[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, hh_195, hh_196, hh_197, hh_198, \
                         hh_199, kh_426, kh_427, kh_428, kh_429, \
                         kh_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * hh_195[k]
                   + f_0 * kh_426[k];

        t_301[k] = -4.0 * hh_196[k]
                   + f_0 * kh_427[k];

        t_302[k] = -4.0 * hh_197[k]
                   + f_0 * kh_428[k];

        t_303[k] = -4.0 * hh_198[k]
                   + f_0 * kh_429[k];

        t_304[k] = -4.0 * hh_199[k]
                   + f_0 * kh_430[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, hh_200, hh_201, hh_202, hh_203, \
                         hh_204, kh_431, kh_432, kh_433, kh_434, \
                         kh_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * hh_200[k]
                   + f_0 * kh_431[k];

        t_306[k] = -4.0 * hh_201[k]
                   + f_0 * kh_432[k];

        t_307[k] = -4.0 * hh_202[k]
                   + f_0 * kh_433[k];

        t_308[k] = -4.0 * hh_203[k]
                   + f_0 * kh_434[k];

        t_309[k] = -4.0 * hh_204[k]
                   + f_0 * kh_435[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, hh_205, hh_206, hh_207, hh_208, \
                         hh_209, kh_436, kh_437, kh_438, kh_439, \
                         kh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -4.0 * hh_205[k]
                   + f_0 * kh_436[k];

        t_311[k] = -4.0 * hh_206[k]
                   + f_0 * kh_437[k];

        t_312[k] = -4.0 * hh_207[k]
                   + f_0 * kh_438[k];

        t_313[k] = -4.0 * hh_208[k]
                   + f_0 * kh_439[k];

        t_314[k] = -4.0 * hh_209[k]
                   + f_0 * kh_440[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, t_320, t_321, t_322, kh_462, \
                         kh_463, kh_464, kh_465, kh_466, kh_467, kh_468, \
                         kh_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_0 * kh_462[k];

        t_316[k] = f_0 * kh_463[k];

        t_317[k] = f_0 * kh_464[k];

        t_318[k] = f_0 * kh_465[k];

        t_319[k] = f_0 * kh_466[k];

        t_320[k] = f_0 * kh_467[k];

        t_321[k] = f_0 * kh_468[k];

        t_322[k] = f_0 * kh_469[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, t_328, t_329, t_330, kh_470, \
                         kh_471, kh_472, kh_473, kh_474, kh_475, kh_476, \
                         kh_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_0 * kh_470[k];

        t_324[k] = f_0 * kh_471[k];

        t_325[k] = f_0 * kh_472[k];

        t_326[k] = f_0 * kh_473[k];

        t_327[k] = f_0 * kh_474[k];

        t_328[k] = f_0 * kh_475[k];

        t_329[k] = f_0 * kh_476[k];

        t_330[k] = f_0 * kh_477[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, t_335, t_336, t_337, hh_210, hh_211, \
                         kh_478, kh_479, kh_480, kh_481, kh_482, kh_483, \
                         kh_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_0 * kh_478[k];

        t_332[k] = f_0 * kh_479[k];

        t_333[k] = f_0 * kh_480[k];

        t_334[k] = f_0 * kh_481[k];

        t_335[k] = f_0 * kh_482[k];

        t_336[k] = -hh_210[k]
                   + f_0 * kh_483[k];

        t_337[k] = -hh_211[k]
                   + f_0 * kh_484[k];
    }
}

static auto
compute_prim_geom_10_ih_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hh, const size_t kh,
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

    const auto *hh_212 = buffer.data(hh + 212);
    const auto *hh_213 = buffer.data(hh + 213);
    const auto *hh_214 = buffer.data(hh + 214);
    const auto *hh_215 = buffer.data(hh + 215);
    const auto *hh_216 = buffer.data(hh + 216);
    const auto *hh_217 = buffer.data(hh + 217);
    const auto *hh_218 = buffer.data(hh + 218);
    const auto *hh_219 = buffer.data(hh + 219);
    const auto *hh_220 = buffer.data(hh + 220);
    const auto *hh_221 = buffer.data(hh + 221);
    const auto *hh_222 = buffer.data(hh + 222);
    const auto *hh_223 = buffer.data(hh + 223);
    const auto *hh_224 = buffer.data(hh + 224);
    const auto *hh_225 = buffer.data(hh + 225);
    const auto *hh_226 = buffer.data(hh + 226);
    const auto *hh_227 = buffer.data(hh + 227);
    const auto *hh_228 = buffer.data(hh + 228);
    const auto *hh_229 = buffer.data(hh + 229);
    const auto *hh_230 = buffer.data(hh + 230);
    const auto *hh_231 = buffer.data(hh + 231);
    const auto *hh_232 = buffer.data(hh + 232);
    const auto *hh_233 = buffer.data(hh + 233);
    const auto *hh_234 = buffer.data(hh + 234);
    const auto *hh_235 = buffer.data(hh + 235);
    const auto *hh_236 = buffer.data(hh + 236);
    const auto *hh_237 = buffer.data(hh + 237);
    const auto *hh_238 = buffer.data(hh + 238);
    const auto *hh_239 = buffer.data(hh + 239);
    const auto *hh_240 = buffer.data(hh + 240);
    const auto *hh_241 = buffer.data(hh + 241);
    const auto *hh_242 = buffer.data(hh + 242);
    const auto *hh_243 = buffer.data(hh + 243);
    const auto *hh_244 = buffer.data(hh + 244);
    const auto *hh_245 = buffer.data(hh + 245);
    const auto *hh_246 = buffer.data(hh + 246);
    const auto *hh_247 = buffer.data(hh + 247);
    const auto *hh_248 = buffer.data(hh + 248);
    const auto *hh_249 = buffer.data(hh + 249);
    const auto *hh_250 = buffer.data(hh + 250);
    const auto *hh_251 = buffer.data(hh + 251);
    const auto *hh_252 = buffer.data(hh + 252);
    const auto *hh_253 = buffer.data(hh + 253);
    const auto *hh_254 = buffer.data(hh + 254);
    const auto *hh_255 = buffer.data(hh + 255);
    const auto *hh_256 = buffer.data(hh + 256);
    const auto *hh_257 = buffer.data(hh + 257);
    const auto *hh_258 = buffer.data(hh + 258);
    const auto *hh_259 = buffer.data(hh + 259);
    const auto *hh_260 = buffer.data(hh + 260);
    const auto *hh_261 = buffer.data(hh + 261);
    const auto *hh_262 = buffer.data(hh + 262);
    const auto *hh_263 = buffer.data(hh + 263);
    const auto *hh_264 = buffer.data(hh + 264);
    const auto *hh_265 = buffer.data(hh + 265);
    const auto *hh_266 = buffer.data(hh + 266);
    const auto *hh_267 = buffer.data(hh + 267);
    const auto *hh_268 = buffer.data(hh + 268);
    const auto *hh_269 = buffer.data(hh + 269);
    const auto *hh_270 = buffer.data(hh + 270);
    const auto *hh_271 = buffer.data(hh + 271);
    const auto *hh_272 = buffer.data(hh + 272);
    const auto *hh_273 = buffer.data(hh + 273);
    const auto *hh_274 = buffer.data(hh + 274);
    const auto *hh_275 = buffer.data(hh + 275);
    const auto *hh_276 = buffer.data(hh + 276);
    const auto *hh_277 = buffer.data(hh + 277);
    const auto *hh_278 = buffer.data(hh + 278);
    const auto *hh_279 = buffer.data(hh + 279);
    const auto *hh_280 = buffer.data(hh + 280);
    const auto *hh_281 = buffer.data(hh + 281);
    const auto *hh_282 = buffer.data(hh + 282);
    const auto *hh_283 = buffer.data(hh + 283);
    const auto *hh_284 = buffer.data(hh + 284);
    const auto *hh_285 = buffer.data(hh + 285);
    const auto *hh_286 = buffer.data(hh + 286);
    const auto *hh_287 = buffer.data(hh + 287);
    const auto *hh_288 = buffer.data(hh + 288);
    const auto *hh_289 = buffer.data(hh + 289);
    const auto *hh_290 = buffer.data(hh + 290);
    const auto *hh_291 = buffer.data(hh + 291);
    const auto *hh_292 = buffer.data(hh + 292);
    const auto *hh_293 = buffer.data(hh + 293);
    const auto *hh_294 = buffer.data(hh + 294);
    const auto *hh_295 = buffer.data(hh + 295);
    const auto *hh_296 = buffer.data(hh + 296);
    const auto *hh_297 = buffer.data(hh + 297);
    const auto *hh_298 = buffer.data(hh + 298);
    const auto *hh_299 = buffer.data(hh + 299);
    const auto *hh_300 = buffer.data(hh + 300);
    const auto *hh_301 = buffer.data(hh + 301);
    const auto *hh_302 = buffer.data(hh + 302);
    const auto *hh_303 = buffer.data(hh + 303);
    const auto *hh_304 = buffer.data(hh + 304);
    const auto *hh_305 = buffer.data(hh + 305);
    const auto *hh_306 = buffer.data(hh + 306);
    const auto *hh_307 = buffer.data(hh + 307);
    const auto *hh_308 = buffer.data(hh + 308);
    const auto *hh_309 = buffer.data(hh + 309);
    const auto *hh_310 = buffer.data(hh + 310);
    const auto *hh_311 = buffer.data(hh + 311);
    const auto *hh_312 = buffer.data(hh + 312);
    const auto *hh_313 = buffer.data(hh + 313);
    const auto *hh_314 = buffer.data(hh + 314);
    const auto *hh_315 = buffer.data(hh + 315);
    const auto *hh_316 = buffer.data(hh + 316);
    const auto *hh_317 = buffer.data(hh + 317);
    const auto *hh_318 = buffer.data(hh + 318);
    const auto *hh_319 = buffer.data(hh + 319);
    const auto *hh_320 = buffer.data(hh + 320);
    const auto *hh_321 = buffer.data(hh + 321);
    const auto *hh_322 = buffer.data(hh + 322);
    const auto *hh_323 = buffer.data(hh + 323);
    const auto *hh_324 = buffer.data(hh + 324);
    const auto *hh_325 = buffer.data(hh + 325);
    const auto *hh_326 = buffer.data(hh + 326);
    const auto *hh_327 = buffer.data(hh + 327);
    const auto *hh_328 = buffer.data(hh + 328);
    const auto *hh_329 = buffer.data(hh + 329);
    const auto *hh_330 = buffer.data(hh + 330);
    const auto *hh_331 = buffer.data(hh + 331);
    const auto *hh_332 = buffer.data(hh + 332);
    const auto *hh_333 = buffer.data(hh + 333);
    const auto *hh_334 = buffer.data(hh + 334);
    const auto *hh_335 = buffer.data(hh + 335);
    const auto *hh_336 = buffer.data(hh + 336);
    const auto *hh_337 = buffer.data(hh + 337);
    const auto *hh_338 = buffer.data(hh + 338);
    const auto *hh_339 = buffer.data(hh + 339);
    const auto *hh_340 = buffer.data(hh + 340);
    const auto *hh_341 = buffer.data(hh + 341);
    const auto *hh_342 = buffer.data(hh + 342);
    const auto *hh_343 = buffer.data(hh + 343);
    const auto *hh_344 = buffer.data(hh + 344);
    const auto *hh_345 = buffer.data(hh + 345);
    const auto *hh_346 = buffer.data(hh + 346);
    const auto *hh_347 = buffer.data(hh + 347);
    const auto *hh_348 = buffer.data(hh + 348);

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

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, hh_212, hh_213, hh_214, hh_215, \
                         hh_216, kh_485, kh_486, kh_487, kh_488, \
                         kh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = -hh_212[k]
                   + f_0 * kh_485[k];

        t_339[k] = -hh_213[k]
                   + f_0 * kh_486[k];

        t_340[k] = -hh_214[k]
                   + f_0 * kh_487[k];

        t_341[k] = -hh_215[k]
                   + f_0 * kh_488[k];

        t_342[k] = -hh_216[k]
                   + f_0 * kh_489[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, hh_217, hh_218, hh_219, hh_220, \
                         hh_221, kh_490, kh_491, kh_492, kh_493, \
                         kh_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = -hh_217[k]
                   + f_0 * kh_490[k];

        t_344[k] = -hh_218[k]
                   + f_0 * kh_491[k];

        t_345[k] = -hh_219[k]
                   + f_0 * kh_492[k];

        t_346[k] = -hh_220[k]
                   + f_0 * kh_493[k];

        t_347[k] = -hh_221[k]
                   + f_0 * kh_494[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, hh_222, hh_223, hh_224, hh_225, \
                         hh_226, kh_495, kh_496, kh_497, kh_498, \
                         kh_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = -hh_222[k]
                   + f_0 * kh_495[k];

        t_349[k] = -hh_223[k]
                   + f_0 * kh_496[k];

        t_350[k] = -hh_224[k]
                   + f_0 * kh_497[k];

        t_351[k] = -hh_225[k]
                   + f_0 * kh_498[k];

        t_352[k] = -hh_226[k]
                   + f_0 * kh_499[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, hh_227, hh_228, hh_229, hh_230, \
                         hh_231, kh_500, kh_501, kh_502, kh_503, \
                         kh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = -hh_227[k]
                   + f_0 * kh_500[k];

        t_354[k] = -hh_228[k]
                   + f_0 * kh_501[k];

        t_355[k] = -hh_229[k]
                   + f_0 * kh_502[k];

        t_356[k] = -hh_230[k]
                   + f_0 * kh_503[k];

        t_357[k] = -2.0 * hh_231[k]
                   + f_0 * kh_504[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, hh_232, hh_233, hh_234, hh_235, \
                         hh_236, kh_505, kh_506, kh_507, kh_508, \
                         kh_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = -2.0 * hh_232[k]
                   + f_0 * kh_505[k];

        t_359[k] = -2.0 * hh_233[k]
                   + f_0 * kh_506[k];

        t_360[k] = -2.0 * hh_234[k]
                   + f_0 * kh_507[k];

        t_361[k] = -2.0 * hh_235[k]
                   + f_0 * kh_508[k];

        t_362[k] = -2.0 * hh_236[k]
                   + f_0 * kh_509[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, hh_237, hh_238, hh_239, hh_240, \
                         hh_241, kh_510, kh_511, kh_512, kh_513, \
                         kh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = -2.0 * hh_237[k]
                   + f_0 * kh_510[k];

        t_364[k] = -2.0 * hh_238[k]
                   + f_0 * kh_511[k];

        t_365[k] = -2.0 * hh_239[k]
                   + f_0 * kh_512[k];

        t_366[k] = -2.0 * hh_240[k]
                   + f_0 * kh_513[k];

        t_367[k] = -2.0 * hh_241[k]
                   + f_0 * kh_514[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, hh_242, hh_243, hh_244, hh_245, \
                         hh_246, kh_515, kh_516, kh_517, kh_518, \
                         kh_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = -2.0 * hh_242[k]
                   + f_0 * kh_515[k];

        t_369[k] = -2.0 * hh_243[k]
                   + f_0 * kh_516[k];

        t_370[k] = -2.0 * hh_244[k]
                   + f_0 * kh_517[k];

        t_371[k] = -2.0 * hh_245[k]
                   + f_0 * kh_518[k];

        t_372[k] = -2.0 * hh_246[k]
                   + f_0 * kh_519[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, hh_247, hh_248, hh_249, hh_250, \
                         hh_251, kh_520, kh_521, kh_522, kh_523, \
                         kh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = -2.0 * hh_247[k]
                   + f_0 * kh_520[k];

        t_374[k] = -2.0 * hh_248[k]
                   + f_0 * kh_521[k];

        t_375[k] = -2.0 * hh_249[k]
                   + f_0 * kh_522[k];

        t_376[k] = -2.0 * hh_250[k]
                   + f_0 * kh_523[k];

        t_377[k] = -2.0 * hh_251[k]
                   + f_0 * kh_524[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, hh_252, hh_253, hh_254, hh_255, \
                         hh_256, kh_525, kh_526, kh_527, kh_528, \
                         kh_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = -3.0 * hh_252[k]
                   + f_0 * kh_525[k];

        t_379[k] = -3.0 * hh_253[k]
                   + f_0 * kh_526[k];

        t_380[k] = -3.0 * hh_254[k]
                   + f_0 * kh_527[k];

        t_381[k] = -3.0 * hh_255[k]
                   + f_0 * kh_528[k];

        t_382[k] = -3.0 * hh_256[k]
                   + f_0 * kh_529[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, hh_257, hh_258, hh_259, hh_260, \
                         hh_261, kh_530, kh_531, kh_532, kh_533, \
                         kh_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = -3.0 * hh_257[k]
                   + f_0 * kh_530[k];

        t_384[k] = -3.0 * hh_258[k]
                   + f_0 * kh_531[k];

        t_385[k] = -3.0 * hh_259[k]
                   + f_0 * kh_532[k];

        t_386[k] = -3.0 * hh_260[k]
                   + f_0 * kh_533[k];

        t_387[k] = -3.0 * hh_261[k]
                   + f_0 * kh_534[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, hh_262, hh_263, hh_264, hh_265, \
                         hh_266, kh_535, kh_536, kh_537, kh_538, \
                         kh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = -3.0 * hh_262[k]
                   + f_0 * kh_535[k];

        t_389[k] = -3.0 * hh_263[k]
                   + f_0 * kh_536[k];

        t_390[k] = -3.0 * hh_264[k]
                   + f_0 * kh_537[k];

        t_391[k] = -3.0 * hh_265[k]
                   + f_0 * kh_538[k];

        t_392[k] = -3.0 * hh_266[k]
                   + f_0 * kh_539[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, t_397, hh_267, hh_268, hh_269, hh_270, \
                         hh_271, kh_540, kh_541, kh_542, kh_543, \
                         kh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = -3.0 * hh_267[k]
                   + f_0 * kh_540[k];

        t_394[k] = -3.0 * hh_268[k]
                   + f_0 * kh_541[k];

        t_395[k] = -3.0 * hh_269[k]
                   + f_0 * kh_542[k];

        t_396[k] = -3.0 * hh_270[k]
                   + f_0 * kh_543[k];

        t_397[k] = -3.0 * hh_271[k]
                   + f_0 * kh_544[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, t_402, hh_272, hh_273, hh_274, hh_275, \
                         hh_276, kh_545, kh_546, kh_547, kh_548, \
                         kh_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = -3.0 * hh_272[k]
                   + f_0 * kh_545[k];

        t_399[k] = -4.0 * hh_273[k]
                   + f_0 * kh_546[k];

        t_400[k] = -4.0 * hh_274[k]
                   + f_0 * kh_547[k];

        t_401[k] = -4.0 * hh_275[k]
                   + f_0 * kh_548[k];

        t_402[k] = -4.0 * hh_276[k]
                   + f_0 * kh_549[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, t_407, hh_277, hh_278, hh_279, hh_280, \
                         hh_281, kh_550, kh_551, kh_552, kh_553, \
                         kh_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = -4.0 * hh_277[k]
                   + f_0 * kh_550[k];

        t_404[k] = -4.0 * hh_278[k]
                   + f_0 * kh_551[k];

        t_405[k] = -4.0 * hh_279[k]
                   + f_0 * kh_552[k];

        t_406[k] = -4.0 * hh_280[k]
                   + f_0 * kh_553[k];

        t_407[k] = -4.0 * hh_281[k]
                   + f_0 * kh_554[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, hh_282, hh_283, hh_284, hh_285, \
                         hh_286, kh_555, kh_556, kh_557, kh_558, \
                         kh_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = -4.0 * hh_282[k]
                   + f_0 * kh_555[k];

        t_409[k] = -4.0 * hh_283[k]
                   + f_0 * kh_556[k];

        t_410[k] = -4.0 * hh_284[k]
                   + f_0 * kh_557[k];

        t_411[k] = -4.0 * hh_285[k]
                   + f_0 * kh_558[k];

        t_412[k] = -4.0 * hh_286[k]
                   + f_0 * kh_559[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, hh_287, hh_288, hh_289, hh_290, \
                         hh_291, kh_560, kh_561, kh_562, kh_563, \
                         kh_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = -4.0 * hh_287[k]
                   + f_0 * kh_560[k];

        t_414[k] = -4.0 * hh_288[k]
                   + f_0 * kh_561[k];

        t_415[k] = -4.0 * hh_289[k]
                   + f_0 * kh_562[k];

        t_416[k] = -4.0 * hh_290[k]
                   + f_0 * kh_563[k];

        t_417[k] = -4.0 * hh_291[k]
                   + f_0 * kh_564[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, t_422, hh_292, hh_293, hh_294, hh_295, \
                         hh_296, kh_565, kh_566, kh_567, kh_568, \
                         kh_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = -4.0 * hh_292[k]
                   + f_0 * kh_565[k];

        t_419[k] = -4.0 * hh_293[k]
                   + f_0 * kh_566[k];

        t_420[k] = -5.0 * hh_294[k]
                   + f_0 * kh_567[k];

        t_421[k] = -5.0 * hh_295[k]
                   + f_0 * kh_568[k];

        t_422[k] = -5.0 * hh_296[k]
                   + f_0 * kh_569[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, hh_297, hh_298, hh_299, hh_300, \
                         hh_301, kh_570, kh_571, kh_572, kh_573, \
                         kh_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = -5.0 * hh_297[k]
                   + f_0 * kh_570[k];

        t_424[k] = -5.0 * hh_298[k]
                   + f_0 * kh_571[k];

        t_425[k] = -5.0 * hh_299[k]
                   + f_0 * kh_572[k];

        t_426[k] = -5.0 * hh_300[k]
                   + f_0 * kh_573[k];

        t_427[k] = -5.0 * hh_301[k]
                   + f_0 * kh_574[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, hh_302, hh_303, hh_304, hh_305, \
                         hh_306, kh_575, kh_576, kh_577, kh_578, \
                         kh_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = -5.0 * hh_302[k]
                   + f_0 * kh_575[k];

        t_429[k] = -5.0 * hh_303[k]
                   + f_0 * kh_576[k];

        t_430[k] = -5.0 * hh_304[k]
                   + f_0 * kh_577[k];

        t_431[k] = -5.0 * hh_305[k]
                   + f_0 * kh_578[k];

        t_432[k] = -5.0 * hh_306[k]
                   + f_0 * kh_579[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, hh_307, hh_308, hh_309, hh_310, \
                         hh_311, kh_580, kh_581, kh_582, kh_583, \
                         kh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = -5.0 * hh_307[k]
                   + f_0 * kh_580[k];

        t_434[k] = -5.0 * hh_308[k]
                   + f_0 * kh_581[k];

        t_435[k] = -5.0 * hh_309[k]
                   + f_0 * kh_582[k];

        t_436[k] = -5.0 * hh_310[k]
                   + f_0 * kh_583[k];

        t_437[k] = -5.0 * hh_311[k]
                   + f_0 * kh_584[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, t_442, t_443, hh_312, hh_313, hh_314, \
                         kh_585, kh_586, kh_587, kh_609, kh_610, \
                         kh_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = -5.0 * hh_312[k]
                   + f_0 * kh_585[k];

        t_439[k] = -5.0 * hh_313[k]
                   + f_0 * kh_586[k];

        t_440[k] = -5.0 * hh_314[k]
                   + f_0 * kh_587[k];

        t_441[k] = f_0 * kh_609[k];

        t_442[k] = f_0 * kh_610[k];

        t_443[k] = f_0 * kh_611[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, t_449, t_450, t_451, kh_612, \
                         kh_613, kh_614, kh_615, kh_616, kh_617, kh_618, \
                         kh_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_0 * kh_612[k];

        t_445[k] = f_0 * kh_613[k];

        t_446[k] = f_0 * kh_614[k];

        t_447[k] = f_0 * kh_615[k];

        t_448[k] = f_0 * kh_616[k];

        t_449[k] = f_0 * kh_617[k];

        t_450[k] = f_0 * kh_618[k];

        t_451[k] = f_0 * kh_619[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, t_457, t_458, t_459, kh_620, \
                         kh_621, kh_622, kh_623, kh_624, kh_625, kh_626, \
                         kh_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_0 * kh_620[k];

        t_453[k] = f_0 * kh_621[k];

        t_454[k] = f_0 * kh_622[k];

        t_455[k] = f_0 * kh_623[k];

        t_456[k] = f_0 * kh_624[k];

        t_457[k] = f_0 * kh_625[k];

        t_458[k] = f_0 * kh_626[k];

        t_459[k] = f_0 * kh_627[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, t_465, hh_315, hh_316, hh_317, \
                         hh_318, kh_628, kh_629, kh_630, kh_631, kh_632, \
                         kh_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_0 * kh_628[k];

        t_461[k] = f_0 * kh_629[k];

        t_462[k] = -hh_315[k]
                   + f_0 * kh_630[k];

        t_463[k] = -hh_316[k]
                   + f_0 * kh_631[k];

        t_464[k] = -hh_317[k]
                   + f_0 * kh_632[k];

        t_465[k] = -hh_318[k]
                   + f_0 * kh_633[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, hh_319, hh_320, hh_321, hh_322, \
                         hh_323, kh_634, kh_635, kh_636, kh_637, \
                         kh_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = -hh_319[k]
                   + f_0 * kh_634[k];

        t_467[k] = -hh_320[k]
                   + f_0 * kh_635[k];

        t_468[k] = -hh_321[k]
                   + f_0 * kh_636[k];

        t_469[k] = -hh_322[k]
                   + f_0 * kh_637[k];

        t_470[k] = -hh_323[k]
                   + f_0 * kh_638[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, t_475, hh_324, hh_325, hh_326, hh_327, \
                         hh_328, kh_639, kh_640, kh_641, kh_642, \
                         kh_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = -hh_324[k]
                   + f_0 * kh_639[k];

        t_472[k] = -hh_325[k]
                   + f_0 * kh_640[k];

        t_473[k] = -hh_326[k]
                   + f_0 * kh_641[k];

        t_474[k] = -hh_327[k]
                   + f_0 * kh_642[k];

        t_475[k] = -hh_328[k]
                   + f_0 * kh_643[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, hh_329, hh_330, hh_331, hh_332, \
                         hh_333, kh_644, kh_645, kh_646, kh_647, \
                         kh_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = -hh_329[k]
                   + f_0 * kh_644[k];

        t_477[k] = -hh_330[k]
                   + f_0 * kh_645[k];

        t_478[k] = -hh_331[k]
                   + f_0 * kh_646[k];

        t_479[k] = -hh_332[k]
                   + f_0 * kh_647[k];

        t_480[k] = -hh_333[k]
                   + f_0 * kh_648[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, hh_334, hh_335, hh_336, hh_337, \
                         hh_338, kh_649, kh_650, kh_651, kh_652, \
                         kh_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = -hh_334[k]
                   + f_0 * kh_649[k];

        t_482[k] = -hh_335[k]
                   + f_0 * kh_650[k];

        t_483[k] = -2.0 * hh_336[k]
                   + f_0 * kh_651[k];

        t_484[k] = -2.0 * hh_337[k]
                   + f_0 * kh_652[k];

        t_485[k] = -2.0 * hh_338[k]
                   + f_0 * kh_653[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, hh_339, hh_340, hh_341, hh_342, \
                         hh_343, kh_654, kh_655, kh_656, kh_657, \
                         kh_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = -2.0 * hh_339[k]
                   + f_0 * kh_654[k];

        t_487[k] = -2.0 * hh_340[k]
                   + f_0 * kh_655[k];

        t_488[k] = -2.0 * hh_341[k]
                   + f_0 * kh_656[k];

        t_489[k] = -2.0 * hh_342[k]
                   + f_0 * kh_657[k];

        t_490[k] = -2.0 * hh_343[k]
                   + f_0 * kh_658[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, hh_344, hh_345, hh_346, hh_347, \
                         hh_348, kh_659, kh_660, kh_661, kh_662, \
                         kh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = -2.0 * hh_344[k]
                   + f_0 * kh_659[k];

        t_492[k] = -2.0 * hh_345[k]
                   + f_0 * kh_660[k];

        t_493[k] = -2.0 * hh_346[k]
                   + f_0 * kh_661[k];

        t_494[k] = -2.0 * hh_347[k]
                   + f_0 * kh_662[k];

        t_495[k] = -2.0 * hh_348[k]
                   + f_0 * kh_663[k];
    }
}

static auto
compute_prim_geom_10_ih_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hh, const size_t kh,
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

    const auto *hh_349 = buffer.data(hh + 349);
    const auto *hh_350 = buffer.data(hh + 350);
    const auto *hh_351 = buffer.data(hh + 351);
    const auto *hh_352 = buffer.data(hh + 352);
    const auto *hh_353 = buffer.data(hh + 353);
    const auto *hh_354 = buffer.data(hh + 354);
    const auto *hh_355 = buffer.data(hh + 355);
    const auto *hh_356 = buffer.data(hh + 356);
    const auto *hh_357 = buffer.data(hh + 357);
    const auto *hh_358 = buffer.data(hh + 358);
    const auto *hh_359 = buffer.data(hh + 359);
    const auto *hh_360 = buffer.data(hh + 360);
    const auto *hh_361 = buffer.data(hh + 361);
    const auto *hh_362 = buffer.data(hh + 362);
    const auto *hh_363 = buffer.data(hh + 363);
    const auto *hh_364 = buffer.data(hh + 364);
    const auto *hh_365 = buffer.data(hh + 365);
    const auto *hh_366 = buffer.data(hh + 366);
    const auto *hh_367 = buffer.data(hh + 367);
    const auto *hh_368 = buffer.data(hh + 368);
    const auto *hh_369 = buffer.data(hh + 369);
    const auto *hh_370 = buffer.data(hh + 370);
    const auto *hh_371 = buffer.data(hh + 371);
    const auto *hh_372 = buffer.data(hh + 372);
    const auto *hh_373 = buffer.data(hh + 373);
    const auto *hh_374 = buffer.data(hh + 374);
    const auto *hh_375 = buffer.data(hh + 375);
    const auto *hh_376 = buffer.data(hh + 376);
    const auto *hh_377 = buffer.data(hh + 377);
    const auto *hh_378 = buffer.data(hh + 378);
    const auto *hh_379 = buffer.data(hh + 379);
    const auto *hh_380 = buffer.data(hh + 380);
    const auto *hh_381 = buffer.data(hh + 381);
    const auto *hh_382 = buffer.data(hh + 382);
    const auto *hh_383 = buffer.data(hh + 383);
    const auto *hh_384 = buffer.data(hh + 384);
    const auto *hh_385 = buffer.data(hh + 385);
    const auto *hh_386 = buffer.data(hh + 386);
    const auto *hh_387 = buffer.data(hh + 387);
    const auto *hh_388 = buffer.data(hh + 388);
    const auto *hh_389 = buffer.data(hh + 389);
    const auto *hh_390 = buffer.data(hh + 390);
    const auto *hh_391 = buffer.data(hh + 391);
    const auto *hh_392 = buffer.data(hh + 392);
    const auto *hh_393 = buffer.data(hh + 393);
    const auto *hh_394 = buffer.data(hh + 394);
    const auto *hh_395 = buffer.data(hh + 395);
    const auto *hh_396 = buffer.data(hh + 396);
    const auto *hh_397 = buffer.data(hh + 397);
    const auto *hh_398 = buffer.data(hh + 398);
    const auto *hh_399 = buffer.data(hh + 399);
    const auto *hh_400 = buffer.data(hh + 400);
    const auto *hh_401 = buffer.data(hh + 401);
    const auto *hh_402 = buffer.data(hh + 402);
    const auto *hh_403 = buffer.data(hh + 403);
    const auto *hh_404 = buffer.data(hh + 404);
    const auto *hh_405 = buffer.data(hh + 405);
    const auto *hh_406 = buffer.data(hh + 406);
    const auto *hh_407 = buffer.data(hh + 407);
    const auto *hh_408 = buffer.data(hh + 408);
    const auto *hh_409 = buffer.data(hh + 409);
    const auto *hh_410 = buffer.data(hh + 410);
    const auto *hh_411 = buffer.data(hh + 411);
    const auto *hh_412 = buffer.data(hh + 412);
    const auto *hh_413 = buffer.data(hh + 413);
    const auto *hh_414 = buffer.data(hh + 414);
    const auto *hh_415 = buffer.data(hh + 415);
    const auto *hh_416 = buffer.data(hh + 416);
    const auto *hh_417 = buffer.data(hh + 417);
    const auto *hh_418 = buffer.data(hh + 418);
    const auto *hh_419 = buffer.data(hh + 419);
    const auto *hh_420 = buffer.data(hh + 420);
    const auto *hh_421 = buffer.data(hh + 421);
    const auto *hh_422 = buffer.data(hh + 422);
    const auto *hh_423 = buffer.data(hh + 423);
    const auto *hh_424 = buffer.data(hh + 424);
    const auto *hh_425 = buffer.data(hh + 425);
    const auto *hh_426 = buffer.data(hh + 426);
    const auto *hh_427 = buffer.data(hh + 427);
    const auto *hh_428 = buffer.data(hh + 428);
    const auto *hh_429 = buffer.data(hh + 429);
    const auto *hh_430 = buffer.data(hh + 430);
    const auto *hh_431 = buffer.data(hh + 431);
    const auto *hh_432 = buffer.data(hh + 432);
    const auto *hh_433 = buffer.data(hh + 433);
    const auto *hh_434 = buffer.data(hh + 434);
    const auto *hh_435 = buffer.data(hh + 435);
    const auto *hh_436 = buffer.data(hh + 436);
    const auto *hh_437 = buffer.data(hh + 437);
    const auto *hh_438 = buffer.data(hh + 438);
    const auto *hh_439 = buffer.data(hh + 439);
    const auto *hh_440 = buffer.data(hh + 440);

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

#pragma omp simd aligned(t_496, t_497, t_498, t_499, t_500, hh_349, hh_350, hh_351, hh_352, \
                         hh_353, kh_664, kh_665, kh_666, kh_667, \
                         kh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = -2.0 * hh_349[k]
                   + f_0 * kh_664[k];

        t_497[k] = -2.0 * hh_350[k]
                   + f_0 * kh_665[k];

        t_498[k] = -2.0 * hh_351[k]
                   + f_0 * kh_666[k];

        t_499[k] = -2.0 * hh_352[k]
                   + f_0 * kh_667[k];

        t_500[k] = -2.0 * hh_353[k]
                   + f_0 * kh_668[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, t_505, hh_354, hh_355, hh_356, hh_357, \
                         hh_358, kh_669, kh_670, kh_671, kh_672, \
                         kh_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = -2.0 * hh_354[k]
                   + f_0 * kh_669[k];

        t_502[k] = -2.0 * hh_355[k]
                   + f_0 * kh_670[k];

        t_503[k] = -2.0 * hh_356[k]
                   + f_0 * kh_671[k];

        t_504[k] = -3.0 * hh_357[k]
                   + f_0 * kh_672[k];

        t_505[k] = -3.0 * hh_358[k]
                   + f_0 * kh_673[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, t_509, t_510, hh_359, hh_360, hh_361, hh_362, \
                         hh_363, kh_674, kh_675, kh_676, kh_677, \
                         kh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = -3.0 * hh_359[k]
                   + f_0 * kh_674[k];

        t_507[k] = -3.0 * hh_360[k]
                   + f_0 * kh_675[k];

        t_508[k] = -3.0 * hh_361[k]
                   + f_0 * kh_676[k];

        t_509[k] = -3.0 * hh_362[k]
                   + f_0 * kh_677[k];

        t_510[k] = -3.0 * hh_363[k]
                   + f_0 * kh_678[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, t_515, hh_364, hh_365, hh_366, hh_367, \
                         hh_368, kh_679, kh_680, kh_681, kh_682, \
                         kh_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = -3.0 * hh_364[k]
                   + f_0 * kh_679[k];

        t_512[k] = -3.0 * hh_365[k]
                   + f_0 * kh_680[k];

        t_513[k] = -3.0 * hh_366[k]
                   + f_0 * kh_681[k];

        t_514[k] = -3.0 * hh_367[k]
                   + f_0 * kh_682[k];

        t_515[k] = -3.0 * hh_368[k]
                   + f_0 * kh_683[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, t_520, hh_369, hh_370, hh_371, hh_372, \
                         hh_373, kh_684, kh_685, kh_686, kh_687, \
                         kh_688 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = -3.0 * hh_369[k]
                   + f_0 * kh_684[k];

        t_517[k] = -3.0 * hh_370[k]
                   + f_0 * kh_685[k];

        t_518[k] = -3.0 * hh_371[k]
                   + f_0 * kh_686[k];

        t_519[k] = -3.0 * hh_372[k]
                   + f_0 * kh_687[k];

        t_520[k] = -3.0 * hh_373[k]
                   + f_0 * kh_688[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, t_525, hh_374, hh_375, hh_376, hh_377, \
                         hh_378, kh_689, kh_690, kh_691, kh_692, \
                         kh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = -3.0 * hh_374[k]
                   + f_0 * kh_689[k];

        t_522[k] = -3.0 * hh_375[k]
                   + f_0 * kh_690[k];

        t_523[k] = -3.0 * hh_376[k]
                   + f_0 * kh_691[k];

        t_524[k] = -3.0 * hh_377[k]
                   + f_0 * kh_692[k];

        t_525[k] = -4.0 * hh_378[k]
                   + f_0 * kh_693[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, hh_379, hh_380, hh_381, hh_382, \
                         hh_383, kh_694, kh_695, kh_696, kh_697, \
                         kh_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = -4.0 * hh_379[k]
                   + f_0 * kh_694[k];

        t_527[k] = -4.0 * hh_380[k]
                   + f_0 * kh_695[k];

        t_528[k] = -4.0 * hh_381[k]
                   + f_0 * kh_696[k];

        t_529[k] = -4.0 * hh_382[k]
                   + f_0 * kh_697[k];

        t_530[k] = -4.0 * hh_383[k]
                   + f_0 * kh_698[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, t_535, hh_384, hh_385, hh_386, hh_387, \
                         hh_388, kh_699, kh_700, kh_701, kh_702, \
                         kh_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = -4.0 * hh_384[k]
                   + f_0 * kh_699[k];

        t_532[k] = -4.0 * hh_385[k]
                   + f_0 * kh_700[k];

        t_533[k] = -4.0 * hh_386[k]
                   + f_0 * kh_701[k];

        t_534[k] = -4.0 * hh_387[k]
                   + f_0 * kh_702[k];

        t_535[k] = -4.0 * hh_388[k]
                   + f_0 * kh_703[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, hh_389, hh_390, hh_391, hh_392, \
                         hh_393, kh_704, kh_705, kh_706, kh_707, \
                         kh_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = -4.0 * hh_389[k]
                   + f_0 * kh_704[k];

        t_537[k] = -4.0 * hh_390[k]
                   + f_0 * kh_705[k];

        t_538[k] = -4.0 * hh_391[k]
                   + f_0 * kh_706[k];

        t_539[k] = -4.0 * hh_392[k]
                   + f_0 * kh_707[k];

        t_540[k] = -4.0 * hh_393[k]
                   + f_0 * kh_708[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, t_545, hh_394, hh_395, hh_396, hh_397, \
                         hh_398, kh_709, kh_710, kh_711, kh_712, \
                         kh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = -4.0 * hh_394[k]
                   + f_0 * kh_709[k];

        t_542[k] = -4.0 * hh_395[k]
                   + f_0 * kh_710[k];

        t_543[k] = -4.0 * hh_396[k]
                   + f_0 * kh_711[k];

        t_544[k] = -4.0 * hh_397[k]
                   + f_0 * kh_712[k];

        t_545[k] = -4.0 * hh_398[k]
                   + f_0 * kh_713[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, t_550, hh_399, hh_400, hh_401, hh_402, \
                         hh_403, kh_714, kh_715, kh_716, kh_717, \
                         kh_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = -5.0 * hh_399[k]
                   + f_0 * kh_714[k];

        t_547[k] = -5.0 * hh_400[k]
                   + f_0 * kh_715[k];

        t_548[k] = -5.0 * hh_401[k]
                   + f_0 * kh_716[k];

        t_549[k] = -5.0 * hh_402[k]
                   + f_0 * kh_717[k];

        t_550[k] = -5.0 * hh_403[k]
                   + f_0 * kh_718[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, t_554, t_555, hh_404, hh_405, hh_406, hh_407, \
                         hh_408, kh_719, kh_720, kh_721, kh_722, \
                         kh_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = -5.0 * hh_404[k]
                   + f_0 * kh_719[k];

        t_552[k] = -5.0 * hh_405[k]
                   + f_0 * kh_720[k];

        t_553[k] = -5.0 * hh_406[k]
                   + f_0 * kh_721[k];

        t_554[k] = -5.0 * hh_407[k]
                   + f_0 * kh_722[k];

        t_555[k] = -5.0 * hh_408[k]
                   + f_0 * kh_723[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, hh_409, hh_410, hh_411, hh_412, \
                         hh_413, kh_724, kh_725, kh_726, kh_727, \
                         kh_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = -5.0 * hh_409[k]
                   + f_0 * kh_724[k];

        t_557[k] = -5.0 * hh_410[k]
                   + f_0 * kh_725[k];

        t_558[k] = -5.0 * hh_411[k]
                   + f_0 * kh_726[k];

        t_559[k] = -5.0 * hh_412[k]
                   + f_0 * kh_727[k];

        t_560[k] = -5.0 * hh_413[k]
                   + f_0 * kh_728[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, t_565, hh_414, hh_415, hh_416, hh_417, \
                         hh_418, kh_729, kh_730, kh_731, kh_732, \
                         kh_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = -5.0 * hh_414[k]
                   + f_0 * kh_729[k];

        t_562[k] = -5.0 * hh_415[k]
                   + f_0 * kh_730[k];

        t_563[k] = -5.0 * hh_416[k]
                   + f_0 * kh_731[k];

        t_564[k] = -5.0 * hh_417[k]
                   + f_0 * kh_732[k];

        t_565[k] = -5.0 * hh_418[k]
                   + f_0 * kh_733[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, t_570, hh_419, hh_420, hh_421, hh_422, \
                         hh_423, kh_734, kh_735, kh_736, kh_737, \
                         kh_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = -5.0 * hh_419[k]
                   + f_0 * kh_734[k];

        t_567[k] = -6.0 * hh_420[k]
                   + f_0 * kh_735[k];

        t_568[k] = -6.0 * hh_421[k]
                   + f_0 * kh_736[k];

        t_569[k] = -6.0 * hh_422[k]
                   + f_0 * kh_737[k];

        t_570[k] = -6.0 * hh_423[k]
                   + f_0 * kh_738[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, t_574, t_575, hh_424, hh_425, hh_426, hh_427, \
                         hh_428, kh_739, kh_740, kh_741, kh_742, \
                         kh_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = -6.0 * hh_424[k]
                   + f_0 * kh_739[k];

        t_572[k] = -6.0 * hh_425[k]
                   + f_0 * kh_740[k];

        t_573[k] = -6.0 * hh_426[k]
                   + f_0 * kh_741[k];

        t_574[k] = -6.0 * hh_427[k]
                   + f_0 * kh_742[k];

        t_575[k] = -6.0 * hh_428[k]
                   + f_0 * kh_743[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, hh_429, hh_430, hh_431, hh_432, \
                         hh_433, kh_744, kh_745, kh_746, kh_747, \
                         kh_748 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = -6.0 * hh_429[k]
                   + f_0 * kh_744[k];

        t_577[k] = -6.0 * hh_430[k]
                   + f_0 * kh_745[k];

        t_578[k] = -6.0 * hh_431[k]
                   + f_0 * kh_746[k];

        t_579[k] = -6.0 * hh_432[k]
                   + f_0 * kh_747[k];

        t_580[k] = -6.0 * hh_433[k]
                   + f_0 * kh_748[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, t_584, t_585, hh_434, hh_435, hh_436, hh_437, \
                         hh_438, kh_749, kh_750, kh_751, kh_752, \
                         kh_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = -6.0 * hh_434[k]
                   + f_0 * kh_749[k];

        t_582[k] = -6.0 * hh_435[k]
                   + f_0 * kh_750[k];

        t_583[k] = -6.0 * hh_436[k]
                   + f_0 * kh_751[k];

        t_584[k] = -6.0 * hh_437[k]
                   + f_0 * kh_752[k];

        t_585[k] = -6.0 * hh_438[k]
                   + f_0 * kh_753[k];
    }

#pragma omp simd aligned(t_586, t_587, hh_439, hh_440, kh_754, kh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = -6.0 * hh_439[k]
                   + f_0 * kh_754[k];

        t_587[k] = -6.0 * hh_440[k]
                   + f_0 * kh_755[k];
    }
}

auto
compute_prim_geom_10_ih_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t hh, const size_t kh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ih_electron_repulsion_2_piece0(buffer, target, hh, kh, ncols, alpha);

    compute_prim_geom_10_ih_electron_repulsion_2_piece1(buffer, target, hh, kh, ncols, alpha);

    compute_prim_geom_10_ih_electron_repulsion_2_piece2(buffer, target, hh, kh, ncols, alpha);

    compute_prim_geom_10_ih_electron_repulsion_2_piece3(buffer, target, hh, kh, ncols, alpha);
}

}  // namespace simdt2ceri
