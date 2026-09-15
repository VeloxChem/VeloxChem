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


#include "SimdElectronRepulsionGeom10VrrRecGG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_gg_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fg, const size_t hg,
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

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_64 = buffer.data(fg + 64);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_74 = buffer.data(fg + 74);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_76 = buffer.data(fg + 76);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_81 = buffer.data(fg + 81);
    const auto *fg_82 = buffer.data(fg + 82);
    const auto *fg_83 = buffer.data(fg + 83);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_88 = buffer.data(fg + 88);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_92 = buffer.data(fg + 92);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_94 = buffer.data(fg + 94);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_96 = buffer.data(fg + 96);
    const auto *fg_97 = buffer.data(fg + 97);
    const auto *fg_98 = buffer.data(fg + 98);
    const auto *fg_99 = buffer.data(fg + 99);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_101 = buffer.data(fg + 101);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_105 = buffer.data(fg + 105);
    const auto *fg_106 = buffer.data(fg + 106);
    const auto *fg_107 = buffer.data(fg + 107);
    const auto *fg_108 = buffer.data(fg + 108);
    const auto *fg_109 = buffer.data(fg + 109);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_111 = buffer.data(fg + 111);
    const auto *fg_112 = buffer.data(fg + 112);
    const auto *fg_113 = buffer.data(fg + 113);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_120 = buffer.data(fg + 120);
    const auto *fg_121 = buffer.data(fg + 121);
    const auto *fg_122 = buffer.data(fg + 122);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_124 = buffer.data(fg + 124);
    const auto *fg_125 = buffer.data(fg + 125);
    const auto *fg_126 = buffer.data(fg + 126);
    const auto *fg_127 = buffer.data(fg + 127);
    const auto *fg_128 = buffer.data(fg + 128);
    const auto *fg_129 = buffer.data(fg + 129);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_136 = buffer.data(fg + 136);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_139 = buffer.data(fg + 139);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_141 = buffer.data(fg + 141);
    const auto *fg_142 = buffer.data(fg + 142);
    const auto *fg_143 = buffer.data(fg + 143);
    const auto *fg_144 = buffer.data(fg + 144);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_148 = buffer.data(fg + 148);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_1 = buffer.data(hg + 1);
    const auto *hg_2 = buffer.data(hg + 2);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_4 = buffer.data(hg + 4);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_6 = buffer.data(hg + 6);
    const auto *hg_7 = buffer.data(hg + 7);
    const auto *hg_8 = buffer.data(hg + 8);
    const auto *hg_9 = buffer.data(hg + 9);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_11 = buffer.data(hg + 11);
    const auto *hg_12 = buffer.data(hg + 12);
    const auto *hg_13 = buffer.data(hg + 13);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_16 = buffer.data(hg + 16);
    const auto *hg_17 = buffer.data(hg + 17);
    const auto *hg_18 = buffer.data(hg + 18);
    const auto *hg_19 = buffer.data(hg + 19);
    const auto *hg_20 = buffer.data(hg + 20);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_22 = buffer.data(hg + 22);
    const auto *hg_23 = buffer.data(hg + 23);
    const auto *hg_24 = buffer.data(hg + 24);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_26 = buffer.data(hg + 26);
    const auto *hg_27 = buffer.data(hg + 27);
    const auto *hg_28 = buffer.data(hg + 28);
    const auto *hg_29 = buffer.data(hg + 29);
    const auto *hg_30 = buffer.data(hg + 30);
    const auto *hg_31 = buffer.data(hg + 31);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_33 = buffer.data(hg + 33);
    const auto *hg_34 = buffer.data(hg + 34);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_36 = buffer.data(hg + 36);
    const auto *hg_37 = buffer.data(hg + 37);
    const auto *hg_38 = buffer.data(hg + 38);
    const auto *hg_39 = buffer.data(hg + 39);
    const auto *hg_40 = buffer.data(hg + 40);
    const auto *hg_41 = buffer.data(hg + 41);
    const auto *hg_42 = buffer.data(hg + 42);
    const auto *hg_43 = buffer.data(hg + 43);
    const auto *hg_44 = buffer.data(hg + 44);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_47 = buffer.data(hg + 47);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_49 = buffer.data(hg + 49);
    const auto *hg_50 = buffer.data(hg + 50);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_52 = buffer.data(hg + 52);
    const auto *hg_53 = buffer.data(hg + 53);
    const auto *hg_54 = buffer.data(hg + 54);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_56 = buffer.data(hg + 56);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_58 = buffer.data(hg + 58);
    const auto *hg_59 = buffer.data(hg + 59);
    const auto *hg_60 = buffer.data(hg + 60);
    const auto *hg_61 = buffer.data(hg + 61);
    const auto *hg_62 = buffer.data(hg + 62);
    const auto *hg_63 = buffer.data(hg + 63);
    const auto *hg_64 = buffer.data(hg + 64);
    const auto *hg_65 = buffer.data(hg + 65);
    const auto *hg_66 = buffer.data(hg + 66);
    const auto *hg_67 = buffer.data(hg + 67);
    const auto *hg_68 = buffer.data(hg + 68);
    const auto *hg_69 = buffer.data(hg + 69);
    const auto *hg_70 = buffer.data(hg + 70);
    const auto *hg_71 = buffer.data(hg + 71);
    const auto *hg_72 = buffer.data(hg + 72);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_74 = buffer.data(hg + 74);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_79 = buffer.data(hg + 79);
    const auto *hg_80 = buffer.data(hg + 80);
    const auto *hg_81 = buffer.data(hg + 81);
    const auto *hg_82 = buffer.data(hg + 82);
    const auto *hg_83 = buffer.data(hg + 83);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_85 = buffer.data(hg + 85);
    const auto *hg_86 = buffer.data(hg + 86);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_88 = buffer.data(hg + 88);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_91 = buffer.data(hg + 91);
    const auto *hg_92 = buffer.data(hg + 92);
    const auto *hg_93 = buffer.data(hg + 93);
    const auto *hg_94 = buffer.data(hg + 94);
    const auto *hg_95 = buffer.data(hg + 95);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_97 = buffer.data(hg + 97);
    const auto *hg_98 = buffer.data(hg + 98);
    const auto *hg_99 = buffer.data(hg + 99);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_101 = buffer.data(hg + 101);
    const auto *hg_102 = buffer.data(hg + 102);
    const auto *hg_103 = buffer.data(hg + 103);
    const auto *hg_104 = buffer.data(hg + 104);
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_106 = buffer.data(hg + 106);
    const auto *hg_107 = buffer.data(hg + 107);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_109 = buffer.data(hg + 109);
    const auto *hg_110 = buffer.data(hg + 110);
    const auto *hg_111 = buffer.data(hg + 111);
    const auto *hg_112 = buffer.data(hg + 112);
    const auto *hg_113 = buffer.data(hg + 113);
    const auto *hg_114 = buffer.data(hg + 114);
    const auto *hg_115 = buffer.data(hg + 115);
    const auto *hg_116 = buffer.data(hg + 116);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_118 = buffer.data(hg + 118);
    const auto *hg_119 = buffer.data(hg + 119);
    const auto *hg_120 = buffer.data(hg + 120);
    const auto *hg_121 = buffer.data(hg + 121);
    const auto *hg_122 = buffer.data(hg + 122);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_124 = buffer.data(hg + 124);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_126 = buffer.data(hg + 126);
    const auto *hg_127 = buffer.data(hg + 127);
    const auto *hg_128 = buffer.data(hg + 128);
    const auto *hg_129 = buffer.data(hg + 129);
    const auto *hg_130 = buffer.data(hg + 130);
    const auto *hg_131 = buffer.data(hg + 131);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_133 = buffer.data(hg + 133);
    const auto *hg_134 = buffer.data(hg + 134);
    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_136 = buffer.data(hg + 136);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_138 = buffer.data(hg + 138);
    const auto *hg_139 = buffer.data(hg + 139);
    const auto *hg_140 = buffer.data(hg + 140);
    const auto *hg_141 = buffer.data(hg + 141);
    const auto *hg_142 = buffer.data(hg + 142);
    const auto *hg_143 = buffer.data(hg + 143);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_145 = buffer.data(hg + 145);
    const auto *hg_146 = buffer.data(hg + 146);
    const auto *hg_147 = buffer.data(hg + 147);
    const auto *hg_148 = buffer.data(hg + 148);
    const auto *hg_149 = buffer.data(hg + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, fg_0, fg_1, fg_2, fg_3, fg_4, hg_0, hg_1, \
                         hg_2, hg_3, hg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -4.0 * fg_0[k]
                 + f_0 * hg_0[k];

        t_1[k] = -4.0 * fg_1[k]
                 + f_0 * hg_1[k];

        t_2[k] = -4.0 * fg_2[k]
                 + f_0 * hg_2[k];

        t_3[k] = -4.0 * fg_3[k]
                 + f_0 * hg_3[k];

        t_4[k] = -4.0 * fg_4[k]
                 + f_0 * hg_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, fg_5, fg_6, fg_7, fg_8, fg_9, hg_5, hg_6, \
                         hg_7, hg_8, hg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -4.0 * fg_5[k]
                 + f_0 * hg_5[k];

        t_6[k] = -4.0 * fg_6[k]
                 + f_0 * hg_6[k];

        t_7[k] = -4.0 * fg_7[k]
                 + f_0 * hg_7[k];

        t_8[k] = -4.0 * fg_8[k]
                 + f_0 * hg_8[k];

        t_9[k] = -4.0 * fg_9[k]
                 + f_0 * hg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, fg_10, fg_11, fg_12, fg_13, fg_14, \
                         hg_10, hg_11, hg_12, hg_13, hg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -4.0 * fg_10[k]
                  + f_0 * hg_10[k];

        t_11[k] = -4.0 * fg_11[k]
                  + f_0 * hg_11[k];

        t_12[k] = -4.0 * fg_12[k]
                  + f_0 * hg_12[k];

        t_13[k] = -4.0 * fg_13[k]
                  + f_0 * hg_13[k];

        t_14[k] = -4.0 * fg_14[k]
                  + f_0 * hg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, fg_15, fg_16, fg_17, fg_18, fg_19, \
                         hg_15, hg_16, hg_17, hg_18, hg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -3.0 * fg_15[k]
                  + f_0 * hg_15[k];

        t_16[k] = -3.0 * fg_16[k]
                  + f_0 * hg_16[k];

        t_17[k] = -3.0 * fg_17[k]
                  + f_0 * hg_17[k];

        t_18[k] = -3.0 * fg_18[k]
                  + f_0 * hg_18[k];

        t_19[k] = -3.0 * fg_19[k]
                  + f_0 * hg_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, fg_20, fg_21, fg_22, fg_23, fg_24, \
                         hg_20, hg_21, hg_22, hg_23, hg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -3.0 * fg_20[k]
                  + f_0 * hg_20[k];

        t_21[k] = -3.0 * fg_21[k]
                  + f_0 * hg_21[k];

        t_22[k] = -3.0 * fg_22[k]
                  + f_0 * hg_22[k];

        t_23[k] = -3.0 * fg_23[k]
                  + f_0 * hg_23[k];

        t_24[k] = -3.0 * fg_24[k]
                  + f_0 * hg_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, fg_25, fg_26, fg_27, fg_28, fg_29, \
                         hg_25, hg_26, hg_27, hg_28, hg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -3.0 * fg_25[k]
                  + f_0 * hg_25[k];

        t_26[k] = -3.0 * fg_26[k]
                  + f_0 * hg_26[k];

        t_27[k] = -3.0 * fg_27[k]
                  + f_0 * hg_27[k];

        t_28[k] = -3.0 * fg_28[k]
                  + f_0 * hg_28[k];

        t_29[k] = -3.0 * fg_29[k]
                  + f_0 * hg_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, fg_30, fg_31, fg_32, fg_33, fg_34, \
                         hg_30, hg_31, hg_32, hg_33, hg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -3.0 * fg_30[k]
                  + f_0 * hg_30[k];

        t_31[k] = -3.0 * fg_31[k]
                  + f_0 * hg_31[k];

        t_32[k] = -3.0 * fg_32[k]
                  + f_0 * hg_32[k];

        t_33[k] = -3.0 * fg_33[k]
                  + f_0 * hg_33[k];

        t_34[k] = -3.0 * fg_34[k]
                  + f_0 * hg_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, fg_35, fg_36, fg_37, fg_38, fg_39, \
                         hg_35, hg_36, hg_37, hg_38, hg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -3.0 * fg_35[k]
                  + f_0 * hg_35[k];

        t_36[k] = -3.0 * fg_36[k]
                  + f_0 * hg_36[k];

        t_37[k] = -3.0 * fg_37[k]
                  + f_0 * hg_37[k];

        t_38[k] = -3.0 * fg_38[k]
                  + f_0 * hg_38[k];

        t_39[k] = -3.0 * fg_39[k]
                  + f_0 * hg_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, fg_40, fg_41, fg_42, fg_43, fg_44, \
                         hg_40, hg_41, hg_42, hg_43, hg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * fg_40[k]
                  + f_0 * hg_40[k];

        t_41[k] = -3.0 * fg_41[k]
                  + f_0 * hg_41[k];

        t_42[k] = -3.0 * fg_42[k]
                  + f_0 * hg_42[k];

        t_43[k] = -3.0 * fg_43[k]
                  + f_0 * hg_43[k];

        t_44[k] = -3.0 * fg_44[k]
                  + f_0 * hg_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, fg_45, fg_46, fg_47, fg_48, fg_49, \
                         hg_45, hg_46, hg_47, hg_48, hg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * fg_45[k]
                  + f_0 * hg_45[k];

        t_46[k] = -2.0 * fg_46[k]
                  + f_0 * hg_46[k];

        t_47[k] = -2.0 * fg_47[k]
                  + f_0 * hg_47[k];

        t_48[k] = -2.0 * fg_48[k]
                  + f_0 * hg_48[k];

        t_49[k] = -2.0 * fg_49[k]
                  + f_0 * hg_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, fg_50, fg_51, fg_52, fg_53, fg_54, \
                         hg_50, hg_51, hg_52, hg_53, hg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * fg_50[k]
                  + f_0 * hg_50[k];

        t_51[k] = -2.0 * fg_51[k]
                  + f_0 * hg_51[k];

        t_52[k] = -2.0 * fg_52[k]
                  + f_0 * hg_52[k];

        t_53[k] = -2.0 * fg_53[k]
                  + f_0 * hg_53[k];

        t_54[k] = -2.0 * fg_54[k]
                  + f_0 * hg_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, fg_55, fg_56, fg_57, fg_58, fg_59, \
                         hg_55, hg_56, hg_57, hg_58, hg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * fg_55[k]
                  + f_0 * hg_55[k];

        t_56[k] = -2.0 * fg_56[k]
                  + f_0 * hg_56[k];

        t_57[k] = -2.0 * fg_57[k]
                  + f_0 * hg_57[k];

        t_58[k] = -2.0 * fg_58[k]
                  + f_0 * hg_58[k];

        t_59[k] = -2.0 * fg_59[k]
                  + f_0 * hg_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, fg_60, fg_61, fg_62, fg_63, fg_64, \
                         hg_60, hg_61, hg_62, hg_63, hg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -2.0 * fg_60[k]
                  + f_0 * hg_60[k];

        t_61[k] = -2.0 * fg_61[k]
                  + f_0 * hg_61[k];

        t_62[k] = -2.0 * fg_62[k]
                  + f_0 * hg_62[k];

        t_63[k] = -2.0 * fg_63[k]
                  + f_0 * hg_63[k];

        t_64[k] = -2.0 * fg_64[k]
                  + f_0 * hg_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, fg_65, fg_66, fg_67, fg_68, fg_69, \
                         hg_65, hg_66, hg_67, hg_68, hg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -2.0 * fg_65[k]
                  + f_0 * hg_65[k];

        t_66[k] = -2.0 * fg_66[k]
                  + f_0 * hg_66[k];

        t_67[k] = -2.0 * fg_67[k]
                  + f_0 * hg_67[k];

        t_68[k] = -2.0 * fg_68[k]
                  + f_0 * hg_68[k];

        t_69[k] = -2.0 * fg_69[k]
                  + f_0 * hg_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, fg_70, fg_71, fg_72, fg_73, fg_74, \
                         hg_70, hg_71, hg_72, hg_73, hg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * fg_70[k]
                  + f_0 * hg_70[k];

        t_71[k] = -2.0 * fg_71[k]
                  + f_0 * hg_71[k];

        t_72[k] = -2.0 * fg_72[k]
                  + f_0 * hg_72[k];

        t_73[k] = -2.0 * fg_73[k]
                  + f_0 * hg_73[k];

        t_74[k] = -2.0 * fg_74[k]
                  + f_0 * hg_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, fg_75, fg_76, fg_77, fg_78, fg_79, \
                         hg_75, hg_76, hg_77, hg_78, hg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * fg_75[k]
                  + f_0 * hg_75[k];

        t_76[k] = -2.0 * fg_76[k]
                  + f_0 * hg_76[k];

        t_77[k] = -2.0 * fg_77[k]
                  + f_0 * hg_77[k];

        t_78[k] = -2.0 * fg_78[k]
                  + f_0 * hg_78[k];

        t_79[k] = -2.0 * fg_79[k]
                  + f_0 * hg_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, fg_80, fg_81, fg_82, fg_83, fg_84, \
                         hg_80, hg_81, hg_82, hg_83, hg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -2.0 * fg_80[k]
                  + f_0 * hg_80[k];

        t_81[k] = -2.0 * fg_81[k]
                  + f_0 * hg_81[k];

        t_82[k] = -2.0 * fg_82[k]
                  + f_0 * hg_82[k];

        t_83[k] = -2.0 * fg_83[k]
                  + f_0 * hg_83[k];

        t_84[k] = -2.0 * fg_84[k]
                  + f_0 * hg_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, fg_85, fg_86, fg_87, fg_88, fg_89, \
                         hg_85, hg_86, hg_87, hg_88, hg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * fg_85[k]
                  + f_0 * hg_85[k];

        t_86[k] = -2.0 * fg_86[k]
                  + f_0 * hg_86[k];

        t_87[k] = -2.0 * fg_87[k]
                  + f_0 * hg_87[k];

        t_88[k] = -2.0 * fg_88[k]
                  + f_0 * hg_88[k];

        t_89[k] = -2.0 * fg_89[k]
                  + f_0 * hg_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, fg_90, fg_91, fg_92, fg_93, fg_94, \
                         hg_90, hg_91, hg_92, hg_93, hg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -fg_90[k]
                  + f_0 * hg_90[k];

        t_91[k] = -fg_91[k]
                  + f_0 * hg_91[k];

        t_92[k] = -fg_92[k]
                  + f_0 * hg_92[k];

        t_93[k] = -fg_93[k]
                  + f_0 * hg_93[k];

        t_94[k] = -fg_94[k]
                  + f_0 * hg_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, fg_95, fg_96, fg_97, fg_98, fg_99, \
                         hg_95, hg_96, hg_97, hg_98, hg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -fg_95[k]
                  + f_0 * hg_95[k];

        t_96[k] = -fg_96[k]
                  + f_0 * hg_96[k];

        t_97[k] = -fg_97[k]
                  + f_0 * hg_97[k];

        t_98[k] = -fg_98[k]
                  + f_0 * hg_98[k];

        t_99[k] = -fg_99[k]
                  + f_0 * hg_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, fg_100, fg_101, fg_102, fg_103, \
                         fg_104, hg_100, hg_101, hg_102, hg_103, \
                         hg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -fg_100[k]
                   + f_0 * hg_100[k];

        t_101[k] = -fg_101[k]
                   + f_0 * hg_101[k];

        t_102[k] = -fg_102[k]
                   + f_0 * hg_102[k];

        t_103[k] = -fg_103[k]
                   + f_0 * hg_103[k];

        t_104[k] = -fg_104[k]
                   + f_0 * hg_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, fg_105, fg_106, fg_107, fg_108, \
                         fg_109, hg_105, hg_106, hg_107, hg_108, \
                         hg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -fg_105[k]
                   + f_0 * hg_105[k];

        t_106[k] = -fg_106[k]
                   + f_0 * hg_106[k];

        t_107[k] = -fg_107[k]
                   + f_0 * hg_107[k];

        t_108[k] = -fg_108[k]
                   + f_0 * hg_108[k];

        t_109[k] = -fg_109[k]
                   + f_0 * hg_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, fg_110, fg_111, fg_112, fg_113, \
                         fg_114, hg_110, hg_111, hg_112, hg_113, \
                         hg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -fg_110[k]
                   + f_0 * hg_110[k];

        t_111[k] = -fg_111[k]
                   + f_0 * hg_111[k];

        t_112[k] = -fg_112[k]
                   + f_0 * hg_112[k];

        t_113[k] = -fg_113[k]
                   + f_0 * hg_113[k];

        t_114[k] = -fg_114[k]
                   + f_0 * hg_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, fg_115, fg_116, fg_117, fg_118, \
                         fg_119, hg_115, hg_116, hg_117, hg_118, \
                         hg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -fg_115[k]
                   + f_0 * hg_115[k];

        t_116[k] = -fg_116[k]
                   + f_0 * hg_116[k];

        t_117[k] = -fg_117[k]
                   + f_0 * hg_117[k];

        t_118[k] = -fg_118[k]
                   + f_0 * hg_118[k];

        t_119[k] = -fg_119[k]
                   + f_0 * hg_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, fg_120, fg_121, fg_122, fg_123, \
                         fg_124, hg_120, hg_121, hg_122, hg_123, \
                         hg_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -fg_120[k]
                   + f_0 * hg_120[k];

        t_121[k] = -fg_121[k]
                   + f_0 * hg_121[k];

        t_122[k] = -fg_122[k]
                   + f_0 * hg_122[k];

        t_123[k] = -fg_123[k]
                   + f_0 * hg_123[k];

        t_124[k] = -fg_124[k]
                   + f_0 * hg_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, fg_125, fg_126, fg_127, fg_128, \
                         fg_129, hg_125, hg_126, hg_127, hg_128, \
                         hg_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -fg_125[k]
                   + f_0 * hg_125[k];

        t_126[k] = -fg_126[k]
                   + f_0 * hg_126[k];

        t_127[k] = -fg_127[k]
                   + f_0 * hg_127[k];

        t_128[k] = -fg_128[k]
                   + f_0 * hg_128[k];

        t_129[k] = -fg_129[k]
                   + f_0 * hg_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, fg_130, fg_131, fg_132, fg_133, \
                         fg_134, hg_130, hg_131, hg_132, hg_133, \
                         hg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -fg_130[k]
                   + f_0 * hg_130[k];

        t_131[k] = -fg_131[k]
                   + f_0 * hg_131[k];

        t_132[k] = -fg_132[k]
                   + f_0 * hg_132[k];

        t_133[k] = -fg_133[k]
                   + f_0 * hg_133[k];

        t_134[k] = -fg_134[k]
                   + f_0 * hg_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, fg_135, fg_136, fg_137, fg_138, \
                         fg_139, hg_135, hg_136, hg_137, hg_138, \
                         hg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -fg_135[k]
                   + f_0 * hg_135[k];

        t_136[k] = -fg_136[k]
                   + f_0 * hg_136[k];

        t_137[k] = -fg_137[k]
                   + f_0 * hg_137[k];

        t_138[k] = -fg_138[k]
                   + f_0 * hg_138[k];

        t_139[k] = -fg_139[k]
                   + f_0 * hg_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, fg_140, fg_141, fg_142, fg_143, \
                         fg_144, hg_140, hg_141, hg_142, hg_143, \
                         hg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -fg_140[k]
                   + f_0 * hg_140[k];

        t_141[k] = -fg_141[k]
                   + f_0 * hg_141[k];

        t_142[k] = -fg_142[k]
                   + f_0 * hg_142[k];

        t_143[k] = -fg_143[k]
                   + f_0 * hg_143[k];

        t_144[k] = -fg_144[k]
                   + f_0 * hg_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, fg_145, fg_146, fg_147, fg_148, \
                         fg_149, hg_145, hg_146, hg_147, hg_148, \
                         hg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -fg_145[k]
                   + f_0 * hg_145[k];

        t_146[k] = -fg_146[k]
                   + f_0 * hg_146[k];

        t_147[k] = -fg_147[k]
                   + f_0 * hg_147[k];

        t_148[k] = -fg_148[k]
                   + f_0 * hg_148[k];

        t_149[k] = -fg_149[k]
                   + f_0 * hg_149[k];
    }
}

static auto
compute_prim_geom_10_gg_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hg, const size_t ncols,
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

    const auto *hg_150 = buffer.data(hg + 150);
    const auto *hg_151 = buffer.data(hg + 151);
    const auto *hg_152 = buffer.data(hg + 152);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_154 = buffer.data(hg + 154);
    const auto *hg_155 = buffer.data(hg + 155);
    const auto *hg_156 = buffer.data(hg + 156);
    const auto *hg_157 = buffer.data(hg + 157);
    const auto *hg_158 = buffer.data(hg + 158);
    const auto *hg_159 = buffer.data(hg + 159);
    const auto *hg_160 = buffer.data(hg + 160);
    const auto *hg_161 = buffer.data(hg + 161);
    const auto *hg_162 = buffer.data(hg + 162);
    const auto *hg_163 = buffer.data(hg + 163);
    const auto *hg_164 = buffer.data(hg + 164);
    const auto *hg_165 = buffer.data(hg + 165);
    const auto *hg_166 = buffer.data(hg + 166);
    const auto *hg_167 = buffer.data(hg + 167);
    const auto *hg_168 = buffer.data(hg + 168);
    const auto *hg_169 = buffer.data(hg + 169);
    const auto *hg_170 = buffer.data(hg + 170);
    const auto *hg_171 = buffer.data(hg + 171);
    const auto *hg_172 = buffer.data(hg + 172);
    const auto *hg_173 = buffer.data(hg + 173);
    const auto *hg_174 = buffer.data(hg + 174);
    const auto *hg_175 = buffer.data(hg + 175);
    const auto *hg_176 = buffer.data(hg + 176);
    const auto *hg_177 = buffer.data(hg + 177);
    const auto *hg_178 = buffer.data(hg + 178);
    const auto *hg_179 = buffer.data(hg + 179);
    const auto *hg_180 = buffer.data(hg + 180);
    const auto *hg_181 = buffer.data(hg + 181);
    const auto *hg_182 = buffer.data(hg + 182);
    const auto *hg_183 = buffer.data(hg + 183);
    const auto *hg_184 = buffer.data(hg + 184);
    const auto *hg_185 = buffer.data(hg + 185);
    const auto *hg_186 = buffer.data(hg + 186);
    const auto *hg_187 = buffer.data(hg + 187);
    const auto *hg_188 = buffer.data(hg + 188);
    const auto *hg_189 = buffer.data(hg + 189);
    const auto *hg_190 = buffer.data(hg + 190);
    const auto *hg_191 = buffer.data(hg + 191);
    const auto *hg_192 = buffer.data(hg + 192);
    const auto *hg_193 = buffer.data(hg + 193);
    const auto *hg_194 = buffer.data(hg + 194);
    const auto *hg_195 = buffer.data(hg + 195);
    const auto *hg_196 = buffer.data(hg + 196);
    const auto *hg_197 = buffer.data(hg + 197);
    const auto *hg_198 = buffer.data(hg + 198);
    const auto *hg_199 = buffer.data(hg + 199);
    const auto *hg_200 = buffer.data(hg + 200);
    const auto *hg_201 = buffer.data(hg + 201);
    const auto *hg_202 = buffer.data(hg + 202);
    const auto *hg_203 = buffer.data(hg + 203);
    const auto *hg_204 = buffer.data(hg + 204);
    const auto *hg_205 = buffer.data(hg + 205);
    const auto *hg_206 = buffer.data(hg + 206);
    const auto *hg_207 = buffer.data(hg + 207);
    const auto *hg_208 = buffer.data(hg + 208);
    const auto *hg_209 = buffer.data(hg + 209);
    const auto *hg_210 = buffer.data(hg + 210);
    const auto *hg_211 = buffer.data(hg + 211);
    const auto *hg_212 = buffer.data(hg + 212);
    const auto *hg_213 = buffer.data(hg + 213);
    const auto *hg_214 = buffer.data(hg + 214);
    const auto *hg_215 = buffer.data(hg + 215);
    const auto *hg_216 = buffer.data(hg + 216);
    const auto *hg_217 = buffer.data(hg + 217);
    const auto *hg_218 = buffer.data(hg + 218);
    const auto *hg_219 = buffer.data(hg + 219);
    const auto *hg_220 = buffer.data(hg + 220);
    const auto *hg_221 = buffer.data(hg + 221);
    const auto *hg_222 = buffer.data(hg + 222);
    const auto *hg_223 = buffer.data(hg + 223);
    const auto *hg_224 = buffer.data(hg + 224);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, t_156, t_157, hg_150, \
                         hg_151, hg_152, hg_153, hg_154, hg_155, hg_156, \
                         hg_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_0 * hg_150[k];

        t_151[k] = f_0 * hg_151[k];

        t_152[k] = f_0 * hg_152[k];

        t_153[k] = f_0 * hg_153[k];

        t_154[k] = f_0 * hg_154[k];

        t_155[k] = f_0 * hg_155[k];

        t_156[k] = f_0 * hg_156[k];

        t_157[k] = f_0 * hg_157[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, t_163, t_164, t_165, hg_158, \
                         hg_159, hg_160, hg_161, hg_162, hg_163, hg_164, \
                         hg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_0 * hg_158[k];

        t_159[k] = f_0 * hg_159[k];

        t_160[k] = f_0 * hg_160[k];

        t_161[k] = f_0 * hg_161[k];

        t_162[k] = f_0 * hg_162[k];

        t_163[k] = f_0 * hg_163[k];

        t_164[k] = f_0 * hg_164[k];

        t_165[k] = f_0 * hg_165[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, t_171, t_172, t_173, hg_166, \
                         hg_167, hg_168, hg_169, hg_170, hg_171, hg_172, \
                         hg_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_0 * hg_166[k];

        t_167[k] = f_0 * hg_167[k];

        t_168[k] = f_0 * hg_168[k];

        t_169[k] = f_0 * hg_169[k];

        t_170[k] = f_0 * hg_170[k];

        t_171[k] = f_0 * hg_171[k];

        t_172[k] = f_0 * hg_172[k];

        t_173[k] = f_0 * hg_173[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, t_179, t_180, t_181, hg_174, \
                         hg_175, hg_176, hg_177, hg_178, hg_179, hg_180, \
                         hg_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_0 * hg_174[k];

        t_175[k] = f_0 * hg_175[k];

        t_176[k] = f_0 * hg_176[k];

        t_177[k] = f_0 * hg_177[k];

        t_178[k] = f_0 * hg_178[k];

        t_179[k] = f_0 * hg_179[k];

        t_180[k] = f_0 * hg_180[k];

        t_181[k] = f_0 * hg_181[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, t_187, t_188, t_189, hg_182, \
                         hg_183, hg_184, hg_185, hg_186, hg_187, hg_188, \
                         hg_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_0 * hg_182[k];

        t_183[k] = f_0 * hg_183[k];

        t_184[k] = f_0 * hg_184[k];

        t_185[k] = f_0 * hg_185[k];

        t_186[k] = f_0 * hg_186[k];

        t_187[k] = f_0 * hg_187[k];

        t_188[k] = f_0 * hg_188[k];

        t_189[k] = f_0 * hg_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, t_195, t_196, t_197, hg_190, \
                         hg_191, hg_192, hg_193, hg_194, hg_195, hg_196, \
                         hg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_0 * hg_190[k];

        t_191[k] = f_0 * hg_191[k];

        t_192[k] = f_0 * hg_192[k];

        t_193[k] = f_0 * hg_193[k];

        t_194[k] = f_0 * hg_194[k];

        t_195[k] = f_0 * hg_195[k];

        t_196[k] = f_0 * hg_196[k];

        t_197[k] = f_0 * hg_197[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, t_203, t_204, t_205, hg_198, \
                         hg_199, hg_200, hg_201, hg_202, hg_203, hg_204, \
                         hg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_0 * hg_198[k];

        t_199[k] = f_0 * hg_199[k];

        t_200[k] = f_0 * hg_200[k];

        t_201[k] = f_0 * hg_201[k];

        t_202[k] = f_0 * hg_202[k];

        t_203[k] = f_0 * hg_203[k];

        t_204[k] = f_0 * hg_204[k];

        t_205[k] = f_0 * hg_205[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, t_211, t_212, t_213, hg_206, \
                         hg_207, hg_208, hg_209, hg_210, hg_211, hg_212, \
                         hg_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_0 * hg_206[k];

        t_207[k] = f_0 * hg_207[k];

        t_208[k] = f_0 * hg_208[k];

        t_209[k] = f_0 * hg_209[k];

        t_210[k] = f_0 * hg_210[k];

        t_211[k] = f_0 * hg_211[k];

        t_212[k] = f_0 * hg_212[k];

        t_213[k] = f_0 * hg_213[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, t_219, t_220, t_221, hg_214, \
                         hg_215, hg_216, hg_217, hg_218, hg_219, hg_220, \
                         hg_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_0 * hg_214[k];

        t_215[k] = f_0 * hg_215[k];

        t_216[k] = f_0 * hg_216[k];

        t_217[k] = f_0 * hg_217[k];

        t_218[k] = f_0 * hg_218[k];

        t_219[k] = f_0 * hg_219[k];

        t_220[k] = f_0 * hg_220[k];

        t_221[k] = f_0 * hg_221[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, hg_222, hg_223, hg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_0 * hg_222[k];

        t_223[k] = f_0 * hg_223[k];

        t_224[k] = f_0 * hg_224[k];
    }
}

auto
compute_prim_geom_10_gg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t fg, const size_t hg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gg_electron_repulsion_0_piece0(buffer, target, fg, hg, ncols, alpha);

    compute_prim_geom_10_gg_electron_repulsion_0_piece1(buffer, target, hg, ncols, alpha);
}

static auto
compute_prim_geom_10_gg_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fg, const size_t hg,
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

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_64 = buffer.data(fg + 64);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_74 = buffer.data(fg + 74);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_76 = buffer.data(fg + 76);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_81 = buffer.data(fg + 81);
    const auto *fg_82 = buffer.data(fg + 82);
    const auto *fg_83 = buffer.data(fg + 83);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_88 = buffer.data(fg + 88);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_92 = buffer.data(fg + 92);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_94 = buffer.data(fg + 94);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_96 = buffer.data(fg + 96);
    const auto *fg_97 = buffer.data(fg + 97);
    const auto *fg_98 = buffer.data(fg + 98);
    const auto *fg_99 = buffer.data(fg + 99);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_101 = buffer.data(fg + 101);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_105 = buffer.data(fg + 105);
    const auto *fg_106 = buffer.data(fg + 106);
    const auto *fg_107 = buffer.data(fg + 107);
    const auto *fg_108 = buffer.data(fg + 108);
    const auto *fg_109 = buffer.data(fg + 109);

    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_16 = buffer.data(hg + 16);
    const auto *hg_17 = buffer.data(hg + 17);
    const auto *hg_18 = buffer.data(hg + 18);
    const auto *hg_19 = buffer.data(hg + 19);
    const auto *hg_20 = buffer.data(hg + 20);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_22 = buffer.data(hg + 22);
    const auto *hg_23 = buffer.data(hg + 23);
    const auto *hg_24 = buffer.data(hg + 24);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_26 = buffer.data(hg + 26);
    const auto *hg_27 = buffer.data(hg + 27);
    const auto *hg_28 = buffer.data(hg + 28);
    const auto *hg_29 = buffer.data(hg + 29);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_47 = buffer.data(hg + 47);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_49 = buffer.data(hg + 49);
    const auto *hg_50 = buffer.data(hg + 50);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_52 = buffer.data(hg + 52);
    const auto *hg_53 = buffer.data(hg + 53);
    const auto *hg_54 = buffer.data(hg + 54);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_56 = buffer.data(hg + 56);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_58 = buffer.data(hg + 58);
    const auto *hg_59 = buffer.data(hg + 59);
    const auto *hg_60 = buffer.data(hg + 60);
    const auto *hg_61 = buffer.data(hg + 61);
    const auto *hg_62 = buffer.data(hg + 62);
    const auto *hg_63 = buffer.data(hg + 63);
    const auto *hg_64 = buffer.data(hg + 64);
    const auto *hg_65 = buffer.data(hg + 65);
    const auto *hg_66 = buffer.data(hg + 66);
    const auto *hg_67 = buffer.data(hg + 67);
    const auto *hg_68 = buffer.data(hg + 68);
    const auto *hg_69 = buffer.data(hg + 69);
    const auto *hg_70 = buffer.data(hg + 70);
    const auto *hg_71 = buffer.data(hg + 71);
    const auto *hg_72 = buffer.data(hg + 72);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_74 = buffer.data(hg + 74);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_91 = buffer.data(hg + 91);
    const auto *hg_92 = buffer.data(hg + 92);
    const auto *hg_93 = buffer.data(hg + 93);
    const auto *hg_94 = buffer.data(hg + 94);
    const auto *hg_95 = buffer.data(hg + 95);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_97 = buffer.data(hg + 97);
    const auto *hg_98 = buffer.data(hg + 98);
    const auto *hg_99 = buffer.data(hg + 99);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_101 = buffer.data(hg + 101);
    const auto *hg_102 = buffer.data(hg + 102);
    const auto *hg_103 = buffer.data(hg + 103);
    const auto *hg_104 = buffer.data(hg + 104);
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_106 = buffer.data(hg + 106);
    const auto *hg_107 = buffer.data(hg + 107);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_109 = buffer.data(hg + 109);
    const auto *hg_110 = buffer.data(hg + 110);
    const auto *hg_111 = buffer.data(hg + 111);
    const auto *hg_112 = buffer.data(hg + 112);
    const auto *hg_113 = buffer.data(hg + 113);
    const auto *hg_114 = buffer.data(hg + 114);
    const auto *hg_115 = buffer.data(hg + 115);
    const auto *hg_116 = buffer.data(hg + 116);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_118 = buffer.data(hg + 118);
    const auto *hg_119 = buffer.data(hg + 119);
    const auto *hg_120 = buffer.data(hg + 120);
    const auto *hg_121 = buffer.data(hg + 121);
    const auto *hg_122 = buffer.data(hg + 122);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_124 = buffer.data(hg + 124);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_126 = buffer.data(hg + 126);
    const auto *hg_127 = buffer.data(hg + 127);
    const auto *hg_128 = buffer.data(hg + 128);
    const auto *hg_129 = buffer.data(hg + 129);
    const auto *hg_130 = buffer.data(hg + 130);
    const auto *hg_131 = buffer.data(hg + 131);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_133 = buffer.data(hg + 133);
    const auto *hg_134 = buffer.data(hg + 134);
    const auto *hg_150 = buffer.data(hg + 150);
    const auto *hg_151 = buffer.data(hg + 151);
    const auto *hg_152 = buffer.data(hg + 152);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_154 = buffer.data(hg + 154);
    const auto *hg_155 = buffer.data(hg + 155);
    const auto *hg_156 = buffer.data(hg + 156);
    const auto *hg_157 = buffer.data(hg + 157);
    const auto *hg_158 = buffer.data(hg + 158);
    const auto *hg_159 = buffer.data(hg + 159);
    const auto *hg_160 = buffer.data(hg + 160);
    const auto *hg_161 = buffer.data(hg + 161);
    const auto *hg_162 = buffer.data(hg + 162);
    const auto *hg_163 = buffer.data(hg + 163);
    const auto *hg_164 = buffer.data(hg + 164);
    const auto *hg_165 = buffer.data(hg + 165);
    const auto *hg_166 = buffer.data(hg + 166);
    const auto *hg_167 = buffer.data(hg + 167);
    const auto *hg_168 = buffer.data(hg + 168);
    const auto *hg_169 = buffer.data(hg + 169);
    const auto *hg_170 = buffer.data(hg + 170);
    const auto *hg_171 = buffer.data(hg + 171);
    const auto *hg_172 = buffer.data(hg + 172);
    const auto *hg_173 = buffer.data(hg + 173);
    const auto *hg_174 = buffer.data(hg + 174);
    const auto *hg_175 = buffer.data(hg + 175);
    const auto *hg_176 = buffer.data(hg + 176);
    const auto *hg_177 = buffer.data(hg + 177);
    const auto *hg_178 = buffer.data(hg + 178);
    const auto *hg_179 = buffer.data(hg + 179);
    const auto *hg_180 = buffer.data(hg + 180);
    const auto *hg_181 = buffer.data(hg + 181);
    const auto *hg_182 = buffer.data(hg + 182);
    const auto *hg_183 = buffer.data(hg + 183);
    const auto *hg_184 = buffer.data(hg + 184);
    const auto *hg_185 = buffer.data(hg + 185);
    const auto *hg_186 = buffer.data(hg + 186);
    const auto *hg_187 = buffer.data(hg + 187);
    const auto *hg_188 = buffer.data(hg + 188);
    const auto *hg_189 = buffer.data(hg + 189);
    const auto *hg_190 = buffer.data(hg + 190);
    const auto *hg_191 = buffer.data(hg + 191);
    const auto *hg_192 = buffer.data(hg + 192);
    const auto *hg_193 = buffer.data(hg + 193);
    const auto *hg_194 = buffer.data(hg + 194);
    const auto *hg_195 = buffer.data(hg + 195);
    const auto *hg_196 = buffer.data(hg + 196);
    const auto *hg_197 = buffer.data(hg + 197);
    const auto *hg_198 = buffer.data(hg + 198);
    const auto *hg_199 = buffer.data(hg + 199);
    const auto *hg_200 = buffer.data(hg + 200);
    const auto *hg_201 = buffer.data(hg + 201);
    const auto *hg_202 = buffer.data(hg + 202);
    const auto *hg_203 = buffer.data(hg + 203);
    const auto *hg_204 = buffer.data(hg + 204);
    const auto *hg_205 = buffer.data(hg + 205);
    const auto *hg_206 = buffer.data(hg + 206);
    const auto *hg_207 = buffer.data(hg + 207);
    const auto *hg_208 = buffer.data(hg + 208);
    const auto *hg_209 = buffer.data(hg + 209);
    const auto *hg_225 = buffer.data(hg + 225);
    const auto *hg_226 = buffer.data(hg + 226);
    const auto *hg_227 = buffer.data(hg + 227);
    const auto *hg_228 = buffer.data(hg + 228);
    const auto *hg_229 = buffer.data(hg + 229);
    const auto *hg_230 = buffer.data(hg + 230);
    const auto *hg_231 = buffer.data(hg + 231);
    const auto *hg_232 = buffer.data(hg + 232);
    const auto *hg_233 = buffer.data(hg + 233);
    const auto *hg_234 = buffer.data(hg + 234);
    const auto *hg_235 = buffer.data(hg + 235);
    const auto *hg_236 = buffer.data(hg + 236);
    const auto *hg_237 = buffer.data(hg + 237);
    const auto *hg_238 = buffer.data(hg + 238);
    const auto *hg_239 = buffer.data(hg + 239);
    const auto *hg_240 = buffer.data(hg + 240);
    const auto *hg_241 = buffer.data(hg + 241);
    const auto *hg_242 = buffer.data(hg + 242);
    const auto *hg_243 = buffer.data(hg + 243);
    const auto *hg_244 = buffer.data(hg + 244);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hg_15, hg_16, hg_17, hg_18, \
                         hg_19, hg_20, hg_21, hg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hg_15[k];

        t_1[k] = f_0 * hg_16[k];

        t_2[k] = f_0 * hg_17[k];

        t_3[k] = f_0 * hg_18[k];

        t_4[k] = f_0 * hg_19[k];

        t_5[k] = f_0 * hg_20[k];

        t_6[k] = f_0 * hg_21[k];

        t_7[k] = f_0 * hg_22[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, hg_23, hg_24, hg_25, hg_26, \
                         hg_27, hg_28, hg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hg_23[k];

        t_9[k] = f_0 * hg_24[k];

        t_10[k] = f_0 * hg_25[k];

        t_11[k] = f_0 * hg_26[k];

        t_12[k] = f_0 * hg_27[k];

        t_13[k] = f_0 * hg_28[k];

        t_14[k] = f_0 * hg_29[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, fg_0, fg_1, fg_2, fg_3, fg_4, hg_45, \
                         hg_46, hg_47, hg_48, hg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -fg_0[k]
                  + f_0 * hg_45[k];

        t_16[k] = -fg_1[k]
                  + f_0 * hg_46[k];

        t_17[k] = -fg_2[k]
                  + f_0 * hg_47[k];

        t_18[k] = -fg_3[k]
                  + f_0 * hg_48[k];

        t_19[k] = -fg_4[k]
                  + f_0 * hg_49[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, fg_5, fg_6, fg_7, fg_8, fg_9, hg_50, \
                         hg_51, hg_52, hg_53, hg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -fg_5[k]
                  + f_0 * hg_50[k];

        t_21[k] = -fg_6[k]
                  + f_0 * hg_51[k];

        t_22[k] = -fg_7[k]
                  + f_0 * hg_52[k];

        t_23[k] = -fg_8[k]
                  + f_0 * hg_53[k];

        t_24[k] = -fg_9[k]
                  + f_0 * hg_54[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, fg_10, fg_11, fg_12, fg_13, fg_14, \
                         hg_55, hg_56, hg_57, hg_58, hg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -fg_10[k]
                  + f_0 * hg_55[k];

        t_26[k] = -fg_11[k]
                  + f_0 * hg_56[k];

        t_27[k] = -fg_12[k]
                  + f_0 * hg_57[k];

        t_28[k] = -fg_13[k]
                  + f_0 * hg_58[k];

        t_29[k] = -fg_14[k]
                  + f_0 * hg_59[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, t_36, t_37, hg_60, hg_61, hg_62, \
                         hg_63, hg_64, hg_65, hg_66, hg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * hg_60[k];

        t_31[k] = f_0 * hg_61[k];

        t_32[k] = f_0 * hg_62[k];

        t_33[k] = f_0 * hg_63[k];

        t_34[k] = f_0 * hg_64[k];

        t_35[k] = f_0 * hg_65[k];

        t_36[k] = f_0 * hg_66[k];

        t_37[k] = f_0 * hg_67[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, hg_68, hg_69, hg_70, hg_71, \
                         hg_72, hg_73, hg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * hg_68[k];

        t_39[k] = f_0 * hg_69[k];

        t_40[k] = f_0 * hg_70[k];

        t_41[k] = f_0 * hg_71[k];

        t_42[k] = f_0 * hg_72[k];

        t_43[k] = f_0 * hg_73[k];

        t_44[k] = f_0 * hg_74[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, fg_15, fg_16, fg_17, fg_18, fg_19, \
                         hg_90, hg_91, hg_92, hg_93, hg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * fg_15[k]
                  + f_0 * hg_90[k];

        t_46[k] = -2.0 * fg_16[k]
                  + f_0 * hg_91[k];

        t_47[k] = -2.0 * fg_17[k]
                  + f_0 * hg_92[k];

        t_48[k] = -2.0 * fg_18[k]
                  + f_0 * hg_93[k];

        t_49[k] = -2.0 * fg_19[k]
                  + f_0 * hg_94[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, fg_20, fg_21, fg_22, fg_23, fg_24, \
                         hg_95, hg_96, hg_97, hg_98, hg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * fg_20[k]
                  + f_0 * hg_95[k];

        t_51[k] = -2.0 * fg_21[k]
                  + f_0 * hg_96[k];

        t_52[k] = -2.0 * fg_22[k]
                  + f_0 * hg_97[k];

        t_53[k] = -2.0 * fg_23[k]
                  + f_0 * hg_98[k];

        t_54[k] = -2.0 * fg_24[k]
                  + f_0 * hg_99[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, fg_25, fg_26, fg_27, fg_28, fg_29, \
                         hg_100, hg_101, hg_102, hg_103, hg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * fg_25[k]
                  + f_0 * hg_100[k];

        t_56[k] = -2.0 * fg_26[k]
                  + f_0 * hg_101[k];

        t_57[k] = -2.0 * fg_27[k]
                  + f_0 * hg_102[k];

        t_58[k] = -2.0 * fg_28[k]
                  + f_0 * hg_103[k];

        t_59[k] = -2.0 * fg_29[k]
                  + f_0 * hg_104[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, fg_30, fg_31, fg_32, fg_33, fg_34, \
                         hg_105, hg_106, hg_107, hg_108, hg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -fg_30[k]
                  + f_0 * hg_105[k];

        t_61[k] = -fg_31[k]
                  + f_0 * hg_106[k];

        t_62[k] = -fg_32[k]
                  + f_0 * hg_107[k];

        t_63[k] = -fg_33[k]
                  + f_0 * hg_108[k];

        t_64[k] = -fg_34[k]
                  + f_0 * hg_109[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, fg_35, fg_36, fg_37, fg_38, fg_39, \
                         hg_110, hg_111, hg_112, hg_113, hg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -fg_35[k]
                  + f_0 * hg_110[k];

        t_66[k] = -fg_36[k]
                  + f_0 * hg_111[k];

        t_67[k] = -fg_37[k]
                  + f_0 * hg_112[k];

        t_68[k] = -fg_38[k]
                  + f_0 * hg_113[k];

        t_69[k] = -fg_39[k]
                  + f_0 * hg_114[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, fg_40, fg_41, fg_42, fg_43, fg_44, \
                         hg_115, hg_116, hg_117, hg_118, hg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -fg_40[k]
                  + f_0 * hg_115[k];

        t_71[k] = -fg_41[k]
                  + f_0 * hg_116[k];

        t_72[k] = -fg_42[k]
                  + f_0 * hg_117[k];

        t_73[k] = -fg_43[k]
                  + f_0 * hg_118[k];

        t_74[k] = -fg_44[k]
                  + f_0 * hg_119[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, t_80, t_81, t_82, hg_120, hg_121, \
                         hg_122, hg_123, hg_124, hg_125, hg_126, \
                         hg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * hg_120[k];

        t_76[k] = f_0 * hg_121[k];

        t_77[k] = f_0 * hg_122[k];

        t_78[k] = f_0 * hg_123[k];

        t_79[k] = f_0 * hg_124[k];

        t_80[k] = f_0 * hg_125[k];

        t_81[k] = f_0 * hg_126[k];

        t_82[k] = f_0 * hg_127[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, hg_128, hg_129, hg_130, \
                         hg_131, hg_132, hg_133, hg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * hg_128[k];

        t_84[k] = f_0 * hg_129[k];

        t_85[k] = f_0 * hg_130[k];

        t_86[k] = f_0 * hg_131[k];

        t_87[k] = f_0 * hg_132[k];

        t_88[k] = f_0 * hg_133[k];

        t_89[k] = f_0 * hg_134[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, fg_45, fg_46, fg_47, fg_48, fg_49, \
                         hg_150, hg_151, hg_152, hg_153, hg_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * fg_45[k]
                  + f_0 * hg_150[k];

        t_91[k] = -3.0 * fg_46[k]
                  + f_0 * hg_151[k];

        t_92[k] = -3.0 * fg_47[k]
                  + f_0 * hg_152[k];

        t_93[k] = -3.0 * fg_48[k]
                  + f_0 * hg_153[k];

        t_94[k] = -3.0 * fg_49[k]
                  + f_0 * hg_154[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, fg_50, fg_51, fg_52, fg_53, fg_54, \
                         hg_155, hg_156, hg_157, hg_158, hg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * fg_50[k]
                  + f_0 * hg_155[k];

        t_96[k] = -3.0 * fg_51[k]
                  + f_0 * hg_156[k];

        t_97[k] = -3.0 * fg_52[k]
                  + f_0 * hg_157[k];

        t_98[k] = -3.0 * fg_53[k]
                  + f_0 * hg_158[k];

        t_99[k] = -3.0 * fg_54[k]
                  + f_0 * hg_159[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, fg_55, fg_56, fg_57, fg_58, fg_59, \
                         hg_160, hg_161, hg_162, hg_163, hg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * fg_55[k]
                   + f_0 * hg_160[k];

        t_101[k] = -3.0 * fg_56[k]
                   + f_0 * hg_161[k];

        t_102[k] = -3.0 * fg_57[k]
                   + f_0 * hg_162[k];

        t_103[k] = -3.0 * fg_58[k]
                   + f_0 * hg_163[k];

        t_104[k] = -3.0 * fg_59[k]
                   + f_0 * hg_164[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, fg_60, fg_61, fg_62, fg_63, fg_64, \
                         hg_165, hg_166, hg_167, hg_168, hg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * fg_60[k]
                   + f_0 * hg_165[k];

        t_106[k] = -2.0 * fg_61[k]
                   + f_0 * hg_166[k];

        t_107[k] = -2.0 * fg_62[k]
                   + f_0 * hg_167[k];

        t_108[k] = -2.0 * fg_63[k]
                   + f_0 * hg_168[k];

        t_109[k] = -2.0 * fg_64[k]
                   + f_0 * hg_169[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, fg_65, fg_66, fg_67, fg_68, fg_69, \
                         hg_170, hg_171, hg_172, hg_173, hg_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * fg_65[k]
                   + f_0 * hg_170[k];

        t_111[k] = -2.0 * fg_66[k]
                   + f_0 * hg_171[k];

        t_112[k] = -2.0 * fg_67[k]
                   + f_0 * hg_172[k];

        t_113[k] = -2.0 * fg_68[k]
                   + f_0 * hg_173[k];

        t_114[k] = -2.0 * fg_69[k]
                   + f_0 * hg_174[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, fg_70, fg_71, fg_72, fg_73, fg_74, \
                         hg_175, hg_176, hg_177, hg_178, hg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * fg_70[k]
                   + f_0 * hg_175[k];

        t_116[k] = -2.0 * fg_71[k]
                   + f_0 * hg_176[k];

        t_117[k] = -2.0 * fg_72[k]
                   + f_0 * hg_177[k];

        t_118[k] = -2.0 * fg_73[k]
                   + f_0 * hg_178[k];

        t_119[k] = -2.0 * fg_74[k]
                   + f_0 * hg_179[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, fg_75, fg_76, fg_77, fg_78, fg_79, \
                         hg_180, hg_181, hg_182, hg_183, hg_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -fg_75[k]
                   + f_0 * hg_180[k];

        t_121[k] = -fg_76[k]
                   + f_0 * hg_181[k];

        t_122[k] = -fg_77[k]
                   + f_0 * hg_182[k];

        t_123[k] = -fg_78[k]
                   + f_0 * hg_183[k];

        t_124[k] = -fg_79[k]
                   + f_0 * hg_184[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, fg_80, fg_81, fg_82, fg_83, fg_84, \
                         hg_185, hg_186, hg_187, hg_188, hg_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -fg_80[k]
                   + f_0 * hg_185[k];

        t_126[k] = -fg_81[k]
                   + f_0 * hg_186[k];

        t_127[k] = -fg_82[k]
                   + f_0 * hg_187[k];

        t_128[k] = -fg_83[k]
                   + f_0 * hg_188[k];

        t_129[k] = -fg_84[k]
                   + f_0 * hg_189[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, fg_85, fg_86, fg_87, fg_88, fg_89, \
                         hg_190, hg_191, hg_192, hg_193, hg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -fg_85[k]
                   + f_0 * hg_190[k];

        t_131[k] = -fg_86[k]
                   + f_0 * hg_191[k];

        t_132[k] = -fg_87[k]
                   + f_0 * hg_192[k];

        t_133[k] = -fg_88[k]
                   + f_0 * hg_193[k];

        t_134[k] = -fg_89[k]
                   + f_0 * hg_194[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, t_140, t_141, t_142, hg_195, \
                         hg_196, hg_197, hg_198, hg_199, hg_200, hg_201, \
                         hg_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * hg_195[k];

        t_136[k] = f_0 * hg_196[k];

        t_137[k] = f_0 * hg_197[k];

        t_138[k] = f_0 * hg_198[k];

        t_139[k] = f_0 * hg_199[k];

        t_140[k] = f_0 * hg_200[k];

        t_141[k] = f_0 * hg_201[k];

        t_142[k] = f_0 * hg_202[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, t_148, t_149, hg_203, hg_204, \
                         hg_205, hg_206, hg_207, hg_208, hg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_0 * hg_203[k];

        t_144[k] = f_0 * hg_204[k];

        t_145[k] = f_0 * hg_205[k];

        t_146[k] = f_0 * hg_206[k];

        t_147[k] = f_0 * hg_207[k];

        t_148[k] = f_0 * hg_208[k];

        t_149[k] = f_0 * hg_209[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, fg_90, fg_91, fg_92, fg_93, fg_94, \
                         hg_225, hg_226, hg_227, hg_228, hg_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -4.0 * fg_90[k]
                   + f_0 * hg_225[k];

        t_151[k] = -4.0 * fg_91[k]
                   + f_0 * hg_226[k];

        t_152[k] = -4.0 * fg_92[k]
                   + f_0 * hg_227[k];

        t_153[k] = -4.0 * fg_93[k]
                   + f_0 * hg_228[k];

        t_154[k] = -4.0 * fg_94[k]
                   + f_0 * hg_229[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, fg_95, fg_96, fg_97, fg_98, fg_99, \
                         hg_230, hg_231, hg_232, hg_233, hg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -4.0 * fg_95[k]
                   + f_0 * hg_230[k];

        t_156[k] = -4.0 * fg_96[k]
                   + f_0 * hg_231[k];

        t_157[k] = -4.0 * fg_97[k]
                   + f_0 * hg_232[k];

        t_158[k] = -4.0 * fg_98[k]
                   + f_0 * hg_233[k];

        t_159[k] = -4.0 * fg_99[k]
                   + f_0 * hg_234[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, fg_100, fg_101, fg_102, fg_103, \
                         fg_104, hg_235, hg_236, hg_237, hg_238, \
                         hg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -4.0 * fg_100[k]
                   + f_0 * hg_235[k];

        t_161[k] = -4.0 * fg_101[k]
                   + f_0 * hg_236[k];

        t_162[k] = -4.0 * fg_102[k]
                   + f_0 * hg_237[k];

        t_163[k] = -4.0 * fg_103[k]
                   + f_0 * hg_238[k];

        t_164[k] = -4.0 * fg_104[k]
                   + f_0 * hg_239[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, fg_105, fg_106, fg_107, fg_108, \
                         fg_109, hg_240, hg_241, hg_242, hg_243, \
                         hg_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -3.0 * fg_105[k]
                   + f_0 * hg_240[k];

        t_166[k] = -3.0 * fg_106[k]
                   + f_0 * hg_241[k];

        t_167[k] = -3.0 * fg_107[k]
                   + f_0 * hg_242[k];

        t_168[k] = -3.0 * fg_108[k]
                   + f_0 * hg_243[k];

        t_169[k] = -3.0 * fg_109[k]
                   + f_0 * hg_244[k];
    }
}

static auto
compute_prim_geom_10_gg_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fg, const size_t hg,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_111 = buffer.data(fg + 111);
    const auto *fg_112 = buffer.data(fg + 112);
    const auto *fg_113 = buffer.data(fg + 113);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_120 = buffer.data(fg + 120);
    const auto *fg_121 = buffer.data(fg + 121);
    const auto *fg_122 = buffer.data(fg + 122);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_124 = buffer.data(fg + 124);
    const auto *fg_125 = buffer.data(fg + 125);
    const auto *fg_126 = buffer.data(fg + 126);
    const auto *fg_127 = buffer.data(fg + 127);
    const auto *fg_128 = buffer.data(fg + 128);
    const auto *fg_129 = buffer.data(fg + 129);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_136 = buffer.data(fg + 136);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_139 = buffer.data(fg + 139);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_141 = buffer.data(fg + 141);
    const auto *fg_142 = buffer.data(fg + 142);
    const auto *fg_143 = buffer.data(fg + 143);
    const auto *fg_144 = buffer.data(fg + 144);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_148 = buffer.data(fg + 148);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *hg_245 = buffer.data(hg + 245);
    const auto *hg_246 = buffer.data(hg + 246);
    const auto *hg_247 = buffer.data(hg + 247);
    const auto *hg_248 = buffer.data(hg + 248);
    const auto *hg_249 = buffer.data(hg + 249);
    const auto *hg_250 = buffer.data(hg + 250);
    const auto *hg_251 = buffer.data(hg + 251);
    const auto *hg_252 = buffer.data(hg + 252);
    const auto *hg_253 = buffer.data(hg + 253);
    const auto *hg_254 = buffer.data(hg + 254);
    const auto *hg_255 = buffer.data(hg + 255);
    const auto *hg_256 = buffer.data(hg + 256);
    const auto *hg_257 = buffer.data(hg + 257);
    const auto *hg_258 = buffer.data(hg + 258);
    const auto *hg_259 = buffer.data(hg + 259);
    const auto *hg_260 = buffer.data(hg + 260);
    const auto *hg_261 = buffer.data(hg + 261);
    const auto *hg_262 = buffer.data(hg + 262);
    const auto *hg_263 = buffer.data(hg + 263);
    const auto *hg_264 = buffer.data(hg + 264);
    const auto *hg_265 = buffer.data(hg + 265);
    const auto *hg_266 = buffer.data(hg + 266);
    const auto *hg_267 = buffer.data(hg + 267);
    const auto *hg_268 = buffer.data(hg + 268);
    const auto *hg_269 = buffer.data(hg + 269);
    const auto *hg_270 = buffer.data(hg + 270);
    const auto *hg_271 = buffer.data(hg + 271);
    const auto *hg_272 = buffer.data(hg + 272);
    const auto *hg_273 = buffer.data(hg + 273);
    const auto *hg_274 = buffer.data(hg + 274);
    const auto *hg_275 = buffer.data(hg + 275);
    const auto *hg_276 = buffer.data(hg + 276);
    const auto *hg_277 = buffer.data(hg + 277);
    const auto *hg_278 = buffer.data(hg + 278);
    const auto *hg_279 = buffer.data(hg + 279);
    const auto *hg_280 = buffer.data(hg + 280);
    const auto *hg_281 = buffer.data(hg + 281);
    const auto *hg_282 = buffer.data(hg + 282);
    const auto *hg_283 = buffer.data(hg + 283);
    const auto *hg_284 = buffer.data(hg + 284);
    const auto *hg_285 = buffer.data(hg + 285);
    const auto *hg_286 = buffer.data(hg + 286);
    const auto *hg_287 = buffer.data(hg + 287);
    const auto *hg_288 = buffer.data(hg + 288);
    const auto *hg_289 = buffer.data(hg + 289);
    const auto *hg_290 = buffer.data(hg + 290);
    const auto *hg_291 = buffer.data(hg + 291);
    const auto *hg_292 = buffer.data(hg + 292);
    const auto *hg_293 = buffer.data(hg + 293);
    const auto *hg_294 = buffer.data(hg + 294);
    const auto *hg_295 = buffer.data(hg + 295);
    const auto *hg_296 = buffer.data(hg + 296);
    const auto *hg_297 = buffer.data(hg + 297);
    const auto *hg_298 = buffer.data(hg + 298);
    const auto *hg_299 = buffer.data(hg + 299);

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, fg_110, fg_111, fg_112, fg_113, \
                         fg_114, hg_245, hg_246, hg_247, hg_248, \
                         hg_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * fg_110[k]
                   + f_0 * hg_245[k];

        t_171[k] = -3.0 * fg_111[k]
                   + f_0 * hg_246[k];

        t_172[k] = -3.0 * fg_112[k]
                   + f_0 * hg_247[k];

        t_173[k] = -3.0 * fg_113[k]
                   + f_0 * hg_248[k];

        t_174[k] = -3.0 * fg_114[k]
                   + f_0 * hg_249[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, fg_115, fg_116, fg_117, fg_118, \
                         fg_119, hg_250, hg_251, hg_252, hg_253, \
                         hg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * fg_115[k]
                   + f_0 * hg_250[k];

        t_176[k] = -3.0 * fg_116[k]
                   + f_0 * hg_251[k];

        t_177[k] = -3.0 * fg_117[k]
                   + f_0 * hg_252[k];

        t_178[k] = -3.0 * fg_118[k]
                   + f_0 * hg_253[k];

        t_179[k] = -3.0 * fg_119[k]
                   + f_0 * hg_254[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, fg_120, fg_121, fg_122, fg_123, \
                         fg_124, hg_255, hg_256, hg_257, hg_258, \
                         hg_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * fg_120[k]
                   + f_0 * hg_255[k];

        t_181[k] = -2.0 * fg_121[k]
                   + f_0 * hg_256[k];

        t_182[k] = -2.0 * fg_122[k]
                   + f_0 * hg_257[k];

        t_183[k] = -2.0 * fg_123[k]
                   + f_0 * hg_258[k];

        t_184[k] = -2.0 * fg_124[k]
                   + f_0 * hg_259[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, fg_125, fg_126, fg_127, fg_128, \
                         fg_129, hg_260, hg_261, hg_262, hg_263, \
                         hg_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * fg_125[k]
                   + f_0 * hg_260[k];

        t_186[k] = -2.0 * fg_126[k]
                   + f_0 * hg_261[k];

        t_187[k] = -2.0 * fg_127[k]
                   + f_0 * hg_262[k];

        t_188[k] = -2.0 * fg_128[k]
                   + f_0 * hg_263[k];

        t_189[k] = -2.0 * fg_129[k]
                   + f_0 * hg_264[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, fg_130, fg_131, fg_132, fg_133, \
                         fg_134, hg_265, hg_266, hg_267, hg_268, \
                         hg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * fg_130[k]
                   + f_0 * hg_265[k];

        t_191[k] = -2.0 * fg_131[k]
                   + f_0 * hg_266[k];

        t_192[k] = -2.0 * fg_132[k]
                   + f_0 * hg_267[k];

        t_193[k] = -2.0 * fg_133[k]
                   + f_0 * hg_268[k];

        t_194[k] = -2.0 * fg_134[k]
                   + f_0 * hg_269[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, fg_135, fg_136, fg_137, fg_138, \
                         fg_139, hg_270, hg_271, hg_272, hg_273, \
                         hg_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -fg_135[k]
                   + f_0 * hg_270[k];

        t_196[k] = -fg_136[k]
                   + f_0 * hg_271[k];

        t_197[k] = -fg_137[k]
                   + f_0 * hg_272[k];

        t_198[k] = -fg_138[k]
                   + f_0 * hg_273[k];

        t_199[k] = -fg_139[k]
                   + f_0 * hg_274[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, fg_140, fg_141, fg_142, fg_143, \
                         fg_144, hg_275, hg_276, hg_277, hg_278, \
                         hg_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -fg_140[k]
                   + f_0 * hg_275[k];

        t_201[k] = -fg_141[k]
                   + f_0 * hg_276[k];

        t_202[k] = -fg_142[k]
                   + f_0 * hg_277[k];

        t_203[k] = -fg_143[k]
                   + f_0 * hg_278[k];

        t_204[k] = -fg_144[k]
                   + f_0 * hg_279[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, fg_145, fg_146, fg_147, fg_148, \
                         fg_149, hg_280, hg_281, hg_282, hg_283, \
                         hg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -fg_145[k]
                   + f_0 * hg_280[k];

        t_206[k] = -fg_146[k]
                   + f_0 * hg_281[k];

        t_207[k] = -fg_147[k]
                   + f_0 * hg_282[k];

        t_208[k] = -fg_148[k]
                   + f_0 * hg_283[k];

        t_209[k] = -fg_149[k]
                   + f_0 * hg_284[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, t_217, hg_285, \
                         hg_286, hg_287, hg_288, hg_289, hg_290, hg_291, \
                         hg_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * hg_285[k];

        t_211[k] = f_0 * hg_286[k];

        t_212[k] = f_0 * hg_287[k];

        t_213[k] = f_0 * hg_288[k];

        t_214[k] = f_0 * hg_289[k];

        t_215[k] = f_0 * hg_290[k];

        t_216[k] = f_0 * hg_291[k];

        t_217[k] = f_0 * hg_292[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, t_224, hg_293, hg_294, \
                         hg_295, hg_296, hg_297, hg_298, hg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * hg_293[k];

        t_219[k] = f_0 * hg_294[k];

        t_220[k] = f_0 * hg_295[k];

        t_221[k] = f_0 * hg_296[k];

        t_222[k] = f_0 * hg_297[k];

        t_223[k] = f_0 * hg_298[k];

        t_224[k] = f_0 * hg_299[k];
    }
}

auto
compute_prim_geom_10_gg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t fg, const size_t hg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gg_electron_repulsion_1_piece0(buffer, target, fg, hg, ncols, alpha);

    compute_prim_geom_10_gg_electron_repulsion_1_piece1(buffer, target, fg, hg, ncols, alpha);
}

static auto
compute_prim_geom_10_gg_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fg, const size_t hg,
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

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_64 = buffer.data(fg + 64);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_74 = buffer.data(fg + 74);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_76 = buffer.data(fg + 76);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_81 = buffer.data(fg + 81);
    const auto *fg_82 = buffer.data(fg + 82);
    const auto *fg_83 = buffer.data(fg + 83);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_88 = buffer.data(fg + 88);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_92 = buffer.data(fg + 92);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_94 = buffer.data(fg + 94);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_96 = buffer.data(fg + 96);

    const auto *hg_30 = buffer.data(hg + 30);
    const auto *hg_31 = buffer.data(hg + 31);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_33 = buffer.data(hg + 33);
    const auto *hg_34 = buffer.data(hg + 34);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_36 = buffer.data(hg + 36);
    const auto *hg_37 = buffer.data(hg + 37);
    const auto *hg_38 = buffer.data(hg + 38);
    const auto *hg_39 = buffer.data(hg + 39);
    const auto *hg_40 = buffer.data(hg + 40);
    const auto *hg_41 = buffer.data(hg + 41);
    const auto *hg_42 = buffer.data(hg + 42);
    const auto *hg_43 = buffer.data(hg + 43);
    const auto *hg_44 = buffer.data(hg + 44);
    const auto *hg_60 = buffer.data(hg + 60);
    const auto *hg_61 = buffer.data(hg + 61);
    const auto *hg_62 = buffer.data(hg + 62);
    const auto *hg_63 = buffer.data(hg + 63);
    const auto *hg_64 = buffer.data(hg + 64);
    const auto *hg_65 = buffer.data(hg + 65);
    const auto *hg_66 = buffer.data(hg + 66);
    const auto *hg_67 = buffer.data(hg + 67);
    const auto *hg_68 = buffer.data(hg + 68);
    const auto *hg_69 = buffer.data(hg + 69);
    const auto *hg_70 = buffer.data(hg + 70);
    const auto *hg_71 = buffer.data(hg + 71);
    const auto *hg_72 = buffer.data(hg + 72);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_74 = buffer.data(hg + 74);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_79 = buffer.data(hg + 79);
    const auto *hg_80 = buffer.data(hg + 80);
    const auto *hg_81 = buffer.data(hg + 81);
    const auto *hg_82 = buffer.data(hg + 82);
    const auto *hg_83 = buffer.data(hg + 83);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_85 = buffer.data(hg + 85);
    const auto *hg_86 = buffer.data(hg + 86);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_88 = buffer.data(hg + 88);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_106 = buffer.data(hg + 106);
    const auto *hg_107 = buffer.data(hg + 107);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_109 = buffer.data(hg + 109);
    const auto *hg_110 = buffer.data(hg + 110);
    const auto *hg_111 = buffer.data(hg + 111);
    const auto *hg_112 = buffer.data(hg + 112);
    const auto *hg_113 = buffer.data(hg + 113);
    const auto *hg_114 = buffer.data(hg + 114);
    const auto *hg_115 = buffer.data(hg + 115);
    const auto *hg_116 = buffer.data(hg + 116);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_118 = buffer.data(hg + 118);
    const auto *hg_119 = buffer.data(hg + 119);
    const auto *hg_120 = buffer.data(hg + 120);
    const auto *hg_121 = buffer.data(hg + 121);
    const auto *hg_122 = buffer.data(hg + 122);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_124 = buffer.data(hg + 124);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_126 = buffer.data(hg + 126);
    const auto *hg_127 = buffer.data(hg + 127);
    const auto *hg_128 = buffer.data(hg + 128);
    const auto *hg_129 = buffer.data(hg + 129);
    const auto *hg_130 = buffer.data(hg + 130);
    const auto *hg_131 = buffer.data(hg + 131);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_133 = buffer.data(hg + 133);
    const auto *hg_134 = buffer.data(hg + 134);
    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_136 = buffer.data(hg + 136);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_138 = buffer.data(hg + 138);
    const auto *hg_139 = buffer.data(hg + 139);
    const auto *hg_140 = buffer.data(hg + 140);
    const auto *hg_141 = buffer.data(hg + 141);
    const auto *hg_142 = buffer.data(hg + 142);
    const auto *hg_143 = buffer.data(hg + 143);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_145 = buffer.data(hg + 145);
    const auto *hg_146 = buffer.data(hg + 146);
    const auto *hg_147 = buffer.data(hg + 147);
    const auto *hg_148 = buffer.data(hg + 148);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_165 = buffer.data(hg + 165);
    const auto *hg_166 = buffer.data(hg + 166);
    const auto *hg_167 = buffer.data(hg + 167);
    const auto *hg_168 = buffer.data(hg + 168);
    const auto *hg_169 = buffer.data(hg + 169);
    const auto *hg_170 = buffer.data(hg + 170);
    const auto *hg_171 = buffer.data(hg + 171);
    const auto *hg_172 = buffer.data(hg + 172);
    const auto *hg_173 = buffer.data(hg + 173);
    const auto *hg_174 = buffer.data(hg + 174);
    const auto *hg_175 = buffer.data(hg + 175);
    const auto *hg_176 = buffer.data(hg + 176);
    const auto *hg_177 = buffer.data(hg + 177);
    const auto *hg_178 = buffer.data(hg + 178);
    const auto *hg_179 = buffer.data(hg + 179);
    const auto *hg_180 = buffer.data(hg + 180);
    const auto *hg_181 = buffer.data(hg + 181);
    const auto *hg_182 = buffer.data(hg + 182);
    const auto *hg_183 = buffer.data(hg + 183);
    const auto *hg_184 = buffer.data(hg + 184);
    const auto *hg_185 = buffer.data(hg + 185);
    const auto *hg_186 = buffer.data(hg + 186);
    const auto *hg_187 = buffer.data(hg + 187);
    const auto *hg_188 = buffer.data(hg + 188);
    const auto *hg_189 = buffer.data(hg + 189);
    const auto *hg_190 = buffer.data(hg + 190);
    const auto *hg_191 = buffer.data(hg + 191);
    const auto *hg_192 = buffer.data(hg + 192);
    const auto *hg_193 = buffer.data(hg + 193);
    const auto *hg_194 = buffer.data(hg + 194);
    const auto *hg_195 = buffer.data(hg + 195);
    const auto *hg_196 = buffer.data(hg + 196);
    const auto *hg_197 = buffer.data(hg + 197);
    const auto *hg_198 = buffer.data(hg + 198);
    const auto *hg_199 = buffer.data(hg + 199);
    const auto *hg_200 = buffer.data(hg + 200);
    const auto *hg_201 = buffer.data(hg + 201);
    const auto *hg_202 = buffer.data(hg + 202);
    const auto *hg_203 = buffer.data(hg + 203);
    const auto *hg_204 = buffer.data(hg + 204);
    const auto *hg_205 = buffer.data(hg + 205);
    const auto *hg_206 = buffer.data(hg + 206);
    const auto *hg_207 = buffer.data(hg + 207);
    const auto *hg_208 = buffer.data(hg + 208);
    const auto *hg_209 = buffer.data(hg + 209);
    const auto *hg_210 = buffer.data(hg + 210);
    const auto *hg_211 = buffer.data(hg + 211);
    const auto *hg_212 = buffer.data(hg + 212);
    const auto *hg_213 = buffer.data(hg + 213);
    const auto *hg_214 = buffer.data(hg + 214);
    const auto *hg_215 = buffer.data(hg + 215);
    const auto *hg_216 = buffer.data(hg + 216);
    const auto *hg_217 = buffer.data(hg + 217);
    const auto *hg_218 = buffer.data(hg + 218);
    const auto *hg_219 = buffer.data(hg + 219);
    const auto *hg_220 = buffer.data(hg + 220);
    const auto *hg_221 = buffer.data(hg + 221);
    const auto *hg_222 = buffer.data(hg + 222);
    const auto *hg_223 = buffer.data(hg + 223);
    const auto *hg_224 = buffer.data(hg + 224);
    const auto *hg_240 = buffer.data(hg + 240);
    const auto *hg_241 = buffer.data(hg + 241);
    const auto *hg_242 = buffer.data(hg + 242);
    const auto *hg_243 = buffer.data(hg + 243);
    const auto *hg_244 = buffer.data(hg + 244);
    const auto *hg_245 = buffer.data(hg + 245);
    const auto *hg_246 = buffer.data(hg + 246);
    const auto *hg_247 = buffer.data(hg + 247);
    const auto *hg_248 = buffer.data(hg + 248);
    const auto *hg_249 = buffer.data(hg + 249);
    const auto *hg_250 = buffer.data(hg + 250);
    const auto *hg_251 = buffer.data(hg + 251);
    const auto *hg_252 = buffer.data(hg + 252);
    const auto *hg_253 = buffer.data(hg + 253);
    const auto *hg_254 = buffer.data(hg + 254);
    const auto *hg_255 = buffer.data(hg + 255);
    const auto *hg_256 = buffer.data(hg + 256);
    const auto *hg_257 = buffer.data(hg + 257);
    const auto *hg_258 = buffer.data(hg + 258);
    const auto *hg_259 = buffer.data(hg + 259);
    const auto *hg_260 = buffer.data(hg + 260);
    const auto *hg_261 = buffer.data(hg + 261);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hg_30, hg_31, hg_32, hg_33, \
                         hg_34, hg_35, hg_36, hg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hg_30[k];

        t_1[k] = f_0 * hg_31[k];

        t_2[k] = f_0 * hg_32[k];

        t_3[k] = f_0 * hg_33[k];

        t_4[k] = f_0 * hg_34[k];

        t_5[k] = f_0 * hg_35[k];

        t_6[k] = f_0 * hg_36[k];

        t_7[k] = f_0 * hg_37[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, hg_38, hg_39, hg_40, \
                         hg_41, hg_42, hg_43, hg_44, hg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hg_38[k];

        t_9[k] = f_0 * hg_39[k];

        t_10[k] = f_0 * hg_40[k];

        t_11[k] = f_0 * hg_41[k];

        t_12[k] = f_0 * hg_42[k];

        t_13[k] = f_0 * hg_43[k];

        t_14[k] = f_0 * hg_44[k];

        t_15[k] = f_0 * hg_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, hg_61, hg_62, hg_63, \
                         hg_64, hg_65, hg_66, hg_67, hg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * hg_61[k];

        t_17[k] = f_0 * hg_62[k];

        t_18[k] = f_0 * hg_63[k];

        t_19[k] = f_0 * hg_64[k];

        t_20[k] = f_0 * hg_65[k];

        t_21[k] = f_0 * hg_66[k];

        t_22[k] = f_0 * hg_67[k];

        t_23[k] = f_0 * hg_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, fg_0, hg_69, hg_70, hg_71, \
                         hg_72, hg_73, hg_74, hg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * hg_69[k];

        t_25[k] = f_0 * hg_70[k];

        t_26[k] = f_0 * hg_71[k];

        t_27[k] = f_0 * hg_72[k];

        t_28[k] = f_0 * hg_73[k];

        t_29[k] = f_0 * hg_74[k];

        t_30[k] = -fg_0[k]
                  + f_0 * hg_75[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, fg_1, fg_2, fg_3, fg_4, fg_5, hg_76, \
                         hg_77, hg_78, hg_79, hg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -fg_1[k]
                  + f_0 * hg_76[k];

        t_32[k] = -fg_2[k]
                  + f_0 * hg_77[k];

        t_33[k] = -fg_3[k]
                  + f_0 * hg_78[k];

        t_34[k] = -fg_4[k]
                  + f_0 * hg_79[k];

        t_35[k] = -fg_5[k]
                  + f_0 * hg_80[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, fg_6, fg_7, fg_8, fg_9, fg_10, hg_81, \
                         hg_82, hg_83, hg_84, hg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -fg_6[k]
                  + f_0 * hg_81[k];

        t_37[k] = -fg_7[k]
                  + f_0 * hg_82[k];

        t_38[k] = -fg_8[k]
                  + f_0 * hg_83[k];

        t_39[k] = -fg_9[k]
                  + f_0 * hg_84[k];

        t_40[k] = -fg_10[k]
                  + f_0 * hg_85[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, t_46, fg_11, fg_12, fg_13, fg_14, \
                         hg_86, hg_87, hg_88, hg_89, hg_105, hg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -fg_11[k]
                  + f_0 * hg_86[k];

        t_42[k] = -fg_12[k]
                  + f_0 * hg_87[k];

        t_43[k] = -fg_13[k]
                  + f_0 * hg_88[k];

        t_44[k] = -fg_14[k]
                  + f_0 * hg_89[k];

        t_45[k] = f_0 * hg_105[k];

        t_46[k] = f_0 * hg_106[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, t_53, t_54, hg_107, hg_108, \
                         hg_109, hg_110, hg_111, hg_112, hg_113, \
                         hg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * hg_107[k];

        t_48[k] = f_0 * hg_108[k];

        t_49[k] = f_0 * hg_109[k];

        t_50[k] = f_0 * hg_110[k];

        t_51[k] = f_0 * hg_111[k];

        t_52[k] = f_0 * hg_112[k];

        t_53[k] = f_0 * hg_113[k];

        t_54[k] = f_0 * hg_114[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, fg_15, fg_16, hg_115, \
                         hg_116, hg_117, hg_118, hg_119, hg_120, \
                         hg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * hg_115[k];

        t_56[k] = f_0 * hg_116[k];

        t_57[k] = f_0 * hg_117[k];

        t_58[k] = f_0 * hg_118[k];

        t_59[k] = f_0 * hg_119[k];

        t_60[k] = -fg_15[k]
                  + f_0 * hg_120[k];

        t_61[k] = -fg_16[k]
                  + f_0 * hg_121[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, fg_17, fg_18, fg_19, fg_20, fg_21, \
                         hg_122, hg_123, hg_124, hg_125, hg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -fg_17[k]
                  + f_0 * hg_122[k];

        t_63[k] = -fg_18[k]
                  + f_0 * hg_123[k];

        t_64[k] = -fg_19[k]
                  + f_0 * hg_124[k];

        t_65[k] = -fg_20[k]
                  + f_0 * hg_125[k];

        t_66[k] = -fg_21[k]
                  + f_0 * hg_126[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, fg_22, fg_23, fg_24, fg_25, fg_26, \
                         hg_127, hg_128, hg_129, hg_130, hg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -fg_22[k]
                  + f_0 * hg_127[k];

        t_68[k] = -fg_23[k]
                  + f_0 * hg_128[k];

        t_69[k] = -fg_24[k]
                  + f_0 * hg_129[k];

        t_70[k] = -fg_25[k]
                  + f_0 * hg_130[k];

        t_71[k] = -fg_26[k]
                  + f_0 * hg_131[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, fg_27, fg_28, fg_29, fg_30, fg_31, \
                         hg_132, hg_133, hg_134, hg_135, hg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -fg_27[k]
                  + f_0 * hg_132[k];

        t_73[k] = -fg_28[k]
                  + f_0 * hg_133[k];

        t_74[k] = -fg_29[k]
                  + f_0 * hg_134[k];

        t_75[k] = -2.0 * fg_30[k]
                  + f_0 * hg_135[k];

        t_76[k] = -2.0 * fg_31[k]
                  + f_0 * hg_136[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, fg_32, fg_33, fg_34, fg_35, fg_36, \
                         hg_137, hg_138, hg_139, hg_140, hg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -2.0 * fg_32[k]
                  + f_0 * hg_137[k];

        t_78[k] = -2.0 * fg_33[k]
                  + f_0 * hg_138[k];

        t_79[k] = -2.0 * fg_34[k]
                  + f_0 * hg_139[k];

        t_80[k] = -2.0 * fg_35[k]
                  + f_0 * hg_140[k];

        t_81[k] = -2.0 * fg_36[k]
                  + f_0 * hg_141[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, fg_37, fg_38, fg_39, fg_40, fg_41, \
                         hg_142, hg_143, hg_144, hg_145, hg_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -2.0 * fg_37[k]
                  + f_0 * hg_142[k];

        t_83[k] = -2.0 * fg_38[k]
                  + f_0 * hg_143[k];

        t_84[k] = -2.0 * fg_39[k]
                  + f_0 * hg_144[k];

        t_85[k] = -2.0 * fg_40[k]
                  + f_0 * hg_145[k];

        t_86[k] = -2.0 * fg_41[k]
                  + f_0 * hg_146[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, fg_42, fg_43, fg_44, hg_147, \
                         hg_148, hg_149, hg_165, hg_166, hg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -2.0 * fg_42[k]
                  + f_0 * hg_147[k];

        t_88[k] = -2.0 * fg_43[k]
                  + f_0 * hg_148[k];

        t_89[k] = -2.0 * fg_44[k]
                  + f_0 * hg_149[k];

        t_90[k] = f_0 * hg_165[k];

        t_91[k] = f_0 * hg_166[k];

        t_92[k] = f_0 * hg_167[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, hg_168, hg_169, \
                         hg_170, hg_171, hg_172, hg_173, hg_174, \
                         hg_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * hg_168[k];

        t_94[k] = f_0 * hg_169[k];

        t_95[k] = f_0 * hg_170[k];

        t_96[k] = f_0 * hg_171[k];

        t_97[k] = f_0 * hg_172[k];

        t_98[k] = f_0 * hg_173[k];

        t_99[k] = f_0 * hg_174[k];

        t_100[k] = f_0 * hg_175[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, fg_45, fg_46, hg_176, \
                         hg_177, hg_178, hg_179, hg_180, hg_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * hg_176[k];

        t_102[k] = f_0 * hg_177[k];

        t_103[k] = f_0 * hg_178[k];

        t_104[k] = f_0 * hg_179[k];

        t_105[k] = -fg_45[k]
                   + f_0 * hg_180[k];

        t_106[k] = -fg_46[k]
                   + f_0 * hg_181[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, fg_47, fg_48, fg_49, fg_50, fg_51, \
                         hg_182, hg_183, hg_184, hg_185, hg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -fg_47[k]
                   + f_0 * hg_182[k];

        t_108[k] = -fg_48[k]
                   + f_0 * hg_183[k];

        t_109[k] = -fg_49[k]
                   + f_0 * hg_184[k];

        t_110[k] = -fg_50[k]
                   + f_0 * hg_185[k];

        t_111[k] = -fg_51[k]
                   + f_0 * hg_186[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, fg_52, fg_53, fg_54, fg_55, fg_56, \
                         hg_187, hg_188, hg_189, hg_190, hg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -fg_52[k]
                   + f_0 * hg_187[k];

        t_113[k] = -fg_53[k]
                   + f_0 * hg_188[k];

        t_114[k] = -fg_54[k]
                   + f_0 * hg_189[k];

        t_115[k] = -fg_55[k]
                   + f_0 * hg_190[k];

        t_116[k] = -fg_56[k]
                   + f_0 * hg_191[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, fg_57, fg_58, fg_59, fg_60, fg_61, \
                         hg_192, hg_193, hg_194, hg_195, hg_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -fg_57[k]
                   + f_0 * hg_192[k];

        t_118[k] = -fg_58[k]
                   + f_0 * hg_193[k];

        t_119[k] = -fg_59[k]
                   + f_0 * hg_194[k];

        t_120[k] = -2.0 * fg_60[k]
                   + f_0 * hg_195[k];

        t_121[k] = -2.0 * fg_61[k]
                   + f_0 * hg_196[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, fg_62, fg_63, fg_64, fg_65, fg_66, \
                         hg_197, hg_198, hg_199, hg_200, hg_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * fg_62[k]
                   + f_0 * hg_197[k];

        t_123[k] = -2.0 * fg_63[k]
                   + f_0 * hg_198[k];

        t_124[k] = -2.0 * fg_64[k]
                   + f_0 * hg_199[k];

        t_125[k] = -2.0 * fg_65[k]
                   + f_0 * hg_200[k];

        t_126[k] = -2.0 * fg_66[k]
                   + f_0 * hg_201[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, fg_67, fg_68, fg_69, fg_70, fg_71, \
                         hg_202, hg_203, hg_204, hg_205, hg_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * fg_67[k]
                   + f_0 * hg_202[k];

        t_128[k] = -2.0 * fg_68[k]
                   + f_0 * hg_203[k];

        t_129[k] = -2.0 * fg_69[k]
                   + f_0 * hg_204[k];

        t_130[k] = -2.0 * fg_70[k]
                   + f_0 * hg_205[k];

        t_131[k] = -2.0 * fg_71[k]
                   + f_0 * hg_206[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, fg_72, fg_73, fg_74, fg_75, fg_76, \
                         hg_207, hg_208, hg_209, hg_210, hg_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * fg_72[k]
                   + f_0 * hg_207[k];

        t_133[k] = -2.0 * fg_73[k]
                   + f_0 * hg_208[k];

        t_134[k] = -2.0 * fg_74[k]
                   + f_0 * hg_209[k];

        t_135[k] = -3.0 * fg_75[k]
                   + f_0 * hg_210[k];

        t_136[k] = -3.0 * fg_76[k]
                   + f_0 * hg_211[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, fg_77, fg_78, fg_79, fg_80, fg_81, \
                         hg_212, hg_213, hg_214, hg_215, hg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -3.0 * fg_77[k]
                   + f_0 * hg_212[k];

        t_138[k] = -3.0 * fg_78[k]
                   + f_0 * hg_213[k];

        t_139[k] = -3.0 * fg_79[k]
                   + f_0 * hg_214[k];

        t_140[k] = -3.0 * fg_80[k]
                   + f_0 * hg_215[k];

        t_141[k] = -3.0 * fg_81[k]
                   + f_0 * hg_216[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, fg_82, fg_83, fg_84, fg_85, fg_86, \
                         hg_217, hg_218, hg_219, hg_220, hg_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -3.0 * fg_82[k]
                   + f_0 * hg_217[k];

        t_143[k] = -3.0 * fg_83[k]
                   + f_0 * hg_218[k];

        t_144[k] = -3.0 * fg_84[k]
                   + f_0 * hg_219[k];

        t_145[k] = -3.0 * fg_85[k]
                   + f_0 * hg_220[k];

        t_146[k] = -3.0 * fg_86[k]
                   + f_0 * hg_221[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, t_152, fg_87, fg_88, fg_89, \
                         hg_222, hg_223, hg_224, hg_240, hg_241, \
                         hg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -3.0 * fg_87[k]
                   + f_0 * hg_222[k];

        t_148[k] = -3.0 * fg_88[k]
                   + f_0 * hg_223[k];

        t_149[k] = -3.0 * fg_89[k]
                   + f_0 * hg_224[k];

        t_150[k] = f_0 * hg_240[k];

        t_151[k] = f_0 * hg_241[k];

        t_152[k] = f_0 * hg_242[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, t_158, t_159, t_160, hg_243, \
                         hg_244, hg_245, hg_246, hg_247, hg_248, hg_249, \
                         hg_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_0 * hg_243[k];

        t_154[k] = f_0 * hg_244[k];

        t_155[k] = f_0 * hg_245[k];

        t_156[k] = f_0 * hg_246[k];

        t_157[k] = f_0 * hg_247[k];

        t_158[k] = f_0 * hg_248[k];

        t_159[k] = f_0 * hg_249[k];

        t_160[k] = f_0 * hg_250[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, fg_90, fg_91, hg_251, \
                         hg_252, hg_253, hg_254, hg_255, hg_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_0 * hg_251[k];

        t_162[k] = f_0 * hg_252[k];

        t_163[k] = f_0 * hg_253[k];

        t_164[k] = f_0 * hg_254[k];

        t_165[k] = -fg_90[k]
                   + f_0 * hg_255[k];

        t_166[k] = -fg_91[k]
                   + f_0 * hg_256[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, fg_92, fg_93, fg_94, fg_95, fg_96, \
                         hg_257, hg_258, hg_259, hg_260, hg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -fg_92[k]
                   + f_0 * hg_257[k];

        t_168[k] = -fg_93[k]
                   + f_0 * hg_258[k];

        t_169[k] = -fg_94[k]
                   + f_0 * hg_259[k];

        t_170[k] = -fg_95[k]
                   + f_0 * hg_260[k];

        t_171[k] = -fg_96[k]
                   + f_0 * hg_261[k];
    }
}

static auto
compute_prim_geom_10_gg_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fg, const size_t hg,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *fg_97 = buffer.data(fg + 97);
    const auto *fg_98 = buffer.data(fg + 98);
    const auto *fg_99 = buffer.data(fg + 99);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_101 = buffer.data(fg + 101);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_105 = buffer.data(fg + 105);
    const auto *fg_106 = buffer.data(fg + 106);
    const auto *fg_107 = buffer.data(fg + 107);
    const auto *fg_108 = buffer.data(fg + 108);
    const auto *fg_109 = buffer.data(fg + 109);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_111 = buffer.data(fg + 111);
    const auto *fg_112 = buffer.data(fg + 112);
    const auto *fg_113 = buffer.data(fg + 113);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_120 = buffer.data(fg + 120);
    const auto *fg_121 = buffer.data(fg + 121);
    const auto *fg_122 = buffer.data(fg + 122);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_124 = buffer.data(fg + 124);
    const auto *fg_125 = buffer.data(fg + 125);
    const auto *fg_126 = buffer.data(fg + 126);
    const auto *fg_127 = buffer.data(fg + 127);
    const auto *fg_128 = buffer.data(fg + 128);
    const auto *fg_129 = buffer.data(fg + 129);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_136 = buffer.data(fg + 136);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_139 = buffer.data(fg + 139);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_141 = buffer.data(fg + 141);
    const auto *fg_142 = buffer.data(fg + 142);
    const auto *fg_143 = buffer.data(fg + 143);
    const auto *fg_144 = buffer.data(fg + 144);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_148 = buffer.data(fg + 148);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *hg_262 = buffer.data(hg + 262);
    const auto *hg_263 = buffer.data(hg + 263);
    const auto *hg_264 = buffer.data(hg + 264);
    const auto *hg_265 = buffer.data(hg + 265);
    const auto *hg_266 = buffer.data(hg + 266);
    const auto *hg_267 = buffer.data(hg + 267);
    const auto *hg_268 = buffer.data(hg + 268);
    const auto *hg_269 = buffer.data(hg + 269);
    const auto *hg_270 = buffer.data(hg + 270);
    const auto *hg_271 = buffer.data(hg + 271);
    const auto *hg_272 = buffer.data(hg + 272);
    const auto *hg_273 = buffer.data(hg + 273);
    const auto *hg_274 = buffer.data(hg + 274);
    const auto *hg_275 = buffer.data(hg + 275);
    const auto *hg_276 = buffer.data(hg + 276);
    const auto *hg_277 = buffer.data(hg + 277);
    const auto *hg_278 = buffer.data(hg + 278);
    const auto *hg_279 = buffer.data(hg + 279);
    const auto *hg_280 = buffer.data(hg + 280);
    const auto *hg_281 = buffer.data(hg + 281);
    const auto *hg_282 = buffer.data(hg + 282);
    const auto *hg_283 = buffer.data(hg + 283);
    const auto *hg_284 = buffer.data(hg + 284);
    const auto *hg_285 = buffer.data(hg + 285);
    const auto *hg_286 = buffer.data(hg + 286);
    const auto *hg_287 = buffer.data(hg + 287);
    const auto *hg_288 = buffer.data(hg + 288);
    const auto *hg_289 = buffer.data(hg + 289);
    const auto *hg_290 = buffer.data(hg + 290);
    const auto *hg_291 = buffer.data(hg + 291);
    const auto *hg_292 = buffer.data(hg + 292);
    const auto *hg_293 = buffer.data(hg + 293);
    const auto *hg_294 = buffer.data(hg + 294);
    const auto *hg_295 = buffer.data(hg + 295);
    const auto *hg_296 = buffer.data(hg + 296);
    const auto *hg_297 = buffer.data(hg + 297);
    const auto *hg_298 = buffer.data(hg + 298);
    const auto *hg_299 = buffer.data(hg + 299);
    const auto *hg_300 = buffer.data(hg + 300);
    const auto *hg_301 = buffer.data(hg + 301);
    const auto *hg_302 = buffer.data(hg + 302);
    const auto *hg_303 = buffer.data(hg + 303);
    const auto *hg_304 = buffer.data(hg + 304);
    const auto *hg_305 = buffer.data(hg + 305);
    const auto *hg_306 = buffer.data(hg + 306);
    const auto *hg_307 = buffer.data(hg + 307);
    const auto *hg_308 = buffer.data(hg + 308);
    const auto *hg_309 = buffer.data(hg + 309);
    const auto *hg_310 = buffer.data(hg + 310);
    const auto *hg_311 = buffer.data(hg + 311);
    const auto *hg_312 = buffer.data(hg + 312);
    const auto *hg_313 = buffer.data(hg + 313);
    const auto *hg_314 = buffer.data(hg + 314);

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, fg_97, fg_98, fg_99, fg_100, \
                         fg_101, hg_262, hg_263, hg_264, hg_265, \
                         hg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -fg_97[k]
                   + f_0 * hg_262[k];

        t_173[k] = -fg_98[k]
                   + f_0 * hg_263[k];

        t_174[k] = -fg_99[k]
                   + f_0 * hg_264[k];

        t_175[k] = -fg_100[k]
                   + f_0 * hg_265[k];

        t_176[k] = -fg_101[k]
                   + f_0 * hg_266[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, fg_102, fg_103, fg_104, fg_105, \
                         fg_106, hg_267, hg_268, hg_269, hg_270, \
                         hg_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -fg_102[k]
                   + f_0 * hg_267[k];

        t_178[k] = -fg_103[k]
                   + f_0 * hg_268[k];

        t_179[k] = -fg_104[k]
                   + f_0 * hg_269[k];

        t_180[k] = -2.0 * fg_105[k]
                   + f_0 * hg_270[k];

        t_181[k] = -2.0 * fg_106[k]
                   + f_0 * hg_271[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, fg_107, fg_108, fg_109, fg_110, \
                         fg_111, hg_272, hg_273, hg_274, hg_275, \
                         hg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * fg_107[k]
                   + f_0 * hg_272[k];

        t_183[k] = -2.0 * fg_108[k]
                   + f_0 * hg_273[k];

        t_184[k] = -2.0 * fg_109[k]
                   + f_0 * hg_274[k];

        t_185[k] = -2.0 * fg_110[k]
                   + f_0 * hg_275[k];

        t_186[k] = -2.0 * fg_111[k]
                   + f_0 * hg_276[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, fg_112, fg_113, fg_114, fg_115, \
                         fg_116, hg_277, hg_278, hg_279, hg_280, \
                         hg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * fg_112[k]
                   + f_0 * hg_277[k];

        t_188[k] = -2.0 * fg_113[k]
                   + f_0 * hg_278[k];

        t_189[k] = -2.0 * fg_114[k]
                   + f_0 * hg_279[k];

        t_190[k] = -2.0 * fg_115[k]
                   + f_0 * hg_280[k];

        t_191[k] = -2.0 * fg_116[k]
                   + f_0 * hg_281[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, fg_117, fg_118, fg_119, fg_120, \
                         fg_121, hg_282, hg_283, hg_284, hg_285, \
                         hg_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -2.0 * fg_117[k]
                   + f_0 * hg_282[k];

        t_193[k] = -2.0 * fg_118[k]
                   + f_0 * hg_283[k];

        t_194[k] = -2.0 * fg_119[k]
                   + f_0 * hg_284[k];

        t_195[k] = -3.0 * fg_120[k]
                   + f_0 * hg_285[k];

        t_196[k] = -3.0 * fg_121[k]
                   + f_0 * hg_286[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, fg_122, fg_123, fg_124, fg_125, \
                         fg_126, hg_287, hg_288, hg_289, hg_290, \
                         hg_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -3.0 * fg_122[k]
                   + f_0 * hg_287[k];

        t_198[k] = -3.0 * fg_123[k]
                   + f_0 * hg_288[k];

        t_199[k] = -3.0 * fg_124[k]
                   + f_0 * hg_289[k];

        t_200[k] = -3.0 * fg_125[k]
                   + f_0 * hg_290[k];

        t_201[k] = -3.0 * fg_126[k]
                   + f_0 * hg_291[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, fg_127, fg_128, fg_129, fg_130, \
                         fg_131, hg_292, hg_293, hg_294, hg_295, \
                         hg_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -3.0 * fg_127[k]
                   + f_0 * hg_292[k];

        t_203[k] = -3.0 * fg_128[k]
                   + f_0 * hg_293[k];

        t_204[k] = -3.0 * fg_129[k]
                   + f_0 * hg_294[k];

        t_205[k] = -3.0 * fg_130[k]
                   + f_0 * hg_295[k];

        t_206[k] = -3.0 * fg_131[k]
                   + f_0 * hg_296[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, fg_132, fg_133, fg_134, fg_135, \
                         fg_136, hg_297, hg_298, hg_299, hg_300, \
                         hg_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -3.0 * fg_132[k]
                   + f_0 * hg_297[k];

        t_208[k] = -3.0 * fg_133[k]
                   + f_0 * hg_298[k];

        t_209[k] = -3.0 * fg_134[k]
                   + f_0 * hg_299[k];

        t_210[k] = -4.0 * fg_135[k]
                   + f_0 * hg_300[k];

        t_211[k] = -4.0 * fg_136[k]
                   + f_0 * hg_301[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, fg_137, fg_138, fg_139, fg_140, \
                         fg_141, hg_302, hg_303, hg_304, hg_305, \
                         hg_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = -4.0 * fg_137[k]
                   + f_0 * hg_302[k];

        t_213[k] = -4.0 * fg_138[k]
                   + f_0 * hg_303[k];

        t_214[k] = -4.0 * fg_139[k]
                   + f_0 * hg_304[k];

        t_215[k] = -4.0 * fg_140[k]
                   + f_0 * hg_305[k];

        t_216[k] = -4.0 * fg_141[k]
                   + f_0 * hg_306[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, fg_142, fg_143, fg_144, fg_145, \
                         fg_146, hg_307, hg_308, hg_309, hg_310, \
                         hg_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = -4.0 * fg_142[k]
                   + f_0 * hg_307[k];

        t_218[k] = -4.0 * fg_143[k]
                   + f_0 * hg_308[k];

        t_219[k] = -4.0 * fg_144[k]
                   + f_0 * hg_309[k];

        t_220[k] = -4.0 * fg_145[k]
                   + f_0 * hg_310[k];

        t_221[k] = -4.0 * fg_146[k]
                   + f_0 * hg_311[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, fg_147, fg_148, fg_149, hg_312, hg_313, \
                         hg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -4.0 * fg_147[k]
                   + f_0 * hg_312[k];

        t_223[k] = -4.0 * fg_148[k]
                   + f_0 * hg_313[k];

        t_224[k] = -4.0 * fg_149[k]
                   + f_0 * hg_314[k];
    }
}

auto
compute_prim_geom_10_gg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t fg, const size_t hg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gg_electron_repulsion_2_piece0(buffer, target, fg, hg, ncols, alpha);

    compute_prim_geom_10_gg_electron_repulsion_2_piece1(buffer, target, fg, hg, ncols, alpha);
}

}  // namespace simdt2ceri
