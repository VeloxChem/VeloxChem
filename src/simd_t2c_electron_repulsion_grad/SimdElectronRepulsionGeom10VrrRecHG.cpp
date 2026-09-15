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


#include "SimdElectronRepulsionGeom10VrrRecHG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_hg_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gg, const size_t ig,
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

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_94 = buffer.data(gg + 94);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_97 = buffer.data(gg + 97);
    const auto *gg_98 = buffer.data(gg + 98);
    const auto *gg_99 = buffer.data(gg + 99);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_101 = buffer.data(gg + 101);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_106 = buffer.data(gg + 106);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_109 = buffer.data(gg + 109);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_113 = buffer.data(gg + 113);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_124 = buffer.data(gg + 124);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_128 = buffer.data(gg + 128);
    const auto *gg_129 = buffer.data(gg + 129);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_134 = buffer.data(gg + 134);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_136 = buffer.data(gg + 136);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_138 = buffer.data(gg + 138);
    const auto *gg_139 = buffer.data(gg + 139);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_141 = buffer.data(gg + 141);
    const auto *gg_142 = buffer.data(gg + 142);
    const auto *gg_143 = buffer.data(gg + 143);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_145 = buffer.data(gg + 145);
    const auto *gg_146 = buffer.data(gg + 146);
    const auto *gg_147 = buffer.data(gg + 147);
    const auto *gg_148 = buffer.data(gg + 148);
    const auto *gg_149 = buffer.data(gg + 149);

    const auto *ig_0 = buffer.data(ig + 0);
    const auto *ig_1 = buffer.data(ig + 1);
    const auto *ig_2 = buffer.data(ig + 2);
    const auto *ig_3 = buffer.data(ig + 3);
    const auto *ig_4 = buffer.data(ig + 4);
    const auto *ig_5 = buffer.data(ig + 5);
    const auto *ig_6 = buffer.data(ig + 6);
    const auto *ig_7 = buffer.data(ig + 7);
    const auto *ig_8 = buffer.data(ig + 8);
    const auto *ig_9 = buffer.data(ig + 9);
    const auto *ig_10 = buffer.data(ig + 10);
    const auto *ig_11 = buffer.data(ig + 11);
    const auto *ig_12 = buffer.data(ig + 12);
    const auto *ig_13 = buffer.data(ig + 13);
    const auto *ig_14 = buffer.data(ig + 14);
    const auto *ig_15 = buffer.data(ig + 15);
    const auto *ig_16 = buffer.data(ig + 16);
    const auto *ig_17 = buffer.data(ig + 17);
    const auto *ig_18 = buffer.data(ig + 18);
    const auto *ig_19 = buffer.data(ig + 19);
    const auto *ig_20 = buffer.data(ig + 20);
    const auto *ig_21 = buffer.data(ig + 21);
    const auto *ig_22 = buffer.data(ig + 22);
    const auto *ig_23 = buffer.data(ig + 23);
    const auto *ig_24 = buffer.data(ig + 24);
    const auto *ig_25 = buffer.data(ig + 25);
    const auto *ig_26 = buffer.data(ig + 26);
    const auto *ig_27 = buffer.data(ig + 27);
    const auto *ig_28 = buffer.data(ig + 28);
    const auto *ig_29 = buffer.data(ig + 29);
    const auto *ig_30 = buffer.data(ig + 30);
    const auto *ig_31 = buffer.data(ig + 31);
    const auto *ig_32 = buffer.data(ig + 32);
    const auto *ig_33 = buffer.data(ig + 33);
    const auto *ig_34 = buffer.data(ig + 34);
    const auto *ig_35 = buffer.data(ig + 35);
    const auto *ig_36 = buffer.data(ig + 36);
    const auto *ig_37 = buffer.data(ig + 37);
    const auto *ig_38 = buffer.data(ig + 38);
    const auto *ig_39 = buffer.data(ig + 39);
    const auto *ig_40 = buffer.data(ig + 40);
    const auto *ig_41 = buffer.data(ig + 41);
    const auto *ig_42 = buffer.data(ig + 42);
    const auto *ig_43 = buffer.data(ig + 43);
    const auto *ig_44 = buffer.data(ig + 44);
    const auto *ig_45 = buffer.data(ig + 45);
    const auto *ig_46 = buffer.data(ig + 46);
    const auto *ig_47 = buffer.data(ig + 47);
    const auto *ig_48 = buffer.data(ig + 48);
    const auto *ig_49 = buffer.data(ig + 49);
    const auto *ig_50 = buffer.data(ig + 50);
    const auto *ig_51 = buffer.data(ig + 51);
    const auto *ig_52 = buffer.data(ig + 52);
    const auto *ig_53 = buffer.data(ig + 53);
    const auto *ig_54 = buffer.data(ig + 54);
    const auto *ig_55 = buffer.data(ig + 55);
    const auto *ig_56 = buffer.data(ig + 56);
    const auto *ig_57 = buffer.data(ig + 57);
    const auto *ig_58 = buffer.data(ig + 58);
    const auto *ig_59 = buffer.data(ig + 59);
    const auto *ig_60 = buffer.data(ig + 60);
    const auto *ig_61 = buffer.data(ig + 61);
    const auto *ig_62 = buffer.data(ig + 62);
    const auto *ig_63 = buffer.data(ig + 63);
    const auto *ig_64 = buffer.data(ig + 64);
    const auto *ig_65 = buffer.data(ig + 65);
    const auto *ig_66 = buffer.data(ig + 66);
    const auto *ig_67 = buffer.data(ig + 67);
    const auto *ig_68 = buffer.data(ig + 68);
    const auto *ig_69 = buffer.data(ig + 69);
    const auto *ig_70 = buffer.data(ig + 70);
    const auto *ig_71 = buffer.data(ig + 71);
    const auto *ig_72 = buffer.data(ig + 72);
    const auto *ig_73 = buffer.data(ig + 73);
    const auto *ig_74 = buffer.data(ig + 74);
    const auto *ig_75 = buffer.data(ig + 75);
    const auto *ig_76 = buffer.data(ig + 76);
    const auto *ig_77 = buffer.data(ig + 77);
    const auto *ig_78 = buffer.data(ig + 78);
    const auto *ig_79 = buffer.data(ig + 79);
    const auto *ig_80 = buffer.data(ig + 80);
    const auto *ig_81 = buffer.data(ig + 81);
    const auto *ig_82 = buffer.data(ig + 82);
    const auto *ig_83 = buffer.data(ig + 83);
    const auto *ig_84 = buffer.data(ig + 84);
    const auto *ig_85 = buffer.data(ig + 85);
    const auto *ig_86 = buffer.data(ig + 86);
    const auto *ig_87 = buffer.data(ig + 87);
    const auto *ig_88 = buffer.data(ig + 88);
    const auto *ig_89 = buffer.data(ig + 89);
    const auto *ig_90 = buffer.data(ig + 90);
    const auto *ig_91 = buffer.data(ig + 91);
    const auto *ig_92 = buffer.data(ig + 92);
    const auto *ig_93 = buffer.data(ig + 93);
    const auto *ig_94 = buffer.data(ig + 94);
    const auto *ig_95 = buffer.data(ig + 95);
    const auto *ig_96 = buffer.data(ig + 96);
    const auto *ig_97 = buffer.data(ig + 97);
    const auto *ig_98 = buffer.data(ig + 98);
    const auto *ig_99 = buffer.data(ig + 99);
    const auto *ig_100 = buffer.data(ig + 100);
    const auto *ig_101 = buffer.data(ig + 101);
    const auto *ig_102 = buffer.data(ig + 102);
    const auto *ig_103 = buffer.data(ig + 103);
    const auto *ig_104 = buffer.data(ig + 104);
    const auto *ig_105 = buffer.data(ig + 105);
    const auto *ig_106 = buffer.data(ig + 106);
    const auto *ig_107 = buffer.data(ig + 107);
    const auto *ig_108 = buffer.data(ig + 108);
    const auto *ig_109 = buffer.data(ig + 109);
    const auto *ig_110 = buffer.data(ig + 110);
    const auto *ig_111 = buffer.data(ig + 111);
    const auto *ig_112 = buffer.data(ig + 112);
    const auto *ig_113 = buffer.data(ig + 113);
    const auto *ig_114 = buffer.data(ig + 114);
    const auto *ig_115 = buffer.data(ig + 115);
    const auto *ig_116 = buffer.data(ig + 116);
    const auto *ig_117 = buffer.data(ig + 117);
    const auto *ig_118 = buffer.data(ig + 118);
    const auto *ig_119 = buffer.data(ig + 119);
    const auto *ig_120 = buffer.data(ig + 120);
    const auto *ig_121 = buffer.data(ig + 121);
    const auto *ig_122 = buffer.data(ig + 122);
    const auto *ig_123 = buffer.data(ig + 123);
    const auto *ig_124 = buffer.data(ig + 124);
    const auto *ig_125 = buffer.data(ig + 125);
    const auto *ig_126 = buffer.data(ig + 126);
    const auto *ig_127 = buffer.data(ig + 127);
    const auto *ig_128 = buffer.data(ig + 128);
    const auto *ig_129 = buffer.data(ig + 129);
    const auto *ig_130 = buffer.data(ig + 130);
    const auto *ig_131 = buffer.data(ig + 131);
    const auto *ig_132 = buffer.data(ig + 132);
    const auto *ig_133 = buffer.data(ig + 133);
    const auto *ig_134 = buffer.data(ig + 134);
    const auto *ig_135 = buffer.data(ig + 135);
    const auto *ig_136 = buffer.data(ig + 136);
    const auto *ig_137 = buffer.data(ig + 137);
    const auto *ig_138 = buffer.data(ig + 138);
    const auto *ig_139 = buffer.data(ig + 139);
    const auto *ig_140 = buffer.data(ig + 140);
    const auto *ig_141 = buffer.data(ig + 141);
    const auto *ig_142 = buffer.data(ig + 142);
    const auto *ig_143 = buffer.data(ig + 143);
    const auto *ig_144 = buffer.data(ig + 144);
    const auto *ig_145 = buffer.data(ig + 145);
    const auto *ig_146 = buffer.data(ig + 146);
    const auto *ig_147 = buffer.data(ig + 147);
    const auto *ig_148 = buffer.data(ig + 148);
    const auto *ig_149 = buffer.data(ig + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, gg_0, gg_1, gg_2, gg_3, gg_4, ig_0, ig_1, \
                         ig_2, ig_3, ig_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -5.0 * gg_0[k]
                 + f_0 * ig_0[k];

        t_1[k] = -5.0 * gg_1[k]
                 + f_0 * ig_1[k];

        t_2[k] = -5.0 * gg_2[k]
                 + f_0 * ig_2[k];

        t_3[k] = -5.0 * gg_3[k]
                 + f_0 * ig_3[k];

        t_4[k] = -5.0 * gg_4[k]
                 + f_0 * ig_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, gg_5, gg_6, gg_7, gg_8, gg_9, ig_5, ig_6, \
                         ig_7, ig_8, ig_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -5.0 * gg_5[k]
                 + f_0 * ig_5[k];

        t_6[k] = -5.0 * gg_6[k]
                 + f_0 * ig_6[k];

        t_7[k] = -5.0 * gg_7[k]
                 + f_0 * ig_7[k];

        t_8[k] = -5.0 * gg_8[k]
                 + f_0 * ig_8[k];

        t_9[k] = -5.0 * gg_9[k]
                 + f_0 * ig_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, gg_10, gg_11, gg_12, gg_13, gg_14, \
                         ig_10, ig_11, ig_12, ig_13, ig_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -5.0 * gg_10[k]
                  + f_0 * ig_10[k];

        t_11[k] = -5.0 * gg_11[k]
                  + f_0 * ig_11[k];

        t_12[k] = -5.0 * gg_12[k]
                  + f_0 * ig_12[k];

        t_13[k] = -5.0 * gg_13[k]
                  + f_0 * ig_13[k];

        t_14[k] = -5.0 * gg_14[k]
                  + f_0 * ig_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, gg_15, gg_16, gg_17, gg_18, gg_19, \
                         ig_15, ig_16, ig_17, ig_18, ig_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -4.0 * gg_15[k]
                  + f_0 * ig_15[k];

        t_16[k] = -4.0 * gg_16[k]
                  + f_0 * ig_16[k];

        t_17[k] = -4.0 * gg_17[k]
                  + f_0 * ig_17[k];

        t_18[k] = -4.0 * gg_18[k]
                  + f_0 * ig_18[k];

        t_19[k] = -4.0 * gg_19[k]
                  + f_0 * ig_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, gg_20, gg_21, gg_22, gg_23, gg_24, \
                         ig_20, ig_21, ig_22, ig_23, ig_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -4.0 * gg_20[k]
                  + f_0 * ig_20[k];

        t_21[k] = -4.0 * gg_21[k]
                  + f_0 * ig_21[k];

        t_22[k] = -4.0 * gg_22[k]
                  + f_0 * ig_22[k];

        t_23[k] = -4.0 * gg_23[k]
                  + f_0 * ig_23[k];

        t_24[k] = -4.0 * gg_24[k]
                  + f_0 * ig_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, gg_25, gg_26, gg_27, gg_28, gg_29, \
                         ig_25, ig_26, ig_27, ig_28, ig_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -4.0 * gg_25[k]
                  + f_0 * ig_25[k];

        t_26[k] = -4.0 * gg_26[k]
                  + f_0 * ig_26[k];

        t_27[k] = -4.0 * gg_27[k]
                  + f_0 * ig_27[k];

        t_28[k] = -4.0 * gg_28[k]
                  + f_0 * ig_28[k];

        t_29[k] = -4.0 * gg_29[k]
                  + f_0 * ig_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, gg_30, gg_31, gg_32, gg_33, gg_34, \
                         ig_30, ig_31, ig_32, ig_33, ig_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -4.0 * gg_30[k]
                  + f_0 * ig_30[k];

        t_31[k] = -4.0 * gg_31[k]
                  + f_0 * ig_31[k];

        t_32[k] = -4.0 * gg_32[k]
                  + f_0 * ig_32[k];

        t_33[k] = -4.0 * gg_33[k]
                  + f_0 * ig_33[k];

        t_34[k] = -4.0 * gg_34[k]
                  + f_0 * ig_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, gg_35, gg_36, gg_37, gg_38, gg_39, \
                         ig_35, ig_36, ig_37, ig_38, ig_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -4.0 * gg_35[k]
                  + f_0 * ig_35[k];

        t_36[k] = -4.0 * gg_36[k]
                  + f_0 * ig_36[k];

        t_37[k] = -4.0 * gg_37[k]
                  + f_0 * ig_37[k];

        t_38[k] = -4.0 * gg_38[k]
                  + f_0 * ig_38[k];

        t_39[k] = -4.0 * gg_39[k]
                  + f_0 * ig_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, gg_40, gg_41, gg_42, gg_43, gg_44, \
                         ig_40, ig_41, ig_42, ig_43, ig_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -4.0 * gg_40[k]
                  + f_0 * ig_40[k];

        t_41[k] = -4.0 * gg_41[k]
                  + f_0 * ig_41[k];

        t_42[k] = -4.0 * gg_42[k]
                  + f_0 * ig_42[k];

        t_43[k] = -4.0 * gg_43[k]
                  + f_0 * ig_43[k];

        t_44[k] = -4.0 * gg_44[k]
                  + f_0 * ig_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, gg_45, gg_46, gg_47, gg_48, gg_49, \
                         ig_45, ig_46, ig_47, ig_48, ig_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -3.0 * gg_45[k]
                  + f_0 * ig_45[k];

        t_46[k] = -3.0 * gg_46[k]
                  + f_0 * ig_46[k];

        t_47[k] = -3.0 * gg_47[k]
                  + f_0 * ig_47[k];

        t_48[k] = -3.0 * gg_48[k]
                  + f_0 * ig_48[k];

        t_49[k] = -3.0 * gg_49[k]
                  + f_0 * ig_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, gg_50, gg_51, gg_52, gg_53, gg_54, \
                         ig_50, ig_51, ig_52, ig_53, ig_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -3.0 * gg_50[k]
                  + f_0 * ig_50[k];

        t_51[k] = -3.0 * gg_51[k]
                  + f_0 * ig_51[k];

        t_52[k] = -3.0 * gg_52[k]
                  + f_0 * ig_52[k];

        t_53[k] = -3.0 * gg_53[k]
                  + f_0 * ig_53[k];

        t_54[k] = -3.0 * gg_54[k]
                  + f_0 * ig_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, gg_55, gg_56, gg_57, gg_58, gg_59, \
                         ig_55, ig_56, ig_57, ig_58, ig_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -3.0 * gg_55[k]
                  + f_0 * ig_55[k];

        t_56[k] = -3.0 * gg_56[k]
                  + f_0 * ig_56[k];

        t_57[k] = -3.0 * gg_57[k]
                  + f_0 * ig_57[k];

        t_58[k] = -3.0 * gg_58[k]
                  + f_0 * ig_58[k];

        t_59[k] = -3.0 * gg_59[k]
                  + f_0 * ig_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, gg_60, gg_61, gg_62, gg_63, gg_64, \
                         ig_60, ig_61, ig_62, ig_63, ig_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * gg_60[k]
                  + f_0 * ig_60[k];

        t_61[k] = -3.0 * gg_61[k]
                  + f_0 * ig_61[k];

        t_62[k] = -3.0 * gg_62[k]
                  + f_0 * ig_62[k];

        t_63[k] = -3.0 * gg_63[k]
                  + f_0 * ig_63[k];

        t_64[k] = -3.0 * gg_64[k]
                  + f_0 * ig_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, gg_65, gg_66, gg_67, gg_68, gg_69, \
                         ig_65, ig_66, ig_67, ig_68, ig_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * gg_65[k]
                  + f_0 * ig_65[k];

        t_66[k] = -3.0 * gg_66[k]
                  + f_0 * ig_66[k];

        t_67[k] = -3.0 * gg_67[k]
                  + f_0 * ig_67[k];

        t_68[k] = -3.0 * gg_68[k]
                  + f_0 * ig_68[k];

        t_69[k] = -3.0 * gg_69[k]
                  + f_0 * ig_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, gg_70, gg_71, gg_72, gg_73, gg_74, \
                         ig_70, ig_71, ig_72, ig_73, ig_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -3.0 * gg_70[k]
                  + f_0 * ig_70[k];

        t_71[k] = -3.0 * gg_71[k]
                  + f_0 * ig_71[k];

        t_72[k] = -3.0 * gg_72[k]
                  + f_0 * ig_72[k];

        t_73[k] = -3.0 * gg_73[k]
                  + f_0 * ig_73[k];

        t_74[k] = -3.0 * gg_74[k]
                  + f_0 * ig_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, gg_75, gg_76, gg_77, gg_78, gg_79, \
                         ig_75, ig_76, ig_77, ig_78, ig_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -3.0 * gg_75[k]
                  + f_0 * ig_75[k];

        t_76[k] = -3.0 * gg_76[k]
                  + f_0 * ig_76[k];

        t_77[k] = -3.0 * gg_77[k]
                  + f_0 * ig_77[k];

        t_78[k] = -3.0 * gg_78[k]
                  + f_0 * ig_78[k];

        t_79[k] = -3.0 * gg_79[k]
                  + f_0 * ig_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, gg_80, gg_81, gg_82, gg_83, gg_84, \
                         ig_80, ig_81, ig_82, ig_83, ig_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * gg_80[k]
                  + f_0 * ig_80[k];

        t_81[k] = -3.0 * gg_81[k]
                  + f_0 * ig_81[k];

        t_82[k] = -3.0 * gg_82[k]
                  + f_0 * ig_82[k];

        t_83[k] = -3.0 * gg_83[k]
                  + f_0 * ig_83[k];

        t_84[k] = -3.0 * gg_84[k]
                  + f_0 * ig_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, gg_85, gg_86, gg_87, gg_88, gg_89, \
                         ig_85, ig_86, ig_87, ig_88, ig_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -3.0 * gg_85[k]
                  + f_0 * ig_85[k];

        t_86[k] = -3.0 * gg_86[k]
                  + f_0 * ig_86[k];

        t_87[k] = -3.0 * gg_87[k]
                  + f_0 * ig_87[k];

        t_88[k] = -3.0 * gg_88[k]
                  + f_0 * ig_88[k];

        t_89[k] = -3.0 * gg_89[k]
                  + f_0 * ig_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, gg_90, gg_91, gg_92, gg_93, gg_94, \
                         ig_90, ig_91, ig_92, ig_93, ig_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * gg_90[k]
                  + f_0 * ig_90[k];

        t_91[k] = -2.0 * gg_91[k]
                  + f_0 * ig_91[k];

        t_92[k] = -2.0 * gg_92[k]
                  + f_0 * ig_92[k];

        t_93[k] = -2.0 * gg_93[k]
                  + f_0 * ig_93[k];

        t_94[k] = -2.0 * gg_94[k]
                  + f_0 * ig_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, gg_95, gg_96, gg_97, gg_98, gg_99, \
                         ig_95, ig_96, ig_97, ig_98, ig_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * gg_95[k]
                  + f_0 * ig_95[k];

        t_96[k] = -2.0 * gg_96[k]
                  + f_0 * ig_96[k];

        t_97[k] = -2.0 * gg_97[k]
                  + f_0 * ig_97[k];

        t_98[k] = -2.0 * gg_98[k]
                  + f_0 * ig_98[k];

        t_99[k] = -2.0 * gg_99[k]
                  + f_0 * ig_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, gg_100, gg_101, gg_102, gg_103, \
                         gg_104, ig_100, ig_101, ig_102, ig_103, \
                         ig_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * gg_100[k]
                   + f_0 * ig_100[k];

        t_101[k] = -2.0 * gg_101[k]
                   + f_0 * ig_101[k];

        t_102[k] = -2.0 * gg_102[k]
                   + f_0 * ig_102[k];

        t_103[k] = -2.0 * gg_103[k]
                   + f_0 * ig_103[k];

        t_104[k] = -2.0 * gg_104[k]
                   + f_0 * ig_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, gg_105, gg_106, gg_107, gg_108, \
                         gg_109, ig_105, ig_106, ig_107, ig_108, \
                         ig_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * gg_105[k]
                   + f_0 * ig_105[k];

        t_106[k] = -2.0 * gg_106[k]
                   + f_0 * ig_106[k];

        t_107[k] = -2.0 * gg_107[k]
                   + f_0 * ig_107[k];

        t_108[k] = -2.0 * gg_108[k]
                   + f_0 * ig_108[k];

        t_109[k] = -2.0 * gg_109[k]
                   + f_0 * ig_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, gg_110, gg_111, gg_112, gg_113, \
                         gg_114, ig_110, ig_111, ig_112, ig_113, \
                         ig_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * gg_110[k]
                   + f_0 * ig_110[k];

        t_111[k] = -2.0 * gg_111[k]
                   + f_0 * ig_111[k];

        t_112[k] = -2.0 * gg_112[k]
                   + f_0 * ig_112[k];

        t_113[k] = -2.0 * gg_113[k]
                   + f_0 * ig_113[k];

        t_114[k] = -2.0 * gg_114[k]
                   + f_0 * ig_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, gg_115, gg_116, gg_117, gg_118, \
                         gg_119, ig_115, ig_116, ig_117, ig_118, \
                         ig_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * gg_115[k]
                   + f_0 * ig_115[k];

        t_116[k] = -2.0 * gg_116[k]
                   + f_0 * ig_116[k];

        t_117[k] = -2.0 * gg_117[k]
                   + f_0 * ig_117[k];

        t_118[k] = -2.0 * gg_118[k]
                   + f_0 * ig_118[k];

        t_119[k] = -2.0 * gg_119[k]
                   + f_0 * ig_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, gg_120, gg_121, gg_122, gg_123, \
                         gg_124, ig_120, ig_121, ig_122, ig_123, \
                         ig_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -2.0 * gg_120[k]
                   + f_0 * ig_120[k];

        t_121[k] = -2.0 * gg_121[k]
                   + f_0 * ig_121[k];

        t_122[k] = -2.0 * gg_122[k]
                   + f_0 * ig_122[k];

        t_123[k] = -2.0 * gg_123[k]
                   + f_0 * ig_123[k];

        t_124[k] = -2.0 * gg_124[k]
                   + f_0 * ig_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, gg_125, gg_126, gg_127, gg_128, \
                         gg_129, ig_125, ig_126, ig_127, ig_128, \
                         ig_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -2.0 * gg_125[k]
                   + f_0 * ig_125[k];

        t_126[k] = -2.0 * gg_126[k]
                   + f_0 * ig_126[k];

        t_127[k] = -2.0 * gg_127[k]
                   + f_0 * ig_127[k];

        t_128[k] = -2.0 * gg_128[k]
                   + f_0 * ig_128[k];

        t_129[k] = -2.0 * gg_129[k]
                   + f_0 * ig_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, gg_130, gg_131, gg_132, gg_133, \
                         gg_134, ig_130, ig_131, ig_132, ig_133, \
                         ig_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -2.0 * gg_130[k]
                   + f_0 * ig_130[k];

        t_131[k] = -2.0 * gg_131[k]
                   + f_0 * ig_131[k];

        t_132[k] = -2.0 * gg_132[k]
                   + f_0 * ig_132[k];

        t_133[k] = -2.0 * gg_133[k]
                   + f_0 * ig_133[k];

        t_134[k] = -2.0 * gg_134[k]
                   + f_0 * ig_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, gg_135, gg_136, gg_137, gg_138, \
                         gg_139, ig_135, ig_136, ig_137, ig_138, \
                         ig_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -2.0 * gg_135[k]
                   + f_0 * ig_135[k];

        t_136[k] = -2.0 * gg_136[k]
                   + f_0 * ig_136[k];

        t_137[k] = -2.0 * gg_137[k]
                   + f_0 * ig_137[k];

        t_138[k] = -2.0 * gg_138[k]
                   + f_0 * ig_138[k];

        t_139[k] = -2.0 * gg_139[k]
                   + f_0 * ig_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, gg_140, gg_141, gg_142, gg_143, \
                         gg_144, ig_140, ig_141, ig_142, ig_143, \
                         ig_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -2.0 * gg_140[k]
                   + f_0 * ig_140[k];

        t_141[k] = -2.0 * gg_141[k]
                   + f_0 * ig_141[k];

        t_142[k] = -2.0 * gg_142[k]
                   + f_0 * ig_142[k];

        t_143[k] = -2.0 * gg_143[k]
                   + f_0 * ig_143[k];

        t_144[k] = -2.0 * gg_144[k]
                   + f_0 * ig_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, gg_145, gg_146, gg_147, gg_148, \
                         gg_149, ig_145, ig_146, ig_147, ig_148, \
                         ig_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -2.0 * gg_145[k]
                   + f_0 * ig_145[k];

        t_146[k] = -2.0 * gg_146[k]
                   + f_0 * ig_146[k];

        t_147[k] = -2.0 * gg_147[k]
                   + f_0 * ig_147[k];

        t_148[k] = -2.0 * gg_148[k]
                   + f_0 * ig_148[k];

        t_149[k] = -2.0 * gg_149[k]
                   + f_0 * ig_149[k];
    }
}

static auto
compute_prim_geom_10_hg_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gg, const size_t ig,
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

    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_151 = buffer.data(gg + 151);
    const auto *gg_152 = buffer.data(gg + 152);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_154 = buffer.data(gg + 154);
    const auto *gg_155 = buffer.data(gg + 155);
    const auto *gg_156 = buffer.data(gg + 156);
    const auto *gg_157 = buffer.data(gg + 157);
    const auto *gg_158 = buffer.data(gg + 158);
    const auto *gg_159 = buffer.data(gg + 159);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_161 = buffer.data(gg + 161);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_163 = buffer.data(gg + 163);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_165 = buffer.data(gg + 165);
    const auto *gg_166 = buffer.data(gg + 166);
    const auto *gg_167 = buffer.data(gg + 167);
    const auto *gg_168 = buffer.data(gg + 168);
    const auto *gg_169 = buffer.data(gg + 169);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_171 = buffer.data(gg + 171);
    const auto *gg_172 = buffer.data(gg + 172);
    const auto *gg_173 = buffer.data(gg + 173);
    const auto *gg_174 = buffer.data(gg + 174);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_181 = buffer.data(gg + 181);
    const auto *gg_182 = buffer.data(gg + 182);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_184 = buffer.data(gg + 184);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_186 = buffer.data(gg + 186);
    const auto *gg_187 = buffer.data(gg + 187);
    const auto *gg_188 = buffer.data(gg + 188);
    const auto *gg_189 = buffer.data(gg + 189);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_195 = buffer.data(gg + 195);
    const auto *gg_196 = buffer.data(gg + 196);
    const auto *gg_197 = buffer.data(gg + 197);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_199 = buffer.data(gg + 199);
    const auto *gg_200 = buffer.data(gg + 200);
    const auto *gg_201 = buffer.data(gg + 201);
    const auto *gg_202 = buffer.data(gg + 202);
    const auto *gg_203 = buffer.data(gg + 203);
    const auto *gg_204 = buffer.data(gg + 204);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_211 = buffer.data(gg + 211);
    const auto *gg_212 = buffer.data(gg + 212);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_214 = buffer.data(gg + 214);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_216 = buffer.data(gg + 216);
    const auto *gg_217 = buffer.data(gg + 217);
    const auto *gg_218 = buffer.data(gg + 218);
    const auto *gg_219 = buffer.data(gg + 219);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_221 = buffer.data(gg + 221);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_223 = buffer.data(gg + 223);
    const auto *gg_224 = buffer.data(gg + 224);

    const auto *ig_150 = buffer.data(ig + 150);
    const auto *ig_151 = buffer.data(ig + 151);
    const auto *ig_152 = buffer.data(ig + 152);
    const auto *ig_153 = buffer.data(ig + 153);
    const auto *ig_154 = buffer.data(ig + 154);
    const auto *ig_155 = buffer.data(ig + 155);
    const auto *ig_156 = buffer.data(ig + 156);
    const auto *ig_157 = buffer.data(ig + 157);
    const auto *ig_158 = buffer.data(ig + 158);
    const auto *ig_159 = buffer.data(ig + 159);
    const auto *ig_160 = buffer.data(ig + 160);
    const auto *ig_161 = buffer.data(ig + 161);
    const auto *ig_162 = buffer.data(ig + 162);
    const auto *ig_163 = buffer.data(ig + 163);
    const auto *ig_164 = buffer.data(ig + 164);
    const auto *ig_165 = buffer.data(ig + 165);
    const auto *ig_166 = buffer.data(ig + 166);
    const auto *ig_167 = buffer.data(ig + 167);
    const auto *ig_168 = buffer.data(ig + 168);
    const auto *ig_169 = buffer.data(ig + 169);
    const auto *ig_170 = buffer.data(ig + 170);
    const auto *ig_171 = buffer.data(ig + 171);
    const auto *ig_172 = buffer.data(ig + 172);
    const auto *ig_173 = buffer.data(ig + 173);
    const auto *ig_174 = buffer.data(ig + 174);
    const auto *ig_175 = buffer.data(ig + 175);
    const auto *ig_176 = buffer.data(ig + 176);
    const auto *ig_177 = buffer.data(ig + 177);
    const auto *ig_178 = buffer.data(ig + 178);
    const auto *ig_179 = buffer.data(ig + 179);
    const auto *ig_180 = buffer.data(ig + 180);
    const auto *ig_181 = buffer.data(ig + 181);
    const auto *ig_182 = buffer.data(ig + 182);
    const auto *ig_183 = buffer.data(ig + 183);
    const auto *ig_184 = buffer.data(ig + 184);
    const auto *ig_185 = buffer.data(ig + 185);
    const auto *ig_186 = buffer.data(ig + 186);
    const auto *ig_187 = buffer.data(ig + 187);
    const auto *ig_188 = buffer.data(ig + 188);
    const auto *ig_189 = buffer.data(ig + 189);
    const auto *ig_190 = buffer.data(ig + 190);
    const auto *ig_191 = buffer.data(ig + 191);
    const auto *ig_192 = buffer.data(ig + 192);
    const auto *ig_193 = buffer.data(ig + 193);
    const auto *ig_194 = buffer.data(ig + 194);
    const auto *ig_195 = buffer.data(ig + 195);
    const auto *ig_196 = buffer.data(ig + 196);
    const auto *ig_197 = buffer.data(ig + 197);
    const auto *ig_198 = buffer.data(ig + 198);
    const auto *ig_199 = buffer.data(ig + 199);
    const auto *ig_200 = buffer.data(ig + 200);
    const auto *ig_201 = buffer.data(ig + 201);
    const auto *ig_202 = buffer.data(ig + 202);
    const auto *ig_203 = buffer.data(ig + 203);
    const auto *ig_204 = buffer.data(ig + 204);
    const auto *ig_205 = buffer.data(ig + 205);
    const auto *ig_206 = buffer.data(ig + 206);
    const auto *ig_207 = buffer.data(ig + 207);
    const auto *ig_208 = buffer.data(ig + 208);
    const auto *ig_209 = buffer.data(ig + 209);
    const auto *ig_210 = buffer.data(ig + 210);
    const auto *ig_211 = buffer.data(ig + 211);
    const auto *ig_212 = buffer.data(ig + 212);
    const auto *ig_213 = buffer.data(ig + 213);
    const auto *ig_214 = buffer.data(ig + 214);
    const auto *ig_215 = buffer.data(ig + 215);
    const auto *ig_216 = buffer.data(ig + 216);
    const auto *ig_217 = buffer.data(ig + 217);
    const auto *ig_218 = buffer.data(ig + 218);
    const auto *ig_219 = buffer.data(ig + 219);
    const auto *ig_220 = buffer.data(ig + 220);
    const auto *ig_221 = buffer.data(ig + 221);
    const auto *ig_222 = buffer.data(ig + 222);
    const auto *ig_223 = buffer.data(ig + 223);
    const auto *ig_224 = buffer.data(ig + 224);
    const auto *ig_225 = buffer.data(ig + 225);
    const auto *ig_226 = buffer.data(ig + 226);
    const auto *ig_227 = buffer.data(ig + 227);
    const auto *ig_228 = buffer.data(ig + 228);
    const auto *ig_229 = buffer.data(ig + 229);
    const auto *ig_230 = buffer.data(ig + 230);
    const auto *ig_231 = buffer.data(ig + 231);
    const auto *ig_232 = buffer.data(ig + 232);
    const auto *ig_233 = buffer.data(ig + 233);
    const auto *ig_234 = buffer.data(ig + 234);
    const auto *ig_235 = buffer.data(ig + 235);
    const auto *ig_236 = buffer.data(ig + 236);
    const auto *ig_237 = buffer.data(ig + 237);
    const auto *ig_238 = buffer.data(ig + 238);
    const auto *ig_239 = buffer.data(ig + 239);
    const auto *ig_240 = buffer.data(ig + 240);
    const auto *ig_241 = buffer.data(ig + 241);
    const auto *ig_242 = buffer.data(ig + 242);
    const auto *ig_243 = buffer.data(ig + 243);
    const auto *ig_244 = buffer.data(ig + 244);
    const auto *ig_245 = buffer.data(ig + 245);
    const auto *ig_246 = buffer.data(ig + 246);
    const auto *ig_247 = buffer.data(ig + 247);
    const auto *ig_248 = buffer.data(ig + 248);
    const auto *ig_249 = buffer.data(ig + 249);
    const auto *ig_250 = buffer.data(ig + 250);
    const auto *ig_251 = buffer.data(ig + 251);
    const auto *ig_252 = buffer.data(ig + 252);
    const auto *ig_253 = buffer.data(ig + 253);
    const auto *ig_254 = buffer.data(ig + 254);
    const auto *ig_255 = buffer.data(ig + 255);
    const auto *ig_256 = buffer.data(ig + 256);
    const auto *ig_257 = buffer.data(ig + 257);
    const auto *ig_258 = buffer.data(ig + 258);
    const auto *ig_259 = buffer.data(ig + 259);
    const auto *ig_260 = buffer.data(ig + 260);
    const auto *ig_261 = buffer.data(ig + 261);
    const auto *ig_262 = buffer.data(ig + 262);
    const auto *ig_263 = buffer.data(ig + 263);
    const auto *ig_264 = buffer.data(ig + 264);
    const auto *ig_265 = buffer.data(ig + 265);
    const auto *ig_266 = buffer.data(ig + 266);
    const auto *ig_267 = buffer.data(ig + 267);
    const auto *ig_268 = buffer.data(ig + 268);
    const auto *ig_269 = buffer.data(ig + 269);
    const auto *ig_270 = buffer.data(ig + 270);
    const auto *ig_271 = buffer.data(ig + 271);
    const auto *ig_272 = buffer.data(ig + 272);
    const auto *ig_273 = buffer.data(ig + 273);
    const auto *ig_274 = buffer.data(ig + 274);
    const auto *ig_275 = buffer.data(ig + 275);
    const auto *ig_276 = buffer.data(ig + 276);
    const auto *ig_277 = buffer.data(ig + 277);
    const auto *ig_278 = buffer.data(ig + 278);
    const auto *ig_279 = buffer.data(ig + 279);
    const auto *ig_280 = buffer.data(ig + 280);
    const auto *ig_281 = buffer.data(ig + 281);
    const auto *ig_282 = buffer.data(ig + 282);
    const auto *ig_283 = buffer.data(ig + 283);
    const auto *ig_284 = buffer.data(ig + 284);
    const auto *ig_285 = buffer.data(ig + 285);
    const auto *ig_286 = buffer.data(ig + 286);
    const auto *ig_287 = buffer.data(ig + 287);
    const auto *ig_288 = buffer.data(ig + 288);
    const auto *ig_289 = buffer.data(ig + 289);
    const auto *ig_290 = buffer.data(ig + 290);
    const auto *ig_291 = buffer.data(ig + 291);
    const auto *ig_292 = buffer.data(ig + 292);
    const auto *ig_293 = buffer.data(ig + 293);
    const auto *ig_294 = buffer.data(ig + 294);
    const auto *ig_295 = buffer.data(ig + 295);
    const auto *ig_296 = buffer.data(ig + 296);
    const auto *ig_297 = buffer.data(ig + 297);
    const auto *ig_298 = buffer.data(ig + 298);
    const auto *ig_299 = buffer.data(ig + 299);
    const auto *ig_300 = buffer.data(ig + 300);
    const auto *ig_301 = buffer.data(ig + 301);
    const auto *ig_302 = buffer.data(ig + 302);
    const auto *ig_303 = buffer.data(ig + 303);
    const auto *ig_304 = buffer.data(ig + 304);
    const auto *ig_305 = buffer.data(ig + 305);
    const auto *ig_306 = buffer.data(ig + 306);
    const auto *ig_307 = buffer.data(ig + 307);
    const auto *ig_308 = buffer.data(ig + 308);
    const auto *ig_309 = buffer.data(ig + 309);
    const auto *ig_310 = buffer.data(ig + 310);
    const auto *ig_311 = buffer.data(ig + 311);
    const auto *ig_312 = buffer.data(ig + 312);
    const auto *ig_313 = buffer.data(ig + 313);
    const auto *ig_314 = buffer.data(ig + 314);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, gg_150, gg_151, gg_152, gg_153, \
                         gg_154, ig_150, ig_151, ig_152, ig_153, \
                         ig_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -gg_150[k]
                   + f_0 * ig_150[k];

        t_151[k] = -gg_151[k]
                   + f_0 * ig_151[k];

        t_152[k] = -gg_152[k]
                   + f_0 * ig_152[k];

        t_153[k] = -gg_153[k]
                   + f_0 * ig_153[k];

        t_154[k] = -gg_154[k]
                   + f_0 * ig_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, gg_155, gg_156, gg_157, gg_158, \
                         gg_159, ig_155, ig_156, ig_157, ig_158, \
                         ig_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -gg_155[k]
                   + f_0 * ig_155[k];

        t_156[k] = -gg_156[k]
                   + f_0 * ig_156[k];

        t_157[k] = -gg_157[k]
                   + f_0 * ig_157[k];

        t_158[k] = -gg_158[k]
                   + f_0 * ig_158[k];

        t_159[k] = -gg_159[k]
                   + f_0 * ig_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, gg_160, gg_161, gg_162, gg_163, \
                         gg_164, ig_160, ig_161, ig_162, ig_163, \
                         ig_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -gg_160[k]
                   + f_0 * ig_160[k];

        t_161[k] = -gg_161[k]
                   + f_0 * ig_161[k];

        t_162[k] = -gg_162[k]
                   + f_0 * ig_162[k];

        t_163[k] = -gg_163[k]
                   + f_0 * ig_163[k];

        t_164[k] = -gg_164[k]
                   + f_0 * ig_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, gg_165, gg_166, gg_167, gg_168, \
                         gg_169, ig_165, ig_166, ig_167, ig_168, \
                         ig_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -gg_165[k]
                   + f_0 * ig_165[k];

        t_166[k] = -gg_166[k]
                   + f_0 * ig_166[k];

        t_167[k] = -gg_167[k]
                   + f_0 * ig_167[k];

        t_168[k] = -gg_168[k]
                   + f_0 * ig_168[k];

        t_169[k] = -gg_169[k]
                   + f_0 * ig_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, gg_170, gg_171, gg_172, gg_173, \
                         gg_174, ig_170, ig_171, ig_172, ig_173, \
                         ig_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -gg_170[k]
                   + f_0 * ig_170[k];

        t_171[k] = -gg_171[k]
                   + f_0 * ig_171[k];

        t_172[k] = -gg_172[k]
                   + f_0 * ig_172[k];

        t_173[k] = -gg_173[k]
                   + f_0 * ig_173[k];

        t_174[k] = -gg_174[k]
                   + f_0 * ig_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, gg_175, gg_176, gg_177, gg_178, \
                         gg_179, ig_175, ig_176, ig_177, ig_178, \
                         ig_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -gg_175[k]
                   + f_0 * ig_175[k];

        t_176[k] = -gg_176[k]
                   + f_0 * ig_176[k];

        t_177[k] = -gg_177[k]
                   + f_0 * ig_177[k];

        t_178[k] = -gg_178[k]
                   + f_0 * ig_178[k];

        t_179[k] = -gg_179[k]
                   + f_0 * ig_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, gg_180, gg_181, gg_182, gg_183, \
                         gg_184, ig_180, ig_181, ig_182, ig_183, \
                         ig_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -gg_180[k]
                   + f_0 * ig_180[k];

        t_181[k] = -gg_181[k]
                   + f_0 * ig_181[k];

        t_182[k] = -gg_182[k]
                   + f_0 * ig_182[k];

        t_183[k] = -gg_183[k]
                   + f_0 * ig_183[k];

        t_184[k] = -gg_184[k]
                   + f_0 * ig_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, gg_185, gg_186, gg_187, gg_188, \
                         gg_189, ig_185, ig_186, ig_187, ig_188, \
                         ig_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -gg_185[k]
                   + f_0 * ig_185[k];

        t_186[k] = -gg_186[k]
                   + f_0 * ig_186[k];

        t_187[k] = -gg_187[k]
                   + f_0 * ig_187[k];

        t_188[k] = -gg_188[k]
                   + f_0 * ig_188[k];

        t_189[k] = -gg_189[k]
                   + f_0 * ig_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, gg_190, gg_191, gg_192, gg_193, \
                         gg_194, ig_190, ig_191, ig_192, ig_193, \
                         ig_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -gg_190[k]
                   + f_0 * ig_190[k];

        t_191[k] = -gg_191[k]
                   + f_0 * ig_191[k];

        t_192[k] = -gg_192[k]
                   + f_0 * ig_192[k];

        t_193[k] = -gg_193[k]
                   + f_0 * ig_193[k];

        t_194[k] = -gg_194[k]
                   + f_0 * ig_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, gg_195, gg_196, gg_197, gg_198, \
                         gg_199, ig_195, ig_196, ig_197, ig_198, \
                         ig_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -gg_195[k]
                   + f_0 * ig_195[k];

        t_196[k] = -gg_196[k]
                   + f_0 * ig_196[k];

        t_197[k] = -gg_197[k]
                   + f_0 * ig_197[k];

        t_198[k] = -gg_198[k]
                   + f_0 * ig_198[k];

        t_199[k] = -gg_199[k]
                   + f_0 * ig_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, gg_200, gg_201, gg_202, gg_203, \
                         gg_204, ig_200, ig_201, ig_202, ig_203, \
                         ig_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -gg_200[k]
                   + f_0 * ig_200[k];

        t_201[k] = -gg_201[k]
                   + f_0 * ig_201[k];

        t_202[k] = -gg_202[k]
                   + f_0 * ig_202[k];

        t_203[k] = -gg_203[k]
                   + f_0 * ig_203[k];

        t_204[k] = -gg_204[k]
                   + f_0 * ig_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, gg_205, gg_206, gg_207, gg_208, \
                         gg_209, ig_205, ig_206, ig_207, ig_208, \
                         ig_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -gg_205[k]
                   + f_0 * ig_205[k];

        t_206[k] = -gg_206[k]
                   + f_0 * ig_206[k];

        t_207[k] = -gg_207[k]
                   + f_0 * ig_207[k];

        t_208[k] = -gg_208[k]
                   + f_0 * ig_208[k];

        t_209[k] = -gg_209[k]
                   + f_0 * ig_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, gg_210, gg_211, gg_212, gg_213, \
                         gg_214, ig_210, ig_211, ig_212, ig_213, \
                         ig_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -gg_210[k]
                   + f_0 * ig_210[k];

        t_211[k] = -gg_211[k]
                   + f_0 * ig_211[k];

        t_212[k] = -gg_212[k]
                   + f_0 * ig_212[k];

        t_213[k] = -gg_213[k]
                   + f_0 * ig_213[k];

        t_214[k] = -gg_214[k]
                   + f_0 * ig_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, gg_215, gg_216, gg_217, gg_218, \
                         gg_219, ig_215, ig_216, ig_217, ig_218, \
                         ig_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -gg_215[k]
                   + f_0 * ig_215[k];

        t_216[k] = -gg_216[k]
                   + f_0 * ig_216[k];

        t_217[k] = -gg_217[k]
                   + f_0 * ig_217[k];

        t_218[k] = -gg_218[k]
                   + f_0 * ig_218[k];

        t_219[k] = -gg_219[k]
                   + f_0 * ig_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, gg_220, gg_221, gg_222, gg_223, \
                         gg_224, ig_220, ig_221, ig_222, ig_223, \
                         ig_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -gg_220[k]
                   + f_0 * ig_220[k];

        t_221[k] = -gg_221[k]
                   + f_0 * ig_221[k];

        t_222[k] = -gg_222[k]
                   + f_0 * ig_222[k];

        t_223[k] = -gg_223[k]
                   + f_0 * ig_223[k];

        t_224[k] = -gg_224[k]
                   + f_0 * ig_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, t_230, t_231, t_232, ig_225, \
                         ig_226, ig_227, ig_228, ig_229, ig_230, ig_231, \
                         ig_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_0 * ig_225[k];

        t_226[k] = f_0 * ig_226[k];

        t_227[k] = f_0 * ig_227[k];

        t_228[k] = f_0 * ig_228[k];

        t_229[k] = f_0 * ig_229[k];

        t_230[k] = f_0 * ig_230[k];

        t_231[k] = f_0 * ig_231[k];

        t_232[k] = f_0 * ig_232[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, t_237, t_238, t_239, t_240, ig_233, \
                         ig_234, ig_235, ig_236, ig_237, ig_238, ig_239, \
                         ig_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_0 * ig_233[k];

        t_234[k] = f_0 * ig_234[k];

        t_235[k] = f_0 * ig_235[k];

        t_236[k] = f_0 * ig_236[k];

        t_237[k] = f_0 * ig_237[k];

        t_238[k] = f_0 * ig_238[k];

        t_239[k] = f_0 * ig_239[k];

        t_240[k] = f_0 * ig_240[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, t_246, t_247, t_248, ig_241, \
                         ig_242, ig_243, ig_244, ig_245, ig_246, ig_247, \
                         ig_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_0 * ig_241[k];

        t_242[k] = f_0 * ig_242[k];

        t_243[k] = f_0 * ig_243[k];

        t_244[k] = f_0 * ig_244[k];

        t_245[k] = f_0 * ig_245[k];

        t_246[k] = f_0 * ig_246[k];

        t_247[k] = f_0 * ig_247[k];

        t_248[k] = f_0 * ig_248[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, t_253, t_254, t_255, t_256, ig_249, \
                         ig_250, ig_251, ig_252, ig_253, ig_254, ig_255, \
                         ig_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_0 * ig_249[k];

        t_250[k] = f_0 * ig_250[k];

        t_251[k] = f_0 * ig_251[k];

        t_252[k] = f_0 * ig_252[k];

        t_253[k] = f_0 * ig_253[k];

        t_254[k] = f_0 * ig_254[k];

        t_255[k] = f_0 * ig_255[k];

        t_256[k] = f_0 * ig_256[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, t_262, t_263, t_264, ig_257, \
                         ig_258, ig_259, ig_260, ig_261, ig_262, ig_263, \
                         ig_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_0 * ig_257[k];

        t_258[k] = f_0 * ig_258[k];

        t_259[k] = f_0 * ig_259[k];

        t_260[k] = f_0 * ig_260[k];

        t_261[k] = f_0 * ig_261[k];

        t_262[k] = f_0 * ig_262[k];

        t_263[k] = f_0 * ig_263[k];

        t_264[k] = f_0 * ig_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, t_270, t_271, t_272, ig_265, \
                         ig_266, ig_267, ig_268, ig_269, ig_270, ig_271, \
                         ig_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_0 * ig_265[k];

        t_266[k] = f_0 * ig_266[k];

        t_267[k] = f_0 * ig_267[k];

        t_268[k] = f_0 * ig_268[k];

        t_269[k] = f_0 * ig_269[k];

        t_270[k] = f_0 * ig_270[k];

        t_271[k] = f_0 * ig_271[k];

        t_272[k] = f_0 * ig_272[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, t_278, t_279, t_280, ig_273, \
                         ig_274, ig_275, ig_276, ig_277, ig_278, ig_279, \
                         ig_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_0 * ig_273[k];

        t_274[k] = f_0 * ig_274[k];

        t_275[k] = f_0 * ig_275[k];

        t_276[k] = f_0 * ig_276[k];

        t_277[k] = f_0 * ig_277[k];

        t_278[k] = f_0 * ig_278[k];

        t_279[k] = f_0 * ig_279[k];

        t_280[k] = f_0 * ig_280[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, t_285, t_286, t_287, t_288, ig_281, \
                         ig_282, ig_283, ig_284, ig_285, ig_286, ig_287, \
                         ig_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_0 * ig_281[k];

        t_282[k] = f_0 * ig_282[k];

        t_283[k] = f_0 * ig_283[k];

        t_284[k] = f_0 * ig_284[k];

        t_285[k] = f_0 * ig_285[k];

        t_286[k] = f_0 * ig_286[k];

        t_287[k] = f_0 * ig_287[k];

        t_288[k] = f_0 * ig_288[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, t_294, t_295, t_296, ig_289, \
                         ig_290, ig_291, ig_292, ig_293, ig_294, ig_295, \
                         ig_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_0 * ig_289[k];

        t_290[k] = f_0 * ig_290[k];

        t_291[k] = f_0 * ig_291[k];

        t_292[k] = f_0 * ig_292[k];

        t_293[k] = f_0 * ig_293[k];

        t_294[k] = f_0 * ig_294[k];

        t_295[k] = f_0 * ig_295[k];

        t_296[k] = f_0 * ig_296[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, t_302, t_303, t_304, ig_297, \
                         ig_298, ig_299, ig_300, ig_301, ig_302, ig_303, \
                         ig_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_0 * ig_297[k];

        t_298[k] = f_0 * ig_298[k];

        t_299[k] = f_0 * ig_299[k];

        t_300[k] = f_0 * ig_300[k];

        t_301[k] = f_0 * ig_301[k];

        t_302[k] = f_0 * ig_302[k];

        t_303[k] = f_0 * ig_303[k];

        t_304[k] = f_0 * ig_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, t_310, t_311, t_312, ig_305, \
                         ig_306, ig_307, ig_308, ig_309, ig_310, ig_311, \
                         ig_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_0 * ig_305[k];

        t_306[k] = f_0 * ig_306[k];

        t_307[k] = f_0 * ig_307[k];

        t_308[k] = f_0 * ig_308[k];

        t_309[k] = f_0 * ig_309[k];

        t_310[k] = f_0 * ig_310[k];

        t_311[k] = f_0 * ig_311[k];

        t_312[k] = f_0 * ig_312[k];
    }

#pragma omp simd aligned(t_313, t_314, ig_313, ig_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_0 * ig_313[k];

        t_314[k] = f_0 * ig_314[k];
    }
}

auto
compute_prim_geom_10_hg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t gg, const size_t ig,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hg_electron_repulsion_0_piece0(buffer, target, gg, ig, ncols, alpha);

    compute_prim_geom_10_hg_electron_repulsion_0_piece1(buffer, target, gg, ig, ncols, alpha);
}

static auto
compute_prim_geom_10_hg_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gg, const size_t ig,
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

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_94 = buffer.data(gg + 94);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_97 = buffer.data(gg + 97);
    const auto *gg_98 = buffer.data(gg + 98);
    const auto *gg_99 = buffer.data(gg + 99);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_101 = buffer.data(gg + 101);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_106 = buffer.data(gg + 106);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_109 = buffer.data(gg + 109);

    const auto *ig_15 = buffer.data(ig + 15);
    const auto *ig_16 = buffer.data(ig + 16);
    const auto *ig_17 = buffer.data(ig + 17);
    const auto *ig_18 = buffer.data(ig + 18);
    const auto *ig_19 = buffer.data(ig + 19);
    const auto *ig_20 = buffer.data(ig + 20);
    const auto *ig_21 = buffer.data(ig + 21);
    const auto *ig_22 = buffer.data(ig + 22);
    const auto *ig_23 = buffer.data(ig + 23);
    const auto *ig_24 = buffer.data(ig + 24);
    const auto *ig_25 = buffer.data(ig + 25);
    const auto *ig_26 = buffer.data(ig + 26);
    const auto *ig_27 = buffer.data(ig + 27);
    const auto *ig_28 = buffer.data(ig + 28);
    const auto *ig_29 = buffer.data(ig + 29);
    const auto *ig_45 = buffer.data(ig + 45);
    const auto *ig_46 = buffer.data(ig + 46);
    const auto *ig_47 = buffer.data(ig + 47);
    const auto *ig_48 = buffer.data(ig + 48);
    const auto *ig_49 = buffer.data(ig + 49);
    const auto *ig_50 = buffer.data(ig + 50);
    const auto *ig_51 = buffer.data(ig + 51);
    const auto *ig_52 = buffer.data(ig + 52);
    const auto *ig_53 = buffer.data(ig + 53);
    const auto *ig_54 = buffer.data(ig + 54);
    const auto *ig_55 = buffer.data(ig + 55);
    const auto *ig_56 = buffer.data(ig + 56);
    const auto *ig_57 = buffer.data(ig + 57);
    const auto *ig_58 = buffer.data(ig + 58);
    const auto *ig_59 = buffer.data(ig + 59);
    const auto *ig_60 = buffer.data(ig + 60);
    const auto *ig_61 = buffer.data(ig + 61);
    const auto *ig_62 = buffer.data(ig + 62);
    const auto *ig_63 = buffer.data(ig + 63);
    const auto *ig_64 = buffer.data(ig + 64);
    const auto *ig_65 = buffer.data(ig + 65);
    const auto *ig_66 = buffer.data(ig + 66);
    const auto *ig_67 = buffer.data(ig + 67);
    const auto *ig_68 = buffer.data(ig + 68);
    const auto *ig_69 = buffer.data(ig + 69);
    const auto *ig_70 = buffer.data(ig + 70);
    const auto *ig_71 = buffer.data(ig + 71);
    const auto *ig_72 = buffer.data(ig + 72);
    const auto *ig_73 = buffer.data(ig + 73);
    const auto *ig_74 = buffer.data(ig + 74);
    const auto *ig_90 = buffer.data(ig + 90);
    const auto *ig_91 = buffer.data(ig + 91);
    const auto *ig_92 = buffer.data(ig + 92);
    const auto *ig_93 = buffer.data(ig + 93);
    const auto *ig_94 = buffer.data(ig + 94);
    const auto *ig_95 = buffer.data(ig + 95);
    const auto *ig_96 = buffer.data(ig + 96);
    const auto *ig_97 = buffer.data(ig + 97);
    const auto *ig_98 = buffer.data(ig + 98);
    const auto *ig_99 = buffer.data(ig + 99);
    const auto *ig_100 = buffer.data(ig + 100);
    const auto *ig_101 = buffer.data(ig + 101);
    const auto *ig_102 = buffer.data(ig + 102);
    const auto *ig_103 = buffer.data(ig + 103);
    const auto *ig_104 = buffer.data(ig + 104);
    const auto *ig_105 = buffer.data(ig + 105);
    const auto *ig_106 = buffer.data(ig + 106);
    const auto *ig_107 = buffer.data(ig + 107);
    const auto *ig_108 = buffer.data(ig + 108);
    const auto *ig_109 = buffer.data(ig + 109);
    const auto *ig_110 = buffer.data(ig + 110);
    const auto *ig_111 = buffer.data(ig + 111);
    const auto *ig_112 = buffer.data(ig + 112);
    const auto *ig_113 = buffer.data(ig + 113);
    const auto *ig_114 = buffer.data(ig + 114);
    const auto *ig_115 = buffer.data(ig + 115);
    const auto *ig_116 = buffer.data(ig + 116);
    const auto *ig_117 = buffer.data(ig + 117);
    const auto *ig_118 = buffer.data(ig + 118);
    const auto *ig_119 = buffer.data(ig + 119);
    const auto *ig_120 = buffer.data(ig + 120);
    const auto *ig_121 = buffer.data(ig + 121);
    const auto *ig_122 = buffer.data(ig + 122);
    const auto *ig_123 = buffer.data(ig + 123);
    const auto *ig_124 = buffer.data(ig + 124);
    const auto *ig_125 = buffer.data(ig + 125);
    const auto *ig_126 = buffer.data(ig + 126);
    const auto *ig_127 = buffer.data(ig + 127);
    const auto *ig_128 = buffer.data(ig + 128);
    const auto *ig_129 = buffer.data(ig + 129);
    const auto *ig_130 = buffer.data(ig + 130);
    const auto *ig_131 = buffer.data(ig + 131);
    const auto *ig_132 = buffer.data(ig + 132);
    const auto *ig_133 = buffer.data(ig + 133);
    const auto *ig_134 = buffer.data(ig + 134);
    const auto *ig_150 = buffer.data(ig + 150);
    const auto *ig_151 = buffer.data(ig + 151);
    const auto *ig_152 = buffer.data(ig + 152);
    const auto *ig_153 = buffer.data(ig + 153);
    const auto *ig_154 = buffer.data(ig + 154);
    const auto *ig_155 = buffer.data(ig + 155);
    const auto *ig_156 = buffer.data(ig + 156);
    const auto *ig_157 = buffer.data(ig + 157);
    const auto *ig_158 = buffer.data(ig + 158);
    const auto *ig_159 = buffer.data(ig + 159);
    const auto *ig_160 = buffer.data(ig + 160);
    const auto *ig_161 = buffer.data(ig + 161);
    const auto *ig_162 = buffer.data(ig + 162);
    const auto *ig_163 = buffer.data(ig + 163);
    const auto *ig_164 = buffer.data(ig + 164);
    const auto *ig_165 = buffer.data(ig + 165);
    const auto *ig_166 = buffer.data(ig + 166);
    const auto *ig_167 = buffer.data(ig + 167);
    const auto *ig_168 = buffer.data(ig + 168);
    const auto *ig_169 = buffer.data(ig + 169);
    const auto *ig_170 = buffer.data(ig + 170);
    const auto *ig_171 = buffer.data(ig + 171);
    const auto *ig_172 = buffer.data(ig + 172);
    const auto *ig_173 = buffer.data(ig + 173);
    const auto *ig_174 = buffer.data(ig + 174);
    const auto *ig_175 = buffer.data(ig + 175);
    const auto *ig_176 = buffer.data(ig + 176);
    const auto *ig_177 = buffer.data(ig + 177);
    const auto *ig_178 = buffer.data(ig + 178);
    const auto *ig_179 = buffer.data(ig + 179);
    const auto *ig_180 = buffer.data(ig + 180);
    const auto *ig_181 = buffer.data(ig + 181);
    const auto *ig_182 = buffer.data(ig + 182);
    const auto *ig_183 = buffer.data(ig + 183);
    const auto *ig_184 = buffer.data(ig + 184);
    const auto *ig_185 = buffer.data(ig + 185);
    const auto *ig_186 = buffer.data(ig + 186);
    const auto *ig_187 = buffer.data(ig + 187);
    const auto *ig_188 = buffer.data(ig + 188);
    const auto *ig_189 = buffer.data(ig + 189);
    const auto *ig_190 = buffer.data(ig + 190);
    const auto *ig_191 = buffer.data(ig + 191);
    const auto *ig_192 = buffer.data(ig + 192);
    const auto *ig_193 = buffer.data(ig + 193);
    const auto *ig_194 = buffer.data(ig + 194);
    const auto *ig_195 = buffer.data(ig + 195);
    const auto *ig_196 = buffer.data(ig + 196);
    const auto *ig_197 = buffer.data(ig + 197);
    const auto *ig_198 = buffer.data(ig + 198);
    const auto *ig_199 = buffer.data(ig + 199);
    const auto *ig_200 = buffer.data(ig + 200);
    const auto *ig_201 = buffer.data(ig + 201);
    const auto *ig_202 = buffer.data(ig + 202);
    const auto *ig_203 = buffer.data(ig + 203);
    const auto *ig_204 = buffer.data(ig + 204);
    const auto *ig_205 = buffer.data(ig + 205);
    const auto *ig_206 = buffer.data(ig + 206);
    const auto *ig_207 = buffer.data(ig + 207);
    const auto *ig_208 = buffer.data(ig + 208);
    const auto *ig_209 = buffer.data(ig + 209);
    const auto *ig_225 = buffer.data(ig + 225);
    const auto *ig_226 = buffer.data(ig + 226);
    const auto *ig_227 = buffer.data(ig + 227);
    const auto *ig_228 = buffer.data(ig + 228);
    const auto *ig_229 = buffer.data(ig + 229);
    const auto *ig_230 = buffer.data(ig + 230);
    const auto *ig_231 = buffer.data(ig + 231);
    const auto *ig_232 = buffer.data(ig + 232);
    const auto *ig_233 = buffer.data(ig + 233);
    const auto *ig_234 = buffer.data(ig + 234);
    const auto *ig_235 = buffer.data(ig + 235);
    const auto *ig_236 = buffer.data(ig + 236);
    const auto *ig_237 = buffer.data(ig + 237);
    const auto *ig_238 = buffer.data(ig + 238);
    const auto *ig_239 = buffer.data(ig + 239);
    const auto *ig_240 = buffer.data(ig + 240);
    const auto *ig_241 = buffer.data(ig + 241);
    const auto *ig_242 = buffer.data(ig + 242);
    const auto *ig_243 = buffer.data(ig + 243);
    const auto *ig_244 = buffer.data(ig + 244);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ig_15, ig_16, ig_17, ig_18, \
                         ig_19, ig_20, ig_21, ig_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ig_15[k];

        t_1[k] = f_0 * ig_16[k];

        t_2[k] = f_0 * ig_17[k];

        t_3[k] = f_0 * ig_18[k];

        t_4[k] = f_0 * ig_19[k];

        t_5[k] = f_0 * ig_20[k];

        t_6[k] = f_0 * ig_21[k];

        t_7[k] = f_0 * ig_22[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, ig_23, ig_24, ig_25, ig_26, \
                         ig_27, ig_28, ig_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ig_23[k];

        t_9[k] = f_0 * ig_24[k];

        t_10[k] = f_0 * ig_25[k];

        t_11[k] = f_0 * ig_26[k];

        t_12[k] = f_0 * ig_27[k];

        t_13[k] = f_0 * ig_28[k];

        t_14[k] = f_0 * ig_29[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, gg_0, gg_1, gg_2, gg_3, gg_4, ig_45, \
                         ig_46, ig_47, ig_48, ig_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -gg_0[k]
                  + f_0 * ig_45[k];

        t_16[k] = -gg_1[k]
                  + f_0 * ig_46[k];

        t_17[k] = -gg_2[k]
                  + f_0 * ig_47[k];

        t_18[k] = -gg_3[k]
                  + f_0 * ig_48[k];

        t_19[k] = -gg_4[k]
                  + f_0 * ig_49[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, gg_5, gg_6, gg_7, gg_8, gg_9, ig_50, \
                         ig_51, ig_52, ig_53, ig_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -gg_5[k]
                  + f_0 * ig_50[k];

        t_21[k] = -gg_6[k]
                  + f_0 * ig_51[k];

        t_22[k] = -gg_7[k]
                  + f_0 * ig_52[k];

        t_23[k] = -gg_8[k]
                  + f_0 * ig_53[k];

        t_24[k] = -gg_9[k]
                  + f_0 * ig_54[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, gg_10, gg_11, gg_12, gg_13, gg_14, \
                         ig_55, ig_56, ig_57, ig_58, ig_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -gg_10[k]
                  + f_0 * ig_55[k];

        t_26[k] = -gg_11[k]
                  + f_0 * ig_56[k];

        t_27[k] = -gg_12[k]
                  + f_0 * ig_57[k];

        t_28[k] = -gg_13[k]
                  + f_0 * ig_58[k];

        t_29[k] = -gg_14[k]
                  + f_0 * ig_59[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, t_36, t_37, ig_60, ig_61, ig_62, \
                         ig_63, ig_64, ig_65, ig_66, ig_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * ig_60[k];

        t_31[k] = f_0 * ig_61[k];

        t_32[k] = f_0 * ig_62[k];

        t_33[k] = f_0 * ig_63[k];

        t_34[k] = f_0 * ig_64[k];

        t_35[k] = f_0 * ig_65[k];

        t_36[k] = f_0 * ig_66[k];

        t_37[k] = f_0 * ig_67[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, ig_68, ig_69, ig_70, ig_71, \
                         ig_72, ig_73, ig_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * ig_68[k];

        t_39[k] = f_0 * ig_69[k];

        t_40[k] = f_0 * ig_70[k];

        t_41[k] = f_0 * ig_71[k];

        t_42[k] = f_0 * ig_72[k];

        t_43[k] = f_0 * ig_73[k];

        t_44[k] = f_0 * ig_74[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, gg_15, gg_16, gg_17, gg_18, gg_19, \
                         ig_90, ig_91, ig_92, ig_93, ig_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * gg_15[k]
                  + f_0 * ig_90[k];

        t_46[k] = -2.0 * gg_16[k]
                  + f_0 * ig_91[k];

        t_47[k] = -2.0 * gg_17[k]
                  + f_0 * ig_92[k];

        t_48[k] = -2.0 * gg_18[k]
                  + f_0 * ig_93[k];

        t_49[k] = -2.0 * gg_19[k]
                  + f_0 * ig_94[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, gg_20, gg_21, gg_22, gg_23, gg_24, \
                         ig_95, ig_96, ig_97, ig_98, ig_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * gg_20[k]
                  + f_0 * ig_95[k];

        t_51[k] = -2.0 * gg_21[k]
                  + f_0 * ig_96[k];

        t_52[k] = -2.0 * gg_22[k]
                  + f_0 * ig_97[k];

        t_53[k] = -2.0 * gg_23[k]
                  + f_0 * ig_98[k];

        t_54[k] = -2.0 * gg_24[k]
                  + f_0 * ig_99[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, gg_25, gg_26, gg_27, gg_28, gg_29, \
                         ig_100, ig_101, ig_102, ig_103, ig_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * gg_25[k]
                  + f_0 * ig_100[k];

        t_56[k] = -2.0 * gg_26[k]
                  + f_0 * ig_101[k];

        t_57[k] = -2.0 * gg_27[k]
                  + f_0 * ig_102[k];

        t_58[k] = -2.0 * gg_28[k]
                  + f_0 * ig_103[k];

        t_59[k] = -2.0 * gg_29[k]
                  + f_0 * ig_104[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, gg_30, gg_31, gg_32, gg_33, gg_34, \
                         ig_105, ig_106, ig_107, ig_108, ig_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -gg_30[k]
                  + f_0 * ig_105[k];

        t_61[k] = -gg_31[k]
                  + f_0 * ig_106[k];

        t_62[k] = -gg_32[k]
                  + f_0 * ig_107[k];

        t_63[k] = -gg_33[k]
                  + f_0 * ig_108[k];

        t_64[k] = -gg_34[k]
                  + f_0 * ig_109[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, gg_35, gg_36, gg_37, gg_38, gg_39, \
                         ig_110, ig_111, ig_112, ig_113, ig_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -gg_35[k]
                  + f_0 * ig_110[k];

        t_66[k] = -gg_36[k]
                  + f_0 * ig_111[k];

        t_67[k] = -gg_37[k]
                  + f_0 * ig_112[k];

        t_68[k] = -gg_38[k]
                  + f_0 * ig_113[k];

        t_69[k] = -gg_39[k]
                  + f_0 * ig_114[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, gg_40, gg_41, gg_42, gg_43, gg_44, \
                         ig_115, ig_116, ig_117, ig_118, ig_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -gg_40[k]
                  + f_0 * ig_115[k];

        t_71[k] = -gg_41[k]
                  + f_0 * ig_116[k];

        t_72[k] = -gg_42[k]
                  + f_0 * ig_117[k];

        t_73[k] = -gg_43[k]
                  + f_0 * ig_118[k];

        t_74[k] = -gg_44[k]
                  + f_0 * ig_119[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, t_80, t_81, t_82, ig_120, ig_121, \
                         ig_122, ig_123, ig_124, ig_125, ig_126, \
                         ig_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * ig_120[k];

        t_76[k] = f_0 * ig_121[k];

        t_77[k] = f_0 * ig_122[k];

        t_78[k] = f_0 * ig_123[k];

        t_79[k] = f_0 * ig_124[k];

        t_80[k] = f_0 * ig_125[k];

        t_81[k] = f_0 * ig_126[k];

        t_82[k] = f_0 * ig_127[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, ig_128, ig_129, ig_130, \
                         ig_131, ig_132, ig_133, ig_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * ig_128[k];

        t_84[k] = f_0 * ig_129[k];

        t_85[k] = f_0 * ig_130[k];

        t_86[k] = f_0 * ig_131[k];

        t_87[k] = f_0 * ig_132[k];

        t_88[k] = f_0 * ig_133[k];

        t_89[k] = f_0 * ig_134[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, gg_45, gg_46, gg_47, gg_48, gg_49, \
                         ig_150, ig_151, ig_152, ig_153, ig_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * gg_45[k]
                  + f_0 * ig_150[k];

        t_91[k] = -3.0 * gg_46[k]
                  + f_0 * ig_151[k];

        t_92[k] = -3.0 * gg_47[k]
                  + f_0 * ig_152[k];

        t_93[k] = -3.0 * gg_48[k]
                  + f_0 * ig_153[k];

        t_94[k] = -3.0 * gg_49[k]
                  + f_0 * ig_154[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, gg_50, gg_51, gg_52, gg_53, gg_54, \
                         ig_155, ig_156, ig_157, ig_158, ig_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * gg_50[k]
                  + f_0 * ig_155[k];

        t_96[k] = -3.0 * gg_51[k]
                  + f_0 * ig_156[k];

        t_97[k] = -3.0 * gg_52[k]
                  + f_0 * ig_157[k];

        t_98[k] = -3.0 * gg_53[k]
                  + f_0 * ig_158[k];

        t_99[k] = -3.0 * gg_54[k]
                  + f_0 * ig_159[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, gg_55, gg_56, gg_57, gg_58, gg_59, \
                         ig_160, ig_161, ig_162, ig_163, ig_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * gg_55[k]
                   + f_0 * ig_160[k];

        t_101[k] = -3.0 * gg_56[k]
                   + f_0 * ig_161[k];

        t_102[k] = -3.0 * gg_57[k]
                   + f_0 * ig_162[k];

        t_103[k] = -3.0 * gg_58[k]
                   + f_0 * ig_163[k];

        t_104[k] = -3.0 * gg_59[k]
                   + f_0 * ig_164[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, gg_60, gg_61, gg_62, gg_63, gg_64, \
                         ig_165, ig_166, ig_167, ig_168, ig_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * gg_60[k]
                   + f_0 * ig_165[k];

        t_106[k] = -2.0 * gg_61[k]
                   + f_0 * ig_166[k];

        t_107[k] = -2.0 * gg_62[k]
                   + f_0 * ig_167[k];

        t_108[k] = -2.0 * gg_63[k]
                   + f_0 * ig_168[k];

        t_109[k] = -2.0 * gg_64[k]
                   + f_0 * ig_169[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, gg_65, gg_66, gg_67, gg_68, gg_69, \
                         ig_170, ig_171, ig_172, ig_173, ig_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * gg_65[k]
                   + f_0 * ig_170[k];

        t_111[k] = -2.0 * gg_66[k]
                   + f_0 * ig_171[k];

        t_112[k] = -2.0 * gg_67[k]
                   + f_0 * ig_172[k];

        t_113[k] = -2.0 * gg_68[k]
                   + f_0 * ig_173[k];

        t_114[k] = -2.0 * gg_69[k]
                   + f_0 * ig_174[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, gg_70, gg_71, gg_72, gg_73, gg_74, \
                         ig_175, ig_176, ig_177, ig_178, ig_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * gg_70[k]
                   + f_0 * ig_175[k];

        t_116[k] = -2.0 * gg_71[k]
                   + f_0 * ig_176[k];

        t_117[k] = -2.0 * gg_72[k]
                   + f_0 * ig_177[k];

        t_118[k] = -2.0 * gg_73[k]
                   + f_0 * ig_178[k];

        t_119[k] = -2.0 * gg_74[k]
                   + f_0 * ig_179[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, gg_75, gg_76, gg_77, gg_78, gg_79, \
                         ig_180, ig_181, ig_182, ig_183, ig_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -gg_75[k]
                   + f_0 * ig_180[k];

        t_121[k] = -gg_76[k]
                   + f_0 * ig_181[k];

        t_122[k] = -gg_77[k]
                   + f_0 * ig_182[k];

        t_123[k] = -gg_78[k]
                   + f_0 * ig_183[k];

        t_124[k] = -gg_79[k]
                   + f_0 * ig_184[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, gg_80, gg_81, gg_82, gg_83, gg_84, \
                         ig_185, ig_186, ig_187, ig_188, ig_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -gg_80[k]
                   + f_0 * ig_185[k];

        t_126[k] = -gg_81[k]
                   + f_0 * ig_186[k];

        t_127[k] = -gg_82[k]
                   + f_0 * ig_187[k];

        t_128[k] = -gg_83[k]
                   + f_0 * ig_188[k];

        t_129[k] = -gg_84[k]
                   + f_0 * ig_189[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, gg_85, gg_86, gg_87, gg_88, gg_89, \
                         ig_190, ig_191, ig_192, ig_193, ig_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -gg_85[k]
                   + f_0 * ig_190[k];

        t_131[k] = -gg_86[k]
                   + f_0 * ig_191[k];

        t_132[k] = -gg_87[k]
                   + f_0 * ig_192[k];

        t_133[k] = -gg_88[k]
                   + f_0 * ig_193[k];

        t_134[k] = -gg_89[k]
                   + f_0 * ig_194[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, t_140, t_141, t_142, ig_195, \
                         ig_196, ig_197, ig_198, ig_199, ig_200, ig_201, \
                         ig_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * ig_195[k];

        t_136[k] = f_0 * ig_196[k];

        t_137[k] = f_0 * ig_197[k];

        t_138[k] = f_0 * ig_198[k];

        t_139[k] = f_0 * ig_199[k];

        t_140[k] = f_0 * ig_200[k];

        t_141[k] = f_0 * ig_201[k];

        t_142[k] = f_0 * ig_202[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, t_148, t_149, ig_203, ig_204, \
                         ig_205, ig_206, ig_207, ig_208, ig_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_0 * ig_203[k];

        t_144[k] = f_0 * ig_204[k];

        t_145[k] = f_0 * ig_205[k];

        t_146[k] = f_0 * ig_206[k];

        t_147[k] = f_0 * ig_207[k];

        t_148[k] = f_0 * ig_208[k];

        t_149[k] = f_0 * ig_209[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, gg_90, gg_91, gg_92, gg_93, gg_94, \
                         ig_225, ig_226, ig_227, ig_228, ig_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -4.0 * gg_90[k]
                   + f_0 * ig_225[k];

        t_151[k] = -4.0 * gg_91[k]
                   + f_0 * ig_226[k];

        t_152[k] = -4.0 * gg_92[k]
                   + f_0 * ig_227[k];

        t_153[k] = -4.0 * gg_93[k]
                   + f_0 * ig_228[k];

        t_154[k] = -4.0 * gg_94[k]
                   + f_0 * ig_229[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, gg_95, gg_96, gg_97, gg_98, gg_99, \
                         ig_230, ig_231, ig_232, ig_233, ig_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -4.0 * gg_95[k]
                   + f_0 * ig_230[k];

        t_156[k] = -4.0 * gg_96[k]
                   + f_0 * ig_231[k];

        t_157[k] = -4.0 * gg_97[k]
                   + f_0 * ig_232[k];

        t_158[k] = -4.0 * gg_98[k]
                   + f_0 * ig_233[k];

        t_159[k] = -4.0 * gg_99[k]
                   + f_0 * ig_234[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, gg_100, gg_101, gg_102, gg_103, \
                         gg_104, ig_235, ig_236, ig_237, ig_238, \
                         ig_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -4.0 * gg_100[k]
                   + f_0 * ig_235[k];

        t_161[k] = -4.0 * gg_101[k]
                   + f_0 * ig_236[k];

        t_162[k] = -4.0 * gg_102[k]
                   + f_0 * ig_237[k];

        t_163[k] = -4.0 * gg_103[k]
                   + f_0 * ig_238[k];

        t_164[k] = -4.0 * gg_104[k]
                   + f_0 * ig_239[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, gg_105, gg_106, gg_107, gg_108, \
                         gg_109, ig_240, ig_241, ig_242, ig_243, \
                         ig_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -3.0 * gg_105[k]
                   + f_0 * ig_240[k];

        t_166[k] = -3.0 * gg_106[k]
                   + f_0 * ig_241[k];

        t_167[k] = -3.0 * gg_107[k]
                   + f_0 * ig_242[k];

        t_168[k] = -3.0 * gg_108[k]
                   + f_0 * ig_243[k];

        t_169[k] = -3.0 * gg_109[k]
                   + f_0 * ig_244[k];
    }
}

static auto
compute_prim_geom_10_hg_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gg, const size_t ig,
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

    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_113 = buffer.data(gg + 113);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_124 = buffer.data(gg + 124);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_128 = buffer.data(gg + 128);
    const auto *gg_129 = buffer.data(gg + 129);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_134 = buffer.data(gg + 134);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_136 = buffer.data(gg + 136);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_138 = buffer.data(gg + 138);
    const auto *gg_139 = buffer.data(gg + 139);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_141 = buffer.data(gg + 141);
    const auto *gg_142 = buffer.data(gg + 142);
    const auto *gg_143 = buffer.data(gg + 143);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_145 = buffer.data(gg + 145);
    const auto *gg_146 = buffer.data(gg + 146);
    const auto *gg_147 = buffer.data(gg + 147);
    const auto *gg_148 = buffer.data(gg + 148);
    const auto *gg_149 = buffer.data(gg + 149);
    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_151 = buffer.data(gg + 151);
    const auto *gg_152 = buffer.data(gg + 152);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_154 = buffer.data(gg + 154);
    const auto *gg_155 = buffer.data(gg + 155);
    const auto *gg_156 = buffer.data(gg + 156);
    const auto *gg_157 = buffer.data(gg + 157);
    const auto *gg_158 = buffer.data(gg + 158);
    const auto *gg_159 = buffer.data(gg + 159);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_161 = buffer.data(gg + 161);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_163 = buffer.data(gg + 163);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_165 = buffer.data(gg + 165);
    const auto *gg_166 = buffer.data(gg + 166);
    const auto *gg_167 = buffer.data(gg + 167);
    const auto *gg_168 = buffer.data(gg + 168);
    const auto *gg_169 = buffer.data(gg + 169);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_171 = buffer.data(gg + 171);
    const auto *gg_172 = buffer.data(gg + 172);
    const auto *gg_173 = buffer.data(gg + 173);
    const auto *gg_174 = buffer.data(gg + 174);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_181 = buffer.data(gg + 181);
    const auto *gg_182 = buffer.data(gg + 182);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_184 = buffer.data(gg + 184);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_186 = buffer.data(gg + 186);
    const auto *gg_187 = buffer.data(gg + 187);
    const auto *gg_188 = buffer.data(gg + 188);
    const auto *gg_189 = buffer.data(gg + 189);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_195 = buffer.data(gg + 195);
    const auto *gg_196 = buffer.data(gg + 196);
    const auto *gg_197 = buffer.data(gg + 197);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_199 = buffer.data(gg + 199);
    const auto *gg_200 = buffer.data(gg + 200);
    const auto *gg_201 = buffer.data(gg + 201);
    const auto *gg_202 = buffer.data(gg + 202);
    const auto *gg_203 = buffer.data(gg + 203);
    const auto *gg_204 = buffer.data(gg + 204);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_211 = buffer.data(gg + 211);
    const auto *gg_212 = buffer.data(gg + 212);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_214 = buffer.data(gg + 214);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_216 = buffer.data(gg + 216);
    const auto *gg_217 = buffer.data(gg + 217);
    const auto *gg_218 = buffer.data(gg + 218);
    const auto *gg_219 = buffer.data(gg + 219);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_221 = buffer.data(gg + 221);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_223 = buffer.data(gg + 223);
    const auto *gg_224 = buffer.data(gg + 224);

    const auto *ig_245 = buffer.data(ig + 245);
    const auto *ig_246 = buffer.data(ig + 246);
    const auto *ig_247 = buffer.data(ig + 247);
    const auto *ig_248 = buffer.data(ig + 248);
    const auto *ig_249 = buffer.data(ig + 249);
    const auto *ig_250 = buffer.data(ig + 250);
    const auto *ig_251 = buffer.data(ig + 251);
    const auto *ig_252 = buffer.data(ig + 252);
    const auto *ig_253 = buffer.data(ig + 253);
    const auto *ig_254 = buffer.data(ig + 254);
    const auto *ig_255 = buffer.data(ig + 255);
    const auto *ig_256 = buffer.data(ig + 256);
    const auto *ig_257 = buffer.data(ig + 257);
    const auto *ig_258 = buffer.data(ig + 258);
    const auto *ig_259 = buffer.data(ig + 259);
    const auto *ig_260 = buffer.data(ig + 260);
    const auto *ig_261 = buffer.data(ig + 261);
    const auto *ig_262 = buffer.data(ig + 262);
    const auto *ig_263 = buffer.data(ig + 263);
    const auto *ig_264 = buffer.data(ig + 264);
    const auto *ig_265 = buffer.data(ig + 265);
    const auto *ig_266 = buffer.data(ig + 266);
    const auto *ig_267 = buffer.data(ig + 267);
    const auto *ig_268 = buffer.data(ig + 268);
    const auto *ig_269 = buffer.data(ig + 269);
    const auto *ig_270 = buffer.data(ig + 270);
    const auto *ig_271 = buffer.data(ig + 271);
    const auto *ig_272 = buffer.data(ig + 272);
    const auto *ig_273 = buffer.data(ig + 273);
    const auto *ig_274 = buffer.data(ig + 274);
    const auto *ig_275 = buffer.data(ig + 275);
    const auto *ig_276 = buffer.data(ig + 276);
    const auto *ig_277 = buffer.data(ig + 277);
    const auto *ig_278 = buffer.data(ig + 278);
    const auto *ig_279 = buffer.data(ig + 279);
    const auto *ig_280 = buffer.data(ig + 280);
    const auto *ig_281 = buffer.data(ig + 281);
    const auto *ig_282 = buffer.data(ig + 282);
    const auto *ig_283 = buffer.data(ig + 283);
    const auto *ig_284 = buffer.data(ig + 284);
    const auto *ig_285 = buffer.data(ig + 285);
    const auto *ig_286 = buffer.data(ig + 286);
    const auto *ig_287 = buffer.data(ig + 287);
    const auto *ig_288 = buffer.data(ig + 288);
    const auto *ig_289 = buffer.data(ig + 289);
    const auto *ig_290 = buffer.data(ig + 290);
    const auto *ig_291 = buffer.data(ig + 291);
    const auto *ig_292 = buffer.data(ig + 292);
    const auto *ig_293 = buffer.data(ig + 293);
    const auto *ig_294 = buffer.data(ig + 294);
    const auto *ig_295 = buffer.data(ig + 295);
    const auto *ig_296 = buffer.data(ig + 296);
    const auto *ig_297 = buffer.data(ig + 297);
    const auto *ig_298 = buffer.data(ig + 298);
    const auto *ig_299 = buffer.data(ig + 299);
    const auto *ig_315 = buffer.data(ig + 315);
    const auto *ig_316 = buffer.data(ig + 316);
    const auto *ig_317 = buffer.data(ig + 317);
    const auto *ig_318 = buffer.data(ig + 318);
    const auto *ig_319 = buffer.data(ig + 319);
    const auto *ig_320 = buffer.data(ig + 320);
    const auto *ig_321 = buffer.data(ig + 321);
    const auto *ig_322 = buffer.data(ig + 322);
    const auto *ig_323 = buffer.data(ig + 323);
    const auto *ig_324 = buffer.data(ig + 324);
    const auto *ig_325 = buffer.data(ig + 325);
    const auto *ig_326 = buffer.data(ig + 326);
    const auto *ig_327 = buffer.data(ig + 327);
    const auto *ig_328 = buffer.data(ig + 328);
    const auto *ig_329 = buffer.data(ig + 329);
    const auto *ig_330 = buffer.data(ig + 330);
    const auto *ig_331 = buffer.data(ig + 331);
    const auto *ig_332 = buffer.data(ig + 332);
    const auto *ig_333 = buffer.data(ig + 333);
    const auto *ig_334 = buffer.data(ig + 334);
    const auto *ig_335 = buffer.data(ig + 335);
    const auto *ig_336 = buffer.data(ig + 336);
    const auto *ig_337 = buffer.data(ig + 337);
    const auto *ig_338 = buffer.data(ig + 338);
    const auto *ig_339 = buffer.data(ig + 339);
    const auto *ig_340 = buffer.data(ig + 340);
    const auto *ig_341 = buffer.data(ig + 341);
    const auto *ig_342 = buffer.data(ig + 342);
    const auto *ig_343 = buffer.data(ig + 343);
    const auto *ig_344 = buffer.data(ig + 344);
    const auto *ig_345 = buffer.data(ig + 345);
    const auto *ig_346 = buffer.data(ig + 346);
    const auto *ig_347 = buffer.data(ig + 347);
    const auto *ig_348 = buffer.data(ig + 348);
    const auto *ig_349 = buffer.data(ig + 349);
    const auto *ig_350 = buffer.data(ig + 350);
    const auto *ig_351 = buffer.data(ig + 351);
    const auto *ig_352 = buffer.data(ig + 352);
    const auto *ig_353 = buffer.data(ig + 353);
    const auto *ig_354 = buffer.data(ig + 354);
    const auto *ig_355 = buffer.data(ig + 355);
    const auto *ig_356 = buffer.data(ig + 356);
    const auto *ig_357 = buffer.data(ig + 357);
    const auto *ig_358 = buffer.data(ig + 358);
    const auto *ig_359 = buffer.data(ig + 359);
    const auto *ig_360 = buffer.data(ig + 360);
    const auto *ig_361 = buffer.data(ig + 361);
    const auto *ig_362 = buffer.data(ig + 362);
    const auto *ig_363 = buffer.data(ig + 363);
    const auto *ig_364 = buffer.data(ig + 364);
    const auto *ig_365 = buffer.data(ig + 365);
    const auto *ig_366 = buffer.data(ig + 366);
    const auto *ig_367 = buffer.data(ig + 367);
    const auto *ig_368 = buffer.data(ig + 368);
    const auto *ig_369 = buffer.data(ig + 369);
    const auto *ig_370 = buffer.data(ig + 370);
    const auto *ig_371 = buffer.data(ig + 371);
    const auto *ig_372 = buffer.data(ig + 372);
    const auto *ig_373 = buffer.data(ig + 373);
    const auto *ig_374 = buffer.data(ig + 374);
    const auto *ig_375 = buffer.data(ig + 375);
    const auto *ig_376 = buffer.data(ig + 376);
    const auto *ig_377 = buffer.data(ig + 377);
    const auto *ig_378 = buffer.data(ig + 378);
    const auto *ig_379 = buffer.data(ig + 379);
    const auto *ig_380 = buffer.data(ig + 380);
    const auto *ig_381 = buffer.data(ig + 381);
    const auto *ig_382 = buffer.data(ig + 382);
    const auto *ig_383 = buffer.data(ig + 383);
    const auto *ig_384 = buffer.data(ig + 384);
    const auto *ig_385 = buffer.data(ig + 385);
    const auto *ig_386 = buffer.data(ig + 386);
    const auto *ig_387 = buffer.data(ig + 387);
    const auto *ig_388 = buffer.data(ig + 388);
    const auto *ig_389 = buffer.data(ig + 389);
    const auto *ig_390 = buffer.data(ig + 390);
    const auto *ig_391 = buffer.data(ig + 391);
    const auto *ig_392 = buffer.data(ig + 392);
    const auto *ig_393 = buffer.data(ig + 393);
    const auto *ig_394 = buffer.data(ig + 394);
    const auto *ig_395 = buffer.data(ig + 395);
    const auto *ig_396 = buffer.data(ig + 396);
    const auto *ig_397 = buffer.data(ig + 397);
    const auto *ig_398 = buffer.data(ig + 398);
    const auto *ig_399 = buffer.data(ig + 399);
    const auto *ig_400 = buffer.data(ig + 400);
    const auto *ig_401 = buffer.data(ig + 401);
    const auto *ig_402 = buffer.data(ig + 402);
    const auto *ig_403 = buffer.data(ig + 403);
    const auto *ig_404 = buffer.data(ig + 404);

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, gg_110, gg_111, gg_112, gg_113, \
                         gg_114, ig_245, ig_246, ig_247, ig_248, \
                         ig_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * gg_110[k]
                   + f_0 * ig_245[k];

        t_171[k] = -3.0 * gg_111[k]
                   + f_0 * ig_246[k];

        t_172[k] = -3.0 * gg_112[k]
                   + f_0 * ig_247[k];

        t_173[k] = -3.0 * gg_113[k]
                   + f_0 * ig_248[k];

        t_174[k] = -3.0 * gg_114[k]
                   + f_0 * ig_249[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, gg_115, gg_116, gg_117, gg_118, \
                         gg_119, ig_250, ig_251, ig_252, ig_253, \
                         ig_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * gg_115[k]
                   + f_0 * ig_250[k];

        t_176[k] = -3.0 * gg_116[k]
                   + f_0 * ig_251[k];

        t_177[k] = -3.0 * gg_117[k]
                   + f_0 * ig_252[k];

        t_178[k] = -3.0 * gg_118[k]
                   + f_0 * ig_253[k];

        t_179[k] = -3.0 * gg_119[k]
                   + f_0 * ig_254[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, gg_120, gg_121, gg_122, gg_123, \
                         gg_124, ig_255, ig_256, ig_257, ig_258, \
                         ig_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * gg_120[k]
                   + f_0 * ig_255[k];

        t_181[k] = -2.0 * gg_121[k]
                   + f_0 * ig_256[k];

        t_182[k] = -2.0 * gg_122[k]
                   + f_0 * ig_257[k];

        t_183[k] = -2.0 * gg_123[k]
                   + f_0 * ig_258[k];

        t_184[k] = -2.0 * gg_124[k]
                   + f_0 * ig_259[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, gg_125, gg_126, gg_127, gg_128, \
                         gg_129, ig_260, ig_261, ig_262, ig_263, \
                         ig_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * gg_125[k]
                   + f_0 * ig_260[k];

        t_186[k] = -2.0 * gg_126[k]
                   + f_0 * ig_261[k];

        t_187[k] = -2.0 * gg_127[k]
                   + f_0 * ig_262[k];

        t_188[k] = -2.0 * gg_128[k]
                   + f_0 * ig_263[k];

        t_189[k] = -2.0 * gg_129[k]
                   + f_0 * ig_264[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, gg_130, gg_131, gg_132, gg_133, \
                         gg_134, ig_265, ig_266, ig_267, ig_268, \
                         ig_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * gg_130[k]
                   + f_0 * ig_265[k];

        t_191[k] = -2.0 * gg_131[k]
                   + f_0 * ig_266[k];

        t_192[k] = -2.0 * gg_132[k]
                   + f_0 * ig_267[k];

        t_193[k] = -2.0 * gg_133[k]
                   + f_0 * ig_268[k];

        t_194[k] = -2.0 * gg_134[k]
                   + f_0 * ig_269[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, gg_135, gg_136, gg_137, gg_138, \
                         gg_139, ig_270, ig_271, ig_272, ig_273, \
                         ig_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -gg_135[k]
                   + f_0 * ig_270[k];

        t_196[k] = -gg_136[k]
                   + f_0 * ig_271[k];

        t_197[k] = -gg_137[k]
                   + f_0 * ig_272[k];

        t_198[k] = -gg_138[k]
                   + f_0 * ig_273[k];

        t_199[k] = -gg_139[k]
                   + f_0 * ig_274[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, gg_140, gg_141, gg_142, gg_143, \
                         gg_144, ig_275, ig_276, ig_277, ig_278, \
                         ig_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -gg_140[k]
                   + f_0 * ig_275[k];

        t_201[k] = -gg_141[k]
                   + f_0 * ig_276[k];

        t_202[k] = -gg_142[k]
                   + f_0 * ig_277[k];

        t_203[k] = -gg_143[k]
                   + f_0 * ig_278[k];

        t_204[k] = -gg_144[k]
                   + f_0 * ig_279[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, gg_145, gg_146, gg_147, gg_148, \
                         gg_149, ig_280, ig_281, ig_282, ig_283, \
                         ig_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -gg_145[k]
                   + f_0 * ig_280[k];

        t_206[k] = -gg_146[k]
                   + f_0 * ig_281[k];

        t_207[k] = -gg_147[k]
                   + f_0 * ig_282[k];

        t_208[k] = -gg_148[k]
                   + f_0 * ig_283[k];

        t_209[k] = -gg_149[k]
                   + f_0 * ig_284[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, t_217, ig_285, \
                         ig_286, ig_287, ig_288, ig_289, ig_290, ig_291, \
                         ig_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * ig_285[k];

        t_211[k] = f_0 * ig_286[k];

        t_212[k] = f_0 * ig_287[k];

        t_213[k] = f_0 * ig_288[k];

        t_214[k] = f_0 * ig_289[k];

        t_215[k] = f_0 * ig_290[k];

        t_216[k] = f_0 * ig_291[k];

        t_217[k] = f_0 * ig_292[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, t_224, ig_293, ig_294, \
                         ig_295, ig_296, ig_297, ig_298, ig_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * ig_293[k];

        t_219[k] = f_0 * ig_294[k];

        t_220[k] = f_0 * ig_295[k];

        t_221[k] = f_0 * ig_296[k];

        t_222[k] = f_0 * ig_297[k];

        t_223[k] = f_0 * ig_298[k];

        t_224[k] = f_0 * ig_299[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, gg_150, gg_151, gg_152, gg_153, \
                         gg_154, ig_315, ig_316, ig_317, ig_318, \
                         ig_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -5.0 * gg_150[k]
                   + f_0 * ig_315[k];

        t_226[k] = -5.0 * gg_151[k]
                   + f_0 * ig_316[k];

        t_227[k] = -5.0 * gg_152[k]
                   + f_0 * ig_317[k];

        t_228[k] = -5.0 * gg_153[k]
                   + f_0 * ig_318[k];

        t_229[k] = -5.0 * gg_154[k]
                   + f_0 * ig_319[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, gg_155, gg_156, gg_157, gg_158, \
                         gg_159, ig_320, ig_321, ig_322, ig_323, \
                         ig_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -5.0 * gg_155[k]
                   + f_0 * ig_320[k];

        t_231[k] = -5.0 * gg_156[k]
                   + f_0 * ig_321[k];

        t_232[k] = -5.0 * gg_157[k]
                   + f_0 * ig_322[k];

        t_233[k] = -5.0 * gg_158[k]
                   + f_0 * ig_323[k];

        t_234[k] = -5.0 * gg_159[k]
                   + f_0 * ig_324[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, gg_160, gg_161, gg_162, gg_163, \
                         gg_164, ig_325, ig_326, ig_327, ig_328, \
                         ig_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -5.0 * gg_160[k]
                   + f_0 * ig_325[k];

        t_236[k] = -5.0 * gg_161[k]
                   + f_0 * ig_326[k];

        t_237[k] = -5.0 * gg_162[k]
                   + f_0 * ig_327[k];

        t_238[k] = -5.0 * gg_163[k]
                   + f_0 * ig_328[k];

        t_239[k] = -5.0 * gg_164[k]
                   + f_0 * ig_329[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, gg_165, gg_166, gg_167, gg_168, \
                         gg_169, ig_330, ig_331, ig_332, ig_333, \
                         ig_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -4.0 * gg_165[k]
                   + f_0 * ig_330[k];

        t_241[k] = -4.0 * gg_166[k]
                   + f_0 * ig_331[k];

        t_242[k] = -4.0 * gg_167[k]
                   + f_0 * ig_332[k];

        t_243[k] = -4.0 * gg_168[k]
                   + f_0 * ig_333[k];

        t_244[k] = -4.0 * gg_169[k]
                   + f_0 * ig_334[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, gg_170, gg_171, gg_172, gg_173, \
                         gg_174, ig_335, ig_336, ig_337, ig_338, \
                         ig_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -4.0 * gg_170[k]
                   + f_0 * ig_335[k];

        t_246[k] = -4.0 * gg_171[k]
                   + f_0 * ig_336[k];

        t_247[k] = -4.0 * gg_172[k]
                   + f_0 * ig_337[k];

        t_248[k] = -4.0 * gg_173[k]
                   + f_0 * ig_338[k];

        t_249[k] = -4.0 * gg_174[k]
                   + f_0 * ig_339[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, gg_175, gg_176, gg_177, gg_178, \
                         gg_179, ig_340, ig_341, ig_342, ig_343, \
                         ig_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -4.0 * gg_175[k]
                   + f_0 * ig_340[k];

        t_251[k] = -4.0 * gg_176[k]
                   + f_0 * ig_341[k];

        t_252[k] = -4.0 * gg_177[k]
                   + f_0 * ig_342[k];

        t_253[k] = -4.0 * gg_178[k]
                   + f_0 * ig_343[k];

        t_254[k] = -4.0 * gg_179[k]
                   + f_0 * ig_344[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, gg_180, gg_181, gg_182, gg_183, \
                         gg_184, ig_345, ig_346, ig_347, ig_348, \
                         ig_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -3.0 * gg_180[k]
                   + f_0 * ig_345[k];

        t_256[k] = -3.0 * gg_181[k]
                   + f_0 * ig_346[k];

        t_257[k] = -3.0 * gg_182[k]
                   + f_0 * ig_347[k];

        t_258[k] = -3.0 * gg_183[k]
                   + f_0 * ig_348[k];

        t_259[k] = -3.0 * gg_184[k]
                   + f_0 * ig_349[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, gg_185, gg_186, gg_187, gg_188, \
                         gg_189, ig_350, ig_351, ig_352, ig_353, \
                         ig_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -3.0 * gg_185[k]
                   + f_0 * ig_350[k];

        t_261[k] = -3.0 * gg_186[k]
                   + f_0 * ig_351[k];

        t_262[k] = -3.0 * gg_187[k]
                   + f_0 * ig_352[k];

        t_263[k] = -3.0 * gg_188[k]
                   + f_0 * ig_353[k];

        t_264[k] = -3.0 * gg_189[k]
                   + f_0 * ig_354[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, gg_190, gg_191, gg_192, gg_193, \
                         gg_194, ig_355, ig_356, ig_357, ig_358, \
                         ig_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -3.0 * gg_190[k]
                   + f_0 * ig_355[k];

        t_266[k] = -3.0 * gg_191[k]
                   + f_0 * ig_356[k];

        t_267[k] = -3.0 * gg_192[k]
                   + f_0 * ig_357[k];

        t_268[k] = -3.0 * gg_193[k]
                   + f_0 * ig_358[k];

        t_269[k] = -3.0 * gg_194[k]
                   + f_0 * ig_359[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, gg_195, gg_196, gg_197, gg_198, \
                         gg_199, ig_360, ig_361, ig_362, ig_363, \
                         ig_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * gg_195[k]
                   + f_0 * ig_360[k];

        t_271[k] = -2.0 * gg_196[k]
                   + f_0 * ig_361[k];

        t_272[k] = -2.0 * gg_197[k]
                   + f_0 * ig_362[k];

        t_273[k] = -2.0 * gg_198[k]
                   + f_0 * ig_363[k];

        t_274[k] = -2.0 * gg_199[k]
                   + f_0 * ig_364[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, gg_200, gg_201, gg_202, gg_203, \
                         gg_204, ig_365, ig_366, ig_367, ig_368, \
                         ig_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * gg_200[k]
                   + f_0 * ig_365[k];

        t_276[k] = -2.0 * gg_201[k]
                   + f_0 * ig_366[k];

        t_277[k] = -2.0 * gg_202[k]
                   + f_0 * ig_367[k];

        t_278[k] = -2.0 * gg_203[k]
                   + f_0 * ig_368[k];

        t_279[k] = -2.0 * gg_204[k]
                   + f_0 * ig_369[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, gg_205, gg_206, gg_207, gg_208, \
                         gg_209, ig_370, ig_371, ig_372, ig_373, \
                         ig_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * gg_205[k]
                   + f_0 * ig_370[k];

        t_281[k] = -2.0 * gg_206[k]
                   + f_0 * ig_371[k];

        t_282[k] = -2.0 * gg_207[k]
                   + f_0 * ig_372[k];

        t_283[k] = -2.0 * gg_208[k]
                   + f_0 * ig_373[k];

        t_284[k] = -2.0 * gg_209[k]
                   + f_0 * ig_374[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, gg_210, gg_211, gg_212, gg_213, \
                         gg_214, ig_375, ig_376, ig_377, ig_378, \
                         ig_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -gg_210[k]
                   + f_0 * ig_375[k];

        t_286[k] = -gg_211[k]
                   + f_0 * ig_376[k];

        t_287[k] = -gg_212[k]
                   + f_0 * ig_377[k];

        t_288[k] = -gg_213[k]
                   + f_0 * ig_378[k];

        t_289[k] = -gg_214[k]
                   + f_0 * ig_379[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, gg_215, gg_216, gg_217, gg_218, \
                         gg_219, ig_380, ig_381, ig_382, ig_383, \
                         ig_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -gg_215[k]
                   + f_0 * ig_380[k];

        t_291[k] = -gg_216[k]
                   + f_0 * ig_381[k];

        t_292[k] = -gg_217[k]
                   + f_0 * ig_382[k];

        t_293[k] = -gg_218[k]
                   + f_0 * ig_383[k];

        t_294[k] = -gg_219[k]
                   + f_0 * ig_384[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, gg_220, gg_221, gg_222, gg_223, \
                         gg_224, ig_385, ig_386, ig_387, ig_388, \
                         ig_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -gg_220[k]
                   + f_0 * ig_385[k];

        t_296[k] = -gg_221[k]
                   + f_0 * ig_386[k];

        t_297[k] = -gg_222[k]
                   + f_0 * ig_387[k];

        t_298[k] = -gg_223[k]
                   + f_0 * ig_388[k];

        t_299[k] = -gg_224[k]
                   + f_0 * ig_389[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, t_305, t_306, t_307, ig_390, \
                         ig_391, ig_392, ig_393, ig_394, ig_395, ig_396, \
                         ig_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_0 * ig_390[k];

        t_301[k] = f_0 * ig_391[k];

        t_302[k] = f_0 * ig_392[k];

        t_303[k] = f_0 * ig_393[k];

        t_304[k] = f_0 * ig_394[k];

        t_305[k] = f_0 * ig_395[k];

        t_306[k] = f_0 * ig_396[k];

        t_307[k] = f_0 * ig_397[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, t_313, t_314, ig_398, ig_399, \
                         ig_400, ig_401, ig_402, ig_403, ig_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_0 * ig_398[k];

        t_309[k] = f_0 * ig_399[k];

        t_310[k] = f_0 * ig_400[k];

        t_311[k] = f_0 * ig_401[k];

        t_312[k] = f_0 * ig_402[k];

        t_313[k] = f_0 * ig_403[k];

        t_314[k] = f_0 * ig_404[k];
    }
}

auto
compute_prim_geom_10_hg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t gg, const size_t ig,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hg_electron_repulsion_1_piece0(buffer, target, gg, ig, ncols, alpha);

    compute_prim_geom_10_hg_electron_repulsion_1_piece1(buffer, target, gg, ig, ncols, alpha);
}

static auto
compute_prim_geom_10_hg_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gg, const size_t ig,
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

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_94 = buffer.data(gg + 94);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);

    const auto *ig_30 = buffer.data(ig + 30);
    const auto *ig_31 = buffer.data(ig + 31);
    const auto *ig_32 = buffer.data(ig + 32);
    const auto *ig_33 = buffer.data(ig + 33);
    const auto *ig_34 = buffer.data(ig + 34);
    const auto *ig_35 = buffer.data(ig + 35);
    const auto *ig_36 = buffer.data(ig + 36);
    const auto *ig_37 = buffer.data(ig + 37);
    const auto *ig_38 = buffer.data(ig + 38);
    const auto *ig_39 = buffer.data(ig + 39);
    const auto *ig_40 = buffer.data(ig + 40);
    const auto *ig_41 = buffer.data(ig + 41);
    const auto *ig_42 = buffer.data(ig + 42);
    const auto *ig_43 = buffer.data(ig + 43);
    const auto *ig_44 = buffer.data(ig + 44);
    const auto *ig_60 = buffer.data(ig + 60);
    const auto *ig_61 = buffer.data(ig + 61);
    const auto *ig_62 = buffer.data(ig + 62);
    const auto *ig_63 = buffer.data(ig + 63);
    const auto *ig_64 = buffer.data(ig + 64);
    const auto *ig_65 = buffer.data(ig + 65);
    const auto *ig_66 = buffer.data(ig + 66);
    const auto *ig_67 = buffer.data(ig + 67);
    const auto *ig_68 = buffer.data(ig + 68);
    const auto *ig_69 = buffer.data(ig + 69);
    const auto *ig_70 = buffer.data(ig + 70);
    const auto *ig_71 = buffer.data(ig + 71);
    const auto *ig_72 = buffer.data(ig + 72);
    const auto *ig_73 = buffer.data(ig + 73);
    const auto *ig_74 = buffer.data(ig + 74);
    const auto *ig_75 = buffer.data(ig + 75);
    const auto *ig_76 = buffer.data(ig + 76);
    const auto *ig_77 = buffer.data(ig + 77);
    const auto *ig_78 = buffer.data(ig + 78);
    const auto *ig_79 = buffer.data(ig + 79);
    const auto *ig_80 = buffer.data(ig + 80);
    const auto *ig_81 = buffer.data(ig + 81);
    const auto *ig_82 = buffer.data(ig + 82);
    const auto *ig_83 = buffer.data(ig + 83);
    const auto *ig_84 = buffer.data(ig + 84);
    const auto *ig_85 = buffer.data(ig + 85);
    const auto *ig_86 = buffer.data(ig + 86);
    const auto *ig_87 = buffer.data(ig + 87);
    const auto *ig_88 = buffer.data(ig + 88);
    const auto *ig_89 = buffer.data(ig + 89);
    const auto *ig_105 = buffer.data(ig + 105);
    const auto *ig_106 = buffer.data(ig + 106);
    const auto *ig_107 = buffer.data(ig + 107);
    const auto *ig_108 = buffer.data(ig + 108);
    const auto *ig_109 = buffer.data(ig + 109);
    const auto *ig_110 = buffer.data(ig + 110);
    const auto *ig_111 = buffer.data(ig + 111);
    const auto *ig_112 = buffer.data(ig + 112);
    const auto *ig_113 = buffer.data(ig + 113);
    const auto *ig_114 = buffer.data(ig + 114);
    const auto *ig_115 = buffer.data(ig + 115);
    const auto *ig_116 = buffer.data(ig + 116);
    const auto *ig_117 = buffer.data(ig + 117);
    const auto *ig_118 = buffer.data(ig + 118);
    const auto *ig_119 = buffer.data(ig + 119);
    const auto *ig_120 = buffer.data(ig + 120);
    const auto *ig_121 = buffer.data(ig + 121);
    const auto *ig_122 = buffer.data(ig + 122);
    const auto *ig_123 = buffer.data(ig + 123);
    const auto *ig_124 = buffer.data(ig + 124);
    const auto *ig_125 = buffer.data(ig + 125);
    const auto *ig_126 = buffer.data(ig + 126);
    const auto *ig_127 = buffer.data(ig + 127);
    const auto *ig_128 = buffer.data(ig + 128);
    const auto *ig_129 = buffer.data(ig + 129);
    const auto *ig_130 = buffer.data(ig + 130);
    const auto *ig_131 = buffer.data(ig + 131);
    const auto *ig_132 = buffer.data(ig + 132);
    const auto *ig_133 = buffer.data(ig + 133);
    const auto *ig_134 = buffer.data(ig + 134);
    const auto *ig_135 = buffer.data(ig + 135);
    const auto *ig_136 = buffer.data(ig + 136);
    const auto *ig_137 = buffer.data(ig + 137);
    const auto *ig_138 = buffer.data(ig + 138);
    const auto *ig_139 = buffer.data(ig + 139);
    const auto *ig_140 = buffer.data(ig + 140);
    const auto *ig_141 = buffer.data(ig + 141);
    const auto *ig_142 = buffer.data(ig + 142);
    const auto *ig_143 = buffer.data(ig + 143);
    const auto *ig_144 = buffer.data(ig + 144);
    const auto *ig_145 = buffer.data(ig + 145);
    const auto *ig_146 = buffer.data(ig + 146);
    const auto *ig_147 = buffer.data(ig + 147);
    const auto *ig_148 = buffer.data(ig + 148);
    const auto *ig_149 = buffer.data(ig + 149);
    const auto *ig_165 = buffer.data(ig + 165);
    const auto *ig_166 = buffer.data(ig + 166);
    const auto *ig_167 = buffer.data(ig + 167);
    const auto *ig_168 = buffer.data(ig + 168);
    const auto *ig_169 = buffer.data(ig + 169);
    const auto *ig_170 = buffer.data(ig + 170);
    const auto *ig_171 = buffer.data(ig + 171);
    const auto *ig_172 = buffer.data(ig + 172);
    const auto *ig_173 = buffer.data(ig + 173);
    const auto *ig_174 = buffer.data(ig + 174);
    const auto *ig_175 = buffer.data(ig + 175);
    const auto *ig_176 = buffer.data(ig + 176);
    const auto *ig_177 = buffer.data(ig + 177);
    const auto *ig_178 = buffer.data(ig + 178);
    const auto *ig_179 = buffer.data(ig + 179);
    const auto *ig_180 = buffer.data(ig + 180);
    const auto *ig_181 = buffer.data(ig + 181);
    const auto *ig_182 = buffer.data(ig + 182);
    const auto *ig_183 = buffer.data(ig + 183);
    const auto *ig_184 = buffer.data(ig + 184);
    const auto *ig_185 = buffer.data(ig + 185);
    const auto *ig_186 = buffer.data(ig + 186);
    const auto *ig_187 = buffer.data(ig + 187);
    const auto *ig_188 = buffer.data(ig + 188);
    const auto *ig_189 = buffer.data(ig + 189);
    const auto *ig_190 = buffer.data(ig + 190);
    const auto *ig_191 = buffer.data(ig + 191);
    const auto *ig_192 = buffer.data(ig + 192);
    const auto *ig_193 = buffer.data(ig + 193);
    const auto *ig_194 = buffer.data(ig + 194);
    const auto *ig_195 = buffer.data(ig + 195);
    const auto *ig_196 = buffer.data(ig + 196);
    const auto *ig_197 = buffer.data(ig + 197);
    const auto *ig_198 = buffer.data(ig + 198);
    const auto *ig_199 = buffer.data(ig + 199);
    const auto *ig_200 = buffer.data(ig + 200);
    const auto *ig_201 = buffer.data(ig + 201);
    const auto *ig_202 = buffer.data(ig + 202);
    const auto *ig_203 = buffer.data(ig + 203);
    const auto *ig_204 = buffer.data(ig + 204);
    const auto *ig_205 = buffer.data(ig + 205);
    const auto *ig_206 = buffer.data(ig + 206);
    const auto *ig_207 = buffer.data(ig + 207);
    const auto *ig_208 = buffer.data(ig + 208);
    const auto *ig_209 = buffer.data(ig + 209);
    const auto *ig_210 = buffer.data(ig + 210);
    const auto *ig_211 = buffer.data(ig + 211);
    const auto *ig_212 = buffer.data(ig + 212);
    const auto *ig_213 = buffer.data(ig + 213);
    const auto *ig_214 = buffer.data(ig + 214);
    const auto *ig_215 = buffer.data(ig + 215);
    const auto *ig_216 = buffer.data(ig + 216);
    const auto *ig_217 = buffer.data(ig + 217);
    const auto *ig_218 = buffer.data(ig + 218);
    const auto *ig_219 = buffer.data(ig + 219);
    const auto *ig_220 = buffer.data(ig + 220);
    const auto *ig_221 = buffer.data(ig + 221);
    const auto *ig_222 = buffer.data(ig + 222);
    const auto *ig_223 = buffer.data(ig + 223);
    const auto *ig_224 = buffer.data(ig + 224);
    const auto *ig_240 = buffer.data(ig + 240);
    const auto *ig_241 = buffer.data(ig + 241);
    const auto *ig_242 = buffer.data(ig + 242);
    const auto *ig_243 = buffer.data(ig + 243);
    const auto *ig_244 = buffer.data(ig + 244);
    const auto *ig_245 = buffer.data(ig + 245);
    const auto *ig_246 = buffer.data(ig + 246);
    const auto *ig_247 = buffer.data(ig + 247);
    const auto *ig_248 = buffer.data(ig + 248);
    const auto *ig_249 = buffer.data(ig + 249);
    const auto *ig_250 = buffer.data(ig + 250);
    const auto *ig_251 = buffer.data(ig + 251);
    const auto *ig_252 = buffer.data(ig + 252);
    const auto *ig_253 = buffer.data(ig + 253);
    const auto *ig_254 = buffer.data(ig + 254);
    const auto *ig_255 = buffer.data(ig + 255);
    const auto *ig_256 = buffer.data(ig + 256);
    const auto *ig_257 = buffer.data(ig + 257);
    const auto *ig_258 = buffer.data(ig + 258);
    const auto *ig_259 = buffer.data(ig + 259);
    const auto *ig_260 = buffer.data(ig + 260);
    const auto *ig_261 = buffer.data(ig + 261);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ig_30, ig_31, ig_32, ig_33, \
                         ig_34, ig_35, ig_36, ig_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ig_30[k];

        t_1[k] = f_0 * ig_31[k];

        t_2[k] = f_0 * ig_32[k];

        t_3[k] = f_0 * ig_33[k];

        t_4[k] = f_0 * ig_34[k];

        t_5[k] = f_0 * ig_35[k];

        t_6[k] = f_0 * ig_36[k];

        t_7[k] = f_0 * ig_37[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, ig_38, ig_39, ig_40, \
                         ig_41, ig_42, ig_43, ig_44, ig_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ig_38[k];

        t_9[k] = f_0 * ig_39[k];

        t_10[k] = f_0 * ig_40[k];

        t_11[k] = f_0 * ig_41[k];

        t_12[k] = f_0 * ig_42[k];

        t_13[k] = f_0 * ig_43[k];

        t_14[k] = f_0 * ig_44[k];

        t_15[k] = f_0 * ig_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, ig_61, ig_62, ig_63, \
                         ig_64, ig_65, ig_66, ig_67, ig_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * ig_61[k];

        t_17[k] = f_0 * ig_62[k];

        t_18[k] = f_0 * ig_63[k];

        t_19[k] = f_0 * ig_64[k];

        t_20[k] = f_0 * ig_65[k];

        t_21[k] = f_0 * ig_66[k];

        t_22[k] = f_0 * ig_67[k];

        t_23[k] = f_0 * ig_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, gg_0, ig_69, ig_70, ig_71, \
                         ig_72, ig_73, ig_74, ig_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * ig_69[k];

        t_25[k] = f_0 * ig_70[k];

        t_26[k] = f_0 * ig_71[k];

        t_27[k] = f_0 * ig_72[k];

        t_28[k] = f_0 * ig_73[k];

        t_29[k] = f_0 * ig_74[k];

        t_30[k] = -gg_0[k]
                  + f_0 * ig_75[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, gg_1, gg_2, gg_3, gg_4, gg_5, ig_76, \
                         ig_77, ig_78, ig_79, ig_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -gg_1[k]
                  + f_0 * ig_76[k];

        t_32[k] = -gg_2[k]
                  + f_0 * ig_77[k];

        t_33[k] = -gg_3[k]
                  + f_0 * ig_78[k];

        t_34[k] = -gg_4[k]
                  + f_0 * ig_79[k];

        t_35[k] = -gg_5[k]
                  + f_0 * ig_80[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, gg_6, gg_7, gg_8, gg_9, gg_10, ig_81, \
                         ig_82, ig_83, ig_84, ig_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -gg_6[k]
                  + f_0 * ig_81[k];

        t_37[k] = -gg_7[k]
                  + f_0 * ig_82[k];

        t_38[k] = -gg_8[k]
                  + f_0 * ig_83[k];

        t_39[k] = -gg_9[k]
                  + f_0 * ig_84[k];

        t_40[k] = -gg_10[k]
                  + f_0 * ig_85[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, t_46, gg_11, gg_12, gg_13, gg_14, \
                         ig_86, ig_87, ig_88, ig_89, ig_105, ig_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -gg_11[k]
                  + f_0 * ig_86[k];

        t_42[k] = -gg_12[k]
                  + f_0 * ig_87[k];

        t_43[k] = -gg_13[k]
                  + f_0 * ig_88[k];

        t_44[k] = -gg_14[k]
                  + f_0 * ig_89[k];

        t_45[k] = f_0 * ig_105[k];

        t_46[k] = f_0 * ig_106[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, t_53, t_54, ig_107, ig_108, \
                         ig_109, ig_110, ig_111, ig_112, ig_113, \
                         ig_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * ig_107[k];

        t_48[k] = f_0 * ig_108[k];

        t_49[k] = f_0 * ig_109[k];

        t_50[k] = f_0 * ig_110[k];

        t_51[k] = f_0 * ig_111[k];

        t_52[k] = f_0 * ig_112[k];

        t_53[k] = f_0 * ig_113[k];

        t_54[k] = f_0 * ig_114[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, gg_15, gg_16, ig_115, \
                         ig_116, ig_117, ig_118, ig_119, ig_120, \
                         ig_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * ig_115[k];

        t_56[k] = f_0 * ig_116[k];

        t_57[k] = f_0 * ig_117[k];

        t_58[k] = f_0 * ig_118[k];

        t_59[k] = f_0 * ig_119[k];

        t_60[k] = -gg_15[k]
                  + f_0 * ig_120[k];

        t_61[k] = -gg_16[k]
                  + f_0 * ig_121[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, gg_17, gg_18, gg_19, gg_20, gg_21, \
                         ig_122, ig_123, ig_124, ig_125, ig_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -gg_17[k]
                  + f_0 * ig_122[k];

        t_63[k] = -gg_18[k]
                  + f_0 * ig_123[k];

        t_64[k] = -gg_19[k]
                  + f_0 * ig_124[k];

        t_65[k] = -gg_20[k]
                  + f_0 * ig_125[k];

        t_66[k] = -gg_21[k]
                  + f_0 * ig_126[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, gg_22, gg_23, gg_24, gg_25, gg_26, \
                         ig_127, ig_128, ig_129, ig_130, ig_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -gg_22[k]
                  + f_0 * ig_127[k];

        t_68[k] = -gg_23[k]
                  + f_0 * ig_128[k];

        t_69[k] = -gg_24[k]
                  + f_0 * ig_129[k];

        t_70[k] = -gg_25[k]
                  + f_0 * ig_130[k];

        t_71[k] = -gg_26[k]
                  + f_0 * ig_131[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, gg_27, gg_28, gg_29, gg_30, gg_31, \
                         ig_132, ig_133, ig_134, ig_135, ig_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -gg_27[k]
                  + f_0 * ig_132[k];

        t_73[k] = -gg_28[k]
                  + f_0 * ig_133[k];

        t_74[k] = -gg_29[k]
                  + f_0 * ig_134[k];

        t_75[k] = -2.0 * gg_30[k]
                  + f_0 * ig_135[k];

        t_76[k] = -2.0 * gg_31[k]
                  + f_0 * ig_136[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, gg_32, gg_33, gg_34, gg_35, gg_36, \
                         ig_137, ig_138, ig_139, ig_140, ig_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -2.0 * gg_32[k]
                  + f_0 * ig_137[k];

        t_78[k] = -2.0 * gg_33[k]
                  + f_0 * ig_138[k];

        t_79[k] = -2.0 * gg_34[k]
                  + f_0 * ig_139[k];

        t_80[k] = -2.0 * gg_35[k]
                  + f_0 * ig_140[k];

        t_81[k] = -2.0 * gg_36[k]
                  + f_0 * ig_141[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, gg_37, gg_38, gg_39, gg_40, gg_41, \
                         ig_142, ig_143, ig_144, ig_145, ig_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -2.0 * gg_37[k]
                  + f_0 * ig_142[k];

        t_83[k] = -2.0 * gg_38[k]
                  + f_0 * ig_143[k];

        t_84[k] = -2.0 * gg_39[k]
                  + f_0 * ig_144[k];

        t_85[k] = -2.0 * gg_40[k]
                  + f_0 * ig_145[k];

        t_86[k] = -2.0 * gg_41[k]
                  + f_0 * ig_146[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, gg_42, gg_43, gg_44, ig_147, \
                         ig_148, ig_149, ig_165, ig_166, ig_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -2.0 * gg_42[k]
                  + f_0 * ig_147[k];

        t_88[k] = -2.0 * gg_43[k]
                  + f_0 * ig_148[k];

        t_89[k] = -2.0 * gg_44[k]
                  + f_0 * ig_149[k];

        t_90[k] = f_0 * ig_165[k];

        t_91[k] = f_0 * ig_166[k];

        t_92[k] = f_0 * ig_167[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, ig_168, ig_169, \
                         ig_170, ig_171, ig_172, ig_173, ig_174, \
                         ig_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * ig_168[k];

        t_94[k] = f_0 * ig_169[k];

        t_95[k] = f_0 * ig_170[k];

        t_96[k] = f_0 * ig_171[k];

        t_97[k] = f_0 * ig_172[k];

        t_98[k] = f_0 * ig_173[k];

        t_99[k] = f_0 * ig_174[k];

        t_100[k] = f_0 * ig_175[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, gg_45, gg_46, ig_176, \
                         ig_177, ig_178, ig_179, ig_180, ig_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * ig_176[k];

        t_102[k] = f_0 * ig_177[k];

        t_103[k] = f_0 * ig_178[k];

        t_104[k] = f_0 * ig_179[k];

        t_105[k] = -gg_45[k]
                   + f_0 * ig_180[k];

        t_106[k] = -gg_46[k]
                   + f_0 * ig_181[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, gg_47, gg_48, gg_49, gg_50, gg_51, \
                         ig_182, ig_183, ig_184, ig_185, ig_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -gg_47[k]
                   + f_0 * ig_182[k];

        t_108[k] = -gg_48[k]
                   + f_0 * ig_183[k];

        t_109[k] = -gg_49[k]
                   + f_0 * ig_184[k];

        t_110[k] = -gg_50[k]
                   + f_0 * ig_185[k];

        t_111[k] = -gg_51[k]
                   + f_0 * ig_186[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, gg_52, gg_53, gg_54, gg_55, gg_56, \
                         ig_187, ig_188, ig_189, ig_190, ig_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -gg_52[k]
                   + f_0 * ig_187[k];

        t_113[k] = -gg_53[k]
                   + f_0 * ig_188[k];

        t_114[k] = -gg_54[k]
                   + f_0 * ig_189[k];

        t_115[k] = -gg_55[k]
                   + f_0 * ig_190[k];

        t_116[k] = -gg_56[k]
                   + f_0 * ig_191[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, gg_57, gg_58, gg_59, gg_60, gg_61, \
                         ig_192, ig_193, ig_194, ig_195, ig_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -gg_57[k]
                   + f_0 * ig_192[k];

        t_118[k] = -gg_58[k]
                   + f_0 * ig_193[k];

        t_119[k] = -gg_59[k]
                   + f_0 * ig_194[k];

        t_120[k] = -2.0 * gg_60[k]
                   + f_0 * ig_195[k];

        t_121[k] = -2.0 * gg_61[k]
                   + f_0 * ig_196[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, gg_62, gg_63, gg_64, gg_65, gg_66, \
                         ig_197, ig_198, ig_199, ig_200, ig_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * gg_62[k]
                   + f_0 * ig_197[k];

        t_123[k] = -2.0 * gg_63[k]
                   + f_0 * ig_198[k];

        t_124[k] = -2.0 * gg_64[k]
                   + f_0 * ig_199[k];

        t_125[k] = -2.0 * gg_65[k]
                   + f_0 * ig_200[k];

        t_126[k] = -2.0 * gg_66[k]
                   + f_0 * ig_201[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, gg_67, gg_68, gg_69, gg_70, gg_71, \
                         ig_202, ig_203, ig_204, ig_205, ig_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * gg_67[k]
                   + f_0 * ig_202[k];

        t_128[k] = -2.0 * gg_68[k]
                   + f_0 * ig_203[k];

        t_129[k] = -2.0 * gg_69[k]
                   + f_0 * ig_204[k];

        t_130[k] = -2.0 * gg_70[k]
                   + f_0 * ig_205[k];

        t_131[k] = -2.0 * gg_71[k]
                   + f_0 * ig_206[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, gg_72, gg_73, gg_74, gg_75, gg_76, \
                         ig_207, ig_208, ig_209, ig_210, ig_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * gg_72[k]
                   + f_0 * ig_207[k];

        t_133[k] = -2.0 * gg_73[k]
                   + f_0 * ig_208[k];

        t_134[k] = -2.0 * gg_74[k]
                   + f_0 * ig_209[k];

        t_135[k] = -3.0 * gg_75[k]
                   + f_0 * ig_210[k];

        t_136[k] = -3.0 * gg_76[k]
                   + f_0 * ig_211[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, gg_77, gg_78, gg_79, gg_80, gg_81, \
                         ig_212, ig_213, ig_214, ig_215, ig_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -3.0 * gg_77[k]
                   + f_0 * ig_212[k];

        t_138[k] = -3.0 * gg_78[k]
                   + f_0 * ig_213[k];

        t_139[k] = -3.0 * gg_79[k]
                   + f_0 * ig_214[k];

        t_140[k] = -3.0 * gg_80[k]
                   + f_0 * ig_215[k];

        t_141[k] = -3.0 * gg_81[k]
                   + f_0 * ig_216[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, gg_82, gg_83, gg_84, gg_85, gg_86, \
                         ig_217, ig_218, ig_219, ig_220, ig_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -3.0 * gg_82[k]
                   + f_0 * ig_217[k];

        t_143[k] = -3.0 * gg_83[k]
                   + f_0 * ig_218[k];

        t_144[k] = -3.0 * gg_84[k]
                   + f_0 * ig_219[k];

        t_145[k] = -3.0 * gg_85[k]
                   + f_0 * ig_220[k];

        t_146[k] = -3.0 * gg_86[k]
                   + f_0 * ig_221[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, t_152, gg_87, gg_88, gg_89, \
                         ig_222, ig_223, ig_224, ig_240, ig_241, \
                         ig_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -3.0 * gg_87[k]
                   + f_0 * ig_222[k];

        t_148[k] = -3.0 * gg_88[k]
                   + f_0 * ig_223[k];

        t_149[k] = -3.0 * gg_89[k]
                   + f_0 * ig_224[k];

        t_150[k] = f_0 * ig_240[k];

        t_151[k] = f_0 * ig_241[k];

        t_152[k] = f_0 * ig_242[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, t_158, t_159, t_160, ig_243, \
                         ig_244, ig_245, ig_246, ig_247, ig_248, ig_249, \
                         ig_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_0 * ig_243[k];

        t_154[k] = f_0 * ig_244[k];

        t_155[k] = f_0 * ig_245[k];

        t_156[k] = f_0 * ig_246[k];

        t_157[k] = f_0 * ig_247[k];

        t_158[k] = f_0 * ig_248[k];

        t_159[k] = f_0 * ig_249[k];

        t_160[k] = f_0 * ig_250[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, gg_90, gg_91, ig_251, \
                         ig_252, ig_253, ig_254, ig_255, ig_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_0 * ig_251[k];

        t_162[k] = f_0 * ig_252[k];

        t_163[k] = f_0 * ig_253[k];

        t_164[k] = f_0 * ig_254[k];

        t_165[k] = -gg_90[k]
                   + f_0 * ig_255[k];

        t_166[k] = -gg_91[k]
                   + f_0 * ig_256[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, gg_92, gg_93, gg_94, gg_95, gg_96, \
                         ig_257, ig_258, ig_259, ig_260, ig_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -gg_92[k]
                   + f_0 * ig_257[k];

        t_168[k] = -gg_93[k]
                   + f_0 * ig_258[k];

        t_169[k] = -gg_94[k]
                   + f_0 * ig_259[k];

        t_170[k] = -gg_95[k]
                   + f_0 * ig_260[k];

        t_171[k] = -gg_96[k]
                   + f_0 * ig_261[k];
    }
}

static auto
compute_prim_geom_10_hg_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gg, const size_t ig,
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

    const auto *gg_97 = buffer.data(gg + 97);
    const auto *gg_98 = buffer.data(gg + 98);
    const auto *gg_99 = buffer.data(gg + 99);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_101 = buffer.data(gg + 101);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_106 = buffer.data(gg + 106);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_109 = buffer.data(gg + 109);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_113 = buffer.data(gg + 113);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_124 = buffer.data(gg + 124);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_128 = buffer.data(gg + 128);
    const auto *gg_129 = buffer.data(gg + 129);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_134 = buffer.data(gg + 134);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_136 = buffer.data(gg + 136);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_138 = buffer.data(gg + 138);
    const auto *gg_139 = buffer.data(gg + 139);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_141 = buffer.data(gg + 141);
    const auto *gg_142 = buffer.data(gg + 142);
    const auto *gg_143 = buffer.data(gg + 143);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_145 = buffer.data(gg + 145);
    const auto *gg_146 = buffer.data(gg + 146);
    const auto *gg_147 = buffer.data(gg + 147);
    const auto *gg_148 = buffer.data(gg + 148);
    const auto *gg_149 = buffer.data(gg + 149);
    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_151 = buffer.data(gg + 151);
    const auto *gg_152 = buffer.data(gg + 152);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_154 = buffer.data(gg + 154);
    const auto *gg_155 = buffer.data(gg + 155);
    const auto *gg_156 = buffer.data(gg + 156);
    const auto *gg_157 = buffer.data(gg + 157);
    const auto *gg_158 = buffer.data(gg + 158);
    const auto *gg_159 = buffer.data(gg + 159);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_161 = buffer.data(gg + 161);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_163 = buffer.data(gg + 163);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_165 = buffer.data(gg + 165);
    const auto *gg_166 = buffer.data(gg + 166);
    const auto *gg_167 = buffer.data(gg + 167);
    const auto *gg_168 = buffer.data(gg + 168);
    const auto *gg_169 = buffer.data(gg + 169);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_171 = buffer.data(gg + 171);
    const auto *gg_172 = buffer.data(gg + 172);
    const auto *gg_173 = buffer.data(gg + 173);
    const auto *gg_174 = buffer.data(gg + 174);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_181 = buffer.data(gg + 181);
    const auto *gg_182 = buffer.data(gg + 182);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_184 = buffer.data(gg + 184);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_186 = buffer.data(gg + 186);
    const auto *gg_187 = buffer.data(gg + 187);
    const auto *gg_188 = buffer.data(gg + 188);
    const auto *gg_189 = buffer.data(gg + 189);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_195 = buffer.data(gg + 195);
    const auto *gg_196 = buffer.data(gg + 196);
    const auto *gg_197 = buffer.data(gg + 197);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_199 = buffer.data(gg + 199);
    const auto *gg_200 = buffer.data(gg + 200);
    const auto *gg_201 = buffer.data(gg + 201);
    const auto *gg_202 = buffer.data(gg + 202);
    const auto *gg_203 = buffer.data(gg + 203);
    const auto *gg_204 = buffer.data(gg + 204);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_211 = buffer.data(gg + 211);
    const auto *gg_212 = buffer.data(gg + 212);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_214 = buffer.data(gg + 214);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_216 = buffer.data(gg + 216);
    const auto *gg_217 = buffer.data(gg + 217);
    const auto *gg_218 = buffer.data(gg + 218);
    const auto *gg_219 = buffer.data(gg + 219);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_221 = buffer.data(gg + 221);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_223 = buffer.data(gg + 223);
    const auto *gg_224 = buffer.data(gg + 224);

    const auto *ig_262 = buffer.data(ig + 262);
    const auto *ig_263 = buffer.data(ig + 263);
    const auto *ig_264 = buffer.data(ig + 264);
    const auto *ig_265 = buffer.data(ig + 265);
    const auto *ig_266 = buffer.data(ig + 266);
    const auto *ig_267 = buffer.data(ig + 267);
    const auto *ig_268 = buffer.data(ig + 268);
    const auto *ig_269 = buffer.data(ig + 269);
    const auto *ig_270 = buffer.data(ig + 270);
    const auto *ig_271 = buffer.data(ig + 271);
    const auto *ig_272 = buffer.data(ig + 272);
    const auto *ig_273 = buffer.data(ig + 273);
    const auto *ig_274 = buffer.data(ig + 274);
    const auto *ig_275 = buffer.data(ig + 275);
    const auto *ig_276 = buffer.data(ig + 276);
    const auto *ig_277 = buffer.data(ig + 277);
    const auto *ig_278 = buffer.data(ig + 278);
    const auto *ig_279 = buffer.data(ig + 279);
    const auto *ig_280 = buffer.data(ig + 280);
    const auto *ig_281 = buffer.data(ig + 281);
    const auto *ig_282 = buffer.data(ig + 282);
    const auto *ig_283 = buffer.data(ig + 283);
    const auto *ig_284 = buffer.data(ig + 284);
    const auto *ig_285 = buffer.data(ig + 285);
    const auto *ig_286 = buffer.data(ig + 286);
    const auto *ig_287 = buffer.data(ig + 287);
    const auto *ig_288 = buffer.data(ig + 288);
    const auto *ig_289 = buffer.data(ig + 289);
    const auto *ig_290 = buffer.data(ig + 290);
    const auto *ig_291 = buffer.data(ig + 291);
    const auto *ig_292 = buffer.data(ig + 292);
    const auto *ig_293 = buffer.data(ig + 293);
    const auto *ig_294 = buffer.data(ig + 294);
    const auto *ig_295 = buffer.data(ig + 295);
    const auto *ig_296 = buffer.data(ig + 296);
    const auto *ig_297 = buffer.data(ig + 297);
    const auto *ig_298 = buffer.data(ig + 298);
    const auto *ig_299 = buffer.data(ig + 299);
    const auto *ig_300 = buffer.data(ig + 300);
    const auto *ig_301 = buffer.data(ig + 301);
    const auto *ig_302 = buffer.data(ig + 302);
    const auto *ig_303 = buffer.data(ig + 303);
    const auto *ig_304 = buffer.data(ig + 304);
    const auto *ig_305 = buffer.data(ig + 305);
    const auto *ig_306 = buffer.data(ig + 306);
    const auto *ig_307 = buffer.data(ig + 307);
    const auto *ig_308 = buffer.data(ig + 308);
    const auto *ig_309 = buffer.data(ig + 309);
    const auto *ig_310 = buffer.data(ig + 310);
    const auto *ig_311 = buffer.data(ig + 311);
    const auto *ig_312 = buffer.data(ig + 312);
    const auto *ig_313 = buffer.data(ig + 313);
    const auto *ig_314 = buffer.data(ig + 314);
    const auto *ig_330 = buffer.data(ig + 330);
    const auto *ig_331 = buffer.data(ig + 331);
    const auto *ig_332 = buffer.data(ig + 332);
    const auto *ig_333 = buffer.data(ig + 333);
    const auto *ig_334 = buffer.data(ig + 334);
    const auto *ig_335 = buffer.data(ig + 335);
    const auto *ig_336 = buffer.data(ig + 336);
    const auto *ig_337 = buffer.data(ig + 337);
    const auto *ig_338 = buffer.data(ig + 338);
    const auto *ig_339 = buffer.data(ig + 339);
    const auto *ig_340 = buffer.data(ig + 340);
    const auto *ig_341 = buffer.data(ig + 341);
    const auto *ig_342 = buffer.data(ig + 342);
    const auto *ig_343 = buffer.data(ig + 343);
    const auto *ig_344 = buffer.data(ig + 344);
    const auto *ig_345 = buffer.data(ig + 345);
    const auto *ig_346 = buffer.data(ig + 346);
    const auto *ig_347 = buffer.data(ig + 347);
    const auto *ig_348 = buffer.data(ig + 348);
    const auto *ig_349 = buffer.data(ig + 349);
    const auto *ig_350 = buffer.data(ig + 350);
    const auto *ig_351 = buffer.data(ig + 351);
    const auto *ig_352 = buffer.data(ig + 352);
    const auto *ig_353 = buffer.data(ig + 353);
    const auto *ig_354 = buffer.data(ig + 354);
    const auto *ig_355 = buffer.data(ig + 355);
    const auto *ig_356 = buffer.data(ig + 356);
    const auto *ig_357 = buffer.data(ig + 357);
    const auto *ig_358 = buffer.data(ig + 358);
    const auto *ig_359 = buffer.data(ig + 359);
    const auto *ig_360 = buffer.data(ig + 360);
    const auto *ig_361 = buffer.data(ig + 361);
    const auto *ig_362 = buffer.data(ig + 362);
    const auto *ig_363 = buffer.data(ig + 363);
    const auto *ig_364 = buffer.data(ig + 364);
    const auto *ig_365 = buffer.data(ig + 365);
    const auto *ig_366 = buffer.data(ig + 366);
    const auto *ig_367 = buffer.data(ig + 367);
    const auto *ig_368 = buffer.data(ig + 368);
    const auto *ig_369 = buffer.data(ig + 369);
    const auto *ig_370 = buffer.data(ig + 370);
    const auto *ig_371 = buffer.data(ig + 371);
    const auto *ig_372 = buffer.data(ig + 372);
    const auto *ig_373 = buffer.data(ig + 373);
    const auto *ig_374 = buffer.data(ig + 374);
    const auto *ig_375 = buffer.data(ig + 375);
    const auto *ig_376 = buffer.data(ig + 376);
    const auto *ig_377 = buffer.data(ig + 377);
    const auto *ig_378 = buffer.data(ig + 378);
    const auto *ig_379 = buffer.data(ig + 379);
    const auto *ig_380 = buffer.data(ig + 380);
    const auto *ig_381 = buffer.data(ig + 381);
    const auto *ig_382 = buffer.data(ig + 382);
    const auto *ig_383 = buffer.data(ig + 383);
    const auto *ig_384 = buffer.data(ig + 384);
    const auto *ig_385 = buffer.data(ig + 385);
    const auto *ig_386 = buffer.data(ig + 386);
    const auto *ig_387 = buffer.data(ig + 387);
    const auto *ig_388 = buffer.data(ig + 388);
    const auto *ig_389 = buffer.data(ig + 389);
    const auto *ig_390 = buffer.data(ig + 390);
    const auto *ig_391 = buffer.data(ig + 391);
    const auto *ig_392 = buffer.data(ig + 392);
    const auto *ig_393 = buffer.data(ig + 393);
    const auto *ig_394 = buffer.data(ig + 394);
    const auto *ig_395 = buffer.data(ig + 395);
    const auto *ig_396 = buffer.data(ig + 396);
    const auto *ig_397 = buffer.data(ig + 397);
    const auto *ig_398 = buffer.data(ig + 398);
    const auto *ig_399 = buffer.data(ig + 399);
    const auto *ig_400 = buffer.data(ig + 400);
    const auto *ig_401 = buffer.data(ig + 401);
    const auto *ig_402 = buffer.data(ig + 402);
    const auto *ig_403 = buffer.data(ig + 403);
    const auto *ig_404 = buffer.data(ig + 404);
    const auto *ig_405 = buffer.data(ig + 405);
    const auto *ig_406 = buffer.data(ig + 406);
    const auto *ig_407 = buffer.data(ig + 407);
    const auto *ig_408 = buffer.data(ig + 408);
    const auto *ig_409 = buffer.data(ig + 409);
    const auto *ig_410 = buffer.data(ig + 410);
    const auto *ig_411 = buffer.data(ig + 411);
    const auto *ig_412 = buffer.data(ig + 412);
    const auto *ig_413 = buffer.data(ig + 413);
    const auto *ig_414 = buffer.data(ig + 414);
    const auto *ig_415 = buffer.data(ig + 415);
    const auto *ig_416 = buffer.data(ig + 416);
    const auto *ig_417 = buffer.data(ig + 417);
    const auto *ig_418 = buffer.data(ig + 418);
    const auto *ig_419 = buffer.data(ig + 419);

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, gg_97, gg_98, gg_99, gg_100, \
                         gg_101, ig_262, ig_263, ig_264, ig_265, \
                         ig_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -gg_97[k]
                   + f_0 * ig_262[k];

        t_173[k] = -gg_98[k]
                   + f_0 * ig_263[k];

        t_174[k] = -gg_99[k]
                   + f_0 * ig_264[k];

        t_175[k] = -gg_100[k]
                   + f_0 * ig_265[k];

        t_176[k] = -gg_101[k]
                   + f_0 * ig_266[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, gg_102, gg_103, gg_104, gg_105, \
                         gg_106, ig_267, ig_268, ig_269, ig_270, \
                         ig_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -gg_102[k]
                   + f_0 * ig_267[k];

        t_178[k] = -gg_103[k]
                   + f_0 * ig_268[k];

        t_179[k] = -gg_104[k]
                   + f_0 * ig_269[k];

        t_180[k] = -2.0 * gg_105[k]
                   + f_0 * ig_270[k];

        t_181[k] = -2.0 * gg_106[k]
                   + f_0 * ig_271[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, gg_107, gg_108, gg_109, gg_110, \
                         gg_111, ig_272, ig_273, ig_274, ig_275, \
                         ig_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * gg_107[k]
                   + f_0 * ig_272[k];

        t_183[k] = -2.0 * gg_108[k]
                   + f_0 * ig_273[k];

        t_184[k] = -2.0 * gg_109[k]
                   + f_0 * ig_274[k];

        t_185[k] = -2.0 * gg_110[k]
                   + f_0 * ig_275[k];

        t_186[k] = -2.0 * gg_111[k]
                   + f_0 * ig_276[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, gg_112, gg_113, gg_114, gg_115, \
                         gg_116, ig_277, ig_278, ig_279, ig_280, \
                         ig_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * gg_112[k]
                   + f_0 * ig_277[k];

        t_188[k] = -2.0 * gg_113[k]
                   + f_0 * ig_278[k];

        t_189[k] = -2.0 * gg_114[k]
                   + f_0 * ig_279[k];

        t_190[k] = -2.0 * gg_115[k]
                   + f_0 * ig_280[k];

        t_191[k] = -2.0 * gg_116[k]
                   + f_0 * ig_281[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, gg_117, gg_118, gg_119, gg_120, \
                         gg_121, ig_282, ig_283, ig_284, ig_285, \
                         ig_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -2.0 * gg_117[k]
                   + f_0 * ig_282[k];

        t_193[k] = -2.0 * gg_118[k]
                   + f_0 * ig_283[k];

        t_194[k] = -2.0 * gg_119[k]
                   + f_0 * ig_284[k];

        t_195[k] = -3.0 * gg_120[k]
                   + f_0 * ig_285[k];

        t_196[k] = -3.0 * gg_121[k]
                   + f_0 * ig_286[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, gg_122, gg_123, gg_124, gg_125, \
                         gg_126, ig_287, ig_288, ig_289, ig_290, \
                         ig_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -3.0 * gg_122[k]
                   + f_0 * ig_287[k];

        t_198[k] = -3.0 * gg_123[k]
                   + f_0 * ig_288[k];

        t_199[k] = -3.0 * gg_124[k]
                   + f_0 * ig_289[k];

        t_200[k] = -3.0 * gg_125[k]
                   + f_0 * ig_290[k];

        t_201[k] = -3.0 * gg_126[k]
                   + f_0 * ig_291[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, gg_127, gg_128, gg_129, gg_130, \
                         gg_131, ig_292, ig_293, ig_294, ig_295, \
                         ig_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -3.0 * gg_127[k]
                   + f_0 * ig_292[k];

        t_203[k] = -3.0 * gg_128[k]
                   + f_0 * ig_293[k];

        t_204[k] = -3.0 * gg_129[k]
                   + f_0 * ig_294[k];

        t_205[k] = -3.0 * gg_130[k]
                   + f_0 * ig_295[k];

        t_206[k] = -3.0 * gg_131[k]
                   + f_0 * ig_296[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, gg_132, gg_133, gg_134, gg_135, \
                         gg_136, ig_297, ig_298, ig_299, ig_300, \
                         ig_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -3.0 * gg_132[k]
                   + f_0 * ig_297[k];

        t_208[k] = -3.0 * gg_133[k]
                   + f_0 * ig_298[k];

        t_209[k] = -3.0 * gg_134[k]
                   + f_0 * ig_299[k];

        t_210[k] = -4.0 * gg_135[k]
                   + f_0 * ig_300[k];

        t_211[k] = -4.0 * gg_136[k]
                   + f_0 * ig_301[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, gg_137, gg_138, gg_139, gg_140, \
                         gg_141, ig_302, ig_303, ig_304, ig_305, \
                         ig_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = -4.0 * gg_137[k]
                   + f_0 * ig_302[k];

        t_213[k] = -4.0 * gg_138[k]
                   + f_0 * ig_303[k];

        t_214[k] = -4.0 * gg_139[k]
                   + f_0 * ig_304[k];

        t_215[k] = -4.0 * gg_140[k]
                   + f_0 * ig_305[k];

        t_216[k] = -4.0 * gg_141[k]
                   + f_0 * ig_306[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, gg_142, gg_143, gg_144, gg_145, \
                         gg_146, ig_307, ig_308, ig_309, ig_310, \
                         ig_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = -4.0 * gg_142[k]
                   + f_0 * ig_307[k];

        t_218[k] = -4.0 * gg_143[k]
                   + f_0 * ig_308[k];

        t_219[k] = -4.0 * gg_144[k]
                   + f_0 * ig_309[k];

        t_220[k] = -4.0 * gg_145[k]
                   + f_0 * ig_310[k];

        t_221[k] = -4.0 * gg_146[k]
                   + f_0 * ig_311[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, gg_147, gg_148, gg_149, \
                         ig_312, ig_313, ig_314, ig_330, ig_331, \
                         ig_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -4.0 * gg_147[k]
                   + f_0 * ig_312[k];

        t_223[k] = -4.0 * gg_148[k]
                   + f_0 * ig_313[k];

        t_224[k] = -4.0 * gg_149[k]
                   + f_0 * ig_314[k];

        t_225[k] = f_0 * ig_330[k];

        t_226[k] = f_0 * ig_331[k];

        t_227[k] = f_0 * ig_332[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, t_233, t_234, t_235, ig_333, \
                         ig_334, ig_335, ig_336, ig_337, ig_338, ig_339, \
                         ig_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_0 * ig_333[k];

        t_229[k] = f_0 * ig_334[k];

        t_230[k] = f_0 * ig_335[k];

        t_231[k] = f_0 * ig_336[k];

        t_232[k] = f_0 * ig_337[k];

        t_233[k] = f_0 * ig_338[k];

        t_234[k] = f_0 * ig_339[k];

        t_235[k] = f_0 * ig_340[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, t_241, gg_150, gg_151, ig_341, \
                         ig_342, ig_343, ig_344, ig_345, ig_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_0 * ig_341[k];

        t_237[k] = f_0 * ig_342[k];

        t_238[k] = f_0 * ig_343[k];

        t_239[k] = f_0 * ig_344[k];

        t_240[k] = -gg_150[k]
                   + f_0 * ig_345[k];

        t_241[k] = -gg_151[k]
                   + f_0 * ig_346[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, gg_152, gg_153, gg_154, gg_155, \
                         gg_156, ig_347, ig_348, ig_349, ig_350, \
                         ig_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = -gg_152[k]
                   + f_0 * ig_347[k];

        t_243[k] = -gg_153[k]
                   + f_0 * ig_348[k];

        t_244[k] = -gg_154[k]
                   + f_0 * ig_349[k];

        t_245[k] = -gg_155[k]
                   + f_0 * ig_350[k];

        t_246[k] = -gg_156[k]
                   + f_0 * ig_351[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, gg_157, gg_158, gg_159, gg_160, \
                         gg_161, ig_352, ig_353, ig_354, ig_355, \
                         ig_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = -gg_157[k]
                   + f_0 * ig_352[k];

        t_248[k] = -gg_158[k]
                   + f_0 * ig_353[k];

        t_249[k] = -gg_159[k]
                   + f_0 * ig_354[k];

        t_250[k] = -gg_160[k]
                   + f_0 * ig_355[k];

        t_251[k] = -gg_161[k]
                   + f_0 * ig_356[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, gg_162, gg_163, gg_164, gg_165, \
                         gg_166, ig_357, ig_358, ig_359, ig_360, \
                         ig_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -gg_162[k]
                   + f_0 * ig_357[k];

        t_253[k] = -gg_163[k]
                   + f_0 * ig_358[k];

        t_254[k] = -gg_164[k]
                   + f_0 * ig_359[k];

        t_255[k] = -2.0 * gg_165[k]
                   + f_0 * ig_360[k];

        t_256[k] = -2.0 * gg_166[k]
                   + f_0 * ig_361[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, gg_167, gg_168, gg_169, gg_170, \
                         gg_171, ig_362, ig_363, ig_364, ig_365, \
                         ig_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = -2.0 * gg_167[k]
                   + f_0 * ig_362[k];

        t_258[k] = -2.0 * gg_168[k]
                   + f_0 * ig_363[k];

        t_259[k] = -2.0 * gg_169[k]
                   + f_0 * ig_364[k];

        t_260[k] = -2.0 * gg_170[k]
                   + f_0 * ig_365[k];

        t_261[k] = -2.0 * gg_171[k]
                   + f_0 * ig_366[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, gg_172, gg_173, gg_174, gg_175, \
                         gg_176, ig_367, ig_368, ig_369, ig_370, \
                         ig_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = -2.0 * gg_172[k]
                   + f_0 * ig_367[k];

        t_263[k] = -2.0 * gg_173[k]
                   + f_0 * ig_368[k];

        t_264[k] = -2.0 * gg_174[k]
                   + f_0 * ig_369[k];

        t_265[k] = -2.0 * gg_175[k]
                   + f_0 * ig_370[k];

        t_266[k] = -2.0 * gg_176[k]
                   + f_0 * ig_371[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, gg_177, gg_178, gg_179, gg_180, \
                         gg_181, ig_372, ig_373, ig_374, ig_375, \
                         ig_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = -2.0 * gg_177[k]
                   + f_0 * ig_372[k];

        t_268[k] = -2.0 * gg_178[k]
                   + f_0 * ig_373[k];

        t_269[k] = -2.0 * gg_179[k]
                   + f_0 * ig_374[k];

        t_270[k] = -3.0 * gg_180[k]
                   + f_0 * ig_375[k];

        t_271[k] = -3.0 * gg_181[k]
                   + f_0 * ig_376[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, gg_182, gg_183, gg_184, gg_185, \
                         gg_186, ig_377, ig_378, ig_379, ig_380, \
                         ig_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -3.0 * gg_182[k]
                   + f_0 * ig_377[k];

        t_273[k] = -3.0 * gg_183[k]
                   + f_0 * ig_378[k];

        t_274[k] = -3.0 * gg_184[k]
                   + f_0 * ig_379[k];

        t_275[k] = -3.0 * gg_185[k]
                   + f_0 * ig_380[k];

        t_276[k] = -3.0 * gg_186[k]
                   + f_0 * ig_381[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, t_281, gg_187, gg_188, gg_189, gg_190, \
                         gg_191, ig_382, ig_383, ig_384, ig_385, \
                         ig_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = -3.0 * gg_187[k]
                   + f_0 * ig_382[k];

        t_278[k] = -3.0 * gg_188[k]
                   + f_0 * ig_383[k];

        t_279[k] = -3.0 * gg_189[k]
                   + f_0 * ig_384[k];

        t_280[k] = -3.0 * gg_190[k]
                   + f_0 * ig_385[k];

        t_281[k] = -3.0 * gg_191[k]
                   + f_0 * ig_386[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, gg_192, gg_193, gg_194, gg_195, \
                         gg_196, ig_387, ig_388, ig_389, ig_390, \
                         ig_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = -3.0 * gg_192[k]
                   + f_0 * ig_387[k];

        t_283[k] = -3.0 * gg_193[k]
                   + f_0 * ig_388[k];

        t_284[k] = -3.0 * gg_194[k]
                   + f_0 * ig_389[k];

        t_285[k] = -4.0 * gg_195[k]
                   + f_0 * ig_390[k];

        t_286[k] = -4.0 * gg_196[k]
                   + f_0 * ig_391[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, gg_197, gg_198, gg_199, gg_200, \
                         gg_201, ig_392, ig_393, ig_394, ig_395, \
                         ig_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = -4.0 * gg_197[k]
                   + f_0 * ig_392[k];

        t_288[k] = -4.0 * gg_198[k]
                   + f_0 * ig_393[k];

        t_289[k] = -4.0 * gg_199[k]
                   + f_0 * ig_394[k];

        t_290[k] = -4.0 * gg_200[k]
                   + f_0 * ig_395[k];

        t_291[k] = -4.0 * gg_201[k]
                   + f_0 * ig_396[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, gg_202, gg_203, gg_204, gg_205, \
                         gg_206, ig_397, ig_398, ig_399, ig_400, \
                         ig_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = -4.0 * gg_202[k]
                   + f_0 * ig_397[k];

        t_293[k] = -4.0 * gg_203[k]
                   + f_0 * ig_398[k];

        t_294[k] = -4.0 * gg_204[k]
                   + f_0 * ig_399[k];

        t_295[k] = -4.0 * gg_205[k]
                   + f_0 * ig_400[k];

        t_296[k] = -4.0 * gg_206[k]
                   + f_0 * ig_401[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, gg_207, gg_208, gg_209, gg_210, \
                         gg_211, ig_402, ig_403, ig_404, ig_405, \
                         ig_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = -4.0 * gg_207[k]
                   + f_0 * ig_402[k];

        t_298[k] = -4.0 * gg_208[k]
                   + f_0 * ig_403[k];

        t_299[k] = -4.0 * gg_209[k]
                   + f_0 * ig_404[k];

        t_300[k] = -5.0 * gg_210[k]
                   + f_0 * ig_405[k];

        t_301[k] = -5.0 * gg_211[k]
                   + f_0 * ig_406[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, gg_212, gg_213, gg_214, gg_215, \
                         gg_216, ig_407, ig_408, ig_409, ig_410, \
                         ig_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = -5.0 * gg_212[k]
                   + f_0 * ig_407[k];

        t_303[k] = -5.0 * gg_213[k]
                   + f_0 * ig_408[k];

        t_304[k] = -5.0 * gg_214[k]
                   + f_0 * ig_409[k];

        t_305[k] = -5.0 * gg_215[k]
                   + f_0 * ig_410[k];

        t_306[k] = -5.0 * gg_216[k]
                   + f_0 * ig_411[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, gg_217, gg_218, gg_219, gg_220, \
                         gg_221, ig_412, ig_413, ig_414, ig_415, \
                         ig_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = -5.0 * gg_217[k]
                   + f_0 * ig_412[k];

        t_308[k] = -5.0 * gg_218[k]
                   + f_0 * ig_413[k];

        t_309[k] = -5.0 * gg_219[k]
                   + f_0 * ig_414[k];

        t_310[k] = -5.0 * gg_220[k]
                   + f_0 * ig_415[k];

        t_311[k] = -5.0 * gg_221[k]
                   + f_0 * ig_416[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, gg_222, gg_223, gg_224, ig_417, ig_418, \
                         ig_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -5.0 * gg_222[k]
                   + f_0 * ig_417[k];

        t_313[k] = -5.0 * gg_223[k]
                   + f_0 * ig_418[k];

        t_314[k] = -5.0 * gg_224[k]
                   + f_0 * ig_419[k];
    }
}

auto
compute_prim_geom_10_hg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t gg, const size_t ig,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hg_electron_repulsion_2_piece0(buffer, target, gg, ig, ncols, alpha);

    compute_prim_geom_10_hg_electron_repulsion_2_piece1(buffer, target, gg, ig, ncols, alpha);
}

}  // namespace simdt2ceri
