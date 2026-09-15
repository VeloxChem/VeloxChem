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


#include "SimdElectronRepulsionGeom10VrrRecIG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_ig_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hg, const size_t kg,
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

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_1 = buffer.data(kg + 1);
    const auto *kg_2 = buffer.data(kg + 2);
    const auto *kg_3 = buffer.data(kg + 3);
    const auto *kg_4 = buffer.data(kg + 4);
    const auto *kg_5 = buffer.data(kg + 5);
    const auto *kg_6 = buffer.data(kg + 6);
    const auto *kg_7 = buffer.data(kg + 7);
    const auto *kg_8 = buffer.data(kg + 8);
    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_11 = buffer.data(kg + 11);
    const auto *kg_12 = buffer.data(kg + 12);
    const auto *kg_13 = buffer.data(kg + 13);
    const auto *kg_14 = buffer.data(kg + 14);
    const auto *kg_15 = buffer.data(kg + 15);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_17 = buffer.data(kg + 17);
    const auto *kg_18 = buffer.data(kg + 18);
    const auto *kg_19 = buffer.data(kg + 19);
    const auto *kg_20 = buffer.data(kg + 20);
    const auto *kg_21 = buffer.data(kg + 21);
    const auto *kg_22 = buffer.data(kg + 22);
    const auto *kg_23 = buffer.data(kg + 23);
    const auto *kg_24 = buffer.data(kg + 24);
    const auto *kg_25 = buffer.data(kg + 25);
    const auto *kg_26 = buffer.data(kg + 26);
    const auto *kg_27 = buffer.data(kg + 27);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_30 = buffer.data(kg + 30);
    const auto *kg_31 = buffer.data(kg + 31);
    const auto *kg_32 = buffer.data(kg + 32);
    const auto *kg_33 = buffer.data(kg + 33);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_35 = buffer.data(kg + 35);
    const auto *kg_36 = buffer.data(kg + 36);
    const auto *kg_37 = buffer.data(kg + 37);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_39 = buffer.data(kg + 39);
    const auto *kg_40 = buffer.data(kg + 40);
    const auto *kg_41 = buffer.data(kg + 41);
    const auto *kg_42 = buffer.data(kg + 42);
    const auto *kg_43 = buffer.data(kg + 43);
    const auto *kg_44 = buffer.data(kg + 44);
    const auto *kg_45 = buffer.data(kg + 45);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_48 = buffer.data(kg + 48);
    const auto *kg_49 = buffer.data(kg + 49);
    const auto *kg_50 = buffer.data(kg + 50);
    const auto *kg_51 = buffer.data(kg + 51);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_53 = buffer.data(kg + 53);
    const auto *kg_54 = buffer.data(kg + 54);
    const auto *kg_55 = buffer.data(kg + 55);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_58 = buffer.data(kg + 58);
    const auto *kg_59 = buffer.data(kg + 59);
    const auto *kg_60 = buffer.data(kg + 60);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_62 = buffer.data(kg + 62);
    const auto *kg_63 = buffer.data(kg + 63);
    const auto *kg_64 = buffer.data(kg + 64);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_67 = buffer.data(kg + 67);
    const auto *kg_68 = buffer.data(kg + 68);
    const auto *kg_69 = buffer.data(kg + 69);
    const auto *kg_70 = buffer.data(kg + 70);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_72 = buffer.data(kg + 72);
    const auto *kg_73 = buffer.data(kg + 73);
    const auto *kg_74 = buffer.data(kg + 74);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_78 = buffer.data(kg + 78);
    const auto *kg_79 = buffer.data(kg + 79);
    const auto *kg_80 = buffer.data(kg + 80);
    const auto *kg_81 = buffer.data(kg + 81);
    const auto *kg_82 = buffer.data(kg + 82);
    const auto *kg_83 = buffer.data(kg + 83);
    const auto *kg_84 = buffer.data(kg + 84);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_86 = buffer.data(kg + 86);
    const auto *kg_87 = buffer.data(kg + 87);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_91 = buffer.data(kg + 91);
    const auto *kg_92 = buffer.data(kg + 92);
    const auto *kg_93 = buffer.data(kg + 93);
    const auto *kg_94 = buffer.data(kg + 94);
    const auto *kg_95 = buffer.data(kg + 95);
    const auto *kg_96 = buffer.data(kg + 96);
    const auto *kg_97 = buffer.data(kg + 97);
    const auto *kg_98 = buffer.data(kg + 98);
    const auto *kg_99 = buffer.data(kg + 99);
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_101 = buffer.data(kg + 101);
    const auto *kg_102 = buffer.data(kg + 102);
    const auto *kg_103 = buffer.data(kg + 103);
    const auto *kg_104 = buffer.data(kg + 104);
    const auto *kg_105 = buffer.data(kg + 105);
    const auto *kg_106 = buffer.data(kg + 106);
    const auto *kg_107 = buffer.data(kg + 107);
    const auto *kg_108 = buffer.data(kg + 108);
    const auto *kg_109 = buffer.data(kg + 109);
    const auto *kg_110 = buffer.data(kg + 110);
    const auto *kg_111 = buffer.data(kg + 111);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_114 = buffer.data(kg + 114);
    const auto *kg_115 = buffer.data(kg + 115);
    const auto *kg_116 = buffer.data(kg + 116);
    const auto *kg_117 = buffer.data(kg + 117);
    const auto *kg_118 = buffer.data(kg + 118);
    const auto *kg_119 = buffer.data(kg + 119);
    const auto *kg_120 = buffer.data(kg + 120);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_123 = buffer.data(kg + 123);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_125 = buffer.data(kg + 125);
    const auto *kg_126 = buffer.data(kg + 126);
    const auto *kg_127 = buffer.data(kg + 127);
    const auto *kg_128 = buffer.data(kg + 128);
    const auto *kg_129 = buffer.data(kg + 129);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_132 = buffer.data(kg + 132);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_134 = buffer.data(kg + 134);
    const auto *kg_135 = buffer.data(kg + 135);
    const auto *kg_136 = buffer.data(kg + 136);
    const auto *kg_137 = buffer.data(kg + 137);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_139 = buffer.data(kg + 139);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_141 = buffer.data(kg + 141);
    const auto *kg_142 = buffer.data(kg + 142);
    const auto *kg_143 = buffer.data(kg + 143);
    const auto *kg_144 = buffer.data(kg + 144);
    const auto *kg_145 = buffer.data(kg + 145);
    const auto *kg_146 = buffer.data(kg + 146);
    const auto *kg_147 = buffer.data(kg + 147);
    const auto *kg_148 = buffer.data(kg + 148);
    const auto *kg_149 = buffer.data(kg + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, hg_0, hg_1, hg_2, hg_3, hg_4, kg_0, kg_1, \
                         kg_2, kg_3, kg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -6.0 * hg_0[k]
                 + f_0 * kg_0[k];

        t_1[k] = -6.0 * hg_1[k]
                 + f_0 * kg_1[k];

        t_2[k] = -6.0 * hg_2[k]
                 + f_0 * kg_2[k];

        t_3[k] = -6.0 * hg_3[k]
                 + f_0 * kg_3[k];

        t_4[k] = -6.0 * hg_4[k]
                 + f_0 * kg_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, hg_5, hg_6, hg_7, hg_8, hg_9, kg_5, kg_6, \
                         kg_7, kg_8, kg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -6.0 * hg_5[k]
                 + f_0 * kg_5[k];

        t_6[k] = -6.0 * hg_6[k]
                 + f_0 * kg_6[k];

        t_7[k] = -6.0 * hg_7[k]
                 + f_0 * kg_7[k];

        t_8[k] = -6.0 * hg_8[k]
                 + f_0 * kg_8[k];

        t_9[k] = -6.0 * hg_9[k]
                 + f_0 * kg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, hg_10, hg_11, hg_12, hg_13, hg_14, \
                         kg_10, kg_11, kg_12, kg_13, kg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -6.0 * hg_10[k]
                  + f_0 * kg_10[k];

        t_11[k] = -6.0 * hg_11[k]
                  + f_0 * kg_11[k];

        t_12[k] = -6.0 * hg_12[k]
                  + f_0 * kg_12[k];

        t_13[k] = -6.0 * hg_13[k]
                  + f_0 * kg_13[k];

        t_14[k] = -6.0 * hg_14[k]
                  + f_0 * kg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, hg_15, hg_16, hg_17, hg_18, hg_19, \
                         kg_15, kg_16, kg_17, kg_18, kg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -5.0 * hg_15[k]
                  + f_0 * kg_15[k];

        t_16[k] = -5.0 * hg_16[k]
                  + f_0 * kg_16[k];

        t_17[k] = -5.0 * hg_17[k]
                  + f_0 * kg_17[k];

        t_18[k] = -5.0 * hg_18[k]
                  + f_0 * kg_18[k];

        t_19[k] = -5.0 * hg_19[k]
                  + f_0 * kg_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, hg_20, hg_21, hg_22, hg_23, hg_24, \
                         kg_20, kg_21, kg_22, kg_23, kg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -5.0 * hg_20[k]
                  + f_0 * kg_20[k];

        t_21[k] = -5.0 * hg_21[k]
                  + f_0 * kg_21[k];

        t_22[k] = -5.0 * hg_22[k]
                  + f_0 * kg_22[k];

        t_23[k] = -5.0 * hg_23[k]
                  + f_0 * kg_23[k];

        t_24[k] = -5.0 * hg_24[k]
                  + f_0 * kg_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, hg_25, hg_26, hg_27, hg_28, hg_29, \
                         kg_25, kg_26, kg_27, kg_28, kg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -5.0 * hg_25[k]
                  + f_0 * kg_25[k];

        t_26[k] = -5.0 * hg_26[k]
                  + f_0 * kg_26[k];

        t_27[k] = -5.0 * hg_27[k]
                  + f_0 * kg_27[k];

        t_28[k] = -5.0 * hg_28[k]
                  + f_0 * kg_28[k];

        t_29[k] = -5.0 * hg_29[k]
                  + f_0 * kg_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, hg_30, hg_31, hg_32, hg_33, hg_34, \
                         kg_30, kg_31, kg_32, kg_33, kg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -5.0 * hg_30[k]
                  + f_0 * kg_30[k];

        t_31[k] = -5.0 * hg_31[k]
                  + f_0 * kg_31[k];

        t_32[k] = -5.0 * hg_32[k]
                  + f_0 * kg_32[k];

        t_33[k] = -5.0 * hg_33[k]
                  + f_0 * kg_33[k];

        t_34[k] = -5.0 * hg_34[k]
                  + f_0 * kg_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, hg_35, hg_36, hg_37, hg_38, hg_39, \
                         kg_35, kg_36, kg_37, kg_38, kg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -5.0 * hg_35[k]
                  + f_0 * kg_35[k];

        t_36[k] = -5.0 * hg_36[k]
                  + f_0 * kg_36[k];

        t_37[k] = -5.0 * hg_37[k]
                  + f_0 * kg_37[k];

        t_38[k] = -5.0 * hg_38[k]
                  + f_0 * kg_38[k];

        t_39[k] = -5.0 * hg_39[k]
                  + f_0 * kg_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, hg_40, hg_41, hg_42, hg_43, hg_44, \
                         kg_40, kg_41, kg_42, kg_43, kg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -5.0 * hg_40[k]
                  + f_0 * kg_40[k];

        t_41[k] = -5.0 * hg_41[k]
                  + f_0 * kg_41[k];

        t_42[k] = -5.0 * hg_42[k]
                  + f_0 * kg_42[k];

        t_43[k] = -5.0 * hg_43[k]
                  + f_0 * kg_43[k];

        t_44[k] = -5.0 * hg_44[k]
                  + f_0 * kg_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hg_45, hg_46, hg_47, hg_48, hg_49, \
                         kg_45, kg_46, kg_47, kg_48, kg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -4.0 * hg_45[k]
                  + f_0 * kg_45[k];

        t_46[k] = -4.0 * hg_46[k]
                  + f_0 * kg_46[k];

        t_47[k] = -4.0 * hg_47[k]
                  + f_0 * kg_47[k];

        t_48[k] = -4.0 * hg_48[k]
                  + f_0 * kg_48[k];

        t_49[k] = -4.0 * hg_49[k]
                  + f_0 * kg_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hg_50, hg_51, hg_52, hg_53, hg_54, \
                         kg_50, kg_51, kg_52, kg_53, kg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -4.0 * hg_50[k]
                  + f_0 * kg_50[k];

        t_51[k] = -4.0 * hg_51[k]
                  + f_0 * kg_51[k];

        t_52[k] = -4.0 * hg_52[k]
                  + f_0 * kg_52[k];

        t_53[k] = -4.0 * hg_53[k]
                  + f_0 * kg_53[k];

        t_54[k] = -4.0 * hg_54[k]
                  + f_0 * kg_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, hg_55, hg_56, hg_57, hg_58, hg_59, \
                         kg_55, kg_56, kg_57, kg_58, kg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -4.0 * hg_55[k]
                  + f_0 * kg_55[k];

        t_56[k] = -4.0 * hg_56[k]
                  + f_0 * kg_56[k];

        t_57[k] = -4.0 * hg_57[k]
                  + f_0 * kg_57[k];

        t_58[k] = -4.0 * hg_58[k]
                  + f_0 * kg_58[k];

        t_59[k] = -4.0 * hg_59[k]
                  + f_0 * kg_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, hg_60, hg_61, hg_62, hg_63, hg_64, \
                         kg_60, kg_61, kg_62, kg_63, kg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -4.0 * hg_60[k]
                  + f_0 * kg_60[k];

        t_61[k] = -4.0 * hg_61[k]
                  + f_0 * kg_61[k];

        t_62[k] = -4.0 * hg_62[k]
                  + f_0 * kg_62[k];

        t_63[k] = -4.0 * hg_63[k]
                  + f_0 * kg_63[k];

        t_64[k] = -4.0 * hg_64[k]
                  + f_0 * kg_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, hg_65, hg_66, hg_67, hg_68, hg_69, \
                         kg_65, kg_66, kg_67, kg_68, kg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -4.0 * hg_65[k]
                  + f_0 * kg_65[k];

        t_66[k] = -4.0 * hg_66[k]
                  + f_0 * kg_66[k];

        t_67[k] = -4.0 * hg_67[k]
                  + f_0 * kg_67[k];

        t_68[k] = -4.0 * hg_68[k]
                  + f_0 * kg_68[k];

        t_69[k] = -4.0 * hg_69[k]
                  + f_0 * kg_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, hg_70, hg_71, hg_72, hg_73, hg_74, \
                         kg_70, kg_71, kg_72, kg_73, kg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -4.0 * hg_70[k]
                  + f_0 * kg_70[k];

        t_71[k] = -4.0 * hg_71[k]
                  + f_0 * kg_71[k];

        t_72[k] = -4.0 * hg_72[k]
                  + f_0 * kg_72[k];

        t_73[k] = -4.0 * hg_73[k]
                  + f_0 * kg_73[k];

        t_74[k] = -4.0 * hg_74[k]
                  + f_0 * kg_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, hg_75, hg_76, hg_77, hg_78, hg_79, \
                         kg_75, kg_76, kg_77, kg_78, kg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -4.0 * hg_75[k]
                  + f_0 * kg_75[k];

        t_76[k] = -4.0 * hg_76[k]
                  + f_0 * kg_76[k];

        t_77[k] = -4.0 * hg_77[k]
                  + f_0 * kg_77[k];

        t_78[k] = -4.0 * hg_78[k]
                  + f_0 * kg_78[k];

        t_79[k] = -4.0 * hg_79[k]
                  + f_0 * kg_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, hg_80, hg_81, hg_82, hg_83, hg_84, \
                         kg_80, kg_81, kg_82, kg_83, kg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -4.0 * hg_80[k]
                  + f_0 * kg_80[k];

        t_81[k] = -4.0 * hg_81[k]
                  + f_0 * kg_81[k];

        t_82[k] = -4.0 * hg_82[k]
                  + f_0 * kg_82[k];

        t_83[k] = -4.0 * hg_83[k]
                  + f_0 * kg_83[k];

        t_84[k] = -4.0 * hg_84[k]
                  + f_0 * kg_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, hg_85, hg_86, hg_87, hg_88, hg_89, \
                         kg_85, kg_86, kg_87, kg_88, kg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * hg_85[k]
                  + f_0 * kg_85[k];

        t_86[k] = -4.0 * hg_86[k]
                  + f_0 * kg_86[k];

        t_87[k] = -4.0 * hg_87[k]
                  + f_0 * kg_87[k];

        t_88[k] = -4.0 * hg_88[k]
                  + f_0 * kg_88[k];

        t_89[k] = -4.0 * hg_89[k]
                  + f_0 * kg_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, hg_90, hg_91, hg_92, hg_93, hg_94, \
                         kg_90, kg_91, kg_92, kg_93, kg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * hg_90[k]
                  + f_0 * kg_90[k];

        t_91[k] = -3.0 * hg_91[k]
                  + f_0 * kg_91[k];

        t_92[k] = -3.0 * hg_92[k]
                  + f_0 * kg_92[k];

        t_93[k] = -3.0 * hg_93[k]
                  + f_0 * kg_93[k];

        t_94[k] = -3.0 * hg_94[k]
                  + f_0 * kg_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, hg_95, hg_96, hg_97, hg_98, hg_99, \
                         kg_95, kg_96, kg_97, kg_98, kg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * hg_95[k]
                  + f_0 * kg_95[k];

        t_96[k] = -3.0 * hg_96[k]
                  + f_0 * kg_96[k];

        t_97[k] = -3.0 * hg_97[k]
                  + f_0 * kg_97[k];

        t_98[k] = -3.0 * hg_98[k]
                  + f_0 * kg_98[k];

        t_99[k] = -3.0 * hg_99[k]
                  + f_0 * kg_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, hg_100, hg_101, hg_102, hg_103, \
                         hg_104, kg_100, kg_101, kg_102, kg_103, \
                         kg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * hg_100[k]
                   + f_0 * kg_100[k];

        t_101[k] = -3.0 * hg_101[k]
                   + f_0 * kg_101[k];

        t_102[k] = -3.0 * hg_102[k]
                   + f_0 * kg_102[k];

        t_103[k] = -3.0 * hg_103[k]
                   + f_0 * kg_103[k];

        t_104[k] = -3.0 * hg_104[k]
                   + f_0 * kg_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, hg_105, hg_106, hg_107, hg_108, \
                         hg_109, kg_105, kg_106, kg_107, kg_108, \
                         kg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -3.0 * hg_105[k]
                   + f_0 * kg_105[k];

        t_106[k] = -3.0 * hg_106[k]
                   + f_0 * kg_106[k];

        t_107[k] = -3.0 * hg_107[k]
                   + f_0 * kg_107[k];

        t_108[k] = -3.0 * hg_108[k]
                   + f_0 * kg_108[k];

        t_109[k] = -3.0 * hg_109[k]
                   + f_0 * kg_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, hg_110, hg_111, hg_112, hg_113, \
                         hg_114, kg_110, kg_111, kg_112, kg_113, \
                         kg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -3.0 * hg_110[k]
                   + f_0 * kg_110[k];

        t_111[k] = -3.0 * hg_111[k]
                   + f_0 * kg_111[k];

        t_112[k] = -3.0 * hg_112[k]
                   + f_0 * kg_112[k];

        t_113[k] = -3.0 * hg_113[k]
                   + f_0 * kg_113[k];

        t_114[k] = -3.0 * hg_114[k]
                   + f_0 * kg_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, hg_115, hg_116, hg_117, hg_118, \
                         hg_119, kg_115, kg_116, kg_117, kg_118, \
                         kg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -3.0 * hg_115[k]
                   + f_0 * kg_115[k];

        t_116[k] = -3.0 * hg_116[k]
                   + f_0 * kg_116[k];

        t_117[k] = -3.0 * hg_117[k]
                   + f_0 * kg_117[k];

        t_118[k] = -3.0 * hg_118[k]
                   + f_0 * kg_118[k];

        t_119[k] = -3.0 * hg_119[k]
                   + f_0 * kg_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, hg_120, hg_121, hg_122, hg_123, \
                         hg_124, kg_120, kg_121, kg_122, kg_123, \
                         kg_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -3.0 * hg_120[k]
                   + f_0 * kg_120[k];

        t_121[k] = -3.0 * hg_121[k]
                   + f_0 * kg_121[k];

        t_122[k] = -3.0 * hg_122[k]
                   + f_0 * kg_122[k];

        t_123[k] = -3.0 * hg_123[k]
                   + f_0 * kg_123[k];

        t_124[k] = -3.0 * hg_124[k]
                   + f_0 * kg_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, hg_125, hg_126, hg_127, hg_128, \
                         hg_129, kg_125, kg_126, kg_127, kg_128, \
                         kg_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -3.0 * hg_125[k]
                   + f_0 * kg_125[k];

        t_126[k] = -3.0 * hg_126[k]
                   + f_0 * kg_126[k];

        t_127[k] = -3.0 * hg_127[k]
                   + f_0 * kg_127[k];

        t_128[k] = -3.0 * hg_128[k]
                   + f_0 * kg_128[k];

        t_129[k] = -3.0 * hg_129[k]
                   + f_0 * kg_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, hg_130, hg_131, hg_132, hg_133, \
                         hg_134, kg_130, kg_131, kg_132, kg_133, \
                         kg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -3.0 * hg_130[k]
                   + f_0 * kg_130[k];

        t_131[k] = -3.0 * hg_131[k]
                   + f_0 * kg_131[k];

        t_132[k] = -3.0 * hg_132[k]
                   + f_0 * kg_132[k];

        t_133[k] = -3.0 * hg_133[k]
                   + f_0 * kg_133[k];

        t_134[k] = -3.0 * hg_134[k]
                   + f_0 * kg_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, hg_135, hg_136, hg_137, hg_138, \
                         hg_139, kg_135, kg_136, kg_137, kg_138, \
                         kg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -3.0 * hg_135[k]
                   + f_0 * kg_135[k];

        t_136[k] = -3.0 * hg_136[k]
                   + f_0 * kg_136[k];

        t_137[k] = -3.0 * hg_137[k]
                   + f_0 * kg_137[k];

        t_138[k] = -3.0 * hg_138[k]
                   + f_0 * kg_138[k];

        t_139[k] = -3.0 * hg_139[k]
                   + f_0 * kg_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, hg_140, hg_141, hg_142, hg_143, \
                         hg_144, kg_140, kg_141, kg_142, kg_143, \
                         kg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -3.0 * hg_140[k]
                   + f_0 * kg_140[k];

        t_141[k] = -3.0 * hg_141[k]
                   + f_0 * kg_141[k];

        t_142[k] = -3.0 * hg_142[k]
                   + f_0 * kg_142[k];

        t_143[k] = -3.0 * hg_143[k]
                   + f_0 * kg_143[k];

        t_144[k] = -3.0 * hg_144[k]
                   + f_0 * kg_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, hg_145, hg_146, hg_147, hg_148, \
                         hg_149, kg_145, kg_146, kg_147, kg_148, \
                         kg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -3.0 * hg_145[k]
                   + f_0 * kg_145[k];

        t_146[k] = -3.0 * hg_146[k]
                   + f_0 * kg_146[k];

        t_147[k] = -3.0 * hg_147[k]
                   + f_0 * kg_147[k];

        t_148[k] = -3.0 * hg_148[k]
                   + f_0 * kg_148[k];

        t_149[k] = -3.0 * hg_149[k]
                   + f_0 * kg_149[k];
    }
}

static auto
compute_prim_geom_10_ig_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hg, const size_t kg,
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

    const auto *kg_150 = buffer.data(kg + 150);
    const auto *kg_151 = buffer.data(kg + 151);
    const auto *kg_152 = buffer.data(kg + 152);
    const auto *kg_153 = buffer.data(kg + 153);
    const auto *kg_154 = buffer.data(kg + 154);
    const auto *kg_155 = buffer.data(kg + 155);
    const auto *kg_156 = buffer.data(kg + 156);
    const auto *kg_157 = buffer.data(kg + 157);
    const auto *kg_158 = buffer.data(kg + 158);
    const auto *kg_159 = buffer.data(kg + 159);
    const auto *kg_160 = buffer.data(kg + 160);
    const auto *kg_161 = buffer.data(kg + 161);
    const auto *kg_162 = buffer.data(kg + 162);
    const auto *kg_163 = buffer.data(kg + 163);
    const auto *kg_164 = buffer.data(kg + 164);
    const auto *kg_165 = buffer.data(kg + 165);
    const auto *kg_166 = buffer.data(kg + 166);
    const auto *kg_167 = buffer.data(kg + 167);
    const auto *kg_168 = buffer.data(kg + 168);
    const auto *kg_169 = buffer.data(kg + 169);
    const auto *kg_170 = buffer.data(kg + 170);
    const auto *kg_171 = buffer.data(kg + 171);
    const auto *kg_172 = buffer.data(kg + 172);
    const auto *kg_173 = buffer.data(kg + 173);
    const auto *kg_174 = buffer.data(kg + 174);
    const auto *kg_175 = buffer.data(kg + 175);
    const auto *kg_176 = buffer.data(kg + 176);
    const auto *kg_177 = buffer.data(kg + 177);
    const auto *kg_178 = buffer.data(kg + 178);
    const auto *kg_179 = buffer.data(kg + 179);
    const auto *kg_180 = buffer.data(kg + 180);
    const auto *kg_181 = buffer.data(kg + 181);
    const auto *kg_182 = buffer.data(kg + 182);
    const auto *kg_183 = buffer.data(kg + 183);
    const auto *kg_184 = buffer.data(kg + 184);
    const auto *kg_185 = buffer.data(kg + 185);
    const auto *kg_186 = buffer.data(kg + 186);
    const auto *kg_187 = buffer.data(kg + 187);
    const auto *kg_188 = buffer.data(kg + 188);
    const auto *kg_189 = buffer.data(kg + 189);
    const auto *kg_190 = buffer.data(kg + 190);
    const auto *kg_191 = buffer.data(kg + 191);
    const auto *kg_192 = buffer.data(kg + 192);
    const auto *kg_193 = buffer.data(kg + 193);
    const auto *kg_194 = buffer.data(kg + 194);
    const auto *kg_195 = buffer.data(kg + 195);
    const auto *kg_196 = buffer.data(kg + 196);
    const auto *kg_197 = buffer.data(kg + 197);
    const auto *kg_198 = buffer.data(kg + 198);
    const auto *kg_199 = buffer.data(kg + 199);
    const auto *kg_200 = buffer.data(kg + 200);
    const auto *kg_201 = buffer.data(kg + 201);
    const auto *kg_202 = buffer.data(kg + 202);
    const auto *kg_203 = buffer.data(kg + 203);
    const auto *kg_204 = buffer.data(kg + 204);
    const auto *kg_205 = buffer.data(kg + 205);
    const auto *kg_206 = buffer.data(kg + 206);
    const auto *kg_207 = buffer.data(kg + 207);
    const auto *kg_208 = buffer.data(kg + 208);
    const auto *kg_209 = buffer.data(kg + 209);
    const auto *kg_210 = buffer.data(kg + 210);
    const auto *kg_211 = buffer.data(kg + 211);
    const auto *kg_212 = buffer.data(kg + 212);
    const auto *kg_213 = buffer.data(kg + 213);
    const auto *kg_214 = buffer.data(kg + 214);
    const auto *kg_215 = buffer.data(kg + 215);
    const auto *kg_216 = buffer.data(kg + 216);
    const auto *kg_217 = buffer.data(kg + 217);
    const auto *kg_218 = buffer.data(kg + 218);
    const auto *kg_219 = buffer.data(kg + 219);
    const auto *kg_220 = buffer.data(kg + 220);
    const auto *kg_221 = buffer.data(kg + 221);
    const auto *kg_222 = buffer.data(kg + 222);
    const auto *kg_223 = buffer.data(kg + 223);
    const auto *kg_224 = buffer.data(kg + 224);
    const auto *kg_225 = buffer.data(kg + 225);
    const auto *kg_226 = buffer.data(kg + 226);
    const auto *kg_227 = buffer.data(kg + 227);
    const auto *kg_228 = buffer.data(kg + 228);
    const auto *kg_229 = buffer.data(kg + 229);
    const auto *kg_230 = buffer.data(kg + 230);
    const auto *kg_231 = buffer.data(kg + 231);
    const auto *kg_232 = buffer.data(kg + 232);
    const auto *kg_233 = buffer.data(kg + 233);
    const auto *kg_234 = buffer.data(kg + 234);
    const auto *kg_235 = buffer.data(kg + 235);
    const auto *kg_236 = buffer.data(kg + 236);
    const auto *kg_237 = buffer.data(kg + 237);
    const auto *kg_238 = buffer.data(kg + 238);
    const auto *kg_239 = buffer.data(kg + 239);
    const auto *kg_240 = buffer.data(kg + 240);
    const auto *kg_241 = buffer.data(kg + 241);
    const auto *kg_242 = buffer.data(kg + 242);
    const auto *kg_243 = buffer.data(kg + 243);
    const auto *kg_244 = buffer.data(kg + 244);
    const auto *kg_245 = buffer.data(kg + 245);
    const auto *kg_246 = buffer.data(kg + 246);
    const auto *kg_247 = buffer.data(kg + 247);
    const auto *kg_248 = buffer.data(kg + 248);
    const auto *kg_249 = buffer.data(kg + 249);
    const auto *kg_250 = buffer.data(kg + 250);
    const auto *kg_251 = buffer.data(kg + 251);
    const auto *kg_252 = buffer.data(kg + 252);
    const auto *kg_253 = buffer.data(kg + 253);
    const auto *kg_254 = buffer.data(kg + 254);
    const auto *kg_255 = buffer.data(kg + 255);
    const auto *kg_256 = buffer.data(kg + 256);
    const auto *kg_257 = buffer.data(kg + 257);
    const auto *kg_258 = buffer.data(kg + 258);
    const auto *kg_259 = buffer.data(kg + 259);
    const auto *kg_260 = buffer.data(kg + 260);
    const auto *kg_261 = buffer.data(kg + 261);
    const auto *kg_262 = buffer.data(kg + 262);
    const auto *kg_263 = buffer.data(kg + 263);
    const auto *kg_264 = buffer.data(kg + 264);
    const auto *kg_265 = buffer.data(kg + 265);
    const auto *kg_266 = buffer.data(kg + 266);
    const auto *kg_267 = buffer.data(kg + 267);
    const auto *kg_268 = buffer.data(kg + 268);
    const auto *kg_269 = buffer.data(kg + 269);
    const auto *kg_270 = buffer.data(kg + 270);
    const auto *kg_271 = buffer.data(kg + 271);
    const auto *kg_272 = buffer.data(kg + 272);
    const auto *kg_273 = buffer.data(kg + 273);
    const auto *kg_274 = buffer.data(kg + 274);
    const auto *kg_275 = buffer.data(kg + 275);
    const auto *kg_276 = buffer.data(kg + 276);
    const auto *kg_277 = buffer.data(kg + 277);
    const auto *kg_278 = buffer.data(kg + 278);
    const auto *kg_279 = buffer.data(kg + 279);
    const auto *kg_280 = buffer.data(kg + 280);
    const auto *kg_281 = buffer.data(kg + 281);
    const auto *kg_282 = buffer.data(kg + 282);
    const auto *kg_283 = buffer.data(kg + 283);
    const auto *kg_284 = buffer.data(kg + 284);
    const auto *kg_285 = buffer.data(kg + 285);
    const auto *kg_286 = buffer.data(kg + 286);
    const auto *kg_287 = buffer.data(kg + 287);
    const auto *kg_288 = buffer.data(kg + 288);
    const auto *kg_289 = buffer.data(kg + 289);
    const auto *kg_290 = buffer.data(kg + 290);
    const auto *kg_291 = buffer.data(kg + 291);
    const auto *kg_292 = buffer.data(kg + 292);
    const auto *kg_293 = buffer.data(kg + 293);
    const auto *kg_294 = buffer.data(kg + 294);
    const auto *kg_295 = buffer.data(kg + 295);
    const auto *kg_296 = buffer.data(kg + 296);
    const auto *kg_297 = buffer.data(kg + 297);
    const auto *kg_298 = buffer.data(kg + 298);
    const auto *kg_299 = buffer.data(kg + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, hg_150, hg_151, hg_152, hg_153, \
                         hg_154, kg_150, kg_151, kg_152, kg_153, \
                         kg_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -2.0 * hg_150[k]
                   + f_0 * kg_150[k];

        t_151[k] = -2.0 * hg_151[k]
                   + f_0 * kg_151[k];

        t_152[k] = -2.0 * hg_152[k]
                   + f_0 * kg_152[k];

        t_153[k] = -2.0 * hg_153[k]
                   + f_0 * kg_153[k];

        t_154[k] = -2.0 * hg_154[k]
                   + f_0 * kg_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, hg_155, hg_156, hg_157, hg_158, \
                         hg_159, kg_155, kg_156, kg_157, kg_158, \
                         kg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -2.0 * hg_155[k]
                   + f_0 * kg_155[k];

        t_156[k] = -2.0 * hg_156[k]
                   + f_0 * kg_156[k];

        t_157[k] = -2.0 * hg_157[k]
                   + f_0 * kg_157[k];

        t_158[k] = -2.0 * hg_158[k]
                   + f_0 * kg_158[k];

        t_159[k] = -2.0 * hg_159[k]
                   + f_0 * kg_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, hg_160, hg_161, hg_162, hg_163, \
                         hg_164, kg_160, kg_161, kg_162, kg_163, \
                         kg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -2.0 * hg_160[k]
                   + f_0 * kg_160[k];

        t_161[k] = -2.0 * hg_161[k]
                   + f_0 * kg_161[k];

        t_162[k] = -2.0 * hg_162[k]
                   + f_0 * kg_162[k];

        t_163[k] = -2.0 * hg_163[k]
                   + f_0 * kg_163[k];

        t_164[k] = -2.0 * hg_164[k]
                   + f_0 * kg_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, hg_165, hg_166, hg_167, hg_168, \
                         hg_169, kg_165, kg_166, kg_167, kg_168, \
                         kg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -2.0 * hg_165[k]
                   + f_0 * kg_165[k];

        t_166[k] = -2.0 * hg_166[k]
                   + f_0 * kg_166[k];

        t_167[k] = -2.0 * hg_167[k]
                   + f_0 * kg_167[k];

        t_168[k] = -2.0 * hg_168[k]
                   + f_0 * kg_168[k];

        t_169[k] = -2.0 * hg_169[k]
                   + f_0 * kg_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, hg_170, hg_171, hg_172, hg_173, \
                         hg_174, kg_170, kg_171, kg_172, kg_173, \
                         kg_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -2.0 * hg_170[k]
                   + f_0 * kg_170[k];

        t_171[k] = -2.0 * hg_171[k]
                   + f_0 * kg_171[k];

        t_172[k] = -2.0 * hg_172[k]
                   + f_0 * kg_172[k];

        t_173[k] = -2.0 * hg_173[k]
                   + f_0 * kg_173[k];

        t_174[k] = -2.0 * hg_174[k]
                   + f_0 * kg_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, hg_175, hg_176, hg_177, hg_178, \
                         hg_179, kg_175, kg_176, kg_177, kg_178, \
                         kg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -2.0 * hg_175[k]
                   + f_0 * kg_175[k];

        t_176[k] = -2.0 * hg_176[k]
                   + f_0 * kg_176[k];

        t_177[k] = -2.0 * hg_177[k]
                   + f_0 * kg_177[k];

        t_178[k] = -2.0 * hg_178[k]
                   + f_0 * kg_178[k];

        t_179[k] = -2.0 * hg_179[k]
                   + f_0 * kg_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, hg_180, hg_181, hg_182, hg_183, \
                         hg_184, kg_180, kg_181, kg_182, kg_183, \
                         kg_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * hg_180[k]
                   + f_0 * kg_180[k];

        t_181[k] = -2.0 * hg_181[k]
                   + f_0 * kg_181[k];

        t_182[k] = -2.0 * hg_182[k]
                   + f_0 * kg_182[k];

        t_183[k] = -2.0 * hg_183[k]
                   + f_0 * kg_183[k];

        t_184[k] = -2.0 * hg_184[k]
                   + f_0 * kg_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, hg_185, hg_186, hg_187, hg_188, \
                         hg_189, kg_185, kg_186, kg_187, kg_188, \
                         kg_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * hg_185[k]
                   + f_0 * kg_185[k];

        t_186[k] = -2.0 * hg_186[k]
                   + f_0 * kg_186[k];

        t_187[k] = -2.0 * hg_187[k]
                   + f_0 * kg_187[k];

        t_188[k] = -2.0 * hg_188[k]
                   + f_0 * kg_188[k];

        t_189[k] = -2.0 * hg_189[k]
                   + f_0 * kg_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, hg_190, hg_191, hg_192, hg_193, \
                         hg_194, kg_190, kg_191, kg_192, kg_193, \
                         kg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * hg_190[k]
                   + f_0 * kg_190[k];

        t_191[k] = -2.0 * hg_191[k]
                   + f_0 * kg_191[k];

        t_192[k] = -2.0 * hg_192[k]
                   + f_0 * kg_192[k];

        t_193[k] = -2.0 * hg_193[k]
                   + f_0 * kg_193[k];

        t_194[k] = -2.0 * hg_194[k]
                   + f_0 * kg_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, hg_195, hg_196, hg_197, hg_198, \
                         hg_199, kg_195, kg_196, kg_197, kg_198, \
                         kg_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * hg_195[k]
                   + f_0 * kg_195[k];

        t_196[k] = -2.0 * hg_196[k]
                   + f_0 * kg_196[k];

        t_197[k] = -2.0 * hg_197[k]
                   + f_0 * kg_197[k];

        t_198[k] = -2.0 * hg_198[k]
                   + f_0 * kg_198[k];

        t_199[k] = -2.0 * hg_199[k]
                   + f_0 * kg_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, hg_200, hg_201, hg_202, hg_203, \
                         hg_204, kg_200, kg_201, kg_202, kg_203, \
                         kg_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * hg_200[k]
                   + f_0 * kg_200[k];

        t_201[k] = -2.0 * hg_201[k]
                   + f_0 * kg_201[k];

        t_202[k] = -2.0 * hg_202[k]
                   + f_0 * kg_202[k];

        t_203[k] = -2.0 * hg_203[k]
                   + f_0 * kg_203[k];

        t_204[k] = -2.0 * hg_204[k]
                   + f_0 * kg_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, hg_205, hg_206, hg_207, hg_208, \
                         hg_209, kg_205, kg_206, kg_207, kg_208, \
                         kg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * hg_205[k]
                   + f_0 * kg_205[k];

        t_206[k] = -2.0 * hg_206[k]
                   + f_0 * kg_206[k];

        t_207[k] = -2.0 * hg_207[k]
                   + f_0 * kg_207[k];

        t_208[k] = -2.0 * hg_208[k]
                   + f_0 * kg_208[k];

        t_209[k] = -2.0 * hg_209[k]
                   + f_0 * kg_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, hg_210, hg_211, hg_212, hg_213, \
                         hg_214, kg_210, kg_211, kg_212, kg_213, \
                         kg_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * hg_210[k]
                   + f_0 * kg_210[k];

        t_211[k] = -2.0 * hg_211[k]
                   + f_0 * kg_211[k];

        t_212[k] = -2.0 * hg_212[k]
                   + f_0 * kg_212[k];

        t_213[k] = -2.0 * hg_213[k]
                   + f_0 * kg_213[k];

        t_214[k] = -2.0 * hg_214[k]
                   + f_0 * kg_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, hg_215, hg_216, hg_217, hg_218, \
                         hg_219, kg_215, kg_216, kg_217, kg_218, \
                         kg_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * hg_215[k]
                   + f_0 * kg_215[k];

        t_216[k] = -2.0 * hg_216[k]
                   + f_0 * kg_216[k];

        t_217[k] = -2.0 * hg_217[k]
                   + f_0 * kg_217[k];

        t_218[k] = -2.0 * hg_218[k]
                   + f_0 * kg_218[k];

        t_219[k] = -2.0 * hg_219[k]
                   + f_0 * kg_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, hg_220, hg_221, hg_222, hg_223, \
                         hg_224, kg_220, kg_221, kg_222, kg_223, \
                         kg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * hg_220[k]
                   + f_0 * kg_220[k];

        t_221[k] = -2.0 * hg_221[k]
                   + f_0 * kg_221[k];

        t_222[k] = -2.0 * hg_222[k]
                   + f_0 * kg_222[k];

        t_223[k] = -2.0 * hg_223[k]
                   + f_0 * kg_223[k];

        t_224[k] = -2.0 * hg_224[k]
                   + f_0 * kg_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, hg_225, hg_226, hg_227, hg_228, \
                         hg_229, kg_225, kg_226, kg_227, kg_228, \
                         kg_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -hg_225[k]
                   + f_0 * kg_225[k];

        t_226[k] = -hg_226[k]
                   + f_0 * kg_226[k];

        t_227[k] = -hg_227[k]
                   + f_0 * kg_227[k];

        t_228[k] = -hg_228[k]
                   + f_0 * kg_228[k];

        t_229[k] = -hg_229[k]
                   + f_0 * kg_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, hg_230, hg_231, hg_232, hg_233, \
                         hg_234, kg_230, kg_231, kg_232, kg_233, \
                         kg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -hg_230[k]
                   + f_0 * kg_230[k];

        t_231[k] = -hg_231[k]
                   + f_0 * kg_231[k];

        t_232[k] = -hg_232[k]
                   + f_0 * kg_232[k];

        t_233[k] = -hg_233[k]
                   + f_0 * kg_233[k];

        t_234[k] = -hg_234[k]
                   + f_0 * kg_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, hg_235, hg_236, hg_237, hg_238, \
                         hg_239, kg_235, kg_236, kg_237, kg_238, \
                         kg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -hg_235[k]
                   + f_0 * kg_235[k];

        t_236[k] = -hg_236[k]
                   + f_0 * kg_236[k];

        t_237[k] = -hg_237[k]
                   + f_0 * kg_237[k];

        t_238[k] = -hg_238[k]
                   + f_0 * kg_238[k];

        t_239[k] = -hg_239[k]
                   + f_0 * kg_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, hg_240, hg_241, hg_242, hg_243, \
                         hg_244, kg_240, kg_241, kg_242, kg_243, \
                         kg_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -hg_240[k]
                   + f_0 * kg_240[k];

        t_241[k] = -hg_241[k]
                   + f_0 * kg_241[k];

        t_242[k] = -hg_242[k]
                   + f_0 * kg_242[k];

        t_243[k] = -hg_243[k]
                   + f_0 * kg_243[k];

        t_244[k] = -hg_244[k]
                   + f_0 * kg_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, hg_245, hg_246, hg_247, hg_248, \
                         hg_249, kg_245, kg_246, kg_247, kg_248, \
                         kg_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -hg_245[k]
                   + f_0 * kg_245[k];

        t_246[k] = -hg_246[k]
                   + f_0 * kg_246[k];

        t_247[k] = -hg_247[k]
                   + f_0 * kg_247[k];

        t_248[k] = -hg_248[k]
                   + f_0 * kg_248[k];

        t_249[k] = -hg_249[k]
                   + f_0 * kg_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, hg_250, hg_251, hg_252, hg_253, \
                         hg_254, kg_250, kg_251, kg_252, kg_253, \
                         kg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -hg_250[k]
                   + f_0 * kg_250[k];

        t_251[k] = -hg_251[k]
                   + f_0 * kg_251[k];

        t_252[k] = -hg_252[k]
                   + f_0 * kg_252[k];

        t_253[k] = -hg_253[k]
                   + f_0 * kg_253[k];

        t_254[k] = -hg_254[k]
                   + f_0 * kg_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, hg_255, hg_256, hg_257, hg_258, \
                         hg_259, kg_255, kg_256, kg_257, kg_258, \
                         kg_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -hg_255[k]
                   + f_0 * kg_255[k];

        t_256[k] = -hg_256[k]
                   + f_0 * kg_256[k];

        t_257[k] = -hg_257[k]
                   + f_0 * kg_257[k];

        t_258[k] = -hg_258[k]
                   + f_0 * kg_258[k];

        t_259[k] = -hg_259[k]
                   + f_0 * kg_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, hg_260, hg_261, hg_262, hg_263, \
                         hg_264, kg_260, kg_261, kg_262, kg_263, \
                         kg_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -hg_260[k]
                   + f_0 * kg_260[k];

        t_261[k] = -hg_261[k]
                   + f_0 * kg_261[k];

        t_262[k] = -hg_262[k]
                   + f_0 * kg_262[k];

        t_263[k] = -hg_263[k]
                   + f_0 * kg_263[k];

        t_264[k] = -hg_264[k]
                   + f_0 * kg_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, hg_265, hg_266, hg_267, hg_268, \
                         hg_269, kg_265, kg_266, kg_267, kg_268, \
                         kg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -hg_265[k]
                   + f_0 * kg_265[k];

        t_266[k] = -hg_266[k]
                   + f_0 * kg_266[k];

        t_267[k] = -hg_267[k]
                   + f_0 * kg_267[k];

        t_268[k] = -hg_268[k]
                   + f_0 * kg_268[k];

        t_269[k] = -hg_269[k]
                   + f_0 * kg_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, hg_270, hg_271, hg_272, hg_273, \
                         hg_274, kg_270, kg_271, kg_272, kg_273, \
                         kg_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -hg_270[k]
                   + f_0 * kg_270[k];

        t_271[k] = -hg_271[k]
                   + f_0 * kg_271[k];

        t_272[k] = -hg_272[k]
                   + f_0 * kg_272[k];

        t_273[k] = -hg_273[k]
                   + f_0 * kg_273[k];

        t_274[k] = -hg_274[k]
                   + f_0 * kg_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, hg_275, hg_276, hg_277, hg_278, \
                         hg_279, kg_275, kg_276, kg_277, kg_278, \
                         kg_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -hg_275[k]
                   + f_0 * kg_275[k];

        t_276[k] = -hg_276[k]
                   + f_0 * kg_276[k];

        t_277[k] = -hg_277[k]
                   + f_0 * kg_277[k];

        t_278[k] = -hg_278[k]
                   + f_0 * kg_278[k];

        t_279[k] = -hg_279[k]
                   + f_0 * kg_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, hg_280, hg_281, hg_282, hg_283, \
                         hg_284, kg_280, kg_281, kg_282, kg_283, \
                         kg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -hg_280[k]
                   + f_0 * kg_280[k];

        t_281[k] = -hg_281[k]
                   + f_0 * kg_281[k];

        t_282[k] = -hg_282[k]
                   + f_0 * kg_282[k];

        t_283[k] = -hg_283[k]
                   + f_0 * kg_283[k];

        t_284[k] = -hg_284[k]
                   + f_0 * kg_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, hg_285, hg_286, hg_287, hg_288, \
                         hg_289, kg_285, kg_286, kg_287, kg_288, \
                         kg_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -hg_285[k]
                   + f_0 * kg_285[k];

        t_286[k] = -hg_286[k]
                   + f_0 * kg_286[k];

        t_287[k] = -hg_287[k]
                   + f_0 * kg_287[k];

        t_288[k] = -hg_288[k]
                   + f_0 * kg_288[k];

        t_289[k] = -hg_289[k]
                   + f_0 * kg_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, hg_290, hg_291, hg_292, hg_293, \
                         hg_294, kg_290, kg_291, kg_292, kg_293, \
                         kg_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -hg_290[k]
                   + f_0 * kg_290[k];

        t_291[k] = -hg_291[k]
                   + f_0 * kg_291[k];

        t_292[k] = -hg_292[k]
                   + f_0 * kg_292[k];

        t_293[k] = -hg_293[k]
                   + f_0 * kg_293[k];

        t_294[k] = -hg_294[k]
                   + f_0 * kg_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, hg_295, hg_296, hg_297, hg_298, \
                         hg_299, kg_295, kg_296, kg_297, kg_298, \
                         kg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -hg_295[k]
                   + f_0 * kg_295[k];

        t_296[k] = -hg_296[k]
                   + f_0 * kg_296[k];

        t_297[k] = -hg_297[k]
                   + f_0 * kg_297[k];

        t_298[k] = -hg_298[k]
                   + f_0 * kg_298[k];

        t_299[k] = -hg_299[k]
                   + f_0 * kg_299[k];
    }
}

static auto
compute_prim_geom_10_ig_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hg, const size_t kg,
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

    const auto *kg_300 = buffer.data(kg + 300);
    const auto *kg_301 = buffer.data(kg + 301);
    const auto *kg_302 = buffer.data(kg + 302);
    const auto *kg_303 = buffer.data(kg + 303);
    const auto *kg_304 = buffer.data(kg + 304);
    const auto *kg_305 = buffer.data(kg + 305);
    const auto *kg_306 = buffer.data(kg + 306);
    const auto *kg_307 = buffer.data(kg + 307);
    const auto *kg_308 = buffer.data(kg + 308);
    const auto *kg_309 = buffer.data(kg + 309);
    const auto *kg_310 = buffer.data(kg + 310);
    const auto *kg_311 = buffer.data(kg + 311);
    const auto *kg_312 = buffer.data(kg + 312);
    const auto *kg_313 = buffer.data(kg + 313);
    const auto *kg_314 = buffer.data(kg + 314);
    const auto *kg_315 = buffer.data(kg + 315);
    const auto *kg_316 = buffer.data(kg + 316);
    const auto *kg_317 = buffer.data(kg + 317);
    const auto *kg_318 = buffer.data(kg + 318);
    const auto *kg_319 = buffer.data(kg + 319);
    const auto *kg_320 = buffer.data(kg + 320);
    const auto *kg_321 = buffer.data(kg + 321);
    const auto *kg_322 = buffer.data(kg + 322);
    const auto *kg_323 = buffer.data(kg + 323);
    const auto *kg_324 = buffer.data(kg + 324);
    const auto *kg_325 = buffer.data(kg + 325);
    const auto *kg_326 = buffer.data(kg + 326);
    const auto *kg_327 = buffer.data(kg + 327);
    const auto *kg_328 = buffer.data(kg + 328);
    const auto *kg_329 = buffer.data(kg + 329);
    const auto *kg_330 = buffer.data(kg + 330);
    const auto *kg_331 = buffer.data(kg + 331);
    const auto *kg_332 = buffer.data(kg + 332);
    const auto *kg_333 = buffer.data(kg + 333);
    const auto *kg_334 = buffer.data(kg + 334);
    const auto *kg_335 = buffer.data(kg + 335);
    const auto *kg_336 = buffer.data(kg + 336);
    const auto *kg_337 = buffer.data(kg + 337);
    const auto *kg_338 = buffer.data(kg + 338);
    const auto *kg_339 = buffer.data(kg + 339);
    const auto *kg_340 = buffer.data(kg + 340);
    const auto *kg_341 = buffer.data(kg + 341);
    const auto *kg_342 = buffer.data(kg + 342);
    const auto *kg_343 = buffer.data(kg + 343);
    const auto *kg_344 = buffer.data(kg + 344);
    const auto *kg_345 = buffer.data(kg + 345);
    const auto *kg_346 = buffer.data(kg + 346);
    const auto *kg_347 = buffer.data(kg + 347);
    const auto *kg_348 = buffer.data(kg + 348);
    const auto *kg_349 = buffer.data(kg + 349);
    const auto *kg_350 = buffer.data(kg + 350);
    const auto *kg_351 = buffer.data(kg + 351);
    const auto *kg_352 = buffer.data(kg + 352);
    const auto *kg_353 = buffer.data(kg + 353);
    const auto *kg_354 = buffer.data(kg + 354);
    const auto *kg_355 = buffer.data(kg + 355);
    const auto *kg_356 = buffer.data(kg + 356);
    const auto *kg_357 = buffer.data(kg + 357);
    const auto *kg_358 = buffer.data(kg + 358);
    const auto *kg_359 = buffer.data(kg + 359);
    const auto *kg_360 = buffer.data(kg + 360);
    const auto *kg_361 = buffer.data(kg + 361);
    const auto *kg_362 = buffer.data(kg + 362);
    const auto *kg_363 = buffer.data(kg + 363);
    const auto *kg_364 = buffer.data(kg + 364);
    const auto *kg_365 = buffer.data(kg + 365);
    const auto *kg_366 = buffer.data(kg + 366);
    const auto *kg_367 = buffer.data(kg + 367);
    const auto *kg_368 = buffer.data(kg + 368);
    const auto *kg_369 = buffer.data(kg + 369);
    const auto *kg_370 = buffer.data(kg + 370);
    const auto *kg_371 = buffer.data(kg + 371);
    const auto *kg_372 = buffer.data(kg + 372);
    const auto *kg_373 = buffer.data(kg + 373);
    const auto *kg_374 = buffer.data(kg + 374);
    const auto *kg_375 = buffer.data(kg + 375);
    const auto *kg_376 = buffer.data(kg + 376);
    const auto *kg_377 = buffer.data(kg + 377);
    const auto *kg_378 = buffer.data(kg + 378);
    const auto *kg_379 = buffer.data(kg + 379);
    const auto *kg_380 = buffer.data(kg + 380);
    const auto *kg_381 = buffer.data(kg + 381);
    const auto *kg_382 = buffer.data(kg + 382);
    const auto *kg_383 = buffer.data(kg + 383);
    const auto *kg_384 = buffer.data(kg + 384);
    const auto *kg_385 = buffer.data(kg + 385);
    const auto *kg_386 = buffer.data(kg + 386);
    const auto *kg_387 = buffer.data(kg + 387);
    const auto *kg_388 = buffer.data(kg + 388);
    const auto *kg_389 = buffer.data(kg + 389);
    const auto *kg_390 = buffer.data(kg + 390);
    const auto *kg_391 = buffer.data(kg + 391);
    const auto *kg_392 = buffer.data(kg + 392);
    const auto *kg_393 = buffer.data(kg + 393);
    const auto *kg_394 = buffer.data(kg + 394);
    const auto *kg_395 = buffer.data(kg + 395);
    const auto *kg_396 = buffer.data(kg + 396);
    const auto *kg_397 = buffer.data(kg + 397);
    const auto *kg_398 = buffer.data(kg + 398);
    const auto *kg_399 = buffer.data(kg + 399);
    const auto *kg_400 = buffer.data(kg + 400);
    const auto *kg_401 = buffer.data(kg + 401);
    const auto *kg_402 = buffer.data(kg + 402);
    const auto *kg_403 = buffer.data(kg + 403);
    const auto *kg_404 = buffer.data(kg + 404);
    const auto *kg_405 = buffer.data(kg + 405);
    const auto *kg_406 = buffer.data(kg + 406);
    const auto *kg_407 = buffer.data(kg + 407);
    const auto *kg_408 = buffer.data(kg + 408);
    const auto *kg_409 = buffer.data(kg + 409);
    const auto *kg_410 = buffer.data(kg + 410);
    const auto *kg_411 = buffer.data(kg + 411);
    const auto *kg_412 = buffer.data(kg + 412);
    const auto *kg_413 = buffer.data(kg + 413);
    const auto *kg_414 = buffer.data(kg + 414);
    const auto *kg_415 = buffer.data(kg + 415);
    const auto *kg_416 = buffer.data(kg + 416);
    const auto *kg_417 = buffer.data(kg + 417);
    const auto *kg_418 = buffer.data(kg + 418);
    const auto *kg_419 = buffer.data(kg + 419);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, hg_300, hg_301, hg_302, hg_303, \
                         hg_304, kg_300, kg_301, kg_302, kg_303, \
                         kg_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -hg_300[k]
                   + f_0 * kg_300[k];

        t_301[k] = -hg_301[k]
                   + f_0 * kg_301[k];

        t_302[k] = -hg_302[k]
                   + f_0 * kg_302[k];

        t_303[k] = -hg_303[k]
                   + f_0 * kg_303[k];

        t_304[k] = -hg_304[k]
                   + f_0 * kg_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, hg_305, hg_306, hg_307, hg_308, \
                         hg_309, kg_305, kg_306, kg_307, kg_308, \
                         kg_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -hg_305[k]
                   + f_0 * kg_305[k];

        t_306[k] = -hg_306[k]
                   + f_0 * kg_306[k];

        t_307[k] = -hg_307[k]
                   + f_0 * kg_307[k];

        t_308[k] = -hg_308[k]
                   + f_0 * kg_308[k];

        t_309[k] = -hg_309[k]
                   + f_0 * kg_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, hg_310, hg_311, hg_312, hg_313, \
                         hg_314, kg_310, kg_311, kg_312, kg_313, \
                         kg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -hg_310[k]
                   + f_0 * kg_310[k];

        t_311[k] = -hg_311[k]
                   + f_0 * kg_311[k];

        t_312[k] = -hg_312[k]
                   + f_0 * kg_312[k];

        t_313[k] = -hg_313[k]
                   + f_0 * kg_313[k];

        t_314[k] = -hg_314[k]
                   + f_0 * kg_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, t_320, t_321, t_322, kg_315, \
                         kg_316, kg_317, kg_318, kg_319, kg_320, kg_321, \
                         kg_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_0 * kg_315[k];

        t_316[k] = f_0 * kg_316[k];

        t_317[k] = f_0 * kg_317[k];

        t_318[k] = f_0 * kg_318[k];

        t_319[k] = f_0 * kg_319[k];

        t_320[k] = f_0 * kg_320[k];

        t_321[k] = f_0 * kg_321[k];

        t_322[k] = f_0 * kg_322[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, t_328, t_329, t_330, kg_323, \
                         kg_324, kg_325, kg_326, kg_327, kg_328, kg_329, \
                         kg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_0 * kg_323[k];

        t_324[k] = f_0 * kg_324[k];

        t_325[k] = f_0 * kg_325[k];

        t_326[k] = f_0 * kg_326[k];

        t_327[k] = f_0 * kg_327[k];

        t_328[k] = f_0 * kg_328[k];

        t_329[k] = f_0 * kg_329[k];

        t_330[k] = f_0 * kg_330[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, t_335, t_336, t_337, t_338, kg_331, \
                         kg_332, kg_333, kg_334, kg_335, kg_336, kg_337, \
                         kg_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_0 * kg_331[k];

        t_332[k] = f_0 * kg_332[k];

        t_333[k] = f_0 * kg_333[k];

        t_334[k] = f_0 * kg_334[k];

        t_335[k] = f_0 * kg_335[k];

        t_336[k] = f_0 * kg_336[k];

        t_337[k] = f_0 * kg_337[k];

        t_338[k] = f_0 * kg_338[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, t_343, t_344, t_345, t_346, kg_339, \
                         kg_340, kg_341, kg_342, kg_343, kg_344, kg_345, \
                         kg_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_0 * kg_339[k];

        t_340[k] = f_0 * kg_340[k];

        t_341[k] = f_0 * kg_341[k];

        t_342[k] = f_0 * kg_342[k];

        t_343[k] = f_0 * kg_343[k];

        t_344[k] = f_0 * kg_344[k];

        t_345[k] = f_0 * kg_345[k];

        t_346[k] = f_0 * kg_346[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, t_352, t_353, t_354, kg_347, \
                         kg_348, kg_349, kg_350, kg_351, kg_352, kg_353, \
                         kg_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_0 * kg_347[k];

        t_348[k] = f_0 * kg_348[k];

        t_349[k] = f_0 * kg_349[k];

        t_350[k] = f_0 * kg_350[k];

        t_351[k] = f_0 * kg_351[k];

        t_352[k] = f_0 * kg_352[k];

        t_353[k] = f_0 * kg_353[k];

        t_354[k] = f_0 * kg_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, t_360, t_361, t_362, kg_355, \
                         kg_356, kg_357, kg_358, kg_359, kg_360, kg_361, \
                         kg_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_0 * kg_355[k];

        t_356[k] = f_0 * kg_356[k];

        t_357[k] = f_0 * kg_357[k];

        t_358[k] = f_0 * kg_358[k];

        t_359[k] = f_0 * kg_359[k];

        t_360[k] = f_0 * kg_360[k];

        t_361[k] = f_0 * kg_361[k];

        t_362[k] = f_0 * kg_362[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, t_368, t_369, t_370, kg_363, \
                         kg_364, kg_365, kg_366, kg_367, kg_368, kg_369, \
                         kg_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_0 * kg_363[k];

        t_364[k] = f_0 * kg_364[k];

        t_365[k] = f_0 * kg_365[k];

        t_366[k] = f_0 * kg_366[k];

        t_367[k] = f_0 * kg_367[k];

        t_368[k] = f_0 * kg_368[k];

        t_369[k] = f_0 * kg_369[k];

        t_370[k] = f_0 * kg_370[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, t_375, t_376, t_377, t_378, kg_371, \
                         kg_372, kg_373, kg_374, kg_375, kg_376, kg_377, \
                         kg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_0 * kg_371[k];

        t_372[k] = f_0 * kg_372[k];

        t_373[k] = f_0 * kg_373[k];

        t_374[k] = f_0 * kg_374[k];

        t_375[k] = f_0 * kg_375[k];

        t_376[k] = f_0 * kg_376[k];

        t_377[k] = f_0 * kg_377[k];

        t_378[k] = f_0 * kg_378[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, t_384, t_385, t_386, kg_379, \
                         kg_380, kg_381, kg_382, kg_383, kg_384, kg_385, \
                         kg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_0 * kg_379[k];

        t_380[k] = f_0 * kg_380[k];

        t_381[k] = f_0 * kg_381[k];

        t_382[k] = f_0 * kg_382[k];

        t_383[k] = f_0 * kg_383[k];

        t_384[k] = f_0 * kg_384[k];

        t_385[k] = f_0 * kg_385[k];

        t_386[k] = f_0 * kg_386[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, t_392, t_393, t_394, kg_387, \
                         kg_388, kg_389, kg_390, kg_391, kg_392, kg_393, \
                         kg_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_0 * kg_387[k];

        t_388[k] = f_0 * kg_388[k];

        t_389[k] = f_0 * kg_389[k];

        t_390[k] = f_0 * kg_390[k];

        t_391[k] = f_0 * kg_391[k];

        t_392[k] = f_0 * kg_392[k];

        t_393[k] = f_0 * kg_393[k];

        t_394[k] = f_0 * kg_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, t_400, t_401, t_402, kg_395, \
                         kg_396, kg_397, kg_398, kg_399, kg_400, kg_401, \
                         kg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_0 * kg_395[k];

        t_396[k] = f_0 * kg_396[k];

        t_397[k] = f_0 * kg_397[k];

        t_398[k] = f_0 * kg_398[k];

        t_399[k] = f_0 * kg_399[k];

        t_400[k] = f_0 * kg_400[k];

        t_401[k] = f_0 * kg_401[k];

        t_402[k] = f_0 * kg_402[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, t_407, t_408, t_409, t_410, kg_403, \
                         kg_404, kg_405, kg_406, kg_407, kg_408, kg_409, \
                         kg_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_0 * kg_403[k];

        t_404[k] = f_0 * kg_404[k];

        t_405[k] = f_0 * kg_405[k];

        t_406[k] = f_0 * kg_406[k];

        t_407[k] = f_0 * kg_407[k];

        t_408[k] = f_0 * kg_408[k];

        t_409[k] = f_0 * kg_409[k];

        t_410[k] = f_0 * kg_410[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, t_415, t_416, t_417, t_418, kg_411, \
                         kg_412, kg_413, kg_414, kg_415, kg_416, kg_417, \
                         kg_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_0 * kg_411[k];

        t_412[k] = f_0 * kg_412[k];

        t_413[k] = f_0 * kg_413[k];

        t_414[k] = f_0 * kg_414[k];

        t_415[k] = f_0 * kg_415[k];

        t_416[k] = f_0 * kg_416[k];

        t_417[k] = f_0 * kg_417[k];

        t_418[k] = f_0 * kg_418[k];
    }

#pragma omp simd aligned(t_419, kg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_0 * kg_419[k];
    }
}

auto
compute_prim_geom_10_ig_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t hg, const size_t kg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ig_electron_repulsion_0_piece0(buffer, target, hg, kg, ncols, alpha);

    compute_prim_geom_10_ig_electron_repulsion_0_piece1(buffer, target, hg, kg, ncols, alpha);

    compute_prim_geom_10_ig_electron_repulsion_0_piece2(buffer, target, hg, kg, ncols, alpha);
}

static auto
compute_prim_geom_10_ig_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hg, const size_t kg,
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

    const auto *kg_15 = buffer.data(kg + 15);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_17 = buffer.data(kg + 17);
    const auto *kg_18 = buffer.data(kg + 18);
    const auto *kg_19 = buffer.data(kg + 19);
    const auto *kg_20 = buffer.data(kg + 20);
    const auto *kg_21 = buffer.data(kg + 21);
    const auto *kg_22 = buffer.data(kg + 22);
    const auto *kg_23 = buffer.data(kg + 23);
    const auto *kg_24 = buffer.data(kg + 24);
    const auto *kg_25 = buffer.data(kg + 25);
    const auto *kg_26 = buffer.data(kg + 26);
    const auto *kg_27 = buffer.data(kg + 27);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_45 = buffer.data(kg + 45);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_48 = buffer.data(kg + 48);
    const auto *kg_49 = buffer.data(kg + 49);
    const auto *kg_50 = buffer.data(kg + 50);
    const auto *kg_51 = buffer.data(kg + 51);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_53 = buffer.data(kg + 53);
    const auto *kg_54 = buffer.data(kg + 54);
    const auto *kg_55 = buffer.data(kg + 55);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_58 = buffer.data(kg + 58);
    const auto *kg_59 = buffer.data(kg + 59);
    const auto *kg_60 = buffer.data(kg + 60);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_62 = buffer.data(kg + 62);
    const auto *kg_63 = buffer.data(kg + 63);
    const auto *kg_64 = buffer.data(kg + 64);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_67 = buffer.data(kg + 67);
    const auto *kg_68 = buffer.data(kg + 68);
    const auto *kg_69 = buffer.data(kg + 69);
    const auto *kg_70 = buffer.data(kg + 70);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_72 = buffer.data(kg + 72);
    const auto *kg_73 = buffer.data(kg + 73);
    const auto *kg_74 = buffer.data(kg + 74);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_91 = buffer.data(kg + 91);
    const auto *kg_92 = buffer.data(kg + 92);
    const auto *kg_93 = buffer.data(kg + 93);
    const auto *kg_94 = buffer.data(kg + 94);
    const auto *kg_95 = buffer.data(kg + 95);
    const auto *kg_96 = buffer.data(kg + 96);
    const auto *kg_97 = buffer.data(kg + 97);
    const auto *kg_98 = buffer.data(kg + 98);
    const auto *kg_99 = buffer.data(kg + 99);
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_101 = buffer.data(kg + 101);
    const auto *kg_102 = buffer.data(kg + 102);
    const auto *kg_103 = buffer.data(kg + 103);
    const auto *kg_104 = buffer.data(kg + 104);
    const auto *kg_105 = buffer.data(kg + 105);
    const auto *kg_106 = buffer.data(kg + 106);
    const auto *kg_107 = buffer.data(kg + 107);
    const auto *kg_108 = buffer.data(kg + 108);
    const auto *kg_109 = buffer.data(kg + 109);
    const auto *kg_110 = buffer.data(kg + 110);
    const auto *kg_111 = buffer.data(kg + 111);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_114 = buffer.data(kg + 114);
    const auto *kg_115 = buffer.data(kg + 115);
    const auto *kg_116 = buffer.data(kg + 116);
    const auto *kg_117 = buffer.data(kg + 117);
    const auto *kg_118 = buffer.data(kg + 118);
    const auto *kg_119 = buffer.data(kg + 119);
    const auto *kg_120 = buffer.data(kg + 120);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_123 = buffer.data(kg + 123);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_125 = buffer.data(kg + 125);
    const auto *kg_126 = buffer.data(kg + 126);
    const auto *kg_127 = buffer.data(kg + 127);
    const auto *kg_128 = buffer.data(kg + 128);
    const auto *kg_129 = buffer.data(kg + 129);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_132 = buffer.data(kg + 132);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_134 = buffer.data(kg + 134);
    const auto *kg_150 = buffer.data(kg + 150);
    const auto *kg_151 = buffer.data(kg + 151);
    const auto *kg_152 = buffer.data(kg + 152);
    const auto *kg_153 = buffer.data(kg + 153);
    const auto *kg_154 = buffer.data(kg + 154);
    const auto *kg_155 = buffer.data(kg + 155);
    const auto *kg_156 = buffer.data(kg + 156);
    const auto *kg_157 = buffer.data(kg + 157);
    const auto *kg_158 = buffer.data(kg + 158);
    const auto *kg_159 = buffer.data(kg + 159);
    const auto *kg_160 = buffer.data(kg + 160);
    const auto *kg_161 = buffer.data(kg + 161);
    const auto *kg_162 = buffer.data(kg + 162);
    const auto *kg_163 = buffer.data(kg + 163);
    const auto *kg_164 = buffer.data(kg + 164);
    const auto *kg_165 = buffer.data(kg + 165);
    const auto *kg_166 = buffer.data(kg + 166);
    const auto *kg_167 = buffer.data(kg + 167);
    const auto *kg_168 = buffer.data(kg + 168);
    const auto *kg_169 = buffer.data(kg + 169);
    const auto *kg_170 = buffer.data(kg + 170);
    const auto *kg_171 = buffer.data(kg + 171);
    const auto *kg_172 = buffer.data(kg + 172);
    const auto *kg_173 = buffer.data(kg + 173);
    const auto *kg_174 = buffer.data(kg + 174);
    const auto *kg_175 = buffer.data(kg + 175);
    const auto *kg_176 = buffer.data(kg + 176);
    const auto *kg_177 = buffer.data(kg + 177);
    const auto *kg_178 = buffer.data(kg + 178);
    const auto *kg_179 = buffer.data(kg + 179);
    const auto *kg_180 = buffer.data(kg + 180);
    const auto *kg_181 = buffer.data(kg + 181);
    const auto *kg_182 = buffer.data(kg + 182);
    const auto *kg_183 = buffer.data(kg + 183);
    const auto *kg_184 = buffer.data(kg + 184);
    const auto *kg_185 = buffer.data(kg + 185);
    const auto *kg_186 = buffer.data(kg + 186);
    const auto *kg_187 = buffer.data(kg + 187);
    const auto *kg_188 = buffer.data(kg + 188);
    const auto *kg_189 = buffer.data(kg + 189);
    const auto *kg_190 = buffer.data(kg + 190);
    const auto *kg_191 = buffer.data(kg + 191);
    const auto *kg_192 = buffer.data(kg + 192);
    const auto *kg_193 = buffer.data(kg + 193);
    const auto *kg_194 = buffer.data(kg + 194);
    const auto *kg_195 = buffer.data(kg + 195);
    const auto *kg_196 = buffer.data(kg + 196);
    const auto *kg_197 = buffer.data(kg + 197);
    const auto *kg_198 = buffer.data(kg + 198);
    const auto *kg_199 = buffer.data(kg + 199);
    const auto *kg_200 = buffer.data(kg + 200);
    const auto *kg_201 = buffer.data(kg + 201);
    const auto *kg_202 = buffer.data(kg + 202);
    const auto *kg_203 = buffer.data(kg + 203);
    const auto *kg_204 = buffer.data(kg + 204);
    const auto *kg_205 = buffer.data(kg + 205);
    const auto *kg_206 = buffer.data(kg + 206);
    const auto *kg_207 = buffer.data(kg + 207);
    const auto *kg_208 = buffer.data(kg + 208);
    const auto *kg_209 = buffer.data(kg + 209);
    const auto *kg_225 = buffer.data(kg + 225);
    const auto *kg_226 = buffer.data(kg + 226);
    const auto *kg_227 = buffer.data(kg + 227);
    const auto *kg_228 = buffer.data(kg + 228);
    const auto *kg_229 = buffer.data(kg + 229);
    const auto *kg_230 = buffer.data(kg + 230);
    const auto *kg_231 = buffer.data(kg + 231);
    const auto *kg_232 = buffer.data(kg + 232);
    const auto *kg_233 = buffer.data(kg + 233);
    const auto *kg_234 = buffer.data(kg + 234);
    const auto *kg_235 = buffer.data(kg + 235);
    const auto *kg_236 = buffer.data(kg + 236);
    const auto *kg_237 = buffer.data(kg + 237);
    const auto *kg_238 = buffer.data(kg + 238);
    const auto *kg_239 = buffer.data(kg + 239);
    const auto *kg_240 = buffer.data(kg + 240);
    const auto *kg_241 = buffer.data(kg + 241);
    const auto *kg_242 = buffer.data(kg + 242);
    const auto *kg_243 = buffer.data(kg + 243);
    const auto *kg_244 = buffer.data(kg + 244);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, kg_15, kg_16, kg_17, kg_18, \
                         kg_19, kg_20, kg_21, kg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kg_15[k];

        t_1[k] = f_0 * kg_16[k];

        t_2[k] = f_0 * kg_17[k];

        t_3[k] = f_0 * kg_18[k];

        t_4[k] = f_0 * kg_19[k];

        t_5[k] = f_0 * kg_20[k];

        t_6[k] = f_0 * kg_21[k];

        t_7[k] = f_0 * kg_22[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, kg_23, kg_24, kg_25, kg_26, \
                         kg_27, kg_28, kg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * kg_23[k];

        t_9[k] = f_0 * kg_24[k];

        t_10[k] = f_0 * kg_25[k];

        t_11[k] = f_0 * kg_26[k];

        t_12[k] = f_0 * kg_27[k];

        t_13[k] = f_0 * kg_28[k];

        t_14[k] = f_0 * kg_29[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, hg_0, hg_1, hg_2, hg_3, hg_4, kg_45, \
                         kg_46, kg_47, kg_48, kg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -hg_0[k]
                  + f_0 * kg_45[k];

        t_16[k] = -hg_1[k]
                  + f_0 * kg_46[k];

        t_17[k] = -hg_2[k]
                  + f_0 * kg_47[k];

        t_18[k] = -hg_3[k]
                  + f_0 * kg_48[k];

        t_19[k] = -hg_4[k]
                  + f_0 * kg_49[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, hg_5, hg_6, hg_7, hg_8, hg_9, kg_50, \
                         kg_51, kg_52, kg_53, kg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -hg_5[k]
                  + f_0 * kg_50[k];

        t_21[k] = -hg_6[k]
                  + f_0 * kg_51[k];

        t_22[k] = -hg_7[k]
                  + f_0 * kg_52[k];

        t_23[k] = -hg_8[k]
                  + f_0 * kg_53[k];

        t_24[k] = -hg_9[k]
                  + f_0 * kg_54[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, hg_10, hg_11, hg_12, hg_13, hg_14, \
                         kg_55, kg_56, kg_57, kg_58, kg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -hg_10[k]
                  + f_0 * kg_55[k];

        t_26[k] = -hg_11[k]
                  + f_0 * kg_56[k];

        t_27[k] = -hg_12[k]
                  + f_0 * kg_57[k];

        t_28[k] = -hg_13[k]
                  + f_0 * kg_58[k];

        t_29[k] = -hg_14[k]
                  + f_0 * kg_59[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, t_36, t_37, kg_60, kg_61, kg_62, \
                         kg_63, kg_64, kg_65, kg_66, kg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * kg_60[k];

        t_31[k] = f_0 * kg_61[k];

        t_32[k] = f_0 * kg_62[k];

        t_33[k] = f_0 * kg_63[k];

        t_34[k] = f_0 * kg_64[k];

        t_35[k] = f_0 * kg_65[k];

        t_36[k] = f_0 * kg_66[k];

        t_37[k] = f_0 * kg_67[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, kg_68, kg_69, kg_70, kg_71, \
                         kg_72, kg_73, kg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * kg_68[k];

        t_39[k] = f_0 * kg_69[k];

        t_40[k] = f_0 * kg_70[k];

        t_41[k] = f_0 * kg_71[k];

        t_42[k] = f_0 * kg_72[k];

        t_43[k] = f_0 * kg_73[k];

        t_44[k] = f_0 * kg_74[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hg_15, hg_16, hg_17, hg_18, hg_19, \
                         kg_90, kg_91, kg_92, kg_93, kg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * hg_15[k]
                  + f_0 * kg_90[k];

        t_46[k] = -2.0 * hg_16[k]
                  + f_0 * kg_91[k];

        t_47[k] = -2.0 * hg_17[k]
                  + f_0 * kg_92[k];

        t_48[k] = -2.0 * hg_18[k]
                  + f_0 * kg_93[k];

        t_49[k] = -2.0 * hg_19[k]
                  + f_0 * kg_94[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hg_20, hg_21, hg_22, hg_23, hg_24, \
                         kg_95, kg_96, kg_97, kg_98, kg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * hg_20[k]
                  + f_0 * kg_95[k];

        t_51[k] = -2.0 * hg_21[k]
                  + f_0 * kg_96[k];

        t_52[k] = -2.0 * hg_22[k]
                  + f_0 * kg_97[k];

        t_53[k] = -2.0 * hg_23[k]
                  + f_0 * kg_98[k];

        t_54[k] = -2.0 * hg_24[k]
                  + f_0 * kg_99[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, hg_25, hg_26, hg_27, hg_28, hg_29, \
                         kg_100, kg_101, kg_102, kg_103, kg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * hg_25[k]
                  + f_0 * kg_100[k];

        t_56[k] = -2.0 * hg_26[k]
                  + f_0 * kg_101[k];

        t_57[k] = -2.0 * hg_27[k]
                  + f_0 * kg_102[k];

        t_58[k] = -2.0 * hg_28[k]
                  + f_0 * kg_103[k];

        t_59[k] = -2.0 * hg_29[k]
                  + f_0 * kg_104[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, hg_30, hg_31, hg_32, hg_33, hg_34, \
                         kg_105, kg_106, kg_107, kg_108, kg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -hg_30[k]
                  + f_0 * kg_105[k];

        t_61[k] = -hg_31[k]
                  + f_0 * kg_106[k];

        t_62[k] = -hg_32[k]
                  + f_0 * kg_107[k];

        t_63[k] = -hg_33[k]
                  + f_0 * kg_108[k];

        t_64[k] = -hg_34[k]
                  + f_0 * kg_109[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, hg_35, hg_36, hg_37, hg_38, hg_39, \
                         kg_110, kg_111, kg_112, kg_113, kg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -hg_35[k]
                  + f_0 * kg_110[k];

        t_66[k] = -hg_36[k]
                  + f_0 * kg_111[k];

        t_67[k] = -hg_37[k]
                  + f_0 * kg_112[k];

        t_68[k] = -hg_38[k]
                  + f_0 * kg_113[k];

        t_69[k] = -hg_39[k]
                  + f_0 * kg_114[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, hg_40, hg_41, hg_42, hg_43, hg_44, \
                         kg_115, kg_116, kg_117, kg_118, kg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -hg_40[k]
                  + f_0 * kg_115[k];

        t_71[k] = -hg_41[k]
                  + f_0 * kg_116[k];

        t_72[k] = -hg_42[k]
                  + f_0 * kg_117[k];

        t_73[k] = -hg_43[k]
                  + f_0 * kg_118[k];

        t_74[k] = -hg_44[k]
                  + f_0 * kg_119[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, t_80, t_81, t_82, kg_120, kg_121, \
                         kg_122, kg_123, kg_124, kg_125, kg_126, \
                         kg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * kg_120[k];

        t_76[k] = f_0 * kg_121[k];

        t_77[k] = f_0 * kg_122[k];

        t_78[k] = f_0 * kg_123[k];

        t_79[k] = f_0 * kg_124[k];

        t_80[k] = f_0 * kg_125[k];

        t_81[k] = f_0 * kg_126[k];

        t_82[k] = f_0 * kg_127[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, kg_128, kg_129, kg_130, \
                         kg_131, kg_132, kg_133, kg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * kg_128[k];

        t_84[k] = f_0 * kg_129[k];

        t_85[k] = f_0 * kg_130[k];

        t_86[k] = f_0 * kg_131[k];

        t_87[k] = f_0 * kg_132[k];

        t_88[k] = f_0 * kg_133[k];

        t_89[k] = f_0 * kg_134[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, hg_45, hg_46, hg_47, hg_48, hg_49, \
                         kg_150, kg_151, kg_152, kg_153, kg_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * hg_45[k]
                  + f_0 * kg_150[k];

        t_91[k] = -3.0 * hg_46[k]
                  + f_0 * kg_151[k];

        t_92[k] = -3.0 * hg_47[k]
                  + f_0 * kg_152[k];

        t_93[k] = -3.0 * hg_48[k]
                  + f_0 * kg_153[k];

        t_94[k] = -3.0 * hg_49[k]
                  + f_0 * kg_154[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, hg_50, hg_51, hg_52, hg_53, hg_54, \
                         kg_155, kg_156, kg_157, kg_158, kg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * hg_50[k]
                  + f_0 * kg_155[k];

        t_96[k] = -3.0 * hg_51[k]
                  + f_0 * kg_156[k];

        t_97[k] = -3.0 * hg_52[k]
                  + f_0 * kg_157[k];

        t_98[k] = -3.0 * hg_53[k]
                  + f_0 * kg_158[k];

        t_99[k] = -3.0 * hg_54[k]
                  + f_0 * kg_159[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, hg_55, hg_56, hg_57, hg_58, hg_59, \
                         kg_160, kg_161, kg_162, kg_163, kg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * hg_55[k]
                   + f_0 * kg_160[k];

        t_101[k] = -3.0 * hg_56[k]
                   + f_0 * kg_161[k];

        t_102[k] = -3.0 * hg_57[k]
                   + f_0 * kg_162[k];

        t_103[k] = -3.0 * hg_58[k]
                   + f_0 * kg_163[k];

        t_104[k] = -3.0 * hg_59[k]
                   + f_0 * kg_164[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, hg_60, hg_61, hg_62, hg_63, hg_64, \
                         kg_165, kg_166, kg_167, kg_168, kg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * hg_60[k]
                   + f_0 * kg_165[k];

        t_106[k] = -2.0 * hg_61[k]
                   + f_0 * kg_166[k];

        t_107[k] = -2.0 * hg_62[k]
                   + f_0 * kg_167[k];

        t_108[k] = -2.0 * hg_63[k]
                   + f_0 * kg_168[k];

        t_109[k] = -2.0 * hg_64[k]
                   + f_0 * kg_169[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, hg_65, hg_66, hg_67, hg_68, hg_69, \
                         kg_170, kg_171, kg_172, kg_173, kg_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * hg_65[k]
                   + f_0 * kg_170[k];

        t_111[k] = -2.0 * hg_66[k]
                   + f_0 * kg_171[k];

        t_112[k] = -2.0 * hg_67[k]
                   + f_0 * kg_172[k];

        t_113[k] = -2.0 * hg_68[k]
                   + f_0 * kg_173[k];

        t_114[k] = -2.0 * hg_69[k]
                   + f_0 * kg_174[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, hg_70, hg_71, hg_72, hg_73, hg_74, \
                         kg_175, kg_176, kg_177, kg_178, kg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * hg_70[k]
                   + f_0 * kg_175[k];

        t_116[k] = -2.0 * hg_71[k]
                   + f_0 * kg_176[k];

        t_117[k] = -2.0 * hg_72[k]
                   + f_0 * kg_177[k];

        t_118[k] = -2.0 * hg_73[k]
                   + f_0 * kg_178[k];

        t_119[k] = -2.0 * hg_74[k]
                   + f_0 * kg_179[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, hg_75, hg_76, hg_77, hg_78, hg_79, \
                         kg_180, kg_181, kg_182, kg_183, kg_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -hg_75[k]
                   + f_0 * kg_180[k];

        t_121[k] = -hg_76[k]
                   + f_0 * kg_181[k];

        t_122[k] = -hg_77[k]
                   + f_0 * kg_182[k];

        t_123[k] = -hg_78[k]
                   + f_0 * kg_183[k];

        t_124[k] = -hg_79[k]
                   + f_0 * kg_184[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, hg_80, hg_81, hg_82, hg_83, hg_84, \
                         kg_185, kg_186, kg_187, kg_188, kg_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -hg_80[k]
                   + f_0 * kg_185[k];

        t_126[k] = -hg_81[k]
                   + f_0 * kg_186[k];

        t_127[k] = -hg_82[k]
                   + f_0 * kg_187[k];

        t_128[k] = -hg_83[k]
                   + f_0 * kg_188[k];

        t_129[k] = -hg_84[k]
                   + f_0 * kg_189[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, hg_85, hg_86, hg_87, hg_88, hg_89, \
                         kg_190, kg_191, kg_192, kg_193, kg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -hg_85[k]
                   + f_0 * kg_190[k];

        t_131[k] = -hg_86[k]
                   + f_0 * kg_191[k];

        t_132[k] = -hg_87[k]
                   + f_0 * kg_192[k];

        t_133[k] = -hg_88[k]
                   + f_0 * kg_193[k];

        t_134[k] = -hg_89[k]
                   + f_0 * kg_194[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, t_140, t_141, t_142, kg_195, \
                         kg_196, kg_197, kg_198, kg_199, kg_200, kg_201, \
                         kg_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * kg_195[k];

        t_136[k] = f_0 * kg_196[k];

        t_137[k] = f_0 * kg_197[k];

        t_138[k] = f_0 * kg_198[k];

        t_139[k] = f_0 * kg_199[k];

        t_140[k] = f_0 * kg_200[k];

        t_141[k] = f_0 * kg_201[k];

        t_142[k] = f_0 * kg_202[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, t_148, t_149, kg_203, kg_204, \
                         kg_205, kg_206, kg_207, kg_208, kg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_0 * kg_203[k];

        t_144[k] = f_0 * kg_204[k];

        t_145[k] = f_0 * kg_205[k];

        t_146[k] = f_0 * kg_206[k];

        t_147[k] = f_0 * kg_207[k];

        t_148[k] = f_0 * kg_208[k];

        t_149[k] = f_0 * kg_209[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, hg_90, hg_91, hg_92, hg_93, hg_94, \
                         kg_225, kg_226, kg_227, kg_228, kg_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -4.0 * hg_90[k]
                   + f_0 * kg_225[k];

        t_151[k] = -4.0 * hg_91[k]
                   + f_0 * kg_226[k];

        t_152[k] = -4.0 * hg_92[k]
                   + f_0 * kg_227[k];

        t_153[k] = -4.0 * hg_93[k]
                   + f_0 * kg_228[k];

        t_154[k] = -4.0 * hg_94[k]
                   + f_0 * kg_229[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, hg_95, hg_96, hg_97, hg_98, hg_99, \
                         kg_230, kg_231, kg_232, kg_233, kg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -4.0 * hg_95[k]
                   + f_0 * kg_230[k];

        t_156[k] = -4.0 * hg_96[k]
                   + f_0 * kg_231[k];

        t_157[k] = -4.0 * hg_97[k]
                   + f_0 * kg_232[k];

        t_158[k] = -4.0 * hg_98[k]
                   + f_0 * kg_233[k];

        t_159[k] = -4.0 * hg_99[k]
                   + f_0 * kg_234[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, hg_100, hg_101, hg_102, hg_103, \
                         hg_104, kg_235, kg_236, kg_237, kg_238, \
                         kg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -4.0 * hg_100[k]
                   + f_0 * kg_235[k];

        t_161[k] = -4.0 * hg_101[k]
                   + f_0 * kg_236[k];

        t_162[k] = -4.0 * hg_102[k]
                   + f_0 * kg_237[k];

        t_163[k] = -4.0 * hg_103[k]
                   + f_0 * kg_238[k];

        t_164[k] = -4.0 * hg_104[k]
                   + f_0 * kg_239[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, hg_105, hg_106, hg_107, hg_108, \
                         hg_109, kg_240, kg_241, kg_242, kg_243, \
                         kg_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -3.0 * hg_105[k]
                   + f_0 * kg_240[k];

        t_166[k] = -3.0 * hg_106[k]
                   + f_0 * kg_241[k];

        t_167[k] = -3.0 * hg_107[k]
                   + f_0 * kg_242[k];

        t_168[k] = -3.0 * hg_108[k]
                   + f_0 * kg_243[k];

        t_169[k] = -3.0 * hg_109[k]
                   + f_0 * kg_244[k];
    }
}

static auto
compute_prim_geom_10_ig_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hg, const size_t kg,
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

    const auto *kg_245 = buffer.data(kg + 245);
    const auto *kg_246 = buffer.data(kg + 246);
    const auto *kg_247 = buffer.data(kg + 247);
    const auto *kg_248 = buffer.data(kg + 248);
    const auto *kg_249 = buffer.data(kg + 249);
    const auto *kg_250 = buffer.data(kg + 250);
    const auto *kg_251 = buffer.data(kg + 251);
    const auto *kg_252 = buffer.data(kg + 252);
    const auto *kg_253 = buffer.data(kg + 253);
    const auto *kg_254 = buffer.data(kg + 254);
    const auto *kg_255 = buffer.data(kg + 255);
    const auto *kg_256 = buffer.data(kg + 256);
    const auto *kg_257 = buffer.data(kg + 257);
    const auto *kg_258 = buffer.data(kg + 258);
    const auto *kg_259 = buffer.data(kg + 259);
    const auto *kg_260 = buffer.data(kg + 260);
    const auto *kg_261 = buffer.data(kg + 261);
    const auto *kg_262 = buffer.data(kg + 262);
    const auto *kg_263 = buffer.data(kg + 263);
    const auto *kg_264 = buffer.data(kg + 264);
    const auto *kg_265 = buffer.data(kg + 265);
    const auto *kg_266 = buffer.data(kg + 266);
    const auto *kg_267 = buffer.data(kg + 267);
    const auto *kg_268 = buffer.data(kg + 268);
    const auto *kg_269 = buffer.data(kg + 269);
    const auto *kg_270 = buffer.data(kg + 270);
    const auto *kg_271 = buffer.data(kg + 271);
    const auto *kg_272 = buffer.data(kg + 272);
    const auto *kg_273 = buffer.data(kg + 273);
    const auto *kg_274 = buffer.data(kg + 274);
    const auto *kg_275 = buffer.data(kg + 275);
    const auto *kg_276 = buffer.data(kg + 276);
    const auto *kg_277 = buffer.data(kg + 277);
    const auto *kg_278 = buffer.data(kg + 278);
    const auto *kg_279 = buffer.data(kg + 279);
    const auto *kg_280 = buffer.data(kg + 280);
    const auto *kg_281 = buffer.data(kg + 281);
    const auto *kg_282 = buffer.data(kg + 282);
    const auto *kg_283 = buffer.data(kg + 283);
    const auto *kg_284 = buffer.data(kg + 284);
    const auto *kg_285 = buffer.data(kg + 285);
    const auto *kg_286 = buffer.data(kg + 286);
    const auto *kg_287 = buffer.data(kg + 287);
    const auto *kg_288 = buffer.data(kg + 288);
    const auto *kg_289 = buffer.data(kg + 289);
    const auto *kg_290 = buffer.data(kg + 290);
    const auto *kg_291 = buffer.data(kg + 291);
    const auto *kg_292 = buffer.data(kg + 292);
    const auto *kg_293 = buffer.data(kg + 293);
    const auto *kg_294 = buffer.data(kg + 294);
    const auto *kg_295 = buffer.data(kg + 295);
    const auto *kg_296 = buffer.data(kg + 296);
    const auto *kg_297 = buffer.data(kg + 297);
    const auto *kg_298 = buffer.data(kg + 298);
    const auto *kg_299 = buffer.data(kg + 299);
    const auto *kg_315 = buffer.data(kg + 315);
    const auto *kg_316 = buffer.data(kg + 316);
    const auto *kg_317 = buffer.data(kg + 317);
    const auto *kg_318 = buffer.data(kg + 318);
    const auto *kg_319 = buffer.data(kg + 319);
    const auto *kg_320 = buffer.data(kg + 320);
    const auto *kg_321 = buffer.data(kg + 321);
    const auto *kg_322 = buffer.data(kg + 322);
    const auto *kg_323 = buffer.data(kg + 323);
    const auto *kg_324 = buffer.data(kg + 324);
    const auto *kg_325 = buffer.data(kg + 325);
    const auto *kg_326 = buffer.data(kg + 326);
    const auto *kg_327 = buffer.data(kg + 327);
    const auto *kg_328 = buffer.data(kg + 328);
    const auto *kg_329 = buffer.data(kg + 329);
    const auto *kg_330 = buffer.data(kg + 330);
    const auto *kg_331 = buffer.data(kg + 331);
    const auto *kg_332 = buffer.data(kg + 332);
    const auto *kg_333 = buffer.data(kg + 333);
    const auto *kg_334 = buffer.data(kg + 334);
    const auto *kg_335 = buffer.data(kg + 335);
    const auto *kg_336 = buffer.data(kg + 336);
    const auto *kg_337 = buffer.data(kg + 337);
    const auto *kg_338 = buffer.data(kg + 338);
    const auto *kg_339 = buffer.data(kg + 339);
    const auto *kg_340 = buffer.data(kg + 340);
    const auto *kg_341 = buffer.data(kg + 341);
    const auto *kg_342 = buffer.data(kg + 342);
    const auto *kg_343 = buffer.data(kg + 343);
    const auto *kg_344 = buffer.data(kg + 344);
    const auto *kg_345 = buffer.data(kg + 345);
    const auto *kg_346 = buffer.data(kg + 346);
    const auto *kg_347 = buffer.data(kg + 347);
    const auto *kg_348 = buffer.data(kg + 348);
    const auto *kg_349 = buffer.data(kg + 349);
    const auto *kg_350 = buffer.data(kg + 350);
    const auto *kg_351 = buffer.data(kg + 351);
    const auto *kg_352 = buffer.data(kg + 352);
    const auto *kg_353 = buffer.data(kg + 353);
    const auto *kg_354 = buffer.data(kg + 354);
    const auto *kg_355 = buffer.data(kg + 355);
    const auto *kg_356 = buffer.data(kg + 356);
    const auto *kg_357 = buffer.data(kg + 357);
    const auto *kg_358 = buffer.data(kg + 358);
    const auto *kg_359 = buffer.data(kg + 359);
    const auto *kg_360 = buffer.data(kg + 360);
    const auto *kg_361 = buffer.data(kg + 361);
    const auto *kg_362 = buffer.data(kg + 362);
    const auto *kg_363 = buffer.data(kg + 363);
    const auto *kg_364 = buffer.data(kg + 364);
    const auto *kg_365 = buffer.data(kg + 365);
    const auto *kg_366 = buffer.data(kg + 366);
    const auto *kg_367 = buffer.data(kg + 367);
    const auto *kg_368 = buffer.data(kg + 368);
    const auto *kg_369 = buffer.data(kg + 369);
    const auto *kg_370 = buffer.data(kg + 370);
    const auto *kg_371 = buffer.data(kg + 371);
    const auto *kg_372 = buffer.data(kg + 372);
    const auto *kg_373 = buffer.data(kg + 373);
    const auto *kg_374 = buffer.data(kg + 374);
    const auto *kg_375 = buffer.data(kg + 375);
    const auto *kg_376 = buffer.data(kg + 376);
    const auto *kg_377 = buffer.data(kg + 377);
    const auto *kg_378 = buffer.data(kg + 378);
    const auto *kg_379 = buffer.data(kg + 379);
    const auto *kg_380 = buffer.data(kg + 380);
    const auto *kg_381 = buffer.data(kg + 381);
    const auto *kg_382 = buffer.data(kg + 382);
    const auto *kg_383 = buffer.data(kg + 383);
    const auto *kg_384 = buffer.data(kg + 384);
    const auto *kg_385 = buffer.data(kg + 385);
    const auto *kg_386 = buffer.data(kg + 386);
    const auto *kg_387 = buffer.data(kg + 387);
    const auto *kg_388 = buffer.data(kg + 388);
    const auto *kg_389 = buffer.data(kg + 389);
    const auto *kg_390 = buffer.data(kg + 390);
    const auto *kg_391 = buffer.data(kg + 391);
    const auto *kg_392 = buffer.data(kg + 392);
    const auto *kg_393 = buffer.data(kg + 393);
    const auto *kg_394 = buffer.data(kg + 394);
    const auto *kg_395 = buffer.data(kg + 395);
    const auto *kg_396 = buffer.data(kg + 396);
    const auto *kg_397 = buffer.data(kg + 397);
    const auto *kg_398 = buffer.data(kg + 398);
    const auto *kg_399 = buffer.data(kg + 399);
    const auto *kg_400 = buffer.data(kg + 400);
    const auto *kg_401 = buffer.data(kg + 401);
    const auto *kg_402 = buffer.data(kg + 402);
    const auto *kg_403 = buffer.data(kg + 403);
    const auto *kg_404 = buffer.data(kg + 404);
    const auto *kg_420 = buffer.data(kg + 420);
    const auto *kg_421 = buffer.data(kg + 421);
    const auto *kg_422 = buffer.data(kg + 422);
    const auto *kg_423 = buffer.data(kg + 423);
    const auto *kg_424 = buffer.data(kg + 424);
    const auto *kg_425 = buffer.data(kg + 425);
    const auto *kg_426 = buffer.data(kg + 426);
    const auto *kg_427 = buffer.data(kg + 427);
    const auto *kg_428 = buffer.data(kg + 428);
    const auto *kg_429 = buffer.data(kg + 429);
    const auto *kg_430 = buffer.data(kg + 430);
    const auto *kg_431 = buffer.data(kg + 431);
    const auto *kg_432 = buffer.data(kg + 432);
    const auto *kg_433 = buffer.data(kg + 433);
    const auto *kg_434 = buffer.data(kg + 434);

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, hg_110, hg_111, hg_112, hg_113, \
                         hg_114, kg_245, kg_246, kg_247, kg_248, \
                         kg_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * hg_110[k]
                   + f_0 * kg_245[k];

        t_171[k] = -3.0 * hg_111[k]
                   + f_0 * kg_246[k];

        t_172[k] = -3.0 * hg_112[k]
                   + f_0 * kg_247[k];

        t_173[k] = -3.0 * hg_113[k]
                   + f_0 * kg_248[k];

        t_174[k] = -3.0 * hg_114[k]
                   + f_0 * kg_249[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, hg_115, hg_116, hg_117, hg_118, \
                         hg_119, kg_250, kg_251, kg_252, kg_253, \
                         kg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * hg_115[k]
                   + f_0 * kg_250[k];

        t_176[k] = -3.0 * hg_116[k]
                   + f_0 * kg_251[k];

        t_177[k] = -3.0 * hg_117[k]
                   + f_0 * kg_252[k];

        t_178[k] = -3.0 * hg_118[k]
                   + f_0 * kg_253[k];

        t_179[k] = -3.0 * hg_119[k]
                   + f_0 * kg_254[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, hg_120, hg_121, hg_122, hg_123, \
                         hg_124, kg_255, kg_256, kg_257, kg_258, \
                         kg_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * hg_120[k]
                   + f_0 * kg_255[k];

        t_181[k] = -2.0 * hg_121[k]
                   + f_0 * kg_256[k];

        t_182[k] = -2.0 * hg_122[k]
                   + f_0 * kg_257[k];

        t_183[k] = -2.0 * hg_123[k]
                   + f_0 * kg_258[k];

        t_184[k] = -2.0 * hg_124[k]
                   + f_0 * kg_259[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, hg_125, hg_126, hg_127, hg_128, \
                         hg_129, kg_260, kg_261, kg_262, kg_263, \
                         kg_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * hg_125[k]
                   + f_0 * kg_260[k];

        t_186[k] = -2.0 * hg_126[k]
                   + f_0 * kg_261[k];

        t_187[k] = -2.0 * hg_127[k]
                   + f_0 * kg_262[k];

        t_188[k] = -2.0 * hg_128[k]
                   + f_0 * kg_263[k];

        t_189[k] = -2.0 * hg_129[k]
                   + f_0 * kg_264[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, hg_130, hg_131, hg_132, hg_133, \
                         hg_134, kg_265, kg_266, kg_267, kg_268, \
                         kg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * hg_130[k]
                   + f_0 * kg_265[k];

        t_191[k] = -2.0 * hg_131[k]
                   + f_0 * kg_266[k];

        t_192[k] = -2.0 * hg_132[k]
                   + f_0 * kg_267[k];

        t_193[k] = -2.0 * hg_133[k]
                   + f_0 * kg_268[k];

        t_194[k] = -2.0 * hg_134[k]
                   + f_0 * kg_269[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, hg_135, hg_136, hg_137, hg_138, \
                         hg_139, kg_270, kg_271, kg_272, kg_273, \
                         kg_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -hg_135[k]
                   + f_0 * kg_270[k];

        t_196[k] = -hg_136[k]
                   + f_0 * kg_271[k];

        t_197[k] = -hg_137[k]
                   + f_0 * kg_272[k];

        t_198[k] = -hg_138[k]
                   + f_0 * kg_273[k];

        t_199[k] = -hg_139[k]
                   + f_0 * kg_274[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, hg_140, hg_141, hg_142, hg_143, \
                         hg_144, kg_275, kg_276, kg_277, kg_278, \
                         kg_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -hg_140[k]
                   + f_0 * kg_275[k];

        t_201[k] = -hg_141[k]
                   + f_0 * kg_276[k];

        t_202[k] = -hg_142[k]
                   + f_0 * kg_277[k];

        t_203[k] = -hg_143[k]
                   + f_0 * kg_278[k];

        t_204[k] = -hg_144[k]
                   + f_0 * kg_279[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, hg_145, hg_146, hg_147, hg_148, \
                         hg_149, kg_280, kg_281, kg_282, kg_283, \
                         kg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -hg_145[k]
                   + f_0 * kg_280[k];

        t_206[k] = -hg_146[k]
                   + f_0 * kg_281[k];

        t_207[k] = -hg_147[k]
                   + f_0 * kg_282[k];

        t_208[k] = -hg_148[k]
                   + f_0 * kg_283[k];

        t_209[k] = -hg_149[k]
                   + f_0 * kg_284[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, t_217, kg_285, \
                         kg_286, kg_287, kg_288, kg_289, kg_290, kg_291, \
                         kg_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * kg_285[k];

        t_211[k] = f_0 * kg_286[k];

        t_212[k] = f_0 * kg_287[k];

        t_213[k] = f_0 * kg_288[k];

        t_214[k] = f_0 * kg_289[k];

        t_215[k] = f_0 * kg_290[k];

        t_216[k] = f_0 * kg_291[k];

        t_217[k] = f_0 * kg_292[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, t_224, kg_293, kg_294, \
                         kg_295, kg_296, kg_297, kg_298, kg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * kg_293[k];

        t_219[k] = f_0 * kg_294[k];

        t_220[k] = f_0 * kg_295[k];

        t_221[k] = f_0 * kg_296[k];

        t_222[k] = f_0 * kg_297[k];

        t_223[k] = f_0 * kg_298[k];

        t_224[k] = f_0 * kg_299[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, hg_150, hg_151, hg_152, hg_153, \
                         hg_154, kg_315, kg_316, kg_317, kg_318, \
                         kg_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -5.0 * hg_150[k]
                   + f_0 * kg_315[k];

        t_226[k] = -5.0 * hg_151[k]
                   + f_0 * kg_316[k];

        t_227[k] = -5.0 * hg_152[k]
                   + f_0 * kg_317[k];

        t_228[k] = -5.0 * hg_153[k]
                   + f_0 * kg_318[k];

        t_229[k] = -5.0 * hg_154[k]
                   + f_0 * kg_319[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, hg_155, hg_156, hg_157, hg_158, \
                         hg_159, kg_320, kg_321, kg_322, kg_323, \
                         kg_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -5.0 * hg_155[k]
                   + f_0 * kg_320[k];

        t_231[k] = -5.0 * hg_156[k]
                   + f_0 * kg_321[k];

        t_232[k] = -5.0 * hg_157[k]
                   + f_0 * kg_322[k];

        t_233[k] = -5.0 * hg_158[k]
                   + f_0 * kg_323[k];

        t_234[k] = -5.0 * hg_159[k]
                   + f_0 * kg_324[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, hg_160, hg_161, hg_162, hg_163, \
                         hg_164, kg_325, kg_326, kg_327, kg_328, \
                         kg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -5.0 * hg_160[k]
                   + f_0 * kg_325[k];

        t_236[k] = -5.0 * hg_161[k]
                   + f_0 * kg_326[k];

        t_237[k] = -5.0 * hg_162[k]
                   + f_0 * kg_327[k];

        t_238[k] = -5.0 * hg_163[k]
                   + f_0 * kg_328[k];

        t_239[k] = -5.0 * hg_164[k]
                   + f_0 * kg_329[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, hg_165, hg_166, hg_167, hg_168, \
                         hg_169, kg_330, kg_331, kg_332, kg_333, \
                         kg_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -4.0 * hg_165[k]
                   + f_0 * kg_330[k];

        t_241[k] = -4.0 * hg_166[k]
                   + f_0 * kg_331[k];

        t_242[k] = -4.0 * hg_167[k]
                   + f_0 * kg_332[k];

        t_243[k] = -4.0 * hg_168[k]
                   + f_0 * kg_333[k];

        t_244[k] = -4.0 * hg_169[k]
                   + f_0 * kg_334[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, hg_170, hg_171, hg_172, hg_173, \
                         hg_174, kg_335, kg_336, kg_337, kg_338, \
                         kg_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -4.0 * hg_170[k]
                   + f_0 * kg_335[k];

        t_246[k] = -4.0 * hg_171[k]
                   + f_0 * kg_336[k];

        t_247[k] = -4.0 * hg_172[k]
                   + f_0 * kg_337[k];

        t_248[k] = -4.0 * hg_173[k]
                   + f_0 * kg_338[k];

        t_249[k] = -4.0 * hg_174[k]
                   + f_0 * kg_339[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, hg_175, hg_176, hg_177, hg_178, \
                         hg_179, kg_340, kg_341, kg_342, kg_343, \
                         kg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -4.0 * hg_175[k]
                   + f_0 * kg_340[k];

        t_251[k] = -4.0 * hg_176[k]
                   + f_0 * kg_341[k];

        t_252[k] = -4.0 * hg_177[k]
                   + f_0 * kg_342[k];

        t_253[k] = -4.0 * hg_178[k]
                   + f_0 * kg_343[k];

        t_254[k] = -4.0 * hg_179[k]
                   + f_0 * kg_344[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, hg_180, hg_181, hg_182, hg_183, \
                         hg_184, kg_345, kg_346, kg_347, kg_348, \
                         kg_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -3.0 * hg_180[k]
                   + f_0 * kg_345[k];

        t_256[k] = -3.0 * hg_181[k]
                   + f_0 * kg_346[k];

        t_257[k] = -3.0 * hg_182[k]
                   + f_0 * kg_347[k];

        t_258[k] = -3.0 * hg_183[k]
                   + f_0 * kg_348[k];

        t_259[k] = -3.0 * hg_184[k]
                   + f_0 * kg_349[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, hg_185, hg_186, hg_187, hg_188, \
                         hg_189, kg_350, kg_351, kg_352, kg_353, \
                         kg_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -3.0 * hg_185[k]
                   + f_0 * kg_350[k];

        t_261[k] = -3.0 * hg_186[k]
                   + f_0 * kg_351[k];

        t_262[k] = -3.0 * hg_187[k]
                   + f_0 * kg_352[k];

        t_263[k] = -3.0 * hg_188[k]
                   + f_0 * kg_353[k];

        t_264[k] = -3.0 * hg_189[k]
                   + f_0 * kg_354[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, hg_190, hg_191, hg_192, hg_193, \
                         hg_194, kg_355, kg_356, kg_357, kg_358, \
                         kg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -3.0 * hg_190[k]
                   + f_0 * kg_355[k];

        t_266[k] = -3.0 * hg_191[k]
                   + f_0 * kg_356[k];

        t_267[k] = -3.0 * hg_192[k]
                   + f_0 * kg_357[k];

        t_268[k] = -3.0 * hg_193[k]
                   + f_0 * kg_358[k];

        t_269[k] = -3.0 * hg_194[k]
                   + f_0 * kg_359[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, hg_195, hg_196, hg_197, hg_198, \
                         hg_199, kg_360, kg_361, kg_362, kg_363, \
                         kg_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * hg_195[k]
                   + f_0 * kg_360[k];

        t_271[k] = -2.0 * hg_196[k]
                   + f_0 * kg_361[k];

        t_272[k] = -2.0 * hg_197[k]
                   + f_0 * kg_362[k];

        t_273[k] = -2.0 * hg_198[k]
                   + f_0 * kg_363[k];

        t_274[k] = -2.0 * hg_199[k]
                   + f_0 * kg_364[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, hg_200, hg_201, hg_202, hg_203, \
                         hg_204, kg_365, kg_366, kg_367, kg_368, \
                         kg_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * hg_200[k]
                   + f_0 * kg_365[k];

        t_276[k] = -2.0 * hg_201[k]
                   + f_0 * kg_366[k];

        t_277[k] = -2.0 * hg_202[k]
                   + f_0 * kg_367[k];

        t_278[k] = -2.0 * hg_203[k]
                   + f_0 * kg_368[k];

        t_279[k] = -2.0 * hg_204[k]
                   + f_0 * kg_369[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, hg_205, hg_206, hg_207, hg_208, \
                         hg_209, kg_370, kg_371, kg_372, kg_373, \
                         kg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * hg_205[k]
                   + f_0 * kg_370[k];

        t_281[k] = -2.0 * hg_206[k]
                   + f_0 * kg_371[k];

        t_282[k] = -2.0 * hg_207[k]
                   + f_0 * kg_372[k];

        t_283[k] = -2.0 * hg_208[k]
                   + f_0 * kg_373[k];

        t_284[k] = -2.0 * hg_209[k]
                   + f_0 * kg_374[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, hg_210, hg_211, hg_212, hg_213, \
                         hg_214, kg_375, kg_376, kg_377, kg_378, \
                         kg_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -hg_210[k]
                   + f_0 * kg_375[k];

        t_286[k] = -hg_211[k]
                   + f_0 * kg_376[k];

        t_287[k] = -hg_212[k]
                   + f_0 * kg_377[k];

        t_288[k] = -hg_213[k]
                   + f_0 * kg_378[k];

        t_289[k] = -hg_214[k]
                   + f_0 * kg_379[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, hg_215, hg_216, hg_217, hg_218, \
                         hg_219, kg_380, kg_381, kg_382, kg_383, \
                         kg_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -hg_215[k]
                   + f_0 * kg_380[k];

        t_291[k] = -hg_216[k]
                   + f_0 * kg_381[k];

        t_292[k] = -hg_217[k]
                   + f_0 * kg_382[k];

        t_293[k] = -hg_218[k]
                   + f_0 * kg_383[k];

        t_294[k] = -hg_219[k]
                   + f_0 * kg_384[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, hg_220, hg_221, hg_222, hg_223, \
                         hg_224, kg_385, kg_386, kg_387, kg_388, \
                         kg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -hg_220[k]
                   + f_0 * kg_385[k];

        t_296[k] = -hg_221[k]
                   + f_0 * kg_386[k];

        t_297[k] = -hg_222[k]
                   + f_0 * kg_387[k];

        t_298[k] = -hg_223[k]
                   + f_0 * kg_388[k];

        t_299[k] = -hg_224[k]
                   + f_0 * kg_389[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, t_305, t_306, t_307, kg_390, \
                         kg_391, kg_392, kg_393, kg_394, kg_395, kg_396, \
                         kg_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_0 * kg_390[k];

        t_301[k] = f_0 * kg_391[k];

        t_302[k] = f_0 * kg_392[k];

        t_303[k] = f_0 * kg_393[k];

        t_304[k] = f_0 * kg_394[k];

        t_305[k] = f_0 * kg_395[k];

        t_306[k] = f_0 * kg_396[k];

        t_307[k] = f_0 * kg_397[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, t_313, t_314, kg_398, kg_399, \
                         kg_400, kg_401, kg_402, kg_403, kg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_0 * kg_398[k];

        t_309[k] = f_0 * kg_399[k];

        t_310[k] = f_0 * kg_400[k];

        t_311[k] = f_0 * kg_401[k];

        t_312[k] = f_0 * kg_402[k];

        t_313[k] = f_0 * kg_403[k];

        t_314[k] = f_0 * kg_404[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, hg_225, hg_226, hg_227, hg_228, \
                         hg_229, kg_420, kg_421, kg_422, kg_423, \
                         kg_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -6.0 * hg_225[k]
                   + f_0 * kg_420[k];

        t_316[k] = -6.0 * hg_226[k]
                   + f_0 * kg_421[k];

        t_317[k] = -6.0 * hg_227[k]
                   + f_0 * kg_422[k];

        t_318[k] = -6.0 * hg_228[k]
                   + f_0 * kg_423[k];

        t_319[k] = -6.0 * hg_229[k]
                   + f_0 * kg_424[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, hg_230, hg_231, hg_232, hg_233, \
                         hg_234, kg_425, kg_426, kg_427, kg_428, \
                         kg_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -6.0 * hg_230[k]
                   + f_0 * kg_425[k];

        t_321[k] = -6.0 * hg_231[k]
                   + f_0 * kg_426[k];

        t_322[k] = -6.0 * hg_232[k]
                   + f_0 * kg_427[k];

        t_323[k] = -6.0 * hg_233[k]
                   + f_0 * kg_428[k];

        t_324[k] = -6.0 * hg_234[k]
                   + f_0 * kg_429[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, hg_235, hg_236, hg_237, hg_238, \
                         hg_239, kg_430, kg_431, kg_432, kg_433, \
                         kg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -6.0 * hg_235[k]
                   + f_0 * kg_430[k];

        t_326[k] = -6.0 * hg_236[k]
                   + f_0 * kg_431[k];

        t_327[k] = -6.0 * hg_237[k]
                   + f_0 * kg_432[k];

        t_328[k] = -6.0 * hg_238[k]
                   + f_0 * kg_433[k];

        t_329[k] = -6.0 * hg_239[k]
                   + f_0 * kg_434[k];
    }
}

static auto
compute_prim_geom_10_ig_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hg, const size_t kg,
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

    const auto *kg_435 = buffer.data(kg + 435);
    const auto *kg_436 = buffer.data(kg + 436);
    const auto *kg_437 = buffer.data(kg + 437);
    const auto *kg_438 = buffer.data(kg + 438);
    const auto *kg_439 = buffer.data(kg + 439);
    const auto *kg_440 = buffer.data(kg + 440);
    const auto *kg_441 = buffer.data(kg + 441);
    const auto *kg_442 = buffer.data(kg + 442);
    const auto *kg_443 = buffer.data(kg + 443);
    const auto *kg_444 = buffer.data(kg + 444);
    const auto *kg_445 = buffer.data(kg + 445);
    const auto *kg_446 = buffer.data(kg + 446);
    const auto *kg_447 = buffer.data(kg + 447);
    const auto *kg_448 = buffer.data(kg + 448);
    const auto *kg_449 = buffer.data(kg + 449);
    const auto *kg_450 = buffer.data(kg + 450);
    const auto *kg_451 = buffer.data(kg + 451);
    const auto *kg_452 = buffer.data(kg + 452);
    const auto *kg_453 = buffer.data(kg + 453);
    const auto *kg_454 = buffer.data(kg + 454);
    const auto *kg_455 = buffer.data(kg + 455);
    const auto *kg_456 = buffer.data(kg + 456);
    const auto *kg_457 = buffer.data(kg + 457);
    const auto *kg_458 = buffer.data(kg + 458);
    const auto *kg_459 = buffer.data(kg + 459);
    const auto *kg_460 = buffer.data(kg + 460);
    const auto *kg_461 = buffer.data(kg + 461);
    const auto *kg_462 = buffer.data(kg + 462);
    const auto *kg_463 = buffer.data(kg + 463);
    const auto *kg_464 = buffer.data(kg + 464);
    const auto *kg_465 = buffer.data(kg + 465);
    const auto *kg_466 = buffer.data(kg + 466);
    const auto *kg_467 = buffer.data(kg + 467);
    const auto *kg_468 = buffer.data(kg + 468);
    const auto *kg_469 = buffer.data(kg + 469);
    const auto *kg_470 = buffer.data(kg + 470);
    const auto *kg_471 = buffer.data(kg + 471);
    const auto *kg_472 = buffer.data(kg + 472);
    const auto *kg_473 = buffer.data(kg + 473);
    const auto *kg_474 = buffer.data(kg + 474);
    const auto *kg_475 = buffer.data(kg + 475);
    const auto *kg_476 = buffer.data(kg + 476);
    const auto *kg_477 = buffer.data(kg + 477);
    const auto *kg_478 = buffer.data(kg + 478);
    const auto *kg_479 = buffer.data(kg + 479);
    const auto *kg_480 = buffer.data(kg + 480);
    const auto *kg_481 = buffer.data(kg + 481);
    const auto *kg_482 = buffer.data(kg + 482);
    const auto *kg_483 = buffer.data(kg + 483);
    const auto *kg_484 = buffer.data(kg + 484);
    const auto *kg_485 = buffer.data(kg + 485);
    const auto *kg_486 = buffer.data(kg + 486);
    const auto *kg_487 = buffer.data(kg + 487);
    const auto *kg_488 = buffer.data(kg + 488);
    const auto *kg_489 = buffer.data(kg + 489);
    const auto *kg_490 = buffer.data(kg + 490);
    const auto *kg_491 = buffer.data(kg + 491);
    const auto *kg_492 = buffer.data(kg + 492);
    const auto *kg_493 = buffer.data(kg + 493);
    const auto *kg_494 = buffer.data(kg + 494);
    const auto *kg_495 = buffer.data(kg + 495);
    const auto *kg_496 = buffer.data(kg + 496);
    const auto *kg_497 = buffer.data(kg + 497);
    const auto *kg_498 = buffer.data(kg + 498);
    const auto *kg_499 = buffer.data(kg + 499);
    const auto *kg_500 = buffer.data(kg + 500);
    const auto *kg_501 = buffer.data(kg + 501);
    const auto *kg_502 = buffer.data(kg + 502);
    const auto *kg_503 = buffer.data(kg + 503);
    const auto *kg_504 = buffer.data(kg + 504);
    const auto *kg_505 = buffer.data(kg + 505);
    const auto *kg_506 = buffer.data(kg + 506);
    const auto *kg_507 = buffer.data(kg + 507);
    const auto *kg_508 = buffer.data(kg + 508);
    const auto *kg_509 = buffer.data(kg + 509);
    const auto *kg_510 = buffer.data(kg + 510);
    const auto *kg_511 = buffer.data(kg + 511);
    const auto *kg_512 = buffer.data(kg + 512);
    const auto *kg_513 = buffer.data(kg + 513);
    const auto *kg_514 = buffer.data(kg + 514);
    const auto *kg_515 = buffer.data(kg + 515);
    const auto *kg_516 = buffer.data(kg + 516);
    const auto *kg_517 = buffer.data(kg + 517);
    const auto *kg_518 = buffer.data(kg + 518);
    const auto *kg_519 = buffer.data(kg + 519);
    const auto *kg_520 = buffer.data(kg + 520);
    const auto *kg_521 = buffer.data(kg + 521);
    const auto *kg_522 = buffer.data(kg + 522);
    const auto *kg_523 = buffer.data(kg + 523);
    const auto *kg_524 = buffer.data(kg + 524);

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, hg_240, hg_241, hg_242, hg_243, \
                         hg_244, kg_435, kg_436, kg_437, kg_438, \
                         kg_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -5.0 * hg_240[k]
                   + f_0 * kg_435[k];

        t_331[k] = -5.0 * hg_241[k]
                   + f_0 * kg_436[k];

        t_332[k] = -5.0 * hg_242[k]
                   + f_0 * kg_437[k];

        t_333[k] = -5.0 * hg_243[k]
                   + f_0 * kg_438[k];

        t_334[k] = -5.0 * hg_244[k]
                   + f_0 * kg_439[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, hg_245, hg_246, hg_247, hg_248, \
                         hg_249, kg_440, kg_441, kg_442, kg_443, \
                         kg_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -5.0 * hg_245[k]
                   + f_0 * kg_440[k];

        t_336[k] = -5.0 * hg_246[k]
                   + f_0 * kg_441[k];

        t_337[k] = -5.0 * hg_247[k]
                   + f_0 * kg_442[k];

        t_338[k] = -5.0 * hg_248[k]
                   + f_0 * kg_443[k];

        t_339[k] = -5.0 * hg_249[k]
                   + f_0 * kg_444[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, hg_250, hg_251, hg_252, hg_253, \
                         hg_254, kg_445, kg_446, kg_447, kg_448, \
                         kg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -5.0 * hg_250[k]
                   + f_0 * kg_445[k];

        t_341[k] = -5.0 * hg_251[k]
                   + f_0 * kg_446[k];

        t_342[k] = -5.0 * hg_252[k]
                   + f_0 * kg_447[k];

        t_343[k] = -5.0 * hg_253[k]
                   + f_0 * kg_448[k];

        t_344[k] = -5.0 * hg_254[k]
                   + f_0 * kg_449[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, hg_255, hg_256, hg_257, hg_258, \
                         hg_259, kg_450, kg_451, kg_452, kg_453, \
                         kg_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -4.0 * hg_255[k]
                   + f_0 * kg_450[k];

        t_346[k] = -4.0 * hg_256[k]
                   + f_0 * kg_451[k];

        t_347[k] = -4.0 * hg_257[k]
                   + f_0 * kg_452[k];

        t_348[k] = -4.0 * hg_258[k]
                   + f_0 * kg_453[k];

        t_349[k] = -4.0 * hg_259[k]
                   + f_0 * kg_454[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, hg_260, hg_261, hg_262, hg_263, \
                         hg_264, kg_455, kg_456, kg_457, kg_458, \
                         kg_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -4.0 * hg_260[k]
                   + f_0 * kg_455[k];

        t_351[k] = -4.0 * hg_261[k]
                   + f_0 * kg_456[k];

        t_352[k] = -4.0 * hg_262[k]
                   + f_0 * kg_457[k];

        t_353[k] = -4.0 * hg_263[k]
                   + f_0 * kg_458[k];

        t_354[k] = -4.0 * hg_264[k]
                   + f_0 * kg_459[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, hg_265, hg_266, hg_267, hg_268, \
                         hg_269, kg_460, kg_461, kg_462, kg_463, \
                         kg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -4.0 * hg_265[k]
                   + f_0 * kg_460[k];

        t_356[k] = -4.0 * hg_266[k]
                   + f_0 * kg_461[k];

        t_357[k] = -4.0 * hg_267[k]
                   + f_0 * kg_462[k];

        t_358[k] = -4.0 * hg_268[k]
                   + f_0 * kg_463[k];

        t_359[k] = -4.0 * hg_269[k]
                   + f_0 * kg_464[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, hg_270, hg_271, hg_272, hg_273, \
                         hg_274, kg_465, kg_466, kg_467, kg_468, \
                         kg_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -3.0 * hg_270[k]
                   + f_0 * kg_465[k];

        t_361[k] = -3.0 * hg_271[k]
                   + f_0 * kg_466[k];

        t_362[k] = -3.0 * hg_272[k]
                   + f_0 * kg_467[k];

        t_363[k] = -3.0 * hg_273[k]
                   + f_0 * kg_468[k];

        t_364[k] = -3.0 * hg_274[k]
                   + f_0 * kg_469[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, hg_275, hg_276, hg_277, hg_278, \
                         hg_279, kg_470, kg_471, kg_472, kg_473, \
                         kg_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -3.0 * hg_275[k]
                   + f_0 * kg_470[k];

        t_366[k] = -3.0 * hg_276[k]
                   + f_0 * kg_471[k];

        t_367[k] = -3.0 * hg_277[k]
                   + f_0 * kg_472[k];

        t_368[k] = -3.0 * hg_278[k]
                   + f_0 * kg_473[k];

        t_369[k] = -3.0 * hg_279[k]
                   + f_0 * kg_474[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, hg_280, hg_281, hg_282, hg_283, \
                         hg_284, kg_475, kg_476, kg_477, kg_478, \
                         kg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -3.0 * hg_280[k]
                   + f_0 * kg_475[k];

        t_371[k] = -3.0 * hg_281[k]
                   + f_0 * kg_476[k];

        t_372[k] = -3.0 * hg_282[k]
                   + f_0 * kg_477[k];

        t_373[k] = -3.0 * hg_283[k]
                   + f_0 * kg_478[k];

        t_374[k] = -3.0 * hg_284[k]
                   + f_0 * kg_479[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, hg_285, hg_286, hg_287, hg_288, \
                         hg_289, kg_480, kg_481, kg_482, kg_483, \
                         kg_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -2.0 * hg_285[k]
                   + f_0 * kg_480[k];

        t_376[k] = -2.0 * hg_286[k]
                   + f_0 * kg_481[k];

        t_377[k] = -2.0 * hg_287[k]
                   + f_0 * kg_482[k];

        t_378[k] = -2.0 * hg_288[k]
                   + f_0 * kg_483[k];

        t_379[k] = -2.0 * hg_289[k]
                   + f_0 * kg_484[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, hg_290, hg_291, hg_292, hg_293, \
                         hg_294, kg_485, kg_486, kg_487, kg_488, \
                         kg_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * hg_290[k]
                   + f_0 * kg_485[k];

        t_381[k] = -2.0 * hg_291[k]
                   + f_0 * kg_486[k];

        t_382[k] = -2.0 * hg_292[k]
                   + f_0 * kg_487[k];

        t_383[k] = -2.0 * hg_293[k]
                   + f_0 * kg_488[k];

        t_384[k] = -2.0 * hg_294[k]
                   + f_0 * kg_489[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, hg_295, hg_296, hg_297, hg_298, \
                         hg_299, kg_490, kg_491, kg_492, kg_493, \
                         kg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * hg_295[k]
                   + f_0 * kg_490[k];

        t_386[k] = -2.0 * hg_296[k]
                   + f_0 * kg_491[k];

        t_387[k] = -2.0 * hg_297[k]
                   + f_0 * kg_492[k];

        t_388[k] = -2.0 * hg_298[k]
                   + f_0 * kg_493[k];

        t_389[k] = -2.0 * hg_299[k]
                   + f_0 * kg_494[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, hg_300, hg_301, hg_302, hg_303, \
                         hg_304, kg_495, kg_496, kg_497, kg_498, \
                         kg_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -hg_300[k]
                   + f_0 * kg_495[k];

        t_391[k] = -hg_301[k]
                   + f_0 * kg_496[k];

        t_392[k] = -hg_302[k]
                   + f_0 * kg_497[k];

        t_393[k] = -hg_303[k]
                   + f_0 * kg_498[k];

        t_394[k] = -hg_304[k]
                   + f_0 * kg_499[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, hg_305, hg_306, hg_307, hg_308, \
                         hg_309, kg_500, kg_501, kg_502, kg_503, \
                         kg_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -hg_305[k]
                   + f_0 * kg_500[k];

        t_396[k] = -hg_306[k]
                   + f_0 * kg_501[k];

        t_397[k] = -hg_307[k]
                   + f_0 * kg_502[k];

        t_398[k] = -hg_308[k]
                   + f_0 * kg_503[k];

        t_399[k] = -hg_309[k]
                   + f_0 * kg_504[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, hg_310, hg_311, hg_312, hg_313, \
                         hg_314, kg_505, kg_506, kg_507, kg_508, \
                         kg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -hg_310[k]
                   + f_0 * kg_505[k];

        t_401[k] = -hg_311[k]
                   + f_0 * kg_506[k];

        t_402[k] = -hg_312[k]
                   + f_0 * kg_507[k];

        t_403[k] = -hg_313[k]
                   + f_0 * kg_508[k];

        t_404[k] = -hg_314[k]
                   + f_0 * kg_509[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, kg_510, \
                         kg_511, kg_512, kg_513, kg_514, kg_515, kg_516, \
                         kg_517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * kg_510[k];

        t_406[k] = f_0 * kg_511[k];

        t_407[k] = f_0 * kg_512[k];

        t_408[k] = f_0 * kg_513[k];

        t_409[k] = f_0 * kg_514[k];

        t_410[k] = f_0 * kg_515[k];

        t_411[k] = f_0 * kg_516[k];

        t_412[k] = f_0 * kg_517[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, kg_518, kg_519, \
                         kg_520, kg_521, kg_522, kg_523, kg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * kg_518[k];

        t_414[k] = f_0 * kg_519[k];

        t_415[k] = f_0 * kg_520[k];

        t_416[k] = f_0 * kg_521[k];

        t_417[k] = f_0 * kg_522[k];

        t_418[k] = f_0 * kg_523[k];

        t_419[k] = f_0 * kg_524[k];
    }
}

auto
compute_prim_geom_10_ig_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t hg, const size_t kg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ig_electron_repulsion_1_piece0(buffer, target, hg, kg, ncols, alpha);

    compute_prim_geom_10_ig_electron_repulsion_1_piece1(buffer, target, hg, kg, ncols, alpha);

    compute_prim_geom_10_ig_electron_repulsion_1_piece2(buffer, target, hg, kg, ncols, alpha);
}

static auto
compute_prim_geom_10_ig_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hg, const size_t kg,
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

    const auto *kg_30 = buffer.data(kg + 30);
    const auto *kg_31 = buffer.data(kg + 31);
    const auto *kg_32 = buffer.data(kg + 32);
    const auto *kg_33 = buffer.data(kg + 33);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_35 = buffer.data(kg + 35);
    const auto *kg_36 = buffer.data(kg + 36);
    const auto *kg_37 = buffer.data(kg + 37);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_39 = buffer.data(kg + 39);
    const auto *kg_40 = buffer.data(kg + 40);
    const auto *kg_41 = buffer.data(kg + 41);
    const auto *kg_42 = buffer.data(kg + 42);
    const auto *kg_43 = buffer.data(kg + 43);
    const auto *kg_44 = buffer.data(kg + 44);
    const auto *kg_60 = buffer.data(kg + 60);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_62 = buffer.data(kg + 62);
    const auto *kg_63 = buffer.data(kg + 63);
    const auto *kg_64 = buffer.data(kg + 64);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_67 = buffer.data(kg + 67);
    const auto *kg_68 = buffer.data(kg + 68);
    const auto *kg_69 = buffer.data(kg + 69);
    const auto *kg_70 = buffer.data(kg + 70);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_72 = buffer.data(kg + 72);
    const auto *kg_73 = buffer.data(kg + 73);
    const auto *kg_74 = buffer.data(kg + 74);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_78 = buffer.data(kg + 78);
    const auto *kg_79 = buffer.data(kg + 79);
    const auto *kg_80 = buffer.data(kg + 80);
    const auto *kg_81 = buffer.data(kg + 81);
    const auto *kg_82 = buffer.data(kg + 82);
    const auto *kg_83 = buffer.data(kg + 83);
    const auto *kg_84 = buffer.data(kg + 84);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_86 = buffer.data(kg + 86);
    const auto *kg_87 = buffer.data(kg + 87);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_105 = buffer.data(kg + 105);
    const auto *kg_106 = buffer.data(kg + 106);
    const auto *kg_107 = buffer.data(kg + 107);
    const auto *kg_108 = buffer.data(kg + 108);
    const auto *kg_109 = buffer.data(kg + 109);
    const auto *kg_110 = buffer.data(kg + 110);
    const auto *kg_111 = buffer.data(kg + 111);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_114 = buffer.data(kg + 114);
    const auto *kg_115 = buffer.data(kg + 115);
    const auto *kg_116 = buffer.data(kg + 116);
    const auto *kg_117 = buffer.data(kg + 117);
    const auto *kg_118 = buffer.data(kg + 118);
    const auto *kg_119 = buffer.data(kg + 119);
    const auto *kg_120 = buffer.data(kg + 120);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_123 = buffer.data(kg + 123);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_125 = buffer.data(kg + 125);
    const auto *kg_126 = buffer.data(kg + 126);
    const auto *kg_127 = buffer.data(kg + 127);
    const auto *kg_128 = buffer.data(kg + 128);
    const auto *kg_129 = buffer.data(kg + 129);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_132 = buffer.data(kg + 132);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_134 = buffer.data(kg + 134);
    const auto *kg_135 = buffer.data(kg + 135);
    const auto *kg_136 = buffer.data(kg + 136);
    const auto *kg_137 = buffer.data(kg + 137);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_139 = buffer.data(kg + 139);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_141 = buffer.data(kg + 141);
    const auto *kg_142 = buffer.data(kg + 142);
    const auto *kg_143 = buffer.data(kg + 143);
    const auto *kg_144 = buffer.data(kg + 144);
    const auto *kg_145 = buffer.data(kg + 145);
    const auto *kg_146 = buffer.data(kg + 146);
    const auto *kg_147 = buffer.data(kg + 147);
    const auto *kg_148 = buffer.data(kg + 148);
    const auto *kg_149 = buffer.data(kg + 149);
    const auto *kg_165 = buffer.data(kg + 165);
    const auto *kg_166 = buffer.data(kg + 166);
    const auto *kg_167 = buffer.data(kg + 167);
    const auto *kg_168 = buffer.data(kg + 168);
    const auto *kg_169 = buffer.data(kg + 169);
    const auto *kg_170 = buffer.data(kg + 170);
    const auto *kg_171 = buffer.data(kg + 171);
    const auto *kg_172 = buffer.data(kg + 172);
    const auto *kg_173 = buffer.data(kg + 173);
    const auto *kg_174 = buffer.data(kg + 174);
    const auto *kg_175 = buffer.data(kg + 175);
    const auto *kg_176 = buffer.data(kg + 176);
    const auto *kg_177 = buffer.data(kg + 177);
    const auto *kg_178 = buffer.data(kg + 178);
    const auto *kg_179 = buffer.data(kg + 179);
    const auto *kg_180 = buffer.data(kg + 180);
    const auto *kg_181 = buffer.data(kg + 181);
    const auto *kg_182 = buffer.data(kg + 182);
    const auto *kg_183 = buffer.data(kg + 183);
    const auto *kg_184 = buffer.data(kg + 184);
    const auto *kg_185 = buffer.data(kg + 185);
    const auto *kg_186 = buffer.data(kg + 186);
    const auto *kg_187 = buffer.data(kg + 187);
    const auto *kg_188 = buffer.data(kg + 188);
    const auto *kg_189 = buffer.data(kg + 189);
    const auto *kg_190 = buffer.data(kg + 190);
    const auto *kg_191 = buffer.data(kg + 191);
    const auto *kg_192 = buffer.data(kg + 192);
    const auto *kg_193 = buffer.data(kg + 193);
    const auto *kg_194 = buffer.data(kg + 194);
    const auto *kg_195 = buffer.data(kg + 195);
    const auto *kg_196 = buffer.data(kg + 196);
    const auto *kg_197 = buffer.data(kg + 197);
    const auto *kg_198 = buffer.data(kg + 198);
    const auto *kg_199 = buffer.data(kg + 199);
    const auto *kg_200 = buffer.data(kg + 200);
    const auto *kg_201 = buffer.data(kg + 201);
    const auto *kg_202 = buffer.data(kg + 202);
    const auto *kg_203 = buffer.data(kg + 203);
    const auto *kg_204 = buffer.data(kg + 204);
    const auto *kg_205 = buffer.data(kg + 205);
    const auto *kg_206 = buffer.data(kg + 206);
    const auto *kg_207 = buffer.data(kg + 207);
    const auto *kg_208 = buffer.data(kg + 208);
    const auto *kg_209 = buffer.data(kg + 209);
    const auto *kg_210 = buffer.data(kg + 210);
    const auto *kg_211 = buffer.data(kg + 211);
    const auto *kg_212 = buffer.data(kg + 212);
    const auto *kg_213 = buffer.data(kg + 213);
    const auto *kg_214 = buffer.data(kg + 214);
    const auto *kg_215 = buffer.data(kg + 215);
    const auto *kg_216 = buffer.data(kg + 216);
    const auto *kg_217 = buffer.data(kg + 217);
    const auto *kg_218 = buffer.data(kg + 218);
    const auto *kg_219 = buffer.data(kg + 219);
    const auto *kg_220 = buffer.data(kg + 220);
    const auto *kg_221 = buffer.data(kg + 221);
    const auto *kg_222 = buffer.data(kg + 222);
    const auto *kg_223 = buffer.data(kg + 223);
    const auto *kg_224 = buffer.data(kg + 224);
    const auto *kg_240 = buffer.data(kg + 240);
    const auto *kg_241 = buffer.data(kg + 241);
    const auto *kg_242 = buffer.data(kg + 242);
    const auto *kg_243 = buffer.data(kg + 243);
    const auto *kg_244 = buffer.data(kg + 244);
    const auto *kg_245 = buffer.data(kg + 245);
    const auto *kg_246 = buffer.data(kg + 246);
    const auto *kg_247 = buffer.data(kg + 247);
    const auto *kg_248 = buffer.data(kg + 248);
    const auto *kg_249 = buffer.data(kg + 249);
    const auto *kg_250 = buffer.data(kg + 250);
    const auto *kg_251 = buffer.data(kg + 251);
    const auto *kg_252 = buffer.data(kg + 252);
    const auto *kg_253 = buffer.data(kg + 253);
    const auto *kg_254 = buffer.data(kg + 254);
    const auto *kg_255 = buffer.data(kg + 255);
    const auto *kg_256 = buffer.data(kg + 256);
    const auto *kg_257 = buffer.data(kg + 257);
    const auto *kg_258 = buffer.data(kg + 258);
    const auto *kg_259 = buffer.data(kg + 259);
    const auto *kg_260 = buffer.data(kg + 260);
    const auto *kg_261 = buffer.data(kg + 261);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, kg_30, kg_31, kg_32, kg_33, \
                         kg_34, kg_35, kg_36, kg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kg_30[k];

        t_1[k] = f_0 * kg_31[k];

        t_2[k] = f_0 * kg_32[k];

        t_3[k] = f_0 * kg_33[k];

        t_4[k] = f_0 * kg_34[k];

        t_5[k] = f_0 * kg_35[k];

        t_6[k] = f_0 * kg_36[k];

        t_7[k] = f_0 * kg_37[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, kg_38, kg_39, kg_40, \
                         kg_41, kg_42, kg_43, kg_44, kg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * kg_38[k];

        t_9[k] = f_0 * kg_39[k];

        t_10[k] = f_0 * kg_40[k];

        t_11[k] = f_0 * kg_41[k];

        t_12[k] = f_0 * kg_42[k];

        t_13[k] = f_0 * kg_43[k];

        t_14[k] = f_0 * kg_44[k];

        t_15[k] = f_0 * kg_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, kg_61, kg_62, kg_63, \
                         kg_64, kg_65, kg_66, kg_67, kg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * kg_61[k];

        t_17[k] = f_0 * kg_62[k];

        t_18[k] = f_0 * kg_63[k];

        t_19[k] = f_0 * kg_64[k];

        t_20[k] = f_0 * kg_65[k];

        t_21[k] = f_0 * kg_66[k];

        t_22[k] = f_0 * kg_67[k];

        t_23[k] = f_0 * kg_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, hg_0, kg_69, kg_70, kg_71, \
                         kg_72, kg_73, kg_74, kg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * kg_69[k];

        t_25[k] = f_0 * kg_70[k];

        t_26[k] = f_0 * kg_71[k];

        t_27[k] = f_0 * kg_72[k];

        t_28[k] = f_0 * kg_73[k];

        t_29[k] = f_0 * kg_74[k];

        t_30[k] = -hg_0[k]
                  + f_0 * kg_75[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, hg_1, hg_2, hg_3, hg_4, hg_5, kg_76, \
                         kg_77, kg_78, kg_79, kg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -hg_1[k]
                  + f_0 * kg_76[k];

        t_32[k] = -hg_2[k]
                  + f_0 * kg_77[k];

        t_33[k] = -hg_3[k]
                  + f_0 * kg_78[k];

        t_34[k] = -hg_4[k]
                  + f_0 * kg_79[k];

        t_35[k] = -hg_5[k]
                  + f_0 * kg_80[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, hg_6, hg_7, hg_8, hg_9, hg_10, kg_81, \
                         kg_82, kg_83, kg_84, kg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -hg_6[k]
                  + f_0 * kg_81[k];

        t_37[k] = -hg_7[k]
                  + f_0 * kg_82[k];

        t_38[k] = -hg_8[k]
                  + f_0 * kg_83[k];

        t_39[k] = -hg_9[k]
                  + f_0 * kg_84[k];

        t_40[k] = -hg_10[k]
                  + f_0 * kg_85[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, t_46, hg_11, hg_12, hg_13, hg_14, \
                         kg_86, kg_87, kg_88, kg_89, kg_105, kg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -hg_11[k]
                  + f_0 * kg_86[k];

        t_42[k] = -hg_12[k]
                  + f_0 * kg_87[k];

        t_43[k] = -hg_13[k]
                  + f_0 * kg_88[k];

        t_44[k] = -hg_14[k]
                  + f_0 * kg_89[k];

        t_45[k] = f_0 * kg_105[k];

        t_46[k] = f_0 * kg_106[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, t_53, t_54, kg_107, kg_108, \
                         kg_109, kg_110, kg_111, kg_112, kg_113, \
                         kg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * kg_107[k];

        t_48[k] = f_0 * kg_108[k];

        t_49[k] = f_0 * kg_109[k];

        t_50[k] = f_0 * kg_110[k];

        t_51[k] = f_0 * kg_111[k];

        t_52[k] = f_0 * kg_112[k];

        t_53[k] = f_0 * kg_113[k];

        t_54[k] = f_0 * kg_114[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, hg_15, hg_16, kg_115, \
                         kg_116, kg_117, kg_118, kg_119, kg_120, \
                         kg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * kg_115[k];

        t_56[k] = f_0 * kg_116[k];

        t_57[k] = f_0 * kg_117[k];

        t_58[k] = f_0 * kg_118[k];

        t_59[k] = f_0 * kg_119[k];

        t_60[k] = -hg_15[k]
                  + f_0 * kg_120[k];

        t_61[k] = -hg_16[k]
                  + f_0 * kg_121[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, hg_17, hg_18, hg_19, hg_20, hg_21, \
                         kg_122, kg_123, kg_124, kg_125, kg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -hg_17[k]
                  + f_0 * kg_122[k];

        t_63[k] = -hg_18[k]
                  + f_0 * kg_123[k];

        t_64[k] = -hg_19[k]
                  + f_0 * kg_124[k];

        t_65[k] = -hg_20[k]
                  + f_0 * kg_125[k];

        t_66[k] = -hg_21[k]
                  + f_0 * kg_126[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, hg_22, hg_23, hg_24, hg_25, hg_26, \
                         kg_127, kg_128, kg_129, kg_130, kg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -hg_22[k]
                  + f_0 * kg_127[k];

        t_68[k] = -hg_23[k]
                  + f_0 * kg_128[k];

        t_69[k] = -hg_24[k]
                  + f_0 * kg_129[k];

        t_70[k] = -hg_25[k]
                  + f_0 * kg_130[k];

        t_71[k] = -hg_26[k]
                  + f_0 * kg_131[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, hg_27, hg_28, hg_29, hg_30, hg_31, \
                         kg_132, kg_133, kg_134, kg_135, kg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -hg_27[k]
                  + f_0 * kg_132[k];

        t_73[k] = -hg_28[k]
                  + f_0 * kg_133[k];

        t_74[k] = -hg_29[k]
                  + f_0 * kg_134[k];

        t_75[k] = -2.0 * hg_30[k]
                  + f_0 * kg_135[k];

        t_76[k] = -2.0 * hg_31[k]
                  + f_0 * kg_136[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, hg_32, hg_33, hg_34, hg_35, hg_36, \
                         kg_137, kg_138, kg_139, kg_140, kg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -2.0 * hg_32[k]
                  + f_0 * kg_137[k];

        t_78[k] = -2.0 * hg_33[k]
                  + f_0 * kg_138[k];

        t_79[k] = -2.0 * hg_34[k]
                  + f_0 * kg_139[k];

        t_80[k] = -2.0 * hg_35[k]
                  + f_0 * kg_140[k];

        t_81[k] = -2.0 * hg_36[k]
                  + f_0 * kg_141[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, hg_37, hg_38, hg_39, hg_40, hg_41, \
                         kg_142, kg_143, kg_144, kg_145, kg_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -2.0 * hg_37[k]
                  + f_0 * kg_142[k];

        t_83[k] = -2.0 * hg_38[k]
                  + f_0 * kg_143[k];

        t_84[k] = -2.0 * hg_39[k]
                  + f_0 * kg_144[k];

        t_85[k] = -2.0 * hg_40[k]
                  + f_0 * kg_145[k];

        t_86[k] = -2.0 * hg_41[k]
                  + f_0 * kg_146[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, hg_42, hg_43, hg_44, kg_147, \
                         kg_148, kg_149, kg_165, kg_166, kg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -2.0 * hg_42[k]
                  + f_0 * kg_147[k];

        t_88[k] = -2.0 * hg_43[k]
                  + f_0 * kg_148[k];

        t_89[k] = -2.0 * hg_44[k]
                  + f_0 * kg_149[k];

        t_90[k] = f_0 * kg_165[k];

        t_91[k] = f_0 * kg_166[k];

        t_92[k] = f_0 * kg_167[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, kg_168, kg_169, \
                         kg_170, kg_171, kg_172, kg_173, kg_174, \
                         kg_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * kg_168[k];

        t_94[k] = f_0 * kg_169[k];

        t_95[k] = f_0 * kg_170[k];

        t_96[k] = f_0 * kg_171[k];

        t_97[k] = f_0 * kg_172[k];

        t_98[k] = f_0 * kg_173[k];

        t_99[k] = f_0 * kg_174[k];

        t_100[k] = f_0 * kg_175[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, hg_45, hg_46, kg_176, \
                         kg_177, kg_178, kg_179, kg_180, kg_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * kg_176[k];

        t_102[k] = f_0 * kg_177[k];

        t_103[k] = f_0 * kg_178[k];

        t_104[k] = f_0 * kg_179[k];

        t_105[k] = -hg_45[k]
                   + f_0 * kg_180[k];

        t_106[k] = -hg_46[k]
                   + f_0 * kg_181[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, hg_47, hg_48, hg_49, hg_50, hg_51, \
                         kg_182, kg_183, kg_184, kg_185, kg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -hg_47[k]
                   + f_0 * kg_182[k];

        t_108[k] = -hg_48[k]
                   + f_0 * kg_183[k];

        t_109[k] = -hg_49[k]
                   + f_0 * kg_184[k];

        t_110[k] = -hg_50[k]
                   + f_0 * kg_185[k];

        t_111[k] = -hg_51[k]
                   + f_0 * kg_186[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, hg_52, hg_53, hg_54, hg_55, hg_56, \
                         kg_187, kg_188, kg_189, kg_190, kg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -hg_52[k]
                   + f_0 * kg_187[k];

        t_113[k] = -hg_53[k]
                   + f_0 * kg_188[k];

        t_114[k] = -hg_54[k]
                   + f_0 * kg_189[k];

        t_115[k] = -hg_55[k]
                   + f_0 * kg_190[k];

        t_116[k] = -hg_56[k]
                   + f_0 * kg_191[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, hg_57, hg_58, hg_59, hg_60, hg_61, \
                         kg_192, kg_193, kg_194, kg_195, kg_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -hg_57[k]
                   + f_0 * kg_192[k];

        t_118[k] = -hg_58[k]
                   + f_0 * kg_193[k];

        t_119[k] = -hg_59[k]
                   + f_0 * kg_194[k];

        t_120[k] = -2.0 * hg_60[k]
                   + f_0 * kg_195[k];

        t_121[k] = -2.0 * hg_61[k]
                   + f_0 * kg_196[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, hg_62, hg_63, hg_64, hg_65, hg_66, \
                         kg_197, kg_198, kg_199, kg_200, kg_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * hg_62[k]
                   + f_0 * kg_197[k];

        t_123[k] = -2.0 * hg_63[k]
                   + f_0 * kg_198[k];

        t_124[k] = -2.0 * hg_64[k]
                   + f_0 * kg_199[k];

        t_125[k] = -2.0 * hg_65[k]
                   + f_0 * kg_200[k];

        t_126[k] = -2.0 * hg_66[k]
                   + f_0 * kg_201[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, hg_67, hg_68, hg_69, hg_70, hg_71, \
                         kg_202, kg_203, kg_204, kg_205, kg_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * hg_67[k]
                   + f_0 * kg_202[k];

        t_128[k] = -2.0 * hg_68[k]
                   + f_0 * kg_203[k];

        t_129[k] = -2.0 * hg_69[k]
                   + f_0 * kg_204[k];

        t_130[k] = -2.0 * hg_70[k]
                   + f_0 * kg_205[k];

        t_131[k] = -2.0 * hg_71[k]
                   + f_0 * kg_206[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, hg_72, hg_73, hg_74, hg_75, hg_76, \
                         kg_207, kg_208, kg_209, kg_210, kg_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * hg_72[k]
                   + f_0 * kg_207[k];

        t_133[k] = -2.0 * hg_73[k]
                   + f_0 * kg_208[k];

        t_134[k] = -2.0 * hg_74[k]
                   + f_0 * kg_209[k];

        t_135[k] = -3.0 * hg_75[k]
                   + f_0 * kg_210[k];

        t_136[k] = -3.0 * hg_76[k]
                   + f_0 * kg_211[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, hg_77, hg_78, hg_79, hg_80, hg_81, \
                         kg_212, kg_213, kg_214, kg_215, kg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -3.0 * hg_77[k]
                   + f_0 * kg_212[k];

        t_138[k] = -3.0 * hg_78[k]
                   + f_0 * kg_213[k];

        t_139[k] = -3.0 * hg_79[k]
                   + f_0 * kg_214[k];

        t_140[k] = -3.0 * hg_80[k]
                   + f_0 * kg_215[k];

        t_141[k] = -3.0 * hg_81[k]
                   + f_0 * kg_216[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, hg_82, hg_83, hg_84, hg_85, hg_86, \
                         kg_217, kg_218, kg_219, kg_220, kg_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -3.0 * hg_82[k]
                   + f_0 * kg_217[k];

        t_143[k] = -3.0 * hg_83[k]
                   + f_0 * kg_218[k];

        t_144[k] = -3.0 * hg_84[k]
                   + f_0 * kg_219[k];

        t_145[k] = -3.0 * hg_85[k]
                   + f_0 * kg_220[k];

        t_146[k] = -3.0 * hg_86[k]
                   + f_0 * kg_221[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, t_152, hg_87, hg_88, hg_89, \
                         kg_222, kg_223, kg_224, kg_240, kg_241, \
                         kg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -3.0 * hg_87[k]
                   + f_0 * kg_222[k];

        t_148[k] = -3.0 * hg_88[k]
                   + f_0 * kg_223[k];

        t_149[k] = -3.0 * hg_89[k]
                   + f_0 * kg_224[k];

        t_150[k] = f_0 * kg_240[k];

        t_151[k] = f_0 * kg_241[k];

        t_152[k] = f_0 * kg_242[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, t_158, t_159, t_160, kg_243, \
                         kg_244, kg_245, kg_246, kg_247, kg_248, kg_249, \
                         kg_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_0 * kg_243[k];

        t_154[k] = f_0 * kg_244[k];

        t_155[k] = f_0 * kg_245[k];

        t_156[k] = f_0 * kg_246[k];

        t_157[k] = f_0 * kg_247[k];

        t_158[k] = f_0 * kg_248[k];

        t_159[k] = f_0 * kg_249[k];

        t_160[k] = f_0 * kg_250[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, hg_90, hg_91, kg_251, \
                         kg_252, kg_253, kg_254, kg_255, kg_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_0 * kg_251[k];

        t_162[k] = f_0 * kg_252[k];

        t_163[k] = f_0 * kg_253[k];

        t_164[k] = f_0 * kg_254[k];

        t_165[k] = -hg_90[k]
                   + f_0 * kg_255[k];

        t_166[k] = -hg_91[k]
                   + f_0 * kg_256[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, hg_92, hg_93, hg_94, hg_95, hg_96, \
                         kg_257, kg_258, kg_259, kg_260, kg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -hg_92[k]
                   + f_0 * kg_257[k];

        t_168[k] = -hg_93[k]
                   + f_0 * kg_258[k];

        t_169[k] = -hg_94[k]
                   + f_0 * kg_259[k];

        t_170[k] = -hg_95[k]
                   + f_0 * kg_260[k];

        t_171[k] = -hg_96[k]
                   + f_0 * kg_261[k];
    }
}

static auto
compute_prim_geom_10_ig_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hg, const size_t kg,
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
    const auto *hg_225 = buffer.data(hg + 225);
    const auto *hg_226 = buffer.data(hg + 226);

    const auto *kg_262 = buffer.data(kg + 262);
    const auto *kg_263 = buffer.data(kg + 263);
    const auto *kg_264 = buffer.data(kg + 264);
    const auto *kg_265 = buffer.data(kg + 265);
    const auto *kg_266 = buffer.data(kg + 266);
    const auto *kg_267 = buffer.data(kg + 267);
    const auto *kg_268 = buffer.data(kg + 268);
    const auto *kg_269 = buffer.data(kg + 269);
    const auto *kg_270 = buffer.data(kg + 270);
    const auto *kg_271 = buffer.data(kg + 271);
    const auto *kg_272 = buffer.data(kg + 272);
    const auto *kg_273 = buffer.data(kg + 273);
    const auto *kg_274 = buffer.data(kg + 274);
    const auto *kg_275 = buffer.data(kg + 275);
    const auto *kg_276 = buffer.data(kg + 276);
    const auto *kg_277 = buffer.data(kg + 277);
    const auto *kg_278 = buffer.data(kg + 278);
    const auto *kg_279 = buffer.data(kg + 279);
    const auto *kg_280 = buffer.data(kg + 280);
    const auto *kg_281 = buffer.data(kg + 281);
    const auto *kg_282 = buffer.data(kg + 282);
    const auto *kg_283 = buffer.data(kg + 283);
    const auto *kg_284 = buffer.data(kg + 284);
    const auto *kg_285 = buffer.data(kg + 285);
    const auto *kg_286 = buffer.data(kg + 286);
    const auto *kg_287 = buffer.data(kg + 287);
    const auto *kg_288 = buffer.data(kg + 288);
    const auto *kg_289 = buffer.data(kg + 289);
    const auto *kg_290 = buffer.data(kg + 290);
    const auto *kg_291 = buffer.data(kg + 291);
    const auto *kg_292 = buffer.data(kg + 292);
    const auto *kg_293 = buffer.data(kg + 293);
    const auto *kg_294 = buffer.data(kg + 294);
    const auto *kg_295 = buffer.data(kg + 295);
    const auto *kg_296 = buffer.data(kg + 296);
    const auto *kg_297 = buffer.data(kg + 297);
    const auto *kg_298 = buffer.data(kg + 298);
    const auto *kg_299 = buffer.data(kg + 299);
    const auto *kg_300 = buffer.data(kg + 300);
    const auto *kg_301 = buffer.data(kg + 301);
    const auto *kg_302 = buffer.data(kg + 302);
    const auto *kg_303 = buffer.data(kg + 303);
    const auto *kg_304 = buffer.data(kg + 304);
    const auto *kg_305 = buffer.data(kg + 305);
    const auto *kg_306 = buffer.data(kg + 306);
    const auto *kg_307 = buffer.data(kg + 307);
    const auto *kg_308 = buffer.data(kg + 308);
    const auto *kg_309 = buffer.data(kg + 309);
    const auto *kg_310 = buffer.data(kg + 310);
    const auto *kg_311 = buffer.data(kg + 311);
    const auto *kg_312 = buffer.data(kg + 312);
    const auto *kg_313 = buffer.data(kg + 313);
    const auto *kg_314 = buffer.data(kg + 314);
    const auto *kg_330 = buffer.data(kg + 330);
    const auto *kg_331 = buffer.data(kg + 331);
    const auto *kg_332 = buffer.data(kg + 332);
    const auto *kg_333 = buffer.data(kg + 333);
    const auto *kg_334 = buffer.data(kg + 334);
    const auto *kg_335 = buffer.data(kg + 335);
    const auto *kg_336 = buffer.data(kg + 336);
    const auto *kg_337 = buffer.data(kg + 337);
    const auto *kg_338 = buffer.data(kg + 338);
    const auto *kg_339 = buffer.data(kg + 339);
    const auto *kg_340 = buffer.data(kg + 340);
    const auto *kg_341 = buffer.data(kg + 341);
    const auto *kg_342 = buffer.data(kg + 342);
    const auto *kg_343 = buffer.data(kg + 343);
    const auto *kg_344 = buffer.data(kg + 344);
    const auto *kg_345 = buffer.data(kg + 345);
    const auto *kg_346 = buffer.data(kg + 346);
    const auto *kg_347 = buffer.data(kg + 347);
    const auto *kg_348 = buffer.data(kg + 348);
    const auto *kg_349 = buffer.data(kg + 349);
    const auto *kg_350 = buffer.data(kg + 350);
    const auto *kg_351 = buffer.data(kg + 351);
    const auto *kg_352 = buffer.data(kg + 352);
    const auto *kg_353 = buffer.data(kg + 353);
    const auto *kg_354 = buffer.data(kg + 354);
    const auto *kg_355 = buffer.data(kg + 355);
    const auto *kg_356 = buffer.data(kg + 356);
    const auto *kg_357 = buffer.data(kg + 357);
    const auto *kg_358 = buffer.data(kg + 358);
    const auto *kg_359 = buffer.data(kg + 359);
    const auto *kg_360 = buffer.data(kg + 360);
    const auto *kg_361 = buffer.data(kg + 361);
    const auto *kg_362 = buffer.data(kg + 362);
    const auto *kg_363 = buffer.data(kg + 363);
    const auto *kg_364 = buffer.data(kg + 364);
    const auto *kg_365 = buffer.data(kg + 365);
    const auto *kg_366 = buffer.data(kg + 366);
    const auto *kg_367 = buffer.data(kg + 367);
    const auto *kg_368 = buffer.data(kg + 368);
    const auto *kg_369 = buffer.data(kg + 369);
    const auto *kg_370 = buffer.data(kg + 370);
    const auto *kg_371 = buffer.data(kg + 371);
    const auto *kg_372 = buffer.data(kg + 372);
    const auto *kg_373 = buffer.data(kg + 373);
    const auto *kg_374 = buffer.data(kg + 374);
    const auto *kg_375 = buffer.data(kg + 375);
    const auto *kg_376 = buffer.data(kg + 376);
    const auto *kg_377 = buffer.data(kg + 377);
    const auto *kg_378 = buffer.data(kg + 378);
    const auto *kg_379 = buffer.data(kg + 379);
    const auto *kg_380 = buffer.data(kg + 380);
    const auto *kg_381 = buffer.data(kg + 381);
    const auto *kg_382 = buffer.data(kg + 382);
    const auto *kg_383 = buffer.data(kg + 383);
    const auto *kg_384 = buffer.data(kg + 384);
    const auto *kg_385 = buffer.data(kg + 385);
    const auto *kg_386 = buffer.data(kg + 386);
    const auto *kg_387 = buffer.data(kg + 387);
    const auto *kg_388 = buffer.data(kg + 388);
    const auto *kg_389 = buffer.data(kg + 389);
    const auto *kg_390 = buffer.data(kg + 390);
    const auto *kg_391 = buffer.data(kg + 391);
    const auto *kg_392 = buffer.data(kg + 392);
    const auto *kg_393 = buffer.data(kg + 393);
    const auto *kg_394 = buffer.data(kg + 394);
    const auto *kg_395 = buffer.data(kg + 395);
    const auto *kg_396 = buffer.data(kg + 396);
    const auto *kg_397 = buffer.data(kg + 397);
    const auto *kg_398 = buffer.data(kg + 398);
    const auto *kg_399 = buffer.data(kg + 399);
    const auto *kg_400 = buffer.data(kg + 400);
    const auto *kg_401 = buffer.data(kg + 401);
    const auto *kg_402 = buffer.data(kg + 402);
    const auto *kg_403 = buffer.data(kg + 403);
    const auto *kg_404 = buffer.data(kg + 404);
    const auto *kg_405 = buffer.data(kg + 405);
    const auto *kg_406 = buffer.data(kg + 406);
    const auto *kg_407 = buffer.data(kg + 407);
    const auto *kg_408 = buffer.data(kg + 408);
    const auto *kg_409 = buffer.data(kg + 409);
    const auto *kg_410 = buffer.data(kg + 410);
    const auto *kg_411 = buffer.data(kg + 411);
    const auto *kg_412 = buffer.data(kg + 412);
    const auto *kg_413 = buffer.data(kg + 413);
    const auto *kg_414 = buffer.data(kg + 414);
    const auto *kg_415 = buffer.data(kg + 415);
    const auto *kg_416 = buffer.data(kg + 416);
    const auto *kg_417 = buffer.data(kg + 417);
    const auto *kg_418 = buffer.data(kg + 418);
    const auto *kg_419 = buffer.data(kg + 419);
    const auto *kg_435 = buffer.data(kg + 435);
    const auto *kg_436 = buffer.data(kg + 436);
    const auto *kg_437 = buffer.data(kg + 437);
    const auto *kg_438 = buffer.data(kg + 438);
    const auto *kg_439 = buffer.data(kg + 439);
    const auto *kg_440 = buffer.data(kg + 440);
    const auto *kg_441 = buffer.data(kg + 441);
    const auto *kg_442 = buffer.data(kg + 442);
    const auto *kg_443 = buffer.data(kg + 443);
    const auto *kg_444 = buffer.data(kg + 444);
    const auto *kg_445 = buffer.data(kg + 445);
    const auto *kg_446 = buffer.data(kg + 446);
    const auto *kg_447 = buffer.data(kg + 447);
    const auto *kg_448 = buffer.data(kg + 448);
    const auto *kg_449 = buffer.data(kg + 449);
    const auto *kg_450 = buffer.data(kg + 450);
    const auto *kg_451 = buffer.data(kg + 451);

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, hg_97, hg_98, hg_99, hg_100, \
                         hg_101, kg_262, kg_263, kg_264, kg_265, \
                         kg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -hg_97[k]
                   + f_0 * kg_262[k];

        t_173[k] = -hg_98[k]
                   + f_0 * kg_263[k];

        t_174[k] = -hg_99[k]
                   + f_0 * kg_264[k];

        t_175[k] = -hg_100[k]
                   + f_0 * kg_265[k];

        t_176[k] = -hg_101[k]
                   + f_0 * kg_266[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, hg_102, hg_103, hg_104, hg_105, \
                         hg_106, kg_267, kg_268, kg_269, kg_270, \
                         kg_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -hg_102[k]
                   + f_0 * kg_267[k];

        t_178[k] = -hg_103[k]
                   + f_0 * kg_268[k];

        t_179[k] = -hg_104[k]
                   + f_0 * kg_269[k];

        t_180[k] = -2.0 * hg_105[k]
                   + f_0 * kg_270[k];

        t_181[k] = -2.0 * hg_106[k]
                   + f_0 * kg_271[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, hg_107, hg_108, hg_109, hg_110, \
                         hg_111, kg_272, kg_273, kg_274, kg_275, \
                         kg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * hg_107[k]
                   + f_0 * kg_272[k];

        t_183[k] = -2.0 * hg_108[k]
                   + f_0 * kg_273[k];

        t_184[k] = -2.0 * hg_109[k]
                   + f_0 * kg_274[k];

        t_185[k] = -2.0 * hg_110[k]
                   + f_0 * kg_275[k];

        t_186[k] = -2.0 * hg_111[k]
                   + f_0 * kg_276[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, hg_112, hg_113, hg_114, hg_115, \
                         hg_116, kg_277, kg_278, kg_279, kg_280, \
                         kg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * hg_112[k]
                   + f_0 * kg_277[k];

        t_188[k] = -2.0 * hg_113[k]
                   + f_0 * kg_278[k];

        t_189[k] = -2.0 * hg_114[k]
                   + f_0 * kg_279[k];

        t_190[k] = -2.0 * hg_115[k]
                   + f_0 * kg_280[k];

        t_191[k] = -2.0 * hg_116[k]
                   + f_0 * kg_281[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, hg_117, hg_118, hg_119, hg_120, \
                         hg_121, kg_282, kg_283, kg_284, kg_285, \
                         kg_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -2.0 * hg_117[k]
                   + f_0 * kg_282[k];

        t_193[k] = -2.0 * hg_118[k]
                   + f_0 * kg_283[k];

        t_194[k] = -2.0 * hg_119[k]
                   + f_0 * kg_284[k];

        t_195[k] = -3.0 * hg_120[k]
                   + f_0 * kg_285[k];

        t_196[k] = -3.0 * hg_121[k]
                   + f_0 * kg_286[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, hg_122, hg_123, hg_124, hg_125, \
                         hg_126, kg_287, kg_288, kg_289, kg_290, \
                         kg_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -3.0 * hg_122[k]
                   + f_0 * kg_287[k];

        t_198[k] = -3.0 * hg_123[k]
                   + f_0 * kg_288[k];

        t_199[k] = -3.0 * hg_124[k]
                   + f_0 * kg_289[k];

        t_200[k] = -3.0 * hg_125[k]
                   + f_0 * kg_290[k];

        t_201[k] = -3.0 * hg_126[k]
                   + f_0 * kg_291[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, hg_127, hg_128, hg_129, hg_130, \
                         hg_131, kg_292, kg_293, kg_294, kg_295, \
                         kg_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -3.0 * hg_127[k]
                   + f_0 * kg_292[k];

        t_203[k] = -3.0 * hg_128[k]
                   + f_0 * kg_293[k];

        t_204[k] = -3.0 * hg_129[k]
                   + f_0 * kg_294[k];

        t_205[k] = -3.0 * hg_130[k]
                   + f_0 * kg_295[k];

        t_206[k] = -3.0 * hg_131[k]
                   + f_0 * kg_296[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, hg_132, hg_133, hg_134, hg_135, \
                         hg_136, kg_297, kg_298, kg_299, kg_300, \
                         kg_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -3.0 * hg_132[k]
                   + f_0 * kg_297[k];

        t_208[k] = -3.0 * hg_133[k]
                   + f_0 * kg_298[k];

        t_209[k] = -3.0 * hg_134[k]
                   + f_0 * kg_299[k];

        t_210[k] = -4.0 * hg_135[k]
                   + f_0 * kg_300[k];

        t_211[k] = -4.0 * hg_136[k]
                   + f_0 * kg_301[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, hg_137, hg_138, hg_139, hg_140, \
                         hg_141, kg_302, kg_303, kg_304, kg_305, \
                         kg_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = -4.0 * hg_137[k]
                   + f_0 * kg_302[k];

        t_213[k] = -4.0 * hg_138[k]
                   + f_0 * kg_303[k];

        t_214[k] = -4.0 * hg_139[k]
                   + f_0 * kg_304[k];

        t_215[k] = -4.0 * hg_140[k]
                   + f_0 * kg_305[k];

        t_216[k] = -4.0 * hg_141[k]
                   + f_0 * kg_306[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, hg_142, hg_143, hg_144, hg_145, \
                         hg_146, kg_307, kg_308, kg_309, kg_310, \
                         kg_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = -4.0 * hg_142[k]
                   + f_0 * kg_307[k];

        t_218[k] = -4.0 * hg_143[k]
                   + f_0 * kg_308[k];

        t_219[k] = -4.0 * hg_144[k]
                   + f_0 * kg_309[k];

        t_220[k] = -4.0 * hg_145[k]
                   + f_0 * kg_310[k];

        t_221[k] = -4.0 * hg_146[k]
                   + f_0 * kg_311[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, hg_147, hg_148, hg_149, \
                         kg_312, kg_313, kg_314, kg_330, kg_331, \
                         kg_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -4.0 * hg_147[k]
                   + f_0 * kg_312[k];

        t_223[k] = -4.0 * hg_148[k]
                   + f_0 * kg_313[k];

        t_224[k] = -4.0 * hg_149[k]
                   + f_0 * kg_314[k];

        t_225[k] = f_0 * kg_330[k];

        t_226[k] = f_0 * kg_331[k];

        t_227[k] = f_0 * kg_332[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, t_233, t_234, t_235, kg_333, \
                         kg_334, kg_335, kg_336, kg_337, kg_338, kg_339, \
                         kg_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_0 * kg_333[k];

        t_229[k] = f_0 * kg_334[k];

        t_230[k] = f_0 * kg_335[k];

        t_231[k] = f_0 * kg_336[k];

        t_232[k] = f_0 * kg_337[k];

        t_233[k] = f_0 * kg_338[k];

        t_234[k] = f_0 * kg_339[k];

        t_235[k] = f_0 * kg_340[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, t_241, hg_150, hg_151, kg_341, \
                         kg_342, kg_343, kg_344, kg_345, kg_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_0 * kg_341[k];

        t_237[k] = f_0 * kg_342[k];

        t_238[k] = f_0 * kg_343[k];

        t_239[k] = f_0 * kg_344[k];

        t_240[k] = -hg_150[k]
                   + f_0 * kg_345[k];

        t_241[k] = -hg_151[k]
                   + f_0 * kg_346[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, hg_152, hg_153, hg_154, hg_155, \
                         hg_156, kg_347, kg_348, kg_349, kg_350, \
                         kg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = -hg_152[k]
                   + f_0 * kg_347[k];

        t_243[k] = -hg_153[k]
                   + f_0 * kg_348[k];

        t_244[k] = -hg_154[k]
                   + f_0 * kg_349[k];

        t_245[k] = -hg_155[k]
                   + f_0 * kg_350[k];

        t_246[k] = -hg_156[k]
                   + f_0 * kg_351[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, hg_157, hg_158, hg_159, hg_160, \
                         hg_161, kg_352, kg_353, kg_354, kg_355, \
                         kg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = -hg_157[k]
                   + f_0 * kg_352[k];

        t_248[k] = -hg_158[k]
                   + f_0 * kg_353[k];

        t_249[k] = -hg_159[k]
                   + f_0 * kg_354[k];

        t_250[k] = -hg_160[k]
                   + f_0 * kg_355[k];

        t_251[k] = -hg_161[k]
                   + f_0 * kg_356[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, hg_162, hg_163, hg_164, hg_165, \
                         hg_166, kg_357, kg_358, kg_359, kg_360, \
                         kg_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -hg_162[k]
                   + f_0 * kg_357[k];

        t_253[k] = -hg_163[k]
                   + f_0 * kg_358[k];

        t_254[k] = -hg_164[k]
                   + f_0 * kg_359[k];

        t_255[k] = -2.0 * hg_165[k]
                   + f_0 * kg_360[k];

        t_256[k] = -2.0 * hg_166[k]
                   + f_0 * kg_361[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, hg_167, hg_168, hg_169, hg_170, \
                         hg_171, kg_362, kg_363, kg_364, kg_365, \
                         kg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = -2.0 * hg_167[k]
                   + f_0 * kg_362[k];

        t_258[k] = -2.0 * hg_168[k]
                   + f_0 * kg_363[k];

        t_259[k] = -2.0 * hg_169[k]
                   + f_0 * kg_364[k];

        t_260[k] = -2.0 * hg_170[k]
                   + f_0 * kg_365[k];

        t_261[k] = -2.0 * hg_171[k]
                   + f_0 * kg_366[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, hg_172, hg_173, hg_174, hg_175, \
                         hg_176, kg_367, kg_368, kg_369, kg_370, \
                         kg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = -2.0 * hg_172[k]
                   + f_0 * kg_367[k];

        t_263[k] = -2.0 * hg_173[k]
                   + f_0 * kg_368[k];

        t_264[k] = -2.0 * hg_174[k]
                   + f_0 * kg_369[k];

        t_265[k] = -2.0 * hg_175[k]
                   + f_0 * kg_370[k];

        t_266[k] = -2.0 * hg_176[k]
                   + f_0 * kg_371[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, hg_177, hg_178, hg_179, hg_180, \
                         hg_181, kg_372, kg_373, kg_374, kg_375, \
                         kg_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = -2.0 * hg_177[k]
                   + f_0 * kg_372[k];

        t_268[k] = -2.0 * hg_178[k]
                   + f_0 * kg_373[k];

        t_269[k] = -2.0 * hg_179[k]
                   + f_0 * kg_374[k];

        t_270[k] = -3.0 * hg_180[k]
                   + f_0 * kg_375[k];

        t_271[k] = -3.0 * hg_181[k]
                   + f_0 * kg_376[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, hg_182, hg_183, hg_184, hg_185, \
                         hg_186, kg_377, kg_378, kg_379, kg_380, \
                         kg_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -3.0 * hg_182[k]
                   + f_0 * kg_377[k];

        t_273[k] = -3.0 * hg_183[k]
                   + f_0 * kg_378[k];

        t_274[k] = -3.0 * hg_184[k]
                   + f_0 * kg_379[k];

        t_275[k] = -3.0 * hg_185[k]
                   + f_0 * kg_380[k];

        t_276[k] = -3.0 * hg_186[k]
                   + f_0 * kg_381[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, t_281, hg_187, hg_188, hg_189, hg_190, \
                         hg_191, kg_382, kg_383, kg_384, kg_385, \
                         kg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = -3.0 * hg_187[k]
                   + f_0 * kg_382[k];

        t_278[k] = -3.0 * hg_188[k]
                   + f_0 * kg_383[k];

        t_279[k] = -3.0 * hg_189[k]
                   + f_0 * kg_384[k];

        t_280[k] = -3.0 * hg_190[k]
                   + f_0 * kg_385[k];

        t_281[k] = -3.0 * hg_191[k]
                   + f_0 * kg_386[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, hg_192, hg_193, hg_194, hg_195, \
                         hg_196, kg_387, kg_388, kg_389, kg_390, \
                         kg_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = -3.0 * hg_192[k]
                   + f_0 * kg_387[k];

        t_283[k] = -3.0 * hg_193[k]
                   + f_0 * kg_388[k];

        t_284[k] = -3.0 * hg_194[k]
                   + f_0 * kg_389[k];

        t_285[k] = -4.0 * hg_195[k]
                   + f_0 * kg_390[k];

        t_286[k] = -4.0 * hg_196[k]
                   + f_0 * kg_391[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, hg_197, hg_198, hg_199, hg_200, \
                         hg_201, kg_392, kg_393, kg_394, kg_395, \
                         kg_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = -4.0 * hg_197[k]
                   + f_0 * kg_392[k];

        t_288[k] = -4.0 * hg_198[k]
                   + f_0 * kg_393[k];

        t_289[k] = -4.0 * hg_199[k]
                   + f_0 * kg_394[k];

        t_290[k] = -4.0 * hg_200[k]
                   + f_0 * kg_395[k];

        t_291[k] = -4.0 * hg_201[k]
                   + f_0 * kg_396[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, hg_202, hg_203, hg_204, hg_205, \
                         hg_206, kg_397, kg_398, kg_399, kg_400, \
                         kg_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = -4.0 * hg_202[k]
                   + f_0 * kg_397[k];

        t_293[k] = -4.0 * hg_203[k]
                   + f_0 * kg_398[k];

        t_294[k] = -4.0 * hg_204[k]
                   + f_0 * kg_399[k];

        t_295[k] = -4.0 * hg_205[k]
                   + f_0 * kg_400[k];

        t_296[k] = -4.0 * hg_206[k]
                   + f_0 * kg_401[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, hg_207, hg_208, hg_209, hg_210, \
                         hg_211, kg_402, kg_403, kg_404, kg_405, \
                         kg_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = -4.0 * hg_207[k]
                   + f_0 * kg_402[k];

        t_298[k] = -4.0 * hg_208[k]
                   + f_0 * kg_403[k];

        t_299[k] = -4.0 * hg_209[k]
                   + f_0 * kg_404[k];

        t_300[k] = -5.0 * hg_210[k]
                   + f_0 * kg_405[k];

        t_301[k] = -5.0 * hg_211[k]
                   + f_0 * kg_406[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, hg_212, hg_213, hg_214, hg_215, \
                         hg_216, kg_407, kg_408, kg_409, kg_410, \
                         kg_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = -5.0 * hg_212[k]
                   + f_0 * kg_407[k];

        t_303[k] = -5.0 * hg_213[k]
                   + f_0 * kg_408[k];

        t_304[k] = -5.0 * hg_214[k]
                   + f_0 * kg_409[k];

        t_305[k] = -5.0 * hg_215[k]
                   + f_0 * kg_410[k];

        t_306[k] = -5.0 * hg_216[k]
                   + f_0 * kg_411[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, hg_217, hg_218, hg_219, hg_220, \
                         hg_221, kg_412, kg_413, kg_414, kg_415, \
                         kg_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = -5.0 * hg_217[k]
                   + f_0 * kg_412[k];

        t_308[k] = -5.0 * hg_218[k]
                   + f_0 * kg_413[k];

        t_309[k] = -5.0 * hg_219[k]
                   + f_0 * kg_414[k];

        t_310[k] = -5.0 * hg_220[k]
                   + f_0 * kg_415[k];

        t_311[k] = -5.0 * hg_221[k]
                   + f_0 * kg_416[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, t_317, hg_222, hg_223, hg_224, \
                         kg_417, kg_418, kg_419, kg_435, kg_436, \
                         kg_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -5.0 * hg_222[k]
                   + f_0 * kg_417[k];

        t_313[k] = -5.0 * hg_223[k]
                   + f_0 * kg_418[k];

        t_314[k] = -5.0 * hg_224[k]
                   + f_0 * kg_419[k];

        t_315[k] = f_0 * kg_435[k];

        t_316[k] = f_0 * kg_436[k];

        t_317[k] = f_0 * kg_437[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, t_323, t_324, t_325, kg_438, \
                         kg_439, kg_440, kg_441, kg_442, kg_443, kg_444, \
                         kg_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_0 * kg_438[k];

        t_319[k] = f_0 * kg_439[k];

        t_320[k] = f_0 * kg_440[k];

        t_321[k] = f_0 * kg_441[k];

        t_322[k] = f_0 * kg_442[k];

        t_323[k] = f_0 * kg_443[k];

        t_324[k] = f_0 * kg_444[k];

        t_325[k] = f_0 * kg_445[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, hg_225, hg_226, kg_446, \
                         kg_447, kg_448, kg_449, kg_450, kg_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_0 * kg_446[k];

        t_327[k] = f_0 * kg_447[k];

        t_328[k] = f_0 * kg_448[k];

        t_329[k] = f_0 * kg_449[k];

        t_330[k] = -hg_225[k]
                   + f_0 * kg_450[k];

        t_331[k] = -hg_226[k]
                   + f_0 * kg_451[k];
    }
}

static auto
compute_prim_geom_10_ig_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hg, const size_t kg,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *kg_452 = buffer.data(kg + 452);
    const auto *kg_453 = buffer.data(kg + 453);
    const auto *kg_454 = buffer.data(kg + 454);
    const auto *kg_455 = buffer.data(kg + 455);
    const auto *kg_456 = buffer.data(kg + 456);
    const auto *kg_457 = buffer.data(kg + 457);
    const auto *kg_458 = buffer.data(kg + 458);
    const auto *kg_459 = buffer.data(kg + 459);
    const auto *kg_460 = buffer.data(kg + 460);
    const auto *kg_461 = buffer.data(kg + 461);
    const auto *kg_462 = buffer.data(kg + 462);
    const auto *kg_463 = buffer.data(kg + 463);
    const auto *kg_464 = buffer.data(kg + 464);
    const auto *kg_465 = buffer.data(kg + 465);
    const auto *kg_466 = buffer.data(kg + 466);
    const auto *kg_467 = buffer.data(kg + 467);
    const auto *kg_468 = buffer.data(kg + 468);
    const auto *kg_469 = buffer.data(kg + 469);
    const auto *kg_470 = buffer.data(kg + 470);
    const auto *kg_471 = buffer.data(kg + 471);
    const auto *kg_472 = buffer.data(kg + 472);
    const auto *kg_473 = buffer.data(kg + 473);
    const auto *kg_474 = buffer.data(kg + 474);
    const auto *kg_475 = buffer.data(kg + 475);
    const auto *kg_476 = buffer.data(kg + 476);
    const auto *kg_477 = buffer.data(kg + 477);
    const auto *kg_478 = buffer.data(kg + 478);
    const auto *kg_479 = buffer.data(kg + 479);
    const auto *kg_480 = buffer.data(kg + 480);
    const auto *kg_481 = buffer.data(kg + 481);
    const auto *kg_482 = buffer.data(kg + 482);
    const auto *kg_483 = buffer.data(kg + 483);
    const auto *kg_484 = buffer.data(kg + 484);
    const auto *kg_485 = buffer.data(kg + 485);
    const auto *kg_486 = buffer.data(kg + 486);
    const auto *kg_487 = buffer.data(kg + 487);
    const auto *kg_488 = buffer.data(kg + 488);
    const auto *kg_489 = buffer.data(kg + 489);
    const auto *kg_490 = buffer.data(kg + 490);
    const auto *kg_491 = buffer.data(kg + 491);
    const auto *kg_492 = buffer.data(kg + 492);
    const auto *kg_493 = buffer.data(kg + 493);
    const auto *kg_494 = buffer.data(kg + 494);
    const auto *kg_495 = buffer.data(kg + 495);
    const auto *kg_496 = buffer.data(kg + 496);
    const auto *kg_497 = buffer.data(kg + 497);
    const auto *kg_498 = buffer.data(kg + 498);
    const auto *kg_499 = buffer.data(kg + 499);
    const auto *kg_500 = buffer.data(kg + 500);
    const auto *kg_501 = buffer.data(kg + 501);
    const auto *kg_502 = buffer.data(kg + 502);
    const auto *kg_503 = buffer.data(kg + 503);
    const auto *kg_504 = buffer.data(kg + 504);
    const auto *kg_505 = buffer.data(kg + 505);
    const auto *kg_506 = buffer.data(kg + 506);
    const auto *kg_507 = buffer.data(kg + 507);
    const auto *kg_508 = buffer.data(kg + 508);
    const auto *kg_509 = buffer.data(kg + 509);
    const auto *kg_510 = buffer.data(kg + 510);
    const auto *kg_511 = buffer.data(kg + 511);
    const auto *kg_512 = buffer.data(kg + 512);
    const auto *kg_513 = buffer.data(kg + 513);
    const auto *kg_514 = buffer.data(kg + 514);
    const auto *kg_515 = buffer.data(kg + 515);
    const auto *kg_516 = buffer.data(kg + 516);
    const auto *kg_517 = buffer.data(kg + 517);
    const auto *kg_518 = buffer.data(kg + 518);
    const auto *kg_519 = buffer.data(kg + 519);
    const auto *kg_520 = buffer.data(kg + 520);
    const auto *kg_521 = buffer.data(kg + 521);
    const auto *kg_522 = buffer.data(kg + 522);
    const auto *kg_523 = buffer.data(kg + 523);
    const auto *kg_524 = buffer.data(kg + 524);
    const auto *kg_525 = buffer.data(kg + 525);
    const auto *kg_526 = buffer.data(kg + 526);
    const auto *kg_527 = buffer.data(kg + 527);
    const auto *kg_528 = buffer.data(kg + 528);
    const auto *kg_529 = buffer.data(kg + 529);
    const auto *kg_530 = buffer.data(kg + 530);
    const auto *kg_531 = buffer.data(kg + 531);
    const auto *kg_532 = buffer.data(kg + 532);
    const auto *kg_533 = buffer.data(kg + 533);
    const auto *kg_534 = buffer.data(kg + 534);
    const auto *kg_535 = buffer.data(kg + 535);
    const auto *kg_536 = buffer.data(kg + 536);
    const auto *kg_537 = buffer.data(kg + 537);
    const auto *kg_538 = buffer.data(kg + 538);
    const auto *kg_539 = buffer.data(kg + 539);

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, hg_227, hg_228, hg_229, hg_230, \
                         hg_231, kg_452, kg_453, kg_454, kg_455, \
                         kg_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -hg_227[k]
                   + f_0 * kg_452[k];

        t_333[k] = -hg_228[k]
                   + f_0 * kg_453[k];

        t_334[k] = -hg_229[k]
                   + f_0 * kg_454[k];

        t_335[k] = -hg_230[k]
                   + f_0 * kg_455[k];

        t_336[k] = -hg_231[k]
                   + f_0 * kg_456[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, hg_232, hg_233, hg_234, hg_235, \
                         hg_236, kg_457, kg_458, kg_459, kg_460, \
                         kg_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -hg_232[k]
                   + f_0 * kg_457[k];

        t_338[k] = -hg_233[k]
                   + f_0 * kg_458[k];

        t_339[k] = -hg_234[k]
                   + f_0 * kg_459[k];

        t_340[k] = -hg_235[k]
                   + f_0 * kg_460[k];

        t_341[k] = -hg_236[k]
                   + f_0 * kg_461[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, hg_237, hg_238, hg_239, hg_240, \
                         hg_241, kg_462, kg_463, kg_464, kg_465, \
                         kg_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -hg_237[k]
                   + f_0 * kg_462[k];

        t_343[k] = -hg_238[k]
                   + f_0 * kg_463[k];

        t_344[k] = -hg_239[k]
                   + f_0 * kg_464[k];

        t_345[k] = -2.0 * hg_240[k]
                   + f_0 * kg_465[k];

        t_346[k] = -2.0 * hg_241[k]
                   + f_0 * kg_466[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, hg_242, hg_243, hg_244, hg_245, \
                         hg_246, kg_467, kg_468, kg_469, kg_470, \
                         kg_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -2.0 * hg_242[k]
                   + f_0 * kg_467[k];

        t_348[k] = -2.0 * hg_243[k]
                   + f_0 * kg_468[k];

        t_349[k] = -2.0 * hg_244[k]
                   + f_0 * kg_469[k];

        t_350[k] = -2.0 * hg_245[k]
                   + f_0 * kg_470[k];

        t_351[k] = -2.0 * hg_246[k]
                   + f_0 * kg_471[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, hg_247, hg_248, hg_249, hg_250, \
                         hg_251, kg_472, kg_473, kg_474, kg_475, \
                         kg_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -2.0 * hg_247[k]
                   + f_0 * kg_472[k];

        t_353[k] = -2.0 * hg_248[k]
                   + f_0 * kg_473[k];

        t_354[k] = -2.0 * hg_249[k]
                   + f_0 * kg_474[k];

        t_355[k] = -2.0 * hg_250[k]
                   + f_0 * kg_475[k];

        t_356[k] = -2.0 * hg_251[k]
                   + f_0 * kg_476[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, hg_252, hg_253, hg_254, hg_255, \
                         hg_256, kg_477, kg_478, kg_479, kg_480, \
                         kg_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -2.0 * hg_252[k]
                   + f_0 * kg_477[k];

        t_358[k] = -2.0 * hg_253[k]
                   + f_0 * kg_478[k];

        t_359[k] = -2.0 * hg_254[k]
                   + f_0 * kg_479[k];

        t_360[k] = -3.0 * hg_255[k]
                   + f_0 * kg_480[k];

        t_361[k] = -3.0 * hg_256[k]
                   + f_0 * kg_481[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, hg_257, hg_258, hg_259, hg_260, \
                         hg_261, kg_482, kg_483, kg_484, kg_485, \
                         kg_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -3.0 * hg_257[k]
                   + f_0 * kg_482[k];

        t_363[k] = -3.0 * hg_258[k]
                   + f_0 * kg_483[k];

        t_364[k] = -3.0 * hg_259[k]
                   + f_0 * kg_484[k];

        t_365[k] = -3.0 * hg_260[k]
                   + f_0 * kg_485[k];

        t_366[k] = -3.0 * hg_261[k]
                   + f_0 * kg_486[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, hg_262, hg_263, hg_264, hg_265, \
                         hg_266, kg_487, kg_488, kg_489, kg_490, \
                         kg_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -3.0 * hg_262[k]
                   + f_0 * kg_487[k];

        t_368[k] = -3.0 * hg_263[k]
                   + f_0 * kg_488[k];

        t_369[k] = -3.0 * hg_264[k]
                   + f_0 * kg_489[k];

        t_370[k] = -3.0 * hg_265[k]
                   + f_0 * kg_490[k];

        t_371[k] = -3.0 * hg_266[k]
                   + f_0 * kg_491[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, hg_267, hg_268, hg_269, hg_270, \
                         hg_271, kg_492, kg_493, kg_494, kg_495, \
                         kg_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -3.0 * hg_267[k]
                   + f_0 * kg_492[k];

        t_373[k] = -3.0 * hg_268[k]
                   + f_0 * kg_493[k];

        t_374[k] = -3.0 * hg_269[k]
                   + f_0 * kg_494[k];

        t_375[k] = -4.0 * hg_270[k]
                   + f_0 * kg_495[k];

        t_376[k] = -4.0 * hg_271[k]
                   + f_0 * kg_496[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, hg_272, hg_273, hg_274, hg_275, \
                         hg_276, kg_497, kg_498, kg_499, kg_500, \
                         kg_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -4.0 * hg_272[k]
                   + f_0 * kg_497[k];

        t_378[k] = -4.0 * hg_273[k]
                   + f_0 * kg_498[k];

        t_379[k] = -4.0 * hg_274[k]
                   + f_0 * kg_499[k];

        t_380[k] = -4.0 * hg_275[k]
                   + f_0 * kg_500[k];

        t_381[k] = -4.0 * hg_276[k]
                   + f_0 * kg_501[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, hg_277, hg_278, hg_279, hg_280, \
                         hg_281, kg_502, kg_503, kg_504, kg_505, \
                         kg_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -4.0 * hg_277[k]
                   + f_0 * kg_502[k];

        t_383[k] = -4.0 * hg_278[k]
                   + f_0 * kg_503[k];

        t_384[k] = -4.0 * hg_279[k]
                   + f_0 * kg_504[k];

        t_385[k] = -4.0 * hg_280[k]
                   + f_0 * kg_505[k];

        t_386[k] = -4.0 * hg_281[k]
                   + f_0 * kg_506[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, hg_282, hg_283, hg_284, hg_285, \
                         hg_286, kg_507, kg_508, kg_509, kg_510, \
                         kg_511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -4.0 * hg_282[k]
                   + f_0 * kg_507[k];

        t_388[k] = -4.0 * hg_283[k]
                   + f_0 * kg_508[k];

        t_389[k] = -4.0 * hg_284[k]
                   + f_0 * kg_509[k];

        t_390[k] = -5.0 * hg_285[k]
                   + f_0 * kg_510[k];

        t_391[k] = -5.0 * hg_286[k]
                   + f_0 * kg_511[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, hg_287, hg_288, hg_289, hg_290, \
                         hg_291, kg_512, kg_513, kg_514, kg_515, \
                         kg_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -5.0 * hg_287[k]
                   + f_0 * kg_512[k];

        t_393[k] = -5.0 * hg_288[k]
                   + f_0 * kg_513[k];

        t_394[k] = -5.0 * hg_289[k]
                   + f_0 * kg_514[k];

        t_395[k] = -5.0 * hg_290[k]
                   + f_0 * kg_515[k];

        t_396[k] = -5.0 * hg_291[k]
                   + f_0 * kg_516[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, hg_292, hg_293, hg_294, hg_295, \
                         hg_296, kg_517, kg_518, kg_519, kg_520, \
                         kg_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -5.0 * hg_292[k]
                   + f_0 * kg_517[k];

        t_398[k] = -5.0 * hg_293[k]
                   + f_0 * kg_518[k];

        t_399[k] = -5.0 * hg_294[k]
                   + f_0 * kg_519[k];

        t_400[k] = -5.0 * hg_295[k]
                   + f_0 * kg_520[k];

        t_401[k] = -5.0 * hg_296[k]
                   + f_0 * kg_521[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, hg_297, hg_298, hg_299, hg_300, \
                         hg_301, kg_522, kg_523, kg_524, kg_525, \
                         kg_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -5.0 * hg_297[k]
                   + f_0 * kg_522[k];

        t_403[k] = -5.0 * hg_298[k]
                   + f_0 * kg_523[k];

        t_404[k] = -5.0 * hg_299[k]
                   + f_0 * kg_524[k];

        t_405[k] = -6.0 * hg_300[k]
                   + f_0 * kg_525[k];

        t_406[k] = -6.0 * hg_301[k]
                   + f_0 * kg_526[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, hg_302, hg_303, hg_304, hg_305, \
                         hg_306, kg_527, kg_528, kg_529, kg_530, \
                         kg_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -6.0 * hg_302[k]
                   + f_0 * kg_527[k];

        t_408[k] = -6.0 * hg_303[k]
                   + f_0 * kg_528[k];

        t_409[k] = -6.0 * hg_304[k]
                   + f_0 * kg_529[k];

        t_410[k] = -6.0 * hg_305[k]
                   + f_0 * kg_530[k];

        t_411[k] = -6.0 * hg_306[k]
                   + f_0 * kg_531[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, hg_307, hg_308, hg_309, hg_310, \
                         hg_311, kg_532, kg_533, kg_534, kg_535, \
                         kg_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -6.0 * hg_307[k]
                   + f_0 * kg_532[k];

        t_413[k] = -6.0 * hg_308[k]
                   + f_0 * kg_533[k];

        t_414[k] = -6.0 * hg_309[k]
                   + f_0 * kg_534[k];

        t_415[k] = -6.0 * hg_310[k]
                   + f_0 * kg_535[k];

        t_416[k] = -6.0 * hg_311[k]
                   + f_0 * kg_536[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, hg_312, hg_313, hg_314, kg_537, kg_538, \
                         kg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -6.0 * hg_312[k]
                   + f_0 * kg_537[k];

        t_418[k] = -6.0 * hg_313[k]
                   + f_0 * kg_538[k];

        t_419[k] = -6.0 * hg_314[k]
                   + f_0 * kg_539[k];
    }
}

auto
compute_prim_geom_10_ig_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t hg, const size_t kg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ig_electron_repulsion_2_piece0(buffer, target, hg, kg, ncols, alpha);

    compute_prim_geom_10_ig_electron_repulsion_2_piece1(buffer, target, hg, kg, ncols, alpha);

    compute_prim_geom_10_ig_electron_repulsion_2_piece2(buffer, target, hg, kg, ncols, alpha);
}

}  // namespace simdt2ceri
