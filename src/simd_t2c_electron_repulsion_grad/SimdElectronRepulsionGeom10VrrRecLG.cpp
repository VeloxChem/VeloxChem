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


#include "SimdElectronRepulsionGeom10VrrRecLG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_lg_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
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

    const auto *mg_0 = buffer.data(mg + 0);
    const auto *mg_1 = buffer.data(mg + 1);
    const auto *mg_2 = buffer.data(mg + 2);
    const auto *mg_3 = buffer.data(mg + 3);
    const auto *mg_4 = buffer.data(mg + 4);
    const auto *mg_5 = buffer.data(mg + 5);
    const auto *mg_6 = buffer.data(mg + 6);
    const auto *mg_7 = buffer.data(mg + 7);
    const auto *mg_8 = buffer.data(mg + 8);
    const auto *mg_9 = buffer.data(mg + 9);
    const auto *mg_10 = buffer.data(mg + 10);
    const auto *mg_11 = buffer.data(mg + 11);
    const auto *mg_12 = buffer.data(mg + 12);
    const auto *mg_13 = buffer.data(mg + 13);
    const auto *mg_14 = buffer.data(mg + 14);
    const auto *mg_15 = buffer.data(mg + 15);
    const auto *mg_16 = buffer.data(mg + 16);
    const auto *mg_17 = buffer.data(mg + 17);
    const auto *mg_18 = buffer.data(mg + 18);
    const auto *mg_19 = buffer.data(mg + 19);
    const auto *mg_20 = buffer.data(mg + 20);
    const auto *mg_21 = buffer.data(mg + 21);
    const auto *mg_22 = buffer.data(mg + 22);
    const auto *mg_23 = buffer.data(mg + 23);
    const auto *mg_24 = buffer.data(mg + 24);
    const auto *mg_25 = buffer.data(mg + 25);
    const auto *mg_26 = buffer.data(mg + 26);
    const auto *mg_27 = buffer.data(mg + 27);
    const auto *mg_28 = buffer.data(mg + 28);
    const auto *mg_29 = buffer.data(mg + 29);
    const auto *mg_30 = buffer.data(mg + 30);
    const auto *mg_31 = buffer.data(mg + 31);
    const auto *mg_32 = buffer.data(mg + 32);
    const auto *mg_33 = buffer.data(mg + 33);
    const auto *mg_34 = buffer.data(mg + 34);
    const auto *mg_35 = buffer.data(mg + 35);
    const auto *mg_36 = buffer.data(mg + 36);
    const auto *mg_37 = buffer.data(mg + 37);
    const auto *mg_38 = buffer.data(mg + 38);
    const auto *mg_39 = buffer.data(mg + 39);
    const auto *mg_40 = buffer.data(mg + 40);
    const auto *mg_41 = buffer.data(mg + 41);
    const auto *mg_42 = buffer.data(mg + 42);
    const auto *mg_43 = buffer.data(mg + 43);
    const auto *mg_44 = buffer.data(mg + 44);
    const auto *mg_45 = buffer.data(mg + 45);
    const auto *mg_46 = buffer.data(mg + 46);
    const auto *mg_47 = buffer.data(mg + 47);
    const auto *mg_48 = buffer.data(mg + 48);
    const auto *mg_49 = buffer.data(mg + 49);
    const auto *mg_50 = buffer.data(mg + 50);
    const auto *mg_51 = buffer.data(mg + 51);
    const auto *mg_52 = buffer.data(mg + 52);
    const auto *mg_53 = buffer.data(mg + 53);
    const auto *mg_54 = buffer.data(mg + 54);
    const auto *mg_55 = buffer.data(mg + 55);
    const auto *mg_56 = buffer.data(mg + 56);
    const auto *mg_57 = buffer.data(mg + 57);
    const auto *mg_58 = buffer.data(mg + 58);
    const auto *mg_59 = buffer.data(mg + 59);
    const auto *mg_60 = buffer.data(mg + 60);
    const auto *mg_61 = buffer.data(mg + 61);
    const auto *mg_62 = buffer.data(mg + 62);
    const auto *mg_63 = buffer.data(mg + 63);
    const auto *mg_64 = buffer.data(mg + 64);
    const auto *mg_65 = buffer.data(mg + 65);
    const auto *mg_66 = buffer.data(mg + 66);
    const auto *mg_67 = buffer.data(mg + 67);
    const auto *mg_68 = buffer.data(mg + 68);
    const auto *mg_69 = buffer.data(mg + 69);
    const auto *mg_70 = buffer.data(mg + 70);
    const auto *mg_71 = buffer.data(mg + 71);
    const auto *mg_72 = buffer.data(mg + 72);
    const auto *mg_73 = buffer.data(mg + 73);
    const auto *mg_74 = buffer.data(mg + 74);
    const auto *mg_75 = buffer.data(mg + 75);
    const auto *mg_76 = buffer.data(mg + 76);
    const auto *mg_77 = buffer.data(mg + 77);
    const auto *mg_78 = buffer.data(mg + 78);
    const auto *mg_79 = buffer.data(mg + 79);
    const auto *mg_80 = buffer.data(mg + 80);
    const auto *mg_81 = buffer.data(mg + 81);
    const auto *mg_82 = buffer.data(mg + 82);
    const auto *mg_83 = buffer.data(mg + 83);
    const auto *mg_84 = buffer.data(mg + 84);
    const auto *mg_85 = buffer.data(mg + 85);
    const auto *mg_86 = buffer.data(mg + 86);
    const auto *mg_87 = buffer.data(mg + 87);
    const auto *mg_88 = buffer.data(mg + 88);
    const auto *mg_89 = buffer.data(mg + 89);
    const auto *mg_90 = buffer.data(mg + 90);
    const auto *mg_91 = buffer.data(mg + 91);
    const auto *mg_92 = buffer.data(mg + 92);
    const auto *mg_93 = buffer.data(mg + 93);
    const auto *mg_94 = buffer.data(mg + 94);
    const auto *mg_95 = buffer.data(mg + 95);
    const auto *mg_96 = buffer.data(mg + 96);
    const auto *mg_97 = buffer.data(mg + 97);
    const auto *mg_98 = buffer.data(mg + 98);
    const auto *mg_99 = buffer.data(mg + 99);
    const auto *mg_100 = buffer.data(mg + 100);
    const auto *mg_101 = buffer.data(mg + 101);
    const auto *mg_102 = buffer.data(mg + 102);
    const auto *mg_103 = buffer.data(mg + 103);
    const auto *mg_104 = buffer.data(mg + 104);
    const auto *mg_105 = buffer.data(mg + 105);
    const auto *mg_106 = buffer.data(mg + 106);
    const auto *mg_107 = buffer.data(mg + 107);
    const auto *mg_108 = buffer.data(mg + 108);
    const auto *mg_109 = buffer.data(mg + 109);
    const auto *mg_110 = buffer.data(mg + 110);
    const auto *mg_111 = buffer.data(mg + 111);
    const auto *mg_112 = buffer.data(mg + 112);
    const auto *mg_113 = buffer.data(mg + 113);
    const auto *mg_114 = buffer.data(mg + 114);
    const auto *mg_115 = buffer.data(mg + 115);
    const auto *mg_116 = buffer.data(mg + 116);
    const auto *mg_117 = buffer.data(mg + 117);
    const auto *mg_118 = buffer.data(mg + 118);
    const auto *mg_119 = buffer.data(mg + 119);
    const auto *mg_120 = buffer.data(mg + 120);
    const auto *mg_121 = buffer.data(mg + 121);
    const auto *mg_122 = buffer.data(mg + 122);
    const auto *mg_123 = buffer.data(mg + 123);
    const auto *mg_124 = buffer.data(mg + 124);
    const auto *mg_125 = buffer.data(mg + 125);
    const auto *mg_126 = buffer.data(mg + 126);
    const auto *mg_127 = buffer.data(mg + 127);
    const auto *mg_128 = buffer.data(mg + 128);
    const auto *mg_129 = buffer.data(mg + 129);
    const auto *mg_130 = buffer.data(mg + 130);
    const auto *mg_131 = buffer.data(mg + 131);
    const auto *mg_132 = buffer.data(mg + 132);
    const auto *mg_133 = buffer.data(mg + 133);
    const auto *mg_134 = buffer.data(mg + 134);
    const auto *mg_135 = buffer.data(mg + 135);
    const auto *mg_136 = buffer.data(mg + 136);
    const auto *mg_137 = buffer.data(mg + 137);
    const auto *mg_138 = buffer.data(mg + 138);
    const auto *mg_139 = buffer.data(mg + 139);
    const auto *mg_140 = buffer.data(mg + 140);
    const auto *mg_141 = buffer.data(mg + 141);
    const auto *mg_142 = buffer.data(mg + 142);
    const auto *mg_143 = buffer.data(mg + 143);
    const auto *mg_144 = buffer.data(mg + 144);
    const auto *mg_145 = buffer.data(mg + 145);
    const auto *mg_146 = buffer.data(mg + 146);
    const auto *mg_147 = buffer.data(mg + 147);
    const auto *mg_148 = buffer.data(mg + 148);
    const auto *mg_149 = buffer.data(mg + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, kg_0, kg_1, kg_2, kg_3, kg_4, mg_0, mg_1, \
                         mg_2, mg_3, mg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -8.0 * kg_0[k]
                 + f_0 * mg_0[k];

        t_1[k] = -8.0 * kg_1[k]
                 + f_0 * mg_1[k];

        t_2[k] = -8.0 * kg_2[k]
                 + f_0 * mg_2[k];

        t_3[k] = -8.0 * kg_3[k]
                 + f_0 * mg_3[k];

        t_4[k] = -8.0 * kg_4[k]
                 + f_0 * mg_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, kg_5, kg_6, kg_7, kg_8, kg_9, mg_5, mg_6, \
                         mg_7, mg_8, mg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -8.0 * kg_5[k]
                 + f_0 * mg_5[k];

        t_6[k] = -8.0 * kg_6[k]
                 + f_0 * mg_6[k];

        t_7[k] = -8.0 * kg_7[k]
                 + f_0 * mg_7[k];

        t_8[k] = -8.0 * kg_8[k]
                 + f_0 * mg_8[k];

        t_9[k] = -8.0 * kg_9[k]
                 + f_0 * mg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, kg_10, kg_11, kg_12, kg_13, kg_14, \
                         mg_10, mg_11, mg_12, mg_13, mg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -8.0 * kg_10[k]
                  + f_0 * mg_10[k];

        t_11[k] = -8.0 * kg_11[k]
                  + f_0 * mg_11[k];

        t_12[k] = -8.0 * kg_12[k]
                  + f_0 * mg_12[k];

        t_13[k] = -8.0 * kg_13[k]
                  + f_0 * mg_13[k];

        t_14[k] = -8.0 * kg_14[k]
                  + f_0 * mg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, kg_15, kg_16, kg_17, kg_18, kg_19, \
                         mg_15, mg_16, mg_17, mg_18, mg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -7.0 * kg_15[k]
                  + f_0 * mg_15[k];

        t_16[k] = -7.0 * kg_16[k]
                  + f_0 * mg_16[k];

        t_17[k] = -7.0 * kg_17[k]
                  + f_0 * mg_17[k];

        t_18[k] = -7.0 * kg_18[k]
                  + f_0 * mg_18[k];

        t_19[k] = -7.0 * kg_19[k]
                  + f_0 * mg_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, kg_20, kg_21, kg_22, kg_23, kg_24, \
                         mg_20, mg_21, mg_22, mg_23, mg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -7.0 * kg_20[k]
                  + f_0 * mg_20[k];

        t_21[k] = -7.0 * kg_21[k]
                  + f_0 * mg_21[k];

        t_22[k] = -7.0 * kg_22[k]
                  + f_0 * mg_22[k];

        t_23[k] = -7.0 * kg_23[k]
                  + f_0 * mg_23[k];

        t_24[k] = -7.0 * kg_24[k]
                  + f_0 * mg_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, kg_25, kg_26, kg_27, kg_28, kg_29, \
                         mg_25, mg_26, mg_27, mg_28, mg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -7.0 * kg_25[k]
                  + f_0 * mg_25[k];

        t_26[k] = -7.0 * kg_26[k]
                  + f_0 * mg_26[k];

        t_27[k] = -7.0 * kg_27[k]
                  + f_0 * mg_27[k];

        t_28[k] = -7.0 * kg_28[k]
                  + f_0 * mg_28[k];

        t_29[k] = -7.0 * kg_29[k]
                  + f_0 * mg_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, kg_30, kg_31, kg_32, kg_33, kg_34, \
                         mg_30, mg_31, mg_32, mg_33, mg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -7.0 * kg_30[k]
                  + f_0 * mg_30[k];

        t_31[k] = -7.0 * kg_31[k]
                  + f_0 * mg_31[k];

        t_32[k] = -7.0 * kg_32[k]
                  + f_0 * mg_32[k];

        t_33[k] = -7.0 * kg_33[k]
                  + f_0 * mg_33[k];

        t_34[k] = -7.0 * kg_34[k]
                  + f_0 * mg_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, kg_35, kg_36, kg_37, kg_38, kg_39, \
                         mg_35, mg_36, mg_37, mg_38, mg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -7.0 * kg_35[k]
                  + f_0 * mg_35[k];

        t_36[k] = -7.0 * kg_36[k]
                  + f_0 * mg_36[k];

        t_37[k] = -7.0 * kg_37[k]
                  + f_0 * mg_37[k];

        t_38[k] = -7.0 * kg_38[k]
                  + f_0 * mg_38[k];

        t_39[k] = -7.0 * kg_39[k]
                  + f_0 * mg_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, kg_40, kg_41, kg_42, kg_43, kg_44, \
                         mg_40, mg_41, mg_42, mg_43, mg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -7.0 * kg_40[k]
                  + f_0 * mg_40[k];

        t_41[k] = -7.0 * kg_41[k]
                  + f_0 * mg_41[k];

        t_42[k] = -7.0 * kg_42[k]
                  + f_0 * mg_42[k];

        t_43[k] = -7.0 * kg_43[k]
                  + f_0 * mg_43[k];

        t_44[k] = -7.0 * kg_44[k]
                  + f_0 * mg_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, kg_45, kg_46, kg_47, kg_48, kg_49, \
                         mg_45, mg_46, mg_47, mg_48, mg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -6.0 * kg_45[k]
                  + f_0 * mg_45[k];

        t_46[k] = -6.0 * kg_46[k]
                  + f_0 * mg_46[k];

        t_47[k] = -6.0 * kg_47[k]
                  + f_0 * mg_47[k];

        t_48[k] = -6.0 * kg_48[k]
                  + f_0 * mg_48[k];

        t_49[k] = -6.0 * kg_49[k]
                  + f_0 * mg_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, kg_50, kg_51, kg_52, kg_53, kg_54, \
                         mg_50, mg_51, mg_52, mg_53, mg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -6.0 * kg_50[k]
                  + f_0 * mg_50[k];

        t_51[k] = -6.0 * kg_51[k]
                  + f_0 * mg_51[k];

        t_52[k] = -6.0 * kg_52[k]
                  + f_0 * mg_52[k];

        t_53[k] = -6.0 * kg_53[k]
                  + f_0 * mg_53[k];

        t_54[k] = -6.0 * kg_54[k]
                  + f_0 * mg_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, kg_55, kg_56, kg_57, kg_58, kg_59, \
                         mg_55, mg_56, mg_57, mg_58, mg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -6.0 * kg_55[k]
                  + f_0 * mg_55[k];

        t_56[k] = -6.0 * kg_56[k]
                  + f_0 * mg_56[k];

        t_57[k] = -6.0 * kg_57[k]
                  + f_0 * mg_57[k];

        t_58[k] = -6.0 * kg_58[k]
                  + f_0 * mg_58[k];

        t_59[k] = -6.0 * kg_59[k]
                  + f_0 * mg_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, kg_60, kg_61, kg_62, kg_63, kg_64, \
                         mg_60, mg_61, mg_62, mg_63, mg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -6.0 * kg_60[k]
                  + f_0 * mg_60[k];

        t_61[k] = -6.0 * kg_61[k]
                  + f_0 * mg_61[k];

        t_62[k] = -6.0 * kg_62[k]
                  + f_0 * mg_62[k];

        t_63[k] = -6.0 * kg_63[k]
                  + f_0 * mg_63[k];

        t_64[k] = -6.0 * kg_64[k]
                  + f_0 * mg_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, kg_65, kg_66, kg_67, kg_68, kg_69, \
                         mg_65, mg_66, mg_67, mg_68, mg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -6.0 * kg_65[k]
                  + f_0 * mg_65[k];

        t_66[k] = -6.0 * kg_66[k]
                  + f_0 * mg_66[k];

        t_67[k] = -6.0 * kg_67[k]
                  + f_0 * mg_67[k];

        t_68[k] = -6.0 * kg_68[k]
                  + f_0 * mg_68[k];

        t_69[k] = -6.0 * kg_69[k]
                  + f_0 * mg_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, kg_70, kg_71, kg_72, kg_73, kg_74, \
                         mg_70, mg_71, mg_72, mg_73, mg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -6.0 * kg_70[k]
                  + f_0 * mg_70[k];

        t_71[k] = -6.0 * kg_71[k]
                  + f_0 * mg_71[k];

        t_72[k] = -6.0 * kg_72[k]
                  + f_0 * mg_72[k];

        t_73[k] = -6.0 * kg_73[k]
                  + f_0 * mg_73[k];

        t_74[k] = -6.0 * kg_74[k]
                  + f_0 * mg_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, kg_75, kg_76, kg_77, kg_78, kg_79, \
                         mg_75, mg_76, mg_77, mg_78, mg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -6.0 * kg_75[k]
                  + f_0 * mg_75[k];

        t_76[k] = -6.0 * kg_76[k]
                  + f_0 * mg_76[k];

        t_77[k] = -6.0 * kg_77[k]
                  + f_0 * mg_77[k];

        t_78[k] = -6.0 * kg_78[k]
                  + f_0 * mg_78[k];

        t_79[k] = -6.0 * kg_79[k]
                  + f_0 * mg_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, kg_80, kg_81, kg_82, kg_83, kg_84, \
                         mg_80, mg_81, mg_82, mg_83, mg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -6.0 * kg_80[k]
                  + f_0 * mg_80[k];

        t_81[k] = -6.0 * kg_81[k]
                  + f_0 * mg_81[k];

        t_82[k] = -6.0 * kg_82[k]
                  + f_0 * mg_82[k];

        t_83[k] = -6.0 * kg_83[k]
                  + f_0 * mg_83[k];

        t_84[k] = -6.0 * kg_84[k]
                  + f_0 * mg_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, kg_85, kg_86, kg_87, kg_88, kg_89, \
                         mg_85, mg_86, mg_87, mg_88, mg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -6.0 * kg_85[k]
                  + f_0 * mg_85[k];

        t_86[k] = -6.0 * kg_86[k]
                  + f_0 * mg_86[k];

        t_87[k] = -6.0 * kg_87[k]
                  + f_0 * mg_87[k];

        t_88[k] = -6.0 * kg_88[k]
                  + f_0 * mg_88[k];

        t_89[k] = -6.0 * kg_89[k]
                  + f_0 * mg_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, kg_90, kg_91, kg_92, kg_93, kg_94, \
                         mg_90, mg_91, mg_92, mg_93, mg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -5.0 * kg_90[k]
                  + f_0 * mg_90[k];

        t_91[k] = -5.0 * kg_91[k]
                  + f_0 * mg_91[k];

        t_92[k] = -5.0 * kg_92[k]
                  + f_0 * mg_92[k];

        t_93[k] = -5.0 * kg_93[k]
                  + f_0 * mg_93[k];

        t_94[k] = -5.0 * kg_94[k]
                  + f_0 * mg_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, kg_95, kg_96, kg_97, kg_98, kg_99, \
                         mg_95, mg_96, mg_97, mg_98, mg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -5.0 * kg_95[k]
                  + f_0 * mg_95[k];

        t_96[k] = -5.0 * kg_96[k]
                  + f_0 * mg_96[k];

        t_97[k] = -5.0 * kg_97[k]
                  + f_0 * mg_97[k];

        t_98[k] = -5.0 * kg_98[k]
                  + f_0 * mg_98[k];

        t_99[k] = -5.0 * kg_99[k]
                  + f_0 * mg_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, kg_100, kg_101, kg_102, kg_103, \
                         kg_104, mg_100, mg_101, mg_102, mg_103, \
                         mg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -5.0 * kg_100[k]
                   + f_0 * mg_100[k];

        t_101[k] = -5.0 * kg_101[k]
                   + f_0 * mg_101[k];

        t_102[k] = -5.0 * kg_102[k]
                   + f_0 * mg_102[k];

        t_103[k] = -5.0 * kg_103[k]
                   + f_0 * mg_103[k];

        t_104[k] = -5.0 * kg_104[k]
                   + f_0 * mg_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, kg_105, kg_106, kg_107, kg_108, \
                         kg_109, mg_105, mg_106, mg_107, mg_108, \
                         mg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -5.0 * kg_105[k]
                   + f_0 * mg_105[k];

        t_106[k] = -5.0 * kg_106[k]
                   + f_0 * mg_106[k];

        t_107[k] = -5.0 * kg_107[k]
                   + f_0 * mg_107[k];

        t_108[k] = -5.0 * kg_108[k]
                   + f_0 * mg_108[k];

        t_109[k] = -5.0 * kg_109[k]
                   + f_0 * mg_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, kg_110, kg_111, kg_112, kg_113, \
                         kg_114, mg_110, mg_111, mg_112, mg_113, \
                         mg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -5.0 * kg_110[k]
                   + f_0 * mg_110[k];

        t_111[k] = -5.0 * kg_111[k]
                   + f_0 * mg_111[k];

        t_112[k] = -5.0 * kg_112[k]
                   + f_0 * mg_112[k];

        t_113[k] = -5.0 * kg_113[k]
                   + f_0 * mg_113[k];

        t_114[k] = -5.0 * kg_114[k]
                   + f_0 * mg_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, kg_115, kg_116, kg_117, kg_118, \
                         kg_119, mg_115, mg_116, mg_117, mg_118, \
                         mg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -5.0 * kg_115[k]
                   + f_0 * mg_115[k];

        t_116[k] = -5.0 * kg_116[k]
                   + f_0 * mg_116[k];

        t_117[k] = -5.0 * kg_117[k]
                   + f_0 * mg_117[k];

        t_118[k] = -5.0 * kg_118[k]
                   + f_0 * mg_118[k];

        t_119[k] = -5.0 * kg_119[k]
                   + f_0 * mg_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, kg_120, kg_121, kg_122, kg_123, \
                         kg_124, mg_120, mg_121, mg_122, mg_123, \
                         mg_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -5.0 * kg_120[k]
                   + f_0 * mg_120[k];

        t_121[k] = -5.0 * kg_121[k]
                   + f_0 * mg_121[k];

        t_122[k] = -5.0 * kg_122[k]
                   + f_0 * mg_122[k];

        t_123[k] = -5.0 * kg_123[k]
                   + f_0 * mg_123[k];

        t_124[k] = -5.0 * kg_124[k]
                   + f_0 * mg_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, kg_125, kg_126, kg_127, kg_128, \
                         kg_129, mg_125, mg_126, mg_127, mg_128, \
                         mg_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -5.0 * kg_125[k]
                   + f_0 * mg_125[k];

        t_126[k] = -5.0 * kg_126[k]
                   + f_0 * mg_126[k];

        t_127[k] = -5.0 * kg_127[k]
                   + f_0 * mg_127[k];

        t_128[k] = -5.0 * kg_128[k]
                   + f_0 * mg_128[k];

        t_129[k] = -5.0 * kg_129[k]
                   + f_0 * mg_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, kg_130, kg_131, kg_132, kg_133, \
                         kg_134, mg_130, mg_131, mg_132, mg_133, \
                         mg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -5.0 * kg_130[k]
                   + f_0 * mg_130[k];

        t_131[k] = -5.0 * kg_131[k]
                   + f_0 * mg_131[k];

        t_132[k] = -5.0 * kg_132[k]
                   + f_0 * mg_132[k];

        t_133[k] = -5.0 * kg_133[k]
                   + f_0 * mg_133[k];

        t_134[k] = -5.0 * kg_134[k]
                   + f_0 * mg_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, kg_135, kg_136, kg_137, kg_138, \
                         kg_139, mg_135, mg_136, mg_137, mg_138, \
                         mg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -5.0 * kg_135[k]
                   + f_0 * mg_135[k];

        t_136[k] = -5.0 * kg_136[k]
                   + f_0 * mg_136[k];

        t_137[k] = -5.0 * kg_137[k]
                   + f_0 * mg_137[k];

        t_138[k] = -5.0 * kg_138[k]
                   + f_0 * mg_138[k];

        t_139[k] = -5.0 * kg_139[k]
                   + f_0 * mg_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, kg_140, kg_141, kg_142, kg_143, \
                         kg_144, mg_140, mg_141, mg_142, mg_143, \
                         mg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -5.0 * kg_140[k]
                   + f_0 * mg_140[k];

        t_141[k] = -5.0 * kg_141[k]
                   + f_0 * mg_141[k];

        t_142[k] = -5.0 * kg_142[k]
                   + f_0 * mg_142[k];

        t_143[k] = -5.0 * kg_143[k]
                   + f_0 * mg_143[k];

        t_144[k] = -5.0 * kg_144[k]
                   + f_0 * mg_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, kg_145, kg_146, kg_147, kg_148, \
                         kg_149, mg_145, mg_146, mg_147, mg_148, \
                         mg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -5.0 * kg_145[k]
                   + f_0 * mg_145[k];

        t_146[k] = -5.0 * kg_146[k]
                   + f_0 * mg_146[k];

        t_147[k] = -5.0 * kg_147[k]
                   + f_0 * mg_147[k];

        t_148[k] = -5.0 * kg_148[k]
                   + f_0 * mg_148[k];

        t_149[k] = -5.0 * kg_149[k]
                   + f_0 * mg_149[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
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

    const auto *mg_150 = buffer.data(mg + 150);
    const auto *mg_151 = buffer.data(mg + 151);
    const auto *mg_152 = buffer.data(mg + 152);
    const auto *mg_153 = buffer.data(mg + 153);
    const auto *mg_154 = buffer.data(mg + 154);
    const auto *mg_155 = buffer.data(mg + 155);
    const auto *mg_156 = buffer.data(mg + 156);
    const auto *mg_157 = buffer.data(mg + 157);
    const auto *mg_158 = buffer.data(mg + 158);
    const auto *mg_159 = buffer.data(mg + 159);
    const auto *mg_160 = buffer.data(mg + 160);
    const auto *mg_161 = buffer.data(mg + 161);
    const auto *mg_162 = buffer.data(mg + 162);
    const auto *mg_163 = buffer.data(mg + 163);
    const auto *mg_164 = buffer.data(mg + 164);
    const auto *mg_165 = buffer.data(mg + 165);
    const auto *mg_166 = buffer.data(mg + 166);
    const auto *mg_167 = buffer.data(mg + 167);
    const auto *mg_168 = buffer.data(mg + 168);
    const auto *mg_169 = buffer.data(mg + 169);
    const auto *mg_170 = buffer.data(mg + 170);
    const auto *mg_171 = buffer.data(mg + 171);
    const auto *mg_172 = buffer.data(mg + 172);
    const auto *mg_173 = buffer.data(mg + 173);
    const auto *mg_174 = buffer.data(mg + 174);
    const auto *mg_175 = buffer.data(mg + 175);
    const auto *mg_176 = buffer.data(mg + 176);
    const auto *mg_177 = buffer.data(mg + 177);
    const auto *mg_178 = buffer.data(mg + 178);
    const auto *mg_179 = buffer.data(mg + 179);
    const auto *mg_180 = buffer.data(mg + 180);
    const auto *mg_181 = buffer.data(mg + 181);
    const auto *mg_182 = buffer.data(mg + 182);
    const auto *mg_183 = buffer.data(mg + 183);
    const auto *mg_184 = buffer.data(mg + 184);
    const auto *mg_185 = buffer.data(mg + 185);
    const auto *mg_186 = buffer.data(mg + 186);
    const auto *mg_187 = buffer.data(mg + 187);
    const auto *mg_188 = buffer.data(mg + 188);
    const auto *mg_189 = buffer.data(mg + 189);
    const auto *mg_190 = buffer.data(mg + 190);
    const auto *mg_191 = buffer.data(mg + 191);
    const auto *mg_192 = buffer.data(mg + 192);
    const auto *mg_193 = buffer.data(mg + 193);
    const auto *mg_194 = buffer.data(mg + 194);
    const auto *mg_195 = buffer.data(mg + 195);
    const auto *mg_196 = buffer.data(mg + 196);
    const auto *mg_197 = buffer.data(mg + 197);
    const auto *mg_198 = buffer.data(mg + 198);
    const auto *mg_199 = buffer.data(mg + 199);
    const auto *mg_200 = buffer.data(mg + 200);
    const auto *mg_201 = buffer.data(mg + 201);
    const auto *mg_202 = buffer.data(mg + 202);
    const auto *mg_203 = buffer.data(mg + 203);
    const auto *mg_204 = buffer.data(mg + 204);
    const auto *mg_205 = buffer.data(mg + 205);
    const auto *mg_206 = buffer.data(mg + 206);
    const auto *mg_207 = buffer.data(mg + 207);
    const auto *mg_208 = buffer.data(mg + 208);
    const auto *mg_209 = buffer.data(mg + 209);
    const auto *mg_210 = buffer.data(mg + 210);
    const auto *mg_211 = buffer.data(mg + 211);
    const auto *mg_212 = buffer.data(mg + 212);
    const auto *mg_213 = buffer.data(mg + 213);
    const auto *mg_214 = buffer.data(mg + 214);
    const auto *mg_215 = buffer.data(mg + 215);
    const auto *mg_216 = buffer.data(mg + 216);
    const auto *mg_217 = buffer.data(mg + 217);
    const auto *mg_218 = buffer.data(mg + 218);
    const auto *mg_219 = buffer.data(mg + 219);
    const auto *mg_220 = buffer.data(mg + 220);
    const auto *mg_221 = buffer.data(mg + 221);
    const auto *mg_222 = buffer.data(mg + 222);
    const auto *mg_223 = buffer.data(mg + 223);
    const auto *mg_224 = buffer.data(mg + 224);
    const auto *mg_225 = buffer.data(mg + 225);
    const auto *mg_226 = buffer.data(mg + 226);
    const auto *mg_227 = buffer.data(mg + 227);
    const auto *mg_228 = buffer.data(mg + 228);
    const auto *mg_229 = buffer.data(mg + 229);
    const auto *mg_230 = buffer.data(mg + 230);
    const auto *mg_231 = buffer.data(mg + 231);
    const auto *mg_232 = buffer.data(mg + 232);
    const auto *mg_233 = buffer.data(mg + 233);
    const auto *mg_234 = buffer.data(mg + 234);
    const auto *mg_235 = buffer.data(mg + 235);
    const auto *mg_236 = buffer.data(mg + 236);
    const auto *mg_237 = buffer.data(mg + 237);
    const auto *mg_238 = buffer.data(mg + 238);
    const auto *mg_239 = buffer.data(mg + 239);
    const auto *mg_240 = buffer.data(mg + 240);
    const auto *mg_241 = buffer.data(mg + 241);
    const auto *mg_242 = buffer.data(mg + 242);
    const auto *mg_243 = buffer.data(mg + 243);
    const auto *mg_244 = buffer.data(mg + 244);
    const auto *mg_245 = buffer.data(mg + 245);
    const auto *mg_246 = buffer.data(mg + 246);
    const auto *mg_247 = buffer.data(mg + 247);
    const auto *mg_248 = buffer.data(mg + 248);
    const auto *mg_249 = buffer.data(mg + 249);
    const auto *mg_250 = buffer.data(mg + 250);
    const auto *mg_251 = buffer.data(mg + 251);
    const auto *mg_252 = buffer.data(mg + 252);
    const auto *mg_253 = buffer.data(mg + 253);
    const auto *mg_254 = buffer.data(mg + 254);
    const auto *mg_255 = buffer.data(mg + 255);
    const auto *mg_256 = buffer.data(mg + 256);
    const auto *mg_257 = buffer.data(mg + 257);
    const auto *mg_258 = buffer.data(mg + 258);
    const auto *mg_259 = buffer.data(mg + 259);
    const auto *mg_260 = buffer.data(mg + 260);
    const auto *mg_261 = buffer.data(mg + 261);
    const auto *mg_262 = buffer.data(mg + 262);
    const auto *mg_263 = buffer.data(mg + 263);
    const auto *mg_264 = buffer.data(mg + 264);
    const auto *mg_265 = buffer.data(mg + 265);
    const auto *mg_266 = buffer.data(mg + 266);
    const auto *mg_267 = buffer.data(mg + 267);
    const auto *mg_268 = buffer.data(mg + 268);
    const auto *mg_269 = buffer.data(mg + 269);
    const auto *mg_270 = buffer.data(mg + 270);
    const auto *mg_271 = buffer.data(mg + 271);
    const auto *mg_272 = buffer.data(mg + 272);
    const auto *mg_273 = buffer.data(mg + 273);
    const auto *mg_274 = buffer.data(mg + 274);
    const auto *mg_275 = buffer.data(mg + 275);
    const auto *mg_276 = buffer.data(mg + 276);
    const auto *mg_277 = buffer.data(mg + 277);
    const auto *mg_278 = buffer.data(mg + 278);
    const auto *mg_279 = buffer.data(mg + 279);
    const auto *mg_280 = buffer.data(mg + 280);
    const auto *mg_281 = buffer.data(mg + 281);
    const auto *mg_282 = buffer.data(mg + 282);
    const auto *mg_283 = buffer.data(mg + 283);
    const auto *mg_284 = buffer.data(mg + 284);
    const auto *mg_285 = buffer.data(mg + 285);
    const auto *mg_286 = buffer.data(mg + 286);
    const auto *mg_287 = buffer.data(mg + 287);
    const auto *mg_288 = buffer.data(mg + 288);
    const auto *mg_289 = buffer.data(mg + 289);
    const auto *mg_290 = buffer.data(mg + 290);
    const auto *mg_291 = buffer.data(mg + 291);
    const auto *mg_292 = buffer.data(mg + 292);
    const auto *mg_293 = buffer.data(mg + 293);
    const auto *mg_294 = buffer.data(mg + 294);
    const auto *mg_295 = buffer.data(mg + 295);
    const auto *mg_296 = buffer.data(mg + 296);
    const auto *mg_297 = buffer.data(mg + 297);
    const auto *mg_298 = buffer.data(mg + 298);
    const auto *mg_299 = buffer.data(mg + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, kg_150, kg_151, kg_152, kg_153, \
                         kg_154, mg_150, mg_151, mg_152, mg_153, \
                         mg_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -4.0 * kg_150[k]
                   + f_0 * mg_150[k];

        t_151[k] = -4.0 * kg_151[k]
                   + f_0 * mg_151[k];

        t_152[k] = -4.0 * kg_152[k]
                   + f_0 * mg_152[k];

        t_153[k] = -4.0 * kg_153[k]
                   + f_0 * mg_153[k];

        t_154[k] = -4.0 * kg_154[k]
                   + f_0 * mg_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, kg_155, kg_156, kg_157, kg_158, \
                         kg_159, mg_155, mg_156, mg_157, mg_158, \
                         mg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -4.0 * kg_155[k]
                   + f_0 * mg_155[k];

        t_156[k] = -4.0 * kg_156[k]
                   + f_0 * mg_156[k];

        t_157[k] = -4.0 * kg_157[k]
                   + f_0 * mg_157[k];

        t_158[k] = -4.0 * kg_158[k]
                   + f_0 * mg_158[k];

        t_159[k] = -4.0 * kg_159[k]
                   + f_0 * mg_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, kg_160, kg_161, kg_162, kg_163, \
                         kg_164, mg_160, mg_161, mg_162, mg_163, \
                         mg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -4.0 * kg_160[k]
                   + f_0 * mg_160[k];

        t_161[k] = -4.0 * kg_161[k]
                   + f_0 * mg_161[k];

        t_162[k] = -4.0 * kg_162[k]
                   + f_0 * mg_162[k];

        t_163[k] = -4.0 * kg_163[k]
                   + f_0 * mg_163[k];

        t_164[k] = -4.0 * kg_164[k]
                   + f_0 * mg_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, kg_165, kg_166, kg_167, kg_168, \
                         kg_169, mg_165, mg_166, mg_167, mg_168, \
                         mg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -4.0 * kg_165[k]
                   + f_0 * mg_165[k];

        t_166[k] = -4.0 * kg_166[k]
                   + f_0 * mg_166[k];

        t_167[k] = -4.0 * kg_167[k]
                   + f_0 * mg_167[k];

        t_168[k] = -4.0 * kg_168[k]
                   + f_0 * mg_168[k];

        t_169[k] = -4.0 * kg_169[k]
                   + f_0 * mg_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, kg_170, kg_171, kg_172, kg_173, \
                         kg_174, mg_170, mg_171, mg_172, mg_173, \
                         mg_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -4.0 * kg_170[k]
                   + f_0 * mg_170[k];

        t_171[k] = -4.0 * kg_171[k]
                   + f_0 * mg_171[k];

        t_172[k] = -4.0 * kg_172[k]
                   + f_0 * mg_172[k];

        t_173[k] = -4.0 * kg_173[k]
                   + f_0 * mg_173[k];

        t_174[k] = -4.0 * kg_174[k]
                   + f_0 * mg_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, kg_175, kg_176, kg_177, kg_178, \
                         kg_179, mg_175, mg_176, mg_177, mg_178, \
                         mg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -4.0 * kg_175[k]
                   + f_0 * mg_175[k];

        t_176[k] = -4.0 * kg_176[k]
                   + f_0 * mg_176[k];

        t_177[k] = -4.0 * kg_177[k]
                   + f_0 * mg_177[k];

        t_178[k] = -4.0 * kg_178[k]
                   + f_0 * mg_178[k];

        t_179[k] = -4.0 * kg_179[k]
                   + f_0 * mg_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, kg_180, kg_181, kg_182, kg_183, \
                         kg_184, mg_180, mg_181, mg_182, mg_183, \
                         mg_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -4.0 * kg_180[k]
                   + f_0 * mg_180[k];

        t_181[k] = -4.0 * kg_181[k]
                   + f_0 * mg_181[k];

        t_182[k] = -4.0 * kg_182[k]
                   + f_0 * mg_182[k];

        t_183[k] = -4.0 * kg_183[k]
                   + f_0 * mg_183[k];

        t_184[k] = -4.0 * kg_184[k]
                   + f_0 * mg_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, kg_185, kg_186, kg_187, kg_188, \
                         kg_189, mg_185, mg_186, mg_187, mg_188, \
                         mg_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -4.0 * kg_185[k]
                   + f_0 * mg_185[k];

        t_186[k] = -4.0 * kg_186[k]
                   + f_0 * mg_186[k];

        t_187[k] = -4.0 * kg_187[k]
                   + f_0 * mg_187[k];

        t_188[k] = -4.0 * kg_188[k]
                   + f_0 * mg_188[k];

        t_189[k] = -4.0 * kg_189[k]
                   + f_0 * mg_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, kg_190, kg_191, kg_192, kg_193, \
                         kg_194, mg_190, mg_191, mg_192, mg_193, \
                         mg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -4.0 * kg_190[k]
                   + f_0 * mg_190[k];

        t_191[k] = -4.0 * kg_191[k]
                   + f_0 * mg_191[k];

        t_192[k] = -4.0 * kg_192[k]
                   + f_0 * mg_192[k];

        t_193[k] = -4.0 * kg_193[k]
                   + f_0 * mg_193[k];

        t_194[k] = -4.0 * kg_194[k]
                   + f_0 * mg_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, kg_195, kg_196, kg_197, kg_198, \
                         kg_199, mg_195, mg_196, mg_197, mg_198, \
                         mg_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -4.0 * kg_195[k]
                   + f_0 * mg_195[k];

        t_196[k] = -4.0 * kg_196[k]
                   + f_0 * mg_196[k];

        t_197[k] = -4.0 * kg_197[k]
                   + f_0 * mg_197[k];

        t_198[k] = -4.0 * kg_198[k]
                   + f_0 * mg_198[k];

        t_199[k] = -4.0 * kg_199[k]
                   + f_0 * mg_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, kg_200, kg_201, kg_202, kg_203, \
                         kg_204, mg_200, mg_201, mg_202, mg_203, \
                         mg_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -4.0 * kg_200[k]
                   + f_0 * mg_200[k];

        t_201[k] = -4.0 * kg_201[k]
                   + f_0 * mg_201[k];

        t_202[k] = -4.0 * kg_202[k]
                   + f_0 * mg_202[k];

        t_203[k] = -4.0 * kg_203[k]
                   + f_0 * mg_203[k];

        t_204[k] = -4.0 * kg_204[k]
                   + f_0 * mg_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, kg_205, kg_206, kg_207, kg_208, \
                         kg_209, mg_205, mg_206, mg_207, mg_208, \
                         mg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -4.0 * kg_205[k]
                   + f_0 * mg_205[k];

        t_206[k] = -4.0 * kg_206[k]
                   + f_0 * mg_206[k];

        t_207[k] = -4.0 * kg_207[k]
                   + f_0 * mg_207[k];

        t_208[k] = -4.0 * kg_208[k]
                   + f_0 * mg_208[k];

        t_209[k] = -4.0 * kg_209[k]
                   + f_0 * mg_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, kg_210, kg_211, kg_212, kg_213, \
                         kg_214, mg_210, mg_211, mg_212, mg_213, \
                         mg_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -4.0 * kg_210[k]
                   + f_0 * mg_210[k];

        t_211[k] = -4.0 * kg_211[k]
                   + f_0 * mg_211[k];

        t_212[k] = -4.0 * kg_212[k]
                   + f_0 * mg_212[k];

        t_213[k] = -4.0 * kg_213[k]
                   + f_0 * mg_213[k];

        t_214[k] = -4.0 * kg_214[k]
                   + f_0 * mg_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, kg_215, kg_216, kg_217, kg_218, \
                         kg_219, mg_215, mg_216, mg_217, mg_218, \
                         mg_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -4.0 * kg_215[k]
                   + f_0 * mg_215[k];

        t_216[k] = -4.0 * kg_216[k]
                   + f_0 * mg_216[k];

        t_217[k] = -4.0 * kg_217[k]
                   + f_0 * mg_217[k];

        t_218[k] = -4.0 * kg_218[k]
                   + f_0 * mg_218[k];

        t_219[k] = -4.0 * kg_219[k]
                   + f_0 * mg_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, kg_220, kg_221, kg_222, kg_223, \
                         kg_224, mg_220, mg_221, mg_222, mg_223, \
                         mg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -4.0 * kg_220[k]
                   + f_0 * mg_220[k];

        t_221[k] = -4.0 * kg_221[k]
                   + f_0 * mg_221[k];

        t_222[k] = -4.0 * kg_222[k]
                   + f_0 * mg_222[k];

        t_223[k] = -4.0 * kg_223[k]
                   + f_0 * mg_223[k];

        t_224[k] = -4.0 * kg_224[k]
                   + f_0 * mg_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, kg_225, kg_226, kg_227, kg_228, \
                         kg_229, mg_225, mg_226, mg_227, mg_228, \
                         mg_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -3.0 * kg_225[k]
                   + f_0 * mg_225[k];

        t_226[k] = -3.0 * kg_226[k]
                   + f_0 * mg_226[k];

        t_227[k] = -3.0 * kg_227[k]
                   + f_0 * mg_227[k];

        t_228[k] = -3.0 * kg_228[k]
                   + f_0 * mg_228[k];

        t_229[k] = -3.0 * kg_229[k]
                   + f_0 * mg_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, kg_230, kg_231, kg_232, kg_233, \
                         kg_234, mg_230, mg_231, mg_232, mg_233, \
                         mg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -3.0 * kg_230[k]
                   + f_0 * mg_230[k];

        t_231[k] = -3.0 * kg_231[k]
                   + f_0 * mg_231[k];

        t_232[k] = -3.0 * kg_232[k]
                   + f_0 * mg_232[k];

        t_233[k] = -3.0 * kg_233[k]
                   + f_0 * mg_233[k];

        t_234[k] = -3.0 * kg_234[k]
                   + f_0 * mg_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, kg_235, kg_236, kg_237, kg_238, \
                         kg_239, mg_235, mg_236, mg_237, mg_238, \
                         mg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -3.0 * kg_235[k]
                   + f_0 * mg_235[k];

        t_236[k] = -3.0 * kg_236[k]
                   + f_0 * mg_236[k];

        t_237[k] = -3.0 * kg_237[k]
                   + f_0 * mg_237[k];

        t_238[k] = -3.0 * kg_238[k]
                   + f_0 * mg_238[k];

        t_239[k] = -3.0 * kg_239[k]
                   + f_0 * mg_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, kg_240, kg_241, kg_242, kg_243, \
                         kg_244, mg_240, mg_241, mg_242, mg_243, \
                         mg_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * kg_240[k]
                   + f_0 * mg_240[k];

        t_241[k] = -3.0 * kg_241[k]
                   + f_0 * mg_241[k];

        t_242[k] = -3.0 * kg_242[k]
                   + f_0 * mg_242[k];

        t_243[k] = -3.0 * kg_243[k]
                   + f_0 * mg_243[k];

        t_244[k] = -3.0 * kg_244[k]
                   + f_0 * mg_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, kg_245, kg_246, kg_247, kg_248, \
                         kg_249, mg_245, mg_246, mg_247, mg_248, \
                         mg_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * kg_245[k]
                   + f_0 * mg_245[k];

        t_246[k] = -3.0 * kg_246[k]
                   + f_0 * mg_246[k];

        t_247[k] = -3.0 * kg_247[k]
                   + f_0 * mg_247[k];

        t_248[k] = -3.0 * kg_248[k]
                   + f_0 * mg_248[k];

        t_249[k] = -3.0 * kg_249[k]
                   + f_0 * mg_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, kg_250, kg_251, kg_252, kg_253, \
                         kg_254, mg_250, mg_251, mg_252, mg_253, \
                         mg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -3.0 * kg_250[k]
                   + f_0 * mg_250[k];

        t_251[k] = -3.0 * kg_251[k]
                   + f_0 * mg_251[k];

        t_252[k] = -3.0 * kg_252[k]
                   + f_0 * mg_252[k];

        t_253[k] = -3.0 * kg_253[k]
                   + f_0 * mg_253[k];

        t_254[k] = -3.0 * kg_254[k]
                   + f_0 * mg_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, kg_255, kg_256, kg_257, kg_258, \
                         kg_259, mg_255, mg_256, mg_257, mg_258, \
                         mg_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -3.0 * kg_255[k]
                   + f_0 * mg_255[k];

        t_256[k] = -3.0 * kg_256[k]
                   + f_0 * mg_256[k];

        t_257[k] = -3.0 * kg_257[k]
                   + f_0 * mg_257[k];

        t_258[k] = -3.0 * kg_258[k]
                   + f_0 * mg_258[k];

        t_259[k] = -3.0 * kg_259[k]
                   + f_0 * mg_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, kg_260, kg_261, kg_262, kg_263, \
                         kg_264, mg_260, mg_261, mg_262, mg_263, \
                         mg_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -3.0 * kg_260[k]
                   + f_0 * mg_260[k];

        t_261[k] = -3.0 * kg_261[k]
                   + f_0 * mg_261[k];

        t_262[k] = -3.0 * kg_262[k]
                   + f_0 * mg_262[k];

        t_263[k] = -3.0 * kg_263[k]
                   + f_0 * mg_263[k];

        t_264[k] = -3.0 * kg_264[k]
                   + f_0 * mg_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, kg_265, kg_266, kg_267, kg_268, \
                         kg_269, mg_265, mg_266, mg_267, mg_268, \
                         mg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -3.0 * kg_265[k]
                   + f_0 * mg_265[k];

        t_266[k] = -3.0 * kg_266[k]
                   + f_0 * mg_266[k];

        t_267[k] = -3.0 * kg_267[k]
                   + f_0 * mg_267[k];

        t_268[k] = -3.0 * kg_268[k]
                   + f_0 * mg_268[k];

        t_269[k] = -3.0 * kg_269[k]
                   + f_0 * mg_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, kg_270, kg_271, kg_272, kg_273, \
                         kg_274, mg_270, mg_271, mg_272, mg_273, \
                         mg_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -3.0 * kg_270[k]
                   + f_0 * mg_270[k];

        t_271[k] = -3.0 * kg_271[k]
                   + f_0 * mg_271[k];

        t_272[k] = -3.0 * kg_272[k]
                   + f_0 * mg_272[k];

        t_273[k] = -3.0 * kg_273[k]
                   + f_0 * mg_273[k];

        t_274[k] = -3.0 * kg_274[k]
                   + f_0 * mg_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, kg_275, kg_276, kg_277, kg_278, \
                         kg_279, mg_275, mg_276, mg_277, mg_278, \
                         mg_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * kg_275[k]
                   + f_0 * mg_275[k];

        t_276[k] = -3.0 * kg_276[k]
                   + f_0 * mg_276[k];

        t_277[k] = -3.0 * kg_277[k]
                   + f_0 * mg_277[k];

        t_278[k] = -3.0 * kg_278[k]
                   + f_0 * mg_278[k];

        t_279[k] = -3.0 * kg_279[k]
                   + f_0 * mg_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, kg_280, kg_281, kg_282, kg_283, \
                         kg_284, mg_280, mg_281, mg_282, mg_283, \
                         mg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * kg_280[k]
                   + f_0 * mg_280[k];

        t_281[k] = -3.0 * kg_281[k]
                   + f_0 * mg_281[k];

        t_282[k] = -3.0 * kg_282[k]
                   + f_0 * mg_282[k];

        t_283[k] = -3.0 * kg_283[k]
                   + f_0 * mg_283[k];

        t_284[k] = -3.0 * kg_284[k]
                   + f_0 * mg_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, kg_285, kg_286, kg_287, kg_288, \
                         kg_289, mg_285, mg_286, mg_287, mg_288, \
                         mg_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * kg_285[k]
                   + f_0 * mg_285[k];

        t_286[k] = -3.0 * kg_286[k]
                   + f_0 * mg_286[k];

        t_287[k] = -3.0 * kg_287[k]
                   + f_0 * mg_287[k];

        t_288[k] = -3.0 * kg_288[k]
                   + f_0 * mg_288[k];

        t_289[k] = -3.0 * kg_289[k]
                   + f_0 * mg_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, kg_290, kg_291, kg_292, kg_293, \
                         kg_294, mg_290, mg_291, mg_292, mg_293, \
                         mg_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * kg_290[k]
                   + f_0 * mg_290[k];

        t_291[k] = -3.0 * kg_291[k]
                   + f_0 * mg_291[k];

        t_292[k] = -3.0 * kg_292[k]
                   + f_0 * mg_292[k];

        t_293[k] = -3.0 * kg_293[k]
                   + f_0 * mg_293[k];

        t_294[k] = -3.0 * kg_294[k]
                   + f_0 * mg_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, kg_295, kg_296, kg_297, kg_298, \
                         kg_299, mg_295, mg_296, mg_297, mg_298, \
                         mg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -3.0 * kg_295[k]
                   + f_0 * mg_295[k];

        t_296[k] = -3.0 * kg_296[k]
                   + f_0 * mg_296[k];

        t_297[k] = -3.0 * kg_297[k]
                   + f_0 * mg_297[k];

        t_298[k] = -3.0 * kg_298[k]
                   + f_0 * mg_298[k];

        t_299[k] = -3.0 * kg_299[k]
                   + f_0 * mg_299[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
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

    const auto *mg_300 = buffer.data(mg + 300);
    const auto *mg_301 = buffer.data(mg + 301);
    const auto *mg_302 = buffer.data(mg + 302);
    const auto *mg_303 = buffer.data(mg + 303);
    const auto *mg_304 = buffer.data(mg + 304);
    const auto *mg_305 = buffer.data(mg + 305);
    const auto *mg_306 = buffer.data(mg + 306);
    const auto *mg_307 = buffer.data(mg + 307);
    const auto *mg_308 = buffer.data(mg + 308);
    const auto *mg_309 = buffer.data(mg + 309);
    const auto *mg_310 = buffer.data(mg + 310);
    const auto *mg_311 = buffer.data(mg + 311);
    const auto *mg_312 = buffer.data(mg + 312);
    const auto *mg_313 = buffer.data(mg + 313);
    const auto *mg_314 = buffer.data(mg + 314);
    const auto *mg_315 = buffer.data(mg + 315);
    const auto *mg_316 = buffer.data(mg + 316);
    const auto *mg_317 = buffer.data(mg + 317);
    const auto *mg_318 = buffer.data(mg + 318);
    const auto *mg_319 = buffer.data(mg + 319);
    const auto *mg_320 = buffer.data(mg + 320);
    const auto *mg_321 = buffer.data(mg + 321);
    const auto *mg_322 = buffer.data(mg + 322);
    const auto *mg_323 = buffer.data(mg + 323);
    const auto *mg_324 = buffer.data(mg + 324);
    const auto *mg_325 = buffer.data(mg + 325);
    const auto *mg_326 = buffer.data(mg + 326);
    const auto *mg_327 = buffer.data(mg + 327);
    const auto *mg_328 = buffer.data(mg + 328);
    const auto *mg_329 = buffer.data(mg + 329);
    const auto *mg_330 = buffer.data(mg + 330);
    const auto *mg_331 = buffer.data(mg + 331);
    const auto *mg_332 = buffer.data(mg + 332);
    const auto *mg_333 = buffer.data(mg + 333);
    const auto *mg_334 = buffer.data(mg + 334);
    const auto *mg_335 = buffer.data(mg + 335);
    const auto *mg_336 = buffer.data(mg + 336);
    const auto *mg_337 = buffer.data(mg + 337);
    const auto *mg_338 = buffer.data(mg + 338);
    const auto *mg_339 = buffer.data(mg + 339);
    const auto *mg_340 = buffer.data(mg + 340);
    const auto *mg_341 = buffer.data(mg + 341);
    const auto *mg_342 = buffer.data(mg + 342);
    const auto *mg_343 = buffer.data(mg + 343);
    const auto *mg_344 = buffer.data(mg + 344);
    const auto *mg_345 = buffer.data(mg + 345);
    const auto *mg_346 = buffer.data(mg + 346);
    const auto *mg_347 = buffer.data(mg + 347);
    const auto *mg_348 = buffer.data(mg + 348);
    const auto *mg_349 = buffer.data(mg + 349);
    const auto *mg_350 = buffer.data(mg + 350);
    const auto *mg_351 = buffer.data(mg + 351);
    const auto *mg_352 = buffer.data(mg + 352);
    const auto *mg_353 = buffer.data(mg + 353);
    const auto *mg_354 = buffer.data(mg + 354);
    const auto *mg_355 = buffer.data(mg + 355);
    const auto *mg_356 = buffer.data(mg + 356);
    const auto *mg_357 = buffer.data(mg + 357);
    const auto *mg_358 = buffer.data(mg + 358);
    const auto *mg_359 = buffer.data(mg + 359);
    const auto *mg_360 = buffer.data(mg + 360);
    const auto *mg_361 = buffer.data(mg + 361);
    const auto *mg_362 = buffer.data(mg + 362);
    const auto *mg_363 = buffer.data(mg + 363);
    const auto *mg_364 = buffer.data(mg + 364);
    const auto *mg_365 = buffer.data(mg + 365);
    const auto *mg_366 = buffer.data(mg + 366);
    const auto *mg_367 = buffer.data(mg + 367);
    const auto *mg_368 = buffer.data(mg + 368);
    const auto *mg_369 = buffer.data(mg + 369);
    const auto *mg_370 = buffer.data(mg + 370);
    const auto *mg_371 = buffer.data(mg + 371);
    const auto *mg_372 = buffer.data(mg + 372);
    const auto *mg_373 = buffer.data(mg + 373);
    const auto *mg_374 = buffer.data(mg + 374);
    const auto *mg_375 = buffer.data(mg + 375);
    const auto *mg_376 = buffer.data(mg + 376);
    const auto *mg_377 = buffer.data(mg + 377);
    const auto *mg_378 = buffer.data(mg + 378);
    const auto *mg_379 = buffer.data(mg + 379);
    const auto *mg_380 = buffer.data(mg + 380);
    const auto *mg_381 = buffer.data(mg + 381);
    const auto *mg_382 = buffer.data(mg + 382);
    const auto *mg_383 = buffer.data(mg + 383);
    const auto *mg_384 = buffer.data(mg + 384);
    const auto *mg_385 = buffer.data(mg + 385);
    const auto *mg_386 = buffer.data(mg + 386);
    const auto *mg_387 = buffer.data(mg + 387);
    const auto *mg_388 = buffer.data(mg + 388);
    const auto *mg_389 = buffer.data(mg + 389);
    const auto *mg_390 = buffer.data(mg + 390);
    const auto *mg_391 = buffer.data(mg + 391);
    const auto *mg_392 = buffer.data(mg + 392);
    const auto *mg_393 = buffer.data(mg + 393);
    const auto *mg_394 = buffer.data(mg + 394);
    const auto *mg_395 = buffer.data(mg + 395);
    const auto *mg_396 = buffer.data(mg + 396);
    const auto *mg_397 = buffer.data(mg + 397);
    const auto *mg_398 = buffer.data(mg + 398);
    const auto *mg_399 = buffer.data(mg + 399);
    const auto *mg_400 = buffer.data(mg + 400);
    const auto *mg_401 = buffer.data(mg + 401);
    const auto *mg_402 = buffer.data(mg + 402);
    const auto *mg_403 = buffer.data(mg + 403);
    const auto *mg_404 = buffer.data(mg + 404);
    const auto *mg_405 = buffer.data(mg + 405);
    const auto *mg_406 = buffer.data(mg + 406);
    const auto *mg_407 = buffer.data(mg + 407);
    const auto *mg_408 = buffer.data(mg + 408);
    const auto *mg_409 = buffer.data(mg + 409);
    const auto *mg_410 = buffer.data(mg + 410);
    const auto *mg_411 = buffer.data(mg + 411);
    const auto *mg_412 = buffer.data(mg + 412);
    const auto *mg_413 = buffer.data(mg + 413);
    const auto *mg_414 = buffer.data(mg + 414);
    const auto *mg_415 = buffer.data(mg + 415);
    const auto *mg_416 = buffer.data(mg + 416);
    const auto *mg_417 = buffer.data(mg + 417);
    const auto *mg_418 = buffer.data(mg + 418);
    const auto *mg_419 = buffer.data(mg + 419);
    const auto *mg_420 = buffer.data(mg + 420);
    const auto *mg_421 = buffer.data(mg + 421);
    const auto *mg_422 = buffer.data(mg + 422);
    const auto *mg_423 = buffer.data(mg + 423);
    const auto *mg_424 = buffer.data(mg + 424);
    const auto *mg_425 = buffer.data(mg + 425);
    const auto *mg_426 = buffer.data(mg + 426);
    const auto *mg_427 = buffer.data(mg + 427);
    const auto *mg_428 = buffer.data(mg + 428);
    const auto *mg_429 = buffer.data(mg + 429);
    const auto *mg_430 = buffer.data(mg + 430);
    const auto *mg_431 = buffer.data(mg + 431);
    const auto *mg_432 = buffer.data(mg + 432);
    const auto *mg_433 = buffer.data(mg + 433);
    const auto *mg_434 = buffer.data(mg + 434);
    const auto *mg_435 = buffer.data(mg + 435);
    const auto *mg_436 = buffer.data(mg + 436);
    const auto *mg_437 = buffer.data(mg + 437);
    const auto *mg_438 = buffer.data(mg + 438);
    const auto *mg_439 = buffer.data(mg + 439);
    const auto *mg_440 = buffer.data(mg + 440);
    const auto *mg_441 = buffer.data(mg + 441);
    const auto *mg_442 = buffer.data(mg + 442);
    const auto *mg_443 = buffer.data(mg + 443);
    const auto *mg_444 = buffer.data(mg + 444);
    const auto *mg_445 = buffer.data(mg + 445);
    const auto *mg_446 = buffer.data(mg + 446);
    const auto *mg_447 = buffer.data(mg + 447);
    const auto *mg_448 = buffer.data(mg + 448);
    const auto *mg_449 = buffer.data(mg + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, kg_300, kg_301, kg_302, kg_303, \
                         kg_304, mg_300, mg_301, mg_302, mg_303, \
                         mg_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -3.0 * kg_300[k]
                   + f_0 * mg_300[k];

        t_301[k] = -3.0 * kg_301[k]
                   + f_0 * mg_301[k];

        t_302[k] = -3.0 * kg_302[k]
                   + f_0 * mg_302[k];

        t_303[k] = -3.0 * kg_303[k]
                   + f_0 * mg_303[k];

        t_304[k] = -3.0 * kg_304[k]
                   + f_0 * mg_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, kg_305, kg_306, kg_307, kg_308, \
                         kg_309, mg_305, mg_306, mg_307, mg_308, \
                         mg_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -3.0 * kg_305[k]
                   + f_0 * mg_305[k];

        t_306[k] = -3.0 * kg_306[k]
                   + f_0 * mg_306[k];

        t_307[k] = -3.0 * kg_307[k]
                   + f_0 * mg_307[k];

        t_308[k] = -3.0 * kg_308[k]
                   + f_0 * mg_308[k];

        t_309[k] = -3.0 * kg_309[k]
                   + f_0 * mg_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, kg_310, kg_311, kg_312, kg_313, \
                         kg_314, mg_310, mg_311, mg_312, mg_313, \
                         mg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * kg_310[k]
                   + f_0 * mg_310[k];

        t_311[k] = -3.0 * kg_311[k]
                   + f_0 * mg_311[k];

        t_312[k] = -3.0 * kg_312[k]
                   + f_0 * mg_312[k];

        t_313[k] = -3.0 * kg_313[k]
                   + f_0 * mg_313[k];

        t_314[k] = -3.0 * kg_314[k]
                   + f_0 * mg_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, kg_315, kg_316, kg_317, kg_318, \
                         kg_319, mg_315, mg_316, mg_317, mg_318, \
                         mg_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -2.0 * kg_315[k]
                   + f_0 * mg_315[k];

        t_316[k] = -2.0 * kg_316[k]
                   + f_0 * mg_316[k];

        t_317[k] = -2.0 * kg_317[k]
                   + f_0 * mg_317[k];

        t_318[k] = -2.0 * kg_318[k]
                   + f_0 * mg_318[k];

        t_319[k] = -2.0 * kg_319[k]
                   + f_0 * mg_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, kg_320, kg_321, kg_322, kg_323, \
                         kg_324, mg_320, mg_321, mg_322, mg_323, \
                         mg_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -2.0 * kg_320[k]
                   + f_0 * mg_320[k];

        t_321[k] = -2.0 * kg_321[k]
                   + f_0 * mg_321[k];

        t_322[k] = -2.0 * kg_322[k]
                   + f_0 * mg_322[k];

        t_323[k] = -2.0 * kg_323[k]
                   + f_0 * mg_323[k];

        t_324[k] = -2.0 * kg_324[k]
                   + f_0 * mg_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, kg_325, kg_326, kg_327, kg_328, \
                         kg_329, mg_325, mg_326, mg_327, mg_328, \
                         mg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -2.0 * kg_325[k]
                   + f_0 * mg_325[k];

        t_326[k] = -2.0 * kg_326[k]
                   + f_0 * mg_326[k];

        t_327[k] = -2.0 * kg_327[k]
                   + f_0 * mg_327[k];

        t_328[k] = -2.0 * kg_328[k]
                   + f_0 * mg_328[k];

        t_329[k] = -2.0 * kg_329[k]
                   + f_0 * mg_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, kg_330, kg_331, kg_332, kg_333, \
                         kg_334, mg_330, mg_331, mg_332, mg_333, \
                         mg_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -2.0 * kg_330[k]
                   + f_0 * mg_330[k];

        t_331[k] = -2.0 * kg_331[k]
                   + f_0 * mg_331[k];

        t_332[k] = -2.0 * kg_332[k]
                   + f_0 * mg_332[k];

        t_333[k] = -2.0 * kg_333[k]
                   + f_0 * mg_333[k];

        t_334[k] = -2.0 * kg_334[k]
                   + f_0 * mg_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, kg_335, kg_336, kg_337, kg_338, \
                         kg_339, mg_335, mg_336, mg_337, mg_338, \
                         mg_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -2.0 * kg_335[k]
                   + f_0 * mg_335[k];

        t_336[k] = -2.0 * kg_336[k]
                   + f_0 * mg_336[k];

        t_337[k] = -2.0 * kg_337[k]
                   + f_0 * mg_337[k];

        t_338[k] = -2.0 * kg_338[k]
                   + f_0 * mg_338[k];

        t_339[k] = -2.0 * kg_339[k]
                   + f_0 * mg_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, kg_340, kg_341, kg_342, kg_343, \
                         kg_344, mg_340, mg_341, mg_342, mg_343, \
                         mg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * kg_340[k]
                   + f_0 * mg_340[k];

        t_341[k] = -2.0 * kg_341[k]
                   + f_0 * mg_341[k];

        t_342[k] = -2.0 * kg_342[k]
                   + f_0 * mg_342[k];

        t_343[k] = -2.0 * kg_343[k]
                   + f_0 * mg_343[k];

        t_344[k] = -2.0 * kg_344[k]
                   + f_0 * mg_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, kg_345, kg_346, kg_347, kg_348, \
                         kg_349, mg_345, mg_346, mg_347, mg_348, \
                         mg_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * kg_345[k]
                   + f_0 * mg_345[k];

        t_346[k] = -2.0 * kg_346[k]
                   + f_0 * mg_346[k];

        t_347[k] = -2.0 * kg_347[k]
                   + f_0 * mg_347[k];

        t_348[k] = -2.0 * kg_348[k]
                   + f_0 * mg_348[k];

        t_349[k] = -2.0 * kg_349[k]
                   + f_0 * mg_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, kg_350, kg_351, kg_352, kg_353, \
                         kg_354, mg_350, mg_351, mg_352, mg_353, \
                         mg_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * kg_350[k]
                   + f_0 * mg_350[k];

        t_351[k] = -2.0 * kg_351[k]
                   + f_0 * mg_351[k];

        t_352[k] = -2.0 * kg_352[k]
                   + f_0 * mg_352[k];

        t_353[k] = -2.0 * kg_353[k]
                   + f_0 * mg_353[k];

        t_354[k] = -2.0 * kg_354[k]
                   + f_0 * mg_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, kg_355, kg_356, kg_357, kg_358, \
                         kg_359, mg_355, mg_356, mg_357, mg_358, \
                         mg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * kg_355[k]
                   + f_0 * mg_355[k];

        t_356[k] = -2.0 * kg_356[k]
                   + f_0 * mg_356[k];

        t_357[k] = -2.0 * kg_357[k]
                   + f_0 * mg_357[k];

        t_358[k] = -2.0 * kg_358[k]
                   + f_0 * mg_358[k];

        t_359[k] = -2.0 * kg_359[k]
                   + f_0 * mg_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, kg_360, kg_361, kg_362, kg_363, \
                         kg_364, mg_360, mg_361, mg_362, mg_363, \
                         mg_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -2.0 * kg_360[k]
                   + f_0 * mg_360[k];

        t_361[k] = -2.0 * kg_361[k]
                   + f_0 * mg_361[k];

        t_362[k] = -2.0 * kg_362[k]
                   + f_0 * mg_362[k];

        t_363[k] = -2.0 * kg_363[k]
                   + f_0 * mg_363[k];

        t_364[k] = -2.0 * kg_364[k]
                   + f_0 * mg_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, kg_365, kg_366, kg_367, kg_368, \
                         kg_369, mg_365, mg_366, mg_367, mg_368, \
                         mg_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -2.0 * kg_365[k]
                   + f_0 * mg_365[k];

        t_366[k] = -2.0 * kg_366[k]
                   + f_0 * mg_366[k];

        t_367[k] = -2.0 * kg_367[k]
                   + f_0 * mg_367[k];

        t_368[k] = -2.0 * kg_368[k]
                   + f_0 * mg_368[k];

        t_369[k] = -2.0 * kg_369[k]
                   + f_0 * mg_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, kg_370, kg_371, kg_372, kg_373, \
                         kg_374, mg_370, mg_371, mg_372, mg_373, \
                         mg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -2.0 * kg_370[k]
                   + f_0 * mg_370[k];

        t_371[k] = -2.0 * kg_371[k]
                   + f_0 * mg_371[k];

        t_372[k] = -2.0 * kg_372[k]
                   + f_0 * mg_372[k];

        t_373[k] = -2.0 * kg_373[k]
                   + f_0 * mg_373[k];

        t_374[k] = -2.0 * kg_374[k]
                   + f_0 * mg_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, kg_375, kg_376, kg_377, kg_378, \
                         kg_379, mg_375, mg_376, mg_377, mg_378, \
                         mg_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -2.0 * kg_375[k]
                   + f_0 * mg_375[k];

        t_376[k] = -2.0 * kg_376[k]
                   + f_0 * mg_376[k];

        t_377[k] = -2.0 * kg_377[k]
                   + f_0 * mg_377[k];

        t_378[k] = -2.0 * kg_378[k]
                   + f_0 * mg_378[k];

        t_379[k] = -2.0 * kg_379[k]
                   + f_0 * mg_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, kg_380, kg_381, kg_382, kg_383, \
                         kg_384, mg_380, mg_381, mg_382, mg_383, \
                         mg_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * kg_380[k]
                   + f_0 * mg_380[k];

        t_381[k] = -2.0 * kg_381[k]
                   + f_0 * mg_381[k];

        t_382[k] = -2.0 * kg_382[k]
                   + f_0 * mg_382[k];

        t_383[k] = -2.0 * kg_383[k]
                   + f_0 * mg_383[k];

        t_384[k] = -2.0 * kg_384[k]
                   + f_0 * mg_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, kg_385, kg_386, kg_387, kg_388, \
                         kg_389, mg_385, mg_386, mg_387, mg_388, \
                         mg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * kg_385[k]
                   + f_0 * mg_385[k];

        t_386[k] = -2.0 * kg_386[k]
                   + f_0 * mg_386[k];

        t_387[k] = -2.0 * kg_387[k]
                   + f_0 * mg_387[k];

        t_388[k] = -2.0 * kg_388[k]
                   + f_0 * mg_388[k];

        t_389[k] = -2.0 * kg_389[k]
                   + f_0 * mg_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, kg_390, kg_391, kg_392, kg_393, \
                         kg_394, mg_390, mg_391, mg_392, mg_393, \
                         mg_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -2.0 * kg_390[k]
                   + f_0 * mg_390[k];

        t_391[k] = -2.0 * kg_391[k]
                   + f_0 * mg_391[k];

        t_392[k] = -2.0 * kg_392[k]
                   + f_0 * mg_392[k];

        t_393[k] = -2.0 * kg_393[k]
                   + f_0 * mg_393[k];

        t_394[k] = -2.0 * kg_394[k]
                   + f_0 * mg_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, kg_395, kg_396, kg_397, kg_398, \
                         kg_399, mg_395, mg_396, mg_397, mg_398, \
                         mg_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -2.0 * kg_395[k]
                   + f_0 * mg_395[k];

        t_396[k] = -2.0 * kg_396[k]
                   + f_0 * mg_396[k];

        t_397[k] = -2.0 * kg_397[k]
                   + f_0 * mg_397[k];

        t_398[k] = -2.0 * kg_398[k]
                   + f_0 * mg_398[k];

        t_399[k] = -2.0 * kg_399[k]
                   + f_0 * mg_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, kg_400, kg_401, kg_402, kg_403, \
                         kg_404, mg_400, mg_401, mg_402, mg_403, \
                         mg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -2.0 * kg_400[k]
                   + f_0 * mg_400[k];

        t_401[k] = -2.0 * kg_401[k]
                   + f_0 * mg_401[k];

        t_402[k] = -2.0 * kg_402[k]
                   + f_0 * mg_402[k];

        t_403[k] = -2.0 * kg_403[k]
                   + f_0 * mg_403[k];

        t_404[k] = -2.0 * kg_404[k]
                   + f_0 * mg_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, kg_405, kg_406, kg_407, kg_408, \
                         kg_409, mg_405, mg_406, mg_407, mg_408, \
                         mg_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -2.0 * kg_405[k]
                   + f_0 * mg_405[k];

        t_406[k] = -2.0 * kg_406[k]
                   + f_0 * mg_406[k];

        t_407[k] = -2.0 * kg_407[k]
                   + f_0 * mg_407[k];

        t_408[k] = -2.0 * kg_408[k]
                   + f_0 * mg_408[k];

        t_409[k] = -2.0 * kg_409[k]
                   + f_0 * mg_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, kg_410, kg_411, kg_412, kg_413, \
                         kg_414, mg_410, mg_411, mg_412, mg_413, \
                         mg_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -2.0 * kg_410[k]
                   + f_0 * mg_410[k];

        t_411[k] = -2.0 * kg_411[k]
                   + f_0 * mg_411[k];

        t_412[k] = -2.0 * kg_412[k]
                   + f_0 * mg_412[k];

        t_413[k] = -2.0 * kg_413[k]
                   + f_0 * mg_413[k];

        t_414[k] = -2.0 * kg_414[k]
                   + f_0 * mg_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, kg_415, kg_416, kg_417, kg_418, \
                         kg_419, mg_415, mg_416, mg_417, mg_418, \
                         mg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -2.0 * kg_415[k]
                   + f_0 * mg_415[k];

        t_416[k] = -2.0 * kg_416[k]
                   + f_0 * mg_416[k];

        t_417[k] = -2.0 * kg_417[k]
                   + f_0 * mg_417[k];

        t_418[k] = -2.0 * kg_418[k]
                   + f_0 * mg_418[k];

        t_419[k] = -2.0 * kg_419[k]
                   + f_0 * mg_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, kg_420, kg_421, kg_422, kg_423, \
                         kg_424, mg_420, mg_421, mg_422, mg_423, \
                         mg_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -kg_420[k]
                   + f_0 * mg_420[k];

        t_421[k] = -kg_421[k]
                   + f_0 * mg_421[k];

        t_422[k] = -kg_422[k]
                   + f_0 * mg_422[k];

        t_423[k] = -kg_423[k]
                   + f_0 * mg_423[k];

        t_424[k] = -kg_424[k]
                   + f_0 * mg_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, kg_425, kg_426, kg_427, kg_428, \
                         kg_429, mg_425, mg_426, mg_427, mg_428, \
                         mg_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -kg_425[k]
                   + f_0 * mg_425[k];

        t_426[k] = -kg_426[k]
                   + f_0 * mg_426[k];

        t_427[k] = -kg_427[k]
                   + f_0 * mg_427[k];

        t_428[k] = -kg_428[k]
                   + f_0 * mg_428[k];

        t_429[k] = -kg_429[k]
                   + f_0 * mg_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, kg_430, kg_431, kg_432, kg_433, \
                         kg_434, mg_430, mg_431, mg_432, mg_433, \
                         mg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -kg_430[k]
                   + f_0 * mg_430[k];

        t_431[k] = -kg_431[k]
                   + f_0 * mg_431[k];

        t_432[k] = -kg_432[k]
                   + f_0 * mg_432[k];

        t_433[k] = -kg_433[k]
                   + f_0 * mg_433[k];

        t_434[k] = -kg_434[k]
                   + f_0 * mg_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, kg_435, kg_436, kg_437, kg_438, \
                         kg_439, mg_435, mg_436, mg_437, mg_438, \
                         mg_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -kg_435[k]
                   + f_0 * mg_435[k];

        t_436[k] = -kg_436[k]
                   + f_0 * mg_436[k];

        t_437[k] = -kg_437[k]
                   + f_0 * mg_437[k];

        t_438[k] = -kg_438[k]
                   + f_0 * mg_438[k];

        t_439[k] = -kg_439[k]
                   + f_0 * mg_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, kg_440, kg_441, kg_442, kg_443, \
                         kg_444, mg_440, mg_441, mg_442, mg_443, \
                         mg_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -kg_440[k]
                   + f_0 * mg_440[k];

        t_441[k] = -kg_441[k]
                   + f_0 * mg_441[k];

        t_442[k] = -kg_442[k]
                   + f_0 * mg_442[k];

        t_443[k] = -kg_443[k]
                   + f_0 * mg_443[k];

        t_444[k] = -kg_444[k]
                   + f_0 * mg_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, kg_445, kg_446, kg_447, kg_448, \
                         kg_449, mg_445, mg_446, mg_447, mg_448, \
                         mg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -kg_445[k]
                   + f_0 * mg_445[k];

        t_446[k] = -kg_446[k]
                   + f_0 * mg_446[k];

        t_447[k] = -kg_447[k]
                   + f_0 * mg_447[k];

        t_448[k] = -kg_448[k]
                   + f_0 * mg_448[k];

        t_449[k] = -kg_449[k]
                   + f_0 * mg_449[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
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

    const auto *mg_450 = buffer.data(mg + 450);
    const auto *mg_451 = buffer.data(mg + 451);
    const auto *mg_452 = buffer.data(mg + 452);
    const auto *mg_453 = buffer.data(mg + 453);
    const auto *mg_454 = buffer.data(mg + 454);
    const auto *mg_455 = buffer.data(mg + 455);
    const auto *mg_456 = buffer.data(mg + 456);
    const auto *mg_457 = buffer.data(mg + 457);
    const auto *mg_458 = buffer.data(mg + 458);
    const auto *mg_459 = buffer.data(mg + 459);
    const auto *mg_460 = buffer.data(mg + 460);
    const auto *mg_461 = buffer.data(mg + 461);
    const auto *mg_462 = buffer.data(mg + 462);
    const auto *mg_463 = buffer.data(mg + 463);
    const auto *mg_464 = buffer.data(mg + 464);
    const auto *mg_465 = buffer.data(mg + 465);
    const auto *mg_466 = buffer.data(mg + 466);
    const auto *mg_467 = buffer.data(mg + 467);
    const auto *mg_468 = buffer.data(mg + 468);
    const auto *mg_469 = buffer.data(mg + 469);
    const auto *mg_470 = buffer.data(mg + 470);
    const auto *mg_471 = buffer.data(mg + 471);
    const auto *mg_472 = buffer.data(mg + 472);
    const auto *mg_473 = buffer.data(mg + 473);
    const auto *mg_474 = buffer.data(mg + 474);
    const auto *mg_475 = buffer.data(mg + 475);
    const auto *mg_476 = buffer.data(mg + 476);
    const auto *mg_477 = buffer.data(mg + 477);
    const auto *mg_478 = buffer.data(mg + 478);
    const auto *mg_479 = buffer.data(mg + 479);
    const auto *mg_480 = buffer.data(mg + 480);
    const auto *mg_481 = buffer.data(mg + 481);
    const auto *mg_482 = buffer.data(mg + 482);
    const auto *mg_483 = buffer.data(mg + 483);
    const auto *mg_484 = buffer.data(mg + 484);
    const auto *mg_485 = buffer.data(mg + 485);
    const auto *mg_486 = buffer.data(mg + 486);
    const auto *mg_487 = buffer.data(mg + 487);
    const auto *mg_488 = buffer.data(mg + 488);
    const auto *mg_489 = buffer.data(mg + 489);
    const auto *mg_490 = buffer.data(mg + 490);
    const auto *mg_491 = buffer.data(mg + 491);
    const auto *mg_492 = buffer.data(mg + 492);
    const auto *mg_493 = buffer.data(mg + 493);
    const auto *mg_494 = buffer.data(mg + 494);
    const auto *mg_495 = buffer.data(mg + 495);
    const auto *mg_496 = buffer.data(mg + 496);
    const auto *mg_497 = buffer.data(mg + 497);
    const auto *mg_498 = buffer.data(mg + 498);
    const auto *mg_499 = buffer.data(mg + 499);
    const auto *mg_500 = buffer.data(mg + 500);
    const auto *mg_501 = buffer.data(mg + 501);
    const auto *mg_502 = buffer.data(mg + 502);
    const auto *mg_503 = buffer.data(mg + 503);
    const auto *mg_504 = buffer.data(mg + 504);
    const auto *mg_505 = buffer.data(mg + 505);
    const auto *mg_506 = buffer.data(mg + 506);
    const auto *mg_507 = buffer.data(mg + 507);
    const auto *mg_508 = buffer.data(mg + 508);
    const auto *mg_509 = buffer.data(mg + 509);
    const auto *mg_510 = buffer.data(mg + 510);
    const auto *mg_511 = buffer.data(mg + 511);
    const auto *mg_512 = buffer.data(mg + 512);
    const auto *mg_513 = buffer.data(mg + 513);
    const auto *mg_514 = buffer.data(mg + 514);
    const auto *mg_515 = buffer.data(mg + 515);
    const auto *mg_516 = buffer.data(mg + 516);
    const auto *mg_517 = buffer.data(mg + 517);
    const auto *mg_518 = buffer.data(mg + 518);
    const auto *mg_519 = buffer.data(mg + 519);
    const auto *mg_520 = buffer.data(mg + 520);
    const auto *mg_521 = buffer.data(mg + 521);
    const auto *mg_522 = buffer.data(mg + 522);
    const auto *mg_523 = buffer.data(mg + 523);
    const auto *mg_524 = buffer.data(mg + 524);
    const auto *mg_525 = buffer.data(mg + 525);
    const auto *mg_526 = buffer.data(mg + 526);
    const auto *mg_527 = buffer.data(mg + 527);
    const auto *mg_528 = buffer.data(mg + 528);
    const auto *mg_529 = buffer.data(mg + 529);
    const auto *mg_530 = buffer.data(mg + 530);
    const auto *mg_531 = buffer.data(mg + 531);
    const auto *mg_532 = buffer.data(mg + 532);
    const auto *mg_533 = buffer.data(mg + 533);
    const auto *mg_534 = buffer.data(mg + 534);
    const auto *mg_535 = buffer.data(mg + 535);
    const auto *mg_536 = buffer.data(mg + 536);
    const auto *mg_537 = buffer.data(mg + 537);
    const auto *mg_538 = buffer.data(mg + 538);
    const auto *mg_539 = buffer.data(mg + 539);
    const auto *mg_540 = buffer.data(mg + 540);
    const auto *mg_541 = buffer.data(mg + 541);
    const auto *mg_542 = buffer.data(mg + 542);
    const auto *mg_543 = buffer.data(mg + 543);
    const auto *mg_544 = buffer.data(mg + 544);
    const auto *mg_545 = buffer.data(mg + 545);
    const auto *mg_546 = buffer.data(mg + 546);
    const auto *mg_547 = buffer.data(mg + 547);
    const auto *mg_548 = buffer.data(mg + 548);
    const auto *mg_549 = buffer.data(mg + 549);
    const auto *mg_550 = buffer.data(mg + 550);
    const auto *mg_551 = buffer.data(mg + 551);
    const auto *mg_552 = buffer.data(mg + 552);
    const auto *mg_553 = buffer.data(mg + 553);
    const auto *mg_554 = buffer.data(mg + 554);
    const auto *mg_555 = buffer.data(mg + 555);
    const auto *mg_556 = buffer.data(mg + 556);
    const auto *mg_557 = buffer.data(mg + 557);
    const auto *mg_558 = buffer.data(mg + 558);
    const auto *mg_559 = buffer.data(mg + 559);
    const auto *mg_560 = buffer.data(mg + 560);
    const auto *mg_561 = buffer.data(mg + 561);
    const auto *mg_562 = buffer.data(mg + 562);
    const auto *mg_563 = buffer.data(mg + 563);
    const auto *mg_564 = buffer.data(mg + 564);
    const auto *mg_565 = buffer.data(mg + 565);
    const auto *mg_566 = buffer.data(mg + 566);
    const auto *mg_567 = buffer.data(mg + 567);
    const auto *mg_568 = buffer.data(mg + 568);
    const auto *mg_569 = buffer.data(mg + 569);
    const auto *mg_570 = buffer.data(mg + 570);
    const auto *mg_571 = buffer.data(mg + 571);
    const auto *mg_572 = buffer.data(mg + 572);
    const auto *mg_573 = buffer.data(mg + 573);
    const auto *mg_574 = buffer.data(mg + 574);
    const auto *mg_575 = buffer.data(mg + 575);
    const auto *mg_576 = buffer.data(mg + 576);
    const auto *mg_577 = buffer.data(mg + 577);
    const auto *mg_578 = buffer.data(mg + 578);
    const auto *mg_579 = buffer.data(mg + 579);
    const auto *mg_580 = buffer.data(mg + 580);
    const auto *mg_581 = buffer.data(mg + 581);
    const auto *mg_582 = buffer.data(mg + 582);
    const auto *mg_583 = buffer.data(mg + 583);
    const auto *mg_584 = buffer.data(mg + 584);
    const auto *mg_585 = buffer.data(mg + 585);
    const auto *mg_586 = buffer.data(mg + 586);
    const auto *mg_587 = buffer.data(mg + 587);
    const auto *mg_588 = buffer.data(mg + 588);
    const auto *mg_589 = buffer.data(mg + 589);
    const auto *mg_590 = buffer.data(mg + 590);
    const auto *mg_591 = buffer.data(mg + 591);
    const auto *mg_592 = buffer.data(mg + 592);
    const auto *mg_593 = buffer.data(mg + 593);
    const auto *mg_594 = buffer.data(mg + 594);
    const auto *mg_595 = buffer.data(mg + 595);
    const auto *mg_596 = buffer.data(mg + 596);
    const auto *mg_597 = buffer.data(mg + 597);
    const auto *mg_598 = buffer.data(mg + 598);
    const auto *mg_599 = buffer.data(mg + 599);
    const auto *mg_600 = buffer.data(mg + 600);
    const auto *mg_601 = buffer.data(mg + 601);
    const auto *mg_602 = buffer.data(mg + 602);
    const auto *mg_603 = buffer.data(mg + 603);
    const auto *mg_604 = buffer.data(mg + 604);
    const auto *mg_605 = buffer.data(mg + 605);
    const auto *mg_606 = buffer.data(mg + 606);
    const auto *mg_607 = buffer.data(mg + 607);
    const auto *mg_608 = buffer.data(mg + 608);
    const auto *mg_609 = buffer.data(mg + 609);
    const auto *mg_610 = buffer.data(mg + 610);
    const auto *mg_611 = buffer.data(mg + 611);
    const auto *mg_612 = buffer.data(mg + 612);
    const auto *mg_613 = buffer.data(mg + 613);
    const auto *mg_614 = buffer.data(mg + 614);
    const auto *mg_615 = buffer.data(mg + 615);
    const auto *mg_616 = buffer.data(mg + 616);
    const auto *mg_617 = buffer.data(mg + 617);
    const auto *mg_618 = buffer.data(mg + 618);
    const auto *mg_619 = buffer.data(mg + 619);
    const auto *mg_620 = buffer.data(mg + 620);
    const auto *mg_621 = buffer.data(mg + 621);
    const auto *mg_622 = buffer.data(mg + 622);
    const auto *mg_623 = buffer.data(mg + 623);
    const auto *mg_624 = buffer.data(mg + 624);
    const auto *mg_625 = buffer.data(mg + 625);
    const auto *mg_626 = buffer.data(mg + 626);
    const auto *mg_627 = buffer.data(mg + 627);

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, kg_450, kg_451, kg_452, kg_453, \
                         kg_454, mg_450, mg_451, mg_452, mg_453, \
                         mg_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -kg_450[k]
                   + f_0 * mg_450[k];

        t_451[k] = -kg_451[k]
                   + f_0 * mg_451[k];

        t_452[k] = -kg_452[k]
                   + f_0 * mg_452[k];

        t_453[k] = -kg_453[k]
                   + f_0 * mg_453[k];

        t_454[k] = -kg_454[k]
                   + f_0 * mg_454[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, kg_455, kg_456, kg_457, kg_458, \
                         kg_459, mg_455, mg_456, mg_457, mg_458, \
                         mg_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -kg_455[k]
                   + f_0 * mg_455[k];

        t_456[k] = -kg_456[k]
                   + f_0 * mg_456[k];

        t_457[k] = -kg_457[k]
                   + f_0 * mg_457[k];

        t_458[k] = -kg_458[k]
                   + f_0 * mg_458[k];

        t_459[k] = -kg_459[k]
                   + f_0 * mg_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, kg_460, kg_461, kg_462, kg_463, \
                         kg_464, mg_460, mg_461, mg_462, mg_463, \
                         mg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -kg_460[k]
                   + f_0 * mg_460[k];

        t_461[k] = -kg_461[k]
                   + f_0 * mg_461[k];

        t_462[k] = -kg_462[k]
                   + f_0 * mg_462[k];

        t_463[k] = -kg_463[k]
                   + f_0 * mg_463[k];

        t_464[k] = -kg_464[k]
                   + f_0 * mg_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, kg_465, kg_466, kg_467, kg_468, \
                         kg_469, mg_465, mg_466, mg_467, mg_468, \
                         mg_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -kg_465[k]
                   + f_0 * mg_465[k];

        t_466[k] = -kg_466[k]
                   + f_0 * mg_466[k];

        t_467[k] = -kg_467[k]
                   + f_0 * mg_467[k];

        t_468[k] = -kg_468[k]
                   + f_0 * mg_468[k];

        t_469[k] = -kg_469[k]
                   + f_0 * mg_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, kg_470, kg_471, kg_472, kg_473, \
                         kg_474, mg_470, mg_471, mg_472, mg_473, \
                         mg_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -kg_470[k]
                   + f_0 * mg_470[k];

        t_471[k] = -kg_471[k]
                   + f_0 * mg_471[k];

        t_472[k] = -kg_472[k]
                   + f_0 * mg_472[k];

        t_473[k] = -kg_473[k]
                   + f_0 * mg_473[k];

        t_474[k] = -kg_474[k]
                   + f_0 * mg_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, kg_475, kg_476, kg_477, kg_478, \
                         kg_479, mg_475, mg_476, mg_477, mg_478, \
                         mg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -kg_475[k]
                   + f_0 * mg_475[k];

        t_476[k] = -kg_476[k]
                   + f_0 * mg_476[k];

        t_477[k] = -kg_477[k]
                   + f_0 * mg_477[k];

        t_478[k] = -kg_478[k]
                   + f_0 * mg_478[k];

        t_479[k] = -kg_479[k]
                   + f_0 * mg_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, kg_480, kg_481, kg_482, kg_483, \
                         kg_484, mg_480, mg_481, mg_482, mg_483, \
                         mg_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -kg_480[k]
                   + f_0 * mg_480[k];

        t_481[k] = -kg_481[k]
                   + f_0 * mg_481[k];

        t_482[k] = -kg_482[k]
                   + f_0 * mg_482[k];

        t_483[k] = -kg_483[k]
                   + f_0 * mg_483[k];

        t_484[k] = -kg_484[k]
                   + f_0 * mg_484[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, kg_485, kg_486, kg_487, kg_488, \
                         kg_489, mg_485, mg_486, mg_487, mg_488, \
                         mg_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -kg_485[k]
                   + f_0 * mg_485[k];

        t_486[k] = -kg_486[k]
                   + f_0 * mg_486[k];

        t_487[k] = -kg_487[k]
                   + f_0 * mg_487[k];

        t_488[k] = -kg_488[k]
                   + f_0 * mg_488[k];

        t_489[k] = -kg_489[k]
                   + f_0 * mg_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, kg_490, kg_491, kg_492, kg_493, \
                         kg_494, mg_490, mg_491, mg_492, mg_493, \
                         mg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -kg_490[k]
                   + f_0 * mg_490[k];

        t_491[k] = -kg_491[k]
                   + f_0 * mg_491[k];

        t_492[k] = -kg_492[k]
                   + f_0 * mg_492[k];

        t_493[k] = -kg_493[k]
                   + f_0 * mg_493[k];

        t_494[k] = -kg_494[k]
                   + f_0 * mg_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, kg_495, kg_496, kg_497, kg_498, \
                         kg_499, mg_495, mg_496, mg_497, mg_498, \
                         mg_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -kg_495[k]
                   + f_0 * mg_495[k];

        t_496[k] = -kg_496[k]
                   + f_0 * mg_496[k];

        t_497[k] = -kg_497[k]
                   + f_0 * mg_497[k];

        t_498[k] = -kg_498[k]
                   + f_0 * mg_498[k];

        t_499[k] = -kg_499[k]
                   + f_0 * mg_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, kg_500, kg_501, kg_502, kg_503, \
                         kg_504, mg_500, mg_501, mg_502, mg_503, \
                         mg_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -kg_500[k]
                   + f_0 * mg_500[k];

        t_501[k] = -kg_501[k]
                   + f_0 * mg_501[k];

        t_502[k] = -kg_502[k]
                   + f_0 * mg_502[k];

        t_503[k] = -kg_503[k]
                   + f_0 * mg_503[k];

        t_504[k] = -kg_504[k]
                   + f_0 * mg_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, kg_505, kg_506, kg_507, kg_508, \
                         kg_509, mg_505, mg_506, mg_507, mg_508, \
                         mg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -kg_505[k]
                   + f_0 * mg_505[k];

        t_506[k] = -kg_506[k]
                   + f_0 * mg_506[k];

        t_507[k] = -kg_507[k]
                   + f_0 * mg_507[k];

        t_508[k] = -kg_508[k]
                   + f_0 * mg_508[k];

        t_509[k] = -kg_509[k]
                   + f_0 * mg_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, kg_510, kg_511, kg_512, kg_513, \
                         kg_514, mg_510, mg_511, mg_512, mg_513, \
                         mg_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -kg_510[k]
                   + f_0 * mg_510[k];

        t_511[k] = -kg_511[k]
                   + f_0 * mg_511[k];

        t_512[k] = -kg_512[k]
                   + f_0 * mg_512[k];

        t_513[k] = -kg_513[k]
                   + f_0 * mg_513[k];

        t_514[k] = -kg_514[k]
                   + f_0 * mg_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, kg_515, kg_516, kg_517, kg_518, \
                         kg_519, mg_515, mg_516, mg_517, mg_518, \
                         mg_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -kg_515[k]
                   + f_0 * mg_515[k];

        t_516[k] = -kg_516[k]
                   + f_0 * mg_516[k];

        t_517[k] = -kg_517[k]
                   + f_0 * mg_517[k];

        t_518[k] = -kg_518[k]
                   + f_0 * mg_518[k];

        t_519[k] = -kg_519[k]
                   + f_0 * mg_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, kg_520, kg_521, kg_522, kg_523, \
                         kg_524, mg_520, mg_521, mg_522, mg_523, \
                         mg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -kg_520[k]
                   + f_0 * mg_520[k];

        t_521[k] = -kg_521[k]
                   + f_0 * mg_521[k];

        t_522[k] = -kg_522[k]
                   + f_0 * mg_522[k];

        t_523[k] = -kg_523[k]
                   + f_0 * mg_523[k];

        t_524[k] = -kg_524[k]
                   + f_0 * mg_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, kg_525, kg_526, kg_527, kg_528, \
                         kg_529, mg_525, mg_526, mg_527, mg_528, \
                         mg_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -kg_525[k]
                   + f_0 * mg_525[k];

        t_526[k] = -kg_526[k]
                   + f_0 * mg_526[k];

        t_527[k] = -kg_527[k]
                   + f_0 * mg_527[k];

        t_528[k] = -kg_528[k]
                   + f_0 * mg_528[k];

        t_529[k] = -kg_529[k]
                   + f_0 * mg_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, kg_530, kg_531, kg_532, kg_533, \
                         kg_534, mg_530, mg_531, mg_532, mg_533, \
                         mg_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -kg_530[k]
                   + f_0 * mg_530[k];

        t_531[k] = -kg_531[k]
                   + f_0 * mg_531[k];

        t_532[k] = -kg_532[k]
                   + f_0 * mg_532[k];

        t_533[k] = -kg_533[k]
                   + f_0 * mg_533[k];

        t_534[k] = -kg_534[k]
                   + f_0 * mg_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, kg_535, kg_536, kg_537, kg_538, \
                         kg_539, mg_535, mg_536, mg_537, mg_538, \
                         mg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -kg_535[k]
                   + f_0 * mg_535[k];

        t_536[k] = -kg_536[k]
                   + f_0 * mg_536[k];

        t_537[k] = -kg_537[k]
                   + f_0 * mg_537[k];

        t_538[k] = -kg_538[k]
                   + f_0 * mg_538[k];

        t_539[k] = -kg_539[k]
                   + f_0 * mg_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, t_545, t_546, t_547, mg_540, \
                         mg_541, mg_542, mg_543, mg_544, mg_545, mg_546, \
                         mg_547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_0 * mg_540[k];

        t_541[k] = f_0 * mg_541[k];

        t_542[k] = f_0 * mg_542[k];

        t_543[k] = f_0 * mg_543[k];

        t_544[k] = f_0 * mg_544[k];

        t_545[k] = f_0 * mg_545[k];

        t_546[k] = f_0 * mg_546[k];

        t_547[k] = f_0 * mg_547[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, t_553, t_554, t_555, mg_548, \
                         mg_549, mg_550, mg_551, mg_552, mg_553, mg_554, \
                         mg_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_0 * mg_548[k];

        t_549[k] = f_0 * mg_549[k];

        t_550[k] = f_0 * mg_550[k];

        t_551[k] = f_0 * mg_551[k];

        t_552[k] = f_0 * mg_552[k];

        t_553[k] = f_0 * mg_553[k];

        t_554[k] = f_0 * mg_554[k];

        t_555[k] = f_0 * mg_555[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, t_561, t_562, t_563, mg_556, \
                         mg_557, mg_558, mg_559, mg_560, mg_561, mg_562, \
                         mg_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_0 * mg_556[k];

        t_557[k] = f_0 * mg_557[k];

        t_558[k] = f_0 * mg_558[k];

        t_559[k] = f_0 * mg_559[k];

        t_560[k] = f_0 * mg_560[k];

        t_561[k] = f_0 * mg_561[k];

        t_562[k] = f_0 * mg_562[k];

        t_563[k] = f_0 * mg_563[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, t_568, t_569, t_570, t_571, mg_564, \
                         mg_565, mg_566, mg_567, mg_568, mg_569, mg_570, \
                         mg_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_0 * mg_564[k];

        t_565[k] = f_0 * mg_565[k];

        t_566[k] = f_0 * mg_566[k];

        t_567[k] = f_0 * mg_567[k];

        t_568[k] = f_0 * mg_568[k];

        t_569[k] = f_0 * mg_569[k];

        t_570[k] = f_0 * mg_570[k];

        t_571[k] = f_0 * mg_571[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, t_577, t_578, t_579, mg_572, \
                         mg_573, mg_574, mg_575, mg_576, mg_577, mg_578, \
                         mg_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_0 * mg_572[k];

        t_573[k] = f_0 * mg_573[k];

        t_574[k] = f_0 * mg_574[k];

        t_575[k] = f_0 * mg_575[k];

        t_576[k] = f_0 * mg_576[k];

        t_577[k] = f_0 * mg_577[k];

        t_578[k] = f_0 * mg_578[k];

        t_579[k] = f_0 * mg_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, t_585, t_586, t_587, mg_580, \
                         mg_581, mg_582, mg_583, mg_584, mg_585, mg_586, \
                         mg_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_0 * mg_580[k];

        t_581[k] = f_0 * mg_581[k];

        t_582[k] = f_0 * mg_582[k];

        t_583[k] = f_0 * mg_583[k];

        t_584[k] = f_0 * mg_584[k];

        t_585[k] = f_0 * mg_585[k];

        t_586[k] = f_0 * mg_586[k];

        t_587[k] = f_0 * mg_587[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, t_592, t_593, t_594, t_595, mg_588, \
                         mg_589, mg_590, mg_591, mg_592, mg_593, mg_594, \
                         mg_595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_0 * mg_588[k];

        t_589[k] = f_0 * mg_589[k];

        t_590[k] = f_0 * mg_590[k];

        t_591[k] = f_0 * mg_591[k];

        t_592[k] = f_0 * mg_592[k];

        t_593[k] = f_0 * mg_593[k];

        t_594[k] = f_0 * mg_594[k];

        t_595[k] = f_0 * mg_595[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, t_600, t_601, t_602, t_603, mg_596, \
                         mg_597, mg_598, mg_599, mg_600, mg_601, mg_602, \
                         mg_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_0 * mg_596[k];

        t_597[k] = f_0 * mg_597[k];

        t_598[k] = f_0 * mg_598[k];

        t_599[k] = f_0 * mg_599[k];

        t_600[k] = f_0 * mg_600[k];

        t_601[k] = f_0 * mg_601[k];

        t_602[k] = f_0 * mg_602[k];

        t_603[k] = f_0 * mg_603[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, t_609, t_610, t_611, mg_604, \
                         mg_605, mg_606, mg_607, mg_608, mg_609, mg_610, \
                         mg_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_0 * mg_604[k];

        t_605[k] = f_0 * mg_605[k];

        t_606[k] = f_0 * mg_606[k];

        t_607[k] = f_0 * mg_607[k];

        t_608[k] = f_0 * mg_608[k];

        t_609[k] = f_0 * mg_609[k];

        t_610[k] = f_0 * mg_610[k];

        t_611[k] = f_0 * mg_611[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, t_617, t_618, t_619, mg_612, \
                         mg_613, mg_614, mg_615, mg_616, mg_617, mg_618, \
                         mg_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_0 * mg_612[k];

        t_613[k] = f_0 * mg_613[k];

        t_614[k] = f_0 * mg_614[k];

        t_615[k] = f_0 * mg_615[k];

        t_616[k] = f_0 * mg_616[k];

        t_617[k] = f_0 * mg_617[k];

        t_618[k] = f_0 * mg_618[k];

        t_619[k] = f_0 * mg_619[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, t_625, t_626, t_627, mg_620, \
                         mg_621, mg_622, mg_623, mg_624, mg_625, mg_626, \
                         mg_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_0 * mg_620[k];

        t_621[k] = f_0 * mg_621[k];

        t_622[k] = f_0 * mg_622[k];

        t_623[k] = f_0 * mg_623[k];

        t_624[k] = f_0 * mg_624[k];

        t_625[k] = f_0 * mg_625[k];

        t_626[k] = f_0 * mg_626[k];

        t_627[k] = f_0 * mg_627[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t mg, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *mg_628 = buffer.data(mg + 628);
    const auto *mg_629 = buffer.data(mg + 629);
    const auto *mg_630 = buffer.data(mg + 630);
    const auto *mg_631 = buffer.data(mg + 631);
    const auto *mg_632 = buffer.data(mg + 632);
    const auto *mg_633 = buffer.data(mg + 633);
    const auto *mg_634 = buffer.data(mg + 634);
    const auto *mg_635 = buffer.data(mg + 635);
    const auto *mg_636 = buffer.data(mg + 636);
    const auto *mg_637 = buffer.data(mg + 637);
    const auto *mg_638 = buffer.data(mg + 638);
    const auto *mg_639 = buffer.data(mg + 639);
    const auto *mg_640 = buffer.data(mg + 640);
    const auto *mg_641 = buffer.data(mg + 641);
    const auto *mg_642 = buffer.data(mg + 642);
    const auto *mg_643 = buffer.data(mg + 643);
    const auto *mg_644 = buffer.data(mg + 644);
    const auto *mg_645 = buffer.data(mg + 645);
    const auto *mg_646 = buffer.data(mg + 646);
    const auto *mg_647 = buffer.data(mg + 647);
    const auto *mg_648 = buffer.data(mg + 648);
    const auto *mg_649 = buffer.data(mg + 649);
    const auto *mg_650 = buffer.data(mg + 650);
    const auto *mg_651 = buffer.data(mg + 651);
    const auto *mg_652 = buffer.data(mg + 652);
    const auto *mg_653 = buffer.data(mg + 653);
    const auto *mg_654 = buffer.data(mg + 654);
    const auto *mg_655 = buffer.data(mg + 655);
    const auto *mg_656 = buffer.data(mg + 656);
    const auto *mg_657 = buffer.data(mg + 657);
    const auto *mg_658 = buffer.data(mg + 658);
    const auto *mg_659 = buffer.data(mg + 659);
    const auto *mg_660 = buffer.data(mg + 660);
    const auto *mg_661 = buffer.data(mg + 661);
    const auto *mg_662 = buffer.data(mg + 662);
    const auto *mg_663 = buffer.data(mg + 663);
    const auto *mg_664 = buffer.data(mg + 664);
    const auto *mg_665 = buffer.data(mg + 665);
    const auto *mg_666 = buffer.data(mg + 666);
    const auto *mg_667 = buffer.data(mg + 667);
    const auto *mg_668 = buffer.data(mg + 668);
    const auto *mg_669 = buffer.data(mg + 669);
    const auto *mg_670 = buffer.data(mg + 670);
    const auto *mg_671 = buffer.data(mg + 671);
    const auto *mg_672 = buffer.data(mg + 672);
    const auto *mg_673 = buffer.data(mg + 673);
    const auto *mg_674 = buffer.data(mg + 674);

#pragma omp simd aligned(t_628, t_629, t_630, t_631, t_632, t_633, t_634, t_635, mg_628, \
                         mg_629, mg_630, mg_631, mg_632, mg_633, mg_634, \
                         mg_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_0 * mg_628[k];

        t_629[k] = f_0 * mg_629[k];

        t_630[k] = f_0 * mg_630[k];

        t_631[k] = f_0 * mg_631[k];

        t_632[k] = f_0 * mg_632[k];

        t_633[k] = f_0 * mg_633[k];

        t_634[k] = f_0 * mg_634[k];

        t_635[k] = f_0 * mg_635[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, t_640, t_641, t_642, t_643, mg_636, \
                         mg_637, mg_638, mg_639, mg_640, mg_641, mg_642, \
                         mg_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_0 * mg_636[k];

        t_637[k] = f_0 * mg_637[k];

        t_638[k] = f_0 * mg_638[k];

        t_639[k] = f_0 * mg_639[k];

        t_640[k] = f_0 * mg_640[k];

        t_641[k] = f_0 * mg_641[k];

        t_642[k] = f_0 * mg_642[k];

        t_643[k] = f_0 * mg_643[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, t_648, t_649, t_650, t_651, mg_644, \
                         mg_645, mg_646, mg_647, mg_648, mg_649, mg_650, \
                         mg_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_0 * mg_644[k];

        t_645[k] = f_0 * mg_645[k];

        t_646[k] = f_0 * mg_646[k];

        t_647[k] = f_0 * mg_647[k];

        t_648[k] = f_0 * mg_648[k];

        t_649[k] = f_0 * mg_649[k];

        t_650[k] = f_0 * mg_650[k];

        t_651[k] = f_0 * mg_651[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, t_656, t_657, t_658, t_659, mg_652, \
                         mg_653, mg_654, mg_655, mg_656, mg_657, mg_658, \
                         mg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_0 * mg_652[k];

        t_653[k] = f_0 * mg_653[k];

        t_654[k] = f_0 * mg_654[k];

        t_655[k] = f_0 * mg_655[k];

        t_656[k] = f_0 * mg_656[k];

        t_657[k] = f_0 * mg_657[k];

        t_658[k] = f_0 * mg_658[k];

        t_659[k] = f_0 * mg_659[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, t_665, t_666, t_667, mg_660, \
                         mg_661, mg_662, mg_663, mg_664, mg_665, mg_666, \
                         mg_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_0 * mg_660[k];

        t_661[k] = f_0 * mg_661[k];

        t_662[k] = f_0 * mg_662[k];

        t_663[k] = f_0 * mg_663[k];

        t_664[k] = f_0 * mg_664[k];

        t_665[k] = f_0 * mg_665[k];

        t_666[k] = f_0 * mg_666[k];

        t_667[k] = f_0 * mg_667[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, t_672, t_673, t_674, mg_668, mg_669, \
                         mg_670, mg_671, mg_672, mg_673, mg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_0 * mg_668[k];

        t_669[k] = f_0 * mg_669[k];

        t_670[k] = f_0 * mg_670[k];

        t_671[k] = f_0 * mg_671[k];

        t_672[k] = f_0 * mg_672[k];

        t_673[k] = f_0 * mg_673[k];

        t_674[k] = f_0 * mg_674[k];
    }
}

auto
compute_prim_geom_10_lg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t kg, const size_t mg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_lg_electron_repulsion_0_piece0(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_0_piece1(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_0_piece2(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_0_piece3(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_0_piece4(buffer, target, mg, ncols, alpha);
}

static auto
compute_prim_geom_10_lg_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
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

    const auto *mg_15 = buffer.data(mg + 15);
    const auto *mg_16 = buffer.data(mg + 16);
    const auto *mg_17 = buffer.data(mg + 17);
    const auto *mg_18 = buffer.data(mg + 18);
    const auto *mg_19 = buffer.data(mg + 19);
    const auto *mg_20 = buffer.data(mg + 20);
    const auto *mg_21 = buffer.data(mg + 21);
    const auto *mg_22 = buffer.data(mg + 22);
    const auto *mg_23 = buffer.data(mg + 23);
    const auto *mg_24 = buffer.data(mg + 24);
    const auto *mg_25 = buffer.data(mg + 25);
    const auto *mg_26 = buffer.data(mg + 26);
    const auto *mg_27 = buffer.data(mg + 27);
    const auto *mg_28 = buffer.data(mg + 28);
    const auto *mg_29 = buffer.data(mg + 29);
    const auto *mg_45 = buffer.data(mg + 45);
    const auto *mg_46 = buffer.data(mg + 46);
    const auto *mg_47 = buffer.data(mg + 47);
    const auto *mg_48 = buffer.data(mg + 48);
    const auto *mg_49 = buffer.data(mg + 49);
    const auto *mg_50 = buffer.data(mg + 50);
    const auto *mg_51 = buffer.data(mg + 51);
    const auto *mg_52 = buffer.data(mg + 52);
    const auto *mg_53 = buffer.data(mg + 53);
    const auto *mg_54 = buffer.data(mg + 54);
    const auto *mg_55 = buffer.data(mg + 55);
    const auto *mg_56 = buffer.data(mg + 56);
    const auto *mg_57 = buffer.data(mg + 57);
    const auto *mg_58 = buffer.data(mg + 58);
    const auto *mg_59 = buffer.data(mg + 59);
    const auto *mg_60 = buffer.data(mg + 60);
    const auto *mg_61 = buffer.data(mg + 61);
    const auto *mg_62 = buffer.data(mg + 62);
    const auto *mg_63 = buffer.data(mg + 63);
    const auto *mg_64 = buffer.data(mg + 64);
    const auto *mg_65 = buffer.data(mg + 65);
    const auto *mg_66 = buffer.data(mg + 66);
    const auto *mg_67 = buffer.data(mg + 67);
    const auto *mg_68 = buffer.data(mg + 68);
    const auto *mg_69 = buffer.data(mg + 69);
    const auto *mg_70 = buffer.data(mg + 70);
    const auto *mg_71 = buffer.data(mg + 71);
    const auto *mg_72 = buffer.data(mg + 72);
    const auto *mg_73 = buffer.data(mg + 73);
    const auto *mg_74 = buffer.data(mg + 74);
    const auto *mg_90 = buffer.data(mg + 90);
    const auto *mg_91 = buffer.data(mg + 91);
    const auto *mg_92 = buffer.data(mg + 92);
    const auto *mg_93 = buffer.data(mg + 93);
    const auto *mg_94 = buffer.data(mg + 94);
    const auto *mg_95 = buffer.data(mg + 95);
    const auto *mg_96 = buffer.data(mg + 96);
    const auto *mg_97 = buffer.data(mg + 97);
    const auto *mg_98 = buffer.data(mg + 98);
    const auto *mg_99 = buffer.data(mg + 99);
    const auto *mg_100 = buffer.data(mg + 100);
    const auto *mg_101 = buffer.data(mg + 101);
    const auto *mg_102 = buffer.data(mg + 102);
    const auto *mg_103 = buffer.data(mg + 103);
    const auto *mg_104 = buffer.data(mg + 104);
    const auto *mg_105 = buffer.data(mg + 105);
    const auto *mg_106 = buffer.data(mg + 106);
    const auto *mg_107 = buffer.data(mg + 107);
    const auto *mg_108 = buffer.data(mg + 108);
    const auto *mg_109 = buffer.data(mg + 109);
    const auto *mg_110 = buffer.data(mg + 110);
    const auto *mg_111 = buffer.data(mg + 111);
    const auto *mg_112 = buffer.data(mg + 112);
    const auto *mg_113 = buffer.data(mg + 113);
    const auto *mg_114 = buffer.data(mg + 114);
    const auto *mg_115 = buffer.data(mg + 115);
    const auto *mg_116 = buffer.data(mg + 116);
    const auto *mg_117 = buffer.data(mg + 117);
    const auto *mg_118 = buffer.data(mg + 118);
    const auto *mg_119 = buffer.data(mg + 119);
    const auto *mg_120 = buffer.data(mg + 120);
    const auto *mg_121 = buffer.data(mg + 121);
    const auto *mg_122 = buffer.data(mg + 122);
    const auto *mg_123 = buffer.data(mg + 123);
    const auto *mg_124 = buffer.data(mg + 124);
    const auto *mg_125 = buffer.data(mg + 125);
    const auto *mg_126 = buffer.data(mg + 126);
    const auto *mg_127 = buffer.data(mg + 127);
    const auto *mg_128 = buffer.data(mg + 128);
    const auto *mg_129 = buffer.data(mg + 129);
    const auto *mg_130 = buffer.data(mg + 130);
    const auto *mg_131 = buffer.data(mg + 131);
    const auto *mg_132 = buffer.data(mg + 132);
    const auto *mg_133 = buffer.data(mg + 133);
    const auto *mg_134 = buffer.data(mg + 134);
    const auto *mg_150 = buffer.data(mg + 150);
    const auto *mg_151 = buffer.data(mg + 151);
    const auto *mg_152 = buffer.data(mg + 152);
    const auto *mg_153 = buffer.data(mg + 153);
    const auto *mg_154 = buffer.data(mg + 154);
    const auto *mg_155 = buffer.data(mg + 155);
    const auto *mg_156 = buffer.data(mg + 156);
    const auto *mg_157 = buffer.data(mg + 157);
    const auto *mg_158 = buffer.data(mg + 158);
    const auto *mg_159 = buffer.data(mg + 159);
    const auto *mg_160 = buffer.data(mg + 160);
    const auto *mg_161 = buffer.data(mg + 161);
    const auto *mg_162 = buffer.data(mg + 162);
    const auto *mg_163 = buffer.data(mg + 163);
    const auto *mg_164 = buffer.data(mg + 164);
    const auto *mg_165 = buffer.data(mg + 165);
    const auto *mg_166 = buffer.data(mg + 166);
    const auto *mg_167 = buffer.data(mg + 167);
    const auto *mg_168 = buffer.data(mg + 168);
    const auto *mg_169 = buffer.data(mg + 169);
    const auto *mg_170 = buffer.data(mg + 170);
    const auto *mg_171 = buffer.data(mg + 171);
    const auto *mg_172 = buffer.data(mg + 172);
    const auto *mg_173 = buffer.data(mg + 173);
    const auto *mg_174 = buffer.data(mg + 174);
    const auto *mg_175 = buffer.data(mg + 175);
    const auto *mg_176 = buffer.data(mg + 176);
    const auto *mg_177 = buffer.data(mg + 177);
    const auto *mg_178 = buffer.data(mg + 178);
    const auto *mg_179 = buffer.data(mg + 179);
    const auto *mg_180 = buffer.data(mg + 180);
    const auto *mg_181 = buffer.data(mg + 181);
    const auto *mg_182 = buffer.data(mg + 182);
    const auto *mg_183 = buffer.data(mg + 183);
    const auto *mg_184 = buffer.data(mg + 184);
    const auto *mg_185 = buffer.data(mg + 185);
    const auto *mg_186 = buffer.data(mg + 186);
    const auto *mg_187 = buffer.data(mg + 187);
    const auto *mg_188 = buffer.data(mg + 188);
    const auto *mg_189 = buffer.data(mg + 189);
    const auto *mg_190 = buffer.data(mg + 190);
    const auto *mg_191 = buffer.data(mg + 191);
    const auto *mg_192 = buffer.data(mg + 192);
    const auto *mg_193 = buffer.data(mg + 193);
    const auto *mg_194 = buffer.data(mg + 194);
    const auto *mg_195 = buffer.data(mg + 195);
    const auto *mg_196 = buffer.data(mg + 196);
    const auto *mg_197 = buffer.data(mg + 197);
    const auto *mg_198 = buffer.data(mg + 198);
    const auto *mg_199 = buffer.data(mg + 199);
    const auto *mg_200 = buffer.data(mg + 200);
    const auto *mg_201 = buffer.data(mg + 201);
    const auto *mg_202 = buffer.data(mg + 202);
    const auto *mg_203 = buffer.data(mg + 203);
    const auto *mg_204 = buffer.data(mg + 204);
    const auto *mg_205 = buffer.data(mg + 205);
    const auto *mg_206 = buffer.data(mg + 206);
    const auto *mg_207 = buffer.data(mg + 207);
    const auto *mg_208 = buffer.data(mg + 208);
    const auto *mg_209 = buffer.data(mg + 209);
    const auto *mg_225 = buffer.data(mg + 225);
    const auto *mg_226 = buffer.data(mg + 226);
    const auto *mg_227 = buffer.data(mg + 227);
    const auto *mg_228 = buffer.data(mg + 228);
    const auto *mg_229 = buffer.data(mg + 229);
    const auto *mg_230 = buffer.data(mg + 230);
    const auto *mg_231 = buffer.data(mg + 231);
    const auto *mg_232 = buffer.data(mg + 232);
    const auto *mg_233 = buffer.data(mg + 233);
    const auto *mg_234 = buffer.data(mg + 234);
    const auto *mg_235 = buffer.data(mg + 235);
    const auto *mg_236 = buffer.data(mg + 236);
    const auto *mg_237 = buffer.data(mg + 237);
    const auto *mg_238 = buffer.data(mg + 238);
    const auto *mg_239 = buffer.data(mg + 239);
    const auto *mg_240 = buffer.data(mg + 240);
    const auto *mg_241 = buffer.data(mg + 241);
    const auto *mg_242 = buffer.data(mg + 242);
    const auto *mg_243 = buffer.data(mg + 243);
    const auto *mg_244 = buffer.data(mg + 244);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, mg_15, mg_16, mg_17, mg_18, \
                         mg_19, mg_20, mg_21, mg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * mg_15[k];

        t_1[k] = f_0 * mg_16[k];

        t_2[k] = f_0 * mg_17[k];

        t_3[k] = f_0 * mg_18[k];

        t_4[k] = f_0 * mg_19[k];

        t_5[k] = f_0 * mg_20[k];

        t_6[k] = f_0 * mg_21[k];

        t_7[k] = f_0 * mg_22[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, mg_23, mg_24, mg_25, mg_26, \
                         mg_27, mg_28, mg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * mg_23[k];

        t_9[k] = f_0 * mg_24[k];

        t_10[k] = f_0 * mg_25[k];

        t_11[k] = f_0 * mg_26[k];

        t_12[k] = f_0 * mg_27[k];

        t_13[k] = f_0 * mg_28[k];

        t_14[k] = f_0 * mg_29[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, kg_0, kg_1, kg_2, kg_3, kg_4, mg_45, \
                         mg_46, mg_47, mg_48, mg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -kg_0[k]
                  + f_0 * mg_45[k];

        t_16[k] = -kg_1[k]
                  + f_0 * mg_46[k];

        t_17[k] = -kg_2[k]
                  + f_0 * mg_47[k];

        t_18[k] = -kg_3[k]
                  + f_0 * mg_48[k];

        t_19[k] = -kg_4[k]
                  + f_0 * mg_49[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, kg_5, kg_6, kg_7, kg_8, kg_9, mg_50, \
                         mg_51, mg_52, mg_53, mg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -kg_5[k]
                  + f_0 * mg_50[k];

        t_21[k] = -kg_6[k]
                  + f_0 * mg_51[k];

        t_22[k] = -kg_7[k]
                  + f_0 * mg_52[k];

        t_23[k] = -kg_8[k]
                  + f_0 * mg_53[k];

        t_24[k] = -kg_9[k]
                  + f_0 * mg_54[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, kg_10, kg_11, kg_12, kg_13, kg_14, \
                         mg_55, mg_56, mg_57, mg_58, mg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -kg_10[k]
                  + f_0 * mg_55[k];

        t_26[k] = -kg_11[k]
                  + f_0 * mg_56[k];

        t_27[k] = -kg_12[k]
                  + f_0 * mg_57[k];

        t_28[k] = -kg_13[k]
                  + f_0 * mg_58[k];

        t_29[k] = -kg_14[k]
                  + f_0 * mg_59[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, t_36, t_37, mg_60, mg_61, mg_62, \
                         mg_63, mg_64, mg_65, mg_66, mg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * mg_60[k];

        t_31[k] = f_0 * mg_61[k];

        t_32[k] = f_0 * mg_62[k];

        t_33[k] = f_0 * mg_63[k];

        t_34[k] = f_0 * mg_64[k];

        t_35[k] = f_0 * mg_65[k];

        t_36[k] = f_0 * mg_66[k];

        t_37[k] = f_0 * mg_67[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, mg_68, mg_69, mg_70, mg_71, \
                         mg_72, mg_73, mg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * mg_68[k];

        t_39[k] = f_0 * mg_69[k];

        t_40[k] = f_0 * mg_70[k];

        t_41[k] = f_0 * mg_71[k];

        t_42[k] = f_0 * mg_72[k];

        t_43[k] = f_0 * mg_73[k];

        t_44[k] = f_0 * mg_74[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, kg_15, kg_16, kg_17, kg_18, kg_19, \
                         mg_90, mg_91, mg_92, mg_93, mg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * kg_15[k]
                  + f_0 * mg_90[k];

        t_46[k] = -2.0 * kg_16[k]
                  + f_0 * mg_91[k];

        t_47[k] = -2.0 * kg_17[k]
                  + f_0 * mg_92[k];

        t_48[k] = -2.0 * kg_18[k]
                  + f_0 * mg_93[k];

        t_49[k] = -2.0 * kg_19[k]
                  + f_0 * mg_94[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, kg_20, kg_21, kg_22, kg_23, kg_24, \
                         mg_95, mg_96, mg_97, mg_98, mg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * kg_20[k]
                  + f_0 * mg_95[k];

        t_51[k] = -2.0 * kg_21[k]
                  + f_0 * mg_96[k];

        t_52[k] = -2.0 * kg_22[k]
                  + f_0 * mg_97[k];

        t_53[k] = -2.0 * kg_23[k]
                  + f_0 * mg_98[k];

        t_54[k] = -2.0 * kg_24[k]
                  + f_0 * mg_99[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, kg_25, kg_26, kg_27, kg_28, kg_29, \
                         mg_100, mg_101, mg_102, mg_103, mg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * kg_25[k]
                  + f_0 * mg_100[k];

        t_56[k] = -2.0 * kg_26[k]
                  + f_0 * mg_101[k];

        t_57[k] = -2.0 * kg_27[k]
                  + f_0 * mg_102[k];

        t_58[k] = -2.0 * kg_28[k]
                  + f_0 * mg_103[k];

        t_59[k] = -2.0 * kg_29[k]
                  + f_0 * mg_104[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, kg_30, kg_31, kg_32, kg_33, kg_34, \
                         mg_105, mg_106, mg_107, mg_108, mg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -kg_30[k]
                  + f_0 * mg_105[k];

        t_61[k] = -kg_31[k]
                  + f_0 * mg_106[k];

        t_62[k] = -kg_32[k]
                  + f_0 * mg_107[k];

        t_63[k] = -kg_33[k]
                  + f_0 * mg_108[k];

        t_64[k] = -kg_34[k]
                  + f_0 * mg_109[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, kg_35, kg_36, kg_37, kg_38, kg_39, \
                         mg_110, mg_111, mg_112, mg_113, mg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -kg_35[k]
                  + f_0 * mg_110[k];

        t_66[k] = -kg_36[k]
                  + f_0 * mg_111[k];

        t_67[k] = -kg_37[k]
                  + f_0 * mg_112[k];

        t_68[k] = -kg_38[k]
                  + f_0 * mg_113[k];

        t_69[k] = -kg_39[k]
                  + f_0 * mg_114[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, kg_40, kg_41, kg_42, kg_43, kg_44, \
                         mg_115, mg_116, mg_117, mg_118, mg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -kg_40[k]
                  + f_0 * mg_115[k];

        t_71[k] = -kg_41[k]
                  + f_0 * mg_116[k];

        t_72[k] = -kg_42[k]
                  + f_0 * mg_117[k];

        t_73[k] = -kg_43[k]
                  + f_0 * mg_118[k];

        t_74[k] = -kg_44[k]
                  + f_0 * mg_119[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, t_80, t_81, t_82, mg_120, mg_121, \
                         mg_122, mg_123, mg_124, mg_125, mg_126, \
                         mg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * mg_120[k];

        t_76[k] = f_0 * mg_121[k];

        t_77[k] = f_0 * mg_122[k];

        t_78[k] = f_0 * mg_123[k];

        t_79[k] = f_0 * mg_124[k];

        t_80[k] = f_0 * mg_125[k];

        t_81[k] = f_0 * mg_126[k];

        t_82[k] = f_0 * mg_127[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, mg_128, mg_129, mg_130, \
                         mg_131, mg_132, mg_133, mg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * mg_128[k];

        t_84[k] = f_0 * mg_129[k];

        t_85[k] = f_0 * mg_130[k];

        t_86[k] = f_0 * mg_131[k];

        t_87[k] = f_0 * mg_132[k];

        t_88[k] = f_0 * mg_133[k];

        t_89[k] = f_0 * mg_134[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, kg_45, kg_46, kg_47, kg_48, kg_49, \
                         mg_150, mg_151, mg_152, mg_153, mg_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * kg_45[k]
                  + f_0 * mg_150[k];

        t_91[k] = -3.0 * kg_46[k]
                  + f_0 * mg_151[k];

        t_92[k] = -3.0 * kg_47[k]
                  + f_0 * mg_152[k];

        t_93[k] = -3.0 * kg_48[k]
                  + f_0 * mg_153[k];

        t_94[k] = -3.0 * kg_49[k]
                  + f_0 * mg_154[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, kg_50, kg_51, kg_52, kg_53, kg_54, \
                         mg_155, mg_156, mg_157, mg_158, mg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * kg_50[k]
                  + f_0 * mg_155[k];

        t_96[k] = -3.0 * kg_51[k]
                  + f_0 * mg_156[k];

        t_97[k] = -3.0 * kg_52[k]
                  + f_0 * mg_157[k];

        t_98[k] = -3.0 * kg_53[k]
                  + f_0 * mg_158[k];

        t_99[k] = -3.0 * kg_54[k]
                  + f_0 * mg_159[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, kg_55, kg_56, kg_57, kg_58, kg_59, \
                         mg_160, mg_161, mg_162, mg_163, mg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * kg_55[k]
                   + f_0 * mg_160[k];

        t_101[k] = -3.0 * kg_56[k]
                   + f_0 * mg_161[k];

        t_102[k] = -3.0 * kg_57[k]
                   + f_0 * mg_162[k];

        t_103[k] = -3.0 * kg_58[k]
                   + f_0 * mg_163[k];

        t_104[k] = -3.0 * kg_59[k]
                   + f_0 * mg_164[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, kg_60, kg_61, kg_62, kg_63, kg_64, \
                         mg_165, mg_166, mg_167, mg_168, mg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * kg_60[k]
                   + f_0 * mg_165[k];

        t_106[k] = -2.0 * kg_61[k]
                   + f_0 * mg_166[k];

        t_107[k] = -2.0 * kg_62[k]
                   + f_0 * mg_167[k];

        t_108[k] = -2.0 * kg_63[k]
                   + f_0 * mg_168[k];

        t_109[k] = -2.0 * kg_64[k]
                   + f_0 * mg_169[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, kg_65, kg_66, kg_67, kg_68, kg_69, \
                         mg_170, mg_171, mg_172, mg_173, mg_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * kg_65[k]
                   + f_0 * mg_170[k];

        t_111[k] = -2.0 * kg_66[k]
                   + f_0 * mg_171[k];

        t_112[k] = -2.0 * kg_67[k]
                   + f_0 * mg_172[k];

        t_113[k] = -2.0 * kg_68[k]
                   + f_0 * mg_173[k];

        t_114[k] = -2.0 * kg_69[k]
                   + f_0 * mg_174[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, kg_70, kg_71, kg_72, kg_73, kg_74, \
                         mg_175, mg_176, mg_177, mg_178, mg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * kg_70[k]
                   + f_0 * mg_175[k];

        t_116[k] = -2.0 * kg_71[k]
                   + f_0 * mg_176[k];

        t_117[k] = -2.0 * kg_72[k]
                   + f_0 * mg_177[k];

        t_118[k] = -2.0 * kg_73[k]
                   + f_0 * mg_178[k];

        t_119[k] = -2.0 * kg_74[k]
                   + f_0 * mg_179[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, kg_75, kg_76, kg_77, kg_78, kg_79, \
                         mg_180, mg_181, mg_182, mg_183, mg_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -kg_75[k]
                   + f_0 * mg_180[k];

        t_121[k] = -kg_76[k]
                   + f_0 * mg_181[k];

        t_122[k] = -kg_77[k]
                   + f_0 * mg_182[k];

        t_123[k] = -kg_78[k]
                   + f_0 * mg_183[k];

        t_124[k] = -kg_79[k]
                   + f_0 * mg_184[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, kg_80, kg_81, kg_82, kg_83, kg_84, \
                         mg_185, mg_186, mg_187, mg_188, mg_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -kg_80[k]
                   + f_0 * mg_185[k];

        t_126[k] = -kg_81[k]
                   + f_0 * mg_186[k];

        t_127[k] = -kg_82[k]
                   + f_0 * mg_187[k];

        t_128[k] = -kg_83[k]
                   + f_0 * mg_188[k];

        t_129[k] = -kg_84[k]
                   + f_0 * mg_189[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, kg_85, kg_86, kg_87, kg_88, kg_89, \
                         mg_190, mg_191, mg_192, mg_193, mg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -kg_85[k]
                   + f_0 * mg_190[k];

        t_131[k] = -kg_86[k]
                   + f_0 * mg_191[k];

        t_132[k] = -kg_87[k]
                   + f_0 * mg_192[k];

        t_133[k] = -kg_88[k]
                   + f_0 * mg_193[k];

        t_134[k] = -kg_89[k]
                   + f_0 * mg_194[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, t_140, t_141, t_142, mg_195, \
                         mg_196, mg_197, mg_198, mg_199, mg_200, mg_201, \
                         mg_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * mg_195[k];

        t_136[k] = f_0 * mg_196[k];

        t_137[k] = f_0 * mg_197[k];

        t_138[k] = f_0 * mg_198[k];

        t_139[k] = f_0 * mg_199[k];

        t_140[k] = f_0 * mg_200[k];

        t_141[k] = f_0 * mg_201[k];

        t_142[k] = f_0 * mg_202[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, t_148, t_149, mg_203, mg_204, \
                         mg_205, mg_206, mg_207, mg_208, mg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_0 * mg_203[k];

        t_144[k] = f_0 * mg_204[k];

        t_145[k] = f_0 * mg_205[k];

        t_146[k] = f_0 * mg_206[k];

        t_147[k] = f_0 * mg_207[k];

        t_148[k] = f_0 * mg_208[k];

        t_149[k] = f_0 * mg_209[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, kg_90, kg_91, kg_92, kg_93, kg_94, \
                         mg_225, mg_226, mg_227, mg_228, mg_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -4.0 * kg_90[k]
                   + f_0 * mg_225[k];

        t_151[k] = -4.0 * kg_91[k]
                   + f_0 * mg_226[k];

        t_152[k] = -4.0 * kg_92[k]
                   + f_0 * mg_227[k];

        t_153[k] = -4.0 * kg_93[k]
                   + f_0 * mg_228[k];

        t_154[k] = -4.0 * kg_94[k]
                   + f_0 * mg_229[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, kg_95, kg_96, kg_97, kg_98, kg_99, \
                         mg_230, mg_231, mg_232, mg_233, mg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -4.0 * kg_95[k]
                   + f_0 * mg_230[k];

        t_156[k] = -4.0 * kg_96[k]
                   + f_0 * mg_231[k];

        t_157[k] = -4.0 * kg_97[k]
                   + f_0 * mg_232[k];

        t_158[k] = -4.0 * kg_98[k]
                   + f_0 * mg_233[k];

        t_159[k] = -4.0 * kg_99[k]
                   + f_0 * mg_234[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, kg_100, kg_101, kg_102, kg_103, \
                         kg_104, mg_235, mg_236, mg_237, mg_238, \
                         mg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -4.0 * kg_100[k]
                   + f_0 * mg_235[k];

        t_161[k] = -4.0 * kg_101[k]
                   + f_0 * mg_236[k];

        t_162[k] = -4.0 * kg_102[k]
                   + f_0 * mg_237[k];

        t_163[k] = -4.0 * kg_103[k]
                   + f_0 * mg_238[k];

        t_164[k] = -4.0 * kg_104[k]
                   + f_0 * mg_239[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, kg_105, kg_106, kg_107, kg_108, \
                         kg_109, mg_240, mg_241, mg_242, mg_243, \
                         mg_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -3.0 * kg_105[k]
                   + f_0 * mg_240[k];

        t_166[k] = -3.0 * kg_106[k]
                   + f_0 * mg_241[k];

        t_167[k] = -3.0 * kg_107[k]
                   + f_0 * mg_242[k];

        t_168[k] = -3.0 * kg_108[k]
                   + f_0 * mg_243[k];

        t_169[k] = -3.0 * kg_109[k]
                   + f_0 * mg_244[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
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

    const auto *mg_245 = buffer.data(mg + 245);
    const auto *mg_246 = buffer.data(mg + 246);
    const auto *mg_247 = buffer.data(mg + 247);
    const auto *mg_248 = buffer.data(mg + 248);
    const auto *mg_249 = buffer.data(mg + 249);
    const auto *mg_250 = buffer.data(mg + 250);
    const auto *mg_251 = buffer.data(mg + 251);
    const auto *mg_252 = buffer.data(mg + 252);
    const auto *mg_253 = buffer.data(mg + 253);
    const auto *mg_254 = buffer.data(mg + 254);
    const auto *mg_255 = buffer.data(mg + 255);
    const auto *mg_256 = buffer.data(mg + 256);
    const auto *mg_257 = buffer.data(mg + 257);
    const auto *mg_258 = buffer.data(mg + 258);
    const auto *mg_259 = buffer.data(mg + 259);
    const auto *mg_260 = buffer.data(mg + 260);
    const auto *mg_261 = buffer.data(mg + 261);
    const auto *mg_262 = buffer.data(mg + 262);
    const auto *mg_263 = buffer.data(mg + 263);
    const auto *mg_264 = buffer.data(mg + 264);
    const auto *mg_265 = buffer.data(mg + 265);
    const auto *mg_266 = buffer.data(mg + 266);
    const auto *mg_267 = buffer.data(mg + 267);
    const auto *mg_268 = buffer.data(mg + 268);
    const auto *mg_269 = buffer.data(mg + 269);
    const auto *mg_270 = buffer.data(mg + 270);
    const auto *mg_271 = buffer.data(mg + 271);
    const auto *mg_272 = buffer.data(mg + 272);
    const auto *mg_273 = buffer.data(mg + 273);
    const auto *mg_274 = buffer.data(mg + 274);
    const auto *mg_275 = buffer.data(mg + 275);
    const auto *mg_276 = buffer.data(mg + 276);
    const auto *mg_277 = buffer.data(mg + 277);
    const auto *mg_278 = buffer.data(mg + 278);
    const auto *mg_279 = buffer.data(mg + 279);
    const auto *mg_280 = buffer.data(mg + 280);
    const auto *mg_281 = buffer.data(mg + 281);
    const auto *mg_282 = buffer.data(mg + 282);
    const auto *mg_283 = buffer.data(mg + 283);
    const auto *mg_284 = buffer.data(mg + 284);
    const auto *mg_285 = buffer.data(mg + 285);
    const auto *mg_286 = buffer.data(mg + 286);
    const auto *mg_287 = buffer.data(mg + 287);
    const auto *mg_288 = buffer.data(mg + 288);
    const auto *mg_289 = buffer.data(mg + 289);
    const auto *mg_290 = buffer.data(mg + 290);
    const auto *mg_291 = buffer.data(mg + 291);
    const auto *mg_292 = buffer.data(mg + 292);
    const auto *mg_293 = buffer.data(mg + 293);
    const auto *mg_294 = buffer.data(mg + 294);
    const auto *mg_295 = buffer.data(mg + 295);
    const auto *mg_296 = buffer.data(mg + 296);
    const auto *mg_297 = buffer.data(mg + 297);
    const auto *mg_298 = buffer.data(mg + 298);
    const auto *mg_299 = buffer.data(mg + 299);
    const auto *mg_315 = buffer.data(mg + 315);
    const auto *mg_316 = buffer.data(mg + 316);
    const auto *mg_317 = buffer.data(mg + 317);
    const auto *mg_318 = buffer.data(mg + 318);
    const auto *mg_319 = buffer.data(mg + 319);
    const auto *mg_320 = buffer.data(mg + 320);
    const auto *mg_321 = buffer.data(mg + 321);
    const auto *mg_322 = buffer.data(mg + 322);
    const auto *mg_323 = buffer.data(mg + 323);
    const auto *mg_324 = buffer.data(mg + 324);
    const auto *mg_325 = buffer.data(mg + 325);
    const auto *mg_326 = buffer.data(mg + 326);
    const auto *mg_327 = buffer.data(mg + 327);
    const auto *mg_328 = buffer.data(mg + 328);
    const auto *mg_329 = buffer.data(mg + 329);
    const auto *mg_330 = buffer.data(mg + 330);
    const auto *mg_331 = buffer.data(mg + 331);
    const auto *mg_332 = buffer.data(mg + 332);
    const auto *mg_333 = buffer.data(mg + 333);
    const auto *mg_334 = buffer.data(mg + 334);
    const auto *mg_335 = buffer.data(mg + 335);
    const auto *mg_336 = buffer.data(mg + 336);
    const auto *mg_337 = buffer.data(mg + 337);
    const auto *mg_338 = buffer.data(mg + 338);
    const auto *mg_339 = buffer.data(mg + 339);
    const auto *mg_340 = buffer.data(mg + 340);
    const auto *mg_341 = buffer.data(mg + 341);
    const auto *mg_342 = buffer.data(mg + 342);
    const auto *mg_343 = buffer.data(mg + 343);
    const auto *mg_344 = buffer.data(mg + 344);
    const auto *mg_345 = buffer.data(mg + 345);
    const auto *mg_346 = buffer.data(mg + 346);
    const auto *mg_347 = buffer.data(mg + 347);
    const auto *mg_348 = buffer.data(mg + 348);
    const auto *mg_349 = buffer.data(mg + 349);
    const auto *mg_350 = buffer.data(mg + 350);
    const auto *mg_351 = buffer.data(mg + 351);
    const auto *mg_352 = buffer.data(mg + 352);
    const auto *mg_353 = buffer.data(mg + 353);
    const auto *mg_354 = buffer.data(mg + 354);
    const auto *mg_355 = buffer.data(mg + 355);
    const auto *mg_356 = buffer.data(mg + 356);
    const auto *mg_357 = buffer.data(mg + 357);
    const auto *mg_358 = buffer.data(mg + 358);
    const auto *mg_359 = buffer.data(mg + 359);
    const auto *mg_360 = buffer.data(mg + 360);
    const auto *mg_361 = buffer.data(mg + 361);
    const auto *mg_362 = buffer.data(mg + 362);
    const auto *mg_363 = buffer.data(mg + 363);
    const auto *mg_364 = buffer.data(mg + 364);
    const auto *mg_365 = buffer.data(mg + 365);
    const auto *mg_366 = buffer.data(mg + 366);
    const auto *mg_367 = buffer.data(mg + 367);
    const auto *mg_368 = buffer.data(mg + 368);
    const auto *mg_369 = buffer.data(mg + 369);
    const auto *mg_370 = buffer.data(mg + 370);
    const auto *mg_371 = buffer.data(mg + 371);
    const auto *mg_372 = buffer.data(mg + 372);
    const auto *mg_373 = buffer.data(mg + 373);
    const auto *mg_374 = buffer.data(mg + 374);
    const auto *mg_375 = buffer.data(mg + 375);
    const auto *mg_376 = buffer.data(mg + 376);
    const auto *mg_377 = buffer.data(mg + 377);
    const auto *mg_378 = buffer.data(mg + 378);
    const auto *mg_379 = buffer.data(mg + 379);
    const auto *mg_380 = buffer.data(mg + 380);
    const auto *mg_381 = buffer.data(mg + 381);
    const auto *mg_382 = buffer.data(mg + 382);
    const auto *mg_383 = buffer.data(mg + 383);
    const auto *mg_384 = buffer.data(mg + 384);
    const auto *mg_385 = buffer.data(mg + 385);
    const auto *mg_386 = buffer.data(mg + 386);
    const auto *mg_387 = buffer.data(mg + 387);
    const auto *mg_388 = buffer.data(mg + 388);
    const auto *mg_389 = buffer.data(mg + 389);
    const auto *mg_390 = buffer.data(mg + 390);
    const auto *mg_391 = buffer.data(mg + 391);
    const auto *mg_392 = buffer.data(mg + 392);
    const auto *mg_393 = buffer.data(mg + 393);
    const auto *mg_394 = buffer.data(mg + 394);
    const auto *mg_395 = buffer.data(mg + 395);
    const auto *mg_396 = buffer.data(mg + 396);
    const auto *mg_397 = buffer.data(mg + 397);
    const auto *mg_398 = buffer.data(mg + 398);
    const auto *mg_399 = buffer.data(mg + 399);
    const auto *mg_400 = buffer.data(mg + 400);
    const auto *mg_401 = buffer.data(mg + 401);
    const auto *mg_402 = buffer.data(mg + 402);
    const auto *mg_403 = buffer.data(mg + 403);
    const auto *mg_404 = buffer.data(mg + 404);
    const auto *mg_420 = buffer.data(mg + 420);
    const auto *mg_421 = buffer.data(mg + 421);
    const auto *mg_422 = buffer.data(mg + 422);
    const auto *mg_423 = buffer.data(mg + 423);
    const auto *mg_424 = buffer.data(mg + 424);
    const auto *mg_425 = buffer.data(mg + 425);
    const auto *mg_426 = buffer.data(mg + 426);
    const auto *mg_427 = buffer.data(mg + 427);
    const auto *mg_428 = buffer.data(mg + 428);
    const auto *mg_429 = buffer.data(mg + 429);
    const auto *mg_430 = buffer.data(mg + 430);
    const auto *mg_431 = buffer.data(mg + 431);
    const auto *mg_432 = buffer.data(mg + 432);
    const auto *mg_433 = buffer.data(mg + 433);
    const auto *mg_434 = buffer.data(mg + 434);

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, kg_110, kg_111, kg_112, kg_113, \
                         kg_114, mg_245, mg_246, mg_247, mg_248, \
                         mg_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * kg_110[k]
                   + f_0 * mg_245[k];

        t_171[k] = -3.0 * kg_111[k]
                   + f_0 * mg_246[k];

        t_172[k] = -3.0 * kg_112[k]
                   + f_0 * mg_247[k];

        t_173[k] = -3.0 * kg_113[k]
                   + f_0 * mg_248[k];

        t_174[k] = -3.0 * kg_114[k]
                   + f_0 * mg_249[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, kg_115, kg_116, kg_117, kg_118, \
                         kg_119, mg_250, mg_251, mg_252, mg_253, \
                         mg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * kg_115[k]
                   + f_0 * mg_250[k];

        t_176[k] = -3.0 * kg_116[k]
                   + f_0 * mg_251[k];

        t_177[k] = -3.0 * kg_117[k]
                   + f_0 * mg_252[k];

        t_178[k] = -3.0 * kg_118[k]
                   + f_0 * mg_253[k];

        t_179[k] = -3.0 * kg_119[k]
                   + f_0 * mg_254[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, kg_120, kg_121, kg_122, kg_123, \
                         kg_124, mg_255, mg_256, mg_257, mg_258, \
                         mg_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * kg_120[k]
                   + f_0 * mg_255[k];

        t_181[k] = -2.0 * kg_121[k]
                   + f_0 * mg_256[k];

        t_182[k] = -2.0 * kg_122[k]
                   + f_0 * mg_257[k];

        t_183[k] = -2.0 * kg_123[k]
                   + f_0 * mg_258[k];

        t_184[k] = -2.0 * kg_124[k]
                   + f_0 * mg_259[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, kg_125, kg_126, kg_127, kg_128, \
                         kg_129, mg_260, mg_261, mg_262, mg_263, \
                         mg_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * kg_125[k]
                   + f_0 * mg_260[k];

        t_186[k] = -2.0 * kg_126[k]
                   + f_0 * mg_261[k];

        t_187[k] = -2.0 * kg_127[k]
                   + f_0 * mg_262[k];

        t_188[k] = -2.0 * kg_128[k]
                   + f_0 * mg_263[k];

        t_189[k] = -2.0 * kg_129[k]
                   + f_0 * mg_264[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, kg_130, kg_131, kg_132, kg_133, \
                         kg_134, mg_265, mg_266, mg_267, mg_268, \
                         mg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * kg_130[k]
                   + f_0 * mg_265[k];

        t_191[k] = -2.0 * kg_131[k]
                   + f_0 * mg_266[k];

        t_192[k] = -2.0 * kg_132[k]
                   + f_0 * mg_267[k];

        t_193[k] = -2.0 * kg_133[k]
                   + f_0 * mg_268[k];

        t_194[k] = -2.0 * kg_134[k]
                   + f_0 * mg_269[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, kg_135, kg_136, kg_137, kg_138, \
                         kg_139, mg_270, mg_271, mg_272, mg_273, \
                         mg_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -kg_135[k]
                   + f_0 * mg_270[k];

        t_196[k] = -kg_136[k]
                   + f_0 * mg_271[k];

        t_197[k] = -kg_137[k]
                   + f_0 * mg_272[k];

        t_198[k] = -kg_138[k]
                   + f_0 * mg_273[k];

        t_199[k] = -kg_139[k]
                   + f_0 * mg_274[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, kg_140, kg_141, kg_142, kg_143, \
                         kg_144, mg_275, mg_276, mg_277, mg_278, \
                         mg_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -kg_140[k]
                   + f_0 * mg_275[k];

        t_201[k] = -kg_141[k]
                   + f_0 * mg_276[k];

        t_202[k] = -kg_142[k]
                   + f_0 * mg_277[k];

        t_203[k] = -kg_143[k]
                   + f_0 * mg_278[k];

        t_204[k] = -kg_144[k]
                   + f_0 * mg_279[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, kg_145, kg_146, kg_147, kg_148, \
                         kg_149, mg_280, mg_281, mg_282, mg_283, \
                         mg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -kg_145[k]
                   + f_0 * mg_280[k];

        t_206[k] = -kg_146[k]
                   + f_0 * mg_281[k];

        t_207[k] = -kg_147[k]
                   + f_0 * mg_282[k];

        t_208[k] = -kg_148[k]
                   + f_0 * mg_283[k];

        t_209[k] = -kg_149[k]
                   + f_0 * mg_284[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, t_217, mg_285, \
                         mg_286, mg_287, mg_288, mg_289, mg_290, mg_291, \
                         mg_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * mg_285[k];

        t_211[k] = f_0 * mg_286[k];

        t_212[k] = f_0 * mg_287[k];

        t_213[k] = f_0 * mg_288[k];

        t_214[k] = f_0 * mg_289[k];

        t_215[k] = f_0 * mg_290[k];

        t_216[k] = f_0 * mg_291[k];

        t_217[k] = f_0 * mg_292[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, t_224, mg_293, mg_294, \
                         mg_295, mg_296, mg_297, mg_298, mg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * mg_293[k];

        t_219[k] = f_0 * mg_294[k];

        t_220[k] = f_0 * mg_295[k];

        t_221[k] = f_0 * mg_296[k];

        t_222[k] = f_0 * mg_297[k];

        t_223[k] = f_0 * mg_298[k];

        t_224[k] = f_0 * mg_299[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, kg_150, kg_151, kg_152, kg_153, \
                         kg_154, mg_315, mg_316, mg_317, mg_318, \
                         mg_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -5.0 * kg_150[k]
                   + f_0 * mg_315[k];

        t_226[k] = -5.0 * kg_151[k]
                   + f_0 * mg_316[k];

        t_227[k] = -5.0 * kg_152[k]
                   + f_0 * mg_317[k];

        t_228[k] = -5.0 * kg_153[k]
                   + f_0 * mg_318[k];

        t_229[k] = -5.0 * kg_154[k]
                   + f_0 * mg_319[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, kg_155, kg_156, kg_157, kg_158, \
                         kg_159, mg_320, mg_321, mg_322, mg_323, \
                         mg_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -5.0 * kg_155[k]
                   + f_0 * mg_320[k];

        t_231[k] = -5.0 * kg_156[k]
                   + f_0 * mg_321[k];

        t_232[k] = -5.0 * kg_157[k]
                   + f_0 * mg_322[k];

        t_233[k] = -5.0 * kg_158[k]
                   + f_0 * mg_323[k];

        t_234[k] = -5.0 * kg_159[k]
                   + f_0 * mg_324[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, kg_160, kg_161, kg_162, kg_163, \
                         kg_164, mg_325, mg_326, mg_327, mg_328, \
                         mg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -5.0 * kg_160[k]
                   + f_0 * mg_325[k];

        t_236[k] = -5.0 * kg_161[k]
                   + f_0 * mg_326[k];

        t_237[k] = -5.0 * kg_162[k]
                   + f_0 * mg_327[k];

        t_238[k] = -5.0 * kg_163[k]
                   + f_0 * mg_328[k];

        t_239[k] = -5.0 * kg_164[k]
                   + f_0 * mg_329[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, kg_165, kg_166, kg_167, kg_168, \
                         kg_169, mg_330, mg_331, mg_332, mg_333, \
                         mg_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -4.0 * kg_165[k]
                   + f_0 * mg_330[k];

        t_241[k] = -4.0 * kg_166[k]
                   + f_0 * mg_331[k];

        t_242[k] = -4.0 * kg_167[k]
                   + f_0 * mg_332[k];

        t_243[k] = -4.0 * kg_168[k]
                   + f_0 * mg_333[k];

        t_244[k] = -4.0 * kg_169[k]
                   + f_0 * mg_334[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, kg_170, kg_171, kg_172, kg_173, \
                         kg_174, mg_335, mg_336, mg_337, mg_338, \
                         mg_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -4.0 * kg_170[k]
                   + f_0 * mg_335[k];

        t_246[k] = -4.0 * kg_171[k]
                   + f_0 * mg_336[k];

        t_247[k] = -4.0 * kg_172[k]
                   + f_0 * mg_337[k];

        t_248[k] = -4.0 * kg_173[k]
                   + f_0 * mg_338[k];

        t_249[k] = -4.0 * kg_174[k]
                   + f_0 * mg_339[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, kg_175, kg_176, kg_177, kg_178, \
                         kg_179, mg_340, mg_341, mg_342, mg_343, \
                         mg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -4.0 * kg_175[k]
                   + f_0 * mg_340[k];

        t_251[k] = -4.0 * kg_176[k]
                   + f_0 * mg_341[k];

        t_252[k] = -4.0 * kg_177[k]
                   + f_0 * mg_342[k];

        t_253[k] = -4.0 * kg_178[k]
                   + f_0 * mg_343[k];

        t_254[k] = -4.0 * kg_179[k]
                   + f_0 * mg_344[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, kg_180, kg_181, kg_182, kg_183, \
                         kg_184, mg_345, mg_346, mg_347, mg_348, \
                         mg_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -3.0 * kg_180[k]
                   + f_0 * mg_345[k];

        t_256[k] = -3.0 * kg_181[k]
                   + f_0 * mg_346[k];

        t_257[k] = -3.0 * kg_182[k]
                   + f_0 * mg_347[k];

        t_258[k] = -3.0 * kg_183[k]
                   + f_0 * mg_348[k];

        t_259[k] = -3.0 * kg_184[k]
                   + f_0 * mg_349[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, kg_185, kg_186, kg_187, kg_188, \
                         kg_189, mg_350, mg_351, mg_352, mg_353, \
                         mg_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -3.0 * kg_185[k]
                   + f_0 * mg_350[k];

        t_261[k] = -3.0 * kg_186[k]
                   + f_0 * mg_351[k];

        t_262[k] = -3.0 * kg_187[k]
                   + f_0 * mg_352[k];

        t_263[k] = -3.0 * kg_188[k]
                   + f_0 * mg_353[k];

        t_264[k] = -3.0 * kg_189[k]
                   + f_0 * mg_354[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, kg_190, kg_191, kg_192, kg_193, \
                         kg_194, mg_355, mg_356, mg_357, mg_358, \
                         mg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -3.0 * kg_190[k]
                   + f_0 * mg_355[k];

        t_266[k] = -3.0 * kg_191[k]
                   + f_0 * mg_356[k];

        t_267[k] = -3.0 * kg_192[k]
                   + f_0 * mg_357[k];

        t_268[k] = -3.0 * kg_193[k]
                   + f_0 * mg_358[k];

        t_269[k] = -3.0 * kg_194[k]
                   + f_0 * mg_359[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, kg_195, kg_196, kg_197, kg_198, \
                         kg_199, mg_360, mg_361, mg_362, mg_363, \
                         mg_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * kg_195[k]
                   + f_0 * mg_360[k];

        t_271[k] = -2.0 * kg_196[k]
                   + f_0 * mg_361[k];

        t_272[k] = -2.0 * kg_197[k]
                   + f_0 * mg_362[k];

        t_273[k] = -2.0 * kg_198[k]
                   + f_0 * mg_363[k];

        t_274[k] = -2.0 * kg_199[k]
                   + f_0 * mg_364[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, kg_200, kg_201, kg_202, kg_203, \
                         kg_204, mg_365, mg_366, mg_367, mg_368, \
                         mg_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * kg_200[k]
                   + f_0 * mg_365[k];

        t_276[k] = -2.0 * kg_201[k]
                   + f_0 * mg_366[k];

        t_277[k] = -2.0 * kg_202[k]
                   + f_0 * mg_367[k];

        t_278[k] = -2.0 * kg_203[k]
                   + f_0 * mg_368[k];

        t_279[k] = -2.0 * kg_204[k]
                   + f_0 * mg_369[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, kg_205, kg_206, kg_207, kg_208, \
                         kg_209, mg_370, mg_371, mg_372, mg_373, \
                         mg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * kg_205[k]
                   + f_0 * mg_370[k];

        t_281[k] = -2.0 * kg_206[k]
                   + f_0 * mg_371[k];

        t_282[k] = -2.0 * kg_207[k]
                   + f_0 * mg_372[k];

        t_283[k] = -2.0 * kg_208[k]
                   + f_0 * mg_373[k];

        t_284[k] = -2.0 * kg_209[k]
                   + f_0 * mg_374[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, kg_210, kg_211, kg_212, kg_213, \
                         kg_214, mg_375, mg_376, mg_377, mg_378, \
                         mg_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -kg_210[k]
                   + f_0 * mg_375[k];

        t_286[k] = -kg_211[k]
                   + f_0 * mg_376[k];

        t_287[k] = -kg_212[k]
                   + f_0 * mg_377[k];

        t_288[k] = -kg_213[k]
                   + f_0 * mg_378[k];

        t_289[k] = -kg_214[k]
                   + f_0 * mg_379[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, kg_215, kg_216, kg_217, kg_218, \
                         kg_219, mg_380, mg_381, mg_382, mg_383, \
                         mg_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -kg_215[k]
                   + f_0 * mg_380[k];

        t_291[k] = -kg_216[k]
                   + f_0 * mg_381[k];

        t_292[k] = -kg_217[k]
                   + f_0 * mg_382[k];

        t_293[k] = -kg_218[k]
                   + f_0 * mg_383[k];

        t_294[k] = -kg_219[k]
                   + f_0 * mg_384[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, kg_220, kg_221, kg_222, kg_223, \
                         kg_224, mg_385, mg_386, mg_387, mg_388, \
                         mg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -kg_220[k]
                   + f_0 * mg_385[k];

        t_296[k] = -kg_221[k]
                   + f_0 * mg_386[k];

        t_297[k] = -kg_222[k]
                   + f_0 * mg_387[k];

        t_298[k] = -kg_223[k]
                   + f_0 * mg_388[k];

        t_299[k] = -kg_224[k]
                   + f_0 * mg_389[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, t_305, t_306, t_307, mg_390, \
                         mg_391, mg_392, mg_393, mg_394, mg_395, mg_396, \
                         mg_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_0 * mg_390[k];

        t_301[k] = f_0 * mg_391[k];

        t_302[k] = f_0 * mg_392[k];

        t_303[k] = f_0 * mg_393[k];

        t_304[k] = f_0 * mg_394[k];

        t_305[k] = f_0 * mg_395[k];

        t_306[k] = f_0 * mg_396[k];

        t_307[k] = f_0 * mg_397[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, t_313, t_314, mg_398, mg_399, \
                         mg_400, mg_401, mg_402, mg_403, mg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_0 * mg_398[k];

        t_309[k] = f_0 * mg_399[k];

        t_310[k] = f_0 * mg_400[k];

        t_311[k] = f_0 * mg_401[k];

        t_312[k] = f_0 * mg_402[k];

        t_313[k] = f_0 * mg_403[k];

        t_314[k] = f_0 * mg_404[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, kg_225, kg_226, kg_227, kg_228, \
                         kg_229, mg_420, mg_421, mg_422, mg_423, \
                         mg_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -6.0 * kg_225[k]
                   + f_0 * mg_420[k];

        t_316[k] = -6.0 * kg_226[k]
                   + f_0 * mg_421[k];

        t_317[k] = -6.0 * kg_227[k]
                   + f_0 * mg_422[k];

        t_318[k] = -6.0 * kg_228[k]
                   + f_0 * mg_423[k];

        t_319[k] = -6.0 * kg_229[k]
                   + f_0 * mg_424[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, kg_230, kg_231, kg_232, kg_233, \
                         kg_234, mg_425, mg_426, mg_427, mg_428, \
                         mg_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -6.0 * kg_230[k]
                   + f_0 * mg_425[k];

        t_321[k] = -6.0 * kg_231[k]
                   + f_0 * mg_426[k];

        t_322[k] = -6.0 * kg_232[k]
                   + f_0 * mg_427[k];

        t_323[k] = -6.0 * kg_233[k]
                   + f_0 * mg_428[k];

        t_324[k] = -6.0 * kg_234[k]
                   + f_0 * mg_429[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, kg_235, kg_236, kg_237, kg_238, \
                         kg_239, mg_430, mg_431, mg_432, mg_433, \
                         mg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -6.0 * kg_235[k]
                   + f_0 * mg_430[k];

        t_326[k] = -6.0 * kg_236[k]
                   + f_0 * mg_431[k];

        t_327[k] = -6.0 * kg_237[k]
                   + f_0 * mg_432[k];

        t_328[k] = -6.0 * kg_238[k]
                   + f_0 * mg_433[k];

        t_329[k] = -6.0 * kg_239[k]
                   + f_0 * mg_434[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
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

    const auto *mg_435 = buffer.data(mg + 435);
    const auto *mg_436 = buffer.data(mg + 436);
    const auto *mg_437 = buffer.data(mg + 437);
    const auto *mg_438 = buffer.data(mg + 438);
    const auto *mg_439 = buffer.data(mg + 439);
    const auto *mg_440 = buffer.data(mg + 440);
    const auto *mg_441 = buffer.data(mg + 441);
    const auto *mg_442 = buffer.data(mg + 442);
    const auto *mg_443 = buffer.data(mg + 443);
    const auto *mg_444 = buffer.data(mg + 444);
    const auto *mg_445 = buffer.data(mg + 445);
    const auto *mg_446 = buffer.data(mg + 446);
    const auto *mg_447 = buffer.data(mg + 447);
    const auto *mg_448 = buffer.data(mg + 448);
    const auto *mg_449 = buffer.data(mg + 449);
    const auto *mg_450 = buffer.data(mg + 450);
    const auto *mg_451 = buffer.data(mg + 451);
    const auto *mg_452 = buffer.data(mg + 452);
    const auto *mg_453 = buffer.data(mg + 453);
    const auto *mg_454 = buffer.data(mg + 454);
    const auto *mg_455 = buffer.data(mg + 455);
    const auto *mg_456 = buffer.data(mg + 456);
    const auto *mg_457 = buffer.data(mg + 457);
    const auto *mg_458 = buffer.data(mg + 458);
    const auto *mg_459 = buffer.data(mg + 459);
    const auto *mg_460 = buffer.data(mg + 460);
    const auto *mg_461 = buffer.data(mg + 461);
    const auto *mg_462 = buffer.data(mg + 462);
    const auto *mg_463 = buffer.data(mg + 463);
    const auto *mg_464 = buffer.data(mg + 464);
    const auto *mg_465 = buffer.data(mg + 465);
    const auto *mg_466 = buffer.data(mg + 466);
    const auto *mg_467 = buffer.data(mg + 467);
    const auto *mg_468 = buffer.data(mg + 468);
    const auto *mg_469 = buffer.data(mg + 469);
    const auto *mg_470 = buffer.data(mg + 470);
    const auto *mg_471 = buffer.data(mg + 471);
    const auto *mg_472 = buffer.data(mg + 472);
    const auto *mg_473 = buffer.data(mg + 473);
    const auto *mg_474 = buffer.data(mg + 474);
    const auto *mg_475 = buffer.data(mg + 475);
    const auto *mg_476 = buffer.data(mg + 476);
    const auto *mg_477 = buffer.data(mg + 477);
    const auto *mg_478 = buffer.data(mg + 478);
    const auto *mg_479 = buffer.data(mg + 479);
    const auto *mg_480 = buffer.data(mg + 480);
    const auto *mg_481 = buffer.data(mg + 481);
    const auto *mg_482 = buffer.data(mg + 482);
    const auto *mg_483 = buffer.data(mg + 483);
    const auto *mg_484 = buffer.data(mg + 484);
    const auto *mg_485 = buffer.data(mg + 485);
    const auto *mg_486 = buffer.data(mg + 486);
    const auto *mg_487 = buffer.data(mg + 487);
    const auto *mg_488 = buffer.data(mg + 488);
    const auto *mg_489 = buffer.data(mg + 489);
    const auto *mg_490 = buffer.data(mg + 490);
    const auto *mg_491 = buffer.data(mg + 491);
    const auto *mg_492 = buffer.data(mg + 492);
    const auto *mg_493 = buffer.data(mg + 493);
    const auto *mg_494 = buffer.data(mg + 494);
    const auto *mg_495 = buffer.data(mg + 495);
    const auto *mg_496 = buffer.data(mg + 496);
    const auto *mg_497 = buffer.data(mg + 497);
    const auto *mg_498 = buffer.data(mg + 498);
    const auto *mg_499 = buffer.data(mg + 499);
    const auto *mg_500 = buffer.data(mg + 500);
    const auto *mg_501 = buffer.data(mg + 501);
    const auto *mg_502 = buffer.data(mg + 502);
    const auto *mg_503 = buffer.data(mg + 503);
    const auto *mg_504 = buffer.data(mg + 504);
    const auto *mg_505 = buffer.data(mg + 505);
    const auto *mg_506 = buffer.data(mg + 506);
    const auto *mg_507 = buffer.data(mg + 507);
    const auto *mg_508 = buffer.data(mg + 508);
    const auto *mg_509 = buffer.data(mg + 509);
    const auto *mg_510 = buffer.data(mg + 510);
    const auto *mg_511 = buffer.data(mg + 511);
    const auto *mg_512 = buffer.data(mg + 512);
    const auto *mg_513 = buffer.data(mg + 513);
    const auto *mg_514 = buffer.data(mg + 514);
    const auto *mg_515 = buffer.data(mg + 515);
    const auto *mg_516 = buffer.data(mg + 516);
    const auto *mg_517 = buffer.data(mg + 517);
    const auto *mg_518 = buffer.data(mg + 518);
    const auto *mg_519 = buffer.data(mg + 519);
    const auto *mg_520 = buffer.data(mg + 520);
    const auto *mg_521 = buffer.data(mg + 521);
    const auto *mg_522 = buffer.data(mg + 522);
    const auto *mg_523 = buffer.data(mg + 523);
    const auto *mg_524 = buffer.data(mg + 524);
    const auto *mg_540 = buffer.data(mg + 540);
    const auto *mg_541 = buffer.data(mg + 541);
    const auto *mg_542 = buffer.data(mg + 542);
    const auto *mg_543 = buffer.data(mg + 543);
    const auto *mg_544 = buffer.data(mg + 544);
    const auto *mg_545 = buffer.data(mg + 545);
    const auto *mg_546 = buffer.data(mg + 546);
    const auto *mg_547 = buffer.data(mg + 547);
    const auto *mg_548 = buffer.data(mg + 548);
    const auto *mg_549 = buffer.data(mg + 549);
    const auto *mg_550 = buffer.data(mg + 550);
    const auto *mg_551 = buffer.data(mg + 551);
    const auto *mg_552 = buffer.data(mg + 552);
    const auto *mg_553 = buffer.data(mg + 553);
    const auto *mg_554 = buffer.data(mg + 554);
    const auto *mg_555 = buffer.data(mg + 555);
    const auto *mg_556 = buffer.data(mg + 556);
    const auto *mg_557 = buffer.data(mg + 557);
    const auto *mg_558 = buffer.data(mg + 558);
    const auto *mg_559 = buffer.data(mg + 559);
    const auto *mg_560 = buffer.data(mg + 560);
    const auto *mg_561 = buffer.data(mg + 561);
    const auto *mg_562 = buffer.data(mg + 562);
    const auto *mg_563 = buffer.data(mg + 563);
    const auto *mg_564 = buffer.data(mg + 564);
    const auto *mg_565 = buffer.data(mg + 565);
    const auto *mg_566 = buffer.data(mg + 566);
    const auto *mg_567 = buffer.data(mg + 567);
    const auto *mg_568 = buffer.data(mg + 568);
    const auto *mg_569 = buffer.data(mg + 569);
    const auto *mg_570 = buffer.data(mg + 570);
    const auto *mg_571 = buffer.data(mg + 571);
    const auto *mg_572 = buffer.data(mg + 572);
    const auto *mg_573 = buffer.data(mg + 573);
    const auto *mg_574 = buffer.data(mg + 574);
    const auto *mg_575 = buffer.data(mg + 575);
    const auto *mg_576 = buffer.data(mg + 576);
    const auto *mg_577 = buffer.data(mg + 577);
    const auto *mg_578 = buffer.data(mg + 578);
    const auto *mg_579 = buffer.data(mg + 579);
    const auto *mg_580 = buffer.data(mg + 580);
    const auto *mg_581 = buffer.data(mg + 581);
    const auto *mg_582 = buffer.data(mg + 582);
    const auto *mg_583 = buffer.data(mg + 583);
    const auto *mg_584 = buffer.data(mg + 584);
    const auto *mg_585 = buffer.data(mg + 585);
    const auto *mg_586 = buffer.data(mg + 586);
    const auto *mg_587 = buffer.data(mg + 587);
    const auto *mg_588 = buffer.data(mg + 588);
    const auto *mg_589 = buffer.data(mg + 589);
    const auto *mg_590 = buffer.data(mg + 590);
    const auto *mg_591 = buffer.data(mg + 591);
    const auto *mg_592 = buffer.data(mg + 592);
    const auto *mg_593 = buffer.data(mg + 593);
    const auto *mg_594 = buffer.data(mg + 594);
    const auto *mg_595 = buffer.data(mg + 595);
    const auto *mg_596 = buffer.data(mg + 596);
    const auto *mg_597 = buffer.data(mg + 597);
    const auto *mg_598 = buffer.data(mg + 598);
    const auto *mg_599 = buffer.data(mg + 599);
    const auto *mg_600 = buffer.data(mg + 600);
    const auto *mg_601 = buffer.data(mg + 601);
    const auto *mg_602 = buffer.data(mg + 602);
    const auto *mg_603 = buffer.data(mg + 603);
    const auto *mg_604 = buffer.data(mg + 604);

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, kg_240, kg_241, kg_242, kg_243, \
                         kg_244, mg_435, mg_436, mg_437, mg_438, \
                         mg_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -5.0 * kg_240[k]
                   + f_0 * mg_435[k];

        t_331[k] = -5.0 * kg_241[k]
                   + f_0 * mg_436[k];

        t_332[k] = -5.0 * kg_242[k]
                   + f_0 * mg_437[k];

        t_333[k] = -5.0 * kg_243[k]
                   + f_0 * mg_438[k];

        t_334[k] = -5.0 * kg_244[k]
                   + f_0 * mg_439[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, kg_245, kg_246, kg_247, kg_248, \
                         kg_249, mg_440, mg_441, mg_442, mg_443, \
                         mg_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -5.0 * kg_245[k]
                   + f_0 * mg_440[k];

        t_336[k] = -5.0 * kg_246[k]
                   + f_0 * mg_441[k];

        t_337[k] = -5.0 * kg_247[k]
                   + f_0 * mg_442[k];

        t_338[k] = -5.0 * kg_248[k]
                   + f_0 * mg_443[k];

        t_339[k] = -5.0 * kg_249[k]
                   + f_0 * mg_444[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, kg_250, kg_251, kg_252, kg_253, \
                         kg_254, mg_445, mg_446, mg_447, mg_448, \
                         mg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -5.0 * kg_250[k]
                   + f_0 * mg_445[k];

        t_341[k] = -5.0 * kg_251[k]
                   + f_0 * mg_446[k];

        t_342[k] = -5.0 * kg_252[k]
                   + f_0 * mg_447[k];

        t_343[k] = -5.0 * kg_253[k]
                   + f_0 * mg_448[k];

        t_344[k] = -5.0 * kg_254[k]
                   + f_0 * mg_449[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, kg_255, kg_256, kg_257, kg_258, \
                         kg_259, mg_450, mg_451, mg_452, mg_453, \
                         mg_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -4.0 * kg_255[k]
                   + f_0 * mg_450[k];

        t_346[k] = -4.0 * kg_256[k]
                   + f_0 * mg_451[k];

        t_347[k] = -4.0 * kg_257[k]
                   + f_0 * mg_452[k];

        t_348[k] = -4.0 * kg_258[k]
                   + f_0 * mg_453[k];

        t_349[k] = -4.0 * kg_259[k]
                   + f_0 * mg_454[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, kg_260, kg_261, kg_262, kg_263, \
                         kg_264, mg_455, mg_456, mg_457, mg_458, \
                         mg_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -4.0 * kg_260[k]
                   + f_0 * mg_455[k];

        t_351[k] = -4.0 * kg_261[k]
                   + f_0 * mg_456[k];

        t_352[k] = -4.0 * kg_262[k]
                   + f_0 * mg_457[k];

        t_353[k] = -4.0 * kg_263[k]
                   + f_0 * mg_458[k];

        t_354[k] = -4.0 * kg_264[k]
                   + f_0 * mg_459[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, kg_265, kg_266, kg_267, kg_268, \
                         kg_269, mg_460, mg_461, mg_462, mg_463, \
                         mg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -4.0 * kg_265[k]
                   + f_0 * mg_460[k];

        t_356[k] = -4.0 * kg_266[k]
                   + f_0 * mg_461[k];

        t_357[k] = -4.0 * kg_267[k]
                   + f_0 * mg_462[k];

        t_358[k] = -4.0 * kg_268[k]
                   + f_0 * mg_463[k];

        t_359[k] = -4.0 * kg_269[k]
                   + f_0 * mg_464[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, kg_270, kg_271, kg_272, kg_273, \
                         kg_274, mg_465, mg_466, mg_467, mg_468, \
                         mg_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -3.0 * kg_270[k]
                   + f_0 * mg_465[k];

        t_361[k] = -3.0 * kg_271[k]
                   + f_0 * mg_466[k];

        t_362[k] = -3.0 * kg_272[k]
                   + f_0 * mg_467[k];

        t_363[k] = -3.0 * kg_273[k]
                   + f_0 * mg_468[k];

        t_364[k] = -3.0 * kg_274[k]
                   + f_0 * mg_469[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, kg_275, kg_276, kg_277, kg_278, \
                         kg_279, mg_470, mg_471, mg_472, mg_473, \
                         mg_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -3.0 * kg_275[k]
                   + f_0 * mg_470[k];

        t_366[k] = -3.0 * kg_276[k]
                   + f_0 * mg_471[k];

        t_367[k] = -3.0 * kg_277[k]
                   + f_0 * mg_472[k];

        t_368[k] = -3.0 * kg_278[k]
                   + f_0 * mg_473[k];

        t_369[k] = -3.0 * kg_279[k]
                   + f_0 * mg_474[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, kg_280, kg_281, kg_282, kg_283, \
                         kg_284, mg_475, mg_476, mg_477, mg_478, \
                         mg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -3.0 * kg_280[k]
                   + f_0 * mg_475[k];

        t_371[k] = -3.0 * kg_281[k]
                   + f_0 * mg_476[k];

        t_372[k] = -3.0 * kg_282[k]
                   + f_0 * mg_477[k];

        t_373[k] = -3.0 * kg_283[k]
                   + f_0 * mg_478[k];

        t_374[k] = -3.0 * kg_284[k]
                   + f_0 * mg_479[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, kg_285, kg_286, kg_287, kg_288, \
                         kg_289, mg_480, mg_481, mg_482, mg_483, \
                         mg_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -2.0 * kg_285[k]
                   + f_0 * mg_480[k];

        t_376[k] = -2.0 * kg_286[k]
                   + f_0 * mg_481[k];

        t_377[k] = -2.0 * kg_287[k]
                   + f_0 * mg_482[k];

        t_378[k] = -2.0 * kg_288[k]
                   + f_0 * mg_483[k];

        t_379[k] = -2.0 * kg_289[k]
                   + f_0 * mg_484[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, kg_290, kg_291, kg_292, kg_293, \
                         kg_294, mg_485, mg_486, mg_487, mg_488, \
                         mg_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * kg_290[k]
                   + f_0 * mg_485[k];

        t_381[k] = -2.0 * kg_291[k]
                   + f_0 * mg_486[k];

        t_382[k] = -2.0 * kg_292[k]
                   + f_0 * mg_487[k];

        t_383[k] = -2.0 * kg_293[k]
                   + f_0 * mg_488[k];

        t_384[k] = -2.0 * kg_294[k]
                   + f_0 * mg_489[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, kg_295, kg_296, kg_297, kg_298, \
                         kg_299, mg_490, mg_491, mg_492, mg_493, \
                         mg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * kg_295[k]
                   + f_0 * mg_490[k];

        t_386[k] = -2.0 * kg_296[k]
                   + f_0 * mg_491[k];

        t_387[k] = -2.0 * kg_297[k]
                   + f_0 * mg_492[k];

        t_388[k] = -2.0 * kg_298[k]
                   + f_0 * mg_493[k];

        t_389[k] = -2.0 * kg_299[k]
                   + f_0 * mg_494[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, kg_300, kg_301, kg_302, kg_303, \
                         kg_304, mg_495, mg_496, mg_497, mg_498, \
                         mg_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -kg_300[k]
                   + f_0 * mg_495[k];

        t_391[k] = -kg_301[k]
                   + f_0 * mg_496[k];

        t_392[k] = -kg_302[k]
                   + f_0 * mg_497[k];

        t_393[k] = -kg_303[k]
                   + f_0 * mg_498[k];

        t_394[k] = -kg_304[k]
                   + f_0 * mg_499[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, kg_305, kg_306, kg_307, kg_308, \
                         kg_309, mg_500, mg_501, mg_502, mg_503, \
                         mg_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -kg_305[k]
                   + f_0 * mg_500[k];

        t_396[k] = -kg_306[k]
                   + f_0 * mg_501[k];

        t_397[k] = -kg_307[k]
                   + f_0 * mg_502[k];

        t_398[k] = -kg_308[k]
                   + f_0 * mg_503[k];

        t_399[k] = -kg_309[k]
                   + f_0 * mg_504[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, kg_310, kg_311, kg_312, kg_313, \
                         kg_314, mg_505, mg_506, mg_507, mg_508, \
                         mg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -kg_310[k]
                   + f_0 * mg_505[k];

        t_401[k] = -kg_311[k]
                   + f_0 * mg_506[k];

        t_402[k] = -kg_312[k]
                   + f_0 * mg_507[k];

        t_403[k] = -kg_313[k]
                   + f_0 * mg_508[k];

        t_404[k] = -kg_314[k]
                   + f_0 * mg_509[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, mg_510, \
                         mg_511, mg_512, mg_513, mg_514, mg_515, mg_516, \
                         mg_517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * mg_510[k];

        t_406[k] = f_0 * mg_511[k];

        t_407[k] = f_0 * mg_512[k];

        t_408[k] = f_0 * mg_513[k];

        t_409[k] = f_0 * mg_514[k];

        t_410[k] = f_0 * mg_515[k];

        t_411[k] = f_0 * mg_516[k];

        t_412[k] = f_0 * mg_517[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, mg_518, mg_519, \
                         mg_520, mg_521, mg_522, mg_523, mg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * mg_518[k];

        t_414[k] = f_0 * mg_519[k];

        t_415[k] = f_0 * mg_520[k];

        t_416[k] = f_0 * mg_521[k];

        t_417[k] = f_0 * mg_522[k];

        t_418[k] = f_0 * mg_523[k];

        t_419[k] = f_0 * mg_524[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, kg_315, kg_316, kg_317, kg_318, \
                         kg_319, mg_540, mg_541, mg_542, mg_543, \
                         mg_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -7.0 * kg_315[k]
                   + f_0 * mg_540[k];

        t_421[k] = -7.0 * kg_316[k]
                   + f_0 * mg_541[k];

        t_422[k] = -7.0 * kg_317[k]
                   + f_0 * mg_542[k];

        t_423[k] = -7.0 * kg_318[k]
                   + f_0 * mg_543[k];

        t_424[k] = -7.0 * kg_319[k]
                   + f_0 * mg_544[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, kg_320, kg_321, kg_322, kg_323, \
                         kg_324, mg_545, mg_546, mg_547, mg_548, \
                         mg_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -7.0 * kg_320[k]
                   + f_0 * mg_545[k];

        t_426[k] = -7.0 * kg_321[k]
                   + f_0 * mg_546[k];

        t_427[k] = -7.0 * kg_322[k]
                   + f_0 * mg_547[k];

        t_428[k] = -7.0 * kg_323[k]
                   + f_0 * mg_548[k];

        t_429[k] = -7.0 * kg_324[k]
                   + f_0 * mg_549[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, kg_325, kg_326, kg_327, kg_328, \
                         kg_329, mg_550, mg_551, mg_552, mg_553, \
                         mg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -7.0 * kg_325[k]
                   + f_0 * mg_550[k];

        t_431[k] = -7.0 * kg_326[k]
                   + f_0 * mg_551[k];

        t_432[k] = -7.0 * kg_327[k]
                   + f_0 * mg_552[k];

        t_433[k] = -7.0 * kg_328[k]
                   + f_0 * mg_553[k];

        t_434[k] = -7.0 * kg_329[k]
                   + f_0 * mg_554[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, kg_330, kg_331, kg_332, kg_333, \
                         kg_334, mg_555, mg_556, mg_557, mg_558, \
                         mg_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -6.0 * kg_330[k]
                   + f_0 * mg_555[k];

        t_436[k] = -6.0 * kg_331[k]
                   + f_0 * mg_556[k];

        t_437[k] = -6.0 * kg_332[k]
                   + f_0 * mg_557[k];

        t_438[k] = -6.0 * kg_333[k]
                   + f_0 * mg_558[k];

        t_439[k] = -6.0 * kg_334[k]
                   + f_0 * mg_559[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, kg_335, kg_336, kg_337, kg_338, \
                         kg_339, mg_560, mg_561, mg_562, mg_563, \
                         mg_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -6.0 * kg_335[k]
                   + f_0 * mg_560[k];

        t_441[k] = -6.0 * kg_336[k]
                   + f_0 * mg_561[k];

        t_442[k] = -6.0 * kg_337[k]
                   + f_0 * mg_562[k];

        t_443[k] = -6.0 * kg_338[k]
                   + f_0 * mg_563[k];

        t_444[k] = -6.0 * kg_339[k]
                   + f_0 * mg_564[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, kg_340, kg_341, kg_342, kg_343, \
                         kg_344, mg_565, mg_566, mg_567, mg_568, \
                         mg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -6.0 * kg_340[k]
                   + f_0 * mg_565[k];

        t_446[k] = -6.0 * kg_341[k]
                   + f_0 * mg_566[k];

        t_447[k] = -6.0 * kg_342[k]
                   + f_0 * mg_567[k];

        t_448[k] = -6.0 * kg_343[k]
                   + f_0 * mg_568[k];

        t_449[k] = -6.0 * kg_344[k]
                   + f_0 * mg_569[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, kg_345, kg_346, kg_347, kg_348, \
                         kg_349, mg_570, mg_571, mg_572, mg_573, \
                         mg_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -5.0 * kg_345[k]
                   + f_0 * mg_570[k];

        t_451[k] = -5.0 * kg_346[k]
                   + f_0 * mg_571[k];

        t_452[k] = -5.0 * kg_347[k]
                   + f_0 * mg_572[k];

        t_453[k] = -5.0 * kg_348[k]
                   + f_0 * mg_573[k];

        t_454[k] = -5.0 * kg_349[k]
                   + f_0 * mg_574[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, kg_350, kg_351, kg_352, kg_353, \
                         kg_354, mg_575, mg_576, mg_577, mg_578, \
                         mg_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -5.0 * kg_350[k]
                   + f_0 * mg_575[k];

        t_456[k] = -5.0 * kg_351[k]
                   + f_0 * mg_576[k];

        t_457[k] = -5.0 * kg_352[k]
                   + f_0 * mg_577[k];

        t_458[k] = -5.0 * kg_353[k]
                   + f_0 * mg_578[k];

        t_459[k] = -5.0 * kg_354[k]
                   + f_0 * mg_579[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, kg_355, kg_356, kg_357, kg_358, \
                         kg_359, mg_580, mg_581, mg_582, mg_583, \
                         mg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -5.0 * kg_355[k]
                   + f_0 * mg_580[k];

        t_461[k] = -5.0 * kg_356[k]
                   + f_0 * mg_581[k];

        t_462[k] = -5.0 * kg_357[k]
                   + f_0 * mg_582[k];

        t_463[k] = -5.0 * kg_358[k]
                   + f_0 * mg_583[k];

        t_464[k] = -5.0 * kg_359[k]
                   + f_0 * mg_584[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, kg_360, kg_361, kg_362, kg_363, \
                         kg_364, mg_585, mg_586, mg_587, mg_588, \
                         mg_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -4.0 * kg_360[k]
                   + f_0 * mg_585[k];

        t_466[k] = -4.0 * kg_361[k]
                   + f_0 * mg_586[k];

        t_467[k] = -4.0 * kg_362[k]
                   + f_0 * mg_587[k];

        t_468[k] = -4.0 * kg_363[k]
                   + f_0 * mg_588[k];

        t_469[k] = -4.0 * kg_364[k]
                   + f_0 * mg_589[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, kg_365, kg_366, kg_367, kg_368, \
                         kg_369, mg_590, mg_591, mg_592, mg_593, \
                         mg_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -4.0 * kg_365[k]
                   + f_0 * mg_590[k];

        t_471[k] = -4.0 * kg_366[k]
                   + f_0 * mg_591[k];

        t_472[k] = -4.0 * kg_367[k]
                   + f_0 * mg_592[k];

        t_473[k] = -4.0 * kg_368[k]
                   + f_0 * mg_593[k];

        t_474[k] = -4.0 * kg_369[k]
                   + f_0 * mg_594[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, kg_370, kg_371, kg_372, kg_373, \
                         kg_374, mg_595, mg_596, mg_597, mg_598, \
                         mg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -4.0 * kg_370[k]
                   + f_0 * mg_595[k];

        t_476[k] = -4.0 * kg_371[k]
                   + f_0 * mg_596[k];

        t_477[k] = -4.0 * kg_372[k]
                   + f_0 * mg_597[k];

        t_478[k] = -4.0 * kg_373[k]
                   + f_0 * mg_598[k];

        t_479[k] = -4.0 * kg_374[k]
                   + f_0 * mg_599[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, kg_375, kg_376, kg_377, kg_378, \
                         kg_379, mg_600, mg_601, mg_602, mg_603, \
                         mg_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -3.0 * kg_375[k]
                   + f_0 * mg_600[k];

        t_481[k] = -3.0 * kg_376[k]
                   + f_0 * mg_601[k];

        t_482[k] = -3.0 * kg_377[k]
                   + f_0 * mg_602[k];

        t_483[k] = -3.0 * kg_378[k]
                   + f_0 * mg_603[k];

        t_484[k] = -3.0 * kg_379[k]
                   + f_0 * mg_604[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *mg_605 = buffer.data(mg + 605);
    const auto *mg_606 = buffer.data(mg + 606);
    const auto *mg_607 = buffer.data(mg + 607);
    const auto *mg_608 = buffer.data(mg + 608);
    const auto *mg_609 = buffer.data(mg + 609);
    const auto *mg_610 = buffer.data(mg + 610);
    const auto *mg_611 = buffer.data(mg + 611);
    const auto *mg_612 = buffer.data(mg + 612);
    const auto *mg_613 = buffer.data(mg + 613);
    const auto *mg_614 = buffer.data(mg + 614);
    const auto *mg_615 = buffer.data(mg + 615);
    const auto *mg_616 = buffer.data(mg + 616);
    const auto *mg_617 = buffer.data(mg + 617);
    const auto *mg_618 = buffer.data(mg + 618);
    const auto *mg_619 = buffer.data(mg + 619);
    const auto *mg_620 = buffer.data(mg + 620);
    const auto *mg_621 = buffer.data(mg + 621);
    const auto *mg_622 = buffer.data(mg + 622);
    const auto *mg_623 = buffer.data(mg + 623);
    const auto *mg_624 = buffer.data(mg + 624);
    const auto *mg_625 = buffer.data(mg + 625);
    const auto *mg_626 = buffer.data(mg + 626);
    const auto *mg_627 = buffer.data(mg + 627);
    const auto *mg_628 = buffer.data(mg + 628);
    const auto *mg_629 = buffer.data(mg + 629);
    const auto *mg_630 = buffer.data(mg + 630);
    const auto *mg_631 = buffer.data(mg + 631);
    const auto *mg_632 = buffer.data(mg + 632);
    const auto *mg_633 = buffer.data(mg + 633);
    const auto *mg_634 = buffer.data(mg + 634);
    const auto *mg_635 = buffer.data(mg + 635);
    const auto *mg_636 = buffer.data(mg + 636);
    const auto *mg_637 = buffer.data(mg + 637);
    const auto *mg_638 = buffer.data(mg + 638);
    const auto *mg_639 = buffer.data(mg + 639);
    const auto *mg_640 = buffer.data(mg + 640);
    const auto *mg_641 = buffer.data(mg + 641);
    const auto *mg_642 = buffer.data(mg + 642);
    const auto *mg_643 = buffer.data(mg + 643);
    const auto *mg_644 = buffer.data(mg + 644);
    const auto *mg_645 = buffer.data(mg + 645);
    const auto *mg_646 = buffer.data(mg + 646);
    const auto *mg_647 = buffer.data(mg + 647);
    const auto *mg_648 = buffer.data(mg + 648);
    const auto *mg_649 = buffer.data(mg + 649);
    const auto *mg_650 = buffer.data(mg + 650);
    const auto *mg_651 = buffer.data(mg + 651);
    const auto *mg_652 = buffer.data(mg + 652);
    const auto *mg_653 = buffer.data(mg + 653);
    const auto *mg_654 = buffer.data(mg + 654);
    const auto *mg_655 = buffer.data(mg + 655);
    const auto *mg_656 = buffer.data(mg + 656);
    const auto *mg_657 = buffer.data(mg + 657);
    const auto *mg_658 = buffer.data(mg + 658);
    const auto *mg_659 = buffer.data(mg + 659);
    const auto *mg_675 = buffer.data(mg + 675);
    const auto *mg_676 = buffer.data(mg + 676);
    const auto *mg_677 = buffer.data(mg + 677);
    const auto *mg_678 = buffer.data(mg + 678);
    const auto *mg_679 = buffer.data(mg + 679);
    const auto *mg_680 = buffer.data(mg + 680);
    const auto *mg_681 = buffer.data(mg + 681);
    const auto *mg_682 = buffer.data(mg + 682);
    const auto *mg_683 = buffer.data(mg + 683);
    const auto *mg_684 = buffer.data(mg + 684);
    const auto *mg_685 = buffer.data(mg + 685);
    const auto *mg_686 = buffer.data(mg + 686);
    const auto *mg_687 = buffer.data(mg + 687);
    const auto *mg_688 = buffer.data(mg + 688);
    const auto *mg_689 = buffer.data(mg + 689);
    const auto *mg_690 = buffer.data(mg + 690);
    const auto *mg_691 = buffer.data(mg + 691);
    const auto *mg_692 = buffer.data(mg + 692);
    const auto *mg_693 = buffer.data(mg + 693);
    const auto *mg_694 = buffer.data(mg + 694);
    const auto *mg_695 = buffer.data(mg + 695);
    const auto *mg_696 = buffer.data(mg + 696);
    const auto *mg_697 = buffer.data(mg + 697);
    const auto *mg_698 = buffer.data(mg + 698);
    const auto *mg_699 = buffer.data(mg + 699);
    const auto *mg_700 = buffer.data(mg + 700);
    const auto *mg_701 = buffer.data(mg + 701);
    const auto *mg_702 = buffer.data(mg + 702);
    const auto *mg_703 = buffer.data(mg + 703);
    const auto *mg_704 = buffer.data(mg + 704);
    const auto *mg_705 = buffer.data(mg + 705);
    const auto *mg_706 = buffer.data(mg + 706);
    const auto *mg_707 = buffer.data(mg + 707);
    const auto *mg_708 = buffer.data(mg + 708);
    const auto *mg_709 = buffer.data(mg + 709);
    const auto *mg_710 = buffer.data(mg + 710);
    const auto *mg_711 = buffer.data(mg + 711);
    const auto *mg_712 = buffer.data(mg + 712);
    const auto *mg_713 = buffer.data(mg + 713);
    const auto *mg_714 = buffer.data(mg + 714);
    const auto *mg_715 = buffer.data(mg + 715);
    const auto *mg_716 = buffer.data(mg + 716);
    const auto *mg_717 = buffer.data(mg + 717);
    const auto *mg_718 = buffer.data(mg + 718);
    const auto *mg_719 = buffer.data(mg + 719);
    const auto *mg_720 = buffer.data(mg + 720);
    const auto *mg_721 = buffer.data(mg + 721);
    const auto *mg_722 = buffer.data(mg + 722);
    const auto *mg_723 = buffer.data(mg + 723);
    const auto *mg_724 = buffer.data(mg + 724);
    const auto *mg_725 = buffer.data(mg + 725);
    const auto *mg_726 = buffer.data(mg + 726);
    const auto *mg_727 = buffer.data(mg + 727);
    const auto *mg_728 = buffer.data(mg + 728);
    const auto *mg_729 = buffer.data(mg + 729);
    const auto *mg_730 = buffer.data(mg + 730);
    const auto *mg_731 = buffer.data(mg + 731);
    const auto *mg_732 = buffer.data(mg + 732);
    const auto *mg_733 = buffer.data(mg + 733);
    const auto *mg_734 = buffer.data(mg + 734);
    const auto *mg_735 = buffer.data(mg + 735);
    const auto *mg_736 = buffer.data(mg + 736);
    const auto *mg_737 = buffer.data(mg + 737);
    const auto *mg_738 = buffer.data(mg + 738);
    const auto *mg_739 = buffer.data(mg + 739);
    const auto *mg_740 = buffer.data(mg + 740);
    const auto *mg_741 = buffer.data(mg + 741);
    const auto *mg_742 = buffer.data(mg + 742);
    const auto *mg_743 = buffer.data(mg + 743);
    const auto *mg_744 = buffer.data(mg + 744);
    const auto *mg_745 = buffer.data(mg + 745);
    const auto *mg_746 = buffer.data(mg + 746);
    const auto *mg_747 = buffer.data(mg + 747);
    const auto *mg_748 = buffer.data(mg + 748);
    const auto *mg_749 = buffer.data(mg + 749);
    const auto *mg_750 = buffer.data(mg + 750);
    const auto *mg_751 = buffer.data(mg + 751);
    const auto *mg_752 = buffer.data(mg + 752);
    const auto *mg_753 = buffer.data(mg + 753);
    const auto *mg_754 = buffer.data(mg + 754);
    const auto *mg_755 = buffer.data(mg + 755);
    const auto *mg_756 = buffer.data(mg + 756);
    const auto *mg_757 = buffer.data(mg + 757);
    const auto *mg_758 = buffer.data(mg + 758);
    const auto *mg_759 = buffer.data(mg + 759);
    const auto *mg_760 = buffer.data(mg + 760);
    const auto *mg_761 = buffer.data(mg + 761);
    const auto *mg_762 = buffer.data(mg + 762);
    const auto *mg_763 = buffer.data(mg + 763);
    const auto *mg_764 = buffer.data(mg + 764);
    const auto *mg_765 = buffer.data(mg + 765);
    const auto *mg_766 = buffer.data(mg + 766);
    const auto *mg_767 = buffer.data(mg + 767);
    const auto *mg_768 = buffer.data(mg + 768);
    const auto *mg_769 = buffer.data(mg + 769);
    const auto *mg_770 = buffer.data(mg + 770);
    const auto *mg_771 = buffer.data(mg + 771);
    const auto *mg_772 = buffer.data(mg + 772);
    const auto *mg_773 = buffer.data(mg + 773);
    const auto *mg_774 = buffer.data(mg + 774);

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, kg_380, kg_381, kg_382, kg_383, \
                         kg_384, mg_605, mg_606, mg_607, mg_608, \
                         mg_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -3.0 * kg_380[k]
                   + f_0 * mg_605[k];

        t_486[k] = -3.0 * kg_381[k]
                   + f_0 * mg_606[k];

        t_487[k] = -3.0 * kg_382[k]
                   + f_0 * mg_607[k];

        t_488[k] = -3.0 * kg_383[k]
                   + f_0 * mg_608[k];

        t_489[k] = -3.0 * kg_384[k]
                   + f_0 * mg_609[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, kg_385, kg_386, kg_387, kg_388, \
                         kg_389, mg_610, mg_611, mg_612, mg_613, \
                         mg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -3.0 * kg_385[k]
                   + f_0 * mg_610[k];

        t_491[k] = -3.0 * kg_386[k]
                   + f_0 * mg_611[k];

        t_492[k] = -3.0 * kg_387[k]
                   + f_0 * mg_612[k];

        t_493[k] = -3.0 * kg_388[k]
                   + f_0 * mg_613[k];

        t_494[k] = -3.0 * kg_389[k]
                   + f_0 * mg_614[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, kg_390, kg_391, kg_392, kg_393, \
                         kg_394, mg_615, mg_616, mg_617, mg_618, \
                         mg_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -2.0 * kg_390[k]
                   + f_0 * mg_615[k];

        t_496[k] = -2.0 * kg_391[k]
                   + f_0 * mg_616[k];

        t_497[k] = -2.0 * kg_392[k]
                   + f_0 * mg_617[k];

        t_498[k] = -2.0 * kg_393[k]
                   + f_0 * mg_618[k];

        t_499[k] = -2.0 * kg_394[k]
                   + f_0 * mg_619[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, kg_395, kg_396, kg_397, kg_398, \
                         kg_399, mg_620, mg_621, mg_622, mg_623, \
                         mg_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -2.0 * kg_395[k]
                   + f_0 * mg_620[k];

        t_501[k] = -2.0 * kg_396[k]
                   + f_0 * mg_621[k];

        t_502[k] = -2.0 * kg_397[k]
                   + f_0 * mg_622[k];

        t_503[k] = -2.0 * kg_398[k]
                   + f_0 * mg_623[k];

        t_504[k] = -2.0 * kg_399[k]
                   + f_0 * mg_624[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, kg_400, kg_401, kg_402, kg_403, \
                         kg_404, mg_625, mg_626, mg_627, mg_628, \
                         mg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -2.0 * kg_400[k]
                   + f_0 * mg_625[k];

        t_506[k] = -2.0 * kg_401[k]
                   + f_0 * mg_626[k];

        t_507[k] = -2.0 * kg_402[k]
                   + f_0 * mg_627[k];

        t_508[k] = -2.0 * kg_403[k]
                   + f_0 * mg_628[k];

        t_509[k] = -2.0 * kg_404[k]
                   + f_0 * mg_629[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, kg_405, kg_406, kg_407, kg_408, \
                         kg_409, mg_630, mg_631, mg_632, mg_633, \
                         mg_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -kg_405[k]
                   + f_0 * mg_630[k];

        t_511[k] = -kg_406[k]
                   + f_0 * mg_631[k];

        t_512[k] = -kg_407[k]
                   + f_0 * mg_632[k];

        t_513[k] = -kg_408[k]
                   + f_0 * mg_633[k];

        t_514[k] = -kg_409[k]
                   + f_0 * mg_634[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, kg_410, kg_411, kg_412, kg_413, \
                         kg_414, mg_635, mg_636, mg_637, mg_638, \
                         mg_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -kg_410[k]
                   + f_0 * mg_635[k];

        t_516[k] = -kg_411[k]
                   + f_0 * mg_636[k];

        t_517[k] = -kg_412[k]
                   + f_0 * mg_637[k];

        t_518[k] = -kg_413[k]
                   + f_0 * mg_638[k];

        t_519[k] = -kg_414[k]
                   + f_0 * mg_639[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, kg_415, kg_416, kg_417, kg_418, \
                         kg_419, mg_640, mg_641, mg_642, mg_643, \
                         mg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -kg_415[k]
                   + f_0 * mg_640[k];

        t_521[k] = -kg_416[k]
                   + f_0 * mg_641[k];

        t_522[k] = -kg_417[k]
                   + f_0 * mg_642[k];

        t_523[k] = -kg_418[k]
                   + f_0 * mg_643[k];

        t_524[k] = -kg_419[k]
                   + f_0 * mg_644[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, t_530, t_531, t_532, mg_645, \
                         mg_646, mg_647, mg_648, mg_649, mg_650, mg_651, \
                         mg_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_0 * mg_645[k];

        t_526[k] = f_0 * mg_646[k];

        t_527[k] = f_0 * mg_647[k];

        t_528[k] = f_0 * mg_648[k];

        t_529[k] = f_0 * mg_649[k];

        t_530[k] = f_0 * mg_650[k];

        t_531[k] = f_0 * mg_651[k];

        t_532[k] = f_0 * mg_652[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, t_537, t_538, t_539, mg_653, mg_654, \
                         mg_655, mg_656, mg_657, mg_658, mg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_0 * mg_653[k];

        t_534[k] = f_0 * mg_654[k];

        t_535[k] = f_0 * mg_655[k];

        t_536[k] = f_0 * mg_656[k];

        t_537[k] = f_0 * mg_657[k];

        t_538[k] = f_0 * mg_658[k];

        t_539[k] = f_0 * mg_659[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, kg_420, kg_421, kg_422, kg_423, \
                         kg_424, mg_675, mg_676, mg_677, mg_678, \
                         mg_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -8.0 * kg_420[k]
                   + f_0 * mg_675[k];

        t_541[k] = -8.0 * kg_421[k]
                   + f_0 * mg_676[k];

        t_542[k] = -8.0 * kg_422[k]
                   + f_0 * mg_677[k];

        t_543[k] = -8.0 * kg_423[k]
                   + f_0 * mg_678[k];

        t_544[k] = -8.0 * kg_424[k]
                   + f_0 * mg_679[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, kg_425, kg_426, kg_427, kg_428, \
                         kg_429, mg_680, mg_681, mg_682, mg_683, \
                         mg_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -8.0 * kg_425[k]
                   + f_0 * mg_680[k];

        t_546[k] = -8.0 * kg_426[k]
                   + f_0 * mg_681[k];

        t_547[k] = -8.0 * kg_427[k]
                   + f_0 * mg_682[k];

        t_548[k] = -8.0 * kg_428[k]
                   + f_0 * mg_683[k];

        t_549[k] = -8.0 * kg_429[k]
                   + f_0 * mg_684[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, kg_430, kg_431, kg_432, kg_433, \
                         kg_434, mg_685, mg_686, mg_687, mg_688, \
                         mg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -8.0 * kg_430[k]
                   + f_0 * mg_685[k];

        t_551[k] = -8.0 * kg_431[k]
                   + f_0 * mg_686[k];

        t_552[k] = -8.0 * kg_432[k]
                   + f_0 * mg_687[k];

        t_553[k] = -8.0 * kg_433[k]
                   + f_0 * mg_688[k];

        t_554[k] = -8.0 * kg_434[k]
                   + f_0 * mg_689[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, kg_435, kg_436, kg_437, kg_438, \
                         kg_439, mg_690, mg_691, mg_692, mg_693, \
                         mg_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -7.0 * kg_435[k]
                   + f_0 * mg_690[k];

        t_556[k] = -7.0 * kg_436[k]
                   + f_0 * mg_691[k];

        t_557[k] = -7.0 * kg_437[k]
                   + f_0 * mg_692[k];

        t_558[k] = -7.0 * kg_438[k]
                   + f_0 * mg_693[k];

        t_559[k] = -7.0 * kg_439[k]
                   + f_0 * mg_694[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, kg_440, kg_441, kg_442, kg_443, \
                         kg_444, mg_695, mg_696, mg_697, mg_698, \
                         mg_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = -7.0 * kg_440[k]
                   + f_0 * mg_695[k];

        t_561[k] = -7.0 * kg_441[k]
                   + f_0 * mg_696[k];

        t_562[k] = -7.0 * kg_442[k]
                   + f_0 * mg_697[k];

        t_563[k] = -7.0 * kg_443[k]
                   + f_0 * mg_698[k];

        t_564[k] = -7.0 * kg_444[k]
                   + f_0 * mg_699[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, kg_445, kg_446, kg_447, kg_448, \
                         kg_449, mg_700, mg_701, mg_702, mg_703, \
                         mg_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = -7.0 * kg_445[k]
                   + f_0 * mg_700[k];

        t_566[k] = -7.0 * kg_446[k]
                   + f_0 * mg_701[k];

        t_567[k] = -7.0 * kg_447[k]
                   + f_0 * mg_702[k];

        t_568[k] = -7.0 * kg_448[k]
                   + f_0 * mg_703[k];

        t_569[k] = -7.0 * kg_449[k]
                   + f_0 * mg_704[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, kg_450, kg_451, kg_452, kg_453, \
                         kg_454, mg_705, mg_706, mg_707, mg_708, \
                         mg_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -6.0 * kg_450[k]
                   + f_0 * mg_705[k];

        t_571[k] = -6.0 * kg_451[k]
                   + f_0 * mg_706[k];

        t_572[k] = -6.0 * kg_452[k]
                   + f_0 * mg_707[k];

        t_573[k] = -6.0 * kg_453[k]
                   + f_0 * mg_708[k];

        t_574[k] = -6.0 * kg_454[k]
                   + f_0 * mg_709[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, kg_455, kg_456, kg_457, kg_458, \
                         kg_459, mg_710, mg_711, mg_712, mg_713, \
                         mg_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = -6.0 * kg_455[k]
                   + f_0 * mg_710[k];

        t_576[k] = -6.0 * kg_456[k]
                   + f_0 * mg_711[k];

        t_577[k] = -6.0 * kg_457[k]
                   + f_0 * mg_712[k];

        t_578[k] = -6.0 * kg_458[k]
                   + f_0 * mg_713[k];

        t_579[k] = -6.0 * kg_459[k]
                   + f_0 * mg_714[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, kg_460, kg_461, kg_462, kg_463, \
                         kg_464, mg_715, mg_716, mg_717, mg_718, \
                         mg_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -6.0 * kg_460[k]
                   + f_0 * mg_715[k];

        t_581[k] = -6.0 * kg_461[k]
                   + f_0 * mg_716[k];

        t_582[k] = -6.0 * kg_462[k]
                   + f_0 * mg_717[k];

        t_583[k] = -6.0 * kg_463[k]
                   + f_0 * mg_718[k];

        t_584[k] = -6.0 * kg_464[k]
                   + f_0 * mg_719[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, kg_465, kg_466, kg_467, kg_468, \
                         kg_469, mg_720, mg_721, mg_722, mg_723, \
                         mg_724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -5.0 * kg_465[k]
                   + f_0 * mg_720[k];

        t_586[k] = -5.0 * kg_466[k]
                   + f_0 * mg_721[k];

        t_587[k] = -5.0 * kg_467[k]
                   + f_0 * mg_722[k];

        t_588[k] = -5.0 * kg_468[k]
                   + f_0 * mg_723[k];

        t_589[k] = -5.0 * kg_469[k]
                   + f_0 * mg_724[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, kg_470, kg_471, kg_472, kg_473, \
                         kg_474, mg_725, mg_726, mg_727, mg_728, \
                         mg_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -5.0 * kg_470[k]
                   + f_0 * mg_725[k];

        t_591[k] = -5.0 * kg_471[k]
                   + f_0 * mg_726[k];

        t_592[k] = -5.0 * kg_472[k]
                   + f_0 * mg_727[k];

        t_593[k] = -5.0 * kg_473[k]
                   + f_0 * mg_728[k];

        t_594[k] = -5.0 * kg_474[k]
                   + f_0 * mg_729[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, kg_475, kg_476, kg_477, kg_478, \
                         kg_479, mg_730, mg_731, mg_732, mg_733, \
                         mg_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -5.0 * kg_475[k]
                   + f_0 * mg_730[k];

        t_596[k] = -5.0 * kg_476[k]
                   + f_0 * mg_731[k];

        t_597[k] = -5.0 * kg_477[k]
                   + f_0 * mg_732[k];

        t_598[k] = -5.0 * kg_478[k]
                   + f_0 * mg_733[k];

        t_599[k] = -5.0 * kg_479[k]
                   + f_0 * mg_734[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, kg_480, kg_481, kg_482, kg_483, \
                         kg_484, mg_735, mg_736, mg_737, mg_738, \
                         mg_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -4.0 * kg_480[k]
                   + f_0 * mg_735[k];

        t_601[k] = -4.0 * kg_481[k]
                   + f_0 * mg_736[k];

        t_602[k] = -4.0 * kg_482[k]
                   + f_0 * mg_737[k];

        t_603[k] = -4.0 * kg_483[k]
                   + f_0 * mg_738[k];

        t_604[k] = -4.0 * kg_484[k]
                   + f_0 * mg_739[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, kg_485, kg_486, kg_487, kg_488, \
                         kg_489, mg_740, mg_741, mg_742, mg_743, \
                         mg_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -4.0 * kg_485[k]
                   + f_0 * mg_740[k];

        t_606[k] = -4.0 * kg_486[k]
                   + f_0 * mg_741[k];

        t_607[k] = -4.0 * kg_487[k]
                   + f_0 * mg_742[k];

        t_608[k] = -4.0 * kg_488[k]
                   + f_0 * mg_743[k];

        t_609[k] = -4.0 * kg_489[k]
                   + f_0 * mg_744[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, kg_490, kg_491, kg_492, kg_493, \
                         kg_494, mg_745, mg_746, mg_747, mg_748, \
                         mg_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -4.0 * kg_490[k]
                   + f_0 * mg_745[k];

        t_611[k] = -4.0 * kg_491[k]
                   + f_0 * mg_746[k];

        t_612[k] = -4.0 * kg_492[k]
                   + f_0 * mg_747[k];

        t_613[k] = -4.0 * kg_493[k]
                   + f_0 * mg_748[k];

        t_614[k] = -4.0 * kg_494[k]
                   + f_0 * mg_749[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, kg_495, kg_496, kg_497, kg_498, \
                         kg_499, mg_750, mg_751, mg_752, mg_753, \
                         mg_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -3.0 * kg_495[k]
                   + f_0 * mg_750[k];

        t_616[k] = -3.0 * kg_496[k]
                   + f_0 * mg_751[k];

        t_617[k] = -3.0 * kg_497[k]
                   + f_0 * mg_752[k];

        t_618[k] = -3.0 * kg_498[k]
                   + f_0 * mg_753[k];

        t_619[k] = -3.0 * kg_499[k]
                   + f_0 * mg_754[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, kg_500, kg_501, kg_502, kg_503, \
                         kg_504, mg_755, mg_756, mg_757, mg_758, \
                         mg_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -3.0 * kg_500[k]
                   + f_0 * mg_755[k];

        t_621[k] = -3.0 * kg_501[k]
                   + f_0 * mg_756[k];

        t_622[k] = -3.0 * kg_502[k]
                   + f_0 * mg_757[k];

        t_623[k] = -3.0 * kg_503[k]
                   + f_0 * mg_758[k];

        t_624[k] = -3.0 * kg_504[k]
                   + f_0 * mg_759[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, kg_505, kg_506, kg_507, kg_508, \
                         kg_509, mg_760, mg_761, mg_762, mg_763, \
                         mg_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -3.0 * kg_505[k]
                   + f_0 * mg_760[k];

        t_626[k] = -3.0 * kg_506[k]
                   + f_0 * mg_761[k];

        t_627[k] = -3.0 * kg_507[k]
                   + f_0 * mg_762[k];

        t_628[k] = -3.0 * kg_508[k]
                   + f_0 * mg_763[k];

        t_629[k] = -3.0 * kg_509[k]
                   + f_0 * mg_764[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, kg_510, kg_511, kg_512, kg_513, \
                         kg_514, mg_765, mg_766, mg_767, mg_768, \
                         mg_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -2.0 * kg_510[k]
                   + f_0 * mg_765[k];

        t_631[k] = -2.0 * kg_511[k]
                   + f_0 * mg_766[k];

        t_632[k] = -2.0 * kg_512[k]
                   + f_0 * mg_767[k];

        t_633[k] = -2.0 * kg_513[k]
                   + f_0 * mg_768[k];

        t_634[k] = -2.0 * kg_514[k]
                   + f_0 * mg_769[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, kg_515, kg_516, kg_517, kg_518, \
                         kg_519, mg_770, mg_771, mg_772, mg_773, \
                         mg_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -2.0 * kg_515[k]
                   + f_0 * mg_770[k];

        t_636[k] = -2.0 * kg_516[k]
                   + f_0 * mg_771[k];

        t_637[k] = -2.0 * kg_517[k]
                   + f_0 * mg_772[k];

        t_638[k] = -2.0 * kg_518[k]
                   + f_0 * mg_773[k];

        t_639[k] = -2.0 * kg_519[k]
                   + f_0 * mg_774[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *mg_775 = buffer.data(mg + 775);
    const auto *mg_776 = buffer.data(mg + 776);
    const auto *mg_777 = buffer.data(mg + 777);
    const auto *mg_778 = buffer.data(mg + 778);
    const auto *mg_779 = buffer.data(mg + 779);
    const auto *mg_780 = buffer.data(mg + 780);
    const auto *mg_781 = buffer.data(mg + 781);
    const auto *mg_782 = buffer.data(mg + 782);
    const auto *mg_783 = buffer.data(mg + 783);
    const auto *mg_784 = buffer.data(mg + 784);
    const auto *mg_785 = buffer.data(mg + 785);
    const auto *mg_786 = buffer.data(mg + 786);
    const auto *mg_787 = buffer.data(mg + 787);
    const auto *mg_788 = buffer.data(mg + 788);
    const auto *mg_789 = buffer.data(mg + 789);
    const auto *mg_790 = buffer.data(mg + 790);
    const auto *mg_791 = buffer.data(mg + 791);
    const auto *mg_792 = buffer.data(mg + 792);
    const auto *mg_793 = buffer.data(mg + 793);
    const auto *mg_794 = buffer.data(mg + 794);
    const auto *mg_795 = buffer.data(mg + 795);
    const auto *mg_796 = buffer.data(mg + 796);
    const auto *mg_797 = buffer.data(mg + 797);
    const auto *mg_798 = buffer.data(mg + 798);
    const auto *mg_799 = buffer.data(mg + 799);
    const auto *mg_800 = buffer.data(mg + 800);
    const auto *mg_801 = buffer.data(mg + 801);
    const auto *mg_802 = buffer.data(mg + 802);
    const auto *mg_803 = buffer.data(mg + 803);
    const auto *mg_804 = buffer.data(mg + 804);
    const auto *mg_805 = buffer.data(mg + 805);
    const auto *mg_806 = buffer.data(mg + 806);
    const auto *mg_807 = buffer.data(mg + 807);
    const auto *mg_808 = buffer.data(mg + 808);
    const auto *mg_809 = buffer.data(mg + 809);

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, kg_520, kg_521, kg_522, kg_523, \
                         kg_524, mg_775, mg_776, mg_777, mg_778, \
                         mg_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -2.0 * kg_520[k]
                   + f_0 * mg_775[k];

        t_641[k] = -2.0 * kg_521[k]
                   + f_0 * mg_776[k];

        t_642[k] = -2.0 * kg_522[k]
                   + f_0 * mg_777[k];

        t_643[k] = -2.0 * kg_523[k]
                   + f_0 * mg_778[k];

        t_644[k] = -2.0 * kg_524[k]
                   + f_0 * mg_779[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, kg_525, kg_526, kg_527, kg_528, \
                         kg_529, mg_780, mg_781, mg_782, mg_783, \
                         mg_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -kg_525[k]
                   + f_0 * mg_780[k];

        t_646[k] = -kg_526[k]
                   + f_0 * mg_781[k];

        t_647[k] = -kg_527[k]
                   + f_0 * mg_782[k];

        t_648[k] = -kg_528[k]
                   + f_0 * mg_783[k];

        t_649[k] = -kg_529[k]
                   + f_0 * mg_784[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, kg_530, kg_531, kg_532, kg_533, \
                         kg_534, mg_785, mg_786, mg_787, mg_788, \
                         mg_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -kg_530[k]
                   + f_0 * mg_785[k];

        t_651[k] = -kg_531[k]
                   + f_0 * mg_786[k];

        t_652[k] = -kg_532[k]
                   + f_0 * mg_787[k];

        t_653[k] = -kg_533[k]
                   + f_0 * mg_788[k];

        t_654[k] = -kg_534[k]
                   + f_0 * mg_789[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, kg_535, kg_536, kg_537, kg_538, \
                         kg_539, mg_790, mg_791, mg_792, mg_793, \
                         mg_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -kg_535[k]
                   + f_0 * mg_790[k];

        t_656[k] = -kg_536[k]
                   + f_0 * mg_791[k];

        t_657[k] = -kg_537[k]
                   + f_0 * mg_792[k];

        t_658[k] = -kg_538[k]
                   + f_0 * mg_793[k];

        t_659[k] = -kg_539[k]
                   + f_0 * mg_794[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, t_665, t_666, t_667, mg_795, \
                         mg_796, mg_797, mg_798, mg_799, mg_800, mg_801, \
                         mg_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_0 * mg_795[k];

        t_661[k] = f_0 * mg_796[k];

        t_662[k] = f_0 * mg_797[k];

        t_663[k] = f_0 * mg_798[k];

        t_664[k] = f_0 * mg_799[k];

        t_665[k] = f_0 * mg_800[k];

        t_666[k] = f_0 * mg_801[k];

        t_667[k] = f_0 * mg_802[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, t_672, t_673, t_674, mg_803, mg_804, \
                         mg_805, mg_806, mg_807, mg_808, mg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_0 * mg_803[k];

        t_669[k] = f_0 * mg_804[k];

        t_670[k] = f_0 * mg_805[k];

        t_671[k] = f_0 * mg_806[k];

        t_672[k] = f_0 * mg_807[k];

        t_673[k] = f_0 * mg_808[k];

        t_674[k] = f_0 * mg_809[k];
    }
}

auto
compute_prim_geom_10_lg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t kg, const size_t mg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_lg_electron_repulsion_1_piece0(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_1_piece1(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_1_piece2(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_1_piece3(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_1_piece4(buffer, target, kg, mg, ncols, alpha);
}

static auto
compute_prim_geom_10_lg_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
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

    const auto *mg_30 = buffer.data(mg + 30);
    const auto *mg_31 = buffer.data(mg + 31);
    const auto *mg_32 = buffer.data(mg + 32);
    const auto *mg_33 = buffer.data(mg + 33);
    const auto *mg_34 = buffer.data(mg + 34);
    const auto *mg_35 = buffer.data(mg + 35);
    const auto *mg_36 = buffer.data(mg + 36);
    const auto *mg_37 = buffer.data(mg + 37);
    const auto *mg_38 = buffer.data(mg + 38);
    const auto *mg_39 = buffer.data(mg + 39);
    const auto *mg_40 = buffer.data(mg + 40);
    const auto *mg_41 = buffer.data(mg + 41);
    const auto *mg_42 = buffer.data(mg + 42);
    const auto *mg_43 = buffer.data(mg + 43);
    const auto *mg_44 = buffer.data(mg + 44);
    const auto *mg_60 = buffer.data(mg + 60);
    const auto *mg_61 = buffer.data(mg + 61);
    const auto *mg_62 = buffer.data(mg + 62);
    const auto *mg_63 = buffer.data(mg + 63);
    const auto *mg_64 = buffer.data(mg + 64);
    const auto *mg_65 = buffer.data(mg + 65);
    const auto *mg_66 = buffer.data(mg + 66);
    const auto *mg_67 = buffer.data(mg + 67);
    const auto *mg_68 = buffer.data(mg + 68);
    const auto *mg_69 = buffer.data(mg + 69);
    const auto *mg_70 = buffer.data(mg + 70);
    const auto *mg_71 = buffer.data(mg + 71);
    const auto *mg_72 = buffer.data(mg + 72);
    const auto *mg_73 = buffer.data(mg + 73);
    const auto *mg_74 = buffer.data(mg + 74);
    const auto *mg_75 = buffer.data(mg + 75);
    const auto *mg_76 = buffer.data(mg + 76);
    const auto *mg_77 = buffer.data(mg + 77);
    const auto *mg_78 = buffer.data(mg + 78);
    const auto *mg_79 = buffer.data(mg + 79);
    const auto *mg_80 = buffer.data(mg + 80);
    const auto *mg_81 = buffer.data(mg + 81);
    const auto *mg_82 = buffer.data(mg + 82);
    const auto *mg_83 = buffer.data(mg + 83);
    const auto *mg_84 = buffer.data(mg + 84);
    const auto *mg_85 = buffer.data(mg + 85);
    const auto *mg_86 = buffer.data(mg + 86);
    const auto *mg_87 = buffer.data(mg + 87);
    const auto *mg_88 = buffer.data(mg + 88);
    const auto *mg_89 = buffer.data(mg + 89);
    const auto *mg_105 = buffer.data(mg + 105);
    const auto *mg_106 = buffer.data(mg + 106);
    const auto *mg_107 = buffer.data(mg + 107);
    const auto *mg_108 = buffer.data(mg + 108);
    const auto *mg_109 = buffer.data(mg + 109);
    const auto *mg_110 = buffer.data(mg + 110);
    const auto *mg_111 = buffer.data(mg + 111);
    const auto *mg_112 = buffer.data(mg + 112);
    const auto *mg_113 = buffer.data(mg + 113);
    const auto *mg_114 = buffer.data(mg + 114);
    const auto *mg_115 = buffer.data(mg + 115);
    const auto *mg_116 = buffer.data(mg + 116);
    const auto *mg_117 = buffer.data(mg + 117);
    const auto *mg_118 = buffer.data(mg + 118);
    const auto *mg_119 = buffer.data(mg + 119);
    const auto *mg_120 = buffer.data(mg + 120);
    const auto *mg_121 = buffer.data(mg + 121);
    const auto *mg_122 = buffer.data(mg + 122);
    const auto *mg_123 = buffer.data(mg + 123);
    const auto *mg_124 = buffer.data(mg + 124);
    const auto *mg_125 = buffer.data(mg + 125);
    const auto *mg_126 = buffer.data(mg + 126);
    const auto *mg_127 = buffer.data(mg + 127);
    const auto *mg_128 = buffer.data(mg + 128);
    const auto *mg_129 = buffer.data(mg + 129);
    const auto *mg_130 = buffer.data(mg + 130);
    const auto *mg_131 = buffer.data(mg + 131);
    const auto *mg_132 = buffer.data(mg + 132);
    const auto *mg_133 = buffer.data(mg + 133);
    const auto *mg_134 = buffer.data(mg + 134);
    const auto *mg_135 = buffer.data(mg + 135);
    const auto *mg_136 = buffer.data(mg + 136);
    const auto *mg_137 = buffer.data(mg + 137);
    const auto *mg_138 = buffer.data(mg + 138);
    const auto *mg_139 = buffer.data(mg + 139);
    const auto *mg_140 = buffer.data(mg + 140);
    const auto *mg_141 = buffer.data(mg + 141);
    const auto *mg_142 = buffer.data(mg + 142);
    const auto *mg_143 = buffer.data(mg + 143);
    const auto *mg_144 = buffer.data(mg + 144);
    const auto *mg_145 = buffer.data(mg + 145);
    const auto *mg_146 = buffer.data(mg + 146);
    const auto *mg_147 = buffer.data(mg + 147);
    const auto *mg_148 = buffer.data(mg + 148);
    const auto *mg_149 = buffer.data(mg + 149);
    const auto *mg_165 = buffer.data(mg + 165);
    const auto *mg_166 = buffer.data(mg + 166);
    const auto *mg_167 = buffer.data(mg + 167);
    const auto *mg_168 = buffer.data(mg + 168);
    const auto *mg_169 = buffer.data(mg + 169);
    const auto *mg_170 = buffer.data(mg + 170);
    const auto *mg_171 = buffer.data(mg + 171);
    const auto *mg_172 = buffer.data(mg + 172);
    const auto *mg_173 = buffer.data(mg + 173);
    const auto *mg_174 = buffer.data(mg + 174);
    const auto *mg_175 = buffer.data(mg + 175);
    const auto *mg_176 = buffer.data(mg + 176);
    const auto *mg_177 = buffer.data(mg + 177);
    const auto *mg_178 = buffer.data(mg + 178);
    const auto *mg_179 = buffer.data(mg + 179);
    const auto *mg_180 = buffer.data(mg + 180);
    const auto *mg_181 = buffer.data(mg + 181);
    const auto *mg_182 = buffer.data(mg + 182);
    const auto *mg_183 = buffer.data(mg + 183);
    const auto *mg_184 = buffer.data(mg + 184);
    const auto *mg_185 = buffer.data(mg + 185);
    const auto *mg_186 = buffer.data(mg + 186);
    const auto *mg_187 = buffer.data(mg + 187);
    const auto *mg_188 = buffer.data(mg + 188);
    const auto *mg_189 = buffer.data(mg + 189);
    const auto *mg_190 = buffer.data(mg + 190);
    const auto *mg_191 = buffer.data(mg + 191);
    const auto *mg_192 = buffer.data(mg + 192);
    const auto *mg_193 = buffer.data(mg + 193);
    const auto *mg_194 = buffer.data(mg + 194);
    const auto *mg_195 = buffer.data(mg + 195);
    const auto *mg_196 = buffer.data(mg + 196);
    const auto *mg_197 = buffer.data(mg + 197);
    const auto *mg_198 = buffer.data(mg + 198);
    const auto *mg_199 = buffer.data(mg + 199);
    const auto *mg_200 = buffer.data(mg + 200);
    const auto *mg_201 = buffer.data(mg + 201);
    const auto *mg_202 = buffer.data(mg + 202);
    const auto *mg_203 = buffer.data(mg + 203);
    const auto *mg_204 = buffer.data(mg + 204);
    const auto *mg_205 = buffer.data(mg + 205);
    const auto *mg_206 = buffer.data(mg + 206);
    const auto *mg_207 = buffer.data(mg + 207);
    const auto *mg_208 = buffer.data(mg + 208);
    const auto *mg_209 = buffer.data(mg + 209);
    const auto *mg_210 = buffer.data(mg + 210);
    const auto *mg_211 = buffer.data(mg + 211);
    const auto *mg_212 = buffer.data(mg + 212);
    const auto *mg_213 = buffer.data(mg + 213);
    const auto *mg_214 = buffer.data(mg + 214);
    const auto *mg_215 = buffer.data(mg + 215);
    const auto *mg_216 = buffer.data(mg + 216);
    const auto *mg_217 = buffer.data(mg + 217);
    const auto *mg_218 = buffer.data(mg + 218);
    const auto *mg_219 = buffer.data(mg + 219);
    const auto *mg_220 = buffer.data(mg + 220);
    const auto *mg_221 = buffer.data(mg + 221);
    const auto *mg_222 = buffer.data(mg + 222);
    const auto *mg_223 = buffer.data(mg + 223);
    const auto *mg_224 = buffer.data(mg + 224);
    const auto *mg_240 = buffer.data(mg + 240);
    const auto *mg_241 = buffer.data(mg + 241);
    const auto *mg_242 = buffer.data(mg + 242);
    const auto *mg_243 = buffer.data(mg + 243);
    const auto *mg_244 = buffer.data(mg + 244);
    const auto *mg_245 = buffer.data(mg + 245);
    const auto *mg_246 = buffer.data(mg + 246);
    const auto *mg_247 = buffer.data(mg + 247);
    const auto *mg_248 = buffer.data(mg + 248);
    const auto *mg_249 = buffer.data(mg + 249);
    const auto *mg_250 = buffer.data(mg + 250);
    const auto *mg_251 = buffer.data(mg + 251);
    const auto *mg_252 = buffer.data(mg + 252);
    const auto *mg_253 = buffer.data(mg + 253);
    const auto *mg_254 = buffer.data(mg + 254);
    const auto *mg_255 = buffer.data(mg + 255);
    const auto *mg_256 = buffer.data(mg + 256);
    const auto *mg_257 = buffer.data(mg + 257);
    const auto *mg_258 = buffer.data(mg + 258);
    const auto *mg_259 = buffer.data(mg + 259);
    const auto *mg_260 = buffer.data(mg + 260);
    const auto *mg_261 = buffer.data(mg + 261);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, mg_30, mg_31, mg_32, mg_33, \
                         mg_34, mg_35, mg_36, mg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * mg_30[k];

        t_1[k] = f_0 * mg_31[k];

        t_2[k] = f_0 * mg_32[k];

        t_3[k] = f_0 * mg_33[k];

        t_4[k] = f_0 * mg_34[k];

        t_5[k] = f_0 * mg_35[k];

        t_6[k] = f_0 * mg_36[k];

        t_7[k] = f_0 * mg_37[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, mg_38, mg_39, mg_40, \
                         mg_41, mg_42, mg_43, mg_44, mg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * mg_38[k];

        t_9[k] = f_0 * mg_39[k];

        t_10[k] = f_0 * mg_40[k];

        t_11[k] = f_0 * mg_41[k];

        t_12[k] = f_0 * mg_42[k];

        t_13[k] = f_0 * mg_43[k];

        t_14[k] = f_0 * mg_44[k];

        t_15[k] = f_0 * mg_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, mg_61, mg_62, mg_63, \
                         mg_64, mg_65, mg_66, mg_67, mg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * mg_61[k];

        t_17[k] = f_0 * mg_62[k];

        t_18[k] = f_0 * mg_63[k];

        t_19[k] = f_0 * mg_64[k];

        t_20[k] = f_0 * mg_65[k];

        t_21[k] = f_0 * mg_66[k];

        t_22[k] = f_0 * mg_67[k];

        t_23[k] = f_0 * mg_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, kg_0, mg_69, mg_70, mg_71, \
                         mg_72, mg_73, mg_74, mg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * mg_69[k];

        t_25[k] = f_0 * mg_70[k];

        t_26[k] = f_0 * mg_71[k];

        t_27[k] = f_0 * mg_72[k];

        t_28[k] = f_0 * mg_73[k];

        t_29[k] = f_0 * mg_74[k];

        t_30[k] = -kg_0[k]
                  + f_0 * mg_75[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, kg_1, kg_2, kg_3, kg_4, kg_5, mg_76, \
                         mg_77, mg_78, mg_79, mg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -kg_1[k]
                  + f_0 * mg_76[k];

        t_32[k] = -kg_2[k]
                  + f_0 * mg_77[k];

        t_33[k] = -kg_3[k]
                  + f_0 * mg_78[k];

        t_34[k] = -kg_4[k]
                  + f_0 * mg_79[k];

        t_35[k] = -kg_5[k]
                  + f_0 * mg_80[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, kg_6, kg_7, kg_8, kg_9, kg_10, mg_81, \
                         mg_82, mg_83, mg_84, mg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -kg_6[k]
                  + f_0 * mg_81[k];

        t_37[k] = -kg_7[k]
                  + f_0 * mg_82[k];

        t_38[k] = -kg_8[k]
                  + f_0 * mg_83[k];

        t_39[k] = -kg_9[k]
                  + f_0 * mg_84[k];

        t_40[k] = -kg_10[k]
                  + f_0 * mg_85[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, t_46, kg_11, kg_12, kg_13, kg_14, \
                         mg_86, mg_87, mg_88, mg_89, mg_105, mg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -kg_11[k]
                  + f_0 * mg_86[k];

        t_42[k] = -kg_12[k]
                  + f_0 * mg_87[k];

        t_43[k] = -kg_13[k]
                  + f_0 * mg_88[k];

        t_44[k] = -kg_14[k]
                  + f_0 * mg_89[k];

        t_45[k] = f_0 * mg_105[k];

        t_46[k] = f_0 * mg_106[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, t_53, t_54, mg_107, mg_108, \
                         mg_109, mg_110, mg_111, mg_112, mg_113, \
                         mg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * mg_107[k];

        t_48[k] = f_0 * mg_108[k];

        t_49[k] = f_0 * mg_109[k];

        t_50[k] = f_0 * mg_110[k];

        t_51[k] = f_0 * mg_111[k];

        t_52[k] = f_0 * mg_112[k];

        t_53[k] = f_0 * mg_113[k];

        t_54[k] = f_0 * mg_114[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, kg_15, kg_16, mg_115, \
                         mg_116, mg_117, mg_118, mg_119, mg_120, \
                         mg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * mg_115[k];

        t_56[k] = f_0 * mg_116[k];

        t_57[k] = f_0 * mg_117[k];

        t_58[k] = f_0 * mg_118[k];

        t_59[k] = f_0 * mg_119[k];

        t_60[k] = -kg_15[k]
                  + f_0 * mg_120[k];

        t_61[k] = -kg_16[k]
                  + f_0 * mg_121[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, kg_17, kg_18, kg_19, kg_20, kg_21, \
                         mg_122, mg_123, mg_124, mg_125, mg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -kg_17[k]
                  + f_0 * mg_122[k];

        t_63[k] = -kg_18[k]
                  + f_0 * mg_123[k];

        t_64[k] = -kg_19[k]
                  + f_0 * mg_124[k];

        t_65[k] = -kg_20[k]
                  + f_0 * mg_125[k];

        t_66[k] = -kg_21[k]
                  + f_0 * mg_126[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, kg_22, kg_23, kg_24, kg_25, kg_26, \
                         mg_127, mg_128, mg_129, mg_130, mg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -kg_22[k]
                  + f_0 * mg_127[k];

        t_68[k] = -kg_23[k]
                  + f_0 * mg_128[k];

        t_69[k] = -kg_24[k]
                  + f_0 * mg_129[k];

        t_70[k] = -kg_25[k]
                  + f_0 * mg_130[k];

        t_71[k] = -kg_26[k]
                  + f_0 * mg_131[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, kg_27, kg_28, kg_29, kg_30, kg_31, \
                         mg_132, mg_133, mg_134, mg_135, mg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -kg_27[k]
                  + f_0 * mg_132[k];

        t_73[k] = -kg_28[k]
                  + f_0 * mg_133[k];

        t_74[k] = -kg_29[k]
                  + f_0 * mg_134[k];

        t_75[k] = -2.0 * kg_30[k]
                  + f_0 * mg_135[k];

        t_76[k] = -2.0 * kg_31[k]
                  + f_0 * mg_136[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, kg_32, kg_33, kg_34, kg_35, kg_36, \
                         mg_137, mg_138, mg_139, mg_140, mg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -2.0 * kg_32[k]
                  + f_0 * mg_137[k];

        t_78[k] = -2.0 * kg_33[k]
                  + f_0 * mg_138[k];

        t_79[k] = -2.0 * kg_34[k]
                  + f_0 * mg_139[k];

        t_80[k] = -2.0 * kg_35[k]
                  + f_0 * mg_140[k];

        t_81[k] = -2.0 * kg_36[k]
                  + f_0 * mg_141[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, kg_37, kg_38, kg_39, kg_40, kg_41, \
                         mg_142, mg_143, mg_144, mg_145, mg_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -2.0 * kg_37[k]
                  + f_0 * mg_142[k];

        t_83[k] = -2.0 * kg_38[k]
                  + f_0 * mg_143[k];

        t_84[k] = -2.0 * kg_39[k]
                  + f_0 * mg_144[k];

        t_85[k] = -2.0 * kg_40[k]
                  + f_0 * mg_145[k];

        t_86[k] = -2.0 * kg_41[k]
                  + f_0 * mg_146[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, kg_42, kg_43, kg_44, mg_147, \
                         mg_148, mg_149, mg_165, mg_166, mg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -2.0 * kg_42[k]
                  + f_0 * mg_147[k];

        t_88[k] = -2.0 * kg_43[k]
                  + f_0 * mg_148[k];

        t_89[k] = -2.0 * kg_44[k]
                  + f_0 * mg_149[k];

        t_90[k] = f_0 * mg_165[k];

        t_91[k] = f_0 * mg_166[k];

        t_92[k] = f_0 * mg_167[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, mg_168, mg_169, \
                         mg_170, mg_171, mg_172, mg_173, mg_174, \
                         mg_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * mg_168[k];

        t_94[k] = f_0 * mg_169[k];

        t_95[k] = f_0 * mg_170[k];

        t_96[k] = f_0 * mg_171[k];

        t_97[k] = f_0 * mg_172[k];

        t_98[k] = f_0 * mg_173[k];

        t_99[k] = f_0 * mg_174[k];

        t_100[k] = f_0 * mg_175[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, kg_45, kg_46, mg_176, \
                         mg_177, mg_178, mg_179, mg_180, mg_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * mg_176[k];

        t_102[k] = f_0 * mg_177[k];

        t_103[k] = f_0 * mg_178[k];

        t_104[k] = f_0 * mg_179[k];

        t_105[k] = -kg_45[k]
                   + f_0 * mg_180[k];

        t_106[k] = -kg_46[k]
                   + f_0 * mg_181[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, kg_47, kg_48, kg_49, kg_50, kg_51, \
                         mg_182, mg_183, mg_184, mg_185, mg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -kg_47[k]
                   + f_0 * mg_182[k];

        t_108[k] = -kg_48[k]
                   + f_0 * mg_183[k];

        t_109[k] = -kg_49[k]
                   + f_0 * mg_184[k];

        t_110[k] = -kg_50[k]
                   + f_0 * mg_185[k];

        t_111[k] = -kg_51[k]
                   + f_0 * mg_186[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, kg_52, kg_53, kg_54, kg_55, kg_56, \
                         mg_187, mg_188, mg_189, mg_190, mg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -kg_52[k]
                   + f_0 * mg_187[k];

        t_113[k] = -kg_53[k]
                   + f_0 * mg_188[k];

        t_114[k] = -kg_54[k]
                   + f_0 * mg_189[k];

        t_115[k] = -kg_55[k]
                   + f_0 * mg_190[k];

        t_116[k] = -kg_56[k]
                   + f_0 * mg_191[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, kg_57, kg_58, kg_59, kg_60, kg_61, \
                         mg_192, mg_193, mg_194, mg_195, mg_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -kg_57[k]
                   + f_0 * mg_192[k];

        t_118[k] = -kg_58[k]
                   + f_0 * mg_193[k];

        t_119[k] = -kg_59[k]
                   + f_0 * mg_194[k];

        t_120[k] = -2.0 * kg_60[k]
                   + f_0 * mg_195[k];

        t_121[k] = -2.0 * kg_61[k]
                   + f_0 * mg_196[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, kg_62, kg_63, kg_64, kg_65, kg_66, \
                         mg_197, mg_198, mg_199, mg_200, mg_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * kg_62[k]
                   + f_0 * mg_197[k];

        t_123[k] = -2.0 * kg_63[k]
                   + f_0 * mg_198[k];

        t_124[k] = -2.0 * kg_64[k]
                   + f_0 * mg_199[k];

        t_125[k] = -2.0 * kg_65[k]
                   + f_0 * mg_200[k];

        t_126[k] = -2.0 * kg_66[k]
                   + f_0 * mg_201[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, kg_67, kg_68, kg_69, kg_70, kg_71, \
                         mg_202, mg_203, mg_204, mg_205, mg_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * kg_67[k]
                   + f_0 * mg_202[k];

        t_128[k] = -2.0 * kg_68[k]
                   + f_0 * mg_203[k];

        t_129[k] = -2.0 * kg_69[k]
                   + f_0 * mg_204[k];

        t_130[k] = -2.0 * kg_70[k]
                   + f_0 * mg_205[k];

        t_131[k] = -2.0 * kg_71[k]
                   + f_0 * mg_206[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, kg_72, kg_73, kg_74, kg_75, kg_76, \
                         mg_207, mg_208, mg_209, mg_210, mg_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * kg_72[k]
                   + f_0 * mg_207[k];

        t_133[k] = -2.0 * kg_73[k]
                   + f_0 * mg_208[k];

        t_134[k] = -2.0 * kg_74[k]
                   + f_0 * mg_209[k];

        t_135[k] = -3.0 * kg_75[k]
                   + f_0 * mg_210[k];

        t_136[k] = -3.0 * kg_76[k]
                   + f_0 * mg_211[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, kg_77, kg_78, kg_79, kg_80, kg_81, \
                         mg_212, mg_213, mg_214, mg_215, mg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -3.0 * kg_77[k]
                   + f_0 * mg_212[k];

        t_138[k] = -3.0 * kg_78[k]
                   + f_0 * mg_213[k];

        t_139[k] = -3.0 * kg_79[k]
                   + f_0 * mg_214[k];

        t_140[k] = -3.0 * kg_80[k]
                   + f_0 * mg_215[k];

        t_141[k] = -3.0 * kg_81[k]
                   + f_0 * mg_216[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, kg_82, kg_83, kg_84, kg_85, kg_86, \
                         mg_217, mg_218, mg_219, mg_220, mg_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -3.0 * kg_82[k]
                   + f_0 * mg_217[k];

        t_143[k] = -3.0 * kg_83[k]
                   + f_0 * mg_218[k];

        t_144[k] = -3.0 * kg_84[k]
                   + f_0 * mg_219[k];

        t_145[k] = -3.0 * kg_85[k]
                   + f_0 * mg_220[k];

        t_146[k] = -3.0 * kg_86[k]
                   + f_0 * mg_221[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, t_152, kg_87, kg_88, kg_89, \
                         mg_222, mg_223, mg_224, mg_240, mg_241, \
                         mg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -3.0 * kg_87[k]
                   + f_0 * mg_222[k];

        t_148[k] = -3.0 * kg_88[k]
                   + f_0 * mg_223[k];

        t_149[k] = -3.0 * kg_89[k]
                   + f_0 * mg_224[k];

        t_150[k] = f_0 * mg_240[k];

        t_151[k] = f_0 * mg_241[k];

        t_152[k] = f_0 * mg_242[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, t_158, t_159, t_160, mg_243, \
                         mg_244, mg_245, mg_246, mg_247, mg_248, mg_249, \
                         mg_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_0 * mg_243[k];

        t_154[k] = f_0 * mg_244[k];

        t_155[k] = f_0 * mg_245[k];

        t_156[k] = f_0 * mg_246[k];

        t_157[k] = f_0 * mg_247[k];

        t_158[k] = f_0 * mg_248[k];

        t_159[k] = f_0 * mg_249[k];

        t_160[k] = f_0 * mg_250[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, kg_90, kg_91, mg_251, \
                         mg_252, mg_253, mg_254, mg_255, mg_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_0 * mg_251[k];

        t_162[k] = f_0 * mg_252[k];

        t_163[k] = f_0 * mg_253[k];

        t_164[k] = f_0 * mg_254[k];

        t_165[k] = -kg_90[k]
                   + f_0 * mg_255[k];

        t_166[k] = -kg_91[k]
                   + f_0 * mg_256[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, kg_92, kg_93, kg_94, kg_95, kg_96, \
                         mg_257, mg_258, mg_259, mg_260, mg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -kg_92[k]
                   + f_0 * mg_257[k];

        t_168[k] = -kg_93[k]
                   + f_0 * mg_258[k];

        t_169[k] = -kg_94[k]
                   + f_0 * mg_259[k];

        t_170[k] = -kg_95[k]
                   + f_0 * mg_260[k];

        t_171[k] = -kg_96[k]
                   + f_0 * mg_261[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
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

    const auto *mg_262 = buffer.data(mg + 262);
    const auto *mg_263 = buffer.data(mg + 263);
    const auto *mg_264 = buffer.data(mg + 264);
    const auto *mg_265 = buffer.data(mg + 265);
    const auto *mg_266 = buffer.data(mg + 266);
    const auto *mg_267 = buffer.data(mg + 267);
    const auto *mg_268 = buffer.data(mg + 268);
    const auto *mg_269 = buffer.data(mg + 269);
    const auto *mg_270 = buffer.data(mg + 270);
    const auto *mg_271 = buffer.data(mg + 271);
    const auto *mg_272 = buffer.data(mg + 272);
    const auto *mg_273 = buffer.data(mg + 273);
    const auto *mg_274 = buffer.data(mg + 274);
    const auto *mg_275 = buffer.data(mg + 275);
    const auto *mg_276 = buffer.data(mg + 276);
    const auto *mg_277 = buffer.data(mg + 277);
    const auto *mg_278 = buffer.data(mg + 278);
    const auto *mg_279 = buffer.data(mg + 279);
    const auto *mg_280 = buffer.data(mg + 280);
    const auto *mg_281 = buffer.data(mg + 281);
    const auto *mg_282 = buffer.data(mg + 282);
    const auto *mg_283 = buffer.data(mg + 283);
    const auto *mg_284 = buffer.data(mg + 284);
    const auto *mg_285 = buffer.data(mg + 285);
    const auto *mg_286 = buffer.data(mg + 286);
    const auto *mg_287 = buffer.data(mg + 287);
    const auto *mg_288 = buffer.data(mg + 288);
    const auto *mg_289 = buffer.data(mg + 289);
    const auto *mg_290 = buffer.data(mg + 290);
    const auto *mg_291 = buffer.data(mg + 291);
    const auto *mg_292 = buffer.data(mg + 292);
    const auto *mg_293 = buffer.data(mg + 293);
    const auto *mg_294 = buffer.data(mg + 294);
    const auto *mg_295 = buffer.data(mg + 295);
    const auto *mg_296 = buffer.data(mg + 296);
    const auto *mg_297 = buffer.data(mg + 297);
    const auto *mg_298 = buffer.data(mg + 298);
    const auto *mg_299 = buffer.data(mg + 299);
    const auto *mg_300 = buffer.data(mg + 300);
    const auto *mg_301 = buffer.data(mg + 301);
    const auto *mg_302 = buffer.data(mg + 302);
    const auto *mg_303 = buffer.data(mg + 303);
    const auto *mg_304 = buffer.data(mg + 304);
    const auto *mg_305 = buffer.data(mg + 305);
    const auto *mg_306 = buffer.data(mg + 306);
    const auto *mg_307 = buffer.data(mg + 307);
    const auto *mg_308 = buffer.data(mg + 308);
    const auto *mg_309 = buffer.data(mg + 309);
    const auto *mg_310 = buffer.data(mg + 310);
    const auto *mg_311 = buffer.data(mg + 311);
    const auto *mg_312 = buffer.data(mg + 312);
    const auto *mg_313 = buffer.data(mg + 313);
    const auto *mg_314 = buffer.data(mg + 314);
    const auto *mg_330 = buffer.data(mg + 330);
    const auto *mg_331 = buffer.data(mg + 331);
    const auto *mg_332 = buffer.data(mg + 332);
    const auto *mg_333 = buffer.data(mg + 333);
    const auto *mg_334 = buffer.data(mg + 334);
    const auto *mg_335 = buffer.data(mg + 335);
    const auto *mg_336 = buffer.data(mg + 336);
    const auto *mg_337 = buffer.data(mg + 337);
    const auto *mg_338 = buffer.data(mg + 338);
    const auto *mg_339 = buffer.data(mg + 339);
    const auto *mg_340 = buffer.data(mg + 340);
    const auto *mg_341 = buffer.data(mg + 341);
    const auto *mg_342 = buffer.data(mg + 342);
    const auto *mg_343 = buffer.data(mg + 343);
    const auto *mg_344 = buffer.data(mg + 344);
    const auto *mg_345 = buffer.data(mg + 345);
    const auto *mg_346 = buffer.data(mg + 346);
    const auto *mg_347 = buffer.data(mg + 347);
    const auto *mg_348 = buffer.data(mg + 348);
    const auto *mg_349 = buffer.data(mg + 349);
    const auto *mg_350 = buffer.data(mg + 350);
    const auto *mg_351 = buffer.data(mg + 351);
    const auto *mg_352 = buffer.data(mg + 352);
    const auto *mg_353 = buffer.data(mg + 353);
    const auto *mg_354 = buffer.data(mg + 354);
    const auto *mg_355 = buffer.data(mg + 355);
    const auto *mg_356 = buffer.data(mg + 356);
    const auto *mg_357 = buffer.data(mg + 357);
    const auto *mg_358 = buffer.data(mg + 358);
    const auto *mg_359 = buffer.data(mg + 359);
    const auto *mg_360 = buffer.data(mg + 360);
    const auto *mg_361 = buffer.data(mg + 361);
    const auto *mg_362 = buffer.data(mg + 362);
    const auto *mg_363 = buffer.data(mg + 363);
    const auto *mg_364 = buffer.data(mg + 364);
    const auto *mg_365 = buffer.data(mg + 365);
    const auto *mg_366 = buffer.data(mg + 366);
    const auto *mg_367 = buffer.data(mg + 367);
    const auto *mg_368 = buffer.data(mg + 368);
    const auto *mg_369 = buffer.data(mg + 369);
    const auto *mg_370 = buffer.data(mg + 370);
    const auto *mg_371 = buffer.data(mg + 371);
    const auto *mg_372 = buffer.data(mg + 372);
    const auto *mg_373 = buffer.data(mg + 373);
    const auto *mg_374 = buffer.data(mg + 374);
    const auto *mg_375 = buffer.data(mg + 375);
    const auto *mg_376 = buffer.data(mg + 376);
    const auto *mg_377 = buffer.data(mg + 377);
    const auto *mg_378 = buffer.data(mg + 378);
    const auto *mg_379 = buffer.data(mg + 379);
    const auto *mg_380 = buffer.data(mg + 380);
    const auto *mg_381 = buffer.data(mg + 381);
    const auto *mg_382 = buffer.data(mg + 382);
    const auto *mg_383 = buffer.data(mg + 383);
    const auto *mg_384 = buffer.data(mg + 384);
    const auto *mg_385 = buffer.data(mg + 385);
    const auto *mg_386 = buffer.data(mg + 386);
    const auto *mg_387 = buffer.data(mg + 387);
    const auto *mg_388 = buffer.data(mg + 388);
    const auto *mg_389 = buffer.data(mg + 389);
    const auto *mg_390 = buffer.data(mg + 390);
    const auto *mg_391 = buffer.data(mg + 391);
    const auto *mg_392 = buffer.data(mg + 392);
    const auto *mg_393 = buffer.data(mg + 393);
    const auto *mg_394 = buffer.data(mg + 394);
    const auto *mg_395 = buffer.data(mg + 395);
    const auto *mg_396 = buffer.data(mg + 396);
    const auto *mg_397 = buffer.data(mg + 397);
    const auto *mg_398 = buffer.data(mg + 398);
    const auto *mg_399 = buffer.data(mg + 399);
    const auto *mg_400 = buffer.data(mg + 400);
    const auto *mg_401 = buffer.data(mg + 401);
    const auto *mg_402 = buffer.data(mg + 402);
    const auto *mg_403 = buffer.data(mg + 403);
    const auto *mg_404 = buffer.data(mg + 404);
    const auto *mg_405 = buffer.data(mg + 405);
    const auto *mg_406 = buffer.data(mg + 406);
    const auto *mg_407 = buffer.data(mg + 407);
    const auto *mg_408 = buffer.data(mg + 408);
    const auto *mg_409 = buffer.data(mg + 409);
    const auto *mg_410 = buffer.data(mg + 410);
    const auto *mg_411 = buffer.data(mg + 411);
    const auto *mg_412 = buffer.data(mg + 412);
    const auto *mg_413 = buffer.data(mg + 413);
    const auto *mg_414 = buffer.data(mg + 414);
    const auto *mg_415 = buffer.data(mg + 415);
    const auto *mg_416 = buffer.data(mg + 416);
    const auto *mg_417 = buffer.data(mg + 417);
    const auto *mg_418 = buffer.data(mg + 418);
    const auto *mg_419 = buffer.data(mg + 419);
    const auto *mg_435 = buffer.data(mg + 435);
    const auto *mg_436 = buffer.data(mg + 436);
    const auto *mg_437 = buffer.data(mg + 437);
    const auto *mg_438 = buffer.data(mg + 438);
    const auto *mg_439 = buffer.data(mg + 439);
    const auto *mg_440 = buffer.data(mg + 440);
    const auto *mg_441 = buffer.data(mg + 441);
    const auto *mg_442 = buffer.data(mg + 442);
    const auto *mg_443 = buffer.data(mg + 443);
    const auto *mg_444 = buffer.data(mg + 444);
    const auto *mg_445 = buffer.data(mg + 445);
    const auto *mg_446 = buffer.data(mg + 446);
    const auto *mg_447 = buffer.data(mg + 447);
    const auto *mg_448 = buffer.data(mg + 448);
    const auto *mg_449 = buffer.data(mg + 449);
    const auto *mg_450 = buffer.data(mg + 450);
    const auto *mg_451 = buffer.data(mg + 451);

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, kg_97, kg_98, kg_99, kg_100, \
                         kg_101, mg_262, mg_263, mg_264, mg_265, \
                         mg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -kg_97[k]
                   + f_0 * mg_262[k];

        t_173[k] = -kg_98[k]
                   + f_0 * mg_263[k];

        t_174[k] = -kg_99[k]
                   + f_0 * mg_264[k];

        t_175[k] = -kg_100[k]
                   + f_0 * mg_265[k];

        t_176[k] = -kg_101[k]
                   + f_0 * mg_266[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, kg_102, kg_103, kg_104, kg_105, \
                         kg_106, mg_267, mg_268, mg_269, mg_270, \
                         mg_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -kg_102[k]
                   + f_0 * mg_267[k];

        t_178[k] = -kg_103[k]
                   + f_0 * mg_268[k];

        t_179[k] = -kg_104[k]
                   + f_0 * mg_269[k];

        t_180[k] = -2.0 * kg_105[k]
                   + f_0 * mg_270[k];

        t_181[k] = -2.0 * kg_106[k]
                   + f_0 * mg_271[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, kg_107, kg_108, kg_109, kg_110, \
                         kg_111, mg_272, mg_273, mg_274, mg_275, \
                         mg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * kg_107[k]
                   + f_0 * mg_272[k];

        t_183[k] = -2.0 * kg_108[k]
                   + f_0 * mg_273[k];

        t_184[k] = -2.0 * kg_109[k]
                   + f_0 * mg_274[k];

        t_185[k] = -2.0 * kg_110[k]
                   + f_0 * mg_275[k];

        t_186[k] = -2.0 * kg_111[k]
                   + f_0 * mg_276[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, kg_112, kg_113, kg_114, kg_115, \
                         kg_116, mg_277, mg_278, mg_279, mg_280, \
                         mg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * kg_112[k]
                   + f_0 * mg_277[k];

        t_188[k] = -2.0 * kg_113[k]
                   + f_0 * mg_278[k];

        t_189[k] = -2.0 * kg_114[k]
                   + f_0 * mg_279[k];

        t_190[k] = -2.0 * kg_115[k]
                   + f_0 * mg_280[k];

        t_191[k] = -2.0 * kg_116[k]
                   + f_0 * mg_281[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, kg_117, kg_118, kg_119, kg_120, \
                         kg_121, mg_282, mg_283, mg_284, mg_285, \
                         mg_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -2.0 * kg_117[k]
                   + f_0 * mg_282[k];

        t_193[k] = -2.0 * kg_118[k]
                   + f_0 * mg_283[k];

        t_194[k] = -2.0 * kg_119[k]
                   + f_0 * mg_284[k];

        t_195[k] = -3.0 * kg_120[k]
                   + f_0 * mg_285[k];

        t_196[k] = -3.0 * kg_121[k]
                   + f_0 * mg_286[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, kg_122, kg_123, kg_124, kg_125, \
                         kg_126, mg_287, mg_288, mg_289, mg_290, \
                         mg_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -3.0 * kg_122[k]
                   + f_0 * mg_287[k];

        t_198[k] = -3.0 * kg_123[k]
                   + f_0 * mg_288[k];

        t_199[k] = -3.0 * kg_124[k]
                   + f_0 * mg_289[k];

        t_200[k] = -3.0 * kg_125[k]
                   + f_0 * mg_290[k];

        t_201[k] = -3.0 * kg_126[k]
                   + f_0 * mg_291[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, kg_127, kg_128, kg_129, kg_130, \
                         kg_131, mg_292, mg_293, mg_294, mg_295, \
                         mg_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -3.0 * kg_127[k]
                   + f_0 * mg_292[k];

        t_203[k] = -3.0 * kg_128[k]
                   + f_0 * mg_293[k];

        t_204[k] = -3.0 * kg_129[k]
                   + f_0 * mg_294[k];

        t_205[k] = -3.0 * kg_130[k]
                   + f_0 * mg_295[k];

        t_206[k] = -3.0 * kg_131[k]
                   + f_0 * mg_296[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, kg_132, kg_133, kg_134, kg_135, \
                         kg_136, mg_297, mg_298, mg_299, mg_300, \
                         mg_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -3.0 * kg_132[k]
                   + f_0 * mg_297[k];

        t_208[k] = -3.0 * kg_133[k]
                   + f_0 * mg_298[k];

        t_209[k] = -3.0 * kg_134[k]
                   + f_0 * mg_299[k];

        t_210[k] = -4.0 * kg_135[k]
                   + f_0 * mg_300[k];

        t_211[k] = -4.0 * kg_136[k]
                   + f_0 * mg_301[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, kg_137, kg_138, kg_139, kg_140, \
                         kg_141, mg_302, mg_303, mg_304, mg_305, \
                         mg_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = -4.0 * kg_137[k]
                   + f_0 * mg_302[k];

        t_213[k] = -4.0 * kg_138[k]
                   + f_0 * mg_303[k];

        t_214[k] = -4.0 * kg_139[k]
                   + f_0 * mg_304[k];

        t_215[k] = -4.0 * kg_140[k]
                   + f_0 * mg_305[k];

        t_216[k] = -4.0 * kg_141[k]
                   + f_0 * mg_306[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, kg_142, kg_143, kg_144, kg_145, \
                         kg_146, mg_307, mg_308, mg_309, mg_310, \
                         mg_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = -4.0 * kg_142[k]
                   + f_0 * mg_307[k];

        t_218[k] = -4.0 * kg_143[k]
                   + f_0 * mg_308[k];

        t_219[k] = -4.0 * kg_144[k]
                   + f_0 * mg_309[k];

        t_220[k] = -4.0 * kg_145[k]
                   + f_0 * mg_310[k];

        t_221[k] = -4.0 * kg_146[k]
                   + f_0 * mg_311[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, kg_147, kg_148, kg_149, \
                         mg_312, mg_313, mg_314, mg_330, mg_331, \
                         mg_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -4.0 * kg_147[k]
                   + f_0 * mg_312[k];

        t_223[k] = -4.0 * kg_148[k]
                   + f_0 * mg_313[k];

        t_224[k] = -4.0 * kg_149[k]
                   + f_0 * mg_314[k];

        t_225[k] = f_0 * mg_330[k];

        t_226[k] = f_0 * mg_331[k];

        t_227[k] = f_0 * mg_332[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, t_233, t_234, t_235, mg_333, \
                         mg_334, mg_335, mg_336, mg_337, mg_338, mg_339, \
                         mg_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_0 * mg_333[k];

        t_229[k] = f_0 * mg_334[k];

        t_230[k] = f_0 * mg_335[k];

        t_231[k] = f_0 * mg_336[k];

        t_232[k] = f_0 * mg_337[k];

        t_233[k] = f_0 * mg_338[k];

        t_234[k] = f_0 * mg_339[k];

        t_235[k] = f_0 * mg_340[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, t_241, kg_150, kg_151, mg_341, \
                         mg_342, mg_343, mg_344, mg_345, mg_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_0 * mg_341[k];

        t_237[k] = f_0 * mg_342[k];

        t_238[k] = f_0 * mg_343[k];

        t_239[k] = f_0 * mg_344[k];

        t_240[k] = -kg_150[k]
                   + f_0 * mg_345[k];

        t_241[k] = -kg_151[k]
                   + f_0 * mg_346[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, kg_152, kg_153, kg_154, kg_155, \
                         kg_156, mg_347, mg_348, mg_349, mg_350, \
                         mg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = -kg_152[k]
                   + f_0 * mg_347[k];

        t_243[k] = -kg_153[k]
                   + f_0 * mg_348[k];

        t_244[k] = -kg_154[k]
                   + f_0 * mg_349[k];

        t_245[k] = -kg_155[k]
                   + f_0 * mg_350[k];

        t_246[k] = -kg_156[k]
                   + f_0 * mg_351[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, kg_157, kg_158, kg_159, kg_160, \
                         kg_161, mg_352, mg_353, mg_354, mg_355, \
                         mg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = -kg_157[k]
                   + f_0 * mg_352[k];

        t_248[k] = -kg_158[k]
                   + f_0 * mg_353[k];

        t_249[k] = -kg_159[k]
                   + f_0 * mg_354[k];

        t_250[k] = -kg_160[k]
                   + f_0 * mg_355[k];

        t_251[k] = -kg_161[k]
                   + f_0 * mg_356[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, kg_162, kg_163, kg_164, kg_165, \
                         kg_166, mg_357, mg_358, mg_359, mg_360, \
                         mg_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -kg_162[k]
                   + f_0 * mg_357[k];

        t_253[k] = -kg_163[k]
                   + f_0 * mg_358[k];

        t_254[k] = -kg_164[k]
                   + f_0 * mg_359[k];

        t_255[k] = -2.0 * kg_165[k]
                   + f_0 * mg_360[k];

        t_256[k] = -2.0 * kg_166[k]
                   + f_0 * mg_361[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, kg_167, kg_168, kg_169, kg_170, \
                         kg_171, mg_362, mg_363, mg_364, mg_365, \
                         mg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = -2.0 * kg_167[k]
                   + f_0 * mg_362[k];

        t_258[k] = -2.0 * kg_168[k]
                   + f_0 * mg_363[k];

        t_259[k] = -2.0 * kg_169[k]
                   + f_0 * mg_364[k];

        t_260[k] = -2.0 * kg_170[k]
                   + f_0 * mg_365[k];

        t_261[k] = -2.0 * kg_171[k]
                   + f_0 * mg_366[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, kg_172, kg_173, kg_174, kg_175, \
                         kg_176, mg_367, mg_368, mg_369, mg_370, \
                         mg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = -2.0 * kg_172[k]
                   + f_0 * mg_367[k];

        t_263[k] = -2.0 * kg_173[k]
                   + f_0 * mg_368[k];

        t_264[k] = -2.0 * kg_174[k]
                   + f_0 * mg_369[k];

        t_265[k] = -2.0 * kg_175[k]
                   + f_0 * mg_370[k];

        t_266[k] = -2.0 * kg_176[k]
                   + f_0 * mg_371[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, kg_177, kg_178, kg_179, kg_180, \
                         kg_181, mg_372, mg_373, mg_374, mg_375, \
                         mg_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = -2.0 * kg_177[k]
                   + f_0 * mg_372[k];

        t_268[k] = -2.0 * kg_178[k]
                   + f_0 * mg_373[k];

        t_269[k] = -2.0 * kg_179[k]
                   + f_0 * mg_374[k];

        t_270[k] = -3.0 * kg_180[k]
                   + f_0 * mg_375[k];

        t_271[k] = -3.0 * kg_181[k]
                   + f_0 * mg_376[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, kg_182, kg_183, kg_184, kg_185, \
                         kg_186, mg_377, mg_378, mg_379, mg_380, \
                         mg_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -3.0 * kg_182[k]
                   + f_0 * mg_377[k];

        t_273[k] = -3.0 * kg_183[k]
                   + f_0 * mg_378[k];

        t_274[k] = -3.0 * kg_184[k]
                   + f_0 * mg_379[k];

        t_275[k] = -3.0 * kg_185[k]
                   + f_0 * mg_380[k];

        t_276[k] = -3.0 * kg_186[k]
                   + f_0 * mg_381[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, t_281, kg_187, kg_188, kg_189, kg_190, \
                         kg_191, mg_382, mg_383, mg_384, mg_385, \
                         mg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = -3.0 * kg_187[k]
                   + f_0 * mg_382[k];

        t_278[k] = -3.0 * kg_188[k]
                   + f_0 * mg_383[k];

        t_279[k] = -3.0 * kg_189[k]
                   + f_0 * mg_384[k];

        t_280[k] = -3.0 * kg_190[k]
                   + f_0 * mg_385[k];

        t_281[k] = -3.0 * kg_191[k]
                   + f_0 * mg_386[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, kg_192, kg_193, kg_194, kg_195, \
                         kg_196, mg_387, mg_388, mg_389, mg_390, \
                         mg_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = -3.0 * kg_192[k]
                   + f_0 * mg_387[k];

        t_283[k] = -3.0 * kg_193[k]
                   + f_0 * mg_388[k];

        t_284[k] = -3.0 * kg_194[k]
                   + f_0 * mg_389[k];

        t_285[k] = -4.0 * kg_195[k]
                   + f_0 * mg_390[k];

        t_286[k] = -4.0 * kg_196[k]
                   + f_0 * mg_391[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, kg_197, kg_198, kg_199, kg_200, \
                         kg_201, mg_392, mg_393, mg_394, mg_395, \
                         mg_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = -4.0 * kg_197[k]
                   + f_0 * mg_392[k];

        t_288[k] = -4.0 * kg_198[k]
                   + f_0 * mg_393[k];

        t_289[k] = -4.0 * kg_199[k]
                   + f_0 * mg_394[k];

        t_290[k] = -4.0 * kg_200[k]
                   + f_0 * mg_395[k];

        t_291[k] = -4.0 * kg_201[k]
                   + f_0 * mg_396[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, kg_202, kg_203, kg_204, kg_205, \
                         kg_206, mg_397, mg_398, mg_399, mg_400, \
                         mg_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = -4.0 * kg_202[k]
                   + f_0 * mg_397[k];

        t_293[k] = -4.0 * kg_203[k]
                   + f_0 * mg_398[k];

        t_294[k] = -4.0 * kg_204[k]
                   + f_0 * mg_399[k];

        t_295[k] = -4.0 * kg_205[k]
                   + f_0 * mg_400[k];

        t_296[k] = -4.0 * kg_206[k]
                   + f_0 * mg_401[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, kg_207, kg_208, kg_209, kg_210, \
                         kg_211, mg_402, mg_403, mg_404, mg_405, \
                         mg_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = -4.0 * kg_207[k]
                   + f_0 * mg_402[k];

        t_298[k] = -4.0 * kg_208[k]
                   + f_0 * mg_403[k];

        t_299[k] = -4.0 * kg_209[k]
                   + f_0 * mg_404[k];

        t_300[k] = -5.0 * kg_210[k]
                   + f_0 * mg_405[k];

        t_301[k] = -5.0 * kg_211[k]
                   + f_0 * mg_406[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, kg_212, kg_213, kg_214, kg_215, \
                         kg_216, mg_407, mg_408, mg_409, mg_410, \
                         mg_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = -5.0 * kg_212[k]
                   + f_0 * mg_407[k];

        t_303[k] = -5.0 * kg_213[k]
                   + f_0 * mg_408[k];

        t_304[k] = -5.0 * kg_214[k]
                   + f_0 * mg_409[k];

        t_305[k] = -5.0 * kg_215[k]
                   + f_0 * mg_410[k];

        t_306[k] = -5.0 * kg_216[k]
                   + f_0 * mg_411[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, kg_217, kg_218, kg_219, kg_220, \
                         kg_221, mg_412, mg_413, mg_414, mg_415, \
                         mg_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = -5.0 * kg_217[k]
                   + f_0 * mg_412[k];

        t_308[k] = -5.0 * kg_218[k]
                   + f_0 * mg_413[k];

        t_309[k] = -5.0 * kg_219[k]
                   + f_0 * mg_414[k];

        t_310[k] = -5.0 * kg_220[k]
                   + f_0 * mg_415[k];

        t_311[k] = -5.0 * kg_221[k]
                   + f_0 * mg_416[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, t_317, kg_222, kg_223, kg_224, \
                         mg_417, mg_418, mg_419, mg_435, mg_436, \
                         mg_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -5.0 * kg_222[k]
                   + f_0 * mg_417[k];

        t_313[k] = -5.0 * kg_223[k]
                   + f_0 * mg_418[k];

        t_314[k] = -5.0 * kg_224[k]
                   + f_0 * mg_419[k];

        t_315[k] = f_0 * mg_435[k];

        t_316[k] = f_0 * mg_436[k];

        t_317[k] = f_0 * mg_437[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, t_323, t_324, t_325, mg_438, \
                         mg_439, mg_440, mg_441, mg_442, mg_443, mg_444, \
                         mg_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_0 * mg_438[k];

        t_319[k] = f_0 * mg_439[k];

        t_320[k] = f_0 * mg_440[k];

        t_321[k] = f_0 * mg_441[k];

        t_322[k] = f_0 * mg_442[k];

        t_323[k] = f_0 * mg_443[k];

        t_324[k] = f_0 * mg_444[k];

        t_325[k] = f_0 * mg_445[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, kg_225, kg_226, mg_446, \
                         mg_447, mg_448, mg_449, mg_450, mg_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_0 * mg_446[k];

        t_327[k] = f_0 * mg_447[k];

        t_328[k] = f_0 * mg_448[k];

        t_329[k] = f_0 * mg_449[k];

        t_330[k] = -kg_225[k]
                   + f_0 * mg_450[k];

        t_331[k] = -kg_226[k]
                   + f_0 * mg_451[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
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

    const auto *mg_452 = buffer.data(mg + 452);
    const auto *mg_453 = buffer.data(mg + 453);
    const auto *mg_454 = buffer.data(mg + 454);
    const auto *mg_455 = buffer.data(mg + 455);
    const auto *mg_456 = buffer.data(mg + 456);
    const auto *mg_457 = buffer.data(mg + 457);
    const auto *mg_458 = buffer.data(mg + 458);
    const auto *mg_459 = buffer.data(mg + 459);
    const auto *mg_460 = buffer.data(mg + 460);
    const auto *mg_461 = buffer.data(mg + 461);
    const auto *mg_462 = buffer.data(mg + 462);
    const auto *mg_463 = buffer.data(mg + 463);
    const auto *mg_464 = buffer.data(mg + 464);
    const auto *mg_465 = buffer.data(mg + 465);
    const auto *mg_466 = buffer.data(mg + 466);
    const auto *mg_467 = buffer.data(mg + 467);
    const auto *mg_468 = buffer.data(mg + 468);
    const auto *mg_469 = buffer.data(mg + 469);
    const auto *mg_470 = buffer.data(mg + 470);
    const auto *mg_471 = buffer.data(mg + 471);
    const auto *mg_472 = buffer.data(mg + 472);
    const auto *mg_473 = buffer.data(mg + 473);
    const auto *mg_474 = buffer.data(mg + 474);
    const auto *mg_475 = buffer.data(mg + 475);
    const auto *mg_476 = buffer.data(mg + 476);
    const auto *mg_477 = buffer.data(mg + 477);
    const auto *mg_478 = buffer.data(mg + 478);
    const auto *mg_479 = buffer.data(mg + 479);
    const auto *mg_480 = buffer.data(mg + 480);
    const auto *mg_481 = buffer.data(mg + 481);
    const auto *mg_482 = buffer.data(mg + 482);
    const auto *mg_483 = buffer.data(mg + 483);
    const auto *mg_484 = buffer.data(mg + 484);
    const auto *mg_485 = buffer.data(mg + 485);
    const auto *mg_486 = buffer.data(mg + 486);
    const auto *mg_487 = buffer.data(mg + 487);
    const auto *mg_488 = buffer.data(mg + 488);
    const auto *mg_489 = buffer.data(mg + 489);
    const auto *mg_490 = buffer.data(mg + 490);
    const auto *mg_491 = buffer.data(mg + 491);
    const auto *mg_492 = buffer.data(mg + 492);
    const auto *mg_493 = buffer.data(mg + 493);
    const auto *mg_494 = buffer.data(mg + 494);
    const auto *mg_495 = buffer.data(mg + 495);
    const auto *mg_496 = buffer.data(mg + 496);
    const auto *mg_497 = buffer.data(mg + 497);
    const auto *mg_498 = buffer.data(mg + 498);
    const auto *mg_499 = buffer.data(mg + 499);
    const auto *mg_500 = buffer.data(mg + 500);
    const auto *mg_501 = buffer.data(mg + 501);
    const auto *mg_502 = buffer.data(mg + 502);
    const auto *mg_503 = buffer.data(mg + 503);
    const auto *mg_504 = buffer.data(mg + 504);
    const auto *mg_505 = buffer.data(mg + 505);
    const auto *mg_506 = buffer.data(mg + 506);
    const auto *mg_507 = buffer.data(mg + 507);
    const auto *mg_508 = buffer.data(mg + 508);
    const auto *mg_509 = buffer.data(mg + 509);
    const auto *mg_510 = buffer.data(mg + 510);
    const auto *mg_511 = buffer.data(mg + 511);
    const auto *mg_512 = buffer.data(mg + 512);
    const auto *mg_513 = buffer.data(mg + 513);
    const auto *mg_514 = buffer.data(mg + 514);
    const auto *mg_515 = buffer.data(mg + 515);
    const auto *mg_516 = buffer.data(mg + 516);
    const auto *mg_517 = buffer.data(mg + 517);
    const auto *mg_518 = buffer.data(mg + 518);
    const auto *mg_519 = buffer.data(mg + 519);
    const auto *mg_520 = buffer.data(mg + 520);
    const auto *mg_521 = buffer.data(mg + 521);
    const auto *mg_522 = buffer.data(mg + 522);
    const auto *mg_523 = buffer.data(mg + 523);
    const auto *mg_524 = buffer.data(mg + 524);
    const auto *mg_525 = buffer.data(mg + 525);
    const auto *mg_526 = buffer.data(mg + 526);
    const auto *mg_527 = buffer.data(mg + 527);
    const auto *mg_528 = buffer.data(mg + 528);
    const auto *mg_529 = buffer.data(mg + 529);
    const auto *mg_530 = buffer.data(mg + 530);
    const auto *mg_531 = buffer.data(mg + 531);
    const auto *mg_532 = buffer.data(mg + 532);
    const auto *mg_533 = buffer.data(mg + 533);
    const auto *mg_534 = buffer.data(mg + 534);
    const auto *mg_535 = buffer.data(mg + 535);
    const auto *mg_536 = buffer.data(mg + 536);
    const auto *mg_537 = buffer.data(mg + 537);
    const auto *mg_538 = buffer.data(mg + 538);
    const auto *mg_539 = buffer.data(mg + 539);
    const auto *mg_555 = buffer.data(mg + 555);
    const auto *mg_556 = buffer.data(mg + 556);
    const auto *mg_557 = buffer.data(mg + 557);
    const auto *mg_558 = buffer.data(mg + 558);
    const auto *mg_559 = buffer.data(mg + 559);
    const auto *mg_560 = buffer.data(mg + 560);
    const auto *mg_561 = buffer.data(mg + 561);
    const auto *mg_562 = buffer.data(mg + 562);
    const auto *mg_563 = buffer.data(mg + 563);
    const auto *mg_564 = buffer.data(mg + 564);
    const auto *mg_565 = buffer.data(mg + 565);
    const auto *mg_566 = buffer.data(mg + 566);
    const auto *mg_567 = buffer.data(mg + 567);
    const auto *mg_568 = buffer.data(mg + 568);
    const auto *mg_569 = buffer.data(mg + 569);
    const auto *mg_570 = buffer.data(mg + 570);
    const auto *mg_571 = buffer.data(mg + 571);
    const auto *mg_572 = buffer.data(mg + 572);
    const auto *mg_573 = buffer.data(mg + 573);
    const auto *mg_574 = buffer.data(mg + 574);
    const auto *mg_575 = buffer.data(mg + 575);
    const auto *mg_576 = buffer.data(mg + 576);
    const auto *mg_577 = buffer.data(mg + 577);
    const auto *mg_578 = buffer.data(mg + 578);
    const auto *mg_579 = buffer.data(mg + 579);
    const auto *mg_580 = buffer.data(mg + 580);
    const auto *mg_581 = buffer.data(mg + 581);
    const auto *mg_582 = buffer.data(mg + 582);
    const auto *mg_583 = buffer.data(mg + 583);
    const auto *mg_584 = buffer.data(mg + 584);
    const auto *mg_585 = buffer.data(mg + 585);
    const auto *mg_586 = buffer.data(mg + 586);
    const auto *mg_587 = buffer.data(mg + 587);
    const auto *mg_588 = buffer.data(mg + 588);
    const auto *mg_589 = buffer.data(mg + 589);
    const auto *mg_590 = buffer.data(mg + 590);
    const auto *mg_591 = buffer.data(mg + 591);
    const auto *mg_592 = buffer.data(mg + 592);
    const auto *mg_593 = buffer.data(mg + 593);
    const auto *mg_594 = buffer.data(mg + 594);
    const auto *mg_595 = buffer.data(mg + 595);
    const auto *mg_596 = buffer.data(mg + 596);
    const auto *mg_597 = buffer.data(mg + 597);
    const auto *mg_598 = buffer.data(mg + 598);
    const auto *mg_599 = buffer.data(mg + 599);
    const auto *mg_600 = buffer.data(mg + 600);
    const auto *mg_601 = buffer.data(mg + 601);
    const auto *mg_602 = buffer.data(mg + 602);
    const auto *mg_603 = buffer.data(mg + 603);
    const auto *mg_604 = buffer.data(mg + 604);
    const auto *mg_605 = buffer.data(mg + 605);
    const auto *mg_606 = buffer.data(mg + 606);
    const auto *mg_607 = buffer.data(mg + 607);
    const auto *mg_608 = buffer.data(mg + 608);
    const auto *mg_609 = buffer.data(mg + 609);
    const auto *mg_610 = buffer.data(mg + 610);
    const auto *mg_611 = buffer.data(mg + 611);
    const auto *mg_612 = buffer.data(mg + 612);
    const auto *mg_613 = buffer.data(mg + 613);
    const auto *mg_614 = buffer.data(mg + 614);
    const auto *mg_615 = buffer.data(mg + 615);
    const auto *mg_616 = buffer.data(mg + 616);
    const auto *mg_617 = buffer.data(mg + 617);
    const auto *mg_618 = buffer.data(mg + 618);
    const auto *mg_619 = buffer.data(mg + 619);
    const auto *mg_620 = buffer.data(mg + 620);
    const auto *mg_621 = buffer.data(mg + 621);

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, kg_227, kg_228, kg_229, kg_230, \
                         kg_231, mg_452, mg_453, mg_454, mg_455, \
                         mg_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -kg_227[k]
                   + f_0 * mg_452[k];

        t_333[k] = -kg_228[k]
                   + f_0 * mg_453[k];

        t_334[k] = -kg_229[k]
                   + f_0 * mg_454[k];

        t_335[k] = -kg_230[k]
                   + f_0 * mg_455[k];

        t_336[k] = -kg_231[k]
                   + f_0 * mg_456[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, kg_232, kg_233, kg_234, kg_235, \
                         kg_236, mg_457, mg_458, mg_459, mg_460, \
                         mg_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -kg_232[k]
                   + f_0 * mg_457[k];

        t_338[k] = -kg_233[k]
                   + f_0 * mg_458[k];

        t_339[k] = -kg_234[k]
                   + f_0 * mg_459[k];

        t_340[k] = -kg_235[k]
                   + f_0 * mg_460[k];

        t_341[k] = -kg_236[k]
                   + f_0 * mg_461[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, kg_237, kg_238, kg_239, kg_240, \
                         kg_241, mg_462, mg_463, mg_464, mg_465, \
                         mg_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -kg_237[k]
                   + f_0 * mg_462[k];

        t_343[k] = -kg_238[k]
                   + f_0 * mg_463[k];

        t_344[k] = -kg_239[k]
                   + f_0 * mg_464[k];

        t_345[k] = -2.0 * kg_240[k]
                   + f_0 * mg_465[k];

        t_346[k] = -2.0 * kg_241[k]
                   + f_0 * mg_466[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, kg_242, kg_243, kg_244, kg_245, \
                         kg_246, mg_467, mg_468, mg_469, mg_470, \
                         mg_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -2.0 * kg_242[k]
                   + f_0 * mg_467[k];

        t_348[k] = -2.0 * kg_243[k]
                   + f_0 * mg_468[k];

        t_349[k] = -2.0 * kg_244[k]
                   + f_0 * mg_469[k];

        t_350[k] = -2.0 * kg_245[k]
                   + f_0 * mg_470[k];

        t_351[k] = -2.0 * kg_246[k]
                   + f_0 * mg_471[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, kg_247, kg_248, kg_249, kg_250, \
                         kg_251, mg_472, mg_473, mg_474, mg_475, \
                         mg_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -2.0 * kg_247[k]
                   + f_0 * mg_472[k];

        t_353[k] = -2.0 * kg_248[k]
                   + f_0 * mg_473[k];

        t_354[k] = -2.0 * kg_249[k]
                   + f_0 * mg_474[k];

        t_355[k] = -2.0 * kg_250[k]
                   + f_0 * mg_475[k];

        t_356[k] = -2.0 * kg_251[k]
                   + f_0 * mg_476[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, kg_252, kg_253, kg_254, kg_255, \
                         kg_256, mg_477, mg_478, mg_479, mg_480, \
                         mg_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -2.0 * kg_252[k]
                   + f_0 * mg_477[k];

        t_358[k] = -2.0 * kg_253[k]
                   + f_0 * mg_478[k];

        t_359[k] = -2.0 * kg_254[k]
                   + f_0 * mg_479[k];

        t_360[k] = -3.0 * kg_255[k]
                   + f_0 * mg_480[k];

        t_361[k] = -3.0 * kg_256[k]
                   + f_0 * mg_481[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, kg_257, kg_258, kg_259, kg_260, \
                         kg_261, mg_482, mg_483, mg_484, mg_485, \
                         mg_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -3.0 * kg_257[k]
                   + f_0 * mg_482[k];

        t_363[k] = -3.0 * kg_258[k]
                   + f_0 * mg_483[k];

        t_364[k] = -3.0 * kg_259[k]
                   + f_0 * mg_484[k];

        t_365[k] = -3.0 * kg_260[k]
                   + f_0 * mg_485[k];

        t_366[k] = -3.0 * kg_261[k]
                   + f_0 * mg_486[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, kg_262, kg_263, kg_264, kg_265, \
                         kg_266, mg_487, mg_488, mg_489, mg_490, \
                         mg_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -3.0 * kg_262[k]
                   + f_0 * mg_487[k];

        t_368[k] = -3.0 * kg_263[k]
                   + f_0 * mg_488[k];

        t_369[k] = -3.0 * kg_264[k]
                   + f_0 * mg_489[k];

        t_370[k] = -3.0 * kg_265[k]
                   + f_0 * mg_490[k];

        t_371[k] = -3.0 * kg_266[k]
                   + f_0 * mg_491[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, kg_267, kg_268, kg_269, kg_270, \
                         kg_271, mg_492, mg_493, mg_494, mg_495, \
                         mg_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -3.0 * kg_267[k]
                   + f_0 * mg_492[k];

        t_373[k] = -3.0 * kg_268[k]
                   + f_0 * mg_493[k];

        t_374[k] = -3.0 * kg_269[k]
                   + f_0 * mg_494[k];

        t_375[k] = -4.0 * kg_270[k]
                   + f_0 * mg_495[k];

        t_376[k] = -4.0 * kg_271[k]
                   + f_0 * mg_496[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, kg_272, kg_273, kg_274, kg_275, \
                         kg_276, mg_497, mg_498, mg_499, mg_500, \
                         mg_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -4.0 * kg_272[k]
                   + f_0 * mg_497[k];

        t_378[k] = -4.0 * kg_273[k]
                   + f_0 * mg_498[k];

        t_379[k] = -4.0 * kg_274[k]
                   + f_0 * mg_499[k];

        t_380[k] = -4.0 * kg_275[k]
                   + f_0 * mg_500[k];

        t_381[k] = -4.0 * kg_276[k]
                   + f_0 * mg_501[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, kg_277, kg_278, kg_279, kg_280, \
                         kg_281, mg_502, mg_503, mg_504, mg_505, \
                         mg_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -4.0 * kg_277[k]
                   + f_0 * mg_502[k];

        t_383[k] = -4.0 * kg_278[k]
                   + f_0 * mg_503[k];

        t_384[k] = -4.0 * kg_279[k]
                   + f_0 * mg_504[k];

        t_385[k] = -4.0 * kg_280[k]
                   + f_0 * mg_505[k];

        t_386[k] = -4.0 * kg_281[k]
                   + f_0 * mg_506[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, kg_282, kg_283, kg_284, kg_285, \
                         kg_286, mg_507, mg_508, mg_509, mg_510, \
                         mg_511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -4.0 * kg_282[k]
                   + f_0 * mg_507[k];

        t_388[k] = -4.0 * kg_283[k]
                   + f_0 * mg_508[k];

        t_389[k] = -4.0 * kg_284[k]
                   + f_0 * mg_509[k];

        t_390[k] = -5.0 * kg_285[k]
                   + f_0 * mg_510[k];

        t_391[k] = -5.0 * kg_286[k]
                   + f_0 * mg_511[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, kg_287, kg_288, kg_289, kg_290, \
                         kg_291, mg_512, mg_513, mg_514, mg_515, \
                         mg_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -5.0 * kg_287[k]
                   + f_0 * mg_512[k];

        t_393[k] = -5.0 * kg_288[k]
                   + f_0 * mg_513[k];

        t_394[k] = -5.0 * kg_289[k]
                   + f_0 * mg_514[k];

        t_395[k] = -5.0 * kg_290[k]
                   + f_0 * mg_515[k];

        t_396[k] = -5.0 * kg_291[k]
                   + f_0 * mg_516[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, kg_292, kg_293, kg_294, kg_295, \
                         kg_296, mg_517, mg_518, mg_519, mg_520, \
                         mg_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -5.0 * kg_292[k]
                   + f_0 * mg_517[k];

        t_398[k] = -5.0 * kg_293[k]
                   + f_0 * mg_518[k];

        t_399[k] = -5.0 * kg_294[k]
                   + f_0 * mg_519[k];

        t_400[k] = -5.0 * kg_295[k]
                   + f_0 * mg_520[k];

        t_401[k] = -5.0 * kg_296[k]
                   + f_0 * mg_521[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, kg_297, kg_298, kg_299, kg_300, \
                         kg_301, mg_522, mg_523, mg_524, mg_525, \
                         mg_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -5.0 * kg_297[k]
                   + f_0 * mg_522[k];

        t_403[k] = -5.0 * kg_298[k]
                   + f_0 * mg_523[k];

        t_404[k] = -5.0 * kg_299[k]
                   + f_0 * mg_524[k];

        t_405[k] = -6.0 * kg_300[k]
                   + f_0 * mg_525[k];

        t_406[k] = -6.0 * kg_301[k]
                   + f_0 * mg_526[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, kg_302, kg_303, kg_304, kg_305, \
                         kg_306, mg_527, mg_528, mg_529, mg_530, \
                         mg_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -6.0 * kg_302[k]
                   + f_0 * mg_527[k];

        t_408[k] = -6.0 * kg_303[k]
                   + f_0 * mg_528[k];

        t_409[k] = -6.0 * kg_304[k]
                   + f_0 * mg_529[k];

        t_410[k] = -6.0 * kg_305[k]
                   + f_0 * mg_530[k];

        t_411[k] = -6.0 * kg_306[k]
                   + f_0 * mg_531[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, kg_307, kg_308, kg_309, kg_310, \
                         kg_311, mg_532, mg_533, mg_534, mg_535, \
                         mg_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -6.0 * kg_307[k]
                   + f_0 * mg_532[k];

        t_413[k] = -6.0 * kg_308[k]
                   + f_0 * mg_533[k];

        t_414[k] = -6.0 * kg_309[k]
                   + f_0 * mg_534[k];

        t_415[k] = -6.0 * kg_310[k]
                   + f_0 * mg_535[k];

        t_416[k] = -6.0 * kg_311[k]
                   + f_0 * mg_536[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, t_422, kg_312, kg_313, kg_314, \
                         mg_537, mg_538, mg_539, mg_555, mg_556, \
                         mg_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -6.0 * kg_312[k]
                   + f_0 * mg_537[k];

        t_418[k] = -6.0 * kg_313[k]
                   + f_0 * mg_538[k];

        t_419[k] = -6.0 * kg_314[k]
                   + f_0 * mg_539[k];

        t_420[k] = f_0 * mg_555[k];

        t_421[k] = f_0 * mg_556[k];

        t_422[k] = f_0 * mg_557[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, t_428, t_429, t_430, mg_558, \
                         mg_559, mg_560, mg_561, mg_562, mg_563, mg_564, \
                         mg_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_0 * mg_558[k];

        t_424[k] = f_0 * mg_559[k];

        t_425[k] = f_0 * mg_560[k];

        t_426[k] = f_0 * mg_561[k];

        t_427[k] = f_0 * mg_562[k];

        t_428[k] = f_0 * mg_563[k];

        t_429[k] = f_0 * mg_564[k];

        t_430[k] = f_0 * mg_565[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, t_436, kg_315, kg_316, mg_566, \
                         mg_567, mg_568, mg_569, mg_570, mg_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_0 * mg_566[k];

        t_432[k] = f_0 * mg_567[k];

        t_433[k] = f_0 * mg_568[k];

        t_434[k] = f_0 * mg_569[k];

        t_435[k] = -kg_315[k]
                   + f_0 * mg_570[k];

        t_436[k] = -kg_316[k]
                   + f_0 * mg_571[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, t_441, kg_317, kg_318, kg_319, kg_320, \
                         kg_321, mg_572, mg_573, mg_574, mg_575, \
                         mg_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = -kg_317[k]
                   + f_0 * mg_572[k];

        t_438[k] = -kg_318[k]
                   + f_0 * mg_573[k];

        t_439[k] = -kg_319[k]
                   + f_0 * mg_574[k];

        t_440[k] = -kg_320[k]
                   + f_0 * mg_575[k];

        t_441[k] = -kg_321[k]
                   + f_0 * mg_576[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, t_446, kg_322, kg_323, kg_324, kg_325, \
                         kg_326, mg_577, mg_578, mg_579, mg_580, \
                         mg_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = -kg_322[k]
                   + f_0 * mg_577[k];

        t_443[k] = -kg_323[k]
                   + f_0 * mg_578[k];

        t_444[k] = -kg_324[k]
                   + f_0 * mg_579[k];

        t_445[k] = -kg_325[k]
                   + f_0 * mg_580[k];

        t_446[k] = -kg_326[k]
                   + f_0 * mg_581[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, t_451, kg_327, kg_328, kg_329, kg_330, \
                         kg_331, mg_582, mg_583, mg_584, mg_585, \
                         mg_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = -kg_327[k]
                   + f_0 * mg_582[k];

        t_448[k] = -kg_328[k]
                   + f_0 * mg_583[k];

        t_449[k] = -kg_329[k]
                   + f_0 * mg_584[k];

        t_450[k] = -2.0 * kg_330[k]
                   + f_0 * mg_585[k];

        t_451[k] = -2.0 * kg_331[k]
                   + f_0 * mg_586[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, kg_332, kg_333, kg_334, kg_335, \
                         kg_336, mg_587, mg_588, mg_589, mg_590, \
                         mg_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = -2.0 * kg_332[k]
                   + f_0 * mg_587[k];

        t_453[k] = -2.0 * kg_333[k]
                   + f_0 * mg_588[k];

        t_454[k] = -2.0 * kg_334[k]
                   + f_0 * mg_589[k];

        t_455[k] = -2.0 * kg_335[k]
                   + f_0 * mg_590[k];

        t_456[k] = -2.0 * kg_336[k]
                   + f_0 * mg_591[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, t_461, kg_337, kg_338, kg_339, kg_340, \
                         kg_341, mg_592, mg_593, mg_594, mg_595, \
                         mg_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = -2.0 * kg_337[k]
                   + f_0 * mg_592[k];

        t_458[k] = -2.0 * kg_338[k]
                   + f_0 * mg_593[k];

        t_459[k] = -2.0 * kg_339[k]
                   + f_0 * mg_594[k];

        t_460[k] = -2.0 * kg_340[k]
                   + f_0 * mg_595[k];

        t_461[k] = -2.0 * kg_341[k]
                   + f_0 * mg_596[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, kg_342, kg_343, kg_344, kg_345, \
                         kg_346, mg_597, mg_598, mg_599, mg_600, \
                         mg_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = -2.0 * kg_342[k]
                   + f_0 * mg_597[k];

        t_463[k] = -2.0 * kg_343[k]
                   + f_0 * mg_598[k];

        t_464[k] = -2.0 * kg_344[k]
                   + f_0 * mg_599[k];

        t_465[k] = -3.0 * kg_345[k]
                   + f_0 * mg_600[k];

        t_466[k] = -3.0 * kg_346[k]
                   + f_0 * mg_601[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, kg_347, kg_348, kg_349, kg_350, \
                         kg_351, mg_602, mg_603, mg_604, mg_605, \
                         mg_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = -3.0 * kg_347[k]
                   + f_0 * mg_602[k];

        t_468[k] = -3.0 * kg_348[k]
                   + f_0 * mg_603[k];

        t_469[k] = -3.0 * kg_349[k]
                   + f_0 * mg_604[k];

        t_470[k] = -3.0 * kg_350[k]
                   + f_0 * mg_605[k];

        t_471[k] = -3.0 * kg_351[k]
                   + f_0 * mg_606[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, kg_352, kg_353, kg_354, kg_355, \
                         kg_356, mg_607, mg_608, mg_609, mg_610, \
                         mg_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = -3.0 * kg_352[k]
                   + f_0 * mg_607[k];

        t_473[k] = -3.0 * kg_353[k]
                   + f_0 * mg_608[k];

        t_474[k] = -3.0 * kg_354[k]
                   + f_0 * mg_609[k];

        t_475[k] = -3.0 * kg_355[k]
                   + f_0 * mg_610[k];

        t_476[k] = -3.0 * kg_356[k]
                   + f_0 * mg_611[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, kg_357, kg_358, kg_359, kg_360, \
                         kg_361, mg_612, mg_613, mg_614, mg_615, \
                         mg_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = -3.0 * kg_357[k]
                   + f_0 * mg_612[k];

        t_478[k] = -3.0 * kg_358[k]
                   + f_0 * mg_613[k];

        t_479[k] = -3.0 * kg_359[k]
                   + f_0 * mg_614[k];

        t_480[k] = -4.0 * kg_360[k]
                   + f_0 * mg_615[k];

        t_481[k] = -4.0 * kg_361[k]
                   + f_0 * mg_616[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, t_486, kg_362, kg_363, kg_364, kg_365, \
                         kg_366, mg_617, mg_618, mg_619, mg_620, \
                         mg_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = -4.0 * kg_362[k]
                   + f_0 * mg_617[k];

        t_483[k] = -4.0 * kg_363[k]
                   + f_0 * mg_618[k];

        t_484[k] = -4.0 * kg_364[k]
                   + f_0 * mg_619[k];

        t_485[k] = -4.0 * kg_365[k]
                   + f_0 * mg_620[k];

        t_486[k] = -4.0 * kg_366[k]
                   + f_0 * mg_621[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *mg_622 = buffer.data(mg + 622);
    const auto *mg_623 = buffer.data(mg + 623);
    const auto *mg_624 = buffer.data(mg + 624);
    const auto *mg_625 = buffer.data(mg + 625);
    const auto *mg_626 = buffer.data(mg + 626);
    const auto *mg_627 = buffer.data(mg + 627);
    const auto *mg_628 = buffer.data(mg + 628);
    const auto *mg_629 = buffer.data(mg + 629);
    const auto *mg_630 = buffer.data(mg + 630);
    const auto *mg_631 = buffer.data(mg + 631);
    const auto *mg_632 = buffer.data(mg + 632);
    const auto *mg_633 = buffer.data(mg + 633);
    const auto *mg_634 = buffer.data(mg + 634);
    const auto *mg_635 = buffer.data(mg + 635);
    const auto *mg_636 = buffer.data(mg + 636);
    const auto *mg_637 = buffer.data(mg + 637);
    const auto *mg_638 = buffer.data(mg + 638);
    const auto *mg_639 = buffer.data(mg + 639);
    const auto *mg_640 = buffer.data(mg + 640);
    const auto *mg_641 = buffer.data(mg + 641);
    const auto *mg_642 = buffer.data(mg + 642);
    const auto *mg_643 = buffer.data(mg + 643);
    const auto *mg_644 = buffer.data(mg + 644);
    const auto *mg_645 = buffer.data(mg + 645);
    const auto *mg_646 = buffer.data(mg + 646);
    const auto *mg_647 = buffer.data(mg + 647);
    const auto *mg_648 = buffer.data(mg + 648);
    const auto *mg_649 = buffer.data(mg + 649);
    const auto *mg_650 = buffer.data(mg + 650);
    const auto *mg_651 = buffer.data(mg + 651);
    const auto *mg_652 = buffer.data(mg + 652);
    const auto *mg_653 = buffer.data(mg + 653);
    const auto *mg_654 = buffer.data(mg + 654);
    const auto *mg_655 = buffer.data(mg + 655);
    const auto *mg_656 = buffer.data(mg + 656);
    const auto *mg_657 = buffer.data(mg + 657);
    const auto *mg_658 = buffer.data(mg + 658);
    const auto *mg_659 = buffer.data(mg + 659);
    const auto *mg_660 = buffer.data(mg + 660);
    const auto *mg_661 = buffer.data(mg + 661);
    const auto *mg_662 = buffer.data(mg + 662);
    const auto *mg_663 = buffer.data(mg + 663);
    const auto *mg_664 = buffer.data(mg + 664);
    const auto *mg_665 = buffer.data(mg + 665);
    const auto *mg_666 = buffer.data(mg + 666);
    const auto *mg_667 = buffer.data(mg + 667);
    const auto *mg_668 = buffer.data(mg + 668);
    const auto *mg_669 = buffer.data(mg + 669);
    const auto *mg_670 = buffer.data(mg + 670);
    const auto *mg_671 = buffer.data(mg + 671);
    const auto *mg_672 = buffer.data(mg + 672);
    const auto *mg_673 = buffer.data(mg + 673);
    const auto *mg_674 = buffer.data(mg + 674);
    const auto *mg_690 = buffer.data(mg + 690);
    const auto *mg_691 = buffer.data(mg + 691);
    const auto *mg_692 = buffer.data(mg + 692);
    const auto *mg_693 = buffer.data(mg + 693);
    const auto *mg_694 = buffer.data(mg + 694);
    const auto *mg_695 = buffer.data(mg + 695);
    const auto *mg_696 = buffer.data(mg + 696);
    const auto *mg_697 = buffer.data(mg + 697);
    const auto *mg_698 = buffer.data(mg + 698);
    const auto *mg_699 = buffer.data(mg + 699);
    const auto *mg_700 = buffer.data(mg + 700);
    const auto *mg_701 = buffer.data(mg + 701);
    const auto *mg_702 = buffer.data(mg + 702);
    const auto *mg_703 = buffer.data(mg + 703);
    const auto *mg_704 = buffer.data(mg + 704);
    const auto *mg_705 = buffer.data(mg + 705);
    const auto *mg_706 = buffer.data(mg + 706);
    const auto *mg_707 = buffer.data(mg + 707);
    const auto *mg_708 = buffer.data(mg + 708);
    const auto *mg_709 = buffer.data(mg + 709);
    const auto *mg_710 = buffer.data(mg + 710);
    const auto *mg_711 = buffer.data(mg + 711);
    const auto *mg_712 = buffer.data(mg + 712);
    const auto *mg_713 = buffer.data(mg + 713);
    const auto *mg_714 = buffer.data(mg + 714);
    const auto *mg_715 = buffer.data(mg + 715);
    const auto *mg_716 = buffer.data(mg + 716);
    const auto *mg_717 = buffer.data(mg + 717);
    const auto *mg_718 = buffer.data(mg + 718);
    const auto *mg_719 = buffer.data(mg + 719);
    const auto *mg_720 = buffer.data(mg + 720);
    const auto *mg_721 = buffer.data(mg + 721);
    const auto *mg_722 = buffer.data(mg + 722);
    const auto *mg_723 = buffer.data(mg + 723);
    const auto *mg_724 = buffer.data(mg + 724);
    const auto *mg_725 = buffer.data(mg + 725);
    const auto *mg_726 = buffer.data(mg + 726);
    const auto *mg_727 = buffer.data(mg + 727);
    const auto *mg_728 = buffer.data(mg + 728);
    const auto *mg_729 = buffer.data(mg + 729);
    const auto *mg_730 = buffer.data(mg + 730);
    const auto *mg_731 = buffer.data(mg + 731);
    const auto *mg_732 = buffer.data(mg + 732);
    const auto *mg_733 = buffer.data(mg + 733);
    const auto *mg_734 = buffer.data(mg + 734);
    const auto *mg_735 = buffer.data(mg + 735);
    const auto *mg_736 = buffer.data(mg + 736);
    const auto *mg_737 = buffer.data(mg + 737);
    const auto *mg_738 = buffer.data(mg + 738);
    const auto *mg_739 = buffer.data(mg + 739);
    const auto *mg_740 = buffer.data(mg + 740);
    const auto *mg_741 = buffer.data(mg + 741);
    const auto *mg_742 = buffer.data(mg + 742);
    const auto *mg_743 = buffer.data(mg + 743);
    const auto *mg_744 = buffer.data(mg + 744);
    const auto *mg_745 = buffer.data(mg + 745);
    const auto *mg_746 = buffer.data(mg + 746);
    const auto *mg_747 = buffer.data(mg + 747);
    const auto *mg_748 = buffer.data(mg + 748);
    const auto *mg_749 = buffer.data(mg + 749);
    const auto *mg_750 = buffer.data(mg + 750);
    const auto *mg_751 = buffer.data(mg + 751);
    const auto *mg_752 = buffer.data(mg + 752);
    const auto *mg_753 = buffer.data(mg + 753);
    const auto *mg_754 = buffer.data(mg + 754);
    const auto *mg_755 = buffer.data(mg + 755);
    const auto *mg_756 = buffer.data(mg + 756);
    const auto *mg_757 = buffer.data(mg + 757);
    const auto *mg_758 = buffer.data(mg + 758);
    const auto *mg_759 = buffer.data(mg + 759);
    const auto *mg_760 = buffer.data(mg + 760);
    const auto *mg_761 = buffer.data(mg + 761);
    const auto *mg_762 = buffer.data(mg + 762);
    const auto *mg_763 = buffer.data(mg + 763);
    const auto *mg_764 = buffer.data(mg + 764);
    const auto *mg_765 = buffer.data(mg + 765);
    const auto *mg_766 = buffer.data(mg + 766);
    const auto *mg_767 = buffer.data(mg + 767);
    const auto *mg_768 = buffer.data(mg + 768);
    const auto *mg_769 = buffer.data(mg + 769);
    const auto *mg_770 = buffer.data(mg + 770);
    const auto *mg_771 = buffer.data(mg + 771);
    const auto *mg_772 = buffer.data(mg + 772);
    const auto *mg_773 = buffer.data(mg + 773);
    const auto *mg_774 = buffer.data(mg + 774);
    const auto *mg_775 = buffer.data(mg + 775);
    const auto *mg_776 = buffer.data(mg + 776);
    const auto *mg_777 = buffer.data(mg + 777);
    const auto *mg_778 = buffer.data(mg + 778);
    const auto *mg_779 = buffer.data(mg + 779);
    const auto *mg_780 = buffer.data(mg + 780);
    const auto *mg_781 = buffer.data(mg + 781);
    const auto *mg_782 = buffer.data(mg + 782);
    const auto *mg_783 = buffer.data(mg + 783);
    const auto *mg_784 = buffer.data(mg + 784);
    const auto *mg_785 = buffer.data(mg + 785);
    const auto *mg_786 = buffer.data(mg + 786);
    const auto *mg_787 = buffer.data(mg + 787);
    const auto *mg_788 = buffer.data(mg + 788);
    const auto *mg_789 = buffer.data(mg + 789);
    const auto *mg_790 = buffer.data(mg + 790);
    const auto *mg_791 = buffer.data(mg + 791);

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, kg_367, kg_368, kg_369, kg_370, \
                         kg_371, mg_622, mg_623, mg_624, mg_625, \
                         mg_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = -4.0 * kg_367[k]
                   + f_0 * mg_622[k];

        t_488[k] = -4.0 * kg_368[k]
                   + f_0 * mg_623[k];

        t_489[k] = -4.0 * kg_369[k]
                   + f_0 * mg_624[k];

        t_490[k] = -4.0 * kg_370[k]
                   + f_0 * mg_625[k];

        t_491[k] = -4.0 * kg_371[k]
                   + f_0 * mg_626[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, kg_372, kg_373, kg_374, kg_375, \
                         kg_376, mg_627, mg_628, mg_629, mg_630, \
                         mg_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = -4.0 * kg_372[k]
                   + f_0 * mg_627[k];

        t_493[k] = -4.0 * kg_373[k]
                   + f_0 * mg_628[k];

        t_494[k] = -4.0 * kg_374[k]
                   + f_0 * mg_629[k];

        t_495[k] = -5.0 * kg_375[k]
                   + f_0 * mg_630[k];

        t_496[k] = -5.0 * kg_376[k]
                   + f_0 * mg_631[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, t_500, t_501, kg_377, kg_378, kg_379, kg_380, \
                         kg_381, mg_632, mg_633, mg_634, mg_635, \
                         mg_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = -5.0 * kg_377[k]
                   + f_0 * mg_632[k];

        t_498[k] = -5.0 * kg_378[k]
                   + f_0 * mg_633[k];

        t_499[k] = -5.0 * kg_379[k]
                   + f_0 * mg_634[k];

        t_500[k] = -5.0 * kg_380[k]
                   + f_0 * mg_635[k];

        t_501[k] = -5.0 * kg_381[k]
                   + f_0 * mg_636[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, kg_382, kg_383, kg_384, kg_385, \
                         kg_386, mg_637, mg_638, mg_639, mg_640, \
                         mg_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = -5.0 * kg_382[k]
                   + f_0 * mg_637[k];

        t_503[k] = -5.0 * kg_383[k]
                   + f_0 * mg_638[k];

        t_504[k] = -5.0 * kg_384[k]
                   + f_0 * mg_639[k];

        t_505[k] = -5.0 * kg_385[k]
                   + f_0 * mg_640[k];

        t_506[k] = -5.0 * kg_386[k]
                   + f_0 * mg_641[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, t_510, t_511, kg_387, kg_388, kg_389, kg_390, \
                         kg_391, mg_642, mg_643, mg_644, mg_645, \
                         mg_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = -5.0 * kg_387[k]
                   + f_0 * mg_642[k];

        t_508[k] = -5.0 * kg_388[k]
                   + f_0 * mg_643[k];

        t_509[k] = -5.0 * kg_389[k]
                   + f_0 * mg_644[k];

        t_510[k] = -6.0 * kg_390[k]
                   + f_0 * mg_645[k];

        t_511[k] = -6.0 * kg_391[k]
                   + f_0 * mg_646[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, kg_392, kg_393, kg_394, kg_395, \
                         kg_396, mg_647, mg_648, mg_649, mg_650, \
                         mg_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = -6.0 * kg_392[k]
                   + f_0 * mg_647[k];

        t_513[k] = -6.0 * kg_393[k]
                   + f_0 * mg_648[k];

        t_514[k] = -6.0 * kg_394[k]
                   + f_0 * mg_649[k];

        t_515[k] = -6.0 * kg_395[k]
                   + f_0 * mg_650[k];

        t_516[k] = -6.0 * kg_396[k]
                   + f_0 * mg_651[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, kg_397, kg_398, kg_399, kg_400, \
                         kg_401, mg_652, mg_653, mg_654, mg_655, \
                         mg_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = -6.0 * kg_397[k]
                   + f_0 * mg_652[k];

        t_518[k] = -6.0 * kg_398[k]
                   + f_0 * mg_653[k];

        t_519[k] = -6.0 * kg_399[k]
                   + f_0 * mg_654[k];

        t_520[k] = -6.0 * kg_400[k]
                   + f_0 * mg_655[k];

        t_521[k] = -6.0 * kg_401[k]
                   + f_0 * mg_656[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, kg_402, kg_403, kg_404, kg_405, \
                         kg_406, mg_657, mg_658, mg_659, mg_660, \
                         mg_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = -6.0 * kg_402[k]
                   + f_0 * mg_657[k];

        t_523[k] = -6.0 * kg_403[k]
                   + f_0 * mg_658[k];

        t_524[k] = -6.0 * kg_404[k]
                   + f_0 * mg_659[k];

        t_525[k] = -7.0 * kg_405[k]
                   + f_0 * mg_660[k];

        t_526[k] = -7.0 * kg_406[k]
                   + f_0 * mg_661[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, kg_407, kg_408, kg_409, kg_410, \
                         kg_411, mg_662, mg_663, mg_664, mg_665, \
                         mg_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = -7.0 * kg_407[k]
                   + f_0 * mg_662[k];

        t_528[k] = -7.0 * kg_408[k]
                   + f_0 * mg_663[k];

        t_529[k] = -7.0 * kg_409[k]
                   + f_0 * mg_664[k];

        t_530[k] = -7.0 * kg_410[k]
                   + f_0 * mg_665[k];

        t_531[k] = -7.0 * kg_411[k]
                   + f_0 * mg_666[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, kg_412, kg_413, kg_414, kg_415, \
                         kg_416, mg_667, mg_668, mg_669, mg_670, \
                         mg_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = -7.0 * kg_412[k]
                   + f_0 * mg_667[k];

        t_533[k] = -7.0 * kg_413[k]
                   + f_0 * mg_668[k];

        t_534[k] = -7.0 * kg_414[k]
                   + f_0 * mg_669[k];

        t_535[k] = -7.0 * kg_415[k]
                   + f_0 * mg_670[k];

        t_536[k] = -7.0 * kg_416[k]
                   + f_0 * mg_671[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, t_541, t_542, kg_417, kg_418, kg_419, \
                         mg_672, mg_673, mg_674, mg_690, mg_691, \
                         mg_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = -7.0 * kg_417[k]
                   + f_0 * mg_672[k];

        t_538[k] = -7.0 * kg_418[k]
                   + f_0 * mg_673[k];

        t_539[k] = -7.0 * kg_419[k]
                   + f_0 * mg_674[k];

        t_540[k] = f_0 * mg_690[k];

        t_541[k] = f_0 * mg_691[k];

        t_542[k] = f_0 * mg_692[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, t_547, t_548, t_549, t_550, mg_693, \
                         mg_694, mg_695, mg_696, mg_697, mg_698, mg_699, \
                         mg_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_0 * mg_693[k];

        t_544[k] = f_0 * mg_694[k];

        t_545[k] = f_0 * mg_695[k];

        t_546[k] = f_0 * mg_696[k];

        t_547[k] = f_0 * mg_697[k];

        t_548[k] = f_0 * mg_698[k];

        t_549[k] = f_0 * mg_699[k];

        t_550[k] = f_0 * mg_700[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, t_554, t_555, t_556, kg_420, kg_421, mg_701, \
                         mg_702, mg_703, mg_704, mg_705, mg_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_0 * mg_701[k];

        t_552[k] = f_0 * mg_702[k];

        t_553[k] = f_0 * mg_703[k];

        t_554[k] = f_0 * mg_704[k];

        t_555[k] = -kg_420[k]
                   + f_0 * mg_705[k];

        t_556[k] = -kg_421[k]
                   + f_0 * mg_706[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, kg_422, kg_423, kg_424, kg_425, \
                         kg_426, mg_707, mg_708, mg_709, mg_710, \
                         mg_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = -kg_422[k]
                   + f_0 * mg_707[k];

        t_558[k] = -kg_423[k]
                   + f_0 * mg_708[k];

        t_559[k] = -kg_424[k]
                   + f_0 * mg_709[k];

        t_560[k] = -kg_425[k]
                   + f_0 * mg_710[k];

        t_561[k] = -kg_426[k]
                   + f_0 * mg_711[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, kg_427, kg_428, kg_429, kg_430, \
                         kg_431, mg_712, mg_713, mg_714, mg_715, \
                         mg_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = -kg_427[k]
                   + f_0 * mg_712[k];

        t_563[k] = -kg_428[k]
                   + f_0 * mg_713[k];

        t_564[k] = -kg_429[k]
                   + f_0 * mg_714[k];

        t_565[k] = -kg_430[k]
                   + f_0 * mg_715[k];

        t_566[k] = -kg_431[k]
                   + f_0 * mg_716[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, kg_432, kg_433, kg_434, kg_435, \
                         kg_436, mg_717, mg_718, mg_719, mg_720, \
                         mg_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -kg_432[k]
                   + f_0 * mg_717[k];

        t_568[k] = -kg_433[k]
                   + f_0 * mg_718[k];

        t_569[k] = -kg_434[k]
                   + f_0 * mg_719[k];

        t_570[k] = -2.0 * kg_435[k]
                   + f_0 * mg_720[k];

        t_571[k] = -2.0 * kg_436[k]
                   + f_0 * mg_721[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, kg_437, kg_438, kg_439, kg_440, \
                         kg_441, mg_722, mg_723, mg_724, mg_725, \
                         mg_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = -2.0 * kg_437[k]
                   + f_0 * mg_722[k];

        t_573[k] = -2.0 * kg_438[k]
                   + f_0 * mg_723[k];

        t_574[k] = -2.0 * kg_439[k]
                   + f_0 * mg_724[k];

        t_575[k] = -2.0 * kg_440[k]
                   + f_0 * mg_725[k];

        t_576[k] = -2.0 * kg_441[k]
                   + f_0 * mg_726[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, kg_442, kg_443, kg_444, kg_445, \
                         kg_446, mg_727, mg_728, mg_729, mg_730, \
                         mg_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = -2.0 * kg_442[k]
                   + f_0 * mg_727[k];

        t_578[k] = -2.0 * kg_443[k]
                   + f_0 * mg_728[k];

        t_579[k] = -2.0 * kg_444[k]
                   + f_0 * mg_729[k];

        t_580[k] = -2.0 * kg_445[k]
                   + f_0 * mg_730[k];

        t_581[k] = -2.0 * kg_446[k]
                   + f_0 * mg_731[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, kg_447, kg_448, kg_449, kg_450, \
                         kg_451, mg_732, mg_733, mg_734, mg_735, \
                         mg_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = -2.0 * kg_447[k]
                   + f_0 * mg_732[k];

        t_583[k] = -2.0 * kg_448[k]
                   + f_0 * mg_733[k];

        t_584[k] = -2.0 * kg_449[k]
                   + f_0 * mg_734[k];

        t_585[k] = -3.0 * kg_450[k]
                   + f_0 * mg_735[k];

        t_586[k] = -3.0 * kg_451[k]
                   + f_0 * mg_736[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, t_591, kg_452, kg_453, kg_454, kg_455, \
                         kg_456, mg_737, mg_738, mg_739, mg_740, \
                         mg_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = -3.0 * kg_452[k]
                   + f_0 * mg_737[k];

        t_588[k] = -3.0 * kg_453[k]
                   + f_0 * mg_738[k];

        t_589[k] = -3.0 * kg_454[k]
                   + f_0 * mg_739[k];

        t_590[k] = -3.0 * kg_455[k]
                   + f_0 * mg_740[k];

        t_591[k] = -3.0 * kg_456[k]
                   + f_0 * mg_741[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, t_596, kg_457, kg_458, kg_459, kg_460, \
                         kg_461, mg_742, mg_743, mg_744, mg_745, \
                         mg_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = -3.0 * kg_457[k]
                   + f_0 * mg_742[k];

        t_593[k] = -3.0 * kg_458[k]
                   + f_0 * mg_743[k];

        t_594[k] = -3.0 * kg_459[k]
                   + f_0 * mg_744[k];

        t_595[k] = -3.0 * kg_460[k]
                   + f_0 * mg_745[k];

        t_596[k] = -3.0 * kg_461[k]
                   + f_0 * mg_746[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, t_600, t_601, kg_462, kg_463, kg_464, kg_465, \
                         kg_466, mg_747, mg_748, mg_749, mg_750, \
                         mg_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = -3.0 * kg_462[k]
                   + f_0 * mg_747[k];

        t_598[k] = -3.0 * kg_463[k]
                   + f_0 * mg_748[k];

        t_599[k] = -3.0 * kg_464[k]
                   + f_0 * mg_749[k];

        t_600[k] = -4.0 * kg_465[k]
                   + f_0 * mg_750[k];

        t_601[k] = -4.0 * kg_466[k]
                   + f_0 * mg_751[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, kg_467, kg_468, kg_469, kg_470, \
                         kg_471, mg_752, mg_753, mg_754, mg_755, \
                         mg_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = -4.0 * kg_467[k]
                   + f_0 * mg_752[k];

        t_603[k] = -4.0 * kg_468[k]
                   + f_0 * mg_753[k];

        t_604[k] = -4.0 * kg_469[k]
                   + f_0 * mg_754[k];

        t_605[k] = -4.0 * kg_470[k]
                   + f_0 * mg_755[k];

        t_606[k] = -4.0 * kg_471[k]
                   + f_0 * mg_756[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, t_611, kg_472, kg_473, kg_474, kg_475, \
                         kg_476, mg_757, mg_758, mg_759, mg_760, \
                         mg_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = -4.0 * kg_472[k]
                   + f_0 * mg_757[k];

        t_608[k] = -4.0 * kg_473[k]
                   + f_0 * mg_758[k];

        t_609[k] = -4.0 * kg_474[k]
                   + f_0 * mg_759[k];

        t_610[k] = -4.0 * kg_475[k]
                   + f_0 * mg_760[k];

        t_611[k] = -4.0 * kg_476[k]
                   + f_0 * mg_761[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, kg_477, kg_478, kg_479, kg_480, \
                         kg_481, mg_762, mg_763, mg_764, mg_765, \
                         mg_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = -4.0 * kg_477[k]
                   + f_0 * mg_762[k];

        t_613[k] = -4.0 * kg_478[k]
                   + f_0 * mg_763[k];

        t_614[k] = -4.0 * kg_479[k]
                   + f_0 * mg_764[k];

        t_615[k] = -5.0 * kg_480[k]
                   + f_0 * mg_765[k];

        t_616[k] = -5.0 * kg_481[k]
                   + f_0 * mg_766[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, kg_482, kg_483, kg_484, kg_485, \
                         kg_486, mg_767, mg_768, mg_769, mg_770, \
                         mg_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = -5.0 * kg_482[k]
                   + f_0 * mg_767[k];

        t_618[k] = -5.0 * kg_483[k]
                   + f_0 * mg_768[k];

        t_619[k] = -5.0 * kg_484[k]
                   + f_0 * mg_769[k];

        t_620[k] = -5.0 * kg_485[k]
                   + f_0 * mg_770[k];

        t_621[k] = -5.0 * kg_486[k]
                   + f_0 * mg_771[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, kg_487, kg_488, kg_489, kg_490, \
                         kg_491, mg_772, mg_773, mg_774, mg_775, \
                         mg_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = -5.0 * kg_487[k]
                   + f_0 * mg_772[k];

        t_623[k] = -5.0 * kg_488[k]
                   + f_0 * mg_773[k];

        t_624[k] = -5.0 * kg_489[k]
                   + f_0 * mg_774[k];

        t_625[k] = -5.0 * kg_490[k]
                   + f_0 * mg_775[k];

        t_626[k] = -5.0 * kg_491[k]
                   + f_0 * mg_776[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, t_631, kg_492, kg_493, kg_494, kg_495, \
                         kg_496, mg_777, mg_778, mg_779, mg_780, \
                         mg_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = -5.0 * kg_492[k]
                   + f_0 * mg_777[k];

        t_628[k] = -5.0 * kg_493[k]
                   + f_0 * mg_778[k];

        t_629[k] = -5.0 * kg_494[k]
                   + f_0 * mg_779[k];

        t_630[k] = -6.0 * kg_495[k]
                   + f_0 * mg_780[k];

        t_631[k] = -6.0 * kg_496[k]
                   + f_0 * mg_781[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, t_636, kg_497, kg_498, kg_499, kg_500, \
                         kg_501, mg_782, mg_783, mg_784, mg_785, \
                         mg_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = -6.0 * kg_497[k]
                   + f_0 * mg_782[k];

        t_633[k] = -6.0 * kg_498[k]
                   + f_0 * mg_783[k];

        t_634[k] = -6.0 * kg_499[k]
                   + f_0 * mg_784[k];

        t_635[k] = -6.0 * kg_500[k]
                   + f_0 * mg_785[k];

        t_636[k] = -6.0 * kg_501[k]
                   + f_0 * mg_786[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, t_641, kg_502, kg_503, kg_504, kg_505, \
                         kg_506, mg_787, mg_788, mg_789, mg_790, \
                         mg_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = -6.0 * kg_502[k]
                   + f_0 * mg_787[k];

        t_638[k] = -6.0 * kg_503[k]
                   + f_0 * mg_788[k];

        t_639[k] = -6.0 * kg_504[k]
                   + f_0 * mg_789[k];

        t_640[k] = -6.0 * kg_505[k]
                   + f_0 * mg_790[k];

        t_641[k] = -6.0 * kg_506[k]
                   + f_0 * mg_791[k];
    }
}

static auto
compute_prim_geom_10_lg_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kg, const size_t mg,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *mg_792 = buffer.data(mg + 792);
    const auto *mg_793 = buffer.data(mg + 793);
    const auto *mg_794 = buffer.data(mg + 794);
    const auto *mg_795 = buffer.data(mg + 795);
    const auto *mg_796 = buffer.data(mg + 796);
    const auto *mg_797 = buffer.data(mg + 797);
    const auto *mg_798 = buffer.data(mg + 798);
    const auto *mg_799 = buffer.data(mg + 799);
    const auto *mg_800 = buffer.data(mg + 800);
    const auto *mg_801 = buffer.data(mg + 801);
    const auto *mg_802 = buffer.data(mg + 802);
    const auto *mg_803 = buffer.data(mg + 803);
    const auto *mg_804 = buffer.data(mg + 804);
    const auto *mg_805 = buffer.data(mg + 805);
    const auto *mg_806 = buffer.data(mg + 806);
    const auto *mg_807 = buffer.data(mg + 807);
    const auto *mg_808 = buffer.data(mg + 808);
    const auto *mg_809 = buffer.data(mg + 809);
    const auto *mg_810 = buffer.data(mg + 810);
    const auto *mg_811 = buffer.data(mg + 811);
    const auto *mg_812 = buffer.data(mg + 812);
    const auto *mg_813 = buffer.data(mg + 813);
    const auto *mg_814 = buffer.data(mg + 814);
    const auto *mg_815 = buffer.data(mg + 815);
    const auto *mg_816 = buffer.data(mg + 816);
    const auto *mg_817 = buffer.data(mg + 817);
    const auto *mg_818 = buffer.data(mg + 818);
    const auto *mg_819 = buffer.data(mg + 819);
    const auto *mg_820 = buffer.data(mg + 820);
    const auto *mg_821 = buffer.data(mg + 821);
    const auto *mg_822 = buffer.data(mg + 822);
    const auto *mg_823 = buffer.data(mg + 823);
    const auto *mg_824 = buffer.data(mg + 824);

#pragma omp simd aligned(t_642, t_643, t_644, t_645, t_646, kg_507, kg_508, kg_509, kg_510, \
                         kg_511, mg_792, mg_793, mg_794, mg_795, \
                         mg_796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = -6.0 * kg_507[k]
                   + f_0 * mg_792[k];

        t_643[k] = -6.0 * kg_508[k]
                   + f_0 * mg_793[k];

        t_644[k] = -6.0 * kg_509[k]
                   + f_0 * mg_794[k];

        t_645[k] = -7.0 * kg_510[k]
                   + f_0 * mg_795[k];

        t_646[k] = -7.0 * kg_511[k]
                   + f_0 * mg_796[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, t_651, kg_512, kg_513, kg_514, kg_515, \
                         kg_516, mg_797, mg_798, mg_799, mg_800, \
                         mg_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = -7.0 * kg_512[k]
                   + f_0 * mg_797[k];

        t_648[k] = -7.0 * kg_513[k]
                   + f_0 * mg_798[k];

        t_649[k] = -7.0 * kg_514[k]
                   + f_0 * mg_799[k];

        t_650[k] = -7.0 * kg_515[k]
                   + f_0 * mg_800[k];

        t_651[k] = -7.0 * kg_516[k]
                   + f_0 * mg_801[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, t_656, kg_517, kg_518, kg_519, kg_520, \
                         kg_521, mg_802, mg_803, mg_804, mg_805, \
                         mg_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = -7.0 * kg_517[k]
                   + f_0 * mg_802[k];

        t_653[k] = -7.0 * kg_518[k]
                   + f_0 * mg_803[k];

        t_654[k] = -7.0 * kg_519[k]
                   + f_0 * mg_804[k];

        t_655[k] = -7.0 * kg_520[k]
                   + f_0 * mg_805[k];

        t_656[k] = -7.0 * kg_521[k]
                   + f_0 * mg_806[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, kg_522, kg_523, kg_524, kg_525, \
                         kg_526, mg_807, mg_808, mg_809, mg_810, \
                         mg_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = -7.0 * kg_522[k]
                   + f_0 * mg_807[k];

        t_658[k] = -7.0 * kg_523[k]
                   + f_0 * mg_808[k];

        t_659[k] = -7.0 * kg_524[k]
                   + f_0 * mg_809[k];

        t_660[k] = -8.0 * kg_525[k]
                   + f_0 * mg_810[k];

        t_661[k] = -8.0 * kg_526[k]
                   + f_0 * mg_811[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, t_666, kg_527, kg_528, kg_529, kg_530, \
                         kg_531, mg_812, mg_813, mg_814, mg_815, \
                         mg_816 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = -8.0 * kg_527[k]
                   + f_0 * mg_812[k];

        t_663[k] = -8.0 * kg_528[k]
                   + f_0 * mg_813[k];

        t_664[k] = -8.0 * kg_529[k]
                   + f_0 * mg_814[k];

        t_665[k] = -8.0 * kg_530[k]
                   + f_0 * mg_815[k];

        t_666[k] = -8.0 * kg_531[k]
                   + f_0 * mg_816[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, t_670, t_671, kg_532, kg_533, kg_534, kg_535, \
                         kg_536, mg_817, mg_818, mg_819, mg_820, \
                         mg_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = -8.0 * kg_532[k]
                   + f_0 * mg_817[k];

        t_668[k] = -8.0 * kg_533[k]
                   + f_0 * mg_818[k];

        t_669[k] = -8.0 * kg_534[k]
                   + f_0 * mg_819[k];

        t_670[k] = -8.0 * kg_535[k]
                   + f_0 * mg_820[k];

        t_671[k] = -8.0 * kg_536[k]
                   + f_0 * mg_821[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, kg_537, kg_538, kg_539, mg_822, mg_823, \
                         mg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = -8.0 * kg_537[k]
                   + f_0 * mg_822[k];

        t_673[k] = -8.0 * kg_538[k]
                   + f_0 * mg_823[k];

        t_674[k] = -8.0 * kg_539[k]
                   + f_0 * mg_824[k];
    }
}

auto
compute_prim_geom_10_lg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t kg, const size_t mg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_lg_electron_repulsion_2_piece0(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_2_piece1(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_2_piece2(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_2_piece3(buffer, target, kg, mg, ncols, alpha);

    compute_prim_geom_10_lg_electron_repulsion_2_piece4(buffer, target, kg, mg, ncols, alpha);
}

}  // namespace simdt2ceri
