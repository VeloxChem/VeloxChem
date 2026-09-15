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


#include "SimdElectronRepulsionGeom10VrrRecKG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_kg_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ig, const size_t lg,
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

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_1 = buffer.data(lg + 1);
    const auto *lg_2 = buffer.data(lg + 2);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_4 = buffer.data(lg + 4);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_6 = buffer.data(lg + 6);
    const auto *lg_7 = buffer.data(lg + 7);
    const auto *lg_8 = buffer.data(lg + 8);
    const auto *lg_9 = buffer.data(lg + 9);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_11 = buffer.data(lg + 11);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_13 = buffer.data(lg + 13);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_17 = buffer.data(lg + 17);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_19 = buffer.data(lg + 19);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_22 = buffer.data(lg + 22);
    const auto *lg_23 = buffer.data(lg + 23);
    const auto *lg_24 = buffer.data(lg + 24);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_26 = buffer.data(lg + 26);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_31 = buffer.data(lg + 31);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_34 = buffer.data(lg + 34);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_36 = buffer.data(lg + 36);
    const auto *lg_37 = buffer.data(lg + 37);
    const auto *lg_38 = buffer.data(lg + 38);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_49 = buffer.data(lg + 49);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_52 = buffer.data(lg + 52);
    const auto *lg_53 = buffer.data(lg + 53);
    const auto *lg_54 = buffer.data(lg + 54);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_79 = buffer.data(lg + 79);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_81 = buffer.data(lg + 81);
    const auto *lg_82 = buffer.data(lg + 82);
    const auto *lg_83 = buffer.data(lg + 83);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_94 = buffer.data(lg + 94);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_97 = buffer.data(lg + 97);
    const auto *lg_98 = buffer.data(lg + 98);
    const auto *lg_99 = buffer.data(lg + 99);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_139 = buffer.data(lg + 139);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_141 = buffer.data(lg + 141);
    const auto *lg_142 = buffer.data(lg + 142);
    const auto *lg_143 = buffer.data(lg + 143);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ig_0, ig_1, ig_2, ig_3, ig_4, lg_0, lg_1, \
                         lg_2, lg_3, lg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -7.0 * ig_0[k]
                 + f_0 * lg_0[k];

        t_1[k] = -7.0 * ig_1[k]
                 + f_0 * lg_1[k];

        t_2[k] = -7.0 * ig_2[k]
                 + f_0 * lg_2[k];

        t_3[k] = -7.0 * ig_3[k]
                 + f_0 * lg_3[k];

        t_4[k] = -7.0 * ig_4[k]
                 + f_0 * lg_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ig_5, ig_6, ig_7, ig_8, ig_9, lg_5, lg_6, \
                         lg_7, lg_8, lg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -7.0 * ig_5[k]
                 + f_0 * lg_5[k];

        t_6[k] = -7.0 * ig_6[k]
                 + f_0 * lg_6[k];

        t_7[k] = -7.0 * ig_7[k]
                 + f_0 * lg_7[k];

        t_8[k] = -7.0 * ig_8[k]
                 + f_0 * lg_8[k];

        t_9[k] = -7.0 * ig_9[k]
                 + f_0 * lg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ig_10, ig_11, ig_12, ig_13, ig_14, \
                         lg_10, lg_11, lg_12, lg_13, lg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -7.0 * ig_10[k]
                  + f_0 * lg_10[k];

        t_11[k] = -7.0 * ig_11[k]
                  + f_0 * lg_11[k];

        t_12[k] = -7.0 * ig_12[k]
                  + f_0 * lg_12[k];

        t_13[k] = -7.0 * ig_13[k]
                  + f_0 * lg_13[k];

        t_14[k] = -7.0 * ig_14[k]
                  + f_0 * lg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ig_15, ig_16, ig_17, ig_18, ig_19, \
                         lg_15, lg_16, lg_17, lg_18, lg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -6.0 * ig_15[k]
                  + f_0 * lg_15[k];

        t_16[k] = -6.0 * ig_16[k]
                  + f_0 * lg_16[k];

        t_17[k] = -6.0 * ig_17[k]
                  + f_0 * lg_17[k];

        t_18[k] = -6.0 * ig_18[k]
                  + f_0 * lg_18[k];

        t_19[k] = -6.0 * ig_19[k]
                  + f_0 * lg_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ig_20, ig_21, ig_22, ig_23, ig_24, \
                         lg_20, lg_21, lg_22, lg_23, lg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -6.0 * ig_20[k]
                  + f_0 * lg_20[k];

        t_21[k] = -6.0 * ig_21[k]
                  + f_0 * lg_21[k];

        t_22[k] = -6.0 * ig_22[k]
                  + f_0 * lg_22[k];

        t_23[k] = -6.0 * ig_23[k]
                  + f_0 * lg_23[k];

        t_24[k] = -6.0 * ig_24[k]
                  + f_0 * lg_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ig_25, ig_26, ig_27, ig_28, ig_29, \
                         lg_25, lg_26, lg_27, lg_28, lg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -6.0 * ig_25[k]
                  + f_0 * lg_25[k];

        t_26[k] = -6.0 * ig_26[k]
                  + f_0 * lg_26[k];

        t_27[k] = -6.0 * ig_27[k]
                  + f_0 * lg_27[k];

        t_28[k] = -6.0 * ig_28[k]
                  + f_0 * lg_28[k];

        t_29[k] = -6.0 * ig_29[k]
                  + f_0 * lg_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ig_30, ig_31, ig_32, ig_33, ig_34, \
                         lg_30, lg_31, lg_32, lg_33, lg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -6.0 * ig_30[k]
                  + f_0 * lg_30[k];

        t_31[k] = -6.0 * ig_31[k]
                  + f_0 * lg_31[k];

        t_32[k] = -6.0 * ig_32[k]
                  + f_0 * lg_32[k];

        t_33[k] = -6.0 * ig_33[k]
                  + f_0 * lg_33[k];

        t_34[k] = -6.0 * ig_34[k]
                  + f_0 * lg_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ig_35, ig_36, ig_37, ig_38, ig_39, \
                         lg_35, lg_36, lg_37, lg_38, lg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -6.0 * ig_35[k]
                  + f_0 * lg_35[k];

        t_36[k] = -6.0 * ig_36[k]
                  + f_0 * lg_36[k];

        t_37[k] = -6.0 * ig_37[k]
                  + f_0 * lg_37[k];

        t_38[k] = -6.0 * ig_38[k]
                  + f_0 * lg_38[k];

        t_39[k] = -6.0 * ig_39[k]
                  + f_0 * lg_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ig_40, ig_41, ig_42, ig_43, ig_44, \
                         lg_40, lg_41, lg_42, lg_43, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -6.0 * ig_40[k]
                  + f_0 * lg_40[k];

        t_41[k] = -6.0 * ig_41[k]
                  + f_0 * lg_41[k];

        t_42[k] = -6.0 * ig_42[k]
                  + f_0 * lg_42[k];

        t_43[k] = -6.0 * ig_43[k]
                  + f_0 * lg_43[k];

        t_44[k] = -6.0 * ig_44[k]
                  + f_0 * lg_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ig_45, ig_46, ig_47, ig_48, ig_49, \
                         lg_45, lg_46, lg_47, lg_48, lg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -5.0 * ig_45[k]
                  + f_0 * lg_45[k];

        t_46[k] = -5.0 * ig_46[k]
                  + f_0 * lg_46[k];

        t_47[k] = -5.0 * ig_47[k]
                  + f_0 * lg_47[k];

        t_48[k] = -5.0 * ig_48[k]
                  + f_0 * lg_48[k];

        t_49[k] = -5.0 * ig_49[k]
                  + f_0 * lg_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ig_50, ig_51, ig_52, ig_53, ig_54, \
                         lg_50, lg_51, lg_52, lg_53, lg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -5.0 * ig_50[k]
                  + f_0 * lg_50[k];

        t_51[k] = -5.0 * ig_51[k]
                  + f_0 * lg_51[k];

        t_52[k] = -5.0 * ig_52[k]
                  + f_0 * lg_52[k];

        t_53[k] = -5.0 * ig_53[k]
                  + f_0 * lg_53[k];

        t_54[k] = -5.0 * ig_54[k]
                  + f_0 * lg_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ig_55, ig_56, ig_57, ig_58, ig_59, \
                         lg_55, lg_56, lg_57, lg_58, lg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -5.0 * ig_55[k]
                  + f_0 * lg_55[k];

        t_56[k] = -5.0 * ig_56[k]
                  + f_0 * lg_56[k];

        t_57[k] = -5.0 * ig_57[k]
                  + f_0 * lg_57[k];

        t_58[k] = -5.0 * ig_58[k]
                  + f_0 * lg_58[k];

        t_59[k] = -5.0 * ig_59[k]
                  + f_0 * lg_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ig_60, ig_61, ig_62, ig_63, ig_64, \
                         lg_60, lg_61, lg_62, lg_63, lg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -5.0 * ig_60[k]
                  + f_0 * lg_60[k];

        t_61[k] = -5.0 * ig_61[k]
                  + f_0 * lg_61[k];

        t_62[k] = -5.0 * ig_62[k]
                  + f_0 * lg_62[k];

        t_63[k] = -5.0 * ig_63[k]
                  + f_0 * lg_63[k];

        t_64[k] = -5.0 * ig_64[k]
                  + f_0 * lg_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ig_65, ig_66, ig_67, ig_68, ig_69, \
                         lg_65, lg_66, lg_67, lg_68, lg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -5.0 * ig_65[k]
                  + f_0 * lg_65[k];

        t_66[k] = -5.0 * ig_66[k]
                  + f_0 * lg_66[k];

        t_67[k] = -5.0 * ig_67[k]
                  + f_0 * lg_67[k];

        t_68[k] = -5.0 * ig_68[k]
                  + f_0 * lg_68[k];

        t_69[k] = -5.0 * ig_69[k]
                  + f_0 * lg_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ig_70, ig_71, ig_72, ig_73, ig_74, \
                         lg_70, lg_71, lg_72, lg_73, lg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -5.0 * ig_70[k]
                  + f_0 * lg_70[k];

        t_71[k] = -5.0 * ig_71[k]
                  + f_0 * lg_71[k];

        t_72[k] = -5.0 * ig_72[k]
                  + f_0 * lg_72[k];

        t_73[k] = -5.0 * ig_73[k]
                  + f_0 * lg_73[k];

        t_74[k] = -5.0 * ig_74[k]
                  + f_0 * lg_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ig_75, ig_76, ig_77, ig_78, ig_79, \
                         lg_75, lg_76, lg_77, lg_78, lg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -5.0 * ig_75[k]
                  + f_0 * lg_75[k];

        t_76[k] = -5.0 * ig_76[k]
                  + f_0 * lg_76[k];

        t_77[k] = -5.0 * ig_77[k]
                  + f_0 * lg_77[k];

        t_78[k] = -5.0 * ig_78[k]
                  + f_0 * lg_78[k];

        t_79[k] = -5.0 * ig_79[k]
                  + f_0 * lg_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ig_80, ig_81, ig_82, ig_83, ig_84, \
                         lg_80, lg_81, lg_82, lg_83, lg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -5.0 * ig_80[k]
                  + f_0 * lg_80[k];

        t_81[k] = -5.0 * ig_81[k]
                  + f_0 * lg_81[k];

        t_82[k] = -5.0 * ig_82[k]
                  + f_0 * lg_82[k];

        t_83[k] = -5.0 * ig_83[k]
                  + f_0 * lg_83[k];

        t_84[k] = -5.0 * ig_84[k]
                  + f_0 * lg_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ig_85, ig_86, ig_87, ig_88, ig_89, \
                         lg_85, lg_86, lg_87, lg_88, lg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -5.0 * ig_85[k]
                  + f_0 * lg_85[k];

        t_86[k] = -5.0 * ig_86[k]
                  + f_0 * lg_86[k];

        t_87[k] = -5.0 * ig_87[k]
                  + f_0 * lg_87[k];

        t_88[k] = -5.0 * ig_88[k]
                  + f_0 * lg_88[k];

        t_89[k] = -5.0 * ig_89[k]
                  + f_0 * lg_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ig_90, ig_91, ig_92, ig_93, ig_94, \
                         lg_90, lg_91, lg_92, lg_93, lg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -4.0 * ig_90[k]
                  + f_0 * lg_90[k];

        t_91[k] = -4.0 * ig_91[k]
                  + f_0 * lg_91[k];

        t_92[k] = -4.0 * ig_92[k]
                  + f_0 * lg_92[k];

        t_93[k] = -4.0 * ig_93[k]
                  + f_0 * lg_93[k];

        t_94[k] = -4.0 * ig_94[k]
                  + f_0 * lg_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ig_95, ig_96, ig_97, ig_98, ig_99, \
                         lg_95, lg_96, lg_97, lg_98, lg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -4.0 * ig_95[k]
                  + f_0 * lg_95[k];

        t_96[k] = -4.0 * ig_96[k]
                  + f_0 * lg_96[k];

        t_97[k] = -4.0 * ig_97[k]
                  + f_0 * lg_97[k];

        t_98[k] = -4.0 * ig_98[k]
                  + f_0 * lg_98[k];

        t_99[k] = -4.0 * ig_99[k]
                  + f_0 * lg_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ig_100, ig_101, ig_102, ig_103, \
                         ig_104, lg_100, lg_101, lg_102, lg_103, \
                         lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -4.0 * ig_100[k]
                   + f_0 * lg_100[k];

        t_101[k] = -4.0 * ig_101[k]
                   + f_0 * lg_101[k];

        t_102[k] = -4.0 * ig_102[k]
                   + f_0 * lg_102[k];

        t_103[k] = -4.0 * ig_103[k]
                   + f_0 * lg_103[k];

        t_104[k] = -4.0 * ig_104[k]
                   + f_0 * lg_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ig_105, ig_106, ig_107, ig_108, \
                         ig_109, lg_105, lg_106, lg_107, lg_108, \
                         lg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -4.0 * ig_105[k]
                   + f_0 * lg_105[k];

        t_106[k] = -4.0 * ig_106[k]
                   + f_0 * lg_106[k];

        t_107[k] = -4.0 * ig_107[k]
                   + f_0 * lg_107[k];

        t_108[k] = -4.0 * ig_108[k]
                   + f_0 * lg_108[k];

        t_109[k] = -4.0 * ig_109[k]
                   + f_0 * lg_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ig_110, ig_111, ig_112, ig_113, \
                         ig_114, lg_110, lg_111, lg_112, lg_113, \
                         lg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -4.0 * ig_110[k]
                   + f_0 * lg_110[k];

        t_111[k] = -4.0 * ig_111[k]
                   + f_0 * lg_111[k];

        t_112[k] = -4.0 * ig_112[k]
                   + f_0 * lg_112[k];

        t_113[k] = -4.0 * ig_113[k]
                   + f_0 * lg_113[k];

        t_114[k] = -4.0 * ig_114[k]
                   + f_0 * lg_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ig_115, ig_116, ig_117, ig_118, \
                         ig_119, lg_115, lg_116, lg_117, lg_118, \
                         lg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -4.0 * ig_115[k]
                   + f_0 * lg_115[k];

        t_116[k] = -4.0 * ig_116[k]
                   + f_0 * lg_116[k];

        t_117[k] = -4.0 * ig_117[k]
                   + f_0 * lg_117[k];

        t_118[k] = -4.0 * ig_118[k]
                   + f_0 * lg_118[k];

        t_119[k] = -4.0 * ig_119[k]
                   + f_0 * lg_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ig_120, ig_121, ig_122, ig_123, \
                         ig_124, lg_120, lg_121, lg_122, lg_123, \
                         lg_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -4.0 * ig_120[k]
                   + f_0 * lg_120[k];

        t_121[k] = -4.0 * ig_121[k]
                   + f_0 * lg_121[k];

        t_122[k] = -4.0 * ig_122[k]
                   + f_0 * lg_122[k];

        t_123[k] = -4.0 * ig_123[k]
                   + f_0 * lg_123[k];

        t_124[k] = -4.0 * ig_124[k]
                   + f_0 * lg_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ig_125, ig_126, ig_127, ig_128, \
                         ig_129, lg_125, lg_126, lg_127, lg_128, \
                         lg_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -4.0 * ig_125[k]
                   + f_0 * lg_125[k];

        t_126[k] = -4.0 * ig_126[k]
                   + f_0 * lg_126[k];

        t_127[k] = -4.0 * ig_127[k]
                   + f_0 * lg_127[k];

        t_128[k] = -4.0 * ig_128[k]
                   + f_0 * lg_128[k];

        t_129[k] = -4.0 * ig_129[k]
                   + f_0 * lg_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ig_130, ig_131, ig_132, ig_133, \
                         ig_134, lg_130, lg_131, lg_132, lg_133, \
                         lg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -4.0 * ig_130[k]
                   + f_0 * lg_130[k];

        t_131[k] = -4.0 * ig_131[k]
                   + f_0 * lg_131[k];

        t_132[k] = -4.0 * ig_132[k]
                   + f_0 * lg_132[k];

        t_133[k] = -4.0 * ig_133[k]
                   + f_0 * lg_133[k];

        t_134[k] = -4.0 * ig_134[k]
                   + f_0 * lg_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ig_135, ig_136, ig_137, ig_138, \
                         ig_139, lg_135, lg_136, lg_137, lg_138, \
                         lg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -4.0 * ig_135[k]
                   + f_0 * lg_135[k];

        t_136[k] = -4.0 * ig_136[k]
                   + f_0 * lg_136[k];

        t_137[k] = -4.0 * ig_137[k]
                   + f_0 * lg_137[k];

        t_138[k] = -4.0 * ig_138[k]
                   + f_0 * lg_138[k];

        t_139[k] = -4.0 * ig_139[k]
                   + f_0 * lg_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ig_140, ig_141, ig_142, ig_143, \
                         ig_144, lg_140, lg_141, lg_142, lg_143, \
                         lg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -4.0 * ig_140[k]
                   + f_0 * lg_140[k];

        t_141[k] = -4.0 * ig_141[k]
                   + f_0 * lg_141[k];

        t_142[k] = -4.0 * ig_142[k]
                   + f_0 * lg_142[k];

        t_143[k] = -4.0 * ig_143[k]
                   + f_0 * lg_143[k];

        t_144[k] = -4.0 * ig_144[k]
                   + f_0 * lg_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ig_145, ig_146, ig_147, ig_148, \
                         ig_149, lg_145, lg_146, lg_147, lg_148, \
                         lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -4.0 * ig_145[k]
                   + f_0 * lg_145[k];

        t_146[k] = -4.0 * ig_146[k]
                   + f_0 * lg_146[k];

        t_147[k] = -4.0 * ig_147[k]
                   + f_0 * lg_147[k];

        t_148[k] = -4.0 * ig_148[k]
                   + f_0 * lg_148[k];

        t_149[k] = -4.0 * ig_149[k]
                   + f_0 * lg_149[k];
    }
}

static auto
compute_prim_geom_10_kg_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ig, const size_t lg,
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

    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_154 = buffer.data(lg + 154);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_157 = buffer.data(lg + 157);
    const auto *lg_158 = buffer.data(lg + 158);
    const auto *lg_159 = buffer.data(lg + 159);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);
    const auto *lg_168 = buffer.data(lg + 168);
    const auto *lg_169 = buffer.data(lg + 169);
    const auto *lg_170 = buffer.data(lg + 170);
    const auto *lg_171 = buffer.data(lg + 171);
    const auto *lg_172 = buffer.data(lg + 172);
    const auto *lg_173 = buffer.data(lg + 173);
    const auto *lg_174 = buffer.data(lg + 174);
    const auto *lg_175 = buffer.data(lg + 175);
    const auto *lg_176 = buffer.data(lg + 176);
    const auto *lg_177 = buffer.data(lg + 177);
    const auto *lg_178 = buffer.data(lg + 178);
    const auto *lg_179 = buffer.data(lg + 179);
    const auto *lg_180 = buffer.data(lg + 180);
    const auto *lg_181 = buffer.data(lg + 181);
    const auto *lg_182 = buffer.data(lg + 182);
    const auto *lg_183 = buffer.data(lg + 183);
    const auto *lg_184 = buffer.data(lg + 184);
    const auto *lg_185 = buffer.data(lg + 185);
    const auto *lg_186 = buffer.data(lg + 186);
    const auto *lg_187 = buffer.data(lg + 187);
    const auto *lg_188 = buffer.data(lg + 188);
    const auto *lg_189 = buffer.data(lg + 189);
    const auto *lg_190 = buffer.data(lg + 190);
    const auto *lg_191 = buffer.data(lg + 191);
    const auto *lg_192 = buffer.data(lg + 192);
    const auto *lg_193 = buffer.data(lg + 193);
    const auto *lg_194 = buffer.data(lg + 194);
    const auto *lg_195 = buffer.data(lg + 195);
    const auto *lg_196 = buffer.data(lg + 196);
    const auto *lg_197 = buffer.data(lg + 197);
    const auto *lg_198 = buffer.data(lg + 198);
    const auto *lg_199 = buffer.data(lg + 199);
    const auto *lg_200 = buffer.data(lg + 200);
    const auto *lg_201 = buffer.data(lg + 201);
    const auto *lg_202 = buffer.data(lg + 202);
    const auto *lg_203 = buffer.data(lg + 203);
    const auto *lg_204 = buffer.data(lg + 204);
    const auto *lg_205 = buffer.data(lg + 205);
    const auto *lg_206 = buffer.data(lg + 206);
    const auto *lg_207 = buffer.data(lg + 207);
    const auto *lg_208 = buffer.data(lg + 208);
    const auto *lg_209 = buffer.data(lg + 209);
    const auto *lg_210 = buffer.data(lg + 210);
    const auto *lg_211 = buffer.data(lg + 211);
    const auto *lg_212 = buffer.data(lg + 212);
    const auto *lg_213 = buffer.data(lg + 213);
    const auto *lg_214 = buffer.data(lg + 214);
    const auto *lg_215 = buffer.data(lg + 215);
    const auto *lg_216 = buffer.data(lg + 216);
    const auto *lg_217 = buffer.data(lg + 217);
    const auto *lg_218 = buffer.data(lg + 218);
    const auto *lg_219 = buffer.data(lg + 219);
    const auto *lg_220 = buffer.data(lg + 220);
    const auto *lg_221 = buffer.data(lg + 221);
    const auto *lg_222 = buffer.data(lg + 222);
    const auto *lg_223 = buffer.data(lg + 223);
    const auto *lg_224 = buffer.data(lg + 224);
    const auto *lg_225 = buffer.data(lg + 225);
    const auto *lg_226 = buffer.data(lg + 226);
    const auto *lg_227 = buffer.data(lg + 227);
    const auto *lg_228 = buffer.data(lg + 228);
    const auto *lg_229 = buffer.data(lg + 229);
    const auto *lg_230 = buffer.data(lg + 230);
    const auto *lg_231 = buffer.data(lg + 231);
    const auto *lg_232 = buffer.data(lg + 232);
    const auto *lg_233 = buffer.data(lg + 233);
    const auto *lg_234 = buffer.data(lg + 234);
    const auto *lg_235 = buffer.data(lg + 235);
    const auto *lg_236 = buffer.data(lg + 236);
    const auto *lg_237 = buffer.data(lg + 237);
    const auto *lg_238 = buffer.data(lg + 238);
    const auto *lg_239 = buffer.data(lg + 239);
    const auto *lg_240 = buffer.data(lg + 240);
    const auto *lg_241 = buffer.data(lg + 241);
    const auto *lg_242 = buffer.data(lg + 242);
    const auto *lg_243 = buffer.data(lg + 243);
    const auto *lg_244 = buffer.data(lg + 244);
    const auto *lg_245 = buffer.data(lg + 245);
    const auto *lg_246 = buffer.data(lg + 246);
    const auto *lg_247 = buffer.data(lg + 247);
    const auto *lg_248 = buffer.data(lg + 248);
    const auto *lg_249 = buffer.data(lg + 249);
    const auto *lg_250 = buffer.data(lg + 250);
    const auto *lg_251 = buffer.data(lg + 251);
    const auto *lg_252 = buffer.data(lg + 252);
    const auto *lg_253 = buffer.data(lg + 253);
    const auto *lg_254 = buffer.data(lg + 254);
    const auto *lg_255 = buffer.data(lg + 255);
    const auto *lg_256 = buffer.data(lg + 256);
    const auto *lg_257 = buffer.data(lg + 257);
    const auto *lg_258 = buffer.data(lg + 258);
    const auto *lg_259 = buffer.data(lg + 259);
    const auto *lg_260 = buffer.data(lg + 260);
    const auto *lg_261 = buffer.data(lg + 261);
    const auto *lg_262 = buffer.data(lg + 262);
    const auto *lg_263 = buffer.data(lg + 263);
    const auto *lg_264 = buffer.data(lg + 264);
    const auto *lg_265 = buffer.data(lg + 265);
    const auto *lg_266 = buffer.data(lg + 266);
    const auto *lg_267 = buffer.data(lg + 267);
    const auto *lg_268 = buffer.data(lg + 268);
    const auto *lg_269 = buffer.data(lg + 269);
    const auto *lg_270 = buffer.data(lg + 270);
    const auto *lg_271 = buffer.data(lg + 271);
    const auto *lg_272 = buffer.data(lg + 272);
    const auto *lg_273 = buffer.data(lg + 273);
    const auto *lg_274 = buffer.data(lg + 274);
    const auto *lg_275 = buffer.data(lg + 275);
    const auto *lg_276 = buffer.data(lg + 276);
    const auto *lg_277 = buffer.data(lg + 277);
    const auto *lg_278 = buffer.data(lg + 278);
    const auto *lg_279 = buffer.data(lg + 279);
    const auto *lg_280 = buffer.data(lg + 280);
    const auto *lg_281 = buffer.data(lg + 281);
    const auto *lg_282 = buffer.data(lg + 282);
    const auto *lg_283 = buffer.data(lg + 283);
    const auto *lg_284 = buffer.data(lg + 284);
    const auto *lg_285 = buffer.data(lg + 285);
    const auto *lg_286 = buffer.data(lg + 286);
    const auto *lg_287 = buffer.data(lg + 287);
    const auto *lg_288 = buffer.data(lg + 288);
    const auto *lg_289 = buffer.data(lg + 289);
    const auto *lg_290 = buffer.data(lg + 290);
    const auto *lg_291 = buffer.data(lg + 291);
    const auto *lg_292 = buffer.data(lg + 292);
    const auto *lg_293 = buffer.data(lg + 293);
    const auto *lg_294 = buffer.data(lg + 294);
    const auto *lg_295 = buffer.data(lg + 295);
    const auto *lg_296 = buffer.data(lg + 296);
    const auto *lg_297 = buffer.data(lg + 297);
    const auto *lg_298 = buffer.data(lg + 298);
    const auto *lg_299 = buffer.data(lg + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ig_150, ig_151, ig_152, ig_153, \
                         ig_154, lg_150, lg_151, lg_152, lg_153, \
                         lg_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -3.0 * ig_150[k]
                   + f_0 * lg_150[k];

        t_151[k] = -3.0 * ig_151[k]
                   + f_0 * lg_151[k];

        t_152[k] = -3.0 * ig_152[k]
                   + f_0 * lg_152[k];

        t_153[k] = -3.0 * ig_153[k]
                   + f_0 * lg_153[k];

        t_154[k] = -3.0 * ig_154[k]
                   + f_0 * lg_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ig_155, ig_156, ig_157, ig_158, \
                         ig_159, lg_155, lg_156, lg_157, lg_158, \
                         lg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -3.0 * ig_155[k]
                   + f_0 * lg_155[k];

        t_156[k] = -3.0 * ig_156[k]
                   + f_0 * lg_156[k];

        t_157[k] = -3.0 * ig_157[k]
                   + f_0 * lg_157[k];

        t_158[k] = -3.0 * ig_158[k]
                   + f_0 * lg_158[k];

        t_159[k] = -3.0 * ig_159[k]
                   + f_0 * lg_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ig_160, ig_161, ig_162, ig_163, \
                         ig_164, lg_160, lg_161, lg_162, lg_163, \
                         lg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -3.0 * ig_160[k]
                   + f_0 * lg_160[k];

        t_161[k] = -3.0 * ig_161[k]
                   + f_0 * lg_161[k];

        t_162[k] = -3.0 * ig_162[k]
                   + f_0 * lg_162[k];

        t_163[k] = -3.0 * ig_163[k]
                   + f_0 * lg_163[k];

        t_164[k] = -3.0 * ig_164[k]
                   + f_0 * lg_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ig_165, ig_166, ig_167, ig_168, \
                         ig_169, lg_165, lg_166, lg_167, lg_168, \
                         lg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -3.0 * ig_165[k]
                   + f_0 * lg_165[k];

        t_166[k] = -3.0 * ig_166[k]
                   + f_0 * lg_166[k];

        t_167[k] = -3.0 * ig_167[k]
                   + f_0 * lg_167[k];

        t_168[k] = -3.0 * ig_168[k]
                   + f_0 * lg_168[k];

        t_169[k] = -3.0 * ig_169[k]
                   + f_0 * lg_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ig_170, ig_171, ig_172, ig_173, \
                         ig_174, lg_170, lg_171, lg_172, lg_173, \
                         lg_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * ig_170[k]
                   + f_0 * lg_170[k];

        t_171[k] = -3.0 * ig_171[k]
                   + f_0 * lg_171[k];

        t_172[k] = -3.0 * ig_172[k]
                   + f_0 * lg_172[k];

        t_173[k] = -3.0 * ig_173[k]
                   + f_0 * lg_173[k];

        t_174[k] = -3.0 * ig_174[k]
                   + f_0 * lg_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ig_175, ig_176, ig_177, ig_178, \
                         ig_179, lg_175, lg_176, lg_177, lg_178, \
                         lg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * ig_175[k]
                   + f_0 * lg_175[k];

        t_176[k] = -3.0 * ig_176[k]
                   + f_0 * lg_176[k];

        t_177[k] = -3.0 * ig_177[k]
                   + f_0 * lg_177[k];

        t_178[k] = -3.0 * ig_178[k]
                   + f_0 * lg_178[k];

        t_179[k] = -3.0 * ig_179[k]
                   + f_0 * lg_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ig_180, ig_181, ig_182, ig_183, \
                         ig_184, lg_180, lg_181, lg_182, lg_183, \
                         lg_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * ig_180[k]
                   + f_0 * lg_180[k];

        t_181[k] = -3.0 * ig_181[k]
                   + f_0 * lg_181[k];

        t_182[k] = -3.0 * ig_182[k]
                   + f_0 * lg_182[k];

        t_183[k] = -3.0 * ig_183[k]
                   + f_0 * lg_183[k];

        t_184[k] = -3.0 * ig_184[k]
                   + f_0 * lg_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ig_185, ig_186, ig_187, ig_188, \
                         ig_189, lg_185, lg_186, lg_187, lg_188, \
                         lg_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * ig_185[k]
                   + f_0 * lg_185[k];

        t_186[k] = -3.0 * ig_186[k]
                   + f_0 * lg_186[k];

        t_187[k] = -3.0 * ig_187[k]
                   + f_0 * lg_187[k];

        t_188[k] = -3.0 * ig_188[k]
                   + f_0 * lg_188[k];

        t_189[k] = -3.0 * ig_189[k]
                   + f_0 * lg_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ig_190, ig_191, ig_192, ig_193, \
                         ig_194, lg_190, lg_191, lg_192, lg_193, \
                         lg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * ig_190[k]
                   + f_0 * lg_190[k];

        t_191[k] = -3.0 * ig_191[k]
                   + f_0 * lg_191[k];

        t_192[k] = -3.0 * ig_192[k]
                   + f_0 * lg_192[k];

        t_193[k] = -3.0 * ig_193[k]
                   + f_0 * lg_193[k];

        t_194[k] = -3.0 * ig_194[k]
                   + f_0 * lg_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ig_195, ig_196, ig_197, ig_198, \
                         ig_199, lg_195, lg_196, lg_197, lg_198, \
                         lg_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * ig_195[k]
                   + f_0 * lg_195[k];

        t_196[k] = -3.0 * ig_196[k]
                   + f_0 * lg_196[k];

        t_197[k] = -3.0 * ig_197[k]
                   + f_0 * lg_197[k];

        t_198[k] = -3.0 * ig_198[k]
                   + f_0 * lg_198[k];

        t_199[k] = -3.0 * ig_199[k]
                   + f_0 * lg_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ig_200, ig_201, ig_202, ig_203, \
                         ig_204, lg_200, lg_201, lg_202, lg_203, \
                         lg_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -3.0 * ig_200[k]
                   + f_0 * lg_200[k];

        t_201[k] = -3.0 * ig_201[k]
                   + f_0 * lg_201[k];

        t_202[k] = -3.0 * ig_202[k]
                   + f_0 * lg_202[k];

        t_203[k] = -3.0 * ig_203[k]
                   + f_0 * lg_203[k];

        t_204[k] = -3.0 * ig_204[k]
                   + f_0 * lg_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ig_205, ig_206, ig_207, ig_208, \
                         ig_209, lg_205, lg_206, lg_207, lg_208, \
                         lg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -3.0 * ig_205[k]
                   + f_0 * lg_205[k];

        t_206[k] = -3.0 * ig_206[k]
                   + f_0 * lg_206[k];

        t_207[k] = -3.0 * ig_207[k]
                   + f_0 * lg_207[k];

        t_208[k] = -3.0 * ig_208[k]
                   + f_0 * lg_208[k];

        t_209[k] = -3.0 * ig_209[k]
                   + f_0 * lg_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ig_210, ig_211, ig_212, ig_213, \
                         ig_214, lg_210, lg_211, lg_212, lg_213, \
                         lg_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -3.0 * ig_210[k]
                   + f_0 * lg_210[k];

        t_211[k] = -3.0 * ig_211[k]
                   + f_0 * lg_211[k];

        t_212[k] = -3.0 * ig_212[k]
                   + f_0 * lg_212[k];

        t_213[k] = -3.0 * ig_213[k]
                   + f_0 * lg_213[k];

        t_214[k] = -3.0 * ig_214[k]
                   + f_0 * lg_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ig_215, ig_216, ig_217, ig_218, \
                         ig_219, lg_215, lg_216, lg_217, lg_218, \
                         lg_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -3.0 * ig_215[k]
                   + f_0 * lg_215[k];

        t_216[k] = -3.0 * ig_216[k]
                   + f_0 * lg_216[k];

        t_217[k] = -3.0 * ig_217[k]
                   + f_0 * lg_217[k];

        t_218[k] = -3.0 * ig_218[k]
                   + f_0 * lg_218[k];

        t_219[k] = -3.0 * ig_219[k]
                   + f_0 * lg_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ig_220, ig_221, ig_222, ig_223, \
                         ig_224, lg_220, lg_221, lg_222, lg_223, \
                         lg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -3.0 * ig_220[k]
                   + f_0 * lg_220[k];

        t_221[k] = -3.0 * ig_221[k]
                   + f_0 * lg_221[k];

        t_222[k] = -3.0 * ig_222[k]
                   + f_0 * lg_222[k];

        t_223[k] = -3.0 * ig_223[k]
                   + f_0 * lg_223[k];

        t_224[k] = -3.0 * ig_224[k]
                   + f_0 * lg_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ig_225, ig_226, ig_227, ig_228, \
                         ig_229, lg_225, lg_226, lg_227, lg_228, \
                         lg_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -2.0 * ig_225[k]
                   + f_0 * lg_225[k];

        t_226[k] = -2.0 * ig_226[k]
                   + f_0 * lg_226[k];

        t_227[k] = -2.0 * ig_227[k]
                   + f_0 * lg_227[k];

        t_228[k] = -2.0 * ig_228[k]
                   + f_0 * lg_228[k];

        t_229[k] = -2.0 * ig_229[k]
                   + f_0 * lg_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ig_230, ig_231, ig_232, ig_233, \
                         ig_234, lg_230, lg_231, lg_232, lg_233, \
                         lg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -2.0 * ig_230[k]
                   + f_0 * lg_230[k];

        t_231[k] = -2.0 * ig_231[k]
                   + f_0 * lg_231[k];

        t_232[k] = -2.0 * ig_232[k]
                   + f_0 * lg_232[k];

        t_233[k] = -2.0 * ig_233[k]
                   + f_0 * lg_233[k];

        t_234[k] = -2.0 * ig_234[k]
                   + f_0 * lg_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ig_235, ig_236, ig_237, ig_238, \
                         ig_239, lg_235, lg_236, lg_237, lg_238, \
                         lg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -2.0 * ig_235[k]
                   + f_0 * lg_235[k];

        t_236[k] = -2.0 * ig_236[k]
                   + f_0 * lg_236[k];

        t_237[k] = -2.0 * ig_237[k]
                   + f_0 * lg_237[k];

        t_238[k] = -2.0 * ig_238[k]
                   + f_0 * lg_238[k];

        t_239[k] = -2.0 * ig_239[k]
                   + f_0 * lg_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ig_240, ig_241, ig_242, ig_243, \
                         ig_244, lg_240, lg_241, lg_242, lg_243, \
                         lg_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -2.0 * ig_240[k]
                   + f_0 * lg_240[k];

        t_241[k] = -2.0 * ig_241[k]
                   + f_0 * lg_241[k];

        t_242[k] = -2.0 * ig_242[k]
                   + f_0 * lg_242[k];

        t_243[k] = -2.0 * ig_243[k]
                   + f_0 * lg_243[k];

        t_244[k] = -2.0 * ig_244[k]
                   + f_0 * lg_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ig_245, ig_246, ig_247, ig_248, \
                         ig_249, lg_245, lg_246, lg_247, lg_248, \
                         lg_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -2.0 * ig_245[k]
                   + f_0 * lg_245[k];

        t_246[k] = -2.0 * ig_246[k]
                   + f_0 * lg_246[k];

        t_247[k] = -2.0 * ig_247[k]
                   + f_0 * lg_247[k];

        t_248[k] = -2.0 * ig_248[k]
                   + f_0 * lg_248[k];

        t_249[k] = -2.0 * ig_249[k]
                   + f_0 * lg_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ig_250, ig_251, ig_252, ig_253, \
                         ig_254, lg_250, lg_251, lg_252, lg_253, \
                         lg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * ig_250[k]
                   + f_0 * lg_250[k];

        t_251[k] = -2.0 * ig_251[k]
                   + f_0 * lg_251[k];

        t_252[k] = -2.0 * ig_252[k]
                   + f_0 * lg_252[k];

        t_253[k] = -2.0 * ig_253[k]
                   + f_0 * lg_253[k];

        t_254[k] = -2.0 * ig_254[k]
                   + f_0 * lg_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ig_255, ig_256, ig_257, ig_258, \
                         ig_259, lg_255, lg_256, lg_257, lg_258, \
                         lg_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * ig_255[k]
                   + f_0 * lg_255[k];

        t_256[k] = -2.0 * ig_256[k]
                   + f_0 * lg_256[k];

        t_257[k] = -2.0 * ig_257[k]
                   + f_0 * lg_257[k];

        t_258[k] = -2.0 * ig_258[k]
                   + f_0 * lg_258[k];

        t_259[k] = -2.0 * ig_259[k]
                   + f_0 * lg_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ig_260, ig_261, ig_262, ig_263, \
                         ig_264, lg_260, lg_261, lg_262, lg_263, \
                         lg_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * ig_260[k]
                   + f_0 * lg_260[k];

        t_261[k] = -2.0 * ig_261[k]
                   + f_0 * lg_261[k];

        t_262[k] = -2.0 * ig_262[k]
                   + f_0 * lg_262[k];

        t_263[k] = -2.0 * ig_263[k]
                   + f_0 * lg_263[k];

        t_264[k] = -2.0 * ig_264[k]
                   + f_0 * lg_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ig_265, ig_266, ig_267, ig_268, \
                         ig_269, lg_265, lg_266, lg_267, lg_268, \
                         lg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * ig_265[k]
                   + f_0 * lg_265[k];

        t_266[k] = -2.0 * ig_266[k]
                   + f_0 * lg_266[k];

        t_267[k] = -2.0 * ig_267[k]
                   + f_0 * lg_267[k];

        t_268[k] = -2.0 * ig_268[k]
                   + f_0 * lg_268[k];

        t_269[k] = -2.0 * ig_269[k]
                   + f_0 * lg_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ig_270, ig_271, ig_272, ig_273, \
                         ig_274, lg_270, lg_271, lg_272, lg_273, \
                         lg_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * ig_270[k]
                   + f_0 * lg_270[k];

        t_271[k] = -2.0 * ig_271[k]
                   + f_0 * lg_271[k];

        t_272[k] = -2.0 * ig_272[k]
                   + f_0 * lg_272[k];

        t_273[k] = -2.0 * ig_273[k]
                   + f_0 * lg_273[k];

        t_274[k] = -2.0 * ig_274[k]
                   + f_0 * lg_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ig_275, ig_276, ig_277, ig_278, \
                         ig_279, lg_275, lg_276, lg_277, lg_278, \
                         lg_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * ig_275[k]
                   + f_0 * lg_275[k];

        t_276[k] = -2.0 * ig_276[k]
                   + f_0 * lg_276[k];

        t_277[k] = -2.0 * ig_277[k]
                   + f_0 * lg_277[k];

        t_278[k] = -2.0 * ig_278[k]
                   + f_0 * lg_278[k];

        t_279[k] = -2.0 * ig_279[k]
                   + f_0 * lg_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ig_280, ig_281, ig_282, ig_283, \
                         ig_284, lg_280, lg_281, lg_282, lg_283, \
                         lg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * ig_280[k]
                   + f_0 * lg_280[k];

        t_281[k] = -2.0 * ig_281[k]
                   + f_0 * lg_281[k];

        t_282[k] = -2.0 * ig_282[k]
                   + f_0 * lg_282[k];

        t_283[k] = -2.0 * ig_283[k]
                   + f_0 * lg_283[k];

        t_284[k] = -2.0 * ig_284[k]
                   + f_0 * lg_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ig_285, ig_286, ig_287, ig_288, \
                         ig_289, lg_285, lg_286, lg_287, lg_288, \
                         lg_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -2.0 * ig_285[k]
                   + f_0 * lg_285[k];

        t_286[k] = -2.0 * ig_286[k]
                   + f_0 * lg_286[k];

        t_287[k] = -2.0 * ig_287[k]
                   + f_0 * lg_287[k];

        t_288[k] = -2.0 * ig_288[k]
                   + f_0 * lg_288[k];

        t_289[k] = -2.0 * ig_289[k]
                   + f_0 * lg_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ig_290, ig_291, ig_292, ig_293, \
                         ig_294, lg_290, lg_291, lg_292, lg_293, \
                         lg_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -2.0 * ig_290[k]
                   + f_0 * lg_290[k];

        t_291[k] = -2.0 * ig_291[k]
                   + f_0 * lg_291[k];

        t_292[k] = -2.0 * ig_292[k]
                   + f_0 * lg_292[k];

        t_293[k] = -2.0 * ig_293[k]
                   + f_0 * lg_293[k];

        t_294[k] = -2.0 * ig_294[k]
                   + f_0 * lg_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ig_295, ig_296, ig_297, ig_298, \
                         ig_299, lg_295, lg_296, lg_297, lg_298, \
                         lg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -2.0 * ig_295[k]
                   + f_0 * lg_295[k];

        t_296[k] = -2.0 * ig_296[k]
                   + f_0 * lg_296[k];

        t_297[k] = -2.0 * ig_297[k]
                   + f_0 * lg_297[k];

        t_298[k] = -2.0 * ig_298[k]
                   + f_0 * lg_298[k];

        t_299[k] = -2.0 * ig_299[k]
                   + f_0 * lg_299[k];
    }
}

static auto
compute_prim_geom_10_kg_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ig, const size_t lg,
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
    auto *t_455 = buffer.data(target + 455);
    auto *t_456 = buffer.data(target + 456);
    auto *t_457 = buffer.data(target + 457);
    auto *t_458 = buffer.data(target + 458);
    auto *t_459 = buffer.data(target + 459);

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

    const auto *lg_300 = buffer.data(lg + 300);
    const auto *lg_301 = buffer.data(lg + 301);
    const auto *lg_302 = buffer.data(lg + 302);
    const auto *lg_303 = buffer.data(lg + 303);
    const auto *lg_304 = buffer.data(lg + 304);
    const auto *lg_305 = buffer.data(lg + 305);
    const auto *lg_306 = buffer.data(lg + 306);
    const auto *lg_307 = buffer.data(lg + 307);
    const auto *lg_308 = buffer.data(lg + 308);
    const auto *lg_309 = buffer.data(lg + 309);
    const auto *lg_310 = buffer.data(lg + 310);
    const auto *lg_311 = buffer.data(lg + 311);
    const auto *lg_312 = buffer.data(lg + 312);
    const auto *lg_313 = buffer.data(lg + 313);
    const auto *lg_314 = buffer.data(lg + 314);
    const auto *lg_315 = buffer.data(lg + 315);
    const auto *lg_316 = buffer.data(lg + 316);
    const auto *lg_317 = buffer.data(lg + 317);
    const auto *lg_318 = buffer.data(lg + 318);
    const auto *lg_319 = buffer.data(lg + 319);
    const auto *lg_320 = buffer.data(lg + 320);
    const auto *lg_321 = buffer.data(lg + 321);
    const auto *lg_322 = buffer.data(lg + 322);
    const auto *lg_323 = buffer.data(lg + 323);
    const auto *lg_324 = buffer.data(lg + 324);
    const auto *lg_325 = buffer.data(lg + 325);
    const auto *lg_326 = buffer.data(lg + 326);
    const auto *lg_327 = buffer.data(lg + 327);
    const auto *lg_328 = buffer.data(lg + 328);
    const auto *lg_329 = buffer.data(lg + 329);
    const auto *lg_330 = buffer.data(lg + 330);
    const auto *lg_331 = buffer.data(lg + 331);
    const auto *lg_332 = buffer.data(lg + 332);
    const auto *lg_333 = buffer.data(lg + 333);
    const auto *lg_334 = buffer.data(lg + 334);
    const auto *lg_335 = buffer.data(lg + 335);
    const auto *lg_336 = buffer.data(lg + 336);
    const auto *lg_337 = buffer.data(lg + 337);
    const auto *lg_338 = buffer.data(lg + 338);
    const auto *lg_339 = buffer.data(lg + 339);
    const auto *lg_340 = buffer.data(lg + 340);
    const auto *lg_341 = buffer.data(lg + 341);
    const auto *lg_342 = buffer.data(lg + 342);
    const auto *lg_343 = buffer.data(lg + 343);
    const auto *lg_344 = buffer.data(lg + 344);
    const auto *lg_345 = buffer.data(lg + 345);
    const auto *lg_346 = buffer.data(lg + 346);
    const auto *lg_347 = buffer.data(lg + 347);
    const auto *lg_348 = buffer.data(lg + 348);
    const auto *lg_349 = buffer.data(lg + 349);
    const auto *lg_350 = buffer.data(lg + 350);
    const auto *lg_351 = buffer.data(lg + 351);
    const auto *lg_352 = buffer.data(lg + 352);
    const auto *lg_353 = buffer.data(lg + 353);
    const auto *lg_354 = buffer.data(lg + 354);
    const auto *lg_355 = buffer.data(lg + 355);
    const auto *lg_356 = buffer.data(lg + 356);
    const auto *lg_357 = buffer.data(lg + 357);
    const auto *lg_358 = buffer.data(lg + 358);
    const auto *lg_359 = buffer.data(lg + 359);
    const auto *lg_360 = buffer.data(lg + 360);
    const auto *lg_361 = buffer.data(lg + 361);
    const auto *lg_362 = buffer.data(lg + 362);
    const auto *lg_363 = buffer.data(lg + 363);
    const auto *lg_364 = buffer.data(lg + 364);
    const auto *lg_365 = buffer.data(lg + 365);
    const auto *lg_366 = buffer.data(lg + 366);
    const auto *lg_367 = buffer.data(lg + 367);
    const auto *lg_368 = buffer.data(lg + 368);
    const auto *lg_369 = buffer.data(lg + 369);
    const auto *lg_370 = buffer.data(lg + 370);
    const auto *lg_371 = buffer.data(lg + 371);
    const auto *lg_372 = buffer.data(lg + 372);
    const auto *lg_373 = buffer.data(lg + 373);
    const auto *lg_374 = buffer.data(lg + 374);
    const auto *lg_375 = buffer.data(lg + 375);
    const auto *lg_376 = buffer.data(lg + 376);
    const auto *lg_377 = buffer.data(lg + 377);
    const auto *lg_378 = buffer.data(lg + 378);
    const auto *lg_379 = buffer.data(lg + 379);
    const auto *lg_380 = buffer.data(lg + 380);
    const auto *lg_381 = buffer.data(lg + 381);
    const auto *lg_382 = buffer.data(lg + 382);
    const auto *lg_383 = buffer.data(lg + 383);
    const auto *lg_384 = buffer.data(lg + 384);
    const auto *lg_385 = buffer.data(lg + 385);
    const auto *lg_386 = buffer.data(lg + 386);
    const auto *lg_387 = buffer.data(lg + 387);
    const auto *lg_388 = buffer.data(lg + 388);
    const auto *lg_389 = buffer.data(lg + 389);
    const auto *lg_390 = buffer.data(lg + 390);
    const auto *lg_391 = buffer.data(lg + 391);
    const auto *lg_392 = buffer.data(lg + 392);
    const auto *lg_393 = buffer.data(lg + 393);
    const auto *lg_394 = buffer.data(lg + 394);
    const auto *lg_395 = buffer.data(lg + 395);
    const auto *lg_396 = buffer.data(lg + 396);
    const auto *lg_397 = buffer.data(lg + 397);
    const auto *lg_398 = buffer.data(lg + 398);
    const auto *lg_399 = buffer.data(lg + 399);
    const auto *lg_400 = buffer.data(lg + 400);
    const auto *lg_401 = buffer.data(lg + 401);
    const auto *lg_402 = buffer.data(lg + 402);
    const auto *lg_403 = buffer.data(lg + 403);
    const auto *lg_404 = buffer.data(lg + 404);
    const auto *lg_405 = buffer.data(lg + 405);
    const auto *lg_406 = buffer.data(lg + 406);
    const auto *lg_407 = buffer.data(lg + 407);
    const auto *lg_408 = buffer.data(lg + 408);
    const auto *lg_409 = buffer.data(lg + 409);
    const auto *lg_410 = buffer.data(lg + 410);
    const auto *lg_411 = buffer.data(lg + 411);
    const auto *lg_412 = buffer.data(lg + 412);
    const auto *lg_413 = buffer.data(lg + 413);
    const auto *lg_414 = buffer.data(lg + 414);
    const auto *lg_415 = buffer.data(lg + 415);
    const auto *lg_416 = buffer.data(lg + 416);
    const auto *lg_417 = buffer.data(lg + 417);
    const auto *lg_418 = buffer.data(lg + 418);
    const auto *lg_419 = buffer.data(lg + 419);
    const auto *lg_420 = buffer.data(lg + 420);
    const auto *lg_421 = buffer.data(lg + 421);
    const auto *lg_422 = buffer.data(lg + 422);
    const auto *lg_423 = buffer.data(lg + 423);
    const auto *lg_424 = buffer.data(lg + 424);
    const auto *lg_425 = buffer.data(lg + 425);
    const auto *lg_426 = buffer.data(lg + 426);
    const auto *lg_427 = buffer.data(lg + 427);
    const auto *lg_428 = buffer.data(lg + 428);
    const auto *lg_429 = buffer.data(lg + 429);
    const auto *lg_430 = buffer.data(lg + 430);
    const auto *lg_431 = buffer.data(lg + 431);
    const auto *lg_432 = buffer.data(lg + 432);
    const auto *lg_433 = buffer.data(lg + 433);
    const auto *lg_434 = buffer.data(lg + 434);
    const auto *lg_435 = buffer.data(lg + 435);
    const auto *lg_436 = buffer.data(lg + 436);
    const auto *lg_437 = buffer.data(lg + 437);
    const auto *lg_438 = buffer.data(lg + 438);
    const auto *lg_439 = buffer.data(lg + 439);
    const auto *lg_440 = buffer.data(lg + 440);
    const auto *lg_441 = buffer.data(lg + 441);
    const auto *lg_442 = buffer.data(lg + 442);
    const auto *lg_443 = buffer.data(lg + 443);
    const auto *lg_444 = buffer.data(lg + 444);
    const auto *lg_445 = buffer.data(lg + 445);
    const auto *lg_446 = buffer.data(lg + 446);
    const auto *lg_447 = buffer.data(lg + 447);
    const auto *lg_448 = buffer.data(lg + 448);
    const auto *lg_449 = buffer.data(lg + 449);
    const auto *lg_450 = buffer.data(lg + 450);
    const auto *lg_451 = buffer.data(lg + 451);
    const auto *lg_452 = buffer.data(lg + 452);
    const auto *lg_453 = buffer.data(lg + 453);
    const auto *lg_454 = buffer.data(lg + 454);
    const auto *lg_455 = buffer.data(lg + 455);
    const auto *lg_456 = buffer.data(lg + 456);
    const auto *lg_457 = buffer.data(lg + 457);
    const auto *lg_458 = buffer.data(lg + 458);
    const auto *lg_459 = buffer.data(lg + 459);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ig_300, ig_301, ig_302, ig_303, \
                         ig_304, lg_300, lg_301, lg_302, lg_303, \
                         lg_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -2.0 * ig_300[k]
                   + f_0 * lg_300[k];

        t_301[k] = -2.0 * ig_301[k]
                   + f_0 * lg_301[k];

        t_302[k] = -2.0 * ig_302[k]
                   + f_0 * lg_302[k];

        t_303[k] = -2.0 * ig_303[k]
                   + f_0 * lg_303[k];

        t_304[k] = -2.0 * ig_304[k]
                   + f_0 * lg_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ig_305, ig_306, ig_307, ig_308, \
                         ig_309, lg_305, lg_306, lg_307, lg_308, \
                         lg_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -2.0 * ig_305[k]
                   + f_0 * lg_305[k];

        t_306[k] = -2.0 * ig_306[k]
                   + f_0 * lg_306[k];

        t_307[k] = -2.0 * ig_307[k]
                   + f_0 * lg_307[k];

        t_308[k] = -2.0 * ig_308[k]
                   + f_0 * lg_308[k];

        t_309[k] = -2.0 * ig_309[k]
                   + f_0 * lg_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ig_310, ig_311, ig_312, ig_313, \
                         ig_314, lg_310, lg_311, lg_312, lg_313, \
                         lg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -2.0 * ig_310[k]
                   + f_0 * lg_310[k];

        t_311[k] = -2.0 * ig_311[k]
                   + f_0 * lg_311[k];

        t_312[k] = -2.0 * ig_312[k]
                   + f_0 * lg_312[k];

        t_313[k] = -2.0 * ig_313[k]
                   + f_0 * lg_313[k];

        t_314[k] = -2.0 * ig_314[k]
                   + f_0 * lg_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ig_315, ig_316, ig_317, ig_318, \
                         ig_319, lg_315, lg_316, lg_317, lg_318, \
                         lg_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -ig_315[k]
                   + f_0 * lg_315[k];

        t_316[k] = -ig_316[k]
                   + f_0 * lg_316[k];

        t_317[k] = -ig_317[k]
                   + f_0 * lg_317[k];

        t_318[k] = -ig_318[k]
                   + f_0 * lg_318[k];

        t_319[k] = -ig_319[k]
                   + f_0 * lg_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ig_320, ig_321, ig_322, ig_323, \
                         ig_324, lg_320, lg_321, lg_322, lg_323, \
                         lg_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -ig_320[k]
                   + f_0 * lg_320[k];

        t_321[k] = -ig_321[k]
                   + f_0 * lg_321[k];

        t_322[k] = -ig_322[k]
                   + f_0 * lg_322[k];

        t_323[k] = -ig_323[k]
                   + f_0 * lg_323[k];

        t_324[k] = -ig_324[k]
                   + f_0 * lg_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ig_325, ig_326, ig_327, ig_328, \
                         ig_329, lg_325, lg_326, lg_327, lg_328, \
                         lg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -ig_325[k]
                   + f_0 * lg_325[k];

        t_326[k] = -ig_326[k]
                   + f_0 * lg_326[k];

        t_327[k] = -ig_327[k]
                   + f_0 * lg_327[k];

        t_328[k] = -ig_328[k]
                   + f_0 * lg_328[k];

        t_329[k] = -ig_329[k]
                   + f_0 * lg_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ig_330, ig_331, ig_332, ig_333, \
                         ig_334, lg_330, lg_331, lg_332, lg_333, \
                         lg_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -ig_330[k]
                   + f_0 * lg_330[k];

        t_331[k] = -ig_331[k]
                   + f_0 * lg_331[k];

        t_332[k] = -ig_332[k]
                   + f_0 * lg_332[k];

        t_333[k] = -ig_333[k]
                   + f_0 * lg_333[k];

        t_334[k] = -ig_334[k]
                   + f_0 * lg_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ig_335, ig_336, ig_337, ig_338, \
                         ig_339, lg_335, lg_336, lg_337, lg_338, \
                         lg_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -ig_335[k]
                   + f_0 * lg_335[k];

        t_336[k] = -ig_336[k]
                   + f_0 * lg_336[k];

        t_337[k] = -ig_337[k]
                   + f_0 * lg_337[k];

        t_338[k] = -ig_338[k]
                   + f_0 * lg_338[k];

        t_339[k] = -ig_339[k]
                   + f_0 * lg_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ig_340, ig_341, ig_342, ig_343, \
                         ig_344, lg_340, lg_341, lg_342, lg_343, \
                         lg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -ig_340[k]
                   + f_0 * lg_340[k];

        t_341[k] = -ig_341[k]
                   + f_0 * lg_341[k];

        t_342[k] = -ig_342[k]
                   + f_0 * lg_342[k];

        t_343[k] = -ig_343[k]
                   + f_0 * lg_343[k];

        t_344[k] = -ig_344[k]
                   + f_0 * lg_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ig_345, ig_346, ig_347, ig_348, \
                         ig_349, lg_345, lg_346, lg_347, lg_348, \
                         lg_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -ig_345[k]
                   + f_0 * lg_345[k];

        t_346[k] = -ig_346[k]
                   + f_0 * lg_346[k];

        t_347[k] = -ig_347[k]
                   + f_0 * lg_347[k];

        t_348[k] = -ig_348[k]
                   + f_0 * lg_348[k];

        t_349[k] = -ig_349[k]
                   + f_0 * lg_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ig_350, ig_351, ig_352, ig_353, \
                         ig_354, lg_350, lg_351, lg_352, lg_353, \
                         lg_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -ig_350[k]
                   + f_0 * lg_350[k];

        t_351[k] = -ig_351[k]
                   + f_0 * lg_351[k];

        t_352[k] = -ig_352[k]
                   + f_0 * lg_352[k];

        t_353[k] = -ig_353[k]
                   + f_0 * lg_353[k];

        t_354[k] = -ig_354[k]
                   + f_0 * lg_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ig_355, ig_356, ig_357, ig_358, \
                         ig_359, lg_355, lg_356, lg_357, lg_358, \
                         lg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -ig_355[k]
                   + f_0 * lg_355[k];

        t_356[k] = -ig_356[k]
                   + f_0 * lg_356[k];

        t_357[k] = -ig_357[k]
                   + f_0 * lg_357[k];

        t_358[k] = -ig_358[k]
                   + f_0 * lg_358[k];

        t_359[k] = -ig_359[k]
                   + f_0 * lg_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ig_360, ig_361, ig_362, ig_363, \
                         ig_364, lg_360, lg_361, lg_362, lg_363, \
                         lg_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -ig_360[k]
                   + f_0 * lg_360[k];

        t_361[k] = -ig_361[k]
                   + f_0 * lg_361[k];

        t_362[k] = -ig_362[k]
                   + f_0 * lg_362[k];

        t_363[k] = -ig_363[k]
                   + f_0 * lg_363[k];

        t_364[k] = -ig_364[k]
                   + f_0 * lg_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ig_365, ig_366, ig_367, ig_368, \
                         ig_369, lg_365, lg_366, lg_367, lg_368, \
                         lg_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -ig_365[k]
                   + f_0 * lg_365[k];

        t_366[k] = -ig_366[k]
                   + f_0 * lg_366[k];

        t_367[k] = -ig_367[k]
                   + f_0 * lg_367[k];

        t_368[k] = -ig_368[k]
                   + f_0 * lg_368[k];

        t_369[k] = -ig_369[k]
                   + f_0 * lg_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ig_370, ig_371, ig_372, ig_373, \
                         ig_374, lg_370, lg_371, lg_372, lg_373, \
                         lg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -ig_370[k]
                   + f_0 * lg_370[k];

        t_371[k] = -ig_371[k]
                   + f_0 * lg_371[k];

        t_372[k] = -ig_372[k]
                   + f_0 * lg_372[k];

        t_373[k] = -ig_373[k]
                   + f_0 * lg_373[k];

        t_374[k] = -ig_374[k]
                   + f_0 * lg_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ig_375, ig_376, ig_377, ig_378, \
                         ig_379, lg_375, lg_376, lg_377, lg_378, \
                         lg_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -ig_375[k]
                   + f_0 * lg_375[k];

        t_376[k] = -ig_376[k]
                   + f_0 * lg_376[k];

        t_377[k] = -ig_377[k]
                   + f_0 * lg_377[k];

        t_378[k] = -ig_378[k]
                   + f_0 * lg_378[k];

        t_379[k] = -ig_379[k]
                   + f_0 * lg_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ig_380, ig_381, ig_382, ig_383, \
                         ig_384, lg_380, lg_381, lg_382, lg_383, \
                         lg_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -ig_380[k]
                   + f_0 * lg_380[k];

        t_381[k] = -ig_381[k]
                   + f_0 * lg_381[k];

        t_382[k] = -ig_382[k]
                   + f_0 * lg_382[k];

        t_383[k] = -ig_383[k]
                   + f_0 * lg_383[k];

        t_384[k] = -ig_384[k]
                   + f_0 * lg_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ig_385, ig_386, ig_387, ig_388, \
                         ig_389, lg_385, lg_386, lg_387, lg_388, \
                         lg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -ig_385[k]
                   + f_0 * lg_385[k];

        t_386[k] = -ig_386[k]
                   + f_0 * lg_386[k];

        t_387[k] = -ig_387[k]
                   + f_0 * lg_387[k];

        t_388[k] = -ig_388[k]
                   + f_0 * lg_388[k];

        t_389[k] = -ig_389[k]
                   + f_0 * lg_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ig_390, ig_391, ig_392, ig_393, \
                         ig_394, lg_390, lg_391, lg_392, lg_393, \
                         lg_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -ig_390[k]
                   + f_0 * lg_390[k];

        t_391[k] = -ig_391[k]
                   + f_0 * lg_391[k];

        t_392[k] = -ig_392[k]
                   + f_0 * lg_392[k];

        t_393[k] = -ig_393[k]
                   + f_0 * lg_393[k];

        t_394[k] = -ig_394[k]
                   + f_0 * lg_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ig_395, ig_396, ig_397, ig_398, \
                         ig_399, lg_395, lg_396, lg_397, lg_398, \
                         lg_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -ig_395[k]
                   + f_0 * lg_395[k];

        t_396[k] = -ig_396[k]
                   + f_0 * lg_396[k];

        t_397[k] = -ig_397[k]
                   + f_0 * lg_397[k];

        t_398[k] = -ig_398[k]
                   + f_0 * lg_398[k];

        t_399[k] = -ig_399[k]
                   + f_0 * lg_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ig_400, ig_401, ig_402, ig_403, \
                         ig_404, lg_400, lg_401, lg_402, lg_403, \
                         lg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -ig_400[k]
                   + f_0 * lg_400[k];

        t_401[k] = -ig_401[k]
                   + f_0 * lg_401[k];

        t_402[k] = -ig_402[k]
                   + f_0 * lg_402[k];

        t_403[k] = -ig_403[k]
                   + f_0 * lg_403[k];

        t_404[k] = -ig_404[k]
                   + f_0 * lg_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ig_405, ig_406, ig_407, ig_408, \
                         ig_409, lg_405, lg_406, lg_407, lg_408, \
                         lg_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -ig_405[k]
                   + f_0 * lg_405[k];

        t_406[k] = -ig_406[k]
                   + f_0 * lg_406[k];

        t_407[k] = -ig_407[k]
                   + f_0 * lg_407[k];

        t_408[k] = -ig_408[k]
                   + f_0 * lg_408[k];

        t_409[k] = -ig_409[k]
                   + f_0 * lg_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ig_410, ig_411, ig_412, ig_413, \
                         ig_414, lg_410, lg_411, lg_412, lg_413, \
                         lg_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -ig_410[k]
                   + f_0 * lg_410[k];

        t_411[k] = -ig_411[k]
                   + f_0 * lg_411[k];

        t_412[k] = -ig_412[k]
                   + f_0 * lg_412[k];

        t_413[k] = -ig_413[k]
                   + f_0 * lg_413[k];

        t_414[k] = -ig_414[k]
                   + f_0 * lg_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ig_415, ig_416, ig_417, ig_418, \
                         ig_419, lg_415, lg_416, lg_417, lg_418, \
                         lg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -ig_415[k]
                   + f_0 * lg_415[k];

        t_416[k] = -ig_416[k]
                   + f_0 * lg_416[k];

        t_417[k] = -ig_417[k]
                   + f_0 * lg_417[k];

        t_418[k] = -ig_418[k]
                   + f_0 * lg_418[k];

        t_419[k] = -ig_419[k]
                   + f_0 * lg_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, t_425, t_426, t_427, lg_420, \
                         lg_421, lg_422, lg_423, lg_424, lg_425, lg_426, \
                         lg_427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_0 * lg_420[k];

        t_421[k] = f_0 * lg_421[k];

        t_422[k] = f_0 * lg_422[k];

        t_423[k] = f_0 * lg_423[k];

        t_424[k] = f_0 * lg_424[k];

        t_425[k] = f_0 * lg_425[k];

        t_426[k] = f_0 * lg_426[k];

        t_427[k] = f_0 * lg_427[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, t_433, t_434, t_435, lg_428, \
                         lg_429, lg_430, lg_431, lg_432, lg_433, lg_434, \
                         lg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_0 * lg_428[k];

        t_429[k] = f_0 * lg_429[k];

        t_430[k] = f_0 * lg_430[k];

        t_431[k] = f_0 * lg_431[k];

        t_432[k] = f_0 * lg_432[k];

        t_433[k] = f_0 * lg_433[k];

        t_434[k] = f_0 * lg_434[k];

        t_435[k] = f_0 * lg_435[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, t_441, t_442, t_443, lg_436, \
                         lg_437, lg_438, lg_439, lg_440, lg_441, lg_442, \
                         lg_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_0 * lg_436[k];

        t_437[k] = f_0 * lg_437[k];

        t_438[k] = f_0 * lg_438[k];

        t_439[k] = f_0 * lg_439[k];

        t_440[k] = f_0 * lg_440[k];

        t_441[k] = f_0 * lg_441[k];

        t_442[k] = f_0 * lg_442[k];

        t_443[k] = f_0 * lg_443[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, t_449, t_450, t_451, lg_444, \
                         lg_445, lg_446, lg_447, lg_448, lg_449, lg_450, \
                         lg_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_0 * lg_444[k];

        t_445[k] = f_0 * lg_445[k];

        t_446[k] = f_0 * lg_446[k];

        t_447[k] = f_0 * lg_447[k];

        t_448[k] = f_0 * lg_448[k];

        t_449[k] = f_0 * lg_449[k];

        t_450[k] = f_0 * lg_450[k];

        t_451[k] = f_0 * lg_451[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, t_457, t_458, t_459, lg_452, \
                         lg_453, lg_454, lg_455, lg_456, lg_457, lg_458, \
                         lg_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_0 * lg_452[k];

        t_453[k] = f_0 * lg_453[k];

        t_454[k] = f_0 * lg_454[k];

        t_455[k] = f_0 * lg_455[k];

        t_456[k] = f_0 * lg_456[k];

        t_457[k] = f_0 * lg_457[k];

        t_458[k] = f_0 * lg_458[k];

        t_459[k] = f_0 * lg_459[k];
    }
}

static auto
compute_prim_geom_10_kg_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t lg, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *lg_460 = buffer.data(lg + 460);
    const auto *lg_461 = buffer.data(lg + 461);
    const auto *lg_462 = buffer.data(lg + 462);
    const auto *lg_463 = buffer.data(lg + 463);
    const auto *lg_464 = buffer.data(lg + 464);
    const auto *lg_465 = buffer.data(lg + 465);
    const auto *lg_466 = buffer.data(lg + 466);
    const auto *lg_467 = buffer.data(lg + 467);
    const auto *lg_468 = buffer.data(lg + 468);
    const auto *lg_469 = buffer.data(lg + 469);
    const auto *lg_470 = buffer.data(lg + 470);
    const auto *lg_471 = buffer.data(lg + 471);
    const auto *lg_472 = buffer.data(lg + 472);
    const auto *lg_473 = buffer.data(lg + 473);
    const auto *lg_474 = buffer.data(lg + 474);
    const auto *lg_475 = buffer.data(lg + 475);
    const auto *lg_476 = buffer.data(lg + 476);
    const auto *lg_477 = buffer.data(lg + 477);
    const auto *lg_478 = buffer.data(lg + 478);
    const auto *lg_479 = buffer.data(lg + 479);
    const auto *lg_480 = buffer.data(lg + 480);
    const auto *lg_481 = buffer.data(lg + 481);
    const auto *lg_482 = buffer.data(lg + 482);
    const auto *lg_483 = buffer.data(lg + 483);
    const auto *lg_484 = buffer.data(lg + 484);
    const auto *lg_485 = buffer.data(lg + 485);
    const auto *lg_486 = buffer.data(lg + 486);
    const auto *lg_487 = buffer.data(lg + 487);
    const auto *lg_488 = buffer.data(lg + 488);
    const auto *lg_489 = buffer.data(lg + 489);
    const auto *lg_490 = buffer.data(lg + 490);
    const auto *lg_491 = buffer.data(lg + 491);
    const auto *lg_492 = buffer.data(lg + 492);
    const auto *lg_493 = buffer.data(lg + 493);
    const auto *lg_494 = buffer.data(lg + 494);
    const auto *lg_495 = buffer.data(lg + 495);
    const auto *lg_496 = buffer.data(lg + 496);
    const auto *lg_497 = buffer.data(lg + 497);
    const auto *lg_498 = buffer.data(lg + 498);
    const auto *lg_499 = buffer.data(lg + 499);
    const auto *lg_500 = buffer.data(lg + 500);
    const auto *lg_501 = buffer.data(lg + 501);
    const auto *lg_502 = buffer.data(lg + 502);
    const auto *lg_503 = buffer.data(lg + 503);
    const auto *lg_504 = buffer.data(lg + 504);
    const auto *lg_505 = buffer.data(lg + 505);
    const auto *lg_506 = buffer.data(lg + 506);
    const auto *lg_507 = buffer.data(lg + 507);
    const auto *lg_508 = buffer.data(lg + 508);
    const auto *lg_509 = buffer.data(lg + 509);
    const auto *lg_510 = buffer.data(lg + 510);
    const auto *lg_511 = buffer.data(lg + 511);
    const auto *lg_512 = buffer.data(lg + 512);
    const auto *lg_513 = buffer.data(lg + 513);
    const auto *lg_514 = buffer.data(lg + 514);
    const auto *lg_515 = buffer.data(lg + 515);
    const auto *lg_516 = buffer.data(lg + 516);
    const auto *lg_517 = buffer.data(lg + 517);
    const auto *lg_518 = buffer.data(lg + 518);
    const auto *lg_519 = buffer.data(lg + 519);
    const auto *lg_520 = buffer.data(lg + 520);
    const auto *lg_521 = buffer.data(lg + 521);
    const auto *lg_522 = buffer.data(lg + 522);
    const auto *lg_523 = buffer.data(lg + 523);
    const auto *lg_524 = buffer.data(lg + 524);
    const auto *lg_525 = buffer.data(lg + 525);
    const auto *lg_526 = buffer.data(lg + 526);
    const auto *lg_527 = buffer.data(lg + 527);
    const auto *lg_528 = buffer.data(lg + 528);
    const auto *lg_529 = buffer.data(lg + 529);
    const auto *lg_530 = buffer.data(lg + 530);
    const auto *lg_531 = buffer.data(lg + 531);
    const auto *lg_532 = buffer.data(lg + 532);
    const auto *lg_533 = buffer.data(lg + 533);
    const auto *lg_534 = buffer.data(lg + 534);
    const auto *lg_535 = buffer.data(lg + 535);
    const auto *lg_536 = buffer.data(lg + 536);
    const auto *lg_537 = buffer.data(lg + 537);
    const auto *lg_538 = buffer.data(lg + 538);
    const auto *lg_539 = buffer.data(lg + 539);

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, t_465, t_466, t_467, lg_460, \
                         lg_461, lg_462, lg_463, lg_464, lg_465, lg_466, \
                         lg_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_0 * lg_460[k];

        t_461[k] = f_0 * lg_461[k];

        t_462[k] = f_0 * lg_462[k];

        t_463[k] = f_0 * lg_463[k];

        t_464[k] = f_0 * lg_464[k];

        t_465[k] = f_0 * lg_465[k];

        t_466[k] = f_0 * lg_466[k];

        t_467[k] = f_0 * lg_467[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, t_471, t_472, t_473, t_474, t_475, lg_468, \
                         lg_469, lg_470, lg_471, lg_472, lg_473, lg_474, \
                         lg_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_0 * lg_468[k];

        t_469[k] = f_0 * lg_469[k];

        t_470[k] = f_0 * lg_470[k];

        t_471[k] = f_0 * lg_471[k];

        t_472[k] = f_0 * lg_472[k];

        t_473[k] = f_0 * lg_473[k];

        t_474[k] = f_0 * lg_474[k];

        t_475[k] = f_0 * lg_475[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, t_481, t_482, t_483, lg_476, \
                         lg_477, lg_478, lg_479, lg_480, lg_481, lg_482, \
                         lg_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_0 * lg_476[k];

        t_477[k] = f_0 * lg_477[k];

        t_478[k] = f_0 * lg_478[k];

        t_479[k] = f_0 * lg_479[k];

        t_480[k] = f_0 * lg_480[k];

        t_481[k] = f_0 * lg_481[k];

        t_482[k] = f_0 * lg_482[k];

        t_483[k] = f_0 * lg_483[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, t_488, t_489, t_490, t_491, lg_484, \
                         lg_485, lg_486, lg_487, lg_488, lg_489, lg_490, \
                         lg_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_0 * lg_484[k];

        t_485[k] = f_0 * lg_485[k];

        t_486[k] = f_0 * lg_486[k];

        t_487[k] = f_0 * lg_487[k];

        t_488[k] = f_0 * lg_488[k];

        t_489[k] = f_0 * lg_489[k];

        t_490[k] = f_0 * lg_490[k];

        t_491[k] = f_0 * lg_491[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, t_497, t_498, t_499, lg_492, \
                         lg_493, lg_494, lg_495, lg_496, lg_497, lg_498, \
                         lg_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_0 * lg_492[k];

        t_493[k] = f_0 * lg_493[k];

        t_494[k] = f_0 * lg_494[k];

        t_495[k] = f_0 * lg_495[k];

        t_496[k] = f_0 * lg_496[k];

        t_497[k] = f_0 * lg_497[k];

        t_498[k] = f_0 * lg_498[k];

        t_499[k] = f_0 * lg_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, t_505, t_506, t_507, lg_500, \
                         lg_501, lg_502, lg_503, lg_504, lg_505, lg_506, \
                         lg_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_0 * lg_500[k];

        t_501[k] = f_0 * lg_501[k];

        t_502[k] = f_0 * lg_502[k];

        t_503[k] = f_0 * lg_503[k];

        t_504[k] = f_0 * lg_504[k];

        t_505[k] = f_0 * lg_505[k];

        t_506[k] = f_0 * lg_506[k];

        t_507[k] = f_0 * lg_507[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, t_512, t_513, t_514, t_515, lg_508, \
                         lg_509, lg_510, lg_511, lg_512, lg_513, lg_514, \
                         lg_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_0 * lg_508[k];

        t_509[k] = f_0 * lg_509[k];

        t_510[k] = f_0 * lg_510[k];

        t_511[k] = f_0 * lg_511[k];

        t_512[k] = f_0 * lg_512[k];

        t_513[k] = f_0 * lg_513[k];

        t_514[k] = f_0 * lg_514[k];

        t_515[k] = f_0 * lg_515[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, t_520, t_521, t_522, t_523, lg_516, \
                         lg_517, lg_518, lg_519, lg_520, lg_521, lg_522, \
                         lg_523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_0 * lg_516[k];

        t_517[k] = f_0 * lg_517[k];

        t_518[k] = f_0 * lg_518[k];

        t_519[k] = f_0 * lg_519[k];

        t_520[k] = f_0 * lg_520[k];

        t_521[k] = f_0 * lg_521[k];

        t_522[k] = f_0 * lg_522[k];

        t_523[k] = f_0 * lg_523[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, t_529, t_530, t_531, lg_524, \
                         lg_525, lg_526, lg_527, lg_528, lg_529, lg_530, \
                         lg_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = f_0 * lg_524[k];

        t_525[k] = f_0 * lg_525[k];

        t_526[k] = f_0 * lg_526[k];

        t_527[k] = f_0 * lg_527[k];

        t_528[k] = f_0 * lg_528[k];

        t_529[k] = f_0 * lg_529[k];

        t_530[k] = f_0 * lg_530[k];

        t_531[k] = f_0 * lg_531[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, t_537, t_538, t_539, lg_532, \
                         lg_533, lg_534, lg_535, lg_536, lg_537, lg_538, \
                         lg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_0 * lg_532[k];

        t_533[k] = f_0 * lg_533[k];

        t_534[k] = f_0 * lg_534[k];

        t_535[k] = f_0 * lg_535[k];

        t_536[k] = f_0 * lg_536[k];

        t_537[k] = f_0 * lg_537[k];

        t_538[k] = f_0 * lg_538[k];

        t_539[k] = f_0 * lg_539[k];
    }
}

auto
compute_prim_geom_10_kg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t ig, const size_t lg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kg_electron_repulsion_0_piece0(buffer, target, ig, lg, ncols, alpha);

    compute_prim_geom_10_kg_electron_repulsion_0_piece1(buffer, target, ig, lg, ncols, alpha);

    compute_prim_geom_10_kg_electron_repulsion_0_piece2(buffer, target, ig, lg, ncols, alpha);

    compute_prim_geom_10_kg_electron_repulsion_0_piece3(buffer, target, lg, ncols, alpha);
}

static auto
compute_prim_geom_10_kg_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ig, const size_t lg,
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

    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_17 = buffer.data(lg + 17);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_19 = buffer.data(lg + 19);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_22 = buffer.data(lg + 22);
    const auto *lg_23 = buffer.data(lg + 23);
    const auto *lg_24 = buffer.data(lg + 24);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_26 = buffer.data(lg + 26);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_49 = buffer.data(lg + 49);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_52 = buffer.data(lg + 52);
    const auto *lg_53 = buffer.data(lg + 53);
    const auto *lg_54 = buffer.data(lg + 54);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_94 = buffer.data(lg + 94);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_97 = buffer.data(lg + 97);
    const auto *lg_98 = buffer.data(lg + 98);
    const auto *lg_99 = buffer.data(lg + 99);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_154 = buffer.data(lg + 154);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_157 = buffer.data(lg + 157);
    const auto *lg_158 = buffer.data(lg + 158);
    const auto *lg_159 = buffer.data(lg + 159);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);
    const auto *lg_168 = buffer.data(lg + 168);
    const auto *lg_169 = buffer.data(lg + 169);
    const auto *lg_170 = buffer.data(lg + 170);
    const auto *lg_171 = buffer.data(lg + 171);
    const auto *lg_172 = buffer.data(lg + 172);
    const auto *lg_173 = buffer.data(lg + 173);
    const auto *lg_174 = buffer.data(lg + 174);
    const auto *lg_175 = buffer.data(lg + 175);
    const auto *lg_176 = buffer.data(lg + 176);
    const auto *lg_177 = buffer.data(lg + 177);
    const auto *lg_178 = buffer.data(lg + 178);
    const auto *lg_179 = buffer.data(lg + 179);
    const auto *lg_180 = buffer.data(lg + 180);
    const auto *lg_181 = buffer.data(lg + 181);
    const auto *lg_182 = buffer.data(lg + 182);
    const auto *lg_183 = buffer.data(lg + 183);
    const auto *lg_184 = buffer.data(lg + 184);
    const auto *lg_185 = buffer.data(lg + 185);
    const auto *lg_186 = buffer.data(lg + 186);
    const auto *lg_187 = buffer.data(lg + 187);
    const auto *lg_188 = buffer.data(lg + 188);
    const auto *lg_189 = buffer.data(lg + 189);
    const auto *lg_190 = buffer.data(lg + 190);
    const auto *lg_191 = buffer.data(lg + 191);
    const auto *lg_192 = buffer.data(lg + 192);
    const auto *lg_193 = buffer.data(lg + 193);
    const auto *lg_194 = buffer.data(lg + 194);
    const auto *lg_195 = buffer.data(lg + 195);
    const auto *lg_196 = buffer.data(lg + 196);
    const auto *lg_197 = buffer.data(lg + 197);
    const auto *lg_198 = buffer.data(lg + 198);
    const auto *lg_199 = buffer.data(lg + 199);
    const auto *lg_200 = buffer.data(lg + 200);
    const auto *lg_201 = buffer.data(lg + 201);
    const auto *lg_202 = buffer.data(lg + 202);
    const auto *lg_203 = buffer.data(lg + 203);
    const auto *lg_204 = buffer.data(lg + 204);
    const auto *lg_205 = buffer.data(lg + 205);
    const auto *lg_206 = buffer.data(lg + 206);
    const auto *lg_207 = buffer.data(lg + 207);
    const auto *lg_208 = buffer.data(lg + 208);
    const auto *lg_209 = buffer.data(lg + 209);
    const auto *lg_225 = buffer.data(lg + 225);
    const auto *lg_226 = buffer.data(lg + 226);
    const auto *lg_227 = buffer.data(lg + 227);
    const auto *lg_228 = buffer.data(lg + 228);
    const auto *lg_229 = buffer.data(lg + 229);
    const auto *lg_230 = buffer.data(lg + 230);
    const auto *lg_231 = buffer.data(lg + 231);
    const auto *lg_232 = buffer.data(lg + 232);
    const auto *lg_233 = buffer.data(lg + 233);
    const auto *lg_234 = buffer.data(lg + 234);
    const auto *lg_235 = buffer.data(lg + 235);
    const auto *lg_236 = buffer.data(lg + 236);
    const auto *lg_237 = buffer.data(lg + 237);
    const auto *lg_238 = buffer.data(lg + 238);
    const auto *lg_239 = buffer.data(lg + 239);
    const auto *lg_240 = buffer.data(lg + 240);
    const auto *lg_241 = buffer.data(lg + 241);
    const auto *lg_242 = buffer.data(lg + 242);
    const auto *lg_243 = buffer.data(lg + 243);
    const auto *lg_244 = buffer.data(lg + 244);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, lg_15, lg_16, lg_17, lg_18, \
                         lg_19, lg_20, lg_21, lg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lg_15[k];

        t_1[k] = f_0 * lg_16[k];

        t_2[k] = f_0 * lg_17[k];

        t_3[k] = f_0 * lg_18[k];

        t_4[k] = f_0 * lg_19[k];

        t_5[k] = f_0 * lg_20[k];

        t_6[k] = f_0 * lg_21[k];

        t_7[k] = f_0 * lg_22[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, lg_23, lg_24, lg_25, lg_26, \
                         lg_27, lg_28, lg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * lg_23[k];

        t_9[k] = f_0 * lg_24[k];

        t_10[k] = f_0 * lg_25[k];

        t_11[k] = f_0 * lg_26[k];

        t_12[k] = f_0 * lg_27[k];

        t_13[k] = f_0 * lg_28[k];

        t_14[k] = f_0 * lg_29[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ig_0, ig_1, ig_2, ig_3, ig_4, lg_45, \
                         lg_46, lg_47, lg_48, lg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -ig_0[k]
                  + f_0 * lg_45[k];

        t_16[k] = -ig_1[k]
                  + f_0 * lg_46[k];

        t_17[k] = -ig_2[k]
                  + f_0 * lg_47[k];

        t_18[k] = -ig_3[k]
                  + f_0 * lg_48[k];

        t_19[k] = -ig_4[k]
                  + f_0 * lg_49[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ig_5, ig_6, ig_7, ig_8, ig_9, lg_50, \
                         lg_51, lg_52, lg_53, lg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -ig_5[k]
                  + f_0 * lg_50[k];

        t_21[k] = -ig_6[k]
                  + f_0 * lg_51[k];

        t_22[k] = -ig_7[k]
                  + f_0 * lg_52[k];

        t_23[k] = -ig_8[k]
                  + f_0 * lg_53[k];

        t_24[k] = -ig_9[k]
                  + f_0 * lg_54[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ig_10, ig_11, ig_12, ig_13, ig_14, \
                         lg_55, lg_56, lg_57, lg_58, lg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -ig_10[k]
                  + f_0 * lg_55[k];

        t_26[k] = -ig_11[k]
                  + f_0 * lg_56[k];

        t_27[k] = -ig_12[k]
                  + f_0 * lg_57[k];

        t_28[k] = -ig_13[k]
                  + f_0 * lg_58[k];

        t_29[k] = -ig_14[k]
                  + f_0 * lg_59[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, t_36, t_37, lg_60, lg_61, lg_62, \
                         lg_63, lg_64, lg_65, lg_66, lg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * lg_60[k];

        t_31[k] = f_0 * lg_61[k];

        t_32[k] = f_0 * lg_62[k];

        t_33[k] = f_0 * lg_63[k];

        t_34[k] = f_0 * lg_64[k];

        t_35[k] = f_0 * lg_65[k];

        t_36[k] = f_0 * lg_66[k];

        t_37[k] = f_0 * lg_67[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, lg_68, lg_69, lg_70, lg_71, \
                         lg_72, lg_73, lg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * lg_68[k];

        t_39[k] = f_0 * lg_69[k];

        t_40[k] = f_0 * lg_70[k];

        t_41[k] = f_0 * lg_71[k];

        t_42[k] = f_0 * lg_72[k];

        t_43[k] = f_0 * lg_73[k];

        t_44[k] = f_0 * lg_74[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ig_15, ig_16, ig_17, ig_18, ig_19, \
                         lg_90, lg_91, lg_92, lg_93, lg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * ig_15[k]
                  + f_0 * lg_90[k];

        t_46[k] = -2.0 * ig_16[k]
                  + f_0 * lg_91[k];

        t_47[k] = -2.0 * ig_17[k]
                  + f_0 * lg_92[k];

        t_48[k] = -2.0 * ig_18[k]
                  + f_0 * lg_93[k];

        t_49[k] = -2.0 * ig_19[k]
                  + f_0 * lg_94[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ig_20, ig_21, ig_22, ig_23, ig_24, \
                         lg_95, lg_96, lg_97, lg_98, lg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * ig_20[k]
                  + f_0 * lg_95[k];

        t_51[k] = -2.0 * ig_21[k]
                  + f_0 * lg_96[k];

        t_52[k] = -2.0 * ig_22[k]
                  + f_0 * lg_97[k];

        t_53[k] = -2.0 * ig_23[k]
                  + f_0 * lg_98[k];

        t_54[k] = -2.0 * ig_24[k]
                  + f_0 * lg_99[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ig_25, ig_26, ig_27, ig_28, ig_29, \
                         lg_100, lg_101, lg_102, lg_103, lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * ig_25[k]
                  + f_0 * lg_100[k];

        t_56[k] = -2.0 * ig_26[k]
                  + f_0 * lg_101[k];

        t_57[k] = -2.0 * ig_27[k]
                  + f_0 * lg_102[k];

        t_58[k] = -2.0 * ig_28[k]
                  + f_0 * lg_103[k];

        t_59[k] = -2.0 * ig_29[k]
                  + f_0 * lg_104[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ig_30, ig_31, ig_32, ig_33, ig_34, \
                         lg_105, lg_106, lg_107, lg_108, lg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -ig_30[k]
                  + f_0 * lg_105[k];

        t_61[k] = -ig_31[k]
                  + f_0 * lg_106[k];

        t_62[k] = -ig_32[k]
                  + f_0 * lg_107[k];

        t_63[k] = -ig_33[k]
                  + f_0 * lg_108[k];

        t_64[k] = -ig_34[k]
                  + f_0 * lg_109[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ig_35, ig_36, ig_37, ig_38, ig_39, \
                         lg_110, lg_111, lg_112, lg_113, lg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -ig_35[k]
                  + f_0 * lg_110[k];

        t_66[k] = -ig_36[k]
                  + f_0 * lg_111[k];

        t_67[k] = -ig_37[k]
                  + f_0 * lg_112[k];

        t_68[k] = -ig_38[k]
                  + f_0 * lg_113[k];

        t_69[k] = -ig_39[k]
                  + f_0 * lg_114[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ig_40, ig_41, ig_42, ig_43, ig_44, \
                         lg_115, lg_116, lg_117, lg_118, lg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -ig_40[k]
                  + f_0 * lg_115[k];

        t_71[k] = -ig_41[k]
                  + f_0 * lg_116[k];

        t_72[k] = -ig_42[k]
                  + f_0 * lg_117[k];

        t_73[k] = -ig_43[k]
                  + f_0 * lg_118[k];

        t_74[k] = -ig_44[k]
                  + f_0 * lg_119[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, t_80, t_81, t_82, lg_120, lg_121, \
                         lg_122, lg_123, lg_124, lg_125, lg_126, \
                         lg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * lg_120[k];

        t_76[k] = f_0 * lg_121[k];

        t_77[k] = f_0 * lg_122[k];

        t_78[k] = f_0 * lg_123[k];

        t_79[k] = f_0 * lg_124[k];

        t_80[k] = f_0 * lg_125[k];

        t_81[k] = f_0 * lg_126[k];

        t_82[k] = f_0 * lg_127[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, lg_128, lg_129, lg_130, \
                         lg_131, lg_132, lg_133, lg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * lg_128[k];

        t_84[k] = f_0 * lg_129[k];

        t_85[k] = f_0 * lg_130[k];

        t_86[k] = f_0 * lg_131[k];

        t_87[k] = f_0 * lg_132[k];

        t_88[k] = f_0 * lg_133[k];

        t_89[k] = f_0 * lg_134[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ig_45, ig_46, ig_47, ig_48, ig_49, \
                         lg_150, lg_151, lg_152, lg_153, lg_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * ig_45[k]
                  + f_0 * lg_150[k];

        t_91[k] = -3.0 * ig_46[k]
                  + f_0 * lg_151[k];

        t_92[k] = -3.0 * ig_47[k]
                  + f_0 * lg_152[k];

        t_93[k] = -3.0 * ig_48[k]
                  + f_0 * lg_153[k];

        t_94[k] = -3.0 * ig_49[k]
                  + f_0 * lg_154[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ig_50, ig_51, ig_52, ig_53, ig_54, \
                         lg_155, lg_156, lg_157, lg_158, lg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * ig_50[k]
                  + f_0 * lg_155[k];

        t_96[k] = -3.0 * ig_51[k]
                  + f_0 * lg_156[k];

        t_97[k] = -3.0 * ig_52[k]
                  + f_0 * lg_157[k];

        t_98[k] = -3.0 * ig_53[k]
                  + f_0 * lg_158[k];

        t_99[k] = -3.0 * ig_54[k]
                  + f_0 * lg_159[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ig_55, ig_56, ig_57, ig_58, ig_59, \
                         lg_160, lg_161, lg_162, lg_163, lg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * ig_55[k]
                   + f_0 * lg_160[k];

        t_101[k] = -3.0 * ig_56[k]
                   + f_0 * lg_161[k];

        t_102[k] = -3.0 * ig_57[k]
                   + f_0 * lg_162[k];

        t_103[k] = -3.0 * ig_58[k]
                   + f_0 * lg_163[k];

        t_104[k] = -3.0 * ig_59[k]
                   + f_0 * lg_164[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ig_60, ig_61, ig_62, ig_63, ig_64, \
                         lg_165, lg_166, lg_167, lg_168, lg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * ig_60[k]
                   + f_0 * lg_165[k];

        t_106[k] = -2.0 * ig_61[k]
                   + f_0 * lg_166[k];

        t_107[k] = -2.0 * ig_62[k]
                   + f_0 * lg_167[k];

        t_108[k] = -2.0 * ig_63[k]
                   + f_0 * lg_168[k];

        t_109[k] = -2.0 * ig_64[k]
                   + f_0 * lg_169[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ig_65, ig_66, ig_67, ig_68, ig_69, \
                         lg_170, lg_171, lg_172, lg_173, lg_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * ig_65[k]
                   + f_0 * lg_170[k];

        t_111[k] = -2.0 * ig_66[k]
                   + f_0 * lg_171[k];

        t_112[k] = -2.0 * ig_67[k]
                   + f_0 * lg_172[k];

        t_113[k] = -2.0 * ig_68[k]
                   + f_0 * lg_173[k];

        t_114[k] = -2.0 * ig_69[k]
                   + f_0 * lg_174[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ig_70, ig_71, ig_72, ig_73, ig_74, \
                         lg_175, lg_176, lg_177, lg_178, lg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * ig_70[k]
                   + f_0 * lg_175[k];

        t_116[k] = -2.0 * ig_71[k]
                   + f_0 * lg_176[k];

        t_117[k] = -2.0 * ig_72[k]
                   + f_0 * lg_177[k];

        t_118[k] = -2.0 * ig_73[k]
                   + f_0 * lg_178[k];

        t_119[k] = -2.0 * ig_74[k]
                   + f_0 * lg_179[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ig_75, ig_76, ig_77, ig_78, ig_79, \
                         lg_180, lg_181, lg_182, lg_183, lg_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -ig_75[k]
                   + f_0 * lg_180[k];

        t_121[k] = -ig_76[k]
                   + f_0 * lg_181[k];

        t_122[k] = -ig_77[k]
                   + f_0 * lg_182[k];

        t_123[k] = -ig_78[k]
                   + f_0 * lg_183[k];

        t_124[k] = -ig_79[k]
                   + f_0 * lg_184[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ig_80, ig_81, ig_82, ig_83, ig_84, \
                         lg_185, lg_186, lg_187, lg_188, lg_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -ig_80[k]
                   + f_0 * lg_185[k];

        t_126[k] = -ig_81[k]
                   + f_0 * lg_186[k];

        t_127[k] = -ig_82[k]
                   + f_0 * lg_187[k];

        t_128[k] = -ig_83[k]
                   + f_0 * lg_188[k];

        t_129[k] = -ig_84[k]
                   + f_0 * lg_189[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ig_85, ig_86, ig_87, ig_88, ig_89, \
                         lg_190, lg_191, lg_192, lg_193, lg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -ig_85[k]
                   + f_0 * lg_190[k];

        t_131[k] = -ig_86[k]
                   + f_0 * lg_191[k];

        t_132[k] = -ig_87[k]
                   + f_0 * lg_192[k];

        t_133[k] = -ig_88[k]
                   + f_0 * lg_193[k];

        t_134[k] = -ig_89[k]
                   + f_0 * lg_194[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, t_140, t_141, t_142, lg_195, \
                         lg_196, lg_197, lg_198, lg_199, lg_200, lg_201, \
                         lg_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * lg_195[k];

        t_136[k] = f_0 * lg_196[k];

        t_137[k] = f_0 * lg_197[k];

        t_138[k] = f_0 * lg_198[k];

        t_139[k] = f_0 * lg_199[k];

        t_140[k] = f_0 * lg_200[k];

        t_141[k] = f_0 * lg_201[k];

        t_142[k] = f_0 * lg_202[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, t_148, t_149, lg_203, lg_204, \
                         lg_205, lg_206, lg_207, lg_208, lg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_0 * lg_203[k];

        t_144[k] = f_0 * lg_204[k];

        t_145[k] = f_0 * lg_205[k];

        t_146[k] = f_0 * lg_206[k];

        t_147[k] = f_0 * lg_207[k];

        t_148[k] = f_0 * lg_208[k];

        t_149[k] = f_0 * lg_209[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ig_90, ig_91, ig_92, ig_93, ig_94, \
                         lg_225, lg_226, lg_227, lg_228, lg_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -4.0 * ig_90[k]
                   + f_0 * lg_225[k];

        t_151[k] = -4.0 * ig_91[k]
                   + f_0 * lg_226[k];

        t_152[k] = -4.0 * ig_92[k]
                   + f_0 * lg_227[k];

        t_153[k] = -4.0 * ig_93[k]
                   + f_0 * lg_228[k];

        t_154[k] = -4.0 * ig_94[k]
                   + f_0 * lg_229[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ig_95, ig_96, ig_97, ig_98, ig_99, \
                         lg_230, lg_231, lg_232, lg_233, lg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -4.0 * ig_95[k]
                   + f_0 * lg_230[k];

        t_156[k] = -4.0 * ig_96[k]
                   + f_0 * lg_231[k];

        t_157[k] = -4.0 * ig_97[k]
                   + f_0 * lg_232[k];

        t_158[k] = -4.0 * ig_98[k]
                   + f_0 * lg_233[k];

        t_159[k] = -4.0 * ig_99[k]
                   + f_0 * lg_234[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ig_100, ig_101, ig_102, ig_103, \
                         ig_104, lg_235, lg_236, lg_237, lg_238, \
                         lg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -4.0 * ig_100[k]
                   + f_0 * lg_235[k];

        t_161[k] = -4.0 * ig_101[k]
                   + f_0 * lg_236[k];

        t_162[k] = -4.0 * ig_102[k]
                   + f_0 * lg_237[k];

        t_163[k] = -4.0 * ig_103[k]
                   + f_0 * lg_238[k];

        t_164[k] = -4.0 * ig_104[k]
                   + f_0 * lg_239[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ig_105, ig_106, ig_107, ig_108, \
                         ig_109, lg_240, lg_241, lg_242, lg_243, \
                         lg_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -3.0 * ig_105[k]
                   + f_0 * lg_240[k];

        t_166[k] = -3.0 * ig_106[k]
                   + f_0 * lg_241[k];

        t_167[k] = -3.0 * ig_107[k]
                   + f_0 * lg_242[k];

        t_168[k] = -3.0 * ig_108[k]
                   + f_0 * lg_243[k];

        t_169[k] = -3.0 * ig_109[k]
                   + f_0 * lg_244[k];
    }
}

static auto
compute_prim_geom_10_kg_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ig, const size_t lg,
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

    const auto *lg_245 = buffer.data(lg + 245);
    const auto *lg_246 = buffer.data(lg + 246);
    const auto *lg_247 = buffer.data(lg + 247);
    const auto *lg_248 = buffer.data(lg + 248);
    const auto *lg_249 = buffer.data(lg + 249);
    const auto *lg_250 = buffer.data(lg + 250);
    const auto *lg_251 = buffer.data(lg + 251);
    const auto *lg_252 = buffer.data(lg + 252);
    const auto *lg_253 = buffer.data(lg + 253);
    const auto *lg_254 = buffer.data(lg + 254);
    const auto *lg_255 = buffer.data(lg + 255);
    const auto *lg_256 = buffer.data(lg + 256);
    const auto *lg_257 = buffer.data(lg + 257);
    const auto *lg_258 = buffer.data(lg + 258);
    const auto *lg_259 = buffer.data(lg + 259);
    const auto *lg_260 = buffer.data(lg + 260);
    const auto *lg_261 = buffer.data(lg + 261);
    const auto *lg_262 = buffer.data(lg + 262);
    const auto *lg_263 = buffer.data(lg + 263);
    const auto *lg_264 = buffer.data(lg + 264);
    const auto *lg_265 = buffer.data(lg + 265);
    const auto *lg_266 = buffer.data(lg + 266);
    const auto *lg_267 = buffer.data(lg + 267);
    const auto *lg_268 = buffer.data(lg + 268);
    const auto *lg_269 = buffer.data(lg + 269);
    const auto *lg_270 = buffer.data(lg + 270);
    const auto *lg_271 = buffer.data(lg + 271);
    const auto *lg_272 = buffer.data(lg + 272);
    const auto *lg_273 = buffer.data(lg + 273);
    const auto *lg_274 = buffer.data(lg + 274);
    const auto *lg_275 = buffer.data(lg + 275);
    const auto *lg_276 = buffer.data(lg + 276);
    const auto *lg_277 = buffer.data(lg + 277);
    const auto *lg_278 = buffer.data(lg + 278);
    const auto *lg_279 = buffer.data(lg + 279);
    const auto *lg_280 = buffer.data(lg + 280);
    const auto *lg_281 = buffer.data(lg + 281);
    const auto *lg_282 = buffer.data(lg + 282);
    const auto *lg_283 = buffer.data(lg + 283);
    const auto *lg_284 = buffer.data(lg + 284);
    const auto *lg_285 = buffer.data(lg + 285);
    const auto *lg_286 = buffer.data(lg + 286);
    const auto *lg_287 = buffer.data(lg + 287);
    const auto *lg_288 = buffer.data(lg + 288);
    const auto *lg_289 = buffer.data(lg + 289);
    const auto *lg_290 = buffer.data(lg + 290);
    const auto *lg_291 = buffer.data(lg + 291);
    const auto *lg_292 = buffer.data(lg + 292);
    const auto *lg_293 = buffer.data(lg + 293);
    const auto *lg_294 = buffer.data(lg + 294);
    const auto *lg_295 = buffer.data(lg + 295);
    const auto *lg_296 = buffer.data(lg + 296);
    const auto *lg_297 = buffer.data(lg + 297);
    const auto *lg_298 = buffer.data(lg + 298);
    const auto *lg_299 = buffer.data(lg + 299);
    const auto *lg_315 = buffer.data(lg + 315);
    const auto *lg_316 = buffer.data(lg + 316);
    const auto *lg_317 = buffer.data(lg + 317);
    const auto *lg_318 = buffer.data(lg + 318);
    const auto *lg_319 = buffer.data(lg + 319);
    const auto *lg_320 = buffer.data(lg + 320);
    const auto *lg_321 = buffer.data(lg + 321);
    const auto *lg_322 = buffer.data(lg + 322);
    const auto *lg_323 = buffer.data(lg + 323);
    const auto *lg_324 = buffer.data(lg + 324);
    const auto *lg_325 = buffer.data(lg + 325);
    const auto *lg_326 = buffer.data(lg + 326);
    const auto *lg_327 = buffer.data(lg + 327);
    const auto *lg_328 = buffer.data(lg + 328);
    const auto *lg_329 = buffer.data(lg + 329);
    const auto *lg_330 = buffer.data(lg + 330);
    const auto *lg_331 = buffer.data(lg + 331);
    const auto *lg_332 = buffer.data(lg + 332);
    const auto *lg_333 = buffer.data(lg + 333);
    const auto *lg_334 = buffer.data(lg + 334);
    const auto *lg_335 = buffer.data(lg + 335);
    const auto *lg_336 = buffer.data(lg + 336);
    const auto *lg_337 = buffer.data(lg + 337);
    const auto *lg_338 = buffer.data(lg + 338);
    const auto *lg_339 = buffer.data(lg + 339);
    const auto *lg_340 = buffer.data(lg + 340);
    const auto *lg_341 = buffer.data(lg + 341);
    const auto *lg_342 = buffer.data(lg + 342);
    const auto *lg_343 = buffer.data(lg + 343);
    const auto *lg_344 = buffer.data(lg + 344);
    const auto *lg_345 = buffer.data(lg + 345);
    const auto *lg_346 = buffer.data(lg + 346);
    const auto *lg_347 = buffer.data(lg + 347);
    const auto *lg_348 = buffer.data(lg + 348);
    const auto *lg_349 = buffer.data(lg + 349);
    const auto *lg_350 = buffer.data(lg + 350);
    const auto *lg_351 = buffer.data(lg + 351);
    const auto *lg_352 = buffer.data(lg + 352);
    const auto *lg_353 = buffer.data(lg + 353);
    const auto *lg_354 = buffer.data(lg + 354);
    const auto *lg_355 = buffer.data(lg + 355);
    const auto *lg_356 = buffer.data(lg + 356);
    const auto *lg_357 = buffer.data(lg + 357);
    const auto *lg_358 = buffer.data(lg + 358);
    const auto *lg_359 = buffer.data(lg + 359);
    const auto *lg_360 = buffer.data(lg + 360);
    const auto *lg_361 = buffer.data(lg + 361);
    const auto *lg_362 = buffer.data(lg + 362);
    const auto *lg_363 = buffer.data(lg + 363);
    const auto *lg_364 = buffer.data(lg + 364);
    const auto *lg_365 = buffer.data(lg + 365);
    const auto *lg_366 = buffer.data(lg + 366);
    const auto *lg_367 = buffer.data(lg + 367);
    const auto *lg_368 = buffer.data(lg + 368);
    const auto *lg_369 = buffer.data(lg + 369);
    const auto *lg_370 = buffer.data(lg + 370);
    const auto *lg_371 = buffer.data(lg + 371);
    const auto *lg_372 = buffer.data(lg + 372);
    const auto *lg_373 = buffer.data(lg + 373);
    const auto *lg_374 = buffer.data(lg + 374);
    const auto *lg_375 = buffer.data(lg + 375);
    const auto *lg_376 = buffer.data(lg + 376);
    const auto *lg_377 = buffer.data(lg + 377);
    const auto *lg_378 = buffer.data(lg + 378);
    const auto *lg_379 = buffer.data(lg + 379);
    const auto *lg_380 = buffer.data(lg + 380);
    const auto *lg_381 = buffer.data(lg + 381);
    const auto *lg_382 = buffer.data(lg + 382);
    const auto *lg_383 = buffer.data(lg + 383);
    const auto *lg_384 = buffer.data(lg + 384);
    const auto *lg_385 = buffer.data(lg + 385);
    const auto *lg_386 = buffer.data(lg + 386);
    const auto *lg_387 = buffer.data(lg + 387);
    const auto *lg_388 = buffer.data(lg + 388);
    const auto *lg_389 = buffer.data(lg + 389);
    const auto *lg_390 = buffer.data(lg + 390);
    const auto *lg_391 = buffer.data(lg + 391);
    const auto *lg_392 = buffer.data(lg + 392);
    const auto *lg_393 = buffer.data(lg + 393);
    const auto *lg_394 = buffer.data(lg + 394);
    const auto *lg_395 = buffer.data(lg + 395);
    const auto *lg_396 = buffer.data(lg + 396);
    const auto *lg_397 = buffer.data(lg + 397);
    const auto *lg_398 = buffer.data(lg + 398);
    const auto *lg_399 = buffer.data(lg + 399);
    const auto *lg_400 = buffer.data(lg + 400);
    const auto *lg_401 = buffer.data(lg + 401);
    const auto *lg_402 = buffer.data(lg + 402);
    const auto *lg_403 = buffer.data(lg + 403);
    const auto *lg_404 = buffer.data(lg + 404);
    const auto *lg_420 = buffer.data(lg + 420);
    const auto *lg_421 = buffer.data(lg + 421);
    const auto *lg_422 = buffer.data(lg + 422);
    const auto *lg_423 = buffer.data(lg + 423);
    const auto *lg_424 = buffer.data(lg + 424);
    const auto *lg_425 = buffer.data(lg + 425);
    const auto *lg_426 = buffer.data(lg + 426);
    const auto *lg_427 = buffer.data(lg + 427);
    const auto *lg_428 = buffer.data(lg + 428);
    const auto *lg_429 = buffer.data(lg + 429);
    const auto *lg_430 = buffer.data(lg + 430);
    const auto *lg_431 = buffer.data(lg + 431);
    const auto *lg_432 = buffer.data(lg + 432);
    const auto *lg_433 = buffer.data(lg + 433);
    const auto *lg_434 = buffer.data(lg + 434);

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ig_110, ig_111, ig_112, ig_113, \
                         ig_114, lg_245, lg_246, lg_247, lg_248, \
                         lg_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * ig_110[k]
                   + f_0 * lg_245[k];

        t_171[k] = -3.0 * ig_111[k]
                   + f_0 * lg_246[k];

        t_172[k] = -3.0 * ig_112[k]
                   + f_0 * lg_247[k];

        t_173[k] = -3.0 * ig_113[k]
                   + f_0 * lg_248[k];

        t_174[k] = -3.0 * ig_114[k]
                   + f_0 * lg_249[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ig_115, ig_116, ig_117, ig_118, \
                         ig_119, lg_250, lg_251, lg_252, lg_253, \
                         lg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * ig_115[k]
                   + f_0 * lg_250[k];

        t_176[k] = -3.0 * ig_116[k]
                   + f_0 * lg_251[k];

        t_177[k] = -3.0 * ig_117[k]
                   + f_0 * lg_252[k];

        t_178[k] = -3.0 * ig_118[k]
                   + f_0 * lg_253[k];

        t_179[k] = -3.0 * ig_119[k]
                   + f_0 * lg_254[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ig_120, ig_121, ig_122, ig_123, \
                         ig_124, lg_255, lg_256, lg_257, lg_258, \
                         lg_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * ig_120[k]
                   + f_0 * lg_255[k];

        t_181[k] = -2.0 * ig_121[k]
                   + f_0 * lg_256[k];

        t_182[k] = -2.0 * ig_122[k]
                   + f_0 * lg_257[k];

        t_183[k] = -2.0 * ig_123[k]
                   + f_0 * lg_258[k];

        t_184[k] = -2.0 * ig_124[k]
                   + f_0 * lg_259[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ig_125, ig_126, ig_127, ig_128, \
                         ig_129, lg_260, lg_261, lg_262, lg_263, \
                         lg_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * ig_125[k]
                   + f_0 * lg_260[k];

        t_186[k] = -2.0 * ig_126[k]
                   + f_0 * lg_261[k];

        t_187[k] = -2.0 * ig_127[k]
                   + f_0 * lg_262[k];

        t_188[k] = -2.0 * ig_128[k]
                   + f_0 * lg_263[k];

        t_189[k] = -2.0 * ig_129[k]
                   + f_0 * lg_264[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ig_130, ig_131, ig_132, ig_133, \
                         ig_134, lg_265, lg_266, lg_267, lg_268, \
                         lg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * ig_130[k]
                   + f_0 * lg_265[k];

        t_191[k] = -2.0 * ig_131[k]
                   + f_0 * lg_266[k];

        t_192[k] = -2.0 * ig_132[k]
                   + f_0 * lg_267[k];

        t_193[k] = -2.0 * ig_133[k]
                   + f_0 * lg_268[k];

        t_194[k] = -2.0 * ig_134[k]
                   + f_0 * lg_269[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ig_135, ig_136, ig_137, ig_138, \
                         ig_139, lg_270, lg_271, lg_272, lg_273, \
                         lg_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -ig_135[k]
                   + f_0 * lg_270[k];

        t_196[k] = -ig_136[k]
                   + f_0 * lg_271[k];

        t_197[k] = -ig_137[k]
                   + f_0 * lg_272[k];

        t_198[k] = -ig_138[k]
                   + f_0 * lg_273[k];

        t_199[k] = -ig_139[k]
                   + f_0 * lg_274[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ig_140, ig_141, ig_142, ig_143, \
                         ig_144, lg_275, lg_276, lg_277, lg_278, \
                         lg_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -ig_140[k]
                   + f_0 * lg_275[k];

        t_201[k] = -ig_141[k]
                   + f_0 * lg_276[k];

        t_202[k] = -ig_142[k]
                   + f_0 * lg_277[k];

        t_203[k] = -ig_143[k]
                   + f_0 * lg_278[k];

        t_204[k] = -ig_144[k]
                   + f_0 * lg_279[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ig_145, ig_146, ig_147, ig_148, \
                         ig_149, lg_280, lg_281, lg_282, lg_283, \
                         lg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -ig_145[k]
                   + f_0 * lg_280[k];

        t_206[k] = -ig_146[k]
                   + f_0 * lg_281[k];

        t_207[k] = -ig_147[k]
                   + f_0 * lg_282[k];

        t_208[k] = -ig_148[k]
                   + f_0 * lg_283[k];

        t_209[k] = -ig_149[k]
                   + f_0 * lg_284[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, t_217, lg_285, \
                         lg_286, lg_287, lg_288, lg_289, lg_290, lg_291, \
                         lg_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * lg_285[k];

        t_211[k] = f_0 * lg_286[k];

        t_212[k] = f_0 * lg_287[k];

        t_213[k] = f_0 * lg_288[k];

        t_214[k] = f_0 * lg_289[k];

        t_215[k] = f_0 * lg_290[k];

        t_216[k] = f_0 * lg_291[k];

        t_217[k] = f_0 * lg_292[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, t_224, lg_293, lg_294, \
                         lg_295, lg_296, lg_297, lg_298, lg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * lg_293[k];

        t_219[k] = f_0 * lg_294[k];

        t_220[k] = f_0 * lg_295[k];

        t_221[k] = f_0 * lg_296[k];

        t_222[k] = f_0 * lg_297[k];

        t_223[k] = f_0 * lg_298[k];

        t_224[k] = f_0 * lg_299[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ig_150, ig_151, ig_152, ig_153, \
                         ig_154, lg_315, lg_316, lg_317, lg_318, \
                         lg_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -5.0 * ig_150[k]
                   + f_0 * lg_315[k];

        t_226[k] = -5.0 * ig_151[k]
                   + f_0 * lg_316[k];

        t_227[k] = -5.0 * ig_152[k]
                   + f_0 * lg_317[k];

        t_228[k] = -5.0 * ig_153[k]
                   + f_0 * lg_318[k];

        t_229[k] = -5.0 * ig_154[k]
                   + f_0 * lg_319[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ig_155, ig_156, ig_157, ig_158, \
                         ig_159, lg_320, lg_321, lg_322, lg_323, \
                         lg_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -5.0 * ig_155[k]
                   + f_0 * lg_320[k];

        t_231[k] = -5.0 * ig_156[k]
                   + f_0 * lg_321[k];

        t_232[k] = -5.0 * ig_157[k]
                   + f_0 * lg_322[k];

        t_233[k] = -5.0 * ig_158[k]
                   + f_0 * lg_323[k];

        t_234[k] = -5.0 * ig_159[k]
                   + f_0 * lg_324[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ig_160, ig_161, ig_162, ig_163, \
                         ig_164, lg_325, lg_326, lg_327, lg_328, \
                         lg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -5.0 * ig_160[k]
                   + f_0 * lg_325[k];

        t_236[k] = -5.0 * ig_161[k]
                   + f_0 * lg_326[k];

        t_237[k] = -5.0 * ig_162[k]
                   + f_0 * lg_327[k];

        t_238[k] = -5.0 * ig_163[k]
                   + f_0 * lg_328[k];

        t_239[k] = -5.0 * ig_164[k]
                   + f_0 * lg_329[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ig_165, ig_166, ig_167, ig_168, \
                         ig_169, lg_330, lg_331, lg_332, lg_333, \
                         lg_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -4.0 * ig_165[k]
                   + f_0 * lg_330[k];

        t_241[k] = -4.0 * ig_166[k]
                   + f_0 * lg_331[k];

        t_242[k] = -4.0 * ig_167[k]
                   + f_0 * lg_332[k];

        t_243[k] = -4.0 * ig_168[k]
                   + f_0 * lg_333[k];

        t_244[k] = -4.0 * ig_169[k]
                   + f_0 * lg_334[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ig_170, ig_171, ig_172, ig_173, \
                         ig_174, lg_335, lg_336, lg_337, lg_338, \
                         lg_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -4.0 * ig_170[k]
                   + f_0 * lg_335[k];

        t_246[k] = -4.0 * ig_171[k]
                   + f_0 * lg_336[k];

        t_247[k] = -4.0 * ig_172[k]
                   + f_0 * lg_337[k];

        t_248[k] = -4.0 * ig_173[k]
                   + f_0 * lg_338[k];

        t_249[k] = -4.0 * ig_174[k]
                   + f_0 * lg_339[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ig_175, ig_176, ig_177, ig_178, \
                         ig_179, lg_340, lg_341, lg_342, lg_343, \
                         lg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -4.0 * ig_175[k]
                   + f_0 * lg_340[k];

        t_251[k] = -4.0 * ig_176[k]
                   + f_0 * lg_341[k];

        t_252[k] = -4.0 * ig_177[k]
                   + f_0 * lg_342[k];

        t_253[k] = -4.0 * ig_178[k]
                   + f_0 * lg_343[k];

        t_254[k] = -4.0 * ig_179[k]
                   + f_0 * lg_344[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ig_180, ig_181, ig_182, ig_183, \
                         ig_184, lg_345, lg_346, lg_347, lg_348, \
                         lg_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -3.0 * ig_180[k]
                   + f_0 * lg_345[k];

        t_256[k] = -3.0 * ig_181[k]
                   + f_0 * lg_346[k];

        t_257[k] = -3.0 * ig_182[k]
                   + f_0 * lg_347[k];

        t_258[k] = -3.0 * ig_183[k]
                   + f_0 * lg_348[k];

        t_259[k] = -3.0 * ig_184[k]
                   + f_0 * lg_349[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ig_185, ig_186, ig_187, ig_188, \
                         ig_189, lg_350, lg_351, lg_352, lg_353, \
                         lg_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -3.0 * ig_185[k]
                   + f_0 * lg_350[k];

        t_261[k] = -3.0 * ig_186[k]
                   + f_0 * lg_351[k];

        t_262[k] = -3.0 * ig_187[k]
                   + f_0 * lg_352[k];

        t_263[k] = -3.0 * ig_188[k]
                   + f_0 * lg_353[k];

        t_264[k] = -3.0 * ig_189[k]
                   + f_0 * lg_354[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ig_190, ig_191, ig_192, ig_193, \
                         ig_194, lg_355, lg_356, lg_357, lg_358, \
                         lg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -3.0 * ig_190[k]
                   + f_0 * lg_355[k];

        t_266[k] = -3.0 * ig_191[k]
                   + f_0 * lg_356[k];

        t_267[k] = -3.0 * ig_192[k]
                   + f_0 * lg_357[k];

        t_268[k] = -3.0 * ig_193[k]
                   + f_0 * lg_358[k];

        t_269[k] = -3.0 * ig_194[k]
                   + f_0 * lg_359[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ig_195, ig_196, ig_197, ig_198, \
                         ig_199, lg_360, lg_361, lg_362, lg_363, \
                         lg_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * ig_195[k]
                   + f_0 * lg_360[k];

        t_271[k] = -2.0 * ig_196[k]
                   + f_0 * lg_361[k];

        t_272[k] = -2.0 * ig_197[k]
                   + f_0 * lg_362[k];

        t_273[k] = -2.0 * ig_198[k]
                   + f_0 * lg_363[k];

        t_274[k] = -2.0 * ig_199[k]
                   + f_0 * lg_364[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ig_200, ig_201, ig_202, ig_203, \
                         ig_204, lg_365, lg_366, lg_367, lg_368, \
                         lg_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * ig_200[k]
                   + f_0 * lg_365[k];

        t_276[k] = -2.0 * ig_201[k]
                   + f_0 * lg_366[k];

        t_277[k] = -2.0 * ig_202[k]
                   + f_0 * lg_367[k];

        t_278[k] = -2.0 * ig_203[k]
                   + f_0 * lg_368[k];

        t_279[k] = -2.0 * ig_204[k]
                   + f_0 * lg_369[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ig_205, ig_206, ig_207, ig_208, \
                         ig_209, lg_370, lg_371, lg_372, lg_373, \
                         lg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * ig_205[k]
                   + f_0 * lg_370[k];

        t_281[k] = -2.0 * ig_206[k]
                   + f_0 * lg_371[k];

        t_282[k] = -2.0 * ig_207[k]
                   + f_0 * lg_372[k];

        t_283[k] = -2.0 * ig_208[k]
                   + f_0 * lg_373[k];

        t_284[k] = -2.0 * ig_209[k]
                   + f_0 * lg_374[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ig_210, ig_211, ig_212, ig_213, \
                         ig_214, lg_375, lg_376, lg_377, lg_378, \
                         lg_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -ig_210[k]
                   + f_0 * lg_375[k];

        t_286[k] = -ig_211[k]
                   + f_0 * lg_376[k];

        t_287[k] = -ig_212[k]
                   + f_0 * lg_377[k];

        t_288[k] = -ig_213[k]
                   + f_0 * lg_378[k];

        t_289[k] = -ig_214[k]
                   + f_0 * lg_379[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ig_215, ig_216, ig_217, ig_218, \
                         ig_219, lg_380, lg_381, lg_382, lg_383, \
                         lg_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -ig_215[k]
                   + f_0 * lg_380[k];

        t_291[k] = -ig_216[k]
                   + f_0 * lg_381[k];

        t_292[k] = -ig_217[k]
                   + f_0 * lg_382[k];

        t_293[k] = -ig_218[k]
                   + f_0 * lg_383[k];

        t_294[k] = -ig_219[k]
                   + f_0 * lg_384[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ig_220, ig_221, ig_222, ig_223, \
                         ig_224, lg_385, lg_386, lg_387, lg_388, \
                         lg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -ig_220[k]
                   + f_0 * lg_385[k];

        t_296[k] = -ig_221[k]
                   + f_0 * lg_386[k];

        t_297[k] = -ig_222[k]
                   + f_0 * lg_387[k];

        t_298[k] = -ig_223[k]
                   + f_0 * lg_388[k];

        t_299[k] = -ig_224[k]
                   + f_0 * lg_389[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, t_305, t_306, t_307, lg_390, \
                         lg_391, lg_392, lg_393, lg_394, lg_395, lg_396, \
                         lg_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_0 * lg_390[k];

        t_301[k] = f_0 * lg_391[k];

        t_302[k] = f_0 * lg_392[k];

        t_303[k] = f_0 * lg_393[k];

        t_304[k] = f_0 * lg_394[k];

        t_305[k] = f_0 * lg_395[k];

        t_306[k] = f_0 * lg_396[k];

        t_307[k] = f_0 * lg_397[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, t_313, t_314, lg_398, lg_399, \
                         lg_400, lg_401, lg_402, lg_403, lg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_0 * lg_398[k];

        t_309[k] = f_0 * lg_399[k];

        t_310[k] = f_0 * lg_400[k];

        t_311[k] = f_0 * lg_401[k];

        t_312[k] = f_0 * lg_402[k];

        t_313[k] = f_0 * lg_403[k];

        t_314[k] = f_0 * lg_404[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ig_225, ig_226, ig_227, ig_228, \
                         ig_229, lg_420, lg_421, lg_422, lg_423, \
                         lg_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -6.0 * ig_225[k]
                   + f_0 * lg_420[k];

        t_316[k] = -6.0 * ig_226[k]
                   + f_0 * lg_421[k];

        t_317[k] = -6.0 * ig_227[k]
                   + f_0 * lg_422[k];

        t_318[k] = -6.0 * ig_228[k]
                   + f_0 * lg_423[k];

        t_319[k] = -6.0 * ig_229[k]
                   + f_0 * lg_424[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ig_230, ig_231, ig_232, ig_233, \
                         ig_234, lg_425, lg_426, lg_427, lg_428, \
                         lg_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -6.0 * ig_230[k]
                   + f_0 * lg_425[k];

        t_321[k] = -6.0 * ig_231[k]
                   + f_0 * lg_426[k];

        t_322[k] = -6.0 * ig_232[k]
                   + f_0 * lg_427[k];

        t_323[k] = -6.0 * ig_233[k]
                   + f_0 * lg_428[k];

        t_324[k] = -6.0 * ig_234[k]
                   + f_0 * lg_429[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ig_235, ig_236, ig_237, ig_238, \
                         ig_239, lg_430, lg_431, lg_432, lg_433, \
                         lg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -6.0 * ig_235[k]
                   + f_0 * lg_430[k];

        t_326[k] = -6.0 * ig_236[k]
                   + f_0 * lg_431[k];

        t_327[k] = -6.0 * ig_237[k]
                   + f_0 * lg_432[k];

        t_328[k] = -6.0 * ig_238[k]
                   + f_0 * lg_433[k];

        t_329[k] = -6.0 * ig_239[k]
                   + f_0 * lg_434[k];
    }
}

static auto
compute_prim_geom_10_kg_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ig, const size_t lg,
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

    const auto *lg_435 = buffer.data(lg + 435);
    const auto *lg_436 = buffer.data(lg + 436);
    const auto *lg_437 = buffer.data(lg + 437);
    const auto *lg_438 = buffer.data(lg + 438);
    const auto *lg_439 = buffer.data(lg + 439);
    const auto *lg_440 = buffer.data(lg + 440);
    const auto *lg_441 = buffer.data(lg + 441);
    const auto *lg_442 = buffer.data(lg + 442);
    const auto *lg_443 = buffer.data(lg + 443);
    const auto *lg_444 = buffer.data(lg + 444);
    const auto *lg_445 = buffer.data(lg + 445);
    const auto *lg_446 = buffer.data(lg + 446);
    const auto *lg_447 = buffer.data(lg + 447);
    const auto *lg_448 = buffer.data(lg + 448);
    const auto *lg_449 = buffer.data(lg + 449);
    const auto *lg_450 = buffer.data(lg + 450);
    const auto *lg_451 = buffer.data(lg + 451);
    const auto *lg_452 = buffer.data(lg + 452);
    const auto *lg_453 = buffer.data(lg + 453);
    const auto *lg_454 = buffer.data(lg + 454);
    const auto *lg_455 = buffer.data(lg + 455);
    const auto *lg_456 = buffer.data(lg + 456);
    const auto *lg_457 = buffer.data(lg + 457);
    const auto *lg_458 = buffer.data(lg + 458);
    const auto *lg_459 = buffer.data(lg + 459);
    const auto *lg_460 = buffer.data(lg + 460);
    const auto *lg_461 = buffer.data(lg + 461);
    const auto *lg_462 = buffer.data(lg + 462);
    const auto *lg_463 = buffer.data(lg + 463);
    const auto *lg_464 = buffer.data(lg + 464);
    const auto *lg_465 = buffer.data(lg + 465);
    const auto *lg_466 = buffer.data(lg + 466);
    const auto *lg_467 = buffer.data(lg + 467);
    const auto *lg_468 = buffer.data(lg + 468);
    const auto *lg_469 = buffer.data(lg + 469);
    const auto *lg_470 = buffer.data(lg + 470);
    const auto *lg_471 = buffer.data(lg + 471);
    const auto *lg_472 = buffer.data(lg + 472);
    const auto *lg_473 = buffer.data(lg + 473);
    const auto *lg_474 = buffer.data(lg + 474);
    const auto *lg_475 = buffer.data(lg + 475);
    const auto *lg_476 = buffer.data(lg + 476);
    const auto *lg_477 = buffer.data(lg + 477);
    const auto *lg_478 = buffer.data(lg + 478);
    const auto *lg_479 = buffer.data(lg + 479);
    const auto *lg_480 = buffer.data(lg + 480);
    const auto *lg_481 = buffer.data(lg + 481);
    const auto *lg_482 = buffer.data(lg + 482);
    const auto *lg_483 = buffer.data(lg + 483);
    const auto *lg_484 = buffer.data(lg + 484);
    const auto *lg_485 = buffer.data(lg + 485);
    const auto *lg_486 = buffer.data(lg + 486);
    const auto *lg_487 = buffer.data(lg + 487);
    const auto *lg_488 = buffer.data(lg + 488);
    const auto *lg_489 = buffer.data(lg + 489);
    const auto *lg_490 = buffer.data(lg + 490);
    const auto *lg_491 = buffer.data(lg + 491);
    const auto *lg_492 = buffer.data(lg + 492);
    const auto *lg_493 = buffer.data(lg + 493);
    const auto *lg_494 = buffer.data(lg + 494);
    const auto *lg_495 = buffer.data(lg + 495);
    const auto *lg_496 = buffer.data(lg + 496);
    const auto *lg_497 = buffer.data(lg + 497);
    const auto *lg_498 = buffer.data(lg + 498);
    const auto *lg_499 = buffer.data(lg + 499);
    const auto *lg_500 = buffer.data(lg + 500);
    const auto *lg_501 = buffer.data(lg + 501);
    const auto *lg_502 = buffer.data(lg + 502);
    const auto *lg_503 = buffer.data(lg + 503);
    const auto *lg_504 = buffer.data(lg + 504);
    const auto *lg_505 = buffer.data(lg + 505);
    const auto *lg_506 = buffer.data(lg + 506);
    const auto *lg_507 = buffer.data(lg + 507);
    const auto *lg_508 = buffer.data(lg + 508);
    const auto *lg_509 = buffer.data(lg + 509);
    const auto *lg_510 = buffer.data(lg + 510);
    const auto *lg_511 = buffer.data(lg + 511);
    const auto *lg_512 = buffer.data(lg + 512);
    const auto *lg_513 = buffer.data(lg + 513);
    const auto *lg_514 = buffer.data(lg + 514);
    const auto *lg_515 = buffer.data(lg + 515);
    const auto *lg_516 = buffer.data(lg + 516);
    const auto *lg_517 = buffer.data(lg + 517);
    const auto *lg_518 = buffer.data(lg + 518);
    const auto *lg_519 = buffer.data(lg + 519);
    const auto *lg_520 = buffer.data(lg + 520);
    const auto *lg_521 = buffer.data(lg + 521);
    const auto *lg_522 = buffer.data(lg + 522);
    const auto *lg_523 = buffer.data(lg + 523);
    const auto *lg_524 = buffer.data(lg + 524);
    const auto *lg_540 = buffer.data(lg + 540);
    const auto *lg_541 = buffer.data(lg + 541);
    const auto *lg_542 = buffer.data(lg + 542);
    const auto *lg_543 = buffer.data(lg + 543);
    const auto *lg_544 = buffer.data(lg + 544);
    const auto *lg_545 = buffer.data(lg + 545);
    const auto *lg_546 = buffer.data(lg + 546);
    const auto *lg_547 = buffer.data(lg + 547);
    const auto *lg_548 = buffer.data(lg + 548);
    const auto *lg_549 = buffer.data(lg + 549);
    const auto *lg_550 = buffer.data(lg + 550);
    const auto *lg_551 = buffer.data(lg + 551);
    const auto *lg_552 = buffer.data(lg + 552);
    const auto *lg_553 = buffer.data(lg + 553);
    const auto *lg_554 = buffer.data(lg + 554);
    const auto *lg_555 = buffer.data(lg + 555);
    const auto *lg_556 = buffer.data(lg + 556);
    const auto *lg_557 = buffer.data(lg + 557);
    const auto *lg_558 = buffer.data(lg + 558);
    const auto *lg_559 = buffer.data(lg + 559);
    const auto *lg_560 = buffer.data(lg + 560);
    const auto *lg_561 = buffer.data(lg + 561);
    const auto *lg_562 = buffer.data(lg + 562);
    const auto *lg_563 = buffer.data(lg + 563);
    const auto *lg_564 = buffer.data(lg + 564);
    const auto *lg_565 = buffer.data(lg + 565);
    const auto *lg_566 = buffer.data(lg + 566);
    const auto *lg_567 = buffer.data(lg + 567);
    const auto *lg_568 = buffer.data(lg + 568);
    const auto *lg_569 = buffer.data(lg + 569);
    const auto *lg_570 = buffer.data(lg + 570);
    const auto *lg_571 = buffer.data(lg + 571);
    const auto *lg_572 = buffer.data(lg + 572);
    const auto *lg_573 = buffer.data(lg + 573);
    const auto *lg_574 = buffer.data(lg + 574);
    const auto *lg_575 = buffer.data(lg + 575);
    const auto *lg_576 = buffer.data(lg + 576);
    const auto *lg_577 = buffer.data(lg + 577);
    const auto *lg_578 = buffer.data(lg + 578);
    const auto *lg_579 = buffer.data(lg + 579);
    const auto *lg_580 = buffer.data(lg + 580);
    const auto *lg_581 = buffer.data(lg + 581);
    const auto *lg_582 = buffer.data(lg + 582);
    const auto *lg_583 = buffer.data(lg + 583);
    const auto *lg_584 = buffer.data(lg + 584);
    const auto *lg_585 = buffer.data(lg + 585);
    const auto *lg_586 = buffer.data(lg + 586);
    const auto *lg_587 = buffer.data(lg + 587);
    const auto *lg_588 = buffer.data(lg + 588);
    const auto *lg_589 = buffer.data(lg + 589);
    const auto *lg_590 = buffer.data(lg + 590);
    const auto *lg_591 = buffer.data(lg + 591);
    const auto *lg_592 = buffer.data(lg + 592);
    const auto *lg_593 = buffer.data(lg + 593);
    const auto *lg_594 = buffer.data(lg + 594);
    const auto *lg_595 = buffer.data(lg + 595);
    const auto *lg_596 = buffer.data(lg + 596);
    const auto *lg_597 = buffer.data(lg + 597);
    const auto *lg_598 = buffer.data(lg + 598);
    const auto *lg_599 = buffer.data(lg + 599);
    const auto *lg_600 = buffer.data(lg + 600);
    const auto *lg_601 = buffer.data(lg + 601);
    const auto *lg_602 = buffer.data(lg + 602);
    const auto *lg_603 = buffer.data(lg + 603);
    const auto *lg_604 = buffer.data(lg + 604);

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ig_240, ig_241, ig_242, ig_243, \
                         ig_244, lg_435, lg_436, lg_437, lg_438, \
                         lg_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -5.0 * ig_240[k]
                   + f_0 * lg_435[k];

        t_331[k] = -5.0 * ig_241[k]
                   + f_0 * lg_436[k];

        t_332[k] = -5.0 * ig_242[k]
                   + f_0 * lg_437[k];

        t_333[k] = -5.0 * ig_243[k]
                   + f_0 * lg_438[k];

        t_334[k] = -5.0 * ig_244[k]
                   + f_0 * lg_439[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ig_245, ig_246, ig_247, ig_248, \
                         ig_249, lg_440, lg_441, lg_442, lg_443, \
                         lg_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -5.0 * ig_245[k]
                   + f_0 * lg_440[k];

        t_336[k] = -5.0 * ig_246[k]
                   + f_0 * lg_441[k];

        t_337[k] = -5.0 * ig_247[k]
                   + f_0 * lg_442[k];

        t_338[k] = -5.0 * ig_248[k]
                   + f_0 * lg_443[k];

        t_339[k] = -5.0 * ig_249[k]
                   + f_0 * lg_444[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ig_250, ig_251, ig_252, ig_253, \
                         ig_254, lg_445, lg_446, lg_447, lg_448, \
                         lg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -5.0 * ig_250[k]
                   + f_0 * lg_445[k];

        t_341[k] = -5.0 * ig_251[k]
                   + f_0 * lg_446[k];

        t_342[k] = -5.0 * ig_252[k]
                   + f_0 * lg_447[k];

        t_343[k] = -5.0 * ig_253[k]
                   + f_0 * lg_448[k];

        t_344[k] = -5.0 * ig_254[k]
                   + f_0 * lg_449[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ig_255, ig_256, ig_257, ig_258, \
                         ig_259, lg_450, lg_451, lg_452, lg_453, \
                         lg_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -4.0 * ig_255[k]
                   + f_0 * lg_450[k];

        t_346[k] = -4.0 * ig_256[k]
                   + f_0 * lg_451[k];

        t_347[k] = -4.0 * ig_257[k]
                   + f_0 * lg_452[k];

        t_348[k] = -4.0 * ig_258[k]
                   + f_0 * lg_453[k];

        t_349[k] = -4.0 * ig_259[k]
                   + f_0 * lg_454[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ig_260, ig_261, ig_262, ig_263, \
                         ig_264, lg_455, lg_456, lg_457, lg_458, \
                         lg_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -4.0 * ig_260[k]
                   + f_0 * lg_455[k];

        t_351[k] = -4.0 * ig_261[k]
                   + f_0 * lg_456[k];

        t_352[k] = -4.0 * ig_262[k]
                   + f_0 * lg_457[k];

        t_353[k] = -4.0 * ig_263[k]
                   + f_0 * lg_458[k];

        t_354[k] = -4.0 * ig_264[k]
                   + f_0 * lg_459[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ig_265, ig_266, ig_267, ig_268, \
                         ig_269, lg_460, lg_461, lg_462, lg_463, \
                         lg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -4.0 * ig_265[k]
                   + f_0 * lg_460[k];

        t_356[k] = -4.0 * ig_266[k]
                   + f_0 * lg_461[k];

        t_357[k] = -4.0 * ig_267[k]
                   + f_0 * lg_462[k];

        t_358[k] = -4.0 * ig_268[k]
                   + f_0 * lg_463[k];

        t_359[k] = -4.0 * ig_269[k]
                   + f_0 * lg_464[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ig_270, ig_271, ig_272, ig_273, \
                         ig_274, lg_465, lg_466, lg_467, lg_468, \
                         lg_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -3.0 * ig_270[k]
                   + f_0 * lg_465[k];

        t_361[k] = -3.0 * ig_271[k]
                   + f_0 * lg_466[k];

        t_362[k] = -3.0 * ig_272[k]
                   + f_0 * lg_467[k];

        t_363[k] = -3.0 * ig_273[k]
                   + f_0 * lg_468[k];

        t_364[k] = -3.0 * ig_274[k]
                   + f_0 * lg_469[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ig_275, ig_276, ig_277, ig_278, \
                         ig_279, lg_470, lg_471, lg_472, lg_473, \
                         lg_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -3.0 * ig_275[k]
                   + f_0 * lg_470[k];

        t_366[k] = -3.0 * ig_276[k]
                   + f_0 * lg_471[k];

        t_367[k] = -3.0 * ig_277[k]
                   + f_0 * lg_472[k];

        t_368[k] = -3.0 * ig_278[k]
                   + f_0 * lg_473[k];

        t_369[k] = -3.0 * ig_279[k]
                   + f_0 * lg_474[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ig_280, ig_281, ig_282, ig_283, \
                         ig_284, lg_475, lg_476, lg_477, lg_478, \
                         lg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -3.0 * ig_280[k]
                   + f_0 * lg_475[k];

        t_371[k] = -3.0 * ig_281[k]
                   + f_0 * lg_476[k];

        t_372[k] = -3.0 * ig_282[k]
                   + f_0 * lg_477[k];

        t_373[k] = -3.0 * ig_283[k]
                   + f_0 * lg_478[k];

        t_374[k] = -3.0 * ig_284[k]
                   + f_0 * lg_479[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ig_285, ig_286, ig_287, ig_288, \
                         ig_289, lg_480, lg_481, lg_482, lg_483, \
                         lg_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -2.0 * ig_285[k]
                   + f_0 * lg_480[k];

        t_376[k] = -2.0 * ig_286[k]
                   + f_0 * lg_481[k];

        t_377[k] = -2.0 * ig_287[k]
                   + f_0 * lg_482[k];

        t_378[k] = -2.0 * ig_288[k]
                   + f_0 * lg_483[k];

        t_379[k] = -2.0 * ig_289[k]
                   + f_0 * lg_484[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ig_290, ig_291, ig_292, ig_293, \
                         ig_294, lg_485, lg_486, lg_487, lg_488, \
                         lg_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * ig_290[k]
                   + f_0 * lg_485[k];

        t_381[k] = -2.0 * ig_291[k]
                   + f_0 * lg_486[k];

        t_382[k] = -2.0 * ig_292[k]
                   + f_0 * lg_487[k];

        t_383[k] = -2.0 * ig_293[k]
                   + f_0 * lg_488[k];

        t_384[k] = -2.0 * ig_294[k]
                   + f_0 * lg_489[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ig_295, ig_296, ig_297, ig_298, \
                         ig_299, lg_490, lg_491, lg_492, lg_493, \
                         lg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * ig_295[k]
                   + f_0 * lg_490[k];

        t_386[k] = -2.0 * ig_296[k]
                   + f_0 * lg_491[k];

        t_387[k] = -2.0 * ig_297[k]
                   + f_0 * lg_492[k];

        t_388[k] = -2.0 * ig_298[k]
                   + f_0 * lg_493[k];

        t_389[k] = -2.0 * ig_299[k]
                   + f_0 * lg_494[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ig_300, ig_301, ig_302, ig_303, \
                         ig_304, lg_495, lg_496, lg_497, lg_498, \
                         lg_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -ig_300[k]
                   + f_0 * lg_495[k];

        t_391[k] = -ig_301[k]
                   + f_0 * lg_496[k];

        t_392[k] = -ig_302[k]
                   + f_0 * lg_497[k];

        t_393[k] = -ig_303[k]
                   + f_0 * lg_498[k];

        t_394[k] = -ig_304[k]
                   + f_0 * lg_499[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ig_305, ig_306, ig_307, ig_308, \
                         ig_309, lg_500, lg_501, lg_502, lg_503, \
                         lg_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -ig_305[k]
                   + f_0 * lg_500[k];

        t_396[k] = -ig_306[k]
                   + f_0 * lg_501[k];

        t_397[k] = -ig_307[k]
                   + f_0 * lg_502[k];

        t_398[k] = -ig_308[k]
                   + f_0 * lg_503[k];

        t_399[k] = -ig_309[k]
                   + f_0 * lg_504[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ig_310, ig_311, ig_312, ig_313, \
                         ig_314, lg_505, lg_506, lg_507, lg_508, \
                         lg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -ig_310[k]
                   + f_0 * lg_505[k];

        t_401[k] = -ig_311[k]
                   + f_0 * lg_506[k];

        t_402[k] = -ig_312[k]
                   + f_0 * lg_507[k];

        t_403[k] = -ig_313[k]
                   + f_0 * lg_508[k];

        t_404[k] = -ig_314[k]
                   + f_0 * lg_509[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, lg_510, \
                         lg_511, lg_512, lg_513, lg_514, lg_515, lg_516, \
                         lg_517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * lg_510[k];

        t_406[k] = f_0 * lg_511[k];

        t_407[k] = f_0 * lg_512[k];

        t_408[k] = f_0 * lg_513[k];

        t_409[k] = f_0 * lg_514[k];

        t_410[k] = f_0 * lg_515[k];

        t_411[k] = f_0 * lg_516[k];

        t_412[k] = f_0 * lg_517[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, lg_518, lg_519, \
                         lg_520, lg_521, lg_522, lg_523, lg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * lg_518[k];

        t_414[k] = f_0 * lg_519[k];

        t_415[k] = f_0 * lg_520[k];

        t_416[k] = f_0 * lg_521[k];

        t_417[k] = f_0 * lg_522[k];

        t_418[k] = f_0 * lg_523[k];

        t_419[k] = f_0 * lg_524[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ig_315, ig_316, ig_317, ig_318, \
                         ig_319, lg_540, lg_541, lg_542, lg_543, \
                         lg_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -7.0 * ig_315[k]
                   + f_0 * lg_540[k];

        t_421[k] = -7.0 * ig_316[k]
                   + f_0 * lg_541[k];

        t_422[k] = -7.0 * ig_317[k]
                   + f_0 * lg_542[k];

        t_423[k] = -7.0 * ig_318[k]
                   + f_0 * lg_543[k];

        t_424[k] = -7.0 * ig_319[k]
                   + f_0 * lg_544[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ig_320, ig_321, ig_322, ig_323, \
                         ig_324, lg_545, lg_546, lg_547, lg_548, \
                         lg_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -7.0 * ig_320[k]
                   + f_0 * lg_545[k];

        t_426[k] = -7.0 * ig_321[k]
                   + f_0 * lg_546[k];

        t_427[k] = -7.0 * ig_322[k]
                   + f_0 * lg_547[k];

        t_428[k] = -7.0 * ig_323[k]
                   + f_0 * lg_548[k];

        t_429[k] = -7.0 * ig_324[k]
                   + f_0 * lg_549[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ig_325, ig_326, ig_327, ig_328, \
                         ig_329, lg_550, lg_551, lg_552, lg_553, \
                         lg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -7.0 * ig_325[k]
                   + f_0 * lg_550[k];

        t_431[k] = -7.0 * ig_326[k]
                   + f_0 * lg_551[k];

        t_432[k] = -7.0 * ig_327[k]
                   + f_0 * lg_552[k];

        t_433[k] = -7.0 * ig_328[k]
                   + f_0 * lg_553[k];

        t_434[k] = -7.0 * ig_329[k]
                   + f_0 * lg_554[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ig_330, ig_331, ig_332, ig_333, \
                         ig_334, lg_555, lg_556, lg_557, lg_558, \
                         lg_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -6.0 * ig_330[k]
                   + f_0 * lg_555[k];

        t_436[k] = -6.0 * ig_331[k]
                   + f_0 * lg_556[k];

        t_437[k] = -6.0 * ig_332[k]
                   + f_0 * lg_557[k];

        t_438[k] = -6.0 * ig_333[k]
                   + f_0 * lg_558[k];

        t_439[k] = -6.0 * ig_334[k]
                   + f_0 * lg_559[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ig_335, ig_336, ig_337, ig_338, \
                         ig_339, lg_560, lg_561, lg_562, lg_563, \
                         lg_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -6.0 * ig_335[k]
                   + f_0 * lg_560[k];

        t_441[k] = -6.0 * ig_336[k]
                   + f_0 * lg_561[k];

        t_442[k] = -6.0 * ig_337[k]
                   + f_0 * lg_562[k];

        t_443[k] = -6.0 * ig_338[k]
                   + f_0 * lg_563[k];

        t_444[k] = -6.0 * ig_339[k]
                   + f_0 * lg_564[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ig_340, ig_341, ig_342, ig_343, \
                         ig_344, lg_565, lg_566, lg_567, lg_568, \
                         lg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -6.0 * ig_340[k]
                   + f_0 * lg_565[k];

        t_446[k] = -6.0 * ig_341[k]
                   + f_0 * lg_566[k];

        t_447[k] = -6.0 * ig_342[k]
                   + f_0 * lg_567[k];

        t_448[k] = -6.0 * ig_343[k]
                   + f_0 * lg_568[k];

        t_449[k] = -6.0 * ig_344[k]
                   + f_0 * lg_569[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ig_345, ig_346, ig_347, ig_348, \
                         ig_349, lg_570, lg_571, lg_572, lg_573, \
                         lg_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -5.0 * ig_345[k]
                   + f_0 * lg_570[k];

        t_451[k] = -5.0 * ig_346[k]
                   + f_0 * lg_571[k];

        t_452[k] = -5.0 * ig_347[k]
                   + f_0 * lg_572[k];

        t_453[k] = -5.0 * ig_348[k]
                   + f_0 * lg_573[k];

        t_454[k] = -5.0 * ig_349[k]
                   + f_0 * lg_574[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ig_350, ig_351, ig_352, ig_353, \
                         ig_354, lg_575, lg_576, lg_577, lg_578, \
                         lg_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -5.0 * ig_350[k]
                   + f_0 * lg_575[k];

        t_456[k] = -5.0 * ig_351[k]
                   + f_0 * lg_576[k];

        t_457[k] = -5.0 * ig_352[k]
                   + f_0 * lg_577[k];

        t_458[k] = -5.0 * ig_353[k]
                   + f_0 * lg_578[k];

        t_459[k] = -5.0 * ig_354[k]
                   + f_0 * lg_579[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ig_355, ig_356, ig_357, ig_358, \
                         ig_359, lg_580, lg_581, lg_582, lg_583, \
                         lg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -5.0 * ig_355[k]
                   + f_0 * lg_580[k];

        t_461[k] = -5.0 * ig_356[k]
                   + f_0 * lg_581[k];

        t_462[k] = -5.0 * ig_357[k]
                   + f_0 * lg_582[k];

        t_463[k] = -5.0 * ig_358[k]
                   + f_0 * lg_583[k];

        t_464[k] = -5.0 * ig_359[k]
                   + f_0 * lg_584[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ig_360, ig_361, ig_362, ig_363, \
                         ig_364, lg_585, lg_586, lg_587, lg_588, \
                         lg_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -4.0 * ig_360[k]
                   + f_0 * lg_585[k];

        t_466[k] = -4.0 * ig_361[k]
                   + f_0 * lg_586[k];

        t_467[k] = -4.0 * ig_362[k]
                   + f_0 * lg_587[k];

        t_468[k] = -4.0 * ig_363[k]
                   + f_0 * lg_588[k];

        t_469[k] = -4.0 * ig_364[k]
                   + f_0 * lg_589[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ig_365, ig_366, ig_367, ig_368, \
                         ig_369, lg_590, lg_591, lg_592, lg_593, \
                         lg_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -4.0 * ig_365[k]
                   + f_0 * lg_590[k];

        t_471[k] = -4.0 * ig_366[k]
                   + f_0 * lg_591[k];

        t_472[k] = -4.0 * ig_367[k]
                   + f_0 * lg_592[k];

        t_473[k] = -4.0 * ig_368[k]
                   + f_0 * lg_593[k];

        t_474[k] = -4.0 * ig_369[k]
                   + f_0 * lg_594[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ig_370, ig_371, ig_372, ig_373, \
                         ig_374, lg_595, lg_596, lg_597, lg_598, \
                         lg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -4.0 * ig_370[k]
                   + f_0 * lg_595[k];

        t_476[k] = -4.0 * ig_371[k]
                   + f_0 * lg_596[k];

        t_477[k] = -4.0 * ig_372[k]
                   + f_0 * lg_597[k];

        t_478[k] = -4.0 * ig_373[k]
                   + f_0 * lg_598[k];

        t_479[k] = -4.0 * ig_374[k]
                   + f_0 * lg_599[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ig_375, ig_376, ig_377, ig_378, \
                         ig_379, lg_600, lg_601, lg_602, lg_603, \
                         lg_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -3.0 * ig_375[k]
                   + f_0 * lg_600[k];

        t_481[k] = -3.0 * ig_376[k]
                   + f_0 * lg_601[k];

        t_482[k] = -3.0 * ig_377[k]
                   + f_0 * lg_602[k];

        t_483[k] = -3.0 * ig_378[k]
                   + f_0 * lg_603[k];

        t_484[k] = -3.0 * ig_379[k]
                   + f_0 * lg_604[k];
    }
}

static auto
compute_prim_geom_10_kg_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ig, const size_t lg,
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

    const auto *lg_605 = buffer.data(lg + 605);
    const auto *lg_606 = buffer.data(lg + 606);
    const auto *lg_607 = buffer.data(lg + 607);
    const auto *lg_608 = buffer.data(lg + 608);
    const auto *lg_609 = buffer.data(lg + 609);
    const auto *lg_610 = buffer.data(lg + 610);
    const auto *lg_611 = buffer.data(lg + 611);
    const auto *lg_612 = buffer.data(lg + 612);
    const auto *lg_613 = buffer.data(lg + 613);
    const auto *lg_614 = buffer.data(lg + 614);
    const auto *lg_615 = buffer.data(lg + 615);
    const auto *lg_616 = buffer.data(lg + 616);
    const auto *lg_617 = buffer.data(lg + 617);
    const auto *lg_618 = buffer.data(lg + 618);
    const auto *lg_619 = buffer.data(lg + 619);
    const auto *lg_620 = buffer.data(lg + 620);
    const auto *lg_621 = buffer.data(lg + 621);
    const auto *lg_622 = buffer.data(lg + 622);
    const auto *lg_623 = buffer.data(lg + 623);
    const auto *lg_624 = buffer.data(lg + 624);
    const auto *lg_625 = buffer.data(lg + 625);
    const auto *lg_626 = buffer.data(lg + 626);
    const auto *lg_627 = buffer.data(lg + 627);
    const auto *lg_628 = buffer.data(lg + 628);
    const auto *lg_629 = buffer.data(lg + 629);
    const auto *lg_630 = buffer.data(lg + 630);
    const auto *lg_631 = buffer.data(lg + 631);
    const auto *lg_632 = buffer.data(lg + 632);
    const auto *lg_633 = buffer.data(lg + 633);
    const auto *lg_634 = buffer.data(lg + 634);
    const auto *lg_635 = buffer.data(lg + 635);
    const auto *lg_636 = buffer.data(lg + 636);
    const auto *lg_637 = buffer.data(lg + 637);
    const auto *lg_638 = buffer.data(lg + 638);
    const auto *lg_639 = buffer.data(lg + 639);
    const auto *lg_640 = buffer.data(lg + 640);
    const auto *lg_641 = buffer.data(lg + 641);
    const auto *lg_642 = buffer.data(lg + 642);
    const auto *lg_643 = buffer.data(lg + 643);
    const auto *lg_644 = buffer.data(lg + 644);
    const auto *lg_645 = buffer.data(lg + 645);
    const auto *lg_646 = buffer.data(lg + 646);
    const auto *lg_647 = buffer.data(lg + 647);
    const auto *lg_648 = buffer.data(lg + 648);
    const auto *lg_649 = buffer.data(lg + 649);
    const auto *lg_650 = buffer.data(lg + 650);
    const auto *lg_651 = buffer.data(lg + 651);
    const auto *lg_652 = buffer.data(lg + 652);
    const auto *lg_653 = buffer.data(lg + 653);
    const auto *lg_654 = buffer.data(lg + 654);
    const auto *lg_655 = buffer.data(lg + 655);
    const auto *lg_656 = buffer.data(lg + 656);
    const auto *lg_657 = buffer.data(lg + 657);
    const auto *lg_658 = buffer.data(lg + 658);
    const auto *lg_659 = buffer.data(lg + 659);

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ig_380, ig_381, ig_382, ig_383, \
                         ig_384, lg_605, lg_606, lg_607, lg_608, \
                         lg_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -3.0 * ig_380[k]
                   + f_0 * lg_605[k];

        t_486[k] = -3.0 * ig_381[k]
                   + f_0 * lg_606[k];

        t_487[k] = -3.0 * ig_382[k]
                   + f_0 * lg_607[k];

        t_488[k] = -3.0 * ig_383[k]
                   + f_0 * lg_608[k];

        t_489[k] = -3.0 * ig_384[k]
                   + f_0 * lg_609[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ig_385, ig_386, ig_387, ig_388, \
                         ig_389, lg_610, lg_611, lg_612, lg_613, \
                         lg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -3.0 * ig_385[k]
                   + f_0 * lg_610[k];

        t_491[k] = -3.0 * ig_386[k]
                   + f_0 * lg_611[k];

        t_492[k] = -3.0 * ig_387[k]
                   + f_0 * lg_612[k];

        t_493[k] = -3.0 * ig_388[k]
                   + f_0 * lg_613[k];

        t_494[k] = -3.0 * ig_389[k]
                   + f_0 * lg_614[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ig_390, ig_391, ig_392, ig_393, \
                         ig_394, lg_615, lg_616, lg_617, lg_618, \
                         lg_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -2.0 * ig_390[k]
                   + f_0 * lg_615[k];

        t_496[k] = -2.0 * ig_391[k]
                   + f_0 * lg_616[k];

        t_497[k] = -2.0 * ig_392[k]
                   + f_0 * lg_617[k];

        t_498[k] = -2.0 * ig_393[k]
                   + f_0 * lg_618[k];

        t_499[k] = -2.0 * ig_394[k]
                   + f_0 * lg_619[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ig_395, ig_396, ig_397, ig_398, \
                         ig_399, lg_620, lg_621, lg_622, lg_623, \
                         lg_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -2.0 * ig_395[k]
                   + f_0 * lg_620[k];

        t_501[k] = -2.0 * ig_396[k]
                   + f_0 * lg_621[k];

        t_502[k] = -2.0 * ig_397[k]
                   + f_0 * lg_622[k];

        t_503[k] = -2.0 * ig_398[k]
                   + f_0 * lg_623[k];

        t_504[k] = -2.0 * ig_399[k]
                   + f_0 * lg_624[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ig_400, ig_401, ig_402, ig_403, \
                         ig_404, lg_625, lg_626, lg_627, lg_628, \
                         lg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -2.0 * ig_400[k]
                   + f_0 * lg_625[k];

        t_506[k] = -2.0 * ig_401[k]
                   + f_0 * lg_626[k];

        t_507[k] = -2.0 * ig_402[k]
                   + f_0 * lg_627[k];

        t_508[k] = -2.0 * ig_403[k]
                   + f_0 * lg_628[k];

        t_509[k] = -2.0 * ig_404[k]
                   + f_0 * lg_629[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ig_405, ig_406, ig_407, ig_408, \
                         ig_409, lg_630, lg_631, lg_632, lg_633, \
                         lg_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -ig_405[k]
                   + f_0 * lg_630[k];

        t_511[k] = -ig_406[k]
                   + f_0 * lg_631[k];

        t_512[k] = -ig_407[k]
                   + f_0 * lg_632[k];

        t_513[k] = -ig_408[k]
                   + f_0 * lg_633[k];

        t_514[k] = -ig_409[k]
                   + f_0 * lg_634[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ig_410, ig_411, ig_412, ig_413, \
                         ig_414, lg_635, lg_636, lg_637, lg_638, \
                         lg_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -ig_410[k]
                   + f_0 * lg_635[k];

        t_516[k] = -ig_411[k]
                   + f_0 * lg_636[k];

        t_517[k] = -ig_412[k]
                   + f_0 * lg_637[k];

        t_518[k] = -ig_413[k]
                   + f_0 * lg_638[k];

        t_519[k] = -ig_414[k]
                   + f_0 * lg_639[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ig_415, ig_416, ig_417, ig_418, \
                         ig_419, lg_640, lg_641, lg_642, lg_643, \
                         lg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -ig_415[k]
                   + f_0 * lg_640[k];

        t_521[k] = -ig_416[k]
                   + f_0 * lg_641[k];

        t_522[k] = -ig_417[k]
                   + f_0 * lg_642[k];

        t_523[k] = -ig_418[k]
                   + f_0 * lg_643[k];

        t_524[k] = -ig_419[k]
                   + f_0 * lg_644[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, t_530, t_531, t_532, lg_645, \
                         lg_646, lg_647, lg_648, lg_649, lg_650, lg_651, \
                         lg_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_0 * lg_645[k];

        t_526[k] = f_0 * lg_646[k];

        t_527[k] = f_0 * lg_647[k];

        t_528[k] = f_0 * lg_648[k];

        t_529[k] = f_0 * lg_649[k];

        t_530[k] = f_0 * lg_650[k];

        t_531[k] = f_0 * lg_651[k];

        t_532[k] = f_0 * lg_652[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, t_537, t_538, t_539, lg_653, lg_654, \
                         lg_655, lg_656, lg_657, lg_658, lg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_0 * lg_653[k];

        t_534[k] = f_0 * lg_654[k];

        t_535[k] = f_0 * lg_655[k];

        t_536[k] = f_0 * lg_656[k];

        t_537[k] = f_0 * lg_657[k];

        t_538[k] = f_0 * lg_658[k];

        t_539[k] = f_0 * lg_659[k];
    }
}

auto
compute_prim_geom_10_kg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t ig, const size_t lg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kg_electron_repulsion_1_piece0(buffer, target, ig, lg, ncols, alpha);

    compute_prim_geom_10_kg_electron_repulsion_1_piece1(buffer, target, ig, lg, ncols, alpha);

    compute_prim_geom_10_kg_electron_repulsion_1_piece2(buffer, target, ig, lg, ncols, alpha);

    compute_prim_geom_10_kg_electron_repulsion_1_piece3(buffer, target, ig, lg, ncols, alpha);
}

static auto
compute_prim_geom_10_kg_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ig, const size_t lg,
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

    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_31 = buffer.data(lg + 31);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_34 = buffer.data(lg + 34);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_36 = buffer.data(lg + 36);
    const auto *lg_37 = buffer.data(lg + 37);
    const auto *lg_38 = buffer.data(lg + 38);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_79 = buffer.data(lg + 79);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_81 = buffer.data(lg + 81);
    const auto *lg_82 = buffer.data(lg + 82);
    const auto *lg_83 = buffer.data(lg + 83);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_139 = buffer.data(lg + 139);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_141 = buffer.data(lg + 141);
    const auto *lg_142 = buffer.data(lg + 142);
    const auto *lg_143 = buffer.data(lg + 143);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);
    const auto *lg_168 = buffer.data(lg + 168);
    const auto *lg_169 = buffer.data(lg + 169);
    const auto *lg_170 = buffer.data(lg + 170);
    const auto *lg_171 = buffer.data(lg + 171);
    const auto *lg_172 = buffer.data(lg + 172);
    const auto *lg_173 = buffer.data(lg + 173);
    const auto *lg_174 = buffer.data(lg + 174);
    const auto *lg_175 = buffer.data(lg + 175);
    const auto *lg_176 = buffer.data(lg + 176);
    const auto *lg_177 = buffer.data(lg + 177);
    const auto *lg_178 = buffer.data(lg + 178);
    const auto *lg_179 = buffer.data(lg + 179);
    const auto *lg_180 = buffer.data(lg + 180);
    const auto *lg_181 = buffer.data(lg + 181);
    const auto *lg_182 = buffer.data(lg + 182);
    const auto *lg_183 = buffer.data(lg + 183);
    const auto *lg_184 = buffer.data(lg + 184);
    const auto *lg_185 = buffer.data(lg + 185);
    const auto *lg_186 = buffer.data(lg + 186);
    const auto *lg_187 = buffer.data(lg + 187);
    const auto *lg_188 = buffer.data(lg + 188);
    const auto *lg_189 = buffer.data(lg + 189);
    const auto *lg_190 = buffer.data(lg + 190);
    const auto *lg_191 = buffer.data(lg + 191);
    const auto *lg_192 = buffer.data(lg + 192);
    const auto *lg_193 = buffer.data(lg + 193);
    const auto *lg_194 = buffer.data(lg + 194);
    const auto *lg_195 = buffer.data(lg + 195);
    const auto *lg_196 = buffer.data(lg + 196);
    const auto *lg_197 = buffer.data(lg + 197);
    const auto *lg_198 = buffer.data(lg + 198);
    const auto *lg_199 = buffer.data(lg + 199);
    const auto *lg_200 = buffer.data(lg + 200);
    const auto *lg_201 = buffer.data(lg + 201);
    const auto *lg_202 = buffer.data(lg + 202);
    const auto *lg_203 = buffer.data(lg + 203);
    const auto *lg_204 = buffer.data(lg + 204);
    const auto *lg_205 = buffer.data(lg + 205);
    const auto *lg_206 = buffer.data(lg + 206);
    const auto *lg_207 = buffer.data(lg + 207);
    const auto *lg_208 = buffer.data(lg + 208);
    const auto *lg_209 = buffer.data(lg + 209);
    const auto *lg_210 = buffer.data(lg + 210);
    const auto *lg_211 = buffer.data(lg + 211);
    const auto *lg_212 = buffer.data(lg + 212);
    const auto *lg_213 = buffer.data(lg + 213);
    const auto *lg_214 = buffer.data(lg + 214);
    const auto *lg_215 = buffer.data(lg + 215);
    const auto *lg_216 = buffer.data(lg + 216);
    const auto *lg_217 = buffer.data(lg + 217);
    const auto *lg_218 = buffer.data(lg + 218);
    const auto *lg_219 = buffer.data(lg + 219);
    const auto *lg_220 = buffer.data(lg + 220);
    const auto *lg_221 = buffer.data(lg + 221);
    const auto *lg_222 = buffer.data(lg + 222);
    const auto *lg_223 = buffer.data(lg + 223);
    const auto *lg_224 = buffer.data(lg + 224);
    const auto *lg_240 = buffer.data(lg + 240);
    const auto *lg_241 = buffer.data(lg + 241);
    const auto *lg_242 = buffer.data(lg + 242);
    const auto *lg_243 = buffer.data(lg + 243);
    const auto *lg_244 = buffer.data(lg + 244);
    const auto *lg_245 = buffer.data(lg + 245);
    const auto *lg_246 = buffer.data(lg + 246);
    const auto *lg_247 = buffer.data(lg + 247);
    const auto *lg_248 = buffer.data(lg + 248);
    const auto *lg_249 = buffer.data(lg + 249);
    const auto *lg_250 = buffer.data(lg + 250);
    const auto *lg_251 = buffer.data(lg + 251);
    const auto *lg_252 = buffer.data(lg + 252);
    const auto *lg_253 = buffer.data(lg + 253);
    const auto *lg_254 = buffer.data(lg + 254);
    const auto *lg_255 = buffer.data(lg + 255);
    const auto *lg_256 = buffer.data(lg + 256);
    const auto *lg_257 = buffer.data(lg + 257);
    const auto *lg_258 = buffer.data(lg + 258);
    const auto *lg_259 = buffer.data(lg + 259);
    const auto *lg_260 = buffer.data(lg + 260);
    const auto *lg_261 = buffer.data(lg + 261);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, lg_30, lg_31, lg_32, lg_33, \
                         lg_34, lg_35, lg_36, lg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lg_30[k];

        t_1[k] = f_0 * lg_31[k];

        t_2[k] = f_0 * lg_32[k];

        t_3[k] = f_0 * lg_33[k];

        t_4[k] = f_0 * lg_34[k];

        t_5[k] = f_0 * lg_35[k];

        t_6[k] = f_0 * lg_36[k];

        t_7[k] = f_0 * lg_37[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, lg_38, lg_39, lg_40, \
                         lg_41, lg_42, lg_43, lg_44, lg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * lg_38[k];

        t_9[k] = f_0 * lg_39[k];

        t_10[k] = f_0 * lg_40[k];

        t_11[k] = f_0 * lg_41[k];

        t_12[k] = f_0 * lg_42[k];

        t_13[k] = f_0 * lg_43[k];

        t_14[k] = f_0 * lg_44[k];

        t_15[k] = f_0 * lg_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, lg_61, lg_62, lg_63, \
                         lg_64, lg_65, lg_66, lg_67, lg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * lg_61[k];

        t_17[k] = f_0 * lg_62[k];

        t_18[k] = f_0 * lg_63[k];

        t_19[k] = f_0 * lg_64[k];

        t_20[k] = f_0 * lg_65[k];

        t_21[k] = f_0 * lg_66[k];

        t_22[k] = f_0 * lg_67[k];

        t_23[k] = f_0 * lg_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, ig_0, lg_69, lg_70, lg_71, \
                         lg_72, lg_73, lg_74, lg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * lg_69[k];

        t_25[k] = f_0 * lg_70[k];

        t_26[k] = f_0 * lg_71[k];

        t_27[k] = f_0 * lg_72[k];

        t_28[k] = f_0 * lg_73[k];

        t_29[k] = f_0 * lg_74[k];

        t_30[k] = -ig_0[k]
                  + f_0 * lg_75[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ig_1, ig_2, ig_3, ig_4, ig_5, lg_76, \
                         lg_77, lg_78, lg_79, lg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -ig_1[k]
                  + f_0 * lg_76[k];

        t_32[k] = -ig_2[k]
                  + f_0 * lg_77[k];

        t_33[k] = -ig_3[k]
                  + f_0 * lg_78[k];

        t_34[k] = -ig_4[k]
                  + f_0 * lg_79[k];

        t_35[k] = -ig_5[k]
                  + f_0 * lg_80[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ig_6, ig_7, ig_8, ig_9, ig_10, lg_81, \
                         lg_82, lg_83, lg_84, lg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -ig_6[k]
                  + f_0 * lg_81[k];

        t_37[k] = -ig_7[k]
                  + f_0 * lg_82[k];

        t_38[k] = -ig_8[k]
                  + f_0 * lg_83[k];

        t_39[k] = -ig_9[k]
                  + f_0 * lg_84[k];

        t_40[k] = -ig_10[k]
                  + f_0 * lg_85[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, t_46, ig_11, ig_12, ig_13, ig_14, \
                         lg_86, lg_87, lg_88, lg_89, lg_105, lg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -ig_11[k]
                  + f_0 * lg_86[k];

        t_42[k] = -ig_12[k]
                  + f_0 * lg_87[k];

        t_43[k] = -ig_13[k]
                  + f_0 * lg_88[k];

        t_44[k] = -ig_14[k]
                  + f_0 * lg_89[k];

        t_45[k] = f_0 * lg_105[k];

        t_46[k] = f_0 * lg_106[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, t_53, t_54, lg_107, lg_108, \
                         lg_109, lg_110, lg_111, lg_112, lg_113, \
                         lg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * lg_107[k];

        t_48[k] = f_0 * lg_108[k];

        t_49[k] = f_0 * lg_109[k];

        t_50[k] = f_0 * lg_110[k];

        t_51[k] = f_0 * lg_111[k];

        t_52[k] = f_0 * lg_112[k];

        t_53[k] = f_0 * lg_113[k];

        t_54[k] = f_0 * lg_114[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, ig_15, ig_16, lg_115, \
                         lg_116, lg_117, lg_118, lg_119, lg_120, \
                         lg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * lg_115[k];

        t_56[k] = f_0 * lg_116[k];

        t_57[k] = f_0 * lg_117[k];

        t_58[k] = f_0 * lg_118[k];

        t_59[k] = f_0 * lg_119[k];

        t_60[k] = -ig_15[k]
                  + f_0 * lg_120[k];

        t_61[k] = -ig_16[k]
                  + f_0 * lg_121[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, ig_17, ig_18, ig_19, ig_20, ig_21, \
                         lg_122, lg_123, lg_124, lg_125, lg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -ig_17[k]
                  + f_0 * lg_122[k];

        t_63[k] = -ig_18[k]
                  + f_0 * lg_123[k];

        t_64[k] = -ig_19[k]
                  + f_0 * lg_124[k];

        t_65[k] = -ig_20[k]
                  + f_0 * lg_125[k];

        t_66[k] = -ig_21[k]
                  + f_0 * lg_126[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ig_22, ig_23, ig_24, ig_25, ig_26, \
                         lg_127, lg_128, lg_129, lg_130, lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -ig_22[k]
                  + f_0 * lg_127[k];

        t_68[k] = -ig_23[k]
                  + f_0 * lg_128[k];

        t_69[k] = -ig_24[k]
                  + f_0 * lg_129[k];

        t_70[k] = -ig_25[k]
                  + f_0 * lg_130[k];

        t_71[k] = -ig_26[k]
                  + f_0 * lg_131[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ig_27, ig_28, ig_29, ig_30, ig_31, \
                         lg_132, lg_133, lg_134, lg_135, lg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -ig_27[k]
                  + f_0 * lg_132[k];

        t_73[k] = -ig_28[k]
                  + f_0 * lg_133[k];

        t_74[k] = -ig_29[k]
                  + f_0 * lg_134[k];

        t_75[k] = -2.0 * ig_30[k]
                  + f_0 * lg_135[k];

        t_76[k] = -2.0 * ig_31[k]
                  + f_0 * lg_136[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, ig_32, ig_33, ig_34, ig_35, ig_36, \
                         lg_137, lg_138, lg_139, lg_140, lg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -2.0 * ig_32[k]
                  + f_0 * lg_137[k];

        t_78[k] = -2.0 * ig_33[k]
                  + f_0 * lg_138[k];

        t_79[k] = -2.0 * ig_34[k]
                  + f_0 * lg_139[k];

        t_80[k] = -2.0 * ig_35[k]
                  + f_0 * lg_140[k];

        t_81[k] = -2.0 * ig_36[k]
                  + f_0 * lg_141[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, ig_37, ig_38, ig_39, ig_40, ig_41, \
                         lg_142, lg_143, lg_144, lg_145, lg_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -2.0 * ig_37[k]
                  + f_0 * lg_142[k];

        t_83[k] = -2.0 * ig_38[k]
                  + f_0 * lg_143[k];

        t_84[k] = -2.0 * ig_39[k]
                  + f_0 * lg_144[k];

        t_85[k] = -2.0 * ig_40[k]
                  + f_0 * lg_145[k];

        t_86[k] = -2.0 * ig_41[k]
                  + f_0 * lg_146[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, ig_42, ig_43, ig_44, lg_147, \
                         lg_148, lg_149, lg_165, lg_166, lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -2.0 * ig_42[k]
                  + f_0 * lg_147[k];

        t_88[k] = -2.0 * ig_43[k]
                  + f_0 * lg_148[k];

        t_89[k] = -2.0 * ig_44[k]
                  + f_0 * lg_149[k];

        t_90[k] = f_0 * lg_165[k];

        t_91[k] = f_0 * lg_166[k];

        t_92[k] = f_0 * lg_167[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, lg_168, lg_169, \
                         lg_170, lg_171, lg_172, lg_173, lg_174, \
                         lg_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * lg_168[k];

        t_94[k] = f_0 * lg_169[k];

        t_95[k] = f_0 * lg_170[k];

        t_96[k] = f_0 * lg_171[k];

        t_97[k] = f_0 * lg_172[k];

        t_98[k] = f_0 * lg_173[k];

        t_99[k] = f_0 * lg_174[k];

        t_100[k] = f_0 * lg_175[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, ig_45, ig_46, lg_176, \
                         lg_177, lg_178, lg_179, lg_180, lg_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * lg_176[k];

        t_102[k] = f_0 * lg_177[k];

        t_103[k] = f_0 * lg_178[k];

        t_104[k] = f_0 * lg_179[k];

        t_105[k] = -ig_45[k]
                   + f_0 * lg_180[k];

        t_106[k] = -ig_46[k]
                   + f_0 * lg_181[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, ig_47, ig_48, ig_49, ig_50, ig_51, \
                         lg_182, lg_183, lg_184, lg_185, lg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -ig_47[k]
                   + f_0 * lg_182[k];

        t_108[k] = -ig_48[k]
                   + f_0 * lg_183[k];

        t_109[k] = -ig_49[k]
                   + f_0 * lg_184[k];

        t_110[k] = -ig_50[k]
                   + f_0 * lg_185[k];

        t_111[k] = -ig_51[k]
                   + f_0 * lg_186[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, ig_52, ig_53, ig_54, ig_55, ig_56, \
                         lg_187, lg_188, lg_189, lg_190, lg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -ig_52[k]
                   + f_0 * lg_187[k];

        t_113[k] = -ig_53[k]
                   + f_0 * lg_188[k];

        t_114[k] = -ig_54[k]
                   + f_0 * lg_189[k];

        t_115[k] = -ig_55[k]
                   + f_0 * lg_190[k];

        t_116[k] = -ig_56[k]
                   + f_0 * lg_191[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, ig_57, ig_58, ig_59, ig_60, ig_61, \
                         lg_192, lg_193, lg_194, lg_195, lg_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -ig_57[k]
                   + f_0 * lg_192[k];

        t_118[k] = -ig_58[k]
                   + f_0 * lg_193[k];

        t_119[k] = -ig_59[k]
                   + f_0 * lg_194[k];

        t_120[k] = -2.0 * ig_60[k]
                   + f_0 * lg_195[k];

        t_121[k] = -2.0 * ig_61[k]
                   + f_0 * lg_196[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, ig_62, ig_63, ig_64, ig_65, ig_66, \
                         lg_197, lg_198, lg_199, lg_200, lg_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * ig_62[k]
                   + f_0 * lg_197[k];

        t_123[k] = -2.0 * ig_63[k]
                   + f_0 * lg_198[k];

        t_124[k] = -2.0 * ig_64[k]
                   + f_0 * lg_199[k];

        t_125[k] = -2.0 * ig_65[k]
                   + f_0 * lg_200[k];

        t_126[k] = -2.0 * ig_66[k]
                   + f_0 * lg_201[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, ig_67, ig_68, ig_69, ig_70, ig_71, \
                         lg_202, lg_203, lg_204, lg_205, lg_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * ig_67[k]
                   + f_0 * lg_202[k];

        t_128[k] = -2.0 * ig_68[k]
                   + f_0 * lg_203[k];

        t_129[k] = -2.0 * ig_69[k]
                   + f_0 * lg_204[k];

        t_130[k] = -2.0 * ig_70[k]
                   + f_0 * lg_205[k];

        t_131[k] = -2.0 * ig_71[k]
                   + f_0 * lg_206[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, ig_72, ig_73, ig_74, ig_75, ig_76, \
                         lg_207, lg_208, lg_209, lg_210, lg_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * ig_72[k]
                   + f_0 * lg_207[k];

        t_133[k] = -2.0 * ig_73[k]
                   + f_0 * lg_208[k];

        t_134[k] = -2.0 * ig_74[k]
                   + f_0 * lg_209[k];

        t_135[k] = -3.0 * ig_75[k]
                   + f_0 * lg_210[k];

        t_136[k] = -3.0 * ig_76[k]
                   + f_0 * lg_211[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, ig_77, ig_78, ig_79, ig_80, ig_81, \
                         lg_212, lg_213, lg_214, lg_215, lg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -3.0 * ig_77[k]
                   + f_0 * lg_212[k];

        t_138[k] = -3.0 * ig_78[k]
                   + f_0 * lg_213[k];

        t_139[k] = -3.0 * ig_79[k]
                   + f_0 * lg_214[k];

        t_140[k] = -3.0 * ig_80[k]
                   + f_0 * lg_215[k];

        t_141[k] = -3.0 * ig_81[k]
                   + f_0 * lg_216[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, ig_82, ig_83, ig_84, ig_85, ig_86, \
                         lg_217, lg_218, lg_219, lg_220, lg_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -3.0 * ig_82[k]
                   + f_0 * lg_217[k];

        t_143[k] = -3.0 * ig_83[k]
                   + f_0 * lg_218[k];

        t_144[k] = -3.0 * ig_84[k]
                   + f_0 * lg_219[k];

        t_145[k] = -3.0 * ig_85[k]
                   + f_0 * lg_220[k];

        t_146[k] = -3.0 * ig_86[k]
                   + f_0 * lg_221[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, t_152, ig_87, ig_88, ig_89, \
                         lg_222, lg_223, lg_224, lg_240, lg_241, \
                         lg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -3.0 * ig_87[k]
                   + f_0 * lg_222[k];

        t_148[k] = -3.0 * ig_88[k]
                   + f_0 * lg_223[k];

        t_149[k] = -3.0 * ig_89[k]
                   + f_0 * lg_224[k];

        t_150[k] = f_0 * lg_240[k];

        t_151[k] = f_0 * lg_241[k];

        t_152[k] = f_0 * lg_242[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, t_158, t_159, t_160, lg_243, \
                         lg_244, lg_245, lg_246, lg_247, lg_248, lg_249, \
                         lg_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_0 * lg_243[k];

        t_154[k] = f_0 * lg_244[k];

        t_155[k] = f_0 * lg_245[k];

        t_156[k] = f_0 * lg_246[k];

        t_157[k] = f_0 * lg_247[k];

        t_158[k] = f_0 * lg_248[k];

        t_159[k] = f_0 * lg_249[k];

        t_160[k] = f_0 * lg_250[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, t_166, ig_90, ig_91, lg_251, \
                         lg_252, lg_253, lg_254, lg_255, lg_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_0 * lg_251[k];

        t_162[k] = f_0 * lg_252[k];

        t_163[k] = f_0 * lg_253[k];

        t_164[k] = f_0 * lg_254[k];

        t_165[k] = -ig_90[k]
                   + f_0 * lg_255[k];

        t_166[k] = -ig_91[k]
                   + f_0 * lg_256[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, ig_92, ig_93, ig_94, ig_95, ig_96, \
                         lg_257, lg_258, lg_259, lg_260, lg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -ig_92[k]
                   + f_0 * lg_257[k];

        t_168[k] = -ig_93[k]
                   + f_0 * lg_258[k];

        t_169[k] = -ig_94[k]
                   + f_0 * lg_259[k];

        t_170[k] = -ig_95[k]
                   + f_0 * lg_260[k];

        t_171[k] = -ig_96[k]
                   + f_0 * lg_261[k];
    }
}

static auto
compute_prim_geom_10_kg_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ig, const size_t lg,
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

    const auto *lg_262 = buffer.data(lg + 262);
    const auto *lg_263 = buffer.data(lg + 263);
    const auto *lg_264 = buffer.data(lg + 264);
    const auto *lg_265 = buffer.data(lg + 265);
    const auto *lg_266 = buffer.data(lg + 266);
    const auto *lg_267 = buffer.data(lg + 267);
    const auto *lg_268 = buffer.data(lg + 268);
    const auto *lg_269 = buffer.data(lg + 269);
    const auto *lg_270 = buffer.data(lg + 270);
    const auto *lg_271 = buffer.data(lg + 271);
    const auto *lg_272 = buffer.data(lg + 272);
    const auto *lg_273 = buffer.data(lg + 273);
    const auto *lg_274 = buffer.data(lg + 274);
    const auto *lg_275 = buffer.data(lg + 275);
    const auto *lg_276 = buffer.data(lg + 276);
    const auto *lg_277 = buffer.data(lg + 277);
    const auto *lg_278 = buffer.data(lg + 278);
    const auto *lg_279 = buffer.data(lg + 279);
    const auto *lg_280 = buffer.data(lg + 280);
    const auto *lg_281 = buffer.data(lg + 281);
    const auto *lg_282 = buffer.data(lg + 282);
    const auto *lg_283 = buffer.data(lg + 283);
    const auto *lg_284 = buffer.data(lg + 284);
    const auto *lg_285 = buffer.data(lg + 285);
    const auto *lg_286 = buffer.data(lg + 286);
    const auto *lg_287 = buffer.data(lg + 287);
    const auto *lg_288 = buffer.data(lg + 288);
    const auto *lg_289 = buffer.data(lg + 289);
    const auto *lg_290 = buffer.data(lg + 290);
    const auto *lg_291 = buffer.data(lg + 291);
    const auto *lg_292 = buffer.data(lg + 292);
    const auto *lg_293 = buffer.data(lg + 293);
    const auto *lg_294 = buffer.data(lg + 294);
    const auto *lg_295 = buffer.data(lg + 295);
    const auto *lg_296 = buffer.data(lg + 296);
    const auto *lg_297 = buffer.data(lg + 297);
    const auto *lg_298 = buffer.data(lg + 298);
    const auto *lg_299 = buffer.data(lg + 299);
    const auto *lg_300 = buffer.data(lg + 300);
    const auto *lg_301 = buffer.data(lg + 301);
    const auto *lg_302 = buffer.data(lg + 302);
    const auto *lg_303 = buffer.data(lg + 303);
    const auto *lg_304 = buffer.data(lg + 304);
    const auto *lg_305 = buffer.data(lg + 305);
    const auto *lg_306 = buffer.data(lg + 306);
    const auto *lg_307 = buffer.data(lg + 307);
    const auto *lg_308 = buffer.data(lg + 308);
    const auto *lg_309 = buffer.data(lg + 309);
    const auto *lg_310 = buffer.data(lg + 310);
    const auto *lg_311 = buffer.data(lg + 311);
    const auto *lg_312 = buffer.data(lg + 312);
    const auto *lg_313 = buffer.data(lg + 313);
    const auto *lg_314 = buffer.data(lg + 314);
    const auto *lg_330 = buffer.data(lg + 330);
    const auto *lg_331 = buffer.data(lg + 331);
    const auto *lg_332 = buffer.data(lg + 332);
    const auto *lg_333 = buffer.data(lg + 333);
    const auto *lg_334 = buffer.data(lg + 334);
    const auto *lg_335 = buffer.data(lg + 335);
    const auto *lg_336 = buffer.data(lg + 336);
    const auto *lg_337 = buffer.data(lg + 337);
    const auto *lg_338 = buffer.data(lg + 338);
    const auto *lg_339 = buffer.data(lg + 339);
    const auto *lg_340 = buffer.data(lg + 340);
    const auto *lg_341 = buffer.data(lg + 341);
    const auto *lg_342 = buffer.data(lg + 342);
    const auto *lg_343 = buffer.data(lg + 343);
    const auto *lg_344 = buffer.data(lg + 344);
    const auto *lg_345 = buffer.data(lg + 345);
    const auto *lg_346 = buffer.data(lg + 346);
    const auto *lg_347 = buffer.data(lg + 347);
    const auto *lg_348 = buffer.data(lg + 348);
    const auto *lg_349 = buffer.data(lg + 349);
    const auto *lg_350 = buffer.data(lg + 350);
    const auto *lg_351 = buffer.data(lg + 351);
    const auto *lg_352 = buffer.data(lg + 352);
    const auto *lg_353 = buffer.data(lg + 353);
    const auto *lg_354 = buffer.data(lg + 354);
    const auto *lg_355 = buffer.data(lg + 355);
    const auto *lg_356 = buffer.data(lg + 356);
    const auto *lg_357 = buffer.data(lg + 357);
    const auto *lg_358 = buffer.data(lg + 358);
    const auto *lg_359 = buffer.data(lg + 359);
    const auto *lg_360 = buffer.data(lg + 360);
    const auto *lg_361 = buffer.data(lg + 361);
    const auto *lg_362 = buffer.data(lg + 362);
    const auto *lg_363 = buffer.data(lg + 363);
    const auto *lg_364 = buffer.data(lg + 364);
    const auto *lg_365 = buffer.data(lg + 365);
    const auto *lg_366 = buffer.data(lg + 366);
    const auto *lg_367 = buffer.data(lg + 367);
    const auto *lg_368 = buffer.data(lg + 368);
    const auto *lg_369 = buffer.data(lg + 369);
    const auto *lg_370 = buffer.data(lg + 370);
    const auto *lg_371 = buffer.data(lg + 371);
    const auto *lg_372 = buffer.data(lg + 372);
    const auto *lg_373 = buffer.data(lg + 373);
    const auto *lg_374 = buffer.data(lg + 374);
    const auto *lg_375 = buffer.data(lg + 375);
    const auto *lg_376 = buffer.data(lg + 376);
    const auto *lg_377 = buffer.data(lg + 377);
    const auto *lg_378 = buffer.data(lg + 378);
    const auto *lg_379 = buffer.data(lg + 379);
    const auto *lg_380 = buffer.data(lg + 380);
    const auto *lg_381 = buffer.data(lg + 381);
    const auto *lg_382 = buffer.data(lg + 382);
    const auto *lg_383 = buffer.data(lg + 383);
    const auto *lg_384 = buffer.data(lg + 384);
    const auto *lg_385 = buffer.data(lg + 385);
    const auto *lg_386 = buffer.data(lg + 386);
    const auto *lg_387 = buffer.data(lg + 387);
    const auto *lg_388 = buffer.data(lg + 388);
    const auto *lg_389 = buffer.data(lg + 389);
    const auto *lg_390 = buffer.data(lg + 390);
    const auto *lg_391 = buffer.data(lg + 391);
    const auto *lg_392 = buffer.data(lg + 392);
    const auto *lg_393 = buffer.data(lg + 393);
    const auto *lg_394 = buffer.data(lg + 394);
    const auto *lg_395 = buffer.data(lg + 395);
    const auto *lg_396 = buffer.data(lg + 396);
    const auto *lg_397 = buffer.data(lg + 397);
    const auto *lg_398 = buffer.data(lg + 398);
    const auto *lg_399 = buffer.data(lg + 399);
    const auto *lg_400 = buffer.data(lg + 400);
    const auto *lg_401 = buffer.data(lg + 401);
    const auto *lg_402 = buffer.data(lg + 402);
    const auto *lg_403 = buffer.data(lg + 403);
    const auto *lg_404 = buffer.data(lg + 404);
    const auto *lg_405 = buffer.data(lg + 405);
    const auto *lg_406 = buffer.data(lg + 406);
    const auto *lg_407 = buffer.data(lg + 407);
    const auto *lg_408 = buffer.data(lg + 408);
    const auto *lg_409 = buffer.data(lg + 409);
    const auto *lg_410 = buffer.data(lg + 410);
    const auto *lg_411 = buffer.data(lg + 411);
    const auto *lg_412 = buffer.data(lg + 412);
    const auto *lg_413 = buffer.data(lg + 413);
    const auto *lg_414 = buffer.data(lg + 414);
    const auto *lg_415 = buffer.data(lg + 415);
    const auto *lg_416 = buffer.data(lg + 416);
    const auto *lg_417 = buffer.data(lg + 417);
    const auto *lg_418 = buffer.data(lg + 418);
    const auto *lg_419 = buffer.data(lg + 419);
    const auto *lg_435 = buffer.data(lg + 435);
    const auto *lg_436 = buffer.data(lg + 436);
    const auto *lg_437 = buffer.data(lg + 437);
    const auto *lg_438 = buffer.data(lg + 438);
    const auto *lg_439 = buffer.data(lg + 439);
    const auto *lg_440 = buffer.data(lg + 440);
    const auto *lg_441 = buffer.data(lg + 441);
    const auto *lg_442 = buffer.data(lg + 442);
    const auto *lg_443 = buffer.data(lg + 443);
    const auto *lg_444 = buffer.data(lg + 444);
    const auto *lg_445 = buffer.data(lg + 445);
    const auto *lg_446 = buffer.data(lg + 446);
    const auto *lg_447 = buffer.data(lg + 447);
    const auto *lg_448 = buffer.data(lg + 448);
    const auto *lg_449 = buffer.data(lg + 449);
    const auto *lg_450 = buffer.data(lg + 450);
    const auto *lg_451 = buffer.data(lg + 451);

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, ig_97, ig_98, ig_99, ig_100, \
                         ig_101, lg_262, lg_263, lg_264, lg_265, \
                         lg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -ig_97[k]
                   + f_0 * lg_262[k];

        t_173[k] = -ig_98[k]
                   + f_0 * lg_263[k];

        t_174[k] = -ig_99[k]
                   + f_0 * lg_264[k];

        t_175[k] = -ig_100[k]
                   + f_0 * lg_265[k];

        t_176[k] = -ig_101[k]
                   + f_0 * lg_266[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, ig_102, ig_103, ig_104, ig_105, \
                         ig_106, lg_267, lg_268, lg_269, lg_270, \
                         lg_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -ig_102[k]
                   + f_0 * lg_267[k];

        t_178[k] = -ig_103[k]
                   + f_0 * lg_268[k];

        t_179[k] = -ig_104[k]
                   + f_0 * lg_269[k];

        t_180[k] = -2.0 * ig_105[k]
                   + f_0 * lg_270[k];

        t_181[k] = -2.0 * ig_106[k]
                   + f_0 * lg_271[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, ig_107, ig_108, ig_109, ig_110, \
                         ig_111, lg_272, lg_273, lg_274, lg_275, \
                         lg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * ig_107[k]
                   + f_0 * lg_272[k];

        t_183[k] = -2.0 * ig_108[k]
                   + f_0 * lg_273[k];

        t_184[k] = -2.0 * ig_109[k]
                   + f_0 * lg_274[k];

        t_185[k] = -2.0 * ig_110[k]
                   + f_0 * lg_275[k];

        t_186[k] = -2.0 * ig_111[k]
                   + f_0 * lg_276[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, ig_112, ig_113, ig_114, ig_115, \
                         ig_116, lg_277, lg_278, lg_279, lg_280, \
                         lg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * ig_112[k]
                   + f_0 * lg_277[k];

        t_188[k] = -2.0 * ig_113[k]
                   + f_0 * lg_278[k];

        t_189[k] = -2.0 * ig_114[k]
                   + f_0 * lg_279[k];

        t_190[k] = -2.0 * ig_115[k]
                   + f_0 * lg_280[k];

        t_191[k] = -2.0 * ig_116[k]
                   + f_0 * lg_281[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, ig_117, ig_118, ig_119, ig_120, \
                         ig_121, lg_282, lg_283, lg_284, lg_285, \
                         lg_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -2.0 * ig_117[k]
                   + f_0 * lg_282[k];

        t_193[k] = -2.0 * ig_118[k]
                   + f_0 * lg_283[k];

        t_194[k] = -2.0 * ig_119[k]
                   + f_0 * lg_284[k];

        t_195[k] = -3.0 * ig_120[k]
                   + f_0 * lg_285[k];

        t_196[k] = -3.0 * ig_121[k]
                   + f_0 * lg_286[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, ig_122, ig_123, ig_124, ig_125, \
                         ig_126, lg_287, lg_288, lg_289, lg_290, \
                         lg_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -3.0 * ig_122[k]
                   + f_0 * lg_287[k];

        t_198[k] = -3.0 * ig_123[k]
                   + f_0 * lg_288[k];

        t_199[k] = -3.0 * ig_124[k]
                   + f_0 * lg_289[k];

        t_200[k] = -3.0 * ig_125[k]
                   + f_0 * lg_290[k];

        t_201[k] = -3.0 * ig_126[k]
                   + f_0 * lg_291[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, ig_127, ig_128, ig_129, ig_130, \
                         ig_131, lg_292, lg_293, lg_294, lg_295, \
                         lg_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -3.0 * ig_127[k]
                   + f_0 * lg_292[k];

        t_203[k] = -3.0 * ig_128[k]
                   + f_0 * lg_293[k];

        t_204[k] = -3.0 * ig_129[k]
                   + f_0 * lg_294[k];

        t_205[k] = -3.0 * ig_130[k]
                   + f_0 * lg_295[k];

        t_206[k] = -3.0 * ig_131[k]
                   + f_0 * lg_296[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, ig_132, ig_133, ig_134, ig_135, \
                         ig_136, lg_297, lg_298, lg_299, lg_300, \
                         lg_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -3.0 * ig_132[k]
                   + f_0 * lg_297[k];

        t_208[k] = -3.0 * ig_133[k]
                   + f_0 * lg_298[k];

        t_209[k] = -3.0 * ig_134[k]
                   + f_0 * lg_299[k];

        t_210[k] = -4.0 * ig_135[k]
                   + f_0 * lg_300[k];

        t_211[k] = -4.0 * ig_136[k]
                   + f_0 * lg_301[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, ig_137, ig_138, ig_139, ig_140, \
                         ig_141, lg_302, lg_303, lg_304, lg_305, \
                         lg_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = -4.0 * ig_137[k]
                   + f_0 * lg_302[k];

        t_213[k] = -4.0 * ig_138[k]
                   + f_0 * lg_303[k];

        t_214[k] = -4.0 * ig_139[k]
                   + f_0 * lg_304[k];

        t_215[k] = -4.0 * ig_140[k]
                   + f_0 * lg_305[k];

        t_216[k] = -4.0 * ig_141[k]
                   + f_0 * lg_306[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, ig_142, ig_143, ig_144, ig_145, \
                         ig_146, lg_307, lg_308, lg_309, lg_310, \
                         lg_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = -4.0 * ig_142[k]
                   + f_0 * lg_307[k];

        t_218[k] = -4.0 * ig_143[k]
                   + f_0 * lg_308[k];

        t_219[k] = -4.0 * ig_144[k]
                   + f_0 * lg_309[k];

        t_220[k] = -4.0 * ig_145[k]
                   + f_0 * lg_310[k];

        t_221[k] = -4.0 * ig_146[k]
                   + f_0 * lg_311[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, ig_147, ig_148, ig_149, \
                         lg_312, lg_313, lg_314, lg_330, lg_331, \
                         lg_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -4.0 * ig_147[k]
                   + f_0 * lg_312[k];

        t_223[k] = -4.0 * ig_148[k]
                   + f_0 * lg_313[k];

        t_224[k] = -4.0 * ig_149[k]
                   + f_0 * lg_314[k];

        t_225[k] = f_0 * lg_330[k];

        t_226[k] = f_0 * lg_331[k];

        t_227[k] = f_0 * lg_332[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, t_233, t_234, t_235, lg_333, \
                         lg_334, lg_335, lg_336, lg_337, lg_338, lg_339, \
                         lg_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_0 * lg_333[k];

        t_229[k] = f_0 * lg_334[k];

        t_230[k] = f_0 * lg_335[k];

        t_231[k] = f_0 * lg_336[k];

        t_232[k] = f_0 * lg_337[k];

        t_233[k] = f_0 * lg_338[k];

        t_234[k] = f_0 * lg_339[k];

        t_235[k] = f_0 * lg_340[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, t_241, ig_150, ig_151, lg_341, \
                         lg_342, lg_343, lg_344, lg_345, lg_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_0 * lg_341[k];

        t_237[k] = f_0 * lg_342[k];

        t_238[k] = f_0 * lg_343[k];

        t_239[k] = f_0 * lg_344[k];

        t_240[k] = -ig_150[k]
                   + f_0 * lg_345[k];

        t_241[k] = -ig_151[k]
                   + f_0 * lg_346[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, ig_152, ig_153, ig_154, ig_155, \
                         ig_156, lg_347, lg_348, lg_349, lg_350, \
                         lg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = -ig_152[k]
                   + f_0 * lg_347[k];

        t_243[k] = -ig_153[k]
                   + f_0 * lg_348[k];

        t_244[k] = -ig_154[k]
                   + f_0 * lg_349[k];

        t_245[k] = -ig_155[k]
                   + f_0 * lg_350[k];

        t_246[k] = -ig_156[k]
                   + f_0 * lg_351[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, ig_157, ig_158, ig_159, ig_160, \
                         ig_161, lg_352, lg_353, lg_354, lg_355, \
                         lg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = -ig_157[k]
                   + f_0 * lg_352[k];

        t_248[k] = -ig_158[k]
                   + f_0 * lg_353[k];

        t_249[k] = -ig_159[k]
                   + f_0 * lg_354[k];

        t_250[k] = -ig_160[k]
                   + f_0 * lg_355[k];

        t_251[k] = -ig_161[k]
                   + f_0 * lg_356[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, ig_162, ig_163, ig_164, ig_165, \
                         ig_166, lg_357, lg_358, lg_359, lg_360, \
                         lg_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -ig_162[k]
                   + f_0 * lg_357[k];

        t_253[k] = -ig_163[k]
                   + f_0 * lg_358[k];

        t_254[k] = -ig_164[k]
                   + f_0 * lg_359[k];

        t_255[k] = -2.0 * ig_165[k]
                   + f_0 * lg_360[k];

        t_256[k] = -2.0 * ig_166[k]
                   + f_0 * lg_361[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, ig_167, ig_168, ig_169, ig_170, \
                         ig_171, lg_362, lg_363, lg_364, lg_365, \
                         lg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = -2.0 * ig_167[k]
                   + f_0 * lg_362[k];

        t_258[k] = -2.0 * ig_168[k]
                   + f_0 * lg_363[k];

        t_259[k] = -2.0 * ig_169[k]
                   + f_0 * lg_364[k];

        t_260[k] = -2.0 * ig_170[k]
                   + f_0 * lg_365[k];

        t_261[k] = -2.0 * ig_171[k]
                   + f_0 * lg_366[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, ig_172, ig_173, ig_174, ig_175, \
                         ig_176, lg_367, lg_368, lg_369, lg_370, \
                         lg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = -2.0 * ig_172[k]
                   + f_0 * lg_367[k];

        t_263[k] = -2.0 * ig_173[k]
                   + f_0 * lg_368[k];

        t_264[k] = -2.0 * ig_174[k]
                   + f_0 * lg_369[k];

        t_265[k] = -2.0 * ig_175[k]
                   + f_0 * lg_370[k];

        t_266[k] = -2.0 * ig_176[k]
                   + f_0 * lg_371[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, ig_177, ig_178, ig_179, ig_180, \
                         ig_181, lg_372, lg_373, lg_374, lg_375, \
                         lg_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = -2.0 * ig_177[k]
                   + f_0 * lg_372[k];

        t_268[k] = -2.0 * ig_178[k]
                   + f_0 * lg_373[k];

        t_269[k] = -2.0 * ig_179[k]
                   + f_0 * lg_374[k];

        t_270[k] = -3.0 * ig_180[k]
                   + f_0 * lg_375[k];

        t_271[k] = -3.0 * ig_181[k]
                   + f_0 * lg_376[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, ig_182, ig_183, ig_184, ig_185, \
                         ig_186, lg_377, lg_378, lg_379, lg_380, \
                         lg_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -3.0 * ig_182[k]
                   + f_0 * lg_377[k];

        t_273[k] = -3.0 * ig_183[k]
                   + f_0 * lg_378[k];

        t_274[k] = -3.0 * ig_184[k]
                   + f_0 * lg_379[k];

        t_275[k] = -3.0 * ig_185[k]
                   + f_0 * lg_380[k];

        t_276[k] = -3.0 * ig_186[k]
                   + f_0 * lg_381[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, t_281, ig_187, ig_188, ig_189, ig_190, \
                         ig_191, lg_382, lg_383, lg_384, lg_385, \
                         lg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = -3.0 * ig_187[k]
                   + f_0 * lg_382[k];

        t_278[k] = -3.0 * ig_188[k]
                   + f_0 * lg_383[k];

        t_279[k] = -3.0 * ig_189[k]
                   + f_0 * lg_384[k];

        t_280[k] = -3.0 * ig_190[k]
                   + f_0 * lg_385[k];

        t_281[k] = -3.0 * ig_191[k]
                   + f_0 * lg_386[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, ig_192, ig_193, ig_194, ig_195, \
                         ig_196, lg_387, lg_388, lg_389, lg_390, \
                         lg_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = -3.0 * ig_192[k]
                   + f_0 * lg_387[k];

        t_283[k] = -3.0 * ig_193[k]
                   + f_0 * lg_388[k];

        t_284[k] = -3.0 * ig_194[k]
                   + f_0 * lg_389[k];

        t_285[k] = -4.0 * ig_195[k]
                   + f_0 * lg_390[k];

        t_286[k] = -4.0 * ig_196[k]
                   + f_0 * lg_391[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, ig_197, ig_198, ig_199, ig_200, \
                         ig_201, lg_392, lg_393, lg_394, lg_395, \
                         lg_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = -4.0 * ig_197[k]
                   + f_0 * lg_392[k];

        t_288[k] = -4.0 * ig_198[k]
                   + f_0 * lg_393[k];

        t_289[k] = -4.0 * ig_199[k]
                   + f_0 * lg_394[k];

        t_290[k] = -4.0 * ig_200[k]
                   + f_0 * lg_395[k];

        t_291[k] = -4.0 * ig_201[k]
                   + f_0 * lg_396[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, ig_202, ig_203, ig_204, ig_205, \
                         ig_206, lg_397, lg_398, lg_399, lg_400, \
                         lg_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = -4.0 * ig_202[k]
                   + f_0 * lg_397[k];

        t_293[k] = -4.0 * ig_203[k]
                   + f_0 * lg_398[k];

        t_294[k] = -4.0 * ig_204[k]
                   + f_0 * lg_399[k];

        t_295[k] = -4.0 * ig_205[k]
                   + f_0 * lg_400[k];

        t_296[k] = -4.0 * ig_206[k]
                   + f_0 * lg_401[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, ig_207, ig_208, ig_209, ig_210, \
                         ig_211, lg_402, lg_403, lg_404, lg_405, \
                         lg_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = -4.0 * ig_207[k]
                   + f_0 * lg_402[k];

        t_298[k] = -4.0 * ig_208[k]
                   + f_0 * lg_403[k];

        t_299[k] = -4.0 * ig_209[k]
                   + f_0 * lg_404[k];

        t_300[k] = -5.0 * ig_210[k]
                   + f_0 * lg_405[k];

        t_301[k] = -5.0 * ig_211[k]
                   + f_0 * lg_406[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, ig_212, ig_213, ig_214, ig_215, \
                         ig_216, lg_407, lg_408, lg_409, lg_410, \
                         lg_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = -5.0 * ig_212[k]
                   + f_0 * lg_407[k];

        t_303[k] = -5.0 * ig_213[k]
                   + f_0 * lg_408[k];

        t_304[k] = -5.0 * ig_214[k]
                   + f_0 * lg_409[k];

        t_305[k] = -5.0 * ig_215[k]
                   + f_0 * lg_410[k];

        t_306[k] = -5.0 * ig_216[k]
                   + f_0 * lg_411[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, ig_217, ig_218, ig_219, ig_220, \
                         ig_221, lg_412, lg_413, lg_414, lg_415, \
                         lg_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = -5.0 * ig_217[k]
                   + f_0 * lg_412[k];

        t_308[k] = -5.0 * ig_218[k]
                   + f_0 * lg_413[k];

        t_309[k] = -5.0 * ig_219[k]
                   + f_0 * lg_414[k];

        t_310[k] = -5.0 * ig_220[k]
                   + f_0 * lg_415[k];

        t_311[k] = -5.0 * ig_221[k]
                   + f_0 * lg_416[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, t_317, ig_222, ig_223, ig_224, \
                         lg_417, lg_418, lg_419, lg_435, lg_436, \
                         lg_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -5.0 * ig_222[k]
                   + f_0 * lg_417[k];

        t_313[k] = -5.0 * ig_223[k]
                   + f_0 * lg_418[k];

        t_314[k] = -5.0 * ig_224[k]
                   + f_0 * lg_419[k];

        t_315[k] = f_0 * lg_435[k];

        t_316[k] = f_0 * lg_436[k];

        t_317[k] = f_0 * lg_437[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, t_323, t_324, t_325, lg_438, \
                         lg_439, lg_440, lg_441, lg_442, lg_443, lg_444, \
                         lg_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_0 * lg_438[k];

        t_319[k] = f_0 * lg_439[k];

        t_320[k] = f_0 * lg_440[k];

        t_321[k] = f_0 * lg_441[k];

        t_322[k] = f_0 * lg_442[k];

        t_323[k] = f_0 * lg_443[k];

        t_324[k] = f_0 * lg_444[k];

        t_325[k] = f_0 * lg_445[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, ig_225, ig_226, lg_446, \
                         lg_447, lg_448, lg_449, lg_450, lg_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_0 * lg_446[k];

        t_327[k] = f_0 * lg_447[k];

        t_328[k] = f_0 * lg_448[k];

        t_329[k] = f_0 * lg_449[k];

        t_330[k] = -ig_225[k]
                   + f_0 * lg_450[k];

        t_331[k] = -ig_226[k]
                   + f_0 * lg_451[k];
    }
}

static auto
compute_prim_geom_10_kg_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ig, const size_t lg,
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

    const auto *lg_452 = buffer.data(lg + 452);
    const auto *lg_453 = buffer.data(lg + 453);
    const auto *lg_454 = buffer.data(lg + 454);
    const auto *lg_455 = buffer.data(lg + 455);
    const auto *lg_456 = buffer.data(lg + 456);
    const auto *lg_457 = buffer.data(lg + 457);
    const auto *lg_458 = buffer.data(lg + 458);
    const auto *lg_459 = buffer.data(lg + 459);
    const auto *lg_460 = buffer.data(lg + 460);
    const auto *lg_461 = buffer.data(lg + 461);
    const auto *lg_462 = buffer.data(lg + 462);
    const auto *lg_463 = buffer.data(lg + 463);
    const auto *lg_464 = buffer.data(lg + 464);
    const auto *lg_465 = buffer.data(lg + 465);
    const auto *lg_466 = buffer.data(lg + 466);
    const auto *lg_467 = buffer.data(lg + 467);
    const auto *lg_468 = buffer.data(lg + 468);
    const auto *lg_469 = buffer.data(lg + 469);
    const auto *lg_470 = buffer.data(lg + 470);
    const auto *lg_471 = buffer.data(lg + 471);
    const auto *lg_472 = buffer.data(lg + 472);
    const auto *lg_473 = buffer.data(lg + 473);
    const auto *lg_474 = buffer.data(lg + 474);
    const auto *lg_475 = buffer.data(lg + 475);
    const auto *lg_476 = buffer.data(lg + 476);
    const auto *lg_477 = buffer.data(lg + 477);
    const auto *lg_478 = buffer.data(lg + 478);
    const auto *lg_479 = buffer.data(lg + 479);
    const auto *lg_480 = buffer.data(lg + 480);
    const auto *lg_481 = buffer.data(lg + 481);
    const auto *lg_482 = buffer.data(lg + 482);
    const auto *lg_483 = buffer.data(lg + 483);
    const auto *lg_484 = buffer.data(lg + 484);
    const auto *lg_485 = buffer.data(lg + 485);
    const auto *lg_486 = buffer.data(lg + 486);
    const auto *lg_487 = buffer.data(lg + 487);
    const auto *lg_488 = buffer.data(lg + 488);
    const auto *lg_489 = buffer.data(lg + 489);
    const auto *lg_490 = buffer.data(lg + 490);
    const auto *lg_491 = buffer.data(lg + 491);
    const auto *lg_492 = buffer.data(lg + 492);
    const auto *lg_493 = buffer.data(lg + 493);
    const auto *lg_494 = buffer.data(lg + 494);
    const auto *lg_495 = buffer.data(lg + 495);
    const auto *lg_496 = buffer.data(lg + 496);
    const auto *lg_497 = buffer.data(lg + 497);
    const auto *lg_498 = buffer.data(lg + 498);
    const auto *lg_499 = buffer.data(lg + 499);
    const auto *lg_500 = buffer.data(lg + 500);
    const auto *lg_501 = buffer.data(lg + 501);
    const auto *lg_502 = buffer.data(lg + 502);
    const auto *lg_503 = buffer.data(lg + 503);
    const auto *lg_504 = buffer.data(lg + 504);
    const auto *lg_505 = buffer.data(lg + 505);
    const auto *lg_506 = buffer.data(lg + 506);
    const auto *lg_507 = buffer.data(lg + 507);
    const auto *lg_508 = buffer.data(lg + 508);
    const auto *lg_509 = buffer.data(lg + 509);
    const auto *lg_510 = buffer.data(lg + 510);
    const auto *lg_511 = buffer.data(lg + 511);
    const auto *lg_512 = buffer.data(lg + 512);
    const auto *lg_513 = buffer.data(lg + 513);
    const auto *lg_514 = buffer.data(lg + 514);
    const auto *lg_515 = buffer.data(lg + 515);
    const auto *lg_516 = buffer.data(lg + 516);
    const auto *lg_517 = buffer.data(lg + 517);
    const auto *lg_518 = buffer.data(lg + 518);
    const auto *lg_519 = buffer.data(lg + 519);
    const auto *lg_520 = buffer.data(lg + 520);
    const auto *lg_521 = buffer.data(lg + 521);
    const auto *lg_522 = buffer.data(lg + 522);
    const auto *lg_523 = buffer.data(lg + 523);
    const auto *lg_524 = buffer.data(lg + 524);
    const auto *lg_525 = buffer.data(lg + 525);
    const auto *lg_526 = buffer.data(lg + 526);
    const auto *lg_527 = buffer.data(lg + 527);
    const auto *lg_528 = buffer.data(lg + 528);
    const auto *lg_529 = buffer.data(lg + 529);
    const auto *lg_530 = buffer.data(lg + 530);
    const auto *lg_531 = buffer.data(lg + 531);
    const auto *lg_532 = buffer.data(lg + 532);
    const auto *lg_533 = buffer.data(lg + 533);
    const auto *lg_534 = buffer.data(lg + 534);
    const auto *lg_535 = buffer.data(lg + 535);
    const auto *lg_536 = buffer.data(lg + 536);
    const auto *lg_537 = buffer.data(lg + 537);
    const auto *lg_538 = buffer.data(lg + 538);
    const auto *lg_539 = buffer.data(lg + 539);
    const auto *lg_555 = buffer.data(lg + 555);
    const auto *lg_556 = buffer.data(lg + 556);
    const auto *lg_557 = buffer.data(lg + 557);
    const auto *lg_558 = buffer.data(lg + 558);
    const auto *lg_559 = buffer.data(lg + 559);
    const auto *lg_560 = buffer.data(lg + 560);
    const auto *lg_561 = buffer.data(lg + 561);
    const auto *lg_562 = buffer.data(lg + 562);
    const auto *lg_563 = buffer.data(lg + 563);
    const auto *lg_564 = buffer.data(lg + 564);
    const auto *lg_565 = buffer.data(lg + 565);
    const auto *lg_566 = buffer.data(lg + 566);
    const auto *lg_567 = buffer.data(lg + 567);
    const auto *lg_568 = buffer.data(lg + 568);
    const auto *lg_569 = buffer.data(lg + 569);
    const auto *lg_570 = buffer.data(lg + 570);
    const auto *lg_571 = buffer.data(lg + 571);
    const auto *lg_572 = buffer.data(lg + 572);
    const auto *lg_573 = buffer.data(lg + 573);
    const auto *lg_574 = buffer.data(lg + 574);
    const auto *lg_575 = buffer.data(lg + 575);
    const auto *lg_576 = buffer.data(lg + 576);
    const auto *lg_577 = buffer.data(lg + 577);
    const auto *lg_578 = buffer.data(lg + 578);
    const auto *lg_579 = buffer.data(lg + 579);
    const auto *lg_580 = buffer.data(lg + 580);
    const auto *lg_581 = buffer.data(lg + 581);
    const auto *lg_582 = buffer.data(lg + 582);
    const auto *lg_583 = buffer.data(lg + 583);
    const auto *lg_584 = buffer.data(lg + 584);
    const auto *lg_585 = buffer.data(lg + 585);
    const auto *lg_586 = buffer.data(lg + 586);
    const auto *lg_587 = buffer.data(lg + 587);
    const auto *lg_588 = buffer.data(lg + 588);
    const auto *lg_589 = buffer.data(lg + 589);
    const auto *lg_590 = buffer.data(lg + 590);
    const auto *lg_591 = buffer.data(lg + 591);
    const auto *lg_592 = buffer.data(lg + 592);
    const auto *lg_593 = buffer.data(lg + 593);
    const auto *lg_594 = buffer.data(lg + 594);
    const auto *lg_595 = buffer.data(lg + 595);
    const auto *lg_596 = buffer.data(lg + 596);
    const auto *lg_597 = buffer.data(lg + 597);
    const auto *lg_598 = buffer.data(lg + 598);
    const auto *lg_599 = buffer.data(lg + 599);
    const auto *lg_600 = buffer.data(lg + 600);
    const auto *lg_601 = buffer.data(lg + 601);
    const auto *lg_602 = buffer.data(lg + 602);
    const auto *lg_603 = buffer.data(lg + 603);
    const auto *lg_604 = buffer.data(lg + 604);
    const auto *lg_605 = buffer.data(lg + 605);
    const auto *lg_606 = buffer.data(lg + 606);
    const auto *lg_607 = buffer.data(lg + 607);
    const auto *lg_608 = buffer.data(lg + 608);
    const auto *lg_609 = buffer.data(lg + 609);
    const auto *lg_610 = buffer.data(lg + 610);
    const auto *lg_611 = buffer.data(lg + 611);
    const auto *lg_612 = buffer.data(lg + 612);
    const auto *lg_613 = buffer.data(lg + 613);
    const auto *lg_614 = buffer.data(lg + 614);
    const auto *lg_615 = buffer.data(lg + 615);
    const auto *lg_616 = buffer.data(lg + 616);
    const auto *lg_617 = buffer.data(lg + 617);
    const auto *lg_618 = buffer.data(lg + 618);
    const auto *lg_619 = buffer.data(lg + 619);
    const auto *lg_620 = buffer.data(lg + 620);
    const auto *lg_621 = buffer.data(lg + 621);

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, ig_227, ig_228, ig_229, ig_230, \
                         ig_231, lg_452, lg_453, lg_454, lg_455, \
                         lg_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -ig_227[k]
                   + f_0 * lg_452[k];

        t_333[k] = -ig_228[k]
                   + f_0 * lg_453[k];

        t_334[k] = -ig_229[k]
                   + f_0 * lg_454[k];

        t_335[k] = -ig_230[k]
                   + f_0 * lg_455[k];

        t_336[k] = -ig_231[k]
                   + f_0 * lg_456[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, ig_232, ig_233, ig_234, ig_235, \
                         ig_236, lg_457, lg_458, lg_459, lg_460, \
                         lg_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -ig_232[k]
                   + f_0 * lg_457[k];

        t_338[k] = -ig_233[k]
                   + f_0 * lg_458[k];

        t_339[k] = -ig_234[k]
                   + f_0 * lg_459[k];

        t_340[k] = -ig_235[k]
                   + f_0 * lg_460[k];

        t_341[k] = -ig_236[k]
                   + f_0 * lg_461[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, ig_237, ig_238, ig_239, ig_240, \
                         ig_241, lg_462, lg_463, lg_464, lg_465, \
                         lg_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -ig_237[k]
                   + f_0 * lg_462[k];

        t_343[k] = -ig_238[k]
                   + f_0 * lg_463[k];

        t_344[k] = -ig_239[k]
                   + f_0 * lg_464[k];

        t_345[k] = -2.0 * ig_240[k]
                   + f_0 * lg_465[k];

        t_346[k] = -2.0 * ig_241[k]
                   + f_0 * lg_466[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, ig_242, ig_243, ig_244, ig_245, \
                         ig_246, lg_467, lg_468, lg_469, lg_470, \
                         lg_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -2.0 * ig_242[k]
                   + f_0 * lg_467[k];

        t_348[k] = -2.0 * ig_243[k]
                   + f_0 * lg_468[k];

        t_349[k] = -2.0 * ig_244[k]
                   + f_0 * lg_469[k];

        t_350[k] = -2.0 * ig_245[k]
                   + f_0 * lg_470[k];

        t_351[k] = -2.0 * ig_246[k]
                   + f_0 * lg_471[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, ig_247, ig_248, ig_249, ig_250, \
                         ig_251, lg_472, lg_473, lg_474, lg_475, \
                         lg_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -2.0 * ig_247[k]
                   + f_0 * lg_472[k];

        t_353[k] = -2.0 * ig_248[k]
                   + f_0 * lg_473[k];

        t_354[k] = -2.0 * ig_249[k]
                   + f_0 * lg_474[k];

        t_355[k] = -2.0 * ig_250[k]
                   + f_0 * lg_475[k];

        t_356[k] = -2.0 * ig_251[k]
                   + f_0 * lg_476[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, ig_252, ig_253, ig_254, ig_255, \
                         ig_256, lg_477, lg_478, lg_479, lg_480, \
                         lg_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -2.0 * ig_252[k]
                   + f_0 * lg_477[k];

        t_358[k] = -2.0 * ig_253[k]
                   + f_0 * lg_478[k];

        t_359[k] = -2.0 * ig_254[k]
                   + f_0 * lg_479[k];

        t_360[k] = -3.0 * ig_255[k]
                   + f_0 * lg_480[k];

        t_361[k] = -3.0 * ig_256[k]
                   + f_0 * lg_481[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, ig_257, ig_258, ig_259, ig_260, \
                         ig_261, lg_482, lg_483, lg_484, lg_485, \
                         lg_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -3.0 * ig_257[k]
                   + f_0 * lg_482[k];

        t_363[k] = -3.0 * ig_258[k]
                   + f_0 * lg_483[k];

        t_364[k] = -3.0 * ig_259[k]
                   + f_0 * lg_484[k];

        t_365[k] = -3.0 * ig_260[k]
                   + f_0 * lg_485[k];

        t_366[k] = -3.0 * ig_261[k]
                   + f_0 * lg_486[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, ig_262, ig_263, ig_264, ig_265, \
                         ig_266, lg_487, lg_488, lg_489, lg_490, \
                         lg_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -3.0 * ig_262[k]
                   + f_0 * lg_487[k];

        t_368[k] = -3.0 * ig_263[k]
                   + f_0 * lg_488[k];

        t_369[k] = -3.0 * ig_264[k]
                   + f_0 * lg_489[k];

        t_370[k] = -3.0 * ig_265[k]
                   + f_0 * lg_490[k];

        t_371[k] = -3.0 * ig_266[k]
                   + f_0 * lg_491[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, ig_267, ig_268, ig_269, ig_270, \
                         ig_271, lg_492, lg_493, lg_494, lg_495, \
                         lg_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -3.0 * ig_267[k]
                   + f_0 * lg_492[k];

        t_373[k] = -3.0 * ig_268[k]
                   + f_0 * lg_493[k];

        t_374[k] = -3.0 * ig_269[k]
                   + f_0 * lg_494[k];

        t_375[k] = -4.0 * ig_270[k]
                   + f_0 * lg_495[k];

        t_376[k] = -4.0 * ig_271[k]
                   + f_0 * lg_496[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, ig_272, ig_273, ig_274, ig_275, \
                         ig_276, lg_497, lg_498, lg_499, lg_500, \
                         lg_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -4.0 * ig_272[k]
                   + f_0 * lg_497[k];

        t_378[k] = -4.0 * ig_273[k]
                   + f_0 * lg_498[k];

        t_379[k] = -4.0 * ig_274[k]
                   + f_0 * lg_499[k];

        t_380[k] = -4.0 * ig_275[k]
                   + f_0 * lg_500[k];

        t_381[k] = -4.0 * ig_276[k]
                   + f_0 * lg_501[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, ig_277, ig_278, ig_279, ig_280, \
                         ig_281, lg_502, lg_503, lg_504, lg_505, \
                         lg_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -4.0 * ig_277[k]
                   + f_0 * lg_502[k];

        t_383[k] = -4.0 * ig_278[k]
                   + f_0 * lg_503[k];

        t_384[k] = -4.0 * ig_279[k]
                   + f_0 * lg_504[k];

        t_385[k] = -4.0 * ig_280[k]
                   + f_0 * lg_505[k];

        t_386[k] = -4.0 * ig_281[k]
                   + f_0 * lg_506[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, ig_282, ig_283, ig_284, ig_285, \
                         ig_286, lg_507, lg_508, lg_509, lg_510, \
                         lg_511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -4.0 * ig_282[k]
                   + f_0 * lg_507[k];

        t_388[k] = -4.0 * ig_283[k]
                   + f_0 * lg_508[k];

        t_389[k] = -4.0 * ig_284[k]
                   + f_0 * lg_509[k];

        t_390[k] = -5.0 * ig_285[k]
                   + f_0 * lg_510[k];

        t_391[k] = -5.0 * ig_286[k]
                   + f_0 * lg_511[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, ig_287, ig_288, ig_289, ig_290, \
                         ig_291, lg_512, lg_513, lg_514, lg_515, \
                         lg_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -5.0 * ig_287[k]
                   + f_0 * lg_512[k];

        t_393[k] = -5.0 * ig_288[k]
                   + f_0 * lg_513[k];

        t_394[k] = -5.0 * ig_289[k]
                   + f_0 * lg_514[k];

        t_395[k] = -5.0 * ig_290[k]
                   + f_0 * lg_515[k];

        t_396[k] = -5.0 * ig_291[k]
                   + f_0 * lg_516[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, ig_292, ig_293, ig_294, ig_295, \
                         ig_296, lg_517, lg_518, lg_519, lg_520, \
                         lg_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -5.0 * ig_292[k]
                   + f_0 * lg_517[k];

        t_398[k] = -5.0 * ig_293[k]
                   + f_0 * lg_518[k];

        t_399[k] = -5.0 * ig_294[k]
                   + f_0 * lg_519[k];

        t_400[k] = -5.0 * ig_295[k]
                   + f_0 * lg_520[k];

        t_401[k] = -5.0 * ig_296[k]
                   + f_0 * lg_521[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, ig_297, ig_298, ig_299, ig_300, \
                         ig_301, lg_522, lg_523, lg_524, lg_525, \
                         lg_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -5.0 * ig_297[k]
                   + f_0 * lg_522[k];

        t_403[k] = -5.0 * ig_298[k]
                   + f_0 * lg_523[k];

        t_404[k] = -5.0 * ig_299[k]
                   + f_0 * lg_524[k];

        t_405[k] = -6.0 * ig_300[k]
                   + f_0 * lg_525[k];

        t_406[k] = -6.0 * ig_301[k]
                   + f_0 * lg_526[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, ig_302, ig_303, ig_304, ig_305, \
                         ig_306, lg_527, lg_528, lg_529, lg_530, \
                         lg_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -6.0 * ig_302[k]
                   + f_0 * lg_527[k];

        t_408[k] = -6.0 * ig_303[k]
                   + f_0 * lg_528[k];

        t_409[k] = -6.0 * ig_304[k]
                   + f_0 * lg_529[k];

        t_410[k] = -6.0 * ig_305[k]
                   + f_0 * lg_530[k];

        t_411[k] = -6.0 * ig_306[k]
                   + f_0 * lg_531[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, ig_307, ig_308, ig_309, ig_310, \
                         ig_311, lg_532, lg_533, lg_534, lg_535, \
                         lg_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -6.0 * ig_307[k]
                   + f_0 * lg_532[k];

        t_413[k] = -6.0 * ig_308[k]
                   + f_0 * lg_533[k];

        t_414[k] = -6.0 * ig_309[k]
                   + f_0 * lg_534[k];

        t_415[k] = -6.0 * ig_310[k]
                   + f_0 * lg_535[k];

        t_416[k] = -6.0 * ig_311[k]
                   + f_0 * lg_536[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, t_422, ig_312, ig_313, ig_314, \
                         lg_537, lg_538, lg_539, lg_555, lg_556, \
                         lg_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -6.0 * ig_312[k]
                   + f_0 * lg_537[k];

        t_418[k] = -6.0 * ig_313[k]
                   + f_0 * lg_538[k];

        t_419[k] = -6.0 * ig_314[k]
                   + f_0 * lg_539[k];

        t_420[k] = f_0 * lg_555[k];

        t_421[k] = f_0 * lg_556[k];

        t_422[k] = f_0 * lg_557[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, t_428, t_429, t_430, lg_558, \
                         lg_559, lg_560, lg_561, lg_562, lg_563, lg_564, \
                         lg_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_0 * lg_558[k];

        t_424[k] = f_0 * lg_559[k];

        t_425[k] = f_0 * lg_560[k];

        t_426[k] = f_0 * lg_561[k];

        t_427[k] = f_0 * lg_562[k];

        t_428[k] = f_0 * lg_563[k];

        t_429[k] = f_0 * lg_564[k];

        t_430[k] = f_0 * lg_565[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, t_436, ig_315, ig_316, lg_566, \
                         lg_567, lg_568, lg_569, lg_570, lg_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_0 * lg_566[k];

        t_432[k] = f_0 * lg_567[k];

        t_433[k] = f_0 * lg_568[k];

        t_434[k] = f_0 * lg_569[k];

        t_435[k] = -ig_315[k]
                   + f_0 * lg_570[k];

        t_436[k] = -ig_316[k]
                   + f_0 * lg_571[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, t_441, ig_317, ig_318, ig_319, ig_320, \
                         ig_321, lg_572, lg_573, lg_574, lg_575, \
                         lg_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = -ig_317[k]
                   + f_0 * lg_572[k];

        t_438[k] = -ig_318[k]
                   + f_0 * lg_573[k];

        t_439[k] = -ig_319[k]
                   + f_0 * lg_574[k];

        t_440[k] = -ig_320[k]
                   + f_0 * lg_575[k];

        t_441[k] = -ig_321[k]
                   + f_0 * lg_576[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, t_446, ig_322, ig_323, ig_324, ig_325, \
                         ig_326, lg_577, lg_578, lg_579, lg_580, \
                         lg_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = -ig_322[k]
                   + f_0 * lg_577[k];

        t_443[k] = -ig_323[k]
                   + f_0 * lg_578[k];

        t_444[k] = -ig_324[k]
                   + f_0 * lg_579[k];

        t_445[k] = -ig_325[k]
                   + f_0 * lg_580[k];

        t_446[k] = -ig_326[k]
                   + f_0 * lg_581[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, t_451, ig_327, ig_328, ig_329, ig_330, \
                         ig_331, lg_582, lg_583, lg_584, lg_585, \
                         lg_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = -ig_327[k]
                   + f_0 * lg_582[k];

        t_448[k] = -ig_328[k]
                   + f_0 * lg_583[k];

        t_449[k] = -ig_329[k]
                   + f_0 * lg_584[k];

        t_450[k] = -2.0 * ig_330[k]
                   + f_0 * lg_585[k];

        t_451[k] = -2.0 * ig_331[k]
                   + f_0 * lg_586[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, ig_332, ig_333, ig_334, ig_335, \
                         ig_336, lg_587, lg_588, lg_589, lg_590, \
                         lg_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = -2.0 * ig_332[k]
                   + f_0 * lg_587[k];

        t_453[k] = -2.0 * ig_333[k]
                   + f_0 * lg_588[k];

        t_454[k] = -2.0 * ig_334[k]
                   + f_0 * lg_589[k];

        t_455[k] = -2.0 * ig_335[k]
                   + f_0 * lg_590[k];

        t_456[k] = -2.0 * ig_336[k]
                   + f_0 * lg_591[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, t_461, ig_337, ig_338, ig_339, ig_340, \
                         ig_341, lg_592, lg_593, lg_594, lg_595, \
                         lg_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = -2.0 * ig_337[k]
                   + f_0 * lg_592[k];

        t_458[k] = -2.0 * ig_338[k]
                   + f_0 * lg_593[k];

        t_459[k] = -2.0 * ig_339[k]
                   + f_0 * lg_594[k];

        t_460[k] = -2.0 * ig_340[k]
                   + f_0 * lg_595[k];

        t_461[k] = -2.0 * ig_341[k]
                   + f_0 * lg_596[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, ig_342, ig_343, ig_344, ig_345, \
                         ig_346, lg_597, lg_598, lg_599, lg_600, \
                         lg_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = -2.0 * ig_342[k]
                   + f_0 * lg_597[k];

        t_463[k] = -2.0 * ig_343[k]
                   + f_0 * lg_598[k];

        t_464[k] = -2.0 * ig_344[k]
                   + f_0 * lg_599[k];

        t_465[k] = -3.0 * ig_345[k]
                   + f_0 * lg_600[k];

        t_466[k] = -3.0 * ig_346[k]
                   + f_0 * lg_601[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, ig_347, ig_348, ig_349, ig_350, \
                         ig_351, lg_602, lg_603, lg_604, lg_605, \
                         lg_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = -3.0 * ig_347[k]
                   + f_0 * lg_602[k];

        t_468[k] = -3.0 * ig_348[k]
                   + f_0 * lg_603[k];

        t_469[k] = -3.0 * ig_349[k]
                   + f_0 * lg_604[k];

        t_470[k] = -3.0 * ig_350[k]
                   + f_0 * lg_605[k];

        t_471[k] = -3.0 * ig_351[k]
                   + f_0 * lg_606[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, ig_352, ig_353, ig_354, ig_355, \
                         ig_356, lg_607, lg_608, lg_609, lg_610, \
                         lg_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = -3.0 * ig_352[k]
                   + f_0 * lg_607[k];

        t_473[k] = -3.0 * ig_353[k]
                   + f_0 * lg_608[k];

        t_474[k] = -3.0 * ig_354[k]
                   + f_0 * lg_609[k];

        t_475[k] = -3.0 * ig_355[k]
                   + f_0 * lg_610[k];

        t_476[k] = -3.0 * ig_356[k]
                   + f_0 * lg_611[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, ig_357, ig_358, ig_359, ig_360, \
                         ig_361, lg_612, lg_613, lg_614, lg_615, \
                         lg_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = -3.0 * ig_357[k]
                   + f_0 * lg_612[k];

        t_478[k] = -3.0 * ig_358[k]
                   + f_0 * lg_613[k];

        t_479[k] = -3.0 * ig_359[k]
                   + f_0 * lg_614[k];

        t_480[k] = -4.0 * ig_360[k]
                   + f_0 * lg_615[k];

        t_481[k] = -4.0 * ig_361[k]
                   + f_0 * lg_616[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, t_486, ig_362, ig_363, ig_364, ig_365, \
                         ig_366, lg_617, lg_618, lg_619, lg_620, \
                         lg_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = -4.0 * ig_362[k]
                   + f_0 * lg_617[k];

        t_483[k] = -4.0 * ig_363[k]
                   + f_0 * lg_618[k];

        t_484[k] = -4.0 * ig_364[k]
                   + f_0 * lg_619[k];

        t_485[k] = -4.0 * ig_365[k]
                   + f_0 * lg_620[k];

        t_486[k] = -4.0 * ig_366[k]
                   + f_0 * lg_621[k];
    }
}

static auto
compute_prim_geom_10_kg_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ig, const size_t lg,
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

    const auto *lg_622 = buffer.data(lg + 622);
    const auto *lg_623 = buffer.data(lg + 623);
    const auto *lg_624 = buffer.data(lg + 624);
    const auto *lg_625 = buffer.data(lg + 625);
    const auto *lg_626 = buffer.data(lg + 626);
    const auto *lg_627 = buffer.data(lg + 627);
    const auto *lg_628 = buffer.data(lg + 628);
    const auto *lg_629 = buffer.data(lg + 629);
    const auto *lg_630 = buffer.data(lg + 630);
    const auto *lg_631 = buffer.data(lg + 631);
    const auto *lg_632 = buffer.data(lg + 632);
    const auto *lg_633 = buffer.data(lg + 633);
    const auto *lg_634 = buffer.data(lg + 634);
    const auto *lg_635 = buffer.data(lg + 635);
    const auto *lg_636 = buffer.data(lg + 636);
    const auto *lg_637 = buffer.data(lg + 637);
    const auto *lg_638 = buffer.data(lg + 638);
    const auto *lg_639 = buffer.data(lg + 639);
    const auto *lg_640 = buffer.data(lg + 640);
    const auto *lg_641 = buffer.data(lg + 641);
    const auto *lg_642 = buffer.data(lg + 642);
    const auto *lg_643 = buffer.data(lg + 643);
    const auto *lg_644 = buffer.data(lg + 644);
    const auto *lg_645 = buffer.data(lg + 645);
    const auto *lg_646 = buffer.data(lg + 646);
    const auto *lg_647 = buffer.data(lg + 647);
    const auto *lg_648 = buffer.data(lg + 648);
    const auto *lg_649 = buffer.data(lg + 649);
    const auto *lg_650 = buffer.data(lg + 650);
    const auto *lg_651 = buffer.data(lg + 651);
    const auto *lg_652 = buffer.data(lg + 652);
    const auto *lg_653 = buffer.data(lg + 653);
    const auto *lg_654 = buffer.data(lg + 654);
    const auto *lg_655 = buffer.data(lg + 655);
    const auto *lg_656 = buffer.data(lg + 656);
    const auto *lg_657 = buffer.data(lg + 657);
    const auto *lg_658 = buffer.data(lg + 658);
    const auto *lg_659 = buffer.data(lg + 659);
    const auto *lg_660 = buffer.data(lg + 660);
    const auto *lg_661 = buffer.data(lg + 661);
    const auto *lg_662 = buffer.data(lg + 662);
    const auto *lg_663 = buffer.data(lg + 663);
    const auto *lg_664 = buffer.data(lg + 664);
    const auto *lg_665 = buffer.data(lg + 665);
    const auto *lg_666 = buffer.data(lg + 666);
    const auto *lg_667 = buffer.data(lg + 667);
    const auto *lg_668 = buffer.data(lg + 668);
    const auto *lg_669 = buffer.data(lg + 669);
    const auto *lg_670 = buffer.data(lg + 670);
    const auto *lg_671 = buffer.data(lg + 671);
    const auto *lg_672 = buffer.data(lg + 672);
    const auto *lg_673 = buffer.data(lg + 673);
    const auto *lg_674 = buffer.data(lg + 674);

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, ig_367, ig_368, ig_369, ig_370, \
                         ig_371, lg_622, lg_623, lg_624, lg_625, \
                         lg_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = -4.0 * ig_367[k]
                   + f_0 * lg_622[k];

        t_488[k] = -4.0 * ig_368[k]
                   + f_0 * lg_623[k];

        t_489[k] = -4.0 * ig_369[k]
                   + f_0 * lg_624[k];

        t_490[k] = -4.0 * ig_370[k]
                   + f_0 * lg_625[k];

        t_491[k] = -4.0 * ig_371[k]
                   + f_0 * lg_626[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, ig_372, ig_373, ig_374, ig_375, \
                         ig_376, lg_627, lg_628, lg_629, lg_630, \
                         lg_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = -4.0 * ig_372[k]
                   + f_0 * lg_627[k];

        t_493[k] = -4.0 * ig_373[k]
                   + f_0 * lg_628[k];

        t_494[k] = -4.0 * ig_374[k]
                   + f_0 * lg_629[k];

        t_495[k] = -5.0 * ig_375[k]
                   + f_0 * lg_630[k];

        t_496[k] = -5.0 * ig_376[k]
                   + f_0 * lg_631[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, t_500, t_501, ig_377, ig_378, ig_379, ig_380, \
                         ig_381, lg_632, lg_633, lg_634, lg_635, \
                         lg_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = -5.0 * ig_377[k]
                   + f_0 * lg_632[k];

        t_498[k] = -5.0 * ig_378[k]
                   + f_0 * lg_633[k];

        t_499[k] = -5.0 * ig_379[k]
                   + f_0 * lg_634[k];

        t_500[k] = -5.0 * ig_380[k]
                   + f_0 * lg_635[k];

        t_501[k] = -5.0 * ig_381[k]
                   + f_0 * lg_636[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, ig_382, ig_383, ig_384, ig_385, \
                         ig_386, lg_637, lg_638, lg_639, lg_640, \
                         lg_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = -5.0 * ig_382[k]
                   + f_0 * lg_637[k];

        t_503[k] = -5.0 * ig_383[k]
                   + f_0 * lg_638[k];

        t_504[k] = -5.0 * ig_384[k]
                   + f_0 * lg_639[k];

        t_505[k] = -5.0 * ig_385[k]
                   + f_0 * lg_640[k];

        t_506[k] = -5.0 * ig_386[k]
                   + f_0 * lg_641[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, t_510, t_511, ig_387, ig_388, ig_389, ig_390, \
                         ig_391, lg_642, lg_643, lg_644, lg_645, \
                         lg_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = -5.0 * ig_387[k]
                   + f_0 * lg_642[k];

        t_508[k] = -5.0 * ig_388[k]
                   + f_0 * lg_643[k];

        t_509[k] = -5.0 * ig_389[k]
                   + f_0 * lg_644[k];

        t_510[k] = -6.0 * ig_390[k]
                   + f_0 * lg_645[k];

        t_511[k] = -6.0 * ig_391[k]
                   + f_0 * lg_646[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, ig_392, ig_393, ig_394, ig_395, \
                         ig_396, lg_647, lg_648, lg_649, lg_650, \
                         lg_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = -6.0 * ig_392[k]
                   + f_0 * lg_647[k];

        t_513[k] = -6.0 * ig_393[k]
                   + f_0 * lg_648[k];

        t_514[k] = -6.0 * ig_394[k]
                   + f_0 * lg_649[k];

        t_515[k] = -6.0 * ig_395[k]
                   + f_0 * lg_650[k];

        t_516[k] = -6.0 * ig_396[k]
                   + f_0 * lg_651[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, ig_397, ig_398, ig_399, ig_400, \
                         ig_401, lg_652, lg_653, lg_654, lg_655, \
                         lg_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = -6.0 * ig_397[k]
                   + f_0 * lg_652[k];

        t_518[k] = -6.0 * ig_398[k]
                   + f_0 * lg_653[k];

        t_519[k] = -6.0 * ig_399[k]
                   + f_0 * lg_654[k];

        t_520[k] = -6.0 * ig_400[k]
                   + f_0 * lg_655[k];

        t_521[k] = -6.0 * ig_401[k]
                   + f_0 * lg_656[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, ig_402, ig_403, ig_404, ig_405, \
                         ig_406, lg_657, lg_658, lg_659, lg_660, \
                         lg_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = -6.0 * ig_402[k]
                   + f_0 * lg_657[k];

        t_523[k] = -6.0 * ig_403[k]
                   + f_0 * lg_658[k];

        t_524[k] = -6.0 * ig_404[k]
                   + f_0 * lg_659[k];

        t_525[k] = -7.0 * ig_405[k]
                   + f_0 * lg_660[k];

        t_526[k] = -7.0 * ig_406[k]
                   + f_0 * lg_661[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, ig_407, ig_408, ig_409, ig_410, \
                         ig_411, lg_662, lg_663, lg_664, lg_665, \
                         lg_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = -7.0 * ig_407[k]
                   + f_0 * lg_662[k];

        t_528[k] = -7.0 * ig_408[k]
                   + f_0 * lg_663[k];

        t_529[k] = -7.0 * ig_409[k]
                   + f_0 * lg_664[k];

        t_530[k] = -7.0 * ig_410[k]
                   + f_0 * lg_665[k];

        t_531[k] = -7.0 * ig_411[k]
                   + f_0 * lg_666[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, ig_412, ig_413, ig_414, ig_415, \
                         ig_416, lg_667, lg_668, lg_669, lg_670, \
                         lg_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = -7.0 * ig_412[k]
                   + f_0 * lg_667[k];

        t_533[k] = -7.0 * ig_413[k]
                   + f_0 * lg_668[k];

        t_534[k] = -7.0 * ig_414[k]
                   + f_0 * lg_669[k];

        t_535[k] = -7.0 * ig_415[k]
                   + f_0 * lg_670[k];

        t_536[k] = -7.0 * ig_416[k]
                   + f_0 * lg_671[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, ig_417, ig_418, ig_419, lg_672, lg_673, \
                         lg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = -7.0 * ig_417[k]
                   + f_0 * lg_672[k];

        t_538[k] = -7.0 * ig_418[k]
                   + f_0 * lg_673[k];

        t_539[k] = -7.0 * ig_419[k]
                   + f_0 * lg_674[k];
    }
}

auto
compute_prim_geom_10_kg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t ig, const size_t lg,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kg_electron_repulsion_2_piece0(buffer, target, ig, lg, ncols, alpha);

    compute_prim_geom_10_kg_electron_repulsion_2_piece1(buffer, target, ig, lg, ncols, alpha);

    compute_prim_geom_10_kg_electron_repulsion_2_piece2(buffer, target, ig, lg, ncols, alpha);

    compute_prim_geom_10_kg_electron_repulsion_2_piece3(buffer, target, ig, lg, ncols, alpha);
}

}  // namespace simdt2ceri
