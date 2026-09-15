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


#include "SimdElectronRepulsionGeom10VrrRecLI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_li_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
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

    const auto *ki_0 = buffer.data(ki + 0);
    const auto *ki_1 = buffer.data(ki + 1);
    const auto *ki_2 = buffer.data(ki + 2);
    const auto *ki_3 = buffer.data(ki + 3);
    const auto *ki_4 = buffer.data(ki + 4);
    const auto *ki_5 = buffer.data(ki + 5);
    const auto *ki_6 = buffer.data(ki + 6);
    const auto *ki_7 = buffer.data(ki + 7);
    const auto *ki_8 = buffer.data(ki + 8);
    const auto *ki_9 = buffer.data(ki + 9);
    const auto *ki_10 = buffer.data(ki + 10);
    const auto *ki_11 = buffer.data(ki + 11);
    const auto *ki_12 = buffer.data(ki + 12);
    const auto *ki_13 = buffer.data(ki + 13);
    const auto *ki_14 = buffer.data(ki + 14);
    const auto *ki_15 = buffer.data(ki + 15);
    const auto *ki_16 = buffer.data(ki + 16);
    const auto *ki_17 = buffer.data(ki + 17);
    const auto *ki_18 = buffer.data(ki + 18);
    const auto *ki_19 = buffer.data(ki + 19);
    const auto *ki_20 = buffer.data(ki + 20);
    const auto *ki_21 = buffer.data(ki + 21);
    const auto *ki_22 = buffer.data(ki + 22);
    const auto *ki_23 = buffer.data(ki + 23);
    const auto *ki_24 = buffer.data(ki + 24);
    const auto *ki_25 = buffer.data(ki + 25);
    const auto *ki_26 = buffer.data(ki + 26);
    const auto *ki_27 = buffer.data(ki + 27);
    const auto *ki_28 = buffer.data(ki + 28);
    const auto *ki_29 = buffer.data(ki + 29);
    const auto *ki_30 = buffer.data(ki + 30);
    const auto *ki_31 = buffer.data(ki + 31);
    const auto *ki_32 = buffer.data(ki + 32);
    const auto *ki_33 = buffer.data(ki + 33);
    const auto *ki_34 = buffer.data(ki + 34);
    const auto *ki_35 = buffer.data(ki + 35);
    const auto *ki_36 = buffer.data(ki + 36);
    const auto *ki_37 = buffer.data(ki + 37);
    const auto *ki_38 = buffer.data(ki + 38);
    const auto *ki_39 = buffer.data(ki + 39);
    const auto *ki_40 = buffer.data(ki + 40);
    const auto *ki_41 = buffer.data(ki + 41);
    const auto *ki_42 = buffer.data(ki + 42);
    const auto *ki_43 = buffer.data(ki + 43);
    const auto *ki_44 = buffer.data(ki + 44);
    const auto *ki_45 = buffer.data(ki + 45);
    const auto *ki_46 = buffer.data(ki + 46);
    const auto *ki_47 = buffer.data(ki + 47);
    const auto *ki_48 = buffer.data(ki + 48);
    const auto *ki_49 = buffer.data(ki + 49);
    const auto *ki_50 = buffer.data(ki + 50);
    const auto *ki_51 = buffer.data(ki + 51);
    const auto *ki_52 = buffer.data(ki + 52);
    const auto *ki_53 = buffer.data(ki + 53);
    const auto *ki_54 = buffer.data(ki + 54);
    const auto *ki_55 = buffer.data(ki + 55);
    const auto *ki_56 = buffer.data(ki + 56);
    const auto *ki_57 = buffer.data(ki + 57);
    const auto *ki_58 = buffer.data(ki + 58);
    const auto *ki_59 = buffer.data(ki + 59);
    const auto *ki_60 = buffer.data(ki + 60);
    const auto *ki_61 = buffer.data(ki + 61);
    const auto *ki_62 = buffer.data(ki + 62);
    const auto *ki_63 = buffer.data(ki + 63);
    const auto *ki_64 = buffer.data(ki + 64);
    const auto *ki_65 = buffer.data(ki + 65);
    const auto *ki_66 = buffer.data(ki + 66);
    const auto *ki_67 = buffer.data(ki + 67);
    const auto *ki_68 = buffer.data(ki + 68);
    const auto *ki_69 = buffer.data(ki + 69);
    const auto *ki_70 = buffer.data(ki + 70);
    const auto *ki_71 = buffer.data(ki + 71);
    const auto *ki_72 = buffer.data(ki + 72);
    const auto *ki_73 = buffer.data(ki + 73);
    const auto *ki_74 = buffer.data(ki + 74);
    const auto *ki_75 = buffer.data(ki + 75);
    const auto *ki_76 = buffer.data(ki + 76);
    const auto *ki_77 = buffer.data(ki + 77);
    const auto *ki_78 = buffer.data(ki + 78);
    const auto *ki_79 = buffer.data(ki + 79);
    const auto *ki_80 = buffer.data(ki + 80);
    const auto *ki_81 = buffer.data(ki + 81);
    const auto *ki_82 = buffer.data(ki + 82);
    const auto *ki_83 = buffer.data(ki + 83);
    const auto *ki_84 = buffer.data(ki + 84);
    const auto *ki_85 = buffer.data(ki + 85);
    const auto *ki_86 = buffer.data(ki + 86);
    const auto *ki_87 = buffer.data(ki + 87);
    const auto *ki_88 = buffer.data(ki + 88);
    const auto *ki_89 = buffer.data(ki + 89);
    const auto *ki_90 = buffer.data(ki + 90);
    const auto *ki_91 = buffer.data(ki + 91);
    const auto *ki_92 = buffer.data(ki + 92);
    const auto *ki_93 = buffer.data(ki + 93);
    const auto *ki_94 = buffer.data(ki + 94);
    const auto *ki_95 = buffer.data(ki + 95);
    const auto *ki_96 = buffer.data(ki + 96);
    const auto *ki_97 = buffer.data(ki + 97);
    const auto *ki_98 = buffer.data(ki + 98);
    const auto *ki_99 = buffer.data(ki + 99);
    const auto *ki_100 = buffer.data(ki + 100);
    const auto *ki_101 = buffer.data(ki + 101);
    const auto *ki_102 = buffer.data(ki + 102);
    const auto *ki_103 = buffer.data(ki + 103);
    const auto *ki_104 = buffer.data(ki + 104);
    const auto *ki_105 = buffer.data(ki + 105);
    const auto *ki_106 = buffer.data(ki + 106);
    const auto *ki_107 = buffer.data(ki + 107);
    const auto *ki_108 = buffer.data(ki + 108);
    const auto *ki_109 = buffer.data(ki + 109);
    const auto *ki_110 = buffer.data(ki + 110);
    const auto *ki_111 = buffer.data(ki + 111);
    const auto *ki_112 = buffer.data(ki + 112);
    const auto *ki_113 = buffer.data(ki + 113);
    const auto *ki_114 = buffer.data(ki + 114);
    const auto *ki_115 = buffer.data(ki + 115);
    const auto *ki_116 = buffer.data(ki + 116);
    const auto *ki_117 = buffer.data(ki + 117);
    const auto *ki_118 = buffer.data(ki + 118);
    const auto *ki_119 = buffer.data(ki + 119);
    const auto *ki_120 = buffer.data(ki + 120);
    const auto *ki_121 = buffer.data(ki + 121);
    const auto *ki_122 = buffer.data(ki + 122);
    const auto *ki_123 = buffer.data(ki + 123);
    const auto *ki_124 = buffer.data(ki + 124);
    const auto *ki_125 = buffer.data(ki + 125);
    const auto *ki_126 = buffer.data(ki + 126);
    const auto *ki_127 = buffer.data(ki + 127);
    const auto *ki_128 = buffer.data(ki + 128);
    const auto *ki_129 = buffer.data(ki + 129);
    const auto *ki_130 = buffer.data(ki + 130);
    const auto *ki_131 = buffer.data(ki + 131);
    const auto *ki_132 = buffer.data(ki + 132);
    const auto *ki_133 = buffer.data(ki + 133);
    const auto *ki_134 = buffer.data(ki + 134);
    const auto *ki_135 = buffer.data(ki + 135);
    const auto *ki_136 = buffer.data(ki + 136);
    const auto *ki_137 = buffer.data(ki + 137);
    const auto *ki_138 = buffer.data(ki + 138);
    const auto *ki_139 = buffer.data(ki + 139);
    const auto *ki_140 = buffer.data(ki + 140);
    const auto *ki_141 = buffer.data(ki + 141);
    const auto *ki_142 = buffer.data(ki + 142);
    const auto *ki_143 = buffer.data(ki + 143);
    const auto *ki_144 = buffer.data(ki + 144);
    const auto *ki_145 = buffer.data(ki + 145);
    const auto *ki_146 = buffer.data(ki + 146);
    const auto *ki_147 = buffer.data(ki + 147);
    const auto *ki_148 = buffer.data(ki + 148);
    const auto *ki_149 = buffer.data(ki + 149);

    const auto *mi_0 = buffer.data(mi + 0);
    const auto *mi_1 = buffer.data(mi + 1);
    const auto *mi_2 = buffer.data(mi + 2);
    const auto *mi_3 = buffer.data(mi + 3);
    const auto *mi_4 = buffer.data(mi + 4);
    const auto *mi_5 = buffer.data(mi + 5);
    const auto *mi_6 = buffer.data(mi + 6);
    const auto *mi_7 = buffer.data(mi + 7);
    const auto *mi_8 = buffer.data(mi + 8);
    const auto *mi_9 = buffer.data(mi + 9);
    const auto *mi_10 = buffer.data(mi + 10);
    const auto *mi_11 = buffer.data(mi + 11);
    const auto *mi_12 = buffer.data(mi + 12);
    const auto *mi_13 = buffer.data(mi + 13);
    const auto *mi_14 = buffer.data(mi + 14);
    const auto *mi_15 = buffer.data(mi + 15);
    const auto *mi_16 = buffer.data(mi + 16);
    const auto *mi_17 = buffer.data(mi + 17);
    const auto *mi_18 = buffer.data(mi + 18);
    const auto *mi_19 = buffer.data(mi + 19);
    const auto *mi_20 = buffer.data(mi + 20);
    const auto *mi_21 = buffer.data(mi + 21);
    const auto *mi_22 = buffer.data(mi + 22);
    const auto *mi_23 = buffer.data(mi + 23);
    const auto *mi_24 = buffer.data(mi + 24);
    const auto *mi_25 = buffer.data(mi + 25);
    const auto *mi_26 = buffer.data(mi + 26);
    const auto *mi_27 = buffer.data(mi + 27);
    const auto *mi_28 = buffer.data(mi + 28);
    const auto *mi_29 = buffer.data(mi + 29);
    const auto *mi_30 = buffer.data(mi + 30);
    const auto *mi_31 = buffer.data(mi + 31);
    const auto *mi_32 = buffer.data(mi + 32);
    const auto *mi_33 = buffer.data(mi + 33);
    const auto *mi_34 = buffer.data(mi + 34);
    const auto *mi_35 = buffer.data(mi + 35);
    const auto *mi_36 = buffer.data(mi + 36);
    const auto *mi_37 = buffer.data(mi + 37);
    const auto *mi_38 = buffer.data(mi + 38);
    const auto *mi_39 = buffer.data(mi + 39);
    const auto *mi_40 = buffer.data(mi + 40);
    const auto *mi_41 = buffer.data(mi + 41);
    const auto *mi_42 = buffer.data(mi + 42);
    const auto *mi_43 = buffer.data(mi + 43);
    const auto *mi_44 = buffer.data(mi + 44);
    const auto *mi_45 = buffer.data(mi + 45);
    const auto *mi_46 = buffer.data(mi + 46);
    const auto *mi_47 = buffer.data(mi + 47);
    const auto *mi_48 = buffer.data(mi + 48);
    const auto *mi_49 = buffer.data(mi + 49);
    const auto *mi_50 = buffer.data(mi + 50);
    const auto *mi_51 = buffer.data(mi + 51);
    const auto *mi_52 = buffer.data(mi + 52);
    const auto *mi_53 = buffer.data(mi + 53);
    const auto *mi_54 = buffer.data(mi + 54);
    const auto *mi_55 = buffer.data(mi + 55);
    const auto *mi_56 = buffer.data(mi + 56);
    const auto *mi_57 = buffer.data(mi + 57);
    const auto *mi_58 = buffer.data(mi + 58);
    const auto *mi_59 = buffer.data(mi + 59);
    const auto *mi_60 = buffer.data(mi + 60);
    const auto *mi_61 = buffer.data(mi + 61);
    const auto *mi_62 = buffer.data(mi + 62);
    const auto *mi_63 = buffer.data(mi + 63);
    const auto *mi_64 = buffer.data(mi + 64);
    const auto *mi_65 = buffer.data(mi + 65);
    const auto *mi_66 = buffer.data(mi + 66);
    const auto *mi_67 = buffer.data(mi + 67);
    const auto *mi_68 = buffer.data(mi + 68);
    const auto *mi_69 = buffer.data(mi + 69);
    const auto *mi_70 = buffer.data(mi + 70);
    const auto *mi_71 = buffer.data(mi + 71);
    const auto *mi_72 = buffer.data(mi + 72);
    const auto *mi_73 = buffer.data(mi + 73);
    const auto *mi_74 = buffer.data(mi + 74);
    const auto *mi_75 = buffer.data(mi + 75);
    const auto *mi_76 = buffer.data(mi + 76);
    const auto *mi_77 = buffer.data(mi + 77);
    const auto *mi_78 = buffer.data(mi + 78);
    const auto *mi_79 = buffer.data(mi + 79);
    const auto *mi_80 = buffer.data(mi + 80);
    const auto *mi_81 = buffer.data(mi + 81);
    const auto *mi_82 = buffer.data(mi + 82);
    const auto *mi_83 = buffer.data(mi + 83);
    const auto *mi_84 = buffer.data(mi + 84);
    const auto *mi_85 = buffer.data(mi + 85);
    const auto *mi_86 = buffer.data(mi + 86);
    const auto *mi_87 = buffer.data(mi + 87);
    const auto *mi_88 = buffer.data(mi + 88);
    const auto *mi_89 = buffer.data(mi + 89);
    const auto *mi_90 = buffer.data(mi + 90);
    const auto *mi_91 = buffer.data(mi + 91);
    const auto *mi_92 = buffer.data(mi + 92);
    const auto *mi_93 = buffer.data(mi + 93);
    const auto *mi_94 = buffer.data(mi + 94);
    const auto *mi_95 = buffer.data(mi + 95);
    const auto *mi_96 = buffer.data(mi + 96);
    const auto *mi_97 = buffer.data(mi + 97);
    const auto *mi_98 = buffer.data(mi + 98);
    const auto *mi_99 = buffer.data(mi + 99);
    const auto *mi_100 = buffer.data(mi + 100);
    const auto *mi_101 = buffer.data(mi + 101);
    const auto *mi_102 = buffer.data(mi + 102);
    const auto *mi_103 = buffer.data(mi + 103);
    const auto *mi_104 = buffer.data(mi + 104);
    const auto *mi_105 = buffer.data(mi + 105);
    const auto *mi_106 = buffer.data(mi + 106);
    const auto *mi_107 = buffer.data(mi + 107);
    const auto *mi_108 = buffer.data(mi + 108);
    const auto *mi_109 = buffer.data(mi + 109);
    const auto *mi_110 = buffer.data(mi + 110);
    const auto *mi_111 = buffer.data(mi + 111);
    const auto *mi_112 = buffer.data(mi + 112);
    const auto *mi_113 = buffer.data(mi + 113);
    const auto *mi_114 = buffer.data(mi + 114);
    const auto *mi_115 = buffer.data(mi + 115);
    const auto *mi_116 = buffer.data(mi + 116);
    const auto *mi_117 = buffer.data(mi + 117);
    const auto *mi_118 = buffer.data(mi + 118);
    const auto *mi_119 = buffer.data(mi + 119);
    const auto *mi_120 = buffer.data(mi + 120);
    const auto *mi_121 = buffer.data(mi + 121);
    const auto *mi_122 = buffer.data(mi + 122);
    const auto *mi_123 = buffer.data(mi + 123);
    const auto *mi_124 = buffer.data(mi + 124);
    const auto *mi_125 = buffer.data(mi + 125);
    const auto *mi_126 = buffer.data(mi + 126);
    const auto *mi_127 = buffer.data(mi + 127);
    const auto *mi_128 = buffer.data(mi + 128);
    const auto *mi_129 = buffer.data(mi + 129);
    const auto *mi_130 = buffer.data(mi + 130);
    const auto *mi_131 = buffer.data(mi + 131);
    const auto *mi_132 = buffer.data(mi + 132);
    const auto *mi_133 = buffer.data(mi + 133);
    const auto *mi_134 = buffer.data(mi + 134);
    const auto *mi_135 = buffer.data(mi + 135);
    const auto *mi_136 = buffer.data(mi + 136);
    const auto *mi_137 = buffer.data(mi + 137);
    const auto *mi_138 = buffer.data(mi + 138);
    const auto *mi_139 = buffer.data(mi + 139);
    const auto *mi_140 = buffer.data(mi + 140);
    const auto *mi_141 = buffer.data(mi + 141);
    const auto *mi_142 = buffer.data(mi + 142);
    const auto *mi_143 = buffer.data(mi + 143);
    const auto *mi_144 = buffer.data(mi + 144);
    const auto *mi_145 = buffer.data(mi + 145);
    const auto *mi_146 = buffer.data(mi + 146);
    const auto *mi_147 = buffer.data(mi + 147);
    const auto *mi_148 = buffer.data(mi + 148);
    const auto *mi_149 = buffer.data(mi + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ki_0, ki_1, ki_2, ki_3, ki_4, mi_0, mi_1, \
                         mi_2, mi_3, mi_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -8.0 * ki_0[k]
                 + f_0 * mi_0[k];

        t_1[k] = -8.0 * ki_1[k]
                 + f_0 * mi_1[k];

        t_2[k] = -8.0 * ki_2[k]
                 + f_0 * mi_2[k];

        t_3[k] = -8.0 * ki_3[k]
                 + f_0 * mi_3[k];

        t_4[k] = -8.0 * ki_4[k]
                 + f_0 * mi_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ki_5, ki_6, ki_7, ki_8, ki_9, mi_5, mi_6, \
                         mi_7, mi_8, mi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -8.0 * ki_5[k]
                 + f_0 * mi_5[k];

        t_6[k] = -8.0 * ki_6[k]
                 + f_0 * mi_6[k];

        t_7[k] = -8.0 * ki_7[k]
                 + f_0 * mi_7[k];

        t_8[k] = -8.0 * ki_8[k]
                 + f_0 * mi_8[k];

        t_9[k] = -8.0 * ki_9[k]
                 + f_0 * mi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ki_10, ki_11, ki_12, ki_13, ki_14, \
                         mi_10, mi_11, mi_12, mi_13, mi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -8.0 * ki_10[k]
                  + f_0 * mi_10[k];

        t_11[k] = -8.0 * ki_11[k]
                  + f_0 * mi_11[k];

        t_12[k] = -8.0 * ki_12[k]
                  + f_0 * mi_12[k];

        t_13[k] = -8.0 * ki_13[k]
                  + f_0 * mi_13[k];

        t_14[k] = -8.0 * ki_14[k]
                  + f_0 * mi_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ki_15, ki_16, ki_17, ki_18, ki_19, \
                         mi_15, mi_16, mi_17, mi_18, mi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -8.0 * ki_15[k]
                  + f_0 * mi_15[k];

        t_16[k] = -8.0 * ki_16[k]
                  + f_0 * mi_16[k];

        t_17[k] = -8.0 * ki_17[k]
                  + f_0 * mi_17[k];

        t_18[k] = -8.0 * ki_18[k]
                  + f_0 * mi_18[k];

        t_19[k] = -8.0 * ki_19[k]
                  + f_0 * mi_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ki_20, ki_21, ki_22, ki_23, ki_24, \
                         mi_20, mi_21, mi_22, mi_23, mi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -8.0 * ki_20[k]
                  + f_0 * mi_20[k];

        t_21[k] = -8.0 * ki_21[k]
                  + f_0 * mi_21[k];

        t_22[k] = -8.0 * ki_22[k]
                  + f_0 * mi_22[k];

        t_23[k] = -8.0 * ki_23[k]
                  + f_0 * mi_23[k];

        t_24[k] = -8.0 * ki_24[k]
                  + f_0 * mi_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ki_25, ki_26, ki_27, ki_28, ki_29, \
                         mi_25, mi_26, mi_27, mi_28, mi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -8.0 * ki_25[k]
                  + f_0 * mi_25[k];

        t_26[k] = -8.0 * ki_26[k]
                  + f_0 * mi_26[k];

        t_27[k] = -8.0 * ki_27[k]
                  + f_0 * mi_27[k];

        t_28[k] = -7.0 * ki_28[k]
                  + f_0 * mi_28[k];

        t_29[k] = -7.0 * ki_29[k]
                  + f_0 * mi_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ki_30, ki_31, ki_32, ki_33, ki_34, \
                         mi_30, mi_31, mi_32, mi_33, mi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -7.0 * ki_30[k]
                  + f_0 * mi_30[k];

        t_31[k] = -7.0 * ki_31[k]
                  + f_0 * mi_31[k];

        t_32[k] = -7.0 * ki_32[k]
                  + f_0 * mi_32[k];

        t_33[k] = -7.0 * ki_33[k]
                  + f_0 * mi_33[k];

        t_34[k] = -7.0 * ki_34[k]
                  + f_0 * mi_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ki_35, ki_36, ki_37, ki_38, ki_39, \
                         mi_35, mi_36, mi_37, mi_38, mi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -7.0 * ki_35[k]
                  + f_0 * mi_35[k];

        t_36[k] = -7.0 * ki_36[k]
                  + f_0 * mi_36[k];

        t_37[k] = -7.0 * ki_37[k]
                  + f_0 * mi_37[k];

        t_38[k] = -7.0 * ki_38[k]
                  + f_0 * mi_38[k];

        t_39[k] = -7.0 * ki_39[k]
                  + f_0 * mi_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ki_40, ki_41, ki_42, ki_43, ki_44, \
                         mi_40, mi_41, mi_42, mi_43, mi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -7.0 * ki_40[k]
                  + f_0 * mi_40[k];

        t_41[k] = -7.0 * ki_41[k]
                  + f_0 * mi_41[k];

        t_42[k] = -7.0 * ki_42[k]
                  + f_0 * mi_42[k];

        t_43[k] = -7.0 * ki_43[k]
                  + f_0 * mi_43[k];

        t_44[k] = -7.0 * ki_44[k]
                  + f_0 * mi_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ki_45, ki_46, ki_47, ki_48, ki_49, \
                         mi_45, mi_46, mi_47, mi_48, mi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -7.0 * ki_45[k]
                  + f_0 * mi_45[k];

        t_46[k] = -7.0 * ki_46[k]
                  + f_0 * mi_46[k];

        t_47[k] = -7.0 * ki_47[k]
                  + f_0 * mi_47[k];

        t_48[k] = -7.0 * ki_48[k]
                  + f_0 * mi_48[k];

        t_49[k] = -7.0 * ki_49[k]
                  + f_0 * mi_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ki_50, ki_51, ki_52, ki_53, ki_54, \
                         mi_50, mi_51, mi_52, mi_53, mi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -7.0 * ki_50[k]
                  + f_0 * mi_50[k];

        t_51[k] = -7.0 * ki_51[k]
                  + f_0 * mi_51[k];

        t_52[k] = -7.0 * ki_52[k]
                  + f_0 * mi_52[k];

        t_53[k] = -7.0 * ki_53[k]
                  + f_0 * mi_53[k];

        t_54[k] = -7.0 * ki_54[k]
                  + f_0 * mi_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ki_55, ki_56, ki_57, ki_58, ki_59, \
                         mi_55, mi_56, mi_57, mi_58, mi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -7.0 * ki_55[k]
                  + f_0 * mi_55[k];

        t_56[k] = -7.0 * ki_56[k]
                  + f_0 * mi_56[k];

        t_57[k] = -7.0 * ki_57[k]
                  + f_0 * mi_57[k];

        t_58[k] = -7.0 * ki_58[k]
                  + f_0 * mi_58[k];

        t_59[k] = -7.0 * ki_59[k]
                  + f_0 * mi_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ki_60, ki_61, ki_62, ki_63, ki_64, \
                         mi_60, mi_61, mi_62, mi_63, mi_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -7.0 * ki_60[k]
                  + f_0 * mi_60[k];

        t_61[k] = -7.0 * ki_61[k]
                  + f_0 * mi_61[k];

        t_62[k] = -7.0 * ki_62[k]
                  + f_0 * mi_62[k];

        t_63[k] = -7.0 * ki_63[k]
                  + f_0 * mi_63[k];

        t_64[k] = -7.0 * ki_64[k]
                  + f_0 * mi_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ki_65, ki_66, ki_67, ki_68, ki_69, \
                         mi_65, mi_66, mi_67, mi_68, mi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -7.0 * ki_65[k]
                  + f_0 * mi_65[k];

        t_66[k] = -7.0 * ki_66[k]
                  + f_0 * mi_66[k];

        t_67[k] = -7.0 * ki_67[k]
                  + f_0 * mi_67[k];

        t_68[k] = -7.0 * ki_68[k]
                  + f_0 * mi_68[k];

        t_69[k] = -7.0 * ki_69[k]
                  + f_0 * mi_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ki_70, ki_71, ki_72, ki_73, ki_74, \
                         mi_70, mi_71, mi_72, mi_73, mi_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -7.0 * ki_70[k]
                  + f_0 * mi_70[k];

        t_71[k] = -7.0 * ki_71[k]
                  + f_0 * mi_71[k];

        t_72[k] = -7.0 * ki_72[k]
                  + f_0 * mi_72[k];

        t_73[k] = -7.0 * ki_73[k]
                  + f_0 * mi_73[k];

        t_74[k] = -7.0 * ki_74[k]
                  + f_0 * mi_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ki_75, ki_76, ki_77, ki_78, ki_79, \
                         mi_75, mi_76, mi_77, mi_78, mi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -7.0 * ki_75[k]
                  + f_0 * mi_75[k];

        t_76[k] = -7.0 * ki_76[k]
                  + f_0 * mi_76[k];

        t_77[k] = -7.0 * ki_77[k]
                  + f_0 * mi_77[k];

        t_78[k] = -7.0 * ki_78[k]
                  + f_0 * mi_78[k];

        t_79[k] = -7.0 * ki_79[k]
                  + f_0 * mi_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ki_80, ki_81, ki_82, ki_83, ki_84, \
                         mi_80, mi_81, mi_82, mi_83, mi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -7.0 * ki_80[k]
                  + f_0 * mi_80[k];

        t_81[k] = -7.0 * ki_81[k]
                  + f_0 * mi_81[k];

        t_82[k] = -7.0 * ki_82[k]
                  + f_0 * mi_82[k];

        t_83[k] = -7.0 * ki_83[k]
                  + f_0 * mi_83[k];

        t_84[k] = -6.0 * ki_84[k]
                  + f_0 * mi_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ki_85, ki_86, ki_87, ki_88, ki_89, \
                         mi_85, mi_86, mi_87, mi_88, mi_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -6.0 * ki_85[k]
                  + f_0 * mi_85[k];

        t_86[k] = -6.0 * ki_86[k]
                  + f_0 * mi_86[k];

        t_87[k] = -6.0 * ki_87[k]
                  + f_0 * mi_87[k];

        t_88[k] = -6.0 * ki_88[k]
                  + f_0 * mi_88[k];

        t_89[k] = -6.0 * ki_89[k]
                  + f_0 * mi_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ki_90, ki_91, ki_92, ki_93, ki_94, \
                         mi_90, mi_91, mi_92, mi_93, mi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -6.0 * ki_90[k]
                  + f_0 * mi_90[k];

        t_91[k] = -6.0 * ki_91[k]
                  + f_0 * mi_91[k];

        t_92[k] = -6.0 * ki_92[k]
                  + f_0 * mi_92[k];

        t_93[k] = -6.0 * ki_93[k]
                  + f_0 * mi_93[k];

        t_94[k] = -6.0 * ki_94[k]
                  + f_0 * mi_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ki_95, ki_96, ki_97, ki_98, ki_99, \
                         mi_95, mi_96, mi_97, mi_98, mi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -6.0 * ki_95[k]
                  + f_0 * mi_95[k];

        t_96[k] = -6.0 * ki_96[k]
                  + f_0 * mi_96[k];

        t_97[k] = -6.0 * ki_97[k]
                  + f_0 * mi_97[k];

        t_98[k] = -6.0 * ki_98[k]
                  + f_0 * mi_98[k];

        t_99[k] = -6.0 * ki_99[k]
                  + f_0 * mi_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ki_100, ki_101, ki_102, ki_103, \
                         ki_104, mi_100, mi_101, mi_102, mi_103, \
                         mi_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -6.0 * ki_100[k]
                   + f_0 * mi_100[k];

        t_101[k] = -6.0 * ki_101[k]
                   + f_0 * mi_101[k];

        t_102[k] = -6.0 * ki_102[k]
                   + f_0 * mi_102[k];

        t_103[k] = -6.0 * ki_103[k]
                   + f_0 * mi_103[k];

        t_104[k] = -6.0 * ki_104[k]
                   + f_0 * mi_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ki_105, ki_106, ki_107, ki_108, \
                         ki_109, mi_105, mi_106, mi_107, mi_108, \
                         mi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -6.0 * ki_105[k]
                   + f_0 * mi_105[k];

        t_106[k] = -6.0 * ki_106[k]
                   + f_0 * mi_106[k];

        t_107[k] = -6.0 * ki_107[k]
                   + f_0 * mi_107[k];

        t_108[k] = -6.0 * ki_108[k]
                   + f_0 * mi_108[k];

        t_109[k] = -6.0 * ki_109[k]
                   + f_0 * mi_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ki_110, ki_111, ki_112, ki_113, \
                         ki_114, mi_110, mi_111, mi_112, mi_113, \
                         mi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -6.0 * ki_110[k]
                   + f_0 * mi_110[k];

        t_111[k] = -6.0 * ki_111[k]
                   + f_0 * mi_111[k];

        t_112[k] = -6.0 * ki_112[k]
                   + f_0 * mi_112[k];

        t_113[k] = -6.0 * ki_113[k]
                   + f_0 * mi_113[k];

        t_114[k] = -6.0 * ki_114[k]
                   + f_0 * mi_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ki_115, ki_116, ki_117, ki_118, \
                         ki_119, mi_115, mi_116, mi_117, mi_118, \
                         mi_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -6.0 * ki_115[k]
                   + f_0 * mi_115[k];

        t_116[k] = -6.0 * ki_116[k]
                   + f_0 * mi_116[k];

        t_117[k] = -6.0 * ki_117[k]
                   + f_0 * mi_117[k];

        t_118[k] = -6.0 * ki_118[k]
                   + f_0 * mi_118[k];

        t_119[k] = -6.0 * ki_119[k]
                   + f_0 * mi_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ki_120, ki_121, ki_122, ki_123, \
                         ki_124, mi_120, mi_121, mi_122, mi_123, \
                         mi_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -6.0 * ki_120[k]
                   + f_0 * mi_120[k];

        t_121[k] = -6.0 * ki_121[k]
                   + f_0 * mi_121[k];

        t_122[k] = -6.0 * ki_122[k]
                   + f_0 * mi_122[k];

        t_123[k] = -6.0 * ki_123[k]
                   + f_0 * mi_123[k];

        t_124[k] = -6.0 * ki_124[k]
                   + f_0 * mi_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ki_125, ki_126, ki_127, ki_128, \
                         ki_129, mi_125, mi_126, mi_127, mi_128, \
                         mi_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -6.0 * ki_125[k]
                   + f_0 * mi_125[k];

        t_126[k] = -6.0 * ki_126[k]
                   + f_0 * mi_126[k];

        t_127[k] = -6.0 * ki_127[k]
                   + f_0 * mi_127[k];

        t_128[k] = -6.0 * ki_128[k]
                   + f_0 * mi_128[k];

        t_129[k] = -6.0 * ki_129[k]
                   + f_0 * mi_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ki_130, ki_131, ki_132, ki_133, \
                         ki_134, mi_130, mi_131, mi_132, mi_133, \
                         mi_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -6.0 * ki_130[k]
                   + f_0 * mi_130[k];

        t_131[k] = -6.0 * ki_131[k]
                   + f_0 * mi_131[k];

        t_132[k] = -6.0 * ki_132[k]
                   + f_0 * mi_132[k];

        t_133[k] = -6.0 * ki_133[k]
                   + f_0 * mi_133[k];

        t_134[k] = -6.0 * ki_134[k]
                   + f_0 * mi_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ki_135, ki_136, ki_137, ki_138, \
                         ki_139, mi_135, mi_136, mi_137, mi_138, \
                         mi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -6.0 * ki_135[k]
                   + f_0 * mi_135[k];

        t_136[k] = -6.0 * ki_136[k]
                   + f_0 * mi_136[k];

        t_137[k] = -6.0 * ki_137[k]
                   + f_0 * mi_137[k];

        t_138[k] = -6.0 * ki_138[k]
                   + f_0 * mi_138[k];

        t_139[k] = -6.0 * ki_139[k]
                   + f_0 * mi_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ki_140, ki_141, ki_142, ki_143, \
                         ki_144, mi_140, mi_141, mi_142, mi_143, \
                         mi_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -6.0 * ki_140[k]
                   + f_0 * mi_140[k];

        t_141[k] = -6.0 * ki_141[k]
                   + f_0 * mi_141[k];

        t_142[k] = -6.0 * ki_142[k]
                   + f_0 * mi_142[k];

        t_143[k] = -6.0 * ki_143[k]
                   + f_0 * mi_143[k];

        t_144[k] = -6.0 * ki_144[k]
                   + f_0 * mi_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ki_145, ki_146, ki_147, ki_148, \
                         ki_149, mi_145, mi_146, mi_147, mi_148, \
                         mi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -6.0 * ki_145[k]
                   + f_0 * mi_145[k];

        t_146[k] = -6.0 * ki_146[k]
                   + f_0 * mi_146[k];

        t_147[k] = -6.0 * ki_147[k]
                   + f_0 * mi_147[k];

        t_148[k] = -6.0 * ki_148[k]
                   + f_0 * mi_148[k];

        t_149[k] = -6.0 * ki_149[k]
                   + f_0 * mi_149[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
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

    const auto *ki_150 = buffer.data(ki + 150);
    const auto *ki_151 = buffer.data(ki + 151);
    const auto *ki_152 = buffer.data(ki + 152);
    const auto *ki_153 = buffer.data(ki + 153);
    const auto *ki_154 = buffer.data(ki + 154);
    const auto *ki_155 = buffer.data(ki + 155);
    const auto *ki_156 = buffer.data(ki + 156);
    const auto *ki_157 = buffer.data(ki + 157);
    const auto *ki_158 = buffer.data(ki + 158);
    const auto *ki_159 = buffer.data(ki + 159);
    const auto *ki_160 = buffer.data(ki + 160);
    const auto *ki_161 = buffer.data(ki + 161);
    const auto *ki_162 = buffer.data(ki + 162);
    const auto *ki_163 = buffer.data(ki + 163);
    const auto *ki_164 = buffer.data(ki + 164);
    const auto *ki_165 = buffer.data(ki + 165);
    const auto *ki_166 = buffer.data(ki + 166);
    const auto *ki_167 = buffer.data(ki + 167);
    const auto *ki_168 = buffer.data(ki + 168);
    const auto *ki_169 = buffer.data(ki + 169);
    const auto *ki_170 = buffer.data(ki + 170);
    const auto *ki_171 = buffer.data(ki + 171);
    const auto *ki_172 = buffer.data(ki + 172);
    const auto *ki_173 = buffer.data(ki + 173);
    const auto *ki_174 = buffer.data(ki + 174);
    const auto *ki_175 = buffer.data(ki + 175);
    const auto *ki_176 = buffer.data(ki + 176);
    const auto *ki_177 = buffer.data(ki + 177);
    const auto *ki_178 = buffer.data(ki + 178);
    const auto *ki_179 = buffer.data(ki + 179);
    const auto *ki_180 = buffer.data(ki + 180);
    const auto *ki_181 = buffer.data(ki + 181);
    const auto *ki_182 = buffer.data(ki + 182);
    const auto *ki_183 = buffer.data(ki + 183);
    const auto *ki_184 = buffer.data(ki + 184);
    const auto *ki_185 = buffer.data(ki + 185);
    const auto *ki_186 = buffer.data(ki + 186);
    const auto *ki_187 = buffer.data(ki + 187);
    const auto *ki_188 = buffer.data(ki + 188);
    const auto *ki_189 = buffer.data(ki + 189);
    const auto *ki_190 = buffer.data(ki + 190);
    const auto *ki_191 = buffer.data(ki + 191);
    const auto *ki_192 = buffer.data(ki + 192);
    const auto *ki_193 = buffer.data(ki + 193);
    const auto *ki_194 = buffer.data(ki + 194);
    const auto *ki_195 = buffer.data(ki + 195);
    const auto *ki_196 = buffer.data(ki + 196);
    const auto *ki_197 = buffer.data(ki + 197);
    const auto *ki_198 = buffer.data(ki + 198);
    const auto *ki_199 = buffer.data(ki + 199);
    const auto *ki_200 = buffer.data(ki + 200);
    const auto *ki_201 = buffer.data(ki + 201);
    const auto *ki_202 = buffer.data(ki + 202);
    const auto *ki_203 = buffer.data(ki + 203);
    const auto *ki_204 = buffer.data(ki + 204);
    const auto *ki_205 = buffer.data(ki + 205);
    const auto *ki_206 = buffer.data(ki + 206);
    const auto *ki_207 = buffer.data(ki + 207);
    const auto *ki_208 = buffer.data(ki + 208);
    const auto *ki_209 = buffer.data(ki + 209);
    const auto *ki_210 = buffer.data(ki + 210);
    const auto *ki_211 = buffer.data(ki + 211);
    const auto *ki_212 = buffer.data(ki + 212);
    const auto *ki_213 = buffer.data(ki + 213);
    const auto *ki_214 = buffer.data(ki + 214);
    const auto *ki_215 = buffer.data(ki + 215);
    const auto *ki_216 = buffer.data(ki + 216);
    const auto *ki_217 = buffer.data(ki + 217);
    const auto *ki_218 = buffer.data(ki + 218);
    const auto *ki_219 = buffer.data(ki + 219);
    const auto *ki_220 = buffer.data(ki + 220);
    const auto *ki_221 = buffer.data(ki + 221);
    const auto *ki_222 = buffer.data(ki + 222);
    const auto *ki_223 = buffer.data(ki + 223);
    const auto *ki_224 = buffer.data(ki + 224);
    const auto *ki_225 = buffer.data(ki + 225);
    const auto *ki_226 = buffer.data(ki + 226);
    const auto *ki_227 = buffer.data(ki + 227);
    const auto *ki_228 = buffer.data(ki + 228);
    const auto *ki_229 = buffer.data(ki + 229);
    const auto *ki_230 = buffer.data(ki + 230);
    const auto *ki_231 = buffer.data(ki + 231);
    const auto *ki_232 = buffer.data(ki + 232);
    const auto *ki_233 = buffer.data(ki + 233);
    const auto *ki_234 = buffer.data(ki + 234);
    const auto *ki_235 = buffer.data(ki + 235);
    const auto *ki_236 = buffer.data(ki + 236);
    const auto *ki_237 = buffer.data(ki + 237);
    const auto *ki_238 = buffer.data(ki + 238);
    const auto *ki_239 = buffer.data(ki + 239);
    const auto *ki_240 = buffer.data(ki + 240);
    const auto *ki_241 = buffer.data(ki + 241);
    const auto *ki_242 = buffer.data(ki + 242);
    const auto *ki_243 = buffer.data(ki + 243);
    const auto *ki_244 = buffer.data(ki + 244);
    const auto *ki_245 = buffer.data(ki + 245);
    const auto *ki_246 = buffer.data(ki + 246);
    const auto *ki_247 = buffer.data(ki + 247);
    const auto *ki_248 = buffer.data(ki + 248);
    const auto *ki_249 = buffer.data(ki + 249);
    const auto *ki_250 = buffer.data(ki + 250);
    const auto *ki_251 = buffer.data(ki + 251);
    const auto *ki_252 = buffer.data(ki + 252);
    const auto *ki_253 = buffer.data(ki + 253);
    const auto *ki_254 = buffer.data(ki + 254);
    const auto *ki_255 = buffer.data(ki + 255);
    const auto *ki_256 = buffer.data(ki + 256);
    const auto *ki_257 = buffer.data(ki + 257);
    const auto *ki_258 = buffer.data(ki + 258);
    const auto *ki_259 = buffer.data(ki + 259);
    const auto *ki_260 = buffer.data(ki + 260);
    const auto *ki_261 = buffer.data(ki + 261);
    const auto *ki_262 = buffer.data(ki + 262);
    const auto *ki_263 = buffer.data(ki + 263);
    const auto *ki_264 = buffer.data(ki + 264);
    const auto *ki_265 = buffer.data(ki + 265);
    const auto *ki_266 = buffer.data(ki + 266);
    const auto *ki_267 = buffer.data(ki + 267);
    const auto *ki_268 = buffer.data(ki + 268);
    const auto *ki_269 = buffer.data(ki + 269);
    const auto *ki_270 = buffer.data(ki + 270);
    const auto *ki_271 = buffer.data(ki + 271);
    const auto *ki_272 = buffer.data(ki + 272);
    const auto *ki_273 = buffer.data(ki + 273);
    const auto *ki_274 = buffer.data(ki + 274);
    const auto *ki_275 = buffer.data(ki + 275);
    const auto *ki_276 = buffer.data(ki + 276);
    const auto *ki_277 = buffer.data(ki + 277);
    const auto *ki_278 = buffer.data(ki + 278);
    const auto *ki_279 = buffer.data(ki + 279);
    const auto *ki_280 = buffer.data(ki + 280);
    const auto *ki_281 = buffer.data(ki + 281);
    const auto *ki_282 = buffer.data(ki + 282);
    const auto *ki_283 = buffer.data(ki + 283);
    const auto *ki_284 = buffer.data(ki + 284);
    const auto *ki_285 = buffer.data(ki + 285);
    const auto *ki_286 = buffer.data(ki + 286);
    const auto *ki_287 = buffer.data(ki + 287);
    const auto *ki_288 = buffer.data(ki + 288);
    const auto *ki_289 = buffer.data(ki + 289);
    const auto *ki_290 = buffer.data(ki + 290);
    const auto *ki_291 = buffer.data(ki + 291);
    const auto *ki_292 = buffer.data(ki + 292);
    const auto *ki_293 = buffer.data(ki + 293);
    const auto *ki_294 = buffer.data(ki + 294);
    const auto *ki_295 = buffer.data(ki + 295);
    const auto *ki_296 = buffer.data(ki + 296);
    const auto *ki_297 = buffer.data(ki + 297);
    const auto *ki_298 = buffer.data(ki + 298);
    const auto *ki_299 = buffer.data(ki + 299);

    const auto *mi_150 = buffer.data(mi + 150);
    const auto *mi_151 = buffer.data(mi + 151);
    const auto *mi_152 = buffer.data(mi + 152);
    const auto *mi_153 = buffer.data(mi + 153);
    const auto *mi_154 = buffer.data(mi + 154);
    const auto *mi_155 = buffer.data(mi + 155);
    const auto *mi_156 = buffer.data(mi + 156);
    const auto *mi_157 = buffer.data(mi + 157);
    const auto *mi_158 = buffer.data(mi + 158);
    const auto *mi_159 = buffer.data(mi + 159);
    const auto *mi_160 = buffer.data(mi + 160);
    const auto *mi_161 = buffer.data(mi + 161);
    const auto *mi_162 = buffer.data(mi + 162);
    const auto *mi_163 = buffer.data(mi + 163);
    const auto *mi_164 = buffer.data(mi + 164);
    const auto *mi_165 = buffer.data(mi + 165);
    const auto *mi_166 = buffer.data(mi + 166);
    const auto *mi_167 = buffer.data(mi + 167);
    const auto *mi_168 = buffer.data(mi + 168);
    const auto *mi_169 = buffer.data(mi + 169);
    const auto *mi_170 = buffer.data(mi + 170);
    const auto *mi_171 = buffer.data(mi + 171);
    const auto *mi_172 = buffer.data(mi + 172);
    const auto *mi_173 = buffer.data(mi + 173);
    const auto *mi_174 = buffer.data(mi + 174);
    const auto *mi_175 = buffer.data(mi + 175);
    const auto *mi_176 = buffer.data(mi + 176);
    const auto *mi_177 = buffer.data(mi + 177);
    const auto *mi_178 = buffer.data(mi + 178);
    const auto *mi_179 = buffer.data(mi + 179);
    const auto *mi_180 = buffer.data(mi + 180);
    const auto *mi_181 = buffer.data(mi + 181);
    const auto *mi_182 = buffer.data(mi + 182);
    const auto *mi_183 = buffer.data(mi + 183);
    const auto *mi_184 = buffer.data(mi + 184);
    const auto *mi_185 = buffer.data(mi + 185);
    const auto *mi_186 = buffer.data(mi + 186);
    const auto *mi_187 = buffer.data(mi + 187);
    const auto *mi_188 = buffer.data(mi + 188);
    const auto *mi_189 = buffer.data(mi + 189);
    const auto *mi_190 = buffer.data(mi + 190);
    const auto *mi_191 = buffer.data(mi + 191);
    const auto *mi_192 = buffer.data(mi + 192);
    const auto *mi_193 = buffer.data(mi + 193);
    const auto *mi_194 = buffer.data(mi + 194);
    const auto *mi_195 = buffer.data(mi + 195);
    const auto *mi_196 = buffer.data(mi + 196);
    const auto *mi_197 = buffer.data(mi + 197);
    const auto *mi_198 = buffer.data(mi + 198);
    const auto *mi_199 = buffer.data(mi + 199);
    const auto *mi_200 = buffer.data(mi + 200);
    const auto *mi_201 = buffer.data(mi + 201);
    const auto *mi_202 = buffer.data(mi + 202);
    const auto *mi_203 = buffer.data(mi + 203);
    const auto *mi_204 = buffer.data(mi + 204);
    const auto *mi_205 = buffer.data(mi + 205);
    const auto *mi_206 = buffer.data(mi + 206);
    const auto *mi_207 = buffer.data(mi + 207);
    const auto *mi_208 = buffer.data(mi + 208);
    const auto *mi_209 = buffer.data(mi + 209);
    const auto *mi_210 = buffer.data(mi + 210);
    const auto *mi_211 = buffer.data(mi + 211);
    const auto *mi_212 = buffer.data(mi + 212);
    const auto *mi_213 = buffer.data(mi + 213);
    const auto *mi_214 = buffer.data(mi + 214);
    const auto *mi_215 = buffer.data(mi + 215);
    const auto *mi_216 = buffer.data(mi + 216);
    const auto *mi_217 = buffer.data(mi + 217);
    const auto *mi_218 = buffer.data(mi + 218);
    const auto *mi_219 = buffer.data(mi + 219);
    const auto *mi_220 = buffer.data(mi + 220);
    const auto *mi_221 = buffer.data(mi + 221);
    const auto *mi_222 = buffer.data(mi + 222);
    const auto *mi_223 = buffer.data(mi + 223);
    const auto *mi_224 = buffer.data(mi + 224);
    const auto *mi_225 = buffer.data(mi + 225);
    const auto *mi_226 = buffer.data(mi + 226);
    const auto *mi_227 = buffer.data(mi + 227);
    const auto *mi_228 = buffer.data(mi + 228);
    const auto *mi_229 = buffer.data(mi + 229);
    const auto *mi_230 = buffer.data(mi + 230);
    const auto *mi_231 = buffer.data(mi + 231);
    const auto *mi_232 = buffer.data(mi + 232);
    const auto *mi_233 = buffer.data(mi + 233);
    const auto *mi_234 = buffer.data(mi + 234);
    const auto *mi_235 = buffer.data(mi + 235);
    const auto *mi_236 = buffer.data(mi + 236);
    const auto *mi_237 = buffer.data(mi + 237);
    const auto *mi_238 = buffer.data(mi + 238);
    const auto *mi_239 = buffer.data(mi + 239);
    const auto *mi_240 = buffer.data(mi + 240);
    const auto *mi_241 = buffer.data(mi + 241);
    const auto *mi_242 = buffer.data(mi + 242);
    const auto *mi_243 = buffer.data(mi + 243);
    const auto *mi_244 = buffer.data(mi + 244);
    const auto *mi_245 = buffer.data(mi + 245);
    const auto *mi_246 = buffer.data(mi + 246);
    const auto *mi_247 = buffer.data(mi + 247);
    const auto *mi_248 = buffer.data(mi + 248);
    const auto *mi_249 = buffer.data(mi + 249);
    const auto *mi_250 = buffer.data(mi + 250);
    const auto *mi_251 = buffer.data(mi + 251);
    const auto *mi_252 = buffer.data(mi + 252);
    const auto *mi_253 = buffer.data(mi + 253);
    const auto *mi_254 = buffer.data(mi + 254);
    const auto *mi_255 = buffer.data(mi + 255);
    const auto *mi_256 = buffer.data(mi + 256);
    const auto *mi_257 = buffer.data(mi + 257);
    const auto *mi_258 = buffer.data(mi + 258);
    const auto *mi_259 = buffer.data(mi + 259);
    const auto *mi_260 = buffer.data(mi + 260);
    const auto *mi_261 = buffer.data(mi + 261);
    const auto *mi_262 = buffer.data(mi + 262);
    const auto *mi_263 = buffer.data(mi + 263);
    const auto *mi_264 = buffer.data(mi + 264);
    const auto *mi_265 = buffer.data(mi + 265);
    const auto *mi_266 = buffer.data(mi + 266);
    const auto *mi_267 = buffer.data(mi + 267);
    const auto *mi_268 = buffer.data(mi + 268);
    const auto *mi_269 = buffer.data(mi + 269);
    const auto *mi_270 = buffer.data(mi + 270);
    const auto *mi_271 = buffer.data(mi + 271);
    const auto *mi_272 = buffer.data(mi + 272);
    const auto *mi_273 = buffer.data(mi + 273);
    const auto *mi_274 = buffer.data(mi + 274);
    const auto *mi_275 = buffer.data(mi + 275);
    const auto *mi_276 = buffer.data(mi + 276);
    const auto *mi_277 = buffer.data(mi + 277);
    const auto *mi_278 = buffer.data(mi + 278);
    const auto *mi_279 = buffer.data(mi + 279);
    const auto *mi_280 = buffer.data(mi + 280);
    const auto *mi_281 = buffer.data(mi + 281);
    const auto *mi_282 = buffer.data(mi + 282);
    const auto *mi_283 = buffer.data(mi + 283);
    const auto *mi_284 = buffer.data(mi + 284);
    const auto *mi_285 = buffer.data(mi + 285);
    const auto *mi_286 = buffer.data(mi + 286);
    const auto *mi_287 = buffer.data(mi + 287);
    const auto *mi_288 = buffer.data(mi + 288);
    const auto *mi_289 = buffer.data(mi + 289);
    const auto *mi_290 = buffer.data(mi + 290);
    const auto *mi_291 = buffer.data(mi + 291);
    const auto *mi_292 = buffer.data(mi + 292);
    const auto *mi_293 = buffer.data(mi + 293);
    const auto *mi_294 = buffer.data(mi + 294);
    const auto *mi_295 = buffer.data(mi + 295);
    const auto *mi_296 = buffer.data(mi + 296);
    const auto *mi_297 = buffer.data(mi + 297);
    const auto *mi_298 = buffer.data(mi + 298);
    const auto *mi_299 = buffer.data(mi + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ki_150, ki_151, ki_152, ki_153, \
                         ki_154, mi_150, mi_151, mi_152, mi_153, \
                         mi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -6.0 * ki_150[k]
                   + f_0 * mi_150[k];

        t_151[k] = -6.0 * ki_151[k]
                   + f_0 * mi_151[k];

        t_152[k] = -6.0 * ki_152[k]
                   + f_0 * mi_152[k];

        t_153[k] = -6.0 * ki_153[k]
                   + f_0 * mi_153[k];

        t_154[k] = -6.0 * ki_154[k]
                   + f_0 * mi_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ki_155, ki_156, ki_157, ki_158, \
                         ki_159, mi_155, mi_156, mi_157, mi_158, \
                         mi_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -6.0 * ki_155[k]
                   + f_0 * mi_155[k];

        t_156[k] = -6.0 * ki_156[k]
                   + f_0 * mi_156[k];

        t_157[k] = -6.0 * ki_157[k]
                   + f_0 * mi_157[k];

        t_158[k] = -6.0 * ki_158[k]
                   + f_0 * mi_158[k];

        t_159[k] = -6.0 * ki_159[k]
                   + f_0 * mi_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ki_160, ki_161, ki_162, ki_163, \
                         ki_164, mi_160, mi_161, mi_162, mi_163, \
                         mi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -6.0 * ki_160[k]
                   + f_0 * mi_160[k];

        t_161[k] = -6.0 * ki_161[k]
                   + f_0 * mi_161[k];

        t_162[k] = -6.0 * ki_162[k]
                   + f_0 * mi_162[k];

        t_163[k] = -6.0 * ki_163[k]
                   + f_0 * mi_163[k];

        t_164[k] = -6.0 * ki_164[k]
                   + f_0 * mi_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ki_165, ki_166, ki_167, ki_168, \
                         ki_169, mi_165, mi_166, mi_167, mi_168, \
                         mi_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -6.0 * ki_165[k]
                   + f_0 * mi_165[k];

        t_166[k] = -6.0 * ki_166[k]
                   + f_0 * mi_166[k];

        t_167[k] = -6.0 * ki_167[k]
                   + f_0 * mi_167[k];

        t_168[k] = -5.0 * ki_168[k]
                   + f_0 * mi_168[k];

        t_169[k] = -5.0 * ki_169[k]
                   + f_0 * mi_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ki_170, ki_171, ki_172, ki_173, \
                         ki_174, mi_170, mi_171, mi_172, mi_173, \
                         mi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -5.0 * ki_170[k]
                   + f_0 * mi_170[k];

        t_171[k] = -5.0 * ki_171[k]
                   + f_0 * mi_171[k];

        t_172[k] = -5.0 * ki_172[k]
                   + f_0 * mi_172[k];

        t_173[k] = -5.0 * ki_173[k]
                   + f_0 * mi_173[k];

        t_174[k] = -5.0 * ki_174[k]
                   + f_0 * mi_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ki_175, ki_176, ki_177, ki_178, \
                         ki_179, mi_175, mi_176, mi_177, mi_178, \
                         mi_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -5.0 * ki_175[k]
                   + f_0 * mi_175[k];

        t_176[k] = -5.0 * ki_176[k]
                   + f_0 * mi_176[k];

        t_177[k] = -5.0 * ki_177[k]
                   + f_0 * mi_177[k];

        t_178[k] = -5.0 * ki_178[k]
                   + f_0 * mi_178[k];

        t_179[k] = -5.0 * ki_179[k]
                   + f_0 * mi_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ki_180, ki_181, ki_182, ki_183, \
                         ki_184, mi_180, mi_181, mi_182, mi_183, \
                         mi_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -5.0 * ki_180[k]
                   + f_0 * mi_180[k];

        t_181[k] = -5.0 * ki_181[k]
                   + f_0 * mi_181[k];

        t_182[k] = -5.0 * ki_182[k]
                   + f_0 * mi_182[k];

        t_183[k] = -5.0 * ki_183[k]
                   + f_0 * mi_183[k];

        t_184[k] = -5.0 * ki_184[k]
                   + f_0 * mi_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ki_185, ki_186, ki_187, ki_188, \
                         ki_189, mi_185, mi_186, mi_187, mi_188, \
                         mi_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -5.0 * ki_185[k]
                   + f_0 * mi_185[k];

        t_186[k] = -5.0 * ki_186[k]
                   + f_0 * mi_186[k];

        t_187[k] = -5.0 * ki_187[k]
                   + f_0 * mi_187[k];

        t_188[k] = -5.0 * ki_188[k]
                   + f_0 * mi_188[k];

        t_189[k] = -5.0 * ki_189[k]
                   + f_0 * mi_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ki_190, ki_191, ki_192, ki_193, \
                         ki_194, mi_190, mi_191, mi_192, mi_193, \
                         mi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -5.0 * ki_190[k]
                   + f_0 * mi_190[k];

        t_191[k] = -5.0 * ki_191[k]
                   + f_0 * mi_191[k];

        t_192[k] = -5.0 * ki_192[k]
                   + f_0 * mi_192[k];

        t_193[k] = -5.0 * ki_193[k]
                   + f_0 * mi_193[k];

        t_194[k] = -5.0 * ki_194[k]
                   + f_0 * mi_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ki_195, ki_196, ki_197, ki_198, \
                         ki_199, mi_195, mi_196, mi_197, mi_198, \
                         mi_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -5.0 * ki_195[k]
                   + f_0 * mi_195[k];

        t_196[k] = -5.0 * ki_196[k]
                   + f_0 * mi_196[k];

        t_197[k] = -5.0 * ki_197[k]
                   + f_0 * mi_197[k];

        t_198[k] = -5.0 * ki_198[k]
                   + f_0 * mi_198[k];

        t_199[k] = -5.0 * ki_199[k]
                   + f_0 * mi_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ki_200, ki_201, ki_202, ki_203, \
                         ki_204, mi_200, mi_201, mi_202, mi_203, \
                         mi_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -5.0 * ki_200[k]
                   + f_0 * mi_200[k];

        t_201[k] = -5.0 * ki_201[k]
                   + f_0 * mi_201[k];

        t_202[k] = -5.0 * ki_202[k]
                   + f_0 * mi_202[k];

        t_203[k] = -5.0 * ki_203[k]
                   + f_0 * mi_203[k];

        t_204[k] = -5.0 * ki_204[k]
                   + f_0 * mi_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ki_205, ki_206, ki_207, ki_208, \
                         ki_209, mi_205, mi_206, mi_207, mi_208, \
                         mi_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -5.0 * ki_205[k]
                   + f_0 * mi_205[k];

        t_206[k] = -5.0 * ki_206[k]
                   + f_0 * mi_206[k];

        t_207[k] = -5.0 * ki_207[k]
                   + f_0 * mi_207[k];

        t_208[k] = -5.0 * ki_208[k]
                   + f_0 * mi_208[k];

        t_209[k] = -5.0 * ki_209[k]
                   + f_0 * mi_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ki_210, ki_211, ki_212, ki_213, \
                         ki_214, mi_210, mi_211, mi_212, mi_213, \
                         mi_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -5.0 * ki_210[k]
                   + f_0 * mi_210[k];

        t_211[k] = -5.0 * ki_211[k]
                   + f_0 * mi_211[k];

        t_212[k] = -5.0 * ki_212[k]
                   + f_0 * mi_212[k];

        t_213[k] = -5.0 * ki_213[k]
                   + f_0 * mi_213[k];

        t_214[k] = -5.0 * ki_214[k]
                   + f_0 * mi_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ki_215, ki_216, ki_217, ki_218, \
                         ki_219, mi_215, mi_216, mi_217, mi_218, \
                         mi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -5.0 * ki_215[k]
                   + f_0 * mi_215[k];

        t_216[k] = -5.0 * ki_216[k]
                   + f_0 * mi_216[k];

        t_217[k] = -5.0 * ki_217[k]
                   + f_0 * mi_217[k];

        t_218[k] = -5.0 * ki_218[k]
                   + f_0 * mi_218[k];

        t_219[k] = -5.0 * ki_219[k]
                   + f_0 * mi_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ki_220, ki_221, ki_222, ki_223, \
                         ki_224, mi_220, mi_221, mi_222, mi_223, \
                         mi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -5.0 * ki_220[k]
                   + f_0 * mi_220[k];

        t_221[k] = -5.0 * ki_221[k]
                   + f_0 * mi_221[k];

        t_222[k] = -5.0 * ki_222[k]
                   + f_0 * mi_222[k];

        t_223[k] = -5.0 * ki_223[k]
                   + f_0 * mi_223[k];

        t_224[k] = -5.0 * ki_224[k]
                   + f_0 * mi_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ki_225, ki_226, ki_227, ki_228, \
                         ki_229, mi_225, mi_226, mi_227, mi_228, \
                         mi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -5.0 * ki_225[k]
                   + f_0 * mi_225[k];

        t_226[k] = -5.0 * ki_226[k]
                   + f_0 * mi_226[k];

        t_227[k] = -5.0 * ki_227[k]
                   + f_0 * mi_227[k];

        t_228[k] = -5.0 * ki_228[k]
                   + f_0 * mi_228[k];

        t_229[k] = -5.0 * ki_229[k]
                   + f_0 * mi_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ki_230, ki_231, ki_232, ki_233, \
                         ki_234, mi_230, mi_231, mi_232, mi_233, \
                         mi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -5.0 * ki_230[k]
                   + f_0 * mi_230[k];

        t_231[k] = -5.0 * ki_231[k]
                   + f_0 * mi_231[k];

        t_232[k] = -5.0 * ki_232[k]
                   + f_0 * mi_232[k];

        t_233[k] = -5.0 * ki_233[k]
                   + f_0 * mi_233[k];

        t_234[k] = -5.0 * ki_234[k]
                   + f_0 * mi_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ki_235, ki_236, ki_237, ki_238, \
                         ki_239, mi_235, mi_236, mi_237, mi_238, \
                         mi_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -5.0 * ki_235[k]
                   + f_0 * mi_235[k];

        t_236[k] = -5.0 * ki_236[k]
                   + f_0 * mi_236[k];

        t_237[k] = -5.0 * ki_237[k]
                   + f_0 * mi_237[k];

        t_238[k] = -5.0 * ki_238[k]
                   + f_0 * mi_238[k];

        t_239[k] = -5.0 * ki_239[k]
                   + f_0 * mi_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ki_240, ki_241, ki_242, ki_243, \
                         ki_244, mi_240, mi_241, mi_242, mi_243, \
                         mi_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -5.0 * ki_240[k]
                   + f_0 * mi_240[k];

        t_241[k] = -5.0 * ki_241[k]
                   + f_0 * mi_241[k];

        t_242[k] = -5.0 * ki_242[k]
                   + f_0 * mi_242[k];

        t_243[k] = -5.0 * ki_243[k]
                   + f_0 * mi_243[k];

        t_244[k] = -5.0 * ki_244[k]
                   + f_0 * mi_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ki_245, ki_246, ki_247, ki_248, \
                         ki_249, mi_245, mi_246, mi_247, mi_248, \
                         mi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -5.0 * ki_245[k]
                   + f_0 * mi_245[k];

        t_246[k] = -5.0 * ki_246[k]
                   + f_0 * mi_246[k];

        t_247[k] = -5.0 * ki_247[k]
                   + f_0 * mi_247[k];

        t_248[k] = -5.0 * ki_248[k]
                   + f_0 * mi_248[k];

        t_249[k] = -5.0 * ki_249[k]
                   + f_0 * mi_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ki_250, ki_251, ki_252, ki_253, \
                         ki_254, mi_250, mi_251, mi_252, mi_253, \
                         mi_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -5.0 * ki_250[k]
                   + f_0 * mi_250[k];

        t_251[k] = -5.0 * ki_251[k]
                   + f_0 * mi_251[k];

        t_252[k] = -5.0 * ki_252[k]
                   + f_0 * mi_252[k];

        t_253[k] = -5.0 * ki_253[k]
                   + f_0 * mi_253[k];

        t_254[k] = -5.0 * ki_254[k]
                   + f_0 * mi_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ki_255, ki_256, ki_257, ki_258, \
                         ki_259, mi_255, mi_256, mi_257, mi_258, \
                         mi_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -5.0 * ki_255[k]
                   + f_0 * mi_255[k];

        t_256[k] = -5.0 * ki_256[k]
                   + f_0 * mi_256[k];

        t_257[k] = -5.0 * ki_257[k]
                   + f_0 * mi_257[k];

        t_258[k] = -5.0 * ki_258[k]
                   + f_0 * mi_258[k];

        t_259[k] = -5.0 * ki_259[k]
                   + f_0 * mi_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ki_260, ki_261, ki_262, ki_263, \
                         ki_264, mi_260, mi_261, mi_262, mi_263, \
                         mi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -5.0 * ki_260[k]
                   + f_0 * mi_260[k];

        t_261[k] = -5.0 * ki_261[k]
                   + f_0 * mi_261[k];

        t_262[k] = -5.0 * ki_262[k]
                   + f_0 * mi_262[k];

        t_263[k] = -5.0 * ki_263[k]
                   + f_0 * mi_263[k];

        t_264[k] = -5.0 * ki_264[k]
                   + f_0 * mi_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ki_265, ki_266, ki_267, ki_268, \
                         ki_269, mi_265, mi_266, mi_267, mi_268, \
                         mi_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -5.0 * ki_265[k]
                   + f_0 * mi_265[k];

        t_266[k] = -5.0 * ki_266[k]
                   + f_0 * mi_266[k];

        t_267[k] = -5.0 * ki_267[k]
                   + f_0 * mi_267[k];

        t_268[k] = -5.0 * ki_268[k]
                   + f_0 * mi_268[k];

        t_269[k] = -5.0 * ki_269[k]
                   + f_0 * mi_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ki_270, ki_271, ki_272, ki_273, \
                         ki_274, mi_270, mi_271, mi_272, mi_273, \
                         mi_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -5.0 * ki_270[k]
                   + f_0 * mi_270[k];

        t_271[k] = -5.0 * ki_271[k]
                   + f_0 * mi_271[k];

        t_272[k] = -5.0 * ki_272[k]
                   + f_0 * mi_272[k];

        t_273[k] = -5.0 * ki_273[k]
                   + f_0 * mi_273[k];

        t_274[k] = -5.0 * ki_274[k]
                   + f_0 * mi_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ki_275, ki_276, ki_277, ki_278, \
                         ki_279, mi_275, mi_276, mi_277, mi_278, \
                         mi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -5.0 * ki_275[k]
                   + f_0 * mi_275[k];

        t_276[k] = -5.0 * ki_276[k]
                   + f_0 * mi_276[k];

        t_277[k] = -5.0 * ki_277[k]
                   + f_0 * mi_277[k];

        t_278[k] = -5.0 * ki_278[k]
                   + f_0 * mi_278[k];

        t_279[k] = -5.0 * ki_279[k]
                   + f_0 * mi_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ki_280, ki_281, ki_282, ki_283, \
                         ki_284, mi_280, mi_281, mi_282, mi_283, \
                         mi_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -4.0 * ki_280[k]
                   + f_0 * mi_280[k];

        t_281[k] = -4.0 * ki_281[k]
                   + f_0 * mi_281[k];

        t_282[k] = -4.0 * ki_282[k]
                   + f_0 * mi_282[k];

        t_283[k] = -4.0 * ki_283[k]
                   + f_0 * mi_283[k];

        t_284[k] = -4.0 * ki_284[k]
                   + f_0 * mi_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ki_285, ki_286, ki_287, ki_288, \
                         ki_289, mi_285, mi_286, mi_287, mi_288, \
                         mi_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -4.0 * ki_285[k]
                   + f_0 * mi_285[k];

        t_286[k] = -4.0 * ki_286[k]
                   + f_0 * mi_286[k];

        t_287[k] = -4.0 * ki_287[k]
                   + f_0 * mi_287[k];

        t_288[k] = -4.0 * ki_288[k]
                   + f_0 * mi_288[k];

        t_289[k] = -4.0 * ki_289[k]
                   + f_0 * mi_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ki_290, ki_291, ki_292, ki_293, \
                         ki_294, mi_290, mi_291, mi_292, mi_293, \
                         mi_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -4.0 * ki_290[k]
                   + f_0 * mi_290[k];

        t_291[k] = -4.0 * ki_291[k]
                   + f_0 * mi_291[k];

        t_292[k] = -4.0 * ki_292[k]
                   + f_0 * mi_292[k];

        t_293[k] = -4.0 * ki_293[k]
                   + f_0 * mi_293[k];

        t_294[k] = -4.0 * ki_294[k]
                   + f_0 * mi_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ki_295, ki_296, ki_297, ki_298, \
                         ki_299, mi_295, mi_296, mi_297, mi_298, \
                         mi_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * ki_295[k]
                   + f_0 * mi_295[k];

        t_296[k] = -4.0 * ki_296[k]
                   + f_0 * mi_296[k];

        t_297[k] = -4.0 * ki_297[k]
                   + f_0 * mi_297[k];

        t_298[k] = -4.0 * ki_298[k]
                   + f_0 * mi_298[k];

        t_299[k] = -4.0 * ki_299[k]
                   + f_0 * mi_299[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
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

    const auto *ki_300 = buffer.data(ki + 300);
    const auto *ki_301 = buffer.data(ki + 301);
    const auto *ki_302 = buffer.data(ki + 302);
    const auto *ki_303 = buffer.data(ki + 303);
    const auto *ki_304 = buffer.data(ki + 304);
    const auto *ki_305 = buffer.data(ki + 305);
    const auto *ki_306 = buffer.data(ki + 306);
    const auto *ki_307 = buffer.data(ki + 307);
    const auto *ki_308 = buffer.data(ki + 308);
    const auto *ki_309 = buffer.data(ki + 309);
    const auto *ki_310 = buffer.data(ki + 310);
    const auto *ki_311 = buffer.data(ki + 311);
    const auto *ki_312 = buffer.data(ki + 312);
    const auto *ki_313 = buffer.data(ki + 313);
    const auto *ki_314 = buffer.data(ki + 314);
    const auto *ki_315 = buffer.data(ki + 315);
    const auto *ki_316 = buffer.data(ki + 316);
    const auto *ki_317 = buffer.data(ki + 317);
    const auto *ki_318 = buffer.data(ki + 318);
    const auto *ki_319 = buffer.data(ki + 319);
    const auto *ki_320 = buffer.data(ki + 320);
    const auto *ki_321 = buffer.data(ki + 321);
    const auto *ki_322 = buffer.data(ki + 322);
    const auto *ki_323 = buffer.data(ki + 323);
    const auto *ki_324 = buffer.data(ki + 324);
    const auto *ki_325 = buffer.data(ki + 325);
    const auto *ki_326 = buffer.data(ki + 326);
    const auto *ki_327 = buffer.data(ki + 327);
    const auto *ki_328 = buffer.data(ki + 328);
    const auto *ki_329 = buffer.data(ki + 329);
    const auto *ki_330 = buffer.data(ki + 330);
    const auto *ki_331 = buffer.data(ki + 331);
    const auto *ki_332 = buffer.data(ki + 332);
    const auto *ki_333 = buffer.data(ki + 333);
    const auto *ki_334 = buffer.data(ki + 334);
    const auto *ki_335 = buffer.data(ki + 335);
    const auto *ki_336 = buffer.data(ki + 336);
    const auto *ki_337 = buffer.data(ki + 337);
    const auto *ki_338 = buffer.data(ki + 338);
    const auto *ki_339 = buffer.data(ki + 339);
    const auto *ki_340 = buffer.data(ki + 340);
    const auto *ki_341 = buffer.data(ki + 341);
    const auto *ki_342 = buffer.data(ki + 342);
    const auto *ki_343 = buffer.data(ki + 343);
    const auto *ki_344 = buffer.data(ki + 344);
    const auto *ki_345 = buffer.data(ki + 345);
    const auto *ki_346 = buffer.data(ki + 346);
    const auto *ki_347 = buffer.data(ki + 347);
    const auto *ki_348 = buffer.data(ki + 348);
    const auto *ki_349 = buffer.data(ki + 349);
    const auto *ki_350 = buffer.data(ki + 350);
    const auto *ki_351 = buffer.data(ki + 351);
    const auto *ki_352 = buffer.data(ki + 352);
    const auto *ki_353 = buffer.data(ki + 353);
    const auto *ki_354 = buffer.data(ki + 354);
    const auto *ki_355 = buffer.data(ki + 355);
    const auto *ki_356 = buffer.data(ki + 356);
    const auto *ki_357 = buffer.data(ki + 357);
    const auto *ki_358 = buffer.data(ki + 358);
    const auto *ki_359 = buffer.data(ki + 359);
    const auto *ki_360 = buffer.data(ki + 360);
    const auto *ki_361 = buffer.data(ki + 361);
    const auto *ki_362 = buffer.data(ki + 362);
    const auto *ki_363 = buffer.data(ki + 363);
    const auto *ki_364 = buffer.data(ki + 364);
    const auto *ki_365 = buffer.data(ki + 365);
    const auto *ki_366 = buffer.data(ki + 366);
    const auto *ki_367 = buffer.data(ki + 367);
    const auto *ki_368 = buffer.data(ki + 368);
    const auto *ki_369 = buffer.data(ki + 369);
    const auto *ki_370 = buffer.data(ki + 370);
    const auto *ki_371 = buffer.data(ki + 371);
    const auto *ki_372 = buffer.data(ki + 372);
    const auto *ki_373 = buffer.data(ki + 373);
    const auto *ki_374 = buffer.data(ki + 374);
    const auto *ki_375 = buffer.data(ki + 375);
    const auto *ki_376 = buffer.data(ki + 376);
    const auto *ki_377 = buffer.data(ki + 377);
    const auto *ki_378 = buffer.data(ki + 378);
    const auto *ki_379 = buffer.data(ki + 379);
    const auto *ki_380 = buffer.data(ki + 380);
    const auto *ki_381 = buffer.data(ki + 381);
    const auto *ki_382 = buffer.data(ki + 382);
    const auto *ki_383 = buffer.data(ki + 383);
    const auto *ki_384 = buffer.data(ki + 384);
    const auto *ki_385 = buffer.data(ki + 385);
    const auto *ki_386 = buffer.data(ki + 386);
    const auto *ki_387 = buffer.data(ki + 387);
    const auto *ki_388 = buffer.data(ki + 388);
    const auto *ki_389 = buffer.data(ki + 389);
    const auto *ki_390 = buffer.data(ki + 390);
    const auto *ki_391 = buffer.data(ki + 391);
    const auto *ki_392 = buffer.data(ki + 392);
    const auto *ki_393 = buffer.data(ki + 393);
    const auto *ki_394 = buffer.data(ki + 394);
    const auto *ki_395 = buffer.data(ki + 395);
    const auto *ki_396 = buffer.data(ki + 396);
    const auto *ki_397 = buffer.data(ki + 397);
    const auto *ki_398 = buffer.data(ki + 398);
    const auto *ki_399 = buffer.data(ki + 399);
    const auto *ki_400 = buffer.data(ki + 400);
    const auto *ki_401 = buffer.data(ki + 401);
    const auto *ki_402 = buffer.data(ki + 402);
    const auto *ki_403 = buffer.data(ki + 403);
    const auto *ki_404 = buffer.data(ki + 404);
    const auto *ki_405 = buffer.data(ki + 405);
    const auto *ki_406 = buffer.data(ki + 406);
    const auto *ki_407 = buffer.data(ki + 407);
    const auto *ki_408 = buffer.data(ki + 408);
    const auto *ki_409 = buffer.data(ki + 409);
    const auto *ki_410 = buffer.data(ki + 410);
    const auto *ki_411 = buffer.data(ki + 411);
    const auto *ki_412 = buffer.data(ki + 412);
    const auto *ki_413 = buffer.data(ki + 413);
    const auto *ki_414 = buffer.data(ki + 414);
    const auto *ki_415 = buffer.data(ki + 415);
    const auto *ki_416 = buffer.data(ki + 416);
    const auto *ki_417 = buffer.data(ki + 417);
    const auto *ki_418 = buffer.data(ki + 418);
    const auto *ki_419 = buffer.data(ki + 419);
    const auto *ki_420 = buffer.data(ki + 420);
    const auto *ki_421 = buffer.data(ki + 421);
    const auto *ki_422 = buffer.data(ki + 422);
    const auto *ki_423 = buffer.data(ki + 423);
    const auto *ki_424 = buffer.data(ki + 424);
    const auto *ki_425 = buffer.data(ki + 425);
    const auto *ki_426 = buffer.data(ki + 426);
    const auto *ki_427 = buffer.data(ki + 427);
    const auto *ki_428 = buffer.data(ki + 428);
    const auto *ki_429 = buffer.data(ki + 429);
    const auto *ki_430 = buffer.data(ki + 430);
    const auto *ki_431 = buffer.data(ki + 431);
    const auto *ki_432 = buffer.data(ki + 432);
    const auto *ki_433 = buffer.data(ki + 433);
    const auto *ki_434 = buffer.data(ki + 434);
    const auto *ki_435 = buffer.data(ki + 435);
    const auto *ki_436 = buffer.data(ki + 436);
    const auto *ki_437 = buffer.data(ki + 437);
    const auto *ki_438 = buffer.data(ki + 438);
    const auto *ki_439 = buffer.data(ki + 439);
    const auto *ki_440 = buffer.data(ki + 440);
    const auto *ki_441 = buffer.data(ki + 441);
    const auto *ki_442 = buffer.data(ki + 442);
    const auto *ki_443 = buffer.data(ki + 443);
    const auto *ki_444 = buffer.data(ki + 444);
    const auto *ki_445 = buffer.data(ki + 445);
    const auto *ki_446 = buffer.data(ki + 446);
    const auto *ki_447 = buffer.data(ki + 447);
    const auto *ki_448 = buffer.data(ki + 448);
    const auto *ki_449 = buffer.data(ki + 449);

    const auto *mi_300 = buffer.data(mi + 300);
    const auto *mi_301 = buffer.data(mi + 301);
    const auto *mi_302 = buffer.data(mi + 302);
    const auto *mi_303 = buffer.data(mi + 303);
    const auto *mi_304 = buffer.data(mi + 304);
    const auto *mi_305 = buffer.data(mi + 305);
    const auto *mi_306 = buffer.data(mi + 306);
    const auto *mi_307 = buffer.data(mi + 307);
    const auto *mi_308 = buffer.data(mi + 308);
    const auto *mi_309 = buffer.data(mi + 309);
    const auto *mi_310 = buffer.data(mi + 310);
    const auto *mi_311 = buffer.data(mi + 311);
    const auto *mi_312 = buffer.data(mi + 312);
    const auto *mi_313 = buffer.data(mi + 313);
    const auto *mi_314 = buffer.data(mi + 314);
    const auto *mi_315 = buffer.data(mi + 315);
    const auto *mi_316 = buffer.data(mi + 316);
    const auto *mi_317 = buffer.data(mi + 317);
    const auto *mi_318 = buffer.data(mi + 318);
    const auto *mi_319 = buffer.data(mi + 319);
    const auto *mi_320 = buffer.data(mi + 320);
    const auto *mi_321 = buffer.data(mi + 321);
    const auto *mi_322 = buffer.data(mi + 322);
    const auto *mi_323 = buffer.data(mi + 323);
    const auto *mi_324 = buffer.data(mi + 324);
    const auto *mi_325 = buffer.data(mi + 325);
    const auto *mi_326 = buffer.data(mi + 326);
    const auto *mi_327 = buffer.data(mi + 327);
    const auto *mi_328 = buffer.data(mi + 328);
    const auto *mi_329 = buffer.data(mi + 329);
    const auto *mi_330 = buffer.data(mi + 330);
    const auto *mi_331 = buffer.data(mi + 331);
    const auto *mi_332 = buffer.data(mi + 332);
    const auto *mi_333 = buffer.data(mi + 333);
    const auto *mi_334 = buffer.data(mi + 334);
    const auto *mi_335 = buffer.data(mi + 335);
    const auto *mi_336 = buffer.data(mi + 336);
    const auto *mi_337 = buffer.data(mi + 337);
    const auto *mi_338 = buffer.data(mi + 338);
    const auto *mi_339 = buffer.data(mi + 339);
    const auto *mi_340 = buffer.data(mi + 340);
    const auto *mi_341 = buffer.data(mi + 341);
    const auto *mi_342 = buffer.data(mi + 342);
    const auto *mi_343 = buffer.data(mi + 343);
    const auto *mi_344 = buffer.data(mi + 344);
    const auto *mi_345 = buffer.data(mi + 345);
    const auto *mi_346 = buffer.data(mi + 346);
    const auto *mi_347 = buffer.data(mi + 347);
    const auto *mi_348 = buffer.data(mi + 348);
    const auto *mi_349 = buffer.data(mi + 349);
    const auto *mi_350 = buffer.data(mi + 350);
    const auto *mi_351 = buffer.data(mi + 351);
    const auto *mi_352 = buffer.data(mi + 352);
    const auto *mi_353 = buffer.data(mi + 353);
    const auto *mi_354 = buffer.data(mi + 354);
    const auto *mi_355 = buffer.data(mi + 355);
    const auto *mi_356 = buffer.data(mi + 356);
    const auto *mi_357 = buffer.data(mi + 357);
    const auto *mi_358 = buffer.data(mi + 358);
    const auto *mi_359 = buffer.data(mi + 359);
    const auto *mi_360 = buffer.data(mi + 360);
    const auto *mi_361 = buffer.data(mi + 361);
    const auto *mi_362 = buffer.data(mi + 362);
    const auto *mi_363 = buffer.data(mi + 363);
    const auto *mi_364 = buffer.data(mi + 364);
    const auto *mi_365 = buffer.data(mi + 365);
    const auto *mi_366 = buffer.data(mi + 366);
    const auto *mi_367 = buffer.data(mi + 367);
    const auto *mi_368 = buffer.data(mi + 368);
    const auto *mi_369 = buffer.data(mi + 369);
    const auto *mi_370 = buffer.data(mi + 370);
    const auto *mi_371 = buffer.data(mi + 371);
    const auto *mi_372 = buffer.data(mi + 372);
    const auto *mi_373 = buffer.data(mi + 373);
    const auto *mi_374 = buffer.data(mi + 374);
    const auto *mi_375 = buffer.data(mi + 375);
    const auto *mi_376 = buffer.data(mi + 376);
    const auto *mi_377 = buffer.data(mi + 377);
    const auto *mi_378 = buffer.data(mi + 378);
    const auto *mi_379 = buffer.data(mi + 379);
    const auto *mi_380 = buffer.data(mi + 380);
    const auto *mi_381 = buffer.data(mi + 381);
    const auto *mi_382 = buffer.data(mi + 382);
    const auto *mi_383 = buffer.data(mi + 383);
    const auto *mi_384 = buffer.data(mi + 384);
    const auto *mi_385 = buffer.data(mi + 385);
    const auto *mi_386 = buffer.data(mi + 386);
    const auto *mi_387 = buffer.data(mi + 387);
    const auto *mi_388 = buffer.data(mi + 388);
    const auto *mi_389 = buffer.data(mi + 389);
    const auto *mi_390 = buffer.data(mi + 390);
    const auto *mi_391 = buffer.data(mi + 391);
    const auto *mi_392 = buffer.data(mi + 392);
    const auto *mi_393 = buffer.data(mi + 393);
    const auto *mi_394 = buffer.data(mi + 394);
    const auto *mi_395 = buffer.data(mi + 395);
    const auto *mi_396 = buffer.data(mi + 396);
    const auto *mi_397 = buffer.data(mi + 397);
    const auto *mi_398 = buffer.data(mi + 398);
    const auto *mi_399 = buffer.data(mi + 399);
    const auto *mi_400 = buffer.data(mi + 400);
    const auto *mi_401 = buffer.data(mi + 401);
    const auto *mi_402 = buffer.data(mi + 402);
    const auto *mi_403 = buffer.data(mi + 403);
    const auto *mi_404 = buffer.data(mi + 404);
    const auto *mi_405 = buffer.data(mi + 405);
    const auto *mi_406 = buffer.data(mi + 406);
    const auto *mi_407 = buffer.data(mi + 407);
    const auto *mi_408 = buffer.data(mi + 408);
    const auto *mi_409 = buffer.data(mi + 409);
    const auto *mi_410 = buffer.data(mi + 410);
    const auto *mi_411 = buffer.data(mi + 411);
    const auto *mi_412 = buffer.data(mi + 412);
    const auto *mi_413 = buffer.data(mi + 413);
    const auto *mi_414 = buffer.data(mi + 414);
    const auto *mi_415 = buffer.data(mi + 415);
    const auto *mi_416 = buffer.data(mi + 416);
    const auto *mi_417 = buffer.data(mi + 417);
    const auto *mi_418 = buffer.data(mi + 418);
    const auto *mi_419 = buffer.data(mi + 419);
    const auto *mi_420 = buffer.data(mi + 420);
    const auto *mi_421 = buffer.data(mi + 421);
    const auto *mi_422 = buffer.data(mi + 422);
    const auto *mi_423 = buffer.data(mi + 423);
    const auto *mi_424 = buffer.data(mi + 424);
    const auto *mi_425 = buffer.data(mi + 425);
    const auto *mi_426 = buffer.data(mi + 426);
    const auto *mi_427 = buffer.data(mi + 427);
    const auto *mi_428 = buffer.data(mi + 428);
    const auto *mi_429 = buffer.data(mi + 429);
    const auto *mi_430 = buffer.data(mi + 430);
    const auto *mi_431 = buffer.data(mi + 431);
    const auto *mi_432 = buffer.data(mi + 432);
    const auto *mi_433 = buffer.data(mi + 433);
    const auto *mi_434 = buffer.data(mi + 434);
    const auto *mi_435 = buffer.data(mi + 435);
    const auto *mi_436 = buffer.data(mi + 436);
    const auto *mi_437 = buffer.data(mi + 437);
    const auto *mi_438 = buffer.data(mi + 438);
    const auto *mi_439 = buffer.data(mi + 439);
    const auto *mi_440 = buffer.data(mi + 440);
    const auto *mi_441 = buffer.data(mi + 441);
    const auto *mi_442 = buffer.data(mi + 442);
    const auto *mi_443 = buffer.data(mi + 443);
    const auto *mi_444 = buffer.data(mi + 444);
    const auto *mi_445 = buffer.data(mi + 445);
    const auto *mi_446 = buffer.data(mi + 446);
    const auto *mi_447 = buffer.data(mi + 447);
    const auto *mi_448 = buffer.data(mi + 448);
    const auto *mi_449 = buffer.data(mi + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ki_300, ki_301, ki_302, ki_303, \
                         ki_304, mi_300, mi_301, mi_302, mi_303, \
                         mi_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * ki_300[k]
                   + f_0 * mi_300[k];

        t_301[k] = -4.0 * ki_301[k]
                   + f_0 * mi_301[k];

        t_302[k] = -4.0 * ki_302[k]
                   + f_0 * mi_302[k];

        t_303[k] = -4.0 * ki_303[k]
                   + f_0 * mi_303[k];

        t_304[k] = -4.0 * ki_304[k]
                   + f_0 * mi_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ki_305, ki_306, ki_307, ki_308, \
                         ki_309, mi_305, mi_306, mi_307, mi_308, \
                         mi_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * ki_305[k]
                   + f_0 * mi_305[k];

        t_306[k] = -4.0 * ki_306[k]
                   + f_0 * mi_306[k];

        t_307[k] = -4.0 * ki_307[k]
                   + f_0 * mi_307[k];

        t_308[k] = -4.0 * ki_308[k]
                   + f_0 * mi_308[k];

        t_309[k] = -4.0 * ki_309[k]
                   + f_0 * mi_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ki_310, ki_311, ki_312, ki_313, \
                         ki_314, mi_310, mi_311, mi_312, mi_313, \
                         mi_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -4.0 * ki_310[k]
                   + f_0 * mi_310[k];

        t_311[k] = -4.0 * ki_311[k]
                   + f_0 * mi_311[k];

        t_312[k] = -4.0 * ki_312[k]
                   + f_0 * mi_312[k];

        t_313[k] = -4.0 * ki_313[k]
                   + f_0 * mi_313[k];

        t_314[k] = -4.0 * ki_314[k]
                   + f_0 * mi_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ki_315, ki_316, ki_317, ki_318, \
                         ki_319, mi_315, mi_316, mi_317, mi_318, \
                         mi_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -4.0 * ki_315[k]
                   + f_0 * mi_315[k];

        t_316[k] = -4.0 * ki_316[k]
                   + f_0 * mi_316[k];

        t_317[k] = -4.0 * ki_317[k]
                   + f_0 * mi_317[k];

        t_318[k] = -4.0 * ki_318[k]
                   + f_0 * mi_318[k];

        t_319[k] = -4.0 * ki_319[k]
                   + f_0 * mi_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ki_320, ki_321, ki_322, ki_323, \
                         ki_324, mi_320, mi_321, mi_322, mi_323, \
                         mi_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -4.0 * ki_320[k]
                   + f_0 * mi_320[k];

        t_321[k] = -4.0 * ki_321[k]
                   + f_0 * mi_321[k];

        t_322[k] = -4.0 * ki_322[k]
                   + f_0 * mi_322[k];

        t_323[k] = -4.0 * ki_323[k]
                   + f_0 * mi_323[k];

        t_324[k] = -4.0 * ki_324[k]
                   + f_0 * mi_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ki_325, ki_326, ki_327, ki_328, \
                         ki_329, mi_325, mi_326, mi_327, mi_328, \
                         mi_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -4.0 * ki_325[k]
                   + f_0 * mi_325[k];

        t_326[k] = -4.0 * ki_326[k]
                   + f_0 * mi_326[k];

        t_327[k] = -4.0 * ki_327[k]
                   + f_0 * mi_327[k];

        t_328[k] = -4.0 * ki_328[k]
                   + f_0 * mi_328[k];

        t_329[k] = -4.0 * ki_329[k]
                   + f_0 * mi_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ki_330, ki_331, ki_332, ki_333, \
                         ki_334, mi_330, mi_331, mi_332, mi_333, \
                         mi_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -4.0 * ki_330[k]
                   + f_0 * mi_330[k];

        t_331[k] = -4.0 * ki_331[k]
                   + f_0 * mi_331[k];

        t_332[k] = -4.0 * ki_332[k]
                   + f_0 * mi_332[k];

        t_333[k] = -4.0 * ki_333[k]
                   + f_0 * mi_333[k];

        t_334[k] = -4.0 * ki_334[k]
                   + f_0 * mi_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ki_335, ki_336, ki_337, ki_338, \
                         ki_339, mi_335, mi_336, mi_337, mi_338, \
                         mi_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -4.0 * ki_335[k]
                   + f_0 * mi_335[k];

        t_336[k] = -4.0 * ki_336[k]
                   + f_0 * mi_336[k];

        t_337[k] = -4.0 * ki_337[k]
                   + f_0 * mi_337[k];

        t_338[k] = -4.0 * ki_338[k]
                   + f_0 * mi_338[k];

        t_339[k] = -4.0 * ki_339[k]
                   + f_0 * mi_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ki_340, ki_341, ki_342, ki_343, \
                         ki_344, mi_340, mi_341, mi_342, mi_343, \
                         mi_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -4.0 * ki_340[k]
                   + f_0 * mi_340[k];

        t_341[k] = -4.0 * ki_341[k]
                   + f_0 * mi_341[k];

        t_342[k] = -4.0 * ki_342[k]
                   + f_0 * mi_342[k];

        t_343[k] = -4.0 * ki_343[k]
                   + f_0 * mi_343[k];

        t_344[k] = -4.0 * ki_344[k]
                   + f_0 * mi_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ki_345, ki_346, ki_347, ki_348, \
                         ki_349, mi_345, mi_346, mi_347, mi_348, \
                         mi_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -4.0 * ki_345[k]
                   + f_0 * mi_345[k];

        t_346[k] = -4.0 * ki_346[k]
                   + f_0 * mi_346[k];

        t_347[k] = -4.0 * ki_347[k]
                   + f_0 * mi_347[k];

        t_348[k] = -4.0 * ki_348[k]
                   + f_0 * mi_348[k];

        t_349[k] = -4.0 * ki_349[k]
                   + f_0 * mi_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ki_350, ki_351, ki_352, ki_353, \
                         ki_354, mi_350, mi_351, mi_352, mi_353, \
                         mi_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -4.0 * ki_350[k]
                   + f_0 * mi_350[k];

        t_351[k] = -4.0 * ki_351[k]
                   + f_0 * mi_351[k];

        t_352[k] = -4.0 * ki_352[k]
                   + f_0 * mi_352[k];

        t_353[k] = -4.0 * ki_353[k]
                   + f_0 * mi_353[k];

        t_354[k] = -4.0 * ki_354[k]
                   + f_0 * mi_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ki_355, ki_356, ki_357, ki_358, \
                         ki_359, mi_355, mi_356, mi_357, mi_358, \
                         mi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -4.0 * ki_355[k]
                   + f_0 * mi_355[k];

        t_356[k] = -4.0 * ki_356[k]
                   + f_0 * mi_356[k];

        t_357[k] = -4.0 * ki_357[k]
                   + f_0 * mi_357[k];

        t_358[k] = -4.0 * ki_358[k]
                   + f_0 * mi_358[k];

        t_359[k] = -4.0 * ki_359[k]
                   + f_0 * mi_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ki_360, ki_361, ki_362, ki_363, \
                         ki_364, mi_360, mi_361, mi_362, mi_363, \
                         mi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -4.0 * ki_360[k]
                   + f_0 * mi_360[k];

        t_361[k] = -4.0 * ki_361[k]
                   + f_0 * mi_361[k];

        t_362[k] = -4.0 * ki_362[k]
                   + f_0 * mi_362[k];

        t_363[k] = -4.0 * ki_363[k]
                   + f_0 * mi_363[k];

        t_364[k] = -4.0 * ki_364[k]
                   + f_0 * mi_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ki_365, ki_366, ki_367, ki_368, \
                         ki_369, mi_365, mi_366, mi_367, mi_368, \
                         mi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -4.0 * ki_365[k]
                   + f_0 * mi_365[k];

        t_366[k] = -4.0 * ki_366[k]
                   + f_0 * mi_366[k];

        t_367[k] = -4.0 * ki_367[k]
                   + f_0 * mi_367[k];

        t_368[k] = -4.0 * ki_368[k]
                   + f_0 * mi_368[k];

        t_369[k] = -4.0 * ki_369[k]
                   + f_0 * mi_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ki_370, ki_371, ki_372, ki_373, \
                         ki_374, mi_370, mi_371, mi_372, mi_373, \
                         mi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -4.0 * ki_370[k]
                   + f_0 * mi_370[k];

        t_371[k] = -4.0 * ki_371[k]
                   + f_0 * mi_371[k];

        t_372[k] = -4.0 * ki_372[k]
                   + f_0 * mi_372[k];

        t_373[k] = -4.0 * ki_373[k]
                   + f_0 * mi_373[k];

        t_374[k] = -4.0 * ki_374[k]
                   + f_0 * mi_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ki_375, ki_376, ki_377, ki_378, \
                         ki_379, mi_375, mi_376, mi_377, mi_378, \
                         mi_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -4.0 * ki_375[k]
                   + f_0 * mi_375[k];

        t_376[k] = -4.0 * ki_376[k]
                   + f_0 * mi_376[k];

        t_377[k] = -4.0 * ki_377[k]
                   + f_0 * mi_377[k];

        t_378[k] = -4.0 * ki_378[k]
                   + f_0 * mi_378[k];

        t_379[k] = -4.0 * ki_379[k]
                   + f_0 * mi_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ki_380, ki_381, ki_382, ki_383, \
                         ki_384, mi_380, mi_381, mi_382, mi_383, \
                         mi_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -4.0 * ki_380[k]
                   + f_0 * mi_380[k];

        t_381[k] = -4.0 * ki_381[k]
                   + f_0 * mi_381[k];

        t_382[k] = -4.0 * ki_382[k]
                   + f_0 * mi_382[k];

        t_383[k] = -4.0 * ki_383[k]
                   + f_0 * mi_383[k];

        t_384[k] = -4.0 * ki_384[k]
                   + f_0 * mi_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ki_385, ki_386, ki_387, ki_388, \
                         ki_389, mi_385, mi_386, mi_387, mi_388, \
                         mi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -4.0 * ki_385[k]
                   + f_0 * mi_385[k];

        t_386[k] = -4.0 * ki_386[k]
                   + f_0 * mi_386[k];

        t_387[k] = -4.0 * ki_387[k]
                   + f_0 * mi_387[k];

        t_388[k] = -4.0 * ki_388[k]
                   + f_0 * mi_388[k];

        t_389[k] = -4.0 * ki_389[k]
                   + f_0 * mi_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ki_390, ki_391, ki_392, ki_393, \
                         ki_394, mi_390, mi_391, mi_392, mi_393, \
                         mi_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -4.0 * ki_390[k]
                   + f_0 * mi_390[k];

        t_391[k] = -4.0 * ki_391[k]
                   + f_0 * mi_391[k];

        t_392[k] = -4.0 * ki_392[k]
                   + f_0 * mi_392[k];

        t_393[k] = -4.0 * ki_393[k]
                   + f_0 * mi_393[k];

        t_394[k] = -4.0 * ki_394[k]
                   + f_0 * mi_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ki_395, ki_396, ki_397, ki_398, \
                         ki_399, mi_395, mi_396, mi_397, mi_398, \
                         mi_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -4.0 * ki_395[k]
                   + f_0 * mi_395[k];

        t_396[k] = -4.0 * ki_396[k]
                   + f_0 * mi_396[k];

        t_397[k] = -4.0 * ki_397[k]
                   + f_0 * mi_397[k];

        t_398[k] = -4.0 * ki_398[k]
                   + f_0 * mi_398[k];

        t_399[k] = -4.0 * ki_399[k]
                   + f_0 * mi_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ki_400, ki_401, ki_402, ki_403, \
                         ki_404, mi_400, mi_401, mi_402, mi_403, \
                         mi_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -4.0 * ki_400[k]
                   + f_0 * mi_400[k];

        t_401[k] = -4.0 * ki_401[k]
                   + f_0 * mi_401[k];

        t_402[k] = -4.0 * ki_402[k]
                   + f_0 * mi_402[k];

        t_403[k] = -4.0 * ki_403[k]
                   + f_0 * mi_403[k];

        t_404[k] = -4.0 * ki_404[k]
                   + f_0 * mi_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ki_405, ki_406, ki_407, ki_408, \
                         ki_409, mi_405, mi_406, mi_407, mi_408, \
                         mi_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -4.0 * ki_405[k]
                   + f_0 * mi_405[k];

        t_406[k] = -4.0 * ki_406[k]
                   + f_0 * mi_406[k];

        t_407[k] = -4.0 * ki_407[k]
                   + f_0 * mi_407[k];

        t_408[k] = -4.0 * ki_408[k]
                   + f_0 * mi_408[k];

        t_409[k] = -4.0 * ki_409[k]
                   + f_0 * mi_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ki_410, ki_411, ki_412, ki_413, \
                         ki_414, mi_410, mi_411, mi_412, mi_413, \
                         mi_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -4.0 * ki_410[k]
                   + f_0 * mi_410[k];

        t_411[k] = -4.0 * ki_411[k]
                   + f_0 * mi_411[k];

        t_412[k] = -4.0 * ki_412[k]
                   + f_0 * mi_412[k];

        t_413[k] = -4.0 * ki_413[k]
                   + f_0 * mi_413[k];

        t_414[k] = -4.0 * ki_414[k]
                   + f_0 * mi_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ki_415, ki_416, ki_417, ki_418, \
                         ki_419, mi_415, mi_416, mi_417, mi_418, \
                         mi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -4.0 * ki_415[k]
                   + f_0 * mi_415[k];

        t_416[k] = -4.0 * ki_416[k]
                   + f_0 * mi_416[k];

        t_417[k] = -4.0 * ki_417[k]
                   + f_0 * mi_417[k];

        t_418[k] = -4.0 * ki_418[k]
                   + f_0 * mi_418[k];

        t_419[k] = -4.0 * ki_419[k]
                   + f_0 * mi_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ki_420, ki_421, ki_422, ki_423, \
                         ki_424, mi_420, mi_421, mi_422, mi_423, \
                         mi_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -3.0 * ki_420[k]
                   + f_0 * mi_420[k];

        t_421[k] = -3.0 * ki_421[k]
                   + f_0 * mi_421[k];

        t_422[k] = -3.0 * ki_422[k]
                   + f_0 * mi_422[k];

        t_423[k] = -3.0 * ki_423[k]
                   + f_0 * mi_423[k];

        t_424[k] = -3.0 * ki_424[k]
                   + f_0 * mi_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ki_425, ki_426, ki_427, ki_428, \
                         ki_429, mi_425, mi_426, mi_427, mi_428, \
                         mi_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -3.0 * ki_425[k]
                   + f_0 * mi_425[k];

        t_426[k] = -3.0 * ki_426[k]
                   + f_0 * mi_426[k];

        t_427[k] = -3.0 * ki_427[k]
                   + f_0 * mi_427[k];

        t_428[k] = -3.0 * ki_428[k]
                   + f_0 * mi_428[k];

        t_429[k] = -3.0 * ki_429[k]
                   + f_0 * mi_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ki_430, ki_431, ki_432, ki_433, \
                         ki_434, mi_430, mi_431, mi_432, mi_433, \
                         mi_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -3.0 * ki_430[k]
                   + f_0 * mi_430[k];

        t_431[k] = -3.0 * ki_431[k]
                   + f_0 * mi_431[k];

        t_432[k] = -3.0 * ki_432[k]
                   + f_0 * mi_432[k];

        t_433[k] = -3.0 * ki_433[k]
                   + f_0 * mi_433[k];

        t_434[k] = -3.0 * ki_434[k]
                   + f_0 * mi_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ki_435, ki_436, ki_437, ki_438, \
                         ki_439, mi_435, mi_436, mi_437, mi_438, \
                         mi_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -3.0 * ki_435[k]
                   + f_0 * mi_435[k];

        t_436[k] = -3.0 * ki_436[k]
                   + f_0 * mi_436[k];

        t_437[k] = -3.0 * ki_437[k]
                   + f_0 * mi_437[k];

        t_438[k] = -3.0 * ki_438[k]
                   + f_0 * mi_438[k];

        t_439[k] = -3.0 * ki_439[k]
                   + f_0 * mi_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ki_440, ki_441, ki_442, ki_443, \
                         ki_444, mi_440, mi_441, mi_442, mi_443, \
                         mi_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -3.0 * ki_440[k]
                   + f_0 * mi_440[k];

        t_441[k] = -3.0 * ki_441[k]
                   + f_0 * mi_441[k];

        t_442[k] = -3.0 * ki_442[k]
                   + f_0 * mi_442[k];

        t_443[k] = -3.0 * ki_443[k]
                   + f_0 * mi_443[k];

        t_444[k] = -3.0 * ki_444[k]
                   + f_0 * mi_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ki_445, ki_446, ki_447, ki_448, \
                         ki_449, mi_445, mi_446, mi_447, mi_448, \
                         mi_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -3.0 * ki_445[k]
                   + f_0 * mi_445[k];

        t_446[k] = -3.0 * ki_446[k]
                   + f_0 * mi_446[k];

        t_447[k] = -3.0 * ki_447[k]
                   + f_0 * mi_447[k];

        t_448[k] = -3.0 * ki_448[k]
                   + f_0 * mi_448[k];

        t_449[k] = -3.0 * ki_449[k]
                   + f_0 * mi_449[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
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

    const auto *ki_450 = buffer.data(ki + 450);
    const auto *ki_451 = buffer.data(ki + 451);
    const auto *ki_452 = buffer.data(ki + 452);
    const auto *ki_453 = buffer.data(ki + 453);
    const auto *ki_454 = buffer.data(ki + 454);
    const auto *ki_455 = buffer.data(ki + 455);
    const auto *ki_456 = buffer.data(ki + 456);
    const auto *ki_457 = buffer.data(ki + 457);
    const auto *ki_458 = buffer.data(ki + 458);
    const auto *ki_459 = buffer.data(ki + 459);
    const auto *ki_460 = buffer.data(ki + 460);
    const auto *ki_461 = buffer.data(ki + 461);
    const auto *ki_462 = buffer.data(ki + 462);
    const auto *ki_463 = buffer.data(ki + 463);
    const auto *ki_464 = buffer.data(ki + 464);
    const auto *ki_465 = buffer.data(ki + 465);
    const auto *ki_466 = buffer.data(ki + 466);
    const auto *ki_467 = buffer.data(ki + 467);
    const auto *ki_468 = buffer.data(ki + 468);
    const auto *ki_469 = buffer.data(ki + 469);
    const auto *ki_470 = buffer.data(ki + 470);
    const auto *ki_471 = buffer.data(ki + 471);
    const auto *ki_472 = buffer.data(ki + 472);
    const auto *ki_473 = buffer.data(ki + 473);
    const auto *ki_474 = buffer.data(ki + 474);
    const auto *ki_475 = buffer.data(ki + 475);
    const auto *ki_476 = buffer.data(ki + 476);
    const auto *ki_477 = buffer.data(ki + 477);
    const auto *ki_478 = buffer.data(ki + 478);
    const auto *ki_479 = buffer.data(ki + 479);
    const auto *ki_480 = buffer.data(ki + 480);
    const auto *ki_481 = buffer.data(ki + 481);
    const auto *ki_482 = buffer.data(ki + 482);
    const auto *ki_483 = buffer.data(ki + 483);
    const auto *ki_484 = buffer.data(ki + 484);
    const auto *ki_485 = buffer.data(ki + 485);
    const auto *ki_486 = buffer.data(ki + 486);
    const auto *ki_487 = buffer.data(ki + 487);
    const auto *ki_488 = buffer.data(ki + 488);
    const auto *ki_489 = buffer.data(ki + 489);
    const auto *ki_490 = buffer.data(ki + 490);
    const auto *ki_491 = buffer.data(ki + 491);
    const auto *ki_492 = buffer.data(ki + 492);
    const auto *ki_493 = buffer.data(ki + 493);
    const auto *ki_494 = buffer.data(ki + 494);
    const auto *ki_495 = buffer.data(ki + 495);
    const auto *ki_496 = buffer.data(ki + 496);
    const auto *ki_497 = buffer.data(ki + 497);
    const auto *ki_498 = buffer.data(ki + 498);
    const auto *ki_499 = buffer.data(ki + 499);
    const auto *ki_500 = buffer.data(ki + 500);
    const auto *ki_501 = buffer.data(ki + 501);
    const auto *ki_502 = buffer.data(ki + 502);
    const auto *ki_503 = buffer.data(ki + 503);
    const auto *ki_504 = buffer.data(ki + 504);
    const auto *ki_505 = buffer.data(ki + 505);
    const auto *ki_506 = buffer.data(ki + 506);
    const auto *ki_507 = buffer.data(ki + 507);
    const auto *ki_508 = buffer.data(ki + 508);
    const auto *ki_509 = buffer.data(ki + 509);
    const auto *ki_510 = buffer.data(ki + 510);
    const auto *ki_511 = buffer.data(ki + 511);
    const auto *ki_512 = buffer.data(ki + 512);
    const auto *ki_513 = buffer.data(ki + 513);
    const auto *ki_514 = buffer.data(ki + 514);
    const auto *ki_515 = buffer.data(ki + 515);
    const auto *ki_516 = buffer.data(ki + 516);
    const auto *ki_517 = buffer.data(ki + 517);
    const auto *ki_518 = buffer.data(ki + 518);
    const auto *ki_519 = buffer.data(ki + 519);
    const auto *ki_520 = buffer.data(ki + 520);
    const auto *ki_521 = buffer.data(ki + 521);
    const auto *ki_522 = buffer.data(ki + 522);
    const auto *ki_523 = buffer.data(ki + 523);
    const auto *ki_524 = buffer.data(ki + 524);
    const auto *ki_525 = buffer.data(ki + 525);
    const auto *ki_526 = buffer.data(ki + 526);
    const auto *ki_527 = buffer.data(ki + 527);
    const auto *ki_528 = buffer.data(ki + 528);
    const auto *ki_529 = buffer.data(ki + 529);
    const auto *ki_530 = buffer.data(ki + 530);
    const auto *ki_531 = buffer.data(ki + 531);
    const auto *ki_532 = buffer.data(ki + 532);
    const auto *ki_533 = buffer.data(ki + 533);
    const auto *ki_534 = buffer.data(ki + 534);
    const auto *ki_535 = buffer.data(ki + 535);
    const auto *ki_536 = buffer.data(ki + 536);
    const auto *ki_537 = buffer.data(ki + 537);
    const auto *ki_538 = buffer.data(ki + 538);
    const auto *ki_539 = buffer.data(ki + 539);
    const auto *ki_540 = buffer.data(ki + 540);
    const auto *ki_541 = buffer.data(ki + 541);
    const auto *ki_542 = buffer.data(ki + 542);
    const auto *ki_543 = buffer.data(ki + 543);
    const auto *ki_544 = buffer.data(ki + 544);
    const auto *ki_545 = buffer.data(ki + 545);
    const auto *ki_546 = buffer.data(ki + 546);
    const auto *ki_547 = buffer.data(ki + 547);
    const auto *ki_548 = buffer.data(ki + 548);
    const auto *ki_549 = buffer.data(ki + 549);
    const auto *ki_550 = buffer.data(ki + 550);
    const auto *ki_551 = buffer.data(ki + 551);
    const auto *ki_552 = buffer.data(ki + 552);
    const auto *ki_553 = buffer.data(ki + 553);
    const auto *ki_554 = buffer.data(ki + 554);
    const auto *ki_555 = buffer.data(ki + 555);
    const auto *ki_556 = buffer.data(ki + 556);
    const auto *ki_557 = buffer.data(ki + 557);
    const auto *ki_558 = buffer.data(ki + 558);
    const auto *ki_559 = buffer.data(ki + 559);
    const auto *ki_560 = buffer.data(ki + 560);
    const auto *ki_561 = buffer.data(ki + 561);
    const auto *ki_562 = buffer.data(ki + 562);
    const auto *ki_563 = buffer.data(ki + 563);
    const auto *ki_564 = buffer.data(ki + 564);
    const auto *ki_565 = buffer.data(ki + 565);
    const auto *ki_566 = buffer.data(ki + 566);
    const auto *ki_567 = buffer.data(ki + 567);
    const auto *ki_568 = buffer.data(ki + 568);
    const auto *ki_569 = buffer.data(ki + 569);
    const auto *ki_570 = buffer.data(ki + 570);
    const auto *ki_571 = buffer.data(ki + 571);
    const auto *ki_572 = buffer.data(ki + 572);
    const auto *ki_573 = buffer.data(ki + 573);
    const auto *ki_574 = buffer.data(ki + 574);
    const auto *ki_575 = buffer.data(ki + 575);
    const auto *ki_576 = buffer.data(ki + 576);
    const auto *ki_577 = buffer.data(ki + 577);
    const auto *ki_578 = buffer.data(ki + 578);
    const auto *ki_579 = buffer.data(ki + 579);
    const auto *ki_580 = buffer.data(ki + 580);
    const auto *ki_581 = buffer.data(ki + 581);
    const auto *ki_582 = buffer.data(ki + 582);
    const auto *ki_583 = buffer.data(ki + 583);
    const auto *ki_584 = buffer.data(ki + 584);
    const auto *ki_585 = buffer.data(ki + 585);
    const auto *ki_586 = buffer.data(ki + 586);
    const auto *ki_587 = buffer.data(ki + 587);
    const auto *ki_588 = buffer.data(ki + 588);
    const auto *ki_589 = buffer.data(ki + 589);
    const auto *ki_590 = buffer.data(ki + 590);
    const auto *ki_591 = buffer.data(ki + 591);
    const auto *ki_592 = buffer.data(ki + 592);
    const auto *ki_593 = buffer.data(ki + 593);
    const auto *ki_594 = buffer.data(ki + 594);
    const auto *ki_595 = buffer.data(ki + 595);
    const auto *ki_596 = buffer.data(ki + 596);
    const auto *ki_597 = buffer.data(ki + 597);
    const auto *ki_598 = buffer.data(ki + 598);
    const auto *ki_599 = buffer.data(ki + 599);

    const auto *mi_450 = buffer.data(mi + 450);
    const auto *mi_451 = buffer.data(mi + 451);
    const auto *mi_452 = buffer.data(mi + 452);
    const auto *mi_453 = buffer.data(mi + 453);
    const auto *mi_454 = buffer.data(mi + 454);
    const auto *mi_455 = buffer.data(mi + 455);
    const auto *mi_456 = buffer.data(mi + 456);
    const auto *mi_457 = buffer.data(mi + 457);
    const auto *mi_458 = buffer.data(mi + 458);
    const auto *mi_459 = buffer.data(mi + 459);
    const auto *mi_460 = buffer.data(mi + 460);
    const auto *mi_461 = buffer.data(mi + 461);
    const auto *mi_462 = buffer.data(mi + 462);
    const auto *mi_463 = buffer.data(mi + 463);
    const auto *mi_464 = buffer.data(mi + 464);
    const auto *mi_465 = buffer.data(mi + 465);
    const auto *mi_466 = buffer.data(mi + 466);
    const auto *mi_467 = buffer.data(mi + 467);
    const auto *mi_468 = buffer.data(mi + 468);
    const auto *mi_469 = buffer.data(mi + 469);
    const auto *mi_470 = buffer.data(mi + 470);
    const auto *mi_471 = buffer.data(mi + 471);
    const auto *mi_472 = buffer.data(mi + 472);
    const auto *mi_473 = buffer.data(mi + 473);
    const auto *mi_474 = buffer.data(mi + 474);
    const auto *mi_475 = buffer.data(mi + 475);
    const auto *mi_476 = buffer.data(mi + 476);
    const auto *mi_477 = buffer.data(mi + 477);
    const auto *mi_478 = buffer.data(mi + 478);
    const auto *mi_479 = buffer.data(mi + 479);
    const auto *mi_480 = buffer.data(mi + 480);
    const auto *mi_481 = buffer.data(mi + 481);
    const auto *mi_482 = buffer.data(mi + 482);
    const auto *mi_483 = buffer.data(mi + 483);
    const auto *mi_484 = buffer.data(mi + 484);
    const auto *mi_485 = buffer.data(mi + 485);
    const auto *mi_486 = buffer.data(mi + 486);
    const auto *mi_487 = buffer.data(mi + 487);
    const auto *mi_488 = buffer.data(mi + 488);
    const auto *mi_489 = buffer.data(mi + 489);
    const auto *mi_490 = buffer.data(mi + 490);
    const auto *mi_491 = buffer.data(mi + 491);
    const auto *mi_492 = buffer.data(mi + 492);
    const auto *mi_493 = buffer.data(mi + 493);
    const auto *mi_494 = buffer.data(mi + 494);
    const auto *mi_495 = buffer.data(mi + 495);
    const auto *mi_496 = buffer.data(mi + 496);
    const auto *mi_497 = buffer.data(mi + 497);
    const auto *mi_498 = buffer.data(mi + 498);
    const auto *mi_499 = buffer.data(mi + 499);
    const auto *mi_500 = buffer.data(mi + 500);
    const auto *mi_501 = buffer.data(mi + 501);
    const auto *mi_502 = buffer.data(mi + 502);
    const auto *mi_503 = buffer.data(mi + 503);
    const auto *mi_504 = buffer.data(mi + 504);
    const auto *mi_505 = buffer.data(mi + 505);
    const auto *mi_506 = buffer.data(mi + 506);
    const auto *mi_507 = buffer.data(mi + 507);
    const auto *mi_508 = buffer.data(mi + 508);
    const auto *mi_509 = buffer.data(mi + 509);
    const auto *mi_510 = buffer.data(mi + 510);
    const auto *mi_511 = buffer.data(mi + 511);
    const auto *mi_512 = buffer.data(mi + 512);
    const auto *mi_513 = buffer.data(mi + 513);
    const auto *mi_514 = buffer.data(mi + 514);
    const auto *mi_515 = buffer.data(mi + 515);
    const auto *mi_516 = buffer.data(mi + 516);
    const auto *mi_517 = buffer.data(mi + 517);
    const auto *mi_518 = buffer.data(mi + 518);
    const auto *mi_519 = buffer.data(mi + 519);
    const auto *mi_520 = buffer.data(mi + 520);
    const auto *mi_521 = buffer.data(mi + 521);
    const auto *mi_522 = buffer.data(mi + 522);
    const auto *mi_523 = buffer.data(mi + 523);
    const auto *mi_524 = buffer.data(mi + 524);
    const auto *mi_525 = buffer.data(mi + 525);
    const auto *mi_526 = buffer.data(mi + 526);
    const auto *mi_527 = buffer.data(mi + 527);
    const auto *mi_528 = buffer.data(mi + 528);
    const auto *mi_529 = buffer.data(mi + 529);
    const auto *mi_530 = buffer.data(mi + 530);
    const auto *mi_531 = buffer.data(mi + 531);
    const auto *mi_532 = buffer.data(mi + 532);
    const auto *mi_533 = buffer.data(mi + 533);
    const auto *mi_534 = buffer.data(mi + 534);
    const auto *mi_535 = buffer.data(mi + 535);
    const auto *mi_536 = buffer.data(mi + 536);
    const auto *mi_537 = buffer.data(mi + 537);
    const auto *mi_538 = buffer.data(mi + 538);
    const auto *mi_539 = buffer.data(mi + 539);
    const auto *mi_540 = buffer.data(mi + 540);
    const auto *mi_541 = buffer.data(mi + 541);
    const auto *mi_542 = buffer.data(mi + 542);
    const auto *mi_543 = buffer.data(mi + 543);
    const auto *mi_544 = buffer.data(mi + 544);
    const auto *mi_545 = buffer.data(mi + 545);
    const auto *mi_546 = buffer.data(mi + 546);
    const auto *mi_547 = buffer.data(mi + 547);
    const auto *mi_548 = buffer.data(mi + 548);
    const auto *mi_549 = buffer.data(mi + 549);
    const auto *mi_550 = buffer.data(mi + 550);
    const auto *mi_551 = buffer.data(mi + 551);
    const auto *mi_552 = buffer.data(mi + 552);
    const auto *mi_553 = buffer.data(mi + 553);
    const auto *mi_554 = buffer.data(mi + 554);
    const auto *mi_555 = buffer.data(mi + 555);
    const auto *mi_556 = buffer.data(mi + 556);
    const auto *mi_557 = buffer.data(mi + 557);
    const auto *mi_558 = buffer.data(mi + 558);
    const auto *mi_559 = buffer.data(mi + 559);
    const auto *mi_560 = buffer.data(mi + 560);
    const auto *mi_561 = buffer.data(mi + 561);
    const auto *mi_562 = buffer.data(mi + 562);
    const auto *mi_563 = buffer.data(mi + 563);
    const auto *mi_564 = buffer.data(mi + 564);
    const auto *mi_565 = buffer.data(mi + 565);
    const auto *mi_566 = buffer.data(mi + 566);
    const auto *mi_567 = buffer.data(mi + 567);
    const auto *mi_568 = buffer.data(mi + 568);
    const auto *mi_569 = buffer.data(mi + 569);
    const auto *mi_570 = buffer.data(mi + 570);
    const auto *mi_571 = buffer.data(mi + 571);
    const auto *mi_572 = buffer.data(mi + 572);
    const auto *mi_573 = buffer.data(mi + 573);
    const auto *mi_574 = buffer.data(mi + 574);
    const auto *mi_575 = buffer.data(mi + 575);
    const auto *mi_576 = buffer.data(mi + 576);
    const auto *mi_577 = buffer.data(mi + 577);
    const auto *mi_578 = buffer.data(mi + 578);
    const auto *mi_579 = buffer.data(mi + 579);
    const auto *mi_580 = buffer.data(mi + 580);
    const auto *mi_581 = buffer.data(mi + 581);
    const auto *mi_582 = buffer.data(mi + 582);
    const auto *mi_583 = buffer.data(mi + 583);
    const auto *mi_584 = buffer.data(mi + 584);
    const auto *mi_585 = buffer.data(mi + 585);
    const auto *mi_586 = buffer.data(mi + 586);
    const auto *mi_587 = buffer.data(mi + 587);
    const auto *mi_588 = buffer.data(mi + 588);
    const auto *mi_589 = buffer.data(mi + 589);
    const auto *mi_590 = buffer.data(mi + 590);
    const auto *mi_591 = buffer.data(mi + 591);
    const auto *mi_592 = buffer.data(mi + 592);
    const auto *mi_593 = buffer.data(mi + 593);
    const auto *mi_594 = buffer.data(mi + 594);
    const auto *mi_595 = buffer.data(mi + 595);
    const auto *mi_596 = buffer.data(mi + 596);
    const auto *mi_597 = buffer.data(mi + 597);
    const auto *mi_598 = buffer.data(mi + 598);
    const auto *mi_599 = buffer.data(mi + 599);

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ki_450, ki_451, ki_452, ki_453, \
                         ki_454, mi_450, mi_451, mi_452, mi_453, \
                         mi_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -3.0 * ki_450[k]
                   + f_0 * mi_450[k];

        t_451[k] = -3.0 * ki_451[k]
                   + f_0 * mi_451[k];

        t_452[k] = -3.0 * ki_452[k]
                   + f_0 * mi_452[k];

        t_453[k] = -3.0 * ki_453[k]
                   + f_0 * mi_453[k];

        t_454[k] = -3.0 * ki_454[k]
                   + f_0 * mi_454[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ki_455, ki_456, ki_457, ki_458, \
                         ki_459, mi_455, mi_456, mi_457, mi_458, \
                         mi_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -3.0 * ki_455[k]
                   + f_0 * mi_455[k];

        t_456[k] = -3.0 * ki_456[k]
                   + f_0 * mi_456[k];

        t_457[k] = -3.0 * ki_457[k]
                   + f_0 * mi_457[k];

        t_458[k] = -3.0 * ki_458[k]
                   + f_0 * mi_458[k];

        t_459[k] = -3.0 * ki_459[k]
                   + f_0 * mi_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ki_460, ki_461, ki_462, ki_463, \
                         ki_464, mi_460, mi_461, mi_462, mi_463, \
                         mi_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -3.0 * ki_460[k]
                   + f_0 * mi_460[k];

        t_461[k] = -3.0 * ki_461[k]
                   + f_0 * mi_461[k];

        t_462[k] = -3.0 * ki_462[k]
                   + f_0 * mi_462[k];

        t_463[k] = -3.0 * ki_463[k]
                   + f_0 * mi_463[k];

        t_464[k] = -3.0 * ki_464[k]
                   + f_0 * mi_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ki_465, ki_466, ki_467, ki_468, \
                         ki_469, mi_465, mi_466, mi_467, mi_468, \
                         mi_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -3.0 * ki_465[k]
                   + f_0 * mi_465[k];

        t_466[k] = -3.0 * ki_466[k]
                   + f_0 * mi_466[k];

        t_467[k] = -3.0 * ki_467[k]
                   + f_0 * mi_467[k];

        t_468[k] = -3.0 * ki_468[k]
                   + f_0 * mi_468[k];

        t_469[k] = -3.0 * ki_469[k]
                   + f_0 * mi_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ki_470, ki_471, ki_472, ki_473, \
                         ki_474, mi_470, mi_471, mi_472, mi_473, \
                         mi_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -3.0 * ki_470[k]
                   + f_0 * mi_470[k];

        t_471[k] = -3.0 * ki_471[k]
                   + f_0 * mi_471[k];

        t_472[k] = -3.0 * ki_472[k]
                   + f_0 * mi_472[k];

        t_473[k] = -3.0 * ki_473[k]
                   + f_0 * mi_473[k];

        t_474[k] = -3.0 * ki_474[k]
                   + f_0 * mi_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ki_475, ki_476, ki_477, ki_478, \
                         ki_479, mi_475, mi_476, mi_477, mi_478, \
                         mi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -3.0 * ki_475[k]
                   + f_0 * mi_475[k];

        t_476[k] = -3.0 * ki_476[k]
                   + f_0 * mi_476[k];

        t_477[k] = -3.0 * ki_477[k]
                   + f_0 * mi_477[k];

        t_478[k] = -3.0 * ki_478[k]
                   + f_0 * mi_478[k];

        t_479[k] = -3.0 * ki_479[k]
                   + f_0 * mi_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ki_480, ki_481, ki_482, ki_483, \
                         ki_484, mi_480, mi_481, mi_482, mi_483, \
                         mi_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -3.0 * ki_480[k]
                   + f_0 * mi_480[k];

        t_481[k] = -3.0 * ki_481[k]
                   + f_0 * mi_481[k];

        t_482[k] = -3.0 * ki_482[k]
                   + f_0 * mi_482[k];

        t_483[k] = -3.0 * ki_483[k]
                   + f_0 * mi_483[k];

        t_484[k] = -3.0 * ki_484[k]
                   + f_0 * mi_484[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ki_485, ki_486, ki_487, ki_488, \
                         ki_489, mi_485, mi_486, mi_487, mi_488, \
                         mi_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -3.0 * ki_485[k]
                   + f_0 * mi_485[k];

        t_486[k] = -3.0 * ki_486[k]
                   + f_0 * mi_486[k];

        t_487[k] = -3.0 * ki_487[k]
                   + f_0 * mi_487[k];

        t_488[k] = -3.0 * ki_488[k]
                   + f_0 * mi_488[k];

        t_489[k] = -3.0 * ki_489[k]
                   + f_0 * mi_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ki_490, ki_491, ki_492, ki_493, \
                         ki_494, mi_490, mi_491, mi_492, mi_493, \
                         mi_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -3.0 * ki_490[k]
                   + f_0 * mi_490[k];

        t_491[k] = -3.0 * ki_491[k]
                   + f_0 * mi_491[k];

        t_492[k] = -3.0 * ki_492[k]
                   + f_0 * mi_492[k];

        t_493[k] = -3.0 * ki_493[k]
                   + f_0 * mi_493[k];

        t_494[k] = -3.0 * ki_494[k]
                   + f_0 * mi_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ki_495, ki_496, ki_497, ki_498, \
                         ki_499, mi_495, mi_496, mi_497, mi_498, \
                         mi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -3.0 * ki_495[k]
                   + f_0 * mi_495[k];

        t_496[k] = -3.0 * ki_496[k]
                   + f_0 * mi_496[k];

        t_497[k] = -3.0 * ki_497[k]
                   + f_0 * mi_497[k];

        t_498[k] = -3.0 * ki_498[k]
                   + f_0 * mi_498[k];

        t_499[k] = -3.0 * ki_499[k]
                   + f_0 * mi_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ki_500, ki_501, ki_502, ki_503, \
                         ki_504, mi_500, mi_501, mi_502, mi_503, \
                         mi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -3.0 * ki_500[k]
                   + f_0 * mi_500[k];

        t_501[k] = -3.0 * ki_501[k]
                   + f_0 * mi_501[k];

        t_502[k] = -3.0 * ki_502[k]
                   + f_0 * mi_502[k];

        t_503[k] = -3.0 * ki_503[k]
                   + f_0 * mi_503[k];

        t_504[k] = -3.0 * ki_504[k]
                   + f_0 * mi_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ki_505, ki_506, ki_507, ki_508, \
                         ki_509, mi_505, mi_506, mi_507, mi_508, \
                         mi_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -3.0 * ki_505[k]
                   + f_0 * mi_505[k];

        t_506[k] = -3.0 * ki_506[k]
                   + f_0 * mi_506[k];

        t_507[k] = -3.0 * ki_507[k]
                   + f_0 * mi_507[k];

        t_508[k] = -3.0 * ki_508[k]
                   + f_0 * mi_508[k];

        t_509[k] = -3.0 * ki_509[k]
                   + f_0 * mi_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ki_510, ki_511, ki_512, ki_513, \
                         ki_514, mi_510, mi_511, mi_512, mi_513, \
                         mi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -3.0 * ki_510[k]
                   + f_0 * mi_510[k];

        t_511[k] = -3.0 * ki_511[k]
                   + f_0 * mi_511[k];

        t_512[k] = -3.0 * ki_512[k]
                   + f_0 * mi_512[k];

        t_513[k] = -3.0 * ki_513[k]
                   + f_0 * mi_513[k];

        t_514[k] = -3.0 * ki_514[k]
                   + f_0 * mi_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ki_515, ki_516, ki_517, ki_518, \
                         ki_519, mi_515, mi_516, mi_517, mi_518, \
                         mi_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -3.0 * ki_515[k]
                   + f_0 * mi_515[k];

        t_516[k] = -3.0 * ki_516[k]
                   + f_0 * mi_516[k];

        t_517[k] = -3.0 * ki_517[k]
                   + f_0 * mi_517[k];

        t_518[k] = -3.0 * ki_518[k]
                   + f_0 * mi_518[k];

        t_519[k] = -3.0 * ki_519[k]
                   + f_0 * mi_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ki_520, ki_521, ki_522, ki_523, \
                         ki_524, mi_520, mi_521, mi_522, mi_523, \
                         mi_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -3.0 * ki_520[k]
                   + f_0 * mi_520[k];

        t_521[k] = -3.0 * ki_521[k]
                   + f_0 * mi_521[k];

        t_522[k] = -3.0 * ki_522[k]
                   + f_0 * mi_522[k];

        t_523[k] = -3.0 * ki_523[k]
                   + f_0 * mi_523[k];

        t_524[k] = -3.0 * ki_524[k]
                   + f_0 * mi_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ki_525, ki_526, ki_527, ki_528, \
                         ki_529, mi_525, mi_526, mi_527, mi_528, \
                         mi_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -3.0 * ki_525[k]
                   + f_0 * mi_525[k];

        t_526[k] = -3.0 * ki_526[k]
                   + f_0 * mi_526[k];

        t_527[k] = -3.0 * ki_527[k]
                   + f_0 * mi_527[k];

        t_528[k] = -3.0 * ki_528[k]
                   + f_0 * mi_528[k];

        t_529[k] = -3.0 * ki_529[k]
                   + f_0 * mi_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ki_530, ki_531, ki_532, ki_533, \
                         ki_534, mi_530, mi_531, mi_532, mi_533, \
                         mi_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -3.0 * ki_530[k]
                   + f_0 * mi_530[k];

        t_531[k] = -3.0 * ki_531[k]
                   + f_0 * mi_531[k];

        t_532[k] = -3.0 * ki_532[k]
                   + f_0 * mi_532[k];

        t_533[k] = -3.0 * ki_533[k]
                   + f_0 * mi_533[k];

        t_534[k] = -3.0 * ki_534[k]
                   + f_0 * mi_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ki_535, ki_536, ki_537, ki_538, \
                         ki_539, mi_535, mi_536, mi_537, mi_538, \
                         mi_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -3.0 * ki_535[k]
                   + f_0 * mi_535[k];

        t_536[k] = -3.0 * ki_536[k]
                   + f_0 * mi_536[k];

        t_537[k] = -3.0 * ki_537[k]
                   + f_0 * mi_537[k];

        t_538[k] = -3.0 * ki_538[k]
                   + f_0 * mi_538[k];

        t_539[k] = -3.0 * ki_539[k]
                   + f_0 * mi_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ki_540, ki_541, ki_542, ki_543, \
                         ki_544, mi_540, mi_541, mi_542, mi_543, \
                         mi_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -3.0 * ki_540[k]
                   + f_0 * mi_540[k];

        t_541[k] = -3.0 * ki_541[k]
                   + f_0 * mi_541[k];

        t_542[k] = -3.0 * ki_542[k]
                   + f_0 * mi_542[k];

        t_543[k] = -3.0 * ki_543[k]
                   + f_0 * mi_543[k];

        t_544[k] = -3.0 * ki_544[k]
                   + f_0 * mi_544[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ki_545, ki_546, ki_547, ki_548, \
                         ki_549, mi_545, mi_546, mi_547, mi_548, \
                         mi_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -3.0 * ki_545[k]
                   + f_0 * mi_545[k];

        t_546[k] = -3.0 * ki_546[k]
                   + f_0 * mi_546[k];

        t_547[k] = -3.0 * ki_547[k]
                   + f_0 * mi_547[k];

        t_548[k] = -3.0 * ki_548[k]
                   + f_0 * mi_548[k];

        t_549[k] = -3.0 * ki_549[k]
                   + f_0 * mi_549[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ki_550, ki_551, ki_552, ki_553, \
                         ki_554, mi_550, mi_551, mi_552, mi_553, \
                         mi_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -3.0 * ki_550[k]
                   + f_0 * mi_550[k];

        t_551[k] = -3.0 * ki_551[k]
                   + f_0 * mi_551[k];

        t_552[k] = -3.0 * ki_552[k]
                   + f_0 * mi_552[k];

        t_553[k] = -3.0 * ki_553[k]
                   + f_0 * mi_553[k];

        t_554[k] = -3.0 * ki_554[k]
                   + f_0 * mi_554[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ki_555, ki_556, ki_557, ki_558, \
                         ki_559, mi_555, mi_556, mi_557, mi_558, \
                         mi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -3.0 * ki_555[k]
                   + f_0 * mi_555[k];

        t_556[k] = -3.0 * ki_556[k]
                   + f_0 * mi_556[k];

        t_557[k] = -3.0 * ki_557[k]
                   + f_0 * mi_557[k];

        t_558[k] = -3.0 * ki_558[k]
                   + f_0 * mi_558[k];

        t_559[k] = -3.0 * ki_559[k]
                   + f_0 * mi_559[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, ki_560, ki_561, ki_562, ki_563, \
                         ki_564, mi_560, mi_561, mi_562, mi_563, \
                         mi_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = -3.0 * ki_560[k]
                   + f_0 * mi_560[k];

        t_561[k] = -3.0 * ki_561[k]
                   + f_0 * mi_561[k];

        t_562[k] = -3.0 * ki_562[k]
                   + f_0 * mi_562[k];

        t_563[k] = -3.0 * ki_563[k]
                   + f_0 * mi_563[k];

        t_564[k] = -3.0 * ki_564[k]
                   + f_0 * mi_564[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, ki_565, ki_566, ki_567, ki_568, \
                         ki_569, mi_565, mi_566, mi_567, mi_568, \
                         mi_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = -3.0 * ki_565[k]
                   + f_0 * mi_565[k];

        t_566[k] = -3.0 * ki_566[k]
                   + f_0 * mi_566[k];

        t_567[k] = -3.0 * ki_567[k]
                   + f_0 * mi_567[k];

        t_568[k] = -3.0 * ki_568[k]
                   + f_0 * mi_568[k];

        t_569[k] = -3.0 * ki_569[k]
                   + f_0 * mi_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ki_570, ki_571, ki_572, ki_573, \
                         ki_574, mi_570, mi_571, mi_572, mi_573, \
                         mi_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -3.0 * ki_570[k]
                   + f_0 * mi_570[k];

        t_571[k] = -3.0 * ki_571[k]
                   + f_0 * mi_571[k];

        t_572[k] = -3.0 * ki_572[k]
                   + f_0 * mi_572[k];

        t_573[k] = -3.0 * ki_573[k]
                   + f_0 * mi_573[k];

        t_574[k] = -3.0 * ki_574[k]
                   + f_0 * mi_574[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ki_575, ki_576, ki_577, ki_578, \
                         ki_579, mi_575, mi_576, mi_577, mi_578, \
                         mi_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = -3.0 * ki_575[k]
                   + f_0 * mi_575[k];

        t_576[k] = -3.0 * ki_576[k]
                   + f_0 * mi_576[k];

        t_577[k] = -3.0 * ki_577[k]
                   + f_0 * mi_577[k];

        t_578[k] = -3.0 * ki_578[k]
                   + f_0 * mi_578[k];

        t_579[k] = -3.0 * ki_579[k]
                   + f_0 * mi_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, ki_580, ki_581, ki_582, ki_583, \
                         ki_584, mi_580, mi_581, mi_582, mi_583, \
                         mi_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -3.0 * ki_580[k]
                   + f_0 * mi_580[k];

        t_581[k] = -3.0 * ki_581[k]
                   + f_0 * mi_581[k];

        t_582[k] = -3.0 * ki_582[k]
                   + f_0 * mi_582[k];

        t_583[k] = -3.0 * ki_583[k]
                   + f_0 * mi_583[k];

        t_584[k] = -3.0 * ki_584[k]
                   + f_0 * mi_584[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, ki_585, ki_586, ki_587, ki_588, \
                         ki_589, mi_585, mi_586, mi_587, mi_588, \
                         mi_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -3.0 * ki_585[k]
                   + f_0 * mi_585[k];

        t_586[k] = -3.0 * ki_586[k]
                   + f_0 * mi_586[k];

        t_587[k] = -3.0 * ki_587[k]
                   + f_0 * mi_587[k];

        t_588[k] = -2.0 * ki_588[k]
                   + f_0 * mi_588[k];

        t_589[k] = -2.0 * ki_589[k]
                   + f_0 * mi_589[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, ki_590, ki_591, ki_592, ki_593, \
                         ki_594, mi_590, mi_591, mi_592, mi_593, \
                         mi_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -2.0 * ki_590[k]
                   + f_0 * mi_590[k];

        t_591[k] = -2.0 * ki_591[k]
                   + f_0 * mi_591[k];

        t_592[k] = -2.0 * ki_592[k]
                   + f_0 * mi_592[k];

        t_593[k] = -2.0 * ki_593[k]
                   + f_0 * mi_593[k];

        t_594[k] = -2.0 * ki_594[k]
                   + f_0 * mi_594[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, ki_595, ki_596, ki_597, ki_598, \
                         ki_599, mi_595, mi_596, mi_597, mi_598, \
                         mi_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -2.0 * ki_595[k]
                   + f_0 * mi_595[k];

        t_596[k] = -2.0 * ki_596[k]
                   + f_0 * mi_596[k];

        t_597[k] = -2.0 * ki_597[k]
                   + f_0 * mi_597[k];

        t_598[k] = -2.0 * ki_598[k]
                   + f_0 * mi_598[k];

        t_599[k] = -2.0 * ki_599[k]
                   + f_0 * mi_599[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
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

    const auto *ki_600 = buffer.data(ki + 600);
    const auto *ki_601 = buffer.data(ki + 601);
    const auto *ki_602 = buffer.data(ki + 602);
    const auto *ki_603 = buffer.data(ki + 603);
    const auto *ki_604 = buffer.data(ki + 604);
    const auto *ki_605 = buffer.data(ki + 605);
    const auto *ki_606 = buffer.data(ki + 606);
    const auto *ki_607 = buffer.data(ki + 607);
    const auto *ki_608 = buffer.data(ki + 608);
    const auto *ki_609 = buffer.data(ki + 609);
    const auto *ki_610 = buffer.data(ki + 610);
    const auto *ki_611 = buffer.data(ki + 611);
    const auto *ki_612 = buffer.data(ki + 612);
    const auto *ki_613 = buffer.data(ki + 613);
    const auto *ki_614 = buffer.data(ki + 614);
    const auto *ki_615 = buffer.data(ki + 615);
    const auto *ki_616 = buffer.data(ki + 616);
    const auto *ki_617 = buffer.data(ki + 617);
    const auto *ki_618 = buffer.data(ki + 618);
    const auto *ki_619 = buffer.data(ki + 619);
    const auto *ki_620 = buffer.data(ki + 620);
    const auto *ki_621 = buffer.data(ki + 621);
    const auto *ki_622 = buffer.data(ki + 622);
    const auto *ki_623 = buffer.data(ki + 623);
    const auto *ki_624 = buffer.data(ki + 624);
    const auto *ki_625 = buffer.data(ki + 625);
    const auto *ki_626 = buffer.data(ki + 626);
    const auto *ki_627 = buffer.data(ki + 627);
    const auto *ki_628 = buffer.data(ki + 628);
    const auto *ki_629 = buffer.data(ki + 629);
    const auto *ki_630 = buffer.data(ki + 630);
    const auto *ki_631 = buffer.data(ki + 631);
    const auto *ki_632 = buffer.data(ki + 632);
    const auto *ki_633 = buffer.data(ki + 633);
    const auto *ki_634 = buffer.data(ki + 634);
    const auto *ki_635 = buffer.data(ki + 635);
    const auto *ki_636 = buffer.data(ki + 636);
    const auto *ki_637 = buffer.data(ki + 637);
    const auto *ki_638 = buffer.data(ki + 638);
    const auto *ki_639 = buffer.data(ki + 639);
    const auto *ki_640 = buffer.data(ki + 640);
    const auto *ki_641 = buffer.data(ki + 641);
    const auto *ki_642 = buffer.data(ki + 642);
    const auto *ki_643 = buffer.data(ki + 643);
    const auto *ki_644 = buffer.data(ki + 644);
    const auto *ki_645 = buffer.data(ki + 645);
    const auto *ki_646 = buffer.data(ki + 646);
    const auto *ki_647 = buffer.data(ki + 647);
    const auto *ki_648 = buffer.data(ki + 648);
    const auto *ki_649 = buffer.data(ki + 649);
    const auto *ki_650 = buffer.data(ki + 650);
    const auto *ki_651 = buffer.data(ki + 651);
    const auto *ki_652 = buffer.data(ki + 652);
    const auto *ki_653 = buffer.data(ki + 653);
    const auto *ki_654 = buffer.data(ki + 654);
    const auto *ki_655 = buffer.data(ki + 655);
    const auto *ki_656 = buffer.data(ki + 656);
    const auto *ki_657 = buffer.data(ki + 657);
    const auto *ki_658 = buffer.data(ki + 658);
    const auto *ki_659 = buffer.data(ki + 659);
    const auto *ki_660 = buffer.data(ki + 660);
    const auto *ki_661 = buffer.data(ki + 661);
    const auto *ki_662 = buffer.data(ki + 662);
    const auto *ki_663 = buffer.data(ki + 663);
    const auto *ki_664 = buffer.data(ki + 664);
    const auto *ki_665 = buffer.data(ki + 665);
    const auto *ki_666 = buffer.data(ki + 666);
    const auto *ki_667 = buffer.data(ki + 667);
    const auto *ki_668 = buffer.data(ki + 668);
    const auto *ki_669 = buffer.data(ki + 669);
    const auto *ki_670 = buffer.data(ki + 670);
    const auto *ki_671 = buffer.data(ki + 671);
    const auto *ki_672 = buffer.data(ki + 672);
    const auto *ki_673 = buffer.data(ki + 673);
    const auto *ki_674 = buffer.data(ki + 674);
    const auto *ki_675 = buffer.data(ki + 675);
    const auto *ki_676 = buffer.data(ki + 676);
    const auto *ki_677 = buffer.data(ki + 677);
    const auto *ki_678 = buffer.data(ki + 678);
    const auto *ki_679 = buffer.data(ki + 679);
    const auto *ki_680 = buffer.data(ki + 680);
    const auto *ki_681 = buffer.data(ki + 681);
    const auto *ki_682 = buffer.data(ki + 682);
    const auto *ki_683 = buffer.data(ki + 683);
    const auto *ki_684 = buffer.data(ki + 684);
    const auto *ki_685 = buffer.data(ki + 685);
    const auto *ki_686 = buffer.data(ki + 686);
    const auto *ki_687 = buffer.data(ki + 687);
    const auto *ki_688 = buffer.data(ki + 688);
    const auto *ki_689 = buffer.data(ki + 689);
    const auto *ki_690 = buffer.data(ki + 690);
    const auto *ki_691 = buffer.data(ki + 691);
    const auto *ki_692 = buffer.data(ki + 692);
    const auto *ki_693 = buffer.data(ki + 693);
    const auto *ki_694 = buffer.data(ki + 694);
    const auto *ki_695 = buffer.data(ki + 695);
    const auto *ki_696 = buffer.data(ki + 696);
    const auto *ki_697 = buffer.data(ki + 697);
    const auto *ki_698 = buffer.data(ki + 698);
    const auto *ki_699 = buffer.data(ki + 699);
    const auto *ki_700 = buffer.data(ki + 700);
    const auto *ki_701 = buffer.data(ki + 701);
    const auto *ki_702 = buffer.data(ki + 702);
    const auto *ki_703 = buffer.data(ki + 703);
    const auto *ki_704 = buffer.data(ki + 704);
    const auto *ki_705 = buffer.data(ki + 705);
    const auto *ki_706 = buffer.data(ki + 706);
    const auto *ki_707 = buffer.data(ki + 707);
    const auto *ki_708 = buffer.data(ki + 708);
    const auto *ki_709 = buffer.data(ki + 709);
    const auto *ki_710 = buffer.data(ki + 710);
    const auto *ki_711 = buffer.data(ki + 711);
    const auto *ki_712 = buffer.data(ki + 712);
    const auto *ki_713 = buffer.data(ki + 713);
    const auto *ki_714 = buffer.data(ki + 714);
    const auto *ki_715 = buffer.data(ki + 715);
    const auto *ki_716 = buffer.data(ki + 716);
    const auto *ki_717 = buffer.data(ki + 717);
    const auto *ki_718 = buffer.data(ki + 718);
    const auto *ki_719 = buffer.data(ki + 719);
    const auto *ki_720 = buffer.data(ki + 720);
    const auto *ki_721 = buffer.data(ki + 721);
    const auto *ki_722 = buffer.data(ki + 722);
    const auto *ki_723 = buffer.data(ki + 723);
    const auto *ki_724 = buffer.data(ki + 724);
    const auto *ki_725 = buffer.data(ki + 725);
    const auto *ki_726 = buffer.data(ki + 726);
    const auto *ki_727 = buffer.data(ki + 727);
    const auto *ki_728 = buffer.data(ki + 728);
    const auto *ki_729 = buffer.data(ki + 729);
    const auto *ki_730 = buffer.data(ki + 730);
    const auto *ki_731 = buffer.data(ki + 731);
    const auto *ki_732 = buffer.data(ki + 732);
    const auto *ki_733 = buffer.data(ki + 733);
    const auto *ki_734 = buffer.data(ki + 734);
    const auto *ki_735 = buffer.data(ki + 735);
    const auto *ki_736 = buffer.data(ki + 736);
    const auto *ki_737 = buffer.data(ki + 737);
    const auto *ki_738 = buffer.data(ki + 738);
    const auto *ki_739 = buffer.data(ki + 739);
    const auto *ki_740 = buffer.data(ki + 740);
    const auto *ki_741 = buffer.data(ki + 741);
    const auto *ki_742 = buffer.data(ki + 742);
    const auto *ki_743 = buffer.data(ki + 743);
    const auto *ki_744 = buffer.data(ki + 744);
    const auto *ki_745 = buffer.data(ki + 745);
    const auto *ki_746 = buffer.data(ki + 746);
    const auto *ki_747 = buffer.data(ki + 747);
    const auto *ki_748 = buffer.data(ki + 748);
    const auto *ki_749 = buffer.data(ki + 749);

    const auto *mi_600 = buffer.data(mi + 600);
    const auto *mi_601 = buffer.data(mi + 601);
    const auto *mi_602 = buffer.data(mi + 602);
    const auto *mi_603 = buffer.data(mi + 603);
    const auto *mi_604 = buffer.data(mi + 604);
    const auto *mi_605 = buffer.data(mi + 605);
    const auto *mi_606 = buffer.data(mi + 606);
    const auto *mi_607 = buffer.data(mi + 607);
    const auto *mi_608 = buffer.data(mi + 608);
    const auto *mi_609 = buffer.data(mi + 609);
    const auto *mi_610 = buffer.data(mi + 610);
    const auto *mi_611 = buffer.data(mi + 611);
    const auto *mi_612 = buffer.data(mi + 612);
    const auto *mi_613 = buffer.data(mi + 613);
    const auto *mi_614 = buffer.data(mi + 614);
    const auto *mi_615 = buffer.data(mi + 615);
    const auto *mi_616 = buffer.data(mi + 616);
    const auto *mi_617 = buffer.data(mi + 617);
    const auto *mi_618 = buffer.data(mi + 618);
    const auto *mi_619 = buffer.data(mi + 619);
    const auto *mi_620 = buffer.data(mi + 620);
    const auto *mi_621 = buffer.data(mi + 621);
    const auto *mi_622 = buffer.data(mi + 622);
    const auto *mi_623 = buffer.data(mi + 623);
    const auto *mi_624 = buffer.data(mi + 624);
    const auto *mi_625 = buffer.data(mi + 625);
    const auto *mi_626 = buffer.data(mi + 626);
    const auto *mi_627 = buffer.data(mi + 627);
    const auto *mi_628 = buffer.data(mi + 628);
    const auto *mi_629 = buffer.data(mi + 629);
    const auto *mi_630 = buffer.data(mi + 630);
    const auto *mi_631 = buffer.data(mi + 631);
    const auto *mi_632 = buffer.data(mi + 632);
    const auto *mi_633 = buffer.data(mi + 633);
    const auto *mi_634 = buffer.data(mi + 634);
    const auto *mi_635 = buffer.data(mi + 635);
    const auto *mi_636 = buffer.data(mi + 636);
    const auto *mi_637 = buffer.data(mi + 637);
    const auto *mi_638 = buffer.data(mi + 638);
    const auto *mi_639 = buffer.data(mi + 639);
    const auto *mi_640 = buffer.data(mi + 640);
    const auto *mi_641 = buffer.data(mi + 641);
    const auto *mi_642 = buffer.data(mi + 642);
    const auto *mi_643 = buffer.data(mi + 643);
    const auto *mi_644 = buffer.data(mi + 644);
    const auto *mi_645 = buffer.data(mi + 645);
    const auto *mi_646 = buffer.data(mi + 646);
    const auto *mi_647 = buffer.data(mi + 647);
    const auto *mi_648 = buffer.data(mi + 648);
    const auto *mi_649 = buffer.data(mi + 649);
    const auto *mi_650 = buffer.data(mi + 650);
    const auto *mi_651 = buffer.data(mi + 651);
    const auto *mi_652 = buffer.data(mi + 652);
    const auto *mi_653 = buffer.data(mi + 653);
    const auto *mi_654 = buffer.data(mi + 654);
    const auto *mi_655 = buffer.data(mi + 655);
    const auto *mi_656 = buffer.data(mi + 656);
    const auto *mi_657 = buffer.data(mi + 657);
    const auto *mi_658 = buffer.data(mi + 658);
    const auto *mi_659 = buffer.data(mi + 659);
    const auto *mi_660 = buffer.data(mi + 660);
    const auto *mi_661 = buffer.data(mi + 661);
    const auto *mi_662 = buffer.data(mi + 662);
    const auto *mi_663 = buffer.data(mi + 663);
    const auto *mi_664 = buffer.data(mi + 664);
    const auto *mi_665 = buffer.data(mi + 665);
    const auto *mi_666 = buffer.data(mi + 666);
    const auto *mi_667 = buffer.data(mi + 667);
    const auto *mi_668 = buffer.data(mi + 668);
    const auto *mi_669 = buffer.data(mi + 669);
    const auto *mi_670 = buffer.data(mi + 670);
    const auto *mi_671 = buffer.data(mi + 671);
    const auto *mi_672 = buffer.data(mi + 672);
    const auto *mi_673 = buffer.data(mi + 673);
    const auto *mi_674 = buffer.data(mi + 674);
    const auto *mi_675 = buffer.data(mi + 675);
    const auto *mi_676 = buffer.data(mi + 676);
    const auto *mi_677 = buffer.data(mi + 677);
    const auto *mi_678 = buffer.data(mi + 678);
    const auto *mi_679 = buffer.data(mi + 679);
    const auto *mi_680 = buffer.data(mi + 680);
    const auto *mi_681 = buffer.data(mi + 681);
    const auto *mi_682 = buffer.data(mi + 682);
    const auto *mi_683 = buffer.data(mi + 683);
    const auto *mi_684 = buffer.data(mi + 684);
    const auto *mi_685 = buffer.data(mi + 685);
    const auto *mi_686 = buffer.data(mi + 686);
    const auto *mi_687 = buffer.data(mi + 687);
    const auto *mi_688 = buffer.data(mi + 688);
    const auto *mi_689 = buffer.data(mi + 689);
    const auto *mi_690 = buffer.data(mi + 690);
    const auto *mi_691 = buffer.data(mi + 691);
    const auto *mi_692 = buffer.data(mi + 692);
    const auto *mi_693 = buffer.data(mi + 693);
    const auto *mi_694 = buffer.data(mi + 694);
    const auto *mi_695 = buffer.data(mi + 695);
    const auto *mi_696 = buffer.data(mi + 696);
    const auto *mi_697 = buffer.data(mi + 697);
    const auto *mi_698 = buffer.data(mi + 698);
    const auto *mi_699 = buffer.data(mi + 699);
    const auto *mi_700 = buffer.data(mi + 700);
    const auto *mi_701 = buffer.data(mi + 701);
    const auto *mi_702 = buffer.data(mi + 702);
    const auto *mi_703 = buffer.data(mi + 703);
    const auto *mi_704 = buffer.data(mi + 704);
    const auto *mi_705 = buffer.data(mi + 705);
    const auto *mi_706 = buffer.data(mi + 706);
    const auto *mi_707 = buffer.data(mi + 707);
    const auto *mi_708 = buffer.data(mi + 708);
    const auto *mi_709 = buffer.data(mi + 709);
    const auto *mi_710 = buffer.data(mi + 710);
    const auto *mi_711 = buffer.data(mi + 711);
    const auto *mi_712 = buffer.data(mi + 712);
    const auto *mi_713 = buffer.data(mi + 713);
    const auto *mi_714 = buffer.data(mi + 714);
    const auto *mi_715 = buffer.data(mi + 715);
    const auto *mi_716 = buffer.data(mi + 716);
    const auto *mi_717 = buffer.data(mi + 717);
    const auto *mi_718 = buffer.data(mi + 718);
    const auto *mi_719 = buffer.data(mi + 719);
    const auto *mi_720 = buffer.data(mi + 720);
    const auto *mi_721 = buffer.data(mi + 721);
    const auto *mi_722 = buffer.data(mi + 722);
    const auto *mi_723 = buffer.data(mi + 723);
    const auto *mi_724 = buffer.data(mi + 724);
    const auto *mi_725 = buffer.data(mi + 725);
    const auto *mi_726 = buffer.data(mi + 726);
    const auto *mi_727 = buffer.data(mi + 727);
    const auto *mi_728 = buffer.data(mi + 728);
    const auto *mi_729 = buffer.data(mi + 729);
    const auto *mi_730 = buffer.data(mi + 730);
    const auto *mi_731 = buffer.data(mi + 731);
    const auto *mi_732 = buffer.data(mi + 732);
    const auto *mi_733 = buffer.data(mi + 733);
    const auto *mi_734 = buffer.data(mi + 734);
    const auto *mi_735 = buffer.data(mi + 735);
    const auto *mi_736 = buffer.data(mi + 736);
    const auto *mi_737 = buffer.data(mi + 737);
    const auto *mi_738 = buffer.data(mi + 738);
    const auto *mi_739 = buffer.data(mi + 739);
    const auto *mi_740 = buffer.data(mi + 740);
    const auto *mi_741 = buffer.data(mi + 741);
    const auto *mi_742 = buffer.data(mi + 742);
    const auto *mi_743 = buffer.data(mi + 743);
    const auto *mi_744 = buffer.data(mi + 744);
    const auto *mi_745 = buffer.data(mi + 745);
    const auto *mi_746 = buffer.data(mi + 746);
    const auto *mi_747 = buffer.data(mi + 747);
    const auto *mi_748 = buffer.data(mi + 748);
    const auto *mi_749 = buffer.data(mi + 749);

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, ki_600, ki_601, ki_602, ki_603, \
                         ki_604, mi_600, mi_601, mi_602, mi_603, \
                         mi_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -2.0 * ki_600[k]
                   + f_0 * mi_600[k];

        t_601[k] = -2.0 * ki_601[k]
                   + f_0 * mi_601[k];

        t_602[k] = -2.0 * ki_602[k]
                   + f_0 * mi_602[k];

        t_603[k] = -2.0 * ki_603[k]
                   + f_0 * mi_603[k];

        t_604[k] = -2.0 * ki_604[k]
                   + f_0 * mi_604[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, ki_605, ki_606, ki_607, ki_608, \
                         ki_609, mi_605, mi_606, mi_607, mi_608, \
                         mi_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -2.0 * ki_605[k]
                   + f_0 * mi_605[k];

        t_606[k] = -2.0 * ki_606[k]
                   + f_0 * mi_606[k];

        t_607[k] = -2.0 * ki_607[k]
                   + f_0 * mi_607[k];

        t_608[k] = -2.0 * ki_608[k]
                   + f_0 * mi_608[k];

        t_609[k] = -2.0 * ki_609[k]
                   + f_0 * mi_609[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, ki_610, ki_611, ki_612, ki_613, \
                         ki_614, mi_610, mi_611, mi_612, mi_613, \
                         mi_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -2.0 * ki_610[k]
                   + f_0 * mi_610[k];

        t_611[k] = -2.0 * ki_611[k]
                   + f_0 * mi_611[k];

        t_612[k] = -2.0 * ki_612[k]
                   + f_0 * mi_612[k];

        t_613[k] = -2.0 * ki_613[k]
                   + f_0 * mi_613[k];

        t_614[k] = -2.0 * ki_614[k]
                   + f_0 * mi_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, ki_615, ki_616, ki_617, ki_618, \
                         ki_619, mi_615, mi_616, mi_617, mi_618, \
                         mi_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -2.0 * ki_615[k]
                   + f_0 * mi_615[k];

        t_616[k] = -2.0 * ki_616[k]
                   + f_0 * mi_616[k];

        t_617[k] = -2.0 * ki_617[k]
                   + f_0 * mi_617[k];

        t_618[k] = -2.0 * ki_618[k]
                   + f_0 * mi_618[k];

        t_619[k] = -2.0 * ki_619[k]
                   + f_0 * mi_619[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, ki_620, ki_621, ki_622, ki_623, \
                         ki_624, mi_620, mi_621, mi_622, mi_623, \
                         mi_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -2.0 * ki_620[k]
                   + f_0 * mi_620[k];

        t_621[k] = -2.0 * ki_621[k]
                   + f_0 * mi_621[k];

        t_622[k] = -2.0 * ki_622[k]
                   + f_0 * mi_622[k];

        t_623[k] = -2.0 * ki_623[k]
                   + f_0 * mi_623[k];

        t_624[k] = -2.0 * ki_624[k]
                   + f_0 * mi_624[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, ki_625, ki_626, ki_627, ki_628, \
                         ki_629, mi_625, mi_626, mi_627, mi_628, \
                         mi_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -2.0 * ki_625[k]
                   + f_0 * mi_625[k];

        t_626[k] = -2.0 * ki_626[k]
                   + f_0 * mi_626[k];

        t_627[k] = -2.0 * ki_627[k]
                   + f_0 * mi_627[k];

        t_628[k] = -2.0 * ki_628[k]
                   + f_0 * mi_628[k];

        t_629[k] = -2.0 * ki_629[k]
                   + f_0 * mi_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, ki_630, ki_631, ki_632, ki_633, \
                         ki_634, mi_630, mi_631, mi_632, mi_633, \
                         mi_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -2.0 * ki_630[k]
                   + f_0 * mi_630[k];

        t_631[k] = -2.0 * ki_631[k]
                   + f_0 * mi_631[k];

        t_632[k] = -2.0 * ki_632[k]
                   + f_0 * mi_632[k];

        t_633[k] = -2.0 * ki_633[k]
                   + f_0 * mi_633[k];

        t_634[k] = -2.0 * ki_634[k]
                   + f_0 * mi_634[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, ki_635, ki_636, ki_637, ki_638, \
                         ki_639, mi_635, mi_636, mi_637, mi_638, \
                         mi_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -2.0 * ki_635[k]
                   + f_0 * mi_635[k];

        t_636[k] = -2.0 * ki_636[k]
                   + f_0 * mi_636[k];

        t_637[k] = -2.0 * ki_637[k]
                   + f_0 * mi_637[k];

        t_638[k] = -2.0 * ki_638[k]
                   + f_0 * mi_638[k];

        t_639[k] = -2.0 * ki_639[k]
                   + f_0 * mi_639[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, ki_640, ki_641, ki_642, ki_643, \
                         ki_644, mi_640, mi_641, mi_642, mi_643, \
                         mi_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -2.0 * ki_640[k]
                   + f_0 * mi_640[k];

        t_641[k] = -2.0 * ki_641[k]
                   + f_0 * mi_641[k];

        t_642[k] = -2.0 * ki_642[k]
                   + f_0 * mi_642[k];

        t_643[k] = -2.0 * ki_643[k]
                   + f_0 * mi_643[k];

        t_644[k] = -2.0 * ki_644[k]
                   + f_0 * mi_644[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, ki_645, ki_646, ki_647, ki_648, \
                         ki_649, mi_645, mi_646, mi_647, mi_648, \
                         mi_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -2.0 * ki_645[k]
                   + f_0 * mi_645[k];

        t_646[k] = -2.0 * ki_646[k]
                   + f_0 * mi_646[k];

        t_647[k] = -2.0 * ki_647[k]
                   + f_0 * mi_647[k];

        t_648[k] = -2.0 * ki_648[k]
                   + f_0 * mi_648[k];

        t_649[k] = -2.0 * ki_649[k]
                   + f_0 * mi_649[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, ki_650, ki_651, ki_652, ki_653, \
                         ki_654, mi_650, mi_651, mi_652, mi_653, \
                         mi_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -2.0 * ki_650[k]
                   + f_0 * mi_650[k];

        t_651[k] = -2.0 * ki_651[k]
                   + f_0 * mi_651[k];

        t_652[k] = -2.0 * ki_652[k]
                   + f_0 * mi_652[k];

        t_653[k] = -2.0 * ki_653[k]
                   + f_0 * mi_653[k];

        t_654[k] = -2.0 * ki_654[k]
                   + f_0 * mi_654[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, ki_655, ki_656, ki_657, ki_658, \
                         ki_659, mi_655, mi_656, mi_657, mi_658, \
                         mi_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -2.0 * ki_655[k]
                   + f_0 * mi_655[k];

        t_656[k] = -2.0 * ki_656[k]
                   + f_0 * mi_656[k];

        t_657[k] = -2.0 * ki_657[k]
                   + f_0 * mi_657[k];

        t_658[k] = -2.0 * ki_658[k]
                   + f_0 * mi_658[k];

        t_659[k] = -2.0 * ki_659[k]
                   + f_0 * mi_659[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, ki_660, ki_661, ki_662, ki_663, \
                         ki_664, mi_660, mi_661, mi_662, mi_663, \
                         mi_664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -2.0 * ki_660[k]
                   + f_0 * mi_660[k];

        t_661[k] = -2.0 * ki_661[k]
                   + f_0 * mi_661[k];

        t_662[k] = -2.0 * ki_662[k]
                   + f_0 * mi_662[k];

        t_663[k] = -2.0 * ki_663[k]
                   + f_0 * mi_663[k];

        t_664[k] = -2.0 * ki_664[k]
                   + f_0 * mi_664[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, ki_665, ki_666, ki_667, ki_668, \
                         ki_669, mi_665, mi_666, mi_667, mi_668, \
                         mi_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -2.0 * ki_665[k]
                   + f_0 * mi_665[k];

        t_666[k] = -2.0 * ki_666[k]
                   + f_0 * mi_666[k];

        t_667[k] = -2.0 * ki_667[k]
                   + f_0 * mi_667[k];

        t_668[k] = -2.0 * ki_668[k]
                   + f_0 * mi_668[k];

        t_669[k] = -2.0 * ki_669[k]
                   + f_0 * mi_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, ki_670, ki_671, ki_672, ki_673, \
                         ki_674, mi_670, mi_671, mi_672, mi_673, \
                         mi_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -2.0 * ki_670[k]
                   + f_0 * mi_670[k];

        t_671[k] = -2.0 * ki_671[k]
                   + f_0 * mi_671[k];

        t_672[k] = -2.0 * ki_672[k]
                   + f_0 * mi_672[k];

        t_673[k] = -2.0 * ki_673[k]
                   + f_0 * mi_673[k];

        t_674[k] = -2.0 * ki_674[k]
                   + f_0 * mi_674[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, ki_675, ki_676, ki_677, ki_678, \
                         ki_679, mi_675, mi_676, mi_677, mi_678, \
                         mi_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -2.0 * ki_675[k]
                   + f_0 * mi_675[k];

        t_676[k] = -2.0 * ki_676[k]
                   + f_0 * mi_676[k];

        t_677[k] = -2.0 * ki_677[k]
                   + f_0 * mi_677[k];

        t_678[k] = -2.0 * ki_678[k]
                   + f_0 * mi_678[k];

        t_679[k] = -2.0 * ki_679[k]
                   + f_0 * mi_679[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, ki_680, ki_681, ki_682, ki_683, \
                         ki_684, mi_680, mi_681, mi_682, mi_683, \
                         mi_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -2.0 * ki_680[k]
                   + f_0 * mi_680[k];

        t_681[k] = -2.0 * ki_681[k]
                   + f_0 * mi_681[k];

        t_682[k] = -2.0 * ki_682[k]
                   + f_0 * mi_682[k];

        t_683[k] = -2.0 * ki_683[k]
                   + f_0 * mi_683[k];

        t_684[k] = -2.0 * ki_684[k]
                   + f_0 * mi_684[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, ki_685, ki_686, ki_687, ki_688, \
                         ki_689, mi_685, mi_686, mi_687, mi_688, \
                         mi_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -2.0 * ki_685[k]
                   + f_0 * mi_685[k];

        t_686[k] = -2.0 * ki_686[k]
                   + f_0 * mi_686[k];

        t_687[k] = -2.0 * ki_687[k]
                   + f_0 * mi_687[k];

        t_688[k] = -2.0 * ki_688[k]
                   + f_0 * mi_688[k];

        t_689[k] = -2.0 * ki_689[k]
                   + f_0 * mi_689[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, ki_690, ki_691, ki_692, ki_693, \
                         ki_694, mi_690, mi_691, mi_692, mi_693, \
                         mi_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -2.0 * ki_690[k]
                   + f_0 * mi_690[k];

        t_691[k] = -2.0 * ki_691[k]
                   + f_0 * mi_691[k];

        t_692[k] = -2.0 * ki_692[k]
                   + f_0 * mi_692[k];

        t_693[k] = -2.0 * ki_693[k]
                   + f_0 * mi_693[k];

        t_694[k] = -2.0 * ki_694[k]
                   + f_0 * mi_694[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, ki_695, ki_696, ki_697, ki_698, \
                         ki_699, mi_695, mi_696, mi_697, mi_698, \
                         mi_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -2.0 * ki_695[k]
                   + f_0 * mi_695[k];

        t_696[k] = -2.0 * ki_696[k]
                   + f_0 * mi_696[k];

        t_697[k] = -2.0 * ki_697[k]
                   + f_0 * mi_697[k];

        t_698[k] = -2.0 * ki_698[k]
                   + f_0 * mi_698[k];

        t_699[k] = -2.0 * ki_699[k]
                   + f_0 * mi_699[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, ki_700, ki_701, ki_702, ki_703, \
                         ki_704, mi_700, mi_701, mi_702, mi_703, \
                         mi_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -2.0 * ki_700[k]
                   + f_0 * mi_700[k];

        t_701[k] = -2.0 * ki_701[k]
                   + f_0 * mi_701[k];

        t_702[k] = -2.0 * ki_702[k]
                   + f_0 * mi_702[k];

        t_703[k] = -2.0 * ki_703[k]
                   + f_0 * mi_703[k];

        t_704[k] = -2.0 * ki_704[k]
                   + f_0 * mi_704[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, ki_705, ki_706, ki_707, ki_708, \
                         ki_709, mi_705, mi_706, mi_707, mi_708, \
                         mi_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -2.0 * ki_705[k]
                   + f_0 * mi_705[k];

        t_706[k] = -2.0 * ki_706[k]
                   + f_0 * mi_706[k];

        t_707[k] = -2.0 * ki_707[k]
                   + f_0 * mi_707[k];

        t_708[k] = -2.0 * ki_708[k]
                   + f_0 * mi_708[k];

        t_709[k] = -2.0 * ki_709[k]
                   + f_0 * mi_709[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, ki_710, ki_711, ki_712, ki_713, \
                         ki_714, mi_710, mi_711, mi_712, mi_713, \
                         mi_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -2.0 * ki_710[k]
                   + f_0 * mi_710[k];

        t_711[k] = -2.0 * ki_711[k]
                   + f_0 * mi_711[k];

        t_712[k] = -2.0 * ki_712[k]
                   + f_0 * mi_712[k];

        t_713[k] = -2.0 * ki_713[k]
                   + f_0 * mi_713[k];

        t_714[k] = -2.0 * ki_714[k]
                   + f_0 * mi_714[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, ki_715, ki_716, ki_717, ki_718, \
                         ki_719, mi_715, mi_716, mi_717, mi_718, \
                         mi_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -2.0 * ki_715[k]
                   + f_0 * mi_715[k];

        t_716[k] = -2.0 * ki_716[k]
                   + f_0 * mi_716[k];

        t_717[k] = -2.0 * ki_717[k]
                   + f_0 * mi_717[k];

        t_718[k] = -2.0 * ki_718[k]
                   + f_0 * mi_718[k];

        t_719[k] = -2.0 * ki_719[k]
                   + f_0 * mi_719[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, ki_720, ki_721, ki_722, ki_723, \
                         ki_724, mi_720, mi_721, mi_722, mi_723, \
                         mi_724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -2.0 * ki_720[k]
                   + f_0 * mi_720[k];

        t_721[k] = -2.0 * ki_721[k]
                   + f_0 * mi_721[k];

        t_722[k] = -2.0 * ki_722[k]
                   + f_0 * mi_722[k];

        t_723[k] = -2.0 * ki_723[k]
                   + f_0 * mi_723[k];

        t_724[k] = -2.0 * ki_724[k]
                   + f_0 * mi_724[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, ki_725, ki_726, ki_727, ki_728, \
                         ki_729, mi_725, mi_726, mi_727, mi_728, \
                         mi_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -2.0 * ki_725[k]
                   + f_0 * mi_725[k];

        t_726[k] = -2.0 * ki_726[k]
                   + f_0 * mi_726[k];

        t_727[k] = -2.0 * ki_727[k]
                   + f_0 * mi_727[k];

        t_728[k] = -2.0 * ki_728[k]
                   + f_0 * mi_728[k];

        t_729[k] = -2.0 * ki_729[k]
                   + f_0 * mi_729[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, ki_730, ki_731, ki_732, ki_733, \
                         ki_734, mi_730, mi_731, mi_732, mi_733, \
                         mi_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -2.0 * ki_730[k]
                   + f_0 * mi_730[k];

        t_731[k] = -2.0 * ki_731[k]
                   + f_0 * mi_731[k];

        t_732[k] = -2.0 * ki_732[k]
                   + f_0 * mi_732[k];

        t_733[k] = -2.0 * ki_733[k]
                   + f_0 * mi_733[k];

        t_734[k] = -2.0 * ki_734[k]
                   + f_0 * mi_734[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, ki_735, ki_736, ki_737, ki_738, \
                         ki_739, mi_735, mi_736, mi_737, mi_738, \
                         mi_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -2.0 * ki_735[k]
                   + f_0 * mi_735[k];

        t_736[k] = -2.0 * ki_736[k]
                   + f_0 * mi_736[k];

        t_737[k] = -2.0 * ki_737[k]
                   + f_0 * mi_737[k];

        t_738[k] = -2.0 * ki_738[k]
                   + f_0 * mi_738[k];

        t_739[k] = -2.0 * ki_739[k]
                   + f_0 * mi_739[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, ki_740, ki_741, ki_742, ki_743, \
                         ki_744, mi_740, mi_741, mi_742, mi_743, \
                         mi_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -2.0 * ki_740[k]
                   + f_0 * mi_740[k];

        t_741[k] = -2.0 * ki_741[k]
                   + f_0 * mi_741[k];

        t_742[k] = -2.0 * ki_742[k]
                   + f_0 * mi_742[k];

        t_743[k] = -2.0 * ki_743[k]
                   + f_0 * mi_743[k];

        t_744[k] = -2.0 * ki_744[k]
                   + f_0 * mi_744[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, ki_745, ki_746, ki_747, ki_748, \
                         ki_749, mi_745, mi_746, mi_747, mi_748, \
                         mi_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -2.0 * ki_745[k]
                   + f_0 * mi_745[k];

        t_746[k] = -2.0 * ki_746[k]
                   + f_0 * mi_746[k];

        t_747[k] = -2.0 * ki_747[k]
                   + f_0 * mi_747[k];

        t_748[k] = -2.0 * ki_748[k]
                   + f_0 * mi_748[k];

        t_749[k] = -2.0 * ki_749[k]
                   + f_0 * mi_749[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
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

    const auto *ki_750 = buffer.data(ki + 750);
    const auto *ki_751 = buffer.data(ki + 751);
    const auto *ki_752 = buffer.data(ki + 752);
    const auto *ki_753 = buffer.data(ki + 753);
    const auto *ki_754 = buffer.data(ki + 754);
    const auto *ki_755 = buffer.data(ki + 755);
    const auto *ki_756 = buffer.data(ki + 756);
    const auto *ki_757 = buffer.data(ki + 757);
    const auto *ki_758 = buffer.data(ki + 758);
    const auto *ki_759 = buffer.data(ki + 759);
    const auto *ki_760 = buffer.data(ki + 760);
    const auto *ki_761 = buffer.data(ki + 761);
    const auto *ki_762 = buffer.data(ki + 762);
    const auto *ki_763 = buffer.data(ki + 763);
    const auto *ki_764 = buffer.data(ki + 764);
    const auto *ki_765 = buffer.data(ki + 765);
    const auto *ki_766 = buffer.data(ki + 766);
    const auto *ki_767 = buffer.data(ki + 767);
    const auto *ki_768 = buffer.data(ki + 768);
    const auto *ki_769 = buffer.data(ki + 769);
    const auto *ki_770 = buffer.data(ki + 770);
    const auto *ki_771 = buffer.data(ki + 771);
    const auto *ki_772 = buffer.data(ki + 772);
    const auto *ki_773 = buffer.data(ki + 773);
    const auto *ki_774 = buffer.data(ki + 774);
    const auto *ki_775 = buffer.data(ki + 775);
    const auto *ki_776 = buffer.data(ki + 776);
    const auto *ki_777 = buffer.data(ki + 777);
    const auto *ki_778 = buffer.data(ki + 778);
    const auto *ki_779 = buffer.data(ki + 779);
    const auto *ki_780 = buffer.data(ki + 780);
    const auto *ki_781 = buffer.data(ki + 781);
    const auto *ki_782 = buffer.data(ki + 782);
    const auto *ki_783 = buffer.data(ki + 783);
    const auto *ki_784 = buffer.data(ki + 784);
    const auto *ki_785 = buffer.data(ki + 785);
    const auto *ki_786 = buffer.data(ki + 786);
    const auto *ki_787 = buffer.data(ki + 787);
    const auto *ki_788 = buffer.data(ki + 788);
    const auto *ki_789 = buffer.data(ki + 789);
    const auto *ki_790 = buffer.data(ki + 790);
    const auto *ki_791 = buffer.data(ki + 791);
    const auto *ki_792 = buffer.data(ki + 792);
    const auto *ki_793 = buffer.data(ki + 793);
    const auto *ki_794 = buffer.data(ki + 794);
    const auto *ki_795 = buffer.data(ki + 795);
    const auto *ki_796 = buffer.data(ki + 796);
    const auto *ki_797 = buffer.data(ki + 797);
    const auto *ki_798 = buffer.data(ki + 798);
    const auto *ki_799 = buffer.data(ki + 799);
    const auto *ki_800 = buffer.data(ki + 800);
    const auto *ki_801 = buffer.data(ki + 801);
    const auto *ki_802 = buffer.data(ki + 802);
    const auto *ki_803 = buffer.data(ki + 803);
    const auto *ki_804 = buffer.data(ki + 804);
    const auto *ki_805 = buffer.data(ki + 805);
    const auto *ki_806 = buffer.data(ki + 806);
    const auto *ki_807 = buffer.data(ki + 807);
    const auto *ki_808 = buffer.data(ki + 808);
    const auto *ki_809 = buffer.data(ki + 809);
    const auto *ki_810 = buffer.data(ki + 810);
    const auto *ki_811 = buffer.data(ki + 811);
    const auto *ki_812 = buffer.data(ki + 812);
    const auto *ki_813 = buffer.data(ki + 813);
    const auto *ki_814 = buffer.data(ki + 814);
    const auto *ki_815 = buffer.data(ki + 815);
    const auto *ki_816 = buffer.data(ki + 816);
    const auto *ki_817 = buffer.data(ki + 817);
    const auto *ki_818 = buffer.data(ki + 818);
    const auto *ki_819 = buffer.data(ki + 819);
    const auto *ki_820 = buffer.data(ki + 820);
    const auto *ki_821 = buffer.data(ki + 821);
    const auto *ki_822 = buffer.data(ki + 822);
    const auto *ki_823 = buffer.data(ki + 823);
    const auto *ki_824 = buffer.data(ki + 824);
    const auto *ki_825 = buffer.data(ki + 825);
    const auto *ki_826 = buffer.data(ki + 826);
    const auto *ki_827 = buffer.data(ki + 827);
    const auto *ki_828 = buffer.data(ki + 828);
    const auto *ki_829 = buffer.data(ki + 829);
    const auto *ki_830 = buffer.data(ki + 830);
    const auto *ki_831 = buffer.data(ki + 831);
    const auto *ki_832 = buffer.data(ki + 832);
    const auto *ki_833 = buffer.data(ki + 833);
    const auto *ki_834 = buffer.data(ki + 834);
    const auto *ki_835 = buffer.data(ki + 835);
    const auto *ki_836 = buffer.data(ki + 836);
    const auto *ki_837 = buffer.data(ki + 837);
    const auto *ki_838 = buffer.data(ki + 838);
    const auto *ki_839 = buffer.data(ki + 839);
    const auto *ki_840 = buffer.data(ki + 840);
    const auto *ki_841 = buffer.data(ki + 841);
    const auto *ki_842 = buffer.data(ki + 842);
    const auto *ki_843 = buffer.data(ki + 843);
    const auto *ki_844 = buffer.data(ki + 844);
    const auto *ki_845 = buffer.data(ki + 845);
    const auto *ki_846 = buffer.data(ki + 846);
    const auto *ki_847 = buffer.data(ki + 847);
    const auto *ki_848 = buffer.data(ki + 848);
    const auto *ki_849 = buffer.data(ki + 849);
    const auto *ki_850 = buffer.data(ki + 850);
    const auto *ki_851 = buffer.data(ki + 851);
    const auto *ki_852 = buffer.data(ki + 852);
    const auto *ki_853 = buffer.data(ki + 853);
    const auto *ki_854 = buffer.data(ki + 854);
    const auto *ki_855 = buffer.data(ki + 855);
    const auto *ki_856 = buffer.data(ki + 856);
    const auto *ki_857 = buffer.data(ki + 857);
    const auto *ki_858 = buffer.data(ki + 858);
    const auto *ki_859 = buffer.data(ki + 859);
    const auto *ki_860 = buffer.data(ki + 860);
    const auto *ki_861 = buffer.data(ki + 861);
    const auto *ki_862 = buffer.data(ki + 862);
    const auto *ki_863 = buffer.data(ki + 863);
    const auto *ki_864 = buffer.data(ki + 864);
    const auto *ki_865 = buffer.data(ki + 865);
    const auto *ki_866 = buffer.data(ki + 866);
    const auto *ki_867 = buffer.data(ki + 867);
    const auto *ki_868 = buffer.data(ki + 868);
    const auto *ki_869 = buffer.data(ki + 869);
    const auto *ki_870 = buffer.data(ki + 870);
    const auto *ki_871 = buffer.data(ki + 871);
    const auto *ki_872 = buffer.data(ki + 872);
    const auto *ki_873 = buffer.data(ki + 873);
    const auto *ki_874 = buffer.data(ki + 874);
    const auto *ki_875 = buffer.data(ki + 875);
    const auto *ki_876 = buffer.data(ki + 876);
    const auto *ki_877 = buffer.data(ki + 877);
    const auto *ki_878 = buffer.data(ki + 878);
    const auto *ki_879 = buffer.data(ki + 879);
    const auto *ki_880 = buffer.data(ki + 880);
    const auto *ki_881 = buffer.data(ki + 881);
    const auto *ki_882 = buffer.data(ki + 882);
    const auto *ki_883 = buffer.data(ki + 883);
    const auto *ki_884 = buffer.data(ki + 884);
    const auto *ki_885 = buffer.data(ki + 885);
    const auto *ki_886 = buffer.data(ki + 886);
    const auto *ki_887 = buffer.data(ki + 887);
    const auto *ki_888 = buffer.data(ki + 888);
    const auto *ki_889 = buffer.data(ki + 889);
    const auto *ki_890 = buffer.data(ki + 890);
    const auto *ki_891 = buffer.data(ki + 891);
    const auto *ki_892 = buffer.data(ki + 892);
    const auto *ki_893 = buffer.data(ki + 893);
    const auto *ki_894 = buffer.data(ki + 894);
    const auto *ki_895 = buffer.data(ki + 895);
    const auto *ki_896 = buffer.data(ki + 896);
    const auto *ki_897 = buffer.data(ki + 897);
    const auto *ki_898 = buffer.data(ki + 898);
    const auto *ki_899 = buffer.data(ki + 899);

    const auto *mi_750 = buffer.data(mi + 750);
    const auto *mi_751 = buffer.data(mi + 751);
    const auto *mi_752 = buffer.data(mi + 752);
    const auto *mi_753 = buffer.data(mi + 753);
    const auto *mi_754 = buffer.data(mi + 754);
    const auto *mi_755 = buffer.data(mi + 755);
    const auto *mi_756 = buffer.data(mi + 756);
    const auto *mi_757 = buffer.data(mi + 757);
    const auto *mi_758 = buffer.data(mi + 758);
    const auto *mi_759 = buffer.data(mi + 759);
    const auto *mi_760 = buffer.data(mi + 760);
    const auto *mi_761 = buffer.data(mi + 761);
    const auto *mi_762 = buffer.data(mi + 762);
    const auto *mi_763 = buffer.data(mi + 763);
    const auto *mi_764 = buffer.data(mi + 764);
    const auto *mi_765 = buffer.data(mi + 765);
    const auto *mi_766 = buffer.data(mi + 766);
    const auto *mi_767 = buffer.data(mi + 767);
    const auto *mi_768 = buffer.data(mi + 768);
    const auto *mi_769 = buffer.data(mi + 769);
    const auto *mi_770 = buffer.data(mi + 770);
    const auto *mi_771 = buffer.data(mi + 771);
    const auto *mi_772 = buffer.data(mi + 772);
    const auto *mi_773 = buffer.data(mi + 773);
    const auto *mi_774 = buffer.data(mi + 774);
    const auto *mi_775 = buffer.data(mi + 775);
    const auto *mi_776 = buffer.data(mi + 776);
    const auto *mi_777 = buffer.data(mi + 777);
    const auto *mi_778 = buffer.data(mi + 778);
    const auto *mi_779 = buffer.data(mi + 779);
    const auto *mi_780 = buffer.data(mi + 780);
    const auto *mi_781 = buffer.data(mi + 781);
    const auto *mi_782 = buffer.data(mi + 782);
    const auto *mi_783 = buffer.data(mi + 783);
    const auto *mi_784 = buffer.data(mi + 784);
    const auto *mi_785 = buffer.data(mi + 785);
    const auto *mi_786 = buffer.data(mi + 786);
    const auto *mi_787 = buffer.data(mi + 787);
    const auto *mi_788 = buffer.data(mi + 788);
    const auto *mi_789 = buffer.data(mi + 789);
    const auto *mi_790 = buffer.data(mi + 790);
    const auto *mi_791 = buffer.data(mi + 791);
    const auto *mi_792 = buffer.data(mi + 792);
    const auto *mi_793 = buffer.data(mi + 793);
    const auto *mi_794 = buffer.data(mi + 794);
    const auto *mi_795 = buffer.data(mi + 795);
    const auto *mi_796 = buffer.data(mi + 796);
    const auto *mi_797 = buffer.data(mi + 797);
    const auto *mi_798 = buffer.data(mi + 798);
    const auto *mi_799 = buffer.data(mi + 799);
    const auto *mi_800 = buffer.data(mi + 800);
    const auto *mi_801 = buffer.data(mi + 801);
    const auto *mi_802 = buffer.data(mi + 802);
    const auto *mi_803 = buffer.data(mi + 803);
    const auto *mi_804 = buffer.data(mi + 804);
    const auto *mi_805 = buffer.data(mi + 805);
    const auto *mi_806 = buffer.data(mi + 806);
    const auto *mi_807 = buffer.data(mi + 807);
    const auto *mi_808 = buffer.data(mi + 808);
    const auto *mi_809 = buffer.data(mi + 809);
    const auto *mi_810 = buffer.data(mi + 810);
    const auto *mi_811 = buffer.data(mi + 811);
    const auto *mi_812 = buffer.data(mi + 812);
    const auto *mi_813 = buffer.data(mi + 813);
    const auto *mi_814 = buffer.data(mi + 814);
    const auto *mi_815 = buffer.data(mi + 815);
    const auto *mi_816 = buffer.data(mi + 816);
    const auto *mi_817 = buffer.data(mi + 817);
    const auto *mi_818 = buffer.data(mi + 818);
    const auto *mi_819 = buffer.data(mi + 819);
    const auto *mi_820 = buffer.data(mi + 820);
    const auto *mi_821 = buffer.data(mi + 821);
    const auto *mi_822 = buffer.data(mi + 822);
    const auto *mi_823 = buffer.data(mi + 823);
    const auto *mi_824 = buffer.data(mi + 824);
    const auto *mi_825 = buffer.data(mi + 825);
    const auto *mi_826 = buffer.data(mi + 826);
    const auto *mi_827 = buffer.data(mi + 827);
    const auto *mi_828 = buffer.data(mi + 828);
    const auto *mi_829 = buffer.data(mi + 829);
    const auto *mi_830 = buffer.data(mi + 830);
    const auto *mi_831 = buffer.data(mi + 831);
    const auto *mi_832 = buffer.data(mi + 832);
    const auto *mi_833 = buffer.data(mi + 833);
    const auto *mi_834 = buffer.data(mi + 834);
    const auto *mi_835 = buffer.data(mi + 835);
    const auto *mi_836 = buffer.data(mi + 836);
    const auto *mi_837 = buffer.data(mi + 837);
    const auto *mi_838 = buffer.data(mi + 838);
    const auto *mi_839 = buffer.data(mi + 839);
    const auto *mi_840 = buffer.data(mi + 840);
    const auto *mi_841 = buffer.data(mi + 841);
    const auto *mi_842 = buffer.data(mi + 842);
    const auto *mi_843 = buffer.data(mi + 843);
    const auto *mi_844 = buffer.data(mi + 844);
    const auto *mi_845 = buffer.data(mi + 845);
    const auto *mi_846 = buffer.data(mi + 846);
    const auto *mi_847 = buffer.data(mi + 847);
    const auto *mi_848 = buffer.data(mi + 848);
    const auto *mi_849 = buffer.data(mi + 849);
    const auto *mi_850 = buffer.data(mi + 850);
    const auto *mi_851 = buffer.data(mi + 851);
    const auto *mi_852 = buffer.data(mi + 852);
    const auto *mi_853 = buffer.data(mi + 853);
    const auto *mi_854 = buffer.data(mi + 854);
    const auto *mi_855 = buffer.data(mi + 855);
    const auto *mi_856 = buffer.data(mi + 856);
    const auto *mi_857 = buffer.data(mi + 857);
    const auto *mi_858 = buffer.data(mi + 858);
    const auto *mi_859 = buffer.data(mi + 859);
    const auto *mi_860 = buffer.data(mi + 860);
    const auto *mi_861 = buffer.data(mi + 861);
    const auto *mi_862 = buffer.data(mi + 862);
    const auto *mi_863 = buffer.data(mi + 863);
    const auto *mi_864 = buffer.data(mi + 864);
    const auto *mi_865 = buffer.data(mi + 865);
    const auto *mi_866 = buffer.data(mi + 866);
    const auto *mi_867 = buffer.data(mi + 867);
    const auto *mi_868 = buffer.data(mi + 868);
    const auto *mi_869 = buffer.data(mi + 869);
    const auto *mi_870 = buffer.data(mi + 870);
    const auto *mi_871 = buffer.data(mi + 871);
    const auto *mi_872 = buffer.data(mi + 872);
    const auto *mi_873 = buffer.data(mi + 873);
    const auto *mi_874 = buffer.data(mi + 874);
    const auto *mi_875 = buffer.data(mi + 875);
    const auto *mi_876 = buffer.data(mi + 876);
    const auto *mi_877 = buffer.data(mi + 877);
    const auto *mi_878 = buffer.data(mi + 878);
    const auto *mi_879 = buffer.data(mi + 879);
    const auto *mi_880 = buffer.data(mi + 880);
    const auto *mi_881 = buffer.data(mi + 881);
    const auto *mi_882 = buffer.data(mi + 882);
    const auto *mi_883 = buffer.data(mi + 883);
    const auto *mi_884 = buffer.data(mi + 884);
    const auto *mi_885 = buffer.data(mi + 885);
    const auto *mi_886 = buffer.data(mi + 886);
    const auto *mi_887 = buffer.data(mi + 887);
    const auto *mi_888 = buffer.data(mi + 888);
    const auto *mi_889 = buffer.data(mi + 889);
    const auto *mi_890 = buffer.data(mi + 890);
    const auto *mi_891 = buffer.data(mi + 891);
    const auto *mi_892 = buffer.data(mi + 892);
    const auto *mi_893 = buffer.data(mi + 893);
    const auto *mi_894 = buffer.data(mi + 894);
    const auto *mi_895 = buffer.data(mi + 895);
    const auto *mi_896 = buffer.data(mi + 896);
    const auto *mi_897 = buffer.data(mi + 897);
    const auto *mi_898 = buffer.data(mi + 898);
    const auto *mi_899 = buffer.data(mi + 899);

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, ki_750, ki_751, ki_752, ki_753, \
                         ki_754, mi_750, mi_751, mi_752, mi_753, \
                         mi_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -2.0 * ki_750[k]
                   + f_0 * mi_750[k];

        t_751[k] = -2.0 * ki_751[k]
                   + f_0 * mi_751[k];

        t_752[k] = -2.0 * ki_752[k]
                   + f_0 * mi_752[k];

        t_753[k] = -2.0 * ki_753[k]
                   + f_0 * mi_753[k];

        t_754[k] = -2.0 * ki_754[k]
                   + f_0 * mi_754[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, ki_755, ki_756, ki_757, ki_758, \
                         ki_759, mi_755, mi_756, mi_757, mi_758, \
                         mi_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -2.0 * ki_755[k]
                   + f_0 * mi_755[k];

        t_756[k] = -2.0 * ki_756[k]
                   + f_0 * mi_756[k];

        t_757[k] = -2.0 * ki_757[k]
                   + f_0 * mi_757[k];

        t_758[k] = -2.0 * ki_758[k]
                   + f_0 * mi_758[k];

        t_759[k] = -2.0 * ki_759[k]
                   + f_0 * mi_759[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, t_764, ki_760, ki_761, ki_762, ki_763, \
                         ki_764, mi_760, mi_761, mi_762, mi_763, \
                         mi_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = -2.0 * ki_760[k]
                   + f_0 * mi_760[k];

        t_761[k] = -2.0 * ki_761[k]
                   + f_0 * mi_761[k];

        t_762[k] = -2.0 * ki_762[k]
                   + f_0 * mi_762[k];

        t_763[k] = -2.0 * ki_763[k]
                   + f_0 * mi_763[k];

        t_764[k] = -2.0 * ki_764[k]
                   + f_0 * mi_764[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, t_769, ki_765, ki_766, ki_767, ki_768, \
                         ki_769, mi_765, mi_766, mi_767, mi_768, \
                         mi_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = -2.0 * ki_765[k]
                   + f_0 * mi_765[k];

        t_766[k] = -2.0 * ki_766[k]
                   + f_0 * mi_766[k];

        t_767[k] = -2.0 * ki_767[k]
                   + f_0 * mi_767[k];

        t_768[k] = -2.0 * ki_768[k]
                   + f_0 * mi_768[k];

        t_769[k] = -2.0 * ki_769[k]
                   + f_0 * mi_769[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, ki_770, ki_771, ki_772, ki_773, \
                         ki_774, mi_770, mi_771, mi_772, mi_773, \
                         mi_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = -2.0 * ki_770[k]
                   + f_0 * mi_770[k];

        t_771[k] = -2.0 * ki_771[k]
                   + f_0 * mi_771[k];

        t_772[k] = -2.0 * ki_772[k]
                   + f_0 * mi_772[k];

        t_773[k] = -2.0 * ki_773[k]
                   + f_0 * mi_773[k];

        t_774[k] = -2.0 * ki_774[k]
                   + f_0 * mi_774[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, ki_775, ki_776, ki_777, ki_778, \
                         ki_779, mi_775, mi_776, mi_777, mi_778, \
                         mi_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = -2.0 * ki_775[k]
                   + f_0 * mi_775[k];

        t_776[k] = -2.0 * ki_776[k]
                   + f_0 * mi_776[k];

        t_777[k] = -2.0 * ki_777[k]
                   + f_0 * mi_777[k];

        t_778[k] = -2.0 * ki_778[k]
                   + f_0 * mi_778[k];

        t_779[k] = -2.0 * ki_779[k]
                   + f_0 * mi_779[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, t_784, ki_780, ki_781, ki_782, ki_783, \
                         ki_784, mi_780, mi_781, mi_782, mi_783, \
                         mi_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = -2.0 * ki_780[k]
                   + f_0 * mi_780[k];

        t_781[k] = -2.0 * ki_781[k]
                   + f_0 * mi_781[k];

        t_782[k] = -2.0 * ki_782[k]
                   + f_0 * mi_782[k];

        t_783[k] = -2.0 * ki_783[k]
                   + f_0 * mi_783[k];

        t_784[k] = -ki_784[k]
                   + f_0 * mi_784[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, t_789, ki_785, ki_786, ki_787, ki_788, \
                         ki_789, mi_785, mi_786, mi_787, mi_788, \
                         mi_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = -ki_785[k]
                   + f_0 * mi_785[k];

        t_786[k] = -ki_786[k]
                   + f_0 * mi_786[k];

        t_787[k] = -ki_787[k]
                   + f_0 * mi_787[k];

        t_788[k] = -ki_788[k]
                   + f_0 * mi_788[k];

        t_789[k] = -ki_789[k]
                   + f_0 * mi_789[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, t_794, ki_790, ki_791, ki_792, ki_793, \
                         ki_794, mi_790, mi_791, mi_792, mi_793, \
                         mi_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = -ki_790[k]
                   + f_0 * mi_790[k];

        t_791[k] = -ki_791[k]
                   + f_0 * mi_791[k];

        t_792[k] = -ki_792[k]
                   + f_0 * mi_792[k];

        t_793[k] = -ki_793[k]
                   + f_0 * mi_793[k];

        t_794[k] = -ki_794[k]
                   + f_0 * mi_794[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, ki_795, ki_796, ki_797, ki_798, \
                         ki_799, mi_795, mi_796, mi_797, mi_798, \
                         mi_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = -ki_795[k]
                   + f_0 * mi_795[k];

        t_796[k] = -ki_796[k]
                   + f_0 * mi_796[k];

        t_797[k] = -ki_797[k]
                   + f_0 * mi_797[k];

        t_798[k] = -ki_798[k]
                   + f_0 * mi_798[k];

        t_799[k] = -ki_799[k]
                   + f_0 * mi_799[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, ki_800, ki_801, ki_802, ki_803, \
                         ki_804, mi_800, mi_801, mi_802, mi_803, \
                         mi_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = -ki_800[k]
                   + f_0 * mi_800[k];

        t_801[k] = -ki_801[k]
                   + f_0 * mi_801[k];

        t_802[k] = -ki_802[k]
                   + f_0 * mi_802[k];

        t_803[k] = -ki_803[k]
                   + f_0 * mi_803[k];

        t_804[k] = -ki_804[k]
                   + f_0 * mi_804[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, ki_805, ki_806, ki_807, ki_808, \
                         ki_809, mi_805, mi_806, mi_807, mi_808, \
                         mi_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = -ki_805[k]
                   + f_0 * mi_805[k];

        t_806[k] = -ki_806[k]
                   + f_0 * mi_806[k];

        t_807[k] = -ki_807[k]
                   + f_0 * mi_807[k];

        t_808[k] = -ki_808[k]
                   + f_0 * mi_808[k];

        t_809[k] = -ki_809[k]
                   + f_0 * mi_809[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, ki_810, ki_811, ki_812, ki_813, \
                         ki_814, mi_810, mi_811, mi_812, mi_813, \
                         mi_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = -ki_810[k]
                   + f_0 * mi_810[k];

        t_811[k] = -ki_811[k]
                   + f_0 * mi_811[k];

        t_812[k] = -ki_812[k]
                   + f_0 * mi_812[k];

        t_813[k] = -ki_813[k]
                   + f_0 * mi_813[k];

        t_814[k] = -ki_814[k]
                   + f_0 * mi_814[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, ki_815, ki_816, ki_817, ki_818, \
                         ki_819, mi_815, mi_816, mi_817, mi_818, \
                         mi_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = -ki_815[k]
                   + f_0 * mi_815[k];

        t_816[k] = -ki_816[k]
                   + f_0 * mi_816[k];

        t_817[k] = -ki_817[k]
                   + f_0 * mi_817[k];

        t_818[k] = -ki_818[k]
                   + f_0 * mi_818[k];

        t_819[k] = -ki_819[k]
                   + f_0 * mi_819[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, t_824, ki_820, ki_821, ki_822, ki_823, \
                         ki_824, mi_820, mi_821, mi_822, mi_823, \
                         mi_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = -ki_820[k]
                   + f_0 * mi_820[k];

        t_821[k] = -ki_821[k]
                   + f_0 * mi_821[k];

        t_822[k] = -ki_822[k]
                   + f_0 * mi_822[k];

        t_823[k] = -ki_823[k]
                   + f_0 * mi_823[k];

        t_824[k] = -ki_824[k]
                   + f_0 * mi_824[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, ki_825, ki_826, ki_827, ki_828, \
                         ki_829, mi_825, mi_826, mi_827, mi_828, \
                         mi_829 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = -ki_825[k]
                   + f_0 * mi_825[k];

        t_826[k] = -ki_826[k]
                   + f_0 * mi_826[k];

        t_827[k] = -ki_827[k]
                   + f_0 * mi_827[k];

        t_828[k] = -ki_828[k]
                   + f_0 * mi_828[k];

        t_829[k] = -ki_829[k]
                   + f_0 * mi_829[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, ki_830, ki_831, ki_832, ki_833, \
                         ki_834, mi_830, mi_831, mi_832, mi_833, \
                         mi_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = -ki_830[k]
                   + f_0 * mi_830[k];

        t_831[k] = -ki_831[k]
                   + f_0 * mi_831[k];

        t_832[k] = -ki_832[k]
                   + f_0 * mi_832[k];

        t_833[k] = -ki_833[k]
                   + f_0 * mi_833[k];

        t_834[k] = -ki_834[k]
                   + f_0 * mi_834[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, t_838, t_839, ki_835, ki_836, ki_837, ki_838, \
                         ki_839, mi_835, mi_836, mi_837, mi_838, \
                         mi_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = -ki_835[k]
                   + f_0 * mi_835[k];

        t_836[k] = -ki_836[k]
                   + f_0 * mi_836[k];

        t_837[k] = -ki_837[k]
                   + f_0 * mi_837[k];

        t_838[k] = -ki_838[k]
                   + f_0 * mi_838[k];

        t_839[k] = -ki_839[k]
                   + f_0 * mi_839[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, ki_840, ki_841, ki_842, ki_843, \
                         ki_844, mi_840, mi_841, mi_842, mi_843, \
                         mi_844 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = -ki_840[k]
                   + f_0 * mi_840[k];

        t_841[k] = -ki_841[k]
                   + f_0 * mi_841[k];

        t_842[k] = -ki_842[k]
                   + f_0 * mi_842[k];

        t_843[k] = -ki_843[k]
                   + f_0 * mi_843[k];

        t_844[k] = -ki_844[k]
                   + f_0 * mi_844[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, t_849, ki_845, ki_846, ki_847, ki_848, \
                         ki_849, mi_845, mi_846, mi_847, mi_848, \
                         mi_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = -ki_845[k]
                   + f_0 * mi_845[k];

        t_846[k] = -ki_846[k]
                   + f_0 * mi_846[k];

        t_847[k] = -ki_847[k]
                   + f_0 * mi_847[k];

        t_848[k] = -ki_848[k]
                   + f_0 * mi_848[k];

        t_849[k] = -ki_849[k]
                   + f_0 * mi_849[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, ki_850, ki_851, ki_852, ki_853, \
                         ki_854, mi_850, mi_851, mi_852, mi_853, \
                         mi_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = -ki_850[k]
                   + f_0 * mi_850[k];

        t_851[k] = -ki_851[k]
                   + f_0 * mi_851[k];

        t_852[k] = -ki_852[k]
                   + f_0 * mi_852[k];

        t_853[k] = -ki_853[k]
                   + f_0 * mi_853[k];

        t_854[k] = -ki_854[k]
                   + f_0 * mi_854[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, ki_855, ki_856, ki_857, ki_858, \
                         ki_859, mi_855, mi_856, mi_857, mi_858, \
                         mi_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = -ki_855[k]
                   + f_0 * mi_855[k];

        t_856[k] = -ki_856[k]
                   + f_0 * mi_856[k];

        t_857[k] = -ki_857[k]
                   + f_0 * mi_857[k];

        t_858[k] = -ki_858[k]
                   + f_0 * mi_858[k];

        t_859[k] = -ki_859[k]
                   + f_0 * mi_859[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, t_864, ki_860, ki_861, ki_862, ki_863, \
                         ki_864, mi_860, mi_861, mi_862, mi_863, \
                         mi_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = -ki_860[k]
                   + f_0 * mi_860[k];

        t_861[k] = -ki_861[k]
                   + f_0 * mi_861[k];

        t_862[k] = -ki_862[k]
                   + f_0 * mi_862[k];

        t_863[k] = -ki_863[k]
                   + f_0 * mi_863[k];

        t_864[k] = -ki_864[k]
                   + f_0 * mi_864[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, t_869, ki_865, ki_866, ki_867, ki_868, \
                         ki_869, mi_865, mi_866, mi_867, mi_868, \
                         mi_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = -ki_865[k]
                   + f_0 * mi_865[k];

        t_866[k] = -ki_866[k]
                   + f_0 * mi_866[k];

        t_867[k] = -ki_867[k]
                   + f_0 * mi_867[k];

        t_868[k] = -ki_868[k]
                   + f_0 * mi_868[k];

        t_869[k] = -ki_869[k]
                   + f_0 * mi_869[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, t_874, ki_870, ki_871, ki_872, ki_873, \
                         ki_874, mi_870, mi_871, mi_872, mi_873, \
                         mi_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = -ki_870[k]
                   + f_0 * mi_870[k];

        t_871[k] = -ki_871[k]
                   + f_0 * mi_871[k];

        t_872[k] = -ki_872[k]
                   + f_0 * mi_872[k];

        t_873[k] = -ki_873[k]
                   + f_0 * mi_873[k];

        t_874[k] = -ki_874[k]
                   + f_0 * mi_874[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, t_879, ki_875, ki_876, ki_877, ki_878, \
                         ki_879, mi_875, mi_876, mi_877, mi_878, \
                         mi_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = -ki_875[k]
                   + f_0 * mi_875[k];

        t_876[k] = -ki_876[k]
                   + f_0 * mi_876[k];

        t_877[k] = -ki_877[k]
                   + f_0 * mi_877[k];

        t_878[k] = -ki_878[k]
                   + f_0 * mi_878[k];

        t_879[k] = -ki_879[k]
                   + f_0 * mi_879[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, t_884, ki_880, ki_881, ki_882, ki_883, \
                         ki_884, mi_880, mi_881, mi_882, mi_883, \
                         mi_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = -ki_880[k]
                   + f_0 * mi_880[k];

        t_881[k] = -ki_881[k]
                   + f_0 * mi_881[k];

        t_882[k] = -ki_882[k]
                   + f_0 * mi_882[k];

        t_883[k] = -ki_883[k]
                   + f_0 * mi_883[k];

        t_884[k] = -ki_884[k]
                   + f_0 * mi_884[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, ki_885, ki_886, ki_887, ki_888, \
                         ki_889, mi_885, mi_886, mi_887, mi_888, \
                         mi_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = -ki_885[k]
                   + f_0 * mi_885[k];

        t_886[k] = -ki_886[k]
                   + f_0 * mi_886[k];

        t_887[k] = -ki_887[k]
                   + f_0 * mi_887[k];

        t_888[k] = -ki_888[k]
                   + f_0 * mi_888[k];

        t_889[k] = -ki_889[k]
                   + f_0 * mi_889[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, ki_890, ki_891, ki_892, ki_893, \
                         ki_894, mi_890, mi_891, mi_892, mi_893, \
                         mi_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = -ki_890[k]
                   + f_0 * mi_890[k];

        t_891[k] = -ki_891[k]
                   + f_0 * mi_891[k];

        t_892[k] = -ki_892[k]
                   + f_0 * mi_892[k];

        t_893[k] = -ki_893[k]
                   + f_0 * mi_893[k];

        t_894[k] = -ki_894[k]
                   + f_0 * mi_894[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, t_899, ki_895, ki_896, ki_897, ki_898, \
                         ki_899, mi_895, mi_896, mi_897, mi_898, \
                         mi_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = -ki_895[k]
                   + f_0 * mi_895[k];

        t_896[k] = -ki_896[k]
                   + f_0 * mi_896[k];

        t_897[k] = -ki_897[k]
                   + f_0 * mi_897[k];

        t_898[k] = -ki_898[k]
                   + f_0 * mi_898[k];

        t_899[k] = -ki_899[k]
                   + f_0 * mi_899[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_0_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);
    auto *t_948 = buffer.data(target + 948);
    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);
    auto *t_955 = buffer.data(target + 955);
    auto *t_956 = buffer.data(target + 956);
    auto *t_957 = buffer.data(target + 957);
    auto *t_958 = buffer.data(target + 958);
    auto *t_959 = buffer.data(target + 959);
    auto *t_960 = buffer.data(target + 960);
    auto *t_961 = buffer.data(target + 961);
    auto *t_962 = buffer.data(target + 962);
    auto *t_963 = buffer.data(target + 963);
    auto *t_964 = buffer.data(target + 964);
    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);
    auto *t_970 = buffer.data(target + 970);
    auto *t_971 = buffer.data(target + 971);
    auto *t_972 = buffer.data(target + 972);
    auto *t_973 = buffer.data(target + 973);
    auto *t_974 = buffer.data(target + 974);
    auto *t_975 = buffer.data(target + 975);
    auto *t_976 = buffer.data(target + 976);
    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);
    auto *t_983 = buffer.data(target + 983);
    auto *t_984 = buffer.data(target + 984);
    auto *t_985 = buffer.data(target + 985);
    auto *t_986 = buffer.data(target + 986);
    auto *t_987 = buffer.data(target + 987);
    auto *t_988 = buffer.data(target + 988);
    auto *t_989 = buffer.data(target + 989);
    auto *t_990 = buffer.data(target + 990);
    auto *t_991 = buffer.data(target + 991);
    auto *t_992 = buffer.data(target + 992);
    auto *t_993 = buffer.data(target + 993);
    auto *t_994 = buffer.data(target + 994);
    auto *t_995 = buffer.data(target + 995);
    auto *t_996 = buffer.data(target + 996);
    auto *t_997 = buffer.data(target + 997);
    auto *t_998 = buffer.data(target + 998);
    auto *t_999 = buffer.data(target + 999);
    auto *t_1000 = buffer.data(target + 1000);
    auto *t_1001 = buffer.data(target + 1001);
    auto *t_1002 = buffer.data(target + 1002);
    auto *t_1003 = buffer.data(target + 1003);
    auto *t_1004 = buffer.data(target + 1004);
    auto *t_1005 = buffer.data(target + 1005);
    auto *t_1006 = buffer.data(target + 1006);
    auto *t_1007 = buffer.data(target + 1007);
    auto *t_1008 = buffer.data(target + 1008);
    auto *t_1009 = buffer.data(target + 1009);
    auto *t_1010 = buffer.data(target + 1010);
    auto *t_1011 = buffer.data(target + 1011);
    auto *t_1012 = buffer.data(target + 1012);
    auto *t_1013 = buffer.data(target + 1013);
    auto *t_1014 = buffer.data(target + 1014);
    auto *t_1015 = buffer.data(target + 1015);
    auto *t_1016 = buffer.data(target + 1016);
    auto *t_1017 = buffer.data(target + 1017);
    auto *t_1018 = buffer.data(target + 1018);
    auto *t_1019 = buffer.data(target + 1019);
    auto *t_1020 = buffer.data(target + 1020);
    auto *t_1021 = buffer.data(target + 1021);
    auto *t_1022 = buffer.data(target + 1022);
    auto *t_1023 = buffer.data(target + 1023);
    auto *t_1024 = buffer.data(target + 1024);
    auto *t_1025 = buffer.data(target + 1025);
    auto *t_1026 = buffer.data(target + 1026);
    auto *t_1027 = buffer.data(target + 1027);
    auto *t_1028 = buffer.data(target + 1028);
    auto *t_1029 = buffer.data(target + 1029);
    auto *t_1030 = buffer.data(target + 1030);
    auto *t_1031 = buffer.data(target + 1031);
    auto *t_1032 = buffer.data(target + 1032);
    auto *t_1033 = buffer.data(target + 1033);
    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);
    auto *t_1050 = buffer.data(target + 1050);
    auto *t_1051 = buffer.data(target + 1051);
    auto *t_1052 = buffer.data(target + 1052);
    auto *t_1053 = buffer.data(target + 1053);
    auto *t_1054 = buffer.data(target + 1054);
    auto *t_1055 = buffer.data(target + 1055);
    auto *t_1056 = buffer.data(target + 1056);
    auto *t_1057 = buffer.data(target + 1057);
    auto *t_1058 = buffer.data(target + 1058);
    auto *t_1059 = buffer.data(target + 1059);
    auto *t_1060 = buffer.data(target + 1060);
    auto *t_1061 = buffer.data(target + 1061);
    auto *t_1062 = buffer.data(target + 1062);
    auto *t_1063 = buffer.data(target + 1063);
    auto *t_1064 = buffer.data(target + 1064);
    auto *t_1065 = buffer.data(target + 1065);
    auto *t_1066 = buffer.data(target + 1066);

    const auto *ki_900 = buffer.data(ki + 900);
    const auto *ki_901 = buffer.data(ki + 901);
    const auto *ki_902 = buffer.data(ki + 902);
    const auto *ki_903 = buffer.data(ki + 903);
    const auto *ki_904 = buffer.data(ki + 904);
    const auto *ki_905 = buffer.data(ki + 905);
    const auto *ki_906 = buffer.data(ki + 906);
    const auto *ki_907 = buffer.data(ki + 907);
    const auto *ki_908 = buffer.data(ki + 908);
    const auto *ki_909 = buffer.data(ki + 909);
    const auto *ki_910 = buffer.data(ki + 910);
    const auto *ki_911 = buffer.data(ki + 911);
    const auto *ki_912 = buffer.data(ki + 912);
    const auto *ki_913 = buffer.data(ki + 913);
    const auto *ki_914 = buffer.data(ki + 914);
    const auto *ki_915 = buffer.data(ki + 915);
    const auto *ki_916 = buffer.data(ki + 916);
    const auto *ki_917 = buffer.data(ki + 917);
    const auto *ki_918 = buffer.data(ki + 918);
    const auto *ki_919 = buffer.data(ki + 919);
    const auto *ki_920 = buffer.data(ki + 920);
    const auto *ki_921 = buffer.data(ki + 921);
    const auto *ki_922 = buffer.data(ki + 922);
    const auto *ki_923 = buffer.data(ki + 923);
    const auto *ki_924 = buffer.data(ki + 924);
    const auto *ki_925 = buffer.data(ki + 925);
    const auto *ki_926 = buffer.data(ki + 926);
    const auto *ki_927 = buffer.data(ki + 927);
    const auto *ki_928 = buffer.data(ki + 928);
    const auto *ki_929 = buffer.data(ki + 929);
    const auto *ki_930 = buffer.data(ki + 930);
    const auto *ki_931 = buffer.data(ki + 931);
    const auto *ki_932 = buffer.data(ki + 932);
    const auto *ki_933 = buffer.data(ki + 933);
    const auto *ki_934 = buffer.data(ki + 934);
    const auto *ki_935 = buffer.data(ki + 935);
    const auto *ki_936 = buffer.data(ki + 936);
    const auto *ki_937 = buffer.data(ki + 937);
    const auto *ki_938 = buffer.data(ki + 938);
    const auto *ki_939 = buffer.data(ki + 939);
    const auto *ki_940 = buffer.data(ki + 940);
    const auto *ki_941 = buffer.data(ki + 941);
    const auto *ki_942 = buffer.data(ki + 942);
    const auto *ki_943 = buffer.data(ki + 943);
    const auto *ki_944 = buffer.data(ki + 944);
    const auto *ki_945 = buffer.data(ki + 945);
    const auto *ki_946 = buffer.data(ki + 946);
    const auto *ki_947 = buffer.data(ki + 947);
    const auto *ki_948 = buffer.data(ki + 948);
    const auto *ki_949 = buffer.data(ki + 949);
    const auto *ki_950 = buffer.data(ki + 950);
    const auto *ki_951 = buffer.data(ki + 951);
    const auto *ki_952 = buffer.data(ki + 952);
    const auto *ki_953 = buffer.data(ki + 953);
    const auto *ki_954 = buffer.data(ki + 954);
    const auto *ki_955 = buffer.data(ki + 955);
    const auto *ki_956 = buffer.data(ki + 956);
    const auto *ki_957 = buffer.data(ki + 957);
    const auto *ki_958 = buffer.data(ki + 958);
    const auto *ki_959 = buffer.data(ki + 959);
    const auto *ki_960 = buffer.data(ki + 960);
    const auto *ki_961 = buffer.data(ki + 961);
    const auto *ki_962 = buffer.data(ki + 962);
    const auto *ki_963 = buffer.data(ki + 963);
    const auto *ki_964 = buffer.data(ki + 964);
    const auto *ki_965 = buffer.data(ki + 965);
    const auto *ki_966 = buffer.data(ki + 966);
    const auto *ki_967 = buffer.data(ki + 967);
    const auto *ki_968 = buffer.data(ki + 968);
    const auto *ki_969 = buffer.data(ki + 969);
    const auto *ki_970 = buffer.data(ki + 970);
    const auto *ki_971 = buffer.data(ki + 971);
    const auto *ki_972 = buffer.data(ki + 972);
    const auto *ki_973 = buffer.data(ki + 973);
    const auto *ki_974 = buffer.data(ki + 974);
    const auto *ki_975 = buffer.data(ki + 975);
    const auto *ki_976 = buffer.data(ki + 976);
    const auto *ki_977 = buffer.data(ki + 977);
    const auto *ki_978 = buffer.data(ki + 978);
    const auto *ki_979 = buffer.data(ki + 979);
    const auto *ki_980 = buffer.data(ki + 980);
    const auto *ki_981 = buffer.data(ki + 981);
    const auto *ki_982 = buffer.data(ki + 982);
    const auto *ki_983 = buffer.data(ki + 983);
    const auto *ki_984 = buffer.data(ki + 984);
    const auto *ki_985 = buffer.data(ki + 985);
    const auto *ki_986 = buffer.data(ki + 986);
    const auto *ki_987 = buffer.data(ki + 987);
    const auto *ki_988 = buffer.data(ki + 988);
    const auto *ki_989 = buffer.data(ki + 989);
    const auto *ki_990 = buffer.data(ki + 990);
    const auto *ki_991 = buffer.data(ki + 991);
    const auto *ki_992 = buffer.data(ki + 992);
    const auto *ki_993 = buffer.data(ki + 993);
    const auto *ki_994 = buffer.data(ki + 994);
    const auto *ki_995 = buffer.data(ki + 995);
    const auto *ki_996 = buffer.data(ki + 996);
    const auto *ki_997 = buffer.data(ki + 997);
    const auto *ki_998 = buffer.data(ki + 998);
    const auto *ki_999 = buffer.data(ki + 999);
    const auto *ki_1000 = buffer.data(ki + 1000);
    const auto *ki_1001 = buffer.data(ki + 1001);
    const auto *ki_1002 = buffer.data(ki + 1002);
    const auto *ki_1003 = buffer.data(ki + 1003);
    const auto *ki_1004 = buffer.data(ki + 1004);
    const auto *ki_1005 = buffer.data(ki + 1005);
    const auto *ki_1006 = buffer.data(ki + 1006);
    const auto *ki_1007 = buffer.data(ki + 1007);

    const auto *mi_900 = buffer.data(mi + 900);
    const auto *mi_901 = buffer.data(mi + 901);
    const auto *mi_902 = buffer.data(mi + 902);
    const auto *mi_903 = buffer.data(mi + 903);
    const auto *mi_904 = buffer.data(mi + 904);
    const auto *mi_905 = buffer.data(mi + 905);
    const auto *mi_906 = buffer.data(mi + 906);
    const auto *mi_907 = buffer.data(mi + 907);
    const auto *mi_908 = buffer.data(mi + 908);
    const auto *mi_909 = buffer.data(mi + 909);
    const auto *mi_910 = buffer.data(mi + 910);
    const auto *mi_911 = buffer.data(mi + 911);
    const auto *mi_912 = buffer.data(mi + 912);
    const auto *mi_913 = buffer.data(mi + 913);
    const auto *mi_914 = buffer.data(mi + 914);
    const auto *mi_915 = buffer.data(mi + 915);
    const auto *mi_916 = buffer.data(mi + 916);
    const auto *mi_917 = buffer.data(mi + 917);
    const auto *mi_918 = buffer.data(mi + 918);
    const auto *mi_919 = buffer.data(mi + 919);
    const auto *mi_920 = buffer.data(mi + 920);
    const auto *mi_921 = buffer.data(mi + 921);
    const auto *mi_922 = buffer.data(mi + 922);
    const auto *mi_923 = buffer.data(mi + 923);
    const auto *mi_924 = buffer.data(mi + 924);
    const auto *mi_925 = buffer.data(mi + 925);
    const auto *mi_926 = buffer.data(mi + 926);
    const auto *mi_927 = buffer.data(mi + 927);
    const auto *mi_928 = buffer.data(mi + 928);
    const auto *mi_929 = buffer.data(mi + 929);
    const auto *mi_930 = buffer.data(mi + 930);
    const auto *mi_931 = buffer.data(mi + 931);
    const auto *mi_932 = buffer.data(mi + 932);
    const auto *mi_933 = buffer.data(mi + 933);
    const auto *mi_934 = buffer.data(mi + 934);
    const auto *mi_935 = buffer.data(mi + 935);
    const auto *mi_936 = buffer.data(mi + 936);
    const auto *mi_937 = buffer.data(mi + 937);
    const auto *mi_938 = buffer.data(mi + 938);
    const auto *mi_939 = buffer.data(mi + 939);
    const auto *mi_940 = buffer.data(mi + 940);
    const auto *mi_941 = buffer.data(mi + 941);
    const auto *mi_942 = buffer.data(mi + 942);
    const auto *mi_943 = buffer.data(mi + 943);
    const auto *mi_944 = buffer.data(mi + 944);
    const auto *mi_945 = buffer.data(mi + 945);
    const auto *mi_946 = buffer.data(mi + 946);
    const auto *mi_947 = buffer.data(mi + 947);
    const auto *mi_948 = buffer.data(mi + 948);
    const auto *mi_949 = buffer.data(mi + 949);
    const auto *mi_950 = buffer.data(mi + 950);
    const auto *mi_951 = buffer.data(mi + 951);
    const auto *mi_952 = buffer.data(mi + 952);
    const auto *mi_953 = buffer.data(mi + 953);
    const auto *mi_954 = buffer.data(mi + 954);
    const auto *mi_955 = buffer.data(mi + 955);
    const auto *mi_956 = buffer.data(mi + 956);
    const auto *mi_957 = buffer.data(mi + 957);
    const auto *mi_958 = buffer.data(mi + 958);
    const auto *mi_959 = buffer.data(mi + 959);
    const auto *mi_960 = buffer.data(mi + 960);
    const auto *mi_961 = buffer.data(mi + 961);
    const auto *mi_962 = buffer.data(mi + 962);
    const auto *mi_963 = buffer.data(mi + 963);
    const auto *mi_964 = buffer.data(mi + 964);
    const auto *mi_965 = buffer.data(mi + 965);
    const auto *mi_966 = buffer.data(mi + 966);
    const auto *mi_967 = buffer.data(mi + 967);
    const auto *mi_968 = buffer.data(mi + 968);
    const auto *mi_969 = buffer.data(mi + 969);
    const auto *mi_970 = buffer.data(mi + 970);
    const auto *mi_971 = buffer.data(mi + 971);
    const auto *mi_972 = buffer.data(mi + 972);
    const auto *mi_973 = buffer.data(mi + 973);
    const auto *mi_974 = buffer.data(mi + 974);
    const auto *mi_975 = buffer.data(mi + 975);
    const auto *mi_976 = buffer.data(mi + 976);
    const auto *mi_977 = buffer.data(mi + 977);
    const auto *mi_978 = buffer.data(mi + 978);
    const auto *mi_979 = buffer.data(mi + 979);
    const auto *mi_980 = buffer.data(mi + 980);
    const auto *mi_981 = buffer.data(mi + 981);
    const auto *mi_982 = buffer.data(mi + 982);
    const auto *mi_983 = buffer.data(mi + 983);
    const auto *mi_984 = buffer.data(mi + 984);
    const auto *mi_985 = buffer.data(mi + 985);
    const auto *mi_986 = buffer.data(mi + 986);
    const auto *mi_987 = buffer.data(mi + 987);
    const auto *mi_988 = buffer.data(mi + 988);
    const auto *mi_989 = buffer.data(mi + 989);
    const auto *mi_990 = buffer.data(mi + 990);
    const auto *mi_991 = buffer.data(mi + 991);
    const auto *mi_992 = buffer.data(mi + 992);
    const auto *mi_993 = buffer.data(mi + 993);
    const auto *mi_994 = buffer.data(mi + 994);
    const auto *mi_995 = buffer.data(mi + 995);
    const auto *mi_996 = buffer.data(mi + 996);
    const auto *mi_997 = buffer.data(mi + 997);
    const auto *mi_998 = buffer.data(mi + 998);
    const auto *mi_999 = buffer.data(mi + 999);
    const auto *mi_1000 = buffer.data(mi + 1000);
    const auto *mi_1001 = buffer.data(mi + 1001);
    const auto *mi_1002 = buffer.data(mi + 1002);
    const auto *mi_1003 = buffer.data(mi + 1003);
    const auto *mi_1004 = buffer.data(mi + 1004);
    const auto *mi_1005 = buffer.data(mi + 1005);
    const auto *mi_1006 = buffer.data(mi + 1006);
    const auto *mi_1007 = buffer.data(mi + 1007);
    const auto *mi_1008 = buffer.data(mi + 1008);
    const auto *mi_1009 = buffer.data(mi + 1009);
    const auto *mi_1010 = buffer.data(mi + 1010);
    const auto *mi_1011 = buffer.data(mi + 1011);
    const auto *mi_1012 = buffer.data(mi + 1012);
    const auto *mi_1013 = buffer.data(mi + 1013);
    const auto *mi_1014 = buffer.data(mi + 1014);
    const auto *mi_1015 = buffer.data(mi + 1015);
    const auto *mi_1016 = buffer.data(mi + 1016);
    const auto *mi_1017 = buffer.data(mi + 1017);
    const auto *mi_1018 = buffer.data(mi + 1018);
    const auto *mi_1019 = buffer.data(mi + 1019);
    const auto *mi_1020 = buffer.data(mi + 1020);
    const auto *mi_1021 = buffer.data(mi + 1021);
    const auto *mi_1022 = buffer.data(mi + 1022);
    const auto *mi_1023 = buffer.data(mi + 1023);
    const auto *mi_1024 = buffer.data(mi + 1024);
    const auto *mi_1025 = buffer.data(mi + 1025);
    const auto *mi_1026 = buffer.data(mi + 1026);
    const auto *mi_1027 = buffer.data(mi + 1027);
    const auto *mi_1028 = buffer.data(mi + 1028);
    const auto *mi_1029 = buffer.data(mi + 1029);
    const auto *mi_1030 = buffer.data(mi + 1030);
    const auto *mi_1031 = buffer.data(mi + 1031);
    const auto *mi_1032 = buffer.data(mi + 1032);
    const auto *mi_1033 = buffer.data(mi + 1033);
    const auto *mi_1034 = buffer.data(mi + 1034);
    const auto *mi_1035 = buffer.data(mi + 1035);
    const auto *mi_1036 = buffer.data(mi + 1036);
    const auto *mi_1037 = buffer.data(mi + 1037);
    const auto *mi_1038 = buffer.data(mi + 1038);
    const auto *mi_1039 = buffer.data(mi + 1039);
    const auto *mi_1040 = buffer.data(mi + 1040);
    const auto *mi_1041 = buffer.data(mi + 1041);
    const auto *mi_1042 = buffer.data(mi + 1042);
    const auto *mi_1043 = buffer.data(mi + 1043);
    const auto *mi_1044 = buffer.data(mi + 1044);
    const auto *mi_1045 = buffer.data(mi + 1045);
    const auto *mi_1046 = buffer.data(mi + 1046);
    const auto *mi_1047 = buffer.data(mi + 1047);
    const auto *mi_1048 = buffer.data(mi + 1048);
    const auto *mi_1049 = buffer.data(mi + 1049);
    const auto *mi_1050 = buffer.data(mi + 1050);
    const auto *mi_1051 = buffer.data(mi + 1051);
    const auto *mi_1052 = buffer.data(mi + 1052);
    const auto *mi_1053 = buffer.data(mi + 1053);
    const auto *mi_1054 = buffer.data(mi + 1054);
    const auto *mi_1055 = buffer.data(mi + 1055);
    const auto *mi_1056 = buffer.data(mi + 1056);
    const auto *mi_1057 = buffer.data(mi + 1057);
    const auto *mi_1058 = buffer.data(mi + 1058);
    const auto *mi_1059 = buffer.data(mi + 1059);
    const auto *mi_1060 = buffer.data(mi + 1060);
    const auto *mi_1061 = buffer.data(mi + 1061);
    const auto *mi_1062 = buffer.data(mi + 1062);
    const auto *mi_1063 = buffer.data(mi + 1063);
    const auto *mi_1064 = buffer.data(mi + 1064);
    const auto *mi_1065 = buffer.data(mi + 1065);
    const auto *mi_1066 = buffer.data(mi + 1066);

#pragma omp simd aligned(t_900, t_901, t_902, t_903, t_904, ki_900, ki_901, ki_902, ki_903, \
                         ki_904, mi_900, mi_901, mi_902, mi_903, \
                         mi_904 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = -ki_900[k]
                   + f_0 * mi_900[k];

        t_901[k] = -ki_901[k]
                   + f_0 * mi_901[k];

        t_902[k] = -ki_902[k]
                   + f_0 * mi_902[k];

        t_903[k] = -ki_903[k]
                   + f_0 * mi_903[k];

        t_904[k] = -ki_904[k]
                   + f_0 * mi_904[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, t_909, ki_905, ki_906, ki_907, ki_908, \
                         ki_909, mi_905, mi_906, mi_907, mi_908, \
                         mi_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = -ki_905[k]
                   + f_0 * mi_905[k];

        t_906[k] = -ki_906[k]
                   + f_0 * mi_906[k];

        t_907[k] = -ki_907[k]
                   + f_0 * mi_907[k];

        t_908[k] = -ki_908[k]
                   + f_0 * mi_908[k];

        t_909[k] = -ki_909[k]
                   + f_0 * mi_909[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, t_914, ki_910, ki_911, ki_912, ki_913, \
                         ki_914, mi_910, mi_911, mi_912, mi_913, \
                         mi_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = -ki_910[k]
                   + f_0 * mi_910[k];

        t_911[k] = -ki_911[k]
                   + f_0 * mi_911[k];

        t_912[k] = -ki_912[k]
                   + f_0 * mi_912[k];

        t_913[k] = -ki_913[k]
                   + f_0 * mi_913[k];

        t_914[k] = -ki_914[k]
                   + f_0 * mi_914[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, t_919, ki_915, ki_916, ki_917, ki_918, \
                         ki_919, mi_915, mi_916, mi_917, mi_918, \
                         mi_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = -ki_915[k]
                   + f_0 * mi_915[k];

        t_916[k] = -ki_916[k]
                   + f_0 * mi_916[k];

        t_917[k] = -ki_917[k]
                   + f_0 * mi_917[k];

        t_918[k] = -ki_918[k]
                   + f_0 * mi_918[k];

        t_919[k] = -ki_919[k]
                   + f_0 * mi_919[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, t_924, ki_920, ki_921, ki_922, ki_923, \
                         ki_924, mi_920, mi_921, mi_922, mi_923, \
                         mi_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = -ki_920[k]
                   + f_0 * mi_920[k];

        t_921[k] = -ki_921[k]
                   + f_0 * mi_921[k];

        t_922[k] = -ki_922[k]
                   + f_0 * mi_922[k];

        t_923[k] = -ki_923[k]
                   + f_0 * mi_923[k];

        t_924[k] = -ki_924[k]
                   + f_0 * mi_924[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, t_929, ki_925, ki_926, ki_927, ki_928, \
                         ki_929, mi_925, mi_926, mi_927, mi_928, \
                         mi_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = -ki_925[k]
                   + f_0 * mi_925[k];

        t_926[k] = -ki_926[k]
                   + f_0 * mi_926[k];

        t_927[k] = -ki_927[k]
                   + f_0 * mi_927[k];

        t_928[k] = -ki_928[k]
                   + f_0 * mi_928[k];

        t_929[k] = -ki_929[k]
                   + f_0 * mi_929[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, ki_930, ki_931, ki_932, ki_933, \
                         ki_934, mi_930, mi_931, mi_932, mi_933, \
                         mi_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = -ki_930[k]
                   + f_0 * mi_930[k];

        t_931[k] = -ki_931[k]
                   + f_0 * mi_931[k];

        t_932[k] = -ki_932[k]
                   + f_0 * mi_932[k];

        t_933[k] = -ki_933[k]
                   + f_0 * mi_933[k];

        t_934[k] = -ki_934[k]
                   + f_0 * mi_934[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, t_939, ki_935, ki_936, ki_937, ki_938, \
                         ki_939, mi_935, mi_936, mi_937, mi_938, \
                         mi_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = -ki_935[k]
                   + f_0 * mi_935[k];

        t_936[k] = -ki_936[k]
                   + f_0 * mi_936[k];

        t_937[k] = -ki_937[k]
                   + f_0 * mi_937[k];

        t_938[k] = -ki_938[k]
                   + f_0 * mi_938[k];

        t_939[k] = -ki_939[k]
                   + f_0 * mi_939[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, t_944, ki_940, ki_941, ki_942, ki_943, \
                         ki_944, mi_940, mi_941, mi_942, mi_943, \
                         mi_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = -ki_940[k]
                   + f_0 * mi_940[k];

        t_941[k] = -ki_941[k]
                   + f_0 * mi_941[k];

        t_942[k] = -ki_942[k]
                   + f_0 * mi_942[k];

        t_943[k] = -ki_943[k]
                   + f_0 * mi_943[k];

        t_944[k] = -ki_944[k]
                   + f_0 * mi_944[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, t_949, ki_945, ki_946, ki_947, ki_948, \
                         ki_949, mi_945, mi_946, mi_947, mi_948, \
                         mi_949 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = -ki_945[k]
                   + f_0 * mi_945[k];

        t_946[k] = -ki_946[k]
                   + f_0 * mi_946[k];

        t_947[k] = -ki_947[k]
                   + f_0 * mi_947[k];

        t_948[k] = -ki_948[k]
                   + f_0 * mi_948[k];

        t_949[k] = -ki_949[k]
                   + f_0 * mi_949[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, t_954, ki_950, ki_951, ki_952, ki_953, \
                         ki_954, mi_950, mi_951, mi_952, mi_953, \
                         mi_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = -ki_950[k]
                   + f_0 * mi_950[k];

        t_951[k] = -ki_951[k]
                   + f_0 * mi_951[k];

        t_952[k] = -ki_952[k]
                   + f_0 * mi_952[k];

        t_953[k] = -ki_953[k]
                   + f_0 * mi_953[k];

        t_954[k] = -ki_954[k]
                   + f_0 * mi_954[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, t_959, ki_955, ki_956, ki_957, ki_958, \
                         ki_959, mi_955, mi_956, mi_957, mi_958, \
                         mi_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = -ki_955[k]
                   + f_0 * mi_955[k];

        t_956[k] = -ki_956[k]
                   + f_0 * mi_956[k];

        t_957[k] = -ki_957[k]
                   + f_0 * mi_957[k];

        t_958[k] = -ki_958[k]
                   + f_0 * mi_958[k];

        t_959[k] = -ki_959[k]
                   + f_0 * mi_959[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, t_963, t_964, ki_960, ki_961, ki_962, ki_963, \
                         ki_964, mi_960, mi_961, mi_962, mi_963, \
                         mi_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = -ki_960[k]
                   + f_0 * mi_960[k];

        t_961[k] = -ki_961[k]
                   + f_0 * mi_961[k];

        t_962[k] = -ki_962[k]
                   + f_0 * mi_962[k];

        t_963[k] = -ki_963[k]
                   + f_0 * mi_963[k];

        t_964[k] = -ki_964[k]
                   + f_0 * mi_964[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, t_969, ki_965, ki_966, ki_967, ki_968, \
                         ki_969, mi_965, mi_966, mi_967, mi_968, \
                         mi_969 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = -ki_965[k]
                   + f_0 * mi_965[k];

        t_966[k] = -ki_966[k]
                   + f_0 * mi_966[k];

        t_967[k] = -ki_967[k]
                   + f_0 * mi_967[k];

        t_968[k] = -ki_968[k]
                   + f_0 * mi_968[k];

        t_969[k] = -ki_969[k]
                   + f_0 * mi_969[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, t_974, ki_970, ki_971, ki_972, ki_973, \
                         ki_974, mi_970, mi_971, mi_972, mi_973, \
                         mi_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = -ki_970[k]
                   + f_0 * mi_970[k];

        t_971[k] = -ki_971[k]
                   + f_0 * mi_971[k];

        t_972[k] = -ki_972[k]
                   + f_0 * mi_972[k];

        t_973[k] = -ki_973[k]
                   + f_0 * mi_973[k];

        t_974[k] = -ki_974[k]
                   + f_0 * mi_974[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, t_979, ki_975, ki_976, ki_977, ki_978, \
                         ki_979, mi_975, mi_976, mi_977, mi_978, \
                         mi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = -ki_975[k]
                   + f_0 * mi_975[k];

        t_976[k] = -ki_976[k]
                   + f_0 * mi_976[k];

        t_977[k] = -ki_977[k]
                   + f_0 * mi_977[k];

        t_978[k] = -ki_978[k]
                   + f_0 * mi_978[k];

        t_979[k] = -ki_979[k]
                   + f_0 * mi_979[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, t_984, ki_980, ki_981, ki_982, ki_983, \
                         ki_984, mi_980, mi_981, mi_982, mi_983, \
                         mi_984 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = -ki_980[k]
                   + f_0 * mi_980[k];

        t_981[k] = -ki_981[k]
                   + f_0 * mi_981[k];

        t_982[k] = -ki_982[k]
                   + f_0 * mi_982[k];

        t_983[k] = -ki_983[k]
                   + f_0 * mi_983[k];

        t_984[k] = -ki_984[k]
                   + f_0 * mi_984[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, t_989, ki_985, ki_986, ki_987, ki_988, \
                         ki_989, mi_985, mi_986, mi_987, mi_988, \
                         mi_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = -ki_985[k]
                   + f_0 * mi_985[k];

        t_986[k] = -ki_986[k]
                   + f_0 * mi_986[k];

        t_987[k] = -ki_987[k]
                   + f_0 * mi_987[k];

        t_988[k] = -ki_988[k]
                   + f_0 * mi_988[k];

        t_989[k] = -ki_989[k]
                   + f_0 * mi_989[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, t_993, t_994, ki_990, ki_991, ki_992, ki_993, \
                         ki_994, mi_990, mi_991, mi_992, mi_993, \
                         mi_994 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = -ki_990[k]
                   + f_0 * mi_990[k];

        t_991[k] = -ki_991[k]
                   + f_0 * mi_991[k];

        t_992[k] = -ki_992[k]
                   + f_0 * mi_992[k];

        t_993[k] = -ki_993[k]
                   + f_0 * mi_993[k];

        t_994[k] = -ki_994[k]
                   + f_0 * mi_994[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, ki_995, ki_996, ki_997, ki_998, \
                         ki_999, mi_995, mi_996, mi_997, mi_998, \
                         mi_999 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = -ki_995[k]
                   + f_0 * mi_995[k];

        t_996[k] = -ki_996[k]
                   + f_0 * mi_996[k];

        t_997[k] = -ki_997[k]
                   + f_0 * mi_997[k];

        t_998[k] = -ki_998[k]
                   + f_0 * mi_998[k];

        t_999[k] = -ki_999[k]
                   + f_0 * mi_999[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, t_1004, ki_1000, ki_1001, ki_1002, \
                         ki_1003, ki_1004, mi_1000, mi_1001, mi_1002, mi_1003, \
                         mi_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = -ki_1000[k]
                    + f_0 * mi_1000[k];

        t_1001[k] = -ki_1001[k]
                    + f_0 * mi_1001[k];

        t_1002[k] = -ki_1002[k]
                    + f_0 * mi_1002[k];

        t_1003[k] = -ki_1003[k]
                    + f_0 * mi_1003[k];

        t_1004[k] = -ki_1004[k]
                    + f_0 * mi_1004[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, t_1008, t_1009, t_1010, ki_1005, ki_1006, \
                         ki_1007, mi_1005, mi_1006, mi_1007, mi_1008, mi_1009, \
                         mi_1010 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = -ki_1005[k]
                    + f_0 * mi_1005[k];

        t_1006[k] = -ki_1006[k]
                    + f_0 * mi_1006[k];

        t_1007[k] = -ki_1007[k]
                    + f_0 * mi_1007[k];

        t_1008[k] = f_0 * mi_1008[k];

        t_1009[k] = f_0 * mi_1009[k];

        t_1010[k] = f_0 * mi_1010[k];
    }

#pragma omp simd aligned(t_1011, t_1012, t_1013, t_1014, t_1015, t_1016, t_1017, t_1018, \
                         mi_1011, mi_1012, mi_1013, mi_1014, mi_1015, mi_1016, mi_1017, \
                         mi_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1011[k] = f_0 * mi_1011[k];

        t_1012[k] = f_0 * mi_1012[k];

        t_1013[k] = f_0 * mi_1013[k];

        t_1014[k] = f_0 * mi_1014[k];

        t_1015[k] = f_0 * mi_1015[k];

        t_1016[k] = f_0 * mi_1016[k];

        t_1017[k] = f_0 * mi_1017[k];

        t_1018[k] = f_0 * mi_1018[k];
    }

#pragma omp simd aligned(t_1019, t_1020, t_1021, t_1022, t_1023, t_1024, t_1025, t_1026, \
                         mi_1019, mi_1020, mi_1021, mi_1022, mi_1023, mi_1024, mi_1025, \
                         mi_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1019[k] = f_0 * mi_1019[k];

        t_1020[k] = f_0 * mi_1020[k];

        t_1021[k] = f_0 * mi_1021[k];

        t_1022[k] = f_0 * mi_1022[k];

        t_1023[k] = f_0 * mi_1023[k];

        t_1024[k] = f_0 * mi_1024[k];

        t_1025[k] = f_0 * mi_1025[k];

        t_1026[k] = f_0 * mi_1026[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, t_1031, t_1032, t_1033, t_1034, \
                         mi_1027, mi_1028, mi_1029, mi_1030, mi_1031, mi_1032, mi_1033, \
                         mi_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_0 * mi_1027[k];

        t_1028[k] = f_0 * mi_1028[k];

        t_1029[k] = f_0 * mi_1029[k];

        t_1030[k] = f_0 * mi_1030[k];

        t_1031[k] = f_0 * mi_1031[k];

        t_1032[k] = f_0 * mi_1032[k];

        t_1033[k] = f_0 * mi_1033[k];

        t_1034[k] = f_0 * mi_1034[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, t_1038, t_1039, t_1040, t_1041, t_1042, \
                         mi_1035, mi_1036, mi_1037, mi_1038, mi_1039, mi_1040, mi_1041, \
                         mi_1042 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = f_0 * mi_1035[k];

        t_1036[k] = f_0 * mi_1036[k];

        t_1037[k] = f_0 * mi_1037[k];

        t_1038[k] = f_0 * mi_1038[k];

        t_1039[k] = f_0 * mi_1039[k];

        t_1040[k] = f_0 * mi_1040[k];

        t_1041[k] = f_0 * mi_1041[k];

        t_1042[k] = f_0 * mi_1042[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, t_1046, t_1047, t_1048, t_1049, t_1050, \
                         mi_1043, mi_1044, mi_1045, mi_1046, mi_1047, mi_1048, mi_1049, \
                         mi_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_0 * mi_1043[k];

        t_1044[k] = f_0 * mi_1044[k];

        t_1045[k] = f_0 * mi_1045[k];

        t_1046[k] = f_0 * mi_1046[k];

        t_1047[k] = f_0 * mi_1047[k];

        t_1048[k] = f_0 * mi_1048[k];

        t_1049[k] = f_0 * mi_1049[k];

        t_1050[k] = f_0 * mi_1050[k];
    }

#pragma omp simd aligned(t_1051, t_1052, t_1053, t_1054, t_1055, t_1056, t_1057, t_1058, \
                         mi_1051, mi_1052, mi_1053, mi_1054, mi_1055, mi_1056, mi_1057, \
                         mi_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1051[k] = f_0 * mi_1051[k];

        t_1052[k] = f_0 * mi_1052[k];

        t_1053[k] = f_0 * mi_1053[k];

        t_1054[k] = f_0 * mi_1054[k];

        t_1055[k] = f_0 * mi_1055[k];

        t_1056[k] = f_0 * mi_1056[k];

        t_1057[k] = f_0 * mi_1057[k];

        t_1058[k] = f_0 * mi_1058[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, t_1062, t_1063, t_1064, t_1065, t_1066, \
                         mi_1059, mi_1060, mi_1061, mi_1062, mi_1063, mi_1064, mi_1065, \
                         mi_1066 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = f_0 * mi_1059[k];

        t_1060[k] = f_0 * mi_1060[k];

        t_1061[k] = f_0 * mi_1061[k];

        t_1062[k] = f_0 * mi_1062[k];

        t_1063[k] = f_0 * mi_1063[k];

        t_1064[k] = f_0 * mi_1064[k];

        t_1065[k] = f_0 * mi_1065[k];

        t_1066[k] = f_0 * mi_1066[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_0_piece7(CSimdMatrix &buffer, const size_t target,
                                                    const size_t mi, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_1067 = buffer.data(target + 1067);
    auto *t_1068 = buffer.data(target + 1068);
    auto *t_1069 = buffer.data(target + 1069);
    auto *t_1070 = buffer.data(target + 1070);
    auto *t_1071 = buffer.data(target + 1071);
    auto *t_1072 = buffer.data(target + 1072);
    auto *t_1073 = buffer.data(target + 1073);
    auto *t_1074 = buffer.data(target + 1074);
    auto *t_1075 = buffer.data(target + 1075);
    auto *t_1076 = buffer.data(target + 1076);
    auto *t_1077 = buffer.data(target + 1077);
    auto *t_1078 = buffer.data(target + 1078);
    auto *t_1079 = buffer.data(target + 1079);
    auto *t_1080 = buffer.data(target + 1080);
    auto *t_1081 = buffer.data(target + 1081);
    auto *t_1082 = buffer.data(target + 1082);
    auto *t_1083 = buffer.data(target + 1083);
    auto *t_1084 = buffer.data(target + 1084);
    auto *t_1085 = buffer.data(target + 1085);
    auto *t_1086 = buffer.data(target + 1086);
    auto *t_1087 = buffer.data(target + 1087);
    auto *t_1088 = buffer.data(target + 1088);
    auto *t_1089 = buffer.data(target + 1089);
    auto *t_1090 = buffer.data(target + 1090);
    auto *t_1091 = buffer.data(target + 1091);
    auto *t_1092 = buffer.data(target + 1092);
    auto *t_1093 = buffer.data(target + 1093);
    auto *t_1094 = buffer.data(target + 1094);
    auto *t_1095 = buffer.data(target + 1095);
    auto *t_1096 = buffer.data(target + 1096);
    auto *t_1097 = buffer.data(target + 1097);
    auto *t_1098 = buffer.data(target + 1098);
    auto *t_1099 = buffer.data(target + 1099);
    auto *t_1100 = buffer.data(target + 1100);
    auto *t_1101 = buffer.data(target + 1101);
    auto *t_1102 = buffer.data(target + 1102);
    auto *t_1103 = buffer.data(target + 1103);
    auto *t_1104 = buffer.data(target + 1104);
    auto *t_1105 = buffer.data(target + 1105);
    auto *t_1106 = buffer.data(target + 1106);
    auto *t_1107 = buffer.data(target + 1107);
    auto *t_1108 = buffer.data(target + 1108);
    auto *t_1109 = buffer.data(target + 1109);
    auto *t_1110 = buffer.data(target + 1110);
    auto *t_1111 = buffer.data(target + 1111);
    auto *t_1112 = buffer.data(target + 1112);
    auto *t_1113 = buffer.data(target + 1113);
    auto *t_1114 = buffer.data(target + 1114);
    auto *t_1115 = buffer.data(target + 1115);
    auto *t_1116 = buffer.data(target + 1116);
    auto *t_1117 = buffer.data(target + 1117);
    auto *t_1118 = buffer.data(target + 1118);
    auto *t_1119 = buffer.data(target + 1119);
    auto *t_1120 = buffer.data(target + 1120);
    auto *t_1121 = buffer.data(target + 1121);
    auto *t_1122 = buffer.data(target + 1122);
    auto *t_1123 = buffer.data(target + 1123);
    auto *t_1124 = buffer.data(target + 1124);
    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);
    auto *t_1133 = buffer.data(target + 1133);
    auto *t_1134 = buffer.data(target + 1134);
    auto *t_1135 = buffer.data(target + 1135);
    auto *t_1136 = buffer.data(target + 1136);
    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);
    auto *t_1169 = buffer.data(target + 1169);
    auto *t_1170 = buffer.data(target + 1170);
    auto *t_1171 = buffer.data(target + 1171);
    auto *t_1172 = buffer.data(target + 1172);
    auto *t_1173 = buffer.data(target + 1173);
    auto *t_1174 = buffer.data(target + 1174);
    auto *t_1175 = buffer.data(target + 1175);
    auto *t_1176 = buffer.data(target + 1176);
    auto *t_1177 = buffer.data(target + 1177);
    auto *t_1178 = buffer.data(target + 1178);
    auto *t_1179 = buffer.data(target + 1179);
    auto *t_1180 = buffer.data(target + 1180);
    auto *t_1181 = buffer.data(target + 1181);
    auto *t_1182 = buffer.data(target + 1182);
    auto *t_1183 = buffer.data(target + 1183);
    auto *t_1184 = buffer.data(target + 1184);
    auto *t_1185 = buffer.data(target + 1185);
    auto *t_1186 = buffer.data(target + 1186);
    auto *t_1187 = buffer.data(target + 1187);
    auto *t_1188 = buffer.data(target + 1188);
    auto *t_1189 = buffer.data(target + 1189);
    auto *t_1190 = buffer.data(target + 1190);
    auto *t_1191 = buffer.data(target + 1191);
    auto *t_1192 = buffer.data(target + 1192);
    auto *t_1193 = buffer.data(target + 1193);
    auto *t_1194 = buffer.data(target + 1194);
    auto *t_1195 = buffer.data(target + 1195);
    auto *t_1196 = buffer.data(target + 1196);
    auto *t_1197 = buffer.data(target + 1197);
    auto *t_1198 = buffer.data(target + 1198);
    auto *t_1199 = buffer.data(target + 1199);
    auto *t_1200 = buffer.data(target + 1200);
    auto *t_1201 = buffer.data(target + 1201);
    auto *t_1202 = buffer.data(target + 1202);
    auto *t_1203 = buffer.data(target + 1203);
    auto *t_1204 = buffer.data(target + 1204);
    auto *t_1205 = buffer.data(target + 1205);
    auto *t_1206 = buffer.data(target + 1206);
    auto *t_1207 = buffer.data(target + 1207);
    auto *t_1208 = buffer.data(target + 1208);
    auto *t_1209 = buffer.data(target + 1209);
    auto *t_1210 = buffer.data(target + 1210);
    auto *t_1211 = buffer.data(target + 1211);
    auto *t_1212 = buffer.data(target + 1212);
    auto *t_1213 = buffer.data(target + 1213);
    auto *t_1214 = buffer.data(target + 1214);
    auto *t_1215 = buffer.data(target + 1215);
    auto *t_1216 = buffer.data(target + 1216);
    auto *t_1217 = buffer.data(target + 1217);
    auto *t_1218 = buffer.data(target + 1218);
    auto *t_1219 = buffer.data(target + 1219);
    auto *t_1220 = buffer.data(target + 1220);
    auto *t_1221 = buffer.data(target + 1221);
    auto *t_1222 = buffer.data(target + 1222);
    auto *t_1223 = buffer.data(target + 1223);
    auto *t_1224 = buffer.data(target + 1224);
    auto *t_1225 = buffer.data(target + 1225);
    auto *t_1226 = buffer.data(target + 1226);
    auto *t_1227 = buffer.data(target + 1227);
    auto *t_1228 = buffer.data(target + 1228);
    auto *t_1229 = buffer.data(target + 1229);
    auto *t_1230 = buffer.data(target + 1230);
    auto *t_1231 = buffer.data(target + 1231);
    auto *t_1232 = buffer.data(target + 1232);
    auto *t_1233 = buffer.data(target + 1233);
    auto *t_1234 = buffer.data(target + 1234);
    auto *t_1235 = buffer.data(target + 1235);
    auto *t_1236 = buffer.data(target + 1236);
    auto *t_1237 = buffer.data(target + 1237);
    auto *t_1238 = buffer.data(target + 1238);
    auto *t_1239 = buffer.data(target + 1239);
    auto *t_1240 = buffer.data(target + 1240);
    auto *t_1241 = buffer.data(target + 1241);
    auto *t_1242 = buffer.data(target + 1242);
    auto *t_1243 = buffer.data(target + 1243);
    auto *t_1244 = buffer.data(target + 1244);
    auto *t_1245 = buffer.data(target + 1245);
    auto *t_1246 = buffer.data(target + 1246);
    auto *t_1247 = buffer.data(target + 1247);
    auto *t_1248 = buffer.data(target + 1248);
    auto *t_1249 = buffer.data(target + 1249);
    auto *t_1250 = buffer.data(target + 1250);
    auto *t_1251 = buffer.data(target + 1251);
    auto *t_1252 = buffer.data(target + 1252);
    auto *t_1253 = buffer.data(target + 1253);
    auto *t_1254 = buffer.data(target + 1254);
    auto *t_1255 = buffer.data(target + 1255);
    auto *t_1256 = buffer.data(target + 1256);
    auto *t_1257 = buffer.data(target + 1257);
    auto *t_1258 = buffer.data(target + 1258);
    auto *t_1259 = buffer.data(target + 1259);

    const auto *mi_1067 = buffer.data(mi + 1067);
    const auto *mi_1068 = buffer.data(mi + 1068);
    const auto *mi_1069 = buffer.data(mi + 1069);
    const auto *mi_1070 = buffer.data(mi + 1070);
    const auto *mi_1071 = buffer.data(mi + 1071);
    const auto *mi_1072 = buffer.data(mi + 1072);
    const auto *mi_1073 = buffer.data(mi + 1073);
    const auto *mi_1074 = buffer.data(mi + 1074);
    const auto *mi_1075 = buffer.data(mi + 1075);
    const auto *mi_1076 = buffer.data(mi + 1076);
    const auto *mi_1077 = buffer.data(mi + 1077);
    const auto *mi_1078 = buffer.data(mi + 1078);
    const auto *mi_1079 = buffer.data(mi + 1079);
    const auto *mi_1080 = buffer.data(mi + 1080);
    const auto *mi_1081 = buffer.data(mi + 1081);
    const auto *mi_1082 = buffer.data(mi + 1082);
    const auto *mi_1083 = buffer.data(mi + 1083);
    const auto *mi_1084 = buffer.data(mi + 1084);
    const auto *mi_1085 = buffer.data(mi + 1085);
    const auto *mi_1086 = buffer.data(mi + 1086);
    const auto *mi_1087 = buffer.data(mi + 1087);
    const auto *mi_1088 = buffer.data(mi + 1088);
    const auto *mi_1089 = buffer.data(mi + 1089);
    const auto *mi_1090 = buffer.data(mi + 1090);
    const auto *mi_1091 = buffer.data(mi + 1091);
    const auto *mi_1092 = buffer.data(mi + 1092);
    const auto *mi_1093 = buffer.data(mi + 1093);
    const auto *mi_1094 = buffer.data(mi + 1094);
    const auto *mi_1095 = buffer.data(mi + 1095);
    const auto *mi_1096 = buffer.data(mi + 1096);
    const auto *mi_1097 = buffer.data(mi + 1097);
    const auto *mi_1098 = buffer.data(mi + 1098);
    const auto *mi_1099 = buffer.data(mi + 1099);
    const auto *mi_1100 = buffer.data(mi + 1100);
    const auto *mi_1101 = buffer.data(mi + 1101);
    const auto *mi_1102 = buffer.data(mi + 1102);
    const auto *mi_1103 = buffer.data(mi + 1103);
    const auto *mi_1104 = buffer.data(mi + 1104);
    const auto *mi_1105 = buffer.data(mi + 1105);
    const auto *mi_1106 = buffer.data(mi + 1106);
    const auto *mi_1107 = buffer.data(mi + 1107);
    const auto *mi_1108 = buffer.data(mi + 1108);
    const auto *mi_1109 = buffer.data(mi + 1109);
    const auto *mi_1110 = buffer.data(mi + 1110);
    const auto *mi_1111 = buffer.data(mi + 1111);
    const auto *mi_1112 = buffer.data(mi + 1112);
    const auto *mi_1113 = buffer.data(mi + 1113);
    const auto *mi_1114 = buffer.data(mi + 1114);
    const auto *mi_1115 = buffer.data(mi + 1115);
    const auto *mi_1116 = buffer.data(mi + 1116);
    const auto *mi_1117 = buffer.data(mi + 1117);
    const auto *mi_1118 = buffer.data(mi + 1118);
    const auto *mi_1119 = buffer.data(mi + 1119);
    const auto *mi_1120 = buffer.data(mi + 1120);
    const auto *mi_1121 = buffer.data(mi + 1121);
    const auto *mi_1122 = buffer.data(mi + 1122);
    const auto *mi_1123 = buffer.data(mi + 1123);
    const auto *mi_1124 = buffer.data(mi + 1124);
    const auto *mi_1125 = buffer.data(mi + 1125);
    const auto *mi_1126 = buffer.data(mi + 1126);
    const auto *mi_1127 = buffer.data(mi + 1127);
    const auto *mi_1128 = buffer.data(mi + 1128);
    const auto *mi_1129 = buffer.data(mi + 1129);
    const auto *mi_1130 = buffer.data(mi + 1130);
    const auto *mi_1131 = buffer.data(mi + 1131);
    const auto *mi_1132 = buffer.data(mi + 1132);
    const auto *mi_1133 = buffer.data(mi + 1133);
    const auto *mi_1134 = buffer.data(mi + 1134);
    const auto *mi_1135 = buffer.data(mi + 1135);
    const auto *mi_1136 = buffer.data(mi + 1136);
    const auto *mi_1137 = buffer.data(mi + 1137);
    const auto *mi_1138 = buffer.data(mi + 1138);
    const auto *mi_1139 = buffer.data(mi + 1139);
    const auto *mi_1140 = buffer.data(mi + 1140);
    const auto *mi_1141 = buffer.data(mi + 1141);
    const auto *mi_1142 = buffer.data(mi + 1142);
    const auto *mi_1143 = buffer.data(mi + 1143);
    const auto *mi_1144 = buffer.data(mi + 1144);
    const auto *mi_1145 = buffer.data(mi + 1145);
    const auto *mi_1146 = buffer.data(mi + 1146);
    const auto *mi_1147 = buffer.data(mi + 1147);
    const auto *mi_1148 = buffer.data(mi + 1148);
    const auto *mi_1149 = buffer.data(mi + 1149);
    const auto *mi_1150 = buffer.data(mi + 1150);
    const auto *mi_1151 = buffer.data(mi + 1151);
    const auto *mi_1152 = buffer.data(mi + 1152);
    const auto *mi_1153 = buffer.data(mi + 1153);
    const auto *mi_1154 = buffer.data(mi + 1154);
    const auto *mi_1155 = buffer.data(mi + 1155);
    const auto *mi_1156 = buffer.data(mi + 1156);
    const auto *mi_1157 = buffer.data(mi + 1157);
    const auto *mi_1158 = buffer.data(mi + 1158);
    const auto *mi_1159 = buffer.data(mi + 1159);
    const auto *mi_1160 = buffer.data(mi + 1160);
    const auto *mi_1161 = buffer.data(mi + 1161);
    const auto *mi_1162 = buffer.data(mi + 1162);
    const auto *mi_1163 = buffer.data(mi + 1163);
    const auto *mi_1164 = buffer.data(mi + 1164);
    const auto *mi_1165 = buffer.data(mi + 1165);
    const auto *mi_1166 = buffer.data(mi + 1166);
    const auto *mi_1167 = buffer.data(mi + 1167);
    const auto *mi_1168 = buffer.data(mi + 1168);
    const auto *mi_1169 = buffer.data(mi + 1169);
    const auto *mi_1170 = buffer.data(mi + 1170);
    const auto *mi_1171 = buffer.data(mi + 1171);
    const auto *mi_1172 = buffer.data(mi + 1172);
    const auto *mi_1173 = buffer.data(mi + 1173);
    const auto *mi_1174 = buffer.data(mi + 1174);
    const auto *mi_1175 = buffer.data(mi + 1175);
    const auto *mi_1176 = buffer.data(mi + 1176);
    const auto *mi_1177 = buffer.data(mi + 1177);
    const auto *mi_1178 = buffer.data(mi + 1178);
    const auto *mi_1179 = buffer.data(mi + 1179);
    const auto *mi_1180 = buffer.data(mi + 1180);
    const auto *mi_1181 = buffer.data(mi + 1181);
    const auto *mi_1182 = buffer.data(mi + 1182);
    const auto *mi_1183 = buffer.data(mi + 1183);
    const auto *mi_1184 = buffer.data(mi + 1184);
    const auto *mi_1185 = buffer.data(mi + 1185);
    const auto *mi_1186 = buffer.data(mi + 1186);
    const auto *mi_1187 = buffer.data(mi + 1187);
    const auto *mi_1188 = buffer.data(mi + 1188);
    const auto *mi_1189 = buffer.data(mi + 1189);
    const auto *mi_1190 = buffer.data(mi + 1190);
    const auto *mi_1191 = buffer.data(mi + 1191);
    const auto *mi_1192 = buffer.data(mi + 1192);
    const auto *mi_1193 = buffer.data(mi + 1193);
    const auto *mi_1194 = buffer.data(mi + 1194);
    const auto *mi_1195 = buffer.data(mi + 1195);
    const auto *mi_1196 = buffer.data(mi + 1196);
    const auto *mi_1197 = buffer.data(mi + 1197);
    const auto *mi_1198 = buffer.data(mi + 1198);
    const auto *mi_1199 = buffer.data(mi + 1199);
    const auto *mi_1200 = buffer.data(mi + 1200);
    const auto *mi_1201 = buffer.data(mi + 1201);
    const auto *mi_1202 = buffer.data(mi + 1202);
    const auto *mi_1203 = buffer.data(mi + 1203);
    const auto *mi_1204 = buffer.data(mi + 1204);
    const auto *mi_1205 = buffer.data(mi + 1205);
    const auto *mi_1206 = buffer.data(mi + 1206);
    const auto *mi_1207 = buffer.data(mi + 1207);
    const auto *mi_1208 = buffer.data(mi + 1208);
    const auto *mi_1209 = buffer.data(mi + 1209);
    const auto *mi_1210 = buffer.data(mi + 1210);
    const auto *mi_1211 = buffer.data(mi + 1211);
    const auto *mi_1212 = buffer.data(mi + 1212);
    const auto *mi_1213 = buffer.data(mi + 1213);
    const auto *mi_1214 = buffer.data(mi + 1214);
    const auto *mi_1215 = buffer.data(mi + 1215);
    const auto *mi_1216 = buffer.data(mi + 1216);
    const auto *mi_1217 = buffer.data(mi + 1217);
    const auto *mi_1218 = buffer.data(mi + 1218);
    const auto *mi_1219 = buffer.data(mi + 1219);
    const auto *mi_1220 = buffer.data(mi + 1220);
    const auto *mi_1221 = buffer.data(mi + 1221);
    const auto *mi_1222 = buffer.data(mi + 1222);
    const auto *mi_1223 = buffer.data(mi + 1223);
    const auto *mi_1224 = buffer.data(mi + 1224);
    const auto *mi_1225 = buffer.data(mi + 1225);
    const auto *mi_1226 = buffer.data(mi + 1226);
    const auto *mi_1227 = buffer.data(mi + 1227);
    const auto *mi_1228 = buffer.data(mi + 1228);
    const auto *mi_1229 = buffer.data(mi + 1229);
    const auto *mi_1230 = buffer.data(mi + 1230);
    const auto *mi_1231 = buffer.data(mi + 1231);
    const auto *mi_1232 = buffer.data(mi + 1232);
    const auto *mi_1233 = buffer.data(mi + 1233);
    const auto *mi_1234 = buffer.data(mi + 1234);
    const auto *mi_1235 = buffer.data(mi + 1235);
    const auto *mi_1236 = buffer.data(mi + 1236);
    const auto *mi_1237 = buffer.data(mi + 1237);
    const auto *mi_1238 = buffer.data(mi + 1238);
    const auto *mi_1239 = buffer.data(mi + 1239);
    const auto *mi_1240 = buffer.data(mi + 1240);
    const auto *mi_1241 = buffer.data(mi + 1241);
    const auto *mi_1242 = buffer.data(mi + 1242);
    const auto *mi_1243 = buffer.data(mi + 1243);
    const auto *mi_1244 = buffer.data(mi + 1244);
    const auto *mi_1245 = buffer.data(mi + 1245);
    const auto *mi_1246 = buffer.data(mi + 1246);
    const auto *mi_1247 = buffer.data(mi + 1247);
    const auto *mi_1248 = buffer.data(mi + 1248);
    const auto *mi_1249 = buffer.data(mi + 1249);
    const auto *mi_1250 = buffer.data(mi + 1250);
    const auto *mi_1251 = buffer.data(mi + 1251);
    const auto *mi_1252 = buffer.data(mi + 1252);
    const auto *mi_1253 = buffer.data(mi + 1253);
    const auto *mi_1254 = buffer.data(mi + 1254);
    const auto *mi_1255 = buffer.data(mi + 1255);
    const auto *mi_1256 = buffer.data(mi + 1256);
    const auto *mi_1257 = buffer.data(mi + 1257);
    const auto *mi_1258 = buffer.data(mi + 1258);
    const auto *mi_1259 = buffer.data(mi + 1259);

#pragma omp simd aligned(t_1067, t_1068, t_1069, t_1070, t_1071, t_1072, t_1073, t_1074, \
                         mi_1067, mi_1068, mi_1069, mi_1070, mi_1071, mi_1072, mi_1073, \
                         mi_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_0 * mi_1067[k];

        t_1068[k] = f_0 * mi_1068[k];

        t_1069[k] = f_0 * mi_1069[k];

        t_1070[k] = f_0 * mi_1070[k];

        t_1071[k] = f_0 * mi_1071[k];

        t_1072[k] = f_0 * mi_1072[k];

        t_1073[k] = f_0 * mi_1073[k];

        t_1074[k] = f_0 * mi_1074[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, t_1078, t_1079, t_1080, t_1081, t_1082, \
                         mi_1075, mi_1076, mi_1077, mi_1078, mi_1079, mi_1080, mi_1081, \
                         mi_1082 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = f_0 * mi_1075[k];

        t_1076[k] = f_0 * mi_1076[k];

        t_1077[k] = f_0 * mi_1077[k];

        t_1078[k] = f_0 * mi_1078[k];

        t_1079[k] = f_0 * mi_1079[k];

        t_1080[k] = f_0 * mi_1080[k];

        t_1081[k] = f_0 * mi_1081[k];

        t_1082[k] = f_0 * mi_1082[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, t_1086, t_1087, t_1088, t_1089, t_1090, \
                         mi_1083, mi_1084, mi_1085, mi_1086, mi_1087, mi_1088, mi_1089, \
                         mi_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = f_0 * mi_1083[k];

        t_1084[k] = f_0 * mi_1084[k];

        t_1085[k] = f_0 * mi_1085[k];

        t_1086[k] = f_0 * mi_1086[k];

        t_1087[k] = f_0 * mi_1087[k];

        t_1088[k] = f_0 * mi_1088[k];

        t_1089[k] = f_0 * mi_1089[k];

        t_1090[k] = f_0 * mi_1090[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, t_1094, t_1095, t_1096, t_1097, t_1098, \
                         mi_1091, mi_1092, mi_1093, mi_1094, mi_1095, mi_1096, mi_1097, \
                         mi_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_0 * mi_1091[k];

        t_1092[k] = f_0 * mi_1092[k];

        t_1093[k] = f_0 * mi_1093[k];

        t_1094[k] = f_0 * mi_1094[k];

        t_1095[k] = f_0 * mi_1095[k];

        t_1096[k] = f_0 * mi_1096[k];

        t_1097[k] = f_0 * mi_1097[k];

        t_1098[k] = f_0 * mi_1098[k];
    }

#pragma omp simd aligned(t_1099, t_1100, t_1101, t_1102, t_1103, t_1104, t_1105, t_1106, \
                         mi_1099, mi_1100, mi_1101, mi_1102, mi_1103, mi_1104, mi_1105, \
                         mi_1106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1099[k] = f_0 * mi_1099[k];

        t_1100[k] = f_0 * mi_1100[k];

        t_1101[k] = f_0 * mi_1101[k];

        t_1102[k] = f_0 * mi_1102[k];

        t_1103[k] = f_0 * mi_1103[k];

        t_1104[k] = f_0 * mi_1104[k];

        t_1105[k] = f_0 * mi_1105[k];

        t_1106[k] = f_0 * mi_1106[k];
    }

#pragma omp simd aligned(t_1107, t_1108, t_1109, t_1110, t_1111, t_1112, t_1113, t_1114, \
                         mi_1107, mi_1108, mi_1109, mi_1110, mi_1111, mi_1112, mi_1113, \
                         mi_1114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1107[k] = f_0 * mi_1107[k];

        t_1108[k] = f_0 * mi_1108[k];

        t_1109[k] = f_0 * mi_1109[k];

        t_1110[k] = f_0 * mi_1110[k];

        t_1111[k] = f_0 * mi_1111[k];

        t_1112[k] = f_0 * mi_1112[k];

        t_1113[k] = f_0 * mi_1113[k];

        t_1114[k] = f_0 * mi_1114[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, t_1118, t_1119, t_1120, t_1121, t_1122, \
                         mi_1115, mi_1116, mi_1117, mi_1118, mi_1119, mi_1120, mi_1121, \
                         mi_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_0 * mi_1115[k];

        t_1116[k] = f_0 * mi_1116[k];

        t_1117[k] = f_0 * mi_1117[k];

        t_1118[k] = f_0 * mi_1118[k];

        t_1119[k] = f_0 * mi_1119[k];

        t_1120[k] = f_0 * mi_1120[k];

        t_1121[k] = f_0 * mi_1121[k];

        t_1122[k] = f_0 * mi_1122[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, t_1126, t_1127, t_1128, t_1129, t_1130, \
                         mi_1123, mi_1124, mi_1125, mi_1126, mi_1127, mi_1128, mi_1129, \
                         mi_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_0 * mi_1123[k];

        t_1124[k] = f_0 * mi_1124[k];

        t_1125[k] = f_0 * mi_1125[k];

        t_1126[k] = f_0 * mi_1126[k];

        t_1127[k] = f_0 * mi_1127[k];

        t_1128[k] = f_0 * mi_1128[k];

        t_1129[k] = f_0 * mi_1129[k];

        t_1130[k] = f_0 * mi_1130[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, t_1134, t_1135, t_1136, t_1137, t_1138, \
                         mi_1131, mi_1132, mi_1133, mi_1134, mi_1135, mi_1136, mi_1137, \
                         mi_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = f_0 * mi_1131[k];

        t_1132[k] = f_0 * mi_1132[k];

        t_1133[k] = f_0 * mi_1133[k];

        t_1134[k] = f_0 * mi_1134[k];

        t_1135[k] = f_0 * mi_1135[k];

        t_1136[k] = f_0 * mi_1136[k];

        t_1137[k] = f_0 * mi_1137[k];

        t_1138[k] = f_0 * mi_1138[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, t_1142, t_1143, t_1144, t_1145, t_1146, \
                         mi_1139, mi_1140, mi_1141, mi_1142, mi_1143, mi_1144, mi_1145, \
                         mi_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_0 * mi_1139[k];

        t_1140[k] = f_0 * mi_1140[k];

        t_1141[k] = f_0 * mi_1141[k];

        t_1142[k] = f_0 * mi_1142[k];

        t_1143[k] = f_0 * mi_1143[k];

        t_1144[k] = f_0 * mi_1144[k];

        t_1145[k] = f_0 * mi_1145[k];

        t_1146[k] = f_0 * mi_1146[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, t_1150, t_1151, t_1152, t_1153, t_1154, \
                         mi_1147, mi_1148, mi_1149, mi_1150, mi_1151, mi_1152, mi_1153, \
                         mi_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_0 * mi_1147[k];

        t_1148[k] = f_0 * mi_1148[k];

        t_1149[k] = f_0 * mi_1149[k];

        t_1150[k] = f_0 * mi_1150[k];

        t_1151[k] = f_0 * mi_1151[k];

        t_1152[k] = f_0 * mi_1152[k];

        t_1153[k] = f_0 * mi_1153[k];

        t_1154[k] = f_0 * mi_1154[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, t_1159, t_1160, t_1161, t_1162, \
                         mi_1155, mi_1156, mi_1157, mi_1158, mi_1159, mi_1160, mi_1161, \
                         mi_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_0 * mi_1155[k];

        t_1156[k] = f_0 * mi_1156[k];

        t_1157[k] = f_0 * mi_1157[k];

        t_1158[k] = f_0 * mi_1158[k];

        t_1159[k] = f_0 * mi_1159[k];

        t_1160[k] = f_0 * mi_1160[k];

        t_1161[k] = f_0 * mi_1161[k];

        t_1162[k] = f_0 * mi_1162[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, t_1166, t_1167, t_1168, t_1169, t_1170, \
                         mi_1163, mi_1164, mi_1165, mi_1166, mi_1167, mi_1168, mi_1169, \
                         mi_1170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_0 * mi_1163[k];

        t_1164[k] = f_0 * mi_1164[k];

        t_1165[k] = f_0 * mi_1165[k];

        t_1166[k] = f_0 * mi_1166[k];

        t_1167[k] = f_0 * mi_1167[k];

        t_1168[k] = f_0 * mi_1168[k];

        t_1169[k] = f_0 * mi_1169[k];

        t_1170[k] = f_0 * mi_1170[k];
    }

#pragma omp simd aligned(t_1171, t_1172, t_1173, t_1174, t_1175, t_1176, t_1177, t_1178, \
                         mi_1171, mi_1172, mi_1173, mi_1174, mi_1175, mi_1176, mi_1177, \
                         mi_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1171[k] = f_0 * mi_1171[k];

        t_1172[k] = f_0 * mi_1172[k];

        t_1173[k] = f_0 * mi_1173[k];

        t_1174[k] = f_0 * mi_1174[k];

        t_1175[k] = f_0 * mi_1175[k];

        t_1176[k] = f_0 * mi_1176[k];

        t_1177[k] = f_0 * mi_1177[k];

        t_1178[k] = f_0 * mi_1178[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, t_1182, t_1183, t_1184, t_1185, t_1186, \
                         mi_1179, mi_1180, mi_1181, mi_1182, mi_1183, mi_1184, mi_1185, \
                         mi_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = f_0 * mi_1179[k];

        t_1180[k] = f_0 * mi_1180[k];

        t_1181[k] = f_0 * mi_1181[k];

        t_1182[k] = f_0 * mi_1182[k];

        t_1183[k] = f_0 * mi_1183[k];

        t_1184[k] = f_0 * mi_1184[k];

        t_1185[k] = f_0 * mi_1185[k];

        t_1186[k] = f_0 * mi_1186[k];
    }

#pragma omp simd aligned(t_1187, t_1188, t_1189, t_1190, t_1191, t_1192, t_1193, t_1194, \
                         mi_1187, mi_1188, mi_1189, mi_1190, mi_1191, mi_1192, mi_1193, \
                         mi_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1187[k] = f_0 * mi_1187[k];

        t_1188[k] = f_0 * mi_1188[k];

        t_1189[k] = f_0 * mi_1189[k];

        t_1190[k] = f_0 * mi_1190[k];

        t_1191[k] = f_0 * mi_1191[k];

        t_1192[k] = f_0 * mi_1192[k];

        t_1193[k] = f_0 * mi_1193[k];

        t_1194[k] = f_0 * mi_1194[k];
    }

#pragma omp simd aligned(t_1195, t_1196, t_1197, t_1198, t_1199, t_1200, t_1201, t_1202, \
                         mi_1195, mi_1196, mi_1197, mi_1198, mi_1199, mi_1200, mi_1201, \
                         mi_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1195[k] = f_0 * mi_1195[k];

        t_1196[k] = f_0 * mi_1196[k];

        t_1197[k] = f_0 * mi_1197[k];

        t_1198[k] = f_0 * mi_1198[k];

        t_1199[k] = f_0 * mi_1199[k];

        t_1200[k] = f_0 * mi_1200[k];

        t_1201[k] = f_0 * mi_1201[k];

        t_1202[k] = f_0 * mi_1202[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, t_1206, t_1207, t_1208, t_1209, t_1210, \
                         mi_1203, mi_1204, mi_1205, mi_1206, mi_1207, mi_1208, mi_1209, \
                         mi_1210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = f_0 * mi_1203[k];

        t_1204[k] = f_0 * mi_1204[k];

        t_1205[k] = f_0 * mi_1205[k];

        t_1206[k] = f_0 * mi_1206[k];

        t_1207[k] = f_0 * mi_1207[k];

        t_1208[k] = f_0 * mi_1208[k];

        t_1209[k] = f_0 * mi_1209[k];

        t_1210[k] = f_0 * mi_1210[k];
    }

#pragma omp simd aligned(t_1211, t_1212, t_1213, t_1214, t_1215, t_1216, t_1217, t_1218, \
                         mi_1211, mi_1212, mi_1213, mi_1214, mi_1215, mi_1216, mi_1217, \
                         mi_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1211[k] = f_0 * mi_1211[k];

        t_1212[k] = f_0 * mi_1212[k];

        t_1213[k] = f_0 * mi_1213[k];

        t_1214[k] = f_0 * mi_1214[k];

        t_1215[k] = f_0 * mi_1215[k];

        t_1216[k] = f_0 * mi_1216[k];

        t_1217[k] = f_0 * mi_1217[k];

        t_1218[k] = f_0 * mi_1218[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, t_1222, t_1223, t_1224, t_1225, t_1226, \
                         mi_1219, mi_1220, mi_1221, mi_1222, mi_1223, mi_1224, mi_1225, \
                         mi_1226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_0 * mi_1219[k];

        t_1220[k] = f_0 * mi_1220[k];

        t_1221[k] = f_0 * mi_1221[k];

        t_1222[k] = f_0 * mi_1222[k];

        t_1223[k] = f_0 * mi_1223[k];

        t_1224[k] = f_0 * mi_1224[k];

        t_1225[k] = f_0 * mi_1225[k];

        t_1226[k] = f_0 * mi_1226[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, t_1230, t_1231, t_1232, t_1233, t_1234, \
                         mi_1227, mi_1228, mi_1229, mi_1230, mi_1231, mi_1232, mi_1233, \
                         mi_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = f_0 * mi_1227[k];

        t_1228[k] = f_0 * mi_1228[k];

        t_1229[k] = f_0 * mi_1229[k];

        t_1230[k] = f_0 * mi_1230[k];

        t_1231[k] = f_0 * mi_1231[k];

        t_1232[k] = f_0 * mi_1232[k];

        t_1233[k] = f_0 * mi_1233[k];

        t_1234[k] = f_0 * mi_1234[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, t_1238, t_1239, t_1240, t_1241, t_1242, \
                         mi_1235, mi_1236, mi_1237, mi_1238, mi_1239, mi_1240, mi_1241, \
                         mi_1242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_0 * mi_1235[k];

        t_1236[k] = f_0 * mi_1236[k];

        t_1237[k] = f_0 * mi_1237[k];

        t_1238[k] = f_0 * mi_1238[k];

        t_1239[k] = f_0 * mi_1239[k];

        t_1240[k] = f_0 * mi_1240[k];

        t_1241[k] = f_0 * mi_1241[k];

        t_1242[k] = f_0 * mi_1242[k];
    }

#pragma omp simd aligned(t_1243, t_1244, t_1245, t_1246, t_1247, t_1248, t_1249, t_1250, \
                         mi_1243, mi_1244, mi_1245, mi_1246, mi_1247, mi_1248, mi_1249, \
                         mi_1250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1243[k] = f_0 * mi_1243[k];

        t_1244[k] = f_0 * mi_1244[k];

        t_1245[k] = f_0 * mi_1245[k];

        t_1246[k] = f_0 * mi_1246[k];

        t_1247[k] = f_0 * mi_1247[k];

        t_1248[k] = f_0 * mi_1248[k];

        t_1249[k] = f_0 * mi_1249[k];

        t_1250[k] = f_0 * mi_1250[k];
    }

#pragma omp simd aligned(t_1251, t_1252, t_1253, t_1254, t_1255, t_1256, t_1257, t_1258, \
                         mi_1251, mi_1252, mi_1253, mi_1254, mi_1255, mi_1256, mi_1257, \
                         mi_1258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1251[k] = f_0 * mi_1251[k];

        t_1252[k] = f_0 * mi_1252[k];

        t_1253[k] = f_0 * mi_1253[k];

        t_1254[k] = f_0 * mi_1254[k];

        t_1255[k] = f_0 * mi_1255[k];

        t_1256[k] = f_0 * mi_1256[k];

        t_1257[k] = f_0 * mi_1257[k];

        t_1258[k] = f_0 * mi_1258[k];
    }

#pragma omp simd aligned(t_1259, mi_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1259[k] = f_0 * mi_1259[k];
    }
}

auto
compute_prim_geom_10_li_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t ki, const size_t mi,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_li_electron_repulsion_0_piece0(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_0_piece1(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_0_piece2(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_0_piece3(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_0_piece4(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_0_piece5(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_0_piece6(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_0_piece7(buffer, target, mi, ncols, alpha);
}

static auto
compute_prim_geom_10_li_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
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

    const auto *ki_0 = buffer.data(ki + 0);
    const auto *ki_1 = buffer.data(ki + 1);
    const auto *ki_2 = buffer.data(ki + 2);
    const auto *ki_3 = buffer.data(ki + 3);
    const auto *ki_4 = buffer.data(ki + 4);
    const auto *ki_5 = buffer.data(ki + 5);
    const auto *ki_6 = buffer.data(ki + 6);
    const auto *ki_7 = buffer.data(ki + 7);
    const auto *ki_8 = buffer.data(ki + 8);
    const auto *ki_9 = buffer.data(ki + 9);
    const auto *ki_10 = buffer.data(ki + 10);
    const auto *ki_11 = buffer.data(ki + 11);
    const auto *ki_12 = buffer.data(ki + 12);
    const auto *ki_13 = buffer.data(ki + 13);
    const auto *ki_14 = buffer.data(ki + 14);
    const auto *ki_15 = buffer.data(ki + 15);
    const auto *ki_16 = buffer.data(ki + 16);
    const auto *ki_17 = buffer.data(ki + 17);
    const auto *ki_18 = buffer.data(ki + 18);
    const auto *ki_19 = buffer.data(ki + 19);
    const auto *ki_20 = buffer.data(ki + 20);
    const auto *ki_21 = buffer.data(ki + 21);
    const auto *ki_22 = buffer.data(ki + 22);
    const auto *ki_23 = buffer.data(ki + 23);
    const auto *ki_24 = buffer.data(ki + 24);
    const auto *ki_25 = buffer.data(ki + 25);
    const auto *ki_26 = buffer.data(ki + 26);
    const auto *ki_27 = buffer.data(ki + 27);
    const auto *ki_28 = buffer.data(ki + 28);
    const auto *ki_29 = buffer.data(ki + 29);
    const auto *ki_30 = buffer.data(ki + 30);
    const auto *ki_31 = buffer.data(ki + 31);
    const auto *ki_32 = buffer.data(ki + 32);
    const auto *ki_33 = buffer.data(ki + 33);
    const auto *ki_34 = buffer.data(ki + 34);
    const auto *ki_35 = buffer.data(ki + 35);
    const auto *ki_36 = buffer.data(ki + 36);
    const auto *ki_37 = buffer.data(ki + 37);
    const auto *ki_38 = buffer.data(ki + 38);
    const auto *ki_39 = buffer.data(ki + 39);
    const auto *ki_40 = buffer.data(ki + 40);
    const auto *ki_41 = buffer.data(ki + 41);
    const auto *ki_42 = buffer.data(ki + 42);
    const auto *ki_43 = buffer.data(ki + 43);
    const auto *ki_44 = buffer.data(ki + 44);
    const auto *ki_45 = buffer.data(ki + 45);
    const auto *ki_46 = buffer.data(ki + 46);
    const auto *ki_47 = buffer.data(ki + 47);
    const auto *ki_48 = buffer.data(ki + 48);
    const auto *ki_49 = buffer.data(ki + 49);
    const auto *ki_50 = buffer.data(ki + 50);
    const auto *ki_51 = buffer.data(ki + 51);
    const auto *ki_52 = buffer.data(ki + 52);
    const auto *ki_53 = buffer.data(ki + 53);
    const auto *ki_54 = buffer.data(ki + 54);
    const auto *ki_55 = buffer.data(ki + 55);
    const auto *ki_56 = buffer.data(ki + 56);
    const auto *ki_57 = buffer.data(ki + 57);
    const auto *ki_58 = buffer.data(ki + 58);
    const auto *ki_59 = buffer.data(ki + 59);
    const auto *ki_60 = buffer.data(ki + 60);
    const auto *ki_61 = buffer.data(ki + 61);
    const auto *ki_62 = buffer.data(ki + 62);
    const auto *ki_63 = buffer.data(ki + 63);
    const auto *ki_64 = buffer.data(ki + 64);
    const auto *ki_65 = buffer.data(ki + 65);
    const auto *ki_66 = buffer.data(ki + 66);
    const auto *ki_67 = buffer.data(ki + 67);
    const auto *ki_68 = buffer.data(ki + 68);
    const auto *ki_69 = buffer.data(ki + 69);
    const auto *ki_70 = buffer.data(ki + 70);
    const auto *ki_71 = buffer.data(ki + 71);
    const auto *ki_72 = buffer.data(ki + 72);
    const auto *ki_73 = buffer.data(ki + 73);
    const auto *ki_74 = buffer.data(ki + 74);
    const auto *ki_75 = buffer.data(ki + 75);
    const auto *ki_76 = buffer.data(ki + 76);
    const auto *ki_77 = buffer.data(ki + 77);
    const auto *ki_78 = buffer.data(ki + 78);
    const auto *ki_79 = buffer.data(ki + 79);
    const auto *ki_80 = buffer.data(ki + 80);
    const auto *ki_81 = buffer.data(ki + 81);
    const auto *ki_82 = buffer.data(ki + 82);
    const auto *ki_83 = buffer.data(ki + 83);
    const auto *ki_84 = buffer.data(ki + 84);
    const auto *ki_85 = buffer.data(ki + 85);
    const auto *ki_86 = buffer.data(ki + 86);
    const auto *ki_87 = buffer.data(ki + 87);
    const auto *ki_88 = buffer.data(ki + 88);
    const auto *ki_89 = buffer.data(ki + 89);
    const auto *ki_90 = buffer.data(ki + 90);

    const auto *mi_28 = buffer.data(mi + 28);
    const auto *mi_29 = buffer.data(mi + 29);
    const auto *mi_30 = buffer.data(mi + 30);
    const auto *mi_31 = buffer.data(mi + 31);
    const auto *mi_32 = buffer.data(mi + 32);
    const auto *mi_33 = buffer.data(mi + 33);
    const auto *mi_34 = buffer.data(mi + 34);
    const auto *mi_35 = buffer.data(mi + 35);
    const auto *mi_36 = buffer.data(mi + 36);
    const auto *mi_37 = buffer.data(mi + 37);
    const auto *mi_38 = buffer.data(mi + 38);
    const auto *mi_39 = buffer.data(mi + 39);
    const auto *mi_40 = buffer.data(mi + 40);
    const auto *mi_41 = buffer.data(mi + 41);
    const auto *mi_42 = buffer.data(mi + 42);
    const auto *mi_43 = buffer.data(mi + 43);
    const auto *mi_44 = buffer.data(mi + 44);
    const auto *mi_45 = buffer.data(mi + 45);
    const auto *mi_46 = buffer.data(mi + 46);
    const auto *mi_47 = buffer.data(mi + 47);
    const auto *mi_48 = buffer.data(mi + 48);
    const auto *mi_49 = buffer.data(mi + 49);
    const auto *mi_50 = buffer.data(mi + 50);
    const auto *mi_51 = buffer.data(mi + 51);
    const auto *mi_52 = buffer.data(mi + 52);
    const auto *mi_53 = buffer.data(mi + 53);
    const auto *mi_54 = buffer.data(mi + 54);
    const auto *mi_55 = buffer.data(mi + 55);
    const auto *mi_84 = buffer.data(mi + 84);
    const auto *mi_85 = buffer.data(mi + 85);
    const auto *mi_86 = buffer.data(mi + 86);
    const auto *mi_87 = buffer.data(mi + 87);
    const auto *mi_88 = buffer.data(mi + 88);
    const auto *mi_89 = buffer.data(mi + 89);
    const auto *mi_90 = buffer.data(mi + 90);
    const auto *mi_91 = buffer.data(mi + 91);
    const auto *mi_92 = buffer.data(mi + 92);
    const auto *mi_93 = buffer.data(mi + 93);
    const auto *mi_94 = buffer.data(mi + 94);
    const auto *mi_95 = buffer.data(mi + 95);
    const auto *mi_96 = buffer.data(mi + 96);
    const auto *mi_97 = buffer.data(mi + 97);
    const auto *mi_98 = buffer.data(mi + 98);
    const auto *mi_99 = buffer.data(mi + 99);
    const auto *mi_100 = buffer.data(mi + 100);
    const auto *mi_101 = buffer.data(mi + 101);
    const auto *mi_102 = buffer.data(mi + 102);
    const auto *mi_103 = buffer.data(mi + 103);
    const auto *mi_104 = buffer.data(mi + 104);
    const auto *mi_105 = buffer.data(mi + 105);
    const auto *mi_106 = buffer.data(mi + 106);
    const auto *mi_107 = buffer.data(mi + 107);
    const auto *mi_108 = buffer.data(mi + 108);
    const auto *mi_109 = buffer.data(mi + 109);
    const auto *mi_110 = buffer.data(mi + 110);
    const auto *mi_111 = buffer.data(mi + 111);
    const auto *mi_112 = buffer.data(mi + 112);
    const auto *mi_113 = buffer.data(mi + 113);
    const auto *mi_114 = buffer.data(mi + 114);
    const auto *mi_115 = buffer.data(mi + 115);
    const auto *mi_116 = buffer.data(mi + 116);
    const auto *mi_117 = buffer.data(mi + 117);
    const auto *mi_118 = buffer.data(mi + 118);
    const auto *mi_119 = buffer.data(mi + 119);
    const auto *mi_120 = buffer.data(mi + 120);
    const auto *mi_121 = buffer.data(mi + 121);
    const auto *mi_122 = buffer.data(mi + 122);
    const auto *mi_123 = buffer.data(mi + 123);
    const auto *mi_124 = buffer.data(mi + 124);
    const auto *mi_125 = buffer.data(mi + 125);
    const auto *mi_126 = buffer.data(mi + 126);
    const auto *mi_127 = buffer.data(mi + 127);
    const auto *mi_128 = buffer.data(mi + 128);
    const auto *mi_129 = buffer.data(mi + 129);
    const auto *mi_130 = buffer.data(mi + 130);
    const auto *mi_131 = buffer.data(mi + 131);
    const auto *mi_132 = buffer.data(mi + 132);
    const auto *mi_133 = buffer.data(mi + 133);
    const auto *mi_134 = buffer.data(mi + 134);
    const auto *mi_135 = buffer.data(mi + 135);
    const auto *mi_136 = buffer.data(mi + 136);
    const auto *mi_137 = buffer.data(mi + 137);
    const auto *mi_138 = buffer.data(mi + 138);
    const auto *mi_139 = buffer.data(mi + 139);
    const auto *mi_168 = buffer.data(mi + 168);
    const auto *mi_169 = buffer.data(mi + 169);
    const auto *mi_170 = buffer.data(mi + 170);
    const auto *mi_171 = buffer.data(mi + 171);
    const auto *mi_172 = buffer.data(mi + 172);
    const auto *mi_173 = buffer.data(mi + 173);
    const auto *mi_174 = buffer.data(mi + 174);
    const auto *mi_175 = buffer.data(mi + 175);
    const auto *mi_176 = buffer.data(mi + 176);
    const auto *mi_177 = buffer.data(mi + 177);
    const auto *mi_178 = buffer.data(mi + 178);
    const auto *mi_179 = buffer.data(mi + 179);
    const auto *mi_180 = buffer.data(mi + 180);
    const auto *mi_181 = buffer.data(mi + 181);
    const auto *mi_182 = buffer.data(mi + 182);
    const auto *mi_183 = buffer.data(mi + 183);
    const auto *mi_184 = buffer.data(mi + 184);
    const auto *mi_185 = buffer.data(mi + 185);
    const auto *mi_186 = buffer.data(mi + 186);
    const auto *mi_187 = buffer.data(mi + 187);
    const auto *mi_188 = buffer.data(mi + 188);
    const auto *mi_189 = buffer.data(mi + 189);
    const auto *mi_190 = buffer.data(mi + 190);
    const auto *mi_191 = buffer.data(mi + 191);
    const auto *mi_192 = buffer.data(mi + 192);
    const auto *mi_193 = buffer.data(mi + 193);
    const auto *mi_194 = buffer.data(mi + 194);
    const auto *mi_195 = buffer.data(mi + 195);
    const auto *mi_196 = buffer.data(mi + 196);
    const auto *mi_197 = buffer.data(mi + 197);
    const auto *mi_198 = buffer.data(mi + 198);
    const auto *mi_199 = buffer.data(mi + 199);
    const auto *mi_200 = buffer.data(mi + 200);
    const auto *mi_201 = buffer.data(mi + 201);
    const auto *mi_202 = buffer.data(mi + 202);
    const auto *mi_203 = buffer.data(mi + 203);
    const auto *mi_204 = buffer.data(mi + 204);
    const auto *mi_205 = buffer.data(mi + 205);
    const auto *mi_206 = buffer.data(mi + 206);
    const auto *mi_207 = buffer.data(mi + 207);
    const auto *mi_208 = buffer.data(mi + 208);
    const auto *mi_209 = buffer.data(mi + 209);
    const auto *mi_210 = buffer.data(mi + 210);
    const auto *mi_211 = buffer.data(mi + 211);
    const auto *mi_212 = buffer.data(mi + 212);
    const auto *mi_213 = buffer.data(mi + 213);
    const auto *mi_214 = buffer.data(mi + 214);
    const auto *mi_215 = buffer.data(mi + 215);
    const auto *mi_216 = buffer.data(mi + 216);
    const auto *mi_217 = buffer.data(mi + 217);
    const auto *mi_218 = buffer.data(mi + 218);
    const auto *mi_219 = buffer.data(mi + 219);
    const auto *mi_220 = buffer.data(mi + 220);
    const auto *mi_221 = buffer.data(mi + 221);
    const auto *mi_222 = buffer.data(mi + 222);
    const auto *mi_223 = buffer.data(mi + 223);
    const auto *mi_224 = buffer.data(mi + 224);
    const auto *mi_225 = buffer.data(mi + 225);
    const auto *mi_226 = buffer.data(mi + 226);
    const auto *mi_227 = buffer.data(mi + 227);
    const auto *mi_228 = buffer.data(mi + 228);
    const auto *mi_229 = buffer.data(mi + 229);
    const auto *mi_230 = buffer.data(mi + 230);
    const auto *mi_231 = buffer.data(mi + 231);
    const auto *mi_232 = buffer.data(mi + 232);
    const auto *mi_233 = buffer.data(mi + 233);
    const auto *mi_234 = buffer.data(mi + 234);
    const auto *mi_235 = buffer.data(mi + 235);
    const auto *mi_236 = buffer.data(mi + 236);
    const auto *mi_237 = buffer.data(mi + 237);
    const auto *mi_238 = buffer.data(mi + 238);
    const auto *mi_239 = buffer.data(mi + 239);
    const auto *mi_240 = buffer.data(mi + 240);
    const auto *mi_241 = buffer.data(mi + 241);
    const auto *mi_242 = buffer.data(mi + 242);
    const auto *mi_243 = buffer.data(mi + 243);
    const auto *mi_244 = buffer.data(mi + 244);
    const auto *mi_245 = buffer.data(mi + 245);
    const auto *mi_246 = buffer.data(mi + 246);
    const auto *mi_247 = buffer.data(mi + 247);
    const auto *mi_248 = buffer.data(mi + 248);
    const auto *mi_249 = buffer.data(mi + 249);
    const auto *mi_250 = buffer.data(mi + 250);
    const auto *mi_251 = buffer.data(mi + 251);
    const auto *mi_280 = buffer.data(mi + 280);
    const auto *mi_281 = buffer.data(mi + 281);
    const auto *mi_282 = buffer.data(mi + 282);
    const auto *mi_283 = buffer.data(mi + 283);
    const auto *mi_284 = buffer.data(mi + 284);
    const auto *mi_285 = buffer.data(mi + 285);
    const auto *mi_286 = buffer.data(mi + 286);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, mi_28, mi_29, mi_30, mi_31, \
                         mi_32, mi_33, mi_34, mi_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * mi_28[k];

        t_1[k] = f_0 * mi_29[k];

        t_2[k] = f_0 * mi_30[k];

        t_3[k] = f_0 * mi_31[k];

        t_4[k] = f_0 * mi_32[k];

        t_5[k] = f_0 * mi_33[k];

        t_6[k] = f_0 * mi_34[k];

        t_7[k] = f_0 * mi_35[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, mi_36, mi_37, mi_38, \
                         mi_39, mi_40, mi_41, mi_42, mi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * mi_36[k];

        t_9[k] = f_0 * mi_37[k];

        t_10[k] = f_0 * mi_38[k];

        t_11[k] = f_0 * mi_39[k];

        t_12[k] = f_0 * mi_40[k];

        t_13[k] = f_0 * mi_41[k];

        t_14[k] = f_0 * mi_42[k];

        t_15[k] = f_0 * mi_43[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, mi_44, mi_45, mi_46, \
                         mi_47, mi_48, mi_49, mi_50, mi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * mi_44[k];

        t_17[k] = f_0 * mi_45[k];

        t_18[k] = f_0 * mi_46[k];

        t_19[k] = f_0 * mi_47[k];

        t_20[k] = f_0 * mi_48[k];

        t_21[k] = f_0 * mi_49[k];

        t_22[k] = f_0 * mi_50[k];

        t_23[k] = f_0 * mi_51[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, ki_0, ki_1, mi_52, mi_53, mi_54, \
                         mi_55, mi_84, mi_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * mi_52[k];

        t_25[k] = f_0 * mi_53[k];

        t_26[k] = f_0 * mi_54[k];

        t_27[k] = f_0 * mi_55[k];

        t_28[k] = -ki_0[k]
                  + f_0 * mi_84[k];

        t_29[k] = -ki_1[k]
                  + f_0 * mi_85[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ki_2, ki_3, ki_4, ki_5, ki_6, mi_86, \
                         mi_87, mi_88, mi_89, mi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -ki_2[k]
                  + f_0 * mi_86[k];

        t_31[k] = -ki_3[k]
                  + f_0 * mi_87[k];

        t_32[k] = -ki_4[k]
                  + f_0 * mi_88[k];

        t_33[k] = -ki_5[k]
                  + f_0 * mi_89[k];

        t_34[k] = -ki_6[k]
                  + f_0 * mi_90[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ki_7, ki_8, ki_9, ki_10, ki_11, mi_91, \
                         mi_92, mi_93, mi_94, mi_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -ki_7[k]
                  + f_0 * mi_91[k];

        t_36[k] = -ki_8[k]
                  + f_0 * mi_92[k];

        t_37[k] = -ki_9[k]
                  + f_0 * mi_93[k];

        t_38[k] = -ki_10[k]
                  + f_0 * mi_94[k];

        t_39[k] = -ki_11[k]
                  + f_0 * mi_95[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ki_12, ki_13, ki_14, ki_15, ki_16, \
                         mi_96, mi_97, mi_98, mi_99, mi_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -ki_12[k]
                  + f_0 * mi_96[k];

        t_41[k] = -ki_13[k]
                  + f_0 * mi_97[k];

        t_42[k] = -ki_14[k]
                  + f_0 * mi_98[k];

        t_43[k] = -ki_15[k]
                  + f_0 * mi_99[k];

        t_44[k] = -ki_16[k]
                  + f_0 * mi_100[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ki_17, ki_18, ki_19, ki_20, ki_21, \
                         mi_101, mi_102, mi_103, mi_104, mi_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -ki_17[k]
                  + f_0 * mi_101[k];

        t_46[k] = -ki_18[k]
                  + f_0 * mi_102[k];

        t_47[k] = -ki_19[k]
                  + f_0 * mi_103[k];

        t_48[k] = -ki_20[k]
                  + f_0 * mi_104[k];

        t_49[k] = -ki_21[k]
                  + f_0 * mi_105[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ki_22, ki_23, ki_24, ki_25, ki_26, \
                         mi_106, mi_107, mi_108, mi_109, mi_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -ki_22[k]
                  + f_0 * mi_106[k];

        t_51[k] = -ki_23[k]
                  + f_0 * mi_107[k];

        t_52[k] = -ki_24[k]
                  + f_0 * mi_108[k];

        t_53[k] = -ki_25[k]
                  + f_0 * mi_109[k];

        t_54[k] = -ki_26[k]
                  + f_0 * mi_110[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, ki_27, mi_111, mi_112, \
                         mi_113, mi_114, mi_115, mi_116, mi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -ki_27[k]
                  + f_0 * mi_111[k];

        t_56[k] = f_0 * mi_112[k];

        t_57[k] = f_0 * mi_113[k];

        t_58[k] = f_0 * mi_114[k];

        t_59[k] = f_0 * mi_115[k];

        t_60[k] = f_0 * mi_116[k];

        t_61[k] = f_0 * mi_117[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, t_67, t_68, t_69, mi_118, mi_119, \
                         mi_120, mi_121, mi_122, mi_123, mi_124, \
                         mi_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * mi_118[k];

        t_63[k] = f_0 * mi_119[k];

        t_64[k] = f_0 * mi_120[k];

        t_65[k] = f_0 * mi_121[k];

        t_66[k] = f_0 * mi_122[k];

        t_67[k] = f_0 * mi_123[k];

        t_68[k] = f_0 * mi_124[k];

        t_69[k] = f_0 * mi_125[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, t_75, t_76, t_77, mi_126, mi_127, \
                         mi_128, mi_129, mi_130, mi_131, mi_132, \
                         mi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_0 * mi_126[k];

        t_71[k] = f_0 * mi_127[k];

        t_72[k] = f_0 * mi_128[k];

        t_73[k] = f_0 * mi_129[k];

        t_74[k] = f_0 * mi_130[k];

        t_75[k] = f_0 * mi_131[k];

        t_76[k] = f_0 * mi_132[k];

        t_77[k] = f_0 * mi_133[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, t_84, ki_28, mi_134, mi_135, \
                         mi_136, mi_137, mi_138, mi_139, mi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * mi_134[k];

        t_79[k] = f_0 * mi_135[k];

        t_80[k] = f_0 * mi_136[k];

        t_81[k] = f_0 * mi_137[k];

        t_82[k] = f_0 * mi_138[k];

        t_83[k] = f_0 * mi_139[k];

        t_84[k] = -2.0 * ki_28[k]
                  + f_0 * mi_168[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ki_29, ki_30, ki_31, ki_32, ki_33, \
                         mi_169, mi_170, mi_171, mi_172, mi_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * ki_29[k]
                  + f_0 * mi_169[k];

        t_86[k] = -2.0 * ki_30[k]
                  + f_0 * mi_170[k];

        t_87[k] = -2.0 * ki_31[k]
                  + f_0 * mi_171[k];

        t_88[k] = -2.0 * ki_32[k]
                  + f_0 * mi_172[k];

        t_89[k] = -2.0 * ki_33[k]
                  + f_0 * mi_173[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ki_34, ki_35, ki_36, ki_37, ki_38, \
                         mi_174, mi_175, mi_176, mi_177, mi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * ki_34[k]
                  + f_0 * mi_174[k];

        t_91[k] = -2.0 * ki_35[k]
                  + f_0 * mi_175[k];

        t_92[k] = -2.0 * ki_36[k]
                  + f_0 * mi_176[k];

        t_93[k] = -2.0 * ki_37[k]
                  + f_0 * mi_177[k];

        t_94[k] = -2.0 * ki_38[k]
                  + f_0 * mi_178[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ki_39, ki_40, ki_41, ki_42, ki_43, \
                         mi_179, mi_180, mi_181, mi_182, mi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * ki_39[k]
                  + f_0 * mi_179[k];

        t_96[k] = -2.0 * ki_40[k]
                  + f_0 * mi_180[k];

        t_97[k] = -2.0 * ki_41[k]
                  + f_0 * mi_181[k];

        t_98[k] = -2.0 * ki_42[k]
                  + f_0 * mi_182[k];

        t_99[k] = -2.0 * ki_43[k]
                  + f_0 * mi_183[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ki_44, ki_45, ki_46, ki_47, ki_48, \
                         mi_184, mi_185, mi_186, mi_187, mi_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * ki_44[k]
                   + f_0 * mi_184[k];

        t_101[k] = -2.0 * ki_45[k]
                   + f_0 * mi_185[k];

        t_102[k] = -2.0 * ki_46[k]
                   + f_0 * mi_186[k];

        t_103[k] = -2.0 * ki_47[k]
                   + f_0 * mi_187[k];

        t_104[k] = -2.0 * ki_48[k]
                   + f_0 * mi_188[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ki_49, ki_50, ki_51, ki_52, ki_53, \
                         mi_189, mi_190, mi_191, mi_192, mi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * ki_49[k]
                   + f_0 * mi_189[k];

        t_106[k] = -2.0 * ki_50[k]
                   + f_0 * mi_190[k];

        t_107[k] = -2.0 * ki_51[k]
                   + f_0 * mi_191[k];

        t_108[k] = -2.0 * ki_52[k]
                   + f_0 * mi_192[k];

        t_109[k] = -2.0 * ki_53[k]
                   + f_0 * mi_193[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ki_54, ki_55, ki_56, ki_57, ki_58, \
                         mi_194, mi_195, mi_196, mi_197, mi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * ki_54[k]
                   + f_0 * mi_194[k];

        t_111[k] = -2.0 * ki_55[k]
                   + f_0 * mi_195[k];

        t_112[k] = -ki_56[k]
                   + f_0 * mi_196[k];

        t_113[k] = -ki_57[k]
                   + f_0 * mi_197[k];

        t_114[k] = -ki_58[k]
                   + f_0 * mi_198[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ki_59, ki_60, ki_61, ki_62, ki_63, \
                         mi_199, mi_200, mi_201, mi_202, mi_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -ki_59[k]
                   + f_0 * mi_199[k];

        t_116[k] = -ki_60[k]
                   + f_0 * mi_200[k];

        t_117[k] = -ki_61[k]
                   + f_0 * mi_201[k];

        t_118[k] = -ki_62[k]
                   + f_0 * mi_202[k];

        t_119[k] = -ki_63[k]
                   + f_0 * mi_203[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ki_64, ki_65, ki_66, ki_67, ki_68, \
                         mi_204, mi_205, mi_206, mi_207, mi_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -ki_64[k]
                   + f_0 * mi_204[k];

        t_121[k] = -ki_65[k]
                   + f_0 * mi_205[k];

        t_122[k] = -ki_66[k]
                   + f_0 * mi_206[k];

        t_123[k] = -ki_67[k]
                   + f_0 * mi_207[k];

        t_124[k] = -ki_68[k]
                   + f_0 * mi_208[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ki_69, ki_70, ki_71, ki_72, ki_73, \
                         mi_209, mi_210, mi_211, mi_212, mi_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -ki_69[k]
                   + f_0 * mi_209[k];

        t_126[k] = -ki_70[k]
                   + f_0 * mi_210[k];

        t_127[k] = -ki_71[k]
                   + f_0 * mi_211[k];

        t_128[k] = -ki_72[k]
                   + f_0 * mi_212[k];

        t_129[k] = -ki_73[k]
                   + f_0 * mi_213[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ki_74, ki_75, ki_76, ki_77, ki_78, \
                         mi_214, mi_215, mi_216, mi_217, mi_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -ki_74[k]
                   + f_0 * mi_214[k];

        t_131[k] = -ki_75[k]
                   + f_0 * mi_215[k];

        t_132[k] = -ki_76[k]
                   + f_0 * mi_216[k];

        t_133[k] = -ki_77[k]
                   + f_0 * mi_217[k];

        t_134[k] = -ki_78[k]
                   + f_0 * mi_218[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ki_79, ki_80, ki_81, ki_82, ki_83, \
                         mi_219, mi_220, mi_221, mi_222, mi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -ki_79[k]
                   + f_0 * mi_219[k];

        t_136[k] = -ki_80[k]
                   + f_0 * mi_220[k];

        t_137[k] = -ki_81[k]
                   + f_0 * mi_221[k];

        t_138[k] = -ki_82[k]
                   + f_0 * mi_222[k];

        t_139[k] = -ki_83[k]
                   + f_0 * mi_223[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, t_146, t_147, mi_224, \
                         mi_225, mi_226, mi_227, mi_228, mi_229, mi_230, \
                         mi_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * mi_224[k];

        t_141[k] = f_0 * mi_225[k];

        t_142[k] = f_0 * mi_226[k];

        t_143[k] = f_0 * mi_227[k];

        t_144[k] = f_0 * mi_228[k];

        t_145[k] = f_0 * mi_229[k];

        t_146[k] = f_0 * mi_230[k];

        t_147[k] = f_0 * mi_231[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, t_153, t_154, t_155, mi_232, \
                         mi_233, mi_234, mi_235, mi_236, mi_237, mi_238, \
                         mi_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * mi_232[k];

        t_149[k] = f_0 * mi_233[k];

        t_150[k] = f_0 * mi_234[k];

        t_151[k] = f_0 * mi_235[k];

        t_152[k] = f_0 * mi_236[k];

        t_153[k] = f_0 * mi_237[k];

        t_154[k] = f_0 * mi_238[k];

        t_155[k] = f_0 * mi_239[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, t_162, t_163, mi_240, \
                         mi_241, mi_242, mi_243, mi_244, mi_245, mi_246, \
                         mi_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * mi_240[k];

        t_157[k] = f_0 * mi_241[k];

        t_158[k] = f_0 * mi_242[k];

        t_159[k] = f_0 * mi_243[k];

        t_160[k] = f_0 * mi_244[k];

        t_161[k] = f_0 * mi_245[k];

        t_162[k] = f_0 * mi_246[k];

        t_163[k] = f_0 * mi_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, t_169, ki_84, ki_85, mi_248, \
                         mi_249, mi_250, mi_251, mi_280, mi_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_0 * mi_248[k];

        t_165[k] = f_0 * mi_249[k];

        t_166[k] = f_0 * mi_250[k];

        t_167[k] = f_0 * mi_251[k];

        t_168[k] = -3.0 * ki_84[k]
                   + f_0 * mi_280[k];

        t_169[k] = -3.0 * ki_85[k]
                   + f_0 * mi_281[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ki_86, ki_87, ki_88, ki_89, ki_90, \
                         mi_282, mi_283, mi_284, mi_285, mi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * ki_86[k]
                   + f_0 * mi_282[k];

        t_171[k] = -3.0 * ki_87[k]
                   + f_0 * mi_283[k];

        t_172[k] = -3.0 * ki_88[k]
                   + f_0 * mi_284[k];

        t_173[k] = -3.0 * ki_89[k]
                   + f_0 * mi_285[k];

        t_174[k] = -3.0 * ki_90[k]
                   + f_0 * mi_286[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);
    auto *t_334 = buffer.data(target + 334);

    const auto *ki_91 = buffer.data(ki + 91);
    const auto *ki_92 = buffer.data(ki + 92);
    const auto *ki_93 = buffer.data(ki + 93);
    const auto *ki_94 = buffer.data(ki + 94);
    const auto *ki_95 = buffer.data(ki + 95);
    const auto *ki_96 = buffer.data(ki + 96);
    const auto *ki_97 = buffer.data(ki + 97);
    const auto *ki_98 = buffer.data(ki + 98);
    const auto *ki_99 = buffer.data(ki + 99);
    const auto *ki_100 = buffer.data(ki + 100);
    const auto *ki_101 = buffer.data(ki + 101);
    const auto *ki_102 = buffer.data(ki + 102);
    const auto *ki_103 = buffer.data(ki + 103);
    const auto *ki_104 = buffer.data(ki + 104);
    const auto *ki_105 = buffer.data(ki + 105);
    const auto *ki_106 = buffer.data(ki + 106);
    const auto *ki_107 = buffer.data(ki + 107);
    const auto *ki_108 = buffer.data(ki + 108);
    const auto *ki_109 = buffer.data(ki + 109);
    const auto *ki_110 = buffer.data(ki + 110);
    const auto *ki_111 = buffer.data(ki + 111);
    const auto *ki_112 = buffer.data(ki + 112);
    const auto *ki_113 = buffer.data(ki + 113);
    const auto *ki_114 = buffer.data(ki + 114);
    const auto *ki_115 = buffer.data(ki + 115);
    const auto *ki_116 = buffer.data(ki + 116);
    const auto *ki_117 = buffer.data(ki + 117);
    const auto *ki_118 = buffer.data(ki + 118);
    const auto *ki_119 = buffer.data(ki + 119);
    const auto *ki_120 = buffer.data(ki + 120);
    const auto *ki_121 = buffer.data(ki + 121);
    const auto *ki_122 = buffer.data(ki + 122);
    const auto *ki_123 = buffer.data(ki + 123);
    const auto *ki_124 = buffer.data(ki + 124);
    const auto *ki_125 = buffer.data(ki + 125);
    const auto *ki_126 = buffer.data(ki + 126);
    const auto *ki_127 = buffer.data(ki + 127);
    const auto *ki_128 = buffer.data(ki + 128);
    const auto *ki_129 = buffer.data(ki + 129);
    const auto *ki_130 = buffer.data(ki + 130);
    const auto *ki_131 = buffer.data(ki + 131);
    const auto *ki_132 = buffer.data(ki + 132);
    const auto *ki_133 = buffer.data(ki + 133);
    const auto *ki_134 = buffer.data(ki + 134);
    const auto *ki_135 = buffer.data(ki + 135);
    const auto *ki_136 = buffer.data(ki + 136);
    const auto *ki_137 = buffer.data(ki + 137);
    const auto *ki_138 = buffer.data(ki + 138);
    const auto *ki_139 = buffer.data(ki + 139);
    const auto *ki_140 = buffer.data(ki + 140);
    const auto *ki_141 = buffer.data(ki + 141);
    const auto *ki_142 = buffer.data(ki + 142);
    const auto *ki_143 = buffer.data(ki + 143);
    const auto *ki_144 = buffer.data(ki + 144);
    const auto *ki_145 = buffer.data(ki + 145);
    const auto *ki_146 = buffer.data(ki + 146);
    const auto *ki_147 = buffer.data(ki + 147);
    const auto *ki_148 = buffer.data(ki + 148);
    const auto *ki_149 = buffer.data(ki + 149);
    const auto *ki_150 = buffer.data(ki + 150);
    const auto *ki_151 = buffer.data(ki + 151);
    const auto *ki_152 = buffer.data(ki + 152);
    const auto *ki_153 = buffer.data(ki + 153);
    const auto *ki_154 = buffer.data(ki + 154);
    const auto *ki_155 = buffer.data(ki + 155);
    const auto *ki_156 = buffer.data(ki + 156);
    const auto *ki_157 = buffer.data(ki + 157);
    const auto *ki_158 = buffer.data(ki + 158);
    const auto *ki_159 = buffer.data(ki + 159);
    const auto *ki_160 = buffer.data(ki + 160);
    const auto *ki_161 = buffer.data(ki + 161);
    const auto *ki_162 = buffer.data(ki + 162);
    const auto *ki_163 = buffer.data(ki + 163);
    const auto *ki_164 = buffer.data(ki + 164);
    const auto *ki_165 = buffer.data(ki + 165);
    const auto *ki_166 = buffer.data(ki + 166);
    const auto *ki_167 = buffer.data(ki + 167);
    const auto *ki_168 = buffer.data(ki + 168);
    const auto *ki_169 = buffer.data(ki + 169);
    const auto *ki_170 = buffer.data(ki + 170);
    const auto *ki_171 = buffer.data(ki + 171);
    const auto *ki_172 = buffer.data(ki + 172);
    const auto *ki_173 = buffer.data(ki + 173);
    const auto *ki_174 = buffer.data(ki + 174);
    const auto *ki_175 = buffer.data(ki + 175);
    const auto *ki_176 = buffer.data(ki + 176);
    const auto *ki_177 = buffer.data(ki + 177);
    const auto *ki_178 = buffer.data(ki + 178);
    const auto *ki_179 = buffer.data(ki + 179);
    const auto *ki_180 = buffer.data(ki + 180);
    const auto *ki_181 = buffer.data(ki + 181);
    const auto *ki_182 = buffer.data(ki + 182);
    const auto *ki_183 = buffer.data(ki + 183);
    const auto *ki_184 = buffer.data(ki + 184);
    const auto *ki_185 = buffer.data(ki + 185);
    const auto *ki_186 = buffer.data(ki + 186);
    const auto *ki_187 = buffer.data(ki + 187);
    const auto *ki_188 = buffer.data(ki + 188);
    const auto *ki_189 = buffer.data(ki + 189);
    const auto *ki_190 = buffer.data(ki + 190);
    const auto *ki_191 = buffer.data(ki + 191);
    const auto *ki_192 = buffer.data(ki + 192);
    const auto *ki_193 = buffer.data(ki + 193);
    const auto *ki_194 = buffer.data(ki + 194);
    const auto *ki_195 = buffer.data(ki + 195);
    const auto *ki_196 = buffer.data(ki + 196);
    const auto *ki_197 = buffer.data(ki + 197);
    const auto *ki_198 = buffer.data(ki + 198);
    const auto *ki_199 = buffer.data(ki + 199);
    const auto *ki_200 = buffer.data(ki + 200);
    const auto *ki_201 = buffer.data(ki + 201);
    const auto *ki_202 = buffer.data(ki + 202);
    const auto *ki_203 = buffer.data(ki + 203);
    const auto *ki_204 = buffer.data(ki + 204);
    const auto *ki_205 = buffer.data(ki + 205);
    const auto *ki_206 = buffer.data(ki + 206);
    const auto *ki_207 = buffer.data(ki + 207);
    const auto *ki_208 = buffer.data(ki + 208);
    const auto *ki_209 = buffer.data(ki + 209);
    const auto *ki_210 = buffer.data(ki + 210);
    const auto *ki_211 = buffer.data(ki + 211);
    const auto *ki_212 = buffer.data(ki + 212);
    const auto *ki_213 = buffer.data(ki + 213);
    const auto *ki_214 = buffer.data(ki + 214);
    const auto *ki_215 = buffer.data(ki + 215);
    const auto *ki_216 = buffer.data(ki + 216);
    const auto *ki_217 = buffer.data(ki + 217);
    const auto *ki_218 = buffer.data(ki + 218);
    const auto *ki_219 = buffer.data(ki + 219);
    const auto *ki_220 = buffer.data(ki + 220);
    const auto *ki_221 = buffer.data(ki + 221);
    const auto *ki_222 = buffer.data(ki + 222);

    const auto *mi_287 = buffer.data(mi + 287);
    const auto *mi_288 = buffer.data(mi + 288);
    const auto *mi_289 = buffer.data(mi + 289);
    const auto *mi_290 = buffer.data(mi + 290);
    const auto *mi_291 = buffer.data(mi + 291);
    const auto *mi_292 = buffer.data(mi + 292);
    const auto *mi_293 = buffer.data(mi + 293);
    const auto *mi_294 = buffer.data(mi + 294);
    const auto *mi_295 = buffer.data(mi + 295);
    const auto *mi_296 = buffer.data(mi + 296);
    const auto *mi_297 = buffer.data(mi + 297);
    const auto *mi_298 = buffer.data(mi + 298);
    const auto *mi_299 = buffer.data(mi + 299);
    const auto *mi_300 = buffer.data(mi + 300);
    const auto *mi_301 = buffer.data(mi + 301);
    const auto *mi_302 = buffer.data(mi + 302);
    const auto *mi_303 = buffer.data(mi + 303);
    const auto *mi_304 = buffer.data(mi + 304);
    const auto *mi_305 = buffer.data(mi + 305);
    const auto *mi_306 = buffer.data(mi + 306);
    const auto *mi_307 = buffer.data(mi + 307);
    const auto *mi_308 = buffer.data(mi + 308);
    const auto *mi_309 = buffer.data(mi + 309);
    const auto *mi_310 = buffer.data(mi + 310);
    const auto *mi_311 = buffer.data(mi + 311);
    const auto *mi_312 = buffer.data(mi + 312);
    const auto *mi_313 = buffer.data(mi + 313);
    const auto *mi_314 = buffer.data(mi + 314);
    const auto *mi_315 = buffer.data(mi + 315);
    const auto *mi_316 = buffer.data(mi + 316);
    const auto *mi_317 = buffer.data(mi + 317);
    const auto *mi_318 = buffer.data(mi + 318);
    const auto *mi_319 = buffer.data(mi + 319);
    const auto *mi_320 = buffer.data(mi + 320);
    const auto *mi_321 = buffer.data(mi + 321);
    const auto *mi_322 = buffer.data(mi + 322);
    const auto *mi_323 = buffer.data(mi + 323);
    const auto *mi_324 = buffer.data(mi + 324);
    const auto *mi_325 = buffer.data(mi + 325);
    const auto *mi_326 = buffer.data(mi + 326);
    const auto *mi_327 = buffer.data(mi + 327);
    const auto *mi_328 = buffer.data(mi + 328);
    const auto *mi_329 = buffer.data(mi + 329);
    const auto *mi_330 = buffer.data(mi + 330);
    const auto *mi_331 = buffer.data(mi + 331);
    const auto *mi_332 = buffer.data(mi + 332);
    const auto *mi_333 = buffer.data(mi + 333);
    const auto *mi_334 = buffer.data(mi + 334);
    const auto *mi_335 = buffer.data(mi + 335);
    const auto *mi_336 = buffer.data(mi + 336);
    const auto *mi_337 = buffer.data(mi + 337);
    const auto *mi_338 = buffer.data(mi + 338);
    const auto *mi_339 = buffer.data(mi + 339);
    const auto *mi_340 = buffer.data(mi + 340);
    const auto *mi_341 = buffer.data(mi + 341);
    const auto *mi_342 = buffer.data(mi + 342);
    const auto *mi_343 = buffer.data(mi + 343);
    const auto *mi_344 = buffer.data(mi + 344);
    const auto *mi_345 = buffer.data(mi + 345);
    const auto *mi_346 = buffer.data(mi + 346);
    const auto *mi_347 = buffer.data(mi + 347);
    const auto *mi_348 = buffer.data(mi + 348);
    const auto *mi_349 = buffer.data(mi + 349);
    const auto *mi_350 = buffer.data(mi + 350);
    const auto *mi_351 = buffer.data(mi + 351);
    const auto *mi_352 = buffer.data(mi + 352);
    const auto *mi_353 = buffer.data(mi + 353);
    const auto *mi_354 = buffer.data(mi + 354);
    const auto *mi_355 = buffer.data(mi + 355);
    const auto *mi_356 = buffer.data(mi + 356);
    const auto *mi_357 = buffer.data(mi + 357);
    const auto *mi_358 = buffer.data(mi + 358);
    const auto *mi_359 = buffer.data(mi + 359);
    const auto *mi_360 = buffer.data(mi + 360);
    const auto *mi_361 = buffer.data(mi + 361);
    const auto *mi_362 = buffer.data(mi + 362);
    const auto *mi_363 = buffer.data(mi + 363);
    const auto *mi_364 = buffer.data(mi + 364);
    const auto *mi_365 = buffer.data(mi + 365);
    const auto *mi_366 = buffer.data(mi + 366);
    const auto *mi_367 = buffer.data(mi + 367);
    const auto *mi_368 = buffer.data(mi + 368);
    const auto *mi_369 = buffer.data(mi + 369);
    const auto *mi_370 = buffer.data(mi + 370);
    const auto *mi_371 = buffer.data(mi + 371);
    const auto *mi_372 = buffer.data(mi + 372);
    const auto *mi_373 = buffer.data(mi + 373);
    const auto *mi_374 = buffer.data(mi + 374);
    const auto *mi_375 = buffer.data(mi + 375);
    const auto *mi_376 = buffer.data(mi + 376);
    const auto *mi_377 = buffer.data(mi + 377);
    const auto *mi_378 = buffer.data(mi + 378);
    const auto *mi_379 = buffer.data(mi + 379);
    const auto *mi_380 = buffer.data(mi + 380);
    const auto *mi_381 = buffer.data(mi + 381);
    const auto *mi_382 = buffer.data(mi + 382);
    const auto *mi_383 = buffer.data(mi + 383);
    const auto *mi_384 = buffer.data(mi + 384);
    const auto *mi_385 = buffer.data(mi + 385);
    const auto *mi_386 = buffer.data(mi + 386);
    const auto *mi_387 = buffer.data(mi + 387);
    const auto *mi_388 = buffer.data(mi + 388);
    const auto *mi_389 = buffer.data(mi + 389);
    const auto *mi_390 = buffer.data(mi + 390);
    const auto *mi_391 = buffer.data(mi + 391);
    const auto *mi_420 = buffer.data(mi + 420);
    const auto *mi_421 = buffer.data(mi + 421);
    const auto *mi_422 = buffer.data(mi + 422);
    const auto *mi_423 = buffer.data(mi + 423);
    const auto *mi_424 = buffer.data(mi + 424);
    const auto *mi_425 = buffer.data(mi + 425);
    const auto *mi_426 = buffer.data(mi + 426);
    const auto *mi_427 = buffer.data(mi + 427);
    const auto *mi_428 = buffer.data(mi + 428);
    const auto *mi_429 = buffer.data(mi + 429);
    const auto *mi_430 = buffer.data(mi + 430);
    const auto *mi_431 = buffer.data(mi + 431);
    const auto *mi_432 = buffer.data(mi + 432);
    const auto *mi_433 = buffer.data(mi + 433);
    const auto *mi_434 = buffer.data(mi + 434);
    const auto *mi_435 = buffer.data(mi + 435);
    const auto *mi_436 = buffer.data(mi + 436);
    const auto *mi_437 = buffer.data(mi + 437);
    const auto *mi_438 = buffer.data(mi + 438);
    const auto *mi_439 = buffer.data(mi + 439);
    const auto *mi_440 = buffer.data(mi + 440);
    const auto *mi_441 = buffer.data(mi + 441);
    const auto *mi_442 = buffer.data(mi + 442);
    const auto *mi_443 = buffer.data(mi + 443);
    const auto *mi_444 = buffer.data(mi + 444);
    const auto *mi_445 = buffer.data(mi + 445);
    const auto *mi_446 = buffer.data(mi + 446);
    const auto *mi_447 = buffer.data(mi + 447);
    const auto *mi_448 = buffer.data(mi + 448);
    const auto *mi_449 = buffer.data(mi + 449);
    const auto *mi_450 = buffer.data(mi + 450);
    const auto *mi_451 = buffer.data(mi + 451);
    const auto *mi_452 = buffer.data(mi + 452);
    const auto *mi_453 = buffer.data(mi + 453);
    const auto *mi_454 = buffer.data(mi + 454);
    const auto *mi_455 = buffer.data(mi + 455);
    const auto *mi_456 = buffer.data(mi + 456);
    const auto *mi_457 = buffer.data(mi + 457);
    const auto *mi_458 = buffer.data(mi + 458);
    const auto *mi_459 = buffer.data(mi + 459);
    const auto *mi_460 = buffer.data(mi + 460);
    const auto *mi_461 = buffer.data(mi + 461);
    const auto *mi_462 = buffer.data(mi + 462);
    const auto *mi_463 = buffer.data(mi + 463);
    const auto *mi_464 = buffer.data(mi + 464);
    const auto *mi_465 = buffer.data(mi + 465);
    const auto *mi_466 = buffer.data(mi + 466);
    const auto *mi_467 = buffer.data(mi + 467);
    const auto *mi_468 = buffer.data(mi + 468);
    const auto *mi_469 = buffer.data(mi + 469);
    const auto *mi_470 = buffer.data(mi + 470);
    const auto *mi_471 = buffer.data(mi + 471);
    const auto *mi_472 = buffer.data(mi + 472);
    const auto *mi_473 = buffer.data(mi + 473);
    const auto *mi_474 = buffer.data(mi + 474);

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ki_91, ki_92, ki_93, ki_94, ki_95, \
                         mi_287, mi_288, mi_289, mi_290, mi_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * ki_91[k]
                   + f_0 * mi_287[k];

        t_176[k] = -3.0 * ki_92[k]
                   + f_0 * mi_288[k];

        t_177[k] = -3.0 * ki_93[k]
                   + f_0 * mi_289[k];

        t_178[k] = -3.0 * ki_94[k]
                   + f_0 * mi_290[k];

        t_179[k] = -3.0 * ki_95[k]
                   + f_0 * mi_291[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ki_96, ki_97, ki_98, ki_99, \
                         ki_100, mi_292, mi_293, mi_294, mi_295, \
                         mi_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * ki_96[k]
                   + f_0 * mi_292[k];

        t_181[k] = -3.0 * ki_97[k]
                   + f_0 * mi_293[k];

        t_182[k] = -3.0 * ki_98[k]
                   + f_0 * mi_294[k];

        t_183[k] = -3.0 * ki_99[k]
                   + f_0 * mi_295[k];

        t_184[k] = -3.0 * ki_100[k]
                   + f_0 * mi_296[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ki_101, ki_102, ki_103, ki_104, \
                         ki_105, mi_297, mi_298, mi_299, mi_300, \
                         mi_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * ki_101[k]
                   + f_0 * mi_297[k];

        t_186[k] = -3.0 * ki_102[k]
                   + f_0 * mi_298[k];

        t_187[k] = -3.0 * ki_103[k]
                   + f_0 * mi_299[k];

        t_188[k] = -3.0 * ki_104[k]
                   + f_0 * mi_300[k];

        t_189[k] = -3.0 * ki_105[k]
                   + f_0 * mi_301[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ki_106, ki_107, ki_108, ki_109, \
                         ki_110, mi_302, mi_303, mi_304, mi_305, \
                         mi_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * ki_106[k]
                   + f_0 * mi_302[k];

        t_191[k] = -3.0 * ki_107[k]
                   + f_0 * mi_303[k];

        t_192[k] = -3.0 * ki_108[k]
                   + f_0 * mi_304[k];

        t_193[k] = -3.0 * ki_109[k]
                   + f_0 * mi_305[k];

        t_194[k] = -3.0 * ki_110[k]
                   + f_0 * mi_306[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ki_111, ki_112, ki_113, ki_114, \
                         ki_115, mi_307, mi_308, mi_309, mi_310, \
                         mi_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * ki_111[k]
                   + f_0 * mi_307[k];

        t_196[k] = -2.0 * ki_112[k]
                   + f_0 * mi_308[k];

        t_197[k] = -2.0 * ki_113[k]
                   + f_0 * mi_309[k];

        t_198[k] = -2.0 * ki_114[k]
                   + f_0 * mi_310[k];

        t_199[k] = -2.0 * ki_115[k]
                   + f_0 * mi_311[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ki_116, ki_117, ki_118, ki_119, \
                         ki_120, mi_312, mi_313, mi_314, mi_315, \
                         mi_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * ki_116[k]
                   + f_0 * mi_312[k];

        t_201[k] = -2.0 * ki_117[k]
                   + f_0 * mi_313[k];

        t_202[k] = -2.0 * ki_118[k]
                   + f_0 * mi_314[k];

        t_203[k] = -2.0 * ki_119[k]
                   + f_0 * mi_315[k];

        t_204[k] = -2.0 * ki_120[k]
                   + f_0 * mi_316[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ki_121, ki_122, ki_123, ki_124, \
                         ki_125, mi_317, mi_318, mi_319, mi_320, \
                         mi_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * ki_121[k]
                   + f_0 * mi_317[k];

        t_206[k] = -2.0 * ki_122[k]
                   + f_0 * mi_318[k];

        t_207[k] = -2.0 * ki_123[k]
                   + f_0 * mi_319[k];

        t_208[k] = -2.0 * ki_124[k]
                   + f_0 * mi_320[k];

        t_209[k] = -2.0 * ki_125[k]
                   + f_0 * mi_321[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ki_126, ki_127, ki_128, ki_129, \
                         ki_130, mi_322, mi_323, mi_324, mi_325, \
                         mi_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * ki_126[k]
                   + f_0 * mi_322[k];

        t_211[k] = -2.0 * ki_127[k]
                   + f_0 * mi_323[k];

        t_212[k] = -2.0 * ki_128[k]
                   + f_0 * mi_324[k];

        t_213[k] = -2.0 * ki_129[k]
                   + f_0 * mi_325[k];

        t_214[k] = -2.0 * ki_130[k]
                   + f_0 * mi_326[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ki_131, ki_132, ki_133, ki_134, \
                         ki_135, mi_327, mi_328, mi_329, mi_330, \
                         mi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * ki_131[k]
                   + f_0 * mi_327[k];

        t_216[k] = -2.0 * ki_132[k]
                   + f_0 * mi_328[k];

        t_217[k] = -2.0 * ki_133[k]
                   + f_0 * mi_329[k];

        t_218[k] = -2.0 * ki_134[k]
                   + f_0 * mi_330[k];

        t_219[k] = -2.0 * ki_135[k]
                   + f_0 * mi_331[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ki_136, ki_137, ki_138, ki_139, \
                         ki_140, mi_332, mi_333, mi_334, mi_335, \
                         mi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * ki_136[k]
                   + f_0 * mi_332[k];

        t_221[k] = -2.0 * ki_137[k]
                   + f_0 * mi_333[k];

        t_222[k] = -2.0 * ki_138[k]
                   + f_0 * mi_334[k];

        t_223[k] = -2.0 * ki_139[k]
                   + f_0 * mi_335[k];

        t_224[k] = -ki_140[k]
                   + f_0 * mi_336[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ki_141, ki_142, ki_143, ki_144, \
                         ki_145, mi_337, mi_338, mi_339, mi_340, \
                         mi_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -ki_141[k]
                   + f_0 * mi_337[k];

        t_226[k] = -ki_142[k]
                   + f_0 * mi_338[k];

        t_227[k] = -ki_143[k]
                   + f_0 * mi_339[k];

        t_228[k] = -ki_144[k]
                   + f_0 * mi_340[k];

        t_229[k] = -ki_145[k]
                   + f_0 * mi_341[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ki_146, ki_147, ki_148, ki_149, \
                         ki_150, mi_342, mi_343, mi_344, mi_345, \
                         mi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -ki_146[k]
                   + f_0 * mi_342[k];

        t_231[k] = -ki_147[k]
                   + f_0 * mi_343[k];

        t_232[k] = -ki_148[k]
                   + f_0 * mi_344[k];

        t_233[k] = -ki_149[k]
                   + f_0 * mi_345[k];

        t_234[k] = -ki_150[k]
                   + f_0 * mi_346[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ki_151, ki_152, ki_153, ki_154, \
                         ki_155, mi_347, mi_348, mi_349, mi_350, \
                         mi_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -ki_151[k]
                   + f_0 * mi_347[k];

        t_236[k] = -ki_152[k]
                   + f_0 * mi_348[k];

        t_237[k] = -ki_153[k]
                   + f_0 * mi_349[k];

        t_238[k] = -ki_154[k]
                   + f_0 * mi_350[k];

        t_239[k] = -ki_155[k]
                   + f_0 * mi_351[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ki_156, ki_157, ki_158, ki_159, \
                         ki_160, mi_352, mi_353, mi_354, mi_355, \
                         mi_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -ki_156[k]
                   + f_0 * mi_352[k];

        t_241[k] = -ki_157[k]
                   + f_0 * mi_353[k];

        t_242[k] = -ki_158[k]
                   + f_0 * mi_354[k];

        t_243[k] = -ki_159[k]
                   + f_0 * mi_355[k];

        t_244[k] = -ki_160[k]
                   + f_0 * mi_356[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ki_161, ki_162, ki_163, ki_164, \
                         ki_165, mi_357, mi_358, mi_359, mi_360, \
                         mi_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -ki_161[k]
                   + f_0 * mi_357[k];

        t_246[k] = -ki_162[k]
                   + f_0 * mi_358[k];

        t_247[k] = -ki_163[k]
                   + f_0 * mi_359[k];

        t_248[k] = -ki_164[k]
                   + f_0 * mi_360[k];

        t_249[k] = -ki_165[k]
                   + f_0 * mi_361[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, t_255, t_256, ki_166, ki_167, \
                         mi_362, mi_363, mi_364, mi_365, mi_366, mi_367, \
                         mi_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -ki_166[k]
                   + f_0 * mi_362[k];

        t_251[k] = -ki_167[k]
                   + f_0 * mi_363[k];

        t_252[k] = f_0 * mi_364[k];

        t_253[k] = f_0 * mi_365[k];

        t_254[k] = f_0 * mi_366[k];

        t_255[k] = f_0 * mi_367[k];

        t_256[k] = f_0 * mi_368[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, t_262, t_263, t_264, mi_369, \
                         mi_370, mi_371, mi_372, mi_373, mi_374, mi_375, \
                         mi_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_0 * mi_369[k];

        t_258[k] = f_0 * mi_370[k];

        t_259[k] = f_0 * mi_371[k];

        t_260[k] = f_0 * mi_372[k];

        t_261[k] = f_0 * mi_373[k];

        t_262[k] = f_0 * mi_374[k];

        t_263[k] = f_0 * mi_375[k];

        t_264[k] = f_0 * mi_376[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, t_270, t_271, t_272, mi_377, \
                         mi_378, mi_379, mi_380, mi_381, mi_382, mi_383, \
                         mi_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_0 * mi_377[k];

        t_266[k] = f_0 * mi_378[k];

        t_267[k] = f_0 * mi_379[k];

        t_268[k] = f_0 * mi_380[k];

        t_269[k] = f_0 * mi_381[k];

        t_270[k] = f_0 * mi_382[k];

        t_271[k] = f_0 * mi_383[k];

        t_272[k] = f_0 * mi_384[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, t_278, t_279, mi_385, mi_386, \
                         mi_387, mi_388, mi_389, mi_390, mi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_0 * mi_385[k];

        t_274[k] = f_0 * mi_386[k];

        t_275[k] = f_0 * mi_387[k];

        t_276[k] = f_0 * mi_388[k];

        t_277[k] = f_0 * mi_389[k];

        t_278[k] = f_0 * mi_390[k];

        t_279[k] = f_0 * mi_391[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ki_168, ki_169, ki_170, ki_171, \
                         ki_172, mi_420, mi_421, mi_422, mi_423, \
                         mi_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -4.0 * ki_168[k]
                   + f_0 * mi_420[k];

        t_281[k] = -4.0 * ki_169[k]
                   + f_0 * mi_421[k];

        t_282[k] = -4.0 * ki_170[k]
                   + f_0 * mi_422[k];

        t_283[k] = -4.0 * ki_171[k]
                   + f_0 * mi_423[k];

        t_284[k] = -4.0 * ki_172[k]
                   + f_0 * mi_424[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ki_173, ki_174, ki_175, ki_176, \
                         ki_177, mi_425, mi_426, mi_427, mi_428, \
                         mi_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -4.0 * ki_173[k]
                   + f_0 * mi_425[k];

        t_286[k] = -4.0 * ki_174[k]
                   + f_0 * mi_426[k];

        t_287[k] = -4.0 * ki_175[k]
                   + f_0 * mi_427[k];

        t_288[k] = -4.0 * ki_176[k]
                   + f_0 * mi_428[k];

        t_289[k] = -4.0 * ki_177[k]
                   + f_0 * mi_429[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ki_178, ki_179, ki_180, ki_181, \
                         ki_182, mi_430, mi_431, mi_432, mi_433, \
                         mi_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -4.0 * ki_178[k]
                   + f_0 * mi_430[k];

        t_291[k] = -4.0 * ki_179[k]
                   + f_0 * mi_431[k];

        t_292[k] = -4.0 * ki_180[k]
                   + f_0 * mi_432[k];

        t_293[k] = -4.0 * ki_181[k]
                   + f_0 * mi_433[k];

        t_294[k] = -4.0 * ki_182[k]
                   + f_0 * mi_434[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ki_183, ki_184, ki_185, ki_186, \
                         ki_187, mi_435, mi_436, mi_437, mi_438, \
                         mi_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * ki_183[k]
                   + f_0 * mi_435[k];

        t_296[k] = -4.0 * ki_184[k]
                   + f_0 * mi_436[k];

        t_297[k] = -4.0 * ki_185[k]
                   + f_0 * mi_437[k];

        t_298[k] = -4.0 * ki_186[k]
                   + f_0 * mi_438[k];

        t_299[k] = -4.0 * ki_187[k]
                   + f_0 * mi_439[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ki_188, ki_189, ki_190, ki_191, \
                         ki_192, mi_440, mi_441, mi_442, mi_443, \
                         mi_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * ki_188[k]
                   + f_0 * mi_440[k];

        t_301[k] = -4.0 * ki_189[k]
                   + f_0 * mi_441[k];

        t_302[k] = -4.0 * ki_190[k]
                   + f_0 * mi_442[k];

        t_303[k] = -4.0 * ki_191[k]
                   + f_0 * mi_443[k];

        t_304[k] = -4.0 * ki_192[k]
                   + f_0 * mi_444[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ki_193, ki_194, ki_195, ki_196, \
                         ki_197, mi_445, mi_446, mi_447, mi_448, \
                         mi_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * ki_193[k]
                   + f_0 * mi_445[k];

        t_306[k] = -4.0 * ki_194[k]
                   + f_0 * mi_446[k];

        t_307[k] = -4.0 * ki_195[k]
                   + f_0 * mi_447[k];

        t_308[k] = -3.0 * ki_196[k]
                   + f_0 * mi_448[k];

        t_309[k] = -3.0 * ki_197[k]
                   + f_0 * mi_449[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ki_198, ki_199, ki_200, ki_201, \
                         ki_202, mi_450, mi_451, mi_452, mi_453, \
                         mi_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * ki_198[k]
                   + f_0 * mi_450[k];

        t_311[k] = -3.0 * ki_199[k]
                   + f_0 * mi_451[k];

        t_312[k] = -3.0 * ki_200[k]
                   + f_0 * mi_452[k];

        t_313[k] = -3.0 * ki_201[k]
                   + f_0 * mi_453[k];

        t_314[k] = -3.0 * ki_202[k]
                   + f_0 * mi_454[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ki_203, ki_204, ki_205, ki_206, \
                         ki_207, mi_455, mi_456, mi_457, mi_458, \
                         mi_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -3.0 * ki_203[k]
                   + f_0 * mi_455[k];

        t_316[k] = -3.0 * ki_204[k]
                   + f_0 * mi_456[k];

        t_317[k] = -3.0 * ki_205[k]
                   + f_0 * mi_457[k];

        t_318[k] = -3.0 * ki_206[k]
                   + f_0 * mi_458[k];

        t_319[k] = -3.0 * ki_207[k]
                   + f_0 * mi_459[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ki_208, ki_209, ki_210, ki_211, \
                         ki_212, mi_460, mi_461, mi_462, mi_463, \
                         mi_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -3.0 * ki_208[k]
                   + f_0 * mi_460[k];

        t_321[k] = -3.0 * ki_209[k]
                   + f_0 * mi_461[k];

        t_322[k] = -3.0 * ki_210[k]
                   + f_0 * mi_462[k];

        t_323[k] = -3.0 * ki_211[k]
                   + f_0 * mi_463[k];

        t_324[k] = -3.0 * ki_212[k]
                   + f_0 * mi_464[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ki_213, ki_214, ki_215, ki_216, \
                         ki_217, mi_465, mi_466, mi_467, mi_468, \
                         mi_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -3.0 * ki_213[k]
                   + f_0 * mi_465[k];

        t_326[k] = -3.0 * ki_214[k]
                   + f_0 * mi_466[k];

        t_327[k] = -3.0 * ki_215[k]
                   + f_0 * mi_467[k];

        t_328[k] = -3.0 * ki_216[k]
                   + f_0 * mi_468[k];

        t_329[k] = -3.0 * ki_217[k]
                   + f_0 * mi_469[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ki_218, ki_219, ki_220, ki_221, \
                         ki_222, mi_470, mi_471, mi_472, mi_473, \
                         mi_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -3.0 * ki_218[k]
                   + f_0 * mi_470[k];

        t_331[k] = -3.0 * ki_219[k]
                   + f_0 * mi_471[k];

        t_332[k] = -3.0 * ki_220[k]
                   + f_0 * mi_472[k];

        t_333[k] = -3.0 * ki_221[k]
                   + f_0 * mi_473[k];

        t_334[k] = -3.0 * ki_222[k]
                   + f_0 * mi_474[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);

    const auto *ki_223 = buffer.data(ki + 223);
    const auto *ki_224 = buffer.data(ki + 224);
    const auto *ki_225 = buffer.data(ki + 225);
    const auto *ki_226 = buffer.data(ki + 226);
    const auto *ki_227 = buffer.data(ki + 227);
    const auto *ki_228 = buffer.data(ki + 228);
    const auto *ki_229 = buffer.data(ki + 229);
    const auto *ki_230 = buffer.data(ki + 230);
    const auto *ki_231 = buffer.data(ki + 231);
    const auto *ki_232 = buffer.data(ki + 232);
    const auto *ki_233 = buffer.data(ki + 233);
    const auto *ki_234 = buffer.data(ki + 234);
    const auto *ki_235 = buffer.data(ki + 235);
    const auto *ki_236 = buffer.data(ki + 236);
    const auto *ki_237 = buffer.data(ki + 237);
    const auto *ki_238 = buffer.data(ki + 238);
    const auto *ki_239 = buffer.data(ki + 239);
    const auto *ki_240 = buffer.data(ki + 240);
    const auto *ki_241 = buffer.data(ki + 241);
    const auto *ki_242 = buffer.data(ki + 242);
    const auto *ki_243 = buffer.data(ki + 243);
    const auto *ki_244 = buffer.data(ki + 244);
    const auto *ki_245 = buffer.data(ki + 245);
    const auto *ki_246 = buffer.data(ki + 246);
    const auto *ki_247 = buffer.data(ki + 247);
    const auto *ki_248 = buffer.data(ki + 248);
    const auto *ki_249 = buffer.data(ki + 249);
    const auto *ki_250 = buffer.data(ki + 250);
    const auto *ki_251 = buffer.data(ki + 251);
    const auto *ki_252 = buffer.data(ki + 252);
    const auto *ki_253 = buffer.data(ki + 253);
    const auto *ki_254 = buffer.data(ki + 254);
    const auto *ki_255 = buffer.data(ki + 255);
    const auto *ki_256 = buffer.data(ki + 256);
    const auto *ki_257 = buffer.data(ki + 257);
    const auto *ki_258 = buffer.data(ki + 258);
    const auto *ki_259 = buffer.data(ki + 259);
    const auto *ki_260 = buffer.data(ki + 260);
    const auto *ki_261 = buffer.data(ki + 261);
    const auto *ki_262 = buffer.data(ki + 262);
    const auto *ki_263 = buffer.data(ki + 263);
    const auto *ki_264 = buffer.data(ki + 264);
    const auto *ki_265 = buffer.data(ki + 265);
    const auto *ki_266 = buffer.data(ki + 266);
    const auto *ki_267 = buffer.data(ki + 267);
    const auto *ki_268 = buffer.data(ki + 268);
    const auto *ki_269 = buffer.data(ki + 269);
    const auto *ki_270 = buffer.data(ki + 270);
    const auto *ki_271 = buffer.data(ki + 271);
    const auto *ki_272 = buffer.data(ki + 272);
    const auto *ki_273 = buffer.data(ki + 273);
    const auto *ki_274 = buffer.data(ki + 274);
    const auto *ki_275 = buffer.data(ki + 275);
    const auto *ki_276 = buffer.data(ki + 276);
    const auto *ki_277 = buffer.data(ki + 277);
    const auto *ki_278 = buffer.data(ki + 278);
    const auto *ki_279 = buffer.data(ki + 279);
    const auto *ki_280 = buffer.data(ki + 280);
    const auto *ki_281 = buffer.data(ki + 281);
    const auto *ki_282 = buffer.data(ki + 282);
    const auto *ki_283 = buffer.data(ki + 283);
    const auto *ki_284 = buffer.data(ki + 284);
    const auto *ki_285 = buffer.data(ki + 285);
    const auto *ki_286 = buffer.data(ki + 286);
    const auto *ki_287 = buffer.data(ki + 287);
    const auto *ki_288 = buffer.data(ki + 288);
    const auto *ki_289 = buffer.data(ki + 289);
    const auto *ki_290 = buffer.data(ki + 290);
    const auto *ki_291 = buffer.data(ki + 291);
    const auto *ki_292 = buffer.data(ki + 292);
    const auto *ki_293 = buffer.data(ki + 293);
    const auto *ki_294 = buffer.data(ki + 294);
    const auto *ki_295 = buffer.data(ki + 295);
    const auto *ki_296 = buffer.data(ki + 296);
    const auto *ki_297 = buffer.data(ki + 297);
    const auto *ki_298 = buffer.data(ki + 298);
    const auto *ki_299 = buffer.data(ki + 299);
    const auto *ki_300 = buffer.data(ki + 300);
    const auto *ki_301 = buffer.data(ki + 301);
    const auto *ki_302 = buffer.data(ki + 302);
    const auto *ki_303 = buffer.data(ki + 303);
    const auto *ki_304 = buffer.data(ki + 304);
    const auto *ki_305 = buffer.data(ki + 305);
    const auto *ki_306 = buffer.data(ki + 306);
    const auto *ki_307 = buffer.data(ki + 307);
    const auto *ki_308 = buffer.data(ki + 308);
    const auto *ki_309 = buffer.data(ki + 309);
    const auto *ki_310 = buffer.data(ki + 310);
    const auto *ki_311 = buffer.data(ki + 311);
    const auto *ki_312 = buffer.data(ki + 312);
    const auto *ki_313 = buffer.data(ki + 313);
    const auto *ki_314 = buffer.data(ki + 314);
    const auto *ki_315 = buffer.data(ki + 315);
    const auto *ki_316 = buffer.data(ki + 316);
    const auto *ki_317 = buffer.data(ki + 317);
    const auto *ki_318 = buffer.data(ki + 318);
    const auto *ki_319 = buffer.data(ki + 319);
    const auto *ki_320 = buffer.data(ki + 320);
    const auto *ki_321 = buffer.data(ki + 321);
    const auto *ki_322 = buffer.data(ki + 322);
    const auto *ki_323 = buffer.data(ki + 323);
    const auto *ki_324 = buffer.data(ki + 324);
    const auto *ki_325 = buffer.data(ki + 325);
    const auto *ki_326 = buffer.data(ki + 326);
    const auto *ki_327 = buffer.data(ki + 327);
    const auto *ki_328 = buffer.data(ki + 328);
    const auto *ki_329 = buffer.data(ki + 329);
    const auto *ki_330 = buffer.data(ki + 330);
    const auto *ki_331 = buffer.data(ki + 331);
    const auto *ki_332 = buffer.data(ki + 332);
    const auto *ki_333 = buffer.data(ki + 333);
    const auto *ki_334 = buffer.data(ki + 334);
    const auto *ki_335 = buffer.data(ki + 335);
    const auto *ki_336 = buffer.data(ki + 336);
    const auto *ki_337 = buffer.data(ki + 337);
    const auto *ki_338 = buffer.data(ki + 338);
    const auto *ki_339 = buffer.data(ki + 339);
    const auto *ki_340 = buffer.data(ki + 340);
    const auto *ki_341 = buffer.data(ki + 341);
    const auto *ki_342 = buffer.data(ki + 342);
    const auto *ki_343 = buffer.data(ki + 343);
    const auto *ki_344 = buffer.data(ki + 344);
    const auto *ki_345 = buffer.data(ki + 345);
    const auto *ki_346 = buffer.data(ki + 346);
    const auto *ki_347 = buffer.data(ki + 347);
    const auto *ki_348 = buffer.data(ki + 348);
    const auto *ki_349 = buffer.data(ki + 349);
    const auto *ki_350 = buffer.data(ki + 350);
    const auto *ki_351 = buffer.data(ki + 351);
    const auto *ki_352 = buffer.data(ki + 352);
    const auto *ki_353 = buffer.data(ki + 353);
    const auto *ki_354 = buffer.data(ki + 354);

    const auto *mi_475 = buffer.data(mi + 475);
    const auto *mi_476 = buffer.data(mi + 476);
    const auto *mi_477 = buffer.data(mi + 477);
    const auto *mi_478 = buffer.data(mi + 478);
    const auto *mi_479 = buffer.data(mi + 479);
    const auto *mi_480 = buffer.data(mi + 480);
    const auto *mi_481 = buffer.data(mi + 481);
    const auto *mi_482 = buffer.data(mi + 482);
    const auto *mi_483 = buffer.data(mi + 483);
    const auto *mi_484 = buffer.data(mi + 484);
    const auto *mi_485 = buffer.data(mi + 485);
    const auto *mi_486 = buffer.data(mi + 486);
    const auto *mi_487 = buffer.data(mi + 487);
    const auto *mi_488 = buffer.data(mi + 488);
    const auto *mi_489 = buffer.data(mi + 489);
    const auto *mi_490 = buffer.data(mi + 490);
    const auto *mi_491 = buffer.data(mi + 491);
    const auto *mi_492 = buffer.data(mi + 492);
    const auto *mi_493 = buffer.data(mi + 493);
    const auto *mi_494 = buffer.data(mi + 494);
    const auto *mi_495 = buffer.data(mi + 495);
    const auto *mi_496 = buffer.data(mi + 496);
    const auto *mi_497 = buffer.data(mi + 497);
    const auto *mi_498 = buffer.data(mi + 498);
    const auto *mi_499 = buffer.data(mi + 499);
    const auto *mi_500 = buffer.data(mi + 500);
    const auto *mi_501 = buffer.data(mi + 501);
    const auto *mi_502 = buffer.data(mi + 502);
    const auto *mi_503 = buffer.data(mi + 503);
    const auto *mi_504 = buffer.data(mi + 504);
    const auto *mi_505 = buffer.data(mi + 505);
    const auto *mi_506 = buffer.data(mi + 506);
    const auto *mi_507 = buffer.data(mi + 507);
    const auto *mi_508 = buffer.data(mi + 508);
    const auto *mi_509 = buffer.data(mi + 509);
    const auto *mi_510 = buffer.data(mi + 510);
    const auto *mi_511 = buffer.data(mi + 511);
    const auto *mi_512 = buffer.data(mi + 512);
    const auto *mi_513 = buffer.data(mi + 513);
    const auto *mi_514 = buffer.data(mi + 514);
    const auto *mi_515 = buffer.data(mi + 515);
    const auto *mi_516 = buffer.data(mi + 516);
    const auto *mi_517 = buffer.data(mi + 517);
    const auto *mi_518 = buffer.data(mi + 518);
    const auto *mi_519 = buffer.data(mi + 519);
    const auto *mi_520 = buffer.data(mi + 520);
    const auto *mi_521 = buffer.data(mi + 521);
    const auto *mi_522 = buffer.data(mi + 522);
    const auto *mi_523 = buffer.data(mi + 523);
    const auto *mi_524 = buffer.data(mi + 524);
    const auto *mi_525 = buffer.data(mi + 525);
    const auto *mi_526 = buffer.data(mi + 526);
    const auto *mi_527 = buffer.data(mi + 527);
    const auto *mi_528 = buffer.data(mi + 528);
    const auto *mi_529 = buffer.data(mi + 529);
    const auto *mi_530 = buffer.data(mi + 530);
    const auto *mi_531 = buffer.data(mi + 531);
    const auto *mi_532 = buffer.data(mi + 532);
    const auto *mi_533 = buffer.data(mi + 533);
    const auto *mi_534 = buffer.data(mi + 534);
    const auto *mi_535 = buffer.data(mi + 535);
    const auto *mi_536 = buffer.data(mi + 536);
    const auto *mi_537 = buffer.data(mi + 537);
    const auto *mi_538 = buffer.data(mi + 538);
    const auto *mi_539 = buffer.data(mi + 539);
    const auto *mi_540 = buffer.data(mi + 540);
    const auto *mi_541 = buffer.data(mi + 541);
    const auto *mi_542 = buffer.data(mi + 542);
    const auto *mi_543 = buffer.data(mi + 543);
    const auto *mi_544 = buffer.data(mi + 544);
    const auto *mi_545 = buffer.data(mi + 545);
    const auto *mi_546 = buffer.data(mi + 546);
    const auto *mi_547 = buffer.data(mi + 547);
    const auto *mi_548 = buffer.data(mi + 548);
    const auto *mi_549 = buffer.data(mi + 549);
    const auto *mi_550 = buffer.data(mi + 550);
    const auto *mi_551 = buffer.data(mi + 551);
    const auto *mi_552 = buffer.data(mi + 552);
    const auto *mi_553 = buffer.data(mi + 553);
    const auto *mi_554 = buffer.data(mi + 554);
    const auto *mi_555 = buffer.data(mi + 555);
    const auto *mi_556 = buffer.data(mi + 556);
    const auto *mi_557 = buffer.data(mi + 557);
    const auto *mi_558 = buffer.data(mi + 558);
    const auto *mi_559 = buffer.data(mi + 559);
    const auto *mi_588 = buffer.data(mi + 588);
    const auto *mi_589 = buffer.data(mi + 589);
    const auto *mi_590 = buffer.data(mi + 590);
    const auto *mi_591 = buffer.data(mi + 591);
    const auto *mi_592 = buffer.data(mi + 592);
    const auto *mi_593 = buffer.data(mi + 593);
    const auto *mi_594 = buffer.data(mi + 594);
    const auto *mi_595 = buffer.data(mi + 595);
    const auto *mi_596 = buffer.data(mi + 596);
    const auto *mi_597 = buffer.data(mi + 597);
    const auto *mi_598 = buffer.data(mi + 598);
    const auto *mi_599 = buffer.data(mi + 599);
    const auto *mi_600 = buffer.data(mi + 600);
    const auto *mi_601 = buffer.data(mi + 601);
    const auto *mi_602 = buffer.data(mi + 602);
    const auto *mi_603 = buffer.data(mi + 603);
    const auto *mi_604 = buffer.data(mi + 604);
    const auto *mi_605 = buffer.data(mi + 605);
    const auto *mi_606 = buffer.data(mi + 606);
    const auto *mi_607 = buffer.data(mi + 607);
    const auto *mi_608 = buffer.data(mi + 608);
    const auto *mi_609 = buffer.data(mi + 609);
    const auto *mi_610 = buffer.data(mi + 610);
    const auto *mi_611 = buffer.data(mi + 611);
    const auto *mi_612 = buffer.data(mi + 612);
    const auto *mi_613 = buffer.data(mi + 613);
    const auto *mi_614 = buffer.data(mi + 614);
    const auto *mi_615 = buffer.data(mi + 615);
    const auto *mi_616 = buffer.data(mi + 616);
    const auto *mi_617 = buffer.data(mi + 617);
    const auto *mi_618 = buffer.data(mi + 618);
    const auto *mi_619 = buffer.data(mi + 619);
    const auto *mi_620 = buffer.data(mi + 620);
    const auto *mi_621 = buffer.data(mi + 621);
    const auto *mi_622 = buffer.data(mi + 622);
    const auto *mi_623 = buffer.data(mi + 623);
    const auto *mi_624 = buffer.data(mi + 624);
    const auto *mi_625 = buffer.data(mi + 625);
    const auto *mi_626 = buffer.data(mi + 626);
    const auto *mi_627 = buffer.data(mi + 627);
    const auto *mi_628 = buffer.data(mi + 628);
    const auto *mi_629 = buffer.data(mi + 629);
    const auto *mi_630 = buffer.data(mi + 630);
    const auto *mi_631 = buffer.data(mi + 631);
    const auto *mi_632 = buffer.data(mi + 632);
    const auto *mi_633 = buffer.data(mi + 633);
    const auto *mi_634 = buffer.data(mi + 634);
    const auto *mi_635 = buffer.data(mi + 635);
    const auto *mi_636 = buffer.data(mi + 636);
    const auto *mi_637 = buffer.data(mi + 637);
    const auto *mi_638 = buffer.data(mi + 638);
    const auto *mi_639 = buffer.data(mi + 639);
    const auto *mi_640 = buffer.data(mi + 640);
    const auto *mi_641 = buffer.data(mi + 641);
    const auto *mi_642 = buffer.data(mi + 642);
    const auto *mi_643 = buffer.data(mi + 643);
    const auto *mi_644 = buffer.data(mi + 644);
    const auto *mi_645 = buffer.data(mi + 645);
    const auto *mi_646 = buffer.data(mi + 646);
    const auto *mi_647 = buffer.data(mi + 647);
    const auto *mi_648 = buffer.data(mi + 648);
    const auto *mi_649 = buffer.data(mi + 649);
    const auto *mi_650 = buffer.data(mi + 650);
    const auto *mi_651 = buffer.data(mi + 651);
    const auto *mi_652 = buffer.data(mi + 652);
    const auto *mi_653 = buffer.data(mi + 653);
    const auto *mi_654 = buffer.data(mi + 654);
    const auto *mi_655 = buffer.data(mi + 655);
    const auto *mi_656 = buffer.data(mi + 656);
    const auto *mi_657 = buffer.data(mi + 657);
    const auto *mi_658 = buffer.data(mi + 658);
    const auto *mi_659 = buffer.data(mi + 659);
    const auto *mi_660 = buffer.data(mi + 660);
    const auto *mi_661 = buffer.data(mi + 661);
    const auto *mi_662 = buffer.data(mi + 662);

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ki_223, ki_224, ki_225, ki_226, \
                         ki_227, mi_475, mi_476, mi_477, mi_478, \
                         mi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -3.0 * ki_223[k]
                   + f_0 * mi_475[k];

        t_336[k] = -2.0 * ki_224[k]
                   + f_0 * mi_476[k];

        t_337[k] = -2.0 * ki_225[k]
                   + f_0 * mi_477[k];

        t_338[k] = -2.0 * ki_226[k]
                   + f_0 * mi_478[k];

        t_339[k] = -2.0 * ki_227[k]
                   + f_0 * mi_479[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ki_228, ki_229, ki_230, ki_231, \
                         ki_232, mi_480, mi_481, mi_482, mi_483, \
                         mi_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * ki_228[k]
                   + f_0 * mi_480[k];

        t_341[k] = -2.0 * ki_229[k]
                   + f_0 * mi_481[k];

        t_342[k] = -2.0 * ki_230[k]
                   + f_0 * mi_482[k];

        t_343[k] = -2.0 * ki_231[k]
                   + f_0 * mi_483[k];

        t_344[k] = -2.0 * ki_232[k]
                   + f_0 * mi_484[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ki_233, ki_234, ki_235, ki_236, \
                         ki_237, mi_485, mi_486, mi_487, mi_488, \
                         mi_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * ki_233[k]
                   + f_0 * mi_485[k];

        t_346[k] = -2.0 * ki_234[k]
                   + f_0 * mi_486[k];

        t_347[k] = -2.0 * ki_235[k]
                   + f_0 * mi_487[k];

        t_348[k] = -2.0 * ki_236[k]
                   + f_0 * mi_488[k];

        t_349[k] = -2.0 * ki_237[k]
                   + f_0 * mi_489[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ki_238, ki_239, ki_240, ki_241, \
                         ki_242, mi_490, mi_491, mi_492, mi_493, \
                         mi_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * ki_238[k]
                   + f_0 * mi_490[k];

        t_351[k] = -2.0 * ki_239[k]
                   + f_0 * mi_491[k];

        t_352[k] = -2.0 * ki_240[k]
                   + f_0 * mi_492[k];

        t_353[k] = -2.0 * ki_241[k]
                   + f_0 * mi_493[k];

        t_354[k] = -2.0 * ki_242[k]
                   + f_0 * mi_494[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ki_243, ki_244, ki_245, ki_246, \
                         ki_247, mi_495, mi_496, mi_497, mi_498, \
                         mi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * ki_243[k]
                   + f_0 * mi_495[k];

        t_356[k] = -2.0 * ki_244[k]
                   + f_0 * mi_496[k];

        t_357[k] = -2.0 * ki_245[k]
                   + f_0 * mi_497[k];

        t_358[k] = -2.0 * ki_246[k]
                   + f_0 * mi_498[k];

        t_359[k] = -2.0 * ki_247[k]
                   + f_0 * mi_499[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ki_248, ki_249, ki_250, ki_251, \
                         ki_252, mi_500, mi_501, mi_502, mi_503, \
                         mi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -2.0 * ki_248[k]
                   + f_0 * mi_500[k];

        t_361[k] = -2.0 * ki_249[k]
                   + f_0 * mi_501[k];

        t_362[k] = -2.0 * ki_250[k]
                   + f_0 * mi_502[k];

        t_363[k] = -2.0 * ki_251[k]
                   + f_0 * mi_503[k];

        t_364[k] = -ki_252[k]
                   + f_0 * mi_504[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ki_253, ki_254, ki_255, ki_256, \
                         ki_257, mi_505, mi_506, mi_507, mi_508, \
                         mi_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -ki_253[k]
                   + f_0 * mi_505[k];

        t_366[k] = -ki_254[k]
                   + f_0 * mi_506[k];

        t_367[k] = -ki_255[k]
                   + f_0 * mi_507[k];

        t_368[k] = -ki_256[k]
                   + f_0 * mi_508[k];

        t_369[k] = -ki_257[k]
                   + f_0 * mi_509[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ki_258, ki_259, ki_260, ki_261, \
                         ki_262, mi_510, mi_511, mi_512, mi_513, \
                         mi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -ki_258[k]
                   + f_0 * mi_510[k];

        t_371[k] = -ki_259[k]
                   + f_0 * mi_511[k];

        t_372[k] = -ki_260[k]
                   + f_0 * mi_512[k];

        t_373[k] = -ki_261[k]
                   + f_0 * mi_513[k];

        t_374[k] = -ki_262[k]
                   + f_0 * mi_514[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ki_263, ki_264, ki_265, ki_266, \
                         ki_267, mi_515, mi_516, mi_517, mi_518, \
                         mi_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -ki_263[k]
                   + f_0 * mi_515[k];

        t_376[k] = -ki_264[k]
                   + f_0 * mi_516[k];

        t_377[k] = -ki_265[k]
                   + f_0 * mi_517[k];

        t_378[k] = -ki_266[k]
                   + f_0 * mi_518[k];

        t_379[k] = -ki_267[k]
                   + f_0 * mi_519[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ki_268, ki_269, ki_270, ki_271, \
                         ki_272, mi_520, mi_521, mi_522, mi_523, \
                         mi_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -ki_268[k]
                   + f_0 * mi_520[k];

        t_381[k] = -ki_269[k]
                   + f_0 * mi_521[k];

        t_382[k] = -ki_270[k]
                   + f_0 * mi_522[k];

        t_383[k] = -ki_271[k]
                   + f_0 * mi_523[k];

        t_384[k] = -ki_272[k]
                   + f_0 * mi_524[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ki_273, ki_274, ki_275, ki_276, \
                         ki_277, mi_525, mi_526, mi_527, mi_528, \
                         mi_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -ki_273[k]
                   + f_0 * mi_525[k];

        t_386[k] = -ki_274[k]
                   + f_0 * mi_526[k];

        t_387[k] = -ki_275[k]
                   + f_0 * mi_527[k];

        t_388[k] = -ki_276[k]
                   + f_0 * mi_528[k];

        t_389[k] = -ki_277[k]
                   + f_0 * mi_529[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, t_395, t_396, ki_278, ki_279, \
                         mi_530, mi_531, mi_532, mi_533, mi_534, mi_535, \
                         mi_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -ki_278[k]
                   + f_0 * mi_530[k];

        t_391[k] = -ki_279[k]
                   + f_0 * mi_531[k];

        t_392[k] = f_0 * mi_532[k];

        t_393[k] = f_0 * mi_533[k];

        t_394[k] = f_0 * mi_534[k];

        t_395[k] = f_0 * mi_535[k];

        t_396[k] = f_0 * mi_536[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, t_402, t_403, t_404, mi_537, \
                         mi_538, mi_539, mi_540, mi_541, mi_542, mi_543, \
                         mi_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_0 * mi_537[k];

        t_398[k] = f_0 * mi_538[k];

        t_399[k] = f_0 * mi_539[k];

        t_400[k] = f_0 * mi_540[k];

        t_401[k] = f_0 * mi_541[k];

        t_402[k] = f_0 * mi_542[k];

        t_403[k] = f_0 * mi_543[k];

        t_404[k] = f_0 * mi_544[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, mi_545, \
                         mi_546, mi_547, mi_548, mi_549, mi_550, mi_551, \
                         mi_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * mi_545[k];

        t_406[k] = f_0 * mi_546[k];

        t_407[k] = f_0 * mi_547[k];

        t_408[k] = f_0 * mi_548[k];

        t_409[k] = f_0 * mi_549[k];

        t_410[k] = f_0 * mi_550[k];

        t_411[k] = f_0 * mi_551[k];

        t_412[k] = f_0 * mi_552[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, mi_553, mi_554, \
                         mi_555, mi_556, mi_557, mi_558, mi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * mi_553[k];

        t_414[k] = f_0 * mi_554[k];

        t_415[k] = f_0 * mi_555[k];

        t_416[k] = f_0 * mi_556[k];

        t_417[k] = f_0 * mi_557[k];

        t_418[k] = f_0 * mi_558[k];

        t_419[k] = f_0 * mi_559[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ki_280, ki_281, ki_282, ki_283, \
                         ki_284, mi_588, mi_589, mi_590, mi_591, \
                         mi_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -5.0 * ki_280[k]
                   + f_0 * mi_588[k];

        t_421[k] = -5.0 * ki_281[k]
                   + f_0 * mi_589[k];

        t_422[k] = -5.0 * ki_282[k]
                   + f_0 * mi_590[k];

        t_423[k] = -5.0 * ki_283[k]
                   + f_0 * mi_591[k];

        t_424[k] = -5.0 * ki_284[k]
                   + f_0 * mi_592[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ki_285, ki_286, ki_287, ki_288, \
                         ki_289, mi_593, mi_594, mi_595, mi_596, \
                         mi_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -5.0 * ki_285[k]
                   + f_0 * mi_593[k];

        t_426[k] = -5.0 * ki_286[k]
                   + f_0 * mi_594[k];

        t_427[k] = -5.0 * ki_287[k]
                   + f_0 * mi_595[k];

        t_428[k] = -5.0 * ki_288[k]
                   + f_0 * mi_596[k];

        t_429[k] = -5.0 * ki_289[k]
                   + f_0 * mi_597[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ki_290, ki_291, ki_292, ki_293, \
                         ki_294, mi_598, mi_599, mi_600, mi_601, \
                         mi_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -5.0 * ki_290[k]
                   + f_0 * mi_598[k];

        t_431[k] = -5.0 * ki_291[k]
                   + f_0 * mi_599[k];

        t_432[k] = -5.0 * ki_292[k]
                   + f_0 * mi_600[k];

        t_433[k] = -5.0 * ki_293[k]
                   + f_0 * mi_601[k];

        t_434[k] = -5.0 * ki_294[k]
                   + f_0 * mi_602[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ki_295, ki_296, ki_297, ki_298, \
                         ki_299, mi_603, mi_604, mi_605, mi_606, \
                         mi_607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -5.0 * ki_295[k]
                   + f_0 * mi_603[k];

        t_436[k] = -5.0 * ki_296[k]
                   + f_0 * mi_604[k];

        t_437[k] = -5.0 * ki_297[k]
                   + f_0 * mi_605[k];

        t_438[k] = -5.0 * ki_298[k]
                   + f_0 * mi_606[k];

        t_439[k] = -5.0 * ki_299[k]
                   + f_0 * mi_607[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ki_300, ki_301, ki_302, ki_303, \
                         ki_304, mi_608, mi_609, mi_610, mi_611, \
                         mi_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -5.0 * ki_300[k]
                   + f_0 * mi_608[k];

        t_441[k] = -5.0 * ki_301[k]
                   + f_0 * mi_609[k];

        t_442[k] = -5.0 * ki_302[k]
                   + f_0 * mi_610[k];

        t_443[k] = -5.0 * ki_303[k]
                   + f_0 * mi_611[k];

        t_444[k] = -5.0 * ki_304[k]
                   + f_0 * mi_612[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ki_305, ki_306, ki_307, ki_308, \
                         ki_309, mi_613, mi_614, mi_615, mi_616, \
                         mi_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -5.0 * ki_305[k]
                   + f_0 * mi_613[k];

        t_446[k] = -5.0 * ki_306[k]
                   + f_0 * mi_614[k];

        t_447[k] = -5.0 * ki_307[k]
                   + f_0 * mi_615[k];

        t_448[k] = -4.0 * ki_308[k]
                   + f_0 * mi_616[k];

        t_449[k] = -4.0 * ki_309[k]
                   + f_0 * mi_617[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ki_310, ki_311, ki_312, ki_313, \
                         ki_314, mi_618, mi_619, mi_620, mi_621, \
                         mi_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -4.0 * ki_310[k]
                   + f_0 * mi_618[k];

        t_451[k] = -4.0 * ki_311[k]
                   + f_0 * mi_619[k];

        t_452[k] = -4.0 * ki_312[k]
                   + f_0 * mi_620[k];

        t_453[k] = -4.0 * ki_313[k]
                   + f_0 * mi_621[k];

        t_454[k] = -4.0 * ki_314[k]
                   + f_0 * mi_622[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ki_315, ki_316, ki_317, ki_318, \
                         ki_319, mi_623, mi_624, mi_625, mi_626, \
                         mi_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -4.0 * ki_315[k]
                   + f_0 * mi_623[k];

        t_456[k] = -4.0 * ki_316[k]
                   + f_0 * mi_624[k];

        t_457[k] = -4.0 * ki_317[k]
                   + f_0 * mi_625[k];

        t_458[k] = -4.0 * ki_318[k]
                   + f_0 * mi_626[k];

        t_459[k] = -4.0 * ki_319[k]
                   + f_0 * mi_627[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ki_320, ki_321, ki_322, ki_323, \
                         ki_324, mi_628, mi_629, mi_630, mi_631, \
                         mi_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -4.0 * ki_320[k]
                   + f_0 * mi_628[k];

        t_461[k] = -4.0 * ki_321[k]
                   + f_0 * mi_629[k];

        t_462[k] = -4.0 * ki_322[k]
                   + f_0 * mi_630[k];

        t_463[k] = -4.0 * ki_323[k]
                   + f_0 * mi_631[k];

        t_464[k] = -4.0 * ki_324[k]
                   + f_0 * mi_632[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ki_325, ki_326, ki_327, ki_328, \
                         ki_329, mi_633, mi_634, mi_635, mi_636, \
                         mi_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -4.0 * ki_325[k]
                   + f_0 * mi_633[k];

        t_466[k] = -4.0 * ki_326[k]
                   + f_0 * mi_634[k];

        t_467[k] = -4.0 * ki_327[k]
                   + f_0 * mi_635[k];

        t_468[k] = -4.0 * ki_328[k]
                   + f_0 * mi_636[k];

        t_469[k] = -4.0 * ki_329[k]
                   + f_0 * mi_637[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ki_330, ki_331, ki_332, ki_333, \
                         ki_334, mi_638, mi_639, mi_640, mi_641, \
                         mi_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -4.0 * ki_330[k]
                   + f_0 * mi_638[k];

        t_471[k] = -4.0 * ki_331[k]
                   + f_0 * mi_639[k];

        t_472[k] = -4.0 * ki_332[k]
                   + f_0 * mi_640[k];

        t_473[k] = -4.0 * ki_333[k]
                   + f_0 * mi_641[k];

        t_474[k] = -4.0 * ki_334[k]
                   + f_0 * mi_642[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ki_335, ki_336, ki_337, ki_338, \
                         ki_339, mi_643, mi_644, mi_645, mi_646, \
                         mi_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -4.0 * ki_335[k]
                   + f_0 * mi_643[k];

        t_476[k] = -3.0 * ki_336[k]
                   + f_0 * mi_644[k];

        t_477[k] = -3.0 * ki_337[k]
                   + f_0 * mi_645[k];

        t_478[k] = -3.0 * ki_338[k]
                   + f_0 * mi_646[k];

        t_479[k] = -3.0 * ki_339[k]
                   + f_0 * mi_647[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ki_340, ki_341, ki_342, ki_343, \
                         ki_344, mi_648, mi_649, mi_650, mi_651, \
                         mi_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -3.0 * ki_340[k]
                   + f_0 * mi_648[k];

        t_481[k] = -3.0 * ki_341[k]
                   + f_0 * mi_649[k];

        t_482[k] = -3.0 * ki_342[k]
                   + f_0 * mi_650[k];

        t_483[k] = -3.0 * ki_343[k]
                   + f_0 * mi_651[k];

        t_484[k] = -3.0 * ki_344[k]
                   + f_0 * mi_652[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ki_345, ki_346, ki_347, ki_348, \
                         ki_349, mi_653, mi_654, mi_655, mi_656, \
                         mi_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -3.0 * ki_345[k]
                   + f_0 * mi_653[k];

        t_486[k] = -3.0 * ki_346[k]
                   + f_0 * mi_654[k];

        t_487[k] = -3.0 * ki_347[k]
                   + f_0 * mi_655[k];

        t_488[k] = -3.0 * ki_348[k]
                   + f_0 * mi_656[k];

        t_489[k] = -3.0 * ki_349[k]
                   + f_0 * mi_657[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ki_350, ki_351, ki_352, ki_353, \
                         ki_354, mi_658, mi_659, mi_660, mi_661, \
                         mi_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -3.0 * ki_350[k]
                   + f_0 * mi_658[k];

        t_491[k] = -3.0 * ki_351[k]
                   + f_0 * mi_659[k];

        t_492[k] = -3.0 * ki_352[k]
                   + f_0 * mi_660[k];

        t_493[k] = -3.0 * ki_353[k]
                   + f_0 * mi_661[k];

        t_494[k] = -3.0 * ki_354[k]
                   + f_0 * mi_662[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);
    auto *t_654 = buffer.data(target + 654);

    const auto *ki_355 = buffer.data(ki + 355);
    const auto *ki_356 = buffer.data(ki + 356);
    const auto *ki_357 = buffer.data(ki + 357);
    const auto *ki_358 = buffer.data(ki + 358);
    const auto *ki_359 = buffer.data(ki + 359);
    const auto *ki_360 = buffer.data(ki + 360);
    const auto *ki_361 = buffer.data(ki + 361);
    const auto *ki_362 = buffer.data(ki + 362);
    const auto *ki_363 = buffer.data(ki + 363);
    const auto *ki_364 = buffer.data(ki + 364);
    const auto *ki_365 = buffer.data(ki + 365);
    const auto *ki_366 = buffer.data(ki + 366);
    const auto *ki_367 = buffer.data(ki + 367);
    const auto *ki_368 = buffer.data(ki + 368);
    const auto *ki_369 = buffer.data(ki + 369);
    const auto *ki_370 = buffer.data(ki + 370);
    const auto *ki_371 = buffer.data(ki + 371);
    const auto *ki_372 = buffer.data(ki + 372);
    const auto *ki_373 = buffer.data(ki + 373);
    const auto *ki_374 = buffer.data(ki + 374);
    const auto *ki_375 = buffer.data(ki + 375);
    const auto *ki_376 = buffer.data(ki + 376);
    const auto *ki_377 = buffer.data(ki + 377);
    const auto *ki_378 = buffer.data(ki + 378);
    const auto *ki_379 = buffer.data(ki + 379);
    const auto *ki_380 = buffer.data(ki + 380);
    const auto *ki_381 = buffer.data(ki + 381);
    const auto *ki_382 = buffer.data(ki + 382);
    const auto *ki_383 = buffer.data(ki + 383);
    const auto *ki_384 = buffer.data(ki + 384);
    const auto *ki_385 = buffer.data(ki + 385);
    const auto *ki_386 = buffer.data(ki + 386);
    const auto *ki_387 = buffer.data(ki + 387);
    const auto *ki_388 = buffer.data(ki + 388);
    const auto *ki_389 = buffer.data(ki + 389);
    const auto *ki_390 = buffer.data(ki + 390);
    const auto *ki_391 = buffer.data(ki + 391);
    const auto *ki_392 = buffer.data(ki + 392);
    const auto *ki_393 = buffer.data(ki + 393);
    const auto *ki_394 = buffer.data(ki + 394);
    const auto *ki_395 = buffer.data(ki + 395);
    const auto *ki_396 = buffer.data(ki + 396);
    const auto *ki_397 = buffer.data(ki + 397);
    const auto *ki_398 = buffer.data(ki + 398);
    const auto *ki_399 = buffer.data(ki + 399);
    const auto *ki_400 = buffer.data(ki + 400);
    const auto *ki_401 = buffer.data(ki + 401);
    const auto *ki_402 = buffer.data(ki + 402);
    const auto *ki_403 = buffer.data(ki + 403);
    const auto *ki_404 = buffer.data(ki + 404);
    const auto *ki_405 = buffer.data(ki + 405);
    const auto *ki_406 = buffer.data(ki + 406);
    const auto *ki_407 = buffer.data(ki + 407);
    const auto *ki_408 = buffer.data(ki + 408);
    const auto *ki_409 = buffer.data(ki + 409);
    const auto *ki_410 = buffer.data(ki + 410);
    const auto *ki_411 = buffer.data(ki + 411);
    const auto *ki_412 = buffer.data(ki + 412);
    const auto *ki_413 = buffer.data(ki + 413);
    const auto *ki_414 = buffer.data(ki + 414);
    const auto *ki_415 = buffer.data(ki + 415);
    const auto *ki_416 = buffer.data(ki + 416);
    const auto *ki_417 = buffer.data(ki + 417);
    const auto *ki_418 = buffer.data(ki + 418);
    const auto *ki_419 = buffer.data(ki + 419);
    const auto *ki_420 = buffer.data(ki + 420);
    const auto *ki_421 = buffer.data(ki + 421);
    const auto *ki_422 = buffer.data(ki + 422);
    const auto *ki_423 = buffer.data(ki + 423);
    const auto *ki_424 = buffer.data(ki + 424);
    const auto *ki_425 = buffer.data(ki + 425);
    const auto *ki_426 = buffer.data(ki + 426);
    const auto *ki_427 = buffer.data(ki + 427);
    const auto *ki_428 = buffer.data(ki + 428);
    const auto *ki_429 = buffer.data(ki + 429);
    const auto *ki_430 = buffer.data(ki + 430);
    const auto *ki_431 = buffer.data(ki + 431);
    const auto *ki_432 = buffer.data(ki + 432);
    const auto *ki_433 = buffer.data(ki + 433);
    const auto *ki_434 = buffer.data(ki + 434);
    const auto *ki_435 = buffer.data(ki + 435);
    const auto *ki_436 = buffer.data(ki + 436);
    const auto *ki_437 = buffer.data(ki + 437);
    const auto *ki_438 = buffer.data(ki + 438);
    const auto *ki_439 = buffer.data(ki + 439);
    const auto *ki_440 = buffer.data(ki + 440);
    const auto *ki_441 = buffer.data(ki + 441);
    const auto *ki_442 = buffer.data(ki + 442);
    const auto *ki_443 = buffer.data(ki + 443);
    const auto *ki_444 = buffer.data(ki + 444);
    const auto *ki_445 = buffer.data(ki + 445);
    const auto *ki_446 = buffer.data(ki + 446);
    const auto *ki_447 = buffer.data(ki + 447);
    const auto *ki_448 = buffer.data(ki + 448);
    const auto *ki_449 = buffer.data(ki + 449);
    const auto *ki_450 = buffer.data(ki + 450);
    const auto *ki_451 = buffer.data(ki + 451);
    const auto *ki_452 = buffer.data(ki + 452);
    const auto *ki_453 = buffer.data(ki + 453);
    const auto *ki_454 = buffer.data(ki + 454);
    const auto *ki_455 = buffer.data(ki + 455);
    const auto *ki_456 = buffer.data(ki + 456);
    const auto *ki_457 = buffer.data(ki + 457);
    const auto *ki_458 = buffer.data(ki + 458);
    const auto *ki_459 = buffer.data(ki + 459);
    const auto *ki_460 = buffer.data(ki + 460);
    const auto *ki_461 = buffer.data(ki + 461);
    const auto *ki_462 = buffer.data(ki + 462);
    const auto *ki_463 = buffer.data(ki + 463);
    const auto *ki_464 = buffer.data(ki + 464);
    const auto *ki_465 = buffer.data(ki + 465);
    const auto *ki_466 = buffer.data(ki + 466);
    const auto *ki_467 = buffer.data(ki + 467);
    const auto *ki_468 = buffer.data(ki + 468);
    const auto *ki_469 = buffer.data(ki + 469);
    const auto *ki_470 = buffer.data(ki + 470);
    const auto *ki_471 = buffer.data(ki + 471);
    const auto *ki_472 = buffer.data(ki + 472);
    const auto *ki_473 = buffer.data(ki + 473);
    const auto *ki_474 = buffer.data(ki + 474);
    const auto *ki_475 = buffer.data(ki + 475);
    const auto *ki_476 = buffer.data(ki + 476);
    const auto *ki_477 = buffer.data(ki + 477);
    const auto *ki_478 = buffer.data(ki + 478);
    const auto *ki_479 = buffer.data(ki + 479);
    const auto *ki_480 = buffer.data(ki + 480);
    const auto *ki_481 = buffer.data(ki + 481);
    const auto *ki_482 = buffer.data(ki + 482);
    const auto *ki_483 = buffer.data(ki + 483);
    const auto *ki_484 = buffer.data(ki + 484);
    const auto *ki_485 = buffer.data(ki + 485);
    const auto *ki_486 = buffer.data(ki + 486);

    const auto *mi_663 = buffer.data(mi + 663);
    const auto *mi_664 = buffer.data(mi + 664);
    const auto *mi_665 = buffer.data(mi + 665);
    const auto *mi_666 = buffer.data(mi + 666);
    const auto *mi_667 = buffer.data(mi + 667);
    const auto *mi_668 = buffer.data(mi + 668);
    const auto *mi_669 = buffer.data(mi + 669);
    const auto *mi_670 = buffer.data(mi + 670);
    const auto *mi_671 = buffer.data(mi + 671);
    const auto *mi_672 = buffer.data(mi + 672);
    const auto *mi_673 = buffer.data(mi + 673);
    const auto *mi_674 = buffer.data(mi + 674);
    const auto *mi_675 = buffer.data(mi + 675);
    const auto *mi_676 = buffer.data(mi + 676);
    const auto *mi_677 = buffer.data(mi + 677);
    const auto *mi_678 = buffer.data(mi + 678);
    const auto *mi_679 = buffer.data(mi + 679);
    const auto *mi_680 = buffer.data(mi + 680);
    const auto *mi_681 = buffer.data(mi + 681);
    const auto *mi_682 = buffer.data(mi + 682);
    const auto *mi_683 = buffer.data(mi + 683);
    const auto *mi_684 = buffer.data(mi + 684);
    const auto *mi_685 = buffer.data(mi + 685);
    const auto *mi_686 = buffer.data(mi + 686);
    const auto *mi_687 = buffer.data(mi + 687);
    const auto *mi_688 = buffer.data(mi + 688);
    const auto *mi_689 = buffer.data(mi + 689);
    const auto *mi_690 = buffer.data(mi + 690);
    const auto *mi_691 = buffer.data(mi + 691);
    const auto *mi_692 = buffer.data(mi + 692);
    const auto *mi_693 = buffer.data(mi + 693);
    const auto *mi_694 = buffer.data(mi + 694);
    const auto *mi_695 = buffer.data(mi + 695);
    const auto *mi_696 = buffer.data(mi + 696);
    const auto *mi_697 = buffer.data(mi + 697);
    const auto *mi_698 = buffer.data(mi + 698);
    const auto *mi_699 = buffer.data(mi + 699);
    const auto *mi_700 = buffer.data(mi + 700);
    const auto *mi_701 = buffer.data(mi + 701);
    const auto *mi_702 = buffer.data(mi + 702);
    const auto *mi_703 = buffer.data(mi + 703);
    const auto *mi_704 = buffer.data(mi + 704);
    const auto *mi_705 = buffer.data(mi + 705);
    const auto *mi_706 = buffer.data(mi + 706);
    const auto *mi_707 = buffer.data(mi + 707);
    const auto *mi_708 = buffer.data(mi + 708);
    const auto *mi_709 = buffer.data(mi + 709);
    const auto *mi_710 = buffer.data(mi + 710);
    const auto *mi_711 = buffer.data(mi + 711);
    const auto *mi_712 = buffer.data(mi + 712);
    const auto *mi_713 = buffer.data(mi + 713);
    const auto *mi_714 = buffer.data(mi + 714);
    const auto *mi_715 = buffer.data(mi + 715);
    const auto *mi_716 = buffer.data(mi + 716);
    const auto *mi_717 = buffer.data(mi + 717);
    const auto *mi_718 = buffer.data(mi + 718);
    const auto *mi_719 = buffer.data(mi + 719);
    const auto *mi_720 = buffer.data(mi + 720);
    const auto *mi_721 = buffer.data(mi + 721);
    const auto *mi_722 = buffer.data(mi + 722);
    const auto *mi_723 = buffer.data(mi + 723);
    const auto *mi_724 = buffer.data(mi + 724);
    const auto *mi_725 = buffer.data(mi + 725);
    const auto *mi_726 = buffer.data(mi + 726);
    const auto *mi_727 = buffer.data(mi + 727);
    const auto *mi_728 = buffer.data(mi + 728);
    const auto *mi_729 = buffer.data(mi + 729);
    const auto *mi_730 = buffer.data(mi + 730);
    const auto *mi_731 = buffer.data(mi + 731);
    const auto *mi_732 = buffer.data(mi + 732);
    const auto *mi_733 = buffer.data(mi + 733);
    const auto *mi_734 = buffer.data(mi + 734);
    const auto *mi_735 = buffer.data(mi + 735);
    const auto *mi_736 = buffer.data(mi + 736);
    const auto *mi_737 = buffer.data(mi + 737);
    const auto *mi_738 = buffer.data(mi + 738);
    const auto *mi_739 = buffer.data(mi + 739);
    const auto *mi_740 = buffer.data(mi + 740);
    const auto *mi_741 = buffer.data(mi + 741);
    const auto *mi_742 = buffer.data(mi + 742);
    const auto *mi_743 = buffer.data(mi + 743);
    const auto *mi_744 = buffer.data(mi + 744);
    const auto *mi_745 = buffer.data(mi + 745);
    const auto *mi_746 = buffer.data(mi + 746);
    const auto *mi_747 = buffer.data(mi + 747);
    const auto *mi_748 = buffer.data(mi + 748);
    const auto *mi_749 = buffer.data(mi + 749);
    const auto *mi_750 = buffer.data(mi + 750);
    const auto *mi_751 = buffer.data(mi + 751);
    const auto *mi_752 = buffer.data(mi + 752);
    const auto *mi_753 = buffer.data(mi + 753);
    const auto *mi_754 = buffer.data(mi + 754);
    const auto *mi_755 = buffer.data(mi + 755);
    const auto *mi_784 = buffer.data(mi + 784);
    const auto *mi_785 = buffer.data(mi + 785);
    const auto *mi_786 = buffer.data(mi + 786);
    const auto *mi_787 = buffer.data(mi + 787);
    const auto *mi_788 = buffer.data(mi + 788);
    const auto *mi_789 = buffer.data(mi + 789);
    const auto *mi_790 = buffer.data(mi + 790);
    const auto *mi_791 = buffer.data(mi + 791);
    const auto *mi_792 = buffer.data(mi + 792);
    const auto *mi_793 = buffer.data(mi + 793);
    const auto *mi_794 = buffer.data(mi + 794);
    const auto *mi_795 = buffer.data(mi + 795);
    const auto *mi_796 = buffer.data(mi + 796);
    const auto *mi_797 = buffer.data(mi + 797);
    const auto *mi_798 = buffer.data(mi + 798);
    const auto *mi_799 = buffer.data(mi + 799);
    const auto *mi_800 = buffer.data(mi + 800);
    const auto *mi_801 = buffer.data(mi + 801);
    const auto *mi_802 = buffer.data(mi + 802);
    const auto *mi_803 = buffer.data(mi + 803);
    const auto *mi_804 = buffer.data(mi + 804);
    const auto *mi_805 = buffer.data(mi + 805);
    const auto *mi_806 = buffer.data(mi + 806);
    const auto *mi_807 = buffer.data(mi + 807);
    const auto *mi_808 = buffer.data(mi + 808);
    const auto *mi_809 = buffer.data(mi + 809);
    const auto *mi_810 = buffer.data(mi + 810);
    const auto *mi_811 = buffer.data(mi + 811);
    const auto *mi_812 = buffer.data(mi + 812);
    const auto *mi_813 = buffer.data(mi + 813);
    const auto *mi_814 = buffer.data(mi + 814);
    const auto *mi_815 = buffer.data(mi + 815);
    const auto *mi_816 = buffer.data(mi + 816);
    const auto *mi_817 = buffer.data(mi + 817);
    const auto *mi_818 = buffer.data(mi + 818);
    const auto *mi_819 = buffer.data(mi + 819);
    const auto *mi_820 = buffer.data(mi + 820);
    const auto *mi_821 = buffer.data(mi + 821);
    const auto *mi_822 = buffer.data(mi + 822);
    const auto *mi_823 = buffer.data(mi + 823);
    const auto *mi_824 = buffer.data(mi + 824);
    const auto *mi_825 = buffer.data(mi + 825);
    const auto *mi_826 = buffer.data(mi + 826);
    const auto *mi_827 = buffer.data(mi + 827);
    const auto *mi_828 = buffer.data(mi + 828);
    const auto *mi_829 = buffer.data(mi + 829);
    const auto *mi_830 = buffer.data(mi + 830);
    const auto *mi_831 = buffer.data(mi + 831);
    const auto *mi_832 = buffer.data(mi + 832);
    const auto *mi_833 = buffer.data(mi + 833);
    const auto *mi_834 = buffer.data(mi + 834);
    const auto *mi_835 = buffer.data(mi + 835);
    const auto *mi_836 = buffer.data(mi + 836);
    const auto *mi_837 = buffer.data(mi + 837);
    const auto *mi_838 = buffer.data(mi + 838);
    const auto *mi_839 = buffer.data(mi + 839);
    const auto *mi_840 = buffer.data(mi + 840);
    const auto *mi_841 = buffer.data(mi + 841);
    const auto *mi_842 = buffer.data(mi + 842);
    const auto *mi_843 = buffer.data(mi + 843);
    const auto *mi_844 = buffer.data(mi + 844);
    const auto *mi_845 = buffer.data(mi + 845);
    const auto *mi_846 = buffer.data(mi + 846);
    const auto *mi_847 = buffer.data(mi + 847);
    const auto *mi_848 = buffer.data(mi + 848);
    const auto *mi_849 = buffer.data(mi + 849);
    const auto *mi_850 = buffer.data(mi + 850);

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ki_355, ki_356, ki_357, ki_358, \
                         ki_359, mi_663, mi_664, mi_665, mi_666, \
                         mi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -3.0 * ki_355[k]
                   + f_0 * mi_663[k];

        t_496[k] = -3.0 * ki_356[k]
                   + f_0 * mi_664[k];

        t_497[k] = -3.0 * ki_357[k]
                   + f_0 * mi_665[k];

        t_498[k] = -3.0 * ki_358[k]
                   + f_0 * mi_666[k];

        t_499[k] = -3.0 * ki_359[k]
                   + f_0 * mi_667[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ki_360, ki_361, ki_362, ki_363, \
                         ki_364, mi_668, mi_669, mi_670, mi_671, \
                         mi_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -3.0 * ki_360[k]
                   + f_0 * mi_668[k];

        t_501[k] = -3.0 * ki_361[k]
                   + f_0 * mi_669[k];

        t_502[k] = -3.0 * ki_362[k]
                   + f_0 * mi_670[k];

        t_503[k] = -3.0 * ki_363[k]
                   + f_0 * mi_671[k];

        t_504[k] = -2.0 * ki_364[k]
                   + f_0 * mi_672[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ki_365, ki_366, ki_367, ki_368, \
                         ki_369, mi_673, mi_674, mi_675, mi_676, \
                         mi_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -2.0 * ki_365[k]
                   + f_0 * mi_673[k];

        t_506[k] = -2.0 * ki_366[k]
                   + f_0 * mi_674[k];

        t_507[k] = -2.0 * ki_367[k]
                   + f_0 * mi_675[k];

        t_508[k] = -2.0 * ki_368[k]
                   + f_0 * mi_676[k];

        t_509[k] = -2.0 * ki_369[k]
                   + f_0 * mi_677[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ki_370, ki_371, ki_372, ki_373, \
                         ki_374, mi_678, mi_679, mi_680, mi_681, \
                         mi_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -2.0 * ki_370[k]
                   + f_0 * mi_678[k];

        t_511[k] = -2.0 * ki_371[k]
                   + f_0 * mi_679[k];

        t_512[k] = -2.0 * ki_372[k]
                   + f_0 * mi_680[k];

        t_513[k] = -2.0 * ki_373[k]
                   + f_0 * mi_681[k];

        t_514[k] = -2.0 * ki_374[k]
                   + f_0 * mi_682[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ki_375, ki_376, ki_377, ki_378, \
                         ki_379, mi_683, mi_684, mi_685, mi_686, \
                         mi_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -2.0 * ki_375[k]
                   + f_0 * mi_683[k];

        t_516[k] = -2.0 * ki_376[k]
                   + f_0 * mi_684[k];

        t_517[k] = -2.0 * ki_377[k]
                   + f_0 * mi_685[k];

        t_518[k] = -2.0 * ki_378[k]
                   + f_0 * mi_686[k];

        t_519[k] = -2.0 * ki_379[k]
                   + f_0 * mi_687[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ki_380, ki_381, ki_382, ki_383, \
                         ki_384, mi_688, mi_689, mi_690, mi_691, \
                         mi_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -2.0 * ki_380[k]
                   + f_0 * mi_688[k];

        t_521[k] = -2.0 * ki_381[k]
                   + f_0 * mi_689[k];

        t_522[k] = -2.0 * ki_382[k]
                   + f_0 * mi_690[k];

        t_523[k] = -2.0 * ki_383[k]
                   + f_0 * mi_691[k];

        t_524[k] = -2.0 * ki_384[k]
                   + f_0 * mi_692[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ki_385, ki_386, ki_387, ki_388, \
                         ki_389, mi_693, mi_694, mi_695, mi_696, \
                         mi_697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -2.0 * ki_385[k]
                   + f_0 * mi_693[k];

        t_526[k] = -2.0 * ki_386[k]
                   + f_0 * mi_694[k];

        t_527[k] = -2.0 * ki_387[k]
                   + f_0 * mi_695[k];

        t_528[k] = -2.0 * ki_388[k]
                   + f_0 * mi_696[k];

        t_529[k] = -2.0 * ki_389[k]
                   + f_0 * mi_697[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ki_390, ki_391, ki_392, ki_393, \
                         ki_394, mi_698, mi_699, mi_700, mi_701, \
                         mi_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -2.0 * ki_390[k]
                   + f_0 * mi_698[k];

        t_531[k] = -2.0 * ki_391[k]
                   + f_0 * mi_699[k];

        t_532[k] = -ki_392[k]
                   + f_0 * mi_700[k];

        t_533[k] = -ki_393[k]
                   + f_0 * mi_701[k];

        t_534[k] = -ki_394[k]
                   + f_0 * mi_702[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ki_395, ki_396, ki_397, ki_398, \
                         ki_399, mi_703, mi_704, mi_705, mi_706, \
                         mi_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -ki_395[k]
                   + f_0 * mi_703[k];

        t_536[k] = -ki_396[k]
                   + f_0 * mi_704[k];

        t_537[k] = -ki_397[k]
                   + f_0 * mi_705[k];

        t_538[k] = -ki_398[k]
                   + f_0 * mi_706[k];

        t_539[k] = -ki_399[k]
                   + f_0 * mi_707[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ki_400, ki_401, ki_402, ki_403, \
                         ki_404, mi_708, mi_709, mi_710, mi_711, \
                         mi_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -ki_400[k]
                   + f_0 * mi_708[k];

        t_541[k] = -ki_401[k]
                   + f_0 * mi_709[k];

        t_542[k] = -ki_402[k]
                   + f_0 * mi_710[k];

        t_543[k] = -ki_403[k]
                   + f_0 * mi_711[k];

        t_544[k] = -ki_404[k]
                   + f_0 * mi_712[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ki_405, ki_406, ki_407, ki_408, \
                         ki_409, mi_713, mi_714, mi_715, mi_716, \
                         mi_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -ki_405[k]
                   + f_0 * mi_713[k];

        t_546[k] = -ki_406[k]
                   + f_0 * mi_714[k];

        t_547[k] = -ki_407[k]
                   + f_0 * mi_715[k];

        t_548[k] = -ki_408[k]
                   + f_0 * mi_716[k];

        t_549[k] = -ki_409[k]
                   + f_0 * mi_717[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ki_410, ki_411, ki_412, ki_413, \
                         ki_414, mi_718, mi_719, mi_720, mi_721, \
                         mi_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -ki_410[k]
                   + f_0 * mi_718[k];

        t_551[k] = -ki_411[k]
                   + f_0 * mi_719[k];

        t_552[k] = -ki_412[k]
                   + f_0 * mi_720[k];

        t_553[k] = -ki_413[k]
                   + f_0 * mi_721[k];

        t_554[k] = -ki_414[k]
                   + f_0 * mi_722[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ki_415, ki_416, ki_417, ki_418, \
                         ki_419, mi_723, mi_724, mi_725, mi_726, \
                         mi_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -ki_415[k]
                   + f_0 * mi_723[k];

        t_556[k] = -ki_416[k]
                   + f_0 * mi_724[k];

        t_557[k] = -ki_417[k]
                   + f_0 * mi_725[k];

        t_558[k] = -ki_418[k]
                   + f_0 * mi_726[k];

        t_559[k] = -ki_419[k]
                   + f_0 * mi_727[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, t_565, t_566, t_567, mi_728, \
                         mi_729, mi_730, mi_731, mi_732, mi_733, mi_734, \
                         mi_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_0 * mi_728[k];

        t_561[k] = f_0 * mi_729[k];

        t_562[k] = f_0 * mi_730[k];

        t_563[k] = f_0 * mi_731[k];

        t_564[k] = f_0 * mi_732[k];

        t_565[k] = f_0 * mi_733[k];

        t_566[k] = f_0 * mi_734[k];

        t_567[k] = f_0 * mi_735[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, t_572, t_573, t_574, t_575, mi_736, \
                         mi_737, mi_738, mi_739, mi_740, mi_741, mi_742, \
                         mi_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_0 * mi_736[k];

        t_569[k] = f_0 * mi_737[k];

        t_570[k] = f_0 * mi_738[k];

        t_571[k] = f_0 * mi_739[k];

        t_572[k] = f_0 * mi_740[k];

        t_573[k] = f_0 * mi_741[k];

        t_574[k] = f_0 * mi_742[k];

        t_575[k] = f_0 * mi_743[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, t_581, t_582, t_583, mi_744, \
                         mi_745, mi_746, mi_747, mi_748, mi_749, mi_750, \
                         mi_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_0 * mi_744[k];

        t_577[k] = f_0 * mi_745[k];

        t_578[k] = f_0 * mi_746[k];

        t_579[k] = f_0 * mi_747[k];

        t_580[k] = f_0 * mi_748[k];

        t_581[k] = f_0 * mi_749[k];

        t_582[k] = f_0 * mi_750[k];

        t_583[k] = f_0 * mi_751[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, t_588, t_589, ki_420, ki_421, mi_752, \
                         mi_753, mi_754, mi_755, mi_784, mi_785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_0 * mi_752[k];

        t_585[k] = f_0 * mi_753[k];

        t_586[k] = f_0 * mi_754[k];

        t_587[k] = f_0 * mi_755[k];

        t_588[k] = -6.0 * ki_420[k]
                   + f_0 * mi_784[k];

        t_589[k] = -6.0 * ki_421[k]
                   + f_0 * mi_785[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, ki_422, ki_423, ki_424, ki_425, \
                         ki_426, mi_786, mi_787, mi_788, mi_789, \
                         mi_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -6.0 * ki_422[k]
                   + f_0 * mi_786[k];

        t_591[k] = -6.0 * ki_423[k]
                   + f_0 * mi_787[k];

        t_592[k] = -6.0 * ki_424[k]
                   + f_0 * mi_788[k];

        t_593[k] = -6.0 * ki_425[k]
                   + f_0 * mi_789[k];

        t_594[k] = -6.0 * ki_426[k]
                   + f_0 * mi_790[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, ki_427, ki_428, ki_429, ki_430, \
                         ki_431, mi_791, mi_792, mi_793, mi_794, \
                         mi_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -6.0 * ki_427[k]
                   + f_0 * mi_791[k];

        t_596[k] = -6.0 * ki_428[k]
                   + f_0 * mi_792[k];

        t_597[k] = -6.0 * ki_429[k]
                   + f_0 * mi_793[k];

        t_598[k] = -6.0 * ki_430[k]
                   + f_0 * mi_794[k];

        t_599[k] = -6.0 * ki_431[k]
                   + f_0 * mi_795[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, ki_432, ki_433, ki_434, ki_435, \
                         ki_436, mi_796, mi_797, mi_798, mi_799, \
                         mi_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -6.0 * ki_432[k]
                   + f_0 * mi_796[k];

        t_601[k] = -6.0 * ki_433[k]
                   + f_0 * mi_797[k];

        t_602[k] = -6.0 * ki_434[k]
                   + f_0 * mi_798[k];

        t_603[k] = -6.0 * ki_435[k]
                   + f_0 * mi_799[k];

        t_604[k] = -6.0 * ki_436[k]
                   + f_0 * mi_800[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, ki_437, ki_438, ki_439, ki_440, \
                         ki_441, mi_801, mi_802, mi_803, mi_804, \
                         mi_805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -6.0 * ki_437[k]
                   + f_0 * mi_801[k];

        t_606[k] = -6.0 * ki_438[k]
                   + f_0 * mi_802[k];

        t_607[k] = -6.0 * ki_439[k]
                   + f_0 * mi_803[k];

        t_608[k] = -6.0 * ki_440[k]
                   + f_0 * mi_804[k];

        t_609[k] = -6.0 * ki_441[k]
                   + f_0 * mi_805[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, ki_442, ki_443, ki_444, ki_445, \
                         ki_446, mi_806, mi_807, mi_808, mi_809, \
                         mi_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -6.0 * ki_442[k]
                   + f_0 * mi_806[k];

        t_611[k] = -6.0 * ki_443[k]
                   + f_0 * mi_807[k];

        t_612[k] = -6.0 * ki_444[k]
                   + f_0 * mi_808[k];

        t_613[k] = -6.0 * ki_445[k]
                   + f_0 * mi_809[k];

        t_614[k] = -6.0 * ki_446[k]
                   + f_0 * mi_810[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, ki_447, ki_448, ki_449, ki_450, \
                         ki_451, mi_811, mi_812, mi_813, mi_814, \
                         mi_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -6.0 * ki_447[k]
                   + f_0 * mi_811[k];

        t_616[k] = -5.0 * ki_448[k]
                   + f_0 * mi_812[k];

        t_617[k] = -5.0 * ki_449[k]
                   + f_0 * mi_813[k];

        t_618[k] = -5.0 * ki_450[k]
                   + f_0 * mi_814[k];

        t_619[k] = -5.0 * ki_451[k]
                   + f_0 * mi_815[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, ki_452, ki_453, ki_454, ki_455, \
                         ki_456, mi_816, mi_817, mi_818, mi_819, \
                         mi_820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -5.0 * ki_452[k]
                   + f_0 * mi_816[k];

        t_621[k] = -5.0 * ki_453[k]
                   + f_0 * mi_817[k];

        t_622[k] = -5.0 * ki_454[k]
                   + f_0 * mi_818[k];

        t_623[k] = -5.0 * ki_455[k]
                   + f_0 * mi_819[k];

        t_624[k] = -5.0 * ki_456[k]
                   + f_0 * mi_820[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, ki_457, ki_458, ki_459, ki_460, \
                         ki_461, mi_821, mi_822, mi_823, mi_824, \
                         mi_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -5.0 * ki_457[k]
                   + f_0 * mi_821[k];

        t_626[k] = -5.0 * ki_458[k]
                   + f_0 * mi_822[k];

        t_627[k] = -5.0 * ki_459[k]
                   + f_0 * mi_823[k];

        t_628[k] = -5.0 * ki_460[k]
                   + f_0 * mi_824[k];

        t_629[k] = -5.0 * ki_461[k]
                   + f_0 * mi_825[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, ki_462, ki_463, ki_464, ki_465, \
                         ki_466, mi_826, mi_827, mi_828, mi_829, \
                         mi_830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -5.0 * ki_462[k]
                   + f_0 * mi_826[k];

        t_631[k] = -5.0 * ki_463[k]
                   + f_0 * mi_827[k];

        t_632[k] = -5.0 * ki_464[k]
                   + f_0 * mi_828[k];

        t_633[k] = -5.0 * ki_465[k]
                   + f_0 * mi_829[k];

        t_634[k] = -5.0 * ki_466[k]
                   + f_0 * mi_830[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, ki_467, ki_468, ki_469, ki_470, \
                         ki_471, mi_831, mi_832, mi_833, mi_834, \
                         mi_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -5.0 * ki_467[k]
                   + f_0 * mi_831[k];

        t_636[k] = -5.0 * ki_468[k]
                   + f_0 * mi_832[k];

        t_637[k] = -5.0 * ki_469[k]
                   + f_0 * mi_833[k];

        t_638[k] = -5.0 * ki_470[k]
                   + f_0 * mi_834[k];

        t_639[k] = -5.0 * ki_471[k]
                   + f_0 * mi_835[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, ki_472, ki_473, ki_474, ki_475, \
                         ki_476, mi_836, mi_837, mi_838, mi_839, \
                         mi_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -5.0 * ki_472[k]
                   + f_0 * mi_836[k];

        t_641[k] = -5.0 * ki_473[k]
                   + f_0 * mi_837[k];

        t_642[k] = -5.0 * ki_474[k]
                   + f_0 * mi_838[k];

        t_643[k] = -5.0 * ki_475[k]
                   + f_0 * mi_839[k];

        t_644[k] = -4.0 * ki_476[k]
                   + f_0 * mi_840[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, ki_477, ki_478, ki_479, ki_480, \
                         ki_481, mi_841, mi_842, mi_843, mi_844, \
                         mi_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -4.0 * ki_477[k]
                   + f_0 * mi_841[k];

        t_646[k] = -4.0 * ki_478[k]
                   + f_0 * mi_842[k];

        t_647[k] = -4.0 * ki_479[k]
                   + f_0 * mi_843[k];

        t_648[k] = -4.0 * ki_480[k]
                   + f_0 * mi_844[k];

        t_649[k] = -4.0 * ki_481[k]
                   + f_0 * mi_845[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, ki_482, ki_483, ki_484, ki_485, \
                         ki_486, mi_846, mi_847, mi_848, mi_849, \
                         mi_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -4.0 * ki_482[k]
                   + f_0 * mi_846[k];

        t_651[k] = -4.0 * ki_483[k]
                   + f_0 * mi_847[k];

        t_652[k] = -4.0 * ki_484[k]
                   + f_0 * mi_848[k];

        t_653[k] = -4.0 * ki_485[k]
                   + f_0 * mi_849[k];

        t_654[k] = -4.0 * ki_486[k]
                   + f_0 * mi_850[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);

    const auto *ki_487 = buffer.data(ki + 487);
    const auto *ki_488 = buffer.data(ki + 488);
    const auto *ki_489 = buffer.data(ki + 489);
    const auto *ki_490 = buffer.data(ki + 490);
    const auto *ki_491 = buffer.data(ki + 491);
    const auto *ki_492 = buffer.data(ki + 492);
    const auto *ki_493 = buffer.data(ki + 493);
    const auto *ki_494 = buffer.data(ki + 494);
    const auto *ki_495 = buffer.data(ki + 495);
    const auto *ki_496 = buffer.data(ki + 496);
    const auto *ki_497 = buffer.data(ki + 497);
    const auto *ki_498 = buffer.data(ki + 498);
    const auto *ki_499 = buffer.data(ki + 499);
    const auto *ki_500 = buffer.data(ki + 500);
    const auto *ki_501 = buffer.data(ki + 501);
    const auto *ki_502 = buffer.data(ki + 502);
    const auto *ki_503 = buffer.data(ki + 503);
    const auto *ki_504 = buffer.data(ki + 504);
    const auto *ki_505 = buffer.data(ki + 505);
    const auto *ki_506 = buffer.data(ki + 506);
    const auto *ki_507 = buffer.data(ki + 507);
    const auto *ki_508 = buffer.data(ki + 508);
    const auto *ki_509 = buffer.data(ki + 509);
    const auto *ki_510 = buffer.data(ki + 510);
    const auto *ki_511 = buffer.data(ki + 511);
    const auto *ki_512 = buffer.data(ki + 512);
    const auto *ki_513 = buffer.data(ki + 513);
    const auto *ki_514 = buffer.data(ki + 514);
    const auto *ki_515 = buffer.data(ki + 515);
    const auto *ki_516 = buffer.data(ki + 516);
    const auto *ki_517 = buffer.data(ki + 517);
    const auto *ki_518 = buffer.data(ki + 518);
    const auto *ki_519 = buffer.data(ki + 519);
    const auto *ki_520 = buffer.data(ki + 520);
    const auto *ki_521 = buffer.data(ki + 521);
    const auto *ki_522 = buffer.data(ki + 522);
    const auto *ki_523 = buffer.data(ki + 523);
    const auto *ki_524 = buffer.data(ki + 524);
    const auto *ki_525 = buffer.data(ki + 525);
    const auto *ki_526 = buffer.data(ki + 526);
    const auto *ki_527 = buffer.data(ki + 527);
    const auto *ki_528 = buffer.data(ki + 528);
    const auto *ki_529 = buffer.data(ki + 529);
    const auto *ki_530 = buffer.data(ki + 530);
    const auto *ki_531 = buffer.data(ki + 531);
    const auto *ki_532 = buffer.data(ki + 532);
    const auto *ki_533 = buffer.data(ki + 533);
    const auto *ki_534 = buffer.data(ki + 534);
    const auto *ki_535 = buffer.data(ki + 535);
    const auto *ki_536 = buffer.data(ki + 536);
    const auto *ki_537 = buffer.data(ki + 537);
    const auto *ki_538 = buffer.data(ki + 538);
    const auto *ki_539 = buffer.data(ki + 539);
    const auto *ki_540 = buffer.data(ki + 540);
    const auto *ki_541 = buffer.data(ki + 541);
    const auto *ki_542 = buffer.data(ki + 542);
    const auto *ki_543 = buffer.data(ki + 543);
    const auto *ki_544 = buffer.data(ki + 544);
    const auto *ki_545 = buffer.data(ki + 545);
    const auto *ki_546 = buffer.data(ki + 546);
    const auto *ki_547 = buffer.data(ki + 547);
    const auto *ki_548 = buffer.data(ki + 548);
    const auto *ki_549 = buffer.data(ki + 549);
    const auto *ki_550 = buffer.data(ki + 550);
    const auto *ki_551 = buffer.data(ki + 551);
    const auto *ki_552 = buffer.data(ki + 552);
    const auto *ki_553 = buffer.data(ki + 553);
    const auto *ki_554 = buffer.data(ki + 554);
    const auto *ki_555 = buffer.data(ki + 555);
    const auto *ki_556 = buffer.data(ki + 556);
    const auto *ki_557 = buffer.data(ki + 557);
    const auto *ki_558 = buffer.data(ki + 558);
    const auto *ki_559 = buffer.data(ki + 559);
    const auto *ki_560 = buffer.data(ki + 560);
    const auto *ki_561 = buffer.data(ki + 561);
    const auto *ki_562 = buffer.data(ki + 562);
    const auto *ki_563 = buffer.data(ki + 563);
    const auto *ki_564 = buffer.data(ki + 564);
    const auto *ki_565 = buffer.data(ki + 565);
    const auto *ki_566 = buffer.data(ki + 566);
    const auto *ki_567 = buffer.data(ki + 567);
    const auto *ki_568 = buffer.data(ki + 568);
    const auto *ki_569 = buffer.data(ki + 569);
    const auto *ki_570 = buffer.data(ki + 570);
    const auto *ki_571 = buffer.data(ki + 571);
    const auto *ki_572 = buffer.data(ki + 572);
    const auto *ki_573 = buffer.data(ki + 573);
    const auto *ki_574 = buffer.data(ki + 574);
    const auto *ki_575 = buffer.data(ki + 575);
    const auto *ki_576 = buffer.data(ki + 576);
    const auto *ki_577 = buffer.data(ki + 577);
    const auto *ki_578 = buffer.data(ki + 578);
    const auto *ki_579 = buffer.data(ki + 579);
    const auto *ki_580 = buffer.data(ki + 580);
    const auto *ki_581 = buffer.data(ki + 581);
    const auto *ki_582 = buffer.data(ki + 582);
    const auto *ki_583 = buffer.data(ki + 583);
    const auto *ki_584 = buffer.data(ki + 584);
    const auto *ki_585 = buffer.data(ki + 585);
    const auto *ki_586 = buffer.data(ki + 586);
    const auto *ki_587 = buffer.data(ki + 587);
    const auto *ki_588 = buffer.data(ki + 588);
    const auto *ki_589 = buffer.data(ki + 589);
    const auto *ki_590 = buffer.data(ki + 590);
    const auto *ki_591 = buffer.data(ki + 591);
    const auto *ki_592 = buffer.data(ki + 592);
    const auto *ki_593 = buffer.data(ki + 593);
    const auto *ki_594 = buffer.data(ki + 594);
    const auto *ki_595 = buffer.data(ki + 595);
    const auto *ki_596 = buffer.data(ki + 596);
    const auto *ki_597 = buffer.data(ki + 597);
    const auto *ki_598 = buffer.data(ki + 598);
    const auto *ki_599 = buffer.data(ki + 599);
    const auto *ki_600 = buffer.data(ki + 600);
    const auto *ki_601 = buffer.data(ki + 601);
    const auto *ki_602 = buffer.data(ki + 602);
    const auto *ki_603 = buffer.data(ki + 603);
    const auto *ki_604 = buffer.data(ki + 604);
    const auto *ki_605 = buffer.data(ki + 605);
    const auto *ki_606 = buffer.data(ki + 606);
    const auto *ki_607 = buffer.data(ki + 607);
    const auto *ki_608 = buffer.data(ki + 608);
    const auto *ki_609 = buffer.data(ki + 609);
    const auto *ki_610 = buffer.data(ki + 610);
    const auto *ki_611 = buffer.data(ki + 611);
    const auto *ki_612 = buffer.data(ki + 612);
    const auto *ki_613 = buffer.data(ki + 613);
    const auto *ki_614 = buffer.data(ki + 614);
    const auto *ki_615 = buffer.data(ki + 615);
    const auto *ki_616 = buffer.data(ki + 616);
    const auto *ki_617 = buffer.data(ki + 617);
    const auto *ki_618 = buffer.data(ki + 618);

    const auto *mi_851 = buffer.data(mi + 851);
    const auto *mi_852 = buffer.data(mi + 852);
    const auto *mi_853 = buffer.data(mi + 853);
    const auto *mi_854 = buffer.data(mi + 854);
    const auto *mi_855 = buffer.data(mi + 855);
    const auto *mi_856 = buffer.data(mi + 856);
    const auto *mi_857 = buffer.data(mi + 857);
    const auto *mi_858 = buffer.data(mi + 858);
    const auto *mi_859 = buffer.data(mi + 859);
    const auto *mi_860 = buffer.data(mi + 860);
    const auto *mi_861 = buffer.data(mi + 861);
    const auto *mi_862 = buffer.data(mi + 862);
    const auto *mi_863 = buffer.data(mi + 863);
    const auto *mi_864 = buffer.data(mi + 864);
    const auto *mi_865 = buffer.data(mi + 865);
    const auto *mi_866 = buffer.data(mi + 866);
    const auto *mi_867 = buffer.data(mi + 867);
    const auto *mi_868 = buffer.data(mi + 868);
    const auto *mi_869 = buffer.data(mi + 869);
    const auto *mi_870 = buffer.data(mi + 870);
    const auto *mi_871 = buffer.data(mi + 871);
    const auto *mi_872 = buffer.data(mi + 872);
    const auto *mi_873 = buffer.data(mi + 873);
    const auto *mi_874 = buffer.data(mi + 874);
    const auto *mi_875 = buffer.data(mi + 875);
    const auto *mi_876 = buffer.data(mi + 876);
    const auto *mi_877 = buffer.data(mi + 877);
    const auto *mi_878 = buffer.data(mi + 878);
    const auto *mi_879 = buffer.data(mi + 879);
    const auto *mi_880 = buffer.data(mi + 880);
    const auto *mi_881 = buffer.data(mi + 881);
    const auto *mi_882 = buffer.data(mi + 882);
    const auto *mi_883 = buffer.data(mi + 883);
    const auto *mi_884 = buffer.data(mi + 884);
    const auto *mi_885 = buffer.data(mi + 885);
    const auto *mi_886 = buffer.data(mi + 886);
    const auto *mi_887 = buffer.data(mi + 887);
    const auto *mi_888 = buffer.data(mi + 888);
    const auto *mi_889 = buffer.data(mi + 889);
    const auto *mi_890 = buffer.data(mi + 890);
    const auto *mi_891 = buffer.data(mi + 891);
    const auto *mi_892 = buffer.data(mi + 892);
    const auto *mi_893 = buffer.data(mi + 893);
    const auto *mi_894 = buffer.data(mi + 894);
    const auto *mi_895 = buffer.data(mi + 895);
    const auto *mi_896 = buffer.data(mi + 896);
    const auto *mi_897 = buffer.data(mi + 897);
    const auto *mi_898 = buffer.data(mi + 898);
    const auto *mi_899 = buffer.data(mi + 899);
    const auto *mi_900 = buffer.data(mi + 900);
    const auto *mi_901 = buffer.data(mi + 901);
    const auto *mi_902 = buffer.data(mi + 902);
    const auto *mi_903 = buffer.data(mi + 903);
    const auto *mi_904 = buffer.data(mi + 904);
    const auto *mi_905 = buffer.data(mi + 905);
    const auto *mi_906 = buffer.data(mi + 906);
    const auto *mi_907 = buffer.data(mi + 907);
    const auto *mi_908 = buffer.data(mi + 908);
    const auto *mi_909 = buffer.data(mi + 909);
    const auto *mi_910 = buffer.data(mi + 910);
    const auto *mi_911 = buffer.data(mi + 911);
    const auto *mi_912 = buffer.data(mi + 912);
    const auto *mi_913 = buffer.data(mi + 913);
    const auto *mi_914 = buffer.data(mi + 914);
    const auto *mi_915 = buffer.data(mi + 915);
    const auto *mi_916 = buffer.data(mi + 916);
    const auto *mi_917 = buffer.data(mi + 917);
    const auto *mi_918 = buffer.data(mi + 918);
    const auto *mi_919 = buffer.data(mi + 919);
    const auto *mi_920 = buffer.data(mi + 920);
    const auto *mi_921 = buffer.data(mi + 921);
    const auto *mi_922 = buffer.data(mi + 922);
    const auto *mi_923 = buffer.data(mi + 923);
    const auto *mi_924 = buffer.data(mi + 924);
    const auto *mi_925 = buffer.data(mi + 925);
    const auto *mi_926 = buffer.data(mi + 926);
    const auto *mi_927 = buffer.data(mi + 927);
    const auto *mi_928 = buffer.data(mi + 928);
    const auto *mi_929 = buffer.data(mi + 929);
    const auto *mi_930 = buffer.data(mi + 930);
    const auto *mi_931 = buffer.data(mi + 931);
    const auto *mi_932 = buffer.data(mi + 932);
    const auto *mi_933 = buffer.data(mi + 933);
    const auto *mi_934 = buffer.data(mi + 934);
    const auto *mi_935 = buffer.data(mi + 935);
    const auto *mi_936 = buffer.data(mi + 936);
    const auto *mi_937 = buffer.data(mi + 937);
    const auto *mi_938 = buffer.data(mi + 938);
    const auto *mi_939 = buffer.data(mi + 939);
    const auto *mi_940 = buffer.data(mi + 940);
    const auto *mi_941 = buffer.data(mi + 941);
    const auto *mi_942 = buffer.data(mi + 942);
    const auto *mi_943 = buffer.data(mi + 943);
    const auto *mi_944 = buffer.data(mi + 944);
    const auto *mi_945 = buffer.data(mi + 945);
    const auto *mi_946 = buffer.data(mi + 946);
    const auto *mi_947 = buffer.data(mi + 947);
    const auto *mi_948 = buffer.data(mi + 948);
    const auto *mi_949 = buffer.data(mi + 949);
    const auto *mi_950 = buffer.data(mi + 950);
    const auto *mi_951 = buffer.data(mi + 951);
    const auto *mi_952 = buffer.data(mi + 952);
    const auto *mi_953 = buffer.data(mi + 953);
    const auto *mi_954 = buffer.data(mi + 954);
    const auto *mi_955 = buffer.data(mi + 955);
    const auto *mi_956 = buffer.data(mi + 956);
    const auto *mi_957 = buffer.data(mi + 957);
    const auto *mi_958 = buffer.data(mi + 958);
    const auto *mi_959 = buffer.data(mi + 959);
    const auto *mi_960 = buffer.data(mi + 960);
    const auto *mi_961 = buffer.data(mi + 961);
    const auto *mi_962 = buffer.data(mi + 962);
    const auto *mi_963 = buffer.data(mi + 963);
    const auto *mi_964 = buffer.data(mi + 964);
    const auto *mi_965 = buffer.data(mi + 965);
    const auto *mi_966 = buffer.data(mi + 966);
    const auto *mi_967 = buffer.data(mi + 967);
    const auto *mi_968 = buffer.data(mi + 968);
    const auto *mi_969 = buffer.data(mi + 969);
    const auto *mi_970 = buffer.data(mi + 970);
    const auto *mi_971 = buffer.data(mi + 971);
    const auto *mi_972 = buffer.data(mi + 972);
    const auto *mi_973 = buffer.data(mi + 973);
    const auto *mi_974 = buffer.data(mi + 974);
    const auto *mi_975 = buffer.data(mi + 975);
    const auto *mi_976 = buffer.data(mi + 976);
    const auto *mi_977 = buffer.data(mi + 977);
    const auto *mi_978 = buffer.data(mi + 978);
    const auto *mi_979 = buffer.data(mi + 979);
    const auto *mi_1008 = buffer.data(mi + 1008);
    const auto *mi_1009 = buffer.data(mi + 1009);
    const auto *mi_1010 = buffer.data(mi + 1010);
    const auto *mi_1011 = buffer.data(mi + 1011);
    const auto *mi_1012 = buffer.data(mi + 1012);
    const auto *mi_1013 = buffer.data(mi + 1013);
    const auto *mi_1014 = buffer.data(mi + 1014);
    const auto *mi_1015 = buffer.data(mi + 1015);
    const auto *mi_1016 = buffer.data(mi + 1016);
    const auto *mi_1017 = buffer.data(mi + 1017);
    const auto *mi_1018 = buffer.data(mi + 1018);
    const auto *mi_1019 = buffer.data(mi + 1019);
    const auto *mi_1020 = buffer.data(mi + 1020);
    const auto *mi_1021 = buffer.data(mi + 1021);
    const auto *mi_1022 = buffer.data(mi + 1022);
    const auto *mi_1023 = buffer.data(mi + 1023);
    const auto *mi_1024 = buffer.data(mi + 1024);
    const auto *mi_1025 = buffer.data(mi + 1025);
    const auto *mi_1026 = buffer.data(mi + 1026);
    const auto *mi_1027 = buffer.data(mi + 1027);
    const auto *mi_1028 = buffer.data(mi + 1028);
    const auto *mi_1029 = buffer.data(mi + 1029);
    const auto *mi_1030 = buffer.data(mi + 1030);
    const auto *mi_1031 = buffer.data(mi + 1031);
    const auto *mi_1032 = buffer.data(mi + 1032);
    const auto *mi_1033 = buffer.data(mi + 1033);
    const auto *mi_1034 = buffer.data(mi + 1034);
    const auto *mi_1035 = buffer.data(mi + 1035);
    const auto *mi_1036 = buffer.data(mi + 1036);
    const auto *mi_1037 = buffer.data(mi + 1037);
    const auto *mi_1038 = buffer.data(mi + 1038);

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, ki_487, ki_488, ki_489, ki_490, \
                         ki_491, mi_851, mi_852, mi_853, mi_854, \
                         mi_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -4.0 * ki_487[k]
                   + f_0 * mi_851[k];

        t_656[k] = -4.0 * ki_488[k]
                   + f_0 * mi_852[k];

        t_657[k] = -4.0 * ki_489[k]
                   + f_0 * mi_853[k];

        t_658[k] = -4.0 * ki_490[k]
                   + f_0 * mi_854[k];

        t_659[k] = -4.0 * ki_491[k]
                   + f_0 * mi_855[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, ki_492, ki_493, ki_494, ki_495, \
                         ki_496, mi_856, mi_857, mi_858, mi_859, \
                         mi_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -4.0 * ki_492[k]
                   + f_0 * mi_856[k];

        t_661[k] = -4.0 * ki_493[k]
                   + f_0 * mi_857[k];

        t_662[k] = -4.0 * ki_494[k]
                   + f_0 * mi_858[k];

        t_663[k] = -4.0 * ki_495[k]
                   + f_0 * mi_859[k];

        t_664[k] = -4.0 * ki_496[k]
                   + f_0 * mi_860[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, ki_497, ki_498, ki_499, ki_500, \
                         ki_501, mi_861, mi_862, mi_863, mi_864, \
                         mi_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -4.0 * ki_497[k]
                   + f_0 * mi_861[k];

        t_666[k] = -4.0 * ki_498[k]
                   + f_0 * mi_862[k];

        t_667[k] = -4.0 * ki_499[k]
                   + f_0 * mi_863[k];

        t_668[k] = -4.0 * ki_500[k]
                   + f_0 * mi_864[k];

        t_669[k] = -4.0 * ki_501[k]
                   + f_0 * mi_865[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, ki_502, ki_503, ki_504, ki_505, \
                         ki_506, mi_866, mi_867, mi_868, mi_869, \
                         mi_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -4.0 * ki_502[k]
                   + f_0 * mi_866[k];

        t_671[k] = -4.0 * ki_503[k]
                   + f_0 * mi_867[k];

        t_672[k] = -3.0 * ki_504[k]
                   + f_0 * mi_868[k];

        t_673[k] = -3.0 * ki_505[k]
                   + f_0 * mi_869[k];

        t_674[k] = -3.0 * ki_506[k]
                   + f_0 * mi_870[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, ki_507, ki_508, ki_509, ki_510, \
                         ki_511, mi_871, mi_872, mi_873, mi_874, \
                         mi_875 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -3.0 * ki_507[k]
                   + f_0 * mi_871[k];

        t_676[k] = -3.0 * ki_508[k]
                   + f_0 * mi_872[k];

        t_677[k] = -3.0 * ki_509[k]
                   + f_0 * mi_873[k];

        t_678[k] = -3.0 * ki_510[k]
                   + f_0 * mi_874[k];

        t_679[k] = -3.0 * ki_511[k]
                   + f_0 * mi_875[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, ki_512, ki_513, ki_514, ki_515, \
                         ki_516, mi_876, mi_877, mi_878, mi_879, \
                         mi_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -3.0 * ki_512[k]
                   + f_0 * mi_876[k];

        t_681[k] = -3.0 * ki_513[k]
                   + f_0 * mi_877[k];

        t_682[k] = -3.0 * ki_514[k]
                   + f_0 * mi_878[k];

        t_683[k] = -3.0 * ki_515[k]
                   + f_0 * mi_879[k];

        t_684[k] = -3.0 * ki_516[k]
                   + f_0 * mi_880[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, ki_517, ki_518, ki_519, ki_520, \
                         ki_521, mi_881, mi_882, mi_883, mi_884, \
                         mi_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -3.0 * ki_517[k]
                   + f_0 * mi_881[k];

        t_686[k] = -3.0 * ki_518[k]
                   + f_0 * mi_882[k];

        t_687[k] = -3.0 * ki_519[k]
                   + f_0 * mi_883[k];

        t_688[k] = -3.0 * ki_520[k]
                   + f_0 * mi_884[k];

        t_689[k] = -3.0 * ki_521[k]
                   + f_0 * mi_885[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, ki_522, ki_523, ki_524, ki_525, \
                         ki_526, mi_886, mi_887, mi_888, mi_889, \
                         mi_890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -3.0 * ki_522[k]
                   + f_0 * mi_886[k];

        t_691[k] = -3.0 * ki_523[k]
                   + f_0 * mi_887[k];

        t_692[k] = -3.0 * ki_524[k]
                   + f_0 * mi_888[k];

        t_693[k] = -3.0 * ki_525[k]
                   + f_0 * mi_889[k];

        t_694[k] = -3.0 * ki_526[k]
                   + f_0 * mi_890[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, ki_527, ki_528, ki_529, ki_530, \
                         ki_531, mi_891, mi_892, mi_893, mi_894, \
                         mi_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -3.0 * ki_527[k]
                   + f_0 * mi_891[k];

        t_696[k] = -3.0 * ki_528[k]
                   + f_0 * mi_892[k];

        t_697[k] = -3.0 * ki_529[k]
                   + f_0 * mi_893[k];

        t_698[k] = -3.0 * ki_530[k]
                   + f_0 * mi_894[k];

        t_699[k] = -3.0 * ki_531[k]
                   + f_0 * mi_895[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, ki_532, ki_533, ki_534, ki_535, \
                         ki_536, mi_896, mi_897, mi_898, mi_899, \
                         mi_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -2.0 * ki_532[k]
                   + f_0 * mi_896[k];

        t_701[k] = -2.0 * ki_533[k]
                   + f_0 * mi_897[k];

        t_702[k] = -2.0 * ki_534[k]
                   + f_0 * mi_898[k];

        t_703[k] = -2.0 * ki_535[k]
                   + f_0 * mi_899[k];

        t_704[k] = -2.0 * ki_536[k]
                   + f_0 * mi_900[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, ki_537, ki_538, ki_539, ki_540, \
                         ki_541, mi_901, mi_902, mi_903, mi_904, \
                         mi_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -2.0 * ki_537[k]
                   + f_0 * mi_901[k];

        t_706[k] = -2.0 * ki_538[k]
                   + f_0 * mi_902[k];

        t_707[k] = -2.0 * ki_539[k]
                   + f_0 * mi_903[k];

        t_708[k] = -2.0 * ki_540[k]
                   + f_0 * mi_904[k];

        t_709[k] = -2.0 * ki_541[k]
                   + f_0 * mi_905[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, ki_542, ki_543, ki_544, ki_545, \
                         ki_546, mi_906, mi_907, mi_908, mi_909, \
                         mi_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -2.0 * ki_542[k]
                   + f_0 * mi_906[k];

        t_711[k] = -2.0 * ki_543[k]
                   + f_0 * mi_907[k];

        t_712[k] = -2.0 * ki_544[k]
                   + f_0 * mi_908[k];

        t_713[k] = -2.0 * ki_545[k]
                   + f_0 * mi_909[k];

        t_714[k] = -2.0 * ki_546[k]
                   + f_0 * mi_910[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, ki_547, ki_548, ki_549, ki_550, \
                         ki_551, mi_911, mi_912, mi_913, mi_914, \
                         mi_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -2.0 * ki_547[k]
                   + f_0 * mi_911[k];

        t_716[k] = -2.0 * ki_548[k]
                   + f_0 * mi_912[k];

        t_717[k] = -2.0 * ki_549[k]
                   + f_0 * mi_913[k];

        t_718[k] = -2.0 * ki_550[k]
                   + f_0 * mi_914[k];

        t_719[k] = -2.0 * ki_551[k]
                   + f_0 * mi_915[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, ki_552, ki_553, ki_554, ki_555, \
                         ki_556, mi_916, mi_917, mi_918, mi_919, \
                         mi_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -2.0 * ki_552[k]
                   + f_0 * mi_916[k];

        t_721[k] = -2.0 * ki_553[k]
                   + f_0 * mi_917[k];

        t_722[k] = -2.0 * ki_554[k]
                   + f_0 * mi_918[k];

        t_723[k] = -2.0 * ki_555[k]
                   + f_0 * mi_919[k];

        t_724[k] = -2.0 * ki_556[k]
                   + f_0 * mi_920[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, ki_557, ki_558, ki_559, ki_560, \
                         ki_561, mi_921, mi_922, mi_923, mi_924, \
                         mi_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -2.0 * ki_557[k]
                   + f_0 * mi_921[k];

        t_726[k] = -2.0 * ki_558[k]
                   + f_0 * mi_922[k];

        t_727[k] = -2.0 * ki_559[k]
                   + f_0 * mi_923[k];

        t_728[k] = -ki_560[k]
                   + f_0 * mi_924[k];

        t_729[k] = -ki_561[k]
                   + f_0 * mi_925[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, ki_562, ki_563, ki_564, ki_565, \
                         ki_566, mi_926, mi_927, mi_928, mi_929, \
                         mi_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -ki_562[k]
                   + f_0 * mi_926[k];

        t_731[k] = -ki_563[k]
                   + f_0 * mi_927[k];

        t_732[k] = -ki_564[k]
                   + f_0 * mi_928[k];

        t_733[k] = -ki_565[k]
                   + f_0 * mi_929[k];

        t_734[k] = -ki_566[k]
                   + f_0 * mi_930[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, ki_567, ki_568, ki_569, ki_570, \
                         ki_571, mi_931, mi_932, mi_933, mi_934, \
                         mi_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -ki_567[k]
                   + f_0 * mi_931[k];

        t_736[k] = -ki_568[k]
                   + f_0 * mi_932[k];

        t_737[k] = -ki_569[k]
                   + f_0 * mi_933[k];

        t_738[k] = -ki_570[k]
                   + f_0 * mi_934[k];

        t_739[k] = -ki_571[k]
                   + f_0 * mi_935[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, ki_572, ki_573, ki_574, ki_575, \
                         ki_576, mi_936, mi_937, mi_938, mi_939, \
                         mi_940 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -ki_572[k]
                   + f_0 * mi_936[k];

        t_741[k] = -ki_573[k]
                   + f_0 * mi_937[k];

        t_742[k] = -ki_574[k]
                   + f_0 * mi_938[k];

        t_743[k] = -ki_575[k]
                   + f_0 * mi_939[k];

        t_744[k] = -ki_576[k]
                   + f_0 * mi_940[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, ki_577, ki_578, ki_579, ki_580, \
                         ki_581, mi_941, mi_942, mi_943, mi_944, \
                         mi_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -ki_577[k]
                   + f_0 * mi_941[k];

        t_746[k] = -ki_578[k]
                   + f_0 * mi_942[k];

        t_747[k] = -ki_579[k]
                   + f_0 * mi_943[k];

        t_748[k] = -ki_580[k]
                   + f_0 * mi_944[k];

        t_749[k] = -ki_581[k]
                   + f_0 * mi_945[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, ki_582, ki_583, ki_584, ki_585, \
                         ki_586, mi_946, mi_947, mi_948, mi_949, \
                         mi_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -ki_582[k]
                   + f_0 * mi_946[k];

        t_751[k] = -ki_583[k]
                   + f_0 * mi_947[k];

        t_752[k] = -ki_584[k]
                   + f_0 * mi_948[k];

        t_753[k] = -ki_585[k]
                   + f_0 * mi_949[k];

        t_754[k] = -ki_586[k]
                   + f_0 * mi_950[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, t_760, t_761, ki_587, mi_951, \
                         mi_952, mi_953, mi_954, mi_955, mi_956, \
                         mi_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -ki_587[k]
                   + f_0 * mi_951[k];

        t_756[k] = f_0 * mi_952[k];

        t_757[k] = f_0 * mi_953[k];

        t_758[k] = f_0 * mi_954[k];

        t_759[k] = f_0 * mi_955[k];

        t_760[k] = f_0 * mi_956[k];

        t_761[k] = f_0 * mi_957[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, t_766, t_767, t_768, t_769, mi_958, \
                         mi_959, mi_960, mi_961, mi_962, mi_963, mi_964, \
                         mi_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_0 * mi_958[k];

        t_763[k] = f_0 * mi_959[k];

        t_764[k] = f_0 * mi_960[k];

        t_765[k] = f_0 * mi_961[k];

        t_766[k] = f_0 * mi_962[k];

        t_767[k] = f_0 * mi_963[k];

        t_768[k] = f_0 * mi_964[k];

        t_769[k] = f_0 * mi_965[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, t_775, t_776, t_777, mi_966, \
                         mi_967, mi_968, mi_969, mi_970, mi_971, mi_972, \
                         mi_973 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_0 * mi_966[k];

        t_771[k] = f_0 * mi_967[k];

        t_772[k] = f_0 * mi_968[k];

        t_773[k] = f_0 * mi_969[k];

        t_774[k] = f_0 * mi_970[k];

        t_775[k] = f_0 * mi_971[k];

        t_776[k] = f_0 * mi_972[k];

        t_777[k] = f_0 * mi_973[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, t_783, t_784, ki_588, mi_974, \
                         mi_975, mi_976, mi_977, mi_978, mi_979, \
                         mi_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_0 * mi_974[k];

        t_779[k] = f_0 * mi_975[k];

        t_780[k] = f_0 * mi_976[k];

        t_781[k] = f_0 * mi_977[k];

        t_782[k] = f_0 * mi_978[k];

        t_783[k] = f_0 * mi_979[k];

        t_784[k] = -7.0 * ki_588[k]
                   + f_0 * mi_1008[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, t_789, ki_589, ki_590, ki_591, ki_592, \
                         ki_593, mi_1009, mi_1010, mi_1011, mi_1012, \
                         mi_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = -7.0 * ki_589[k]
                   + f_0 * mi_1009[k];

        t_786[k] = -7.0 * ki_590[k]
                   + f_0 * mi_1010[k];

        t_787[k] = -7.0 * ki_591[k]
                   + f_0 * mi_1011[k];

        t_788[k] = -7.0 * ki_592[k]
                   + f_0 * mi_1012[k];

        t_789[k] = -7.0 * ki_593[k]
                   + f_0 * mi_1013[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, t_794, ki_594, ki_595, ki_596, ki_597, \
                         ki_598, mi_1014, mi_1015, mi_1016, mi_1017, \
                         mi_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = -7.0 * ki_594[k]
                   + f_0 * mi_1014[k];

        t_791[k] = -7.0 * ki_595[k]
                   + f_0 * mi_1015[k];

        t_792[k] = -7.0 * ki_596[k]
                   + f_0 * mi_1016[k];

        t_793[k] = -7.0 * ki_597[k]
                   + f_0 * mi_1017[k];

        t_794[k] = -7.0 * ki_598[k]
                   + f_0 * mi_1018[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, ki_599, ki_600, ki_601, ki_602, \
                         ki_603, mi_1019, mi_1020, mi_1021, mi_1022, \
                         mi_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = -7.0 * ki_599[k]
                   + f_0 * mi_1019[k];

        t_796[k] = -7.0 * ki_600[k]
                   + f_0 * mi_1020[k];

        t_797[k] = -7.0 * ki_601[k]
                   + f_0 * mi_1021[k];

        t_798[k] = -7.0 * ki_602[k]
                   + f_0 * mi_1022[k];

        t_799[k] = -7.0 * ki_603[k]
                   + f_0 * mi_1023[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, ki_604, ki_605, ki_606, ki_607, \
                         ki_608, mi_1024, mi_1025, mi_1026, mi_1027, \
                         mi_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = -7.0 * ki_604[k]
                   + f_0 * mi_1024[k];

        t_801[k] = -7.0 * ki_605[k]
                   + f_0 * mi_1025[k];

        t_802[k] = -7.0 * ki_606[k]
                   + f_0 * mi_1026[k];

        t_803[k] = -7.0 * ki_607[k]
                   + f_0 * mi_1027[k];

        t_804[k] = -7.0 * ki_608[k]
                   + f_0 * mi_1028[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, ki_609, ki_610, ki_611, ki_612, \
                         ki_613, mi_1029, mi_1030, mi_1031, mi_1032, \
                         mi_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = -7.0 * ki_609[k]
                   + f_0 * mi_1029[k];

        t_806[k] = -7.0 * ki_610[k]
                   + f_0 * mi_1030[k];

        t_807[k] = -7.0 * ki_611[k]
                   + f_0 * mi_1031[k];

        t_808[k] = -7.0 * ki_612[k]
                   + f_0 * mi_1032[k];

        t_809[k] = -7.0 * ki_613[k]
                   + f_0 * mi_1033[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, ki_614, ki_615, ki_616, ki_617, \
                         ki_618, mi_1034, mi_1035, mi_1036, mi_1037, \
                         mi_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = -7.0 * ki_614[k]
                   + f_0 * mi_1034[k];

        t_811[k] = -7.0 * ki_615[k]
                   + f_0 * mi_1035[k];

        t_812[k] = -6.0 * ki_616[k]
                   + f_0 * mi_1036[k];

        t_813[k] = -6.0 * ki_617[k]
                   + f_0 * mi_1037[k];

        t_814[k] = -6.0 * ki_618[k]
                   + f_0 * mi_1038[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_1_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);
    auto *t_948 = buffer.data(target + 948);
    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);
    auto *t_955 = buffer.data(target + 955);
    auto *t_956 = buffer.data(target + 956);
    auto *t_957 = buffer.data(target + 957);
    auto *t_958 = buffer.data(target + 958);
    auto *t_959 = buffer.data(target + 959);
    auto *t_960 = buffer.data(target + 960);
    auto *t_961 = buffer.data(target + 961);
    auto *t_962 = buffer.data(target + 962);
    auto *t_963 = buffer.data(target + 963);
    auto *t_964 = buffer.data(target + 964);

    const auto *ki_619 = buffer.data(ki + 619);
    const auto *ki_620 = buffer.data(ki + 620);
    const auto *ki_621 = buffer.data(ki + 621);
    const auto *ki_622 = buffer.data(ki + 622);
    const auto *ki_623 = buffer.data(ki + 623);
    const auto *ki_624 = buffer.data(ki + 624);
    const auto *ki_625 = buffer.data(ki + 625);
    const auto *ki_626 = buffer.data(ki + 626);
    const auto *ki_627 = buffer.data(ki + 627);
    const auto *ki_628 = buffer.data(ki + 628);
    const auto *ki_629 = buffer.data(ki + 629);
    const auto *ki_630 = buffer.data(ki + 630);
    const auto *ki_631 = buffer.data(ki + 631);
    const auto *ki_632 = buffer.data(ki + 632);
    const auto *ki_633 = buffer.data(ki + 633);
    const auto *ki_634 = buffer.data(ki + 634);
    const auto *ki_635 = buffer.data(ki + 635);
    const auto *ki_636 = buffer.data(ki + 636);
    const auto *ki_637 = buffer.data(ki + 637);
    const auto *ki_638 = buffer.data(ki + 638);
    const auto *ki_639 = buffer.data(ki + 639);
    const auto *ki_640 = buffer.data(ki + 640);
    const auto *ki_641 = buffer.data(ki + 641);
    const auto *ki_642 = buffer.data(ki + 642);
    const auto *ki_643 = buffer.data(ki + 643);
    const auto *ki_644 = buffer.data(ki + 644);
    const auto *ki_645 = buffer.data(ki + 645);
    const auto *ki_646 = buffer.data(ki + 646);
    const auto *ki_647 = buffer.data(ki + 647);
    const auto *ki_648 = buffer.data(ki + 648);
    const auto *ki_649 = buffer.data(ki + 649);
    const auto *ki_650 = buffer.data(ki + 650);
    const auto *ki_651 = buffer.data(ki + 651);
    const auto *ki_652 = buffer.data(ki + 652);
    const auto *ki_653 = buffer.data(ki + 653);
    const auto *ki_654 = buffer.data(ki + 654);
    const auto *ki_655 = buffer.data(ki + 655);
    const auto *ki_656 = buffer.data(ki + 656);
    const auto *ki_657 = buffer.data(ki + 657);
    const auto *ki_658 = buffer.data(ki + 658);
    const auto *ki_659 = buffer.data(ki + 659);
    const auto *ki_660 = buffer.data(ki + 660);
    const auto *ki_661 = buffer.data(ki + 661);
    const auto *ki_662 = buffer.data(ki + 662);
    const auto *ki_663 = buffer.data(ki + 663);
    const auto *ki_664 = buffer.data(ki + 664);
    const auto *ki_665 = buffer.data(ki + 665);
    const auto *ki_666 = buffer.data(ki + 666);
    const auto *ki_667 = buffer.data(ki + 667);
    const auto *ki_668 = buffer.data(ki + 668);
    const auto *ki_669 = buffer.data(ki + 669);
    const auto *ki_670 = buffer.data(ki + 670);
    const auto *ki_671 = buffer.data(ki + 671);
    const auto *ki_672 = buffer.data(ki + 672);
    const auto *ki_673 = buffer.data(ki + 673);
    const auto *ki_674 = buffer.data(ki + 674);
    const auto *ki_675 = buffer.data(ki + 675);
    const auto *ki_676 = buffer.data(ki + 676);
    const auto *ki_677 = buffer.data(ki + 677);
    const auto *ki_678 = buffer.data(ki + 678);
    const auto *ki_679 = buffer.data(ki + 679);
    const auto *ki_680 = buffer.data(ki + 680);
    const auto *ki_681 = buffer.data(ki + 681);
    const auto *ki_682 = buffer.data(ki + 682);
    const auto *ki_683 = buffer.data(ki + 683);
    const auto *ki_684 = buffer.data(ki + 684);
    const auto *ki_685 = buffer.data(ki + 685);
    const auto *ki_686 = buffer.data(ki + 686);
    const auto *ki_687 = buffer.data(ki + 687);
    const auto *ki_688 = buffer.data(ki + 688);
    const auto *ki_689 = buffer.data(ki + 689);
    const auto *ki_690 = buffer.data(ki + 690);
    const auto *ki_691 = buffer.data(ki + 691);
    const auto *ki_692 = buffer.data(ki + 692);
    const auto *ki_693 = buffer.data(ki + 693);
    const auto *ki_694 = buffer.data(ki + 694);
    const auto *ki_695 = buffer.data(ki + 695);
    const auto *ki_696 = buffer.data(ki + 696);
    const auto *ki_697 = buffer.data(ki + 697);
    const auto *ki_698 = buffer.data(ki + 698);
    const auto *ki_699 = buffer.data(ki + 699);
    const auto *ki_700 = buffer.data(ki + 700);
    const auto *ki_701 = buffer.data(ki + 701);
    const auto *ki_702 = buffer.data(ki + 702);
    const auto *ki_703 = buffer.data(ki + 703);
    const auto *ki_704 = buffer.data(ki + 704);
    const auto *ki_705 = buffer.data(ki + 705);
    const auto *ki_706 = buffer.data(ki + 706);
    const auto *ki_707 = buffer.data(ki + 707);
    const auto *ki_708 = buffer.data(ki + 708);
    const auto *ki_709 = buffer.data(ki + 709);
    const auto *ki_710 = buffer.data(ki + 710);
    const auto *ki_711 = buffer.data(ki + 711);
    const auto *ki_712 = buffer.data(ki + 712);
    const auto *ki_713 = buffer.data(ki + 713);
    const auto *ki_714 = buffer.data(ki + 714);
    const auto *ki_715 = buffer.data(ki + 715);
    const auto *ki_716 = buffer.data(ki + 716);
    const auto *ki_717 = buffer.data(ki + 717);
    const auto *ki_718 = buffer.data(ki + 718);
    const auto *ki_719 = buffer.data(ki + 719);
    const auto *ki_720 = buffer.data(ki + 720);
    const auto *ki_721 = buffer.data(ki + 721);
    const auto *ki_722 = buffer.data(ki + 722);
    const auto *ki_723 = buffer.data(ki + 723);
    const auto *ki_724 = buffer.data(ki + 724);
    const auto *ki_725 = buffer.data(ki + 725);
    const auto *ki_726 = buffer.data(ki + 726);
    const auto *ki_727 = buffer.data(ki + 727);
    const auto *ki_728 = buffer.data(ki + 728);
    const auto *ki_729 = buffer.data(ki + 729);
    const auto *ki_730 = buffer.data(ki + 730);
    const auto *ki_731 = buffer.data(ki + 731);
    const auto *ki_732 = buffer.data(ki + 732);
    const auto *ki_733 = buffer.data(ki + 733);
    const auto *ki_734 = buffer.data(ki + 734);
    const auto *ki_735 = buffer.data(ki + 735);
    const auto *ki_736 = buffer.data(ki + 736);
    const auto *ki_737 = buffer.data(ki + 737);
    const auto *ki_738 = buffer.data(ki + 738);
    const auto *ki_739 = buffer.data(ki + 739);
    const auto *ki_740 = buffer.data(ki + 740);
    const auto *ki_741 = buffer.data(ki + 741);
    const auto *ki_742 = buffer.data(ki + 742);
    const auto *ki_743 = buffer.data(ki + 743);
    const auto *ki_744 = buffer.data(ki + 744);
    const auto *ki_745 = buffer.data(ki + 745);
    const auto *ki_746 = buffer.data(ki + 746);
    const auto *ki_747 = buffer.data(ki + 747);
    const auto *ki_748 = buffer.data(ki + 748);
    const auto *ki_749 = buffer.data(ki + 749);
    const auto *ki_750 = buffer.data(ki + 750);
    const auto *ki_751 = buffer.data(ki + 751);
    const auto *ki_752 = buffer.data(ki + 752);
    const auto *ki_753 = buffer.data(ki + 753);
    const auto *ki_754 = buffer.data(ki + 754);
    const auto *ki_755 = buffer.data(ki + 755);
    const auto *ki_756 = buffer.data(ki + 756);
    const auto *ki_757 = buffer.data(ki + 757);
    const auto *ki_758 = buffer.data(ki + 758);
    const auto *ki_759 = buffer.data(ki + 759);
    const auto *ki_760 = buffer.data(ki + 760);
    const auto *ki_761 = buffer.data(ki + 761);
    const auto *ki_762 = buffer.data(ki + 762);
    const auto *ki_763 = buffer.data(ki + 763);
    const auto *ki_764 = buffer.data(ki + 764);
    const auto *ki_765 = buffer.data(ki + 765);
    const auto *ki_766 = buffer.data(ki + 766);
    const auto *ki_767 = buffer.data(ki + 767);
    const auto *ki_768 = buffer.data(ki + 768);

    const auto *mi_1039 = buffer.data(mi + 1039);
    const auto *mi_1040 = buffer.data(mi + 1040);
    const auto *mi_1041 = buffer.data(mi + 1041);
    const auto *mi_1042 = buffer.data(mi + 1042);
    const auto *mi_1043 = buffer.data(mi + 1043);
    const auto *mi_1044 = buffer.data(mi + 1044);
    const auto *mi_1045 = buffer.data(mi + 1045);
    const auto *mi_1046 = buffer.data(mi + 1046);
    const auto *mi_1047 = buffer.data(mi + 1047);
    const auto *mi_1048 = buffer.data(mi + 1048);
    const auto *mi_1049 = buffer.data(mi + 1049);
    const auto *mi_1050 = buffer.data(mi + 1050);
    const auto *mi_1051 = buffer.data(mi + 1051);
    const auto *mi_1052 = buffer.data(mi + 1052);
    const auto *mi_1053 = buffer.data(mi + 1053);
    const auto *mi_1054 = buffer.data(mi + 1054);
    const auto *mi_1055 = buffer.data(mi + 1055);
    const auto *mi_1056 = buffer.data(mi + 1056);
    const auto *mi_1057 = buffer.data(mi + 1057);
    const auto *mi_1058 = buffer.data(mi + 1058);
    const auto *mi_1059 = buffer.data(mi + 1059);
    const auto *mi_1060 = buffer.data(mi + 1060);
    const auto *mi_1061 = buffer.data(mi + 1061);
    const auto *mi_1062 = buffer.data(mi + 1062);
    const auto *mi_1063 = buffer.data(mi + 1063);
    const auto *mi_1064 = buffer.data(mi + 1064);
    const auto *mi_1065 = buffer.data(mi + 1065);
    const auto *mi_1066 = buffer.data(mi + 1066);
    const auto *mi_1067 = buffer.data(mi + 1067);
    const auto *mi_1068 = buffer.data(mi + 1068);
    const auto *mi_1069 = buffer.data(mi + 1069);
    const auto *mi_1070 = buffer.data(mi + 1070);
    const auto *mi_1071 = buffer.data(mi + 1071);
    const auto *mi_1072 = buffer.data(mi + 1072);
    const auto *mi_1073 = buffer.data(mi + 1073);
    const auto *mi_1074 = buffer.data(mi + 1074);
    const auto *mi_1075 = buffer.data(mi + 1075);
    const auto *mi_1076 = buffer.data(mi + 1076);
    const auto *mi_1077 = buffer.data(mi + 1077);
    const auto *mi_1078 = buffer.data(mi + 1078);
    const auto *mi_1079 = buffer.data(mi + 1079);
    const auto *mi_1080 = buffer.data(mi + 1080);
    const auto *mi_1081 = buffer.data(mi + 1081);
    const auto *mi_1082 = buffer.data(mi + 1082);
    const auto *mi_1083 = buffer.data(mi + 1083);
    const auto *mi_1084 = buffer.data(mi + 1084);
    const auto *mi_1085 = buffer.data(mi + 1085);
    const auto *mi_1086 = buffer.data(mi + 1086);
    const auto *mi_1087 = buffer.data(mi + 1087);
    const auto *mi_1088 = buffer.data(mi + 1088);
    const auto *mi_1089 = buffer.data(mi + 1089);
    const auto *mi_1090 = buffer.data(mi + 1090);
    const auto *mi_1091 = buffer.data(mi + 1091);
    const auto *mi_1092 = buffer.data(mi + 1092);
    const auto *mi_1093 = buffer.data(mi + 1093);
    const auto *mi_1094 = buffer.data(mi + 1094);
    const auto *mi_1095 = buffer.data(mi + 1095);
    const auto *mi_1096 = buffer.data(mi + 1096);
    const auto *mi_1097 = buffer.data(mi + 1097);
    const auto *mi_1098 = buffer.data(mi + 1098);
    const auto *mi_1099 = buffer.data(mi + 1099);
    const auto *mi_1100 = buffer.data(mi + 1100);
    const auto *mi_1101 = buffer.data(mi + 1101);
    const auto *mi_1102 = buffer.data(mi + 1102);
    const auto *mi_1103 = buffer.data(mi + 1103);
    const auto *mi_1104 = buffer.data(mi + 1104);
    const auto *mi_1105 = buffer.data(mi + 1105);
    const auto *mi_1106 = buffer.data(mi + 1106);
    const auto *mi_1107 = buffer.data(mi + 1107);
    const auto *mi_1108 = buffer.data(mi + 1108);
    const auto *mi_1109 = buffer.data(mi + 1109);
    const auto *mi_1110 = buffer.data(mi + 1110);
    const auto *mi_1111 = buffer.data(mi + 1111);
    const auto *mi_1112 = buffer.data(mi + 1112);
    const auto *mi_1113 = buffer.data(mi + 1113);
    const auto *mi_1114 = buffer.data(mi + 1114);
    const auto *mi_1115 = buffer.data(mi + 1115);
    const auto *mi_1116 = buffer.data(mi + 1116);
    const auto *mi_1117 = buffer.data(mi + 1117);
    const auto *mi_1118 = buffer.data(mi + 1118);
    const auto *mi_1119 = buffer.data(mi + 1119);
    const auto *mi_1120 = buffer.data(mi + 1120);
    const auto *mi_1121 = buffer.data(mi + 1121);
    const auto *mi_1122 = buffer.data(mi + 1122);
    const auto *mi_1123 = buffer.data(mi + 1123);
    const auto *mi_1124 = buffer.data(mi + 1124);
    const auto *mi_1125 = buffer.data(mi + 1125);
    const auto *mi_1126 = buffer.data(mi + 1126);
    const auto *mi_1127 = buffer.data(mi + 1127);
    const auto *mi_1128 = buffer.data(mi + 1128);
    const auto *mi_1129 = buffer.data(mi + 1129);
    const auto *mi_1130 = buffer.data(mi + 1130);
    const auto *mi_1131 = buffer.data(mi + 1131);
    const auto *mi_1132 = buffer.data(mi + 1132);
    const auto *mi_1133 = buffer.data(mi + 1133);
    const auto *mi_1134 = buffer.data(mi + 1134);
    const auto *mi_1135 = buffer.data(mi + 1135);
    const auto *mi_1136 = buffer.data(mi + 1136);
    const auto *mi_1137 = buffer.data(mi + 1137);
    const auto *mi_1138 = buffer.data(mi + 1138);
    const auto *mi_1139 = buffer.data(mi + 1139);
    const auto *mi_1140 = buffer.data(mi + 1140);
    const auto *mi_1141 = buffer.data(mi + 1141);
    const auto *mi_1142 = buffer.data(mi + 1142);
    const auto *mi_1143 = buffer.data(mi + 1143);
    const auto *mi_1144 = buffer.data(mi + 1144);
    const auto *mi_1145 = buffer.data(mi + 1145);
    const auto *mi_1146 = buffer.data(mi + 1146);
    const auto *mi_1147 = buffer.data(mi + 1147);
    const auto *mi_1148 = buffer.data(mi + 1148);
    const auto *mi_1149 = buffer.data(mi + 1149);
    const auto *mi_1150 = buffer.data(mi + 1150);
    const auto *mi_1151 = buffer.data(mi + 1151);
    const auto *mi_1152 = buffer.data(mi + 1152);
    const auto *mi_1153 = buffer.data(mi + 1153);
    const auto *mi_1154 = buffer.data(mi + 1154);
    const auto *mi_1155 = buffer.data(mi + 1155);
    const auto *mi_1156 = buffer.data(mi + 1156);
    const auto *mi_1157 = buffer.data(mi + 1157);
    const auto *mi_1158 = buffer.data(mi + 1158);
    const auto *mi_1159 = buffer.data(mi + 1159);
    const auto *mi_1160 = buffer.data(mi + 1160);
    const auto *mi_1161 = buffer.data(mi + 1161);
    const auto *mi_1162 = buffer.data(mi + 1162);
    const auto *mi_1163 = buffer.data(mi + 1163);
    const auto *mi_1164 = buffer.data(mi + 1164);
    const auto *mi_1165 = buffer.data(mi + 1165);
    const auto *mi_1166 = buffer.data(mi + 1166);
    const auto *mi_1167 = buffer.data(mi + 1167);
    const auto *mi_1168 = buffer.data(mi + 1168);
    const auto *mi_1169 = buffer.data(mi + 1169);
    const auto *mi_1170 = buffer.data(mi + 1170);
    const auto *mi_1171 = buffer.data(mi + 1171);
    const auto *mi_1172 = buffer.data(mi + 1172);
    const auto *mi_1173 = buffer.data(mi + 1173);
    const auto *mi_1174 = buffer.data(mi + 1174);
    const auto *mi_1175 = buffer.data(mi + 1175);
    const auto *mi_1176 = buffer.data(mi + 1176);
    const auto *mi_1177 = buffer.data(mi + 1177);
    const auto *mi_1178 = buffer.data(mi + 1178);
    const auto *mi_1179 = buffer.data(mi + 1179);
    const auto *mi_1180 = buffer.data(mi + 1180);
    const auto *mi_1181 = buffer.data(mi + 1181);
    const auto *mi_1182 = buffer.data(mi + 1182);
    const auto *mi_1183 = buffer.data(mi + 1183);
    const auto *mi_1184 = buffer.data(mi + 1184);
    const auto *mi_1185 = buffer.data(mi + 1185);
    const auto *mi_1186 = buffer.data(mi + 1186);
    const auto *mi_1187 = buffer.data(mi + 1187);
    const auto *mi_1188 = buffer.data(mi + 1188);

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, ki_619, ki_620, ki_621, ki_622, \
                         ki_623, mi_1039, mi_1040, mi_1041, mi_1042, \
                         mi_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = -6.0 * ki_619[k]
                   + f_0 * mi_1039[k];

        t_816[k] = -6.0 * ki_620[k]
                   + f_0 * mi_1040[k];

        t_817[k] = -6.0 * ki_621[k]
                   + f_0 * mi_1041[k];

        t_818[k] = -6.0 * ki_622[k]
                   + f_0 * mi_1042[k];

        t_819[k] = -6.0 * ki_623[k]
                   + f_0 * mi_1043[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, t_824, ki_624, ki_625, ki_626, ki_627, \
                         ki_628, mi_1044, mi_1045, mi_1046, mi_1047, \
                         mi_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = -6.0 * ki_624[k]
                   + f_0 * mi_1044[k];

        t_821[k] = -6.0 * ki_625[k]
                   + f_0 * mi_1045[k];

        t_822[k] = -6.0 * ki_626[k]
                   + f_0 * mi_1046[k];

        t_823[k] = -6.0 * ki_627[k]
                   + f_0 * mi_1047[k];

        t_824[k] = -6.0 * ki_628[k]
                   + f_0 * mi_1048[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, ki_629, ki_630, ki_631, ki_632, \
                         ki_633, mi_1049, mi_1050, mi_1051, mi_1052, \
                         mi_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = -6.0 * ki_629[k]
                   + f_0 * mi_1049[k];

        t_826[k] = -6.0 * ki_630[k]
                   + f_0 * mi_1050[k];

        t_827[k] = -6.0 * ki_631[k]
                   + f_0 * mi_1051[k];

        t_828[k] = -6.0 * ki_632[k]
                   + f_0 * mi_1052[k];

        t_829[k] = -6.0 * ki_633[k]
                   + f_0 * mi_1053[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, ki_634, ki_635, ki_636, ki_637, \
                         ki_638, mi_1054, mi_1055, mi_1056, mi_1057, \
                         mi_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = -6.0 * ki_634[k]
                   + f_0 * mi_1054[k];

        t_831[k] = -6.0 * ki_635[k]
                   + f_0 * mi_1055[k];

        t_832[k] = -6.0 * ki_636[k]
                   + f_0 * mi_1056[k];

        t_833[k] = -6.0 * ki_637[k]
                   + f_0 * mi_1057[k];

        t_834[k] = -6.0 * ki_638[k]
                   + f_0 * mi_1058[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, t_838, t_839, ki_639, ki_640, ki_641, ki_642, \
                         ki_643, mi_1059, mi_1060, mi_1061, mi_1062, \
                         mi_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = -6.0 * ki_639[k]
                   + f_0 * mi_1059[k];

        t_836[k] = -6.0 * ki_640[k]
                   + f_0 * mi_1060[k];

        t_837[k] = -6.0 * ki_641[k]
                   + f_0 * mi_1061[k];

        t_838[k] = -6.0 * ki_642[k]
                   + f_0 * mi_1062[k];

        t_839[k] = -6.0 * ki_643[k]
                   + f_0 * mi_1063[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, ki_644, ki_645, ki_646, ki_647, \
                         ki_648, mi_1064, mi_1065, mi_1066, mi_1067, \
                         mi_1068 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = -5.0 * ki_644[k]
                   + f_0 * mi_1064[k];

        t_841[k] = -5.0 * ki_645[k]
                   + f_0 * mi_1065[k];

        t_842[k] = -5.0 * ki_646[k]
                   + f_0 * mi_1066[k];

        t_843[k] = -5.0 * ki_647[k]
                   + f_0 * mi_1067[k];

        t_844[k] = -5.0 * ki_648[k]
                   + f_0 * mi_1068[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, t_849, ki_649, ki_650, ki_651, ki_652, \
                         ki_653, mi_1069, mi_1070, mi_1071, mi_1072, \
                         mi_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = -5.0 * ki_649[k]
                   + f_0 * mi_1069[k];

        t_846[k] = -5.0 * ki_650[k]
                   + f_0 * mi_1070[k];

        t_847[k] = -5.0 * ki_651[k]
                   + f_0 * mi_1071[k];

        t_848[k] = -5.0 * ki_652[k]
                   + f_0 * mi_1072[k];

        t_849[k] = -5.0 * ki_653[k]
                   + f_0 * mi_1073[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, ki_654, ki_655, ki_656, ki_657, \
                         ki_658, mi_1074, mi_1075, mi_1076, mi_1077, \
                         mi_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = -5.0 * ki_654[k]
                   + f_0 * mi_1074[k];

        t_851[k] = -5.0 * ki_655[k]
                   + f_0 * mi_1075[k];

        t_852[k] = -5.0 * ki_656[k]
                   + f_0 * mi_1076[k];

        t_853[k] = -5.0 * ki_657[k]
                   + f_0 * mi_1077[k];

        t_854[k] = -5.0 * ki_658[k]
                   + f_0 * mi_1078[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, ki_659, ki_660, ki_661, ki_662, \
                         ki_663, mi_1079, mi_1080, mi_1081, mi_1082, \
                         mi_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = -5.0 * ki_659[k]
                   + f_0 * mi_1079[k];

        t_856[k] = -5.0 * ki_660[k]
                   + f_0 * mi_1080[k];

        t_857[k] = -5.0 * ki_661[k]
                   + f_0 * mi_1081[k];

        t_858[k] = -5.0 * ki_662[k]
                   + f_0 * mi_1082[k];

        t_859[k] = -5.0 * ki_663[k]
                   + f_0 * mi_1083[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, t_864, ki_664, ki_665, ki_666, ki_667, \
                         ki_668, mi_1084, mi_1085, mi_1086, mi_1087, \
                         mi_1088 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = -5.0 * ki_664[k]
                   + f_0 * mi_1084[k];

        t_861[k] = -5.0 * ki_665[k]
                   + f_0 * mi_1085[k];

        t_862[k] = -5.0 * ki_666[k]
                   + f_0 * mi_1086[k];

        t_863[k] = -5.0 * ki_667[k]
                   + f_0 * mi_1087[k];

        t_864[k] = -5.0 * ki_668[k]
                   + f_0 * mi_1088[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, t_869, ki_669, ki_670, ki_671, ki_672, \
                         ki_673, mi_1089, mi_1090, mi_1091, mi_1092, \
                         mi_1093 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = -5.0 * ki_669[k]
                   + f_0 * mi_1089[k];

        t_866[k] = -5.0 * ki_670[k]
                   + f_0 * mi_1090[k];

        t_867[k] = -5.0 * ki_671[k]
                   + f_0 * mi_1091[k];

        t_868[k] = -4.0 * ki_672[k]
                   + f_0 * mi_1092[k];

        t_869[k] = -4.0 * ki_673[k]
                   + f_0 * mi_1093[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, t_874, ki_674, ki_675, ki_676, ki_677, \
                         ki_678, mi_1094, mi_1095, mi_1096, mi_1097, \
                         mi_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = -4.0 * ki_674[k]
                   + f_0 * mi_1094[k];

        t_871[k] = -4.0 * ki_675[k]
                   + f_0 * mi_1095[k];

        t_872[k] = -4.0 * ki_676[k]
                   + f_0 * mi_1096[k];

        t_873[k] = -4.0 * ki_677[k]
                   + f_0 * mi_1097[k];

        t_874[k] = -4.0 * ki_678[k]
                   + f_0 * mi_1098[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, t_879, ki_679, ki_680, ki_681, ki_682, \
                         ki_683, mi_1099, mi_1100, mi_1101, mi_1102, \
                         mi_1103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = -4.0 * ki_679[k]
                   + f_0 * mi_1099[k];

        t_876[k] = -4.0 * ki_680[k]
                   + f_0 * mi_1100[k];

        t_877[k] = -4.0 * ki_681[k]
                   + f_0 * mi_1101[k];

        t_878[k] = -4.0 * ki_682[k]
                   + f_0 * mi_1102[k];

        t_879[k] = -4.0 * ki_683[k]
                   + f_0 * mi_1103[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, t_884, ki_684, ki_685, ki_686, ki_687, \
                         ki_688, mi_1104, mi_1105, mi_1106, mi_1107, \
                         mi_1108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = -4.0 * ki_684[k]
                   + f_0 * mi_1104[k];

        t_881[k] = -4.0 * ki_685[k]
                   + f_0 * mi_1105[k];

        t_882[k] = -4.0 * ki_686[k]
                   + f_0 * mi_1106[k];

        t_883[k] = -4.0 * ki_687[k]
                   + f_0 * mi_1107[k];

        t_884[k] = -4.0 * ki_688[k]
                   + f_0 * mi_1108[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, ki_689, ki_690, ki_691, ki_692, \
                         ki_693, mi_1109, mi_1110, mi_1111, mi_1112, \
                         mi_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = -4.0 * ki_689[k]
                   + f_0 * mi_1109[k];

        t_886[k] = -4.0 * ki_690[k]
                   + f_0 * mi_1110[k];

        t_887[k] = -4.0 * ki_691[k]
                   + f_0 * mi_1111[k];

        t_888[k] = -4.0 * ki_692[k]
                   + f_0 * mi_1112[k];

        t_889[k] = -4.0 * ki_693[k]
                   + f_0 * mi_1113[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, ki_694, ki_695, ki_696, ki_697, \
                         ki_698, mi_1114, mi_1115, mi_1116, mi_1117, \
                         mi_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = -4.0 * ki_694[k]
                   + f_0 * mi_1114[k];

        t_891[k] = -4.0 * ki_695[k]
                   + f_0 * mi_1115[k];

        t_892[k] = -4.0 * ki_696[k]
                   + f_0 * mi_1116[k];

        t_893[k] = -4.0 * ki_697[k]
                   + f_0 * mi_1117[k];

        t_894[k] = -4.0 * ki_698[k]
                   + f_0 * mi_1118[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, t_899, ki_699, ki_700, ki_701, ki_702, \
                         ki_703, mi_1119, mi_1120, mi_1121, mi_1122, \
                         mi_1123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = -4.0 * ki_699[k]
                   + f_0 * mi_1119[k];

        t_896[k] = -3.0 * ki_700[k]
                   + f_0 * mi_1120[k];

        t_897[k] = -3.0 * ki_701[k]
                   + f_0 * mi_1121[k];

        t_898[k] = -3.0 * ki_702[k]
                   + f_0 * mi_1122[k];

        t_899[k] = -3.0 * ki_703[k]
                   + f_0 * mi_1123[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, t_904, ki_704, ki_705, ki_706, ki_707, \
                         ki_708, mi_1124, mi_1125, mi_1126, mi_1127, \
                         mi_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = -3.0 * ki_704[k]
                   + f_0 * mi_1124[k];

        t_901[k] = -3.0 * ki_705[k]
                   + f_0 * mi_1125[k];

        t_902[k] = -3.0 * ki_706[k]
                   + f_0 * mi_1126[k];

        t_903[k] = -3.0 * ki_707[k]
                   + f_0 * mi_1127[k];

        t_904[k] = -3.0 * ki_708[k]
                   + f_0 * mi_1128[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, t_909, ki_709, ki_710, ki_711, ki_712, \
                         ki_713, mi_1129, mi_1130, mi_1131, mi_1132, \
                         mi_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = -3.0 * ki_709[k]
                   + f_0 * mi_1129[k];

        t_906[k] = -3.0 * ki_710[k]
                   + f_0 * mi_1130[k];

        t_907[k] = -3.0 * ki_711[k]
                   + f_0 * mi_1131[k];

        t_908[k] = -3.0 * ki_712[k]
                   + f_0 * mi_1132[k];

        t_909[k] = -3.0 * ki_713[k]
                   + f_0 * mi_1133[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, t_914, ki_714, ki_715, ki_716, ki_717, \
                         ki_718, mi_1134, mi_1135, mi_1136, mi_1137, \
                         mi_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = -3.0 * ki_714[k]
                   + f_0 * mi_1134[k];

        t_911[k] = -3.0 * ki_715[k]
                   + f_0 * mi_1135[k];

        t_912[k] = -3.0 * ki_716[k]
                   + f_0 * mi_1136[k];

        t_913[k] = -3.0 * ki_717[k]
                   + f_0 * mi_1137[k];

        t_914[k] = -3.0 * ki_718[k]
                   + f_0 * mi_1138[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, t_919, ki_719, ki_720, ki_721, ki_722, \
                         ki_723, mi_1139, mi_1140, mi_1141, mi_1142, \
                         mi_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = -3.0 * ki_719[k]
                   + f_0 * mi_1139[k];

        t_916[k] = -3.0 * ki_720[k]
                   + f_0 * mi_1140[k];

        t_917[k] = -3.0 * ki_721[k]
                   + f_0 * mi_1141[k];

        t_918[k] = -3.0 * ki_722[k]
                   + f_0 * mi_1142[k];

        t_919[k] = -3.0 * ki_723[k]
                   + f_0 * mi_1143[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, t_924, ki_724, ki_725, ki_726, ki_727, \
                         ki_728, mi_1144, mi_1145, mi_1146, mi_1147, \
                         mi_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = -3.0 * ki_724[k]
                   + f_0 * mi_1144[k];

        t_921[k] = -3.0 * ki_725[k]
                   + f_0 * mi_1145[k];

        t_922[k] = -3.0 * ki_726[k]
                   + f_0 * mi_1146[k];

        t_923[k] = -3.0 * ki_727[k]
                   + f_0 * mi_1147[k];

        t_924[k] = -2.0 * ki_728[k]
                   + f_0 * mi_1148[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, t_929, ki_729, ki_730, ki_731, ki_732, \
                         ki_733, mi_1149, mi_1150, mi_1151, mi_1152, \
                         mi_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = -2.0 * ki_729[k]
                   + f_0 * mi_1149[k];

        t_926[k] = -2.0 * ki_730[k]
                   + f_0 * mi_1150[k];

        t_927[k] = -2.0 * ki_731[k]
                   + f_0 * mi_1151[k];

        t_928[k] = -2.0 * ki_732[k]
                   + f_0 * mi_1152[k];

        t_929[k] = -2.0 * ki_733[k]
                   + f_0 * mi_1153[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, ki_734, ki_735, ki_736, ki_737, \
                         ki_738, mi_1154, mi_1155, mi_1156, mi_1157, \
                         mi_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = -2.0 * ki_734[k]
                   + f_0 * mi_1154[k];

        t_931[k] = -2.0 * ki_735[k]
                   + f_0 * mi_1155[k];

        t_932[k] = -2.0 * ki_736[k]
                   + f_0 * mi_1156[k];

        t_933[k] = -2.0 * ki_737[k]
                   + f_0 * mi_1157[k];

        t_934[k] = -2.0 * ki_738[k]
                   + f_0 * mi_1158[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, t_939, ki_739, ki_740, ki_741, ki_742, \
                         ki_743, mi_1159, mi_1160, mi_1161, mi_1162, \
                         mi_1163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = -2.0 * ki_739[k]
                   + f_0 * mi_1159[k];

        t_936[k] = -2.0 * ki_740[k]
                   + f_0 * mi_1160[k];

        t_937[k] = -2.0 * ki_741[k]
                   + f_0 * mi_1161[k];

        t_938[k] = -2.0 * ki_742[k]
                   + f_0 * mi_1162[k];

        t_939[k] = -2.0 * ki_743[k]
                   + f_0 * mi_1163[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, t_944, ki_744, ki_745, ki_746, ki_747, \
                         ki_748, mi_1164, mi_1165, mi_1166, mi_1167, \
                         mi_1168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = -2.0 * ki_744[k]
                   + f_0 * mi_1164[k];

        t_941[k] = -2.0 * ki_745[k]
                   + f_0 * mi_1165[k];

        t_942[k] = -2.0 * ki_746[k]
                   + f_0 * mi_1166[k];

        t_943[k] = -2.0 * ki_747[k]
                   + f_0 * mi_1167[k];

        t_944[k] = -2.0 * ki_748[k]
                   + f_0 * mi_1168[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, t_949, ki_749, ki_750, ki_751, ki_752, \
                         ki_753, mi_1169, mi_1170, mi_1171, mi_1172, \
                         mi_1173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = -2.0 * ki_749[k]
                   + f_0 * mi_1169[k];

        t_946[k] = -2.0 * ki_750[k]
                   + f_0 * mi_1170[k];

        t_947[k] = -2.0 * ki_751[k]
                   + f_0 * mi_1171[k];

        t_948[k] = -2.0 * ki_752[k]
                   + f_0 * mi_1172[k];

        t_949[k] = -2.0 * ki_753[k]
                   + f_0 * mi_1173[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, t_954, ki_754, ki_755, ki_756, ki_757, \
                         ki_758, mi_1174, mi_1175, mi_1176, mi_1177, \
                         mi_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = -2.0 * ki_754[k]
                   + f_0 * mi_1174[k];

        t_951[k] = -2.0 * ki_755[k]
                   + f_0 * mi_1175[k];

        t_952[k] = -ki_756[k]
                   + f_0 * mi_1176[k];

        t_953[k] = -ki_757[k]
                   + f_0 * mi_1177[k];

        t_954[k] = -ki_758[k]
                   + f_0 * mi_1178[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, t_959, ki_759, ki_760, ki_761, ki_762, \
                         ki_763, mi_1179, mi_1180, mi_1181, mi_1182, \
                         mi_1183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = -ki_759[k]
                   + f_0 * mi_1179[k];

        t_956[k] = -ki_760[k]
                   + f_0 * mi_1180[k];

        t_957[k] = -ki_761[k]
                   + f_0 * mi_1181[k];

        t_958[k] = -ki_762[k]
                   + f_0 * mi_1182[k];

        t_959[k] = -ki_763[k]
                   + f_0 * mi_1183[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, t_963, t_964, ki_764, ki_765, ki_766, ki_767, \
                         ki_768, mi_1184, mi_1185, mi_1186, mi_1187, \
                         mi_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = -ki_764[k]
                   + f_0 * mi_1184[k];

        t_961[k] = -ki_765[k]
                   + f_0 * mi_1185[k];

        t_962[k] = -ki_766[k]
                   + f_0 * mi_1186[k];

        t_963[k] = -ki_767[k]
                   + f_0 * mi_1187[k];

        t_964[k] = -ki_768[k]
                   + f_0 * mi_1188[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_1_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);
    auto *t_970 = buffer.data(target + 970);
    auto *t_971 = buffer.data(target + 971);
    auto *t_972 = buffer.data(target + 972);
    auto *t_973 = buffer.data(target + 973);
    auto *t_974 = buffer.data(target + 974);
    auto *t_975 = buffer.data(target + 975);
    auto *t_976 = buffer.data(target + 976);
    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);
    auto *t_983 = buffer.data(target + 983);
    auto *t_984 = buffer.data(target + 984);
    auto *t_985 = buffer.data(target + 985);
    auto *t_986 = buffer.data(target + 986);
    auto *t_987 = buffer.data(target + 987);
    auto *t_988 = buffer.data(target + 988);
    auto *t_989 = buffer.data(target + 989);
    auto *t_990 = buffer.data(target + 990);
    auto *t_991 = buffer.data(target + 991);
    auto *t_992 = buffer.data(target + 992);
    auto *t_993 = buffer.data(target + 993);
    auto *t_994 = buffer.data(target + 994);
    auto *t_995 = buffer.data(target + 995);
    auto *t_996 = buffer.data(target + 996);
    auto *t_997 = buffer.data(target + 997);
    auto *t_998 = buffer.data(target + 998);
    auto *t_999 = buffer.data(target + 999);
    auto *t_1000 = buffer.data(target + 1000);
    auto *t_1001 = buffer.data(target + 1001);
    auto *t_1002 = buffer.data(target + 1002);
    auto *t_1003 = buffer.data(target + 1003);
    auto *t_1004 = buffer.data(target + 1004);
    auto *t_1005 = buffer.data(target + 1005);
    auto *t_1006 = buffer.data(target + 1006);
    auto *t_1007 = buffer.data(target + 1007);
    auto *t_1008 = buffer.data(target + 1008);
    auto *t_1009 = buffer.data(target + 1009);
    auto *t_1010 = buffer.data(target + 1010);
    auto *t_1011 = buffer.data(target + 1011);
    auto *t_1012 = buffer.data(target + 1012);
    auto *t_1013 = buffer.data(target + 1013);
    auto *t_1014 = buffer.data(target + 1014);
    auto *t_1015 = buffer.data(target + 1015);
    auto *t_1016 = buffer.data(target + 1016);
    auto *t_1017 = buffer.data(target + 1017);
    auto *t_1018 = buffer.data(target + 1018);
    auto *t_1019 = buffer.data(target + 1019);
    auto *t_1020 = buffer.data(target + 1020);
    auto *t_1021 = buffer.data(target + 1021);
    auto *t_1022 = buffer.data(target + 1022);
    auto *t_1023 = buffer.data(target + 1023);
    auto *t_1024 = buffer.data(target + 1024);
    auto *t_1025 = buffer.data(target + 1025);
    auto *t_1026 = buffer.data(target + 1026);
    auto *t_1027 = buffer.data(target + 1027);
    auto *t_1028 = buffer.data(target + 1028);
    auto *t_1029 = buffer.data(target + 1029);
    auto *t_1030 = buffer.data(target + 1030);
    auto *t_1031 = buffer.data(target + 1031);
    auto *t_1032 = buffer.data(target + 1032);
    auto *t_1033 = buffer.data(target + 1033);
    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);
    auto *t_1050 = buffer.data(target + 1050);
    auto *t_1051 = buffer.data(target + 1051);
    auto *t_1052 = buffer.data(target + 1052);
    auto *t_1053 = buffer.data(target + 1053);
    auto *t_1054 = buffer.data(target + 1054);
    auto *t_1055 = buffer.data(target + 1055);
    auto *t_1056 = buffer.data(target + 1056);
    auto *t_1057 = buffer.data(target + 1057);
    auto *t_1058 = buffer.data(target + 1058);
    auto *t_1059 = buffer.data(target + 1059);
    auto *t_1060 = buffer.data(target + 1060);
    auto *t_1061 = buffer.data(target + 1061);
    auto *t_1062 = buffer.data(target + 1062);
    auto *t_1063 = buffer.data(target + 1063);
    auto *t_1064 = buffer.data(target + 1064);
    auto *t_1065 = buffer.data(target + 1065);
    auto *t_1066 = buffer.data(target + 1066);
    auto *t_1067 = buffer.data(target + 1067);
    auto *t_1068 = buffer.data(target + 1068);
    auto *t_1069 = buffer.data(target + 1069);
    auto *t_1070 = buffer.data(target + 1070);
    auto *t_1071 = buffer.data(target + 1071);
    auto *t_1072 = buffer.data(target + 1072);
    auto *t_1073 = buffer.data(target + 1073);
    auto *t_1074 = buffer.data(target + 1074);
    auto *t_1075 = buffer.data(target + 1075);
    auto *t_1076 = buffer.data(target + 1076);
    auto *t_1077 = buffer.data(target + 1077);
    auto *t_1078 = buffer.data(target + 1078);
    auto *t_1079 = buffer.data(target + 1079);
    auto *t_1080 = buffer.data(target + 1080);
    auto *t_1081 = buffer.data(target + 1081);
    auto *t_1082 = buffer.data(target + 1082);
    auto *t_1083 = buffer.data(target + 1083);
    auto *t_1084 = buffer.data(target + 1084);
    auto *t_1085 = buffer.data(target + 1085);
    auto *t_1086 = buffer.data(target + 1086);
    auto *t_1087 = buffer.data(target + 1087);
    auto *t_1088 = buffer.data(target + 1088);
    auto *t_1089 = buffer.data(target + 1089);
    auto *t_1090 = buffer.data(target + 1090);
    auto *t_1091 = buffer.data(target + 1091);
    auto *t_1092 = buffer.data(target + 1092);
    auto *t_1093 = buffer.data(target + 1093);
    auto *t_1094 = buffer.data(target + 1094);
    auto *t_1095 = buffer.data(target + 1095);
    auto *t_1096 = buffer.data(target + 1096);
    auto *t_1097 = buffer.data(target + 1097);
    auto *t_1098 = buffer.data(target + 1098);
    auto *t_1099 = buffer.data(target + 1099);
    auto *t_1100 = buffer.data(target + 1100);
    auto *t_1101 = buffer.data(target + 1101);
    auto *t_1102 = buffer.data(target + 1102);
    auto *t_1103 = buffer.data(target + 1103);
    auto *t_1104 = buffer.data(target + 1104);
    auto *t_1105 = buffer.data(target + 1105);
    auto *t_1106 = buffer.data(target + 1106);
    auto *t_1107 = buffer.data(target + 1107);
    auto *t_1108 = buffer.data(target + 1108);
    auto *t_1109 = buffer.data(target + 1109);
    auto *t_1110 = buffer.data(target + 1110);
    auto *t_1111 = buffer.data(target + 1111);
    auto *t_1112 = buffer.data(target + 1112);
    auto *t_1113 = buffer.data(target + 1113);
    auto *t_1114 = buffer.data(target + 1114);
    auto *t_1115 = buffer.data(target + 1115);
    auto *t_1116 = buffer.data(target + 1116);
    auto *t_1117 = buffer.data(target + 1117);
    auto *t_1118 = buffer.data(target + 1118);
    auto *t_1119 = buffer.data(target + 1119);
    auto *t_1120 = buffer.data(target + 1120);
    auto *t_1121 = buffer.data(target + 1121);
    auto *t_1122 = buffer.data(target + 1122);
    auto *t_1123 = buffer.data(target + 1123);
    auto *t_1124 = buffer.data(target + 1124);

    const auto *ki_769 = buffer.data(ki + 769);
    const auto *ki_770 = buffer.data(ki + 770);
    const auto *ki_771 = buffer.data(ki + 771);
    const auto *ki_772 = buffer.data(ki + 772);
    const auto *ki_773 = buffer.data(ki + 773);
    const auto *ki_774 = buffer.data(ki + 774);
    const auto *ki_775 = buffer.data(ki + 775);
    const auto *ki_776 = buffer.data(ki + 776);
    const auto *ki_777 = buffer.data(ki + 777);
    const auto *ki_778 = buffer.data(ki + 778);
    const auto *ki_779 = buffer.data(ki + 779);
    const auto *ki_780 = buffer.data(ki + 780);
    const auto *ki_781 = buffer.data(ki + 781);
    const auto *ki_782 = buffer.data(ki + 782);
    const auto *ki_783 = buffer.data(ki + 783);
    const auto *ki_784 = buffer.data(ki + 784);
    const auto *ki_785 = buffer.data(ki + 785);
    const auto *ki_786 = buffer.data(ki + 786);
    const auto *ki_787 = buffer.data(ki + 787);
    const auto *ki_788 = buffer.data(ki + 788);
    const auto *ki_789 = buffer.data(ki + 789);
    const auto *ki_790 = buffer.data(ki + 790);
    const auto *ki_791 = buffer.data(ki + 791);
    const auto *ki_792 = buffer.data(ki + 792);
    const auto *ki_793 = buffer.data(ki + 793);
    const auto *ki_794 = buffer.data(ki + 794);
    const auto *ki_795 = buffer.data(ki + 795);
    const auto *ki_796 = buffer.data(ki + 796);
    const auto *ki_797 = buffer.data(ki + 797);
    const auto *ki_798 = buffer.data(ki + 798);
    const auto *ki_799 = buffer.data(ki + 799);
    const auto *ki_800 = buffer.data(ki + 800);
    const auto *ki_801 = buffer.data(ki + 801);
    const auto *ki_802 = buffer.data(ki + 802);
    const auto *ki_803 = buffer.data(ki + 803);
    const auto *ki_804 = buffer.data(ki + 804);
    const auto *ki_805 = buffer.data(ki + 805);
    const auto *ki_806 = buffer.data(ki + 806);
    const auto *ki_807 = buffer.data(ki + 807);
    const auto *ki_808 = buffer.data(ki + 808);
    const auto *ki_809 = buffer.data(ki + 809);
    const auto *ki_810 = buffer.data(ki + 810);
    const auto *ki_811 = buffer.data(ki + 811);
    const auto *ki_812 = buffer.data(ki + 812);
    const auto *ki_813 = buffer.data(ki + 813);
    const auto *ki_814 = buffer.data(ki + 814);
    const auto *ki_815 = buffer.data(ki + 815);
    const auto *ki_816 = buffer.data(ki + 816);
    const auto *ki_817 = buffer.data(ki + 817);
    const auto *ki_818 = buffer.data(ki + 818);
    const auto *ki_819 = buffer.data(ki + 819);
    const auto *ki_820 = buffer.data(ki + 820);
    const auto *ki_821 = buffer.data(ki + 821);
    const auto *ki_822 = buffer.data(ki + 822);
    const auto *ki_823 = buffer.data(ki + 823);
    const auto *ki_824 = buffer.data(ki + 824);
    const auto *ki_825 = buffer.data(ki + 825);
    const auto *ki_826 = buffer.data(ki + 826);
    const auto *ki_827 = buffer.data(ki + 827);
    const auto *ki_828 = buffer.data(ki + 828);
    const auto *ki_829 = buffer.data(ki + 829);
    const auto *ki_830 = buffer.data(ki + 830);
    const auto *ki_831 = buffer.data(ki + 831);
    const auto *ki_832 = buffer.data(ki + 832);
    const auto *ki_833 = buffer.data(ki + 833);
    const auto *ki_834 = buffer.data(ki + 834);
    const auto *ki_835 = buffer.data(ki + 835);
    const auto *ki_836 = buffer.data(ki + 836);
    const auto *ki_837 = buffer.data(ki + 837);
    const auto *ki_838 = buffer.data(ki + 838);
    const auto *ki_839 = buffer.data(ki + 839);
    const auto *ki_840 = buffer.data(ki + 840);
    const auto *ki_841 = buffer.data(ki + 841);
    const auto *ki_842 = buffer.data(ki + 842);
    const auto *ki_843 = buffer.data(ki + 843);
    const auto *ki_844 = buffer.data(ki + 844);
    const auto *ki_845 = buffer.data(ki + 845);
    const auto *ki_846 = buffer.data(ki + 846);
    const auto *ki_847 = buffer.data(ki + 847);
    const auto *ki_848 = buffer.data(ki + 848);
    const auto *ki_849 = buffer.data(ki + 849);
    const auto *ki_850 = buffer.data(ki + 850);
    const auto *ki_851 = buffer.data(ki + 851);
    const auto *ki_852 = buffer.data(ki + 852);
    const auto *ki_853 = buffer.data(ki + 853);
    const auto *ki_854 = buffer.data(ki + 854);
    const auto *ki_855 = buffer.data(ki + 855);
    const auto *ki_856 = buffer.data(ki + 856);
    const auto *ki_857 = buffer.data(ki + 857);
    const auto *ki_858 = buffer.data(ki + 858);
    const auto *ki_859 = buffer.data(ki + 859);
    const auto *ki_860 = buffer.data(ki + 860);
    const auto *ki_861 = buffer.data(ki + 861);
    const auto *ki_862 = buffer.data(ki + 862);
    const auto *ki_863 = buffer.data(ki + 863);
    const auto *ki_864 = buffer.data(ki + 864);
    const auto *ki_865 = buffer.data(ki + 865);
    const auto *ki_866 = buffer.data(ki + 866);
    const auto *ki_867 = buffer.data(ki + 867);
    const auto *ki_868 = buffer.data(ki + 868);
    const auto *ki_869 = buffer.data(ki + 869);
    const auto *ki_870 = buffer.data(ki + 870);
    const auto *ki_871 = buffer.data(ki + 871);
    const auto *ki_872 = buffer.data(ki + 872);
    const auto *ki_873 = buffer.data(ki + 873);
    const auto *ki_874 = buffer.data(ki + 874);
    const auto *ki_875 = buffer.data(ki + 875);
    const auto *ki_876 = buffer.data(ki + 876);
    const auto *ki_877 = buffer.data(ki + 877);
    const auto *ki_878 = buffer.data(ki + 878);
    const auto *ki_879 = buffer.data(ki + 879);
    const auto *ki_880 = buffer.data(ki + 880);
    const auto *ki_881 = buffer.data(ki + 881);
    const auto *ki_882 = buffer.data(ki + 882);
    const auto *ki_883 = buffer.data(ki + 883);
    const auto *ki_884 = buffer.data(ki + 884);
    const auto *ki_885 = buffer.data(ki + 885);
    const auto *ki_886 = buffer.data(ki + 886);
    const auto *ki_887 = buffer.data(ki + 887);
    const auto *ki_888 = buffer.data(ki + 888);
    const auto *ki_889 = buffer.data(ki + 889);
    const auto *ki_890 = buffer.data(ki + 890);
    const auto *ki_891 = buffer.data(ki + 891);
    const auto *ki_892 = buffer.data(ki + 892);
    const auto *ki_893 = buffer.data(ki + 893);
    const auto *ki_894 = buffer.data(ki + 894);
    const auto *ki_895 = buffer.data(ki + 895);
    const auto *ki_896 = buffer.data(ki + 896);
    const auto *ki_897 = buffer.data(ki + 897);
    const auto *ki_898 = buffer.data(ki + 898);
    const auto *ki_899 = buffer.data(ki + 899);
    const auto *ki_900 = buffer.data(ki + 900);

    const auto *mi_1189 = buffer.data(mi + 1189);
    const auto *mi_1190 = buffer.data(mi + 1190);
    const auto *mi_1191 = buffer.data(mi + 1191);
    const auto *mi_1192 = buffer.data(mi + 1192);
    const auto *mi_1193 = buffer.data(mi + 1193);
    const auto *mi_1194 = buffer.data(mi + 1194);
    const auto *mi_1195 = buffer.data(mi + 1195);
    const auto *mi_1196 = buffer.data(mi + 1196);
    const auto *mi_1197 = buffer.data(mi + 1197);
    const auto *mi_1198 = buffer.data(mi + 1198);
    const auto *mi_1199 = buffer.data(mi + 1199);
    const auto *mi_1200 = buffer.data(mi + 1200);
    const auto *mi_1201 = buffer.data(mi + 1201);
    const auto *mi_1202 = buffer.data(mi + 1202);
    const auto *mi_1203 = buffer.data(mi + 1203);
    const auto *mi_1204 = buffer.data(mi + 1204);
    const auto *mi_1205 = buffer.data(mi + 1205);
    const auto *mi_1206 = buffer.data(mi + 1206);
    const auto *mi_1207 = buffer.data(mi + 1207);
    const auto *mi_1208 = buffer.data(mi + 1208);
    const auto *mi_1209 = buffer.data(mi + 1209);
    const auto *mi_1210 = buffer.data(mi + 1210);
    const auto *mi_1211 = buffer.data(mi + 1211);
    const auto *mi_1212 = buffer.data(mi + 1212);
    const auto *mi_1213 = buffer.data(mi + 1213);
    const auto *mi_1214 = buffer.data(mi + 1214);
    const auto *mi_1215 = buffer.data(mi + 1215);
    const auto *mi_1216 = buffer.data(mi + 1216);
    const auto *mi_1217 = buffer.data(mi + 1217);
    const auto *mi_1218 = buffer.data(mi + 1218);
    const auto *mi_1219 = buffer.data(mi + 1219);
    const auto *mi_1220 = buffer.data(mi + 1220);
    const auto *mi_1221 = buffer.data(mi + 1221);
    const auto *mi_1222 = buffer.data(mi + 1222);
    const auto *mi_1223 = buffer.data(mi + 1223);
    const auto *mi_1224 = buffer.data(mi + 1224);
    const auto *mi_1225 = buffer.data(mi + 1225);
    const auto *mi_1226 = buffer.data(mi + 1226);
    const auto *mi_1227 = buffer.data(mi + 1227);
    const auto *mi_1228 = buffer.data(mi + 1228);
    const auto *mi_1229 = buffer.data(mi + 1229);
    const auto *mi_1230 = buffer.data(mi + 1230);
    const auto *mi_1231 = buffer.data(mi + 1231);
    const auto *mi_1260 = buffer.data(mi + 1260);
    const auto *mi_1261 = buffer.data(mi + 1261);
    const auto *mi_1262 = buffer.data(mi + 1262);
    const auto *mi_1263 = buffer.data(mi + 1263);
    const auto *mi_1264 = buffer.data(mi + 1264);
    const auto *mi_1265 = buffer.data(mi + 1265);
    const auto *mi_1266 = buffer.data(mi + 1266);
    const auto *mi_1267 = buffer.data(mi + 1267);
    const auto *mi_1268 = buffer.data(mi + 1268);
    const auto *mi_1269 = buffer.data(mi + 1269);
    const auto *mi_1270 = buffer.data(mi + 1270);
    const auto *mi_1271 = buffer.data(mi + 1271);
    const auto *mi_1272 = buffer.data(mi + 1272);
    const auto *mi_1273 = buffer.data(mi + 1273);
    const auto *mi_1274 = buffer.data(mi + 1274);
    const auto *mi_1275 = buffer.data(mi + 1275);
    const auto *mi_1276 = buffer.data(mi + 1276);
    const auto *mi_1277 = buffer.data(mi + 1277);
    const auto *mi_1278 = buffer.data(mi + 1278);
    const auto *mi_1279 = buffer.data(mi + 1279);
    const auto *mi_1280 = buffer.data(mi + 1280);
    const auto *mi_1281 = buffer.data(mi + 1281);
    const auto *mi_1282 = buffer.data(mi + 1282);
    const auto *mi_1283 = buffer.data(mi + 1283);
    const auto *mi_1284 = buffer.data(mi + 1284);
    const auto *mi_1285 = buffer.data(mi + 1285);
    const auto *mi_1286 = buffer.data(mi + 1286);
    const auto *mi_1287 = buffer.data(mi + 1287);
    const auto *mi_1288 = buffer.data(mi + 1288);
    const auto *mi_1289 = buffer.data(mi + 1289);
    const auto *mi_1290 = buffer.data(mi + 1290);
    const auto *mi_1291 = buffer.data(mi + 1291);
    const auto *mi_1292 = buffer.data(mi + 1292);
    const auto *mi_1293 = buffer.data(mi + 1293);
    const auto *mi_1294 = buffer.data(mi + 1294);
    const auto *mi_1295 = buffer.data(mi + 1295);
    const auto *mi_1296 = buffer.data(mi + 1296);
    const auto *mi_1297 = buffer.data(mi + 1297);
    const auto *mi_1298 = buffer.data(mi + 1298);
    const auto *mi_1299 = buffer.data(mi + 1299);
    const auto *mi_1300 = buffer.data(mi + 1300);
    const auto *mi_1301 = buffer.data(mi + 1301);
    const auto *mi_1302 = buffer.data(mi + 1302);
    const auto *mi_1303 = buffer.data(mi + 1303);
    const auto *mi_1304 = buffer.data(mi + 1304);
    const auto *mi_1305 = buffer.data(mi + 1305);
    const auto *mi_1306 = buffer.data(mi + 1306);
    const auto *mi_1307 = buffer.data(mi + 1307);
    const auto *mi_1308 = buffer.data(mi + 1308);
    const auto *mi_1309 = buffer.data(mi + 1309);
    const auto *mi_1310 = buffer.data(mi + 1310);
    const auto *mi_1311 = buffer.data(mi + 1311);
    const auto *mi_1312 = buffer.data(mi + 1312);
    const auto *mi_1313 = buffer.data(mi + 1313);
    const auto *mi_1314 = buffer.data(mi + 1314);
    const auto *mi_1315 = buffer.data(mi + 1315);
    const auto *mi_1316 = buffer.data(mi + 1316);
    const auto *mi_1317 = buffer.data(mi + 1317);
    const auto *mi_1318 = buffer.data(mi + 1318);
    const auto *mi_1319 = buffer.data(mi + 1319);
    const auto *mi_1320 = buffer.data(mi + 1320);
    const auto *mi_1321 = buffer.data(mi + 1321);
    const auto *mi_1322 = buffer.data(mi + 1322);
    const auto *mi_1323 = buffer.data(mi + 1323);
    const auto *mi_1324 = buffer.data(mi + 1324);
    const auto *mi_1325 = buffer.data(mi + 1325);
    const auto *mi_1326 = buffer.data(mi + 1326);
    const auto *mi_1327 = buffer.data(mi + 1327);
    const auto *mi_1328 = buffer.data(mi + 1328);
    const auto *mi_1329 = buffer.data(mi + 1329);
    const auto *mi_1330 = buffer.data(mi + 1330);
    const auto *mi_1331 = buffer.data(mi + 1331);
    const auto *mi_1332 = buffer.data(mi + 1332);
    const auto *mi_1333 = buffer.data(mi + 1333);
    const auto *mi_1334 = buffer.data(mi + 1334);
    const auto *mi_1335 = buffer.data(mi + 1335);
    const auto *mi_1336 = buffer.data(mi + 1336);
    const auto *mi_1337 = buffer.data(mi + 1337);
    const auto *mi_1338 = buffer.data(mi + 1338);
    const auto *mi_1339 = buffer.data(mi + 1339);
    const auto *mi_1340 = buffer.data(mi + 1340);
    const auto *mi_1341 = buffer.data(mi + 1341);
    const auto *mi_1342 = buffer.data(mi + 1342);
    const auto *mi_1343 = buffer.data(mi + 1343);
    const auto *mi_1344 = buffer.data(mi + 1344);
    const auto *mi_1345 = buffer.data(mi + 1345);
    const auto *mi_1346 = buffer.data(mi + 1346);
    const auto *mi_1347 = buffer.data(mi + 1347);
    const auto *mi_1348 = buffer.data(mi + 1348);
    const auto *mi_1349 = buffer.data(mi + 1349);
    const auto *mi_1350 = buffer.data(mi + 1350);
    const auto *mi_1351 = buffer.data(mi + 1351);
    const auto *mi_1352 = buffer.data(mi + 1352);
    const auto *mi_1353 = buffer.data(mi + 1353);
    const auto *mi_1354 = buffer.data(mi + 1354);
    const auto *mi_1355 = buffer.data(mi + 1355);
    const auto *mi_1356 = buffer.data(mi + 1356);
    const auto *mi_1357 = buffer.data(mi + 1357);
    const auto *mi_1358 = buffer.data(mi + 1358);
    const auto *mi_1359 = buffer.data(mi + 1359);
    const auto *mi_1360 = buffer.data(mi + 1360);
    const auto *mi_1361 = buffer.data(mi + 1361);
    const auto *mi_1362 = buffer.data(mi + 1362);
    const auto *mi_1363 = buffer.data(mi + 1363);
    const auto *mi_1364 = buffer.data(mi + 1364);
    const auto *mi_1365 = buffer.data(mi + 1365);
    const auto *mi_1366 = buffer.data(mi + 1366);
    const auto *mi_1367 = buffer.data(mi + 1367);
    const auto *mi_1368 = buffer.data(mi + 1368);
    const auto *mi_1369 = buffer.data(mi + 1369);
    const auto *mi_1370 = buffer.data(mi + 1370);
    const auto *mi_1371 = buffer.data(mi + 1371);
    const auto *mi_1372 = buffer.data(mi + 1372);
    const auto *mi_1373 = buffer.data(mi + 1373);
    const auto *mi_1374 = buffer.data(mi + 1374);
    const auto *mi_1375 = buffer.data(mi + 1375);
    const auto *mi_1376 = buffer.data(mi + 1376);

#pragma omp simd aligned(t_965, t_966, t_967, t_968, t_969, ki_769, ki_770, ki_771, ki_772, \
                         ki_773, mi_1189, mi_1190, mi_1191, mi_1192, \
                         mi_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = -ki_769[k]
                   + f_0 * mi_1189[k];

        t_966[k] = -ki_770[k]
                   + f_0 * mi_1190[k];

        t_967[k] = -ki_771[k]
                   + f_0 * mi_1191[k];

        t_968[k] = -ki_772[k]
                   + f_0 * mi_1192[k];

        t_969[k] = -ki_773[k]
                   + f_0 * mi_1193[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, t_974, ki_774, ki_775, ki_776, ki_777, \
                         ki_778, mi_1194, mi_1195, mi_1196, mi_1197, \
                         mi_1198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = -ki_774[k]
                   + f_0 * mi_1194[k];

        t_971[k] = -ki_775[k]
                   + f_0 * mi_1195[k];

        t_972[k] = -ki_776[k]
                   + f_0 * mi_1196[k];

        t_973[k] = -ki_777[k]
                   + f_0 * mi_1197[k];

        t_974[k] = -ki_778[k]
                   + f_0 * mi_1198[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, t_979, ki_779, ki_780, ki_781, ki_782, \
                         ki_783, mi_1199, mi_1200, mi_1201, mi_1202, \
                         mi_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = -ki_779[k]
                   + f_0 * mi_1199[k];

        t_976[k] = -ki_780[k]
                   + f_0 * mi_1200[k];

        t_977[k] = -ki_781[k]
                   + f_0 * mi_1201[k];

        t_978[k] = -ki_782[k]
                   + f_0 * mi_1202[k];

        t_979[k] = -ki_783[k]
                   + f_0 * mi_1203[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, t_984, t_985, t_986, t_987, mi_1204, \
                         mi_1205, mi_1206, mi_1207, mi_1208, mi_1209, mi_1210, \
                         mi_1211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = f_0 * mi_1204[k];

        t_981[k] = f_0 * mi_1205[k];

        t_982[k] = f_0 * mi_1206[k];

        t_983[k] = f_0 * mi_1207[k];

        t_984[k] = f_0 * mi_1208[k];

        t_985[k] = f_0 * mi_1209[k];

        t_986[k] = f_0 * mi_1210[k];

        t_987[k] = f_0 * mi_1211[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, t_991, t_992, t_993, t_994, t_995, mi_1212, \
                         mi_1213, mi_1214, mi_1215, mi_1216, mi_1217, mi_1218, \
                         mi_1219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_0 * mi_1212[k];

        t_989[k] = f_0 * mi_1213[k];

        t_990[k] = f_0 * mi_1214[k];

        t_991[k] = f_0 * mi_1215[k];

        t_992[k] = f_0 * mi_1216[k];

        t_993[k] = f_0 * mi_1217[k];

        t_994[k] = f_0 * mi_1218[k];

        t_995[k] = f_0 * mi_1219[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, t_1000, t_1001, t_1002, t_1003, mi_1220, \
                         mi_1221, mi_1222, mi_1223, mi_1224, mi_1225, mi_1226, \
                         mi_1227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = f_0 * mi_1220[k];

        t_997[k] = f_0 * mi_1221[k];

        t_998[k] = f_0 * mi_1222[k];

        t_999[k] = f_0 * mi_1223[k];

        t_1000[k] = f_0 * mi_1224[k];

        t_1001[k] = f_0 * mi_1225[k];

        t_1002[k] = f_0 * mi_1226[k];

        t_1003[k] = f_0 * mi_1227[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, t_1008, t_1009, ki_784, ki_785, \
                         mi_1228, mi_1229, mi_1230, mi_1231, mi_1260, \
                         mi_1261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_0 * mi_1228[k];

        t_1005[k] = f_0 * mi_1229[k];

        t_1006[k] = f_0 * mi_1230[k];

        t_1007[k] = f_0 * mi_1231[k];

        t_1008[k] = -8.0 * ki_784[k]
                    + f_0 * mi_1260[k];

        t_1009[k] = -8.0 * ki_785[k]
                    + f_0 * mi_1261[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, t_1013, t_1014, ki_786, ki_787, ki_788, \
                         ki_789, ki_790, mi_1262, mi_1263, mi_1264, mi_1265, \
                         mi_1266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = -8.0 * ki_786[k]
                    + f_0 * mi_1262[k];

        t_1011[k] = -8.0 * ki_787[k]
                    + f_0 * mi_1263[k];

        t_1012[k] = -8.0 * ki_788[k]
                    + f_0 * mi_1264[k];

        t_1013[k] = -8.0 * ki_789[k]
                    + f_0 * mi_1265[k];

        t_1014[k] = -8.0 * ki_790[k]
                    + f_0 * mi_1266[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, t_1018, t_1019, ki_791, ki_792, ki_793, \
                         ki_794, ki_795, mi_1267, mi_1268, mi_1269, mi_1270, \
                         mi_1271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = -8.0 * ki_791[k]
                    + f_0 * mi_1267[k];

        t_1016[k] = -8.0 * ki_792[k]
                    + f_0 * mi_1268[k];

        t_1017[k] = -8.0 * ki_793[k]
                    + f_0 * mi_1269[k];

        t_1018[k] = -8.0 * ki_794[k]
                    + f_0 * mi_1270[k];

        t_1019[k] = -8.0 * ki_795[k]
                    + f_0 * mi_1271[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, t_1024, ki_796, ki_797, ki_798, \
                         ki_799, ki_800, mi_1272, mi_1273, mi_1274, mi_1275, \
                         mi_1276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = -8.0 * ki_796[k]
                    + f_0 * mi_1272[k];

        t_1021[k] = -8.0 * ki_797[k]
                    + f_0 * mi_1273[k];

        t_1022[k] = -8.0 * ki_798[k]
                    + f_0 * mi_1274[k];

        t_1023[k] = -8.0 * ki_799[k]
                    + f_0 * mi_1275[k];

        t_1024[k] = -8.0 * ki_800[k]
                    + f_0 * mi_1276[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, t_1029, ki_801, ki_802, ki_803, \
                         ki_804, ki_805, mi_1277, mi_1278, mi_1279, mi_1280, \
                         mi_1281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = -8.0 * ki_801[k]
                    + f_0 * mi_1277[k];

        t_1026[k] = -8.0 * ki_802[k]
                    + f_0 * mi_1278[k];

        t_1027[k] = -8.0 * ki_803[k]
                    + f_0 * mi_1279[k];

        t_1028[k] = -8.0 * ki_804[k]
                    + f_0 * mi_1280[k];

        t_1029[k] = -8.0 * ki_805[k]
                    + f_0 * mi_1281[k];
    }

#pragma omp simd aligned(t_1030, t_1031, t_1032, t_1033, t_1034, ki_806, ki_807, ki_808, \
                         ki_809, ki_810, mi_1282, mi_1283, mi_1284, mi_1285, \
                         mi_1286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1030[k] = -8.0 * ki_806[k]
                    + f_0 * mi_1282[k];

        t_1031[k] = -8.0 * ki_807[k]
                    + f_0 * mi_1283[k];

        t_1032[k] = -8.0 * ki_808[k]
                    + f_0 * mi_1284[k];

        t_1033[k] = -8.0 * ki_809[k]
                    + f_0 * mi_1285[k];

        t_1034[k] = -8.0 * ki_810[k]
                    + f_0 * mi_1286[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, t_1038, t_1039, ki_811, ki_812, ki_813, \
                         ki_814, ki_815, mi_1287, mi_1288, mi_1289, mi_1290, \
                         mi_1291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = -8.0 * ki_811[k]
                    + f_0 * mi_1287[k];

        t_1036[k] = -7.0 * ki_812[k]
                    + f_0 * mi_1288[k];

        t_1037[k] = -7.0 * ki_813[k]
                    + f_0 * mi_1289[k];

        t_1038[k] = -7.0 * ki_814[k]
                    + f_0 * mi_1290[k];

        t_1039[k] = -7.0 * ki_815[k]
                    + f_0 * mi_1291[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, t_1043, t_1044, ki_816, ki_817, ki_818, \
                         ki_819, ki_820, mi_1292, mi_1293, mi_1294, mi_1295, \
                         mi_1296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = -7.0 * ki_816[k]
                    + f_0 * mi_1292[k];

        t_1041[k] = -7.0 * ki_817[k]
                    + f_0 * mi_1293[k];

        t_1042[k] = -7.0 * ki_818[k]
                    + f_0 * mi_1294[k];

        t_1043[k] = -7.0 * ki_819[k]
                    + f_0 * mi_1295[k];

        t_1044[k] = -7.0 * ki_820[k]
                    + f_0 * mi_1296[k];
    }

#pragma omp simd aligned(t_1045, t_1046, t_1047, t_1048, t_1049, ki_821, ki_822, ki_823, \
                         ki_824, ki_825, mi_1297, mi_1298, mi_1299, mi_1300, \
                         mi_1301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1045[k] = -7.0 * ki_821[k]
                    + f_0 * mi_1297[k];

        t_1046[k] = -7.0 * ki_822[k]
                    + f_0 * mi_1298[k];

        t_1047[k] = -7.0 * ki_823[k]
                    + f_0 * mi_1299[k];

        t_1048[k] = -7.0 * ki_824[k]
                    + f_0 * mi_1300[k];

        t_1049[k] = -7.0 * ki_825[k]
                    + f_0 * mi_1301[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, t_1054, ki_826, ki_827, ki_828, \
                         ki_829, ki_830, mi_1302, mi_1303, mi_1304, mi_1305, \
                         mi_1306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = -7.0 * ki_826[k]
                    + f_0 * mi_1302[k];

        t_1051[k] = -7.0 * ki_827[k]
                    + f_0 * mi_1303[k];

        t_1052[k] = -7.0 * ki_828[k]
                    + f_0 * mi_1304[k];

        t_1053[k] = -7.0 * ki_829[k]
                    + f_0 * mi_1305[k];

        t_1054[k] = -7.0 * ki_830[k]
                    + f_0 * mi_1306[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, t_1058, t_1059, ki_831, ki_832, ki_833, \
                         ki_834, ki_835, mi_1307, mi_1308, mi_1309, mi_1310, \
                         mi_1311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = -7.0 * ki_831[k]
                    + f_0 * mi_1307[k];

        t_1056[k] = -7.0 * ki_832[k]
                    + f_0 * mi_1308[k];

        t_1057[k] = -7.0 * ki_833[k]
                    + f_0 * mi_1309[k];

        t_1058[k] = -7.0 * ki_834[k]
                    + f_0 * mi_1310[k];

        t_1059[k] = -7.0 * ki_835[k]
                    + f_0 * mi_1311[k];
    }

#pragma omp simd aligned(t_1060, t_1061, t_1062, t_1063, t_1064, ki_836, ki_837, ki_838, \
                         ki_839, ki_840, mi_1312, mi_1313, mi_1314, mi_1315, \
                         mi_1316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1060[k] = -7.0 * ki_836[k]
                    + f_0 * mi_1312[k];

        t_1061[k] = -7.0 * ki_837[k]
                    + f_0 * mi_1313[k];

        t_1062[k] = -7.0 * ki_838[k]
                    + f_0 * mi_1314[k];

        t_1063[k] = -7.0 * ki_839[k]
                    + f_0 * mi_1315[k];

        t_1064[k] = -6.0 * ki_840[k]
                    + f_0 * mi_1316[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, t_1069, ki_841, ki_842, ki_843, \
                         ki_844, ki_845, mi_1317, mi_1318, mi_1319, mi_1320, \
                         mi_1321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = -6.0 * ki_841[k]
                    + f_0 * mi_1317[k];

        t_1066[k] = -6.0 * ki_842[k]
                    + f_0 * mi_1318[k];

        t_1067[k] = -6.0 * ki_843[k]
                    + f_0 * mi_1319[k];

        t_1068[k] = -6.0 * ki_844[k]
                    + f_0 * mi_1320[k];

        t_1069[k] = -6.0 * ki_845[k]
                    + f_0 * mi_1321[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, t_1073, t_1074, ki_846, ki_847, ki_848, \
                         ki_849, ki_850, mi_1322, mi_1323, mi_1324, mi_1325, \
                         mi_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = -6.0 * ki_846[k]
                    + f_0 * mi_1322[k];

        t_1071[k] = -6.0 * ki_847[k]
                    + f_0 * mi_1323[k];

        t_1072[k] = -6.0 * ki_848[k]
                    + f_0 * mi_1324[k];

        t_1073[k] = -6.0 * ki_849[k]
                    + f_0 * mi_1325[k];

        t_1074[k] = -6.0 * ki_850[k]
                    + f_0 * mi_1326[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, t_1078, t_1079, ki_851, ki_852, ki_853, \
                         ki_854, ki_855, mi_1327, mi_1328, mi_1329, mi_1330, \
                         mi_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = -6.0 * ki_851[k]
                    + f_0 * mi_1327[k];

        t_1076[k] = -6.0 * ki_852[k]
                    + f_0 * mi_1328[k];

        t_1077[k] = -6.0 * ki_853[k]
                    + f_0 * mi_1329[k];

        t_1078[k] = -6.0 * ki_854[k]
                    + f_0 * mi_1330[k];

        t_1079[k] = -6.0 * ki_855[k]
                    + f_0 * mi_1331[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, t_1083, t_1084, ki_856, ki_857, ki_858, \
                         ki_859, ki_860, mi_1332, mi_1333, mi_1334, mi_1335, \
                         mi_1336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = -6.0 * ki_856[k]
                    + f_0 * mi_1332[k];

        t_1081[k] = -6.0 * ki_857[k]
                    + f_0 * mi_1333[k];

        t_1082[k] = -6.0 * ki_858[k]
                    + f_0 * mi_1334[k];

        t_1083[k] = -6.0 * ki_859[k]
                    + f_0 * mi_1335[k];

        t_1084[k] = -6.0 * ki_860[k]
                    + f_0 * mi_1336[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, t_1088, t_1089, ki_861, ki_862, ki_863, \
                         ki_864, ki_865, mi_1337, mi_1338, mi_1339, mi_1340, \
                         mi_1341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = -6.0 * ki_861[k]
                    + f_0 * mi_1337[k];

        t_1086[k] = -6.0 * ki_862[k]
                    + f_0 * mi_1338[k];

        t_1087[k] = -6.0 * ki_863[k]
                    + f_0 * mi_1339[k];

        t_1088[k] = -6.0 * ki_864[k]
                    + f_0 * mi_1340[k];

        t_1089[k] = -6.0 * ki_865[k]
                    + f_0 * mi_1341[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, t_1093, t_1094, ki_866, ki_867, ki_868, \
                         ki_869, ki_870, mi_1342, mi_1343, mi_1344, mi_1345, \
                         mi_1346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = -6.0 * ki_866[k]
                    + f_0 * mi_1342[k];

        t_1091[k] = -6.0 * ki_867[k]
                    + f_0 * mi_1343[k];

        t_1092[k] = -5.0 * ki_868[k]
                    + f_0 * mi_1344[k];

        t_1093[k] = -5.0 * ki_869[k]
                    + f_0 * mi_1345[k];

        t_1094[k] = -5.0 * ki_870[k]
                    + f_0 * mi_1346[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, t_1098, t_1099, ki_871, ki_872, ki_873, \
                         ki_874, ki_875, mi_1347, mi_1348, mi_1349, mi_1350, \
                         mi_1351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = -5.0 * ki_871[k]
                    + f_0 * mi_1347[k];

        t_1096[k] = -5.0 * ki_872[k]
                    + f_0 * mi_1348[k];

        t_1097[k] = -5.0 * ki_873[k]
                    + f_0 * mi_1349[k];

        t_1098[k] = -5.0 * ki_874[k]
                    + f_0 * mi_1350[k];

        t_1099[k] = -5.0 * ki_875[k]
                    + f_0 * mi_1351[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, t_1104, ki_876, ki_877, ki_878, \
                         ki_879, ki_880, mi_1352, mi_1353, mi_1354, mi_1355, \
                         mi_1356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = -5.0 * ki_876[k]
                    + f_0 * mi_1352[k];

        t_1101[k] = -5.0 * ki_877[k]
                    + f_0 * mi_1353[k];

        t_1102[k] = -5.0 * ki_878[k]
                    + f_0 * mi_1354[k];

        t_1103[k] = -5.0 * ki_879[k]
                    + f_0 * mi_1355[k];

        t_1104[k] = -5.0 * ki_880[k]
                    + f_0 * mi_1356[k];
    }

#pragma omp simd aligned(t_1105, t_1106, t_1107, t_1108, t_1109, ki_881, ki_882, ki_883, \
                         ki_884, ki_885, mi_1357, mi_1358, mi_1359, mi_1360, \
                         mi_1361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1105[k] = -5.0 * ki_881[k]
                    + f_0 * mi_1357[k];

        t_1106[k] = -5.0 * ki_882[k]
                    + f_0 * mi_1358[k];

        t_1107[k] = -5.0 * ki_883[k]
                    + f_0 * mi_1359[k];

        t_1108[k] = -5.0 * ki_884[k]
                    + f_0 * mi_1360[k];

        t_1109[k] = -5.0 * ki_885[k]
                    + f_0 * mi_1361[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, t_1114, ki_886, ki_887, ki_888, \
                         ki_889, ki_890, mi_1362, mi_1363, mi_1364, mi_1365, \
                         mi_1366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = -5.0 * ki_886[k]
                    + f_0 * mi_1362[k];

        t_1111[k] = -5.0 * ki_887[k]
                    + f_0 * mi_1363[k];

        t_1112[k] = -5.0 * ki_888[k]
                    + f_0 * mi_1364[k];

        t_1113[k] = -5.0 * ki_889[k]
                    + f_0 * mi_1365[k];

        t_1114[k] = -5.0 * ki_890[k]
                    + f_0 * mi_1366[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, t_1118, t_1119, ki_891, ki_892, ki_893, \
                         ki_894, ki_895, mi_1367, mi_1368, mi_1369, mi_1370, \
                         mi_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = -5.0 * ki_891[k]
                    + f_0 * mi_1367[k];

        t_1116[k] = -5.0 * ki_892[k]
                    + f_0 * mi_1368[k];

        t_1117[k] = -5.0 * ki_893[k]
                    + f_0 * mi_1369[k];

        t_1118[k] = -5.0 * ki_894[k]
                    + f_0 * mi_1370[k];

        t_1119[k] = -5.0 * ki_895[k]
                    + f_0 * mi_1371[k];
    }

#pragma omp simd aligned(t_1120, t_1121, t_1122, t_1123, t_1124, ki_896, ki_897, ki_898, \
                         ki_899, ki_900, mi_1372, mi_1373, mi_1374, mi_1375, \
                         mi_1376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1120[k] = -4.0 * ki_896[k]
                    + f_0 * mi_1372[k];

        t_1121[k] = -4.0 * ki_897[k]
                    + f_0 * mi_1373[k];

        t_1122[k] = -4.0 * ki_898[k]
                    + f_0 * mi_1374[k];

        t_1123[k] = -4.0 * ki_899[k]
                    + f_0 * mi_1375[k];

        t_1124[k] = -4.0 * ki_900[k]
                    + f_0 * mi_1376[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_1_piece7(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);
    auto *t_1133 = buffer.data(target + 1133);
    auto *t_1134 = buffer.data(target + 1134);
    auto *t_1135 = buffer.data(target + 1135);
    auto *t_1136 = buffer.data(target + 1136);
    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);
    auto *t_1169 = buffer.data(target + 1169);
    auto *t_1170 = buffer.data(target + 1170);
    auto *t_1171 = buffer.data(target + 1171);
    auto *t_1172 = buffer.data(target + 1172);
    auto *t_1173 = buffer.data(target + 1173);
    auto *t_1174 = buffer.data(target + 1174);
    auto *t_1175 = buffer.data(target + 1175);
    auto *t_1176 = buffer.data(target + 1176);
    auto *t_1177 = buffer.data(target + 1177);
    auto *t_1178 = buffer.data(target + 1178);
    auto *t_1179 = buffer.data(target + 1179);
    auto *t_1180 = buffer.data(target + 1180);
    auto *t_1181 = buffer.data(target + 1181);
    auto *t_1182 = buffer.data(target + 1182);
    auto *t_1183 = buffer.data(target + 1183);
    auto *t_1184 = buffer.data(target + 1184);
    auto *t_1185 = buffer.data(target + 1185);
    auto *t_1186 = buffer.data(target + 1186);
    auto *t_1187 = buffer.data(target + 1187);
    auto *t_1188 = buffer.data(target + 1188);
    auto *t_1189 = buffer.data(target + 1189);
    auto *t_1190 = buffer.data(target + 1190);
    auto *t_1191 = buffer.data(target + 1191);
    auto *t_1192 = buffer.data(target + 1192);
    auto *t_1193 = buffer.data(target + 1193);
    auto *t_1194 = buffer.data(target + 1194);
    auto *t_1195 = buffer.data(target + 1195);
    auto *t_1196 = buffer.data(target + 1196);
    auto *t_1197 = buffer.data(target + 1197);
    auto *t_1198 = buffer.data(target + 1198);
    auto *t_1199 = buffer.data(target + 1199);
    auto *t_1200 = buffer.data(target + 1200);
    auto *t_1201 = buffer.data(target + 1201);
    auto *t_1202 = buffer.data(target + 1202);
    auto *t_1203 = buffer.data(target + 1203);
    auto *t_1204 = buffer.data(target + 1204);
    auto *t_1205 = buffer.data(target + 1205);
    auto *t_1206 = buffer.data(target + 1206);
    auto *t_1207 = buffer.data(target + 1207);
    auto *t_1208 = buffer.data(target + 1208);
    auto *t_1209 = buffer.data(target + 1209);
    auto *t_1210 = buffer.data(target + 1210);
    auto *t_1211 = buffer.data(target + 1211);
    auto *t_1212 = buffer.data(target + 1212);
    auto *t_1213 = buffer.data(target + 1213);
    auto *t_1214 = buffer.data(target + 1214);
    auto *t_1215 = buffer.data(target + 1215);
    auto *t_1216 = buffer.data(target + 1216);
    auto *t_1217 = buffer.data(target + 1217);
    auto *t_1218 = buffer.data(target + 1218);
    auto *t_1219 = buffer.data(target + 1219);
    auto *t_1220 = buffer.data(target + 1220);
    auto *t_1221 = buffer.data(target + 1221);
    auto *t_1222 = buffer.data(target + 1222);
    auto *t_1223 = buffer.data(target + 1223);
    auto *t_1224 = buffer.data(target + 1224);
    auto *t_1225 = buffer.data(target + 1225);
    auto *t_1226 = buffer.data(target + 1226);
    auto *t_1227 = buffer.data(target + 1227);
    auto *t_1228 = buffer.data(target + 1228);
    auto *t_1229 = buffer.data(target + 1229);
    auto *t_1230 = buffer.data(target + 1230);
    auto *t_1231 = buffer.data(target + 1231);
    auto *t_1232 = buffer.data(target + 1232);
    auto *t_1233 = buffer.data(target + 1233);
    auto *t_1234 = buffer.data(target + 1234);
    auto *t_1235 = buffer.data(target + 1235);
    auto *t_1236 = buffer.data(target + 1236);
    auto *t_1237 = buffer.data(target + 1237);
    auto *t_1238 = buffer.data(target + 1238);
    auto *t_1239 = buffer.data(target + 1239);
    auto *t_1240 = buffer.data(target + 1240);
    auto *t_1241 = buffer.data(target + 1241);
    auto *t_1242 = buffer.data(target + 1242);
    auto *t_1243 = buffer.data(target + 1243);
    auto *t_1244 = buffer.data(target + 1244);
    auto *t_1245 = buffer.data(target + 1245);
    auto *t_1246 = buffer.data(target + 1246);
    auto *t_1247 = buffer.data(target + 1247);
    auto *t_1248 = buffer.data(target + 1248);
    auto *t_1249 = buffer.data(target + 1249);
    auto *t_1250 = buffer.data(target + 1250);
    auto *t_1251 = buffer.data(target + 1251);
    auto *t_1252 = buffer.data(target + 1252);
    auto *t_1253 = buffer.data(target + 1253);
    auto *t_1254 = buffer.data(target + 1254);
    auto *t_1255 = buffer.data(target + 1255);
    auto *t_1256 = buffer.data(target + 1256);
    auto *t_1257 = buffer.data(target + 1257);
    auto *t_1258 = buffer.data(target + 1258);
    auto *t_1259 = buffer.data(target + 1259);

    const auto *ki_901 = buffer.data(ki + 901);
    const auto *ki_902 = buffer.data(ki + 902);
    const auto *ki_903 = buffer.data(ki + 903);
    const auto *ki_904 = buffer.data(ki + 904);
    const auto *ki_905 = buffer.data(ki + 905);
    const auto *ki_906 = buffer.data(ki + 906);
    const auto *ki_907 = buffer.data(ki + 907);
    const auto *ki_908 = buffer.data(ki + 908);
    const auto *ki_909 = buffer.data(ki + 909);
    const auto *ki_910 = buffer.data(ki + 910);
    const auto *ki_911 = buffer.data(ki + 911);
    const auto *ki_912 = buffer.data(ki + 912);
    const auto *ki_913 = buffer.data(ki + 913);
    const auto *ki_914 = buffer.data(ki + 914);
    const auto *ki_915 = buffer.data(ki + 915);
    const auto *ki_916 = buffer.data(ki + 916);
    const auto *ki_917 = buffer.data(ki + 917);
    const auto *ki_918 = buffer.data(ki + 918);
    const auto *ki_919 = buffer.data(ki + 919);
    const auto *ki_920 = buffer.data(ki + 920);
    const auto *ki_921 = buffer.data(ki + 921);
    const auto *ki_922 = buffer.data(ki + 922);
    const auto *ki_923 = buffer.data(ki + 923);
    const auto *ki_924 = buffer.data(ki + 924);
    const auto *ki_925 = buffer.data(ki + 925);
    const auto *ki_926 = buffer.data(ki + 926);
    const auto *ki_927 = buffer.data(ki + 927);
    const auto *ki_928 = buffer.data(ki + 928);
    const auto *ki_929 = buffer.data(ki + 929);
    const auto *ki_930 = buffer.data(ki + 930);
    const auto *ki_931 = buffer.data(ki + 931);
    const auto *ki_932 = buffer.data(ki + 932);
    const auto *ki_933 = buffer.data(ki + 933);
    const auto *ki_934 = buffer.data(ki + 934);
    const auto *ki_935 = buffer.data(ki + 935);
    const auto *ki_936 = buffer.data(ki + 936);
    const auto *ki_937 = buffer.data(ki + 937);
    const auto *ki_938 = buffer.data(ki + 938);
    const auto *ki_939 = buffer.data(ki + 939);
    const auto *ki_940 = buffer.data(ki + 940);
    const auto *ki_941 = buffer.data(ki + 941);
    const auto *ki_942 = buffer.data(ki + 942);
    const auto *ki_943 = buffer.data(ki + 943);
    const auto *ki_944 = buffer.data(ki + 944);
    const auto *ki_945 = buffer.data(ki + 945);
    const auto *ki_946 = buffer.data(ki + 946);
    const auto *ki_947 = buffer.data(ki + 947);
    const auto *ki_948 = buffer.data(ki + 948);
    const auto *ki_949 = buffer.data(ki + 949);
    const auto *ki_950 = buffer.data(ki + 950);
    const auto *ki_951 = buffer.data(ki + 951);
    const auto *ki_952 = buffer.data(ki + 952);
    const auto *ki_953 = buffer.data(ki + 953);
    const auto *ki_954 = buffer.data(ki + 954);
    const auto *ki_955 = buffer.data(ki + 955);
    const auto *ki_956 = buffer.data(ki + 956);
    const auto *ki_957 = buffer.data(ki + 957);
    const auto *ki_958 = buffer.data(ki + 958);
    const auto *ki_959 = buffer.data(ki + 959);
    const auto *ki_960 = buffer.data(ki + 960);
    const auto *ki_961 = buffer.data(ki + 961);
    const auto *ki_962 = buffer.data(ki + 962);
    const auto *ki_963 = buffer.data(ki + 963);
    const auto *ki_964 = buffer.data(ki + 964);
    const auto *ki_965 = buffer.data(ki + 965);
    const auto *ki_966 = buffer.data(ki + 966);
    const auto *ki_967 = buffer.data(ki + 967);
    const auto *ki_968 = buffer.data(ki + 968);
    const auto *ki_969 = buffer.data(ki + 969);
    const auto *ki_970 = buffer.data(ki + 970);
    const auto *ki_971 = buffer.data(ki + 971);
    const auto *ki_972 = buffer.data(ki + 972);
    const auto *ki_973 = buffer.data(ki + 973);
    const auto *ki_974 = buffer.data(ki + 974);
    const auto *ki_975 = buffer.data(ki + 975);
    const auto *ki_976 = buffer.data(ki + 976);
    const auto *ki_977 = buffer.data(ki + 977);
    const auto *ki_978 = buffer.data(ki + 978);
    const auto *ki_979 = buffer.data(ki + 979);
    const auto *ki_980 = buffer.data(ki + 980);
    const auto *ki_981 = buffer.data(ki + 981);
    const auto *ki_982 = buffer.data(ki + 982);
    const auto *ki_983 = buffer.data(ki + 983);
    const auto *ki_984 = buffer.data(ki + 984);
    const auto *ki_985 = buffer.data(ki + 985);
    const auto *ki_986 = buffer.data(ki + 986);
    const auto *ki_987 = buffer.data(ki + 987);
    const auto *ki_988 = buffer.data(ki + 988);
    const auto *ki_989 = buffer.data(ki + 989);
    const auto *ki_990 = buffer.data(ki + 990);
    const auto *ki_991 = buffer.data(ki + 991);
    const auto *ki_992 = buffer.data(ki + 992);
    const auto *ki_993 = buffer.data(ki + 993);
    const auto *ki_994 = buffer.data(ki + 994);
    const auto *ki_995 = buffer.data(ki + 995);
    const auto *ki_996 = buffer.data(ki + 996);
    const auto *ki_997 = buffer.data(ki + 997);
    const auto *ki_998 = buffer.data(ki + 998);
    const auto *ki_999 = buffer.data(ki + 999);
    const auto *ki_1000 = buffer.data(ki + 1000);
    const auto *ki_1001 = buffer.data(ki + 1001);
    const auto *ki_1002 = buffer.data(ki + 1002);
    const auto *ki_1003 = buffer.data(ki + 1003);
    const auto *ki_1004 = buffer.data(ki + 1004);
    const auto *ki_1005 = buffer.data(ki + 1005);
    const auto *ki_1006 = buffer.data(ki + 1006);
    const auto *ki_1007 = buffer.data(ki + 1007);

    const auto *mi_1377 = buffer.data(mi + 1377);
    const auto *mi_1378 = buffer.data(mi + 1378);
    const auto *mi_1379 = buffer.data(mi + 1379);
    const auto *mi_1380 = buffer.data(mi + 1380);
    const auto *mi_1381 = buffer.data(mi + 1381);
    const auto *mi_1382 = buffer.data(mi + 1382);
    const auto *mi_1383 = buffer.data(mi + 1383);
    const auto *mi_1384 = buffer.data(mi + 1384);
    const auto *mi_1385 = buffer.data(mi + 1385);
    const auto *mi_1386 = buffer.data(mi + 1386);
    const auto *mi_1387 = buffer.data(mi + 1387);
    const auto *mi_1388 = buffer.data(mi + 1388);
    const auto *mi_1389 = buffer.data(mi + 1389);
    const auto *mi_1390 = buffer.data(mi + 1390);
    const auto *mi_1391 = buffer.data(mi + 1391);
    const auto *mi_1392 = buffer.data(mi + 1392);
    const auto *mi_1393 = buffer.data(mi + 1393);
    const auto *mi_1394 = buffer.data(mi + 1394);
    const auto *mi_1395 = buffer.data(mi + 1395);
    const auto *mi_1396 = buffer.data(mi + 1396);
    const auto *mi_1397 = buffer.data(mi + 1397);
    const auto *mi_1398 = buffer.data(mi + 1398);
    const auto *mi_1399 = buffer.data(mi + 1399);
    const auto *mi_1400 = buffer.data(mi + 1400);
    const auto *mi_1401 = buffer.data(mi + 1401);
    const auto *mi_1402 = buffer.data(mi + 1402);
    const auto *mi_1403 = buffer.data(mi + 1403);
    const auto *mi_1404 = buffer.data(mi + 1404);
    const auto *mi_1405 = buffer.data(mi + 1405);
    const auto *mi_1406 = buffer.data(mi + 1406);
    const auto *mi_1407 = buffer.data(mi + 1407);
    const auto *mi_1408 = buffer.data(mi + 1408);
    const auto *mi_1409 = buffer.data(mi + 1409);
    const auto *mi_1410 = buffer.data(mi + 1410);
    const auto *mi_1411 = buffer.data(mi + 1411);
    const auto *mi_1412 = buffer.data(mi + 1412);
    const auto *mi_1413 = buffer.data(mi + 1413);
    const auto *mi_1414 = buffer.data(mi + 1414);
    const auto *mi_1415 = buffer.data(mi + 1415);
    const auto *mi_1416 = buffer.data(mi + 1416);
    const auto *mi_1417 = buffer.data(mi + 1417);
    const auto *mi_1418 = buffer.data(mi + 1418);
    const auto *mi_1419 = buffer.data(mi + 1419);
    const auto *mi_1420 = buffer.data(mi + 1420);
    const auto *mi_1421 = buffer.data(mi + 1421);
    const auto *mi_1422 = buffer.data(mi + 1422);
    const auto *mi_1423 = buffer.data(mi + 1423);
    const auto *mi_1424 = buffer.data(mi + 1424);
    const auto *mi_1425 = buffer.data(mi + 1425);
    const auto *mi_1426 = buffer.data(mi + 1426);
    const auto *mi_1427 = buffer.data(mi + 1427);
    const auto *mi_1428 = buffer.data(mi + 1428);
    const auto *mi_1429 = buffer.data(mi + 1429);
    const auto *mi_1430 = buffer.data(mi + 1430);
    const auto *mi_1431 = buffer.data(mi + 1431);
    const auto *mi_1432 = buffer.data(mi + 1432);
    const auto *mi_1433 = buffer.data(mi + 1433);
    const auto *mi_1434 = buffer.data(mi + 1434);
    const auto *mi_1435 = buffer.data(mi + 1435);
    const auto *mi_1436 = buffer.data(mi + 1436);
    const auto *mi_1437 = buffer.data(mi + 1437);
    const auto *mi_1438 = buffer.data(mi + 1438);
    const auto *mi_1439 = buffer.data(mi + 1439);
    const auto *mi_1440 = buffer.data(mi + 1440);
    const auto *mi_1441 = buffer.data(mi + 1441);
    const auto *mi_1442 = buffer.data(mi + 1442);
    const auto *mi_1443 = buffer.data(mi + 1443);
    const auto *mi_1444 = buffer.data(mi + 1444);
    const auto *mi_1445 = buffer.data(mi + 1445);
    const auto *mi_1446 = buffer.data(mi + 1446);
    const auto *mi_1447 = buffer.data(mi + 1447);
    const auto *mi_1448 = buffer.data(mi + 1448);
    const auto *mi_1449 = buffer.data(mi + 1449);
    const auto *mi_1450 = buffer.data(mi + 1450);
    const auto *mi_1451 = buffer.data(mi + 1451);
    const auto *mi_1452 = buffer.data(mi + 1452);
    const auto *mi_1453 = buffer.data(mi + 1453);
    const auto *mi_1454 = buffer.data(mi + 1454);
    const auto *mi_1455 = buffer.data(mi + 1455);
    const auto *mi_1456 = buffer.data(mi + 1456);
    const auto *mi_1457 = buffer.data(mi + 1457);
    const auto *mi_1458 = buffer.data(mi + 1458);
    const auto *mi_1459 = buffer.data(mi + 1459);
    const auto *mi_1460 = buffer.data(mi + 1460);
    const auto *mi_1461 = buffer.data(mi + 1461);
    const auto *mi_1462 = buffer.data(mi + 1462);
    const auto *mi_1463 = buffer.data(mi + 1463);
    const auto *mi_1464 = buffer.data(mi + 1464);
    const auto *mi_1465 = buffer.data(mi + 1465);
    const auto *mi_1466 = buffer.data(mi + 1466);
    const auto *mi_1467 = buffer.data(mi + 1467);
    const auto *mi_1468 = buffer.data(mi + 1468);
    const auto *mi_1469 = buffer.data(mi + 1469);
    const auto *mi_1470 = buffer.data(mi + 1470);
    const auto *mi_1471 = buffer.data(mi + 1471);
    const auto *mi_1472 = buffer.data(mi + 1472);
    const auto *mi_1473 = buffer.data(mi + 1473);
    const auto *mi_1474 = buffer.data(mi + 1474);
    const auto *mi_1475 = buffer.data(mi + 1475);
    const auto *mi_1476 = buffer.data(mi + 1476);
    const auto *mi_1477 = buffer.data(mi + 1477);
    const auto *mi_1478 = buffer.data(mi + 1478);
    const auto *mi_1479 = buffer.data(mi + 1479);
    const auto *mi_1480 = buffer.data(mi + 1480);
    const auto *mi_1481 = buffer.data(mi + 1481);
    const auto *mi_1482 = buffer.data(mi + 1482);
    const auto *mi_1483 = buffer.data(mi + 1483);
    const auto *mi_1484 = buffer.data(mi + 1484);
    const auto *mi_1485 = buffer.data(mi + 1485);
    const auto *mi_1486 = buffer.data(mi + 1486);
    const auto *mi_1487 = buffer.data(mi + 1487);
    const auto *mi_1488 = buffer.data(mi + 1488);
    const auto *mi_1489 = buffer.data(mi + 1489);
    const auto *mi_1490 = buffer.data(mi + 1490);
    const auto *mi_1491 = buffer.data(mi + 1491);
    const auto *mi_1492 = buffer.data(mi + 1492);
    const auto *mi_1493 = buffer.data(mi + 1493);
    const auto *mi_1494 = buffer.data(mi + 1494);
    const auto *mi_1495 = buffer.data(mi + 1495);
    const auto *mi_1496 = buffer.data(mi + 1496);
    const auto *mi_1497 = buffer.data(mi + 1497);
    const auto *mi_1498 = buffer.data(mi + 1498);
    const auto *mi_1499 = buffer.data(mi + 1499);
    const auto *mi_1500 = buffer.data(mi + 1500);
    const auto *mi_1501 = buffer.data(mi + 1501);
    const auto *mi_1502 = buffer.data(mi + 1502);
    const auto *mi_1503 = buffer.data(mi + 1503);
    const auto *mi_1504 = buffer.data(mi + 1504);
    const auto *mi_1505 = buffer.data(mi + 1505);
    const auto *mi_1506 = buffer.data(mi + 1506);
    const auto *mi_1507 = buffer.data(mi + 1507);
    const auto *mi_1508 = buffer.data(mi + 1508);
    const auto *mi_1509 = buffer.data(mi + 1509);
    const auto *mi_1510 = buffer.data(mi + 1510);
    const auto *mi_1511 = buffer.data(mi + 1511);

#pragma omp simd aligned(t_1125, t_1126, t_1127, t_1128, t_1129, ki_901, ki_902, ki_903, \
                         ki_904, ki_905, mi_1377, mi_1378, mi_1379, mi_1380, \
                         mi_1381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = -4.0 * ki_901[k]
                    + f_0 * mi_1377[k];

        t_1126[k] = -4.0 * ki_902[k]
                    + f_0 * mi_1378[k];

        t_1127[k] = -4.0 * ki_903[k]
                    + f_0 * mi_1379[k];

        t_1128[k] = -4.0 * ki_904[k]
                    + f_0 * mi_1380[k];

        t_1129[k] = -4.0 * ki_905[k]
                    + f_0 * mi_1381[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, t_1133, t_1134, ki_906, ki_907, ki_908, \
                         ki_909, ki_910, mi_1382, mi_1383, mi_1384, mi_1385, \
                         mi_1386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = -4.0 * ki_906[k]
                    + f_0 * mi_1382[k];

        t_1131[k] = -4.0 * ki_907[k]
                    + f_0 * mi_1383[k];

        t_1132[k] = -4.0 * ki_908[k]
                    + f_0 * mi_1384[k];

        t_1133[k] = -4.0 * ki_909[k]
                    + f_0 * mi_1385[k];

        t_1134[k] = -4.0 * ki_910[k]
                    + f_0 * mi_1386[k];
    }

#pragma omp simd aligned(t_1135, t_1136, t_1137, t_1138, t_1139, ki_911, ki_912, ki_913, \
                         ki_914, ki_915, mi_1387, mi_1388, mi_1389, mi_1390, \
                         mi_1391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1135[k] = -4.0 * ki_911[k]
                    + f_0 * mi_1387[k];

        t_1136[k] = -4.0 * ki_912[k]
                    + f_0 * mi_1388[k];

        t_1137[k] = -4.0 * ki_913[k]
                    + f_0 * mi_1389[k];

        t_1138[k] = -4.0 * ki_914[k]
                    + f_0 * mi_1390[k];

        t_1139[k] = -4.0 * ki_915[k]
                    + f_0 * mi_1391[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, t_1143, t_1144, ki_916, ki_917, ki_918, \
                         ki_919, ki_920, mi_1392, mi_1393, mi_1394, mi_1395, \
                         mi_1396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = -4.0 * ki_916[k]
                    + f_0 * mi_1392[k];

        t_1141[k] = -4.0 * ki_917[k]
                    + f_0 * mi_1393[k];

        t_1142[k] = -4.0 * ki_918[k]
                    + f_0 * mi_1394[k];

        t_1143[k] = -4.0 * ki_919[k]
                    + f_0 * mi_1395[k];

        t_1144[k] = -4.0 * ki_920[k]
                    + f_0 * mi_1396[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, t_1148, t_1149, ki_921, ki_922, ki_923, \
                         ki_924, ki_925, mi_1397, mi_1398, mi_1399, mi_1400, \
                         mi_1401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = -4.0 * ki_921[k]
                    + f_0 * mi_1397[k];

        t_1146[k] = -4.0 * ki_922[k]
                    + f_0 * mi_1398[k];

        t_1147[k] = -4.0 * ki_923[k]
                    + f_0 * mi_1399[k];

        t_1148[k] = -3.0 * ki_924[k]
                    + f_0 * mi_1400[k];

        t_1149[k] = -3.0 * ki_925[k]
                    + f_0 * mi_1401[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, t_1153, t_1154, ki_926, ki_927, ki_928, \
                         ki_929, ki_930, mi_1402, mi_1403, mi_1404, mi_1405, \
                         mi_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = -3.0 * ki_926[k]
                    + f_0 * mi_1402[k];

        t_1151[k] = -3.0 * ki_927[k]
                    + f_0 * mi_1403[k];

        t_1152[k] = -3.0 * ki_928[k]
                    + f_0 * mi_1404[k];

        t_1153[k] = -3.0 * ki_929[k]
                    + f_0 * mi_1405[k];

        t_1154[k] = -3.0 * ki_930[k]
                    + f_0 * mi_1406[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, t_1159, ki_931, ki_932, ki_933, \
                         ki_934, ki_935, mi_1407, mi_1408, mi_1409, mi_1410, \
                         mi_1411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = -3.0 * ki_931[k]
                    + f_0 * mi_1407[k];

        t_1156[k] = -3.0 * ki_932[k]
                    + f_0 * mi_1408[k];

        t_1157[k] = -3.0 * ki_933[k]
                    + f_0 * mi_1409[k];

        t_1158[k] = -3.0 * ki_934[k]
                    + f_0 * mi_1410[k];

        t_1159[k] = -3.0 * ki_935[k]
                    + f_0 * mi_1411[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, t_1163, t_1164, ki_936, ki_937, ki_938, \
                         ki_939, ki_940, mi_1412, mi_1413, mi_1414, mi_1415, \
                         mi_1416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = -3.0 * ki_936[k]
                    + f_0 * mi_1412[k];

        t_1161[k] = -3.0 * ki_937[k]
                    + f_0 * mi_1413[k];

        t_1162[k] = -3.0 * ki_938[k]
                    + f_0 * mi_1414[k];

        t_1163[k] = -3.0 * ki_939[k]
                    + f_0 * mi_1415[k];

        t_1164[k] = -3.0 * ki_940[k]
                    + f_0 * mi_1416[k];
    }

#pragma omp simd aligned(t_1165, t_1166, t_1167, t_1168, t_1169, ki_941, ki_942, ki_943, \
                         ki_944, ki_945, mi_1417, mi_1418, mi_1419, mi_1420, \
                         mi_1421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1165[k] = -3.0 * ki_941[k]
                    + f_0 * mi_1417[k];

        t_1166[k] = -3.0 * ki_942[k]
                    + f_0 * mi_1418[k];

        t_1167[k] = -3.0 * ki_943[k]
                    + f_0 * mi_1419[k];

        t_1168[k] = -3.0 * ki_944[k]
                    + f_0 * mi_1420[k];

        t_1169[k] = -3.0 * ki_945[k]
                    + f_0 * mi_1421[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, t_1173, t_1174, ki_946, ki_947, ki_948, \
                         ki_949, ki_950, mi_1422, mi_1423, mi_1424, mi_1425, \
                         mi_1426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = -3.0 * ki_946[k]
                    + f_0 * mi_1422[k];

        t_1171[k] = -3.0 * ki_947[k]
                    + f_0 * mi_1423[k];

        t_1172[k] = -3.0 * ki_948[k]
                    + f_0 * mi_1424[k];

        t_1173[k] = -3.0 * ki_949[k]
                    + f_0 * mi_1425[k];

        t_1174[k] = -3.0 * ki_950[k]
                    + f_0 * mi_1426[k];
    }

#pragma omp simd aligned(t_1175, t_1176, t_1177, t_1178, t_1179, ki_951, ki_952, ki_953, \
                         ki_954, ki_955, mi_1427, mi_1428, mi_1429, mi_1430, \
                         mi_1431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1175[k] = -3.0 * ki_951[k]
                    + f_0 * mi_1427[k];

        t_1176[k] = -2.0 * ki_952[k]
                    + f_0 * mi_1428[k];

        t_1177[k] = -2.0 * ki_953[k]
                    + f_0 * mi_1429[k];

        t_1178[k] = -2.0 * ki_954[k]
                    + f_0 * mi_1430[k];

        t_1179[k] = -2.0 * ki_955[k]
                    + f_0 * mi_1431[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, t_1183, t_1184, ki_956, ki_957, ki_958, \
                         ki_959, ki_960, mi_1432, mi_1433, mi_1434, mi_1435, \
                         mi_1436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = -2.0 * ki_956[k]
                    + f_0 * mi_1432[k];

        t_1181[k] = -2.0 * ki_957[k]
                    + f_0 * mi_1433[k];

        t_1182[k] = -2.0 * ki_958[k]
                    + f_0 * mi_1434[k];

        t_1183[k] = -2.0 * ki_959[k]
                    + f_0 * mi_1435[k];

        t_1184[k] = -2.0 * ki_960[k]
                    + f_0 * mi_1436[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, t_1188, t_1189, ki_961, ki_962, ki_963, \
                         ki_964, ki_965, mi_1437, mi_1438, mi_1439, mi_1440, \
                         mi_1441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = -2.0 * ki_961[k]
                    + f_0 * mi_1437[k];

        t_1186[k] = -2.0 * ki_962[k]
                    + f_0 * mi_1438[k];

        t_1187[k] = -2.0 * ki_963[k]
                    + f_0 * mi_1439[k];

        t_1188[k] = -2.0 * ki_964[k]
                    + f_0 * mi_1440[k];

        t_1189[k] = -2.0 * ki_965[k]
                    + f_0 * mi_1441[k];
    }

#pragma omp simd aligned(t_1190, t_1191, t_1192, t_1193, t_1194, ki_966, ki_967, ki_968, \
                         ki_969, ki_970, mi_1442, mi_1443, mi_1444, mi_1445, \
                         mi_1446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1190[k] = -2.0 * ki_966[k]
                    + f_0 * mi_1442[k];

        t_1191[k] = -2.0 * ki_967[k]
                    + f_0 * mi_1443[k];

        t_1192[k] = -2.0 * ki_968[k]
                    + f_0 * mi_1444[k];

        t_1193[k] = -2.0 * ki_969[k]
                    + f_0 * mi_1445[k];

        t_1194[k] = -2.0 * ki_970[k]
                    + f_0 * mi_1446[k];
    }

#pragma omp simd aligned(t_1195, t_1196, t_1197, t_1198, t_1199, ki_971, ki_972, ki_973, \
                         ki_974, ki_975, mi_1447, mi_1448, mi_1449, mi_1450, \
                         mi_1451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1195[k] = -2.0 * ki_971[k]
                    + f_0 * mi_1447[k];

        t_1196[k] = -2.0 * ki_972[k]
                    + f_0 * mi_1448[k];

        t_1197[k] = -2.0 * ki_973[k]
                    + f_0 * mi_1449[k];

        t_1198[k] = -2.0 * ki_974[k]
                    + f_0 * mi_1450[k];

        t_1199[k] = -2.0 * ki_975[k]
                    + f_0 * mi_1451[k];
    }

#pragma omp simd aligned(t_1200, t_1201, t_1202, t_1203, t_1204, ki_976, ki_977, ki_978, \
                         ki_979, ki_980, mi_1452, mi_1453, mi_1454, mi_1455, \
                         mi_1456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1200[k] = -2.0 * ki_976[k]
                    + f_0 * mi_1452[k];

        t_1201[k] = -2.0 * ki_977[k]
                    + f_0 * mi_1453[k];

        t_1202[k] = -2.0 * ki_978[k]
                    + f_0 * mi_1454[k];

        t_1203[k] = -2.0 * ki_979[k]
                    + f_0 * mi_1455[k];

        t_1204[k] = -ki_980[k]
                    + f_0 * mi_1456[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, t_1208, t_1209, ki_981, ki_982, ki_983, \
                         ki_984, ki_985, mi_1457, mi_1458, mi_1459, mi_1460, \
                         mi_1461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = -ki_981[k]
                    + f_0 * mi_1457[k];

        t_1206[k] = -ki_982[k]
                    + f_0 * mi_1458[k];

        t_1207[k] = -ki_983[k]
                    + f_0 * mi_1459[k];

        t_1208[k] = -ki_984[k]
                    + f_0 * mi_1460[k];

        t_1209[k] = -ki_985[k]
                    + f_0 * mi_1461[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, t_1214, ki_986, ki_987, ki_988, \
                         ki_989, ki_990, mi_1462, mi_1463, mi_1464, mi_1465, \
                         mi_1466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = -ki_986[k]
                    + f_0 * mi_1462[k];

        t_1211[k] = -ki_987[k]
                    + f_0 * mi_1463[k];

        t_1212[k] = -ki_988[k]
                    + f_0 * mi_1464[k];

        t_1213[k] = -ki_989[k]
                    + f_0 * mi_1465[k];

        t_1214[k] = -ki_990[k]
                    + f_0 * mi_1466[k];
    }

#pragma omp simd aligned(t_1215, t_1216, t_1217, t_1218, t_1219, ki_991, ki_992, ki_993, \
                         ki_994, ki_995, mi_1467, mi_1468, mi_1469, mi_1470, \
                         mi_1471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1215[k] = -ki_991[k]
                    + f_0 * mi_1467[k];

        t_1216[k] = -ki_992[k]
                    + f_0 * mi_1468[k];

        t_1217[k] = -ki_993[k]
                    + f_0 * mi_1469[k];

        t_1218[k] = -ki_994[k]
                    + f_0 * mi_1470[k];

        t_1219[k] = -ki_995[k]
                    + f_0 * mi_1471[k];
    }

#pragma omp simd aligned(t_1220, t_1221, t_1222, t_1223, t_1224, ki_996, ki_997, ki_998, \
                         ki_999, ki_1000, mi_1472, mi_1473, mi_1474, mi_1475, \
                         mi_1476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1220[k] = -ki_996[k]
                    + f_0 * mi_1472[k];

        t_1221[k] = -ki_997[k]
                    + f_0 * mi_1473[k];

        t_1222[k] = -ki_998[k]
                    + f_0 * mi_1474[k];

        t_1223[k] = -ki_999[k]
                    + f_0 * mi_1475[k];

        t_1224[k] = -ki_1000[k]
                    + f_0 * mi_1476[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, t_1228, t_1229, ki_1001, ki_1002, ki_1003, \
                         ki_1004, ki_1005, mi_1477, mi_1478, mi_1479, mi_1480, \
                         mi_1481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = -ki_1001[k]
                    + f_0 * mi_1477[k];

        t_1226[k] = -ki_1002[k]
                    + f_0 * mi_1478[k];

        t_1227[k] = -ki_1003[k]
                    + f_0 * mi_1479[k];

        t_1228[k] = -ki_1004[k]
                    + f_0 * mi_1480[k];

        t_1229[k] = -ki_1005[k]
                    + f_0 * mi_1481[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, t_1233, t_1234, t_1235, t_1236, ki_1006, \
                         ki_1007, mi_1482, mi_1483, mi_1484, mi_1485, mi_1486, mi_1487, \
                         mi_1488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = -ki_1006[k]
                    + f_0 * mi_1482[k];

        t_1231[k] = -ki_1007[k]
                    + f_0 * mi_1483[k];

        t_1232[k] = f_0 * mi_1484[k];

        t_1233[k] = f_0 * mi_1485[k];

        t_1234[k] = f_0 * mi_1486[k];

        t_1235[k] = f_0 * mi_1487[k];

        t_1236[k] = f_0 * mi_1488[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, t_1241, t_1242, t_1243, t_1244, \
                         mi_1489, mi_1490, mi_1491, mi_1492, mi_1493, mi_1494, mi_1495, \
                         mi_1496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_0 * mi_1489[k];

        t_1238[k] = f_0 * mi_1490[k];

        t_1239[k] = f_0 * mi_1491[k];

        t_1240[k] = f_0 * mi_1492[k];

        t_1241[k] = f_0 * mi_1493[k];

        t_1242[k] = f_0 * mi_1494[k];

        t_1243[k] = f_0 * mi_1495[k];

        t_1244[k] = f_0 * mi_1496[k];
    }

#pragma omp simd aligned(t_1245, t_1246, t_1247, t_1248, t_1249, t_1250, t_1251, t_1252, \
                         mi_1497, mi_1498, mi_1499, mi_1500, mi_1501, mi_1502, mi_1503, \
                         mi_1504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1245[k] = f_0 * mi_1497[k];

        t_1246[k] = f_0 * mi_1498[k];

        t_1247[k] = f_0 * mi_1499[k];

        t_1248[k] = f_0 * mi_1500[k];

        t_1249[k] = f_0 * mi_1501[k];

        t_1250[k] = f_0 * mi_1502[k];

        t_1251[k] = f_0 * mi_1503[k];

        t_1252[k] = f_0 * mi_1504[k];
    }

#pragma omp simd aligned(t_1253, t_1254, t_1255, t_1256, t_1257, t_1258, t_1259, mi_1505, \
                         mi_1506, mi_1507, mi_1508, mi_1509, mi_1510, \
                         mi_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1253[k] = f_0 * mi_1505[k];

        t_1254[k] = f_0 * mi_1506[k];

        t_1255[k] = f_0 * mi_1507[k];

        t_1256[k] = f_0 * mi_1508[k];

        t_1257[k] = f_0 * mi_1509[k];

        t_1258[k] = f_0 * mi_1510[k];

        t_1259[k] = f_0 * mi_1511[k];
    }
}

auto
compute_prim_geom_10_li_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t ki, const size_t mi,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_li_electron_repulsion_1_piece0(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_1_piece1(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_1_piece2(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_1_piece3(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_1_piece4(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_1_piece5(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_1_piece6(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_1_piece7(buffer, target, ki, mi, ncols, alpha);
}

static auto
compute_prim_geom_10_li_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
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
    auto *t_177 = buffer.data(target + 177);
    auto *t_178 = buffer.data(target + 178);
    auto *t_179 = buffer.data(target + 179);
    auto *t_180 = buffer.data(target + 180);

    const auto *ki_0 = buffer.data(ki + 0);
    const auto *ki_1 = buffer.data(ki + 1);
    const auto *ki_2 = buffer.data(ki + 2);
    const auto *ki_3 = buffer.data(ki + 3);
    const auto *ki_4 = buffer.data(ki + 4);
    const auto *ki_5 = buffer.data(ki + 5);
    const auto *ki_6 = buffer.data(ki + 6);
    const auto *ki_7 = buffer.data(ki + 7);
    const auto *ki_8 = buffer.data(ki + 8);
    const auto *ki_9 = buffer.data(ki + 9);
    const auto *ki_10 = buffer.data(ki + 10);
    const auto *ki_11 = buffer.data(ki + 11);
    const auto *ki_12 = buffer.data(ki + 12);
    const auto *ki_13 = buffer.data(ki + 13);
    const auto *ki_14 = buffer.data(ki + 14);
    const auto *ki_15 = buffer.data(ki + 15);
    const auto *ki_16 = buffer.data(ki + 16);
    const auto *ki_17 = buffer.data(ki + 17);
    const auto *ki_18 = buffer.data(ki + 18);
    const auto *ki_19 = buffer.data(ki + 19);
    const auto *ki_20 = buffer.data(ki + 20);
    const auto *ki_21 = buffer.data(ki + 21);
    const auto *ki_22 = buffer.data(ki + 22);
    const auto *ki_23 = buffer.data(ki + 23);
    const auto *ki_24 = buffer.data(ki + 24);
    const auto *ki_25 = buffer.data(ki + 25);
    const auto *ki_26 = buffer.data(ki + 26);
    const auto *ki_27 = buffer.data(ki + 27);
    const auto *ki_28 = buffer.data(ki + 28);
    const auto *ki_29 = buffer.data(ki + 29);
    const auto *ki_30 = buffer.data(ki + 30);
    const auto *ki_31 = buffer.data(ki + 31);
    const auto *ki_32 = buffer.data(ki + 32);
    const auto *ki_33 = buffer.data(ki + 33);
    const auto *ki_34 = buffer.data(ki + 34);
    const auto *ki_35 = buffer.data(ki + 35);
    const auto *ki_36 = buffer.data(ki + 36);
    const auto *ki_37 = buffer.data(ki + 37);
    const auto *ki_38 = buffer.data(ki + 38);
    const auto *ki_39 = buffer.data(ki + 39);
    const auto *ki_40 = buffer.data(ki + 40);
    const auto *ki_41 = buffer.data(ki + 41);
    const auto *ki_42 = buffer.data(ki + 42);
    const auto *ki_43 = buffer.data(ki + 43);
    const auto *ki_44 = buffer.data(ki + 44);
    const auto *ki_45 = buffer.data(ki + 45);
    const auto *ki_46 = buffer.data(ki + 46);
    const auto *ki_47 = buffer.data(ki + 47);
    const auto *ki_48 = buffer.data(ki + 48);
    const auto *ki_49 = buffer.data(ki + 49);
    const auto *ki_50 = buffer.data(ki + 50);
    const auto *ki_51 = buffer.data(ki + 51);
    const auto *ki_52 = buffer.data(ki + 52);
    const auto *ki_53 = buffer.data(ki + 53);
    const auto *ki_54 = buffer.data(ki + 54);
    const auto *ki_55 = buffer.data(ki + 55);
    const auto *ki_56 = buffer.data(ki + 56);
    const auto *ki_57 = buffer.data(ki + 57);
    const auto *ki_58 = buffer.data(ki + 58);
    const auto *ki_59 = buffer.data(ki + 59);
    const auto *ki_60 = buffer.data(ki + 60);
    const auto *ki_61 = buffer.data(ki + 61);
    const auto *ki_62 = buffer.data(ki + 62);
    const auto *ki_63 = buffer.data(ki + 63);
    const auto *ki_64 = buffer.data(ki + 64);
    const auto *ki_65 = buffer.data(ki + 65);
    const auto *ki_66 = buffer.data(ki + 66);
    const auto *ki_67 = buffer.data(ki + 67);
    const auto *ki_68 = buffer.data(ki + 68);
    const auto *ki_69 = buffer.data(ki + 69);
    const auto *ki_70 = buffer.data(ki + 70);
    const auto *ki_71 = buffer.data(ki + 71);
    const auto *ki_72 = buffer.data(ki + 72);
    const auto *ki_73 = buffer.data(ki + 73);
    const auto *ki_74 = buffer.data(ki + 74);
    const auto *ki_75 = buffer.data(ki + 75);
    const auto *ki_76 = buffer.data(ki + 76);
    const auto *ki_77 = buffer.data(ki + 77);
    const auto *ki_78 = buffer.data(ki + 78);
    const auto *ki_79 = buffer.data(ki + 79);
    const auto *ki_80 = buffer.data(ki + 80);
    const auto *ki_81 = buffer.data(ki + 81);
    const auto *ki_82 = buffer.data(ki + 82);
    const auto *ki_83 = buffer.data(ki + 83);

    const auto *mi_56 = buffer.data(mi + 56);
    const auto *mi_57 = buffer.data(mi + 57);
    const auto *mi_58 = buffer.data(mi + 58);
    const auto *mi_59 = buffer.data(mi + 59);
    const auto *mi_60 = buffer.data(mi + 60);
    const auto *mi_61 = buffer.data(mi + 61);
    const auto *mi_62 = buffer.data(mi + 62);
    const auto *mi_63 = buffer.data(mi + 63);
    const auto *mi_64 = buffer.data(mi + 64);
    const auto *mi_65 = buffer.data(mi + 65);
    const auto *mi_66 = buffer.data(mi + 66);
    const auto *mi_67 = buffer.data(mi + 67);
    const auto *mi_68 = buffer.data(mi + 68);
    const auto *mi_69 = buffer.data(mi + 69);
    const auto *mi_70 = buffer.data(mi + 70);
    const auto *mi_71 = buffer.data(mi + 71);
    const auto *mi_72 = buffer.data(mi + 72);
    const auto *mi_73 = buffer.data(mi + 73);
    const auto *mi_74 = buffer.data(mi + 74);
    const auto *mi_75 = buffer.data(mi + 75);
    const auto *mi_76 = buffer.data(mi + 76);
    const auto *mi_77 = buffer.data(mi + 77);
    const auto *mi_78 = buffer.data(mi + 78);
    const auto *mi_79 = buffer.data(mi + 79);
    const auto *mi_80 = buffer.data(mi + 80);
    const auto *mi_81 = buffer.data(mi + 81);
    const auto *mi_82 = buffer.data(mi + 82);
    const auto *mi_83 = buffer.data(mi + 83);
    const auto *mi_112 = buffer.data(mi + 112);
    const auto *mi_113 = buffer.data(mi + 113);
    const auto *mi_114 = buffer.data(mi + 114);
    const auto *mi_115 = buffer.data(mi + 115);
    const auto *mi_116 = buffer.data(mi + 116);
    const auto *mi_117 = buffer.data(mi + 117);
    const auto *mi_118 = buffer.data(mi + 118);
    const auto *mi_119 = buffer.data(mi + 119);
    const auto *mi_120 = buffer.data(mi + 120);
    const auto *mi_121 = buffer.data(mi + 121);
    const auto *mi_122 = buffer.data(mi + 122);
    const auto *mi_123 = buffer.data(mi + 123);
    const auto *mi_124 = buffer.data(mi + 124);
    const auto *mi_125 = buffer.data(mi + 125);
    const auto *mi_126 = buffer.data(mi + 126);
    const auto *mi_127 = buffer.data(mi + 127);
    const auto *mi_128 = buffer.data(mi + 128);
    const auto *mi_129 = buffer.data(mi + 129);
    const auto *mi_130 = buffer.data(mi + 130);
    const auto *mi_131 = buffer.data(mi + 131);
    const auto *mi_132 = buffer.data(mi + 132);
    const auto *mi_133 = buffer.data(mi + 133);
    const auto *mi_134 = buffer.data(mi + 134);
    const auto *mi_135 = buffer.data(mi + 135);
    const auto *mi_136 = buffer.data(mi + 136);
    const auto *mi_137 = buffer.data(mi + 137);
    const auto *mi_138 = buffer.data(mi + 138);
    const auto *mi_139 = buffer.data(mi + 139);
    const auto *mi_140 = buffer.data(mi + 140);
    const auto *mi_141 = buffer.data(mi + 141);
    const auto *mi_142 = buffer.data(mi + 142);
    const auto *mi_143 = buffer.data(mi + 143);
    const auto *mi_144 = buffer.data(mi + 144);
    const auto *mi_145 = buffer.data(mi + 145);
    const auto *mi_146 = buffer.data(mi + 146);
    const auto *mi_147 = buffer.data(mi + 147);
    const auto *mi_148 = buffer.data(mi + 148);
    const auto *mi_149 = buffer.data(mi + 149);
    const auto *mi_150 = buffer.data(mi + 150);
    const auto *mi_151 = buffer.data(mi + 151);
    const auto *mi_152 = buffer.data(mi + 152);
    const auto *mi_153 = buffer.data(mi + 153);
    const auto *mi_154 = buffer.data(mi + 154);
    const auto *mi_155 = buffer.data(mi + 155);
    const auto *mi_156 = buffer.data(mi + 156);
    const auto *mi_157 = buffer.data(mi + 157);
    const auto *mi_158 = buffer.data(mi + 158);
    const auto *mi_159 = buffer.data(mi + 159);
    const auto *mi_160 = buffer.data(mi + 160);
    const auto *mi_161 = buffer.data(mi + 161);
    const auto *mi_162 = buffer.data(mi + 162);
    const auto *mi_163 = buffer.data(mi + 163);
    const auto *mi_164 = buffer.data(mi + 164);
    const auto *mi_165 = buffer.data(mi + 165);
    const auto *mi_166 = buffer.data(mi + 166);
    const auto *mi_167 = buffer.data(mi + 167);
    const auto *mi_196 = buffer.data(mi + 196);
    const auto *mi_197 = buffer.data(mi + 197);
    const auto *mi_198 = buffer.data(mi + 198);
    const auto *mi_199 = buffer.data(mi + 199);
    const auto *mi_200 = buffer.data(mi + 200);
    const auto *mi_201 = buffer.data(mi + 201);
    const auto *mi_202 = buffer.data(mi + 202);
    const auto *mi_203 = buffer.data(mi + 203);
    const auto *mi_204 = buffer.data(mi + 204);
    const auto *mi_205 = buffer.data(mi + 205);
    const auto *mi_206 = buffer.data(mi + 206);
    const auto *mi_207 = buffer.data(mi + 207);
    const auto *mi_208 = buffer.data(mi + 208);
    const auto *mi_209 = buffer.data(mi + 209);
    const auto *mi_210 = buffer.data(mi + 210);
    const auto *mi_211 = buffer.data(mi + 211);
    const auto *mi_212 = buffer.data(mi + 212);
    const auto *mi_213 = buffer.data(mi + 213);
    const auto *mi_214 = buffer.data(mi + 214);
    const auto *mi_215 = buffer.data(mi + 215);
    const auto *mi_216 = buffer.data(mi + 216);
    const auto *mi_217 = buffer.data(mi + 217);
    const auto *mi_218 = buffer.data(mi + 218);
    const auto *mi_219 = buffer.data(mi + 219);
    const auto *mi_220 = buffer.data(mi + 220);
    const auto *mi_221 = buffer.data(mi + 221);
    const auto *mi_222 = buffer.data(mi + 222);
    const auto *mi_223 = buffer.data(mi + 223);
    const auto *mi_224 = buffer.data(mi + 224);
    const auto *mi_225 = buffer.data(mi + 225);
    const auto *mi_226 = buffer.data(mi + 226);
    const auto *mi_227 = buffer.data(mi + 227);
    const auto *mi_228 = buffer.data(mi + 228);
    const auto *mi_229 = buffer.data(mi + 229);
    const auto *mi_230 = buffer.data(mi + 230);
    const auto *mi_231 = buffer.data(mi + 231);
    const auto *mi_232 = buffer.data(mi + 232);
    const auto *mi_233 = buffer.data(mi + 233);
    const auto *mi_234 = buffer.data(mi + 234);
    const auto *mi_235 = buffer.data(mi + 235);
    const auto *mi_236 = buffer.data(mi + 236);
    const auto *mi_237 = buffer.data(mi + 237);
    const auto *mi_238 = buffer.data(mi + 238);
    const auto *mi_239 = buffer.data(mi + 239);
    const auto *mi_240 = buffer.data(mi + 240);
    const auto *mi_241 = buffer.data(mi + 241);
    const auto *mi_242 = buffer.data(mi + 242);
    const auto *mi_243 = buffer.data(mi + 243);
    const auto *mi_244 = buffer.data(mi + 244);
    const auto *mi_245 = buffer.data(mi + 245);
    const auto *mi_246 = buffer.data(mi + 246);
    const auto *mi_247 = buffer.data(mi + 247);
    const auto *mi_248 = buffer.data(mi + 248);
    const auto *mi_249 = buffer.data(mi + 249);
    const auto *mi_250 = buffer.data(mi + 250);
    const auto *mi_251 = buffer.data(mi + 251);
    const auto *mi_252 = buffer.data(mi + 252);
    const auto *mi_253 = buffer.data(mi + 253);
    const auto *mi_254 = buffer.data(mi + 254);
    const auto *mi_255 = buffer.data(mi + 255);
    const auto *mi_256 = buffer.data(mi + 256);
    const auto *mi_257 = buffer.data(mi + 257);
    const auto *mi_258 = buffer.data(mi + 258);
    const auto *mi_259 = buffer.data(mi + 259);
    const auto *mi_260 = buffer.data(mi + 260);
    const auto *mi_261 = buffer.data(mi + 261);
    const auto *mi_262 = buffer.data(mi + 262);
    const auto *mi_263 = buffer.data(mi + 263);
    const auto *mi_264 = buffer.data(mi + 264);
    const auto *mi_265 = buffer.data(mi + 265);
    const auto *mi_266 = buffer.data(mi + 266);
    const auto *mi_267 = buffer.data(mi + 267);
    const auto *mi_268 = buffer.data(mi + 268);
    const auto *mi_269 = buffer.data(mi + 269);
    const auto *mi_270 = buffer.data(mi + 270);
    const auto *mi_271 = buffer.data(mi + 271);
    const auto *mi_272 = buffer.data(mi + 272);
    const auto *mi_273 = buffer.data(mi + 273);
    const auto *mi_274 = buffer.data(mi + 274);
    const auto *mi_275 = buffer.data(mi + 275);
    const auto *mi_276 = buffer.data(mi + 276);
    const auto *mi_277 = buffer.data(mi + 277);
    const auto *mi_278 = buffer.data(mi + 278);
    const auto *mi_279 = buffer.data(mi + 279);
    const auto *mi_308 = buffer.data(mi + 308);
    const auto *mi_309 = buffer.data(mi + 309);
    const auto *mi_310 = buffer.data(mi + 310);
    const auto *mi_311 = buffer.data(mi + 311);
    const auto *mi_312 = buffer.data(mi + 312);
    const auto *mi_313 = buffer.data(mi + 313);
    const auto *mi_314 = buffer.data(mi + 314);
    const auto *mi_315 = buffer.data(mi + 315);
    const auto *mi_316 = buffer.data(mi + 316);
    const auto *mi_317 = buffer.data(mi + 317);
    const auto *mi_318 = buffer.data(mi + 318);
    const auto *mi_319 = buffer.data(mi + 319);
    const auto *mi_320 = buffer.data(mi + 320);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, mi_56, mi_57, mi_58, mi_59, \
                         mi_60, mi_61, mi_62, mi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * mi_56[k];

        t_1[k] = f_0 * mi_57[k];

        t_2[k] = f_0 * mi_58[k];

        t_3[k] = f_0 * mi_59[k];

        t_4[k] = f_0 * mi_60[k];

        t_5[k] = f_0 * mi_61[k];

        t_6[k] = f_0 * mi_62[k];

        t_7[k] = f_0 * mi_63[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, mi_64, mi_65, mi_66, \
                         mi_67, mi_68, mi_69, mi_70, mi_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * mi_64[k];

        t_9[k] = f_0 * mi_65[k];

        t_10[k] = f_0 * mi_66[k];

        t_11[k] = f_0 * mi_67[k];

        t_12[k] = f_0 * mi_68[k];

        t_13[k] = f_0 * mi_69[k];

        t_14[k] = f_0 * mi_70[k];

        t_15[k] = f_0 * mi_71[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, mi_72, mi_73, mi_74, \
                         mi_75, mi_76, mi_77, mi_78, mi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * mi_72[k];

        t_17[k] = f_0 * mi_73[k];

        t_18[k] = f_0 * mi_74[k];

        t_19[k] = f_0 * mi_75[k];

        t_20[k] = f_0 * mi_76[k];

        t_21[k] = f_0 * mi_77[k];

        t_22[k] = f_0 * mi_78[k];

        t_23[k] = f_0 * mi_79[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, mi_80, mi_81, mi_82, \
                         mi_83, mi_112, mi_113, mi_114, mi_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * mi_80[k];

        t_25[k] = f_0 * mi_81[k];

        t_26[k] = f_0 * mi_82[k];

        t_27[k] = f_0 * mi_83[k];

        t_28[k] = f_0 * mi_112[k];

        t_29[k] = f_0 * mi_113[k];

        t_30[k] = f_0 * mi_114[k];

        t_31[k] = f_0 * mi_115[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, mi_116, mi_117, \
                         mi_118, mi_119, mi_120, mi_121, mi_122, \
                         mi_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * mi_116[k];

        t_33[k] = f_0 * mi_117[k];

        t_34[k] = f_0 * mi_118[k];

        t_35[k] = f_0 * mi_119[k];

        t_36[k] = f_0 * mi_120[k];

        t_37[k] = f_0 * mi_121[k];

        t_38[k] = f_0 * mi_122[k];

        t_39[k] = f_0 * mi_123[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, mi_124, mi_125, \
                         mi_126, mi_127, mi_128, mi_129, mi_130, \
                         mi_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * mi_124[k];

        t_41[k] = f_0 * mi_125[k];

        t_42[k] = f_0 * mi_126[k];

        t_43[k] = f_0 * mi_127[k];

        t_44[k] = f_0 * mi_128[k];

        t_45[k] = f_0 * mi_129[k];

        t_46[k] = f_0 * mi_130[k];

        t_47[k] = f_0 * mi_131[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, mi_132, mi_133, \
                         mi_134, mi_135, mi_136, mi_137, mi_138, \
                         mi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * mi_132[k];

        t_49[k] = f_0 * mi_133[k];

        t_50[k] = f_0 * mi_134[k];

        t_51[k] = f_0 * mi_135[k];

        t_52[k] = f_0 * mi_136[k];

        t_53[k] = f_0 * mi_137[k];

        t_54[k] = f_0 * mi_138[k];

        t_55[k] = f_0 * mi_139[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, ki_0, ki_1, ki_2, ki_3, ki_4, mi_140, \
                         mi_141, mi_142, mi_143, mi_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -ki_0[k]
                  + f_0 * mi_140[k];

        t_57[k] = -ki_1[k]
                  + f_0 * mi_141[k];

        t_58[k] = -ki_2[k]
                  + f_0 * mi_142[k];

        t_59[k] = -ki_3[k]
                  + f_0 * mi_143[k];

        t_60[k] = -ki_4[k]
                  + f_0 * mi_144[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, ki_5, ki_6, ki_7, ki_8, ki_9, mi_145, \
                         mi_146, mi_147, mi_148, mi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -ki_5[k]
                  + f_0 * mi_145[k];

        t_62[k] = -ki_6[k]
                  + f_0 * mi_146[k];

        t_63[k] = -ki_7[k]
                  + f_0 * mi_147[k];

        t_64[k] = -ki_8[k]
                  + f_0 * mi_148[k];

        t_65[k] = -ki_9[k]
                  + f_0 * mi_149[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, ki_10, ki_11, ki_12, ki_13, ki_14, \
                         mi_150, mi_151, mi_152, mi_153, mi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -ki_10[k]
                  + f_0 * mi_150[k];

        t_67[k] = -ki_11[k]
                  + f_0 * mi_151[k];

        t_68[k] = -ki_12[k]
                  + f_0 * mi_152[k];

        t_69[k] = -ki_13[k]
                  + f_0 * mi_153[k];

        t_70[k] = -ki_14[k]
                  + f_0 * mi_154[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, ki_15, ki_16, ki_17, ki_18, ki_19, \
                         mi_155, mi_156, mi_157, mi_158, mi_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -ki_15[k]
                  + f_0 * mi_155[k];

        t_72[k] = -ki_16[k]
                  + f_0 * mi_156[k];

        t_73[k] = -ki_17[k]
                  + f_0 * mi_157[k];

        t_74[k] = -ki_18[k]
                  + f_0 * mi_158[k];

        t_75[k] = -ki_19[k]
                  + f_0 * mi_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, ki_20, ki_21, ki_22, ki_23, ki_24, \
                         mi_160, mi_161, mi_162, mi_163, mi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -ki_20[k]
                  + f_0 * mi_160[k];

        t_77[k] = -ki_21[k]
                  + f_0 * mi_161[k];

        t_78[k] = -ki_22[k]
                  + f_0 * mi_162[k];

        t_79[k] = -ki_23[k]
                  + f_0 * mi_163[k];

        t_80[k] = -ki_24[k]
                  + f_0 * mi_164[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, ki_25, ki_26, ki_27, mi_165, \
                         mi_166, mi_167, mi_196, mi_197, mi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -ki_25[k]
                  + f_0 * mi_165[k];

        t_82[k] = -ki_26[k]
                  + f_0 * mi_166[k];

        t_83[k] = -ki_27[k]
                  + f_0 * mi_167[k];

        t_84[k] = f_0 * mi_196[k];

        t_85[k] = f_0 * mi_197[k];

        t_86[k] = f_0 * mi_198[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, t_93, t_94, mi_199, mi_200, \
                         mi_201, mi_202, mi_203, mi_204, mi_205, \
                         mi_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_0 * mi_199[k];

        t_88[k] = f_0 * mi_200[k];

        t_89[k] = f_0 * mi_201[k];

        t_90[k] = f_0 * mi_202[k];

        t_91[k] = f_0 * mi_203[k];

        t_92[k] = f_0 * mi_204[k];

        t_93[k] = f_0 * mi_205[k];

        t_94[k] = f_0 * mi_206[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, t_100, t_101, t_102, mi_207, mi_208, \
                         mi_209, mi_210, mi_211, mi_212, mi_213, \
                         mi_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_0 * mi_207[k];

        t_96[k] = f_0 * mi_208[k];

        t_97[k] = f_0 * mi_209[k];

        t_98[k] = f_0 * mi_210[k];

        t_99[k] = f_0 * mi_211[k];

        t_100[k] = f_0 * mi_212[k];

        t_101[k] = f_0 * mi_213[k];

        t_102[k] = f_0 * mi_214[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, t_109, t_110, mi_215, \
                         mi_216, mi_217, mi_218, mi_219, mi_220, mi_221, \
                         mi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_0 * mi_215[k];

        t_104[k] = f_0 * mi_216[k];

        t_105[k] = f_0 * mi_217[k];

        t_106[k] = f_0 * mi_218[k];

        t_107[k] = f_0 * mi_219[k];

        t_108[k] = f_0 * mi_220[k];

        t_109[k] = f_0 * mi_221[k];

        t_110[k] = f_0 * mi_222[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, ki_28, ki_29, ki_30, ki_31, \
                         mi_223, mi_224, mi_225, mi_226, mi_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_0 * mi_223[k];

        t_112[k] = -ki_28[k]
                   + f_0 * mi_224[k];

        t_113[k] = -ki_29[k]
                   + f_0 * mi_225[k];

        t_114[k] = -ki_30[k]
                   + f_0 * mi_226[k];

        t_115[k] = -ki_31[k]
                   + f_0 * mi_227[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, ki_32, ki_33, ki_34, ki_35, ki_36, \
                         mi_228, mi_229, mi_230, mi_231, mi_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = -ki_32[k]
                   + f_0 * mi_228[k];

        t_117[k] = -ki_33[k]
                   + f_0 * mi_229[k];

        t_118[k] = -ki_34[k]
                   + f_0 * mi_230[k];

        t_119[k] = -ki_35[k]
                   + f_0 * mi_231[k];

        t_120[k] = -ki_36[k]
                   + f_0 * mi_232[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, ki_37, ki_38, ki_39, ki_40, ki_41, \
                         mi_233, mi_234, mi_235, mi_236, mi_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -ki_37[k]
                   + f_0 * mi_233[k];

        t_122[k] = -ki_38[k]
                   + f_0 * mi_234[k];

        t_123[k] = -ki_39[k]
                   + f_0 * mi_235[k];

        t_124[k] = -ki_40[k]
                   + f_0 * mi_236[k];

        t_125[k] = -ki_41[k]
                   + f_0 * mi_237[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, ki_42, ki_43, ki_44, ki_45, ki_46, \
                         mi_238, mi_239, mi_240, mi_241, mi_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = -ki_42[k]
                   + f_0 * mi_238[k];

        t_127[k] = -ki_43[k]
                   + f_0 * mi_239[k];

        t_128[k] = -ki_44[k]
                   + f_0 * mi_240[k];

        t_129[k] = -ki_45[k]
                   + f_0 * mi_241[k];

        t_130[k] = -ki_46[k]
                   + f_0 * mi_242[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, ki_47, ki_48, ki_49, ki_50, ki_51, \
                         mi_243, mi_244, mi_245, mi_246, mi_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = -ki_47[k]
                   + f_0 * mi_243[k];

        t_132[k] = -ki_48[k]
                   + f_0 * mi_244[k];

        t_133[k] = -ki_49[k]
                   + f_0 * mi_245[k];

        t_134[k] = -ki_50[k]
                   + f_0 * mi_246[k];

        t_135[k] = -ki_51[k]
                   + f_0 * mi_247[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, ki_52, ki_53, ki_54, ki_55, ki_56, \
                         mi_248, mi_249, mi_250, mi_251, mi_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = -ki_52[k]
                   + f_0 * mi_248[k];

        t_137[k] = -ki_53[k]
                   + f_0 * mi_249[k];

        t_138[k] = -ki_54[k]
                   + f_0 * mi_250[k];

        t_139[k] = -ki_55[k]
                   + f_0 * mi_251[k];

        t_140[k] = -2.0 * ki_56[k]
                   + f_0 * mi_252[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, ki_57, ki_58, ki_59, ki_60, ki_61, \
                         mi_253, mi_254, mi_255, mi_256, mi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = -2.0 * ki_57[k]
                   + f_0 * mi_253[k];

        t_142[k] = -2.0 * ki_58[k]
                   + f_0 * mi_254[k];

        t_143[k] = -2.0 * ki_59[k]
                   + f_0 * mi_255[k];

        t_144[k] = -2.0 * ki_60[k]
                   + f_0 * mi_256[k];

        t_145[k] = -2.0 * ki_61[k]
                   + f_0 * mi_257[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, ki_62, ki_63, ki_64, ki_65, ki_66, \
                         mi_258, mi_259, mi_260, mi_261, mi_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = -2.0 * ki_62[k]
                   + f_0 * mi_258[k];

        t_147[k] = -2.0 * ki_63[k]
                   + f_0 * mi_259[k];

        t_148[k] = -2.0 * ki_64[k]
                   + f_0 * mi_260[k];

        t_149[k] = -2.0 * ki_65[k]
                   + f_0 * mi_261[k];

        t_150[k] = -2.0 * ki_66[k]
                   + f_0 * mi_262[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, ki_67, ki_68, ki_69, ki_70, ki_71, \
                         mi_263, mi_264, mi_265, mi_266, mi_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = -2.0 * ki_67[k]
                   + f_0 * mi_263[k];

        t_152[k] = -2.0 * ki_68[k]
                   + f_0 * mi_264[k];

        t_153[k] = -2.0 * ki_69[k]
                   + f_0 * mi_265[k];

        t_154[k] = -2.0 * ki_70[k]
                   + f_0 * mi_266[k];

        t_155[k] = -2.0 * ki_71[k]
                   + f_0 * mi_267[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, ki_72, ki_73, ki_74, ki_75, ki_76, \
                         mi_268, mi_269, mi_270, mi_271, mi_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = -2.0 * ki_72[k]
                   + f_0 * mi_268[k];

        t_157[k] = -2.0 * ki_73[k]
                   + f_0 * mi_269[k];

        t_158[k] = -2.0 * ki_74[k]
                   + f_0 * mi_270[k];

        t_159[k] = -2.0 * ki_75[k]
                   + f_0 * mi_271[k];

        t_160[k] = -2.0 * ki_76[k]
                   + f_0 * mi_272[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, ki_77, ki_78, ki_79, ki_80, ki_81, \
                         mi_273, mi_274, mi_275, mi_276, mi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = -2.0 * ki_77[k]
                   + f_0 * mi_273[k];

        t_162[k] = -2.0 * ki_78[k]
                   + f_0 * mi_274[k];

        t_163[k] = -2.0 * ki_79[k]
                   + f_0 * mi_275[k];

        t_164[k] = -2.0 * ki_80[k]
                   + f_0 * mi_276[k];

        t_165[k] = -2.0 * ki_81[k]
                   + f_0 * mi_277[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, t_171, t_172, ki_82, ki_83, \
                         mi_278, mi_279, mi_308, mi_309, mi_310, mi_311, \
                         mi_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = -2.0 * ki_82[k]
                   + f_0 * mi_278[k];

        t_167[k] = -2.0 * ki_83[k]
                   + f_0 * mi_279[k];

        t_168[k] = f_0 * mi_308[k];

        t_169[k] = f_0 * mi_309[k];

        t_170[k] = f_0 * mi_310[k];

        t_171[k] = f_0 * mi_311[k];

        t_172[k] = f_0 * mi_312[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, t_179, t_180, mi_313, \
                         mi_314, mi_315, mi_316, mi_317, mi_318, mi_319, \
                         mi_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_0 * mi_313[k];

        t_174[k] = f_0 * mi_314[k];

        t_175[k] = f_0 * mi_315[k];

        t_176[k] = f_0 * mi_316[k];

        t_177[k] = f_0 * mi_317[k];

        t_178[k] = f_0 * mi_318[k];

        t_179[k] = f_0 * mi_319[k];

        t_180[k] = f_0 * mi_320[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_338 = buffer.data(target + 338);
    auto *t_339 = buffer.data(target + 339);
    auto *t_340 = buffer.data(target + 340);
    auto *t_341 = buffer.data(target + 341);
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);

    const auto *ki_84 = buffer.data(ki + 84);
    const auto *ki_85 = buffer.data(ki + 85);
    const auto *ki_86 = buffer.data(ki + 86);
    const auto *ki_87 = buffer.data(ki + 87);
    const auto *ki_88 = buffer.data(ki + 88);
    const auto *ki_89 = buffer.data(ki + 89);
    const auto *ki_90 = buffer.data(ki + 90);
    const auto *ki_91 = buffer.data(ki + 91);
    const auto *ki_92 = buffer.data(ki + 92);
    const auto *ki_93 = buffer.data(ki + 93);
    const auto *ki_94 = buffer.data(ki + 94);
    const auto *ki_95 = buffer.data(ki + 95);
    const auto *ki_96 = buffer.data(ki + 96);
    const auto *ki_97 = buffer.data(ki + 97);
    const auto *ki_98 = buffer.data(ki + 98);
    const auto *ki_99 = buffer.data(ki + 99);
    const auto *ki_100 = buffer.data(ki + 100);
    const auto *ki_101 = buffer.data(ki + 101);
    const auto *ki_102 = buffer.data(ki + 102);
    const auto *ki_103 = buffer.data(ki + 103);
    const auto *ki_104 = buffer.data(ki + 104);
    const auto *ki_105 = buffer.data(ki + 105);
    const auto *ki_106 = buffer.data(ki + 106);
    const auto *ki_107 = buffer.data(ki + 107);
    const auto *ki_108 = buffer.data(ki + 108);
    const auto *ki_109 = buffer.data(ki + 109);
    const auto *ki_110 = buffer.data(ki + 110);
    const auto *ki_111 = buffer.data(ki + 111);
    const auto *ki_112 = buffer.data(ki + 112);
    const auto *ki_113 = buffer.data(ki + 113);
    const auto *ki_114 = buffer.data(ki + 114);
    const auto *ki_115 = buffer.data(ki + 115);
    const auto *ki_116 = buffer.data(ki + 116);
    const auto *ki_117 = buffer.data(ki + 117);
    const auto *ki_118 = buffer.data(ki + 118);
    const auto *ki_119 = buffer.data(ki + 119);
    const auto *ki_120 = buffer.data(ki + 120);
    const auto *ki_121 = buffer.data(ki + 121);
    const auto *ki_122 = buffer.data(ki + 122);
    const auto *ki_123 = buffer.data(ki + 123);
    const auto *ki_124 = buffer.data(ki + 124);
    const auto *ki_125 = buffer.data(ki + 125);
    const auto *ki_126 = buffer.data(ki + 126);
    const auto *ki_127 = buffer.data(ki + 127);
    const auto *ki_128 = buffer.data(ki + 128);
    const auto *ki_129 = buffer.data(ki + 129);
    const auto *ki_130 = buffer.data(ki + 130);
    const auto *ki_131 = buffer.data(ki + 131);
    const auto *ki_132 = buffer.data(ki + 132);
    const auto *ki_133 = buffer.data(ki + 133);
    const auto *ki_134 = buffer.data(ki + 134);
    const auto *ki_135 = buffer.data(ki + 135);
    const auto *ki_136 = buffer.data(ki + 136);
    const auto *ki_137 = buffer.data(ki + 137);
    const auto *ki_138 = buffer.data(ki + 138);
    const auto *ki_139 = buffer.data(ki + 139);
    const auto *ki_140 = buffer.data(ki + 140);
    const auto *ki_141 = buffer.data(ki + 141);
    const auto *ki_142 = buffer.data(ki + 142);
    const auto *ki_143 = buffer.data(ki + 143);
    const auto *ki_144 = buffer.data(ki + 144);
    const auto *ki_145 = buffer.data(ki + 145);
    const auto *ki_146 = buffer.data(ki + 146);
    const auto *ki_147 = buffer.data(ki + 147);
    const auto *ki_148 = buffer.data(ki + 148);
    const auto *ki_149 = buffer.data(ki + 149);
    const auto *ki_150 = buffer.data(ki + 150);
    const auto *ki_151 = buffer.data(ki + 151);
    const auto *ki_152 = buffer.data(ki + 152);
    const auto *ki_153 = buffer.data(ki + 153);
    const auto *ki_154 = buffer.data(ki + 154);
    const auto *ki_155 = buffer.data(ki + 155);
    const auto *ki_156 = buffer.data(ki + 156);
    const auto *ki_157 = buffer.data(ki + 157);
    const auto *ki_158 = buffer.data(ki + 158);
    const auto *ki_159 = buffer.data(ki + 159);
    const auto *ki_160 = buffer.data(ki + 160);
    const auto *ki_161 = buffer.data(ki + 161);
    const auto *ki_162 = buffer.data(ki + 162);
    const auto *ki_163 = buffer.data(ki + 163);
    const auto *ki_164 = buffer.data(ki + 164);
    const auto *ki_165 = buffer.data(ki + 165);
    const auto *ki_166 = buffer.data(ki + 166);
    const auto *ki_167 = buffer.data(ki + 167);
    const auto *ki_168 = buffer.data(ki + 168);
    const auto *ki_169 = buffer.data(ki + 169);
    const auto *ki_170 = buffer.data(ki + 170);
    const auto *ki_171 = buffer.data(ki + 171);
    const auto *ki_172 = buffer.data(ki + 172);
    const auto *ki_173 = buffer.data(ki + 173);
    const auto *ki_174 = buffer.data(ki + 174);
    const auto *ki_175 = buffer.data(ki + 175);
    const auto *ki_176 = buffer.data(ki + 176);
    const auto *ki_177 = buffer.data(ki + 177);
    const auto *ki_178 = buffer.data(ki + 178);
    const auto *ki_179 = buffer.data(ki + 179);
    const auto *ki_180 = buffer.data(ki + 180);
    const auto *ki_181 = buffer.data(ki + 181);
    const auto *ki_182 = buffer.data(ki + 182);
    const auto *ki_183 = buffer.data(ki + 183);
    const auto *ki_184 = buffer.data(ki + 184);
    const auto *ki_185 = buffer.data(ki + 185);
    const auto *ki_186 = buffer.data(ki + 186);
    const auto *ki_187 = buffer.data(ki + 187);
    const auto *ki_188 = buffer.data(ki + 188);
    const auto *ki_189 = buffer.data(ki + 189);
    const auto *ki_190 = buffer.data(ki + 190);
    const auto *ki_191 = buffer.data(ki + 191);
    const auto *ki_192 = buffer.data(ki + 192);
    const auto *ki_193 = buffer.data(ki + 193);
    const auto *ki_194 = buffer.data(ki + 194);
    const auto *ki_195 = buffer.data(ki + 195);
    const auto *ki_196 = buffer.data(ki + 196);
    const auto *ki_197 = buffer.data(ki + 197);
    const auto *ki_198 = buffer.data(ki + 198);
    const auto *ki_199 = buffer.data(ki + 199);
    const auto *ki_200 = buffer.data(ki + 200);
    const auto *ki_201 = buffer.data(ki + 201);
    const auto *ki_202 = buffer.data(ki + 202);
    const auto *ki_203 = buffer.data(ki + 203);
    const auto *ki_204 = buffer.data(ki + 204);
    const auto *ki_205 = buffer.data(ki + 205);
    const auto *ki_206 = buffer.data(ki + 206);

    const auto *mi_321 = buffer.data(mi + 321);
    const auto *mi_322 = buffer.data(mi + 322);
    const auto *mi_323 = buffer.data(mi + 323);
    const auto *mi_324 = buffer.data(mi + 324);
    const auto *mi_325 = buffer.data(mi + 325);
    const auto *mi_326 = buffer.data(mi + 326);
    const auto *mi_327 = buffer.data(mi + 327);
    const auto *mi_328 = buffer.data(mi + 328);
    const auto *mi_329 = buffer.data(mi + 329);
    const auto *mi_330 = buffer.data(mi + 330);
    const auto *mi_331 = buffer.data(mi + 331);
    const auto *mi_332 = buffer.data(mi + 332);
    const auto *mi_333 = buffer.data(mi + 333);
    const auto *mi_334 = buffer.data(mi + 334);
    const auto *mi_335 = buffer.data(mi + 335);
    const auto *mi_336 = buffer.data(mi + 336);
    const auto *mi_337 = buffer.data(mi + 337);
    const auto *mi_338 = buffer.data(mi + 338);
    const auto *mi_339 = buffer.data(mi + 339);
    const auto *mi_340 = buffer.data(mi + 340);
    const auto *mi_341 = buffer.data(mi + 341);
    const auto *mi_342 = buffer.data(mi + 342);
    const auto *mi_343 = buffer.data(mi + 343);
    const auto *mi_344 = buffer.data(mi + 344);
    const auto *mi_345 = buffer.data(mi + 345);
    const auto *mi_346 = buffer.data(mi + 346);
    const auto *mi_347 = buffer.data(mi + 347);
    const auto *mi_348 = buffer.data(mi + 348);
    const auto *mi_349 = buffer.data(mi + 349);
    const auto *mi_350 = buffer.data(mi + 350);
    const auto *mi_351 = buffer.data(mi + 351);
    const auto *mi_352 = buffer.data(mi + 352);
    const auto *mi_353 = buffer.data(mi + 353);
    const auto *mi_354 = buffer.data(mi + 354);
    const auto *mi_355 = buffer.data(mi + 355);
    const auto *mi_356 = buffer.data(mi + 356);
    const auto *mi_357 = buffer.data(mi + 357);
    const auto *mi_358 = buffer.data(mi + 358);
    const auto *mi_359 = buffer.data(mi + 359);
    const auto *mi_360 = buffer.data(mi + 360);
    const auto *mi_361 = buffer.data(mi + 361);
    const auto *mi_362 = buffer.data(mi + 362);
    const auto *mi_363 = buffer.data(mi + 363);
    const auto *mi_364 = buffer.data(mi + 364);
    const auto *mi_365 = buffer.data(mi + 365);
    const auto *mi_366 = buffer.data(mi + 366);
    const auto *mi_367 = buffer.data(mi + 367);
    const auto *mi_368 = buffer.data(mi + 368);
    const auto *mi_369 = buffer.data(mi + 369);
    const auto *mi_370 = buffer.data(mi + 370);
    const auto *mi_371 = buffer.data(mi + 371);
    const auto *mi_372 = buffer.data(mi + 372);
    const auto *mi_373 = buffer.data(mi + 373);
    const auto *mi_374 = buffer.data(mi + 374);
    const auto *mi_375 = buffer.data(mi + 375);
    const auto *mi_376 = buffer.data(mi + 376);
    const auto *mi_377 = buffer.data(mi + 377);
    const auto *mi_378 = buffer.data(mi + 378);
    const auto *mi_379 = buffer.data(mi + 379);
    const auto *mi_380 = buffer.data(mi + 380);
    const auto *mi_381 = buffer.data(mi + 381);
    const auto *mi_382 = buffer.data(mi + 382);
    const auto *mi_383 = buffer.data(mi + 383);
    const auto *mi_384 = buffer.data(mi + 384);
    const auto *mi_385 = buffer.data(mi + 385);
    const auto *mi_386 = buffer.data(mi + 386);
    const auto *mi_387 = buffer.data(mi + 387);
    const auto *mi_388 = buffer.data(mi + 388);
    const auto *mi_389 = buffer.data(mi + 389);
    const auto *mi_390 = buffer.data(mi + 390);
    const auto *mi_391 = buffer.data(mi + 391);
    const auto *mi_392 = buffer.data(mi + 392);
    const auto *mi_393 = buffer.data(mi + 393);
    const auto *mi_394 = buffer.data(mi + 394);
    const auto *mi_395 = buffer.data(mi + 395);
    const auto *mi_396 = buffer.data(mi + 396);
    const auto *mi_397 = buffer.data(mi + 397);
    const auto *mi_398 = buffer.data(mi + 398);
    const auto *mi_399 = buffer.data(mi + 399);
    const auto *mi_400 = buffer.data(mi + 400);
    const auto *mi_401 = buffer.data(mi + 401);
    const auto *mi_402 = buffer.data(mi + 402);
    const auto *mi_403 = buffer.data(mi + 403);
    const auto *mi_404 = buffer.data(mi + 404);
    const auto *mi_405 = buffer.data(mi + 405);
    const auto *mi_406 = buffer.data(mi + 406);
    const auto *mi_407 = buffer.data(mi + 407);
    const auto *mi_408 = buffer.data(mi + 408);
    const auto *mi_409 = buffer.data(mi + 409);
    const auto *mi_410 = buffer.data(mi + 410);
    const auto *mi_411 = buffer.data(mi + 411);
    const auto *mi_412 = buffer.data(mi + 412);
    const auto *mi_413 = buffer.data(mi + 413);
    const auto *mi_414 = buffer.data(mi + 414);
    const auto *mi_415 = buffer.data(mi + 415);
    const auto *mi_416 = buffer.data(mi + 416);
    const auto *mi_417 = buffer.data(mi + 417);
    const auto *mi_418 = buffer.data(mi + 418);
    const auto *mi_419 = buffer.data(mi + 419);
    const auto *mi_448 = buffer.data(mi + 448);
    const auto *mi_449 = buffer.data(mi + 449);
    const auto *mi_450 = buffer.data(mi + 450);
    const auto *mi_451 = buffer.data(mi + 451);
    const auto *mi_452 = buffer.data(mi + 452);
    const auto *mi_453 = buffer.data(mi + 453);
    const auto *mi_454 = buffer.data(mi + 454);
    const auto *mi_455 = buffer.data(mi + 455);
    const auto *mi_456 = buffer.data(mi + 456);
    const auto *mi_457 = buffer.data(mi + 457);
    const auto *mi_458 = buffer.data(mi + 458);
    const auto *mi_459 = buffer.data(mi + 459);
    const auto *mi_460 = buffer.data(mi + 460);
    const auto *mi_461 = buffer.data(mi + 461);
    const auto *mi_462 = buffer.data(mi + 462);
    const auto *mi_463 = buffer.data(mi + 463);
    const auto *mi_464 = buffer.data(mi + 464);
    const auto *mi_465 = buffer.data(mi + 465);
    const auto *mi_466 = buffer.data(mi + 466);
    const auto *mi_467 = buffer.data(mi + 467);
    const auto *mi_468 = buffer.data(mi + 468);
    const auto *mi_469 = buffer.data(mi + 469);
    const auto *mi_470 = buffer.data(mi + 470);
    const auto *mi_471 = buffer.data(mi + 471);
    const auto *mi_472 = buffer.data(mi + 472);
    const auto *mi_473 = buffer.data(mi + 473);
    const auto *mi_474 = buffer.data(mi + 474);
    const auto *mi_475 = buffer.data(mi + 475);
    const auto *mi_476 = buffer.data(mi + 476);
    const auto *mi_477 = buffer.data(mi + 477);
    const auto *mi_478 = buffer.data(mi + 478);
    const auto *mi_479 = buffer.data(mi + 479);
    const auto *mi_480 = buffer.data(mi + 480);
    const auto *mi_481 = buffer.data(mi + 481);
    const auto *mi_482 = buffer.data(mi + 482);
    const auto *mi_483 = buffer.data(mi + 483);
    const auto *mi_484 = buffer.data(mi + 484);
    const auto *mi_485 = buffer.data(mi + 485);
    const auto *mi_486 = buffer.data(mi + 486);
    const auto *mi_487 = buffer.data(mi + 487);
    const auto *mi_488 = buffer.data(mi + 488);
    const auto *mi_489 = buffer.data(mi + 489);
    const auto *mi_490 = buffer.data(mi + 490);
    const auto *mi_491 = buffer.data(mi + 491);
    const auto *mi_492 = buffer.data(mi + 492);
    const auto *mi_493 = buffer.data(mi + 493);
    const auto *mi_494 = buffer.data(mi + 494);
    const auto *mi_495 = buffer.data(mi + 495);
    const auto *mi_496 = buffer.data(mi + 496);
    const auto *mi_497 = buffer.data(mi + 497);
    const auto *mi_498 = buffer.data(mi + 498);
    const auto *mi_499 = buffer.data(mi + 499);
    const auto *mi_500 = buffer.data(mi + 500);
    const auto *mi_501 = buffer.data(mi + 501);
    const auto *mi_502 = buffer.data(mi + 502);
    const auto *mi_503 = buffer.data(mi + 503);
    const auto *mi_504 = buffer.data(mi + 504);
    const auto *mi_505 = buffer.data(mi + 505);
    const auto *mi_506 = buffer.data(mi + 506);
    const auto *mi_507 = buffer.data(mi + 507);
    const auto *mi_508 = buffer.data(mi + 508);
    const auto *mi_509 = buffer.data(mi + 509);
    const auto *mi_510 = buffer.data(mi + 510);
    const auto *mi_511 = buffer.data(mi + 511);
    const auto *mi_512 = buffer.data(mi + 512);
    const auto *mi_513 = buffer.data(mi + 513);
    const auto *mi_514 = buffer.data(mi + 514);

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, t_186, t_187, t_188, mi_321, \
                         mi_322, mi_323, mi_324, mi_325, mi_326, mi_327, \
                         mi_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_0 * mi_321[k];

        t_182[k] = f_0 * mi_322[k];

        t_183[k] = f_0 * mi_323[k];

        t_184[k] = f_0 * mi_324[k];

        t_185[k] = f_0 * mi_325[k];

        t_186[k] = f_0 * mi_326[k];

        t_187[k] = f_0 * mi_327[k];

        t_188[k] = f_0 * mi_328[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, mi_329, mi_330, \
                         mi_331, mi_332, mi_333, mi_334, mi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * mi_329[k];

        t_190[k] = f_0 * mi_330[k];

        t_191[k] = f_0 * mi_331[k];

        t_192[k] = f_0 * mi_332[k];

        t_193[k] = f_0 * mi_333[k];

        t_194[k] = f_0 * mi_334[k];

        t_195[k] = f_0 * mi_335[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, ki_84, ki_85, ki_86, ki_87, ki_88, \
                         mi_336, mi_337, mi_338, mi_339, mi_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -ki_84[k]
                   + f_0 * mi_336[k];

        t_197[k] = -ki_85[k]
                   + f_0 * mi_337[k];

        t_198[k] = -ki_86[k]
                   + f_0 * mi_338[k];

        t_199[k] = -ki_87[k]
                   + f_0 * mi_339[k];

        t_200[k] = -ki_88[k]
                   + f_0 * mi_340[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, ki_89, ki_90, ki_91, ki_92, ki_93, \
                         mi_341, mi_342, mi_343, mi_344, mi_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = -ki_89[k]
                   + f_0 * mi_341[k];

        t_202[k] = -ki_90[k]
                   + f_0 * mi_342[k];

        t_203[k] = -ki_91[k]
                   + f_0 * mi_343[k];

        t_204[k] = -ki_92[k]
                   + f_0 * mi_344[k];

        t_205[k] = -ki_93[k]
                   + f_0 * mi_345[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, ki_94, ki_95, ki_96, ki_97, ki_98, \
                         mi_346, mi_347, mi_348, mi_349, mi_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = -ki_94[k]
                   + f_0 * mi_346[k];

        t_207[k] = -ki_95[k]
                   + f_0 * mi_347[k];

        t_208[k] = -ki_96[k]
                   + f_0 * mi_348[k];

        t_209[k] = -ki_97[k]
                   + f_0 * mi_349[k];

        t_210[k] = -ki_98[k]
                   + f_0 * mi_350[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, ki_99, ki_100, ki_101, ki_102, \
                         ki_103, mi_351, mi_352, mi_353, mi_354, \
                         mi_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = -ki_99[k]
                   + f_0 * mi_351[k];

        t_212[k] = -ki_100[k]
                   + f_0 * mi_352[k];

        t_213[k] = -ki_101[k]
                   + f_0 * mi_353[k];

        t_214[k] = -ki_102[k]
                   + f_0 * mi_354[k];

        t_215[k] = -ki_103[k]
                   + f_0 * mi_355[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, ki_104, ki_105, ki_106, ki_107, \
                         ki_108, mi_356, mi_357, mi_358, mi_359, \
                         mi_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = -ki_104[k]
                   + f_0 * mi_356[k];

        t_217[k] = -ki_105[k]
                   + f_0 * mi_357[k];

        t_218[k] = -ki_106[k]
                   + f_0 * mi_358[k];

        t_219[k] = -ki_107[k]
                   + f_0 * mi_359[k];

        t_220[k] = -ki_108[k]
                   + f_0 * mi_360[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, ki_109, ki_110, ki_111, ki_112, \
                         ki_113, mi_361, mi_362, mi_363, mi_364, \
                         mi_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = -ki_109[k]
                   + f_0 * mi_361[k];

        t_222[k] = -ki_110[k]
                   + f_0 * mi_362[k];

        t_223[k] = -ki_111[k]
                   + f_0 * mi_363[k];

        t_224[k] = -2.0 * ki_112[k]
                   + f_0 * mi_364[k];

        t_225[k] = -2.0 * ki_113[k]
                   + f_0 * mi_365[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, t_230, ki_114, ki_115, ki_116, ki_117, \
                         ki_118, mi_366, mi_367, mi_368, mi_369, \
                         mi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = -2.0 * ki_114[k]
                   + f_0 * mi_366[k];

        t_227[k] = -2.0 * ki_115[k]
                   + f_0 * mi_367[k];

        t_228[k] = -2.0 * ki_116[k]
                   + f_0 * mi_368[k];

        t_229[k] = -2.0 * ki_117[k]
                   + f_0 * mi_369[k];

        t_230[k] = -2.0 * ki_118[k]
                   + f_0 * mi_370[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, ki_119, ki_120, ki_121, ki_122, \
                         ki_123, mi_371, mi_372, mi_373, mi_374, \
                         mi_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = -2.0 * ki_119[k]
                   + f_0 * mi_371[k];

        t_232[k] = -2.0 * ki_120[k]
                   + f_0 * mi_372[k];

        t_233[k] = -2.0 * ki_121[k]
                   + f_0 * mi_373[k];

        t_234[k] = -2.0 * ki_122[k]
                   + f_0 * mi_374[k];

        t_235[k] = -2.0 * ki_123[k]
                   + f_0 * mi_375[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, ki_124, ki_125, ki_126, ki_127, \
                         ki_128, mi_376, mi_377, mi_378, mi_379, \
                         mi_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = -2.0 * ki_124[k]
                   + f_0 * mi_376[k];

        t_237[k] = -2.0 * ki_125[k]
                   + f_0 * mi_377[k];

        t_238[k] = -2.0 * ki_126[k]
                   + f_0 * mi_378[k];

        t_239[k] = -2.0 * ki_127[k]
                   + f_0 * mi_379[k];

        t_240[k] = -2.0 * ki_128[k]
                   + f_0 * mi_380[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, ki_129, ki_130, ki_131, ki_132, \
                         ki_133, mi_381, mi_382, mi_383, mi_384, \
                         mi_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = -2.0 * ki_129[k]
                   + f_0 * mi_381[k];

        t_242[k] = -2.0 * ki_130[k]
                   + f_0 * mi_382[k];

        t_243[k] = -2.0 * ki_131[k]
                   + f_0 * mi_383[k];

        t_244[k] = -2.0 * ki_132[k]
                   + f_0 * mi_384[k];

        t_245[k] = -2.0 * ki_133[k]
                   + f_0 * mi_385[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, ki_134, ki_135, ki_136, ki_137, \
                         ki_138, mi_386, mi_387, mi_388, mi_389, \
                         mi_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = -2.0 * ki_134[k]
                   + f_0 * mi_386[k];

        t_247[k] = -2.0 * ki_135[k]
                   + f_0 * mi_387[k];

        t_248[k] = -2.0 * ki_136[k]
                   + f_0 * mi_388[k];

        t_249[k] = -2.0 * ki_137[k]
                   + f_0 * mi_389[k];

        t_250[k] = -2.0 * ki_138[k]
                   + f_0 * mi_390[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, ki_139, ki_140, ki_141, ki_142, \
                         ki_143, mi_391, mi_392, mi_393, mi_394, \
                         mi_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = -2.0 * ki_139[k]
                   + f_0 * mi_391[k];

        t_252[k] = -3.0 * ki_140[k]
                   + f_0 * mi_392[k];

        t_253[k] = -3.0 * ki_141[k]
                   + f_0 * mi_393[k];

        t_254[k] = -3.0 * ki_142[k]
                   + f_0 * mi_394[k];

        t_255[k] = -3.0 * ki_143[k]
                   + f_0 * mi_395[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, ki_144, ki_145, ki_146, ki_147, \
                         ki_148, mi_396, mi_397, mi_398, mi_399, \
                         mi_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = -3.0 * ki_144[k]
                   + f_0 * mi_396[k];

        t_257[k] = -3.0 * ki_145[k]
                   + f_0 * mi_397[k];

        t_258[k] = -3.0 * ki_146[k]
                   + f_0 * mi_398[k];

        t_259[k] = -3.0 * ki_147[k]
                   + f_0 * mi_399[k];

        t_260[k] = -3.0 * ki_148[k]
                   + f_0 * mi_400[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, ki_149, ki_150, ki_151, ki_152, \
                         ki_153, mi_401, mi_402, mi_403, mi_404, \
                         mi_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = -3.0 * ki_149[k]
                   + f_0 * mi_401[k];

        t_262[k] = -3.0 * ki_150[k]
                   + f_0 * mi_402[k];

        t_263[k] = -3.0 * ki_151[k]
                   + f_0 * mi_403[k];

        t_264[k] = -3.0 * ki_152[k]
                   + f_0 * mi_404[k];

        t_265[k] = -3.0 * ki_153[k]
                   + f_0 * mi_405[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, ki_154, ki_155, ki_156, ki_157, \
                         ki_158, mi_406, mi_407, mi_408, mi_409, \
                         mi_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = -3.0 * ki_154[k]
                   + f_0 * mi_406[k];

        t_267[k] = -3.0 * ki_155[k]
                   + f_0 * mi_407[k];

        t_268[k] = -3.0 * ki_156[k]
                   + f_0 * mi_408[k];

        t_269[k] = -3.0 * ki_157[k]
                   + f_0 * mi_409[k];

        t_270[k] = -3.0 * ki_158[k]
                   + f_0 * mi_410[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, ki_159, ki_160, ki_161, ki_162, \
                         ki_163, mi_411, mi_412, mi_413, mi_414, \
                         mi_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = -3.0 * ki_159[k]
                   + f_0 * mi_411[k];

        t_272[k] = -3.0 * ki_160[k]
                   + f_0 * mi_412[k];

        t_273[k] = -3.0 * ki_161[k]
                   + f_0 * mi_413[k];

        t_274[k] = -3.0 * ki_162[k]
                   + f_0 * mi_414[k];

        t_275[k] = -3.0 * ki_163[k]
                   + f_0 * mi_415[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, t_281, ki_164, ki_165, ki_166, \
                         ki_167, mi_416, mi_417, mi_418, mi_419, mi_448, \
                         mi_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = -3.0 * ki_164[k]
                   + f_0 * mi_416[k];

        t_277[k] = -3.0 * ki_165[k]
                   + f_0 * mi_417[k];

        t_278[k] = -3.0 * ki_166[k]
                   + f_0 * mi_418[k];

        t_279[k] = -3.0 * ki_167[k]
                   + f_0 * mi_419[k];

        t_280[k] = f_0 * mi_448[k];

        t_281[k] = f_0 * mi_449[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, t_287, t_288, t_289, mi_450, \
                         mi_451, mi_452, mi_453, mi_454, mi_455, mi_456, \
                         mi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_0 * mi_450[k];

        t_283[k] = f_0 * mi_451[k];

        t_284[k] = f_0 * mi_452[k];

        t_285[k] = f_0 * mi_453[k];

        t_286[k] = f_0 * mi_454[k];

        t_287[k] = f_0 * mi_455[k];

        t_288[k] = f_0 * mi_456[k];

        t_289[k] = f_0 * mi_457[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, t_295, t_296, t_297, mi_458, \
                         mi_459, mi_460, mi_461, mi_462, mi_463, mi_464, \
                         mi_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_0 * mi_458[k];

        t_291[k] = f_0 * mi_459[k];

        t_292[k] = f_0 * mi_460[k];

        t_293[k] = f_0 * mi_461[k];

        t_294[k] = f_0 * mi_462[k];

        t_295[k] = f_0 * mi_463[k];

        t_296[k] = f_0 * mi_464[k];

        t_297[k] = f_0 * mi_465[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, t_303, t_304, t_305, mi_466, \
                         mi_467, mi_468, mi_469, mi_470, mi_471, mi_472, \
                         mi_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_0 * mi_466[k];

        t_299[k] = f_0 * mi_467[k];

        t_300[k] = f_0 * mi_468[k];

        t_301[k] = f_0 * mi_469[k];

        t_302[k] = f_0 * mi_470[k];

        t_303[k] = f_0 * mi_471[k];

        t_304[k] = f_0 * mi_472[k];

        t_305[k] = f_0 * mi_473[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, t_311, ki_168, ki_169, ki_170, \
                         ki_171, mi_474, mi_475, mi_476, mi_477, mi_478, \
                         mi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_0 * mi_474[k];

        t_307[k] = f_0 * mi_475[k];

        t_308[k] = -ki_168[k]
                   + f_0 * mi_476[k];

        t_309[k] = -ki_169[k]
                   + f_0 * mi_477[k];

        t_310[k] = -ki_170[k]
                   + f_0 * mi_478[k];

        t_311[k] = -ki_171[k]
                   + f_0 * mi_479[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, ki_172, ki_173, ki_174, ki_175, \
                         ki_176, mi_480, mi_481, mi_482, mi_483, \
                         mi_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -ki_172[k]
                   + f_0 * mi_480[k];

        t_313[k] = -ki_173[k]
                   + f_0 * mi_481[k];

        t_314[k] = -ki_174[k]
                   + f_0 * mi_482[k];

        t_315[k] = -ki_175[k]
                   + f_0 * mi_483[k];

        t_316[k] = -ki_176[k]
                   + f_0 * mi_484[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, ki_177, ki_178, ki_179, ki_180, \
                         ki_181, mi_485, mi_486, mi_487, mi_488, \
                         mi_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = -ki_177[k]
                   + f_0 * mi_485[k];

        t_318[k] = -ki_178[k]
                   + f_0 * mi_486[k];

        t_319[k] = -ki_179[k]
                   + f_0 * mi_487[k];

        t_320[k] = -ki_180[k]
                   + f_0 * mi_488[k];

        t_321[k] = -ki_181[k]
                   + f_0 * mi_489[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, ki_182, ki_183, ki_184, ki_185, \
                         ki_186, mi_490, mi_491, mi_492, mi_493, \
                         mi_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = -ki_182[k]
                   + f_0 * mi_490[k];

        t_323[k] = -ki_183[k]
                   + f_0 * mi_491[k];

        t_324[k] = -ki_184[k]
                   + f_0 * mi_492[k];

        t_325[k] = -ki_185[k]
                   + f_0 * mi_493[k];

        t_326[k] = -ki_186[k]
                   + f_0 * mi_494[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, ki_187, ki_188, ki_189, ki_190, \
                         ki_191, mi_495, mi_496, mi_497, mi_498, \
                         mi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = -ki_187[k]
                   + f_0 * mi_495[k];

        t_328[k] = -ki_188[k]
                   + f_0 * mi_496[k];

        t_329[k] = -ki_189[k]
                   + f_0 * mi_497[k];

        t_330[k] = -ki_190[k]
                   + f_0 * mi_498[k];

        t_331[k] = -ki_191[k]
                   + f_0 * mi_499[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, ki_192, ki_193, ki_194, ki_195, \
                         ki_196, mi_500, mi_501, mi_502, mi_503, \
                         mi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -ki_192[k]
                   + f_0 * mi_500[k];

        t_333[k] = -ki_193[k]
                   + f_0 * mi_501[k];

        t_334[k] = -ki_194[k]
                   + f_0 * mi_502[k];

        t_335[k] = -ki_195[k]
                   + f_0 * mi_503[k];

        t_336[k] = -2.0 * ki_196[k]
                   + f_0 * mi_504[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, ki_197, ki_198, ki_199, ki_200, \
                         ki_201, mi_505, mi_506, mi_507, mi_508, \
                         mi_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -2.0 * ki_197[k]
                   + f_0 * mi_505[k];

        t_338[k] = -2.0 * ki_198[k]
                   + f_0 * mi_506[k];

        t_339[k] = -2.0 * ki_199[k]
                   + f_0 * mi_507[k];

        t_340[k] = -2.0 * ki_200[k]
                   + f_0 * mi_508[k];

        t_341[k] = -2.0 * ki_201[k]
                   + f_0 * mi_509[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, ki_202, ki_203, ki_204, ki_205, \
                         ki_206, mi_510, mi_511, mi_512, mi_513, \
                         mi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -2.0 * ki_202[k]
                   + f_0 * mi_510[k];

        t_343[k] = -2.0 * ki_203[k]
                   + f_0 * mi_511[k];

        t_344[k] = -2.0 * ki_204[k]
                   + f_0 * mi_512[k];

        t_345[k] = -2.0 * ki_205[k]
                   + f_0 * mi_513[k];

        t_346[k] = -2.0 * ki_206[k]
                   + f_0 * mi_514[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ki_207 = buffer.data(ki + 207);
    const auto *ki_208 = buffer.data(ki + 208);
    const auto *ki_209 = buffer.data(ki + 209);
    const auto *ki_210 = buffer.data(ki + 210);
    const auto *ki_211 = buffer.data(ki + 211);
    const auto *ki_212 = buffer.data(ki + 212);
    const auto *ki_213 = buffer.data(ki + 213);
    const auto *ki_214 = buffer.data(ki + 214);
    const auto *ki_215 = buffer.data(ki + 215);
    const auto *ki_216 = buffer.data(ki + 216);
    const auto *ki_217 = buffer.data(ki + 217);
    const auto *ki_218 = buffer.data(ki + 218);
    const auto *ki_219 = buffer.data(ki + 219);
    const auto *ki_220 = buffer.data(ki + 220);
    const auto *ki_221 = buffer.data(ki + 221);
    const auto *ki_222 = buffer.data(ki + 222);
    const auto *ki_223 = buffer.data(ki + 223);
    const auto *ki_224 = buffer.data(ki + 224);
    const auto *ki_225 = buffer.data(ki + 225);
    const auto *ki_226 = buffer.data(ki + 226);
    const auto *ki_227 = buffer.data(ki + 227);
    const auto *ki_228 = buffer.data(ki + 228);
    const auto *ki_229 = buffer.data(ki + 229);
    const auto *ki_230 = buffer.data(ki + 230);
    const auto *ki_231 = buffer.data(ki + 231);
    const auto *ki_232 = buffer.data(ki + 232);
    const auto *ki_233 = buffer.data(ki + 233);
    const auto *ki_234 = buffer.data(ki + 234);
    const auto *ki_235 = buffer.data(ki + 235);
    const auto *ki_236 = buffer.data(ki + 236);
    const auto *ki_237 = buffer.data(ki + 237);
    const auto *ki_238 = buffer.data(ki + 238);
    const auto *ki_239 = buffer.data(ki + 239);
    const auto *ki_240 = buffer.data(ki + 240);
    const auto *ki_241 = buffer.data(ki + 241);
    const auto *ki_242 = buffer.data(ki + 242);
    const auto *ki_243 = buffer.data(ki + 243);
    const auto *ki_244 = buffer.data(ki + 244);
    const auto *ki_245 = buffer.data(ki + 245);
    const auto *ki_246 = buffer.data(ki + 246);
    const auto *ki_247 = buffer.data(ki + 247);
    const auto *ki_248 = buffer.data(ki + 248);
    const auto *ki_249 = buffer.data(ki + 249);
    const auto *ki_250 = buffer.data(ki + 250);
    const auto *ki_251 = buffer.data(ki + 251);
    const auto *ki_252 = buffer.data(ki + 252);
    const auto *ki_253 = buffer.data(ki + 253);
    const auto *ki_254 = buffer.data(ki + 254);
    const auto *ki_255 = buffer.data(ki + 255);
    const auto *ki_256 = buffer.data(ki + 256);
    const auto *ki_257 = buffer.data(ki + 257);
    const auto *ki_258 = buffer.data(ki + 258);
    const auto *ki_259 = buffer.data(ki + 259);
    const auto *ki_260 = buffer.data(ki + 260);
    const auto *ki_261 = buffer.data(ki + 261);
    const auto *ki_262 = buffer.data(ki + 262);
    const auto *ki_263 = buffer.data(ki + 263);
    const auto *ki_264 = buffer.data(ki + 264);
    const auto *ki_265 = buffer.data(ki + 265);
    const auto *ki_266 = buffer.data(ki + 266);
    const auto *ki_267 = buffer.data(ki + 267);
    const auto *ki_268 = buffer.data(ki + 268);
    const auto *ki_269 = buffer.data(ki + 269);
    const auto *ki_270 = buffer.data(ki + 270);
    const auto *ki_271 = buffer.data(ki + 271);
    const auto *ki_272 = buffer.data(ki + 272);
    const auto *ki_273 = buffer.data(ki + 273);
    const auto *ki_274 = buffer.data(ki + 274);
    const auto *ki_275 = buffer.data(ki + 275);
    const auto *ki_276 = buffer.data(ki + 276);
    const auto *ki_277 = buffer.data(ki + 277);
    const auto *ki_278 = buffer.data(ki + 278);
    const auto *ki_279 = buffer.data(ki + 279);
    const auto *ki_280 = buffer.data(ki + 280);
    const auto *ki_281 = buffer.data(ki + 281);
    const auto *ki_282 = buffer.data(ki + 282);
    const auto *ki_283 = buffer.data(ki + 283);
    const auto *ki_284 = buffer.data(ki + 284);
    const auto *ki_285 = buffer.data(ki + 285);
    const auto *ki_286 = buffer.data(ki + 286);
    const auto *ki_287 = buffer.data(ki + 287);
    const auto *ki_288 = buffer.data(ki + 288);
    const auto *ki_289 = buffer.data(ki + 289);
    const auto *ki_290 = buffer.data(ki + 290);
    const auto *ki_291 = buffer.data(ki + 291);
    const auto *ki_292 = buffer.data(ki + 292);
    const auto *ki_293 = buffer.data(ki + 293);
    const auto *ki_294 = buffer.data(ki + 294);
    const auto *ki_295 = buffer.data(ki + 295);
    const auto *ki_296 = buffer.data(ki + 296);
    const auto *ki_297 = buffer.data(ki + 297);
    const auto *ki_298 = buffer.data(ki + 298);
    const auto *ki_299 = buffer.data(ki + 299);
    const auto *ki_300 = buffer.data(ki + 300);
    const auto *ki_301 = buffer.data(ki + 301);
    const auto *ki_302 = buffer.data(ki + 302);
    const auto *ki_303 = buffer.data(ki + 303);
    const auto *ki_304 = buffer.data(ki + 304);
    const auto *ki_305 = buffer.data(ki + 305);
    const auto *ki_306 = buffer.data(ki + 306);
    const auto *ki_307 = buffer.data(ki + 307);
    const auto *ki_308 = buffer.data(ki + 308);
    const auto *ki_309 = buffer.data(ki + 309);
    const auto *ki_310 = buffer.data(ki + 310);
    const auto *ki_311 = buffer.data(ki + 311);
    const auto *ki_312 = buffer.data(ki + 312);
    const auto *ki_313 = buffer.data(ki + 313);
    const auto *ki_314 = buffer.data(ki + 314);
    const auto *ki_315 = buffer.data(ki + 315);
    const auto *ki_316 = buffer.data(ki + 316);
    const auto *ki_317 = buffer.data(ki + 317);
    const auto *ki_318 = buffer.data(ki + 318);
    const auto *ki_319 = buffer.data(ki + 319);
    const auto *ki_320 = buffer.data(ki + 320);
    const auto *ki_321 = buffer.data(ki + 321);
    const auto *ki_322 = buffer.data(ki + 322);
    const auto *ki_323 = buffer.data(ki + 323);
    const auto *ki_324 = buffer.data(ki + 324);
    const auto *ki_325 = buffer.data(ki + 325);
    const auto *ki_326 = buffer.data(ki + 326);
    const auto *ki_327 = buffer.data(ki + 327);
    const auto *ki_328 = buffer.data(ki + 328);
    const auto *ki_329 = buffer.data(ki + 329);
    const auto *ki_330 = buffer.data(ki + 330);
    const auto *ki_331 = buffer.data(ki + 331);
    const auto *ki_332 = buffer.data(ki + 332);
    const auto *ki_333 = buffer.data(ki + 333);
    const auto *ki_334 = buffer.data(ki + 334);
    const auto *ki_335 = buffer.data(ki + 335);
    const auto *ki_336 = buffer.data(ki + 336);
    const auto *ki_337 = buffer.data(ki + 337);
    const auto *ki_338 = buffer.data(ki + 338);

    const auto *mi_515 = buffer.data(mi + 515);
    const auto *mi_516 = buffer.data(mi + 516);
    const auto *mi_517 = buffer.data(mi + 517);
    const auto *mi_518 = buffer.data(mi + 518);
    const auto *mi_519 = buffer.data(mi + 519);
    const auto *mi_520 = buffer.data(mi + 520);
    const auto *mi_521 = buffer.data(mi + 521);
    const auto *mi_522 = buffer.data(mi + 522);
    const auto *mi_523 = buffer.data(mi + 523);
    const auto *mi_524 = buffer.data(mi + 524);
    const auto *mi_525 = buffer.data(mi + 525);
    const auto *mi_526 = buffer.data(mi + 526);
    const auto *mi_527 = buffer.data(mi + 527);
    const auto *mi_528 = buffer.data(mi + 528);
    const auto *mi_529 = buffer.data(mi + 529);
    const auto *mi_530 = buffer.data(mi + 530);
    const auto *mi_531 = buffer.data(mi + 531);
    const auto *mi_532 = buffer.data(mi + 532);
    const auto *mi_533 = buffer.data(mi + 533);
    const auto *mi_534 = buffer.data(mi + 534);
    const auto *mi_535 = buffer.data(mi + 535);
    const auto *mi_536 = buffer.data(mi + 536);
    const auto *mi_537 = buffer.data(mi + 537);
    const auto *mi_538 = buffer.data(mi + 538);
    const auto *mi_539 = buffer.data(mi + 539);
    const auto *mi_540 = buffer.data(mi + 540);
    const auto *mi_541 = buffer.data(mi + 541);
    const auto *mi_542 = buffer.data(mi + 542);
    const auto *mi_543 = buffer.data(mi + 543);
    const auto *mi_544 = buffer.data(mi + 544);
    const auto *mi_545 = buffer.data(mi + 545);
    const auto *mi_546 = buffer.data(mi + 546);
    const auto *mi_547 = buffer.data(mi + 547);
    const auto *mi_548 = buffer.data(mi + 548);
    const auto *mi_549 = buffer.data(mi + 549);
    const auto *mi_550 = buffer.data(mi + 550);
    const auto *mi_551 = buffer.data(mi + 551);
    const auto *mi_552 = buffer.data(mi + 552);
    const auto *mi_553 = buffer.data(mi + 553);
    const auto *mi_554 = buffer.data(mi + 554);
    const auto *mi_555 = buffer.data(mi + 555);
    const auto *mi_556 = buffer.data(mi + 556);
    const auto *mi_557 = buffer.data(mi + 557);
    const auto *mi_558 = buffer.data(mi + 558);
    const auto *mi_559 = buffer.data(mi + 559);
    const auto *mi_560 = buffer.data(mi + 560);
    const auto *mi_561 = buffer.data(mi + 561);
    const auto *mi_562 = buffer.data(mi + 562);
    const auto *mi_563 = buffer.data(mi + 563);
    const auto *mi_564 = buffer.data(mi + 564);
    const auto *mi_565 = buffer.data(mi + 565);
    const auto *mi_566 = buffer.data(mi + 566);
    const auto *mi_567 = buffer.data(mi + 567);
    const auto *mi_568 = buffer.data(mi + 568);
    const auto *mi_569 = buffer.data(mi + 569);
    const auto *mi_570 = buffer.data(mi + 570);
    const auto *mi_571 = buffer.data(mi + 571);
    const auto *mi_572 = buffer.data(mi + 572);
    const auto *mi_573 = buffer.data(mi + 573);
    const auto *mi_574 = buffer.data(mi + 574);
    const auto *mi_575 = buffer.data(mi + 575);
    const auto *mi_576 = buffer.data(mi + 576);
    const auto *mi_577 = buffer.data(mi + 577);
    const auto *mi_578 = buffer.data(mi + 578);
    const auto *mi_579 = buffer.data(mi + 579);
    const auto *mi_580 = buffer.data(mi + 580);
    const auto *mi_581 = buffer.data(mi + 581);
    const auto *mi_582 = buffer.data(mi + 582);
    const auto *mi_583 = buffer.data(mi + 583);
    const auto *mi_584 = buffer.data(mi + 584);
    const auto *mi_585 = buffer.data(mi + 585);
    const auto *mi_586 = buffer.data(mi + 586);
    const auto *mi_587 = buffer.data(mi + 587);
    const auto *mi_616 = buffer.data(mi + 616);
    const auto *mi_617 = buffer.data(mi + 617);
    const auto *mi_618 = buffer.data(mi + 618);
    const auto *mi_619 = buffer.data(mi + 619);
    const auto *mi_620 = buffer.data(mi + 620);
    const auto *mi_621 = buffer.data(mi + 621);
    const auto *mi_622 = buffer.data(mi + 622);
    const auto *mi_623 = buffer.data(mi + 623);
    const auto *mi_624 = buffer.data(mi + 624);
    const auto *mi_625 = buffer.data(mi + 625);
    const auto *mi_626 = buffer.data(mi + 626);
    const auto *mi_627 = buffer.data(mi + 627);
    const auto *mi_628 = buffer.data(mi + 628);
    const auto *mi_629 = buffer.data(mi + 629);
    const auto *mi_630 = buffer.data(mi + 630);
    const auto *mi_631 = buffer.data(mi + 631);
    const auto *mi_632 = buffer.data(mi + 632);
    const auto *mi_633 = buffer.data(mi + 633);
    const auto *mi_634 = buffer.data(mi + 634);
    const auto *mi_635 = buffer.data(mi + 635);
    const auto *mi_636 = buffer.data(mi + 636);
    const auto *mi_637 = buffer.data(mi + 637);
    const auto *mi_638 = buffer.data(mi + 638);
    const auto *mi_639 = buffer.data(mi + 639);
    const auto *mi_640 = buffer.data(mi + 640);
    const auto *mi_641 = buffer.data(mi + 641);
    const auto *mi_642 = buffer.data(mi + 642);
    const auto *mi_643 = buffer.data(mi + 643);
    const auto *mi_644 = buffer.data(mi + 644);
    const auto *mi_645 = buffer.data(mi + 645);
    const auto *mi_646 = buffer.data(mi + 646);
    const auto *mi_647 = buffer.data(mi + 647);
    const auto *mi_648 = buffer.data(mi + 648);
    const auto *mi_649 = buffer.data(mi + 649);
    const auto *mi_650 = buffer.data(mi + 650);
    const auto *mi_651 = buffer.data(mi + 651);
    const auto *mi_652 = buffer.data(mi + 652);
    const auto *mi_653 = buffer.data(mi + 653);
    const auto *mi_654 = buffer.data(mi + 654);
    const auto *mi_655 = buffer.data(mi + 655);
    const auto *mi_656 = buffer.data(mi + 656);
    const auto *mi_657 = buffer.data(mi + 657);
    const auto *mi_658 = buffer.data(mi + 658);
    const auto *mi_659 = buffer.data(mi + 659);
    const auto *mi_660 = buffer.data(mi + 660);
    const auto *mi_661 = buffer.data(mi + 661);
    const auto *mi_662 = buffer.data(mi + 662);
    const auto *mi_663 = buffer.data(mi + 663);
    const auto *mi_664 = buffer.data(mi + 664);
    const auto *mi_665 = buffer.data(mi + 665);
    const auto *mi_666 = buffer.data(mi + 666);
    const auto *mi_667 = buffer.data(mi + 667);
    const auto *mi_668 = buffer.data(mi + 668);
    const auto *mi_669 = buffer.data(mi + 669);
    const auto *mi_670 = buffer.data(mi + 670);
    const auto *mi_671 = buffer.data(mi + 671);
    const auto *mi_672 = buffer.data(mi + 672);
    const auto *mi_673 = buffer.data(mi + 673);
    const auto *mi_674 = buffer.data(mi + 674);
    const auto *mi_675 = buffer.data(mi + 675);
    const auto *mi_676 = buffer.data(mi + 676);
    const auto *mi_677 = buffer.data(mi + 677);
    const auto *mi_678 = buffer.data(mi + 678);
    const auto *mi_679 = buffer.data(mi + 679);
    const auto *mi_680 = buffer.data(mi + 680);
    const auto *mi_681 = buffer.data(mi + 681);
    const auto *mi_682 = buffer.data(mi + 682);
    const auto *mi_683 = buffer.data(mi + 683);
    const auto *mi_684 = buffer.data(mi + 684);
    const auto *mi_685 = buffer.data(mi + 685);
    const auto *mi_686 = buffer.data(mi + 686);
    const auto *mi_687 = buffer.data(mi + 687);
    const auto *mi_688 = buffer.data(mi + 688);
    const auto *mi_689 = buffer.data(mi + 689);
    const auto *mi_690 = buffer.data(mi + 690);
    const auto *mi_691 = buffer.data(mi + 691);
    const auto *mi_692 = buffer.data(mi + 692);
    const auto *mi_693 = buffer.data(mi + 693);
    const auto *mi_694 = buffer.data(mi + 694);
    const auto *mi_695 = buffer.data(mi + 695);
    const auto *mi_696 = buffer.data(mi + 696);
    const auto *mi_697 = buffer.data(mi + 697);
    const auto *mi_698 = buffer.data(mi + 698);
    const auto *mi_699 = buffer.data(mi + 699);
    const auto *mi_700 = buffer.data(mi + 700);
    const auto *mi_701 = buffer.data(mi + 701);
    const auto *mi_702 = buffer.data(mi + 702);

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, ki_207, ki_208, ki_209, ki_210, \
                         ki_211, mi_515, mi_516, mi_517, mi_518, \
                         mi_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -2.0 * ki_207[k]
                   + f_0 * mi_515[k];

        t_348[k] = -2.0 * ki_208[k]
                   + f_0 * mi_516[k];

        t_349[k] = -2.0 * ki_209[k]
                   + f_0 * mi_517[k];

        t_350[k] = -2.0 * ki_210[k]
                   + f_0 * mi_518[k];

        t_351[k] = -2.0 * ki_211[k]
                   + f_0 * mi_519[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, ki_212, ki_213, ki_214, ki_215, \
                         ki_216, mi_520, mi_521, mi_522, mi_523, \
                         mi_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -2.0 * ki_212[k]
                   + f_0 * mi_520[k];

        t_353[k] = -2.0 * ki_213[k]
                   + f_0 * mi_521[k];

        t_354[k] = -2.0 * ki_214[k]
                   + f_0 * mi_522[k];

        t_355[k] = -2.0 * ki_215[k]
                   + f_0 * mi_523[k];

        t_356[k] = -2.0 * ki_216[k]
                   + f_0 * mi_524[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, ki_217, ki_218, ki_219, ki_220, \
                         ki_221, mi_525, mi_526, mi_527, mi_528, \
                         mi_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -2.0 * ki_217[k]
                   + f_0 * mi_525[k];

        t_358[k] = -2.0 * ki_218[k]
                   + f_0 * mi_526[k];

        t_359[k] = -2.0 * ki_219[k]
                   + f_0 * mi_527[k];

        t_360[k] = -2.0 * ki_220[k]
                   + f_0 * mi_528[k];

        t_361[k] = -2.0 * ki_221[k]
                   + f_0 * mi_529[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, ki_222, ki_223, ki_224, ki_225, \
                         ki_226, mi_530, mi_531, mi_532, mi_533, \
                         mi_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -2.0 * ki_222[k]
                   + f_0 * mi_530[k];

        t_363[k] = -2.0 * ki_223[k]
                   + f_0 * mi_531[k];

        t_364[k] = -3.0 * ki_224[k]
                   + f_0 * mi_532[k];

        t_365[k] = -3.0 * ki_225[k]
                   + f_0 * mi_533[k];

        t_366[k] = -3.0 * ki_226[k]
                   + f_0 * mi_534[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, ki_227, ki_228, ki_229, ki_230, \
                         ki_231, mi_535, mi_536, mi_537, mi_538, \
                         mi_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -3.0 * ki_227[k]
                   + f_0 * mi_535[k];

        t_368[k] = -3.0 * ki_228[k]
                   + f_0 * mi_536[k];

        t_369[k] = -3.0 * ki_229[k]
                   + f_0 * mi_537[k];

        t_370[k] = -3.0 * ki_230[k]
                   + f_0 * mi_538[k];

        t_371[k] = -3.0 * ki_231[k]
                   + f_0 * mi_539[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, ki_232, ki_233, ki_234, ki_235, \
                         ki_236, mi_540, mi_541, mi_542, mi_543, \
                         mi_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -3.0 * ki_232[k]
                   + f_0 * mi_540[k];

        t_373[k] = -3.0 * ki_233[k]
                   + f_0 * mi_541[k];

        t_374[k] = -3.0 * ki_234[k]
                   + f_0 * mi_542[k];

        t_375[k] = -3.0 * ki_235[k]
                   + f_0 * mi_543[k];

        t_376[k] = -3.0 * ki_236[k]
                   + f_0 * mi_544[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, ki_237, ki_238, ki_239, ki_240, \
                         ki_241, mi_545, mi_546, mi_547, mi_548, \
                         mi_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -3.0 * ki_237[k]
                   + f_0 * mi_545[k];

        t_378[k] = -3.0 * ki_238[k]
                   + f_0 * mi_546[k];

        t_379[k] = -3.0 * ki_239[k]
                   + f_0 * mi_547[k];

        t_380[k] = -3.0 * ki_240[k]
                   + f_0 * mi_548[k];

        t_381[k] = -3.0 * ki_241[k]
                   + f_0 * mi_549[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, ki_242, ki_243, ki_244, ki_245, \
                         ki_246, mi_550, mi_551, mi_552, mi_553, \
                         mi_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -3.0 * ki_242[k]
                   + f_0 * mi_550[k];

        t_383[k] = -3.0 * ki_243[k]
                   + f_0 * mi_551[k];

        t_384[k] = -3.0 * ki_244[k]
                   + f_0 * mi_552[k];

        t_385[k] = -3.0 * ki_245[k]
                   + f_0 * mi_553[k];

        t_386[k] = -3.0 * ki_246[k]
                   + f_0 * mi_554[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, ki_247, ki_248, ki_249, ki_250, \
                         ki_251, mi_555, mi_556, mi_557, mi_558, \
                         mi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -3.0 * ki_247[k]
                   + f_0 * mi_555[k];

        t_388[k] = -3.0 * ki_248[k]
                   + f_0 * mi_556[k];

        t_389[k] = -3.0 * ki_249[k]
                   + f_0 * mi_557[k];

        t_390[k] = -3.0 * ki_250[k]
                   + f_0 * mi_558[k];

        t_391[k] = -3.0 * ki_251[k]
                   + f_0 * mi_559[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, ki_252, ki_253, ki_254, ki_255, \
                         ki_256, mi_560, mi_561, mi_562, mi_563, \
                         mi_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -4.0 * ki_252[k]
                   + f_0 * mi_560[k];

        t_393[k] = -4.0 * ki_253[k]
                   + f_0 * mi_561[k];

        t_394[k] = -4.0 * ki_254[k]
                   + f_0 * mi_562[k];

        t_395[k] = -4.0 * ki_255[k]
                   + f_0 * mi_563[k];

        t_396[k] = -4.0 * ki_256[k]
                   + f_0 * mi_564[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, ki_257, ki_258, ki_259, ki_260, \
                         ki_261, mi_565, mi_566, mi_567, mi_568, \
                         mi_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -4.0 * ki_257[k]
                   + f_0 * mi_565[k];

        t_398[k] = -4.0 * ki_258[k]
                   + f_0 * mi_566[k];

        t_399[k] = -4.0 * ki_259[k]
                   + f_0 * mi_567[k];

        t_400[k] = -4.0 * ki_260[k]
                   + f_0 * mi_568[k];

        t_401[k] = -4.0 * ki_261[k]
                   + f_0 * mi_569[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, ki_262, ki_263, ki_264, ki_265, \
                         ki_266, mi_570, mi_571, mi_572, mi_573, \
                         mi_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -4.0 * ki_262[k]
                   + f_0 * mi_570[k];

        t_403[k] = -4.0 * ki_263[k]
                   + f_0 * mi_571[k];

        t_404[k] = -4.0 * ki_264[k]
                   + f_0 * mi_572[k];

        t_405[k] = -4.0 * ki_265[k]
                   + f_0 * mi_573[k];

        t_406[k] = -4.0 * ki_266[k]
                   + f_0 * mi_574[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, ki_267, ki_268, ki_269, ki_270, \
                         ki_271, mi_575, mi_576, mi_577, mi_578, \
                         mi_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -4.0 * ki_267[k]
                   + f_0 * mi_575[k];

        t_408[k] = -4.0 * ki_268[k]
                   + f_0 * mi_576[k];

        t_409[k] = -4.0 * ki_269[k]
                   + f_0 * mi_577[k];

        t_410[k] = -4.0 * ki_270[k]
                   + f_0 * mi_578[k];

        t_411[k] = -4.0 * ki_271[k]
                   + f_0 * mi_579[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, ki_272, ki_273, ki_274, ki_275, \
                         ki_276, mi_580, mi_581, mi_582, mi_583, \
                         mi_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -4.0 * ki_272[k]
                   + f_0 * mi_580[k];

        t_413[k] = -4.0 * ki_273[k]
                   + f_0 * mi_581[k];

        t_414[k] = -4.0 * ki_274[k]
                   + f_0 * mi_582[k];

        t_415[k] = -4.0 * ki_275[k]
                   + f_0 * mi_583[k];

        t_416[k] = -4.0 * ki_276[k]
                   + f_0 * mi_584[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, t_422, ki_277, ki_278, ki_279, \
                         mi_585, mi_586, mi_587, mi_616, mi_617, \
                         mi_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -4.0 * ki_277[k]
                   + f_0 * mi_585[k];

        t_418[k] = -4.0 * ki_278[k]
                   + f_0 * mi_586[k];

        t_419[k] = -4.0 * ki_279[k]
                   + f_0 * mi_587[k];

        t_420[k] = f_0 * mi_616[k];

        t_421[k] = f_0 * mi_617[k];

        t_422[k] = f_0 * mi_618[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, t_428, t_429, t_430, mi_619, \
                         mi_620, mi_621, mi_622, mi_623, mi_624, mi_625, \
                         mi_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_0 * mi_619[k];

        t_424[k] = f_0 * mi_620[k];

        t_425[k] = f_0 * mi_621[k];

        t_426[k] = f_0 * mi_622[k];

        t_427[k] = f_0 * mi_623[k];

        t_428[k] = f_0 * mi_624[k];

        t_429[k] = f_0 * mi_625[k];

        t_430[k] = f_0 * mi_626[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, t_436, t_437, t_438, mi_627, \
                         mi_628, mi_629, mi_630, mi_631, mi_632, mi_633, \
                         mi_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_0 * mi_627[k];

        t_432[k] = f_0 * mi_628[k];

        t_433[k] = f_0 * mi_629[k];

        t_434[k] = f_0 * mi_630[k];

        t_435[k] = f_0 * mi_631[k];

        t_436[k] = f_0 * mi_632[k];

        t_437[k] = f_0 * mi_633[k];

        t_438[k] = f_0 * mi_634[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, t_443, t_444, t_445, t_446, mi_635, \
                         mi_636, mi_637, mi_638, mi_639, mi_640, mi_641, \
                         mi_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_0 * mi_635[k];

        t_440[k] = f_0 * mi_636[k];

        t_441[k] = f_0 * mi_637[k];

        t_442[k] = f_0 * mi_638[k];

        t_443[k] = f_0 * mi_639[k];

        t_444[k] = f_0 * mi_640[k];

        t_445[k] = f_0 * mi_641[k];

        t_446[k] = f_0 * mi_642[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, t_451, ki_280, ki_281, ki_282, ki_283, \
                         mi_643, mi_644, mi_645, mi_646, mi_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_0 * mi_643[k];

        t_448[k] = -ki_280[k]
                   + f_0 * mi_644[k];

        t_449[k] = -ki_281[k]
                   + f_0 * mi_645[k];

        t_450[k] = -ki_282[k]
                   + f_0 * mi_646[k];

        t_451[k] = -ki_283[k]
                   + f_0 * mi_647[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, ki_284, ki_285, ki_286, ki_287, \
                         ki_288, mi_648, mi_649, mi_650, mi_651, \
                         mi_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = -ki_284[k]
                   + f_0 * mi_648[k];

        t_453[k] = -ki_285[k]
                   + f_0 * mi_649[k];

        t_454[k] = -ki_286[k]
                   + f_0 * mi_650[k];

        t_455[k] = -ki_287[k]
                   + f_0 * mi_651[k];

        t_456[k] = -ki_288[k]
                   + f_0 * mi_652[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, t_461, ki_289, ki_290, ki_291, ki_292, \
                         ki_293, mi_653, mi_654, mi_655, mi_656, \
                         mi_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = -ki_289[k]
                   + f_0 * mi_653[k];

        t_458[k] = -ki_290[k]
                   + f_0 * mi_654[k];

        t_459[k] = -ki_291[k]
                   + f_0 * mi_655[k];

        t_460[k] = -ki_292[k]
                   + f_0 * mi_656[k];

        t_461[k] = -ki_293[k]
                   + f_0 * mi_657[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, ki_294, ki_295, ki_296, ki_297, \
                         ki_298, mi_658, mi_659, mi_660, mi_661, \
                         mi_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = -ki_294[k]
                   + f_0 * mi_658[k];

        t_463[k] = -ki_295[k]
                   + f_0 * mi_659[k];

        t_464[k] = -ki_296[k]
                   + f_0 * mi_660[k];

        t_465[k] = -ki_297[k]
                   + f_0 * mi_661[k];

        t_466[k] = -ki_298[k]
                   + f_0 * mi_662[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, ki_299, ki_300, ki_301, ki_302, \
                         ki_303, mi_663, mi_664, mi_665, mi_666, \
                         mi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = -ki_299[k]
                   + f_0 * mi_663[k];

        t_468[k] = -ki_300[k]
                   + f_0 * mi_664[k];

        t_469[k] = -ki_301[k]
                   + f_0 * mi_665[k];

        t_470[k] = -ki_302[k]
                   + f_0 * mi_666[k];

        t_471[k] = -ki_303[k]
                   + f_0 * mi_667[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, ki_304, ki_305, ki_306, ki_307, \
                         ki_308, mi_668, mi_669, mi_670, mi_671, \
                         mi_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = -ki_304[k]
                   + f_0 * mi_668[k];

        t_473[k] = -ki_305[k]
                   + f_0 * mi_669[k];

        t_474[k] = -ki_306[k]
                   + f_0 * mi_670[k];

        t_475[k] = -ki_307[k]
                   + f_0 * mi_671[k];

        t_476[k] = -2.0 * ki_308[k]
                   + f_0 * mi_672[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, ki_309, ki_310, ki_311, ki_312, \
                         ki_313, mi_673, mi_674, mi_675, mi_676, \
                         mi_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = -2.0 * ki_309[k]
                   + f_0 * mi_673[k];

        t_478[k] = -2.0 * ki_310[k]
                   + f_0 * mi_674[k];

        t_479[k] = -2.0 * ki_311[k]
                   + f_0 * mi_675[k];

        t_480[k] = -2.0 * ki_312[k]
                   + f_0 * mi_676[k];

        t_481[k] = -2.0 * ki_313[k]
                   + f_0 * mi_677[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, t_486, ki_314, ki_315, ki_316, ki_317, \
                         ki_318, mi_678, mi_679, mi_680, mi_681, \
                         mi_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = -2.0 * ki_314[k]
                   + f_0 * mi_678[k];

        t_483[k] = -2.0 * ki_315[k]
                   + f_0 * mi_679[k];

        t_484[k] = -2.0 * ki_316[k]
                   + f_0 * mi_680[k];

        t_485[k] = -2.0 * ki_317[k]
                   + f_0 * mi_681[k];

        t_486[k] = -2.0 * ki_318[k]
                   + f_0 * mi_682[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, ki_319, ki_320, ki_321, ki_322, \
                         ki_323, mi_683, mi_684, mi_685, mi_686, \
                         mi_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = -2.0 * ki_319[k]
                   + f_0 * mi_683[k];

        t_488[k] = -2.0 * ki_320[k]
                   + f_0 * mi_684[k];

        t_489[k] = -2.0 * ki_321[k]
                   + f_0 * mi_685[k];

        t_490[k] = -2.0 * ki_322[k]
                   + f_0 * mi_686[k];

        t_491[k] = -2.0 * ki_323[k]
                   + f_0 * mi_687[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, ki_324, ki_325, ki_326, ki_327, \
                         ki_328, mi_688, mi_689, mi_690, mi_691, \
                         mi_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = -2.0 * ki_324[k]
                   + f_0 * mi_688[k];

        t_493[k] = -2.0 * ki_325[k]
                   + f_0 * mi_689[k];

        t_494[k] = -2.0 * ki_326[k]
                   + f_0 * mi_690[k];

        t_495[k] = -2.0 * ki_327[k]
                   + f_0 * mi_691[k];

        t_496[k] = -2.0 * ki_328[k]
                   + f_0 * mi_692[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, t_500, t_501, ki_329, ki_330, ki_331, ki_332, \
                         ki_333, mi_693, mi_694, mi_695, mi_696, \
                         mi_697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = -2.0 * ki_329[k]
                   + f_0 * mi_693[k];

        t_498[k] = -2.0 * ki_330[k]
                   + f_0 * mi_694[k];

        t_499[k] = -2.0 * ki_331[k]
                   + f_0 * mi_695[k];

        t_500[k] = -2.0 * ki_332[k]
                   + f_0 * mi_696[k];

        t_501[k] = -2.0 * ki_333[k]
                   + f_0 * mi_697[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, ki_334, ki_335, ki_336, ki_337, \
                         ki_338, mi_698, mi_699, mi_700, mi_701, \
                         mi_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = -2.0 * ki_334[k]
                   + f_0 * mi_698[k];

        t_503[k] = -2.0 * ki_335[k]
                   + f_0 * mi_699[k];

        t_504[k] = -3.0 * ki_336[k]
                   + f_0 * mi_700[k];

        t_505[k] = -3.0 * ki_337[k]
                   + f_0 * mi_701[k];

        t_506[k] = -3.0 * ki_338[k]
                   + f_0 * mi_702[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ki_339 = buffer.data(ki + 339);
    const auto *ki_340 = buffer.data(ki + 340);
    const auto *ki_341 = buffer.data(ki + 341);
    const auto *ki_342 = buffer.data(ki + 342);
    const auto *ki_343 = buffer.data(ki + 343);
    const auto *ki_344 = buffer.data(ki + 344);
    const auto *ki_345 = buffer.data(ki + 345);
    const auto *ki_346 = buffer.data(ki + 346);
    const auto *ki_347 = buffer.data(ki + 347);
    const auto *ki_348 = buffer.data(ki + 348);
    const auto *ki_349 = buffer.data(ki + 349);
    const auto *ki_350 = buffer.data(ki + 350);
    const auto *ki_351 = buffer.data(ki + 351);
    const auto *ki_352 = buffer.data(ki + 352);
    const auto *ki_353 = buffer.data(ki + 353);
    const auto *ki_354 = buffer.data(ki + 354);
    const auto *ki_355 = buffer.data(ki + 355);
    const auto *ki_356 = buffer.data(ki + 356);
    const auto *ki_357 = buffer.data(ki + 357);
    const auto *ki_358 = buffer.data(ki + 358);
    const auto *ki_359 = buffer.data(ki + 359);
    const auto *ki_360 = buffer.data(ki + 360);
    const auto *ki_361 = buffer.data(ki + 361);
    const auto *ki_362 = buffer.data(ki + 362);
    const auto *ki_363 = buffer.data(ki + 363);
    const auto *ki_364 = buffer.data(ki + 364);
    const auto *ki_365 = buffer.data(ki + 365);
    const auto *ki_366 = buffer.data(ki + 366);
    const auto *ki_367 = buffer.data(ki + 367);
    const auto *ki_368 = buffer.data(ki + 368);
    const auto *ki_369 = buffer.data(ki + 369);
    const auto *ki_370 = buffer.data(ki + 370);
    const auto *ki_371 = buffer.data(ki + 371);
    const auto *ki_372 = buffer.data(ki + 372);
    const auto *ki_373 = buffer.data(ki + 373);
    const auto *ki_374 = buffer.data(ki + 374);
    const auto *ki_375 = buffer.data(ki + 375);
    const auto *ki_376 = buffer.data(ki + 376);
    const auto *ki_377 = buffer.data(ki + 377);
    const auto *ki_378 = buffer.data(ki + 378);
    const auto *ki_379 = buffer.data(ki + 379);
    const auto *ki_380 = buffer.data(ki + 380);
    const auto *ki_381 = buffer.data(ki + 381);
    const auto *ki_382 = buffer.data(ki + 382);
    const auto *ki_383 = buffer.data(ki + 383);
    const auto *ki_384 = buffer.data(ki + 384);
    const auto *ki_385 = buffer.data(ki + 385);
    const auto *ki_386 = buffer.data(ki + 386);
    const auto *ki_387 = buffer.data(ki + 387);
    const auto *ki_388 = buffer.data(ki + 388);
    const auto *ki_389 = buffer.data(ki + 389);
    const auto *ki_390 = buffer.data(ki + 390);
    const auto *ki_391 = buffer.data(ki + 391);
    const auto *ki_392 = buffer.data(ki + 392);
    const auto *ki_393 = buffer.data(ki + 393);
    const auto *ki_394 = buffer.data(ki + 394);
    const auto *ki_395 = buffer.data(ki + 395);
    const auto *ki_396 = buffer.data(ki + 396);
    const auto *ki_397 = buffer.data(ki + 397);
    const auto *ki_398 = buffer.data(ki + 398);
    const auto *ki_399 = buffer.data(ki + 399);
    const auto *ki_400 = buffer.data(ki + 400);
    const auto *ki_401 = buffer.data(ki + 401);
    const auto *ki_402 = buffer.data(ki + 402);
    const auto *ki_403 = buffer.data(ki + 403);
    const auto *ki_404 = buffer.data(ki + 404);
    const auto *ki_405 = buffer.data(ki + 405);
    const auto *ki_406 = buffer.data(ki + 406);
    const auto *ki_407 = buffer.data(ki + 407);
    const auto *ki_408 = buffer.data(ki + 408);
    const auto *ki_409 = buffer.data(ki + 409);
    const auto *ki_410 = buffer.data(ki + 410);
    const auto *ki_411 = buffer.data(ki + 411);
    const auto *ki_412 = buffer.data(ki + 412);
    const auto *ki_413 = buffer.data(ki + 413);
    const auto *ki_414 = buffer.data(ki + 414);
    const auto *ki_415 = buffer.data(ki + 415);
    const auto *ki_416 = buffer.data(ki + 416);
    const auto *ki_417 = buffer.data(ki + 417);
    const auto *ki_418 = buffer.data(ki + 418);
    const auto *ki_419 = buffer.data(ki + 419);
    const auto *ki_420 = buffer.data(ki + 420);
    const auto *ki_421 = buffer.data(ki + 421);
    const auto *ki_422 = buffer.data(ki + 422);
    const auto *ki_423 = buffer.data(ki + 423);
    const auto *ki_424 = buffer.data(ki + 424);
    const auto *ki_425 = buffer.data(ki + 425);
    const auto *ki_426 = buffer.data(ki + 426);
    const auto *ki_427 = buffer.data(ki + 427);
    const auto *ki_428 = buffer.data(ki + 428);
    const auto *ki_429 = buffer.data(ki + 429);
    const auto *ki_430 = buffer.data(ki + 430);
    const auto *ki_431 = buffer.data(ki + 431);
    const auto *ki_432 = buffer.data(ki + 432);
    const auto *ki_433 = buffer.data(ki + 433);
    const auto *ki_434 = buffer.data(ki + 434);
    const auto *ki_435 = buffer.data(ki + 435);
    const auto *ki_436 = buffer.data(ki + 436);
    const auto *ki_437 = buffer.data(ki + 437);
    const auto *ki_438 = buffer.data(ki + 438);
    const auto *ki_439 = buffer.data(ki + 439);
    const auto *ki_440 = buffer.data(ki + 440);
    const auto *ki_441 = buffer.data(ki + 441);
    const auto *ki_442 = buffer.data(ki + 442);
    const auto *ki_443 = buffer.data(ki + 443);
    const auto *ki_444 = buffer.data(ki + 444);
    const auto *ki_445 = buffer.data(ki + 445);
    const auto *ki_446 = buffer.data(ki + 446);
    const auto *ki_447 = buffer.data(ki + 447);
    const auto *ki_448 = buffer.data(ki + 448);
    const auto *ki_449 = buffer.data(ki + 449);
    const auto *ki_450 = buffer.data(ki + 450);
    const auto *ki_451 = buffer.data(ki + 451);
    const auto *ki_452 = buffer.data(ki + 452);
    const auto *ki_453 = buffer.data(ki + 453);
    const auto *ki_454 = buffer.data(ki + 454);
    const auto *ki_455 = buffer.data(ki + 455);
    const auto *ki_456 = buffer.data(ki + 456);
    const auto *ki_457 = buffer.data(ki + 457);
    const auto *ki_458 = buffer.data(ki + 458);
    const auto *ki_459 = buffer.data(ki + 459);
    const auto *ki_460 = buffer.data(ki + 460);
    const auto *ki_461 = buffer.data(ki + 461);
    const auto *ki_462 = buffer.data(ki + 462);
    const auto *ki_463 = buffer.data(ki + 463);
    const auto *ki_464 = buffer.data(ki + 464);
    const auto *ki_465 = buffer.data(ki + 465);
    const auto *ki_466 = buffer.data(ki + 466);
    const auto *ki_467 = buffer.data(ki + 467);
    const auto *ki_468 = buffer.data(ki + 468);
    const auto *ki_469 = buffer.data(ki + 469);
    const auto *ki_470 = buffer.data(ki + 470);

    const auto *mi_703 = buffer.data(mi + 703);
    const auto *mi_704 = buffer.data(mi + 704);
    const auto *mi_705 = buffer.data(mi + 705);
    const auto *mi_706 = buffer.data(mi + 706);
    const auto *mi_707 = buffer.data(mi + 707);
    const auto *mi_708 = buffer.data(mi + 708);
    const auto *mi_709 = buffer.data(mi + 709);
    const auto *mi_710 = buffer.data(mi + 710);
    const auto *mi_711 = buffer.data(mi + 711);
    const auto *mi_712 = buffer.data(mi + 712);
    const auto *mi_713 = buffer.data(mi + 713);
    const auto *mi_714 = buffer.data(mi + 714);
    const auto *mi_715 = buffer.data(mi + 715);
    const auto *mi_716 = buffer.data(mi + 716);
    const auto *mi_717 = buffer.data(mi + 717);
    const auto *mi_718 = buffer.data(mi + 718);
    const auto *mi_719 = buffer.data(mi + 719);
    const auto *mi_720 = buffer.data(mi + 720);
    const auto *mi_721 = buffer.data(mi + 721);
    const auto *mi_722 = buffer.data(mi + 722);
    const auto *mi_723 = buffer.data(mi + 723);
    const auto *mi_724 = buffer.data(mi + 724);
    const auto *mi_725 = buffer.data(mi + 725);
    const auto *mi_726 = buffer.data(mi + 726);
    const auto *mi_727 = buffer.data(mi + 727);
    const auto *mi_728 = buffer.data(mi + 728);
    const auto *mi_729 = buffer.data(mi + 729);
    const auto *mi_730 = buffer.data(mi + 730);
    const auto *mi_731 = buffer.data(mi + 731);
    const auto *mi_732 = buffer.data(mi + 732);
    const auto *mi_733 = buffer.data(mi + 733);
    const auto *mi_734 = buffer.data(mi + 734);
    const auto *mi_735 = buffer.data(mi + 735);
    const auto *mi_736 = buffer.data(mi + 736);
    const auto *mi_737 = buffer.data(mi + 737);
    const auto *mi_738 = buffer.data(mi + 738);
    const auto *mi_739 = buffer.data(mi + 739);
    const auto *mi_740 = buffer.data(mi + 740);
    const auto *mi_741 = buffer.data(mi + 741);
    const auto *mi_742 = buffer.data(mi + 742);
    const auto *mi_743 = buffer.data(mi + 743);
    const auto *mi_744 = buffer.data(mi + 744);
    const auto *mi_745 = buffer.data(mi + 745);
    const auto *mi_746 = buffer.data(mi + 746);
    const auto *mi_747 = buffer.data(mi + 747);
    const auto *mi_748 = buffer.data(mi + 748);
    const auto *mi_749 = buffer.data(mi + 749);
    const auto *mi_750 = buffer.data(mi + 750);
    const auto *mi_751 = buffer.data(mi + 751);
    const auto *mi_752 = buffer.data(mi + 752);
    const auto *mi_753 = buffer.data(mi + 753);
    const auto *mi_754 = buffer.data(mi + 754);
    const auto *mi_755 = buffer.data(mi + 755);
    const auto *mi_756 = buffer.data(mi + 756);
    const auto *mi_757 = buffer.data(mi + 757);
    const auto *mi_758 = buffer.data(mi + 758);
    const auto *mi_759 = buffer.data(mi + 759);
    const auto *mi_760 = buffer.data(mi + 760);
    const auto *mi_761 = buffer.data(mi + 761);
    const auto *mi_762 = buffer.data(mi + 762);
    const auto *mi_763 = buffer.data(mi + 763);
    const auto *mi_764 = buffer.data(mi + 764);
    const auto *mi_765 = buffer.data(mi + 765);
    const auto *mi_766 = buffer.data(mi + 766);
    const auto *mi_767 = buffer.data(mi + 767);
    const auto *mi_768 = buffer.data(mi + 768);
    const auto *mi_769 = buffer.data(mi + 769);
    const auto *mi_770 = buffer.data(mi + 770);
    const auto *mi_771 = buffer.data(mi + 771);
    const auto *mi_772 = buffer.data(mi + 772);
    const auto *mi_773 = buffer.data(mi + 773);
    const auto *mi_774 = buffer.data(mi + 774);
    const auto *mi_775 = buffer.data(mi + 775);
    const auto *mi_776 = buffer.data(mi + 776);
    const auto *mi_777 = buffer.data(mi + 777);
    const auto *mi_778 = buffer.data(mi + 778);
    const auto *mi_779 = buffer.data(mi + 779);
    const auto *mi_780 = buffer.data(mi + 780);
    const auto *mi_781 = buffer.data(mi + 781);
    const auto *mi_782 = buffer.data(mi + 782);
    const auto *mi_783 = buffer.data(mi + 783);
    const auto *mi_812 = buffer.data(mi + 812);
    const auto *mi_813 = buffer.data(mi + 813);
    const auto *mi_814 = buffer.data(mi + 814);
    const auto *mi_815 = buffer.data(mi + 815);
    const auto *mi_816 = buffer.data(mi + 816);
    const auto *mi_817 = buffer.data(mi + 817);
    const auto *mi_818 = buffer.data(mi + 818);
    const auto *mi_819 = buffer.data(mi + 819);
    const auto *mi_820 = buffer.data(mi + 820);
    const auto *mi_821 = buffer.data(mi + 821);
    const auto *mi_822 = buffer.data(mi + 822);
    const auto *mi_823 = buffer.data(mi + 823);
    const auto *mi_824 = buffer.data(mi + 824);
    const auto *mi_825 = buffer.data(mi + 825);
    const auto *mi_826 = buffer.data(mi + 826);
    const auto *mi_827 = buffer.data(mi + 827);
    const auto *mi_828 = buffer.data(mi + 828);
    const auto *mi_829 = buffer.data(mi + 829);
    const auto *mi_830 = buffer.data(mi + 830);
    const auto *mi_831 = buffer.data(mi + 831);
    const auto *mi_832 = buffer.data(mi + 832);
    const auto *mi_833 = buffer.data(mi + 833);
    const auto *mi_834 = buffer.data(mi + 834);
    const auto *mi_835 = buffer.data(mi + 835);
    const auto *mi_836 = buffer.data(mi + 836);
    const auto *mi_837 = buffer.data(mi + 837);
    const auto *mi_838 = buffer.data(mi + 838);
    const auto *mi_839 = buffer.data(mi + 839);
    const auto *mi_840 = buffer.data(mi + 840);
    const auto *mi_841 = buffer.data(mi + 841);
    const auto *mi_842 = buffer.data(mi + 842);
    const auto *mi_843 = buffer.data(mi + 843);
    const auto *mi_844 = buffer.data(mi + 844);
    const auto *mi_845 = buffer.data(mi + 845);
    const auto *mi_846 = buffer.data(mi + 846);
    const auto *mi_847 = buffer.data(mi + 847);
    const auto *mi_848 = buffer.data(mi + 848);
    const auto *mi_849 = buffer.data(mi + 849);
    const auto *mi_850 = buffer.data(mi + 850);
    const auto *mi_851 = buffer.data(mi + 851);
    const auto *mi_852 = buffer.data(mi + 852);
    const auto *mi_853 = buffer.data(mi + 853);
    const auto *mi_854 = buffer.data(mi + 854);
    const auto *mi_855 = buffer.data(mi + 855);
    const auto *mi_856 = buffer.data(mi + 856);
    const auto *mi_857 = buffer.data(mi + 857);
    const auto *mi_858 = buffer.data(mi + 858);
    const auto *mi_859 = buffer.data(mi + 859);
    const auto *mi_860 = buffer.data(mi + 860);
    const auto *mi_861 = buffer.data(mi + 861);
    const auto *mi_862 = buffer.data(mi + 862);
    const auto *mi_863 = buffer.data(mi + 863);
    const auto *mi_864 = buffer.data(mi + 864);
    const auto *mi_865 = buffer.data(mi + 865);
    const auto *mi_866 = buffer.data(mi + 866);
    const auto *mi_867 = buffer.data(mi + 867);
    const auto *mi_868 = buffer.data(mi + 868);
    const auto *mi_869 = buffer.data(mi + 869);
    const auto *mi_870 = buffer.data(mi + 870);
    const auto *mi_871 = buffer.data(mi + 871);
    const auto *mi_872 = buffer.data(mi + 872);
    const auto *mi_873 = buffer.data(mi + 873);
    const auto *mi_874 = buffer.data(mi + 874);
    const auto *mi_875 = buffer.data(mi + 875);
    const auto *mi_876 = buffer.data(mi + 876);
    const auto *mi_877 = buffer.data(mi + 877);
    const auto *mi_878 = buffer.data(mi + 878);
    const auto *mi_879 = buffer.data(mi + 879);
    const auto *mi_880 = buffer.data(mi + 880);
    const auto *mi_881 = buffer.data(mi + 881);
    const auto *mi_882 = buffer.data(mi + 882);
    const auto *mi_883 = buffer.data(mi + 883);
    const auto *mi_884 = buffer.data(mi + 884);
    const auto *mi_885 = buffer.data(mi + 885);
    const auto *mi_886 = buffer.data(mi + 886);
    const auto *mi_887 = buffer.data(mi + 887);
    const auto *mi_888 = buffer.data(mi + 888);
    const auto *mi_889 = buffer.data(mi + 889);
    const auto *mi_890 = buffer.data(mi + 890);

#pragma omp simd aligned(t_507, t_508, t_509, t_510, t_511, ki_339, ki_340, ki_341, ki_342, \
                         ki_343, mi_703, mi_704, mi_705, mi_706, \
                         mi_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = -3.0 * ki_339[k]
                   + f_0 * mi_703[k];

        t_508[k] = -3.0 * ki_340[k]
                   + f_0 * mi_704[k];

        t_509[k] = -3.0 * ki_341[k]
                   + f_0 * mi_705[k];

        t_510[k] = -3.0 * ki_342[k]
                   + f_0 * mi_706[k];

        t_511[k] = -3.0 * ki_343[k]
                   + f_0 * mi_707[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, ki_344, ki_345, ki_346, ki_347, \
                         ki_348, mi_708, mi_709, mi_710, mi_711, \
                         mi_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = -3.0 * ki_344[k]
                   + f_0 * mi_708[k];

        t_513[k] = -3.0 * ki_345[k]
                   + f_0 * mi_709[k];

        t_514[k] = -3.0 * ki_346[k]
                   + f_0 * mi_710[k];

        t_515[k] = -3.0 * ki_347[k]
                   + f_0 * mi_711[k];

        t_516[k] = -3.0 * ki_348[k]
                   + f_0 * mi_712[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, ki_349, ki_350, ki_351, ki_352, \
                         ki_353, mi_713, mi_714, mi_715, mi_716, \
                         mi_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = -3.0 * ki_349[k]
                   + f_0 * mi_713[k];

        t_518[k] = -3.0 * ki_350[k]
                   + f_0 * mi_714[k];

        t_519[k] = -3.0 * ki_351[k]
                   + f_0 * mi_715[k];

        t_520[k] = -3.0 * ki_352[k]
                   + f_0 * mi_716[k];

        t_521[k] = -3.0 * ki_353[k]
                   + f_0 * mi_717[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, ki_354, ki_355, ki_356, ki_357, \
                         ki_358, mi_718, mi_719, mi_720, mi_721, \
                         mi_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = -3.0 * ki_354[k]
                   + f_0 * mi_718[k];

        t_523[k] = -3.0 * ki_355[k]
                   + f_0 * mi_719[k];

        t_524[k] = -3.0 * ki_356[k]
                   + f_0 * mi_720[k];

        t_525[k] = -3.0 * ki_357[k]
                   + f_0 * mi_721[k];

        t_526[k] = -3.0 * ki_358[k]
                   + f_0 * mi_722[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, ki_359, ki_360, ki_361, ki_362, \
                         ki_363, mi_723, mi_724, mi_725, mi_726, \
                         mi_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = -3.0 * ki_359[k]
                   + f_0 * mi_723[k];

        t_528[k] = -3.0 * ki_360[k]
                   + f_0 * mi_724[k];

        t_529[k] = -3.0 * ki_361[k]
                   + f_0 * mi_725[k];

        t_530[k] = -3.0 * ki_362[k]
                   + f_0 * mi_726[k];

        t_531[k] = -3.0 * ki_363[k]
                   + f_0 * mi_727[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, ki_364, ki_365, ki_366, ki_367, \
                         ki_368, mi_728, mi_729, mi_730, mi_731, \
                         mi_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = -4.0 * ki_364[k]
                   + f_0 * mi_728[k];

        t_533[k] = -4.0 * ki_365[k]
                   + f_0 * mi_729[k];

        t_534[k] = -4.0 * ki_366[k]
                   + f_0 * mi_730[k];

        t_535[k] = -4.0 * ki_367[k]
                   + f_0 * mi_731[k];

        t_536[k] = -4.0 * ki_368[k]
                   + f_0 * mi_732[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, t_541, ki_369, ki_370, ki_371, ki_372, \
                         ki_373, mi_733, mi_734, mi_735, mi_736, \
                         mi_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = -4.0 * ki_369[k]
                   + f_0 * mi_733[k];

        t_538[k] = -4.0 * ki_370[k]
                   + f_0 * mi_734[k];

        t_539[k] = -4.0 * ki_371[k]
                   + f_0 * mi_735[k];

        t_540[k] = -4.0 * ki_372[k]
                   + f_0 * mi_736[k];

        t_541[k] = -4.0 * ki_373[k]
                   + f_0 * mi_737[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, ki_374, ki_375, ki_376, ki_377, \
                         ki_378, mi_738, mi_739, mi_740, mi_741, \
                         mi_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = -4.0 * ki_374[k]
                   + f_0 * mi_738[k];

        t_543[k] = -4.0 * ki_375[k]
                   + f_0 * mi_739[k];

        t_544[k] = -4.0 * ki_376[k]
                   + f_0 * mi_740[k];

        t_545[k] = -4.0 * ki_377[k]
                   + f_0 * mi_741[k];

        t_546[k] = -4.0 * ki_378[k]
                   + f_0 * mi_742[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, ki_379, ki_380, ki_381, ki_382, \
                         ki_383, mi_743, mi_744, mi_745, mi_746, \
                         mi_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = -4.0 * ki_379[k]
                   + f_0 * mi_743[k];

        t_548[k] = -4.0 * ki_380[k]
                   + f_0 * mi_744[k];

        t_549[k] = -4.0 * ki_381[k]
                   + f_0 * mi_745[k];

        t_550[k] = -4.0 * ki_382[k]
                   + f_0 * mi_746[k];

        t_551[k] = -4.0 * ki_383[k]
                   + f_0 * mi_747[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, t_556, ki_384, ki_385, ki_386, ki_387, \
                         ki_388, mi_748, mi_749, mi_750, mi_751, \
                         mi_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = -4.0 * ki_384[k]
                   + f_0 * mi_748[k];

        t_553[k] = -4.0 * ki_385[k]
                   + f_0 * mi_749[k];

        t_554[k] = -4.0 * ki_386[k]
                   + f_0 * mi_750[k];

        t_555[k] = -4.0 * ki_387[k]
                   + f_0 * mi_751[k];

        t_556[k] = -4.0 * ki_388[k]
                   + f_0 * mi_752[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, ki_389, ki_390, ki_391, ki_392, \
                         ki_393, mi_753, mi_754, mi_755, mi_756, \
                         mi_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = -4.0 * ki_389[k]
                   + f_0 * mi_753[k];

        t_558[k] = -4.0 * ki_390[k]
                   + f_0 * mi_754[k];

        t_559[k] = -4.0 * ki_391[k]
                   + f_0 * mi_755[k];

        t_560[k] = -5.0 * ki_392[k]
                   + f_0 * mi_756[k];

        t_561[k] = -5.0 * ki_393[k]
                   + f_0 * mi_757[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, ki_394, ki_395, ki_396, ki_397, \
                         ki_398, mi_758, mi_759, mi_760, mi_761, \
                         mi_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = -5.0 * ki_394[k]
                   + f_0 * mi_758[k];

        t_563[k] = -5.0 * ki_395[k]
                   + f_0 * mi_759[k];

        t_564[k] = -5.0 * ki_396[k]
                   + f_0 * mi_760[k];

        t_565[k] = -5.0 * ki_397[k]
                   + f_0 * mi_761[k];

        t_566[k] = -5.0 * ki_398[k]
                   + f_0 * mi_762[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, ki_399, ki_400, ki_401, ki_402, \
                         ki_403, mi_763, mi_764, mi_765, mi_766, \
                         mi_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -5.0 * ki_399[k]
                   + f_0 * mi_763[k];

        t_568[k] = -5.0 * ki_400[k]
                   + f_0 * mi_764[k];

        t_569[k] = -5.0 * ki_401[k]
                   + f_0 * mi_765[k];

        t_570[k] = -5.0 * ki_402[k]
                   + f_0 * mi_766[k];

        t_571[k] = -5.0 * ki_403[k]
                   + f_0 * mi_767[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, ki_404, ki_405, ki_406, ki_407, \
                         ki_408, mi_768, mi_769, mi_770, mi_771, \
                         mi_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = -5.0 * ki_404[k]
                   + f_0 * mi_768[k];

        t_573[k] = -5.0 * ki_405[k]
                   + f_0 * mi_769[k];

        t_574[k] = -5.0 * ki_406[k]
                   + f_0 * mi_770[k];

        t_575[k] = -5.0 * ki_407[k]
                   + f_0 * mi_771[k];

        t_576[k] = -5.0 * ki_408[k]
                   + f_0 * mi_772[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, ki_409, ki_410, ki_411, ki_412, \
                         ki_413, mi_773, mi_774, mi_775, mi_776, \
                         mi_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = -5.0 * ki_409[k]
                   + f_0 * mi_773[k];

        t_578[k] = -5.0 * ki_410[k]
                   + f_0 * mi_774[k];

        t_579[k] = -5.0 * ki_411[k]
                   + f_0 * mi_775[k];

        t_580[k] = -5.0 * ki_412[k]
                   + f_0 * mi_776[k];

        t_581[k] = -5.0 * ki_413[k]
                   + f_0 * mi_777[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, ki_414, ki_415, ki_416, ki_417, \
                         ki_418, mi_778, mi_779, mi_780, mi_781, \
                         mi_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = -5.0 * ki_414[k]
                   + f_0 * mi_778[k];

        t_583[k] = -5.0 * ki_415[k]
                   + f_0 * mi_779[k];

        t_584[k] = -5.0 * ki_416[k]
                   + f_0 * mi_780[k];

        t_585[k] = -5.0 * ki_417[k]
                   + f_0 * mi_781[k];

        t_586[k] = -5.0 * ki_418[k]
                   + f_0 * mi_782[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, t_591, t_592, t_593, ki_419, mi_783, \
                         mi_812, mi_813, mi_814, mi_815, mi_816, \
                         mi_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = -5.0 * ki_419[k]
                   + f_0 * mi_783[k];

        t_588[k] = f_0 * mi_812[k];

        t_589[k] = f_0 * mi_813[k];

        t_590[k] = f_0 * mi_814[k];

        t_591[k] = f_0 * mi_815[k];

        t_592[k] = f_0 * mi_816[k];

        t_593[k] = f_0 * mi_817[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, t_598, t_599, t_600, t_601, mi_818, \
                         mi_819, mi_820, mi_821, mi_822, mi_823, mi_824, \
                         mi_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_0 * mi_818[k];

        t_595[k] = f_0 * mi_819[k];

        t_596[k] = f_0 * mi_820[k];

        t_597[k] = f_0 * mi_821[k];

        t_598[k] = f_0 * mi_822[k];

        t_599[k] = f_0 * mi_823[k];

        t_600[k] = f_0 * mi_824[k];

        t_601[k] = f_0 * mi_825[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, t_607, t_608, t_609, mi_826, \
                         mi_827, mi_828, mi_829, mi_830, mi_831, mi_832, \
                         mi_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_0 * mi_826[k];

        t_603[k] = f_0 * mi_827[k];

        t_604[k] = f_0 * mi_828[k];

        t_605[k] = f_0 * mi_829[k];

        t_606[k] = f_0 * mi_830[k];

        t_607[k] = f_0 * mi_831[k];

        t_608[k] = f_0 * mi_832[k];

        t_609[k] = f_0 * mi_833[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, t_615, t_616, ki_420, mi_834, \
                         mi_835, mi_836, mi_837, mi_838, mi_839, \
                         mi_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_0 * mi_834[k];

        t_611[k] = f_0 * mi_835[k];

        t_612[k] = f_0 * mi_836[k];

        t_613[k] = f_0 * mi_837[k];

        t_614[k] = f_0 * mi_838[k];

        t_615[k] = f_0 * mi_839[k];

        t_616[k] = -ki_420[k]
                   + f_0 * mi_840[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, ki_421, ki_422, ki_423, ki_424, \
                         ki_425, mi_841, mi_842, mi_843, mi_844, \
                         mi_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = -ki_421[k]
                   + f_0 * mi_841[k];

        t_618[k] = -ki_422[k]
                   + f_0 * mi_842[k];

        t_619[k] = -ki_423[k]
                   + f_0 * mi_843[k];

        t_620[k] = -ki_424[k]
                   + f_0 * mi_844[k];

        t_621[k] = -ki_425[k]
                   + f_0 * mi_845[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, ki_426, ki_427, ki_428, ki_429, \
                         ki_430, mi_846, mi_847, mi_848, mi_849, \
                         mi_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = -ki_426[k]
                   + f_0 * mi_846[k];

        t_623[k] = -ki_427[k]
                   + f_0 * mi_847[k];

        t_624[k] = -ki_428[k]
                   + f_0 * mi_848[k];

        t_625[k] = -ki_429[k]
                   + f_0 * mi_849[k];

        t_626[k] = -ki_430[k]
                   + f_0 * mi_850[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, t_631, ki_431, ki_432, ki_433, ki_434, \
                         ki_435, mi_851, mi_852, mi_853, mi_854, \
                         mi_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = -ki_431[k]
                   + f_0 * mi_851[k];

        t_628[k] = -ki_432[k]
                   + f_0 * mi_852[k];

        t_629[k] = -ki_433[k]
                   + f_0 * mi_853[k];

        t_630[k] = -ki_434[k]
                   + f_0 * mi_854[k];

        t_631[k] = -ki_435[k]
                   + f_0 * mi_855[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, t_636, ki_436, ki_437, ki_438, ki_439, \
                         ki_440, mi_856, mi_857, mi_858, mi_859, \
                         mi_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = -ki_436[k]
                   + f_0 * mi_856[k];

        t_633[k] = -ki_437[k]
                   + f_0 * mi_857[k];

        t_634[k] = -ki_438[k]
                   + f_0 * mi_858[k];

        t_635[k] = -ki_439[k]
                   + f_0 * mi_859[k];

        t_636[k] = -ki_440[k]
                   + f_0 * mi_860[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, t_641, ki_441, ki_442, ki_443, ki_444, \
                         ki_445, mi_861, mi_862, mi_863, mi_864, \
                         mi_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = -ki_441[k]
                   + f_0 * mi_861[k];

        t_638[k] = -ki_442[k]
                   + f_0 * mi_862[k];

        t_639[k] = -ki_443[k]
                   + f_0 * mi_863[k];

        t_640[k] = -ki_444[k]
                   + f_0 * mi_864[k];

        t_641[k] = -ki_445[k]
                   + f_0 * mi_865[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, t_646, ki_446, ki_447, ki_448, ki_449, \
                         ki_450, mi_866, mi_867, mi_868, mi_869, \
                         mi_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = -ki_446[k]
                   + f_0 * mi_866[k];

        t_643[k] = -ki_447[k]
                   + f_0 * mi_867[k];

        t_644[k] = -2.0 * ki_448[k]
                   + f_0 * mi_868[k];

        t_645[k] = -2.0 * ki_449[k]
                   + f_0 * mi_869[k];

        t_646[k] = -2.0 * ki_450[k]
                   + f_0 * mi_870[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, t_651, ki_451, ki_452, ki_453, ki_454, \
                         ki_455, mi_871, mi_872, mi_873, mi_874, \
                         mi_875 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = -2.0 * ki_451[k]
                   + f_0 * mi_871[k];

        t_648[k] = -2.0 * ki_452[k]
                   + f_0 * mi_872[k];

        t_649[k] = -2.0 * ki_453[k]
                   + f_0 * mi_873[k];

        t_650[k] = -2.0 * ki_454[k]
                   + f_0 * mi_874[k];

        t_651[k] = -2.0 * ki_455[k]
                   + f_0 * mi_875[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, t_656, ki_456, ki_457, ki_458, ki_459, \
                         ki_460, mi_876, mi_877, mi_878, mi_879, \
                         mi_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = -2.0 * ki_456[k]
                   + f_0 * mi_876[k];

        t_653[k] = -2.0 * ki_457[k]
                   + f_0 * mi_877[k];

        t_654[k] = -2.0 * ki_458[k]
                   + f_0 * mi_878[k];

        t_655[k] = -2.0 * ki_459[k]
                   + f_0 * mi_879[k];

        t_656[k] = -2.0 * ki_460[k]
                   + f_0 * mi_880[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, ki_461, ki_462, ki_463, ki_464, \
                         ki_465, mi_881, mi_882, mi_883, mi_884, \
                         mi_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = -2.0 * ki_461[k]
                   + f_0 * mi_881[k];

        t_658[k] = -2.0 * ki_462[k]
                   + f_0 * mi_882[k];

        t_659[k] = -2.0 * ki_463[k]
                   + f_0 * mi_883[k];

        t_660[k] = -2.0 * ki_464[k]
                   + f_0 * mi_884[k];

        t_661[k] = -2.0 * ki_465[k]
                   + f_0 * mi_885[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, t_666, ki_466, ki_467, ki_468, ki_469, \
                         ki_470, mi_886, mi_887, mi_888, mi_889, \
                         mi_890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = -2.0 * ki_466[k]
                   + f_0 * mi_886[k];

        t_663[k] = -2.0 * ki_467[k]
                   + f_0 * mi_887[k];

        t_664[k] = -2.0 * ki_468[k]
                   + f_0 * mi_888[k];

        t_665[k] = -2.0 * ki_469[k]
                   + f_0 * mi_889[k];

        t_666[k] = -2.0 * ki_470[k]
                   + f_0 * mi_890[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ki_471 = buffer.data(ki + 471);
    const auto *ki_472 = buffer.data(ki + 472);
    const auto *ki_473 = buffer.data(ki + 473);
    const auto *ki_474 = buffer.data(ki + 474);
    const auto *ki_475 = buffer.data(ki + 475);
    const auto *ki_476 = buffer.data(ki + 476);
    const auto *ki_477 = buffer.data(ki + 477);
    const auto *ki_478 = buffer.data(ki + 478);
    const auto *ki_479 = buffer.data(ki + 479);
    const auto *ki_480 = buffer.data(ki + 480);
    const auto *ki_481 = buffer.data(ki + 481);
    const auto *ki_482 = buffer.data(ki + 482);
    const auto *ki_483 = buffer.data(ki + 483);
    const auto *ki_484 = buffer.data(ki + 484);
    const auto *ki_485 = buffer.data(ki + 485);
    const auto *ki_486 = buffer.data(ki + 486);
    const auto *ki_487 = buffer.data(ki + 487);
    const auto *ki_488 = buffer.data(ki + 488);
    const auto *ki_489 = buffer.data(ki + 489);
    const auto *ki_490 = buffer.data(ki + 490);
    const auto *ki_491 = buffer.data(ki + 491);
    const auto *ki_492 = buffer.data(ki + 492);
    const auto *ki_493 = buffer.data(ki + 493);
    const auto *ki_494 = buffer.data(ki + 494);
    const auto *ki_495 = buffer.data(ki + 495);
    const auto *ki_496 = buffer.data(ki + 496);
    const auto *ki_497 = buffer.data(ki + 497);
    const auto *ki_498 = buffer.data(ki + 498);
    const auto *ki_499 = buffer.data(ki + 499);
    const auto *ki_500 = buffer.data(ki + 500);
    const auto *ki_501 = buffer.data(ki + 501);
    const auto *ki_502 = buffer.data(ki + 502);
    const auto *ki_503 = buffer.data(ki + 503);
    const auto *ki_504 = buffer.data(ki + 504);
    const auto *ki_505 = buffer.data(ki + 505);
    const auto *ki_506 = buffer.data(ki + 506);
    const auto *ki_507 = buffer.data(ki + 507);
    const auto *ki_508 = buffer.data(ki + 508);
    const auto *ki_509 = buffer.data(ki + 509);
    const auto *ki_510 = buffer.data(ki + 510);
    const auto *ki_511 = buffer.data(ki + 511);
    const auto *ki_512 = buffer.data(ki + 512);
    const auto *ki_513 = buffer.data(ki + 513);
    const auto *ki_514 = buffer.data(ki + 514);
    const auto *ki_515 = buffer.data(ki + 515);
    const auto *ki_516 = buffer.data(ki + 516);
    const auto *ki_517 = buffer.data(ki + 517);
    const auto *ki_518 = buffer.data(ki + 518);
    const auto *ki_519 = buffer.data(ki + 519);
    const auto *ki_520 = buffer.data(ki + 520);
    const auto *ki_521 = buffer.data(ki + 521);
    const auto *ki_522 = buffer.data(ki + 522);
    const auto *ki_523 = buffer.data(ki + 523);
    const auto *ki_524 = buffer.data(ki + 524);
    const auto *ki_525 = buffer.data(ki + 525);
    const auto *ki_526 = buffer.data(ki + 526);
    const auto *ki_527 = buffer.data(ki + 527);
    const auto *ki_528 = buffer.data(ki + 528);
    const auto *ki_529 = buffer.data(ki + 529);
    const auto *ki_530 = buffer.data(ki + 530);
    const auto *ki_531 = buffer.data(ki + 531);
    const auto *ki_532 = buffer.data(ki + 532);
    const auto *ki_533 = buffer.data(ki + 533);
    const auto *ki_534 = buffer.data(ki + 534);
    const auto *ki_535 = buffer.data(ki + 535);
    const auto *ki_536 = buffer.data(ki + 536);
    const auto *ki_537 = buffer.data(ki + 537);
    const auto *ki_538 = buffer.data(ki + 538);
    const auto *ki_539 = buffer.data(ki + 539);
    const auto *ki_540 = buffer.data(ki + 540);
    const auto *ki_541 = buffer.data(ki + 541);
    const auto *ki_542 = buffer.data(ki + 542);
    const auto *ki_543 = buffer.data(ki + 543);
    const auto *ki_544 = buffer.data(ki + 544);
    const auto *ki_545 = buffer.data(ki + 545);
    const auto *ki_546 = buffer.data(ki + 546);
    const auto *ki_547 = buffer.data(ki + 547);
    const auto *ki_548 = buffer.data(ki + 548);
    const auto *ki_549 = buffer.data(ki + 549);
    const auto *ki_550 = buffer.data(ki + 550);
    const auto *ki_551 = buffer.data(ki + 551);
    const auto *ki_552 = buffer.data(ki + 552);
    const auto *ki_553 = buffer.data(ki + 553);
    const auto *ki_554 = buffer.data(ki + 554);
    const auto *ki_555 = buffer.data(ki + 555);
    const auto *ki_556 = buffer.data(ki + 556);
    const auto *ki_557 = buffer.data(ki + 557);
    const auto *ki_558 = buffer.data(ki + 558);
    const auto *ki_559 = buffer.data(ki + 559);
    const auto *ki_560 = buffer.data(ki + 560);
    const auto *ki_561 = buffer.data(ki + 561);
    const auto *ki_562 = buffer.data(ki + 562);
    const auto *ki_563 = buffer.data(ki + 563);
    const auto *ki_564 = buffer.data(ki + 564);
    const auto *ki_565 = buffer.data(ki + 565);
    const auto *ki_566 = buffer.data(ki + 566);
    const auto *ki_567 = buffer.data(ki + 567);
    const auto *ki_568 = buffer.data(ki + 568);
    const auto *ki_569 = buffer.data(ki + 569);
    const auto *ki_570 = buffer.data(ki + 570);
    const auto *ki_571 = buffer.data(ki + 571);
    const auto *ki_572 = buffer.data(ki + 572);
    const auto *ki_573 = buffer.data(ki + 573);
    const auto *ki_574 = buffer.data(ki + 574);
    const auto *ki_575 = buffer.data(ki + 575);
    const auto *ki_576 = buffer.data(ki + 576);
    const auto *ki_577 = buffer.data(ki + 577);
    const auto *ki_578 = buffer.data(ki + 578);
    const auto *ki_579 = buffer.data(ki + 579);
    const auto *ki_580 = buffer.data(ki + 580);
    const auto *ki_581 = buffer.data(ki + 581);
    const auto *ki_582 = buffer.data(ki + 582);
    const auto *ki_583 = buffer.data(ki + 583);
    const auto *ki_584 = buffer.data(ki + 584);
    const auto *ki_585 = buffer.data(ki + 585);
    const auto *ki_586 = buffer.data(ki + 586);
    const auto *ki_587 = buffer.data(ki + 587);
    const auto *ki_588 = buffer.data(ki + 588);
    const auto *ki_589 = buffer.data(ki + 589);
    const auto *ki_590 = buffer.data(ki + 590);
    const auto *ki_591 = buffer.data(ki + 591);
    const auto *ki_592 = buffer.data(ki + 592);
    const auto *ki_593 = buffer.data(ki + 593);
    const auto *ki_594 = buffer.data(ki + 594);
    const auto *ki_595 = buffer.data(ki + 595);
    const auto *ki_596 = buffer.data(ki + 596);
    const auto *ki_597 = buffer.data(ki + 597);
    const auto *ki_598 = buffer.data(ki + 598);
    const auto *ki_599 = buffer.data(ki + 599);
    const auto *ki_600 = buffer.data(ki + 600);
    const auto *ki_601 = buffer.data(ki + 601);
    const auto *ki_602 = buffer.data(ki + 602);

    const auto *mi_891 = buffer.data(mi + 891);
    const auto *mi_892 = buffer.data(mi + 892);
    const auto *mi_893 = buffer.data(mi + 893);
    const auto *mi_894 = buffer.data(mi + 894);
    const auto *mi_895 = buffer.data(mi + 895);
    const auto *mi_896 = buffer.data(mi + 896);
    const auto *mi_897 = buffer.data(mi + 897);
    const auto *mi_898 = buffer.data(mi + 898);
    const auto *mi_899 = buffer.data(mi + 899);
    const auto *mi_900 = buffer.data(mi + 900);
    const auto *mi_901 = buffer.data(mi + 901);
    const auto *mi_902 = buffer.data(mi + 902);
    const auto *mi_903 = buffer.data(mi + 903);
    const auto *mi_904 = buffer.data(mi + 904);
    const auto *mi_905 = buffer.data(mi + 905);
    const auto *mi_906 = buffer.data(mi + 906);
    const auto *mi_907 = buffer.data(mi + 907);
    const auto *mi_908 = buffer.data(mi + 908);
    const auto *mi_909 = buffer.data(mi + 909);
    const auto *mi_910 = buffer.data(mi + 910);
    const auto *mi_911 = buffer.data(mi + 911);
    const auto *mi_912 = buffer.data(mi + 912);
    const auto *mi_913 = buffer.data(mi + 913);
    const auto *mi_914 = buffer.data(mi + 914);
    const auto *mi_915 = buffer.data(mi + 915);
    const auto *mi_916 = buffer.data(mi + 916);
    const auto *mi_917 = buffer.data(mi + 917);
    const auto *mi_918 = buffer.data(mi + 918);
    const auto *mi_919 = buffer.data(mi + 919);
    const auto *mi_920 = buffer.data(mi + 920);
    const auto *mi_921 = buffer.data(mi + 921);
    const auto *mi_922 = buffer.data(mi + 922);
    const auto *mi_923 = buffer.data(mi + 923);
    const auto *mi_924 = buffer.data(mi + 924);
    const auto *mi_925 = buffer.data(mi + 925);
    const auto *mi_926 = buffer.data(mi + 926);
    const auto *mi_927 = buffer.data(mi + 927);
    const auto *mi_928 = buffer.data(mi + 928);
    const auto *mi_929 = buffer.data(mi + 929);
    const auto *mi_930 = buffer.data(mi + 930);
    const auto *mi_931 = buffer.data(mi + 931);
    const auto *mi_932 = buffer.data(mi + 932);
    const auto *mi_933 = buffer.data(mi + 933);
    const auto *mi_934 = buffer.data(mi + 934);
    const auto *mi_935 = buffer.data(mi + 935);
    const auto *mi_936 = buffer.data(mi + 936);
    const auto *mi_937 = buffer.data(mi + 937);
    const auto *mi_938 = buffer.data(mi + 938);
    const auto *mi_939 = buffer.data(mi + 939);
    const auto *mi_940 = buffer.data(mi + 940);
    const auto *mi_941 = buffer.data(mi + 941);
    const auto *mi_942 = buffer.data(mi + 942);
    const auto *mi_943 = buffer.data(mi + 943);
    const auto *mi_944 = buffer.data(mi + 944);
    const auto *mi_945 = buffer.data(mi + 945);
    const auto *mi_946 = buffer.data(mi + 946);
    const auto *mi_947 = buffer.data(mi + 947);
    const auto *mi_948 = buffer.data(mi + 948);
    const auto *mi_949 = buffer.data(mi + 949);
    const auto *mi_950 = buffer.data(mi + 950);
    const auto *mi_951 = buffer.data(mi + 951);
    const auto *mi_952 = buffer.data(mi + 952);
    const auto *mi_953 = buffer.data(mi + 953);
    const auto *mi_954 = buffer.data(mi + 954);
    const auto *mi_955 = buffer.data(mi + 955);
    const auto *mi_956 = buffer.data(mi + 956);
    const auto *mi_957 = buffer.data(mi + 957);
    const auto *mi_958 = buffer.data(mi + 958);
    const auto *mi_959 = buffer.data(mi + 959);
    const auto *mi_960 = buffer.data(mi + 960);
    const auto *mi_961 = buffer.data(mi + 961);
    const auto *mi_962 = buffer.data(mi + 962);
    const auto *mi_963 = buffer.data(mi + 963);
    const auto *mi_964 = buffer.data(mi + 964);
    const auto *mi_965 = buffer.data(mi + 965);
    const auto *mi_966 = buffer.data(mi + 966);
    const auto *mi_967 = buffer.data(mi + 967);
    const auto *mi_968 = buffer.data(mi + 968);
    const auto *mi_969 = buffer.data(mi + 969);
    const auto *mi_970 = buffer.data(mi + 970);
    const auto *mi_971 = buffer.data(mi + 971);
    const auto *mi_972 = buffer.data(mi + 972);
    const auto *mi_973 = buffer.data(mi + 973);
    const auto *mi_974 = buffer.data(mi + 974);
    const auto *mi_975 = buffer.data(mi + 975);
    const auto *mi_976 = buffer.data(mi + 976);
    const auto *mi_977 = buffer.data(mi + 977);
    const auto *mi_978 = buffer.data(mi + 978);
    const auto *mi_979 = buffer.data(mi + 979);
    const auto *mi_980 = buffer.data(mi + 980);
    const auto *mi_981 = buffer.data(mi + 981);
    const auto *mi_982 = buffer.data(mi + 982);
    const auto *mi_983 = buffer.data(mi + 983);
    const auto *mi_984 = buffer.data(mi + 984);
    const auto *mi_985 = buffer.data(mi + 985);
    const auto *mi_986 = buffer.data(mi + 986);
    const auto *mi_987 = buffer.data(mi + 987);
    const auto *mi_988 = buffer.data(mi + 988);
    const auto *mi_989 = buffer.data(mi + 989);
    const auto *mi_990 = buffer.data(mi + 990);
    const auto *mi_991 = buffer.data(mi + 991);
    const auto *mi_992 = buffer.data(mi + 992);
    const auto *mi_993 = buffer.data(mi + 993);
    const auto *mi_994 = buffer.data(mi + 994);
    const auto *mi_995 = buffer.data(mi + 995);
    const auto *mi_996 = buffer.data(mi + 996);
    const auto *mi_997 = buffer.data(mi + 997);
    const auto *mi_998 = buffer.data(mi + 998);
    const auto *mi_999 = buffer.data(mi + 999);
    const auto *mi_1000 = buffer.data(mi + 1000);
    const auto *mi_1001 = buffer.data(mi + 1001);
    const auto *mi_1002 = buffer.data(mi + 1002);
    const auto *mi_1003 = buffer.data(mi + 1003);
    const auto *mi_1004 = buffer.data(mi + 1004);
    const auto *mi_1005 = buffer.data(mi + 1005);
    const auto *mi_1006 = buffer.data(mi + 1006);
    const auto *mi_1007 = buffer.data(mi + 1007);
    const auto *mi_1036 = buffer.data(mi + 1036);
    const auto *mi_1037 = buffer.data(mi + 1037);
    const auto *mi_1038 = buffer.data(mi + 1038);
    const auto *mi_1039 = buffer.data(mi + 1039);
    const auto *mi_1040 = buffer.data(mi + 1040);
    const auto *mi_1041 = buffer.data(mi + 1041);
    const auto *mi_1042 = buffer.data(mi + 1042);
    const auto *mi_1043 = buffer.data(mi + 1043);
    const auto *mi_1044 = buffer.data(mi + 1044);
    const auto *mi_1045 = buffer.data(mi + 1045);
    const auto *mi_1046 = buffer.data(mi + 1046);
    const auto *mi_1047 = buffer.data(mi + 1047);
    const auto *mi_1048 = buffer.data(mi + 1048);
    const auto *mi_1049 = buffer.data(mi + 1049);
    const auto *mi_1050 = buffer.data(mi + 1050);
    const auto *mi_1051 = buffer.data(mi + 1051);
    const auto *mi_1052 = buffer.data(mi + 1052);
    const auto *mi_1053 = buffer.data(mi + 1053);
    const auto *mi_1054 = buffer.data(mi + 1054);
    const auto *mi_1055 = buffer.data(mi + 1055);
    const auto *mi_1056 = buffer.data(mi + 1056);
    const auto *mi_1057 = buffer.data(mi + 1057);
    const auto *mi_1058 = buffer.data(mi + 1058);
    const auto *mi_1059 = buffer.data(mi + 1059);
    const auto *mi_1060 = buffer.data(mi + 1060);
    const auto *mi_1061 = buffer.data(mi + 1061);
    const auto *mi_1062 = buffer.data(mi + 1062);
    const auto *mi_1063 = buffer.data(mi + 1063);
    const auto *mi_1064 = buffer.data(mi + 1064);
    const auto *mi_1065 = buffer.data(mi + 1065);
    const auto *mi_1066 = buffer.data(mi + 1066);
    const auto *mi_1067 = buffer.data(mi + 1067);
    const auto *mi_1068 = buffer.data(mi + 1068);
    const auto *mi_1069 = buffer.data(mi + 1069);
    const auto *mi_1070 = buffer.data(mi + 1070);
    const auto *mi_1071 = buffer.data(mi + 1071);
    const auto *mi_1072 = buffer.data(mi + 1072);
    const auto *mi_1073 = buffer.data(mi + 1073);
    const auto *mi_1074 = buffer.data(mi + 1074);
    const auto *mi_1075 = buffer.data(mi + 1075);
    const auto *mi_1076 = buffer.data(mi + 1076);
    const auto *mi_1077 = buffer.data(mi + 1077);
    const auto *mi_1078 = buffer.data(mi + 1078);

#pragma omp simd aligned(t_667, t_668, t_669, t_670, t_671, ki_471, ki_472, ki_473, ki_474, \
                         ki_475, mi_891, mi_892, mi_893, mi_894, \
                         mi_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = -2.0 * ki_471[k]
                   + f_0 * mi_891[k];

        t_668[k] = -2.0 * ki_472[k]
                   + f_0 * mi_892[k];

        t_669[k] = -2.0 * ki_473[k]
                   + f_0 * mi_893[k];

        t_670[k] = -2.0 * ki_474[k]
                   + f_0 * mi_894[k];

        t_671[k] = -2.0 * ki_475[k]
                   + f_0 * mi_895[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, t_676, ki_476, ki_477, ki_478, ki_479, \
                         ki_480, mi_896, mi_897, mi_898, mi_899, \
                         mi_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = -3.0 * ki_476[k]
                   + f_0 * mi_896[k];

        t_673[k] = -3.0 * ki_477[k]
                   + f_0 * mi_897[k];

        t_674[k] = -3.0 * ki_478[k]
                   + f_0 * mi_898[k];

        t_675[k] = -3.0 * ki_479[k]
                   + f_0 * mi_899[k];

        t_676[k] = -3.0 * ki_480[k]
                   + f_0 * mi_900[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, t_680, t_681, ki_481, ki_482, ki_483, ki_484, \
                         ki_485, mi_901, mi_902, mi_903, mi_904, \
                         mi_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = -3.0 * ki_481[k]
                   + f_0 * mi_901[k];

        t_678[k] = -3.0 * ki_482[k]
                   + f_0 * mi_902[k];

        t_679[k] = -3.0 * ki_483[k]
                   + f_0 * mi_903[k];

        t_680[k] = -3.0 * ki_484[k]
                   + f_0 * mi_904[k];

        t_681[k] = -3.0 * ki_485[k]
                   + f_0 * mi_905[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, t_686, ki_486, ki_487, ki_488, ki_489, \
                         ki_490, mi_906, mi_907, mi_908, mi_909, \
                         mi_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = -3.0 * ki_486[k]
                   + f_0 * mi_906[k];

        t_683[k] = -3.0 * ki_487[k]
                   + f_0 * mi_907[k];

        t_684[k] = -3.0 * ki_488[k]
                   + f_0 * mi_908[k];

        t_685[k] = -3.0 * ki_489[k]
                   + f_0 * mi_909[k];

        t_686[k] = -3.0 * ki_490[k]
                   + f_0 * mi_910[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, t_691, ki_491, ki_492, ki_493, ki_494, \
                         ki_495, mi_911, mi_912, mi_913, mi_914, \
                         mi_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = -3.0 * ki_491[k]
                   + f_0 * mi_911[k];

        t_688[k] = -3.0 * ki_492[k]
                   + f_0 * mi_912[k];

        t_689[k] = -3.0 * ki_493[k]
                   + f_0 * mi_913[k];

        t_690[k] = -3.0 * ki_494[k]
                   + f_0 * mi_914[k];

        t_691[k] = -3.0 * ki_495[k]
                   + f_0 * mi_915[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, t_696, ki_496, ki_497, ki_498, ki_499, \
                         ki_500, mi_916, mi_917, mi_918, mi_919, \
                         mi_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = -3.0 * ki_496[k]
                   + f_0 * mi_916[k];

        t_693[k] = -3.0 * ki_497[k]
                   + f_0 * mi_917[k];

        t_694[k] = -3.0 * ki_498[k]
                   + f_0 * mi_918[k];

        t_695[k] = -3.0 * ki_499[k]
                   + f_0 * mi_919[k];

        t_696[k] = -3.0 * ki_500[k]
                   + f_0 * mi_920[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, t_701, ki_501, ki_502, ki_503, ki_504, \
                         ki_505, mi_921, mi_922, mi_923, mi_924, \
                         mi_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = -3.0 * ki_501[k]
                   + f_0 * mi_921[k];

        t_698[k] = -3.0 * ki_502[k]
                   + f_0 * mi_922[k];

        t_699[k] = -3.0 * ki_503[k]
                   + f_0 * mi_923[k];

        t_700[k] = -4.0 * ki_504[k]
                   + f_0 * mi_924[k];

        t_701[k] = -4.0 * ki_505[k]
                   + f_0 * mi_925[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, t_706, ki_506, ki_507, ki_508, ki_509, \
                         ki_510, mi_926, mi_927, mi_928, mi_929, \
                         mi_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = -4.0 * ki_506[k]
                   + f_0 * mi_926[k];

        t_703[k] = -4.0 * ki_507[k]
                   + f_0 * mi_927[k];

        t_704[k] = -4.0 * ki_508[k]
                   + f_0 * mi_928[k];

        t_705[k] = -4.0 * ki_509[k]
                   + f_0 * mi_929[k];

        t_706[k] = -4.0 * ki_510[k]
                   + f_0 * mi_930[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, t_710, t_711, ki_511, ki_512, ki_513, ki_514, \
                         ki_515, mi_931, mi_932, mi_933, mi_934, \
                         mi_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = -4.0 * ki_511[k]
                   + f_0 * mi_931[k];

        t_708[k] = -4.0 * ki_512[k]
                   + f_0 * mi_932[k];

        t_709[k] = -4.0 * ki_513[k]
                   + f_0 * mi_933[k];

        t_710[k] = -4.0 * ki_514[k]
                   + f_0 * mi_934[k];

        t_711[k] = -4.0 * ki_515[k]
                   + f_0 * mi_935[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, t_716, ki_516, ki_517, ki_518, ki_519, \
                         ki_520, mi_936, mi_937, mi_938, mi_939, \
                         mi_940 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = -4.0 * ki_516[k]
                   + f_0 * mi_936[k];

        t_713[k] = -4.0 * ki_517[k]
                   + f_0 * mi_937[k];

        t_714[k] = -4.0 * ki_518[k]
                   + f_0 * mi_938[k];

        t_715[k] = -4.0 * ki_519[k]
                   + f_0 * mi_939[k];

        t_716[k] = -4.0 * ki_520[k]
                   + f_0 * mi_940[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, t_720, t_721, ki_521, ki_522, ki_523, ki_524, \
                         ki_525, mi_941, mi_942, mi_943, mi_944, \
                         mi_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = -4.0 * ki_521[k]
                   + f_0 * mi_941[k];

        t_718[k] = -4.0 * ki_522[k]
                   + f_0 * mi_942[k];

        t_719[k] = -4.0 * ki_523[k]
                   + f_0 * mi_943[k];

        t_720[k] = -4.0 * ki_524[k]
                   + f_0 * mi_944[k];

        t_721[k] = -4.0 * ki_525[k]
                   + f_0 * mi_945[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, t_725, t_726, ki_526, ki_527, ki_528, ki_529, \
                         ki_530, mi_946, mi_947, mi_948, mi_949, \
                         mi_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = -4.0 * ki_526[k]
                   + f_0 * mi_946[k];

        t_723[k] = -4.0 * ki_527[k]
                   + f_0 * mi_947[k];

        t_724[k] = -4.0 * ki_528[k]
                   + f_0 * mi_948[k];

        t_725[k] = -4.0 * ki_529[k]
                   + f_0 * mi_949[k];

        t_726[k] = -4.0 * ki_530[k]
                   + f_0 * mi_950[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, t_730, t_731, ki_531, ki_532, ki_533, ki_534, \
                         ki_535, mi_951, mi_952, mi_953, mi_954, \
                         mi_955 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = -4.0 * ki_531[k]
                   + f_0 * mi_951[k];

        t_728[k] = -5.0 * ki_532[k]
                   + f_0 * mi_952[k];

        t_729[k] = -5.0 * ki_533[k]
                   + f_0 * mi_953[k];

        t_730[k] = -5.0 * ki_534[k]
                   + f_0 * mi_954[k];

        t_731[k] = -5.0 * ki_535[k]
                   + f_0 * mi_955[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, t_735, t_736, ki_536, ki_537, ki_538, ki_539, \
                         ki_540, mi_956, mi_957, mi_958, mi_959, \
                         mi_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = -5.0 * ki_536[k]
                   + f_0 * mi_956[k];

        t_733[k] = -5.0 * ki_537[k]
                   + f_0 * mi_957[k];

        t_734[k] = -5.0 * ki_538[k]
                   + f_0 * mi_958[k];

        t_735[k] = -5.0 * ki_539[k]
                   + f_0 * mi_959[k];

        t_736[k] = -5.0 * ki_540[k]
                   + f_0 * mi_960[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, t_740, t_741, ki_541, ki_542, ki_543, ki_544, \
                         ki_545, mi_961, mi_962, mi_963, mi_964, \
                         mi_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = -5.0 * ki_541[k]
                   + f_0 * mi_961[k];

        t_738[k] = -5.0 * ki_542[k]
                   + f_0 * mi_962[k];

        t_739[k] = -5.0 * ki_543[k]
                   + f_0 * mi_963[k];

        t_740[k] = -5.0 * ki_544[k]
                   + f_0 * mi_964[k];

        t_741[k] = -5.0 * ki_545[k]
                   + f_0 * mi_965[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, t_746, ki_546, ki_547, ki_548, ki_549, \
                         ki_550, mi_966, mi_967, mi_968, mi_969, \
                         mi_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = -5.0 * ki_546[k]
                   + f_0 * mi_966[k];

        t_743[k] = -5.0 * ki_547[k]
                   + f_0 * mi_967[k];

        t_744[k] = -5.0 * ki_548[k]
                   + f_0 * mi_968[k];

        t_745[k] = -5.0 * ki_549[k]
                   + f_0 * mi_969[k];

        t_746[k] = -5.0 * ki_550[k]
                   + f_0 * mi_970[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, t_750, t_751, ki_551, ki_552, ki_553, ki_554, \
                         ki_555, mi_971, mi_972, mi_973, mi_974, \
                         mi_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = -5.0 * ki_551[k]
                   + f_0 * mi_971[k];

        t_748[k] = -5.0 * ki_552[k]
                   + f_0 * mi_972[k];

        t_749[k] = -5.0 * ki_553[k]
                   + f_0 * mi_973[k];

        t_750[k] = -5.0 * ki_554[k]
                   + f_0 * mi_974[k];

        t_751[k] = -5.0 * ki_555[k]
                   + f_0 * mi_975[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, t_756, ki_556, ki_557, ki_558, ki_559, \
                         ki_560, mi_976, mi_977, mi_978, mi_979, \
                         mi_980 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = -5.0 * ki_556[k]
                   + f_0 * mi_976[k];

        t_753[k] = -5.0 * ki_557[k]
                   + f_0 * mi_977[k];

        t_754[k] = -5.0 * ki_558[k]
                   + f_0 * mi_978[k];

        t_755[k] = -5.0 * ki_559[k]
                   + f_0 * mi_979[k];

        t_756[k] = -6.0 * ki_560[k]
                   + f_0 * mi_980[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, t_761, ki_561, ki_562, ki_563, ki_564, \
                         ki_565, mi_981, mi_982, mi_983, mi_984, \
                         mi_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = -6.0 * ki_561[k]
                   + f_0 * mi_981[k];

        t_758[k] = -6.0 * ki_562[k]
                   + f_0 * mi_982[k];

        t_759[k] = -6.0 * ki_563[k]
                   + f_0 * mi_983[k];

        t_760[k] = -6.0 * ki_564[k]
                   + f_0 * mi_984[k];

        t_761[k] = -6.0 * ki_565[k]
                   + f_0 * mi_985[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, t_766, ki_566, ki_567, ki_568, ki_569, \
                         ki_570, mi_986, mi_987, mi_988, mi_989, \
                         mi_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = -6.0 * ki_566[k]
                   + f_0 * mi_986[k];

        t_763[k] = -6.0 * ki_567[k]
                   + f_0 * mi_987[k];

        t_764[k] = -6.0 * ki_568[k]
                   + f_0 * mi_988[k];

        t_765[k] = -6.0 * ki_569[k]
                   + f_0 * mi_989[k];

        t_766[k] = -6.0 * ki_570[k]
                   + f_0 * mi_990[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, t_770, t_771, ki_571, ki_572, ki_573, ki_574, \
                         ki_575, mi_991, mi_992, mi_993, mi_994, \
                         mi_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = -6.0 * ki_571[k]
                   + f_0 * mi_991[k];

        t_768[k] = -6.0 * ki_572[k]
                   + f_0 * mi_992[k];

        t_769[k] = -6.0 * ki_573[k]
                   + f_0 * mi_993[k];

        t_770[k] = -6.0 * ki_574[k]
                   + f_0 * mi_994[k];

        t_771[k] = -6.0 * ki_575[k]
                   + f_0 * mi_995[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, ki_576, ki_577, ki_578, ki_579, \
                         ki_580, mi_996, mi_997, mi_998, mi_999, \
                         mi_1000 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = -6.0 * ki_576[k]
                   + f_0 * mi_996[k];

        t_773[k] = -6.0 * ki_577[k]
                   + f_0 * mi_997[k];

        t_774[k] = -6.0 * ki_578[k]
                   + f_0 * mi_998[k];

        t_775[k] = -6.0 * ki_579[k]
                   + f_0 * mi_999[k];

        t_776[k] = -6.0 * ki_580[k]
                   + f_0 * mi_1000[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, t_780, t_781, ki_581, ki_582, ki_583, ki_584, \
                         ki_585, mi_1001, mi_1002, mi_1003, mi_1004, \
                         mi_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = -6.0 * ki_581[k]
                   + f_0 * mi_1001[k];

        t_778[k] = -6.0 * ki_582[k]
                   + f_0 * mi_1002[k];

        t_779[k] = -6.0 * ki_583[k]
                   + f_0 * mi_1003[k];

        t_780[k] = -6.0 * ki_584[k]
                   + f_0 * mi_1004[k];

        t_781[k] = -6.0 * ki_585[k]
                   + f_0 * mi_1005[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, t_785, t_786, t_787, t_788, ki_586, ki_587, \
                         mi_1006, mi_1007, mi_1036, mi_1037, mi_1038, mi_1039, \
                         mi_1040 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = -6.0 * ki_586[k]
                   + f_0 * mi_1006[k];

        t_783[k] = -6.0 * ki_587[k]
                   + f_0 * mi_1007[k];

        t_784[k] = f_0 * mi_1036[k];

        t_785[k] = f_0 * mi_1037[k];

        t_786[k] = f_0 * mi_1038[k];

        t_787[k] = f_0 * mi_1039[k];

        t_788[k] = f_0 * mi_1040[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, t_792, t_793, t_794, t_795, t_796, mi_1041, \
                         mi_1042, mi_1043, mi_1044, mi_1045, mi_1046, mi_1047, \
                         mi_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_0 * mi_1041[k];

        t_790[k] = f_0 * mi_1042[k];

        t_791[k] = f_0 * mi_1043[k];

        t_792[k] = f_0 * mi_1044[k];

        t_793[k] = f_0 * mi_1045[k];

        t_794[k] = f_0 * mi_1046[k];

        t_795[k] = f_0 * mi_1047[k];

        t_796[k] = f_0 * mi_1048[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, t_800, t_801, t_802, t_803, t_804, mi_1049, \
                         mi_1050, mi_1051, mi_1052, mi_1053, mi_1054, mi_1055, \
                         mi_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_0 * mi_1049[k];

        t_798[k] = f_0 * mi_1050[k];

        t_799[k] = f_0 * mi_1051[k];

        t_800[k] = f_0 * mi_1052[k];

        t_801[k] = f_0 * mi_1053[k];

        t_802[k] = f_0 * mi_1054[k];

        t_803[k] = f_0 * mi_1055[k];

        t_804[k] = f_0 * mi_1056[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, t_810, t_811, mi_1057, mi_1058, \
                         mi_1059, mi_1060, mi_1061, mi_1062, mi_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_0 * mi_1057[k];

        t_806[k] = f_0 * mi_1058[k];

        t_807[k] = f_0 * mi_1059[k];

        t_808[k] = f_0 * mi_1060[k];

        t_809[k] = f_0 * mi_1061[k];

        t_810[k] = f_0 * mi_1062[k];

        t_811[k] = f_0 * mi_1063[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, t_816, ki_588, ki_589, ki_590, ki_591, \
                         ki_592, mi_1064, mi_1065, mi_1066, mi_1067, \
                         mi_1068 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = -ki_588[k]
                   + f_0 * mi_1064[k];

        t_813[k] = -ki_589[k]
                   + f_0 * mi_1065[k];

        t_814[k] = -ki_590[k]
                   + f_0 * mi_1066[k];

        t_815[k] = -ki_591[k]
                   + f_0 * mi_1067[k];

        t_816[k] = -ki_592[k]
                   + f_0 * mi_1068[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, t_820, t_821, ki_593, ki_594, ki_595, ki_596, \
                         ki_597, mi_1069, mi_1070, mi_1071, mi_1072, \
                         mi_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = -ki_593[k]
                   + f_0 * mi_1069[k];

        t_818[k] = -ki_594[k]
                   + f_0 * mi_1070[k];

        t_819[k] = -ki_595[k]
                   + f_0 * mi_1071[k];

        t_820[k] = -ki_596[k]
                   + f_0 * mi_1072[k];

        t_821[k] = -ki_597[k]
                   + f_0 * mi_1073[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, t_825, t_826, ki_598, ki_599, ki_600, ki_601, \
                         ki_602, mi_1074, mi_1075, mi_1076, mi_1077, \
                         mi_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = -ki_598[k]
                   + f_0 * mi_1074[k];

        t_823[k] = -ki_599[k]
                   + f_0 * mi_1075[k];

        t_824[k] = -ki_600[k]
                   + f_0 * mi_1076[k];

        t_825[k] = -ki_601[k]
                   + f_0 * mi_1077[k];

        t_826[k] = -ki_602[k]
                   + f_0 * mi_1078[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_2_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);
    auto *t_948 = buffer.data(target + 948);
    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);
    auto *t_955 = buffer.data(target + 955);
    auto *t_956 = buffer.data(target + 956);
    auto *t_957 = buffer.data(target + 957);
    auto *t_958 = buffer.data(target + 958);
    auto *t_959 = buffer.data(target + 959);
    auto *t_960 = buffer.data(target + 960);
    auto *t_961 = buffer.data(target + 961);
    auto *t_962 = buffer.data(target + 962);
    auto *t_963 = buffer.data(target + 963);
    auto *t_964 = buffer.data(target + 964);
    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);
    auto *t_970 = buffer.data(target + 970);
    auto *t_971 = buffer.data(target + 971);
    auto *t_972 = buffer.data(target + 972);
    auto *t_973 = buffer.data(target + 973);
    auto *t_974 = buffer.data(target + 974);
    auto *t_975 = buffer.data(target + 975);
    auto *t_976 = buffer.data(target + 976);

    const auto *ki_603 = buffer.data(ki + 603);
    const auto *ki_604 = buffer.data(ki + 604);
    const auto *ki_605 = buffer.data(ki + 605);
    const auto *ki_606 = buffer.data(ki + 606);
    const auto *ki_607 = buffer.data(ki + 607);
    const auto *ki_608 = buffer.data(ki + 608);
    const auto *ki_609 = buffer.data(ki + 609);
    const auto *ki_610 = buffer.data(ki + 610);
    const auto *ki_611 = buffer.data(ki + 611);
    const auto *ki_612 = buffer.data(ki + 612);
    const auto *ki_613 = buffer.data(ki + 613);
    const auto *ki_614 = buffer.data(ki + 614);
    const auto *ki_615 = buffer.data(ki + 615);
    const auto *ki_616 = buffer.data(ki + 616);
    const auto *ki_617 = buffer.data(ki + 617);
    const auto *ki_618 = buffer.data(ki + 618);
    const auto *ki_619 = buffer.data(ki + 619);
    const auto *ki_620 = buffer.data(ki + 620);
    const auto *ki_621 = buffer.data(ki + 621);
    const auto *ki_622 = buffer.data(ki + 622);
    const auto *ki_623 = buffer.data(ki + 623);
    const auto *ki_624 = buffer.data(ki + 624);
    const auto *ki_625 = buffer.data(ki + 625);
    const auto *ki_626 = buffer.data(ki + 626);
    const auto *ki_627 = buffer.data(ki + 627);
    const auto *ki_628 = buffer.data(ki + 628);
    const auto *ki_629 = buffer.data(ki + 629);
    const auto *ki_630 = buffer.data(ki + 630);
    const auto *ki_631 = buffer.data(ki + 631);
    const auto *ki_632 = buffer.data(ki + 632);
    const auto *ki_633 = buffer.data(ki + 633);
    const auto *ki_634 = buffer.data(ki + 634);
    const auto *ki_635 = buffer.data(ki + 635);
    const auto *ki_636 = buffer.data(ki + 636);
    const auto *ki_637 = buffer.data(ki + 637);
    const auto *ki_638 = buffer.data(ki + 638);
    const auto *ki_639 = buffer.data(ki + 639);
    const auto *ki_640 = buffer.data(ki + 640);
    const auto *ki_641 = buffer.data(ki + 641);
    const auto *ki_642 = buffer.data(ki + 642);
    const auto *ki_643 = buffer.data(ki + 643);
    const auto *ki_644 = buffer.data(ki + 644);
    const auto *ki_645 = buffer.data(ki + 645);
    const auto *ki_646 = buffer.data(ki + 646);
    const auto *ki_647 = buffer.data(ki + 647);
    const auto *ki_648 = buffer.data(ki + 648);
    const auto *ki_649 = buffer.data(ki + 649);
    const auto *ki_650 = buffer.data(ki + 650);
    const auto *ki_651 = buffer.data(ki + 651);
    const auto *ki_652 = buffer.data(ki + 652);
    const auto *ki_653 = buffer.data(ki + 653);
    const auto *ki_654 = buffer.data(ki + 654);
    const auto *ki_655 = buffer.data(ki + 655);
    const auto *ki_656 = buffer.data(ki + 656);
    const auto *ki_657 = buffer.data(ki + 657);
    const auto *ki_658 = buffer.data(ki + 658);
    const auto *ki_659 = buffer.data(ki + 659);
    const auto *ki_660 = buffer.data(ki + 660);
    const auto *ki_661 = buffer.data(ki + 661);
    const auto *ki_662 = buffer.data(ki + 662);
    const auto *ki_663 = buffer.data(ki + 663);
    const auto *ki_664 = buffer.data(ki + 664);
    const auto *ki_665 = buffer.data(ki + 665);
    const auto *ki_666 = buffer.data(ki + 666);
    const auto *ki_667 = buffer.data(ki + 667);
    const auto *ki_668 = buffer.data(ki + 668);
    const auto *ki_669 = buffer.data(ki + 669);
    const auto *ki_670 = buffer.data(ki + 670);
    const auto *ki_671 = buffer.data(ki + 671);
    const auto *ki_672 = buffer.data(ki + 672);
    const auto *ki_673 = buffer.data(ki + 673);
    const auto *ki_674 = buffer.data(ki + 674);
    const auto *ki_675 = buffer.data(ki + 675);
    const auto *ki_676 = buffer.data(ki + 676);
    const auto *ki_677 = buffer.data(ki + 677);
    const auto *ki_678 = buffer.data(ki + 678);
    const auto *ki_679 = buffer.data(ki + 679);
    const auto *ki_680 = buffer.data(ki + 680);
    const auto *ki_681 = buffer.data(ki + 681);
    const auto *ki_682 = buffer.data(ki + 682);
    const auto *ki_683 = buffer.data(ki + 683);
    const auto *ki_684 = buffer.data(ki + 684);
    const auto *ki_685 = buffer.data(ki + 685);
    const auto *ki_686 = buffer.data(ki + 686);
    const auto *ki_687 = buffer.data(ki + 687);
    const auto *ki_688 = buffer.data(ki + 688);
    const auto *ki_689 = buffer.data(ki + 689);
    const auto *ki_690 = buffer.data(ki + 690);
    const auto *ki_691 = buffer.data(ki + 691);
    const auto *ki_692 = buffer.data(ki + 692);
    const auto *ki_693 = buffer.data(ki + 693);
    const auto *ki_694 = buffer.data(ki + 694);
    const auto *ki_695 = buffer.data(ki + 695);
    const auto *ki_696 = buffer.data(ki + 696);
    const auto *ki_697 = buffer.data(ki + 697);
    const auto *ki_698 = buffer.data(ki + 698);
    const auto *ki_699 = buffer.data(ki + 699);
    const auto *ki_700 = buffer.data(ki + 700);
    const auto *ki_701 = buffer.data(ki + 701);
    const auto *ki_702 = buffer.data(ki + 702);
    const auto *ki_703 = buffer.data(ki + 703);
    const auto *ki_704 = buffer.data(ki + 704);
    const auto *ki_705 = buffer.data(ki + 705);
    const auto *ki_706 = buffer.data(ki + 706);
    const auto *ki_707 = buffer.data(ki + 707);
    const auto *ki_708 = buffer.data(ki + 708);
    const auto *ki_709 = buffer.data(ki + 709);
    const auto *ki_710 = buffer.data(ki + 710);
    const auto *ki_711 = buffer.data(ki + 711);
    const auto *ki_712 = buffer.data(ki + 712);
    const auto *ki_713 = buffer.data(ki + 713);
    const auto *ki_714 = buffer.data(ki + 714);
    const auto *ki_715 = buffer.data(ki + 715);
    const auto *ki_716 = buffer.data(ki + 716);
    const auto *ki_717 = buffer.data(ki + 717);
    const auto *ki_718 = buffer.data(ki + 718);
    const auto *ki_719 = buffer.data(ki + 719);
    const auto *ki_720 = buffer.data(ki + 720);
    const auto *ki_721 = buffer.data(ki + 721);
    const auto *ki_722 = buffer.data(ki + 722);
    const auto *ki_723 = buffer.data(ki + 723);
    const auto *ki_724 = buffer.data(ki + 724);
    const auto *ki_725 = buffer.data(ki + 725);
    const auto *ki_726 = buffer.data(ki + 726);
    const auto *ki_727 = buffer.data(ki + 727);
    const auto *ki_728 = buffer.data(ki + 728);
    const auto *ki_729 = buffer.data(ki + 729);
    const auto *ki_730 = buffer.data(ki + 730);
    const auto *ki_731 = buffer.data(ki + 731);
    const auto *ki_732 = buffer.data(ki + 732);
    const auto *ki_733 = buffer.data(ki + 733);
    const auto *ki_734 = buffer.data(ki + 734);
    const auto *ki_735 = buffer.data(ki + 735);
    const auto *ki_736 = buffer.data(ki + 736);
    const auto *ki_737 = buffer.data(ki + 737);
    const auto *ki_738 = buffer.data(ki + 738);
    const auto *ki_739 = buffer.data(ki + 739);
    const auto *ki_740 = buffer.data(ki + 740);
    const auto *ki_741 = buffer.data(ki + 741);
    const auto *ki_742 = buffer.data(ki + 742);
    const auto *ki_743 = buffer.data(ki + 743);
    const auto *ki_744 = buffer.data(ki + 744);
    const auto *ki_745 = buffer.data(ki + 745);
    const auto *ki_746 = buffer.data(ki + 746);
    const auto *ki_747 = buffer.data(ki + 747);
    const auto *ki_748 = buffer.data(ki + 748);
    const auto *ki_749 = buffer.data(ki + 749);
    const auto *ki_750 = buffer.data(ki + 750);
    const auto *ki_751 = buffer.data(ki + 751);
    const auto *ki_752 = buffer.data(ki + 752);

    const auto *mi_1079 = buffer.data(mi + 1079);
    const auto *mi_1080 = buffer.data(mi + 1080);
    const auto *mi_1081 = buffer.data(mi + 1081);
    const auto *mi_1082 = buffer.data(mi + 1082);
    const auto *mi_1083 = buffer.data(mi + 1083);
    const auto *mi_1084 = buffer.data(mi + 1084);
    const auto *mi_1085 = buffer.data(mi + 1085);
    const auto *mi_1086 = buffer.data(mi + 1086);
    const auto *mi_1087 = buffer.data(mi + 1087);
    const auto *mi_1088 = buffer.data(mi + 1088);
    const auto *mi_1089 = buffer.data(mi + 1089);
    const auto *mi_1090 = buffer.data(mi + 1090);
    const auto *mi_1091 = buffer.data(mi + 1091);
    const auto *mi_1092 = buffer.data(mi + 1092);
    const auto *mi_1093 = buffer.data(mi + 1093);
    const auto *mi_1094 = buffer.data(mi + 1094);
    const auto *mi_1095 = buffer.data(mi + 1095);
    const auto *mi_1096 = buffer.data(mi + 1096);
    const auto *mi_1097 = buffer.data(mi + 1097);
    const auto *mi_1098 = buffer.data(mi + 1098);
    const auto *mi_1099 = buffer.data(mi + 1099);
    const auto *mi_1100 = buffer.data(mi + 1100);
    const auto *mi_1101 = buffer.data(mi + 1101);
    const auto *mi_1102 = buffer.data(mi + 1102);
    const auto *mi_1103 = buffer.data(mi + 1103);
    const auto *mi_1104 = buffer.data(mi + 1104);
    const auto *mi_1105 = buffer.data(mi + 1105);
    const auto *mi_1106 = buffer.data(mi + 1106);
    const auto *mi_1107 = buffer.data(mi + 1107);
    const auto *mi_1108 = buffer.data(mi + 1108);
    const auto *mi_1109 = buffer.data(mi + 1109);
    const auto *mi_1110 = buffer.data(mi + 1110);
    const auto *mi_1111 = buffer.data(mi + 1111);
    const auto *mi_1112 = buffer.data(mi + 1112);
    const auto *mi_1113 = buffer.data(mi + 1113);
    const auto *mi_1114 = buffer.data(mi + 1114);
    const auto *mi_1115 = buffer.data(mi + 1115);
    const auto *mi_1116 = buffer.data(mi + 1116);
    const auto *mi_1117 = buffer.data(mi + 1117);
    const auto *mi_1118 = buffer.data(mi + 1118);
    const auto *mi_1119 = buffer.data(mi + 1119);
    const auto *mi_1120 = buffer.data(mi + 1120);
    const auto *mi_1121 = buffer.data(mi + 1121);
    const auto *mi_1122 = buffer.data(mi + 1122);
    const auto *mi_1123 = buffer.data(mi + 1123);
    const auto *mi_1124 = buffer.data(mi + 1124);
    const auto *mi_1125 = buffer.data(mi + 1125);
    const auto *mi_1126 = buffer.data(mi + 1126);
    const auto *mi_1127 = buffer.data(mi + 1127);
    const auto *mi_1128 = buffer.data(mi + 1128);
    const auto *mi_1129 = buffer.data(mi + 1129);
    const auto *mi_1130 = buffer.data(mi + 1130);
    const auto *mi_1131 = buffer.data(mi + 1131);
    const auto *mi_1132 = buffer.data(mi + 1132);
    const auto *mi_1133 = buffer.data(mi + 1133);
    const auto *mi_1134 = buffer.data(mi + 1134);
    const auto *mi_1135 = buffer.data(mi + 1135);
    const auto *mi_1136 = buffer.data(mi + 1136);
    const auto *mi_1137 = buffer.data(mi + 1137);
    const auto *mi_1138 = buffer.data(mi + 1138);
    const auto *mi_1139 = buffer.data(mi + 1139);
    const auto *mi_1140 = buffer.data(mi + 1140);
    const auto *mi_1141 = buffer.data(mi + 1141);
    const auto *mi_1142 = buffer.data(mi + 1142);
    const auto *mi_1143 = buffer.data(mi + 1143);
    const auto *mi_1144 = buffer.data(mi + 1144);
    const auto *mi_1145 = buffer.data(mi + 1145);
    const auto *mi_1146 = buffer.data(mi + 1146);
    const auto *mi_1147 = buffer.data(mi + 1147);
    const auto *mi_1148 = buffer.data(mi + 1148);
    const auto *mi_1149 = buffer.data(mi + 1149);
    const auto *mi_1150 = buffer.data(mi + 1150);
    const auto *mi_1151 = buffer.data(mi + 1151);
    const auto *mi_1152 = buffer.data(mi + 1152);
    const auto *mi_1153 = buffer.data(mi + 1153);
    const auto *mi_1154 = buffer.data(mi + 1154);
    const auto *mi_1155 = buffer.data(mi + 1155);
    const auto *mi_1156 = buffer.data(mi + 1156);
    const auto *mi_1157 = buffer.data(mi + 1157);
    const auto *mi_1158 = buffer.data(mi + 1158);
    const auto *mi_1159 = buffer.data(mi + 1159);
    const auto *mi_1160 = buffer.data(mi + 1160);
    const auto *mi_1161 = buffer.data(mi + 1161);
    const auto *mi_1162 = buffer.data(mi + 1162);
    const auto *mi_1163 = buffer.data(mi + 1163);
    const auto *mi_1164 = buffer.data(mi + 1164);
    const auto *mi_1165 = buffer.data(mi + 1165);
    const auto *mi_1166 = buffer.data(mi + 1166);
    const auto *mi_1167 = buffer.data(mi + 1167);
    const auto *mi_1168 = buffer.data(mi + 1168);
    const auto *mi_1169 = buffer.data(mi + 1169);
    const auto *mi_1170 = buffer.data(mi + 1170);
    const auto *mi_1171 = buffer.data(mi + 1171);
    const auto *mi_1172 = buffer.data(mi + 1172);
    const auto *mi_1173 = buffer.data(mi + 1173);
    const auto *mi_1174 = buffer.data(mi + 1174);
    const auto *mi_1175 = buffer.data(mi + 1175);
    const auto *mi_1176 = buffer.data(mi + 1176);
    const auto *mi_1177 = buffer.data(mi + 1177);
    const auto *mi_1178 = buffer.data(mi + 1178);
    const auto *mi_1179 = buffer.data(mi + 1179);
    const auto *mi_1180 = buffer.data(mi + 1180);
    const auto *mi_1181 = buffer.data(mi + 1181);
    const auto *mi_1182 = buffer.data(mi + 1182);
    const auto *mi_1183 = buffer.data(mi + 1183);
    const auto *mi_1184 = buffer.data(mi + 1184);
    const auto *mi_1185 = buffer.data(mi + 1185);
    const auto *mi_1186 = buffer.data(mi + 1186);
    const auto *mi_1187 = buffer.data(mi + 1187);
    const auto *mi_1188 = buffer.data(mi + 1188);
    const auto *mi_1189 = buffer.data(mi + 1189);
    const auto *mi_1190 = buffer.data(mi + 1190);
    const auto *mi_1191 = buffer.data(mi + 1191);
    const auto *mi_1192 = buffer.data(mi + 1192);
    const auto *mi_1193 = buffer.data(mi + 1193);
    const auto *mi_1194 = buffer.data(mi + 1194);
    const auto *mi_1195 = buffer.data(mi + 1195);
    const auto *mi_1196 = buffer.data(mi + 1196);
    const auto *mi_1197 = buffer.data(mi + 1197);
    const auto *mi_1198 = buffer.data(mi + 1198);
    const auto *mi_1199 = buffer.data(mi + 1199);
    const auto *mi_1200 = buffer.data(mi + 1200);
    const auto *mi_1201 = buffer.data(mi + 1201);
    const auto *mi_1202 = buffer.data(mi + 1202);
    const auto *mi_1203 = buffer.data(mi + 1203);
    const auto *mi_1204 = buffer.data(mi + 1204);
    const auto *mi_1205 = buffer.data(mi + 1205);
    const auto *mi_1206 = buffer.data(mi + 1206);
    const auto *mi_1207 = buffer.data(mi + 1207);
    const auto *mi_1208 = buffer.data(mi + 1208);
    const auto *mi_1209 = buffer.data(mi + 1209);
    const auto *mi_1210 = buffer.data(mi + 1210);
    const auto *mi_1211 = buffer.data(mi + 1211);
    const auto *mi_1212 = buffer.data(mi + 1212);
    const auto *mi_1213 = buffer.data(mi + 1213);
    const auto *mi_1214 = buffer.data(mi + 1214);
    const auto *mi_1215 = buffer.data(mi + 1215);
    const auto *mi_1216 = buffer.data(mi + 1216);
    const auto *mi_1217 = buffer.data(mi + 1217);
    const auto *mi_1218 = buffer.data(mi + 1218);
    const auto *mi_1219 = buffer.data(mi + 1219);
    const auto *mi_1220 = buffer.data(mi + 1220);
    const auto *mi_1221 = buffer.data(mi + 1221);
    const auto *mi_1222 = buffer.data(mi + 1222);
    const auto *mi_1223 = buffer.data(mi + 1223);
    const auto *mi_1224 = buffer.data(mi + 1224);
    const auto *mi_1225 = buffer.data(mi + 1225);
    const auto *mi_1226 = buffer.data(mi + 1226);
    const auto *mi_1227 = buffer.data(mi + 1227);
    const auto *mi_1228 = buffer.data(mi + 1228);

#pragma omp simd aligned(t_827, t_828, t_829, t_830, t_831, ki_603, ki_604, ki_605, ki_606, \
                         ki_607, mi_1079, mi_1080, mi_1081, mi_1082, \
                         mi_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = -ki_603[k]
                   + f_0 * mi_1079[k];

        t_828[k] = -ki_604[k]
                   + f_0 * mi_1080[k];

        t_829[k] = -ki_605[k]
                   + f_0 * mi_1081[k];

        t_830[k] = -ki_606[k]
                   + f_0 * mi_1082[k];

        t_831[k] = -ki_607[k]
                   + f_0 * mi_1083[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, t_835, t_836, ki_608, ki_609, ki_610, ki_611, \
                         ki_612, mi_1084, mi_1085, mi_1086, mi_1087, \
                         mi_1088 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = -ki_608[k]
                   + f_0 * mi_1084[k];

        t_833[k] = -ki_609[k]
                   + f_0 * mi_1085[k];

        t_834[k] = -ki_610[k]
                   + f_0 * mi_1086[k];

        t_835[k] = -ki_611[k]
                   + f_0 * mi_1087[k];

        t_836[k] = -ki_612[k]
                   + f_0 * mi_1088[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, t_840, t_841, ki_613, ki_614, ki_615, ki_616, \
                         ki_617, mi_1089, mi_1090, mi_1091, mi_1092, \
                         mi_1093 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = -ki_613[k]
                   + f_0 * mi_1089[k];

        t_838[k] = -ki_614[k]
                   + f_0 * mi_1090[k];

        t_839[k] = -ki_615[k]
                   + f_0 * mi_1091[k];

        t_840[k] = -2.0 * ki_616[k]
                   + f_0 * mi_1092[k];

        t_841[k] = -2.0 * ki_617[k]
                   + f_0 * mi_1093[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, t_845, t_846, ki_618, ki_619, ki_620, ki_621, \
                         ki_622, mi_1094, mi_1095, mi_1096, mi_1097, \
                         mi_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = -2.0 * ki_618[k]
                   + f_0 * mi_1094[k];

        t_843[k] = -2.0 * ki_619[k]
                   + f_0 * mi_1095[k];

        t_844[k] = -2.0 * ki_620[k]
                   + f_0 * mi_1096[k];

        t_845[k] = -2.0 * ki_621[k]
                   + f_0 * mi_1097[k];

        t_846[k] = -2.0 * ki_622[k]
                   + f_0 * mi_1098[k];
    }

#pragma omp simd aligned(t_847, t_848, t_849, t_850, t_851, ki_623, ki_624, ki_625, ki_626, \
                         ki_627, mi_1099, mi_1100, mi_1101, mi_1102, \
                         mi_1103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_847[k] = -2.0 * ki_623[k]
                   + f_0 * mi_1099[k];

        t_848[k] = -2.0 * ki_624[k]
                   + f_0 * mi_1100[k];

        t_849[k] = -2.0 * ki_625[k]
                   + f_0 * mi_1101[k];

        t_850[k] = -2.0 * ki_626[k]
                   + f_0 * mi_1102[k];

        t_851[k] = -2.0 * ki_627[k]
                   + f_0 * mi_1103[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, t_856, ki_628, ki_629, ki_630, ki_631, \
                         ki_632, mi_1104, mi_1105, mi_1106, mi_1107, \
                         mi_1108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = -2.0 * ki_628[k]
                   + f_0 * mi_1104[k];

        t_853[k] = -2.0 * ki_629[k]
                   + f_0 * mi_1105[k];

        t_854[k] = -2.0 * ki_630[k]
                   + f_0 * mi_1106[k];

        t_855[k] = -2.0 * ki_631[k]
                   + f_0 * mi_1107[k];

        t_856[k] = -2.0 * ki_632[k]
                   + f_0 * mi_1108[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, t_860, t_861, ki_633, ki_634, ki_635, ki_636, \
                         ki_637, mi_1109, mi_1110, mi_1111, mi_1112, \
                         mi_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = -2.0 * ki_633[k]
                   + f_0 * mi_1109[k];

        t_858[k] = -2.0 * ki_634[k]
                   + f_0 * mi_1110[k];

        t_859[k] = -2.0 * ki_635[k]
                   + f_0 * mi_1111[k];

        t_860[k] = -2.0 * ki_636[k]
                   + f_0 * mi_1112[k];

        t_861[k] = -2.0 * ki_637[k]
                   + f_0 * mi_1113[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, t_866, ki_638, ki_639, ki_640, ki_641, \
                         ki_642, mi_1114, mi_1115, mi_1116, mi_1117, \
                         mi_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = -2.0 * ki_638[k]
                   + f_0 * mi_1114[k];

        t_863[k] = -2.0 * ki_639[k]
                   + f_0 * mi_1115[k];

        t_864[k] = -2.0 * ki_640[k]
                   + f_0 * mi_1116[k];

        t_865[k] = -2.0 * ki_641[k]
                   + f_0 * mi_1117[k];

        t_866[k] = -2.0 * ki_642[k]
                   + f_0 * mi_1118[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, t_870, t_871, ki_643, ki_644, ki_645, ki_646, \
                         ki_647, mi_1119, mi_1120, mi_1121, mi_1122, \
                         mi_1123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = -2.0 * ki_643[k]
                   + f_0 * mi_1119[k];

        t_868[k] = -3.0 * ki_644[k]
                   + f_0 * mi_1120[k];

        t_869[k] = -3.0 * ki_645[k]
                   + f_0 * mi_1121[k];

        t_870[k] = -3.0 * ki_646[k]
                   + f_0 * mi_1122[k];

        t_871[k] = -3.0 * ki_647[k]
                   + f_0 * mi_1123[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, t_875, t_876, ki_648, ki_649, ki_650, ki_651, \
                         ki_652, mi_1124, mi_1125, mi_1126, mi_1127, \
                         mi_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = -3.0 * ki_648[k]
                   + f_0 * mi_1124[k];

        t_873[k] = -3.0 * ki_649[k]
                   + f_0 * mi_1125[k];

        t_874[k] = -3.0 * ki_650[k]
                   + f_0 * mi_1126[k];

        t_875[k] = -3.0 * ki_651[k]
                   + f_0 * mi_1127[k];

        t_876[k] = -3.0 * ki_652[k]
                   + f_0 * mi_1128[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, t_880, t_881, ki_653, ki_654, ki_655, ki_656, \
                         ki_657, mi_1129, mi_1130, mi_1131, mi_1132, \
                         mi_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = -3.0 * ki_653[k]
                   + f_0 * mi_1129[k];

        t_878[k] = -3.0 * ki_654[k]
                   + f_0 * mi_1130[k];

        t_879[k] = -3.0 * ki_655[k]
                   + f_0 * mi_1131[k];

        t_880[k] = -3.0 * ki_656[k]
                   + f_0 * mi_1132[k];

        t_881[k] = -3.0 * ki_657[k]
                   + f_0 * mi_1133[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, t_886, ki_658, ki_659, ki_660, ki_661, \
                         ki_662, mi_1134, mi_1135, mi_1136, mi_1137, \
                         mi_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = -3.0 * ki_658[k]
                   + f_0 * mi_1134[k];

        t_883[k] = -3.0 * ki_659[k]
                   + f_0 * mi_1135[k];

        t_884[k] = -3.0 * ki_660[k]
                   + f_0 * mi_1136[k];

        t_885[k] = -3.0 * ki_661[k]
                   + f_0 * mi_1137[k];

        t_886[k] = -3.0 * ki_662[k]
                   + f_0 * mi_1138[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, t_890, t_891, ki_663, ki_664, ki_665, ki_666, \
                         ki_667, mi_1139, mi_1140, mi_1141, mi_1142, \
                         mi_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = -3.0 * ki_663[k]
                   + f_0 * mi_1139[k];

        t_888[k] = -3.0 * ki_664[k]
                   + f_0 * mi_1140[k];

        t_889[k] = -3.0 * ki_665[k]
                   + f_0 * mi_1141[k];

        t_890[k] = -3.0 * ki_666[k]
                   + f_0 * mi_1142[k];

        t_891[k] = -3.0 * ki_667[k]
                   + f_0 * mi_1143[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, t_895, t_896, ki_668, ki_669, ki_670, ki_671, \
                         ki_672, mi_1144, mi_1145, mi_1146, mi_1147, \
                         mi_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = -3.0 * ki_668[k]
                   + f_0 * mi_1144[k];

        t_893[k] = -3.0 * ki_669[k]
                   + f_0 * mi_1145[k];

        t_894[k] = -3.0 * ki_670[k]
                   + f_0 * mi_1146[k];

        t_895[k] = -3.0 * ki_671[k]
                   + f_0 * mi_1147[k];

        t_896[k] = -4.0 * ki_672[k]
                   + f_0 * mi_1148[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, t_900, t_901, ki_673, ki_674, ki_675, ki_676, \
                         ki_677, mi_1149, mi_1150, mi_1151, mi_1152, \
                         mi_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = -4.0 * ki_673[k]
                   + f_0 * mi_1149[k];

        t_898[k] = -4.0 * ki_674[k]
                   + f_0 * mi_1150[k];

        t_899[k] = -4.0 * ki_675[k]
                   + f_0 * mi_1151[k];

        t_900[k] = -4.0 * ki_676[k]
                   + f_0 * mi_1152[k];

        t_901[k] = -4.0 * ki_677[k]
                   + f_0 * mi_1153[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, t_905, t_906, ki_678, ki_679, ki_680, ki_681, \
                         ki_682, mi_1154, mi_1155, mi_1156, mi_1157, \
                         mi_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = -4.0 * ki_678[k]
                   + f_0 * mi_1154[k];

        t_903[k] = -4.0 * ki_679[k]
                   + f_0 * mi_1155[k];

        t_904[k] = -4.0 * ki_680[k]
                   + f_0 * mi_1156[k];

        t_905[k] = -4.0 * ki_681[k]
                   + f_0 * mi_1157[k];

        t_906[k] = -4.0 * ki_682[k]
                   + f_0 * mi_1158[k];
    }

#pragma omp simd aligned(t_907, t_908, t_909, t_910, t_911, ki_683, ki_684, ki_685, ki_686, \
                         ki_687, mi_1159, mi_1160, mi_1161, mi_1162, \
                         mi_1163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_907[k] = -4.0 * ki_683[k]
                   + f_0 * mi_1159[k];

        t_908[k] = -4.0 * ki_684[k]
                   + f_0 * mi_1160[k];

        t_909[k] = -4.0 * ki_685[k]
                   + f_0 * mi_1161[k];

        t_910[k] = -4.0 * ki_686[k]
                   + f_0 * mi_1162[k];

        t_911[k] = -4.0 * ki_687[k]
                   + f_0 * mi_1163[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, t_915, t_916, ki_688, ki_689, ki_690, ki_691, \
                         ki_692, mi_1164, mi_1165, mi_1166, mi_1167, \
                         mi_1168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = -4.0 * ki_688[k]
                   + f_0 * mi_1164[k];

        t_913[k] = -4.0 * ki_689[k]
                   + f_0 * mi_1165[k];

        t_914[k] = -4.0 * ki_690[k]
                   + f_0 * mi_1166[k];

        t_915[k] = -4.0 * ki_691[k]
                   + f_0 * mi_1167[k];

        t_916[k] = -4.0 * ki_692[k]
                   + f_0 * mi_1168[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, t_920, t_921, ki_693, ki_694, ki_695, ki_696, \
                         ki_697, mi_1169, mi_1170, mi_1171, mi_1172, \
                         mi_1173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = -4.0 * ki_693[k]
                   + f_0 * mi_1169[k];

        t_918[k] = -4.0 * ki_694[k]
                   + f_0 * mi_1170[k];

        t_919[k] = -4.0 * ki_695[k]
                   + f_0 * mi_1171[k];

        t_920[k] = -4.0 * ki_696[k]
                   + f_0 * mi_1172[k];

        t_921[k] = -4.0 * ki_697[k]
                   + f_0 * mi_1173[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, t_926, ki_698, ki_699, ki_700, ki_701, \
                         ki_702, mi_1174, mi_1175, mi_1176, mi_1177, \
                         mi_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = -4.0 * ki_698[k]
                   + f_0 * mi_1174[k];

        t_923[k] = -4.0 * ki_699[k]
                   + f_0 * mi_1175[k];

        t_924[k] = -5.0 * ki_700[k]
                   + f_0 * mi_1176[k];

        t_925[k] = -5.0 * ki_701[k]
                   + f_0 * mi_1177[k];

        t_926[k] = -5.0 * ki_702[k]
                   + f_0 * mi_1178[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, t_931, ki_703, ki_704, ki_705, ki_706, \
                         ki_707, mi_1179, mi_1180, mi_1181, mi_1182, \
                         mi_1183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = -5.0 * ki_703[k]
                   + f_0 * mi_1179[k];

        t_928[k] = -5.0 * ki_704[k]
                   + f_0 * mi_1180[k];

        t_929[k] = -5.0 * ki_705[k]
                   + f_0 * mi_1181[k];

        t_930[k] = -5.0 * ki_706[k]
                   + f_0 * mi_1182[k];

        t_931[k] = -5.0 * ki_707[k]
                   + f_0 * mi_1183[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, t_935, t_936, ki_708, ki_709, ki_710, ki_711, \
                         ki_712, mi_1184, mi_1185, mi_1186, mi_1187, \
                         mi_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = -5.0 * ki_708[k]
                   + f_0 * mi_1184[k];

        t_933[k] = -5.0 * ki_709[k]
                   + f_0 * mi_1185[k];

        t_934[k] = -5.0 * ki_710[k]
                   + f_0 * mi_1186[k];

        t_935[k] = -5.0 * ki_711[k]
                   + f_0 * mi_1187[k];

        t_936[k] = -5.0 * ki_712[k]
                   + f_0 * mi_1188[k];
    }

#pragma omp simd aligned(t_937, t_938, t_939, t_940, t_941, ki_713, ki_714, ki_715, ki_716, \
                         ki_717, mi_1189, mi_1190, mi_1191, mi_1192, \
                         mi_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_937[k] = -5.0 * ki_713[k]
                   + f_0 * mi_1189[k];

        t_938[k] = -5.0 * ki_714[k]
                   + f_0 * mi_1190[k];

        t_939[k] = -5.0 * ki_715[k]
                   + f_0 * mi_1191[k];

        t_940[k] = -5.0 * ki_716[k]
                   + f_0 * mi_1192[k];

        t_941[k] = -5.0 * ki_717[k]
                   + f_0 * mi_1193[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, t_946, ki_718, ki_719, ki_720, ki_721, \
                         ki_722, mi_1194, mi_1195, mi_1196, mi_1197, \
                         mi_1198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = -5.0 * ki_718[k]
                   + f_0 * mi_1194[k];

        t_943[k] = -5.0 * ki_719[k]
                   + f_0 * mi_1195[k];

        t_944[k] = -5.0 * ki_720[k]
                   + f_0 * mi_1196[k];

        t_945[k] = -5.0 * ki_721[k]
                   + f_0 * mi_1197[k];

        t_946[k] = -5.0 * ki_722[k]
                   + f_0 * mi_1198[k];
    }

#pragma omp simd aligned(t_947, t_948, t_949, t_950, t_951, ki_723, ki_724, ki_725, ki_726, \
                         ki_727, mi_1199, mi_1200, mi_1201, mi_1202, \
                         mi_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = -5.0 * ki_723[k]
                   + f_0 * mi_1199[k];

        t_948[k] = -5.0 * ki_724[k]
                   + f_0 * mi_1200[k];

        t_949[k] = -5.0 * ki_725[k]
                   + f_0 * mi_1201[k];

        t_950[k] = -5.0 * ki_726[k]
                   + f_0 * mi_1202[k];

        t_951[k] = -5.0 * ki_727[k]
                   + f_0 * mi_1203[k];
    }

#pragma omp simd aligned(t_952, t_953, t_954, t_955, t_956, ki_728, ki_729, ki_730, ki_731, \
                         ki_732, mi_1204, mi_1205, mi_1206, mi_1207, \
                         mi_1208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_952[k] = -6.0 * ki_728[k]
                   + f_0 * mi_1204[k];

        t_953[k] = -6.0 * ki_729[k]
                   + f_0 * mi_1205[k];

        t_954[k] = -6.0 * ki_730[k]
                   + f_0 * mi_1206[k];

        t_955[k] = -6.0 * ki_731[k]
                   + f_0 * mi_1207[k];

        t_956[k] = -6.0 * ki_732[k]
                   + f_0 * mi_1208[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, t_961, ki_733, ki_734, ki_735, ki_736, \
                         ki_737, mi_1209, mi_1210, mi_1211, mi_1212, \
                         mi_1213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = -6.0 * ki_733[k]
                   + f_0 * mi_1209[k];

        t_958[k] = -6.0 * ki_734[k]
                   + f_0 * mi_1210[k];

        t_959[k] = -6.0 * ki_735[k]
                   + f_0 * mi_1211[k];

        t_960[k] = -6.0 * ki_736[k]
                   + f_0 * mi_1212[k];

        t_961[k] = -6.0 * ki_737[k]
                   + f_0 * mi_1213[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, t_966, ki_738, ki_739, ki_740, ki_741, \
                         ki_742, mi_1214, mi_1215, mi_1216, mi_1217, \
                         mi_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = -6.0 * ki_738[k]
                   + f_0 * mi_1214[k];

        t_963[k] = -6.0 * ki_739[k]
                   + f_0 * mi_1215[k];

        t_964[k] = -6.0 * ki_740[k]
                   + f_0 * mi_1216[k];

        t_965[k] = -6.0 * ki_741[k]
                   + f_0 * mi_1217[k];

        t_966[k] = -6.0 * ki_742[k]
                   + f_0 * mi_1218[k];
    }

#pragma omp simd aligned(t_967, t_968, t_969, t_970, t_971, ki_743, ki_744, ki_745, ki_746, \
                         ki_747, mi_1219, mi_1220, mi_1221, mi_1222, \
                         mi_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_967[k] = -6.0 * ki_743[k]
                   + f_0 * mi_1219[k];

        t_968[k] = -6.0 * ki_744[k]
                   + f_0 * mi_1220[k];

        t_969[k] = -6.0 * ki_745[k]
                   + f_0 * mi_1221[k];

        t_970[k] = -6.0 * ki_746[k]
                   + f_0 * mi_1222[k];

        t_971[k] = -6.0 * ki_747[k]
                   + f_0 * mi_1223[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, t_976, ki_748, ki_749, ki_750, ki_751, \
                         ki_752, mi_1224, mi_1225, mi_1226, mi_1227, \
                         mi_1228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = -6.0 * ki_748[k]
                   + f_0 * mi_1224[k];

        t_973[k] = -6.0 * ki_749[k]
                   + f_0 * mi_1225[k];

        t_974[k] = -6.0 * ki_750[k]
                   + f_0 * mi_1226[k];

        t_975[k] = -6.0 * ki_751[k]
                   + f_0 * mi_1227[k];

        t_976[k] = -6.0 * ki_752[k]
                   + f_0 * mi_1228[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_2_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);
    auto *t_983 = buffer.data(target + 983);
    auto *t_984 = buffer.data(target + 984);
    auto *t_985 = buffer.data(target + 985);
    auto *t_986 = buffer.data(target + 986);
    auto *t_987 = buffer.data(target + 987);
    auto *t_988 = buffer.data(target + 988);
    auto *t_989 = buffer.data(target + 989);
    auto *t_990 = buffer.data(target + 990);
    auto *t_991 = buffer.data(target + 991);
    auto *t_992 = buffer.data(target + 992);
    auto *t_993 = buffer.data(target + 993);
    auto *t_994 = buffer.data(target + 994);
    auto *t_995 = buffer.data(target + 995);
    auto *t_996 = buffer.data(target + 996);
    auto *t_997 = buffer.data(target + 997);
    auto *t_998 = buffer.data(target + 998);
    auto *t_999 = buffer.data(target + 999);
    auto *t_1000 = buffer.data(target + 1000);
    auto *t_1001 = buffer.data(target + 1001);
    auto *t_1002 = buffer.data(target + 1002);
    auto *t_1003 = buffer.data(target + 1003);
    auto *t_1004 = buffer.data(target + 1004);
    auto *t_1005 = buffer.data(target + 1005);
    auto *t_1006 = buffer.data(target + 1006);
    auto *t_1007 = buffer.data(target + 1007);
    auto *t_1008 = buffer.data(target + 1008);
    auto *t_1009 = buffer.data(target + 1009);
    auto *t_1010 = buffer.data(target + 1010);
    auto *t_1011 = buffer.data(target + 1011);
    auto *t_1012 = buffer.data(target + 1012);
    auto *t_1013 = buffer.data(target + 1013);
    auto *t_1014 = buffer.data(target + 1014);
    auto *t_1015 = buffer.data(target + 1015);
    auto *t_1016 = buffer.data(target + 1016);
    auto *t_1017 = buffer.data(target + 1017);
    auto *t_1018 = buffer.data(target + 1018);
    auto *t_1019 = buffer.data(target + 1019);
    auto *t_1020 = buffer.data(target + 1020);
    auto *t_1021 = buffer.data(target + 1021);
    auto *t_1022 = buffer.data(target + 1022);
    auto *t_1023 = buffer.data(target + 1023);
    auto *t_1024 = buffer.data(target + 1024);
    auto *t_1025 = buffer.data(target + 1025);
    auto *t_1026 = buffer.data(target + 1026);
    auto *t_1027 = buffer.data(target + 1027);
    auto *t_1028 = buffer.data(target + 1028);
    auto *t_1029 = buffer.data(target + 1029);
    auto *t_1030 = buffer.data(target + 1030);
    auto *t_1031 = buffer.data(target + 1031);
    auto *t_1032 = buffer.data(target + 1032);
    auto *t_1033 = buffer.data(target + 1033);
    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);
    auto *t_1050 = buffer.data(target + 1050);
    auto *t_1051 = buffer.data(target + 1051);
    auto *t_1052 = buffer.data(target + 1052);
    auto *t_1053 = buffer.data(target + 1053);
    auto *t_1054 = buffer.data(target + 1054);
    auto *t_1055 = buffer.data(target + 1055);
    auto *t_1056 = buffer.data(target + 1056);
    auto *t_1057 = buffer.data(target + 1057);
    auto *t_1058 = buffer.data(target + 1058);
    auto *t_1059 = buffer.data(target + 1059);
    auto *t_1060 = buffer.data(target + 1060);
    auto *t_1061 = buffer.data(target + 1061);
    auto *t_1062 = buffer.data(target + 1062);
    auto *t_1063 = buffer.data(target + 1063);
    auto *t_1064 = buffer.data(target + 1064);
    auto *t_1065 = buffer.data(target + 1065);
    auto *t_1066 = buffer.data(target + 1066);
    auto *t_1067 = buffer.data(target + 1067);
    auto *t_1068 = buffer.data(target + 1068);
    auto *t_1069 = buffer.data(target + 1069);
    auto *t_1070 = buffer.data(target + 1070);
    auto *t_1071 = buffer.data(target + 1071);
    auto *t_1072 = buffer.data(target + 1072);
    auto *t_1073 = buffer.data(target + 1073);
    auto *t_1074 = buffer.data(target + 1074);
    auto *t_1075 = buffer.data(target + 1075);
    auto *t_1076 = buffer.data(target + 1076);
    auto *t_1077 = buffer.data(target + 1077);
    auto *t_1078 = buffer.data(target + 1078);
    auto *t_1079 = buffer.data(target + 1079);
    auto *t_1080 = buffer.data(target + 1080);
    auto *t_1081 = buffer.data(target + 1081);
    auto *t_1082 = buffer.data(target + 1082);
    auto *t_1083 = buffer.data(target + 1083);
    auto *t_1084 = buffer.data(target + 1084);
    auto *t_1085 = buffer.data(target + 1085);
    auto *t_1086 = buffer.data(target + 1086);
    auto *t_1087 = buffer.data(target + 1087);
    auto *t_1088 = buffer.data(target + 1088);
    auto *t_1089 = buffer.data(target + 1089);
    auto *t_1090 = buffer.data(target + 1090);
    auto *t_1091 = buffer.data(target + 1091);
    auto *t_1092 = buffer.data(target + 1092);
    auto *t_1093 = buffer.data(target + 1093);
    auto *t_1094 = buffer.data(target + 1094);
    auto *t_1095 = buffer.data(target + 1095);
    auto *t_1096 = buffer.data(target + 1096);
    auto *t_1097 = buffer.data(target + 1097);
    auto *t_1098 = buffer.data(target + 1098);
    auto *t_1099 = buffer.data(target + 1099);
    auto *t_1100 = buffer.data(target + 1100);
    auto *t_1101 = buffer.data(target + 1101);
    auto *t_1102 = buffer.data(target + 1102);
    auto *t_1103 = buffer.data(target + 1103);
    auto *t_1104 = buffer.data(target + 1104);
    auto *t_1105 = buffer.data(target + 1105);
    auto *t_1106 = buffer.data(target + 1106);
    auto *t_1107 = buffer.data(target + 1107);
    auto *t_1108 = buffer.data(target + 1108);
    auto *t_1109 = buffer.data(target + 1109);
    auto *t_1110 = buffer.data(target + 1110);
    auto *t_1111 = buffer.data(target + 1111);
    auto *t_1112 = buffer.data(target + 1112);
    auto *t_1113 = buffer.data(target + 1113);
    auto *t_1114 = buffer.data(target + 1114);
    auto *t_1115 = buffer.data(target + 1115);
    auto *t_1116 = buffer.data(target + 1116);
    auto *t_1117 = buffer.data(target + 1117);
    auto *t_1118 = buffer.data(target + 1118);
    auto *t_1119 = buffer.data(target + 1119);
    auto *t_1120 = buffer.data(target + 1120);
    auto *t_1121 = buffer.data(target + 1121);
    auto *t_1122 = buffer.data(target + 1122);
    auto *t_1123 = buffer.data(target + 1123);
    auto *t_1124 = buffer.data(target + 1124);
    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);
    auto *t_1133 = buffer.data(target + 1133);
    auto *t_1134 = buffer.data(target + 1134);
    auto *t_1135 = buffer.data(target + 1135);
    auto *t_1136 = buffer.data(target + 1136);

    const auto *ki_753 = buffer.data(ki + 753);
    const auto *ki_754 = buffer.data(ki + 754);
    const auto *ki_755 = buffer.data(ki + 755);
    const auto *ki_756 = buffer.data(ki + 756);
    const auto *ki_757 = buffer.data(ki + 757);
    const auto *ki_758 = buffer.data(ki + 758);
    const auto *ki_759 = buffer.data(ki + 759);
    const auto *ki_760 = buffer.data(ki + 760);
    const auto *ki_761 = buffer.data(ki + 761);
    const auto *ki_762 = buffer.data(ki + 762);
    const auto *ki_763 = buffer.data(ki + 763);
    const auto *ki_764 = buffer.data(ki + 764);
    const auto *ki_765 = buffer.data(ki + 765);
    const auto *ki_766 = buffer.data(ki + 766);
    const auto *ki_767 = buffer.data(ki + 767);
    const auto *ki_768 = buffer.data(ki + 768);
    const auto *ki_769 = buffer.data(ki + 769);
    const auto *ki_770 = buffer.data(ki + 770);
    const auto *ki_771 = buffer.data(ki + 771);
    const auto *ki_772 = buffer.data(ki + 772);
    const auto *ki_773 = buffer.data(ki + 773);
    const auto *ki_774 = buffer.data(ki + 774);
    const auto *ki_775 = buffer.data(ki + 775);
    const auto *ki_776 = buffer.data(ki + 776);
    const auto *ki_777 = buffer.data(ki + 777);
    const auto *ki_778 = buffer.data(ki + 778);
    const auto *ki_779 = buffer.data(ki + 779);
    const auto *ki_780 = buffer.data(ki + 780);
    const auto *ki_781 = buffer.data(ki + 781);
    const auto *ki_782 = buffer.data(ki + 782);
    const auto *ki_783 = buffer.data(ki + 783);
    const auto *ki_784 = buffer.data(ki + 784);
    const auto *ki_785 = buffer.data(ki + 785);
    const auto *ki_786 = buffer.data(ki + 786);
    const auto *ki_787 = buffer.data(ki + 787);
    const auto *ki_788 = buffer.data(ki + 788);
    const auto *ki_789 = buffer.data(ki + 789);
    const auto *ki_790 = buffer.data(ki + 790);
    const auto *ki_791 = buffer.data(ki + 791);
    const auto *ki_792 = buffer.data(ki + 792);
    const auto *ki_793 = buffer.data(ki + 793);
    const auto *ki_794 = buffer.data(ki + 794);
    const auto *ki_795 = buffer.data(ki + 795);
    const auto *ki_796 = buffer.data(ki + 796);
    const auto *ki_797 = buffer.data(ki + 797);
    const auto *ki_798 = buffer.data(ki + 798);
    const auto *ki_799 = buffer.data(ki + 799);
    const auto *ki_800 = buffer.data(ki + 800);
    const auto *ki_801 = buffer.data(ki + 801);
    const auto *ki_802 = buffer.data(ki + 802);
    const auto *ki_803 = buffer.data(ki + 803);
    const auto *ki_804 = buffer.data(ki + 804);
    const auto *ki_805 = buffer.data(ki + 805);
    const auto *ki_806 = buffer.data(ki + 806);
    const auto *ki_807 = buffer.data(ki + 807);
    const auto *ki_808 = buffer.data(ki + 808);
    const auto *ki_809 = buffer.data(ki + 809);
    const auto *ki_810 = buffer.data(ki + 810);
    const auto *ki_811 = buffer.data(ki + 811);
    const auto *ki_812 = buffer.data(ki + 812);
    const auto *ki_813 = buffer.data(ki + 813);
    const auto *ki_814 = buffer.data(ki + 814);
    const auto *ki_815 = buffer.data(ki + 815);
    const auto *ki_816 = buffer.data(ki + 816);
    const auto *ki_817 = buffer.data(ki + 817);
    const auto *ki_818 = buffer.data(ki + 818);
    const auto *ki_819 = buffer.data(ki + 819);
    const auto *ki_820 = buffer.data(ki + 820);
    const auto *ki_821 = buffer.data(ki + 821);
    const auto *ki_822 = buffer.data(ki + 822);
    const auto *ki_823 = buffer.data(ki + 823);
    const auto *ki_824 = buffer.data(ki + 824);
    const auto *ki_825 = buffer.data(ki + 825);
    const auto *ki_826 = buffer.data(ki + 826);
    const auto *ki_827 = buffer.data(ki + 827);
    const auto *ki_828 = buffer.data(ki + 828);
    const auto *ki_829 = buffer.data(ki + 829);
    const auto *ki_830 = buffer.data(ki + 830);
    const auto *ki_831 = buffer.data(ki + 831);
    const auto *ki_832 = buffer.data(ki + 832);
    const auto *ki_833 = buffer.data(ki + 833);
    const auto *ki_834 = buffer.data(ki + 834);
    const auto *ki_835 = buffer.data(ki + 835);
    const auto *ki_836 = buffer.data(ki + 836);
    const auto *ki_837 = buffer.data(ki + 837);
    const auto *ki_838 = buffer.data(ki + 838);
    const auto *ki_839 = buffer.data(ki + 839);
    const auto *ki_840 = buffer.data(ki + 840);
    const auto *ki_841 = buffer.data(ki + 841);
    const auto *ki_842 = buffer.data(ki + 842);
    const auto *ki_843 = buffer.data(ki + 843);
    const auto *ki_844 = buffer.data(ki + 844);
    const auto *ki_845 = buffer.data(ki + 845);
    const auto *ki_846 = buffer.data(ki + 846);
    const auto *ki_847 = buffer.data(ki + 847);
    const auto *ki_848 = buffer.data(ki + 848);
    const auto *ki_849 = buffer.data(ki + 849);
    const auto *ki_850 = buffer.data(ki + 850);
    const auto *ki_851 = buffer.data(ki + 851);
    const auto *ki_852 = buffer.data(ki + 852);
    const auto *ki_853 = buffer.data(ki + 853);
    const auto *ki_854 = buffer.data(ki + 854);
    const auto *ki_855 = buffer.data(ki + 855);
    const auto *ki_856 = buffer.data(ki + 856);
    const auto *ki_857 = buffer.data(ki + 857);
    const auto *ki_858 = buffer.data(ki + 858);
    const auto *ki_859 = buffer.data(ki + 859);
    const auto *ki_860 = buffer.data(ki + 860);
    const auto *ki_861 = buffer.data(ki + 861);
    const auto *ki_862 = buffer.data(ki + 862);
    const auto *ki_863 = buffer.data(ki + 863);
    const auto *ki_864 = buffer.data(ki + 864);
    const auto *ki_865 = buffer.data(ki + 865);
    const auto *ki_866 = buffer.data(ki + 866);
    const auto *ki_867 = buffer.data(ki + 867);
    const auto *ki_868 = buffer.data(ki + 868);
    const auto *ki_869 = buffer.data(ki + 869);
    const auto *ki_870 = buffer.data(ki + 870);
    const auto *ki_871 = buffer.data(ki + 871);
    const auto *ki_872 = buffer.data(ki + 872);
    const auto *ki_873 = buffer.data(ki + 873);
    const auto *ki_874 = buffer.data(ki + 874);
    const auto *ki_875 = buffer.data(ki + 875);
    const auto *ki_876 = buffer.data(ki + 876);
    const auto *ki_877 = buffer.data(ki + 877);
    const auto *ki_878 = buffer.data(ki + 878);
    const auto *ki_879 = buffer.data(ki + 879);
    const auto *ki_880 = buffer.data(ki + 880);
    const auto *ki_881 = buffer.data(ki + 881);
    const auto *ki_882 = buffer.data(ki + 882);
    const auto *ki_883 = buffer.data(ki + 883);
    const auto *ki_884 = buffer.data(ki + 884);

    const auto *mi_1229 = buffer.data(mi + 1229);
    const auto *mi_1230 = buffer.data(mi + 1230);
    const auto *mi_1231 = buffer.data(mi + 1231);
    const auto *mi_1232 = buffer.data(mi + 1232);
    const auto *mi_1233 = buffer.data(mi + 1233);
    const auto *mi_1234 = buffer.data(mi + 1234);
    const auto *mi_1235 = buffer.data(mi + 1235);
    const auto *mi_1236 = buffer.data(mi + 1236);
    const auto *mi_1237 = buffer.data(mi + 1237);
    const auto *mi_1238 = buffer.data(mi + 1238);
    const auto *mi_1239 = buffer.data(mi + 1239);
    const auto *mi_1240 = buffer.data(mi + 1240);
    const auto *mi_1241 = buffer.data(mi + 1241);
    const auto *mi_1242 = buffer.data(mi + 1242);
    const auto *mi_1243 = buffer.data(mi + 1243);
    const auto *mi_1244 = buffer.data(mi + 1244);
    const auto *mi_1245 = buffer.data(mi + 1245);
    const auto *mi_1246 = buffer.data(mi + 1246);
    const auto *mi_1247 = buffer.data(mi + 1247);
    const auto *mi_1248 = buffer.data(mi + 1248);
    const auto *mi_1249 = buffer.data(mi + 1249);
    const auto *mi_1250 = buffer.data(mi + 1250);
    const auto *mi_1251 = buffer.data(mi + 1251);
    const auto *mi_1252 = buffer.data(mi + 1252);
    const auto *mi_1253 = buffer.data(mi + 1253);
    const auto *mi_1254 = buffer.data(mi + 1254);
    const auto *mi_1255 = buffer.data(mi + 1255);
    const auto *mi_1256 = buffer.data(mi + 1256);
    const auto *mi_1257 = buffer.data(mi + 1257);
    const auto *mi_1258 = buffer.data(mi + 1258);
    const auto *mi_1259 = buffer.data(mi + 1259);
    const auto *mi_1288 = buffer.data(mi + 1288);
    const auto *mi_1289 = buffer.data(mi + 1289);
    const auto *mi_1290 = buffer.data(mi + 1290);
    const auto *mi_1291 = buffer.data(mi + 1291);
    const auto *mi_1292 = buffer.data(mi + 1292);
    const auto *mi_1293 = buffer.data(mi + 1293);
    const auto *mi_1294 = buffer.data(mi + 1294);
    const auto *mi_1295 = buffer.data(mi + 1295);
    const auto *mi_1296 = buffer.data(mi + 1296);
    const auto *mi_1297 = buffer.data(mi + 1297);
    const auto *mi_1298 = buffer.data(mi + 1298);
    const auto *mi_1299 = buffer.data(mi + 1299);
    const auto *mi_1300 = buffer.data(mi + 1300);
    const auto *mi_1301 = buffer.data(mi + 1301);
    const auto *mi_1302 = buffer.data(mi + 1302);
    const auto *mi_1303 = buffer.data(mi + 1303);
    const auto *mi_1304 = buffer.data(mi + 1304);
    const auto *mi_1305 = buffer.data(mi + 1305);
    const auto *mi_1306 = buffer.data(mi + 1306);
    const auto *mi_1307 = buffer.data(mi + 1307);
    const auto *mi_1308 = buffer.data(mi + 1308);
    const auto *mi_1309 = buffer.data(mi + 1309);
    const auto *mi_1310 = buffer.data(mi + 1310);
    const auto *mi_1311 = buffer.data(mi + 1311);
    const auto *mi_1312 = buffer.data(mi + 1312);
    const auto *mi_1313 = buffer.data(mi + 1313);
    const auto *mi_1314 = buffer.data(mi + 1314);
    const auto *mi_1315 = buffer.data(mi + 1315);
    const auto *mi_1316 = buffer.data(mi + 1316);
    const auto *mi_1317 = buffer.data(mi + 1317);
    const auto *mi_1318 = buffer.data(mi + 1318);
    const auto *mi_1319 = buffer.data(mi + 1319);
    const auto *mi_1320 = buffer.data(mi + 1320);
    const auto *mi_1321 = buffer.data(mi + 1321);
    const auto *mi_1322 = buffer.data(mi + 1322);
    const auto *mi_1323 = buffer.data(mi + 1323);
    const auto *mi_1324 = buffer.data(mi + 1324);
    const auto *mi_1325 = buffer.data(mi + 1325);
    const auto *mi_1326 = buffer.data(mi + 1326);
    const auto *mi_1327 = buffer.data(mi + 1327);
    const auto *mi_1328 = buffer.data(mi + 1328);
    const auto *mi_1329 = buffer.data(mi + 1329);
    const auto *mi_1330 = buffer.data(mi + 1330);
    const auto *mi_1331 = buffer.data(mi + 1331);
    const auto *mi_1332 = buffer.data(mi + 1332);
    const auto *mi_1333 = buffer.data(mi + 1333);
    const auto *mi_1334 = buffer.data(mi + 1334);
    const auto *mi_1335 = buffer.data(mi + 1335);
    const auto *mi_1336 = buffer.data(mi + 1336);
    const auto *mi_1337 = buffer.data(mi + 1337);
    const auto *mi_1338 = buffer.data(mi + 1338);
    const auto *mi_1339 = buffer.data(mi + 1339);
    const auto *mi_1340 = buffer.data(mi + 1340);
    const auto *mi_1341 = buffer.data(mi + 1341);
    const auto *mi_1342 = buffer.data(mi + 1342);
    const auto *mi_1343 = buffer.data(mi + 1343);
    const auto *mi_1344 = buffer.data(mi + 1344);
    const auto *mi_1345 = buffer.data(mi + 1345);
    const auto *mi_1346 = buffer.data(mi + 1346);
    const auto *mi_1347 = buffer.data(mi + 1347);
    const auto *mi_1348 = buffer.data(mi + 1348);
    const auto *mi_1349 = buffer.data(mi + 1349);
    const auto *mi_1350 = buffer.data(mi + 1350);
    const auto *mi_1351 = buffer.data(mi + 1351);
    const auto *mi_1352 = buffer.data(mi + 1352);
    const auto *mi_1353 = buffer.data(mi + 1353);
    const auto *mi_1354 = buffer.data(mi + 1354);
    const auto *mi_1355 = buffer.data(mi + 1355);
    const auto *mi_1356 = buffer.data(mi + 1356);
    const auto *mi_1357 = buffer.data(mi + 1357);
    const auto *mi_1358 = buffer.data(mi + 1358);
    const auto *mi_1359 = buffer.data(mi + 1359);
    const auto *mi_1360 = buffer.data(mi + 1360);
    const auto *mi_1361 = buffer.data(mi + 1361);
    const auto *mi_1362 = buffer.data(mi + 1362);
    const auto *mi_1363 = buffer.data(mi + 1363);
    const auto *mi_1364 = buffer.data(mi + 1364);
    const auto *mi_1365 = buffer.data(mi + 1365);
    const auto *mi_1366 = buffer.data(mi + 1366);
    const auto *mi_1367 = buffer.data(mi + 1367);
    const auto *mi_1368 = buffer.data(mi + 1368);
    const auto *mi_1369 = buffer.data(mi + 1369);
    const auto *mi_1370 = buffer.data(mi + 1370);
    const auto *mi_1371 = buffer.data(mi + 1371);
    const auto *mi_1372 = buffer.data(mi + 1372);
    const auto *mi_1373 = buffer.data(mi + 1373);
    const auto *mi_1374 = buffer.data(mi + 1374);
    const auto *mi_1375 = buffer.data(mi + 1375);
    const auto *mi_1376 = buffer.data(mi + 1376);
    const auto *mi_1377 = buffer.data(mi + 1377);
    const auto *mi_1378 = buffer.data(mi + 1378);
    const auto *mi_1379 = buffer.data(mi + 1379);
    const auto *mi_1380 = buffer.data(mi + 1380);
    const auto *mi_1381 = buffer.data(mi + 1381);
    const auto *mi_1382 = buffer.data(mi + 1382);
    const auto *mi_1383 = buffer.data(mi + 1383);
    const auto *mi_1384 = buffer.data(mi + 1384);
    const auto *mi_1385 = buffer.data(mi + 1385);
    const auto *mi_1386 = buffer.data(mi + 1386);
    const auto *mi_1387 = buffer.data(mi + 1387);
    const auto *mi_1388 = buffer.data(mi + 1388);
    const auto *mi_1389 = buffer.data(mi + 1389);
    const auto *mi_1390 = buffer.data(mi + 1390);
    const auto *mi_1391 = buffer.data(mi + 1391);
    const auto *mi_1392 = buffer.data(mi + 1392);
    const auto *mi_1393 = buffer.data(mi + 1393);
    const auto *mi_1394 = buffer.data(mi + 1394);
    const auto *mi_1395 = buffer.data(mi + 1395);
    const auto *mi_1396 = buffer.data(mi + 1396);
    const auto *mi_1397 = buffer.data(mi + 1397);
    const auto *mi_1398 = buffer.data(mi + 1398);
    const auto *mi_1399 = buffer.data(mi + 1399);
    const auto *mi_1400 = buffer.data(mi + 1400);
    const auto *mi_1401 = buffer.data(mi + 1401);
    const auto *mi_1402 = buffer.data(mi + 1402);
    const auto *mi_1403 = buffer.data(mi + 1403);
    const auto *mi_1404 = buffer.data(mi + 1404);
    const auto *mi_1405 = buffer.data(mi + 1405);
    const auto *mi_1406 = buffer.data(mi + 1406);
    const auto *mi_1407 = buffer.data(mi + 1407);
    const auto *mi_1408 = buffer.data(mi + 1408);
    const auto *mi_1409 = buffer.data(mi + 1409);
    const auto *mi_1410 = buffer.data(mi + 1410);
    const auto *mi_1411 = buffer.data(mi + 1411);
    const auto *mi_1412 = buffer.data(mi + 1412);
    const auto *mi_1413 = buffer.data(mi + 1413);
    const auto *mi_1414 = buffer.data(mi + 1414);
    const auto *mi_1415 = buffer.data(mi + 1415);
    const auto *mi_1416 = buffer.data(mi + 1416);

#pragma omp simd aligned(t_977, t_978, t_979, t_980, t_981, ki_753, ki_754, ki_755, ki_756, \
                         ki_757, mi_1229, mi_1230, mi_1231, mi_1232, \
                         mi_1233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = -6.0 * ki_753[k]
                   + f_0 * mi_1229[k];

        t_978[k] = -6.0 * ki_754[k]
                   + f_0 * mi_1230[k];

        t_979[k] = -6.0 * ki_755[k]
                   + f_0 * mi_1231[k];

        t_980[k] = -7.0 * ki_756[k]
                   + f_0 * mi_1232[k];

        t_981[k] = -7.0 * ki_757[k]
                   + f_0 * mi_1233[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, t_985, t_986, ki_758, ki_759, ki_760, ki_761, \
                         ki_762, mi_1234, mi_1235, mi_1236, mi_1237, \
                         mi_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = -7.0 * ki_758[k]
                   + f_0 * mi_1234[k];

        t_983[k] = -7.0 * ki_759[k]
                   + f_0 * mi_1235[k];

        t_984[k] = -7.0 * ki_760[k]
                   + f_0 * mi_1236[k];

        t_985[k] = -7.0 * ki_761[k]
                   + f_0 * mi_1237[k];

        t_986[k] = -7.0 * ki_762[k]
                   + f_0 * mi_1238[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, t_990, t_991, ki_763, ki_764, ki_765, ki_766, \
                         ki_767, mi_1239, mi_1240, mi_1241, mi_1242, \
                         mi_1243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = -7.0 * ki_763[k]
                   + f_0 * mi_1239[k];

        t_988[k] = -7.0 * ki_764[k]
                   + f_0 * mi_1240[k];

        t_989[k] = -7.0 * ki_765[k]
                   + f_0 * mi_1241[k];

        t_990[k] = -7.0 * ki_766[k]
                   + f_0 * mi_1242[k];

        t_991[k] = -7.0 * ki_767[k]
                   + f_0 * mi_1243[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, t_996, ki_768, ki_769, ki_770, ki_771, \
                         ki_772, mi_1244, mi_1245, mi_1246, mi_1247, \
                         mi_1248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = -7.0 * ki_768[k]
                   + f_0 * mi_1244[k];

        t_993[k] = -7.0 * ki_769[k]
                   + f_0 * mi_1245[k];

        t_994[k] = -7.0 * ki_770[k]
                   + f_0 * mi_1246[k];

        t_995[k] = -7.0 * ki_771[k]
                   + f_0 * mi_1247[k];

        t_996[k] = -7.0 * ki_772[k]
                   + f_0 * mi_1248[k];
    }

#pragma omp simd aligned(t_997, t_998, t_999, t_1000, t_1001, ki_773, ki_774, ki_775, ki_776, \
                         ki_777, mi_1249, mi_1250, mi_1251, mi_1252, \
                         mi_1253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_997[k] = -7.0 * ki_773[k]
                   + f_0 * mi_1249[k];

        t_998[k] = -7.0 * ki_774[k]
                   + f_0 * mi_1250[k];

        t_999[k] = -7.0 * ki_775[k]
                   + f_0 * mi_1251[k];

        t_1000[k] = -7.0 * ki_776[k]
                    + f_0 * mi_1252[k];

        t_1001[k] = -7.0 * ki_777[k]
                    + f_0 * mi_1253[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, t_1005, t_1006, ki_778, ki_779, ki_780, \
                         ki_781, ki_782, mi_1254, mi_1255, mi_1256, mi_1257, \
                         mi_1258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = -7.0 * ki_778[k]
                    + f_0 * mi_1254[k];

        t_1003[k] = -7.0 * ki_779[k]
                    + f_0 * mi_1255[k];

        t_1004[k] = -7.0 * ki_780[k]
                    + f_0 * mi_1256[k];

        t_1005[k] = -7.0 * ki_781[k]
                    + f_0 * mi_1257[k];

        t_1006[k] = -7.0 * ki_782[k]
                    + f_0 * mi_1258[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, t_1010, t_1011, t_1012, t_1013, ki_783, \
                         mi_1259, mi_1288, mi_1289, mi_1290, mi_1291, mi_1292, \
                         mi_1293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = -7.0 * ki_783[k]
                    + f_0 * mi_1259[k];

        t_1008[k] = f_0 * mi_1288[k];

        t_1009[k] = f_0 * mi_1289[k];

        t_1010[k] = f_0 * mi_1290[k];

        t_1011[k] = f_0 * mi_1291[k];

        t_1012[k] = f_0 * mi_1292[k];

        t_1013[k] = f_0 * mi_1293[k];
    }

#pragma omp simd aligned(t_1014, t_1015, t_1016, t_1017, t_1018, t_1019, t_1020, t_1021, \
                         mi_1294, mi_1295, mi_1296, mi_1297, mi_1298, mi_1299, mi_1300, \
                         mi_1301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1014[k] = f_0 * mi_1294[k];

        t_1015[k] = f_0 * mi_1295[k];

        t_1016[k] = f_0 * mi_1296[k];

        t_1017[k] = f_0 * mi_1297[k];

        t_1018[k] = f_0 * mi_1298[k];

        t_1019[k] = f_0 * mi_1299[k];

        t_1020[k] = f_0 * mi_1300[k];

        t_1021[k] = f_0 * mi_1301[k];
    }

#pragma omp simd aligned(t_1022, t_1023, t_1024, t_1025, t_1026, t_1027, t_1028, t_1029, \
                         mi_1302, mi_1303, mi_1304, mi_1305, mi_1306, mi_1307, mi_1308, \
                         mi_1309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1022[k] = f_0 * mi_1302[k];

        t_1023[k] = f_0 * mi_1303[k];

        t_1024[k] = f_0 * mi_1304[k];

        t_1025[k] = f_0 * mi_1305[k];

        t_1026[k] = f_0 * mi_1306[k];

        t_1027[k] = f_0 * mi_1307[k];

        t_1028[k] = f_0 * mi_1308[k];

        t_1029[k] = f_0 * mi_1309[k];
    }

#pragma omp simd aligned(t_1030, t_1031, t_1032, t_1033, t_1034, t_1035, t_1036, ki_784, \
                         mi_1310, mi_1311, mi_1312, mi_1313, mi_1314, mi_1315, \
                         mi_1316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1030[k] = f_0 * mi_1310[k];

        t_1031[k] = f_0 * mi_1311[k];

        t_1032[k] = f_0 * mi_1312[k];

        t_1033[k] = f_0 * mi_1313[k];

        t_1034[k] = f_0 * mi_1314[k];

        t_1035[k] = f_0 * mi_1315[k];

        t_1036[k] = -ki_784[k]
                    + f_0 * mi_1316[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, t_1040, t_1041, ki_785, ki_786, ki_787, \
                         ki_788, ki_789, mi_1317, mi_1318, mi_1319, mi_1320, \
                         mi_1321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = -ki_785[k]
                    + f_0 * mi_1317[k];

        t_1038[k] = -ki_786[k]
                    + f_0 * mi_1318[k];

        t_1039[k] = -ki_787[k]
                    + f_0 * mi_1319[k];

        t_1040[k] = -ki_788[k]
                    + f_0 * mi_1320[k];

        t_1041[k] = -ki_789[k]
                    + f_0 * mi_1321[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, t_1046, ki_790, ki_791, ki_792, \
                         ki_793, ki_794, mi_1322, mi_1323, mi_1324, mi_1325, \
                         mi_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = -ki_790[k]
                    + f_0 * mi_1322[k];

        t_1043[k] = -ki_791[k]
                    + f_0 * mi_1323[k];

        t_1044[k] = -ki_792[k]
                    + f_0 * mi_1324[k];

        t_1045[k] = -ki_793[k]
                    + f_0 * mi_1325[k];

        t_1046[k] = -ki_794[k]
                    + f_0 * mi_1326[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, t_1050, t_1051, ki_795, ki_796, ki_797, \
                         ki_798, ki_799, mi_1327, mi_1328, mi_1329, mi_1330, \
                         mi_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = -ki_795[k]
                    + f_0 * mi_1327[k];

        t_1048[k] = -ki_796[k]
                    + f_0 * mi_1328[k];

        t_1049[k] = -ki_797[k]
                    + f_0 * mi_1329[k];

        t_1050[k] = -ki_798[k]
                    + f_0 * mi_1330[k];

        t_1051[k] = -ki_799[k]
                    + f_0 * mi_1331[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, t_1055, t_1056, ki_800, ki_801, ki_802, \
                         ki_803, ki_804, mi_1332, mi_1333, mi_1334, mi_1335, \
                         mi_1336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = -ki_800[k]
                    + f_0 * mi_1332[k];

        t_1053[k] = -ki_801[k]
                    + f_0 * mi_1333[k];

        t_1054[k] = -ki_802[k]
                    + f_0 * mi_1334[k];

        t_1055[k] = -ki_803[k]
                    + f_0 * mi_1335[k];

        t_1056[k] = -ki_804[k]
                    + f_0 * mi_1336[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, t_1060, t_1061, ki_805, ki_806, ki_807, \
                         ki_808, ki_809, mi_1337, mi_1338, mi_1339, mi_1340, \
                         mi_1341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = -ki_805[k]
                    + f_0 * mi_1337[k];

        t_1058[k] = -ki_806[k]
                    + f_0 * mi_1338[k];

        t_1059[k] = -ki_807[k]
                    + f_0 * mi_1339[k];

        t_1060[k] = -ki_808[k]
                    + f_0 * mi_1340[k];

        t_1061[k] = -ki_809[k]
                    + f_0 * mi_1341[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, t_1065, t_1066, ki_810, ki_811, ki_812, \
                         ki_813, ki_814, mi_1342, mi_1343, mi_1344, mi_1345, \
                         mi_1346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = -ki_810[k]
                    + f_0 * mi_1342[k];

        t_1063[k] = -ki_811[k]
                    + f_0 * mi_1343[k];

        t_1064[k] = -2.0 * ki_812[k]
                    + f_0 * mi_1344[k];

        t_1065[k] = -2.0 * ki_813[k]
                    + f_0 * mi_1345[k];

        t_1066[k] = -2.0 * ki_814[k]
                    + f_0 * mi_1346[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, t_1070, t_1071, ki_815, ki_816, ki_817, \
                         ki_818, ki_819, mi_1347, mi_1348, mi_1349, mi_1350, \
                         mi_1351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = -2.0 * ki_815[k]
                    + f_0 * mi_1347[k];

        t_1068[k] = -2.0 * ki_816[k]
                    + f_0 * mi_1348[k];

        t_1069[k] = -2.0 * ki_817[k]
                    + f_0 * mi_1349[k];

        t_1070[k] = -2.0 * ki_818[k]
                    + f_0 * mi_1350[k];

        t_1071[k] = -2.0 * ki_819[k]
                    + f_0 * mi_1351[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, t_1075, t_1076, ki_820, ki_821, ki_822, \
                         ki_823, ki_824, mi_1352, mi_1353, mi_1354, mi_1355, \
                         mi_1356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = -2.0 * ki_820[k]
                    + f_0 * mi_1352[k];

        t_1073[k] = -2.0 * ki_821[k]
                    + f_0 * mi_1353[k];

        t_1074[k] = -2.0 * ki_822[k]
                    + f_0 * mi_1354[k];

        t_1075[k] = -2.0 * ki_823[k]
                    + f_0 * mi_1355[k];

        t_1076[k] = -2.0 * ki_824[k]
                    + f_0 * mi_1356[k];
    }

#pragma omp simd aligned(t_1077, t_1078, t_1079, t_1080, t_1081, ki_825, ki_826, ki_827, \
                         ki_828, ki_829, mi_1357, mi_1358, mi_1359, mi_1360, \
                         mi_1361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = -2.0 * ki_825[k]
                    + f_0 * mi_1357[k];

        t_1078[k] = -2.0 * ki_826[k]
                    + f_0 * mi_1358[k];

        t_1079[k] = -2.0 * ki_827[k]
                    + f_0 * mi_1359[k];

        t_1080[k] = -2.0 * ki_828[k]
                    + f_0 * mi_1360[k];

        t_1081[k] = -2.0 * ki_829[k]
                    + f_0 * mi_1361[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, t_1085, t_1086, ki_830, ki_831, ki_832, \
                         ki_833, ki_834, mi_1362, mi_1363, mi_1364, mi_1365, \
                         mi_1366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = -2.0 * ki_830[k]
                    + f_0 * mi_1362[k];

        t_1083[k] = -2.0 * ki_831[k]
                    + f_0 * mi_1363[k];

        t_1084[k] = -2.0 * ki_832[k]
                    + f_0 * mi_1364[k];

        t_1085[k] = -2.0 * ki_833[k]
                    + f_0 * mi_1365[k];

        t_1086[k] = -2.0 * ki_834[k]
                    + f_0 * mi_1366[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, t_1090, t_1091, ki_835, ki_836, ki_837, \
                         ki_838, ki_839, mi_1367, mi_1368, mi_1369, mi_1370, \
                         mi_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = -2.0 * ki_835[k]
                    + f_0 * mi_1367[k];

        t_1088[k] = -2.0 * ki_836[k]
                    + f_0 * mi_1368[k];

        t_1089[k] = -2.0 * ki_837[k]
                    + f_0 * mi_1369[k];

        t_1090[k] = -2.0 * ki_838[k]
                    + f_0 * mi_1370[k];

        t_1091[k] = -2.0 * ki_839[k]
                    + f_0 * mi_1371[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, t_1095, t_1096, ki_840, ki_841, ki_842, \
                         ki_843, ki_844, mi_1372, mi_1373, mi_1374, mi_1375, \
                         mi_1376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = -3.0 * ki_840[k]
                    + f_0 * mi_1372[k];

        t_1093[k] = -3.0 * ki_841[k]
                    + f_0 * mi_1373[k];

        t_1094[k] = -3.0 * ki_842[k]
                    + f_0 * mi_1374[k];

        t_1095[k] = -3.0 * ki_843[k]
                    + f_0 * mi_1375[k];

        t_1096[k] = -3.0 * ki_844[k]
                    + f_0 * mi_1376[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, t_1100, t_1101, ki_845, ki_846, ki_847, \
                         ki_848, ki_849, mi_1377, mi_1378, mi_1379, mi_1380, \
                         mi_1381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = -3.0 * ki_845[k]
                    + f_0 * mi_1377[k];

        t_1098[k] = -3.0 * ki_846[k]
                    + f_0 * mi_1378[k];

        t_1099[k] = -3.0 * ki_847[k]
                    + f_0 * mi_1379[k];

        t_1100[k] = -3.0 * ki_848[k]
                    + f_0 * mi_1380[k];

        t_1101[k] = -3.0 * ki_849[k]
                    + f_0 * mi_1381[k];
    }

#pragma omp simd aligned(t_1102, t_1103, t_1104, t_1105, t_1106, ki_850, ki_851, ki_852, \
                         ki_853, ki_854, mi_1382, mi_1383, mi_1384, mi_1385, \
                         mi_1386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1102[k] = -3.0 * ki_850[k]
                    + f_0 * mi_1382[k];

        t_1103[k] = -3.0 * ki_851[k]
                    + f_0 * mi_1383[k];

        t_1104[k] = -3.0 * ki_852[k]
                    + f_0 * mi_1384[k];

        t_1105[k] = -3.0 * ki_853[k]
                    + f_0 * mi_1385[k];

        t_1106[k] = -3.0 * ki_854[k]
                    + f_0 * mi_1386[k];
    }

#pragma omp simd aligned(t_1107, t_1108, t_1109, t_1110, t_1111, ki_855, ki_856, ki_857, \
                         ki_858, ki_859, mi_1387, mi_1388, mi_1389, mi_1390, \
                         mi_1391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1107[k] = -3.0 * ki_855[k]
                    + f_0 * mi_1387[k];

        t_1108[k] = -3.0 * ki_856[k]
                    + f_0 * mi_1388[k];

        t_1109[k] = -3.0 * ki_857[k]
                    + f_0 * mi_1389[k];

        t_1110[k] = -3.0 * ki_858[k]
                    + f_0 * mi_1390[k];

        t_1111[k] = -3.0 * ki_859[k]
                    + f_0 * mi_1391[k];
    }

#pragma omp simd aligned(t_1112, t_1113, t_1114, t_1115, t_1116, ki_860, ki_861, ki_862, \
                         ki_863, ki_864, mi_1392, mi_1393, mi_1394, mi_1395, \
                         mi_1396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1112[k] = -3.0 * ki_860[k]
                    + f_0 * mi_1392[k];

        t_1113[k] = -3.0 * ki_861[k]
                    + f_0 * mi_1393[k];

        t_1114[k] = -3.0 * ki_862[k]
                    + f_0 * mi_1394[k];

        t_1115[k] = -3.0 * ki_863[k]
                    + f_0 * mi_1395[k];

        t_1116[k] = -3.0 * ki_864[k]
                    + f_0 * mi_1396[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, t_1120, t_1121, ki_865, ki_866, ki_867, \
                         ki_868, ki_869, mi_1397, mi_1398, mi_1399, mi_1400, \
                         mi_1401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = -3.0 * ki_865[k]
                    + f_0 * mi_1397[k];

        t_1118[k] = -3.0 * ki_866[k]
                    + f_0 * mi_1398[k];

        t_1119[k] = -3.0 * ki_867[k]
                    + f_0 * mi_1399[k];

        t_1120[k] = -4.0 * ki_868[k]
                    + f_0 * mi_1400[k];

        t_1121[k] = -4.0 * ki_869[k]
                    + f_0 * mi_1401[k];
    }

#pragma omp simd aligned(t_1122, t_1123, t_1124, t_1125, t_1126, ki_870, ki_871, ki_872, \
                         ki_873, ki_874, mi_1402, mi_1403, mi_1404, mi_1405, \
                         mi_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1122[k] = -4.0 * ki_870[k]
                    + f_0 * mi_1402[k];

        t_1123[k] = -4.0 * ki_871[k]
                    + f_0 * mi_1403[k];

        t_1124[k] = -4.0 * ki_872[k]
                    + f_0 * mi_1404[k];

        t_1125[k] = -4.0 * ki_873[k]
                    + f_0 * mi_1405[k];

        t_1126[k] = -4.0 * ki_874[k]
                    + f_0 * mi_1406[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, t_1130, t_1131, ki_875, ki_876, ki_877, \
                         ki_878, ki_879, mi_1407, mi_1408, mi_1409, mi_1410, \
                         mi_1411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = -4.0 * ki_875[k]
                    + f_0 * mi_1407[k];

        t_1128[k] = -4.0 * ki_876[k]
                    + f_0 * mi_1408[k];

        t_1129[k] = -4.0 * ki_877[k]
                    + f_0 * mi_1409[k];

        t_1130[k] = -4.0 * ki_878[k]
                    + f_0 * mi_1410[k];

        t_1131[k] = -4.0 * ki_879[k]
                    + f_0 * mi_1411[k];
    }

#pragma omp simd aligned(t_1132, t_1133, t_1134, t_1135, t_1136, ki_880, ki_881, ki_882, \
                         ki_883, ki_884, mi_1412, mi_1413, mi_1414, mi_1415, \
                         mi_1416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = -4.0 * ki_880[k]
                    + f_0 * mi_1412[k];

        t_1133[k] = -4.0 * ki_881[k]
                    + f_0 * mi_1413[k];

        t_1134[k] = -4.0 * ki_882[k]
                    + f_0 * mi_1414[k];

        t_1135[k] = -4.0 * ki_883[k]
                    + f_0 * mi_1415[k];

        t_1136[k] = -4.0 * ki_884[k]
                    + f_0 * mi_1416[k];
    }
}

static auto
compute_prim_geom_10_li_electron_repulsion_2_piece7(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t mi,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);
    auto *t_1169 = buffer.data(target + 1169);
    auto *t_1170 = buffer.data(target + 1170);
    auto *t_1171 = buffer.data(target + 1171);
    auto *t_1172 = buffer.data(target + 1172);
    auto *t_1173 = buffer.data(target + 1173);
    auto *t_1174 = buffer.data(target + 1174);
    auto *t_1175 = buffer.data(target + 1175);
    auto *t_1176 = buffer.data(target + 1176);
    auto *t_1177 = buffer.data(target + 1177);
    auto *t_1178 = buffer.data(target + 1178);
    auto *t_1179 = buffer.data(target + 1179);
    auto *t_1180 = buffer.data(target + 1180);
    auto *t_1181 = buffer.data(target + 1181);
    auto *t_1182 = buffer.data(target + 1182);
    auto *t_1183 = buffer.data(target + 1183);
    auto *t_1184 = buffer.data(target + 1184);
    auto *t_1185 = buffer.data(target + 1185);
    auto *t_1186 = buffer.data(target + 1186);
    auto *t_1187 = buffer.data(target + 1187);
    auto *t_1188 = buffer.data(target + 1188);
    auto *t_1189 = buffer.data(target + 1189);
    auto *t_1190 = buffer.data(target + 1190);
    auto *t_1191 = buffer.data(target + 1191);
    auto *t_1192 = buffer.data(target + 1192);
    auto *t_1193 = buffer.data(target + 1193);
    auto *t_1194 = buffer.data(target + 1194);
    auto *t_1195 = buffer.data(target + 1195);
    auto *t_1196 = buffer.data(target + 1196);
    auto *t_1197 = buffer.data(target + 1197);
    auto *t_1198 = buffer.data(target + 1198);
    auto *t_1199 = buffer.data(target + 1199);
    auto *t_1200 = buffer.data(target + 1200);
    auto *t_1201 = buffer.data(target + 1201);
    auto *t_1202 = buffer.data(target + 1202);
    auto *t_1203 = buffer.data(target + 1203);
    auto *t_1204 = buffer.data(target + 1204);
    auto *t_1205 = buffer.data(target + 1205);
    auto *t_1206 = buffer.data(target + 1206);
    auto *t_1207 = buffer.data(target + 1207);
    auto *t_1208 = buffer.data(target + 1208);
    auto *t_1209 = buffer.data(target + 1209);
    auto *t_1210 = buffer.data(target + 1210);
    auto *t_1211 = buffer.data(target + 1211);
    auto *t_1212 = buffer.data(target + 1212);
    auto *t_1213 = buffer.data(target + 1213);
    auto *t_1214 = buffer.data(target + 1214);
    auto *t_1215 = buffer.data(target + 1215);
    auto *t_1216 = buffer.data(target + 1216);
    auto *t_1217 = buffer.data(target + 1217);
    auto *t_1218 = buffer.data(target + 1218);
    auto *t_1219 = buffer.data(target + 1219);
    auto *t_1220 = buffer.data(target + 1220);
    auto *t_1221 = buffer.data(target + 1221);
    auto *t_1222 = buffer.data(target + 1222);
    auto *t_1223 = buffer.data(target + 1223);
    auto *t_1224 = buffer.data(target + 1224);
    auto *t_1225 = buffer.data(target + 1225);
    auto *t_1226 = buffer.data(target + 1226);
    auto *t_1227 = buffer.data(target + 1227);
    auto *t_1228 = buffer.data(target + 1228);
    auto *t_1229 = buffer.data(target + 1229);
    auto *t_1230 = buffer.data(target + 1230);
    auto *t_1231 = buffer.data(target + 1231);
    auto *t_1232 = buffer.data(target + 1232);
    auto *t_1233 = buffer.data(target + 1233);
    auto *t_1234 = buffer.data(target + 1234);
    auto *t_1235 = buffer.data(target + 1235);
    auto *t_1236 = buffer.data(target + 1236);
    auto *t_1237 = buffer.data(target + 1237);
    auto *t_1238 = buffer.data(target + 1238);
    auto *t_1239 = buffer.data(target + 1239);
    auto *t_1240 = buffer.data(target + 1240);
    auto *t_1241 = buffer.data(target + 1241);
    auto *t_1242 = buffer.data(target + 1242);
    auto *t_1243 = buffer.data(target + 1243);
    auto *t_1244 = buffer.data(target + 1244);
    auto *t_1245 = buffer.data(target + 1245);
    auto *t_1246 = buffer.data(target + 1246);
    auto *t_1247 = buffer.data(target + 1247);
    auto *t_1248 = buffer.data(target + 1248);
    auto *t_1249 = buffer.data(target + 1249);
    auto *t_1250 = buffer.data(target + 1250);
    auto *t_1251 = buffer.data(target + 1251);
    auto *t_1252 = buffer.data(target + 1252);
    auto *t_1253 = buffer.data(target + 1253);
    auto *t_1254 = buffer.data(target + 1254);
    auto *t_1255 = buffer.data(target + 1255);
    auto *t_1256 = buffer.data(target + 1256);
    auto *t_1257 = buffer.data(target + 1257);
    auto *t_1258 = buffer.data(target + 1258);
    auto *t_1259 = buffer.data(target + 1259);

    const auto *ki_885 = buffer.data(ki + 885);
    const auto *ki_886 = buffer.data(ki + 886);
    const auto *ki_887 = buffer.data(ki + 887);
    const auto *ki_888 = buffer.data(ki + 888);
    const auto *ki_889 = buffer.data(ki + 889);
    const auto *ki_890 = buffer.data(ki + 890);
    const auto *ki_891 = buffer.data(ki + 891);
    const auto *ki_892 = buffer.data(ki + 892);
    const auto *ki_893 = buffer.data(ki + 893);
    const auto *ki_894 = buffer.data(ki + 894);
    const auto *ki_895 = buffer.data(ki + 895);
    const auto *ki_896 = buffer.data(ki + 896);
    const auto *ki_897 = buffer.data(ki + 897);
    const auto *ki_898 = buffer.data(ki + 898);
    const auto *ki_899 = buffer.data(ki + 899);
    const auto *ki_900 = buffer.data(ki + 900);
    const auto *ki_901 = buffer.data(ki + 901);
    const auto *ki_902 = buffer.data(ki + 902);
    const auto *ki_903 = buffer.data(ki + 903);
    const auto *ki_904 = buffer.data(ki + 904);
    const auto *ki_905 = buffer.data(ki + 905);
    const auto *ki_906 = buffer.data(ki + 906);
    const auto *ki_907 = buffer.data(ki + 907);
    const auto *ki_908 = buffer.data(ki + 908);
    const auto *ki_909 = buffer.data(ki + 909);
    const auto *ki_910 = buffer.data(ki + 910);
    const auto *ki_911 = buffer.data(ki + 911);
    const auto *ki_912 = buffer.data(ki + 912);
    const auto *ki_913 = buffer.data(ki + 913);
    const auto *ki_914 = buffer.data(ki + 914);
    const auto *ki_915 = buffer.data(ki + 915);
    const auto *ki_916 = buffer.data(ki + 916);
    const auto *ki_917 = buffer.data(ki + 917);
    const auto *ki_918 = buffer.data(ki + 918);
    const auto *ki_919 = buffer.data(ki + 919);
    const auto *ki_920 = buffer.data(ki + 920);
    const auto *ki_921 = buffer.data(ki + 921);
    const auto *ki_922 = buffer.data(ki + 922);
    const auto *ki_923 = buffer.data(ki + 923);
    const auto *ki_924 = buffer.data(ki + 924);
    const auto *ki_925 = buffer.data(ki + 925);
    const auto *ki_926 = buffer.data(ki + 926);
    const auto *ki_927 = buffer.data(ki + 927);
    const auto *ki_928 = buffer.data(ki + 928);
    const auto *ki_929 = buffer.data(ki + 929);
    const auto *ki_930 = buffer.data(ki + 930);
    const auto *ki_931 = buffer.data(ki + 931);
    const auto *ki_932 = buffer.data(ki + 932);
    const auto *ki_933 = buffer.data(ki + 933);
    const auto *ki_934 = buffer.data(ki + 934);
    const auto *ki_935 = buffer.data(ki + 935);
    const auto *ki_936 = buffer.data(ki + 936);
    const auto *ki_937 = buffer.data(ki + 937);
    const auto *ki_938 = buffer.data(ki + 938);
    const auto *ki_939 = buffer.data(ki + 939);
    const auto *ki_940 = buffer.data(ki + 940);
    const auto *ki_941 = buffer.data(ki + 941);
    const auto *ki_942 = buffer.data(ki + 942);
    const auto *ki_943 = buffer.data(ki + 943);
    const auto *ki_944 = buffer.data(ki + 944);
    const auto *ki_945 = buffer.data(ki + 945);
    const auto *ki_946 = buffer.data(ki + 946);
    const auto *ki_947 = buffer.data(ki + 947);
    const auto *ki_948 = buffer.data(ki + 948);
    const auto *ki_949 = buffer.data(ki + 949);
    const auto *ki_950 = buffer.data(ki + 950);
    const auto *ki_951 = buffer.data(ki + 951);
    const auto *ki_952 = buffer.data(ki + 952);
    const auto *ki_953 = buffer.data(ki + 953);
    const auto *ki_954 = buffer.data(ki + 954);
    const auto *ki_955 = buffer.data(ki + 955);
    const auto *ki_956 = buffer.data(ki + 956);
    const auto *ki_957 = buffer.data(ki + 957);
    const auto *ki_958 = buffer.data(ki + 958);
    const auto *ki_959 = buffer.data(ki + 959);
    const auto *ki_960 = buffer.data(ki + 960);
    const auto *ki_961 = buffer.data(ki + 961);
    const auto *ki_962 = buffer.data(ki + 962);
    const auto *ki_963 = buffer.data(ki + 963);
    const auto *ki_964 = buffer.data(ki + 964);
    const auto *ki_965 = buffer.data(ki + 965);
    const auto *ki_966 = buffer.data(ki + 966);
    const auto *ki_967 = buffer.data(ki + 967);
    const auto *ki_968 = buffer.data(ki + 968);
    const auto *ki_969 = buffer.data(ki + 969);
    const auto *ki_970 = buffer.data(ki + 970);
    const auto *ki_971 = buffer.data(ki + 971);
    const auto *ki_972 = buffer.data(ki + 972);
    const auto *ki_973 = buffer.data(ki + 973);
    const auto *ki_974 = buffer.data(ki + 974);
    const auto *ki_975 = buffer.data(ki + 975);
    const auto *ki_976 = buffer.data(ki + 976);
    const auto *ki_977 = buffer.data(ki + 977);
    const auto *ki_978 = buffer.data(ki + 978);
    const auto *ki_979 = buffer.data(ki + 979);
    const auto *ki_980 = buffer.data(ki + 980);
    const auto *ki_981 = buffer.data(ki + 981);
    const auto *ki_982 = buffer.data(ki + 982);
    const auto *ki_983 = buffer.data(ki + 983);
    const auto *ki_984 = buffer.data(ki + 984);
    const auto *ki_985 = buffer.data(ki + 985);
    const auto *ki_986 = buffer.data(ki + 986);
    const auto *ki_987 = buffer.data(ki + 987);
    const auto *ki_988 = buffer.data(ki + 988);
    const auto *ki_989 = buffer.data(ki + 989);
    const auto *ki_990 = buffer.data(ki + 990);
    const auto *ki_991 = buffer.data(ki + 991);
    const auto *ki_992 = buffer.data(ki + 992);
    const auto *ki_993 = buffer.data(ki + 993);
    const auto *ki_994 = buffer.data(ki + 994);
    const auto *ki_995 = buffer.data(ki + 995);
    const auto *ki_996 = buffer.data(ki + 996);
    const auto *ki_997 = buffer.data(ki + 997);
    const auto *ki_998 = buffer.data(ki + 998);
    const auto *ki_999 = buffer.data(ki + 999);
    const auto *ki_1000 = buffer.data(ki + 1000);
    const auto *ki_1001 = buffer.data(ki + 1001);
    const auto *ki_1002 = buffer.data(ki + 1002);
    const auto *ki_1003 = buffer.data(ki + 1003);
    const auto *ki_1004 = buffer.data(ki + 1004);
    const auto *ki_1005 = buffer.data(ki + 1005);
    const auto *ki_1006 = buffer.data(ki + 1006);
    const auto *ki_1007 = buffer.data(ki + 1007);

    const auto *mi_1417 = buffer.data(mi + 1417);
    const auto *mi_1418 = buffer.data(mi + 1418);
    const auto *mi_1419 = buffer.data(mi + 1419);
    const auto *mi_1420 = buffer.data(mi + 1420);
    const auto *mi_1421 = buffer.data(mi + 1421);
    const auto *mi_1422 = buffer.data(mi + 1422);
    const auto *mi_1423 = buffer.data(mi + 1423);
    const auto *mi_1424 = buffer.data(mi + 1424);
    const auto *mi_1425 = buffer.data(mi + 1425);
    const auto *mi_1426 = buffer.data(mi + 1426);
    const auto *mi_1427 = buffer.data(mi + 1427);
    const auto *mi_1428 = buffer.data(mi + 1428);
    const auto *mi_1429 = buffer.data(mi + 1429);
    const auto *mi_1430 = buffer.data(mi + 1430);
    const auto *mi_1431 = buffer.data(mi + 1431);
    const auto *mi_1432 = buffer.data(mi + 1432);
    const auto *mi_1433 = buffer.data(mi + 1433);
    const auto *mi_1434 = buffer.data(mi + 1434);
    const auto *mi_1435 = buffer.data(mi + 1435);
    const auto *mi_1436 = buffer.data(mi + 1436);
    const auto *mi_1437 = buffer.data(mi + 1437);
    const auto *mi_1438 = buffer.data(mi + 1438);
    const auto *mi_1439 = buffer.data(mi + 1439);
    const auto *mi_1440 = buffer.data(mi + 1440);
    const auto *mi_1441 = buffer.data(mi + 1441);
    const auto *mi_1442 = buffer.data(mi + 1442);
    const auto *mi_1443 = buffer.data(mi + 1443);
    const auto *mi_1444 = buffer.data(mi + 1444);
    const auto *mi_1445 = buffer.data(mi + 1445);
    const auto *mi_1446 = buffer.data(mi + 1446);
    const auto *mi_1447 = buffer.data(mi + 1447);
    const auto *mi_1448 = buffer.data(mi + 1448);
    const auto *mi_1449 = buffer.data(mi + 1449);
    const auto *mi_1450 = buffer.data(mi + 1450);
    const auto *mi_1451 = buffer.data(mi + 1451);
    const auto *mi_1452 = buffer.data(mi + 1452);
    const auto *mi_1453 = buffer.data(mi + 1453);
    const auto *mi_1454 = buffer.data(mi + 1454);
    const auto *mi_1455 = buffer.data(mi + 1455);
    const auto *mi_1456 = buffer.data(mi + 1456);
    const auto *mi_1457 = buffer.data(mi + 1457);
    const auto *mi_1458 = buffer.data(mi + 1458);
    const auto *mi_1459 = buffer.data(mi + 1459);
    const auto *mi_1460 = buffer.data(mi + 1460);
    const auto *mi_1461 = buffer.data(mi + 1461);
    const auto *mi_1462 = buffer.data(mi + 1462);
    const auto *mi_1463 = buffer.data(mi + 1463);
    const auto *mi_1464 = buffer.data(mi + 1464);
    const auto *mi_1465 = buffer.data(mi + 1465);
    const auto *mi_1466 = buffer.data(mi + 1466);
    const auto *mi_1467 = buffer.data(mi + 1467);
    const auto *mi_1468 = buffer.data(mi + 1468);
    const auto *mi_1469 = buffer.data(mi + 1469);
    const auto *mi_1470 = buffer.data(mi + 1470);
    const auto *mi_1471 = buffer.data(mi + 1471);
    const auto *mi_1472 = buffer.data(mi + 1472);
    const auto *mi_1473 = buffer.data(mi + 1473);
    const auto *mi_1474 = buffer.data(mi + 1474);
    const auto *mi_1475 = buffer.data(mi + 1475);
    const auto *mi_1476 = buffer.data(mi + 1476);
    const auto *mi_1477 = buffer.data(mi + 1477);
    const auto *mi_1478 = buffer.data(mi + 1478);
    const auto *mi_1479 = buffer.data(mi + 1479);
    const auto *mi_1480 = buffer.data(mi + 1480);
    const auto *mi_1481 = buffer.data(mi + 1481);
    const auto *mi_1482 = buffer.data(mi + 1482);
    const auto *mi_1483 = buffer.data(mi + 1483);
    const auto *mi_1484 = buffer.data(mi + 1484);
    const auto *mi_1485 = buffer.data(mi + 1485);
    const auto *mi_1486 = buffer.data(mi + 1486);
    const auto *mi_1487 = buffer.data(mi + 1487);
    const auto *mi_1488 = buffer.data(mi + 1488);
    const auto *mi_1489 = buffer.data(mi + 1489);
    const auto *mi_1490 = buffer.data(mi + 1490);
    const auto *mi_1491 = buffer.data(mi + 1491);
    const auto *mi_1492 = buffer.data(mi + 1492);
    const auto *mi_1493 = buffer.data(mi + 1493);
    const auto *mi_1494 = buffer.data(mi + 1494);
    const auto *mi_1495 = buffer.data(mi + 1495);
    const auto *mi_1496 = buffer.data(mi + 1496);
    const auto *mi_1497 = buffer.data(mi + 1497);
    const auto *mi_1498 = buffer.data(mi + 1498);
    const auto *mi_1499 = buffer.data(mi + 1499);
    const auto *mi_1500 = buffer.data(mi + 1500);
    const auto *mi_1501 = buffer.data(mi + 1501);
    const auto *mi_1502 = buffer.data(mi + 1502);
    const auto *mi_1503 = buffer.data(mi + 1503);
    const auto *mi_1504 = buffer.data(mi + 1504);
    const auto *mi_1505 = buffer.data(mi + 1505);
    const auto *mi_1506 = buffer.data(mi + 1506);
    const auto *mi_1507 = buffer.data(mi + 1507);
    const auto *mi_1508 = buffer.data(mi + 1508);
    const auto *mi_1509 = buffer.data(mi + 1509);
    const auto *mi_1510 = buffer.data(mi + 1510);
    const auto *mi_1511 = buffer.data(mi + 1511);
    const auto *mi_1512 = buffer.data(mi + 1512);
    const auto *mi_1513 = buffer.data(mi + 1513);
    const auto *mi_1514 = buffer.data(mi + 1514);
    const auto *mi_1515 = buffer.data(mi + 1515);
    const auto *mi_1516 = buffer.data(mi + 1516);
    const auto *mi_1517 = buffer.data(mi + 1517);
    const auto *mi_1518 = buffer.data(mi + 1518);
    const auto *mi_1519 = buffer.data(mi + 1519);
    const auto *mi_1520 = buffer.data(mi + 1520);
    const auto *mi_1521 = buffer.data(mi + 1521);
    const auto *mi_1522 = buffer.data(mi + 1522);
    const auto *mi_1523 = buffer.data(mi + 1523);
    const auto *mi_1524 = buffer.data(mi + 1524);
    const auto *mi_1525 = buffer.data(mi + 1525);
    const auto *mi_1526 = buffer.data(mi + 1526);
    const auto *mi_1527 = buffer.data(mi + 1527);
    const auto *mi_1528 = buffer.data(mi + 1528);
    const auto *mi_1529 = buffer.data(mi + 1529);
    const auto *mi_1530 = buffer.data(mi + 1530);
    const auto *mi_1531 = buffer.data(mi + 1531);
    const auto *mi_1532 = buffer.data(mi + 1532);
    const auto *mi_1533 = buffer.data(mi + 1533);
    const auto *mi_1534 = buffer.data(mi + 1534);
    const auto *mi_1535 = buffer.data(mi + 1535);
    const auto *mi_1536 = buffer.data(mi + 1536);
    const auto *mi_1537 = buffer.data(mi + 1537);
    const auto *mi_1538 = buffer.data(mi + 1538);
    const auto *mi_1539 = buffer.data(mi + 1539);

#pragma omp simd aligned(t_1137, t_1138, t_1139, t_1140, t_1141, ki_885, ki_886, ki_887, \
                         ki_888, ki_889, mi_1417, mi_1418, mi_1419, mi_1420, \
                         mi_1421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1137[k] = -4.0 * ki_885[k]
                    + f_0 * mi_1417[k];

        t_1138[k] = -4.0 * ki_886[k]
                    + f_0 * mi_1418[k];

        t_1139[k] = -4.0 * ki_887[k]
                    + f_0 * mi_1419[k];

        t_1140[k] = -4.0 * ki_888[k]
                    + f_0 * mi_1420[k];

        t_1141[k] = -4.0 * ki_889[k]
                    + f_0 * mi_1421[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, t_1145, t_1146, ki_890, ki_891, ki_892, \
                         ki_893, ki_894, mi_1422, mi_1423, mi_1424, mi_1425, \
                         mi_1426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = -4.0 * ki_890[k]
                    + f_0 * mi_1422[k];

        t_1143[k] = -4.0 * ki_891[k]
                    + f_0 * mi_1423[k];

        t_1144[k] = -4.0 * ki_892[k]
                    + f_0 * mi_1424[k];

        t_1145[k] = -4.0 * ki_893[k]
                    + f_0 * mi_1425[k];

        t_1146[k] = -4.0 * ki_894[k]
                    + f_0 * mi_1426[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, t_1150, t_1151, ki_895, ki_896, ki_897, \
                         ki_898, ki_899, mi_1427, mi_1428, mi_1429, mi_1430, \
                         mi_1431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = -4.0 * ki_895[k]
                    + f_0 * mi_1427[k];

        t_1148[k] = -5.0 * ki_896[k]
                    + f_0 * mi_1428[k];

        t_1149[k] = -5.0 * ki_897[k]
                    + f_0 * mi_1429[k];

        t_1150[k] = -5.0 * ki_898[k]
                    + f_0 * mi_1430[k];

        t_1151[k] = -5.0 * ki_899[k]
                    + f_0 * mi_1431[k];
    }

#pragma omp simd aligned(t_1152, t_1153, t_1154, t_1155, t_1156, ki_900, ki_901, ki_902, \
                         ki_903, ki_904, mi_1432, mi_1433, mi_1434, mi_1435, \
                         mi_1436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1152[k] = -5.0 * ki_900[k]
                    + f_0 * mi_1432[k];

        t_1153[k] = -5.0 * ki_901[k]
                    + f_0 * mi_1433[k];

        t_1154[k] = -5.0 * ki_902[k]
                    + f_0 * mi_1434[k];

        t_1155[k] = -5.0 * ki_903[k]
                    + f_0 * mi_1435[k];

        t_1156[k] = -5.0 * ki_904[k]
                    + f_0 * mi_1436[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, t_1160, t_1161, ki_905, ki_906, ki_907, \
                         ki_908, ki_909, mi_1437, mi_1438, mi_1439, mi_1440, \
                         mi_1441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = -5.0 * ki_905[k]
                    + f_0 * mi_1437[k];

        t_1158[k] = -5.0 * ki_906[k]
                    + f_0 * mi_1438[k];

        t_1159[k] = -5.0 * ki_907[k]
                    + f_0 * mi_1439[k];

        t_1160[k] = -5.0 * ki_908[k]
                    + f_0 * mi_1440[k];

        t_1161[k] = -5.0 * ki_909[k]
                    + f_0 * mi_1441[k];
    }

#pragma omp simd aligned(t_1162, t_1163, t_1164, t_1165, t_1166, ki_910, ki_911, ki_912, \
                         ki_913, ki_914, mi_1442, mi_1443, mi_1444, mi_1445, \
                         mi_1446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1162[k] = -5.0 * ki_910[k]
                    + f_0 * mi_1442[k];

        t_1163[k] = -5.0 * ki_911[k]
                    + f_0 * mi_1443[k];

        t_1164[k] = -5.0 * ki_912[k]
                    + f_0 * mi_1444[k];

        t_1165[k] = -5.0 * ki_913[k]
                    + f_0 * mi_1445[k];

        t_1166[k] = -5.0 * ki_914[k]
                    + f_0 * mi_1446[k];
    }

#pragma omp simd aligned(t_1167, t_1168, t_1169, t_1170, t_1171, ki_915, ki_916, ki_917, \
                         ki_918, ki_919, mi_1447, mi_1448, mi_1449, mi_1450, \
                         mi_1451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1167[k] = -5.0 * ki_915[k]
                    + f_0 * mi_1447[k];

        t_1168[k] = -5.0 * ki_916[k]
                    + f_0 * mi_1448[k];

        t_1169[k] = -5.0 * ki_917[k]
                    + f_0 * mi_1449[k];

        t_1170[k] = -5.0 * ki_918[k]
                    + f_0 * mi_1450[k];

        t_1171[k] = -5.0 * ki_919[k]
                    + f_0 * mi_1451[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, t_1175, t_1176, ki_920, ki_921, ki_922, \
                         ki_923, ki_924, mi_1452, mi_1453, mi_1454, mi_1455, \
                         mi_1456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = -5.0 * ki_920[k]
                    + f_0 * mi_1452[k];

        t_1173[k] = -5.0 * ki_921[k]
                    + f_0 * mi_1453[k];

        t_1174[k] = -5.0 * ki_922[k]
                    + f_0 * mi_1454[k];

        t_1175[k] = -5.0 * ki_923[k]
                    + f_0 * mi_1455[k];

        t_1176[k] = -6.0 * ki_924[k]
                    + f_0 * mi_1456[k];
    }

#pragma omp simd aligned(t_1177, t_1178, t_1179, t_1180, t_1181, ki_925, ki_926, ki_927, \
                         ki_928, ki_929, mi_1457, mi_1458, mi_1459, mi_1460, \
                         mi_1461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1177[k] = -6.0 * ki_925[k]
                    + f_0 * mi_1457[k];

        t_1178[k] = -6.0 * ki_926[k]
                    + f_0 * mi_1458[k];

        t_1179[k] = -6.0 * ki_927[k]
                    + f_0 * mi_1459[k];

        t_1180[k] = -6.0 * ki_928[k]
                    + f_0 * mi_1460[k];

        t_1181[k] = -6.0 * ki_929[k]
                    + f_0 * mi_1461[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, t_1185, t_1186, ki_930, ki_931, ki_932, \
                         ki_933, ki_934, mi_1462, mi_1463, mi_1464, mi_1465, \
                         mi_1466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = -6.0 * ki_930[k]
                    + f_0 * mi_1462[k];

        t_1183[k] = -6.0 * ki_931[k]
                    + f_0 * mi_1463[k];

        t_1184[k] = -6.0 * ki_932[k]
                    + f_0 * mi_1464[k];

        t_1185[k] = -6.0 * ki_933[k]
                    + f_0 * mi_1465[k];

        t_1186[k] = -6.0 * ki_934[k]
                    + f_0 * mi_1466[k];
    }

#pragma omp simd aligned(t_1187, t_1188, t_1189, t_1190, t_1191, ki_935, ki_936, ki_937, \
                         ki_938, ki_939, mi_1467, mi_1468, mi_1469, mi_1470, \
                         mi_1471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1187[k] = -6.0 * ki_935[k]
                    + f_0 * mi_1467[k];

        t_1188[k] = -6.0 * ki_936[k]
                    + f_0 * mi_1468[k];

        t_1189[k] = -6.0 * ki_937[k]
                    + f_0 * mi_1469[k];

        t_1190[k] = -6.0 * ki_938[k]
                    + f_0 * mi_1470[k];

        t_1191[k] = -6.0 * ki_939[k]
                    + f_0 * mi_1471[k];
    }

#pragma omp simd aligned(t_1192, t_1193, t_1194, t_1195, t_1196, ki_940, ki_941, ki_942, \
                         ki_943, ki_944, mi_1472, mi_1473, mi_1474, mi_1475, \
                         mi_1476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1192[k] = -6.0 * ki_940[k]
                    + f_0 * mi_1472[k];

        t_1193[k] = -6.0 * ki_941[k]
                    + f_0 * mi_1473[k];

        t_1194[k] = -6.0 * ki_942[k]
                    + f_0 * mi_1474[k];

        t_1195[k] = -6.0 * ki_943[k]
                    + f_0 * mi_1475[k];

        t_1196[k] = -6.0 * ki_944[k]
                    + f_0 * mi_1476[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, t_1200, t_1201, ki_945, ki_946, ki_947, \
                         ki_948, ki_949, mi_1477, mi_1478, mi_1479, mi_1480, \
                         mi_1481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = -6.0 * ki_945[k]
                    + f_0 * mi_1477[k];

        t_1198[k] = -6.0 * ki_946[k]
                    + f_0 * mi_1478[k];

        t_1199[k] = -6.0 * ki_947[k]
                    + f_0 * mi_1479[k];

        t_1200[k] = -6.0 * ki_948[k]
                    + f_0 * mi_1480[k];

        t_1201[k] = -6.0 * ki_949[k]
                    + f_0 * mi_1481[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, t_1205, t_1206, ki_950, ki_951, ki_952, \
                         ki_953, ki_954, mi_1482, mi_1483, mi_1484, mi_1485, \
                         mi_1486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = -6.0 * ki_950[k]
                    + f_0 * mi_1482[k];

        t_1203[k] = -6.0 * ki_951[k]
                    + f_0 * mi_1483[k];

        t_1204[k] = -7.0 * ki_952[k]
                    + f_0 * mi_1484[k];

        t_1205[k] = -7.0 * ki_953[k]
                    + f_0 * mi_1485[k];

        t_1206[k] = -7.0 * ki_954[k]
                    + f_0 * mi_1486[k];
    }

#pragma omp simd aligned(t_1207, t_1208, t_1209, t_1210, t_1211, ki_955, ki_956, ki_957, \
                         ki_958, ki_959, mi_1487, mi_1488, mi_1489, mi_1490, \
                         mi_1491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1207[k] = -7.0 * ki_955[k]
                    + f_0 * mi_1487[k];

        t_1208[k] = -7.0 * ki_956[k]
                    + f_0 * mi_1488[k];

        t_1209[k] = -7.0 * ki_957[k]
                    + f_0 * mi_1489[k];

        t_1210[k] = -7.0 * ki_958[k]
                    + f_0 * mi_1490[k];

        t_1211[k] = -7.0 * ki_959[k]
                    + f_0 * mi_1491[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, t_1215, t_1216, ki_960, ki_961, ki_962, \
                         ki_963, ki_964, mi_1492, mi_1493, mi_1494, mi_1495, \
                         mi_1496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = -7.0 * ki_960[k]
                    + f_0 * mi_1492[k];

        t_1213[k] = -7.0 * ki_961[k]
                    + f_0 * mi_1493[k];

        t_1214[k] = -7.0 * ki_962[k]
                    + f_0 * mi_1494[k];

        t_1215[k] = -7.0 * ki_963[k]
                    + f_0 * mi_1495[k];

        t_1216[k] = -7.0 * ki_964[k]
                    + f_0 * mi_1496[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, t_1220, t_1221, ki_965, ki_966, ki_967, \
                         ki_968, ki_969, mi_1497, mi_1498, mi_1499, mi_1500, \
                         mi_1501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = -7.0 * ki_965[k]
                    + f_0 * mi_1497[k];

        t_1218[k] = -7.0 * ki_966[k]
                    + f_0 * mi_1498[k];

        t_1219[k] = -7.0 * ki_967[k]
                    + f_0 * mi_1499[k];

        t_1220[k] = -7.0 * ki_968[k]
                    + f_0 * mi_1500[k];

        t_1221[k] = -7.0 * ki_969[k]
                    + f_0 * mi_1501[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, t_1226, ki_970, ki_971, ki_972, \
                         ki_973, ki_974, mi_1502, mi_1503, mi_1504, mi_1505, \
                         mi_1506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = -7.0 * ki_970[k]
                    + f_0 * mi_1502[k];

        t_1223[k] = -7.0 * ki_971[k]
                    + f_0 * mi_1503[k];

        t_1224[k] = -7.0 * ki_972[k]
                    + f_0 * mi_1504[k];

        t_1225[k] = -7.0 * ki_973[k]
                    + f_0 * mi_1505[k];

        t_1226[k] = -7.0 * ki_974[k]
                    + f_0 * mi_1506[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, t_1230, t_1231, ki_975, ki_976, ki_977, \
                         ki_978, ki_979, mi_1507, mi_1508, mi_1509, mi_1510, \
                         mi_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = -7.0 * ki_975[k]
                    + f_0 * mi_1507[k];

        t_1228[k] = -7.0 * ki_976[k]
                    + f_0 * mi_1508[k];

        t_1229[k] = -7.0 * ki_977[k]
                    + f_0 * mi_1509[k];

        t_1230[k] = -7.0 * ki_978[k]
                    + f_0 * mi_1510[k];

        t_1231[k] = -7.0 * ki_979[k]
                    + f_0 * mi_1511[k];
    }

#pragma omp simd aligned(t_1232, t_1233, t_1234, t_1235, t_1236, ki_980, ki_981, ki_982, \
                         ki_983, ki_984, mi_1512, mi_1513, mi_1514, mi_1515, \
                         mi_1516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1232[k] = -8.0 * ki_980[k]
                    + f_0 * mi_1512[k];

        t_1233[k] = -8.0 * ki_981[k]
                    + f_0 * mi_1513[k];

        t_1234[k] = -8.0 * ki_982[k]
                    + f_0 * mi_1514[k];

        t_1235[k] = -8.0 * ki_983[k]
                    + f_0 * mi_1515[k];

        t_1236[k] = -8.0 * ki_984[k]
                    + f_0 * mi_1516[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, t_1241, ki_985, ki_986, ki_987, \
                         ki_988, ki_989, mi_1517, mi_1518, mi_1519, mi_1520, \
                         mi_1521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = -8.0 * ki_985[k]
                    + f_0 * mi_1517[k];

        t_1238[k] = -8.0 * ki_986[k]
                    + f_0 * mi_1518[k];

        t_1239[k] = -8.0 * ki_987[k]
                    + f_0 * mi_1519[k];

        t_1240[k] = -8.0 * ki_988[k]
                    + f_0 * mi_1520[k];

        t_1241[k] = -8.0 * ki_989[k]
                    + f_0 * mi_1521[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, t_1245, t_1246, ki_990, ki_991, ki_992, \
                         ki_993, ki_994, mi_1522, mi_1523, mi_1524, mi_1525, \
                         mi_1526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = -8.0 * ki_990[k]
                    + f_0 * mi_1522[k];

        t_1243[k] = -8.0 * ki_991[k]
                    + f_0 * mi_1523[k];

        t_1244[k] = -8.0 * ki_992[k]
                    + f_0 * mi_1524[k];

        t_1245[k] = -8.0 * ki_993[k]
                    + f_0 * mi_1525[k];

        t_1246[k] = -8.0 * ki_994[k]
                    + f_0 * mi_1526[k];
    }

#pragma omp simd aligned(t_1247, t_1248, t_1249, t_1250, t_1251, ki_995, ki_996, ki_997, \
                         ki_998, ki_999, mi_1527, mi_1528, mi_1529, mi_1530, \
                         mi_1531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1247[k] = -8.0 * ki_995[k]
                    + f_0 * mi_1527[k];

        t_1248[k] = -8.0 * ki_996[k]
                    + f_0 * mi_1528[k];

        t_1249[k] = -8.0 * ki_997[k]
                    + f_0 * mi_1529[k];

        t_1250[k] = -8.0 * ki_998[k]
                    + f_0 * mi_1530[k];

        t_1251[k] = -8.0 * ki_999[k]
                    + f_0 * mi_1531[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, t_1256, ki_1000, ki_1001, ki_1002, \
                         ki_1003, ki_1004, mi_1532, mi_1533, mi_1534, mi_1535, \
                         mi_1536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = -8.0 * ki_1000[k]
                    + f_0 * mi_1532[k];

        t_1253[k] = -8.0 * ki_1001[k]
                    + f_0 * mi_1533[k];

        t_1254[k] = -8.0 * ki_1002[k]
                    + f_0 * mi_1534[k];

        t_1255[k] = -8.0 * ki_1003[k]
                    + f_0 * mi_1535[k];

        t_1256[k] = -8.0 * ki_1004[k]
                    + f_0 * mi_1536[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, ki_1005, ki_1006, ki_1007, mi_1537, mi_1538, \
                         mi_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = -8.0 * ki_1005[k]
                    + f_0 * mi_1537[k];

        t_1258[k] = -8.0 * ki_1006[k]
                    + f_0 * mi_1538[k];

        t_1259[k] = -8.0 * ki_1007[k]
                    + f_0 * mi_1539[k];
    }
}

auto
compute_prim_geom_10_li_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t ki, const size_t mi,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_li_electron_repulsion_2_piece0(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_2_piece1(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_2_piece2(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_2_piece3(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_2_piece4(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_2_piece5(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_2_piece6(buffer, target, ki, mi, ncols, alpha);

    compute_prim_geom_10_li_electron_repulsion_2_piece7(buffer, target, ki, mi, ncols, alpha);
}

}  // namespace simdt2ceri
