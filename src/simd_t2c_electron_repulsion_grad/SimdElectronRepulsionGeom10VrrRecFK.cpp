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


#include "SimdElectronRepulsionGeom10VrrRecFK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_fk_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dk, const size_t gk,
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

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_6 = buffer.data(dk + 6);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_9 = buffer.data(dk + 9);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_13 = buffer.data(dk + 13);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_22 = buffer.data(dk + 22);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_25 = buffer.data(dk + 25);
    const auto *dk_26 = buffer.data(dk + 26);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_29 = buffer.data(dk + 29);
    const auto *dk_30 = buffer.data(dk + 30);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_33 = buffer.data(dk + 33);
    const auto *dk_34 = buffer.data(dk + 34);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_36 = buffer.data(dk + 36);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_49 = buffer.data(dk + 49);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_53 = buffer.data(dk + 53);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_60 = buffer.data(dk + 60);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_62 = buffer.data(dk + 62);
    const auto *dk_63 = buffer.data(dk + 63);
    const auto *dk_64 = buffer.data(dk + 64);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_71 = buffer.data(dk + 71);
    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_76 = buffer.data(dk + 76);
    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_80 = buffer.data(dk + 80);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_83 = buffer.data(dk + 83);
    const auto *dk_84 = buffer.data(dk + 84);
    const auto *dk_85 = buffer.data(dk + 85);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_87 = buffer.data(dk + 87);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_89 = buffer.data(dk + 89);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_94 = buffer.data(dk + 94);
    const auto *dk_95 = buffer.data(dk + 95);
    const auto *dk_96 = buffer.data(dk + 96);
    const auto *dk_97 = buffer.data(dk + 97);
    const auto *dk_98 = buffer.data(dk + 98);
    const auto *dk_99 = buffer.data(dk + 99);
    const auto *dk_100 = buffer.data(dk + 100);
    const auto *dk_101 = buffer.data(dk + 101);
    const auto *dk_102 = buffer.data(dk + 102);
    const auto *dk_103 = buffer.data(dk + 103);
    const auto *dk_104 = buffer.data(dk + 104);
    const auto *dk_105 = buffer.data(dk + 105);
    const auto *dk_106 = buffer.data(dk + 106);
    const auto *dk_107 = buffer.data(dk + 107);
    const auto *dk_108 = buffer.data(dk + 108);
    const auto *dk_109 = buffer.data(dk + 109);
    const auto *dk_110 = buffer.data(dk + 110);
    const auto *dk_111 = buffer.data(dk + 111);
    const auto *dk_112 = buffer.data(dk + 112);
    const auto *dk_113 = buffer.data(dk + 113);
    const auto *dk_114 = buffer.data(dk + 114);
    const auto *dk_115 = buffer.data(dk + 115);
    const auto *dk_116 = buffer.data(dk + 116);
    const auto *dk_117 = buffer.data(dk + 117);
    const auto *dk_118 = buffer.data(dk + 118);
    const auto *dk_119 = buffer.data(dk + 119);
    const auto *dk_120 = buffer.data(dk + 120);
    const auto *dk_121 = buffer.data(dk + 121);
    const auto *dk_122 = buffer.data(dk + 122);
    const auto *dk_123 = buffer.data(dk + 123);
    const auto *dk_124 = buffer.data(dk + 124);
    const auto *dk_125 = buffer.data(dk + 125);
    const auto *dk_126 = buffer.data(dk + 126);
    const auto *dk_127 = buffer.data(dk + 127);
    const auto *dk_128 = buffer.data(dk + 128);
    const auto *dk_129 = buffer.data(dk + 129);
    const auto *dk_130 = buffer.data(dk + 130);
    const auto *dk_131 = buffer.data(dk + 131);
    const auto *dk_132 = buffer.data(dk + 132);
    const auto *dk_133 = buffer.data(dk + 133);
    const auto *dk_134 = buffer.data(dk + 134);
    const auto *dk_135 = buffer.data(dk + 135);
    const auto *dk_136 = buffer.data(dk + 136);
    const auto *dk_137 = buffer.data(dk + 137);
    const auto *dk_138 = buffer.data(dk + 138);
    const auto *dk_139 = buffer.data(dk + 139);
    const auto *dk_140 = buffer.data(dk + 140);
    const auto *dk_141 = buffer.data(dk + 141);
    const auto *dk_142 = buffer.data(dk + 142);
    const auto *dk_143 = buffer.data(dk + 143);
    const auto *dk_144 = buffer.data(dk + 144);
    const auto *dk_145 = buffer.data(dk + 145);
    const auto *dk_146 = buffer.data(dk + 146);
    const auto *dk_147 = buffer.data(dk + 147);
    const auto *dk_148 = buffer.data(dk + 148);
    const auto *dk_149 = buffer.data(dk + 149);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_1 = buffer.data(gk + 1);
    const auto *gk_2 = buffer.data(gk + 2);
    const auto *gk_3 = buffer.data(gk + 3);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_9 = buffer.data(gk + 9);
    const auto *gk_10 = buffer.data(gk + 10);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_16 = buffer.data(gk + 16);
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_22 = buffer.data(gk + 22);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_28 = buffer.data(gk + 28);
    const auto *gk_29 = buffer.data(gk + 29);
    const auto *gk_30 = buffer.data(gk + 30);
    const auto *gk_31 = buffer.data(gk + 31);
    const auto *gk_32 = buffer.data(gk + 32);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_34 = buffer.data(gk + 34);
    const auto *gk_35 = buffer.data(gk + 35);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_41 = buffer.data(gk + 41);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_43 = buffer.data(gk + 43);
    const auto *gk_44 = buffer.data(gk + 44);
    const auto *gk_45 = buffer.data(gk + 45);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_47 = buffer.data(gk + 47);
    const auto *gk_48 = buffer.data(gk + 48);
    const auto *gk_49 = buffer.data(gk + 49);
    const auto *gk_50 = buffer.data(gk + 50);
    const auto *gk_51 = buffer.data(gk + 51);
    const auto *gk_52 = buffer.data(gk + 52);
    const auto *gk_53 = buffer.data(gk + 53);
    const auto *gk_54 = buffer.data(gk + 54);
    const auto *gk_55 = buffer.data(gk + 55);
    const auto *gk_56 = buffer.data(gk + 56);
    const auto *gk_57 = buffer.data(gk + 57);
    const auto *gk_58 = buffer.data(gk + 58);
    const auto *gk_59 = buffer.data(gk + 59);
    const auto *gk_60 = buffer.data(gk + 60);
    const auto *gk_61 = buffer.data(gk + 61);
    const auto *gk_62 = buffer.data(gk + 62);
    const auto *gk_63 = buffer.data(gk + 63);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_65 = buffer.data(gk + 65);
    const auto *gk_66 = buffer.data(gk + 66);
    const auto *gk_67 = buffer.data(gk + 67);
    const auto *gk_68 = buffer.data(gk + 68);
    const auto *gk_69 = buffer.data(gk + 69);
    const auto *gk_70 = buffer.data(gk + 70);
    const auto *gk_71 = buffer.data(gk + 71);
    const auto *gk_72 = buffer.data(gk + 72);
    const auto *gk_73 = buffer.data(gk + 73);
    const auto *gk_74 = buffer.data(gk + 74);
    const auto *gk_75 = buffer.data(gk + 75);
    const auto *gk_76 = buffer.data(gk + 76);
    const auto *gk_77 = buffer.data(gk + 77);
    const auto *gk_78 = buffer.data(gk + 78);
    const auto *gk_79 = buffer.data(gk + 79);
    const auto *gk_80 = buffer.data(gk + 80);
    const auto *gk_81 = buffer.data(gk + 81);
    const auto *gk_82 = buffer.data(gk + 82);
    const auto *gk_83 = buffer.data(gk + 83);
    const auto *gk_84 = buffer.data(gk + 84);
    const auto *gk_85 = buffer.data(gk + 85);
    const auto *gk_86 = buffer.data(gk + 86);
    const auto *gk_87 = buffer.data(gk + 87);
    const auto *gk_88 = buffer.data(gk + 88);
    const auto *gk_89 = buffer.data(gk + 89);
    const auto *gk_90 = buffer.data(gk + 90);
    const auto *gk_91 = buffer.data(gk + 91);
    const auto *gk_92 = buffer.data(gk + 92);
    const auto *gk_93 = buffer.data(gk + 93);
    const auto *gk_94 = buffer.data(gk + 94);
    const auto *gk_95 = buffer.data(gk + 95);
    const auto *gk_96 = buffer.data(gk + 96);
    const auto *gk_97 = buffer.data(gk + 97);
    const auto *gk_98 = buffer.data(gk + 98);
    const auto *gk_99 = buffer.data(gk + 99);
    const auto *gk_100 = buffer.data(gk + 100);
    const auto *gk_101 = buffer.data(gk + 101);
    const auto *gk_102 = buffer.data(gk + 102);
    const auto *gk_103 = buffer.data(gk + 103);
    const auto *gk_104 = buffer.data(gk + 104);
    const auto *gk_105 = buffer.data(gk + 105);
    const auto *gk_106 = buffer.data(gk + 106);
    const auto *gk_107 = buffer.data(gk + 107);
    const auto *gk_108 = buffer.data(gk + 108);
    const auto *gk_109 = buffer.data(gk + 109);
    const auto *gk_110 = buffer.data(gk + 110);
    const auto *gk_111 = buffer.data(gk + 111);
    const auto *gk_112 = buffer.data(gk + 112);
    const auto *gk_113 = buffer.data(gk + 113);
    const auto *gk_114 = buffer.data(gk + 114);
    const auto *gk_115 = buffer.data(gk + 115);
    const auto *gk_116 = buffer.data(gk + 116);
    const auto *gk_117 = buffer.data(gk + 117);
    const auto *gk_118 = buffer.data(gk + 118);
    const auto *gk_119 = buffer.data(gk + 119);
    const auto *gk_120 = buffer.data(gk + 120);
    const auto *gk_121 = buffer.data(gk + 121);
    const auto *gk_122 = buffer.data(gk + 122);
    const auto *gk_123 = buffer.data(gk + 123);
    const auto *gk_124 = buffer.data(gk + 124);
    const auto *gk_125 = buffer.data(gk + 125);
    const auto *gk_126 = buffer.data(gk + 126);
    const auto *gk_127 = buffer.data(gk + 127);
    const auto *gk_128 = buffer.data(gk + 128);
    const auto *gk_129 = buffer.data(gk + 129);
    const auto *gk_130 = buffer.data(gk + 130);
    const auto *gk_131 = buffer.data(gk + 131);
    const auto *gk_132 = buffer.data(gk + 132);
    const auto *gk_133 = buffer.data(gk + 133);
    const auto *gk_134 = buffer.data(gk + 134);
    const auto *gk_135 = buffer.data(gk + 135);
    const auto *gk_136 = buffer.data(gk + 136);
    const auto *gk_137 = buffer.data(gk + 137);
    const auto *gk_138 = buffer.data(gk + 138);
    const auto *gk_139 = buffer.data(gk + 139);
    const auto *gk_140 = buffer.data(gk + 140);
    const auto *gk_141 = buffer.data(gk + 141);
    const auto *gk_142 = buffer.data(gk + 142);
    const auto *gk_143 = buffer.data(gk + 143);
    const auto *gk_144 = buffer.data(gk + 144);
    const auto *gk_145 = buffer.data(gk + 145);
    const auto *gk_146 = buffer.data(gk + 146);
    const auto *gk_147 = buffer.data(gk + 147);
    const auto *gk_148 = buffer.data(gk + 148);
    const auto *gk_149 = buffer.data(gk + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, dk_0, dk_1, dk_2, dk_3, dk_4, gk_0, gk_1, \
                         gk_2, gk_3, gk_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -3.0 * dk_0[k]
                 + f_0 * gk_0[k];

        t_1[k] = -3.0 * dk_1[k]
                 + f_0 * gk_1[k];

        t_2[k] = -3.0 * dk_2[k]
                 + f_0 * gk_2[k];

        t_3[k] = -3.0 * dk_3[k]
                 + f_0 * gk_3[k];

        t_4[k] = -3.0 * dk_4[k]
                 + f_0 * gk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, dk_5, dk_6, dk_7, dk_8, dk_9, gk_5, gk_6, \
                         gk_7, gk_8, gk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -3.0 * dk_5[k]
                 + f_0 * gk_5[k];

        t_6[k] = -3.0 * dk_6[k]
                 + f_0 * gk_6[k];

        t_7[k] = -3.0 * dk_7[k]
                 + f_0 * gk_7[k];

        t_8[k] = -3.0 * dk_8[k]
                 + f_0 * gk_8[k];

        t_9[k] = -3.0 * dk_9[k]
                 + f_0 * gk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, dk_10, dk_11, dk_12, dk_13, dk_14, \
                         gk_10, gk_11, gk_12, gk_13, gk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -3.0 * dk_10[k]
                  + f_0 * gk_10[k];

        t_11[k] = -3.0 * dk_11[k]
                  + f_0 * gk_11[k];

        t_12[k] = -3.0 * dk_12[k]
                  + f_0 * gk_12[k];

        t_13[k] = -3.0 * dk_13[k]
                  + f_0 * gk_13[k];

        t_14[k] = -3.0 * dk_14[k]
                  + f_0 * gk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, dk_15, dk_16, dk_17, dk_18, dk_19, \
                         gk_15, gk_16, gk_17, gk_18, gk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -3.0 * dk_15[k]
                  + f_0 * gk_15[k];

        t_16[k] = -3.0 * dk_16[k]
                  + f_0 * gk_16[k];

        t_17[k] = -3.0 * dk_17[k]
                  + f_0 * gk_17[k];

        t_18[k] = -3.0 * dk_18[k]
                  + f_0 * gk_18[k];

        t_19[k] = -3.0 * dk_19[k]
                  + f_0 * gk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, dk_20, dk_21, dk_22, dk_23, dk_24, \
                         gk_20, gk_21, gk_22, gk_23, gk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -3.0 * dk_20[k]
                  + f_0 * gk_20[k];

        t_21[k] = -3.0 * dk_21[k]
                  + f_0 * gk_21[k];

        t_22[k] = -3.0 * dk_22[k]
                  + f_0 * gk_22[k];

        t_23[k] = -3.0 * dk_23[k]
                  + f_0 * gk_23[k];

        t_24[k] = -3.0 * dk_24[k]
                  + f_0 * gk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, dk_25, dk_26, dk_27, dk_28, dk_29, \
                         gk_25, gk_26, gk_27, gk_28, gk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -3.0 * dk_25[k]
                  + f_0 * gk_25[k];

        t_26[k] = -3.0 * dk_26[k]
                  + f_0 * gk_26[k];

        t_27[k] = -3.0 * dk_27[k]
                  + f_0 * gk_27[k];

        t_28[k] = -3.0 * dk_28[k]
                  + f_0 * gk_28[k];

        t_29[k] = -3.0 * dk_29[k]
                  + f_0 * gk_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, dk_30, dk_31, dk_32, dk_33, dk_34, \
                         gk_30, gk_31, gk_32, gk_33, gk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -3.0 * dk_30[k]
                  + f_0 * gk_30[k];

        t_31[k] = -3.0 * dk_31[k]
                  + f_0 * gk_31[k];

        t_32[k] = -3.0 * dk_32[k]
                  + f_0 * gk_32[k];

        t_33[k] = -3.0 * dk_33[k]
                  + f_0 * gk_33[k];

        t_34[k] = -3.0 * dk_34[k]
                  + f_0 * gk_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, dk_35, dk_36, dk_37, dk_38, dk_39, \
                         gk_35, gk_36, gk_37, gk_38, gk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -3.0 * dk_35[k]
                  + f_0 * gk_35[k];

        t_36[k] = -2.0 * dk_36[k]
                  + f_0 * gk_36[k];

        t_37[k] = -2.0 * dk_37[k]
                  + f_0 * gk_37[k];

        t_38[k] = -2.0 * dk_38[k]
                  + f_0 * gk_38[k];

        t_39[k] = -2.0 * dk_39[k]
                  + f_0 * gk_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, dk_40, dk_41, dk_42, dk_43, dk_44, \
                         gk_40, gk_41, gk_42, gk_43, gk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -2.0 * dk_40[k]
                  + f_0 * gk_40[k];

        t_41[k] = -2.0 * dk_41[k]
                  + f_0 * gk_41[k];

        t_42[k] = -2.0 * dk_42[k]
                  + f_0 * gk_42[k];

        t_43[k] = -2.0 * dk_43[k]
                  + f_0 * gk_43[k];

        t_44[k] = -2.0 * dk_44[k]
                  + f_0 * gk_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, dk_45, dk_46, dk_47, dk_48, dk_49, \
                         gk_45, gk_46, gk_47, gk_48, gk_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * dk_45[k]
                  + f_0 * gk_45[k];

        t_46[k] = -2.0 * dk_46[k]
                  + f_0 * gk_46[k];

        t_47[k] = -2.0 * dk_47[k]
                  + f_0 * gk_47[k];

        t_48[k] = -2.0 * dk_48[k]
                  + f_0 * gk_48[k];

        t_49[k] = -2.0 * dk_49[k]
                  + f_0 * gk_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, dk_50, dk_51, dk_52, dk_53, dk_54, \
                         gk_50, gk_51, gk_52, gk_53, gk_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * dk_50[k]
                  + f_0 * gk_50[k];

        t_51[k] = -2.0 * dk_51[k]
                  + f_0 * gk_51[k];

        t_52[k] = -2.0 * dk_52[k]
                  + f_0 * gk_52[k];

        t_53[k] = -2.0 * dk_53[k]
                  + f_0 * gk_53[k];

        t_54[k] = -2.0 * dk_54[k]
                  + f_0 * gk_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, dk_55, dk_56, dk_57, dk_58, dk_59, \
                         gk_55, gk_56, gk_57, gk_58, gk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * dk_55[k]
                  + f_0 * gk_55[k];

        t_56[k] = -2.0 * dk_56[k]
                  + f_0 * gk_56[k];

        t_57[k] = -2.0 * dk_57[k]
                  + f_0 * gk_57[k];

        t_58[k] = -2.0 * dk_58[k]
                  + f_0 * gk_58[k];

        t_59[k] = -2.0 * dk_59[k]
                  + f_0 * gk_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, dk_60, dk_61, dk_62, dk_63, dk_64, \
                         gk_60, gk_61, gk_62, gk_63, gk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -2.0 * dk_60[k]
                  + f_0 * gk_60[k];

        t_61[k] = -2.0 * dk_61[k]
                  + f_0 * gk_61[k];

        t_62[k] = -2.0 * dk_62[k]
                  + f_0 * gk_62[k];

        t_63[k] = -2.0 * dk_63[k]
                  + f_0 * gk_63[k];

        t_64[k] = -2.0 * dk_64[k]
                  + f_0 * gk_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, dk_65, dk_66, dk_67, dk_68, dk_69, \
                         gk_65, gk_66, gk_67, gk_68, gk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -2.0 * dk_65[k]
                  + f_0 * gk_65[k];

        t_66[k] = -2.0 * dk_66[k]
                  + f_0 * gk_66[k];

        t_67[k] = -2.0 * dk_67[k]
                  + f_0 * gk_67[k];

        t_68[k] = -2.0 * dk_68[k]
                  + f_0 * gk_68[k];

        t_69[k] = -2.0 * dk_69[k]
                  + f_0 * gk_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, dk_70, dk_71, dk_72, dk_73, dk_74, \
                         gk_70, gk_71, gk_72, gk_73, gk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * dk_70[k]
                  + f_0 * gk_70[k];

        t_71[k] = -2.0 * dk_71[k]
                  + f_0 * gk_71[k];

        t_72[k] = -2.0 * dk_72[k]
                  + f_0 * gk_72[k];

        t_73[k] = -2.0 * dk_73[k]
                  + f_0 * gk_73[k];

        t_74[k] = -2.0 * dk_74[k]
                  + f_0 * gk_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, dk_75, dk_76, dk_77, dk_78, dk_79, \
                         gk_75, gk_76, gk_77, gk_78, gk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * dk_75[k]
                  + f_0 * gk_75[k];

        t_76[k] = -2.0 * dk_76[k]
                  + f_0 * gk_76[k];

        t_77[k] = -2.0 * dk_77[k]
                  + f_0 * gk_77[k];

        t_78[k] = -2.0 * dk_78[k]
                  + f_0 * gk_78[k];

        t_79[k] = -2.0 * dk_79[k]
                  + f_0 * gk_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, dk_80, dk_81, dk_82, dk_83, dk_84, \
                         gk_80, gk_81, gk_82, gk_83, gk_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -2.0 * dk_80[k]
                  + f_0 * gk_80[k];

        t_81[k] = -2.0 * dk_81[k]
                  + f_0 * gk_81[k];

        t_82[k] = -2.0 * dk_82[k]
                  + f_0 * gk_82[k];

        t_83[k] = -2.0 * dk_83[k]
                  + f_0 * gk_83[k];

        t_84[k] = -2.0 * dk_84[k]
                  + f_0 * gk_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, dk_85, dk_86, dk_87, dk_88, dk_89, \
                         gk_85, gk_86, gk_87, gk_88, gk_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * dk_85[k]
                  + f_0 * gk_85[k];

        t_86[k] = -2.0 * dk_86[k]
                  + f_0 * gk_86[k];

        t_87[k] = -2.0 * dk_87[k]
                  + f_0 * gk_87[k];

        t_88[k] = -2.0 * dk_88[k]
                  + f_0 * gk_88[k];

        t_89[k] = -2.0 * dk_89[k]
                  + f_0 * gk_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, dk_90, dk_91, dk_92, dk_93, dk_94, \
                         gk_90, gk_91, gk_92, gk_93, gk_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * dk_90[k]
                  + f_0 * gk_90[k];

        t_91[k] = -2.0 * dk_91[k]
                  + f_0 * gk_91[k];

        t_92[k] = -2.0 * dk_92[k]
                  + f_0 * gk_92[k];

        t_93[k] = -2.0 * dk_93[k]
                  + f_0 * gk_93[k];

        t_94[k] = -2.0 * dk_94[k]
                  + f_0 * gk_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, dk_95, dk_96, dk_97, dk_98, dk_99, \
                         gk_95, gk_96, gk_97, gk_98, gk_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * dk_95[k]
                  + f_0 * gk_95[k];

        t_96[k] = -2.0 * dk_96[k]
                  + f_0 * gk_96[k];

        t_97[k] = -2.0 * dk_97[k]
                  + f_0 * gk_97[k];

        t_98[k] = -2.0 * dk_98[k]
                  + f_0 * gk_98[k];

        t_99[k] = -2.0 * dk_99[k]
                  + f_0 * gk_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, dk_100, dk_101, dk_102, dk_103, \
                         dk_104, gk_100, gk_101, gk_102, gk_103, \
                         gk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * dk_100[k]
                   + f_0 * gk_100[k];

        t_101[k] = -2.0 * dk_101[k]
                   + f_0 * gk_101[k];

        t_102[k] = -2.0 * dk_102[k]
                   + f_0 * gk_102[k];

        t_103[k] = -2.0 * dk_103[k]
                   + f_0 * gk_103[k];

        t_104[k] = -2.0 * dk_104[k]
                   + f_0 * gk_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, dk_105, dk_106, dk_107, dk_108, \
                         dk_109, gk_105, gk_106, gk_107, gk_108, \
                         gk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * dk_105[k]
                   + f_0 * gk_105[k];

        t_106[k] = -2.0 * dk_106[k]
                   + f_0 * gk_106[k];

        t_107[k] = -2.0 * dk_107[k]
                   + f_0 * gk_107[k];

        t_108[k] = -dk_108[k]
                   + f_0 * gk_108[k];

        t_109[k] = -dk_109[k]
                   + f_0 * gk_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, dk_110, dk_111, dk_112, dk_113, \
                         dk_114, gk_110, gk_111, gk_112, gk_113, \
                         gk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -dk_110[k]
                   + f_0 * gk_110[k];

        t_111[k] = -dk_111[k]
                   + f_0 * gk_111[k];

        t_112[k] = -dk_112[k]
                   + f_0 * gk_112[k];

        t_113[k] = -dk_113[k]
                   + f_0 * gk_113[k];

        t_114[k] = -dk_114[k]
                   + f_0 * gk_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, dk_115, dk_116, dk_117, dk_118, \
                         dk_119, gk_115, gk_116, gk_117, gk_118, \
                         gk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -dk_115[k]
                   + f_0 * gk_115[k];

        t_116[k] = -dk_116[k]
                   + f_0 * gk_116[k];

        t_117[k] = -dk_117[k]
                   + f_0 * gk_117[k];

        t_118[k] = -dk_118[k]
                   + f_0 * gk_118[k];

        t_119[k] = -dk_119[k]
                   + f_0 * gk_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, dk_120, dk_121, dk_122, dk_123, \
                         dk_124, gk_120, gk_121, gk_122, gk_123, \
                         gk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -dk_120[k]
                   + f_0 * gk_120[k];

        t_121[k] = -dk_121[k]
                   + f_0 * gk_121[k];

        t_122[k] = -dk_122[k]
                   + f_0 * gk_122[k];

        t_123[k] = -dk_123[k]
                   + f_0 * gk_123[k];

        t_124[k] = -dk_124[k]
                   + f_0 * gk_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, dk_125, dk_126, dk_127, dk_128, \
                         dk_129, gk_125, gk_126, gk_127, gk_128, \
                         gk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -dk_125[k]
                   + f_0 * gk_125[k];

        t_126[k] = -dk_126[k]
                   + f_0 * gk_126[k];

        t_127[k] = -dk_127[k]
                   + f_0 * gk_127[k];

        t_128[k] = -dk_128[k]
                   + f_0 * gk_128[k];

        t_129[k] = -dk_129[k]
                   + f_0 * gk_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, dk_130, dk_131, dk_132, dk_133, \
                         dk_134, gk_130, gk_131, gk_132, gk_133, \
                         gk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -dk_130[k]
                   + f_0 * gk_130[k];

        t_131[k] = -dk_131[k]
                   + f_0 * gk_131[k];

        t_132[k] = -dk_132[k]
                   + f_0 * gk_132[k];

        t_133[k] = -dk_133[k]
                   + f_0 * gk_133[k];

        t_134[k] = -dk_134[k]
                   + f_0 * gk_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, dk_135, dk_136, dk_137, dk_138, \
                         dk_139, gk_135, gk_136, gk_137, gk_138, \
                         gk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -dk_135[k]
                   + f_0 * gk_135[k];

        t_136[k] = -dk_136[k]
                   + f_0 * gk_136[k];

        t_137[k] = -dk_137[k]
                   + f_0 * gk_137[k];

        t_138[k] = -dk_138[k]
                   + f_0 * gk_138[k];

        t_139[k] = -dk_139[k]
                   + f_0 * gk_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, dk_140, dk_141, dk_142, dk_143, \
                         dk_144, gk_140, gk_141, gk_142, gk_143, \
                         gk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -dk_140[k]
                   + f_0 * gk_140[k];

        t_141[k] = -dk_141[k]
                   + f_0 * gk_141[k];

        t_142[k] = -dk_142[k]
                   + f_0 * gk_142[k];

        t_143[k] = -dk_143[k]
                   + f_0 * gk_143[k];

        t_144[k] = -dk_144[k]
                   + f_0 * gk_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, dk_145, dk_146, dk_147, dk_148, \
                         dk_149, gk_145, gk_146, gk_147, gk_148, \
                         gk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -dk_145[k]
                   + f_0 * gk_145[k];

        t_146[k] = -dk_146[k]
                   + f_0 * gk_146[k];

        t_147[k] = -dk_147[k]
                   + f_0 * gk_147[k];

        t_148[k] = -dk_148[k]
                   + f_0 * gk_148[k];

        t_149[k] = -dk_149[k]
                   + f_0 * gk_149[k];
    }
}

static auto
compute_prim_geom_10_fk_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dk, const size_t gk,
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

    const auto *dk_150 = buffer.data(dk + 150);
    const auto *dk_151 = buffer.data(dk + 151);
    const auto *dk_152 = buffer.data(dk + 152);
    const auto *dk_153 = buffer.data(dk + 153);
    const auto *dk_154 = buffer.data(dk + 154);
    const auto *dk_155 = buffer.data(dk + 155);
    const auto *dk_156 = buffer.data(dk + 156);
    const auto *dk_157 = buffer.data(dk + 157);
    const auto *dk_158 = buffer.data(dk + 158);
    const auto *dk_159 = buffer.data(dk + 159);
    const auto *dk_160 = buffer.data(dk + 160);
    const auto *dk_161 = buffer.data(dk + 161);
    const auto *dk_162 = buffer.data(dk + 162);
    const auto *dk_163 = buffer.data(dk + 163);
    const auto *dk_164 = buffer.data(dk + 164);
    const auto *dk_165 = buffer.data(dk + 165);
    const auto *dk_166 = buffer.data(dk + 166);
    const auto *dk_167 = buffer.data(dk + 167);
    const auto *dk_168 = buffer.data(dk + 168);
    const auto *dk_169 = buffer.data(dk + 169);
    const auto *dk_170 = buffer.data(dk + 170);
    const auto *dk_171 = buffer.data(dk + 171);
    const auto *dk_172 = buffer.data(dk + 172);
    const auto *dk_173 = buffer.data(dk + 173);
    const auto *dk_174 = buffer.data(dk + 174);
    const auto *dk_175 = buffer.data(dk + 175);
    const auto *dk_176 = buffer.data(dk + 176);
    const auto *dk_177 = buffer.data(dk + 177);
    const auto *dk_178 = buffer.data(dk + 178);
    const auto *dk_179 = buffer.data(dk + 179);
    const auto *dk_180 = buffer.data(dk + 180);
    const auto *dk_181 = buffer.data(dk + 181);
    const auto *dk_182 = buffer.data(dk + 182);
    const auto *dk_183 = buffer.data(dk + 183);
    const auto *dk_184 = buffer.data(dk + 184);
    const auto *dk_185 = buffer.data(dk + 185);
    const auto *dk_186 = buffer.data(dk + 186);
    const auto *dk_187 = buffer.data(dk + 187);
    const auto *dk_188 = buffer.data(dk + 188);
    const auto *dk_189 = buffer.data(dk + 189);
    const auto *dk_190 = buffer.data(dk + 190);
    const auto *dk_191 = buffer.data(dk + 191);
    const auto *dk_192 = buffer.data(dk + 192);
    const auto *dk_193 = buffer.data(dk + 193);
    const auto *dk_194 = buffer.data(dk + 194);
    const auto *dk_195 = buffer.data(dk + 195);
    const auto *dk_196 = buffer.data(dk + 196);
    const auto *dk_197 = buffer.data(dk + 197);
    const auto *dk_198 = buffer.data(dk + 198);
    const auto *dk_199 = buffer.data(dk + 199);
    const auto *dk_200 = buffer.data(dk + 200);
    const auto *dk_201 = buffer.data(dk + 201);
    const auto *dk_202 = buffer.data(dk + 202);
    const auto *dk_203 = buffer.data(dk + 203);
    const auto *dk_204 = buffer.data(dk + 204);
    const auto *dk_205 = buffer.data(dk + 205);
    const auto *dk_206 = buffer.data(dk + 206);
    const auto *dk_207 = buffer.data(dk + 207);
    const auto *dk_208 = buffer.data(dk + 208);
    const auto *dk_209 = buffer.data(dk + 209);
    const auto *dk_210 = buffer.data(dk + 210);
    const auto *dk_211 = buffer.data(dk + 211);
    const auto *dk_212 = buffer.data(dk + 212);
    const auto *dk_213 = buffer.data(dk + 213);
    const auto *dk_214 = buffer.data(dk + 214);
    const auto *dk_215 = buffer.data(dk + 215);

    const auto *gk_150 = buffer.data(gk + 150);
    const auto *gk_151 = buffer.data(gk + 151);
    const auto *gk_152 = buffer.data(gk + 152);
    const auto *gk_153 = buffer.data(gk + 153);
    const auto *gk_154 = buffer.data(gk + 154);
    const auto *gk_155 = buffer.data(gk + 155);
    const auto *gk_156 = buffer.data(gk + 156);
    const auto *gk_157 = buffer.data(gk + 157);
    const auto *gk_158 = buffer.data(gk + 158);
    const auto *gk_159 = buffer.data(gk + 159);
    const auto *gk_160 = buffer.data(gk + 160);
    const auto *gk_161 = buffer.data(gk + 161);
    const auto *gk_162 = buffer.data(gk + 162);
    const auto *gk_163 = buffer.data(gk + 163);
    const auto *gk_164 = buffer.data(gk + 164);
    const auto *gk_165 = buffer.data(gk + 165);
    const auto *gk_166 = buffer.data(gk + 166);
    const auto *gk_167 = buffer.data(gk + 167);
    const auto *gk_168 = buffer.data(gk + 168);
    const auto *gk_169 = buffer.data(gk + 169);
    const auto *gk_170 = buffer.data(gk + 170);
    const auto *gk_171 = buffer.data(gk + 171);
    const auto *gk_172 = buffer.data(gk + 172);
    const auto *gk_173 = buffer.data(gk + 173);
    const auto *gk_174 = buffer.data(gk + 174);
    const auto *gk_175 = buffer.data(gk + 175);
    const auto *gk_176 = buffer.data(gk + 176);
    const auto *gk_177 = buffer.data(gk + 177);
    const auto *gk_178 = buffer.data(gk + 178);
    const auto *gk_179 = buffer.data(gk + 179);
    const auto *gk_180 = buffer.data(gk + 180);
    const auto *gk_181 = buffer.data(gk + 181);
    const auto *gk_182 = buffer.data(gk + 182);
    const auto *gk_183 = buffer.data(gk + 183);
    const auto *gk_184 = buffer.data(gk + 184);
    const auto *gk_185 = buffer.data(gk + 185);
    const auto *gk_186 = buffer.data(gk + 186);
    const auto *gk_187 = buffer.data(gk + 187);
    const auto *gk_188 = buffer.data(gk + 188);
    const auto *gk_189 = buffer.data(gk + 189);
    const auto *gk_190 = buffer.data(gk + 190);
    const auto *gk_191 = buffer.data(gk + 191);
    const auto *gk_192 = buffer.data(gk + 192);
    const auto *gk_193 = buffer.data(gk + 193);
    const auto *gk_194 = buffer.data(gk + 194);
    const auto *gk_195 = buffer.data(gk + 195);
    const auto *gk_196 = buffer.data(gk + 196);
    const auto *gk_197 = buffer.data(gk + 197);
    const auto *gk_198 = buffer.data(gk + 198);
    const auto *gk_199 = buffer.data(gk + 199);
    const auto *gk_200 = buffer.data(gk + 200);
    const auto *gk_201 = buffer.data(gk + 201);
    const auto *gk_202 = buffer.data(gk + 202);
    const auto *gk_203 = buffer.data(gk + 203);
    const auto *gk_204 = buffer.data(gk + 204);
    const auto *gk_205 = buffer.data(gk + 205);
    const auto *gk_206 = buffer.data(gk + 206);
    const auto *gk_207 = buffer.data(gk + 207);
    const auto *gk_208 = buffer.data(gk + 208);
    const auto *gk_209 = buffer.data(gk + 209);
    const auto *gk_210 = buffer.data(gk + 210);
    const auto *gk_211 = buffer.data(gk + 211);
    const auto *gk_212 = buffer.data(gk + 212);
    const auto *gk_213 = buffer.data(gk + 213);
    const auto *gk_214 = buffer.data(gk + 214);
    const auto *gk_215 = buffer.data(gk + 215);
    const auto *gk_216 = buffer.data(gk + 216);
    const auto *gk_217 = buffer.data(gk + 217);
    const auto *gk_218 = buffer.data(gk + 218);
    const auto *gk_219 = buffer.data(gk + 219);
    const auto *gk_220 = buffer.data(gk + 220);
    const auto *gk_221 = buffer.data(gk + 221);
    const auto *gk_222 = buffer.data(gk + 222);
    const auto *gk_223 = buffer.data(gk + 223);
    const auto *gk_224 = buffer.data(gk + 224);
    const auto *gk_225 = buffer.data(gk + 225);
    const auto *gk_226 = buffer.data(gk + 226);
    const auto *gk_227 = buffer.data(gk + 227);
    const auto *gk_228 = buffer.data(gk + 228);
    const auto *gk_229 = buffer.data(gk + 229);
    const auto *gk_230 = buffer.data(gk + 230);
    const auto *gk_231 = buffer.data(gk + 231);
    const auto *gk_232 = buffer.data(gk + 232);
    const auto *gk_233 = buffer.data(gk + 233);
    const auto *gk_234 = buffer.data(gk + 234);
    const auto *gk_235 = buffer.data(gk + 235);
    const auto *gk_236 = buffer.data(gk + 236);
    const auto *gk_237 = buffer.data(gk + 237);
    const auto *gk_238 = buffer.data(gk + 238);
    const auto *gk_239 = buffer.data(gk + 239);
    const auto *gk_240 = buffer.data(gk + 240);
    const auto *gk_241 = buffer.data(gk + 241);
    const auto *gk_242 = buffer.data(gk + 242);
    const auto *gk_243 = buffer.data(gk + 243);
    const auto *gk_244 = buffer.data(gk + 244);
    const auto *gk_245 = buffer.data(gk + 245);
    const auto *gk_246 = buffer.data(gk + 246);
    const auto *gk_247 = buffer.data(gk + 247);
    const auto *gk_248 = buffer.data(gk + 248);
    const auto *gk_249 = buffer.data(gk + 249);
    const auto *gk_250 = buffer.data(gk + 250);
    const auto *gk_251 = buffer.data(gk + 251);
    const auto *gk_252 = buffer.data(gk + 252);
    const auto *gk_253 = buffer.data(gk + 253);
    const auto *gk_254 = buffer.data(gk + 254);
    const auto *gk_255 = buffer.data(gk + 255);
    const auto *gk_256 = buffer.data(gk + 256);
    const auto *gk_257 = buffer.data(gk + 257);
    const auto *gk_258 = buffer.data(gk + 258);
    const auto *gk_259 = buffer.data(gk + 259);
    const auto *gk_260 = buffer.data(gk + 260);
    const auto *gk_261 = buffer.data(gk + 261);
    const auto *gk_262 = buffer.data(gk + 262);
    const auto *gk_263 = buffer.data(gk + 263);
    const auto *gk_264 = buffer.data(gk + 264);
    const auto *gk_265 = buffer.data(gk + 265);
    const auto *gk_266 = buffer.data(gk + 266);
    const auto *gk_267 = buffer.data(gk + 267);
    const auto *gk_268 = buffer.data(gk + 268);
    const auto *gk_269 = buffer.data(gk + 269);
    const auto *gk_270 = buffer.data(gk + 270);
    const auto *gk_271 = buffer.data(gk + 271);
    const auto *gk_272 = buffer.data(gk + 272);
    const auto *gk_273 = buffer.data(gk + 273);
    const auto *gk_274 = buffer.data(gk + 274);
    const auto *gk_275 = buffer.data(gk + 275);
    const auto *gk_276 = buffer.data(gk + 276);
    const auto *gk_277 = buffer.data(gk + 277);
    const auto *gk_278 = buffer.data(gk + 278);
    const auto *gk_279 = buffer.data(gk + 279);
    const auto *gk_280 = buffer.data(gk + 280);
    const auto *gk_281 = buffer.data(gk + 281);
    const auto *gk_282 = buffer.data(gk + 282);
    const auto *gk_283 = buffer.data(gk + 283);
    const auto *gk_284 = buffer.data(gk + 284);
    const auto *gk_285 = buffer.data(gk + 285);
    const auto *gk_286 = buffer.data(gk + 286);
    const auto *gk_287 = buffer.data(gk + 287);
    const auto *gk_288 = buffer.data(gk + 288);
    const auto *gk_289 = buffer.data(gk + 289);
    const auto *gk_290 = buffer.data(gk + 290);
    const auto *gk_291 = buffer.data(gk + 291);
    const auto *gk_292 = buffer.data(gk + 292);
    const auto *gk_293 = buffer.data(gk + 293);
    const auto *gk_294 = buffer.data(gk + 294);
    const auto *gk_295 = buffer.data(gk + 295);
    const auto *gk_296 = buffer.data(gk + 296);
    const auto *gk_297 = buffer.data(gk + 297);
    const auto *gk_298 = buffer.data(gk + 298);
    const auto *gk_299 = buffer.data(gk + 299);
    const auto *gk_300 = buffer.data(gk + 300);
    const auto *gk_301 = buffer.data(gk + 301);
    const auto *gk_302 = buffer.data(gk + 302);
    const auto *gk_303 = buffer.data(gk + 303);
    const auto *gk_304 = buffer.data(gk + 304);
    const auto *gk_305 = buffer.data(gk + 305);
    const auto *gk_306 = buffer.data(gk + 306);
    const auto *gk_307 = buffer.data(gk + 307);
    const auto *gk_308 = buffer.data(gk + 308);
    const auto *gk_309 = buffer.data(gk + 309);
    const auto *gk_310 = buffer.data(gk + 310);
    const auto *gk_311 = buffer.data(gk + 311);
    const auto *gk_312 = buffer.data(gk + 312);
    const auto *gk_313 = buffer.data(gk + 313);
    const auto *gk_314 = buffer.data(gk + 314);
    const auto *gk_315 = buffer.data(gk + 315);
    const auto *gk_316 = buffer.data(gk + 316);
    const auto *gk_317 = buffer.data(gk + 317);
    const auto *gk_318 = buffer.data(gk + 318);
    const auto *gk_319 = buffer.data(gk + 319);
    const auto *gk_320 = buffer.data(gk + 320);
    const auto *gk_321 = buffer.data(gk + 321);
    const auto *gk_322 = buffer.data(gk + 322);
    const auto *gk_323 = buffer.data(gk + 323);
    const auto *gk_324 = buffer.data(gk + 324);
    const auto *gk_325 = buffer.data(gk + 325);
    const auto *gk_326 = buffer.data(gk + 326);
    const auto *gk_327 = buffer.data(gk + 327);
    const auto *gk_328 = buffer.data(gk + 328);
    const auto *gk_329 = buffer.data(gk + 329);
    const auto *gk_330 = buffer.data(gk + 330);
    const auto *gk_331 = buffer.data(gk + 331);
    const auto *gk_332 = buffer.data(gk + 332);
    const auto *gk_333 = buffer.data(gk + 333);
    const auto *gk_334 = buffer.data(gk + 334);
    const auto *gk_335 = buffer.data(gk + 335);
    const auto *gk_336 = buffer.data(gk + 336);
    const auto *gk_337 = buffer.data(gk + 337);
    const auto *gk_338 = buffer.data(gk + 338);
    const auto *gk_339 = buffer.data(gk + 339);
    const auto *gk_340 = buffer.data(gk + 340);
    const auto *gk_341 = buffer.data(gk + 341);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, dk_150, dk_151, dk_152, dk_153, \
                         dk_154, gk_150, gk_151, gk_152, gk_153, \
                         gk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -dk_150[k]
                   + f_0 * gk_150[k];

        t_151[k] = -dk_151[k]
                   + f_0 * gk_151[k];

        t_152[k] = -dk_152[k]
                   + f_0 * gk_152[k];

        t_153[k] = -dk_153[k]
                   + f_0 * gk_153[k];

        t_154[k] = -dk_154[k]
                   + f_0 * gk_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, dk_155, dk_156, dk_157, dk_158, \
                         dk_159, gk_155, gk_156, gk_157, gk_158, \
                         gk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -dk_155[k]
                   + f_0 * gk_155[k];

        t_156[k] = -dk_156[k]
                   + f_0 * gk_156[k];

        t_157[k] = -dk_157[k]
                   + f_0 * gk_157[k];

        t_158[k] = -dk_158[k]
                   + f_0 * gk_158[k];

        t_159[k] = -dk_159[k]
                   + f_0 * gk_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, dk_160, dk_161, dk_162, dk_163, \
                         dk_164, gk_160, gk_161, gk_162, gk_163, \
                         gk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -dk_160[k]
                   + f_0 * gk_160[k];

        t_161[k] = -dk_161[k]
                   + f_0 * gk_161[k];

        t_162[k] = -dk_162[k]
                   + f_0 * gk_162[k];

        t_163[k] = -dk_163[k]
                   + f_0 * gk_163[k];

        t_164[k] = -dk_164[k]
                   + f_0 * gk_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, dk_165, dk_166, dk_167, dk_168, \
                         dk_169, gk_165, gk_166, gk_167, gk_168, \
                         gk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -dk_165[k]
                   + f_0 * gk_165[k];

        t_166[k] = -dk_166[k]
                   + f_0 * gk_166[k];

        t_167[k] = -dk_167[k]
                   + f_0 * gk_167[k];

        t_168[k] = -dk_168[k]
                   + f_0 * gk_168[k];

        t_169[k] = -dk_169[k]
                   + f_0 * gk_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, dk_170, dk_171, dk_172, dk_173, \
                         dk_174, gk_170, gk_171, gk_172, gk_173, \
                         gk_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -dk_170[k]
                   + f_0 * gk_170[k];

        t_171[k] = -dk_171[k]
                   + f_0 * gk_171[k];

        t_172[k] = -dk_172[k]
                   + f_0 * gk_172[k];

        t_173[k] = -dk_173[k]
                   + f_0 * gk_173[k];

        t_174[k] = -dk_174[k]
                   + f_0 * gk_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, dk_175, dk_176, dk_177, dk_178, \
                         dk_179, gk_175, gk_176, gk_177, gk_178, \
                         gk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -dk_175[k]
                   + f_0 * gk_175[k];

        t_176[k] = -dk_176[k]
                   + f_0 * gk_176[k];

        t_177[k] = -dk_177[k]
                   + f_0 * gk_177[k];

        t_178[k] = -dk_178[k]
                   + f_0 * gk_178[k];

        t_179[k] = -dk_179[k]
                   + f_0 * gk_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, dk_180, dk_181, dk_182, dk_183, \
                         dk_184, gk_180, gk_181, gk_182, gk_183, \
                         gk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -dk_180[k]
                   + f_0 * gk_180[k];

        t_181[k] = -dk_181[k]
                   + f_0 * gk_181[k];

        t_182[k] = -dk_182[k]
                   + f_0 * gk_182[k];

        t_183[k] = -dk_183[k]
                   + f_0 * gk_183[k];

        t_184[k] = -dk_184[k]
                   + f_0 * gk_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, dk_185, dk_186, dk_187, dk_188, \
                         dk_189, gk_185, gk_186, gk_187, gk_188, \
                         gk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -dk_185[k]
                   + f_0 * gk_185[k];

        t_186[k] = -dk_186[k]
                   + f_0 * gk_186[k];

        t_187[k] = -dk_187[k]
                   + f_0 * gk_187[k];

        t_188[k] = -dk_188[k]
                   + f_0 * gk_188[k];

        t_189[k] = -dk_189[k]
                   + f_0 * gk_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, dk_190, dk_191, dk_192, dk_193, \
                         dk_194, gk_190, gk_191, gk_192, gk_193, \
                         gk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -dk_190[k]
                   + f_0 * gk_190[k];

        t_191[k] = -dk_191[k]
                   + f_0 * gk_191[k];

        t_192[k] = -dk_192[k]
                   + f_0 * gk_192[k];

        t_193[k] = -dk_193[k]
                   + f_0 * gk_193[k];

        t_194[k] = -dk_194[k]
                   + f_0 * gk_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, dk_195, dk_196, dk_197, dk_198, \
                         dk_199, gk_195, gk_196, gk_197, gk_198, \
                         gk_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -dk_195[k]
                   + f_0 * gk_195[k];

        t_196[k] = -dk_196[k]
                   + f_0 * gk_196[k];

        t_197[k] = -dk_197[k]
                   + f_0 * gk_197[k];

        t_198[k] = -dk_198[k]
                   + f_0 * gk_198[k];

        t_199[k] = -dk_199[k]
                   + f_0 * gk_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, dk_200, dk_201, dk_202, dk_203, \
                         dk_204, gk_200, gk_201, gk_202, gk_203, \
                         gk_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -dk_200[k]
                   + f_0 * gk_200[k];

        t_201[k] = -dk_201[k]
                   + f_0 * gk_201[k];

        t_202[k] = -dk_202[k]
                   + f_0 * gk_202[k];

        t_203[k] = -dk_203[k]
                   + f_0 * gk_203[k];

        t_204[k] = -dk_204[k]
                   + f_0 * gk_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, dk_205, dk_206, dk_207, dk_208, \
                         dk_209, gk_205, gk_206, gk_207, gk_208, \
                         gk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -dk_205[k]
                   + f_0 * gk_205[k];

        t_206[k] = -dk_206[k]
                   + f_0 * gk_206[k];

        t_207[k] = -dk_207[k]
                   + f_0 * gk_207[k];

        t_208[k] = -dk_208[k]
                   + f_0 * gk_208[k];

        t_209[k] = -dk_209[k]
                   + f_0 * gk_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, dk_210, dk_211, dk_212, dk_213, \
                         dk_214, gk_210, gk_211, gk_212, gk_213, \
                         gk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -dk_210[k]
                   + f_0 * gk_210[k];

        t_211[k] = -dk_211[k]
                   + f_0 * gk_211[k];

        t_212[k] = -dk_212[k]
                   + f_0 * gk_212[k];

        t_213[k] = -dk_213[k]
                   + f_0 * gk_213[k];

        t_214[k] = -dk_214[k]
                   + f_0 * gk_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, t_220, t_221, dk_215, gk_215, \
                         gk_216, gk_217, gk_218, gk_219, gk_220, \
                         gk_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -dk_215[k]
                   + f_0 * gk_215[k];

        t_216[k] = f_0 * gk_216[k];

        t_217[k] = f_0 * gk_217[k];

        t_218[k] = f_0 * gk_218[k];

        t_219[k] = f_0 * gk_219[k];

        t_220[k] = f_0 * gk_220[k];

        t_221[k] = f_0 * gk_221[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, t_228, t_229, gk_222, \
                         gk_223, gk_224, gk_225, gk_226, gk_227, gk_228, \
                         gk_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_0 * gk_222[k];

        t_223[k] = f_0 * gk_223[k];

        t_224[k] = f_0 * gk_224[k];

        t_225[k] = f_0 * gk_225[k];

        t_226[k] = f_0 * gk_226[k];

        t_227[k] = f_0 * gk_227[k];

        t_228[k] = f_0 * gk_228[k];

        t_229[k] = f_0 * gk_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, t_235, t_236, t_237, gk_230, \
                         gk_231, gk_232, gk_233, gk_234, gk_235, gk_236, \
                         gk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_0 * gk_230[k];

        t_231[k] = f_0 * gk_231[k];

        t_232[k] = f_0 * gk_232[k];

        t_233[k] = f_0 * gk_233[k];

        t_234[k] = f_0 * gk_234[k];

        t_235[k] = f_0 * gk_235[k];

        t_236[k] = f_0 * gk_236[k];

        t_237[k] = f_0 * gk_237[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, t_243, t_244, t_245, gk_238, \
                         gk_239, gk_240, gk_241, gk_242, gk_243, gk_244, \
                         gk_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_0 * gk_238[k];

        t_239[k] = f_0 * gk_239[k];

        t_240[k] = f_0 * gk_240[k];

        t_241[k] = f_0 * gk_241[k];

        t_242[k] = f_0 * gk_242[k];

        t_243[k] = f_0 * gk_243[k];

        t_244[k] = f_0 * gk_244[k];

        t_245[k] = f_0 * gk_245[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, t_251, t_252, t_253, gk_246, \
                         gk_247, gk_248, gk_249, gk_250, gk_251, gk_252, \
                         gk_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_0 * gk_246[k];

        t_247[k] = f_0 * gk_247[k];

        t_248[k] = f_0 * gk_248[k];

        t_249[k] = f_0 * gk_249[k];

        t_250[k] = f_0 * gk_250[k];

        t_251[k] = f_0 * gk_251[k];

        t_252[k] = f_0 * gk_252[k];

        t_253[k] = f_0 * gk_253[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, t_258, t_259, t_260, t_261, gk_254, \
                         gk_255, gk_256, gk_257, gk_258, gk_259, gk_260, \
                         gk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_0 * gk_254[k];

        t_255[k] = f_0 * gk_255[k];

        t_256[k] = f_0 * gk_256[k];

        t_257[k] = f_0 * gk_257[k];

        t_258[k] = f_0 * gk_258[k];

        t_259[k] = f_0 * gk_259[k];

        t_260[k] = f_0 * gk_260[k];

        t_261[k] = f_0 * gk_261[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, t_267, t_268, t_269, gk_262, \
                         gk_263, gk_264, gk_265, gk_266, gk_267, gk_268, \
                         gk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_0 * gk_262[k];

        t_263[k] = f_0 * gk_263[k];

        t_264[k] = f_0 * gk_264[k];

        t_265[k] = f_0 * gk_265[k];

        t_266[k] = f_0 * gk_266[k];

        t_267[k] = f_0 * gk_267[k];

        t_268[k] = f_0 * gk_268[k];

        t_269[k] = f_0 * gk_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, t_276, t_277, gk_270, \
                         gk_271, gk_272, gk_273, gk_274, gk_275, gk_276, \
                         gk_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_0 * gk_270[k];

        t_271[k] = f_0 * gk_271[k];

        t_272[k] = f_0 * gk_272[k];

        t_273[k] = f_0 * gk_273[k];

        t_274[k] = f_0 * gk_274[k];

        t_275[k] = f_0 * gk_275[k];

        t_276[k] = f_0 * gk_276[k];

        t_277[k] = f_0 * gk_277[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, t_283, t_284, t_285, gk_278, \
                         gk_279, gk_280, gk_281, gk_282, gk_283, gk_284, \
                         gk_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_0 * gk_278[k];

        t_279[k] = f_0 * gk_279[k];

        t_280[k] = f_0 * gk_280[k];

        t_281[k] = f_0 * gk_281[k];

        t_282[k] = f_0 * gk_282[k];

        t_283[k] = f_0 * gk_283[k];

        t_284[k] = f_0 * gk_284[k];

        t_285[k] = f_0 * gk_285[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, t_291, t_292, t_293, gk_286, \
                         gk_287, gk_288, gk_289, gk_290, gk_291, gk_292, \
                         gk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * gk_286[k];

        t_287[k] = f_0 * gk_287[k];

        t_288[k] = f_0 * gk_288[k];

        t_289[k] = f_0 * gk_289[k];

        t_290[k] = f_0 * gk_290[k];

        t_291[k] = f_0 * gk_291[k];

        t_292[k] = f_0 * gk_292[k];

        t_293[k] = f_0 * gk_293[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, t_299, t_300, t_301, gk_294, \
                         gk_295, gk_296, gk_297, gk_298, gk_299, gk_300, \
                         gk_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_0 * gk_294[k];

        t_295[k] = f_0 * gk_295[k];

        t_296[k] = f_0 * gk_296[k];

        t_297[k] = f_0 * gk_297[k];

        t_298[k] = f_0 * gk_298[k];

        t_299[k] = f_0 * gk_299[k];

        t_300[k] = f_0 * gk_300[k];

        t_301[k] = f_0 * gk_301[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, t_307, t_308, t_309, gk_302, \
                         gk_303, gk_304, gk_305, gk_306, gk_307, gk_308, \
                         gk_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_0 * gk_302[k];

        t_303[k] = f_0 * gk_303[k];

        t_304[k] = f_0 * gk_304[k];

        t_305[k] = f_0 * gk_305[k];

        t_306[k] = f_0 * gk_306[k];

        t_307[k] = f_0 * gk_307[k];

        t_308[k] = f_0 * gk_308[k];

        t_309[k] = f_0 * gk_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, t_315, t_316, t_317, gk_310, \
                         gk_311, gk_312, gk_313, gk_314, gk_315, gk_316, \
                         gk_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_0 * gk_310[k];

        t_311[k] = f_0 * gk_311[k];

        t_312[k] = f_0 * gk_312[k];

        t_313[k] = f_0 * gk_313[k];

        t_314[k] = f_0 * gk_314[k];

        t_315[k] = f_0 * gk_315[k];

        t_316[k] = f_0 * gk_316[k];

        t_317[k] = f_0 * gk_317[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, t_323, t_324, t_325, gk_318, \
                         gk_319, gk_320, gk_321, gk_322, gk_323, gk_324, \
                         gk_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_0 * gk_318[k];

        t_319[k] = f_0 * gk_319[k];

        t_320[k] = f_0 * gk_320[k];

        t_321[k] = f_0 * gk_321[k];

        t_322[k] = f_0 * gk_322[k];

        t_323[k] = f_0 * gk_323[k];

        t_324[k] = f_0 * gk_324[k];

        t_325[k] = f_0 * gk_325[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, t_332, t_333, gk_326, \
                         gk_327, gk_328, gk_329, gk_330, gk_331, gk_332, \
                         gk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_0 * gk_326[k];

        t_327[k] = f_0 * gk_327[k];

        t_328[k] = f_0 * gk_328[k];

        t_329[k] = f_0 * gk_329[k];

        t_330[k] = f_0 * gk_330[k];

        t_331[k] = f_0 * gk_331[k];

        t_332[k] = f_0 * gk_332[k];

        t_333[k] = f_0 * gk_333[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, t_339, t_340, t_341, gk_334, \
                         gk_335, gk_336, gk_337, gk_338, gk_339, gk_340, \
                         gk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_0 * gk_334[k];

        t_335[k] = f_0 * gk_335[k];

        t_336[k] = f_0 * gk_336[k];

        t_337[k] = f_0 * gk_337[k];

        t_338[k] = f_0 * gk_338[k];

        t_339[k] = f_0 * gk_339[k];

        t_340[k] = f_0 * gk_340[k];

        t_341[k] = f_0 * gk_341[k];
    }
}

static auto
compute_prim_geom_10_fk_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gk_342 = buffer.data(gk + 342);
    const auto *gk_343 = buffer.data(gk + 343);
    const auto *gk_344 = buffer.data(gk + 344);
    const auto *gk_345 = buffer.data(gk + 345);
    const auto *gk_346 = buffer.data(gk + 346);
    const auto *gk_347 = buffer.data(gk + 347);
    const auto *gk_348 = buffer.data(gk + 348);
    const auto *gk_349 = buffer.data(gk + 349);
    const auto *gk_350 = buffer.data(gk + 350);
    const auto *gk_351 = buffer.data(gk + 351);
    const auto *gk_352 = buffer.data(gk + 352);
    const auto *gk_353 = buffer.data(gk + 353);
    const auto *gk_354 = buffer.data(gk + 354);
    const auto *gk_355 = buffer.data(gk + 355);
    const auto *gk_356 = buffer.data(gk + 356);
    const auto *gk_357 = buffer.data(gk + 357);
    const auto *gk_358 = buffer.data(gk + 358);
    const auto *gk_359 = buffer.data(gk + 359);

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, t_347, t_348, t_349, gk_342, \
                         gk_343, gk_344, gk_345, gk_346, gk_347, gk_348, \
                         gk_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_0 * gk_342[k];

        t_343[k] = f_0 * gk_343[k];

        t_344[k] = f_0 * gk_344[k];

        t_345[k] = f_0 * gk_345[k];

        t_346[k] = f_0 * gk_346[k];

        t_347[k] = f_0 * gk_347[k];

        t_348[k] = f_0 * gk_348[k];

        t_349[k] = f_0 * gk_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, t_355, t_356, t_357, gk_350, \
                         gk_351, gk_352, gk_353, gk_354, gk_355, gk_356, \
                         gk_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_0 * gk_350[k];

        t_351[k] = f_0 * gk_351[k];

        t_352[k] = f_0 * gk_352[k];

        t_353[k] = f_0 * gk_353[k];

        t_354[k] = f_0 * gk_354[k];

        t_355[k] = f_0 * gk_355[k];

        t_356[k] = f_0 * gk_356[k];

        t_357[k] = f_0 * gk_357[k];
    }

#pragma omp simd aligned(t_358, t_359, gk_358, gk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_0 * gk_358[k];

        t_359[k] = f_0 * gk_359[k];
    }
}

auto
compute_prim_geom_10_fk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t dk, const size_t gk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fk_electron_repulsion_0_piece0(buffer, target, dk, gk, ncols, alpha);

    compute_prim_geom_10_fk_electron_repulsion_0_piece1(buffer, target, dk, gk, ncols, alpha);

    compute_prim_geom_10_fk_electron_repulsion_0_piece2(buffer, target, gk, ncols, alpha);
}

static auto
compute_prim_geom_10_fk_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dk, const size_t gk,
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

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_6 = buffer.data(dk + 6);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_9 = buffer.data(dk + 9);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_13 = buffer.data(dk + 13);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_22 = buffer.data(dk + 22);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_25 = buffer.data(dk + 25);
    const auto *dk_26 = buffer.data(dk + 26);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_29 = buffer.data(dk + 29);
    const auto *dk_30 = buffer.data(dk + 30);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_33 = buffer.data(dk + 33);
    const auto *dk_34 = buffer.data(dk + 34);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_36 = buffer.data(dk + 36);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_49 = buffer.data(dk + 49);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_53 = buffer.data(dk + 53);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_60 = buffer.data(dk + 60);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_62 = buffer.data(dk + 62);
    const auto *dk_63 = buffer.data(dk + 63);
    const auto *dk_64 = buffer.data(dk + 64);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_71 = buffer.data(dk + 71);
    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_76 = buffer.data(dk + 76);
    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_80 = buffer.data(dk + 80);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_83 = buffer.data(dk + 83);
    const auto *dk_84 = buffer.data(dk + 84);
    const auto *dk_85 = buffer.data(dk + 85);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_87 = buffer.data(dk + 87);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_89 = buffer.data(dk + 89);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_94 = buffer.data(dk + 94);
    const auto *dk_95 = buffer.data(dk + 95);
    const auto *dk_96 = buffer.data(dk + 96);
    const auto *dk_97 = buffer.data(dk + 97);
    const auto *dk_98 = buffer.data(dk + 98);
    const auto *dk_99 = buffer.data(dk + 99);

    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_41 = buffer.data(gk + 41);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_43 = buffer.data(gk + 43);
    const auto *gk_44 = buffer.data(gk + 44);
    const auto *gk_45 = buffer.data(gk + 45);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_47 = buffer.data(gk + 47);
    const auto *gk_48 = buffer.data(gk + 48);
    const auto *gk_49 = buffer.data(gk + 49);
    const auto *gk_50 = buffer.data(gk + 50);
    const auto *gk_51 = buffer.data(gk + 51);
    const auto *gk_52 = buffer.data(gk + 52);
    const auto *gk_53 = buffer.data(gk + 53);
    const auto *gk_54 = buffer.data(gk + 54);
    const auto *gk_55 = buffer.data(gk + 55);
    const auto *gk_56 = buffer.data(gk + 56);
    const auto *gk_57 = buffer.data(gk + 57);
    const auto *gk_58 = buffer.data(gk + 58);
    const auto *gk_59 = buffer.data(gk + 59);
    const auto *gk_60 = buffer.data(gk + 60);
    const auto *gk_61 = buffer.data(gk + 61);
    const auto *gk_62 = buffer.data(gk + 62);
    const auto *gk_63 = buffer.data(gk + 63);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_65 = buffer.data(gk + 65);
    const auto *gk_66 = buffer.data(gk + 66);
    const auto *gk_67 = buffer.data(gk + 67);
    const auto *gk_68 = buffer.data(gk + 68);
    const auto *gk_69 = buffer.data(gk + 69);
    const auto *gk_70 = buffer.data(gk + 70);
    const auto *gk_71 = buffer.data(gk + 71);
    const auto *gk_108 = buffer.data(gk + 108);
    const auto *gk_109 = buffer.data(gk + 109);
    const auto *gk_110 = buffer.data(gk + 110);
    const auto *gk_111 = buffer.data(gk + 111);
    const auto *gk_112 = buffer.data(gk + 112);
    const auto *gk_113 = buffer.data(gk + 113);
    const auto *gk_114 = buffer.data(gk + 114);
    const auto *gk_115 = buffer.data(gk + 115);
    const auto *gk_116 = buffer.data(gk + 116);
    const auto *gk_117 = buffer.data(gk + 117);
    const auto *gk_118 = buffer.data(gk + 118);
    const auto *gk_119 = buffer.data(gk + 119);
    const auto *gk_120 = buffer.data(gk + 120);
    const auto *gk_121 = buffer.data(gk + 121);
    const auto *gk_122 = buffer.data(gk + 122);
    const auto *gk_123 = buffer.data(gk + 123);
    const auto *gk_124 = buffer.data(gk + 124);
    const auto *gk_125 = buffer.data(gk + 125);
    const auto *gk_126 = buffer.data(gk + 126);
    const auto *gk_127 = buffer.data(gk + 127);
    const auto *gk_128 = buffer.data(gk + 128);
    const auto *gk_129 = buffer.data(gk + 129);
    const auto *gk_130 = buffer.data(gk + 130);
    const auto *gk_131 = buffer.data(gk + 131);
    const auto *gk_132 = buffer.data(gk + 132);
    const auto *gk_133 = buffer.data(gk + 133);
    const auto *gk_134 = buffer.data(gk + 134);
    const auto *gk_135 = buffer.data(gk + 135);
    const auto *gk_136 = buffer.data(gk + 136);
    const auto *gk_137 = buffer.data(gk + 137);
    const auto *gk_138 = buffer.data(gk + 138);
    const auto *gk_139 = buffer.data(gk + 139);
    const auto *gk_140 = buffer.data(gk + 140);
    const auto *gk_141 = buffer.data(gk + 141);
    const auto *gk_142 = buffer.data(gk + 142);
    const auto *gk_143 = buffer.data(gk + 143);
    const auto *gk_144 = buffer.data(gk + 144);
    const auto *gk_145 = buffer.data(gk + 145);
    const auto *gk_146 = buffer.data(gk + 146);
    const auto *gk_147 = buffer.data(gk + 147);
    const auto *gk_148 = buffer.data(gk + 148);
    const auto *gk_149 = buffer.data(gk + 149);
    const auto *gk_150 = buffer.data(gk + 150);
    const auto *gk_151 = buffer.data(gk + 151);
    const auto *gk_152 = buffer.data(gk + 152);
    const auto *gk_153 = buffer.data(gk + 153);
    const auto *gk_154 = buffer.data(gk + 154);
    const auto *gk_155 = buffer.data(gk + 155);
    const auto *gk_156 = buffer.data(gk + 156);
    const auto *gk_157 = buffer.data(gk + 157);
    const auto *gk_158 = buffer.data(gk + 158);
    const auto *gk_159 = buffer.data(gk + 159);
    const auto *gk_160 = buffer.data(gk + 160);
    const auto *gk_161 = buffer.data(gk + 161);
    const auto *gk_162 = buffer.data(gk + 162);
    const auto *gk_163 = buffer.data(gk + 163);
    const auto *gk_164 = buffer.data(gk + 164);
    const auto *gk_165 = buffer.data(gk + 165);
    const auto *gk_166 = buffer.data(gk + 166);
    const auto *gk_167 = buffer.data(gk + 167);
    const auto *gk_168 = buffer.data(gk + 168);
    const auto *gk_169 = buffer.data(gk + 169);
    const auto *gk_170 = buffer.data(gk + 170);
    const auto *gk_171 = buffer.data(gk + 171);
    const auto *gk_172 = buffer.data(gk + 172);
    const auto *gk_173 = buffer.data(gk + 173);
    const auto *gk_174 = buffer.data(gk + 174);
    const auto *gk_175 = buffer.data(gk + 175);
    const auto *gk_176 = buffer.data(gk + 176);
    const auto *gk_177 = buffer.data(gk + 177);
    const auto *gk_178 = buffer.data(gk + 178);
    const auto *gk_179 = buffer.data(gk + 179);
    const auto *gk_216 = buffer.data(gk + 216);
    const auto *gk_217 = buffer.data(gk + 217);
    const auto *gk_218 = buffer.data(gk + 218);
    const auto *gk_219 = buffer.data(gk + 219);
    const auto *gk_220 = buffer.data(gk + 220);
    const auto *gk_221 = buffer.data(gk + 221);
    const auto *gk_222 = buffer.data(gk + 222);
    const auto *gk_223 = buffer.data(gk + 223);
    const auto *gk_224 = buffer.data(gk + 224);
    const auto *gk_225 = buffer.data(gk + 225);
    const auto *gk_226 = buffer.data(gk + 226);
    const auto *gk_227 = buffer.data(gk + 227);
    const auto *gk_228 = buffer.data(gk + 228);
    const auto *gk_229 = buffer.data(gk + 229);
    const auto *gk_230 = buffer.data(gk + 230);
    const auto *gk_231 = buffer.data(gk + 231);
    const auto *gk_232 = buffer.data(gk + 232);
    const auto *gk_233 = buffer.data(gk + 233);
    const auto *gk_234 = buffer.data(gk + 234);
    const auto *gk_235 = buffer.data(gk + 235);
    const auto *gk_236 = buffer.data(gk + 236);
    const auto *gk_237 = buffer.data(gk + 237);
    const auto *gk_238 = buffer.data(gk + 238);
    const auto *gk_239 = buffer.data(gk + 239);
    const auto *gk_240 = buffer.data(gk + 240);
    const auto *gk_241 = buffer.data(gk + 241);
    const auto *gk_242 = buffer.data(gk + 242);
    const auto *gk_243 = buffer.data(gk + 243);
    const auto *gk_244 = buffer.data(gk + 244);
    const auto *gk_245 = buffer.data(gk + 245);
    const auto *gk_246 = buffer.data(gk + 246);
    const auto *gk_247 = buffer.data(gk + 247);
    const auto *gk_248 = buffer.data(gk + 248);
    const auto *gk_249 = buffer.data(gk + 249);
    const auto *gk_250 = buffer.data(gk + 250);
    const auto *gk_251 = buffer.data(gk + 251);
    const auto *gk_252 = buffer.data(gk + 252);
    const auto *gk_253 = buffer.data(gk + 253);
    const auto *gk_254 = buffer.data(gk + 254);
    const auto *gk_255 = buffer.data(gk + 255);
    const auto *gk_256 = buffer.data(gk + 256);
    const auto *gk_257 = buffer.data(gk + 257);
    const auto *gk_258 = buffer.data(gk + 258);
    const auto *gk_259 = buffer.data(gk + 259);
    const auto *gk_260 = buffer.data(gk + 260);
    const auto *gk_261 = buffer.data(gk + 261);
    const auto *gk_262 = buffer.data(gk + 262);
    const auto *gk_263 = buffer.data(gk + 263);
    const auto *gk_264 = buffer.data(gk + 264);
    const auto *gk_265 = buffer.data(gk + 265);
    const auto *gk_266 = buffer.data(gk + 266);
    const auto *gk_267 = buffer.data(gk + 267);
    const auto *gk_268 = buffer.data(gk + 268);
    const auto *gk_269 = buffer.data(gk + 269);
    const auto *gk_270 = buffer.data(gk + 270);
    const auto *gk_271 = buffer.data(gk + 271);
    const auto *gk_272 = buffer.data(gk + 272);
    const auto *gk_273 = buffer.data(gk + 273);
    const auto *gk_274 = buffer.data(gk + 274);
    const auto *gk_275 = buffer.data(gk + 275);
    const auto *gk_276 = buffer.data(gk + 276);
    const auto *gk_277 = buffer.data(gk + 277);
    const auto *gk_278 = buffer.data(gk + 278);
    const auto *gk_279 = buffer.data(gk + 279);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gk_36, gk_37, gk_38, gk_39, \
                         gk_40, gk_41, gk_42, gk_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gk_36[k];

        t_1[k] = f_0 * gk_37[k];

        t_2[k] = f_0 * gk_38[k];

        t_3[k] = f_0 * gk_39[k];

        t_4[k] = f_0 * gk_40[k];

        t_5[k] = f_0 * gk_41[k];

        t_6[k] = f_0 * gk_42[k];

        t_7[k] = f_0 * gk_43[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, gk_44, gk_45, gk_46, \
                         gk_47, gk_48, gk_49, gk_50, gk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gk_44[k];

        t_9[k] = f_0 * gk_45[k];

        t_10[k] = f_0 * gk_46[k];

        t_11[k] = f_0 * gk_47[k];

        t_12[k] = f_0 * gk_48[k];

        t_13[k] = f_0 * gk_49[k];

        t_14[k] = f_0 * gk_50[k];

        t_15[k] = f_0 * gk_51[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, gk_52, gk_53, gk_54, \
                         gk_55, gk_56, gk_57, gk_58, gk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * gk_52[k];

        t_17[k] = f_0 * gk_53[k];

        t_18[k] = f_0 * gk_54[k];

        t_19[k] = f_0 * gk_55[k];

        t_20[k] = f_0 * gk_56[k];

        t_21[k] = f_0 * gk_57[k];

        t_22[k] = f_0 * gk_58[k];

        t_23[k] = f_0 * gk_59[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, gk_60, gk_61, gk_62, \
                         gk_63, gk_64, gk_65, gk_66, gk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * gk_60[k];

        t_25[k] = f_0 * gk_61[k];

        t_26[k] = f_0 * gk_62[k];

        t_27[k] = f_0 * gk_63[k];

        t_28[k] = f_0 * gk_64[k];

        t_29[k] = f_0 * gk_65[k];

        t_30[k] = f_0 * gk_66[k];

        t_31[k] = f_0 * gk_67[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, dk_0, dk_1, gk_68, gk_69, gk_70, \
                         gk_71, gk_108, gk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * gk_68[k];

        t_33[k] = f_0 * gk_69[k];

        t_34[k] = f_0 * gk_70[k];

        t_35[k] = f_0 * gk_71[k];

        t_36[k] = -dk_0[k]
                  + f_0 * gk_108[k];

        t_37[k] = -dk_1[k]
                  + f_0 * gk_109[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, dk_2, dk_3, dk_4, dk_5, dk_6, gk_110, \
                         gk_111, gk_112, gk_113, gk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -dk_2[k]
                  + f_0 * gk_110[k];

        t_39[k] = -dk_3[k]
                  + f_0 * gk_111[k];

        t_40[k] = -dk_4[k]
                  + f_0 * gk_112[k];

        t_41[k] = -dk_5[k]
                  + f_0 * gk_113[k];

        t_42[k] = -dk_6[k]
                  + f_0 * gk_114[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, dk_7, dk_8, dk_9, dk_10, dk_11, gk_115, \
                         gk_116, gk_117, gk_118, gk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -dk_7[k]
                  + f_0 * gk_115[k];

        t_44[k] = -dk_8[k]
                  + f_0 * gk_116[k];

        t_45[k] = -dk_9[k]
                  + f_0 * gk_117[k];

        t_46[k] = -dk_10[k]
                  + f_0 * gk_118[k];

        t_47[k] = -dk_11[k]
                  + f_0 * gk_119[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, dk_12, dk_13, dk_14, dk_15, dk_16, \
                         gk_120, gk_121, gk_122, gk_123, gk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -dk_12[k]
                  + f_0 * gk_120[k];

        t_49[k] = -dk_13[k]
                  + f_0 * gk_121[k];

        t_50[k] = -dk_14[k]
                  + f_0 * gk_122[k];

        t_51[k] = -dk_15[k]
                  + f_0 * gk_123[k];

        t_52[k] = -dk_16[k]
                  + f_0 * gk_124[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, dk_17, dk_18, dk_19, dk_20, dk_21, \
                         gk_125, gk_126, gk_127, gk_128, gk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -dk_17[k]
                  + f_0 * gk_125[k];

        t_54[k] = -dk_18[k]
                  + f_0 * gk_126[k];

        t_55[k] = -dk_19[k]
                  + f_0 * gk_127[k];

        t_56[k] = -dk_20[k]
                  + f_0 * gk_128[k];

        t_57[k] = -dk_21[k]
                  + f_0 * gk_129[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, dk_22, dk_23, dk_24, dk_25, dk_26, \
                         gk_130, gk_131, gk_132, gk_133, gk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -dk_22[k]
                  + f_0 * gk_130[k];

        t_59[k] = -dk_23[k]
                  + f_0 * gk_131[k];

        t_60[k] = -dk_24[k]
                  + f_0 * gk_132[k];

        t_61[k] = -dk_25[k]
                  + f_0 * gk_133[k];

        t_62[k] = -dk_26[k]
                  + f_0 * gk_134[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, dk_27, dk_28, dk_29, dk_30, dk_31, \
                         gk_135, gk_136, gk_137, gk_138, gk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -dk_27[k]
                  + f_0 * gk_135[k];

        t_64[k] = -dk_28[k]
                  + f_0 * gk_136[k];

        t_65[k] = -dk_29[k]
                  + f_0 * gk_137[k];

        t_66[k] = -dk_30[k]
                  + f_0 * gk_138[k];

        t_67[k] = -dk_31[k]
                  + f_0 * gk_139[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, dk_32, dk_33, dk_34, dk_35, \
                         gk_140, gk_141, gk_142, gk_143, gk_144, \
                         gk_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -dk_32[k]
                  + f_0 * gk_140[k];

        t_69[k] = -dk_33[k]
                  + f_0 * gk_141[k];

        t_70[k] = -dk_34[k]
                  + f_0 * gk_142[k];

        t_71[k] = -dk_35[k]
                  + f_0 * gk_143[k];

        t_72[k] = f_0 * gk_144[k];

        t_73[k] = f_0 * gk_145[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, t_80, t_81, gk_146, gk_147, \
                         gk_148, gk_149, gk_150, gk_151, gk_152, \
                         gk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * gk_146[k];

        t_75[k] = f_0 * gk_147[k];

        t_76[k] = f_0 * gk_148[k];

        t_77[k] = f_0 * gk_149[k];

        t_78[k] = f_0 * gk_150[k];

        t_79[k] = f_0 * gk_151[k];

        t_80[k] = f_0 * gk_152[k];

        t_81[k] = f_0 * gk_153[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, t_87, t_88, t_89, gk_154, gk_155, \
                         gk_156, gk_157, gk_158, gk_159, gk_160, \
                         gk_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * gk_154[k];

        t_83[k] = f_0 * gk_155[k];

        t_84[k] = f_0 * gk_156[k];

        t_85[k] = f_0 * gk_157[k];

        t_86[k] = f_0 * gk_158[k];

        t_87[k] = f_0 * gk_159[k];

        t_88[k] = f_0 * gk_160[k];

        t_89[k] = f_0 * gk_161[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, gk_162, gk_163, \
                         gk_164, gk_165, gk_166, gk_167, gk_168, \
                         gk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * gk_162[k];

        t_91[k] = f_0 * gk_163[k];

        t_92[k] = f_0 * gk_164[k];

        t_93[k] = f_0 * gk_165[k];

        t_94[k] = f_0 * gk_166[k];

        t_95[k] = f_0 * gk_167[k];

        t_96[k] = f_0 * gk_168[k];

        t_97[k] = f_0 * gk_169[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, t_104, t_105, gk_170, gk_171, \
                         gk_172, gk_173, gk_174, gk_175, gk_176, \
                         gk_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * gk_170[k];

        t_99[k] = f_0 * gk_171[k];

        t_100[k] = f_0 * gk_172[k];

        t_101[k] = f_0 * gk_173[k];

        t_102[k] = f_0 * gk_174[k];

        t_103[k] = f_0 * gk_175[k];

        t_104[k] = f_0 * gk_176[k];

        t_105[k] = f_0 * gk_177[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, dk_36, dk_37, dk_38, dk_39, \
                         gk_178, gk_179, gk_216, gk_217, gk_218, \
                         gk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * gk_178[k];

        t_107[k] = f_0 * gk_179[k];

        t_108[k] = -2.0 * dk_36[k]
                   + f_0 * gk_216[k];

        t_109[k] = -2.0 * dk_37[k]
                   + f_0 * gk_217[k];

        t_110[k] = -2.0 * dk_38[k]
                   + f_0 * gk_218[k];

        t_111[k] = -2.0 * dk_39[k]
                   + f_0 * gk_219[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, dk_40, dk_41, dk_42, dk_43, dk_44, \
                         gk_220, gk_221, gk_222, gk_223, gk_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -2.0 * dk_40[k]
                   + f_0 * gk_220[k];

        t_113[k] = -2.0 * dk_41[k]
                   + f_0 * gk_221[k];

        t_114[k] = -2.0 * dk_42[k]
                   + f_0 * gk_222[k];

        t_115[k] = -2.0 * dk_43[k]
                   + f_0 * gk_223[k];

        t_116[k] = -2.0 * dk_44[k]
                   + f_0 * gk_224[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, dk_45, dk_46, dk_47, dk_48, dk_49, \
                         gk_225, gk_226, gk_227, gk_228, gk_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -2.0 * dk_45[k]
                   + f_0 * gk_225[k];

        t_118[k] = -2.0 * dk_46[k]
                   + f_0 * gk_226[k];

        t_119[k] = -2.0 * dk_47[k]
                   + f_0 * gk_227[k];

        t_120[k] = -2.0 * dk_48[k]
                   + f_0 * gk_228[k];

        t_121[k] = -2.0 * dk_49[k]
                   + f_0 * gk_229[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, dk_50, dk_51, dk_52, dk_53, dk_54, \
                         gk_230, gk_231, gk_232, gk_233, gk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * dk_50[k]
                   + f_0 * gk_230[k];

        t_123[k] = -2.0 * dk_51[k]
                   + f_0 * gk_231[k];

        t_124[k] = -2.0 * dk_52[k]
                   + f_0 * gk_232[k];

        t_125[k] = -2.0 * dk_53[k]
                   + f_0 * gk_233[k];

        t_126[k] = -2.0 * dk_54[k]
                   + f_0 * gk_234[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, dk_55, dk_56, dk_57, dk_58, dk_59, \
                         gk_235, gk_236, gk_237, gk_238, gk_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * dk_55[k]
                   + f_0 * gk_235[k];

        t_128[k] = -2.0 * dk_56[k]
                   + f_0 * gk_236[k];

        t_129[k] = -2.0 * dk_57[k]
                   + f_0 * gk_237[k];

        t_130[k] = -2.0 * dk_58[k]
                   + f_0 * gk_238[k];

        t_131[k] = -2.0 * dk_59[k]
                   + f_0 * gk_239[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, dk_60, dk_61, dk_62, dk_63, dk_64, \
                         gk_240, gk_241, gk_242, gk_243, gk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * dk_60[k]
                   + f_0 * gk_240[k];

        t_133[k] = -2.0 * dk_61[k]
                   + f_0 * gk_241[k];

        t_134[k] = -2.0 * dk_62[k]
                   + f_0 * gk_242[k];

        t_135[k] = -2.0 * dk_63[k]
                   + f_0 * gk_243[k];

        t_136[k] = -2.0 * dk_64[k]
                   + f_0 * gk_244[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, dk_65, dk_66, dk_67, dk_68, dk_69, \
                         gk_245, gk_246, gk_247, gk_248, gk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -2.0 * dk_65[k]
                   + f_0 * gk_245[k];

        t_138[k] = -2.0 * dk_66[k]
                   + f_0 * gk_246[k];

        t_139[k] = -2.0 * dk_67[k]
                   + f_0 * gk_247[k];

        t_140[k] = -2.0 * dk_68[k]
                   + f_0 * gk_248[k];

        t_141[k] = -2.0 * dk_69[k]
                   + f_0 * gk_249[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, dk_70, dk_71, dk_72, dk_73, dk_74, \
                         gk_250, gk_251, gk_252, gk_253, gk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -2.0 * dk_70[k]
                   + f_0 * gk_250[k];

        t_143[k] = -2.0 * dk_71[k]
                   + f_0 * gk_251[k];

        t_144[k] = -dk_72[k]
                   + f_0 * gk_252[k];

        t_145[k] = -dk_73[k]
                   + f_0 * gk_253[k];

        t_146[k] = -dk_74[k]
                   + f_0 * gk_254[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, dk_75, dk_76, dk_77, dk_78, dk_79, \
                         gk_255, gk_256, gk_257, gk_258, gk_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -dk_75[k]
                   + f_0 * gk_255[k];

        t_148[k] = -dk_76[k]
                   + f_0 * gk_256[k];

        t_149[k] = -dk_77[k]
                   + f_0 * gk_257[k];

        t_150[k] = -dk_78[k]
                   + f_0 * gk_258[k];

        t_151[k] = -dk_79[k]
                   + f_0 * gk_259[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, dk_80, dk_81, dk_82, dk_83, dk_84, \
                         gk_260, gk_261, gk_262, gk_263, gk_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -dk_80[k]
                   + f_0 * gk_260[k];

        t_153[k] = -dk_81[k]
                   + f_0 * gk_261[k];

        t_154[k] = -dk_82[k]
                   + f_0 * gk_262[k];

        t_155[k] = -dk_83[k]
                   + f_0 * gk_263[k];

        t_156[k] = -dk_84[k]
                   + f_0 * gk_264[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, dk_85, dk_86, dk_87, dk_88, dk_89, \
                         gk_265, gk_266, gk_267, gk_268, gk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -dk_85[k]
                   + f_0 * gk_265[k];

        t_158[k] = -dk_86[k]
                   + f_0 * gk_266[k];

        t_159[k] = -dk_87[k]
                   + f_0 * gk_267[k];

        t_160[k] = -dk_88[k]
                   + f_0 * gk_268[k];

        t_161[k] = -dk_89[k]
                   + f_0 * gk_269[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, dk_90, dk_91, dk_92, dk_93, dk_94, \
                         gk_270, gk_271, gk_272, gk_273, gk_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -dk_90[k]
                   + f_0 * gk_270[k];

        t_163[k] = -dk_91[k]
                   + f_0 * gk_271[k];

        t_164[k] = -dk_92[k]
                   + f_0 * gk_272[k];

        t_165[k] = -dk_93[k]
                   + f_0 * gk_273[k];

        t_166[k] = -dk_94[k]
                   + f_0 * gk_274[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, dk_95, dk_96, dk_97, dk_98, dk_99, \
                         gk_275, gk_276, gk_277, gk_278, gk_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -dk_95[k]
                   + f_0 * gk_275[k];

        t_168[k] = -dk_96[k]
                   + f_0 * gk_276[k];

        t_169[k] = -dk_97[k]
                   + f_0 * gk_277[k];

        t_170[k] = -dk_98[k]
                   + f_0 * gk_278[k];

        t_171[k] = -dk_99[k]
                   + f_0 * gk_279[k];
    }
}

static auto
compute_prim_geom_10_fk_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dk, const size_t gk,
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
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);

    const auto *dk_100 = buffer.data(dk + 100);
    const auto *dk_101 = buffer.data(dk + 101);
    const auto *dk_102 = buffer.data(dk + 102);
    const auto *dk_103 = buffer.data(dk + 103);
    const auto *dk_104 = buffer.data(dk + 104);
    const auto *dk_105 = buffer.data(dk + 105);
    const auto *dk_106 = buffer.data(dk + 106);
    const auto *dk_107 = buffer.data(dk + 107);
    const auto *dk_108 = buffer.data(dk + 108);
    const auto *dk_109 = buffer.data(dk + 109);
    const auto *dk_110 = buffer.data(dk + 110);
    const auto *dk_111 = buffer.data(dk + 111);
    const auto *dk_112 = buffer.data(dk + 112);
    const auto *dk_113 = buffer.data(dk + 113);
    const auto *dk_114 = buffer.data(dk + 114);
    const auto *dk_115 = buffer.data(dk + 115);
    const auto *dk_116 = buffer.data(dk + 116);
    const auto *dk_117 = buffer.data(dk + 117);
    const auto *dk_118 = buffer.data(dk + 118);
    const auto *dk_119 = buffer.data(dk + 119);
    const auto *dk_120 = buffer.data(dk + 120);
    const auto *dk_121 = buffer.data(dk + 121);
    const auto *dk_122 = buffer.data(dk + 122);
    const auto *dk_123 = buffer.data(dk + 123);
    const auto *dk_124 = buffer.data(dk + 124);
    const auto *dk_125 = buffer.data(dk + 125);
    const auto *dk_126 = buffer.data(dk + 126);
    const auto *dk_127 = buffer.data(dk + 127);
    const auto *dk_128 = buffer.data(dk + 128);
    const auto *dk_129 = buffer.data(dk + 129);
    const auto *dk_130 = buffer.data(dk + 130);
    const auto *dk_131 = buffer.data(dk + 131);
    const auto *dk_132 = buffer.data(dk + 132);
    const auto *dk_133 = buffer.data(dk + 133);
    const auto *dk_134 = buffer.data(dk + 134);
    const auto *dk_135 = buffer.data(dk + 135);
    const auto *dk_136 = buffer.data(dk + 136);
    const auto *dk_137 = buffer.data(dk + 137);
    const auto *dk_138 = buffer.data(dk + 138);
    const auto *dk_139 = buffer.data(dk + 139);
    const auto *dk_140 = buffer.data(dk + 140);
    const auto *dk_141 = buffer.data(dk + 141);
    const auto *dk_142 = buffer.data(dk + 142);
    const auto *dk_143 = buffer.data(dk + 143);
    const auto *dk_144 = buffer.data(dk + 144);
    const auto *dk_145 = buffer.data(dk + 145);
    const auto *dk_146 = buffer.data(dk + 146);
    const auto *dk_147 = buffer.data(dk + 147);
    const auto *dk_148 = buffer.data(dk + 148);
    const auto *dk_149 = buffer.data(dk + 149);
    const auto *dk_150 = buffer.data(dk + 150);
    const auto *dk_151 = buffer.data(dk + 151);
    const auto *dk_152 = buffer.data(dk + 152);
    const auto *dk_153 = buffer.data(dk + 153);
    const auto *dk_154 = buffer.data(dk + 154);
    const auto *dk_155 = buffer.data(dk + 155);
    const auto *dk_156 = buffer.data(dk + 156);
    const auto *dk_157 = buffer.data(dk + 157);
    const auto *dk_158 = buffer.data(dk + 158);
    const auto *dk_159 = buffer.data(dk + 159);
    const auto *dk_160 = buffer.data(dk + 160);
    const auto *dk_161 = buffer.data(dk + 161);
    const auto *dk_162 = buffer.data(dk + 162);
    const auto *dk_163 = buffer.data(dk + 163);
    const auto *dk_164 = buffer.data(dk + 164);
    const auto *dk_165 = buffer.data(dk + 165);
    const auto *dk_166 = buffer.data(dk + 166);
    const auto *dk_167 = buffer.data(dk + 167);
    const auto *dk_168 = buffer.data(dk + 168);
    const auto *dk_169 = buffer.data(dk + 169);
    const auto *dk_170 = buffer.data(dk + 170);
    const auto *dk_171 = buffer.data(dk + 171);
    const auto *dk_172 = buffer.data(dk + 172);
    const auto *dk_173 = buffer.data(dk + 173);
    const auto *dk_174 = buffer.data(dk + 174);
    const auto *dk_175 = buffer.data(dk + 175);
    const auto *dk_176 = buffer.data(dk + 176);
    const auto *dk_177 = buffer.data(dk + 177);
    const auto *dk_178 = buffer.data(dk + 178);
    const auto *dk_179 = buffer.data(dk + 179);
    const auto *dk_180 = buffer.data(dk + 180);
    const auto *dk_181 = buffer.data(dk + 181);
    const auto *dk_182 = buffer.data(dk + 182);
    const auto *dk_183 = buffer.data(dk + 183);
    const auto *dk_184 = buffer.data(dk + 184);
    const auto *dk_185 = buffer.data(dk + 185);
    const auto *dk_186 = buffer.data(dk + 186);
    const auto *dk_187 = buffer.data(dk + 187);
    const auto *dk_188 = buffer.data(dk + 188);
    const auto *dk_189 = buffer.data(dk + 189);
    const auto *dk_190 = buffer.data(dk + 190);
    const auto *dk_191 = buffer.data(dk + 191);
    const auto *dk_192 = buffer.data(dk + 192);
    const auto *dk_193 = buffer.data(dk + 193);
    const auto *dk_194 = buffer.data(dk + 194);
    const auto *dk_195 = buffer.data(dk + 195);
    const auto *dk_196 = buffer.data(dk + 196);
    const auto *dk_197 = buffer.data(dk + 197);
    const auto *dk_198 = buffer.data(dk + 198);
    const auto *dk_199 = buffer.data(dk + 199);
    const auto *dk_200 = buffer.data(dk + 200);
    const auto *dk_201 = buffer.data(dk + 201);
    const auto *dk_202 = buffer.data(dk + 202);
    const auto *dk_203 = buffer.data(dk + 203);
    const auto *dk_204 = buffer.data(dk + 204);
    const auto *dk_205 = buffer.data(dk + 205);
    const auto *dk_206 = buffer.data(dk + 206);
    const auto *dk_207 = buffer.data(dk + 207);
    const auto *dk_208 = buffer.data(dk + 208);
    const auto *dk_209 = buffer.data(dk + 209);
    const auto *dk_210 = buffer.data(dk + 210);
    const auto *dk_211 = buffer.data(dk + 211);
    const auto *dk_212 = buffer.data(dk + 212);
    const auto *dk_213 = buffer.data(dk + 213);
    const auto *dk_214 = buffer.data(dk + 214);
    const auto *dk_215 = buffer.data(dk + 215);

    const auto *gk_280 = buffer.data(gk + 280);
    const auto *gk_281 = buffer.data(gk + 281);
    const auto *gk_282 = buffer.data(gk + 282);
    const auto *gk_283 = buffer.data(gk + 283);
    const auto *gk_284 = buffer.data(gk + 284);
    const auto *gk_285 = buffer.data(gk + 285);
    const auto *gk_286 = buffer.data(gk + 286);
    const auto *gk_287 = buffer.data(gk + 287);
    const auto *gk_288 = buffer.data(gk + 288);
    const auto *gk_289 = buffer.data(gk + 289);
    const auto *gk_290 = buffer.data(gk + 290);
    const auto *gk_291 = buffer.data(gk + 291);
    const auto *gk_292 = buffer.data(gk + 292);
    const auto *gk_293 = buffer.data(gk + 293);
    const auto *gk_294 = buffer.data(gk + 294);
    const auto *gk_295 = buffer.data(gk + 295);
    const auto *gk_296 = buffer.data(gk + 296);
    const auto *gk_297 = buffer.data(gk + 297);
    const auto *gk_298 = buffer.data(gk + 298);
    const auto *gk_299 = buffer.data(gk + 299);
    const auto *gk_300 = buffer.data(gk + 300);
    const auto *gk_301 = buffer.data(gk + 301);
    const auto *gk_302 = buffer.data(gk + 302);
    const auto *gk_303 = buffer.data(gk + 303);
    const auto *gk_304 = buffer.data(gk + 304);
    const auto *gk_305 = buffer.data(gk + 305);
    const auto *gk_306 = buffer.data(gk + 306);
    const auto *gk_307 = buffer.data(gk + 307);
    const auto *gk_308 = buffer.data(gk + 308);
    const auto *gk_309 = buffer.data(gk + 309);
    const auto *gk_310 = buffer.data(gk + 310);
    const auto *gk_311 = buffer.data(gk + 311);
    const auto *gk_312 = buffer.data(gk + 312);
    const auto *gk_313 = buffer.data(gk + 313);
    const auto *gk_314 = buffer.data(gk + 314);
    const auto *gk_315 = buffer.data(gk + 315);
    const auto *gk_316 = buffer.data(gk + 316);
    const auto *gk_317 = buffer.data(gk + 317);
    const auto *gk_318 = buffer.data(gk + 318);
    const auto *gk_319 = buffer.data(gk + 319);
    const auto *gk_320 = buffer.data(gk + 320);
    const auto *gk_321 = buffer.data(gk + 321);
    const auto *gk_322 = buffer.data(gk + 322);
    const auto *gk_323 = buffer.data(gk + 323);
    const auto *gk_360 = buffer.data(gk + 360);
    const auto *gk_361 = buffer.data(gk + 361);
    const auto *gk_362 = buffer.data(gk + 362);
    const auto *gk_363 = buffer.data(gk + 363);
    const auto *gk_364 = buffer.data(gk + 364);
    const auto *gk_365 = buffer.data(gk + 365);
    const auto *gk_366 = buffer.data(gk + 366);
    const auto *gk_367 = buffer.data(gk + 367);
    const auto *gk_368 = buffer.data(gk + 368);
    const auto *gk_369 = buffer.data(gk + 369);
    const auto *gk_370 = buffer.data(gk + 370);
    const auto *gk_371 = buffer.data(gk + 371);
    const auto *gk_372 = buffer.data(gk + 372);
    const auto *gk_373 = buffer.data(gk + 373);
    const auto *gk_374 = buffer.data(gk + 374);
    const auto *gk_375 = buffer.data(gk + 375);
    const auto *gk_376 = buffer.data(gk + 376);
    const auto *gk_377 = buffer.data(gk + 377);
    const auto *gk_378 = buffer.data(gk + 378);
    const auto *gk_379 = buffer.data(gk + 379);
    const auto *gk_380 = buffer.data(gk + 380);
    const auto *gk_381 = buffer.data(gk + 381);
    const auto *gk_382 = buffer.data(gk + 382);
    const auto *gk_383 = buffer.data(gk + 383);
    const auto *gk_384 = buffer.data(gk + 384);
    const auto *gk_385 = buffer.data(gk + 385);
    const auto *gk_386 = buffer.data(gk + 386);
    const auto *gk_387 = buffer.data(gk + 387);
    const auto *gk_388 = buffer.data(gk + 388);
    const auto *gk_389 = buffer.data(gk + 389);
    const auto *gk_390 = buffer.data(gk + 390);
    const auto *gk_391 = buffer.data(gk + 391);
    const auto *gk_392 = buffer.data(gk + 392);
    const auto *gk_393 = buffer.data(gk + 393);
    const auto *gk_394 = buffer.data(gk + 394);
    const auto *gk_395 = buffer.data(gk + 395);
    const auto *gk_396 = buffer.data(gk + 396);
    const auto *gk_397 = buffer.data(gk + 397);
    const auto *gk_398 = buffer.data(gk + 398);
    const auto *gk_399 = buffer.data(gk + 399);
    const auto *gk_400 = buffer.data(gk + 400);
    const auto *gk_401 = buffer.data(gk + 401);
    const auto *gk_402 = buffer.data(gk + 402);
    const auto *gk_403 = buffer.data(gk + 403);
    const auto *gk_404 = buffer.data(gk + 404);
    const auto *gk_405 = buffer.data(gk + 405);
    const auto *gk_406 = buffer.data(gk + 406);
    const auto *gk_407 = buffer.data(gk + 407);
    const auto *gk_408 = buffer.data(gk + 408);
    const auto *gk_409 = buffer.data(gk + 409);
    const auto *gk_410 = buffer.data(gk + 410);
    const auto *gk_411 = buffer.data(gk + 411);
    const auto *gk_412 = buffer.data(gk + 412);
    const auto *gk_413 = buffer.data(gk + 413);
    const auto *gk_414 = buffer.data(gk + 414);
    const auto *gk_415 = buffer.data(gk + 415);
    const auto *gk_416 = buffer.data(gk + 416);
    const auto *gk_417 = buffer.data(gk + 417);
    const auto *gk_418 = buffer.data(gk + 418);
    const auto *gk_419 = buffer.data(gk + 419);
    const auto *gk_420 = buffer.data(gk + 420);
    const auto *gk_421 = buffer.data(gk + 421);
    const auto *gk_422 = buffer.data(gk + 422);
    const auto *gk_423 = buffer.data(gk + 423);
    const auto *gk_424 = buffer.data(gk + 424);
    const auto *gk_425 = buffer.data(gk + 425);
    const auto *gk_426 = buffer.data(gk + 426);
    const auto *gk_427 = buffer.data(gk + 427);
    const auto *gk_428 = buffer.data(gk + 428);
    const auto *gk_429 = buffer.data(gk + 429);
    const auto *gk_430 = buffer.data(gk + 430);
    const auto *gk_431 = buffer.data(gk + 431);
    const auto *gk_432 = buffer.data(gk + 432);
    const auto *gk_433 = buffer.data(gk + 433);
    const auto *gk_434 = buffer.data(gk + 434);
    const auto *gk_435 = buffer.data(gk + 435);
    const auto *gk_436 = buffer.data(gk + 436);
    const auto *gk_437 = buffer.data(gk + 437);
    const auto *gk_438 = buffer.data(gk + 438);
    const auto *gk_439 = buffer.data(gk + 439);
    const auto *gk_440 = buffer.data(gk + 440);
    const auto *gk_441 = buffer.data(gk + 441);
    const auto *gk_442 = buffer.data(gk + 442);
    const auto *gk_443 = buffer.data(gk + 443);
    const auto *gk_444 = buffer.data(gk + 444);
    const auto *gk_445 = buffer.data(gk + 445);
    const auto *gk_446 = buffer.data(gk + 446);
    const auto *gk_447 = buffer.data(gk + 447);
    const auto *gk_448 = buffer.data(gk + 448);
    const auto *gk_449 = buffer.data(gk + 449);
    const auto *gk_450 = buffer.data(gk + 450);
    const auto *gk_451 = buffer.data(gk + 451);
    const auto *gk_452 = buffer.data(gk + 452);
    const auto *gk_453 = buffer.data(gk + 453);
    const auto *gk_454 = buffer.data(gk + 454);
    const auto *gk_455 = buffer.data(gk + 455);
    const auto *gk_456 = buffer.data(gk + 456);
    const auto *gk_457 = buffer.data(gk + 457);
    const auto *gk_458 = buffer.data(gk + 458);
    const auto *gk_459 = buffer.data(gk + 459);
    const auto *gk_460 = buffer.data(gk + 460);
    const auto *gk_461 = buffer.data(gk + 461);
    const auto *gk_462 = buffer.data(gk + 462);
    const auto *gk_463 = buffer.data(gk + 463);
    const auto *gk_464 = buffer.data(gk + 464);
    const auto *gk_465 = buffer.data(gk + 465);
    const auto *gk_466 = buffer.data(gk + 466);
    const auto *gk_467 = buffer.data(gk + 467);
    const auto *gk_468 = buffer.data(gk + 468);
    const auto *gk_469 = buffer.data(gk + 469);
    const auto *gk_470 = buffer.data(gk + 470);
    const auto *gk_471 = buffer.data(gk + 471);
    const auto *gk_472 = buffer.data(gk + 472);
    const auto *gk_473 = buffer.data(gk + 473);
    const auto *gk_474 = buffer.data(gk + 474);
    const auto *gk_475 = buffer.data(gk + 475);
    const auto *gk_476 = buffer.data(gk + 476);
    const auto *gk_477 = buffer.data(gk + 477);

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, dk_100, dk_101, dk_102, dk_103, \
                         dk_104, gk_280, gk_281, gk_282, gk_283, \
                         gk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -dk_100[k]
                   + f_0 * gk_280[k];

        t_173[k] = -dk_101[k]
                   + f_0 * gk_281[k];

        t_174[k] = -dk_102[k]
                   + f_0 * gk_282[k];

        t_175[k] = -dk_103[k]
                   + f_0 * gk_283[k];

        t_176[k] = -dk_104[k]
                   + f_0 * gk_284[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, t_182, dk_105, dk_106, dk_107, \
                         gk_285, gk_286, gk_287, gk_288, gk_289, \
                         gk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -dk_105[k]
                   + f_0 * gk_285[k];

        t_178[k] = -dk_106[k]
                   + f_0 * gk_286[k];

        t_179[k] = -dk_107[k]
                   + f_0 * gk_287[k];

        t_180[k] = f_0 * gk_288[k];

        t_181[k] = f_0 * gk_289[k];

        t_182[k] = f_0 * gk_290[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, t_188, t_189, t_190, gk_291, \
                         gk_292, gk_293, gk_294, gk_295, gk_296, gk_297, \
                         gk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_0 * gk_291[k];

        t_184[k] = f_0 * gk_292[k];

        t_185[k] = f_0 * gk_293[k];

        t_186[k] = f_0 * gk_294[k];

        t_187[k] = f_0 * gk_295[k];

        t_188[k] = f_0 * gk_296[k];

        t_189[k] = f_0 * gk_297[k];

        t_190[k] = f_0 * gk_298[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, t_197, t_198, gk_299, \
                         gk_300, gk_301, gk_302, gk_303, gk_304, gk_305, \
                         gk_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_0 * gk_299[k];

        t_192[k] = f_0 * gk_300[k];

        t_193[k] = f_0 * gk_301[k];

        t_194[k] = f_0 * gk_302[k];

        t_195[k] = f_0 * gk_303[k];

        t_196[k] = f_0 * gk_304[k];

        t_197[k] = f_0 * gk_305[k];

        t_198[k] = f_0 * gk_306[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, t_204, t_205, t_206, gk_307, \
                         gk_308, gk_309, gk_310, gk_311, gk_312, gk_313, \
                         gk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_0 * gk_307[k];

        t_200[k] = f_0 * gk_308[k];

        t_201[k] = f_0 * gk_309[k];

        t_202[k] = f_0 * gk_310[k];

        t_203[k] = f_0 * gk_311[k];

        t_204[k] = f_0 * gk_312[k];

        t_205[k] = f_0 * gk_313[k];

        t_206[k] = f_0 * gk_314[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, t_213, t_214, gk_315, \
                         gk_316, gk_317, gk_318, gk_319, gk_320, gk_321, \
                         gk_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_0 * gk_315[k];

        t_208[k] = f_0 * gk_316[k];

        t_209[k] = f_0 * gk_317[k];

        t_210[k] = f_0 * gk_318[k];

        t_211[k] = f_0 * gk_319[k];

        t_212[k] = f_0 * gk_320[k];

        t_213[k] = f_0 * gk_321[k];

        t_214[k] = f_0 * gk_322[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, dk_108, dk_109, dk_110, dk_111, \
                         gk_323, gk_360, gk_361, gk_362, gk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_0 * gk_323[k];

        t_216[k] = -3.0 * dk_108[k]
                   + f_0 * gk_360[k];

        t_217[k] = -3.0 * dk_109[k]
                   + f_0 * gk_361[k];

        t_218[k] = -3.0 * dk_110[k]
                   + f_0 * gk_362[k];

        t_219[k] = -3.0 * dk_111[k]
                   + f_0 * gk_363[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, dk_112, dk_113, dk_114, dk_115, \
                         dk_116, gk_364, gk_365, gk_366, gk_367, \
                         gk_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -3.0 * dk_112[k]
                   + f_0 * gk_364[k];

        t_221[k] = -3.0 * dk_113[k]
                   + f_0 * gk_365[k];

        t_222[k] = -3.0 * dk_114[k]
                   + f_0 * gk_366[k];

        t_223[k] = -3.0 * dk_115[k]
                   + f_0 * gk_367[k];

        t_224[k] = -3.0 * dk_116[k]
                   + f_0 * gk_368[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, dk_117, dk_118, dk_119, dk_120, \
                         dk_121, gk_369, gk_370, gk_371, gk_372, \
                         gk_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -3.0 * dk_117[k]
                   + f_0 * gk_369[k];

        t_226[k] = -3.0 * dk_118[k]
                   + f_0 * gk_370[k];

        t_227[k] = -3.0 * dk_119[k]
                   + f_0 * gk_371[k];

        t_228[k] = -3.0 * dk_120[k]
                   + f_0 * gk_372[k];

        t_229[k] = -3.0 * dk_121[k]
                   + f_0 * gk_373[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, dk_122, dk_123, dk_124, dk_125, \
                         dk_126, gk_374, gk_375, gk_376, gk_377, \
                         gk_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -3.0 * dk_122[k]
                   + f_0 * gk_374[k];

        t_231[k] = -3.0 * dk_123[k]
                   + f_0 * gk_375[k];

        t_232[k] = -3.0 * dk_124[k]
                   + f_0 * gk_376[k];

        t_233[k] = -3.0 * dk_125[k]
                   + f_0 * gk_377[k];

        t_234[k] = -3.0 * dk_126[k]
                   + f_0 * gk_378[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, dk_127, dk_128, dk_129, dk_130, \
                         dk_131, gk_379, gk_380, gk_381, gk_382, \
                         gk_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -3.0 * dk_127[k]
                   + f_0 * gk_379[k];

        t_236[k] = -3.0 * dk_128[k]
                   + f_0 * gk_380[k];

        t_237[k] = -3.0 * dk_129[k]
                   + f_0 * gk_381[k];

        t_238[k] = -3.0 * dk_130[k]
                   + f_0 * gk_382[k];

        t_239[k] = -3.0 * dk_131[k]
                   + f_0 * gk_383[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, dk_132, dk_133, dk_134, dk_135, \
                         dk_136, gk_384, gk_385, gk_386, gk_387, \
                         gk_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * dk_132[k]
                   + f_0 * gk_384[k];

        t_241[k] = -3.0 * dk_133[k]
                   + f_0 * gk_385[k];

        t_242[k] = -3.0 * dk_134[k]
                   + f_0 * gk_386[k];

        t_243[k] = -3.0 * dk_135[k]
                   + f_0 * gk_387[k];

        t_244[k] = -3.0 * dk_136[k]
                   + f_0 * gk_388[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, dk_137, dk_138, dk_139, dk_140, \
                         dk_141, gk_389, gk_390, gk_391, gk_392, \
                         gk_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * dk_137[k]
                   + f_0 * gk_389[k];

        t_246[k] = -3.0 * dk_138[k]
                   + f_0 * gk_390[k];

        t_247[k] = -3.0 * dk_139[k]
                   + f_0 * gk_391[k];

        t_248[k] = -3.0 * dk_140[k]
                   + f_0 * gk_392[k];

        t_249[k] = -3.0 * dk_141[k]
                   + f_0 * gk_393[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, dk_142, dk_143, dk_144, dk_145, \
                         dk_146, gk_394, gk_395, gk_396, gk_397, \
                         gk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -3.0 * dk_142[k]
                   + f_0 * gk_394[k];

        t_251[k] = -3.0 * dk_143[k]
                   + f_0 * gk_395[k];

        t_252[k] = -2.0 * dk_144[k]
                   + f_0 * gk_396[k];

        t_253[k] = -2.0 * dk_145[k]
                   + f_0 * gk_397[k];

        t_254[k] = -2.0 * dk_146[k]
                   + f_0 * gk_398[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, dk_147, dk_148, dk_149, dk_150, \
                         dk_151, gk_399, gk_400, gk_401, gk_402, \
                         gk_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * dk_147[k]
                   + f_0 * gk_399[k];

        t_256[k] = -2.0 * dk_148[k]
                   + f_0 * gk_400[k];

        t_257[k] = -2.0 * dk_149[k]
                   + f_0 * gk_401[k];

        t_258[k] = -2.0 * dk_150[k]
                   + f_0 * gk_402[k];

        t_259[k] = -2.0 * dk_151[k]
                   + f_0 * gk_403[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, dk_152, dk_153, dk_154, dk_155, \
                         dk_156, gk_404, gk_405, gk_406, gk_407, \
                         gk_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * dk_152[k]
                   + f_0 * gk_404[k];

        t_261[k] = -2.0 * dk_153[k]
                   + f_0 * gk_405[k];

        t_262[k] = -2.0 * dk_154[k]
                   + f_0 * gk_406[k];

        t_263[k] = -2.0 * dk_155[k]
                   + f_0 * gk_407[k];

        t_264[k] = -2.0 * dk_156[k]
                   + f_0 * gk_408[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, dk_157, dk_158, dk_159, dk_160, \
                         dk_161, gk_409, gk_410, gk_411, gk_412, \
                         gk_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * dk_157[k]
                   + f_0 * gk_409[k];

        t_266[k] = -2.0 * dk_158[k]
                   + f_0 * gk_410[k];

        t_267[k] = -2.0 * dk_159[k]
                   + f_0 * gk_411[k];

        t_268[k] = -2.0 * dk_160[k]
                   + f_0 * gk_412[k];

        t_269[k] = -2.0 * dk_161[k]
                   + f_0 * gk_413[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, dk_162, dk_163, dk_164, dk_165, \
                         dk_166, gk_414, gk_415, gk_416, gk_417, \
                         gk_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * dk_162[k]
                   + f_0 * gk_414[k];

        t_271[k] = -2.0 * dk_163[k]
                   + f_0 * gk_415[k];

        t_272[k] = -2.0 * dk_164[k]
                   + f_0 * gk_416[k];

        t_273[k] = -2.0 * dk_165[k]
                   + f_0 * gk_417[k];

        t_274[k] = -2.0 * dk_166[k]
                   + f_0 * gk_418[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, dk_167, dk_168, dk_169, dk_170, \
                         dk_171, gk_419, gk_420, gk_421, gk_422, \
                         gk_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * dk_167[k]
                   + f_0 * gk_419[k];

        t_276[k] = -2.0 * dk_168[k]
                   + f_0 * gk_420[k];

        t_277[k] = -2.0 * dk_169[k]
                   + f_0 * gk_421[k];

        t_278[k] = -2.0 * dk_170[k]
                   + f_0 * gk_422[k];

        t_279[k] = -2.0 * dk_171[k]
                   + f_0 * gk_423[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, dk_172, dk_173, dk_174, dk_175, \
                         dk_176, gk_424, gk_425, gk_426, gk_427, \
                         gk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * dk_172[k]
                   + f_0 * gk_424[k];

        t_281[k] = -2.0 * dk_173[k]
                   + f_0 * gk_425[k];

        t_282[k] = -2.0 * dk_174[k]
                   + f_0 * gk_426[k];

        t_283[k] = -2.0 * dk_175[k]
                   + f_0 * gk_427[k];

        t_284[k] = -2.0 * dk_176[k]
                   + f_0 * gk_428[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, dk_177, dk_178, dk_179, dk_180, \
                         dk_181, gk_429, gk_430, gk_431, gk_432, \
                         gk_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -2.0 * dk_177[k]
                   + f_0 * gk_429[k];

        t_286[k] = -2.0 * dk_178[k]
                   + f_0 * gk_430[k];

        t_287[k] = -2.0 * dk_179[k]
                   + f_0 * gk_431[k];

        t_288[k] = -dk_180[k]
                   + f_0 * gk_432[k];

        t_289[k] = -dk_181[k]
                   + f_0 * gk_433[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, dk_182, dk_183, dk_184, dk_185, \
                         dk_186, gk_434, gk_435, gk_436, gk_437, \
                         gk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -dk_182[k]
                   + f_0 * gk_434[k];

        t_291[k] = -dk_183[k]
                   + f_0 * gk_435[k];

        t_292[k] = -dk_184[k]
                   + f_0 * gk_436[k];

        t_293[k] = -dk_185[k]
                   + f_0 * gk_437[k];

        t_294[k] = -dk_186[k]
                   + f_0 * gk_438[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, dk_187, dk_188, dk_189, dk_190, \
                         dk_191, gk_439, gk_440, gk_441, gk_442, \
                         gk_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -dk_187[k]
                   + f_0 * gk_439[k];

        t_296[k] = -dk_188[k]
                   + f_0 * gk_440[k];

        t_297[k] = -dk_189[k]
                   + f_0 * gk_441[k];

        t_298[k] = -dk_190[k]
                   + f_0 * gk_442[k];

        t_299[k] = -dk_191[k]
                   + f_0 * gk_443[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, dk_192, dk_193, dk_194, dk_195, \
                         dk_196, gk_444, gk_445, gk_446, gk_447, \
                         gk_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -dk_192[k]
                   + f_0 * gk_444[k];

        t_301[k] = -dk_193[k]
                   + f_0 * gk_445[k];

        t_302[k] = -dk_194[k]
                   + f_0 * gk_446[k];

        t_303[k] = -dk_195[k]
                   + f_0 * gk_447[k];

        t_304[k] = -dk_196[k]
                   + f_0 * gk_448[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, dk_197, dk_198, dk_199, dk_200, \
                         dk_201, gk_449, gk_450, gk_451, gk_452, \
                         gk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -dk_197[k]
                   + f_0 * gk_449[k];

        t_306[k] = -dk_198[k]
                   + f_0 * gk_450[k];

        t_307[k] = -dk_199[k]
                   + f_0 * gk_451[k];

        t_308[k] = -dk_200[k]
                   + f_0 * gk_452[k];

        t_309[k] = -dk_201[k]
                   + f_0 * gk_453[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, dk_202, dk_203, dk_204, dk_205, \
                         dk_206, gk_454, gk_455, gk_456, gk_457, \
                         gk_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -dk_202[k]
                   + f_0 * gk_454[k];

        t_311[k] = -dk_203[k]
                   + f_0 * gk_455[k];

        t_312[k] = -dk_204[k]
                   + f_0 * gk_456[k];

        t_313[k] = -dk_205[k]
                   + f_0 * gk_457[k];

        t_314[k] = -dk_206[k]
                   + f_0 * gk_458[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, dk_207, dk_208, dk_209, dk_210, \
                         dk_211, gk_459, gk_460, gk_461, gk_462, \
                         gk_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -dk_207[k]
                   + f_0 * gk_459[k];

        t_316[k] = -dk_208[k]
                   + f_0 * gk_460[k];

        t_317[k] = -dk_209[k]
                   + f_0 * gk_461[k];

        t_318[k] = -dk_210[k]
                   + f_0 * gk_462[k];

        t_319[k] = -dk_211[k]
                   + f_0 * gk_463[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, t_325, dk_212, dk_213, dk_214, \
                         dk_215, gk_464, gk_465, gk_466, gk_467, gk_468, \
                         gk_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -dk_212[k]
                   + f_0 * gk_464[k];

        t_321[k] = -dk_213[k]
                   + f_0 * gk_465[k];

        t_322[k] = -dk_214[k]
                   + f_0 * gk_466[k];

        t_323[k] = -dk_215[k]
                   + f_0 * gk_467[k];

        t_324[k] = f_0 * gk_468[k];

        t_325[k] = f_0 * gk_469[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, t_332, t_333, gk_470, \
                         gk_471, gk_472, gk_473, gk_474, gk_475, gk_476, \
                         gk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_0 * gk_470[k];

        t_327[k] = f_0 * gk_471[k];

        t_328[k] = f_0 * gk_472[k];

        t_329[k] = f_0 * gk_473[k];

        t_330[k] = f_0 * gk_474[k];

        t_331[k] = f_0 * gk_475[k];

        t_332[k] = f_0 * gk_476[k];

        t_333[k] = f_0 * gk_477[k];
    }
}

static auto
compute_prim_geom_10_fk_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gk_478 = buffer.data(gk + 478);
    const auto *gk_479 = buffer.data(gk + 479);
    const auto *gk_480 = buffer.data(gk + 480);
    const auto *gk_481 = buffer.data(gk + 481);
    const auto *gk_482 = buffer.data(gk + 482);
    const auto *gk_483 = buffer.data(gk + 483);
    const auto *gk_484 = buffer.data(gk + 484);
    const auto *gk_485 = buffer.data(gk + 485);
    const auto *gk_486 = buffer.data(gk + 486);
    const auto *gk_487 = buffer.data(gk + 487);
    const auto *gk_488 = buffer.data(gk + 488);
    const auto *gk_489 = buffer.data(gk + 489);
    const auto *gk_490 = buffer.data(gk + 490);
    const auto *gk_491 = buffer.data(gk + 491);
    const auto *gk_492 = buffer.data(gk + 492);
    const auto *gk_493 = buffer.data(gk + 493);
    const auto *gk_494 = buffer.data(gk + 494);
    const auto *gk_495 = buffer.data(gk + 495);
    const auto *gk_496 = buffer.data(gk + 496);
    const auto *gk_497 = buffer.data(gk + 497);
    const auto *gk_498 = buffer.data(gk + 498);
    const auto *gk_499 = buffer.data(gk + 499);
    const auto *gk_500 = buffer.data(gk + 500);
    const auto *gk_501 = buffer.data(gk + 501);
    const auto *gk_502 = buffer.data(gk + 502);
    const auto *gk_503 = buffer.data(gk + 503);

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, t_339, t_340, t_341, gk_478, \
                         gk_479, gk_480, gk_481, gk_482, gk_483, gk_484, \
                         gk_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_0 * gk_478[k];

        t_335[k] = f_0 * gk_479[k];

        t_336[k] = f_0 * gk_480[k];

        t_337[k] = f_0 * gk_481[k];

        t_338[k] = f_0 * gk_482[k];

        t_339[k] = f_0 * gk_483[k];

        t_340[k] = f_0 * gk_484[k];

        t_341[k] = f_0 * gk_485[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, t_347, t_348, t_349, gk_486, \
                         gk_487, gk_488, gk_489, gk_490, gk_491, gk_492, \
                         gk_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_0 * gk_486[k];

        t_343[k] = f_0 * gk_487[k];

        t_344[k] = f_0 * gk_488[k];

        t_345[k] = f_0 * gk_489[k];

        t_346[k] = f_0 * gk_490[k];

        t_347[k] = f_0 * gk_491[k];

        t_348[k] = f_0 * gk_492[k];

        t_349[k] = f_0 * gk_493[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, t_355, t_356, t_357, gk_494, \
                         gk_495, gk_496, gk_497, gk_498, gk_499, gk_500, \
                         gk_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_0 * gk_494[k];

        t_351[k] = f_0 * gk_495[k];

        t_352[k] = f_0 * gk_496[k];

        t_353[k] = f_0 * gk_497[k];

        t_354[k] = f_0 * gk_498[k];

        t_355[k] = f_0 * gk_499[k];

        t_356[k] = f_0 * gk_500[k];

        t_357[k] = f_0 * gk_501[k];
    }

#pragma omp simd aligned(t_358, t_359, gk_502, gk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_0 * gk_502[k];

        t_359[k] = f_0 * gk_503[k];
    }
}

auto
compute_prim_geom_10_fk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t dk, const size_t gk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fk_electron_repulsion_1_piece0(buffer, target, dk, gk, ncols, alpha);

    compute_prim_geom_10_fk_electron_repulsion_1_piece1(buffer, target, dk, gk, ncols, alpha);

    compute_prim_geom_10_fk_electron_repulsion_1_piece2(buffer, target, gk, ncols, alpha);
}

static auto
compute_prim_geom_10_fk_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dk, const size_t gk,
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
    auto *t_181 = buffer.data(target + 181);
    auto *t_182 = buffer.data(target + 182);
    auto *t_183 = buffer.data(target + 183);
    auto *t_184 = buffer.data(target + 184);

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_6 = buffer.data(dk + 6);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_9 = buffer.data(dk + 9);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_13 = buffer.data(dk + 13);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_22 = buffer.data(dk + 22);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_25 = buffer.data(dk + 25);
    const auto *dk_26 = buffer.data(dk + 26);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_29 = buffer.data(dk + 29);
    const auto *dk_30 = buffer.data(dk + 30);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_33 = buffer.data(dk + 33);
    const auto *dk_34 = buffer.data(dk + 34);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_36 = buffer.data(dk + 36);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_49 = buffer.data(dk + 49);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_53 = buffer.data(dk + 53);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_60 = buffer.data(dk + 60);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_62 = buffer.data(dk + 62);
    const auto *dk_63 = buffer.data(dk + 63);
    const auto *dk_64 = buffer.data(dk + 64);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_71 = buffer.data(dk + 71);
    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_76 = buffer.data(dk + 76);

    const auto *gk_72 = buffer.data(gk + 72);
    const auto *gk_73 = buffer.data(gk + 73);
    const auto *gk_74 = buffer.data(gk + 74);
    const auto *gk_75 = buffer.data(gk + 75);
    const auto *gk_76 = buffer.data(gk + 76);
    const auto *gk_77 = buffer.data(gk + 77);
    const auto *gk_78 = buffer.data(gk + 78);
    const auto *gk_79 = buffer.data(gk + 79);
    const auto *gk_80 = buffer.data(gk + 80);
    const auto *gk_81 = buffer.data(gk + 81);
    const auto *gk_82 = buffer.data(gk + 82);
    const auto *gk_83 = buffer.data(gk + 83);
    const auto *gk_84 = buffer.data(gk + 84);
    const auto *gk_85 = buffer.data(gk + 85);
    const auto *gk_86 = buffer.data(gk + 86);
    const auto *gk_87 = buffer.data(gk + 87);
    const auto *gk_88 = buffer.data(gk + 88);
    const auto *gk_89 = buffer.data(gk + 89);
    const auto *gk_90 = buffer.data(gk + 90);
    const auto *gk_91 = buffer.data(gk + 91);
    const auto *gk_92 = buffer.data(gk + 92);
    const auto *gk_93 = buffer.data(gk + 93);
    const auto *gk_94 = buffer.data(gk + 94);
    const auto *gk_95 = buffer.data(gk + 95);
    const auto *gk_96 = buffer.data(gk + 96);
    const auto *gk_97 = buffer.data(gk + 97);
    const auto *gk_98 = buffer.data(gk + 98);
    const auto *gk_99 = buffer.data(gk + 99);
    const auto *gk_100 = buffer.data(gk + 100);
    const auto *gk_101 = buffer.data(gk + 101);
    const auto *gk_102 = buffer.data(gk + 102);
    const auto *gk_103 = buffer.data(gk + 103);
    const auto *gk_104 = buffer.data(gk + 104);
    const auto *gk_105 = buffer.data(gk + 105);
    const auto *gk_106 = buffer.data(gk + 106);
    const auto *gk_107 = buffer.data(gk + 107);
    const auto *gk_144 = buffer.data(gk + 144);
    const auto *gk_145 = buffer.data(gk + 145);
    const auto *gk_146 = buffer.data(gk + 146);
    const auto *gk_147 = buffer.data(gk + 147);
    const auto *gk_148 = buffer.data(gk + 148);
    const auto *gk_149 = buffer.data(gk + 149);
    const auto *gk_150 = buffer.data(gk + 150);
    const auto *gk_151 = buffer.data(gk + 151);
    const auto *gk_152 = buffer.data(gk + 152);
    const auto *gk_153 = buffer.data(gk + 153);
    const auto *gk_154 = buffer.data(gk + 154);
    const auto *gk_155 = buffer.data(gk + 155);
    const auto *gk_156 = buffer.data(gk + 156);
    const auto *gk_157 = buffer.data(gk + 157);
    const auto *gk_158 = buffer.data(gk + 158);
    const auto *gk_159 = buffer.data(gk + 159);
    const auto *gk_160 = buffer.data(gk + 160);
    const auto *gk_161 = buffer.data(gk + 161);
    const auto *gk_162 = buffer.data(gk + 162);
    const auto *gk_163 = buffer.data(gk + 163);
    const auto *gk_164 = buffer.data(gk + 164);
    const auto *gk_165 = buffer.data(gk + 165);
    const auto *gk_166 = buffer.data(gk + 166);
    const auto *gk_167 = buffer.data(gk + 167);
    const auto *gk_168 = buffer.data(gk + 168);
    const auto *gk_169 = buffer.data(gk + 169);
    const auto *gk_170 = buffer.data(gk + 170);
    const auto *gk_171 = buffer.data(gk + 171);
    const auto *gk_172 = buffer.data(gk + 172);
    const auto *gk_173 = buffer.data(gk + 173);
    const auto *gk_174 = buffer.data(gk + 174);
    const auto *gk_175 = buffer.data(gk + 175);
    const auto *gk_176 = buffer.data(gk + 176);
    const auto *gk_177 = buffer.data(gk + 177);
    const auto *gk_178 = buffer.data(gk + 178);
    const auto *gk_179 = buffer.data(gk + 179);
    const auto *gk_180 = buffer.data(gk + 180);
    const auto *gk_181 = buffer.data(gk + 181);
    const auto *gk_182 = buffer.data(gk + 182);
    const auto *gk_183 = buffer.data(gk + 183);
    const auto *gk_184 = buffer.data(gk + 184);
    const auto *gk_185 = buffer.data(gk + 185);
    const auto *gk_186 = buffer.data(gk + 186);
    const auto *gk_187 = buffer.data(gk + 187);
    const auto *gk_188 = buffer.data(gk + 188);
    const auto *gk_189 = buffer.data(gk + 189);
    const auto *gk_190 = buffer.data(gk + 190);
    const auto *gk_191 = buffer.data(gk + 191);
    const auto *gk_192 = buffer.data(gk + 192);
    const auto *gk_193 = buffer.data(gk + 193);
    const auto *gk_194 = buffer.data(gk + 194);
    const auto *gk_195 = buffer.data(gk + 195);
    const auto *gk_196 = buffer.data(gk + 196);
    const auto *gk_197 = buffer.data(gk + 197);
    const auto *gk_198 = buffer.data(gk + 198);
    const auto *gk_199 = buffer.data(gk + 199);
    const auto *gk_200 = buffer.data(gk + 200);
    const auto *gk_201 = buffer.data(gk + 201);
    const auto *gk_202 = buffer.data(gk + 202);
    const auto *gk_203 = buffer.data(gk + 203);
    const auto *gk_204 = buffer.data(gk + 204);
    const auto *gk_205 = buffer.data(gk + 205);
    const auto *gk_206 = buffer.data(gk + 206);
    const auto *gk_207 = buffer.data(gk + 207);
    const auto *gk_208 = buffer.data(gk + 208);
    const auto *gk_209 = buffer.data(gk + 209);
    const auto *gk_210 = buffer.data(gk + 210);
    const auto *gk_211 = buffer.data(gk + 211);
    const auto *gk_212 = buffer.data(gk + 212);
    const auto *gk_213 = buffer.data(gk + 213);
    const auto *gk_214 = buffer.data(gk + 214);
    const auto *gk_215 = buffer.data(gk + 215);
    const auto *gk_252 = buffer.data(gk + 252);
    const auto *gk_253 = buffer.data(gk + 253);
    const auto *gk_254 = buffer.data(gk + 254);
    const auto *gk_255 = buffer.data(gk + 255);
    const auto *gk_256 = buffer.data(gk + 256);
    const auto *gk_257 = buffer.data(gk + 257);
    const auto *gk_258 = buffer.data(gk + 258);
    const auto *gk_259 = buffer.data(gk + 259);
    const auto *gk_260 = buffer.data(gk + 260);
    const auto *gk_261 = buffer.data(gk + 261);
    const auto *gk_262 = buffer.data(gk + 262);
    const auto *gk_263 = buffer.data(gk + 263);
    const auto *gk_264 = buffer.data(gk + 264);
    const auto *gk_265 = buffer.data(gk + 265);
    const auto *gk_266 = buffer.data(gk + 266);
    const auto *gk_267 = buffer.data(gk + 267);
    const auto *gk_268 = buffer.data(gk + 268);
    const auto *gk_269 = buffer.data(gk + 269);
    const auto *gk_270 = buffer.data(gk + 270);
    const auto *gk_271 = buffer.data(gk + 271);
    const auto *gk_272 = buffer.data(gk + 272);
    const auto *gk_273 = buffer.data(gk + 273);
    const auto *gk_274 = buffer.data(gk + 274);
    const auto *gk_275 = buffer.data(gk + 275);
    const auto *gk_276 = buffer.data(gk + 276);
    const auto *gk_277 = buffer.data(gk + 277);
    const auto *gk_278 = buffer.data(gk + 278);
    const auto *gk_279 = buffer.data(gk + 279);
    const auto *gk_280 = buffer.data(gk + 280);
    const auto *gk_281 = buffer.data(gk + 281);
    const auto *gk_282 = buffer.data(gk + 282);
    const auto *gk_283 = buffer.data(gk + 283);
    const auto *gk_284 = buffer.data(gk + 284);
    const auto *gk_285 = buffer.data(gk + 285);
    const auto *gk_286 = buffer.data(gk + 286);
    const auto *gk_287 = buffer.data(gk + 287);
    const auto *gk_288 = buffer.data(gk + 288);
    const auto *gk_289 = buffer.data(gk + 289);
    const auto *gk_290 = buffer.data(gk + 290);
    const auto *gk_291 = buffer.data(gk + 291);
    const auto *gk_292 = buffer.data(gk + 292);
    const auto *gk_293 = buffer.data(gk + 293);
    const auto *gk_294 = buffer.data(gk + 294);
    const auto *gk_295 = buffer.data(gk + 295);
    const auto *gk_296 = buffer.data(gk + 296);
    const auto *gk_297 = buffer.data(gk + 297);
    const auto *gk_298 = buffer.data(gk + 298);
    const auto *gk_299 = buffer.data(gk + 299);
    const auto *gk_300 = buffer.data(gk + 300);
    const auto *gk_301 = buffer.data(gk + 301);
    const auto *gk_302 = buffer.data(gk + 302);
    const auto *gk_303 = buffer.data(gk + 303);
    const auto *gk_304 = buffer.data(gk + 304);
    const auto *gk_305 = buffer.data(gk + 305);
    const auto *gk_306 = buffer.data(gk + 306);
    const auto *gk_307 = buffer.data(gk + 307);
    const auto *gk_308 = buffer.data(gk + 308);
    const auto *gk_309 = buffer.data(gk + 309);
    const auto *gk_310 = buffer.data(gk + 310);
    const auto *gk_311 = buffer.data(gk + 311);
    const auto *gk_312 = buffer.data(gk + 312);
    const auto *gk_313 = buffer.data(gk + 313);
    const auto *gk_314 = buffer.data(gk + 314);
    const auto *gk_315 = buffer.data(gk + 315);
    const auto *gk_316 = buffer.data(gk + 316);
    const auto *gk_317 = buffer.data(gk + 317);
    const auto *gk_318 = buffer.data(gk + 318);
    const auto *gk_319 = buffer.data(gk + 319);
    const auto *gk_320 = buffer.data(gk + 320);
    const auto *gk_321 = buffer.data(gk + 321);
    const auto *gk_322 = buffer.data(gk + 322);
    const auto *gk_323 = buffer.data(gk + 323);
    const auto *gk_324 = buffer.data(gk + 324);
    const auto *gk_325 = buffer.data(gk + 325);
    const auto *gk_326 = buffer.data(gk + 326);
    const auto *gk_327 = buffer.data(gk + 327);
    const auto *gk_328 = buffer.data(gk + 328);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gk_72, gk_73, gk_74, gk_75, \
                         gk_76, gk_77, gk_78, gk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gk_72[k];

        t_1[k] = f_0 * gk_73[k];

        t_2[k] = f_0 * gk_74[k];

        t_3[k] = f_0 * gk_75[k];

        t_4[k] = f_0 * gk_76[k];

        t_5[k] = f_0 * gk_77[k];

        t_6[k] = f_0 * gk_78[k];

        t_7[k] = f_0 * gk_79[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, gk_80, gk_81, gk_82, \
                         gk_83, gk_84, gk_85, gk_86, gk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gk_80[k];

        t_9[k] = f_0 * gk_81[k];

        t_10[k] = f_0 * gk_82[k];

        t_11[k] = f_0 * gk_83[k];

        t_12[k] = f_0 * gk_84[k];

        t_13[k] = f_0 * gk_85[k];

        t_14[k] = f_0 * gk_86[k];

        t_15[k] = f_0 * gk_87[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, gk_88, gk_89, gk_90, \
                         gk_91, gk_92, gk_93, gk_94, gk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * gk_88[k];

        t_17[k] = f_0 * gk_89[k];

        t_18[k] = f_0 * gk_90[k];

        t_19[k] = f_0 * gk_91[k];

        t_20[k] = f_0 * gk_92[k];

        t_21[k] = f_0 * gk_93[k];

        t_22[k] = f_0 * gk_94[k];

        t_23[k] = f_0 * gk_95[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, gk_96, gk_97, gk_98, \
                         gk_99, gk_100, gk_101, gk_102, gk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * gk_96[k];

        t_25[k] = f_0 * gk_97[k];

        t_26[k] = f_0 * gk_98[k];

        t_27[k] = f_0 * gk_99[k];

        t_28[k] = f_0 * gk_100[k];

        t_29[k] = f_0 * gk_101[k];

        t_30[k] = f_0 * gk_102[k];

        t_31[k] = f_0 * gk_103[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, gk_104, gk_105, \
                         gk_106, gk_107, gk_144, gk_145, gk_146, \
                         gk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * gk_104[k];

        t_33[k] = f_0 * gk_105[k];

        t_34[k] = f_0 * gk_106[k];

        t_35[k] = f_0 * gk_107[k];

        t_36[k] = f_0 * gk_144[k];

        t_37[k] = f_0 * gk_145[k];

        t_38[k] = f_0 * gk_146[k];

        t_39[k] = f_0 * gk_147[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, gk_148, gk_149, \
                         gk_150, gk_151, gk_152, gk_153, gk_154, \
                         gk_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * gk_148[k];

        t_41[k] = f_0 * gk_149[k];

        t_42[k] = f_0 * gk_150[k];

        t_43[k] = f_0 * gk_151[k];

        t_44[k] = f_0 * gk_152[k];

        t_45[k] = f_0 * gk_153[k];

        t_46[k] = f_0 * gk_154[k];

        t_47[k] = f_0 * gk_155[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, gk_156, gk_157, \
                         gk_158, gk_159, gk_160, gk_161, gk_162, \
                         gk_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * gk_156[k];

        t_49[k] = f_0 * gk_157[k];

        t_50[k] = f_0 * gk_158[k];

        t_51[k] = f_0 * gk_159[k];

        t_52[k] = f_0 * gk_160[k];

        t_53[k] = f_0 * gk_161[k];

        t_54[k] = f_0 * gk_162[k];

        t_55[k] = f_0 * gk_163[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, gk_164, gk_165, \
                         gk_166, gk_167, gk_168, gk_169, gk_170, \
                         gk_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * gk_164[k];

        t_57[k] = f_0 * gk_165[k];

        t_58[k] = f_0 * gk_166[k];

        t_59[k] = f_0 * gk_167[k];

        t_60[k] = f_0 * gk_168[k];

        t_61[k] = f_0 * gk_169[k];

        t_62[k] = f_0 * gk_170[k];

        t_63[k] = f_0 * gk_171[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, gk_172, gk_173, \
                         gk_174, gk_175, gk_176, gk_177, gk_178, \
                         gk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * gk_172[k];

        t_65[k] = f_0 * gk_173[k];

        t_66[k] = f_0 * gk_174[k];

        t_67[k] = f_0 * gk_175[k];

        t_68[k] = f_0 * gk_176[k];

        t_69[k] = f_0 * gk_177[k];

        t_70[k] = f_0 * gk_178[k];

        t_71[k] = f_0 * gk_179[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, dk_0, dk_1, dk_2, dk_3, dk_4, gk_180, \
                         gk_181, gk_182, gk_183, gk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -dk_0[k]
                  + f_0 * gk_180[k];

        t_73[k] = -dk_1[k]
                  + f_0 * gk_181[k];

        t_74[k] = -dk_2[k]
                  + f_0 * gk_182[k];

        t_75[k] = -dk_3[k]
                  + f_0 * gk_183[k];

        t_76[k] = -dk_4[k]
                  + f_0 * gk_184[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, dk_5, dk_6, dk_7, dk_8, dk_9, gk_185, \
                         gk_186, gk_187, gk_188, gk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -dk_5[k]
                  + f_0 * gk_185[k];

        t_78[k] = -dk_6[k]
                  + f_0 * gk_186[k];

        t_79[k] = -dk_7[k]
                  + f_0 * gk_187[k];

        t_80[k] = -dk_8[k]
                  + f_0 * gk_188[k];

        t_81[k] = -dk_9[k]
                  + f_0 * gk_189[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, dk_10, dk_11, dk_12, dk_13, dk_14, \
                         gk_190, gk_191, gk_192, gk_193, gk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -dk_10[k]
                  + f_0 * gk_190[k];

        t_83[k] = -dk_11[k]
                  + f_0 * gk_191[k];

        t_84[k] = -dk_12[k]
                  + f_0 * gk_192[k];

        t_85[k] = -dk_13[k]
                  + f_0 * gk_193[k];

        t_86[k] = -dk_14[k]
                  + f_0 * gk_194[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, dk_15, dk_16, dk_17, dk_18, dk_19, \
                         gk_195, gk_196, gk_197, gk_198, gk_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -dk_15[k]
                  + f_0 * gk_195[k];

        t_88[k] = -dk_16[k]
                  + f_0 * gk_196[k];

        t_89[k] = -dk_17[k]
                  + f_0 * gk_197[k];

        t_90[k] = -dk_18[k]
                  + f_0 * gk_198[k];

        t_91[k] = -dk_19[k]
                  + f_0 * gk_199[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, dk_20, dk_21, dk_22, dk_23, dk_24, \
                         gk_200, gk_201, gk_202, gk_203, gk_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = -dk_20[k]
                  + f_0 * gk_200[k];

        t_93[k] = -dk_21[k]
                  + f_0 * gk_201[k];

        t_94[k] = -dk_22[k]
                  + f_0 * gk_202[k];

        t_95[k] = -dk_23[k]
                  + f_0 * gk_203[k];

        t_96[k] = -dk_24[k]
                  + f_0 * gk_204[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, dk_25, dk_26, dk_27, dk_28, dk_29, \
                         gk_205, gk_206, gk_207, gk_208, gk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -dk_25[k]
                  + f_0 * gk_205[k];

        t_98[k] = -dk_26[k]
                  + f_0 * gk_206[k];

        t_99[k] = -dk_27[k]
                  + f_0 * gk_207[k];

        t_100[k] = -dk_28[k]
                   + f_0 * gk_208[k];

        t_101[k] = -dk_29[k]
                   + f_0 * gk_209[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, dk_30, dk_31, dk_32, dk_33, dk_34, \
                         gk_210, gk_211, gk_212, gk_213, gk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -dk_30[k]
                   + f_0 * gk_210[k];

        t_103[k] = -dk_31[k]
                   + f_0 * gk_211[k];

        t_104[k] = -dk_32[k]
                   + f_0 * gk_212[k];

        t_105[k] = -dk_33[k]
                   + f_0 * gk_213[k];

        t_106[k] = -dk_34[k]
                   + f_0 * gk_214[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, dk_35, gk_215, \
                         gk_252, gk_253, gk_254, gk_255, gk_256, \
                         gk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -dk_35[k]
                   + f_0 * gk_215[k];

        t_108[k] = f_0 * gk_252[k];

        t_109[k] = f_0 * gk_253[k];

        t_110[k] = f_0 * gk_254[k];

        t_111[k] = f_0 * gk_255[k];

        t_112[k] = f_0 * gk_256[k];

        t_113[k] = f_0 * gk_257[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, t_120, t_121, gk_258, \
                         gk_259, gk_260, gk_261, gk_262, gk_263, gk_264, \
                         gk_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * gk_258[k];

        t_115[k] = f_0 * gk_259[k];

        t_116[k] = f_0 * gk_260[k];

        t_117[k] = f_0 * gk_261[k];

        t_118[k] = f_0 * gk_262[k];

        t_119[k] = f_0 * gk_263[k];

        t_120[k] = f_0 * gk_264[k];

        t_121[k] = f_0 * gk_265[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, t_128, t_129, gk_266, \
                         gk_267, gk_268, gk_269, gk_270, gk_271, gk_272, \
                         gk_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * gk_266[k];

        t_123[k] = f_0 * gk_267[k];

        t_124[k] = f_0 * gk_268[k];

        t_125[k] = f_0 * gk_269[k];

        t_126[k] = f_0 * gk_270[k];

        t_127[k] = f_0 * gk_271[k];

        t_128[k] = f_0 * gk_272[k];

        t_129[k] = f_0 * gk_273[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, t_135, t_136, t_137, gk_274, \
                         gk_275, gk_276, gk_277, gk_278, gk_279, gk_280, \
                         gk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * gk_274[k];

        t_131[k] = f_0 * gk_275[k];

        t_132[k] = f_0 * gk_276[k];

        t_133[k] = f_0 * gk_277[k];

        t_134[k] = f_0 * gk_278[k];

        t_135[k] = f_0 * gk_279[k];

        t_136[k] = f_0 * gk_280[k];

        t_137[k] = f_0 * gk_281[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, t_144, dk_36, gk_282, \
                         gk_283, gk_284, gk_285, gk_286, gk_287, \
                         gk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * gk_282[k];

        t_139[k] = f_0 * gk_283[k];

        t_140[k] = f_0 * gk_284[k];

        t_141[k] = f_0 * gk_285[k];

        t_142[k] = f_0 * gk_286[k];

        t_143[k] = f_0 * gk_287[k];

        t_144[k] = -dk_36[k]
                   + f_0 * gk_288[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, dk_37, dk_38, dk_39, dk_40, dk_41, \
                         gk_289, gk_290, gk_291, gk_292, gk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -dk_37[k]
                   + f_0 * gk_289[k];

        t_146[k] = -dk_38[k]
                   + f_0 * gk_290[k];

        t_147[k] = -dk_39[k]
                   + f_0 * gk_291[k];

        t_148[k] = -dk_40[k]
                   + f_0 * gk_292[k];

        t_149[k] = -dk_41[k]
                   + f_0 * gk_293[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, dk_42, dk_43, dk_44, dk_45, dk_46, \
                         gk_294, gk_295, gk_296, gk_297, gk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -dk_42[k]
                   + f_0 * gk_294[k];

        t_151[k] = -dk_43[k]
                   + f_0 * gk_295[k];

        t_152[k] = -dk_44[k]
                   + f_0 * gk_296[k];

        t_153[k] = -dk_45[k]
                   + f_0 * gk_297[k];

        t_154[k] = -dk_46[k]
                   + f_0 * gk_298[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, dk_47, dk_48, dk_49, dk_50, dk_51, \
                         gk_299, gk_300, gk_301, gk_302, gk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -dk_47[k]
                   + f_0 * gk_299[k];

        t_156[k] = -dk_48[k]
                   + f_0 * gk_300[k];

        t_157[k] = -dk_49[k]
                   + f_0 * gk_301[k];

        t_158[k] = -dk_50[k]
                   + f_0 * gk_302[k];

        t_159[k] = -dk_51[k]
                   + f_0 * gk_303[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, dk_52, dk_53, dk_54, dk_55, dk_56, \
                         gk_304, gk_305, gk_306, gk_307, gk_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -dk_52[k]
                   + f_0 * gk_304[k];

        t_161[k] = -dk_53[k]
                   + f_0 * gk_305[k];

        t_162[k] = -dk_54[k]
                   + f_0 * gk_306[k];

        t_163[k] = -dk_55[k]
                   + f_0 * gk_307[k];

        t_164[k] = -dk_56[k]
                   + f_0 * gk_308[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, dk_57, dk_58, dk_59, dk_60, dk_61, \
                         gk_309, gk_310, gk_311, gk_312, gk_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -dk_57[k]
                   + f_0 * gk_309[k];

        t_166[k] = -dk_58[k]
                   + f_0 * gk_310[k];

        t_167[k] = -dk_59[k]
                   + f_0 * gk_311[k];

        t_168[k] = -dk_60[k]
                   + f_0 * gk_312[k];

        t_169[k] = -dk_61[k]
                   + f_0 * gk_313[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, dk_62, dk_63, dk_64, dk_65, dk_66, \
                         gk_314, gk_315, gk_316, gk_317, gk_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -dk_62[k]
                   + f_0 * gk_314[k];

        t_171[k] = -dk_63[k]
                   + f_0 * gk_315[k];

        t_172[k] = -dk_64[k]
                   + f_0 * gk_316[k];

        t_173[k] = -dk_65[k]
                   + f_0 * gk_317[k];

        t_174[k] = -dk_66[k]
                   + f_0 * gk_318[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, dk_67, dk_68, dk_69, dk_70, dk_71, \
                         gk_319, gk_320, gk_321, gk_322, gk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -dk_67[k]
                   + f_0 * gk_319[k];

        t_176[k] = -dk_68[k]
                   + f_0 * gk_320[k];

        t_177[k] = -dk_69[k]
                   + f_0 * gk_321[k];

        t_178[k] = -dk_70[k]
                   + f_0 * gk_322[k];

        t_179[k] = -dk_71[k]
                   + f_0 * gk_323[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, dk_72, dk_73, dk_74, dk_75, dk_76, \
                         gk_324, gk_325, gk_326, gk_327, gk_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * dk_72[k]
                   + f_0 * gk_324[k];

        t_181[k] = -2.0 * dk_73[k]
                   + f_0 * gk_325[k];

        t_182[k] = -2.0 * dk_74[k]
                   + f_0 * gk_326[k];

        t_183[k] = -2.0 * dk_75[k]
                   + f_0 * gk_327[k];

        t_184[k] = -2.0 * dk_76[k]
                   + f_0 * gk_328[k];
    }
}

static auto
compute_prim_geom_10_fk_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dk, const size_t gk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_347 = buffer.data(target + 347);

    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_80 = buffer.data(dk + 80);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_83 = buffer.data(dk + 83);
    const auto *dk_84 = buffer.data(dk + 84);
    const auto *dk_85 = buffer.data(dk + 85);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_87 = buffer.data(dk + 87);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_89 = buffer.data(dk + 89);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_94 = buffer.data(dk + 94);
    const auto *dk_95 = buffer.data(dk + 95);
    const auto *dk_96 = buffer.data(dk + 96);
    const auto *dk_97 = buffer.data(dk + 97);
    const auto *dk_98 = buffer.data(dk + 98);
    const auto *dk_99 = buffer.data(dk + 99);
    const auto *dk_100 = buffer.data(dk + 100);
    const auto *dk_101 = buffer.data(dk + 101);
    const auto *dk_102 = buffer.data(dk + 102);
    const auto *dk_103 = buffer.data(dk + 103);
    const auto *dk_104 = buffer.data(dk + 104);
    const auto *dk_105 = buffer.data(dk + 105);
    const auto *dk_106 = buffer.data(dk + 106);
    const auto *dk_107 = buffer.data(dk + 107);
    const auto *dk_108 = buffer.data(dk + 108);
    const auto *dk_109 = buffer.data(dk + 109);
    const auto *dk_110 = buffer.data(dk + 110);
    const auto *dk_111 = buffer.data(dk + 111);
    const auto *dk_112 = buffer.data(dk + 112);
    const auto *dk_113 = buffer.data(dk + 113);
    const auto *dk_114 = buffer.data(dk + 114);
    const auto *dk_115 = buffer.data(dk + 115);
    const auto *dk_116 = buffer.data(dk + 116);
    const auto *dk_117 = buffer.data(dk + 117);
    const auto *dk_118 = buffer.data(dk + 118);
    const auto *dk_119 = buffer.data(dk + 119);
    const auto *dk_120 = buffer.data(dk + 120);
    const auto *dk_121 = buffer.data(dk + 121);
    const auto *dk_122 = buffer.data(dk + 122);
    const auto *dk_123 = buffer.data(dk + 123);
    const auto *dk_124 = buffer.data(dk + 124);
    const auto *dk_125 = buffer.data(dk + 125);
    const auto *dk_126 = buffer.data(dk + 126);
    const auto *dk_127 = buffer.data(dk + 127);
    const auto *dk_128 = buffer.data(dk + 128);
    const auto *dk_129 = buffer.data(dk + 129);
    const auto *dk_130 = buffer.data(dk + 130);
    const auto *dk_131 = buffer.data(dk + 131);
    const auto *dk_132 = buffer.data(dk + 132);
    const auto *dk_133 = buffer.data(dk + 133);
    const auto *dk_134 = buffer.data(dk + 134);
    const auto *dk_135 = buffer.data(dk + 135);
    const auto *dk_136 = buffer.data(dk + 136);
    const auto *dk_137 = buffer.data(dk + 137);
    const auto *dk_138 = buffer.data(dk + 138);
    const auto *dk_139 = buffer.data(dk + 139);
    const auto *dk_140 = buffer.data(dk + 140);
    const auto *dk_141 = buffer.data(dk + 141);
    const auto *dk_142 = buffer.data(dk + 142);
    const auto *dk_143 = buffer.data(dk + 143);
    const auto *dk_144 = buffer.data(dk + 144);
    const auto *dk_145 = buffer.data(dk + 145);
    const auto *dk_146 = buffer.data(dk + 146);
    const auto *dk_147 = buffer.data(dk + 147);
    const auto *dk_148 = buffer.data(dk + 148);
    const auto *dk_149 = buffer.data(dk + 149);
    const auto *dk_150 = buffer.data(dk + 150);
    const auto *dk_151 = buffer.data(dk + 151);
    const auto *dk_152 = buffer.data(dk + 152);
    const auto *dk_153 = buffer.data(dk + 153);
    const auto *dk_154 = buffer.data(dk + 154);
    const auto *dk_155 = buffer.data(dk + 155);
    const auto *dk_156 = buffer.data(dk + 156);
    const auto *dk_157 = buffer.data(dk + 157);
    const auto *dk_158 = buffer.data(dk + 158);
    const auto *dk_159 = buffer.data(dk + 159);
    const auto *dk_160 = buffer.data(dk + 160);
    const auto *dk_161 = buffer.data(dk + 161);
    const auto *dk_162 = buffer.data(dk + 162);
    const auto *dk_163 = buffer.data(dk + 163);
    const auto *dk_164 = buffer.data(dk + 164);
    const auto *dk_165 = buffer.data(dk + 165);
    const auto *dk_166 = buffer.data(dk + 166);
    const auto *dk_167 = buffer.data(dk + 167);
    const auto *dk_168 = buffer.data(dk + 168);
    const auto *dk_169 = buffer.data(dk + 169);
    const auto *dk_170 = buffer.data(dk + 170);
    const auto *dk_171 = buffer.data(dk + 171);
    const auto *dk_172 = buffer.data(dk + 172);
    const auto *dk_173 = buffer.data(dk + 173);
    const auto *dk_174 = buffer.data(dk + 174);
    const auto *dk_175 = buffer.data(dk + 175);
    const auto *dk_176 = buffer.data(dk + 176);
    const auto *dk_177 = buffer.data(dk + 177);
    const auto *dk_178 = buffer.data(dk + 178);
    const auto *dk_179 = buffer.data(dk + 179);
    const auto *dk_180 = buffer.data(dk + 180);
    const auto *dk_181 = buffer.data(dk + 181);
    const auto *dk_182 = buffer.data(dk + 182);
    const auto *dk_183 = buffer.data(dk + 183);
    const auto *dk_184 = buffer.data(dk + 184);
    const auto *dk_185 = buffer.data(dk + 185);
    const auto *dk_186 = buffer.data(dk + 186);
    const auto *dk_187 = buffer.data(dk + 187);
    const auto *dk_188 = buffer.data(dk + 188);
    const auto *dk_189 = buffer.data(dk + 189);
    const auto *dk_190 = buffer.data(dk + 190);
    const auto *dk_191 = buffer.data(dk + 191);
    const auto *dk_192 = buffer.data(dk + 192);
    const auto *dk_193 = buffer.data(dk + 193);
    const auto *dk_194 = buffer.data(dk + 194);
    const auto *dk_195 = buffer.data(dk + 195);
    const auto *dk_196 = buffer.data(dk + 196);
    const auto *dk_197 = buffer.data(dk + 197);
    const auto *dk_198 = buffer.data(dk + 198);
    const auto *dk_199 = buffer.data(dk + 199);
    const auto *dk_200 = buffer.data(dk + 200);
    const auto *dk_201 = buffer.data(dk + 201);
    const auto *dk_202 = buffer.data(dk + 202);
    const auto *dk_203 = buffer.data(dk + 203);

    const auto *gk_329 = buffer.data(gk + 329);
    const auto *gk_330 = buffer.data(gk + 330);
    const auto *gk_331 = buffer.data(gk + 331);
    const auto *gk_332 = buffer.data(gk + 332);
    const auto *gk_333 = buffer.data(gk + 333);
    const auto *gk_334 = buffer.data(gk + 334);
    const auto *gk_335 = buffer.data(gk + 335);
    const auto *gk_336 = buffer.data(gk + 336);
    const auto *gk_337 = buffer.data(gk + 337);
    const auto *gk_338 = buffer.data(gk + 338);
    const auto *gk_339 = buffer.data(gk + 339);
    const auto *gk_340 = buffer.data(gk + 340);
    const auto *gk_341 = buffer.data(gk + 341);
    const auto *gk_342 = buffer.data(gk + 342);
    const auto *gk_343 = buffer.data(gk + 343);
    const auto *gk_344 = buffer.data(gk + 344);
    const auto *gk_345 = buffer.data(gk + 345);
    const auto *gk_346 = buffer.data(gk + 346);
    const auto *gk_347 = buffer.data(gk + 347);
    const auto *gk_348 = buffer.data(gk + 348);
    const auto *gk_349 = buffer.data(gk + 349);
    const auto *gk_350 = buffer.data(gk + 350);
    const auto *gk_351 = buffer.data(gk + 351);
    const auto *gk_352 = buffer.data(gk + 352);
    const auto *gk_353 = buffer.data(gk + 353);
    const auto *gk_354 = buffer.data(gk + 354);
    const auto *gk_355 = buffer.data(gk + 355);
    const auto *gk_356 = buffer.data(gk + 356);
    const auto *gk_357 = buffer.data(gk + 357);
    const auto *gk_358 = buffer.data(gk + 358);
    const auto *gk_359 = buffer.data(gk + 359);
    const auto *gk_396 = buffer.data(gk + 396);
    const auto *gk_397 = buffer.data(gk + 397);
    const auto *gk_398 = buffer.data(gk + 398);
    const auto *gk_399 = buffer.data(gk + 399);
    const auto *gk_400 = buffer.data(gk + 400);
    const auto *gk_401 = buffer.data(gk + 401);
    const auto *gk_402 = buffer.data(gk + 402);
    const auto *gk_403 = buffer.data(gk + 403);
    const auto *gk_404 = buffer.data(gk + 404);
    const auto *gk_405 = buffer.data(gk + 405);
    const auto *gk_406 = buffer.data(gk + 406);
    const auto *gk_407 = buffer.data(gk + 407);
    const auto *gk_408 = buffer.data(gk + 408);
    const auto *gk_409 = buffer.data(gk + 409);
    const auto *gk_410 = buffer.data(gk + 410);
    const auto *gk_411 = buffer.data(gk + 411);
    const auto *gk_412 = buffer.data(gk + 412);
    const auto *gk_413 = buffer.data(gk + 413);
    const auto *gk_414 = buffer.data(gk + 414);
    const auto *gk_415 = buffer.data(gk + 415);
    const auto *gk_416 = buffer.data(gk + 416);
    const auto *gk_417 = buffer.data(gk + 417);
    const auto *gk_418 = buffer.data(gk + 418);
    const auto *gk_419 = buffer.data(gk + 419);
    const auto *gk_420 = buffer.data(gk + 420);
    const auto *gk_421 = buffer.data(gk + 421);
    const auto *gk_422 = buffer.data(gk + 422);
    const auto *gk_423 = buffer.data(gk + 423);
    const auto *gk_424 = buffer.data(gk + 424);
    const auto *gk_425 = buffer.data(gk + 425);
    const auto *gk_426 = buffer.data(gk + 426);
    const auto *gk_427 = buffer.data(gk + 427);
    const auto *gk_428 = buffer.data(gk + 428);
    const auto *gk_429 = buffer.data(gk + 429);
    const auto *gk_430 = buffer.data(gk + 430);
    const auto *gk_431 = buffer.data(gk + 431);
    const auto *gk_432 = buffer.data(gk + 432);
    const auto *gk_433 = buffer.data(gk + 433);
    const auto *gk_434 = buffer.data(gk + 434);
    const auto *gk_435 = buffer.data(gk + 435);
    const auto *gk_436 = buffer.data(gk + 436);
    const auto *gk_437 = buffer.data(gk + 437);
    const auto *gk_438 = buffer.data(gk + 438);
    const auto *gk_439 = buffer.data(gk + 439);
    const auto *gk_440 = buffer.data(gk + 440);
    const auto *gk_441 = buffer.data(gk + 441);
    const auto *gk_442 = buffer.data(gk + 442);
    const auto *gk_443 = buffer.data(gk + 443);
    const auto *gk_444 = buffer.data(gk + 444);
    const auto *gk_445 = buffer.data(gk + 445);
    const auto *gk_446 = buffer.data(gk + 446);
    const auto *gk_447 = buffer.data(gk + 447);
    const auto *gk_448 = buffer.data(gk + 448);
    const auto *gk_449 = buffer.data(gk + 449);
    const auto *gk_450 = buffer.data(gk + 450);
    const auto *gk_451 = buffer.data(gk + 451);
    const auto *gk_452 = buffer.data(gk + 452);
    const auto *gk_453 = buffer.data(gk + 453);
    const auto *gk_454 = buffer.data(gk + 454);
    const auto *gk_455 = buffer.data(gk + 455);
    const auto *gk_456 = buffer.data(gk + 456);
    const auto *gk_457 = buffer.data(gk + 457);
    const auto *gk_458 = buffer.data(gk + 458);
    const auto *gk_459 = buffer.data(gk + 459);
    const auto *gk_460 = buffer.data(gk + 460);
    const auto *gk_461 = buffer.data(gk + 461);
    const auto *gk_462 = buffer.data(gk + 462);
    const auto *gk_463 = buffer.data(gk + 463);
    const auto *gk_464 = buffer.data(gk + 464);
    const auto *gk_465 = buffer.data(gk + 465);
    const auto *gk_466 = buffer.data(gk + 466);
    const auto *gk_467 = buffer.data(gk + 467);
    const auto *gk_468 = buffer.data(gk + 468);
    const auto *gk_469 = buffer.data(gk + 469);
    const auto *gk_470 = buffer.data(gk + 470);
    const auto *gk_471 = buffer.data(gk + 471);
    const auto *gk_472 = buffer.data(gk + 472);
    const auto *gk_473 = buffer.data(gk + 473);
    const auto *gk_474 = buffer.data(gk + 474);
    const auto *gk_475 = buffer.data(gk + 475);
    const auto *gk_476 = buffer.data(gk + 476);
    const auto *gk_477 = buffer.data(gk + 477);
    const auto *gk_478 = buffer.data(gk + 478);
    const auto *gk_479 = buffer.data(gk + 479);
    const auto *gk_480 = buffer.data(gk + 480);
    const auto *gk_481 = buffer.data(gk + 481);
    const auto *gk_482 = buffer.data(gk + 482);
    const auto *gk_483 = buffer.data(gk + 483);
    const auto *gk_484 = buffer.data(gk + 484);
    const auto *gk_485 = buffer.data(gk + 485);
    const auto *gk_486 = buffer.data(gk + 486);
    const auto *gk_487 = buffer.data(gk + 487);
    const auto *gk_488 = buffer.data(gk + 488);
    const auto *gk_489 = buffer.data(gk + 489);
    const auto *gk_490 = buffer.data(gk + 490);
    const auto *gk_491 = buffer.data(gk + 491);
    const auto *gk_492 = buffer.data(gk + 492);
    const auto *gk_493 = buffer.data(gk + 493);
    const auto *gk_494 = buffer.data(gk + 494);
    const auto *gk_495 = buffer.data(gk + 495);
    const auto *gk_496 = buffer.data(gk + 496);
    const auto *gk_497 = buffer.data(gk + 497);
    const auto *gk_498 = buffer.data(gk + 498);
    const auto *gk_499 = buffer.data(gk + 499);
    const auto *gk_500 = buffer.data(gk + 500);
    const auto *gk_501 = buffer.data(gk + 501);
    const auto *gk_502 = buffer.data(gk + 502);
    const auto *gk_503 = buffer.data(gk + 503);
    const auto *gk_504 = buffer.data(gk + 504);
    const auto *gk_505 = buffer.data(gk + 505);
    const auto *gk_506 = buffer.data(gk + 506);
    const auto *gk_507 = buffer.data(gk + 507);
    const auto *gk_508 = buffer.data(gk + 508);
    const auto *gk_509 = buffer.data(gk + 509);
    const auto *gk_510 = buffer.data(gk + 510);
    const auto *gk_511 = buffer.data(gk + 511);
    const auto *gk_512 = buffer.data(gk + 512);
    const auto *gk_513 = buffer.data(gk + 513);
    const auto *gk_514 = buffer.data(gk + 514);
    const auto *gk_515 = buffer.data(gk + 515);
    const auto *gk_516 = buffer.data(gk + 516);
    const auto *gk_517 = buffer.data(gk + 517);
    const auto *gk_518 = buffer.data(gk + 518);
    const auto *gk_519 = buffer.data(gk + 519);
    const auto *gk_520 = buffer.data(gk + 520);
    const auto *gk_521 = buffer.data(gk + 521);
    const auto *gk_522 = buffer.data(gk + 522);
    const auto *gk_523 = buffer.data(gk + 523);
    const auto *gk_524 = buffer.data(gk + 524);
    const auto *gk_525 = buffer.data(gk + 525);
    const auto *gk_526 = buffer.data(gk + 526);
    const auto *gk_527 = buffer.data(gk + 527);

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, dk_77, dk_78, dk_79, dk_80, dk_81, \
                         gk_329, gk_330, gk_331, gk_332, gk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * dk_77[k]
                   + f_0 * gk_329[k];

        t_186[k] = -2.0 * dk_78[k]
                   + f_0 * gk_330[k];

        t_187[k] = -2.0 * dk_79[k]
                   + f_0 * gk_331[k];

        t_188[k] = -2.0 * dk_80[k]
                   + f_0 * gk_332[k];

        t_189[k] = -2.0 * dk_81[k]
                   + f_0 * gk_333[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, dk_82, dk_83, dk_84, dk_85, dk_86, \
                         gk_334, gk_335, gk_336, gk_337, gk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * dk_82[k]
                   + f_0 * gk_334[k];

        t_191[k] = -2.0 * dk_83[k]
                   + f_0 * gk_335[k];

        t_192[k] = -2.0 * dk_84[k]
                   + f_0 * gk_336[k];

        t_193[k] = -2.0 * dk_85[k]
                   + f_0 * gk_337[k];

        t_194[k] = -2.0 * dk_86[k]
                   + f_0 * gk_338[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, dk_87, dk_88, dk_89, dk_90, dk_91, \
                         gk_339, gk_340, gk_341, gk_342, gk_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * dk_87[k]
                   + f_0 * gk_339[k];

        t_196[k] = -2.0 * dk_88[k]
                   + f_0 * gk_340[k];

        t_197[k] = -2.0 * dk_89[k]
                   + f_0 * gk_341[k];

        t_198[k] = -2.0 * dk_90[k]
                   + f_0 * gk_342[k];

        t_199[k] = -2.0 * dk_91[k]
                   + f_0 * gk_343[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, dk_92, dk_93, dk_94, dk_95, dk_96, \
                         gk_344, gk_345, gk_346, gk_347, gk_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * dk_92[k]
                   + f_0 * gk_344[k];

        t_201[k] = -2.0 * dk_93[k]
                   + f_0 * gk_345[k];

        t_202[k] = -2.0 * dk_94[k]
                   + f_0 * gk_346[k];

        t_203[k] = -2.0 * dk_95[k]
                   + f_0 * gk_347[k];

        t_204[k] = -2.0 * dk_96[k]
                   + f_0 * gk_348[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, dk_97, dk_98, dk_99, dk_100, \
                         dk_101, gk_349, gk_350, gk_351, gk_352, \
                         gk_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * dk_97[k]
                   + f_0 * gk_349[k];

        t_206[k] = -2.0 * dk_98[k]
                   + f_0 * gk_350[k];

        t_207[k] = -2.0 * dk_99[k]
                   + f_0 * gk_351[k];

        t_208[k] = -2.0 * dk_100[k]
                   + f_0 * gk_352[k];

        t_209[k] = -2.0 * dk_101[k]
                   + f_0 * gk_353[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, dk_102, dk_103, dk_104, dk_105, \
                         dk_106, gk_354, gk_355, gk_356, gk_357, \
                         gk_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * dk_102[k]
                   + f_0 * gk_354[k];

        t_211[k] = -2.0 * dk_103[k]
                   + f_0 * gk_355[k];

        t_212[k] = -2.0 * dk_104[k]
                   + f_0 * gk_356[k];

        t_213[k] = -2.0 * dk_105[k]
                   + f_0 * gk_357[k];

        t_214[k] = -2.0 * dk_106[k]
                   + f_0 * gk_358[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, t_220, t_221, dk_107, gk_359, \
                         gk_396, gk_397, gk_398, gk_399, gk_400, \
                         gk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * dk_107[k]
                   + f_0 * gk_359[k];

        t_216[k] = f_0 * gk_396[k];

        t_217[k] = f_0 * gk_397[k];

        t_218[k] = f_0 * gk_398[k];

        t_219[k] = f_0 * gk_399[k];

        t_220[k] = f_0 * gk_400[k];

        t_221[k] = f_0 * gk_401[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, t_228, t_229, gk_402, \
                         gk_403, gk_404, gk_405, gk_406, gk_407, gk_408, \
                         gk_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_0 * gk_402[k];

        t_223[k] = f_0 * gk_403[k];

        t_224[k] = f_0 * gk_404[k];

        t_225[k] = f_0 * gk_405[k];

        t_226[k] = f_0 * gk_406[k];

        t_227[k] = f_0 * gk_407[k];

        t_228[k] = f_0 * gk_408[k];

        t_229[k] = f_0 * gk_409[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, t_235, t_236, t_237, gk_410, \
                         gk_411, gk_412, gk_413, gk_414, gk_415, gk_416, \
                         gk_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_0 * gk_410[k];

        t_231[k] = f_0 * gk_411[k];

        t_232[k] = f_0 * gk_412[k];

        t_233[k] = f_0 * gk_413[k];

        t_234[k] = f_0 * gk_414[k];

        t_235[k] = f_0 * gk_415[k];

        t_236[k] = f_0 * gk_416[k];

        t_237[k] = f_0 * gk_417[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, t_243, t_244, t_245, gk_418, \
                         gk_419, gk_420, gk_421, gk_422, gk_423, gk_424, \
                         gk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_0 * gk_418[k];

        t_239[k] = f_0 * gk_419[k];

        t_240[k] = f_0 * gk_420[k];

        t_241[k] = f_0 * gk_421[k];

        t_242[k] = f_0 * gk_422[k];

        t_243[k] = f_0 * gk_423[k];

        t_244[k] = f_0 * gk_424[k];

        t_245[k] = f_0 * gk_425[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, t_251, t_252, dk_108, gk_426, \
                         gk_427, gk_428, gk_429, gk_430, gk_431, \
                         gk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_0 * gk_426[k];

        t_247[k] = f_0 * gk_427[k];

        t_248[k] = f_0 * gk_428[k];

        t_249[k] = f_0 * gk_429[k];

        t_250[k] = f_0 * gk_430[k];

        t_251[k] = f_0 * gk_431[k];

        t_252[k] = -dk_108[k]
                   + f_0 * gk_432[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, dk_109, dk_110, dk_111, dk_112, \
                         dk_113, gk_433, gk_434, gk_435, gk_436, \
                         gk_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = -dk_109[k]
                   + f_0 * gk_433[k];

        t_254[k] = -dk_110[k]
                   + f_0 * gk_434[k];

        t_255[k] = -dk_111[k]
                   + f_0 * gk_435[k];

        t_256[k] = -dk_112[k]
                   + f_0 * gk_436[k];

        t_257[k] = -dk_113[k]
                   + f_0 * gk_437[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, dk_114, dk_115, dk_116, dk_117, \
                         dk_118, gk_438, gk_439, gk_440, gk_441, \
                         gk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = -dk_114[k]
                   + f_0 * gk_438[k];

        t_259[k] = -dk_115[k]
                   + f_0 * gk_439[k];

        t_260[k] = -dk_116[k]
                   + f_0 * gk_440[k];

        t_261[k] = -dk_117[k]
                   + f_0 * gk_441[k];

        t_262[k] = -dk_118[k]
                   + f_0 * gk_442[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, dk_119, dk_120, dk_121, dk_122, \
                         dk_123, gk_443, gk_444, gk_445, gk_446, \
                         gk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = -dk_119[k]
                   + f_0 * gk_443[k];

        t_264[k] = -dk_120[k]
                   + f_0 * gk_444[k];

        t_265[k] = -dk_121[k]
                   + f_0 * gk_445[k];

        t_266[k] = -dk_122[k]
                   + f_0 * gk_446[k];

        t_267[k] = -dk_123[k]
                   + f_0 * gk_447[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, dk_124, dk_125, dk_126, dk_127, \
                         dk_128, gk_448, gk_449, gk_450, gk_451, \
                         gk_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = -dk_124[k]
                   + f_0 * gk_448[k];

        t_269[k] = -dk_125[k]
                   + f_0 * gk_449[k];

        t_270[k] = -dk_126[k]
                   + f_0 * gk_450[k];

        t_271[k] = -dk_127[k]
                   + f_0 * gk_451[k];

        t_272[k] = -dk_128[k]
                   + f_0 * gk_452[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, dk_129, dk_130, dk_131, dk_132, \
                         dk_133, gk_453, gk_454, gk_455, gk_456, \
                         gk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = -dk_129[k]
                   + f_0 * gk_453[k];

        t_274[k] = -dk_130[k]
                   + f_0 * gk_454[k];

        t_275[k] = -dk_131[k]
                   + f_0 * gk_455[k];

        t_276[k] = -dk_132[k]
                   + f_0 * gk_456[k];

        t_277[k] = -dk_133[k]
                   + f_0 * gk_457[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, dk_134, dk_135, dk_136, dk_137, \
                         dk_138, gk_458, gk_459, gk_460, gk_461, \
                         gk_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = -dk_134[k]
                   + f_0 * gk_458[k];

        t_279[k] = -dk_135[k]
                   + f_0 * gk_459[k];

        t_280[k] = -dk_136[k]
                   + f_0 * gk_460[k];

        t_281[k] = -dk_137[k]
                   + f_0 * gk_461[k];

        t_282[k] = -dk_138[k]
                   + f_0 * gk_462[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, dk_139, dk_140, dk_141, dk_142, \
                         dk_143, gk_463, gk_464, gk_465, gk_466, \
                         gk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = -dk_139[k]
                   + f_0 * gk_463[k];

        t_284[k] = -dk_140[k]
                   + f_0 * gk_464[k];

        t_285[k] = -dk_141[k]
                   + f_0 * gk_465[k];

        t_286[k] = -dk_142[k]
                   + f_0 * gk_466[k];

        t_287[k] = -dk_143[k]
                   + f_0 * gk_467[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, dk_144, dk_145, dk_146, dk_147, \
                         dk_148, gk_468, gk_469, gk_470, gk_471, \
                         gk_472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = -2.0 * dk_144[k]
                   + f_0 * gk_468[k];

        t_289[k] = -2.0 * dk_145[k]
                   + f_0 * gk_469[k];

        t_290[k] = -2.0 * dk_146[k]
                   + f_0 * gk_470[k];

        t_291[k] = -2.0 * dk_147[k]
                   + f_0 * gk_471[k];

        t_292[k] = -2.0 * dk_148[k]
                   + f_0 * gk_472[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, t_297, dk_149, dk_150, dk_151, dk_152, \
                         dk_153, gk_473, gk_474, gk_475, gk_476, \
                         gk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = -2.0 * dk_149[k]
                   + f_0 * gk_473[k];

        t_294[k] = -2.0 * dk_150[k]
                   + f_0 * gk_474[k];

        t_295[k] = -2.0 * dk_151[k]
                   + f_0 * gk_475[k];

        t_296[k] = -2.0 * dk_152[k]
                   + f_0 * gk_476[k];

        t_297[k] = -2.0 * dk_153[k]
                   + f_0 * gk_477[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, dk_154, dk_155, dk_156, dk_157, \
                         dk_158, gk_478, gk_479, gk_480, gk_481, \
                         gk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = -2.0 * dk_154[k]
                   + f_0 * gk_478[k];

        t_299[k] = -2.0 * dk_155[k]
                   + f_0 * gk_479[k];

        t_300[k] = -2.0 * dk_156[k]
                   + f_0 * gk_480[k];

        t_301[k] = -2.0 * dk_157[k]
                   + f_0 * gk_481[k];

        t_302[k] = -2.0 * dk_158[k]
                   + f_0 * gk_482[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, dk_159, dk_160, dk_161, dk_162, \
                         dk_163, gk_483, gk_484, gk_485, gk_486, \
                         gk_487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = -2.0 * dk_159[k]
                   + f_0 * gk_483[k];

        t_304[k] = -2.0 * dk_160[k]
                   + f_0 * gk_484[k];

        t_305[k] = -2.0 * dk_161[k]
                   + f_0 * gk_485[k];

        t_306[k] = -2.0 * dk_162[k]
                   + f_0 * gk_486[k];

        t_307[k] = -2.0 * dk_163[k]
                   + f_0 * gk_487[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, dk_164, dk_165, dk_166, dk_167, \
                         dk_168, gk_488, gk_489, gk_490, gk_491, \
                         gk_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = -2.0 * dk_164[k]
                   + f_0 * gk_488[k];

        t_309[k] = -2.0 * dk_165[k]
                   + f_0 * gk_489[k];

        t_310[k] = -2.0 * dk_166[k]
                   + f_0 * gk_490[k];

        t_311[k] = -2.0 * dk_167[k]
                   + f_0 * gk_491[k];

        t_312[k] = -2.0 * dk_168[k]
                   + f_0 * gk_492[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, dk_169, dk_170, dk_171, dk_172, \
                         dk_173, gk_493, gk_494, gk_495, gk_496, \
                         gk_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = -2.0 * dk_169[k]
                   + f_0 * gk_493[k];

        t_314[k] = -2.0 * dk_170[k]
                   + f_0 * gk_494[k];

        t_315[k] = -2.0 * dk_171[k]
                   + f_0 * gk_495[k];

        t_316[k] = -2.0 * dk_172[k]
                   + f_0 * gk_496[k];

        t_317[k] = -2.0 * dk_173[k]
                   + f_0 * gk_497[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, dk_174, dk_175, dk_176, dk_177, \
                         dk_178, gk_498, gk_499, gk_500, gk_501, \
                         gk_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = -2.0 * dk_174[k]
                   + f_0 * gk_498[k];

        t_319[k] = -2.0 * dk_175[k]
                   + f_0 * gk_499[k];

        t_320[k] = -2.0 * dk_176[k]
                   + f_0 * gk_500[k];

        t_321[k] = -2.0 * dk_177[k]
                   + f_0 * gk_501[k];

        t_322[k] = -2.0 * dk_178[k]
                   + f_0 * gk_502[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, dk_179, dk_180, dk_181, dk_182, \
                         dk_183, gk_503, gk_504, gk_505, gk_506, \
                         gk_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = -2.0 * dk_179[k]
                   + f_0 * gk_503[k];

        t_324[k] = -3.0 * dk_180[k]
                   + f_0 * gk_504[k];

        t_325[k] = -3.0 * dk_181[k]
                   + f_0 * gk_505[k];

        t_326[k] = -3.0 * dk_182[k]
                   + f_0 * gk_506[k];

        t_327[k] = -3.0 * dk_183[k]
                   + f_0 * gk_507[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, dk_184, dk_185, dk_186, dk_187, \
                         dk_188, gk_508, gk_509, gk_510, gk_511, \
                         gk_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = -3.0 * dk_184[k]
                   + f_0 * gk_508[k];

        t_329[k] = -3.0 * dk_185[k]
                   + f_0 * gk_509[k];

        t_330[k] = -3.0 * dk_186[k]
                   + f_0 * gk_510[k];

        t_331[k] = -3.0 * dk_187[k]
                   + f_0 * gk_511[k];

        t_332[k] = -3.0 * dk_188[k]
                   + f_0 * gk_512[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, t_337, dk_189, dk_190, dk_191, dk_192, \
                         dk_193, gk_513, gk_514, gk_515, gk_516, \
                         gk_517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = -3.0 * dk_189[k]
                   + f_0 * gk_513[k];

        t_334[k] = -3.0 * dk_190[k]
                   + f_0 * gk_514[k];

        t_335[k] = -3.0 * dk_191[k]
                   + f_0 * gk_515[k];

        t_336[k] = -3.0 * dk_192[k]
                   + f_0 * gk_516[k];

        t_337[k] = -3.0 * dk_193[k]
                   + f_0 * gk_517[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, dk_194, dk_195, dk_196, dk_197, \
                         dk_198, gk_518, gk_519, gk_520, gk_521, \
                         gk_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = -3.0 * dk_194[k]
                   + f_0 * gk_518[k];

        t_339[k] = -3.0 * dk_195[k]
                   + f_0 * gk_519[k];

        t_340[k] = -3.0 * dk_196[k]
                   + f_0 * gk_520[k];

        t_341[k] = -3.0 * dk_197[k]
                   + f_0 * gk_521[k];

        t_342[k] = -3.0 * dk_198[k]
                   + f_0 * gk_522[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, dk_199, dk_200, dk_201, dk_202, \
                         dk_203, gk_523, gk_524, gk_525, gk_526, \
                         gk_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = -3.0 * dk_199[k]
                   + f_0 * gk_523[k];

        t_344[k] = -3.0 * dk_200[k]
                   + f_0 * gk_524[k];

        t_345[k] = -3.0 * dk_201[k]
                   + f_0 * gk_525[k];

        t_346[k] = -3.0 * dk_202[k]
                   + f_0 * gk_526[k];

        t_347[k] = -3.0 * dk_203[k]
                   + f_0 * gk_527[k];
    }
}

static auto
compute_prim_geom_10_fk_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dk, const size_t gk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *dk_204 = buffer.data(dk + 204);
    const auto *dk_205 = buffer.data(dk + 205);
    const auto *dk_206 = buffer.data(dk + 206);
    const auto *dk_207 = buffer.data(dk + 207);
    const auto *dk_208 = buffer.data(dk + 208);
    const auto *dk_209 = buffer.data(dk + 209);
    const auto *dk_210 = buffer.data(dk + 210);
    const auto *dk_211 = buffer.data(dk + 211);
    const auto *dk_212 = buffer.data(dk + 212);
    const auto *dk_213 = buffer.data(dk + 213);
    const auto *dk_214 = buffer.data(dk + 214);
    const auto *dk_215 = buffer.data(dk + 215);

    const auto *gk_528 = buffer.data(gk + 528);
    const auto *gk_529 = buffer.data(gk + 529);
    const auto *gk_530 = buffer.data(gk + 530);
    const auto *gk_531 = buffer.data(gk + 531);
    const auto *gk_532 = buffer.data(gk + 532);
    const auto *gk_533 = buffer.data(gk + 533);
    const auto *gk_534 = buffer.data(gk + 534);
    const auto *gk_535 = buffer.data(gk + 535);
    const auto *gk_536 = buffer.data(gk + 536);
    const auto *gk_537 = buffer.data(gk + 537);
    const auto *gk_538 = buffer.data(gk + 538);
    const auto *gk_539 = buffer.data(gk + 539);

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, dk_204, dk_205, dk_206, dk_207, \
                         dk_208, gk_528, gk_529, gk_530, gk_531, \
                         gk_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = -3.0 * dk_204[k]
                   + f_0 * gk_528[k];

        t_349[k] = -3.0 * dk_205[k]
                   + f_0 * gk_529[k];

        t_350[k] = -3.0 * dk_206[k]
                   + f_0 * gk_530[k];

        t_351[k] = -3.0 * dk_207[k]
                   + f_0 * gk_531[k];

        t_352[k] = -3.0 * dk_208[k]
                   + f_0 * gk_532[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, dk_209, dk_210, dk_211, dk_212, \
                         dk_213, gk_533, gk_534, gk_535, gk_536, \
                         gk_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = -3.0 * dk_209[k]
                   + f_0 * gk_533[k];

        t_354[k] = -3.0 * dk_210[k]
                   + f_0 * gk_534[k];

        t_355[k] = -3.0 * dk_211[k]
                   + f_0 * gk_535[k];

        t_356[k] = -3.0 * dk_212[k]
                   + f_0 * gk_536[k];

        t_357[k] = -3.0 * dk_213[k]
                   + f_0 * gk_537[k];
    }

#pragma omp simd aligned(t_358, t_359, dk_214, dk_215, gk_538, gk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = -3.0 * dk_214[k]
                   + f_0 * gk_538[k];

        t_359[k] = -3.0 * dk_215[k]
                   + f_0 * gk_539[k];
    }
}

auto
compute_prim_geom_10_fk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t dk, const size_t gk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fk_electron_repulsion_2_piece0(buffer, target, dk, gk, ncols, alpha);

    compute_prim_geom_10_fk_electron_repulsion_2_piece1(buffer, target, dk, gk, ncols, alpha);

    compute_prim_geom_10_fk_electron_repulsion_2_piece2(buffer, target, dk, gk, ncols, alpha);
}

}  // namespace simdt2ceri
