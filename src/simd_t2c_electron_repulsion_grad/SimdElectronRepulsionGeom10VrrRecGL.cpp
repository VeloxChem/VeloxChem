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


#include "SimdElectronRepulsionGeom10VrrRecGL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_gl_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
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

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);
    const auto *fl_9 = buffer.data(fl + 9);
    const auto *fl_10 = buffer.data(fl + 10);
    const auto *fl_11 = buffer.data(fl + 11);
    const auto *fl_12 = buffer.data(fl + 12);
    const auto *fl_13 = buffer.data(fl + 13);
    const auto *fl_14 = buffer.data(fl + 14);
    const auto *fl_15 = buffer.data(fl + 15);
    const auto *fl_16 = buffer.data(fl + 16);
    const auto *fl_17 = buffer.data(fl + 17);
    const auto *fl_18 = buffer.data(fl + 18);
    const auto *fl_19 = buffer.data(fl + 19);
    const auto *fl_20 = buffer.data(fl + 20);
    const auto *fl_21 = buffer.data(fl + 21);
    const auto *fl_22 = buffer.data(fl + 22);
    const auto *fl_23 = buffer.data(fl + 23);
    const auto *fl_24 = buffer.data(fl + 24);
    const auto *fl_25 = buffer.data(fl + 25);
    const auto *fl_26 = buffer.data(fl + 26);
    const auto *fl_27 = buffer.data(fl + 27);
    const auto *fl_28 = buffer.data(fl + 28);
    const auto *fl_29 = buffer.data(fl + 29);
    const auto *fl_30 = buffer.data(fl + 30);
    const auto *fl_31 = buffer.data(fl + 31);
    const auto *fl_32 = buffer.data(fl + 32);
    const auto *fl_33 = buffer.data(fl + 33);
    const auto *fl_34 = buffer.data(fl + 34);
    const auto *fl_35 = buffer.data(fl + 35);
    const auto *fl_36 = buffer.data(fl + 36);
    const auto *fl_37 = buffer.data(fl + 37);
    const auto *fl_38 = buffer.data(fl + 38);
    const auto *fl_39 = buffer.data(fl + 39);
    const auto *fl_40 = buffer.data(fl + 40);
    const auto *fl_41 = buffer.data(fl + 41);
    const auto *fl_42 = buffer.data(fl + 42);
    const auto *fl_43 = buffer.data(fl + 43);
    const auto *fl_44 = buffer.data(fl + 44);
    const auto *fl_45 = buffer.data(fl + 45);
    const auto *fl_46 = buffer.data(fl + 46);
    const auto *fl_47 = buffer.data(fl + 47);
    const auto *fl_48 = buffer.data(fl + 48);
    const auto *fl_49 = buffer.data(fl + 49);
    const auto *fl_50 = buffer.data(fl + 50);
    const auto *fl_51 = buffer.data(fl + 51);
    const auto *fl_52 = buffer.data(fl + 52);
    const auto *fl_53 = buffer.data(fl + 53);
    const auto *fl_54 = buffer.data(fl + 54);
    const auto *fl_55 = buffer.data(fl + 55);
    const auto *fl_56 = buffer.data(fl + 56);
    const auto *fl_57 = buffer.data(fl + 57);
    const auto *fl_58 = buffer.data(fl + 58);
    const auto *fl_59 = buffer.data(fl + 59);
    const auto *fl_60 = buffer.data(fl + 60);
    const auto *fl_61 = buffer.data(fl + 61);
    const auto *fl_62 = buffer.data(fl + 62);
    const auto *fl_63 = buffer.data(fl + 63);
    const auto *fl_64 = buffer.data(fl + 64);
    const auto *fl_65 = buffer.data(fl + 65);
    const auto *fl_66 = buffer.data(fl + 66);
    const auto *fl_67 = buffer.data(fl + 67);
    const auto *fl_68 = buffer.data(fl + 68);
    const auto *fl_69 = buffer.data(fl + 69);
    const auto *fl_70 = buffer.data(fl + 70);
    const auto *fl_71 = buffer.data(fl + 71);
    const auto *fl_72 = buffer.data(fl + 72);
    const auto *fl_73 = buffer.data(fl + 73);
    const auto *fl_74 = buffer.data(fl + 74);
    const auto *fl_75 = buffer.data(fl + 75);
    const auto *fl_76 = buffer.data(fl + 76);
    const auto *fl_77 = buffer.data(fl + 77);
    const auto *fl_78 = buffer.data(fl + 78);
    const auto *fl_79 = buffer.data(fl + 79);
    const auto *fl_80 = buffer.data(fl + 80);
    const auto *fl_81 = buffer.data(fl + 81);
    const auto *fl_82 = buffer.data(fl + 82);
    const auto *fl_83 = buffer.data(fl + 83);
    const auto *fl_84 = buffer.data(fl + 84);
    const auto *fl_85 = buffer.data(fl + 85);
    const auto *fl_86 = buffer.data(fl + 86);
    const auto *fl_87 = buffer.data(fl + 87);
    const auto *fl_88 = buffer.data(fl + 88);
    const auto *fl_89 = buffer.data(fl + 89);
    const auto *fl_90 = buffer.data(fl + 90);
    const auto *fl_91 = buffer.data(fl + 91);
    const auto *fl_92 = buffer.data(fl + 92);
    const auto *fl_93 = buffer.data(fl + 93);
    const auto *fl_94 = buffer.data(fl + 94);
    const auto *fl_95 = buffer.data(fl + 95);
    const auto *fl_96 = buffer.data(fl + 96);
    const auto *fl_97 = buffer.data(fl + 97);
    const auto *fl_98 = buffer.data(fl + 98);
    const auto *fl_99 = buffer.data(fl + 99);
    const auto *fl_100 = buffer.data(fl + 100);
    const auto *fl_101 = buffer.data(fl + 101);
    const auto *fl_102 = buffer.data(fl + 102);
    const auto *fl_103 = buffer.data(fl + 103);
    const auto *fl_104 = buffer.data(fl + 104);
    const auto *fl_105 = buffer.data(fl + 105);
    const auto *fl_106 = buffer.data(fl + 106);
    const auto *fl_107 = buffer.data(fl + 107);
    const auto *fl_108 = buffer.data(fl + 108);
    const auto *fl_109 = buffer.data(fl + 109);
    const auto *fl_110 = buffer.data(fl + 110);
    const auto *fl_111 = buffer.data(fl + 111);
    const auto *fl_112 = buffer.data(fl + 112);
    const auto *fl_113 = buffer.data(fl + 113);
    const auto *fl_114 = buffer.data(fl + 114);
    const auto *fl_115 = buffer.data(fl + 115);
    const auto *fl_116 = buffer.data(fl + 116);
    const auto *fl_117 = buffer.data(fl + 117);
    const auto *fl_118 = buffer.data(fl + 118);
    const auto *fl_119 = buffer.data(fl + 119);
    const auto *fl_120 = buffer.data(fl + 120);
    const auto *fl_121 = buffer.data(fl + 121);
    const auto *fl_122 = buffer.data(fl + 122);
    const auto *fl_123 = buffer.data(fl + 123);
    const auto *fl_124 = buffer.data(fl + 124);
    const auto *fl_125 = buffer.data(fl + 125);
    const auto *fl_126 = buffer.data(fl + 126);
    const auto *fl_127 = buffer.data(fl + 127);
    const auto *fl_128 = buffer.data(fl + 128);
    const auto *fl_129 = buffer.data(fl + 129);
    const auto *fl_130 = buffer.data(fl + 130);
    const auto *fl_131 = buffer.data(fl + 131);
    const auto *fl_132 = buffer.data(fl + 132);
    const auto *fl_133 = buffer.data(fl + 133);
    const auto *fl_134 = buffer.data(fl + 134);
    const auto *fl_135 = buffer.data(fl + 135);
    const auto *fl_136 = buffer.data(fl + 136);
    const auto *fl_137 = buffer.data(fl + 137);
    const auto *fl_138 = buffer.data(fl + 138);
    const auto *fl_139 = buffer.data(fl + 139);
    const auto *fl_140 = buffer.data(fl + 140);
    const auto *fl_141 = buffer.data(fl + 141);
    const auto *fl_142 = buffer.data(fl + 142);
    const auto *fl_143 = buffer.data(fl + 143);
    const auto *fl_144 = buffer.data(fl + 144);
    const auto *fl_145 = buffer.data(fl + 145);
    const auto *fl_146 = buffer.data(fl + 146);
    const auto *fl_147 = buffer.data(fl + 147);
    const auto *fl_148 = buffer.data(fl + 148);
    const auto *fl_149 = buffer.data(fl + 149);

    const auto *hl_0 = buffer.data(hl + 0);
    const auto *hl_1 = buffer.data(hl + 1);
    const auto *hl_2 = buffer.data(hl + 2);
    const auto *hl_3 = buffer.data(hl + 3);
    const auto *hl_4 = buffer.data(hl + 4);
    const auto *hl_5 = buffer.data(hl + 5);
    const auto *hl_6 = buffer.data(hl + 6);
    const auto *hl_7 = buffer.data(hl + 7);
    const auto *hl_8 = buffer.data(hl + 8);
    const auto *hl_9 = buffer.data(hl + 9);
    const auto *hl_10 = buffer.data(hl + 10);
    const auto *hl_11 = buffer.data(hl + 11);
    const auto *hl_12 = buffer.data(hl + 12);
    const auto *hl_13 = buffer.data(hl + 13);
    const auto *hl_14 = buffer.data(hl + 14);
    const auto *hl_15 = buffer.data(hl + 15);
    const auto *hl_16 = buffer.data(hl + 16);
    const auto *hl_17 = buffer.data(hl + 17);
    const auto *hl_18 = buffer.data(hl + 18);
    const auto *hl_19 = buffer.data(hl + 19);
    const auto *hl_20 = buffer.data(hl + 20);
    const auto *hl_21 = buffer.data(hl + 21);
    const auto *hl_22 = buffer.data(hl + 22);
    const auto *hl_23 = buffer.data(hl + 23);
    const auto *hl_24 = buffer.data(hl + 24);
    const auto *hl_25 = buffer.data(hl + 25);
    const auto *hl_26 = buffer.data(hl + 26);
    const auto *hl_27 = buffer.data(hl + 27);
    const auto *hl_28 = buffer.data(hl + 28);
    const auto *hl_29 = buffer.data(hl + 29);
    const auto *hl_30 = buffer.data(hl + 30);
    const auto *hl_31 = buffer.data(hl + 31);
    const auto *hl_32 = buffer.data(hl + 32);
    const auto *hl_33 = buffer.data(hl + 33);
    const auto *hl_34 = buffer.data(hl + 34);
    const auto *hl_35 = buffer.data(hl + 35);
    const auto *hl_36 = buffer.data(hl + 36);
    const auto *hl_37 = buffer.data(hl + 37);
    const auto *hl_38 = buffer.data(hl + 38);
    const auto *hl_39 = buffer.data(hl + 39);
    const auto *hl_40 = buffer.data(hl + 40);
    const auto *hl_41 = buffer.data(hl + 41);
    const auto *hl_42 = buffer.data(hl + 42);
    const auto *hl_43 = buffer.data(hl + 43);
    const auto *hl_44 = buffer.data(hl + 44);
    const auto *hl_45 = buffer.data(hl + 45);
    const auto *hl_46 = buffer.data(hl + 46);
    const auto *hl_47 = buffer.data(hl + 47);
    const auto *hl_48 = buffer.data(hl + 48);
    const auto *hl_49 = buffer.data(hl + 49);
    const auto *hl_50 = buffer.data(hl + 50);
    const auto *hl_51 = buffer.data(hl + 51);
    const auto *hl_52 = buffer.data(hl + 52);
    const auto *hl_53 = buffer.data(hl + 53);
    const auto *hl_54 = buffer.data(hl + 54);
    const auto *hl_55 = buffer.data(hl + 55);
    const auto *hl_56 = buffer.data(hl + 56);
    const auto *hl_57 = buffer.data(hl + 57);
    const auto *hl_58 = buffer.data(hl + 58);
    const auto *hl_59 = buffer.data(hl + 59);
    const auto *hl_60 = buffer.data(hl + 60);
    const auto *hl_61 = buffer.data(hl + 61);
    const auto *hl_62 = buffer.data(hl + 62);
    const auto *hl_63 = buffer.data(hl + 63);
    const auto *hl_64 = buffer.data(hl + 64);
    const auto *hl_65 = buffer.data(hl + 65);
    const auto *hl_66 = buffer.data(hl + 66);
    const auto *hl_67 = buffer.data(hl + 67);
    const auto *hl_68 = buffer.data(hl + 68);
    const auto *hl_69 = buffer.data(hl + 69);
    const auto *hl_70 = buffer.data(hl + 70);
    const auto *hl_71 = buffer.data(hl + 71);
    const auto *hl_72 = buffer.data(hl + 72);
    const auto *hl_73 = buffer.data(hl + 73);
    const auto *hl_74 = buffer.data(hl + 74);
    const auto *hl_75 = buffer.data(hl + 75);
    const auto *hl_76 = buffer.data(hl + 76);
    const auto *hl_77 = buffer.data(hl + 77);
    const auto *hl_78 = buffer.data(hl + 78);
    const auto *hl_79 = buffer.data(hl + 79);
    const auto *hl_80 = buffer.data(hl + 80);
    const auto *hl_81 = buffer.data(hl + 81);
    const auto *hl_82 = buffer.data(hl + 82);
    const auto *hl_83 = buffer.data(hl + 83);
    const auto *hl_84 = buffer.data(hl + 84);
    const auto *hl_85 = buffer.data(hl + 85);
    const auto *hl_86 = buffer.data(hl + 86);
    const auto *hl_87 = buffer.data(hl + 87);
    const auto *hl_88 = buffer.data(hl + 88);
    const auto *hl_89 = buffer.data(hl + 89);
    const auto *hl_90 = buffer.data(hl + 90);
    const auto *hl_91 = buffer.data(hl + 91);
    const auto *hl_92 = buffer.data(hl + 92);
    const auto *hl_93 = buffer.data(hl + 93);
    const auto *hl_94 = buffer.data(hl + 94);
    const auto *hl_95 = buffer.data(hl + 95);
    const auto *hl_96 = buffer.data(hl + 96);
    const auto *hl_97 = buffer.data(hl + 97);
    const auto *hl_98 = buffer.data(hl + 98);
    const auto *hl_99 = buffer.data(hl + 99);
    const auto *hl_100 = buffer.data(hl + 100);
    const auto *hl_101 = buffer.data(hl + 101);
    const auto *hl_102 = buffer.data(hl + 102);
    const auto *hl_103 = buffer.data(hl + 103);
    const auto *hl_104 = buffer.data(hl + 104);
    const auto *hl_105 = buffer.data(hl + 105);
    const auto *hl_106 = buffer.data(hl + 106);
    const auto *hl_107 = buffer.data(hl + 107);
    const auto *hl_108 = buffer.data(hl + 108);
    const auto *hl_109 = buffer.data(hl + 109);
    const auto *hl_110 = buffer.data(hl + 110);
    const auto *hl_111 = buffer.data(hl + 111);
    const auto *hl_112 = buffer.data(hl + 112);
    const auto *hl_113 = buffer.data(hl + 113);
    const auto *hl_114 = buffer.data(hl + 114);
    const auto *hl_115 = buffer.data(hl + 115);
    const auto *hl_116 = buffer.data(hl + 116);
    const auto *hl_117 = buffer.data(hl + 117);
    const auto *hl_118 = buffer.data(hl + 118);
    const auto *hl_119 = buffer.data(hl + 119);
    const auto *hl_120 = buffer.data(hl + 120);
    const auto *hl_121 = buffer.data(hl + 121);
    const auto *hl_122 = buffer.data(hl + 122);
    const auto *hl_123 = buffer.data(hl + 123);
    const auto *hl_124 = buffer.data(hl + 124);
    const auto *hl_125 = buffer.data(hl + 125);
    const auto *hl_126 = buffer.data(hl + 126);
    const auto *hl_127 = buffer.data(hl + 127);
    const auto *hl_128 = buffer.data(hl + 128);
    const auto *hl_129 = buffer.data(hl + 129);
    const auto *hl_130 = buffer.data(hl + 130);
    const auto *hl_131 = buffer.data(hl + 131);
    const auto *hl_132 = buffer.data(hl + 132);
    const auto *hl_133 = buffer.data(hl + 133);
    const auto *hl_134 = buffer.data(hl + 134);
    const auto *hl_135 = buffer.data(hl + 135);
    const auto *hl_136 = buffer.data(hl + 136);
    const auto *hl_137 = buffer.data(hl + 137);
    const auto *hl_138 = buffer.data(hl + 138);
    const auto *hl_139 = buffer.data(hl + 139);
    const auto *hl_140 = buffer.data(hl + 140);
    const auto *hl_141 = buffer.data(hl + 141);
    const auto *hl_142 = buffer.data(hl + 142);
    const auto *hl_143 = buffer.data(hl + 143);
    const auto *hl_144 = buffer.data(hl + 144);
    const auto *hl_145 = buffer.data(hl + 145);
    const auto *hl_146 = buffer.data(hl + 146);
    const auto *hl_147 = buffer.data(hl + 147);
    const auto *hl_148 = buffer.data(hl + 148);
    const auto *hl_149 = buffer.data(hl + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, fl_0, fl_1, fl_2, fl_3, fl_4, hl_0, hl_1, \
                         hl_2, hl_3, hl_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -4.0 * fl_0[k]
                 + f_0 * hl_0[k];

        t_1[k] = -4.0 * fl_1[k]
                 + f_0 * hl_1[k];

        t_2[k] = -4.0 * fl_2[k]
                 + f_0 * hl_2[k];

        t_3[k] = -4.0 * fl_3[k]
                 + f_0 * hl_3[k];

        t_4[k] = -4.0 * fl_4[k]
                 + f_0 * hl_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, fl_5, fl_6, fl_7, fl_8, fl_9, hl_5, hl_6, \
                         hl_7, hl_8, hl_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -4.0 * fl_5[k]
                 + f_0 * hl_5[k];

        t_6[k] = -4.0 * fl_6[k]
                 + f_0 * hl_6[k];

        t_7[k] = -4.0 * fl_7[k]
                 + f_0 * hl_7[k];

        t_8[k] = -4.0 * fl_8[k]
                 + f_0 * hl_8[k];

        t_9[k] = -4.0 * fl_9[k]
                 + f_0 * hl_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, fl_10, fl_11, fl_12, fl_13, fl_14, \
                         hl_10, hl_11, hl_12, hl_13, hl_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -4.0 * fl_10[k]
                  + f_0 * hl_10[k];

        t_11[k] = -4.0 * fl_11[k]
                  + f_0 * hl_11[k];

        t_12[k] = -4.0 * fl_12[k]
                  + f_0 * hl_12[k];

        t_13[k] = -4.0 * fl_13[k]
                  + f_0 * hl_13[k];

        t_14[k] = -4.0 * fl_14[k]
                  + f_0 * hl_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, fl_15, fl_16, fl_17, fl_18, fl_19, \
                         hl_15, hl_16, hl_17, hl_18, hl_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -4.0 * fl_15[k]
                  + f_0 * hl_15[k];

        t_16[k] = -4.0 * fl_16[k]
                  + f_0 * hl_16[k];

        t_17[k] = -4.0 * fl_17[k]
                  + f_0 * hl_17[k];

        t_18[k] = -4.0 * fl_18[k]
                  + f_0 * hl_18[k];

        t_19[k] = -4.0 * fl_19[k]
                  + f_0 * hl_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, fl_20, fl_21, fl_22, fl_23, fl_24, \
                         hl_20, hl_21, hl_22, hl_23, hl_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -4.0 * fl_20[k]
                  + f_0 * hl_20[k];

        t_21[k] = -4.0 * fl_21[k]
                  + f_0 * hl_21[k];

        t_22[k] = -4.0 * fl_22[k]
                  + f_0 * hl_22[k];

        t_23[k] = -4.0 * fl_23[k]
                  + f_0 * hl_23[k];

        t_24[k] = -4.0 * fl_24[k]
                  + f_0 * hl_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, fl_25, fl_26, fl_27, fl_28, fl_29, \
                         hl_25, hl_26, hl_27, hl_28, hl_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -4.0 * fl_25[k]
                  + f_0 * hl_25[k];

        t_26[k] = -4.0 * fl_26[k]
                  + f_0 * hl_26[k];

        t_27[k] = -4.0 * fl_27[k]
                  + f_0 * hl_27[k];

        t_28[k] = -4.0 * fl_28[k]
                  + f_0 * hl_28[k];

        t_29[k] = -4.0 * fl_29[k]
                  + f_0 * hl_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, fl_30, fl_31, fl_32, fl_33, fl_34, \
                         hl_30, hl_31, hl_32, hl_33, hl_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -4.0 * fl_30[k]
                  + f_0 * hl_30[k];

        t_31[k] = -4.0 * fl_31[k]
                  + f_0 * hl_31[k];

        t_32[k] = -4.0 * fl_32[k]
                  + f_0 * hl_32[k];

        t_33[k] = -4.0 * fl_33[k]
                  + f_0 * hl_33[k];

        t_34[k] = -4.0 * fl_34[k]
                  + f_0 * hl_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, fl_35, fl_36, fl_37, fl_38, fl_39, \
                         hl_35, hl_36, hl_37, hl_38, hl_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -4.0 * fl_35[k]
                  + f_0 * hl_35[k];

        t_36[k] = -4.0 * fl_36[k]
                  + f_0 * hl_36[k];

        t_37[k] = -4.0 * fl_37[k]
                  + f_0 * hl_37[k];

        t_38[k] = -4.0 * fl_38[k]
                  + f_0 * hl_38[k];

        t_39[k] = -4.0 * fl_39[k]
                  + f_0 * hl_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, fl_40, fl_41, fl_42, fl_43, fl_44, \
                         hl_40, hl_41, hl_42, hl_43, hl_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -4.0 * fl_40[k]
                  + f_0 * hl_40[k];

        t_41[k] = -4.0 * fl_41[k]
                  + f_0 * hl_41[k];

        t_42[k] = -4.0 * fl_42[k]
                  + f_0 * hl_42[k];

        t_43[k] = -4.0 * fl_43[k]
                  + f_0 * hl_43[k];

        t_44[k] = -4.0 * fl_44[k]
                  + f_0 * hl_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, fl_45, fl_46, fl_47, fl_48, fl_49, \
                         hl_45, hl_46, hl_47, hl_48, hl_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -3.0 * fl_45[k]
                  + f_0 * hl_45[k];

        t_46[k] = -3.0 * fl_46[k]
                  + f_0 * hl_46[k];

        t_47[k] = -3.0 * fl_47[k]
                  + f_0 * hl_47[k];

        t_48[k] = -3.0 * fl_48[k]
                  + f_0 * hl_48[k];

        t_49[k] = -3.0 * fl_49[k]
                  + f_0 * hl_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, fl_50, fl_51, fl_52, fl_53, fl_54, \
                         hl_50, hl_51, hl_52, hl_53, hl_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -3.0 * fl_50[k]
                  + f_0 * hl_50[k];

        t_51[k] = -3.0 * fl_51[k]
                  + f_0 * hl_51[k];

        t_52[k] = -3.0 * fl_52[k]
                  + f_0 * hl_52[k];

        t_53[k] = -3.0 * fl_53[k]
                  + f_0 * hl_53[k];

        t_54[k] = -3.0 * fl_54[k]
                  + f_0 * hl_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, fl_55, fl_56, fl_57, fl_58, fl_59, \
                         hl_55, hl_56, hl_57, hl_58, hl_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -3.0 * fl_55[k]
                  + f_0 * hl_55[k];

        t_56[k] = -3.0 * fl_56[k]
                  + f_0 * hl_56[k];

        t_57[k] = -3.0 * fl_57[k]
                  + f_0 * hl_57[k];

        t_58[k] = -3.0 * fl_58[k]
                  + f_0 * hl_58[k];

        t_59[k] = -3.0 * fl_59[k]
                  + f_0 * hl_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, fl_60, fl_61, fl_62, fl_63, fl_64, \
                         hl_60, hl_61, hl_62, hl_63, hl_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * fl_60[k]
                  + f_0 * hl_60[k];

        t_61[k] = -3.0 * fl_61[k]
                  + f_0 * hl_61[k];

        t_62[k] = -3.0 * fl_62[k]
                  + f_0 * hl_62[k];

        t_63[k] = -3.0 * fl_63[k]
                  + f_0 * hl_63[k];

        t_64[k] = -3.0 * fl_64[k]
                  + f_0 * hl_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, fl_65, fl_66, fl_67, fl_68, fl_69, \
                         hl_65, hl_66, hl_67, hl_68, hl_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * fl_65[k]
                  + f_0 * hl_65[k];

        t_66[k] = -3.0 * fl_66[k]
                  + f_0 * hl_66[k];

        t_67[k] = -3.0 * fl_67[k]
                  + f_0 * hl_67[k];

        t_68[k] = -3.0 * fl_68[k]
                  + f_0 * hl_68[k];

        t_69[k] = -3.0 * fl_69[k]
                  + f_0 * hl_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, fl_70, fl_71, fl_72, fl_73, fl_74, \
                         hl_70, hl_71, hl_72, hl_73, hl_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -3.0 * fl_70[k]
                  + f_0 * hl_70[k];

        t_71[k] = -3.0 * fl_71[k]
                  + f_0 * hl_71[k];

        t_72[k] = -3.0 * fl_72[k]
                  + f_0 * hl_72[k];

        t_73[k] = -3.0 * fl_73[k]
                  + f_0 * hl_73[k];

        t_74[k] = -3.0 * fl_74[k]
                  + f_0 * hl_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, fl_75, fl_76, fl_77, fl_78, fl_79, \
                         hl_75, hl_76, hl_77, hl_78, hl_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -3.0 * fl_75[k]
                  + f_0 * hl_75[k];

        t_76[k] = -3.0 * fl_76[k]
                  + f_0 * hl_76[k];

        t_77[k] = -3.0 * fl_77[k]
                  + f_0 * hl_77[k];

        t_78[k] = -3.0 * fl_78[k]
                  + f_0 * hl_78[k];

        t_79[k] = -3.0 * fl_79[k]
                  + f_0 * hl_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, fl_80, fl_81, fl_82, fl_83, fl_84, \
                         hl_80, hl_81, hl_82, hl_83, hl_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * fl_80[k]
                  + f_0 * hl_80[k];

        t_81[k] = -3.0 * fl_81[k]
                  + f_0 * hl_81[k];

        t_82[k] = -3.0 * fl_82[k]
                  + f_0 * hl_82[k];

        t_83[k] = -3.0 * fl_83[k]
                  + f_0 * hl_83[k];

        t_84[k] = -3.0 * fl_84[k]
                  + f_0 * hl_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, fl_85, fl_86, fl_87, fl_88, fl_89, \
                         hl_85, hl_86, hl_87, hl_88, hl_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -3.0 * fl_85[k]
                  + f_0 * hl_85[k];

        t_86[k] = -3.0 * fl_86[k]
                  + f_0 * hl_86[k];

        t_87[k] = -3.0 * fl_87[k]
                  + f_0 * hl_87[k];

        t_88[k] = -3.0 * fl_88[k]
                  + f_0 * hl_88[k];

        t_89[k] = -3.0 * fl_89[k]
                  + f_0 * hl_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, fl_90, fl_91, fl_92, fl_93, fl_94, \
                         hl_90, hl_91, hl_92, hl_93, hl_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * fl_90[k]
                  + f_0 * hl_90[k];

        t_91[k] = -3.0 * fl_91[k]
                  + f_0 * hl_91[k];

        t_92[k] = -3.0 * fl_92[k]
                  + f_0 * hl_92[k];

        t_93[k] = -3.0 * fl_93[k]
                  + f_0 * hl_93[k];

        t_94[k] = -3.0 * fl_94[k]
                  + f_0 * hl_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, fl_95, fl_96, fl_97, fl_98, fl_99, \
                         hl_95, hl_96, hl_97, hl_98, hl_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * fl_95[k]
                  + f_0 * hl_95[k];

        t_96[k] = -3.0 * fl_96[k]
                  + f_0 * hl_96[k];

        t_97[k] = -3.0 * fl_97[k]
                  + f_0 * hl_97[k];

        t_98[k] = -3.0 * fl_98[k]
                  + f_0 * hl_98[k];

        t_99[k] = -3.0 * fl_99[k]
                  + f_0 * hl_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, fl_100, fl_101, fl_102, fl_103, \
                         fl_104, hl_100, hl_101, hl_102, hl_103, \
                         hl_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * fl_100[k]
                   + f_0 * hl_100[k];

        t_101[k] = -3.0 * fl_101[k]
                   + f_0 * hl_101[k];

        t_102[k] = -3.0 * fl_102[k]
                   + f_0 * hl_102[k];

        t_103[k] = -3.0 * fl_103[k]
                   + f_0 * hl_103[k];

        t_104[k] = -3.0 * fl_104[k]
                   + f_0 * hl_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, fl_105, fl_106, fl_107, fl_108, \
                         fl_109, hl_105, hl_106, hl_107, hl_108, \
                         hl_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -3.0 * fl_105[k]
                   + f_0 * hl_105[k];

        t_106[k] = -3.0 * fl_106[k]
                   + f_0 * hl_106[k];

        t_107[k] = -3.0 * fl_107[k]
                   + f_0 * hl_107[k];

        t_108[k] = -3.0 * fl_108[k]
                   + f_0 * hl_108[k];

        t_109[k] = -3.0 * fl_109[k]
                   + f_0 * hl_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, fl_110, fl_111, fl_112, fl_113, \
                         fl_114, hl_110, hl_111, hl_112, hl_113, \
                         hl_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -3.0 * fl_110[k]
                   + f_0 * hl_110[k];

        t_111[k] = -3.0 * fl_111[k]
                   + f_0 * hl_111[k];

        t_112[k] = -3.0 * fl_112[k]
                   + f_0 * hl_112[k];

        t_113[k] = -3.0 * fl_113[k]
                   + f_0 * hl_113[k];

        t_114[k] = -3.0 * fl_114[k]
                   + f_0 * hl_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, fl_115, fl_116, fl_117, fl_118, \
                         fl_119, hl_115, hl_116, hl_117, hl_118, \
                         hl_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -3.0 * fl_115[k]
                   + f_0 * hl_115[k];

        t_116[k] = -3.0 * fl_116[k]
                   + f_0 * hl_116[k];

        t_117[k] = -3.0 * fl_117[k]
                   + f_0 * hl_117[k];

        t_118[k] = -3.0 * fl_118[k]
                   + f_0 * hl_118[k];

        t_119[k] = -3.0 * fl_119[k]
                   + f_0 * hl_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, fl_120, fl_121, fl_122, fl_123, \
                         fl_124, hl_120, hl_121, hl_122, hl_123, \
                         hl_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -3.0 * fl_120[k]
                   + f_0 * hl_120[k];

        t_121[k] = -3.0 * fl_121[k]
                   + f_0 * hl_121[k];

        t_122[k] = -3.0 * fl_122[k]
                   + f_0 * hl_122[k];

        t_123[k] = -3.0 * fl_123[k]
                   + f_0 * hl_123[k];

        t_124[k] = -3.0 * fl_124[k]
                   + f_0 * hl_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, fl_125, fl_126, fl_127, fl_128, \
                         fl_129, hl_125, hl_126, hl_127, hl_128, \
                         hl_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -3.0 * fl_125[k]
                   + f_0 * hl_125[k];

        t_126[k] = -3.0 * fl_126[k]
                   + f_0 * hl_126[k];

        t_127[k] = -3.0 * fl_127[k]
                   + f_0 * hl_127[k];

        t_128[k] = -3.0 * fl_128[k]
                   + f_0 * hl_128[k];

        t_129[k] = -3.0 * fl_129[k]
                   + f_0 * hl_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, fl_130, fl_131, fl_132, fl_133, \
                         fl_134, hl_130, hl_131, hl_132, hl_133, \
                         hl_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -3.0 * fl_130[k]
                   + f_0 * hl_130[k];

        t_131[k] = -3.0 * fl_131[k]
                   + f_0 * hl_131[k];

        t_132[k] = -3.0 * fl_132[k]
                   + f_0 * hl_132[k];

        t_133[k] = -3.0 * fl_133[k]
                   + f_0 * hl_133[k];

        t_134[k] = -3.0 * fl_134[k]
                   + f_0 * hl_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, fl_135, fl_136, fl_137, fl_138, \
                         fl_139, hl_135, hl_136, hl_137, hl_138, \
                         hl_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -2.0 * fl_135[k]
                   + f_0 * hl_135[k];

        t_136[k] = -2.0 * fl_136[k]
                   + f_0 * hl_136[k];

        t_137[k] = -2.0 * fl_137[k]
                   + f_0 * hl_137[k];

        t_138[k] = -2.0 * fl_138[k]
                   + f_0 * hl_138[k];

        t_139[k] = -2.0 * fl_139[k]
                   + f_0 * hl_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, fl_140, fl_141, fl_142, fl_143, \
                         fl_144, hl_140, hl_141, hl_142, hl_143, \
                         hl_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -2.0 * fl_140[k]
                   + f_0 * hl_140[k];

        t_141[k] = -2.0 * fl_141[k]
                   + f_0 * hl_141[k];

        t_142[k] = -2.0 * fl_142[k]
                   + f_0 * hl_142[k];

        t_143[k] = -2.0 * fl_143[k]
                   + f_0 * hl_143[k];

        t_144[k] = -2.0 * fl_144[k]
                   + f_0 * hl_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, fl_145, fl_146, fl_147, fl_148, \
                         fl_149, hl_145, hl_146, hl_147, hl_148, \
                         hl_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -2.0 * fl_145[k]
                   + f_0 * hl_145[k];

        t_146[k] = -2.0 * fl_146[k]
                   + f_0 * hl_146[k];

        t_147[k] = -2.0 * fl_147[k]
                   + f_0 * hl_147[k];

        t_148[k] = -2.0 * fl_148[k]
                   + f_0 * hl_148[k];

        t_149[k] = -2.0 * fl_149[k]
                   + f_0 * hl_149[k];
    }
}

static auto
compute_prim_geom_10_gl_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
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

    const auto *fl_150 = buffer.data(fl + 150);
    const auto *fl_151 = buffer.data(fl + 151);
    const auto *fl_152 = buffer.data(fl + 152);
    const auto *fl_153 = buffer.data(fl + 153);
    const auto *fl_154 = buffer.data(fl + 154);
    const auto *fl_155 = buffer.data(fl + 155);
    const auto *fl_156 = buffer.data(fl + 156);
    const auto *fl_157 = buffer.data(fl + 157);
    const auto *fl_158 = buffer.data(fl + 158);
    const auto *fl_159 = buffer.data(fl + 159);
    const auto *fl_160 = buffer.data(fl + 160);
    const auto *fl_161 = buffer.data(fl + 161);
    const auto *fl_162 = buffer.data(fl + 162);
    const auto *fl_163 = buffer.data(fl + 163);
    const auto *fl_164 = buffer.data(fl + 164);
    const auto *fl_165 = buffer.data(fl + 165);
    const auto *fl_166 = buffer.data(fl + 166);
    const auto *fl_167 = buffer.data(fl + 167);
    const auto *fl_168 = buffer.data(fl + 168);
    const auto *fl_169 = buffer.data(fl + 169);
    const auto *fl_170 = buffer.data(fl + 170);
    const auto *fl_171 = buffer.data(fl + 171);
    const auto *fl_172 = buffer.data(fl + 172);
    const auto *fl_173 = buffer.data(fl + 173);
    const auto *fl_174 = buffer.data(fl + 174);
    const auto *fl_175 = buffer.data(fl + 175);
    const auto *fl_176 = buffer.data(fl + 176);
    const auto *fl_177 = buffer.data(fl + 177);
    const auto *fl_178 = buffer.data(fl + 178);
    const auto *fl_179 = buffer.data(fl + 179);
    const auto *fl_180 = buffer.data(fl + 180);
    const auto *fl_181 = buffer.data(fl + 181);
    const auto *fl_182 = buffer.data(fl + 182);
    const auto *fl_183 = buffer.data(fl + 183);
    const auto *fl_184 = buffer.data(fl + 184);
    const auto *fl_185 = buffer.data(fl + 185);
    const auto *fl_186 = buffer.data(fl + 186);
    const auto *fl_187 = buffer.data(fl + 187);
    const auto *fl_188 = buffer.data(fl + 188);
    const auto *fl_189 = buffer.data(fl + 189);
    const auto *fl_190 = buffer.data(fl + 190);
    const auto *fl_191 = buffer.data(fl + 191);
    const auto *fl_192 = buffer.data(fl + 192);
    const auto *fl_193 = buffer.data(fl + 193);
    const auto *fl_194 = buffer.data(fl + 194);
    const auto *fl_195 = buffer.data(fl + 195);
    const auto *fl_196 = buffer.data(fl + 196);
    const auto *fl_197 = buffer.data(fl + 197);
    const auto *fl_198 = buffer.data(fl + 198);
    const auto *fl_199 = buffer.data(fl + 199);
    const auto *fl_200 = buffer.data(fl + 200);
    const auto *fl_201 = buffer.data(fl + 201);
    const auto *fl_202 = buffer.data(fl + 202);
    const auto *fl_203 = buffer.data(fl + 203);
    const auto *fl_204 = buffer.data(fl + 204);
    const auto *fl_205 = buffer.data(fl + 205);
    const auto *fl_206 = buffer.data(fl + 206);
    const auto *fl_207 = buffer.data(fl + 207);
    const auto *fl_208 = buffer.data(fl + 208);
    const auto *fl_209 = buffer.data(fl + 209);
    const auto *fl_210 = buffer.data(fl + 210);
    const auto *fl_211 = buffer.data(fl + 211);
    const auto *fl_212 = buffer.data(fl + 212);
    const auto *fl_213 = buffer.data(fl + 213);
    const auto *fl_214 = buffer.data(fl + 214);
    const auto *fl_215 = buffer.data(fl + 215);
    const auto *fl_216 = buffer.data(fl + 216);
    const auto *fl_217 = buffer.data(fl + 217);
    const auto *fl_218 = buffer.data(fl + 218);
    const auto *fl_219 = buffer.data(fl + 219);
    const auto *fl_220 = buffer.data(fl + 220);
    const auto *fl_221 = buffer.data(fl + 221);
    const auto *fl_222 = buffer.data(fl + 222);
    const auto *fl_223 = buffer.data(fl + 223);
    const auto *fl_224 = buffer.data(fl + 224);
    const auto *fl_225 = buffer.data(fl + 225);
    const auto *fl_226 = buffer.data(fl + 226);
    const auto *fl_227 = buffer.data(fl + 227);
    const auto *fl_228 = buffer.data(fl + 228);
    const auto *fl_229 = buffer.data(fl + 229);
    const auto *fl_230 = buffer.data(fl + 230);
    const auto *fl_231 = buffer.data(fl + 231);
    const auto *fl_232 = buffer.data(fl + 232);
    const auto *fl_233 = buffer.data(fl + 233);
    const auto *fl_234 = buffer.data(fl + 234);
    const auto *fl_235 = buffer.data(fl + 235);
    const auto *fl_236 = buffer.data(fl + 236);
    const auto *fl_237 = buffer.data(fl + 237);
    const auto *fl_238 = buffer.data(fl + 238);
    const auto *fl_239 = buffer.data(fl + 239);
    const auto *fl_240 = buffer.data(fl + 240);
    const auto *fl_241 = buffer.data(fl + 241);
    const auto *fl_242 = buffer.data(fl + 242);
    const auto *fl_243 = buffer.data(fl + 243);
    const auto *fl_244 = buffer.data(fl + 244);
    const auto *fl_245 = buffer.data(fl + 245);
    const auto *fl_246 = buffer.data(fl + 246);
    const auto *fl_247 = buffer.data(fl + 247);
    const auto *fl_248 = buffer.data(fl + 248);
    const auto *fl_249 = buffer.data(fl + 249);
    const auto *fl_250 = buffer.data(fl + 250);
    const auto *fl_251 = buffer.data(fl + 251);
    const auto *fl_252 = buffer.data(fl + 252);
    const auto *fl_253 = buffer.data(fl + 253);
    const auto *fl_254 = buffer.data(fl + 254);
    const auto *fl_255 = buffer.data(fl + 255);
    const auto *fl_256 = buffer.data(fl + 256);
    const auto *fl_257 = buffer.data(fl + 257);
    const auto *fl_258 = buffer.data(fl + 258);
    const auto *fl_259 = buffer.data(fl + 259);
    const auto *fl_260 = buffer.data(fl + 260);
    const auto *fl_261 = buffer.data(fl + 261);
    const auto *fl_262 = buffer.data(fl + 262);
    const auto *fl_263 = buffer.data(fl + 263);
    const auto *fl_264 = buffer.data(fl + 264);
    const auto *fl_265 = buffer.data(fl + 265);
    const auto *fl_266 = buffer.data(fl + 266);
    const auto *fl_267 = buffer.data(fl + 267);
    const auto *fl_268 = buffer.data(fl + 268);
    const auto *fl_269 = buffer.data(fl + 269);
    const auto *fl_270 = buffer.data(fl + 270);
    const auto *fl_271 = buffer.data(fl + 271);
    const auto *fl_272 = buffer.data(fl + 272);
    const auto *fl_273 = buffer.data(fl + 273);
    const auto *fl_274 = buffer.data(fl + 274);
    const auto *fl_275 = buffer.data(fl + 275);
    const auto *fl_276 = buffer.data(fl + 276);
    const auto *fl_277 = buffer.data(fl + 277);
    const auto *fl_278 = buffer.data(fl + 278);
    const auto *fl_279 = buffer.data(fl + 279);
    const auto *fl_280 = buffer.data(fl + 280);
    const auto *fl_281 = buffer.data(fl + 281);
    const auto *fl_282 = buffer.data(fl + 282);
    const auto *fl_283 = buffer.data(fl + 283);
    const auto *fl_284 = buffer.data(fl + 284);
    const auto *fl_285 = buffer.data(fl + 285);
    const auto *fl_286 = buffer.data(fl + 286);
    const auto *fl_287 = buffer.data(fl + 287);
    const auto *fl_288 = buffer.data(fl + 288);
    const auto *fl_289 = buffer.data(fl + 289);
    const auto *fl_290 = buffer.data(fl + 290);
    const auto *fl_291 = buffer.data(fl + 291);
    const auto *fl_292 = buffer.data(fl + 292);
    const auto *fl_293 = buffer.data(fl + 293);
    const auto *fl_294 = buffer.data(fl + 294);
    const auto *fl_295 = buffer.data(fl + 295);
    const auto *fl_296 = buffer.data(fl + 296);
    const auto *fl_297 = buffer.data(fl + 297);
    const auto *fl_298 = buffer.data(fl + 298);
    const auto *fl_299 = buffer.data(fl + 299);

    const auto *hl_150 = buffer.data(hl + 150);
    const auto *hl_151 = buffer.data(hl + 151);
    const auto *hl_152 = buffer.data(hl + 152);
    const auto *hl_153 = buffer.data(hl + 153);
    const auto *hl_154 = buffer.data(hl + 154);
    const auto *hl_155 = buffer.data(hl + 155);
    const auto *hl_156 = buffer.data(hl + 156);
    const auto *hl_157 = buffer.data(hl + 157);
    const auto *hl_158 = buffer.data(hl + 158);
    const auto *hl_159 = buffer.data(hl + 159);
    const auto *hl_160 = buffer.data(hl + 160);
    const auto *hl_161 = buffer.data(hl + 161);
    const auto *hl_162 = buffer.data(hl + 162);
    const auto *hl_163 = buffer.data(hl + 163);
    const auto *hl_164 = buffer.data(hl + 164);
    const auto *hl_165 = buffer.data(hl + 165);
    const auto *hl_166 = buffer.data(hl + 166);
    const auto *hl_167 = buffer.data(hl + 167);
    const auto *hl_168 = buffer.data(hl + 168);
    const auto *hl_169 = buffer.data(hl + 169);
    const auto *hl_170 = buffer.data(hl + 170);
    const auto *hl_171 = buffer.data(hl + 171);
    const auto *hl_172 = buffer.data(hl + 172);
    const auto *hl_173 = buffer.data(hl + 173);
    const auto *hl_174 = buffer.data(hl + 174);
    const auto *hl_175 = buffer.data(hl + 175);
    const auto *hl_176 = buffer.data(hl + 176);
    const auto *hl_177 = buffer.data(hl + 177);
    const auto *hl_178 = buffer.data(hl + 178);
    const auto *hl_179 = buffer.data(hl + 179);
    const auto *hl_180 = buffer.data(hl + 180);
    const auto *hl_181 = buffer.data(hl + 181);
    const auto *hl_182 = buffer.data(hl + 182);
    const auto *hl_183 = buffer.data(hl + 183);
    const auto *hl_184 = buffer.data(hl + 184);
    const auto *hl_185 = buffer.data(hl + 185);
    const auto *hl_186 = buffer.data(hl + 186);
    const auto *hl_187 = buffer.data(hl + 187);
    const auto *hl_188 = buffer.data(hl + 188);
    const auto *hl_189 = buffer.data(hl + 189);
    const auto *hl_190 = buffer.data(hl + 190);
    const auto *hl_191 = buffer.data(hl + 191);
    const auto *hl_192 = buffer.data(hl + 192);
    const auto *hl_193 = buffer.data(hl + 193);
    const auto *hl_194 = buffer.data(hl + 194);
    const auto *hl_195 = buffer.data(hl + 195);
    const auto *hl_196 = buffer.data(hl + 196);
    const auto *hl_197 = buffer.data(hl + 197);
    const auto *hl_198 = buffer.data(hl + 198);
    const auto *hl_199 = buffer.data(hl + 199);
    const auto *hl_200 = buffer.data(hl + 200);
    const auto *hl_201 = buffer.data(hl + 201);
    const auto *hl_202 = buffer.data(hl + 202);
    const auto *hl_203 = buffer.data(hl + 203);
    const auto *hl_204 = buffer.data(hl + 204);
    const auto *hl_205 = buffer.data(hl + 205);
    const auto *hl_206 = buffer.data(hl + 206);
    const auto *hl_207 = buffer.data(hl + 207);
    const auto *hl_208 = buffer.data(hl + 208);
    const auto *hl_209 = buffer.data(hl + 209);
    const auto *hl_210 = buffer.data(hl + 210);
    const auto *hl_211 = buffer.data(hl + 211);
    const auto *hl_212 = buffer.data(hl + 212);
    const auto *hl_213 = buffer.data(hl + 213);
    const auto *hl_214 = buffer.data(hl + 214);
    const auto *hl_215 = buffer.data(hl + 215);
    const auto *hl_216 = buffer.data(hl + 216);
    const auto *hl_217 = buffer.data(hl + 217);
    const auto *hl_218 = buffer.data(hl + 218);
    const auto *hl_219 = buffer.data(hl + 219);
    const auto *hl_220 = buffer.data(hl + 220);
    const auto *hl_221 = buffer.data(hl + 221);
    const auto *hl_222 = buffer.data(hl + 222);
    const auto *hl_223 = buffer.data(hl + 223);
    const auto *hl_224 = buffer.data(hl + 224);
    const auto *hl_225 = buffer.data(hl + 225);
    const auto *hl_226 = buffer.data(hl + 226);
    const auto *hl_227 = buffer.data(hl + 227);
    const auto *hl_228 = buffer.data(hl + 228);
    const auto *hl_229 = buffer.data(hl + 229);
    const auto *hl_230 = buffer.data(hl + 230);
    const auto *hl_231 = buffer.data(hl + 231);
    const auto *hl_232 = buffer.data(hl + 232);
    const auto *hl_233 = buffer.data(hl + 233);
    const auto *hl_234 = buffer.data(hl + 234);
    const auto *hl_235 = buffer.data(hl + 235);
    const auto *hl_236 = buffer.data(hl + 236);
    const auto *hl_237 = buffer.data(hl + 237);
    const auto *hl_238 = buffer.data(hl + 238);
    const auto *hl_239 = buffer.data(hl + 239);
    const auto *hl_240 = buffer.data(hl + 240);
    const auto *hl_241 = buffer.data(hl + 241);
    const auto *hl_242 = buffer.data(hl + 242);
    const auto *hl_243 = buffer.data(hl + 243);
    const auto *hl_244 = buffer.data(hl + 244);
    const auto *hl_245 = buffer.data(hl + 245);
    const auto *hl_246 = buffer.data(hl + 246);
    const auto *hl_247 = buffer.data(hl + 247);
    const auto *hl_248 = buffer.data(hl + 248);
    const auto *hl_249 = buffer.data(hl + 249);
    const auto *hl_250 = buffer.data(hl + 250);
    const auto *hl_251 = buffer.data(hl + 251);
    const auto *hl_252 = buffer.data(hl + 252);
    const auto *hl_253 = buffer.data(hl + 253);
    const auto *hl_254 = buffer.data(hl + 254);
    const auto *hl_255 = buffer.data(hl + 255);
    const auto *hl_256 = buffer.data(hl + 256);
    const auto *hl_257 = buffer.data(hl + 257);
    const auto *hl_258 = buffer.data(hl + 258);
    const auto *hl_259 = buffer.data(hl + 259);
    const auto *hl_260 = buffer.data(hl + 260);
    const auto *hl_261 = buffer.data(hl + 261);
    const auto *hl_262 = buffer.data(hl + 262);
    const auto *hl_263 = buffer.data(hl + 263);
    const auto *hl_264 = buffer.data(hl + 264);
    const auto *hl_265 = buffer.data(hl + 265);
    const auto *hl_266 = buffer.data(hl + 266);
    const auto *hl_267 = buffer.data(hl + 267);
    const auto *hl_268 = buffer.data(hl + 268);
    const auto *hl_269 = buffer.data(hl + 269);
    const auto *hl_270 = buffer.data(hl + 270);
    const auto *hl_271 = buffer.data(hl + 271);
    const auto *hl_272 = buffer.data(hl + 272);
    const auto *hl_273 = buffer.data(hl + 273);
    const auto *hl_274 = buffer.data(hl + 274);
    const auto *hl_275 = buffer.data(hl + 275);
    const auto *hl_276 = buffer.data(hl + 276);
    const auto *hl_277 = buffer.data(hl + 277);
    const auto *hl_278 = buffer.data(hl + 278);
    const auto *hl_279 = buffer.data(hl + 279);
    const auto *hl_280 = buffer.data(hl + 280);
    const auto *hl_281 = buffer.data(hl + 281);
    const auto *hl_282 = buffer.data(hl + 282);
    const auto *hl_283 = buffer.data(hl + 283);
    const auto *hl_284 = buffer.data(hl + 284);
    const auto *hl_285 = buffer.data(hl + 285);
    const auto *hl_286 = buffer.data(hl + 286);
    const auto *hl_287 = buffer.data(hl + 287);
    const auto *hl_288 = buffer.data(hl + 288);
    const auto *hl_289 = buffer.data(hl + 289);
    const auto *hl_290 = buffer.data(hl + 290);
    const auto *hl_291 = buffer.data(hl + 291);
    const auto *hl_292 = buffer.data(hl + 292);
    const auto *hl_293 = buffer.data(hl + 293);
    const auto *hl_294 = buffer.data(hl + 294);
    const auto *hl_295 = buffer.data(hl + 295);
    const auto *hl_296 = buffer.data(hl + 296);
    const auto *hl_297 = buffer.data(hl + 297);
    const auto *hl_298 = buffer.data(hl + 298);
    const auto *hl_299 = buffer.data(hl + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, fl_150, fl_151, fl_152, fl_153, \
                         fl_154, hl_150, hl_151, hl_152, hl_153, \
                         hl_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -2.0 * fl_150[k]
                   + f_0 * hl_150[k];

        t_151[k] = -2.0 * fl_151[k]
                   + f_0 * hl_151[k];

        t_152[k] = -2.0 * fl_152[k]
                   + f_0 * hl_152[k];

        t_153[k] = -2.0 * fl_153[k]
                   + f_0 * hl_153[k];

        t_154[k] = -2.0 * fl_154[k]
                   + f_0 * hl_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, fl_155, fl_156, fl_157, fl_158, \
                         fl_159, hl_155, hl_156, hl_157, hl_158, \
                         hl_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -2.0 * fl_155[k]
                   + f_0 * hl_155[k];

        t_156[k] = -2.0 * fl_156[k]
                   + f_0 * hl_156[k];

        t_157[k] = -2.0 * fl_157[k]
                   + f_0 * hl_157[k];

        t_158[k] = -2.0 * fl_158[k]
                   + f_0 * hl_158[k];

        t_159[k] = -2.0 * fl_159[k]
                   + f_0 * hl_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, fl_160, fl_161, fl_162, fl_163, \
                         fl_164, hl_160, hl_161, hl_162, hl_163, \
                         hl_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -2.0 * fl_160[k]
                   + f_0 * hl_160[k];

        t_161[k] = -2.0 * fl_161[k]
                   + f_0 * hl_161[k];

        t_162[k] = -2.0 * fl_162[k]
                   + f_0 * hl_162[k];

        t_163[k] = -2.0 * fl_163[k]
                   + f_0 * hl_163[k];

        t_164[k] = -2.0 * fl_164[k]
                   + f_0 * hl_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, fl_165, fl_166, fl_167, fl_168, \
                         fl_169, hl_165, hl_166, hl_167, hl_168, \
                         hl_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -2.0 * fl_165[k]
                   + f_0 * hl_165[k];

        t_166[k] = -2.0 * fl_166[k]
                   + f_0 * hl_166[k];

        t_167[k] = -2.0 * fl_167[k]
                   + f_0 * hl_167[k];

        t_168[k] = -2.0 * fl_168[k]
                   + f_0 * hl_168[k];

        t_169[k] = -2.0 * fl_169[k]
                   + f_0 * hl_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, fl_170, fl_171, fl_172, fl_173, \
                         fl_174, hl_170, hl_171, hl_172, hl_173, \
                         hl_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -2.0 * fl_170[k]
                   + f_0 * hl_170[k];

        t_171[k] = -2.0 * fl_171[k]
                   + f_0 * hl_171[k];

        t_172[k] = -2.0 * fl_172[k]
                   + f_0 * hl_172[k];

        t_173[k] = -2.0 * fl_173[k]
                   + f_0 * hl_173[k];

        t_174[k] = -2.0 * fl_174[k]
                   + f_0 * hl_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, fl_175, fl_176, fl_177, fl_178, \
                         fl_179, hl_175, hl_176, hl_177, hl_178, \
                         hl_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -2.0 * fl_175[k]
                   + f_0 * hl_175[k];

        t_176[k] = -2.0 * fl_176[k]
                   + f_0 * hl_176[k];

        t_177[k] = -2.0 * fl_177[k]
                   + f_0 * hl_177[k];

        t_178[k] = -2.0 * fl_178[k]
                   + f_0 * hl_178[k];

        t_179[k] = -2.0 * fl_179[k]
                   + f_0 * hl_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, fl_180, fl_181, fl_182, fl_183, \
                         fl_184, hl_180, hl_181, hl_182, hl_183, \
                         hl_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * fl_180[k]
                   + f_0 * hl_180[k];

        t_181[k] = -2.0 * fl_181[k]
                   + f_0 * hl_181[k];

        t_182[k] = -2.0 * fl_182[k]
                   + f_0 * hl_182[k];

        t_183[k] = -2.0 * fl_183[k]
                   + f_0 * hl_183[k];

        t_184[k] = -2.0 * fl_184[k]
                   + f_0 * hl_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, fl_185, fl_186, fl_187, fl_188, \
                         fl_189, hl_185, hl_186, hl_187, hl_188, \
                         hl_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * fl_185[k]
                   + f_0 * hl_185[k];

        t_186[k] = -2.0 * fl_186[k]
                   + f_0 * hl_186[k];

        t_187[k] = -2.0 * fl_187[k]
                   + f_0 * hl_187[k];

        t_188[k] = -2.0 * fl_188[k]
                   + f_0 * hl_188[k];

        t_189[k] = -2.0 * fl_189[k]
                   + f_0 * hl_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, fl_190, fl_191, fl_192, fl_193, \
                         fl_194, hl_190, hl_191, hl_192, hl_193, \
                         hl_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * fl_190[k]
                   + f_0 * hl_190[k];

        t_191[k] = -2.0 * fl_191[k]
                   + f_0 * hl_191[k];

        t_192[k] = -2.0 * fl_192[k]
                   + f_0 * hl_192[k];

        t_193[k] = -2.0 * fl_193[k]
                   + f_0 * hl_193[k];

        t_194[k] = -2.0 * fl_194[k]
                   + f_0 * hl_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, fl_195, fl_196, fl_197, fl_198, \
                         fl_199, hl_195, hl_196, hl_197, hl_198, \
                         hl_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * fl_195[k]
                   + f_0 * hl_195[k];

        t_196[k] = -2.0 * fl_196[k]
                   + f_0 * hl_196[k];

        t_197[k] = -2.0 * fl_197[k]
                   + f_0 * hl_197[k];

        t_198[k] = -2.0 * fl_198[k]
                   + f_0 * hl_198[k];

        t_199[k] = -2.0 * fl_199[k]
                   + f_0 * hl_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, fl_200, fl_201, fl_202, fl_203, \
                         fl_204, hl_200, hl_201, hl_202, hl_203, \
                         hl_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * fl_200[k]
                   + f_0 * hl_200[k];

        t_201[k] = -2.0 * fl_201[k]
                   + f_0 * hl_201[k];

        t_202[k] = -2.0 * fl_202[k]
                   + f_0 * hl_202[k];

        t_203[k] = -2.0 * fl_203[k]
                   + f_0 * hl_203[k];

        t_204[k] = -2.0 * fl_204[k]
                   + f_0 * hl_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, fl_205, fl_206, fl_207, fl_208, \
                         fl_209, hl_205, hl_206, hl_207, hl_208, \
                         hl_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * fl_205[k]
                   + f_0 * hl_205[k];

        t_206[k] = -2.0 * fl_206[k]
                   + f_0 * hl_206[k];

        t_207[k] = -2.0 * fl_207[k]
                   + f_0 * hl_207[k];

        t_208[k] = -2.0 * fl_208[k]
                   + f_0 * hl_208[k];

        t_209[k] = -2.0 * fl_209[k]
                   + f_0 * hl_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, fl_210, fl_211, fl_212, fl_213, \
                         fl_214, hl_210, hl_211, hl_212, hl_213, \
                         hl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * fl_210[k]
                   + f_0 * hl_210[k];

        t_211[k] = -2.0 * fl_211[k]
                   + f_0 * hl_211[k];

        t_212[k] = -2.0 * fl_212[k]
                   + f_0 * hl_212[k];

        t_213[k] = -2.0 * fl_213[k]
                   + f_0 * hl_213[k];

        t_214[k] = -2.0 * fl_214[k]
                   + f_0 * hl_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, fl_215, fl_216, fl_217, fl_218, \
                         fl_219, hl_215, hl_216, hl_217, hl_218, \
                         hl_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * fl_215[k]
                   + f_0 * hl_215[k];

        t_216[k] = -2.0 * fl_216[k]
                   + f_0 * hl_216[k];

        t_217[k] = -2.0 * fl_217[k]
                   + f_0 * hl_217[k];

        t_218[k] = -2.0 * fl_218[k]
                   + f_0 * hl_218[k];

        t_219[k] = -2.0 * fl_219[k]
                   + f_0 * hl_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, fl_220, fl_221, fl_222, fl_223, \
                         fl_224, hl_220, hl_221, hl_222, hl_223, \
                         hl_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * fl_220[k]
                   + f_0 * hl_220[k];

        t_221[k] = -2.0 * fl_221[k]
                   + f_0 * hl_221[k];

        t_222[k] = -2.0 * fl_222[k]
                   + f_0 * hl_222[k];

        t_223[k] = -2.0 * fl_223[k]
                   + f_0 * hl_223[k];

        t_224[k] = -2.0 * fl_224[k]
                   + f_0 * hl_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, fl_225, fl_226, fl_227, fl_228, \
                         fl_229, hl_225, hl_226, hl_227, hl_228, \
                         hl_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -2.0 * fl_225[k]
                   + f_0 * hl_225[k];

        t_226[k] = -2.0 * fl_226[k]
                   + f_0 * hl_226[k];

        t_227[k] = -2.0 * fl_227[k]
                   + f_0 * hl_227[k];

        t_228[k] = -2.0 * fl_228[k]
                   + f_0 * hl_228[k];

        t_229[k] = -2.0 * fl_229[k]
                   + f_0 * hl_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, fl_230, fl_231, fl_232, fl_233, \
                         fl_234, hl_230, hl_231, hl_232, hl_233, \
                         hl_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -2.0 * fl_230[k]
                   + f_0 * hl_230[k];

        t_231[k] = -2.0 * fl_231[k]
                   + f_0 * hl_231[k];

        t_232[k] = -2.0 * fl_232[k]
                   + f_0 * hl_232[k];

        t_233[k] = -2.0 * fl_233[k]
                   + f_0 * hl_233[k];

        t_234[k] = -2.0 * fl_234[k]
                   + f_0 * hl_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, fl_235, fl_236, fl_237, fl_238, \
                         fl_239, hl_235, hl_236, hl_237, hl_238, \
                         hl_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -2.0 * fl_235[k]
                   + f_0 * hl_235[k];

        t_236[k] = -2.0 * fl_236[k]
                   + f_0 * hl_236[k];

        t_237[k] = -2.0 * fl_237[k]
                   + f_0 * hl_237[k];

        t_238[k] = -2.0 * fl_238[k]
                   + f_0 * hl_238[k];

        t_239[k] = -2.0 * fl_239[k]
                   + f_0 * hl_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, fl_240, fl_241, fl_242, fl_243, \
                         fl_244, hl_240, hl_241, hl_242, hl_243, \
                         hl_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -2.0 * fl_240[k]
                   + f_0 * hl_240[k];

        t_241[k] = -2.0 * fl_241[k]
                   + f_0 * hl_241[k];

        t_242[k] = -2.0 * fl_242[k]
                   + f_0 * hl_242[k];

        t_243[k] = -2.0 * fl_243[k]
                   + f_0 * hl_243[k];

        t_244[k] = -2.0 * fl_244[k]
                   + f_0 * hl_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, fl_245, fl_246, fl_247, fl_248, \
                         fl_249, hl_245, hl_246, hl_247, hl_248, \
                         hl_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -2.0 * fl_245[k]
                   + f_0 * hl_245[k];

        t_246[k] = -2.0 * fl_246[k]
                   + f_0 * hl_246[k];

        t_247[k] = -2.0 * fl_247[k]
                   + f_0 * hl_247[k];

        t_248[k] = -2.0 * fl_248[k]
                   + f_0 * hl_248[k];

        t_249[k] = -2.0 * fl_249[k]
                   + f_0 * hl_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, fl_250, fl_251, fl_252, fl_253, \
                         fl_254, hl_250, hl_251, hl_252, hl_253, \
                         hl_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * fl_250[k]
                   + f_0 * hl_250[k];

        t_251[k] = -2.0 * fl_251[k]
                   + f_0 * hl_251[k];

        t_252[k] = -2.0 * fl_252[k]
                   + f_0 * hl_252[k];

        t_253[k] = -2.0 * fl_253[k]
                   + f_0 * hl_253[k];

        t_254[k] = -2.0 * fl_254[k]
                   + f_0 * hl_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, fl_255, fl_256, fl_257, fl_258, \
                         fl_259, hl_255, hl_256, hl_257, hl_258, \
                         hl_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * fl_255[k]
                   + f_0 * hl_255[k];

        t_256[k] = -2.0 * fl_256[k]
                   + f_0 * hl_256[k];

        t_257[k] = -2.0 * fl_257[k]
                   + f_0 * hl_257[k];

        t_258[k] = -2.0 * fl_258[k]
                   + f_0 * hl_258[k];

        t_259[k] = -2.0 * fl_259[k]
                   + f_0 * hl_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, fl_260, fl_261, fl_262, fl_263, \
                         fl_264, hl_260, hl_261, hl_262, hl_263, \
                         hl_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * fl_260[k]
                   + f_0 * hl_260[k];

        t_261[k] = -2.0 * fl_261[k]
                   + f_0 * hl_261[k];

        t_262[k] = -2.0 * fl_262[k]
                   + f_0 * hl_262[k];

        t_263[k] = -2.0 * fl_263[k]
                   + f_0 * hl_263[k];

        t_264[k] = -2.0 * fl_264[k]
                   + f_0 * hl_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, fl_265, fl_266, fl_267, fl_268, \
                         fl_269, hl_265, hl_266, hl_267, hl_268, \
                         hl_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * fl_265[k]
                   + f_0 * hl_265[k];

        t_266[k] = -2.0 * fl_266[k]
                   + f_0 * hl_266[k];

        t_267[k] = -2.0 * fl_267[k]
                   + f_0 * hl_267[k];

        t_268[k] = -2.0 * fl_268[k]
                   + f_0 * hl_268[k];

        t_269[k] = -2.0 * fl_269[k]
                   + f_0 * hl_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, fl_270, fl_271, fl_272, fl_273, \
                         fl_274, hl_270, hl_271, hl_272, hl_273, \
                         hl_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -fl_270[k]
                   + f_0 * hl_270[k];

        t_271[k] = -fl_271[k]
                   + f_0 * hl_271[k];

        t_272[k] = -fl_272[k]
                   + f_0 * hl_272[k];

        t_273[k] = -fl_273[k]
                   + f_0 * hl_273[k];

        t_274[k] = -fl_274[k]
                   + f_0 * hl_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, fl_275, fl_276, fl_277, fl_278, \
                         fl_279, hl_275, hl_276, hl_277, hl_278, \
                         hl_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -fl_275[k]
                   + f_0 * hl_275[k];

        t_276[k] = -fl_276[k]
                   + f_0 * hl_276[k];

        t_277[k] = -fl_277[k]
                   + f_0 * hl_277[k];

        t_278[k] = -fl_278[k]
                   + f_0 * hl_278[k];

        t_279[k] = -fl_279[k]
                   + f_0 * hl_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, fl_280, fl_281, fl_282, fl_283, \
                         fl_284, hl_280, hl_281, hl_282, hl_283, \
                         hl_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -fl_280[k]
                   + f_0 * hl_280[k];

        t_281[k] = -fl_281[k]
                   + f_0 * hl_281[k];

        t_282[k] = -fl_282[k]
                   + f_0 * hl_282[k];

        t_283[k] = -fl_283[k]
                   + f_0 * hl_283[k];

        t_284[k] = -fl_284[k]
                   + f_0 * hl_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, fl_285, fl_286, fl_287, fl_288, \
                         fl_289, hl_285, hl_286, hl_287, hl_288, \
                         hl_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -fl_285[k]
                   + f_0 * hl_285[k];

        t_286[k] = -fl_286[k]
                   + f_0 * hl_286[k];

        t_287[k] = -fl_287[k]
                   + f_0 * hl_287[k];

        t_288[k] = -fl_288[k]
                   + f_0 * hl_288[k];

        t_289[k] = -fl_289[k]
                   + f_0 * hl_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, fl_290, fl_291, fl_292, fl_293, \
                         fl_294, hl_290, hl_291, hl_292, hl_293, \
                         hl_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -fl_290[k]
                   + f_0 * hl_290[k];

        t_291[k] = -fl_291[k]
                   + f_0 * hl_291[k];

        t_292[k] = -fl_292[k]
                   + f_0 * hl_292[k];

        t_293[k] = -fl_293[k]
                   + f_0 * hl_293[k];

        t_294[k] = -fl_294[k]
                   + f_0 * hl_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, fl_295, fl_296, fl_297, fl_298, \
                         fl_299, hl_295, hl_296, hl_297, hl_298, \
                         hl_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -fl_295[k]
                   + f_0 * hl_295[k];

        t_296[k] = -fl_296[k]
                   + f_0 * hl_296[k];

        t_297[k] = -fl_297[k]
                   + f_0 * hl_297[k];

        t_298[k] = -fl_298[k]
                   + f_0 * hl_298[k];

        t_299[k] = -fl_299[k]
                   + f_0 * hl_299[k];
    }
}

static auto
compute_prim_geom_10_gl_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
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

    const auto *fl_300 = buffer.data(fl + 300);
    const auto *fl_301 = buffer.data(fl + 301);
    const auto *fl_302 = buffer.data(fl + 302);
    const auto *fl_303 = buffer.data(fl + 303);
    const auto *fl_304 = buffer.data(fl + 304);
    const auto *fl_305 = buffer.data(fl + 305);
    const auto *fl_306 = buffer.data(fl + 306);
    const auto *fl_307 = buffer.data(fl + 307);
    const auto *fl_308 = buffer.data(fl + 308);
    const auto *fl_309 = buffer.data(fl + 309);
    const auto *fl_310 = buffer.data(fl + 310);
    const auto *fl_311 = buffer.data(fl + 311);
    const auto *fl_312 = buffer.data(fl + 312);
    const auto *fl_313 = buffer.data(fl + 313);
    const auto *fl_314 = buffer.data(fl + 314);
    const auto *fl_315 = buffer.data(fl + 315);
    const auto *fl_316 = buffer.data(fl + 316);
    const auto *fl_317 = buffer.data(fl + 317);
    const auto *fl_318 = buffer.data(fl + 318);
    const auto *fl_319 = buffer.data(fl + 319);
    const auto *fl_320 = buffer.data(fl + 320);
    const auto *fl_321 = buffer.data(fl + 321);
    const auto *fl_322 = buffer.data(fl + 322);
    const auto *fl_323 = buffer.data(fl + 323);
    const auto *fl_324 = buffer.data(fl + 324);
    const auto *fl_325 = buffer.data(fl + 325);
    const auto *fl_326 = buffer.data(fl + 326);
    const auto *fl_327 = buffer.data(fl + 327);
    const auto *fl_328 = buffer.data(fl + 328);
    const auto *fl_329 = buffer.data(fl + 329);
    const auto *fl_330 = buffer.data(fl + 330);
    const auto *fl_331 = buffer.data(fl + 331);
    const auto *fl_332 = buffer.data(fl + 332);
    const auto *fl_333 = buffer.data(fl + 333);
    const auto *fl_334 = buffer.data(fl + 334);
    const auto *fl_335 = buffer.data(fl + 335);
    const auto *fl_336 = buffer.data(fl + 336);
    const auto *fl_337 = buffer.data(fl + 337);
    const auto *fl_338 = buffer.data(fl + 338);
    const auto *fl_339 = buffer.data(fl + 339);
    const auto *fl_340 = buffer.data(fl + 340);
    const auto *fl_341 = buffer.data(fl + 341);
    const auto *fl_342 = buffer.data(fl + 342);
    const auto *fl_343 = buffer.data(fl + 343);
    const auto *fl_344 = buffer.data(fl + 344);
    const auto *fl_345 = buffer.data(fl + 345);
    const auto *fl_346 = buffer.data(fl + 346);
    const auto *fl_347 = buffer.data(fl + 347);
    const auto *fl_348 = buffer.data(fl + 348);
    const auto *fl_349 = buffer.data(fl + 349);
    const auto *fl_350 = buffer.data(fl + 350);
    const auto *fl_351 = buffer.data(fl + 351);
    const auto *fl_352 = buffer.data(fl + 352);
    const auto *fl_353 = buffer.data(fl + 353);
    const auto *fl_354 = buffer.data(fl + 354);
    const auto *fl_355 = buffer.data(fl + 355);
    const auto *fl_356 = buffer.data(fl + 356);
    const auto *fl_357 = buffer.data(fl + 357);
    const auto *fl_358 = buffer.data(fl + 358);
    const auto *fl_359 = buffer.data(fl + 359);
    const auto *fl_360 = buffer.data(fl + 360);
    const auto *fl_361 = buffer.data(fl + 361);
    const auto *fl_362 = buffer.data(fl + 362);
    const auto *fl_363 = buffer.data(fl + 363);
    const auto *fl_364 = buffer.data(fl + 364);
    const auto *fl_365 = buffer.data(fl + 365);
    const auto *fl_366 = buffer.data(fl + 366);
    const auto *fl_367 = buffer.data(fl + 367);
    const auto *fl_368 = buffer.data(fl + 368);
    const auto *fl_369 = buffer.data(fl + 369);
    const auto *fl_370 = buffer.data(fl + 370);
    const auto *fl_371 = buffer.data(fl + 371);
    const auto *fl_372 = buffer.data(fl + 372);
    const auto *fl_373 = buffer.data(fl + 373);
    const auto *fl_374 = buffer.data(fl + 374);
    const auto *fl_375 = buffer.data(fl + 375);
    const auto *fl_376 = buffer.data(fl + 376);
    const auto *fl_377 = buffer.data(fl + 377);
    const auto *fl_378 = buffer.data(fl + 378);
    const auto *fl_379 = buffer.data(fl + 379);
    const auto *fl_380 = buffer.data(fl + 380);
    const auto *fl_381 = buffer.data(fl + 381);
    const auto *fl_382 = buffer.data(fl + 382);
    const auto *fl_383 = buffer.data(fl + 383);
    const auto *fl_384 = buffer.data(fl + 384);
    const auto *fl_385 = buffer.data(fl + 385);
    const auto *fl_386 = buffer.data(fl + 386);
    const auto *fl_387 = buffer.data(fl + 387);
    const auto *fl_388 = buffer.data(fl + 388);
    const auto *fl_389 = buffer.data(fl + 389);
    const auto *fl_390 = buffer.data(fl + 390);
    const auto *fl_391 = buffer.data(fl + 391);
    const auto *fl_392 = buffer.data(fl + 392);
    const auto *fl_393 = buffer.data(fl + 393);
    const auto *fl_394 = buffer.data(fl + 394);
    const auto *fl_395 = buffer.data(fl + 395);
    const auto *fl_396 = buffer.data(fl + 396);
    const auto *fl_397 = buffer.data(fl + 397);
    const auto *fl_398 = buffer.data(fl + 398);
    const auto *fl_399 = buffer.data(fl + 399);
    const auto *fl_400 = buffer.data(fl + 400);
    const auto *fl_401 = buffer.data(fl + 401);
    const auto *fl_402 = buffer.data(fl + 402);
    const auto *fl_403 = buffer.data(fl + 403);
    const auto *fl_404 = buffer.data(fl + 404);
    const auto *fl_405 = buffer.data(fl + 405);
    const auto *fl_406 = buffer.data(fl + 406);
    const auto *fl_407 = buffer.data(fl + 407);
    const auto *fl_408 = buffer.data(fl + 408);
    const auto *fl_409 = buffer.data(fl + 409);
    const auto *fl_410 = buffer.data(fl + 410);
    const auto *fl_411 = buffer.data(fl + 411);
    const auto *fl_412 = buffer.data(fl + 412);
    const auto *fl_413 = buffer.data(fl + 413);
    const auto *fl_414 = buffer.data(fl + 414);
    const auto *fl_415 = buffer.data(fl + 415);
    const auto *fl_416 = buffer.data(fl + 416);
    const auto *fl_417 = buffer.data(fl + 417);
    const auto *fl_418 = buffer.data(fl + 418);
    const auto *fl_419 = buffer.data(fl + 419);
    const auto *fl_420 = buffer.data(fl + 420);
    const auto *fl_421 = buffer.data(fl + 421);
    const auto *fl_422 = buffer.data(fl + 422);
    const auto *fl_423 = buffer.data(fl + 423);
    const auto *fl_424 = buffer.data(fl + 424);
    const auto *fl_425 = buffer.data(fl + 425);
    const auto *fl_426 = buffer.data(fl + 426);
    const auto *fl_427 = buffer.data(fl + 427);
    const auto *fl_428 = buffer.data(fl + 428);
    const auto *fl_429 = buffer.data(fl + 429);
    const auto *fl_430 = buffer.data(fl + 430);
    const auto *fl_431 = buffer.data(fl + 431);
    const auto *fl_432 = buffer.data(fl + 432);
    const auto *fl_433 = buffer.data(fl + 433);
    const auto *fl_434 = buffer.data(fl + 434);
    const auto *fl_435 = buffer.data(fl + 435);
    const auto *fl_436 = buffer.data(fl + 436);
    const auto *fl_437 = buffer.data(fl + 437);
    const auto *fl_438 = buffer.data(fl + 438);
    const auto *fl_439 = buffer.data(fl + 439);
    const auto *fl_440 = buffer.data(fl + 440);
    const auto *fl_441 = buffer.data(fl + 441);
    const auto *fl_442 = buffer.data(fl + 442);
    const auto *fl_443 = buffer.data(fl + 443);
    const auto *fl_444 = buffer.data(fl + 444);
    const auto *fl_445 = buffer.data(fl + 445);
    const auto *fl_446 = buffer.data(fl + 446);
    const auto *fl_447 = buffer.data(fl + 447);
    const auto *fl_448 = buffer.data(fl + 448);
    const auto *fl_449 = buffer.data(fl + 449);

    const auto *hl_300 = buffer.data(hl + 300);
    const auto *hl_301 = buffer.data(hl + 301);
    const auto *hl_302 = buffer.data(hl + 302);
    const auto *hl_303 = buffer.data(hl + 303);
    const auto *hl_304 = buffer.data(hl + 304);
    const auto *hl_305 = buffer.data(hl + 305);
    const auto *hl_306 = buffer.data(hl + 306);
    const auto *hl_307 = buffer.data(hl + 307);
    const auto *hl_308 = buffer.data(hl + 308);
    const auto *hl_309 = buffer.data(hl + 309);
    const auto *hl_310 = buffer.data(hl + 310);
    const auto *hl_311 = buffer.data(hl + 311);
    const auto *hl_312 = buffer.data(hl + 312);
    const auto *hl_313 = buffer.data(hl + 313);
    const auto *hl_314 = buffer.data(hl + 314);
    const auto *hl_315 = buffer.data(hl + 315);
    const auto *hl_316 = buffer.data(hl + 316);
    const auto *hl_317 = buffer.data(hl + 317);
    const auto *hl_318 = buffer.data(hl + 318);
    const auto *hl_319 = buffer.data(hl + 319);
    const auto *hl_320 = buffer.data(hl + 320);
    const auto *hl_321 = buffer.data(hl + 321);
    const auto *hl_322 = buffer.data(hl + 322);
    const auto *hl_323 = buffer.data(hl + 323);
    const auto *hl_324 = buffer.data(hl + 324);
    const auto *hl_325 = buffer.data(hl + 325);
    const auto *hl_326 = buffer.data(hl + 326);
    const auto *hl_327 = buffer.data(hl + 327);
    const auto *hl_328 = buffer.data(hl + 328);
    const auto *hl_329 = buffer.data(hl + 329);
    const auto *hl_330 = buffer.data(hl + 330);
    const auto *hl_331 = buffer.data(hl + 331);
    const auto *hl_332 = buffer.data(hl + 332);
    const auto *hl_333 = buffer.data(hl + 333);
    const auto *hl_334 = buffer.data(hl + 334);
    const auto *hl_335 = buffer.data(hl + 335);
    const auto *hl_336 = buffer.data(hl + 336);
    const auto *hl_337 = buffer.data(hl + 337);
    const auto *hl_338 = buffer.data(hl + 338);
    const auto *hl_339 = buffer.data(hl + 339);
    const auto *hl_340 = buffer.data(hl + 340);
    const auto *hl_341 = buffer.data(hl + 341);
    const auto *hl_342 = buffer.data(hl + 342);
    const auto *hl_343 = buffer.data(hl + 343);
    const auto *hl_344 = buffer.data(hl + 344);
    const auto *hl_345 = buffer.data(hl + 345);
    const auto *hl_346 = buffer.data(hl + 346);
    const auto *hl_347 = buffer.data(hl + 347);
    const auto *hl_348 = buffer.data(hl + 348);
    const auto *hl_349 = buffer.data(hl + 349);
    const auto *hl_350 = buffer.data(hl + 350);
    const auto *hl_351 = buffer.data(hl + 351);
    const auto *hl_352 = buffer.data(hl + 352);
    const auto *hl_353 = buffer.data(hl + 353);
    const auto *hl_354 = buffer.data(hl + 354);
    const auto *hl_355 = buffer.data(hl + 355);
    const auto *hl_356 = buffer.data(hl + 356);
    const auto *hl_357 = buffer.data(hl + 357);
    const auto *hl_358 = buffer.data(hl + 358);
    const auto *hl_359 = buffer.data(hl + 359);
    const auto *hl_360 = buffer.data(hl + 360);
    const auto *hl_361 = buffer.data(hl + 361);
    const auto *hl_362 = buffer.data(hl + 362);
    const auto *hl_363 = buffer.data(hl + 363);
    const auto *hl_364 = buffer.data(hl + 364);
    const auto *hl_365 = buffer.data(hl + 365);
    const auto *hl_366 = buffer.data(hl + 366);
    const auto *hl_367 = buffer.data(hl + 367);
    const auto *hl_368 = buffer.data(hl + 368);
    const auto *hl_369 = buffer.data(hl + 369);
    const auto *hl_370 = buffer.data(hl + 370);
    const auto *hl_371 = buffer.data(hl + 371);
    const auto *hl_372 = buffer.data(hl + 372);
    const auto *hl_373 = buffer.data(hl + 373);
    const auto *hl_374 = buffer.data(hl + 374);
    const auto *hl_375 = buffer.data(hl + 375);
    const auto *hl_376 = buffer.data(hl + 376);
    const auto *hl_377 = buffer.data(hl + 377);
    const auto *hl_378 = buffer.data(hl + 378);
    const auto *hl_379 = buffer.data(hl + 379);
    const auto *hl_380 = buffer.data(hl + 380);
    const auto *hl_381 = buffer.data(hl + 381);
    const auto *hl_382 = buffer.data(hl + 382);
    const auto *hl_383 = buffer.data(hl + 383);
    const auto *hl_384 = buffer.data(hl + 384);
    const auto *hl_385 = buffer.data(hl + 385);
    const auto *hl_386 = buffer.data(hl + 386);
    const auto *hl_387 = buffer.data(hl + 387);
    const auto *hl_388 = buffer.data(hl + 388);
    const auto *hl_389 = buffer.data(hl + 389);
    const auto *hl_390 = buffer.data(hl + 390);
    const auto *hl_391 = buffer.data(hl + 391);
    const auto *hl_392 = buffer.data(hl + 392);
    const auto *hl_393 = buffer.data(hl + 393);
    const auto *hl_394 = buffer.data(hl + 394);
    const auto *hl_395 = buffer.data(hl + 395);
    const auto *hl_396 = buffer.data(hl + 396);
    const auto *hl_397 = buffer.data(hl + 397);
    const auto *hl_398 = buffer.data(hl + 398);
    const auto *hl_399 = buffer.data(hl + 399);
    const auto *hl_400 = buffer.data(hl + 400);
    const auto *hl_401 = buffer.data(hl + 401);
    const auto *hl_402 = buffer.data(hl + 402);
    const auto *hl_403 = buffer.data(hl + 403);
    const auto *hl_404 = buffer.data(hl + 404);
    const auto *hl_405 = buffer.data(hl + 405);
    const auto *hl_406 = buffer.data(hl + 406);
    const auto *hl_407 = buffer.data(hl + 407);
    const auto *hl_408 = buffer.data(hl + 408);
    const auto *hl_409 = buffer.data(hl + 409);
    const auto *hl_410 = buffer.data(hl + 410);
    const auto *hl_411 = buffer.data(hl + 411);
    const auto *hl_412 = buffer.data(hl + 412);
    const auto *hl_413 = buffer.data(hl + 413);
    const auto *hl_414 = buffer.data(hl + 414);
    const auto *hl_415 = buffer.data(hl + 415);
    const auto *hl_416 = buffer.data(hl + 416);
    const auto *hl_417 = buffer.data(hl + 417);
    const auto *hl_418 = buffer.data(hl + 418);
    const auto *hl_419 = buffer.data(hl + 419);
    const auto *hl_420 = buffer.data(hl + 420);
    const auto *hl_421 = buffer.data(hl + 421);
    const auto *hl_422 = buffer.data(hl + 422);
    const auto *hl_423 = buffer.data(hl + 423);
    const auto *hl_424 = buffer.data(hl + 424);
    const auto *hl_425 = buffer.data(hl + 425);
    const auto *hl_426 = buffer.data(hl + 426);
    const auto *hl_427 = buffer.data(hl + 427);
    const auto *hl_428 = buffer.data(hl + 428);
    const auto *hl_429 = buffer.data(hl + 429);
    const auto *hl_430 = buffer.data(hl + 430);
    const auto *hl_431 = buffer.data(hl + 431);
    const auto *hl_432 = buffer.data(hl + 432);
    const auto *hl_433 = buffer.data(hl + 433);
    const auto *hl_434 = buffer.data(hl + 434);
    const auto *hl_435 = buffer.data(hl + 435);
    const auto *hl_436 = buffer.data(hl + 436);
    const auto *hl_437 = buffer.data(hl + 437);
    const auto *hl_438 = buffer.data(hl + 438);
    const auto *hl_439 = buffer.data(hl + 439);
    const auto *hl_440 = buffer.data(hl + 440);
    const auto *hl_441 = buffer.data(hl + 441);
    const auto *hl_442 = buffer.data(hl + 442);
    const auto *hl_443 = buffer.data(hl + 443);
    const auto *hl_444 = buffer.data(hl + 444);
    const auto *hl_445 = buffer.data(hl + 445);
    const auto *hl_446 = buffer.data(hl + 446);
    const auto *hl_447 = buffer.data(hl + 447);
    const auto *hl_448 = buffer.data(hl + 448);
    const auto *hl_449 = buffer.data(hl + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, fl_300, fl_301, fl_302, fl_303, \
                         fl_304, hl_300, hl_301, hl_302, hl_303, \
                         hl_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -fl_300[k]
                   + f_0 * hl_300[k];

        t_301[k] = -fl_301[k]
                   + f_0 * hl_301[k];

        t_302[k] = -fl_302[k]
                   + f_0 * hl_302[k];

        t_303[k] = -fl_303[k]
                   + f_0 * hl_303[k];

        t_304[k] = -fl_304[k]
                   + f_0 * hl_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, fl_305, fl_306, fl_307, fl_308, \
                         fl_309, hl_305, hl_306, hl_307, hl_308, \
                         hl_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -fl_305[k]
                   + f_0 * hl_305[k];

        t_306[k] = -fl_306[k]
                   + f_0 * hl_306[k];

        t_307[k] = -fl_307[k]
                   + f_0 * hl_307[k];

        t_308[k] = -fl_308[k]
                   + f_0 * hl_308[k];

        t_309[k] = -fl_309[k]
                   + f_0 * hl_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, fl_310, fl_311, fl_312, fl_313, \
                         fl_314, hl_310, hl_311, hl_312, hl_313, \
                         hl_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -fl_310[k]
                   + f_0 * hl_310[k];

        t_311[k] = -fl_311[k]
                   + f_0 * hl_311[k];

        t_312[k] = -fl_312[k]
                   + f_0 * hl_312[k];

        t_313[k] = -fl_313[k]
                   + f_0 * hl_313[k];

        t_314[k] = -fl_314[k]
                   + f_0 * hl_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, fl_315, fl_316, fl_317, fl_318, \
                         fl_319, hl_315, hl_316, hl_317, hl_318, \
                         hl_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -fl_315[k]
                   + f_0 * hl_315[k];

        t_316[k] = -fl_316[k]
                   + f_0 * hl_316[k];

        t_317[k] = -fl_317[k]
                   + f_0 * hl_317[k];

        t_318[k] = -fl_318[k]
                   + f_0 * hl_318[k];

        t_319[k] = -fl_319[k]
                   + f_0 * hl_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, fl_320, fl_321, fl_322, fl_323, \
                         fl_324, hl_320, hl_321, hl_322, hl_323, \
                         hl_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -fl_320[k]
                   + f_0 * hl_320[k];

        t_321[k] = -fl_321[k]
                   + f_0 * hl_321[k];

        t_322[k] = -fl_322[k]
                   + f_0 * hl_322[k];

        t_323[k] = -fl_323[k]
                   + f_0 * hl_323[k];

        t_324[k] = -fl_324[k]
                   + f_0 * hl_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, fl_325, fl_326, fl_327, fl_328, \
                         fl_329, hl_325, hl_326, hl_327, hl_328, \
                         hl_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -fl_325[k]
                   + f_0 * hl_325[k];

        t_326[k] = -fl_326[k]
                   + f_0 * hl_326[k];

        t_327[k] = -fl_327[k]
                   + f_0 * hl_327[k];

        t_328[k] = -fl_328[k]
                   + f_0 * hl_328[k];

        t_329[k] = -fl_329[k]
                   + f_0 * hl_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, fl_330, fl_331, fl_332, fl_333, \
                         fl_334, hl_330, hl_331, hl_332, hl_333, \
                         hl_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -fl_330[k]
                   + f_0 * hl_330[k];

        t_331[k] = -fl_331[k]
                   + f_0 * hl_331[k];

        t_332[k] = -fl_332[k]
                   + f_0 * hl_332[k];

        t_333[k] = -fl_333[k]
                   + f_0 * hl_333[k];

        t_334[k] = -fl_334[k]
                   + f_0 * hl_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, fl_335, fl_336, fl_337, fl_338, \
                         fl_339, hl_335, hl_336, hl_337, hl_338, \
                         hl_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -fl_335[k]
                   + f_0 * hl_335[k];

        t_336[k] = -fl_336[k]
                   + f_0 * hl_336[k];

        t_337[k] = -fl_337[k]
                   + f_0 * hl_337[k];

        t_338[k] = -fl_338[k]
                   + f_0 * hl_338[k];

        t_339[k] = -fl_339[k]
                   + f_0 * hl_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, fl_340, fl_341, fl_342, fl_343, \
                         fl_344, hl_340, hl_341, hl_342, hl_343, \
                         hl_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -fl_340[k]
                   + f_0 * hl_340[k];

        t_341[k] = -fl_341[k]
                   + f_0 * hl_341[k];

        t_342[k] = -fl_342[k]
                   + f_0 * hl_342[k];

        t_343[k] = -fl_343[k]
                   + f_0 * hl_343[k];

        t_344[k] = -fl_344[k]
                   + f_0 * hl_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, fl_345, fl_346, fl_347, fl_348, \
                         fl_349, hl_345, hl_346, hl_347, hl_348, \
                         hl_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -fl_345[k]
                   + f_0 * hl_345[k];

        t_346[k] = -fl_346[k]
                   + f_0 * hl_346[k];

        t_347[k] = -fl_347[k]
                   + f_0 * hl_347[k];

        t_348[k] = -fl_348[k]
                   + f_0 * hl_348[k];

        t_349[k] = -fl_349[k]
                   + f_0 * hl_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, fl_350, fl_351, fl_352, fl_353, \
                         fl_354, hl_350, hl_351, hl_352, hl_353, \
                         hl_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -fl_350[k]
                   + f_0 * hl_350[k];

        t_351[k] = -fl_351[k]
                   + f_0 * hl_351[k];

        t_352[k] = -fl_352[k]
                   + f_0 * hl_352[k];

        t_353[k] = -fl_353[k]
                   + f_0 * hl_353[k];

        t_354[k] = -fl_354[k]
                   + f_0 * hl_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, fl_355, fl_356, fl_357, fl_358, \
                         fl_359, hl_355, hl_356, hl_357, hl_358, \
                         hl_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -fl_355[k]
                   + f_0 * hl_355[k];

        t_356[k] = -fl_356[k]
                   + f_0 * hl_356[k];

        t_357[k] = -fl_357[k]
                   + f_0 * hl_357[k];

        t_358[k] = -fl_358[k]
                   + f_0 * hl_358[k];

        t_359[k] = -fl_359[k]
                   + f_0 * hl_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, fl_360, fl_361, fl_362, fl_363, \
                         fl_364, hl_360, hl_361, hl_362, hl_363, \
                         hl_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -fl_360[k]
                   + f_0 * hl_360[k];

        t_361[k] = -fl_361[k]
                   + f_0 * hl_361[k];

        t_362[k] = -fl_362[k]
                   + f_0 * hl_362[k];

        t_363[k] = -fl_363[k]
                   + f_0 * hl_363[k];

        t_364[k] = -fl_364[k]
                   + f_0 * hl_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, fl_365, fl_366, fl_367, fl_368, \
                         fl_369, hl_365, hl_366, hl_367, hl_368, \
                         hl_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -fl_365[k]
                   + f_0 * hl_365[k];

        t_366[k] = -fl_366[k]
                   + f_0 * hl_366[k];

        t_367[k] = -fl_367[k]
                   + f_0 * hl_367[k];

        t_368[k] = -fl_368[k]
                   + f_0 * hl_368[k];

        t_369[k] = -fl_369[k]
                   + f_0 * hl_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, fl_370, fl_371, fl_372, fl_373, \
                         fl_374, hl_370, hl_371, hl_372, hl_373, \
                         hl_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -fl_370[k]
                   + f_0 * hl_370[k];

        t_371[k] = -fl_371[k]
                   + f_0 * hl_371[k];

        t_372[k] = -fl_372[k]
                   + f_0 * hl_372[k];

        t_373[k] = -fl_373[k]
                   + f_0 * hl_373[k];

        t_374[k] = -fl_374[k]
                   + f_0 * hl_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, fl_375, fl_376, fl_377, fl_378, \
                         fl_379, hl_375, hl_376, hl_377, hl_378, \
                         hl_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -fl_375[k]
                   + f_0 * hl_375[k];

        t_376[k] = -fl_376[k]
                   + f_0 * hl_376[k];

        t_377[k] = -fl_377[k]
                   + f_0 * hl_377[k];

        t_378[k] = -fl_378[k]
                   + f_0 * hl_378[k];

        t_379[k] = -fl_379[k]
                   + f_0 * hl_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, fl_380, fl_381, fl_382, fl_383, \
                         fl_384, hl_380, hl_381, hl_382, hl_383, \
                         hl_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -fl_380[k]
                   + f_0 * hl_380[k];

        t_381[k] = -fl_381[k]
                   + f_0 * hl_381[k];

        t_382[k] = -fl_382[k]
                   + f_0 * hl_382[k];

        t_383[k] = -fl_383[k]
                   + f_0 * hl_383[k];

        t_384[k] = -fl_384[k]
                   + f_0 * hl_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, fl_385, fl_386, fl_387, fl_388, \
                         fl_389, hl_385, hl_386, hl_387, hl_388, \
                         hl_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -fl_385[k]
                   + f_0 * hl_385[k];

        t_386[k] = -fl_386[k]
                   + f_0 * hl_386[k];

        t_387[k] = -fl_387[k]
                   + f_0 * hl_387[k];

        t_388[k] = -fl_388[k]
                   + f_0 * hl_388[k];

        t_389[k] = -fl_389[k]
                   + f_0 * hl_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, fl_390, fl_391, fl_392, fl_393, \
                         fl_394, hl_390, hl_391, hl_392, hl_393, \
                         hl_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -fl_390[k]
                   + f_0 * hl_390[k];

        t_391[k] = -fl_391[k]
                   + f_0 * hl_391[k];

        t_392[k] = -fl_392[k]
                   + f_0 * hl_392[k];

        t_393[k] = -fl_393[k]
                   + f_0 * hl_393[k];

        t_394[k] = -fl_394[k]
                   + f_0 * hl_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, fl_395, fl_396, fl_397, fl_398, \
                         fl_399, hl_395, hl_396, hl_397, hl_398, \
                         hl_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -fl_395[k]
                   + f_0 * hl_395[k];

        t_396[k] = -fl_396[k]
                   + f_0 * hl_396[k];

        t_397[k] = -fl_397[k]
                   + f_0 * hl_397[k];

        t_398[k] = -fl_398[k]
                   + f_0 * hl_398[k];

        t_399[k] = -fl_399[k]
                   + f_0 * hl_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, fl_400, fl_401, fl_402, fl_403, \
                         fl_404, hl_400, hl_401, hl_402, hl_403, \
                         hl_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -fl_400[k]
                   + f_0 * hl_400[k];

        t_401[k] = -fl_401[k]
                   + f_0 * hl_401[k];

        t_402[k] = -fl_402[k]
                   + f_0 * hl_402[k];

        t_403[k] = -fl_403[k]
                   + f_0 * hl_403[k];

        t_404[k] = -fl_404[k]
                   + f_0 * hl_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, fl_405, fl_406, fl_407, fl_408, \
                         fl_409, hl_405, hl_406, hl_407, hl_408, \
                         hl_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -fl_405[k]
                   + f_0 * hl_405[k];

        t_406[k] = -fl_406[k]
                   + f_0 * hl_406[k];

        t_407[k] = -fl_407[k]
                   + f_0 * hl_407[k];

        t_408[k] = -fl_408[k]
                   + f_0 * hl_408[k];

        t_409[k] = -fl_409[k]
                   + f_0 * hl_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, fl_410, fl_411, fl_412, fl_413, \
                         fl_414, hl_410, hl_411, hl_412, hl_413, \
                         hl_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -fl_410[k]
                   + f_0 * hl_410[k];

        t_411[k] = -fl_411[k]
                   + f_0 * hl_411[k];

        t_412[k] = -fl_412[k]
                   + f_0 * hl_412[k];

        t_413[k] = -fl_413[k]
                   + f_0 * hl_413[k];

        t_414[k] = -fl_414[k]
                   + f_0 * hl_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, fl_415, fl_416, fl_417, fl_418, \
                         fl_419, hl_415, hl_416, hl_417, hl_418, \
                         hl_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -fl_415[k]
                   + f_0 * hl_415[k];

        t_416[k] = -fl_416[k]
                   + f_0 * hl_416[k];

        t_417[k] = -fl_417[k]
                   + f_0 * hl_417[k];

        t_418[k] = -fl_418[k]
                   + f_0 * hl_418[k];

        t_419[k] = -fl_419[k]
                   + f_0 * hl_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, fl_420, fl_421, fl_422, fl_423, \
                         fl_424, hl_420, hl_421, hl_422, hl_423, \
                         hl_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -fl_420[k]
                   + f_0 * hl_420[k];

        t_421[k] = -fl_421[k]
                   + f_0 * hl_421[k];

        t_422[k] = -fl_422[k]
                   + f_0 * hl_422[k];

        t_423[k] = -fl_423[k]
                   + f_0 * hl_423[k];

        t_424[k] = -fl_424[k]
                   + f_0 * hl_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, fl_425, fl_426, fl_427, fl_428, \
                         fl_429, hl_425, hl_426, hl_427, hl_428, \
                         hl_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -fl_425[k]
                   + f_0 * hl_425[k];

        t_426[k] = -fl_426[k]
                   + f_0 * hl_426[k];

        t_427[k] = -fl_427[k]
                   + f_0 * hl_427[k];

        t_428[k] = -fl_428[k]
                   + f_0 * hl_428[k];

        t_429[k] = -fl_429[k]
                   + f_0 * hl_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, fl_430, fl_431, fl_432, fl_433, \
                         fl_434, hl_430, hl_431, hl_432, hl_433, \
                         hl_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -fl_430[k]
                   + f_0 * hl_430[k];

        t_431[k] = -fl_431[k]
                   + f_0 * hl_431[k];

        t_432[k] = -fl_432[k]
                   + f_0 * hl_432[k];

        t_433[k] = -fl_433[k]
                   + f_0 * hl_433[k];

        t_434[k] = -fl_434[k]
                   + f_0 * hl_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, fl_435, fl_436, fl_437, fl_438, \
                         fl_439, hl_435, hl_436, hl_437, hl_438, \
                         hl_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -fl_435[k]
                   + f_0 * hl_435[k];

        t_436[k] = -fl_436[k]
                   + f_0 * hl_436[k];

        t_437[k] = -fl_437[k]
                   + f_0 * hl_437[k];

        t_438[k] = -fl_438[k]
                   + f_0 * hl_438[k];

        t_439[k] = -fl_439[k]
                   + f_0 * hl_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, fl_440, fl_441, fl_442, fl_443, \
                         fl_444, hl_440, hl_441, hl_442, hl_443, \
                         hl_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -fl_440[k]
                   + f_0 * hl_440[k];

        t_441[k] = -fl_441[k]
                   + f_0 * hl_441[k];

        t_442[k] = -fl_442[k]
                   + f_0 * hl_442[k];

        t_443[k] = -fl_443[k]
                   + f_0 * hl_443[k];

        t_444[k] = -fl_444[k]
                   + f_0 * hl_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, fl_445, fl_446, fl_447, fl_448, \
                         fl_449, hl_445, hl_446, hl_447, hl_448, \
                         hl_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -fl_445[k]
                   + f_0 * hl_445[k];

        t_446[k] = -fl_446[k]
                   + f_0 * hl_446[k];

        t_447[k] = -fl_447[k]
                   + f_0 * hl_447[k];

        t_448[k] = -fl_448[k]
                   + f_0 * hl_448[k];

        t_449[k] = -fl_449[k]
                   + f_0 * hl_449[k];
    }
}

static auto
compute_prim_geom_10_gl_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t ncols,
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

    const auto *hl_450 = buffer.data(hl + 450);
    const auto *hl_451 = buffer.data(hl + 451);
    const auto *hl_452 = buffer.data(hl + 452);
    const auto *hl_453 = buffer.data(hl + 453);
    const auto *hl_454 = buffer.data(hl + 454);
    const auto *hl_455 = buffer.data(hl + 455);
    const auto *hl_456 = buffer.data(hl + 456);
    const auto *hl_457 = buffer.data(hl + 457);
    const auto *hl_458 = buffer.data(hl + 458);
    const auto *hl_459 = buffer.data(hl + 459);
    const auto *hl_460 = buffer.data(hl + 460);
    const auto *hl_461 = buffer.data(hl + 461);
    const auto *hl_462 = buffer.data(hl + 462);
    const auto *hl_463 = buffer.data(hl + 463);
    const auto *hl_464 = buffer.data(hl + 464);
    const auto *hl_465 = buffer.data(hl + 465);
    const auto *hl_466 = buffer.data(hl + 466);
    const auto *hl_467 = buffer.data(hl + 467);
    const auto *hl_468 = buffer.data(hl + 468);
    const auto *hl_469 = buffer.data(hl + 469);
    const auto *hl_470 = buffer.data(hl + 470);
    const auto *hl_471 = buffer.data(hl + 471);
    const auto *hl_472 = buffer.data(hl + 472);
    const auto *hl_473 = buffer.data(hl + 473);
    const auto *hl_474 = buffer.data(hl + 474);
    const auto *hl_475 = buffer.data(hl + 475);
    const auto *hl_476 = buffer.data(hl + 476);
    const auto *hl_477 = buffer.data(hl + 477);
    const auto *hl_478 = buffer.data(hl + 478);
    const auto *hl_479 = buffer.data(hl + 479);
    const auto *hl_480 = buffer.data(hl + 480);
    const auto *hl_481 = buffer.data(hl + 481);
    const auto *hl_482 = buffer.data(hl + 482);
    const auto *hl_483 = buffer.data(hl + 483);
    const auto *hl_484 = buffer.data(hl + 484);
    const auto *hl_485 = buffer.data(hl + 485);
    const auto *hl_486 = buffer.data(hl + 486);
    const auto *hl_487 = buffer.data(hl + 487);
    const auto *hl_488 = buffer.data(hl + 488);
    const auto *hl_489 = buffer.data(hl + 489);
    const auto *hl_490 = buffer.data(hl + 490);
    const auto *hl_491 = buffer.data(hl + 491);
    const auto *hl_492 = buffer.data(hl + 492);
    const auto *hl_493 = buffer.data(hl + 493);
    const auto *hl_494 = buffer.data(hl + 494);
    const auto *hl_495 = buffer.data(hl + 495);
    const auto *hl_496 = buffer.data(hl + 496);
    const auto *hl_497 = buffer.data(hl + 497);
    const auto *hl_498 = buffer.data(hl + 498);
    const auto *hl_499 = buffer.data(hl + 499);
    const auto *hl_500 = buffer.data(hl + 500);
    const auto *hl_501 = buffer.data(hl + 501);
    const auto *hl_502 = buffer.data(hl + 502);
    const auto *hl_503 = buffer.data(hl + 503);
    const auto *hl_504 = buffer.data(hl + 504);
    const auto *hl_505 = buffer.data(hl + 505);
    const auto *hl_506 = buffer.data(hl + 506);
    const auto *hl_507 = buffer.data(hl + 507);
    const auto *hl_508 = buffer.data(hl + 508);
    const auto *hl_509 = buffer.data(hl + 509);
    const auto *hl_510 = buffer.data(hl + 510);
    const auto *hl_511 = buffer.data(hl + 511);
    const auto *hl_512 = buffer.data(hl + 512);
    const auto *hl_513 = buffer.data(hl + 513);
    const auto *hl_514 = buffer.data(hl + 514);
    const auto *hl_515 = buffer.data(hl + 515);
    const auto *hl_516 = buffer.data(hl + 516);
    const auto *hl_517 = buffer.data(hl + 517);
    const auto *hl_518 = buffer.data(hl + 518);
    const auto *hl_519 = buffer.data(hl + 519);
    const auto *hl_520 = buffer.data(hl + 520);
    const auto *hl_521 = buffer.data(hl + 521);
    const auto *hl_522 = buffer.data(hl + 522);
    const auto *hl_523 = buffer.data(hl + 523);
    const auto *hl_524 = buffer.data(hl + 524);
    const auto *hl_525 = buffer.data(hl + 525);
    const auto *hl_526 = buffer.data(hl + 526);
    const auto *hl_527 = buffer.data(hl + 527);
    const auto *hl_528 = buffer.data(hl + 528);
    const auto *hl_529 = buffer.data(hl + 529);
    const auto *hl_530 = buffer.data(hl + 530);
    const auto *hl_531 = buffer.data(hl + 531);
    const auto *hl_532 = buffer.data(hl + 532);
    const auto *hl_533 = buffer.data(hl + 533);
    const auto *hl_534 = buffer.data(hl + 534);
    const auto *hl_535 = buffer.data(hl + 535);
    const auto *hl_536 = buffer.data(hl + 536);
    const auto *hl_537 = buffer.data(hl + 537);
    const auto *hl_538 = buffer.data(hl + 538);
    const auto *hl_539 = buffer.data(hl + 539);
    const auto *hl_540 = buffer.data(hl + 540);
    const auto *hl_541 = buffer.data(hl + 541);
    const auto *hl_542 = buffer.data(hl + 542);
    const auto *hl_543 = buffer.data(hl + 543);
    const auto *hl_544 = buffer.data(hl + 544);
    const auto *hl_545 = buffer.data(hl + 545);
    const auto *hl_546 = buffer.data(hl + 546);
    const auto *hl_547 = buffer.data(hl + 547);
    const auto *hl_548 = buffer.data(hl + 548);
    const auto *hl_549 = buffer.data(hl + 549);
    const auto *hl_550 = buffer.data(hl + 550);
    const auto *hl_551 = buffer.data(hl + 551);
    const auto *hl_552 = buffer.data(hl + 552);
    const auto *hl_553 = buffer.data(hl + 553);
    const auto *hl_554 = buffer.data(hl + 554);
    const auto *hl_555 = buffer.data(hl + 555);
    const auto *hl_556 = buffer.data(hl + 556);
    const auto *hl_557 = buffer.data(hl + 557);
    const auto *hl_558 = buffer.data(hl + 558);
    const auto *hl_559 = buffer.data(hl + 559);
    const auto *hl_560 = buffer.data(hl + 560);
    const auto *hl_561 = buffer.data(hl + 561);
    const auto *hl_562 = buffer.data(hl + 562);
    const auto *hl_563 = buffer.data(hl + 563);
    const auto *hl_564 = buffer.data(hl + 564);
    const auto *hl_565 = buffer.data(hl + 565);
    const auto *hl_566 = buffer.data(hl + 566);
    const auto *hl_567 = buffer.data(hl + 567);
    const auto *hl_568 = buffer.data(hl + 568);
    const auto *hl_569 = buffer.data(hl + 569);
    const auto *hl_570 = buffer.data(hl + 570);
    const auto *hl_571 = buffer.data(hl + 571);
    const auto *hl_572 = buffer.data(hl + 572);
    const auto *hl_573 = buffer.data(hl + 573);
    const auto *hl_574 = buffer.data(hl + 574);
    const auto *hl_575 = buffer.data(hl + 575);
    const auto *hl_576 = buffer.data(hl + 576);
    const auto *hl_577 = buffer.data(hl + 577);
    const auto *hl_578 = buffer.data(hl + 578);
    const auto *hl_579 = buffer.data(hl + 579);
    const auto *hl_580 = buffer.data(hl + 580);
    const auto *hl_581 = buffer.data(hl + 581);
    const auto *hl_582 = buffer.data(hl + 582);
    const auto *hl_583 = buffer.data(hl + 583);
    const auto *hl_584 = buffer.data(hl + 584);
    const auto *hl_585 = buffer.data(hl + 585);
    const auto *hl_586 = buffer.data(hl + 586);
    const auto *hl_587 = buffer.data(hl + 587);
    const auto *hl_588 = buffer.data(hl + 588);
    const auto *hl_589 = buffer.data(hl + 589);
    const auto *hl_590 = buffer.data(hl + 590);
    const auto *hl_591 = buffer.data(hl + 591);
    const auto *hl_592 = buffer.data(hl + 592);
    const auto *hl_593 = buffer.data(hl + 593);
    const auto *hl_594 = buffer.data(hl + 594);
    const auto *hl_595 = buffer.data(hl + 595);
    const auto *hl_596 = buffer.data(hl + 596);
    const auto *hl_597 = buffer.data(hl + 597);
    const auto *hl_598 = buffer.data(hl + 598);
    const auto *hl_599 = buffer.data(hl + 599);
    const auto *hl_600 = buffer.data(hl + 600);
    const auto *hl_601 = buffer.data(hl + 601);
    const auto *hl_602 = buffer.data(hl + 602);
    const auto *hl_603 = buffer.data(hl + 603);
    const auto *hl_604 = buffer.data(hl + 604);
    const auto *hl_605 = buffer.data(hl + 605);
    const auto *hl_606 = buffer.data(hl + 606);
    const auto *hl_607 = buffer.data(hl + 607);
    const auto *hl_608 = buffer.data(hl + 608);
    const auto *hl_609 = buffer.data(hl + 609);
    const auto *hl_610 = buffer.data(hl + 610);
    const auto *hl_611 = buffer.data(hl + 611);
    const auto *hl_612 = buffer.data(hl + 612);
    const auto *hl_613 = buffer.data(hl + 613);
    const auto *hl_614 = buffer.data(hl + 614);
    const auto *hl_615 = buffer.data(hl + 615);
    const auto *hl_616 = buffer.data(hl + 616);
    const auto *hl_617 = buffer.data(hl + 617);
    const auto *hl_618 = buffer.data(hl + 618);
    const auto *hl_619 = buffer.data(hl + 619);
    const auto *hl_620 = buffer.data(hl + 620);
    const auto *hl_621 = buffer.data(hl + 621);
    const auto *hl_622 = buffer.data(hl + 622);
    const auto *hl_623 = buffer.data(hl + 623);
    const auto *hl_624 = buffer.data(hl + 624);
    const auto *hl_625 = buffer.data(hl + 625);
    const auto *hl_626 = buffer.data(hl + 626);
    const auto *hl_627 = buffer.data(hl + 627);
    const auto *hl_628 = buffer.data(hl + 628);
    const auto *hl_629 = buffer.data(hl + 629);
    const auto *hl_630 = buffer.data(hl + 630);
    const auto *hl_631 = buffer.data(hl + 631);
    const auto *hl_632 = buffer.data(hl + 632);
    const auto *hl_633 = buffer.data(hl + 633);
    const auto *hl_634 = buffer.data(hl + 634);
    const auto *hl_635 = buffer.data(hl + 635);
    const auto *hl_636 = buffer.data(hl + 636);
    const auto *hl_637 = buffer.data(hl + 637);
    const auto *hl_638 = buffer.data(hl + 638);
    const auto *hl_639 = buffer.data(hl + 639);
    const auto *hl_640 = buffer.data(hl + 640);
    const auto *hl_641 = buffer.data(hl + 641);
    const auto *hl_642 = buffer.data(hl + 642);
    const auto *hl_643 = buffer.data(hl + 643);
    const auto *hl_644 = buffer.data(hl + 644);
    const auto *hl_645 = buffer.data(hl + 645);
    const auto *hl_646 = buffer.data(hl + 646);
    const auto *hl_647 = buffer.data(hl + 647);
    const auto *hl_648 = buffer.data(hl + 648);
    const auto *hl_649 = buffer.data(hl + 649);
    const auto *hl_650 = buffer.data(hl + 650);
    const auto *hl_651 = buffer.data(hl + 651);
    const auto *hl_652 = buffer.data(hl + 652);
    const auto *hl_653 = buffer.data(hl + 653);
    const auto *hl_654 = buffer.data(hl + 654);
    const auto *hl_655 = buffer.data(hl + 655);
    const auto *hl_656 = buffer.data(hl + 656);
    const auto *hl_657 = buffer.data(hl + 657);
    const auto *hl_658 = buffer.data(hl + 658);
    const auto *hl_659 = buffer.data(hl + 659);
    const auto *hl_660 = buffer.data(hl + 660);
    const auto *hl_661 = buffer.data(hl + 661);
    const auto *hl_662 = buffer.data(hl + 662);
    const auto *hl_663 = buffer.data(hl + 663);
    const auto *hl_664 = buffer.data(hl + 664);
    const auto *hl_665 = buffer.data(hl + 665);
    const auto *hl_666 = buffer.data(hl + 666);
    const auto *hl_667 = buffer.data(hl + 667);
    const auto *hl_668 = buffer.data(hl + 668);
    const auto *hl_669 = buffer.data(hl + 669);
    const auto *hl_670 = buffer.data(hl + 670);
    const auto *hl_671 = buffer.data(hl + 671);
    const auto *hl_672 = buffer.data(hl + 672);
    const auto *hl_673 = buffer.data(hl + 673);
    const auto *hl_674 = buffer.data(hl + 674);

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, t_455, t_456, t_457, hl_450, \
                         hl_451, hl_452, hl_453, hl_454, hl_455, hl_456, \
                         hl_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_0 * hl_450[k];

        t_451[k] = f_0 * hl_451[k];

        t_452[k] = f_0 * hl_452[k];

        t_453[k] = f_0 * hl_453[k];

        t_454[k] = f_0 * hl_454[k];

        t_455[k] = f_0 * hl_455[k];

        t_456[k] = f_0 * hl_456[k];

        t_457[k] = f_0 * hl_457[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, t_462, t_463, t_464, t_465, hl_458, \
                         hl_459, hl_460, hl_461, hl_462, hl_463, hl_464, \
                         hl_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_0 * hl_458[k];

        t_459[k] = f_0 * hl_459[k];

        t_460[k] = f_0 * hl_460[k];

        t_461[k] = f_0 * hl_461[k];

        t_462[k] = f_0 * hl_462[k];

        t_463[k] = f_0 * hl_463[k];

        t_464[k] = f_0 * hl_464[k];

        t_465[k] = f_0 * hl_465[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, t_471, t_472, t_473, hl_466, \
                         hl_467, hl_468, hl_469, hl_470, hl_471, hl_472, \
                         hl_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_0 * hl_466[k];

        t_467[k] = f_0 * hl_467[k];

        t_468[k] = f_0 * hl_468[k];

        t_469[k] = f_0 * hl_469[k];

        t_470[k] = f_0 * hl_470[k];

        t_471[k] = f_0 * hl_471[k];

        t_472[k] = f_0 * hl_472[k];

        t_473[k] = f_0 * hl_473[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, t_478, t_479, t_480, t_481, hl_474, \
                         hl_475, hl_476, hl_477, hl_478, hl_479, hl_480, \
                         hl_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_0 * hl_474[k];

        t_475[k] = f_0 * hl_475[k];

        t_476[k] = f_0 * hl_476[k];

        t_477[k] = f_0 * hl_477[k];

        t_478[k] = f_0 * hl_478[k];

        t_479[k] = f_0 * hl_479[k];

        t_480[k] = f_0 * hl_480[k];

        t_481[k] = f_0 * hl_481[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, t_486, t_487, t_488, t_489, hl_482, \
                         hl_483, hl_484, hl_485, hl_486, hl_487, hl_488, \
                         hl_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_0 * hl_482[k];

        t_483[k] = f_0 * hl_483[k];

        t_484[k] = f_0 * hl_484[k];

        t_485[k] = f_0 * hl_485[k];

        t_486[k] = f_0 * hl_486[k];

        t_487[k] = f_0 * hl_487[k];

        t_488[k] = f_0 * hl_488[k];

        t_489[k] = f_0 * hl_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, t_495, t_496, t_497, hl_490, \
                         hl_491, hl_492, hl_493, hl_494, hl_495, hl_496, \
                         hl_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = f_0 * hl_490[k];

        t_491[k] = f_0 * hl_491[k];

        t_492[k] = f_0 * hl_492[k];

        t_493[k] = f_0 * hl_493[k];

        t_494[k] = f_0 * hl_494[k];

        t_495[k] = f_0 * hl_495[k];

        t_496[k] = f_0 * hl_496[k];

        t_497[k] = f_0 * hl_497[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, t_501, t_502, t_503, t_504, t_505, hl_498, \
                         hl_499, hl_500, hl_501, hl_502, hl_503, hl_504, \
                         hl_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_0 * hl_498[k];

        t_499[k] = f_0 * hl_499[k];

        t_500[k] = f_0 * hl_500[k];

        t_501[k] = f_0 * hl_501[k];

        t_502[k] = f_0 * hl_502[k];

        t_503[k] = f_0 * hl_503[k];

        t_504[k] = f_0 * hl_504[k];

        t_505[k] = f_0 * hl_505[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, t_509, t_510, t_511, t_512, t_513, hl_506, \
                         hl_507, hl_508, hl_509, hl_510, hl_511, hl_512, \
                         hl_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_0 * hl_506[k];

        t_507[k] = f_0 * hl_507[k];

        t_508[k] = f_0 * hl_508[k];

        t_509[k] = f_0 * hl_509[k];

        t_510[k] = f_0 * hl_510[k];

        t_511[k] = f_0 * hl_511[k];

        t_512[k] = f_0 * hl_512[k];

        t_513[k] = f_0 * hl_513[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, t_519, t_520, t_521, hl_514, \
                         hl_515, hl_516, hl_517, hl_518, hl_519, hl_520, \
                         hl_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_0 * hl_514[k];

        t_515[k] = f_0 * hl_515[k];

        t_516[k] = f_0 * hl_516[k];

        t_517[k] = f_0 * hl_517[k];

        t_518[k] = f_0 * hl_518[k];

        t_519[k] = f_0 * hl_519[k];

        t_520[k] = f_0 * hl_520[k];

        t_521[k] = f_0 * hl_521[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, t_527, t_528, t_529, hl_522, \
                         hl_523, hl_524, hl_525, hl_526, hl_527, hl_528, \
                         hl_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_0 * hl_522[k];

        t_523[k] = f_0 * hl_523[k];

        t_524[k] = f_0 * hl_524[k];

        t_525[k] = f_0 * hl_525[k];

        t_526[k] = f_0 * hl_526[k];

        t_527[k] = f_0 * hl_527[k];

        t_528[k] = f_0 * hl_528[k];

        t_529[k] = f_0 * hl_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, t_535, t_536, t_537, hl_530, \
                         hl_531, hl_532, hl_533, hl_534, hl_535, hl_536, \
                         hl_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_0 * hl_530[k];

        t_531[k] = f_0 * hl_531[k];

        t_532[k] = f_0 * hl_532[k];

        t_533[k] = f_0 * hl_533[k];

        t_534[k] = f_0 * hl_534[k];

        t_535[k] = f_0 * hl_535[k];

        t_536[k] = f_0 * hl_536[k];

        t_537[k] = f_0 * hl_537[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, t_542, t_543, t_544, t_545, hl_538, \
                         hl_539, hl_540, hl_541, hl_542, hl_543, hl_544, \
                         hl_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_0 * hl_538[k];

        t_539[k] = f_0 * hl_539[k];

        t_540[k] = f_0 * hl_540[k];

        t_541[k] = f_0 * hl_541[k];

        t_542[k] = f_0 * hl_542[k];

        t_543[k] = f_0 * hl_543[k];

        t_544[k] = f_0 * hl_544[k];

        t_545[k] = f_0 * hl_545[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, t_550, t_551, t_552, t_553, hl_546, \
                         hl_547, hl_548, hl_549, hl_550, hl_551, hl_552, \
                         hl_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = f_0 * hl_546[k];

        t_547[k] = f_0 * hl_547[k];

        t_548[k] = f_0 * hl_548[k];

        t_549[k] = f_0 * hl_549[k];

        t_550[k] = f_0 * hl_550[k];

        t_551[k] = f_0 * hl_551[k];

        t_552[k] = f_0 * hl_552[k];

        t_553[k] = f_0 * hl_553[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, t_558, t_559, t_560, t_561, hl_554, \
                         hl_555, hl_556, hl_557, hl_558, hl_559, hl_560, \
                         hl_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_0 * hl_554[k];

        t_555[k] = f_0 * hl_555[k];

        t_556[k] = f_0 * hl_556[k];

        t_557[k] = f_0 * hl_557[k];

        t_558[k] = f_0 * hl_558[k];

        t_559[k] = f_0 * hl_559[k];

        t_560[k] = f_0 * hl_560[k];

        t_561[k] = f_0 * hl_561[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, t_567, t_568, t_569, hl_562, \
                         hl_563, hl_564, hl_565, hl_566, hl_567, hl_568, \
                         hl_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_0 * hl_562[k];

        t_563[k] = f_0 * hl_563[k];

        t_564[k] = f_0 * hl_564[k];

        t_565[k] = f_0 * hl_565[k];

        t_566[k] = f_0 * hl_566[k];

        t_567[k] = f_0 * hl_567[k];

        t_568[k] = f_0 * hl_568[k];

        t_569[k] = f_0 * hl_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, t_575, t_576, t_577, hl_570, \
                         hl_571, hl_572, hl_573, hl_574, hl_575, hl_576, \
                         hl_577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_0 * hl_570[k];

        t_571[k] = f_0 * hl_571[k];

        t_572[k] = f_0 * hl_572[k];

        t_573[k] = f_0 * hl_573[k];

        t_574[k] = f_0 * hl_574[k];

        t_575[k] = f_0 * hl_575[k];

        t_576[k] = f_0 * hl_576[k];

        t_577[k] = f_0 * hl_577[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, t_582, t_583, t_584, t_585, hl_578, \
                         hl_579, hl_580, hl_581, hl_582, hl_583, hl_584, \
                         hl_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_0 * hl_578[k];

        t_579[k] = f_0 * hl_579[k];

        t_580[k] = f_0 * hl_580[k];

        t_581[k] = f_0 * hl_581[k];

        t_582[k] = f_0 * hl_582[k];

        t_583[k] = f_0 * hl_583[k];

        t_584[k] = f_0 * hl_584[k];

        t_585[k] = f_0 * hl_585[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, t_590, t_591, t_592, t_593, hl_586, \
                         hl_587, hl_588, hl_589, hl_590, hl_591, hl_592, \
                         hl_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_0 * hl_586[k];

        t_587[k] = f_0 * hl_587[k];

        t_588[k] = f_0 * hl_588[k];

        t_589[k] = f_0 * hl_589[k];

        t_590[k] = f_0 * hl_590[k];

        t_591[k] = f_0 * hl_591[k];

        t_592[k] = f_0 * hl_592[k];

        t_593[k] = f_0 * hl_593[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, t_598, t_599, t_600, t_601, hl_594, \
                         hl_595, hl_596, hl_597, hl_598, hl_599, hl_600, \
                         hl_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_0 * hl_594[k];

        t_595[k] = f_0 * hl_595[k];

        t_596[k] = f_0 * hl_596[k];

        t_597[k] = f_0 * hl_597[k];

        t_598[k] = f_0 * hl_598[k];

        t_599[k] = f_0 * hl_599[k];

        t_600[k] = f_0 * hl_600[k];

        t_601[k] = f_0 * hl_601[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, t_607, t_608, t_609, hl_602, \
                         hl_603, hl_604, hl_605, hl_606, hl_607, hl_608, \
                         hl_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_0 * hl_602[k];

        t_603[k] = f_0 * hl_603[k];

        t_604[k] = f_0 * hl_604[k];

        t_605[k] = f_0 * hl_605[k];

        t_606[k] = f_0 * hl_606[k];

        t_607[k] = f_0 * hl_607[k];

        t_608[k] = f_0 * hl_608[k];

        t_609[k] = f_0 * hl_609[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, t_615, t_616, t_617, hl_610, \
                         hl_611, hl_612, hl_613, hl_614, hl_615, hl_616, \
                         hl_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_0 * hl_610[k];

        t_611[k] = f_0 * hl_611[k];

        t_612[k] = f_0 * hl_612[k];

        t_613[k] = f_0 * hl_613[k];

        t_614[k] = f_0 * hl_614[k];

        t_615[k] = f_0 * hl_615[k];

        t_616[k] = f_0 * hl_616[k];

        t_617[k] = f_0 * hl_617[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, t_622, t_623, t_624, t_625, hl_618, \
                         hl_619, hl_620, hl_621, hl_622, hl_623, hl_624, \
                         hl_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_0 * hl_618[k];

        t_619[k] = f_0 * hl_619[k];

        t_620[k] = f_0 * hl_620[k];

        t_621[k] = f_0 * hl_621[k];

        t_622[k] = f_0 * hl_622[k];

        t_623[k] = f_0 * hl_623[k];

        t_624[k] = f_0 * hl_624[k];

        t_625[k] = f_0 * hl_625[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, t_630, t_631, t_632, t_633, hl_626, \
                         hl_627, hl_628, hl_629, hl_630, hl_631, hl_632, \
                         hl_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_0 * hl_626[k];

        t_627[k] = f_0 * hl_627[k];

        t_628[k] = f_0 * hl_628[k];

        t_629[k] = f_0 * hl_629[k];

        t_630[k] = f_0 * hl_630[k];

        t_631[k] = f_0 * hl_631[k];

        t_632[k] = f_0 * hl_632[k];

        t_633[k] = f_0 * hl_633[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, t_639, t_640, t_641, hl_634, \
                         hl_635, hl_636, hl_637, hl_638, hl_639, hl_640, \
                         hl_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_0 * hl_634[k];

        t_635[k] = f_0 * hl_635[k];

        t_636[k] = f_0 * hl_636[k];

        t_637[k] = f_0 * hl_637[k];

        t_638[k] = f_0 * hl_638[k];

        t_639[k] = f_0 * hl_639[k];

        t_640[k] = f_0 * hl_640[k];

        t_641[k] = f_0 * hl_641[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, t_646, t_647, t_648, t_649, hl_642, \
                         hl_643, hl_644, hl_645, hl_646, hl_647, hl_648, \
                         hl_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_0 * hl_642[k];

        t_643[k] = f_0 * hl_643[k];

        t_644[k] = f_0 * hl_644[k];

        t_645[k] = f_0 * hl_645[k];

        t_646[k] = f_0 * hl_646[k];

        t_647[k] = f_0 * hl_647[k];

        t_648[k] = f_0 * hl_648[k];

        t_649[k] = f_0 * hl_649[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, t_655, t_656, t_657, hl_650, \
                         hl_651, hl_652, hl_653, hl_654, hl_655, hl_656, \
                         hl_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_0 * hl_650[k];

        t_651[k] = f_0 * hl_651[k];

        t_652[k] = f_0 * hl_652[k];

        t_653[k] = f_0 * hl_653[k];

        t_654[k] = f_0 * hl_654[k];

        t_655[k] = f_0 * hl_655[k];

        t_656[k] = f_0 * hl_656[k];

        t_657[k] = f_0 * hl_657[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, t_662, t_663, t_664, t_665, hl_658, \
                         hl_659, hl_660, hl_661, hl_662, hl_663, hl_664, \
                         hl_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = f_0 * hl_658[k];

        t_659[k] = f_0 * hl_659[k];

        t_660[k] = f_0 * hl_660[k];

        t_661[k] = f_0 * hl_661[k];

        t_662[k] = f_0 * hl_662[k];

        t_663[k] = f_0 * hl_663[k];

        t_664[k] = f_0 * hl_664[k];

        t_665[k] = f_0 * hl_665[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, t_670, t_671, t_672, t_673, hl_666, \
                         hl_667, hl_668, hl_669, hl_670, hl_671, hl_672, \
                         hl_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_0 * hl_666[k];

        t_667[k] = f_0 * hl_667[k];

        t_668[k] = f_0 * hl_668[k];

        t_669[k] = f_0 * hl_669[k];

        t_670[k] = f_0 * hl_670[k];

        t_671[k] = f_0 * hl_671[k];

        t_672[k] = f_0 * hl_672[k];

        t_673[k] = f_0 * hl_673[k];
    }

#pragma omp simd aligned(t_674, hl_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_0 * hl_674[k];
    }
}

auto
compute_prim_geom_10_gl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t fl, const size_t hl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gl_electron_repulsion_0_piece0(buffer, target, fl, hl, ncols, alpha);

    compute_prim_geom_10_gl_electron_repulsion_0_piece1(buffer, target, fl, hl, ncols, alpha);

    compute_prim_geom_10_gl_electron_repulsion_0_piece2(buffer, target, fl, hl, ncols, alpha);

    compute_prim_geom_10_gl_electron_repulsion_0_piece3(buffer, target, hl, ncols, alpha);
}

static auto
compute_prim_geom_10_gl_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
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

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);
    const auto *fl_9 = buffer.data(fl + 9);
    const auto *fl_10 = buffer.data(fl + 10);
    const auto *fl_11 = buffer.data(fl + 11);
    const auto *fl_12 = buffer.data(fl + 12);
    const auto *fl_13 = buffer.data(fl + 13);
    const auto *fl_14 = buffer.data(fl + 14);
    const auto *fl_15 = buffer.data(fl + 15);
    const auto *fl_16 = buffer.data(fl + 16);
    const auto *fl_17 = buffer.data(fl + 17);
    const auto *fl_18 = buffer.data(fl + 18);
    const auto *fl_19 = buffer.data(fl + 19);
    const auto *fl_20 = buffer.data(fl + 20);
    const auto *fl_21 = buffer.data(fl + 21);
    const auto *fl_22 = buffer.data(fl + 22);
    const auto *fl_23 = buffer.data(fl + 23);
    const auto *fl_24 = buffer.data(fl + 24);
    const auto *fl_25 = buffer.data(fl + 25);
    const auto *fl_26 = buffer.data(fl + 26);
    const auto *fl_27 = buffer.data(fl + 27);
    const auto *fl_28 = buffer.data(fl + 28);
    const auto *fl_29 = buffer.data(fl + 29);
    const auto *fl_30 = buffer.data(fl + 30);
    const auto *fl_31 = buffer.data(fl + 31);
    const auto *fl_32 = buffer.data(fl + 32);
    const auto *fl_33 = buffer.data(fl + 33);
    const auto *fl_34 = buffer.data(fl + 34);
    const auto *fl_35 = buffer.data(fl + 35);
    const auto *fl_36 = buffer.data(fl + 36);
    const auto *fl_37 = buffer.data(fl + 37);
    const auto *fl_38 = buffer.data(fl + 38);
    const auto *fl_39 = buffer.data(fl + 39);
    const auto *fl_40 = buffer.data(fl + 40);
    const auto *fl_41 = buffer.data(fl + 41);
    const auto *fl_42 = buffer.data(fl + 42);
    const auto *fl_43 = buffer.data(fl + 43);
    const auto *fl_44 = buffer.data(fl + 44);
    const auto *fl_45 = buffer.data(fl + 45);
    const auto *fl_46 = buffer.data(fl + 46);
    const auto *fl_47 = buffer.data(fl + 47);
    const auto *fl_48 = buffer.data(fl + 48);
    const auto *fl_49 = buffer.data(fl + 49);
    const auto *fl_50 = buffer.data(fl + 50);
    const auto *fl_51 = buffer.data(fl + 51);
    const auto *fl_52 = buffer.data(fl + 52);
    const auto *fl_53 = buffer.data(fl + 53);
    const auto *fl_54 = buffer.data(fl + 54);
    const auto *fl_55 = buffer.data(fl + 55);
    const auto *fl_56 = buffer.data(fl + 56);
    const auto *fl_57 = buffer.data(fl + 57);
    const auto *fl_58 = buffer.data(fl + 58);
    const auto *fl_59 = buffer.data(fl + 59);
    const auto *fl_60 = buffer.data(fl + 60);
    const auto *fl_61 = buffer.data(fl + 61);
    const auto *fl_62 = buffer.data(fl + 62);
    const auto *fl_63 = buffer.data(fl + 63);
    const auto *fl_64 = buffer.data(fl + 64);
    const auto *fl_65 = buffer.data(fl + 65);
    const auto *fl_66 = buffer.data(fl + 66);
    const auto *fl_67 = buffer.data(fl + 67);
    const auto *fl_68 = buffer.data(fl + 68);
    const auto *fl_69 = buffer.data(fl + 69);
    const auto *fl_70 = buffer.data(fl + 70);
    const auto *fl_71 = buffer.data(fl + 71);
    const auto *fl_72 = buffer.data(fl + 72);
    const auto *fl_73 = buffer.data(fl + 73);
    const auto *fl_74 = buffer.data(fl + 74);
    const auto *fl_75 = buffer.data(fl + 75);
    const auto *fl_76 = buffer.data(fl + 76);
    const auto *fl_77 = buffer.data(fl + 77);
    const auto *fl_78 = buffer.data(fl + 78);
    const auto *fl_79 = buffer.data(fl + 79);
    const auto *fl_80 = buffer.data(fl + 80);
    const auto *fl_81 = buffer.data(fl + 81);
    const auto *fl_82 = buffer.data(fl + 82);
    const auto *fl_83 = buffer.data(fl + 83);
    const auto *fl_84 = buffer.data(fl + 84);
    const auto *fl_85 = buffer.data(fl + 85);
    const auto *fl_86 = buffer.data(fl + 86);
    const auto *fl_87 = buffer.data(fl + 87);
    const auto *fl_88 = buffer.data(fl + 88);

    const auto *hl_45 = buffer.data(hl + 45);
    const auto *hl_46 = buffer.data(hl + 46);
    const auto *hl_47 = buffer.data(hl + 47);
    const auto *hl_48 = buffer.data(hl + 48);
    const auto *hl_49 = buffer.data(hl + 49);
    const auto *hl_50 = buffer.data(hl + 50);
    const auto *hl_51 = buffer.data(hl + 51);
    const auto *hl_52 = buffer.data(hl + 52);
    const auto *hl_53 = buffer.data(hl + 53);
    const auto *hl_54 = buffer.data(hl + 54);
    const auto *hl_55 = buffer.data(hl + 55);
    const auto *hl_56 = buffer.data(hl + 56);
    const auto *hl_57 = buffer.data(hl + 57);
    const auto *hl_58 = buffer.data(hl + 58);
    const auto *hl_59 = buffer.data(hl + 59);
    const auto *hl_60 = buffer.data(hl + 60);
    const auto *hl_61 = buffer.data(hl + 61);
    const auto *hl_62 = buffer.data(hl + 62);
    const auto *hl_63 = buffer.data(hl + 63);
    const auto *hl_64 = buffer.data(hl + 64);
    const auto *hl_65 = buffer.data(hl + 65);
    const auto *hl_66 = buffer.data(hl + 66);
    const auto *hl_67 = buffer.data(hl + 67);
    const auto *hl_68 = buffer.data(hl + 68);
    const auto *hl_69 = buffer.data(hl + 69);
    const auto *hl_70 = buffer.data(hl + 70);
    const auto *hl_71 = buffer.data(hl + 71);
    const auto *hl_72 = buffer.data(hl + 72);
    const auto *hl_73 = buffer.data(hl + 73);
    const auto *hl_74 = buffer.data(hl + 74);
    const auto *hl_75 = buffer.data(hl + 75);
    const auto *hl_76 = buffer.data(hl + 76);
    const auto *hl_77 = buffer.data(hl + 77);
    const auto *hl_78 = buffer.data(hl + 78);
    const auto *hl_79 = buffer.data(hl + 79);
    const auto *hl_80 = buffer.data(hl + 80);
    const auto *hl_81 = buffer.data(hl + 81);
    const auto *hl_82 = buffer.data(hl + 82);
    const auto *hl_83 = buffer.data(hl + 83);
    const auto *hl_84 = buffer.data(hl + 84);
    const auto *hl_85 = buffer.data(hl + 85);
    const auto *hl_86 = buffer.data(hl + 86);
    const auto *hl_87 = buffer.data(hl + 87);
    const auto *hl_88 = buffer.data(hl + 88);
    const auto *hl_89 = buffer.data(hl + 89);
    const auto *hl_135 = buffer.data(hl + 135);
    const auto *hl_136 = buffer.data(hl + 136);
    const auto *hl_137 = buffer.data(hl + 137);
    const auto *hl_138 = buffer.data(hl + 138);
    const auto *hl_139 = buffer.data(hl + 139);
    const auto *hl_140 = buffer.data(hl + 140);
    const auto *hl_141 = buffer.data(hl + 141);
    const auto *hl_142 = buffer.data(hl + 142);
    const auto *hl_143 = buffer.data(hl + 143);
    const auto *hl_144 = buffer.data(hl + 144);
    const auto *hl_145 = buffer.data(hl + 145);
    const auto *hl_146 = buffer.data(hl + 146);
    const auto *hl_147 = buffer.data(hl + 147);
    const auto *hl_148 = buffer.data(hl + 148);
    const auto *hl_149 = buffer.data(hl + 149);
    const auto *hl_150 = buffer.data(hl + 150);
    const auto *hl_151 = buffer.data(hl + 151);
    const auto *hl_152 = buffer.data(hl + 152);
    const auto *hl_153 = buffer.data(hl + 153);
    const auto *hl_154 = buffer.data(hl + 154);
    const auto *hl_155 = buffer.data(hl + 155);
    const auto *hl_156 = buffer.data(hl + 156);
    const auto *hl_157 = buffer.data(hl + 157);
    const auto *hl_158 = buffer.data(hl + 158);
    const auto *hl_159 = buffer.data(hl + 159);
    const auto *hl_160 = buffer.data(hl + 160);
    const auto *hl_161 = buffer.data(hl + 161);
    const auto *hl_162 = buffer.data(hl + 162);
    const auto *hl_163 = buffer.data(hl + 163);
    const auto *hl_164 = buffer.data(hl + 164);
    const auto *hl_165 = buffer.data(hl + 165);
    const auto *hl_166 = buffer.data(hl + 166);
    const auto *hl_167 = buffer.data(hl + 167);
    const auto *hl_168 = buffer.data(hl + 168);
    const auto *hl_169 = buffer.data(hl + 169);
    const auto *hl_170 = buffer.data(hl + 170);
    const auto *hl_171 = buffer.data(hl + 171);
    const auto *hl_172 = buffer.data(hl + 172);
    const auto *hl_173 = buffer.data(hl + 173);
    const auto *hl_174 = buffer.data(hl + 174);
    const auto *hl_175 = buffer.data(hl + 175);
    const auto *hl_176 = buffer.data(hl + 176);
    const auto *hl_177 = buffer.data(hl + 177);
    const auto *hl_178 = buffer.data(hl + 178);
    const auto *hl_179 = buffer.data(hl + 179);
    const auto *hl_180 = buffer.data(hl + 180);
    const auto *hl_181 = buffer.data(hl + 181);
    const auto *hl_182 = buffer.data(hl + 182);
    const auto *hl_183 = buffer.data(hl + 183);
    const auto *hl_184 = buffer.data(hl + 184);
    const auto *hl_185 = buffer.data(hl + 185);
    const auto *hl_186 = buffer.data(hl + 186);
    const auto *hl_187 = buffer.data(hl + 187);
    const auto *hl_188 = buffer.data(hl + 188);
    const auto *hl_189 = buffer.data(hl + 189);
    const auto *hl_190 = buffer.data(hl + 190);
    const auto *hl_191 = buffer.data(hl + 191);
    const auto *hl_192 = buffer.data(hl + 192);
    const auto *hl_193 = buffer.data(hl + 193);
    const auto *hl_194 = buffer.data(hl + 194);
    const auto *hl_195 = buffer.data(hl + 195);
    const auto *hl_196 = buffer.data(hl + 196);
    const auto *hl_197 = buffer.data(hl + 197);
    const auto *hl_198 = buffer.data(hl + 198);
    const auto *hl_199 = buffer.data(hl + 199);
    const auto *hl_200 = buffer.data(hl + 200);
    const auto *hl_201 = buffer.data(hl + 201);
    const auto *hl_202 = buffer.data(hl + 202);
    const auto *hl_203 = buffer.data(hl + 203);
    const auto *hl_204 = buffer.data(hl + 204);
    const auto *hl_205 = buffer.data(hl + 205);
    const auto *hl_206 = buffer.data(hl + 206);
    const auto *hl_207 = buffer.data(hl + 207);
    const auto *hl_208 = buffer.data(hl + 208);
    const auto *hl_209 = buffer.data(hl + 209);
    const auto *hl_210 = buffer.data(hl + 210);
    const auto *hl_211 = buffer.data(hl + 211);
    const auto *hl_212 = buffer.data(hl + 212);
    const auto *hl_213 = buffer.data(hl + 213);
    const auto *hl_214 = buffer.data(hl + 214);
    const auto *hl_215 = buffer.data(hl + 215);
    const auto *hl_216 = buffer.data(hl + 216);
    const auto *hl_217 = buffer.data(hl + 217);
    const auto *hl_218 = buffer.data(hl + 218);
    const auto *hl_219 = buffer.data(hl + 219);
    const auto *hl_220 = buffer.data(hl + 220);
    const auto *hl_221 = buffer.data(hl + 221);
    const auto *hl_222 = buffer.data(hl + 222);
    const auto *hl_223 = buffer.data(hl + 223);
    const auto *hl_224 = buffer.data(hl + 224);
    const auto *hl_270 = buffer.data(hl + 270);
    const auto *hl_271 = buffer.data(hl + 271);
    const auto *hl_272 = buffer.data(hl + 272);
    const auto *hl_273 = buffer.data(hl + 273);
    const auto *hl_274 = buffer.data(hl + 274);
    const auto *hl_275 = buffer.data(hl + 275);
    const auto *hl_276 = buffer.data(hl + 276);
    const auto *hl_277 = buffer.data(hl + 277);
    const auto *hl_278 = buffer.data(hl + 278);
    const auto *hl_279 = buffer.data(hl + 279);
    const auto *hl_280 = buffer.data(hl + 280);
    const auto *hl_281 = buffer.data(hl + 281);
    const auto *hl_282 = buffer.data(hl + 282);
    const auto *hl_283 = buffer.data(hl + 283);
    const auto *hl_284 = buffer.data(hl + 284);
    const auto *hl_285 = buffer.data(hl + 285);
    const auto *hl_286 = buffer.data(hl + 286);
    const auto *hl_287 = buffer.data(hl + 287);
    const auto *hl_288 = buffer.data(hl + 288);
    const auto *hl_289 = buffer.data(hl + 289);
    const auto *hl_290 = buffer.data(hl + 290);
    const auto *hl_291 = buffer.data(hl + 291);
    const auto *hl_292 = buffer.data(hl + 292);
    const auto *hl_293 = buffer.data(hl + 293);
    const auto *hl_294 = buffer.data(hl + 294);
    const auto *hl_295 = buffer.data(hl + 295);
    const auto *hl_296 = buffer.data(hl + 296);
    const auto *hl_297 = buffer.data(hl + 297);
    const auto *hl_298 = buffer.data(hl + 298);
    const auto *hl_299 = buffer.data(hl + 299);
    const auto *hl_300 = buffer.data(hl + 300);
    const auto *hl_301 = buffer.data(hl + 301);
    const auto *hl_302 = buffer.data(hl + 302);
    const auto *hl_303 = buffer.data(hl + 303);
    const auto *hl_304 = buffer.data(hl + 304);
    const auto *hl_305 = buffer.data(hl + 305);
    const auto *hl_306 = buffer.data(hl + 306);
    const auto *hl_307 = buffer.data(hl + 307);
    const auto *hl_308 = buffer.data(hl + 308);
    const auto *hl_309 = buffer.data(hl + 309);
    const auto *hl_310 = buffer.data(hl + 310);
    const auto *hl_311 = buffer.data(hl + 311);
    const auto *hl_312 = buffer.data(hl + 312);
    const auto *hl_313 = buffer.data(hl + 313);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hl_45, hl_46, hl_47, hl_48, \
                         hl_49, hl_50, hl_51, hl_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hl_45[k];

        t_1[k] = f_0 * hl_46[k];

        t_2[k] = f_0 * hl_47[k];

        t_3[k] = f_0 * hl_48[k];

        t_4[k] = f_0 * hl_49[k];

        t_5[k] = f_0 * hl_50[k];

        t_6[k] = f_0 * hl_51[k];

        t_7[k] = f_0 * hl_52[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, hl_53, hl_54, hl_55, \
                         hl_56, hl_57, hl_58, hl_59, hl_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hl_53[k];

        t_9[k] = f_0 * hl_54[k];

        t_10[k] = f_0 * hl_55[k];

        t_11[k] = f_0 * hl_56[k];

        t_12[k] = f_0 * hl_57[k];

        t_13[k] = f_0 * hl_58[k];

        t_14[k] = f_0 * hl_59[k];

        t_15[k] = f_0 * hl_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, hl_61, hl_62, hl_63, \
                         hl_64, hl_65, hl_66, hl_67, hl_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * hl_61[k];

        t_17[k] = f_0 * hl_62[k];

        t_18[k] = f_0 * hl_63[k];

        t_19[k] = f_0 * hl_64[k];

        t_20[k] = f_0 * hl_65[k];

        t_21[k] = f_0 * hl_66[k];

        t_22[k] = f_0 * hl_67[k];

        t_23[k] = f_0 * hl_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, hl_69, hl_70, hl_71, \
                         hl_72, hl_73, hl_74, hl_75, hl_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * hl_69[k];

        t_25[k] = f_0 * hl_70[k];

        t_26[k] = f_0 * hl_71[k];

        t_27[k] = f_0 * hl_72[k];

        t_28[k] = f_0 * hl_73[k];

        t_29[k] = f_0 * hl_74[k];

        t_30[k] = f_0 * hl_75[k];

        t_31[k] = f_0 * hl_76[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, hl_77, hl_78, hl_79, \
                         hl_80, hl_81, hl_82, hl_83, hl_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * hl_77[k];

        t_33[k] = f_0 * hl_78[k];

        t_34[k] = f_0 * hl_79[k];

        t_35[k] = f_0 * hl_80[k];

        t_36[k] = f_0 * hl_81[k];

        t_37[k] = f_0 * hl_82[k];

        t_38[k] = f_0 * hl_83[k];

        t_39[k] = f_0 * hl_84[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, fl_0, fl_1, hl_85, hl_86, \
                         hl_87, hl_88, hl_89, hl_135, hl_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * hl_85[k];

        t_41[k] = f_0 * hl_86[k];

        t_42[k] = f_0 * hl_87[k];

        t_43[k] = f_0 * hl_88[k];

        t_44[k] = f_0 * hl_89[k];

        t_45[k] = -fl_0[k]
                  + f_0 * hl_135[k];

        t_46[k] = -fl_1[k]
                  + f_0 * hl_136[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, fl_2, fl_3, fl_4, fl_5, fl_6, hl_137, \
                         hl_138, hl_139, hl_140, hl_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -fl_2[k]
                  + f_0 * hl_137[k];

        t_48[k] = -fl_3[k]
                  + f_0 * hl_138[k];

        t_49[k] = -fl_4[k]
                  + f_0 * hl_139[k];

        t_50[k] = -fl_5[k]
                  + f_0 * hl_140[k];

        t_51[k] = -fl_6[k]
                  + f_0 * hl_141[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, fl_7, fl_8, fl_9, fl_10, fl_11, hl_142, \
                         hl_143, hl_144, hl_145, hl_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = -fl_7[k]
                  + f_0 * hl_142[k];

        t_53[k] = -fl_8[k]
                  + f_0 * hl_143[k];

        t_54[k] = -fl_9[k]
                  + f_0 * hl_144[k];

        t_55[k] = -fl_10[k]
                  + f_0 * hl_145[k];

        t_56[k] = -fl_11[k]
                  + f_0 * hl_146[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, fl_12, fl_13, fl_14, fl_15, fl_16, \
                         hl_147, hl_148, hl_149, hl_150, hl_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -fl_12[k]
                  + f_0 * hl_147[k];

        t_58[k] = -fl_13[k]
                  + f_0 * hl_148[k];

        t_59[k] = -fl_14[k]
                  + f_0 * hl_149[k];

        t_60[k] = -fl_15[k]
                  + f_0 * hl_150[k];

        t_61[k] = -fl_16[k]
                  + f_0 * hl_151[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, fl_17, fl_18, fl_19, fl_20, fl_21, \
                         hl_152, hl_153, hl_154, hl_155, hl_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -fl_17[k]
                  + f_0 * hl_152[k];

        t_63[k] = -fl_18[k]
                  + f_0 * hl_153[k];

        t_64[k] = -fl_19[k]
                  + f_0 * hl_154[k];

        t_65[k] = -fl_20[k]
                  + f_0 * hl_155[k];

        t_66[k] = -fl_21[k]
                  + f_0 * hl_156[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, fl_22, fl_23, fl_24, fl_25, fl_26, \
                         hl_157, hl_158, hl_159, hl_160, hl_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -fl_22[k]
                  + f_0 * hl_157[k];

        t_68[k] = -fl_23[k]
                  + f_0 * hl_158[k];

        t_69[k] = -fl_24[k]
                  + f_0 * hl_159[k];

        t_70[k] = -fl_25[k]
                  + f_0 * hl_160[k];

        t_71[k] = -fl_26[k]
                  + f_0 * hl_161[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, fl_27, fl_28, fl_29, fl_30, fl_31, \
                         hl_162, hl_163, hl_164, hl_165, hl_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -fl_27[k]
                  + f_0 * hl_162[k];

        t_73[k] = -fl_28[k]
                  + f_0 * hl_163[k];

        t_74[k] = -fl_29[k]
                  + f_0 * hl_164[k];

        t_75[k] = -fl_30[k]
                  + f_0 * hl_165[k];

        t_76[k] = -fl_31[k]
                  + f_0 * hl_166[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, fl_32, fl_33, fl_34, fl_35, fl_36, \
                         hl_167, hl_168, hl_169, hl_170, hl_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -fl_32[k]
                  + f_0 * hl_167[k];

        t_78[k] = -fl_33[k]
                  + f_0 * hl_168[k];

        t_79[k] = -fl_34[k]
                  + f_0 * hl_169[k];

        t_80[k] = -fl_35[k]
                  + f_0 * hl_170[k];

        t_81[k] = -fl_36[k]
                  + f_0 * hl_171[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, fl_37, fl_38, fl_39, fl_40, fl_41, \
                         hl_172, hl_173, hl_174, hl_175, hl_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -fl_37[k]
                  + f_0 * hl_172[k];

        t_83[k] = -fl_38[k]
                  + f_0 * hl_173[k];

        t_84[k] = -fl_39[k]
                  + f_0 * hl_174[k];

        t_85[k] = -fl_40[k]
                  + f_0 * hl_175[k];

        t_86[k] = -fl_41[k]
                  + f_0 * hl_176[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, fl_42, fl_43, fl_44, hl_177, \
                         hl_178, hl_179, hl_180, hl_181, hl_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -fl_42[k]
                  + f_0 * hl_177[k];

        t_88[k] = -fl_43[k]
                  + f_0 * hl_178[k];

        t_89[k] = -fl_44[k]
                  + f_0 * hl_179[k];

        t_90[k] = f_0 * hl_180[k];

        t_91[k] = f_0 * hl_181[k];

        t_92[k] = f_0 * hl_182[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, hl_183, hl_184, \
                         hl_185, hl_186, hl_187, hl_188, hl_189, \
                         hl_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * hl_183[k];

        t_94[k] = f_0 * hl_184[k];

        t_95[k] = f_0 * hl_185[k];

        t_96[k] = f_0 * hl_186[k];

        t_97[k] = f_0 * hl_187[k];

        t_98[k] = f_0 * hl_188[k];

        t_99[k] = f_0 * hl_189[k];

        t_100[k] = f_0 * hl_190[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, t_107, t_108, hl_191, \
                         hl_192, hl_193, hl_194, hl_195, hl_196, hl_197, \
                         hl_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * hl_191[k];

        t_102[k] = f_0 * hl_192[k];

        t_103[k] = f_0 * hl_193[k];

        t_104[k] = f_0 * hl_194[k];

        t_105[k] = f_0 * hl_195[k];

        t_106[k] = f_0 * hl_196[k];

        t_107[k] = f_0 * hl_197[k];

        t_108[k] = f_0 * hl_198[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, t_114, t_115, t_116, hl_199, \
                         hl_200, hl_201, hl_202, hl_203, hl_204, hl_205, \
                         hl_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_0 * hl_199[k];

        t_110[k] = f_0 * hl_200[k];

        t_111[k] = f_0 * hl_201[k];

        t_112[k] = f_0 * hl_202[k];

        t_113[k] = f_0 * hl_203[k];

        t_114[k] = f_0 * hl_204[k];

        t_115[k] = f_0 * hl_205[k];

        t_116[k] = f_0 * hl_206[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, t_122, t_123, t_124, hl_207, \
                         hl_208, hl_209, hl_210, hl_211, hl_212, hl_213, \
                         hl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_0 * hl_207[k];

        t_118[k] = f_0 * hl_208[k];

        t_119[k] = f_0 * hl_209[k];

        t_120[k] = f_0 * hl_210[k];

        t_121[k] = f_0 * hl_211[k];

        t_122[k] = f_0 * hl_212[k];

        t_123[k] = f_0 * hl_213[k];

        t_124[k] = f_0 * hl_214[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, t_131, t_132, hl_215, \
                         hl_216, hl_217, hl_218, hl_219, hl_220, hl_221, \
                         hl_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * hl_215[k];

        t_126[k] = f_0 * hl_216[k];

        t_127[k] = f_0 * hl_217[k];

        t_128[k] = f_0 * hl_218[k];

        t_129[k] = f_0 * hl_219[k];

        t_130[k] = f_0 * hl_220[k];

        t_131[k] = f_0 * hl_221[k];

        t_132[k] = f_0 * hl_222[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, t_138, fl_45, fl_46, fl_47, fl_48, \
                         hl_223, hl_224, hl_270, hl_271, hl_272, \
                         hl_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_0 * hl_223[k];

        t_134[k] = f_0 * hl_224[k];

        t_135[k] = -2.0 * fl_45[k]
                   + f_0 * hl_270[k];

        t_136[k] = -2.0 * fl_46[k]
                   + f_0 * hl_271[k];

        t_137[k] = -2.0 * fl_47[k]
                   + f_0 * hl_272[k];

        t_138[k] = -2.0 * fl_48[k]
                   + f_0 * hl_273[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, fl_49, fl_50, fl_51, fl_52, fl_53, \
                         hl_274, hl_275, hl_276, hl_277, hl_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -2.0 * fl_49[k]
                   + f_0 * hl_274[k];

        t_140[k] = -2.0 * fl_50[k]
                   + f_0 * hl_275[k];

        t_141[k] = -2.0 * fl_51[k]
                   + f_0 * hl_276[k];

        t_142[k] = -2.0 * fl_52[k]
                   + f_0 * hl_277[k];

        t_143[k] = -2.0 * fl_53[k]
                   + f_0 * hl_278[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, fl_54, fl_55, fl_56, fl_57, fl_58, \
                         hl_279, hl_280, hl_281, hl_282, hl_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -2.0 * fl_54[k]
                   + f_0 * hl_279[k];

        t_145[k] = -2.0 * fl_55[k]
                   + f_0 * hl_280[k];

        t_146[k] = -2.0 * fl_56[k]
                   + f_0 * hl_281[k];

        t_147[k] = -2.0 * fl_57[k]
                   + f_0 * hl_282[k];

        t_148[k] = -2.0 * fl_58[k]
                   + f_0 * hl_283[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, fl_59, fl_60, fl_61, fl_62, fl_63, \
                         hl_284, hl_285, hl_286, hl_287, hl_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -2.0 * fl_59[k]
                   + f_0 * hl_284[k];

        t_150[k] = -2.0 * fl_60[k]
                   + f_0 * hl_285[k];

        t_151[k] = -2.0 * fl_61[k]
                   + f_0 * hl_286[k];

        t_152[k] = -2.0 * fl_62[k]
                   + f_0 * hl_287[k];

        t_153[k] = -2.0 * fl_63[k]
                   + f_0 * hl_288[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, fl_64, fl_65, fl_66, fl_67, fl_68, \
                         hl_289, hl_290, hl_291, hl_292, hl_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -2.0 * fl_64[k]
                   + f_0 * hl_289[k];

        t_155[k] = -2.0 * fl_65[k]
                   + f_0 * hl_290[k];

        t_156[k] = -2.0 * fl_66[k]
                   + f_0 * hl_291[k];

        t_157[k] = -2.0 * fl_67[k]
                   + f_0 * hl_292[k];

        t_158[k] = -2.0 * fl_68[k]
                   + f_0 * hl_293[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, fl_69, fl_70, fl_71, fl_72, fl_73, \
                         hl_294, hl_295, hl_296, hl_297, hl_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -2.0 * fl_69[k]
                   + f_0 * hl_294[k];

        t_160[k] = -2.0 * fl_70[k]
                   + f_0 * hl_295[k];

        t_161[k] = -2.0 * fl_71[k]
                   + f_0 * hl_296[k];

        t_162[k] = -2.0 * fl_72[k]
                   + f_0 * hl_297[k];

        t_163[k] = -2.0 * fl_73[k]
                   + f_0 * hl_298[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, fl_74, fl_75, fl_76, fl_77, fl_78, \
                         hl_299, hl_300, hl_301, hl_302, hl_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -2.0 * fl_74[k]
                   + f_0 * hl_299[k];

        t_165[k] = -2.0 * fl_75[k]
                   + f_0 * hl_300[k];

        t_166[k] = -2.0 * fl_76[k]
                   + f_0 * hl_301[k];

        t_167[k] = -2.0 * fl_77[k]
                   + f_0 * hl_302[k];

        t_168[k] = -2.0 * fl_78[k]
                   + f_0 * hl_303[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, fl_79, fl_80, fl_81, fl_82, fl_83, \
                         hl_304, hl_305, hl_306, hl_307, hl_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -2.0 * fl_79[k]
                   + f_0 * hl_304[k];

        t_170[k] = -2.0 * fl_80[k]
                   + f_0 * hl_305[k];

        t_171[k] = -2.0 * fl_81[k]
                   + f_0 * hl_306[k];

        t_172[k] = -2.0 * fl_82[k]
                   + f_0 * hl_307[k];

        t_173[k] = -2.0 * fl_83[k]
                   + f_0 * hl_308[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, fl_84, fl_85, fl_86, fl_87, fl_88, \
                         hl_309, hl_310, hl_311, hl_312, hl_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -2.0 * fl_84[k]
                   + f_0 * hl_309[k];

        t_175[k] = -2.0 * fl_85[k]
                   + f_0 * hl_310[k];

        t_176[k] = -2.0 * fl_86[k]
                   + f_0 * hl_311[k];

        t_177[k] = -2.0 * fl_87[k]
                   + f_0 * hl_312[k];

        t_178[k] = -2.0 * fl_88[k]
                   + f_0 * hl_313[k];
    }
}

static auto
compute_prim_geom_10_gl_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);

    const auto *fl_89 = buffer.data(fl + 89);
    const auto *fl_90 = buffer.data(fl + 90);
    const auto *fl_91 = buffer.data(fl + 91);
    const auto *fl_92 = buffer.data(fl + 92);
    const auto *fl_93 = buffer.data(fl + 93);
    const auto *fl_94 = buffer.data(fl + 94);
    const auto *fl_95 = buffer.data(fl + 95);
    const auto *fl_96 = buffer.data(fl + 96);
    const auto *fl_97 = buffer.data(fl + 97);
    const auto *fl_98 = buffer.data(fl + 98);
    const auto *fl_99 = buffer.data(fl + 99);
    const auto *fl_100 = buffer.data(fl + 100);
    const auto *fl_101 = buffer.data(fl + 101);
    const auto *fl_102 = buffer.data(fl + 102);
    const auto *fl_103 = buffer.data(fl + 103);
    const auto *fl_104 = buffer.data(fl + 104);
    const auto *fl_105 = buffer.data(fl + 105);
    const auto *fl_106 = buffer.data(fl + 106);
    const auto *fl_107 = buffer.data(fl + 107);
    const auto *fl_108 = buffer.data(fl + 108);
    const auto *fl_109 = buffer.data(fl + 109);
    const auto *fl_110 = buffer.data(fl + 110);
    const auto *fl_111 = buffer.data(fl + 111);
    const auto *fl_112 = buffer.data(fl + 112);
    const auto *fl_113 = buffer.data(fl + 113);
    const auto *fl_114 = buffer.data(fl + 114);
    const auto *fl_115 = buffer.data(fl + 115);
    const auto *fl_116 = buffer.data(fl + 116);
    const auto *fl_117 = buffer.data(fl + 117);
    const auto *fl_118 = buffer.data(fl + 118);
    const auto *fl_119 = buffer.data(fl + 119);
    const auto *fl_120 = buffer.data(fl + 120);
    const auto *fl_121 = buffer.data(fl + 121);
    const auto *fl_122 = buffer.data(fl + 122);
    const auto *fl_123 = buffer.data(fl + 123);
    const auto *fl_124 = buffer.data(fl + 124);
    const auto *fl_125 = buffer.data(fl + 125);
    const auto *fl_126 = buffer.data(fl + 126);
    const auto *fl_127 = buffer.data(fl + 127);
    const auto *fl_128 = buffer.data(fl + 128);
    const auto *fl_129 = buffer.data(fl + 129);
    const auto *fl_130 = buffer.data(fl + 130);
    const auto *fl_131 = buffer.data(fl + 131);
    const auto *fl_132 = buffer.data(fl + 132);
    const auto *fl_133 = buffer.data(fl + 133);
    const auto *fl_134 = buffer.data(fl + 134);
    const auto *fl_135 = buffer.data(fl + 135);
    const auto *fl_136 = buffer.data(fl + 136);
    const auto *fl_137 = buffer.data(fl + 137);
    const auto *fl_138 = buffer.data(fl + 138);
    const auto *fl_139 = buffer.data(fl + 139);
    const auto *fl_140 = buffer.data(fl + 140);
    const auto *fl_141 = buffer.data(fl + 141);
    const auto *fl_142 = buffer.data(fl + 142);
    const auto *fl_143 = buffer.data(fl + 143);
    const auto *fl_144 = buffer.data(fl + 144);
    const auto *fl_145 = buffer.data(fl + 145);
    const auto *fl_146 = buffer.data(fl + 146);
    const auto *fl_147 = buffer.data(fl + 147);
    const auto *fl_148 = buffer.data(fl + 148);
    const auto *fl_149 = buffer.data(fl + 149);
    const auto *fl_150 = buffer.data(fl + 150);
    const auto *fl_151 = buffer.data(fl + 151);
    const auto *fl_152 = buffer.data(fl + 152);
    const auto *fl_153 = buffer.data(fl + 153);
    const auto *fl_154 = buffer.data(fl + 154);
    const auto *fl_155 = buffer.data(fl + 155);
    const auto *fl_156 = buffer.data(fl + 156);
    const auto *fl_157 = buffer.data(fl + 157);
    const auto *fl_158 = buffer.data(fl + 158);
    const auto *fl_159 = buffer.data(fl + 159);
    const auto *fl_160 = buffer.data(fl + 160);
    const auto *fl_161 = buffer.data(fl + 161);
    const auto *fl_162 = buffer.data(fl + 162);
    const auto *fl_163 = buffer.data(fl + 163);
    const auto *fl_164 = buffer.data(fl + 164);
    const auto *fl_165 = buffer.data(fl + 165);
    const auto *fl_166 = buffer.data(fl + 166);
    const auto *fl_167 = buffer.data(fl + 167);
    const auto *fl_168 = buffer.data(fl + 168);
    const auto *fl_169 = buffer.data(fl + 169);
    const auto *fl_170 = buffer.data(fl + 170);
    const auto *fl_171 = buffer.data(fl + 171);
    const auto *fl_172 = buffer.data(fl + 172);
    const auto *fl_173 = buffer.data(fl + 173);
    const auto *fl_174 = buffer.data(fl + 174);
    const auto *fl_175 = buffer.data(fl + 175);
    const auto *fl_176 = buffer.data(fl + 176);
    const auto *fl_177 = buffer.data(fl + 177);
    const auto *fl_178 = buffer.data(fl + 178);
    const auto *fl_179 = buffer.data(fl + 179);
    const auto *fl_180 = buffer.data(fl + 180);
    const auto *fl_181 = buffer.data(fl + 181);
    const auto *fl_182 = buffer.data(fl + 182);
    const auto *fl_183 = buffer.data(fl + 183);
    const auto *fl_184 = buffer.data(fl + 184);
    const auto *fl_185 = buffer.data(fl + 185);
    const auto *fl_186 = buffer.data(fl + 186);
    const auto *fl_187 = buffer.data(fl + 187);
    const auto *fl_188 = buffer.data(fl + 188);
    const auto *fl_189 = buffer.data(fl + 189);
    const auto *fl_190 = buffer.data(fl + 190);
    const auto *fl_191 = buffer.data(fl + 191);
    const auto *fl_192 = buffer.data(fl + 192);
    const auto *fl_193 = buffer.data(fl + 193);
    const auto *fl_194 = buffer.data(fl + 194);
    const auto *fl_195 = buffer.data(fl + 195);
    const auto *fl_196 = buffer.data(fl + 196);
    const auto *fl_197 = buffer.data(fl + 197);
    const auto *fl_198 = buffer.data(fl + 198);
    const auto *fl_199 = buffer.data(fl + 199);
    const auto *fl_200 = buffer.data(fl + 200);
    const auto *fl_201 = buffer.data(fl + 201);
    const auto *fl_202 = buffer.data(fl + 202);
    const auto *fl_203 = buffer.data(fl + 203);
    const auto *fl_204 = buffer.data(fl + 204);
    const auto *fl_205 = buffer.data(fl + 205);
    const auto *fl_206 = buffer.data(fl + 206);
    const auto *fl_207 = buffer.data(fl + 207);
    const auto *fl_208 = buffer.data(fl + 208);
    const auto *fl_209 = buffer.data(fl + 209);

    const auto *hl_314 = buffer.data(hl + 314);
    const auto *hl_315 = buffer.data(hl + 315);
    const auto *hl_316 = buffer.data(hl + 316);
    const auto *hl_317 = buffer.data(hl + 317);
    const auto *hl_318 = buffer.data(hl + 318);
    const auto *hl_319 = buffer.data(hl + 319);
    const auto *hl_320 = buffer.data(hl + 320);
    const auto *hl_321 = buffer.data(hl + 321);
    const auto *hl_322 = buffer.data(hl + 322);
    const auto *hl_323 = buffer.data(hl + 323);
    const auto *hl_324 = buffer.data(hl + 324);
    const auto *hl_325 = buffer.data(hl + 325);
    const auto *hl_326 = buffer.data(hl + 326);
    const auto *hl_327 = buffer.data(hl + 327);
    const auto *hl_328 = buffer.data(hl + 328);
    const auto *hl_329 = buffer.data(hl + 329);
    const auto *hl_330 = buffer.data(hl + 330);
    const auto *hl_331 = buffer.data(hl + 331);
    const auto *hl_332 = buffer.data(hl + 332);
    const auto *hl_333 = buffer.data(hl + 333);
    const auto *hl_334 = buffer.data(hl + 334);
    const auto *hl_335 = buffer.data(hl + 335);
    const auto *hl_336 = buffer.data(hl + 336);
    const auto *hl_337 = buffer.data(hl + 337);
    const auto *hl_338 = buffer.data(hl + 338);
    const auto *hl_339 = buffer.data(hl + 339);
    const auto *hl_340 = buffer.data(hl + 340);
    const auto *hl_341 = buffer.data(hl + 341);
    const auto *hl_342 = buffer.data(hl + 342);
    const auto *hl_343 = buffer.data(hl + 343);
    const auto *hl_344 = buffer.data(hl + 344);
    const auto *hl_345 = buffer.data(hl + 345);
    const auto *hl_346 = buffer.data(hl + 346);
    const auto *hl_347 = buffer.data(hl + 347);
    const auto *hl_348 = buffer.data(hl + 348);
    const auto *hl_349 = buffer.data(hl + 349);
    const auto *hl_350 = buffer.data(hl + 350);
    const auto *hl_351 = buffer.data(hl + 351);
    const auto *hl_352 = buffer.data(hl + 352);
    const auto *hl_353 = buffer.data(hl + 353);
    const auto *hl_354 = buffer.data(hl + 354);
    const auto *hl_355 = buffer.data(hl + 355);
    const auto *hl_356 = buffer.data(hl + 356);
    const auto *hl_357 = buffer.data(hl + 357);
    const auto *hl_358 = buffer.data(hl + 358);
    const auto *hl_359 = buffer.data(hl + 359);
    const auto *hl_360 = buffer.data(hl + 360);
    const auto *hl_361 = buffer.data(hl + 361);
    const auto *hl_362 = buffer.data(hl + 362);
    const auto *hl_363 = buffer.data(hl + 363);
    const auto *hl_364 = buffer.data(hl + 364);
    const auto *hl_365 = buffer.data(hl + 365);
    const auto *hl_366 = buffer.data(hl + 366);
    const auto *hl_367 = buffer.data(hl + 367);
    const auto *hl_368 = buffer.data(hl + 368);
    const auto *hl_369 = buffer.data(hl + 369);
    const auto *hl_370 = buffer.data(hl + 370);
    const auto *hl_371 = buffer.data(hl + 371);
    const auto *hl_372 = buffer.data(hl + 372);
    const auto *hl_373 = buffer.data(hl + 373);
    const auto *hl_374 = buffer.data(hl + 374);
    const auto *hl_375 = buffer.data(hl + 375);
    const auto *hl_376 = buffer.data(hl + 376);
    const auto *hl_377 = buffer.data(hl + 377);
    const auto *hl_378 = buffer.data(hl + 378);
    const auto *hl_379 = buffer.data(hl + 379);
    const auto *hl_380 = buffer.data(hl + 380);
    const auto *hl_381 = buffer.data(hl + 381);
    const auto *hl_382 = buffer.data(hl + 382);
    const auto *hl_383 = buffer.data(hl + 383);
    const auto *hl_384 = buffer.data(hl + 384);
    const auto *hl_385 = buffer.data(hl + 385);
    const auto *hl_386 = buffer.data(hl + 386);
    const auto *hl_387 = buffer.data(hl + 387);
    const auto *hl_388 = buffer.data(hl + 388);
    const auto *hl_389 = buffer.data(hl + 389);
    const auto *hl_390 = buffer.data(hl + 390);
    const auto *hl_391 = buffer.data(hl + 391);
    const auto *hl_392 = buffer.data(hl + 392);
    const auto *hl_393 = buffer.data(hl + 393);
    const auto *hl_394 = buffer.data(hl + 394);
    const auto *hl_395 = buffer.data(hl + 395);
    const auto *hl_396 = buffer.data(hl + 396);
    const auto *hl_397 = buffer.data(hl + 397);
    const auto *hl_398 = buffer.data(hl + 398);
    const auto *hl_399 = buffer.data(hl + 399);
    const auto *hl_400 = buffer.data(hl + 400);
    const auto *hl_401 = buffer.data(hl + 401);
    const auto *hl_402 = buffer.data(hl + 402);
    const auto *hl_403 = buffer.data(hl + 403);
    const auto *hl_404 = buffer.data(hl + 404);
    const auto *hl_450 = buffer.data(hl + 450);
    const auto *hl_451 = buffer.data(hl + 451);
    const auto *hl_452 = buffer.data(hl + 452);
    const auto *hl_453 = buffer.data(hl + 453);
    const auto *hl_454 = buffer.data(hl + 454);
    const auto *hl_455 = buffer.data(hl + 455);
    const auto *hl_456 = buffer.data(hl + 456);
    const auto *hl_457 = buffer.data(hl + 457);
    const auto *hl_458 = buffer.data(hl + 458);
    const auto *hl_459 = buffer.data(hl + 459);
    const auto *hl_460 = buffer.data(hl + 460);
    const auto *hl_461 = buffer.data(hl + 461);
    const auto *hl_462 = buffer.data(hl + 462);
    const auto *hl_463 = buffer.data(hl + 463);
    const auto *hl_464 = buffer.data(hl + 464);
    const auto *hl_465 = buffer.data(hl + 465);
    const auto *hl_466 = buffer.data(hl + 466);
    const auto *hl_467 = buffer.data(hl + 467);
    const auto *hl_468 = buffer.data(hl + 468);
    const auto *hl_469 = buffer.data(hl + 469);
    const auto *hl_470 = buffer.data(hl + 470);
    const auto *hl_471 = buffer.data(hl + 471);
    const auto *hl_472 = buffer.data(hl + 472);
    const auto *hl_473 = buffer.data(hl + 473);
    const auto *hl_474 = buffer.data(hl + 474);
    const auto *hl_475 = buffer.data(hl + 475);
    const auto *hl_476 = buffer.data(hl + 476);
    const auto *hl_477 = buffer.data(hl + 477);
    const auto *hl_478 = buffer.data(hl + 478);
    const auto *hl_479 = buffer.data(hl + 479);
    const auto *hl_480 = buffer.data(hl + 480);
    const auto *hl_481 = buffer.data(hl + 481);
    const auto *hl_482 = buffer.data(hl + 482);
    const auto *hl_483 = buffer.data(hl + 483);
    const auto *hl_484 = buffer.data(hl + 484);
    const auto *hl_485 = buffer.data(hl + 485);
    const auto *hl_486 = buffer.data(hl + 486);
    const auto *hl_487 = buffer.data(hl + 487);
    const auto *hl_488 = buffer.data(hl + 488);
    const auto *hl_489 = buffer.data(hl + 489);
    const auto *hl_490 = buffer.data(hl + 490);
    const auto *hl_491 = buffer.data(hl + 491);
    const auto *hl_492 = buffer.data(hl + 492);
    const auto *hl_493 = buffer.data(hl + 493);
    const auto *hl_494 = buffer.data(hl + 494);
    const auto *hl_495 = buffer.data(hl + 495);
    const auto *hl_496 = buffer.data(hl + 496);
    const auto *hl_497 = buffer.data(hl + 497);
    const auto *hl_498 = buffer.data(hl + 498);
    const auto *hl_499 = buffer.data(hl + 499);
    const auto *hl_500 = buffer.data(hl + 500);
    const auto *hl_501 = buffer.data(hl + 501);
    const auto *hl_502 = buffer.data(hl + 502);
    const auto *hl_503 = buffer.data(hl + 503);
    const auto *hl_504 = buffer.data(hl + 504);
    const auto *hl_505 = buffer.data(hl + 505);
    const auto *hl_506 = buffer.data(hl + 506);
    const auto *hl_507 = buffer.data(hl + 507);
    const auto *hl_508 = buffer.data(hl + 508);
    const auto *hl_509 = buffer.data(hl + 509);
    const auto *hl_510 = buffer.data(hl + 510);
    const auto *hl_511 = buffer.data(hl + 511);
    const auto *hl_512 = buffer.data(hl + 512);
    const auto *hl_513 = buffer.data(hl + 513);
    const auto *hl_514 = buffer.data(hl + 514);
    const auto *hl_515 = buffer.data(hl + 515);
    const auto *hl_516 = buffer.data(hl + 516);
    const auto *hl_517 = buffer.data(hl + 517);
    const auto *hl_518 = buffer.data(hl + 518);
    const auto *hl_519 = buffer.data(hl + 519);
    const auto *hl_520 = buffer.data(hl + 520);
    const auto *hl_521 = buffer.data(hl + 521);
    const auto *hl_522 = buffer.data(hl + 522);
    const auto *hl_523 = buffer.data(hl + 523);
    const auto *hl_524 = buffer.data(hl + 524);

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, fl_89, fl_90, fl_91, fl_92, fl_93, \
                         hl_314, hl_315, hl_316, hl_317, hl_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -2.0 * fl_89[k]
                   + f_0 * hl_314[k];

        t_180[k] = -fl_90[k]
                   + f_0 * hl_315[k];

        t_181[k] = -fl_91[k]
                   + f_0 * hl_316[k];

        t_182[k] = -fl_92[k]
                   + f_0 * hl_317[k];

        t_183[k] = -fl_93[k]
                   + f_0 * hl_318[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, fl_94, fl_95, fl_96, fl_97, fl_98, \
                         hl_319, hl_320, hl_321, hl_322, hl_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -fl_94[k]
                   + f_0 * hl_319[k];

        t_185[k] = -fl_95[k]
                   + f_0 * hl_320[k];

        t_186[k] = -fl_96[k]
                   + f_0 * hl_321[k];

        t_187[k] = -fl_97[k]
                   + f_0 * hl_322[k];

        t_188[k] = -fl_98[k]
                   + f_0 * hl_323[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, fl_99, fl_100, fl_101, fl_102, \
                         fl_103, hl_324, hl_325, hl_326, hl_327, \
                         hl_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = -fl_99[k]
                   + f_0 * hl_324[k];

        t_190[k] = -fl_100[k]
                   + f_0 * hl_325[k];

        t_191[k] = -fl_101[k]
                   + f_0 * hl_326[k];

        t_192[k] = -fl_102[k]
                   + f_0 * hl_327[k];

        t_193[k] = -fl_103[k]
                   + f_0 * hl_328[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, fl_104, fl_105, fl_106, fl_107, \
                         fl_108, hl_329, hl_330, hl_331, hl_332, \
                         hl_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = -fl_104[k]
                   + f_0 * hl_329[k];

        t_195[k] = -fl_105[k]
                   + f_0 * hl_330[k];

        t_196[k] = -fl_106[k]
                   + f_0 * hl_331[k];

        t_197[k] = -fl_107[k]
                   + f_0 * hl_332[k];

        t_198[k] = -fl_108[k]
                   + f_0 * hl_333[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, fl_109, fl_110, fl_111, fl_112, \
                         fl_113, hl_334, hl_335, hl_336, hl_337, \
                         hl_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = -fl_109[k]
                   + f_0 * hl_334[k];

        t_200[k] = -fl_110[k]
                   + f_0 * hl_335[k];

        t_201[k] = -fl_111[k]
                   + f_0 * hl_336[k];

        t_202[k] = -fl_112[k]
                   + f_0 * hl_337[k];

        t_203[k] = -fl_113[k]
                   + f_0 * hl_338[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, fl_114, fl_115, fl_116, fl_117, \
                         fl_118, hl_339, hl_340, hl_341, hl_342, \
                         hl_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = -fl_114[k]
                   + f_0 * hl_339[k];

        t_205[k] = -fl_115[k]
                   + f_0 * hl_340[k];

        t_206[k] = -fl_116[k]
                   + f_0 * hl_341[k];

        t_207[k] = -fl_117[k]
                   + f_0 * hl_342[k];

        t_208[k] = -fl_118[k]
                   + f_0 * hl_343[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, fl_119, fl_120, fl_121, fl_122, \
                         fl_123, hl_344, hl_345, hl_346, hl_347, \
                         hl_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = -fl_119[k]
                   + f_0 * hl_344[k];

        t_210[k] = -fl_120[k]
                   + f_0 * hl_345[k];

        t_211[k] = -fl_121[k]
                   + f_0 * hl_346[k];

        t_212[k] = -fl_122[k]
                   + f_0 * hl_347[k];

        t_213[k] = -fl_123[k]
                   + f_0 * hl_348[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, fl_124, fl_125, fl_126, fl_127, \
                         fl_128, hl_349, hl_350, hl_351, hl_352, \
                         hl_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = -fl_124[k]
                   + f_0 * hl_349[k];

        t_215[k] = -fl_125[k]
                   + f_0 * hl_350[k];

        t_216[k] = -fl_126[k]
                   + f_0 * hl_351[k];

        t_217[k] = -fl_127[k]
                   + f_0 * hl_352[k];

        t_218[k] = -fl_128[k]
                   + f_0 * hl_353[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, fl_129, fl_130, fl_131, fl_132, \
                         fl_133, hl_354, hl_355, hl_356, hl_357, \
                         hl_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = -fl_129[k]
                   + f_0 * hl_354[k];

        t_220[k] = -fl_130[k]
                   + f_0 * hl_355[k];

        t_221[k] = -fl_131[k]
                   + f_0 * hl_356[k];

        t_222[k] = -fl_132[k]
                   + f_0 * hl_357[k];

        t_223[k] = -fl_133[k]
                   + f_0 * hl_358[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, t_229, t_230, fl_134, hl_359, \
                         hl_360, hl_361, hl_362, hl_363, hl_364, \
                         hl_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = -fl_134[k]
                   + f_0 * hl_359[k];

        t_225[k] = f_0 * hl_360[k];

        t_226[k] = f_0 * hl_361[k];

        t_227[k] = f_0 * hl_362[k];

        t_228[k] = f_0 * hl_363[k];

        t_229[k] = f_0 * hl_364[k];

        t_230[k] = f_0 * hl_365[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, t_236, t_237, t_238, hl_366, \
                         hl_367, hl_368, hl_369, hl_370, hl_371, hl_372, \
                         hl_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_0 * hl_366[k];

        t_232[k] = f_0 * hl_367[k];

        t_233[k] = f_0 * hl_368[k];

        t_234[k] = f_0 * hl_369[k];

        t_235[k] = f_0 * hl_370[k];

        t_236[k] = f_0 * hl_371[k];

        t_237[k] = f_0 * hl_372[k];

        t_238[k] = f_0 * hl_373[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, t_244, t_245, t_246, hl_374, \
                         hl_375, hl_376, hl_377, hl_378, hl_379, hl_380, \
                         hl_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_0 * hl_374[k];

        t_240[k] = f_0 * hl_375[k];

        t_241[k] = f_0 * hl_376[k];

        t_242[k] = f_0 * hl_377[k];

        t_243[k] = f_0 * hl_378[k];

        t_244[k] = f_0 * hl_379[k];

        t_245[k] = f_0 * hl_380[k];

        t_246[k] = f_0 * hl_381[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, t_252, t_253, t_254, hl_382, \
                         hl_383, hl_384, hl_385, hl_386, hl_387, hl_388, \
                         hl_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_0 * hl_382[k];

        t_248[k] = f_0 * hl_383[k];

        t_249[k] = f_0 * hl_384[k];

        t_250[k] = f_0 * hl_385[k];

        t_251[k] = f_0 * hl_386[k];

        t_252[k] = f_0 * hl_387[k];

        t_253[k] = f_0 * hl_388[k];

        t_254[k] = f_0 * hl_389[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, t_260, t_261, t_262, hl_390, \
                         hl_391, hl_392, hl_393, hl_394, hl_395, hl_396, \
                         hl_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_0 * hl_390[k];

        t_256[k] = f_0 * hl_391[k];

        t_257[k] = f_0 * hl_392[k];

        t_258[k] = f_0 * hl_393[k];

        t_259[k] = f_0 * hl_394[k];

        t_260[k] = f_0 * hl_395[k];

        t_261[k] = f_0 * hl_396[k];

        t_262[k] = f_0 * hl_397[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, t_268, t_269, hl_398, hl_399, \
                         hl_400, hl_401, hl_402, hl_403, hl_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_0 * hl_398[k];

        t_264[k] = f_0 * hl_399[k];

        t_265[k] = f_0 * hl_400[k];

        t_266[k] = f_0 * hl_401[k];

        t_267[k] = f_0 * hl_402[k];

        t_268[k] = f_0 * hl_403[k];

        t_269[k] = f_0 * hl_404[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, fl_135, fl_136, fl_137, fl_138, \
                         fl_139, hl_450, hl_451, hl_452, hl_453, \
                         hl_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -3.0 * fl_135[k]
                   + f_0 * hl_450[k];

        t_271[k] = -3.0 * fl_136[k]
                   + f_0 * hl_451[k];

        t_272[k] = -3.0 * fl_137[k]
                   + f_0 * hl_452[k];

        t_273[k] = -3.0 * fl_138[k]
                   + f_0 * hl_453[k];

        t_274[k] = -3.0 * fl_139[k]
                   + f_0 * hl_454[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, fl_140, fl_141, fl_142, fl_143, \
                         fl_144, hl_455, hl_456, hl_457, hl_458, \
                         hl_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * fl_140[k]
                   + f_0 * hl_455[k];

        t_276[k] = -3.0 * fl_141[k]
                   + f_0 * hl_456[k];

        t_277[k] = -3.0 * fl_142[k]
                   + f_0 * hl_457[k];

        t_278[k] = -3.0 * fl_143[k]
                   + f_0 * hl_458[k];

        t_279[k] = -3.0 * fl_144[k]
                   + f_0 * hl_459[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, fl_145, fl_146, fl_147, fl_148, \
                         fl_149, hl_460, hl_461, hl_462, hl_463, \
                         hl_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * fl_145[k]
                   + f_0 * hl_460[k];

        t_281[k] = -3.0 * fl_146[k]
                   + f_0 * hl_461[k];

        t_282[k] = -3.0 * fl_147[k]
                   + f_0 * hl_462[k];

        t_283[k] = -3.0 * fl_148[k]
                   + f_0 * hl_463[k];

        t_284[k] = -3.0 * fl_149[k]
                   + f_0 * hl_464[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, fl_150, fl_151, fl_152, fl_153, \
                         fl_154, hl_465, hl_466, hl_467, hl_468, \
                         hl_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * fl_150[k]
                   + f_0 * hl_465[k];

        t_286[k] = -3.0 * fl_151[k]
                   + f_0 * hl_466[k];

        t_287[k] = -3.0 * fl_152[k]
                   + f_0 * hl_467[k];

        t_288[k] = -3.0 * fl_153[k]
                   + f_0 * hl_468[k];

        t_289[k] = -3.0 * fl_154[k]
                   + f_0 * hl_469[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, fl_155, fl_156, fl_157, fl_158, \
                         fl_159, hl_470, hl_471, hl_472, hl_473, \
                         hl_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * fl_155[k]
                   + f_0 * hl_470[k];

        t_291[k] = -3.0 * fl_156[k]
                   + f_0 * hl_471[k];

        t_292[k] = -3.0 * fl_157[k]
                   + f_0 * hl_472[k];

        t_293[k] = -3.0 * fl_158[k]
                   + f_0 * hl_473[k];

        t_294[k] = -3.0 * fl_159[k]
                   + f_0 * hl_474[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, fl_160, fl_161, fl_162, fl_163, \
                         fl_164, hl_475, hl_476, hl_477, hl_478, \
                         hl_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -3.0 * fl_160[k]
                   + f_0 * hl_475[k];

        t_296[k] = -3.0 * fl_161[k]
                   + f_0 * hl_476[k];

        t_297[k] = -3.0 * fl_162[k]
                   + f_0 * hl_477[k];

        t_298[k] = -3.0 * fl_163[k]
                   + f_0 * hl_478[k];

        t_299[k] = -3.0 * fl_164[k]
                   + f_0 * hl_479[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, fl_165, fl_166, fl_167, fl_168, \
                         fl_169, hl_480, hl_481, hl_482, hl_483, \
                         hl_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -3.0 * fl_165[k]
                   + f_0 * hl_480[k];

        t_301[k] = -3.0 * fl_166[k]
                   + f_0 * hl_481[k];

        t_302[k] = -3.0 * fl_167[k]
                   + f_0 * hl_482[k];

        t_303[k] = -3.0 * fl_168[k]
                   + f_0 * hl_483[k];

        t_304[k] = -3.0 * fl_169[k]
                   + f_0 * hl_484[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, fl_170, fl_171, fl_172, fl_173, \
                         fl_174, hl_485, hl_486, hl_487, hl_488, \
                         hl_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -3.0 * fl_170[k]
                   + f_0 * hl_485[k];

        t_306[k] = -3.0 * fl_171[k]
                   + f_0 * hl_486[k];

        t_307[k] = -3.0 * fl_172[k]
                   + f_0 * hl_487[k];

        t_308[k] = -3.0 * fl_173[k]
                   + f_0 * hl_488[k];

        t_309[k] = -3.0 * fl_174[k]
                   + f_0 * hl_489[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, fl_175, fl_176, fl_177, fl_178, \
                         fl_179, hl_490, hl_491, hl_492, hl_493, \
                         hl_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * fl_175[k]
                   + f_0 * hl_490[k];

        t_311[k] = -3.0 * fl_176[k]
                   + f_0 * hl_491[k];

        t_312[k] = -3.0 * fl_177[k]
                   + f_0 * hl_492[k];

        t_313[k] = -3.0 * fl_178[k]
                   + f_0 * hl_493[k];

        t_314[k] = -3.0 * fl_179[k]
                   + f_0 * hl_494[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, fl_180, fl_181, fl_182, fl_183, \
                         fl_184, hl_495, hl_496, hl_497, hl_498, \
                         hl_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -2.0 * fl_180[k]
                   + f_0 * hl_495[k];

        t_316[k] = -2.0 * fl_181[k]
                   + f_0 * hl_496[k];

        t_317[k] = -2.0 * fl_182[k]
                   + f_0 * hl_497[k];

        t_318[k] = -2.0 * fl_183[k]
                   + f_0 * hl_498[k];

        t_319[k] = -2.0 * fl_184[k]
                   + f_0 * hl_499[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, fl_185, fl_186, fl_187, fl_188, \
                         fl_189, hl_500, hl_501, hl_502, hl_503, \
                         hl_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -2.0 * fl_185[k]
                   + f_0 * hl_500[k];

        t_321[k] = -2.0 * fl_186[k]
                   + f_0 * hl_501[k];

        t_322[k] = -2.0 * fl_187[k]
                   + f_0 * hl_502[k];

        t_323[k] = -2.0 * fl_188[k]
                   + f_0 * hl_503[k];

        t_324[k] = -2.0 * fl_189[k]
                   + f_0 * hl_504[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, fl_190, fl_191, fl_192, fl_193, \
                         fl_194, hl_505, hl_506, hl_507, hl_508, \
                         hl_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -2.0 * fl_190[k]
                   + f_0 * hl_505[k];

        t_326[k] = -2.0 * fl_191[k]
                   + f_0 * hl_506[k];

        t_327[k] = -2.0 * fl_192[k]
                   + f_0 * hl_507[k];

        t_328[k] = -2.0 * fl_193[k]
                   + f_0 * hl_508[k];

        t_329[k] = -2.0 * fl_194[k]
                   + f_0 * hl_509[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, fl_195, fl_196, fl_197, fl_198, \
                         fl_199, hl_510, hl_511, hl_512, hl_513, \
                         hl_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -2.0 * fl_195[k]
                   + f_0 * hl_510[k];

        t_331[k] = -2.0 * fl_196[k]
                   + f_0 * hl_511[k];

        t_332[k] = -2.0 * fl_197[k]
                   + f_0 * hl_512[k];

        t_333[k] = -2.0 * fl_198[k]
                   + f_0 * hl_513[k];

        t_334[k] = -2.0 * fl_199[k]
                   + f_0 * hl_514[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, fl_200, fl_201, fl_202, fl_203, \
                         fl_204, hl_515, hl_516, hl_517, hl_518, \
                         hl_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -2.0 * fl_200[k]
                   + f_0 * hl_515[k];

        t_336[k] = -2.0 * fl_201[k]
                   + f_0 * hl_516[k];

        t_337[k] = -2.0 * fl_202[k]
                   + f_0 * hl_517[k];

        t_338[k] = -2.0 * fl_203[k]
                   + f_0 * hl_518[k];

        t_339[k] = -2.0 * fl_204[k]
                   + f_0 * hl_519[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, fl_205, fl_206, fl_207, fl_208, \
                         fl_209, hl_520, hl_521, hl_522, hl_523, \
                         hl_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * fl_205[k]
                   + f_0 * hl_520[k];

        t_341[k] = -2.0 * fl_206[k]
                   + f_0 * hl_521[k];

        t_342[k] = -2.0 * fl_207[k]
                   + f_0 * hl_522[k];

        t_343[k] = -2.0 * fl_208[k]
                   + f_0 * hl_523[k];

        t_344[k] = -2.0 * fl_209[k]
                   + f_0 * hl_524[k];
    }
}

static auto
compute_prim_geom_10_gl_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *fl_210 = buffer.data(fl + 210);
    const auto *fl_211 = buffer.data(fl + 211);
    const auto *fl_212 = buffer.data(fl + 212);
    const auto *fl_213 = buffer.data(fl + 213);
    const auto *fl_214 = buffer.data(fl + 214);
    const auto *fl_215 = buffer.data(fl + 215);
    const auto *fl_216 = buffer.data(fl + 216);
    const auto *fl_217 = buffer.data(fl + 217);
    const auto *fl_218 = buffer.data(fl + 218);
    const auto *fl_219 = buffer.data(fl + 219);
    const auto *fl_220 = buffer.data(fl + 220);
    const auto *fl_221 = buffer.data(fl + 221);
    const auto *fl_222 = buffer.data(fl + 222);
    const auto *fl_223 = buffer.data(fl + 223);
    const auto *fl_224 = buffer.data(fl + 224);
    const auto *fl_225 = buffer.data(fl + 225);
    const auto *fl_226 = buffer.data(fl + 226);
    const auto *fl_227 = buffer.data(fl + 227);
    const auto *fl_228 = buffer.data(fl + 228);
    const auto *fl_229 = buffer.data(fl + 229);
    const auto *fl_230 = buffer.data(fl + 230);
    const auto *fl_231 = buffer.data(fl + 231);
    const auto *fl_232 = buffer.data(fl + 232);
    const auto *fl_233 = buffer.data(fl + 233);
    const auto *fl_234 = buffer.data(fl + 234);
    const auto *fl_235 = buffer.data(fl + 235);
    const auto *fl_236 = buffer.data(fl + 236);
    const auto *fl_237 = buffer.data(fl + 237);
    const auto *fl_238 = buffer.data(fl + 238);
    const auto *fl_239 = buffer.data(fl + 239);
    const auto *fl_240 = buffer.data(fl + 240);
    const auto *fl_241 = buffer.data(fl + 241);
    const auto *fl_242 = buffer.data(fl + 242);
    const auto *fl_243 = buffer.data(fl + 243);
    const auto *fl_244 = buffer.data(fl + 244);
    const auto *fl_245 = buffer.data(fl + 245);
    const auto *fl_246 = buffer.data(fl + 246);
    const auto *fl_247 = buffer.data(fl + 247);
    const auto *fl_248 = buffer.data(fl + 248);
    const auto *fl_249 = buffer.data(fl + 249);
    const auto *fl_250 = buffer.data(fl + 250);
    const auto *fl_251 = buffer.data(fl + 251);
    const auto *fl_252 = buffer.data(fl + 252);
    const auto *fl_253 = buffer.data(fl + 253);
    const auto *fl_254 = buffer.data(fl + 254);
    const auto *fl_255 = buffer.data(fl + 255);
    const auto *fl_256 = buffer.data(fl + 256);
    const auto *fl_257 = buffer.data(fl + 257);
    const auto *fl_258 = buffer.data(fl + 258);
    const auto *fl_259 = buffer.data(fl + 259);
    const auto *fl_260 = buffer.data(fl + 260);
    const auto *fl_261 = buffer.data(fl + 261);
    const auto *fl_262 = buffer.data(fl + 262);
    const auto *fl_263 = buffer.data(fl + 263);
    const auto *fl_264 = buffer.data(fl + 264);
    const auto *fl_265 = buffer.data(fl + 265);
    const auto *fl_266 = buffer.data(fl + 266);
    const auto *fl_267 = buffer.data(fl + 267);
    const auto *fl_268 = buffer.data(fl + 268);
    const auto *fl_269 = buffer.data(fl + 269);
    const auto *fl_270 = buffer.data(fl + 270);
    const auto *fl_271 = buffer.data(fl + 271);
    const auto *fl_272 = buffer.data(fl + 272);
    const auto *fl_273 = buffer.data(fl + 273);
    const auto *fl_274 = buffer.data(fl + 274);
    const auto *fl_275 = buffer.data(fl + 275);
    const auto *fl_276 = buffer.data(fl + 276);
    const auto *fl_277 = buffer.data(fl + 277);
    const auto *fl_278 = buffer.data(fl + 278);
    const auto *fl_279 = buffer.data(fl + 279);
    const auto *fl_280 = buffer.data(fl + 280);
    const auto *fl_281 = buffer.data(fl + 281);
    const auto *fl_282 = buffer.data(fl + 282);
    const auto *fl_283 = buffer.data(fl + 283);
    const auto *fl_284 = buffer.data(fl + 284);
    const auto *fl_285 = buffer.data(fl + 285);
    const auto *fl_286 = buffer.data(fl + 286);
    const auto *fl_287 = buffer.data(fl + 287);
    const auto *fl_288 = buffer.data(fl + 288);
    const auto *fl_289 = buffer.data(fl + 289);
    const auto *fl_290 = buffer.data(fl + 290);
    const auto *fl_291 = buffer.data(fl + 291);
    const auto *fl_292 = buffer.data(fl + 292);
    const auto *fl_293 = buffer.data(fl + 293);
    const auto *fl_294 = buffer.data(fl + 294);
    const auto *fl_295 = buffer.data(fl + 295);
    const auto *fl_296 = buffer.data(fl + 296);
    const auto *fl_297 = buffer.data(fl + 297);
    const auto *fl_298 = buffer.data(fl + 298);
    const auto *fl_299 = buffer.data(fl + 299);
    const auto *fl_300 = buffer.data(fl + 300);
    const auto *fl_301 = buffer.data(fl + 301);
    const auto *fl_302 = buffer.data(fl + 302);
    const auto *fl_303 = buffer.data(fl + 303);
    const auto *fl_304 = buffer.data(fl + 304);
    const auto *fl_305 = buffer.data(fl + 305);
    const auto *fl_306 = buffer.data(fl + 306);
    const auto *fl_307 = buffer.data(fl + 307);
    const auto *fl_308 = buffer.data(fl + 308);
    const auto *fl_309 = buffer.data(fl + 309);
    const auto *fl_310 = buffer.data(fl + 310);
    const auto *fl_311 = buffer.data(fl + 311);
    const auto *fl_312 = buffer.data(fl + 312);
    const auto *fl_313 = buffer.data(fl + 313);
    const auto *fl_314 = buffer.data(fl + 314);
    const auto *fl_315 = buffer.data(fl + 315);
    const auto *fl_316 = buffer.data(fl + 316);
    const auto *fl_317 = buffer.data(fl + 317);
    const auto *fl_318 = buffer.data(fl + 318);
    const auto *fl_319 = buffer.data(fl + 319);
    const auto *fl_320 = buffer.data(fl + 320);
    const auto *fl_321 = buffer.data(fl + 321);
    const auto *fl_322 = buffer.data(fl + 322);
    const auto *fl_323 = buffer.data(fl + 323);
    const auto *fl_324 = buffer.data(fl + 324);
    const auto *fl_325 = buffer.data(fl + 325);
    const auto *fl_326 = buffer.data(fl + 326);

    const auto *hl_525 = buffer.data(hl + 525);
    const auto *hl_526 = buffer.data(hl + 526);
    const auto *hl_527 = buffer.data(hl + 527);
    const auto *hl_528 = buffer.data(hl + 528);
    const auto *hl_529 = buffer.data(hl + 529);
    const auto *hl_530 = buffer.data(hl + 530);
    const auto *hl_531 = buffer.data(hl + 531);
    const auto *hl_532 = buffer.data(hl + 532);
    const auto *hl_533 = buffer.data(hl + 533);
    const auto *hl_534 = buffer.data(hl + 534);
    const auto *hl_535 = buffer.data(hl + 535);
    const auto *hl_536 = buffer.data(hl + 536);
    const auto *hl_537 = buffer.data(hl + 537);
    const auto *hl_538 = buffer.data(hl + 538);
    const auto *hl_539 = buffer.data(hl + 539);
    const auto *hl_540 = buffer.data(hl + 540);
    const auto *hl_541 = buffer.data(hl + 541);
    const auto *hl_542 = buffer.data(hl + 542);
    const auto *hl_543 = buffer.data(hl + 543);
    const auto *hl_544 = buffer.data(hl + 544);
    const auto *hl_545 = buffer.data(hl + 545);
    const auto *hl_546 = buffer.data(hl + 546);
    const auto *hl_547 = buffer.data(hl + 547);
    const auto *hl_548 = buffer.data(hl + 548);
    const auto *hl_549 = buffer.data(hl + 549);
    const auto *hl_550 = buffer.data(hl + 550);
    const auto *hl_551 = buffer.data(hl + 551);
    const auto *hl_552 = buffer.data(hl + 552);
    const auto *hl_553 = buffer.data(hl + 553);
    const auto *hl_554 = buffer.data(hl + 554);
    const auto *hl_555 = buffer.data(hl + 555);
    const auto *hl_556 = buffer.data(hl + 556);
    const auto *hl_557 = buffer.data(hl + 557);
    const auto *hl_558 = buffer.data(hl + 558);
    const auto *hl_559 = buffer.data(hl + 559);
    const auto *hl_560 = buffer.data(hl + 560);
    const auto *hl_561 = buffer.data(hl + 561);
    const auto *hl_562 = buffer.data(hl + 562);
    const auto *hl_563 = buffer.data(hl + 563);
    const auto *hl_564 = buffer.data(hl + 564);
    const auto *hl_565 = buffer.data(hl + 565);
    const auto *hl_566 = buffer.data(hl + 566);
    const auto *hl_567 = buffer.data(hl + 567);
    const auto *hl_568 = buffer.data(hl + 568);
    const auto *hl_569 = buffer.data(hl + 569);
    const auto *hl_570 = buffer.data(hl + 570);
    const auto *hl_571 = buffer.data(hl + 571);
    const auto *hl_572 = buffer.data(hl + 572);
    const auto *hl_573 = buffer.data(hl + 573);
    const auto *hl_574 = buffer.data(hl + 574);
    const auto *hl_575 = buffer.data(hl + 575);
    const auto *hl_576 = buffer.data(hl + 576);
    const auto *hl_577 = buffer.data(hl + 577);
    const auto *hl_578 = buffer.data(hl + 578);
    const auto *hl_579 = buffer.data(hl + 579);
    const auto *hl_580 = buffer.data(hl + 580);
    const auto *hl_581 = buffer.data(hl + 581);
    const auto *hl_582 = buffer.data(hl + 582);
    const auto *hl_583 = buffer.data(hl + 583);
    const auto *hl_584 = buffer.data(hl + 584);
    const auto *hl_585 = buffer.data(hl + 585);
    const auto *hl_586 = buffer.data(hl + 586);
    const auto *hl_587 = buffer.data(hl + 587);
    const auto *hl_588 = buffer.data(hl + 588);
    const auto *hl_589 = buffer.data(hl + 589);
    const auto *hl_590 = buffer.data(hl + 590);
    const auto *hl_591 = buffer.data(hl + 591);
    const auto *hl_592 = buffer.data(hl + 592);
    const auto *hl_593 = buffer.data(hl + 593);
    const auto *hl_594 = buffer.data(hl + 594);
    const auto *hl_595 = buffer.data(hl + 595);
    const auto *hl_596 = buffer.data(hl + 596);
    const auto *hl_597 = buffer.data(hl + 597);
    const auto *hl_598 = buffer.data(hl + 598);
    const auto *hl_599 = buffer.data(hl + 599);
    const auto *hl_600 = buffer.data(hl + 600);
    const auto *hl_601 = buffer.data(hl + 601);
    const auto *hl_602 = buffer.data(hl + 602);
    const auto *hl_603 = buffer.data(hl + 603);
    const auto *hl_604 = buffer.data(hl + 604);
    const auto *hl_605 = buffer.data(hl + 605);
    const auto *hl_606 = buffer.data(hl + 606);
    const auto *hl_607 = buffer.data(hl + 607);
    const auto *hl_608 = buffer.data(hl + 608);
    const auto *hl_609 = buffer.data(hl + 609);
    const auto *hl_610 = buffer.data(hl + 610);
    const auto *hl_611 = buffer.data(hl + 611);
    const auto *hl_612 = buffer.data(hl + 612);
    const auto *hl_613 = buffer.data(hl + 613);
    const auto *hl_614 = buffer.data(hl + 614);
    const auto *hl_615 = buffer.data(hl + 615);
    const auto *hl_616 = buffer.data(hl + 616);
    const auto *hl_617 = buffer.data(hl + 617);
    const auto *hl_618 = buffer.data(hl + 618);
    const auto *hl_619 = buffer.data(hl + 619);
    const auto *hl_620 = buffer.data(hl + 620);
    const auto *hl_621 = buffer.data(hl + 621);
    const auto *hl_622 = buffer.data(hl + 622);
    const auto *hl_623 = buffer.data(hl + 623);
    const auto *hl_624 = buffer.data(hl + 624);
    const auto *hl_625 = buffer.data(hl + 625);
    const auto *hl_626 = buffer.data(hl + 626);
    const auto *hl_627 = buffer.data(hl + 627);
    const auto *hl_628 = buffer.data(hl + 628);
    const auto *hl_629 = buffer.data(hl + 629);
    const auto *hl_675 = buffer.data(hl + 675);
    const auto *hl_676 = buffer.data(hl + 676);
    const auto *hl_677 = buffer.data(hl + 677);
    const auto *hl_678 = buffer.data(hl + 678);
    const auto *hl_679 = buffer.data(hl + 679);
    const auto *hl_680 = buffer.data(hl + 680);
    const auto *hl_681 = buffer.data(hl + 681);
    const auto *hl_682 = buffer.data(hl + 682);
    const auto *hl_683 = buffer.data(hl + 683);
    const auto *hl_684 = buffer.data(hl + 684);
    const auto *hl_685 = buffer.data(hl + 685);
    const auto *hl_686 = buffer.data(hl + 686);
    const auto *hl_687 = buffer.data(hl + 687);
    const auto *hl_688 = buffer.data(hl + 688);
    const auto *hl_689 = buffer.data(hl + 689);
    const auto *hl_690 = buffer.data(hl + 690);
    const auto *hl_691 = buffer.data(hl + 691);
    const auto *hl_692 = buffer.data(hl + 692);
    const auto *hl_693 = buffer.data(hl + 693);
    const auto *hl_694 = buffer.data(hl + 694);
    const auto *hl_695 = buffer.data(hl + 695);
    const auto *hl_696 = buffer.data(hl + 696);
    const auto *hl_697 = buffer.data(hl + 697);
    const auto *hl_698 = buffer.data(hl + 698);
    const auto *hl_699 = buffer.data(hl + 699);
    const auto *hl_700 = buffer.data(hl + 700);
    const auto *hl_701 = buffer.data(hl + 701);
    const auto *hl_702 = buffer.data(hl + 702);
    const auto *hl_703 = buffer.data(hl + 703);
    const auto *hl_704 = buffer.data(hl + 704);
    const auto *hl_705 = buffer.data(hl + 705);
    const auto *hl_706 = buffer.data(hl + 706);
    const auto *hl_707 = buffer.data(hl + 707);
    const auto *hl_708 = buffer.data(hl + 708);
    const auto *hl_709 = buffer.data(hl + 709);
    const auto *hl_710 = buffer.data(hl + 710);
    const auto *hl_711 = buffer.data(hl + 711);
    const auto *hl_712 = buffer.data(hl + 712);
    const auto *hl_713 = buffer.data(hl + 713);
    const auto *hl_714 = buffer.data(hl + 714);
    const auto *hl_715 = buffer.data(hl + 715);
    const auto *hl_716 = buffer.data(hl + 716);
    const auto *hl_717 = buffer.data(hl + 717);
    const auto *hl_718 = buffer.data(hl + 718);
    const auto *hl_719 = buffer.data(hl + 719);
    const auto *hl_720 = buffer.data(hl + 720);
    const auto *hl_721 = buffer.data(hl + 721);
    const auto *hl_722 = buffer.data(hl + 722);
    const auto *hl_723 = buffer.data(hl + 723);
    const auto *hl_724 = buffer.data(hl + 724);
    const auto *hl_725 = buffer.data(hl + 725);
    const auto *hl_726 = buffer.data(hl + 726);
    const auto *hl_727 = buffer.data(hl + 727);
    const auto *hl_728 = buffer.data(hl + 728);
    const auto *hl_729 = buffer.data(hl + 729);
    const auto *hl_730 = buffer.data(hl + 730);
    const auto *hl_731 = buffer.data(hl + 731);

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, fl_210, fl_211, fl_212, fl_213, \
                         fl_214, hl_525, hl_526, hl_527, hl_528, \
                         hl_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * fl_210[k]
                   + f_0 * hl_525[k];

        t_346[k] = -2.0 * fl_211[k]
                   + f_0 * hl_526[k];

        t_347[k] = -2.0 * fl_212[k]
                   + f_0 * hl_527[k];

        t_348[k] = -2.0 * fl_213[k]
                   + f_0 * hl_528[k];

        t_349[k] = -2.0 * fl_214[k]
                   + f_0 * hl_529[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, fl_215, fl_216, fl_217, fl_218, \
                         fl_219, hl_530, hl_531, hl_532, hl_533, \
                         hl_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * fl_215[k]
                   + f_0 * hl_530[k];

        t_351[k] = -2.0 * fl_216[k]
                   + f_0 * hl_531[k];

        t_352[k] = -2.0 * fl_217[k]
                   + f_0 * hl_532[k];

        t_353[k] = -2.0 * fl_218[k]
                   + f_0 * hl_533[k];

        t_354[k] = -2.0 * fl_219[k]
                   + f_0 * hl_534[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, fl_220, fl_221, fl_222, fl_223, \
                         fl_224, hl_535, hl_536, hl_537, hl_538, \
                         hl_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * fl_220[k]
                   + f_0 * hl_535[k];

        t_356[k] = -2.0 * fl_221[k]
                   + f_0 * hl_536[k];

        t_357[k] = -2.0 * fl_222[k]
                   + f_0 * hl_537[k];

        t_358[k] = -2.0 * fl_223[k]
                   + f_0 * hl_538[k];

        t_359[k] = -2.0 * fl_224[k]
                   + f_0 * hl_539[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, fl_225, fl_226, fl_227, fl_228, \
                         fl_229, hl_540, hl_541, hl_542, hl_543, \
                         hl_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -fl_225[k]
                   + f_0 * hl_540[k];

        t_361[k] = -fl_226[k]
                   + f_0 * hl_541[k];

        t_362[k] = -fl_227[k]
                   + f_0 * hl_542[k];

        t_363[k] = -fl_228[k]
                   + f_0 * hl_543[k];

        t_364[k] = -fl_229[k]
                   + f_0 * hl_544[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, fl_230, fl_231, fl_232, fl_233, \
                         fl_234, hl_545, hl_546, hl_547, hl_548, \
                         hl_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -fl_230[k]
                   + f_0 * hl_545[k];

        t_366[k] = -fl_231[k]
                   + f_0 * hl_546[k];

        t_367[k] = -fl_232[k]
                   + f_0 * hl_547[k];

        t_368[k] = -fl_233[k]
                   + f_0 * hl_548[k];

        t_369[k] = -fl_234[k]
                   + f_0 * hl_549[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, fl_235, fl_236, fl_237, fl_238, \
                         fl_239, hl_550, hl_551, hl_552, hl_553, \
                         hl_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -fl_235[k]
                   + f_0 * hl_550[k];

        t_371[k] = -fl_236[k]
                   + f_0 * hl_551[k];

        t_372[k] = -fl_237[k]
                   + f_0 * hl_552[k];

        t_373[k] = -fl_238[k]
                   + f_0 * hl_553[k];

        t_374[k] = -fl_239[k]
                   + f_0 * hl_554[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, fl_240, fl_241, fl_242, fl_243, \
                         fl_244, hl_555, hl_556, hl_557, hl_558, \
                         hl_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -fl_240[k]
                   + f_0 * hl_555[k];

        t_376[k] = -fl_241[k]
                   + f_0 * hl_556[k];

        t_377[k] = -fl_242[k]
                   + f_0 * hl_557[k];

        t_378[k] = -fl_243[k]
                   + f_0 * hl_558[k];

        t_379[k] = -fl_244[k]
                   + f_0 * hl_559[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, fl_245, fl_246, fl_247, fl_248, \
                         fl_249, hl_560, hl_561, hl_562, hl_563, \
                         hl_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -fl_245[k]
                   + f_0 * hl_560[k];

        t_381[k] = -fl_246[k]
                   + f_0 * hl_561[k];

        t_382[k] = -fl_247[k]
                   + f_0 * hl_562[k];

        t_383[k] = -fl_248[k]
                   + f_0 * hl_563[k];

        t_384[k] = -fl_249[k]
                   + f_0 * hl_564[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, fl_250, fl_251, fl_252, fl_253, \
                         fl_254, hl_565, hl_566, hl_567, hl_568, \
                         hl_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -fl_250[k]
                   + f_0 * hl_565[k];

        t_386[k] = -fl_251[k]
                   + f_0 * hl_566[k];

        t_387[k] = -fl_252[k]
                   + f_0 * hl_567[k];

        t_388[k] = -fl_253[k]
                   + f_0 * hl_568[k];

        t_389[k] = -fl_254[k]
                   + f_0 * hl_569[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, fl_255, fl_256, fl_257, fl_258, \
                         fl_259, hl_570, hl_571, hl_572, hl_573, \
                         hl_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -fl_255[k]
                   + f_0 * hl_570[k];

        t_391[k] = -fl_256[k]
                   + f_0 * hl_571[k];

        t_392[k] = -fl_257[k]
                   + f_0 * hl_572[k];

        t_393[k] = -fl_258[k]
                   + f_0 * hl_573[k];

        t_394[k] = -fl_259[k]
                   + f_0 * hl_574[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, fl_260, fl_261, fl_262, fl_263, \
                         fl_264, hl_575, hl_576, hl_577, hl_578, \
                         hl_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -fl_260[k]
                   + f_0 * hl_575[k];

        t_396[k] = -fl_261[k]
                   + f_0 * hl_576[k];

        t_397[k] = -fl_262[k]
                   + f_0 * hl_577[k];

        t_398[k] = -fl_263[k]
                   + f_0 * hl_578[k];

        t_399[k] = -fl_264[k]
                   + f_0 * hl_579[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, fl_265, fl_266, fl_267, fl_268, \
                         fl_269, hl_580, hl_581, hl_582, hl_583, \
                         hl_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -fl_265[k]
                   + f_0 * hl_580[k];

        t_401[k] = -fl_266[k]
                   + f_0 * hl_581[k];

        t_402[k] = -fl_267[k]
                   + f_0 * hl_582[k];

        t_403[k] = -fl_268[k]
                   + f_0 * hl_583[k];

        t_404[k] = -fl_269[k]
                   + f_0 * hl_584[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, hl_585, \
                         hl_586, hl_587, hl_588, hl_589, hl_590, hl_591, \
                         hl_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * hl_585[k];

        t_406[k] = f_0 * hl_586[k];

        t_407[k] = f_0 * hl_587[k];

        t_408[k] = f_0 * hl_588[k];

        t_409[k] = f_0 * hl_589[k];

        t_410[k] = f_0 * hl_590[k];

        t_411[k] = f_0 * hl_591[k];

        t_412[k] = f_0 * hl_592[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, t_420, hl_593, \
                         hl_594, hl_595, hl_596, hl_597, hl_598, hl_599, \
                         hl_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * hl_593[k];

        t_414[k] = f_0 * hl_594[k];

        t_415[k] = f_0 * hl_595[k];

        t_416[k] = f_0 * hl_596[k];

        t_417[k] = f_0 * hl_597[k];

        t_418[k] = f_0 * hl_598[k];

        t_419[k] = f_0 * hl_599[k];

        t_420[k] = f_0 * hl_600[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, t_425, t_426, t_427, t_428, hl_601, \
                         hl_602, hl_603, hl_604, hl_605, hl_606, hl_607, \
                         hl_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_0 * hl_601[k];

        t_422[k] = f_0 * hl_602[k];

        t_423[k] = f_0 * hl_603[k];

        t_424[k] = f_0 * hl_604[k];

        t_425[k] = f_0 * hl_605[k];

        t_426[k] = f_0 * hl_606[k];

        t_427[k] = f_0 * hl_607[k];

        t_428[k] = f_0 * hl_608[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, t_434, t_435, t_436, hl_609, \
                         hl_610, hl_611, hl_612, hl_613, hl_614, hl_615, \
                         hl_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_0 * hl_609[k];

        t_430[k] = f_0 * hl_610[k];

        t_431[k] = f_0 * hl_611[k];

        t_432[k] = f_0 * hl_612[k];

        t_433[k] = f_0 * hl_613[k];

        t_434[k] = f_0 * hl_614[k];

        t_435[k] = f_0 * hl_615[k];

        t_436[k] = f_0 * hl_616[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, t_441, t_442, t_443, t_444, hl_617, \
                         hl_618, hl_619, hl_620, hl_621, hl_622, hl_623, \
                         hl_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_0 * hl_617[k];

        t_438[k] = f_0 * hl_618[k];

        t_439[k] = f_0 * hl_619[k];

        t_440[k] = f_0 * hl_620[k];

        t_441[k] = f_0 * hl_621[k];

        t_442[k] = f_0 * hl_622[k];

        t_443[k] = f_0 * hl_623[k];

        t_444[k] = f_0 * hl_624[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, t_450, t_451, fl_270, fl_271, \
                         hl_625, hl_626, hl_627, hl_628, hl_629, hl_675, \
                         hl_676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_0 * hl_625[k];

        t_446[k] = f_0 * hl_626[k];

        t_447[k] = f_0 * hl_627[k];

        t_448[k] = f_0 * hl_628[k];

        t_449[k] = f_0 * hl_629[k];

        t_450[k] = -4.0 * fl_270[k]
                   + f_0 * hl_675[k];

        t_451[k] = -4.0 * fl_271[k]
                   + f_0 * hl_676[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, fl_272, fl_273, fl_274, fl_275, \
                         fl_276, hl_677, hl_678, hl_679, hl_680, \
                         hl_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = -4.0 * fl_272[k]
                   + f_0 * hl_677[k];

        t_453[k] = -4.0 * fl_273[k]
                   + f_0 * hl_678[k];

        t_454[k] = -4.0 * fl_274[k]
                   + f_0 * hl_679[k];

        t_455[k] = -4.0 * fl_275[k]
                   + f_0 * hl_680[k];

        t_456[k] = -4.0 * fl_276[k]
                   + f_0 * hl_681[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, t_461, fl_277, fl_278, fl_279, fl_280, \
                         fl_281, hl_682, hl_683, hl_684, hl_685, \
                         hl_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = -4.0 * fl_277[k]
                   + f_0 * hl_682[k];

        t_458[k] = -4.0 * fl_278[k]
                   + f_0 * hl_683[k];

        t_459[k] = -4.0 * fl_279[k]
                   + f_0 * hl_684[k];

        t_460[k] = -4.0 * fl_280[k]
                   + f_0 * hl_685[k];

        t_461[k] = -4.0 * fl_281[k]
                   + f_0 * hl_686[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, fl_282, fl_283, fl_284, fl_285, \
                         fl_286, hl_687, hl_688, hl_689, hl_690, \
                         hl_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = -4.0 * fl_282[k]
                   + f_0 * hl_687[k];

        t_463[k] = -4.0 * fl_283[k]
                   + f_0 * hl_688[k];

        t_464[k] = -4.0 * fl_284[k]
                   + f_0 * hl_689[k];

        t_465[k] = -4.0 * fl_285[k]
                   + f_0 * hl_690[k];

        t_466[k] = -4.0 * fl_286[k]
                   + f_0 * hl_691[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, fl_287, fl_288, fl_289, fl_290, \
                         fl_291, hl_692, hl_693, hl_694, hl_695, \
                         hl_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = -4.0 * fl_287[k]
                   + f_0 * hl_692[k];

        t_468[k] = -4.0 * fl_288[k]
                   + f_0 * hl_693[k];

        t_469[k] = -4.0 * fl_289[k]
                   + f_0 * hl_694[k];

        t_470[k] = -4.0 * fl_290[k]
                   + f_0 * hl_695[k];

        t_471[k] = -4.0 * fl_291[k]
                   + f_0 * hl_696[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, fl_292, fl_293, fl_294, fl_295, \
                         fl_296, hl_697, hl_698, hl_699, hl_700, \
                         hl_701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = -4.0 * fl_292[k]
                   + f_0 * hl_697[k];

        t_473[k] = -4.0 * fl_293[k]
                   + f_0 * hl_698[k];

        t_474[k] = -4.0 * fl_294[k]
                   + f_0 * hl_699[k];

        t_475[k] = -4.0 * fl_295[k]
                   + f_0 * hl_700[k];

        t_476[k] = -4.0 * fl_296[k]
                   + f_0 * hl_701[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, fl_297, fl_298, fl_299, fl_300, \
                         fl_301, hl_702, hl_703, hl_704, hl_705, \
                         hl_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = -4.0 * fl_297[k]
                   + f_0 * hl_702[k];

        t_478[k] = -4.0 * fl_298[k]
                   + f_0 * hl_703[k];

        t_479[k] = -4.0 * fl_299[k]
                   + f_0 * hl_704[k];

        t_480[k] = -4.0 * fl_300[k]
                   + f_0 * hl_705[k];

        t_481[k] = -4.0 * fl_301[k]
                   + f_0 * hl_706[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, t_486, fl_302, fl_303, fl_304, fl_305, \
                         fl_306, hl_707, hl_708, hl_709, hl_710, \
                         hl_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = -4.0 * fl_302[k]
                   + f_0 * hl_707[k];

        t_483[k] = -4.0 * fl_303[k]
                   + f_0 * hl_708[k];

        t_484[k] = -4.0 * fl_304[k]
                   + f_0 * hl_709[k];

        t_485[k] = -4.0 * fl_305[k]
                   + f_0 * hl_710[k];

        t_486[k] = -4.0 * fl_306[k]
                   + f_0 * hl_711[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, fl_307, fl_308, fl_309, fl_310, \
                         fl_311, hl_712, hl_713, hl_714, hl_715, \
                         hl_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = -4.0 * fl_307[k]
                   + f_0 * hl_712[k];

        t_488[k] = -4.0 * fl_308[k]
                   + f_0 * hl_713[k];

        t_489[k] = -4.0 * fl_309[k]
                   + f_0 * hl_714[k];

        t_490[k] = -4.0 * fl_310[k]
                   + f_0 * hl_715[k];

        t_491[k] = -4.0 * fl_311[k]
                   + f_0 * hl_716[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, fl_312, fl_313, fl_314, fl_315, \
                         fl_316, hl_717, hl_718, hl_719, hl_720, \
                         hl_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = -4.0 * fl_312[k]
                   + f_0 * hl_717[k];

        t_493[k] = -4.0 * fl_313[k]
                   + f_0 * hl_718[k];

        t_494[k] = -4.0 * fl_314[k]
                   + f_0 * hl_719[k];

        t_495[k] = -3.0 * fl_315[k]
                   + f_0 * hl_720[k];

        t_496[k] = -3.0 * fl_316[k]
                   + f_0 * hl_721[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, t_500, t_501, fl_317, fl_318, fl_319, fl_320, \
                         fl_321, hl_722, hl_723, hl_724, hl_725, \
                         hl_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = -3.0 * fl_317[k]
                   + f_0 * hl_722[k];

        t_498[k] = -3.0 * fl_318[k]
                   + f_0 * hl_723[k];

        t_499[k] = -3.0 * fl_319[k]
                   + f_0 * hl_724[k];

        t_500[k] = -3.0 * fl_320[k]
                   + f_0 * hl_725[k];

        t_501[k] = -3.0 * fl_321[k]
                   + f_0 * hl_726[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, fl_322, fl_323, fl_324, fl_325, \
                         fl_326, hl_727, hl_728, hl_729, hl_730, \
                         hl_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = -3.0 * fl_322[k]
                   + f_0 * hl_727[k];

        t_503[k] = -3.0 * fl_323[k]
                   + f_0 * hl_728[k];

        t_504[k] = -3.0 * fl_324[k]
                   + f_0 * hl_729[k];

        t_505[k] = -3.0 * fl_325[k]
                   + f_0 * hl_730[k];

        t_506[k] = -3.0 * fl_326[k]
                   + f_0 * hl_731[k];
    }
}

static auto
compute_prim_geom_10_gl_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
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
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);

    const auto *fl_327 = buffer.data(fl + 327);
    const auto *fl_328 = buffer.data(fl + 328);
    const auto *fl_329 = buffer.data(fl + 329);
    const auto *fl_330 = buffer.data(fl + 330);
    const auto *fl_331 = buffer.data(fl + 331);
    const auto *fl_332 = buffer.data(fl + 332);
    const auto *fl_333 = buffer.data(fl + 333);
    const auto *fl_334 = buffer.data(fl + 334);
    const auto *fl_335 = buffer.data(fl + 335);
    const auto *fl_336 = buffer.data(fl + 336);
    const auto *fl_337 = buffer.data(fl + 337);
    const auto *fl_338 = buffer.data(fl + 338);
    const auto *fl_339 = buffer.data(fl + 339);
    const auto *fl_340 = buffer.data(fl + 340);
    const auto *fl_341 = buffer.data(fl + 341);
    const auto *fl_342 = buffer.data(fl + 342);
    const auto *fl_343 = buffer.data(fl + 343);
    const auto *fl_344 = buffer.data(fl + 344);
    const auto *fl_345 = buffer.data(fl + 345);
    const auto *fl_346 = buffer.data(fl + 346);
    const auto *fl_347 = buffer.data(fl + 347);
    const auto *fl_348 = buffer.data(fl + 348);
    const auto *fl_349 = buffer.data(fl + 349);
    const auto *fl_350 = buffer.data(fl + 350);
    const auto *fl_351 = buffer.data(fl + 351);
    const auto *fl_352 = buffer.data(fl + 352);
    const auto *fl_353 = buffer.data(fl + 353);
    const auto *fl_354 = buffer.data(fl + 354);
    const auto *fl_355 = buffer.data(fl + 355);
    const auto *fl_356 = buffer.data(fl + 356);
    const auto *fl_357 = buffer.data(fl + 357);
    const auto *fl_358 = buffer.data(fl + 358);
    const auto *fl_359 = buffer.data(fl + 359);
    const auto *fl_360 = buffer.data(fl + 360);
    const auto *fl_361 = buffer.data(fl + 361);
    const auto *fl_362 = buffer.data(fl + 362);
    const auto *fl_363 = buffer.data(fl + 363);
    const auto *fl_364 = buffer.data(fl + 364);
    const auto *fl_365 = buffer.data(fl + 365);
    const auto *fl_366 = buffer.data(fl + 366);
    const auto *fl_367 = buffer.data(fl + 367);
    const auto *fl_368 = buffer.data(fl + 368);
    const auto *fl_369 = buffer.data(fl + 369);
    const auto *fl_370 = buffer.data(fl + 370);
    const auto *fl_371 = buffer.data(fl + 371);
    const auto *fl_372 = buffer.data(fl + 372);
    const auto *fl_373 = buffer.data(fl + 373);
    const auto *fl_374 = buffer.data(fl + 374);
    const auto *fl_375 = buffer.data(fl + 375);
    const auto *fl_376 = buffer.data(fl + 376);
    const auto *fl_377 = buffer.data(fl + 377);
    const auto *fl_378 = buffer.data(fl + 378);
    const auto *fl_379 = buffer.data(fl + 379);
    const auto *fl_380 = buffer.data(fl + 380);
    const auto *fl_381 = buffer.data(fl + 381);
    const auto *fl_382 = buffer.data(fl + 382);
    const auto *fl_383 = buffer.data(fl + 383);
    const auto *fl_384 = buffer.data(fl + 384);
    const auto *fl_385 = buffer.data(fl + 385);
    const auto *fl_386 = buffer.data(fl + 386);
    const auto *fl_387 = buffer.data(fl + 387);
    const auto *fl_388 = buffer.data(fl + 388);
    const auto *fl_389 = buffer.data(fl + 389);
    const auto *fl_390 = buffer.data(fl + 390);
    const auto *fl_391 = buffer.data(fl + 391);
    const auto *fl_392 = buffer.data(fl + 392);
    const auto *fl_393 = buffer.data(fl + 393);
    const auto *fl_394 = buffer.data(fl + 394);
    const auto *fl_395 = buffer.data(fl + 395);
    const auto *fl_396 = buffer.data(fl + 396);
    const auto *fl_397 = buffer.data(fl + 397);
    const auto *fl_398 = buffer.data(fl + 398);
    const auto *fl_399 = buffer.data(fl + 399);
    const auto *fl_400 = buffer.data(fl + 400);
    const auto *fl_401 = buffer.data(fl + 401);
    const auto *fl_402 = buffer.data(fl + 402);
    const auto *fl_403 = buffer.data(fl + 403);
    const auto *fl_404 = buffer.data(fl + 404);
    const auto *fl_405 = buffer.data(fl + 405);
    const auto *fl_406 = buffer.data(fl + 406);
    const auto *fl_407 = buffer.data(fl + 407);
    const auto *fl_408 = buffer.data(fl + 408);
    const auto *fl_409 = buffer.data(fl + 409);
    const auto *fl_410 = buffer.data(fl + 410);
    const auto *fl_411 = buffer.data(fl + 411);
    const auto *fl_412 = buffer.data(fl + 412);
    const auto *fl_413 = buffer.data(fl + 413);
    const auto *fl_414 = buffer.data(fl + 414);
    const auto *fl_415 = buffer.data(fl + 415);
    const auto *fl_416 = buffer.data(fl + 416);
    const auto *fl_417 = buffer.data(fl + 417);
    const auto *fl_418 = buffer.data(fl + 418);
    const auto *fl_419 = buffer.data(fl + 419);
    const auto *fl_420 = buffer.data(fl + 420);
    const auto *fl_421 = buffer.data(fl + 421);
    const auto *fl_422 = buffer.data(fl + 422);
    const auto *fl_423 = buffer.data(fl + 423);
    const auto *fl_424 = buffer.data(fl + 424);
    const auto *fl_425 = buffer.data(fl + 425);
    const auto *fl_426 = buffer.data(fl + 426);
    const auto *fl_427 = buffer.data(fl + 427);
    const auto *fl_428 = buffer.data(fl + 428);
    const auto *fl_429 = buffer.data(fl + 429);
    const auto *fl_430 = buffer.data(fl + 430);
    const auto *fl_431 = buffer.data(fl + 431);
    const auto *fl_432 = buffer.data(fl + 432);
    const auto *fl_433 = buffer.data(fl + 433);
    const auto *fl_434 = buffer.data(fl + 434);
    const auto *fl_435 = buffer.data(fl + 435);
    const auto *fl_436 = buffer.data(fl + 436);
    const auto *fl_437 = buffer.data(fl + 437);
    const auto *fl_438 = buffer.data(fl + 438);
    const auto *fl_439 = buffer.data(fl + 439);
    const auto *fl_440 = buffer.data(fl + 440);
    const auto *fl_441 = buffer.data(fl + 441);
    const auto *fl_442 = buffer.data(fl + 442);
    const auto *fl_443 = buffer.data(fl + 443);
    const auto *fl_444 = buffer.data(fl + 444);
    const auto *fl_445 = buffer.data(fl + 445);
    const auto *fl_446 = buffer.data(fl + 446);
    const auto *fl_447 = buffer.data(fl + 447);
    const auto *fl_448 = buffer.data(fl + 448);
    const auto *fl_449 = buffer.data(fl + 449);

    const auto *hl_732 = buffer.data(hl + 732);
    const auto *hl_733 = buffer.data(hl + 733);
    const auto *hl_734 = buffer.data(hl + 734);
    const auto *hl_735 = buffer.data(hl + 735);
    const auto *hl_736 = buffer.data(hl + 736);
    const auto *hl_737 = buffer.data(hl + 737);
    const auto *hl_738 = buffer.data(hl + 738);
    const auto *hl_739 = buffer.data(hl + 739);
    const auto *hl_740 = buffer.data(hl + 740);
    const auto *hl_741 = buffer.data(hl + 741);
    const auto *hl_742 = buffer.data(hl + 742);
    const auto *hl_743 = buffer.data(hl + 743);
    const auto *hl_744 = buffer.data(hl + 744);
    const auto *hl_745 = buffer.data(hl + 745);
    const auto *hl_746 = buffer.data(hl + 746);
    const auto *hl_747 = buffer.data(hl + 747);
    const auto *hl_748 = buffer.data(hl + 748);
    const auto *hl_749 = buffer.data(hl + 749);
    const auto *hl_750 = buffer.data(hl + 750);
    const auto *hl_751 = buffer.data(hl + 751);
    const auto *hl_752 = buffer.data(hl + 752);
    const auto *hl_753 = buffer.data(hl + 753);
    const auto *hl_754 = buffer.data(hl + 754);
    const auto *hl_755 = buffer.data(hl + 755);
    const auto *hl_756 = buffer.data(hl + 756);
    const auto *hl_757 = buffer.data(hl + 757);
    const auto *hl_758 = buffer.data(hl + 758);
    const auto *hl_759 = buffer.data(hl + 759);
    const auto *hl_760 = buffer.data(hl + 760);
    const auto *hl_761 = buffer.data(hl + 761);
    const auto *hl_762 = buffer.data(hl + 762);
    const auto *hl_763 = buffer.data(hl + 763);
    const auto *hl_764 = buffer.data(hl + 764);
    const auto *hl_765 = buffer.data(hl + 765);
    const auto *hl_766 = buffer.data(hl + 766);
    const auto *hl_767 = buffer.data(hl + 767);
    const auto *hl_768 = buffer.data(hl + 768);
    const auto *hl_769 = buffer.data(hl + 769);
    const auto *hl_770 = buffer.data(hl + 770);
    const auto *hl_771 = buffer.data(hl + 771);
    const auto *hl_772 = buffer.data(hl + 772);
    const auto *hl_773 = buffer.data(hl + 773);
    const auto *hl_774 = buffer.data(hl + 774);
    const auto *hl_775 = buffer.data(hl + 775);
    const auto *hl_776 = buffer.data(hl + 776);
    const auto *hl_777 = buffer.data(hl + 777);
    const auto *hl_778 = buffer.data(hl + 778);
    const auto *hl_779 = buffer.data(hl + 779);
    const auto *hl_780 = buffer.data(hl + 780);
    const auto *hl_781 = buffer.data(hl + 781);
    const auto *hl_782 = buffer.data(hl + 782);
    const auto *hl_783 = buffer.data(hl + 783);
    const auto *hl_784 = buffer.data(hl + 784);
    const auto *hl_785 = buffer.data(hl + 785);
    const auto *hl_786 = buffer.data(hl + 786);
    const auto *hl_787 = buffer.data(hl + 787);
    const auto *hl_788 = buffer.data(hl + 788);
    const auto *hl_789 = buffer.data(hl + 789);
    const auto *hl_790 = buffer.data(hl + 790);
    const auto *hl_791 = buffer.data(hl + 791);
    const auto *hl_792 = buffer.data(hl + 792);
    const auto *hl_793 = buffer.data(hl + 793);
    const auto *hl_794 = buffer.data(hl + 794);
    const auto *hl_795 = buffer.data(hl + 795);
    const auto *hl_796 = buffer.data(hl + 796);
    const auto *hl_797 = buffer.data(hl + 797);
    const auto *hl_798 = buffer.data(hl + 798);
    const auto *hl_799 = buffer.data(hl + 799);
    const auto *hl_800 = buffer.data(hl + 800);
    const auto *hl_801 = buffer.data(hl + 801);
    const auto *hl_802 = buffer.data(hl + 802);
    const auto *hl_803 = buffer.data(hl + 803);
    const auto *hl_804 = buffer.data(hl + 804);
    const auto *hl_805 = buffer.data(hl + 805);
    const auto *hl_806 = buffer.data(hl + 806);
    const auto *hl_807 = buffer.data(hl + 807);
    const auto *hl_808 = buffer.data(hl + 808);
    const auto *hl_809 = buffer.data(hl + 809);
    const auto *hl_810 = buffer.data(hl + 810);
    const auto *hl_811 = buffer.data(hl + 811);
    const auto *hl_812 = buffer.data(hl + 812);
    const auto *hl_813 = buffer.data(hl + 813);
    const auto *hl_814 = buffer.data(hl + 814);
    const auto *hl_815 = buffer.data(hl + 815);
    const auto *hl_816 = buffer.data(hl + 816);
    const auto *hl_817 = buffer.data(hl + 817);
    const auto *hl_818 = buffer.data(hl + 818);
    const auto *hl_819 = buffer.data(hl + 819);
    const auto *hl_820 = buffer.data(hl + 820);
    const auto *hl_821 = buffer.data(hl + 821);
    const auto *hl_822 = buffer.data(hl + 822);
    const auto *hl_823 = buffer.data(hl + 823);
    const auto *hl_824 = buffer.data(hl + 824);
    const auto *hl_825 = buffer.data(hl + 825);
    const auto *hl_826 = buffer.data(hl + 826);
    const auto *hl_827 = buffer.data(hl + 827);
    const auto *hl_828 = buffer.data(hl + 828);
    const auto *hl_829 = buffer.data(hl + 829);
    const auto *hl_830 = buffer.data(hl + 830);
    const auto *hl_831 = buffer.data(hl + 831);
    const auto *hl_832 = buffer.data(hl + 832);
    const auto *hl_833 = buffer.data(hl + 833);
    const auto *hl_834 = buffer.data(hl + 834);
    const auto *hl_835 = buffer.data(hl + 835);
    const auto *hl_836 = buffer.data(hl + 836);
    const auto *hl_837 = buffer.data(hl + 837);
    const auto *hl_838 = buffer.data(hl + 838);
    const auto *hl_839 = buffer.data(hl + 839);
    const auto *hl_840 = buffer.data(hl + 840);
    const auto *hl_841 = buffer.data(hl + 841);
    const auto *hl_842 = buffer.data(hl + 842);
    const auto *hl_843 = buffer.data(hl + 843);
    const auto *hl_844 = buffer.data(hl + 844);
    const auto *hl_845 = buffer.data(hl + 845);
    const auto *hl_846 = buffer.data(hl + 846);
    const auto *hl_847 = buffer.data(hl + 847);
    const auto *hl_848 = buffer.data(hl + 848);
    const auto *hl_849 = buffer.data(hl + 849);
    const auto *hl_850 = buffer.data(hl + 850);
    const auto *hl_851 = buffer.data(hl + 851);
    const auto *hl_852 = buffer.data(hl + 852);
    const auto *hl_853 = buffer.data(hl + 853);
    const auto *hl_854 = buffer.data(hl + 854);
    const auto *hl_855 = buffer.data(hl + 855);
    const auto *hl_856 = buffer.data(hl + 856);
    const auto *hl_857 = buffer.data(hl + 857);
    const auto *hl_858 = buffer.data(hl + 858);
    const auto *hl_859 = buffer.data(hl + 859);
    const auto *hl_860 = buffer.data(hl + 860);
    const auto *hl_861 = buffer.data(hl + 861);
    const auto *hl_862 = buffer.data(hl + 862);
    const auto *hl_863 = buffer.data(hl + 863);
    const auto *hl_864 = buffer.data(hl + 864);
    const auto *hl_865 = buffer.data(hl + 865);
    const auto *hl_866 = buffer.data(hl + 866);
    const auto *hl_867 = buffer.data(hl + 867);
    const auto *hl_868 = buffer.data(hl + 868);
    const auto *hl_869 = buffer.data(hl + 869);
    const auto *hl_870 = buffer.data(hl + 870);
    const auto *hl_871 = buffer.data(hl + 871);
    const auto *hl_872 = buffer.data(hl + 872);
    const auto *hl_873 = buffer.data(hl + 873);
    const auto *hl_874 = buffer.data(hl + 874);
    const auto *hl_875 = buffer.data(hl + 875);
    const auto *hl_876 = buffer.data(hl + 876);
    const auto *hl_877 = buffer.data(hl + 877);
    const auto *hl_878 = buffer.data(hl + 878);
    const auto *hl_879 = buffer.data(hl + 879);
    const auto *hl_880 = buffer.data(hl + 880);
    const auto *hl_881 = buffer.data(hl + 881);
    const auto *hl_882 = buffer.data(hl + 882);
    const auto *hl_883 = buffer.data(hl + 883);
    const auto *hl_884 = buffer.data(hl + 884);
    const auto *hl_885 = buffer.data(hl + 885);
    const auto *hl_886 = buffer.data(hl + 886);
    const auto *hl_887 = buffer.data(hl + 887);
    const auto *hl_888 = buffer.data(hl + 888);
    const auto *hl_889 = buffer.data(hl + 889);
    const auto *hl_890 = buffer.data(hl + 890);
    const auto *hl_891 = buffer.data(hl + 891);
    const auto *hl_892 = buffer.data(hl + 892);
    const auto *hl_893 = buffer.data(hl + 893);
    const auto *hl_894 = buffer.data(hl + 894);
    const auto *hl_895 = buffer.data(hl + 895);
    const auto *hl_896 = buffer.data(hl + 896);
    const auto *hl_897 = buffer.data(hl + 897);

#pragma omp simd aligned(t_507, t_508, t_509, t_510, t_511, fl_327, fl_328, fl_329, fl_330, \
                         fl_331, hl_732, hl_733, hl_734, hl_735, \
                         hl_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = -3.0 * fl_327[k]
                   + f_0 * hl_732[k];

        t_508[k] = -3.0 * fl_328[k]
                   + f_0 * hl_733[k];

        t_509[k] = -3.0 * fl_329[k]
                   + f_0 * hl_734[k];

        t_510[k] = -3.0 * fl_330[k]
                   + f_0 * hl_735[k];

        t_511[k] = -3.0 * fl_331[k]
                   + f_0 * hl_736[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, fl_332, fl_333, fl_334, fl_335, \
                         fl_336, hl_737, hl_738, hl_739, hl_740, \
                         hl_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = -3.0 * fl_332[k]
                   + f_0 * hl_737[k];

        t_513[k] = -3.0 * fl_333[k]
                   + f_0 * hl_738[k];

        t_514[k] = -3.0 * fl_334[k]
                   + f_0 * hl_739[k];

        t_515[k] = -3.0 * fl_335[k]
                   + f_0 * hl_740[k];

        t_516[k] = -3.0 * fl_336[k]
                   + f_0 * hl_741[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, fl_337, fl_338, fl_339, fl_340, \
                         fl_341, hl_742, hl_743, hl_744, hl_745, \
                         hl_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = -3.0 * fl_337[k]
                   + f_0 * hl_742[k];

        t_518[k] = -3.0 * fl_338[k]
                   + f_0 * hl_743[k];

        t_519[k] = -3.0 * fl_339[k]
                   + f_0 * hl_744[k];

        t_520[k] = -3.0 * fl_340[k]
                   + f_0 * hl_745[k];

        t_521[k] = -3.0 * fl_341[k]
                   + f_0 * hl_746[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, fl_342, fl_343, fl_344, fl_345, \
                         fl_346, hl_747, hl_748, hl_749, hl_750, \
                         hl_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = -3.0 * fl_342[k]
                   + f_0 * hl_747[k];

        t_523[k] = -3.0 * fl_343[k]
                   + f_0 * hl_748[k];

        t_524[k] = -3.0 * fl_344[k]
                   + f_0 * hl_749[k];

        t_525[k] = -3.0 * fl_345[k]
                   + f_0 * hl_750[k];

        t_526[k] = -3.0 * fl_346[k]
                   + f_0 * hl_751[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, fl_347, fl_348, fl_349, fl_350, \
                         fl_351, hl_752, hl_753, hl_754, hl_755, \
                         hl_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = -3.0 * fl_347[k]
                   + f_0 * hl_752[k];

        t_528[k] = -3.0 * fl_348[k]
                   + f_0 * hl_753[k];

        t_529[k] = -3.0 * fl_349[k]
                   + f_0 * hl_754[k];

        t_530[k] = -3.0 * fl_350[k]
                   + f_0 * hl_755[k];

        t_531[k] = -3.0 * fl_351[k]
                   + f_0 * hl_756[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, fl_352, fl_353, fl_354, fl_355, \
                         fl_356, hl_757, hl_758, hl_759, hl_760, \
                         hl_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = -3.0 * fl_352[k]
                   + f_0 * hl_757[k];

        t_533[k] = -3.0 * fl_353[k]
                   + f_0 * hl_758[k];

        t_534[k] = -3.0 * fl_354[k]
                   + f_0 * hl_759[k];

        t_535[k] = -3.0 * fl_355[k]
                   + f_0 * hl_760[k];

        t_536[k] = -3.0 * fl_356[k]
                   + f_0 * hl_761[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, t_541, fl_357, fl_358, fl_359, fl_360, \
                         fl_361, hl_762, hl_763, hl_764, hl_765, \
                         hl_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = -3.0 * fl_357[k]
                   + f_0 * hl_762[k];

        t_538[k] = -3.0 * fl_358[k]
                   + f_0 * hl_763[k];

        t_539[k] = -3.0 * fl_359[k]
                   + f_0 * hl_764[k];

        t_540[k] = -2.0 * fl_360[k]
                   + f_0 * hl_765[k];

        t_541[k] = -2.0 * fl_361[k]
                   + f_0 * hl_766[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, fl_362, fl_363, fl_364, fl_365, \
                         fl_366, hl_767, hl_768, hl_769, hl_770, \
                         hl_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = -2.0 * fl_362[k]
                   + f_0 * hl_767[k];

        t_543[k] = -2.0 * fl_363[k]
                   + f_0 * hl_768[k];

        t_544[k] = -2.0 * fl_364[k]
                   + f_0 * hl_769[k];

        t_545[k] = -2.0 * fl_365[k]
                   + f_0 * hl_770[k];

        t_546[k] = -2.0 * fl_366[k]
                   + f_0 * hl_771[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, fl_367, fl_368, fl_369, fl_370, \
                         fl_371, hl_772, hl_773, hl_774, hl_775, \
                         hl_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = -2.0 * fl_367[k]
                   + f_0 * hl_772[k];

        t_548[k] = -2.0 * fl_368[k]
                   + f_0 * hl_773[k];

        t_549[k] = -2.0 * fl_369[k]
                   + f_0 * hl_774[k];

        t_550[k] = -2.0 * fl_370[k]
                   + f_0 * hl_775[k];

        t_551[k] = -2.0 * fl_371[k]
                   + f_0 * hl_776[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, t_556, fl_372, fl_373, fl_374, fl_375, \
                         fl_376, hl_777, hl_778, hl_779, hl_780, \
                         hl_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = -2.0 * fl_372[k]
                   + f_0 * hl_777[k];

        t_553[k] = -2.0 * fl_373[k]
                   + f_0 * hl_778[k];

        t_554[k] = -2.0 * fl_374[k]
                   + f_0 * hl_779[k];

        t_555[k] = -2.0 * fl_375[k]
                   + f_0 * hl_780[k];

        t_556[k] = -2.0 * fl_376[k]
                   + f_0 * hl_781[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, fl_377, fl_378, fl_379, fl_380, \
                         fl_381, hl_782, hl_783, hl_784, hl_785, \
                         hl_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = -2.0 * fl_377[k]
                   + f_0 * hl_782[k];

        t_558[k] = -2.0 * fl_378[k]
                   + f_0 * hl_783[k];

        t_559[k] = -2.0 * fl_379[k]
                   + f_0 * hl_784[k];

        t_560[k] = -2.0 * fl_380[k]
                   + f_0 * hl_785[k];

        t_561[k] = -2.0 * fl_381[k]
                   + f_0 * hl_786[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, fl_382, fl_383, fl_384, fl_385, \
                         fl_386, hl_787, hl_788, hl_789, hl_790, \
                         hl_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = -2.0 * fl_382[k]
                   + f_0 * hl_787[k];

        t_563[k] = -2.0 * fl_383[k]
                   + f_0 * hl_788[k];

        t_564[k] = -2.0 * fl_384[k]
                   + f_0 * hl_789[k];

        t_565[k] = -2.0 * fl_385[k]
                   + f_0 * hl_790[k];

        t_566[k] = -2.0 * fl_386[k]
                   + f_0 * hl_791[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, fl_387, fl_388, fl_389, fl_390, \
                         fl_391, hl_792, hl_793, hl_794, hl_795, \
                         hl_796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -2.0 * fl_387[k]
                   + f_0 * hl_792[k];

        t_568[k] = -2.0 * fl_388[k]
                   + f_0 * hl_793[k];

        t_569[k] = -2.0 * fl_389[k]
                   + f_0 * hl_794[k];

        t_570[k] = -2.0 * fl_390[k]
                   + f_0 * hl_795[k];

        t_571[k] = -2.0 * fl_391[k]
                   + f_0 * hl_796[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, fl_392, fl_393, fl_394, fl_395, \
                         fl_396, hl_797, hl_798, hl_799, hl_800, \
                         hl_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = -2.0 * fl_392[k]
                   + f_0 * hl_797[k];

        t_573[k] = -2.0 * fl_393[k]
                   + f_0 * hl_798[k];

        t_574[k] = -2.0 * fl_394[k]
                   + f_0 * hl_799[k];

        t_575[k] = -2.0 * fl_395[k]
                   + f_0 * hl_800[k];

        t_576[k] = -2.0 * fl_396[k]
                   + f_0 * hl_801[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, fl_397, fl_398, fl_399, fl_400, \
                         fl_401, hl_802, hl_803, hl_804, hl_805, \
                         hl_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = -2.0 * fl_397[k]
                   + f_0 * hl_802[k];

        t_578[k] = -2.0 * fl_398[k]
                   + f_0 * hl_803[k];

        t_579[k] = -2.0 * fl_399[k]
                   + f_0 * hl_804[k];

        t_580[k] = -2.0 * fl_400[k]
                   + f_0 * hl_805[k];

        t_581[k] = -2.0 * fl_401[k]
                   + f_0 * hl_806[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, fl_402, fl_403, fl_404, fl_405, \
                         fl_406, hl_807, hl_808, hl_809, hl_810, \
                         hl_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = -2.0 * fl_402[k]
                   + f_0 * hl_807[k];

        t_583[k] = -2.0 * fl_403[k]
                   + f_0 * hl_808[k];

        t_584[k] = -2.0 * fl_404[k]
                   + f_0 * hl_809[k];

        t_585[k] = -fl_405[k]
                   + f_0 * hl_810[k];

        t_586[k] = -fl_406[k]
                   + f_0 * hl_811[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, t_591, fl_407, fl_408, fl_409, fl_410, \
                         fl_411, hl_812, hl_813, hl_814, hl_815, \
                         hl_816 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = -fl_407[k]
                   + f_0 * hl_812[k];

        t_588[k] = -fl_408[k]
                   + f_0 * hl_813[k];

        t_589[k] = -fl_409[k]
                   + f_0 * hl_814[k];

        t_590[k] = -fl_410[k]
                   + f_0 * hl_815[k];

        t_591[k] = -fl_411[k]
                   + f_0 * hl_816[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, t_596, fl_412, fl_413, fl_414, fl_415, \
                         fl_416, hl_817, hl_818, hl_819, hl_820, \
                         hl_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = -fl_412[k]
                   + f_0 * hl_817[k];

        t_593[k] = -fl_413[k]
                   + f_0 * hl_818[k];

        t_594[k] = -fl_414[k]
                   + f_0 * hl_819[k];

        t_595[k] = -fl_415[k]
                   + f_0 * hl_820[k];

        t_596[k] = -fl_416[k]
                   + f_0 * hl_821[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, t_600, t_601, fl_417, fl_418, fl_419, fl_420, \
                         fl_421, hl_822, hl_823, hl_824, hl_825, \
                         hl_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = -fl_417[k]
                   + f_0 * hl_822[k];

        t_598[k] = -fl_418[k]
                   + f_0 * hl_823[k];

        t_599[k] = -fl_419[k]
                   + f_0 * hl_824[k];

        t_600[k] = -fl_420[k]
                   + f_0 * hl_825[k];

        t_601[k] = -fl_421[k]
                   + f_0 * hl_826[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, fl_422, fl_423, fl_424, fl_425, \
                         fl_426, hl_827, hl_828, hl_829, hl_830, \
                         hl_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = -fl_422[k]
                   + f_0 * hl_827[k];

        t_603[k] = -fl_423[k]
                   + f_0 * hl_828[k];

        t_604[k] = -fl_424[k]
                   + f_0 * hl_829[k];

        t_605[k] = -fl_425[k]
                   + f_0 * hl_830[k];

        t_606[k] = -fl_426[k]
                   + f_0 * hl_831[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, t_611, fl_427, fl_428, fl_429, fl_430, \
                         fl_431, hl_832, hl_833, hl_834, hl_835, \
                         hl_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = -fl_427[k]
                   + f_0 * hl_832[k];

        t_608[k] = -fl_428[k]
                   + f_0 * hl_833[k];

        t_609[k] = -fl_429[k]
                   + f_0 * hl_834[k];

        t_610[k] = -fl_430[k]
                   + f_0 * hl_835[k];

        t_611[k] = -fl_431[k]
                   + f_0 * hl_836[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, fl_432, fl_433, fl_434, fl_435, \
                         fl_436, hl_837, hl_838, hl_839, hl_840, \
                         hl_841 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = -fl_432[k]
                   + f_0 * hl_837[k];

        t_613[k] = -fl_433[k]
                   + f_0 * hl_838[k];

        t_614[k] = -fl_434[k]
                   + f_0 * hl_839[k];

        t_615[k] = -fl_435[k]
                   + f_0 * hl_840[k];

        t_616[k] = -fl_436[k]
                   + f_0 * hl_841[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, fl_437, fl_438, fl_439, fl_440, \
                         fl_441, hl_842, hl_843, hl_844, hl_845, \
                         hl_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = -fl_437[k]
                   + f_0 * hl_842[k];

        t_618[k] = -fl_438[k]
                   + f_0 * hl_843[k];

        t_619[k] = -fl_439[k]
                   + f_0 * hl_844[k];

        t_620[k] = -fl_440[k]
                   + f_0 * hl_845[k];

        t_621[k] = -fl_441[k]
                   + f_0 * hl_846[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, fl_442, fl_443, fl_444, fl_445, \
                         fl_446, hl_847, hl_848, hl_849, hl_850, \
                         hl_851 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = -fl_442[k]
                   + f_0 * hl_847[k];

        t_623[k] = -fl_443[k]
                   + f_0 * hl_848[k];

        t_624[k] = -fl_444[k]
                   + f_0 * hl_849[k];

        t_625[k] = -fl_445[k]
                   + f_0 * hl_850[k];

        t_626[k] = -fl_446[k]
                   + f_0 * hl_851[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, t_631, t_632, fl_447, fl_448, fl_449, \
                         hl_852, hl_853, hl_854, hl_855, hl_856, \
                         hl_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = -fl_447[k]
                   + f_0 * hl_852[k];

        t_628[k] = -fl_448[k]
                   + f_0 * hl_853[k];

        t_629[k] = -fl_449[k]
                   + f_0 * hl_854[k];

        t_630[k] = f_0 * hl_855[k];

        t_631[k] = f_0 * hl_856[k];

        t_632[k] = f_0 * hl_857[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, t_636, t_637, t_638, t_639, t_640, hl_858, \
                         hl_859, hl_860, hl_861, hl_862, hl_863, hl_864, \
                         hl_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_0 * hl_858[k];

        t_634[k] = f_0 * hl_859[k];

        t_635[k] = f_0 * hl_860[k];

        t_636[k] = f_0 * hl_861[k];

        t_637[k] = f_0 * hl_862[k];

        t_638[k] = f_0 * hl_863[k];

        t_639[k] = f_0 * hl_864[k];

        t_640[k] = f_0 * hl_865[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, t_645, t_646, t_647, t_648, hl_866, \
                         hl_867, hl_868, hl_869, hl_870, hl_871, hl_872, \
                         hl_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_0 * hl_866[k];

        t_642[k] = f_0 * hl_867[k];

        t_643[k] = f_0 * hl_868[k];

        t_644[k] = f_0 * hl_869[k];

        t_645[k] = f_0 * hl_870[k];

        t_646[k] = f_0 * hl_871[k];

        t_647[k] = f_0 * hl_872[k];

        t_648[k] = f_0 * hl_873[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, t_654, t_655, t_656, hl_874, \
                         hl_875, hl_876, hl_877, hl_878, hl_879, hl_880, \
                         hl_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_0 * hl_874[k];

        t_650[k] = f_0 * hl_875[k];

        t_651[k] = f_0 * hl_876[k];

        t_652[k] = f_0 * hl_877[k];

        t_653[k] = f_0 * hl_878[k];

        t_654[k] = f_0 * hl_879[k];

        t_655[k] = f_0 * hl_880[k];

        t_656[k] = f_0 * hl_881[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, t_662, t_663, t_664, hl_882, \
                         hl_883, hl_884, hl_885, hl_886, hl_887, hl_888, \
                         hl_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_0 * hl_882[k];

        t_658[k] = f_0 * hl_883[k];

        t_659[k] = f_0 * hl_884[k];

        t_660[k] = f_0 * hl_885[k];

        t_661[k] = f_0 * hl_886[k];

        t_662[k] = f_0 * hl_887[k];

        t_663[k] = f_0 * hl_888[k];

        t_664[k] = f_0 * hl_889[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, t_670, t_671, t_672, hl_890, \
                         hl_891, hl_892, hl_893, hl_894, hl_895, hl_896, \
                         hl_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_0 * hl_890[k];

        t_666[k] = f_0 * hl_891[k];

        t_667[k] = f_0 * hl_892[k];

        t_668[k] = f_0 * hl_893[k];

        t_669[k] = f_0 * hl_894[k];

        t_670[k] = f_0 * hl_895[k];

        t_671[k] = f_0 * hl_896[k];

        t_672[k] = f_0 * hl_897[k];
    }
}

static auto
compute_prim_geom_10_gl_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);

    const auto *hl_898 = buffer.data(hl + 898);
    const auto *hl_899 = buffer.data(hl + 899);

#pragma omp simd aligned(t_673, t_674, hl_898, hl_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_0 * hl_898[k];

        t_674[k] = f_0 * hl_899[k];
    }
}

auto
compute_prim_geom_10_gl_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t fl, const size_t hl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gl_electron_repulsion_1_piece0(buffer, target, fl, hl, ncols, alpha);

    compute_prim_geom_10_gl_electron_repulsion_1_piece1(buffer, target, fl, hl, ncols, alpha);

    compute_prim_geom_10_gl_electron_repulsion_1_piece2(buffer, target, fl, hl, ncols, alpha);

    compute_prim_geom_10_gl_electron_repulsion_1_piece3(buffer, target, fl, hl, ncols, alpha);

    compute_prim_geom_10_gl_electron_repulsion_1_piece4(buffer, target, hl, ncols, alpha);
}

static auto
compute_prim_geom_10_gl_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
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

    const auto *fl_0 = buffer.data(fl + 0);
    const auto *fl_1 = buffer.data(fl + 1);
    const auto *fl_2 = buffer.data(fl + 2);
    const auto *fl_3 = buffer.data(fl + 3);
    const auto *fl_4 = buffer.data(fl + 4);
    const auto *fl_5 = buffer.data(fl + 5);
    const auto *fl_6 = buffer.data(fl + 6);
    const auto *fl_7 = buffer.data(fl + 7);
    const auto *fl_8 = buffer.data(fl + 8);
    const auto *fl_9 = buffer.data(fl + 9);
    const auto *fl_10 = buffer.data(fl + 10);
    const auto *fl_11 = buffer.data(fl + 11);
    const auto *fl_12 = buffer.data(fl + 12);
    const auto *fl_13 = buffer.data(fl + 13);
    const auto *fl_14 = buffer.data(fl + 14);
    const auto *fl_15 = buffer.data(fl + 15);
    const auto *fl_16 = buffer.data(fl + 16);
    const auto *fl_17 = buffer.data(fl + 17);
    const auto *fl_18 = buffer.data(fl + 18);
    const auto *fl_19 = buffer.data(fl + 19);
    const auto *fl_20 = buffer.data(fl + 20);
    const auto *fl_21 = buffer.data(fl + 21);
    const auto *fl_22 = buffer.data(fl + 22);
    const auto *fl_23 = buffer.data(fl + 23);
    const auto *fl_24 = buffer.data(fl + 24);
    const auto *fl_25 = buffer.data(fl + 25);
    const auto *fl_26 = buffer.data(fl + 26);
    const auto *fl_27 = buffer.data(fl + 27);
    const auto *fl_28 = buffer.data(fl + 28);
    const auto *fl_29 = buffer.data(fl + 29);
    const auto *fl_30 = buffer.data(fl + 30);
    const auto *fl_31 = buffer.data(fl + 31);
    const auto *fl_32 = buffer.data(fl + 32);
    const auto *fl_33 = buffer.data(fl + 33);
    const auto *fl_34 = buffer.data(fl + 34);
    const auto *fl_35 = buffer.data(fl + 35);
    const auto *fl_36 = buffer.data(fl + 36);
    const auto *fl_37 = buffer.data(fl + 37);
    const auto *fl_38 = buffer.data(fl + 38);
    const auto *fl_39 = buffer.data(fl + 39);
    const auto *fl_40 = buffer.data(fl + 40);
    const auto *fl_41 = buffer.data(fl + 41);
    const auto *fl_42 = buffer.data(fl + 42);
    const auto *fl_43 = buffer.data(fl + 43);
    const auto *fl_44 = buffer.data(fl + 44);
    const auto *fl_45 = buffer.data(fl + 45);
    const auto *fl_46 = buffer.data(fl + 46);
    const auto *fl_47 = buffer.data(fl + 47);
    const auto *fl_48 = buffer.data(fl + 48);
    const auto *fl_49 = buffer.data(fl + 49);
    const auto *fl_50 = buffer.data(fl + 50);
    const auto *fl_51 = buffer.data(fl + 51);
    const auto *fl_52 = buffer.data(fl + 52);
    const auto *fl_53 = buffer.data(fl + 53);
    const auto *fl_54 = buffer.data(fl + 54);
    const auto *fl_55 = buffer.data(fl + 55);
    const auto *fl_56 = buffer.data(fl + 56);
    const auto *fl_57 = buffer.data(fl + 57);
    const auto *fl_58 = buffer.data(fl + 58);
    const auto *fl_59 = buffer.data(fl + 59);

    const auto *hl_90 = buffer.data(hl + 90);
    const auto *hl_91 = buffer.data(hl + 91);
    const auto *hl_92 = buffer.data(hl + 92);
    const auto *hl_93 = buffer.data(hl + 93);
    const auto *hl_94 = buffer.data(hl + 94);
    const auto *hl_95 = buffer.data(hl + 95);
    const auto *hl_96 = buffer.data(hl + 96);
    const auto *hl_97 = buffer.data(hl + 97);
    const auto *hl_98 = buffer.data(hl + 98);
    const auto *hl_99 = buffer.data(hl + 99);
    const auto *hl_100 = buffer.data(hl + 100);
    const auto *hl_101 = buffer.data(hl + 101);
    const auto *hl_102 = buffer.data(hl + 102);
    const auto *hl_103 = buffer.data(hl + 103);
    const auto *hl_104 = buffer.data(hl + 104);
    const auto *hl_105 = buffer.data(hl + 105);
    const auto *hl_106 = buffer.data(hl + 106);
    const auto *hl_107 = buffer.data(hl + 107);
    const auto *hl_108 = buffer.data(hl + 108);
    const auto *hl_109 = buffer.data(hl + 109);
    const auto *hl_110 = buffer.data(hl + 110);
    const auto *hl_111 = buffer.data(hl + 111);
    const auto *hl_112 = buffer.data(hl + 112);
    const auto *hl_113 = buffer.data(hl + 113);
    const auto *hl_114 = buffer.data(hl + 114);
    const auto *hl_115 = buffer.data(hl + 115);
    const auto *hl_116 = buffer.data(hl + 116);
    const auto *hl_117 = buffer.data(hl + 117);
    const auto *hl_118 = buffer.data(hl + 118);
    const auto *hl_119 = buffer.data(hl + 119);
    const auto *hl_120 = buffer.data(hl + 120);
    const auto *hl_121 = buffer.data(hl + 121);
    const auto *hl_122 = buffer.data(hl + 122);
    const auto *hl_123 = buffer.data(hl + 123);
    const auto *hl_124 = buffer.data(hl + 124);
    const auto *hl_125 = buffer.data(hl + 125);
    const auto *hl_126 = buffer.data(hl + 126);
    const auto *hl_127 = buffer.data(hl + 127);
    const auto *hl_128 = buffer.data(hl + 128);
    const auto *hl_129 = buffer.data(hl + 129);
    const auto *hl_130 = buffer.data(hl + 130);
    const auto *hl_131 = buffer.data(hl + 131);
    const auto *hl_132 = buffer.data(hl + 132);
    const auto *hl_133 = buffer.data(hl + 133);
    const auto *hl_134 = buffer.data(hl + 134);
    const auto *hl_180 = buffer.data(hl + 180);
    const auto *hl_181 = buffer.data(hl + 181);
    const auto *hl_182 = buffer.data(hl + 182);
    const auto *hl_183 = buffer.data(hl + 183);
    const auto *hl_184 = buffer.data(hl + 184);
    const auto *hl_185 = buffer.data(hl + 185);
    const auto *hl_186 = buffer.data(hl + 186);
    const auto *hl_187 = buffer.data(hl + 187);
    const auto *hl_188 = buffer.data(hl + 188);
    const auto *hl_189 = buffer.data(hl + 189);
    const auto *hl_190 = buffer.data(hl + 190);
    const auto *hl_191 = buffer.data(hl + 191);
    const auto *hl_192 = buffer.data(hl + 192);
    const auto *hl_193 = buffer.data(hl + 193);
    const auto *hl_194 = buffer.data(hl + 194);
    const auto *hl_195 = buffer.data(hl + 195);
    const auto *hl_196 = buffer.data(hl + 196);
    const auto *hl_197 = buffer.data(hl + 197);
    const auto *hl_198 = buffer.data(hl + 198);
    const auto *hl_199 = buffer.data(hl + 199);
    const auto *hl_200 = buffer.data(hl + 200);
    const auto *hl_201 = buffer.data(hl + 201);
    const auto *hl_202 = buffer.data(hl + 202);
    const auto *hl_203 = buffer.data(hl + 203);
    const auto *hl_204 = buffer.data(hl + 204);
    const auto *hl_205 = buffer.data(hl + 205);
    const auto *hl_206 = buffer.data(hl + 206);
    const auto *hl_207 = buffer.data(hl + 207);
    const auto *hl_208 = buffer.data(hl + 208);
    const auto *hl_209 = buffer.data(hl + 209);
    const auto *hl_210 = buffer.data(hl + 210);
    const auto *hl_211 = buffer.data(hl + 211);
    const auto *hl_212 = buffer.data(hl + 212);
    const auto *hl_213 = buffer.data(hl + 213);
    const auto *hl_214 = buffer.data(hl + 214);
    const auto *hl_215 = buffer.data(hl + 215);
    const auto *hl_216 = buffer.data(hl + 216);
    const auto *hl_217 = buffer.data(hl + 217);
    const auto *hl_218 = buffer.data(hl + 218);
    const auto *hl_219 = buffer.data(hl + 219);
    const auto *hl_220 = buffer.data(hl + 220);
    const auto *hl_221 = buffer.data(hl + 221);
    const auto *hl_222 = buffer.data(hl + 222);
    const auto *hl_223 = buffer.data(hl + 223);
    const auto *hl_224 = buffer.data(hl + 224);
    const auto *hl_225 = buffer.data(hl + 225);
    const auto *hl_226 = buffer.data(hl + 226);
    const auto *hl_227 = buffer.data(hl + 227);
    const auto *hl_228 = buffer.data(hl + 228);
    const auto *hl_229 = buffer.data(hl + 229);
    const auto *hl_230 = buffer.data(hl + 230);
    const auto *hl_231 = buffer.data(hl + 231);
    const auto *hl_232 = buffer.data(hl + 232);
    const auto *hl_233 = buffer.data(hl + 233);
    const auto *hl_234 = buffer.data(hl + 234);
    const auto *hl_235 = buffer.data(hl + 235);
    const auto *hl_236 = buffer.data(hl + 236);
    const auto *hl_237 = buffer.data(hl + 237);
    const auto *hl_238 = buffer.data(hl + 238);
    const auto *hl_239 = buffer.data(hl + 239);
    const auto *hl_240 = buffer.data(hl + 240);
    const auto *hl_241 = buffer.data(hl + 241);
    const auto *hl_242 = buffer.data(hl + 242);
    const auto *hl_243 = buffer.data(hl + 243);
    const auto *hl_244 = buffer.data(hl + 244);
    const auto *hl_245 = buffer.data(hl + 245);
    const auto *hl_246 = buffer.data(hl + 246);
    const auto *hl_247 = buffer.data(hl + 247);
    const auto *hl_248 = buffer.data(hl + 248);
    const auto *hl_249 = buffer.data(hl + 249);
    const auto *hl_250 = buffer.data(hl + 250);
    const auto *hl_251 = buffer.data(hl + 251);
    const auto *hl_252 = buffer.data(hl + 252);
    const auto *hl_253 = buffer.data(hl + 253);
    const auto *hl_254 = buffer.data(hl + 254);
    const auto *hl_255 = buffer.data(hl + 255);
    const auto *hl_256 = buffer.data(hl + 256);
    const auto *hl_257 = buffer.data(hl + 257);
    const auto *hl_258 = buffer.data(hl + 258);
    const auto *hl_259 = buffer.data(hl + 259);
    const auto *hl_260 = buffer.data(hl + 260);
    const auto *hl_261 = buffer.data(hl + 261);
    const auto *hl_262 = buffer.data(hl + 262);
    const auto *hl_263 = buffer.data(hl + 263);
    const auto *hl_264 = buffer.data(hl + 264);
    const auto *hl_265 = buffer.data(hl + 265);
    const auto *hl_266 = buffer.data(hl + 266);
    const auto *hl_267 = buffer.data(hl + 267);
    const auto *hl_268 = buffer.data(hl + 268);
    const auto *hl_269 = buffer.data(hl + 269);
    const auto *hl_315 = buffer.data(hl + 315);
    const auto *hl_316 = buffer.data(hl + 316);
    const auto *hl_317 = buffer.data(hl + 317);
    const auto *hl_318 = buffer.data(hl + 318);
    const auto *hl_319 = buffer.data(hl + 319);
    const auto *hl_320 = buffer.data(hl + 320);
    const auto *hl_321 = buffer.data(hl + 321);
    const auto *hl_322 = buffer.data(hl + 322);
    const auto *hl_323 = buffer.data(hl + 323);
    const auto *hl_324 = buffer.data(hl + 324);
    const auto *hl_325 = buffer.data(hl + 325);
    const auto *hl_326 = buffer.data(hl + 326);
    const auto *hl_327 = buffer.data(hl + 327);
    const auto *hl_328 = buffer.data(hl + 328);
    const auto *hl_329 = buffer.data(hl + 329);
    const auto *hl_330 = buffer.data(hl + 330);
    const auto *hl_331 = buffer.data(hl + 331);
    const auto *hl_332 = buffer.data(hl + 332);
    const auto *hl_333 = buffer.data(hl + 333);
    const auto *hl_334 = buffer.data(hl + 334);
    const auto *hl_335 = buffer.data(hl + 335);
    const auto *hl_336 = buffer.data(hl + 336);
    const auto *hl_337 = buffer.data(hl + 337);
    const auto *hl_338 = buffer.data(hl + 338);
    const auto *hl_339 = buffer.data(hl + 339);
    const auto *hl_340 = buffer.data(hl + 340);
    const auto *hl_341 = buffer.data(hl + 341);
    const auto *hl_342 = buffer.data(hl + 342);
    const auto *hl_343 = buffer.data(hl + 343);
    const auto *hl_344 = buffer.data(hl + 344);
    const auto *hl_345 = buffer.data(hl + 345);
    const auto *hl_346 = buffer.data(hl + 346);
    const auto *hl_347 = buffer.data(hl + 347);
    const auto *hl_348 = buffer.data(hl + 348);
    const auto *hl_349 = buffer.data(hl + 349);
    const auto *hl_350 = buffer.data(hl + 350);
    const auto *hl_351 = buffer.data(hl + 351);
    const auto *hl_352 = buffer.data(hl + 352);
    const auto *hl_353 = buffer.data(hl + 353);
    const auto *hl_354 = buffer.data(hl + 354);
    const auto *hl_355 = buffer.data(hl + 355);
    const auto *hl_356 = buffer.data(hl + 356);
    const auto *hl_357 = buffer.data(hl + 357);
    const auto *hl_358 = buffer.data(hl + 358);
    const auto *hl_359 = buffer.data(hl + 359);
    const auto *hl_360 = buffer.data(hl + 360);
    const auto *hl_361 = buffer.data(hl + 361);
    const auto *hl_362 = buffer.data(hl + 362);
    const auto *hl_363 = buffer.data(hl + 363);
    const auto *hl_364 = buffer.data(hl + 364);
    const auto *hl_365 = buffer.data(hl + 365);
    const auto *hl_366 = buffer.data(hl + 366);
    const auto *hl_367 = buffer.data(hl + 367);
    const auto *hl_368 = buffer.data(hl + 368);
    const auto *hl_369 = buffer.data(hl + 369);
    const auto *hl_370 = buffer.data(hl + 370);
    const auto *hl_371 = buffer.data(hl + 371);
    const auto *hl_372 = buffer.data(hl + 372);
    const auto *hl_373 = buffer.data(hl + 373);
    const auto *hl_374 = buffer.data(hl + 374);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hl_90, hl_91, hl_92, hl_93, \
                         hl_94, hl_95, hl_96, hl_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hl_90[k];

        t_1[k] = f_0 * hl_91[k];

        t_2[k] = f_0 * hl_92[k];

        t_3[k] = f_0 * hl_93[k];

        t_4[k] = f_0 * hl_94[k];

        t_5[k] = f_0 * hl_95[k];

        t_6[k] = f_0 * hl_96[k];

        t_7[k] = f_0 * hl_97[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, hl_98, hl_99, hl_100, \
                         hl_101, hl_102, hl_103, hl_104, hl_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hl_98[k];

        t_9[k] = f_0 * hl_99[k];

        t_10[k] = f_0 * hl_100[k];

        t_11[k] = f_0 * hl_101[k];

        t_12[k] = f_0 * hl_102[k];

        t_13[k] = f_0 * hl_103[k];

        t_14[k] = f_0 * hl_104[k];

        t_15[k] = f_0 * hl_105[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, hl_106, hl_107, \
                         hl_108, hl_109, hl_110, hl_111, hl_112, \
                         hl_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * hl_106[k];

        t_17[k] = f_0 * hl_107[k];

        t_18[k] = f_0 * hl_108[k];

        t_19[k] = f_0 * hl_109[k];

        t_20[k] = f_0 * hl_110[k];

        t_21[k] = f_0 * hl_111[k];

        t_22[k] = f_0 * hl_112[k];

        t_23[k] = f_0 * hl_113[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, hl_114, hl_115, \
                         hl_116, hl_117, hl_118, hl_119, hl_120, \
                         hl_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * hl_114[k];

        t_25[k] = f_0 * hl_115[k];

        t_26[k] = f_0 * hl_116[k];

        t_27[k] = f_0 * hl_117[k];

        t_28[k] = f_0 * hl_118[k];

        t_29[k] = f_0 * hl_119[k];

        t_30[k] = f_0 * hl_120[k];

        t_31[k] = f_0 * hl_121[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, hl_122, hl_123, \
                         hl_124, hl_125, hl_126, hl_127, hl_128, \
                         hl_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * hl_122[k];

        t_33[k] = f_0 * hl_123[k];

        t_34[k] = f_0 * hl_124[k];

        t_35[k] = f_0 * hl_125[k];

        t_36[k] = f_0 * hl_126[k];

        t_37[k] = f_0 * hl_127[k];

        t_38[k] = f_0 * hl_128[k];

        t_39[k] = f_0 * hl_129[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, hl_130, hl_131, \
                         hl_132, hl_133, hl_134, hl_180, hl_181, \
                         hl_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * hl_130[k];

        t_41[k] = f_0 * hl_131[k];

        t_42[k] = f_0 * hl_132[k];

        t_43[k] = f_0 * hl_133[k];

        t_44[k] = f_0 * hl_134[k];

        t_45[k] = f_0 * hl_180[k];

        t_46[k] = f_0 * hl_181[k];

        t_47[k] = f_0 * hl_182[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, hl_183, hl_184, \
                         hl_185, hl_186, hl_187, hl_188, hl_189, \
                         hl_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * hl_183[k];

        t_49[k] = f_0 * hl_184[k];

        t_50[k] = f_0 * hl_185[k];

        t_51[k] = f_0 * hl_186[k];

        t_52[k] = f_0 * hl_187[k];

        t_53[k] = f_0 * hl_188[k];

        t_54[k] = f_0 * hl_189[k];

        t_55[k] = f_0 * hl_190[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, hl_191, hl_192, \
                         hl_193, hl_194, hl_195, hl_196, hl_197, \
                         hl_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * hl_191[k];

        t_57[k] = f_0 * hl_192[k];

        t_58[k] = f_0 * hl_193[k];

        t_59[k] = f_0 * hl_194[k];

        t_60[k] = f_0 * hl_195[k];

        t_61[k] = f_0 * hl_196[k];

        t_62[k] = f_0 * hl_197[k];

        t_63[k] = f_0 * hl_198[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, hl_199, hl_200, \
                         hl_201, hl_202, hl_203, hl_204, hl_205, \
                         hl_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * hl_199[k];

        t_65[k] = f_0 * hl_200[k];

        t_66[k] = f_0 * hl_201[k];

        t_67[k] = f_0 * hl_202[k];

        t_68[k] = f_0 * hl_203[k];

        t_69[k] = f_0 * hl_204[k];

        t_70[k] = f_0 * hl_205[k];

        t_71[k] = f_0 * hl_206[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, t_78, t_79, hl_207, hl_208, \
                         hl_209, hl_210, hl_211, hl_212, hl_213, \
                         hl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * hl_207[k];

        t_73[k] = f_0 * hl_208[k];

        t_74[k] = f_0 * hl_209[k];

        t_75[k] = f_0 * hl_210[k];

        t_76[k] = f_0 * hl_211[k];

        t_77[k] = f_0 * hl_212[k];

        t_78[k] = f_0 * hl_213[k];

        t_79[k] = f_0 * hl_214[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, t_86, t_87, hl_215, hl_216, \
                         hl_217, hl_218, hl_219, hl_220, hl_221, \
                         hl_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * hl_215[k];

        t_81[k] = f_0 * hl_216[k];

        t_82[k] = f_0 * hl_217[k];

        t_83[k] = f_0 * hl_218[k];

        t_84[k] = f_0 * hl_219[k];

        t_85[k] = f_0 * hl_220[k];

        t_86[k] = f_0 * hl_221[k];

        t_87[k] = f_0 * hl_222[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, t_93, fl_0, fl_1, fl_2, fl_3, hl_223, \
                         hl_224, hl_225, hl_226, hl_227, hl_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_0 * hl_223[k];

        t_89[k] = f_0 * hl_224[k];

        t_90[k] = -fl_0[k]
                  + f_0 * hl_225[k];

        t_91[k] = -fl_1[k]
                  + f_0 * hl_226[k];

        t_92[k] = -fl_2[k]
                  + f_0 * hl_227[k];

        t_93[k] = -fl_3[k]
                  + f_0 * hl_228[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, fl_4, fl_5, fl_6, fl_7, fl_8, hl_229, \
                         hl_230, hl_231, hl_232, hl_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -fl_4[k]
                  + f_0 * hl_229[k];

        t_95[k] = -fl_5[k]
                  + f_0 * hl_230[k];

        t_96[k] = -fl_6[k]
                  + f_0 * hl_231[k];

        t_97[k] = -fl_7[k]
                  + f_0 * hl_232[k];

        t_98[k] = -fl_8[k]
                  + f_0 * hl_233[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, fl_9, fl_10, fl_11, fl_12, fl_13, \
                         hl_234, hl_235, hl_236, hl_237, hl_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -fl_9[k]
                  + f_0 * hl_234[k];

        t_100[k] = -fl_10[k]
                   + f_0 * hl_235[k];

        t_101[k] = -fl_11[k]
                   + f_0 * hl_236[k];

        t_102[k] = -fl_12[k]
                   + f_0 * hl_237[k];

        t_103[k] = -fl_13[k]
                   + f_0 * hl_238[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, fl_14, fl_15, fl_16, fl_17, fl_18, \
                         hl_239, hl_240, hl_241, hl_242, hl_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -fl_14[k]
                   + f_0 * hl_239[k];

        t_105[k] = -fl_15[k]
                   + f_0 * hl_240[k];

        t_106[k] = -fl_16[k]
                   + f_0 * hl_241[k];

        t_107[k] = -fl_17[k]
                   + f_0 * hl_242[k];

        t_108[k] = -fl_18[k]
                   + f_0 * hl_243[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, fl_19, fl_20, fl_21, fl_22, fl_23, \
                         hl_244, hl_245, hl_246, hl_247, hl_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -fl_19[k]
                   + f_0 * hl_244[k];

        t_110[k] = -fl_20[k]
                   + f_0 * hl_245[k];

        t_111[k] = -fl_21[k]
                   + f_0 * hl_246[k];

        t_112[k] = -fl_22[k]
                   + f_0 * hl_247[k];

        t_113[k] = -fl_23[k]
                   + f_0 * hl_248[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, fl_24, fl_25, fl_26, fl_27, fl_28, \
                         hl_249, hl_250, hl_251, hl_252, hl_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -fl_24[k]
                   + f_0 * hl_249[k];

        t_115[k] = -fl_25[k]
                   + f_0 * hl_250[k];

        t_116[k] = -fl_26[k]
                   + f_0 * hl_251[k];

        t_117[k] = -fl_27[k]
                   + f_0 * hl_252[k];

        t_118[k] = -fl_28[k]
                   + f_0 * hl_253[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, fl_29, fl_30, fl_31, fl_32, fl_33, \
                         hl_254, hl_255, hl_256, hl_257, hl_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -fl_29[k]
                   + f_0 * hl_254[k];

        t_120[k] = -fl_30[k]
                   + f_0 * hl_255[k];

        t_121[k] = -fl_31[k]
                   + f_0 * hl_256[k];

        t_122[k] = -fl_32[k]
                   + f_0 * hl_257[k];

        t_123[k] = -fl_33[k]
                   + f_0 * hl_258[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, fl_34, fl_35, fl_36, fl_37, fl_38, \
                         hl_259, hl_260, hl_261, hl_262, hl_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -fl_34[k]
                   + f_0 * hl_259[k];

        t_125[k] = -fl_35[k]
                   + f_0 * hl_260[k];

        t_126[k] = -fl_36[k]
                   + f_0 * hl_261[k];

        t_127[k] = -fl_37[k]
                   + f_0 * hl_262[k];

        t_128[k] = -fl_38[k]
                   + f_0 * hl_263[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, fl_39, fl_40, fl_41, fl_42, fl_43, \
                         hl_264, hl_265, hl_266, hl_267, hl_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -fl_39[k]
                   + f_0 * hl_264[k];

        t_130[k] = -fl_40[k]
                   + f_0 * hl_265[k];

        t_131[k] = -fl_41[k]
                   + f_0 * hl_266[k];

        t_132[k] = -fl_42[k]
                   + f_0 * hl_267[k];

        t_133[k] = -fl_43[k]
                   + f_0 * hl_268[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, t_139, t_140, fl_44, hl_269, \
                         hl_315, hl_316, hl_317, hl_318, hl_319, \
                         hl_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -fl_44[k]
                   + f_0 * hl_269[k];

        t_135[k] = f_0 * hl_315[k];

        t_136[k] = f_0 * hl_316[k];

        t_137[k] = f_0 * hl_317[k];

        t_138[k] = f_0 * hl_318[k];

        t_139[k] = f_0 * hl_319[k];

        t_140[k] = f_0 * hl_320[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, t_146, t_147, t_148, hl_321, \
                         hl_322, hl_323, hl_324, hl_325, hl_326, hl_327, \
                         hl_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_0 * hl_321[k];

        t_142[k] = f_0 * hl_322[k];

        t_143[k] = f_0 * hl_323[k];

        t_144[k] = f_0 * hl_324[k];

        t_145[k] = f_0 * hl_325[k];

        t_146[k] = f_0 * hl_326[k];

        t_147[k] = f_0 * hl_327[k];

        t_148[k] = f_0 * hl_328[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, t_154, t_155, t_156, hl_329, \
                         hl_330, hl_331, hl_332, hl_333, hl_334, hl_335, \
                         hl_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * hl_329[k];

        t_150[k] = f_0 * hl_330[k];

        t_151[k] = f_0 * hl_331[k];

        t_152[k] = f_0 * hl_332[k];

        t_153[k] = f_0 * hl_333[k];

        t_154[k] = f_0 * hl_334[k];

        t_155[k] = f_0 * hl_335[k];

        t_156[k] = f_0 * hl_336[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, t_162, t_163, t_164, hl_337, \
                         hl_338, hl_339, hl_340, hl_341, hl_342, hl_343, \
                         hl_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_0 * hl_337[k];

        t_158[k] = f_0 * hl_338[k];

        t_159[k] = f_0 * hl_339[k];

        t_160[k] = f_0 * hl_340[k];

        t_161[k] = f_0 * hl_341[k];

        t_162[k] = f_0 * hl_342[k];

        t_163[k] = f_0 * hl_343[k];

        t_164[k] = f_0 * hl_344[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, t_170, t_171, t_172, hl_345, \
                         hl_346, hl_347, hl_348, hl_349, hl_350, hl_351, \
                         hl_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_0 * hl_345[k];

        t_166[k] = f_0 * hl_346[k];

        t_167[k] = f_0 * hl_347[k];

        t_168[k] = f_0 * hl_348[k];

        t_169[k] = f_0 * hl_349[k];

        t_170[k] = f_0 * hl_350[k];

        t_171[k] = f_0 * hl_351[k];

        t_172[k] = f_0 * hl_352[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, t_179, hl_353, hl_354, \
                         hl_355, hl_356, hl_357, hl_358, hl_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_0 * hl_353[k];

        t_174[k] = f_0 * hl_354[k];

        t_175[k] = f_0 * hl_355[k];

        t_176[k] = f_0 * hl_356[k];

        t_177[k] = f_0 * hl_357[k];

        t_178[k] = f_0 * hl_358[k];

        t_179[k] = f_0 * hl_359[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, fl_45, fl_46, fl_47, fl_48, fl_49, \
                         hl_360, hl_361, hl_362, hl_363, hl_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -fl_45[k]
                   + f_0 * hl_360[k];

        t_181[k] = -fl_46[k]
                   + f_0 * hl_361[k];

        t_182[k] = -fl_47[k]
                   + f_0 * hl_362[k];

        t_183[k] = -fl_48[k]
                   + f_0 * hl_363[k];

        t_184[k] = -fl_49[k]
                   + f_0 * hl_364[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, fl_50, fl_51, fl_52, fl_53, fl_54, \
                         hl_365, hl_366, hl_367, hl_368, hl_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -fl_50[k]
                   + f_0 * hl_365[k];

        t_186[k] = -fl_51[k]
                   + f_0 * hl_366[k];

        t_187[k] = -fl_52[k]
                   + f_0 * hl_367[k];

        t_188[k] = -fl_53[k]
                   + f_0 * hl_368[k];

        t_189[k] = -fl_54[k]
                   + f_0 * hl_369[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, fl_55, fl_56, fl_57, fl_58, fl_59, \
                         hl_370, hl_371, hl_372, hl_373, hl_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -fl_55[k]
                   + f_0 * hl_370[k];

        t_191[k] = -fl_56[k]
                   + f_0 * hl_371[k];

        t_192[k] = -fl_57[k]
                   + f_0 * hl_372[k];

        t_193[k] = -fl_58[k]
                   + f_0 * hl_373[k];

        t_194[k] = -fl_59[k]
                   + f_0 * hl_374[k];
    }
}

static auto
compute_prim_geom_10_gl_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);

    const auto *fl_60 = buffer.data(fl + 60);
    const auto *fl_61 = buffer.data(fl + 61);
    const auto *fl_62 = buffer.data(fl + 62);
    const auto *fl_63 = buffer.data(fl + 63);
    const auto *fl_64 = buffer.data(fl + 64);
    const auto *fl_65 = buffer.data(fl + 65);
    const auto *fl_66 = buffer.data(fl + 66);
    const auto *fl_67 = buffer.data(fl + 67);
    const auto *fl_68 = buffer.data(fl + 68);
    const auto *fl_69 = buffer.data(fl + 69);
    const auto *fl_70 = buffer.data(fl + 70);
    const auto *fl_71 = buffer.data(fl + 71);
    const auto *fl_72 = buffer.data(fl + 72);
    const auto *fl_73 = buffer.data(fl + 73);
    const auto *fl_74 = buffer.data(fl + 74);
    const auto *fl_75 = buffer.data(fl + 75);
    const auto *fl_76 = buffer.data(fl + 76);
    const auto *fl_77 = buffer.data(fl + 77);
    const auto *fl_78 = buffer.data(fl + 78);
    const auto *fl_79 = buffer.data(fl + 79);
    const auto *fl_80 = buffer.data(fl + 80);
    const auto *fl_81 = buffer.data(fl + 81);
    const auto *fl_82 = buffer.data(fl + 82);
    const auto *fl_83 = buffer.data(fl + 83);
    const auto *fl_84 = buffer.data(fl + 84);
    const auto *fl_85 = buffer.data(fl + 85);
    const auto *fl_86 = buffer.data(fl + 86);
    const auto *fl_87 = buffer.data(fl + 87);
    const auto *fl_88 = buffer.data(fl + 88);
    const auto *fl_89 = buffer.data(fl + 89);
    const auto *fl_90 = buffer.data(fl + 90);
    const auto *fl_91 = buffer.data(fl + 91);
    const auto *fl_92 = buffer.data(fl + 92);
    const auto *fl_93 = buffer.data(fl + 93);
    const auto *fl_94 = buffer.data(fl + 94);
    const auto *fl_95 = buffer.data(fl + 95);
    const auto *fl_96 = buffer.data(fl + 96);
    const auto *fl_97 = buffer.data(fl + 97);
    const auto *fl_98 = buffer.data(fl + 98);
    const auto *fl_99 = buffer.data(fl + 99);
    const auto *fl_100 = buffer.data(fl + 100);
    const auto *fl_101 = buffer.data(fl + 101);
    const auto *fl_102 = buffer.data(fl + 102);
    const auto *fl_103 = buffer.data(fl + 103);
    const auto *fl_104 = buffer.data(fl + 104);
    const auto *fl_105 = buffer.data(fl + 105);
    const auto *fl_106 = buffer.data(fl + 106);
    const auto *fl_107 = buffer.data(fl + 107);
    const auto *fl_108 = buffer.data(fl + 108);
    const auto *fl_109 = buffer.data(fl + 109);
    const auto *fl_110 = buffer.data(fl + 110);
    const auto *fl_111 = buffer.data(fl + 111);
    const auto *fl_112 = buffer.data(fl + 112);
    const auto *fl_113 = buffer.data(fl + 113);
    const auto *fl_114 = buffer.data(fl + 114);
    const auto *fl_115 = buffer.data(fl + 115);
    const auto *fl_116 = buffer.data(fl + 116);
    const auto *fl_117 = buffer.data(fl + 117);
    const auto *fl_118 = buffer.data(fl + 118);
    const auto *fl_119 = buffer.data(fl + 119);
    const auto *fl_120 = buffer.data(fl + 120);
    const auto *fl_121 = buffer.data(fl + 121);
    const auto *fl_122 = buffer.data(fl + 122);
    const auto *fl_123 = buffer.data(fl + 123);
    const auto *fl_124 = buffer.data(fl + 124);
    const auto *fl_125 = buffer.data(fl + 125);
    const auto *fl_126 = buffer.data(fl + 126);
    const auto *fl_127 = buffer.data(fl + 127);
    const auto *fl_128 = buffer.data(fl + 128);
    const auto *fl_129 = buffer.data(fl + 129);
    const auto *fl_130 = buffer.data(fl + 130);
    const auto *fl_131 = buffer.data(fl + 131);
    const auto *fl_132 = buffer.data(fl + 132);
    const auto *fl_133 = buffer.data(fl + 133);
    const auto *fl_134 = buffer.data(fl + 134);
    const auto *fl_135 = buffer.data(fl + 135);
    const auto *fl_136 = buffer.data(fl + 136);
    const auto *fl_137 = buffer.data(fl + 137);
    const auto *fl_138 = buffer.data(fl + 138);
    const auto *fl_139 = buffer.data(fl + 139);
    const auto *fl_140 = buffer.data(fl + 140);
    const auto *fl_141 = buffer.data(fl + 141);
    const auto *fl_142 = buffer.data(fl + 142);
    const auto *fl_143 = buffer.data(fl + 143);
    const auto *fl_144 = buffer.data(fl + 144);
    const auto *fl_145 = buffer.data(fl + 145);
    const auto *fl_146 = buffer.data(fl + 146);
    const auto *fl_147 = buffer.data(fl + 147);
    const auto *fl_148 = buffer.data(fl + 148);
    const auto *fl_149 = buffer.data(fl + 149);
    const auto *fl_150 = buffer.data(fl + 150);
    const auto *fl_151 = buffer.data(fl + 151);
    const auto *fl_152 = buffer.data(fl + 152);
    const auto *fl_153 = buffer.data(fl + 153);
    const auto *fl_154 = buffer.data(fl + 154);
    const auto *fl_155 = buffer.data(fl + 155);
    const auto *fl_156 = buffer.data(fl + 156);
    const auto *fl_157 = buffer.data(fl + 157);
    const auto *fl_158 = buffer.data(fl + 158);
    const auto *fl_159 = buffer.data(fl + 159);
    const auto *fl_160 = buffer.data(fl + 160);
    const auto *fl_161 = buffer.data(fl + 161);
    const auto *fl_162 = buffer.data(fl + 162);
    const auto *fl_163 = buffer.data(fl + 163);
    const auto *fl_164 = buffer.data(fl + 164);
    const auto *fl_165 = buffer.data(fl + 165);
    const auto *fl_166 = buffer.data(fl + 166);
    const auto *fl_167 = buffer.data(fl + 167);
    const auto *fl_168 = buffer.data(fl + 168);
    const auto *fl_169 = buffer.data(fl + 169);
    const auto *fl_170 = buffer.data(fl + 170);
    const auto *fl_171 = buffer.data(fl + 171);
    const auto *fl_172 = buffer.data(fl + 172);
    const auto *fl_173 = buffer.data(fl + 173);
    const auto *fl_174 = buffer.data(fl + 174);
    const auto *fl_175 = buffer.data(fl + 175);
    const auto *fl_176 = buffer.data(fl + 176);

    const auto *hl_375 = buffer.data(hl + 375);
    const auto *hl_376 = buffer.data(hl + 376);
    const auto *hl_377 = buffer.data(hl + 377);
    const auto *hl_378 = buffer.data(hl + 378);
    const auto *hl_379 = buffer.data(hl + 379);
    const auto *hl_380 = buffer.data(hl + 380);
    const auto *hl_381 = buffer.data(hl + 381);
    const auto *hl_382 = buffer.data(hl + 382);
    const auto *hl_383 = buffer.data(hl + 383);
    const auto *hl_384 = buffer.data(hl + 384);
    const auto *hl_385 = buffer.data(hl + 385);
    const auto *hl_386 = buffer.data(hl + 386);
    const auto *hl_387 = buffer.data(hl + 387);
    const auto *hl_388 = buffer.data(hl + 388);
    const auto *hl_389 = buffer.data(hl + 389);
    const auto *hl_390 = buffer.data(hl + 390);
    const auto *hl_391 = buffer.data(hl + 391);
    const auto *hl_392 = buffer.data(hl + 392);
    const auto *hl_393 = buffer.data(hl + 393);
    const auto *hl_394 = buffer.data(hl + 394);
    const auto *hl_395 = buffer.data(hl + 395);
    const auto *hl_396 = buffer.data(hl + 396);
    const auto *hl_397 = buffer.data(hl + 397);
    const auto *hl_398 = buffer.data(hl + 398);
    const auto *hl_399 = buffer.data(hl + 399);
    const auto *hl_400 = buffer.data(hl + 400);
    const auto *hl_401 = buffer.data(hl + 401);
    const auto *hl_402 = buffer.data(hl + 402);
    const auto *hl_403 = buffer.data(hl + 403);
    const auto *hl_404 = buffer.data(hl + 404);
    const auto *hl_405 = buffer.data(hl + 405);
    const auto *hl_406 = buffer.data(hl + 406);
    const auto *hl_407 = buffer.data(hl + 407);
    const auto *hl_408 = buffer.data(hl + 408);
    const auto *hl_409 = buffer.data(hl + 409);
    const auto *hl_410 = buffer.data(hl + 410);
    const auto *hl_411 = buffer.data(hl + 411);
    const auto *hl_412 = buffer.data(hl + 412);
    const auto *hl_413 = buffer.data(hl + 413);
    const auto *hl_414 = buffer.data(hl + 414);
    const auto *hl_415 = buffer.data(hl + 415);
    const auto *hl_416 = buffer.data(hl + 416);
    const auto *hl_417 = buffer.data(hl + 417);
    const auto *hl_418 = buffer.data(hl + 418);
    const auto *hl_419 = buffer.data(hl + 419);
    const auto *hl_420 = buffer.data(hl + 420);
    const auto *hl_421 = buffer.data(hl + 421);
    const auto *hl_422 = buffer.data(hl + 422);
    const auto *hl_423 = buffer.data(hl + 423);
    const auto *hl_424 = buffer.data(hl + 424);
    const auto *hl_425 = buffer.data(hl + 425);
    const auto *hl_426 = buffer.data(hl + 426);
    const auto *hl_427 = buffer.data(hl + 427);
    const auto *hl_428 = buffer.data(hl + 428);
    const auto *hl_429 = buffer.data(hl + 429);
    const auto *hl_430 = buffer.data(hl + 430);
    const auto *hl_431 = buffer.data(hl + 431);
    const auto *hl_432 = buffer.data(hl + 432);
    const auto *hl_433 = buffer.data(hl + 433);
    const auto *hl_434 = buffer.data(hl + 434);
    const auto *hl_435 = buffer.data(hl + 435);
    const auto *hl_436 = buffer.data(hl + 436);
    const auto *hl_437 = buffer.data(hl + 437);
    const auto *hl_438 = buffer.data(hl + 438);
    const auto *hl_439 = buffer.data(hl + 439);
    const auto *hl_440 = buffer.data(hl + 440);
    const auto *hl_441 = buffer.data(hl + 441);
    const auto *hl_442 = buffer.data(hl + 442);
    const auto *hl_443 = buffer.data(hl + 443);
    const auto *hl_444 = buffer.data(hl + 444);
    const auto *hl_445 = buffer.data(hl + 445);
    const auto *hl_446 = buffer.data(hl + 446);
    const auto *hl_447 = buffer.data(hl + 447);
    const auto *hl_448 = buffer.data(hl + 448);
    const auto *hl_449 = buffer.data(hl + 449);
    const auto *hl_495 = buffer.data(hl + 495);
    const auto *hl_496 = buffer.data(hl + 496);
    const auto *hl_497 = buffer.data(hl + 497);
    const auto *hl_498 = buffer.data(hl + 498);
    const auto *hl_499 = buffer.data(hl + 499);
    const auto *hl_500 = buffer.data(hl + 500);
    const auto *hl_501 = buffer.data(hl + 501);
    const auto *hl_502 = buffer.data(hl + 502);
    const auto *hl_503 = buffer.data(hl + 503);
    const auto *hl_504 = buffer.data(hl + 504);
    const auto *hl_505 = buffer.data(hl + 505);
    const auto *hl_506 = buffer.data(hl + 506);
    const auto *hl_507 = buffer.data(hl + 507);
    const auto *hl_508 = buffer.data(hl + 508);
    const auto *hl_509 = buffer.data(hl + 509);
    const auto *hl_510 = buffer.data(hl + 510);
    const auto *hl_511 = buffer.data(hl + 511);
    const auto *hl_512 = buffer.data(hl + 512);
    const auto *hl_513 = buffer.data(hl + 513);
    const auto *hl_514 = buffer.data(hl + 514);
    const auto *hl_515 = buffer.data(hl + 515);
    const auto *hl_516 = buffer.data(hl + 516);
    const auto *hl_517 = buffer.data(hl + 517);
    const auto *hl_518 = buffer.data(hl + 518);
    const auto *hl_519 = buffer.data(hl + 519);
    const auto *hl_520 = buffer.data(hl + 520);
    const auto *hl_521 = buffer.data(hl + 521);
    const auto *hl_522 = buffer.data(hl + 522);
    const auto *hl_523 = buffer.data(hl + 523);
    const auto *hl_524 = buffer.data(hl + 524);
    const auto *hl_525 = buffer.data(hl + 525);
    const auto *hl_526 = buffer.data(hl + 526);
    const auto *hl_527 = buffer.data(hl + 527);
    const auto *hl_528 = buffer.data(hl + 528);
    const auto *hl_529 = buffer.data(hl + 529);
    const auto *hl_530 = buffer.data(hl + 530);
    const auto *hl_531 = buffer.data(hl + 531);
    const auto *hl_532 = buffer.data(hl + 532);
    const auto *hl_533 = buffer.data(hl + 533);
    const auto *hl_534 = buffer.data(hl + 534);
    const auto *hl_535 = buffer.data(hl + 535);
    const auto *hl_536 = buffer.data(hl + 536);
    const auto *hl_537 = buffer.data(hl + 537);
    const auto *hl_538 = buffer.data(hl + 538);
    const auto *hl_539 = buffer.data(hl + 539);
    const auto *hl_540 = buffer.data(hl + 540);
    const auto *hl_541 = buffer.data(hl + 541);
    const auto *hl_542 = buffer.data(hl + 542);
    const auto *hl_543 = buffer.data(hl + 543);
    const auto *hl_544 = buffer.data(hl + 544);
    const auto *hl_545 = buffer.data(hl + 545);
    const auto *hl_546 = buffer.data(hl + 546);
    const auto *hl_547 = buffer.data(hl + 547);
    const auto *hl_548 = buffer.data(hl + 548);
    const auto *hl_549 = buffer.data(hl + 549);
    const auto *hl_550 = buffer.data(hl + 550);
    const auto *hl_551 = buffer.data(hl + 551);
    const auto *hl_552 = buffer.data(hl + 552);
    const auto *hl_553 = buffer.data(hl + 553);
    const auto *hl_554 = buffer.data(hl + 554);
    const auto *hl_555 = buffer.data(hl + 555);
    const auto *hl_556 = buffer.data(hl + 556);
    const auto *hl_557 = buffer.data(hl + 557);
    const auto *hl_558 = buffer.data(hl + 558);
    const auto *hl_559 = buffer.data(hl + 559);
    const auto *hl_560 = buffer.data(hl + 560);
    const auto *hl_561 = buffer.data(hl + 561);
    const auto *hl_562 = buffer.data(hl + 562);
    const auto *hl_563 = buffer.data(hl + 563);
    const auto *hl_564 = buffer.data(hl + 564);
    const auto *hl_565 = buffer.data(hl + 565);
    const auto *hl_566 = buffer.data(hl + 566);
    const auto *hl_567 = buffer.data(hl + 567);
    const auto *hl_568 = buffer.data(hl + 568);
    const auto *hl_569 = buffer.data(hl + 569);
    const auto *hl_570 = buffer.data(hl + 570);
    const auto *hl_571 = buffer.data(hl + 571);
    const auto *hl_572 = buffer.data(hl + 572);
    const auto *hl_573 = buffer.data(hl + 573);
    const auto *hl_574 = buffer.data(hl + 574);
    const auto *hl_575 = buffer.data(hl + 575);
    const auto *hl_576 = buffer.data(hl + 576);
    const auto *hl_577 = buffer.data(hl + 577);
    const auto *hl_578 = buffer.data(hl + 578);
    const auto *hl_579 = buffer.data(hl + 579);
    const auto *hl_580 = buffer.data(hl + 580);
    const auto *hl_581 = buffer.data(hl + 581);

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, fl_60, fl_61, fl_62, fl_63, fl_64, \
                         hl_375, hl_376, hl_377, hl_378, hl_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -fl_60[k]
                   + f_0 * hl_375[k];

        t_196[k] = -fl_61[k]
                   + f_0 * hl_376[k];

        t_197[k] = -fl_62[k]
                   + f_0 * hl_377[k];

        t_198[k] = -fl_63[k]
                   + f_0 * hl_378[k];

        t_199[k] = -fl_64[k]
                   + f_0 * hl_379[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, fl_65, fl_66, fl_67, fl_68, fl_69, \
                         hl_380, hl_381, hl_382, hl_383, hl_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -fl_65[k]
                   + f_0 * hl_380[k];

        t_201[k] = -fl_66[k]
                   + f_0 * hl_381[k];

        t_202[k] = -fl_67[k]
                   + f_0 * hl_382[k];

        t_203[k] = -fl_68[k]
                   + f_0 * hl_383[k];

        t_204[k] = -fl_69[k]
                   + f_0 * hl_384[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, fl_70, fl_71, fl_72, fl_73, fl_74, \
                         hl_385, hl_386, hl_387, hl_388, hl_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -fl_70[k]
                   + f_0 * hl_385[k];

        t_206[k] = -fl_71[k]
                   + f_0 * hl_386[k];

        t_207[k] = -fl_72[k]
                   + f_0 * hl_387[k];

        t_208[k] = -fl_73[k]
                   + f_0 * hl_388[k];

        t_209[k] = -fl_74[k]
                   + f_0 * hl_389[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, fl_75, fl_76, fl_77, fl_78, fl_79, \
                         hl_390, hl_391, hl_392, hl_393, hl_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -fl_75[k]
                   + f_0 * hl_390[k];

        t_211[k] = -fl_76[k]
                   + f_0 * hl_391[k];

        t_212[k] = -fl_77[k]
                   + f_0 * hl_392[k];

        t_213[k] = -fl_78[k]
                   + f_0 * hl_393[k];

        t_214[k] = -fl_79[k]
                   + f_0 * hl_394[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, fl_80, fl_81, fl_82, fl_83, fl_84, \
                         hl_395, hl_396, hl_397, hl_398, hl_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -fl_80[k]
                   + f_0 * hl_395[k];

        t_216[k] = -fl_81[k]
                   + f_0 * hl_396[k];

        t_217[k] = -fl_82[k]
                   + f_0 * hl_397[k];

        t_218[k] = -fl_83[k]
                   + f_0 * hl_398[k];

        t_219[k] = -fl_84[k]
                   + f_0 * hl_399[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, fl_85, fl_86, fl_87, fl_88, fl_89, \
                         hl_400, hl_401, hl_402, hl_403, hl_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -fl_85[k]
                   + f_0 * hl_400[k];

        t_221[k] = -fl_86[k]
                   + f_0 * hl_401[k];

        t_222[k] = -fl_87[k]
                   + f_0 * hl_402[k];

        t_223[k] = -fl_88[k]
                   + f_0 * hl_403[k];

        t_224[k] = -fl_89[k]
                   + f_0 * hl_404[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, fl_90, fl_91, fl_92, fl_93, fl_94, \
                         hl_405, hl_406, hl_407, hl_408, hl_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -2.0 * fl_90[k]
                   + f_0 * hl_405[k];

        t_226[k] = -2.0 * fl_91[k]
                   + f_0 * hl_406[k];

        t_227[k] = -2.0 * fl_92[k]
                   + f_0 * hl_407[k];

        t_228[k] = -2.0 * fl_93[k]
                   + f_0 * hl_408[k];

        t_229[k] = -2.0 * fl_94[k]
                   + f_0 * hl_409[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, fl_95, fl_96, fl_97, fl_98, fl_99, \
                         hl_410, hl_411, hl_412, hl_413, hl_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -2.0 * fl_95[k]
                   + f_0 * hl_410[k];

        t_231[k] = -2.0 * fl_96[k]
                   + f_0 * hl_411[k];

        t_232[k] = -2.0 * fl_97[k]
                   + f_0 * hl_412[k];

        t_233[k] = -2.0 * fl_98[k]
                   + f_0 * hl_413[k];

        t_234[k] = -2.0 * fl_99[k]
                   + f_0 * hl_414[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, fl_100, fl_101, fl_102, fl_103, \
                         fl_104, hl_415, hl_416, hl_417, hl_418, \
                         hl_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -2.0 * fl_100[k]
                   + f_0 * hl_415[k];

        t_236[k] = -2.0 * fl_101[k]
                   + f_0 * hl_416[k];

        t_237[k] = -2.0 * fl_102[k]
                   + f_0 * hl_417[k];

        t_238[k] = -2.0 * fl_103[k]
                   + f_0 * hl_418[k];

        t_239[k] = -2.0 * fl_104[k]
                   + f_0 * hl_419[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, fl_105, fl_106, fl_107, fl_108, \
                         fl_109, hl_420, hl_421, hl_422, hl_423, \
                         hl_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -2.0 * fl_105[k]
                   + f_0 * hl_420[k];

        t_241[k] = -2.0 * fl_106[k]
                   + f_0 * hl_421[k];

        t_242[k] = -2.0 * fl_107[k]
                   + f_0 * hl_422[k];

        t_243[k] = -2.0 * fl_108[k]
                   + f_0 * hl_423[k];

        t_244[k] = -2.0 * fl_109[k]
                   + f_0 * hl_424[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, fl_110, fl_111, fl_112, fl_113, \
                         fl_114, hl_425, hl_426, hl_427, hl_428, \
                         hl_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -2.0 * fl_110[k]
                   + f_0 * hl_425[k];

        t_246[k] = -2.0 * fl_111[k]
                   + f_0 * hl_426[k];

        t_247[k] = -2.0 * fl_112[k]
                   + f_0 * hl_427[k];

        t_248[k] = -2.0 * fl_113[k]
                   + f_0 * hl_428[k];

        t_249[k] = -2.0 * fl_114[k]
                   + f_0 * hl_429[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, fl_115, fl_116, fl_117, fl_118, \
                         fl_119, hl_430, hl_431, hl_432, hl_433, \
                         hl_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * fl_115[k]
                   + f_0 * hl_430[k];

        t_251[k] = -2.0 * fl_116[k]
                   + f_0 * hl_431[k];

        t_252[k] = -2.0 * fl_117[k]
                   + f_0 * hl_432[k];

        t_253[k] = -2.0 * fl_118[k]
                   + f_0 * hl_433[k];

        t_254[k] = -2.0 * fl_119[k]
                   + f_0 * hl_434[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, fl_120, fl_121, fl_122, fl_123, \
                         fl_124, hl_435, hl_436, hl_437, hl_438, \
                         hl_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * fl_120[k]
                   + f_0 * hl_435[k];

        t_256[k] = -2.0 * fl_121[k]
                   + f_0 * hl_436[k];

        t_257[k] = -2.0 * fl_122[k]
                   + f_0 * hl_437[k];

        t_258[k] = -2.0 * fl_123[k]
                   + f_0 * hl_438[k];

        t_259[k] = -2.0 * fl_124[k]
                   + f_0 * hl_439[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, fl_125, fl_126, fl_127, fl_128, \
                         fl_129, hl_440, hl_441, hl_442, hl_443, \
                         hl_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * fl_125[k]
                   + f_0 * hl_440[k];

        t_261[k] = -2.0 * fl_126[k]
                   + f_0 * hl_441[k];

        t_262[k] = -2.0 * fl_127[k]
                   + f_0 * hl_442[k];

        t_263[k] = -2.0 * fl_128[k]
                   + f_0 * hl_443[k];

        t_264[k] = -2.0 * fl_129[k]
                   + f_0 * hl_444[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, fl_130, fl_131, fl_132, fl_133, \
                         fl_134, hl_445, hl_446, hl_447, hl_448, \
                         hl_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * fl_130[k]
                   + f_0 * hl_445[k];

        t_266[k] = -2.0 * fl_131[k]
                   + f_0 * hl_446[k];

        t_267[k] = -2.0 * fl_132[k]
                   + f_0 * hl_447[k];

        t_268[k] = -2.0 * fl_133[k]
                   + f_0 * hl_448[k];

        t_269[k] = -2.0 * fl_134[k]
                   + f_0 * hl_449[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, t_276, t_277, hl_495, \
                         hl_496, hl_497, hl_498, hl_499, hl_500, hl_501, \
                         hl_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_0 * hl_495[k];

        t_271[k] = f_0 * hl_496[k];

        t_272[k] = f_0 * hl_497[k];

        t_273[k] = f_0 * hl_498[k];

        t_274[k] = f_0 * hl_499[k];

        t_275[k] = f_0 * hl_500[k];

        t_276[k] = f_0 * hl_501[k];

        t_277[k] = f_0 * hl_502[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, t_283, t_284, t_285, hl_503, \
                         hl_504, hl_505, hl_506, hl_507, hl_508, hl_509, \
                         hl_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_0 * hl_503[k];

        t_279[k] = f_0 * hl_504[k];

        t_280[k] = f_0 * hl_505[k];

        t_281[k] = f_0 * hl_506[k];

        t_282[k] = f_0 * hl_507[k];

        t_283[k] = f_0 * hl_508[k];

        t_284[k] = f_0 * hl_509[k];

        t_285[k] = f_0 * hl_510[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, t_291, t_292, t_293, hl_511, \
                         hl_512, hl_513, hl_514, hl_515, hl_516, hl_517, \
                         hl_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * hl_511[k];

        t_287[k] = f_0 * hl_512[k];

        t_288[k] = f_0 * hl_513[k];

        t_289[k] = f_0 * hl_514[k];

        t_290[k] = f_0 * hl_515[k];

        t_291[k] = f_0 * hl_516[k];

        t_292[k] = f_0 * hl_517[k];

        t_293[k] = f_0 * hl_518[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, t_299, t_300, t_301, hl_519, \
                         hl_520, hl_521, hl_522, hl_523, hl_524, hl_525, \
                         hl_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_0 * hl_519[k];

        t_295[k] = f_0 * hl_520[k];

        t_296[k] = f_0 * hl_521[k];

        t_297[k] = f_0 * hl_522[k];

        t_298[k] = f_0 * hl_523[k];

        t_299[k] = f_0 * hl_524[k];

        t_300[k] = f_0 * hl_525[k];

        t_301[k] = f_0 * hl_526[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, t_307, t_308, t_309, hl_527, \
                         hl_528, hl_529, hl_530, hl_531, hl_532, hl_533, \
                         hl_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_0 * hl_527[k];

        t_303[k] = f_0 * hl_528[k];

        t_304[k] = f_0 * hl_529[k];

        t_305[k] = f_0 * hl_530[k];

        t_306[k] = f_0 * hl_531[k];

        t_307[k] = f_0 * hl_532[k];

        t_308[k] = f_0 * hl_533[k];

        t_309[k] = f_0 * hl_534[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, t_315, t_316, fl_135, fl_136, \
                         hl_535, hl_536, hl_537, hl_538, hl_539, hl_540, \
                         hl_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_0 * hl_535[k];

        t_311[k] = f_0 * hl_536[k];

        t_312[k] = f_0 * hl_537[k];

        t_313[k] = f_0 * hl_538[k];

        t_314[k] = f_0 * hl_539[k];

        t_315[k] = -fl_135[k]
                   + f_0 * hl_540[k];

        t_316[k] = -fl_136[k]
                   + f_0 * hl_541[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, fl_137, fl_138, fl_139, fl_140, \
                         fl_141, hl_542, hl_543, hl_544, hl_545, \
                         hl_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = -fl_137[k]
                   + f_0 * hl_542[k];

        t_318[k] = -fl_138[k]
                   + f_0 * hl_543[k];

        t_319[k] = -fl_139[k]
                   + f_0 * hl_544[k];

        t_320[k] = -fl_140[k]
                   + f_0 * hl_545[k];

        t_321[k] = -fl_141[k]
                   + f_0 * hl_546[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, fl_142, fl_143, fl_144, fl_145, \
                         fl_146, hl_547, hl_548, hl_549, hl_550, \
                         hl_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = -fl_142[k]
                   + f_0 * hl_547[k];

        t_323[k] = -fl_143[k]
                   + f_0 * hl_548[k];

        t_324[k] = -fl_144[k]
                   + f_0 * hl_549[k];

        t_325[k] = -fl_145[k]
                   + f_0 * hl_550[k];

        t_326[k] = -fl_146[k]
                   + f_0 * hl_551[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, fl_147, fl_148, fl_149, fl_150, \
                         fl_151, hl_552, hl_553, hl_554, hl_555, \
                         hl_556 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = -fl_147[k]
                   + f_0 * hl_552[k];

        t_328[k] = -fl_148[k]
                   + f_0 * hl_553[k];

        t_329[k] = -fl_149[k]
                   + f_0 * hl_554[k];

        t_330[k] = -fl_150[k]
                   + f_0 * hl_555[k];

        t_331[k] = -fl_151[k]
                   + f_0 * hl_556[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, fl_152, fl_153, fl_154, fl_155, \
                         fl_156, hl_557, hl_558, hl_559, hl_560, \
                         hl_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -fl_152[k]
                   + f_0 * hl_557[k];

        t_333[k] = -fl_153[k]
                   + f_0 * hl_558[k];

        t_334[k] = -fl_154[k]
                   + f_0 * hl_559[k];

        t_335[k] = -fl_155[k]
                   + f_0 * hl_560[k];

        t_336[k] = -fl_156[k]
                   + f_0 * hl_561[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, fl_157, fl_158, fl_159, fl_160, \
                         fl_161, hl_562, hl_563, hl_564, hl_565, \
                         hl_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -fl_157[k]
                   + f_0 * hl_562[k];

        t_338[k] = -fl_158[k]
                   + f_0 * hl_563[k];

        t_339[k] = -fl_159[k]
                   + f_0 * hl_564[k];

        t_340[k] = -fl_160[k]
                   + f_0 * hl_565[k];

        t_341[k] = -fl_161[k]
                   + f_0 * hl_566[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, fl_162, fl_163, fl_164, fl_165, \
                         fl_166, hl_567, hl_568, hl_569, hl_570, \
                         hl_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -fl_162[k]
                   + f_0 * hl_567[k];

        t_343[k] = -fl_163[k]
                   + f_0 * hl_568[k];

        t_344[k] = -fl_164[k]
                   + f_0 * hl_569[k];

        t_345[k] = -fl_165[k]
                   + f_0 * hl_570[k];

        t_346[k] = -fl_166[k]
                   + f_0 * hl_571[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, fl_167, fl_168, fl_169, fl_170, \
                         fl_171, hl_572, hl_573, hl_574, hl_575, \
                         hl_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -fl_167[k]
                   + f_0 * hl_572[k];

        t_348[k] = -fl_168[k]
                   + f_0 * hl_573[k];

        t_349[k] = -fl_169[k]
                   + f_0 * hl_574[k];

        t_350[k] = -fl_170[k]
                   + f_0 * hl_575[k];

        t_351[k] = -fl_171[k]
                   + f_0 * hl_576[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, fl_172, fl_173, fl_174, fl_175, \
                         fl_176, hl_577, hl_578, hl_579, hl_580, \
                         hl_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -fl_172[k]
                   + f_0 * hl_577[k];

        t_353[k] = -fl_173[k]
                   + f_0 * hl_578[k];

        t_354[k] = -fl_174[k]
                   + f_0 * hl_579[k];

        t_355[k] = -fl_175[k]
                   + f_0 * hl_580[k];

        t_356[k] = -fl_176[k]
                   + f_0 * hl_581[k];
    }
}

static auto
compute_prim_geom_10_gl_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *fl_177 = buffer.data(fl + 177);
    const auto *fl_178 = buffer.data(fl + 178);
    const auto *fl_179 = buffer.data(fl + 179);
    const auto *fl_180 = buffer.data(fl + 180);
    const auto *fl_181 = buffer.data(fl + 181);
    const auto *fl_182 = buffer.data(fl + 182);
    const auto *fl_183 = buffer.data(fl + 183);
    const auto *fl_184 = buffer.data(fl + 184);
    const auto *fl_185 = buffer.data(fl + 185);
    const auto *fl_186 = buffer.data(fl + 186);
    const auto *fl_187 = buffer.data(fl + 187);
    const auto *fl_188 = buffer.data(fl + 188);
    const auto *fl_189 = buffer.data(fl + 189);
    const auto *fl_190 = buffer.data(fl + 190);
    const auto *fl_191 = buffer.data(fl + 191);
    const auto *fl_192 = buffer.data(fl + 192);
    const auto *fl_193 = buffer.data(fl + 193);
    const auto *fl_194 = buffer.data(fl + 194);
    const auto *fl_195 = buffer.data(fl + 195);
    const auto *fl_196 = buffer.data(fl + 196);
    const auto *fl_197 = buffer.data(fl + 197);
    const auto *fl_198 = buffer.data(fl + 198);
    const auto *fl_199 = buffer.data(fl + 199);
    const auto *fl_200 = buffer.data(fl + 200);
    const auto *fl_201 = buffer.data(fl + 201);
    const auto *fl_202 = buffer.data(fl + 202);
    const auto *fl_203 = buffer.data(fl + 203);
    const auto *fl_204 = buffer.data(fl + 204);
    const auto *fl_205 = buffer.data(fl + 205);
    const auto *fl_206 = buffer.data(fl + 206);
    const auto *fl_207 = buffer.data(fl + 207);
    const auto *fl_208 = buffer.data(fl + 208);
    const auto *fl_209 = buffer.data(fl + 209);
    const auto *fl_210 = buffer.data(fl + 210);
    const auto *fl_211 = buffer.data(fl + 211);
    const auto *fl_212 = buffer.data(fl + 212);
    const auto *fl_213 = buffer.data(fl + 213);
    const auto *fl_214 = buffer.data(fl + 214);
    const auto *fl_215 = buffer.data(fl + 215);
    const auto *fl_216 = buffer.data(fl + 216);
    const auto *fl_217 = buffer.data(fl + 217);
    const auto *fl_218 = buffer.data(fl + 218);
    const auto *fl_219 = buffer.data(fl + 219);
    const auto *fl_220 = buffer.data(fl + 220);
    const auto *fl_221 = buffer.data(fl + 221);
    const auto *fl_222 = buffer.data(fl + 222);
    const auto *fl_223 = buffer.data(fl + 223);
    const auto *fl_224 = buffer.data(fl + 224);
    const auto *fl_225 = buffer.data(fl + 225);
    const auto *fl_226 = buffer.data(fl + 226);
    const auto *fl_227 = buffer.data(fl + 227);
    const auto *fl_228 = buffer.data(fl + 228);
    const auto *fl_229 = buffer.data(fl + 229);
    const auto *fl_230 = buffer.data(fl + 230);
    const auto *fl_231 = buffer.data(fl + 231);
    const auto *fl_232 = buffer.data(fl + 232);
    const auto *fl_233 = buffer.data(fl + 233);
    const auto *fl_234 = buffer.data(fl + 234);
    const auto *fl_235 = buffer.data(fl + 235);
    const auto *fl_236 = buffer.data(fl + 236);
    const auto *fl_237 = buffer.data(fl + 237);
    const auto *fl_238 = buffer.data(fl + 238);
    const auto *fl_239 = buffer.data(fl + 239);
    const auto *fl_240 = buffer.data(fl + 240);
    const auto *fl_241 = buffer.data(fl + 241);
    const auto *fl_242 = buffer.data(fl + 242);
    const auto *fl_243 = buffer.data(fl + 243);
    const auto *fl_244 = buffer.data(fl + 244);
    const auto *fl_245 = buffer.data(fl + 245);
    const auto *fl_246 = buffer.data(fl + 246);
    const auto *fl_247 = buffer.data(fl + 247);
    const auto *fl_248 = buffer.data(fl + 248);
    const auto *fl_249 = buffer.data(fl + 249);
    const auto *fl_250 = buffer.data(fl + 250);
    const auto *fl_251 = buffer.data(fl + 251);
    const auto *fl_252 = buffer.data(fl + 252);
    const auto *fl_253 = buffer.data(fl + 253);
    const auto *fl_254 = buffer.data(fl + 254);
    const auto *fl_255 = buffer.data(fl + 255);
    const auto *fl_256 = buffer.data(fl + 256);
    const auto *fl_257 = buffer.data(fl + 257);
    const auto *fl_258 = buffer.data(fl + 258);
    const auto *fl_259 = buffer.data(fl + 259);
    const auto *fl_260 = buffer.data(fl + 260);
    const auto *fl_261 = buffer.data(fl + 261);
    const auto *fl_262 = buffer.data(fl + 262);
    const auto *fl_263 = buffer.data(fl + 263);
    const auto *fl_264 = buffer.data(fl + 264);
    const auto *fl_265 = buffer.data(fl + 265);
    const auto *fl_266 = buffer.data(fl + 266);
    const auto *fl_267 = buffer.data(fl + 267);
    const auto *fl_268 = buffer.data(fl + 268);
    const auto *fl_269 = buffer.data(fl + 269);
    const auto *fl_270 = buffer.data(fl + 270);
    const auto *fl_271 = buffer.data(fl + 271);
    const auto *fl_272 = buffer.data(fl + 272);
    const auto *fl_273 = buffer.data(fl + 273);
    const auto *fl_274 = buffer.data(fl + 274);
    const auto *fl_275 = buffer.data(fl + 275);
    const auto *fl_276 = buffer.data(fl + 276);
    const auto *fl_277 = buffer.data(fl + 277);
    const auto *fl_278 = buffer.data(fl + 278);
    const auto *fl_279 = buffer.data(fl + 279);
    const auto *fl_280 = buffer.data(fl + 280);
    const auto *fl_281 = buffer.data(fl + 281);
    const auto *fl_282 = buffer.data(fl + 282);
    const auto *fl_283 = buffer.data(fl + 283);
    const auto *fl_284 = buffer.data(fl + 284);
    const auto *fl_285 = buffer.data(fl + 285);
    const auto *fl_286 = buffer.data(fl + 286);
    const auto *fl_287 = buffer.data(fl + 287);
    const auto *fl_288 = buffer.data(fl + 288);
    const auto *fl_289 = buffer.data(fl + 289);
    const auto *fl_290 = buffer.data(fl + 290);
    const auto *fl_291 = buffer.data(fl + 291);
    const auto *fl_292 = buffer.data(fl + 292);
    const auto *fl_293 = buffer.data(fl + 293);

    const auto *hl_582 = buffer.data(hl + 582);
    const auto *hl_583 = buffer.data(hl + 583);
    const auto *hl_584 = buffer.data(hl + 584);
    const auto *hl_585 = buffer.data(hl + 585);
    const auto *hl_586 = buffer.data(hl + 586);
    const auto *hl_587 = buffer.data(hl + 587);
    const auto *hl_588 = buffer.data(hl + 588);
    const auto *hl_589 = buffer.data(hl + 589);
    const auto *hl_590 = buffer.data(hl + 590);
    const auto *hl_591 = buffer.data(hl + 591);
    const auto *hl_592 = buffer.data(hl + 592);
    const auto *hl_593 = buffer.data(hl + 593);
    const auto *hl_594 = buffer.data(hl + 594);
    const auto *hl_595 = buffer.data(hl + 595);
    const auto *hl_596 = buffer.data(hl + 596);
    const auto *hl_597 = buffer.data(hl + 597);
    const auto *hl_598 = buffer.data(hl + 598);
    const auto *hl_599 = buffer.data(hl + 599);
    const auto *hl_600 = buffer.data(hl + 600);
    const auto *hl_601 = buffer.data(hl + 601);
    const auto *hl_602 = buffer.data(hl + 602);
    const auto *hl_603 = buffer.data(hl + 603);
    const auto *hl_604 = buffer.data(hl + 604);
    const auto *hl_605 = buffer.data(hl + 605);
    const auto *hl_606 = buffer.data(hl + 606);
    const auto *hl_607 = buffer.data(hl + 607);
    const auto *hl_608 = buffer.data(hl + 608);
    const auto *hl_609 = buffer.data(hl + 609);
    const auto *hl_610 = buffer.data(hl + 610);
    const auto *hl_611 = buffer.data(hl + 611);
    const auto *hl_612 = buffer.data(hl + 612);
    const auto *hl_613 = buffer.data(hl + 613);
    const auto *hl_614 = buffer.data(hl + 614);
    const auto *hl_615 = buffer.data(hl + 615);
    const auto *hl_616 = buffer.data(hl + 616);
    const auto *hl_617 = buffer.data(hl + 617);
    const auto *hl_618 = buffer.data(hl + 618);
    const auto *hl_619 = buffer.data(hl + 619);
    const auto *hl_620 = buffer.data(hl + 620);
    const auto *hl_621 = buffer.data(hl + 621);
    const auto *hl_622 = buffer.data(hl + 622);
    const auto *hl_623 = buffer.data(hl + 623);
    const auto *hl_624 = buffer.data(hl + 624);
    const auto *hl_625 = buffer.data(hl + 625);
    const auto *hl_626 = buffer.data(hl + 626);
    const auto *hl_627 = buffer.data(hl + 627);
    const auto *hl_628 = buffer.data(hl + 628);
    const auto *hl_629 = buffer.data(hl + 629);
    const auto *hl_630 = buffer.data(hl + 630);
    const auto *hl_631 = buffer.data(hl + 631);
    const auto *hl_632 = buffer.data(hl + 632);
    const auto *hl_633 = buffer.data(hl + 633);
    const auto *hl_634 = buffer.data(hl + 634);
    const auto *hl_635 = buffer.data(hl + 635);
    const auto *hl_636 = buffer.data(hl + 636);
    const auto *hl_637 = buffer.data(hl + 637);
    const auto *hl_638 = buffer.data(hl + 638);
    const auto *hl_639 = buffer.data(hl + 639);
    const auto *hl_640 = buffer.data(hl + 640);
    const auto *hl_641 = buffer.data(hl + 641);
    const auto *hl_642 = buffer.data(hl + 642);
    const auto *hl_643 = buffer.data(hl + 643);
    const auto *hl_644 = buffer.data(hl + 644);
    const auto *hl_645 = buffer.data(hl + 645);
    const auto *hl_646 = buffer.data(hl + 646);
    const auto *hl_647 = buffer.data(hl + 647);
    const auto *hl_648 = buffer.data(hl + 648);
    const auto *hl_649 = buffer.data(hl + 649);
    const auto *hl_650 = buffer.data(hl + 650);
    const auto *hl_651 = buffer.data(hl + 651);
    const auto *hl_652 = buffer.data(hl + 652);
    const auto *hl_653 = buffer.data(hl + 653);
    const auto *hl_654 = buffer.data(hl + 654);
    const auto *hl_655 = buffer.data(hl + 655);
    const auto *hl_656 = buffer.data(hl + 656);
    const auto *hl_657 = buffer.data(hl + 657);
    const auto *hl_658 = buffer.data(hl + 658);
    const auto *hl_659 = buffer.data(hl + 659);
    const auto *hl_660 = buffer.data(hl + 660);
    const auto *hl_661 = buffer.data(hl + 661);
    const auto *hl_662 = buffer.data(hl + 662);
    const auto *hl_663 = buffer.data(hl + 663);
    const auto *hl_664 = buffer.data(hl + 664);
    const auto *hl_665 = buffer.data(hl + 665);
    const auto *hl_666 = buffer.data(hl + 666);
    const auto *hl_667 = buffer.data(hl + 667);
    const auto *hl_668 = buffer.data(hl + 668);
    const auto *hl_669 = buffer.data(hl + 669);
    const auto *hl_670 = buffer.data(hl + 670);
    const auto *hl_671 = buffer.data(hl + 671);
    const auto *hl_672 = buffer.data(hl + 672);
    const auto *hl_673 = buffer.data(hl + 673);
    const auto *hl_674 = buffer.data(hl + 674);
    const auto *hl_720 = buffer.data(hl + 720);
    const auto *hl_721 = buffer.data(hl + 721);
    const auto *hl_722 = buffer.data(hl + 722);
    const auto *hl_723 = buffer.data(hl + 723);
    const auto *hl_724 = buffer.data(hl + 724);
    const auto *hl_725 = buffer.data(hl + 725);
    const auto *hl_726 = buffer.data(hl + 726);
    const auto *hl_727 = buffer.data(hl + 727);
    const auto *hl_728 = buffer.data(hl + 728);
    const auto *hl_729 = buffer.data(hl + 729);
    const auto *hl_730 = buffer.data(hl + 730);
    const auto *hl_731 = buffer.data(hl + 731);
    const auto *hl_732 = buffer.data(hl + 732);
    const auto *hl_733 = buffer.data(hl + 733);
    const auto *hl_734 = buffer.data(hl + 734);
    const auto *hl_735 = buffer.data(hl + 735);
    const auto *hl_736 = buffer.data(hl + 736);
    const auto *hl_737 = buffer.data(hl + 737);
    const auto *hl_738 = buffer.data(hl + 738);
    const auto *hl_739 = buffer.data(hl + 739);
    const auto *hl_740 = buffer.data(hl + 740);
    const auto *hl_741 = buffer.data(hl + 741);
    const auto *hl_742 = buffer.data(hl + 742);
    const auto *hl_743 = buffer.data(hl + 743);
    const auto *hl_744 = buffer.data(hl + 744);
    const auto *hl_745 = buffer.data(hl + 745);
    const auto *hl_746 = buffer.data(hl + 746);
    const auto *hl_747 = buffer.data(hl + 747);
    const auto *hl_748 = buffer.data(hl + 748);
    const auto *hl_749 = buffer.data(hl + 749);
    const auto *hl_750 = buffer.data(hl + 750);
    const auto *hl_751 = buffer.data(hl + 751);
    const auto *hl_752 = buffer.data(hl + 752);
    const auto *hl_753 = buffer.data(hl + 753);
    const auto *hl_754 = buffer.data(hl + 754);
    const auto *hl_755 = buffer.data(hl + 755);
    const auto *hl_756 = buffer.data(hl + 756);
    const auto *hl_757 = buffer.data(hl + 757);
    const auto *hl_758 = buffer.data(hl + 758);
    const auto *hl_759 = buffer.data(hl + 759);
    const auto *hl_760 = buffer.data(hl + 760);
    const auto *hl_761 = buffer.data(hl + 761);
    const auto *hl_762 = buffer.data(hl + 762);
    const auto *hl_763 = buffer.data(hl + 763);
    const auto *hl_764 = buffer.data(hl + 764);
    const auto *hl_765 = buffer.data(hl + 765);
    const auto *hl_766 = buffer.data(hl + 766);
    const auto *hl_767 = buffer.data(hl + 767);
    const auto *hl_768 = buffer.data(hl + 768);
    const auto *hl_769 = buffer.data(hl + 769);
    const auto *hl_770 = buffer.data(hl + 770);
    const auto *hl_771 = buffer.data(hl + 771);
    const auto *hl_772 = buffer.data(hl + 772);
    const auto *hl_773 = buffer.data(hl + 773);
    const auto *hl_774 = buffer.data(hl + 774);
    const auto *hl_775 = buffer.data(hl + 775);
    const auto *hl_776 = buffer.data(hl + 776);
    const auto *hl_777 = buffer.data(hl + 777);
    const auto *hl_778 = buffer.data(hl + 778);
    const auto *hl_779 = buffer.data(hl + 779);
    const auto *hl_780 = buffer.data(hl + 780);
    const auto *hl_781 = buffer.data(hl + 781);
    const auto *hl_782 = buffer.data(hl + 782);
    const auto *hl_783 = buffer.data(hl + 783);
    const auto *hl_784 = buffer.data(hl + 784);
    const auto *hl_785 = buffer.data(hl + 785);
    const auto *hl_786 = buffer.data(hl + 786);
    const auto *hl_787 = buffer.data(hl + 787);
    const auto *hl_788 = buffer.data(hl + 788);

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, fl_177, fl_178, fl_179, fl_180, \
                         fl_181, hl_582, hl_583, hl_584, hl_585, \
                         hl_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -fl_177[k]
                   + f_0 * hl_582[k];

        t_358[k] = -fl_178[k]
                   + f_0 * hl_583[k];

        t_359[k] = -fl_179[k]
                   + f_0 * hl_584[k];

        t_360[k] = -2.0 * fl_180[k]
                   + f_0 * hl_585[k];

        t_361[k] = -2.0 * fl_181[k]
                   + f_0 * hl_586[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, fl_182, fl_183, fl_184, fl_185, \
                         fl_186, hl_587, hl_588, hl_589, hl_590, \
                         hl_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -2.0 * fl_182[k]
                   + f_0 * hl_587[k];

        t_363[k] = -2.0 * fl_183[k]
                   + f_0 * hl_588[k];

        t_364[k] = -2.0 * fl_184[k]
                   + f_0 * hl_589[k];

        t_365[k] = -2.0 * fl_185[k]
                   + f_0 * hl_590[k];

        t_366[k] = -2.0 * fl_186[k]
                   + f_0 * hl_591[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, fl_187, fl_188, fl_189, fl_190, \
                         fl_191, hl_592, hl_593, hl_594, hl_595, \
                         hl_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -2.0 * fl_187[k]
                   + f_0 * hl_592[k];

        t_368[k] = -2.0 * fl_188[k]
                   + f_0 * hl_593[k];

        t_369[k] = -2.0 * fl_189[k]
                   + f_0 * hl_594[k];

        t_370[k] = -2.0 * fl_190[k]
                   + f_0 * hl_595[k];

        t_371[k] = -2.0 * fl_191[k]
                   + f_0 * hl_596[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, fl_192, fl_193, fl_194, fl_195, \
                         fl_196, hl_597, hl_598, hl_599, hl_600, \
                         hl_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -2.0 * fl_192[k]
                   + f_0 * hl_597[k];

        t_373[k] = -2.0 * fl_193[k]
                   + f_0 * hl_598[k];

        t_374[k] = -2.0 * fl_194[k]
                   + f_0 * hl_599[k];

        t_375[k] = -2.0 * fl_195[k]
                   + f_0 * hl_600[k];

        t_376[k] = -2.0 * fl_196[k]
                   + f_0 * hl_601[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, fl_197, fl_198, fl_199, fl_200, \
                         fl_201, hl_602, hl_603, hl_604, hl_605, \
                         hl_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -2.0 * fl_197[k]
                   + f_0 * hl_602[k];

        t_378[k] = -2.0 * fl_198[k]
                   + f_0 * hl_603[k];

        t_379[k] = -2.0 * fl_199[k]
                   + f_0 * hl_604[k];

        t_380[k] = -2.0 * fl_200[k]
                   + f_0 * hl_605[k];

        t_381[k] = -2.0 * fl_201[k]
                   + f_0 * hl_606[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, fl_202, fl_203, fl_204, fl_205, \
                         fl_206, hl_607, hl_608, hl_609, hl_610, \
                         hl_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -2.0 * fl_202[k]
                   + f_0 * hl_607[k];

        t_383[k] = -2.0 * fl_203[k]
                   + f_0 * hl_608[k];

        t_384[k] = -2.0 * fl_204[k]
                   + f_0 * hl_609[k];

        t_385[k] = -2.0 * fl_205[k]
                   + f_0 * hl_610[k];

        t_386[k] = -2.0 * fl_206[k]
                   + f_0 * hl_611[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, fl_207, fl_208, fl_209, fl_210, \
                         fl_211, hl_612, hl_613, hl_614, hl_615, \
                         hl_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -2.0 * fl_207[k]
                   + f_0 * hl_612[k];

        t_388[k] = -2.0 * fl_208[k]
                   + f_0 * hl_613[k];

        t_389[k] = -2.0 * fl_209[k]
                   + f_0 * hl_614[k];

        t_390[k] = -2.0 * fl_210[k]
                   + f_0 * hl_615[k];

        t_391[k] = -2.0 * fl_211[k]
                   + f_0 * hl_616[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, fl_212, fl_213, fl_214, fl_215, \
                         fl_216, hl_617, hl_618, hl_619, hl_620, \
                         hl_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -2.0 * fl_212[k]
                   + f_0 * hl_617[k];

        t_393[k] = -2.0 * fl_213[k]
                   + f_0 * hl_618[k];

        t_394[k] = -2.0 * fl_214[k]
                   + f_0 * hl_619[k];

        t_395[k] = -2.0 * fl_215[k]
                   + f_0 * hl_620[k];

        t_396[k] = -2.0 * fl_216[k]
                   + f_0 * hl_621[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, fl_217, fl_218, fl_219, fl_220, \
                         fl_221, hl_622, hl_623, hl_624, hl_625, \
                         hl_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -2.0 * fl_217[k]
                   + f_0 * hl_622[k];

        t_398[k] = -2.0 * fl_218[k]
                   + f_0 * hl_623[k];

        t_399[k] = -2.0 * fl_219[k]
                   + f_0 * hl_624[k];

        t_400[k] = -2.0 * fl_220[k]
                   + f_0 * hl_625[k];

        t_401[k] = -2.0 * fl_221[k]
                   + f_0 * hl_626[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, fl_222, fl_223, fl_224, fl_225, \
                         fl_226, hl_627, hl_628, hl_629, hl_630, \
                         hl_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -2.0 * fl_222[k]
                   + f_0 * hl_627[k];

        t_403[k] = -2.0 * fl_223[k]
                   + f_0 * hl_628[k];

        t_404[k] = -2.0 * fl_224[k]
                   + f_0 * hl_629[k];

        t_405[k] = -3.0 * fl_225[k]
                   + f_0 * hl_630[k];

        t_406[k] = -3.0 * fl_226[k]
                   + f_0 * hl_631[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, fl_227, fl_228, fl_229, fl_230, \
                         fl_231, hl_632, hl_633, hl_634, hl_635, \
                         hl_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -3.0 * fl_227[k]
                   + f_0 * hl_632[k];

        t_408[k] = -3.0 * fl_228[k]
                   + f_0 * hl_633[k];

        t_409[k] = -3.0 * fl_229[k]
                   + f_0 * hl_634[k];

        t_410[k] = -3.0 * fl_230[k]
                   + f_0 * hl_635[k];

        t_411[k] = -3.0 * fl_231[k]
                   + f_0 * hl_636[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, fl_232, fl_233, fl_234, fl_235, \
                         fl_236, hl_637, hl_638, hl_639, hl_640, \
                         hl_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -3.0 * fl_232[k]
                   + f_0 * hl_637[k];

        t_413[k] = -3.0 * fl_233[k]
                   + f_0 * hl_638[k];

        t_414[k] = -3.0 * fl_234[k]
                   + f_0 * hl_639[k];

        t_415[k] = -3.0 * fl_235[k]
                   + f_0 * hl_640[k];

        t_416[k] = -3.0 * fl_236[k]
                   + f_0 * hl_641[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, fl_237, fl_238, fl_239, fl_240, \
                         fl_241, hl_642, hl_643, hl_644, hl_645, \
                         hl_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -3.0 * fl_237[k]
                   + f_0 * hl_642[k];

        t_418[k] = -3.0 * fl_238[k]
                   + f_0 * hl_643[k];

        t_419[k] = -3.0 * fl_239[k]
                   + f_0 * hl_644[k];

        t_420[k] = -3.0 * fl_240[k]
                   + f_0 * hl_645[k];

        t_421[k] = -3.0 * fl_241[k]
                   + f_0 * hl_646[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, t_426, fl_242, fl_243, fl_244, fl_245, \
                         fl_246, hl_647, hl_648, hl_649, hl_650, \
                         hl_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = -3.0 * fl_242[k]
                   + f_0 * hl_647[k];

        t_423[k] = -3.0 * fl_243[k]
                   + f_0 * hl_648[k];

        t_424[k] = -3.0 * fl_244[k]
                   + f_0 * hl_649[k];

        t_425[k] = -3.0 * fl_245[k]
                   + f_0 * hl_650[k];

        t_426[k] = -3.0 * fl_246[k]
                   + f_0 * hl_651[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, t_431, fl_247, fl_248, fl_249, fl_250, \
                         fl_251, hl_652, hl_653, hl_654, hl_655, \
                         hl_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = -3.0 * fl_247[k]
                   + f_0 * hl_652[k];

        t_428[k] = -3.0 * fl_248[k]
                   + f_0 * hl_653[k];

        t_429[k] = -3.0 * fl_249[k]
                   + f_0 * hl_654[k];

        t_430[k] = -3.0 * fl_250[k]
                   + f_0 * hl_655[k];

        t_431[k] = -3.0 * fl_251[k]
                   + f_0 * hl_656[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, t_436, fl_252, fl_253, fl_254, fl_255, \
                         fl_256, hl_657, hl_658, hl_659, hl_660, \
                         hl_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = -3.0 * fl_252[k]
                   + f_0 * hl_657[k];

        t_433[k] = -3.0 * fl_253[k]
                   + f_0 * hl_658[k];

        t_434[k] = -3.0 * fl_254[k]
                   + f_0 * hl_659[k];

        t_435[k] = -3.0 * fl_255[k]
                   + f_0 * hl_660[k];

        t_436[k] = -3.0 * fl_256[k]
                   + f_0 * hl_661[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, t_441, fl_257, fl_258, fl_259, fl_260, \
                         fl_261, hl_662, hl_663, hl_664, hl_665, \
                         hl_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = -3.0 * fl_257[k]
                   + f_0 * hl_662[k];

        t_438[k] = -3.0 * fl_258[k]
                   + f_0 * hl_663[k];

        t_439[k] = -3.0 * fl_259[k]
                   + f_0 * hl_664[k];

        t_440[k] = -3.0 * fl_260[k]
                   + f_0 * hl_665[k];

        t_441[k] = -3.0 * fl_261[k]
                   + f_0 * hl_666[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, t_446, fl_262, fl_263, fl_264, fl_265, \
                         fl_266, hl_667, hl_668, hl_669, hl_670, \
                         hl_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = -3.0 * fl_262[k]
                   + f_0 * hl_667[k];

        t_443[k] = -3.0 * fl_263[k]
                   + f_0 * hl_668[k];

        t_444[k] = -3.0 * fl_264[k]
                   + f_0 * hl_669[k];

        t_445[k] = -3.0 * fl_265[k]
                   + f_0 * hl_670[k];

        t_446[k] = -3.0 * fl_266[k]
                   + f_0 * hl_671[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, t_451, t_452, fl_267, fl_268, fl_269, \
                         hl_672, hl_673, hl_674, hl_720, hl_721, \
                         hl_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = -3.0 * fl_267[k]
                   + f_0 * hl_672[k];

        t_448[k] = -3.0 * fl_268[k]
                   + f_0 * hl_673[k];

        t_449[k] = -3.0 * fl_269[k]
                   + f_0 * hl_674[k];

        t_450[k] = f_0 * hl_720[k];

        t_451[k] = f_0 * hl_721[k];

        t_452[k] = f_0 * hl_722[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, t_457, t_458, t_459, t_460, hl_723, \
                         hl_724, hl_725, hl_726, hl_727, hl_728, hl_729, \
                         hl_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_0 * hl_723[k];

        t_454[k] = f_0 * hl_724[k];

        t_455[k] = f_0 * hl_725[k];

        t_456[k] = f_0 * hl_726[k];

        t_457[k] = f_0 * hl_727[k];

        t_458[k] = f_0 * hl_728[k];

        t_459[k] = f_0 * hl_729[k];

        t_460[k] = f_0 * hl_730[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, t_465, t_466, t_467, t_468, hl_731, \
                         hl_732, hl_733, hl_734, hl_735, hl_736, hl_737, \
                         hl_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_0 * hl_731[k];

        t_462[k] = f_0 * hl_732[k];

        t_463[k] = f_0 * hl_733[k];

        t_464[k] = f_0 * hl_734[k];

        t_465[k] = f_0 * hl_735[k];

        t_466[k] = f_0 * hl_736[k];

        t_467[k] = f_0 * hl_737[k];

        t_468[k] = f_0 * hl_738[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, t_473, t_474, t_475, t_476, hl_739, \
                         hl_740, hl_741, hl_742, hl_743, hl_744, hl_745, \
                         hl_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_0 * hl_739[k];

        t_470[k] = f_0 * hl_740[k];

        t_471[k] = f_0 * hl_741[k];

        t_472[k] = f_0 * hl_742[k];

        t_473[k] = f_0 * hl_743[k];

        t_474[k] = f_0 * hl_744[k];

        t_475[k] = f_0 * hl_745[k];

        t_476[k] = f_0 * hl_746[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, t_482, t_483, t_484, hl_747, \
                         hl_748, hl_749, hl_750, hl_751, hl_752, hl_753, \
                         hl_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_0 * hl_747[k];

        t_478[k] = f_0 * hl_748[k];

        t_479[k] = f_0 * hl_749[k];

        t_480[k] = f_0 * hl_750[k];

        t_481[k] = f_0 * hl_751[k];

        t_482[k] = f_0 * hl_752[k];

        t_483[k] = f_0 * hl_753[k];

        t_484[k] = f_0 * hl_754[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, t_490, t_491, t_492, hl_755, \
                         hl_756, hl_757, hl_758, hl_759, hl_760, hl_761, \
                         hl_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_0 * hl_755[k];

        t_486[k] = f_0 * hl_756[k];

        t_487[k] = f_0 * hl_757[k];

        t_488[k] = f_0 * hl_758[k];

        t_489[k] = f_0 * hl_759[k];

        t_490[k] = f_0 * hl_760[k];

        t_491[k] = f_0 * hl_761[k];

        t_492[k] = f_0 * hl_762[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, t_496, t_497, t_498, fl_270, fl_271, fl_272, \
                         fl_273, hl_763, hl_764, hl_765, hl_766, hl_767, \
                         hl_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = f_0 * hl_763[k];

        t_494[k] = f_0 * hl_764[k];

        t_495[k] = -fl_270[k]
                   + f_0 * hl_765[k];

        t_496[k] = -fl_271[k]
                   + f_0 * hl_766[k];

        t_497[k] = -fl_272[k]
                   + f_0 * hl_767[k];

        t_498[k] = -fl_273[k]
                   + f_0 * hl_768[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, t_503, fl_274, fl_275, fl_276, fl_277, \
                         fl_278, hl_769, hl_770, hl_771, hl_772, \
                         hl_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = -fl_274[k]
                   + f_0 * hl_769[k];

        t_500[k] = -fl_275[k]
                   + f_0 * hl_770[k];

        t_501[k] = -fl_276[k]
                   + f_0 * hl_771[k];

        t_502[k] = -fl_277[k]
                   + f_0 * hl_772[k];

        t_503[k] = -fl_278[k]
                   + f_0 * hl_773[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, fl_279, fl_280, fl_281, fl_282, \
                         fl_283, hl_774, hl_775, hl_776, hl_777, \
                         hl_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = -fl_279[k]
                   + f_0 * hl_774[k];

        t_505[k] = -fl_280[k]
                   + f_0 * hl_775[k];

        t_506[k] = -fl_281[k]
                   + f_0 * hl_776[k];

        t_507[k] = -fl_282[k]
                   + f_0 * hl_777[k];

        t_508[k] = -fl_283[k]
                   + f_0 * hl_778[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, t_513, fl_284, fl_285, fl_286, fl_287, \
                         fl_288, hl_779, hl_780, hl_781, hl_782, \
                         hl_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = -fl_284[k]
                   + f_0 * hl_779[k];

        t_510[k] = -fl_285[k]
                   + f_0 * hl_780[k];

        t_511[k] = -fl_286[k]
                   + f_0 * hl_781[k];

        t_512[k] = -fl_287[k]
                   + f_0 * hl_782[k];

        t_513[k] = -fl_288[k]
                   + f_0 * hl_783[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, fl_289, fl_290, fl_291, fl_292, \
                         fl_293, hl_784, hl_785, hl_786, hl_787, \
                         hl_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = -fl_289[k]
                   + f_0 * hl_784[k];

        t_515[k] = -fl_290[k]
                   + f_0 * hl_785[k];

        t_516[k] = -fl_291[k]
                   + f_0 * hl_786[k];

        t_517[k] = -fl_292[k]
                   + f_0 * hl_787[k];

        t_518[k] = -fl_293[k]
                   + f_0 * hl_788[k];
    }
}

static auto
compute_prim_geom_10_gl_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);

    const auto *fl_294 = buffer.data(fl + 294);
    const auto *fl_295 = buffer.data(fl + 295);
    const auto *fl_296 = buffer.data(fl + 296);
    const auto *fl_297 = buffer.data(fl + 297);
    const auto *fl_298 = buffer.data(fl + 298);
    const auto *fl_299 = buffer.data(fl + 299);
    const auto *fl_300 = buffer.data(fl + 300);
    const auto *fl_301 = buffer.data(fl + 301);
    const auto *fl_302 = buffer.data(fl + 302);
    const auto *fl_303 = buffer.data(fl + 303);
    const auto *fl_304 = buffer.data(fl + 304);
    const auto *fl_305 = buffer.data(fl + 305);
    const auto *fl_306 = buffer.data(fl + 306);
    const auto *fl_307 = buffer.data(fl + 307);
    const auto *fl_308 = buffer.data(fl + 308);
    const auto *fl_309 = buffer.data(fl + 309);
    const auto *fl_310 = buffer.data(fl + 310);
    const auto *fl_311 = buffer.data(fl + 311);
    const auto *fl_312 = buffer.data(fl + 312);
    const auto *fl_313 = buffer.data(fl + 313);
    const auto *fl_314 = buffer.data(fl + 314);
    const auto *fl_315 = buffer.data(fl + 315);
    const auto *fl_316 = buffer.data(fl + 316);
    const auto *fl_317 = buffer.data(fl + 317);
    const auto *fl_318 = buffer.data(fl + 318);
    const auto *fl_319 = buffer.data(fl + 319);
    const auto *fl_320 = buffer.data(fl + 320);
    const auto *fl_321 = buffer.data(fl + 321);
    const auto *fl_322 = buffer.data(fl + 322);
    const auto *fl_323 = buffer.data(fl + 323);
    const auto *fl_324 = buffer.data(fl + 324);
    const auto *fl_325 = buffer.data(fl + 325);
    const auto *fl_326 = buffer.data(fl + 326);
    const auto *fl_327 = buffer.data(fl + 327);
    const auto *fl_328 = buffer.data(fl + 328);
    const auto *fl_329 = buffer.data(fl + 329);
    const auto *fl_330 = buffer.data(fl + 330);
    const auto *fl_331 = buffer.data(fl + 331);
    const auto *fl_332 = buffer.data(fl + 332);
    const auto *fl_333 = buffer.data(fl + 333);
    const auto *fl_334 = buffer.data(fl + 334);
    const auto *fl_335 = buffer.data(fl + 335);
    const auto *fl_336 = buffer.data(fl + 336);
    const auto *fl_337 = buffer.data(fl + 337);
    const auto *fl_338 = buffer.data(fl + 338);
    const auto *fl_339 = buffer.data(fl + 339);
    const auto *fl_340 = buffer.data(fl + 340);
    const auto *fl_341 = buffer.data(fl + 341);
    const auto *fl_342 = buffer.data(fl + 342);
    const auto *fl_343 = buffer.data(fl + 343);
    const auto *fl_344 = buffer.data(fl + 344);
    const auto *fl_345 = buffer.data(fl + 345);
    const auto *fl_346 = buffer.data(fl + 346);
    const auto *fl_347 = buffer.data(fl + 347);
    const auto *fl_348 = buffer.data(fl + 348);
    const auto *fl_349 = buffer.data(fl + 349);
    const auto *fl_350 = buffer.data(fl + 350);
    const auto *fl_351 = buffer.data(fl + 351);
    const auto *fl_352 = buffer.data(fl + 352);
    const auto *fl_353 = buffer.data(fl + 353);
    const auto *fl_354 = buffer.data(fl + 354);
    const auto *fl_355 = buffer.data(fl + 355);
    const auto *fl_356 = buffer.data(fl + 356);
    const auto *fl_357 = buffer.data(fl + 357);
    const auto *fl_358 = buffer.data(fl + 358);
    const auto *fl_359 = buffer.data(fl + 359);
    const auto *fl_360 = buffer.data(fl + 360);
    const auto *fl_361 = buffer.data(fl + 361);
    const auto *fl_362 = buffer.data(fl + 362);
    const auto *fl_363 = buffer.data(fl + 363);
    const auto *fl_364 = buffer.data(fl + 364);
    const auto *fl_365 = buffer.data(fl + 365);
    const auto *fl_366 = buffer.data(fl + 366);
    const auto *fl_367 = buffer.data(fl + 367);
    const auto *fl_368 = buffer.data(fl + 368);
    const auto *fl_369 = buffer.data(fl + 369);
    const auto *fl_370 = buffer.data(fl + 370);
    const auto *fl_371 = buffer.data(fl + 371);
    const auto *fl_372 = buffer.data(fl + 372);
    const auto *fl_373 = buffer.data(fl + 373);
    const auto *fl_374 = buffer.data(fl + 374);
    const auto *fl_375 = buffer.data(fl + 375);
    const auto *fl_376 = buffer.data(fl + 376);
    const auto *fl_377 = buffer.data(fl + 377);
    const auto *fl_378 = buffer.data(fl + 378);
    const auto *fl_379 = buffer.data(fl + 379);
    const auto *fl_380 = buffer.data(fl + 380);
    const auto *fl_381 = buffer.data(fl + 381);
    const auto *fl_382 = buffer.data(fl + 382);
    const auto *fl_383 = buffer.data(fl + 383);
    const auto *fl_384 = buffer.data(fl + 384);
    const auto *fl_385 = buffer.data(fl + 385);
    const auto *fl_386 = buffer.data(fl + 386);
    const auto *fl_387 = buffer.data(fl + 387);
    const auto *fl_388 = buffer.data(fl + 388);
    const auto *fl_389 = buffer.data(fl + 389);
    const auto *fl_390 = buffer.data(fl + 390);
    const auto *fl_391 = buffer.data(fl + 391);
    const auto *fl_392 = buffer.data(fl + 392);
    const auto *fl_393 = buffer.data(fl + 393);
    const auto *fl_394 = buffer.data(fl + 394);
    const auto *fl_395 = buffer.data(fl + 395);
    const auto *fl_396 = buffer.data(fl + 396);
    const auto *fl_397 = buffer.data(fl + 397);
    const auto *fl_398 = buffer.data(fl + 398);
    const auto *fl_399 = buffer.data(fl + 399);
    const auto *fl_400 = buffer.data(fl + 400);
    const auto *fl_401 = buffer.data(fl + 401);
    const auto *fl_402 = buffer.data(fl + 402);
    const auto *fl_403 = buffer.data(fl + 403);
    const auto *fl_404 = buffer.data(fl + 404);
    const auto *fl_405 = buffer.data(fl + 405);
    const auto *fl_406 = buffer.data(fl + 406);
    const auto *fl_407 = buffer.data(fl + 407);
    const auto *fl_408 = buffer.data(fl + 408);
    const auto *fl_409 = buffer.data(fl + 409);
    const auto *fl_410 = buffer.data(fl + 410);
    const auto *fl_411 = buffer.data(fl + 411);
    const auto *fl_412 = buffer.data(fl + 412);
    const auto *fl_413 = buffer.data(fl + 413);
    const auto *fl_414 = buffer.data(fl + 414);
    const auto *fl_415 = buffer.data(fl + 415);
    const auto *fl_416 = buffer.data(fl + 416);
    const auto *fl_417 = buffer.data(fl + 417);
    const auto *fl_418 = buffer.data(fl + 418);
    const auto *fl_419 = buffer.data(fl + 419);
    const auto *fl_420 = buffer.data(fl + 420);
    const auto *fl_421 = buffer.data(fl + 421);
    const auto *fl_422 = buffer.data(fl + 422);
    const auto *fl_423 = buffer.data(fl + 423);
    const auto *fl_424 = buffer.data(fl + 424);
    const auto *fl_425 = buffer.data(fl + 425);
    const auto *fl_426 = buffer.data(fl + 426);
    const auto *fl_427 = buffer.data(fl + 427);
    const auto *fl_428 = buffer.data(fl + 428);
    const auto *fl_429 = buffer.data(fl + 429);
    const auto *fl_430 = buffer.data(fl + 430);
    const auto *fl_431 = buffer.data(fl + 431);
    const auto *fl_432 = buffer.data(fl + 432);
    const auto *fl_433 = buffer.data(fl + 433);
    const auto *fl_434 = buffer.data(fl + 434);
    const auto *fl_435 = buffer.data(fl + 435);
    const auto *fl_436 = buffer.data(fl + 436);
    const auto *fl_437 = buffer.data(fl + 437);
    const auto *fl_438 = buffer.data(fl + 438);
    const auto *fl_439 = buffer.data(fl + 439);
    const auto *fl_440 = buffer.data(fl + 440);
    const auto *fl_441 = buffer.data(fl + 441);
    const auto *fl_442 = buffer.data(fl + 442);
    const auto *fl_443 = buffer.data(fl + 443);

    const auto *hl_789 = buffer.data(hl + 789);
    const auto *hl_790 = buffer.data(hl + 790);
    const auto *hl_791 = buffer.data(hl + 791);
    const auto *hl_792 = buffer.data(hl + 792);
    const auto *hl_793 = buffer.data(hl + 793);
    const auto *hl_794 = buffer.data(hl + 794);
    const auto *hl_795 = buffer.data(hl + 795);
    const auto *hl_796 = buffer.data(hl + 796);
    const auto *hl_797 = buffer.data(hl + 797);
    const auto *hl_798 = buffer.data(hl + 798);
    const auto *hl_799 = buffer.data(hl + 799);
    const auto *hl_800 = buffer.data(hl + 800);
    const auto *hl_801 = buffer.data(hl + 801);
    const auto *hl_802 = buffer.data(hl + 802);
    const auto *hl_803 = buffer.data(hl + 803);
    const auto *hl_804 = buffer.data(hl + 804);
    const auto *hl_805 = buffer.data(hl + 805);
    const auto *hl_806 = buffer.data(hl + 806);
    const auto *hl_807 = buffer.data(hl + 807);
    const auto *hl_808 = buffer.data(hl + 808);
    const auto *hl_809 = buffer.data(hl + 809);
    const auto *hl_810 = buffer.data(hl + 810);
    const auto *hl_811 = buffer.data(hl + 811);
    const auto *hl_812 = buffer.data(hl + 812);
    const auto *hl_813 = buffer.data(hl + 813);
    const auto *hl_814 = buffer.data(hl + 814);
    const auto *hl_815 = buffer.data(hl + 815);
    const auto *hl_816 = buffer.data(hl + 816);
    const auto *hl_817 = buffer.data(hl + 817);
    const auto *hl_818 = buffer.data(hl + 818);
    const auto *hl_819 = buffer.data(hl + 819);
    const auto *hl_820 = buffer.data(hl + 820);
    const auto *hl_821 = buffer.data(hl + 821);
    const auto *hl_822 = buffer.data(hl + 822);
    const auto *hl_823 = buffer.data(hl + 823);
    const auto *hl_824 = buffer.data(hl + 824);
    const auto *hl_825 = buffer.data(hl + 825);
    const auto *hl_826 = buffer.data(hl + 826);
    const auto *hl_827 = buffer.data(hl + 827);
    const auto *hl_828 = buffer.data(hl + 828);
    const auto *hl_829 = buffer.data(hl + 829);
    const auto *hl_830 = buffer.data(hl + 830);
    const auto *hl_831 = buffer.data(hl + 831);
    const auto *hl_832 = buffer.data(hl + 832);
    const auto *hl_833 = buffer.data(hl + 833);
    const auto *hl_834 = buffer.data(hl + 834);
    const auto *hl_835 = buffer.data(hl + 835);
    const auto *hl_836 = buffer.data(hl + 836);
    const auto *hl_837 = buffer.data(hl + 837);
    const auto *hl_838 = buffer.data(hl + 838);
    const auto *hl_839 = buffer.data(hl + 839);
    const auto *hl_840 = buffer.data(hl + 840);
    const auto *hl_841 = buffer.data(hl + 841);
    const auto *hl_842 = buffer.data(hl + 842);
    const auto *hl_843 = buffer.data(hl + 843);
    const auto *hl_844 = buffer.data(hl + 844);
    const auto *hl_845 = buffer.data(hl + 845);
    const auto *hl_846 = buffer.data(hl + 846);
    const auto *hl_847 = buffer.data(hl + 847);
    const auto *hl_848 = buffer.data(hl + 848);
    const auto *hl_849 = buffer.data(hl + 849);
    const auto *hl_850 = buffer.data(hl + 850);
    const auto *hl_851 = buffer.data(hl + 851);
    const auto *hl_852 = buffer.data(hl + 852);
    const auto *hl_853 = buffer.data(hl + 853);
    const auto *hl_854 = buffer.data(hl + 854);
    const auto *hl_855 = buffer.data(hl + 855);
    const auto *hl_856 = buffer.data(hl + 856);
    const auto *hl_857 = buffer.data(hl + 857);
    const auto *hl_858 = buffer.data(hl + 858);
    const auto *hl_859 = buffer.data(hl + 859);
    const auto *hl_860 = buffer.data(hl + 860);
    const auto *hl_861 = buffer.data(hl + 861);
    const auto *hl_862 = buffer.data(hl + 862);
    const auto *hl_863 = buffer.data(hl + 863);
    const auto *hl_864 = buffer.data(hl + 864);
    const auto *hl_865 = buffer.data(hl + 865);
    const auto *hl_866 = buffer.data(hl + 866);
    const auto *hl_867 = buffer.data(hl + 867);
    const auto *hl_868 = buffer.data(hl + 868);
    const auto *hl_869 = buffer.data(hl + 869);
    const auto *hl_870 = buffer.data(hl + 870);
    const auto *hl_871 = buffer.data(hl + 871);
    const auto *hl_872 = buffer.data(hl + 872);
    const auto *hl_873 = buffer.data(hl + 873);
    const auto *hl_874 = buffer.data(hl + 874);
    const auto *hl_875 = buffer.data(hl + 875);
    const auto *hl_876 = buffer.data(hl + 876);
    const auto *hl_877 = buffer.data(hl + 877);
    const auto *hl_878 = buffer.data(hl + 878);
    const auto *hl_879 = buffer.data(hl + 879);
    const auto *hl_880 = buffer.data(hl + 880);
    const auto *hl_881 = buffer.data(hl + 881);
    const auto *hl_882 = buffer.data(hl + 882);
    const auto *hl_883 = buffer.data(hl + 883);
    const auto *hl_884 = buffer.data(hl + 884);
    const auto *hl_885 = buffer.data(hl + 885);
    const auto *hl_886 = buffer.data(hl + 886);
    const auto *hl_887 = buffer.data(hl + 887);
    const auto *hl_888 = buffer.data(hl + 888);
    const auto *hl_889 = buffer.data(hl + 889);
    const auto *hl_890 = buffer.data(hl + 890);
    const auto *hl_891 = buffer.data(hl + 891);
    const auto *hl_892 = buffer.data(hl + 892);
    const auto *hl_893 = buffer.data(hl + 893);
    const auto *hl_894 = buffer.data(hl + 894);
    const auto *hl_895 = buffer.data(hl + 895);
    const auto *hl_896 = buffer.data(hl + 896);
    const auto *hl_897 = buffer.data(hl + 897);
    const auto *hl_898 = buffer.data(hl + 898);
    const auto *hl_899 = buffer.data(hl + 899);
    const auto *hl_900 = buffer.data(hl + 900);
    const auto *hl_901 = buffer.data(hl + 901);
    const auto *hl_902 = buffer.data(hl + 902);
    const auto *hl_903 = buffer.data(hl + 903);
    const auto *hl_904 = buffer.data(hl + 904);
    const auto *hl_905 = buffer.data(hl + 905);
    const auto *hl_906 = buffer.data(hl + 906);
    const auto *hl_907 = buffer.data(hl + 907);
    const auto *hl_908 = buffer.data(hl + 908);
    const auto *hl_909 = buffer.data(hl + 909);
    const auto *hl_910 = buffer.data(hl + 910);
    const auto *hl_911 = buffer.data(hl + 911);
    const auto *hl_912 = buffer.data(hl + 912);
    const auto *hl_913 = buffer.data(hl + 913);
    const auto *hl_914 = buffer.data(hl + 914);
    const auto *hl_915 = buffer.data(hl + 915);
    const auto *hl_916 = buffer.data(hl + 916);
    const auto *hl_917 = buffer.data(hl + 917);
    const auto *hl_918 = buffer.data(hl + 918);
    const auto *hl_919 = buffer.data(hl + 919);
    const auto *hl_920 = buffer.data(hl + 920);
    const auto *hl_921 = buffer.data(hl + 921);
    const auto *hl_922 = buffer.data(hl + 922);
    const auto *hl_923 = buffer.data(hl + 923);
    const auto *hl_924 = buffer.data(hl + 924);
    const auto *hl_925 = buffer.data(hl + 925);
    const auto *hl_926 = buffer.data(hl + 926);
    const auto *hl_927 = buffer.data(hl + 927);
    const auto *hl_928 = buffer.data(hl + 928);
    const auto *hl_929 = buffer.data(hl + 929);
    const auto *hl_930 = buffer.data(hl + 930);
    const auto *hl_931 = buffer.data(hl + 931);
    const auto *hl_932 = buffer.data(hl + 932);
    const auto *hl_933 = buffer.data(hl + 933);
    const auto *hl_934 = buffer.data(hl + 934);
    const auto *hl_935 = buffer.data(hl + 935);
    const auto *hl_936 = buffer.data(hl + 936);
    const auto *hl_937 = buffer.data(hl + 937);
    const auto *hl_938 = buffer.data(hl + 938);

#pragma omp simd aligned(t_519, t_520, t_521, t_522, t_523, fl_294, fl_295, fl_296, fl_297, \
                         fl_298, hl_789, hl_790, hl_791, hl_792, \
                         hl_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = -fl_294[k]
                   + f_0 * hl_789[k];

        t_520[k] = -fl_295[k]
                   + f_0 * hl_790[k];

        t_521[k] = -fl_296[k]
                   + f_0 * hl_791[k];

        t_522[k] = -fl_297[k]
                   + f_0 * hl_792[k];

        t_523[k] = -fl_298[k]
                   + f_0 * hl_793[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, fl_299, fl_300, fl_301, fl_302, \
                         fl_303, hl_794, hl_795, hl_796, hl_797, \
                         hl_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = -fl_299[k]
                   + f_0 * hl_794[k];

        t_525[k] = -fl_300[k]
                   + f_0 * hl_795[k];

        t_526[k] = -fl_301[k]
                   + f_0 * hl_796[k];

        t_527[k] = -fl_302[k]
                   + f_0 * hl_797[k];

        t_528[k] = -fl_303[k]
                   + f_0 * hl_798[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, t_533, fl_304, fl_305, fl_306, fl_307, \
                         fl_308, hl_799, hl_800, hl_801, hl_802, \
                         hl_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = -fl_304[k]
                   + f_0 * hl_799[k];

        t_530[k] = -fl_305[k]
                   + f_0 * hl_800[k];

        t_531[k] = -fl_306[k]
                   + f_0 * hl_801[k];

        t_532[k] = -fl_307[k]
                   + f_0 * hl_802[k];

        t_533[k] = -fl_308[k]
                   + f_0 * hl_803[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, t_538, fl_309, fl_310, fl_311, fl_312, \
                         fl_313, hl_804, hl_805, hl_806, hl_807, \
                         hl_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = -fl_309[k]
                   + f_0 * hl_804[k];

        t_535[k] = -fl_310[k]
                   + f_0 * hl_805[k];

        t_536[k] = -fl_311[k]
                   + f_0 * hl_806[k];

        t_537[k] = -fl_312[k]
                   + f_0 * hl_807[k];

        t_538[k] = -fl_313[k]
                   + f_0 * hl_808[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, t_542, t_543, fl_314, fl_315, fl_316, fl_317, \
                         fl_318, hl_809, hl_810, hl_811, hl_812, \
                         hl_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = -fl_314[k]
                   + f_0 * hl_809[k];

        t_540[k] = -2.0 * fl_315[k]
                   + f_0 * hl_810[k];

        t_541[k] = -2.0 * fl_316[k]
                   + f_0 * hl_811[k];

        t_542[k] = -2.0 * fl_317[k]
                   + f_0 * hl_812[k];

        t_543[k] = -2.0 * fl_318[k]
                   + f_0 * hl_813[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, t_548, fl_319, fl_320, fl_321, fl_322, \
                         fl_323, hl_814, hl_815, hl_816, hl_817, \
                         hl_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = -2.0 * fl_319[k]
                   + f_0 * hl_814[k];

        t_545[k] = -2.0 * fl_320[k]
                   + f_0 * hl_815[k];

        t_546[k] = -2.0 * fl_321[k]
                   + f_0 * hl_816[k];

        t_547[k] = -2.0 * fl_322[k]
                   + f_0 * hl_817[k];

        t_548[k] = -2.0 * fl_323[k]
                   + f_0 * hl_818[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, t_553, fl_324, fl_325, fl_326, fl_327, \
                         fl_328, hl_819, hl_820, hl_821, hl_822, \
                         hl_823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = -2.0 * fl_324[k]
                   + f_0 * hl_819[k];

        t_550[k] = -2.0 * fl_325[k]
                   + f_0 * hl_820[k];

        t_551[k] = -2.0 * fl_326[k]
                   + f_0 * hl_821[k];

        t_552[k] = -2.0 * fl_327[k]
                   + f_0 * hl_822[k];

        t_553[k] = -2.0 * fl_328[k]
                   + f_0 * hl_823[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, t_558, fl_329, fl_330, fl_331, fl_332, \
                         fl_333, hl_824, hl_825, hl_826, hl_827, \
                         hl_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = -2.0 * fl_329[k]
                   + f_0 * hl_824[k];

        t_555[k] = -2.0 * fl_330[k]
                   + f_0 * hl_825[k];

        t_556[k] = -2.0 * fl_331[k]
                   + f_0 * hl_826[k];

        t_557[k] = -2.0 * fl_332[k]
                   + f_0 * hl_827[k];

        t_558[k] = -2.0 * fl_333[k]
                   + f_0 * hl_828[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, t_563, fl_334, fl_335, fl_336, fl_337, \
                         fl_338, hl_829, hl_830, hl_831, hl_832, \
                         hl_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = -2.0 * fl_334[k]
                   + f_0 * hl_829[k];

        t_560[k] = -2.0 * fl_335[k]
                   + f_0 * hl_830[k];

        t_561[k] = -2.0 * fl_336[k]
                   + f_0 * hl_831[k];

        t_562[k] = -2.0 * fl_337[k]
                   + f_0 * hl_832[k];

        t_563[k] = -2.0 * fl_338[k]
                   + f_0 * hl_833[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, t_568, fl_339, fl_340, fl_341, fl_342, \
                         fl_343, hl_834, hl_835, hl_836, hl_837, \
                         hl_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = -2.0 * fl_339[k]
                   + f_0 * hl_834[k];

        t_565[k] = -2.0 * fl_340[k]
                   + f_0 * hl_835[k];

        t_566[k] = -2.0 * fl_341[k]
                   + f_0 * hl_836[k];

        t_567[k] = -2.0 * fl_342[k]
                   + f_0 * hl_837[k];

        t_568[k] = -2.0 * fl_343[k]
                   + f_0 * hl_838[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, t_573, fl_344, fl_345, fl_346, fl_347, \
                         fl_348, hl_839, hl_840, hl_841, hl_842, \
                         hl_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = -2.0 * fl_344[k]
                   + f_0 * hl_839[k];

        t_570[k] = -2.0 * fl_345[k]
                   + f_0 * hl_840[k];

        t_571[k] = -2.0 * fl_346[k]
                   + f_0 * hl_841[k];

        t_572[k] = -2.0 * fl_347[k]
                   + f_0 * hl_842[k];

        t_573[k] = -2.0 * fl_348[k]
                   + f_0 * hl_843[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, fl_349, fl_350, fl_351, fl_352, \
                         fl_353, hl_844, hl_845, hl_846, hl_847, \
                         hl_848 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = -2.0 * fl_349[k]
                   + f_0 * hl_844[k];

        t_575[k] = -2.0 * fl_350[k]
                   + f_0 * hl_845[k];

        t_576[k] = -2.0 * fl_351[k]
                   + f_0 * hl_846[k];

        t_577[k] = -2.0 * fl_352[k]
                   + f_0 * hl_847[k];

        t_578[k] = -2.0 * fl_353[k]
                   + f_0 * hl_848[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, t_582, t_583, fl_354, fl_355, fl_356, fl_357, \
                         fl_358, hl_849, hl_850, hl_851, hl_852, \
                         hl_853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = -2.0 * fl_354[k]
                   + f_0 * hl_849[k];

        t_580[k] = -2.0 * fl_355[k]
                   + f_0 * hl_850[k];

        t_581[k] = -2.0 * fl_356[k]
                   + f_0 * hl_851[k];

        t_582[k] = -2.0 * fl_357[k]
                   + f_0 * hl_852[k];

        t_583[k] = -2.0 * fl_358[k]
                   + f_0 * hl_853[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, t_588, fl_359, fl_360, fl_361, fl_362, \
                         fl_363, hl_854, hl_855, hl_856, hl_857, \
                         hl_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = -2.0 * fl_359[k]
                   + f_0 * hl_854[k];

        t_585[k] = -3.0 * fl_360[k]
                   + f_0 * hl_855[k];

        t_586[k] = -3.0 * fl_361[k]
                   + f_0 * hl_856[k];

        t_587[k] = -3.0 * fl_362[k]
                   + f_0 * hl_857[k];

        t_588[k] = -3.0 * fl_363[k]
                   + f_0 * hl_858[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, t_593, fl_364, fl_365, fl_366, fl_367, \
                         fl_368, hl_859, hl_860, hl_861, hl_862, \
                         hl_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = -3.0 * fl_364[k]
                   + f_0 * hl_859[k];

        t_590[k] = -3.0 * fl_365[k]
                   + f_0 * hl_860[k];

        t_591[k] = -3.0 * fl_366[k]
                   + f_0 * hl_861[k];

        t_592[k] = -3.0 * fl_367[k]
                   + f_0 * hl_862[k];

        t_593[k] = -3.0 * fl_368[k]
                   + f_0 * hl_863[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, t_598, fl_369, fl_370, fl_371, fl_372, \
                         fl_373, hl_864, hl_865, hl_866, hl_867, \
                         hl_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = -3.0 * fl_369[k]
                   + f_0 * hl_864[k];

        t_595[k] = -3.0 * fl_370[k]
                   + f_0 * hl_865[k];

        t_596[k] = -3.0 * fl_371[k]
                   + f_0 * hl_866[k];

        t_597[k] = -3.0 * fl_372[k]
                   + f_0 * hl_867[k];

        t_598[k] = -3.0 * fl_373[k]
                   + f_0 * hl_868[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, fl_374, fl_375, fl_376, fl_377, \
                         fl_378, hl_869, hl_870, hl_871, hl_872, \
                         hl_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = -3.0 * fl_374[k]
                   + f_0 * hl_869[k];

        t_600[k] = -3.0 * fl_375[k]
                   + f_0 * hl_870[k];

        t_601[k] = -3.0 * fl_376[k]
                   + f_0 * hl_871[k];

        t_602[k] = -3.0 * fl_377[k]
                   + f_0 * hl_872[k];

        t_603[k] = -3.0 * fl_378[k]
                   + f_0 * hl_873[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, fl_379, fl_380, fl_381, fl_382, \
                         fl_383, hl_874, hl_875, hl_876, hl_877, \
                         hl_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = -3.0 * fl_379[k]
                   + f_0 * hl_874[k];

        t_605[k] = -3.0 * fl_380[k]
                   + f_0 * hl_875[k];

        t_606[k] = -3.0 * fl_381[k]
                   + f_0 * hl_876[k];

        t_607[k] = -3.0 * fl_382[k]
                   + f_0 * hl_877[k];

        t_608[k] = -3.0 * fl_383[k]
                   + f_0 * hl_878[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, t_613, fl_384, fl_385, fl_386, fl_387, \
                         fl_388, hl_879, hl_880, hl_881, hl_882, \
                         hl_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = -3.0 * fl_384[k]
                   + f_0 * hl_879[k];

        t_610[k] = -3.0 * fl_385[k]
                   + f_0 * hl_880[k];

        t_611[k] = -3.0 * fl_386[k]
                   + f_0 * hl_881[k];

        t_612[k] = -3.0 * fl_387[k]
                   + f_0 * hl_882[k];

        t_613[k] = -3.0 * fl_388[k]
                   + f_0 * hl_883[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, fl_389, fl_390, fl_391, fl_392, \
                         fl_393, hl_884, hl_885, hl_886, hl_887, \
                         hl_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = -3.0 * fl_389[k]
                   + f_0 * hl_884[k];

        t_615[k] = -3.0 * fl_390[k]
                   + f_0 * hl_885[k];

        t_616[k] = -3.0 * fl_391[k]
                   + f_0 * hl_886[k];

        t_617[k] = -3.0 * fl_392[k]
                   + f_0 * hl_887[k];

        t_618[k] = -3.0 * fl_393[k]
                   + f_0 * hl_888[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, fl_394, fl_395, fl_396, fl_397, \
                         fl_398, hl_889, hl_890, hl_891, hl_892, \
                         hl_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = -3.0 * fl_394[k]
                   + f_0 * hl_889[k];

        t_620[k] = -3.0 * fl_395[k]
                   + f_0 * hl_890[k];

        t_621[k] = -3.0 * fl_396[k]
                   + f_0 * hl_891[k];

        t_622[k] = -3.0 * fl_397[k]
                   + f_0 * hl_892[k];

        t_623[k] = -3.0 * fl_398[k]
                   + f_0 * hl_893[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, t_628, fl_399, fl_400, fl_401, fl_402, \
                         fl_403, hl_894, hl_895, hl_896, hl_897, \
                         hl_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = -3.0 * fl_399[k]
                   + f_0 * hl_894[k];

        t_625[k] = -3.0 * fl_400[k]
                   + f_0 * hl_895[k];

        t_626[k] = -3.0 * fl_401[k]
                   + f_0 * hl_896[k];

        t_627[k] = -3.0 * fl_402[k]
                   + f_0 * hl_897[k];

        t_628[k] = -3.0 * fl_403[k]
                   + f_0 * hl_898[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, t_633, fl_404, fl_405, fl_406, fl_407, \
                         fl_408, hl_899, hl_900, hl_901, hl_902, \
                         hl_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = -3.0 * fl_404[k]
                   + f_0 * hl_899[k];

        t_630[k] = -4.0 * fl_405[k]
                   + f_0 * hl_900[k];

        t_631[k] = -4.0 * fl_406[k]
                   + f_0 * hl_901[k];

        t_632[k] = -4.0 * fl_407[k]
                   + f_0 * hl_902[k];

        t_633[k] = -4.0 * fl_408[k]
                   + f_0 * hl_903[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, fl_409, fl_410, fl_411, fl_412, \
                         fl_413, hl_904, hl_905, hl_906, hl_907, \
                         hl_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = -4.0 * fl_409[k]
                   + f_0 * hl_904[k];

        t_635[k] = -4.0 * fl_410[k]
                   + f_0 * hl_905[k];

        t_636[k] = -4.0 * fl_411[k]
                   + f_0 * hl_906[k];

        t_637[k] = -4.0 * fl_412[k]
                   + f_0 * hl_907[k];

        t_638[k] = -4.0 * fl_413[k]
                   + f_0 * hl_908[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, t_643, fl_414, fl_415, fl_416, fl_417, \
                         fl_418, hl_909, hl_910, hl_911, hl_912, \
                         hl_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = -4.0 * fl_414[k]
                   + f_0 * hl_909[k];

        t_640[k] = -4.0 * fl_415[k]
                   + f_0 * hl_910[k];

        t_641[k] = -4.0 * fl_416[k]
                   + f_0 * hl_911[k];

        t_642[k] = -4.0 * fl_417[k]
                   + f_0 * hl_912[k];

        t_643[k] = -4.0 * fl_418[k]
                   + f_0 * hl_913[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, t_648, fl_419, fl_420, fl_421, fl_422, \
                         fl_423, hl_914, hl_915, hl_916, hl_917, \
                         hl_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = -4.0 * fl_419[k]
                   + f_0 * hl_914[k];

        t_645[k] = -4.0 * fl_420[k]
                   + f_0 * hl_915[k];

        t_646[k] = -4.0 * fl_421[k]
                   + f_0 * hl_916[k];

        t_647[k] = -4.0 * fl_422[k]
                   + f_0 * hl_917[k];

        t_648[k] = -4.0 * fl_423[k]
                   + f_0 * hl_918[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, fl_424, fl_425, fl_426, fl_427, \
                         fl_428, hl_919, hl_920, hl_921, hl_922, \
                         hl_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = -4.0 * fl_424[k]
                   + f_0 * hl_919[k];

        t_650[k] = -4.0 * fl_425[k]
                   + f_0 * hl_920[k];

        t_651[k] = -4.0 * fl_426[k]
                   + f_0 * hl_921[k];

        t_652[k] = -4.0 * fl_427[k]
                   + f_0 * hl_922[k];

        t_653[k] = -4.0 * fl_428[k]
                   + f_0 * hl_923[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, t_657, t_658, fl_429, fl_430, fl_431, fl_432, \
                         fl_433, hl_924, hl_925, hl_926, hl_927, \
                         hl_928 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = -4.0 * fl_429[k]
                   + f_0 * hl_924[k];

        t_655[k] = -4.0 * fl_430[k]
                   + f_0 * hl_925[k];

        t_656[k] = -4.0 * fl_431[k]
                   + f_0 * hl_926[k];

        t_657[k] = -4.0 * fl_432[k]
                   + f_0 * hl_927[k];

        t_658[k] = -4.0 * fl_433[k]
                   + f_0 * hl_928[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, fl_434, fl_435, fl_436, fl_437, \
                         fl_438, hl_929, hl_930, hl_931, hl_932, \
                         hl_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = -4.0 * fl_434[k]
                   + f_0 * hl_929[k];

        t_660[k] = -4.0 * fl_435[k]
                   + f_0 * hl_930[k];

        t_661[k] = -4.0 * fl_436[k]
                   + f_0 * hl_931[k];

        t_662[k] = -4.0 * fl_437[k]
                   + f_0 * hl_932[k];

        t_663[k] = -4.0 * fl_438[k]
                   + f_0 * hl_933[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, t_668, fl_439, fl_440, fl_441, fl_442, \
                         fl_443, hl_934, hl_935, hl_936, hl_937, \
                         hl_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = -4.0 * fl_439[k]
                   + f_0 * hl_934[k];

        t_665[k] = -4.0 * fl_440[k]
                   + f_0 * hl_935[k];

        t_666[k] = -4.0 * fl_441[k]
                   + f_0 * hl_936[k];

        t_667[k] = -4.0 * fl_442[k]
                   + f_0 * hl_937[k];

        t_668[k] = -4.0 * fl_443[k]
                   + f_0 * hl_938[k];
    }
}

static auto
compute_prim_geom_10_gl_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t hl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);

    const auto *fl_444 = buffer.data(fl + 444);
    const auto *fl_445 = buffer.data(fl + 445);
    const auto *fl_446 = buffer.data(fl + 446);
    const auto *fl_447 = buffer.data(fl + 447);
    const auto *fl_448 = buffer.data(fl + 448);
    const auto *fl_449 = buffer.data(fl + 449);

    const auto *hl_939 = buffer.data(hl + 939);
    const auto *hl_940 = buffer.data(hl + 940);
    const auto *hl_941 = buffer.data(hl + 941);
    const auto *hl_942 = buffer.data(hl + 942);
    const auto *hl_943 = buffer.data(hl + 943);
    const auto *hl_944 = buffer.data(hl + 944);

#pragma omp simd aligned(t_669, t_670, t_671, t_672, t_673, fl_444, fl_445, fl_446, fl_447, \
                         fl_448, hl_939, hl_940, hl_941, hl_942, \
                         hl_943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = -4.0 * fl_444[k]
                   + f_0 * hl_939[k];

        t_670[k] = -4.0 * fl_445[k]
                   + f_0 * hl_940[k];

        t_671[k] = -4.0 * fl_446[k]
                   + f_0 * hl_941[k];

        t_672[k] = -4.0 * fl_447[k]
                   + f_0 * hl_942[k];

        t_673[k] = -4.0 * fl_448[k]
                   + f_0 * hl_943[k];
    }

#pragma omp simd aligned(t_674, fl_449, hl_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = -4.0 * fl_449[k]
                   + f_0 * hl_944[k];
    }
}

auto
compute_prim_geom_10_gl_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t fl, const size_t hl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gl_electron_repulsion_2_piece0(buffer, target, fl, hl, ncols, alpha);

    compute_prim_geom_10_gl_electron_repulsion_2_piece1(buffer, target, fl, hl, ncols, alpha);

    compute_prim_geom_10_gl_electron_repulsion_2_piece2(buffer, target, fl, hl, ncols, alpha);

    compute_prim_geom_10_gl_electron_repulsion_2_piece3(buffer, target, fl, hl, ncols, alpha);

    compute_prim_geom_10_gl_electron_repulsion_2_piece4(buffer, target, fl, hl, ncols, alpha);
}

}  // namespace simdt2ceri
