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


#include "SimdElectronRepulsionGeom10VrrRecHL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_hl_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);
    const auto *gl_30 = buffer.data(gl + 30);
    const auto *gl_31 = buffer.data(gl + 31);
    const auto *gl_32 = buffer.data(gl + 32);
    const auto *gl_33 = buffer.data(gl + 33);
    const auto *gl_34 = buffer.data(gl + 34);
    const auto *gl_35 = buffer.data(gl + 35);
    const auto *gl_36 = buffer.data(gl + 36);
    const auto *gl_37 = buffer.data(gl + 37);
    const auto *gl_38 = buffer.data(gl + 38);
    const auto *gl_39 = buffer.data(gl + 39);
    const auto *gl_40 = buffer.data(gl + 40);
    const auto *gl_41 = buffer.data(gl + 41);
    const auto *gl_42 = buffer.data(gl + 42);
    const auto *gl_43 = buffer.data(gl + 43);
    const auto *gl_44 = buffer.data(gl + 44);
    const auto *gl_45 = buffer.data(gl + 45);
    const auto *gl_46 = buffer.data(gl + 46);
    const auto *gl_47 = buffer.data(gl + 47);
    const auto *gl_48 = buffer.data(gl + 48);
    const auto *gl_49 = buffer.data(gl + 49);
    const auto *gl_50 = buffer.data(gl + 50);
    const auto *gl_51 = buffer.data(gl + 51);
    const auto *gl_52 = buffer.data(gl + 52);
    const auto *gl_53 = buffer.data(gl + 53);
    const auto *gl_54 = buffer.data(gl + 54);
    const auto *gl_55 = buffer.data(gl + 55);
    const auto *gl_56 = buffer.data(gl + 56);
    const auto *gl_57 = buffer.data(gl + 57);
    const auto *gl_58 = buffer.data(gl + 58);
    const auto *gl_59 = buffer.data(gl + 59);
    const auto *gl_60 = buffer.data(gl + 60);
    const auto *gl_61 = buffer.data(gl + 61);
    const auto *gl_62 = buffer.data(gl + 62);
    const auto *gl_63 = buffer.data(gl + 63);
    const auto *gl_64 = buffer.data(gl + 64);
    const auto *gl_65 = buffer.data(gl + 65);
    const auto *gl_66 = buffer.data(gl + 66);
    const auto *gl_67 = buffer.data(gl + 67);
    const auto *gl_68 = buffer.data(gl + 68);
    const auto *gl_69 = buffer.data(gl + 69);
    const auto *gl_70 = buffer.data(gl + 70);
    const auto *gl_71 = buffer.data(gl + 71);
    const auto *gl_72 = buffer.data(gl + 72);
    const auto *gl_73 = buffer.data(gl + 73);
    const auto *gl_74 = buffer.data(gl + 74);
    const auto *gl_75 = buffer.data(gl + 75);
    const auto *gl_76 = buffer.data(gl + 76);
    const auto *gl_77 = buffer.data(gl + 77);
    const auto *gl_78 = buffer.data(gl + 78);
    const auto *gl_79 = buffer.data(gl + 79);
    const auto *gl_80 = buffer.data(gl + 80);
    const auto *gl_81 = buffer.data(gl + 81);
    const auto *gl_82 = buffer.data(gl + 82);
    const auto *gl_83 = buffer.data(gl + 83);
    const auto *gl_84 = buffer.data(gl + 84);
    const auto *gl_85 = buffer.data(gl + 85);
    const auto *gl_86 = buffer.data(gl + 86);
    const auto *gl_87 = buffer.data(gl + 87);
    const auto *gl_88 = buffer.data(gl + 88);
    const auto *gl_89 = buffer.data(gl + 89);
    const auto *gl_90 = buffer.data(gl + 90);
    const auto *gl_91 = buffer.data(gl + 91);
    const auto *gl_92 = buffer.data(gl + 92);
    const auto *gl_93 = buffer.data(gl + 93);
    const auto *gl_94 = buffer.data(gl + 94);
    const auto *gl_95 = buffer.data(gl + 95);
    const auto *gl_96 = buffer.data(gl + 96);
    const auto *gl_97 = buffer.data(gl + 97);
    const auto *gl_98 = buffer.data(gl + 98);
    const auto *gl_99 = buffer.data(gl + 99);
    const auto *gl_100 = buffer.data(gl + 100);
    const auto *gl_101 = buffer.data(gl + 101);
    const auto *gl_102 = buffer.data(gl + 102);
    const auto *gl_103 = buffer.data(gl + 103);
    const auto *gl_104 = buffer.data(gl + 104);
    const auto *gl_105 = buffer.data(gl + 105);
    const auto *gl_106 = buffer.data(gl + 106);
    const auto *gl_107 = buffer.data(gl + 107);
    const auto *gl_108 = buffer.data(gl + 108);
    const auto *gl_109 = buffer.data(gl + 109);
    const auto *gl_110 = buffer.data(gl + 110);
    const auto *gl_111 = buffer.data(gl + 111);
    const auto *gl_112 = buffer.data(gl + 112);
    const auto *gl_113 = buffer.data(gl + 113);
    const auto *gl_114 = buffer.data(gl + 114);
    const auto *gl_115 = buffer.data(gl + 115);
    const auto *gl_116 = buffer.data(gl + 116);
    const auto *gl_117 = buffer.data(gl + 117);
    const auto *gl_118 = buffer.data(gl + 118);
    const auto *gl_119 = buffer.data(gl + 119);
    const auto *gl_120 = buffer.data(gl + 120);
    const auto *gl_121 = buffer.data(gl + 121);
    const auto *gl_122 = buffer.data(gl + 122);
    const auto *gl_123 = buffer.data(gl + 123);
    const auto *gl_124 = buffer.data(gl + 124);
    const auto *gl_125 = buffer.data(gl + 125);
    const auto *gl_126 = buffer.data(gl + 126);
    const auto *gl_127 = buffer.data(gl + 127);
    const auto *gl_128 = buffer.data(gl + 128);
    const auto *gl_129 = buffer.data(gl + 129);
    const auto *gl_130 = buffer.data(gl + 130);
    const auto *gl_131 = buffer.data(gl + 131);
    const auto *gl_132 = buffer.data(gl + 132);
    const auto *gl_133 = buffer.data(gl + 133);
    const auto *gl_134 = buffer.data(gl + 134);
    const auto *gl_135 = buffer.data(gl + 135);
    const auto *gl_136 = buffer.data(gl + 136);
    const auto *gl_137 = buffer.data(gl + 137);
    const auto *gl_138 = buffer.data(gl + 138);
    const auto *gl_139 = buffer.data(gl + 139);
    const auto *gl_140 = buffer.data(gl + 140);
    const auto *gl_141 = buffer.data(gl + 141);
    const auto *gl_142 = buffer.data(gl + 142);
    const auto *gl_143 = buffer.data(gl + 143);
    const auto *gl_144 = buffer.data(gl + 144);
    const auto *gl_145 = buffer.data(gl + 145);
    const auto *gl_146 = buffer.data(gl + 146);
    const auto *gl_147 = buffer.data(gl + 147);
    const auto *gl_148 = buffer.data(gl + 148);
    const auto *gl_149 = buffer.data(gl + 149);

    const auto *il_0 = buffer.data(il + 0);
    const auto *il_1 = buffer.data(il + 1);
    const auto *il_2 = buffer.data(il + 2);
    const auto *il_3 = buffer.data(il + 3);
    const auto *il_4 = buffer.data(il + 4);
    const auto *il_5 = buffer.data(il + 5);
    const auto *il_6 = buffer.data(il + 6);
    const auto *il_7 = buffer.data(il + 7);
    const auto *il_8 = buffer.data(il + 8);
    const auto *il_9 = buffer.data(il + 9);
    const auto *il_10 = buffer.data(il + 10);
    const auto *il_11 = buffer.data(il + 11);
    const auto *il_12 = buffer.data(il + 12);
    const auto *il_13 = buffer.data(il + 13);
    const auto *il_14 = buffer.data(il + 14);
    const auto *il_15 = buffer.data(il + 15);
    const auto *il_16 = buffer.data(il + 16);
    const auto *il_17 = buffer.data(il + 17);
    const auto *il_18 = buffer.data(il + 18);
    const auto *il_19 = buffer.data(il + 19);
    const auto *il_20 = buffer.data(il + 20);
    const auto *il_21 = buffer.data(il + 21);
    const auto *il_22 = buffer.data(il + 22);
    const auto *il_23 = buffer.data(il + 23);
    const auto *il_24 = buffer.data(il + 24);
    const auto *il_25 = buffer.data(il + 25);
    const auto *il_26 = buffer.data(il + 26);
    const auto *il_27 = buffer.data(il + 27);
    const auto *il_28 = buffer.data(il + 28);
    const auto *il_29 = buffer.data(il + 29);
    const auto *il_30 = buffer.data(il + 30);
    const auto *il_31 = buffer.data(il + 31);
    const auto *il_32 = buffer.data(il + 32);
    const auto *il_33 = buffer.data(il + 33);
    const auto *il_34 = buffer.data(il + 34);
    const auto *il_35 = buffer.data(il + 35);
    const auto *il_36 = buffer.data(il + 36);
    const auto *il_37 = buffer.data(il + 37);
    const auto *il_38 = buffer.data(il + 38);
    const auto *il_39 = buffer.data(il + 39);
    const auto *il_40 = buffer.data(il + 40);
    const auto *il_41 = buffer.data(il + 41);
    const auto *il_42 = buffer.data(il + 42);
    const auto *il_43 = buffer.data(il + 43);
    const auto *il_44 = buffer.data(il + 44);
    const auto *il_45 = buffer.data(il + 45);
    const auto *il_46 = buffer.data(il + 46);
    const auto *il_47 = buffer.data(il + 47);
    const auto *il_48 = buffer.data(il + 48);
    const auto *il_49 = buffer.data(il + 49);
    const auto *il_50 = buffer.data(il + 50);
    const auto *il_51 = buffer.data(il + 51);
    const auto *il_52 = buffer.data(il + 52);
    const auto *il_53 = buffer.data(il + 53);
    const auto *il_54 = buffer.data(il + 54);
    const auto *il_55 = buffer.data(il + 55);
    const auto *il_56 = buffer.data(il + 56);
    const auto *il_57 = buffer.data(il + 57);
    const auto *il_58 = buffer.data(il + 58);
    const auto *il_59 = buffer.data(il + 59);
    const auto *il_60 = buffer.data(il + 60);
    const auto *il_61 = buffer.data(il + 61);
    const auto *il_62 = buffer.data(il + 62);
    const auto *il_63 = buffer.data(il + 63);
    const auto *il_64 = buffer.data(il + 64);
    const auto *il_65 = buffer.data(il + 65);
    const auto *il_66 = buffer.data(il + 66);
    const auto *il_67 = buffer.data(il + 67);
    const auto *il_68 = buffer.data(il + 68);
    const auto *il_69 = buffer.data(il + 69);
    const auto *il_70 = buffer.data(il + 70);
    const auto *il_71 = buffer.data(il + 71);
    const auto *il_72 = buffer.data(il + 72);
    const auto *il_73 = buffer.data(il + 73);
    const auto *il_74 = buffer.data(il + 74);
    const auto *il_75 = buffer.data(il + 75);
    const auto *il_76 = buffer.data(il + 76);
    const auto *il_77 = buffer.data(il + 77);
    const auto *il_78 = buffer.data(il + 78);
    const auto *il_79 = buffer.data(il + 79);
    const auto *il_80 = buffer.data(il + 80);
    const auto *il_81 = buffer.data(il + 81);
    const auto *il_82 = buffer.data(il + 82);
    const auto *il_83 = buffer.data(il + 83);
    const auto *il_84 = buffer.data(il + 84);
    const auto *il_85 = buffer.data(il + 85);
    const auto *il_86 = buffer.data(il + 86);
    const auto *il_87 = buffer.data(il + 87);
    const auto *il_88 = buffer.data(il + 88);
    const auto *il_89 = buffer.data(il + 89);
    const auto *il_90 = buffer.data(il + 90);
    const auto *il_91 = buffer.data(il + 91);
    const auto *il_92 = buffer.data(il + 92);
    const auto *il_93 = buffer.data(il + 93);
    const auto *il_94 = buffer.data(il + 94);
    const auto *il_95 = buffer.data(il + 95);
    const auto *il_96 = buffer.data(il + 96);
    const auto *il_97 = buffer.data(il + 97);
    const auto *il_98 = buffer.data(il + 98);
    const auto *il_99 = buffer.data(il + 99);
    const auto *il_100 = buffer.data(il + 100);
    const auto *il_101 = buffer.data(il + 101);
    const auto *il_102 = buffer.data(il + 102);
    const auto *il_103 = buffer.data(il + 103);
    const auto *il_104 = buffer.data(il + 104);
    const auto *il_105 = buffer.data(il + 105);
    const auto *il_106 = buffer.data(il + 106);
    const auto *il_107 = buffer.data(il + 107);
    const auto *il_108 = buffer.data(il + 108);
    const auto *il_109 = buffer.data(il + 109);
    const auto *il_110 = buffer.data(il + 110);
    const auto *il_111 = buffer.data(il + 111);
    const auto *il_112 = buffer.data(il + 112);
    const auto *il_113 = buffer.data(il + 113);
    const auto *il_114 = buffer.data(il + 114);
    const auto *il_115 = buffer.data(il + 115);
    const auto *il_116 = buffer.data(il + 116);
    const auto *il_117 = buffer.data(il + 117);
    const auto *il_118 = buffer.data(il + 118);
    const auto *il_119 = buffer.data(il + 119);
    const auto *il_120 = buffer.data(il + 120);
    const auto *il_121 = buffer.data(il + 121);
    const auto *il_122 = buffer.data(il + 122);
    const auto *il_123 = buffer.data(il + 123);
    const auto *il_124 = buffer.data(il + 124);
    const auto *il_125 = buffer.data(il + 125);
    const auto *il_126 = buffer.data(il + 126);
    const auto *il_127 = buffer.data(il + 127);
    const auto *il_128 = buffer.data(il + 128);
    const auto *il_129 = buffer.data(il + 129);
    const auto *il_130 = buffer.data(il + 130);
    const auto *il_131 = buffer.data(il + 131);
    const auto *il_132 = buffer.data(il + 132);
    const auto *il_133 = buffer.data(il + 133);
    const auto *il_134 = buffer.data(il + 134);
    const auto *il_135 = buffer.data(il + 135);
    const auto *il_136 = buffer.data(il + 136);
    const auto *il_137 = buffer.data(il + 137);
    const auto *il_138 = buffer.data(il + 138);
    const auto *il_139 = buffer.data(il + 139);
    const auto *il_140 = buffer.data(il + 140);
    const auto *il_141 = buffer.data(il + 141);
    const auto *il_142 = buffer.data(il + 142);
    const auto *il_143 = buffer.data(il + 143);
    const auto *il_144 = buffer.data(il + 144);
    const auto *il_145 = buffer.data(il + 145);
    const auto *il_146 = buffer.data(il + 146);
    const auto *il_147 = buffer.data(il + 147);
    const auto *il_148 = buffer.data(il + 148);
    const auto *il_149 = buffer.data(il + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, gl_0, gl_1, gl_2, gl_3, gl_4, il_0, il_1, \
                         il_2, il_3, il_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -5.0 * gl_0[k]
                 + f_0 * il_0[k];

        t_1[k] = -5.0 * gl_1[k]
                 + f_0 * il_1[k];

        t_2[k] = -5.0 * gl_2[k]
                 + f_0 * il_2[k];

        t_3[k] = -5.0 * gl_3[k]
                 + f_0 * il_3[k];

        t_4[k] = -5.0 * gl_4[k]
                 + f_0 * il_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, gl_5, gl_6, gl_7, gl_8, gl_9, il_5, il_6, \
                         il_7, il_8, il_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -5.0 * gl_5[k]
                 + f_0 * il_5[k];

        t_6[k] = -5.0 * gl_6[k]
                 + f_0 * il_6[k];

        t_7[k] = -5.0 * gl_7[k]
                 + f_0 * il_7[k];

        t_8[k] = -5.0 * gl_8[k]
                 + f_0 * il_8[k];

        t_9[k] = -5.0 * gl_9[k]
                 + f_0 * il_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, gl_10, gl_11, gl_12, gl_13, gl_14, \
                         il_10, il_11, il_12, il_13, il_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -5.0 * gl_10[k]
                  + f_0 * il_10[k];

        t_11[k] = -5.0 * gl_11[k]
                  + f_0 * il_11[k];

        t_12[k] = -5.0 * gl_12[k]
                  + f_0 * il_12[k];

        t_13[k] = -5.0 * gl_13[k]
                  + f_0 * il_13[k];

        t_14[k] = -5.0 * gl_14[k]
                  + f_0 * il_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, gl_15, gl_16, gl_17, gl_18, gl_19, \
                         il_15, il_16, il_17, il_18, il_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -5.0 * gl_15[k]
                  + f_0 * il_15[k];

        t_16[k] = -5.0 * gl_16[k]
                  + f_0 * il_16[k];

        t_17[k] = -5.0 * gl_17[k]
                  + f_0 * il_17[k];

        t_18[k] = -5.0 * gl_18[k]
                  + f_0 * il_18[k];

        t_19[k] = -5.0 * gl_19[k]
                  + f_0 * il_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, gl_20, gl_21, gl_22, gl_23, gl_24, \
                         il_20, il_21, il_22, il_23, il_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -5.0 * gl_20[k]
                  + f_0 * il_20[k];

        t_21[k] = -5.0 * gl_21[k]
                  + f_0 * il_21[k];

        t_22[k] = -5.0 * gl_22[k]
                  + f_0 * il_22[k];

        t_23[k] = -5.0 * gl_23[k]
                  + f_0 * il_23[k];

        t_24[k] = -5.0 * gl_24[k]
                  + f_0 * il_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, gl_25, gl_26, gl_27, gl_28, gl_29, \
                         il_25, il_26, il_27, il_28, il_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -5.0 * gl_25[k]
                  + f_0 * il_25[k];

        t_26[k] = -5.0 * gl_26[k]
                  + f_0 * il_26[k];

        t_27[k] = -5.0 * gl_27[k]
                  + f_0 * il_27[k];

        t_28[k] = -5.0 * gl_28[k]
                  + f_0 * il_28[k];

        t_29[k] = -5.0 * gl_29[k]
                  + f_0 * il_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, gl_30, gl_31, gl_32, gl_33, gl_34, \
                         il_30, il_31, il_32, il_33, il_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -5.0 * gl_30[k]
                  + f_0 * il_30[k];

        t_31[k] = -5.0 * gl_31[k]
                  + f_0 * il_31[k];

        t_32[k] = -5.0 * gl_32[k]
                  + f_0 * il_32[k];

        t_33[k] = -5.0 * gl_33[k]
                  + f_0 * il_33[k];

        t_34[k] = -5.0 * gl_34[k]
                  + f_0 * il_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, gl_35, gl_36, gl_37, gl_38, gl_39, \
                         il_35, il_36, il_37, il_38, il_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -5.0 * gl_35[k]
                  + f_0 * il_35[k];

        t_36[k] = -5.0 * gl_36[k]
                  + f_0 * il_36[k];

        t_37[k] = -5.0 * gl_37[k]
                  + f_0 * il_37[k];

        t_38[k] = -5.0 * gl_38[k]
                  + f_0 * il_38[k];

        t_39[k] = -5.0 * gl_39[k]
                  + f_0 * il_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, gl_40, gl_41, gl_42, gl_43, gl_44, \
                         il_40, il_41, il_42, il_43, il_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -5.0 * gl_40[k]
                  + f_0 * il_40[k];

        t_41[k] = -5.0 * gl_41[k]
                  + f_0 * il_41[k];

        t_42[k] = -5.0 * gl_42[k]
                  + f_0 * il_42[k];

        t_43[k] = -5.0 * gl_43[k]
                  + f_0 * il_43[k];

        t_44[k] = -5.0 * gl_44[k]
                  + f_0 * il_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, gl_45, gl_46, gl_47, gl_48, gl_49, \
                         il_45, il_46, il_47, il_48, il_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -4.0 * gl_45[k]
                  + f_0 * il_45[k];

        t_46[k] = -4.0 * gl_46[k]
                  + f_0 * il_46[k];

        t_47[k] = -4.0 * gl_47[k]
                  + f_0 * il_47[k];

        t_48[k] = -4.0 * gl_48[k]
                  + f_0 * il_48[k];

        t_49[k] = -4.0 * gl_49[k]
                  + f_0 * il_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, gl_50, gl_51, gl_52, gl_53, gl_54, \
                         il_50, il_51, il_52, il_53, il_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -4.0 * gl_50[k]
                  + f_0 * il_50[k];

        t_51[k] = -4.0 * gl_51[k]
                  + f_0 * il_51[k];

        t_52[k] = -4.0 * gl_52[k]
                  + f_0 * il_52[k];

        t_53[k] = -4.0 * gl_53[k]
                  + f_0 * il_53[k];

        t_54[k] = -4.0 * gl_54[k]
                  + f_0 * il_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, gl_55, gl_56, gl_57, gl_58, gl_59, \
                         il_55, il_56, il_57, il_58, il_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -4.0 * gl_55[k]
                  + f_0 * il_55[k];

        t_56[k] = -4.0 * gl_56[k]
                  + f_0 * il_56[k];

        t_57[k] = -4.0 * gl_57[k]
                  + f_0 * il_57[k];

        t_58[k] = -4.0 * gl_58[k]
                  + f_0 * il_58[k];

        t_59[k] = -4.0 * gl_59[k]
                  + f_0 * il_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, gl_60, gl_61, gl_62, gl_63, gl_64, \
                         il_60, il_61, il_62, il_63, il_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -4.0 * gl_60[k]
                  + f_0 * il_60[k];

        t_61[k] = -4.0 * gl_61[k]
                  + f_0 * il_61[k];

        t_62[k] = -4.0 * gl_62[k]
                  + f_0 * il_62[k];

        t_63[k] = -4.0 * gl_63[k]
                  + f_0 * il_63[k];

        t_64[k] = -4.0 * gl_64[k]
                  + f_0 * il_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, gl_65, gl_66, gl_67, gl_68, gl_69, \
                         il_65, il_66, il_67, il_68, il_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -4.0 * gl_65[k]
                  + f_0 * il_65[k];

        t_66[k] = -4.0 * gl_66[k]
                  + f_0 * il_66[k];

        t_67[k] = -4.0 * gl_67[k]
                  + f_0 * il_67[k];

        t_68[k] = -4.0 * gl_68[k]
                  + f_0 * il_68[k];

        t_69[k] = -4.0 * gl_69[k]
                  + f_0 * il_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, gl_70, gl_71, gl_72, gl_73, gl_74, \
                         il_70, il_71, il_72, il_73, il_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -4.0 * gl_70[k]
                  + f_0 * il_70[k];

        t_71[k] = -4.0 * gl_71[k]
                  + f_0 * il_71[k];

        t_72[k] = -4.0 * gl_72[k]
                  + f_0 * il_72[k];

        t_73[k] = -4.0 * gl_73[k]
                  + f_0 * il_73[k];

        t_74[k] = -4.0 * gl_74[k]
                  + f_0 * il_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, gl_75, gl_76, gl_77, gl_78, gl_79, \
                         il_75, il_76, il_77, il_78, il_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -4.0 * gl_75[k]
                  + f_0 * il_75[k];

        t_76[k] = -4.0 * gl_76[k]
                  + f_0 * il_76[k];

        t_77[k] = -4.0 * gl_77[k]
                  + f_0 * il_77[k];

        t_78[k] = -4.0 * gl_78[k]
                  + f_0 * il_78[k];

        t_79[k] = -4.0 * gl_79[k]
                  + f_0 * il_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, gl_80, gl_81, gl_82, gl_83, gl_84, \
                         il_80, il_81, il_82, il_83, il_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -4.0 * gl_80[k]
                  + f_0 * il_80[k];

        t_81[k] = -4.0 * gl_81[k]
                  + f_0 * il_81[k];

        t_82[k] = -4.0 * gl_82[k]
                  + f_0 * il_82[k];

        t_83[k] = -4.0 * gl_83[k]
                  + f_0 * il_83[k];

        t_84[k] = -4.0 * gl_84[k]
                  + f_0 * il_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, gl_85, gl_86, gl_87, gl_88, gl_89, \
                         il_85, il_86, il_87, il_88, il_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * gl_85[k]
                  + f_0 * il_85[k];

        t_86[k] = -4.0 * gl_86[k]
                  + f_0 * il_86[k];

        t_87[k] = -4.0 * gl_87[k]
                  + f_0 * il_87[k];

        t_88[k] = -4.0 * gl_88[k]
                  + f_0 * il_88[k];

        t_89[k] = -4.0 * gl_89[k]
                  + f_0 * il_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, gl_90, gl_91, gl_92, gl_93, gl_94, \
                         il_90, il_91, il_92, il_93, il_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -4.0 * gl_90[k]
                  + f_0 * il_90[k];

        t_91[k] = -4.0 * gl_91[k]
                  + f_0 * il_91[k];

        t_92[k] = -4.0 * gl_92[k]
                  + f_0 * il_92[k];

        t_93[k] = -4.0 * gl_93[k]
                  + f_0 * il_93[k];

        t_94[k] = -4.0 * gl_94[k]
                  + f_0 * il_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, gl_95, gl_96, gl_97, gl_98, gl_99, \
                         il_95, il_96, il_97, il_98, il_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -4.0 * gl_95[k]
                  + f_0 * il_95[k];

        t_96[k] = -4.0 * gl_96[k]
                  + f_0 * il_96[k];

        t_97[k] = -4.0 * gl_97[k]
                  + f_0 * il_97[k];

        t_98[k] = -4.0 * gl_98[k]
                  + f_0 * il_98[k];

        t_99[k] = -4.0 * gl_99[k]
                  + f_0 * il_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, gl_100, gl_101, gl_102, gl_103, \
                         gl_104, il_100, il_101, il_102, il_103, \
                         il_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -4.0 * gl_100[k]
                   + f_0 * il_100[k];

        t_101[k] = -4.0 * gl_101[k]
                   + f_0 * il_101[k];

        t_102[k] = -4.0 * gl_102[k]
                   + f_0 * il_102[k];

        t_103[k] = -4.0 * gl_103[k]
                   + f_0 * il_103[k];

        t_104[k] = -4.0 * gl_104[k]
                   + f_0 * il_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, gl_105, gl_106, gl_107, gl_108, \
                         gl_109, il_105, il_106, il_107, il_108, \
                         il_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -4.0 * gl_105[k]
                   + f_0 * il_105[k];

        t_106[k] = -4.0 * gl_106[k]
                   + f_0 * il_106[k];

        t_107[k] = -4.0 * gl_107[k]
                   + f_0 * il_107[k];

        t_108[k] = -4.0 * gl_108[k]
                   + f_0 * il_108[k];

        t_109[k] = -4.0 * gl_109[k]
                   + f_0 * il_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, gl_110, gl_111, gl_112, gl_113, \
                         gl_114, il_110, il_111, il_112, il_113, \
                         il_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -4.0 * gl_110[k]
                   + f_0 * il_110[k];

        t_111[k] = -4.0 * gl_111[k]
                   + f_0 * il_111[k];

        t_112[k] = -4.0 * gl_112[k]
                   + f_0 * il_112[k];

        t_113[k] = -4.0 * gl_113[k]
                   + f_0 * il_113[k];

        t_114[k] = -4.0 * gl_114[k]
                   + f_0 * il_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, gl_115, gl_116, gl_117, gl_118, \
                         gl_119, il_115, il_116, il_117, il_118, \
                         il_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -4.0 * gl_115[k]
                   + f_0 * il_115[k];

        t_116[k] = -4.0 * gl_116[k]
                   + f_0 * il_116[k];

        t_117[k] = -4.0 * gl_117[k]
                   + f_0 * il_117[k];

        t_118[k] = -4.0 * gl_118[k]
                   + f_0 * il_118[k];

        t_119[k] = -4.0 * gl_119[k]
                   + f_0 * il_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, gl_120, gl_121, gl_122, gl_123, \
                         gl_124, il_120, il_121, il_122, il_123, \
                         il_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -4.0 * gl_120[k]
                   + f_0 * il_120[k];

        t_121[k] = -4.0 * gl_121[k]
                   + f_0 * il_121[k];

        t_122[k] = -4.0 * gl_122[k]
                   + f_0 * il_122[k];

        t_123[k] = -4.0 * gl_123[k]
                   + f_0 * il_123[k];

        t_124[k] = -4.0 * gl_124[k]
                   + f_0 * il_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, gl_125, gl_126, gl_127, gl_128, \
                         gl_129, il_125, il_126, il_127, il_128, \
                         il_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -4.0 * gl_125[k]
                   + f_0 * il_125[k];

        t_126[k] = -4.0 * gl_126[k]
                   + f_0 * il_126[k];

        t_127[k] = -4.0 * gl_127[k]
                   + f_0 * il_127[k];

        t_128[k] = -4.0 * gl_128[k]
                   + f_0 * il_128[k];

        t_129[k] = -4.0 * gl_129[k]
                   + f_0 * il_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, gl_130, gl_131, gl_132, gl_133, \
                         gl_134, il_130, il_131, il_132, il_133, \
                         il_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -4.0 * gl_130[k]
                   + f_0 * il_130[k];

        t_131[k] = -4.0 * gl_131[k]
                   + f_0 * il_131[k];

        t_132[k] = -4.0 * gl_132[k]
                   + f_0 * il_132[k];

        t_133[k] = -4.0 * gl_133[k]
                   + f_0 * il_133[k];

        t_134[k] = -4.0 * gl_134[k]
                   + f_0 * il_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, gl_135, gl_136, gl_137, gl_138, \
                         gl_139, il_135, il_136, il_137, il_138, \
                         il_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -3.0 * gl_135[k]
                   + f_0 * il_135[k];

        t_136[k] = -3.0 * gl_136[k]
                   + f_0 * il_136[k];

        t_137[k] = -3.0 * gl_137[k]
                   + f_0 * il_137[k];

        t_138[k] = -3.0 * gl_138[k]
                   + f_0 * il_138[k];

        t_139[k] = -3.0 * gl_139[k]
                   + f_0 * il_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, gl_140, gl_141, gl_142, gl_143, \
                         gl_144, il_140, il_141, il_142, il_143, \
                         il_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -3.0 * gl_140[k]
                   + f_0 * il_140[k];

        t_141[k] = -3.0 * gl_141[k]
                   + f_0 * il_141[k];

        t_142[k] = -3.0 * gl_142[k]
                   + f_0 * il_142[k];

        t_143[k] = -3.0 * gl_143[k]
                   + f_0 * il_143[k];

        t_144[k] = -3.0 * gl_144[k]
                   + f_0 * il_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, gl_145, gl_146, gl_147, gl_148, \
                         gl_149, il_145, il_146, il_147, il_148, \
                         il_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -3.0 * gl_145[k]
                   + f_0 * il_145[k];

        t_146[k] = -3.0 * gl_146[k]
                   + f_0 * il_146[k];

        t_147[k] = -3.0 * gl_147[k]
                   + f_0 * il_147[k];

        t_148[k] = -3.0 * gl_148[k]
                   + f_0 * il_148[k];

        t_149[k] = -3.0 * gl_149[k]
                   + f_0 * il_149[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_150 = buffer.data(gl + 150);
    const auto *gl_151 = buffer.data(gl + 151);
    const auto *gl_152 = buffer.data(gl + 152);
    const auto *gl_153 = buffer.data(gl + 153);
    const auto *gl_154 = buffer.data(gl + 154);
    const auto *gl_155 = buffer.data(gl + 155);
    const auto *gl_156 = buffer.data(gl + 156);
    const auto *gl_157 = buffer.data(gl + 157);
    const auto *gl_158 = buffer.data(gl + 158);
    const auto *gl_159 = buffer.data(gl + 159);
    const auto *gl_160 = buffer.data(gl + 160);
    const auto *gl_161 = buffer.data(gl + 161);
    const auto *gl_162 = buffer.data(gl + 162);
    const auto *gl_163 = buffer.data(gl + 163);
    const auto *gl_164 = buffer.data(gl + 164);
    const auto *gl_165 = buffer.data(gl + 165);
    const auto *gl_166 = buffer.data(gl + 166);
    const auto *gl_167 = buffer.data(gl + 167);
    const auto *gl_168 = buffer.data(gl + 168);
    const auto *gl_169 = buffer.data(gl + 169);
    const auto *gl_170 = buffer.data(gl + 170);
    const auto *gl_171 = buffer.data(gl + 171);
    const auto *gl_172 = buffer.data(gl + 172);
    const auto *gl_173 = buffer.data(gl + 173);
    const auto *gl_174 = buffer.data(gl + 174);
    const auto *gl_175 = buffer.data(gl + 175);
    const auto *gl_176 = buffer.data(gl + 176);
    const auto *gl_177 = buffer.data(gl + 177);
    const auto *gl_178 = buffer.data(gl + 178);
    const auto *gl_179 = buffer.data(gl + 179);
    const auto *gl_180 = buffer.data(gl + 180);
    const auto *gl_181 = buffer.data(gl + 181);
    const auto *gl_182 = buffer.data(gl + 182);
    const auto *gl_183 = buffer.data(gl + 183);
    const auto *gl_184 = buffer.data(gl + 184);
    const auto *gl_185 = buffer.data(gl + 185);
    const auto *gl_186 = buffer.data(gl + 186);
    const auto *gl_187 = buffer.data(gl + 187);
    const auto *gl_188 = buffer.data(gl + 188);
    const auto *gl_189 = buffer.data(gl + 189);
    const auto *gl_190 = buffer.data(gl + 190);
    const auto *gl_191 = buffer.data(gl + 191);
    const auto *gl_192 = buffer.data(gl + 192);
    const auto *gl_193 = buffer.data(gl + 193);
    const auto *gl_194 = buffer.data(gl + 194);
    const auto *gl_195 = buffer.data(gl + 195);
    const auto *gl_196 = buffer.data(gl + 196);
    const auto *gl_197 = buffer.data(gl + 197);
    const auto *gl_198 = buffer.data(gl + 198);
    const auto *gl_199 = buffer.data(gl + 199);
    const auto *gl_200 = buffer.data(gl + 200);
    const auto *gl_201 = buffer.data(gl + 201);
    const auto *gl_202 = buffer.data(gl + 202);
    const auto *gl_203 = buffer.data(gl + 203);
    const auto *gl_204 = buffer.data(gl + 204);
    const auto *gl_205 = buffer.data(gl + 205);
    const auto *gl_206 = buffer.data(gl + 206);
    const auto *gl_207 = buffer.data(gl + 207);
    const auto *gl_208 = buffer.data(gl + 208);
    const auto *gl_209 = buffer.data(gl + 209);
    const auto *gl_210 = buffer.data(gl + 210);
    const auto *gl_211 = buffer.data(gl + 211);
    const auto *gl_212 = buffer.data(gl + 212);
    const auto *gl_213 = buffer.data(gl + 213);
    const auto *gl_214 = buffer.data(gl + 214);
    const auto *gl_215 = buffer.data(gl + 215);
    const auto *gl_216 = buffer.data(gl + 216);
    const auto *gl_217 = buffer.data(gl + 217);
    const auto *gl_218 = buffer.data(gl + 218);
    const auto *gl_219 = buffer.data(gl + 219);
    const auto *gl_220 = buffer.data(gl + 220);
    const auto *gl_221 = buffer.data(gl + 221);
    const auto *gl_222 = buffer.data(gl + 222);
    const auto *gl_223 = buffer.data(gl + 223);
    const auto *gl_224 = buffer.data(gl + 224);
    const auto *gl_225 = buffer.data(gl + 225);
    const auto *gl_226 = buffer.data(gl + 226);
    const auto *gl_227 = buffer.data(gl + 227);
    const auto *gl_228 = buffer.data(gl + 228);
    const auto *gl_229 = buffer.data(gl + 229);
    const auto *gl_230 = buffer.data(gl + 230);
    const auto *gl_231 = buffer.data(gl + 231);
    const auto *gl_232 = buffer.data(gl + 232);
    const auto *gl_233 = buffer.data(gl + 233);
    const auto *gl_234 = buffer.data(gl + 234);
    const auto *gl_235 = buffer.data(gl + 235);
    const auto *gl_236 = buffer.data(gl + 236);
    const auto *gl_237 = buffer.data(gl + 237);
    const auto *gl_238 = buffer.data(gl + 238);
    const auto *gl_239 = buffer.data(gl + 239);
    const auto *gl_240 = buffer.data(gl + 240);
    const auto *gl_241 = buffer.data(gl + 241);
    const auto *gl_242 = buffer.data(gl + 242);
    const auto *gl_243 = buffer.data(gl + 243);
    const auto *gl_244 = buffer.data(gl + 244);
    const auto *gl_245 = buffer.data(gl + 245);
    const auto *gl_246 = buffer.data(gl + 246);
    const auto *gl_247 = buffer.data(gl + 247);
    const auto *gl_248 = buffer.data(gl + 248);
    const auto *gl_249 = buffer.data(gl + 249);
    const auto *gl_250 = buffer.data(gl + 250);
    const auto *gl_251 = buffer.data(gl + 251);
    const auto *gl_252 = buffer.data(gl + 252);
    const auto *gl_253 = buffer.data(gl + 253);
    const auto *gl_254 = buffer.data(gl + 254);
    const auto *gl_255 = buffer.data(gl + 255);
    const auto *gl_256 = buffer.data(gl + 256);
    const auto *gl_257 = buffer.data(gl + 257);
    const auto *gl_258 = buffer.data(gl + 258);
    const auto *gl_259 = buffer.data(gl + 259);
    const auto *gl_260 = buffer.data(gl + 260);
    const auto *gl_261 = buffer.data(gl + 261);
    const auto *gl_262 = buffer.data(gl + 262);
    const auto *gl_263 = buffer.data(gl + 263);
    const auto *gl_264 = buffer.data(gl + 264);
    const auto *gl_265 = buffer.data(gl + 265);
    const auto *gl_266 = buffer.data(gl + 266);
    const auto *gl_267 = buffer.data(gl + 267);
    const auto *gl_268 = buffer.data(gl + 268);
    const auto *gl_269 = buffer.data(gl + 269);
    const auto *gl_270 = buffer.data(gl + 270);
    const auto *gl_271 = buffer.data(gl + 271);
    const auto *gl_272 = buffer.data(gl + 272);
    const auto *gl_273 = buffer.data(gl + 273);
    const auto *gl_274 = buffer.data(gl + 274);
    const auto *gl_275 = buffer.data(gl + 275);
    const auto *gl_276 = buffer.data(gl + 276);
    const auto *gl_277 = buffer.data(gl + 277);
    const auto *gl_278 = buffer.data(gl + 278);
    const auto *gl_279 = buffer.data(gl + 279);
    const auto *gl_280 = buffer.data(gl + 280);
    const auto *gl_281 = buffer.data(gl + 281);
    const auto *gl_282 = buffer.data(gl + 282);
    const auto *gl_283 = buffer.data(gl + 283);
    const auto *gl_284 = buffer.data(gl + 284);
    const auto *gl_285 = buffer.data(gl + 285);
    const auto *gl_286 = buffer.data(gl + 286);
    const auto *gl_287 = buffer.data(gl + 287);
    const auto *gl_288 = buffer.data(gl + 288);
    const auto *gl_289 = buffer.data(gl + 289);
    const auto *gl_290 = buffer.data(gl + 290);
    const auto *gl_291 = buffer.data(gl + 291);
    const auto *gl_292 = buffer.data(gl + 292);
    const auto *gl_293 = buffer.data(gl + 293);
    const auto *gl_294 = buffer.data(gl + 294);
    const auto *gl_295 = buffer.data(gl + 295);
    const auto *gl_296 = buffer.data(gl + 296);
    const auto *gl_297 = buffer.data(gl + 297);
    const auto *gl_298 = buffer.data(gl + 298);
    const auto *gl_299 = buffer.data(gl + 299);

    const auto *il_150 = buffer.data(il + 150);
    const auto *il_151 = buffer.data(il + 151);
    const auto *il_152 = buffer.data(il + 152);
    const auto *il_153 = buffer.data(il + 153);
    const auto *il_154 = buffer.data(il + 154);
    const auto *il_155 = buffer.data(il + 155);
    const auto *il_156 = buffer.data(il + 156);
    const auto *il_157 = buffer.data(il + 157);
    const auto *il_158 = buffer.data(il + 158);
    const auto *il_159 = buffer.data(il + 159);
    const auto *il_160 = buffer.data(il + 160);
    const auto *il_161 = buffer.data(il + 161);
    const auto *il_162 = buffer.data(il + 162);
    const auto *il_163 = buffer.data(il + 163);
    const auto *il_164 = buffer.data(il + 164);
    const auto *il_165 = buffer.data(il + 165);
    const auto *il_166 = buffer.data(il + 166);
    const auto *il_167 = buffer.data(il + 167);
    const auto *il_168 = buffer.data(il + 168);
    const auto *il_169 = buffer.data(il + 169);
    const auto *il_170 = buffer.data(il + 170);
    const auto *il_171 = buffer.data(il + 171);
    const auto *il_172 = buffer.data(il + 172);
    const auto *il_173 = buffer.data(il + 173);
    const auto *il_174 = buffer.data(il + 174);
    const auto *il_175 = buffer.data(il + 175);
    const auto *il_176 = buffer.data(il + 176);
    const auto *il_177 = buffer.data(il + 177);
    const auto *il_178 = buffer.data(il + 178);
    const auto *il_179 = buffer.data(il + 179);
    const auto *il_180 = buffer.data(il + 180);
    const auto *il_181 = buffer.data(il + 181);
    const auto *il_182 = buffer.data(il + 182);
    const auto *il_183 = buffer.data(il + 183);
    const auto *il_184 = buffer.data(il + 184);
    const auto *il_185 = buffer.data(il + 185);
    const auto *il_186 = buffer.data(il + 186);
    const auto *il_187 = buffer.data(il + 187);
    const auto *il_188 = buffer.data(il + 188);
    const auto *il_189 = buffer.data(il + 189);
    const auto *il_190 = buffer.data(il + 190);
    const auto *il_191 = buffer.data(il + 191);
    const auto *il_192 = buffer.data(il + 192);
    const auto *il_193 = buffer.data(il + 193);
    const auto *il_194 = buffer.data(il + 194);
    const auto *il_195 = buffer.data(il + 195);
    const auto *il_196 = buffer.data(il + 196);
    const auto *il_197 = buffer.data(il + 197);
    const auto *il_198 = buffer.data(il + 198);
    const auto *il_199 = buffer.data(il + 199);
    const auto *il_200 = buffer.data(il + 200);
    const auto *il_201 = buffer.data(il + 201);
    const auto *il_202 = buffer.data(il + 202);
    const auto *il_203 = buffer.data(il + 203);
    const auto *il_204 = buffer.data(il + 204);
    const auto *il_205 = buffer.data(il + 205);
    const auto *il_206 = buffer.data(il + 206);
    const auto *il_207 = buffer.data(il + 207);
    const auto *il_208 = buffer.data(il + 208);
    const auto *il_209 = buffer.data(il + 209);
    const auto *il_210 = buffer.data(il + 210);
    const auto *il_211 = buffer.data(il + 211);
    const auto *il_212 = buffer.data(il + 212);
    const auto *il_213 = buffer.data(il + 213);
    const auto *il_214 = buffer.data(il + 214);
    const auto *il_215 = buffer.data(il + 215);
    const auto *il_216 = buffer.data(il + 216);
    const auto *il_217 = buffer.data(il + 217);
    const auto *il_218 = buffer.data(il + 218);
    const auto *il_219 = buffer.data(il + 219);
    const auto *il_220 = buffer.data(il + 220);
    const auto *il_221 = buffer.data(il + 221);
    const auto *il_222 = buffer.data(il + 222);
    const auto *il_223 = buffer.data(il + 223);
    const auto *il_224 = buffer.data(il + 224);
    const auto *il_225 = buffer.data(il + 225);
    const auto *il_226 = buffer.data(il + 226);
    const auto *il_227 = buffer.data(il + 227);
    const auto *il_228 = buffer.data(il + 228);
    const auto *il_229 = buffer.data(il + 229);
    const auto *il_230 = buffer.data(il + 230);
    const auto *il_231 = buffer.data(il + 231);
    const auto *il_232 = buffer.data(il + 232);
    const auto *il_233 = buffer.data(il + 233);
    const auto *il_234 = buffer.data(il + 234);
    const auto *il_235 = buffer.data(il + 235);
    const auto *il_236 = buffer.data(il + 236);
    const auto *il_237 = buffer.data(il + 237);
    const auto *il_238 = buffer.data(il + 238);
    const auto *il_239 = buffer.data(il + 239);
    const auto *il_240 = buffer.data(il + 240);
    const auto *il_241 = buffer.data(il + 241);
    const auto *il_242 = buffer.data(il + 242);
    const auto *il_243 = buffer.data(il + 243);
    const auto *il_244 = buffer.data(il + 244);
    const auto *il_245 = buffer.data(il + 245);
    const auto *il_246 = buffer.data(il + 246);
    const auto *il_247 = buffer.data(il + 247);
    const auto *il_248 = buffer.data(il + 248);
    const auto *il_249 = buffer.data(il + 249);
    const auto *il_250 = buffer.data(il + 250);
    const auto *il_251 = buffer.data(il + 251);
    const auto *il_252 = buffer.data(il + 252);
    const auto *il_253 = buffer.data(il + 253);
    const auto *il_254 = buffer.data(il + 254);
    const auto *il_255 = buffer.data(il + 255);
    const auto *il_256 = buffer.data(il + 256);
    const auto *il_257 = buffer.data(il + 257);
    const auto *il_258 = buffer.data(il + 258);
    const auto *il_259 = buffer.data(il + 259);
    const auto *il_260 = buffer.data(il + 260);
    const auto *il_261 = buffer.data(il + 261);
    const auto *il_262 = buffer.data(il + 262);
    const auto *il_263 = buffer.data(il + 263);
    const auto *il_264 = buffer.data(il + 264);
    const auto *il_265 = buffer.data(il + 265);
    const auto *il_266 = buffer.data(il + 266);
    const auto *il_267 = buffer.data(il + 267);
    const auto *il_268 = buffer.data(il + 268);
    const auto *il_269 = buffer.data(il + 269);
    const auto *il_270 = buffer.data(il + 270);
    const auto *il_271 = buffer.data(il + 271);
    const auto *il_272 = buffer.data(il + 272);
    const auto *il_273 = buffer.data(il + 273);
    const auto *il_274 = buffer.data(il + 274);
    const auto *il_275 = buffer.data(il + 275);
    const auto *il_276 = buffer.data(il + 276);
    const auto *il_277 = buffer.data(il + 277);
    const auto *il_278 = buffer.data(il + 278);
    const auto *il_279 = buffer.data(il + 279);
    const auto *il_280 = buffer.data(il + 280);
    const auto *il_281 = buffer.data(il + 281);
    const auto *il_282 = buffer.data(il + 282);
    const auto *il_283 = buffer.data(il + 283);
    const auto *il_284 = buffer.data(il + 284);
    const auto *il_285 = buffer.data(il + 285);
    const auto *il_286 = buffer.data(il + 286);
    const auto *il_287 = buffer.data(il + 287);
    const auto *il_288 = buffer.data(il + 288);
    const auto *il_289 = buffer.data(il + 289);
    const auto *il_290 = buffer.data(il + 290);
    const auto *il_291 = buffer.data(il + 291);
    const auto *il_292 = buffer.data(il + 292);
    const auto *il_293 = buffer.data(il + 293);
    const auto *il_294 = buffer.data(il + 294);
    const auto *il_295 = buffer.data(il + 295);
    const auto *il_296 = buffer.data(il + 296);
    const auto *il_297 = buffer.data(il + 297);
    const auto *il_298 = buffer.data(il + 298);
    const auto *il_299 = buffer.data(il + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, gl_150, gl_151, gl_152, gl_153, \
                         gl_154, il_150, il_151, il_152, il_153, \
                         il_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -3.0 * gl_150[k]
                   + f_0 * il_150[k];

        t_151[k] = -3.0 * gl_151[k]
                   + f_0 * il_151[k];

        t_152[k] = -3.0 * gl_152[k]
                   + f_0 * il_152[k];

        t_153[k] = -3.0 * gl_153[k]
                   + f_0 * il_153[k];

        t_154[k] = -3.0 * gl_154[k]
                   + f_0 * il_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, gl_155, gl_156, gl_157, gl_158, \
                         gl_159, il_155, il_156, il_157, il_158, \
                         il_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -3.0 * gl_155[k]
                   + f_0 * il_155[k];

        t_156[k] = -3.0 * gl_156[k]
                   + f_0 * il_156[k];

        t_157[k] = -3.0 * gl_157[k]
                   + f_0 * il_157[k];

        t_158[k] = -3.0 * gl_158[k]
                   + f_0 * il_158[k];

        t_159[k] = -3.0 * gl_159[k]
                   + f_0 * il_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, gl_160, gl_161, gl_162, gl_163, \
                         gl_164, il_160, il_161, il_162, il_163, \
                         il_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -3.0 * gl_160[k]
                   + f_0 * il_160[k];

        t_161[k] = -3.0 * gl_161[k]
                   + f_0 * il_161[k];

        t_162[k] = -3.0 * gl_162[k]
                   + f_0 * il_162[k];

        t_163[k] = -3.0 * gl_163[k]
                   + f_0 * il_163[k];

        t_164[k] = -3.0 * gl_164[k]
                   + f_0 * il_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, gl_165, gl_166, gl_167, gl_168, \
                         gl_169, il_165, il_166, il_167, il_168, \
                         il_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -3.0 * gl_165[k]
                   + f_0 * il_165[k];

        t_166[k] = -3.0 * gl_166[k]
                   + f_0 * il_166[k];

        t_167[k] = -3.0 * gl_167[k]
                   + f_0 * il_167[k];

        t_168[k] = -3.0 * gl_168[k]
                   + f_0 * il_168[k];

        t_169[k] = -3.0 * gl_169[k]
                   + f_0 * il_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, gl_170, gl_171, gl_172, gl_173, \
                         gl_174, il_170, il_171, il_172, il_173, \
                         il_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * gl_170[k]
                   + f_0 * il_170[k];

        t_171[k] = -3.0 * gl_171[k]
                   + f_0 * il_171[k];

        t_172[k] = -3.0 * gl_172[k]
                   + f_0 * il_172[k];

        t_173[k] = -3.0 * gl_173[k]
                   + f_0 * il_173[k];

        t_174[k] = -3.0 * gl_174[k]
                   + f_0 * il_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, gl_175, gl_176, gl_177, gl_178, \
                         gl_179, il_175, il_176, il_177, il_178, \
                         il_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * gl_175[k]
                   + f_0 * il_175[k];

        t_176[k] = -3.0 * gl_176[k]
                   + f_0 * il_176[k];

        t_177[k] = -3.0 * gl_177[k]
                   + f_0 * il_177[k];

        t_178[k] = -3.0 * gl_178[k]
                   + f_0 * il_178[k];

        t_179[k] = -3.0 * gl_179[k]
                   + f_0 * il_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, gl_180, gl_181, gl_182, gl_183, \
                         gl_184, il_180, il_181, il_182, il_183, \
                         il_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * gl_180[k]
                   + f_0 * il_180[k];

        t_181[k] = -3.0 * gl_181[k]
                   + f_0 * il_181[k];

        t_182[k] = -3.0 * gl_182[k]
                   + f_0 * il_182[k];

        t_183[k] = -3.0 * gl_183[k]
                   + f_0 * il_183[k];

        t_184[k] = -3.0 * gl_184[k]
                   + f_0 * il_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, gl_185, gl_186, gl_187, gl_188, \
                         gl_189, il_185, il_186, il_187, il_188, \
                         il_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * gl_185[k]
                   + f_0 * il_185[k];

        t_186[k] = -3.0 * gl_186[k]
                   + f_0 * il_186[k];

        t_187[k] = -3.0 * gl_187[k]
                   + f_0 * il_187[k];

        t_188[k] = -3.0 * gl_188[k]
                   + f_0 * il_188[k];

        t_189[k] = -3.0 * gl_189[k]
                   + f_0 * il_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, gl_190, gl_191, gl_192, gl_193, \
                         gl_194, il_190, il_191, il_192, il_193, \
                         il_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * gl_190[k]
                   + f_0 * il_190[k];

        t_191[k] = -3.0 * gl_191[k]
                   + f_0 * il_191[k];

        t_192[k] = -3.0 * gl_192[k]
                   + f_0 * il_192[k];

        t_193[k] = -3.0 * gl_193[k]
                   + f_0 * il_193[k];

        t_194[k] = -3.0 * gl_194[k]
                   + f_0 * il_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, gl_195, gl_196, gl_197, gl_198, \
                         gl_199, il_195, il_196, il_197, il_198, \
                         il_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * gl_195[k]
                   + f_0 * il_195[k];

        t_196[k] = -3.0 * gl_196[k]
                   + f_0 * il_196[k];

        t_197[k] = -3.0 * gl_197[k]
                   + f_0 * il_197[k];

        t_198[k] = -3.0 * gl_198[k]
                   + f_0 * il_198[k];

        t_199[k] = -3.0 * gl_199[k]
                   + f_0 * il_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, gl_200, gl_201, gl_202, gl_203, \
                         gl_204, il_200, il_201, il_202, il_203, \
                         il_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -3.0 * gl_200[k]
                   + f_0 * il_200[k];

        t_201[k] = -3.0 * gl_201[k]
                   + f_0 * il_201[k];

        t_202[k] = -3.0 * gl_202[k]
                   + f_0 * il_202[k];

        t_203[k] = -3.0 * gl_203[k]
                   + f_0 * il_203[k];

        t_204[k] = -3.0 * gl_204[k]
                   + f_0 * il_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, gl_205, gl_206, gl_207, gl_208, \
                         gl_209, il_205, il_206, il_207, il_208, \
                         il_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -3.0 * gl_205[k]
                   + f_0 * il_205[k];

        t_206[k] = -3.0 * gl_206[k]
                   + f_0 * il_206[k];

        t_207[k] = -3.0 * gl_207[k]
                   + f_0 * il_207[k];

        t_208[k] = -3.0 * gl_208[k]
                   + f_0 * il_208[k];

        t_209[k] = -3.0 * gl_209[k]
                   + f_0 * il_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, gl_210, gl_211, gl_212, gl_213, \
                         gl_214, il_210, il_211, il_212, il_213, \
                         il_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -3.0 * gl_210[k]
                   + f_0 * il_210[k];

        t_211[k] = -3.0 * gl_211[k]
                   + f_0 * il_211[k];

        t_212[k] = -3.0 * gl_212[k]
                   + f_0 * il_212[k];

        t_213[k] = -3.0 * gl_213[k]
                   + f_0 * il_213[k];

        t_214[k] = -3.0 * gl_214[k]
                   + f_0 * il_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, gl_215, gl_216, gl_217, gl_218, \
                         gl_219, il_215, il_216, il_217, il_218, \
                         il_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -3.0 * gl_215[k]
                   + f_0 * il_215[k];

        t_216[k] = -3.0 * gl_216[k]
                   + f_0 * il_216[k];

        t_217[k] = -3.0 * gl_217[k]
                   + f_0 * il_217[k];

        t_218[k] = -3.0 * gl_218[k]
                   + f_0 * il_218[k];

        t_219[k] = -3.0 * gl_219[k]
                   + f_0 * il_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, gl_220, gl_221, gl_222, gl_223, \
                         gl_224, il_220, il_221, il_222, il_223, \
                         il_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -3.0 * gl_220[k]
                   + f_0 * il_220[k];

        t_221[k] = -3.0 * gl_221[k]
                   + f_0 * il_221[k];

        t_222[k] = -3.0 * gl_222[k]
                   + f_0 * il_222[k];

        t_223[k] = -3.0 * gl_223[k]
                   + f_0 * il_223[k];

        t_224[k] = -3.0 * gl_224[k]
                   + f_0 * il_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, gl_225, gl_226, gl_227, gl_228, \
                         gl_229, il_225, il_226, il_227, il_228, \
                         il_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -3.0 * gl_225[k]
                   + f_0 * il_225[k];

        t_226[k] = -3.0 * gl_226[k]
                   + f_0 * il_226[k];

        t_227[k] = -3.0 * gl_227[k]
                   + f_0 * il_227[k];

        t_228[k] = -3.0 * gl_228[k]
                   + f_0 * il_228[k];

        t_229[k] = -3.0 * gl_229[k]
                   + f_0 * il_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, gl_230, gl_231, gl_232, gl_233, \
                         gl_234, il_230, il_231, il_232, il_233, \
                         il_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -3.0 * gl_230[k]
                   + f_0 * il_230[k];

        t_231[k] = -3.0 * gl_231[k]
                   + f_0 * il_231[k];

        t_232[k] = -3.0 * gl_232[k]
                   + f_0 * il_232[k];

        t_233[k] = -3.0 * gl_233[k]
                   + f_0 * il_233[k];

        t_234[k] = -3.0 * gl_234[k]
                   + f_0 * il_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, gl_235, gl_236, gl_237, gl_238, \
                         gl_239, il_235, il_236, il_237, il_238, \
                         il_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -3.0 * gl_235[k]
                   + f_0 * il_235[k];

        t_236[k] = -3.0 * gl_236[k]
                   + f_0 * il_236[k];

        t_237[k] = -3.0 * gl_237[k]
                   + f_0 * il_237[k];

        t_238[k] = -3.0 * gl_238[k]
                   + f_0 * il_238[k];

        t_239[k] = -3.0 * gl_239[k]
                   + f_0 * il_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, gl_240, gl_241, gl_242, gl_243, \
                         gl_244, il_240, il_241, il_242, il_243, \
                         il_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * gl_240[k]
                   + f_0 * il_240[k];

        t_241[k] = -3.0 * gl_241[k]
                   + f_0 * il_241[k];

        t_242[k] = -3.0 * gl_242[k]
                   + f_0 * il_242[k];

        t_243[k] = -3.0 * gl_243[k]
                   + f_0 * il_243[k];

        t_244[k] = -3.0 * gl_244[k]
                   + f_0 * il_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, gl_245, gl_246, gl_247, gl_248, \
                         gl_249, il_245, il_246, il_247, il_248, \
                         il_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * gl_245[k]
                   + f_0 * il_245[k];

        t_246[k] = -3.0 * gl_246[k]
                   + f_0 * il_246[k];

        t_247[k] = -3.0 * gl_247[k]
                   + f_0 * il_247[k];

        t_248[k] = -3.0 * gl_248[k]
                   + f_0 * il_248[k];

        t_249[k] = -3.0 * gl_249[k]
                   + f_0 * il_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, gl_250, gl_251, gl_252, gl_253, \
                         gl_254, il_250, il_251, il_252, il_253, \
                         il_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -3.0 * gl_250[k]
                   + f_0 * il_250[k];

        t_251[k] = -3.0 * gl_251[k]
                   + f_0 * il_251[k];

        t_252[k] = -3.0 * gl_252[k]
                   + f_0 * il_252[k];

        t_253[k] = -3.0 * gl_253[k]
                   + f_0 * il_253[k];

        t_254[k] = -3.0 * gl_254[k]
                   + f_0 * il_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, gl_255, gl_256, gl_257, gl_258, \
                         gl_259, il_255, il_256, il_257, il_258, \
                         il_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -3.0 * gl_255[k]
                   + f_0 * il_255[k];

        t_256[k] = -3.0 * gl_256[k]
                   + f_0 * il_256[k];

        t_257[k] = -3.0 * gl_257[k]
                   + f_0 * il_257[k];

        t_258[k] = -3.0 * gl_258[k]
                   + f_0 * il_258[k];

        t_259[k] = -3.0 * gl_259[k]
                   + f_0 * il_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, gl_260, gl_261, gl_262, gl_263, \
                         gl_264, il_260, il_261, il_262, il_263, \
                         il_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -3.0 * gl_260[k]
                   + f_0 * il_260[k];

        t_261[k] = -3.0 * gl_261[k]
                   + f_0 * il_261[k];

        t_262[k] = -3.0 * gl_262[k]
                   + f_0 * il_262[k];

        t_263[k] = -3.0 * gl_263[k]
                   + f_0 * il_263[k];

        t_264[k] = -3.0 * gl_264[k]
                   + f_0 * il_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, gl_265, gl_266, gl_267, gl_268, \
                         gl_269, il_265, il_266, il_267, il_268, \
                         il_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -3.0 * gl_265[k]
                   + f_0 * il_265[k];

        t_266[k] = -3.0 * gl_266[k]
                   + f_0 * il_266[k];

        t_267[k] = -3.0 * gl_267[k]
                   + f_0 * il_267[k];

        t_268[k] = -3.0 * gl_268[k]
                   + f_0 * il_268[k];

        t_269[k] = -3.0 * gl_269[k]
                   + f_0 * il_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, gl_270, gl_271, gl_272, gl_273, \
                         gl_274, il_270, il_271, il_272, il_273, \
                         il_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * gl_270[k]
                   + f_0 * il_270[k];

        t_271[k] = -2.0 * gl_271[k]
                   + f_0 * il_271[k];

        t_272[k] = -2.0 * gl_272[k]
                   + f_0 * il_272[k];

        t_273[k] = -2.0 * gl_273[k]
                   + f_0 * il_273[k];

        t_274[k] = -2.0 * gl_274[k]
                   + f_0 * il_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, gl_275, gl_276, gl_277, gl_278, \
                         gl_279, il_275, il_276, il_277, il_278, \
                         il_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * gl_275[k]
                   + f_0 * il_275[k];

        t_276[k] = -2.0 * gl_276[k]
                   + f_0 * il_276[k];

        t_277[k] = -2.0 * gl_277[k]
                   + f_0 * il_277[k];

        t_278[k] = -2.0 * gl_278[k]
                   + f_0 * il_278[k];

        t_279[k] = -2.0 * gl_279[k]
                   + f_0 * il_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, gl_280, gl_281, gl_282, gl_283, \
                         gl_284, il_280, il_281, il_282, il_283, \
                         il_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * gl_280[k]
                   + f_0 * il_280[k];

        t_281[k] = -2.0 * gl_281[k]
                   + f_0 * il_281[k];

        t_282[k] = -2.0 * gl_282[k]
                   + f_0 * il_282[k];

        t_283[k] = -2.0 * gl_283[k]
                   + f_0 * il_283[k];

        t_284[k] = -2.0 * gl_284[k]
                   + f_0 * il_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, gl_285, gl_286, gl_287, gl_288, \
                         gl_289, il_285, il_286, il_287, il_288, \
                         il_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -2.0 * gl_285[k]
                   + f_0 * il_285[k];

        t_286[k] = -2.0 * gl_286[k]
                   + f_0 * il_286[k];

        t_287[k] = -2.0 * gl_287[k]
                   + f_0 * il_287[k];

        t_288[k] = -2.0 * gl_288[k]
                   + f_0 * il_288[k];

        t_289[k] = -2.0 * gl_289[k]
                   + f_0 * il_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, gl_290, gl_291, gl_292, gl_293, \
                         gl_294, il_290, il_291, il_292, il_293, \
                         il_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -2.0 * gl_290[k]
                   + f_0 * il_290[k];

        t_291[k] = -2.0 * gl_291[k]
                   + f_0 * il_291[k];

        t_292[k] = -2.0 * gl_292[k]
                   + f_0 * il_292[k];

        t_293[k] = -2.0 * gl_293[k]
                   + f_0 * il_293[k];

        t_294[k] = -2.0 * gl_294[k]
                   + f_0 * il_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, gl_295, gl_296, gl_297, gl_298, \
                         gl_299, il_295, il_296, il_297, il_298, \
                         il_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -2.0 * gl_295[k]
                   + f_0 * il_295[k];

        t_296[k] = -2.0 * gl_296[k]
                   + f_0 * il_296[k];

        t_297[k] = -2.0 * gl_297[k]
                   + f_0 * il_297[k];

        t_298[k] = -2.0 * gl_298[k]
                   + f_0 * il_298[k];

        t_299[k] = -2.0 * gl_299[k]
                   + f_0 * il_299[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_300 = buffer.data(gl + 300);
    const auto *gl_301 = buffer.data(gl + 301);
    const auto *gl_302 = buffer.data(gl + 302);
    const auto *gl_303 = buffer.data(gl + 303);
    const auto *gl_304 = buffer.data(gl + 304);
    const auto *gl_305 = buffer.data(gl + 305);
    const auto *gl_306 = buffer.data(gl + 306);
    const auto *gl_307 = buffer.data(gl + 307);
    const auto *gl_308 = buffer.data(gl + 308);
    const auto *gl_309 = buffer.data(gl + 309);
    const auto *gl_310 = buffer.data(gl + 310);
    const auto *gl_311 = buffer.data(gl + 311);
    const auto *gl_312 = buffer.data(gl + 312);
    const auto *gl_313 = buffer.data(gl + 313);
    const auto *gl_314 = buffer.data(gl + 314);
    const auto *gl_315 = buffer.data(gl + 315);
    const auto *gl_316 = buffer.data(gl + 316);
    const auto *gl_317 = buffer.data(gl + 317);
    const auto *gl_318 = buffer.data(gl + 318);
    const auto *gl_319 = buffer.data(gl + 319);
    const auto *gl_320 = buffer.data(gl + 320);
    const auto *gl_321 = buffer.data(gl + 321);
    const auto *gl_322 = buffer.data(gl + 322);
    const auto *gl_323 = buffer.data(gl + 323);
    const auto *gl_324 = buffer.data(gl + 324);
    const auto *gl_325 = buffer.data(gl + 325);
    const auto *gl_326 = buffer.data(gl + 326);
    const auto *gl_327 = buffer.data(gl + 327);
    const auto *gl_328 = buffer.data(gl + 328);
    const auto *gl_329 = buffer.data(gl + 329);
    const auto *gl_330 = buffer.data(gl + 330);
    const auto *gl_331 = buffer.data(gl + 331);
    const auto *gl_332 = buffer.data(gl + 332);
    const auto *gl_333 = buffer.data(gl + 333);
    const auto *gl_334 = buffer.data(gl + 334);
    const auto *gl_335 = buffer.data(gl + 335);
    const auto *gl_336 = buffer.data(gl + 336);
    const auto *gl_337 = buffer.data(gl + 337);
    const auto *gl_338 = buffer.data(gl + 338);
    const auto *gl_339 = buffer.data(gl + 339);
    const auto *gl_340 = buffer.data(gl + 340);
    const auto *gl_341 = buffer.data(gl + 341);
    const auto *gl_342 = buffer.data(gl + 342);
    const auto *gl_343 = buffer.data(gl + 343);
    const auto *gl_344 = buffer.data(gl + 344);
    const auto *gl_345 = buffer.data(gl + 345);
    const auto *gl_346 = buffer.data(gl + 346);
    const auto *gl_347 = buffer.data(gl + 347);
    const auto *gl_348 = buffer.data(gl + 348);
    const auto *gl_349 = buffer.data(gl + 349);
    const auto *gl_350 = buffer.data(gl + 350);
    const auto *gl_351 = buffer.data(gl + 351);
    const auto *gl_352 = buffer.data(gl + 352);
    const auto *gl_353 = buffer.data(gl + 353);
    const auto *gl_354 = buffer.data(gl + 354);
    const auto *gl_355 = buffer.data(gl + 355);
    const auto *gl_356 = buffer.data(gl + 356);
    const auto *gl_357 = buffer.data(gl + 357);
    const auto *gl_358 = buffer.data(gl + 358);
    const auto *gl_359 = buffer.data(gl + 359);
    const auto *gl_360 = buffer.data(gl + 360);
    const auto *gl_361 = buffer.data(gl + 361);
    const auto *gl_362 = buffer.data(gl + 362);
    const auto *gl_363 = buffer.data(gl + 363);
    const auto *gl_364 = buffer.data(gl + 364);
    const auto *gl_365 = buffer.data(gl + 365);
    const auto *gl_366 = buffer.data(gl + 366);
    const auto *gl_367 = buffer.data(gl + 367);
    const auto *gl_368 = buffer.data(gl + 368);
    const auto *gl_369 = buffer.data(gl + 369);
    const auto *gl_370 = buffer.data(gl + 370);
    const auto *gl_371 = buffer.data(gl + 371);
    const auto *gl_372 = buffer.data(gl + 372);
    const auto *gl_373 = buffer.data(gl + 373);
    const auto *gl_374 = buffer.data(gl + 374);
    const auto *gl_375 = buffer.data(gl + 375);
    const auto *gl_376 = buffer.data(gl + 376);
    const auto *gl_377 = buffer.data(gl + 377);
    const auto *gl_378 = buffer.data(gl + 378);
    const auto *gl_379 = buffer.data(gl + 379);
    const auto *gl_380 = buffer.data(gl + 380);
    const auto *gl_381 = buffer.data(gl + 381);
    const auto *gl_382 = buffer.data(gl + 382);
    const auto *gl_383 = buffer.data(gl + 383);
    const auto *gl_384 = buffer.data(gl + 384);
    const auto *gl_385 = buffer.data(gl + 385);
    const auto *gl_386 = buffer.data(gl + 386);
    const auto *gl_387 = buffer.data(gl + 387);
    const auto *gl_388 = buffer.data(gl + 388);
    const auto *gl_389 = buffer.data(gl + 389);
    const auto *gl_390 = buffer.data(gl + 390);
    const auto *gl_391 = buffer.data(gl + 391);
    const auto *gl_392 = buffer.data(gl + 392);
    const auto *gl_393 = buffer.data(gl + 393);
    const auto *gl_394 = buffer.data(gl + 394);
    const auto *gl_395 = buffer.data(gl + 395);
    const auto *gl_396 = buffer.data(gl + 396);
    const auto *gl_397 = buffer.data(gl + 397);
    const auto *gl_398 = buffer.data(gl + 398);
    const auto *gl_399 = buffer.data(gl + 399);
    const auto *gl_400 = buffer.data(gl + 400);
    const auto *gl_401 = buffer.data(gl + 401);
    const auto *gl_402 = buffer.data(gl + 402);
    const auto *gl_403 = buffer.data(gl + 403);
    const auto *gl_404 = buffer.data(gl + 404);
    const auto *gl_405 = buffer.data(gl + 405);
    const auto *gl_406 = buffer.data(gl + 406);
    const auto *gl_407 = buffer.data(gl + 407);
    const auto *gl_408 = buffer.data(gl + 408);
    const auto *gl_409 = buffer.data(gl + 409);
    const auto *gl_410 = buffer.data(gl + 410);
    const auto *gl_411 = buffer.data(gl + 411);
    const auto *gl_412 = buffer.data(gl + 412);
    const auto *gl_413 = buffer.data(gl + 413);
    const auto *gl_414 = buffer.data(gl + 414);
    const auto *gl_415 = buffer.data(gl + 415);
    const auto *gl_416 = buffer.data(gl + 416);
    const auto *gl_417 = buffer.data(gl + 417);
    const auto *gl_418 = buffer.data(gl + 418);
    const auto *gl_419 = buffer.data(gl + 419);
    const auto *gl_420 = buffer.data(gl + 420);
    const auto *gl_421 = buffer.data(gl + 421);
    const auto *gl_422 = buffer.data(gl + 422);
    const auto *gl_423 = buffer.data(gl + 423);
    const auto *gl_424 = buffer.data(gl + 424);
    const auto *gl_425 = buffer.data(gl + 425);
    const auto *gl_426 = buffer.data(gl + 426);
    const auto *gl_427 = buffer.data(gl + 427);
    const auto *gl_428 = buffer.data(gl + 428);
    const auto *gl_429 = buffer.data(gl + 429);
    const auto *gl_430 = buffer.data(gl + 430);
    const auto *gl_431 = buffer.data(gl + 431);
    const auto *gl_432 = buffer.data(gl + 432);
    const auto *gl_433 = buffer.data(gl + 433);
    const auto *gl_434 = buffer.data(gl + 434);
    const auto *gl_435 = buffer.data(gl + 435);
    const auto *gl_436 = buffer.data(gl + 436);
    const auto *gl_437 = buffer.data(gl + 437);
    const auto *gl_438 = buffer.data(gl + 438);
    const auto *gl_439 = buffer.data(gl + 439);
    const auto *gl_440 = buffer.data(gl + 440);
    const auto *gl_441 = buffer.data(gl + 441);
    const auto *gl_442 = buffer.data(gl + 442);
    const auto *gl_443 = buffer.data(gl + 443);
    const auto *gl_444 = buffer.data(gl + 444);
    const auto *gl_445 = buffer.data(gl + 445);
    const auto *gl_446 = buffer.data(gl + 446);
    const auto *gl_447 = buffer.data(gl + 447);
    const auto *gl_448 = buffer.data(gl + 448);
    const auto *gl_449 = buffer.data(gl + 449);

    const auto *il_300 = buffer.data(il + 300);
    const auto *il_301 = buffer.data(il + 301);
    const auto *il_302 = buffer.data(il + 302);
    const auto *il_303 = buffer.data(il + 303);
    const auto *il_304 = buffer.data(il + 304);
    const auto *il_305 = buffer.data(il + 305);
    const auto *il_306 = buffer.data(il + 306);
    const auto *il_307 = buffer.data(il + 307);
    const auto *il_308 = buffer.data(il + 308);
    const auto *il_309 = buffer.data(il + 309);
    const auto *il_310 = buffer.data(il + 310);
    const auto *il_311 = buffer.data(il + 311);
    const auto *il_312 = buffer.data(il + 312);
    const auto *il_313 = buffer.data(il + 313);
    const auto *il_314 = buffer.data(il + 314);
    const auto *il_315 = buffer.data(il + 315);
    const auto *il_316 = buffer.data(il + 316);
    const auto *il_317 = buffer.data(il + 317);
    const auto *il_318 = buffer.data(il + 318);
    const auto *il_319 = buffer.data(il + 319);
    const auto *il_320 = buffer.data(il + 320);
    const auto *il_321 = buffer.data(il + 321);
    const auto *il_322 = buffer.data(il + 322);
    const auto *il_323 = buffer.data(il + 323);
    const auto *il_324 = buffer.data(il + 324);
    const auto *il_325 = buffer.data(il + 325);
    const auto *il_326 = buffer.data(il + 326);
    const auto *il_327 = buffer.data(il + 327);
    const auto *il_328 = buffer.data(il + 328);
    const auto *il_329 = buffer.data(il + 329);
    const auto *il_330 = buffer.data(il + 330);
    const auto *il_331 = buffer.data(il + 331);
    const auto *il_332 = buffer.data(il + 332);
    const auto *il_333 = buffer.data(il + 333);
    const auto *il_334 = buffer.data(il + 334);
    const auto *il_335 = buffer.data(il + 335);
    const auto *il_336 = buffer.data(il + 336);
    const auto *il_337 = buffer.data(il + 337);
    const auto *il_338 = buffer.data(il + 338);
    const auto *il_339 = buffer.data(il + 339);
    const auto *il_340 = buffer.data(il + 340);
    const auto *il_341 = buffer.data(il + 341);
    const auto *il_342 = buffer.data(il + 342);
    const auto *il_343 = buffer.data(il + 343);
    const auto *il_344 = buffer.data(il + 344);
    const auto *il_345 = buffer.data(il + 345);
    const auto *il_346 = buffer.data(il + 346);
    const auto *il_347 = buffer.data(il + 347);
    const auto *il_348 = buffer.data(il + 348);
    const auto *il_349 = buffer.data(il + 349);
    const auto *il_350 = buffer.data(il + 350);
    const auto *il_351 = buffer.data(il + 351);
    const auto *il_352 = buffer.data(il + 352);
    const auto *il_353 = buffer.data(il + 353);
    const auto *il_354 = buffer.data(il + 354);
    const auto *il_355 = buffer.data(il + 355);
    const auto *il_356 = buffer.data(il + 356);
    const auto *il_357 = buffer.data(il + 357);
    const auto *il_358 = buffer.data(il + 358);
    const auto *il_359 = buffer.data(il + 359);
    const auto *il_360 = buffer.data(il + 360);
    const auto *il_361 = buffer.data(il + 361);
    const auto *il_362 = buffer.data(il + 362);
    const auto *il_363 = buffer.data(il + 363);
    const auto *il_364 = buffer.data(il + 364);
    const auto *il_365 = buffer.data(il + 365);
    const auto *il_366 = buffer.data(il + 366);
    const auto *il_367 = buffer.data(il + 367);
    const auto *il_368 = buffer.data(il + 368);
    const auto *il_369 = buffer.data(il + 369);
    const auto *il_370 = buffer.data(il + 370);
    const auto *il_371 = buffer.data(il + 371);
    const auto *il_372 = buffer.data(il + 372);
    const auto *il_373 = buffer.data(il + 373);
    const auto *il_374 = buffer.data(il + 374);
    const auto *il_375 = buffer.data(il + 375);
    const auto *il_376 = buffer.data(il + 376);
    const auto *il_377 = buffer.data(il + 377);
    const auto *il_378 = buffer.data(il + 378);
    const auto *il_379 = buffer.data(il + 379);
    const auto *il_380 = buffer.data(il + 380);
    const auto *il_381 = buffer.data(il + 381);
    const auto *il_382 = buffer.data(il + 382);
    const auto *il_383 = buffer.data(il + 383);
    const auto *il_384 = buffer.data(il + 384);
    const auto *il_385 = buffer.data(il + 385);
    const auto *il_386 = buffer.data(il + 386);
    const auto *il_387 = buffer.data(il + 387);
    const auto *il_388 = buffer.data(il + 388);
    const auto *il_389 = buffer.data(il + 389);
    const auto *il_390 = buffer.data(il + 390);
    const auto *il_391 = buffer.data(il + 391);
    const auto *il_392 = buffer.data(il + 392);
    const auto *il_393 = buffer.data(il + 393);
    const auto *il_394 = buffer.data(il + 394);
    const auto *il_395 = buffer.data(il + 395);
    const auto *il_396 = buffer.data(il + 396);
    const auto *il_397 = buffer.data(il + 397);
    const auto *il_398 = buffer.data(il + 398);
    const auto *il_399 = buffer.data(il + 399);
    const auto *il_400 = buffer.data(il + 400);
    const auto *il_401 = buffer.data(il + 401);
    const auto *il_402 = buffer.data(il + 402);
    const auto *il_403 = buffer.data(il + 403);
    const auto *il_404 = buffer.data(il + 404);
    const auto *il_405 = buffer.data(il + 405);
    const auto *il_406 = buffer.data(il + 406);
    const auto *il_407 = buffer.data(il + 407);
    const auto *il_408 = buffer.data(il + 408);
    const auto *il_409 = buffer.data(il + 409);
    const auto *il_410 = buffer.data(il + 410);
    const auto *il_411 = buffer.data(il + 411);
    const auto *il_412 = buffer.data(il + 412);
    const auto *il_413 = buffer.data(il + 413);
    const auto *il_414 = buffer.data(il + 414);
    const auto *il_415 = buffer.data(il + 415);
    const auto *il_416 = buffer.data(il + 416);
    const auto *il_417 = buffer.data(il + 417);
    const auto *il_418 = buffer.data(il + 418);
    const auto *il_419 = buffer.data(il + 419);
    const auto *il_420 = buffer.data(il + 420);
    const auto *il_421 = buffer.data(il + 421);
    const auto *il_422 = buffer.data(il + 422);
    const auto *il_423 = buffer.data(il + 423);
    const auto *il_424 = buffer.data(il + 424);
    const auto *il_425 = buffer.data(il + 425);
    const auto *il_426 = buffer.data(il + 426);
    const auto *il_427 = buffer.data(il + 427);
    const auto *il_428 = buffer.data(il + 428);
    const auto *il_429 = buffer.data(il + 429);
    const auto *il_430 = buffer.data(il + 430);
    const auto *il_431 = buffer.data(il + 431);
    const auto *il_432 = buffer.data(il + 432);
    const auto *il_433 = buffer.data(il + 433);
    const auto *il_434 = buffer.data(il + 434);
    const auto *il_435 = buffer.data(il + 435);
    const auto *il_436 = buffer.data(il + 436);
    const auto *il_437 = buffer.data(il + 437);
    const auto *il_438 = buffer.data(il + 438);
    const auto *il_439 = buffer.data(il + 439);
    const auto *il_440 = buffer.data(il + 440);
    const auto *il_441 = buffer.data(il + 441);
    const auto *il_442 = buffer.data(il + 442);
    const auto *il_443 = buffer.data(il + 443);
    const auto *il_444 = buffer.data(il + 444);
    const auto *il_445 = buffer.data(il + 445);
    const auto *il_446 = buffer.data(il + 446);
    const auto *il_447 = buffer.data(il + 447);
    const auto *il_448 = buffer.data(il + 448);
    const auto *il_449 = buffer.data(il + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, gl_300, gl_301, gl_302, gl_303, \
                         gl_304, il_300, il_301, il_302, il_303, \
                         il_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -2.0 * gl_300[k]
                   + f_0 * il_300[k];

        t_301[k] = -2.0 * gl_301[k]
                   + f_0 * il_301[k];

        t_302[k] = -2.0 * gl_302[k]
                   + f_0 * il_302[k];

        t_303[k] = -2.0 * gl_303[k]
                   + f_0 * il_303[k];

        t_304[k] = -2.0 * gl_304[k]
                   + f_0 * il_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, gl_305, gl_306, gl_307, gl_308, \
                         gl_309, il_305, il_306, il_307, il_308, \
                         il_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -2.0 * gl_305[k]
                   + f_0 * il_305[k];

        t_306[k] = -2.0 * gl_306[k]
                   + f_0 * il_306[k];

        t_307[k] = -2.0 * gl_307[k]
                   + f_0 * il_307[k];

        t_308[k] = -2.0 * gl_308[k]
                   + f_0 * il_308[k];

        t_309[k] = -2.0 * gl_309[k]
                   + f_0 * il_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, gl_310, gl_311, gl_312, gl_313, \
                         gl_314, il_310, il_311, il_312, il_313, \
                         il_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -2.0 * gl_310[k]
                   + f_0 * il_310[k];

        t_311[k] = -2.0 * gl_311[k]
                   + f_0 * il_311[k];

        t_312[k] = -2.0 * gl_312[k]
                   + f_0 * il_312[k];

        t_313[k] = -2.0 * gl_313[k]
                   + f_0 * il_313[k];

        t_314[k] = -2.0 * gl_314[k]
                   + f_0 * il_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, gl_315, gl_316, gl_317, gl_318, \
                         gl_319, il_315, il_316, il_317, il_318, \
                         il_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -2.0 * gl_315[k]
                   + f_0 * il_315[k];

        t_316[k] = -2.0 * gl_316[k]
                   + f_0 * il_316[k];

        t_317[k] = -2.0 * gl_317[k]
                   + f_0 * il_317[k];

        t_318[k] = -2.0 * gl_318[k]
                   + f_0 * il_318[k];

        t_319[k] = -2.0 * gl_319[k]
                   + f_0 * il_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, gl_320, gl_321, gl_322, gl_323, \
                         gl_324, il_320, il_321, il_322, il_323, \
                         il_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -2.0 * gl_320[k]
                   + f_0 * il_320[k];

        t_321[k] = -2.0 * gl_321[k]
                   + f_0 * il_321[k];

        t_322[k] = -2.0 * gl_322[k]
                   + f_0 * il_322[k];

        t_323[k] = -2.0 * gl_323[k]
                   + f_0 * il_323[k];

        t_324[k] = -2.0 * gl_324[k]
                   + f_0 * il_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, gl_325, gl_326, gl_327, gl_328, \
                         gl_329, il_325, il_326, il_327, il_328, \
                         il_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -2.0 * gl_325[k]
                   + f_0 * il_325[k];

        t_326[k] = -2.0 * gl_326[k]
                   + f_0 * il_326[k];

        t_327[k] = -2.0 * gl_327[k]
                   + f_0 * il_327[k];

        t_328[k] = -2.0 * gl_328[k]
                   + f_0 * il_328[k];

        t_329[k] = -2.0 * gl_329[k]
                   + f_0 * il_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, gl_330, gl_331, gl_332, gl_333, \
                         gl_334, il_330, il_331, il_332, il_333, \
                         il_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -2.0 * gl_330[k]
                   + f_0 * il_330[k];

        t_331[k] = -2.0 * gl_331[k]
                   + f_0 * il_331[k];

        t_332[k] = -2.0 * gl_332[k]
                   + f_0 * il_332[k];

        t_333[k] = -2.0 * gl_333[k]
                   + f_0 * il_333[k];

        t_334[k] = -2.0 * gl_334[k]
                   + f_0 * il_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, gl_335, gl_336, gl_337, gl_338, \
                         gl_339, il_335, il_336, il_337, il_338, \
                         il_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -2.0 * gl_335[k]
                   + f_0 * il_335[k];

        t_336[k] = -2.0 * gl_336[k]
                   + f_0 * il_336[k];

        t_337[k] = -2.0 * gl_337[k]
                   + f_0 * il_337[k];

        t_338[k] = -2.0 * gl_338[k]
                   + f_0 * il_338[k];

        t_339[k] = -2.0 * gl_339[k]
                   + f_0 * il_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, gl_340, gl_341, gl_342, gl_343, \
                         gl_344, il_340, il_341, il_342, il_343, \
                         il_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * gl_340[k]
                   + f_0 * il_340[k];

        t_341[k] = -2.0 * gl_341[k]
                   + f_0 * il_341[k];

        t_342[k] = -2.0 * gl_342[k]
                   + f_0 * il_342[k];

        t_343[k] = -2.0 * gl_343[k]
                   + f_0 * il_343[k];

        t_344[k] = -2.0 * gl_344[k]
                   + f_0 * il_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, gl_345, gl_346, gl_347, gl_348, \
                         gl_349, il_345, il_346, il_347, il_348, \
                         il_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * gl_345[k]
                   + f_0 * il_345[k];

        t_346[k] = -2.0 * gl_346[k]
                   + f_0 * il_346[k];

        t_347[k] = -2.0 * gl_347[k]
                   + f_0 * il_347[k];

        t_348[k] = -2.0 * gl_348[k]
                   + f_0 * il_348[k];

        t_349[k] = -2.0 * gl_349[k]
                   + f_0 * il_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, gl_350, gl_351, gl_352, gl_353, \
                         gl_354, il_350, il_351, il_352, il_353, \
                         il_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * gl_350[k]
                   + f_0 * il_350[k];

        t_351[k] = -2.0 * gl_351[k]
                   + f_0 * il_351[k];

        t_352[k] = -2.0 * gl_352[k]
                   + f_0 * il_352[k];

        t_353[k] = -2.0 * gl_353[k]
                   + f_0 * il_353[k];

        t_354[k] = -2.0 * gl_354[k]
                   + f_0 * il_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, gl_355, gl_356, gl_357, gl_358, \
                         gl_359, il_355, il_356, il_357, il_358, \
                         il_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * gl_355[k]
                   + f_0 * il_355[k];

        t_356[k] = -2.0 * gl_356[k]
                   + f_0 * il_356[k];

        t_357[k] = -2.0 * gl_357[k]
                   + f_0 * il_357[k];

        t_358[k] = -2.0 * gl_358[k]
                   + f_0 * il_358[k];

        t_359[k] = -2.0 * gl_359[k]
                   + f_0 * il_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, gl_360, gl_361, gl_362, gl_363, \
                         gl_364, il_360, il_361, il_362, il_363, \
                         il_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -2.0 * gl_360[k]
                   + f_0 * il_360[k];

        t_361[k] = -2.0 * gl_361[k]
                   + f_0 * il_361[k];

        t_362[k] = -2.0 * gl_362[k]
                   + f_0 * il_362[k];

        t_363[k] = -2.0 * gl_363[k]
                   + f_0 * il_363[k];

        t_364[k] = -2.0 * gl_364[k]
                   + f_0 * il_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, gl_365, gl_366, gl_367, gl_368, \
                         gl_369, il_365, il_366, il_367, il_368, \
                         il_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -2.0 * gl_365[k]
                   + f_0 * il_365[k];

        t_366[k] = -2.0 * gl_366[k]
                   + f_0 * il_366[k];

        t_367[k] = -2.0 * gl_367[k]
                   + f_0 * il_367[k];

        t_368[k] = -2.0 * gl_368[k]
                   + f_0 * il_368[k];

        t_369[k] = -2.0 * gl_369[k]
                   + f_0 * il_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, gl_370, gl_371, gl_372, gl_373, \
                         gl_374, il_370, il_371, il_372, il_373, \
                         il_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -2.0 * gl_370[k]
                   + f_0 * il_370[k];

        t_371[k] = -2.0 * gl_371[k]
                   + f_0 * il_371[k];

        t_372[k] = -2.0 * gl_372[k]
                   + f_0 * il_372[k];

        t_373[k] = -2.0 * gl_373[k]
                   + f_0 * il_373[k];

        t_374[k] = -2.0 * gl_374[k]
                   + f_0 * il_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, gl_375, gl_376, gl_377, gl_378, \
                         gl_379, il_375, il_376, il_377, il_378, \
                         il_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -2.0 * gl_375[k]
                   + f_0 * il_375[k];

        t_376[k] = -2.0 * gl_376[k]
                   + f_0 * il_376[k];

        t_377[k] = -2.0 * gl_377[k]
                   + f_0 * il_377[k];

        t_378[k] = -2.0 * gl_378[k]
                   + f_0 * il_378[k];

        t_379[k] = -2.0 * gl_379[k]
                   + f_0 * il_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, gl_380, gl_381, gl_382, gl_383, \
                         gl_384, il_380, il_381, il_382, il_383, \
                         il_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * gl_380[k]
                   + f_0 * il_380[k];

        t_381[k] = -2.0 * gl_381[k]
                   + f_0 * il_381[k];

        t_382[k] = -2.0 * gl_382[k]
                   + f_0 * il_382[k];

        t_383[k] = -2.0 * gl_383[k]
                   + f_0 * il_383[k];

        t_384[k] = -2.0 * gl_384[k]
                   + f_0 * il_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, gl_385, gl_386, gl_387, gl_388, \
                         gl_389, il_385, il_386, il_387, il_388, \
                         il_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * gl_385[k]
                   + f_0 * il_385[k];

        t_386[k] = -2.0 * gl_386[k]
                   + f_0 * il_386[k];

        t_387[k] = -2.0 * gl_387[k]
                   + f_0 * il_387[k];

        t_388[k] = -2.0 * gl_388[k]
                   + f_0 * il_388[k];

        t_389[k] = -2.0 * gl_389[k]
                   + f_0 * il_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, gl_390, gl_391, gl_392, gl_393, \
                         gl_394, il_390, il_391, il_392, il_393, \
                         il_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -2.0 * gl_390[k]
                   + f_0 * il_390[k];

        t_391[k] = -2.0 * gl_391[k]
                   + f_0 * il_391[k];

        t_392[k] = -2.0 * gl_392[k]
                   + f_0 * il_392[k];

        t_393[k] = -2.0 * gl_393[k]
                   + f_0 * il_393[k];

        t_394[k] = -2.0 * gl_394[k]
                   + f_0 * il_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, gl_395, gl_396, gl_397, gl_398, \
                         gl_399, il_395, il_396, il_397, il_398, \
                         il_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -2.0 * gl_395[k]
                   + f_0 * il_395[k];

        t_396[k] = -2.0 * gl_396[k]
                   + f_0 * il_396[k];

        t_397[k] = -2.0 * gl_397[k]
                   + f_0 * il_397[k];

        t_398[k] = -2.0 * gl_398[k]
                   + f_0 * il_398[k];

        t_399[k] = -2.0 * gl_399[k]
                   + f_0 * il_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, gl_400, gl_401, gl_402, gl_403, \
                         gl_404, il_400, il_401, il_402, il_403, \
                         il_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -2.0 * gl_400[k]
                   + f_0 * il_400[k];

        t_401[k] = -2.0 * gl_401[k]
                   + f_0 * il_401[k];

        t_402[k] = -2.0 * gl_402[k]
                   + f_0 * il_402[k];

        t_403[k] = -2.0 * gl_403[k]
                   + f_0 * il_403[k];

        t_404[k] = -2.0 * gl_404[k]
                   + f_0 * il_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, gl_405, gl_406, gl_407, gl_408, \
                         gl_409, il_405, il_406, il_407, il_408, \
                         il_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -2.0 * gl_405[k]
                   + f_0 * il_405[k];

        t_406[k] = -2.0 * gl_406[k]
                   + f_0 * il_406[k];

        t_407[k] = -2.0 * gl_407[k]
                   + f_0 * il_407[k];

        t_408[k] = -2.0 * gl_408[k]
                   + f_0 * il_408[k];

        t_409[k] = -2.0 * gl_409[k]
                   + f_0 * il_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, gl_410, gl_411, gl_412, gl_413, \
                         gl_414, il_410, il_411, il_412, il_413, \
                         il_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -2.0 * gl_410[k]
                   + f_0 * il_410[k];

        t_411[k] = -2.0 * gl_411[k]
                   + f_0 * il_411[k];

        t_412[k] = -2.0 * gl_412[k]
                   + f_0 * il_412[k];

        t_413[k] = -2.0 * gl_413[k]
                   + f_0 * il_413[k];

        t_414[k] = -2.0 * gl_414[k]
                   + f_0 * il_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, gl_415, gl_416, gl_417, gl_418, \
                         gl_419, il_415, il_416, il_417, il_418, \
                         il_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -2.0 * gl_415[k]
                   + f_0 * il_415[k];

        t_416[k] = -2.0 * gl_416[k]
                   + f_0 * il_416[k];

        t_417[k] = -2.0 * gl_417[k]
                   + f_0 * il_417[k];

        t_418[k] = -2.0 * gl_418[k]
                   + f_0 * il_418[k];

        t_419[k] = -2.0 * gl_419[k]
                   + f_0 * il_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, gl_420, gl_421, gl_422, gl_423, \
                         gl_424, il_420, il_421, il_422, il_423, \
                         il_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -2.0 * gl_420[k]
                   + f_0 * il_420[k];

        t_421[k] = -2.0 * gl_421[k]
                   + f_0 * il_421[k];

        t_422[k] = -2.0 * gl_422[k]
                   + f_0 * il_422[k];

        t_423[k] = -2.0 * gl_423[k]
                   + f_0 * il_423[k];

        t_424[k] = -2.0 * gl_424[k]
                   + f_0 * il_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, gl_425, gl_426, gl_427, gl_428, \
                         gl_429, il_425, il_426, il_427, il_428, \
                         il_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -2.0 * gl_425[k]
                   + f_0 * il_425[k];

        t_426[k] = -2.0 * gl_426[k]
                   + f_0 * il_426[k];

        t_427[k] = -2.0 * gl_427[k]
                   + f_0 * il_427[k];

        t_428[k] = -2.0 * gl_428[k]
                   + f_0 * il_428[k];

        t_429[k] = -2.0 * gl_429[k]
                   + f_0 * il_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, gl_430, gl_431, gl_432, gl_433, \
                         gl_434, il_430, il_431, il_432, il_433, \
                         il_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -2.0 * gl_430[k]
                   + f_0 * il_430[k];

        t_431[k] = -2.0 * gl_431[k]
                   + f_0 * il_431[k];

        t_432[k] = -2.0 * gl_432[k]
                   + f_0 * il_432[k];

        t_433[k] = -2.0 * gl_433[k]
                   + f_0 * il_433[k];

        t_434[k] = -2.0 * gl_434[k]
                   + f_0 * il_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, gl_435, gl_436, gl_437, gl_438, \
                         gl_439, il_435, il_436, il_437, il_438, \
                         il_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -2.0 * gl_435[k]
                   + f_0 * il_435[k];

        t_436[k] = -2.0 * gl_436[k]
                   + f_0 * il_436[k];

        t_437[k] = -2.0 * gl_437[k]
                   + f_0 * il_437[k];

        t_438[k] = -2.0 * gl_438[k]
                   + f_0 * il_438[k];

        t_439[k] = -2.0 * gl_439[k]
                   + f_0 * il_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, gl_440, gl_441, gl_442, gl_443, \
                         gl_444, il_440, il_441, il_442, il_443, \
                         il_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -2.0 * gl_440[k]
                   + f_0 * il_440[k];

        t_441[k] = -2.0 * gl_441[k]
                   + f_0 * il_441[k];

        t_442[k] = -2.0 * gl_442[k]
                   + f_0 * il_442[k];

        t_443[k] = -2.0 * gl_443[k]
                   + f_0 * il_443[k];

        t_444[k] = -2.0 * gl_444[k]
                   + f_0 * il_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, gl_445, gl_446, gl_447, gl_448, \
                         gl_449, il_445, il_446, il_447, il_448, \
                         il_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -2.0 * gl_445[k]
                   + f_0 * il_445[k];

        t_446[k] = -2.0 * gl_446[k]
                   + f_0 * il_446[k];

        t_447[k] = -2.0 * gl_447[k]
                   + f_0 * il_447[k];

        t_448[k] = -2.0 * gl_448[k]
                   + f_0 * il_448[k];

        t_449[k] = -2.0 * gl_449[k]
                   + f_0 * il_449[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_450 = buffer.data(gl + 450);
    const auto *gl_451 = buffer.data(gl + 451);
    const auto *gl_452 = buffer.data(gl + 452);
    const auto *gl_453 = buffer.data(gl + 453);
    const auto *gl_454 = buffer.data(gl + 454);
    const auto *gl_455 = buffer.data(gl + 455);
    const auto *gl_456 = buffer.data(gl + 456);
    const auto *gl_457 = buffer.data(gl + 457);
    const auto *gl_458 = buffer.data(gl + 458);
    const auto *gl_459 = buffer.data(gl + 459);
    const auto *gl_460 = buffer.data(gl + 460);
    const auto *gl_461 = buffer.data(gl + 461);
    const auto *gl_462 = buffer.data(gl + 462);
    const auto *gl_463 = buffer.data(gl + 463);
    const auto *gl_464 = buffer.data(gl + 464);
    const auto *gl_465 = buffer.data(gl + 465);
    const auto *gl_466 = buffer.data(gl + 466);
    const auto *gl_467 = buffer.data(gl + 467);
    const auto *gl_468 = buffer.data(gl + 468);
    const auto *gl_469 = buffer.data(gl + 469);
    const auto *gl_470 = buffer.data(gl + 470);
    const auto *gl_471 = buffer.data(gl + 471);
    const auto *gl_472 = buffer.data(gl + 472);
    const auto *gl_473 = buffer.data(gl + 473);
    const auto *gl_474 = buffer.data(gl + 474);
    const auto *gl_475 = buffer.data(gl + 475);
    const auto *gl_476 = buffer.data(gl + 476);
    const auto *gl_477 = buffer.data(gl + 477);
    const auto *gl_478 = buffer.data(gl + 478);
    const auto *gl_479 = buffer.data(gl + 479);
    const auto *gl_480 = buffer.data(gl + 480);
    const auto *gl_481 = buffer.data(gl + 481);
    const auto *gl_482 = buffer.data(gl + 482);
    const auto *gl_483 = buffer.data(gl + 483);
    const auto *gl_484 = buffer.data(gl + 484);
    const auto *gl_485 = buffer.data(gl + 485);
    const auto *gl_486 = buffer.data(gl + 486);
    const auto *gl_487 = buffer.data(gl + 487);
    const auto *gl_488 = buffer.data(gl + 488);
    const auto *gl_489 = buffer.data(gl + 489);
    const auto *gl_490 = buffer.data(gl + 490);
    const auto *gl_491 = buffer.data(gl + 491);
    const auto *gl_492 = buffer.data(gl + 492);
    const auto *gl_493 = buffer.data(gl + 493);
    const auto *gl_494 = buffer.data(gl + 494);
    const auto *gl_495 = buffer.data(gl + 495);
    const auto *gl_496 = buffer.data(gl + 496);
    const auto *gl_497 = buffer.data(gl + 497);
    const auto *gl_498 = buffer.data(gl + 498);
    const auto *gl_499 = buffer.data(gl + 499);
    const auto *gl_500 = buffer.data(gl + 500);
    const auto *gl_501 = buffer.data(gl + 501);
    const auto *gl_502 = buffer.data(gl + 502);
    const auto *gl_503 = buffer.data(gl + 503);
    const auto *gl_504 = buffer.data(gl + 504);
    const auto *gl_505 = buffer.data(gl + 505);
    const auto *gl_506 = buffer.data(gl + 506);
    const auto *gl_507 = buffer.data(gl + 507);
    const auto *gl_508 = buffer.data(gl + 508);
    const auto *gl_509 = buffer.data(gl + 509);
    const auto *gl_510 = buffer.data(gl + 510);
    const auto *gl_511 = buffer.data(gl + 511);
    const auto *gl_512 = buffer.data(gl + 512);
    const auto *gl_513 = buffer.data(gl + 513);
    const auto *gl_514 = buffer.data(gl + 514);
    const auto *gl_515 = buffer.data(gl + 515);
    const auto *gl_516 = buffer.data(gl + 516);
    const auto *gl_517 = buffer.data(gl + 517);
    const auto *gl_518 = buffer.data(gl + 518);
    const auto *gl_519 = buffer.data(gl + 519);
    const auto *gl_520 = buffer.data(gl + 520);
    const auto *gl_521 = buffer.data(gl + 521);
    const auto *gl_522 = buffer.data(gl + 522);
    const auto *gl_523 = buffer.data(gl + 523);
    const auto *gl_524 = buffer.data(gl + 524);
    const auto *gl_525 = buffer.data(gl + 525);
    const auto *gl_526 = buffer.data(gl + 526);
    const auto *gl_527 = buffer.data(gl + 527);
    const auto *gl_528 = buffer.data(gl + 528);
    const auto *gl_529 = buffer.data(gl + 529);
    const auto *gl_530 = buffer.data(gl + 530);
    const auto *gl_531 = buffer.data(gl + 531);
    const auto *gl_532 = buffer.data(gl + 532);
    const auto *gl_533 = buffer.data(gl + 533);
    const auto *gl_534 = buffer.data(gl + 534);
    const auto *gl_535 = buffer.data(gl + 535);
    const auto *gl_536 = buffer.data(gl + 536);
    const auto *gl_537 = buffer.data(gl + 537);
    const auto *gl_538 = buffer.data(gl + 538);
    const auto *gl_539 = buffer.data(gl + 539);
    const auto *gl_540 = buffer.data(gl + 540);
    const auto *gl_541 = buffer.data(gl + 541);
    const auto *gl_542 = buffer.data(gl + 542);
    const auto *gl_543 = buffer.data(gl + 543);
    const auto *gl_544 = buffer.data(gl + 544);
    const auto *gl_545 = buffer.data(gl + 545);
    const auto *gl_546 = buffer.data(gl + 546);
    const auto *gl_547 = buffer.data(gl + 547);
    const auto *gl_548 = buffer.data(gl + 548);
    const auto *gl_549 = buffer.data(gl + 549);
    const auto *gl_550 = buffer.data(gl + 550);
    const auto *gl_551 = buffer.data(gl + 551);
    const auto *gl_552 = buffer.data(gl + 552);
    const auto *gl_553 = buffer.data(gl + 553);
    const auto *gl_554 = buffer.data(gl + 554);
    const auto *gl_555 = buffer.data(gl + 555);
    const auto *gl_556 = buffer.data(gl + 556);
    const auto *gl_557 = buffer.data(gl + 557);
    const auto *gl_558 = buffer.data(gl + 558);
    const auto *gl_559 = buffer.data(gl + 559);
    const auto *gl_560 = buffer.data(gl + 560);
    const auto *gl_561 = buffer.data(gl + 561);
    const auto *gl_562 = buffer.data(gl + 562);
    const auto *gl_563 = buffer.data(gl + 563);
    const auto *gl_564 = buffer.data(gl + 564);
    const auto *gl_565 = buffer.data(gl + 565);
    const auto *gl_566 = buffer.data(gl + 566);
    const auto *gl_567 = buffer.data(gl + 567);
    const auto *gl_568 = buffer.data(gl + 568);
    const auto *gl_569 = buffer.data(gl + 569);
    const auto *gl_570 = buffer.data(gl + 570);
    const auto *gl_571 = buffer.data(gl + 571);
    const auto *gl_572 = buffer.data(gl + 572);
    const auto *gl_573 = buffer.data(gl + 573);
    const auto *gl_574 = buffer.data(gl + 574);
    const auto *gl_575 = buffer.data(gl + 575);
    const auto *gl_576 = buffer.data(gl + 576);
    const auto *gl_577 = buffer.data(gl + 577);
    const auto *gl_578 = buffer.data(gl + 578);
    const auto *gl_579 = buffer.data(gl + 579);
    const auto *gl_580 = buffer.data(gl + 580);
    const auto *gl_581 = buffer.data(gl + 581);
    const auto *gl_582 = buffer.data(gl + 582);
    const auto *gl_583 = buffer.data(gl + 583);
    const auto *gl_584 = buffer.data(gl + 584);
    const auto *gl_585 = buffer.data(gl + 585);
    const auto *gl_586 = buffer.data(gl + 586);
    const auto *gl_587 = buffer.data(gl + 587);
    const auto *gl_588 = buffer.data(gl + 588);
    const auto *gl_589 = buffer.data(gl + 589);
    const auto *gl_590 = buffer.data(gl + 590);
    const auto *gl_591 = buffer.data(gl + 591);
    const auto *gl_592 = buffer.data(gl + 592);
    const auto *gl_593 = buffer.data(gl + 593);
    const auto *gl_594 = buffer.data(gl + 594);
    const auto *gl_595 = buffer.data(gl + 595);
    const auto *gl_596 = buffer.data(gl + 596);
    const auto *gl_597 = buffer.data(gl + 597);
    const auto *gl_598 = buffer.data(gl + 598);
    const auto *gl_599 = buffer.data(gl + 599);

    const auto *il_450 = buffer.data(il + 450);
    const auto *il_451 = buffer.data(il + 451);
    const auto *il_452 = buffer.data(il + 452);
    const auto *il_453 = buffer.data(il + 453);
    const auto *il_454 = buffer.data(il + 454);
    const auto *il_455 = buffer.data(il + 455);
    const auto *il_456 = buffer.data(il + 456);
    const auto *il_457 = buffer.data(il + 457);
    const auto *il_458 = buffer.data(il + 458);
    const auto *il_459 = buffer.data(il + 459);
    const auto *il_460 = buffer.data(il + 460);
    const auto *il_461 = buffer.data(il + 461);
    const auto *il_462 = buffer.data(il + 462);
    const auto *il_463 = buffer.data(il + 463);
    const auto *il_464 = buffer.data(il + 464);
    const auto *il_465 = buffer.data(il + 465);
    const auto *il_466 = buffer.data(il + 466);
    const auto *il_467 = buffer.data(il + 467);
    const auto *il_468 = buffer.data(il + 468);
    const auto *il_469 = buffer.data(il + 469);
    const auto *il_470 = buffer.data(il + 470);
    const auto *il_471 = buffer.data(il + 471);
    const auto *il_472 = buffer.data(il + 472);
    const auto *il_473 = buffer.data(il + 473);
    const auto *il_474 = buffer.data(il + 474);
    const auto *il_475 = buffer.data(il + 475);
    const auto *il_476 = buffer.data(il + 476);
    const auto *il_477 = buffer.data(il + 477);
    const auto *il_478 = buffer.data(il + 478);
    const auto *il_479 = buffer.data(il + 479);
    const auto *il_480 = buffer.data(il + 480);
    const auto *il_481 = buffer.data(il + 481);
    const auto *il_482 = buffer.data(il + 482);
    const auto *il_483 = buffer.data(il + 483);
    const auto *il_484 = buffer.data(il + 484);
    const auto *il_485 = buffer.data(il + 485);
    const auto *il_486 = buffer.data(il + 486);
    const auto *il_487 = buffer.data(il + 487);
    const auto *il_488 = buffer.data(il + 488);
    const auto *il_489 = buffer.data(il + 489);
    const auto *il_490 = buffer.data(il + 490);
    const auto *il_491 = buffer.data(il + 491);
    const auto *il_492 = buffer.data(il + 492);
    const auto *il_493 = buffer.data(il + 493);
    const auto *il_494 = buffer.data(il + 494);
    const auto *il_495 = buffer.data(il + 495);
    const auto *il_496 = buffer.data(il + 496);
    const auto *il_497 = buffer.data(il + 497);
    const auto *il_498 = buffer.data(il + 498);
    const auto *il_499 = buffer.data(il + 499);
    const auto *il_500 = buffer.data(il + 500);
    const auto *il_501 = buffer.data(il + 501);
    const auto *il_502 = buffer.data(il + 502);
    const auto *il_503 = buffer.data(il + 503);
    const auto *il_504 = buffer.data(il + 504);
    const auto *il_505 = buffer.data(il + 505);
    const auto *il_506 = buffer.data(il + 506);
    const auto *il_507 = buffer.data(il + 507);
    const auto *il_508 = buffer.data(il + 508);
    const auto *il_509 = buffer.data(il + 509);
    const auto *il_510 = buffer.data(il + 510);
    const auto *il_511 = buffer.data(il + 511);
    const auto *il_512 = buffer.data(il + 512);
    const auto *il_513 = buffer.data(il + 513);
    const auto *il_514 = buffer.data(il + 514);
    const auto *il_515 = buffer.data(il + 515);
    const auto *il_516 = buffer.data(il + 516);
    const auto *il_517 = buffer.data(il + 517);
    const auto *il_518 = buffer.data(il + 518);
    const auto *il_519 = buffer.data(il + 519);
    const auto *il_520 = buffer.data(il + 520);
    const auto *il_521 = buffer.data(il + 521);
    const auto *il_522 = buffer.data(il + 522);
    const auto *il_523 = buffer.data(il + 523);
    const auto *il_524 = buffer.data(il + 524);
    const auto *il_525 = buffer.data(il + 525);
    const auto *il_526 = buffer.data(il + 526);
    const auto *il_527 = buffer.data(il + 527);
    const auto *il_528 = buffer.data(il + 528);
    const auto *il_529 = buffer.data(il + 529);
    const auto *il_530 = buffer.data(il + 530);
    const auto *il_531 = buffer.data(il + 531);
    const auto *il_532 = buffer.data(il + 532);
    const auto *il_533 = buffer.data(il + 533);
    const auto *il_534 = buffer.data(il + 534);
    const auto *il_535 = buffer.data(il + 535);
    const auto *il_536 = buffer.data(il + 536);
    const auto *il_537 = buffer.data(il + 537);
    const auto *il_538 = buffer.data(il + 538);
    const auto *il_539 = buffer.data(il + 539);
    const auto *il_540 = buffer.data(il + 540);
    const auto *il_541 = buffer.data(il + 541);
    const auto *il_542 = buffer.data(il + 542);
    const auto *il_543 = buffer.data(il + 543);
    const auto *il_544 = buffer.data(il + 544);
    const auto *il_545 = buffer.data(il + 545);
    const auto *il_546 = buffer.data(il + 546);
    const auto *il_547 = buffer.data(il + 547);
    const auto *il_548 = buffer.data(il + 548);
    const auto *il_549 = buffer.data(il + 549);
    const auto *il_550 = buffer.data(il + 550);
    const auto *il_551 = buffer.data(il + 551);
    const auto *il_552 = buffer.data(il + 552);
    const auto *il_553 = buffer.data(il + 553);
    const auto *il_554 = buffer.data(il + 554);
    const auto *il_555 = buffer.data(il + 555);
    const auto *il_556 = buffer.data(il + 556);
    const auto *il_557 = buffer.data(il + 557);
    const auto *il_558 = buffer.data(il + 558);
    const auto *il_559 = buffer.data(il + 559);
    const auto *il_560 = buffer.data(il + 560);
    const auto *il_561 = buffer.data(il + 561);
    const auto *il_562 = buffer.data(il + 562);
    const auto *il_563 = buffer.data(il + 563);
    const auto *il_564 = buffer.data(il + 564);
    const auto *il_565 = buffer.data(il + 565);
    const auto *il_566 = buffer.data(il + 566);
    const auto *il_567 = buffer.data(il + 567);
    const auto *il_568 = buffer.data(il + 568);
    const auto *il_569 = buffer.data(il + 569);
    const auto *il_570 = buffer.data(il + 570);
    const auto *il_571 = buffer.data(il + 571);
    const auto *il_572 = buffer.data(il + 572);
    const auto *il_573 = buffer.data(il + 573);
    const auto *il_574 = buffer.data(il + 574);
    const auto *il_575 = buffer.data(il + 575);
    const auto *il_576 = buffer.data(il + 576);
    const auto *il_577 = buffer.data(il + 577);
    const auto *il_578 = buffer.data(il + 578);
    const auto *il_579 = buffer.data(il + 579);
    const auto *il_580 = buffer.data(il + 580);
    const auto *il_581 = buffer.data(il + 581);
    const auto *il_582 = buffer.data(il + 582);
    const auto *il_583 = buffer.data(il + 583);
    const auto *il_584 = buffer.data(il + 584);
    const auto *il_585 = buffer.data(il + 585);
    const auto *il_586 = buffer.data(il + 586);
    const auto *il_587 = buffer.data(il + 587);
    const auto *il_588 = buffer.data(il + 588);
    const auto *il_589 = buffer.data(il + 589);
    const auto *il_590 = buffer.data(il + 590);
    const auto *il_591 = buffer.data(il + 591);
    const auto *il_592 = buffer.data(il + 592);
    const auto *il_593 = buffer.data(il + 593);
    const auto *il_594 = buffer.data(il + 594);
    const auto *il_595 = buffer.data(il + 595);
    const auto *il_596 = buffer.data(il + 596);
    const auto *il_597 = buffer.data(il + 597);
    const auto *il_598 = buffer.data(il + 598);
    const auto *il_599 = buffer.data(il + 599);

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, gl_450, gl_451, gl_452, gl_453, \
                         gl_454, il_450, il_451, il_452, il_453, \
                         il_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -gl_450[k]
                   + f_0 * il_450[k];

        t_451[k] = -gl_451[k]
                   + f_0 * il_451[k];

        t_452[k] = -gl_452[k]
                   + f_0 * il_452[k];

        t_453[k] = -gl_453[k]
                   + f_0 * il_453[k];

        t_454[k] = -gl_454[k]
                   + f_0 * il_454[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, gl_455, gl_456, gl_457, gl_458, \
                         gl_459, il_455, il_456, il_457, il_458, \
                         il_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -gl_455[k]
                   + f_0 * il_455[k];

        t_456[k] = -gl_456[k]
                   + f_0 * il_456[k];

        t_457[k] = -gl_457[k]
                   + f_0 * il_457[k];

        t_458[k] = -gl_458[k]
                   + f_0 * il_458[k];

        t_459[k] = -gl_459[k]
                   + f_0 * il_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, gl_460, gl_461, gl_462, gl_463, \
                         gl_464, il_460, il_461, il_462, il_463, \
                         il_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -gl_460[k]
                   + f_0 * il_460[k];

        t_461[k] = -gl_461[k]
                   + f_0 * il_461[k];

        t_462[k] = -gl_462[k]
                   + f_0 * il_462[k];

        t_463[k] = -gl_463[k]
                   + f_0 * il_463[k];

        t_464[k] = -gl_464[k]
                   + f_0 * il_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, gl_465, gl_466, gl_467, gl_468, \
                         gl_469, il_465, il_466, il_467, il_468, \
                         il_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -gl_465[k]
                   + f_0 * il_465[k];

        t_466[k] = -gl_466[k]
                   + f_0 * il_466[k];

        t_467[k] = -gl_467[k]
                   + f_0 * il_467[k];

        t_468[k] = -gl_468[k]
                   + f_0 * il_468[k];

        t_469[k] = -gl_469[k]
                   + f_0 * il_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, gl_470, gl_471, gl_472, gl_473, \
                         gl_474, il_470, il_471, il_472, il_473, \
                         il_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -gl_470[k]
                   + f_0 * il_470[k];

        t_471[k] = -gl_471[k]
                   + f_0 * il_471[k];

        t_472[k] = -gl_472[k]
                   + f_0 * il_472[k];

        t_473[k] = -gl_473[k]
                   + f_0 * il_473[k];

        t_474[k] = -gl_474[k]
                   + f_0 * il_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, gl_475, gl_476, gl_477, gl_478, \
                         gl_479, il_475, il_476, il_477, il_478, \
                         il_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -gl_475[k]
                   + f_0 * il_475[k];

        t_476[k] = -gl_476[k]
                   + f_0 * il_476[k];

        t_477[k] = -gl_477[k]
                   + f_0 * il_477[k];

        t_478[k] = -gl_478[k]
                   + f_0 * il_478[k];

        t_479[k] = -gl_479[k]
                   + f_0 * il_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, gl_480, gl_481, gl_482, gl_483, \
                         gl_484, il_480, il_481, il_482, il_483, \
                         il_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -gl_480[k]
                   + f_0 * il_480[k];

        t_481[k] = -gl_481[k]
                   + f_0 * il_481[k];

        t_482[k] = -gl_482[k]
                   + f_0 * il_482[k];

        t_483[k] = -gl_483[k]
                   + f_0 * il_483[k];

        t_484[k] = -gl_484[k]
                   + f_0 * il_484[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, gl_485, gl_486, gl_487, gl_488, \
                         gl_489, il_485, il_486, il_487, il_488, \
                         il_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -gl_485[k]
                   + f_0 * il_485[k];

        t_486[k] = -gl_486[k]
                   + f_0 * il_486[k];

        t_487[k] = -gl_487[k]
                   + f_0 * il_487[k];

        t_488[k] = -gl_488[k]
                   + f_0 * il_488[k];

        t_489[k] = -gl_489[k]
                   + f_0 * il_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, gl_490, gl_491, gl_492, gl_493, \
                         gl_494, il_490, il_491, il_492, il_493, \
                         il_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -gl_490[k]
                   + f_0 * il_490[k];

        t_491[k] = -gl_491[k]
                   + f_0 * il_491[k];

        t_492[k] = -gl_492[k]
                   + f_0 * il_492[k];

        t_493[k] = -gl_493[k]
                   + f_0 * il_493[k];

        t_494[k] = -gl_494[k]
                   + f_0 * il_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, gl_495, gl_496, gl_497, gl_498, \
                         gl_499, il_495, il_496, il_497, il_498, \
                         il_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -gl_495[k]
                   + f_0 * il_495[k];

        t_496[k] = -gl_496[k]
                   + f_0 * il_496[k];

        t_497[k] = -gl_497[k]
                   + f_0 * il_497[k];

        t_498[k] = -gl_498[k]
                   + f_0 * il_498[k];

        t_499[k] = -gl_499[k]
                   + f_0 * il_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, gl_500, gl_501, gl_502, gl_503, \
                         gl_504, il_500, il_501, il_502, il_503, \
                         il_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -gl_500[k]
                   + f_0 * il_500[k];

        t_501[k] = -gl_501[k]
                   + f_0 * il_501[k];

        t_502[k] = -gl_502[k]
                   + f_0 * il_502[k];

        t_503[k] = -gl_503[k]
                   + f_0 * il_503[k];

        t_504[k] = -gl_504[k]
                   + f_0 * il_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, gl_505, gl_506, gl_507, gl_508, \
                         gl_509, il_505, il_506, il_507, il_508, \
                         il_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -gl_505[k]
                   + f_0 * il_505[k];

        t_506[k] = -gl_506[k]
                   + f_0 * il_506[k];

        t_507[k] = -gl_507[k]
                   + f_0 * il_507[k];

        t_508[k] = -gl_508[k]
                   + f_0 * il_508[k];

        t_509[k] = -gl_509[k]
                   + f_0 * il_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, gl_510, gl_511, gl_512, gl_513, \
                         gl_514, il_510, il_511, il_512, il_513, \
                         il_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -gl_510[k]
                   + f_0 * il_510[k];

        t_511[k] = -gl_511[k]
                   + f_0 * il_511[k];

        t_512[k] = -gl_512[k]
                   + f_0 * il_512[k];

        t_513[k] = -gl_513[k]
                   + f_0 * il_513[k];

        t_514[k] = -gl_514[k]
                   + f_0 * il_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, gl_515, gl_516, gl_517, gl_518, \
                         gl_519, il_515, il_516, il_517, il_518, \
                         il_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -gl_515[k]
                   + f_0 * il_515[k];

        t_516[k] = -gl_516[k]
                   + f_0 * il_516[k];

        t_517[k] = -gl_517[k]
                   + f_0 * il_517[k];

        t_518[k] = -gl_518[k]
                   + f_0 * il_518[k];

        t_519[k] = -gl_519[k]
                   + f_0 * il_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, gl_520, gl_521, gl_522, gl_523, \
                         gl_524, il_520, il_521, il_522, il_523, \
                         il_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -gl_520[k]
                   + f_0 * il_520[k];

        t_521[k] = -gl_521[k]
                   + f_0 * il_521[k];

        t_522[k] = -gl_522[k]
                   + f_0 * il_522[k];

        t_523[k] = -gl_523[k]
                   + f_0 * il_523[k];

        t_524[k] = -gl_524[k]
                   + f_0 * il_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, gl_525, gl_526, gl_527, gl_528, \
                         gl_529, il_525, il_526, il_527, il_528, \
                         il_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -gl_525[k]
                   + f_0 * il_525[k];

        t_526[k] = -gl_526[k]
                   + f_0 * il_526[k];

        t_527[k] = -gl_527[k]
                   + f_0 * il_527[k];

        t_528[k] = -gl_528[k]
                   + f_0 * il_528[k];

        t_529[k] = -gl_529[k]
                   + f_0 * il_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, gl_530, gl_531, gl_532, gl_533, \
                         gl_534, il_530, il_531, il_532, il_533, \
                         il_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -gl_530[k]
                   + f_0 * il_530[k];

        t_531[k] = -gl_531[k]
                   + f_0 * il_531[k];

        t_532[k] = -gl_532[k]
                   + f_0 * il_532[k];

        t_533[k] = -gl_533[k]
                   + f_0 * il_533[k];

        t_534[k] = -gl_534[k]
                   + f_0 * il_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, gl_535, gl_536, gl_537, gl_538, \
                         gl_539, il_535, il_536, il_537, il_538, \
                         il_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -gl_535[k]
                   + f_0 * il_535[k];

        t_536[k] = -gl_536[k]
                   + f_0 * il_536[k];

        t_537[k] = -gl_537[k]
                   + f_0 * il_537[k];

        t_538[k] = -gl_538[k]
                   + f_0 * il_538[k];

        t_539[k] = -gl_539[k]
                   + f_0 * il_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, gl_540, gl_541, gl_542, gl_543, \
                         gl_544, il_540, il_541, il_542, il_543, \
                         il_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -gl_540[k]
                   + f_0 * il_540[k];

        t_541[k] = -gl_541[k]
                   + f_0 * il_541[k];

        t_542[k] = -gl_542[k]
                   + f_0 * il_542[k];

        t_543[k] = -gl_543[k]
                   + f_0 * il_543[k];

        t_544[k] = -gl_544[k]
                   + f_0 * il_544[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, gl_545, gl_546, gl_547, gl_548, \
                         gl_549, il_545, il_546, il_547, il_548, \
                         il_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -gl_545[k]
                   + f_0 * il_545[k];

        t_546[k] = -gl_546[k]
                   + f_0 * il_546[k];

        t_547[k] = -gl_547[k]
                   + f_0 * il_547[k];

        t_548[k] = -gl_548[k]
                   + f_0 * il_548[k];

        t_549[k] = -gl_549[k]
                   + f_0 * il_549[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, gl_550, gl_551, gl_552, gl_553, \
                         gl_554, il_550, il_551, il_552, il_553, \
                         il_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -gl_550[k]
                   + f_0 * il_550[k];

        t_551[k] = -gl_551[k]
                   + f_0 * il_551[k];

        t_552[k] = -gl_552[k]
                   + f_0 * il_552[k];

        t_553[k] = -gl_553[k]
                   + f_0 * il_553[k];

        t_554[k] = -gl_554[k]
                   + f_0 * il_554[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, gl_555, gl_556, gl_557, gl_558, \
                         gl_559, il_555, il_556, il_557, il_558, \
                         il_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -gl_555[k]
                   + f_0 * il_555[k];

        t_556[k] = -gl_556[k]
                   + f_0 * il_556[k];

        t_557[k] = -gl_557[k]
                   + f_0 * il_557[k];

        t_558[k] = -gl_558[k]
                   + f_0 * il_558[k];

        t_559[k] = -gl_559[k]
                   + f_0 * il_559[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, gl_560, gl_561, gl_562, gl_563, \
                         gl_564, il_560, il_561, il_562, il_563, \
                         il_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = -gl_560[k]
                   + f_0 * il_560[k];

        t_561[k] = -gl_561[k]
                   + f_0 * il_561[k];

        t_562[k] = -gl_562[k]
                   + f_0 * il_562[k];

        t_563[k] = -gl_563[k]
                   + f_0 * il_563[k];

        t_564[k] = -gl_564[k]
                   + f_0 * il_564[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, gl_565, gl_566, gl_567, gl_568, \
                         gl_569, il_565, il_566, il_567, il_568, \
                         il_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = -gl_565[k]
                   + f_0 * il_565[k];

        t_566[k] = -gl_566[k]
                   + f_0 * il_566[k];

        t_567[k] = -gl_567[k]
                   + f_0 * il_567[k];

        t_568[k] = -gl_568[k]
                   + f_0 * il_568[k];

        t_569[k] = -gl_569[k]
                   + f_0 * il_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, gl_570, gl_571, gl_572, gl_573, \
                         gl_574, il_570, il_571, il_572, il_573, \
                         il_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -gl_570[k]
                   + f_0 * il_570[k];

        t_571[k] = -gl_571[k]
                   + f_0 * il_571[k];

        t_572[k] = -gl_572[k]
                   + f_0 * il_572[k];

        t_573[k] = -gl_573[k]
                   + f_0 * il_573[k];

        t_574[k] = -gl_574[k]
                   + f_0 * il_574[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, gl_575, gl_576, gl_577, gl_578, \
                         gl_579, il_575, il_576, il_577, il_578, \
                         il_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = -gl_575[k]
                   + f_0 * il_575[k];

        t_576[k] = -gl_576[k]
                   + f_0 * il_576[k];

        t_577[k] = -gl_577[k]
                   + f_0 * il_577[k];

        t_578[k] = -gl_578[k]
                   + f_0 * il_578[k];

        t_579[k] = -gl_579[k]
                   + f_0 * il_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, gl_580, gl_581, gl_582, gl_583, \
                         gl_584, il_580, il_581, il_582, il_583, \
                         il_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -gl_580[k]
                   + f_0 * il_580[k];

        t_581[k] = -gl_581[k]
                   + f_0 * il_581[k];

        t_582[k] = -gl_582[k]
                   + f_0 * il_582[k];

        t_583[k] = -gl_583[k]
                   + f_0 * il_583[k];

        t_584[k] = -gl_584[k]
                   + f_0 * il_584[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, gl_585, gl_586, gl_587, gl_588, \
                         gl_589, il_585, il_586, il_587, il_588, \
                         il_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -gl_585[k]
                   + f_0 * il_585[k];

        t_586[k] = -gl_586[k]
                   + f_0 * il_586[k];

        t_587[k] = -gl_587[k]
                   + f_0 * il_587[k];

        t_588[k] = -gl_588[k]
                   + f_0 * il_588[k];

        t_589[k] = -gl_589[k]
                   + f_0 * il_589[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, gl_590, gl_591, gl_592, gl_593, \
                         gl_594, il_590, il_591, il_592, il_593, \
                         il_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -gl_590[k]
                   + f_0 * il_590[k];

        t_591[k] = -gl_591[k]
                   + f_0 * il_591[k];

        t_592[k] = -gl_592[k]
                   + f_0 * il_592[k];

        t_593[k] = -gl_593[k]
                   + f_0 * il_593[k];

        t_594[k] = -gl_594[k]
                   + f_0 * il_594[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, gl_595, gl_596, gl_597, gl_598, \
                         gl_599, il_595, il_596, il_597, il_598, \
                         il_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -gl_595[k]
                   + f_0 * il_595[k];

        t_596[k] = -gl_596[k]
                   + f_0 * il_596[k];

        t_597[k] = -gl_597[k]
                   + f_0 * il_597[k];

        t_598[k] = -gl_598[k]
                   + f_0 * il_598[k];

        t_599[k] = -gl_599[k]
                   + f_0 * il_599[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_600 = buffer.data(gl + 600);
    const auto *gl_601 = buffer.data(gl + 601);
    const auto *gl_602 = buffer.data(gl + 602);
    const auto *gl_603 = buffer.data(gl + 603);
    const auto *gl_604 = buffer.data(gl + 604);
    const auto *gl_605 = buffer.data(gl + 605);
    const auto *gl_606 = buffer.data(gl + 606);
    const auto *gl_607 = buffer.data(gl + 607);
    const auto *gl_608 = buffer.data(gl + 608);
    const auto *gl_609 = buffer.data(gl + 609);
    const auto *gl_610 = buffer.data(gl + 610);
    const auto *gl_611 = buffer.data(gl + 611);
    const auto *gl_612 = buffer.data(gl + 612);
    const auto *gl_613 = buffer.data(gl + 613);
    const auto *gl_614 = buffer.data(gl + 614);
    const auto *gl_615 = buffer.data(gl + 615);
    const auto *gl_616 = buffer.data(gl + 616);
    const auto *gl_617 = buffer.data(gl + 617);
    const auto *gl_618 = buffer.data(gl + 618);
    const auto *gl_619 = buffer.data(gl + 619);
    const auto *gl_620 = buffer.data(gl + 620);
    const auto *gl_621 = buffer.data(gl + 621);
    const auto *gl_622 = buffer.data(gl + 622);
    const auto *gl_623 = buffer.data(gl + 623);
    const auto *gl_624 = buffer.data(gl + 624);
    const auto *gl_625 = buffer.data(gl + 625);
    const auto *gl_626 = buffer.data(gl + 626);
    const auto *gl_627 = buffer.data(gl + 627);
    const auto *gl_628 = buffer.data(gl + 628);
    const auto *gl_629 = buffer.data(gl + 629);
    const auto *gl_630 = buffer.data(gl + 630);
    const auto *gl_631 = buffer.data(gl + 631);
    const auto *gl_632 = buffer.data(gl + 632);
    const auto *gl_633 = buffer.data(gl + 633);
    const auto *gl_634 = buffer.data(gl + 634);
    const auto *gl_635 = buffer.data(gl + 635);
    const auto *gl_636 = buffer.data(gl + 636);
    const auto *gl_637 = buffer.data(gl + 637);
    const auto *gl_638 = buffer.data(gl + 638);
    const auto *gl_639 = buffer.data(gl + 639);
    const auto *gl_640 = buffer.data(gl + 640);
    const auto *gl_641 = buffer.data(gl + 641);
    const auto *gl_642 = buffer.data(gl + 642);
    const auto *gl_643 = buffer.data(gl + 643);
    const auto *gl_644 = buffer.data(gl + 644);
    const auto *gl_645 = buffer.data(gl + 645);
    const auto *gl_646 = buffer.data(gl + 646);
    const auto *gl_647 = buffer.data(gl + 647);
    const auto *gl_648 = buffer.data(gl + 648);
    const auto *gl_649 = buffer.data(gl + 649);
    const auto *gl_650 = buffer.data(gl + 650);
    const auto *gl_651 = buffer.data(gl + 651);
    const auto *gl_652 = buffer.data(gl + 652);
    const auto *gl_653 = buffer.data(gl + 653);
    const auto *gl_654 = buffer.data(gl + 654);
    const auto *gl_655 = buffer.data(gl + 655);
    const auto *gl_656 = buffer.data(gl + 656);
    const auto *gl_657 = buffer.data(gl + 657);
    const auto *gl_658 = buffer.data(gl + 658);
    const auto *gl_659 = buffer.data(gl + 659);
    const auto *gl_660 = buffer.data(gl + 660);
    const auto *gl_661 = buffer.data(gl + 661);
    const auto *gl_662 = buffer.data(gl + 662);
    const auto *gl_663 = buffer.data(gl + 663);
    const auto *gl_664 = buffer.data(gl + 664);
    const auto *gl_665 = buffer.data(gl + 665);
    const auto *gl_666 = buffer.data(gl + 666);
    const auto *gl_667 = buffer.data(gl + 667);
    const auto *gl_668 = buffer.data(gl + 668);
    const auto *gl_669 = buffer.data(gl + 669);
    const auto *gl_670 = buffer.data(gl + 670);
    const auto *gl_671 = buffer.data(gl + 671);
    const auto *gl_672 = buffer.data(gl + 672);
    const auto *gl_673 = buffer.data(gl + 673);
    const auto *gl_674 = buffer.data(gl + 674);

    const auto *il_600 = buffer.data(il + 600);
    const auto *il_601 = buffer.data(il + 601);
    const auto *il_602 = buffer.data(il + 602);
    const auto *il_603 = buffer.data(il + 603);
    const auto *il_604 = buffer.data(il + 604);
    const auto *il_605 = buffer.data(il + 605);
    const auto *il_606 = buffer.data(il + 606);
    const auto *il_607 = buffer.data(il + 607);
    const auto *il_608 = buffer.data(il + 608);
    const auto *il_609 = buffer.data(il + 609);
    const auto *il_610 = buffer.data(il + 610);
    const auto *il_611 = buffer.data(il + 611);
    const auto *il_612 = buffer.data(il + 612);
    const auto *il_613 = buffer.data(il + 613);
    const auto *il_614 = buffer.data(il + 614);
    const auto *il_615 = buffer.data(il + 615);
    const auto *il_616 = buffer.data(il + 616);
    const auto *il_617 = buffer.data(il + 617);
    const auto *il_618 = buffer.data(il + 618);
    const auto *il_619 = buffer.data(il + 619);
    const auto *il_620 = buffer.data(il + 620);
    const auto *il_621 = buffer.data(il + 621);
    const auto *il_622 = buffer.data(il + 622);
    const auto *il_623 = buffer.data(il + 623);
    const auto *il_624 = buffer.data(il + 624);
    const auto *il_625 = buffer.data(il + 625);
    const auto *il_626 = buffer.data(il + 626);
    const auto *il_627 = buffer.data(il + 627);
    const auto *il_628 = buffer.data(il + 628);
    const auto *il_629 = buffer.data(il + 629);
    const auto *il_630 = buffer.data(il + 630);
    const auto *il_631 = buffer.data(il + 631);
    const auto *il_632 = buffer.data(il + 632);
    const auto *il_633 = buffer.data(il + 633);
    const auto *il_634 = buffer.data(il + 634);
    const auto *il_635 = buffer.data(il + 635);
    const auto *il_636 = buffer.data(il + 636);
    const auto *il_637 = buffer.data(il + 637);
    const auto *il_638 = buffer.data(il + 638);
    const auto *il_639 = buffer.data(il + 639);
    const auto *il_640 = buffer.data(il + 640);
    const auto *il_641 = buffer.data(il + 641);
    const auto *il_642 = buffer.data(il + 642);
    const auto *il_643 = buffer.data(il + 643);
    const auto *il_644 = buffer.data(il + 644);
    const auto *il_645 = buffer.data(il + 645);
    const auto *il_646 = buffer.data(il + 646);
    const auto *il_647 = buffer.data(il + 647);
    const auto *il_648 = buffer.data(il + 648);
    const auto *il_649 = buffer.data(il + 649);
    const auto *il_650 = buffer.data(il + 650);
    const auto *il_651 = buffer.data(il + 651);
    const auto *il_652 = buffer.data(il + 652);
    const auto *il_653 = buffer.data(il + 653);
    const auto *il_654 = buffer.data(il + 654);
    const auto *il_655 = buffer.data(il + 655);
    const auto *il_656 = buffer.data(il + 656);
    const auto *il_657 = buffer.data(il + 657);
    const auto *il_658 = buffer.data(il + 658);
    const auto *il_659 = buffer.data(il + 659);
    const auto *il_660 = buffer.data(il + 660);
    const auto *il_661 = buffer.data(il + 661);
    const auto *il_662 = buffer.data(il + 662);
    const auto *il_663 = buffer.data(il + 663);
    const auto *il_664 = buffer.data(il + 664);
    const auto *il_665 = buffer.data(il + 665);
    const auto *il_666 = buffer.data(il + 666);
    const auto *il_667 = buffer.data(il + 667);
    const auto *il_668 = buffer.data(il + 668);
    const auto *il_669 = buffer.data(il + 669);
    const auto *il_670 = buffer.data(il + 670);
    const auto *il_671 = buffer.data(il + 671);
    const auto *il_672 = buffer.data(il + 672);
    const auto *il_673 = buffer.data(il + 673);
    const auto *il_674 = buffer.data(il + 674);
    const auto *il_675 = buffer.data(il + 675);
    const auto *il_676 = buffer.data(il + 676);
    const auto *il_677 = buffer.data(il + 677);
    const auto *il_678 = buffer.data(il + 678);
    const auto *il_679 = buffer.data(il + 679);
    const auto *il_680 = buffer.data(il + 680);
    const auto *il_681 = buffer.data(il + 681);
    const auto *il_682 = buffer.data(il + 682);
    const auto *il_683 = buffer.data(il + 683);
    const auto *il_684 = buffer.data(il + 684);
    const auto *il_685 = buffer.data(il + 685);
    const auto *il_686 = buffer.data(il + 686);
    const auto *il_687 = buffer.data(il + 687);
    const auto *il_688 = buffer.data(il + 688);
    const auto *il_689 = buffer.data(il + 689);
    const auto *il_690 = buffer.data(il + 690);
    const auto *il_691 = buffer.data(il + 691);
    const auto *il_692 = buffer.data(il + 692);
    const auto *il_693 = buffer.data(il + 693);
    const auto *il_694 = buffer.data(il + 694);
    const auto *il_695 = buffer.data(il + 695);
    const auto *il_696 = buffer.data(il + 696);
    const auto *il_697 = buffer.data(il + 697);
    const auto *il_698 = buffer.data(il + 698);
    const auto *il_699 = buffer.data(il + 699);
    const auto *il_700 = buffer.data(il + 700);
    const auto *il_701 = buffer.data(il + 701);
    const auto *il_702 = buffer.data(il + 702);
    const auto *il_703 = buffer.data(il + 703);
    const auto *il_704 = buffer.data(il + 704);
    const auto *il_705 = buffer.data(il + 705);
    const auto *il_706 = buffer.data(il + 706);
    const auto *il_707 = buffer.data(il + 707);
    const auto *il_708 = buffer.data(il + 708);
    const auto *il_709 = buffer.data(il + 709);
    const auto *il_710 = buffer.data(il + 710);
    const auto *il_711 = buffer.data(il + 711);
    const auto *il_712 = buffer.data(il + 712);
    const auto *il_713 = buffer.data(il + 713);
    const auto *il_714 = buffer.data(il + 714);
    const auto *il_715 = buffer.data(il + 715);
    const auto *il_716 = buffer.data(il + 716);
    const auto *il_717 = buffer.data(il + 717);
    const auto *il_718 = buffer.data(il + 718);
    const auto *il_719 = buffer.data(il + 719);
    const auto *il_720 = buffer.data(il + 720);
    const auto *il_721 = buffer.data(il + 721);
    const auto *il_722 = buffer.data(il + 722);
    const auto *il_723 = buffer.data(il + 723);
    const auto *il_724 = buffer.data(il + 724);
    const auto *il_725 = buffer.data(il + 725);
    const auto *il_726 = buffer.data(il + 726);
    const auto *il_727 = buffer.data(il + 727);
    const auto *il_728 = buffer.data(il + 728);
    const auto *il_729 = buffer.data(il + 729);
    const auto *il_730 = buffer.data(il + 730);
    const auto *il_731 = buffer.data(il + 731);
    const auto *il_732 = buffer.data(il + 732);
    const auto *il_733 = buffer.data(il + 733);
    const auto *il_734 = buffer.data(il + 734);
    const auto *il_735 = buffer.data(il + 735);
    const auto *il_736 = buffer.data(il + 736);
    const auto *il_737 = buffer.data(il + 737);
    const auto *il_738 = buffer.data(il + 738);
    const auto *il_739 = buffer.data(il + 739);
    const auto *il_740 = buffer.data(il + 740);
    const auto *il_741 = buffer.data(il + 741);
    const auto *il_742 = buffer.data(il + 742);
    const auto *il_743 = buffer.data(il + 743);
    const auto *il_744 = buffer.data(il + 744);
    const auto *il_745 = buffer.data(il + 745);
    const auto *il_746 = buffer.data(il + 746);
    const auto *il_747 = buffer.data(il + 747);
    const auto *il_748 = buffer.data(il + 748);
    const auto *il_749 = buffer.data(il + 749);
    const auto *il_750 = buffer.data(il + 750);
    const auto *il_751 = buffer.data(il + 751);
    const auto *il_752 = buffer.data(il + 752);
    const auto *il_753 = buffer.data(il + 753);
    const auto *il_754 = buffer.data(il + 754);
    const auto *il_755 = buffer.data(il + 755);
    const auto *il_756 = buffer.data(il + 756);
    const auto *il_757 = buffer.data(il + 757);
    const auto *il_758 = buffer.data(il + 758);
    const auto *il_759 = buffer.data(il + 759);
    const auto *il_760 = buffer.data(il + 760);
    const auto *il_761 = buffer.data(il + 761);
    const auto *il_762 = buffer.data(il + 762);
    const auto *il_763 = buffer.data(il + 763);
    const auto *il_764 = buffer.data(il + 764);
    const auto *il_765 = buffer.data(il + 765);
    const auto *il_766 = buffer.data(il + 766);
    const auto *il_767 = buffer.data(il + 767);
    const auto *il_768 = buffer.data(il + 768);
    const auto *il_769 = buffer.data(il + 769);
    const auto *il_770 = buffer.data(il + 770);
    const auto *il_771 = buffer.data(il + 771);
    const auto *il_772 = buffer.data(il + 772);
    const auto *il_773 = buffer.data(il + 773);
    const auto *il_774 = buffer.data(il + 774);
    const auto *il_775 = buffer.data(il + 775);
    const auto *il_776 = buffer.data(il + 776);
    const auto *il_777 = buffer.data(il + 777);
    const auto *il_778 = buffer.data(il + 778);
    const auto *il_779 = buffer.data(il + 779);
    const auto *il_780 = buffer.data(il + 780);
    const auto *il_781 = buffer.data(il + 781);
    const auto *il_782 = buffer.data(il + 782);
    const auto *il_783 = buffer.data(il + 783);
    const auto *il_784 = buffer.data(il + 784);
    const auto *il_785 = buffer.data(il + 785);
    const auto *il_786 = buffer.data(il + 786);

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, gl_600, gl_601, gl_602, gl_603, \
                         gl_604, il_600, il_601, il_602, il_603, \
                         il_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -gl_600[k]
                   + f_0 * il_600[k];

        t_601[k] = -gl_601[k]
                   + f_0 * il_601[k];

        t_602[k] = -gl_602[k]
                   + f_0 * il_602[k];

        t_603[k] = -gl_603[k]
                   + f_0 * il_603[k];

        t_604[k] = -gl_604[k]
                   + f_0 * il_604[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, gl_605, gl_606, gl_607, gl_608, \
                         gl_609, il_605, il_606, il_607, il_608, \
                         il_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -gl_605[k]
                   + f_0 * il_605[k];

        t_606[k] = -gl_606[k]
                   + f_0 * il_606[k];

        t_607[k] = -gl_607[k]
                   + f_0 * il_607[k];

        t_608[k] = -gl_608[k]
                   + f_0 * il_608[k];

        t_609[k] = -gl_609[k]
                   + f_0 * il_609[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, gl_610, gl_611, gl_612, gl_613, \
                         gl_614, il_610, il_611, il_612, il_613, \
                         il_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -gl_610[k]
                   + f_0 * il_610[k];

        t_611[k] = -gl_611[k]
                   + f_0 * il_611[k];

        t_612[k] = -gl_612[k]
                   + f_0 * il_612[k];

        t_613[k] = -gl_613[k]
                   + f_0 * il_613[k];

        t_614[k] = -gl_614[k]
                   + f_0 * il_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, gl_615, gl_616, gl_617, gl_618, \
                         gl_619, il_615, il_616, il_617, il_618, \
                         il_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -gl_615[k]
                   + f_0 * il_615[k];

        t_616[k] = -gl_616[k]
                   + f_0 * il_616[k];

        t_617[k] = -gl_617[k]
                   + f_0 * il_617[k];

        t_618[k] = -gl_618[k]
                   + f_0 * il_618[k];

        t_619[k] = -gl_619[k]
                   + f_0 * il_619[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, gl_620, gl_621, gl_622, gl_623, \
                         gl_624, il_620, il_621, il_622, il_623, \
                         il_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -gl_620[k]
                   + f_0 * il_620[k];

        t_621[k] = -gl_621[k]
                   + f_0 * il_621[k];

        t_622[k] = -gl_622[k]
                   + f_0 * il_622[k];

        t_623[k] = -gl_623[k]
                   + f_0 * il_623[k];

        t_624[k] = -gl_624[k]
                   + f_0 * il_624[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, gl_625, gl_626, gl_627, gl_628, \
                         gl_629, il_625, il_626, il_627, il_628, \
                         il_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -gl_625[k]
                   + f_0 * il_625[k];

        t_626[k] = -gl_626[k]
                   + f_0 * il_626[k];

        t_627[k] = -gl_627[k]
                   + f_0 * il_627[k];

        t_628[k] = -gl_628[k]
                   + f_0 * il_628[k];

        t_629[k] = -gl_629[k]
                   + f_0 * il_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, gl_630, gl_631, gl_632, gl_633, \
                         gl_634, il_630, il_631, il_632, il_633, \
                         il_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -gl_630[k]
                   + f_0 * il_630[k];

        t_631[k] = -gl_631[k]
                   + f_0 * il_631[k];

        t_632[k] = -gl_632[k]
                   + f_0 * il_632[k];

        t_633[k] = -gl_633[k]
                   + f_0 * il_633[k];

        t_634[k] = -gl_634[k]
                   + f_0 * il_634[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, gl_635, gl_636, gl_637, gl_638, \
                         gl_639, il_635, il_636, il_637, il_638, \
                         il_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -gl_635[k]
                   + f_0 * il_635[k];

        t_636[k] = -gl_636[k]
                   + f_0 * il_636[k];

        t_637[k] = -gl_637[k]
                   + f_0 * il_637[k];

        t_638[k] = -gl_638[k]
                   + f_0 * il_638[k];

        t_639[k] = -gl_639[k]
                   + f_0 * il_639[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, gl_640, gl_641, gl_642, gl_643, \
                         gl_644, il_640, il_641, il_642, il_643, \
                         il_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -gl_640[k]
                   + f_0 * il_640[k];

        t_641[k] = -gl_641[k]
                   + f_0 * il_641[k];

        t_642[k] = -gl_642[k]
                   + f_0 * il_642[k];

        t_643[k] = -gl_643[k]
                   + f_0 * il_643[k];

        t_644[k] = -gl_644[k]
                   + f_0 * il_644[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, gl_645, gl_646, gl_647, gl_648, \
                         gl_649, il_645, il_646, il_647, il_648, \
                         il_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -gl_645[k]
                   + f_0 * il_645[k];

        t_646[k] = -gl_646[k]
                   + f_0 * il_646[k];

        t_647[k] = -gl_647[k]
                   + f_0 * il_647[k];

        t_648[k] = -gl_648[k]
                   + f_0 * il_648[k];

        t_649[k] = -gl_649[k]
                   + f_0 * il_649[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, gl_650, gl_651, gl_652, gl_653, \
                         gl_654, il_650, il_651, il_652, il_653, \
                         il_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -gl_650[k]
                   + f_0 * il_650[k];

        t_651[k] = -gl_651[k]
                   + f_0 * il_651[k];

        t_652[k] = -gl_652[k]
                   + f_0 * il_652[k];

        t_653[k] = -gl_653[k]
                   + f_0 * il_653[k];

        t_654[k] = -gl_654[k]
                   + f_0 * il_654[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, gl_655, gl_656, gl_657, gl_658, \
                         gl_659, il_655, il_656, il_657, il_658, \
                         il_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -gl_655[k]
                   + f_0 * il_655[k];

        t_656[k] = -gl_656[k]
                   + f_0 * il_656[k];

        t_657[k] = -gl_657[k]
                   + f_0 * il_657[k];

        t_658[k] = -gl_658[k]
                   + f_0 * il_658[k];

        t_659[k] = -gl_659[k]
                   + f_0 * il_659[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, gl_660, gl_661, gl_662, gl_663, \
                         gl_664, il_660, il_661, il_662, il_663, \
                         il_664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -gl_660[k]
                   + f_0 * il_660[k];

        t_661[k] = -gl_661[k]
                   + f_0 * il_661[k];

        t_662[k] = -gl_662[k]
                   + f_0 * il_662[k];

        t_663[k] = -gl_663[k]
                   + f_0 * il_663[k];

        t_664[k] = -gl_664[k]
                   + f_0 * il_664[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, gl_665, gl_666, gl_667, gl_668, \
                         gl_669, il_665, il_666, il_667, il_668, \
                         il_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -gl_665[k]
                   + f_0 * il_665[k];

        t_666[k] = -gl_666[k]
                   + f_0 * il_666[k];

        t_667[k] = -gl_667[k]
                   + f_0 * il_667[k];

        t_668[k] = -gl_668[k]
                   + f_0 * il_668[k];

        t_669[k] = -gl_669[k]
                   + f_0 * il_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, gl_670, gl_671, gl_672, gl_673, \
                         gl_674, il_670, il_671, il_672, il_673, \
                         il_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -gl_670[k]
                   + f_0 * il_670[k];

        t_671[k] = -gl_671[k]
                   + f_0 * il_671[k];

        t_672[k] = -gl_672[k]
                   + f_0 * il_672[k];

        t_673[k] = -gl_673[k]
                   + f_0 * il_673[k];

        t_674[k] = -gl_674[k]
                   + f_0 * il_674[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, t_680, t_681, t_682, il_675, \
                         il_676, il_677, il_678, il_679, il_680, il_681, \
                         il_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = f_0 * il_675[k];

        t_676[k] = f_0 * il_676[k];

        t_677[k] = f_0 * il_677[k];

        t_678[k] = f_0 * il_678[k];

        t_679[k] = f_0 * il_679[k];

        t_680[k] = f_0 * il_680[k];

        t_681[k] = f_0 * il_681[k];

        t_682[k] = f_0 * il_682[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, t_686, t_687, t_688, t_689, t_690, il_683, \
                         il_684, il_685, il_686, il_687, il_688, il_689, \
                         il_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_0 * il_683[k];

        t_684[k] = f_0 * il_684[k];

        t_685[k] = f_0 * il_685[k];

        t_686[k] = f_0 * il_686[k];

        t_687[k] = f_0 * il_687[k];

        t_688[k] = f_0 * il_688[k];

        t_689[k] = f_0 * il_689[k];

        t_690[k] = f_0 * il_690[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, t_695, t_696, t_697, t_698, il_691, \
                         il_692, il_693, il_694, il_695, il_696, il_697, \
                         il_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_0 * il_691[k];

        t_692[k] = f_0 * il_692[k];

        t_693[k] = f_0 * il_693[k];

        t_694[k] = f_0 * il_694[k];

        t_695[k] = f_0 * il_695[k];

        t_696[k] = f_0 * il_696[k];

        t_697[k] = f_0 * il_697[k];

        t_698[k] = f_0 * il_698[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, t_703, t_704, t_705, t_706, il_699, \
                         il_700, il_701, il_702, il_703, il_704, il_705, \
                         il_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_0 * il_699[k];

        t_700[k] = f_0 * il_700[k];

        t_701[k] = f_0 * il_701[k];

        t_702[k] = f_0 * il_702[k];

        t_703[k] = f_0 * il_703[k];

        t_704[k] = f_0 * il_704[k];

        t_705[k] = f_0 * il_705[k];

        t_706[k] = f_0 * il_706[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, t_710, t_711, t_712, t_713, t_714, il_707, \
                         il_708, il_709, il_710, il_711, il_712, il_713, \
                         il_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = f_0 * il_707[k];

        t_708[k] = f_0 * il_708[k];

        t_709[k] = f_0 * il_709[k];

        t_710[k] = f_0 * il_710[k];

        t_711[k] = f_0 * il_711[k];

        t_712[k] = f_0 * il_712[k];

        t_713[k] = f_0 * il_713[k];

        t_714[k] = f_0 * il_714[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, t_720, t_721, t_722, il_715, \
                         il_716, il_717, il_718, il_719, il_720, il_721, \
                         il_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_0 * il_715[k];

        t_716[k] = f_0 * il_716[k];

        t_717[k] = f_0 * il_717[k];

        t_718[k] = f_0 * il_718[k];

        t_719[k] = f_0 * il_719[k];

        t_720[k] = f_0 * il_720[k];

        t_721[k] = f_0 * il_721[k];

        t_722[k] = f_0 * il_722[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, t_726, t_727, t_728, t_729, t_730, il_723, \
                         il_724, il_725, il_726, il_727, il_728, il_729, \
                         il_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_0 * il_723[k];

        t_724[k] = f_0 * il_724[k];

        t_725[k] = f_0 * il_725[k];

        t_726[k] = f_0 * il_726[k];

        t_727[k] = f_0 * il_727[k];

        t_728[k] = f_0 * il_728[k];

        t_729[k] = f_0 * il_729[k];

        t_730[k] = f_0 * il_730[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, t_734, t_735, t_736, t_737, t_738, il_731, \
                         il_732, il_733, il_734, il_735, il_736, il_737, \
                         il_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_0 * il_731[k];

        t_732[k] = f_0 * il_732[k];

        t_733[k] = f_0 * il_733[k];

        t_734[k] = f_0 * il_734[k];

        t_735[k] = f_0 * il_735[k];

        t_736[k] = f_0 * il_736[k];

        t_737[k] = f_0 * il_737[k];

        t_738[k] = f_0 * il_738[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, t_743, t_744, t_745, t_746, il_739, \
                         il_740, il_741, il_742, il_743, il_744, il_745, \
                         il_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_0 * il_739[k];

        t_740[k] = f_0 * il_740[k];

        t_741[k] = f_0 * il_741[k];

        t_742[k] = f_0 * il_742[k];

        t_743[k] = f_0 * il_743[k];

        t_744[k] = f_0 * il_744[k];

        t_745[k] = f_0 * il_745[k];

        t_746[k] = f_0 * il_746[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, t_750, t_751, t_752, t_753, t_754, il_747, \
                         il_748, il_749, il_750, il_751, il_752, il_753, \
                         il_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_0 * il_747[k];

        t_748[k] = f_0 * il_748[k];

        t_749[k] = f_0 * il_749[k];

        t_750[k] = f_0 * il_750[k];

        t_751[k] = f_0 * il_751[k];

        t_752[k] = f_0 * il_752[k];

        t_753[k] = f_0 * il_753[k];

        t_754[k] = f_0 * il_754[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, t_760, t_761, t_762, il_755, \
                         il_756, il_757, il_758, il_759, il_760, il_761, \
                         il_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_0 * il_755[k];

        t_756[k] = f_0 * il_756[k];

        t_757[k] = f_0 * il_757[k];

        t_758[k] = f_0 * il_758[k];

        t_759[k] = f_0 * il_759[k];

        t_760[k] = f_0 * il_760[k];

        t_761[k] = f_0 * il_761[k];

        t_762[k] = f_0 * il_762[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, t_766, t_767, t_768, t_769, t_770, il_763, \
                         il_764, il_765, il_766, il_767, il_768, il_769, \
                         il_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_0 * il_763[k];

        t_764[k] = f_0 * il_764[k];

        t_765[k] = f_0 * il_765[k];

        t_766[k] = f_0 * il_766[k];

        t_767[k] = f_0 * il_767[k];

        t_768[k] = f_0 * il_768[k];

        t_769[k] = f_0 * il_769[k];

        t_770[k] = f_0 * il_770[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, t_774, t_775, t_776, t_777, t_778, il_771, \
                         il_772, il_773, il_774, il_775, il_776, il_777, \
                         il_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_0 * il_771[k];

        t_772[k] = f_0 * il_772[k];

        t_773[k] = f_0 * il_773[k];

        t_774[k] = f_0 * il_774[k];

        t_775[k] = f_0 * il_775[k];

        t_776[k] = f_0 * il_776[k];

        t_777[k] = f_0 * il_777[k];

        t_778[k] = f_0 * il_778[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, t_782, t_783, t_784, t_785, t_786, il_779, \
                         il_780, il_781, il_782, il_783, il_784, il_785, \
                         il_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_0 * il_779[k];

        t_780[k] = f_0 * il_780[k];

        t_781[k] = f_0 * il_781[k];

        t_782[k] = f_0 * il_782[k];

        t_783[k] = f_0 * il_783[k];

        t_784[k] = f_0 * il_784[k];

        t_785[k] = f_0 * il_785[k];

        t_786[k] = f_0 * il_786[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t il, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *il_787 = buffer.data(il + 787);
    const auto *il_788 = buffer.data(il + 788);
    const auto *il_789 = buffer.data(il + 789);
    const auto *il_790 = buffer.data(il + 790);
    const auto *il_791 = buffer.data(il + 791);
    const auto *il_792 = buffer.data(il + 792);
    const auto *il_793 = buffer.data(il + 793);
    const auto *il_794 = buffer.data(il + 794);
    const auto *il_795 = buffer.data(il + 795);
    const auto *il_796 = buffer.data(il + 796);
    const auto *il_797 = buffer.data(il + 797);
    const auto *il_798 = buffer.data(il + 798);
    const auto *il_799 = buffer.data(il + 799);
    const auto *il_800 = buffer.data(il + 800);
    const auto *il_801 = buffer.data(il + 801);
    const auto *il_802 = buffer.data(il + 802);
    const auto *il_803 = buffer.data(il + 803);
    const auto *il_804 = buffer.data(il + 804);
    const auto *il_805 = buffer.data(il + 805);
    const auto *il_806 = buffer.data(il + 806);
    const auto *il_807 = buffer.data(il + 807);
    const auto *il_808 = buffer.data(il + 808);
    const auto *il_809 = buffer.data(il + 809);
    const auto *il_810 = buffer.data(il + 810);
    const auto *il_811 = buffer.data(il + 811);
    const auto *il_812 = buffer.data(il + 812);
    const auto *il_813 = buffer.data(il + 813);
    const auto *il_814 = buffer.data(il + 814);
    const auto *il_815 = buffer.data(il + 815);
    const auto *il_816 = buffer.data(il + 816);
    const auto *il_817 = buffer.data(il + 817);
    const auto *il_818 = buffer.data(il + 818);
    const auto *il_819 = buffer.data(il + 819);
    const auto *il_820 = buffer.data(il + 820);
    const auto *il_821 = buffer.data(il + 821);
    const auto *il_822 = buffer.data(il + 822);
    const auto *il_823 = buffer.data(il + 823);
    const auto *il_824 = buffer.data(il + 824);
    const auto *il_825 = buffer.data(il + 825);
    const auto *il_826 = buffer.data(il + 826);
    const auto *il_827 = buffer.data(il + 827);
    const auto *il_828 = buffer.data(il + 828);
    const auto *il_829 = buffer.data(il + 829);
    const auto *il_830 = buffer.data(il + 830);
    const auto *il_831 = buffer.data(il + 831);
    const auto *il_832 = buffer.data(il + 832);
    const auto *il_833 = buffer.data(il + 833);
    const auto *il_834 = buffer.data(il + 834);
    const auto *il_835 = buffer.data(il + 835);
    const auto *il_836 = buffer.data(il + 836);
    const auto *il_837 = buffer.data(il + 837);
    const auto *il_838 = buffer.data(il + 838);
    const auto *il_839 = buffer.data(il + 839);
    const auto *il_840 = buffer.data(il + 840);
    const auto *il_841 = buffer.data(il + 841);
    const auto *il_842 = buffer.data(il + 842);
    const auto *il_843 = buffer.data(il + 843);
    const auto *il_844 = buffer.data(il + 844);
    const auto *il_845 = buffer.data(il + 845);
    const auto *il_846 = buffer.data(il + 846);
    const auto *il_847 = buffer.data(il + 847);
    const auto *il_848 = buffer.data(il + 848);
    const auto *il_849 = buffer.data(il + 849);
    const auto *il_850 = buffer.data(il + 850);
    const auto *il_851 = buffer.data(il + 851);
    const auto *il_852 = buffer.data(il + 852);
    const auto *il_853 = buffer.data(il + 853);
    const auto *il_854 = buffer.data(il + 854);
    const auto *il_855 = buffer.data(il + 855);
    const auto *il_856 = buffer.data(il + 856);
    const auto *il_857 = buffer.data(il + 857);
    const auto *il_858 = buffer.data(il + 858);
    const auto *il_859 = buffer.data(il + 859);
    const auto *il_860 = buffer.data(il + 860);
    const auto *il_861 = buffer.data(il + 861);
    const auto *il_862 = buffer.data(il + 862);
    const auto *il_863 = buffer.data(il + 863);
    const auto *il_864 = buffer.data(il + 864);
    const auto *il_865 = buffer.data(il + 865);
    const auto *il_866 = buffer.data(il + 866);
    const auto *il_867 = buffer.data(il + 867);
    const auto *il_868 = buffer.data(il + 868);
    const auto *il_869 = buffer.data(il + 869);
    const auto *il_870 = buffer.data(il + 870);
    const auto *il_871 = buffer.data(il + 871);
    const auto *il_872 = buffer.data(il + 872);
    const auto *il_873 = buffer.data(il + 873);
    const auto *il_874 = buffer.data(il + 874);
    const auto *il_875 = buffer.data(il + 875);
    const auto *il_876 = buffer.data(il + 876);
    const auto *il_877 = buffer.data(il + 877);
    const auto *il_878 = buffer.data(il + 878);
    const auto *il_879 = buffer.data(il + 879);
    const auto *il_880 = buffer.data(il + 880);
    const auto *il_881 = buffer.data(il + 881);
    const auto *il_882 = buffer.data(il + 882);
    const auto *il_883 = buffer.data(il + 883);
    const auto *il_884 = buffer.data(il + 884);
    const auto *il_885 = buffer.data(il + 885);
    const auto *il_886 = buffer.data(il + 886);
    const auto *il_887 = buffer.data(il + 887);
    const auto *il_888 = buffer.data(il + 888);
    const auto *il_889 = buffer.data(il + 889);
    const auto *il_890 = buffer.data(il + 890);
    const auto *il_891 = buffer.data(il + 891);
    const auto *il_892 = buffer.data(il + 892);
    const auto *il_893 = buffer.data(il + 893);
    const auto *il_894 = buffer.data(il + 894);
    const auto *il_895 = buffer.data(il + 895);
    const auto *il_896 = buffer.data(il + 896);
    const auto *il_897 = buffer.data(il + 897);
    const auto *il_898 = buffer.data(il + 898);
    const auto *il_899 = buffer.data(il + 899);
    const auto *il_900 = buffer.data(il + 900);
    const auto *il_901 = buffer.data(il + 901);
    const auto *il_902 = buffer.data(il + 902);
    const auto *il_903 = buffer.data(il + 903);
    const auto *il_904 = buffer.data(il + 904);
    const auto *il_905 = buffer.data(il + 905);
    const auto *il_906 = buffer.data(il + 906);
    const auto *il_907 = buffer.data(il + 907);
    const auto *il_908 = buffer.data(il + 908);
    const auto *il_909 = buffer.data(il + 909);
    const auto *il_910 = buffer.data(il + 910);
    const auto *il_911 = buffer.data(il + 911);
    const auto *il_912 = buffer.data(il + 912);
    const auto *il_913 = buffer.data(il + 913);
    const auto *il_914 = buffer.data(il + 914);
    const auto *il_915 = buffer.data(il + 915);
    const auto *il_916 = buffer.data(il + 916);
    const auto *il_917 = buffer.data(il + 917);
    const auto *il_918 = buffer.data(il + 918);
    const auto *il_919 = buffer.data(il + 919);
    const auto *il_920 = buffer.data(il + 920);
    const auto *il_921 = buffer.data(il + 921);
    const auto *il_922 = buffer.data(il + 922);
    const auto *il_923 = buffer.data(il + 923);
    const auto *il_924 = buffer.data(il + 924);
    const auto *il_925 = buffer.data(il + 925);
    const auto *il_926 = buffer.data(il + 926);
    const auto *il_927 = buffer.data(il + 927);
    const auto *il_928 = buffer.data(il + 928);
    const auto *il_929 = buffer.data(il + 929);
    const auto *il_930 = buffer.data(il + 930);
    const auto *il_931 = buffer.data(il + 931);
    const auto *il_932 = buffer.data(il + 932);
    const auto *il_933 = buffer.data(il + 933);
    const auto *il_934 = buffer.data(il + 934);
    const auto *il_935 = buffer.data(il + 935);
    const auto *il_936 = buffer.data(il + 936);
    const auto *il_937 = buffer.data(il + 937);
    const auto *il_938 = buffer.data(il + 938);
    const auto *il_939 = buffer.data(il + 939);
    const auto *il_940 = buffer.data(il + 940);
    const auto *il_941 = buffer.data(il + 941);
    const auto *il_942 = buffer.data(il + 942);
    const auto *il_943 = buffer.data(il + 943);
    const auto *il_944 = buffer.data(il + 944);

#pragma omp simd aligned(t_787, t_788, t_789, t_790, t_791, t_792, t_793, t_794, il_787, \
                         il_788, il_789, il_790, il_791, il_792, il_793, \
                         il_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_0 * il_787[k];

        t_788[k] = f_0 * il_788[k];

        t_789[k] = f_0 * il_789[k];

        t_790[k] = f_0 * il_790[k];

        t_791[k] = f_0 * il_791[k];

        t_792[k] = f_0 * il_792[k];

        t_793[k] = f_0 * il_793[k];

        t_794[k] = f_0 * il_794[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, t_800, t_801, t_802, il_795, \
                         il_796, il_797, il_798, il_799, il_800, il_801, \
                         il_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_0 * il_795[k];

        t_796[k] = f_0 * il_796[k];

        t_797[k] = f_0 * il_797[k];

        t_798[k] = f_0 * il_798[k];

        t_799[k] = f_0 * il_799[k];

        t_800[k] = f_0 * il_800[k];

        t_801[k] = f_0 * il_801[k];

        t_802[k] = f_0 * il_802[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, t_806, t_807, t_808, t_809, t_810, il_803, \
                         il_804, il_805, il_806, il_807, il_808, il_809, \
                         il_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_0 * il_803[k];

        t_804[k] = f_0 * il_804[k];

        t_805[k] = f_0 * il_805[k];

        t_806[k] = f_0 * il_806[k];

        t_807[k] = f_0 * il_807[k];

        t_808[k] = f_0 * il_808[k];

        t_809[k] = f_0 * il_809[k];

        t_810[k] = f_0 * il_810[k];
    }

#pragma omp simd aligned(t_811, t_812, t_813, t_814, t_815, t_816, t_817, t_818, il_811, \
                         il_812, il_813, il_814, il_815, il_816, il_817, \
                         il_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = f_0 * il_811[k];

        t_812[k] = f_0 * il_812[k];

        t_813[k] = f_0 * il_813[k];

        t_814[k] = f_0 * il_814[k];

        t_815[k] = f_0 * il_815[k];

        t_816[k] = f_0 * il_816[k];

        t_817[k] = f_0 * il_817[k];

        t_818[k] = f_0 * il_818[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, t_822, t_823, t_824, t_825, t_826, il_819, \
                         il_820, il_821, il_822, il_823, il_824, il_825, \
                         il_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_0 * il_819[k];

        t_820[k] = f_0 * il_820[k];

        t_821[k] = f_0 * il_821[k];

        t_822[k] = f_0 * il_822[k];

        t_823[k] = f_0 * il_823[k];

        t_824[k] = f_0 * il_824[k];

        t_825[k] = f_0 * il_825[k];

        t_826[k] = f_0 * il_826[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, t_830, t_831, t_832, t_833, t_834, il_827, \
                         il_828, il_829, il_830, il_831, il_832, il_833, \
                         il_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_0 * il_827[k];

        t_828[k] = f_0 * il_828[k];

        t_829[k] = f_0 * il_829[k];

        t_830[k] = f_0 * il_830[k];

        t_831[k] = f_0 * il_831[k];

        t_832[k] = f_0 * il_832[k];

        t_833[k] = f_0 * il_833[k];

        t_834[k] = f_0 * il_834[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, t_838, t_839, t_840, t_841, t_842, il_835, \
                         il_836, il_837, il_838, il_839, il_840, il_841, \
                         il_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = f_0 * il_835[k];

        t_836[k] = f_0 * il_836[k];

        t_837[k] = f_0 * il_837[k];

        t_838[k] = f_0 * il_838[k];

        t_839[k] = f_0 * il_839[k];

        t_840[k] = f_0 * il_840[k];

        t_841[k] = f_0 * il_841[k];

        t_842[k] = f_0 * il_842[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, t_846, t_847, t_848, t_849, t_850, il_843, \
                         il_844, il_845, il_846, il_847, il_848, il_849, \
                         il_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_0 * il_843[k];

        t_844[k] = f_0 * il_844[k];

        t_845[k] = f_0 * il_845[k];

        t_846[k] = f_0 * il_846[k];

        t_847[k] = f_0 * il_847[k];

        t_848[k] = f_0 * il_848[k];

        t_849[k] = f_0 * il_849[k];

        t_850[k] = f_0 * il_850[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, t_854, t_855, t_856, t_857, t_858, il_851, \
                         il_852, il_853, il_854, il_855, il_856, il_857, \
                         il_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_0 * il_851[k];

        t_852[k] = f_0 * il_852[k];

        t_853[k] = f_0 * il_853[k];

        t_854[k] = f_0 * il_854[k];

        t_855[k] = f_0 * il_855[k];

        t_856[k] = f_0 * il_856[k];

        t_857[k] = f_0 * il_857[k];

        t_858[k] = f_0 * il_858[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, t_862, t_863, t_864, t_865, t_866, il_859, \
                         il_860, il_861, il_862, il_863, il_864, il_865, \
                         il_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_0 * il_859[k];

        t_860[k] = f_0 * il_860[k];

        t_861[k] = f_0 * il_861[k];

        t_862[k] = f_0 * il_862[k];

        t_863[k] = f_0 * il_863[k];

        t_864[k] = f_0 * il_864[k];

        t_865[k] = f_0 * il_865[k];

        t_866[k] = f_0 * il_866[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, t_870, t_871, t_872, t_873, t_874, il_867, \
                         il_868, il_869, il_870, il_871, il_872, il_873, \
                         il_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = f_0 * il_867[k];

        t_868[k] = f_0 * il_868[k];

        t_869[k] = f_0 * il_869[k];

        t_870[k] = f_0 * il_870[k];

        t_871[k] = f_0 * il_871[k];

        t_872[k] = f_0 * il_872[k];

        t_873[k] = f_0 * il_873[k];

        t_874[k] = f_0 * il_874[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, t_879, t_880, t_881, t_882, il_875, \
                         il_876, il_877, il_878, il_879, il_880, il_881, \
                         il_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = f_0 * il_875[k];

        t_876[k] = f_0 * il_876[k];

        t_877[k] = f_0 * il_877[k];

        t_878[k] = f_0 * il_878[k];

        t_879[k] = f_0 * il_879[k];

        t_880[k] = f_0 * il_880[k];

        t_881[k] = f_0 * il_881[k];

        t_882[k] = f_0 * il_882[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, t_887, t_888, t_889, t_890, il_883, \
                         il_884, il_885, il_886, il_887, il_888, il_889, \
                         il_890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_0 * il_883[k];

        t_884[k] = f_0 * il_884[k];

        t_885[k] = f_0 * il_885[k];

        t_886[k] = f_0 * il_886[k];

        t_887[k] = f_0 * il_887[k];

        t_888[k] = f_0 * il_888[k];

        t_889[k] = f_0 * il_889[k];

        t_890[k] = f_0 * il_890[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, t_894, t_895, t_896, t_897, t_898, il_891, \
                         il_892, il_893, il_894, il_895, il_896, il_897, \
                         il_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = f_0 * il_891[k];

        t_892[k] = f_0 * il_892[k];

        t_893[k] = f_0 * il_893[k];

        t_894[k] = f_0 * il_894[k];

        t_895[k] = f_0 * il_895[k];

        t_896[k] = f_0 * il_896[k];

        t_897[k] = f_0 * il_897[k];

        t_898[k] = f_0 * il_898[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, t_903, t_904, t_905, t_906, il_899, \
                         il_900, il_901, il_902, il_903, il_904, il_905, \
                         il_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_0 * il_899[k];

        t_900[k] = f_0 * il_900[k];

        t_901[k] = f_0 * il_901[k];

        t_902[k] = f_0 * il_902[k];

        t_903[k] = f_0 * il_903[k];

        t_904[k] = f_0 * il_904[k];

        t_905[k] = f_0 * il_905[k];

        t_906[k] = f_0 * il_906[k];
    }

#pragma omp simd aligned(t_907, t_908, t_909, t_910, t_911, t_912, t_913, t_914, il_907, \
                         il_908, il_909, il_910, il_911, il_912, il_913, \
                         il_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_907[k] = f_0 * il_907[k];

        t_908[k] = f_0 * il_908[k];

        t_909[k] = f_0 * il_909[k];

        t_910[k] = f_0 * il_910[k];

        t_911[k] = f_0 * il_911[k];

        t_912[k] = f_0 * il_912[k];

        t_913[k] = f_0 * il_913[k];

        t_914[k] = f_0 * il_914[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, t_919, t_920, t_921, t_922, il_915, \
                         il_916, il_917, il_918, il_919, il_920, il_921, \
                         il_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_0 * il_915[k];

        t_916[k] = f_0 * il_916[k];

        t_917[k] = f_0 * il_917[k];

        t_918[k] = f_0 * il_918[k];

        t_919[k] = f_0 * il_919[k];

        t_920[k] = f_0 * il_920[k];

        t_921[k] = f_0 * il_921[k];

        t_922[k] = f_0 * il_922[k];
    }

#pragma omp simd aligned(t_923, t_924, t_925, t_926, t_927, t_928, t_929, t_930, il_923, \
                         il_924, il_925, il_926, il_927, il_928, il_929, \
                         il_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_923[k] = f_0 * il_923[k];

        t_924[k] = f_0 * il_924[k];

        t_925[k] = f_0 * il_925[k];

        t_926[k] = f_0 * il_926[k];

        t_927[k] = f_0 * il_927[k];

        t_928[k] = f_0 * il_928[k];

        t_929[k] = f_0 * il_929[k];

        t_930[k] = f_0 * il_930[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, t_935, t_936, t_937, t_938, il_931, \
                         il_932, il_933, il_934, il_935, il_936, il_937, \
                         il_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_0 * il_931[k];

        t_932[k] = f_0 * il_932[k];

        t_933[k] = f_0 * il_933[k];

        t_934[k] = f_0 * il_934[k];

        t_935[k] = f_0 * il_935[k];

        t_936[k] = f_0 * il_936[k];

        t_937[k] = f_0 * il_937[k];

        t_938[k] = f_0 * il_938[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, t_942, t_943, t_944, il_939, il_940, il_941, \
                         il_942, il_943, il_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_0 * il_939[k];

        t_940[k] = f_0 * il_940[k];

        t_941[k] = f_0 * il_941[k];

        t_942[k] = f_0 * il_942[k];

        t_943[k] = f_0 * il_943[k];

        t_944[k] = f_0 * il_944[k];
    }
}

auto
compute_prim_geom_10_hl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t gl, const size_t il,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hl_electron_repulsion_0_piece0(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_0_piece1(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_0_piece2(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_0_piece3(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_0_piece4(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_0_piece5(buffer, target, il, ncols, alpha);
}

static auto
compute_prim_geom_10_hl_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);
    const auto *gl_30 = buffer.data(gl + 30);
    const auto *gl_31 = buffer.data(gl + 31);
    const auto *gl_32 = buffer.data(gl + 32);
    const auto *gl_33 = buffer.data(gl + 33);
    const auto *gl_34 = buffer.data(gl + 34);
    const auto *gl_35 = buffer.data(gl + 35);
    const auto *gl_36 = buffer.data(gl + 36);
    const auto *gl_37 = buffer.data(gl + 37);
    const auto *gl_38 = buffer.data(gl + 38);
    const auto *gl_39 = buffer.data(gl + 39);
    const auto *gl_40 = buffer.data(gl + 40);
    const auto *gl_41 = buffer.data(gl + 41);
    const auto *gl_42 = buffer.data(gl + 42);
    const auto *gl_43 = buffer.data(gl + 43);
    const auto *gl_44 = buffer.data(gl + 44);
    const auto *gl_45 = buffer.data(gl + 45);
    const auto *gl_46 = buffer.data(gl + 46);
    const auto *gl_47 = buffer.data(gl + 47);
    const auto *gl_48 = buffer.data(gl + 48);
    const auto *gl_49 = buffer.data(gl + 49);
    const auto *gl_50 = buffer.data(gl + 50);
    const auto *gl_51 = buffer.data(gl + 51);
    const auto *gl_52 = buffer.data(gl + 52);
    const auto *gl_53 = buffer.data(gl + 53);
    const auto *gl_54 = buffer.data(gl + 54);
    const auto *gl_55 = buffer.data(gl + 55);
    const auto *gl_56 = buffer.data(gl + 56);
    const auto *gl_57 = buffer.data(gl + 57);
    const auto *gl_58 = buffer.data(gl + 58);
    const auto *gl_59 = buffer.data(gl + 59);
    const auto *gl_60 = buffer.data(gl + 60);
    const auto *gl_61 = buffer.data(gl + 61);
    const auto *gl_62 = buffer.data(gl + 62);
    const auto *gl_63 = buffer.data(gl + 63);
    const auto *gl_64 = buffer.data(gl + 64);
    const auto *gl_65 = buffer.data(gl + 65);
    const auto *gl_66 = buffer.data(gl + 66);
    const auto *gl_67 = buffer.data(gl + 67);
    const auto *gl_68 = buffer.data(gl + 68);
    const auto *gl_69 = buffer.data(gl + 69);
    const auto *gl_70 = buffer.data(gl + 70);
    const auto *gl_71 = buffer.data(gl + 71);
    const auto *gl_72 = buffer.data(gl + 72);
    const auto *gl_73 = buffer.data(gl + 73);
    const auto *gl_74 = buffer.data(gl + 74);
    const auto *gl_75 = buffer.data(gl + 75);
    const auto *gl_76 = buffer.data(gl + 76);
    const auto *gl_77 = buffer.data(gl + 77);
    const auto *gl_78 = buffer.data(gl + 78);
    const auto *gl_79 = buffer.data(gl + 79);
    const auto *gl_80 = buffer.data(gl + 80);
    const auto *gl_81 = buffer.data(gl + 81);
    const auto *gl_82 = buffer.data(gl + 82);
    const auto *gl_83 = buffer.data(gl + 83);
    const auto *gl_84 = buffer.data(gl + 84);
    const auto *gl_85 = buffer.data(gl + 85);
    const auto *gl_86 = buffer.data(gl + 86);
    const auto *gl_87 = buffer.data(gl + 87);
    const auto *gl_88 = buffer.data(gl + 88);

    const auto *il_45 = buffer.data(il + 45);
    const auto *il_46 = buffer.data(il + 46);
    const auto *il_47 = buffer.data(il + 47);
    const auto *il_48 = buffer.data(il + 48);
    const auto *il_49 = buffer.data(il + 49);
    const auto *il_50 = buffer.data(il + 50);
    const auto *il_51 = buffer.data(il + 51);
    const auto *il_52 = buffer.data(il + 52);
    const auto *il_53 = buffer.data(il + 53);
    const auto *il_54 = buffer.data(il + 54);
    const auto *il_55 = buffer.data(il + 55);
    const auto *il_56 = buffer.data(il + 56);
    const auto *il_57 = buffer.data(il + 57);
    const auto *il_58 = buffer.data(il + 58);
    const auto *il_59 = buffer.data(il + 59);
    const auto *il_60 = buffer.data(il + 60);
    const auto *il_61 = buffer.data(il + 61);
    const auto *il_62 = buffer.data(il + 62);
    const auto *il_63 = buffer.data(il + 63);
    const auto *il_64 = buffer.data(il + 64);
    const auto *il_65 = buffer.data(il + 65);
    const auto *il_66 = buffer.data(il + 66);
    const auto *il_67 = buffer.data(il + 67);
    const auto *il_68 = buffer.data(il + 68);
    const auto *il_69 = buffer.data(il + 69);
    const auto *il_70 = buffer.data(il + 70);
    const auto *il_71 = buffer.data(il + 71);
    const auto *il_72 = buffer.data(il + 72);
    const auto *il_73 = buffer.data(il + 73);
    const auto *il_74 = buffer.data(il + 74);
    const auto *il_75 = buffer.data(il + 75);
    const auto *il_76 = buffer.data(il + 76);
    const auto *il_77 = buffer.data(il + 77);
    const auto *il_78 = buffer.data(il + 78);
    const auto *il_79 = buffer.data(il + 79);
    const auto *il_80 = buffer.data(il + 80);
    const auto *il_81 = buffer.data(il + 81);
    const auto *il_82 = buffer.data(il + 82);
    const auto *il_83 = buffer.data(il + 83);
    const auto *il_84 = buffer.data(il + 84);
    const auto *il_85 = buffer.data(il + 85);
    const auto *il_86 = buffer.data(il + 86);
    const auto *il_87 = buffer.data(il + 87);
    const auto *il_88 = buffer.data(il + 88);
    const auto *il_89 = buffer.data(il + 89);
    const auto *il_135 = buffer.data(il + 135);
    const auto *il_136 = buffer.data(il + 136);
    const auto *il_137 = buffer.data(il + 137);
    const auto *il_138 = buffer.data(il + 138);
    const auto *il_139 = buffer.data(il + 139);
    const auto *il_140 = buffer.data(il + 140);
    const auto *il_141 = buffer.data(il + 141);
    const auto *il_142 = buffer.data(il + 142);
    const auto *il_143 = buffer.data(il + 143);
    const auto *il_144 = buffer.data(il + 144);
    const auto *il_145 = buffer.data(il + 145);
    const auto *il_146 = buffer.data(il + 146);
    const auto *il_147 = buffer.data(il + 147);
    const auto *il_148 = buffer.data(il + 148);
    const auto *il_149 = buffer.data(il + 149);
    const auto *il_150 = buffer.data(il + 150);
    const auto *il_151 = buffer.data(il + 151);
    const auto *il_152 = buffer.data(il + 152);
    const auto *il_153 = buffer.data(il + 153);
    const auto *il_154 = buffer.data(il + 154);
    const auto *il_155 = buffer.data(il + 155);
    const auto *il_156 = buffer.data(il + 156);
    const auto *il_157 = buffer.data(il + 157);
    const auto *il_158 = buffer.data(il + 158);
    const auto *il_159 = buffer.data(il + 159);
    const auto *il_160 = buffer.data(il + 160);
    const auto *il_161 = buffer.data(il + 161);
    const auto *il_162 = buffer.data(il + 162);
    const auto *il_163 = buffer.data(il + 163);
    const auto *il_164 = buffer.data(il + 164);
    const auto *il_165 = buffer.data(il + 165);
    const auto *il_166 = buffer.data(il + 166);
    const auto *il_167 = buffer.data(il + 167);
    const auto *il_168 = buffer.data(il + 168);
    const auto *il_169 = buffer.data(il + 169);
    const auto *il_170 = buffer.data(il + 170);
    const auto *il_171 = buffer.data(il + 171);
    const auto *il_172 = buffer.data(il + 172);
    const auto *il_173 = buffer.data(il + 173);
    const auto *il_174 = buffer.data(il + 174);
    const auto *il_175 = buffer.data(il + 175);
    const auto *il_176 = buffer.data(il + 176);
    const auto *il_177 = buffer.data(il + 177);
    const auto *il_178 = buffer.data(il + 178);
    const auto *il_179 = buffer.data(il + 179);
    const auto *il_180 = buffer.data(il + 180);
    const auto *il_181 = buffer.data(il + 181);
    const auto *il_182 = buffer.data(il + 182);
    const auto *il_183 = buffer.data(il + 183);
    const auto *il_184 = buffer.data(il + 184);
    const auto *il_185 = buffer.data(il + 185);
    const auto *il_186 = buffer.data(il + 186);
    const auto *il_187 = buffer.data(il + 187);
    const auto *il_188 = buffer.data(il + 188);
    const auto *il_189 = buffer.data(il + 189);
    const auto *il_190 = buffer.data(il + 190);
    const auto *il_191 = buffer.data(il + 191);
    const auto *il_192 = buffer.data(il + 192);
    const auto *il_193 = buffer.data(il + 193);
    const auto *il_194 = buffer.data(il + 194);
    const auto *il_195 = buffer.data(il + 195);
    const auto *il_196 = buffer.data(il + 196);
    const auto *il_197 = buffer.data(il + 197);
    const auto *il_198 = buffer.data(il + 198);
    const auto *il_199 = buffer.data(il + 199);
    const auto *il_200 = buffer.data(il + 200);
    const auto *il_201 = buffer.data(il + 201);
    const auto *il_202 = buffer.data(il + 202);
    const auto *il_203 = buffer.data(il + 203);
    const auto *il_204 = buffer.data(il + 204);
    const auto *il_205 = buffer.data(il + 205);
    const auto *il_206 = buffer.data(il + 206);
    const auto *il_207 = buffer.data(il + 207);
    const auto *il_208 = buffer.data(il + 208);
    const auto *il_209 = buffer.data(il + 209);
    const auto *il_210 = buffer.data(il + 210);
    const auto *il_211 = buffer.data(il + 211);
    const auto *il_212 = buffer.data(il + 212);
    const auto *il_213 = buffer.data(il + 213);
    const auto *il_214 = buffer.data(il + 214);
    const auto *il_215 = buffer.data(il + 215);
    const auto *il_216 = buffer.data(il + 216);
    const auto *il_217 = buffer.data(il + 217);
    const auto *il_218 = buffer.data(il + 218);
    const auto *il_219 = buffer.data(il + 219);
    const auto *il_220 = buffer.data(il + 220);
    const auto *il_221 = buffer.data(il + 221);
    const auto *il_222 = buffer.data(il + 222);
    const auto *il_223 = buffer.data(il + 223);
    const auto *il_224 = buffer.data(il + 224);
    const auto *il_270 = buffer.data(il + 270);
    const auto *il_271 = buffer.data(il + 271);
    const auto *il_272 = buffer.data(il + 272);
    const auto *il_273 = buffer.data(il + 273);
    const auto *il_274 = buffer.data(il + 274);
    const auto *il_275 = buffer.data(il + 275);
    const auto *il_276 = buffer.data(il + 276);
    const auto *il_277 = buffer.data(il + 277);
    const auto *il_278 = buffer.data(il + 278);
    const auto *il_279 = buffer.data(il + 279);
    const auto *il_280 = buffer.data(il + 280);
    const auto *il_281 = buffer.data(il + 281);
    const auto *il_282 = buffer.data(il + 282);
    const auto *il_283 = buffer.data(il + 283);
    const auto *il_284 = buffer.data(il + 284);
    const auto *il_285 = buffer.data(il + 285);
    const auto *il_286 = buffer.data(il + 286);
    const auto *il_287 = buffer.data(il + 287);
    const auto *il_288 = buffer.data(il + 288);
    const auto *il_289 = buffer.data(il + 289);
    const auto *il_290 = buffer.data(il + 290);
    const auto *il_291 = buffer.data(il + 291);
    const auto *il_292 = buffer.data(il + 292);
    const auto *il_293 = buffer.data(il + 293);
    const auto *il_294 = buffer.data(il + 294);
    const auto *il_295 = buffer.data(il + 295);
    const auto *il_296 = buffer.data(il + 296);
    const auto *il_297 = buffer.data(il + 297);
    const auto *il_298 = buffer.data(il + 298);
    const auto *il_299 = buffer.data(il + 299);
    const auto *il_300 = buffer.data(il + 300);
    const auto *il_301 = buffer.data(il + 301);
    const auto *il_302 = buffer.data(il + 302);
    const auto *il_303 = buffer.data(il + 303);
    const auto *il_304 = buffer.data(il + 304);
    const auto *il_305 = buffer.data(il + 305);
    const auto *il_306 = buffer.data(il + 306);
    const auto *il_307 = buffer.data(il + 307);
    const auto *il_308 = buffer.data(il + 308);
    const auto *il_309 = buffer.data(il + 309);
    const auto *il_310 = buffer.data(il + 310);
    const auto *il_311 = buffer.data(il + 311);
    const auto *il_312 = buffer.data(il + 312);
    const auto *il_313 = buffer.data(il + 313);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, il_45, il_46, il_47, il_48, \
                         il_49, il_50, il_51, il_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * il_45[k];

        t_1[k] = f_0 * il_46[k];

        t_2[k] = f_0 * il_47[k];

        t_3[k] = f_0 * il_48[k];

        t_4[k] = f_0 * il_49[k];

        t_5[k] = f_0 * il_50[k];

        t_6[k] = f_0 * il_51[k];

        t_7[k] = f_0 * il_52[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, il_53, il_54, il_55, \
                         il_56, il_57, il_58, il_59, il_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * il_53[k];

        t_9[k] = f_0 * il_54[k];

        t_10[k] = f_0 * il_55[k];

        t_11[k] = f_0 * il_56[k];

        t_12[k] = f_0 * il_57[k];

        t_13[k] = f_0 * il_58[k];

        t_14[k] = f_0 * il_59[k];

        t_15[k] = f_0 * il_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, il_61, il_62, il_63, \
                         il_64, il_65, il_66, il_67, il_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * il_61[k];

        t_17[k] = f_0 * il_62[k];

        t_18[k] = f_0 * il_63[k];

        t_19[k] = f_0 * il_64[k];

        t_20[k] = f_0 * il_65[k];

        t_21[k] = f_0 * il_66[k];

        t_22[k] = f_0 * il_67[k];

        t_23[k] = f_0 * il_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, il_69, il_70, il_71, \
                         il_72, il_73, il_74, il_75, il_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * il_69[k];

        t_25[k] = f_0 * il_70[k];

        t_26[k] = f_0 * il_71[k];

        t_27[k] = f_0 * il_72[k];

        t_28[k] = f_0 * il_73[k];

        t_29[k] = f_0 * il_74[k];

        t_30[k] = f_0 * il_75[k];

        t_31[k] = f_0 * il_76[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, il_77, il_78, il_79, \
                         il_80, il_81, il_82, il_83, il_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * il_77[k];

        t_33[k] = f_0 * il_78[k];

        t_34[k] = f_0 * il_79[k];

        t_35[k] = f_0 * il_80[k];

        t_36[k] = f_0 * il_81[k];

        t_37[k] = f_0 * il_82[k];

        t_38[k] = f_0 * il_83[k];

        t_39[k] = f_0 * il_84[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, gl_0, gl_1, il_85, il_86, \
                         il_87, il_88, il_89, il_135, il_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * il_85[k];

        t_41[k] = f_0 * il_86[k];

        t_42[k] = f_0 * il_87[k];

        t_43[k] = f_0 * il_88[k];

        t_44[k] = f_0 * il_89[k];

        t_45[k] = -gl_0[k]
                  + f_0 * il_135[k];

        t_46[k] = -gl_1[k]
                  + f_0 * il_136[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, gl_2, gl_3, gl_4, gl_5, gl_6, il_137, \
                         il_138, il_139, il_140, il_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -gl_2[k]
                  + f_0 * il_137[k];

        t_48[k] = -gl_3[k]
                  + f_0 * il_138[k];

        t_49[k] = -gl_4[k]
                  + f_0 * il_139[k];

        t_50[k] = -gl_5[k]
                  + f_0 * il_140[k];

        t_51[k] = -gl_6[k]
                  + f_0 * il_141[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, gl_7, gl_8, gl_9, gl_10, gl_11, il_142, \
                         il_143, il_144, il_145, il_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = -gl_7[k]
                  + f_0 * il_142[k];

        t_53[k] = -gl_8[k]
                  + f_0 * il_143[k];

        t_54[k] = -gl_9[k]
                  + f_0 * il_144[k];

        t_55[k] = -gl_10[k]
                  + f_0 * il_145[k];

        t_56[k] = -gl_11[k]
                  + f_0 * il_146[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, gl_12, gl_13, gl_14, gl_15, gl_16, \
                         il_147, il_148, il_149, il_150, il_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -gl_12[k]
                  + f_0 * il_147[k];

        t_58[k] = -gl_13[k]
                  + f_0 * il_148[k];

        t_59[k] = -gl_14[k]
                  + f_0 * il_149[k];

        t_60[k] = -gl_15[k]
                  + f_0 * il_150[k];

        t_61[k] = -gl_16[k]
                  + f_0 * il_151[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, gl_17, gl_18, gl_19, gl_20, gl_21, \
                         il_152, il_153, il_154, il_155, il_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -gl_17[k]
                  + f_0 * il_152[k];

        t_63[k] = -gl_18[k]
                  + f_0 * il_153[k];

        t_64[k] = -gl_19[k]
                  + f_0 * il_154[k];

        t_65[k] = -gl_20[k]
                  + f_0 * il_155[k];

        t_66[k] = -gl_21[k]
                  + f_0 * il_156[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, gl_22, gl_23, gl_24, gl_25, gl_26, \
                         il_157, il_158, il_159, il_160, il_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -gl_22[k]
                  + f_0 * il_157[k];

        t_68[k] = -gl_23[k]
                  + f_0 * il_158[k];

        t_69[k] = -gl_24[k]
                  + f_0 * il_159[k];

        t_70[k] = -gl_25[k]
                  + f_0 * il_160[k];

        t_71[k] = -gl_26[k]
                  + f_0 * il_161[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, gl_27, gl_28, gl_29, gl_30, gl_31, \
                         il_162, il_163, il_164, il_165, il_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -gl_27[k]
                  + f_0 * il_162[k];

        t_73[k] = -gl_28[k]
                  + f_0 * il_163[k];

        t_74[k] = -gl_29[k]
                  + f_0 * il_164[k];

        t_75[k] = -gl_30[k]
                  + f_0 * il_165[k];

        t_76[k] = -gl_31[k]
                  + f_0 * il_166[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, gl_32, gl_33, gl_34, gl_35, gl_36, \
                         il_167, il_168, il_169, il_170, il_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -gl_32[k]
                  + f_0 * il_167[k];

        t_78[k] = -gl_33[k]
                  + f_0 * il_168[k];

        t_79[k] = -gl_34[k]
                  + f_0 * il_169[k];

        t_80[k] = -gl_35[k]
                  + f_0 * il_170[k];

        t_81[k] = -gl_36[k]
                  + f_0 * il_171[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, gl_37, gl_38, gl_39, gl_40, gl_41, \
                         il_172, il_173, il_174, il_175, il_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -gl_37[k]
                  + f_0 * il_172[k];

        t_83[k] = -gl_38[k]
                  + f_0 * il_173[k];

        t_84[k] = -gl_39[k]
                  + f_0 * il_174[k];

        t_85[k] = -gl_40[k]
                  + f_0 * il_175[k];

        t_86[k] = -gl_41[k]
                  + f_0 * il_176[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, gl_42, gl_43, gl_44, il_177, \
                         il_178, il_179, il_180, il_181, il_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -gl_42[k]
                  + f_0 * il_177[k];

        t_88[k] = -gl_43[k]
                  + f_0 * il_178[k];

        t_89[k] = -gl_44[k]
                  + f_0 * il_179[k];

        t_90[k] = f_0 * il_180[k];

        t_91[k] = f_0 * il_181[k];

        t_92[k] = f_0 * il_182[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, il_183, il_184, \
                         il_185, il_186, il_187, il_188, il_189, \
                         il_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * il_183[k];

        t_94[k] = f_0 * il_184[k];

        t_95[k] = f_0 * il_185[k];

        t_96[k] = f_0 * il_186[k];

        t_97[k] = f_0 * il_187[k];

        t_98[k] = f_0 * il_188[k];

        t_99[k] = f_0 * il_189[k];

        t_100[k] = f_0 * il_190[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, t_107, t_108, il_191, \
                         il_192, il_193, il_194, il_195, il_196, il_197, \
                         il_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * il_191[k];

        t_102[k] = f_0 * il_192[k];

        t_103[k] = f_0 * il_193[k];

        t_104[k] = f_0 * il_194[k];

        t_105[k] = f_0 * il_195[k];

        t_106[k] = f_0 * il_196[k];

        t_107[k] = f_0 * il_197[k];

        t_108[k] = f_0 * il_198[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, t_114, t_115, t_116, il_199, \
                         il_200, il_201, il_202, il_203, il_204, il_205, \
                         il_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_0 * il_199[k];

        t_110[k] = f_0 * il_200[k];

        t_111[k] = f_0 * il_201[k];

        t_112[k] = f_0 * il_202[k];

        t_113[k] = f_0 * il_203[k];

        t_114[k] = f_0 * il_204[k];

        t_115[k] = f_0 * il_205[k];

        t_116[k] = f_0 * il_206[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, t_122, t_123, t_124, il_207, \
                         il_208, il_209, il_210, il_211, il_212, il_213, \
                         il_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_0 * il_207[k];

        t_118[k] = f_0 * il_208[k];

        t_119[k] = f_0 * il_209[k];

        t_120[k] = f_0 * il_210[k];

        t_121[k] = f_0 * il_211[k];

        t_122[k] = f_0 * il_212[k];

        t_123[k] = f_0 * il_213[k];

        t_124[k] = f_0 * il_214[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, t_131, t_132, il_215, \
                         il_216, il_217, il_218, il_219, il_220, il_221, \
                         il_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * il_215[k];

        t_126[k] = f_0 * il_216[k];

        t_127[k] = f_0 * il_217[k];

        t_128[k] = f_0 * il_218[k];

        t_129[k] = f_0 * il_219[k];

        t_130[k] = f_0 * il_220[k];

        t_131[k] = f_0 * il_221[k];

        t_132[k] = f_0 * il_222[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, t_138, gl_45, gl_46, gl_47, gl_48, \
                         il_223, il_224, il_270, il_271, il_272, \
                         il_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_0 * il_223[k];

        t_134[k] = f_0 * il_224[k];

        t_135[k] = -2.0 * gl_45[k]
                   + f_0 * il_270[k];

        t_136[k] = -2.0 * gl_46[k]
                   + f_0 * il_271[k];

        t_137[k] = -2.0 * gl_47[k]
                   + f_0 * il_272[k];

        t_138[k] = -2.0 * gl_48[k]
                   + f_0 * il_273[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, gl_49, gl_50, gl_51, gl_52, gl_53, \
                         il_274, il_275, il_276, il_277, il_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -2.0 * gl_49[k]
                   + f_0 * il_274[k];

        t_140[k] = -2.0 * gl_50[k]
                   + f_0 * il_275[k];

        t_141[k] = -2.0 * gl_51[k]
                   + f_0 * il_276[k];

        t_142[k] = -2.0 * gl_52[k]
                   + f_0 * il_277[k];

        t_143[k] = -2.0 * gl_53[k]
                   + f_0 * il_278[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, gl_54, gl_55, gl_56, gl_57, gl_58, \
                         il_279, il_280, il_281, il_282, il_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -2.0 * gl_54[k]
                   + f_0 * il_279[k];

        t_145[k] = -2.0 * gl_55[k]
                   + f_0 * il_280[k];

        t_146[k] = -2.0 * gl_56[k]
                   + f_0 * il_281[k];

        t_147[k] = -2.0 * gl_57[k]
                   + f_0 * il_282[k];

        t_148[k] = -2.0 * gl_58[k]
                   + f_0 * il_283[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, gl_59, gl_60, gl_61, gl_62, gl_63, \
                         il_284, il_285, il_286, il_287, il_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -2.0 * gl_59[k]
                   + f_0 * il_284[k];

        t_150[k] = -2.0 * gl_60[k]
                   + f_0 * il_285[k];

        t_151[k] = -2.0 * gl_61[k]
                   + f_0 * il_286[k];

        t_152[k] = -2.0 * gl_62[k]
                   + f_0 * il_287[k];

        t_153[k] = -2.0 * gl_63[k]
                   + f_0 * il_288[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, gl_64, gl_65, gl_66, gl_67, gl_68, \
                         il_289, il_290, il_291, il_292, il_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -2.0 * gl_64[k]
                   + f_0 * il_289[k];

        t_155[k] = -2.0 * gl_65[k]
                   + f_0 * il_290[k];

        t_156[k] = -2.0 * gl_66[k]
                   + f_0 * il_291[k];

        t_157[k] = -2.0 * gl_67[k]
                   + f_0 * il_292[k];

        t_158[k] = -2.0 * gl_68[k]
                   + f_0 * il_293[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, gl_69, gl_70, gl_71, gl_72, gl_73, \
                         il_294, il_295, il_296, il_297, il_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -2.0 * gl_69[k]
                   + f_0 * il_294[k];

        t_160[k] = -2.0 * gl_70[k]
                   + f_0 * il_295[k];

        t_161[k] = -2.0 * gl_71[k]
                   + f_0 * il_296[k];

        t_162[k] = -2.0 * gl_72[k]
                   + f_0 * il_297[k];

        t_163[k] = -2.0 * gl_73[k]
                   + f_0 * il_298[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, gl_74, gl_75, gl_76, gl_77, gl_78, \
                         il_299, il_300, il_301, il_302, il_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -2.0 * gl_74[k]
                   + f_0 * il_299[k];

        t_165[k] = -2.0 * gl_75[k]
                   + f_0 * il_300[k];

        t_166[k] = -2.0 * gl_76[k]
                   + f_0 * il_301[k];

        t_167[k] = -2.0 * gl_77[k]
                   + f_0 * il_302[k];

        t_168[k] = -2.0 * gl_78[k]
                   + f_0 * il_303[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, gl_79, gl_80, gl_81, gl_82, gl_83, \
                         il_304, il_305, il_306, il_307, il_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -2.0 * gl_79[k]
                   + f_0 * il_304[k];

        t_170[k] = -2.0 * gl_80[k]
                   + f_0 * il_305[k];

        t_171[k] = -2.0 * gl_81[k]
                   + f_0 * il_306[k];

        t_172[k] = -2.0 * gl_82[k]
                   + f_0 * il_307[k];

        t_173[k] = -2.0 * gl_83[k]
                   + f_0 * il_308[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, gl_84, gl_85, gl_86, gl_87, gl_88, \
                         il_309, il_310, il_311, il_312, il_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -2.0 * gl_84[k]
                   + f_0 * il_309[k];

        t_175[k] = -2.0 * gl_85[k]
                   + f_0 * il_310[k];

        t_176[k] = -2.0 * gl_86[k]
                   + f_0 * il_311[k];

        t_177[k] = -2.0 * gl_87[k]
                   + f_0 * il_312[k];

        t_178[k] = -2.0 * gl_88[k]
                   + f_0 * il_313[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_89 = buffer.data(gl + 89);
    const auto *gl_90 = buffer.data(gl + 90);
    const auto *gl_91 = buffer.data(gl + 91);
    const auto *gl_92 = buffer.data(gl + 92);
    const auto *gl_93 = buffer.data(gl + 93);
    const auto *gl_94 = buffer.data(gl + 94);
    const auto *gl_95 = buffer.data(gl + 95);
    const auto *gl_96 = buffer.data(gl + 96);
    const auto *gl_97 = buffer.data(gl + 97);
    const auto *gl_98 = buffer.data(gl + 98);
    const auto *gl_99 = buffer.data(gl + 99);
    const auto *gl_100 = buffer.data(gl + 100);
    const auto *gl_101 = buffer.data(gl + 101);
    const auto *gl_102 = buffer.data(gl + 102);
    const auto *gl_103 = buffer.data(gl + 103);
    const auto *gl_104 = buffer.data(gl + 104);
    const auto *gl_105 = buffer.data(gl + 105);
    const auto *gl_106 = buffer.data(gl + 106);
    const auto *gl_107 = buffer.data(gl + 107);
    const auto *gl_108 = buffer.data(gl + 108);
    const auto *gl_109 = buffer.data(gl + 109);
    const auto *gl_110 = buffer.data(gl + 110);
    const auto *gl_111 = buffer.data(gl + 111);
    const auto *gl_112 = buffer.data(gl + 112);
    const auto *gl_113 = buffer.data(gl + 113);
    const auto *gl_114 = buffer.data(gl + 114);
    const auto *gl_115 = buffer.data(gl + 115);
    const auto *gl_116 = buffer.data(gl + 116);
    const auto *gl_117 = buffer.data(gl + 117);
    const auto *gl_118 = buffer.data(gl + 118);
    const auto *gl_119 = buffer.data(gl + 119);
    const auto *gl_120 = buffer.data(gl + 120);
    const auto *gl_121 = buffer.data(gl + 121);
    const auto *gl_122 = buffer.data(gl + 122);
    const auto *gl_123 = buffer.data(gl + 123);
    const auto *gl_124 = buffer.data(gl + 124);
    const auto *gl_125 = buffer.data(gl + 125);
    const auto *gl_126 = buffer.data(gl + 126);
    const auto *gl_127 = buffer.data(gl + 127);
    const auto *gl_128 = buffer.data(gl + 128);
    const auto *gl_129 = buffer.data(gl + 129);
    const auto *gl_130 = buffer.data(gl + 130);
    const auto *gl_131 = buffer.data(gl + 131);
    const auto *gl_132 = buffer.data(gl + 132);
    const auto *gl_133 = buffer.data(gl + 133);
    const auto *gl_134 = buffer.data(gl + 134);
    const auto *gl_135 = buffer.data(gl + 135);
    const auto *gl_136 = buffer.data(gl + 136);
    const auto *gl_137 = buffer.data(gl + 137);
    const auto *gl_138 = buffer.data(gl + 138);
    const auto *gl_139 = buffer.data(gl + 139);
    const auto *gl_140 = buffer.data(gl + 140);
    const auto *gl_141 = buffer.data(gl + 141);
    const auto *gl_142 = buffer.data(gl + 142);
    const auto *gl_143 = buffer.data(gl + 143);
    const auto *gl_144 = buffer.data(gl + 144);
    const auto *gl_145 = buffer.data(gl + 145);
    const auto *gl_146 = buffer.data(gl + 146);
    const auto *gl_147 = buffer.data(gl + 147);
    const auto *gl_148 = buffer.data(gl + 148);
    const auto *gl_149 = buffer.data(gl + 149);
    const auto *gl_150 = buffer.data(gl + 150);
    const auto *gl_151 = buffer.data(gl + 151);
    const auto *gl_152 = buffer.data(gl + 152);
    const auto *gl_153 = buffer.data(gl + 153);
    const auto *gl_154 = buffer.data(gl + 154);
    const auto *gl_155 = buffer.data(gl + 155);
    const auto *gl_156 = buffer.data(gl + 156);
    const auto *gl_157 = buffer.data(gl + 157);
    const auto *gl_158 = buffer.data(gl + 158);
    const auto *gl_159 = buffer.data(gl + 159);
    const auto *gl_160 = buffer.data(gl + 160);
    const auto *gl_161 = buffer.data(gl + 161);
    const auto *gl_162 = buffer.data(gl + 162);
    const auto *gl_163 = buffer.data(gl + 163);
    const auto *gl_164 = buffer.data(gl + 164);
    const auto *gl_165 = buffer.data(gl + 165);
    const auto *gl_166 = buffer.data(gl + 166);
    const auto *gl_167 = buffer.data(gl + 167);
    const auto *gl_168 = buffer.data(gl + 168);
    const auto *gl_169 = buffer.data(gl + 169);
    const auto *gl_170 = buffer.data(gl + 170);
    const auto *gl_171 = buffer.data(gl + 171);
    const auto *gl_172 = buffer.data(gl + 172);
    const auto *gl_173 = buffer.data(gl + 173);
    const auto *gl_174 = buffer.data(gl + 174);
    const auto *gl_175 = buffer.data(gl + 175);
    const auto *gl_176 = buffer.data(gl + 176);
    const auto *gl_177 = buffer.data(gl + 177);
    const auto *gl_178 = buffer.data(gl + 178);
    const auto *gl_179 = buffer.data(gl + 179);
    const auto *gl_180 = buffer.data(gl + 180);
    const auto *gl_181 = buffer.data(gl + 181);
    const auto *gl_182 = buffer.data(gl + 182);
    const auto *gl_183 = buffer.data(gl + 183);
    const auto *gl_184 = buffer.data(gl + 184);
    const auto *gl_185 = buffer.data(gl + 185);
    const auto *gl_186 = buffer.data(gl + 186);
    const auto *gl_187 = buffer.data(gl + 187);
    const auto *gl_188 = buffer.data(gl + 188);
    const auto *gl_189 = buffer.data(gl + 189);
    const auto *gl_190 = buffer.data(gl + 190);
    const auto *gl_191 = buffer.data(gl + 191);
    const auto *gl_192 = buffer.data(gl + 192);
    const auto *gl_193 = buffer.data(gl + 193);
    const auto *gl_194 = buffer.data(gl + 194);
    const auto *gl_195 = buffer.data(gl + 195);
    const auto *gl_196 = buffer.data(gl + 196);
    const auto *gl_197 = buffer.data(gl + 197);
    const auto *gl_198 = buffer.data(gl + 198);
    const auto *gl_199 = buffer.data(gl + 199);
    const auto *gl_200 = buffer.data(gl + 200);
    const auto *gl_201 = buffer.data(gl + 201);
    const auto *gl_202 = buffer.data(gl + 202);
    const auto *gl_203 = buffer.data(gl + 203);
    const auto *gl_204 = buffer.data(gl + 204);
    const auto *gl_205 = buffer.data(gl + 205);
    const auto *gl_206 = buffer.data(gl + 206);
    const auto *gl_207 = buffer.data(gl + 207);
    const auto *gl_208 = buffer.data(gl + 208);
    const auto *gl_209 = buffer.data(gl + 209);

    const auto *il_314 = buffer.data(il + 314);
    const auto *il_315 = buffer.data(il + 315);
    const auto *il_316 = buffer.data(il + 316);
    const auto *il_317 = buffer.data(il + 317);
    const auto *il_318 = buffer.data(il + 318);
    const auto *il_319 = buffer.data(il + 319);
    const auto *il_320 = buffer.data(il + 320);
    const auto *il_321 = buffer.data(il + 321);
    const auto *il_322 = buffer.data(il + 322);
    const auto *il_323 = buffer.data(il + 323);
    const auto *il_324 = buffer.data(il + 324);
    const auto *il_325 = buffer.data(il + 325);
    const auto *il_326 = buffer.data(il + 326);
    const auto *il_327 = buffer.data(il + 327);
    const auto *il_328 = buffer.data(il + 328);
    const auto *il_329 = buffer.data(il + 329);
    const auto *il_330 = buffer.data(il + 330);
    const auto *il_331 = buffer.data(il + 331);
    const auto *il_332 = buffer.data(il + 332);
    const auto *il_333 = buffer.data(il + 333);
    const auto *il_334 = buffer.data(il + 334);
    const auto *il_335 = buffer.data(il + 335);
    const auto *il_336 = buffer.data(il + 336);
    const auto *il_337 = buffer.data(il + 337);
    const auto *il_338 = buffer.data(il + 338);
    const auto *il_339 = buffer.data(il + 339);
    const auto *il_340 = buffer.data(il + 340);
    const auto *il_341 = buffer.data(il + 341);
    const auto *il_342 = buffer.data(il + 342);
    const auto *il_343 = buffer.data(il + 343);
    const auto *il_344 = buffer.data(il + 344);
    const auto *il_345 = buffer.data(il + 345);
    const auto *il_346 = buffer.data(il + 346);
    const auto *il_347 = buffer.data(il + 347);
    const auto *il_348 = buffer.data(il + 348);
    const auto *il_349 = buffer.data(il + 349);
    const auto *il_350 = buffer.data(il + 350);
    const auto *il_351 = buffer.data(il + 351);
    const auto *il_352 = buffer.data(il + 352);
    const auto *il_353 = buffer.data(il + 353);
    const auto *il_354 = buffer.data(il + 354);
    const auto *il_355 = buffer.data(il + 355);
    const auto *il_356 = buffer.data(il + 356);
    const auto *il_357 = buffer.data(il + 357);
    const auto *il_358 = buffer.data(il + 358);
    const auto *il_359 = buffer.data(il + 359);
    const auto *il_360 = buffer.data(il + 360);
    const auto *il_361 = buffer.data(il + 361);
    const auto *il_362 = buffer.data(il + 362);
    const auto *il_363 = buffer.data(il + 363);
    const auto *il_364 = buffer.data(il + 364);
    const auto *il_365 = buffer.data(il + 365);
    const auto *il_366 = buffer.data(il + 366);
    const auto *il_367 = buffer.data(il + 367);
    const auto *il_368 = buffer.data(il + 368);
    const auto *il_369 = buffer.data(il + 369);
    const auto *il_370 = buffer.data(il + 370);
    const auto *il_371 = buffer.data(il + 371);
    const auto *il_372 = buffer.data(il + 372);
    const auto *il_373 = buffer.data(il + 373);
    const auto *il_374 = buffer.data(il + 374);
    const auto *il_375 = buffer.data(il + 375);
    const auto *il_376 = buffer.data(il + 376);
    const auto *il_377 = buffer.data(il + 377);
    const auto *il_378 = buffer.data(il + 378);
    const auto *il_379 = buffer.data(il + 379);
    const auto *il_380 = buffer.data(il + 380);
    const auto *il_381 = buffer.data(il + 381);
    const auto *il_382 = buffer.data(il + 382);
    const auto *il_383 = buffer.data(il + 383);
    const auto *il_384 = buffer.data(il + 384);
    const auto *il_385 = buffer.data(il + 385);
    const auto *il_386 = buffer.data(il + 386);
    const auto *il_387 = buffer.data(il + 387);
    const auto *il_388 = buffer.data(il + 388);
    const auto *il_389 = buffer.data(il + 389);
    const auto *il_390 = buffer.data(il + 390);
    const auto *il_391 = buffer.data(il + 391);
    const auto *il_392 = buffer.data(il + 392);
    const auto *il_393 = buffer.data(il + 393);
    const auto *il_394 = buffer.data(il + 394);
    const auto *il_395 = buffer.data(il + 395);
    const auto *il_396 = buffer.data(il + 396);
    const auto *il_397 = buffer.data(il + 397);
    const auto *il_398 = buffer.data(il + 398);
    const auto *il_399 = buffer.data(il + 399);
    const auto *il_400 = buffer.data(il + 400);
    const auto *il_401 = buffer.data(il + 401);
    const auto *il_402 = buffer.data(il + 402);
    const auto *il_403 = buffer.data(il + 403);
    const auto *il_404 = buffer.data(il + 404);
    const auto *il_450 = buffer.data(il + 450);
    const auto *il_451 = buffer.data(il + 451);
    const auto *il_452 = buffer.data(il + 452);
    const auto *il_453 = buffer.data(il + 453);
    const auto *il_454 = buffer.data(il + 454);
    const auto *il_455 = buffer.data(il + 455);
    const auto *il_456 = buffer.data(il + 456);
    const auto *il_457 = buffer.data(il + 457);
    const auto *il_458 = buffer.data(il + 458);
    const auto *il_459 = buffer.data(il + 459);
    const auto *il_460 = buffer.data(il + 460);
    const auto *il_461 = buffer.data(il + 461);
    const auto *il_462 = buffer.data(il + 462);
    const auto *il_463 = buffer.data(il + 463);
    const auto *il_464 = buffer.data(il + 464);
    const auto *il_465 = buffer.data(il + 465);
    const auto *il_466 = buffer.data(il + 466);
    const auto *il_467 = buffer.data(il + 467);
    const auto *il_468 = buffer.data(il + 468);
    const auto *il_469 = buffer.data(il + 469);
    const auto *il_470 = buffer.data(il + 470);
    const auto *il_471 = buffer.data(il + 471);
    const auto *il_472 = buffer.data(il + 472);
    const auto *il_473 = buffer.data(il + 473);
    const auto *il_474 = buffer.data(il + 474);
    const auto *il_475 = buffer.data(il + 475);
    const auto *il_476 = buffer.data(il + 476);
    const auto *il_477 = buffer.data(il + 477);
    const auto *il_478 = buffer.data(il + 478);
    const auto *il_479 = buffer.data(il + 479);
    const auto *il_480 = buffer.data(il + 480);
    const auto *il_481 = buffer.data(il + 481);
    const auto *il_482 = buffer.data(il + 482);
    const auto *il_483 = buffer.data(il + 483);
    const auto *il_484 = buffer.data(il + 484);
    const auto *il_485 = buffer.data(il + 485);
    const auto *il_486 = buffer.data(il + 486);
    const auto *il_487 = buffer.data(il + 487);
    const auto *il_488 = buffer.data(il + 488);
    const auto *il_489 = buffer.data(il + 489);
    const auto *il_490 = buffer.data(il + 490);
    const auto *il_491 = buffer.data(il + 491);
    const auto *il_492 = buffer.data(il + 492);
    const auto *il_493 = buffer.data(il + 493);
    const auto *il_494 = buffer.data(il + 494);
    const auto *il_495 = buffer.data(il + 495);
    const auto *il_496 = buffer.data(il + 496);
    const auto *il_497 = buffer.data(il + 497);
    const auto *il_498 = buffer.data(il + 498);
    const auto *il_499 = buffer.data(il + 499);
    const auto *il_500 = buffer.data(il + 500);
    const auto *il_501 = buffer.data(il + 501);
    const auto *il_502 = buffer.data(il + 502);
    const auto *il_503 = buffer.data(il + 503);
    const auto *il_504 = buffer.data(il + 504);
    const auto *il_505 = buffer.data(il + 505);
    const auto *il_506 = buffer.data(il + 506);
    const auto *il_507 = buffer.data(il + 507);
    const auto *il_508 = buffer.data(il + 508);
    const auto *il_509 = buffer.data(il + 509);
    const auto *il_510 = buffer.data(il + 510);
    const auto *il_511 = buffer.data(il + 511);
    const auto *il_512 = buffer.data(il + 512);
    const auto *il_513 = buffer.data(il + 513);
    const auto *il_514 = buffer.data(il + 514);
    const auto *il_515 = buffer.data(il + 515);
    const auto *il_516 = buffer.data(il + 516);
    const auto *il_517 = buffer.data(il + 517);
    const auto *il_518 = buffer.data(il + 518);
    const auto *il_519 = buffer.data(il + 519);
    const auto *il_520 = buffer.data(il + 520);
    const auto *il_521 = buffer.data(il + 521);
    const auto *il_522 = buffer.data(il + 522);
    const auto *il_523 = buffer.data(il + 523);
    const auto *il_524 = buffer.data(il + 524);

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, gl_89, gl_90, gl_91, gl_92, gl_93, \
                         il_314, il_315, il_316, il_317, il_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -2.0 * gl_89[k]
                   + f_0 * il_314[k];

        t_180[k] = -gl_90[k]
                   + f_0 * il_315[k];

        t_181[k] = -gl_91[k]
                   + f_0 * il_316[k];

        t_182[k] = -gl_92[k]
                   + f_0 * il_317[k];

        t_183[k] = -gl_93[k]
                   + f_0 * il_318[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, gl_94, gl_95, gl_96, gl_97, gl_98, \
                         il_319, il_320, il_321, il_322, il_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -gl_94[k]
                   + f_0 * il_319[k];

        t_185[k] = -gl_95[k]
                   + f_0 * il_320[k];

        t_186[k] = -gl_96[k]
                   + f_0 * il_321[k];

        t_187[k] = -gl_97[k]
                   + f_0 * il_322[k];

        t_188[k] = -gl_98[k]
                   + f_0 * il_323[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, gl_99, gl_100, gl_101, gl_102, \
                         gl_103, il_324, il_325, il_326, il_327, \
                         il_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = -gl_99[k]
                   + f_0 * il_324[k];

        t_190[k] = -gl_100[k]
                   + f_0 * il_325[k];

        t_191[k] = -gl_101[k]
                   + f_0 * il_326[k];

        t_192[k] = -gl_102[k]
                   + f_0 * il_327[k];

        t_193[k] = -gl_103[k]
                   + f_0 * il_328[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, gl_104, gl_105, gl_106, gl_107, \
                         gl_108, il_329, il_330, il_331, il_332, \
                         il_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = -gl_104[k]
                   + f_0 * il_329[k];

        t_195[k] = -gl_105[k]
                   + f_0 * il_330[k];

        t_196[k] = -gl_106[k]
                   + f_0 * il_331[k];

        t_197[k] = -gl_107[k]
                   + f_0 * il_332[k];

        t_198[k] = -gl_108[k]
                   + f_0 * il_333[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, gl_109, gl_110, gl_111, gl_112, \
                         gl_113, il_334, il_335, il_336, il_337, \
                         il_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = -gl_109[k]
                   + f_0 * il_334[k];

        t_200[k] = -gl_110[k]
                   + f_0 * il_335[k];

        t_201[k] = -gl_111[k]
                   + f_0 * il_336[k];

        t_202[k] = -gl_112[k]
                   + f_0 * il_337[k];

        t_203[k] = -gl_113[k]
                   + f_0 * il_338[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, gl_114, gl_115, gl_116, gl_117, \
                         gl_118, il_339, il_340, il_341, il_342, \
                         il_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = -gl_114[k]
                   + f_0 * il_339[k];

        t_205[k] = -gl_115[k]
                   + f_0 * il_340[k];

        t_206[k] = -gl_116[k]
                   + f_0 * il_341[k];

        t_207[k] = -gl_117[k]
                   + f_0 * il_342[k];

        t_208[k] = -gl_118[k]
                   + f_0 * il_343[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, gl_119, gl_120, gl_121, gl_122, \
                         gl_123, il_344, il_345, il_346, il_347, \
                         il_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = -gl_119[k]
                   + f_0 * il_344[k];

        t_210[k] = -gl_120[k]
                   + f_0 * il_345[k];

        t_211[k] = -gl_121[k]
                   + f_0 * il_346[k];

        t_212[k] = -gl_122[k]
                   + f_0 * il_347[k];

        t_213[k] = -gl_123[k]
                   + f_0 * il_348[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, gl_124, gl_125, gl_126, gl_127, \
                         gl_128, il_349, il_350, il_351, il_352, \
                         il_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = -gl_124[k]
                   + f_0 * il_349[k];

        t_215[k] = -gl_125[k]
                   + f_0 * il_350[k];

        t_216[k] = -gl_126[k]
                   + f_0 * il_351[k];

        t_217[k] = -gl_127[k]
                   + f_0 * il_352[k];

        t_218[k] = -gl_128[k]
                   + f_0 * il_353[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, gl_129, gl_130, gl_131, gl_132, \
                         gl_133, il_354, il_355, il_356, il_357, \
                         il_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = -gl_129[k]
                   + f_0 * il_354[k];

        t_220[k] = -gl_130[k]
                   + f_0 * il_355[k];

        t_221[k] = -gl_131[k]
                   + f_0 * il_356[k];

        t_222[k] = -gl_132[k]
                   + f_0 * il_357[k];

        t_223[k] = -gl_133[k]
                   + f_0 * il_358[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, t_229, t_230, gl_134, il_359, \
                         il_360, il_361, il_362, il_363, il_364, \
                         il_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = -gl_134[k]
                   + f_0 * il_359[k];

        t_225[k] = f_0 * il_360[k];

        t_226[k] = f_0 * il_361[k];

        t_227[k] = f_0 * il_362[k];

        t_228[k] = f_0 * il_363[k];

        t_229[k] = f_0 * il_364[k];

        t_230[k] = f_0 * il_365[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, t_236, t_237, t_238, il_366, \
                         il_367, il_368, il_369, il_370, il_371, il_372, \
                         il_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_0 * il_366[k];

        t_232[k] = f_0 * il_367[k];

        t_233[k] = f_0 * il_368[k];

        t_234[k] = f_0 * il_369[k];

        t_235[k] = f_0 * il_370[k];

        t_236[k] = f_0 * il_371[k];

        t_237[k] = f_0 * il_372[k];

        t_238[k] = f_0 * il_373[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, t_244, t_245, t_246, il_374, \
                         il_375, il_376, il_377, il_378, il_379, il_380, \
                         il_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_0 * il_374[k];

        t_240[k] = f_0 * il_375[k];

        t_241[k] = f_0 * il_376[k];

        t_242[k] = f_0 * il_377[k];

        t_243[k] = f_0 * il_378[k];

        t_244[k] = f_0 * il_379[k];

        t_245[k] = f_0 * il_380[k];

        t_246[k] = f_0 * il_381[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, t_252, t_253, t_254, il_382, \
                         il_383, il_384, il_385, il_386, il_387, il_388, \
                         il_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_0 * il_382[k];

        t_248[k] = f_0 * il_383[k];

        t_249[k] = f_0 * il_384[k];

        t_250[k] = f_0 * il_385[k];

        t_251[k] = f_0 * il_386[k];

        t_252[k] = f_0 * il_387[k];

        t_253[k] = f_0 * il_388[k];

        t_254[k] = f_0 * il_389[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, t_260, t_261, t_262, il_390, \
                         il_391, il_392, il_393, il_394, il_395, il_396, \
                         il_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_0 * il_390[k];

        t_256[k] = f_0 * il_391[k];

        t_257[k] = f_0 * il_392[k];

        t_258[k] = f_0 * il_393[k];

        t_259[k] = f_0 * il_394[k];

        t_260[k] = f_0 * il_395[k];

        t_261[k] = f_0 * il_396[k];

        t_262[k] = f_0 * il_397[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, t_268, t_269, il_398, il_399, \
                         il_400, il_401, il_402, il_403, il_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_0 * il_398[k];

        t_264[k] = f_0 * il_399[k];

        t_265[k] = f_0 * il_400[k];

        t_266[k] = f_0 * il_401[k];

        t_267[k] = f_0 * il_402[k];

        t_268[k] = f_0 * il_403[k];

        t_269[k] = f_0 * il_404[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, gl_135, gl_136, gl_137, gl_138, \
                         gl_139, il_450, il_451, il_452, il_453, \
                         il_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -3.0 * gl_135[k]
                   + f_0 * il_450[k];

        t_271[k] = -3.0 * gl_136[k]
                   + f_0 * il_451[k];

        t_272[k] = -3.0 * gl_137[k]
                   + f_0 * il_452[k];

        t_273[k] = -3.0 * gl_138[k]
                   + f_0 * il_453[k];

        t_274[k] = -3.0 * gl_139[k]
                   + f_0 * il_454[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, gl_140, gl_141, gl_142, gl_143, \
                         gl_144, il_455, il_456, il_457, il_458, \
                         il_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * gl_140[k]
                   + f_0 * il_455[k];

        t_276[k] = -3.0 * gl_141[k]
                   + f_0 * il_456[k];

        t_277[k] = -3.0 * gl_142[k]
                   + f_0 * il_457[k];

        t_278[k] = -3.0 * gl_143[k]
                   + f_0 * il_458[k];

        t_279[k] = -3.0 * gl_144[k]
                   + f_0 * il_459[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, gl_145, gl_146, gl_147, gl_148, \
                         gl_149, il_460, il_461, il_462, il_463, \
                         il_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * gl_145[k]
                   + f_0 * il_460[k];

        t_281[k] = -3.0 * gl_146[k]
                   + f_0 * il_461[k];

        t_282[k] = -3.0 * gl_147[k]
                   + f_0 * il_462[k];

        t_283[k] = -3.0 * gl_148[k]
                   + f_0 * il_463[k];

        t_284[k] = -3.0 * gl_149[k]
                   + f_0 * il_464[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, gl_150, gl_151, gl_152, gl_153, \
                         gl_154, il_465, il_466, il_467, il_468, \
                         il_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * gl_150[k]
                   + f_0 * il_465[k];

        t_286[k] = -3.0 * gl_151[k]
                   + f_0 * il_466[k];

        t_287[k] = -3.0 * gl_152[k]
                   + f_0 * il_467[k];

        t_288[k] = -3.0 * gl_153[k]
                   + f_0 * il_468[k];

        t_289[k] = -3.0 * gl_154[k]
                   + f_0 * il_469[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, gl_155, gl_156, gl_157, gl_158, \
                         gl_159, il_470, il_471, il_472, il_473, \
                         il_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * gl_155[k]
                   + f_0 * il_470[k];

        t_291[k] = -3.0 * gl_156[k]
                   + f_0 * il_471[k];

        t_292[k] = -3.0 * gl_157[k]
                   + f_0 * il_472[k];

        t_293[k] = -3.0 * gl_158[k]
                   + f_0 * il_473[k];

        t_294[k] = -3.0 * gl_159[k]
                   + f_0 * il_474[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, gl_160, gl_161, gl_162, gl_163, \
                         gl_164, il_475, il_476, il_477, il_478, \
                         il_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -3.0 * gl_160[k]
                   + f_0 * il_475[k];

        t_296[k] = -3.0 * gl_161[k]
                   + f_0 * il_476[k];

        t_297[k] = -3.0 * gl_162[k]
                   + f_0 * il_477[k];

        t_298[k] = -3.0 * gl_163[k]
                   + f_0 * il_478[k];

        t_299[k] = -3.0 * gl_164[k]
                   + f_0 * il_479[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, gl_165, gl_166, gl_167, gl_168, \
                         gl_169, il_480, il_481, il_482, il_483, \
                         il_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -3.0 * gl_165[k]
                   + f_0 * il_480[k];

        t_301[k] = -3.0 * gl_166[k]
                   + f_0 * il_481[k];

        t_302[k] = -3.0 * gl_167[k]
                   + f_0 * il_482[k];

        t_303[k] = -3.0 * gl_168[k]
                   + f_0 * il_483[k];

        t_304[k] = -3.0 * gl_169[k]
                   + f_0 * il_484[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, gl_170, gl_171, gl_172, gl_173, \
                         gl_174, il_485, il_486, il_487, il_488, \
                         il_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -3.0 * gl_170[k]
                   + f_0 * il_485[k];

        t_306[k] = -3.0 * gl_171[k]
                   + f_0 * il_486[k];

        t_307[k] = -3.0 * gl_172[k]
                   + f_0 * il_487[k];

        t_308[k] = -3.0 * gl_173[k]
                   + f_0 * il_488[k];

        t_309[k] = -3.0 * gl_174[k]
                   + f_0 * il_489[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, gl_175, gl_176, gl_177, gl_178, \
                         gl_179, il_490, il_491, il_492, il_493, \
                         il_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * gl_175[k]
                   + f_0 * il_490[k];

        t_311[k] = -3.0 * gl_176[k]
                   + f_0 * il_491[k];

        t_312[k] = -3.0 * gl_177[k]
                   + f_0 * il_492[k];

        t_313[k] = -3.0 * gl_178[k]
                   + f_0 * il_493[k];

        t_314[k] = -3.0 * gl_179[k]
                   + f_0 * il_494[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, gl_180, gl_181, gl_182, gl_183, \
                         gl_184, il_495, il_496, il_497, il_498, \
                         il_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -2.0 * gl_180[k]
                   + f_0 * il_495[k];

        t_316[k] = -2.0 * gl_181[k]
                   + f_0 * il_496[k];

        t_317[k] = -2.0 * gl_182[k]
                   + f_0 * il_497[k];

        t_318[k] = -2.0 * gl_183[k]
                   + f_0 * il_498[k];

        t_319[k] = -2.0 * gl_184[k]
                   + f_0 * il_499[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, gl_185, gl_186, gl_187, gl_188, \
                         gl_189, il_500, il_501, il_502, il_503, \
                         il_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -2.0 * gl_185[k]
                   + f_0 * il_500[k];

        t_321[k] = -2.0 * gl_186[k]
                   + f_0 * il_501[k];

        t_322[k] = -2.0 * gl_187[k]
                   + f_0 * il_502[k];

        t_323[k] = -2.0 * gl_188[k]
                   + f_0 * il_503[k];

        t_324[k] = -2.0 * gl_189[k]
                   + f_0 * il_504[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, gl_190, gl_191, gl_192, gl_193, \
                         gl_194, il_505, il_506, il_507, il_508, \
                         il_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -2.0 * gl_190[k]
                   + f_0 * il_505[k];

        t_326[k] = -2.0 * gl_191[k]
                   + f_0 * il_506[k];

        t_327[k] = -2.0 * gl_192[k]
                   + f_0 * il_507[k];

        t_328[k] = -2.0 * gl_193[k]
                   + f_0 * il_508[k];

        t_329[k] = -2.0 * gl_194[k]
                   + f_0 * il_509[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, gl_195, gl_196, gl_197, gl_198, \
                         gl_199, il_510, il_511, il_512, il_513, \
                         il_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -2.0 * gl_195[k]
                   + f_0 * il_510[k];

        t_331[k] = -2.0 * gl_196[k]
                   + f_0 * il_511[k];

        t_332[k] = -2.0 * gl_197[k]
                   + f_0 * il_512[k];

        t_333[k] = -2.0 * gl_198[k]
                   + f_0 * il_513[k];

        t_334[k] = -2.0 * gl_199[k]
                   + f_0 * il_514[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, gl_200, gl_201, gl_202, gl_203, \
                         gl_204, il_515, il_516, il_517, il_518, \
                         il_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -2.0 * gl_200[k]
                   + f_0 * il_515[k];

        t_336[k] = -2.0 * gl_201[k]
                   + f_0 * il_516[k];

        t_337[k] = -2.0 * gl_202[k]
                   + f_0 * il_517[k];

        t_338[k] = -2.0 * gl_203[k]
                   + f_0 * il_518[k];

        t_339[k] = -2.0 * gl_204[k]
                   + f_0 * il_519[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, gl_205, gl_206, gl_207, gl_208, \
                         gl_209, il_520, il_521, il_522, il_523, \
                         il_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * gl_205[k]
                   + f_0 * il_520[k];

        t_341[k] = -2.0 * gl_206[k]
                   + f_0 * il_521[k];

        t_342[k] = -2.0 * gl_207[k]
                   + f_0 * il_522[k];

        t_343[k] = -2.0 * gl_208[k]
                   + f_0 * il_523[k];

        t_344[k] = -2.0 * gl_209[k]
                   + f_0 * il_524[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_210 = buffer.data(gl + 210);
    const auto *gl_211 = buffer.data(gl + 211);
    const auto *gl_212 = buffer.data(gl + 212);
    const auto *gl_213 = buffer.data(gl + 213);
    const auto *gl_214 = buffer.data(gl + 214);
    const auto *gl_215 = buffer.data(gl + 215);
    const auto *gl_216 = buffer.data(gl + 216);
    const auto *gl_217 = buffer.data(gl + 217);
    const auto *gl_218 = buffer.data(gl + 218);
    const auto *gl_219 = buffer.data(gl + 219);
    const auto *gl_220 = buffer.data(gl + 220);
    const auto *gl_221 = buffer.data(gl + 221);
    const auto *gl_222 = buffer.data(gl + 222);
    const auto *gl_223 = buffer.data(gl + 223);
    const auto *gl_224 = buffer.data(gl + 224);
    const auto *gl_225 = buffer.data(gl + 225);
    const auto *gl_226 = buffer.data(gl + 226);
    const auto *gl_227 = buffer.data(gl + 227);
    const auto *gl_228 = buffer.data(gl + 228);
    const auto *gl_229 = buffer.data(gl + 229);
    const auto *gl_230 = buffer.data(gl + 230);
    const auto *gl_231 = buffer.data(gl + 231);
    const auto *gl_232 = buffer.data(gl + 232);
    const auto *gl_233 = buffer.data(gl + 233);
    const auto *gl_234 = buffer.data(gl + 234);
    const auto *gl_235 = buffer.data(gl + 235);
    const auto *gl_236 = buffer.data(gl + 236);
    const auto *gl_237 = buffer.data(gl + 237);
    const auto *gl_238 = buffer.data(gl + 238);
    const auto *gl_239 = buffer.data(gl + 239);
    const auto *gl_240 = buffer.data(gl + 240);
    const auto *gl_241 = buffer.data(gl + 241);
    const auto *gl_242 = buffer.data(gl + 242);
    const auto *gl_243 = buffer.data(gl + 243);
    const auto *gl_244 = buffer.data(gl + 244);
    const auto *gl_245 = buffer.data(gl + 245);
    const auto *gl_246 = buffer.data(gl + 246);
    const auto *gl_247 = buffer.data(gl + 247);
    const auto *gl_248 = buffer.data(gl + 248);
    const auto *gl_249 = buffer.data(gl + 249);
    const auto *gl_250 = buffer.data(gl + 250);
    const auto *gl_251 = buffer.data(gl + 251);
    const auto *gl_252 = buffer.data(gl + 252);
    const auto *gl_253 = buffer.data(gl + 253);
    const auto *gl_254 = buffer.data(gl + 254);
    const auto *gl_255 = buffer.data(gl + 255);
    const auto *gl_256 = buffer.data(gl + 256);
    const auto *gl_257 = buffer.data(gl + 257);
    const auto *gl_258 = buffer.data(gl + 258);
    const auto *gl_259 = buffer.data(gl + 259);
    const auto *gl_260 = buffer.data(gl + 260);
    const auto *gl_261 = buffer.data(gl + 261);
    const auto *gl_262 = buffer.data(gl + 262);
    const auto *gl_263 = buffer.data(gl + 263);
    const auto *gl_264 = buffer.data(gl + 264);
    const auto *gl_265 = buffer.data(gl + 265);
    const auto *gl_266 = buffer.data(gl + 266);
    const auto *gl_267 = buffer.data(gl + 267);
    const auto *gl_268 = buffer.data(gl + 268);
    const auto *gl_269 = buffer.data(gl + 269);
    const auto *gl_270 = buffer.data(gl + 270);
    const auto *gl_271 = buffer.data(gl + 271);
    const auto *gl_272 = buffer.data(gl + 272);
    const auto *gl_273 = buffer.data(gl + 273);
    const auto *gl_274 = buffer.data(gl + 274);
    const auto *gl_275 = buffer.data(gl + 275);
    const auto *gl_276 = buffer.data(gl + 276);
    const auto *gl_277 = buffer.data(gl + 277);
    const auto *gl_278 = buffer.data(gl + 278);
    const auto *gl_279 = buffer.data(gl + 279);
    const auto *gl_280 = buffer.data(gl + 280);
    const auto *gl_281 = buffer.data(gl + 281);
    const auto *gl_282 = buffer.data(gl + 282);
    const auto *gl_283 = buffer.data(gl + 283);
    const auto *gl_284 = buffer.data(gl + 284);
    const auto *gl_285 = buffer.data(gl + 285);
    const auto *gl_286 = buffer.data(gl + 286);
    const auto *gl_287 = buffer.data(gl + 287);
    const auto *gl_288 = buffer.data(gl + 288);
    const auto *gl_289 = buffer.data(gl + 289);
    const auto *gl_290 = buffer.data(gl + 290);
    const auto *gl_291 = buffer.data(gl + 291);
    const auto *gl_292 = buffer.data(gl + 292);
    const auto *gl_293 = buffer.data(gl + 293);
    const auto *gl_294 = buffer.data(gl + 294);
    const auto *gl_295 = buffer.data(gl + 295);
    const auto *gl_296 = buffer.data(gl + 296);
    const auto *gl_297 = buffer.data(gl + 297);
    const auto *gl_298 = buffer.data(gl + 298);
    const auto *gl_299 = buffer.data(gl + 299);
    const auto *gl_300 = buffer.data(gl + 300);
    const auto *gl_301 = buffer.data(gl + 301);
    const auto *gl_302 = buffer.data(gl + 302);
    const auto *gl_303 = buffer.data(gl + 303);
    const auto *gl_304 = buffer.data(gl + 304);
    const auto *gl_305 = buffer.data(gl + 305);
    const auto *gl_306 = buffer.data(gl + 306);
    const auto *gl_307 = buffer.data(gl + 307);
    const auto *gl_308 = buffer.data(gl + 308);
    const auto *gl_309 = buffer.data(gl + 309);
    const auto *gl_310 = buffer.data(gl + 310);
    const auto *gl_311 = buffer.data(gl + 311);
    const auto *gl_312 = buffer.data(gl + 312);
    const auto *gl_313 = buffer.data(gl + 313);
    const auto *gl_314 = buffer.data(gl + 314);
    const auto *gl_315 = buffer.data(gl + 315);
    const auto *gl_316 = buffer.data(gl + 316);
    const auto *gl_317 = buffer.data(gl + 317);
    const auto *gl_318 = buffer.data(gl + 318);
    const auto *gl_319 = buffer.data(gl + 319);
    const auto *gl_320 = buffer.data(gl + 320);
    const auto *gl_321 = buffer.data(gl + 321);
    const auto *gl_322 = buffer.data(gl + 322);
    const auto *gl_323 = buffer.data(gl + 323);
    const auto *gl_324 = buffer.data(gl + 324);
    const auto *gl_325 = buffer.data(gl + 325);
    const auto *gl_326 = buffer.data(gl + 326);

    const auto *il_525 = buffer.data(il + 525);
    const auto *il_526 = buffer.data(il + 526);
    const auto *il_527 = buffer.data(il + 527);
    const auto *il_528 = buffer.data(il + 528);
    const auto *il_529 = buffer.data(il + 529);
    const auto *il_530 = buffer.data(il + 530);
    const auto *il_531 = buffer.data(il + 531);
    const auto *il_532 = buffer.data(il + 532);
    const auto *il_533 = buffer.data(il + 533);
    const auto *il_534 = buffer.data(il + 534);
    const auto *il_535 = buffer.data(il + 535);
    const auto *il_536 = buffer.data(il + 536);
    const auto *il_537 = buffer.data(il + 537);
    const auto *il_538 = buffer.data(il + 538);
    const auto *il_539 = buffer.data(il + 539);
    const auto *il_540 = buffer.data(il + 540);
    const auto *il_541 = buffer.data(il + 541);
    const auto *il_542 = buffer.data(il + 542);
    const auto *il_543 = buffer.data(il + 543);
    const auto *il_544 = buffer.data(il + 544);
    const auto *il_545 = buffer.data(il + 545);
    const auto *il_546 = buffer.data(il + 546);
    const auto *il_547 = buffer.data(il + 547);
    const auto *il_548 = buffer.data(il + 548);
    const auto *il_549 = buffer.data(il + 549);
    const auto *il_550 = buffer.data(il + 550);
    const auto *il_551 = buffer.data(il + 551);
    const auto *il_552 = buffer.data(il + 552);
    const auto *il_553 = buffer.data(il + 553);
    const auto *il_554 = buffer.data(il + 554);
    const auto *il_555 = buffer.data(il + 555);
    const auto *il_556 = buffer.data(il + 556);
    const auto *il_557 = buffer.data(il + 557);
    const auto *il_558 = buffer.data(il + 558);
    const auto *il_559 = buffer.data(il + 559);
    const auto *il_560 = buffer.data(il + 560);
    const auto *il_561 = buffer.data(il + 561);
    const auto *il_562 = buffer.data(il + 562);
    const auto *il_563 = buffer.data(il + 563);
    const auto *il_564 = buffer.data(il + 564);
    const auto *il_565 = buffer.data(il + 565);
    const auto *il_566 = buffer.data(il + 566);
    const auto *il_567 = buffer.data(il + 567);
    const auto *il_568 = buffer.data(il + 568);
    const auto *il_569 = buffer.data(il + 569);
    const auto *il_570 = buffer.data(il + 570);
    const auto *il_571 = buffer.data(il + 571);
    const auto *il_572 = buffer.data(il + 572);
    const auto *il_573 = buffer.data(il + 573);
    const auto *il_574 = buffer.data(il + 574);
    const auto *il_575 = buffer.data(il + 575);
    const auto *il_576 = buffer.data(il + 576);
    const auto *il_577 = buffer.data(il + 577);
    const auto *il_578 = buffer.data(il + 578);
    const auto *il_579 = buffer.data(il + 579);
    const auto *il_580 = buffer.data(il + 580);
    const auto *il_581 = buffer.data(il + 581);
    const auto *il_582 = buffer.data(il + 582);
    const auto *il_583 = buffer.data(il + 583);
    const auto *il_584 = buffer.data(il + 584);
    const auto *il_585 = buffer.data(il + 585);
    const auto *il_586 = buffer.data(il + 586);
    const auto *il_587 = buffer.data(il + 587);
    const auto *il_588 = buffer.data(il + 588);
    const auto *il_589 = buffer.data(il + 589);
    const auto *il_590 = buffer.data(il + 590);
    const auto *il_591 = buffer.data(il + 591);
    const auto *il_592 = buffer.data(il + 592);
    const auto *il_593 = buffer.data(il + 593);
    const auto *il_594 = buffer.data(il + 594);
    const auto *il_595 = buffer.data(il + 595);
    const auto *il_596 = buffer.data(il + 596);
    const auto *il_597 = buffer.data(il + 597);
    const auto *il_598 = buffer.data(il + 598);
    const auto *il_599 = buffer.data(il + 599);
    const auto *il_600 = buffer.data(il + 600);
    const auto *il_601 = buffer.data(il + 601);
    const auto *il_602 = buffer.data(il + 602);
    const auto *il_603 = buffer.data(il + 603);
    const auto *il_604 = buffer.data(il + 604);
    const auto *il_605 = buffer.data(il + 605);
    const auto *il_606 = buffer.data(il + 606);
    const auto *il_607 = buffer.data(il + 607);
    const auto *il_608 = buffer.data(il + 608);
    const auto *il_609 = buffer.data(il + 609);
    const auto *il_610 = buffer.data(il + 610);
    const auto *il_611 = buffer.data(il + 611);
    const auto *il_612 = buffer.data(il + 612);
    const auto *il_613 = buffer.data(il + 613);
    const auto *il_614 = buffer.data(il + 614);
    const auto *il_615 = buffer.data(il + 615);
    const auto *il_616 = buffer.data(il + 616);
    const auto *il_617 = buffer.data(il + 617);
    const auto *il_618 = buffer.data(il + 618);
    const auto *il_619 = buffer.data(il + 619);
    const auto *il_620 = buffer.data(il + 620);
    const auto *il_621 = buffer.data(il + 621);
    const auto *il_622 = buffer.data(il + 622);
    const auto *il_623 = buffer.data(il + 623);
    const auto *il_624 = buffer.data(il + 624);
    const auto *il_625 = buffer.data(il + 625);
    const auto *il_626 = buffer.data(il + 626);
    const auto *il_627 = buffer.data(il + 627);
    const auto *il_628 = buffer.data(il + 628);
    const auto *il_629 = buffer.data(il + 629);
    const auto *il_675 = buffer.data(il + 675);
    const auto *il_676 = buffer.data(il + 676);
    const auto *il_677 = buffer.data(il + 677);
    const auto *il_678 = buffer.data(il + 678);
    const auto *il_679 = buffer.data(il + 679);
    const auto *il_680 = buffer.data(il + 680);
    const auto *il_681 = buffer.data(il + 681);
    const auto *il_682 = buffer.data(il + 682);
    const auto *il_683 = buffer.data(il + 683);
    const auto *il_684 = buffer.data(il + 684);
    const auto *il_685 = buffer.data(il + 685);
    const auto *il_686 = buffer.data(il + 686);
    const auto *il_687 = buffer.data(il + 687);
    const auto *il_688 = buffer.data(il + 688);
    const auto *il_689 = buffer.data(il + 689);
    const auto *il_690 = buffer.data(il + 690);
    const auto *il_691 = buffer.data(il + 691);
    const auto *il_692 = buffer.data(il + 692);
    const auto *il_693 = buffer.data(il + 693);
    const auto *il_694 = buffer.data(il + 694);
    const auto *il_695 = buffer.data(il + 695);
    const auto *il_696 = buffer.data(il + 696);
    const auto *il_697 = buffer.data(il + 697);
    const auto *il_698 = buffer.data(il + 698);
    const auto *il_699 = buffer.data(il + 699);
    const auto *il_700 = buffer.data(il + 700);
    const auto *il_701 = buffer.data(il + 701);
    const auto *il_702 = buffer.data(il + 702);
    const auto *il_703 = buffer.data(il + 703);
    const auto *il_704 = buffer.data(il + 704);
    const auto *il_705 = buffer.data(il + 705);
    const auto *il_706 = buffer.data(il + 706);
    const auto *il_707 = buffer.data(il + 707);
    const auto *il_708 = buffer.data(il + 708);
    const auto *il_709 = buffer.data(il + 709);
    const auto *il_710 = buffer.data(il + 710);
    const auto *il_711 = buffer.data(il + 711);
    const auto *il_712 = buffer.data(il + 712);
    const auto *il_713 = buffer.data(il + 713);
    const auto *il_714 = buffer.data(il + 714);
    const auto *il_715 = buffer.data(il + 715);
    const auto *il_716 = buffer.data(il + 716);
    const auto *il_717 = buffer.data(il + 717);
    const auto *il_718 = buffer.data(il + 718);
    const auto *il_719 = buffer.data(il + 719);
    const auto *il_720 = buffer.data(il + 720);
    const auto *il_721 = buffer.data(il + 721);
    const auto *il_722 = buffer.data(il + 722);
    const auto *il_723 = buffer.data(il + 723);
    const auto *il_724 = buffer.data(il + 724);
    const auto *il_725 = buffer.data(il + 725);
    const auto *il_726 = buffer.data(il + 726);
    const auto *il_727 = buffer.data(il + 727);
    const auto *il_728 = buffer.data(il + 728);
    const auto *il_729 = buffer.data(il + 729);
    const auto *il_730 = buffer.data(il + 730);
    const auto *il_731 = buffer.data(il + 731);

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, gl_210, gl_211, gl_212, gl_213, \
                         gl_214, il_525, il_526, il_527, il_528, \
                         il_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * gl_210[k]
                   + f_0 * il_525[k];

        t_346[k] = -2.0 * gl_211[k]
                   + f_0 * il_526[k];

        t_347[k] = -2.0 * gl_212[k]
                   + f_0 * il_527[k];

        t_348[k] = -2.0 * gl_213[k]
                   + f_0 * il_528[k];

        t_349[k] = -2.0 * gl_214[k]
                   + f_0 * il_529[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, gl_215, gl_216, gl_217, gl_218, \
                         gl_219, il_530, il_531, il_532, il_533, \
                         il_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * gl_215[k]
                   + f_0 * il_530[k];

        t_351[k] = -2.0 * gl_216[k]
                   + f_0 * il_531[k];

        t_352[k] = -2.0 * gl_217[k]
                   + f_0 * il_532[k];

        t_353[k] = -2.0 * gl_218[k]
                   + f_0 * il_533[k];

        t_354[k] = -2.0 * gl_219[k]
                   + f_0 * il_534[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, gl_220, gl_221, gl_222, gl_223, \
                         gl_224, il_535, il_536, il_537, il_538, \
                         il_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * gl_220[k]
                   + f_0 * il_535[k];

        t_356[k] = -2.0 * gl_221[k]
                   + f_0 * il_536[k];

        t_357[k] = -2.0 * gl_222[k]
                   + f_0 * il_537[k];

        t_358[k] = -2.0 * gl_223[k]
                   + f_0 * il_538[k];

        t_359[k] = -2.0 * gl_224[k]
                   + f_0 * il_539[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, gl_225, gl_226, gl_227, gl_228, \
                         gl_229, il_540, il_541, il_542, il_543, \
                         il_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -gl_225[k]
                   + f_0 * il_540[k];

        t_361[k] = -gl_226[k]
                   + f_0 * il_541[k];

        t_362[k] = -gl_227[k]
                   + f_0 * il_542[k];

        t_363[k] = -gl_228[k]
                   + f_0 * il_543[k];

        t_364[k] = -gl_229[k]
                   + f_0 * il_544[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, gl_230, gl_231, gl_232, gl_233, \
                         gl_234, il_545, il_546, il_547, il_548, \
                         il_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -gl_230[k]
                   + f_0 * il_545[k];

        t_366[k] = -gl_231[k]
                   + f_0 * il_546[k];

        t_367[k] = -gl_232[k]
                   + f_0 * il_547[k];

        t_368[k] = -gl_233[k]
                   + f_0 * il_548[k];

        t_369[k] = -gl_234[k]
                   + f_0 * il_549[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, gl_235, gl_236, gl_237, gl_238, \
                         gl_239, il_550, il_551, il_552, il_553, \
                         il_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -gl_235[k]
                   + f_0 * il_550[k];

        t_371[k] = -gl_236[k]
                   + f_0 * il_551[k];

        t_372[k] = -gl_237[k]
                   + f_0 * il_552[k];

        t_373[k] = -gl_238[k]
                   + f_0 * il_553[k];

        t_374[k] = -gl_239[k]
                   + f_0 * il_554[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, gl_240, gl_241, gl_242, gl_243, \
                         gl_244, il_555, il_556, il_557, il_558, \
                         il_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -gl_240[k]
                   + f_0 * il_555[k];

        t_376[k] = -gl_241[k]
                   + f_0 * il_556[k];

        t_377[k] = -gl_242[k]
                   + f_0 * il_557[k];

        t_378[k] = -gl_243[k]
                   + f_0 * il_558[k];

        t_379[k] = -gl_244[k]
                   + f_0 * il_559[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, gl_245, gl_246, gl_247, gl_248, \
                         gl_249, il_560, il_561, il_562, il_563, \
                         il_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -gl_245[k]
                   + f_0 * il_560[k];

        t_381[k] = -gl_246[k]
                   + f_0 * il_561[k];

        t_382[k] = -gl_247[k]
                   + f_0 * il_562[k];

        t_383[k] = -gl_248[k]
                   + f_0 * il_563[k];

        t_384[k] = -gl_249[k]
                   + f_0 * il_564[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, gl_250, gl_251, gl_252, gl_253, \
                         gl_254, il_565, il_566, il_567, il_568, \
                         il_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -gl_250[k]
                   + f_0 * il_565[k];

        t_386[k] = -gl_251[k]
                   + f_0 * il_566[k];

        t_387[k] = -gl_252[k]
                   + f_0 * il_567[k];

        t_388[k] = -gl_253[k]
                   + f_0 * il_568[k];

        t_389[k] = -gl_254[k]
                   + f_0 * il_569[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, gl_255, gl_256, gl_257, gl_258, \
                         gl_259, il_570, il_571, il_572, il_573, \
                         il_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -gl_255[k]
                   + f_0 * il_570[k];

        t_391[k] = -gl_256[k]
                   + f_0 * il_571[k];

        t_392[k] = -gl_257[k]
                   + f_0 * il_572[k];

        t_393[k] = -gl_258[k]
                   + f_0 * il_573[k];

        t_394[k] = -gl_259[k]
                   + f_0 * il_574[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, gl_260, gl_261, gl_262, gl_263, \
                         gl_264, il_575, il_576, il_577, il_578, \
                         il_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -gl_260[k]
                   + f_0 * il_575[k];

        t_396[k] = -gl_261[k]
                   + f_0 * il_576[k];

        t_397[k] = -gl_262[k]
                   + f_0 * il_577[k];

        t_398[k] = -gl_263[k]
                   + f_0 * il_578[k];

        t_399[k] = -gl_264[k]
                   + f_0 * il_579[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, gl_265, gl_266, gl_267, gl_268, \
                         gl_269, il_580, il_581, il_582, il_583, \
                         il_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -gl_265[k]
                   + f_0 * il_580[k];

        t_401[k] = -gl_266[k]
                   + f_0 * il_581[k];

        t_402[k] = -gl_267[k]
                   + f_0 * il_582[k];

        t_403[k] = -gl_268[k]
                   + f_0 * il_583[k];

        t_404[k] = -gl_269[k]
                   + f_0 * il_584[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, il_585, \
                         il_586, il_587, il_588, il_589, il_590, il_591, \
                         il_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * il_585[k];

        t_406[k] = f_0 * il_586[k];

        t_407[k] = f_0 * il_587[k];

        t_408[k] = f_0 * il_588[k];

        t_409[k] = f_0 * il_589[k];

        t_410[k] = f_0 * il_590[k];

        t_411[k] = f_0 * il_591[k];

        t_412[k] = f_0 * il_592[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, t_420, il_593, \
                         il_594, il_595, il_596, il_597, il_598, il_599, \
                         il_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * il_593[k];

        t_414[k] = f_0 * il_594[k];

        t_415[k] = f_0 * il_595[k];

        t_416[k] = f_0 * il_596[k];

        t_417[k] = f_0 * il_597[k];

        t_418[k] = f_0 * il_598[k];

        t_419[k] = f_0 * il_599[k];

        t_420[k] = f_0 * il_600[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, t_425, t_426, t_427, t_428, il_601, \
                         il_602, il_603, il_604, il_605, il_606, il_607, \
                         il_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_0 * il_601[k];

        t_422[k] = f_0 * il_602[k];

        t_423[k] = f_0 * il_603[k];

        t_424[k] = f_0 * il_604[k];

        t_425[k] = f_0 * il_605[k];

        t_426[k] = f_0 * il_606[k];

        t_427[k] = f_0 * il_607[k];

        t_428[k] = f_0 * il_608[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, t_434, t_435, t_436, il_609, \
                         il_610, il_611, il_612, il_613, il_614, il_615, \
                         il_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_0 * il_609[k];

        t_430[k] = f_0 * il_610[k];

        t_431[k] = f_0 * il_611[k];

        t_432[k] = f_0 * il_612[k];

        t_433[k] = f_0 * il_613[k];

        t_434[k] = f_0 * il_614[k];

        t_435[k] = f_0 * il_615[k];

        t_436[k] = f_0 * il_616[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, t_441, t_442, t_443, t_444, il_617, \
                         il_618, il_619, il_620, il_621, il_622, il_623, \
                         il_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_0 * il_617[k];

        t_438[k] = f_0 * il_618[k];

        t_439[k] = f_0 * il_619[k];

        t_440[k] = f_0 * il_620[k];

        t_441[k] = f_0 * il_621[k];

        t_442[k] = f_0 * il_622[k];

        t_443[k] = f_0 * il_623[k];

        t_444[k] = f_0 * il_624[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, t_450, t_451, gl_270, gl_271, \
                         il_625, il_626, il_627, il_628, il_629, il_675, \
                         il_676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_0 * il_625[k];

        t_446[k] = f_0 * il_626[k];

        t_447[k] = f_0 * il_627[k];

        t_448[k] = f_0 * il_628[k];

        t_449[k] = f_0 * il_629[k];

        t_450[k] = -4.0 * gl_270[k]
                   + f_0 * il_675[k];

        t_451[k] = -4.0 * gl_271[k]
                   + f_0 * il_676[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, gl_272, gl_273, gl_274, gl_275, \
                         gl_276, il_677, il_678, il_679, il_680, \
                         il_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = -4.0 * gl_272[k]
                   + f_0 * il_677[k];

        t_453[k] = -4.0 * gl_273[k]
                   + f_0 * il_678[k];

        t_454[k] = -4.0 * gl_274[k]
                   + f_0 * il_679[k];

        t_455[k] = -4.0 * gl_275[k]
                   + f_0 * il_680[k];

        t_456[k] = -4.0 * gl_276[k]
                   + f_0 * il_681[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, t_461, gl_277, gl_278, gl_279, gl_280, \
                         gl_281, il_682, il_683, il_684, il_685, \
                         il_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = -4.0 * gl_277[k]
                   + f_0 * il_682[k];

        t_458[k] = -4.0 * gl_278[k]
                   + f_0 * il_683[k];

        t_459[k] = -4.0 * gl_279[k]
                   + f_0 * il_684[k];

        t_460[k] = -4.0 * gl_280[k]
                   + f_0 * il_685[k];

        t_461[k] = -4.0 * gl_281[k]
                   + f_0 * il_686[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, gl_282, gl_283, gl_284, gl_285, \
                         gl_286, il_687, il_688, il_689, il_690, \
                         il_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = -4.0 * gl_282[k]
                   + f_0 * il_687[k];

        t_463[k] = -4.0 * gl_283[k]
                   + f_0 * il_688[k];

        t_464[k] = -4.0 * gl_284[k]
                   + f_0 * il_689[k];

        t_465[k] = -4.0 * gl_285[k]
                   + f_0 * il_690[k];

        t_466[k] = -4.0 * gl_286[k]
                   + f_0 * il_691[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, gl_287, gl_288, gl_289, gl_290, \
                         gl_291, il_692, il_693, il_694, il_695, \
                         il_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = -4.0 * gl_287[k]
                   + f_0 * il_692[k];

        t_468[k] = -4.0 * gl_288[k]
                   + f_0 * il_693[k];

        t_469[k] = -4.0 * gl_289[k]
                   + f_0 * il_694[k];

        t_470[k] = -4.0 * gl_290[k]
                   + f_0 * il_695[k];

        t_471[k] = -4.0 * gl_291[k]
                   + f_0 * il_696[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, gl_292, gl_293, gl_294, gl_295, \
                         gl_296, il_697, il_698, il_699, il_700, \
                         il_701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = -4.0 * gl_292[k]
                   + f_0 * il_697[k];

        t_473[k] = -4.0 * gl_293[k]
                   + f_0 * il_698[k];

        t_474[k] = -4.0 * gl_294[k]
                   + f_0 * il_699[k];

        t_475[k] = -4.0 * gl_295[k]
                   + f_0 * il_700[k];

        t_476[k] = -4.0 * gl_296[k]
                   + f_0 * il_701[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, gl_297, gl_298, gl_299, gl_300, \
                         gl_301, il_702, il_703, il_704, il_705, \
                         il_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = -4.0 * gl_297[k]
                   + f_0 * il_702[k];

        t_478[k] = -4.0 * gl_298[k]
                   + f_0 * il_703[k];

        t_479[k] = -4.0 * gl_299[k]
                   + f_0 * il_704[k];

        t_480[k] = -4.0 * gl_300[k]
                   + f_0 * il_705[k];

        t_481[k] = -4.0 * gl_301[k]
                   + f_0 * il_706[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, t_486, gl_302, gl_303, gl_304, gl_305, \
                         gl_306, il_707, il_708, il_709, il_710, \
                         il_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = -4.0 * gl_302[k]
                   + f_0 * il_707[k];

        t_483[k] = -4.0 * gl_303[k]
                   + f_0 * il_708[k];

        t_484[k] = -4.0 * gl_304[k]
                   + f_0 * il_709[k];

        t_485[k] = -4.0 * gl_305[k]
                   + f_0 * il_710[k];

        t_486[k] = -4.0 * gl_306[k]
                   + f_0 * il_711[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, gl_307, gl_308, gl_309, gl_310, \
                         gl_311, il_712, il_713, il_714, il_715, \
                         il_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = -4.0 * gl_307[k]
                   + f_0 * il_712[k];

        t_488[k] = -4.0 * gl_308[k]
                   + f_0 * il_713[k];

        t_489[k] = -4.0 * gl_309[k]
                   + f_0 * il_714[k];

        t_490[k] = -4.0 * gl_310[k]
                   + f_0 * il_715[k];

        t_491[k] = -4.0 * gl_311[k]
                   + f_0 * il_716[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, gl_312, gl_313, gl_314, gl_315, \
                         gl_316, il_717, il_718, il_719, il_720, \
                         il_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = -4.0 * gl_312[k]
                   + f_0 * il_717[k];

        t_493[k] = -4.0 * gl_313[k]
                   + f_0 * il_718[k];

        t_494[k] = -4.0 * gl_314[k]
                   + f_0 * il_719[k];

        t_495[k] = -3.0 * gl_315[k]
                   + f_0 * il_720[k];

        t_496[k] = -3.0 * gl_316[k]
                   + f_0 * il_721[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, t_500, t_501, gl_317, gl_318, gl_319, gl_320, \
                         gl_321, il_722, il_723, il_724, il_725, \
                         il_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = -3.0 * gl_317[k]
                   + f_0 * il_722[k];

        t_498[k] = -3.0 * gl_318[k]
                   + f_0 * il_723[k];

        t_499[k] = -3.0 * gl_319[k]
                   + f_0 * il_724[k];

        t_500[k] = -3.0 * gl_320[k]
                   + f_0 * il_725[k];

        t_501[k] = -3.0 * gl_321[k]
                   + f_0 * il_726[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, gl_322, gl_323, gl_324, gl_325, \
                         gl_326, il_727, il_728, il_729, il_730, \
                         il_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = -3.0 * gl_322[k]
                   + f_0 * il_727[k];

        t_503[k] = -3.0 * gl_323[k]
                   + f_0 * il_728[k];

        t_504[k] = -3.0 * gl_324[k]
                   + f_0 * il_729[k];

        t_505[k] = -3.0 * gl_325[k]
                   + f_0 * il_730[k];

        t_506[k] = -3.0 * gl_326[k]
                   + f_0 * il_731[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_327 = buffer.data(gl + 327);
    const auto *gl_328 = buffer.data(gl + 328);
    const auto *gl_329 = buffer.data(gl + 329);
    const auto *gl_330 = buffer.data(gl + 330);
    const auto *gl_331 = buffer.data(gl + 331);
    const auto *gl_332 = buffer.data(gl + 332);
    const auto *gl_333 = buffer.data(gl + 333);
    const auto *gl_334 = buffer.data(gl + 334);
    const auto *gl_335 = buffer.data(gl + 335);
    const auto *gl_336 = buffer.data(gl + 336);
    const auto *gl_337 = buffer.data(gl + 337);
    const auto *gl_338 = buffer.data(gl + 338);
    const auto *gl_339 = buffer.data(gl + 339);
    const auto *gl_340 = buffer.data(gl + 340);
    const auto *gl_341 = buffer.data(gl + 341);
    const auto *gl_342 = buffer.data(gl + 342);
    const auto *gl_343 = buffer.data(gl + 343);
    const auto *gl_344 = buffer.data(gl + 344);
    const auto *gl_345 = buffer.data(gl + 345);
    const auto *gl_346 = buffer.data(gl + 346);
    const auto *gl_347 = buffer.data(gl + 347);
    const auto *gl_348 = buffer.data(gl + 348);
    const auto *gl_349 = buffer.data(gl + 349);
    const auto *gl_350 = buffer.data(gl + 350);
    const auto *gl_351 = buffer.data(gl + 351);
    const auto *gl_352 = buffer.data(gl + 352);
    const auto *gl_353 = buffer.data(gl + 353);
    const auto *gl_354 = buffer.data(gl + 354);
    const auto *gl_355 = buffer.data(gl + 355);
    const auto *gl_356 = buffer.data(gl + 356);
    const auto *gl_357 = buffer.data(gl + 357);
    const auto *gl_358 = buffer.data(gl + 358);
    const auto *gl_359 = buffer.data(gl + 359);
    const auto *gl_360 = buffer.data(gl + 360);
    const auto *gl_361 = buffer.data(gl + 361);
    const auto *gl_362 = buffer.data(gl + 362);
    const auto *gl_363 = buffer.data(gl + 363);
    const auto *gl_364 = buffer.data(gl + 364);
    const auto *gl_365 = buffer.data(gl + 365);
    const auto *gl_366 = buffer.data(gl + 366);
    const auto *gl_367 = buffer.data(gl + 367);
    const auto *gl_368 = buffer.data(gl + 368);
    const auto *gl_369 = buffer.data(gl + 369);
    const auto *gl_370 = buffer.data(gl + 370);
    const auto *gl_371 = buffer.data(gl + 371);
    const auto *gl_372 = buffer.data(gl + 372);
    const auto *gl_373 = buffer.data(gl + 373);
    const auto *gl_374 = buffer.data(gl + 374);
    const auto *gl_375 = buffer.data(gl + 375);
    const auto *gl_376 = buffer.data(gl + 376);
    const auto *gl_377 = buffer.data(gl + 377);
    const auto *gl_378 = buffer.data(gl + 378);
    const auto *gl_379 = buffer.data(gl + 379);
    const auto *gl_380 = buffer.data(gl + 380);
    const auto *gl_381 = buffer.data(gl + 381);
    const auto *gl_382 = buffer.data(gl + 382);
    const auto *gl_383 = buffer.data(gl + 383);
    const auto *gl_384 = buffer.data(gl + 384);
    const auto *gl_385 = buffer.data(gl + 385);
    const auto *gl_386 = buffer.data(gl + 386);
    const auto *gl_387 = buffer.data(gl + 387);
    const auto *gl_388 = buffer.data(gl + 388);
    const auto *gl_389 = buffer.data(gl + 389);
    const auto *gl_390 = buffer.data(gl + 390);
    const auto *gl_391 = buffer.data(gl + 391);
    const auto *gl_392 = buffer.data(gl + 392);
    const auto *gl_393 = buffer.data(gl + 393);
    const auto *gl_394 = buffer.data(gl + 394);
    const auto *gl_395 = buffer.data(gl + 395);
    const auto *gl_396 = buffer.data(gl + 396);
    const auto *gl_397 = buffer.data(gl + 397);
    const auto *gl_398 = buffer.data(gl + 398);
    const auto *gl_399 = buffer.data(gl + 399);
    const auto *gl_400 = buffer.data(gl + 400);
    const auto *gl_401 = buffer.data(gl + 401);
    const auto *gl_402 = buffer.data(gl + 402);
    const auto *gl_403 = buffer.data(gl + 403);
    const auto *gl_404 = buffer.data(gl + 404);
    const auto *gl_405 = buffer.data(gl + 405);
    const auto *gl_406 = buffer.data(gl + 406);
    const auto *gl_407 = buffer.data(gl + 407);
    const auto *gl_408 = buffer.data(gl + 408);
    const auto *gl_409 = buffer.data(gl + 409);
    const auto *gl_410 = buffer.data(gl + 410);
    const auto *gl_411 = buffer.data(gl + 411);
    const auto *gl_412 = buffer.data(gl + 412);
    const auto *gl_413 = buffer.data(gl + 413);
    const auto *gl_414 = buffer.data(gl + 414);
    const auto *gl_415 = buffer.data(gl + 415);
    const auto *gl_416 = buffer.data(gl + 416);
    const auto *gl_417 = buffer.data(gl + 417);
    const auto *gl_418 = buffer.data(gl + 418);
    const auto *gl_419 = buffer.data(gl + 419);
    const auto *gl_420 = buffer.data(gl + 420);
    const auto *gl_421 = buffer.data(gl + 421);
    const auto *gl_422 = buffer.data(gl + 422);
    const auto *gl_423 = buffer.data(gl + 423);
    const auto *gl_424 = buffer.data(gl + 424);
    const auto *gl_425 = buffer.data(gl + 425);
    const auto *gl_426 = buffer.data(gl + 426);
    const auto *gl_427 = buffer.data(gl + 427);
    const auto *gl_428 = buffer.data(gl + 428);
    const auto *gl_429 = buffer.data(gl + 429);
    const auto *gl_430 = buffer.data(gl + 430);
    const auto *gl_431 = buffer.data(gl + 431);
    const auto *gl_432 = buffer.data(gl + 432);
    const auto *gl_433 = buffer.data(gl + 433);
    const auto *gl_434 = buffer.data(gl + 434);
    const auto *gl_435 = buffer.data(gl + 435);
    const auto *gl_436 = buffer.data(gl + 436);
    const auto *gl_437 = buffer.data(gl + 437);
    const auto *gl_438 = buffer.data(gl + 438);
    const auto *gl_439 = buffer.data(gl + 439);
    const auto *gl_440 = buffer.data(gl + 440);
    const auto *gl_441 = buffer.data(gl + 441);
    const auto *gl_442 = buffer.data(gl + 442);
    const auto *gl_443 = buffer.data(gl + 443);
    const auto *gl_444 = buffer.data(gl + 444);
    const auto *gl_445 = buffer.data(gl + 445);
    const auto *gl_446 = buffer.data(gl + 446);
    const auto *gl_447 = buffer.data(gl + 447);
    const auto *gl_448 = buffer.data(gl + 448);
    const auto *gl_449 = buffer.data(gl + 449);

    const auto *il_732 = buffer.data(il + 732);
    const auto *il_733 = buffer.data(il + 733);
    const auto *il_734 = buffer.data(il + 734);
    const auto *il_735 = buffer.data(il + 735);
    const auto *il_736 = buffer.data(il + 736);
    const auto *il_737 = buffer.data(il + 737);
    const auto *il_738 = buffer.data(il + 738);
    const auto *il_739 = buffer.data(il + 739);
    const auto *il_740 = buffer.data(il + 740);
    const auto *il_741 = buffer.data(il + 741);
    const auto *il_742 = buffer.data(il + 742);
    const auto *il_743 = buffer.data(il + 743);
    const auto *il_744 = buffer.data(il + 744);
    const auto *il_745 = buffer.data(il + 745);
    const auto *il_746 = buffer.data(il + 746);
    const auto *il_747 = buffer.data(il + 747);
    const auto *il_748 = buffer.data(il + 748);
    const auto *il_749 = buffer.data(il + 749);
    const auto *il_750 = buffer.data(il + 750);
    const auto *il_751 = buffer.data(il + 751);
    const auto *il_752 = buffer.data(il + 752);
    const auto *il_753 = buffer.data(il + 753);
    const auto *il_754 = buffer.data(il + 754);
    const auto *il_755 = buffer.data(il + 755);
    const auto *il_756 = buffer.data(il + 756);
    const auto *il_757 = buffer.data(il + 757);
    const auto *il_758 = buffer.data(il + 758);
    const auto *il_759 = buffer.data(il + 759);
    const auto *il_760 = buffer.data(il + 760);
    const auto *il_761 = buffer.data(il + 761);
    const auto *il_762 = buffer.data(il + 762);
    const auto *il_763 = buffer.data(il + 763);
    const auto *il_764 = buffer.data(il + 764);
    const auto *il_765 = buffer.data(il + 765);
    const auto *il_766 = buffer.data(il + 766);
    const auto *il_767 = buffer.data(il + 767);
    const auto *il_768 = buffer.data(il + 768);
    const auto *il_769 = buffer.data(il + 769);
    const auto *il_770 = buffer.data(il + 770);
    const auto *il_771 = buffer.data(il + 771);
    const auto *il_772 = buffer.data(il + 772);
    const auto *il_773 = buffer.data(il + 773);
    const auto *il_774 = buffer.data(il + 774);
    const auto *il_775 = buffer.data(il + 775);
    const auto *il_776 = buffer.data(il + 776);
    const auto *il_777 = buffer.data(il + 777);
    const auto *il_778 = buffer.data(il + 778);
    const auto *il_779 = buffer.data(il + 779);
    const auto *il_780 = buffer.data(il + 780);
    const auto *il_781 = buffer.data(il + 781);
    const auto *il_782 = buffer.data(il + 782);
    const auto *il_783 = buffer.data(il + 783);
    const auto *il_784 = buffer.data(il + 784);
    const auto *il_785 = buffer.data(il + 785);
    const auto *il_786 = buffer.data(il + 786);
    const auto *il_787 = buffer.data(il + 787);
    const auto *il_788 = buffer.data(il + 788);
    const auto *il_789 = buffer.data(il + 789);
    const auto *il_790 = buffer.data(il + 790);
    const auto *il_791 = buffer.data(il + 791);
    const auto *il_792 = buffer.data(il + 792);
    const auto *il_793 = buffer.data(il + 793);
    const auto *il_794 = buffer.data(il + 794);
    const auto *il_795 = buffer.data(il + 795);
    const auto *il_796 = buffer.data(il + 796);
    const auto *il_797 = buffer.data(il + 797);
    const auto *il_798 = buffer.data(il + 798);
    const auto *il_799 = buffer.data(il + 799);
    const auto *il_800 = buffer.data(il + 800);
    const auto *il_801 = buffer.data(il + 801);
    const auto *il_802 = buffer.data(il + 802);
    const auto *il_803 = buffer.data(il + 803);
    const auto *il_804 = buffer.data(il + 804);
    const auto *il_805 = buffer.data(il + 805);
    const auto *il_806 = buffer.data(il + 806);
    const auto *il_807 = buffer.data(il + 807);
    const auto *il_808 = buffer.data(il + 808);
    const auto *il_809 = buffer.data(il + 809);
    const auto *il_810 = buffer.data(il + 810);
    const auto *il_811 = buffer.data(il + 811);
    const auto *il_812 = buffer.data(il + 812);
    const auto *il_813 = buffer.data(il + 813);
    const auto *il_814 = buffer.data(il + 814);
    const auto *il_815 = buffer.data(il + 815);
    const auto *il_816 = buffer.data(il + 816);
    const auto *il_817 = buffer.data(il + 817);
    const auto *il_818 = buffer.data(il + 818);
    const auto *il_819 = buffer.data(il + 819);
    const auto *il_820 = buffer.data(il + 820);
    const auto *il_821 = buffer.data(il + 821);
    const auto *il_822 = buffer.data(il + 822);
    const auto *il_823 = buffer.data(il + 823);
    const auto *il_824 = buffer.data(il + 824);
    const auto *il_825 = buffer.data(il + 825);
    const auto *il_826 = buffer.data(il + 826);
    const auto *il_827 = buffer.data(il + 827);
    const auto *il_828 = buffer.data(il + 828);
    const auto *il_829 = buffer.data(il + 829);
    const auto *il_830 = buffer.data(il + 830);
    const auto *il_831 = buffer.data(il + 831);
    const auto *il_832 = buffer.data(il + 832);
    const auto *il_833 = buffer.data(il + 833);
    const auto *il_834 = buffer.data(il + 834);
    const auto *il_835 = buffer.data(il + 835);
    const auto *il_836 = buffer.data(il + 836);
    const auto *il_837 = buffer.data(il + 837);
    const auto *il_838 = buffer.data(il + 838);
    const auto *il_839 = buffer.data(il + 839);
    const auto *il_840 = buffer.data(il + 840);
    const auto *il_841 = buffer.data(il + 841);
    const auto *il_842 = buffer.data(il + 842);
    const auto *il_843 = buffer.data(il + 843);
    const auto *il_844 = buffer.data(il + 844);
    const auto *il_845 = buffer.data(il + 845);
    const auto *il_846 = buffer.data(il + 846);
    const auto *il_847 = buffer.data(il + 847);
    const auto *il_848 = buffer.data(il + 848);
    const auto *il_849 = buffer.data(il + 849);
    const auto *il_850 = buffer.data(il + 850);
    const auto *il_851 = buffer.data(il + 851);
    const auto *il_852 = buffer.data(il + 852);
    const auto *il_853 = buffer.data(il + 853);
    const auto *il_854 = buffer.data(il + 854);
    const auto *il_855 = buffer.data(il + 855);
    const auto *il_856 = buffer.data(il + 856);
    const auto *il_857 = buffer.data(il + 857);
    const auto *il_858 = buffer.data(il + 858);
    const auto *il_859 = buffer.data(il + 859);
    const auto *il_860 = buffer.data(il + 860);
    const auto *il_861 = buffer.data(il + 861);
    const auto *il_862 = buffer.data(il + 862);
    const auto *il_863 = buffer.data(il + 863);
    const auto *il_864 = buffer.data(il + 864);
    const auto *il_865 = buffer.data(il + 865);
    const auto *il_866 = buffer.data(il + 866);
    const auto *il_867 = buffer.data(il + 867);
    const auto *il_868 = buffer.data(il + 868);
    const auto *il_869 = buffer.data(il + 869);
    const auto *il_870 = buffer.data(il + 870);
    const auto *il_871 = buffer.data(il + 871);
    const auto *il_872 = buffer.data(il + 872);
    const auto *il_873 = buffer.data(il + 873);
    const auto *il_874 = buffer.data(il + 874);
    const auto *il_875 = buffer.data(il + 875);
    const auto *il_876 = buffer.data(il + 876);
    const auto *il_877 = buffer.data(il + 877);
    const auto *il_878 = buffer.data(il + 878);
    const auto *il_879 = buffer.data(il + 879);
    const auto *il_880 = buffer.data(il + 880);
    const auto *il_881 = buffer.data(il + 881);
    const auto *il_882 = buffer.data(il + 882);
    const auto *il_883 = buffer.data(il + 883);
    const auto *il_884 = buffer.data(il + 884);
    const auto *il_885 = buffer.data(il + 885);
    const auto *il_886 = buffer.data(il + 886);
    const auto *il_887 = buffer.data(il + 887);
    const auto *il_888 = buffer.data(il + 888);
    const auto *il_889 = buffer.data(il + 889);
    const auto *il_890 = buffer.data(il + 890);
    const auto *il_891 = buffer.data(il + 891);
    const auto *il_892 = buffer.data(il + 892);
    const auto *il_893 = buffer.data(il + 893);
    const auto *il_894 = buffer.data(il + 894);
    const auto *il_895 = buffer.data(il + 895);
    const auto *il_896 = buffer.data(il + 896);
    const auto *il_897 = buffer.data(il + 897);

#pragma omp simd aligned(t_507, t_508, t_509, t_510, t_511, gl_327, gl_328, gl_329, gl_330, \
                         gl_331, il_732, il_733, il_734, il_735, \
                         il_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = -3.0 * gl_327[k]
                   + f_0 * il_732[k];

        t_508[k] = -3.0 * gl_328[k]
                   + f_0 * il_733[k];

        t_509[k] = -3.0 * gl_329[k]
                   + f_0 * il_734[k];

        t_510[k] = -3.0 * gl_330[k]
                   + f_0 * il_735[k];

        t_511[k] = -3.0 * gl_331[k]
                   + f_0 * il_736[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, gl_332, gl_333, gl_334, gl_335, \
                         gl_336, il_737, il_738, il_739, il_740, \
                         il_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = -3.0 * gl_332[k]
                   + f_0 * il_737[k];

        t_513[k] = -3.0 * gl_333[k]
                   + f_0 * il_738[k];

        t_514[k] = -3.0 * gl_334[k]
                   + f_0 * il_739[k];

        t_515[k] = -3.0 * gl_335[k]
                   + f_0 * il_740[k];

        t_516[k] = -3.0 * gl_336[k]
                   + f_0 * il_741[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, gl_337, gl_338, gl_339, gl_340, \
                         gl_341, il_742, il_743, il_744, il_745, \
                         il_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = -3.0 * gl_337[k]
                   + f_0 * il_742[k];

        t_518[k] = -3.0 * gl_338[k]
                   + f_0 * il_743[k];

        t_519[k] = -3.0 * gl_339[k]
                   + f_0 * il_744[k];

        t_520[k] = -3.0 * gl_340[k]
                   + f_0 * il_745[k];

        t_521[k] = -3.0 * gl_341[k]
                   + f_0 * il_746[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, gl_342, gl_343, gl_344, gl_345, \
                         gl_346, il_747, il_748, il_749, il_750, \
                         il_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = -3.0 * gl_342[k]
                   + f_0 * il_747[k];

        t_523[k] = -3.0 * gl_343[k]
                   + f_0 * il_748[k];

        t_524[k] = -3.0 * gl_344[k]
                   + f_0 * il_749[k];

        t_525[k] = -3.0 * gl_345[k]
                   + f_0 * il_750[k];

        t_526[k] = -3.0 * gl_346[k]
                   + f_0 * il_751[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, gl_347, gl_348, gl_349, gl_350, \
                         gl_351, il_752, il_753, il_754, il_755, \
                         il_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = -3.0 * gl_347[k]
                   + f_0 * il_752[k];

        t_528[k] = -3.0 * gl_348[k]
                   + f_0 * il_753[k];

        t_529[k] = -3.0 * gl_349[k]
                   + f_0 * il_754[k];

        t_530[k] = -3.0 * gl_350[k]
                   + f_0 * il_755[k];

        t_531[k] = -3.0 * gl_351[k]
                   + f_0 * il_756[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, gl_352, gl_353, gl_354, gl_355, \
                         gl_356, il_757, il_758, il_759, il_760, \
                         il_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = -3.0 * gl_352[k]
                   + f_0 * il_757[k];

        t_533[k] = -3.0 * gl_353[k]
                   + f_0 * il_758[k];

        t_534[k] = -3.0 * gl_354[k]
                   + f_0 * il_759[k];

        t_535[k] = -3.0 * gl_355[k]
                   + f_0 * il_760[k];

        t_536[k] = -3.0 * gl_356[k]
                   + f_0 * il_761[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, t_541, gl_357, gl_358, gl_359, gl_360, \
                         gl_361, il_762, il_763, il_764, il_765, \
                         il_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = -3.0 * gl_357[k]
                   + f_0 * il_762[k];

        t_538[k] = -3.0 * gl_358[k]
                   + f_0 * il_763[k];

        t_539[k] = -3.0 * gl_359[k]
                   + f_0 * il_764[k];

        t_540[k] = -2.0 * gl_360[k]
                   + f_0 * il_765[k];

        t_541[k] = -2.0 * gl_361[k]
                   + f_0 * il_766[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, gl_362, gl_363, gl_364, gl_365, \
                         gl_366, il_767, il_768, il_769, il_770, \
                         il_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = -2.0 * gl_362[k]
                   + f_0 * il_767[k];

        t_543[k] = -2.0 * gl_363[k]
                   + f_0 * il_768[k];

        t_544[k] = -2.0 * gl_364[k]
                   + f_0 * il_769[k];

        t_545[k] = -2.0 * gl_365[k]
                   + f_0 * il_770[k];

        t_546[k] = -2.0 * gl_366[k]
                   + f_0 * il_771[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, gl_367, gl_368, gl_369, gl_370, \
                         gl_371, il_772, il_773, il_774, il_775, \
                         il_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = -2.0 * gl_367[k]
                   + f_0 * il_772[k];

        t_548[k] = -2.0 * gl_368[k]
                   + f_0 * il_773[k];

        t_549[k] = -2.0 * gl_369[k]
                   + f_0 * il_774[k];

        t_550[k] = -2.0 * gl_370[k]
                   + f_0 * il_775[k];

        t_551[k] = -2.0 * gl_371[k]
                   + f_0 * il_776[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, t_556, gl_372, gl_373, gl_374, gl_375, \
                         gl_376, il_777, il_778, il_779, il_780, \
                         il_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = -2.0 * gl_372[k]
                   + f_0 * il_777[k];

        t_553[k] = -2.0 * gl_373[k]
                   + f_0 * il_778[k];

        t_554[k] = -2.0 * gl_374[k]
                   + f_0 * il_779[k];

        t_555[k] = -2.0 * gl_375[k]
                   + f_0 * il_780[k];

        t_556[k] = -2.0 * gl_376[k]
                   + f_0 * il_781[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, gl_377, gl_378, gl_379, gl_380, \
                         gl_381, il_782, il_783, il_784, il_785, \
                         il_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = -2.0 * gl_377[k]
                   + f_0 * il_782[k];

        t_558[k] = -2.0 * gl_378[k]
                   + f_0 * il_783[k];

        t_559[k] = -2.0 * gl_379[k]
                   + f_0 * il_784[k];

        t_560[k] = -2.0 * gl_380[k]
                   + f_0 * il_785[k];

        t_561[k] = -2.0 * gl_381[k]
                   + f_0 * il_786[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, gl_382, gl_383, gl_384, gl_385, \
                         gl_386, il_787, il_788, il_789, il_790, \
                         il_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = -2.0 * gl_382[k]
                   + f_0 * il_787[k];

        t_563[k] = -2.0 * gl_383[k]
                   + f_0 * il_788[k];

        t_564[k] = -2.0 * gl_384[k]
                   + f_0 * il_789[k];

        t_565[k] = -2.0 * gl_385[k]
                   + f_0 * il_790[k];

        t_566[k] = -2.0 * gl_386[k]
                   + f_0 * il_791[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, gl_387, gl_388, gl_389, gl_390, \
                         gl_391, il_792, il_793, il_794, il_795, \
                         il_796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -2.0 * gl_387[k]
                   + f_0 * il_792[k];

        t_568[k] = -2.0 * gl_388[k]
                   + f_0 * il_793[k];

        t_569[k] = -2.0 * gl_389[k]
                   + f_0 * il_794[k];

        t_570[k] = -2.0 * gl_390[k]
                   + f_0 * il_795[k];

        t_571[k] = -2.0 * gl_391[k]
                   + f_0 * il_796[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, gl_392, gl_393, gl_394, gl_395, \
                         gl_396, il_797, il_798, il_799, il_800, \
                         il_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = -2.0 * gl_392[k]
                   + f_0 * il_797[k];

        t_573[k] = -2.0 * gl_393[k]
                   + f_0 * il_798[k];

        t_574[k] = -2.0 * gl_394[k]
                   + f_0 * il_799[k];

        t_575[k] = -2.0 * gl_395[k]
                   + f_0 * il_800[k];

        t_576[k] = -2.0 * gl_396[k]
                   + f_0 * il_801[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, gl_397, gl_398, gl_399, gl_400, \
                         gl_401, il_802, il_803, il_804, il_805, \
                         il_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = -2.0 * gl_397[k]
                   + f_0 * il_802[k];

        t_578[k] = -2.0 * gl_398[k]
                   + f_0 * il_803[k];

        t_579[k] = -2.0 * gl_399[k]
                   + f_0 * il_804[k];

        t_580[k] = -2.0 * gl_400[k]
                   + f_0 * il_805[k];

        t_581[k] = -2.0 * gl_401[k]
                   + f_0 * il_806[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, gl_402, gl_403, gl_404, gl_405, \
                         gl_406, il_807, il_808, il_809, il_810, \
                         il_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = -2.0 * gl_402[k]
                   + f_0 * il_807[k];

        t_583[k] = -2.0 * gl_403[k]
                   + f_0 * il_808[k];

        t_584[k] = -2.0 * gl_404[k]
                   + f_0 * il_809[k];

        t_585[k] = -gl_405[k]
                   + f_0 * il_810[k];

        t_586[k] = -gl_406[k]
                   + f_0 * il_811[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, t_591, gl_407, gl_408, gl_409, gl_410, \
                         gl_411, il_812, il_813, il_814, il_815, \
                         il_816 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = -gl_407[k]
                   + f_0 * il_812[k];

        t_588[k] = -gl_408[k]
                   + f_0 * il_813[k];

        t_589[k] = -gl_409[k]
                   + f_0 * il_814[k];

        t_590[k] = -gl_410[k]
                   + f_0 * il_815[k];

        t_591[k] = -gl_411[k]
                   + f_0 * il_816[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, t_596, gl_412, gl_413, gl_414, gl_415, \
                         gl_416, il_817, il_818, il_819, il_820, \
                         il_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = -gl_412[k]
                   + f_0 * il_817[k];

        t_593[k] = -gl_413[k]
                   + f_0 * il_818[k];

        t_594[k] = -gl_414[k]
                   + f_0 * il_819[k];

        t_595[k] = -gl_415[k]
                   + f_0 * il_820[k];

        t_596[k] = -gl_416[k]
                   + f_0 * il_821[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, t_600, t_601, gl_417, gl_418, gl_419, gl_420, \
                         gl_421, il_822, il_823, il_824, il_825, \
                         il_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = -gl_417[k]
                   + f_0 * il_822[k];

        t_598[k] = -gl_418[k]
                   + f_0 * il_823[k];

        t_599[k] = -gl_419[k]
                   + f_0 * il_824[k];

        t_600[k] = -gl_420[k]
                   + f_0 * il_825[k];

        t_601[k] = -gl_421[k]
                   + f_0 * il_826[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, gl_422, gl_423, gl_424, gl_425, \
                         gl_426, il_827, il_828, il_829, il_830, \
                         il_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = -gl_422[k]
                   + f_0 * il_827[k];

        t_603[k] = -gl_423[k]
                   + f_0 * il_828[k];

        t_604[k] = -gl_424[k]
                   + f_0 * il_829[k];

        t_605[k] = -gl_425[k]
                   + f_0 * il_830[k];

        t_606[k] = -gl_426[k]
                   + f_0 * il_831[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, t_611, gl_427, gl_428, gl_429, gl_430, \
                         gl_431, il_832, il_833, il_834, il_835, \
                         il_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = -gl_427[k]
                   + f_0 * il_832[k];

        t_608[k] = -gl_428[k]
                   + f_0 * il_833[k];

        t_609[k] = -gl_429[k]
                   + f_0 * il_834[k];

        t_610[k] = -gl_430[k]
                   + f_0 * il_835[k];

        t_611[k] = -gl_431[k]
                   + f_0 * il_836[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, gl_432, gl_433, gl_434, gl_435, \
                         gl_436, il_837, il_838, il_839, il_840, \
                         il_841 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = -gl_432[k]
                   + f_0 * il_837[k];

        t_613[k] = -gl_433[k]
                   + f_0 * il_838[k];

        t_614[k] = -gl_434[k]
                   + f_0 * il_839[k];

        t_615[k] = -gl_435[k]
                   + f_0 * il_840[k];

        t_616[k] = -gl_436[k]
                   + f_0 * il_841[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, gl_437, gl_438, gl_439, gl_440, \
                         gl_441, il_842, il_843, il_844, il_845, \
                         il_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = -gl_437[k]
                   + f_0 * il_842[k];

        t_618[k] = -gl_438[k]
                   + f_0 * il_843[k];

        t_619[k] = -gl_439[k]
                   + f_0 * il_844[k];

        t_620[k] = -gl_440[k]
                   + f_0 * il_845[k];

        t_621[k] = -gl_441[k]
                   + f_0 * il_846[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, gl_442, gl_443, gl_444, gl_445, \
                         gl_446, il_847, il_848, il_849, il_850, \
                         il_851 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = -gl_442[k]
                   + f_0 * il_847[k];

        t_623[k] = -gl_443[k]
                   + f_0 * il_848[k];

        t_624[k] = -gl_444[k]
                   + f_0 * il_849[k];

        t_625[k] = -gl_445[k]
                   + f_0 * il_850[k];

        t_626[k] = -gl_446[k]
                   + f_0 * il_851[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, t_631, t_632, gl_447, gl_448, gl_449, \
                         il_852, il_853, il_854, il_855, il_856, \
                         il_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = -gl_447[k]
                   + f_0 * il_852[k];

        t_628[k] = -gl_448[k]
                   + f_0 * il_853[k];

        t_629[k] = -gl_449[k]
                   + f_0 * il_854[k];

        t_630[k] = f_0 * il_855[k];

        t_631[k] = f_0 * il_856[k];

        t_632[k] = f_0 * il_857[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, t_636, t_637, t_638, t_639, t_640, il_858, \
                         il_859, il_860, il_861, il_862, il_863, il_864, \
                         il_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_0 * il_858[k];

        t_634[k] = f_0 * il_859[k];

        t_635[k] = f_0 * il_860[k];

        t_636[k] = f_0 * il_861[k];

        t_637[k] = f_0 * il_862[k];

        t_638[k] = f_0 * il_863[k];

        t_639[k] = f_0 * il_864[k];

        t_640[k] = f_0 * il_865[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, t_645, t_646, t_647, t_648, il_866, \
                         il_867, il_868, il_869, il_870, il_871, il_872, \
                         il_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_0 * il_866[k];

        t_642[k] = f_0 * il_867[k];

        t_643[k] = f_0 * il_868[k];

        t_644[k] = f_0 * il_869[k];

        t_645[k] = f_0 * il_870[k];

        t_646[k] = f_0 * il_871[k];

        t_647[k] = f_0 * il_872[k];

        t_648[k] = f_0 * il_873[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, t_654, t_655, t_656, il_874, \
                         il_875, il_876, il_877, il_878, il_879, il_880, \
                         il_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_0 * il_874[k];

        t_650[k] = f_0 * il_875[k];

        t_651[k] = f_0 * il_876[k];

        t_652[k] = f_0 * il_877[k];

        t_653[k] = f_0 * il_878[k];

        t_654[k] = f_0 * il_879[k];

        t_655[k] = f_0 * il_880[k];

        t_656[k] = f_0 * il_881[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, t_662, t_663, t_664, il_882, \
                         il_883, il_884, il_885, il_886, il_887, il_888, \
                         il_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_0 * il_882[k];

        t_658[k] = f_0 * il_883[k];

        t_659[k] = f_0 * il_884[k];

        t_660[k] = f_0 * il_885[k];

        t_661[k] = f_0 * il_886[k];

        t_662[k] = f_0 * il_887[k];

        t_663[k] = f_0 * il_888[k];

        t_664[k] = f_0 * il_889[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, t_670, t_671, t_672, il_890, \
                         il_891, il_892, il_893, il_894, il_895, il_896, \
                         il_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_0 * il_890[k];

        t_666[k] = f_0 * il_891[k];

        t_667[k] = f_0 * il_892[k];

        t_668[k] = f_0 * il_893[k];

        t_669[k] = f_0 * il_894[k];

        t_670[k] = f_0 * il_895[k];

        t_671[k] = f_0 * il_896[k];

        t_672[k] = f_0 * il_897[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gl_450 = buffer.data(gl + 450);
    const auto *gl_451 = buffer.data(gl + 451);
    const auto *gl_452 = buffer.data(gl + 452);
    const auto *gl_453 = buffer.data(gl + 453);
    const auto *gl_454 = buffer.data(gl + 454);
    const auto *gl_455 = buffer.data(gl + 455);
    const auto *gl_456 = buffer.data(gl + 456);
    const auto *gl_457 = buffer.data(gl + 457);
    const auto *gl_458 = buffer.data(gl + 458);
    const auto *gl_459 = buffer.data(gl + 459);
    const auto *gl_460 = buffer.data(gl + 460);
    const auto *gl_461 = buffer.data(gl + 461);
    const auto *gl_462 = buffer.data(gl + 462);
    const auto *gl_463 = buffer.data(gl + 463);
    const auto *gl_464 = buffer.data(gl + 464);
    const auto *gl_465 = buffer.data(gl + 465);
    const auto *gl_466 = buffer.data(gl + 466);
    const auto *gl_467 = buffer.data(gl + 467);
    const auto *gl_468 = buffer.data(gl + 468);
    const auto *gl_469 = buffer.data(gl + 469);
    const auto *gl_470 = buffer.data(gl + 470);
    const auto *gl_471 = buffer.data(gl + 471);
    const auto *gl_472 = buffer.data(gl + 472);
    const auto *gl_473 = buffer.data(gl + 473);
    const auto *gl_474 = buffer.data(gl + 474);
    const auto *gl_475 = buffer.data(gl + 475);
    const auto *gl_476 = buffer.data(gl + 476);
    const auto *gl_477 = buffer.data(gl + 477);
    const auto *gl_478 = buffer.data(gl + 478);
    const auto *gl_479 = buffer.data(gl + 479);
    const auto *gl_480 = buffer.data(gl + 480);
    const auto *gl_481 = buffer.data(gl + 481);
    const auto *gl_482 = buffer.data(gl + 482);
    const auto *gl_483 = buffer.data(gl + 483);
    const auto *gl_484 = buffer.data(gl + 484);
    const auto *gl_485 = buffer.data(gl + 485);
    const auto *gl_486 = buffer.data(gl + 486);
    const auto *gl_487 = buffer.data(gl + 487);
    const auto *gl_488 = buffer.data(gl + 488);
    const auto *gl_489 = buffer.data(gl + 489);
    const auto *gl_490 = buffer.data(gl + 490);
    const auto *gl_491 = buffer.data(gl + 491);
    const auto *gl_492 = buffer.data(gl + 492);
    const auto *gl_493 = buffer.data(gl + 493);
    const auto *gl_494 = buffer.data(gl + 494);
    const auto *gl_495 = buffer.data(gl + 495);
    const auto *gl_496 = buffer.data(gl + 496);
    const auto *gl_497 = buffer.data(gl + 497);
    const auto *gl_498 = buffer.data(gl + 498);
    const auto *gl_499 = buffer.data(gl + 499);
    const auto *gl_500 = buffer.data(gl + 500);
    const auto *gl_501 = buffer.data(gl + 501);
    const auto *gl_502 = buffer.data(gl + 502);
    const auto *gl_503 = buffer.data(gl + 503);
    const auto *gl_504 = buffer.data(gl + 504);
    const auto *gl_505 = buffer.data(gl + 505);
    const auto *gl_506 = buffer.data(gl + 506);
    const auto *gl_507 = buffer.data(gl + 507);
    const auto *gl_508 = buffer.data(gl + 508);
    const auto *gl_509 = buffer.data(gl + 509);
    const auto *gl_510 = buffer.data(gl + 510);
    const auto *gl_511 = buffer.data(gl + 511);
    const auto *gl_512 = buffer.data(gl + 512);
    const auto *gl_513 = buffer.data(gl + 513);
    const auto *gl_514 = buffer.data(gl + 514);
    const auto *gl_515 = buffer.data(gl + 515);
    const auto *gl_516 = buffer.data(gl + 516);
    const auto *gl_517 = buffer.data(gl + 517);
    const auto *gl_518 = buffer.data(gl + 518);
    const auto *gl_519 = buffer.data(gl + 519);
    const auto *gl_520 = buffer.data(gl + 520);
    const auto *gl_521 = buffer.data(gl + 521);
    const auto *gl_522 = buffer.data(gl + 522);
    const auto *gl_523 = buffer.data(gl + 523);
    const auto *gl_524 = buffer.data(gl + 524);
    const auto *gl_525 = buffer.data(gl + 525);
    const auto *gl_526 = buffer.data(gl + 526);
    const auto *gl_527 = buffer.data(gl + 527);
    const auto *gl_528 = buffer.data(gl + 528);
    const auto *gl_529 = buffer.data(gl + 529);
    const auto *gl_530 = buffer.data(gl + 530);
    const auto *gl_531 = buffer.data(gl + 531);
    const auto *gl_532 = buffer.data(gl + 532);
    const auto *gl_533 = buffer.data(gl + 533);
    const auto *gl_534 = buffer.data(gl + 534);
    const auto *gl_535 = buffer.data(gl + 535);
    const auto *gl_536 = buffer.data(gl + 536);
    const auto *gl_537 = buffer.data(gl + 537);
    const auto *gl_538 = buffer.data(gl + 538);
    const auto *gl_539 = buffer.data(gl + 539);
    const auto *gl_540 = buffer.data(gl + 540);
    const auto *gl_541 = buffer.data(gl + 541);
    const auto *gl_542 = buffer.data(gl + 542);
    const auto *gl_543 = buffer.data(gl + 543);
    const auto *gl_544 = buffer.data(gl + 544);
    const auto *gl_545 = buffer.data(gl + 545);
    const auto *gl_546 = buffer.data(gl + 546);
    const auto *gl_547 = buffer.data(gl + 547);
    const auto *gl_548 = buffer.data(gl + 548);
    const auto *gl_549 = buffer.data(gl + 549);
    const auto *gl_550 = buffer.data(gl + 550);
    const auto *gl_551 = buffer.data(gl + 551);
    const auto *gl_552 = buffer.data(gl + 552);
    const auto *gl_553 = buffer.data(gl + 553);
    const auto *gl_554 = buffer.data(gl + 554);
    const auto *gl_555 = buffer.data(gl + 555);
    const auto *gl_556 = buffer.data(gl + 556);
    const auto *gl_557 = buffer.data(gl + 557);
    const auto *gl_558 = buffer.data(gl + 558);
    const auto *gl_559 = buffer.data(gl + 559);
    const auto *gl_560 = buffer.data(gl + 560);
    const auto *gl_561 = buffer.data(gl + 561);
    const auto *gl_562 = buffer.data(gl + 562);
    const auto *gl_563 = buffer.data(gl + 563);
    const auto *gl_564 = buffer.data(gl + 564);
    const auto *gl_565 = buffer.data(gl + 565);
    const auto *gl_566 = buffer.data(gl + 566);
    const auto *gl_567 = buffer.data(gl + 567);
    const auto *gl_568 = buffer.data(gl + 568);
    const auto *gl_569 = buffer.data(gl + 569);
    const auto *gl_570 = buffer.data(gl + 570);
    const auto *gl_571 = buffer.data(gl + 571);
    const auto *gl_572 = buffer.data(gl + 572);
    const auto *gl_573 = buffer.data(gl + 573);
    const auto *gl_574 = buffer.data(gl + 574);
    const auto *gl_575 = buffer.data(gl + 575);
    const auto *gl_576 = buffer.data(gl + 576);
    const auto *gl_577 = buffer.data(gl + 577);
    const auto *gl_578 = buffer.data(gl + 578);
    const auto *gl_579 = buffer.data(gl + 579);
    const auto *gl_580 = buffer.data(gl + 580);
    const auto *gl_581 = buffer.data(gl + 581);
    const auto *gl_582 = buffer.data(gl + 582);
    const auto *gl_583 = buffer.data(gl + 583);
    const auto *gl_584 = buffer.data(gl + 584);
    const auto *gl_585 = buffer.data(gl + 585);
    const auto *gl_586 = buffer.data(gl + 586);
    const auto *gl_587 = buffer.data(gl + 587);
    const auto *gl_588 = buffer.data(gl + 588);
    const auto *gl_589 = buffer.data(gl + 589);
    const auto *gl_590 = buffer.data(gl + 590);
    const auto *gl_591 = buffer.data(gl + 591);
    const auto *gl_592 = buffer.data(gl + 592);
    const auto *gl_593 = buffer.data(gl + 593);
    const auto *gl_594 = buffer.data(gl + 594);
    const auto *gl_595 = buffer.data(gl + 595);
    const auto *gl_596 = buffer.data(gl + 596);
    const auto *gl_597 = buffer.data(gl + 597);
    const auto *gl_598 = buffer.data(gl + 598);

    const auto *il_898 = buffer.data(il + 898);
    const auto *il_899 = buffer.data(il + 899);
    const auto *il_945 = buffer.data(il + 945);
    const auto *il_946 = buffer.data(il + 946);
    const auto *il_947 = buffer.data(il + 947);
    const auto *il_948 = buffer.data(il + 948);
    const auto *il_949 = buffer.data(il + 949);
    const auto *il_950 = buffer.data(il + 950);
    const auto *il_951 = buffer.data(il + 951);
    const auto *il_952 = buffer.data(il + 952);
    const auto *il_953 = buffer.data(il + 953);
    const auto *il_954 = buffer.data(il + 954);
    const auto *il_955 = buffer.data(il + 955);
    const auto *il_956 = buffer.data(il + 956);
    const auto *il_957 = buffer.data(il + 957);
    const auto *il_958 = buffer.data(il + 958);
    const auto *il_959 = buffer.data(il + 959);
    const auto *il_960 = buffer.data(il + 960);
    const auto *il_961 = buffer.data(il + 961);
    const auto *il_962 = buffer.data(il + 962);
    const auto *il_963 = buffer.data(il + 963);
    const auto *il_964 = buffer.data(il + 964);
    const auto *il_965 = buffer.data(il + 965);
    const auto *il_966 = buffer.data(il + 966);
    const auto *il_967 = buffer.data(il + 967);
    const auto *il_968 = buffer.data(il + 968);
    const auto *il_969 = buffer.data(il + 969);
    const auto *il_970 = buffer.data(il + 970);
    const auto *il_971 = buffer.data(il + 971);
    const auto *il_972 = buffer.data(il + 972);
    const auto *il_973 = buffer.data(il + 973);
    const auto *il_974 = buffer.data(il + 974);
    const auto *il_975 = buffer.data(il + 975);
    const auto *il_976 = buffer.data(il + 976);
    const auto *il_977 = buffer.data(il + 977);
    const auto *il_978 = buffer.data(il + 978);
    const auto *il_979 = buffer.data(il + 979);
    const auto *il_980 = buffer.data(il + 980);
    const auto *il_981 = buffer.data(il + 981);
    const auto *il_982 = buffer.data(il + 982);
    const auto *il_983 = buffer.data(il + 983);
    const auto *il_984 = buffer.data(il + 984);
    const auto *il_985 = buffer.data(il + 985);
    const auto *il_986 = buffer.data(il + 986);
    const auto *il_987 = buffer.data(il + 987);
    const auto *il_988 = buffer.data(il + 988);
    const auto *il_989 = buffer.data(il + 989);
    const auto *il_990 = buffer.data(il + 990);
    const auto *il_991 = buffer.data(il + 991);
    const auto *il_992 = buffer.data(il + 992);
    const auto *il_993 = buffer.data(il + 993);
    const auto *il_994 = buffer.data(il + 994);
    const auto *il_995 = buffer.data(il + 995);
    const auto *il_996 = buffer.data(il + 996);
    const auto *il_997 = buffer.data(il + 997);
    const auto *il_998 = buffer.data(il + 998);
    const auto *il_999 = buffer.data(il + 999);
    const auto *il_1000 = buffer.data(il + 1000);
    const auto *il_1001 = buffer.data(il + 1001);
    const auto *il_1002 = buffer.data(il + 1002);
    const auto *il_1003 = buffer.data(il + 1003);
    const auto *il_1004 = buffer.data(il + 1004);
    const auto *il_1005 = buffer.data(il + 1005);
    const auto *il_1006 = buffer.data(il + 1006);
    const auto *il_1007 = buffer.data(il + 1007);
    const auto *il_1008 = buffer.data(il + 1008);
    const auto *il_1009 = buffer.data(il + 1009);
    const auto *il_1010 = buffer.data(il + 1010);
    const auto *il_1011 = buffer.data(il + 1011);
    const auto *il_1012 = buffer.data(il + 1012);
    const auto *il_1013 = buffer.data(il + 1013);
    const auto *il_1014 = buffer.data(il + 1014);
    const auto *il_1015 = buffer.data(il + 1015);
    const auto *il_1016 = buffer.data(il + 1016);
    const auto *il_1017 = buffer.data(il + 1017);
    const auto *il_1018 = buffer.data(il + 1018);
    const auto *il_1019 = buffer.data(il + 1019);
    const auto *il_1020 = buffer.data(il + 1020);
    const auto *il_1021 = buffer.data(il + 1021);
    const auto *il_1022 = buffer.data(il + 1022);
    const auto *il_1023 = buffer.data(il + 1023);
    const auto *il_1024 = buffer.data(il + 1024);
    const auto *il_1025 = buffer.data(il + 1025);
    const auto *il_1026 = buffer.data(il + 1026);
    const auto *il_1027 = buffer.data(il + 1027);
    const auto *il_1028 = buffer.data(il + 1028);
    const auto *il_1029 = buffer.data(il + 1029);
    const auto *il_1030 = buffer.data(il + 1030);
    const auto *il_1031 = buffer.data(il + 1031);
    const auto *il_1032 = buffer.data(il + 1032);
    const auto *il_1033 = buffer.data(il + 1033);
    const auto *il_1034 = buffer.data(il + 1034);
    const auto *il_1035 = buffer.data(il + 1035);
    const auto *il_1036 = buffer.data(il + 1036);
    const auto *il_1037 = buffer.data(il + 1037);
    const auto *il_1038 = buffer.data(il + 1038);
    const auto *il_1039 = buffer.data(il + 1039);
    const auto *il_1040 = buffer.data(il + 1040);
    const auto *il_1041 = buffer.data(il + 1041);
    const auto *il_1042 = buffer.data(il + 1042);
    const auto *il_1043 = buffer.data(il + 1043);
    const auto *il_1044 = buffer.data(il + 1044);
    const auto *il_1045 = buffer.data(il + 1045);
    const auto *il_1046 = buffer.data(il + 1046);
    const auto *il_1047 = buffer.data(il + 1047);
    const auto *il_1048 = buffer.data(il + 1048);
    const auto *il_1049 = buffer.data(il + 1049);
    const auto *il_1050 = buffer.data(il + 1050);
    const auto *il_1051 = buffer.data(il + 1051);
    const auto *il_1052 = buffer.data(il + 1052);
    const auto *il_1053 = buffer.data(il + 1053);
    const auto *il_1054 = buffer.data(il + 1054);
    const auto *il_1055 = buffer.data(il + 1055);
    const auto *il_1056 = buffer.data(il + 1056);
    const auto *il_1057 = buffer.data(il + 1057);
    const auto *il_1058 = buffer.data(il + 1058);
    const auto *il_1059 = buffer.data(il + 1059);
    const auto *il_1060 = buffer.data(il + 1060);
    const auto *il_1061 = buffer.data(il + 1061);
    const auto *il_1062 = buffer.data(il + 1062);
    const auto *il_1063 = buffer.data(il + 1063);
    const auto *il_1064 = buffer.data(il + 1064);
    const auto *il_1065 = buffer.data(il + 1065);
    const auto *il_1066 = buffer.data(il + 1066);
    const auto *il_1067 = buffer.data(il + 1067);
    const auto *il_1068 = buffer.data(il + 1068);
    const auto *il_1069 = buffer.data(il + 1069);
    const auto *il_1070 = buffer.data(il + 1070);
    const auto *il_1071 = buffer.data(il + 1071);
    const auto *il_1072 = buffer.data(il + 1072);
    const auto *il_1073 = buffer.data(il + 1073);
    const auto *il_1074 = buffer.data(il + 1074);
    const auto *il_1075 = buffer.data(il + 1075);
    const auto *il_1076 = buffer.data(il + 1076);
    const auto *il_1077 = buffer.data(il + 1077);
    const auto *il_1078 = buffer.data(il + 1078);
    const auto *il_1079 = buffer.data(il + 1079);
    const auto *il_1080 = buffer.data(il + 1080);
    const auto *il_1081 = buffer.data(il + 1081);
    const auto *il_1082 = buffer.data(il + 1082);
    const auto *il_1083 = buffer.data(il + 1083);
    const auto *il_1084 = buffer.data(il + 1084);
    const auto *il_1085 = buffer.data(il + 1085);
    const auto *il_1086 = buffer.data(il + 1086);
    const auto *il_1087 = buffer.data(il + 1087);
    const auto *il_1088 = buffer.data(il + 1088);
    const auto *il_1089 = buffer.data(il + 1089);
    const auto *il_1090 = buffer.data(il + 1090);
    const auto *il_1091 = buffer.data(il + 1091);
    const auto *il_1092 = buffer.data(il + 1092);
    const auto *il_1093 = buffer.data(il + 1093);

#pragma omp simd aligned(t_673, t_674, t_675, t_676, t_677, t_678, gl_450, gl_451, gl_452, \
                         gl_453, il_898, il_899, il_945, il_946, il_947, \
                         il_948 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_0 * il_898[k];

        t_674[k] = f_0 * il_899[k];

        t_675[k] = -5.0 * gl_450[k]
                   + f_0 * il_945[k];

        t_676[k] = -5.0 * gl_451[k]
                   + f_0 * il_946[k];

        t_677[k] = -5.0 * gl_452[k]
                   + f_0 * il_947[k];

        t_678[k] = -5.0 * gl_453[k]
                   + f_0 * il_948[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, t_682, t_683, gl_454, gl_455, gl_456, gl_457, \
                         gl_458, il_949, il_950, il_951, il_952, \
                         il_953 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = -5.0 * gl_454[k]
                   + f_0 * il_949[k];

        t_680[k] = -5.0 * gl_455[k]
                   + f_0 * il_950[k];

        t_681[k] = -5.0 * gl_456[k]
                   + f_0 * il_951[k];

        t_682[k] = -5.0 * gl_457[k]
                   + f_0 * il_952[k];

        t_683[k] = -5.0 * gl_458[k]
                   + f_0 * il_953[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, t_688, gl_459, gl_460, gl_461, gl_462, \
                         gl_463, il_954, il_955, il_956, il_957, \
                         il_958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = -5.0 * gl_459[k]
                   + f_0 * il_954[k];

        t_685[k] = -5.0 * gl_460[k]
                   + f_0 * il_955[k];

        t_686[k] = -5.0 * gl_461[k]
                   + f_0 * il_956[k];

        t_687[k] = -5.0 * gl_462[k]
                   + f_0 * il_957[k];

        t_688[k] = -5.0 * gl_463[k]
                   + f_0 * il_958[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, t_693, gl_464, gl_465, gl_466, gl_467, \
                         gl_468, il_959, il_960, il_961, il_962, \
                         il_963 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = -5.0 * gl_464[k]
                   + f_0 * il_959[k];

        t_690[k] = -5.0 * gl_465[k]
                   + f_0 * il_960[k];

        t_691[k] = -5.0 * gl_466[k]
                   + f_0 * il_961[k];

        t_692[k] = -5.0 * gl_467[k]
                   + f_0 * il_962[k];

        t_693[k] = -5.0 * gl_468[k]
                   + f_0 * il_963[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, t_698, gl_469, gl_470, gl_471, gl_472, \
                         gl_473, il_964, il_965, il_966, il_967, \
                         il_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = -5.0 * gl_469[k]
                   + f_0 * il_964[k];

        t_695[k] = -5.0 * gl_470[k]
                   + f_0 * il_965[k];

        t_696[k] = -5.0 * gl_471[k]
                   + f_0 * il_966[k];

        t_697[k] = -5.0 * gl_472[k]
                   + f_0 * il_967[k];

        t_698[k] = -5.0 * gl_473[k]
                   + f_0 * il_968[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, t_703, gl_474, gl_475, gl_476, gl_477, \
                         gl_478, il_969, il_970, il_971, il_972, \
                         il_973 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = -5.0 * gl_474[k]
                   + f_0 * il_969[k];

        t_700[k] = -5.0 * gl_475[k]
                   + f_0 * il_970[k];

        t_701[k] = -5.0 * gl_476[k]
                   + f_0 * il_971[k];

        t_702[k] = -5.0 * gl_477[k]
                   + f_0 * il_972[k];

        t_703[k] = -5.0 * gl_478[k]
                   + f_0 * il_973[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, gl_479, gl_480, gl_481, gl_482, \
                         gl_483, il_974, il_975, il_976, il_977, \
                         il_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = -5.0 * gl_479[k]
                   + f_0 * il_974[k];

        t_705[k] = -5.0 * gl_480[k]
                   + f_0 * il_975[k];

        t_706[k] = -5.0 * gl_481[k]
                   + f_0 * il_976[k];

        t_707[k] = -5.0 * gl_482[k]
                   + f_0 * il_977[k];

        t_708[k] = -5.0 * gl_483[k]
                   + f_0 * il_978[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, t_712, t_713, gl_484, gl_485, gl_486, gl_487, \
                         gl_488, il_979, il_980, il_981, il_982, \
                         il_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = -5.0 * gl_484[k]
                   + f_0 * il_979[k];

        t_710[k] = -5.0 * gl_485[k]
                   + f_0 * il_980[k];

        t_711[k] = -5.0 * gl_486[k]
                   + f_0 * il_981[k];

        t_712[k] = -5.0 * gl_487[k]
                   + f_0 * il_982[k];

        t_713[k] = -5.0 * gl_488[k]
                   + f_0 * il_983[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, t_718, gl_489, gl_490, gl_491, gl_492, \
                         gl_493, il_984, il_985, il_986, il_987, \
                         il_988 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = -5.0 * gl_489[k]
                   + f_0 * il_984[k];

        t_715[k] = -5.0 * gl_490[k]
                   + f_0 * il_985[k];

        t_716[k] = -5.0 * gl_491[k]
                   + f_0 * il_986[k];

        t_717[k] = -5.0 * gl_492[k]
                   + f_0 * il_987[k];

        t_718[k] = -5.0 * gl_493[k]
                   + f_0 * il_988[k];
    }

#pragma omp simd aligned(t_719, t_720, t_721, t_722, t_723, gl_494, gl_495, gl_496, gl_497, \
                         gl_498, il_989, il_990, il_991, il_992, \
                         il_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_719[k] = -5.0 * gl_494[k]
                   + f_0 * il_989[k];

        t_720[k] = -4.0 * gl_495[k]
                   + f_0 * il_990[k];

        t_721[k] = -4.0 * gl_496[k]
                   + f_0 * il_991[k];

        t_722[k] = -4.0 * gl_497[k]
                   + f_0 * il_992[k];

        t_723[k] = -4.0 * gl_498[k]
                   + f_0 * il_993[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, gl_499, gl_500, gl_501, gl_502, \
                         gl_503, il_994, il_995, il_996, il_997, \
                         il_998 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = -4.0 * gl_499[k]
                   + f_0 * il_994[k];

        t_725[k] = -4.0 * gl_500[k]
                   + f_0 * il_995[k];

        t_726[k] = -4.0 * gl_501[k]
                   + f_0 * il_996[k];

        t_727[k] = -4.0 * gl_502[k]
                   + f_0 * il_997[k];

        t_728[k] = -4.0 * gl_503[k]
                   + f_0 * il_998[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, t_733, gl_504, gl_505, gl_506, gl_507, \
                         gl_508, il_999, il_1000, il_1001, il_1002, \
                         il_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = -4.0 * gl_504[k]
                   + f_0 * il_999[k];

        t_730[k] = -4.0 * gl_505[k]
                   + f_0 * il_1000[k];

        t_731[k] = -4.0 * gl_506[k]
                   + f_0 * il_1001[k];

        t_732[k] = -4.0 * gl_507[k]
                   + f_0 * il_1002[k];

        t_733[k] = -4.0 * gl_508[k]
                   + f_0 * il_1003[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, t_738, gl_509, gl_510, gl_511, gl_512, \
                         gl_513, il_1004, il_1005, il_1006, il_1007, \
                         il_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = -4.0 * gl_509[k]
                   + f_0 * il_1004[k];

        t_735[k] = -4.0 * gl_510[k]
                   + f_0 * il_1005[k];

        t_736[k] = -4.0 * gl_511[k]
                   + f_0 * il_1006[k];

        t_737[k] = -4.0 * gl_512[k]
                   + f_0 * il_1007[k];

        t_738[k] = -4.0 * gl_513[k]
                   + f_0 * il_1008[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, t_743, gl_514, gl_515, gl_516, gl_517, \
                         gl_518, il_1009, il_1010, il_1011, il_1012, \
                         il_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = -4.0 * gl_514[k]
                   + f_0 * il_1009[k];

        t_740[k] = -4.0 * gl_515[k]
                   + f_0 * il_1010[k];

        t_741[k] = -4.0 * gl_516[k]
                   + f_0 * il_1011[k];

        t_742[k] = -4.0 * gl_517[k]
                   + f_0 * il_1012[k];

        t_743[k] = -4.0 * gl_518[k]
                   + f_0 * il_1013[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, t_748, gl_519, gl_520, gl_521, gl_522, \
                         gl_523, il_1014, il_1015, il_1016, il_1017, \
                         il_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = -4.0 * gl_519[k]
                   + f_0 * il_1014[k];

        t_745[k] = -4.0 * gl_520[k]
                   + f_0 * il_1015[k];

        t_746[k] = -4.0 * gl_521[k]
                   + f_0 * il_1016[k];

        t_747[k] = -4.0 * gl_522[k]
                   + f_0 * il_1017[k];

        t_748[k] = -4.0 * gl_523[k]
                   + f_0 * il_1018[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, t_752, t_753, gl_524, gl_525, gl_526, gl_527, \
                         gl_528, il_1019, il_1020, il_1021, il_1022, \
                         il_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = -4.0 * gl_524[k]
                   + f_0 * il_1019[k];

        t_750[k] = -4.0 * gl_525[k]
                   + f_0 * il_1020[k];

        t_751[k] = -4.0 * gl_526[k]
                   + f_0 * il_1021[k];

        t_752[k] = -4.0 * gl_527[k]
                   + f_0 * il_1022[k];

        t_753[k] = -4.0 * gl_528[k]
                   + f_0 * il_1023[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, t_758, gl_529, gl_530, gl_531, gl_532, \
                         gl_533, il_1024, il_1025, il_1026, il_1027, \
                         il_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = -4.0 * gl_529[k]
                   + f_0 * il_1024[k];

        t_755[k] = -4.0 * gl_530[k]
                   + f_0 * il_1025[k];

        t_756[k] = -4.0 * gl_531[k]
                   + f_0 * il_1026[k];

        t_757[k] = -4.0 * gl_532[k]
                   + f_0 * il_1027[k];

        t_758[k] = -4.0 * gl_533[k]
                   + f_0 * il_1028[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, t_762, t_763, gl_534, gl_535, gl_536, gl_537, \
                         gl_538, il_1029, il_1030, il_1031, il_1032, \
                         il_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = -4.0 * gl_534[k]
                   + f_0 * il_1029[k];

        t_760[k] = -4.0 * gl_535[k]
                   + f_0 * il_1030[k];

        t_761[k] = -4.0 * gl_536[k]
                   + f_0 * il_1031[k];

        t_762[k] = -4.0 * gl_537[k]
                   + f_0 * il_1032[k];

        t_763[k] = -4.0 * gl_538[k]
                   + f_0 * il_1033[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, t_768, gl_539, gl_540, gl_541, gl_542, \
                         gl_543, il_1034, il_1035, il_1036, il_1037, \
                         il_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = -4.0 * gl_539[k]
                   + f_0 * il_1034[k];

        t_765[k] = -3.0 * gl_540[k]
                   + f_0 * il_1035[k];

        t_766[k] = -3.0 * gl_541[k]
                   + f_0 * il_1036[k];

        t_767[k] = -3.0 * gl_542[k]
                   + f_0 * il_1037[k];

        t_768[k] = -3.0 * gl_543[k]
                   + f_0 * il_1038[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, t_772, t_773, gl_544, gl_545, gl_546, gl_547, \
                         gl_548, il_1039, il_1040, il_1041, il_1042, \
                         il_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = -3.0 * gl_544[k]
                   + f_0 * il_1039[k];

        t_770[k] = -3.0 * gl_545[k]
                   + f_0 * il_1040[k];

        t_771[k] = -3.0 * gl_546[k]
                   + f_0 * il_1041[k];

        t_772[k] = -3.0 * gl_547[k]
                   + f_0 * il_1042[k];

        t_773[k] = -3.0 * gl_548[k]
                   + f_0 * il_1043[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, t_778, gl_549, gl_550, gl_551, gl_552, \
                         gl_553, il_1044, il_1045, il_1046, il_1047, \
                         il_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = -3.0 * gl_549[k]
                   + f_0 * il_1044[k];

        t_775[k] = -3.0 * gl_550[k]
                   + f_0 * il_1045[k];

        t_776[k] = -3.0 * gl_551[k]
                   + f_0 * il_1046[k];

        t_777[k] = -3.0 * gl_552[k]
                   + f_0 * il_1047[k];

        t_778[k] = -3.0 * gl_553[k]
                   + f_0 * il_1048[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, t_782, t_783, gl_554, gl_555, gl_556, gl_557, \
                         gl_558, il_1049, il_1050, il_1051, il_1052, \
                         il_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = -3.0 * gl_554[k]
                   + f_0 * il_1049[k];

        t_780[k] = -3.0 * gl_555[k]
                   + f_0 * il_1050[k];

        t_781[k] = -3.0 * gl_556[k]
                   + f_0 * il_1051[k];

        t_782[k] = -3.0 * gl_557[k]
                   + f_0 * il_1052[k];

        t_783[k] = -3.0 * gl_558[k]
                   + f_0 * il_1053[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, t_788, gl_559, gl_560, gl_561, gl_562, \
                         gl_563, il_1054, il_1055, il_1056, il_1057, \
                         il_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = -3.0 * gl_559[k]
                   + f_0 * il_1054[k];

        t_785[k] = -3.0 * gl_560[k]
                   + f_0 * il_1055[k];

        t_786[k] = -3.0 * gl_561[k]
                   + f_0 * il_1056[k];

        t_787[k] = -3.0 * gl_562[k]
                   + f_0 * il_1057[k];

        t_788[k] = -3.0 * gl_563[k]
                   + f_0 * il_1058[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, t_792, t_793, gl_564, gl_565, gl_566, gl_567, \
                         gl_568, il_1059, il_1060, il_1061, il_1062, \
                         il_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = -3.0 * gl_564[k]
                   + f_0 * il_1059[k];

        t_790[k] = -3.0 * gl_565[k]
                   + f_0 * il_1060[k];

        t_791[k] = -3.0 * gl_566[k]
                   + f_0 * il_1061[k];

        t_792[k] = -3.0 * gl_567[k]
                   + f_0 * il_1062[k];

        t_793[k] = -3.0 * gl_568[k]
                   + f_0 * il_1063[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, t_798, gl_569, gl_570, gl_571, gl_572, \
                         gl_573, il_1064, il_1065, il_1066, il_1067, \
                         il_1068 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = -3.0 * gl_569[k]
                   + f_0 * il_1064[k];

        t_795[k] = -3.0 * gl_570[k]
                   + f_0 * il_1065[k];

        t_796[k] = -3.0 * gl_571[k]
                   + f_0 * il_1066[k];

        t_797[k] = -3.0 * gl_572[k]
                   + f_0 * il_1067[k];

        t_798[k] = -3.0 * gl_573[k]
                   + f_0 * il_1068[k];
    }

#pragma omp simd aligned(t_799, t_800, t_801, t_802, t_803, gl_574, gl_575, gl_576, gl_577, \
                         gl_578, il_1069, il_1070, il_1071, il_1072, \
                         il_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = -3.0 * gl_574[k]
                   + f_0 * il_1069[k];

        t_800[k] = -3.0 * gl_575[k]
                   + f_0 * il_1070[k];

        t_801[k] = -3.0 * gl_576[k]
                   + f_0 * il_1071[k];

        t_802[k] = -3.0 * gl_577[k]
                   + f_0 * il_1072[k];

        t_803[k] = -3.0 * gl_578[k]
                   + f_0 * il_1073[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, t_807, t_808, gl_579, gl_580, gl_581, gl_582, \
                         gl_583, il_1074, il_1075, il_1076, il_1077, \
                         il_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = -3.0 * gl_579[k]
                   + f_0 * il_1074[k];

        t_805[k] = -3.0 * gl_580[k]
                   + f_0 * il_1075[k];

        t_806[k] = -3.0 * gl_581[k]
                   + f_0 * il_1076[k];

        t_807[k] = -3.0 * gl_582[k]
                   + f_0 * il_1077[k];

        t_808[k] = -3.0 * gl_583[k]
                   + f_0 * il_1078[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, t_812, t_813, gl_584, gl_585, gl_586, gl_587, \
                         gl_588, il_1079, il_1080, il_1081, il_1082, \
                         il_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = -3.0 * gl_584[k]
                   + f_0 * il_1079[k];

        t_810[k] = -2.0 * gl_585[k]
                   + f_0 * il_1080[k];

        t_811[k] = -2.0 * gl_586[k]
                   + f_0 * il_1081[k];

        t_812[k] = -2.0 * gl_587[k]
                   + f_0 * il_1082[k];

        t_813[k] = -2.0 * gl_588[k]
                   + f_0 * il_1083[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, t_817, t_818, gl_589, gl_590, gl_591, gl_592, \
                         gl_593, il_1084, il_1085, il_1086, il_1087, \
                         il_1088 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = -2.0 * gl_589[k]
                   + f_0 * il_1084[k];

        t_815[k] = -2.0 * gl_590[k]
                   + f_0 * il_1085[k];

        t_816[k] = -2.0 * gl_591[k]
                   + f_0 * il_1086[k];

        t_817[k] = -2.0 * gl_592[k]
                   + f_0 * il_1087[k];

        t_818[k] = -2.0 * gl_593[k]
                   + f_0 * il_1088[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, t_822, t_823, gl_594, gl_595, gl_596, gl_597, \
                         gl_598, il_1089, il_1090, il_1091, il_1092, \
                         il_1093 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = -2.0 * gl_594[k]
                   + f_0 * il_1089[k];

        t_820[k] = -2.0 * gl_595[k]
                   + f_0 * il_1090[k];

        t_821[k] = -2.0 * gl_596[k]
                   + f_0 * il_1091[k];

        t_822[k] = -2.0 * gl_597[k]
                   + f_0 * il_1092[k];

        t_823[k] = -2.0 * gl_598[k]
                   + f_0 * il_1093[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_1_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gl_599 = buffer.data(gl + 599);
    const auto *gl_600 = buffer.data(gl + 600);
    const auto *gl_601 = buffer.data(gl + 601);
    const auto *gl_602 = buffer.data(gl + 602);
    const auto *gl_603 = buffer.data(gl + 603);
    const auto *gl_604 = buffer.data(gl + 604);
    const auto *gl_605 = buffer.data(gl + 605);
    const auto *gl_606 = buffer.data(gl + 606);
    const auto *gl_607 = buffer.data(gl + 607);
    const auto *gl_608 = buffer.data(gl + 608);
    const auto *gl_609 = buffer.data(gl + 609);
    const auto *gl_610 = buffer.data(gl + 610);
    const auto *gl_611 = buffer.data(gl + 611);
    const auto *gl_612 = buffer.data(gl + 612);
    const auto *gl_613 = buffer.data(gl + 613);
    const auto *gl_614 = buffer.data(gl + 614);
    const auto *gl_615 = buffer.data(gl + 615);
    const auto *gl_616 = buffer.data(gl + 616);
    const auto *gl_617 = buffer.data(gl + 617);
    const auto *gl_618 = buffer.data(gl + 618);
    const auto *gl_619 = buffer.data(gl + 619);
    const auto *gl_620 = buffer.data(gl + 620);
    const auto *gl_621 = buffer.data(gl + 621);
    const auto *gl_622 = buffer.data(gl + 622);
    const auto *gl_623 = buffer.data(gl + 623);
    const auto *gl_624 = buffer.data(gl + 624);
    const auto *gl_625 = buffer.data(gl + 625);
    const auto *gl_626 = buffer.data(gl + 626);
    const auto *gl_627 = buffer.data(gl + 627);
    const auto *gl_628 = buffer.data(gl + 628);
    const auto *gl_629 = buffer.data(gl + 629);
    const auto *gl_630 = buffer.data(gl + 630);
    const auto *gl_631 = buffer.data(gl + 631);
    const auto *gl_632 = buffer.data(gl + 632);
    const auto *gl_633 = buffer.data(gl + 633);
    const auto *gl_634 = buffer.data(gl + 634);
    const auto *gl_635 = buffer.data(gl + 635);
    const auto *gl_636 = buffer.data(gl + 636);
    const auto *gl_637 = buffer.data(gl + 637);
    const auto *gl_638 = buffer.data(gl + 638);
    const auto *gl_639 = buffer.data(gl + 639);
    const auto *gl_640 = buffer.data(gl + 640);
    const auto *gl_641 = buffer.data(gl + 641);
    const auto *gl_642 = buffer.data(gl + 642);
    const auto *gl_643 = buffer.data(gl + 643);
    const auto *gl_644 = buffer.data(gl + 644);
    const auto *gl_645 = buffer.data(gl + 645);
    const auto *gl_646 = buffer.data(gl + 646);
    const auto *gl_647 = buffer.data(gl + 647);
    const auto *gl_648 = buffer.data(gl + 648);
    const auto *gl_649 = buffer.data(gl + 649);
    const auto *gl_650 = buffer.data(gl + 650);
    const auto *gl_651 = buffer.data(gl + 651);
    const auto *gl_652 = buffer.data(gl + 652);
    const auto *gl_653 = buffer.data(gl + 653);
    const auto *gl_654 = buffer.data(gl + 654);
    const auto *gl_655 = buffer.data(gl + 655);
    const auto *gl_656 = buffer.data(gl + 656);
    const auto *gl_657 = buffer.data(gl + 657);
    const auto *gl_658 = buffer.data(gl + 658);
    const auto *gl_659 = buffer.data(gl + 659);
    const auto *gl_660 = buffer.data(gl + 660);
    const auto *gl_661 = buffer.data(gl + 661);
    const auto *gl_662 = buffer.data(gl + 662);
    const auto *gl_663 = buffer.data(gl + 663);
    const auto *gl_664 = buffer.data(gl + 664);
    const auto *gl_665 = buffer.data(gl + 665);
    const auto *gl_666 = buffer.data(gl + 666);
    const auto *gl_667 = buffer.data(gl + 667);
    const auto *gl_668 = buffer.data(gl + 668);
    const auto *gl_669 = buffer.data(gl + 669);
    const auto *gl_670 = buffer.data(gl + 670);
    const auto *gl_671 = buffer.data(gl + 671);
    const auto *gl_672 = buffer.data(gl + 672);
    const auto *gl_673 = buffer.data(gl + 673);
    const auto *gl_674 = buffer.data(gl + 674);

    const auto *il_1094 = buffer.data(il + 1094);
    const auto *il_1095 = buffer.data(il + 1095);
    const auto *il_1096 = buffer.data(il + 1096);
    const auto *il_1097 = buffer.data(il + 1097);
    const auto *il_1098 = buffer.data(il + 1098);
    const auto *il_1099 = buffer.data(il + 1099);
    const auto *il_1100 = buffer.data(il + 1100);
    const auto *il_1101 = buffer.data(il + 1101);
    const auto *il_1102 = buffer.data(il + 1102);
    const auto *il_1103 = buffer.data(il + 1103);
    const auto *il_1104 = buffer.data(il + 1104);
    const auto *il_1105 = buffer.data(il + 1105);
    const auto *il_1106 = buffer.data(il + 1106);
    const auto *il_1107 = buffer.data(il + 1107);
    const auto *il_1108 = buffer.data(il + 1108);
    const auto *il_1109 = buffer.data(il + 1109);
    const auto *il_1110 = buffer.data(il + 1110);
    const auto *il_1111 = buffer.data(il + 1111);
    const auto *il_1112 = buffer.data(il + 1112);
    const auto *il_1113 = buffer.data(il + 1113);
    const auto *il_1114 = buffer.data(il + 1114);
    const auto *il_1115 = buffer.data(il + 1115);
    const auto *il_1116 = buffer.data(il + 1116);
    const auto *il_1117 = buffer.data(il + 1117);
    const auto *il_1118 = buffer.data(il + 1118);
    const auto *il_1119 = buffer.data(il + 1119);
    const auto *il_1120 = buffer.data(il + 1120);
    const auto *il_1121 = buffer.data(il + 1121);
    const auto *il_1122 = buffer.data(il + 1122);
    const auto *il_1123 = buffer.data(il + 1123);
    const auto *il_1124 = buffer.data(il + 1124);
    const auto *il_1125 = buffer.data(il + 1125);
    const auto *il_1126 = buffer.data(il + 1126);
    const auto *il_1127 = buffer.data(il + 1127);
    const auto *il_1128 = buffer.data(il + 1128);
    const auto *il_1129 = buffer.data(il + 1129);
    const auto *il_1130 = buffer.data(il + 1130);
    const auto *il_1131 = buffer.data(il + 1131);
    const auto *il_1132 = buffer.data(il + 1132);
    const auto *il_1133 = buffer.data(il + 1133);
    const auto *il_1134 = buffer.data(il + 1134);
    const auto *il_1135 = buffer.data(il + 1135);
    const auto *il_1136 = buffer.data(il + 1136);
    const auto *il_1137 = buffer.data(il + 1137);
    const auto *il_1138 = buffer.data(il + 1138);
    const auto *il_1139 = buffer.data(il + 1139);
    const auto *il_1140 = buffer.data(il + 1140);
    const auto *il_1141 = buffer.data(il + 1141);
    const auto *il_1142 = buffer.data(il + 1142);
    const auto *il_1143 = buffer.data(il + 1143);
    const auto *il_1144 = buffer.data(il + 1144);
    const auto *il_1145 = buffer.data(il + 1145);
    const auto *il_1146 = buffer.data(il + 1146);
    const auto *il_1147 = buffer.data(il + 1147);
    const auto *il_1148 = buffer.data(il + 1148);
    const auto *il_1149 = buffer.data(il + 1149);
    const auto *il_1150 = buffer.data(il + 1150);
    const auto *il_1151 = buffer.data(il + 1151);
    const auto *il_1152 = buffer.data(il + 1152);
    const auto *il_1153 = buffer.data(il + 1153);
    const auto *il_1154 = buffer.data(il + 1154);
    const auto *il_1155 = buffer.data(il + 1155);
    const auto *il_1156 = buffer.data(il + 1156);
    const auto *il_1157 = buffer.data(il + 1157);
    const auto *il_1158 = buffer.data(il + 1158);
    const auto *il_1159 = buffer.data(il + 1159);
    const auto *il_1160 = buffer.data(il + 1160);
    const auto *il_1161 = buffer.data(il + 1161);
    const auto *il_1162 = buffer.data(il + 1162);
    const auto *il_1163 = buffer.data(il + 1163);
    const auto *il_1164 = buffer.data(il + 1164);
    const auto *il_1165 = buffer.data(il + 1165);
    const auto *il_1166 = buffer.data(il + 1166);
    const auto *il_1167 = buffer.data(il + 1167);
    const auto *il_1168 = buffer.data(il + 1168);
    const auto *il_1169 = buffer.data(il + 1169);
    const auto *il_1170 = buffer.data(il + 1170);
    const auto *il_1171 = buffer.data(il + 1171);
    const auto *il_1172 = buffer.data(il + 1172);
    const auto *il_1173 = buffer.data(il + 1173);
    const auto *il_1174 = buffer.data(il + 1174);
    const auto *il_1175 = buffer.data(il + 1175);
    const auto *il_1176 = buffer.data(il + 1176);
    const auto *il_1177 = buffer.data(il + 1177);
    const auto *il_1178 = buffer.data(il + 1178);
    const auto *il_1179 = buffer.data(il + 1179);
    const auto *il_1180 = buffer.data(il + 1180);
    const auto *il_1181 = buffer.data(il + 1181);
    const auto *il_1182 = buffer.data(il + 1182);
    const auto *il_1183 = buffer.data(il + 1183);
    const auto *il_1184 = buffer.data(il + 1184);
    const auto *il_1185 = buffer.data(il + 1185);
    const auto *il_1186 = buffer.data(il + 1186);
    const auto *il_1187 = buffer.data(il + 1187);
    const auto *il_1188 = buffer.data(il + 1188);
    const auto *il_1189 = buffer.data(il + 1189);
    const auto *il_1190 = buffer.data(il + 1190);
    const auto *il_1191 = buffer.data(il + 1191);
    const auto *il_1192 = buffer.data(il + 1192);
    const auto *il_1193 = buffer.data(il + 1193);
    const auto *il_1194 = buffer.data(il + 1194);
    const auto *il_1195 = buffer.data(il + 1195);
    const auto *il_1196 = buffer.data(il + 1196);
    const auto *il_1197 = buffer.data(il + 1197);
    const auto *il_1198 = buffer.data(il + 1198);
    const auto *il_1199 = buffer.data(il + 1199);
    const auto *il_1200 = buffer.data(il + 1200);
    const auto *il_1201 = buffer.data(il + 1201);
    const auto *il_1202 = buffer.data(il + 1202);
    const auto *il_1203 = buffer.data(il + 1203);
    const auto *il_1204 = buffer.data(il + 1204);
    const auto *il_1205 = buffer.data(il + 1205);
    const auto *il_1206 = buffer.data(il + 1206);
    const auto *il_1207 = buffer.data(il + 1207);
    const auto *il_1208 = buffer.data(il + 1208);
    const auto *il_1209 = buffer.data(il + 1209);
    const auto *il_1210 = buffer.data(il + 1210);
    const auto *il_1211 = buffer.data(il + 1211);
    const auto *il_1212 = buffer.data(il + 1212);
    const auto *il_1213 = buffer.data(il + 1213);
    const auto *il_1214 = buffer.data(il + 1214);

#pragma omp simd aligned(t_824, t_825, t_826, t_827, t_828, gl_599, gl_600, gl_601, gl_602, \
                         gl_603, il_1094, il_1095, il_1096, il_1097, \
                         il_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = -2.0 * gl_599[k]
                   + f_0 * il_1094[k];

        t_825[k] = -2.0 * gl_600[k]
                   + f_0 * il_1095[k];

        t_826[k] = -2.0 * gl_601[k]
                   + f_0 * il_1096[k];

        t_827[k] = -2.0 * gl_602[k]
                   + f_0 * il_1097[k];

        t_828[k] = -2.0 * gl_603[k]
                   + f_0 * il_1098[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, t_833, gl_604, gl_605, gl_606, gl_607, \
                         gl_608, il_1099, il_1100, il_1101, il_1102, \
                         il_1103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = -2.0 * gl_604[k]
                   + f_0 * il_1099[k];

        t_830[k] = -2.0 * gl_605[k]
                   + f_0 * il_1100[k];

        t_831[k] = -2.0 * gl_606[k]
                   + f_0 * il_1101[k];

        t_832[k] = -2.0 * gl_607[k]
                   + f_0 * il_1102[k];

        t_833[k] = -2.0 * gl_608[k]
                   + f_0 * il_1103[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, t_837, t_838, gl_609, gl_610, gl_611, gl_612, \
                         gl_613, il_1104, il_1105, il_1106, il_1107, \
                         il_1108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = -2.0 * gl_609[k]
                   + f_0 * il_1104[k];

        t_835[k] = -2.0 * gl_610[k]
                   + f_0 * il_1105[k];

        t_836[k] = -2.0 * gl_611[k]
                   + f_0 * il_1106[k];

        t_837[k] = -2.0 * gl_612[k]
                   + f_0 * il_1107[k];

        t_838[k] = -2.0 * gl_613[k]
                   + f_0 * il_1108[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, t_842, t_843, gl_614, gl_615, gl_616, gl_617, \
                         gl_618, il_1109, il_1110, il_1111, il_1112, \
                         il_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = -2.0 * gl_614[k]
                   + f_0 * il_1109[k];

        t_840[k] = -2.0 * gl_615[k]
                   + f_0 * il_1110[k];

        t_841[k] = -2.0 * gl_616[k]
                   + f_0 * il_1111[k];

        t_842[k] = -2.0 * gl_617[k]
                   + f_0 * il_1112[k];

        t_843[k] = -2.0 * gl_618[k]
                   + f_0 * il_1113[k];
    }

#pragma omp simd aligned(t_844, t_845, t_846, t_847, t_848, gl_619, gl_620, gl_621, gl_622, \
                         gl_623, il_1114, il_1115, il_1116, il_1117, \
                         il_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_844[k] = -2.0 * gl_619[k]
                   + f_0 * il_1114[k];

        t_845[k] = -2.0 * gl_620[k]
                   + f_0 * il_1115[k];

        t_846[k] = -2.0 * gl_621[k]
                   + f_0 * il_1116[k];

        t_847[k] = -2.0 * gl_622[k]
                   + f_0 * il_1117[k];

        t_848[k] = -2.0 * gl_623[k]
                   + f_0 * il_1118[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, t_853, gl_624, gl_625, gl_626, gl_627, \
                         gl_628, il_1119, il_1120, il_1121, il_1122, \
                         il_1123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = -2.0 * gl_624[k]
                   + f_0 * il_1119[k];

        t_850[k] = -2.0 * gl_625[k]
                   + f_0 * il_1120[k];

        t_851[k] = -2.0 * gl_626[k]
                   + f_0 * il_1121[k];

        t_852[k] = -2.0 * gl_627[k]
                   + f_0 * il_1122[k];

        t_853[k] = -2.0 * gl_628[k]
                   + f_0 * il_1123[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, t_858, gl_629, gl_630, gl_631, gl_632, \
                         gl_633, il_1124, il_1125, il_1126, il_1127, \
                         il_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = -2.0 * gl_629[k]
                   + f_0 * il_1124[k];

        t_855[k] = -gl_630[k]
                   + f_0 * il_1125[k];

        t_856[k] = -gl_631[k]
                   + f_0 * il_1126[k];

        t_857[k] = -gl_632[k]
                   + f_0 * il_1127[k];

        t_858[k] = -gl_633[k]
                   + f_0 * il_1128[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, t_862, t_863, gl_634, gl_635, gl_636, gl_637, \
                         gl_638, il_1129, il_1130, il_1131, il_1132, \
                         il_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = -gl_634[k]
                   + f_0 * il_1129[k];

        t_860[k] = -gl_635[k]
                   + f_0 * il_1130[k];

        t_861[k] = -gl_636[k]
                   + f_0 * il_1131[k];

        t_862[k] = -gl_637[k]
                   + f_0 * il_1132[k];

        t_863[k] = -gl_638[k]
                   + f_0 * il_1133[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, t_867, t_868, gl_639, gl_640, gl_641, gl_642, \
                         gl_643, il_1134, il_1135, il_1136, il_1137, \
                         il_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = -gl_639[k]
                   + f_0 * il_1134[k];

        t_865[k] = -gl_640[k]
                   + f_0 * il_1135[k];

        t_866[k] = -gl_641[k]
                   + f_0 * il_1136[k];

        t_867[k] = -gl_642[k]
                   + f_0 * il_1137[k];

        t_868[k] = -gl_643[k]
                   + f_0 * il_1138[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, t_873, gl_644, gl_645, gl_646, gl_647, \
                         gl_648, il_1139, il_1140, il_1141, il_1142, \
                         il_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = -gl_644[k]
                   + f_0 * il_1139[k];

        t_870[k] = -gl_645[k]
                   + f_0 * il_1140[k];

        t_871[k] = -gl_646[k]
                   + f_0 * il_1141[k];

        t_872[k] = -gl_647[k]
                   + f_0 * il_1142[k];

        t_873[k] = -gl_648[k]
                   + f_0 * il_1143[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, t_878, gl_649, gl_650, gl_651, gl_652, \
                         gl_653, il_1144, il_1145, il_1146, il_1147, \
                         il_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = -gl_649[k]
                   + f_0 * il_1144[k];

        t_875[k] = -gl_650[k]
                   + f_0 * il_1145[k];

        t_876[k] = -gl_651[k]
                   + f_0 * il_1146[k];

        t_877[k] = -gl_652[k]
                   + f_0 * il_1147[k];

        t_878[k] = -gl_653[k]
                   + f_0 * il_1148[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, t_882, t_883, gl_654, gl_655, gl_656, gl_657, \
                         gl_658, il_1149, il_1150, il_1151, il_1152, \
                         il_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = -gl_654[k]
                   + f_0 * il_1149[k];

        t_880[k] = -gl_655[k]
                   + f_0 * il_1150[k];

        t_881[k] = -gl_656[k]
                   + f_0 * il_1151[k];

        t_882[k] = -gl_657[k]
                   + f_0 * il_1152[k];

        t_883[k] = -gl_658[k]
                   + f_0 * il_1153[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, t_888, gl_659, gl_660, gl_661, gl_662, \
                         gl_663, il_1154, il_1155, il_1156, il_1157, \
                         il_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = -gl_659[k]
                   + f_0 * il_1154[k];

        t_885[k] = -gl_660[k]
                   + f_0 * il_1155[k];

        t_886[k] = -gl_661[k]
                   + f_0 * il_1156[k];

        t_887[k] = -gl_662[k]
                   + f_0 * il_1157[k];

        t_888[k] = -gl_663[k]
                   + f_0 * il_1158[k];
    }

#pragma omp simd aligned(t_889, t_890, t_891, t_892, t_893, gl_664, gl_665, gl_666, gl_667, \
                         gl_668, il_1159, il_1160, il_1161, il_1162, \
                         il_1163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = -gl_664[k]
                   + f_0 * il_1159[k];

        t_890[k] = -gl_665[k]
                   + f_0 * il_1160[k];

        t_891[k] = -gl_666[k]
                   + f_0 * il_1161[k];

        t_892[k] = -gl_667[k]
                   + f_0 * il_1162[k];

        t_893[k] = -gl_668[k]
                   + f_0 * il_1163[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, t_898, gl_669, gl_670, gl_671, gl_672, \
                         gl_673, il_1164, il_1165, il_1166, il_1167, \
                         il_1168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = -gl_669[k]
                   + f_0 * il_1164[k];

        t_895[k] = -gl_670[k]
                   + f_0 * il_1165[k];

        t_896[k] = -gl_671[k]
                   + f_0 * il_1166[k];

        t_897[k] = -gl_672[k]
                   + f_0 * il_1167[k];

        t_898[k] = -gl_673[k]
                   + f_0 * il_1168[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, t_903, t_904, t_905, gl_674, il_1169, \
                         il_1170, il_1171, il_1172, il_1173, il_1174, \
                         il_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = -gl_674[k]
                   + f_0 * il_1169[k];

        t_900[k] = f_0 * il_1170[k];

        t_901[k] = f_0 * il_1171[k];

        t_902[k] = f_0 * il_1172[k];

        t_903[k] = f_0 * il_1173[k];

        t_904[k] = f_0 * il_1174[k];

        t_905[k] = f_0 * il_1175[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, t_909, t_910, t_911, t_912, t_913, il_1176, \
                         il_1177, il_1178, il_1179, il_1180, il_1181, il_1182, \
                         il_1183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_0 * il_1176[k];

        t_907[k] = f_0 * il_1177[k];

        t_908[k] = f_0 * il_1178[k];

        t_909[k] = f_0 * il_1179[k];

        t_910[k] = f_0 * il_1180[k];

        t_911[k] = f_0 * il_1181[k];

        t_912[k] = f_0 * il_1182[k];

        t_913[k] = f_0 * il_1183[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, t_918, t_919, t_920, t_921, il_1184, \
                         il_1185, il_1186, il_1187, il_1188, il_1189, il_1190, \
                         il_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_0 * il_1184[k];

        t_915[k] = f_0 * il_1185[k];

        t_916[k] = f_0 * il_1186[k];

        t_917[k] = f_0 * il_1187[k];

        t_918[k] = f_0 * il_1188[k];

        t_919[k] = f_0 * il_1189[k];

        t_920[k] = f_0 * il_1190[k];

        t_921[k] = f_0 * il_1191[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, t_926, t_927, t_928, t_929, il_1192, \
                         il_1193, il_1194, il_1195, il_1196, il_1197, il_1198, \
                         il_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_0 * il_1192[k];

        t_923[k] = f_0 * il_1193[k];

        t_924[k] = f_0 * il_1194[k];

        t_925[k] = f_0 * il_1195[k];

        t_926[k] = f_0 * il_1196[k];

        t_927[k] = f_0 * il_1197[k];

        t_928[k] = f_0 * il_1198[k];

        t_929[k] = f_0 * il_1199[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, t_935, t_936, t_937, il_1200, \
                         il_1201, il_1202, il_1203, il_1204, il_1205, il_1206, \
                         il_1207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_0 * il_1200[k];

        t_931[k] = f_0 * il_1201[k];

        t_932[k] = f_0 * il_1202[k];

        t_933[k] = f_0 * il_1203[k];

        t_934[k] = f_0 * il_1204[k];

        t_935[k] = f_0 * il_1205[k];

        t_936[k] = f_0 * il_1206[k];

        t_937[k] = f_0 * il_1207[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, t_942, t_943, t_944, il_1208, il_1209, \
                         il_1210, il_1211, il_1212, il_1213, il_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_0 * il_1208[k];

        t_939[k] = f_0 * il_1209[k];

        t_940[k] = f_0 * il_1210[k];

        t_941[k] = f_0 * il_1211[k];

        t_942[k] = f_0 * il_1212[k];

        t_943[k] = f_0 * il_1213[k];

        t_944[k] = f_0 * il_1214[k];
    }
}

auto
compute_prim_geom_10_hl_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t gl, const size_t il,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hl_electron_repulsion_1_piece0(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_1_piece1(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_1_piece2(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_1_piece3(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_1_piece4(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_1_piece5(buffer, target, gl, il, ncols, alpha);
}

static auto
compute_prim_geom_10_hl_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);
    const auto *gl_30 = buffer.data(gl + 30);
    const auto *gl_31 = buffer.data(gl + 31);
    const auto *gl_32 = buffer.data(gl + 32);
    const auto *gl_33 = buffer.data(gl + 33);
    const auto *gl_34 = buffer.data(gl + 34);
    const auto *gl_35 = buffer.data(gl + 35);
    const auto *gl_36 = buffer.data(gl + 36);
    const auto *gl_37 = buffer.data(gl + 37);
    const auto *gl_38 = buffer.data(gl + 38);
    const auto *gl_39 = buffer.data(gl + 39);
    const auto *gl_40 = buffer.data(gl + 40);
    const auto *gl_41 = buffer.data(gl + 41);
    const auto *gl_42 = buffer.data(gl + 42);
    const auto *gl_43 = buffer.data(gl + 43);
    const auto *gl_44 = buffer.data(gl + 44);
    const auto *gl_45 = buffer.data(gl + 45);
    const auto *gl_46 = buffer.data(gl + 46);
    const auto *gl_47 = buffer.data(gl + 47);
    const auto *gl_48 = buffer.data(gl + 48);
    const auto *gl_49 = buffer.data(gl + 49);
    const auto *gl_50 = buffer.data(gl + 50);
    const auto *gl_51 = buffer.data(gl + 51);
    const auto *gl_52 = buffer.data(gl + 52);
    const auto *gl_53 = buffer.data(gl + 53);
    const auto *gl_54 = buffer.data(gl + 54);
    const auto *gl_55 = buffer.data(gl + 55);
    const auto *gl_56 = buffer.data(gl + 56);
    const auto *gl_57 = buffer.data(gl + 57);
    const auto *gl_58 = buffer.data(gl + 58);
    const auto *gl_59 = buffer.data(gl + 59);

    const auto *il_90 = buffer.data(il + 90);
    const auto *il_91 = buffer.data(il + 91);
    const auto *il_92 = buffer.data(il + 92);
    const auto *il_93 = buffer.data(il + 93);
    const auto *il_94 = buffer.data(il + 94);
    const auto *il_95 = buffer.data(il + 95);
    const auto *il_96 = buffer.data(il + 96);
    const auto *il_97 = buffer.data(il + 97);
    const auto *il_98 = buffer.data(il + 98);
    const auto *il_99 = buffer.data(il + 99);
    const auto *il_100 = buffer.data(il + 100);
    const auto *il_101 = buffer.data(il + 101);
    const auto *il_102 = buffer.data(il + 102);
    const auto *il_103 = buffer.data(il + 103);
    const auto *il_104 = buffer.data(il + 104);
    const auto *il_105 = buffer.data(il + 105);
    const auto *il_106 = buffer.data(il + 106);
    const auto *il_107 = buffer.data(il + 107);
    const auto *il_108 = buffer.data(il + 108);
    const auto *il_109 = buffer.data(il + 109);
    const auto *il_110 = buffer.data(il + 110);
    const auto *il_111 = buffer.data(il + 111);
    const auto *il_112 = buffer.data(il + 112);
    const auto *il_113 = buffer.data(il + 113);
    const auto *il_114 = buffer.data(il + 114);
    const auto *il_115 = buffer.data(il + 115);
    const auto *il_116 = buffer.data(il + 116);
    const auto *il_117 = buffer.data(il + 117);
    const auto *il_118 = buffer.data(il + 118);
    const auto *il_119 = buffer.data(il + 119);
    const auto *il_120 = buffer.data(il + 120);
    const auto *il_121 = buffer.data(il + 121);
    const auto *il_122 = buffer.data(il + 122);
    const auto *il_123 = buffer.data(il + 123);
    const auto *il_124 = buffer.data(il + 124);
    const auto *il_125 = buffer.data(il + 125);
    const auto *il_126 = buffer.data(il + 126);
    const auto *il_127 = buffer.data(il + 127);
    const auto *il_128 = buffer.data(il + 128);
    const auto *il_129 = buffer.data(il + 129);
    const auto *il_130 = buffer.data(il + 130);
    const auto *il_131 = buffer.data(il + 131);
    const auto *il_132 = buffer.data(il + 132);
    const auto *il_133 = buffer.data(il + 133);
    const auto *il_134 = buffer.data(il + 134);
    const auto *il_180 = buffer.data(il + 180);
    const auto *il_181 = buffer.data(il + 181);
    const auto *il_182 = buffer.data(il + 182);
    const auto *il_183 = buffer.data(il + 183);
    const auto *il_184 = buffer.data(il + 184);
    const auto *il_185 = buffer.data(il + 185);
    const auto *il_186 = buffer.data(il + 186);
    const auto *il_187 = buffer.data(il + 187);
    const auto *il_188 = buffer.data(il + 188);
    const auto *il_189 = buffer.data(il + 189);
    const auto *il_190 = buffer.data(il + 190);
    const auto *il_191 = buffer.data(il + 191);
    const auto *il_192 = buffer.data(il + 192);
    const auto *il_193 = buffer.data(il + 193);
    const auto *il_194 = buffer.data(il + 194);
    const auto *il_195 = buffer.data(il + 195);
    const auto *il_196 = buffer.data(il + 196);
    const auto *il_197 = buffer.data(il + 197);
    const auto *il_198 = buffer.data(il + 198);
    const auto *il_199 = buffer.data(il + 199);
    const auto *il_200 = buffer.data(il + 200);
    const auto *il_201 = buffer.data(il + 201);
    const auto *il_202 = buffer.data(il + 202);
    const auto *il_203 = buffer.data(il + 203);
    const auto *il_204 = buffer.data(il + 204);
    const auto *il_205 = buffer.data(il + 205);
    const auto *il_206 = buffer.data(il + 206);
    const auto *il_207 = buffer.data(il + 207);
    const auto *il_208 = buffer.data(il + 208);
    const auto *il_209 = buffer.data(il + 209);
    const auto *il_210 = buffer.data(il + 210);
    const auto *il_211 = buffer.data(il + 211);
    const auto *il_212 = buffer.data(il + 212);
    const auto *il_213 = buffer.data(il + 213);
    const auto *il_214 = buffer.data(il + 214);
    const auto *il_215 = buffer.data(il + 215);
    const auto *il_216 = buffer.data(il + 216);
    const auto *il_217 = buffer.data(il + 217);
    const auto *il_218 = buffer.data(il + 218);
    const auto *il_219 = buffer.data(il + 219);
    const auto *il_220 = buffer.data(il + 220);
    const auto *il_221 = buffer.data(il + 221);
    const auto *il_222 = buffer.data(il + 222);
    const auto *il_223 = buffer.data(il + 223);
    const auto *il_224 = buffer.data(il + 224);
    const auto *il_225 = buffer.data(il + 225);
    const auto *il_226 = buffer.data(il + 226);
    const auto *il_227 = buffer.data(il + 227);
    const auto *il_228 = buffer.data(il + 228);
    const auto *il_229 = buffer.data(il + 229);
    const auto *il_230 = buffer.data(il + 230);
    const auto *il_231 = buffer.data(il + 231);
    const auto *il_232 = buffer.data(il + 232);
    const auto *il_233 = buffer.data(il + 233);
    const auto *il_234 = buffer.data(il + 234);
    const auto *il_235 = buffer.data(il + 235);
    const auto *il_236 = buffer.data(il + 236);
    const auto *il_237 = buffer.data(il + 237);
    const auto *il_238 = buffer.data(il + 238);
    const auto *il_239 = buffer.data(il + 239);
    const auto *il_240 = buffer.data(il + 240);
    const auto *il_241 = buffer.data(il + 241);
    const auto *il_242 = buffer.data(il + 242);
    const auto *il_243 = buffer.data(il + 243);
    const auto *il_244 = buffer.data(il + 244);
    const auto *il_245 = buffer.data(il + 245);
    const auto *il_246 = buffer.data(il + 246);
    const auto *il_247 = buffer.data(il + 247);
    const auto *il_248 = buffer.data(il + 248);
    const auto *il_249 = buffer.data(il + 249);
    const auto *il_250 = buffer.data(il + 250);
    const auto *il_251 = buffer.data(il + 251);
    const auto *il_252 = buffer.data(il + 252);
    const auto *il_253 = buffer.data(il + 253);
    const auto *il_254 = buffer.data(il + 254);
    const auto *il_255 = buffer.data(il + 255);
    const auto *il_256 = buffer.data(il + 256);
    const auto *il_257 = buffer.data(il + 257);
    const auto *il_258 = buffer.data(il + 258);
    const auto *il_259 = buffer.data(il + 259);
    const auto *il_260 = buffer.data(il + 260);
    const auto *il_261 = buffer.data(il + 261);
    const auto *il_262 = buffer.data(il + 262);
    const auto *il_263 = buffer.data(il + 263);
    const auto *il_264 = buffer.data(il + 264);
    const auto *il_265 = buffer.data(il + 265);
    const auto *il_266 = buffer.data(il + 266);
    const auto *il_267 = buffer.data(il + 267);
    const auto *il_268 = buffer.data(il + 268);
    const auto *il_269 = buffer.data(il + 269);
    const auto *il_315 = buffer.data(il + 315);
    const auto *il_316 = buffer.data(il + 316);
    const auto *il_317 = buffer.data(il + 317);
    const auto *il_318 = buffer.data(il + 318);
    const auto *il_319 = buffer.data(il + 319);
    const auto *il_320 = buffer.data(il + 320);
    const auto *il_321 = buffer.data(il + 321);
    const auto *il_322 = buffer.data(il + 322);
    const auto *il_323 = buffer.data(il + 323);
    const auto *il_324 = buffer.data(il + 324);
    const auto *il_325 = buffer.data(il + 325);
    const auto *il_326 = buffer.data(il + 326);
    const auto *il_327 = buffer.data(il + 327);
    const auto *il_328 = buffer.data(il + 328);
    const auto *il_329 = buffer.data(il + 329);
    const auto *il_330 = buffer.data(il + 330);
    const auto *il_331 = buffer.data(il + 331);
    const auto *il_332 = buffer.data(il + 332);
    const auto *il_333 = buffer.data(il + 333);
    const auto *il_334 = buffer.data(il + 334);
    const auto *il_335 = buffer.data(il + 335);
    const auto *il_336 = buffer.data(il + 336);
    const auto *il_337 = buffer.data(il + 337);
    const auto *il_338 = buffer.data(il + 338);
    const auto *il_339 = buffer.data(il + 339);
    const auto *il_340 = buffer.data(il + 340);
    const auto *il_341 = buffer.data(il + 341);
    const auto *il_342 = buffer.data(il + 342);
    const auto *il_343 = buffer.data(il + 343);
    const auto *il_344 = buffer.data(il + 344);
    const auto *il_345 = buffer.data(il + 345);
    const auto *il_346 = buffer.data(il + 346);
    const auto *il_347 = buffer.data(il + 347);
    const auto *il_348 = buffer.data(il + 348);
    const auto *il_349 = buffer.data(il + 349);
    const auto *il_350 = buffer.data(il + 350);
    const auto *il_351 = buffer.data(il + 351);
    const auto *il_352 = buffer.data(il + 352);
    const auto *il_353 = buffer.data(il + 353);
    const auto *il_354 = buffer.data(il + 354);
    const auto *il_355 = buffer.data(il + 355);
    const auto *il_356 = buffer.data(il + 356);
    const auto *il_357 = buffer.data(il + 357);
    const auto *il_358 = buffer.data(il + 358);
    const auto *il_359 = buffer.data(il + 359);
    const auto *il_360 = buffer.data(il + 360);
    const auto *il_361 = buffer.data(il + 361);
    const auto *il_362 = buffer.data(il + 362);
    const auto *il_363 = buffer.data(il + 363);
    const auto *il_364 = buffer.data(il + 364);
    const auto *il_365 = buffer.data(il + 365);
    const auto *il_366 = buffer.data(il + 366);
    const auto *il_367 = buffer.data(il + 367);
    const auto *il_368 = buffer.data(il + 368);
    const auto *il_369 = buffer.data(il + 369);
    const auto *il_370 = buffer.data(il + 370);
    const auto *il_371 = buffer.data(il + 371);
    const auto *il_372 = buffer.data(il + 372);
    const auto *il_373 = buffer.data(il + 373);
    const auto *il_374 = buffer.data(il + 374);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, il_90, il_91, il_92, il_93, \
                         il_94, il_95, il_96, il_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * il_90[k];

        t_1[k] = f_0 * il_91[k];

        t_2[k] = f_0 * il_92[k];

        t_3[k] = f_0 * il_93[k];

        t_4[k] = f_0 * il_94[k];

        t_5[k] = f_0 * il_95[k];

        t_6[k] = f_0 * il_96[k];

        t_7[k] = f_0 * il_97[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, il_98, il_99, il_100, \
                         il_101, il_102, il_103, il_104, il_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * il_98[k];

        t_9[k] = f_0 * il_99[k];

        t_10[k] = f_0 * il_100[k];

        t_11[k] = f_0 * il_101[k];

        t_12[k] = f_0 * il_102[k];

        t_13[k] = f_0 * il_103[k];

        t_14[k] = f_0 * il_104[k];

        t_15[k] = f_0 * il_105[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, il_106, il_107, \
                         il_108, il_109, il_110, il_111, il_112, \
                         il_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * il_106[k];

        t_17[k] = f_0 * il_107[k];

        t_18[k] = f_0 * il_108[k];

        t_19[k] = f_0 * il_109[k];

        t_20[k] = f_0 * il_110[k];

        t_21[k] = f_0 * il_111[k];

        t_22[k] = f_0 * il_112[k];

        t_23[k] = f_0 * il_113[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, il_114, il_115, \
                         il_116, il_117, il_118, il_119, il_120, \
                         il_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * il_114[k];

        t_25[k] = f_0 * il_115[k];

        t_26[k] = f_0 * il_116[k];

        t_27[k] = f_0 * il_117[k];

        t_28[k] = f_0 * il_118[k];

        t_29[k] = f_0 * il_119[k];

        t_30[k] = f_0 * il_120[k];

        t_31[k] = f_0 * il_121[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, il_122, il_123, \
                         il_124, il_125, il_126, il_127, il_128, \
                         il_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * il_122[k];

        t_33[k] = f_0 * il_123[k];

        t_34[k] = f_0 * il_124[k];

        t_35[k] = f_0 * il_125[k];

        t_36[k] = f_0 * il_126[k];

        t_37[k] = f_0 * il_127[k];

        t_38[k] = f_0 * il_128[k];

        t_39[k] = f_0 * il_129[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, il_130, il_131, \
                         il_132, il_133, il_134, il_180, il_181, \
                         il_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * il_130[k];

        t_41[k] = f_0 * il_131[k];

        t_42[k] = f_0 * il_132[k];

        t_43[k] = f_0 * il_133[k];

        t_44[k] = f_0 * il_134[k];

        t_45[k] = f_0 * il_180[k];

        t_46[k] = f_0 * il_181[k];

        t_47[k] = f_0 * il_182[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, il_183, il_184, \
                         il_185, il_186, il_187, il_188, il_189, \
                         il_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * il_183[k];

        t_49[k] = f_0 * il_184[k];

        t_50[k] = f_0 * il_185[k];

        t_51[k] = f_0 * il_186[k];

        t_52[k] = f_0 * il_187[k];

        t_53[k] = f_0 * il_188[k];

        t_54[k] = f_0 * il_189[k];

        t_55[k] = f_0 * il_190[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, il_191, il_192, \
                         il_193, il_194, il_195, il_196, il_197, \
                         il_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * il_191[k];

        t_57[k] = f_0 * il_192[k];

        t_58[k] = f_0 * il_193[k];

        t_59[k] = f_0 * il_194[k];

        t_60[k] = f_0 * il_195[k];

        t_61[k] = f_0 * il_196[k];

        t_62[k] = f_0 * il_197[k];

        t_63[k] = f_0 * il_198[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, il_199, il_200, \
                         il_201, il_202, il_203, il_204, il_205, \
                         il_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * il_199[k];

        t_65[k] = f_0 * il_200[k];

        t_66[k] = f_0 * il_201[k];

        t_67[k] = f_0 * il_202[k];

        t_68[k] = f_0 * il_203[k];

        t_69[k] = f_0 * il_204[k];

        t_70[k] = f_0 * il_205[k];

        t_71[k] = f_0 * il_206[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, t_78, t_79, il_207, il_208, \
                         il_209, il_210, il_211, il_212, il_213, \
                         il_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * il_207[k];

        t_73[k] = f_0 * il_208[k];

        t_74[k] = f_0 * il_209[k];

        t_75[k] = f_0 * il_210[k];

        t_76[k] = f_0 * il_211[k];

        t_77[k] = f_0 * il_212[k];

        t_78[k] = f_0 * il_213[k];

        t_79[k] = f_0 * il_214[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, t_86, t_87, il_215, il_216, \
                         il_217, il_218, il_219, il_220, il_221, \
                         il_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * il_215[k];

        t_81[k] = f_0 * il_216[k];

        t_82[k] = f_0 * il_217[k];

        t_83[k] = f_0 * il_218[k];

        t_84[k] = f_0 * il_219[k];

        t_85[k] = f_0 * il_220[k];

        t_86[k] = f_0 * il_221[k];

        t_87[k] = f_0 * il_222[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, t_93, gl_0, gl_1, gl_2, gl_3, il_223, \
                         il_224, il_225, il_226, il_227, il_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_0 * il_223[k];

        t_89[k] = f_0 * il_224[k];

        t_90[k] = -gl_0[k]
                  + f_0 * il_225[k];

        t_91[k] = -gl_1[k]
                  + f_0 * il_226[k];

        t_92[k] = -gl_2[k]
                  + f_0 * il_227[k];

        t_93[k] = -gl_3[k]
                  + f_0 * il_228[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, gl_4, gl_5, gl_6, gl_7, gl_8, il_229, \
                         il_230, il_231, il_232, il_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -gl_4[k]
                  + f_0 * il_229[k];

        t_95[k] = -gl_5[k]
                  + f_0 * il_230[k];

        t_96[k] = -gl_6[k]
                  + f_0 * il_231[k];

        t_97[k] = -gl_7[k]
                  + f_0 * il_232[k];

        t_98[k] = -gl_8[k]
                  + f_0 * il_233[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, gl_9, gl_10, gl_11, gl_12, gl_13, \
                         il_234, il_235, il_236, il_237, il_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -gl_9[k]
                  + f_0 * il_234[k];

        t_100[k] = -gl_10[k]
                   + f_0 * il_235[k];

        t_101[k] = -gl_11[k]
                   + f_0 * il_236[k];

        t_102[k] = -gl_12[k]
                   + f_0 * il_237[k];

        t_103[k] = -gl_13[k]
                   + f_0 * il_238[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, gl_14, gl_15, gl_16, gl_17, gl_18, \
                         il_239, il_240, il_241, il_242, il_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -gl_14[k]
                   + f_0 * il_239[k];

        t_105[k] = -gl_15[k]
                   + f_0 * il_240[k];

        t_106[k] = -gl_16[k]
                   + f_0 * il_241[k];

        t_107[k] = -gl_17[k]
                   + f_0 * il_242[k];

        t_108[k] = -gl_18[k]
                   + f_0 * il_243[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, gl_19, gl_20, gl_21, gl_22, gl_23, \
                         il_244, il_245, il_246, il_247, il_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -gl_19[k]
                   + f_0 * il_244[k];

        t_110[k] = -gl_20[k]
                   + f_0 * il_245[k];

        t_111[k] = -gl_21[k]
                   + f_0 * il_246[k];

        t_112[k] = -gl_22[k]
                   + f_0 * il_247[k];

        t_113[k] = -gl_23[k]
                   + f_0 * il_248[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, gl_24, gl_25, gl_26, gl_27, gl_28, \
                         il_249, il_250, il_251, il_252, il_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -gl_24[k]
                   + f_0 * il_249[k];

        t_115[k] = -gl_25[k]
                   + f_0 * il_250[k];

        t_116[k] = -gl_26[k]
                   + f_0 * il_251[k];

        t_117[k] = -gl_27[k]
                   + f_0 * il_252[k];

        t_118[k] = -gl_28[k]
                   + f_0 * il_253[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, gl_29, gl_30, gl_31, gl_32, gl_33, \
                         il_254, il_255, il_256, il_257, il_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -gl_29[k]
                   + f_0 * il_254[k];

        t_120[k] = -gl_30[k]
                   + f_0 * il_255[k];

        t_121[k] = -gl_31[k]
                   + f_0 * il_256[k];

        t_122[k] = -gl_32[k]
                   + f_0 * il_257[k];

        t_123[k] = -gl_33[k]
                   + f_0 * il_258[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, gl_34, gl_35, gl_36, gl_37, gl_38, \
                         il_259, il_260, il_261, il_262, il_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -gl_34[k]
                   + f_0 * il_259[k];

        t_125[k] = -gl_35[k]
                   + f_0 * il_260[k];

        t_126[k] = -gl_36[k]
                   + f_0 * il_261[k];

        t_127[k] = -gl_37[k]
                   + f_0 * il_262[k];

        t_128[k] = -gl_38[k]
                   + f_0 * il_263[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, gl_39, gl_40, gl_41, gl_42, gl_43, \
                         il_264, il_265, il_266, il_267, il_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -gl_39[k]
                   + f_0 * il_264[k];

        t_130[k] = -gl_40[k]
                   + f_0 * il_265[k];

        t_131[k] = -gl_41[k]
                   + f_0 * il_266[k];

        t_132[k] = -gl_42[k]
                   + f_0 * il_267[k];

        t_133[k] = -gl_43[k]
                   + f_0 * il_268[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, t_139, t_140, gl_44, il_269, \
                         il_315, il_316, il_317, il_318, il_319, \
                         il_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -gl_44[k]
                   + f_0 * il_269[k];

        t_135[k] = f_0 * il_315[k];

        t_136[k] = f_0 * il_316[k];

        t_137[k] = f_0 * il_317[k];

        t_138[k] = f_0 * il_318[k];

        t_139[k] = f_0 * il_319[k];

        t_140[k] = f_0 * il_320[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, t_146, t_147, t_148, il_321, \
                         il_322, il_323, il_324, il_325, il_326, il_327, \
                         il_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_0 * il_321[k];

        t_142[k] = f_0 * il_322[k];

        t_143[k] = f_0 * il_323[k];

        t_144[k] = f_0 * il_324[k];

        t_145[k] = f_0 * il_325[k];

        t_146[k] = f_0 * il_326[k];

        t_147[k] = f_0 * il_327[k];

        t_148[k] = f_0 * il_328[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, t_154, t_155, t_156, il_329, \
                         il_330, il_331, il_332, il_333, il_334, il_335, \
                         il_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * il_329[k];

        t_150[k] = f_0 * il_330[k];

        t_151[k] = f_0 * il_331[k];

        t_152[k] = f_0 * il_332[k];

        t_153[k] = f_0 * il_333[k];

        t_154[k] = f_0 * il_334[k];

        t_155[k] = f_0 * il_335[k];

        t_156[k] = f_0 * il_336[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, t_162, t_163, t_164, il_337, \
                         il_338, il_339, il_340, il_341, il_342, il_343, \
                         il_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_0 * il_337[k];

        t_158[k] = f_0 * il_338[k];

        t_159[k] = f_0 * il_339[k];

        t_160[k] = f_0 * il_340[k];

        t_161[k] = f_0 * il_341[k];

        t_162[k] = f_0 * il_342[k];

        t_163[k] = f_0 * il_343[k];

        t_164[k] = f_0 * il_344[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, t_170, t_171, t_172, il_345, \
                         il_346, il_347, il_348, il_349, il_350, il_351, \
                         il_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_0 * il_345[k];

        t_166[k] = f_0 * il_346[k];

        t_167[k] = f_0 * il_347[k];

        t_168[k] = f_0 * il_348[k];

        t_169[k] = f_0 * il_349[k];

        t_170[k] = f_0 * il_350[k];

        t_171[k] = f_0 * il_351[k];

        t_172[k] = f_0 * il_352[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, t_179, il_353, il_354, \
                         il_355, il_356, il_357, il_358, il_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_0 * il_353[k];

        t_174[k] = f_0 * il_354[k];

        t_175[k] = f_0 * il_355[k];

        t_176[k] = f_0 * il_356[k];

        t_177[k] = f_0 * il_357[k];

        t_178[k] = f_0 * il_358[k];

        t_179[k] = f_0 * il_359[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, gl_45, gl_46, gl_47, gl_48, gl_49, \
                         il_360, il_361, il_362, il_363, il_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -gl_45[k]
                   + f_0 * il_360[k];

        t_181[k] = -gl_46[k]
                   + f_0 * il_361[k];

        t_182[k] = -gl_47[k]
                   + f_0 * il_362[k];

        t_183[k] = -gl_48[k]
                   + f_0 * il_363[k];

        t_184[k] = -gl_49[k]
                   + f_0 * il_364[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, gl_50, gl_51, gl_52, gl_53, gl_54, \
                         il_365, il_366, il_367, il_368, il_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -gl_50[k]
                   + f_0 * il_365[k];

        t_186[k] = -gl_51[k]
                   + f_0 * il_366[k];

        t_187[k] = -gl_52[k]
                   + f_0 * il_367[k];

        t_188[k] = -gl_53[k]
                   + f_0 * il_368[k];

        t_189[k] = -gl_54[k]
                   + f_0 * il_369[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, gl_55, gl_56, gl_57, gl_58, gl_59, \
                         il_370, il_371, il_372, il_373, il_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -gl_55[k]
                   + f_0 * il_370[k];

        t_191[k] = -gl_56[k]
                   + f_0 * il_371[k];

        t_192[k] = -gl_57[k]
                   + f_0 * il_372[k];

        t_193[k] = -gl_58[k]
                   + f_0 * il_373[k];

        t_194[k] = -gl_59[k]
                   + f_0 * il_374[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_60 = buffer.data(gl + 60);
    const auto *gl_61 = buffer.data(gl + 61);
    const auto *gl_62 = buffer.data(gl + 62);
    const auto *gl_63 = buffer.data(gl + 63);
    const auto *gl_64 = buffer.data(gl + 64);
    const auto *gl_65 = buffer.data(gl + 65);
    const auto *gl_66 = buffer.data(gl + 66);
    const auto *gl_67 = buffer.data(gl + 67);
    const auto *gl_68 = buffer.data(gl + 68);
    const auto *gl_69 = buffer.data(gl + 69);
    const auto *gl_70 = buffer.data(gl + 70);
    const auto *gl_71 = buffer.data(gl + 71);
    const auto *gl_72 = buffer.data(gl + 72);
    const auto *gl_73 = buffer.data(gl + 73);
    const auto *gl_74 = buffer.data(gl + 74);
    const auto *gl_75 = buffer.data(gl + 75);
    const auto *gl_76 = buffer.data(gl + 76);
    const auto *gl_77 = buffer.data(gl + 77);
    const auto *gl_78 = buffer.data(gl + 78);
    const auto *gl_79 = buffer.data(gl + 79);
    const auto *gl_80 = buffer.data(gl + 80);
    const auto *gl_81 = buffer.data(gl + 81);
    const auto *gl_82 = buffer.data(gl + 82);
    const auto *gl_83 = buffer.data(gl + 83);
    const auto *gl_84 = buffer.data(gl + 84);
    const auto *gl_85 = buffer.data(gl + 85);
    const auto *gl_86 = buffer.data(gl + 86);
    const auto *gl_87 = buffer.data(gl + 87);
    const auto *gl_88 = buffer.data(gl + 88);
    const auto *gl_89 = buffer.data(gl + 89);
    const auto *gl_90 = buffer.data(gl + 90);
    const auto *gl_91 = buffer.data(gl + 91);
    const auto *gl_92 = buffer.data(gl + 92);
    const auto *gl_93 = buffer.data(gl + 93);
    const auto *gl_94 = buffer.data(gl + 94);
    const auto *gl_95 = buffer.data(gl + 95);
    const auto *gl_96 = buffer.data(gl + 96);
    const auto *gl_97 = buffer.data(gl + 97);
    const auto *gl_98 = buffer.data(gl + 98);
    const auto *gl_99 = buffer.data(gl + 99);
    const auto *gl_100 = buffer.data(gl + 100);
    const auto *gl_101 = buffer.data(gl + 101);
    const auto *gl_102 = buffer.data(gl + 102);
    const auto *gl_103 = buffer.data(gl + 103);
    const auto *gl_104 = buffer.data(gl + 104);
    const auto *gl_105 = buffer.data(gl + 105);
    const auto *gl_106 = buffer.data(gl + 106);
    const auto *gl_107 = buffer.data(gl + 107);
    const auto *gl_108 = buffer.data(gl + 108);
    const auto *gl_109 = buffer.data(gl + 109);
    const auto *gl_110 = buffer.data(gl + 110);
    const auto *gl_111 = buffer.data(gl + 111);
    const auto *gl_112 = buffer.data(gl + 112);
    const auto *gl_113 = buffer.data(gl + 113);
    const auto *gl_114 = buffer.data(gl + 114);
    const auto *gl_115 = buffer.data(gl + 115);
    const auto *gl_116 = buffer.data(gl + 116);
    const auto *gl_117 = buffer.data(gl + 117);
    const auto *gl_118 = buffer.data(gl + 118);
    const auto *gl_119 = buffer.data(gl + 119);
    const auto *gl_120 = buffer.data(gl + 120);
    const auto *gl_121 = buffer.data(gl + 121);
    const auto *gl_122 = buffer.data(gl + 122);
    const auto *gl_123 = buffer.data(gl + 123);
    const auto *gl_124 = buffer.data(gl + 124);
    const auto *gl_125 = buffer.data(gl + 125);
    const auto *gl_126 = buffer.data(gl + 126);
    const auto *gl_127 = buffer.data(gl + 127);
    const auto *gl_128 = buffer.data(gl + 128);
    const auto *gl_129 = buffer.data(gl + 129);
    const auto *gl_130 = buffer.data(gl + 130);
    const auto *gl_131 = buffer.data(gl + 131);
    const auto *gl_132 = buffer.data(gl + 132);
    const auto *gl_133 = buffer.data(gl + 133);
    const auto *gl_134 = buffer.data(gl + 134);
    const auto *gl_135 = buffer.data(gl + 135);
    const auto *gl_136 = buffer.data(gl + 136);
    const auto *gl_137 = buffer.data(gl + 137);
    const auto *gl_138 = buffer.data(gl + 138);
    const auto *gl_139 = buffer.data(gl + 139);
    const auto *gl_140 = buffer.data(gl + 140);
    const auto *gl_141 = buffer.data(gl + 141);
    const auto *gl_142 = buffer.data(gl + 142);
    const auto *gl_143 = buffer.data(gl + 143);
    const auto *gl_144 = buffer.data(gl + 144);
    const auto *gl_145 = buffer.data(gl + 145);
    const auto *gl_146 = buffer.data(gl + 146);
    const auto *gl_147 = buffer.data(gl + 147);
    const auto *gl_148 = buffer.data(gl + 148);
    const auto *gl_149 = buffer.data(gl + 149);
    const auto *gl_150 = buffer.data(gl + 150);
    const auto *gl_151 = buffer.data(gl + 151);
    const auto *gl_152 = buffer.data(gl + 152);
    const auto *gl_153 = buffer.data(gl + 153);
    const auto *gl_154 = buffer.data(gl + 154);
    const auto *gl_155 = buffer.data(gl + 155);
    const auto *gl_156 = buffer.data(gl + 156);
    const auto *gl_157 = buffer.data(gl + 157);
    const auto *gl_158 = buffer.data(gl + 158);
    const auto *gl_159 = buffer.data(gl + 159);
    const auto *gl_160 = buffer.data(gl + 160);
    const auto *gl_161 = buffer.data(gl + 161);
    const auto *gl_162 = buffer.data(gl + 162);
    const auto *gl_163 = buffer.data(gl + 163);
    const auto *gl_164 = buffer.data(gl + 164);
    const auto *gl_165 = buffer.data(gl + 165);
    const auto *gl_166 = buffer.data(gl + 166);
    const auto *gl_167 = buffer.data(gl + 167);
    const auto *gl_168 = buffer.data(gl + 168);
    const auto *gl_169 = buffer.data(gl + 169);
    const auto *gl_170 = buffer.data(gl + 170);
    const auto *gl_171 = buffer.data(gl + 171);
    const auto *gl_172 = buffer.data(gl + 172);
    const auto *gl_173 = buffer.data(gl + 173);
    const auto *gl_174 = buffer.data(gl + 174);
    const auto *gl_175 = buffer.data(gl + 175);
    const auto *gl_176 = buffer.data(gl + 176);

    const auto *il_375 = buffer.data(il + 375);
    const auto *il_376 = buffer.data(il + 376);
    const auto *il_377 = buffer.data(il + 377);
    const auto *il_378 = buffer.data(il + 378);
    const auto *il_379 = buffer.data(il + 379);
    const auto *il_380 = buffer.data(il + 380);
    const auto *il_381 = buffer.data(il + 381);
    const auto *il_382 = buffer.data(il + 382);
    const auto *il_383 = buffer.data(il + 383);
    const auto *il_384 = buffer.data(il + 384);
    const auto *il_385 = buffer.data(il + 385);
    const auto *il_386 = buffer.data(il + 386);
    const auto *il_387 = buffer.data(il + 387);
    const auto *il_388 = buffer.data(il + 388);
    const auto *il_389 = buffer.data(il + 389);
    const auto *il_390 = buffer.data(il + 390);
    const auto *il_391 = buffer.data(il + 391);
    const auto *il_392 = buffer.data(il + 392);
    const auto *il_393 = buffer.data(il + 393);
    const auto *il_394 = buffer.data(il + 394);
    const auto *il_395 = buffer.data(il + 395);
    const auto *il_396 = buffer.data(il + 396);
    const auto *il_397 = buffer.data(il + 397);
    const auto *il_398 = buffer.data(il + 398);
    const auto *il_399 = buffer.data(il + 399);
    const auto *il_400 = buffer.data(il + 400);
    const auto *il_401 = buffer.data(il + 401);
    const auto *il_402 = buffer.data(il + 402);
    const auto *il_403 = buffer.data(il + 403);
    const auto *il_404 = buffer.data(il + 404);
    const auto *il_405 = buffer.data(il + 405);
    const auto *il_406 = buffer.data(il + 406);
    const auto *il_407 = buffer.data(il + 407);
    const auto *il_408 = buffer.data(il + 408);
    const auto *il_409 = buffer.data(il + 409);
    const auto *il_410 = buffer.data(il + 410);
    const auto *il_411 = buffer.data(il + 411);
    const auto *il_412 = buffer.data(il + 412);
    const auto *il_413 = buffer.data(il + 413);
    const auto *il_414 = buffer.data(il + 414);
    const auto *il_415 = buffer.data(il + 415);
    const auto *il_416 = buffer.data(il + 416);
    const auto *il_417 = buffer.data(il + 417);
    const auto *il_418 = buffer.data(il + 418);
    const auto *il_419 = buffer.data(il + 419);
    const auto *il_420 = buffer.data(il + 420);
    const auto *il_421 = buffer.data(il + 421);
    const auto *il_422 = buffer.data(il + 422);
    const auto *il_423 = buffer.data(il + 423);
    const auto *il_424 = buffer.data(il + 424);
    const auto *il_425 = buffer.data(il + 425);
    const auto *il_426 = buffer.data(il + 426);
    const auto *il_427 = buffer.data(il + 427);
    const auto *il_428 = buffer.data(il + 428);
    const auto *il_429 = buffer.data(il + 429);
    const auto *il_430 = buffer.data(il + 430);
    const auto *il_431 = buffer.data(il + 431);
    const auto *il_432 = buffer.data(il + 432);
    const auto *il_433 = buffer.data(il + 433);
    const auto *il_434 = buffer.data(il + 434);
    const auto *il_435 = buffer.data(il + 435);
    const auto *il_436 = buffer.data(il + 436);
    const auto *il_437 = buffer.data(il + 437);
    const auto *il_438 = buffer.data(il + 438);
    const auto *il_439 = buffer.data(il + 439);
    const auto *il_440 = buffer.data(il + 440);
    const auto *il_441 = buffer.data(il + 441);
    const auto *il_442 = buffer.data(il + 442);
    const auto *il_443 = buffer.data(il + 443);
    const auto *il_444 = buffer.data(il + 444);
    const auto *il_445 = buffer.data(il + 445);
    const auto *il_446 = buffer.data(il + 446);
    const auto *il_447 = buffer.data(il + 447);
    const auto *il_448 = buffer.data(il + 448);
    const auto *il_449 = buffer.data(il + 449);
    const auto *il_495 = buffer.data(il + 495);
    const auto *il_496 = buffer.data(il + 496);
    const auto *il_497 = buffer.data(il + 497);
    const auto *il_498 = buffer.data(il + 498);
    const auto *il_499 = buffer.data(il + 499);
    const auto *il_500 = buffer.data(il + 500);
    const auto *il_501 = buffer.data(il + 501);
    const auto *il_502 = buffer.data(il + 502);
    const auto *il_503 = buffer.data(il + 503);
    const auto *il_504 = buffer.data(il + 504);
    const auto *il_505 = buffer.data(il + 505);
    const auto *il_506 = buffer.data(il + 506);
    const auto *il_507 = buffer.data(il + 507);
    const auto *il_508 = buffer.data(il + 508);
    const auto *il_509 = buffer.data(il + 509);
    const auto *il_510 = buffer.data(il + 510);
    const auto *il_511 = buffer.data(il + 511);
    const auto *il_512 = buffer.data(il + 512);
    const auto *il_513 = buffer.data(il + 513);
    const auto *il_514 = buffer.data(il + 514);
    const auto *il_515 = buffer.data(il + 515);
    const auto *il_516 = buffer.data(il + 516);
    const auto *il_517 = buffer.data(il + 517);
    const auto *il_518 = buffer.data(il + 518);
    const auto *il_519 = buffer.data(il + 519);
    const auto *il_520 = buffer.data(il + 520);
    const auto *il_521 = buffer.data(il + 521);
    const auto *il_522 = buffer.data(il + 522);
    const auto *il_523 = buffer.data(il + 523);
    const auto *il_524 = buffer.data(il + 524);
    const auto *il_525 = buffer.data(il + 525);
    const auto *il_526 = buffer.data(il + 526);
    const auto *il_527 = buffer.data(il + 527);
    const auto *il_528 = buffer.data(il + 528);
    const auto *il_529 = buffer.data(il + 529);
    const auto *il_530 = buffer.data(il + 530);
    const auto *il_531 = buffer.data(il + 531);
    const auto *il_532 = buffer.data(il + 532);
    const auto *il_533 = buffer.data(il + 533);
    const auto *il_534 = buffer.data(il + 534);
    const auto *il_535 = buffer.data(il + 535);
    const auto *il_536 = buffer.data(il + 536);
    const auto *il_537 = buffer.data(il + 537);
    const auto *il_538 = buffer.data(il + 538);
    const auto *il_539 = buffer.data(il + 539);
    const auto *il_540 = buffer.data(il + 540);
    const auto *il_541 = buffer.data(il + 541);
    const auto *il_542 = buffer.data(il + 542);
    const auto *il_543 = buffer.data(il + 543);
    const auto *il_544 = buffer.data(il + 544);
    const auto *il_545 = buffer.data(il + 545);
    const auto *il_546 = buffer.data(il + 546);
    const auto *il_547 = buffer.data(il + 547);
    const auto *il_548 = buffer.data(il + 548);
    const auto *il_549 = buffer.data(il + 549);
    const auto *il_550 = buffer.data(il + 550);
    const auto *il_551 = buffer.data(il + 551);
    const auto *il_552 = buffer.data(il + 552);
    const auto *il_553 = buffer.data(il + 553);
    const auto *il_554 = buffer.data(il + 554);
    const auto *il_555 = buffer.data(il + 555);
    const auto *il_556 = buffer.data(il + 556);
    const auto *il_557 = buffer.data(il + 557);
    const auto *il_558 = buffer.data(il + 558);
    const auto *il_559 = buffer.data(il + 559);
    const auto *il_560 = buffer.data(il + 560);
    const auto *il_561 = buffer.data(il + 561);
    const auto *il_562 = buffer.data(il + 562);
    const auto *il_563 = buffer.data(il + 563);
    const auto *il_564 = buffer.data(il + 564);
    const auto *il_565 = buffer.data(il + 565);
    const auto *il_566 = buffer.data(il + 566);
    const auto *il_567 = buffer.data(il + 567);
    const auto *il_568 = buffer.data(il + 568);
    const auto *il_569 = buffer.data(il + 569);
    const auto *il_570 = buffer.data(il + 570);
    const auto *il_571 = buffer.data(il + 571);
    const auto *il_572 = buffer.data(il + 572);
    const auto *il_573 = buffer.data(il + 573);
    const auto *il_574 = buffer.data(il + 574);
    const auto *il_575 = buffer.data(il + 575);
    const auto *il_576 = buffer.data(il + 576);
    const auto *il_577 = buffer.data(il + 577);
    const auto *il_578 = buffer.data(il + 578);
    const auto *il_579 = buffer.data(il + 579);
    const auto *il_580 = buffer.data(il + 580);
    const auto *il_581 = buffer.data(il + 581);

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, gl_60, gl_61, gl_62, gl_63, gl_64, \
                         il_375, il_376, il_377, il_378, il_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -gl_60[k]
                   + f_0 * il_375[k];

        t_196[k] = -gl_61[k]
                   + f_0 * il_376[k];

        t_197[k] = -gl_62[k]
                   + f_0 * il_377[k];

        t_198[k] = -gl_63[k]
                   + f_0 * il_378[k];

        t_199[k] = -gl_64[k]
                   + f_0 * il_379[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, gl_65, gl_66, gl_67, gl_68, gl_69, \
                         il_380, il_381, il_382, il_383, il_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -gl_65[k]
                   + f_0 * il_380[k];

        t_201[k] = -gl_66[k]
                   + f_0 * il_381[k];

        t_202[k] = -gl_67[k]
                   + f_0 * il_382[k];

        t_203[k] = -gl_68[k]
                   + f_0 * il_383[k];

        t_204[k] = -gl_69[k]
                   + f_0 * il_384[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, gl_70, gl_71, gl_72, gl_73, gl_74, \
                         il_385, il_386, il_387, il_388, il_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -gl_70[k]
                   + f_0 * il_385[k];

        t_206[k] = -gl_71[k]
                   + f_0 * il_386[k];

        t_207[k] = -gl_72[k]
                   + f_0 * il_387[k];

        t_208[k] = -gl_73[k]
                   + f_0 * il_388[k];

        t_209[k] = -gl_74[k]
                   + f_0 * il_389[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, gl_75, gl_76, gl_77, gl_78, gl_79, \
                         il_390, il_391, il_392, il_393, il_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -gl_75[k]
                   + f_0 * il_390[k];

        t_211[k] = -gl_76[k]
                   + f_0 * il_391[k];

        t_212[k] = -gl_77[k]
                   + f_0 * il_392[k];

        t_213[k] = -gl_78[k]
                   + f_0 * il_393[k];

        t_214[k] = -gl_79[k]
                   + f_0 * il_394[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, gl_80, gl_81, gl_82, gl_83, gl_84, \
                         il_395, il_396, il_397, il_398, il_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -gl_80[k]
                   + f_0 * il_395[k];

        t_216[k] = -gl_81[k]
                   + f_0 * il_396[k];

        t_217[k] = -gl_82[k]
                   + f_0 * il_397[k];

        t_218[k] = -gl_83[k]
                   + f_0 * il_398[k];

        t_219[k] = -gl_84[k]
                   + f_0 * il_399[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, gl_85, gl_86, gl_87, gl_88, gl_89, \
                         il_400, il_401, il_402, il_403, il_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -gl_85[k]
                   + f_0 * il_400[k];

        t_221[k] = -gl_86[k]
                   + f_0 * il_401[k];

        t_222[k] = -gl_87[k]
                   + f_0 * il_402[k];

        t_223[k] = -gl_88[k]
                   + f_0 * il_403[k];

        t_224[k] = -gl_89[k]
                   + f_0 * il_404[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, gl_90, gl_91, gl_92, gl_93, gl_94, \
                         il_405, il_406, il_407, il_408, il_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -2.0 * gl_90[k]
                   + f_0 * il_405[k];

        t_226[k] = -2.0 * gl_91[k]
                   + f_0 * il_406[k];

        t_227[k] = -2.0 * gl_92[k]
                   + f_0 * il_407[k];

        t_228[k] = -2.0 * gl_93[k]
                   + f_0 * il_408[k];

        t_229[k] = -2.0 * gl_94[k]
                   + f_0 * il_409[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, gl_95, gl_96, gl_97, gl_98, gl_99, \
                         il_410, il_411, il_412, il_413, il_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -2.0 * gl_95[k]
                   + f_0 * il_410[k];

        t_231[k] = -2.0 * gl_96[k]
                   + f_0 * il_411[k];

        t_232[k] = -2.0 * gl_97[k]
                   + f_0 * il_412[k];

        t_233[k] = -2.0 * gl_98[k]
                   + f_0 * il_413[k];

        t_234[k] = -2.0 * gl_99[k]
                   + f_0 * il_414[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, gl_100, gl_101, gl_102, gl_103, \
                         gl_104, il_415, il_416, il_417, il_418, \
                         il_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -2.0 * gl_100[k]
                   + f_0 * il_415[k];

        t_236[k] = -2.0 * gl_101[k]
                   + f_0 * il_416[k];

        t_237[k] = -2.0 * gl_102[k]
                   + f_0 * il_417[k];

        t_238[k] = -2.0 * gl_103[k]
                   + f_0 * il_418[k];

        t_239[k] = -2.0 * gl_104[k]
                   + f_0 * il_419[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, gl_105, gl_106, gl_107, gl_108, \
                         gl_109, il_420, il_421, il_422, il_423, \
                         il_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -2.0 * gl_105[k]
                   + f_0 * il_420[k];

        t_241[k] = -2.0 * gl_106[k]
                   + f_0 * il_421[k];

        t_242[k] = -2.0 * gl_107[k]
                   + f_0 * il_422[k];

        t_243[k] = -2.0 * gl_108[k]
                   + f_0 * il_423[k];

        t_244[k] = -2.0 * gl_109[k]
                   + f_0 * il_424[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, gl_110, gl_111, gl_112, gl_113, \
                         gl_114, il_425, il_426, il_427, il_428, \
                         il_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -2.0 * gl_110[k]
                   + f_0 * il_425[k];

        t_246[k] = -2.0 * gl_111[k]
                   + f_0 * il_426[k];

        t_247[k] = -2.0 * gl_112[k]
                   + f_0 * il_427[k];

        t_248[k] = -2.0 * gl_113[k]
                   + f_0 * il_428[k];

        t_249[k] = -2.0 * gl_114[k]
                   + f_0 * il_429[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, gl_115, gl_116, gl_117, gl_118, \
                         gl_119, il_430, il_431, il_432, il_433, \
                         il_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * gl_115[k]
                   + f_0 * il_430[k];

        t_251[k] = -2.0 * gl_116[k]
                   + f_0 * il_431[k];

        t_252[k] = -2.0 * gl_117[k]
                   + f_0 * il_432[k];

        t_253[k] = -2.0 * gl_118[k]
                   + f_0 * il_433[k];

        t_254[k] = -2.0 * gl_119[k]
                   + f_0 * il_434[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, gl_120, gl_121, gl_122, gl_123, \
                         gl_124, il_435, il_436, il_437, il_438, \
                         il_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * gl_120[k]
                   + f_0 * il_435[k];

        t_256[k] = -2.0 * gl_121[k]
                   + f_0 * il_436[k];

        t_257[k] = -2.0 * gl_122[k]
                   + f_0 * il_437[k];

        t_258[k] = -2.0 * gl_123[k]
                   + f_0 * il_438[k];

        t_259[k] = -2.0 * gl_124[k]
                   + f_0 * il_439[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, gl_125, gl_126, gl_127, gl_128, \
                         gl_129, il_440, il_441, il_442, il_443, \
                         il_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * gl_125[k]
                   + f_0 * il_440[k];

        t_261[k] = -2.0 * gl_126[k]
                   + f_0 * il_441[k];

        t_262[k] = -2.0 * gl_127[k]
                   + f_0 * il_442[k];

        t_263[k] = -2.0 * gl_128[k]
                   + f_0 * il_443[k];

        t_264[k] = -2.0 * gl_129[k]
                   + f_0 * il_444[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, gl_130, gl_131, gl_132, gl_133, \
                         gl_134, il_445, il_446, il_447, il_448, \
                         il_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * gl_130[k]
                   + f_0 * il_445[k];

        t_266[k] = -2.0 * gl_131[k]
                   + f_0 * il_446[k];

        t_267[k] = -2.0 * gl_132[k]
                   + f_0 * il_447[k];

        t_268[k] = -2.0 * gl_133[k]
                   + f_0 * il_448[k];

        t_269[k] = -2.0 * gl_134[k]
                   + f_0 * il_449[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, t_276, t_277, il_495, \
                         il_496, il_497, il_498, il_499, il_500, il_501, \
                         il_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_0 * il_495[k];

        t_271[k] = f_0 * il_496[k];

        t_272[k] = f_0 * il_497[k];

        t_273[k] = f_0 * il_498[k];

        t_274[k] = f_0 * il_499[k];

        t_275[k] = f_0 * il_500[k];

        t_276[k] = f_0 * il_501[k];

        t_277[k] = f_0 * il_502[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, t_283, t_284, t_285, il_503, \
                         il_504, il_505, il_506, il_507, il_508, il_509, \
                         il_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_0 * il_503[k];

        t_279[k] = f_0 * il_504[k];

        t_280[k] = f_0 * il_505[k];

        t_281[k] = f_0 * il_506[k];

        t_282[k] = f_0 * il_507[k];

        t_283[k] = f_0 * il_508[k];

        t_284[k] = f_0 * il_509[k];

        t_285[k] = f_0 * il_510[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, t_291, t_292, t_293, il_511, \
                         il_512, il_513, il_514, il_515, il_516, il_517, \
                         il_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * il_511[k];

        t_287[k] = f_0 * il_512[k];

        t_288[k] = f_0 * il_513[k];

        t_289[k] = f_0 * il_514[k];

        t_290[k] = f_0 * il_515[k];

        t_291[k] = f_0 * il_516[k];

        t_292[k] = f_0 * il_517[k];

        t_293[k] = f_0 * il_518[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, t_299, t_300, t_301, il_519, \
                         il_520, il_521, il_522, il_523, il_524, il_525, \
                         il_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_0 * il_519[k];

        t_295[k] = f_0 * il_520[k];

        t_296[k] = f_0 * il_521[k];

        t_297[k] = f_0 * il_522[k];

        t_298[k] = f_0 * il_523[k];

        t_299[k] = f_0 * il_524[k];

        t_300[k] = f_0 * il_525[k];

        t_301[k] = f_0 * il_526[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, t_307, t_308, t_309, il_527, \
                         il_528, il_529, il_530, il_531, il_532, il_533, \
                         il_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_0 * il_527[k];

        t_303[k] = f_0 * il_528[k];

        t_304[k] = f_0 * il_529[k];

        t_305[k] = f_0 * il_530[k];

        t_306[k] = f_0 * il_531[k];

        t_307[k] = f_0 * il_532[k];

        t_308[k] = f_0 * il_533[k];

        t_309[k] = f_0 * il_534[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, t_315, t_316, gl_135, gl_136, \
                         il_535, il_536, il_537, il_538, il_539, il_540, \
                         il_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_0 * il_535[k];

        t_311[k] = f_0 * il_536[k];

        t_312[k] = f_0 * il_537[k];

        t_313[k] = f_0 * il_538[k];

        t_314[k] = f_0 * il_539[k];

        t_315[k] = -gl_135[k]
                   + f_0 * il_540[k];

        t_316[k] = -gl_136[k]
                   + f_0 * il_541[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, gl_137, gl_138, gl_139, gl_140, \
                         gl_141, il_542, il_543, il_544, il_545, \
                         il_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = -gl_137[k]
                   + f_0 * il_542[k];

        t_318[k] = -gl_138[k]
                   + f_0 * il_543[k];

        t_319[k] = -gl_139[k]
                   + f_0 * il_544[k];

        t_320[k] = -gl_140[k]
                   + f_0 * il_545[k];

        t_321[k] = -gl_141[k]
                   + f_0 * il_546[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, gl_142, gl_143, gl_144, gl_145, \
                         gl_146, il_547, il_548, il_549, il_550, \
                         il_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = -gl_142[k]
                   + f_0 * il_547[k];

        t_323[k] = -gl_143[k]
                   + f_0 * il_548[k];

        t_324[k] = -gl_144[k]
                   + f_0 * il_549[k];

        t_325[k] = -gl_145[k]
                   + f_0 * il_550[k];

        t_326[k] = -gl_146[k]
                   + f_0 * il_551[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, gl_147, gl_148, gl_149, gl_150, \
                         gl_151, il_552, il_553, il_554, il_555, \
                         il_556 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = -gl_147[k]
                   + f_0 * il_552[k];

        t_328[k] = -gl_148[k]
                   + f_0 * il_553[k];

        t_329[k] = -gl_149[k]
                   + f_0 * il_554[k];

        t_330[k] = -gl_150[k]
                   + f_0 * il_555[k];

        t_331[k] = -gl_151[k]
                   + f_0 * il_556[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, gl_152, gl_153, gl_154, gl_155, \
                         gl_156, il_557, il_558, il_559, il_560, \
                         il_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -gl_152[k]
                   + f_0 * il_557[k];

        t_333[k] = -gl_153[k]
                   + f_0 * il_558[k];

        t_334[k] = -gl_154[k]
                   + f_0 * il_559[k];

        t_335[k] = -gl_155[k]
                   + f_0 * il_560[k];

        t_336[k] = -gl_156[k]
                   + f_0 * il_561[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, gl_157, gl_158, gl_159, gl_160, \
                         gl_161, il_562, il_563, il_564, il_565, \
                         il_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -gl_157[k]
                   + f_0 * il_562[k];

        t_338[k] = -gl_158[k]
                   + f_0 * il_563[k];

        t_339[k] = -gl_159[k]
                   + f_0 * il_564[k];

        t_340[k] = -gl_160[k]
                   + f_0 * il_565[k];

        t_341[k] = -gl_161[k]
                   + f_0 * il_566[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, gl_162, gl_163, gl_164, gl_165, \
                         gl_166, il_567, il_568, il_569, il_570, \
                         il_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -gl_162[k]
                   + f_0 * il_567[k];

        t_343[k] = -gl_163[k]
                   + f_0 * il_568[k];

        t_344[k] = -gl_164[k]
                   + f_0 * il_569[k];

        t_345[k] = -gl_165[k]
                   + f_0 * il_570[k];

        t_346[k] = -gl_166[k]
                   + f_0 * il_571[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, gl_167, gl_168, gl_169, gl_170, \
                         gl_171, il_572, il_573, il_574, il_575, \
                         il_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -gl_167[k]
                   + f_0 * il_572[k];

        t_348[k] = -gl_168[k]
                   + f_0 * il_573[k];

        t_349[k] = -gl_169[k]
                   + f_0 * il_574[k];

        t_350[k] = -gl_170[k]
                   + f_0 * il_575[k];

        t_351[k] = -gl_171[k]
                   + f_0 * il_576[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, gl_172, gl_173, gl_174, gl_175, \
                         gl_176, il_577, il_578, il_579, il_580, \
                         il_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -gl_172[k]
                   + f_0 * il_577[k];

        t_353[k] = -gl_173[k]
                   + f_0 * il_578[k];

        t_354[k] = -gl_174[k]
                   + f_0 * il_579[k];

        t_355[k] = -gl_175[k]
                   + f_0 * il_580[k];

        t_356[k] = -gl_176[k]
                   + f_0 * il_581[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_177 = buffer.data(gl + 177);
    const auto *gl_178 = buffer.data(gl + 178);
    const auto *gl_179 = buffer.data(gl + 179);
    const auto *gl_180 = buffer.data(gl + 180);
    const auto *gl_181 = buffer.data(gl + 181);
    const auto *gl_182 = buffer.data(gl + 182);
    const auto *gl_183 = buffer.data(gl + 183);
    const auto *gl_184 = buffer.data(gl + 184);
    const auto *gl_185 = buffer.data(gl + 185);
    const auto *gl_186 = buffer.data(gl + 186);
    const auto *gl_187 = buffer.data(gl + 187);
    const auto *gl_188 = buffer.data(gl + 188);
    const auto *gl_189 = buffer.data(gl + 189);
    const auto *gl_190 = buffer.data(gl + 190);
    const auto *gl_191 = buffer.data(gl + 191);
    const auto *gl_192 = buffer.data(gl + 192);
    const auto *gl_193 = buffer.data(gl + 193);
    const auto *gl_194 = buffer.data(gl + 194);
    const auto *gl_195 = buffer.data(gl + 195);
    const auto *gl_196 = buffer.data(gl + 196);
    const auto *gl_197 = buffer.data(gl + 197);
    const auto *gl_198 = buffer.data(gl + 198);
    const auto *gl_199 = buffer.data(gl + 199);
    const auto *gl_200 = buffer.data(gl + 200);
    const auto *gl_201 = buffer.data(gl + 201);
    const auto *gl_202 = buffer.data(gl + 202);
    const auto *gl_203 = buffer.data(gl + 203);
    const auto *gl_204 = buffer.data(gl + 204);
    const auto *gl_205 = buffer.data(gl + 205);
    const auto *gl_206 = buffer.data(gl + 206);
    const auto *gl_207 = buffer.data(gl + 207);
    const auto *gl_208 = buffer.data(gl + 208);
    const auto *gl_209 = buffer.data(gl + 209);
    const auto *gl_210 = buffer.data(gl + 210);
    const auto *gl_211 = buffer.data(gl + 211);
    const auto *gl_212 = buffer.data(gl + 212);
    const auto *gl_213 = buffer.data(gl + 213);
    const auto *gl_214 = buffer.data(gl + 214);
    const auto *gl_215 = buffer.data(gl + 215);
    const auto *gl_216 = buffer.data(gl + 216);
    const auto *gl_217 = buffer.data(gl + 217);
    const auto *gl_218 = buffer.data(gl + 218);
    const auto *gl_219 = buffer.data(gl + 219);
    const auto *gl_220 = buffer.data(gl + 220);
    const auto *gl_221 = buffer.data(gl + 221);
    const auto *gl_222 = buffer.data(gl + 222);
    const auto *gl_223 = buffer.data(gl + 223);
    const auto *gl_224 = buffer.data(gl + 224);
    const auto *gl_225 = buffer.data(gl + 225);
    const auto *gl_226 = buffer.data(gl + 226);
    const auto *gl_227 = buffer.data(gl + 227);
    const auto *gl_228 = buffer.data(gl + 228);
    const auto *gl_229 = buffer.data(gl + 229);
    const auto *gl_230 = buffer.data(gl + 230);
    const auto *gl_231 = buffer.data(gl + 231);
    const auto *gl_232 = buffer.data(gl + 232);
    const auto *gl_233 = buffer.data(gl + 233);
    const auto *gl_234 = buffer.data(gl + 234);
    const auto *gl_235 = buffer.data(gl + 235);
    const auto *gl_236 = buffer.data(gl + 236);
    const auto *gl_237 = buffer.data(gl + 237);
    const auto *gl_238 = buffer.data(gl + 238);
    const auto *gl_239 = buffer.data(gl + 239);
    const auto *gl_240 = buffer.data(gl + 240);
    const auto *gl_241 = buffer.data(gl + 241);
    const auto *gl_242 = buffer.data(gl + 242);
    const auto *gl_243 = buffer.data(gl + 243);
    const auto *gl_244 = buffer.data(gl + 244);
    const auto *gl_245 = buffer.data(gl + 245);
    const auto *gl_246 = buffer.data(gl + 246);
    const auto *gl_247 = buffer.data(gl + 247);
    const auto *gl_248 = buffer.data(gl + 248);
    const auto *gl_249 = buffer.data(gl + 249);
    const auto *gl_250 = buffer.data(gl + 250);
    const auto *gl_251 = buffer.data(gl + 251);
    const auto *gl_252 = buffer.data(gl + 252);
    const auto *gl_253 = buffer.data(gl + 253);
    const auto *gl_254 = buffer.data(gl + 254);
    const auto *gl_255 = buffer.data(gl + 255);
    const auto *gl_256 = buffer.data(gl + 256);
    const auto *gl_257 = buffer.data(gl + 257);
    const auto *gl_258 = buffer.data(gl + 258);
    const auto *gl_259 = buffer.data(gl + 259);
    const auto *gl_260 = buffer.data(gl + 260);
    const auto *gl_261 = buffer.data(gl + 261);
    const auto *gl_262 = buffer.data(gl + 262);
    const auto *gl_263 = buffer.data(gl + 263);
    const auto *gl_264 = buffer.data(gl + 264);
    const auto *gl_265 = buffer.data(gl + 265);
    const auto *gl_266 = buffer.data(gl + 266);
    const auto *gl_267 = buffer.data(gl + 267);
    const auto *gl_268 = buffer.data(gl + 268);
    const auto *gl_269 = buffer.data(gl + 269);
    const auto *gl_270 = buffer.data(gl + 270);
    const auto *gl_271 = buffer.data(gl + 271);
    const auto *gl_272 = buffer.data(gl + 272);
    const auto *gl_273 = buffer.data(gl + 273);
    const auto *gl_274 = buffer.data(gl + 274);
    const auto *gl_275 = buffer.data(gl + 275);
    const auto *gl_276 = buffer.data(gl + 276);
    const auto *gl_277 = buffer.data(gl + 277);
    const auto *gl_278 = buffer.data(gl + 278);
    const auto *gl_279 = buffer.data(gl + 279);
    const auto *gl_280 = buffer.data(gl + 280);
    const auto *gl_281 = buffer.data(gl + 281);
    const auto *gl_282 = buffer.data(gl + 282);
    const auto *gl_283 = buffer.data(gl + 283);
    const auto *gl_284 = buffer.data(gl + 284);
    const auto *gl_285 = buffer.data(gl + 285);
    const auto *gl_286 = buffer.data(gl + 286);
    const auto *gl_287 = buffer.data(gl + 287);
    const auto *gl_288 = buffer.data(gl + 288);
    const auto *gl_289 = buffer.data(gl + 289);
    const auto *gl_290 = buffer.data(gl + 290);
    const auto *gl_291 = buffer.data(gl + 291);
    const auto *gl_292 = buffer.data(gl + 292);
    const auto *gl_293 = buffer.data(gl + 293);

    const auto *il_582 = buffer.data(il + 582);
    const auto *il_583 = buffer.data(il + 583);
    const auto *il_584 = buffer.data(il + 584);
    const auto *il_585 = buffer.data(il + 585);
    const auto *il_586 = buffer.data(il + 586);
    const auto *il_587 = buffer.data(il + 587);
    const auto *il_588 = buffer.data(il + 588);
    const auto *il_589 = buffer.data(il + 589);
    const auto *il_590 = buffer.data(il + 590);
    const auto *il_591 = buffer.data(il + 591);
    const auto *il_592 = buffer.data(il + 592);
    const auto *il_593 = buffer.data(il + 593);
    const auto *il_594 = buffer.data(il + 594);
    const auto *il_595 = buffer.data(il + 595);
    const auto *il_596 = buffer.data(il + 596);
    const auto *il_597 = buffer.data(il + 597);
    const auto *il_598 = buffer.data(il + 598);
    const auto *il_599 = buffer.data(il + 599);
    const auto *il_600 = buffer.data(il + 600);
    const auto *il_601 = buffer.data(il + 601);
    const auto *il_602 = buffer.data(il + 602);
    const auto *il_603 = buffer.data(il + 603);
    const auto *il_604 = buffer.data(il + 604);
    const auto *il_605 = buffer.data(il + 605);
    const auto *il_606 = buffer.data(il + 606);
    const auto *il_607 = buffer.data(il + 607);
    const auto *il_608 = buffer.data(il + 608);
    const auto *il_609 = buffer.data(il + 609);
    const auto *il_610 = buffer.data(il + 610);
    const auto *il_611 = buffer.data(il + 611);
    const auto *il_612 = buffer.data(il + 612);
    const auto *il_613 = buffer.data(il + 613);
    const auto *il_614 = buffer.data(il + 614);
    const auto *il_615 = buffer.data(il + 615);
    const auto *il_616 = buffer.data(il + 616);
    const auto *il_617 = buffer.data(il + 617);
    const auto *il_618 = buffer.data(il + 618);
    const auto *il_619 = buffer.data(il + 619);
    const auto *il_620 = buffer.data(il + 620);
    const auto *il_621 = buffer.data(il + 621);
    const auto *il_622 = buffer.data(il + 622);
    const auto *il_623 = buffer.data(il + 623);
    const auto *il_624 = buffer.data(il + 624);
    const auto *il_625 = buffer.data(il + 625);
    const auto *il_626 = buffer.data(il + 626);
    const auto *il_627 = buffer.data(il + 627);
    const auto *il_628 = buffer.data(il + 628);
    const auto *il_629 = buffer.data(il + 629);
    const auto *il_630 = buffer.data(il + 630);
    const auto *il_631 = buffer.data(il + 631);
    const auto *il_632 = buffer.data(il + 632);
    const auto *il_633 = buffer.data(il + 633);
    const auto *il_634 = buffer.data(il + 634);
    const auto *il_635 = buffer.data(il + 635);
    const auto *il_636 = buffer.data(il + 636);
    const auto *il_637 = buffer.data(il + 637);
    const auto *il_638 = buffer.data(il + 638);
    const auto *il_639 = buffer.data(il + 639);
    const auto *il_640 = buffer.data(il + 640);
    const auto *il_641 = buffer.data(il + 641);
    const auto *il_642 = buffer.data(il + 642);
    const auto *il_643 = buffer.data(il + 643);
    const auto *il_644 = buffer.data(il + 644);
    const auto *il_645 = buffer.data(il + 645);
    const auto *il_646 = buffer.data(il + 646);
    const auto *il_647 = buffer.data(il + 647);
    const auto *il_648 = buffer.data(il + 648);
    const auto *il_649 = buffer.data(il + 649);
    const auto *il_650 = buffer.data(il + 650);
    const auto *il_651 = buffer.data(il + 651);
    const auto *il_652 = buffer.data(il + 652);
    const auto *il_653 = buffer.data(il + 653);
    const auto *il_654 = buffer.data(il + 654);
    const auto *il_655 = buffer.data(il + 655);
    const auto *il_656 = buffer.data(il + 656);
    const auto *il_657 = buffer.data(il + 657);
    const auto *il_658 = buffer.data(il + 658);
    const auto *il_659 = buffer.data(il + 659);
    const auto *il_660 = buffer.data(il + 660);
    const auto *il_661 = buffer.data(il + 661);
    const auto *il_662 = buffer.data(il + 662);
    const auto *il_663 = buffer.data(il + 663);
    const auto *il_664 = buffer.data(il + 664);
    const auto *il_665 = buffer.data(il + 665);
    const auto *il_666 = buffer.data(il + 666);
    const auto *il_667 = buffer.data(il + 667);
    const auto *il_668 = buffer.data(il + 668);
    const auto *il_669 = buffer.data(il + 669);
    const auto *il_670 = buffer.data(il + 670);
    const auto *il_671 = buffer.data(il + 671);
    const auto *il_672 = buffer.data(il + 672);
    const auto *il_673 = buffer.data(il + 673);
    const auto *il_674 = buffer.data(il + 674);
    const auto *il_720 = buffer.data(il + 720);
    const auto *il_721 = buffer.data(il + 721);
    const auto *il_722 = buffer.data(il + 722);
    const auto *il_723 = buffer.data(il + 723);
    const auto *il_724 = buffer.data(il + 724);
    const auto *il_725 = buffer.data(il + 725);
    const auto *il_726 = buffer.data(il + 726);
    const auto *il_727 = buffer.data(il + 727);
    const auto *il_728 = buffer.data(il + 728);
    const auto *il_729 = buffer.data(il + 729);
    const auto *il_730 = buffer.data(il + 730);
    const auto *il_731 = buffer.data(il + 731);
    const auto *il_732 = buffer.data(il + 732);
    const auto *il_733 = buffer.data(il + 733);
    const auto *il_734 = buffer.data(il + 734);
    const auto *il_735 = buffer.data(il + 735);
    const auto *il_736 = buffer.data(il + 736);
    const auto *il_737 = buffer.data(il + 737);
    const auto *il_738 = buffer.data(il + 738);
    const auto *il_739 = buffer.data(il + 739);
    const auto *il_740 = buffer.data(il + 740);
    const auto *il_741 = buffer.data(il + 741);
    const auto *il_742 = buffer.data(il + 742);
    const auto *il_743 = buffer.data(il + 743);
    const auto *il_744 = buffer.data(il + 744);
    const auto *il_745 = buffer.data(il + 745);
    const auto *il_746 = buffer.data(il + 746);
    const auto *il_747 = buffer.data(il + 747);
    const auto *il_748 = buffer.data(il + 748);
    const auto *il_749 = buffer.data(il + 749);
    const auto *il_750 = buffer.data(il + 750);
    const auto *il_751 = buffer.data(il + 751);
    const auto *il_752 = buffer.data(il + 752);
    const auto *il_753 = buffer.data(il + 753);
    const auto *il_754 = buffer.data(il + 754);
    const auto *il_755 = buffer.data(il + 755);
    const auto *il_756 = buffer.data(il + 756);
    const auto *il_757 = buffer.data(il + 757);
    const auto *il_758 = buffer.data(il + 758);
    const auto *il_759 = buffer.data(il + 759);
    const auto *il_760 = buffer.data(il + 760);
    const auto *il_761 = buffer.data(il + 761);
    const auto *il_762 = buffer.data(il + 762);
    const auto *il_763 = buffer.data(il + 763);
    const auto *il_764 = buffer.data(il + 764);
    const auto *il_765 = buffer.data(il + 765);
    const auto *il_766 = buffer.data(il + 766);
    const auto *il_767 = buffer.data(il + 767);
    const auto *il_768 = buffer.data(il + 768);
    const auto *il_769 = buffer.data(il + 769);
    const auto *il_770 = buffer.data(il + 770);
    const auto *il_771 = buffer.data(il + 771);
    const auto *il_772 = buffer.data(il + 772);
    const auto *il_773 = buffer.data(il + 773);
    const auto *il_774 = buffer.data(il + 774);
    const auto *il_775 = buffer.data(il + 775);
    const auto *il_776 = buffer.data(il + 776);
    const auto *il_777 = buffer.data(il + 777);
    const auto *il_778 = buffer.data(il + 778);
    const auto *il_779 = buffer.data(il + 779);
    const auto *il_780 = buffer.data(il + 780);
    const auto *il_781 = buffer.data(il + 781);
    const auto *il_782 = buffer.data(il + 782);
    const auto *il_783 = buffer.data(il + 783);
    const auto *il_784 = buffer.data(il + 784);
    const auto *il_785 = buffer.data(il + 785);
    const auto *il_786 = buffer.data(il + 786);
    const auto *il_787 = buffer.data(il + 787);
    const auto *il_788 = buffer.data(il + 788);

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, gl_177, gl_178, gl_179, gl_180, \
                         gl_181, il_582, il_583, il_584, il_585, \
                         il_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -gl_177[k]
                   + f_0 * il_582[k];

        t_358[k] = -gl_178[k]
                   + f_0 * il_583[k];

        t_359[k] = -gl_179[k]
                   + f_0 * il_584[k];

        t_360[k] = -2.0 * gl_180[k]
                   + f_0 * il_585[k];

        t_361[k] = -2.0 * gl_181[k]
                   + f_0 * il_586[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, gl_182, gl_183, gl_184, gl_185, \
                         gl_186, il_587, il_588, il_589, il_590, \
                         il_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -2.0 * gl_182[k]
                   + f_0 * il_587[k];

        t_363[k] = -2.0 * gl_183[k]
                   + f_0 * il_588[k];

        t_364[k] = -2.0 * gl_184[k]
                   + f_0 * il_589[k];

        t_365[k] = -2.0 * gl_185[k]
                   + f_0 * il_590[k];

        t_366[k] = -2.0 * gl_186[k]
                   + f_0 * il_591[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, gl_187, gl_188, gl_189, gl_190, \
                         gl_191, il_592, il_593, il_594, il_595, \
                         il_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -2.0 * gl_187[k]
                   + f_0 * il_592[k];

        t_368[k] = -2.0 * gl_188[k]
                   + f_0 * il_593[k];

        t_369[k] = -2.0 * gl_189[k]
                   + f_0 * il_594[k];

        t_370[k] = -2.0 * gl_190[k]
                   + f_0 * il_595[k];

        t_371[k] = -2.0 * gl_191[k]
                   + f_0 * il_596[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, gl_192, gl_193, gl_194, gl_195, \
                         gl_196, il_597, il_598, il_599, il_600, \
                         il_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -2.0 * gl_192[k]
                   + f_0 * il_597[k];

        t_373[k] = -2.0 * gl_193[k]
                   + f_0 * il_598[k];

        t_374[k] = -2.0 * gl_194[k]
                   + f_0 * il_599[k];

        t_375[k] = -2.0 * gl_195[k]
                   + f_0 * il_600[k];

        t_376[k] = -2.0 * gl_196[k]
                   + f_0 * il_601[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, gl_197, gl_198, gl_199, gl_200, \
                         gl_201, il_602, il_603, il_604, il_605, \
                         il_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -2.0 * gl_197[k]
                   + f_0 * il_602[k];

        t_378[k] = -2.0 * gl_198[k]
                   + f_0 * il_603[k];

        t_379[k] = -2.0 * gl_199[k]
                   + f_0 * il_604[k];

        t_380[k] = -2.0 * gl_200[k]
                   + f_0 * il_605[k];

        t_381[k] = -2.0 * gl_201[k]
                   + f_0 * il_606[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, gl_202, gl_203, gl_204, gl_205, \
                         gl_206, il_607, il_608, il_609, il_610, \
                         il_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -2.0 * gl_202[k]
                   + f_0 * il_607[k];

        t_383[k] = -2.0 * gl_203[k]
                   + f_0 * il_608[k];

        t_384[k] = -2.0 * gl_204[k]
                   + f_0 * il_609[k];

        t_385[k] = -2.0 * gl_205[k]
                   + f_0 * il_610[k];

        t_386[k] = -2.0 * gl_206[k]
                   + f_0 * il_611[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, gl_207, gl_208, gl_209, gl_210, \
                         gl_211, il_612, il_613, il_614, il_615, \
                         il_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -2.0 * gl_207[k]
                   + f_0 * il_612[k];

        t_388[k] = -2.0 * gl_208[k]
                   + f_0 * il_613[k];

        t_389[k] = -2.0 * gl_209[k]
                   + f_0 * il_614[k];

        t_390[k] = -2.0 * gl_210[k]
                   + f_0 * il_615[k];

        t_391[k] = -2.0 * gl_211[k]
                   + f_0 * il_616[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, gl_212, gl_213, gl_214, gl_215, \
                         gl_216, il_617, il_618, il_619, il_620, \
                         il_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -2.0 * gl_212[k]
                   + f_0 * il_617[k];

        t_393[k] = -2.0 * gl_213[k]
                   + f_0 * il_618[k];

        t_394[k] = -2.0 * gl_214[k]
                   + f_0 * il_619[k];

        t_395[k] = -2.0 * gl_215[k]
                   + f_0 * il_620[k];

        t_396[k] = -2.0 * gl_216[k]
                   + f_0 * il_621[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, gl_217, gl_218, gl_219, gl_220, \
                         gl_221, il_622, il_623, il_624, il_625, \
                         il_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -2.0 * gl_217[k]
                   + f_0 * il_622[k];

        t_398[k] = -2.0 * gl_218[k]
                   + f_0 * il_623[k];

        t_399[k] = -2.0 * gl_219[k]
                   + f_0 * il_624[k];

        t_400[k] = -2.0 * gl_220[k]
                   + f_0 * il_625[k];

        t_401[k] = -2.0 * gl_221[k]
                   + f_0 * il_626[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, gl_222, gl_223, gl_224, gl_225, \
                         gl_226, il_627, il_628, il_629, il_630, \
                         il_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -2.0 * gl_222[k]
                   + f_0 * il_627[k];

        t_403[k] = -2.0 * gl_223[k]
                   + f_0 * il_628[k];

        t_404[k] = -2.0 * gl_224[k]
                   + f_0 * il_629[k];

        t_405[k] = -3.0 * gl_225[k]
                   + f_0 * il_630[k];

        t_406[k] = -3.0 * gl_226[k]
                   + f_0 * il_631[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, gl_227, gl_228, gl_229, gl_230, \
                         gl_231, il_632, il_633, il_634, il_635, \
                         il_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -3.0 * gl_227[k]
                   + f_0 * il_632[k];

        t_408[k] = -3.0 * gl_228[k]
                   + f_0 * il_633[k];

        t_409[k] = -3.0 * gl_229[k]
                   + f_0 * il_634[k];

        t_410[k] = -3.0 * gl_230[k]
                   + f_0 * il_635[k];

        t_411[k] = -3.0 * gl_231[k]
                   + f_0 * il_636[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, gl_232, gl_233, gl_234, gl_235, \
                         gl_236, il_637, il_638, il_639, il_640, \
                         il_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -3.0 * gl_232[k]
                   + f_0 * il_637[k];

        t_413[k] = -3.0 * gl_233[k]
                   + f_0 * il_638[k];

        t_414[k] = -3.0 * gl_234[k]
                   + f_0 * il_639[k];

        t_415[k] = -3.0 * gl_235[k]
                   + f_0 * il_640[k];

        t_416[k] = -3.0 * gl_236[k]
                   + f_0 * il_641[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, gl_237, gl_238, gl_239, gl_240, \
                         gl_241, il_642, il_643, il_644, il_645, \
                         il_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -3.0 * gl_237[k]
                   + f_0 * il_642[k];

        t_418[k] = -3.0 * gl_238[k]
                   + f_0 * il_643[k];

        t_419[k] = -3.0 * gl_239[k]
                   + f_0 * il_644[k];

        t_420[k] = -3.0 * gl_240[k]
                   + f_0 * il_645[k];

        t_421[k] = -3.0 * gl_241[k]
                   + f_0 * il_646[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, t_426, gl_242, gl_243, gl_244, gl_245, \
                         gl_246, il_647, il_648, il_649, il_650, \
                         il_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = -3.0 * gl_242[k]
                   + f_0 * il_647[k];

        t_423[k] = -3.0 * gl_243[k]
                   + f_0 * il_648[k];

        t_424[k] = -3.0 * gl_244[k]
                   + f_0 * il_649[k];

        t_425[k] = -3.0 * gl_245[k]
                   + f_0 * il_650[k];

        t_426[k] = -3.0 * gl_246[k]
                   + f_0 * il_651[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, t_431, gl_247, gl_248, gl_249, gl_250, \
                         gl_251, il_652, il_653, il_654, il_655, \
                         il_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = -3.0 * gl_247[k]
                   + f_0 * il_652[k];

        t_428[k] = -3.0 * gl_248[k]
                   + f_0 * il_653[k];

        t_429[k] = -3.0 * gl_249[k]
                   + f_0 * il_654[k];

        t_430[k] = -3.0 * gl_250[k]
                   + f_0 * il_655[k];

        t_431[k] = -3.0 * gl_251[k]
                   + f_0 * il_656[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, t_436, gl_252, gl_253, gl_254, gl_255, \
                         gl_256, il_657, il_658, il_659, il_660, \
                         il_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = -3.0 * gl_252[k]
                   + f_0 * il_657[k];

        t_433[k] = -3.0 * gl_253[k]
                   + f_0 * il_658[k];

        t_434[k] = -3.0 * gl_254[k]
                   + f_0 * il_659[k];

        t_435[k] = -3.0 * gl_255[k]
                   + f_0 * il_660[k];

        t_436[k] = -3.0 * gl_256[k]
                   + f_0 * il_661[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, t_441, gl_257, gl_258, gl_259, gl_260, \
                         gl_261, il_662, il_663, il_664, il_665, \
                         il_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = -3.0 * gl_257[k]
                   + f_0 * il_662[k];

        t_438[k] = -3.0 * gl_258[k]
                   + f_0 * il_663[k];

        t_439[k] = -3.0 * gl_259[k]
                   + f_0 * il_664[k];

        t_440[k] = -3.0 * gl_260[k]
                   + f_0 * il_665[k];

        t_441[k] = -3.0 * gl_261[k]
                   + f_0 * il_666[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, t_446, gl_262, gl_263, gl_264, gl_265, \
                         gl_266, il_667, il_668, il_669, il_670, \
                         il_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = -3.0 * gl_262[k]
                   + f_0 * il_667[k];

        t_443[k] = -3.0 * gl_263[k]
                   + f_0 * il_668[k];

        t_444[k] = -3.0 * gl_264[k]
                   + f_0 * il_669[k];

        t_445[k] = -3.0 * gl_265[k]
                   + f_0 * il_670[k];

        t_446[k] = -3.0 * gl_266[k]
                   + f_0 * il_671[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, t_451, t_452, gl_267, gl_268, gl_269, \
                         il_672, il_673, il_674, il_720, il_721, \
                         il_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = -3.0 * gl_267[k]
                   + f_0 * il_672[k];

        t_448[k] = -3.0 * gl_268[k]
                   + f_0 * il_673[k];

        t_449[k] = -3.0 * gl_269[k]
                   + f_0 * il_674[k];

        t_450[k] = f_0 * il_720[k];

        t_451[k] = f_0 * il_721[k];

        t_452[k] = f_0 * il_722[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, t_457, t_458, t_459, t_460, il_723, \
                         il_724, il_725, il_726, il_727, il_728, il_729, \
                         il_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_0 * il_723[k];

        t_454[k] = f_0 * il_724[k];

        t_455[k] = f_0 * il_725[k];

        t_456[k] = f_0 * il_726[k];

        t_457[k] = f_0 * il_727[k];

        t_458[k] = f_0 * il_728[k];

        t_459[k] = f_0 * il_729[k];

        t_460[k] = f_0 * il_730[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, t_465, t_466, t_467, t_468, il_731, \
                         il_732, il_733, il_734, il_735, il_736, il_737, \
                         il_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_0 * il_731[k];

        t_462[k] = f_0 * il_732[k];

        t_463[k] = f_0 * il_733[k];

        t_464[k] = f_0 * il_734[k];

        t_465[k] = f_0 * il_735[k];

        t_466[k] = f_0 * il_736[k];

        t_467[k] = f_0 * il_737[k];

        t_468[k] = f_0 * il_738[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, t_473, t_474, t_475, t_476, il_739, \
                         il_740, il_741, il_742, il_743, il_744, il_745, \
                         il_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_0 * il_739[k];

        t_470[k] = f_0 * il_740[k];

        t_471[k] = f_0 * il_741[k];

        t_472[k] = f_0 * il_742[k];

        t_473[k] = f_0 * il_743[k];

        t_474[k] = f_0 * il_744[k];

        t_475[k] = f_0 * il_745[k];

        t_476[k] = f_0 * il_746[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, t_482, t_483, t_484, il_747, \
                         il_748, il_749, il_750, il_751, il_752, il_753, \
                         il_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_0 * il_747[k];

        t_478[k] = f_0 * il_748[k];

        t_479[k] = f_0 * il_749[k];

        t_480[k] = f_0 * il_750[k];

        t_481[k] = f_0 * il_751[k];

        t_482[k] = f_0 * il_752[k];

        t_483[k] = f_0 * il_753[k];

        t_484[k] = f_0 * il_754[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, t_490, t_491, t_492, il_755, \
                         il_756, il_757, il_758, il_759, il_760, il_761, \
                         il_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_0 * il_755[k];

        t_486[k] = f_0 * il_756[k];

        t_487[k] = f_0 * il_757[k];

        t_488[k] = f_0 * il_758[k];

        t_489[k] = f_0 * il_759[k];

        t_490[k] = f_0 * il_760[k];

        t_491[k] = f_0 * il_761[k];

        t_492[k] = f_0 * il_762[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, t_496, t_497, t_498, gl_270, gl_271, gl_272, \
                         gl_273, il_763, il_764, il_765, il_766, il_767, \
                         il_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = f_0 * il_763[k];

        t_494[k] = f_0 * il_764[k];

        t_495[k] = -gl_270[k]
                   + f_0 * il_765[k];

        t_496[k] = -gl_271[k]
                   + f_0 * il_766[k];

        t_497[k] = -gl_272[k]
                   + f_0 * il_767[k];

        t_498[k] = -gl_273[k]
                   + f_0 * il_768[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, t_503, gl_274, gl_275, gl_276, gl_277, \
                         gl_278, il_769, il_770, il_771, il_772, \
                         il_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = -gl_274[k]
                   + f_0 * il_769[k];

        t_500[k] = -gl_275[k]
                   + f_0 * il_770[k];

        t_501[k] = -gl_276[k]
                   + f_0 * il_771[k];

        t_502[k] = -gl_277[k]
                   + f_0 * il_772[k];

        t_503[k] = -gl_278[k]
                   + f_0 * il_773[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, gl_279, gl_280, gl_281, gl_282, \
                         gl_283, il_774, il_775, il_776, il_777, \
                         il_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = -gl_279[k]
                   + f_0 * il_774[k];

        t_505[k] = -gl_280[k]
                   + f_0 * il_775[k];

        t_506[k] = -gl_281[k]
                   + f_0 * il_776[k];

        t_507[k] = -gl_282[k]
                   + f_0 * il_777[k];

        t_508[k] = -gl_283[k]
                   + f_0 * il_778[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, t_513, gl_284, gl_285, gl_286, gl_287, \
                         gl_288, il_779, il_780, il_781, il_782, \
                         il_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = -gl_284[k]
                   + f_0 * il_779[k];

        t_510[k] = -gl_285[k]
                   + f_0 * il_780[k];

        t_511[k] = -gl_286[k]
                   + f_0 * il_781[k];

        t_512[k] = -gl_287[k]
                   + f_0 * il_782[k];

        t_513[k] = -gl_288[k]
                   + f_0 * il_783[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, gl_289, gl_290, gl_291, gl_292, \
                         gl_293, il_784, il_785, il_786, il_787, \
                         il_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = -gl_289[k]
                   + f_0 * il_784[k];

        t_515[k] = -gl_290[k]
                   + f_0 * il_785[k];

        t_516[k] = -gl_291[k]
                   + f_0 * il_786[k];

        t_517[k] = -gl_292[k]
                   + f_0 * il_787[k];

        t_518[k] = -gl_293[k]
                   + f_0 * il_788[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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

    const auto *gl_294 = buffer.data(gl + 294);
    const auto *gl_295 = buffer.data(gl + 295);
    const auto *gl_296 = buffer.data(gl + 296);
    const auto *gl_297 = buffer.data(gl + 297);
    const auto *gl_298 = buffer.data(gl + 298);
    const auto *gl_299 = buffer.data(gl + 299);
    const auto *gl_300 = buffer.data(gl + 300);
    const auto *gl_301 = buffer.data(gl + 301);
    const auto *gl_302 = buffer.data(gl + 302);
    const auto *gl_303 = buffer.data(gl + 303);
    const auto *gl_304 = buffer.data(gl + 304);
    const auto *gl_305 = buffer.data(gl + 305);
    const auto *gl_306 = buffer.data(gl + 306);
    const auto *gl_307 = buffer.data(gl + 307);
    const auto *gl_308 = buffer.data(gl + 308);
    const auto *gl_309 = buffer.data(gl + 309);
    const auto *gl_310 = buffer.data(gl + 310);
    const auto *gl_311 = buffer.data(gl + 311);
    const auto *gl_312 = buffer.data(gl + 312);
    const auto *gl_313 = buffer.data(gl + 313);
    const auto *gl_314 = buffer.data(gl + 314);
    const auto *gl_315 = buffer.data(gl + 315);
    const auto *gl_316 = buffer.data(gl + 316);
    const auto *gl_317 = buffer.data(gl + 317);
    const auto *gl_318 = buffer.data(gl + 318);
    const auto *gl_319 = buffer.data(gl + 319);
    const auto *gl_320 = buffer.data(gl + 320);
    const auto *gl_321 = buffer.data(gl + 321);
    const auto *gl_322 = buffer.data(gl + 322);
    const auto *gl_323 = buffer.data(gl + 323);
    const auto *gl_324 = buffer.data(gl + 324);
    const auto *gl_325 = buffer.data(gl + 325);
    const auto *gl_326 = buffer.data(gl + 326);
    const auto *gl_327 = buffer.data(gl + 327);
    const auto *gl_328 = buffer.data(gl + 328);
    const auto *gl_329 = buffer.data(gl + 329);
    const auto *gl_330 = buffer.data(gl + 330);
    const auto *gl_331 = buffer.data(gl + 331);
    const auto *gl_332 = buffer.data(gl + 332);
    const auto *gl_333 = buffer.data(gl + 333);
    const auto *gl_334 = buffer.data(gl + 334);
    const auto *gl_335 = buffer.data(gl + 335);
    const auto *gl_336 = buffer.data(gl + 336);
    const auto *gl_337 = buffer.data(gl + 337);
    const auto *gl_338 = buffer.data(gl + 338);
    const auto *gl_339 = buffer.data(gl + 339);
    const auto *gl_340 = buffer.data(gl + 340);
    const auto *gl_341 = buffer.data(gl + 341);
    const auto *gl_342 = buffer.data(gl + 342);
    const auto *gl_343 = buffer.data(gl + 343);
    const auto *gl_344 = buffer.data(gl + 344);
    const auto *gl_345 = buffer.data(gl + 345);
    const auto *gl_346 = buffer.data(gl + 346);
    const auto *gl_347 = buffer.data(gl + 347);
    const auto *gl_348 = buffer.data(gl + 348);
    const auto *gl_349 = buffer.data(gl + 349);
    const auto *gl_350 = buffer.data(gl + 350);
    const auto *gl_351 = buffer.data(gl + 351);
    const auto *gl_352 = buffer.data(gl + 352);
    const auto *gl_353 = buffer.data(gl + 353);
    const auto *gl_354 = buffer.data(gl + 354);
    const auto *gl_355 = buffer.data(gl + 355);
    const auto *gl_356 = buffer.data(gl + 356);
    const auto *gl_357 = buffer.data(gl + 357);
    const auto *gl_358 = buffer.data(gl + 358);
    const auto *gl_359 = buffer.data(gl + 359);
    const auto *gl_360 = buffer.data(gl + 360);
    const auto *gl_361 = buffer.data(gl + 361);
    const auto *gl_362 = buffer.data(gl + 362);
    const auto *gl_363 = buffer.data(gl + 363);
    const auto *gl_364 = buffer.data(gl + 364);
    const auto *gl_365 = buffer.data(gl + 365);
    const auto *gl_366 = buffer.data(gl + 366);
    const auto *gl_367 = buffer.data(gl + 367);
    const auto *gl_368 = buffer.data(gl + 368);
    const auto *gl_369 = buffer.data(gl + 369);
    const auto *gl_370 = buffer.data(gl + 370);
    const auto *gl_371 = buffer.data(gl + 371);
    const auto *gl_372 = buffer.data(gl + 372);
    const auto *gl_373 = buffer.data(gl + 373);
    const auto *gl_374 = buffer.data(gl + 374);
    const auto *gl_375 = buffer.data(gl + 375);
    const auto *gl_376 = buffer.data(gl + 376);
    const auto *gl_377 = buffer.data(gl + 377);
    const auto *gl_378 = buffer.data(gl + 378);
    const auto *gl_379 = buffer.data(gl + 379);
    const auto *gl_380 = buffer.data(gl + 380);
    const auto *gl_381 = buffer.data(gl + 381);
    const auto *gl_382 = buffer.data(gl + 382);
    const auto *gl_383 = buffer.data(gl + 383);
    const auto *gl_384 = buffer.data(gl + 384);
    const auto *gl_385 = buffer.data(gl + 385);
    const auto *gl_386 = buffer.data(gl + 386);
    const auto *gl_387 = buffer.data(gl + 387);
    const auto *gl_388 = buffer.data(gl + 388);
    const auto *gl_389 = buffer.data(gl + 389);
    const auto *gl_390 = buffer.data(gl + 390);
    const auto *gl_391 = buffer.data(gl + 391);
    const auto *gl_392 = buffer.data(gl + 392);
    const auto *gl_393 = buffer.data(gl + 393);
    const auto *gl_394 = buffer.data(gl + 394);
    const auto *gl_395 = buffer.data(gl + 395);
    const auto *gl_396 = buffer.data(gl + 396);
    const auto *gl_397 = buffer.data(gl + 397);
    const auto *gl_398 = buffer.data(gl + 398);
    const auto *gl_399 = buffer.data(gl + 399);
    const auto *gl_400 = buffer.data(gl + 400);
    const auto *gl_401 = buffer.data(gl + 401);
    const auto *gl_402 = buffer.data(gl + 402);
    const auto *gl_403 = buffer.data(gl + 403);
    const auto *gl_404 = buffer.data(gl + 404);
    const auto *gl_405 = buffer.data(gl + 405);
    const auto *gl_406 = buffer.data(gl + 406);
    const auto *gl_407 = buffer.data(gl + 407);
    const auto *gl_408 = buffer.data(gl + 408);
    const auto *gl_409 = buffer.data(gl + 409);
    const auto *gl_410 = buffer.data(gl + 410);
    const auto *gl_411 = buffer.data(gl + 411);
    const auto *gl_412 = buffer.data(gl + 412);
    const auto *gl_413 = buffer.data(gl + 413);
    const auto *gl_414 = buffer.data(gl + 414);
    const auto *gl_415 = buffer.data(gl + 415);
    const auto *gl_416 = buffer.data(gl + 416);
    const auto *gl_417 = buffer.data(gl + 417);
    const auto *gl_418 = buffer.data(gl + 418);
    const auto *gl_419 = buffer.data(gl + 419);
    const auto *gl_420 = buffer.data(gl + 420);
    const auto *gl_421 = buffer.data(gl + 421);
    const auto *gl_422 = buffer.data(gl + 422);
    const auto *gl_423 = buffer.data(gl + 423);
    const auto *gl_424 = buffer.data(gl + 424);
    const auto *gl_425 = buffer.data(gl + 425);
    const auto *gl_426 = buffer.data(gl + 426);
    const auto *gl_427 = buffer.data(gl + 427);
    const auto *gl_428 = buffer.data(gl + 428);
    const auto *gl_429 = buffer.data(gl + 429);
    const auto *gl_430 = buffer.data(gl + 430);
    const auto *gl_431 = buffer.data(gl + 431);
    const auto *gl_432 = buffer.data(gl + 432);
    const auto *gl_433 = buffer.data(gl + 433);
    const auto *gl_434 = buffer.data(gl + 434);
    const auto *gl_435 = buffer.data(gl + 435);
    const auto *gl_436 = buffer.data(gl + 436);
    const auto *gl_437 = buffer.data(gl + 437);
    const auto *gl_438 = buffer.data(gl + 438);
    const auto *gl_439 = buffer.data(gl + 439);
    const auto *gl_440 = buffer.data(gl + 440);
    const auto *gl_441 = buffer.data(gl + 441);
    const auto *gl_442 = buffer.data(gl + 442);
    const auto *gl_443 = buffer.data(gl + 443);

    const auto *il_789 = buffer.data(il + 789);
    const auto *il_790 = buffer.data(il + 790);
    const auto *il_791 = buffer.data(il + 791);
    const auto *il_792 = buffer.data(il + 792);
    const auto *il_793 = buffer.data(il + 793);
    const auto *il_794 = buffer.data(il + 794);
    const auto *il_795 = buffer.data(il + 795);
    const auto *il_796 = buffer.data(il + 796);
    const auto *il_797 = buffer.data(il + 797);
    const auto *il_798 = buffer.data(il + 798);
    const auto *il_799 = buffer.data(il + 799);
    const auto *il_800 = buffer.data(il + 800);
    const auto *il_801 = buffer.data(il + 801);
    const auto *il_802 = buffer.data(il + 802);
    const auto *il_803 = buffer.data(il + 803);
    const auto *il_804 = buffer.data(il + 804);
    const auto *il_805 = buffer.data(il + 805);
    const auto *il_806 = buffer.data(il + 806);
    const auto *il_807 = buffer.data(il + 807);
    const auto *il_808 = buffer.data(il + 808);
    const auto *il_809 = buffer.data(il + 809);
    const auto *il_810 = buffer.data(il + 810);
    const auto *il_811 = buffer.data(il + 811);
    const auto *il_812 = buffer.data(il + 812);
    const auto *il_813 = buffer.data(il + 813);
    const auto *il_814 = buffer.data(il + 814);
    const auto *il_815 = buffer.data(il + 815);
    const auto *il_816 = buffer.data(il + 816);
    const auto *il_817 = buffer.data(il + 817);
    const auto *il_818 = buffer.data(il + 818);
    const auto *il_819 = buffer.data(il + 819);
    const auto *il_820 = buffer.data(il + 820);
    const auto *il_821 = buffer.data(il + 821);
    const auto *il_822 = buffer.data(il + 822);
    const auto *il_823 = buffer.data(il + 823);
    const auto *il_824 = buffer.data(il + 824);
    const auto *il_825 = buffer.data(il + 825);
    const auto *il_826 = buffer.data(il + 826);
    const auto *il_827 = buffer.data(il + 827);
    const auto *il_828 = buffer.data(il + 828);
    const auto *il_829 = buffer.data(il + 829);
    const auto *il_830 = buffer.data(il + 830);
    const auto *il_831 = buffer.data(il + 831);
    const auto *il_832 = buffer.data(il + 832);
    const auto *il_833 = buffer.data(il + 833);
    const auto *il_834 = buffer.data(il + 834);
    const auto *il_835 = buffer.data(il + 835);
    const auto *il_836 = buffer.data(il + 836);
    const auto *il_837 = buffer.data(il + 837);
    const auto *il_838 = buffer.data(il + 838);
    const auto *il_839 = buffer.data(il + 839);
    const auto *il_840 = buffer.data(il + 840);
    const auto *il_841 = buffer.data(il + 841);
    const auto *il_842 = buffer.data(il + 842);
    const auto *il_843 = buffer.data(il + 843);
    const auto *il_844 = buffer.data(il + 844);
    const auto *il_845 = buffer.data(il + 845);
    const auto *il_846 = buffer.data(il + 846);
    const auto *il_847 = buffer.data(il + 847);
    const auto *il_848 = buffer.data(il + 848);
    const auto *il_849 = buffer.data(il + 849);
    const auto *il_850 = buffer.data(il + 850);
    const auto *il_851 = buffer.data(il + 851);
    const auto *il_852 = buffer.data(il + 852);
    const auto *il_853 = buffer.data(il + 853);
    const auto *il_854 = buffer.data(il + 854);
    const auto *il_855 = buffer.data(il + 855);
    const auto *il_856 = buffer.data(il + 856);
    const auto *il_857 = buffer.data(il + 857);
    const auto *il_858 = buffer.data(il + 858);
    const auto *il_859 = buffer.data(il + 859);
    const auto *il_860 = buffer.data(il + 860);
    const auto *il_861 = buffer.data(il + 861);
    const auto *il_862 = buffer.data(il + 862);
    const auto *il_863 = buffer.data(il + 863);
    const auto *il_864 = buffer.data(il + 864);
    const auto *il_865 = buffer.data(il + 865);
    const auto *il_866 = buffer.data(il + 866);
    const auto *il_867 = buffer.data(il + 867);
    const auto *il_868 = buffer.data(il + 868);
    const auto *il_869 = buffer.data(il + 869);
    const auto *il_870 = buffer.data(il + 870);
    const auto *il_871 = buffer.data(il + 871);
    const auto *il_872 = buffer.data(il + 872);
    const auto *il_873 = buffer.data(il + 873);
    const auto *il_874 = buffer.data(il + 874);
    const auto *il_875 = buffer.data(il + 875);
    const auto *il_876 = buffer.data(il + 876);
    const auto *il_877 = buffer.data(il + 877);
    const auto *il_878 = buffer.data(il + 878);
    const auto *il_879 = buffer.data(il + 879);
    const auto *il_880 = buffer.data(il + 880);
    const auto *il_881 = buffer.data(il + 881);
    const auto *il_882 = buffer.data(il + 882);
    const auto *il_883 = buffer.data(il + 883);
    const auto *il_884 = buffer.data(il + 884);
    const auto *il_885 = buffer.data(il + 885);
    const auto *il_886 = buffer.data(il + 886);
    const auto *il_887 = buffer.data(il + 887);
    const auto *il_888 = buffer.data(il + 888);
    const auto *il_889 = buffer.data(il + 889);
    const auto *il_890 = buffer.data(il + 890);
    const auto *il_891 = buffer.data(il + 891);
    const auto *il_892 = buffer.data(il + 892);
    const auto *il_893 = buffer.data(il + 893);
    const auto *il_894 = buffer.data(il + 894);
    const auto *il_895 = buffer.data(il + 895);
    const auto *il_896 = buffer.data(il + 896);
    const auto *il_897 = buffer.data(il + 897);
    const auto *il_898 = buffer.data(il + 898);
    const auto *il_899 = buffer.data(il + 899);
    const auto *il_900 = buffer.data(il + 900);
    const auto *il_901 = buffer.data(il + 901);
    const auto *il_902 = buffer.data(il + 902);
    const auto *il_903 = buffer.data(il + 903);
    const auto *il_904 = buffer.data(il + 904);
    const auto *il_905 = buffer.data(il + 905);
    const auto *il_906 = buffer.data(il + 906);
    const auto *il_907 = buffer.data(il + 907);
    const auto *il_908 = buffer.data(il + 908);
    const auto *il_909 = buffer.data(il + 909);
    const auto *il_910 = buffer.data(il + 910);
    const auto *il_911 = buffer.data(il + 911);
    const auto *il_912 = buffer.data(il + 912);
    const auto *il_913 = buffer.data(il + 913);
    const auto *il_914 = buffer.data(il + 914);
    const auto *il_915 = buffer.data(il + 915);
    const auto *il_916 = buffer.data(il + 916);
    const auto *il_917 = buffer.data(il + 917);
    const auto *il_918 = buffer.data(il + 918);
    const auto *il_919 = buffer.data(il + 919);
    const auto *il_920 = buffer.data(il + 920);
    const auto *il_921 = buffer.data(il + 921);
    const auto *il_922 = buffer.data(il + 922);
    const auto *il_923 = buffer.data(il + 923);
    const auto *il_924 = buffer.data(il + 924);
    const auto *il_925 = buffer.data(il + 925);
    const auto *il_926 = buffer.data(il + 926);
    const auto *il_927 = buffer.data(il + 927);
    const auto *il_928 = buffer.data(il + 928);
    const auto *il_929 = buffer.data(il + 929);
    const auto *il_930 = buffer.data(il + 930);
    const auto *il_931 = buffer.data(il + 931);
    const auto *il_932 = buffer.data(il + 932);
    const auto *il_933 = buffer.data(il + 933);
    const auto *il_934 = buffer.data(il + 934);
    const auto *il_935 = buffer.data(il + 935);
    const auto *il_936 = buffer.data(il + 936);
    const auto *il_937 = buffer.data(il + 937);
    const auto *il_938 = buffer.data(il + 938);

#pragma omp simd aligned(t_519, t_520, t_521, t_522, t_523, gl_294, gl_295, gl_296, gl_297, \
                         gl_298, il_789, il_790, il_791, il_792, \
                         il_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = -gl_294[k]
                   + f_0 * il_789[k];

        t_520[k] = -gl_295[k]
                   + f_0 * il_790[k];

        t_521[k] = -gl_296[k]
                   + f_0 * il_791[k];

        t_522[k] = -gl_297[k]
                   + f_0 * il_792[k];

        t_523[k] = -gl_298[k]
                   + f_0 * il_793[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, gl_299, gl_300, gl_301, gl_302, \
                         gl_303, il_794, il_795, il_796, il_797, \
                         il_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = -gl_299[k]
                   + f_0 * il_794[k];

        t_525[k] = -gl_300[k]
                   + f_0 * il_795[k];

        t_526[k] = -gl_301[k]
                   + f_0 * il_796[k];

        t_527[k] = -gl_302[k]
                   + f_0 * il_797[k];

        t_528[k] = -gl_303[k]
                   + f_0 * il_798[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, t_533, gl_304, gl_305, gl_306, gl_307, \
                         gl_308, il_799, il_800, il_801, il_802, \
                         il_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = -gl_304[k]
                   + f_0 * il_799[k];

        t_530[k] = -gl_305[k]
                   + f_0 * il_800[k];

        t_531[k] = -gl_306[k]
                   + f_0 * il_801[k];

        t_532[k] = -gl_307[k]
                   + f_0 * il_802[k];

        t_533[k] = -gl_308[k]
                   + f_0 * il_803[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, t_538, gl_309, gl_310, gl_311, gl_312, \
                         gl_313, il_804, il_805, il_806, il_807, \
                         il_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = -gl_309[k]
                   + f_0 * il_804[k];

        t_535[k] = -gl_310[k]
                   + f_0 * il_805[k];

        t_536[k] = -gl_311[k]
                   + f_0 * il_806[k];

        t_537[k] = -gl_312[k]
                   + f_0 * il_807[k];

        t_538[k] = -gl_313[k]
                   + f_0 * il_808[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, t_542, t_543, gl_314, gl_315, gl_316, gl_317, \
                         gl_318, il_809, il_810, il_811, il_812, \
                         il_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = -gl_314[k]
                   + f_0 * il_809[k];

        t_540[k] = -2.0 * gl_315[k]
                   + f_0 * il_810[k];

        t_541[k] = -2.0 * gl_316[k]
                   + f_0 * il_811[k];

        t_542[k] = -2.0 * gl_317[k]
                   + f_0 * il_812[k];

        t_543[k] = -2.0 * gl_318[k]
                   + f_0 * il_813[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, t_548, gl_319, gl_320, gl_321, gl_322, \
                         gl_323, il_814, il_815, il_816, il_817, \
                         il_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = -2.0 * gl_319[k]
                   + f_0 * il_814[k];

        t_545[k] = -2.0 * gl_320[k]
                   + f_0 * il_815[k];

        t_546[k] = -2.0 * gl_321[k]
                   + f_0 * il_816[k];

        t_547[k] = -2.0 * gl_322[k]
                   + f_0 * il_817[k];

        t_548[k] = -2.0 * gl_323[k]
                   + f_0 * il_818[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, t_553, gl_324, gl_325, gl_326, gl_327, \
                         gl_328, il_819, il_820, il_821, il_822, \
                         il_823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = -2.0 * gl_324[k]
                   + f_0 * il_819[k];

        t_550[k] = -2.0 * gl_325[k]
                   + f_0 * il_820[k];

        t_551[k] = -2.0 * gl_326[k]
                   + f_0 * il_821[k];

        t_552[k] = -2.0 * gl_327[k]
                   + f_0 * il_822[k];

        t_553[k] = -2.0 * gl_328[k]
                   + f_0 * il_823[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, t_558, gl_329, gl_330, gl_331, gl_332, \
                         gl_333, il_824, il_825, il_826, il_827, \
                         il_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = -2.0 * gl_329[k]
                   + f_0 * il_824[k];

        t_555[k] = -2.0 * gl_330[k]
                   + f_0 * il_825[k];

        t_556[k] = -2.0 * gl_331[k]
                   + f_0 * il_826[k];

        t_557[k] = -2.0 * gl_332[k]
                   + f_0 * il_827[k];

        t_558[k] = -2.0 * gl_333[k]
                   + f_0 * il_828[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, t_563, gl_334, gl_335, gl_336, gl_337, \
                         gl_338, il_829, il_830, il_831, il_832, \
                         il_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = -2.0 * gl_334[k]
                   + f_0 * il_829[k];

        t_560[k] = -2.0 * gl_335[k]
                   + f_0 * il_830[k];

        t_561[k] = -2.0 * gl_336[k]
                   + f_0 * il_831[k];

        t_562[k] = -2.0 * gl_337[k]
                   + f_0 * il_832[k];

        t_563[k] = -2.0 * gl_338[k]
                   + f_0 * il_833[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, t_568, gl_339, gl_340, gl_341, gl_342, \
                         gl_343, il_834, il_835, il_836, il_837, \
                         il_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = -2.0 * gl_339[k]
                   + f_0 * il_834[k];

        t_565[k] = -2.0 * gl_340[k]
                   + f_0 * il_835[k];

        t_566[k] = -2.0 * gl_341[k]
                   + f_0 * il_836[k];

        t_567[k] = -2.0 * gl_342[k]
                   + f_0 * il_837[k];

        t_568[k] = -2.0 * gl_343[k]
                   + f_0 * il_838[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, t_573, gl_344, gl_345, gl_346, gl_347, \
                         gl_348, il_839, il_840, il_841, il_842, \
                         il_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = -2.0 * gl_344[k]
                   + f_0 * il_839[k];

        t_570[k] = -2.0 * gl_345[k]
                   + f_0 * il_840[k];

        t_571[k] = -2.0 * gl_346[k]
                   + f_0 * il_841[k];

        t_572[k] = -2.0 * gl_347[k]
                   + f_0 * il_842[k];

        t_573[k] = -2.0 * gl_348[k]
                   + f_0 * il_843[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, gl_349, gl_350, gl_351, gl_352, \
                         gl_353, il_844, il_845, il_846, il_847, \
                         il_848 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = -2.0 * gl_349[k]
                   + f_0 * il_844[k];

        t_575[k] = -2.0 * gl_350[k]
                   + f_0 * il_845[k];

        t_576[k] = -2.0 * gl_351[k]
                   + f_0 * il_846[k];

        t_577[k] = -2.0 * gl_352[k]
                   + f_0 * il_847[k];

        t_578[k] = -2.0 * gl_353[k]
                   + f_0 * il_848[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, t_582, t_583, gl_354, gl_355, gl_356, gl_357, \
                         gl_358, il_849, il_850, il_851, il_852, \
                         il_853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = -2.0 * gl_354[k]
                   + f_0 * il_849[k];

        t_580[k] = -2.0 * gl_355[k]
                   + f_0 * il_850[k];

        t_581[k] = -2.0 * gl_356[k]
                   + f_0 * il_851[k];

        t_582[k] = -2.0 * gl_357[k]
                   + f_0 * il_852[k];

        t_583[k] = -2.0 * gl_358[k]
                   + f_0 * il_853[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, t_588, gl_359, gl_360, gl_361, gl_362, \
                         gl_363, il_854, il_855, il_856, il_857, \
                         il_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = -2.0 * gl_359[k]
                   + f_0 * il_854[k];

        t_585[k] = -3.0 * gl_360[k]
                   + f_0 * il_855[k];

        t_586[k] = -3.0 * gl_361[k]
                   + f_0 * il_856[k];

        t_587[k] = -3.0 * gl_362[k]
                   + f_0 * il_857[k];

        t_588[k] = -3.0 * gl_363[k]
                   + f_0 * il_858[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, t_593, gl_364, gl_365, gl_366, gl_367, \
                         gl_368, il_859, il_860, il_861, il_862, \
                         il_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = -3.0 * gl_364[k]
                   + f_0 * il_859[k];

        t_590[k] = -3.0 * gl_365[k]
                   + f_0 * il_860[k];

        t_591[k] = -3.0 * gl_366[k]
                   + f_0 * il_861[k];

        t_592[k] = -3.0 * gl_367[k]
                   + f_0 * il_862[k];

        t_593[k] = -3.0 * gl_368[k]
                   + f_0 * il_863[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, t_598, gl_369, gl_370, gl_371, gl_372, \
                         gl_373, il_864, il_865, il_866, il_867, \
                         il_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = -3.0 * gl_369[k]
                   + f_0 * il_864[k];

        t_595[k] = -3.0 * gl_370[k]
                   + f_0 * il_865[k];

        t_596[k] = -3.0 * gl_371[k]
                   + f_0 * il_866[k];

        t_597[k] = -3.0 * gl_372[k]
                   + f_0 * il_867[k];

        t_598[k] = -3.0 * gl_373[k]
                   + f_0 * il_868[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, gl_374, gl_375, gl_376, gl_377, \
                         gl_378, il_869, il_870, il_871, il_872, \
                         il_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = -3.0 * gl_374[k]
                   + f_0 * il_869[k];

        t_600[k] = -3.0 * gl_375[k]
                   + f_0 * il_870[k];

        t_601[k] = -3.0 * gl_376[k]
                   + f_0 * il_871[k];

        t_602[k] = -3.0 * gl_377[k]
                   + f_0 * il_872[k];

        t_603[k] = -3.0 * gl_378[k]
                   + f_0 * il_873[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, gl_379, gl_380, gl_381, gl_382, \
                         gl_383, il_874, il_875, il_876, il_877, \
                         il_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = -3.0 * gl_379[k]
                   + f_0 * il_874[k];

        t_605[k] = -3.0 * gl_380[k]
                   + f_0 * il_875[k];

        t_606[k] = -3.0 * gl_381[k]
                   + f_0 * il_876[k];

        t_607[k] = -3.0 * gl_382[k]
                   + f_0 * il_877[k];

        t_608[k] = -3.0 * gl_383[k]
                   + f_0 * il_878[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, t_613, gl_384, gl_385, gl_386, gl_387, \
                         gl_388, il_879, il_880, il_881, il_882, \
                         il_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = -3.0 * gl_384[k]
                   + f_0 * il_879[k];

        t_610[k] = -3.0 * gl_385[k]
                   + f_0 * il_880[k];

        t_611[k] = -3.0 * gl_386[k]
                   + f_0 * il_881[k];

        t_612[k] = -3.0 * gl_387[k]
                   + f_0 * il_882[k];

        t_613[k] = -3.0 * gl_388[k]
                   + f_0 * il_883[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, gl_389, gl_390, gl_391, gl_392, \
                         gl_393, il_884, il_885, il_886, il_887, \
                         il_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = -3.0 * gl_389[k]
                   + f_0 * il_884[k];

        t_615[k] = -3.0 * gl_390[k]
                   + f_0 * il_885[k];

        t_616[k] = -3.0 * gl_391[k]
                   + f_0 * il_886[k];

        t_617[k] = -3.0 * gl_392[k]
                   + f_0 * il_887[k];

        t_618[k] = -3.0 * gl_393[k]
                   + f_0 * il_888[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, gl_394, gl_395, gl_396, gl_397, \
                         gl_398, il_889, il_890, il_891, il_892, \
                         il_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = -3.0 * gl_394[k]
                   + f_0 * il_889[k];

        t_620[k] = -3.0 * gl_395[k]
                   + f_0 * il_890[k];

        t_621[k] = -3.0 * gl_396[k]
                   + f_0 * il_891[k];

        t_622[k] = -3.0 * gl_397[k]
                   + f_0 * il_892[k];

        t_623[k] = -3.0 * gl_398[k]
                   + f_0 * il_893[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, t_628, gl_399, gl_400, gl_401, gl_402, \
                         gl_403, il_894, il_895, il_896, il_897, \
                         il_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = -3.0 * gl_399[k]
                   + f_0 * il_894[k];

        t_625[k] = -3.0 * gl_400[k]
                   + f_0 * il_895[k];

        t_626[k] = -3.0 * gl_401[k]
                   + f_0 * il_896[k];

        t_627[k] = -3.0 * gl_402[k]
                   + f_0 * il_897[k];

        t_628[k] = -3.0 * gl_403[k]
                   + f_0 * il_898[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, t_633, gl_404, gl_405, gl_406, gl_407, \
                         gl_408, il_899, il_900, il_901, il_902, \
                         il_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = -3.0 * gl_404[k]
                   + f_0 * il_899[k];

        t_630[k] = -4.0 * gl_405[k]
                   + f_0 * il_900[k];

        t_631[k] = -4.0 * gl_406[k]
                   + f_0 * il_901[k];

        t_632[k] = -4.0 * gl_407[k]
                   + f_0 * il_902[k];

        t_633[k] = -4.0 * gl_408[k]
                   + f_0 * il_903[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, gl_409, gl_410, gl_411, gl_412, \
                         gl_413, il_904, il_905, il_906, il_907, \
                         il_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = -4.0 * gl_409[k]
                   + f_0 * il_904[k];

        t_635[k] = -4.0 * gl_410[k]
                   + f_0 * il_905[k];

        t_636[k] = -4.0 * gl_411[k]
                   + f_0 * il_906[k];

        t_637[k] = -4.0 * gl_412[k]
                   + f_0 * il_907[k];

        t_638[k] = -4.0 * gl_413[k]
                   + f_0 * il_908[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, t_643, gl_414, gl_415, gl_416, gl_417, \
                         gl_418, il_909, il_910, il_911, il_912, \
                         il_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = -4.0 * gl_414[k]
                   + f_0 * il_909[k];

        t_640[k] = -4.0 * gl_415[k]
                   + f_0 * il_910[k];

        t_641[k] = -4.0 * gl_416[k]
                   + f_0 * il_911[k];

        t_642[k] = -4.0 * gl_417[k]
                   + f_0 * il_912[k];

        t_643[k] = -4.0 * gl_418[k]
                   + f_0 * il_913[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, t_648, gl_419, gl_420, gl_421, gl_422, \
                         gl_423, il_914, il_915, il_916, il_917, \
                         il_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = -4.0 * gl_419[k]
                   + f_0 * il_914[k];

        t_645[k] = -4.0 * gl_420[k]
                   + f_0 * il_915[k];

        t_646[k] = -4.0 * gl_421[k]
                   + f_0 * il_916[k];

        t_647[k] = -4.0 * gl_422[k]
                   + f_0 * il_917[k];

        t_648[k] = -4.0 * gl_423[k]
                   + f_0 * il_918[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, gl_424, gl_425, gl_426, gl_427, \
                         gl_428, il_919, il_920, il_921, il_922, \
                         il_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = -4.0 * gl_424[k]
                   + f_0 * il_919[k];

        t_650[k] = -4.0 * gl_425[k]
                   + f_0 * il_920[k];

        t_651[k] = -4.0 * gl_426[k]
                   + f_0 * il_921[k];

        t_652[k] = -4.0 * gl_427[k]
                   + f_0 * il_922[k];

        t_653[k] = -4.0 * gl_428[k]
                   + f_0 * il_923[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, t_657, t_658, gl_429, gl_430, gl_431, gl_432, \
                         gl_433, il_924, il_925, il_926, il_927, \
                         il_928 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = -4.0 * gl_429[k]
                   + f_0 * il_924[k];

        t_655[k] = -4.0 * gl_430[k]
                   + f_0 * il_925[k];

        t_656[k] = -4.0 * gl_431[k]
                   + f_0 * il_926[k];

        t_657[k] = -4.0 * gl_432[k]
                   + f_0 * il_927[k];

        t_658[k] = -4.0 * gl_433[k]
                   + f_0 * il_928[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, gl_434, gl_435, gl_436, gl_437, \
                         gl_438, il_929, il_930, il_931, il_932, \
                         il_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = -4.0 * gl_434[k]
                   + f_0 * il_929[k];

        t_660[k] = -4.0 * gl_435[k]
                   + f_0 * il_930[k];

        t_661[k] = -4.0 * gl_436[k]
                   + f_0 * il_931[k];

        t_662[k] = -4.0 * gl_437[k]
                   + f_0 * il_932[k];

        t_663[k] = -4.0 * gl_438[k]
                   + f_0 * il_933[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, t_668, gl_439, gl_440, gl_441, gl_442, \
                         gl_443, il_934, il_935, il_936, il_937, \
                         il_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = -4.0 * gl_439[k]
                   + f_0 * il_934[k];

        t_665[k] = -4.0 * gl_440[k]
                   + f_0 * il_935[k];

        t_666[k] = -4.0 * gl_441[k]
                   + f_0 * il_936[k];

        t_667[k] = -4.0 * gl_442[k]
                   + f_0 * il_937[k];

        t_668[k] = -4.0 * gl_443[k]
                   + f_0 * il_938[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
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
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);

    const auto *gl_444 = buffer.data(gl + 444);
    const auto *gl_445 = buffer.data(gl + 445);
    const auto *gl_446 = buffer.data(gl + 446);
    const auto *gl_447 = buffer.data(gl + 447);
    const auto *gl_448 = buffer.data(gl + 448);
    const auto *gl_449 = buffer.data(gl + 449);
    const auto *gl_450 = buffer.data(gl + 450);
    const auto *gl_451 = buffer.data(gl + 451);
    const auto *gl_452 = buffer.data(gl + 452);
    const auto *gl_453 = buffer.data(gl + 453);
    const auto *gl_454 = buffer.data(gl + 454);
    const auto *gl_455 = buffer.data(gl + 455);
    const auto *gl_456 = buffer.data(gl + 456);
    const auto *gl_457 = buffer.data(gl + 457);
    const auto *gl_458 = buffer.data(gl + 458);
    const auto *gl_459 = buffer.data(gl + 459);
    const auto *gl_460 = buffer.data(gl + 460);
    const auto *gl_461 = buffer.data(gl + 461);
    const auto *gl_462 = buffer.data(gl + 462);
    const auto *gl_463 = buffer.data(gl + 463);
    const auto *gl_464 = buffer.data(gl + 464);
    const auto *gl_465 = buffer.data(gl + 465);
    const auto *gl_466 = buffer.data(gl + 466);
    const auto *gl_467 = buffer.data(gl + 467);
    const auto *gl_468 = buffer.data(gl + 468);
    const auto *gl_469 = buffer.data(gl + 469);
    const auto *gl_470 = buffer.data(gl + 470);
    const auto *gl_471 = buffer.data(gl + 471);
    const auto *gl_472 = buffer.data(gl + 472);
    const auto *gl_473 = buffer.data(gl + 473);
    const auto *gl_474 = buffer.data(gl + 474);
    const auto *gl_475 = buffer.data(gl + 475);
    const auto *gl_476 = buffer.data(gl + 476);
    const auto *gl_477 = buffer.data(gl + 477);
    const auto *gl_478 = buffer.data(gl + 478);
    const auto *gl_479 = buffer.data(gl + 479);
    const auto *gl_480 = buffer.data(gl + 480);
    const auto *gl_481 = buffer.data(gl + 481);
    const auto *gl_482 = buffer.data(gl + 482);
    const auto *gl_483 = buffer.data(gl + 483);
    const auto *gl_484 = buffer.data(gl + 484);
    const auto *gl_485 = buffer.data(gl + 485);
    const auto *gl_486 = buffer.data(gl + 486);
    const auto *gl_487 = buffer.data(gl + 487);
    const auto *gl_488 = buffer.data(gl + 488);
    const auto *gl_489 = buffer.data(gl + 489);
    const auto *gl_490 = buffer.data(gl + 490);
    const auto *gl_491 = buffer.data(gl + 491);
    const auto *gl_492 = buffer.data(gl + 492);
    const auto *gl_493 = buffer.data(gl + 493);
    const auto *gl_494 = buffer.data(gl + 494);
    const auto *gl_495 = buffer.data(gl + 495);
    const auto *gl_496 = buffer.data(gl + 496);
    const auto *gl_497 = buffer.data(gl + 497);
    const auto *gl_498 = buffer.data(gl + 498);
    const auto *gl_499 = buffer.data(gl + 499);
    const auto *gl_500 = buffer.data(gl + 500);
    const auto *gl_501 = buffer.data(gl + 501);
    const auto *gl_502 = buffer.data(gl + 502);
    const auto *gl_503 = buffer.data(gl + 503);
    const auto *gl_504 = buffer.data(gl + 504);
    const auto *gl_505 = buffer.data(gl + 505);
    const auto *gl_506 = buffer.data(gl + 506);
    const auto *gl_507 = buffer.data(gl + 507);
    const auto *gl_508 = buffer.data(gl + 508);
    const auto *gl_509 = buffer.data(gl + 509);
    const auto *gl_510 = buffer.data(gl + 510);
    const auto *gl_511 = buffer.data(gl + 511);
    const auto *gl_512 = buffer.data(gl + 512);
    const auto *gl_513 = buffer.data(gl + 513);
    const auto *gl_514 = buffer.data(gl + 514);
    const auto *gl_515 = buffer.data(gl + 515);
    const auto *gl_516 = buffer.data(gl + 516);
    const auto *gl_517 = buffer.data(gl + 517);
    const auto *gl_518 = buffer.data(gl + 518);
    const auto *gl_519 = buffer.data(gl + 519);
    const auto *gl_520 = buffer.data(gl + 520);
    const auto *gl_521 = buffer.data(gl + 521);
    const auto *gl_522 = buffer.data(gl + 522);
    const auto *gl_523 = buffer.data(gl + 523);
    const auto *gl_524 = buffer.data(gl + 524);
    const auto *gl_525 = buffer.data(gl + 525);
    const auto *gl_526 = buffer.data(gl + 526);
    const auto *gl_527 = buffer.data(gl + 527);
    const auto *gl_528 = buffer.data(gl + 528);
    const auto *gl_529 = buffer.data(gl + 529);
    const auto *gl_530 = buffer.data(gl + 530);
    const auto *gl_531 = buffer.data(gl + 531);
    const auto *gl_532 = buffer.data(gl + 532);
    const auto *gl_533 = buffer.data(gl + 533);
    const auto *gl_534 = buffer.data(gl + 534);
    const auto *gl_535 = buffer.data(gl + 535);
    const auto *gl_536 = buffer.data(gl + 536);
    const auto *gl_537 = buffer.data(gl + 537);
    const auto *gl_538 = buffer.data(gl + 538);
    const auto *gl_539 = buffer.data(gl + 539);
    const auto *gl_540 = buffer.data(gl + 540);
    const auto *gl_541 = buffer.data(gl + 541);
    const auto *gl_542 = buffer.data(gl + 542);
    const auto *gl_543 = buffer.data(gl + 543);
    const auto *gl_544 = buffer.data(gl + 544);
    const auto *gl_545 = buffer.data(gl + 545);
    const auto *gl_546 = buffer.data(gl + 546);
    const auto *gl_547 = buffer.data(gl + 547);
    const auto *gl_548 = buffer.data(gl + 548);
    const auto *gl_549 = buffer.data(gl + 549);
    const auto *gl_550 = buffer.data(gl + 550);
    const auto *gl_551 = buffer.data(gl + 551);
    const auto *gl_552 = buffer.data(gl + 552);
    const auto *gl_553 = buffer.data(gl + 553);
    const auto *gl_554 = buffer.data(gl + 554);
    const auto *gl_555 = buffer.data(gl + 555);
    const auto *gl_556 = buffer.data(gl + 556);
    const auto *gl_557 = buffer.data(gl + 557);
    const auto *gl_558 = buffer.data(gl + 558);
    const auto *gl_559 = buffer.data(gl + 559);
    const auto *gl_560 = buffer.data(gl + 560);
    const auto *gl_561 = buffer.data(gl + 561);
    const auto *gl_562 = buffer.data(gl + 562);
    const auto *gl_563 = buffer.data(gl + 563);
    const auto *gl_564 = buffer.data(gl + 564);

    const auto *il_939 = buffer.data(il + 939);
    const auto *il_940 = buffer.data(il + 940);
    const auto *il_941 = buffer.data(il + 941);
    const auto *il_942 = buffer.data(il + 942);
    const auto *il_943 = buffer.data(il + 943);
    const auto *il_944 = buffer.data(il + 944);
    const auto *il_990 = buffer.data(il + 990);
    const auto *il_991 = buffer.data(il + 991);
    const auto *il_992 = buffer.data(il + 992);
    const auto *il_993 = buffer.data(il + 993);
    const auto *il_994 = buffer.data(il + 994);
    const auto *il_995 = buffer.data(il + 995);
    const auto *il_996 = buffer.data(il + 996);
    const auto *il_997 = buffer.data(il + 997);
    const auto *il_998 = buffer.data(il + 998);
    const auto *il_999 = buffer.data(il + 999);
    const auto *il_1000 = buffer.data(il + 1000);
    const auto *il_1001 = buffer.data(il + 1001);
    const auto *il_1002 = buffer.data(il + 1002);
    const auto *il_1003 = buffer.data(il + 1003);
    const auto *il_1004 = buffer.data(il + 1004);
    const auto *il_1005 = buffer.data(il + 1005);
    const auto *il_1006 = buffer.data(il + 1006);
    const auto *il_1007 = buffer.data(il + 1007);
    const auto *il_1008 = buffer.data(il + 1008);
    const auto *il_1009 = buffer.data(il + 1009);
    const auto *il_1010 = buffer.data(il + 1010);
    const auto *il_1011 = buffer.data(il + 1011);
    const auto *il_1012 = buffer.data(il + 1012);
    const auto *il_1013 = buffer.data(il + 1013);
    const auto *il_1014 = buffer.data(il + 1014);
    const auto *il_1015 = buffer.data(il + 1015);
    const auto *il_1016 = buffer.data(il + 1016);
    const auto *il_1017 = buffer.data(il + 1017);
    const auto *il_1018 = buffer.data(il + 1018);
    const auto *il_1019 = buffer.data(il + 1019);
    const auto *il_1020 = buffer.data(il + 1020);
    const auto *il_1021 = buffer.data(il + 1021);
    const auto *il_1022 = buffer.data(il + 1022);
    const auto *il_1023 = buffer.data(il + 1023);
    const auto *il_1024 = buffer.data(il + 1024);
    const auto *il_1025 = buffer.data(il + 1025);
    const auto *il_1026 = buffer.data(il + 1026);
    const auto *il_1027 = buffer.data(il + 1027);
    const auto *il_1028 = buffer.data(il + 1028);
    const auto *il_1029 = buffer.data(il + 1029);
    const auto *il_1030 = buffer.data(il + 1030);
    const auto *il_1031 = buffer.data(il + 1031);
    const auto *il_1032 = buffer.data(il + 1032);
    const auto *il_1033 = buffer.data(il + 1033);
    const auto *il_1034 = buffer.data(il + 1034);
    const auto *il_1035 = buffer.data(il + 1035);
    const auto *il_1036 = buffer.data(il + 1036);
    const auto *il_1037 = buffer.data(il + 1037);
    const auto *il_1038 = buffer.data(il + 1038);
    const auto *il_1039 = buffer.data(il + 1039);
    const auto *il_1040 = buffer.data(il + 1040);
    const auto *il_1041 = buffer.data(il + 1041);
    const auto *il_1042 = buffer.data(il + 1042);
    const auto *il_1043 = buffer.data(il + 1043);
    const auto *il_1044 = buffer.data(il + 1044);
    const auto *il_1045 = buffer.data(il + 1045);
    const auto *il_1046 = buffer.data(il + 1046);
    const auto *il_1047 = buffer.data(il + 1047);
    const auto *il_1048 = buffer.data(il + 1048);
    const auto *il_1049 = buffer.data(il + 1049);
    const auto *il_1050 = buffer.data(il + 1050);
    const auto *il_1051 = buffer.data(il + 1051);
    const auto *il_1052 = buffer.data(il + 1052);
    const auto *il_1053 = buffer.data(il + 1053);
    const auto *il_1054 = buffer.data(il + 1054);
    const auto *il_1055 = buffer.data(il + 1055);
    const auto *il_1056 = buffer.data(il + 1056);
    const auto *il_1057 = buffer.data(il + 1057);
    const auto *il_1058 = buffer.data(il + 1058);
    const auto *il_1059 = buffer.data(il + 1059);
    const auto *il_1060 = buffer.data(il + 1060);
    const auto *il_1061 = buffer.data(il + 1061);
    const auto *il_1062 = buffer.data(il + 1062);
    const auto *il_1063 = buffer.data(il + 1063);
    const auto *il_1064 = buffer.data(il + 1064);
    const auto *il_1065 = buffer.data(il + 1065);
    const auto *il_1066 = buffer.data(il + 1066);
    const auto *il_1067 = buffer.data(il + 1067);
    const auto *il_1068 = buffer.data(il + 1068);
    const auto *il_1069 = buffer.data(il + 1069);
    const auto *il_1070 = buffer.data(il + 1070);
    const auto *il_1071 = buffer.data(il + 1071);
    const auto *il_1072 = buffer.data(il + 1072);
    const auto *il_1073 = buffer.data(il + 1073);
    const auto *il_1074 = buffer.data(il + 1074);
    const auto *il_1075 = buffer.data(il + 1075);
    const auto *il_1076 = buffer.data(il + 1076);
    const auto *il_1077 = buffer.data(il + 1077);
    const auto *il_1078 = buffer.data(il + 1078);
    const auto *il_1079 = buffer.data(il + 1079);
    const auto *il_1080 = buffer.data(il + 1080);
    const auto *il_1081 = buffer.data(il + 1081);
    const auto *il_1082 = buffer.data(il + 1082);
    const auto *il_1083 = buffer.data(il + 1083);
    const auto *il_1084 = buffer.data(il + 1084);
    const auto *il_1085 = buffer.data(il + 1085);
    const auto *il_1086 = buffer.data(il + 1086);
    const auto *il_1087 = buffer.data(il + 1087);
    const auto *il_1088 = buffer.data(il + 1088);
    const auto *il_1089 = buffer.data(il + 1089);
    const auto *il_1090 = buffer.data(il + 1090);
    const auto *il_1091 = buffer.data(il + 1091);
    const auto *il_1092 = buffer.data(il + 1092);
    const auto *il_1093 = buffer.data(il + 1093);
    const auto *il_1094 = buffer.data(il + 1094);
    const auto *il_1095 = buffer.data(il + 1095);
    const auto *il_1096 = buffer.data(il + 1096);
    const auto *il_1097 = buffer.data(il + 1097);
    const auto *il_1098 = buffer.data(il + 1098);
    const auto *il_1099 = buffer.data(il + 1099);
    const auto *il_1100 = buffer.data(il + 1100);
    const auto *il_1101 = buffer.data(il + 1101);
    const auto *il_1102 = buffer.data(il + 1102);
    const auto *il_1103 = buffer.data(il + 1103);
    const auto *il_1104 = buffer.data(il + 1104);
    const auto *il_1105 = buffer.data(il + 1105);
    const auto *il_1106 = buffer.data(il + 1106);
    const auto *il_1107 = buffer.data(il + 1107);
    const auto *il_1108 = buffer.data(il + 1108);
    const auto *il_1109 = buffer.data(il + 1109);
    const auto *il_1110 = buffer.data(il + 1110);
    const auto *il_1111 = buffer.data(il + 1111);
    const auto *il_1112 = buffer.data(il + 1112);
    const auto *il_1113 = buffer.data(il + 1113);
    const auto *il_1114 = buffer.data(il + 1114);
    const auto *il_1115 = buffer.data(il + 1115);
    const auto *il_1116 = buffer.data(il + 1116);
    const auto *il_1117 = buffer.data(il + 1117);
    const auto *il_1118 = buffer.data(il + 1118);
    const auto *il_1119 = buffer.data(il + 1119);
    const auto *il_1120 = buffer.data(il + 1120);
    const auto *il_1121 = buffer.data(il + 1121);
    const auto *il_1122 = buffer.data(il + 1122);
    const auto *il_1123 = buffer.data(il + 1123);
    const auto *il_1124 = buffer.data(il + 1124);
    const auto *il_1125 = buffer.data(il + 1125);
    const auto *il_1126 = buffer.data(il + 1126);
    const auto *il_1127 = buffer.data(il + 1127);
    const auto *il_1128 = buffer.data(il + 1128);
    const auto *il_1129 = buffer.data(il + 1129);
    const auto *il_1130 = buffer.data(il + 1130);
    const auto *il_1131 = buffer.data(il + 1131);
    const auto *il_1132 = buffer.data(il + 1132);
    const auto *il_1133 = buffer.data(il + 1133);
    const auto *il_1134 = buffer.data(il + 1134);
    const auto *il_1135 = buffer.data(il + 1135);
    const auto *il_1136 = buffer.data(il + 1136);
    const auto *il_1137 = buffer.data(il + 1137);
    const auto *il_1138 = buffer.data(il + 1138);
    const auto *il_1139 = buffer.data(il + 1139);
    const auto *il_1140 = buffer.data(il + 1140);
    const auto *il_1141 = buffer.data(il + 1141);
    const auto *il_1142 = buffer.data(il + 1142);
    const auto *il_1143 = buffer.data(il + 1143);
    const auto *il_1144 = buffer.data(il + 1144);
    const auto *il_1145 = buffer.data(il + 1145);
    const auto *il_1146 = buffer.data(il + 1146);
    const auto *il_1147 = buffer.data(il + 1147);
    const auto *il_1148 = buffer.data(il + 1148);
    const auto *il_1149 = buffer.data(il + 1149);

#pragma omp simd aligned(t_669, t_670, t_671, t_672, t_673, gl_444, gl_445, gl_446, gl_447, \
                         gl_448, il_939, il_940, il_941, il_942, \
                         il_943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = -4.0 * gl_444[k]
                   + f_0 * il_939[k];

        t_670[k] = -4.0 * gl_445[k]
                   + f_0 * il_940[k];

        t_671[k] = -4.0 * gl_446[k]
                   + f_0 * il_941[k];

        t_672[k] = -4.0 * gl_447[k]
                   + f_0 * il_942[k];

        t_673[k] = -4.0 * gl_448[k]
                   + f_0 * il_943[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, t_677, t_678, t_679, t_680, gl_449, il_944, \
                         il_990, il_991, il_992, il_993, il_994, \
                         il_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = -4.0 * gl_449[k]
                   + f_0 * il_944[k];

        t_675[k] = f_0 * il_990[k];

        t_676[k] = f_0 * il_991[k];

        t_677[k] = f_0 * il_992[k];

        t_678[k] = f_0 * il_993[k];

        t_679[k] = f_0 * il_994[k];

        t_680[k] = f_0 * il_995[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, t_684, t_685, t_686, t_687, t_688, il_996, \
                         il_997, il_998, il_999, il_1000, il_1001, il_1002, \
                         il_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = f_0 * il_996[k];

        t_682[k] = f_0 * il_997[k];

        t_683[k] = f_0 * il_998[k];

        t_684[k] = f_0 * il_999[k];

        t_685[k] = f_0 * il_1000[k];

        t_686[k] = f_0 * il_1001[k];

        t_687[k] = f_0 * il_1002[k];

        t_688[k] = f_0 * il_1003[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, t_693, t_694, t_695, t_696, il_1004, \
                         il_1005, il_1006, il_1007, il_1008, il_1009, il_1010, \
                         il_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = f_0 * il_1004[k];

        t_690[k] = f_0 * il_1005[k];

        t_691[k] = f_0 * il_1006[k];

        t_692[k] = f_0 * il_1007[k];

        t_693[k] = f_0 * il_1008[k];

        t_694[k] = f_0 * il_1009[k];

        t_695[k] = f_0 * il_1010[k];

        t_696[k] = f_0 * il_1011[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, t_701, t_702, t_703, t_704, il_1012, \
                         il_1013, il_1014, il_1015, il_1016, il_1017, il_1018, \
                         il_1019 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_0 * il_1012[k];

        t_698[k] = f_0 * il_1013[k];

        t_699[k] = f_0 * il_1014[k];

        t_700[k] = f_0 * il_1015[k];

        t_701[k] = f_0 * il_1016[k];

        t_702[k] = f_0 * il_1017[k];

        t_703[k] = f_0 * il_1018[k];

        t_704[k] = f_0 * il_1019[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, t_710, t_711, t_712, il_1020, \
                         il_1021, il_1022, il_1023, il_1024, il_1025, il_1026, \
                         il_1027 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_0 * il_1020[k];

        t_706[k] = f_0 * il_1021[k];

        t_707[k] = f_0 * il_1022[k];

        t_708[k] = f_0 * il_1023[k];

        t_709[k] = f_0 * il_1024[k];

        t_710[k] = f_0 * il_1025[k];

        t_711[k] = f_0 * il_1026[k];

        t_712[k] = f_0 * il_1027[k];
    }

#pragma omp simd aligned(t_713, t_714, t_715, t_716, t_717, t_718, t_719, il_1028, il_1029, \
                         il_1030, il_1031, il_1032, il_1033, il_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_713[k] = f_0 * il_1028[k];

        t_714[k] = f_0 * il_1029[k];

        t_715[k] = f_0 * il_1030[k];

        t_716[k] = f_0 * il_1031[k];

        t_717[k] = f_0 * il_1032[k];

        t_718[k] = f_0 * il_1033[k];

        t_719[k] = f_0 * il_1034[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, gl_450, gl_451, gl_452, gl_453, \
                         gl_454, il_1035, il_1036, il_1037, il_1038, \
                         il_1039 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -gl_450[k]
                   + f_0 * il_1035[k];

        t_721[k] = -gl_451[k]
                   + f_0 * il_1036[k];

        t_722[k] = -gl_452[k]
                   + f_0 * il_1037[k];

        t_723[k] = -gl_453[k]
                   + f_0 * il_1038[k];

        t_724[k] = -gl_454[k]
                   + f_0 * il_1039[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, gl_455, gl_456, gl_457, gl_458, \
                         gl_459, il_1040, il_1041, il_1042, il_1043, \
                         il_1044 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -gl_455[k]
                   + f_0 * il_1040[k];

        t_726[k] = -gl_456[k]
                   + f_0 * il_1041[k];

        t_727[k] = -gl_457[k]
                   + f_0 * il_1042[k];

        t_728[k] = -gl_458[k]
                   + f_0 * il_1043[k];

        t_729[k] = -gl_459[k]
                   + f_0 * il_1044[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, gl_460, gl_461, gl_462, gl_463, \
                         gl_464, il_1045, il_1046, il_1047, il_1048, \
                         il_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -gl_460[k]
                   + f_0 * il_1045[k];

        t_731[k] = -gl_461[k]
                   + f_0 * il_1046[k];

        t_732[k] = -gl_462[k]
                   + f_0 * il_1047[k];

        t_733[k] = -gl_463[k]
                   + f_0 * il_1048[k];

        t_734[k] = -gl_464[k]
                   + f_0 * il_1049[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, gl_465, gl_466, gl_467, gl_468, \
                         gl_469, il_1050, il_1051, il_1052, il_1053, \
                         il_1054 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -gl_465[k]
                   + f_0 * il_1050[k];

        t_736[k] = -gl_466[k]
                   + f_0 * il_1051[k];

        t_737[k] = -gl_467[k]
                   + f_0 * il_1052[k];

        t_738[k] = -gl_468[k]
                   + f_0 * il_1053[k];

        t_739[k] = -gl_469[k]
                   + f_0 * il_1054[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, gl_470, gl_471, gl_472, gl_473, \
                         gl_474, il_1055, il_1056, il_1057, il_1058, \
                         il_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -gl_470[k]
                   + f_0 * il_1055[k];

        t_741[k] = -gl_471[k]
                   + f_0 * il_1056[k];

        t_742[k] = -gl_472[k]
                   + f_0 * il_1057[k];

        t_743[k] = -gl_473[k]
                   + f_0 * il_1058[k];

        t_744[k] = -gl_474[k]
                   + f_0 * il_1059[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, gl_475, gl_476, gl_477, gl_478, \
                         gl_479, il_1060, il_1061, il_1062, il_1063, \
                         il_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -gl_475[k]
                   + f_0 * il_1060[k];

        t_746[k] = -gl_476[k]
                   + f_0 * il_1061[k];

        t_747[k] = -gl_477[k]
                   + f_0 * il_1062[k];

        t_748[k] = -gl_478[k]
                   + f_0 * il_1063[k];

        t_749[k] = -gl_479[k]
                   + f_0 * il_1064[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, gl_480, gl_481, gl_482, gl_483, \
                         gl_484, il_1065, il_1066, il_1067, il_1068, \
                         il_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -gl_480[k]
                   + f_0 * il_1065[k];

        t_751[k] = -gl_481[k]
                   + f_0 * il_1066[k];

        t_752[k] = -gl_482[k]
                   + f_0 * il_1067[k];

        t_753[k] = -gl_483[k]
                   + f_0 * il_1068[k];

        t_754[k] = -gl_484[k]
                   + f_0 * il_1069[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, gl_485, gl_486, gl_487, gl_488, \
                         gl_489, il_1070, il_1071, il_1072, il_1073, \
                         il_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -gl_485[k]
                   + f_0 * il_1070[k];

        t_756[k] = -gl_486[k]
                   + f_0 * il_1071[k];

        t_757[k] = -gl_487[k]
                   + f_0 * il_1072[k];

        t_758[k] = -gl_488[k]
                   + f_0 * il_1073[k];

        t_759[k] = -gl_489[k]
                   + f_0 * il_1074[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, t_764, gl_490, gl_491, gl_492, gl_493, \
                         gl_494, il_1075, il_1076, il_1077, il_1078, \
                         il_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = -gl_490[k]
                   + f_0 * il_1075[k];

        t_761[k] = -gl_491[k]
                   + f_0 * il_1076[k];

        t_762[k] = -gl_492[k]
                   + f_0 * il_1077[k];

        t_763[k] = -gl_493[k]
                   + f_0 * il_1078[k];

        t_764[k] = -gl_494[k]
                   + f_0 * il_1079[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, t_769, gl_495, gl_496, gl_497, gl_498, \
                         gl_499, il_1080, il_1081, il_1082, il_1083, \
                         il_1084 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = -2.0 * gl_495[k]
                   + f_0 * il_1080[k];

        t_766[k] = -2.0 * gl_496[k]
                   + f_0 * il_1081[k];

        t_767[k] = -2.0 * gl_497[k]
                   + f_0 * il_1082[k];

        t_768[k] = -2.0 * gl_498[k]
                   + f_0 * il_1083[k];

        t_769[k] = -2.0 * gl_499[k]
                   + f_0 * il_1084[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, gl_500, gl_501, gl_502, gl_503, \
                         gl_504, il_1085, il_1086, il_1087, il_1088, \
                         il_1089 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = -2.0 * gl_500[k]
                   + f_0 * il_1085[k];

        t_771[k] = -2.0 * gl_501[k]
                   + f_0 * il_1086[k];

        t_772[k] = -2.0 * gl_502[k]
                   + f_0 * il_1087[k];

        t_773[k] = -2.0 * gl_503[k]
                   + f_0 * il_1088[k];

        t_774[k] = -2.0 * gl_504[k]
                   + f_0 * il_1089[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, gl_505, gl_506, gl_507, gl_508, \
                         gl_509, il_1090, il_1091, il_1092, il_1093, \
                         il_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = -2.0 * gl_505[k]
                   + f_0 * il_1090[k];

        t_776[k] = -2.0 * gl_506[k]
                   + f_0 * il_1091[k];

        t_777[k] = -2.0 * gl_507[k]
                   + f_0 * il_1092[k];

        t_778[k] = -2.0 * gl_508[k]
                   + f_0 * il_1093[k];

        t_779[k] = -2.0 * gl_509[k]
                   + f_0 * il_1094[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, t_784, gl_510, gl_511, gl_512, gl_513, \
                         gl_514, il_1095, il_1096, il_1097, il_1098, \
                         il_1099 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = -2.0 * gl_510[k]
                   + f_0 * il_1095[k];

        t_781[k] = -2.0 * gl_511[k]
                   + f_0 * il_1096[k];

        t_782[k] = -2.0 * gl_512[k]
                   + f_0 * il_1097[k];

        t_783[k] = -2.0 * gl_513[k]
                   + f_0 * il_1098[k];

        t_784[k] = -2.0 * gl_514[k]
                   + f_0 * il_1099[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, t_789, gl_515, gl_516, gl_517, gl_518, \
                         gl_519, il_1100, il_1101, il_1102, il_1103, \
                         il_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = -2.0 * gl_515[k]
                   + f_0 * il_1100[k];

        t_786[k] = -2.0 * gl_516[k]
                   + f_0 * il_1101[k];

        t_787[k] = -2.0 * gl_517[k]
                   + f_0 * il_1102[k];

        t_788[k] = -2.0 * gl_518[k]
                   + f_0 * il_1103[k];

        t_789[k] = -2.0 * gl_519[k]
                   + f_0 * il_1104[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, t_794, gl_520, gl_521, gl_522, gl_523, \
                         gl_524, il_1105, il_1106, il_1107, il_1108, \
                         il_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = -2.0 * gl_520[k]
                   + f_0 * il_1105[k];

        t_791[k] = -2.0 * gl_521[k]
                   + f_0 * il_1106[k];

        t_792[k] = -2.0 * gl_522[k]
                   + f_0 * il_1107[k];

        t_793[k] = -2.0 * gl_523[k]
                   + f_0 * il_1108[k];

        t_794[k] = -2.0 * gl_524[k]
                   + f_0 * il_1109[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, gl_525, gl_526, gl_527, gl_528, \
                         gl_529, il_1110, il_1111, il_1112, il_1113, \
                         il_1114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = -2.0 * gl_525[k]
                   + f_0 * il_1110[k];

        t_796[k] = -2.0 * gl_526[k]
                   + f_0 * il_1111[k];

        t_797[k] = -2.0 * gl_527[k]
                   + f_0 * il_1112[k];

        t_798[k] = -2.0 * gl_528[k]
                   + f_0 * il_1113[k];

        t_799[k] = -2.0 * gl_529[k]
                   + f_0 * il_1114[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, gl_530, gl_531, gl_532, gl_533, \
                         gl_534, il_1115, il_1116, il_1117, il_1118, \
                         il_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = -2.0 * gl_530[k]
                   + f_0 * il_1115[k];

        t_801[k] = -2.0 * gl_531[k]
                   + f_0 * il_1116[k];

        t_802[k] = -2.0 * gl_532[k]
                   + f_0 * il_1117[k];

        t_803[k] = -2.0 * gl_533[k]
                   + f_0 * il_1118[k];

        t_804[k] = -2.0 * gl_534[k]
                   + f_0 * il_1119[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, gl_535, gl_536, gl_537, gl_538, \
                         gl_539, il_1120, il_1121, il_1122, il_1123, \
                         il_1124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = -2.0 * gl_535[k]
                   + f_0 * il_1120[k];

        t_806[k] = -2.0 * gl_536[k]
                   + f_0 * il_1121[k];

        t_807[k] = -2.0 * gl_537[k]
                   + f_0 * il_1122[k];

        t_808[k] = -2.0 * gl_538[k]
                   + f_0 * il_1123[k];

        t_809[k] = -2.0 * gl_539[k]
                   + f_0 * il_1124[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, gl_540, gl_541, gl_542, gl_543, \
                         gl_544, il_1125, il_1126, il_1127, il_1128, \
                         il_1129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = -3.0 * gl_540[k]
                   + f_0 * il_1125[k];

        t_811[k] = -3.0 * gl_541[k]
                   + f_0 * il_1126[k];

        t_812[k] = -3.0 * gl_542[k]
                   + f_0 * il_1127[k];

        t_813[k] = -3.0 * gl_543[k]
                   + f_0 * il_1128[k];

        t_814[k] = -3.0 * gl_544[k]
                   + f_0 * il_1129[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, gl_545, gl_546, gl_547, gl_548, \
                         gl_549, il_1130, il_1131, il_1132, il_1133, \
                         il_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = -3.0 * gl_545[k]
                   + f_0 * il_1130[k];

        t_816[k] = -3.0 * gl_546[k]
                   + f_0 * il_1131[k];

        t_817[k] = -3.0 * gl_547[k]
                   + f_0 * il_1132[k];

        t_818[k] = -3.0 * gl_548[k]
                   + f_0 * il_1133[k];

        t_819[k] = -3.0 * gl_549[k]
                   + f_0 * il_1134[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, t_824, gl_550, gl_551, gl_552, gl_553, \
                         gl_554, il_1135, il_1136, il_1137, il_1138, \
                         il_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = -3.0 * gl_550[k]
                   + f_0 * il_1135[k];

        t_821[k] = -3.0 * gl_551[k]
                   + f_0 * il_1136[k];

        t_822[k] = -3.0 * gl_552[k]
                   + f_0 * il_1137[k];

        t_823[k] = -3.0 * gl_553[k]
                   + f_0 * il_1138[k];

        t_824[k] = -3.0 * gl_554[k]
                   + f_0 * il_1139[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, gl_555, gl_556, gl_557, gl_558, \
                         gl_559, il_1140, il_1141, il_1142, il_1143, \
                         il_1144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = -3.0 * gl_555[k]
                   + f_0 * il_1140[k];

        t_826[k] = -3.0 * gl_556[k]
                   + f_0 * il_1141[k];

        t_827[k] = -3.0 * gl_557[k]
                   + f_0 * il_1142[k];

        t_828[k] = -3.0 * gl_558[k]
                   + f_0 * il_1143[k];

        t_829[k] = -3.0 * gl_559[k]
                   + f_0 * il_1144[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, gl_560, gl_561, gl_562, gl_563, \
                         gl_564, il_1145, il_1146, il_1147, il_1148, \
                         il_1149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = -3.0 * gl_560[k]
                   + f_0 * il_1145[k];

        t_831[k] = -3.0 * gl_561[k]
                   + f_0 * il_1146[k];

        t_832[k] = -3.0 * gl_562[k]
                   + f_0 * il_1147[k];

        t_833[k] = -3.0 * gl_563[k]
                   + f_0 * il_1148[k];

        t_834[k] = -3.0 * gl_564[k]
                   + f_0 * il_1149[k];
    }
}

static auto
compute_prim_geom_10_hl_electron_repulsion_2_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t il,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gl_565 = buffer.data(gl + 565);
    const auto *gl_566 = buffer.data(gl + 566);
    const auto *gl_567 = buffer.data(gl + 567);
    const auto *gl_568 = buffer.data(gl + 568);
    const auto *gl_569 = buffer.data(gl + 569);
    const auto *gl_570 = buffer.data(gl + 570);
    const auto *gl_571 = buffer.data(gl + 571);
    const auto *gl_572 = buffer.data(gl + 572);
    const auto *gl_573 = buffer.data(gl + 573);
    const auto *gl_574 = buffer.data(gl + 574);
    const auto *gl_575 = buffer.data(gl + 575);
    const auto *gl_576 = buffer.data(gl + 576);
    const auto *gl_577 = buffer.data(gl + 577);
    const auto *gl_578 = buffer.data(gl + 578);
    const auto *gl_579 = buffer.data(gl + 579);
    const auto *gl_580 = buffer.data(gl + 580);
    const auto *gl_581 = buffer.data(gl + 581);
    const auto *gl_582 = buffer.data(gl + 582);
    const auto *gl_583 = buffer.data(gl + 583);
    const auto *gl_584 = buffer.data(gl + 584);
    const auto *gl_585 = buffer.data(gl + 585);
    const auto *gl_586 = buffer.data(gl + 586);
    const auto *gl_587 = buffer.data(gl + 587);
    const auto *gl_588 = buffer.data(gl + 588);
    const auto *gl_589 = buffer.data(gl + 589);
    const auto *gl_590 = buffer.data(gl + 590);
    const auto *gl_591 = buffer.data(gl + 591);
    const auto *gl_592 = buffer.data(gl + 592);
    const auto *gl_593 = buffer.data(gl + 593);
    const auto *gl_594 = buffer.data(gl + 594);
    const auto *gl_595 = buffer.data(gl + 595);
    const auto *gl_596 = buffer.data(gl + 596);
    const auto *gl_597 = buffer.data(gl + 597);
    const auto *gl_598 = buffer.data(gl + 598);
    const auto *gl_599 = buffer.data(gl + 599);
    const auto *gl_600 = buffer.data(gl + 600);
    const auto *gl_601 = buffer.data(gl + 601);
    const auto *gl_602 = buffer.data(gl + 602);
    const auto *gl_603 = buffer.data(gl + 603);
    const auto *gl_604 = buffer.data(gl + 604);
    const auto *gl_605 = buffer.data(gl + 605);
    const auto *gl_606 = buffer.data(gl + 606);
    const auto *gl_607 = buffer.data(gl + 607);
    const auto *gl_608 = buffer.data(gl + 608);
    const auto *gl_609 = buffer.data(gl + 609);
    const auto *gl_610 = buffer.data(gl + 610);
    const auto *gl_611 = buffer.data(gl + 611);
    const auto *gl_612 = buffer.data(gl + 612);
    const auto *gl_613 = buffer.data(gl + 613);
    const auto *gl_614 = buffer.data(gl + 614);
    const auto *gl_615 = buffer.data(gl + 615);
    const auto *gl_616 = buffer.data(gl + 616);
    const auto *gl_617 = buffer.data(gl + 617);
    const auto *gl_618 = buffer.data(gl + 618);
    const auto *gl_619 = buffer.data(gl + 619);
    const auto *gl_620 = buffer.data(gl + 620);
    const auto *gl_621 = buffer.data(gl + 621);
    const auto *gl_622 = buffer.data(gl + 622);
    const auto *gl_623 = buffer.data(gl + 623);
    const auto *gl_624 = buffer.data(gl + 624);
    const auto *gl_625 = buffer.data(gl + 625);
    const auto *gl_626 = buffer.data(gl + 626);
    const auto *gl_627 = buffer.data(gl + 627);
    const auto *gl_628 = buffer.data(gl + 628);
    const auto *gl_629 = buffer.data(gl + 629);
    const auto *gl_630 = buffer.data(gl + 630);
    const auto *gl_631 = buffer.data(gl + 631);
    const auto *gl_632 = buffer.data(gl + 632);
    const auto *gl_633 = buffer.data(gl + 633);
    const auto *gl_634 = buffer.data(gl + 634);
    const auto *gl_635 = buffer.data(gl + 635);
    const auto *gl_636 = buffer.data(gl + 636);
    const auto *gl_637 = buffer.data(gl + 637);
    const auto *gl_638 = buffer.data(gl + 638);
    const auto *gl_639 = buffer.data(gl + 639);
    const auto *gl_640 = buffer.data(gl + 640);
    const auto *gl_641 = buffer.data(gl + 641);
    const auto *gl_642 = buffer.data(gl + 642);
    const auto *gl_643 = buffer.data(gl + 643);
    const auto *gl_644 = buffer.data(gl + 644);
    const auto *gl_645 = buffer.data(gl + 645);
    const auto *gl_646 = buffer.data(gl + 646);
    const auto *gl_647 = buffer.data(gl + 647);
    const auto *gl_648 = buffer.data(gl + 648);
    const auto *gl_649 = buffer.data(gl + 649);
    const auto *gl_650 = buffer.data(gl + 650);
    const auto *gl_651 = buffer.data(gl + 651);
    const auto *gl_652 = buffer.data(gl + 652);
    const auto *gl_653 = buffer.data(gl + 653);
    const auto *gl_654 = buffer.data(gl + 654);
    const auto *gl_655 = buffer.data(gl + 655);
    const auto *gl_656 = buffer.data(gl + 656);
    const auto *gl_657 = buffer.data(gl + 657);
    const auto *gl_658 = buffer.data(gl + 658);
    const auto *gl_659 = buffer.data(gl + 659);
    const auto *gl_660 = buffer.data(gl + 660);
    const auto *gl_661 = buffer.data(gl + 661);
    const auto *gl_662 = buffer.data(gl + 662);
    const auto *gl_663 = buffer.data(gl + 663);
    const auto *gl_664 = buffer.data(gl + 664);
    const auto *gl_665 = buffer.data(gl + 665);
    const auto *gl_666 = buffer.data(gl + 666);
    const auto *gl_667 = buffer.data(gl + 667);
    const auto *gl_668 = buffer.data(gl + 668);
    const auto *gl_669 = buffer.data(gl + 669);
    const auto *gl_670 = buffer.data(gl + 670);
    const auto *gl_671 = buffer.data(gl + 671);
    const auto *gl_672 = buffer.data(gl + 672);
    const auto *gl_673 = buffer.data(gl + 673);
    const auto *gl_674 = buffer.data(gl + 674);

    const auto *il_1150 = buffer.data(il + 1150);
    const auto *il_1151 = buffer.data(il + 1151);
    const auto *il_1152 = buffer.data(il + 1152);
    const auto *il_1153 = buffer.data(il + 1153);
    const auto *il_1154 = buffer.data(il + 1154);
    const auto *il_1155 = buffer.data(il + 1155);
    const auto *il_1156 = buffer.data(il + 1156);
    const auto *il_1157 = buffer.data(il + 1157);
    const auto *il_1158 = buffer.data(il + 1158);
    const auto *il_1159 = buffer.data(il + 1159);
    const auto *il_1160 = buffer.data(il + 1160);
    const auto *il_1161 = buffer.data(il + 1161);
    const auto *il_1162 = buffer.data(il + 1162);
    const auto *il_1163 = buffer.data(il + 1163);
    const auto *il_1164 = buffer.data(il + 1164);
    const auto *il_1165 = buffer.data(il + 1165);
    const auto *il_1166 = buffer.data(il + 1166);
    const auto *il_1167 = buffer.data(il + 1167);
    const auto *il_1168 = buffer.data(il + 1168);
    const auto *il_1169 = buffer.data(il + 1169);
    const auto *il_1170 = buffer.data(il + 1170);
    const auto *il_1171 = buffer.data(il + 1171);
    const auto *il_1172 = buffer.data(il + 1172);
    const auto *il_1173 = buffer.data(il + 1173);
    const auto *il_1174 = buffer.data(il + 1174);
    const auto *il_1175 = buffer.data(il + 1175);
    const auto *il_1176 = buffer.data(il + 1176);
    const auto *il_1177 = buffer.data(il + 1177);
    const auto *il_1178 = buffer.data(il + 1178);
    const auto *il_1179 = buffer.data(il + 1179);
    const auto *il_1180 = buffer.data(il + 1180);
    const auto *il_1181 = buffer.data(il + 1181);
    const auto *il_1182 = buffer.data(il + 1182);
    const auto *il_1183 = buffer.data(il + 1183);
    const auto *il_1184 = buffer.data(il + 1184);
    const auto *il_1185 = buffer.data(il + 1185);
    const auto *il_1186 = buffer.data(il + 1186);
    const auto *il_1187 = buffer.data(il + 1187);
    const auto *il_1188 = buffer.data(il + 1188);
    const auto *il_1189 = buffer.data(il + 1189);
    const auto *il_1190 = buffer.data(il + 1190);
    const auto *il_1191 = buffer.data(il + 1191);
    const auto *il_1192 = buffer.data(il + 1192);
    const auto *il_1193 = buffer.data(il + 1193);
    const auto *il_1194 = buffer.data(il + 1194);
    const auto *il_1195 = buffer.data(il + 1195);
    const auto *il_1196 = buffer.data(il + 1196);
    const auto *il_1197 = buffer.data(il + 1197);
    const auto *il_1198 = buffer.data(il + 1198);
    const auto *il_1199 = buffer.data(il + 1199);
    const auto *il_1200 = buffer.data(il + 1200);
    const auto *il_1201 = buffer.data(il + 1201);
    const auto *il_1202 = buffer.data(il + 1202);
    const auto *il_1203 = buffer.data(il + 1203);
    const auto *il_1204 = buffer.data(il + 1204);
    const auto *il_1205 = buffer.data(il + 1205);
    const auto *il_1206 = buffer.data(il + 1206);
    const auto *il_1207 = buffer.data(il + 1207);
    const auto *il_1208 = buffer.data(il + 1208);
    const auto *il_1209 = buffer.data(il + 1209);
    const auto *il_1210 = buffer.data(il + 1210);
    const auto *il_1211 = buffer.data(il + 1211);
    const auto *il_1212 = buffer.data(il + 1212);
    const auto *il_1213 = buffer.data(il + 1213);
    const auto *il_1214 = buffer.data(il + 1214);
    const auto *il_1215 = buffer.data(il + 1215);
    const auto *il_1216 = buffer.data(il + 1216);
    const auto *il_1217 = buffer.data(il + 1217);
    const auto *il_1218 = buffer.data(il + 1218);
    const auto *il_1219 = buffer.data(il + 1219);
    const auto *il_1220 = buffer.data(il + 1220);
    const auto *il_1221 = buffer.data(il + 1221);
    const auto *il_1222 = buffer.data(il + 1222);
    const auto *il_1223 = buffer.data(il + 1223);
    const auto *il_1224 = buffer.data(il + 1224);
    const auto *il_1225 = buffer.data(il + 1225);
    const auto *il_1226 = buffer.data(il + 1226);
    const auto *il_1227 = buffer.data(il + 1227);
    const auto *il_1228 = buffer.data(il + 1228);
    const auto *il_1229 = buffer.data(il + 1229);
    const auto *il_1230 = buffer.data(il + 1230);
    const auto *il_1231 = buffer.data(il + 1231);
    const auto *il_1232 = buffer.data(il + 1232);
    const auto *il_1233 = buffer.data(il + 1233);
    const auto *il_1234 = buffer.data(il + 1234);
    const auto *il_1235 = buffer.data(il + 1235);
    const auto *il_1236 = buffer.data(il + 1236);
    const auto *il_1237 = buffer.data(il + 1237);
    const auto *il_1238 = buffer.data(il + 1238);
    const auto *il_1239 = buffer.data(il + 1239);
    const auto *il_1240 = buffer.data(il + 1240);
    const auto *il_1241 = buffer.data(il + 1241);
    const auto *il_1242 = buffer.data(il + 1242);
    const auto *il_1243 = buffer.data(il + 1243);
    const auto *il_1244 = buffer.data(il + 1244);
    const auto *il_1245 = buffer.data(il + 1245);
    const auto *il_1246 = buffer.data(il + 1246);
    const auto *il_1247 = buffer.data(il + 1247);
    const auto *il_1248 = buffer.data(il + 1248);
    const auto *il_1249 = buffer.data(il + 1249);
    const auto *il_1250 = buffer.data(il + 1250);
    const auto *il_1251 = buffer.data(il + 1251);
    const auto *il_1252 = buffer.data(il + 1252);
    const auto *il_1253 = buffer.data(il + 1253);
    const auto *il_1254 = buffer.data(il + 1254);
    const auto *il_1255 = buffer.data(il + 1255);
    const auto *il_1256 = buffer.data(il + 1256);
    const auto *il_1257 = buffer.data(il + 1257);
    const auto *il_1258 = buffer.data(il + 1258);
    const auto *il_1259 = buffer.data(il + 1259);

#pragma omp simd aligned(t_835, t_836, t_837, t_838, t_839, gl_565, gl_566, gl_567, gl_568, \
                         gl_569, il_1150, il_1151, il_1152, il_1153, \
                         il_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = -3.0 * gl_565[k]
                   + f_0 * il_1150[k];

        t_836[k] = -3.0 * gl_566[k]
                   + f_0 * il_1151[k];

        t_837[k] = -3.0 * gl_567[k]
                   + f_0 * il_1152[k];

        t_838[k] = -3.0 * gl_568[k]
                   + f_0 * il_1153[k];

        t_839[k] = -3.0 * gl_569[k]
                   + f_0 * il_1154[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, gl_570, gl_571, gl_572, gl_573, \
                         gl_574, il_1155, il_1156, il_1157, il_1158, \
                         il_1159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = -3.0 * gl_570[k]
                   + f_0 * il_1155[k];

        t_841[k] = -3.0 * gl_571[k]
                   + f_0 * il_1156[k];

        t_842[k] = -3.0 * gl_572[k]
                   + f_0 * il_1157[k];

        t_843[k] = -3.0 * gl_573[k]
                   + f_0 * il_1158[k];

        t_844[k] = -3.0 * gl_574[k]
                   + f_0 * il_1159[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, t_849, gl_575, gl_576, gl_577, gl_578, \
                         gl_579, il_1160, il_1161, il_1162, il_1163, \
                         il_1164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = -3.0 * gl_575[k]
                   + f_0 * il_1160[k];

        t_846[k] = -3.0 * gl_576[k]
                   + f_0 * il_1161[k];

        t_847[k] = -3.0 * gl_577[k]
                   + f_0 * il_1162[k];

        t_848[k] = -3.0 * gl_578[k]
                   + f_0 * il_1163[k];

        t_849[k] = -3.0 * gl_579[k]
                   + f_0 * il_1164[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, gl_580, gl_581, gl_582, gl_583, \
                         gl_584, il_1165, il_1166, il_1167, il_1168, \
                         il_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = -3.0 * gl_580[k]
                   + f_0 * il_1165[k];

        t_851[k] = -3.0 * gl_581[k]
                   + f_0 * il_1166[k];

        t_852[k] = -3.0 * gl_582[k]
                   + f_0 * il_1167[k];

        t_853[k] = -3.0 * gl_583[k]
                   + f_0 * il_1168[k];

        t_854[k] = -3.0 * gl_584[k]
                   + f_0 * il_1169[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, gl_585, gl_586, gl_587, gl_588, \
                         gl_589, il_1170, il_1171, il_1172, il_1173, \
                         il_1174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = -4.0 * gl_585[k]
                   + f_0 * il_1170[k];

        t_856[k] = -4.0 * gl_586[k]
                   + f_0 * il_1171[k];

        t_857[k] = -4.0 * gl_587[k]
                   + f_0 * il_1172[k];

        t_858[k] = -4.0 * gl_588[k]
                   + f_0 * il_1173[k];

        t_859[k] = -4.0 * gl_589[k]
                   + f_0 * il_1174[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, t_864, gl_590, gl_591, gl_592, gl_593, \
                         gl_594, il_1175, il_1176, il_1177, il_1178, \
                         il_1179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = -4.0 * gl_590[k]
                   + f_0 * il_1175[k];

        t_861[k] = -4.0 * gl_591[k]
                   + f_0 * il_1176[k];

        t_862[k] = -4.0 * gl_592[k]
                   + f_0 * il_1177[k];

        t_863[k] = -4.0 * gl_593[k]
                   + f_0 * il_1178[k];

        t_864[k] = -4.0 * gl_594[k]
                   + f_0 * il_1179[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, t_869, gl_595, gl_596, gl_597, gl_598, \
                         gl_599, il_1180, il_1181, il_1182, il_1183, \
                         il_1184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = -4.0 * gl_595[k]
                   + f_0 * il_1180[k];

        t_866[k] = -4.0 * gl_596[k]
                   + f_0 * il_1181[k];

        t_867[k] = -4.0 * gl_597[k]
                   + f_0 * il_1182[k];

        t_868[k] = -4.0 * gl_598[k]
                   + f_0 * il_1183[k];

        t_869[k] = -4.0 * gl_599[k]
                   + f_0 * il_1184[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, t_874, gl_600, gl_601, gl_602, gl_603, \
                         gl_604, il_1185, il_1186, il_1187, il_1188, \
                         il_1189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = -4.0 * gl_600[k]
                   + f_0 * il_1185[k];

        t_871[k] = -4.0 * gl_601[k]
                   + f_0 * il_1186[k];

        t_872[k] = -4.0 * gl_602[k]
                   + f_0 * il_1187[k];

        t_873[k] = -4.0 * gl_603[k]
                   + f_0 * il_1188[k];

        t_874[k] = -4.0 * gl_604[k]
                   + f_0 * il_1189[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, t_879, gl_605, gl_606, gl_607, gl_608, \
                         gl_609, il_1190, il_1191, il_1192, il_1193, \
                         il_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = -4.0 * gl_605[k]
                   + f_0 * il_1190[k];

        t_876[k] = -4.0 * gl_606[k]
                   + f_0 * il_1191[k];

        t_877[k] = -4.0 * gl_607[k]
                   + f_0 * il_1192[k];

        t_878[k] = -4.0 * gl_608[k]
                   + f_0 * il_1193[k];

        t_879[k] = -4.0 * gl_609[k]
                   + f_0 * il_1194[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, t_884, gl_610, gl_611, gl_612, gl_613, \
                         gl_614, il_1195, il_1196, il_1197, il_1198, \
                         il_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = -4.0 * gl_610[k]
                   + f_0 * il_1195[k];

        t_881[k] = -4.0 * gl_611[k]
                   + f_0 * il_1196[k];

        t_882[k] = -4.0 * gl_612[k]
                   + f_0 * il_1197[k];

        t_883[k] = -4.0 * gl_613[k]
                   + f_0 * il_1198[k];

        t_884[k] = -4.0 * gl_614[k]
                   + f_0 * il_1199[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, gl_615, gl_616, gl_617, gl_618, \
                         gl_619, il_1200, il_1201, il_1202, il_1203, \
                         il_1204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = -4.0 * gl_615[k]
                   + f_0 * il_1200[k];

        t_886[k] = -4.0 * gl_616[k]
                   + f_0 * il_1201[k];

        t_887[k] = -4.0 * gl_617[k]
                   + f_0 * il_1202[k];

        t_888[k] = -4.0 * gl_618[k]
                   + f_0 * il_1203[k];

        t_889[k] = -4.0 * gl_619[k]
                   + f_0 * il_1204[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, gl_620, gl_621, gl_622, gl_623, \
                         gl_624, il_1205, il_1206, il_1207, il_1208, \
                         il_1209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = -4.0 * gl_620[k]
                   + f_0 * il_1205[k];

        t_891[k] = -4.0 * gl_621[k]
                   + f_0 * il_1206[k];

        t_892[k] = -4.0 * gl_622[k]
                   + f_0 * il_1207[k];

        t_893[k] = -4.0 * gl_623[k]
                   + f_0 * il_1208[k];

        t_894[k] = -4.0 * gl_624[k]
                   + f_0 * il_1209[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, t_899, gl_625, gl_626, gl_627, gl_628, \
                         gl_629, il_1210, il_1211, il_1212, il_1213, \
                         il_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = -4.0 * gl_625[k]
                   + f_0 * il_1210[k];

        t_896[k] = -4.0 * gl_626[k]
                   + f_0 * il_1211[k];

        t_897[k] = -4.0 * gl_627[k]
                   + f_0 * il_1212[k];

        t_898[k] = -4.0 * gl_628[k]
                   + f_0 * il_1213[k];

        t_899[k] = -4.0 * gl_629[k]
                   + f_0 * il_1214[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, t_904, gl_630, gl_631, gl_632, gl_633, \
                         gl_634, il_1215, il_1216, il_1217, il_1218, \
                         il_1219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = -5.0 * gl_630[k]
                   + f_0 * il_1215[k];

        t_901[k] = -5.0 * gl_631[k]
                   + f_0 * il_1216[k];

        t_902[k] = -5.0 * gl_632[k]
                   + f_0 * il_1217[k];

        t_903[k] = -5.0 * gl_633[k]
                   + f_0 * il_1218[k];

        t_904[k] = -5.0 * gl_634[k]
                   + f_0 * il_1219[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, t_909, gl_635, gl_636, gl_637, gl_638, \
                         gl_639, il_1220, il_1221, il_1222, il_1223, \
                         il_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = -5.0 * gl_635[k]
                   + f_0 * il_1220[k];

        t_906[k] = -5.0 * gl_636[k]
                   + f_0 * il_1221[k];

        t_907[k] = -5.0 * gl_637[k]
                   + f_0 * il_1222[k];

        t_908[k] = -5.0 * gl_638[k]
                   + f_0 * il_1223[k];

        t_909[k] = -5.0 * gl_639[k]
                   + f_0 * il_1224[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, t_914, gl_640, gl_641, gl_642, gl_643, \
                         gl_644, il_1225, il_1226, il_1227, il_1228, \
                         il_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = -5.0 * gl_640[k]
                   + f_0 * il_1225[k];

        t_911[k] = -5.0 * gl_641[k]
                   + f_0 * il_1226[k];

        t_912[k] = -5.0 * gl_642[k]
                   + f_0 * il_1227[k];

        t_913[k] = -5.0 * gl_643[k]
                   + f_0 * il_1228[k];

        t_914[k] = -5.0 * gl_644[k]
                   + f_0 * il_1229[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, t_919, gl_645, gl_646, gl_647, gl_648, \
                         gl_649, il_1230, il_1231, il_1232, il_1233, \
                         il_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = -5.0 * gl_645[k]
                   + f_0 * il_1230[k];

        t_916[k] = -5.0 * gl_646[k]
                   + f_0 * il_1231[k];

        t_917[k] = -5.0 * gl_647[k]
                   + f_0 * il_1232[k];

        t_918[k] = -5.0 * gl_648[k]
                   + f_0 * il_1233[k];

        t_919[k] = -5.0 * gl_649[k]
                   + f_0 * il_1234[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, t_924, gl_650, gl_651, gl_652, gl_653, \
                         gl_654, il_1235, il_1236, il_1237, il_1238, \
                         il_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = -5.0 * gl_650[k]
                   + f_0 * il_1235[k];

        t_921[k] = -5.0 * gl_651[k]
                   + f_0 * il_1236[k];

        t_922[k] = -5.0 * gl_652[k]
                   + f_0 * il_1237[k];

        t_923[k] = -5.0 * gl_653[k]
                   + f_0 * il_1238[k];

        t_924[k] = -5.0 * gl_654[k]
                   + f_0 * il_1239[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, t_929, gl_655, gl_656, gl_657, gl_658, \
                         gl_659, il_1240, il_1241, il_1242, il_1243, \
                         il_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = -5.0 * gl_655[k]
                   + f_0 * il_1240[k];

        t_926[k] = -5.0 * gl_656[k]
                   + f_0 * il_1241[k];

        t_927[k] = -5.0 * gl_657[k]
                   + f_0 * il_1242[k];

        t_928[k] = -5.0 * gl_658[k]
                   + f_0 * il_1243[k];

        t_929[k] = -5.0 * gl_659[k]
                   + f_0 * il_1244[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, gl_660, gl_661, gl_662, gl_663, \
                         gl_664, il_1245, il_1246, il_1247, il_1248, \
                         il_1249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = -5.0 * gl_660[k]
                   + f_0 * il_1245[k];

        t_931[k] = -5.0 * gl_661[k]
                   + f_0 * il_1246[k];

        t_932[k] = -5.0 * gl_662[k]
                   + f_0 * il_1247[k];

        t_933[k] = -5.0 * gl_663[k]
                   + f_0 * il_1248[k];

        t_934[k] = -5.0 * gl_664[k]
                   + f_0 * il_1249[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, t_939, gl_665, gl_666, gl_667, gl_668, \
                         gl_669, il_1250, il_1251, il_1252, il_1253, \
                         il_1254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = -5.0 * gl_665[k]
                   + f_0 * il_1250[k];

        t_936[k] = -5.0 * gl_666[k]
                   + f_0 * il_1251[k];

        t_937[k] = -5.0 * gl_667[k]
                   + f_0 * il_1252[k];

        t_938[k] = -5.0 * gl_668[k]
                   + f_0 * il_1253[k];

        t_939[k] = -5.0 * gl_669[k]
                   + f_0 * il_1254[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, t_944, gl_670, gl_671, gl_672, gl_673, \
                         gl_674, il_1255, il_1256, il_1257, il_1258, \
                         il_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = -5.0 * gl_670[k]
                   + f_0 * il_1255[k];

        t_941[k] = -5.0 * gl_671[k]
                   + f_0 * il_1256[k];

        t_942[k] = -5.0 * gl_672[k]
                   + f_0 * il_1257[k];

        t_943[k] = -5.0 * gl_673[k]
                   + f_0 * il_1258[k];

        t_944[k] = -5.0 * gl_674[k]
                   + f_0 * il_1259[k];
    }
}

auto
compute_prim_geom_10_hl_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t gl, const size_t il,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hl_electron_repulsion_2_piece0(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_2_piece1(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_2_piece2(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_2_piece3(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_2_piece4(buffer, target, gl, il, ncols, alpha);

    compute_prim_geom_10_hl_electron_repulsion_2_piece5(buffer, target, gl, il, ncols, alpha);
}

}  // namespace simdt2ceri
