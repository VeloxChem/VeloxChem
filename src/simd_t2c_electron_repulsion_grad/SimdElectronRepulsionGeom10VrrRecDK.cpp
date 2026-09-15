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


#include "SimdElectronRepulsionGeom10VrrRecDK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_dk_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t pk, const size_t fk,
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

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_11 = buffer.data(pk + 11);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_16 = buffer.data(pk + 16);
    const auto *pk_17 = buffer.data(pk + 17);
    const auto *pk_18 = buffer.data(pk + 18);
    const auto *pk_19 = buffer.data(pk + 19);
    const auto *pk_20 = buffer.data(pk + 20);
    const auto *pk_21 = buffer.data(pk + 21);
    const auto *pk_22 = buffer.data(pk + 22);
    const auto *pk_23 = buffer.data(pk + 23);
    const auto *pk_24 = buffer.data(pk + 24);
    const auto *pk_25 = buffer.data(pk + 25);
    const auto *pk_26 = buffer.data(pk + 26);
    const auto *pk_27 = buffer.data(pk + 27);
    const auto *pk_28 = buffer.data(pk + 28);
    const auto *pk_29 = buffer.data(pk + 29);
    const auto *pk_30 = buffer.data(pk + 30);
    const auto *pk_31 = buffer.data(pk + 31);
    const auto *pk_32 = buffer.data(pk + 32);
    const auto *pk_33 = buffer.data(pk + 33);
    const auto *pk_34 = buffer.data(pk + 34);
    const auto *pk_35 = buffer.data(pk + 35);
    const auto *pk_36 = buffer.data(pk + 36);
    const auto *pk_37 = buffer.data(pk + 37);
    const auto *pk_38 = buffer.data(pk + 38);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_40 = buffer.data(pk + 40);
    const auto *pk_41 = buffer.data(pk + 41);
    const auto *pk_42 = buffer.data(pk + 42);
    const auto *pk_43 = buffer.data(pk + 43);
    const auto *pk_44 = buffer.data(pk + 44);
    const auto *pk_45 = buffer.data(pk + 45);
    const auto *pk_46 = buffer.data(pk + 46);
    const auto *pk_47 = buffer.data(pk + 47);
    const auto *pk_48 = buffer.data(pk + 48);
    const auto *pk_49 = buffer.data(pk + 49);
    const auto *pk_50 = buffer.data(pk + 50);
    const auto *pk_51 = buffer.data(pk + 51);
    const auto *pk_52 = buffer.data(pk + 52);
    const auto *pk_53 = buffer.data(pk + 53);
    const auto *pk_54 = buffer.data(pk + 54);
    const auto *pk_55 = buffer.data(pk + 55);
    const auto *pk_56 = buffer.data(pk + 56);
    const auto *pk_57 = buffer.data(pk + 57);
    const auto *pk_58 = buffer.data(pk + 58);
    const auto *pk_59 = buffer.data(pk + 59);
    const auto *pk_60 = buffer.data(pk + 60);
    const auto *pk_61 = buffer.data(pk + 61);
    const auto *pk_62 = buffer.data(pk + 62);
    const auto *pk_63 = buffer.data(pk + 63);
    const auto *pk_64 = buffer.data(pk + 64);
    const auto *pk_65 = buffer.data(pk + 65);
    const auto *pk_66 = buffer.data(pk + 66);
    const auto *pk_67 = buffer.data(pk + 67);
    const auto *pk_68 = buffer.data(pk + 68);
    const auto *pk_69 = buffer.data(pk + 69);
    const auto *pk_70 = buffer.data(pk + 70);
    const auto *pk_71 = buffer.data(pk + 71);
    const auto *pk_72 = buffer.data(pk + 72);
    const auto *pk_73 = buffer.data(pk + 73);
    const auto *pk_74 = buffer.data(pk + 74);
    const auto *pk_75 = buffer.data(pk + 75);
    const auto *pk_76 = buffer.data(pk + 76);
    const auto *pk_77 = buffer.data(pk + 77);
    const auto *pk_78 = buffer.data(pk + 78);
    const auto *pk_79 = buffer.data(pk + 79);
    const auto *pk_80 = buffer.data(pk + 80);
    const auto *pk_81 = buffer.data(pk + 81);
    const auto *pk_82 = buffer.data(pk + 82);
    const auto *pk_83 = buffer.data(pk + 83);
    const auto *pk_84 = buffer.data(pk + 84);
    const auto *pk_85 = buffer.data(pk + 85);
    const auto *pk_86 = buffer.data(pk + 86);
    const auto *pk_87 = buffer.data(pk + 87);
    const auto *pk_88 = buffer.data(pk + 88);
    const auto *pk_89 = buffer.data(pk + 89);
    const auto *pk_90 = buffer.data(pk + 90);
    const auto *pk_91 = buffer.data(pk + 91);
    const auto *pk_92 = buffer.data(pk + 92);
    const auto *pk_93 = buffer.data(pk + 93);
    const auto *pk_94 = buffer.data(pk + 94);
    const auto *pk_95 = buffer.data(pk + 95);
    const auto *pk_96 = buffer.data(pk + 96);
    const auto *pk_97 = buffer.data(pk + 97);
    const auto *pk_98 = buffer.data(pk + 98);
    const auto *pk_99 = buffer.data(pk + 99);
    const auto *pk_100 = buffer.data(pk + 100);
    const auto *pk_101 = buffer.data(pk + 101);
    const auto *pk_102 = buffer.data(pk + 102);
    const auto *pk_103 = buffer.data(pk + 103);
    const auto *pk_104 = buffer.data(pk + 104);
    const auto *pk_105 = buffer.data(pk + 105);
    const auto *pk_106 = buffer.data(pk + 106);
    const auto *pk_107 = buffer.data(pk + 107);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);
    const auto *fk_36 = buffer.data(fk + 36);
    const auto *fk_37 = buffer.data(fk + 37);
    const auto *fk_38 = buffer.data(fk + 38);
    const auto *fk_39 = buffer.data(fk + 39);
    const auto *fk_40 = buffer.data(fk + 40);
    const auto *fk_41 = buffer.data(fk + 41);
    const auto *fk_42 = buffer.data(fk + 42);
    const auto *fk_43 = buffer.data(fk + 43);
    const auto *fk_44 = buffer.data(fk + 44);
    const auto *fk_45 = buffer.data(fk + 45);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_47 = buffer.data(fk + 47);
    const auto *fk_48 = buffer.data(fk + 48);
    const auto *fk_49 = buffer.data(fk + 49);
    const auto *fk_50 = buffer.data(fk + 50);
    const auto *fk_51 = buffer.data(fk + 51);
    const auto *fk_52 = buffer.data(fk + 52);
    const auto *fk_53 = buffer.data(fk + 53);
    const auto *fk_54 = buffer.data(fk + 54);
    const auto *fk_55 = buffer.data(fk + 55);
    const auto *fk_56 = buffer.data(fk + 56);
    const auto *fk_57 = buffer.data(fk + 57);
    const auto *fk_58 = buffer.data(fk + 58);
    const auto *fk_59 = buffer.data(fk + 59);
    const auto *fk_60 = buffer.data(fk + 60);
    const auto *fk_61 = buffer.data(fk + 61);
    const auto *fk_62 = buffer.data(fk + 62);
    const auto *fk_63 = buffer.data(fk + 63);
    const auto *fk_64 = buffer.data(fk + 64);
    const auto *fk_65 = buffer.data(fk + 65);
    const auto *fk_66 = buffer.data(fk + 66);
    const auto *fk_67 = buffer.data(fk + 67);
    const auto *fk_68 = buffer.data(fk + 68);
    const auto *fk_69 = buffer.data(fk + 69);
    const auto *fk_70 = buffer.data(fk + 70);
    const auto *fk_71 = buffer.data(fk + 71);
    const auto *fk_72 = buffer.data(fk + 72);
    const auto *fk_73 = buffer.data(fk + 73);
    const auto *fk_74 = buffer.data(fk + 74);
    const auto *fk_75 = buffer.data(fk + 75);
    const auto *fk_76 = buffer.data(fk + 76);
    const auto *fk_77 = buffer.data(fk + 77);
    const auto *fk_78 = buffer.data(fk + 78);
    const auto *fk_79 = buffer.data(fk + 79);
    const auto *fk_80 = buffer.data(fk + 80);
    const auto *fk_81 = buffer.data(fk + 81);
    const auto *fk_82 = buffer.data(fk + 82);
    const auto *fk_83 = buffer.data(fk + 83);
    const auto *fk_84 = buffer.data(fk + 84);
    const auto *fk_85 = buffer.data(fk + 85);
    const auto *fk_86 = buffer.data(fk + 86);
    const auto *fk_87 = buffer.data(fk + 87);
    const auto *fk_88 = buffer.data(fk + 88);
    const auto *fk_89 = buffer.data(fk + 89);
    const auto *fk_90 = buffer.data(fk + 90);
    const auto *fk_91 = buffer.data(fk + 91);
    const auto *fk_92 = buffer.data(fk + 92);
    const auto *fk_93 = buffer.data(fk + 93);
    const auto *fk_94 = buffer.data(fk + 94);
    const auto *fk_95 = buffer.data(fk + 95);
    const auto *fk_96 = buffer.data(fk + 96);
    const auto *fk_97 = buffer.data(fk + 97);
    const auto *fk_98 = buffer.data(fk + 98);
    const auto *fk_99 = buffer.data(fk + 99);
    const auto *fk_100 = buffer.data(fk + 100);
    const auto *fk_101 = buffer.data(fk + 101);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_106 = buffer.data(fk + 106);
    const auto *fk_107 = buffer.data(fk + 107);
    const auto *fk_108 = buffer.data(fk + 108);
    const auto *fk_109 = buffer.data(fk + 109);
    const auto *fk_110 = buffer.data(fk + 110);
    const auto *fk_111 = buffer.data(fk + 111);
    const auto *fk_112 = buffer.data(fk + 112);
    const auto *fk_113 = buffer.data(fk + 113);
    const auto *fk_114 = buffer.data(fk + 114);
    const auto *fk_115 = buffer.data(fk + 115);
    const auto *fk_116 = buffer.data(fk + 116);
    const auto *fk_117 = buffer.data(fk + 117);
    const auto *fk_118 = buffer.data(fk + 118);
    const auto *fk_119 = buffer.data(fk + 119);
    const auto *fk_120 = buffer.data(fk + 120);
    const auto *fk_121 = buffer.data(fk + 121);
    const auto *fk_122 = buffer.data(fk + 122);
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_124 = buffer.data(fk + 124);
    const auto *fk_125 = buffer.data(fk + 125);
    const auto *fk_126 = buffer.data(fk + 126);
    const auto *fk_127 = buffer.data(fk + 127);
    const auto *fk_128 = buffer.data(fk + 128);
    const auto *fk_129 = buffer.data(fk + 129);
    const auto *fk_130 = buffer.data(fk + 130);
    const auto *fk_131 = buffer.data(fk + 131);
    const auto *fk_132 = buffer.data(fk + 132);
    const auto *fk_133 = buffer.data(fk + 133);
    const auto *fk_134 = buffer.data(fk + 134);
    const auto *fk_135 = buffer.data(fk + 135);
    const auto *fk_136 = buffer.data(fk + 136);
    const auto *fk_137 = buffer.data(fk + 137);
    const auto *fk_138 = buffer.data(fk + 138);
    const auto *fk_139 = buffer.data(fk + 139);
    const auto *fk_140 = buffer.data(fk + 140);
    const auto *fk_141 = buffer.data(fk + 141);
    const auto *fk_142 = buffer.data(fk + 142);
    const auto *fk_143 = buffer.data(fk + 143);
    const auto *fk_144 = buffer.data(fk + 144);
    const auto *fk_145 = buffer.data(fk + 145);
    const auto *fk_146 = buffer.data(fk + 146);
    const auto *fk_147 = buffer.data(fk + 147);
    const auto *fk_148 = buffer.data(fk + 148);
    const auto *fk_149 = buffer.data(fk + 149);
    const auto *fk_150 = buffer.data(fk + 150);
    const auto *fk_151 = buffer.data(fk + 151);
    const auto *fk_152 = buffer.data(fk + 152);
    const auto *fk_153 = buffer.data(fk + 153);
    const auto *fk_154 = buffer.data(fk + 154);
    const auto *fk_155 = buffer.data(fk + 155);
    const auto *fk_156 = buffer.data(fk + 156);
    const auto *fk_157 = buffer.data(fk + 157);
    const auto *fk_158 = buffer.data(fk + 158);
    const auto *fk_159 = buffer.data(fk + 159);
    const auto *fk_160 = buffer.data(fk + 160);
    const auto *fk_161 = buffer.data(fk + 161);
    const auto *fk_162 = buffer.data(fk + 162);
    const auto *fk_163 = buffer.data(fk + 163);
    const auto *fk_164 = buffer.data(fk + 164);
    const auto *fk_165 = buffer.data(fk + 165);
    const auto *fk_166 = buffer.data(fk + 166);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pk_0, pk_1, pk_2, pk_3, pk_4, fk_0, fk_1, \
                         fk_2, fk_3, fk_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -2.0 * pk_0[k]
                 + f_0 * fk_0[k];

        t_1[k] = -2.0 * pk_1[k]
                 + f_0 * fk_1[k];

        t_2[k] = -2.0 * pk_2[k]
                 + f_0 * fk_2[k];

        t_3[k] = -2.0 * pk_3[k]
                 + f_0 * fk_3[k];

        t_4[k] = -2.0 * pk_4[k]
                 + f_0 * fk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pk_5, pk_6, pk_7, pk_8, pk_9, fk_5, fk_6, \
                         fk_7, fk_8, fk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -2.0 * pk_5[k]
                 + f_0 * fk_5[k];

        t_6[k] = -2.0 * pk_6[k]
                 + f_0 * fk_6[k];

        t_7[k] = -2.0 * pk_7[k]
                 + f_0 * fk_7[k];

        t_8[k] = -2.0 * pk_8[k]
                 + f_0 * fk_8[k];

        t_9[k] = -2.0 * pk_9[k]
                 + f_0 * fk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pk_10, pk_11, pk_12, pk_13, pk_14, \
                         fk_10, fk_11, fk_12, fk_13, fk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -2.0 * pk_10[k]
                  + f_0 * fk_10[k];

        t_11[k] = -2.0 * pk_11[k]
                  + f_0 * fk_11[k];

        t_12[k] = -2.0 * pk_12[k]
                  + f_0 * fk_12[k];

        t_13[k] = -2.0 * pk_13[k]
                  + f_0 * fk_13[k];

        t_14[k] = -2.0 * pk_14[k]
                  + f_0 * fk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pk_15, pk_16, pk_17, pk_18, pk_19, \
                         fk_15, fk_16, fk_17, fk_18, fk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -2.0 * pk_15[k]
                  + f_0 * fk_15[k];

        t_16[k] = -2.0 * pk_16[k]
                  + f_0 * fk_16[k];

        t_17[k] = -2.0 * pk_17[k]
                  + f_0 * fk_17[k];

        t_18[k] = -2.0 * pk_18[k]
                  + f_0 * fk_18[k];

        t_19[k] = -2.0 * pk_19[k]
                  + f_0 * fk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pk_20, pk_21, pk_22, pk_23, pk_24, \
                         fk_20, fk_21, fk_22, fk_23, fk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -2.0 * pk_20[k]
                  + f_0 * fk_20[k];

        t_21[k] = -2.0 * pk_21[k]
                  + f_0 * fk_21[k];

        t_22[k] = -2.0 * pk_22[k]
                  + f_0 * fk_22[k];

        t_23[k] = -2.0 * pk_23[k]
                  + f_0 * fk_23[k];

        t_24[k] = -2.0 * pk_24[k]
                  + f_0 * fk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pk_25, pk_26, pk_27, pk_28, pk_29, \
                         fk_25, fk_26, fk_27, fk_28, fk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -2.0 * pk_25[k]
                  + f_0 * fk_25[k];

        t_26[k] = -2.0 * pk_26[k]
                  + f_0 * fk_26[k];

        t_27[k] = -2.0 * pk_27[k]
                  + f_0 * fk_27[k];

        t_28[k] = -2.0 * pk_28[k]
                  + f_0 * fk_28[k];

        t_29[k] = -2.0 * pk_29[k]
                  + f_0 * fk_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pk_30, pk_31, pk_32, pk_33, pk_34, \
                         fk_30, fk_31, fk_32, fk_33, fk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -2.0 * pk_30[k]
                  + f_0 * fk_30[k];

        t_31[k] = -2.0 * pk_31[k]
                  + f_0 * fk_31[k];

        t_32[k] = -2.0 * pk_32[k]
                  + f_0 * fk_32[k];

        t_33[k] = -2.0 * pk_33[k]
                  + f_0 * fk_33[k];

        t_34[k] = -2.0 * pk_34[k]
                  + f_0 * fk_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pk_35, pk_36, pk_37, pk_38, pk_39, \
                         fk_35, fk_36, fk_37, fk_38, fk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -2.0 * pk_35[k]
                  + f_0 * fk_35[k];

        t_36[k] = -pk_36[k]
                  + f_0 * fk_36[k];

        t_37[k] = -pk_37[k]
                  + f_0 * fk_37[k];

        t_38[k] = -pk_38[k]
                  + f_0 * fk_38[k];

        t_39[k] = -pk_39[k]
                  + f_0 * fk_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pk_40, pk_41, pk_42, pk_43, pk_44, \
                         fk_40, fk_41, fk_42, fk_43, fk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -pk_40[k]
                  + f_0 * fk_40[k];

        t_41[k] = -pk_41[k]
                  + f_0 * fk_41[k];

        t_42[k] = -pk_42[k]
                  + f_0 * fk_42[k];

        t_43[k] = -pk_43[k]
                  + f_0 * fk_43[k];

        t_44[k] = -pk_44[k]
                  + f_0 * fk_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pk_45, pk_46, pk_47, pk_48, pk_49, \
                         fk_45, fk_46, fk_47, fk_48, fk_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -pk_45[k]
                  + f_0 * fk_45[k];

        t_46[k] = -pk_46[k]
                  + f_0 * fk_46[k];

        t_47[k] = -pk_47[k]
                  + f_0 * fk_47[k];

        t_48[k] = -pk_48[k]
                  + f_0 * fk_48[k];

        t_49[k] = -pk_49[k]
                  + f_0 * fk_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pk_50, pk_51, pk_52, pk_53, pk_54, \
                         fk_50, fk_51, fk_52, fk_53, fk_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -pk_50[k]
                  + f_0 * fk_50[k];

        t_51[k] = -pk_51[k]
                  + f_0 * fk_51[k];

        t_52[k] = -pk_52[k]
                  + f_0 * fk_52[k];

        t_53[k] = -pk_53[k]
                  + f_0 * fk_53[k];

        t_54[k] = -pk_54[k]
                  + f_0 * fk_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pk_55, pk_56, pk_57, pk_58, pk_59, \
                         fk_55, fk_56, fk_57, fk_58, fk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -pk_55[k]
                  + f_0 * fk_55[k];

        t_56[k] = -pk_56[k]
                  + f_0 * fk_56[k];

        t_57[k] = -pk_57[k]
                  + f_0 * fk_57[k];

        t_58[k] = -pk_58[k]
                  + f_0 * fk_58[k];

        t_59[k] = -pk_59[k]
                  + f_0 * fk_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pk_60, pk_61, pk_62, pk_63, pk_64, \
                         fk_60, fk_61, fk_62, fk_63, fk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -pk_60[k]
                  + f_0 * fk_60[k];

        t_61[k] = -pk_61[k]
                  + f_0 * fk_61[k];

        t_62[k] = -pk_62[k]
                  + f_0 * fk_62[k];

        t_63[k] = -pk_63[k]
                  + f_0 * fk_63[k];

        t_64[k] = -pk_64[k]
                  + f_0 * fk_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pk_65, pk_66, pk_67, pk_68, pk_69, \
                         fk_65, fk_66, fk_67, fk_68, fk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -pk_65[k]
                  + f_0 * fk_65[k];

        t_66[k] = -pk_66[k]
                  + f_0 * fk_66[k];

        t_67[k] = -pk_67[k]
                  + f_0 * fk_67[k];

        t_68[k] = -pk_68[k]
                  + f_0 * fk_68[k];

        t_69[k] = -pk_69[k]
                  + f_0 * fk_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pk_70, pk_71, pk_72, pk_73, pk_74, \
                         fk_70, fk_71, fk_72, fk_73, fk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -pk_70[k]
                  + f_0 * fk_70[k];

        t_71[k] = -pk_71[k]
                  + f_0 * fk_71[k];

        t_72[k] = -pk_72[k]
                  + f_0 * fk_72[k];

        t_73[k] = -pk_73[k]
                  + f_0 * fk_73[k];

        t_74[k] = -pk_74[k]
                  + f_0 * fk_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pk_75, pk_76, pk_77, pk_78, pk_79, \
                         fk_75, fk_76, fk_77, fk_78, fk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -pk_75[k]
                  + f_0 * fk_75[k];

        t_76[k] = -pk_76[k]
                  + f_0 * fk_76[k];

        t_77[k] = -pk_77[k]
                  + f_0 * fk_77[k];

        t_78[k] = -pk_78[k]
                  + f_0 * fk_78[k];

        t_79[k] = -pk_79[k]
                  + f_0 * fk_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pk_80, pk_81, pk_82, pk_83, pk_84, \
                         fk_80, fk_81, fk_82, fk_83, fk_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -pk_80[k]
                  + f_0 * fk_80[k];

        t_81[k] = -pk_81[k]
                  + f_0 * fk_81[k];

        t_82[k] = -pk_82[k]
                  + f_0 * fk_82[k];

        t_83[k] = -pk_83[k]
                  + f_0 * fk_83[k];

        t_84[k] = -pk_84[k]
                  + f_0 * fk_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pk_85, pk_86, pk_87, pk_88, pk_89, \
                         fk_85, fk_86, fk_87, fk_88, fk_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -pk_85[k]
                  + f_0 * fk_85[k];

        t_86[k] = -pk_86[k]
                  + f_0 * fk_86[k];

        t_87[k] = -pk_87[k]
                  + f_0 * fk_87[k];

        t_88[k] = -pk_88[k]
                  + f_0 * fk_88[k];

        t_89[k] = -pk_89[k]
                  + f_0 * fk_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pk_90, pk_91, pk_92, pk_93, pk_94, \
                         fk_90, fk_91, fk_92, fk_93, fk_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -pk_90[k]
                  + f_0 * fk_90[k];

        t_91[k] = -pk_91[k]
                  + f_0 * fk_91[k];

        t_92[k] = -pk_92[k]
                  + f_0 * fk_92[k];

        t_93[k] = -pk_93[k]
                  + f_0 * fk_93[k];

        t_94[k] = -pk_94[k]
                  + f_0 * fk_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pk_95, pk_96, pk_97, pk_98, pk_99, \
                         fk_95, fk_96, fk_97, fk_98, fk_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -pk_95[k]
                  + f_0 * fk_95[k];

        t_96[k] = -pk_96[k]
                  + f_0 * fk_96[k];

        t_97[k] = -pk_97[k]
                  + f_0 * fk_97[k];

        t_98[k] = -pk_98[k]
                  + f_0 * fk_98[k];

        t_99[k] = -pk_99[k]
                  + f_0 * fk_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pk_100, pk_101, pk_102, pk_103, \
                         pk_104, fk_100, fk_101, fk_102, fk_103, \
                         fk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -pk_100[k]
                   + f_0 * fk_100[k];

        t_101[k] = -pk_101[k]
                   + f_0 * fk_101[k];

        t_102[k] = -pk_102[k]
                   + f_0 * fk_102[k];

        t_103[k] = -pk_103[k]
                   + f_0 * fk_103[k];

        t_104[k] = -pk_104[k]
                   + f_0 * fk_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, t_110, pk_105, pk_106, pk_107, \
                         fk_105, fk_106, fk_107, fk_108, fk_109, \
                         fk_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -pk_105[k]
                   + f_0 * fk_105[k];

        t_106[k] = -pk_106[k]
                   + f_0 * fk_106[k];

        t_107[k] = -pk_107[k]
                   + f_0 * fk_107[k];

        t_108[k] = f_0 * fk_108[k];

        t_109[k] = f_0 * fk_109[k];

        t_110[k] = f_0 * fk_110[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, t_116, t_117, t_118, fk_111, \
                         fk_112, fk_113, fk_114, fk_115, fk_116, fk_117, \
                         fk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_0 * fk_111[k];

        t_112[k] = f_0 * fk_112[k];

        t_113[k] = f_0 * fk_113[k];

        t_114[k] = f_0 * fk_114[k];

        t_115[k] = f_0 * fk_115[k];

        t_116[k] = f_0 * fk_116[k];

        t_117[k] = f_0 * fk_117[k];

        t_118[k] = f_0 * fk_118[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, t_124, t_125, t_126, fk_119, \
                         fk_120, fk_121, fk_122, fk_123, fk_124, fk_125, \
                         fk_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_0 * fk_119[k];

        t_120[k] = f_0 * fk_120[k];

        t_121[k] = f_0 * fk_121[k];

        t_122[k] = f_0 * fk_122[k];

        t_123[k] = f_0 * fk_123[k];

        t_124[k] = f_0 * fk_124[k];

        t_125[k] = f_0 * fk_125[k];

        t_126[k] = f_0 * fk_126[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, t_132, t_133, t_134, fk_127, \
                         fk_128, fk_129, fk_130, fk_131, fk_132, fk_133, \
                         fk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_0 * fk_127[k];

        t_128[k] = f_0 * fk_128[k];

        t_129[k] = f_0 * fk_129[k];

        t_130[k] = f_0 * fk_130[k];

        t_131[k] = f_0 * fk_131[k];

        t_132[k] = f_0 * fk_132[k];

        t_133[k] = f_0 * fk_133[k];

        t_134[k] = f_0 * fk_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, t_140, t_141, t_142, fk_135, \
                         fk_136, fk_137, fk_138, fk_139, fk_140, fk_141, \
                         fk_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * fk_135[k];

        t_136[k] = f_0 * fk_136[k];

        t_137[k] = f_0 * fk_137[k];

        t_138[k] = f_0 * fk_138[k];

        t_139[k] = f_0 * fk_139[k];

        t_140[k] = f_0 * fk_140[k];

        t_141[k] = f_0 * fk_141[k];

        t_142[k] = f_0 * fk_142[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, t_148, t_149, t_150, fk_143, \
                         fk_144, fk_145, fk_146, fk_147, fk_148, fk_149, \
                         fk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_0 * fk_143[k];

        t_144[k] = f_0 * fk_144[k];

        t_145[k] = f_0 * fk_145[k];

        t_146[k] = f_0 * fk_146[k];

        t_147[k] = f_0 * fk_147[k];

        t_148[k] = f_0 * fk_148[k];

        t_149[k] = f_0 * fk_149[k];

        t_150[k] = f_0 * fk_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, t_156, t_157, t_158, fk_151, \
                         fk_152, fk_153, fk_154, fk_155, fk_156, fk_157, \
                         fk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_0 * fk_151[k];

        t_152[k] = f_0 * fk_152[k];

        t_153[k] = f_0 * fk_153[k];

        t_154[k] = f_0 * fk_154[k];

        t_155[k] = f_0 * fk_155[k];

        t_156[k] = f_0 * fk_156[k];

        t_157[k] = f_0 * fk_157[k];

        t_158[k] = f_0 * fk_158[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, t_164, t_165, t_166, fk_159, \
                         fk_160, fk_161, fk_162, fk_163, fk_164, fk_165, \
                         fk_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_0 * fk_159[k];

        t_160[k] = f_0 * fk_160[k];

        t_161[k] = f_0 * fk_161[k];

        t_162[k] = f_0 * fk_162[k];

        t_163[k] = f_0 * fk_163[k];

        t_164[k] = f_0 * fk_164[k];

        t_165[k] = f_0 * fk_165[k];

        t_166[k] = f_0 * fk_166[k];
    }
}

static auto
compute_prim_geom_10_dk_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *fk_167 = buffer.data(fk + 167);
    const auto *fk_168 = buffer.data(fk + 168);
    const auto *fk_169 = buffer.data(fk + 169);
    const auto *fk_170 = buffer.data(fk + 170);
    const auto *fk_171 = buffer.data(fk + 171);
    const auto *fk_172 = buffer.data(fk + 172);
    const auto *fk_173 = buffer.data(fk + 173);
    const auto *fk_174 = buffer.data(fk + 174);
    const auto *fk_175 = buffer.data(fk + 175);
    const auto *fk_176 = buffer.data(fk + 176);
    const auto *fk_177 = buffer.data(fk + 177);
    const auto *fk_178 = buffer.data(fk + 178);
    const auto *fk_179 = buffer.data(fk + 179);
    const auto *fk_180 = buffer.data(fk + 180);
    const auto *fk_181 = buffer.data(fk + 181);
    const auto *fk_182 = buffer.data(fk + 182);
    const auto *fk_183 = buffer.data(fk + 183);
    const auto *fk_184 = buffer.data(fk + 184);
    const auto *fk_185 = buffer.data(fk + 185);
    const auto *fk_186 = buffer.data(fk + 186);
    const auto *fk_187 = buffer.data(fk + 187);
    const auto *fk_188 = buffer.data(fk + 188);
    const auto *fk_189 = buffer.data(fk + 189);
    const auto *fk_190 = buffer.data(fk + 190);
    const auto *fk_191 = buffer.data(fk + 191);
    const auto *fk_192 = buffer.data(fk + 192);
    const auto *fk_193 = buffer.data(fk + 193);
    const auto *fk_194 = buffer.data(fk + 194);
    const auto *fk_195 = buffer.data(fk + 195);
    const auto *fk_196 = buffer.data(fk + 196);
    const auto *fk_197 = buffer.data(fk + 197);
    const auto *fk_198 = buffer.data(fk + 198);
    const auto *fk_199 = buffer.data(fk + 199);
    const auto *fk_200 = buffer.data(fk + 200);
    const auto *fk_201 = buffer.data(fk + 201);
    const auto *fk_202 = buffer.data(fk + 202);
    const auto *fk_203 = buffer.data(fk + 203);
    const auto *fk_204 = buffer.data(fk + 204);
    const auto *fk_205 = buffer.data(fk + 205);
    const auto *fk_206 = buffer.data(fk + 206);
    const auto *fk_207 = buffer.data(fk + 207);
    const auto *fk_208 = buffer.data(fk + 208);
    const auto *fk_209 = buffer.data(fk + 209);
    const auto *fk_210 = buffer.data(fk + 210);
    const auto *fk_211 = buffer.data(fk + 211);
    const auto *fk_212 = buffer.data(fk + 212);
    const auto *fk_213 = buffer.data(fk + 213);
    const auto *fk_214 = buffer.data(fk + 214);
    const auto *fk_215 = buffer.data(fk + 215);

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, t_172, t_173, t_174, fk_167, \
                         fk_168, fk_169, fk_170, fk_171, fk_172, fk_173, \
                         fk_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_0 * fk_167[k];

        t_168[k] = f_0 * fk_168[k];

        t_169[k] = f_0 * fk_169[k];

        t_170[k] = f_0 * fk_170[k];

        t_171[k] = f_0 * fk_171[k];

        t_172[k] = f_0 * fk_172[k];

        t_173[k] = f_0 * fk_173[k];

        t_174[k] = f_0 * fk_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, t_180, t_181, t_182, fk_175, \
                         fk_176, fk_177, fk_178, fk_179, fk_180, fk_181, \
                         fk_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_0 * fk_175[k];

        t_176[k] = f_0 * fk_176[k];

        t_177[k] = f_0 * fk_177[k];

        t_178[k] = f_0 * fk_178[k];

        t_179[k] = f_0 * fk_179[k];

        t_180[k] = f_0 * fk_180[k];

        t_181[k] = f_0 * fk_181[k];

        t_182[k] = f_0 * fk_182[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, t_188, t_189, t_190, fk_183, \
                         fk_184, fk_185, fk_186, fk_187, fk_188, fk_189, \
                         fk_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_0 * fk_183[k];

        t_184[k] = f_0 * fk_184[k];

        t_185[k] = f_0 * fk_185[k];

        t_186[k] = f_0 * fk_186[k];

        t_187[k] = f_0 * fk_187[k];

        t_188[k] = f_0 * fk_188[k];

        t_189[k] = f_0 * fk_189[k];

        t_190[k] = f_0 * fk_190[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, t_197, t_198, fk_191, \
                         fk_192, fk_193, fk_194, fk_195, fk_196, fk_197, \
                         fk_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_0 * fk_191[k];

        t_192[k] = f_0 * fk_192[k];

        t_193[k] = f_0 * fk_193[k];

        t_194[k] = f_0 * fk_194[k];

        t_195[k] = f_0 * fk_195[k];

        t_196[k] = f_0 * fk_196[k];

        t_197[k] = f_0 * fk_197[k];

        t_198[k] = f_0 * fk_198[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, t_204, t_205, t_206, fk_199, \
                         fk_200, fk_201, fk_202, fk_203, fk_204, fk_205, \
                         fk_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_0 * fk_199[k];

        t_200[k] = f_0 * fk_200[k];

        t_201[k] = f_0 * fk_201[k];

        t_202[k] = f_0 * fk_202[k];

        t_203[k] = f_0 * fk_203[k];

        t_204[k] = f_0 * fk_204[k];

        t_205[k] = f_0 * fk_205[k];

        t_206[k] = f_0 * fk_206[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, t_213, t_214, fk_207, \
                         fk_208, fk_209, fk_210, fk_211, fk_212, fk_213, \
                         fk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_0 * fk_207[k];

        t_208[k] = f_0 * fk_208[k];

        t_209[k] = f_0 * fk_209[k];

        t_210[k] = f_0 * fk_210[k];

        t_211[k] = f_0 * fk_211[k];

        t_212[k] = f_0 * fk_212[k];

        t_213[k] = f_0 * fk_213[k];

        t_214[k] = f_0 * fk_214[k];
    }

#pragma omp simd aligned(t_215, fk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_0 * fk_215[k];
    }
}

auto
compute_prim_geom_10_dk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t pk, const size_t fk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_dk_electron_repulsion_0_piece0(buffer, target, pk, fk, ncols, alpha);

    compute_prim_geom_10_dk_electron_repulsion_0_piece1(buffer, target, fk, ncols, alpha);
}

static auto
compute_prim_geom_10_dk_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t pk, const size_t fk,
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

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_11 = buffer.data(pk + 11);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_16 = buffer.data(pk + 16);
    const auto *pk_17 = buffer.data(pk + 17);
    const auto *pk_18 = buffer.data(pk + 18);
    const auto *pk_19 = buffer.data(pk + 19);
    const auto *pk_20 = buffer.data(pk + 20);
    const auto *pk_21 = buffer.data(pk + 21);
    const auto *pk_22 = buffer.data(pk + 22);
    const auto *pk_23 = buffer.data(pk + 23);
    const auto *pk_24 = buffer.data(pk + 24);
    const auto *pk_25 = buffer.data(pk + 25);
    const auto *pk_26 = buffer.data(pk + 26);
    const auto *pk_27 = buffer.data(pk + 27);
    const auto *pk_28 = buffer.data(pk + 28);
    const auto *pk_29 = buffer.data(pk + 29);
    const auto *pk_30 = buffer.data(pk + 30);
    const auto *pk_31 = buffer.data(pk + 31);
    const auto *pk_32 = buffer.data(pk + 32);
    const auto *pk_33 = buffer.data(pk + 33);
    const auto *pk_34 = buffer.data(pk + 34);
    const auto *pk_35 = buffer.data(pk + 35);
    const auto *pk_36 = buffer.data(pk + 36);
    const auto *pk_37 = buffer.data(pk + 37);
    const auto *pk_38 = buffer.data(pk + 38);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_40 = buffer.data(pk + 40);
    const auto *pk_41 = buffer.data(pk + 41);
    const auto *pk_42 = buffer.data(pk + 42);
    const auto *pk_43 = buffer.data(pk + 43);
    const auto *pk_44 = buffer.data(pk + 44);
    const auto *pk_45 = buffer.data(pk + 45);
    const auto *pk_46 = buffer.data(pk + 46);
    const auto *pk_47 = buffer.data(pk + 47);
    const auto *pk_48 = buffer.data(pk + 48);
    const auto *pk_49 = buffer.data(pk + 49);
    const auto *pk_50 = buffer.data(pk + 50);
    const auto *pk_51 = buffer.data(pk + 51);
    const auto *pk_52 = buffer.data(pk + 52);
    const auto *pk_53 = buffer.data(pk + 53);
    const auto *pk_54 = buffer.data(pk + 54);
    const auto *pk_55 = buffer.data(pk + 55);
    const auto *pk_56 = buffer.data(pk + 56);
    const auto *pk_57 = buffer.data(pk + 57);
    const auto *pk_58 = buffer.data(pk + 58);
    const auto *pk_59 = buffer.data(pk + 59);
    const auto *pk_60 = buffer.data(pk + 60);
    const auto *pk_61 = buffer.data(pk + 61);
    const auto *pk_62 = buffer.data(pk + 62);
    const auto *pk_63 = buffer.data(pk + 63);
    const auto *pk_64 = buffer.data(pk + 64);
    const auto *pk_65 = buffer.data(pk + 65);
    const auto *pk_66 = buffer.data(pk + 66);
    const auto *pk_67 = buffer.data(pk + 67);
    const auto *pk_68 = buffer.data(pk + 68);
    const auto *pk_69 = buffer.data(pk + 69);
    const auto *pk_70 = buffer.data(pk + 70);
    const auto *pk_71 = buffer.data(pk + 71);
    const auto *pk_72 = buffer.data(pk + 72);
    const auto *pk_73 = buffer.data(pk + 73);
    const auto *pk_74 = buffer.data(pk + 74);
    const auto *pk_75 = buffer.data(pk + 75);
    const auto *pk_76 = buffer.data(pk + 76);
    const auto *pk_77 = buffer.data(pk + 77);
    const auto *pk_78 = buffer.data(pk + 78);
    const auto *pk_79 = buffer.data(pk + 79);
    const auto *pk_80 = buffer.data(pk + 80);
    const auto *pk_81 = buffer.data(pk + 81);
    const auto *pk_82 = buffer.data(pk + 82);
    const auto *pk_83 = buffer.data(pk + 83);
    const auto *pk_84 = buffer.data(pk + 84);
    const auto *pk_85 = buffer.data(pk + 85);
    const auto *pk_86 = buffer.data(pk + 86);
    const auto *pk_87 = buffer.data(pk + 87);
    const auto *pk_88 = buffer.data(pk + 88);
    const auto *pk_89 = buffer.data(pk + 89);
    const auto *pk_90 = buffer.data(pk + 90);
    const auto *pk_91 = buffer.data(pk + 91);
    const auto *pk_92 = buffer.data(pk + 92);
    const auto *pk_93 = buffer.data(pk + 93);
    const auto *pk_94 = buffer.data(pk + 94);
    const auto *pk_95 = buffer.data(pk + 95);
    const auto *pk_96 = buffer.data(pk + 96);
    const auto *pk_97 = buffer.data(pk + 97);
    const auto *pk_98 = buffer.data(pk + 98);
    const auto *pk_99 = buffer.data(pk + 99);

    const auto *fk_36 = buffer.data(fk + 36);
    const auto *fk_37 = buffer.data(fk + 37);
    const auto *fk_38 = buffer.data(fk + 38);
    const auto *fk_39 = buffer.data(fk + 39);
    const auto *fk_40 = buffer.data(fk + 40);
    const auto *fk_41 = buffer.data(fk + 41);
    const auto *fk_42 = buffer.data(fk + 42);
    const auto *fk_43 = buffer.data(fk + 43);
    const auto *fk_44 = buffer.data(fk + 44);
    const auto *fk_45 = buffer.data(fk + 45);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_47 = buffer.data(fk + 47);
    const auto *fk_48 = buffer.data(fk + 48);
    const auto *fk_49 = buffer.data(fk + 49);
    const auto *fk_50 = buffer.data(fk + 50);
    const auto *fk_51 = buffer.data(fk + 51);
    const auto *fk_52 = buffer.data(fk + 52);
    const auto *fk_53 = buffer.data(fk + 53);
    const auto *fk_54 = buffer.data(fk + 54);
    const auto *fk_55 = buffer.data(fk + 55);
    const auto *fk_56 = buffer.data(fk + 56);
    const auto *fk_57 = buffer.data(fk + 57);
    const auto *fk_58 = buffer.data(fk + 58);
    const auto *fk_59 = buffer.data(fk + 59);
    const auto *fk_60 = buffer.data(fk + 60);
    const auto *fk_61 = buffer.data(fk + 61);
    const auto *fk_62 = buffer.data(fk + 62);
    const auto *fk_63 = buffer.data(fk + 63);
    const auto *fk_64 = buffer.data(fk + 64);
    const auto *fk_65 = buffer.data(fk + 65);
    const auto *fk_66 = buffer.data(fk + 66);
    const auto *fk_67 = buffer.data(fk + 67);
    const auto *fk_68 = buffer.data(fk + 68);
    const auto *fk_69 = buffer.data(fk + 69);
    const auto *fk_70 = buffer.data(fk + 70);
    const auto *fk_71 = buffer.data(fk + 71);
    const auto *fk_108 = buffer.data(fk + 108);
    const auto *fk_109 = buffer.data(fk + 109);
    const auto *fk_110 = buffer.data(fk + 110);
    const auto *fk_111 = buffer.data(fk + 111);
    const auto *fk_112 = buffer.data(fk + 112);
    const auto *fk_113 = buffer.data(fk + 113);
    const auto *fk_114 = buffer.data(fk + 114);
    const auto *fk_115 = buffer.data(fk + 115);
    const auto *fk_116 = buffer.data(fk + 116);
    const auto *fk_117 = buffer.data(fk + 117);
    const auto *fk_118 = buffer.data(fk + 118);
    const auto *fk_119 = buffer.data(fk + 119);
    const auto *fk_120 = buffer.data(fk + 120);
    const auto *fk_121 = buffer.data(fk + 121);
    const auto *fk_122 = buffer.data(fk + 122);
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_124 = buffer.data(fk + 124);
    const auto *fk_125 = buffer.data(fk + 125);
    const auto *fk_126 = buffer.data(fk + 126);
    const auto *fk_127 = buffer.data(fk + 127);
    const auto *fk_128 = buffer.data(fk + 128);
    const auto *fk_129 = buffer.data(fk + 129);
    const auto *fk_130 = buffer.data(fk + 130);
    const auto *fk_131 = buffer.data(fk + 131);
    const auto *fk_132 = buffer.data(fk + 132);
    const auto *fk_133 = buffer.data(fk + 133);
    const auto *fk_134 = buffer.data(fk + 134);
    const auto *fk_135 = buffer.data(fk + 135);
    const auto *fk_136 = buffer.data(fk + 136);
    const auto *fk_137 = buffer.data(fk + 137);
    const auto *fk_138 = buffer.data(fk + 138);
    const auto *fk_139 = buffer.data(fk + 139);
    const auto *fk_140 = buffer.data(fk + 140);
    const auto *fk_141 = buffer.data(fk + 141);
    const auto *fk_142 = buffer.data(fk + 142);
    const auto *fk_143 = buffer.data(fk + 143);
    const auto *fk_144 = buffer.data(fk + 144);
    const auto *fk_145 = buffer.data(fk + 145);
    const auto *fk_146 = buffer.data(fk + 146);
    const auto *fk_147 = buffer.data(fk + 147);
    const auto *fk_148 = buffer.data(fk + 148);
    const auto *fk_149 = buffer.data(fk + 149);
    const auto *fk_150 = buffer.data(fk + 150);
    const auto *fk_151 = buffer.data(fk + 151);
    const auto *fk_152 = buffer.data(fk + 152);
    const auto *fk_153 = buffer.data(fk + 153);
    const auto *fk_154 = buffer.data(fk + 154);
    const auto *fk_155 = buffer.data(fk + 155);
    const auto *fk_156 = buffer.data(fk + 156);
    const auto *fk_157 = buffer.data(fk + 157);
    const auto *fk_158 = buffer.data(fk + 158);
    const auto *fk_159 = buffer.data(fk + 159);
    const auto *fk_160 = buffer.data(fk + 160);
    const auto *fk_161 = buffer.data(fk + 161);
    const auto *fk_162 = buffer.data(fk + 162);
    const auto *fk_163 = buffer.data(fk + 163);
    const auto *fk_164 = buffer.data(fk + 164);
    const auto *fk_165 = buffer.data(fk + 165);
    const auto *fk_166 = buffer.data(fk + 166);
    const auto *fk_167 = buffer.data(fk + 167);
    const auto *fk_168 = buffer.data(fk + 168);
    const auto *fk_169 = buffer.data(fk + 169);
    const auto *fk_170 = buffer.data(fk + 170);
    const auto *fk_171 = buffer.data(fk + 171);
    const auto *fk_172 = buffer.data(fk + 172);
    const auto *fk_173 = buffer.data(fk + 173);
    const auto *fk_174 = buffer.data(fk + 174);
    const auto *fk_175 = buffer.data(fk + 175);
    const auto *fk_176 = buffer.data(fk + 176);
    const auto *fk_177 = buffer.data(fk + 177);
    const auto *fk_178 = buffer.data(fk + 178);
    const auto *fk_179 = buffer.data(fk + 179);
    const auto *fk_216 = buffer.data(fk + 216);
    const auto *fk_217 = buffer.data(fk + 217);
    const auto *fk_218 = buffer.data(fk + 218);
    const auto *fk_219 = buffer.data(fk + 219);
    const auto *fk_220 = buffer.data(fk + 220);
    const auto *fk_221 = buffer.data(fk + 221);
    const auto *fk_222 = buffer.data(fk + 222);
    const auto *fk_223 = buffer.data(fk + 223);
    const auto *fk_224 = buffer.data(fk + 224);
    const auto *fk_225 = buffer.data(fk + 225);
    const auto *fk_226 = buffer.data(fk + 226);
    const auto *fk_227 = buffer.data(fk + 227);
    const auto *fk_228 = buffer.data(fk + 228);
    const auto *fk_229 = buffer.data(fk + 229);
    const auto *fk_230 = buffer.data(fk + 230);
    const auto *fk_231 = buffer.data(fk + 231);
    const auto *fk_232 = buffer.data(fk + 232);
    const auto *fk_233 = buffer.data(fk + 233);
    const auto *fk_234 = buffer.data(fk + 234);
    const auto *fk_235 = buffer.data(fk + 235);
    const auto *fk_236 = buffer.data(fk + 236);
    const auto *fk_237 = buffer.data(fk + 237);
    const auto *fk_238 = buffer.data(fk + 238);
    const auto *fk_239 = buffer.data(fk + 239);
    const auto *fk_240 = buffer.data(fk + 240);
    const auto *fk_241 = buffer.data(fk + 241);
    const auto *fk_242 = buffer.data(fk + 242);
    const auto *fk_243 = buffer.data(fk + 243);
    const auto *fk_244 = buffer.data(fk + 244);
    const auto *fk_245 = buffer.data(fk + 245);
    const auto *fk_246 = buffer.data(fk + 246);
    const auto *fk_247 = buffer.data(fk + 247);
    const auto *fk_248 = buffer.data(fk + 248);
    const auto *fk_249 = buffer.data(fk + 249);
    const auto *fk_250 = buffer.data(fk + 250);
    const auto *fk_251 = buffer.data(fk + 251);
    const auto *fk_252 = buffer.data(fk + 252);
    const auto *fk_253 = buffer.data(fk + 253);
    const auto *fk_254 = buffer.data(fk + 254);
    const auto *fk_255 = buffer.data(fk + 255);
    const auto *fk_256 = buffer.data(fk + 256);
    const auto *fk_257 = buffer.data(fk + 257);
    const auto *fk_258 = buffer.data(fk + 258);
    const auto *fk_259 = buffer.data(fk + 259);
    const auto *fk_260 = buffer.data(fk + 260);
    const auto *fk_261 = buffer.data(fk + 261);
    const auto *fk_262 = buffer.data(fk + 262);
    const auto *fk_263 = buffer.data(fk + 263);
    const auto *fk_264 = buffer.data(fk + 264);
    const auto *fk_265 = buffer.data(fk + 265);
    const auto *fk_266 = buffer.data(fk + 266);
    const auto *fk_267 = buffer.data(fk + 267);
    const auto *fk_268 = buffer.data(fk + 268);
    const auto *fk_269 = buffer.data(fk + 269);
    const auto *fk_270 = buffer.data(fk + 270);
    const auto *fk_271 = buffer.data(fk + 271);
    const auto *fk_272 = buffer.data(fk + 272);
    const auto *fk_273 = buffer.data(fk + 273);
    const auto *fk_274 = buffer.data(fk + 274);
    const auto *fk_275 = buffer.data(fk + 275);
    const auto *fk_276 = buffer.data(fk + 276);
    const auto *fk_277 = buffer.data(fk + 277);
    const auto *fk_278 = buffer.data(fk + 278);
    const auto *fk_279 = buffer.data(fk + 279);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, fk_36, fk_37, fk_38, fk_39, \
                         fk_40, fk_41, fk_42, fk_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_36[k];

        t_1[k] = f_0 * fk_37[k];

        t_2[k] = f_0 * fk_38[k];

        t_3[k] = f_0 * fk_39[k];

        t_4[k] = f_0 * fk_40[k];

        t_5[k] = f_0 * fk_41[k];

        t_6[k] = f_0 * fk_42[k];

        t_7[k] = f_0 * fk_43[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, fk_44, fk_45, fk_46, \
                         fk_47, fk_48, fk_49, fk_50, fk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * fk_44[k];

        t_9[k] = f_0 * fk_45[k];

        t_10[k] = f_0 * fk_46[k];

        t_11[k] = f_0 * fk_47[k];

        t_12[k] = f_0 * fk_48[k];

        t_13[k] = f_0 * fk_49[k];

        t_14[k] = f_0 * fk_50[k];

        t_15[k] = f_0 * fk_51[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, fk_52, fk_53, fk_54, \
                         fk_55, fk_56, fk_57, fk_58, fk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * fk_52[k];

        t_17[k] = f_0 * fk_53[k];

        t_18[k] = f_0 * fk_54[k];

        t_19[k] = f_0 * fk_55[k];

        t_20[k] = f_0 * fk_56[k];

        t_21[k] = f_0 * fk_57[k];

        t_22[k] = f_0 * fk_58[k];

        t_23[k] = f_0 * fk_59[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, fk_60, fk_61, fk_62, \
                         fk_63, fk_64, fk_65, fk_66, fk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * fk_60[k];

        t_25[k] = f_0 * fk_61[k];

        t_26[k] = f_0 * fk_62[k];

        t_27[k] = f_0 * fk_63[k];

        t_28[k] = f_0 * fk_64[k];

        t_29[k] = f_0 * fk_65[k];

        t_30[k] = f_0 * fk_66[k];

        t_31[k] = f_0 * fk_67[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, pk_0, pk_1, fk_68, fk_69, fk_70, \
                         fk_71, fk_108, fk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * fk_68[k];

        t_33[k] = f_0 * fk_69[k];

        t_34[k] = f_0 * fk_70[k];

        t_35[k] = f_0 * fk_71[k];

        t_36[k] = -pk_0[k]
                  + f_0 * fk_108[k];

        t_37[k] = -pk_1[k]
                  + f_0 * fk_109[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pk_2, pk_3, pk_4, pk_5, pk_6, fk_110, \
                         fk_111, fk_112, fk_113, fk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -pk_2[k]
                  + f_0 * fk_110[k];

        t_39[k] = -pk_3[k]
                  + f_0 * fk_111[k];

        t_40[k] = -pk_4[k]
                  + f_0 * fk_112[k];

        t_41[k] = -pk_5[k]
                  + f_0 * fk_113[k];

        t_42[k] = -pk_6[k]
                  + f_0 * fk_114[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pk_7, pk_8, pk_9, pk_10, pk_11, fk_115, \
                         fk_116, fk_117, fk_118, fk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -pk_7[k]
                  + f_0 * fk_115[k];

        t_44[k] = -pk_8[k]
                  + f_0 * fk_116[k];

        t_45[k] = -pk_9[k]
                  + f_0 * fk_117[k];

        t_46[k] = -pk_10[k]
                  + f_0 * fk_118[k];

        t_47[k] = -pk_11[k]
                  + f_0 * fk_119[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pk_12, pk_13, pk_14, pk_15, pk_16, \
                         fk_120, fk_121, fk_122, fk_123, fk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -pk_12[k]
                  + f_0 * fk_120[k];

        t_49[k] = -pk_13[k]
                  + f_0 * fk_121[k];

        t_50[k] = -pk_14[k]
                  + f_0 * fk_122[k];

        t_51[k] = -pk_15[k]
                  + f_0 * fk_123[k];

        t_52[k] = -pk_16[k]
                  + f_0 * fk_124[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pk_17, pk_18, pk_19, pk_20, pk_21, \
                         fk_125, fk_126, fk_127, fk_128, fk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -pk_17[k]
                  + f_0 * fk_125[k];

        t_54[k] = -pk_18[k]
                  + f_0 * fk_126[k];

        t_55[k] = -pk_19[k]
                  + f_0 * fk_127[k];

        t_56[k] = -pk_20[k]
                  + f_0 * fk_128[k];

        t_57[k] = -pk_21[k]
                  + f_0 * fk_129[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pk_22, pk_23, pk_24, pk_25, pk_26, \
                         fk_130, fk_131, fk_132, fk_133, fk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -pk_22[k]
                  + f_0 * fk_130[k];

        t_59[k] = -pk_23[k]
                  + f_0 * fk_131[k];

        t_60[k] = -pk_24[k]
                  + f_0 * fk_132[k];

        t_61[k] = -pk_25[k]
                  + f_0 * fk_133[k];

        t_62[k] = -pk_26[k]
                  + f_0 * fk_134[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pk_27, pk_28, pk_29, pk_30, pk_31, \
                         fk_135, fk_136, fk_137, fk_138, fk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -pk_27[k]
                  + f_0 * fk_135[k];

        t_64[k] = -pk_28[k]
                  + f_0 * fk_136[k];

        t_65[k] = -pk_29[k]
                  + f_0 * fk_137[k];

        t_66[k] = -pk_30[k]
                  + f_0 * fk_138[k];

        t_67[k] = -pk_31[k]
                  + f_0 * fk_139[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, pk_32, pk_33, pk_34, pk_35, \
                         fk_140, fk_141, fk_142, fk_143, fk_144, \
                         fk_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -pk_32[k]
                  + f_0 * fk_140[k];

        t_69[k] = -pk_33[k]
                  + f_0 * fk_141[k];

        t_70[k] = -pk_34[k]
                  + f_0 * fk_142[k];

        t_71[k] = -pk_35[k]
                  + f_0 * fk_143[k];

        t_72[k] = f_0 * fk_144[k];

        t_73[k] = f_0 * fk_145[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, t_80, t_81, fk_146, fk_147, \
                         fk_148, fk_149, fk_150, fk_151, fk_152, \
                         fk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * fk_146[k];

        t_75[k] = f_0 * fk_147[k];

        t_76[k] = f_0 * fk_148[k];

        t_77[k] = f_0 * fk_149[k];

        t_78[k] = f_0 * fk_150[k];

        t_79[k] = f_0 * fk_151[k];

        t_80[k] = f_0 * fk_152[k];

        t_81[k] = f_0 * fk_153[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, t_87, t_88, t_89, fk_154, fk_155, \
                         fk_156, fk_157, fk_158, fk_159, fk_160, \
                         fk_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * fk_154[k];

        t_83[k] = f_0 * fk_155[k];

        t_84[k] = f_0 * fk_156[k];

        t_85[k] = f_0 * fk_157[k];

        t_86[k] = f_0 * fk_158[k];

        t_87[k] = f_0 * fk_159[k];

        t_88[k] = f_0 * fk_160[k];

        t_89[k] = f_0 * fk_161[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, fk_162, fk_163, \
                         fk_164, fk_165, fk_166, fk_167, fk_168, \
                         fk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * fk_162[k];

        t_91[k] = f_0 * fk_163[k];

        t_92[k] = f_0 * fk_164[k];

        t_93[k] = f_0 * fk_165[k];

        t_94[k] = f_0 * fk_166[k];

        t_95[k] = f_0 * fk_167[k];

        t_96[k] = f_0 * fk_168[k];

        t_97[k] = f_0 * fk_169[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, t_104, t_105, fk_170, fk_171, \
                         fk_172, fk_173, fk_174, fk_175, fk_176, \
                         fk_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * fk_170[k];

        t_99[k] = f_0 * fk_171[k];

        t_100[k] = f_0 * fk_172[k];

        t_101[k] = f_0 * fk_173[k];

        t_102[k] = f_0 * fk_174[k];

        t_103[k] = f_0 * fk_175[k];

        t_104[k] = f_0 * fk_176[k];

        t_105[k] = f_0 * fk_177[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, pk_36, pk_37, pk_38, pk_39, \
                         fk_178, fk_179, fk_216, fk_217, fk_218, \
                         fk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * fk_178[k];

        t_107[k] = f_0 * fk_179[k];

        t_108[k] = -2.0 * pk_36[k]
                   + f_0 * fk_216[k];

        t_109[k] = -2.0 * pk_37[k]
                   + f_0 * fk_217[k];

        t_110[k] = -2.0 * pk_38[k]
                   + f_0 * fk_218[k];

        t_111[k] = -2.0 * pk_39[k]
                   + f_0 * fk_219[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, pk_40, pk_41, pk_42, pk_43, pk_44, \
                         fk_220, fk_221, fk_222, fk_223, fk_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -2.0 * pk_40[k]
                   + f_0 * fk_220[k];

        t_113[k] = -2.0 * pk_41[k]
                   + f_0 * fk_221[k];

        t_114[k] = -2.0 * pk_42[k]
                   + f_0 * fk_222[k];

        t_115[k] = -2.0 * pk_43[k]
                   + f_0 * fk_223[k];

        t_116[k] = -2.0 * pk_44[k]
                   + f_0 * fk_224[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, pk_45, pk_46, pk_47, pk_48, pk_49, \
                         fk_225, fk_226, fk_227, fk_228, fk_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -2.0 * pk_45[k]
                   + f_0 * fk_225[k];

        t_118[k] = -2.0 * pk_46[k]
                   + f_0 * fk_226[k];

        t_119[k] = -2.0 * pk_47[k]
                   + f_0 * fk_227[k];

        t_120[k] = -2.0 * pk_48[k]
                   + f_0 * fk_228[k];

        t_121[k] = -2.0 * pk_49[k]
                   + f_0 * fk_229[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, pk_50, pk_51, pk_52, pk_53, pk_54, \
                         fk_230, fk_231, fk_232, fk_233, fk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * pk_50[k]
                   + f_0 * fk_230[k];

        t_123[k] = -2.0 * pk_51[k]
                   + f_0 * fk_231[k];

        t_124[k] = -2.0 * pk_52[k]
                   + f_0 * fk_232[k];

        t_125[k] = -2.0 * pk_53[k]
                   + f_0 * fk_233[k];

        t_126[k] = -2.0 * pk_54[k]
                   + f_0 * fk_234[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pk_55, pk_56, pk_57, pk_58, pk_59, \
                         fk_235, fk_236, fk_237, fk_238, fk_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * pk_55[k]
                   + f_0 * fk_235[k];

        t_128[k] = -2.0 * pk_56[k]
                   + f_0 * fk_236[k];

        t_129[k] = -2.0 * pk_57[k]
                   + f_0 * fk_237[k];

        t_130[k] = -2.0 * pk_58[k]
                   + f_0 * fk_238[k];

        t_131[k] = -2.0 * pk_59[k]
                   + f_0 * fk_239[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, pk_60, pk_61, pk_62, pk_63, pk_64, \
                         fk_240, fk_241, fk_242, fk_243, fk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * pk_60[k]
                   + f_0 * fk_240[k];

        t_133[k] = -2.0 * pk_61[k]
                   + f_0 * fk_241[k];

        t_134[k] = -2.0 * pk_62[k]
                   + f_0 * fk_242[k];

        t_135[k] = -2.0 * pk_63[k]
                   + f_0 * fk_243[k];

        t_136[k] = -2.0 * pk_64[k]
                   + f_0 * fk_244[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, pk_65, pk_66, pk_67, pk_68, pk_69, \
                         fk_245, fk_246, fk_247, fk_248, fk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -2.0 * pk_65[k]
                   + f_0 * fk_245[k];

        t_138[k] = -2.0 * pk_66[k]
                   + f_0 * fk_246[k];

        t_139[k] = -2.0 * pk_67[k]
                   + f_0 * fk_247[k];

        t_140[k] = -2.0 * pk_68[k]
                   + f_0 * fk_248[k];

        t_141[k] = -2.0 * pk_69[k]
                   + f_0 * fk_249[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, pk_70, pk_71, pk_72, pk_73, pk_74, \
                         fk_250, fk_251, fk_252, fk_253, fk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -2.0 * pk_70[k]
                   + f_0 * fk_250[k];

        t_143[k] = -2.0 * pk_71[k]
                   + f_0 * fk_251[k];

        t_144[k] = -pk_72[k]
                   + f_0 * fk_252[k];

        t_145[k] = -pk_73[k]
                   + f_0 * fk_253[k];

        t_146[k] = -pk_74[k]
                   + f_0 * fk_254[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, pk_75, pk_76, pk_77, pk_78, pk_79, \
                         fk_255, fk_256, fk_257, fk_258, fk_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -pk_75[k]
                   + f_0 * fk_255[k];

        t_148[k] = -pk_76[k]
                   + f_0 * fk_256[k];

        t_149[k] = -pk_77[k]
                   + f_0 * fk_257[k];

        t_150[k] = -pk_78[k]
                   + f_0 * fk_258[k];

        t_151[k] = -pk_79[k]
                   + f_0 * fk_259[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, pk_80, pk_81, pk_82, pk_83, pk_84, \
                         fk_260, fk_261, fk_262, fk_263, fk_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -pk_80[k]
                   + f_0 * fk_260[k];

        t_153[k] = -pk_81[k]
                   + f_0 * fk_261[k];

        t_154[k] = -pk_82[k]
                   + f_0 * fk_262[k];

        t_155[k] = -pk_83[k]
                   + f_0 * fk_263[k];

        t_156[k] = -pk_84[k]
                   + f_0 * fk_264[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pk_85, pk_86, pk_87, pk_88, pk_89, \
                         fk_265, fk_266, fk_267, fk_268, fk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -pk_85[k]
                   + f_0 * fk_265[k];

        t_158[k] = -pk_86[k]
                   + f_0 * fk_266[k];

        t_159[k] = -pk_87[k]
                   + f_0 * fk_267[k];

        t_160[k] = -pk_88[k]
                   + f_0 * fk_268[k];

        t_161[k] = -pk_89[k]
                   + f_0 * fk_269[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, pk_90, pk_91, pk_92, pk_93, pk_94, \
                         fk_270, fk_271, fk_272, fk_273, fk_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -pk_90[k]
                   + f_0 * fk_270[k];

        t_163[k] = -pk_91[k]
                   + f_0 * fk_271[k];

        t_164[k] = -pk_92[k]
                   + f_0 * fk_272[k];

        t_165[k] = -pk_93[k]
                   + f_0 * fk_273[k];

        t_166[k] = -pk_94[k]
                   + f_0 * fk_274[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pk_95, pk_96, pk_97, pk_98, pk_99, \
                         fk_275, fk_276, fk_277, fk_278, fk_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -pk_95[k]
                   + f_0 * fk_275[k];

        t_168[k] = -pk_96[k]
                   + f_0 * fk_276[k];

        t_169[k] = -pk_97[k]
                   + f_0 * fk_277[k];

        t_170[k] = -pk_98[k]
                   + f_0 * fk_278[k];

        t_171[k] = -pk_99[k]
                   + f_0 * fk_279[k];
    }
}

static auto
compute_prim_geom_10_dk_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t pk, const size_t fk,
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

    const auto *pk_100 = buffer.data(pk + 100);
    const auto *pk_101 = buffer.data(pk + 101);
    const auto *pk_102 = buffer.data(pk + 102);
    const auto *pk_103 = buffer.data(pk + 103);
    const auto *pk_104 = buffer.data(pk + 104);
    const auto *pk_105 = buffer.data(pk + 105);
    const auto *pk_106 = buffer.data(pk + 106);
    const auto *pk_107 = buffer.data(pk + 107);

    const auto *fk_280 = buffer.data(fk + 280);
    const auto *fk_281 = buffer.data(fk + 281);
    const auto *fk_282 = buffer.data(fk + 282);
    const auto *fk_283 = buffer.data(fk + 283);
    const auto *fk_284 = buffer.data(fk + 284);
    const auto *fk_285 = buffer.data(fk + 285);
    const auto *fk_286 = buffer.data(fk + 286);
    const auto *fk_287 = buffer.data(fk + 287);
    const auto *fk_288 = buffer.data(fk + 288);
    const auto *fk_289 = buffer.data(fk + 289);
    const auto *fk_290 = buffer.data(fk + 290);
    const auto *fk_291 = buffer.data(fk + 291);
    const auto *fk_292 = buffer.data(fk + 292);
    const auto *fk_293 = buffer.data(fk + 293);
    const auto *fk_294 = buffer.data(fk + 294);
    const auto *fk_295 = buffer.data(fk + 295);
    const auto *fk_296 = buffer.data(fk + 296);
    const auto *fk_297 = buffer.data(fk + 297);
    const auto *fk_298 = buffer.data(fk + 298);
    const auto *fk_299 = buffer.data(fk + 299);
    const auto *fk_300 = buffer.data(fk + 300);
    const auto *fk_301 = buffer.data(fk + 301);
    const auto *fk_302 = buffer.data(fk + 302);
    const auto *fk_303 = buffer.data(fk + 303);
    const auto *fk_304 = buffer.data(fk + 304);
    const auto *fk_305 = buffer.data(fk + 305);
    const auto *fk_306 = buffer.data(fk + 306);
    const auto *fk_307 = buffer.data(fk + 307);
    const auto *fk_308 = buffer.data(fk + 308);
    const auto *fk_309 = buffer.data(fk + 309);
    const auto *fk_310 = buffer.data(fk + 310);
    const auto *fk_311 = buffer.data(fk + 311);
    const auto *fk_312 = buffer.data(fk + 312);
    const auto *fk_313 = buffer.data(fk + 313);
    const auto *fk_314 = buffer.data(fk + 314);
    const auto *fk_315 = buffer.data(fk + 315);
    const auto *fk_316 = buffer.data(fk + 316);
    const auto *fk_317 = buffer.data(fk + 317);
    const auto *fk_318 = buffer.data(fk + 318);
    const auto *fk_319 = buffer.data(fk + 319);
    const auto *fk_320 = buffer.data(fk + 320);
    const auto *fk_321 = buffer.data(fk + 321);
    const auto *fk_322 = buffer.data(fk + 322);
    const auto *fk_323 = buffer.data(fk + 323);

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, pk_100, pk_101, pk_102, pk_103, \
                         pk_104, fk_280, fk_281, fk_282, fk_283, \
                         fk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -pk_100[k]
                   + f_0 * fk_280[k];

        t_173[k] = -pk_101[k]
                   + f_0 * fk_281[k];

        t_174[k] = -pk_102[k]
                   + f_0 * fk_282[k];

        t_175[k] = -pk_103[k]
                   + f_0 * fk_283[k];

        t_176[k] = -pk_104[k]
                   + f_0 * fk_284[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, t_182, pk_105, pk_106, pk_107, \
                         fk_285, fk_286, fk_287, fk_288, fk_289, \
                         fk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -pk_105[k]
                   + f_0 * fk_285[k];

        t_178[k] = -pk_106[k]
                   + f_0 * fk_286[k];

        t_179[k] = -pk_107[k]
                   + f_0 * fk_287[k];

        t_180[k] = f_0 * fk_288[k];

        t_181[k] = f_0 * fk_289[k];

        t_182[k] = f_0 * fk_290[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, t_188, t_189, t_190, fk_291, \
                         fk_292, fk_293, fk_294, fk_295, fk_296, fk_297, \
                         fk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_0 * fk_291[k];

        t_184[k] = f_0 * fk_292[k];

        t_185[k] = f_0 * fk_293[k];

        t_186[k] = f_0 * fk_294[k];

        t_187[k] = f_0 * fk_295[k];

        t_188[k] = f_0 * fk_296[k];

        t_189[k] = f_0 * fk_297[k];

        t_190[k] = f_0 * fk_298[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, t_197, t_198, fk_299, \
                         fk_300, fk_301, fk_302, fk_303, fk_304, fk_305, \
                         fk_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_0 * fk_299[k];

        t_192[k] = f_0 * fk_300[k];

        t_193[k] = f_0 * fk_301[k];

        t_194[k] = f_0 * fk_302[k];

        t_195[k] = f_0 * fk_303[k];

        t_196[k] = f_0 * fk_304[k];

        t_197[k] = f_0 * fk_305[k];

        t_198[k] = f_0 * fk_306[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, t_204, t_205, t_206, fk_307, \
                         fk_308, fk_309, fk_310, fk_311, fk_312, fk_313, \
                         fk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_0 * fk_307[k];

        t_200[k] = f_0 * fk_308[k];

        t_201[k] = f_0 * fk_309[k];

        t_202[k] = f_0 * fk_310[k];

        t_203[k] = f_0 * fk_311[k];

        t_204[k] = f_0 * fk_312[k];

        t_205[k] = f_0 * fk_313[k];

        t_206[k] = f_0 * fk_314[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, t_213, t_214, fk_315, \
                         fk_316, fk_317, fk_318, fk_319, fk_320, fk_321, \
                         fk_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_0 * fk_315[k];

        t_208[k] = f_0 * fk_316[k];

        t_209[k] = f_0 * fk_317[k];

        t_210[k] = f_0 * fk_318[k];

        t_211[k] = f_0 * fk_319[k];

        t_212[k] = f_0 * fk_320[k];

        t_213[k] = f_0 * fk_321[k];

        t_214[k] = f_0 * fk_322[k];
    }

#pragma omp simd aligned(t_215, fk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_0 * fk_323[k];
    }
}

auto
compute_prim_geom_10_dk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t pk, const size_t fk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_dk_electron_repulsion_1_piece0(buffer, target, pk, fk, ncols, alpha);

    compute_prim_geom_10_dk_electron_repulsion_1_piece1(buffer, target, pk, fk, ncols, alpha);
}

static auto
compute_prim_geom_10_dk_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t pk, const size_t fk,
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

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_11 = buffer.data(pk + 11);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_16 = buffer.data(pk + 16);
    const auto *pk_17 = buffer.data(pk + 17);
    const auto *pk_18 = buffer.data(pk + 18);
    const auto *pk_19 = buffer.data(pk + 19);
    const auto *pk_20 = buffer.data(pk + 20);
    const auto *pk_21 = buffer.data(pk + 21);
    const auto *pk_22 = buffer.data(pk + 22);
    const auto *pk_23 = buffer.data(pk + 23);
    const auto *pk_24 = buffer.data(pk + 24);
    const auto *pk_25 = buffer.data(pk + 25);
    const auto *pk_26 = buffer.data(pk + 26);
    const auto *pk_27 = buffer.data(pk + 27);
    const auto *pk_28 = buffer.data(pk + 28);
    const auto *pk_29 = buffer.data(pk + 29);
    const auto *pk_30 = buffer.data(pk + 30);
    const auto *pk_31 = buffer.data(pk + 31);
    const auto *pk_32 = buffer.data(pk + 32);
    const auto *pk_33 = buffer.data(pk + 33);
    const auto *pk_34 = buffer.data(pk + 34);
    const auto *pk_35 = buffer.data(pk + 35);
    const auto *pk_36 = buffer.data(pk + 36);
    const auto *pk_37 = buffer.data(pk + 37);
    const auto *pk_38 = buffer.data(pk + 38);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_40 = buffer.data(pk + 40);
    const auto *pk_41 = buffer.data(pk + 41);
    const auto *pk_42 = buffer.data(pk + 42);
    const auto *pk_43 = buffer.data(pk + 43);
    const auto *pk_44 = buffer.data(pk + 44);
    const auto *pk_45 = buffer.data(pk + 45);
    const auto *pk_46 = buffer.data(pk + 46);
    const auto *pk_47 = buffer.data(pk + 47);
    const auto *pk_48 = buffer.data(pk + 48);
    const auto *pk_49 = buffer.data(pk + 49);
    const auto *pk_50 = buffer.data(pk + 50);
    const auto *pk_51 = buffer.data(pk + 51);
    const auto *pk_52 = buffer.data(pk + 52);
    const auto *pk_53 = buffer.data(pk + 53);
    const auto *pk_54 = buffer.data(pk + 54);
    const auto *pk_55 = buffer.data(pk + 55);
    const auto *pk_56 = buffer.data(pk + 56);
    const auto *pk_57 = buffer.data(pk + 57);
    const auto *pk_58 = buffer.data(pk + 58);
    const auto *pk_59 = buffer.data(pk + 59);
    const auto *pk_60 = buffer.data(pk + 60);
    const auto *pk_61 = buffer.data(pk + 61);
    const auto *pk_62 = buffer.data(pk + 62);
    const auto *pk_63 = buffer.data(pk + 63);
    const auto *pk_64 = buffer.data(pk + 64);
    const auto *pk_65 = buffer.data(pk + 65);
    const auto *pk_66 = buffer.data(pk + 66);
    const auto *pk_67 = buffer.data(pk + 67);
    const auto *pk_68 = buffer.data(pk + 68);
    const auto *pk_69 = buffer.data(pk + 69);
    const auto *pk_70 = buffer.data(pk + 70);
    const auto *pk_71 = buffer.data(pk + 71);
    const auto *pk_72 = buffer.data(pk + 72);
    const auto *pk_73 = buffer.data(pk + 73);
    const auto *pk_74 = buffer.data(pk + 74);
    const auto *pk_75 = buffer.data(pk + 75);
    const auto *pk_76 = buffer.data(pk + 76);

    const auto *fk_72 = buffer.data(fk + 72);
    const auto *fk_73 = buffer.data(fk + 73);
    const auto *fk_74 = buffer.data(fk + 74);
    const auto *fk_75 = buffer.data(fk + 75);
    const auto *fk_76 = buffer.data(fk + 76);
    const auto *fk_77 = buffer.data(fk + 77);
    const auto *fk_78 = buffer.data(fk + 78);
    const auto *fk_79 = buffer.data(fk + 79);
    const auto *fk_80 = buffer.data(fk + 80);
    const auto *fk_81 = buffer.data(fk + 81);
    const auto *fk_82 = buffer.data(fk + 82);
    const auto *fk_83 = buffer.data(fk + 83);
    const auto *fk_84 = buffer.data(fk + 84);
    const auto *fk_85 = buffer.data(fk + 85);
    const auto *fk_86 = buffer.data(fk + 86);
    const auto *fk_87 = buffer.data(fk + 87);
    const auto *fk_88 = buffer.data(fk + 88);
    const auto *fk_89 = buffer.data(fk + 89);
    const auto *fk_90 = buffer.data(fk + 90);
    const auto *fk_91 = buffer.data(fk + 91);
    const auto *fk_92 = buffer.data(fk + 92);
    const auto *fk_93 = buffer.data(fk + 93);
    const auto *fk_94 = buffer.data(fk + 94);
    const auto *fk_95 = buffer.data(fk + 95);
    const auto *fk_96 = buffer.data(fk + 96);
    const auto *fk_97 = buffer.data(fk + 97);
    const auto *fk_98 = buffer.data(fk + 98);
    const auto *fk_99 = buffer.data(fk + 99);
    const auto *fk_100 = buffer.data(fk + 100);
    const auto *fk_101 = buffer.data(fk + 101);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_106 = buffer.data(fk + 106);
    const auto *fk_107 = buffer.data(fk + 107);
    const auto *fk_144 = buffer.data(fk + 144);
    const auto *fk_145 = buffer.data(fk + 145);
    const auto *fk_146 = buffer.data(fk + 146);
    const auto *fk_147 = buffer.data(fk + 147);
    const auto *fk_148 = buffer.data(fk + 148);
    const auto *fk_149 = buffer.data(fk + 149);
    const auto *fk_150 = buffer.data(fk + 150);
    const auto *fk_151 = buffer.data(fk + 151);
    const auto *fk_152 = buffer.data(fk + 152);
    const auto *fk_153 = buffer.data(fk + 153);
    const auto *fk_154 = buffer.data(fk + 154);
    const auto *fk_155 = buffer.data(fk + 155);
    const auto *fk_156 = buffer.data(fk + 156);
    const auto *fk_157 = buffer.data(fk + 157);
    const auto *fk_158 = buffer.data(fk + 158);
    const auto *fk_159 = buffer.data(fk + 159);
    const auto *fk_160 = buffer.data(fk + 160);
    const auto *fk_161 = buffer.data(fk + 161);
    const auto *fk_162 = buffer.data(fk + 162);
    const auto *fk_163 = buffer.data(fk + 163);
    const auto *fk_164 = buffer.data(fk + 164);
    const auto *fk_165 = buffer.data(fk + 165);
    const auto *fk_166 = buffer.data(fk + 166);
    const auto *fk_167 = buffer.data(fk + 167);
    const auto *fk_168 = buffer.data(fk + 168);
    const auto *fk_169 = buffer.data(fk + 169);
    const auto *fk_170 = buffer.data(fk + 170);
    const auto *fk_171 = buffer.data(fk + 171);
    const auto *fk_172 = buffer.data(fk + 172);
    const auto *fk_173 = buffer.data(fk + 173);
    const auto *fk_174 = buffer.data(fk + 174);
    const auto *fk_175 = buffer.data(fk + 175);
    const auto *fk_176 = buffer.data(fk + 176);
    const auto *fk_177 = buffer.data(fk + 177);
    const auto *fk_178 = buffer.data(fk + 178);
    const auto *fk_179 = buffer.data(fk + 179);
    const auto *fk_180 = buffer.data(fk + 180);
    const auto *fk_181 = buffer.data(fk + 181);
    const auto *fk_182 = buffer.data(fk + 182);
    const auto *fk_183 = buffer.data(fk + 183);
    const auto *fk_184 = buffer.data(fk + 184);
    const auto *fk_185 = buffer.data(fk + 185);
    const auto *fk_186 = buffer.data(fk + 186);
    const auto *fk_187 = buffer.data(fk + 187);
    const auto *fk_188 = buffer.data(fk + 188);
    const auto *fk_189 = buffer.data(fk + 189);
    const auto *fk_190 = buffer.data(fk + 190);
    const auto *fk_191 = buffer.data(fk + 191);
    const auto *fk_192 = buffer.data(fk + 192);
    const auto *fk_193 = buffer.data(fk + 193);
    const auto *fk_194 = buffer.data(fk + 194);
    const auto *fk_195 = buffer.data(fk + 195);
    const auto *fk_196 = buffer.data(fk + 196);
    const auto *fk_197 = buffer.data(fk + 197);
    const auto *fk_198 = buffer.data(fk + 198);
    const auto *fk_199 = buffer.data(fk + 199);
    const auto *fk_200 = buffer.data(fk + 200);
    const auto *fk_201 = buffer.data(fk + 201);
    const auto *fk_202 = buffer.data(fk + 202);
    const auto *fk_203 = buffer.data(fk + 203);
    const auto *fk_204 = buffer.data(fk + 204);
    const auto *fk_205 = buffer.data(fk + 205);
    const auto *fk_206 = buffer.data(fk + 206);
    const auto *fk_207 = buffer.data(fk + 207);
    const auto *fk_208 = buffer.data(fk + 208);
    const auto *fk_209 = buffer.data(fk + 209);
    const auto *fk_210 = buffer.data(fk + 210);
    const auto *fk_211 = buffer.data(fk + 211);
    const auto *fk_212 = buffer.data(fk + 212);
    const auto *fk_213 = buffer.data(fk + 213);
    const auto *fk_214 = buffer.data(fk + 214);
    const auto *fk_215 = buffer.data(fk + 215);
    const auto *fk_252 = buffer.data(fk + 252);
    const auto *fk_253 = buffer.data(fk + 253);
    const auto *fk_254 = buffer.data(fk + 254);
    const auto *fk_255 = buffer.data(fk + 255);
    const auto *fk_256 = buffer.data(fk + 256);
    const auto *fk_257 = buffer.data(fk + 257);
    const auto *fk_258 = buffer.data(fk + 258);
    const auto *fk_259 = buffer.data(fk + 259);
    const auto *fk_260 = buffer.data(fk + 260);
    const auto *fk_261 = buffer.data(fk + 261);
    const auto *fk_262 = buffer.data(fk + 262);
    const auto *fk_263 = buffer.data(fk + 263);
    const auto *fk_264 = buffer.data(fk + 264);
    const auto *fk_265 = buffer.data(fk + 265);
    const auto *fk_266 = buffer.data(fk + 266);
    const auto *fk_267 = buffer.data(fk + 267);
    const auto *fk_268 = buffer.data(fk + 268);
    const auto *fk_269 = buffer.data(fk + 269);
    const auto *fk_270 = buffer.data(fk + 270);
    const auto *fk_271 = buffer.data(fk + 271);
    const auto *fk_272 = buffer.data(fk + 272);
    const auto *fk_273 = buffer.data(fk + 273);
    const auto *fk_274 = buffer.data(fk + 274);
    const auto *fk_275 = buffer.data(fk + 275);
    const auto *fk_276 = buffer.data(fk + 276);
    const auto *fk_277 = buffer.data(fk + 277);
    const auto *fk_278 = buffer.data(fk + 278);
    const auto *fk_279 = buffer.data(fk + 279);
    const auto *fk_280 = buffer.data(fk + 280);
    const auto *fk_281 = buffer.data(fk + 281);
    const auto *fk_282 = buffer.data(fk + 282);
    const auto *fk_283 = buffer.data(fk + 283);
    const auto *fk_284 = buffer.data(fk + 284);
    const auto *fk_285 = buffer.data(fk + 285);
    const auto *fk_286 = buffer.data(fk + 286);
    const auto *fk_287 = buffer.data(fk + 287);
    const auto *fk_288 = buffer.data(fk + 288);
    const auto *fk_289 = buffer.data(fk + 289);
    const auto *fk_290 = buffer.data(fk + 290);
    const auto *fk_291 = buffer.data(fk + 291);
    const auto *fk_292 = buffer.data(fk + 292);
    const auto *fk_293 = buffer.data(fk + 293);
    const auto *fk_294 = buffer.data(fk + 294);
    const auto *fk_295 = buffer.data(fk + 295);
    const auto *fk_296 = buffer.data(fk + 296);
    const auto *fk_297 = buffer.data(fk + 297);
    const auto *fk_298 = buffer.data(fk + 298);
    const auto *fk_299 = buffer.data(fk + 299);
    const auto *fk_300 = buffer.data(fk + 300);
    const auto *fk_301 = buffer.data(fk + 301);
    const auto *fk_302 = buffer.data(fk + 302);
    const auto *fk_303 = buffer.data(fk + 303);
    const auto *fk_304 = buffer.data(fk + 304);
    const auto *fk_305 = buffer.data(fk + 305);
    const auto *fk_306 = buffer.data(fk + 306);
    const auto *fk_307 = buffer.data(fk + 307);
    const auto *fk_308 = buffer.data(fk + 308);
    const auto *fk_309 = buffer.data(fk + 309);
    const auto *fk_310 = buffer.data(fk + 310);
    const auto *fk_311 = buffer.data(fk + 311);
    const auto *fk_312 = buffer.data(fk + 312);
    const auto *fk_313 = buffer.data(fk + 313);
    const auto *fk_314 = buffer.data(fk + 314);
    const auto *fk_315 = buffer.data(fk + 315);
    const auto *fk_316 = buffer.data(fk + 316);
    const auto *fk_317 = buffer.data(fk + 317);
    const auto *fk_318 = buffer.data(fk + 318);
    const auto *fk_319 = buffer.data(fk + 319);
    const auto *fk_320 = buffer.data(fk + 320);
    const auto *fk_321 = buffer.data(fk + 321);
    const auto *fk_322 = buffer.data(fk + 322);
    const auto *fk_323 = buffer.data(fk + 323);
    const auto *fk_324 = buffer.data(fk + 324);
    const auto *fk_325 = buffer.data(fk + 325);
    const auto *fk_326 = buffer.data(fk + 326);
    const auto *fk_327 = buffer.data(fk + 327);
    const auto *fk_328 = buffer.data(fk + 328);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, fk_72, fk_73, fk_74, fk_75, \
                         fk_76, fk_77, fk_78, fk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fk_72[k];

        t_1[k] = f_0 * fk_73[k];

        t_2[k] = f_0 * fk_74[k];

        t_3[k] = f_0 * fk_75[k];

        t_4[k] = f_0 * fk_76[k];

        t_5[k] = f_0 * fk_77[k];

        t_6[k] = f_0 * fk_78[k];

        t_7[k] = f_0 * fk_79[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, fk_80, fk_81, fk_82, \
                         fk_83, fk_84, fk_85, fk_86, fk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * fk_80[k];

        t_9[k] = f_0 * fk_81[k];

        t_10[k] = f_0 * fk_82[k];

        t_11[k] = f_0 * fk_83[k];

        t_12[k] = f_0 * fk_84[k];

        t_13[k] = f_0 * fk_85[k];

        t_14[k] = f_0 * fk_86[k];

        t_15[k] = f_0 * fk_87[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, fk_88, fk_89, fk_90, \
                         fk_91, fk_92, fk_93, fk_94, fk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * fk_88[k];

        t_17[k] = f_0 * fk_89[k];

        t_18[k] = f_0 * fk_90[k];

        t_19[k] = f_0 * fk_91[k];

        t_20[k] = f_0 * fk_92[k];

        t_21[k] = f_0 * fk_93[k];

        t_22[k] = f_0 * fk_94[k];

        t_23[k] = f_0 * fk_95[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, fk_96, fk_97, fk_98, \
                         fk_99, fk_100, fk_101, fk_102, fk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * fk_96[k];

        t_25[k] = f_0 * fk_97[k];

        t_26[k] = f_0 * fk_98[k];

        t_27[k] = f_0 * fk_99[k];

        t_28[k] = f_0 * fk_100[k];

        t_29[k] = f_0 * fk_101[k];

        t_30[k] = f_0 * fk_102[k];

        t_31[k] = f_0 * fk_103[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, fk_104, fk_105, \
                         fk_106, fk_107, fk_144, fk_145, fk_146, \
                         fk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * fk_104[k];

        t_33[k] = f_0 * fk_105[k];

        t_34[k] = f_0 * fk_106[k];

        t_35[k] = f_0 * fk_107[k];

        t_36[k] = f_0 * fk_144[k];

        t_37[k] = f_0 * fk_145[k];

        t_38[k] = f_0 * fk_146[k];

        t_39[k] = f_0 * fk_147[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, fk_148, fk_149, \
                         fk_150, fk_151, fk_152, fk_153, fk_154, \
                         fk_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * fk_148[k];

        t_41[k] = f_0 * fk_149[k];

        t_42[k] = f_0 * fk_150[k];

        t_43[k] = f_0 * fk_151[k];

        t_44[k] = f_0 * fk_152[k];

        t_45[k] = f_0 * fk_153[k];

        t_46[k] = f_0 * fk_154[k];

        t_47[k] = f_0 * fk_155[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, fk_156, fk_157, \
                         fk_158, fk_159, fk_160, fk_161, fk_162, \
                         fk_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * fk_156[k];

        t_49[k] = f_0 * fk_157[k];

        t_50[k] = f_0 * fk_158[k];

        t_51[k] = f_0 * fk_159[k];

        t_52[k] = f_0 * fk_160[k];

        t_53[k] = f_0 * fk_161[k];

        t_54[k] = f_0 * fk_162[k];

        t_55[k] = f_0 * fk_163[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, fk_164, fk_165, \
                         fk_166, fk_167, fk_168, fk_169, fk_170, \
                         fk_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * fk_164[k];

        t_57[k] = f_0 * fk_165[k];

        t_58[k] = f_0 * fk_166[k];

        t_59[k] = f_0 * fk_167[k];

        t_60[k] = f_0 * fk_168[k];

        t_61[k] = f_0 * fk_169[k];

        t_62[k] = f_0 * fk_170[k];

        t_63[k] = f_0 * fk_171[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, fk_172, fk_173, \
                         fk_174, fk_175, fk_176, fk_177, fk_178, \
                         fk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * fk_172[k];

        t_65[k] = f_0 * fk_173[k];

        t_66[k] = f_0 * fk_174[k];

        t_67[k] = f_0 * fk_175[k];

        t_68[k] = f_0 * fk_176[k];

        t_69[k] = f_0 * fk_177[k];

        t_70[k] = f_0 * fk_178[k];

        t_71[k] = f_0 * fk_179[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pk_0, pk_1, pk_2, pk_3, pk_4, fk_180, \
                         fk_181, fk_182, fk_183, fk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -pk_0[k]
                  + f_0 * fk_180[k];

        t_73[k] = -pk_1[k]
                  + f_0 * fk_181[k];

        t_74[k] = -pk_2[k]
                  + f_0 * fk_182[k];

        t_75[k] = -pk_3[k]
                  + f_0 * fk_183[k];

        t_76[k] = -pk_4[k]
                  + f_0 * fk_184[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, pk_5, pk_6, pk_7, pk_8, pk_9, fk_185, \
                         fk_186, fk_187, fk_188, fk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -pk_5[k]
                  + f_0 * fk_185[k];

        t_78[k] = -pk_6[k]
                  + f_0 * fk_186[k];

        t_79[k] = -pk_7[k]
                  + f_0 * fk_187[k];

        t_80[k] = -pk_8[k]
                  + f_0 * fk_188[k];

        t_81[k] = -pk_9[k]
                  + f_0 * fk_189[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pk_10, pk_11, pk_12, pk_13, pk_14, \
                         fk_190, fk_191, fk_192, fk_193, fk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -pk_10[k]
                  + f_0 * fk_190[k];

        t_83[k] = -pk_11[k]
                  + f_0 * fk_191[k];

        t_84[k] = -pk_12[k]
                  + f_0 * fk_192[k];

        t_85[k] = -pk_13[k]
                  + f_0 * fk_193[k];

        t_86[k] = -pk_14[k]
                  + f_0 * fk_194[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pk_15, pk_16, pk_17, pk_18, pk_19, \
                         fk_195, fk_196, fk_197, fk_198, fk_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -pk_15[k]
                  + f_0 * fk_195[k];

        t_88[k] = -pk_16[k]
                  + f_0 * fk_196[k];

        t_89[k] = -pk_17[k]
                  + f_0 * fk_197[k];

        t_90[k] = -pk_18[k]
                  + f_0 * fk_198[k];

        t_91[k] = -pk_19[k]
                  + f_0 * fk_199[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, pk_20, pk_21, pk_22, pk_23, pk_24, \
                         fk_200, fk_201, fk_202, fk_203, fk_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = -pk_20[k]
                  + f_0 * fk_200[k];

        t_93[k] = -pk_21[k]
                  + f_0 * fk_201[k];

        t_94[k] = -pk_22[k]
                  + f_0 * fk_202[k];

        t_95[k] = -pk_23[k]
                  + f_0 * fk_203[k];

        t_96[k] = -pk_24[k]
                  + f_0 * fk_204[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, pk_25, pk_26, pk_27, pk_28, pk_29, \
                         fk_205, fk_206, fk_207, fk_208, fk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -pk_25[k]
                  + f_0 * fk_205[k];

        t_98[k] = -pk_26[k]
                  + f_0 * fk_206[k];

        t_99[k] = -pk_27[k]
                  + f_0 * fk_207[k];

        t_100[k] = -pk_28[k]
                   + f_0 * fk_208[k];

        t_101[k] = -pk_29[k]
                   + f_0 * fk_209[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, pk_30, pk_31, pk_32, pk_33, pk_34, \
                         fk_210, fk_211, fk_212, fk_213, fk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -pk_30[k]
                   + f_0 * fk_210[k];

        t_103[k] = -pk_31[k]
                   + f_0 * fk_211[k];

        t_104[k] = -pk_32[k]
                   + f_0 * fk_212[k];

        t_105[k] = -pk_33[k]
                   + f_0 * fk_213[k];

        t_106[k] = -pk_34[k]
                   + f_0 * fk_214[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, pk_35, fk_215, \
                         fk_252, fk_253, fk_254, fk_255, fk_256, \
                         fk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -pk_35[k]
                   + f_0 * fk_215[k];

        t_108[k] = f_0 * fk_252[k];

        t_109[k] = f_0 * fk_253[k];

        t_110[k] = f_0 * fk_254[k];

        t_111[k] = f_0 * fk_255[k];

        t_112[k] = f_0 * fk_256[k];

        t_113[k] = f_0 * fk_257[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, t_120, t_121, fk_258, \
                         fk_259, fk_260, fk_261, fk_262, fk_263, fk_264, \
                         fk_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * fk_258[k];

        t_115[k] = f_0 * fk_259[k];

        t_116[k] = f_0 * fk_260[k];

        t_117[k] = f_0 * fk_261[k];

        t_118[k] = f_0 * fk_262[k];

        t_119[k] = f_0 * fk_263[k];

        t_120[k] = f_0 * fk_264[k];

        t_121[k] = f_0 * fk_265[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, t_128, t_129, fk_266, \
                         fk_267, fk_268, fk_269, fk_270, fk_271, fk_272, \
                         fk_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * fk_266[k];

        t_123[k] = f_0 * fk_267[k];

        t_124[k] = f_0 * fk_268[k];

        t_125[k] = f_0 * fk_269[k];

        t_126[k] = f_0 * fk_270[k];

        t_127[k] = f_0 * fk_271[k];

        t_128[k] = f_0 * fk_272[k];

        t_129[k] = f_0 * fk_273[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, t_135, t_136, t_137, fk_274, \
                         fk_275, fk_276, fk_277, fk_278, fk_279, fk_280, \
                         fk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * fk_274[k];

        t_131[k] = f_0 * fk_275[k];

        t_132[k] = f_0 * fk_276[k];

        t_133[k] = f_0 * fk_277[k];

        t_134[k] = f_0 * fk_278[k];

        t_135[k] = f_0 * fk_279[k];

        t_136[k] = f_0 * fk_280[k];

        t_137[k] = f_0 * fk_281[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, t_144, pk_36, fk_282, \
                         fk_283, fk_284, fk_285, fk_286, fk_287, \
                         fk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * fk_282[k];

        t_139[k] = f_0 * fk_283[k];

        t_140[k] = f_0 * fk_284[k];

        t_141[k] = f_0 * fk_285[k];

        t_142[k] = f_0 * fk_286[k];

        t_143[k] = f_0 * fk_287[k];

        t_144[k] = -pk_36[k]
                   + f_0 * fk_288[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pk_37, pk_38, pk_39, pk_40, pk_41, \
                         fk_289, fk_290, fk_291, fk_292, fk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -pk_37[k]
                   + f_0 * fk_289[k];

        t_146[k] = -pk_38[k]
                   + f_0 * fk_290[k];

        t_147[k] = -pk_39[k]
                   + f_0 * fk_291[k];

        t_148[k] = -pk_40[k]
                   + f_0 * fk_292[k];

        t_149[k] = -pk_41[k]
                   + f_0 * fk_293[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pk_42, pk_43, pk_44, pk_45, pk_46, \
                         fk_294, fk_295, fk_296, fk_297, fk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -pk_42[k]
                   + f_0 * fk_294[k];

        t_151[k] = -pk_43[k]
                   + f_0 * fk_295[k];

        t_152[k] = -pk_44[k]
                   + f_0 * fk_296[k];

        t_153[k] = -pk_45[k]
                   + f_0 * fk_297[k];

        t_154[k] = -pk_46[k]
                   + f_0 * fk_298[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pk_47, pk_48, pk_49, pk_50, pk_51, \
                         fk_299, fk_300, fk_301, fk_302, fk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -pk_47[k]
                   + f_0 * fk_299[k];

        t_156[k] = -pk_48[k]
                   + f_0 * fk_300[k];

        t_157[k] = -pk_49[k]
                   + f_0 * fk_301[k];

        t_158[k] = -pk_50[k]
                   + f_0 * fk_302[k];

        t_159[k] = -pk_51[k]
                   + f_0 * fk_303[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pk_52, pk_53, pk_54, pk_55, pk_56, \
                         fk_304, fk_305, fk_306, fk_307, fk_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -pk_52[k]
                   + f_0 * fk_304[k];

        t_161[k] = -pk_53[k]
                   + f_0 * fk_305[k];

        t_162[k] = -pk_54[k]
                   + f_0 * fk_306[k];

        t_163[k] = -pk_55[k]
                   + f_0 * fk_307[k];

        t_164[k] = -pk_56[k]
                   + f_0 * fk_308[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pk_57, pk_58, pk_59, pk_60, pk_61, \
                         fk_309, fk_310, fk_311, fk_312, fk_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -pk_57[k]
                   + f_0 * fk_309[k];

        t_166[k] = -pk_58[k]
                   + f_0 * fk_310[k];

        t_167[k] = -pk_59[k]
                   + f_0 * fk_311[k];

        t_168[k] = -pk_60[k]
                   + f_0 * fk_312[k];

        t_169[k] = -pk_61[k]
                   + f_0 * fk_313[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, pk_62, pk_63, pk_64, pk_65, pk_66, \
                         fk_314, fk_315, fk_316, fk_317, fk_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -pk_62[k]
                   + f_0 * fk_314[k];

        t_171[k] = -pk_63[k]
                   + f_0 * fk_315[k];

        t_172[k] = -pk_64[k]
                   + f_0 * fk_316[k];

        t_173[k] = -pk_65[k]
                   + f_0 * fk_317[k];

        t_174[k] = -pk_66[k]
                   + f_0 * fk_318[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, pk_67, pk_68, pk_69, pk_70, pk_71, \
                         fk_319, fk_320, fk_321, fk_322, fk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -pk_67[k]
                   + f_0 * fk_319[k];

        t_176[k] = -pk_68[k]
                   + f_0 * fk_320[k];

        t_177[k] = -pk_69[k]
                   + f_0 * fk_321[k];

        t_178[k] = -pk_70[k]
                   + f_0 * fk_322[k];

        t_179[k] = -pk_71[k]
                   + f_0 * fk_323[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pk_72, pk_73, pk_74, pk_75, pk_76, \
                         fk_324, fk_325, fk_326, fk_327, fk_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * pk_72[k]
                   + f_0 * fk_324[k];

        t_181[k] = -2.0 * pk_73[k]
                   + f_0 * fk_325[k];

        t_182[k] = -2.0 * pk_74[k]
                   + f_0 * fk_326[k];

        t_183[k] = -2.0 * pk_75[k]
                   + f_0 * fk_327[k];

        t_184[k] = -2.0 * pk_76[k]
                   + f_0 * fk_328[k];
    }
}

static auto
compute_prim_geom_10_dk_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t pk, const size_t fk,
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

    const auto *pk_77 = buffer.data(pk + 77);
    const auto *pk_78 = buffer.data(pk + 78);
    const auto *pk_79 = buffer.data(pk + 79);
    const auto *pk_80 = buffer.data(pk + 80);
    const auto *pk_81 = buffer.data(pk + 81);
    const auto *pk_82 = buffer.data(pk + 82);
    const auto *pk_83 = buffer.data(pk + 83);
    const auto *pk_84 = buffer.data(pk + 84);
    const auto *pk_85 = buffer.data(pk + 85);
    const auto *pk_86 = buffer.data(pk + 86);
    const auto *pk_87 = buffer.data(pk + 87);
    const auto *pk_88 = buffer.data(pk + 88);
    const auto *pk_89 = buffer.data(pk + 89);
    const auto *pk_90 = buffer.data(pk + 90);
    const auto *pk_91 = buffer.data(pk + 91);
    const auto *pk_92 = buffer.data(pk + 92);
    const auto *pk_93 = buffer.data(pk + 93);
    const auto *pk_94 = buffer.data(pk + 94);
    const auto *pk_95 = buffer.data(pk + 95);
    const auto *pk_96 = buffer.data(pk + 96);
    const auto *pk_97 = buffer.data(pk + 97);
    const auto *pk_98 = buffer.data(pk + 98);
    const auto *pk_99 = buffer.data(pk + 99);
    const auto *pk_100 = buffer.data(pk + 100);
    const auto *pk_101 = buffer.data(pk + 101);
    const auto *pk_102 = buffer.data(pk + 102);
    const auto *pk_103 = buffer.data(pk + 103);
    const auto *pk_104 = buffer.data(pk + 104);
    const auto *pk_105 = buffer.data(pk + 105);
    const auto *pk_106 = buffer.data(pk + 106);
    const auto *pk_107 = buffer.data(pk + 107);

    const auto *fk_329 = buffer.data(fk + 329);
    const auto *fk_330 = buffer.data(fk + 330);
    const auto *fk_331 = buffer.data(fk + 331);
    const auto *fk_332 = buffer.data(fk + 332);
    const auto *fk_333 = buffer.data(fk + 333);
    const auto *fk_334 = buffer.data(fk + 334);
    const auto *fk_335 = buffer.data(fk + 335);
    const auto *fk_336 = buffer.data(fk + 336);
    const auto *fk_337 = buffer.data(fk + 337);
    const auto *fk_338 = buffer.data(fk + 338);
    const auto *fk_339 = buffer.data(fk + 339);
    const auto *fk_340 = buffer.data(fk + 340);
    const auto *fk_341 = buffer.data(fk + 341);
    const auto *fk_342 = buffer.data(fk + 342);
    const auto *fk_343 = buffer.data(fk + 343);
    const auto *fk_344 = buffer.data(fk + 344);
    const auto *fk_345 = buffer.data(fk + 345);
    const auto *fk_346 = buffer.data(fk + 346);
    const auto *fk_347 = buffer.data(fk + 347);
    const auto *fk_348 = buffer.data(fk + 348);
    const auto *fk_349 = buffer.data(fk + 349);
    const auto *fk_350 = buffer.data(fk + 350);
    const auto *fk_351 = buffer.data(fk + 351);
    const auto *fk_352 = buffer.data(fk + 352);
    const auto *fk_353 = buffer.data(fk + 353);
    const auto *fk_354 = buffer.data(fk + 354);
    const auto *fk_355 = buffer.data(fk + 355);
    const auto *fk_356 = buffer.data(fk + 356);
    const auto *fk_357 = buffer.data(fk + 357);
    const auto *fk_358 = buffer.data(fk + 358);
    const auto *fk_359 = buffer.data(fk + 359);

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, pk_77, pk_78, pk_79, pk_80, pk_81, \
                         fk_329, fk_330, fk_331, fk_332, fk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * pk_77[k]
                   + f_0 * fk_329[k];

        t_186[k] = -2.0 * pk_78[k]
                   + f_0 * fk_330[k];

        t_187[k] = -2.0 * pk_79[k]
                   + f_0 * fk_331[k];

        t_188[k] = -2.0 * pk_80[k]
                   + f_0 * fk_332[k];

        t_189[k] = -2.0 * pk_81[k]
                   + f_0 * fk_333[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, pk_82, pk_83, pk_84, pk_85, pk_86, \
                         fk_334, fk_335, fk_336, fk_337, fk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * pk_82[k]
                   + f_0 * fk_334[k];

        t_191[k] = -2.0 * pk_83[k]
                   + f_0 * fk_335[k];

        t_192[k] = -2.0 * pk_84[k]
                   + f_0 * fk_336[k];

        t_193[k] = -2.0 * pk_85[k]
                   + f_0 * fk_337[k];

        t_194[k] = -2.0 * pk_86[k]
                   + f_0 * fk_338[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, pk_87, pk_88, pk_89, pk_90, pk_91, \
                         fk_339, fk_340, fk_341, fk_342, fk_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * pk_87[k]
                   + f_0 * fk_339[k];

        t_196[k] = -2.0 * pk_88[k]
                   + f_0 * fk_340[k];

        t_197[k] = -2.0 * pk_89[k]
                   + f_0 * fk_341[k];

        t_198[k] = -2.0 * pk_90[k]
                   + f_0 * fk_342[k];

        t_199[k] = -2.0 * pk_91[k]
                   + f_0 * fk_343[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pk_92, pk_93, pk_94, pk_95, pk_96, \
                         fk_344, fk_345, fk_346, fk_347, fk_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * pk_92[k]
                   + f_0 * fk_344[k];

        t_201[k] = -2.0 * pk_93[k]
                   + f_0 * fk_345[k];

        t_202[k] = -2.0 * pk_94[k]
                   + f_0 * fk_346[k];

        t_203[k] = -2.0 * pk_95[k]
                   + f_0 * fk_347[k];

        t_204[k] = -2.0 * pk_96[k]
                   + f_0 * fk_348[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, pk_97, pk_98, pk_99, pk_100, \
                         pk_101, fk_349, fk_350, fk_351, fk_352, \
                         fk_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * pk_97[k]
                   + f_0 * fk_349[k];

        t_206[k] = -2.0 * pk_98[k]
                   + f_0 * fk_350[k];

        t_207[k] = -2.0 * pk_99[k]
                   + f_0 * fk_351[k];

        t_208[k] = -2.0 * pk_100[k]
                   + f_0 * fk_352[k];

        t_209[k] = -2.0 * pk_101[k]
                   + f_0 * fk_353[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pk_102, pk_103, pk_104, pk_105, \
                         pk_106, fk_354, fk_355, fk_356, fk_357, \
                         fk_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * pk_102[k]
                   + f_0 * fk_354[k];

        t_211[k] = -2.0 * pk_103[k]
                   + f_0 * fk_355[k];

        t_212[k] = -2.0 * pk_104[k]
                   + f_0 * fk_356[k];

        t_213[k] = -2.0 * pk_105[k]
                   + f_0 * fk_357[k];

        t_214[k] = -2.0 * pk_106[k]
                   + f_0 * fk_358[k];
    }

#pragma omp simd aligned(t_215, pk_107, fk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * pk_107[k]
                   + f_0 * fk_359[k];
    }
}

auto
compute_prim_geom_10_dk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t pk, const size_t fk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_dk_electron_repulsion_2_piece0(buffer, target, pk, fk, ncols, alpha);

    compute_prim_geom_10_dk_electron_repulsion_2_piece1(buffer, target, pk, fk, ncols, alpha);
}

}  // namespace simdt2ceri
