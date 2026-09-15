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


#include "SimdElectronRepulsionGeom10VrrRecID.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_id_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hd, const size_t kd,
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

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_73 = buffer.data(hd + 73);
    const auto *hd_74 = buffer.data(hd + 74);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_79 = buffer.data(hd + 79);
    const auto *hd_80 = buffer.data(hd + 80);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_83 = buffer.data(hd + 83);
    const auto *hd_84 = buffer.data(hd + 84);
    const auto *hd_85 = buffer.data(hd + 85);
    const auto *hd_86 = buffer.data(hd + 86);
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_88 = buffer.data(hd + 88);
    const auto *hd_89 = buffer.data(hd + 89);
    const auto *hd_90 = buffer.data(hd + 90);
    const auto *hd_91 = buffer.data(hd + 91);
    const auto *hd_92 = buffer.data(hd + 92);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_94 = buffer.data(hd + 94);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_96 = buffer.data(hd + 96);
    const auto *hd_97 = buffer.data(hd + 97);
    const auto *hd_98 = buffer.data(hd + 98);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_102 = buffer.data(hd + 102);
    const auto *hd_103 = buffer.data(hd + 103);
    const auto *hd_104 = buffer.data(hd + 104);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_108 = buffer.data(hd + 108);
    const auto *hd_109 = buffer.data(hd + 109);
    const auto *hd_110 = buffer.data(hd + 110);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_114 = buffer.data(hd + 114);
    const auto *hd_115 = buffer.data(hd + 115);
    const auto *hd_116 = buffer.data(hd + 116);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_119 = buffer.data(hd + 119);
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_121 = buffer.data(hd + 121);
    const auto *hd_122 = buffer.data(hd + 122);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_124 = buffer.data(hd + 124);
    const auto *hd_125 = buffer.data(hd + 125);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_91 = buffer.data(kd + 91);
    const auto *kd_92 = buffer.data(kd + 92);
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_94 = buffer.data(kd + 94);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_97 = buffer.data(kd + 97);
    const auto *kd_98 = buffer.data(kd + 98);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_103 = buffer.data(kd + 103);
    const auto *kd_104 = buffer.data(kd + 104);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_109 = buffer.data(kd + 109);
    const auto *kd_110 = buffer.data(kd + 110);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_115 = buffer.data(kd + 115);
    const auto *kd_116 = buffer.data(kd + 116);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_121 = buffer.data(kd + 121);
    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_124 = buffer.data(kd + 124);
    const auto *kd_125 = buffer.data(kd + 125);
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_128 = buffer.data(kd + 128);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_130 = buffer.data(kd + 130);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_133 = buffer.data(kd + 133);
    const auto *kd_134 = buffer.data(kd + 134);
    const auto *kd_135 = buffer.data(kd + 135);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_139 = buffer.data(kd + 139);
    const auto *kd_140 = buffer.data(kd + 140);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_145 = buffer.data(kd + 145);
    const auto *kd_146 = buffer.data(kd + 146);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_151 = buffer.data(kd + 151);
    const auto *kd_152 = buffer.data(kd + 152);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_157 = buffer.data(kd + 157);
    const auto *kd_158 = buffer.data(kd + 158);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_161 = buffer.data(kd + 161);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_163 = buffer.data(kd + 163);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, hd_0, hd_1, hd_2, hd_3, hd_4, kd_0, kd_1, \
                         kd_2, kd_3, kd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -6.0 * hd_0[k]
                 + f_0 * kd_0[k];

        t_1[k] = -6.0 * hd_1[k]
                 + f_0 * kd_1[k];

        t_2[k] = -6.0 * hd_2[k]
                 + f_0 * kd_2[k];

        t_3[k] = -6.0 * hd_3[k]
                 + f_0 * kd_3[k];

        t_4[k] = -6.0 * hd_4[k]
                 + f_0 * kd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, hd_5, hd_6, hd_7, hd_8, hd_9, kd_5, kd_6, \
                         kd_7, kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -6.0 * hd_5[k]
                 + f_0 * kd_5[k];

        t_6[k] = -5.0 * hd_6[k]
                 + f_0 * kd_6[k];

        t_7[k] = -5.0 * hd_7[k]
                 + f_0 * kd_7[k];

        t_8[k] = -5.0 * hd_8[k]
                 + f_0 * kd_8[k];

        t_9[k] = -5.0 * hd_9[k]
                 + f_0 * kd_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, hd_10, hd_11, hd_12, hd_13, hd_14, \
                         kd_10, kd_11, kd_12, kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -5.0 * hd_10[k]
                  + f_0 * kd_10[k];

        t_11[k] = -5.0 * hd_11[k]
                  + f_0 * kd_11[k];

        t_12[k] = -5.0 * hd_12[k]
                  + f_0 * kd_12[k];

        t_13[k] = -5.0 * hd_13[k]
                  + f_0 * kd_13[k];

        t_14[k] = -5.0 * hd_14[k]
                  + f_0 * kd_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, hd_15, hd_16, hd_17, hd_18, hd_19, \
                         kd_15, kd_16, kd_17, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -5.0 * hd_15[k]
                  + f_0 * kd_15[k];

        t_16[k] = -5.0 * hd_16[k]
                  + f_0 * kd_16[k];

        t_17[k] = -5.0 * hd_17[k]
                  + f_0 * kd_17[k];

        t_18[k] = -4.0 * hd_18[k]
                  + f_0 * kd_18[k];

        t_19[k] = -4.0 * hd_19[k]
                  + f_0 * kd_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, hd_20, hd_21, hd_22, hd_23, hd_24, \
                         kd_20, kd_21, kd_22, kd_23, kd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -4.0 * hd_20[k]
                  + f_0 * kd_20[k];

        t_21[k] = -4.0 * hd_21[k]
                  + f_0 * kd_21[k];

        t_22[k] = -4.0 * hd_22[k]
                  + f_0 * kd_22[k];

        t_23[k] = -4.0 * hd_23[k]
                  + f_0 * kd_23[k];

        t_24[k] = -4.0 * hd_24[k]
                  + f_0 * kd_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, hd_25, hd_26, hd_27, hd_28, hd_29, \
                         kd_25, kd_26, kd_27, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -4.0 * hd_25[k]
                  + f_0 * kd_25[k];

        t_26[k] = -4.0 * hd_26[k]
                  + f_0 * kd_26[k];

        t_27[k] = -4.0 * hd_27[k]
                  + f_0 * kd_27[k];

        t_28[k] = -4.0 * hd_28[k]
                  + f_0 * kd_28[k];

        t_29[k] = -4.0 * hd_29[k]
                  + f_0 * kd_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, hd_30, hd_31, hd_32, hd_33, hd_34, \
                         kd_30, kd_31, kd_32, kd_33, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -4.0 * hd_30[k]
                  + f_0 * kd_30[k];

        t_31[k] = -4.0 * hd_31[k]
                  + f_0 * kd_31[k];

        t_32[k] = -4.0 * hd_32[k]
                  + f_0 * kd_32[k];

        t_33[k] = -4.0 * hd_33[k]
                  + f_0 * kd_33[k];

        t_34[k] = -4.0 * hd_34[k]
                  + f_0 * kd_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, hd_35, hd_36, hd_37, hd_38, hd_39, \
                         kd_35, kd_36, kd_37, kd_38, kd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -4.0 * hd_35[k]
                  + f_0 * kd_35[k];

        t_36[k] = -3.0 * hd_36[k]
                  + f_0 * kd_36[k];

        t_37[k] = -3.0 * hd_37[k]
                  + f_0 * kd_37[k];

        t_38[k] = -3.0 * hd_38[k]
                  + f_0 * kd_38[k];

        t_39[k] = -3.0 * hd_39[k]
                  + f_0 * kd_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, hd_40, hd_41, hd_42, hd_43, hd_44, \
                         kd_40, kd_41, kd_42, kd_43, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * hd_40[k]
                  + f_0 * kd_40[k];

        t_41[k] = -3.0 * hd_41[k]
                  + f_0 * kd_41[k];

        t_42[k] = -3.0 * hd_42[k]
                  + f_0 * kd_42[k];

        t_43[k] = -3.0 * hd_43[k]
                  + f_0 * kd_43[k];

        t_44[k] = -3.0 * hd_44[k]
                  + f_0 * kd_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hd_45, hd_46, hd_47, hd_48, hd_49, \
                         kd_45, kd_46, kd_47, kd_48, kd_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -3.0 * hd_45[k]
                  + f_0 * kd_45[k];

        t_46[k] = -3.0 * hd_46[k]
                  + f_0 * kd_46[k];

        t_47[k] = -3.0 * hd_47[k]
                  + f_0 * kd_47[k];

        t_48[k] = -3.0 * hd_48[k]
                  + f_0 * kd_48[k];

        t_49[k] = -3.0 * hd_49[k]
                  + f_0 * kd_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hd_50, hd_51, hd_52, hd_53, hd_54, \
                         kd_50, kd_51, kd_52, kd_53, kd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -3.0 * hd_50[k]
                  + f_0 * kd_50[k];

        t_51[k] = -3.0 * hd_51[k]
                  + f_0 * kd_51[k];

        t_52[k] = -3.0 * hd_52[k]
                  + f_0 * kd_52[k];

        t_53[k] = -3.0 * hd_53[k]
                  + f_0 * kd_53[k];

        t_54[k] = -3.0 * hd_54[k]
                  + f_0 * kd_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, hd_55, hd_56, hd_57, hd_58, hd_59, \
                         kd_55, kd_56, kd_57, kd_58, kd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -3.0 * hd_55[k]
                  + f_0 * kd_55[k];

        t_56[k] = -3.0 * hd_56[k]
                  + f_0 * kd_56[k];

        t_57[k] = -3.0 * hd_57[k]
                  + f_0 * kd_57[k];

        t_58[k] = -3.0 * hd_58[k]
                  + f_0 * kd_58[k];

        t_59[k] = -3.0 * hd_59[k]
                  + f_0 * kd_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, hd_60, hd_61, hd_62, hd_63, hd_64, \
                         kd_60, kd_61, kd_62, kd_63, kd_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -2.0 * hd_60[k]
                  + f_0 * kd_60[k];

        t_61[k] = -2.0 * hd_61[k]
                  + f_0 * kd_61[k];

        t_62[k] = -2.0 * hd_62[k]
                  + f_0 * kd_62[k];

        t_63[k] = -2.0 * hd_63[k]
                  + f_0 * kd_63[k];

        t_64[k] = -2.0 * hd_64[k]
                  + f_0 * kd_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, hd_65, hd_66, hd_67, hd_68, hd_69, \
                         kd_65, kd_66, kd_67, kd_68, kd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -2.0 * hd_65[k]
                  + f_0 * kd_65[k];

        t_66[k] = -2.0 * hd_66[k]
                  + f_0 * kd_66[k];

        t_67[k] = -2.0 * hd_67[k]
                  + f_0 * kd_67[k];

        t_68[k] = -2.0 * hd_68[k]
                  + f_0 * kd_68[k];

        t_69[k] = -2.0 * hd_69[k]
                  + f_0 * kd_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, hd_70, hd_71, hd_72, hd_73, hd_74, \
                         kd_70, kd_71, kd_72, kd_73, kd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * hd_70[k]
                  + f_0 * kd_70[k];

        t_71[k] = -2.0 * hd_71[k]
                  + f_0 * kd_71[k];

        t_72[k] = -2.0 * hd_72[k]
                  + f_0 * kd_72[k];

        t_73[k] = -2.0 * hd_73[k]
                  + f_0 * kd_73[k];

        t_74[k] = -2.0 * hd_74[k]
                  + f_0 * kd_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, hd_75, hd_76, hd_77, hd_78, hd_79, \
                         kd_75, kd_76, kd_77, kd_78, kd_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * hd_75[k]
                  + f_0 * kd_75[k];

        t_76[k] = -2.0 * hd_76[k]
                  + f_0 * kd_76[k];

        t_77[k] = -2.0 * hd_77[k]
                  + f_0 * kd_77[k];

        t_78[k] = -2.0 * hd_78[k]
                  + f_0 * kd_78[k];

        t_79[k] = -2.0 * hd_79[k]
                  + f_0 * kd_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, hd_80, hd_81, hd_82, hd_83, hd_84, \
                         kd_80, kd_81, kd_82, kd_83, kd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -2.0 * hd_80[k]
                  + f_0 * kd_80[k];

        t_81[k] = -2.0 * hd_81[k]
                  + f_0 * kd_81[k];

        t_82[k] = -2.0 * hd_82[k]
                  + f_0 * kd_82[k];

        t_83[k] = -2.0 * hd_83[k]
                  + f_0 * kd_83[k];

        t_84[k] = -2.0 * hd_84[k]
                  + f_0 * kd_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, hd_85, hd_86, hd_87, hd_88, hd_89, \
                         kd_85, kd_86, kd_87, kd_88, kd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * hd_85[k]
                  + f_0 * kd_85[k];

        t_86[k] = -2.0 * hd_86[k]
                  + f_0 * kd_86[k];

        t_87[k] = -2.0 * hd_87[k]
                  + f_0 * kd_87[k];

        t_88[k] = -2.0 * hd_88[k]
                  + f_0 * kd_88[k];

        t_89[k] = -2.0 * hd_89[k]
                  + f_0 * kd_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, hd_90, hd_91, hd_92, hd_93, hd_94, \
                         kd_90, kd_91, kd_92, kd_93, kd_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -hd_90[k]
                  + f_0 * kd_90[k];

        t_91[k] = -hd_91[k]
                  + f_0 * kd_91[k];

        t_92[k] = -hd_92[k]
                  + f_0 * kd_92[k];

        t_93[k] = -hd_93[k]
                  + f_0 * kd_93[k];

        t_94[k] = -hd_94[k]
                  + f_0 * kd_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, hd_95, hd_96, hd_97, hd_98, hd_99, \
                         kd_95, kd_96, kd_97, kd_98, kd_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -hd_95[k]
                  + f_0 * kd_95[k];

        t_96[k] = -hd_96[k]
                  + f_0 * kd_96[k];

        t_97[k] = -hd_97[k]
                  + f_0 * kd_97[k];

        t_98[k] = -hd_98[k]
                  + f_0 * kd_98[k];

        t_99[k] = -hd_99[k]
                  + f_0 * kd_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, hd_100, hd_101, hd_102, hd_103, \
                         hd_104, kd_100, kd_101, kd_102, kd_103, \
                         kd_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -hd_100[k]
                   + f_0 * kd_100[k];

        t_101[k] = -hd_101[k]
                   + f_0 * kd_101[k];

        t_102[k] = -hd_102[k]
                   + f_0 * kd_102[k];

        t_103[k] = -hd_103[k]
                   + f_0 * kd_103[k];

        t_104[k] = -hd_104[k]
                   + f_0 * kd_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, hd_105, hd_106, hd_107, hd_108, \
                         hd_109, kd_105, kd_106, kd_107, kd_108, \
                         kd_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -hd_105[k]
                   + f_0 * kd_105[k];

        t_106[k] = -hd_106[k]
                   + f_0 * kd_106[k];

        t_107[k] = -hd_107[k]
                   + f_0 * kd_107[k];

        t_108[k] = -hd_108[k]
                   + f_0 * kd_108[k];

        t_109[k] = -hd_109[k]
                   + f_0 * kd_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, hd_110, hd_111, hd_112, hd_113, \
                         hd_114, kd_110, kd_111, kd_112, kd_113, \
                         kd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -hd_110[k]
                   + f_0 * kd_110[k];

        t_111[k] = -hd_111[k]
                   + f_0 * kd_111[k];

        t_112[k] = -hd_112[k]
                   + f_0 * kd_112[k];

        t_113[k] = -hd_113[k]
                   + f_0 * kd_113[k];

        t_114[k] = -hd_114[k]
                   + f_0 * kd_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, hd_115, hd_116, hd_117, hd_118, \
                         hd_119, kd_115, kd_116, kd_117, kd_118, \
                         kd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -hd_115[k]
                   + f_0 * kd_115[k];

        t_116[k] = -hd_116[k]
                   + f_0 * kd_116[k];

        t_117[k] = -hd_117[k]
                   + f_0 * kd_117[k];

        t_118[k] = -hd_118[k]
                   + f_0 * kd_118[k];

        t_119[k] = -hd_119[k]
                   + f_0 * kd_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, hd_120, hd_121, hd_122, hd_123, \
                         hd_124, kd_120, kd_121, kd_122, kd_123, \
                         kd_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -hd_120[k]
                   + f_0 * kd_120[k];

        t_121[k] = -hd_121[k]
                   + f_0 * kd_121[k];

        t_122[k] = -hd_122[k]
                   + f_0 * kd_122[k];

        t_123[k] = -hd_123[k]
                   + f_0 * kd_123[k];

        t_124[k] = -hd_124[k]
                   + f_0 * kd_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, t_131, hd_125, kd_125, \
                         kd_126, kd_127, kd_128, kd_129, kd_130, \
                         kd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -hd_125[k]
                   + f_0 * kd_125[k];

        t_126[k] = f_0 * kd_126[k];

        t_127[k] = f_0 * kd_127[k];

        t_128[k] = f_0 * kd_128[k];

        t_129[k] = f_0 * kd_129[k];

        t_130[k] = f_0 * kd_130[k];

        t_131[k] = f_0 * kd_131[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, t_137, t_138, t_139, kd_132, \
                         kd_133, kd_134, kd_135, kd_136, kd_137, kd_138, \
                         kd_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_0 * kd_132[k];

        t_133[k] = f_0 * kd_133[k];

        t_134[k] = f_0 * kd_134[k];

        t_135[k] = f_0 * kd_135[k];

        t_136[k] = f_0 * kd_136[k];

        t_137[k] = f_0 * kd_137[k];

        t_138[k] = f_0 * kd_138[k];

        t_139[k] = f_0 * kd_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, t_146, t_147, kd_140, \
                         kd_141, kd_142, kd_143, kd_144, kd_145, kd_146, \
                         kd_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * kd_140[k];

        t_141[k] = f_0 * kd_141[k];

        t_142[k] = f_0 * kd_142[k];

        t_143[k] = f_0 * kd_143[k];

        t_144[k] = f_0 * kd_144[k];

        t_145[k] = f_0 * kd_145[k];

        t_146[k] = f_0 * kd_146[k];

        t_147[k] = f_0 * kd_147[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, t_153, t_154, t_155, kd_148, \
                         kd_149, kd_150, kd_151, kd_152, kd_153, kd_154, \
                         kd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * kd_148[k];

        t_149[k] = f_0 * kd_149[k];

        t_150[k] = f_0 * kd_150[k];

        t_151[k] = f_0 * kd_151[k];

        t_152[k] = f_0 * kd_152[k];

        t_153[k] = f_0 * kd_153[k];

        t_154[k] = f_0 * kd_154[k];

        t_155[k] = f_0 * kd_155[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, t_162, t_163, kd_156, \
                         kd_157, kd_158, kd_159, kd_160, kd_161, kd_162, \
                         kd_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * kd_156[k];

        t_157[k] = f_0 * kd_157[k];

        t_158[k] = f_0 * kd_158[k];

        t_159[k] = f_0 * kd_159[k];

        t_160[k] = f_0 * kd_160[k];

        t_161[k] = f_0 * kd_161[k];

        t_162[k] = f_0 * kd_162[k];

        t_163[k] = f_0 * kd_163[k];
    }
}

static auto
compute_prim_geom_10_id_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kd, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);

    const auto *kd_164 = buffer.data(kd + 164);
    const auto *kd_165 = buffer.data(kd + 165);
    const auto *kd_166 = buffer.data(kd + 166);
    const auto *kd_167 = buffer.data(kd + 167);

#pragma omp simd aligned(t_164, t_165, t_166, t_167, kd_164, kd_165, kd_166, \
                         kd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_0 * kd_164[k];

        t_165[k] = f_0 * kd_165[k];

        t_166[k] = f_0 * kd_166[k];

        t_167[k] = f_0 * kd_167[k];
    }
}

auto
compute_prim_geom_10_id_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t hd, const size_t kd,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_id_electron_repulsion_0_piece0(buffer, target, hd, kd, ncols, alpha);

    compute_prim_geom_10_id_electron_repulsion_0_piece1(buffer, target, kd, ncols, alpha);
}

static auto
compute_prim_geom_10_id_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hd, const size_t kd,
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

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_73 = buffer.data(hd + 73);
    const auto *hd_74 = buffer.data(hd + 74);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_79 = buffer.data(hd + 79);
    const auto *hd_80 = buffer.data(hd + 80);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_83 = buffer.data(hd + 83);
    const auto *hd_84 = buffer.data(hd + 84);
    const auto *hd_85 = buffer.data(hd + 85);
    const auto *hd_86 = buffer.data(hd + 86);
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_88 = buffer.data(hd + 88);
    const auto *hd_89 = buffer.data(hd + 89);
    const auto *hd_90 = buffer.data(hd + 90);
    const auto *hd_91 = buffer.data(hd + 91);
    const auto *hd_92 = buffer.data(hd + 92);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_94 = buffer.data(hd + 94);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_96 = buffer.data(hd + 96);
    const auto *hd_97 = buffer.data(hd + 97);
    const auto *hd_98 = buffer.data(hd + 98);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_102 = buffer.data(hd + 102);
    const auto *hd_103 = buffer.data(hd + 103);
    const auto *hd_104 = buffer.data(hd + 104);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_108 = buffer.data(hd + 108);
    const auto *hd_109 = buffer.data(hd + 109);
    const auto *hd_110 = buffer.data(hd + 110);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_114 = buffer.data(hd + 114);
    const auto *hd_115 = buffer.data(hd + 115);
    const auto *hd_116 = buffer.data(hd + 116);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_119 = buffer.data(hd + 119);
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_121 = buffer.data(hd + 121);
    const auto *hd_122 = buffer.data(hd + 122);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_124 = buffer.data(hd + 124);
    const auto *hd_125 = buffer.data(hd + 125);

    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_91 = buffer.data(kd + 91);
    const auto *kd_92 = buffer.data(kd + 92);
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_94 = buffer.data(kd + 94);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_97 = buffer.data(kd + 97);
    const auto *kd_98 = buffer.data(kd + 98);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_103 = buffer.data(kd + 103);
    const auto *kd_104 = buffer.data(kd + 104);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_109 = buffer.data(kd + 109);
    const auto *kd_110 = buffer.data(kd + 110);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_115 = buffer.data(kd + 115);
    const auto *kd_116 = buffer.data(kd + 116);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_128 = buffer.data(kd + 128);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_130 = buffer.data(kd + 130);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_133 = buffer.data(kd + 133);
    const auto *kd_134 = buffer.data(kd + 134);
    const auto *kd_135 = buffer.data(kd + 135);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_139 = buffer.data(kd + 139);
    const auto *kd_140 = buffer.data(kd + 140);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_145 = buffer.data(kd + 145);
    const auto *kd_146 = buffer.data(kd + 146);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_151 = buffer.data(kd + 151);
    const auto *kd_152 = buffer.data(kd + 152);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_157 = buffer.data(kd + 157);
    const auto *kd_158 = buffer.data(kd + 158);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_161 = buffer.data(kd + 161);
    const auto *kd_168 = buffer.data(kd + 168);
    const auto *kd_169 = buffer.data(kd + 169);
    const auto *kd_170 = buffer.data(kd + 170);
    const auto *kd_171 = buffer.data(kd + 171);
    const auto *kd_172 = buffer.data(kd + 172);
    const auto *kd_173 = buffer.data(kd + 173);
    const auto *kd_174 = buffer.data(kd + 174);
    const auto *kd_175 = buffer.data(kd + 175);
    const auto *kd_176 = buffer.data(kd + 176);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_180 = buffer.data(kd + 180);
    const auto *kd_181 = buffer.data(kd + 181);
    const auto *kd_182 = buffer.data(kd + 182);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_186 = buffer.data(kd + 186);
    const auto *kd_187 = buffer.data(kd + 187);
    const auto *kd_188 = buffer.data(kd + 188);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_192 = buffer.data(kd + 192);
    const auto *kd_193 = buffer.data(kd + 193);
    const auto *kd_194 = buffer.data(kd + 194);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_198 = buffer.data(kd + 198);
    const auto *kd_199 = buffer.data(kd + 199);
    const auto *kd_200 = buffer.data(kd + 200);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, hd_0, kd_6, kd_7, kd_8, kd_9, \
                         kd_10, kd_11, kd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_6[k];

        t_1[k] = f_0 * kd_7[k];

        t_2[k] = f_0 * kd_8[k];

        t_3[k] = f_0 * kd_9[k];

        t_4[k] = f_0 * kd_10[k];

        t_5[k] = f_0 * kd_11[k];

        t_6[k] = -hd_0[k]
                 + f_0 * kd_18[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, hd_1, hd_2, hd_3, hd_4, hd_5, kd_19, \
                         kd_20, kd_21, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -hd_1[k]
                 + f_0 * kd_19[k];

        t_8[k] = -hd_2[k]
                 + f_0 * kd_20[k];

        t_9[k] = -hd_3[k]
                 + f_0 * kd_21[k];

        t_10[k] = -hd_4[k]
                  + f_0 * kd_22[k];

        t_11[k] = -hd_5[k]
                  + f_0 * kd_23[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, t_18, hd_6, kd_24, kd_25, kd_26, \
                         kd_27, kd_28, kd_29, kd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * kd_24[k];

        t_13[k] = f_0 * kd_25[k];

        t_14[k] = f_0 * kd_26[k];

        t_15[k] = f_0 * kd_27[k];

        t_16[k] = f_0 * kd_28[k];

        t_17[k] = f_0 * kd_29[k];

        t_18[k] = -2.0 * hd_6[k]
                  + f_0 * kd_36[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, hd_7, hd_8, hd_9, hd_10, hd_11, kd_37, \
                         kd_38, kd_39, kd_40, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -2.0 * hd_7[k]
                  + f_0 * kd_37[k];

        t_20[k] = -2.0 * hd_8[k]
                  + f_0 * kd_38[k];

        t_21[k] = -2.0 * hd_9[k]
                  + f_0 * kd_39[k];

        t_22[k] = -2.0 * hd_10[k]
                  + f_0 * kd_40[k];

        t_23[k] = -2.0 * hd_11[k]
                  + f_0 * kd_41[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, hd_12, hd_13, hd_14, hd_15, hd_16, \
                         kd_42, kd_43, kd_44, kd_45, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -hd_12[k]
                  + f_0 * kd_42[k];

        t_25[k] = -hd_13[k]
                  + f_0 * kd_43[k];

        t_26[k] = -hd_14[k]
                  + f_0 * kd_44[k];

        t_27[k] = -hd_15[k]
                  + f_0 * kd_45[k];

        t_28[k] = -hd_16[k]
                  + f_0 * kd_46[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, t_35, hd_17, kd_47, kd_48, kd_49, \
                         kd_50, kd_51, kd_52, kd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -hd_17[k]
                  + f_0 * kd_47[k];

        t_30[k] = f_0 * kd_48[k];

        t_31[k] = f_0 * kd_49[k];

        t_32[k] = f_0 * kd_50[k];

        t_33[k] = f_0 * kd_51[k];

        t_34[k] = f_0 * kd_52[k];

        t_35[k] = f_0 * kd_53[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, hd_18, hd_19, hd_20, hd_21, hd_22, \
                         kd_60, kd_61, kd_62, kd_63, kd_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -3.0 * hd_18[k]
                  + f_0 * kd_60[k];

        t_37[k] = -3.0 * hd_19[k]
                  + f_0 * kd_61[k];

        t_38[k] = -3.0 * hd_20[k]
                  + f_0 * kd_62[k];

        t_39[k] = -3.0 * hd_21[k]
                  + f_0 * kd_63[k];

        t_40[k] = -3.0 * hd_22[k]
                  + f_0 * kd_64[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, hd_23, hd_24, hd_25, hd_26, hd_27, \
                         kd_65, kd_66, kd_67, kd_68, kd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -3.0 * hd_23[k]
                  + f_0 * kd_65[k];

        t_42[k] = -2.0 * hd_24[k]
                  + f_0 * kd_66[k];

        t_43[k] = -2.0 * hd_25[k]
                  + f_0 * kd_67[k];

        t_44[k] = -2.0 * hd_26[k]
                  + f_0 * kd_68[k];

        t_45[k] = -2.0 * hd_27[k]
                  + f_0 * kd_69[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, hd_28, hd_29, hd_30, hd_31, hd_32, \
                         kd_70, kd_71, kd_72, kd_73, kd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -2.0 * hd_28[k]
                  + f_0 * kd_70[k];

        t_47[k] = -2.0 * hd_29[k]
                  + f_0 * kd_71[k];

        t_48[k] = -hd_30[k]
                  + f_0 * kd_72[k];

        t_49[k] = -hd_31[k]
                  + f_0 * kd_73[k];

        t_50[k] = -hd_32[k]
                  + f_0 * kd_74[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, hd_33, hd_34, hd_35, kd_75, \
                         kd_76, kd_77, kd_78, kd_79, kd_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -hd_33[k]
                  + f_0 * kd_75[k];

        t_52[k] = -hd_34[k]
                  + f_0 * kd_76[k];

        t_53[k] = -hd_35[k]
                  + f_0 * kd_77[k];

        t_54[k] = f_0 * kd_78[k];

        t_55[k] = f_0 * kd_79[k];

        t_56[k] = f_0 * kd_80[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, t_62, hd_36, hd_37, hd_38, kd_81, \
                         kd_82, kd_83, kd_90, kd_91, kd_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_0 * kd_81[k];

        t_58[k] = f_0 * kd_82[k];

        t_59[k] = f_0 * kd_83[k];

        t_60[k] = -4.0 * hd_36[k]
                  + f_0 * kd_90[k];

        t_61[k] = -4.0 * hd_37[k]
                  + f_0 * kd_91[k];

        t_62[k] = -4.0 * hd_38[k]
                  + f_0 * kd_92[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, hd_39, hd_40, hd_41, hd_42, hd_43, \
                         kd_93, kd_94, kd_95, kd_96, kd_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -4.0 * hd_39[k]
                  + f_0 * kd_93[k];

        t_64[k] = -4.0 * hd_40[k]
                  + f_0 * kd_94[k];

        t_65[k] = -4.0 * hd_41[k]
                  + f_0 * kd_95[k];

        t_66[k] = -3.0 * hd_42[k]
                  + f_0 * kd_96[k];

        t_67[k] = -3.0 * hd_43[k]
                  + f_0 * kd_97[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, hd_44, hd_45, hd_46, hd_47, hd_48, \
                         kd_98, kd_99, kd_100, kd_101, kd_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -3.0 * hd_44[k]
                  + f_0 * kd_98[k];

        t_69[k] = -3.0 * hd_45[k]
                  + f_0 * kd_99[k];

        t_70[k] = -3.0 * hd_46[k]
                  + f_0 * kd_100[k];

        t_71[k] = -3.0 * hd_47[k]
                  + f_0 * kd_101[k];

        t_72[k] = -2.0 * hd_48[k]
                  + f_0 * kd_102[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, hd_49, hd_50, hd_51, hd_52, hd_53, \
                         kd_103, kd_104, kd_105, kd_106, kd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = -2.0 * hd_49[k]
                  + f_0 * kd_103[k];

        t_74[k] = -2.0 * hd_50[k]
                  + f_0 * kd_104[k];

        t_75[k] = -2.0 * hd_51[k]
                  + f_0 * kd_105[k];

        t_76[k] = -2.0 * hd_52[k]
                  + f_0 * kd_106[k];

        t_77[k] = -2.0 * hd_53[k]
                  + f_0 * kd_107[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, hd_54, hd_55, hd_56, hd_57, hd_58, \
                         kd_108, kd_109, kd_110, kd_111, kd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = -hd_54[k]
                  + f_0 * kd_108[k];

        t_79[k] = -hd_55[k]
                  + f_0 * kd_109[k];

        t_80[k] = -hd_56[k]
                  + f_0 * kd_110[k];

        t_81[k] = -hd_57[k]
                  + f_0 * kd_111[k];

        t_82[k] = -hd_58[k]
                  + f_0 * kd_112[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, hd_59, kd_113, kd_114, \
                         kd_115, kd_116, kd_117, kd_118, kd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = -hd_59[k]
                  + f_0 * kd_113[k];

        t_84[k] = f_0 * kd_114[k];

        t_85[k] = f_0 * kd_115[k];

        t_86[k] = f_0 * kd_116[k];

        t_87[k] = f_0 * kd_117[k];

        t_88[k] = f_0 * kd_118[k];

        t_89[k] = f_0 * kd_119[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, hd_60, hd_61, hd_62, hd_63, hd_64, \
                         kd_126, kd_127, kd_128, kd_129, kd_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -5.0 * hd_60[k]
                  + f_0 * kd_126[k];

        t_91[k] = -5.0 * hd_61[k]
                  + f_0 * kd_127[k];

        t_92[k] = -5.0 * hd_62[k]
                  + f_0 * kd_128[k];

        t_93[k] = -5.0 * hd_63[k]
                  + f_0 * kd_129[k];

        t_94[k] = -5.0 * hd_64[k]
                  + f_0 * kd_130[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, hd_65, hd_66, hd_67, hd_68, hd_69, \
                         kd_131, kd_132, kd_133, kd_134, kd_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -5.0 * hd_65[k]
                  + f_0 * kd_131[k];

        t_96[k] = -4.0 * hd_66[k]
                  + f_0 * kd_132[k];

        t_97[k] = -4.0 * hd_67[k]
                  + f_0 * kd_133[k];

        t_98[k] = -4.0 * hd_68[k]
                  + f_0 * kd_134[k];

        t_99[k] = -4.0 * hd_69[k]
                  + f_0 * kd_135[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, hd_70, hd_71, hd_72, hd_73, hd_74, \
                         kd_136, kd_137, kd_138, kd_139, kd_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -4.0 * hd_70[k]
                   + f_0 * kd_136[k];

        t_101[k] = -4.0 * hd_71[k]
                   + f_0 * kd_137[k];

        t_102[k] = -3.0 * hd_72[k]
                   + f_0 * kd_138[k];

        t_103[k] = -3.0 * hd_73[k]
                   + f_0 * kd_139[k];

        t_104[k] = -3.0 * hd_74[k]
                   + f_0 * kd_140[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, hd_75, hd_76, hd_77, hd_78, hd_79, \
                         kd_141, kd_142, kd_143, kd_144, kd_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -3.0 * hd_75[k]
                   + f_0 * kd_141[k];

        t_106[k] = -3.0 * hd_76[k]
                   + f_0 * kd_142[k];

        t_107[k] = -3.0 * hd_77[k]
                   + f_0 * kd_143[k];

        t_108[k] = -2.0 * hd_78[k]
                   + f_0 * kd_144[k];

        t_109[k] = -2.0 * hd_79[k]
                   + f_0 * kd_145[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, hd_80, hd_81, hd_82, hd_83, hd_84, \
                         kd_146, kd_147, kd_148, kd_149, kd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * hd_80[k]
                   + f_0 * kd_146[k];

        t_111[k] = -2.0 * hd_81[k]
                   + f_0 * kd_147[k];

        t_112[k] = -2.0 * hd_82[k]
                   + f_0 * kd_148[k];

        t_113[k] = -2.0 * hd_83[k]
                   + f_0 * kd_149[k];

        t_114[k] = -hd_84[k]
                   + f_0 * kd_150[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, hd_85, hd_86, hd_87, hd_88, hd_89, \
                         kd_151, kd_152, kd_153, kd_154, kd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -hd_85[k]
                   + f_0 * kd_151[k];

        t_116[k] = -hd_86[k]
                   + f_0 * kd_152[k];

        t_117[k] = -hd_87[k]
                   + f_0 * kd_153[k];

        t_118[k] = -hd_88[k]
                   + f_0 * kd_154[k];

        t_119[k] = -hd_89[k]
                   + f_0 * kd_155[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, t_126, hd_90, kd_156, \
                         kd_157, kd_158, kd_159, kd_160, kd_161, \
                         kd_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_0 * kd_156[k];

        t_121[k] = f_0 * kd_157[k];

        t_122[k] = f_0 * kd_158[k];

        t_123[k] = f_0 * kd_159[k];

        t_124[k] = f_0 * kd_160[k];

        t_125[k] = f_0 * kd_161[k];

        t_126[k] = -6.0 * hd_90[k]
                   + f_0 * kd_168[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, hd_91, hd_92, hd_93, hd_94, hd_95, \
                         kd_169, kd_170, kd_171, kd_172, kd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -6.0 * hd_91[k]
                   + f_0 * kd_169[k];

        t_128[k] = -6.0 * hd_92[k]
                   + f_0 * kd_170[k];

        t_129[k] = -6.0 * hd_93[k]
                   + f_0 * kd_171[k];

        t_130[k] = -6.0 * hd_94[k]
                   + f_0 * kd_172[k];

        t_131[k] = -6.0 * hd_95[k]
                   + f_0 * kd_173[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, hd_96, hd_97, hd_98, hd_99, \
                         hd_100, kd_174, kd_175, kd_176, kd_177, \
                         kd_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -5.0 * hd_96[k]
                   + f_0 * kd_174[k];

        t_133[k] = -5.0 * hd_97[k]
                   + f_0 * kd_175[k];

        t_134[k] = -5.0 * hd_98[k]
                   + f_0 * kd_176[k];

        t_135[k] = -5.0 * hd_99[k]
                   + f_0 * kd_177[k];

        t_136[k] = -5.0 * hd_100[k]
                   + f_0 * kd_178[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, hd_101, hd_102, hd_103, hd_104, \
                         hd_105, kd_179, kd_180, kd_181, kd_182, \
                         kd_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -5.0 * hd_101[k]
                   + f_0 * kd_179[k];

        t_138[k] = -4.0 * hd_102[k]
                   + f_0 * kd_180[k];

        t_139[k] = -4.0 * hd_103[k]
                   + f_0 * kd_181[k];

        t_140[k] = -4.0 * hd_104[k]
                   + f_0 * kd_182[k];

        t_141[k] = -4.0 * hd_105[k]
                   + f_0 * kd_183[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, hd_106, hd_107, hd_108, hd_109, \
                         hd_110, kd_184, kd_185, kd_186, kd_187, \
                         kd_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -4.0 * hd_106[k]
                   + f_0 * kd_184[k];

        t_143[k] = -4.0 * hd_107[k]
                   + f_0 * kd_185[k];

        t_144[k] = -3.0 * hd_108[k]
                   + f_0 * kd_186[k];

        t_145[k] = -3.0 * hd_109[k]
                   + f_0 * kd_187[k];

        t_146[k] = -3.0 * hd_110[k]
                   + f_0 * kd_188[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, hd_111, hd_112, hd_113, hd_114, \
                         hd_115, kd_189, kd_190, kd_191, kd_192, \
                         kd_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -3.0 * hd_111[k]
                   + f_0 * kd_189[k];

        t_148[k] = -3.0 * hd_112[k]
                   + f_0 * kd_190[k];

        t_149[k] = -3.0 * hd_113[k]
                   + f_0 * kd_191[k];

        t_150[k] = -2.0 * hd_114[k]
                   + f_0 * kd_192[k];

        t_151[k] = -2.0 * hd_115[k]
                   + f_0 * kd_193[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, hd_116, hd_117, hd_118, hd_119, \
                         hd_120, kd_194, kd_195, kd_196, kd_197, \
                         kd_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -2.0 * hd_116[k]
                   + f_0 * kd_194[k];

        t_153[k] = -2.0 * hd_117[k]
                   + f_0 * kd_195[k];

        t_154[k] = -2.0 * hd_118[k]
                   + f_0 * kd_196[k];

        t_155[k] = -2.0 * hd_119[k]
                   + f_0 * kd_197[k];

        t_156[k] = -hd_120[k]
                   + f_0 * kd_198[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, hd_121, hd_122, hd_123, hd_124, \
                         hd_125, kd_199, kd_200, kd_201, kd_202, \
                         kd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -hd_121[k]
                   + f_0 * kd_199[k];

        t_158[k] = -hd_122[k]
                   + f_0 * kd_200[k];

        t_159[k] = -hd_123[k]
                   + f_0 * kd_201[k];

        t_160[k] = -hd_124[k]
                   + f_0 * kd_202[k];

        t_161[k] = -hd_125[k]
                   + f_0 * kd_203[k];
    }
}

static auto
compute_prim_geom_10_id_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kd, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);

    const auto *kd_204 = buffer.data(kd + 204);
    const auto *kd_205 = buffer.data(kd + 205);
    const auto *kd_206 = buffer.data(kd + 206);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, t_167, kd_204, kd_205, kd_206, \
                         kd_207, kd_208, kd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_0 * kd_204[k];

        t_163[k] = f_0 * kd_205[k];

        t_164[k] = f_0 * kd_206[k];

        t_165[k] = f_0 * kd_207[k];

        t_166[k] = f_0 * kd_208[k];

        t_167[k] = f_0 * kd_209[k];
    }
}

auto
compute_prim_geom_10_id_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t hd, const size_t kd,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_id_electron_repulsion_1_piece0(buffer, target, hd, kd, ncols, alpha);

    compute_prim_geom_10_id_electron_repulsion_1_piece1(buffer, target, kd, ncols, alpha);
}

static auto
compute_prim_geom_10_id_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hd, const size_t kd,
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

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_73 = buffer.data(hd + 73);
    const auto *hd_74 = buffer.data(hd + 74);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_79 = buffer.data(hd + 79);
    const auto *hd_80 = buffer.data(hd + 80);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_83 = buffer.data(hd + 83);
    const auto *hd_84 = buffer.data(hd + 84);
    const auto *hd_85 = buffer.data(hd + 85);
    const auto *hd_86 = buffer.data(hd + 86);
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_88 = buffer.data(hd + 88);
    const auto *hd_89 = buffer.data(hd + 89);
    const auto *hd_90 = buffer.data(hd + 90);
    const auto *hd_91 = buffer.data(hd + 91);
    const auto *hd_92 = buffer.data(hd + 92);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_94 = buffer.data(hd + 94);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_96 = buffer.data(hd + 96);
    const auto *hd_97 = buffer.data(hd + 97);
    const auto *hd_98 = buffer.data(hd + 98);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_102 = buffer.data(hd + 102);
    const auto *hd_103 = buffer.data(hd + 103);
    const auto *hd_104 = buffer.data(hd + 104);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_108 = buffer.data(hd + 108);
    const auto *hd_109 = buffer.data(hd + 109);
    const auto *hd_110 = buffer.data(hd + 110);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_114 = buffer.data(hd + 114);
    const auto *hd_115 = buffer.data(hd + 115);
    const auto *hd_116 = buffer.data(hd + 116);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_119 = buffer.data(hd + 119);
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_121 = buffer.data(hd + 121);

    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_97 = buffer.data(kd + 97);
    const auto *kd_98 = buffer.data(kd + 98);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_103 = buffer.data(kd + 103);
    const auto *kd_104 = buffer.data(kd + 104);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_109 = buffer.data(kd + 109);
    const auto *kd_110 = buffer.data(kd + 110);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_115 = buffer.data(kd + 115);
    const auto *kd_116 = buffer.data(kd + 116);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_121 = buffer.data(kd + 121);
    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_124 = buffer.data(kd + 124);
    const auto *kd_125 = buffer.data(kd + 125);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_133 = buffer.data(kd + 133);
    const auto *kd_134 = buffer.data(kd + 134);
    const auto *kd_135 = buffer.data(kd + 135);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_139 = buffer.data(kd + 139);
    const auto *kd_140 = buffer.data(kd + 140);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_145 = buffer.data(kd + 145);
    const auto *kd_146 = buffer.data(kd + 146);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_151 = buffer.data(kd + 151);
    const auto *kd_152 = buffer.data(kd + 152);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_157 = buffer.data(kd + 157);
    const auto *kd_158 = buffer.data(kd + 158);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_161 = buffer.data(kd + 161);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_163 = buffer.data(kd + 163);
    const auto *kd_164 = buffer.data(kd + 164);
    const auto *kd_165 = buffer.data(kd + 165);
    const auto *kd_166 = buffer.data(kd + 166);
    const auto *kd_167 = buffer.data(kd + 167);
    const auto *kd_174 = buffer.data(kd + 174);
    const auto *kd_175 = buffer.data(kd + 175);
    const auto *kd_176 = buffer.data(kd + 176);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_180 = buffer.data(kd + 180);
    const auto *kd_181 = buffer.data(kd + 181);
    const auto *kd_182 = buffer.data(kd + 182);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_186 = buffer.data(kd + 186);
    const auto *kd_187 = buffer.data(kd + 187);
    const auto *kd_188 = buffer.data(kd + 188);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_192 = buffer.data(kd + 192);
    const auto *kd_193 = buffer.data(kd + 193);
    const auto *kd_194 = buffer.data(kd + 194);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_198 = buffer.data(kd + 198);
    const auto *kd_199 = buffer.data(kd + 199);
    const auto *kd_200 = buffer.data(kd + 200);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);
    const auto *kd_204 = buffer.data(kd + 204);
    const auto *kd_205 = buffer.data(kd + 205);
    const auto *kd_206 = buffer.data(kd + 206);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);
    const auto *kd_210 = buffer.data(kd + 210);
    const auto *kd_211 = buffer.data(kd + 211);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, kd_12, kd_13, kd_14, kd_15, \
                         kd_16, kd_17, kd_24, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kd_12[k];

        t_1[k] = f_0 * kd_13[k];

        t_2[k] = f_0 * kd_14[k];

        t_3[k] = f_0 * kd_15[k];

        t_4[k] = f_0 * kd_16[k];

        t_5[k] = f_0 * kd_17[k];

        t_6[k] = f_0 * kd_24[k];

        t_7[k] = f_0 * kd_25[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, hd_0, hd_1, kd_26, kd_27, kd_28, \
                         kd_29, kd_30, kd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * kd_26[k];

        t_9[k] = f_0 * kd_27[k];

        t_10[k] = f_0 * kd_28[k];

        t_11[k] = f_0 * kd_29[k];

        t_12[k] = -hd_0[k]
                  + f_0 * kd_30[k];

        t_13[k] = -hd_1[k]
                  + f_0 * kd_31[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, hd_2, hd_3, hd_4, hd_5, kd_32, \
                         kd_33, kd_34, kd_35, kd_42, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -hd_2[k]
                  + f_0 * kd_32[k];

        t_15[k] = -hd_3[k]
                  + f_0 * kd_33[k];

        t_16[k] = -hd_4[k]
                  + f_0 * kd_34[k];

        t_17[k] = -hd_5[k]
                  + f_0 * kd_35[k];

        t_18[k] = f_0 * kd_42[k];

        t_19[k] = f_0 * kd_43[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, hd_6, hd_7, kd_44, kd_45, kd_46, \
                         kd_47, kd_48, kd_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * kd_44[k];

        t_21[k] = f_0 * kd_45[k];

        t_22[k] = f_0 * kd_46[k];

        t_23[k] = f_0 * kd_47[k];

        t_24[k] = -hd_6[k]
                  + f_0 * kd_48[k];

        t_25[k] = -hd_7[k]
                  + f_0 * kd_49[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, hd_8, hd_9, hd_10, hd_11, hd_12, kd_50, \
                         kd_51, kd_52, kd_53, kd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -hd_8[k]
                  + f_0 * kd_50[k];

        t_27[k] = -hd_9[k]
                  + f_0 * kd_51[k];

        t_28[k] = -hd_10[k]
                  + f_0 * kd_52[k];

        t_29[k] = -hd_11[k]
                  + f_0 * kd_53[k];

        t_30[k] = -2.0 * hd_12[k]
                  + f_0 * kd_54[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, hd_13, hd_14, hd_15, hd_16, hd_17, \
                         kd_55, kd_56, kd_57, kd_58, kd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -2.0 * hd_13[k]
                  + f_0 * kd_55[k];

        t_32[k] = -2.0 * hd_14[k]
                  + f_0 * kd_56[k];

        t_33[k] = -2.0 * hd_15[k]
                  + f_0 * kd_57[k];

        t_34[k] = -2.0 * hd_16[k]
                  + f_0 * kd_58[k];

        t_35[k] = -2.0 * hd_17[k]
                  + f_0 * kd_59[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, t_42, hd_18, kd_66, kd_67, kd_68, \
                         kd_69, kd_70, kd_71, kd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * kd_66[k];

        t_37[k] = f_0 * kd_67[k];

        t_38[k] = f_0 * kd_68[k];

        t_39[k] = f_0 * kd_69[k];

        t_40[k] = f_0 * kd_70[k];

        t_41[k] = f_0 * kd_71[k];

        t_42[k] = -hd_18[k]
                  + f_0 * kd_72[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, hd_19, hd_20, hd_21, hd_22, hd_23, \
                         kd_73, kd_74, kd_75, kd_76, kd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -hd_19[k]
                  + f_0 * kd_73[k];

        t_44[k] = -hd_20[k]
                  + f_0 * kd_74[k];

        t_45[k] = -hd_21[k]
                  + f_0 * kd_75[k];

        t_46[k] = -hd_22[k]
                  + f_0 * kd_76[k];

        t_47[k] = -hd_23[k]
                  + f_0 * kd_77[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, hd_24, hd_25, hd_26, hd_27, hd_28, \
                         kd_78, kd_79, kd_80, kd_81, kd_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -2.0 * hd_24[k]
                  + f_0 * kd_78[k];

        t_49[k] = -2.0 * hd_25[k]
                  + f_0 * kd_79[k];

        t_50[k] = -2.0 * hd_26[k]
                  + f_0 * kd_80[k];

        t_51[k] = -2.0 * hd_27[k]
                  + f_0 * kd_81[k];

        t_52[k] = -2.0 * hd_28[k]
                  + f_0 * kd_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, hd_29, hd_30, hd_31, hd_32, hd_33, \
                         kd_83, kd_84, kd_85, kd_86, kd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -2.0 * hd_29[k]
                  + f_0 * kd_83[k];

        t_54[k] = -3.0 * hd_30[k]
                  + f_0 * kd_84[k];

        t_55[k] = -3.0 * hd_31[k]
                  + f_0 * kd_85[k];

        t_56[k] = -3.0 * hd_32[k]
                  + f_0 * kd_86[k];

        t_57[k] = -3.0 * hd_33[k]
                  + f_0 * kd_87[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, t_63, t_64, hd_34, hd_35, kd_88, kd_89, \
                         kd_96, kd_97, kd_98, kd_99, kd_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -3.0 * hd_34[k]
                  + f_0 * kd_88[k];

        t_59[k] = -3.0 * hd_35[k]
                  + f_0 * kd_89[k];

        t_60[k] = f_0 * kd_96[k];

        t_61[k] = f_0 * kd_97[k];

        t_62[k] = f_0 * kd_98[k];

        t_63[k] = f_0 * kd_99[k];

        t_64[k] = f_0 * kd_100[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, hd_36, hd_37, hd_38, hd_39, kd_101, \
                         kd_102, kd_103, kd_104, kd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_0 * kd_101[k];

        t_66[k] = -hd_36[k]
                  + f_0 * kd_102[k];

        t_67[k] = -hd_37[k]
                  + f_0 * kd_103[k];

        t_68[k] = -hd_38[k]
                  + f_0 * kd_104[k];

        t_69[k] = -hd_39[k]
                  + f_0 * kd_105[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, hd_40, hd_41, hd_42, hd_43, hd_44, \
                         kd_106, kd_107, kd_108, kd_109, kd_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -hd_40[k]
                  + f_0 * kd_106[k];

        t_71[k] = -hd_41[k]
                  + f_0 * kd_107[k];

        t_72[k] = -2.0 * hd_42[k]
                  + f_0 * kd_108[k];

        t_73[k] = -2.0 * hd_43[k]
                  + f_0 * kd_109[k];

        t_74[k] = -2.0 * hd_44[k]
                  + f_0 * kd_110[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, hd_45, hd_46, hd_47, hd_48, hd_49, \
                         kd_111, kd_112, kd_113, kd_114, kd_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * hd_45[k]
                  + f_0 * kd_111[k];

        t_76[k] = -2.0 * hd_46[k]
                  + f_0 * kd_112[k];

        t_77[k] = -2.0 * hd_47[k]
                  + f_0 * kd_113[k];

        t_78[k] = -3.0 * hd_48[k]
                  + f_0 * kd_114[k];

        t_79[k] = -3.0 * hd_49[k]
                  + f_0 * kd_115[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, hd_50, hd_51, hd_52, hd_53, hd_54, \
                         kd_116, kd_117, kd_118, kd_119, kd_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * hd_50[k]
                  + f_0 * kd_116[k];

        t_81[k] = -3.0 * hd_51[k]
                  + f_0 * kd_117[k];

        t_82[k] = -3.0 * hd_52[k]
                  + f_0 * kd_118[k];

        t_83[k] = -3.0 * hd_53[k]
                  + f_0 * kd_119[k];

        t_84[k] = -4.0 * hd_54[k]
                  + f_0 * kd_120[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, hd_55, hd_56, hd_57, hd_58, hd_59, \
                         kd_121, kd_122, kd_123, kd_124, kd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * hd_55[k]
                  + f_0 * kd_121[k];

        t_86[k] = -4.0 * hd_56[k]
                  + f_0 * kd_122[k];

        t_87[k] = -4.0 * hd_57[k]
                  + f_0 * kd_123[k];

        t_88[k] = -4.0 * hd_58[k]
                  + f_0 * kd_124[k];

        t_89[k] = -4.0 * hd_59[k]
                  + f_0 * kd_125[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, hd_60, kd_132, kd_133, \
                         kd_134, kd_135, kd_136, kd_137, kd_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * kd_132[k];

        t_91[k] = f_0 * kd_133[k];

        t_92[k] = f_0 * kd_134[k];

        t_93[k] = f_0 * kd_135[k];

        t_94[k] = f_0 * kd_136[k];

        t_95[k] = f_0 * kd_137[k];

        t_96[k] = -hd_60[k]
                  + f_0 * kd_138[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, hd_61, hd_62, hd_63, hd_64, hd_65, \
                         kd_139, kd_140, kd_141, kd_142, kd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -hd_61[k]
                  + f_0 * kd_139[k];

        t_98[k] = -hd_62[k]
                  + f_0 * kd_140[k];

        t_99[k] = -hd_63[k]
                  + f_0 * kd_141[k];

        t_100[k] = -hd_64[k]
                   + f_0 * kd_142[k];

        t_101[k] = -hd_65[k]
                   + f_0 * kd_143[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, hd_66, hd_67, hd_68, hd_69, hd_70, \
                         kd_144, kd_145, kd_146, kd_147, kd_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -2.0 * hd_66[k]
                   + f_0 * kd_144[k];

        t_103[k] = -2.0 * hd_67[k]
                   + f_0 * kd_145[k];

        t_104[k] = -2.0 * hd_68[k]
                   + f_0 * kd_146[k];

        t_105[k] = -2.0 * hd_69[k]
                   + f_0 * kd_147[k];

        t_106[k] = -2.0 * hd_70[k]
                   + f_0 * kd_148[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, hd_71, hd_72, hd_73, hd_74, hd_75, \
                         kd_149, kd_150, kd_151, kd_152, kd_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -2.0 * hd_71[k]
                   + f_0 * kd_149[k];

        t_108[k] = -3.0 * hd_72[k]
                   + f_0 * kd_150[k];

        t_109[k] = -3.0 * hd_73[k]
                   + f_0 * kd_151[k];

        t_110[k] = -3.0 * hd_74[k]
                   + f_0 * kd_152[k];

        t_111[k] = -3.0 * hd_75[k]
                   + f_0 * kd_153[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, hd_76, hd_77, hd_78, hd_79, hd_80, \
                         kd_154, kd_155, kd_156, kd_157, kd_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -3.0 * hd_76[k]
                   + f_0 * kd_154[k];

        t_113[k] = -3.0 * hd_77[k]
                   + f_0 * kd_155[k];

        t_114[k] = -4.0 * hd_78[k]
                   + f_0 * kd_156[k];

        t_115[k] = -4.0 * hd_79[k]
                   + f_0 * kd_157[k];

        t_116[k] = -4.0 * hd_80[k]
                   + f_0 * kd_158[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, hd_81, hd_82, hd_83, hd_84, hd_85, \
                         kd_159, kd_160, kd_161, kd_162, kd_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -4.0 * hd_81[k]
                   + f_0 * kd_159[k];

        t_118[k] = -4.0 * hd_82[k]
                   + f_0 * kd_160[k];

        t_119[k] = -4.0 * hd_83[k]
                   + f_0 * kd_161[k];

        t_120[k] = -5.0 * hd_84[k]
                   + f_0 * kd_162[k];

        t_121[k] = -5.0 * hd_85[k]
                   + f_0 * kd_163[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, hd_86, hd_87, hd_88, hd_89, \
                         kd_164, kd_165, kd_166, kd_167, kd_174, \
                         kd_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -5.0 * hd_86[k]
                   + f_0 * kd_164[k];

        t_123[k] = -5.0 * hd_87[k]
                   + f_0 * kd_165[k];

        t_124[k] = -5.0 * hd_88[k]
                   + f_0 * kd_166[k];

        t_125[k] = -5.0 * hd_89[k]
                   + f_0 * kd_167[k];

        t_126[k] = f_0 * kd_174[k];

        t_127[k] = f_0 * kd_175[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, hd_90, hd_91, kd_176, \
                         kd_177, kd_178, kd_179, kd_180, kd_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_0 * kd_176[k];

        t_129[k] = f_0 * kd_177[k];

        t_130[k] = f_0 * kd_178[k];

        t_131[k] = f_0 * kd_179[k];

        t_132[k] = -hd_90[k]
                   + f_0 * kd_180[k];

        t_133[k] = -hd_91[k]
                   + f_0 * kd_181[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, hd_92, hd_93, hd_94, hd_95, hd_96, \
                         kd_182, kd_183, kd_184, kd_185, kd_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -hd_92[k]
                   + f_0 * kd_182[k];

        t_135[k] = -hd_93[k]
                   + f_0 * kd_183[k];

        t_136[k] = -hd_94[k]
                   + f_0 * kd_184[k];

        t_137[k] = -hd_95[k]
                   + f_0 * kd_185[k];

        t_138[k] = -2.0 * hd_96[k]
                   + f_0 * kd_186[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, hd_97, hd_98, hd_99, hd_100, \
                         hd_101, kd_187, kd_188, kd_189, kd_190, \
                         kd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -2.0 * hd_97[k]
                   + f_0 * kd_187[k];

        t_140[k] = -2.0 * hd_98[k]
                   + f_0 * kd_188[k];

        t_141[k] = -2.0 * hd_99[k]
                   + f_0 * kd_189[k];

        t_142[k] = -2.0 * hd_100[k]
                   + f_0 * kd_190[k];

        t_143[k] = -2.0 * hd_101[k]
                   + f_0 * kd_191[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, hd_102, hd_103, hd_104, hd_105, \
                         hd_106, kd_192, kd_193, kd_194, kd_195, \
                         kd_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -3.0 * hd_102[k]
                   + f_0 * kd_192[k];

        t_145[k] = -3.0 * hd_103[k]
                   + f_0 * kd_193[k];

        t_146[k] = -3.0 * hd_104[k]
                   + f_0 * kd_194[k];

        t_147[k] = -3.0 * hd_105[k]
                   + f_0 * kd_195[k];

        t_148[k] = -3.0 * hd_106[k]
                   + f_0 * kd_196[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, hd_107, hd_108, hd_109, hd_110, \
                         hd_111, kd_197, kd_198, kd_199, kd_200, \
                         kd_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -3.0 * hd_107[k]
                   + f_0 * kd_197[k];

        t_150[k] = -4.0 * hd_108[k]
                   + f_0 * kd_198[k];

        t_151[k] = -4.0 * hd_109[k]
                   + f_0 * kd_199[k];

        t_152[k] = -4.0 * hd_110[k]
                   + f_0 * kd_200[k];

        t_153[k] = -4.0 * hd_111[k]
                   + f_0 * kd_201[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, hd_112, hd_113, hd_114, hd_115, \
                         hd_116, kd_202, kd_203, kd_204, kd_205, \
                         kd_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -4.0 * hd_112[k]
                   + f_0 * kd_202[k];

        t_155[k] = -4.0 * hd_113[k]
                   + f_0 * kd_203[k];

        t_156[k] = -5.0 * hd_114[k]
                   + f_0 * kd_204[k];

        t_157[k] = -5.0 * hd_115[k]
                   + f_0 * kd_205[k];

        t_158[k] = -5.0 * hd_116[k]
                   + f_0 * kd_206[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, hd_117, hd_118, hd_119, hd_120, \
                         hd_121, kd_207, kd_208, kd_209, kd_210, \
                         kd_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -5.0 * hd_117[k]
                   + f_0 * kd_207[k];

        t_160[k] = -5.0 * hd_118[k]
                   + f_0 * kd_208[k];

        t_161[k] = -5.0 * hd_119[k]
                   + f_0 * kd_209[k];

        t_162[k] = -6.0 * hd_120[k]
                   + f_0 * kd_210[k];

        t_163[k] = -6.0 * hd_121[k]
                   + f_0 * kd_211[k];
    }
}

static auto
compute_prim_geom_10_id_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hd, const size_t kd,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);

    const auto *hd_122 = buffer.data(hd + 122);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_124 = buffer.data(hd + 124);
    const auto *hd_125 = buffer.data(hd + 125);

    const auto *kd_212 = buffer.data(kd + 212);
    const auto *kd_213 = buffer.data(kd + 213);
    const auto *kd_214 = buffer.data(kd + 214);
    const auto *kd_215 = buffer.data(kd + 215);

#pragma omp simd aligned(t_164, t_165, t_166, t_167, hd_122, hd_123, hd_124, hd_125, kd_212, \
                         kd_213, kd_214, kd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -6.0 * hd_122[k]
                   + f_0 * kd_212[k];

        t_165[k] = -6.0 * hd_123[k]
                   + f_0 * kd_213[k];

        t_166[k] = -6.0 * hd_124[k]
                   + f_0 * kd_214[k];

        t_167[k] = -6.0 * hd_125[k]
                   + f_0 * kd_215[k];
    }
}

auto
compute_prim_geom_10_id_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t hd, const size_t kd,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_id_electron_repulsion_2_piece0(buffer, target, hd, kd, ncols, alpha);

    compute_prim_geom_10_id_electron_repulsion_2_piece1(buffer, target, hd, kd, ncols, alpha);
}

}  // namespace simdt2ceri
