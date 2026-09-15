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


#include "SimdElectronRepulsionGeom10VrrRecFH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_fh_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dh, const size_t gh,
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

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_41 = buffer.data(dh + 41);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_57 = buffer.data(dh + 57);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_61 = buffer.data(dh + 61);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_65 = buffer.data(dh + 65);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_67 = buffer.data(dh + 67);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_70 = buffer.data(dh + 70);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_74 = buffer.data(dh + 74);
    const auto *dh_75 = buffer.data(dh + 75);
    const auto *dh_76 = buffer.data(dh + 76);
    const auto *dh_77 = buffer.data(dh + 77);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_79 = buffer.data(dh + 79);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_84 = buffer.data(dh + 84);
    const auto *dh_85 = buffer.data(dh + 85);
    const auto *dh_86 = buffer.data(dh + 86);
    const auto *dh_87 = buffer.data(dh + 87);
    const auto *dh_88 = buffer.data(dh + 88);
    const auto *dh_89 = buffer.data(dh + 89);
    const auto *dh_90 = buffer.data(dh + 90);
    const auto *dh_91 = buffer.data(dh + 91);
    const auto *dh_92 = buffer.data(dh + 92);
    const auto *dh_93 = buffer.data(dh + 93);
    const auto *dh_94 = buffer.data(dh + 94);
    const auto *dh_95 = buffer.data(dh + 95);
    const auto *dh_96 = buffer.data(dh + 96);
    const auto *dh_97 = buffer.data(dh + 97);
    const auto *dh_98 = buffer.data(dh + 98);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_106 = buffer.data(dh + 106);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_109 = buffer.data(dh + 109);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_112 = buffer.data(dh + 112);
    const auto *dh_113 = buffer.data(dh + 113);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_115 = buffer.data(dh + 115);
    const auto *dh_116 = buffer.data(dh + 116);
    const auto *dh_117 = buffer.data(dh + 117);
    const auto *dh_118 = buffer.data(dh + 118);
    const auto *dh_119 = buffer.data(dh + 119);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_11 = buffer.data(gh + 11);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_13 = buffer.data(gh + 13);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_19 = buffer.data(gh + 19);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_43 = buffer.data(gh + 43);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_52 = buffer.data(gh + 52);
    const auto *gh_53 = buffer.data(gh + 53);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_61 = buffer.data(gh + 61);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_67 = buffer.data(gh + 67);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_71 = buffer.data(gh + 71);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_74 = buffer.data(gh + 74);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_76 = buffer.data(gh + 76);
    const auto *gh_77 = buffer.data(gh + 77);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_79 = buffer.data(gh + 79);
    const auto *gh_80 = buffer.data(gh + 80);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_82 = buffer.data(gh + 82);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_92 = buffer.data(gh + 92);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_106 = buffer.data(gh + 106);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_113 = buffer.data(gh + 113);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_115 = buffer.data(gh + 115);
    const auto *gh_116 = buffer.data(gh + 116);
    const auto *gh_117 = buffer.data(gh + 117);
    const auto *gh_118 = buffer.data(gh + 118);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_120 = buffer.data(gh + 120);
    const auto *gh_121 = buffer.data(gh + 121);
    const auto *gh_122 = buffer.data(gh + 122);
    const auto *gh_123 = buffer.data(gh + 123);
    const auto *gh_124 = buffer.data(gh + 124);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_128 = buffer.data(gh + 128);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_130 = buffer.data(gh + 130);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_133 = buffer.data(gh + 133);
    const auto *gh_134 = buffer.data(gh + 134);
    const auto *gh_135 = buffer.data(gh + 135);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_137 = buffer.data(gh + 137);
    const auto *gh_138 = buffer.data(gh + 138);
    const auto *gh_139 = buffer.data(gh + 139);
    const auto *gh_140 = buffer.data(gh + 140);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_142 = buffer.data(gh + 142);
    const auto *gh_143 = buffer.data(gh + 143);
    const auto *gh_144 = buffer.data(gh + 144);
    const auto *gh_145 = buffer.data(gh + 145);
    const auto *gh_146 = buffer.data(gh + 146);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_148 = buffer.data(gh + 148);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_151 = buffer.data(gh + 151);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_155 = buffer.data(gh + 155);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_157 = buffer.data(gh + 157);
    const auto *gh_158 = buffer.data(gh + 158);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_160 = buffer.data(gh + 160);
    const auto *gh_161 = buffer.data(gh + 161);
    const auto *gh_162 = buffer.data(gh + 162);
    const auto *gh_163 = buffer.data(gh + 163);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, dh_0, dh_1, dh_2, dh_3, dh_4, gh_0, gh_1, \
                         gh_2, gh_3, gh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -3.0 * dh_0[k]
                 + f_0 * gh_0[k];

        t_1[k] = -3.0 * dh_1[k]
                 + f_0 * gh_1[k];

        t_2[k] = -3.0 * dh_2[k]
                 + f_0 * gh_2[k];

        t_3[k] = -3.0 * dh_3[k]
                 + f_0 * gh_3[k];

        t_4[k] = -3.0 * dh_4[k]
                 + f_0 * gh_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, dh_5, dh_6, dh_7, dh_8, dh_9, gh_5, gh_6, \
                         gh_7, gh_8, gh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -3.0 * dh_5[k]
                 + f_0 * gh_5[k];

        t_6[k] = -3.0 * dh_6[k]
                 + f_0 * gh_6[k];

        t_7[k] = -3.0 * dh_7[k]
                 + f_0 * gh_7[k];

        t_8[k] = -3.0 * dh_8[k]
                 + f_0 * gh_8[k];

        t_9[k] = -3.0 * dh_9[k]
                 + f_0 * gh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, dh_10, dh_11, dh_12, dh_13, dh_14, \
                         gh_10, gh_11, gh_12, gh_13, gh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -3.0 * dh_10[k]
                  + f_0 * gh_10[k];

        t_11[k] = -3.0 * dh_11[k]
                  + f_0 * gh_11[k];

        t_12[k] = -3.0 * dh_12[k]
                  + f_0 * gh_12[k];

        t_13[k] = -3.0 * dh_13[k]
                  + f_0 * gh_13[k];

        t_14[k] = -3.0 * dh_14[k]
                  + f_0 * gh_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, dh_15, dh_16, dh_17, dh_18, dh_19, \
                         gh_15, gh_16, gh_17, gh_18, gh_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -3.0 * dh_15[k]
                  + f_0 * gh_15[k];

        t_16[k] = -3.0 * dh_16[k]
                  + f_0 * gh_16[k];

        t_17[k] = -3.0 * dh_17[k]
                  + f_0 * gh_17[k];

        t_18[k] = -3.0 * dh_18[k]
                  + f_0 * gh_18[k];

        t_19[k] = -3.0 * dh_19[k]
                  + f_0 * gh_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, dh_20, dh_21, dh_22, dh_23, dh_24, \
                         gh_20, gh_21, gh_22, gh_23, gh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -3.0 * dh_20[k]
                  + f_0 * gh_20[k];

        t_21[k] = -2.0 * dh_21[k]
                  + f_0 * gh_21[k];

        t_22[k] = -2.0 * dh_22[k]
                  + f_0 * gh_22[k];

        t_23[k] = -2.0 * dh_23[k]
                  + f_0 * gh_23[k];

        t_24[k] = -2.0 * dh_24[k]
                  + f_0 * gh_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, dh_25, dh_26, dh_27, dh_28, dh_29, \
                         gh_25, gh_26, gh_27, gh_28, gh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -2.0 * dh_25[k]
                  + f_0 * gh_25[k];

        t_26[k] = -2.0 * dh_26[k]
                  + f_0 * gh_26[k];

        t_27[k] = -2.0 * dh_27[k]
                  + f_0 * gh_27[k];

        t_28[k] = -2.0 * dh_28[k]
                  + f_0 * gh_28[k];

        t_29[k] = -2.0 * dh_29[k]
                  + f_0 * gh_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, dh_30, dh_31, dh_32, dh_33, dh_34, \
                         gh_30, gh_31, gh_32, gh_33, gh_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -2.0 * dh_30[k]
                  + f_0 * gh_30[k];

        t_31[k] = -2.0 * dh_31[k]
                  + f_0 * gh_31[k];

        t_32[k] = -2.0 * dh_32[k]
                  + f_0 * gh_32[k];

        t_33[k] = -2.0 * dh_33[k]
                  + f_0 * gh_33[k];

        t_34[k] = -2.0 * dh_34[k]
                  + f_0 * gh_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, dh_35, dh_36, dh_37, dh_38, dh_39, \
                         gh_35, gh_36, gh_37, gh_38, gh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -2.0 * dh_35[k]
                  + f_0 * gh_35[k];

        t_36[k] = -2.0 * dh_36[k]
                  + f_0 * gh_36[k];

        t_37[k] = -2.0 * dh_37[k]
                  + f_0 * gh_37[k];

        t_38[k] = -2.0 * dh_38[k]
                  + f_0 * gh_38[k];

        t_39[k] = -2.0 * dh_39[k]
                  + f_0 * gh_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, dh_40, dh_41, dh_42, dh_43, dh_44, \
                         gh_40, gh_41, gh_42, gh_43, gh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -2.0 * dh_40[k]
                  + f_0 * gh_40[k];

        t_41[k] = -2.0 * dh_41[k]
                  + f_0 * gh_41[k];

        t_42[k] = -2.0 * dh_42[k]
                  + f_0 * gh_42[k];

        t_43[k] = -2.0 * dh_43[k]
                  + f_0 * gh_43[k];

        t_44[k] = -2.0 * dh_44[k]
                  + f_0 * gh_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, dh_45, dh_46, dh_47, dh_48, dh_49, \
                         gh_45, gh_46, gh_47, gh_48, gh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * dh_45[k]
                  + f_0 * gh_45[k];

        t_46[k] = -2.0 * dh_46[k]
                  + f_0 * gh_46[k];

        t_47[k] = -2.0 * dh_47[k]
                  + f_0 * gh_47[k];

        t_48[k] = -2.0 * dh_48[k]
                  + f_0 * gh_48[k];

        t_49[k] = -2.0 * dh_49[k]
                  + f_0 * gh_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, dh_50, dh_51, dh_52, dh_53, dh_54, \
                         gh_50, gh_51, gh_52, gh_53, gh_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * dh_50[k]
                  + f_0 * gh_50[k];

        t_51[k] = -2.0 * dh_51[k]
                  + f_0 * gh_51[k];

        t_52[k] = -2.0 * dh_52[k]
                  + f_0 * gh_52[k];

        t_53[k] = -2.0 * dh_53[k]
                  + f_0 * gh_53[k];

        t_54[k] = -2.0 * dh_54[k]
                  + f_0 * gh_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, dh_55, dh_56, dh_57, dh_58, dh_59, \
                         gh_55, gh_56, gh_57, gh_58, gh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * dh_55[k]
                  + f_0 * gh_55[k];

        t_56[k] = -2.0 * dh_56[k]
                  + f_0 * gh_56[k];

        t_57[k] = -2.0 * dh_57[k]
                  + f_0 * gh_57[k];

        t_58[k] = -2.0 * dh_58[k]
                  + f_0 * gh_58[k];

        t_59[k] = -2.0 * dh_59[k]
                  + f_0 * gh_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, dh_60, dh_61, dh_62, dh_63, dh_64, \
                         gh_60, gh_61, gh_62, gh_63, gh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -2.0 * dh_60[k]
                  + f_0 * gh_60[k];

        t_61[k] = -2.0 * dh_61[k]
                  + f_0 * gh_61[k];

        t_62[k] = -2.0 * dh_62[k]
                  + f_0 * gh_62[k];

        t_63[k] = -dh_63[k]
                  + f_0 * gh_63[k];

        t_64[k] = -dh_64[k]
                  + f_0 * gh_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, dh_65, dh_66, dh_67, dh_68, dh_69, \
                         gh_65, gh_66, gh_67, gh_68, gh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -dh_65[k]
                  + f_0 * gh_65[k];

        t_66[k] = -dh_66[k]
                  + f_0 * gh_66[k];

        t_67[k] = -dh_67[k]
                  + f_0 * gh_67[k];

        t_68[k] = -dh_68[k]
                  + f_0 * gh_68[k];

        t_69[k] = -dh_69[k]
                  + f_0 * gh_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, dh_70, dh_71, dh_72, dh_73, dh_74, \
                         gh_70, gh_71, gh_72, gh_73, gh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -dh_70[k]
                  + f_0 * gh_70[k];

        t_71[k] = -dh_71[k]
                  + f_0 * gh_71[k];

        t_72[k] = -dh_72[k]
                  + f_0 * gh_72[k];

        t_73[k] = -dh_73[k]
                  + f_0 * gh_73[k];

        t_74[k] = -dh_74[k]
                  + f_0 * gh_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, dh_75, dh_76, dh_77, dh_78, dh_79, \
                         gh_75, gh_76, gh_77, gh_78, gh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -dh_75[k]
                  + f_0 * gh_75[k];

        t_76[k] = -dh_76[k]
                  + f_0 * gh_76[k];

        t_77[k] = -dh_77[k]
                  + f_0 * gh_77[k];

        t_78[k] = -dh_78[k]
                  + f_0 * gh_78[k];

        t_79[k] = -dh_79[k]
                  + f_0 * gh_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, dh_80, dh_81, dh_82, dh_83, dh_84, \
                         gh_80, gh_81, gh_82, gh_83, gh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -dh_80[k]
                  + f_0 * gh_80[k];

        t_81[k] = -dh_81[k]
                  + f_0 * gh_81[k];

        t_82[k] = -dh_82[k]
                  + f_0 * gh_82[k];

        t_83[k] = -dh_83[k]
                  + f_0 * gh_83[k];

        t_84[k] = -dh_84[k]
                  + f_0 * gh_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, dh_85, dh_86, dh_87, dh_88, dh_89, \
                         gh_85, gh_86, gh_87, gh_88, gh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -dh_85[k]
                  + f_0 * gh_85[k];

        t_86[k] = -dh_86[k]
                  + f_0 * gh_86[k];

        t_87[k] = -dh_87[k]
                  + f_0 * gh_87[k];

        t_88[k] = -dh_88[k]
                  + f_0 * gh_88[k];

        t_89[k] = -dh_89[k]
                  + f_0 * gh_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, dh_90, dh_91, dh_92, dh_93, dh_94, \
                         gh_90, gh_91, gh_92, gh_93, gh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -dh_90[k]
                  + f_0 * gh_90[k];

        t_91[k] = -dh_91[k]
                  + f_0 * gh_91[k];

        t_92[k] = -dh_92[k]
                  + f_0 * gh_92[k];

        t_93[k] = -dh_93[k]
                  + f_0 * gh_93[k];

        t_94[k] = -dh_94[k]
                  + f_0 * gh_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, dh_95, dh_96, dh_97, dh_98, dh_99, \
                         gh_95, gh_96, gh_97, gh_98, gh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -dh_95[k]
                  + f_0 * gh_95[k];

        t_96[k] = -dh_96[k]
                  + f_0 * gh_96[k];

        t_97[k] = -dh_97[k]
                  + f_0 * gh_97[k];

        t_98[k] = -dh_98[k]
                  + f_0 * gh_98[k];

        t_99[k] = -dh_99[k]
                  + f_0 * gh_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, dh_100, dh_101, dh_102, dh_103, \
                         dh_104, gh_100, gh_101, gh_102, gh_103, \
                         gh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -dh_100[k]
                   + f_0 * gh_100[k];

        t_101[k] = -dh_101[k]
                   + f_0 * gh_101[k];

        t_102[k] = -dh_102[k]
                   + f_0 * gh_102[k];

        t_103[k] = -dh_103[k]
                   + f_0 * gh_103[k];

        t_104[k] = -dh_104[k]
                   + f_0 * gh_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, dh_105, dh_106, dh_107, dh_108, \
                         dh_109, gh_105, gh_106, gh_107, gh_108, \
                         gh_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -dh_105[k]
                   + f_0 * gh_105[k];

        t_106[k] = -dh_106[k]
                   + f_0 * gh_106[k];

        t_107[k] = -dh_107[k]
                   + f_0 * gh_107[k];

        t_108[k] = -dh_108[k]
                   + f_0 * gh_108[k];

        t_109[k] = -dh_109[k]
                   + f_0 * gh_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, dh_110, dh_111, dh_112, dh_113, \
                         dh_114, gh_110, gh_111, gh_112, gh_113, \
                         gh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -dh_110[k]
                   + f_0 * gh_110[k];

        t_111[k] = -dh_111[k]
                   + f_0 * gh_111[k];

        t_112[k] = -dh_112[k]
                   + f_0 * gh_112[k];

        t_113[k] = -dh_113[k]
                   + f_0 * gh_113[k];

        t_114[k] = -dh_114[k]
                   + f_0 * gh_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, dh_115, dh_116, dh_117, dh_118, \
                         dh_119, gh_115, gh_116, gh_117, gh_118, \
                         gh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -dh_115[k]
                   + f_0 * gh_115[k];

        t_116[k] = -dh_116[k]
                   + f_0 * gh_116[k];

        t_117[k] = -dh_117[k]
                   + f_0 * gh_117[k];

        t_118[k] = -dh_118[k]
                   + f_0 * gh_118[k];

        t_119[k] = -dh_119[k]
                   + f_0 * gh_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, dh_120, dh_121, dh_122, dh_123, \
                         dh_124, gh_120, gh_121, gh_122, gh_123, \
                         gh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -dh_120[k]
                   + f_0 * gh_120[k];

        t_121[k] = -dh_121[k]
                   + f_0 * gh_121[k];

        t_122[k] = -dh_122[k]
                   + f_0 * gh_122[k];

        t_123[k] = -dh_123[k]
                   + f_0 * gh_123[k];

        t_124[k] = -dh_124[k]
                   + f_0 * gh_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, t_131, dh_125, gh_125, \
                         gh_126, gh_127, gh_128, gh_129, gh_130, \
                         gh_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -dh_125[k]
                   + f_0 * gh_125[k];

        t_126[k] = f_0 * gh_126[k];

        t_127[k] = f_0 * gh_127[k];

        t_128[k] = f_0 * gh_128[k];

        t_129[k] = f_0 * gh_129[k];

        t_130[k] = f_0 * gh_130[k];

        t_131[k] = f_0 * gh_131[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, t_137, t_138, t_139, gh_132, \
                         gh_133, gh_134, gh_135, gh_136, gh_137, gh_138, \
                         gh_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_0 * gh_132[k];

        t_133[k] = f_0 * gh_133[k];

        t_134[k] = f_0 * gh_134[k];

        t_135[k] = f_0 * gh_135[k];

        t_136[k] = f_0 * gh_136[k];

        t_137[k] = f_0 * gh_137[k];

        t_138[k] = f_0 * gh_138[k];

        t_139[k] = f_0 * gh_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, t_146, t_147, gh_140, \
                         gh_141, gh_142, gh_143, gh_144, gh_145, gh_146, \
                         gh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * gh_140[k];

        t_141[k] = f_0 * gh_141[k];

        t_142[k] = f_0 * gh_142[k];

        t_143[k] = f_0 * gh_143[k];

        t_144[k] = f_0 * gh_144[k];

        t_145[k] = f_0 * gh_145[k];

        t_146[k] = f_0 * gh_146[k];

        t_147[k] = f_0 * gh_147[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, t_153, t_154, t_155, gh_148, \
                         gh_149, gh_150, gh_151, gh_152, gh_153, gh_154, \
                         gh_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * gh_148[k];

        t_149[k] = f_0 * gh_149[k];

        t_150[k] = f_0 * gh_150[k];

        t_151[k] = f_0 * gh_151[k];

        t_152[k] = f_0 * gh_152[k];

        t_153[k] = f_0 * gh_153[k];

        t_154[k] = f_0 * gh_154[k];

        t_155[k] = f_0 * gh_155[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, t_162, t_163, gh_156, \
                         gh_157, gh_158, gh_159, gh_160, gh_161, gh_162, \
                         gh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * gh_156[k];

        t_157[k] = f_0 * gh_157[k];

        t_158[k] = f_0 * gh_158[k];

        t_159[k] = f_0 * gh_159[k];

        t_160[k] = f_0 * gh_160[k];

        t_161[k] = f_0 * gh_161[k];

        t_162[k] = f_0 * gh_162[k];

        t_163[k] = f_0 * gh_163[k];
    }
}

static auto
compute_prim_geom_10_fh_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gh, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_166 = buffer.data(gh + 166);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_169 = buffer.data(gh + 169);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_172 = buffer.data(gh + 172);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_175 = buffer.data(gh + 175);
    const auto *gh_176 = buffer.data(gh + 176);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_178 = buffer.data(gh + 178);
    const auto *gh_179 = buffer.data(gh + 179);
    const auto *gh_180 = buffer.data(gh + 180);
    const auto *gh_181 = buffer.data(gh + 181);
    const auto *gh_182 = buffer.data(gh + 182);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_188 = buffer.data(gh + 188);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_190 = buffer.data(gh + 190);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_192 = buffer.data(gh + 192);
    const auto *gh_193 = buffer.data(gh + 193);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_195 = buffer.data(gh + 195);
    const auto *gh_196 = buffer.data(gh + 196);
    const auto *gh_197 = buffer.data(gh + 197);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_199 = buffer.data(gh + 199);
    const auto *gh_200 = buffer.data(gh + 200);
    const auto *gh_201 = buffer.data(gh + 201);
    const auto *gh_202 = buffer.data(gh + 202);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_204 = buffer.data(gh + 204);
    const auto *gh_205 = buffer.data(gh + 205);
    const auto *gh_206 = buffer.data(gh + 206);
    const auto *gh_207 = buffer.data(gh + 207);
    const auto *gh_208 = buffer.data(gh + 208);
    const auto *gh_209 = buffer.data(gh + 209);

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, t_169, t_170, t_171, gh_164, \
                         gh_165, gh_166, gh_167, gh_168, gh_169, gh_170, \
                         gh_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_0 * gh_164[k];

        t_165[k] = f_0 * gh_165[k];

        t_166[k] = f_0 * gh_166[k];

        t_167[k] = f_0 * gh_167[k];

        t_168[k] = f_0 * gh_168[k];

        t_169[k] = f_0 * gh_169[k];

        t_170[k] = f_0 * gh_170[k];

        t_171[k] = f_0 * gh_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, t_177, t_178, t_179, gh_172, \
                         gh_173, gh_174, gh_175, gh_176, gh_177, gh_178, \
                         gh_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_0 * gh_172[k];

        t_173[k] = f_0 * gh_173[k];

        t_174[k] = f_0 * gh_174[k];

        t_175[k] = f_0 * gh_175[k];

        t_176[k] = f_0 * gh_176[k];

        t_177[k] = f_0 * gh_177[k];

        t_178[k] = f_0 * gh_178[k];

        t_179[k] = f_0 * gh_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, t_185, t_186, t_187, gh_180, \
                         gh_181, gh_182, gh_183, gh_184, gh_185, gh_186, \
                         gh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_0 * gh_180[k];

        t_181[k] = f_0 * gh_181[k];

        t_182[k] = f_0 * gh_182[k];

        t_183[k] = f_0 * gh_183[k];

        t_184[k] = f_0 * gh_184[k];

        t_185[k] = f_0 * gh_185[k];

        t_186[k] = f_0 * gh_186[k];

        t_187[k] = f_0 * gh_187[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, t_193, t_194, t_195, gh_188, \
                         gh_189, gh_190, gh_191, gh_192, gh_193, gh_194, \
                         gh_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_0 * gh_188[k];

        t_189[k] = f_0 * gh_189[k];

        t_190[k] = f_0 * gh_190[k];

        t_191[k] = f_0 * gh_191[k];

        t_192[k] = f_0 * gh_192[k];

        t_193[k] = f_0 * gh_193[k];

        t_194[k] = f_0 * gh_194[k];

        t_195[k] = f_0 * gh_195[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, t_201, t_202, t_203, gh_196, \
                         gh_197, gh_198, gh_199, gh_200, gh_201, gh_202, \
                         gh_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_0 * gh_196[k];

        t_197[k] = f_0 * gh_197[k];

        t_198[k] = f_0 * gh_198[k];

        t_199[k] = f_0 * gh_199[k];

        t_200[k] = f_0 * gh_200[k];

        t_201[k] = f_0 * gh_201[k];

        t_202[k] = f_0 * gh_202[k];

        t_203[k] = f_0 * gh_203[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, t_209, gh_204, gh_205, gh_206, \
                         gh_207, gh_208, gh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_0 * gh_204[k];

        t_205[k] = f_0 * gh_205[k];

        t_206[k] = f_0 * gh_206[k];

        t_207[k] = f_0 * gh_207[k];

        t_208[k] = f_0 * gh_208[k];

        t_209[k] = f_0 * gh_209[k];
    }
}

auto
compute_prim_geom_10_fh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t dh, const size_t gh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fh_electron_repulsion_0_piece0(buffer, target, dh, gh, ncols, alpha);

    compute_prim_geom_10_fh_electron_repulsion_0_piece1(buffer, target, gh, ncols, alpha);
}

static auto
compute_prim_geom_10_fh_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dh, const size_t gh,
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

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_41 = buffer.data(dh + 41);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_57 = buffer.data(dh + 57);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_61 = buffer.data(dh + 61);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_65 = buffer.data(dh + 65);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_67 = buffer.data(dh + 67);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_70 = buffer.data(dh + 70);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_74 = buffer.data(dh + 74);
    const auto *dh_75 = buffer.data(dh + 75);
    const auto *dh_76 = buffer.data(dh + 76);
    const auto *dh_77 = buffer.data(dh + 77);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_79 = buffer.data(dh + 79);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_84 = buffer.data(dh + 84);
    const auto *dh_85 = buffer.data(dh + 85);
    const auto *dh_86 = buffer.data(dh + 86);
    const auto *dh_87 = buffer.data(dh + 87);
    const auto *dh_88 = buffer.data(dh + 88);
    const auto *dh_89 = buffer.data(dh + 89);
    const auto *dh_90 = buffer.data(dh + 90);
    const auto *dh_91 = buffer.data(dh + 91);
    const auto *dh_92 = buffer.data(dh + 92);
    const auto *dh_93 = buffer.data(dh + 93);
    const auto *dh_94 = buffer.data(dh + 94);
    const auto *dh_95 = buffer.data(dh + 95);
    const auto *dh_96 = buffer.data(dh + 96);
    const auto *dh_97 = buffer.data(dh + 97);
    const auto *dh_98 = buffer.data(dh + 98);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_105 = buffer.data(dh + 105);

    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_67 = buffer.data(gh + 67);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_71 = buffer.data(gh + 71);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_74 = buffer.data(gh + 74);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_76 = buffer.data(gh + 76);
    const auto *gh_77 = buffer.data(gh + 77);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_79 = buffer.data(gh + 79);
    const auto *gh_80 = buffer.data(gh + 80);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_82 = buffer.data(gh + 82);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_92 = buffer.data(gh + 92);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_128 = buffer.data(gh + 128);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_130 = buffer.data(gh + 130);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_133 = buffer.data(gh + 133);
    const auto *gh_134 = buffer.data(gh + 134);
    const auto *gh_135 = buffer.data(gh + 135);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_137 = buffer.data(gh + 137);
    const auto *gh_138 = buffer.data(gh + 138);
    const auto *gh_139 = buffer.data(gh + 139);
    const auto *gh_140 = buffer.data(gh + 140);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_142 = buffer.data(gh + 142);
    const auto *gh_143 = buffer.data(gh + 143);
    const auto *gh_144 = buffer.data(gh + 144);
    const auto *gh_145 = buffer.data(gh + 145);
    const auto *gh_146 = buffer.data(gh + 146);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_148 = buffer.data(gh + 148);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_151 = buffer.data(gh + 151);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_155 = buffer.data(gh + 155);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_157 = buffer.data(gh + 157);
    const auto *gh_158 = buffer.data(gh + 158);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_160 = buffer.data(gh + 160);
    const auto *gh_161 = buffer.data(gh + 161);
    const auto *gh_162 = buffer.data(gh + 162);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_166 = buffer.data(gh + 166);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_169 = buffer.data(gh + 169);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_172 = buffer.data(gh + 172);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_175 = buffer.data(gh + 175);
    const auto *gh_176 = buffer.data(gh + 176);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_178 = buffer.data(gh + 178);
    const auto *gh_179 = buffer.data(gh + 179);
    const auto *gh_180 = buffer.data(gh + 180);
    const auto *gh_181 = buffer.data(gh + 181);
    const auto *gh_182 = buffer.data(gh + 182);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_188 = buffer.data(gh + 188);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_211 = buffer.data(gh + 211);
    const auto *gh_212 = buffer.data(gh + 212);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_214 = buffer.data(gh + 214);
    const auto *gh_215 = buffer.data(gh + 215);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_217 = buffer.data(gh + 217);
    const auto *gh_218 = buffer.data(gh + 218);
    const auto *gh_219 = buffer.data(gh + 219);
    const auto *gh_220 = buffer.data(gh + 220);
    const auto *gh_221 = buffer.data(gh + 221);
    const auto *gh_222 = buffer.data(gh + 222);
    const auto *gh_223 = buffer.data(gh + 223);
    const auto *gh_224 = buffer.data(gh + 224);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_226 = buffer.data(gh + 226);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_231 = buffer.data(gh + 231);
    const auto *gh_232 = buffer.data(gh + 232);
    const auto *gh_233 = buffer.data(gh + 233);
    const auto *gh_234 = buffer.data(gh + 234);
    const auto *gh_235 = buffer.data(gh + 235);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_237 = buffer.data(gh + 237);
    const auto *gh_238 = buffer.data(gh + 238);
    const auto *gh_239 = buffer.data(gh + 239);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_241 = buffer.data(gh + 241);
    const auto *gh_242 = buffer.data(gh + 242);
    const auto *gh_243 = buffer.data(gh + 243);
    const auto *gh_244 = buffer.data(gh + 244);
    const auto *gh_245 = buffer.data(gh + 245);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gh_21, gh_22, gh_23, gh_24, \
                         gh_25, gh_26, gh_27, gh_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gh_21[k];

        t_1[k] = f_0 * gh_22[k];

        t_2[k] = f_0 * gh_23[k];

        t_3[k] = f_0 * gh_24[k];

        t_4[k] = f_0 * gh_25[k];

        t_5[k] = f_0 * gh_26[k];

        t_6[k] = f_0 * gh_27[k];

        t_7[k] = f_0 * gh_28[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, gh_29, gh_30, gh_31, \
                         gh_32, gh_33, gh_34, gh_35, gh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gh_29[k];

        t_9[k] = f_0 * gh_30[k];

        t_10[k] = f_0 * gh_31[k];

        t_11[k] = f_0 * gh_32[k];

        t_12[k] = f_0 * gh_33[k];

        t_13[k] = f_0 * gh_34[k];

        t_14[k] = f_0 * gh_35[k];

        t_15[k] = f_0 * gh_36[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, dh_0, dh_1, gh_37, gh_38, \
                         gh_39, gh_40, gh_41, gh_63, gh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * gh_37[k];

        t_17[k] = f_0 * gh_38[k];

        t_18[k] = f_0 * gh_39[k];

        t_19[k] = f_0 * gh_40[k];

        t_20[k] = f_0 * gh_41[k];

        t_21[k] = -dh_0[k]
                  + f_0 * gh_63[k];

        t_22[k] = -dh_1[k]
                  + f_0 * gh_64[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, dh_2, dh_3, dh_4, dh_5, dh_6, gh_65, \
                         gh_66, gh_67, gh_68, gh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -dh_2[k]
                  + f_0 * gh_65[k];

        t_24[k] = -dh_3[k]
                  + f_0 * gh_66[k];

        t_25[k] = -dh_4[k]
                  + f_0 * gh_67[k];

        t_26[k] = -dh_5[k]
                  + f_0 * gh_68[k];

        t_27[k] = -dh_6[k]
                  + f_0 * gh_69[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, dh_7, dh_8, dh_9, dh_10, dh_11, gh_70, \
                         gh_71, gh_72, gh_73, gh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -dh_7[k]
                  + f_0 * gh_70[k];

        t_29[k] = -dh_8[k]
                  + f_0 * gh_71[k];

        t_30[k] = -dh_9[k]
                  + f_0 * gh_72[k];

        t_31[k] = -dh_10[k]
                  + f_0 * gh_73[k];

        t_32[k] = -dh_11[k]
                  + f_0 * gh_74[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, dh_12, dh_13, dh_14, dh_15, dh_16, \
                         gh_75, gh_76, gh_77, gh_78, gh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -dh_12[k]
                  + f_0 * gh_75[k];

        t_34[k] = -dh_13[k]
                  + f_0 * gh_76[k];

        t_35[k] = -dh_14[k]
                  + f_0 * gh_77[k];

        t_36[k] = -dh_15[k]
                  + f_0 * gh_78[k];

        t_37[k] = -dh_16[k]
                  + f_0 * gh_79[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, dh_17, dh_18, dh_19, dh_20, \
                         gh_80, gh_81, gh_82, gh_83, gh_84, gh_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -dh_17[k]
                  + f_0 * gh_80[k];

        t_39[k] = -dh_18[k]
                  + f_0 * gh_81[k];

        t_40[k] = -dh_19[k]
                  + f_0 * gh_82[k];

        t_41[k] = -dh_20[k]
                  + f_0 * gh_83[k];

        t_42[k] = f_0 * gh_84[k];

        t_43[k] = f_0 * gh_85[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, t_51, gh_86, gh_87, gh_88, \
                         gh_89, gh_90, gh_91, gh_92, gh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * gh_86[k];

        t_45[k] = f_0 * gh_87[k];

        t_46[k] = f_0 * gh_88[k];

        t_47[k] = f_0 * gh_89[k];

        t_48[k] = f_0 * gh_90[k];

        t_49[k] = f_0 * gh_91[k];

        t_50[k] = f_0 * gh_92[k];

        t_51[k] = f_0 * gh_93[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, t_58, t_59, gh_94, gh_95, gh_96, \
                         gh_97, gh_98, gh_99, gh_100, gh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * gh_94[k];

        t_53[k] = f_0 * gh_95[k];

        t_54[k] = f_0 * gh_96[k];

        t_55[k] = f_0 * gh_97[k];

        t_56[k] = f_0 * gh_98[k];

        t_57[k] = f_0 * gh_99[k];

        t_58[k] = f_0 * gh_100[k];

        t_59[k] = f_0 * gh_101[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, dh_21, dh_22, dh_23, gh_102, \
                         gh_103, gh_104, gh_126, gh_127, gh_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * gh_102[k];

        t_61[k] = f_0 * gh_103[k];

        t_62[k] = f_0 * gh_104[k];

        t_63[k] = -2.0 * dh_21[k]
                  + f_0 * gh_126[k];

        t_64[k] = -2.0 * dh_22[k]
                  + f_0 * gh_127[k];

        t_65[k] = -2.0 * dh_23[k]
                  + f_0 * gh_128[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, dh_24, dh_25, dh_26, dh_27, dh_28, \
                         gh_129, gh_130, gh_131, gh_132, gh_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -2.0 * dh_24[k]
                  + f_0 * gh_129[k];

        t_67[k] = -2.0 * dh_25[k]
                  + f_0 * gh_130[k];

        t_68[k] = -2.0 * dh_26[k]
                  + f_0 * gh_131[k];

        t_69[k] = -2.0 * dh_27[k]
                  + f_0 * gh_132[k];

        t_70[k] = -2.0 * dh_28[k]
                  + f_0 * gh_133[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, dh_29, dh_30, dh_31, dh_32, dh_33, \
                         gh_134, gh_135, gh_136, gh_137, gh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -2.0 * dh_29[k]
                  + f_0 * gh_134[k];

        t_72[k] = -2.0 * dh_30[k]
                  + f_0 * gh_135[k];

        t_73[k] = -2.0 * dh_31[k]
                  + f_0 * gh_136[k];

        t_74[k] = -2.0 * dh_32[k]
                  + f_0 * gh_137[k];

        t_75[k] = -2.0 * dh_33[k]
                  + f_0 * gh_138[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, dh_34, dh_35, dh_36, dh_37, dh_38, \
                         gh_139, gh_140, gh_141, gh_142, gh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -2.0 * dh_34[k]
                  + f_0 * gh_139[k];

        t_77[k] = -2.0 * dh_35[k]
                  + f_0 * gh_140[k];

        t_78[k] = -2.0 * dh_36[k]
                  + f_0 * gh_141[k];

        t_79[k] = -2.0 * dh_37[k]
                  + f_0 * gh_142[k];

        t_80[k] = -2.0 * dh_38[k]
                  + f_0 * gh_143[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, dh_39, dh_40, dh_41, dh_42, dh_43, \
                         gh_144, gh_145, gh_146, gh_147, gh_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -2.0 * dh_39[k]
                  + f_0 * gh_144[k];

        t_82[k] = -2.0 * dh_40[k]
                  + f_0 * gh_145[k];

        t_83[k] = -2.0 * dh_41[k]
                  + f_0 * gh_146[k];

        t_84[k] = -dh_42[k]
                  + f_0 * gh_147[k];

        t_85[k] = -dh_43[k]
                  + f_0 * gh_148[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, dh_44, dh_45, dh_46, dh_47, dh_48, \
                         gh_149, gh_150, gh_151, gh_152, gh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -dh_44[k]
                  + f_0 * gh_149[k];

        t_87[k] = -dh_45[k]
                  + f_0 * gh_150[k];

        t_88[k] = -dh_46[k]
                  + f_0 * gh_151[k];

        t_89[k] = -dh_47[k]
                  + f_0 * gh_152[k];

        t_90[k] = -dh_48[k]
                  + f_0 * gh_153[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, dh_49, dh_50, dh_51, dh_52, dh_53, \
                         gh_154, gh_155, gh_156, gh_157, gh_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = -dh_49[k]
                  + f_0 * gh_154[k];

        t_92[k] = -dh_50[k]
                  + f_0 * gh_155[k];

        t_93[k] = -dh_51[k]
                  + f_0 * gh_156[k];

        t_94[k] = -dh_52[k]
                  + f_0 * gh_157[k];

        t_95[k] = -dh_53[k]
                  + f_0 * gh_158[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, dh_54, dh_55, dh_56, dh_57, dh_58, \
                         gh_159, gh_160, gh_161, gh_162, gh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = -dh_54[k]
                  + f_0 * gh_159[k];

        t_97[k] = -dh_55[k]
                  + f_0 * gh_160[k];

        t_98[k] = -dh_56[k]
                  + f_0 * gh_161[k];

        t_99[k] = -dh_57[k]
                  + f_0 * gh_162[k];

        t_100[k] = -dh_58[k]
                   + f_0 * gh_163[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, dh_59, dh_60, dh_61, dh_62, \
                         gh_164, gh_165, gh_166, gh_167, gh_168, \
                         gh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = -dh_59[k]
                   + f_0 * gh_164[k];

        t_102[k] = -dh_60[k]
                   + f_0 * gh_165[k];

        t_103[k] = -dh_61[k]
                   + f_0 * gh_166[k];

        t_104[k] = -dh_62[k]
                   + f_0 * gh_167[k];

        t_105[k] = f_0 * gh_168[k];

        t_106[k] = f_0 * gh_169[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, t_114, gh_170, \
                         gh_171, gh_172, gh_173, gh_174, gh_175, gh_176, \
                         gh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * gh_170[k];

        t_108[k] = f_0 * gh_171[k];

        t_109[k] = f_0 * gh_172[k];

        t_110[k] = f_0 * gh_173[k];

        t_111[k] = f_0 * gh_174[k];

        t_112[k] = f_0 * gh_175[k];

        t_113[k] = f_0 * gh_176[k];

        t_114[k] = f_0 * gh_177[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, t_121, t_122, gh_178, \
                         gh_179, gh_180, gh_181, gh_182, gh_183, gh_184, \
                         gh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_0 * gh_178[k];

        t_116[k] = f_0 * gh_179[k];

        t_117[k] = f_0 * gh_180[k];

        t_118[k] = f_0 * gh_181[k];

        t_119[k] = f_0 * gh_182[k];

        t_120[k] = f_0 * gh_183[k];

        t_121[k] = f_0 * gh_184[k];

        t_122[k] = f_0 * gh_185[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, t_127, t_128, dh_63, dh_64, dh_65, \
                         gh_186, gh_187, gh_188, gh_210, gh_211, \
                         gh_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_0 * gh_186[k];

        t_124[k] = f_0 * gh_187[k];

        t_125[k] = f_0 * gh_188[k];

        t_126[k] = -3.0 * dh_63[k]
                   + f_0 * gh_210[k];

        t_127[k] = -3.0 * dh_64[k]
                   + f_0 * gh_211[k];

        t_128[k] = -3.0 * dh_65[k]
                   + f_0 * gh_212[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, dh_66, dh_67, dh_68, dh_69, dh_70, \
                         gh_213, gh_214, gh_215, gh_216, gh_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -3.0 * dh_66[k]
                   + f_0 * gh_213[k];

        t_130[k] = -3.0 * dh_67[k]
                   + f_0 * gh_214[k];

        t_131[k] = -3.0 * dh_68[k]
                   + f_0 * gh_215[k];

        t_132[k] = -3.0 * dh_69[k]
                   + f_0 * gh_216[k];

        t_133[k] = -3.0 * dh_70[k]
                   + f_0 * gh_217[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, dh_71, dh_72, dh_73, dh_74, dh_75, \
                         gh_218, gh_219, gh_220, gh_221, gh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -3.0 * dh_71[k]
                   + f_0 * gh_218[k];

        t_135[k] = -3.0 * dh_72[k]
                   + f_0 * gh_219[k];

        t_136[k] = -3.0 * dh_73[k]
                   + f_0 * gh_220[k];

        t_137[k] = -3.0 * dh_74[k]
                   + f_0 * gh_221[k];

        t_138[k] = -3.0 * dh_75[k]
                   + f_0 * gh_222[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, dh_76, dh_77, dh_78, dh_79, dh_80, \
                         gh_223, gh_224, gh_225, gh_226, gh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -3.0 * dh_76[k]
                   + f_0 * gh_223[k];

        t_140[k] = -3.0 * dh_77[k]
                   + f_0 * gh_224[k];

        t_141[k] = -3.0 * dh_78[k]
                   + f_0 * gh_225[k];

        t_142[k] = -3.0 * dh_79[k]
                   + f_0 * gh_226[k];

        t_143[k] = -3.0 * dh_80[k]
                   + f_0 * gh_227[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, dh_81, dh_82, dh_83, dh_84, dh_85, \
                         gh_228, gh_229, gh_230, gh_231, gh_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -3.0 * dh_81[k]
                   + f_0 * gh_228[k];

        t_145[k] = -3.0 * dh_82[k]
                   + f_0 * gh_229[k];

        t_146[k] = -3.0 * dh_83[k]
                   + f_0 * gh_230[k];

        t_147[k] = -2.0 * dh_84[k]
                   + f_0 * gh_231[k];

        t_148[k] = -2.0 * dh_85[k]
                   + f_0 * gh_232[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, dh_86, dh_87, dh_88, dh_89, dh_90, \
                         gh_233, gh_234, gh_235, gh_236, gh_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -2.0 * dh_86[k]
                   + f_0 * gh_233[k];

        t_150[k] = -2.0 * dh_87[k]
                   + f_0 * gh_234[k];

        t_151[k] = -2.0 * dh_88[k]
                   + f_0 * gh_235[k];

        t_152[k] = -2.0 * dh_89[k]
                   + f_0 * gh_236[k];

        t_153[k] = -2.0 * dh_90[k]
                   + f_0 * gh_237[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, dh_91, dh_92, dh_93, dh_94, dh_95, \
                         gh_238, gh_239, gh_240, gh_241, gh_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -2.0 * dh_91[k]
                   + f_0 * gh_238[k];

        t_155[k] = -2.0 * dh_92[k]
                   + f_0 * gh_239[k];

        t_156[k] = -2.0 * dh_93[k]
                   + f_0 * gh_240[k];

        t_157[k] = -2.0 * dh_94[k]
                   + f_0 * gh_241[k];

        t_158[k] = -2.0 * dh_95[k]
                   + f_0 * gh_242[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, dh_96, dh_97, dh_98, dh_99, \
                         dh_100, gh_243, gh_244, gh_245, gh_246, \
                         gh_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -2.0 * dh_96[k]
                   + f_0 * gh_243[k];

        t_160[k] = -2.0 * dh_97[k]
                   + f_0 * gh_244[k];

        t_161[k] = -2.0 * dh_98[k]
                   + f_0 * gh_245[k];

        t_162[k] = -2.0 * dh_99[k]
                   + f_0 * gh_246[k];

        t_163[k] = -2.0 * dh_100[k]
                   + f_0 * gh_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, dh_101, dh_102, dh_103, dh_104, \
                         dh_105, gh_248, gh_249, gh_250, gh_251, \
                         gh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -2.0 * dh_101[k]
                   + f_0 * gh_248[k];

        t_165[k] = -2.0 * dh_102[k]
                   + f_0 * gh_249[k];

        t_166[k] = -2.0 * dh_103[k]
                   + f_0 * gh_250[k];

        t_167[k] = -2.0 * dh_104[k]
                   + f_0 * gh_251[k];

        t_168[k] = -dh_105[k]
                   + f_0 * gh_252[k];
    }
}

static auto
compute_prim_geom_10_fh_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dh, const size_t gh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *dh_106 = buffer.data(dh + 106);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_109 = buffer.data(dh + 109);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_112 = buffer.data(dh + 112);
    const auto *dh_113 = buffer.data(dh + 113);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_115 = buffer.data(dh + 115);
    const auto *dh_116 = buffer.data(dh + 116);
    const auto *dh_117 = buffer.data(dh + 117);
    const auto *dh_118 = buffer.data(dh + 118);
    const auto *dh_119 = buffer.data(dh + 119);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

    const auto *gh_253 = buffer.data(gh + 253);
    const auto *gh_254 = buffer.data(gh + 254);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_256 = buffer.data(gh + 256);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_259 = buffer.data(gh + 259);
    const auto *gh_260 = buffer.data(gh + 260);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_262 = buffer.data(gh + 262);
    const auto *gh_263 = buffer.data(gh + 263);
    const auto *gh_264 = buffer.data(gh + 264);
    const auto *gh_265 = buffer.data(gh + 265);
    const auto *gh_266 = buffer.data(gh + 266);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_273 = buffer.data(gh + 273);
    const auto *gh_274 = buffer.data(gh + 274);
    const auto *gh_275 = buffer.data(gh + 275);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_277 = buffer.data(gh + 277);
    const auto *gh_278 = buffer.data(gh + 278);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_280 = buffer.data(gh + 280);
    const auto *gh_281 = buffer.data(gh + 281);
    const auto *gh_282 = buffer.data(gh + 282);
    const auto *gh_283 = buffer.data(gh + 283);
    const auto *gh_284 = buffer.data(gh + 284);
    const auto *gh_285 = buffer.data(gh + 285);
    const auto *gh_286 = buffer.data(gh + 286);
    const auto *gh_287 = buffer.data(gh + 287);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, dh_106, dh_107, dh_108, dh_109, \
                         dh_110, gh_253, gh_254, gh_255, gh_256, \
                         gh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -dh_106[k]
                   + f_0 * gh_253[k];

        t_170[k] = -dh_107[k]
                   + f_0 * gh_254[k];

        t_171[k] = -dh_108[k]
                   + f_0 * gh_255[k];

        t_172[k] = -dh_109[k]
                   + f_0 * gh_256[k];

        t_173[k] = -dh_110[k]
                   + f_0 * gh_257[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, dh_111, dh_112, dh_113, dh_114, \
                         dh_115, gh_258, gh_259, gh_260, gh_261, \
                         gh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -dh_111[k]
                   + f_0 * gh_258[k];

        t_175[k] = -dh_112[k]
                   + f_0 * gh_259[k];

        t_176[k] = -dh_113[k]
                   + f_0 * gh_260[k];

        t_177[k] = -dh_114[k]
                   + f_0 * gh_261[k];

        t_178[k] = -dh_115[k]
                   + f_0 * gh_262[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, dh_116, dh_117, dh_118, dh_119, \
                         dh_120, gh_263, gh_264, gh_265, gh_266, \
                         gh_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -dh_116[k]
                   + f_0 * gh_263[k];

        t_180[k] = -dh_117[k]
                   + f_0 * gh_264[k];

        t_181[k] = -dh_118[k]
                   + f_0 * gh_265[k];

        t_182[k] = -dh_119[k]
                   + f_0 * gh_266[k];

        t_183[k] = -dh_120[k]
                   + f_0 * gh_267[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, dh_121, dh_122, dh_123, dh_124, \
                         dh_125, gh_268, gh_269, gh_270, gh_271, \
                         gh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -dh_121[k]
                   + f_0 * gh_268[k];

        t_185[k] = -dh_122[k]
                   + f_0 * gh_269[k];

        t_186[k] = -dh_123[k]
                   + f_0 * gh_270[k];

        t_187[k] = -dh_124[k]
                   + f_0 * gh_271[k];

        t_188[k] = -dh_125[k]
                   + f_0 * gh_272[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, t_196, gh_273, \
                         gh_274, gh_275, gh_276, gh_277, gh_278, gh_279, \
                         gh_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * gh_273[k];

        t_190[k] = f_0 * gh_274[k];

        t_191[k] = f_0 * gh_275[k];

        t_192[k] = f_0 * gh_276[k];

        t_193[k] = f_0 * gh_277[k];

        t_194[k] = f_0 * gh_278[k];

        t_195[k] = f_0 * gh_279[k];

        t_196[k] = f_0 * gh_280[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, t_202, t_203, t_204, gh_281, \
                         gh_282, gh_283, gh_284, gh_285, gh_286, gh_287, \
                         gh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_0 * gh_281[k];

        t_198[k] = f_0 * gh_282[k];

        t_199[k] = f_0 * gh_283[k];

        t_200[k] = f_0 * gh_284[k];

        t_201[k] = f_0 * gh_285[k];

        t_202[k] = f_0 * gh_286[k];

        t_203[k] = f_0 * gh_287[k];

        t_204[k] = f_0 * gh_288[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, gh_289, gh_290, gh_291, gh_292, \
                         gh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_0 * gh_289[k];

        t_206[k] = f_0 * gh_290[k];

        t_207[k] = f_0 * gh_291[k];

        t_208[k] = f_0 * gh_292[k];

        t_209[k] = f_0 * gh_293[k];
    }
}

auto
compute_prim_geom_10_fh_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t dh, const size_t gh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fh_electron_repulsion_1_piece0(buffer, target, dh, gh, ncols, alpha);

    compute_prim_geom_10_fh_electron_repulsion_1_piece1(buffer, target, dh, gh, ncols, alpha);
}

static auto
compute_prim_geom_10_fh_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dh, const size_t gh,
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

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_41 = buffer.data(dh + 41);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_57 = buffer.data(dh + 57);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_61 = buffer.data(dh + 61);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_65 = buffer.data(dh + 65);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_67 = buffer.data(dh + 67);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_70 = buffer.data(dh + 70);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_74 = buffer.data(dh + 74);
    const auto *dh_75 = buffer.data(dh + 75);
    const auto *dh_76 = buffer.data(dh + 76);
    const auto *dh_77 = buffer.data(dh + 77);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_79 = buffer.data(dh + 79);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_84 = buffer.data(dh + 84);
    const auto *dh_85 = buffer.data(dh + 85);
    const auto *dh_86 = buffer.data(dh + 86);
    const auto *dh_87 = buffer.data(dh + 87);
    const auto *dh_88 = buffer.data(dh + 88);
    const auto *dh_89 = buffer.data(dh + 89);
    const auto *dh_90 = buffer.data(dh + 90);
    const auto *dh_91 = buffer.data(dh + 91);
    const auto *dh_92 = buffer.data(dh + 92);

    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_43 = buffer.data(gh + 43);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_52 = buffer.data(gh + 52);
    const auto *gh_53 = buffer.data(gh + 53);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_61 = buffer.data(gh + 61);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_92 = buffer.data(gh + 92);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_106 = buffer.data(gh + 106);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_113 = buffer.data(gh + 113);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_115 = buffer.data(gh + 115);
    const auto *gh_116 = buffer.data(gh + 116);
    const auto *gh_117 = buffer.data(gh + 117);
    const auto *gh_118 = buffer.data(gh + 118);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_120 = buffer.data(gh + 120);
    const auto *gh_121 = buffer.data(gh + 121);
    const auto *gh_122 = buffer.data(gh + 122);
    const auto *gh_123 = buffer.data(gh + 123);
    const auto *gh_124 = buffer.data(gh + 124);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_148 = buffer.data(gh + 148);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_151 = buffer.data(gh + 151);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_155 = buffer.data(gh + 155);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_157 = buffer.data(gh + 157);
    const auto *gh_158 = buffer.data(gh + 158);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_160 = buffer.data(gh + 160);
    const auto *gh_161 = buffer.data(gh + 161);
    const auto *gh_162 = buffer.data(gh + 162);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_166 = buffer.data(gh + 166);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_169 = buffer.data(gh + 169);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_172 = buffer.data(gh + 172);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_175 = buffer.data(gh + 175);
    const auto *gh_176 = buffer.data(gh + 176);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_178 = buffer.data(gh + 178);
    const auto *gh_179 = buffer.data(gh + 179);
    const auto *gh_180 = buffer.data(gh + 180);
    const auto *gh_181 = buffer.data(gh + 181);
    const auto *gh_182 = buffer.data(gh + 182);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_188 = buffer.data(gh + 188);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_190 = buffer.data(gh + 190);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_192 = buffer.data(gh + 192);
    const auto *gh_193 = buffer.data(gh + 193);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_195 = buffer.data(gh + 195);
    const auto *gh_196 = buffer.data(gh + 196);
    const auto *gh_197 = buffer.data(gh + 197);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_199 = buffer.data(gh + 199);
    const auto *gh_200 = buffer.data(gh + 200);
    const auto *gh_201 = buffer.data(gh + 201);
    const auto *gh_202 = buffer.data(gh + 202);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_204 = buffer.data(gh + 204);
    const auto *gh_205 = buffer.data(gh + 205);
    const auto *gh_206 = buffer.data(gh + 206);
    const auto *gh_207 = buffer.data(gh + 207);
    const auto *gh_208 = buffer.data(gh + 208);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_231 = buffer.data(gh + 231);
    const auto *gh_232 = buffer.data(gh + 232);
    const auto *gh_233 = buffer.data(gh + 233);
    const auto *gh_234 = buffer.data(gh + 234);
    const auto *gh_235 = buffer.data(gh + 235);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_237 = buffer.data(gh + 237);
    const auto *gh_238 = buffer.data(gh + 238);
    const auto *gh_239 = buffer.data(gh + 239);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_241 = buffer.data(gh + 241);
    const auto *gh_242 = buffer.data(gh + 242);
    const auto *gh_243 = buffer.data(gh + 243);
    const auto *gh_244 = buffer.data(gh + 244);
    const auto *gh_245 = buffer.data(gh + 245);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_253 = buffer.data(gh + 253);
    const auto *gh_254 = buffer.data(gh + 254);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_256 = buffer.data(gh + 256);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_259 = buffer.data(gh + 259);
    const auto *gh_260 = buffer.data(gh + 260);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_262 = buffer.data(gh + 262);
    const auto *gh_263 = buffer.data(gh + 263);
    const auto *gh_264 = buffer.data(gh + 264);
    const auto *gh_265 = buffer.data(gh + 265);
    const auto *gh_266 = buffer.data(gh + 266);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_273 = buffer.data(gh + 273);
    const auto *gh_274 = buffer.data(gh + 274);
    const auto *gh_275 = buffer.data(gh + 275);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_277 = buffer.data(gh + 277);
    const auto *gh_278 = buffer.data(gh + 278);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_280 = buffer.data(gh + 280);
    const auto *gh_281 = buffer.data(gh + 281);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gh_42, gh_43, gh_44, gh_45, \
                         gh_46, gh_47, gh_48, gh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gh_42[k];

        t_1[k] = f_0 * gh_43[k];

        t_2[k] = f_0 * gh_44[k];

        t_3[k] = f_0 * gh_45[k];

        t_4[k] = f_0 * gh_46[k];

        t_5[k] = f_0 * gh_47[k];

        t_6[k] = f_0 * gh_48[k];

        t_7[k] = f_0 * gh_49[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, gh_50, gh_51, gh_52, \
                         gh_53, gh_54, gh_55, gh_56, gh_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gh_50[k];

        t_9[k] = f_0 * gh_51[k];

        t_10[k] = f_0 * gh_52[k];

        t_11[k] = f_0 * gh_53[k];

        t_12[k] = f_0 * gh_54[k];

        t_13[k] = f_0 * gh_55[k];

        t_14[k] = f_0 * gh_56[k];

        t_15[k] = f_0 * gh_57[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, gh_58, gh_59, gh_60, \
                         gh_61, gh_62, gh_84, gh_85, gh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * gh_58[k];

        t_17[k] = f_0 * gh_59[k];

        t_18[k] = f_0 * gh_60[k];

        t_19[k] = f_0 * gh_61[k];

        t_20[k] = f_0 * gh_62[k];

        t_21[k] = f_0 * gh_84[k];

        t_22[k] = f_0 * gh_85[k];

        t_23[k] = f_0 * gh_86[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, gh_87, gh_88, gh_89, \
                         gh_90, gh_91, gh_92, gh_93, gh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * gh_87[k];

        t_25[k] = f_0 * gh_88[k];

        t_26[k] = f_0 * gh_89[k];

        t_27[k] = f_0 * gh_90[k];

        t_28[k] = f_0 * gh_91[k];

        t_29[k] = f_0 * gh_92[k];

        t_30[k] = f_0 * gh_93[k];

        t_31[k] = f_0 * gh_94[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, gh_95, gh_96, gh_97, \
                         gh_98, gh_99, gh_100, gh_101, gh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * gh_95[k];

        t_33[k] = f_0 * gh_96[k];

        t_34[k] = f_0 * gh_97[k];

        t_35[k] = f_0 * gh_98[k];

        t_36[k] = f_0 * gh_99[k];

        t_37[k] = f_0 * gh_100[k];

        t_38[k] = f_0 * gh_101[k];

        t_39[k] = f_0 * gh_102[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, dh_0, dh_1, dh_2, dh_3, gh_103, \
                         gh_104, gh_105, gh_106, gh_107, gh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * gh_103[k];

        t_41[k] = f_0 * gh_104[k];

        t_42[k] = -dh_0[k]
                  + f_0 * gh_105[k];

        t_43[k] = -dh_1[k]
                  + f_0 * gh_106[k];

        t_44[k] = -dh_2[k]
                  + f_0 * gh_107[k];

        t_45[k] = -dh_3[k]
                  + f_0 * gh_108[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, dh_4, dh_5, dh_6, dh_7, dh_8, gh_109, \
                         gh_110, gh_111, gh_112, gh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -dh_4[k]
                  + f_0 * gh_109[k];

        t_47[k] = -dh_5[k]
                  + f_0 * gh_110[k];

        t_48[k] = -dh_6[k]
                  + f_0 * gh_111[k];

        t_49[k] = -dh_7[k]
                  + f_0 * gh_112[k];

        t_50[k] = -dh_8[k]
                  + f_0 * gh_113[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, dh_9, dh_10, dh_11, dh_12, dh_13, \
                         gh_114, gh_115, gh_116, gh_117, gh_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -dh_9[k]
                  + f_0 * gh_114[k];

        t_52[k] = -dh_10[k]
                  + f_0 * gh_115[k];

        t_53[k] = -dh_11[k]
                  + f_0 * gh_116[k];

        t_54[k] = -dh_12[k]
                  + f_0 * gh_117[k];

        t_55[k] = -dh_13[k]
                  + f_0 * gh_118[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, dh_14, dh_15, dh_16, dh_17, dh_18, \
                         gh_119, gh_120, gh_121, gh_122, gh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -dh_14[k]
                  + f_0 * gh_119[k];

        t_57[k] = -dh_15[k]
                  + f_0 * gh_120[k];

        t_58[k] = -dh_16[k]
                  + f_0 * gh_121[k];

        t_59[k] = -dh_17[k]
                  + f_0 * gh_122[k];

        t_60[k] = -dh_18[k]
                  + f_0 * gh_123[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, t_67, dh_19, dh_20, gh_124, \
                         gh_125, gh_147, gh_148, gh_149, gh_150, \
                         gh_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -dh_19[k]
                  + f_0 * gh_124[k];

        t_62[k] = -dh_20[k]
                  + f_0 * gh_125[k];

        t_63[k] = f_0 * gh_147[k];

        t_64[k] = f_0 * gh_148[k];

        t_65[k] = f_0 * gh_149[k];

        t_66[k] = f_0 * gh_150[k];

        t_67[k] = f_0 * gh_151[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, t_74, t_75, gh_152, gh_153, \
                         gh_154, gh_155, gh_156, gh_157, gh_158, \
                         gh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * gh_152[k];

        t_69[k] = f_0 * gh_153[k];

        t_70[k] = f_0 * gh_154[k];

        t_71[k] = f_0 * gh_155[k];

        t_72[k] = f_0 * gh_156[k];

        t_73[k] = f_0 * gh_157[k];

        t_74[k] = f_0 * gh_158[k];

        t_75[k] = f_0 * gh_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, t_82, t_83, gh_160, gh_161, \
                         gh_162, gh_163, gh_164, gh_165, gh_166, \
                         gh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_0 * gh_160[k];

        t_77[k] = f_0 * gh_161[k];

        t_78[k] = f_0 * gh_162[k];

        t_79[k] = f_0 * gh_163[k];

        t_80[k] = f_0 * gh_164[k];

        t_81[k] = f_0 * gh_165[k];

        t_82[k] = f_0 * gh_166[k];

        t_83[k] = f_0 * gh_167[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, dh_21, dh_22, dh_23, dh_24, dh_25, \
                         gh_168, gh_169, gh_170, gh_171, gh_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -dh_21[k]
                  + f_0 * gh_168[k];

        t_85[k] = -dh_22[k]
                  + f_0 * gh_169[k];

        t_86[k] = -dh_23[k]
                  + f_0 * gh_170[k];

        t_87[k] = -dh_24[k]
                  + f_0 * gh_171[k];

        t_88[k] = -dh_25[k]
                  + f_0 * gh_172[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, dh_26, dh_27, dh_28, dh_29, dh_30, \
                         gh_173, gh_174, gh_175, gh_176, gh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -dh_26[k]
                  + f_0 * gh_173[k];

        t_90[k] = -dh_27[k]
                  + f_0 * gh_174[k];

        t_91[k] = -dh_28[k]
                  + f_0 * gh_175[k];

        t_92[k] = -dh_29[k]
                  + f_0 * gh_176[k];

        t_93[k] = -dh_30[k]
                  + f_0 * gh_177[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, dh_31, dh_32, dh_33, dh_34, dh_35, \
                         gh_178, gh_179, gh_180, gh_181, gh_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -dh_31[k]
                  + f_0 * gh_178[k];

        t_95[k] = -dh_32[k]
                  + f_0 * gh_179[k];

        t_96[k] = -dh_33[k]
                  + f_0 * gh_180[k];

        t_97[k] = -dh_34[k]
                  + f_0 * gh_181[k];

        t_98[k] = -dh_35[k]
                  + f_0 * gh_182[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, dh_36, dh_37, dh_38, dh_39, dh_40, \
                         gh_183, gh_184, gh_185, gh_186, gh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -dh_36[k]
                  + f_0 * gh_183[k];

        t_100[k] = -dh_37[k]
                   + f_0 * gh_184[k];

        t_101[k] = -dh_38[k]
                   + f_0 * gh_185[k];

        t_102[k] = -dh_39[k]
                   + f_0 * gh_186[k];

        t_103[k] = -dh_40[k]
                   + f_0 * gh_187[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, dh_41, dh_42, dh_43, dh_44, dh_45, \
                         gh_188, gh_189, gh_190, gh_191, gh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -dh_41[k]
                   + f_0 * gh_188[k];

        t_105[k] = -2.0 * dh_42[k]
                   + f_0 * gh_189[k];

        t_106[k] = -2.0 * dh_43[k]
                   + f_0 * gh_190[k];

        t_107[k] = -2.0 * dh_44[k]
                   + f_0 * gh_191[k];

        t_108[k] = -2.0 * dh_45[k]
                   + f_0 * gh_192[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, dh_46, dh_47, dh_48, dh_49, dh_50, \
                         gh_193, gh_194, gh_195, gh_196, gh_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -2.0 * dh_46[k]
                   + f_0 * gh_193[k];

        t_110[k] = -2.0 * dh_47[k]
                   + f_0 * gh_194[k];

        t_111[k] = -2.0 * dh_48[k]
                   + f_0 * gh_195[k];

        t_112[k] = -2.0 * dh_49[k]
                   + f_0 * gh_196[k];

        t_113[k] = -2.0 * dh_50[k]
                   + f_0 * gh_197[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, dh_51, dh_52, dh_53, dh_54, dh_55, \
                         gh_198, gh_199, gh_200, gh_201, gh_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -2.0 * dh_51[k]
                   + f_0 * gh_198[k];

        t_115[k] = -2.0 * dh_52[k]
                   + f_0 * gh_199[k];

        t_116[k] = -2.0 * dh_53[k]
                   + f_0 * gh_200[k];

        t_117[k] = -2.0 * dh_54[k]
                   + f_0 * gh_201[k];

        t_118[k] = -2.0 * dh_55[k]
                   + f_0 * gh_202[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, dh_56, dh_57, dh_58, dh_59, dh_60, \
                         gh_203, gh_204, gh_205, gh_206, gh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -2.0 * dh_56[k]
                   + f_0 * gh_203[k];

        t_120[k] = -2.0 * dh_57[k]
                   + f_0 * gh_204[k];

        t_121[k] = -2.0 * dh_58[k]
                   + f_0 * gh_205[k];

        t_122[k] = -2.0 * dh_59[k]
                   + f_0 * gh_206[k];

        t_123[k] = -2.0 * dh_60[k]
                   + f_0 * gh_207[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, t_129, t_130, dh_61, dh_62, \
                         gh_208, gh_209, gh_231, gh_232, gh_233, gh_234, \
                         gh_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * dh_61[k]
                   + f_0 * gh_208[k];

        t_125[k] = -2.0 * dh_62[k]
                   + f_0 * gh_209[k];

        t_126[k] = f_0 * gh_231[k];

        t_127[k] = f_0 * gh_232[k];

        t_128[k] = f_0 * gh_233[k];

        t_129[k] = f_0 * gh_234[k];

        t_130[k] = f_0 * gh_235[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, t_136, t_137, t_138, gh_236, \
                         gh_237, gh_238, gh_239, gh_240, gh_241, gh_242, \
                         gh_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_0 * gh_236[k];

        t_132[k] = f_0 * gh_237[k];

        t_133[k] = f_0 * gh_238[k];

        t_134[k] = f_0 * gh_239[k];

        t_135[k] = f_0 * gh_240[k];

        t_136[k] = f_0 * gh_241[k];

        t_137[k] = f_0 * gh_242[k];

        t_138[k] = f_0 * gh_243[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, t_145, t_146, gh_244, \
                         gh_245, gh_246, gh_247, gh_248, gh_249, gh_250, \
                         gh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_0 * gh_244[k];

        t_140[k] = f_0 * gh_245[k];

        t_141[k] = f_0 * gh_246[k];

        t_142[k] = f_0 * gh_247[k];

        t_143[k] = f_0 * gh_248[k];

        t_144[k] = f_0 * gh_249[k];

        t_145[k] = f_0 * gh_250[k];

        t_146[k] = f_0 * gh_251[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, dh_63, dh_64, dh_65, dh_66, dh_67, \
                         gh_252, gh_253, gh_254, gh_255, gh_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -dh_63[k]
                   + f_0 * gh_252[k];

        t_148[k] = -dh_64[k]
                   + f_0 * gh_253[k];

        t_149[k] = -dh_65[k]
                   + f_0 * gh_254[k];

        t_150[k] = -dh_66[k]
                   + f_0 * gh_255[k];

        t_151[k] = -dh_67[k]
                   + f_0 * gh_256[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, dh_68, dh_69, dh_70, dh_71, dh_72, \
                         gh_257, gh_258, gh_259, gh_260, gh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -dh_68[k]
                   + f_0 * gh_257[k];

        t_153[k] = -dh_69[k]
                   + f_0 * gh_258[k];

        t_154[k] = -dh_70[k]
                   + f_0 * gh_259[k];

        t_155[k] = -dh_71[k]
                   + f_0 * gh_260[k];

        t_156[k] = -dh_72[k]
                   + f_0 * gh_261[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, dh_73, dh_74, dh_75, dh_76, dh_77, \
                         gh_262, gh_263, gh_264, gh_265, gh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -dh_73[k]
                   + f_0 * gh_262[k];

        t_158[k] = -dh_74[k]
                   + f_0 * gh_263[k];

        t_159[k] = -dh_75[k]
                   + f_0 * gh_264[k];

        t_160[k] = -dh_76[k]
                   + f_0 * gh_265[k];

        t_161[k] = -dh_77[k]
                   + f_0 * gh_266[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, dh_78, dh_79, dh_80, dh_81, dh_82, \
                         gh_267, gh_268, gh_269, gh_270, gh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -dh_78[k]
                   + f_0 * gh_267[k];

        t_163[k] = -dh_79[k]
                   + f_0 * gh_268[k];

        t_164[k] = -dh_80[k]
                   + f_0 * gh_269[k];

        t_165[k] = -dh_81[k]
                   + f_0 * gh_270[k];

        t_166[k] = -dh_82[k]
                   + f_0 * gh_271[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, dh_83, dh_84, dh_85, dh_86, dh_87, \
                         gh_272, gh_273, gh_274, gh_275, gh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -dh_83[k]
                   + f_0 * gh_272[k];

        t_168[k] = -2.0 * dh_84[k]
                   + f_0 * gh_273[k];

        t_169[k] = -2.0 * dh_85[k]
                   + f_0 * gh_274[k];

        t_170[k] = -2.0 * dh_86[k]
                   + f_0 * gh_275[k];

        t_171[k] = -2.0 * dh_87[k]
                   + f_0 * gh_276[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, dh_88, dh_89, dh_90, dh_91, dh_92, \
                         gh_277, gh_278, gh_279, gh_280, gh_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -2.0 * dh_88[k]
                   + f_0 * gh_277[k];

        t_173[k] = -2.0 * dh_89[k]
                   + f_0 * gh_278[k];

        t_174[k] = -2.0 * dh_90[k]
                   + f_0 * gh_279[k];

        t_175[k] = -2.0 * dh_91[k]
                   + f_0 * gh_280[k];

        t_176[k] = -2.0 * dh_92[k]
                   + f_0 * gh_281[k];
    }
}

static auto
compute_prim_geom_10_fh_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dh, const size_t gh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *dh_93 = buffer.data(dh + 93);
    const auto *dh_94 = buffer.data(dh + 94);
    const auto *dh_95 = buffer.data(dh + 95);
    const auto *dh_96 = buffer.data(dh + 96);
    const auto *dh_97 = buffer.data(dh + 97);
    const auto *dh_98 = buffer.data(dh + 98);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_106 = buffer.data(dh + 106);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_109 = buffer.data(dh + 109);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_112 = buffer.data(dh + 112);
    const auto *dh_113 = buffer.data(dh + 113);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_115 = buffer.data(dh + 115);
    const auto *dh_116 = buffer.data(dh + 116);
    const auto *dh_117 = buffer.data(dh + 117);
    const auto *dh_118 = buffer.data(dh + 118);
    const auto *dh_119 = buffer.data(dh + 119);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

    const auto *gh_282 = buffer.data(gh + 282);
    const auto *gh_283 = buffer.data(gh + 283);
    const auto *gh_284 = buffer.data(gh + 284);
    const auto *gh_285 = buffer.data(gh + 285);
    const auto *gh_286 = buffer.data(gh + 286);
    const auto *gh_287 = buffer.data(gh + 287);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_295 = buffer.data(gh + 295);
    const auto *gh_296 = buffer.data(gh + 296);
    const auto *gh_297 = buffer.data(gh + 297);
    const auto *gh_298 = buffer.data(gh + 298);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_300 = buffer.data(gh + 300);
    const auto *gh_301 = buffer.data(gh + 301);
    const auto *gh_302 = buffer.data(gh + 302);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_304 = buffer.data(gh + 304);
    const auto *gh_305 = buffer.data(gh + 305);
    const auto *gh_306 = buffer.data(gh + 306);
    const auto *gh_307 = buffer.data(gh + 307);
    const auto *gh_308 = buffer.data(gh + 308);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_313 = buffer.data(gh + 313);
    const auto *gh_314 = buffer.data(gh + 314);

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, dh_93, dh_94, dh_95, dh_96, dh_97, \
                         gh_282, gh_283, gh_284, gh_285, gh_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -2.0 * dh_93[k]
                   + f_0 * gh_282[k];

        t_178[k] = -2.0 * dh_94[k]
                   + f_0 * gh_283[k];

        t_179[k] = -2.0 * dh_95[k]
                   + f_0 * gh_284[k];

        t_180[k] = -2.0 * dh_96[k]
                   + f_0 * gh_285[k];

        t_181[k] = -2.0 * dh_97[k]
                   + f_0 * gh_286[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, dh_98, dh_99, dh_100, dh_101, \
                         dh_102, gh_287, gh_288, gh_289, gh_290, \
                         gh_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * dh_98[k]
                   + f_0 * gh_287[k];

        t_183[k] = -2.0 * dh_99[k]
                   + f_0 * gh_288[k];

        t_184[k] = -2.0 * dh_100[k]
                   + f_0 * gh_289[k];

        t_185[k] = -2.0 * dh_101[k]
                   + f_0 * gh_290[k];

        t_186[k] = -2.0 * dh_102[k]
                   + f_0 * gh_291[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, dh_103, dh_104, dh_105, dh_106, \
                         dh_107, gh_292, gh_293, gh_294, gh_295, \
                         gh_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * dh_103[k]
                   + f_0 * gh_292[k];

        t_188[k] = -2.0 * dh_104[k]
                   + f_0 * gh_293[k];

        t_189[k] = -3.0 * dh_105[k]
                   + f_0 * gh_294[k];

        t_190[k] = -3.0 * dh_106[k]
                   + f_0 * gh_295[k];

        t_191[k] = -3.0 * dh_107[k]
                   + f_0 * gh_296[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, dh_108, dh_109, dh_110, dh_111, \
                         dh_112, gh_297, gh_298, gh_299, gh_300, \
                         gh_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -3.0 * dh_108[k]
                   + f_0 * gh_297[k];

        t_193[k] = -3.0 * dh_109[k]
                   + f_0 * gh_298[k];

        t_194[k] = -3.0 * dh_110[k]
                   + f_0 * gh_299[k];

        t_195[k] = -3.0 * dh_111[k]
                   + f_0 * gh_300[k];

        t_196[k] = -3.0 * dh_112[k]
                   + f_0 * gh_301[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, dh_113, dh_114, dh_115, dh_116, \
                         dh_117, gh_302, gh_303, gh_304, gh_305, \
                         gh_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -3.0 * dh_113[k]
                   + f_0 * gh_302[k];

        t_198[k] = -3.0 * dh_114[k]
                   + f_0 * gh_303[k];

        t_199[k] = -3.0 * dh_115[k]
                   + f_0 * gh_304[k];

        t_200[k] = -3.0 * dh_116[k]
                   + f_0 * gh_305[k];

        t_201[k] = -3.0 * dh_117[k]
                   + f_0 * gh_306[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, dh_118, dh_119, dh_120, dh_121, \
                         dh_122, gh_307, gh_308, gh_309, gh_310, \
                         gh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -3.0 * dh_118[k]
                   + f_0 * gh_307[k];

        t_203[k] = -3.0 * dh_119[k]
                   + f_0 * gh_308[k];

        t_204[k] = -3.0 * dh_120[k]
                   + f_0 * gh_309[k];

        t_205[k] = -3.0 * dh_121[k]
                   + f_0 * gh_310[k];

        t_206[k] = -3.0 * dh_122[k]
                   + f_0 * gh_311[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, dh_123, dh_124, dh_125, gh_312, gh_313, \
                         gh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -3.0 * dh_123[k]
                   + f_0 * gh_312[k];

        t_208[k] = -3.0 * dh_124[k]
                   + f_0 * gh_313[k];

        t_209[k] = -3.0 * dh_125[k]
                   + f_0 * gh_314[k];
    }
}

auto
compute_prim_geom_10_fh_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t dh, const size_t gh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fh_electron_repulsion_2_piece0(buffer, target, dh, gh, ncols, alpha);

    compute_prim_geom_10_fh_electron_repulsion_2_piece1(buffer, target, dh, gh, ncols, alpha);
}

}  // namespace simdt2ceri
