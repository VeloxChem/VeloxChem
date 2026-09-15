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


#include "SimdElectronRepulsionGeom10VrrRecDL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_dl_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t pl, const size_t fl,
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

    const auto *pl_0 = buffer.data(pl + 0);
    const auto *pl_1 = buffer.data(pl + 1);
    const auto *pl_2 = buffer.data(pl + 2);
    const auto *pl_3 = buffer.data(pl + 3);
    const auto *pl_4 = buffer.data(pl + 4);
    const auto *pl_5 = buffer.data(pl + 5);
    const auto *pl_6 = buffer.data(pl + 6);
    const auto *pl_7 = buffer.data(pl + 7);
    const auto *pl_8 = buffer.data(pl + 8);
    const auto *pl_9 = buffer.data(pl + 9);
    const auto *pl_10 = buffer.data(pl + 10);
    const auto *pl_11 = buffer.data(pl + 11);
    const auto *pl_12 = buffer.data(pl + 12);
    const auto *pl_13 = buffer.data(pl + 13);
    const auto *pl_14 = buffer.data(pl + 14);
    const auto *pl_15 = buffer.data(pl + 15);
    const auto *pl_16 = buffer.data(pl + 16);
    const auto *pl_17 = buffer.data(pl + 17);
    const auto *pl_18 = buffer.data(pl + 18);
    const auto *pl_19 = buffer.data(pl + 19);
    const auto *pl_20 = buffer.data(pl + 20);
    const auto *pl_21 = buffer.data(pl + 21);
    const auto *pl_22 = buffer.data(pl + 22);
    const auto *pl_23 = buffer.data(pl + 23);
    const auto *pl_24 = buffer.data(pl + 24);
    const auto *pl_25 = buffer.data(pl + 25);
    const auto *pl_26 = buffer.data(pl + 26);
    const auto *pl_27 = buffer.data(pl + 27);
    const auto *pl_28 = buffer.data(pl + 28);
    const auto *pl_29 = buffer.data(pl + 29);
    const auto *pl_30 = buffer.data(pl + 30);
    const auto *pl_31 = buffer.data(pl + 31);
    const auto *pl_32 = buffer.data(pl + 32);
    const auto *pl_33 = buffer.data(pl + 33);
    const auto *pl_34 = buffer.data(pl + 34);
    const auto *pl_35 = buffer.data(pl + 35);
    const auto *pl_36 = buffer.data(pl + 36);
    const auto *pl_37 = buffer.data(pl + 37);
    const auto *pl_38 = buffer.data(pl + 38);
    const auto *pl_39 = buffer.data(pl + 39);
    const auto *pl_40 = buffer.data(pl + 40);
    const auto *pl_41 = buffer.data(pl + 41);
    const auto *pl_42 = buffer.data(pl + 42);
    const auto *pl_43 = buffer.data(pl + 43);
    const auto *pl_44 = buffer.data(pl + 44);
    const auto *pl_45 = buffer.data(pl + 45);
    const auto *pl_46 = buffer.data(pl + 46);
    const auto *pl_47 = buffer.data(pl + 47);
    const auto *pl_48 = buffer.data(pl + 48);
    const auto *pl_49 = buffer.data(pl + 49);
    const auto *pl_50 = buffer.data(pl + 50);
    const auto *pl_51 = buffer.data(pl + 51);
    const auto *pl_52 = buffer.data(pl + 52);
    const auto *pl_53 = buffer.data(pl + 53);
    const auto *pl_54 = buffer.data(pl + 54);
    const auto *pl_55 = buffer.data(pl + 55);
    const auto *pl_56 = buffer.data(pl + 56);
    const auto *pl_57 = buffer.data(pl + 57);
    const auto *pl_58 = buffer.data(pl + 58);
    const auto *pl_59 = buffer.data(pl + 59);
    const auto *pl_60 = buffer.data(pl + 60);
    const auto *pl_61 = buffer.data(pl + 61);
    const auto *pl_62 = buffer.data(pl + 62);
    const auto *pl_63 = buffer.data(pl + 63);
    const auto *pl_64 = buffer.data(pl + 64);
    const auto *pl_65 = buffer.data(pl + 65);
    const auto *pl_66 = buffer.data(pl + 66);
    const auto *pl_67 = buffer.data(pl + 67);
    const auto *pl_68 = buffer.data(pl + 68);
    const auto *pl_69 = buffer.data(pl + 69);
    const auto *pl_70 = buffer.data(pl + 70);
    const auto *pl_71 = buffer.data(pl + 71);
    const auto *pl_72 = buffer.data(pl + 72);
    const auto *pl_73 = buffer.data(pl + 73);
    const auto *pl_74 = buffer.data(pl + 74);
    const auto *pl_75 = buffer.data(pl + 75);
    const auto *pl_76 = buffer.data(pl + 76);
    const auto *pl_77 = buffer.data(pl + 77);
    const auto *pl_78 = buffer.data(pl + 78);
    const auto *pl_79 = buffer.data(pl + 79);
    const auto *pl_80 = buffer.data(pl + 80);
    const auto *pl_81 = buffer.data(pl + 81);
    const auto *pl_82 = buffer.data(pl + 82);
    const auto *pl_83 = buffer.data(pl + 83);
    const auto *pl_84 = buffer.data(pl + 84);
    const auto *pl_85 = buffer.data(pl + 85);
    const auto *pl_86 = buffer.data(pl + 86);
    const auto *pl_87 = buffer.data(pl + 87);
    const auto *pl_88 = buffer.data(pl + 88);
    const auto *pl_89 = buffer.data(pl + 89);
    const auto *pl_90 = buffer.data(pl + 90);
    const auto *pl_91 = buffer.data(pl + 91);
    const auto *pl_92 = buffer.data(pl + 92);
    const auto *pl_93 = buffer.data(pl + 93);
    const auto *pl_94 = buffer.data(pl + 94);
    const auto *pl_95 = buffer.data(pl + 95);
    const auto *pl_96 = buffer.data(pl + 96);
    const auto *pl_97 = buffer.data(pl + 97);
    const auto *pl_98 = buffer.data(pl + 98);
    const auto *pl_99 = buffer.data(pl + 99);
    const auto *pl_100 = buffer.data(pl + 100);
    const auto *pl_101 = buffer.data(pl + 101);
    const auto *pl_102 = buffer.data(pl + 102);
    const auto *pl_103 = buffer.data(pl + 103);
    const auto *pl_104 = buffer.data(pl + 104);
    const auto *pl_105 = buffer.data(pl + 105);
    const auto *pl_106 = buffer.data(pl + 106);
    const auto *pl_107 = buffer.data(pl + 107);
    const auto *pl_108 = buffer.data(pl + 108);
    const auto *pl_109 = buffer.data(pl + 109);
    const auto *pl_110 = buffer.data(pl + 110);
    const auto *pl_111 = buffer.data(pl + 111);
    const auto *pl_112 = buffer.data(pl + 112);
    const auto *pl_113 = buffer.data(pl + 113);
    const auto *pl_114 = buffer.data(pl + 114);
    const auto *pl_115 = buffer.data(pl + 115);
    const auto *pl_116 = buffer.data(pl + 116);
    const auto *pl_117 = buffer.data(pl + 117);
    const auto *pl_118 = buffer.data(pl + 118);
    const auto *pl_119 = buffer.data(pl + 119);
    const auto *pl_120 = buffer.data(pl + 120);
    const auto *pl_121 = buffer.data(pl + 121);
    const auto *pl_122 = buffer.data(pl + 122);
    const auto *pl_123 = buffer.data(pl + 123);
    const auto *pl_124 = buffer.data(pl + 124);
    const auto *pl_125 = buffer.data(pl + 125);
    const auto *pl_126 = buffer.data(pl + 126);
    const auto *pl_127 = buffer.data(pl + 127);
    const auto *pl_128 = buffer.data(pl + 128);
    const auto *pl_129 = buffer.data(pl + 129);
    const auto *pl_130 = buffer.data(pl + 130);
    const auto *pl_131 = buffer.data(pl + 131);
    const auto *pl_132 = buffer.data(pl + 132);
    const auto *pl_133 = buffer.data(pl + 133);
    const auto *pl_134 = buffer.data(pl + 134);

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
    const auto *fl_150 = buffer.data(fl + 150);
    const auto *fl_151 = buffer.data(fl + 151);
    const auto *fl_152 = buffer.data(fl + 152);
    const auto *fl_153 = buffer.data(fl + 153);
    const auto *fl_154 = buffer.data(fl + 154);
    const auto *fl_155 = buffer.data(fl + 155);
    const auto *fl_156 = buffer.data(fl + 156);
    const auto *fl_157 = buffer.data(fl + 157);
    const auto *fl_158 = buffer.data(fl + 158);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pl_0, pl_1, pl_2, pl_3, pl_4, fl_0, fl_1, \
                         fl_2, fl_3, fl_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -2.0 * pl_0[k]
                 + f_0 * fl_0[k];

        t_1[k] = -2.0 * pl_1[k]
                 + f_0 * fl_1[k];

        t_2[k] = -2.0 * pl_2[k]
                 + f_0 * fl_2[k];

        t_3[k] = -2.0 * pl_3[k]
                 + f_0 * fl_3[k];

        t_4[k] = -2.0 * pl_4[k]
                 + f_0 * fl_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pl_5, pl_6, pl_7, pl_8, pl_9, fl_5, fl_6, \
                         fl_7, fl_8, fl_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -2.0 * pl_5[k]
                 + f_0 * fl_5[k];

        t_6[k] = -2.0 * pl_6[k]
                 + f_0 * fl_6[k];

        t_7[k] = -2.0 * pl_7[k]
                 + f_0 * fl_7[k];

        t_8[k] = -2.0 * pl_8[k]
                 + f_0 * fl_8[k];

        t_9[k] = -2.0 * pl_9[k]
                 + f_0 * fl_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pl_10, pl_11, pl_12, pl_13, pl_14, \
                         fl_10, fl_11, fl_12, fl_13, fl_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -2.0 * pl_10[k]
                  + f_0 * fl_10[k];

        t_11[k] = -2.0 * pl_11[k]
                  + f_0 * fl_11[k];

        t_12[k] = -2.0 * pl_12[k]
                  + f_0 * fl_12[k];

        t_13[k] = -2.0 * pl_13[k]
                  + f_0 * fl_13[k];

        t_14[k] = -2.0 * pl_14[k]
                  + f_0 * fl_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pl_15, pl_16, pl_17, pl_18, pl_19, \
                         fl_15, fl_16, fl_17, fl_18, fl_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -2.0 * pl_15[k]
                  + f_0 * fl_15[k];

        t_16[k] = -2.0 * pl_16[k]
                  + f_0 * fl_16[k];

        t_17[k] = -2.0 * pl_17[k]
                  + f_0 * fl_17[k];

        t_18[k] = -2.0 * pl_18[k]
                  + f_0 * fl_18[k];

        t_19[k] = -2.0 * pl_19[k]
                  + f_0 * fl_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pl_20, pl_21, pl_22, pl_23, pl_24, \
                         fl_20, fl_21, fl_22, fl_23, fl_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -2.0 * pl_20[k]
                  + f_0 * fl_20[k];

        t_21[k] = -2.0 * pl_21[k]
                  + f_0 * fl_21[k];

        t_22[k] = -2.0 * pl_22[k]
                  + f_0 * fl_22[k];

        t_23[k] = -2.0 * pl_23[k]
                  + f_0 * fl_23[k];

        t_24[k] = -2.0 * pl_24[k]
                  + f_0 * fl_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pl_25, pl_26, pl_27, pl_28, pl_29, \
                         fl_25, fl_26, fl_27, fl_28, fl_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -2.0 * pl_25[k]
                  + f_0 * fl_25[k];

        t_26[k] = -2.0 * pl_26[k]
                  + f_0 * fl_26[k];

        t_27[k] = -2.0 * pl_27[k]
                  + f_0 * fl_27[k];

        t_28[k] = -2.0 * pl_28[k]
                  + f_0 * fl_28[k];

        t_29[k] = -2.0 * pl_29[k]
                  + f_0 * fl_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pl_30, pl_31, pl_32, pl_33, pl_34, \
                         fl_30, fl_31, fl_32, fl_33, fl_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -2.0 * pl_30[k]
                  + f_0 * fl_30[k];

        t_31[k] = -2.0 * pl_31[k]
                  + f_0 * fl_31[k];

        t_32[k] = -2.0 * pl_32[k]
                  + f_0 * fl_32[k];

        t_33[k] = -2.0 * pl_33[k]
                  + f_0 * fl_33[k];

        t_34[k] = -2.0 * pl_34[k]
                  + f_0 * fl_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pl_35, pl_36, pl_37, pl_38, pl_39, \
                         fl_35, fl_36, fl_37, fl_38, fl_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -2.0 * pl_35[k]
                  + f_0 * fl_35[k];

        t_36[k] = -2.0 * pl_36[k]
                  + f_0 * fl_36[k];

        t_37[k] = -2.0 * pl_37[k]
                  + f_0 * fl_37[k];

        t_38[k] = -2.0 * pl_38[k]
                  + f_0 * fl_38[k];

        t_39[k] = -2.0 * pl_39[k]
                  + f_0 * fl_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pl_40, pl_41, pl_42, pl_43, pl_44, \
                         fl_40, fl_41, fl_42, fl_43, fl_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -2.0 * pl_40[k]
                  + f_0 * fl_40[k];

        t_41[k] = -2.0 * pl_41[k]
                  + f_0 * fl_41[k];

        t_42[k] = -2.0 * pl_42[k]
                  + f_0 * fl_42[k];

        t_43[k] = -2.0 * pl_43[k]
                  + f_0 * fl_43[k];

        t_44[k] = -2.0 * pl_44[k]
                  + f_0 * fl_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pl_45, pl_46, pl_47, pl_48, pl_49, \
                         fl_45, fl_46, fl_47, fl_48, fl_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -pl_45[k]
                  + f_0 * fl_45[k];

        t_46[k] = -pl_46[k]
                  + f_0 * fl_46[k];

        t_47[k] = -pl_47[k]
                  + f_0 * fl_47[k];

        t_48[k] = -pl_48[k]
                  + f_0 * fl_48[k];

        t_49[k] = -pl_49[k]
                  + f_0 * fl_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pl_50, pl_51, pl_52, pl_53, pl_54, \
                         fl_50, fl_51, fl_52, fl_53, fl_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -pl_50[k]
                  + f_0 * fl_50[k];

        t_51[k] = -pl_51[k]
                  + f_0 * fl_51[k];

        t_52[k] = -pl_52[k]
                  + f_0 * fl_52[k];

        t_53[k] = -pl_53[k]
                  + f_0 * fl_53[k];

        t_54[k] = -pl_54[k]
                  + f_0 * fl_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pl_55, pl_56, pl_57, pl_58, pl_59, \
                         fl_55, fl_56, fl_57, fl_58, fl_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -pl_55[k]
                  + f_0 * fl_55[k];

        t_56[k] = -pl_56[k]
                  + f_0 * fl_56[k];

        t_57[k] = -pl_57[k]
                  + f_0 * fl_57[k];

        t_58[k] = -pl_58[k]
                  + f_0 * fl_58[k];

        t_59[k] = -pl_59[k]
                  + f_0 * fl_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pl_60, pl_61, pl_62, pl_63, pl_64, \
                         fl_60, fl_61, fl_62, fl_63, fl_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -pl_60[k]
                  + f_0 * fl_60[k];

        t_61[k] = -pl_61[k]
                  + f_0 * fl_61[k];

        t_62[k] = -pl_62[k]
                  + f_0 * fl_62[k];

        t_63[k] = -pl_63[k]
                  + f_0 * fl_63[k];

        t_64[k] = -pl_64[k]
                  + f_0 * fl_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pl_65, pl_66, pl_67, pl_68, pl_69, \
                         fl_65, fl_66, fl_67, fl_68, fl_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -pl_65[k]
                  + f_0 * fl_65[k];

        t_66[k] = -pl_66[k]
                  + f_0 * fl_66[k];

        t_67[k] = -pl_67[k]
                  + f_0 * fl_67[k];

        t_68[k] = -pl_68[k]
                  + f_0 * fl_68[k];

        t_69[k] = -pl_69[k]
                  + f_0 * fl_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pl_70, pl_71, pl_72, pl_73, pl_74, \
                         fl_70, fl_71, fl_72, fl_73, fl_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -pl_70[k]
                  + f_0 * fl_70[k];

        t_71[k] = -pl_71[k]
                  + f_0 * fl_71[k];

        t_72[k] = -pl_72[k]
                  + f_0 * fl_72[k];

        t_73[k] = -pl_73[k]
                  + f_0 * fl_73[k];

        t_74[k] = -pl_74[k]
                  + f_0 * fl_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pl_75, pl_76, pl_77, pl_78, pl_79, \
                         fl_75, fl_76, fl_77, fl_78, fl_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -pl_75[k]
                  + f_0 * fl_75[k];

        t_76[k] = -pl_76[k]
                  + f_0 * fl_76[k];

        t_77[k] = -pl_77[k]
                  + f_0 * fl_77[k];

        t_78[k] = -pl_78[k]
                  + f_0 * fl_78[k];

        t_79[k] = -pl_79[k]
                  + f_0 * fl_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pl_80, pl_81, pl_82, pl_83, pl_84, \
                         fl_80, fl_81, fl_82, fl_83, fl_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -pl_80[k]
                  + f_0 * fl_80[k];

        t_81[k] = -pl_81[k]
                  + f_0 * fl_81[k];

        t_82[k] = -pl_82[k]
                  + f_0 * fl_82[k];

        t_83[k] = -pl_83[k]
                  + f_0 * fl_83[k];

        t_84[k] = -pl_84[k]
                  + f_0 * fl_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pl_85, pl_86, pl_87, pl_88, pl_89, \
                         fl_85, fl_86, fl_87, fl_88, fl_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -pl_85[k]
                  + f_0 * fl_85[k];

        t_86[k] = -pl_86[k]
                  + f_0 * fl_86[k];

        t_87[k] = -pl_87[k]
                  + f_0 * fl_87[k];

        t_88[k] = -pl_88[k]
                  + f_0 * fl_88[k];

        t_89[k] = -pl_89[k]
                  + f_0 * fl_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pl_90, pl_91, pl_92, pl_93, pl_94, \
                         fl_90, fl_91, fl_92, fl_93, fl_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -pl_90[k]
                  + f_0 * fl_90[k];

        t_91[k] = -pl_91[k]
                  + f_0 * fl_91[k];

        t_92[k] = -pl_92[k]
                  + f_0 * fl_92[k];

        t_93[k] = -pl_93[k]
                  + f_0 * fl_93[k];

        t_94[k] = -pl_94[k]
                  + f_0 * fl_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pl_95, pl_96, pl_97, pl_98, pl_99, \
                         fl_95, fl_96, fl_97, fl_98, fl_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -pl_95[k]
                  + f_0 * fl_95[k];

        t_96[k] = -pl_96[k]
                  + f_0 * fl_96[k];

        t_97[k] = -pl_97[k]
                  + f_0 * fl_97[k];

        t_98[k] = -pl_98[k]
                  + f_0 * fl_98[k];

        t_99[k] = -pl_99[k]
                  + f_0 * fl_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pl_100, pl_101, pl_102, pl_103, \
                         pl_104, fl_100, fl_101, fl_102, fl_103, \
                         fl_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -pl_100[k]
                   + f_0 * fl_100[k];

        t_101[k] = -pl_101[k]
                   + f_0 * fl_101[k];

        t_102[k] = -pl_102[k]
                   + f_0 * fl_102[k];

        t_103[k] = -pl_103[k]
                   + f_0 * fl_103[k];

        t_104[k] = -pl_104[k]
                   + f_0 * fl_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pl_105, pl_106, pl_107, pl_108, \
                         pl_109, fl_105, fl_106, fl_107, fl_108, \
                         fl_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -pl_105[k]
                   + f_0 * fl_105[k];

        t_106[k] = -pl_106[k]
                   + f_0 * fl_106[k];

        t_107[k] = -pl_107[k]
                   + f_0 * fl_107[k];

        t_108[k] = -pl_108[k]
                   + f_0 * fl_108[k];

        t_109[k] = -pl_109[k]
                   + f_0 * fl_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pl_110, pl_111, pl_112, pl_113, \
                         pl_114, fl_110, fl_111, fl_112, fl_113, \
                         fl_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -pl_110[k]
                   + f_0 * fl_110[k];

        t_111[k] = -pl_111[k]
                   + f_0 * fl_111[k];

        t_112[k] = -pl_112[k]
                   + f_0 * fl_112[k];

        t_113[k] = -pl_113[k]
                   + f_0 * fl_113[k];

        t_114[k] = -pl_114[k]
                   + f_0 * fl_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pl_115, pl_116, pl_117, pl_118, \
                         pl_119, fl_115, fl_116, fl_117, fl_118, \
                         fl_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -pl_115[k]
                   + f_0 * fl_115[k];

        t_116[k] = -pl_116[k]
                   + f_0 * fl_116[k];

        t_117[k] = -pl_117[k]
                   + f_0 * fl_117[k];

        t_118[k] = -pl_118[k]
                   + f_0 * fl_118[k];

        t_119[k] = -pl_119[k]
                   + f_0 * fl_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pl_120, pl_121, pl_122, pl_123, \
                         pl_124, fl_120, fl_121, fl_122, fl_123, \
                         fl_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -pl_120[k]
                   + f_0 * fl_120[k];

        t_121[k] = -pl_121[k]
                   + f_0 * fl_121[k];

        t_122[k] = -pl_122[k]
                   + f_0 * fl_122[k];

        t_123[k] = -pl_123[k]
                   + f_0 * fl_123[k];

        t_124[k] = -pl_124[k]
                   + f_0 * fl_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pl_125, pl_126, pl_127, pl_128, \
                         pl_129, fl_125, fl_126, fl_127, fl_128, \
                         fl_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -pl_125[k]
                   + f_0 * fl_125[k];

        t_126[k] = -pl_126[k]
                   + f_0 * fl_126[k];

        t_127[k] = -pl_127[k]
                   + f_0 * fl_127[k];

        t_128[k] = -pl_128[k]
                   + f_0 * fl_128[k];

        t_129[k] = -pl_129[k]
                   + f_0 * fl_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, pl_130, pl_131, pl_132, pl_133, \
                         pl_134, fl_130, fl_131, fl_132, fl_133, \
                         fl_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -pl_130[k]
                   + f_0 * fl_130[k];

        t_131[k] = -pl_131[k]
                   + f_0 * fl_131[k];

        t_132[k] = -pl_132[k]
                   + f_0 * fl_132[k];

        t_133[k] = -pl_133[k]
                   + f_0 * fl_133[k];

        t_134[k] = -pl_134[k]
                   + f_0 * fl_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, t_140, t_141, t_142, fl_135, \
                         fl_136, fl_137, fl_138, fl_139, fl_140, fl_141, \
                         fl_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * fl_135[k];

        t_136[k] = f_0 * fl_136[k];

        t_137[k] = f_0 * fl_137[k];

        t_138[k] = f_0 * fl_138[k];

        t_139[k] = f_0 * fl_139[k];

        t_140[k] = f_0 * fl_140[k];

        t_141[k] = f_0 * fl_141[k];

        t_142[k] = f_0 * fl_142[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, t_148, t_149, t_150, fl_143, \
                         fl_144, fl_145, fl_146, fl_147, fl_148, fl_149, \
                         fl_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_0 * fl_143[k];

        t_144[k] = f_0 * fl_144[k];

        t_145[k] = f_0 * fl_145[k];

        t_146[k] = f_0 * fl_146[k];

        t_147[k] = f_0 * fl_147[k];

        t_148[k] = f_0 * fl_148[k];

        t_149[k] = f_0 * fl_149[k];

        t_150[k] = f_0 * fl_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, t_156, t_157, t_158, fl_151, \
                         fl_152, fl_153, fl_154, fl_155, fl_156, fl_157, \
                         fl_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_0 * fl_151[k];

        t_152[k] = f_0 * fl_152[k];

        t_153[k] = f_0 * fl_153[k];

        t_154[k] = f_0 * fl_154[k];

        t_155[k] = f_0 * fl_155[k];

        t_156[k] = f_0 * fl_156[k];

        t_157[k] = f_0 * fl_157[k];

        t_158[k] = f_0 * fl_158[k];
    }
}

static auto
compute_prim_geom_10_dl_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fl, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, t_164, t_165, t_166, fl_159, \
                         fl_160, fl_161, fl_162, fl_163, fl_164, fl_165, \
                         fl_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_0 * fl_159[k];

        t_160[k] = f_0 * fl_160[k];

        t_161[k] = f_0 * fl_161[k];

        t_162[k] = f_0 * fl_162[k];

        t_163[k] = f_0 * fl_163[k];

        t_164[k] = f_0 * fl_164[k];

        t_165[k] = f_0 * fl_165[k];

        t_166[k] = f_0 * fl_166[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, t_172, t_173, t_174, fl_167, \
                         fl_168, fl_169, fl_170, fl_171, fl_172, fl_173, \
                         fl_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_0 * fl_167[k];

        t_168[k] = f_0 * fl_168[k];

        t_169[k] = f_0 * fl_169[k];

        t_170[k] = f_0 * fl_170[k];

        t_171[k] = f_0 * fl_171[k];

        t_172[k] = f_0 * fl_172[k];

        t_173[k] = f_0 * fl_173[k];

        t_174[k] = f_0 * fl_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, t_180, t_181, t_182, fl_175, \
                         fl_176, fl_177, fl_178, fl_179, fl_180, fl_181, \
                         fl_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_0 * fl_175[k];

        t_176[k] = f_0 * fl_176[k];

        t_177[k] = f_0 * fl_177[k];

        t_178[k] = f_0 * fl_178[k];

        t_179[k] = f_0 * fl_179[k];

        t_180[k] = f_0 * fl_180[k];

        t_181[k] = f_0 * fl_181[k];

        t_182[k] = f_0 * fl_182[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, t_188, t_189, t_190, fl_183, \
                         fl_184, fl_185, fl_186, fl_187, fl_188, fl_189, \
                         fl_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_0 * fl_183[k];

        t_184[k] = f_0 * fl_184[k];

        t_185[k] = f_0 * fl_185[k];

        t_186[k] = f_0 * fl_186[k];

        t_187[k] = f_0 * fl_187[k];

        t_188[k] = f_0 * fl_188[k];

        t_189[k] = f_0 * fl_189[k];

        t_190[k] = f_0 * fl_190[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, t_197, t_198, fl_191, \
                         fl_192, fl_193, fl_194, fl_195, fl_196, fl_197, \
                         fl_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_0 * fl_191[k];

        t_192[k] = f_0 * fl_192[k];

        t_193[k] = f_0 * fl_193[k];

        t_194[k] = f_0 * fl_194[k];

        t_195[k] = f_0 * fl_195[k];

        t_196[k] = f_0 * fl_196[k];

        t_197[k] = f_0 * fl_197[k];

        t_198[k] = f_0 * fl_198[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, t_204, t_205, t_206, fl_199, \
                         fl_200, fl_201, fl_202, fl_203, fl_204, fl_205, \
                         fl_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_0 * fl_199[k];

        t_200[k] = f_0 * fl_200[k];

        t_201[k] = f_0 * fl_201[k];

        t_202[k] = f_0 * fl_202[k];

        t_203[k] = f_0 * fl_203[k];

        t_204[k] = f_0 * fl_204[k];

        t_205[k] = f_0 * fl_205[k];

        t_206[k] = f_0 * fl_206[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, t_213, t_214, fl_207, \
                         fl_208, fl_209, fl_210, fl_211, fl_212, fl_213, \
                         fl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_0 * fl_207[k];

        t_208[k] = f_0 * fl_208[k];

        t_209[k] = f_0 * fl_209[k];

        t_210[k] = f_0 * fl_210[k];

        t_211[k] = f_0 * fl_211[k];

        t_212[k] = f_0 * fl_212[k];

        t_213[k] = f_0 * fl_213[k];

        t_214[k] = f_0 * fl_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, t_220, t_221, t_222, fl_215, \
                         fl_216, fl_217, fl_218, fl_219, fl_220, fl_221, \
                         fl_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_0 * fl_215[k];

        t_216[k] = f_0 * fl_216[k];

        t_217[k] = f_0 * fl_217[k];

        t_218[k] = f_0 * fl_218[k];

        t_219[k] = f_0 * fl_219[k];

        t_220[k] = f_0 * fl_220[k];

        t_221[k] = f_0 * fl_221[k];

        t_222[k] = f_0 * fl_222[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, t_228, t_229, t_230, fl_223, \
                         fl_224, fl_225, fl_226, fl_227, fl_228, fl_229, \
                         fl_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_0 * fl_223[k];

        t_224[k] = f_0 * fl_224[k];

        t_225[k] = f_0 * fl_225[k];

        t_226[k] = f_0 * fl_226[k];

        t_227[k] = f_0 * fl_227[k];

        t_228[k] = f_0 * fl_228[k];

        t_229[k] = f_0 * fl_229[k];

        t_230[k] = f_0 * fl_230[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, t_236, t_237, t_238, fl_231, \
                         fl_232, fl_233, fl_234, fl_235, fl_236, fl_237, \
                         fl_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_0 * fl_231[k];

        t_232[k] = f_0 * fl_232[k];

        t_233[k] = f_0 * fl_233[k];

        t_234[k] = f_0 * fl_234[k];

        t_235[k] = f_0 * fl_235[k];

        t_236[k] = f_0 * fl_236[k];

        t_237[k] = f_0 * fl_237[k];

        t_238[k] = f_0 * fl_238[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, t_244, t_245, t_246, fl_239, \
                         fl_240, fl_241, fl_242, fl_243, fl_244, fl_245, \
                         fl_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_0 * fl_239[k];

        t_240[k] = f_0 * fl_240[k];

        t_241[k] = f_0 * fl_241[k];

        t_242[k] = f_0 * fl_242[k];

        t_243[k] = f_0 * fl_243[k];

        t_244[k] = f_0 * fl_244[k];

        t_245[k] = f_0 * fl_245[k];

        t_246[k] = f_0 * fl_246[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, t_252, t_253, t_254, fl_247, \
                         fl_248, fl_249, fl_250, fl_251, fl_252, fl_253, \
                         fl_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_0 * fl_247[k];

        t_248[k] = f_0 * fl_248[k];

        t_249[k] = f_0 * fl_249[k];

        t_250[k] = f_0 * fl_250[k];

        t_251[k] = f_0 * fl_251[k];

        t_252[k] = f_0 * fl_252[k];

        t_253[k] = f_0 * fl_253[k];

        t_254[k] = f_0 * fl_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, t_260, t_261, t_262, fl_255, \
                         fl_256, fl_257, fl_258, fl_259, fl_260, fl_261, \
                         fl_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_0 * fl_255[k];

        t_256[k] = f_0 * fl_256[k];

        t_257[k] = f_0 * fl_257[k];

        t_258[k] = f_0 * fl_258[k];

        t_259[k] = f_0 * fl_259[k];

        t_260[k] = f_0 * fl_260[k];

        t_261[k] = f_0 * fl_261[k];

        t_262[k] = f_0 * fl_262[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, t_268, t_269, fl_263, fl_264, \
                         fl_265, fl_266, fl_267, fl_268, fl_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_0 * fl_263[k];

        t_264[k] = f_0 * fl_264[k];

        t_265[k] = f_0 * fl_265[k];

        t_266[k] = f_0 * fl_266[k];

        t_267[k] = f_0 * fl_267[k];

        t_268[k] = f_0 * fl_268[k];

        t_269[k] = f_0 * fl_269[k];
    }
}

auto
compute_prim_geom_10_dl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t pl, const size_t fl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_dl_electron_repulsion_0_piece0(buffer, target, pl, fl, ncols, alpha);

    compute_prim_geom_10_dl_electron_repulsion_0_piece1(buffer, target, fl, ncols, alpha);
}

static auto
compute_prim_geom_10_dl_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t pl, const size_t fl,
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

    const auto *pl_0 = buffer.data(pl + 0);
    const auto *pl_1 = buffer.data(pl + 1);
    const auto *pl_2 = buffer.data(pl + 2);
    const auto *pl_3 = buffer.data(pl + 3);
    const auto *pl_4 = buffer.data(pl + 4);
    const auto *pl_5 = buffer.data(pl + 5);
    const auto *pl_6 = buffer.data(pl + 6);
    const auto *pl_7 = buffer.data(pl + 7);
    const auto *pl_8 = buffer.data(pl + 8);
    const auto *pl_9 = buffer.data(pl + 9);
    const auto *pl_10 = buffer.data(pl + 10);
    const auto *pl_11 = buffer.data(pl + 11);
    const auto *pl_12 = buffer.data(pl + 12);
    const auto *pl_13 = buffer.data(pl + 13);
    const auto *pl_14 = buffer.data(pl + 14);
    const auto *pl_15 = buffer.data(pl + 15);
    const auto *pl_16 = buffer.data(pl + 16);
    const auto *pl_17 = buffer.data(pl + 17);
    const auto *pl_18 = buffer.data(pl + 18);
    const auto *pl_19 = buffer.data(pl + 19);
    const auto *pl_20 = buffer.data(pl + 20);
    const auto *pl_21 = buffer.data(pl + 21);
    const auto *pl_22 = buffer.data(pl + 22);
    const auto *pl_23 = buffer.data(pl + 23);
    const auto *pl_24 = buffer.data(pl + 24);
    const auto *pl_25 = buffer.data(pl + 25);
    const auto *pl_26 = buffer.data(pl + 26);
    const auto *pl_27 = buffer.data(pl + 27);
    const auto *pl_28 = buffer.data(pl + 28);
    const auto *pl_29 = buffer.data(pl + 29);
    const auto *pl_30 = buffer.data(pl + 30);
    const auto *pl_31 = buffer.data(pl + 31);
    const auto *pl_32 = buffer.data(pl + 32);
    const auto *pl_33 = buffer.data(pl + 33);
    const auto *pl_34 = buffer.data(pl + 34);
    const auto *pl_35 = buffer.data(pl + 35);
    const auto *pl_36 = buffer.data(pl + 36);
    const auto *pl_37 = buffer.data(pl + 37);
    const auto *pl_38 = buffer.data(pl + 38);
    const auto *pl_39 = buffer.data(pl + 39);
    const auto *pl_40 = buffer.data(pl + 40);
    const auto *pl_41 = buffer.data(pl + 41);
    const auto *pl_42 = buffer.data(pl + 42);
    const auto *pl_43 = buffer.data(pl + 43);
    const auto *pl_44 = buffer.data(pl + 44);
    const auto *pl_45 = buffer.data(pl + 45);
    const auto *pl_46 = buffer.data(pl + 46);
    const auto *pl_47 = buffer.data(pl + 47);
    const auto *pl_48 = buffer.data(pl + 48);
    const auto *pl_49 = buffer.data(pl + 49);
    const auto *pl_50 = buffer.data(pl + 50);
    const auto *pl_51 = buffer.data(pl + 51);
    const auto *pl_52 = buffer.data(pl + 52);
    const auto *pl_53 = buffer.data(pl + 53);
    const auto *pl_54 = buffer.data(pl + 54);
    const auto *pl_55 = buffer.data(pl + 55);
    const auto *pl_56 = buffer.data(pl + 56);
    const auto *pl_57 = buffer.data(pl + 57);
    const auto *pl_58 = buffer.data(pl + 58);
    const auto *pl_59 = buffer.data(pl + 59);
    const auto *pl_60 = buffer.data(pl + 60);
    const auto *pl_61 = buffer.data(pl + 61);
    const auto *pl_62 = buffer.data(pl + 62);
    const auto *pl_63 = buffer.data(pl + 63);
    const auto *pl_64 = buffer.data(pl + 64);
    const auto *pl_65 = buffer.data(pl + 65);
    const auto *pl_66 = buffer.data(pl + 66);
    const auto *pl_67 = buffer.data(pl + 67);
    const auto *pl_68 = buffer.data(pl + 68);
    const auto *pl_69 = buffer.data(pl + 69);
    const auto *pl_70 = buffer.data(pl + 70);
    const auto *pl_71 = buffer.data(pl + 71);
    const auto *pl_72 = buffer.data(pl + 72);
    const auto *pl_73 = buffer.data(pl + 73);
    const auto *pl_74 = buffer.data(pl + 74);
    const auto *pl_75 = buffer.data(pl + 75);
    const auto *pl_76 = buffer.data(pl + 76);
    const auto *pl_77 = buffer.data(pl + 77);
    const auto *pl_78 = buffer.data(pl + 78);
    const auto *pl_79 = buffer.data(pl + 79);
    const auto *pl_80 = buffer.data(pl + 80);
    const auto *pl_81 = buffer.data(pl + 81);
    const auto *pl_82 = buffer.data(pl + 82);
    const auto *pl_83 = buffer.data(pl + 83);
    const auto *pl_84 = buffer.data(pl + 84);
    const auto *pl_85 = buffer.data(pl + 85);
    const auto *pl_86 = buffer.data(pl + 86);
    const auto *pl_87 = buffer.data(pl + 87);
    const auto *pl_88 = buffer.data(pl + 88);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, fl_45, fl_46, fl_47, fl_48, \
                         fl_49, fl_50, fl_51, fl_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fl_45[k];

        t_1[k] = f_0 * fl_46[k];

        t_2[k] = f_0 * fl_47[k];

        t_3[k] = f_0 * fl_48[k];

        t_4[k] = f_0 * fl_49[k];

        t_5[k] = f_0 * fl_50[k];

        t_6[k] = f_0 * fl_51[k];

        t_7[k] = f_0 * fl_52[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, fl_53, fl_54, fl_55, \
                         fl_56, fl_57, fl_58, fl_59, fl_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * fl_53[k];

        t_9[k] = f_0 * fl_54[k];

        t_10[k] = f_0 * fl_55[k];

        t_11[k] = f_0 * fl_56[k];

        t_12[k] = f_0 * fl_57[k];

        t_13[k] = f_0 * fl_58[k];

        t_14[k] = f_0 * fl_59[k];

        t_15[k] = f_0 * fl_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, fl_61, fl_62, fl_63, \
                         fl_64, fl_65, fl_66, fl_67, fl_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * fl_61[k];

        t_17[k] = f_0 * fl_62[k];

        t_18[k] = f_0 * fl_63[k];

        t_19[k] = f_0 * fl_64[k];

        t_20[k] = f_0 * fl_65[k];

        t_21[k] = f_0 * fl_66[k];

        t_22[k] = f_0 * fl_67[k];

        t_23[k] = f_0 * fl_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, fl_69, fl_70, fl_71, \
                         fl_72, fl_73, fl_74, fl_75, fl_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * fl_69[k];

        t_25[k] = f_0 * fl_70[k];

        t_26[k] = f_0 * fl_71[k];

        t_27[k] = f_0 * fl_72[k];

        t_28[k] = f_0 * fl_73[k];

        t_29[k] = f_0 * fl_74[k];

        t_30[k] = f_0 * fl_75[k];

        t_31[k] = f_0 * fl_76[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, fl_77, fl_78, fl_79, \
                         fl_80, fl_81, fl_82, fl_83, fl_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * fl_77[k];

        t_33[k] = f_0 * fl_78[k];

        t_34[k] = f_0 * fl_79[k];

        t_35[k] = f_0 * fl_80[k];

        t_36[k] = f_0 * fl_81[k];

        t_37[k] = f_0 * fl_82[k];

        t_38[k] = f_0 * fl_83[k];

        t_39[k] = f_0 * fl_84[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, pl_0, pl_1, fl_85, fl_86, \
                         fl_87, fl_88, fl_89, fl_135, fl_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * fl_85[k];

        t_41[k] = f_0 * fl_86[k];

        t_42[k] = f_0 * fl_87[k];

        t_43[k] = f_0 * fl_88[k];

        t_44[k] = f_0 * fl_89[k];

        t_45[k] = -pl_0[k]
                  + f_0 * fl_135[k];

        t_46[k] = -pl_1[k]
                  + f_0 * fl_136[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pl_2, pl_3, pl_4, pl_5, pl_6, fl_137, \
                         fl_138, fl_139, fl_140, fl_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -pl_2[k]
                  + f_0 * fl_137[k];

        t_48[k] = -pl_3[k]
                  + f_0 * fl_138[k];

        t_49[k] = -pl_4[k]
                  + f_0 * fl_139[k];

        t_50[k] = -pl_5[k]
                  + f_0 * fl_140[k];

        t_51[k] = -pl_6[k]
                  + f_0 * fl_141[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pl_7, pl_8, pl_9, pl_10, pl_11, fl_142, \
                         fl_143, fl_144, fl_145, fl_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = -pl_7[k]
                  + f_0 * fl_142[k];

        t_53[k] = -pl_8[k]
                  + f_0 * fl_143[k];

        t_54[k] = -pl_9[k]
                  + f_0 * fl_144[k];

        t_55[k] = -pl_10[k]
                  + f_0 * fl_145[k];

        t_56[k] = -pl_11[k]
                  + f_0 * fl_146[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pl_12, pl_13, pl_14, pl_15, pl_16, \
                         fl_147, fl_148, fl_149, fl_150, fl_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -pl_12[k]
                  + f_0 * fl_147[k];

        t_58[k] = -pl_13[k]
                  + f_0 * fl_148[k];

        t_59[k] = -pl_14[k]
                  + f_0 * fl_149[k];

        t_60[k] = -pl_15[k]
                  + f_0 * fl_150[k];

        t_61[k] = -pl_16[k]
                  + f_0 * fl_151[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, pl_17, pl_18, pl_19, pl_20, pl_21, \
                         fl_152, fl_153, fl_154, fl_155, fl_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -pl_17[k]
                  + f_0 * fl_152[k];

        t_63[k] = -pl_18[k]
                  + f_0 * fl_153[k];

        t_64[k] = -pl_19[k]
                  + f_0 * fl_154[k];

        t_65[k] = -pl_20[k]
                  + f_0 * fl_155[k];

        t_66[k] = -pl_21[k]
                  + f_0 * fl_156[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pl_22, pl_23, pl_24, pl_25, pl_26, \
                         fl_157, fl_158, fl_159, fl_160, fl_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -pl_22[k]
                  + f_0 * fl_157[k];

        t_68[k] = -pl_23[k]
                  + f_0 * fl_158[k];

        t_69[k] = -pl_24[k]
                  + f_0 * fl_159[k];

        t_70[k] = -pl_25[k]
                  + f_0 * fl_160[k];

        t_71[k] = -pl_26[k]
                  + f_0 * fl_161[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pl_27, pl_28, pl_29, pl_30, pl_31, \
                         fl_162, fl_163, fl_164, fl_165, fl_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -pl_27[k]
                  + f_0 * fl_162[k];

        t_73[k] = -pl_28[k]
                  + f_0 * fl_163[k];

        t_74[k] = -pl_29[k]
                  + f_0 * fl_164[k];

        t_75[k] = -pl_30[k]
                  + f_0 * fl_165[k];

        t_76[k] = -pl_31[k]
                  + f_0 * fl_166[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, pl_32, pl_33, pl_34, pl_35, pl_36, \
                         fl_167, fl_168, fl_169, fl_170, fl_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -pl_32[k]
                  + f_0 * fl_167[k];

        t_78[k] = -pl_33[k]
                  + f_0 * fl_168[k];

        t_79[k] = -pl_34[k]
                  + f_0 * fl_169[k];

        t_80[k] = -pl_35[k]
                  + f_0 * fl_170[k];

        t_81[k] = -pl_36[k]
                  + f_0 * fl_171[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pl_37, pl_38, pl_39, pl_40, pl_41, \
                         fl_172, fl_173, fl_174, fl_175, fl_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -pl_37[k]
                  + f_0 * fl_172[k];

        t_83[k] = -pl_38[k]
                  + f_0 * fl_173[k];

        t_84[k] = -pl_39[k]
                  + f_0 * fl_174[k];

        t_85[k] = -pl_40[k]
                  + f_0 * fl_175[k];

        t_86[k] = -pl_41[k]
                  + f_0 * fl_176[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, pl_42, pl_43, pl_44, fl_177, \
                         fl_178, fl_179, fl_180, fl_181, fl_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -pl_42[k]
                  + f_0 * fl_177[k];

        t_88[k] = -pl_43[k]
                  + f_0 * fl_178[k];

        t_89[k] = -pl_44[k]
                  + f_0 * fl_179[k];

        t_90[k] = f_0 * fl_180[k];

        t_91[k] = f_0 * fl_181[k];

        t_92[k] = f_0 * fl_182[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, fl_183, fl_184, \
                         fl_185, fl_186, fl_187, fl_188, fl_189, \
                         fl_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * fl_183[k];

        t_94[k] = f_0 * fl_184[k];

        t_95[k] = f_0 * fl_185[k];

        t_96[k] = f_0 * fl_186[k];

        t_97[k] = f_0 * fl_187[k];

        t_98[k] = f_0 * fl_188[k];

        t_99[k] = f_0 * fl_189[k];

        t_100[k] = f_0 * fl_190[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, t_107, t_108, fl_191, \
                         fl_192, fl_193, fl_194, fl_195, fl_196, fl_197, \
                         fl_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * fl_191[k];

        t_102[k] = f_0 * fl_192[k];

        t_103[k] = f_0 * fl_193[k];

        t_104[k] = f_0 * fl_194[k];

        t_105[k] = f_0 * fl_195[k];

        t_106[k] = f_0 * fl_196[k];

        t_107[k] = f_0 * fl_197[k];

        t_108[k] = f_0 * fl_198[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, t_114, t_115, t_116, fl_199, \
                         fl_200, fl_201, fl_202, fl_203, fl_204, fl_205, \
                         fl_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_0 * fl_199[k];

        t_110[k] = f_0 * fl_200[k];

        t_111[k] = f_0 * fl_201[k];

        t_112[k] = f_0 * fl_202[k];

        t_113[k] = f_0 * fl_203[k];

        t_114[k] = f_0 * fl_204[k];

        t_115[k] = f_0 * fl_205[k];

        t_116[k] = f_0 * fl_206[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, t_122, t_123, t_124, fl_207, \
                         fl_208, fl_209, fl_210, fl_211, fl_212, fl_213, \
                         fl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_0 * fl_207[k];

        t_118[k] = f_0 * fl_208[k];

        t_119[k] = f_0 * fl_209[k];

        t_120[k] = f_0 * fl_210[k];

        t_121[k] = f_0 * fl_211[k];

        t_122[k] = f_0 * fl_212[k];

        t_123[k] = f_0 * fl_213[k];

        t_124[k] = f_0 * fl_214[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, t_131, t_132, fl_215, \
                         fl_216, fl_217, fl_218, fl_219, fl_220, fl_221, \
                         fl_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * fl_215[k];

        t_126[k] = f_0 * fl_216[k];

        t_127[k] = f_0 * fl_217[k];

        t_128[k] = f_0 * fl_218[k];

        t_129[k] = f_0 * fl_219[k];

        t_130[k] = f_0 * fl_220[k];

        t_131[k] = f_0 * fl_221[k];

        t_132[k] = f_0 * fl_222[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, t_138, pl_45, pl_46, pl_47, pl_48, \
                         fl_223, fl_224, fl_270, fl_271, fl_272, \
                         fl_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_0 * fl_223[k];

        t_134[k] = f_0 * fl_224[k];

        t_135[k] = -2.0 * pl_45[k]
                   + f_0 * fl_270[k];

        t_136[k] = -2.0 * pl_46[k]
                   + f_0 * fl_271[k];

        t_137[k] = -2.0 * pl_47[k]
                   + f_0 * fl_272[k];

        t_138[k] = -2.0 * pl_48[k]
                   + f_0 * fl_273[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pl_49, pl_50, pl_51, pl_52, pl_53, \
                         fl_274, fl_275, fl_276, fl_277, fl_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -2.0 * pl_49[k]
                   + f_0 * fl_274[k];

        t_140[k] = -2.0 * pl_50[k]
                   + f_0 * fl_275[k];

        t_141[k] = -2.0 * pl_51[k]
                   + f_0 * fl_276[k];

        t_142[k] = -2.0 * pl_52[k]
                   + f_0 * fl_277[k];

        t_143[k] = -2.0 * pl_53[k]
                   + f_0 * fl_278[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pl_54, pl_55, pl_56, pl_57, pl_58, \
                         fl_279, fl_280, fl_281, fl_282, fl_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -2.0 * pl_54[k]
                   + f_0 * fl_279[k];

        t_145[k] = -2.0 * pl_55[k]
                   + f_0 * fl_280[k];

        t_146[k] = -2.0 * pl_56[k]
                   + f_0 * fl_281[k];

        t_147[k] = -2.0 * pl_57[k]
                   + f_0 * fl_282[k];

        t_148[k] = -2.0 * pl_58[k]
                   + f_0 * fl_283[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, pl_59, pl_60, pl_61, pl_62, pl_63, \
                         fl_284, fl_285, fl_286, fl_287, fl_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -2.0 * pl_59[k]
                   + f_0 * fl_284[k];

        t_150[k] = -2.0 * pl_60[k]
                   + f_0 * fl_285[k];

        t_151[k] = -2.0 * pl_61[k]
                   + f_0 * fl_286[k];

        t_152[k] = -2.0 * pl_62[k]
                   + f_0 * fl_287[k];

        t_153[k] = -2.0 * pl_63[k]
                   + f_0 * fl_288[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, pl_64, pl_65, pl_66, pl_67, pl_68, \
                         fl_289, fl_290, fl_291, fl_292, fl_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -2.0 * pl_64[k]
                   + f_0 * fl_289[k];

        t_155[k] = -2.0 * pl_65[k]
                   + f_0 * fl_290[k];

        t_156[k] = -2.0 * pl_66[k]
                   + f_0 * fl_291[k];

        t_157[k] = -2.0 * pl_67[k]
                   + f_0 * fl_292[k];

        t_158[k] = -2.0 * pl_68[k]
                   + f_0 * fl_293[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, pl_69, pl_70, pl_71, pl_72, pl_73, \
                         fl_294, fl_295, fl_296, fl_297, fl_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -2.0 * pl_69[k]
                   + f_0 * fl_294[k];

        t_160[k] = -2.0 * pl_70[k]
                   + f_0 * fl_295[k];

        t_161[k] = -2.0 * pl_71[k]
                   + f_0 * fl_296[k];

        t_162[k] = -2.0 * pl_72[k]
                   + f_0 * fl_297[k];

        t_163[k] = -2.0 * pl_73[k]
                   + f_0 * fl_298[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, pl_74, pl_75, pl_76, pl_77, pl_78, \
                         fl_299, fl_300, fl_301, fl_302, fl_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -2.0 * pl_74[k]
                   + f_0 * fl_299[k];

        t_165[k] = -2.0 * pl_75[k]
                   + f_0 * fl_300[k];

        t_166[k] = -2.0 * pl_76[k]
                   + f_0 * fl_301[k];

        t_167[k] = -2.0 * pl_77[k]
                   + f_0 * fl_302[k];

        t_168[k] = -2.0 * pl_78[k]
                   + f_0 * fl_303[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, pl_79, pl_80, pl_81, pl_82, pl_83, \
                         fl_304, fl_305, fl_306, fl_307, fl_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -2.0 * pl_79[k]
                   + f_0 * fl_304[k];

        t_170[k] = -2.0 * pl_80[k]
                   + f_0 * fl_305[k];

        t_171[k] = -2.0 * pl_81[k]
                   + f_0 * fl_306[k];

        t_172[k] = -2.0 * pl_82[k]
                   + f_0 * fl_307[k];

        t_173[k] = -2.0 * pl_83[k]
                   + f_0 * fl_308[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, pl_84, pl_85, pl_86, pl_87, pl_88, \
                         fl_309, fl_310, fl_311, fl_312, fl_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -2.0 * pl_84[k]
                   + f_0 * fl_309[k];

        t_175[k] = -2.0 * pl_85[k]
                   + f_0 * fl_310[k];

        t_176[k] = -2.0 * pl_86[k]
                   + f_0 * fl_311[k];

        t_177[k] = -2.0 * pl_87[k]
                   + f_0 * fl_312[k];

        t_178[k] = -2.0 * pl_88[k]
                   + f_0 * fl_313[k];
    }
}

static auto
compute_prim_geom_10_dl_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t pl, const size_t fl,
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

    const auto *pl_89 = buffer.data(pl + 89);
    const auto *pl_90 = buffer.data(pl + 90);
    const auto *pl_91 = buffer.data(pl + 91);
    const auto *pl_92 = buffer.data(pl + 92);
    const auto *pl_93 = buffer.data(pl + 93);
    const auto *pl_94 = buffer.data(pl + 94);
    const auto *pl_95 = buffer.data(pl + 95);
    const auto *pl_96 = buffer.data(pl + 96);
    const auto *pl_97 = buffer.data(pl + 97);
    const auto *pl_98 = buffer.data(pl + 98);
    const auto *pl_99 = buffer.data(pl + 99);
    const auto *pl_100 = buffer.data(pl + 100);
    const auto *pl_101 = buffer.data(pl + 101);
    const auto *pl_102 = buffer.data(pl + 102);
    const auto *pl_103 = buffer.data(pl + 103);
    const auto *pl_104 = buffer.data(pl + 104);
    const auto *pl_105 = buffer.data(pl + 105);
    const auto *pl_106 = buffer.data(pl + 106);
    const auto *pl_107 = buffer.data(pl + 107);
    const auto *pl_108 = buffer.data(pl + 108);
    const auto *pl_109 = buffer.data(pl + 109);
    const auto *pl_110 = buffer.data(pl + 110);
    const auto *pl_111 = buffer.data(pl + 111);
    const auto *pl_112 = buffer.data(pl + 112);
    const auto *pl_113 = buffer.data(pl + 113);
    const auto *pl_114 = buffer.data(pl + 114);
    const auto *pl_115 = buffer.data(pl + 115);
    const auto *pl_116 = buffer.data(pl + 116);
    const auto *pl_117 = buffer.data(pl + 117);
    const auto *pl_118 = buffer.data(pl + 118);
    const auto *pl_119 = buffer.data(pl + 119);
    const auto *pl_120 = buffer.data(pl + 120);
    const auto *pl_121 = buffer.data(pl + 121);
    const auto *pl_122 = buffer.data(pl + 122);
    const auto *pl_123 = buffer.data(pl + 123);
    const auto *pl_124 = buffer.data(pl + 124);
    const auto *pl_125 = buffer.data(pl + 125);
    const auto *pl_126 = buffer.data(pl + 126);
    const auto *pl_127 = buffer.data(pl + 127);
    const auto *pl_128 = buffer.data(pl + 128);
    const auto *pl_129 = buffer.data(pl + 129);
    const auto *pl_130 = buffer.data(pl + 130);
    const auto *pl_131 = buffer.data(pl + 131);
    const auto *pl_132 = buffer.data(pl + 132);
    const auto *pl_133 = buffer.data(pl + 133);
    const auto *pl_134 = buffer.data(pl + 134);

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

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, pl_89, pl_90, pl_91, pl_92, pl_93, \
                         fl_314, fl_315, fl_316, fl_317, fl_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -2.0 * pl_89[k]
                   + f_0 * fl_314[k];

        t_180[k] = -pl_90[k]
                   + f_0 * fl_315[k];

        t_181[k] = -pl_91[k]
                   + f_0 * fl_316[k];

        t_182[k] = -pl_92[k]
                   + f_0 * fl_317[k];

        t_183[k] = -pl_93[k]
                   + f_0 * fl_318[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pl_94, pl_95, pl_96, pl_97, pl_98, \
                         fl_319, fl_320, fl_321, fl_322, fl_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -pl_94[k]
                   + f_0 * fl_319[k];

        t_185[k] = -pl_95[k]
                   + f_0 * fl_320[k];

        t_186[k] = -pl_96[k]
                   + f_0 * fl_321[k];

        t_187[k] = -pl_97[k]
                   + f_0 * fl_322[k];

        t_188[k] = -pl_98[k]
                   + f_0 * fl_323[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pl_99, pl_100, pl_101, pl_102, \
                         pl_103, fl_324, fl_325, fl_326, fl_327, \
                         fl_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = -pl_99[k]
                   + f_0 * fl_324[k];

        t_190[k] = -pl_100[k]
                   + f_0 * fl_325[k];

        t_191[k] = -pl_101[k]
                   + f_0 * fl_326[k];

        t_192[k] = -pl_102[k]
                   + f_0 * fl_327[k];

        t_193[k] = -pl_103[k]
                   + f_0 * fl_328[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, pl_104, pl_105, pl_106, pl_107, \
                         pl_108, fl_329, fl_330, fl_331, fl_332, \
                         fl_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = -pl_104[k]
                   + f_0 * fl_329[k];

        t_195[k] = -pl_105[k]
                   + f_0 * fl_330[k];

        t_196[k] = -pl_106[k]
                   + f_0 * fl_331[k];

        t_197[k] = -pl_107[k]
                   + f_0 * fl_332[k];

        t_198[k] = -pl_108[k]
                   + f_0 * fl_333[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, pl_109, pl_110, pl_111, pl_112, \
                         pl_113, fl_334, fl_335, fl_336, fl_337, \
                         fl_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = -pl_109[k]
                   + f_0 * fl_334[k];

        t_200[k] = -pl_110[k]
                   + f_0 * fl_335[k];

        t_201[k] = -pl_111[k]
                   + f_0 * fl_336[k];

        t_202[k] = -pl_112[k]
                   + f_0 * fl_337[k];

        t_203[k] = -pl_113[k]
                   + f_0 * fl_338[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, pl_114, pl_115, pl_116, pl_117, \
                         pl_118, fl_339, fl_340, fl_341, fl_342, \
                         fl_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = -pl_114[k]
                   + f_0 * fl_339[k];

        t_205[k] = -pl_115[k]
                   + f_0 * fl_340[k];

        t_206[k] = -pl_116[k]
                   + f_0 * fl_341[k];

        t_207[k] = -pl_117[k]
                   + f_0 * fl_342[k];

        t_208[k] = -pl_118[k]
                   + f_0 * fl_343[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, pl_119, pl_120, pl_121, pl_122, \
                         pl_123, fl_344, fl_345, fl_346, fl_347, \
                         fl_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = -pl_119[k]
                   + f_0 * fl_344[k];

        t_210[k] = -pl_120[k]
                   + f_0 * fl_345[k];

        t_211[k] = -pl_121[k]
                   + f_0 * fl_346[k];

        t_212[k] = -pl_122[k]
                   + f_0 * fl_347[k];

        t_213[k] = -pl_123[k]
                   + f_0 * fl_348[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pl_124, pl_125, pl_126, pl_127, \
                         pl_128, fl_349, fl_350, fl_351, fl_352, \
                         fl_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = -pl_124[k]
                   + f_0 * fl_349[k];

        t_215[k] = -pl_125[k]
                   + f_0 * fl_350[k];

        t_216[k] = -pl_126[k]
                   + f_0 * fl_351[k];

        t_217[k] = -pl_127[k]
                   + f_0 * fl_352[k];

        t_218[k] = -pl_128[k]
                   + f_0 * fl_353[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, pl_129, pl_130, pl_131, pl_132, \
                         pl_133, fl_354, fl_355, fl_356, fl_357, \
                         fl_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = -pl_129[k]
                   + f_0 * fl_354[k];

        t_220[k] = -pl_130[k]
                   + f_0 * fl_355[k];

        t_221[k] = -pl_131[k]
                   + f_0 * fl_356[k];

        t_222[k] = -pl_132[k]
                   + f_0 * fl_357[k];

        t_223[k] = -pl_133[k]
                   + f_0 * fl_358[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, t_229, t_230, pl_134, fl_359, \
                         fl_360, fl_361, fl_362, fl_363, fl_364, \
                         fl_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = -pl_134[k]
                   + f_0 * fl_359[k];

        t_225[k] = f_0 * fl_360[k];

        t_226[k] = f_0 * fl_361[k];

        t_227[k] = f_0 * fl_362[k];

        t_228[k] = f_0 * fl_363[k];

        t_229[k] = f_0 * fl_364[k];

        t_230[k] = f_0 * fl_365[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, t_236, t_237, t_238, fl_366, \
                         fl_367, fl_368, fl_369, fl_370, fl_371, fl_372, \
                         fl_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_0 * fl_366[k];

        t_232[k] = f_0 * fl_367[k];

        t_233[k] = f_0 * fl_368[k];

        t_234[k] = f_0 * fl_369[k];

        t_235[k] = f_0 * fl_370[k];

        t_236[k] = f_0 * fl_371[k];

        t_237[k] = f_0 * fl_372[k];

        t_238[k] = f_0 * fl_373[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, t_244, t_245, t_246, fl_374, \
                         fl_375, fl_376, fl_377, fl_378, fl_379, fl_380, \
                         fl_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_0 * fl_374[k];

        t_240[k] = f_0 * fl_375[k];

        t_241[k] = f_0 * fl_376[k];

        t_242[k] = f_0 * fl_377[k];

        t_243[k] = f_0 * fl_378[k];

        t_244[k] = f_0 * fl_379[k];

        t_245[k] = f_0 * fl_380[k];

        t_246[k] = f_0 * fl_381[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, t_252, t_253, t_254, fl_382, \
                         fl_383, fl_384, fl_385, fl_386, fl_387, fl_388, \
                         fl_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_0 * fl_382[k];

        t_248[k] = f_0 * fl_383[k];

        t_249[k] = f_0 * fl_384[k];

        t_250[k] = f_0 * fl_385[k];

        t_251[k] = f_0 * fl_386[k];

        t_252[k] = f_0 * fl_387[k];

        t_253[k] = f_0 * fl_388[k];

        t_254[k] = f_0 * fl_389[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, t_260, t_261, t_262, fl_390, \
                         fl_391, fl_392, fl_393, fl_394, fl_395, fl_396, \
                         fl_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_0 * fl_390[k];

        t_256[k] = f_0 * fl_391[k];

        t_257[k] = f_0 * fl_392[k];

        t_258[k] = f_0 * fl_393[k];

        t_259[k] = f_0 * fl_394[k];

        t_260[k] = f_0 * fl_395[k];

        t_261[k] = f_0 * fl_396[k];

        t_262[k] = f_0 * fl_397[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, t_268, t_269, fl_398, fl_399, \
                         fl_400, fl_401, fl_402, fl_403, fl_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_0 * fl_398[k];

        t_264[k] = f_0 * fl_399[k];

        t_265[k] = f_0 * fl_400[k];

        t_266[k] = f_0 * fl_401[k];

        t_267[k] = f_0 * fl_402[k];

        t_268[k] = f_0 * fl_403[k];

        t_269[k] = f_0 * fl_404[k];
    }
}

auto
compute_prim_geom_10_dl_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t pl, const size_t fl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_dl_electron_repulsion_1_piece0(buffer, target, pl, fl, ncols, alpha);

    compute_prim_geom_10_dl_electron_repulsion_1_piece1(buffer, target, pl, fl, ncols, alpha);
}

static auto
compute_prim_geom_10_dl_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t pl, const size_t fl,
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

    const auto *pl_0 = buffer.data(pl + 0);
    const auto *pl_1 = buffer.data(pl + 1);
    const auto *pl_2 = buffer.data(pl + 2);
    const auto *pl_3 = buffer.data(pl + 3);
    const auto *pl_4 = buffer.data(pl + 4);
    const auto *pl_5 = buffer.data(pl + 5);
    const auto *pl_6 = buffer.data(pl + 6);
    const auto *pl_7 = buffer.data(pl + 7);
    const auto *pl_8 = buffer.data(pl + 8);
    const auto *pl_9 = buffer.data(pl + 9);
    const auto *pl_10 = buffer.data(pl + 10);
    const auto *pl_11 = buffer.data(pl + 11);
    const auto *pl_12 = buffer.data(pl + 12);
    const auto *pl_13 = buffer.data(pl + 13);
    const auto *pl_14 = buffer.data(pl + 14);
    const auto *pl_15 = buffer.data(pl + 15);
    const auto *pl_16 = buffer.data(pl + 16);
    const auto *pl_17 = buffer.data(pl + 17);
    const auto *pl_18 = buffer.data(pl + 18);
    const auto *pl_19 = buffer.data(pl + 19);
    const auto *pl_20 = buffer.data(pl + 20);
    const auto *pl_21 = buffer.data(pl + 21);
    const auto *pl_22 = buffer.data(pl + 22);
    const auto *pl_23 = buffer.data(pl + 23);
    const auto *pl_24 = buffer.data(pl + 24);
    const auto *pl_25 = buffer.data(pl + 25);
    const auto *pl_26 = buffer.data(pl + 26);
    const auto *pl_27 = buffer.data(pl + 27);
    const auto *pl_28 = buffer.data(pl + 28);
    const auto *pl_29 = buffer.data(pl + 29);
    const auto *pl_30 = buffer.data(pl + 30);
    const auto *pl_31 = buffer.data(pl + 31);
    const auto *pl_32 = buffer.data(pl + 32);
    const auto *pl_33 = buffer.data(pl + 33);
    const auto *pl_34 = buffer.data(pl + 34);
    const auto *pl_35 = buffer.data(pl + 35);
    const auto *pl_36 = buffer.data(pl + 36);
    const auto *pl_37 = buffer.data(pl + 37);
    const auto *pl_38 = buffer.data(pl + 38);
    const auto *pl_39 = buffer.data(pl + 39);
    const auto *pl_40 = buffer.data(pl + 40);
    const auto *pl_41 = buffer.data(pl + 41);
    const auto *pl_42 = buffer.data(pl + 42);
    const auto *pl_43 = buffer.data(pl + 43);
    const auto *pl_44 = buffer.data(pl + 44);
    const auto *pl_45 = buffer.data(pl + 45);
    const auto *pl_46 = buffer.data(pl + 46);
    const auto *pl_47 = buffer.data(pl + 47);
    const auto *pl_48 = buffer.data(pl + 48);
    const auto *pl_49 = buffer.data(pl + 49);
    const auto *pl_50 = buffer.data(pl + 50);
    const auto *pl_51 = buffer.data(pl + 51);
    const auto *pl_52 = buffer.data(pl + 52);
    const auto *pl_53 = buffer.data(pl + 53);
    const auto *pl_54 = buffer.data(pl + 54);
    const auto *pl_55 = buffer.data(pl + 55);
    const auto *pl_56 = buffer.data(pl + 56);
    const auto *pl_57 = buffer.data(pl + 57);
    const auto *pl_58 = buffer.data(pl + 58);
    const auto *pl_59 = buffer.data(pl + 59);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, fl_90, fl_91, fl_92, fl_93, \
                         fl_94, fl_95, fl_96, fl_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fl_90[k];

        t_1[k] = f_0 * fl_91[k];

        t_2[k] = f_0 * fl_92[k];

        t_3[k] = f_0 * fl_93[k];

        t_4[k] = f_0 * fl_94[k];

        t_5[k] = f_0 * fl_95[k];

        t_6[k] = f_0 * fl_96[k];

        t_7[k] = f_0 * fl_97[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, fl_98, fl_99, fl_100, \
                         fl_101, fl_102, fl_103, fl_104, fl_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * fl_98[k];

        t_9[k] = f_0 * fl_99[k];

        t_10[k] = f_0 * fl_100[k];

        t_11[k] = f_0 * fl_101[k];

        t_12[k] = f_0 * fl_102[k];

        t_13[k] = f_0 * fl_103[k];

        t_14[k] = f_0 * fl_104[k];

        t_15[k] = f_0 * fl_105[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, fl_106, fl_107, \
                         fl_108, fl_109, fl_110, fl_111, fl_112, \
                         fl_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * fl_106[k];

        t_17[k] = f_0 * fl_107[k];

        t_18[k] = f_0 * fl_108[k];

        t_19[k] = f_0 * fl_109[k];

        t_20[k] = f_0 * fl_110[k];

        t_21[k] = f_0 * fl_111[k];

        t_22[k] = f_0 * fl_112[k];

        t_23[k] = f_0 * fl_113[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, fl_114, fl_115, \
                         fl_116, fl_117, fl_118, fl_119, fl_120, \
                         fl_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * fl_114[k];

        t_25[k] = f_0 * fl_115[k];

        t_26[k] = f_0 * fl_116[k];

        t_27[k] = f_0 * fl_117[k];

        t_28[k] = f_0 * fl_118[k];

        t_29[k] = f_0 * fl_119[k];

        t_30[k] = f_0 * fl_120[k];

        t_31[k] = f_0 * fl_121[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, fl_122, fl_123, \
                         fl_124, fl_125, fl_126, fl_127, fl_128, \
                         fl_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * fl_122[k];

        t_33[k] = f_0 * fl_123[k];

        t_34[k] = f_0 * fl_124[k];

        t_35[k] = f_0 * fl_125[k];

        t_36[k] = f_0 * fl_126[k];

        t_37[k] = f_0 * fl_127[k];

        t_38[k] = f_0 * fl_128[k];

        t_39[k] = f_0 * fl_129[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, fl_130, fl_131, \
                         fl_132, fl_133, fl_134, fl_180, fl_181, \
                         fl_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * fl_130[k];

        t_41[k] = f_0 * fl_131[k];

        t_42[k] = f_0 * fl_132[k];

        t_43[k] = f_0 * fl_133[k];

        t_44[k] = f_0 * fl_134[k];

        t_45[k] = f_0 * fl_180[k];

        t_46[k] = f_0 * fl_181[k];

        t_47[k] = f_0 * fl_182[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, fl_183, fl_184, \
                         fl_185, fl_186, fl_187, fl_188, fl_189, \
                         fl_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * fl_183[k];

        t_49[k] = f_0 * fl_184[k];

        t_50[k] = f_0 * fl_185[k];

        t_51[k] = f_0 * fl_186[k];

        t_52[k] = f_0 * fl_187[k];

        t_53[k] = f_0 * fl_188[k];

        t_54[k] = f_0 * fl_189[k];

        t_55[k] = f_0 * fl_190[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, fl_191, fl_192, \
                         fl_193, fl_194, fl_195, fl_196, fl_197, \
                         fl_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * fl_191[k];

        t_57[k] = f_0 * fl_192[k];

        t_58[k] = f_0 * fl_193[k];

        t_59[k] = f_0 * fl_194[k];

        t_60[k] = f_0 * fl_195[k];

        t_61[k] = f_0 * fl_196[k];

        t_62[k] = f_0 * fl_197[k];

        t_63[k] = f_0 * fl_198[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, fl_199, fl_200, \
                         fl_201, fl_202, fl_203, fl_204, fl_205, \
                         fl_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * fl_199[k];

        t_65[k] = f_0 * fl_200[k];

        t_66[k] = f_0 * fl_201[k];

        t_67[k] = f_0 * fl_202[k];

        t_68[k] = f_0 * fl_203[k];

        t_69[k] = f_0 * fl_204[k];

        t_70[k] = f_0 * fl_205[k];

        t_71[k] = f_0 * fl_206[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, t_78, t_79, fl_207, fl_208, \
                         fl_209, fl_210, fl_211, fl_212, fl_213, \
                         fl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * fl_207[k];

        t_73[k] = f_0 * fl_208[k];

        t_74[k] = f_0 * fl_209[k];

        t_75[k] = f_0 * fl_210[k];

        t_76[k] = f_0 * fl_211[k];

        t_77[k] = f_0 * fl_212[k];

        t_78[k] = f_0 * fl_213[k];

        t_79[k] = f_0 * fl_214[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, t_86, t_87, fl_215, fl_216, \
                         fl_217, fl_218, fl_219, fl_220, fl_221, \
                         fl_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * fl_215[k];

        t_81[k] = f_0 * fl_216[k];

        t_82[k] = f_0 * fl_217[k];

        t_83[k] = f_0 * fl_218[k];

        t_84[k] = f_0 * fl_219[k];

        t_85[k] = f_0 * fl_220[k];

        t_86[k] = f_0 * fl_221[k];

        t_87[k] = f_0 * fl_222[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, t_93, pl_0, pl_1, pl_2, pl_3, fl_223, \
                         fl_224, fl_225, fl_226, fl_227, fl_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_0 * fl_223[k];

        t_89[k] = f_0 * fl_224[k];

        t_90[k] = -pl_0[k]
                  + f_0 * fl_225[k];

        t_91[k] = -pl_1[k]
                  + f_0 * fl_226[k];

        t_92[k] = -pl_2[k]
                  + f_0 * fl_227[k];

        t_93[k] = -pl_3[k]
                  + f_0 * fl_228[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pl_4, pl_5, pl_6, pl_7, pl_8, fl_229, \
                         fl_230, fl_231, fl_232, fl_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -pl_4[k]
                  + f_0 * fl_229[k];

        t_95[k] = -pl_5[k]
                  + f_0 * fl_230[k];

        t_96[k] = -pl_6[k]
                  + f_0 * fl_231[k];

        t_97[k] = -pl_7[k]
                  + f_0 * fl_232[k];

        t_98[k] = -pl_8[k]
                  + f_0 * fl_233[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pl_9, pl_10, pl_11, pl_12, pl_13, \
                         fl_234, fl_235, fl_236, fl_237, fl_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -pl_9[k]
                  + f_0 * fl_234[k];

        t_100[k] = -pl_10[k]
                   + f_0 * fl_235[k];

        t_101[k] = -pl_11[k]
                   + f_0 * fl_236[k];

        t_102[k] = -pl_12[k]
                   + f_0 * fl_237[k];

        t_103[k] = -pl_13[k]
                   + f_0 * fl_238[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pl_14, pl_15, pl_16, pl_17, pl_18, \
                         fl_239, fl_240, fl_241, fl_242, fl_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -pl_14[k]
                   + f_0 * fl_239[k];

        t_105[k] = -pl_15[k]
                   + f_0 * fl_240[k];

        t_106[k] = -pl_16[k]
                   + f_0 * fl_241[k];

        t_107[k] = -pl_17[k]
                   + f_0 * fl_242[k];

        t_108[k] = -pl_18[k]
                   + f_0 * fl_243[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pl_19, pl_20, pl_21, pl_22, pl_23, \
                         fl_244, fl_245, fl_246, fl_247, fl_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -pl_19[k]
                   + f_0 * fl_244[k];

        t_110[k] = -pl_20[k]
                   + f_0 * fl_245[k];

        t_111[k] = -pl_21[k]
                   + f_0 * fl_246[k];

        t_112[k] = -pl_22[k]
                   + f_0 * fl_247[k];

        t_113[k] = -pl_23[k]
                   + f_0 * fl_248[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pl_24, pl_25, pl_26, pl_27, pl_28, \
                         fl_249, fl_250, fl_251, fl_252, fl_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -pl_24[k]
                   + f_0 * fl_249[k];

        t_115[k] = -pl_25[k]
                   + f_0 * fl_250[k];

        t_116[k] = -pl_26[k]
                   + f_0 * fl_251[k];

        t_117[k] = -pl_27[k]
                   + f_0 * fl_252[k];

        t_118[k] = -pl_28[k]
                   + f_0 * fl_253[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, pl_29, pl_30, pl_31, pl_32, pl_33, \
                         fl_254, fl_255, fl_256, fl_257, fl_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -pl_29[k]
                   + f_0 * fl_254[k];

        t_120[k] = -pl_30[k]
                   + f_0 * fl_255[k];

        t_121[k] = -pl_31[k]
                   + f_0 * fl_256[k];

        t_122[k] = -pl_32[k]
                   + f_0 * fl_257[k];

        t_123[k] = -pl_33[k]
                   + f_0 * fl_258[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, pl_34, pl_35, pl_36, pl_37, pl_38, \
                         fl_259, fl_260, fl_261, fl_262, fl_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -pl_34[k]
                   + f_0 * fl_259[k];

        t_125[k] = -pl_35[k]
                   + f_0 * fl_260[k];

        t_126[k] = -pl_36[k]
                   + f_0 * fl_261[k];

        t_127[k] = -pl_37[k]
                   + f_0 * fl_262[k];

        t_128[k] = -pl_38[k]
                   + f_0 * fl_263[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pl_39, pl_40, pl_41, pl_42, pl_43, \
                         fl_264, fl_265, fl_266, fl_267, fl_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -pl_39[k]
                   + f_0 * fl_264[k];

        t_130[k] = -pl_40[k]
                   + f_0 * fl_265[k];

        t_131[k] = -pl_41[k]
                   + f_0 * fl_266[k];

        t_132[k] = -pl_42[k]
                   + f_0 * fl_267[k];

        t_133[k] = -pl_43[k]
                   + f_0 * fl_268[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, t_139, t_140, pl_44, fl_269, \
                         fl_315, fl_316, fl_317, fl_318, fl_319, \
                         fl_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -pl_44[k]
                   + f_0 * fl_269[k];

        t_135[k] = f_0 * fl_315[k];

        t_136[k] = f_0 * fl_316[k];

        t_137[k] = f_0 * fl_317[k];

        t_138[k] = f_0 * fl_318[k];

        t_139[k] = f_0 * fl_319[k];

        t_140[k] = f_0 * fl_320[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, t_146, t_147, t_148, fl_321, \
                         fl_322, fl_323, fl_324, fl_325, fl_326, fl_327, \
                         fl_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_0 * fl_321[k];

        t_142[k] = f_0 * fl_322[k];

        t_143[k] = f_0 * fl_323[k];

        t_144[k] = f_0 * fl_324[k];

        t_145[k] = f_0 * fl_325[k];

        t_146[k] = f_0 * fl_326[k];

        t_147[k] = f_0 * fl_327[k];

        t_148[k] = f_0 * fl_328[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, t_154, t_155, t_156, fl_329, \
                         fl_330, fl_331, fl_332, fl_333, fl_334, fl_335, \
                         fl_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * fl_329[k];

        t_150[k] = f_0 * fl_330[k];

        t_151[k] = f_0 * fl_331[k];

        t_152[k] = f_0 * fl_332[k];

        t_153[k] = f_0 * fl_333[k];

        t_154[k] = f_0 * fl_334[k];

        t_155[k] = f_0 * fl_335[k];

        t_156[k] = f_0 * fl_336[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, t_162, t_163, t_164, fl_337, \
                         fl_338, fl_339, fl_340, fl_341, fl_342, fl_343, \
                         fl_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_0 * fl_337[k];

        t_158[k] = f_0 * fl_338[k];

        t_159[k] = f_0 * fl_339[k];

        t_160[k] = f_0 * fl_340[k];

        t_161[k] = f_0 * fl_341[k];

        t_162[k] = f_0 * fl_342[k];

        t_163[k] = f_0 * fl_343[k];

        t_164[k] = f_0 * fl_344[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, t_170, t_171, t_172, fl_345, \
                         fl_346, fl_347, fl_348, fl_349, fl_350, fl_351, \
                         fl_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_0 * fl_345[k];

        t_166[k] = f_0 * fl_346[k];

        t_167[k] = f_0 * fl_347[k];

        t_168[k] = f_0 * fl_348[k];

        t_169[k] = f_0 * fl_349[k];

        t_170[k] = f_0 * fl_350[k];

        t_171[k] = f_0 * fl_351[k];

        t_172[k] = f_0 * fl_352[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, t_179, fl_353, fl_354, \
                         fl_355, fl_356, fl_357, fl_358, fl_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_0 * fl_353[k];

        t_174[k] = f_0 * fl_354[k];

        t_175[k] = f_0 * fl_355[k];

        t_176[k] = f_0 * fl_356[k];

        t_177[k] = f_0 * fl_357[k];

        t_178[k] = f_0 * fl_358[k];

        t_179[k] = f_0 * fl_359[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pl_45, pl_46, pl_47, pl_48, pl_49, \
                         fl_360, fl_361, fl_362, fl_363, fl_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -pl_45[k]
                   + f_0 * fl_360[k];

        t_181[k] = -pl_46[k]
                   + f_0 * fl_361[k];

        t_182[k] = -pl_47[k]
                   + f_0 * fl_362[k];

        t_183[k] = -pl_48[k]
                   + f_0 * fl_363[k];

        t_184[k] = -pl_49[k]
                   + f_0 * fl_364[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, pl_50, pl_51, pl_52, pl_53, pl_54, \
                         fl_365, fl_366, fl_367, fl_368, fl_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -pl_50[k]
                   + f_0 * fl_365[k];

        t_186[k] = -pl_51[k]
                   + f_0 * fl_366[k];

        t_187[k] = -pl_52[k]
                   + f_0 * fl_367[k];

        t_188[k] = -pl_53[k]
                   + f_0 * fl_368[k];

        t_189[k] = -pl_54[k]
                   + f_0 * fl_369[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, pl_55, pl_56, pl_57, pl_58, pl_59, \
                         fl_370, fl_371, fl_372, fl_373, fl_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -pl_55[k]
                   + f_0 * fl_370[k];

        t_191[k] = -pl_56[k]
                   + f_0 * fl_371[k];

        t_192[k] = -pl_57[k]
                   + f_0 * fl_372[k];

        t_193[k] = -pl_58[k]
                   + f_0 * fl_373[k];

        t_194[k] = -pl_59[k]
                   + f_0 * fl_374[k];
    }
}

static auto
compute_prim_geom_10_dl_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t pl, const size_t fl,
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

    const auto *pl_60 = buffer.data(pl + 60);
    const auto *pl_61 = buffer.data(pl + 61);
    const auto *pl_62 = buffer.data(pl + 62);
    const auto *pl_63 = buffer.data(pl + 63);
    const auto *pl_64 = buffer.data(pl + 64);
    const auto *pl_65 = buffer.data(pl + 65);
    const auto *pl_66 = buffer.data(pl + 66);
    const auto *pl_67 = buffer.data(pl + 67);
    const auto *pl_68 = buffer.data(pl + 68);
    const auto *pl_69 = buffer.data(pl + 69);
    const auto *pl_70 = buffer.data(pl + 70);
    const auto *pl_71 = buffer.data(pl + 71);
    const auto *pl_72 = buffer.data(pl + 72);
    const auto *pl_73 = buffer.data(pl + 73);
    const auto *pl_74 = buffer.data(pl + 74);
    const auto *pl_75 = buffer.data(pl + 75);
    const auto *pl_76 = buffer.data(pl + 76);
    const auto *pl_77 = buffer.data(pl + 77);
    const auto *pl_78 = buffer.data(pl + 78);
    const auto *pl_79 = buffer.data(pl + 79);
    const auto *pl_80 = buffer.data(pl + 80);
    const auto *pl_81 = buffer.data(pl + 81);
    const auto *pl_82 = buffer.data(pl + 82);
    const auto *pl_83 = buffer.data(pl + 83);
    const auto *pl_84 = buffer.data(pl + 84);
    const auto *pl_85 = buffer.data(pl + 85);
    const auto *pl_86 = buffer.data(pl + 86);
    const auto *pl_87 = buffer.data(pl + 87);
    const auto *pl_88 = buffer.data(pl + 88);
    const auto *pl_89 = buffer.data(pl + 89);
    const auto *pl_90 = buffer.data(pl + 90);
    const auto *pl_91 = buffer.data(pl + 91);
    const auto *pl_92 = buffer.data(pl + 92);
    const auto *pl_93 = buffer.data(pl + 93);
    const auto *pl_94 = buffer.data(pl + 94);
    const auto *pl_95 = buffer.data(pl + 95);
    const auto *pl_96 = buffer.data(pl + 96);
    const auto *pl_97 = buffer.data(pl + 97);
    const auto *pl_98 = buffer.data(pl + 98);
    const auto *pl_99 = buffer.data(pl + 99);
    const auto *pl_100 = buffer.data(pl + 100);
    const auto *pl_101 = buffer.data(pl + 101);
    const auto *pl_102 = buffer.data(pl + 102);
    const auto *pl_103 = buffer.data(pl + 103);
    const auto *pl_104 = buffer.data(pl + 104);
    const auto *pl_105 = buffer.data(pl + 105);
    const auto *pl_106 = buffer.data(pl + 106);
    const auto *pl_107 = buffer.data(pl + 107);
    const auto *pl_108 = buffer.data(pl + 108);
    const auto *pl_109 = buffer.data(pl + 109);
    const auto *pl_110 = buffer.data(pl + 110);
    const auto *pl_111 = buffer.data(pl + 111);
    const auto *pl_112 = buffer.data(pl + 112);
    const auto *pl_113 = buffer.data(pl + 113);
    const auto *pl_114 = buffer.data(pl + 114);
    const auto *pl_115 = buffer.data(pl + 115);
    const auto *pl_116 = buffer.data(pl + 116);
    const auto *pl_117 = buffer.data(pl + 117);
    const auto *pl_118 = buffer.data(pl + 118);
    const auto *pl_119 = buffer.data(pl + 119);
    const auto *pl_120 = buffer.data(pl + 120);
    const auto *pl_121 = buffer.data(pl + 121);
    const auto *pl_122 = buffer.data(pl + 122);
    const auto *pl_123 = buffer.data(pl + 123);
    const auto *pl_124 = buffer.data(pl + 124);
    const auto *pl_125 = buffer.data(pl + 125);
    const auto *pl_126 = buffer.data(pl + 126);
    const auto *pl_127 = buffer.data(pl + 127);
    const auto *pl_128 = buffer.data(pl + 128);
    const auto *pl_129 = buffer.data(pl + 129);
    const auto *pl_130 = buffer.data(pl + 130);
    const auto *pl_131 = buffer.data(pl + 131);
    const auto *pl_132 = buffer.data(pl + 132);
    const auto *pl_133 = buffer.data(pl + 133);
    const auto *pl_134 = buffer.data(pl + 134);

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

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, pl_60, pl_61, pl_62, pl_63, pl_64, \
                         fl_375, fl_376, fl_377, fl_378, fl_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -pl_60[k]
                   + f_0 * fl_375[k];

        t_196[k] = -pl_61[k]
                   + f_0 * fl_376[k];

        t_197[k] = -pl_62[k]
                   + f_0 * fl_377[k];

        t_198[k] = -pl_63[k]
                   + f_0 * fl_378[k];

        t_199[k] = -pl_64[k]
                   + f_0 * fl_379[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pl_65, pl_66, pl_67, pl_68, pl_69, \
                         fl_380, fl_381, fl_382, fl_383, fl_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -pl_65[k]
                   + f_0 * fl_380[k];

        t_201[k] = -pl_66[k]
                   + f_0 * fl_381[k];

        t_202[k] = -pl_67[k]
                   + f_0 * fl_382[k];

        t_203[k] = -pl_68[k]
                   + f_0 * fl_383[k];

        t_204[k] = -pl_69[k]
                   + f_0 * fl_384[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, pl_70, pl_71, pl_72, pl_73, pl_74, \
                         fl_385, fl_386, fl_387, fl_388, fl_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -pl_70[k]
                   + f_0 * fl_385[k];

        t_206[k] = -pl_71[k]
                   + f_0 * fl_386[k];

        t_207[k] = -pl_72[k]
                   + f_0 * fl_387[k];

        t_208[k] = -pl_73[k]
                   + f_0 * fl_388[k];

        t_209[k] = -pl_74[k]
                   + f_0 * fl_389[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pl_75, pl_76, pl_77, pl_78, pl_79, \
                         fl_390, fl_391, fl_392, fl_393, fl_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -pl_75[k]
                   + f_0 * fl_390[k];

        t_211[k] = -pl_76[k]
                   + f_0 * fl_391[k];

        t_212[k] = -pl_77[k]
                   + f_0 * fl_392[k];

        t_213[k] = -pl_78[k]
                   + f_0 * fl_393[k];

        t_214[k] = -pl_79[k]
                   + f_0 * fl_394[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, pl_80, pl_81, pl_82, pl_83, pl_84, \
                         fl_395, fl_396, fl_397, fl_398, fl_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -pl_80[k]
                   + f_0 * fl_395[k];

        t_216[k] = -pl_81[k]
                   + f_0 * fl_396[k];

        t_217[k] = -pl_82[k]
                   + f_0 * fl_397[k];

        t_218[k] = -pl_83[k]
                   + f_0 * fl_398[k];

        t_219[k] = -pl_84[k]
                   + f_0 * fl_399[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, pl_85, pl_86, pl_87, pl_88, pl_89, \
                         fl_400, fl_401, fl_402, fl_403, fl_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -pl_85[k]
                   + f_0 * fl_400[k];

        t_221[k] = -pl_86[k]
                   + f_0 * fl_401[k];

        t_222[k] = -pl_87[k]
                   + f_0 * fl_402[k];

        t_223[k] = -pl_88[k]
                   + f_0 * fl_403[k];

        t_224[k] = -pl_89[k]
                   + f_0 * fl_404[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, pl_90, pl_91, pl_92, pl_93, pl_94, \
                         fl_405, fl_406, fl_407, fl_408, fl_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -2.0 * pl_90[k]
                   + f_0 * fl_405[k];

        t_226[k] = -2.0 * pl_91[k]
                   + f_0 * fl_406[k];

        t_227[k] = -2.0 * pl_92[k]
                   + f_0 * fl_407[k];

        t_228[k] = -2.0 * pl_93[k]
                   + f_0 * fl_408[k];

        t_229[k] = -2.0 * pl_94[k]
                   + f_0 * fl_409[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, pl_95, pl_96, pl_97, pl_98, pl_99, \
                         fl_410, fl_411, fl_412, fl_413, fl_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -2.0 * pl_95[k]
                   + f_0 * fl_410[k];

        t_231[k] = -2.0 * pl_96[k]
                   + f_0 * fl_411[k];

        t_232[k] = -2.0 * pl_97[k]
                   + f_0 * fl_412[k];

        t_233[k] = -2.0 * pl_98[k]
                   + f_0 * fl_413[k];

        t_234[k] = -2.0 * pl_99[k]
                   + f_0 * fl_414[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, pl_100, pl_101, pl_102, pl_103, \
                         pl_104, fl_415, fl_416, fl_417, fl_418, \
                         fl_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -2.0 * pl_100[k]
                   + f_0 * fl_415[k];

        t_236[k] = -2.0 * pl_101[k]
                   + f_0 * fl_416[k];

        t_237[k] = -2.0 * pl_102[k]
                   + f_0 * fl_417[k];

        t_238[k] = -2.0 * pl_103[k]
                   + f_0 * fl_418[k];

        t_239[k] = -2.0 * pl_104[k]
                   + f_0 * fl_419[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pl_105, pl_106, pl_107, pl_108, \
                         pl_109, fl_420, fl_421, fl_422, fl_423, \
                         fl_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -2.0 * pl_105[k]
                   + f_0 * fl_420[k];

        t_241[k] = -2.0 * pl_106[k]
                   + f_0 * fl_421[k];

        t_242[k] = -2.0 * pl_107[k]
                   + f_0 * fl_422[k];

        t_243[k] = -2.0 * pl_108[k]
                   + f_0 * fl_423[k];

        t_244[k] = -2.0 * pl_109[k]
                   + f_0 * fl_424[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, pl_110, pl_111, pl_112, pl_113, \
                         pl_114, fl_425, fl_426, fl_427, fl_428, \
                         fl_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -2.0 * pl_110[k]
                   + f_0 * fl_425[k];

        t_246[k] = -2.0 * pl_111[k]
                   + f_0 * fl_426[k];

        t_247[k] = -2.0 * pl_112[k]
                   + f_0 * fl_427[k];

        t_248[k] = -2.0 * pl_113[k]
                   + f_0 * fl_428[k];

        t_249[k] = -2.0 * pl_114[k]
                   + f_0 * fl_429[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, pl_115, pl_116, pl_117, pl_118, \
                         pl_119, fl_430, fl_431, fl_432, fl_433, \
                         fl_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * pl_115[k]
                   + f_0 * fl_430[k];

        t_251[k] = -2.0 * pl_116[k]
                   + f_0 * fl_431[k];

        t_252[k] = -2.0 * pl_117[k]
                   + f_0 * fl_432[k];

        t_253[k] = -2.0 * pl_118[k]
                   + f_0 * fl_433[k];

        t_254[k] = -2.0 * pl_119[k]
                   + f_0 * fl_434[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pl_120, pl_121, pl_122, pl_123, \
                         pl_124, fl_435, fl_436, fl_437, fl_438, \
                         fl_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * pl_120[k]
                   + f_0 * fl_435[k];

        t_256[k] = -2.0 * pl_121[k]
                   + f_0 * fl_436[k];

        t_257[k] = -2.0 * pl_122[k]
                   + f_0 * fl_437[k];

        t_258[k] = -2.0 * pl_123[k]
                   + f_0 * fl_438[k];

        t_259[k] = -2.0 * pl_124[k]
                   + f_0 * fl_439[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, pl_125, pl_126, pl_127, pl_128, \
                         pl_129, fl_440, fl_441, fl_442, fl_443, \
                         fl_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * pl_125[k]
                   + f_0 * fl_440[k];

        t_261[k] = -2.0 * pl_126[k]
                   + f_0 * fl_441[k];

        t_262[k] = -2.0 * pl_127[k]
                   + f_0 * fl_442[k];

        t_263[k] = -2.0 * pl_128[k]
                   + f_0 * fl_443[k];

        t_264[k] = -2.0 * pl_129[k]
                   + f_0 * fl_444[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, pl_130, pl_131, pl_132, pl_133, \
                         pl_134, fl_445, fl_446, fl_447, fl_448, \
                         fl_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * pl_130[k]
                   + f_0 * fl_445[k];

        t_266[k] = -2.0 * pl_131[k]
                   + f_0 * fl_446[k];

        t_267[k] = -2.0 * pl_132[k]
                   + f_0 * fl_447[k];

        t_268[k] = -2.0 * pl_133[k]
                   + f_0 * fl_448[k];

        t_269[k] = -2.0 * pl_134[k]
                   + f_0 * fl_449[k];
    }
}

auto
compute_prim_geom_10_dl_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t pl, const size_t fl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_dl_electron_repulsion_2_piece0(buffer, target, pl, fl, ncols, alpha);

    compute_prim_geom_10_dl_electron_repulsion_2_piece1(buffer, target, pl, fl, ncols, alpha);
}

}  // namespace simdt2ceri
