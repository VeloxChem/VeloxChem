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


#include "SimdElectronRepulsionGeom10VrrRecKH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_kh_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
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

    const auto *ih_0 = buffer.data(ih + 0);
    const auto *ih_1 = buffer.data(ih + 1);
    const auto *ih_2 = buffer.data(ih + 2);
    const auto *ih_3 = buffer.data(ih + 3);
    const auto *ih_4 = buffer.data(ih + 4);
    const auto *ih_5 = buffer.data(ih + 5);
    const auto *ih_6 = buffer.data(ih + 6);
    const auto *ih_7 = buffer.data(ih + 7);
    const auto *ih_8 = buffer.data(ih + 8);
    const auto *ih_9 = buffer.data(ih + 9);
    const auto *ih_10 = buffer.data(ih + 10);
    const auto *ih_11 = buffer.data(ih + 11);
    const auto *ih_12 = buffer.data(ih + 12);
    const auto *ih_13 = buffer.data(ih + 13);
    const auto *ih_14 = buffer.data(ih + 14);
    const auto *ih_15 = buffer.data(ih + 15);
    const auto *ih_16 = buffer.data(ih + 16);
    const auto *ih_17 = buffer.data(ih + 17);
    const auto *ih_18 = buffer.data(ih + 18);
    const auto *ih_19 = buffer.data(ih + 19);
    const auto *ih_20 = buffer.data(ih + 20);
    const auto *ih_21 = buffer.data(ih + 21);
    const auto *ih_22 = buffer.data(ih + 22);
    const auto *ih_23 = buffer.data(ih + 23);
    const auto *ih_24 = buffer.data(ih + 24);
    const auto *ih_25 = buffer.data(ih + 25);
    const auto *ih_26 = buffer.data(ih + 26);
    const auto *ih_27 = buffer.data(ih + 27);
    const auto *ih_28 = buffer.data(ih + 28);
    const auto *ih_29 = buffer.data(ih + 29);
    const auto *ih_30 = buffer.data(ih + 30);
    const auto *ih_31 = buffer.data(ih + 31);
    const auto *ih_32 = buffer.data(ih + 32);
    const auto *ih_33 = buffer.data(ih + 33);
    const auto *ih_34 = buffer.data(ih + 34);
    const auto *ih_35 = buffer.data(ih + 35);
    const auto *ih_36 = buffer.data(ih + 36);
    const auto *ih_37 = buffer.data(ih + 37);
    const auto *ih_38 = buffer.data(ih + 38);
    const auto *ih_39 = buffer.data(ih + 39);
    const auto *ih_40 = buffer.data(ih + 40);
    const auto *ih_41 = buffer.data(ih + 41);
    const auto *ih_42 = buffer.data(ih + 42);
    const auto *ih_43 = buffer.data(ih + 43);
    const auto *ih_44 = buffer.data(ih + 44);
    const auto *ih_45 = buffer.data(ih + 45);
    const auto *ih_46 = buffer.data(ih + 46);
    const auto *ih_47 = buffer.data(ih + 47);
    const auto *ih_48 = buffer.data(ih + 48);
    const auto *ih_49 = buffer.data(ih + 49);
    const auto *ih_50 = buffer.data(ih + 50);
    const auto *ih_51 = buffer.data(ih + 51);
    const auto *ih_52 = buffer.data(ih + 52);
    const auto *ih_53 = buffer.data(ih + 53);
    const auto *ih_54 = buffer.data(ih + 54);
    const auto *ih_55 = buffer.data(ih + 55);
    const auto *ih_56 = buffer.data(ih + 56);
    const auto *ih_57 = buffer.data(ih + 57);
    const auto *ih_58 = buffer.data(ih + 58);
    const auto *ih_59 = buffer.data(ih + 59);
    const auto *ih_60 = buffer.data(ih + 60);
    const auto *ih_61 = buffer.data(ih + 61);
    const auto *ih_62 = buffer.data(ih + 62);
    const auto *ih_63 = buffer.data(ih + 63);
    const auto *ih_64 = buffer.data(ih + 64);
    const auto *ih_65 = buffer.data(ih + 65);
    const auto *ih_66 = buffer.data(ih + 66);
    const auto *ih_67 = buffer.data(ih + 67);
    const auto *ih_68 = buffer.data(ih + 68);
    const auto *ih_69 = buffer.data(ih + 69);
    const auto *ih_70 = buffer.data(ih + 70);
    const auto *ih_71 = buffer.data(ih + 71);
    const auto *ih_72 = buffer.data(ih + 72);
    const auto *ih_73 = buffer.data(ih + 73);
    const auto *ih_74 = buffer.data(ih + 74);
    const auto *ih_75 = buffer.data(ih + 75);
    const auto *ih_76 = buffer.data(ih + 76);
    const auto *ih_77 = buffer.data(ih + 77);
    const auto *ih_78 = buffer.data(ih + 78);
    const auto *ih_79 = buffer.data(ih + 79);
    const auto *ih_80 = buffer.data(ih + 80);
    const auto *ih_81 = buffer.data(ih + 81);
    const auto *ih_82 = buffer.data(ih + 82);
    const auto *ih_83 = buffer.data(ih + 83);
    const auto *ih_84 = buffer.data(ih + 84);
    const auto *ih_85 = buffer.data(ih + 85);
    const auto *ih_86 = buffer.data(ih + 86);
    const auto *ih_87 = buffer.data(ih + 87);
    const auto *ih_88 = buffer.data(ih + 88);
    const auto *ih_89 = buffer.data(ih + 89);
    const auto *ih_90 = buffer.data(ih + 90);
    const auto *ih_91 = buffer.data(ih + 91);
    const auto *ih_92 = buffer.data(ih + 92);
    const auto *ih_93 = buffer.data(ih + 93);
    const auto *ih_94 = buffer.data(ih + 94);
    const auto *ih_95 = buffer.data(ih + 95);
    const auto *ih_96 = buffer.data(ih + 96);
    const auto *ih_97 = buffer.data(ih + 97);
    const auto *ih_98 = buffer.data(ih + 98);
    const auto *ih_99 = buffer.data(ih + 99);
    const auto *ih_100 = buffer.data(ih + 100);
    const auto *ih_101 = buffer.data(ih + 101);
    const auto *ih_102 = buffer.data(ih + 102);
    const auto *ih_103 = buffer.data(ih + 103);
    const auto *ih_104 = buffer.data(ih + 104);
    const auto *ih_105 = buffer.data(ih + 105);
    const auto *ih_106 = buffer.data(ih + 106);
    const auto *ih_107 = buffer.data(ih + 107);
    const auto *ih_108 = buffer.data(ih + 108);
    const auto *ih_109 = buffer.data(ih + 109);
    const auto *ih_110 = buffer.data(ih + 110);
    const auto *ih_111 = buffer.data(ih + 111);
    const auto *ih_112 = buffer.data(ih + 112);
    const auto *ih_113 = buffer.data(ih + 113);
    const auto *ih_114 = buffer.data(ih + 114);
    const auto *ih_115 = buffer.data(ih + 115);
    const auto *ih_116 = buffer.data(ih + 116);
    const auto *ih_117 = buffer.data(ih + 117);
    const auto *ih_118 = buffer.data(ih + 118);
    const auto *ih_119 = buffer.data(ih + 119);
    const auto *ih_120 = buffer.data(ih + 120);
    const auto *ih_121 = buffer.data(ih + 121);
    const auto *ih_122 = buffer.data(ih + 122);
    const auto *ih_123 = buffer.data(ih + 123);
    const auto *ih_124 = buffer.data(ih + 124);
    const auto *ih_125 = buffer.data(ih + 125);
    const auto *ih_126 = buffer.data(ih + 126);
    const auto *ih_127 = buffer.data(ih + 127);
    const auto *ih_128 = buffer.data(ih + 128);
    const auto *ih_129 = buffer.data(ih + 129);
    const auto *ih_130 = buffer.data(ih + 130);
    const auto *ih_131 = buffer.data(ih + 131);
    const auto *ih_132 = buffer.data(ih + 132);
    const auto *ih_133 = buffer.data(ih + 133);
    const auto *ih_134 = buffer.data(ih + 134);
    const auto *ih_135 = buffer.data(ih + 135);
    const auto *ih_136 = buffer.data(ih + 136);
    const auto *ih_137 = buffer.data(ih + 137);
    const auto *ih_138 = buffer.data(ih + 138);
    const auto *ih_139 = buffer.data(ih + 139);
    const auto *ih_140 = buffer.data(ih + 140);
    const auto *ih_141 = buffer.data(ih + 141);
    const auto *ih_142 = buffer.data(ih + 142);
    const auto *ih_143 = buffer.data(ih + 143);
    const auto *ih_144 = buffer.data(ih + 144);
    const auto *ih_145 = buffer.data(ih + 145);
    const auto *ih_146 = buffer.data(ih + 146);
    const auto *ih_147 = buffer.data(ih + 147);
    const auto *ih_148 = buffer.data(ih + 148);
    const auto *ih_149 = buffer.data(ih + 149);

    const auto *lh_0 = buffer.data(lh + 0);
    const auto *lh_1 = buffer.data(lh + 1);
    const auto *lh_2 = buffer.data(lh + 2);
    const auto *lh_3 = buffer.data(lh + 3);
    const auto *lh_4 = buffer.data(lh + 4);
    const auto *lh_5 = buffer.data(lh + 5);
    const auto *lh_6 = buffer.data(lh + 6);
    const auto *lh_7 = buffer.data(lh + 7);
    const auto *lh_8 = buffer.data(lh + 8);
    const auto *lh_9 = buffer.data(lh + 9);
    const auto *lh_10 = buffer.data(lh + 10);
    const auto *lh_11 = buffer.data(lh + 11);
    const auto *lh_12 = buffer.data(lh + 12);
    const auto *lh_13 = buffer.data(lh + 13);
    const auto *lh_14 = buffer.data(lh + 14);
    const auto *lh_15 = buffer.data(lh + 15);
    const auto *lh_16 = buffer.data(lh + 16);
    const auto *lh_17 = buffer.data(lh + 17);
    const auto *lh_18 = buffer.data(lh + 18);
    const auto *lh_19 = buffer.data(lh + 19);
    const auto *lh_20 = buffer.data(lh + 20);
    const auto *lh_21 = buffer.data(lh + 21);
    const auto *lh_22 = buffer.data(lh + 22);
    const auto *lh_23 = buffer.data(lh + 23);
    const auto *lh_24 = buffer.data(lh + 24);
    const auto *lh_25 = buffer.data(lh + 25);
    const auto *lh_26 = buffer.data(lh + 26);
    const auto *lh_27 = buffer.data(lh + 27);
    const auto *lh_28 = buffer.data(lh + 28);
    const auto *lh_29 = buffer.data(lh + 29);
    const auto *lh_30 = buffer.data(lh + 30);
    const auto *lh_31 = buffer.data(lh + 31);
    const auto *lh_32 = buffer.data(lh + 32);
    const auto *lh_33 = buffer.data(lh + 33);
    const auto *lh_34 = buffer.data(lh + 34);
    const auto *lh_35 = buffer.data(lh + 35);
    const auto *lh_36 = buffer.data(lh + 36);
    const auto *lh_37 = buffer.data(lh + 37);
    const auto *lh_38 = buffer.data(lh + 38);
    const auto *lh_39 = buffer.data(lh + 39);
    const auto *lh_40 = buffer.data(lh + 40);
    const auto *lh_41 = buffer.data(lh + 41);
    const auto *lh_42 = buffer.data(lh + 42);
    const auto *lh_43 = buffer.data(lh + 43);
    const auto *lh_44 = buffer.data(lh + 44);
    const auto *lh_45 = buffer.data(lh + 45);
    const auto *lh_46 = buffer.data(lh + 46);
    const auto *lh_47 = buffer.data(lh + 47);
    const auto *lh_48 = buffer.data(lh + 48);
    const auto *lh_49 = buffer.data(lh + 49);
    const auto *lh_50 = buffer.data(lh + 50);
    const auto *lh_51 = buffer.data(lh + 51);
    const auto *lh_52 = buffer.data(lh + 52);
    const auto *lh_53 = buffer.data(lh + 53);
    const auto *lh_54 = buffer.data(lh + 54);
    const auto *lh_55 = buffer.data(lh + 55);
    const auto *lh_56 = buffer.data(lh + 56);
    const auto *lh_57 = buffer.data(lh + 57);
    const auto *lh_58 = buffer.data(lh + 58);
    const auto *lh_59 = buffer.data(lh + 59);
    const auto *lh_60 = buffer.data(lh + 60);
    const auto *lh_61 = buffer.data(lh + 61);
    const auto *lh_62 = buffer.data(lh + 62);
    const auto *lh_63 = buffer.data(lh + 63);
    const auto *lh_64 = buffer.data(lh + 64);
    const auto *lh_65 = buffer.data(lh + 65);
    const auto *lh_66 = buffer.data(lh + 66);
    const auto *lh_67 = buffer.data(lh + 67);
    const auto *lh_68 = buffer.data(lh + 68);
    const auto *lh_69 = buffer.data(lh + 69);
    const auto *lh_70 = buffer.data(lh + 70);
    const auto *lh_71 = buffer.data(lh + 71);
    const auto *lh_72 = buffer.data(lh + 72);
    const auto *lh_73 = buffer.data(lh + 73);
    const auto *lh_74 = buffer.data(lh + 74);
    const auto *lh_75 = buffer.data(lh + 75);
    const auto *lh_76 = buffer.data(lh + 76);
    const auto *lh_77 = buffer.data(lh + 77);
    const auto *lh_78 = buffer.data(lh + 78);
    const auto *lh_79 = buffer.data(lh + 79);
    const auto *lh_80 = buffer.data(lh + 80);
    const auto *lh_81 = buffer.data(lh + 81);
    const auto *lh_82 = buffer.data(lh + 82);
    const auto *lh_83 = buffer.data(lh + 83);
    const auto *lh_84 = buffer.data(lh + 84);
    const auto *lh_85 = buffer.data(lh + 85);
    const auto *lh_86 = buffer.data(lh + 86);
    const auto *lh_87 = buffer.data(lh + 87);
    const auto *lh_88 = buffer.data(lh + 88);
    const auto *lh_89 = buffer.data(lh + 89);
    const auto *lh_90 = buffer.data(lh + 90);
    const auto *lh_91 = buffer.data(lh + 91);
    const auto *lh_92 = buffer.data(lh + 92);
    const auto *lh_93 = buffer.data(lh + 93);
    const auto *lh_94 = buffer.data(lh + 94);
    const auto *lh_95 = buffer.data(lh + 95);
    const auto *lh_96 = buffer.data(lh + 96);
    const auto *lh_97 = buffer.data(lh + 97);
    const auto *lh_98 = buffer.data(lh + 98);
    const auto *lh_99 = buffer.data(lh + 99);
    const auto *lh_100 = buffer.data(lh + 100);
    const auto *lh_101 = buffer.data(lh + 101);
    const auto *lh_102 = buffer.data(lh + 102);
    const auto *lh_103 = buffer.data(lh + 103);
    const auto *lh_104 = buffer.data(lh + 104);
    const auto *lh_105 = buffer.data(lh + 105);
    const auto *lh_106 = buffer.data(lh + 106);
    const auto *lh_107 = buffer.data(lh + 107);
    const auto *lh_108 = buffer.data(lh + 108);
    const auto *lh_109 = buffer.data(lh + 109);
    const auto *lh_110 = buffer.data(lh + 110);
    const auto *lh_111 = buffer.data(lh + 111);
    const auto *lh_112 = buffer.data(lh + 112);
    const auto *lh_113 = buffer.data(lh + 113);
    const auto *lh_114 = buffer.data(lh + 114);
    const auto *lh_115 = buffer.data(lh + 115);
    const auto *lh_116 = buffer.data(lh + 116);
    const auto *lh_117 = buffer.data(lh + 117);
    const auto *lh_118 = buffer.data(lh + 118);
    const auto *lh_119 = buffer.data(lh + 119);
    const auto *lh_120 = buffer.data(lh + 120);
    const auto *lh_121 = buffer.data(lh + 121);
    const auto *lh_122 = buffer.data(lh + 122);
    const auto *lh_123 = buffer.data(lh + 123);
    const auto *lh_124 = buffer.data(lh + 124);
    const auto *lh_125 = buffer.data(lh + 125);
    const auto *lh_126 = buffer.data(lh + 126);
    const auto *lh_127 = buffer.data(lh + 127);
    const auto *lh_128 = buffer.data(lh + 128);
    const auto *lh_129 = buffer.data(lh + 129);
    const auto *lh_130 = buffer.data(lh + 130);
    const auto *lh_131 = buffer.data(lh + 131);
    const auto *lh_132 = buffer.data(lh + 132);
    const auto *lh_133 = buffer.data(lh + 133);
    const auto *lh_134 = buffer.data(lh + 134);
    const auto *lh_135 = buffer.data(lh + 135);
    const auto *lh_136 = buffer.data(lh + 136);
    const auto *lh_137 = buffer.data(lh + 137);
    const auto *lh_138 = buffer.data(lh + 138);
    const auto *lh_139 = buffer.data(lh + 139);
    const auto *lh_140 = buffer.data(lh + 140);
    const auto *lh_141 = buffer.data(lh + 141);
    const auto *lh_142 = buffer.data(lh + 142);
    const auto *lh_143 = buffer.data(lh + 143);
    const auto *lh_144 = buffer.data(lh + 144);
    const auto *lh_145 = buffer.data(lh + 145);
    const auto *lh_146 = buffer.data(lh + 146);
    const auto *lh_147 = buffer.data(lh + 147);
    const auto *lh_148 = buffer.data(lh + 148);
    const auto *lh_149 = buffer.data(lh + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ih_0, ih_1, ih_2, ih_3, ih_4, lh_0, lh_1, \
                         lh_2, lh_3, lh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -7.0 * ih_0[k]
                 + f_0 * lh_0[k];

        t_1[k] = -7.0 * ih_1[k]
                 + f_0 * lh_1[k];

        t_2[k] = -7.0 * ih_2[k]
                 + f_0 * lh_2[k];

        t_3[k] = -7.0 * ih_3[k]
                 + f_0 * lh_3[k];

        t_4[k] = -7.0 * ih_4[k]
                 + f_0 * lh_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ih_5, ih_6, ih_7, ih_8, ih_9, lh_5, lh_6, \
                         lh_7, lh_8, lh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -7.0 * ih_5[k]
                 + f_0 * lh_5[k];

        t_6[k] = -7.0 * ih_6[k]
                 + f_0 * lh_6[k];

        t_7[k] = -7.0 * ih_7[k]
                 + f_0 * lh_7[k];

        t_8[k] = -7.0 * ih_8[k]
                 + f_0 * lh_8[k];

        t_9[k] = -7.0 * ih_9[k]
                 + f_0 * lh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ih_10, ih_11, ih_12, ih_13, ih_14, \
                         lh_10, lh_11, lh_12, lh_13, lh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -7.0 * ih_10[k]
                  + f_0 * lh_10[k];

        t_11[k] = -7.0 * ih_11[k]
                  + f_0 * lh_11[k];

        t_12[k] = -7.0 * ih_12[k]
                  + f_0 * lh_12[k];

        t_13[k] = -7.0 * ih_13[k]
                  + f_0 * lh_13[k];

        t_14[k] = -7.0 * ih_14[k]
                  + f_0 * lh_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ih_15, ih_16, ih_17, ih_18, ih_19, \
                         lh_15, lh_16, lh_17, lh_18, lh_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -7.0 * ih_15[k]
                  + f_0 * lh_15[k];

        t_16[k] = -7.0 * ih_16[k]
                  + f_0 * lh_16[k];

        t_17[k] = -7.0 * ih_17[k]
                  + f_0 * lh_17[k];

        t_18[k] = -7.0 * ih_18[k]
                  + f_0 * lh_18[k];

        t_19[k] = -7.0 * ih_19[k]
                  + f_0 * lh_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ih_20, ih_21, ih_22, ih_23, ih_24, \
                         lh_20, lh_21, lh_22, lh_23, lh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -7.0 * ih_20[k]
                  + f_0 * lh_20[k];

        t_21[k] = -6.0 * ih_21[k]
                  + f_0 * lh_21[k];

        t_22[k] = -6.0 * ih_22[k]
                  + f_0 * lh_22[k];

        t_23[k] = -6.0 * ih_23[k]
                  + f_0 * lh_23[k];

        t_24[k] = -6.0 * ih_24[k]
                  + f_0 * lh_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ih_25, ih_26, ih_27, ih_28, ih_29, \
                         lh_25, lh_26, lh_27, lh_28, lh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -6.0 * ih_25[k]
                  + f_0 * lh_25[k];

        t_26[k] = -6.0 * ih_26[k]
                  + f_0 * lh_26[k];

        t_27[k] = -6.0 * ih_27[k]
                  + f_0 * lh_27[k];

        t_28[k] = -6.0 * ih_28[k]
                  + f_0 * lh_28[k];

        t_29[k] = -6.0 * ih_29[k]
                  + f_0 * lh_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ih_30, ih_31, ih_32, ih_33, ih_34, \
                         lh_30, lh_31, lh_32, lh_33, lh_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -6.0 * ih_30[k]
                  + f_0 * lh_30[k];

        t_31[k] = -6.0 * ih_31[k]
                  + f_0 * lh_31[k];

        t_32[k] = -6.0 * ih_32[k]
                  + f_0 * lh_32[k];

        t_33[k] = -6.0 * ih_33[k]
                  + f_0 * lh_33[k];

        t_34[k] = -6.0 * ih_34[k]
                  + f_0 * lh_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ih_35, ih_36, ih_37, ih_38, ih_39, \
                         lh_35, lh_36, lh_37, lh_38, lh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -6.0 * ih_35[k]
                  + f_0 * lh_35[k];

        t_36[k] = -6.0 * ih_36[k]
                  + f_0 * lh_36[k];

        t_37[k] = -6.0 * ih_37[k]
                  + f_0 * lh_37[k];

        t_38[k] = -6.0 * ih_38[k]
                  + f_0 * lh_38[k];

        t_39[k] = -6.0 * ih_39[k]
                  + f_0 * lh_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ih_40, ih_41, ih_42, ih_43, ih_44, \
                         lh_40, lh_41, lh_42, lh_43, lh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -6.0 * ih_40[k]
                  + f_0 * lh_40[k];

        t_41[k] = -6.0 * ih_41[k]
                  + f_0 * lh_41[k];

        t_42[k] = -6.0 * ih_42[k]
                  + f_0 * lh_42[k];

        t_43[k] = -6.0 * ih_43[k]
                  + f_0 * lh_43[k];

        t_44[k] = -6.0 * ih_44[k]
                  + f_0 * lh_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ih_45, ih_46, ih_47, ih_48, ih_49, \
                         lh_45, lh_46, lh_47, lh_48, lh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -6.0 * ih_45[k]
                  + f_0 * lh_45[k];

        t_46[k] = -6.0 * ih_46[k]
                  + f_0 * lh_46[k];

        t_47[k] = -6.0 * ih_47[k]
                  + f_0 * lh_47[k];

        t_48[k] = -6.0 * ih_48[k]
                  + f_0 * lh_48[k];

        t_49[k] = -6.0 * ih_49[k]
                  + f_0 * lh_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ih_50, ih_51, ih_52, ih_53, ih_54, \
                         lh_50, lh_51, lh_52, lh_53, lh_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -6.0 * ih_50[k]
                  + f_0 * lh_50[k];

        t_51[k] = -6.0 * ih_51[k]
                  + f_0 * lh_51[k];

        t_52[k] = -6.0 * ih_52[k]
                  + f_0 * lh_52[k];

        t_53[k] = -6.0 * ih_53[k]
                  + f_0 * lh_53[k];

        t_54[k] = -6.0 * ih_54[k]
                  + f_0 * lh_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ih_55, ih_56, ih_57, ih_58, ih_59, \
                         lh_55, lh_56, lh_57, lh_58, lh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -6.0 * ih_55[k]
                  + f_0 * lh_55[k];

        t_56[k] = -6.0 * ih_56[k]
                  + f_0 * lh_56[k];

        t_57[k] = -6.0 * ih_57[k]
                  + f_0 * lh_57[k];

        t_58[k] = -6.0 * ih_58[k]
                  + f_0 * lh_58[k];

        t_59[k] = -6.0 * ih_59[k]
                  + f_0 * lh_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ih_60, ih_61, ih_62, ih_63, ih_64, \
                         lh_60, lh_61, lh_62, lh_63, lh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -6.0 * ih_60[k]
                  + f_0 * lh_60[k];

        t_61[k] = -6.0 * ih_61[k]
                  + f_0 * lh_61[k];

        t_62[k] = -6.0 * ih_62[k]
                  + f_0 * lh_62[k];

        t_63[k] = -5.0 * ih_63[k]
                  + f_0 * lh_63[k];

        t_64[k] = -5.0 * ih_64[k]
                  + f_0 * lh_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ih_65, ih_66, ih_67, ih_68, ih_69, \
                         lh_65, lh_66, lh_67, lh_68, lh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -5.0 * ih_65[k]
                  + f_0 * lh_65[k];

        t_66[k] = -5.0 * ih_66[k]
                  + f_0 * lh_66[k];

        t_67[k] = -5.0 * ih_67[k]
                  + f_0 * lh_67[k];

        t_68[k] = -5.0 * ih_68[k]
                  + f_0 * lh_68[k];

        t_69[k] = -5.0 * ih_69[k]
                  + f_0 * lh_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ih_70, ih_71, ih_72, ih_73, ih_74, \
                         lh_70, lh_71, lh_72, lh_73, lh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -5.0 * ih_70[k]
                  + f_0 * lh_70[k];

        t_71[k] = -5.0 * ih_71[k]
                  + f_0 * lh_71[k];

        t_72[k] = -5.0 * ih_72[k]
                  + f_0 * lh_72[k];

        t_73[k] = -5.0 * ih_73[k]
                  + f_0 * lh_73[k];

        t_74[k] = -5.0 * ih_74[k]
                  + f_0 * lh_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ih_75, ih_76, ih_77, ih_78, ih_79, \
                         lh_75, lh_76, lh_77, lh_78, lh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -5.0 * ih_75[k]
                  + f_0 * lh_75[k];

        t_76[k] = -5.0 * ih_76[k]
                  + f_0 * lh_76[k];

        t_77[k] = -5.0 * ih_77[k]
                  + f_0 * lh_77[k];

        t_78[k] = -5.0 * ih_78[k]
                  + f_0 * lh_78[k];

        t_79[k] = -5.0 * ih_79[k]
                  + f_0 * lh_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ih_80, ih_81, ih_82, ih_83, ih_84, \
                         lh_80, lh_81, lh_82, lh_83, lh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -5.0 * ih_80[k]
                  + f_0 * lh_80[k];

        t_81[k] = -5.0 * ih_81[k]
                  + f_0 * lh_81[k];

        t_82[k] = -5.0 * ih_82[k]
                  + f_0 * lh_82[k];

        t_83[k] = -5.0 * ih_83[k]
                  + f_0 * lh_83[k];

        t_84[k] = -5.0 * ih_84[k]
                  + f_0 * lh_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ih_85, ih_86, ih_87, ih_88, ih_89, \
                         lh_85, lh_86, lh_87, lh_88, lh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -5.0 * ih_85[k]
                  + f_0 * lh_85[k];

        t_86[k] = -5.0 * ih_86[k]
                  + f_0 * lh_86[k];

        t_87[k] = -5.0 * ih_87[k]
                  + f_0 * lh_87[k];

        t_88[k] = -5.0 * ih_88[k]
                  + f_0 * lh_88[k];

        t_89[k] = -5.0 * ih_89[k]
                  + f_0 * lh_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ih_90, ih_91, ih_92, ih_93, ih_94, \
                         lh_90, lh_91, lh_92, lh_93, lh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -5.0 * ih_90[k]
                  + f_0 * lh_90[k];

        t_91[k] = -5.0 * ih_91[k]
                  + f_0 * lh_91[k];

        t_92[k] = -5.0 * ih_92[k]
                  + f_0 * lh_92[k];

        t_93[k] = -5.0 * ih_93[k]
                  + f_0 * lh_93[k];

        t_94[k] = -5.0 * ih_94[k]
                  + f_0 * lh_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ih_95, ih_96, ih_97, ih_98, ih_99, \
                         lh_95, lh_96, lh_97, lh_98, lh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -5.0 * ih_95[k]
                  + f_0 * lh_95[k];

        t_96[k] = -5.0 * ih_96[k]
                  + f_0 * lh_96[k];

        t_97[k] = -5.0 * ih_97[k]
                  + f_0 * lh_97[k];

        t_98[k] = -5.0 * ih_98[k]
                  + f_0 * lh_98[k];

        t_99[k] = -5.0 * ih_99[k]
                  + f_0 * lh_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ih_100, ih_101, ih_102, ih_103, \
                         ih_104, lh_100, lh_101, lh_102, lh_103, \
                         lh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -5.0 * ih_100[k]
                   + f_0 * lh_100[k];

        t_101[k] = -5.0 * ih_101[k]
                   + f_0 * lh_101[k];

        t_102[k] = -5.0 * ih_102[k]
                   + f_0 * lh_102[k];

        t_103[k] = -5.0 * ih_103[k]
                   + f_0 * lh_103[k];

        t_104[k] = -5.0 * ih_104[k]
                   + f_0 * lh_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ih_105, ih_106, ih_107, ih_108, \
                         ih_109, lh_105, lh_106, lh_107, lh_108, \
                         lh_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -5.0 * ih_105[k]
                   + f_0 * lh_105[k];

        t_106[k] = -5.0 * ih_106[k]
                   + f_0 * lh_106[k];

        t_107[k] = -5.0 * ih_107[k]
                   + f_0 * lh_107[k];

        t_108[k] = -5.0 * ih_108[k]
                   + f_0 * lh_108[k];

        t_109[k] = -5.0 * ih_109[k]
                   + f_0 * lh_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ih_110, ih_111, ih_112, ih_113, \
                         ih_114, lh_110, lh_111, lh_112, lh_113, \
                         lh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -5.0 * ih_110[k]
                   + f_0 * lh_110[k];

        t_111[k] = -5.0 * ih_111[k]
                   + f_0 * lh_111[k];

        t_112[k] = -5.0 * ih_112[k]
                   + f_0 * lh_112[k];

        t_113[k] = -5.0 * ih_113[k]
                   + f_0 * lh_113[k];

        t_114[k] = -5.0 * ih_114[k]
                   + f_0 * lh_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ih_115, ih_116, ih_117, ih_118, \
                         ih_119, lh_115, lh_116, lh_117, lh_118, \
                         lh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -5.0 * ih_115[k]
                   + f_0 * lh_115[k];

        t_116[k] = -5.0 * ih_116[k]
                   + f_0 * lh_116[k];

        t_117[k] = -5.0 * ih_117[k]
                   + f_0 * lh_117[k];

        t_118[k] = -5.0 * ih_118[k]
                   + f_0 * lh_118[k];

        t_119[k] = -5.0 * ih_119[k]
                   + f_0 * lh_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ih_120, ih_121, ih_122, ih_123, \
                         ih_124, lh_120, lh_121, lh_122, lh_123, \
                         lh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -5.0 * ih_120[k]
                   + f_0 * lh_120[k];

        t_121[k] = -5.0 * ih_121[k]
                   + f_0 * lh_121[k];

        t_122[k] = -5.0 * ih_122[k]
                   + f_0 * lh_122[k];

        t_123[k] = -5.0 * ih_123[k]
                   + f_0 * lh_123[k];

        t_124[k] = -5.0 * ih_124[k]
                   + f_0 * lh_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ih_125, ih_126, ih_127, ih_128, \
                         ih_129, lh_125, lh_126, lh_127, lh_128, \
                         lh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -5.0 * ih_125[k]
                   + f_0 * lh_125[k];

        t_126[k] = -4.0 * ih_126[k]
                   + f_0 * lh_126[k];

        t_127[k] = -4.0 * ih_127[k]
                   + f_0 * lh_127[k];

        t_128[k] = -4.0 * ih_128[k]
                   + f_0 * lh_128[k];

        t_129[k] = -4.0 * ih_129[k]
                   + f_0 * lh_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ih_130, ih_131, ih_132, ih_133, \
                         ih_134, lh_130, lh_131, lh_132, lh_133, \
                         lh_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -4.0 * ih_130[k]
                   + f_0 * lh_130[k];

        t_131[k] = -4.0 * ih_131[k]
                   + f_0 * lh_131[k];

        t_132[k] = -4.0 * ih_132[k]
                   + f_0 * lh_132[k];

        t_133[k] = -4.0 * ih_133[k]
                   + f_0 * lh_133[k];

        t_134[k] = -4.0 * ih_134[k]
                   + f_0 * lh_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ih_135, ih_136, ih_137, ih_138, \
                         ih_139, lh_135, lh_136, lh_137, lh_138, \
                         lh_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -4.0 * ih_135[k]
                   + f_0 * lh_135[k];

        t_136[k] = -4.0 * ih_136[k]
                   + f_0 * lh_136[k];

        t_137[k] = -4.0 * ih_137[k]
                   + f_0 * lh_137[k];

        t_138[k] = -4.0 * ih_138[k]
                   + f_0 * lh_138[k];

        t_139[k] = -4.0 * ih_139[k]
                   + f_0 * lh_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ih_140, ih_141, ih_142, ih_143, \
                         ih_144, lh_140, lh_141, lh_142, lh_143, \
                         lh_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -4.0 * ih_140[k]
                   + f_0 * lh_140[k];

        t_141[k] = -4.0 * ih_141[k]
                   + f_0 * lh_141[k];

        t_142[k] = -4.0 * ih_142[k]
                   + f_0 * lh_142[k];

        t_143[k] = -4.0 * ih_143[k]
                   + f_0 * lh_143[k];

        t_144[k] = -4.0 * ih_144[k]
                   + f_0 * lh_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ih_145, ih_146, ih_147, ih_148, \
                         ih_149, lh_145, lh_146, lh_147, lh_148, \
                         lh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -4.0 * ih_145[k]
                   + f_0 * lh_145[k];

        t_146[k] = -4.0 * ih_146[k]
                   + f_0 * lh_146[k];

        t_147[k] = -4.0 * ih_147[k]
                   + f_0 * lh_147[k];

        t_148[k] = -4.0 * ih_148[k]
                   + f_0 * lh_148[k];

        t_149[k] = -4.0 * ih_149[k]
                   + f_0 * lh_149[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
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

    const auto *ih_150 = buffer.data(ih + 150);
    const auto *ih_151 = buffer.data(ih + 151);
    const auto *ih_152 = buffer.data(ih + 152);
    const auto *ih_153 = buffer.data(ih + 153);
    const auto *ih_154 = buffer.data(ih + 154);
    const auto *ih_155 = buffer.data(ih + 155);
    const auto *ih_156 = buffer.data(ih + 156);
    const auto *ih_157 = buffer.data(ih + 157);
    const auto *ih_158 = buffer.data(ih + 158);
    const auto *ih_159 = buffer.data(ih + 159);
    const auto *ih_160 = buffer.data(ih + 160);
    const auto *ih_161 = buffer.data(ih + 161);
    const auto *ih_162 = buffer.data(ih + 162);
    const auto *ih_163 = buffer.data(ih + 163);
    const auto *ih_164 = buffer.data(ih + 164);
    const auto *ih_165 = buffer.data(ih + 165);
    const auto *ih_166 = buffer.data(ih + 166);
    const auto *ih_167 = buffer.data(ih + 167);
    const auto *ih_168 = buffer.data(ih + 168);
    const auto *ih_169 = buffer.data(ih + 169);
    const auto *ih_170 = buffer.data(ih + 170);
    const auto *ih_171 = buffer.data(ih + 171);
    const auto *ih_172 = buffer.data(ih + 172);
    const auto *ih_173 = buffer.data(ih + 173);
    const auto *ih_174 = buffer.data(ih + 174);
    const auto *ih_175 = buffer.data(ih + 175);
    const auto *ih_176 = buffer.data(ih + 176);
    const auto *ih_177 = buffer.data(ih + 177);
    const auto *ih_178 = buffer.data(ih + 178);
    const auto *ih_179 = buffer.data(ih + 179);
    const auto *ih_180 = buffer.data(ih + 180);
    const auto *ih_181 = buffer.data(ih + 181);
    const auto *ih_182 = buffer.data(ih + 182);
    const auto *ih_183 = buffer.data(ih + 183);
    const auto *ih_184 = buffer.data(ih + 184);
    const auto *ih_185 = buffer.data(ih + 185);
    const auto *ih_186 = buffer.data(ih + 186);
    const auto *ih_187 = buffer.data(ih + 187);
    const auto *ih_188 = buffer.data(ih + 188);
    const auto *ih_189 = buffer.data(ih + 189);
    const auto *ih_190 = buffer.data(ih + 190);
    const auto *ih_191 = buffer.data(ih + 191);
    const auto *ih_192 = buffer.data(ih + 192);
    const auto *ih_193 = buffer.data(ih + 193);
    const auto *ih_194 = buffer.data(ih + 194);
    const auto *ih_195 = buffer.data(ih + 195);
    const auto *ih_196 = buffer.data(ih + 196);
    const auto *ih_197 = buffer.data(ih + 197);
    const auto *ih_198 = buffer.data(ih + 198);
    const auto *ih_199 = buffer.data(ih + 199);
    const auto *ih_200 = buffer.data(ih + 200);
    const auto *ih_201 = buffer.data(ih + 201);
    const auto *ih_202 = buffer.data(ih + 202);
    const auto *ih_203 = buffer.data(ih + 203);
    const auto *ih_204 = buffer.data(ih + 204);
    const auto *ih_205 = buffer.data(ih + 205);
    const auto *ih_206 = buffer.data(ih + 206);
    const auto *ih_207 = buffer.data(ih + 207);
    const auto *ih_208 = buffer.data(ih + 208);
    const auto *ih_209 = buffer.data(ih + 209);
    const auto *ih_210 = buffer.data(ih + 210);
    const auto *ih_211 = buffer.data(ih + 211);
    const auto *ih_212 = buffer.data(ih + 212);
    const auto *ih_213 = buffer.data(ih + 213);
    const auto *ih_214 = buffer.data(ih + 214);
    const auto *ih_215 = buffer.data(ih + 215);
    const auto *ih_216 = buffer.data(ih + 216);
    const auto *ih_217 = buffer.data(ih + 217);
    const auto *ih_218 = buffer.data(ih + 218);
    const auto *ih_219 = buffer.data(ih + 219);
    const auto *ih_220 = buffer.data(ih + 220);
    const auto *ih_221 = buffer.data(ih + 221);
    const auto *ih_222 = buffer.data(ih + 222);
    const auto *ih_223 = buffer.data(ih + 223);
    const auto *ih_224 = buffer.data(ih + 224);
    const auto *ih_225 = buffer.data(ih + 225);
    const auto *ih_226 = buffer.data(ih + 226);
    const auto *ih_227 = buffer.data(ih + 227);
    const auto *ih_228 = buffer.data(ih + 228);
    const auto *ih_229 = buffer.data(ih + 229);
    const auto *ih_230 = buffer.data(ih + 230);
    const auto *ih_231 = buffer.data(ih + 231);
    const auto *ih_232 = buffer.data(ih + 232);
    const auto *ih_233 = buffer.data(ih + 233);
    const auto *ih_234 = buffer.data(ih + 234);
    const auto *ih_235 = buffer.data(ih + 235);
    const auto *ih_236 = buffer.data(ih + 236);
    const auto *ih_237 = buffer.data(ih + 237);
    const auto *ih_238 = buffer.data(ih + 238);
    const auto *ih_239 = buffer.data(ih + 239);
    const auto *ih_240 = buffer.data(ih + 240);
    const auto *ih_241 = buffer.data(ih + 241);
    const auto *ih_242 = buffer.data(ih + 242);
    const auto *ih_243 = buffer.data(ih + 243);
    const auto *ih_244 = buffer.data(ih + 244);
    const auto *ih_245 = buffer.data(ih + 245);
    const auto *ih_246 = buffer.data(ih + 246);
    const auto *ih_247 = buffer.data(ih + 247);
    const auto *ih_248 = buffer.data(ih + 248);
    const auto *ih_249 = buffer.data(ih + 249);
    const auto *ih_250 = buffer.data(ih + 250);
    const auto *ih_251 = buffer.data(ih + 251);
    const auto *ih_252 = buffer.data(ih + 252);
    const auto *ih_253 = buffer.data(ih + 253);
    const auto *ih_254 = buffer.data(ih + 254);
    const auto *ih_255 = buffer.data(ih + 255);
    const auto *ih_256 = buffer.data(ih + 256);
    const auto *ih_257 = buffer.data(ih + 257);
    const auto *ih_258 = buffer.data(ih + 258);
    const auto *ih_259 = buffer.data(ih + 259);
    const auto *ih_260 = buffer.data(ih + 260);
    const auto *ih_261 = buffer.data(ih + 261);
    const auto *ih_262 = buffer.data(ih + 262);
    const auto *ih_263 = buffer.data(ih + 263);
    const auto *ih_264 = buffer.data(ih + 264);
    const auto *ih_265 = buffer.data(ih + 265);
    const auto *ih_266 = buffer.data(ih + 266);
    const auto *ih_267 = buffer.data(ih + 267);
    const auto *ih_268 = buffer.data(ih + 268);
    const auto *ih_269 = buffer.data(ih + 269);
    const auto *ih_270 = buffer.data(ih + 270);
    const auto *ih_271 = buffer.data(ih + 271);
    const auto *ih_272 = buffer.data(ih + 272);
    const auto *ih_273 = buffer.data(ih + 273);
    const auto *ih_274 = buffer.data(ih + 274);
    const auto *ih_275 = buffer.data(ih + 275);
    const auto *ih_276 = buffer.data(ih + 276);
    const auto *ih_277 = buffer.data(ih + 277);
    const auto *ih_278 = buffer.data(ih + 278);
    const auto *ih_279 = buffer.data(ih + 279);
    const auto *ih_280 = buffer.data(ih + 280);
    const auto *ih_281 = buffer.data(ih + 281);
    const auto *ih_282 = buffer.data(ih + 282);
    const auto *ih_283 = buffer.data(ih + 283);
    const auto *ih_284 = buffer.data(ih + 284);
    const auto *ih_285 = buffer.data(ih + 285);
    const auto *ih_286 = buffer.data(ih + 286);
    const auto *ih_287 = buffer.data(ih + 287);
    const auto *ih_288 = buffer.data(ih + 288);
    const auto *ih_289 = buffer.data(ih + 289);
    const auto *ih_290 = buffer.data(ih + 290);
    const auto *ih_291 = buffer.data(ih + 291);
    const auto *ih_292 = buffer.data(ih + 292);
    const auto *ih_293 = buffer.data(ih + 293);
    const auto *ih_294 = buffer.data(ih + 294);
    const auto *ih_295 = buffer.data(ih + 295);
    const auto *ih_296 = buffer.data(ih + 296);
    const auto *ih_297 = buffer.data(ih + 297);
    const auto *ih_298 = buffer.data(ih + 298);
    const auto *ih_299 = buffer.data(ih + 299);

    const auto *lh_150 = buffer.data(lh + 150);
    const auto *lh_151 = buffer.data(lh + 151);
    const auto *lh_152 = buffer.data(lh + 152);
    const auto *lh_153 = buffer.data(lh + 153);
    const auto *lh_154 = buffer.data(lh + 154);
    const auto *lh_155 = buffer.data(lh + 155);
    const auto *lh_156 = buffer.data(lh + 156);
    const auto *lh_157 = buffer.data(lh + 157);
    const auto *lh_158 = buffer.data(lh + 158);
    const auto *lh_159 = buffer.data(lh + 159);
    const auto *lh_160 = buffer.data(lh + 160);
    const auto *lh_161 = buffer.data(lh + 161);
    const auto *lh_162 = buffer.data(lh + 162);
    const auto *lh_163 = buffer.data(lh + 163);
    const auto *lh_164 = buffer.data(lh + 164);
    const auto *lh_165 = buffer.data(lh + 165);
    const auto *lh_166 = buffer.data(lh + 166);
    const auto *lh_167 = buffer.data(lh + 167);
    const auto *lh_168 = buffer.data(lh + 168);
    const auto *lh_169 = buffer.data(lh + 169);
    const auto *lh_170 = buffer.data(lh + 170);
    const auto *lh_171 = buffer.data(lh + 171);
    const auto *lh_172 = buffer.data(lh + 172);
    const auto *lh_173 = buffer.data(lh + 173);
    const auto *lh_174 = buffer.data(lh + 174);
    const auto *lh_175 = buffer.data(lh + 175);
    const auto *lh_176 = buffer.data(lh + 176);
    const auto *lh_177 = buffer.data(lh + 177);
    const auto *lh_178 = buffer.data(lh + 178);
    const auto *lh_179 = buffer.data(lh + 179);
    const auto *lh_180 = buffer.data(lh + 180);
    const auto *lh_181 = buffer.data(lh + 181);
    const auto *lh_182 = buffer.data(lh + 182);
    const auto *lh_183 = buffer.data(lh + 183);
    const auto *lh_184 = buffer.data(lh + 184);
    const auto *lh_185 = buffer.data(lh + 185);
    const auto *lh_186 = buffer.data(lh + 186);
    const auto *lh_187 = buffer.data(lh + 187);
    const auto *lh_188 = buffer.data(lh + 188);
    const auto *lh_189 = buffer.data(lh + 189);
    const auto *lh_190 = buffer.data(lh + 190);
    const auto *lh_191 = buffer.data(lh + 191);
    const auto *lh_192 = buffer.data(lh + 192);
    const auto *lh_193 = buffer.data(lh + 193);
    const auto *lh_194 = buffer.data(lh + 194);
    const auto *lh_195 = buffer.data(lh + 195);
    const auto *lh_196 = buffer.data(lh + 196);
    const auto *lh_197 = buffer.data(lh + 197);
    const auto *lh_198 = buffer.data(lh + 198);
    const auto *lh_199 = buffer.data(lh + 199);
    const auto *lh_200 = buffer.data(lh + 200);
    const auto *lh_201 = buffer.data(lh + 201);
    const auto *lh_202 = buffer.data(lh + 202);
    const auto *lh_203 = buffer.data(lh + 203);
    const auto *lh_204 = buffer.data(lh + 204);
    const auto *lh_205 = buffer.data(lh + 205);
    const auto *lh_206 = buffer.data(lh + 206);
    const auto *lh_207 = buffer.data(lh + 207);
    const auto *lh_208 = buffer.data(lh + 208);
    const auto *lh_209 = buffer.data(lh + 209);
    const auto *lh_210 = buffer.data(lh + 210);
    const auto *lh_211 = buffer.data(lh + 211);
    const auto *lh_212 = buffer.data(lh + 212);
    const auto *lh_213 = buffer.data(lh + 213);
    const auto *lh_214 = buffer.data(lh + 214);
    const auto *lh_215 = buffer.data(lh + 215);
    const auto *lh_216 = buffer.data(lh + 216);
    const auto *lh_217 = buffer.data(lh + 217);
    const auto *lh_218 = buffer.data(lh + 218);
    const auto *lh_219 = buffer.data(lh + 219);
    const auto *lh_220 = buffer.data(lh + 220);
    const auto *lh_221 = buffer.data(lh + 221);
    const auto *lh_222 = buffer.data(lh + 222);
    const auto *lh_223 = buffer.data(lh + 223);
    const auto *lh_224 = buffer.data(lh + 224);
    const auto *lh_225 = buffer.data(lh + 225);
    const auto *lh_226 = buffer.data(lh + 226);
    const auto *lh_227 = buffer.data(lh + 227);
    const auto *lh_228 = buffer.data(lh + 228);
    const auto *lh_229 = buffer.data(lh + 229);
    const auto *lh_230 = buffer.data(lh + 230);
    const auto *lh_231 = buffer.data(lh + 231);
    const auto *lh_232 = buffer.data(lh + 232);
    const auto *lh_233 = buffer.data(lh + 233);
    const auto *lh_234 = buffer.data(lh + 234);
    const auto *lh_235 = buffer.data(lh + 235);
    const auto *lh_236 = buffer.data(lh + 236);
    const auto *lh_237 = buffer.data(lh + 237);
    const auto *lh_238 = buffer.data(lh + 238);
    const auto *lh_239 = buffer.data(lh + 239);
    const auto *lh_240 = buffer.data(lh + 240);
    const auto *lh_241 = buffer.data(lh + 241);
    const auto *lh_242 = buffer.data(lh + 242);
    const auto *lh_243 = buffer.data(lh + 243);
    const auto *lh_244 = buffer.data(lh + 244);
    const auto *lh_245 = buffer.data(lh + 245);
    const auto *lh_246 = buffer.data(lh + 246);
    const auto *lh_247 = buffer.data(lh + 247);
    const auto *lh_248 = buffer.data(lh + 248);
    const auto *lh_249 = buffer.data(lh + 249);
    const auto *lh_250 = buffer.data(lh + 250);
    const auto *lh_251 = buffer.data(lh + 251);
    const auto *lh_252 = buffer.data(lh + 252);
    const auto *lh_253 = buffer.data(lh + 253);
    const auto *lh_254 = buffer.data(lh + 254);
    const auto *lh_255 = buffer.data(lh + 255);
    const auto *lh_256 = buffer.data(lh + 256);
    const auto *lh_257 = buffer.data(lh + 257);
    const auto *lh_258 = buffer.data(lh + 258);
    const auto *lh_259 = buffer.data(lh + 259);
    const auto *lh_260 = buffer.data(lh + 260);
    const auto *lh_261 = buffer.data(lh + 261);
    const auto *lh_262 = buffer.data(lh + 262);
    const auto *lh_263 = buffer.data(lh + 263);
    const auto *lh_264 = buffer.data(lh + 264);
    const auto *lh_265 = buffer.data(lh + 265);
    const auto *lh_266 = buffer.data(lh + 266);
    const auto *lh_267 = buffer.data(lh + 267);
    const auto *lh_268 = buffer.data(lh + 268);
    const auto *lh_269 = buffer.data(lh + 269);
    const auto *lh_270 = buffer.data(lh + 270);
    const auto *lh_271 = buffer.data(lh + 271);
    const auto *lh_272 = buffer.data(lh + 272);
    const auto *lh_273 = buffer.data(lh + 273);
    const auto *lh_274 = buffer.data(lh + 274);
    const auto *lh_275 = buffer.data(lh + 275);
    const auto *lh_276 = buffer.data(lh + 276);
    const auto *lh_277 = buffer.data(lh + 277);
    const auto *lh_278 = buffer.data(lh + 278);
    const auto *lh_279 = buffer.data(lh + 279);
    const auto *lh_280 = buffer.data(lh + 280);
    const auto *lh_281 = buffer.data(lh + 281);
    const auto *lh_282 = buffer.data(lh + 282);
    const auto *lh_283 = buffer.data(lh + 283);
    const auto *lh_284 = buffer.data(lh + 284);
    const auto *lh_285 = buffer.data(lh + 285);
    const auto *lh_286 = buffer.data(lh + 286);
    const auto *lh_287 = buffer.data(lh + 287);
    const auto *lh_288 = buffer.data(lh + 288);
    const auto *lh_289 = buffer.data(lh + 289);
    const auto *lh_290 = buffer.data(lh + 290);
    const auto *lh_291 = buffer.data(lh + 291);
    const auto *lh_292 = buffer.data(lh + 292);
    const auto *lh_293 = buffer.data(lh + 293);
    const auto *lh_294 = buffer.data(lh + 294);
    const auto *lh_295 = buffer.data(lh + 295);
    const auto *lh_296 = buffer.data(lh + 296);
    const auto *lh_297 = buffer.data(lh + 297);
    const auto *lh_298 = buffer.data(lh + 298);
    const auto *lh_299 = buffer.data(lh + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ih_150, ih_151, ih_152, ih_153, \
                         ih_154, lh_150, lh_151, lh_152, lh_153, \
                         lh_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -4.0 * ih_150[k]
                   + f_0 * lh_150[k];

        t_151[k] = -4.0 * ih_151[k]
                   + f_0 * lh_151[k];

        t_152[k] = -4.0 * ih_152[k]
                   + f_0 * lh_152[k];

        t_153[k] = -4.0 * ih_153[k]
                   + f_0 * lh_153[k];

        t_154[k] = -4.0 * ih_154[k]
                   + f_0 * lh_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ih_155, ih_156, ih_157, ih_158, \
                         ih_159, lh_155, lh_156, lh_157, lh_158, \
                         lh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -4.0 * ih_155[k]
                   + f_0 * lh_155[k];

        t_156[k] = -4.0 * ih_156[k]
                   + f_0 * lh_156[k];

        t_157[k] = -4.0 * ih_157[k]
                   + f_0 * lh_157[k];

        t_158[k] = -4.0 * ih_158[k]
                   + f_0 * lh_158[k];

        t_159[k] = -4.0 * ih_159[k]
                   + f_0 * lh_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ih_160, ih_161, ih_162, ih_163, \
                         ih_164, lh_160, lh_161, lh_162, lh_163, \
                         lh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -4.0 * ih_160[k]
                   + f_0 * lh_160[k];

        t_161[k] = -4.0 * ih_161[k]
                   + f_0 * lh_161[k];

        t_162[k] = -4.0 * ih_162[k]
                   + f_0 * lh_162[k];

        t_163[k] = -4.0 * ih_163[k]
                   + f_0 * lh_163[k];

        t_164[k] = -4.0 * ih_164[k]
                   + f_0 * lh_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ih_165, ih_166, ih_167, ih_168, \
                         ih_169, lh_165, lh_166, lh_167, lh_168, \
                         lh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -4.0 * ih_165[k]
                   + f_0 * lh_165[k];

        t_166[k] = -4.0 * ih_166[k]
                   + f_0 * lh_166[k];

        t_167[k] = -4.0 * ih_167[k]
                   + f_0 * lh_167[k];

        t_168[k] = -4.0 * ih_168[k]
                   + f_0 * lh_168[k];

        t_169[k] = -4.0 * ih_169[k]
                   + f_0 * lh_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ih_170, ih_171, ih_172, ih_173, \
                         ih_174, lh_170, lh_171, lh_172, lh_173, \
                         lh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -4.0 * ih_170[k]
                   + f_0 * lh_170[k];

        t_171[k] = -4.0 * ih_171[k]
                   + f_0 * lh_171[k];

        t_172[k] = -4.0 * ih_172[k]
                   + f_0 * lh_172[k];

        t_173[k] = -4.0 * ih_173[k]
                   + f_0 * lh_173[k];

        t_174[k] = -4.0 * ih_174[k]
                   + f_0 * lh_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ih_175, ih_176, ih_177, ih_178, \
                         ih_179, lh_175, lh_176, lh_177, lh_178, \
                         lh_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -4.0 * ih_175[k]
                   + f_0 * lh_175[k];

        t_176[k] = -4.0 * ih_176[k]
                   + f_0 * lh_176[k];

        t_177[k] = -4.0 * ih_177[k]
                   + f_0 * lh_177[k];

        t_178[k] = -4.0 * ih_178[k]
                   + f_0 * lh_178[k];

        t_179[k] = -4.0 * ih_179[k]
                   + f_0 * lh_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ih_180, ih_181, ih_182, ih_183, \
                         ih_184, lh_180, lh_181, lh_182, lh_183, \
                         lh_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -4.0 * ih_180[k]
                   + f_0 * lh_180[k];

        t_181[k] = -4.0 * ih_181[k]
                   + f_0 * lh_181[k];

        t_182[k] = -4.0 * ih_182[k]
                   + f_0 * lh_182[k];

        t_183[k] = -4.0 * ih_183[k]
                   + f_0 * lh_183[k];

        t_184[k] = -4.0 * ih_184[k]
                   + f_0 * lh_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ih_185, ih_186, ih_187, ih_188, \
                         ih_189, lh_185, lh_186, lh_187, lh_188, \
                         lh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -4.0 * ih_185[k]
                   + f_0 * lh_185[k];

        t_186[k] = -4.0 * ih_186[k]
                   + f_0 * lh_186[k];

        t_187[k] = -4.0 * ih_187[k]
                   + f_0 * lh_187[k];

        t_188[k] = -4.0 * ih_188[k]
                   + f_0 * lh_188[k];

        t_189[k] = -4.0 * ih_189[k]
                   + f_0 * lh_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ih_190, ih_191, ih_192, ih_193, \
                         ih_194, lh_190, lh_191, lh_192, lh_193, \
                         lh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -4.0 * ih_190[k]
                   + f_0 * lh_190[k];

        t_191[k] = -4.0 * ih_191[k]
                   + f_0 * lh_191[k];

        t_192[k] = -4.0 * ih_192[k]
                   + f_0 * lh_192[k];

        t_193[k] = -4.0 * ih_193[k]
                   + f_0 * lh_193[k];

        t_194[k] = -4.0 * ih_194[k]
                   + f_0 * lh_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ih_195, ih_196, ih_197, ih_198, \
                         ih_199, lh_195, lh_196, lh_197, lh_198, \
                         lh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -4.0 * ih_195[k]
                   + f_0 * lh_195[k];

        t_196[k] = -4.0 * ih_196[k]
                   + f_0 * lh_196[k];

        t_197[k] = -4.0 * ih_197[k]
                   + f_0 * lh_197[k];

        t_198[k] = -4.0 * ih_198[k]
                   + f_0 * lh_198[k];

        t_199[k] = -4.0 * ih_199[k]
                   + f_0 * lh_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ih_200, ih_201, ih_202, ih_203, \
                         ih_204, lh_200, lh_201, lh_202, lh_203, \
                         lh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -4.0 * ih_200[k]
                   + f_0 * lh_200[k];

        t_201[k] = -4.0 * ih_201[k]
                   + f_0 * lh_201[k];

        t_202[k] = -4.0 * ih_202[k]
                   + f_0 * lh_202[k];

        t_203[k] = -4.0 * ih_203[k]
                   + f_0 * lh_203[k];

        t_204[k] = -4.0 * ih_204[k]
                   + f_0 * lh_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ih_205, ih_206, ih_207, ih_208, \
                         ih_209, lh_205, lh_206, lh_207, lh_208, \
                         lh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -4.0 * ih_205[k]
                   + f_0 * lh_205[k];

        t_206[k] = -4.0 * ih_206[k]
                   + f_0 * lh_206[k];

        t_207[k] = -4.0 * ih_207[k]
                   + f_0 * lh_207[k];

        t_208[k] = -4.0 * ih_208[k]
                   + f_0 * lh_208[k];

        t_209[k] = -4.0 * ih_209[k]
                   + f_0 * lh_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ih_210, ih_211, ih_212, ih_213, \
                         ih_214, lh_210, lh_211, lh_212, lh_213, \
                         lh_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -3.0 * ih_210[k]
                   + f_0 * lh_210[k];

        t_211[k] = -3.0 * ih_211[k]
                   + f_0 * lh_211[k];

        t_212[k] = -3.0 * ih_212[k]
                   + f_0 * lh_212[k];

        t_213[k] = -3.0 * ih_213[k]
                   + f_0 * lh_213[k];

        t_214[k] = -3.0 * ih_214[k]
                   + f_0 * lh_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ih_215, ih_216, ih_217, ih_218, \
                         ih_219, lh_215, lh_216, lh_217, lh_218, \
                         lh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -3.0 * ih_215[k]
                   + f_0 * lh_215[k];

        t_216[k] = -3.0 * ih_216[k]
                   + f_0 * lh_216[k];

        t_217[k] = -3.0 * ih_217[k]
                   + f_0 * lh_217[k];

        t_218[k] = -3.0 * ih_218[k]
                   + f_0 * lh_218[k];

        t_219[k] = -3.0 * ih_219[k]
                   + f_0 * lh_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ih_220, ih_221, ih_222, ih_223, \
                         ih_224, lh_220, lh_221, lh_222, lh_223, \
                         lh_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -3.0 * ih_220[k]
                   + f_0 * lh_220[k];

        t_221[k] = -3.0 * ih_221[k]
                   + f_0 * lh_221[k];

        t_222[k] = -3.0 * ih_222[k]
                   + f_0 * lh_222[k];

        t_223[k] = -3.0 * ih_223[k]
                   + f_0 * lh_223[k];

        t_224[k] = -3.0 * ih_224[k]
                   + f_0 * lh_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ih_225, ih_226, ih_227, ih_228, \
                         ih_229, lh_225, lh_226, lh_227, lh_228, \
                         lh_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -3.0 * ih_225[k]
                   + f_0 * lh_225[k];

        t_226[k] = -3.0 * ih_226[k]
                   + f_0 * lh_226[k];

        t_227[k] = -3.0 * ih_227[k]
                   + f_0 * lh_227[k];

        t_228[k] = -3.0 * ih_228[k]
                   + f_0 * lh_228[k];

        t_229[k] = -3.0 * ih_229[k]
                   + f_0 * lh_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ih_230, ih_231, ih_232, ih_233, \
                         ih_234, lh_230, lh_231, lh_232, lh_233, \
                         lh_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -3.0 * ih_230[k]
                   + f_0 * lh_230[k];

        t_231[k] = -3.0 * ih_231[k]
                   + f_0 * lh_231[k];

        t_232[k] = -3.0 * ih_232[k]
                   + f_0 * lh_232[k];

        t_233[k] = -3.0 * ih_233[k]
                   + f_0 * lh_233[k];

        t_234[k] = -3.0 * ih_234[k]
                   + f_0 * lh_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ih_235, ih_236, ih_237, ih_238, \
                         ih_239, lh_235, lh_236, lh_237, lh_238, \
                         lh_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -3.0 * ih_235[k]
                   + f_0 * lh_235[k];

        t_236[k] = -3.0 * ih_236[k]
                   + f_0 * lh_236[k];

        t_237[k] = -3.0 * ih_237[k]
                   + f_0 * lh_237[k];

        t_238[k] = -3.0 * ih_238[k]
                   + f_0 * lh_238[k];

        t_239[k] = -3.0 * ih_239[k]
                   + f_0 * lh_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ih_240, ih_241, ih_242, ih_243, \
                         ih_244, lh_240, lh_241, lh_242, lh_243, \
                         lh_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * ih_240[k]
                   + f_0 * lh_240[k];

        t_241[k] = -3.0 * ih_241[k]
                   + f_0 * lh_241[k];

        t_242[k] = -3.0 * ih_242[k]
                   + f_0 * lh_242[k];

        t_243[k] = -3.0 * ih_243[k]
                   + f_0 * lh_243[k];

        t_244[k] = -3.0 * ih_244[k]
                   + f_0 * lh_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ih_245, ih_246, ih_247, ih_248, \
                         ih_249, lh_245, lh_246, lh_247, lh_248, \
                         lh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * ih_245[k]
                   + f_0 * lh_245[k];

        t_246[k] = -3.0 * ih_246[k]
                   + f_0 * lh_246[k];

        t_247[k] = -3.0 * ih_247[k]
                   + f_0 * lh_247[k];

        t_248[k] = -3.0 * ih_248[k]
                   + f_0 * lh_248[k];

        t_249[k] = -3.0 * ih_249[k]
                   + f_0 * lh_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ih_250, ih_251, ih_252, ih_253, \
                         ih_254, lh_250, lh_251, lh_252, lh_253, \
                         lh_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -3.0 * ih_250[k]
                   + f_0 * lh_250[k];

        t_251[k] = -3.0 * ih_251[k]
                   + f_0 * lh_251[k];

        t_252[k] = -3.0 * ih_252[k]
                   + f_0 * lh_252[k];

        t_253[k] = -3.0 * ih_253[k]
                   + f_0 * lh_253[k];

        t_254[k] = -3.0 * ih_254[k]
                   + f_0 * lh_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ih_255, ih_256, ih_257, ih_258, \
                         ih_259, lh_255, lh_256, lh_257, lh_258, \
                         lh_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -3.0 * ih_255[k]
                   + f_0 * lh_255[k];

        t_256[k] = -3.0 * ih_256[k]
                   + f_0 * lh_256[k];

        t_257[k] = -3.0 * ih_257[k]
                   + f_0 * lh_257[k];

        t_258[k] = -3.0 * ih_258[k]
                   + f_0 * lh_258[k];

        t_259[k] = -3.0 * ih_259[k]
                   + f_0 * lh_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ih_260, ih_261, ih_262, ih_263, \
                         ih_264, lh_260, lh_261, lh_262, lh_263, \
                         lh_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -3.0 * ih_260[k]
                   + f_0 * lh_260[k];

        t_261[k] = -3.0 * ih_261[k]
                   + f_0 * lh_261[k];

        t_262[k] = -3.0 * ih_262[k]
                   + f_0 * lh_262[k];

        t_263[k] = -3.0 * ih_263[k]
                   + f_0 * lh_263[k];

        t_264[k] = -3.0 * ih_264[k]
                   + f_0 * lh_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ih_265, ih_266, ih_267, ih_268, \
                         ih_269, lh_265, lh_266, lh_267, lh_268, \
                         lh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -3.0 * ih_265[k]
                   + f_0 * lh_265[k];

        t_266[k] = -3.0 * ih_266[k]
                   + f_0 * lh_266[k];

        t_267[k] = -3.0 * ih_267[k]
                   + f_0 * lh_267[k];

        t_268[k] = -3.0 * ih_268[k]
                   + f_0 * lh_268[k];

        t_269[k] = -3.0 * ih_269[k]
                   + f_0 * lh_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ih_270, ih_271, ih_272, ih_273, \
                         ih_274, lh_270, lh_271, lh_272, lh_273, \
                         lh_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -3.0 * ih_270[k]
                   + f_0 * lh_270[k];

        t_271[k] = -3.0 * ih_271[k]
                   + f_0 * lh_271[k];

        t_272[k] = -3.0 * ih_272[k]
                   + f_0 * lh_272[k];

        t_273[k] = -3.0 * ih_273[k]
                   + f_0 * lh_273[k];

        t_274[k] = -3.0 * ih_274[k]
                   + f_0 * lh_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ih_275, ih_276, ih_277, ih_278, \
                         ih_279, lh_275, lh_276, lh_277, lh_278, \
                         lh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * ih_275[k]
                   + f_0 * lh_275[k];

        t_276[k] = -3.0 * ih_276[k]
                   + f_0 * lh_276[k];

        t_277[k] = -3.0 * ih_277[k]
                   + f_0 * lh_277[k];

        t_278[k] = -3.0 * ih_278[k]
                   + f_0 * lh_278[k];

        t_279[k] = -3.0 * ih_279[k]
                   + f_0 * lh_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ih_280, ih_281, ih_282, ih_283, \
                         ih_284, lh_280, lh_281, lh_282, lh_283, \
                         lh_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * ih_280[k]
                   + f_0 * lh_280[k];

        t_281[k] = -3.0 * ih_281[k]
                   + f_0 * lh_281[k];

        t_282[k] = -3.0 * ih_282[k]
                   + f_0 * lh_282[k];

        t_283[k] = -3.0 * ih_283[k]
                   + f_0 * lh_283[k];

        t_284[k] = -3.0 * ih_284[k]
                   + f_0 * lh_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ih_285, ih_286, ih_287, ih_288, \
                         ih_289, lh_285, lh_286, lh_287, lh_288, \
                         lh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * ih_285[k]
                   + f_0 * lh_285[k];

        t_286[k] = -3.0 * ih_286[k]
                   + f_0 * lh_286[k];

        t_287[k] = -3.0 * ih_287[k]
                   + f_0 * lh_287[k];

        t_288[k] = -3.0 * ih_288[k]
                   + f_0 * lh_288[k];

        t_289[k] = -3.0 * ih_289[k]
                   + f_0 * lh_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ih_290, ih_291, ih_292, ih_293, \
                         ih_294, lh_290, lh_291, lh_292, lh_293, \
                         lh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * ih_290[k]
                   + f_0 * lh_290[k];

        t_291[k] = -3.0 * ih_291[k]
                   + f_0 * lh_291[k];

        t_292[k] = -3.0 * ih_292[k]
                   + f_0 * lh_292[k];

        t_293[k] = -3.0 * ih_293[k]
                   + f_0 * lh_293[k];

        t_294[k] = -3.0 * ih_294[k]
                   + f_0 * lh_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ih_295, ih_296, ih_297, ih_298, \
                         ih_299, lh_295, lh_296, lh_297, lh_298, \
                         lh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -3.0 * ih_295[k]
                   + f_0 * lh_295[k];

        t_296[k] = -3.0 * ih_296[k]
                   + f_0 * lh_296[k];

        t_297[k] = -3.0 * ih_297[k]
                   + f_0 * lh_297[k];

        t_298[k] = -3.0 * ih_298[k]
                   + f_0 * lh_298[k];

        t_299[k] = -3.0 * ih_299[k]
                   + f_0 * lh_299[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
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

    const auto *ih_300 = buffer.data(ih + 300);
    const auto *ih_301 = buffer.data(ih + 301);
    const auto *ih_302 = buffer.data(ih + 302);
    const auto *ih_303 = buffer.data(ih + 303);
    const auto *ih_304 = buffer.data(ih + 304);
    const auto *ih_305 = buffer.data(ih + 305);
    const auto *ih_306 = buffer.data(ih + 306);
    const auto *ih_307 = buffer.data(ih + 307);
    const auto *ih_308 = buffer.data(ih + 308);
    const auto *ih_309 = buffer.data(ih + 309);
    const auto *ih_310 = buffer.data(ih + 310);
    const auto *ih_311 = buffer.data(ih + 311);
    const auto *ih_312 = buffer.data(ih + 312);
    const auto *ih_313 = buffer.data(ih + 313);
    const auto *ih_314 = buffer.data(ih + 314);
    const auto *ih_315 = buffer.data(ih + 315);
    const auto *ih_316 = buffer.data(ih + 316);
    const auto *ih_317 = buffer.data(ih + 317);
    const auto *ih_318 = buffer.data(ih + 318);
    const auto *ih_319 = buffer.data(ih + 319);
    const auto *ih_320 = buffer.data(ih + 320);
    const auto *ih_321 = buffer.data(ih + 321);
    const auto *ih_322 = buffer.data(ih + 322);
    const auto *ih_323 = buffer.data(ih + 323);
    const auto *ih_324 = buffer.data(ih + 324);
    const auto *ih_325 = buffer.data(ih + 325);
    const auto *ih_326 = buffer.data(ih + 326);
    const auto *ih_327 = buffer.data(ih + 327);
    const auto *ih_328 = buffer.data(ih + 328);
    const auto *ih_329 = buffer.data(ih + 329);
    const auto *ih_330 = buffer.data(ih + 330);
    const auto *ih_331 = buffer.data(ih + 331);
    const auto *ih_332 = buffer.data(ih + 332);
    const auto *ih_333 = buffer.data(ih + 333);
    const auto *ih_334 = buffer.data(ih + 334);
    const auto *ih_335 = buffer.data(ih + 335);
    const auto *ih_336 = buffer.data(ih + 336);
    const auto *ih_337 = buffer.data(ih + 337);
    const auto *ih_338 = buffer.data(ih + 338);
    const auto *ih_339 = buffer.data(ih + 339);
    const auto *ih_340 = buffer.data(ih + 340);
    const auto *ih_341 = buffer.data(ih + 341);
    const auto *ih_342 = buffer.data(ih + 342);
    const auto *ih_343 = buffer.data(ih + 343);
    const auto *ih_344 = buffer.data(ih + 344);
    const auto *ih_345 = buffer.data(ih + 345);
    const auto *ih_346 = buffer.data(ih + 346);
    const auto *ih_347 = buffer.data(ih + 347);
    const auto *ih_348 = buffer.data(ih + 348);
    const auto *ih_349 = buffer.data(ih + 349);
    const auto *ih_350 = buffer.data(ih + 350);
    const auto *ih_351 = buffer.data(ih + 351);
    const auto *ih_352 = buffer.data(ih + 352);
    const auto *ih_353 = buffer.data(ih + 353);
    const auto *ih_354 = buffer.data(ih + 354);
    const auto *ih_355 = buffer.data(ih + 355);
    const auto *ih_356 = buffer.data(ih + 356);
    const auto *ih_357 = buffer.data(ih + 357);
    const auto *ih_358 = buffer.data(ih + 358);
    const auto *ih_359 = buffer.data(ih + 359);
    const auto *ih_360 = buffer.data(ih + 360);
    const auto *ih_361 = buffer.data(ih + 361);
    const auto *ih_362 = buffer.data(ih + 362);
    const auto *ih_363 = buffer.data(ih + 363);
    const auto *ih_364 = buffer.data(ih + 364);
    const auto *ih_365 = buffer.data(ih + 365);
    const auto *ih_366 = buffer.data(ih + 366);
    const auto *ih_367 = buffer.data(ih + 367);
    const auto *ih_368 = buffer.data(ih + 368);
    const auto *ih_369 = buffer.data(ih + 369);
    const auto *ih_370 = buffer.data(ih + 370);
    const auto *ih_371 = buffer.data(ih + 371);
    const auto *ih_372 = buffer.data(ih + 372);
    const auto *ih_373 = buffer.data(ih + 373);
    const auto *ih_374 = buffer.data(ih + 374);
    const auto *ih_375 = buffer.data(ih + 375);
    const auto *ih_376 = buffer.data(ih + 376);
    const auto *ih_377 = buffer.data(ih + 377);
    const auto *ih_378 = buffer.data(ih + 378);
    const auto *ih_379 = buffer.data(ih + 379);
    const auto *ih_380 = buffer.data(ih + 380);
    const auto *ih_381 = buffer.data(ih + 381);
    const auto *ih_382 = buffer.data(ih + 382);
    const auto *ih_383 = buffer.data(ih + 383);
    const auto *ih_384 = buffer.data(ih + 384);
    const auto *ih_385 = buffer.data(ih + 385);
    const auto *ih_386 = buffer.data(ih + 386);
    const auto *ih_387 = buffer.data(ih + 387);
    const auto *ih_388 = buffer.data(ih + 388);
    const auto *ih_389 = buffer.data(ih + 389);
    const auto *ih_390 = buffer.data(ih + 390);
    const auto *ih_391 = buffer.data(ih + 391);
    const auto *ih_392 = buffer.data(ih + 392);
    const auto *ih_393 = buffer.data(ih + 393);
    const auto *ih_394 = buffer.data(ih + 394);
    const auto *ih_395 = buffer.data(ih + 395);
    const auto *ih_396 = buffer.data(ih + 396);
    const auto *ih_397 = buffer.data(ih + 397);
    const auto *ih_398 = buffer.data(ih + 398);
    const auto *ih_399 = buffer.data(ih + 399);
    const auto *ih_400 = buffer.data(ih + 400);
    const auto *ih_401 = buffer.data(ih + 401);
    const auto *ih_402 = buffer.data(ih + 402);
    const auto *ih_403 = buffer.data(ih + 403);
    const auto *ih_404 = buffer.data(ih + 404);
    const auto *ih_405 = buffer.data(ih + 405);
    const auto *ih_406 = buffer.data(ih + 406);
    const auto *ih_407 = buffer.data(ih + 407);
    const auto *ih_408 = buffer.data(ih + 408);
    const auto *ih_409 = buffer.data(ih + 409);
    const auto *ih_410 = buffer.data(ih + 410);
    const auto *ih_411 = buffer.data(ih + 411);
    const auto *ih_412 = buffer.data(ih + 412);
    const auto *ih_413 = buffer.data(ih + 413);
    const auto *ih_414 = buffer.data(ih + 414);
    const auto *ih_415 = buffer.data(ih + 415);
    const auto *ih_416 = buffer.data(ih + 416);
    const auto *ih_417 = buffer.data(ih + 417);
    const auto *ih_418 = buffer.data(ih + 418);
    const auto *ih_419 = buffer.data(ih + 419);
    const auto *ih_420 = buffer.data(ih + 420);
    const auto *ih_421 = buffer.data(ih + 421);
    const auto *ih_422 = buffer.data(ih + 422);
    const auto *ih_423 = buffer.data(ih + 423);
    const auto *ih_424 = buffer.data(ih + 424);
    const auto *ih_425 = buffer.data(ih + 425);
    const auto *ih_426 = buffer.data(ih + 426);
    const auto *ih_427 = buffer.data(ih + 427);
    const auto *ih_428 = buffer.data(ih + 428);
    const auto *ih_429 = buffer.data(ih + 429);
    const auto *ih_430 = buffer.data(ih + 430);
    const auto *ih_431 = buffer.data(ih + 431);
    const auto *ih_432 = buffer.data(ih + 432);
    const auto *ih_433 = buffer.data(ih + 433);
    const auto *ih_434 = buffer.data(ih + 434);
    const auto *ih_435 = buffer.data(ih + 435);
    const auto *ih_436 = buffer.data(ih + 436);
    const auto *ih_437 = buffer.data(ih + 437);
    const auto *ih_438 = buffer.data(ih + 438);
    const auto *ih_439 = buffer.data(ih + 439);
    const auto *ih_440 = buffer.data(ih + 440);
    const auto *ih_441 = buffer.data(ih + 441);
    const auto *ih_442 = buffer.data(ih + 442);
    const auto *ih_443 = buffer.data(ih + 443);
    const auto *ih_444 = buffer.data(ih + 444);
    const auto *ih_445 = buffer.data(ih + 445);
    const auto *ih_446 = buffer.data(ih + 446);
    const auto *ih_447 = buffer.data(ih + 447);
    const auto *ih_448 = buffer.data(ih + 448);
    const auto *ih_449 = buffer.data(ih + 449);

    const auto *lh_300 = buffer.data(lh + 300);
    const auto *lh_301 = buffer.data(lh + 301);
    const auto *lh_302 = buffer.data(lh + 302);
    const auto *lh_303 = buffer.data(lh + 303);
    const auto *lh_304 = buffer.data(lh + 304);
    const auto *lh_305 = buffer.data(lh + 305);
    const auto *lh_306 = buffer.data(lh + 306);
    const auto *lh_307 = buffer.data(lh + 307);
    const auto *lh_308 = buffer.data(lh + 308);
    const auto *lh_309 = buffer.data(lh + 309);
    const auto *lh_310 = buffer.data(lh + 310);
    const auto *lh_311 = buffer.data(lh + 311);
    const auto *lh_312 = buffer.data(lh + 312);
    const auto *lh_313 = buffer.data(lh + 313);
    const auto *lh_314 = buffer.data(lh + 314);
    const auto *lh_315 = buffer.data(lh + 315);
    const auto *lh_316 = buffer.data(lh + 316);
    const auto *lh_317 = buffer.data(lh + 317);
    const auto *lh_318 = buffer.data(lh + 318);
    const auto *lh_319 = buffer.data(lh + 319);
    const auto *lh_320 = buffer.data(lh + 320);
    const auto *lh_321 = buffer.data(lh + 321);
    const auto *lh_322 = buffer.data(lh + 322);
    const auto *lh_323 = buffer.data(lh + 323);
    const auto *lh_324 = buffer.data(lh + 324);
    const auto *lh_325 = buffer.data(lh + 325);
    const auto *lh_326 = buffer.data(lh + 326);
    const auto *lh_327 = buffer.data(lh + 327);
    const auto *lh_328 = buffer.data(lh + 328);
    const auto *lh_329 = buffer.data(lh + 329);
    const auto *lh_330 = buffer.data(lh + 330);
    const auto *lh_331 = buffer.data(lh + 331);
    const auto *lh_332 = buffer.data(lh + 332);
    const auto *lh_333 = buffer.data(lh + 333);
    const auto *lh_334 = buffer.data(lh + 334);
    const auto *lh_335 = buffer.data(lh + 335);
    const auto *lh_336 = buffer.data(lh + 336);
    const auto *lh_337 = buffer.data(lh + 337);
    const auto *lh_338 = buffer.data(lh + 338);
    const auto *lh_339 = buffer.data(lh + 339);
    const auto *lh_340 = buffer.data(lh + 340);
    const auto *lh_341 = buffer.data(lh + 341);
    const auto *lh_342 = buffer.data(lh + 342);
    const auto *lh_343 = buffer.data(lh + 343);
    const auto *lh_344 = buffer.data(lh + 344);
    const auto *lh_345 = buffer.data(lh + 345);
    const auto *lh_346 = buffer.data(lh + 346);
    const auto *lh_347 = buffer.data(lh + 347);
    const auto *lh_348 = buffer.data(lh + 348);
    const auto *lh_349 = buffer.data(lh + 349);
    const auto *lh_350 = buffer.data(lh + 350);
    const auto *lh_351 = buffer.data(lh + 351);
    const auto *lh_352 = buffer.data(lh + 352);
    const auto *lh_353 = buffer.data(lh + 353);
    const auto *lh_354 = buffer.data(lh + 354);
    const auto *lh_355 = buffer.data(lh + 355);
    const auto *lh_356 = buffer.data(lh + 356);
    const auto *lh_357 = buffer.data(lh + 357);
    const auto *lh_358 = buffer.data(lh + 358);
    const auto *lh_359 = buffer.data(lh + 359);
    const auto *lh_360 = buffer.data(lh + 360);
    const auto *lh_361 = buffer.data(lh + 361);
    const auto *lh_362 = buffer.data(lh + 362);
    const auto *lh_363 = buffer.data(lh + 363);
    const auto *lh_364 = buffer.data(lh + 364);
    const auto *lh_365 = buffer.data(lh + 365);
    const auto *lh_366 = buffer.data(lh + 366);
    const auto *lh_367 = buffer.data(lh + 367);
    const auto *lh_368 = buffer.data(lh + 368);
    const auto *lh_369 = buffer.data(lh + 369);
    const auto *lh_370 = buffer.data(lh + 370);
    const auto *lh_371 = buffer.data(lh + 371);
    const auto *lh_372 = buffer.data(lh + 372);
    const auto *lh_373 = buffer.data(lh + 373);
    const auto *lh_374 = buffer.data(lh + 374);
    const auto *lh_375 = buffer.data(lh + 375);
    const auto *lh_376 = buffer.data(lh + 376);
    const auto *lh_377 = buffer.data(lh + 377);
    const auto *lh_378 = buffer.data(lh + 378);
    const auto *lh_379 = buffer.data(lh + 379);
    const auto *lh_380 = buffer.data(lh + 380);
    const auto *lh_381 = buffer.data(lh + 381);
    const auto *lh_382 = buffer.data(lh + 382);
    const auto *lh_383 = buffer.data(lh + 383);
    const auto *lh_384 = buffer.data(lh + 384);
    const auto *lh_385 = buffer.data(lh + 385);
    const auto *lh_386 = buffer.data(lh + 386);
    const auto *lh_387 = buffer.data(lh + 387);
    const auto *lh_388 = buffer.data(lh + 388);
    const auto *lh_389 = buffer.data(lh + 389);
    const auto *lh_390 = buffer.data(lh + 390);
    const auto *lh_391 = buffer.data(lh + 391);
    const auto *lh_392 = buffer.data(lh + 392);
    const auto *lh_393 = buffer.data(lh + 393);
    const auto *lh_394 = buffer.data(lh + 394);
    const auto *lh_395 = buffer.data(lh + 395);
    const auto *lh_396 = buffer.data(lh + 396);
    const auto *lh_397 = buffer.data(lh + 397);
    const auto *lh_398 = buffer.data(lh + 398);
    const auto *lh_399 = buffer.data(lh + 399);
    const auto *lh_400 = buffer.data(lh + 400);
    const auto *lh_401 = buffer.data(lh + 401);
    const auto *lh_402 = buffer.data(lh + 402);
    const auto *lh_403 = buffer.data(lh + 403);
    const auto *lh_404 = buffer.data(lh + 404);
    const auto *lh_405 = buffer.data(lh + 405);
    const auto *lh_406 = buffer.data(lh + 406);
    const auto *lh_407 = buffer.data(lh + 407);
    const auto *lh_408 = buffer.data(lh + 408);
    const auto *lh_409 = buffer.data(lh + 409);
    const auto *lh_410 = buffer.data(lh + 410);
    const auto *lh_411 = buffer.data(lh + 411);
    const auto *lh_412 = buffer.data(lh + 412);
    const auto *lh_413 = buffer.data(lh + 413);
    const auto *lh_414 = buffer.data(lh + 414);
    const auto *lh_415 = buffer.data(lh + 415);
    const auto *lh_416 = buffer.data(lh + 416);
    const auto *lh_417 = buffer.data(lh + 417);
    const auto *lh_418 = buffer.data(lh + 418);
    const auto *lh_419 = buffer.data(lh + 419);
    const auto *lh_420 = buffer.data(lh + 420);
    const auto *lh_421 = buffer.data(lh + 421);
    const auto *lh_422 = buffer.data(lh + 422);
    const auto *lh_423 = buffer.data(lh + 423);
    const auto *lh_424 = buffer.data(lh + 424);
    const auto *lh_425 = buffer.data(lh + 425);
    const auto *lh_426 = buffer.data(lh + 426);
    const auto *lh_427 = buffer.data(lh + 427);
    const auto *lh_428 = buffer.data(lh + 428);
    const auto *lh_429 = buffer.data(lh + 429);
    const auto *lh_430 = buffer.data(lh + 430);
    const auto *lh_431 = buffer.data(lh + 431);
    const auto *lh_432 = buffer.data(lh + 432);
    const auto *lh_433 = buffer.data(lh + 433);
    const auto *lh_434 = buffer.data(lh + 434);
    const auto *lh_435 = buffer.data(lh + 435);
    const auto *lh_436 = buffer.data(lh + 436);
    const auto *lh_437 = buffer.data(lh + 437);
    const auto *lh_438 = buffer.data(lh + 438);
    const auto *lh_439 = buffer.data(lh + 439);
    const auto *lh_440 = buffer.data(lh + 440);
    const auto *lh_441 = buffer.data(lh + 441);
    const auto *lh_442 = buffer.data(lh + 442);
    const auto *lh_443 = buffer.data(lh + 443);
    const auto *lh_444 = buffer.data(lh + 444);
    const auto *lh_445 = buffer.data(lh + 445);
    const auto *lh_446 = buffer.data(lh + 446);
    const auto *lh_447 = buffer.data(lh + 447);
    const auto *lh_448 = buffer.data(lh + 448);
    const auto *lh_449 = buffer.data(lh + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ih_300, ih_301, ih_302, ih_303, \
                         ih_304, lh_300, lh_301, lh_302, lh_303, \
                         lh_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -3.0 * ih_300[k]
                   + f_0 * lh_300[k];

        t_301[k] = -3.0 * ih_301[k]
                   + f_0 * lh_301[k];

        t_302[k] = -3.0 * ih_302[k]
                   + f_0 * lh_302[k];

        t_303[k] = -3.0 * ih_303[k]
                   + f_0 * lh_303[k];

        t_304[k] = -3.0 * ih_304[k]
                   + f_0 * lh_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ih_305, ih_306, ih_307, ih_308, \
                         ih_309, lh_305, lh_306, lh_307, lh_308, \
                         lh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -3.0 * ih_305[k]
                   + f_0 * lh_305[k];

        t_306[k] = -3.0 * ih_306[k]
                   + f_0 * lh_306[k];

        t_307[k] = -3.0 * ih_307[k]
                   + f_0 * lh_307[k];

        t_308[k] = -3.0 * ih_308[k]
                   + f_0 * lh_308[k];

        t_309[k] = -3.0 * ih_309[k]
                   + f_0 * lh_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ih_310, ih_311, ih_312, ih_313, \
                         ih_314, lh_310, lh_311, lh_312, lh_313, \
                         lh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * ih_310[k]
                   + f_0 * lh_310[k];

        t_311[k] = -3.0 * ih_311[k]
                   + f_0 * lh_311[k];

        t_312[k] = -3.0 * ih_312[k]
                   + f_0 * lh_312[k];

        t_313[k] = -3.0 * ih_313[k]
                   + f_0 * lh_313[k];

        t_314[k] = -3.0 * ih_314[k]
                   + f_0 * lh_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ih_315, ih_316, ih_317, ih_318, \
                         ih_319, lh_315, lh_316, lh_317, lh_318, \
                         lh_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -2.0 * ih_315[k]
                   + f_0 * lh_315[k];

        t_316[k] = -2.0 * ih_316[k]
                   + f_0 * lh_316[k];

        t_317[k] = -2.0 * ih_317[k]
                   + f_0 * lh_317[k];

        t_318[k] = -2.0 * ih_318[k]
                   + f_0 * lh_318[k];

        t_319[k] = -2.0 * ih_319[k]
                   + f_0 * lh_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ih_320, ih_321, ih_322, ih_323, \
                         ih_324, lh_320, lh_321, lh_322, lh_323, \
                         lh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -2.0 * ih_320[k]
                   + f_0 * lh_320[k];

        t_321[k] = -2.0 * ih_321[k]
                   + f_0 * lh_321[k];

        t_322[k] = -2.0 * ih_322[k]
                   + f_0 * lh_322[k];

        t_323[k] = -2.0 * ih_323[k]
                   + f_0 * lh_323[k];

        t_324[k] = -2.0 * ih_324[k]
                   + f_0 * lh_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ih_325, ih_326, ih_327, ih_328, \
                         ih_329, lh_325, lh_326, lh_327, lh_328, \
                         lh_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -2.0 * ih_325[k]
                   + f_0 * lh_325[k];

        t_326[k] = -2.0 * ih_326[k]
                   + f_0 * lh_326[k];

        t_327[k] = -2.0 * ih_327[k]
                   + f_0 * lh_327[k];

        t_328[k] = -2.0 * ih_328[k]
                   + f_0 * lh_328[k];

        t_329[k] = -2.0 * ih_329[k]
                   + f_0 * lh_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ih_330, ih_331, ih_332, ih_333, \
                         ih_334, lh_330, lh_331, lh_332, lh_333, \
                         lh_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -2.0 * ih_330[k]
                   + f_0 * lh_330[k];

        t_331[k] = -2.0 * ih_331[k]
                   + f_0 * lh_331[k];

        t_332[k] = -2.0 * ih_332[k]
                   + f_0 * lh_332[k];

        t_333[k] = -2.0 * ih_333[k]
                   + f_0 * lh_333[k];

        t_334[k] = -2.0 * ih_334[k]
                   + f_0 * lh_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ih_335, ih_336, ih_337, ih_338, \
                         ih_339, lh_335, lh_336, lh_337, lh_338, \
                         lh_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -2.0 * ih_335[k]
                   + f_0 * lh_335[k];

        t_336[k] = -2.0 * ih_336[k]
                   + f_0 * lh_336[k];

        t_337[k] = -2.0 * ih_337[k]
                   + f_0 * lh_337[k];

        t_338[k] = -2.0 * ih_338[k]
                   + f_0 * lh_338[k];

        t_339[k] = -2.0 * ih_339[k]
                   + f_0 * lh_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ih_340, ih_341, ih_342, ih_343, \
                         ih_344, lh_340, lh_341, lh_342, lh_343, \
                         lh_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * ih_340[k]
                   + f_0 * lh_340[k];

        t_341[k] = -2.0 * ih_341[k]
                   + f_0 * lh_341[k];

        t_342[k] = -2.0 * ih_342[k]
                   + f_0 * lh_342[k];

        t_343[k] = -2.0 * ih_343[k]
                   + f_0 * lh_343[k];

        t_344[k] = -2.0 * ih_344[k]
                   + f_0 * lh_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ih_345, ih_346, ih_347, ih_348, \
                         ih_349, lh_345, lh_346, lh_347, lh_348, \
                         lh_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * ih_345[k]
                   + f_0 * lh_345[k];

        t_346[k] = -2.0 * ih_346[k]
                   + f_0 * lh_346[k];

        t_347[k] = -2.0 * ih_347[k]
                   + f_0 * lh_347[k];

        t_348[k] = -2.0 * ih_348[k]
                   + f_0 * lh_348[k];

        t_349[k] = -2.0 * ih_349[k]
                   + f_0 * lh_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ih_350, ih_351, ih_352, ih_353, \
                         ih_354, lh_350, lh_351, lh_352, lh_353, \
                         lh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * ih_350[k]
                   + f_0 * lh_350[k];

        t_351[k] = -2.0 * ih_351[k]
                   + f_0 * lh_351[k];

        t_352[k] = -2.0 * ih_352[k]
                   + f_0 * lh_352[k];

        t_353[k] = -2.0 * ih_353[k]
                   + f_0 * lh_353[k];

        t_354[k] = -2.0 * ih_354[k]
                   + f_0 * lh_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ih_355, ih_356, ih_357, ih_358, \
                         ih_359, lh_355, lh_356, lh_357, lh_358, \
                         lh_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * ih_355[k]
                   + f_0 * lh_355[k];

        t_356[k] = -2.0 * ih_356[k]
                   + f_0 * lh_356[k];

        t_357[k] = -2.0 * ih_357[k]
                   + f_0 * lh_357[k];

        t_358[k] = -2.0 * ih_358[k]
                   + f_0 * lh_358[k];

        t_359[k] = -2.0 * ih_359[k]
                   + f_0 * lh_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ih_360, ih_361, ih_362, ih_363, \
                         ih_364, lh_360, lh_361, lh_362, lh_363, \
                         lh_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -2.0 * ih_360[k]
                   + f_0 * lh_360[k];

        t_361[k] = -2.0 * ih_361[k]
                   + f_0 * lh_361[k];

        t_362[k] = -2.0 * ih_362[k]
                   + f_0 * lh_362[k];

        t_363[k] = -2.0 * ih_363[k]
                   + f_0 * lh_363[k];

        t_364[k] = -2.0 * ih_364[k]
                   + f_0 * lh_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ih_365, ih_366, ih_367, ih_368, \
                         ih_369, lh_365, lh_366, lh_367, lh_368, \
                         lh_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -2.0 * ih_365[k]
                   + f_0 * lh_365[k];

        t_366[k] = -2.0 * ih_366[k]
                   + f_0 * lh_366[k];

        t_367[k] = -2.0 * ih_367[k]
                   + f_0 * lh_367[k];

        t_368[k] = -2.0 * ih_368[k]
                   + f_0 * lh_368[k];

        t_369[k] = -2.0 * ih_369[k]
                   + f_0 * lh_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ih_370, ih_371, ih_372, ih_373, \
                         ih_374, lh_370, lh_371, lh_372, lh_373, \
                         lh_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -2.0 * ih_370[k]
                   + f_0 * lh_370[k];

        t_371[k] = -2.0 * ih_371[k]
                   + f_0 * lh_371[k];

        t_372[k] = -2.0 * ih_372[k]
                   + f_0 * lh_372[k];

        t_373[k] = -2.0 * ih_373[k]
                   + f_0 * lh_373[k];

        t_374[k] = -2.0 * ih_374[k]
                   + f_0 * lh_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ih_375, ih_376, ih_377, ih_378, \
                         ih_379, lh_375, lh_376, lh_377, lh_378, \
                         lh_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -2.0 * ih_375[k]
                   + f_0 * lh_375[k];

        t_376[k] = -2.0 * ih_376[k]
                   + f_0 * lh_376[k];

        t_377[k] = -2.0 * ih_377[k]
                   + f_0 * lh_377[k];

        t_378[k] = -2.0 * ih_378[k]
                   + f_0 * lh_378[k];

        t_379[k] = -2.0 * ih_379[k]
                   + f_0 * lh_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ih_380, ih_381, ih_382, ih_383, \
                         ih_384, lh_380, lh_381, lh_382, lh_383, \
                         lh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * ih_380[k]
                   + f_0 * lh_380[k];

        t_381[k] = -2.0 * ih_381[k]
                   + f_0 * lh_381[k];

        t_382[k] = -2.0 * ih_382[k]
                   + f_0 * lh_382[k];

        t_383[k] = -2.0 * ih_383[k]
                   + f_0 * lh_383[k];

        t_384[k] = -2.0 * ih_384[k]
                   + f_0 * lh_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ih_385, ih_386, ih_387, ih_388, \
                         ih_389, lh_385, lh_386, lh_387, lh_388, \
                         lh_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * ih_385[k]
                   + f_0 * lh_385[k];

        t_386[k] = -2.0 * ih_386[k]
                   + f_0 * lh_386[k];

        t_387[k] = -2.0 * ih_387[k]
                   + f_0 * lh_387[k];

        t_388[k] = -2.0 * ih_388[k]
                   + f_0 * lh_388[k];

        t_389[k] = -2.0 * ih_389[k]
                   + f_0 * lh_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ih_390, ih_391, ih_392, ih_393, \
                         ih_394, lh_390, lh_391, lh_392, lh_393, \
                         lh_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -2.0 * ih_390[k]
                   + f_0 * lh_390[k];

        t_391[k] = -2.0 * ih_391[k]
                   + f_0 * lh_391[k];

        t_392[k] = -2.0 * ih_392[k]
                   + f_0 * lh_392[k];

        t_393[k] = -2.0 * ih_393[k]
                   + f_0 * lh_393[k];

        t_394[k] = -2.0 * ih_394[k]
                   + f_0 * lh_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ih_395, ih_396, ih_397, ih_398, \
                         ih_399, lh_395, lh_396, lh_397, lh_398, \
                         lh_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -2.0 * ih_395[k]
                   + f_0 * lh_395[k];

        t_396[k] = -2.0 * ih_396[k]
                   + f_0 * lh_396[k];

        t_397[k] = -2.0 * ih_397[k]
                   + f_0 * lh_397[k];

        t_398[k] = -2.0 * ih_398[k]
                   + f_0 * lh_398[k];

        t_399[k] = -2.0 * ih_399[k]
                   + f_0 * lh_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ih_400, ih_401, ih_402, ih_403, \
                         ih_404, lh_400, lh_401, lh_402, lh_403, \
                         lh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -2.0 * ih_400[k]
                   + f_0 * lh_400[k];

        t_401[k] = -2.0 * ih_401[k]
                   + f_0 * lh_401[k];

        t_402[k] = -2.0 * ih_402[k]
                   + f_0 * lh_402[k];

        t_403[k] = -2.0 * ih_403[k]
                   + f_0 * lh_403[k];

        t_404[k] = -2.0 * ih_404[k]
                   + f_0 * lh_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ih_405, ih_406, ih_407, ih_408, \
                         ih_409, lh_405, lh_406, lh_407, lh_408, \
                         lh_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -2.0 * ih_405[k]
                   + f_0 * lh_405[k];

        t_406[k] = -2.0 * ih_406[k]
                   + f_0 * lh_406[k];

        t_407[k] = -2.0 * ih_407[k]
                   + f_0 * lh_407[k];

        t_408[k] = -2.0 * ih_408[k]
                   + f_0 * lh_408[k];

        t_409[k] = -2.0 * ih_409[k]
                   + f_0 * lh_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ih_410, ih_411, ih_412, ih_413, \
                         ih_414, lh_410, lh_411, lh_412, lh_413, \
                         lh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -2.0 * ih_410[k]
                   + f_0 * lh_410[k];

        t_411[k] = -2.0 * ih_411[k]
                   + f_0 * lh_411[k];

        t_412[k] = -2.0 * ih_412[k]
                   + f_0 * lh_412[k];

        t_413[k] = -2.0 * ih_413[k]
                   + f_0 * lh_413[k];

        t_414[k] = -2.0 * ih_414[k]
                   + f_0 * lh_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ih_415, ih_416, ih_417, ih_418, \
                         ih_419, lh_415, lh_416, lh_417, lh_418, \
                         lh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -2.0 * ih_415[k]
                   + f_0 * lh_415[k];

        t_416[k] = -2.0 * ih_416[k]
                   + f_0 * lh_416[k];

        t_417[k] = -2.0 * ih_417[k]
                   + f_0 * lh_417[k];

        t_418[k] = -2.0 * ih_418[k]
                   + f_0 * lh_418[k];

        t_419[k] = -2.0 * ih_419[k]
                   + f_0 * lh_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ih_420, ih_421, ih_422, ih_423, \
                         ih_424, lh_420, lh_421, lh_422, lh_423, \
                         lh_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -2.0 * ih_420[k]
                   + f_0 * lh_420[k];

        t_421[k] = -2.0 * ih_421[k]
                   + f_0 * lh_421[k];

        t_422[k] = -2.0 * ih_422[k]
                   + f_0 * lh_422[k];

        t_423[k] = -2.0 * ih_423[k]
                   + f_0 * lh_423[k];

        t_424[k] = -2.0 * ih_424[k]
                   + f_0 * lh_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ih_425, ih_426, ih_427, ih_428, \
                         ih_429, lh_425, lh_426, lh_427, lh_428, \
                         lh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -2.0 * ih_425[k]
                   + f_0 * lh_425[k];

        t_426[k] = -2.0 * ih_426[k]
                   + f_0 * lh_426[k];

        t_427[k] = -2.0 * ih_427[k]
                   + f_0 * lh_427[k];

        t_428[k] = -2.0 * ih_428[k]
                   + f_0 * lh_428[k];

        t_429[k] = -2.0 * ih_429[k]
                   + f_0 * lh_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ih_430, ih_431, ih_432, ih_433, \
                         ih_434, lh_430, lh_431, lh_432, lh_433, \
                         lh_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -2.0 * ih_430[k]
                   + f_0 * lh_430[k];

        t_431[k] = -2.0 * ih_431[k]
                   + f_0 * lh_431[k];

        t_432[k] = -2.0 * ih_432[k]
                   + f_0 * lh_432[k];

        t_433[k] = -2.0 * ih_433[k]
                   + f_0 * lh_433[k];

        t_434[k] = -2.0 * ih_434[k]
                   + f_0 * lh_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ih_435, ih_436, ih_437, ih_438, \
                         ih_439, lh_435, lh_436, lh_437, lh_438, \
                         lh_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -2.0 * ih_435[k]
                   + f_0 * lh_435[k];

        t_436[k] = -2.0 * ih_436[k]
                   + f_0 * lh_436[k];

        t_437[k] = -2.0 * ih_437[k]
                   + f_0 * lh_437[k];

        t_438[k] = -2.0 * ih_438[k]
                   + f_0 * lh_438[k];

        t_439[k] = -2.0 * ih_439[k]
                   + f_0 * lh_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ih_440, ih_441, ih_442, ih_443, \
                         ih_444, lh_440, lh_441, lh_442, lh_443, \
                         lh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -2.0 * ih_440[k]
                   + f_0 * lh_440[k];

        t_441[k] = -ih_441[k]
                   + f_0 * lh_441[k];

        t_442[k] = -ih_442[k]
                   + f_0 * lh_442[k];

        t_443[k] = -ih_443[k]
                   + f_0 * lh_443[k];

        t_444[k] = -ih_444[k]
                   + f_0 * lh_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ih_445, ih_446, ih_447, ih_448, \
                         ih_449, lh_445, lh_446, lh_447, lh_448, \
                         lh_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -ih_445[k]
                   + f_0 * lh_445[k];

        t_446[k] = -ih_446[k]
                   + f_0 * lh_446[k];

        t_447[k] = -ih_447[k]
                   + f_0 * lh_447[k];

        t_448[k] = -ih_448[k]
                   + f_0 * lh_448[k];

        t_449[k] = -ih_449[k]
                   + f_0 * lh_449[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
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

    const auto *ih_450 = buffer.data(ih + 450);
    const auto *ih_451 = buffer.data(ih + 451);
    const auto *ih_452 = buffer.data(ih + 452);
    const auto *ih_453 = buffer.data(ih + 453);
    const auto *ih_454 = buffer.data(ih + 454);
    const auto *ih_455 = buffer.data(ih + 455);
    const auto *ih_456 = buffer.data(ih + 456);
    const auto *ih_457 = buffer.data(ih + 457);
    const auto *ih_458 = buffer.data(ih + 458);
    const auto *ih_459 = buffer.data(ih + 459);
    const auto *ih_460 = buffer.data(ih + 460);
    const auto *ih_461 = buffer.data(ih + 461);
    const auto *ih_462 = buffer.data(ih + 462);
    const auto *ih_463 = buffer.data(ih + 463);
    const auto *ih_464 = buffer.data(ih + 464);
    const auto *ih_465 = buffer.data(ih + 465);
    const auto *ih_466 = buffer.data(ih + 466);
    const auto *ih_467 = buffer.data(ih + 467);
    const auto *ih_468 = buffer.data(ih + 468);
    const auto *ih_469 = buffer.data(ih + 469);
    const auto *ih_470 = buffer.data(ih + 470);
    const auto *ih_471 = buffer.data(ih + 471);
    const auto *ih_472 = buffer.data(ih + 472);
    const auto *ih_473 = buffer.data(ih + 473);
    const auto *ih_474 = buffer.data(ih + 474);
    const auto *ih_475 = buffer.data(ih + 475);
    const auto *ih_476 = buffer.data(ih + 476);
    const auto *ih_477 = buffer.data(ih + 477);
    const auto *ih_478 = buffer.data(ih + 478);
    const auto *ih_479 = buffer.data(ih + 479);
    const auto *ih_480 = buffer.data(ih + 480);
    const auto *ih_481 = buffer.data(ih + 481);
    const auto *ih_482 = buffer.data(ih + 482);
    const auto *ih_483 = buffer.data(ih + 483);
    const auto *ih_484 = buffer.data(ih + 484);
    const auto *ih_485 = buffer.data(ih + 485);
    const auto *ih_486 = buffer.data(ih + 486);
    const auto *ih_487 = buffer.data(ih + 487);
    const auto *ih_488 = buffer.data(ih + 488);
    const auto *ih_489 = buffer.data(ih + 489);
    const auto *ih_490 = buffer.data(ih + 490);
    const auto *ih_491 = buffer.data(ih + 491);
    const auto *ih_492 = buffer.data(ih + 492);
    const auto *ih_493 = buffer.data(ih + 493);
    const auto *ih_494 = buffer.data(ih + 494);
    const auto *ih_495 = buffer.data(ih + 495);
    const auto *ih_496 = buffer.data(ih + 496);
    const auto *ih_497 = buffer.data(ih + 497);
    const auto *ih_498 = buffer.data(ih + 498);
    const auto *ih_499 = buffer.data(ih + 499);
    const auto *ih_500 = buffer.data(ih + 500);
    const auto *ih_501 = buffer.data(ih + 501);
    const auto *ih_502 = buffer.data(ih + 502);
    const auto *ih_503 = buffer.data(ih + 503);
    const auto *ih_504 = buffer.data(ih + 504);
    const auto *ih_505 = buffer.data(ih + 505);
    const auto *ih_506 = buffer.data(ih + 506);
    const auto *ih_507 = buffer.data(ih + 507);
    const auto *ih_508 = buffer.data(ih + 508);
    const auto *ih_509 = buffer.data(ih + 509);
    const auto *ih_510 = buffer.data(ih + 510);
    const auto *ih_511 = buffer.data(ih + 511);
    const auto *ih_512 = buffer.data(ih + 512);
    const auto *ih_513 = buffer.data(ih + 513);
    const auto *ih_514 = buffer.data(ih + 514);
    const auto *ih_515 = buffer.data(ih + 515);
    const auto *ih_516 = buffer.data(ih + 516);
    const auto *ih_517 = buffer.data(ih + 517);
    const auto *ih_518 = buffer.data(ih + 518);
    const auto *ih_519 = buffer.data(ih + 519);
    const auto *ih_520 = buffer.data(ih + 520);
    const auto *ih_521 = buffer.data(ih + 521);
    const auto *ih_522 = buffer.data(ih + 522);
    const auto *ih_523 = buffer.data(ih + 523);
    const auto *ih_524 = buffer.data(ih + 524);
    const auto *ih_525 = buffer.data(ih + 525);
    const auto *ih_526 = buffer.data(ih + 526);
    const auto *ih_527 = buffer.data(ih + 527);
    const auto *ih_528 = buffer.data(ih + 528);
    const auto *ih_529 = buffer.data(ih + 529);
    const auto *ih_530 = buffer.data(ih + 530);
    const auto *ih_531 = buffer.data(ih + 531);
    const auto *ih_532 = buffer.data(ih + 532);
    const auto *ih_533 = buffer.data(ih + 533);
    const auto *ih_534 = buffer.data(ih + 534);
    const auto *ih_535 = buffer.data(ih + 535);
    const auto *ih_536 = buffer.data(ih + 536);
    const auto *ih_537 = buffer.data(ih + 537);
    const auto *ih_538 = buffer.data(ih + 538);
    const auto *ih_539 = buffer.data(ih + 539);
    const auto *ih_540 = buffer.data(ih + 540);
    const auto *ih_541 = buffer.data(ih + 541);
    const auto *ih_542 = buffer.data(ih + 542);
    const auto *ih_543 = buffer.data(ih + 543);
    const auto *ih_544 = buffer.data(ih + 544);
    const auto *ih_545 = buffer.data(ih + 545);
    const auto *ih_546 = buffer.data(ih + 546);
    const auto *ih_547 = buffer.data(ih + 547);
    const auto *ih_548 = buffer.data(ih + 548);
    const auto *ih_549 = buffer.data(ih + 549);
    const auto *ih_550 = buffer.data(ih + 550);
    const auto *ih_551 = buffer.data(ih + 551);
    const auto *ih_552 = buffer.data(ih + 552);
    const auto *ih_553 = buffer.data(ih + 553);
    const auto *ih_554 = buffer.data(ih + 554);
    const auto *ih_555 = buffer.data(ih + 555);
    const auto *ih_556 = buffer.data(ih + 556);
    const auto *ih_557 = buffer.data(ih + 557);
    const auto *ih_558 = buffer.data(ih + 558);
    const auto *ih_559 = buffer.data(ih + 559);
    const auto *ih_560 = buffer.data(ih + 560);
    const auto *ih_561 = buffer.data(ih + 561);
    const auto *ih_562 = buffer.data(ih + 562);
    const auto *ih_563 = buffer.data(ih + 563);
    const auto *ih_564 = buffer.data(ih + 564);
    const auto *ih_565 = buffer.data(ih + 565);
    const auto *ih_566 = buffer.data(ih + 566);
    const auto *ih_567 = buffer.data(ih + 567);
    const auto *ih_568 = buffer.data(ih + 568);
    const auto *ih_569 = buffer.data(ih + 569);
    const auto *ih_570 = buffer.data(ih + 570);
    const auto *ih_571 = buffer.data(ih + 571);
    const auto *ih_572 = buffer.data(ih + 572);
    const auto *ih_573 = buffer.data(ih + 573);
    const auto *ih_574 = buffer.data(ih + 574);
    const auto *ih_575 = buffer.data(ih + 575);
    const auto *ih_576 = buffer.data(ih + 576);
    const auto *ih_577 = buffer.data(ih + 577);
    const auto *ih_578 = buffer.data(ih + 578);
    const auto *ih_579 = buffer.data(ih + 579);
    const auto *ih_580 = buffer.data(ih + 580);
    const auto *ih_581 = buffer.data(ih + 581);
    const auto *ih_582 = buffer.data(ih + 582);
    const auto *ih_583 = buffer.data(ih + 583);
    const auto *ih_584 = buffer.data(ih + 584);
    const auto *ih_585 = buffer.data(ih + 585);
    const auto *ih_586 = buffer.data(ih + 586);
    const auto *ih_587 = buffer.data(ih + 587);

    const auto *lh_450 = buffer.data(lh + 450);
    const auto *lh_451 = buffer.data(lh + 451);
    const auto *lh_452 = buffer.data(lh + 452);
    const auto *lh_453 = buffer.data(lh + 453);
    const auto *lh_454 = buffer.data(lh + 454);
    const auto *lh_455 = buffer.data(lh + 455);
    const auto *lh_456 = buffer.data(lh + 456);
    const auto *lh_457 = buffer.data(lh + 457);
    const auto *lh_458 = buffer.data(lh + 458);
    const auto *lh_459 = buffer.data(lh + 459);
    const auto *lh_460 = buffer.data(lh + 460);
    const auto *lh_461 = buffer.data(lh + 461);
    const auto *lh_462 = buffer.data(lh + 462);
    const auto *lh_463 = buffer.data(lh + 463);
    const auto *lh_464 = buffer.data(lh + 464);
    const auto *lh_465 = buffer.data(lh + 465);
    const auto *lh_466 = buffer.data(lh + 466);
    const auto *lh_467 = buffer.data(lh + 467);
    const auto *lh_468 = buffer.data(lh + 468);
    const auto *lh_469 = buffer.data(lh + 469);
    const auto *lh_470 = buffer.data(lh + 470);
    const auto *lh_471 = buffer.data(lh + 471);
    const auto *lh_472 = buffer.data(lh + 472);
    const auto *lh_473 = buffer.data(lh + 473);
    const auto *lh_474 = buffer.data(lh + 474);
    const auto *lh_475 = buffer.data(lh + 475);
    const auto *lh_476 = buffer.data(lh + 476);
    const auto *lh_477 = buffer.data(lh + 477);
    const auto *lh_478 = buffer.data(lh + 478);
    const auto *lh_479 = buffer.data(lh + 479);
    const auto *lh_480 = buffer.data(lh + 480);
    const auto *lh_481 = buffer.data(lh + 481);
    const auto *lh_482 = buffer.data(lh + 482);
    const auto *lh_483 = buffer.data(lh + 483);
    const auto *lh_484 = buffer.data(lh + 484);
    const auto *lh_485 = buffer.data(lh + 485);
    const auto *lh_486 = buffer.data(lh + 486);
    const auto *lh_487 = buffer.data(lh + 487);
    const auto *lh_488 = buffer.data(lh + 488);
    const auto *lh_489 = buffer.data(lh + 489);
    const auto *lh_490 = buffer.data(lh + 490);
    const auto *lh_491 = buffer.data(lh + 491);
    const auto *lh_492 = buffer.data(lh + 492);
    const auto *lh_493 = buffer.data(lh + 493);
    const auto *lh_494 = buffer.data(lh + 494);
    const auto *lh_495 = buffer.data(lh + 495);
    const auto *lh_496 = buffer.data(lh + 496);
    const auto *lh_497 = buffer.data(lh + 497);
    const auto *lh_498 = buffer.data(lh + 498);
    const auto *lh_499 = buffer.data(lh + 499);
    const auto *lh_500 = buffer.data(lh + 500);
    const auto *lh_501 = buffer.data(lh + 501);
    const auto *lh_502 = buffer.data(lh + 502);
    const auto *lh_503 = buffer.data(lh + 503);
    const auto *lh_504 = buffer.data(lh + 504);
    const auto *lh_505 = buffer.data(lh + 505);
    const auto *lh_506 = buffer.data(lh + 506);
    const auto *lh_507 = buffer.data(lh + 507);
    const auto *lh_508 = buffer.data(lh + 508);
    const auto *lh_509 = buffer.data(lh + 509);
    const auto *lh_510 = buffer.data(lh + 510);
    const auto *lh_511 = buffer.data(lh + 511);
    const auto *lh_512 = buffer.data(lh + 512);
    const auto *lh_513 = buffer.data(lh + 513);
    const auto *lh_514 = buffer.data(lh + 514);
    const auto *lh_515 = buffer.data(lh + 515);
    const auto *lh_516 = buffer.data(lh + 516);
    const auto *lh_517 = buffer.data(lh + 517);
    const auto *lh_518 = buffer.data(lh + 518);
    const auto *lh_519 = buffer.data(lh + 519);
    const auto *lh_520 = buffer.data(lh + 520);
    const auto *lh_521 = buffer.data(lh + 521);
    const auto *lh_522 = buffer.data(lh + 522);
    const auto *lh_523 = buffer.data(lh + 523);
    const auto *lh_524 = buffer.data(lh + 524);
    const auto *lh_525 = buffer.data(lh + 525);
    const auto *lh_526 = buffer.data(lh + 526);
    const auto *lh_527 = buffer.data(lh + 527);
    const auto *lh_528 = buffer.data(lh + 528);
    const auto *lh_529 = buffer.data(lh + 529);
    const auto *lh_530 = buffer.data(lh + 530);
    const auto *lh_531 = buffer.data(lh + 531);
    const auto *lh_532 = buffer.data(lh + 532);
    const auto *lh_533 = buffer.data(lh + 533);
    const auto *lh_534 = buffer.data(lh + 534);
    const auto *lh_535 = buffer.data(lh + 535);
    const auto *lh_536 = buffer.data(lh + 536);
    const auto *lh_537 = buffer.data(lh + 537);
    const auto *lh_538 = buffer.data(lh + 538);
    const auto *lh_539 = buffer.data(lh + 539);
    const auto *lh_540 = buffer.data(lh + 540);
    const auto *lh_541 = buffer.data(lh + 541);
    const auto *lh_542 = buffer.data(lh + 542);
    const auto *lh_543 = buffer.data(lh + 543);
    const auto *lh_544 = buffer.data(lh + 544);
    const auto *lh_545 = buffer.data(lh + 545);
    const auto *lh_546 = buffer.data(lh + 546);
    const auto *lh_547 = buffer.data(lh + 547);
    const auto *lh_548 = buffer.data(lh + 548);
    const auto *lh_549 = buffer.data(lh + 549);
    const auto *lh_550 = buffer.data(lh + 550);
    const auto *lh_551 = buffer.data(lh + 551);
    const auto *lh_552 = buffer.data(lh + 552);
    const auto *lh_553 = buffer.data(lh + 553);
    const auto *lh_554 = buffer.data(lh + 554);
    const auto *lh_555 = buffer.data(lh + 555);
    const auto *lh_556 = buffer.data(lh + 556);
    const auto *lh_557 = buffer.data(lh + 557);
    const auto *lh_558 = buffer.data(lh + 558);
    const auto *lh_559 = buffer.data(lh + 559);
    const auto *lh_560 = buffer.data(lh + 560);
    const auto *lh_561 = buffer.data(lh + 561);
    const auto *lh_562 = buffer.data(lh + 562);
    const auto *lh_563 = buffer.data(lh + 563);
    const auto *lh_564 = buffer.data(lh + 564);
    const auto *lh_565 = buffer.data(lh + 565);
    const auto *lh_566 = buffer.data(lh + 566);
    const auto *lh_567 = buffer.data(lh + 567);
    const auto *lh_568 = buffer.data(lh + 568);
    const auto *lh_569 = buffer.data(lh + 569);
    const auto *lh_570 = buffer.data(lh + 570);
    const auto *lh_571 = buffer.data(lh + 571);
    const auto *lh_572 = buffer.data(lh + 572);
    const auto *lh_573 = buffer.data(lh + 573);
    const auto *lh_574 = buffer.data(lh + 574);
    const auto *lh_575 = buffer.data(lh + 575);
    const auto *lh_576 = buffer.data(lh + 576);
    const auto *lh_577 = buffer.data(lh + 577);
    const auto *lh_578 = buffer.data(lh + 578);
    const auto *lh_579 = buffer.data(lh + 579);
    const auto *lh_580 = buffer.data(lh + 580);
    const auto *lh_581 = buffer.data(lh + 581);
    const auto *lh_582 = buffer.data(lh + 582);
    const auto *lh_583 = buffer.data(lh + 583);
    const auto *lh_584 = buffer.data(lh + 584);
    const auto *lh_585 = buffer.data(lh + 585);
    const auto *lh_586 = buffer.data(lh + 586);
    const auto *lh_587 = buffer.data(lh + 587);
    const auto *lh_588 = buffer.data(lh + 588);
    const auto *lh_589 = buffer.data(lh + 589);
    const auto *lh_590 = buffer.data(lh + 590);
    const auto *lh_591 = buffer.data(lh + 591);
    const auto *lh_592 = buffer.data(lh + 592);
    const auto *lh_593 = buffer.data(lh + 593);
    const auto *lh_594 = buffer.data(lh + 594);
    const auto *lh_595 = buffer.data(lh + 595);
    const auto *lh_596 = buffer.data(lh + 596);
    const auto *lh_597 = buffer.data(lh + 597);
    const auto *lh_598 = buffer.data(lh + 598);
    const auto *lh_599 = buffer.data(lh + 599);
    const auto *lh_600 = buffer.data(lh + 600);
    const auto *lh_601 = buffer.data(lh + 601);
    const auto *lh_602 = buffer.data(lh + 602);
    const auto *lh_603 = buffer.data(lh + 603);
    const auto *lh_604 = buffer.data(lh + 604);
    const auto *lh_605 = buffer.data(lh + 605);
    const auto *lh_606 = buffer.data(lh + 606);

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ih_450, ih_451, ih_452, ih_453, \
                         ih_454, lh_450, lh_451, lh_452, lh_453, \
                         lh_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -ih_450[k]
                   + f_0 * lh_450[k];

        t_451[k] = -ih_451[k]
                   + f_0 * lh_451[k];

        t_452[k] = -ih_452[k]
                   + f_0 * lh_452[k];

        t_453[k] = -ih_453[k]
                   + f_0 * lh_453[k];

        t_454[k] = -ih_454[k]
                   + f_0 * lh_454[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ih_455, ih_456, ih_457, ih_458, \
                         ih_459, lh_455, lh_456, lh_457, lh_458, \
                         lh_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -ih_455[k]
                   + f_0 * lh_455[k];

        t_456[k] = -ih_456[k]
                   + f_0 * lh_456[k];

        t_457[k] = -ih_457[k]
                   + f_0 * lh_457[k];

        t_458[k] = -ih_458[k]
                   + f_0 * lh_458[k];

        t_459[k] = -ih_459[k]
                   + f_0 * lh_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ih_460, ih_461, ih_462, ih_463, \
                         ih_464, lh_460, lh_461, lh_462, lh_463, \
                         lh_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -ih_460[k]
                   + f_0 * lh_460[k];

        t_461[k] = -ih_461[k]
                   + f_0 * lh_461[k];

        t_462[k] = -ih_462[k]
                   + f_0 * lh_462[k];

        t_463[k] = -ih_463[k]
                   + f_0 * lh_463[k];

        t_464[k] = -ih_464[k]
                   + f_0 * lh_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ih_465, ih_466, ih_467, ih_468, \
                         ih_469, lh_465, lh_466, lh_467, lh_468, \
                         lh_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -ih_465[k]
                   + f_0 * lh_465[k];

        t_466[k] = -ih_466[k]
                   + f_0 * lh_466[k];

        t_467[k] = -ih_467[k]
                   + f_0 * lh_467[k];

        t_468[k] = -ih_468[k]
                   + f_0 * lh_468[k];

        t_469[k] = -ih_469[k]
                   + f_0 * lh_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ih_470, ih_471, ih_472, ih_473, \
                         ih_474, lh_470, lh_471, lh_472, lh_473, \
                         lh_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -ih_470[k]
                   + f_0 * lh_470[k];

        t_471[k] = -ih_471[k]
                   + f_0 * lh_471[k];

        t_472[k] = -ih_472[k]
                   + f_0 * lh_472[k];

        t_473[k] = -ih_473[k]
                   + f_0 * lh_473[k];

        t_474[k] = -ih_474[k]
                   + f_0 * lh_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ih_475, ih_476, ih_477, ih_478, \
                         ih_479, lh_475, lh_476, lh_477, lh_478, \
                         lh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -ih_475[k]
                   + f_0 * lh_475[k];

        t_476[k] = -ih_476[k]
                   + f_0 * lh_476[k];

        t_477[k] = -ih_477[k]
                   + f_0 * lh_477[k];

        t_478[k] = -ih_478[k]
                   + f_0 * lh_478[k];

        t_479[k] = -ih_479[k]
                   + f_0 * lh_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ih_480, ih_481, ih_482, ih_483, \
                         ih_484, lh_480, lh_481, lh_482, lh_483, \
                         lh_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -ih_480[k]
                   + f_0 * lh_480[k];

        t_481[k] = -ih_481[k]
                   + f_0 * lh_481[k];

        t_482[k] = -ih_482[k]
                   + f_0 * lh_482[k];

        t_483[k] = -ih_483[k]
                   + f_0 * lh_483[k];

        t_484[k] = -ih_484[k]
                   + f_0 * lh_484[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ih_485, ih_486, ih_487, ih_488, \
                         ih_489, lh_485, lh_486, lh_487, lh_488, \
                         lh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -ih_485[k]
                   + f_0 * lh_485[k];

        t_486[k] = -ih_486[k]
                   + f_0 * lh_486[k];

        t_487[k] = -ih_487[k]
                   + f_0 * lh_487[k];

        t_488[k] = -ih_488[k]
                   + f_0 * lh_488[k];

        t_489[k] = -ih_489[k]
                   + f_0 * lh_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ih_490, ih_491, ih_492, ih_493, \
                         ih_494, lh_490, lh_491, lh_492, lh_493, \
                         lh_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -ih_490[k]
                   + f_0 * lh_490[k];

        t_491[k] = -ih_491[k]
                   + f_0 * lh_491[k];

        t_492[k] = -ih_492[k]
                   + f_0 * lh_492[k];

        t_493[k] = -ih_493[k]
                   + f_0 * lh_493[k];

        t_494[k] = -ih_494[k]
                   + f_0 * lh_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ih_495, ih_496, ih_497, ih_498, \
                         ih_499, lh_495, lh_496, lh_497, lh_498, \
                         lh_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -ih_495[k]
                   + f_0 * lh_495[k];

        t_496[k] = -ih_496[k]
                   + f_0 * lh_496[k];

        t_497[k] = -ih_497[k]
                   + f_0 * lh_497[k];

        t_498[k] = -ih_498[k]
                   + f_0 * lh_498[k];

        t_499[k] = -ih_499[k]
                   + f_0 * lh_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ih_500, ih_501, ih_502, ih_503, \
                         ih_504, lh_500, lh_501, lh_502, lh_503, \
                         lh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -ih_500[k]
                   + f_0 * lh_500[k];

        t_501[k] = -ih_501[k]
                   + f_0 * lh_501[k];

        t_502[k] = -ih_502[k]
                   + f_0 * lh_502[k];

        t_503[k] = -ih_503[k]
                   + f_0 * lh_503[k];

        t_504[k] = -ih_504[k]
                   + f_0 * lh_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ih_505, ih_506, ih_507, ih_508, \
                         ih_509, lh_505, lh_506, lh_507, lh_508, \
                         lh_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -ih_505[k]
                   + f_0 * lh_505[k];

        t_506[k] = -ih_506[k]
                   + f_0 * lh_506[k];

        t_507[k] = -ih_507[k]
                   + f_0 * lh_507[k];

        t_508[k] = -ih_508[k]
                   + f_0 * lh_508[k];

        t_509[k] = -ih_509[k]
                   + f_0 * lh_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ih_510, ih_511, ih_512, ih_513, \
                         ih_514, lh_510, lh_511, lh_512, lh_513, \
                         lh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -ih_510[k]
                   + f_0 * lh_510[k];

        t_511[k] = -ih_511[k]
                   + f_0 * lh_511[k];

        t_512[k] = -ih_512[k]
                   + f_0 * lh_512[k];

        t_513[k] = -ih_513[k]
                   + f_0 * lh_513[k];

        t_514[k] = -ih_514[k]
                   + f_0 * lh_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ih_515, ih_516, ih_517, ih_518, \
                         ih_519, lh_515, lh_516, lh_517, lh_518, \
                         lh_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -ih_515[k]
                   + f_0 * lh_515[k];

        t_516[k] = -ih_516[k]
                   + f_0 * lh_516[k];

        t_517[k] = -ih_517[k]
                   + f_0 * lh_517[k];

        t_518[k] = -ih_518[k]
                   + f_0 * lh_518[k];

        t_519[k] = -ih_519[k]
                   + f_0 * lh_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ih_520, ih_521, ih_522, ih_523, \
                         ih_524, lh_520, lh_521, lh_522, lh_523, \
                         lh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -ih_520[k]
                   + f_0 * lh_520[k];

        t_521[k] = -ih_521[k]
                   + f_0 * lh_521[k];

        t_522[k] = -ih_522[k]
                   + f_0 * lh_522[k];

        t_523[k] = -ih_523[k]
                   + f_0 * lh_523[k];

        t_524[k] = -ih_524[k]
                   + f_0 * lh_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ih_525, ih_526, ih_527, ih_528, \
                         ih_529, lh_525, lh_526, lh_527, lh_528, \
                         lh_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -ih_525[k]
                   + f_0 * lh_525[k];

        t_526[k] = -ih_526[k]
                   + f_0 * lh_526[k];

        t_527[k] = -ih_527[k]
                   + f_0 * lh_527[k];

        t_528[k] = -ih_528[k]
                   + f_0 * lh_528[k];

        t_529[k] = -ih_529[k]
                   + f_0 * lh_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ih_530, ih_531, ih_532, ih_533, \
                         ih_534, lh_530, lh_531, lh_532, lh_533, \
                         lh_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -ih_530[k]
                   + f_0 * lh_530[k];

        t_531[k] = -ih_531[k]
                   + f_0 * lh_531[k];

        t_532[k] = -ih_532[k]
                   + f_0 * lh_532[k];

        t_533[k] = -ih_533[k]
                   + f_0 * lh_533[k];

        t_534[k] = -ih_534[k]
                   + f_0 * lh_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ih_535, ih_536, ih_537, ih_538, \
                         ih_539, lh_535, lh_536, lh_537, lh_538, \
                         lh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -ih_535[k]
                   + f_0 * lh_535[k];

        t_536[k] = -ih_536[k]
                   + f_0 * lh_536[k];

        t_537[k] = -ih_537[k]
                   + f_0 * lh_537[k];

        t_538[k] = -ih_538[k]
                   + f_0 * lh_538[k];

        t_539[k] = -ih_539[k]
                   + f_0 * lh_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ih_540, ih_541, ih_542, ih_543, \
                         ih_544, lh_540, lh_541, lh_542, lh_543, \
                         lh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -ih_540[k]
                   + f_0 * lh_540[k];

        t_541[k] = -ih_541[k]
                   + f_0 * lh_541[k];

        t_542[k] = -ih_542[k]
                   + f_0 * lh_542[k];

        t_543[k] = -ih_543[k]
                   + f_0 * lh_543[k];

        t_544[k] = -ih_544[k]
                   + f_0 * lh_544[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ih_545, ih_546, ih_547, ih_548, \
                         ih_549, lh_545, lh_546, lh_547, lh_548, \
                         lh_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -ih_545[k]
                   + f_0 * lh_545[k];

        t_546[k] = -ih_546[k]
                   + f_0 * lh_546[k];

        t_547[k] = -ih_547[k]
                   + f_0 * lh_547[k];

        t_548[k] = -ih_548[k]
                   + f_0 * lh_548[k];

        t_549[k] = -ih_549[k]
                   + f_0 * lh_549[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ih_550, ih_551, ih_552, ih_553, \
                         ih_554, lh_550, lh_551, lh_552, lh_553, \
                         lh_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -ih_550[k]
                   + f_0 * lh_550[k];

        t_551[k] = -ih_551[k]
                   + f_0 * lh_551[k];

        t_552[k] = -ih_552[k]
                   + f_0 * lh_552[k];

        t_553[k] = -ih_553[k]
                   + f_0 * lh_553[k];

        t_554[k] = -ih_554[k]
                   + f_0 * lh_554[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ih_555, ih_556, ih_557, ih_558, \
                         ih_559, lh_555, lh_556, lh_557, lh_558, \
                         lh_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -ih_555[k]
                   + f_0 * lh_555[k];

        t_556[k] = -ih_556[k]
                   + f_0 * lh_556[k];

        t_557[k] = -ih_557[k]
                   + f_0 * lh_557[k];

        t_558[k] = -ih_558[k]
                   + f_0 * lh_558[k];

        t_559[k] = -ih_559[k]
                   + f_0 * lh_559[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, ih_560, ih_561, ih_562, ih_563, \
                         ih_564, lh_560, lh_561, lh_562, lh_563, \
                         lh_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = -ih_560[k]
                   + f_0 * lh_560[k];

        t_561[k] = -ih_561[k]
                   + f_0 * lh_561[k];

        t_562[k] = -ih_562[k]
                   + f_0 * lh_562[k];

        t_563[k] = -ih_563[k]
                   + f_0 * lh_563[k];

        t_564[k] = -ih_564[k]
                   + f_0 * lh_564[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, ih_565, ih_566, ih_567, ih_568, \
                         ih_569, lh_565, lh_566, lh_567, lh_568, \
                         lh_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = -ih_565[k]
                   + f_0 * lh_565[k];

        t_566[k] = -ih_566[k]
                   + f_0 * lh_566[k];

        t_567[k] = -ih_567[k]
                   + f_0 * lh_567[k];

        t_568[k] = -ih_568[k]
                   + f_0 * lh_568[k];

        t_569[k] = -ih_569[k]
                   + f_0 * lh_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ih_570, ih_571, ih_572, ih_573, \
                         ih_574, lh_570, lh_571, lh_572, lh_573, \
                         lh_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -ih_570[k]
                   + f_0 * lh_570[k];

        t_571[k] = -ih_571[k]
                   + f_0 * lh_571[k];

        t_572[k] = -ih_572[k]
                   + f_0 * lh_572[k];

        t_573[k] = -ih_573[k]
                   + f_0 * lh_573[k];

        t_574[k] = -ih_574[k]
                   + f_0 * lh_574[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ih_575, ih_576, ih_577, ih_578, \
                         ih_579, lh_575, lh_576, lh_577, lh_578, \
                         lh_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = -ih_575[k]
                   + f_0 * lh_575[k];

        t_576[k] = -ih_576[k]
                   + f_0 * lh_576[k];

        t_577[k] = -ih_577[k]
                   + f_0 * lh_577[k];

        t_578[k] = -ih_578[k]
                   + f_0 * lh_578[k];

        t_579[k] = -ih_579[k]
                   + f_0 * lh_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, ih_580, ih_581, ih_582, ih_583, \
                         ih_584, lh_580, lh_581, lh_582, lh_583, \
                         lh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -ih_580[k]
                   + f_0 * lh_580[k];

        t_581[k] = -ih_581[k]
                   + f_0 * lh_581[k];

        t_582[k] = -ih_582[k]
                   + f_0 * lh_582[k];

        t_583[k] = -ih_583[k]
                   + f_0 * lh_583[k];

        t_584[k] = -ih_584[k]
                   + f_0 * lh_584[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, t_590, ih_585, ih_586, ih_587, \
                         lh_585, lh_586, lh_587, lh_588, lh_589, \
                         lh_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -ih_585[k]
                   + f_0 * lh_585[k];

        t_586[k] = -ih_586[k]
                   + f_0 * lh_586[k];

        t_587[k] = -ih_587[k]
                   + f_0 * lh_587[k];

        t_588[k] = f_0 * lh_588[k];

        t_589[k] = f_0 * lh_589[k];

        t_590[k] = f_0 * lh_590[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, t_595, t_596, t_597, t_598, lh_591, \
                         lh_592, lh_593, lh_594, lh_595, lh_596, lh_597, \
                         lh_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_0 * lh_591[k];

        t_592[k] = f_0 * lh_592[k];

        t_593[k] = f_0 * lh_593[k];

        t_594[k] = f_0 * lh_594[k];

        t_595[k] = f_0 * lh_595[k];

        t_596[k] = f_0 * lh_596[k];

        t_597[k] = f_0 * lh_597[k];

        t_598[k] = f_0 * lh_598[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, t_604, t_605, t_606, lh_599, \
                         lh_600, lh_601, lh_602, lh_603, lh_604, lh_605, \
                         lh_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = f_0 * lh_599[k];

        t_600[k] = f_0 * lh_600[k];

        t_601[k] = f_0 * lh_601[k];

        t_602[k] = f_0 * lh_602[k];

        t_603[k] = f_0 * lh_603[k];

        t_604[k] = f_0 * lh_604[k];

        t_605[k] = f_0 * lh_605[k];

        t_606[k] = f_0 * lh_606[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t lh, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *lh_607 = buffer.data(lh + 607);
    const auto *lh_608 = buffer.data(lh + 608);
    const auto *lh_609 = buffer.data(lh + 609);
    const auto *lh_610 = buffer.data(lh + 610);
    const auto *lh_611 = buffer.data(lh + 611);
    const auto *lh_612 = buffer.data(lh + 612);
    const auto *lh_613 = buffer.data(lh + 613);
    const auto *lh_614 = buffer.data(lh + 614);
    const auto *lh_615 = buffer.data(lh + 615);
    const auto *lh_616 = buffer.data(lh + 616);
    const auto *lh_617 = buffer.data(lh + 617);
    const auto *lh_618 = buffer.data(lh + 618);
    const auto *lh_619 = buffer.data(lh + 619);
    const auto *lh_620 = buffer.data(lh + 620);
    const auto *lh_621 = buffer.data(lh + 621);
    const auto *lh_622 = buffer.data(lh + 622);
    const auto *lh_623 = buffer.data(lh + 623);
    const auto *lh_624 = buffer.data(lh + 624);
    const auto *lh_625 = buffer.data(lh + 625);
    const auto *lh_626 = buffer.data(lh + 626);
    const auto *lh_627 = buffer.data(lh + 627);
    const auto *lh_628 = buffer.data(lh + 628);
    const auto *lh_629 = buffer.data(lh + 629);
    const auto *lh_630 = buffer.data(lh + 630);
    const auto *lh_631 = buffer.data(lh + 631);
    const auto *lh_632 = buffer.data(lh + 632);
    const auto *lh_633 = buffer.data(lh + 633);
    const auto *lh_634 = buffer.data(lh + 634);
    const auto *lh_635 = buffer.data(lh + 635);
    const auto *lh_636 = buffer.data(lh + 636);
    const auto *lh_637 = buffer.data(lh + 637);
    const auto *lh_638 = buffer.data(lh + 638);
    const auto *lh_639 = buffer.data(lh + 639);
    const auto *lh_640 = buffer.data(lh + 640);
    const auto *lh_641 = buffer.data(lh + 641);
    const auto *lh_642 = buffer.data(lh + 642);
    const auto *lh_643 = buffer.data(lh + 643);
    const auto *lh_644 = buffer.data(lh + 644);
    const auto *lh_645 = buffer.data(lh + 645);
    const auto *lh_646 = buffer.data(lh + 646);
    const auto *lh_647 = buffer.data(lh + 647);
    const auto *lh_648 = buffer.data(lh + 648);
    const auto *lh_649 = buffer.data(lh + 649);
    const auto *lh_650 = buffer.data(lh + 650);
    const auto *lh_651 = buffer.data(lh + 651);
    const auto *lh_652 = buffer.data(lh + 652);
    const auto *lh_653 = buffer.data(lh + 653);
    const auto *lh_654 = buffer.data(lh + 654);
    const auto *lh_655 = buffer.data(lh + 655);
    const auto *lh_656 = buffer.data(lh + 656);
    const auto *lh_657 = buffer.data(lh + 657);
    const auto *lh_658 = buffer.data(lh + 658);
    const auto *lh_659 = buffer.data(lh + 659);
    const auto *lh_660 = buffer.data(lh + 660);
    const auto *lh_661 = buffer.data(lh + 661);
    const auto *lh_662 = buffer.data(lh + 662);
    const auto *lh_663 = buffer.data(lh + 663);
    const auto *lh_664 = buffer.data(lh + 664);
    const auto *lh_665 = buffer.data(lh + 665);
    const auto *lh_666 = buffer.data(lh + 666);
    const auto *lh_667 = buffer.data(lh + 667);
    const auto *lh_668 = buffer.data(lh + 668);
    const auto *lh_669 = buffer.data(lh + 669);
    const auto *lh_670 = buffer.data(lh + 670);
    const auto *lh_671 = buffer.data(lh + 671);
    const auto *lh_672 = buffer.data(lh + 672);
    const auto *lh_673 = buffer.data(lh + 673);
    const auto *lh_674 = buffer.data(lh + 674);
    const auto *lh_675 = buffer.data(lh + 675);
    const auto *lh_676 = buffer.data(lh + 676);
    const auto *lh_677 = buffer.data(lh + 677);
    const auto *lh_678 = buffer.data(lh + 678);
    const auto *lh_679 = buffer.data(lh + 679);
    const auto *lh_680 = buffer.data(lh + 680);
    const auto *lh_681 = buffer.data(lh + 681);
    const auto *lh_682 = buffer.data(lh + 682);
    const auto *lh_683 = buffer.data(lh + 683);
    const auto *lh_684 = buffer.data(lh + 684);
    const auto *lh_685 = buffer.data(lh + 685);
    const auto *lh_686 = buffer.data(lh + 686);
    const auto *lh_687 = buffer.data(lh + 687);
    const auto *lh_688 = buffer.data(lh + 688);
    const auto *lh_689 = buffer.data(lh + 689);
    const auto *lh_690 = buffer.data(lh + 690);
    const auto *lh_691 = buffer.data(lh + 691);
    const auto *lh_692 = buffer.data(lh + 692);
    const auto *lh_693 = buffer.data(lh + 693);
    const auto *lh_694 = buffer.data(lh + 694);
    const auto *lh_695 = buffer.data(lh + 695);
    const auto *lh_696 = buffer.data(lh + 696);
    const auto *lh_697 = buffer.data(lh + 697);
    const auto *lh_698 = buffer.data(lh + 698);
    const auto *lh_699 = buffer.data(lh + 699);
    const auto *lh_700 = buffer.data(lh + 700);
    const auto *lh_701 = buffer.data(lh + 701);
    const auto *lh_702 = buffer.data(lh + 702);
    const auto *lh_703 = buffer.data(lh + 703);
    const auto *lh_704 = buffer.data(lh + 704);
    const auto *lh_705 = buffer.data(lh + 705);
    const auto *lh_706 = buffer.data(lh + 706);
    const auto *lh_707 = buffer.data(lh + 707);
    const auto *lh_708 = buffer.data(lh + 708);
    const auto *lh_709 = buffer.data(lh + 709);
    const auto *lh_710 = buffer.data(lh + 710);
    const auto *lh_711 = buffer.data(lh + 711);
    const auto *lh_712 = buffer.data(lh + 712);
    const auto *lh_713 = buffer.data(lh + 713);
    const auto *lh_714 = buffer.data(lh + 714);
    const auto *lh_715 = buffer.data(lh + 715);
    const auto *lh_716 = buffer.data(lh + 716);
    const auto *lh_717 = buffer.data(lh + 717);
    const auto *lh_718 = buffer.data(lh + 718);
    const auto *lh_719 = buffer.data(lh + 719);
    const auto *lh_720 = buffer.data(lh + 720);
    const auto *lh_721 = buffer.data(lh + 721);
    const auto *lh_722 = buffer.data(lh + 722);
    const auto *lh_723 = buffer.data(lh + 723);
    const auto *lh_724 = buffer.data(lh + 724);
    const auto *lh_725 = buffer.data(lh + 725);
    const auto *lh_726 = buffer.data(lh + 726);
    const auto *lh_727 = buffer.data(lh + 727);
    const auto *lh_728 = buffer.data(lh + 728);
    const auto *lh_729 = buffer.data(lh + 729);
    const auto *lh_730 = buffer.data(lh + 730);
    const auto *lh_731 = buffer.data(lh + 731);
    const auto *lh_732 = buffer.data(lh + 732);
    const auto *lh_733 = buffer.data(lh + 733);
    const auto *lh_734 = buffer.data(lh + 734);
    const auto *lh_735 = buffer.data(lh + 735);
    const auto *lh_736 = buffer.data(lh + 736);
    const auto *lh_737 = buffer.data(lh + 737);
    const auto *lh_738 = buffer.data(lh + 738);
    const auto *lh_739 = buffer.data(lh + 739);
    const auto *lh_740 = buffer.data(lh + 740);
    const auto *lh_741 = buffer.data(lh + 741);
    const auto *lh_742 = buffer.data(lh + 742);
    const auto *lh_743 = buffer.data(lh + 743);
    const auto *lh_744 = buffer.data(lh + 744);
    const auto *lh_745 = buffer.data(lh + 745);
    const auto *lh_746 = buffer.data(lh + 746);
    const auto *lh_747 = buffer.data(lh + 747);
    const auto *lh_748 = buffer.data(lh + 748);
    const auto *lh_749 = buffer.data(lh + 749);
    const auto *lh_750 = buffer.data(lh + 750);
    const auto *lh_751 = buffer.data(lh + 751);
    const auto *lh_752 = buffer.data(lh + 752);
    const auto *lh_753 = buffer.data(lh + 753);
    const auto *lh_754 = buffer.data(lh + 754);
    const auto *lh_755 = buffer.data(lh + 755);

#pragma omp simd aligned(t_607, t_608, t_609, t_610, t_611, t_612, t_613, t_614, lh_607, \
                         lh_608, lh_609, lh_610, lh_611, lh_612, lh_613, \
                         lh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_0 * lh_607[k];

        t_608[k] = f_0 * lh_608[k];

        t_609[k] = f_0 * lh_609[k];

        t_610[k] = f_0 * lh_610[k];

        t_611[k] = f_0 * lh_611[k];

        t_612[k] = f_0 * lh_612[k];

        t_613[k] = f_0 * lh_613[k];

        t_614[k] = f_0 * lh_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, t_620, t_621, t_622, lh_615, \
                         lh_616, lh_617, lh_618, lh_619, lh_620, lh_621, \
                         lh_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_0 * lh_615[k];

        t_616[k] = f_0 * lh_616[k];

        t_617[k] = f_0 * lh_617[k];

        t_618[k] = f_0 * lh_618[k];

        t_619[k] = f_0 * lh_619[k];

        t_620[k] = f_0 * lh_620[k];

        t_621[k] = f_0 * lh_621[k];

        t_622[k] = f_0 * lh_622[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, t_626, t_627, t_628, t_629, t_630, lh_623, \
                         lh_624, lh_625, lh_626, lh_627, lh_628, lh_629, \
                         lh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_0 * lh_623[k];

        t_624[k] = f_0 * lh_624[k];

        t_625[k] = f_0 * lh_625[k];

        t_626[k] = f_0 * lh_626[k];

        t_627[k] = f_0 * lh_627[k];

        t_628[k] = f_0 * lh_628[k];

        t_629[k] = f_0 * lh_629[k];

        t_630[k] = f_0 * lh_630[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, t_635, t_636, t_637, t_638, lh_631, \
                         lh_632, lh_633, lh_634, lh_635, lh_636, lh_637, \
                         lh_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = f_0 * lh_631[k];

        t_632[k] = f_0 * lh_632[k];

        t_633[k] = f_0 * lh_633[k];

        t_634[k] = f_0 * lh_634[k];

        t_635[k] = f_0 * lh_635[k];

        t_636[k] = f_0 * lh_636[k];

        t_637[k] = f_0 * lh_637[k];

        t_638[k] = f_0 * lh_638[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, t_643, t_644, t_645, t_646, lh_639, \
                         lh_640, lh_641, lh_642, lh_643, lh_644, lh_645, \
                         lh_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_0 * lh_639[k];

        t_640[k] = f_0 * lh_640[k];

        t_641[k] = f_0 * lh_641[k];

        t_642[k] = f_0 * lh_642[k];

        t_643[k] = f_0 * lh_643[k];

        t_644[k] = f_0 * lh_644[k];

        t_645[k] = f_0 * lh_645[k];

        t_646[k] = f_0 * lh_646[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, t_651, t_652, t_653, t_654, lh_647, \
                         lh_648, lh_649, lh_650, lh_651, lh_652, lh_653, \
                         lh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_0 * lh_647[k];

        t_648[k] = f_0 * lh_648[k];

        t_649[k] = f_0 * lh_649[k];

        t_650[k] = f_0 * lh_650[k];

        t_651[k] = f_0 * lh_651[k];

        t_652[k] = f_0 * lh_652[k];

        t_653[k] = f_0 * lh_653[k];

        t_654[k] = f_0 * lh_654[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, t_660, t_661, t_662, lh_655, \
                         lh_656, lh_657, lh_658, lh_659, lh_660, lh_661, \
                         lh_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_0 * lh_655[k];

        t_656[k] = f_0 * lh_656[k];

        t_657[k] = f_0 * lh_657[k];

        t_658[k] = f_0 * lh_658[k];

        t_659[k] = f_0 * lh_659[k];

        t_660[k] = f_0 * lh_660[k];

        t_661[k] = f_0 * lh_661[k];

        t_662[k] = f_0 * lh_662[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, t_667, t_668, t_669, t_670, lh_663, \
                         lh_664, lh_665, lh_666, lh_667, lh_668, lh_669, \
                         lh_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_0 * lh_663[k];

        t_664[k] = f_0 * lh_664[k];

        t_665[k] = f_0 * lh_665[k];

        t_666[k] = f_0 * lh_666[k];

        t_667[k] = f_0 * lh_667[k];

        t_668[k] = f_0 * lh_668[k];

        t_669[k] = f_0 * lh_669[k];

        t_670[k] = f_0 * lh_670[k];
    }

#pragma omp simd aligned(t_671, t_672, t_673, t_674, t_675, t_676, t_677, t_678, lh_671, \
                         lh_672, lh_673, lh_674, lh_675, lh_676, lh_677, \
                         lh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_671[k] = f_0 * lh_671[k];

        t_672[k] = f_0 * lh_672[k];

        t_673[k] = f_0 * lh_673[k];

        t_674[k] = f_0 * lh_674[k];

        t_675[k] = f_0 * lh_675[k];

        t_676[k] = f_0 * lh_676[k];

        t_677[k] = f_0 * lh_677[k];

        t_678[k] = f_0 * lh_678[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, t_682, t_683, t_684, t_685, t_686, lh_679, \
                         lh_680, lh_681, lh_682, lh_683, lh_684, lh_685, \
                         lh_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_0 * lh_679[k];

        t_680[k] = f_0 * lh_680[k];

        t_681[k] = f_0 * lh_681[k];

        t_682[k] = f_0 * lh_682[k];

        t_683[k] = f_0 * lh_683[k];

        t_684[k] = f_0 * lh_684[k];

        t_685[k] = f_0 * lh_685[k];

        t_686[k] = f_0 * lh_686[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, t_691, t_692, t_693, t_694, lh_687, \
                         lh_688, lh_689, lh_690, lh_691, lh_692, lh_693, \
                         lh_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_0 * lh_687[k];

        t_688[k] = f_0 * lh_688[k];

        t_689[k] = f_0 * lh_689[k];

        t_690[k] = f_0 * lh_690[k];

        t_691[k] = f_0 * lh_691[k];

        t_692[k] = f_0 * lh_692[k];

        t_693[k] = f_0 * lh_693[k];

        t_694[k] = f_0 * lh_694[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, t_700, t_701, t_702, lh_695, \
                         lh_696, lh_697, lh_698, lh_699, lh_700, lh_701, \
                         lh_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_0 * lh_695[k];

        t_696[k] = f_0 * lh_696[k];

        t_697[k] = f_0 * lh_697[k];

        t_698[k] = f_0 * lh_698[k];

        t_699[k] = f_0 * lh_699[k];

        t_700[k] = f_0 * lh_700[k];

        t_701[k] = f_0 * lh_701[k];

        t_702[k] = f_0 * lh_702[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, t_706, t_707, t_708, t_709, t_710, lh_703, \
                         lh_704, lh_705, lh_706, lh_707, lh_708, lh_709, \
                         lh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_0 * lh_703[k];

        t_704[k] = f_0 * lh_704[k];

        t_705[k] = f_0 * lh_705[k];

        t_706[k] = f_0 * lh_706[k];

        t_707[k] = f_0 * lh_707[k];

        t_708[k] = f_0 * lh_708[k];

        t_709[k] = f_0 * lh_709[k];

        t_710[k] = f_0 * lh_710[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, t_714, t_715, t_716, t_717, t_718, lh_711, \
                         lh_712, lh_713, lh_714, lh_715, lh_716, lh_717, \
                         lh_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = f_0 * lh_711[k];

        t_712[k] = f_0 * lh_712[k];

        t_713[k] = f_0 * lh_713[k];

        t_714[k] = f_0 * lh_714[k];

        t_715[k] = f_0 * lh_715[k];

        t_716[k] = f_0 * lh_716[k];

        t_717[k] = f_0 * lh_717[k];

        t_718[k] = f_0 * lh_718[k];
    }

#pragma omp simd aligned(t_719, t_720, t_721, t_722, t_723, t_724, t_725, t_726, lh_719, \
                         lh_720, lh_721, lh_722, lh_723, lh_724, lh_725, \
                         lh_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_719[k] = f_0 * lh_719[k];

        t_720[k] = f_0 * lh_720[k];

        t_721[k] = f_0 * lh_721[k];

        t_722[k] = f_0 * lh_722[k];

        t_723[k] = f_0 * lh_723[k];

        t_724[k] = f_0 * lh_724[k];

        t_725[k] = f_0 * lh_725[k];

        t_726[k] = f_0 * lh_726[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, t_730, t_731, t_732, t_733, t_734, lh_727, \
                         lh_728, lh_729, lh_730, lh_731, lh_732, lh_733, \
                         lh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_0 * lh_727[k];

        t_728[k] = f_0 * lh_728[k];

        t_729[k] = f_0 * lh_729[k];

        t_730[k] = f_0 * lh_730[k];

        t_731[k] = f_0 * lh_731[k];

        t_732[k] = f_0 * lh_732[k];

        t_733[k] = f_0 * lh_733[k];

        t_734[k] = f_0 * lh_734[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, t_740, t_741, t_742, lh_735, \
                         lh_736, lh_737, lh_738, lh_739, lh_740, lh_741, \
                         lh_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_0 * lh_735[k];

        t_736[k] = f_0 * lh_736[k];

        t_737[k] = f_0 * lh_737[k];

        t_738[k] = f_0 * lh_738[k];

        t_739[k] = f_0 * lh_739[k];

        t_740[k] = f_0 * lh_740[k];

        t_741[k] = f_0 * lh_741[k];

        t_742[k] = f_0 * lh_742[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, t_748, t_749, t_750, lh_743, \
                         lh_744, lh_745, lh_746, lh_747, lh_748, lh_749, \
                         lh_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_0 * lh_743[k];

        t_744[k] = f_0 * lh_744[k];

        t_745[k] = f_0 * lh_745[k];

        t_746[k] = f_0 * lh_746[k];

        t_747[k] = f_0 * lh_747[k];

        t_748[k] = f_0 * lh_748[k];

        t_749[k] = f_0 * lh_749[k];

        t_750[k] = f_0 * lh_750[k];
    }

#pragma omp simd aligned(t_751, t_752, t_753, t_754, t_755, lh_751, lh_752, lh_753, lh_754, \
                         lh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_751[k] = f_0 * lh_751[k];

        t_752[k] = f_0 * lh_752[k];

        t_753[k] = f_0 * lh_753[k];

        t_754[k] = f_0 * lh_754[k];

        t_755[k] = f_0 * lh_755[k];
    }
}

auto
compute_prim_geom_10_kh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t ih, const size_t lh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kh_electron_repulsion_0_piece0(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_0_piece1(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_0_piece2(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_0_piece3(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_0_piece4(buffer, target, lh, ncols, alpha);
}

static auto
compute_prim_geom_10_kh_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
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

    const auto *ih_0 = buffer.data(ih + 0);
    const auto *ih_1 = buffer.data(ih + 1);
    const auto *ih_2 = buffer.data(ih + 2);
    const auto *ih_3 = buffer.data(ih + 3);
    const auto *ih_4 = buffer.data(ih + 4);
    const auto *ih_5 = buffer.data(ih + 5);
    const auto *ih_6 = buffer.data(ih + 6);
    const auto *ih_7 = buffer.data(ih + 7);
    const auto *ih_8 = buffer.data(ih + 8);
    const auto *ih_9 = buffer.data(ih + 9);
    const auto *ih_10 = buffer.data(ih + 10);
    const auto *ih_11 = buffer.data(ih + 11);
    const auto *ih_12 = buffer.data(ih + 12);
    const auto *ih_13 = buffer.data(ih + 13);
    const auto *ih_14 = buffer.data(ih + 14);
    const auto *ih_15 = buffer.data(ih + 15);
    const auto *ih_16 = buffer.data(ih + 16);
    const auto *ih_17 = buffer.data(ih + 17);
    const auto *ih_18 = buffer.data(ih + 18);
    const auto *ih_19 = buffer.data(ih + 19);
    const auto *ih_20 = buffer.data(ih + 20);
    const auto *ih_21 = buffer.data(ih + 21);
    const auto *ih_22 = buffer.data(ih + 22);
    const auto *ih_23 = buffer.data(ih + 23);
    const auto *ih_24 = buffer.data(ih + 24);
    const auto *ih_25 = buffer.data(ih + 25);
    const auto *ih_26 = buffer.data(ih + 26);
    const auto *ih_27 = buffer.data(ih + 27);
    const auto *ih_28 = buffer.data(ih + 28);
    const auto *ih_29 = buffer.data(ih + 29);
    const auto *ih_30 = buffer.data(ih + 30);
    const auto *ih_31 = buffer.data(ih + 31);
    const auto *ih_32 = buffer.data(ih + 32);
    const auto *ih_33 = buffer.data(ih + 33);
    const auto *ih_34 = buffer.data(ih + 34);
    const auto *ih_35 = buffer.data(ih + 35);
    const auto *ih_36 = buffer.data(ih + 36);
    const auto *ih_37 = buffer.data(ih + 37);
    const auto *ih_38 = buffer.data(ih + 38);
    const auto *ih_39 = buffer.data(ih + 39);
    const auto *ih_40 = buffer.data(ih + 40);
    const auto *ih_41 = buffer.data(ih + 41);
    const auto *ih_42 = buffer.data(ih + 42);
    const auto *ih_43 = buffer.data(ih + 43);
    const auto *ih_44 = buffer.data(ih + 44);
    const auto *ih_45 = buffer.data(ih + 45);
    const auto *ih_46 = buffer.data(ih + 46);
    const auto *ih_47 = buffer.data(ih + 47);
    const auto *ih_48 = buffer.data(ih + 48);
    const auto *ih_49 = buffer.data(ih + 49);
    const auto *ih_50 = buffer.data(ih + 50);
    const auto *ih_51 = buffer.data(ih + 51);
    const auto *ih_52 = buffer.data(ih + 52);
    const auto *ih_53 = buffer.data(ih + 53);
    const auto *ih_54 = buffer.data(ih + 54);
    const auto *ih_55 = buffer.data(ih + 55);
    const auto *ih_56 = buffer.data(ih + 56);
    const auto *ih_57 = buffer.data(ih + 57);
    const auto *ih_58 = buffer.data(ih + 58);
    const auto *ih_59 = buffer.data(ih + 59);
    const auto *ih_60 = buffer.data(ih + 60);
    const auto *ih_61 = buffer.data(ih + 61);
    const auto *ih_62 = buffer.data(ih + 62);
    const auto *ih_63 = buffer.data(ih + 63);
    const auto *ih_64 = buffer.data(ih + 64);
    const auto *ih_65 = buffer.data(ih + 65);
    const auto *ih_66 = buffer.data(ih + 66);
    const auto *ih_67 = buffer.data(ih + 67);
    const auto *ih_68 = buffer.data(ih + 68);
    const auto *ih_69 = buffer.data(ih + 69);
    const auto *ih_70 = buffer.data(ih + 70);
    const auto *ih_71 = buffer.data(ih + 71);
    const auto *ih_72 = buffer.data(ih + 72);
    const auto *ih_73 = buffer.data(ih + 73);
    const auto *ih_74 = buffer.data(ih + 74);
    const auto *ih_75 = buffer.data(ih + 75);
    const auto *ih_76 = buffer.data(ih + 76);
    const auto *ih_77 = buffer.data(ih + 77);
    const auto *ih_78 = buffer.data(ih + 78);
    const auto *ih_79 = buffer.data(ih + 79);
    const auto *ih_80 = buffer.data(ih + 80);
    const auto *ih_81 = buffer.data(ih + 81);
    const auto *ih_82 = buffer.data(ih + 82);
    const auto *ih_83 = buffer.data(ih + 83);
    const auto *ih_84 = buffer.data(ih + 84);
    const auto *ih_85 = buffer.data(ih + 85);
    const auto *ih_86 = buffer.data(ih + 86);
    const auto *ih_87 = buffer.data(ih + 87);
    const auto *ih_88 = buffer.data(ih + 88);
    const auto *ih_89 = buffer.data(ih + 89);
    const auto *ih_90 = buffer.data(ih + 90);
    const auto *ih_91 = buffer.data(ih + 91);
    const auto *ih_92 = buffer.data(ih + 92);
    const auto *ih_93 = buffer.data(ih + 93);
    const auto *ih_94 = buffer.data(ih + 94);
    const auto *ih_95 = buffer.data(ih + 95);
    const auto *ih_96 = buffer.data(ih + 96);
    const auto *ih_97 = buffer.data(ih + 97);
    const auto *ih_98 = buffer.data(ih + 98);
    const auto *ih_99 = buffer.data(ih + 99);
    const auto *ih_100 = buffer.data(ih + 100);
    const auto *ih_101 = buffer.data(ih + 101);
    const auto *ih_102 = buffer.data(ih + 102);
    const auto *ih_103 = buffer.data(ih + 103);
    const auto *ih_104 = buffer.data(ih + 104);
    const auto *ih_105 = buffer.data(ih + 105);

    const auto *lh_21 = buffer.data(lh + 21);
    const auto *lh_22 = buffer.data(lh + 22);
    const auto *lh_23 = buffer.data(lh + 23);
    const auto *lh_24 = buffer.data(lh + 24);
    const auto *lh_25 = buffer.data(lh + 25);
    const auto *lh_26 = buffer.data(lh + 26);
    const auto *lh_27 = buffer.data(lh + 27);
    const auto *lh_28 = buffer.data(lh + 28);
    const auto *lh_29 = buffer.data(lh + 29);
    const auto *lh_30 = buffer.data(lh + 30);
    const auto *lh_31 = buffer.data(lh + 31);
    const auto *lh_32 = buffer.data(lh + 32);
    const auto *lh_33 = buffer.data(lh + 33);
    const auto *lh_34 = buffer.data(lh + 34);
    const auto *lh_35 = buffer.data(lh + 35);
    const auto *lh_36 = buffer.data(lh + 36);
    const auto *lh_37 = buffer.data(lh + 37);
    const auto *lh_38 = buffer.data(lh + 38);
    const auto *lh_39 = buffer.data(lh + 39);
    const auto *lh_40 = buffer.data(lh + 40);
    const auto *lh_41 = buffer.data(lh + 41);
    const auto *lh_63 = buffer.data(lh + 63);
    const auto *lh_64 = buffer.data(lh + 64);
    const auto *lh_65 = buffer.data(lh + 65);
    const auto *lh_66 = buffer.data(lh + 66);
    const auto *lh_67 = buffer.data(lh + 67);
    const auto *lh_68 = buffer.data(lh + 68);
    const auto *lh_69 = buffer.data(lh + 69);
    const auto *lh_70 = buffer.data(lh + 70);
    const auto *lh_71 = buffer.data(lh + 71);
    const auto *lh_72 = buffer.data(lh + 72);
    const auto *lh_73 = buffer.data(lh + 73);
    const auto *lh_74 = buffer.data(lh + 74);
    const auto *lh_75 = buffer.data(lh + 75);
    const auto *lh_76 = buffer.data(lh + 76);
    const auto *lh_77 = buffer.data(lh + 77);
    const auto *lh_78 = buffer.data(lh + 78);
    const auto *lh_79 = buffer.data(lh + 79);
    const auto *lh_80 = buffer.data(lh + 80);
    const auto *lh_81 = buffer.data(lh + 81);
    const auto *lh_82 = buffer.data(lh + 82);
    const auto *lh_83 = buffer.data(lh + 83);
    const auto *lh_84 = buffer.data(lh + 84);
    const auto *lh_85 = buffer.data(lh + 85);
    const auto *lh_86 = buffer.data(lh + 86);
    const auto *lh_87 = buffer.data(lh + 87);
    const auto *lh_88 = buffer.data(lh + 88);
    const auto *lh_89 = buffer.data(lh + 89);
    const auto *lh_90 = buffer.data(lh + 90);
    const auto *lh_91 = buffer.data(lh + 91);
    const auto *lh_92 = buffer.data(lh + 92);
    const auto *lh_93 = buffer.data(lh + 93);
    const auto *lh_94 = buffer.data(lh + 94);
    const auto *lh_95 = buffer.data(lh + 95);
    const auto *lh_96 = buffer.data(lh + 96);
    const auto *lh_97 = buffer.data(lh + 97);
    const auto *lh_98 = buffer.data(lh + 98);
    const auto *lh_99 = buffer.data(lh + 99);
    const auto *lh_100 = buffer.data(lh + 100);
    const auto *lh_101 = buffer.data(lh + 101);
    const auto *lh_102 = buffer.data(lh + 102);
    const auto *lh_103 = buffer.data(lh + 103);
    const auto *lh_104 = buffer.data(lh + 104);
    const auto *lh_126 = buffer.data(lh + 126);
    const auto *lh_127 = buffer.data(lh + 127);
    const auto *lh_128 = buffer.data(lh + 128);
    const auto *lh_129 = buffer.data(lh + 129);
    const auto *lh_130 = buffer.data(lh + 130);
    const auto *lh_131 = buffer.data(lh + 131);
    const auto *lh_132 = buffer.data(lh + 132);
    const auto *lh_133 = buffer.data(lh + 133);
    const auto *lh_134 = buffer.data(lh + 134);
    const auto *lh_135 = buffer.data(lh + 135);
    const auto *lh_136 = buffer.data(lh + 136);
    const auto *lh_137 = buffer.data(lh + 137);
    const auto *lh_138 = buffer.data(lh + 138);
    const auto *lh_139 = buffer.data(lh + 139);
    const auto *lh_140 = buffer.data(lh + 140);
    const auto *lh_141 = buffer.data(lh + 141);
    const auto *lh_142 = buffer.data(lh + 142);
    const auto *lh_143 = buffer.data(lh + 143);
    const auto *lh_144 = buffer.data(lh + 144);
    const auto *lh_145 = buffer.data(lh + 145);
    const auto *lh_146 = buffer.data(lh + 146);
    const auto *lh_147 = buffer.data(lh + 147);
    const auto *lh_148 = buffer.data(lh + 148);
    const auto *lh_149 = buffer.data(lh + 149);
    const auto *lh_150 = buffer.data(lh + 150);
    const auto *lh_151 = buffer.data(lh + 151);
    const auto *lh_152 = buffer.data(lh + 152);
    const auto *lh_153 = buffer.data(lh + 153);
    const auto *lh_154 = buffer.data(lh + 154);
    const auto *lh_155 = buffer.data(lh + 155);
    const auto *lh_156 = buffer.data(lh + 156);
    const auto *lh_157 = buffer.data(lh + 157);
    const auto *lh_158 = buffer.data(lh + 158);
    const auto *lh_159 = buffer.data(lh + 159);
    const auto *lh_160 = buffer.data(lh + 160);
    const auto *lh_161 = buffer.data(lh + 161);
    const auto *lh_162 = buffer.data(lh + 162);
    const auto *lh_163 = buffer.data(lh + 163);
    const auto *lh_164 = buffer.data(lh + 164);
    const auto *lh_165 = buffer.data(lh + 165);
    const auto *lh_166 = buffer.data(lh + 166);
    const auto *lh_167 = buffer.data(lh + 167);
    const auto *lh_168 = buffer.data(lh + 168);
    const auto *lh_169 = buffer.data(lh + 169);
    const auto *lh_170 = buffer.data(lh + 170);
    const auto *lh_171 = buffer.data(lh + 171);
    const auto *lh_172 = buffer.data(lh + 172);
    const auto *lh_173 = buffer.data(lh + 173);
    const auto *lh_174 = buffer.data(lh + 174);
    const auto *lh_175 = buffer.data(lh + 175);
    const auto *lh_176 = buffer.data(lh + 176);
    const auto *lh_177 = buffer.data(lh + 177);
    const auto *lh_178 = buffer.data(lh + 178);
    const auto *lh_179 = buffer.data(lh + 179);
    const auto *lh_180 = buffer.data(lh + 180);
    const auto *lh_181 = buffer.data(lh + 181);
    const auto *lh_182 = buffer.data(lh + 182);
    const auto *lh_183 = buffer.data(lh + 183);
    const auto *lh_184 = buffer.data(lh + 184);
    const auto *lh_185 = buffer.data(lh + 185);
    const auto *lh_186 = buffer.data(lh + 186);
    const auto *lh_187 = buffer.data(lh + 187);
    const auto *lh_188 = buffer.data(lh + 188);
    const auto *lh_210 = buffer.data(lh + 210);
    const auto *lh_211 = buffer.data(lh + 211);
    const auto *lh_212 = buffer.data(lh + 212);
    const auto *lh_213 = buffer.data(lh + 213);
    const auto *lh_214 = buffer.data(lh + 214);
    const auto *lh_215 = buffer.data(lh + 215);
    const auto *lh_216 = buffer.data(lh + 216);
    const auto *lh_217 = buffer.data(lh + 217);
    const auto *lh_218 = buffer.data(lh + 218);
    const auto *lh_219 = buffer.data(lh + 219);
    const auto *lh_220 = buffer.data(lh + 220);
    const auto *lh_221 = buffer.data(lh + 221);
    const auto *lh_222 = buffer.data(lh + 222);
    const auto *lh_223 = buffer.data(lh + 223);
    const auto *lh_224 = buffer.data(lh + 224);
    const auto *lh_225 = buffer.data(lh + 225);
    const auto *lh_226 = buffer.data(lh + 226);
    const auto *lh_227 = buffer.data(lh + 227);
    const auto *lh_228 = buffer.data(lh + 228);
    const auto *lh_229 = buffer.data(lh + 229);
    const auto *lh_230 = buffer.data(lh + 230);
    const auto *lh_231 = buffer.data(lh + 231);
    const auto *lh_232 = buffer.data(lh + 232);
    const auto *lh_233 = buffer.data(lh + 233);
    const auto *lh_234 = buffer.data(lh + 234);
    const auto *lh_235 = buffer.data(lh + 235);
    const auto *lh_236 = buffer.data(lh + 236);
    const auto *lh_237 = buffer.data(lh + 237);
    const auto *lh_238 = buffer.data(lh + 238);
    const auto *lh_239 = buffer.data(lh + 239);
    const auto *lh_240 = buffer.data(lh + 240);
    const auto *lh_241 = buffer.data(lh + 241);
    const auto *lh_242 = buffer.data(lh + 242);
    const auto *lh_243 = buffer.data(lh + 243);
    const auto *lh_244 = buffer.data(lh + 244);
    const auto *lh_245 = buffer.data(lh + 245);
    const auto *lh_246 = buffer.data(lh + 246);
    const auto *lh_247 = buffer.data(lh + 247);
    const auto *lh_248 = buffer.data(lh + 248);
    const auto *lh_249 = buffer.data(lh + 249);
    const auto *lh_250 = buffer.data(lh + 250);
    const auto *lh_251 = buffer.data(lh + 251);
    const auto *lh_252 = buffer.data(lh + 252);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, lh_21, lh_22, lh_23, lh_24, \
                         lh_25, lh_26, lh_27, lh_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lh_21[k];

        t_1[k] = f_0 * lh_22[k];

        t_2[k] = f_0 * lh_23[k];

        t_3[k] = f_0 * lh_24[k];

        t_4[k] = f_0 * lh_25[k];

        t_5[k] = f_0 * lh_26[k];

        t_6[k] = f_0 * lh_27[k];

        t_7[k] = f_0 * lh_28[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, lh_29, lh_30, lh_31, \
                         lh_32, lh_33, lh_34, lh_35, lh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * lh_29[k];

        t_9[k] = f_0 * lh_30[k];

        t_10[k] = f_0 * lh_31[k];

        t_11[k] = f_0 * lh_32[k];

        t_12[k] = f_0 * lh_33[k];

        t_13[k] = f_0 * lh_34[k];

        t_14[k] = f_0 * lh_35[k];

        t_15[k] = f_0 * lh_36[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, ih_0, ih_1, lh_37, lh_38, \
                         lh_39, lh_40, lh_41, lh_63, lh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * lh_37[k];

        t_17[k] = f_0 * lh_38[k];

        t_18[k] = f_0 * lh_39[k];

        t_19[k] = f_0 * lh_40[k];

        t_20[k] = f_0 * lh_41[k];

        t_21[k] = -ih_0[k]
                  + f_0 * lh_63[k];

        t_22[k] = -ih_1[k]
                  + f_0 * lh_64[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, ih_2, ih_3, ih_4, ih_5, ih_6, lh_65, \
                         lh_66, lh_67, lh_68, lh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -ih_2[k]
                  + f_0 * lh_65[k];

        t_24[k] = -ih_3[k]
                  + f_0 * lh_66[k];

        t_25[k] = -ih_4[k]
                  + f_0 * lh_67[k];

        t_26[k] = -ih_5[k]
                  + f_0 * lh_68[k];

        t_27[k] = -ih_6[k]
                  + f_0 * lh_69[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, ih_7, ih_8, ih_9, ih_10, ih_11, lh_70, \
                         lh_71, lh_72, lh_73, lh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -ih_7[k]
                  + f_0 * lh_70[k];

        t_29[k] = -ih_8[k]
                  + f_0 * lh_71[k];

        t_30[k] = -ih_9[k]
                  + f_0 * lh_72[k];

        t_31[k] = -ih_10[k]
                  + f_0 * lh_73[k];

        t_32[k] = -ih_11[k]
                  + f_0 * lh_74[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, ih_12, ih_13, ih_14, ih_15, ih_16, \
                         lh_75, lh_76, lh_77, lh_78, lh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -ih_12[k]
                  + f_0 * lh_75[k];

        t_34[k] = -ih_13[k]
                  + f_0 * lh_76[k];

        t_35[k] = -ih_14[k]
                  + f_0 * lh_77[k];

        t_36[k] = -ih_15[k]
                  + f_0 * lh_78[k];

        t_37[k] = -ih_16[k]
                  + f_0 * lh_79[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ih_17, ih_18, ih_19, ih_20, \
                         lh_80, lh_81, lh_82, lh_83, lh_84, lh_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -ih_17[k]
                  + f_0 * lh_80[k];

        t_39[k] = -ih_18[k]
                  + f_0 * lh_81[k];

        t_40[k] = -ih_19[k]
                  + f_0 * lh_82[k];

        t_41[k] = -ih_20[k]
                  + f_0 * lh_83[k];

        t_42[k] = f_0 * lh_84[k];

        t_43[k] = f_0 * lh_85[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, t_51, lh_86, lh_87, lh_88, \
                         lh_89, lh_90, lh_91, lh_92, lh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * lh_86[k];

        t_45[k] = f_0 * lh_87[k];

        t_46[k] = f_0 * lh_88[k];

        t_47[k] = f_0 * lh_89[k];

        t_48[k] = f_0 * lh_90[k];

        t_49[k] = f_0 * lh_91[k];

        t_50[k] = f_0 * lh_92[k];

        t_51[k] = f_0 * lh_93[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, t_58, t_59, lh_94, lh_95, lh_96, \
                         lh_97, lh_98, lh_99, lh_100, lh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * lh_94[k];

        t_53[k] = f_0 * lh_95[k];

        t_54[k] = f_0 * lh_96[k];

        t_55[k] = f_0 * lh_97[k];

        t_56[k] = f_0 * lh_98[k];

        t_57[k] = f_0 * lh_99[k];

        t_58[k] = f_0 * lh_100[k];

        t_59[k] = f_0 * lh_101[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, ih_21, ih_22, ih_23, lh_102, \
                         lh_103, lh_104, lh_126, lh_127, lh_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * lh_102[k];

        t_61[k] = f_0 * lh_103[k];

        t_62[k] = f_0 * lh_104[k];

        t_63[k] = -2.0 * ih_21[k]
                  + f_0 * lh_126[k];

        t_64[k] = -2.0 * ih_22[k]
                  + f_0 * lh_127[k];

        t_65[k] = -2.0 * ih_23[k]
                  + f_0 * lh_128[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, ih_24, ih_25, ih_26, ih_27, ih_28, \
                         lh_129, lh_130, lh_131, lh_132, lh_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -2.0 * ih_24[k]
                  + f_0 * lh_129[k];

        t_67[k] = -2.0 * ih_25[k]
                  + f_0 * lh_130[k];

        t_68[k] = -2.0 * ih_26[k]
                  + f_0 * lh_131[k];

        t_69[k] = -2.0 * ih_27[k]
                  + f_0 * lh_132[k];

        t_70[k] = -2.0 * ih_28[k]
                  + f_0 * lh_133[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, ih_29, ih_30, ih_31, ih_32, ih_33, \
                         lh_134, lh_135, lh_136, lh_137, lh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -2.0 * ih_29[k]
                  + f_0 * lh_134[k];

        t_72[k] = -2.0 * ih_30[k]
                  + f_0 * lh_135[k];

        t_73[k] = -2.0 * ih_31[k]
                  + f_0 * lh_136[k];

        t_74[k] = -2.0 * ih_32[k]
                  + f_0 * lh_137[k];

        t_75[k] = -2.0 * ih_33[k]
                  + f_0 * lh_138[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, ih_34, ih_35, ih_36, ih_37, ih_38, \
                         lh_139, lh_140, lh_141, lh_142, lh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -2.0 * ih_34[k]
                  + f_0 * lh_139[k];

        t_77[k] = -2.0 * ih_35[k]
                  + f_0 * lh_140[k];

        t_78[k] = -2.0 * ih_36[k]
                  + f_0 * lh_141[k];

        t_79[k] = -2.0 * ih_37[k]
                  + f_0 * lh_142[k];

        t_80[k] = -2.0 * ih_38[k]
                  + f_0 * lh_143[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, ih_39, ih_40, ih_41, ih_42, ih_43, \
                         lh_144, lh_145, lh_146, lh_147, lh_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -2.0 * ih_39[k]
                  + f_0 * lh_144[k];

        t_82[k] = -2.0 * ih_40[k]
                  + f_0 * lh_145[k];

        t_83[k] = -2.0 * ih_41[k]
                  + f_0 * lh_146[k];

        t_84[k] = -ih_42[k]
                  + f_0 * lh_147[k];

        t_85[k] = -ih_43[k]
                  + f_0 * lh_148[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, ih_44, ih_45, ih_46, ih_47, ih_48, \
                         lh_149, lh_150, lh_151, lh_152, lh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -ih_44[k]
                  + f_0 * lh_149[k];

        t_87[k] = -ih_45[k]
                  + f_0 * lh_150[k];

        t_88[k] = -ih_46[k]
                  + f_0 * lh_151[k];

        t_89[k] = -ih_47[k]
                  + f_0 * lh_152[k];

        t_90[k] = -ih_48[k]
                  + f_0 * lh_153[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, ih_49, ih_50, ih_51, ih_52, ih_53, \
                         lh_154, lh_155, lh_156, lh_157, lh_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = -ih_49[k]
                  + f_0 * lh_154[k];

        t_92[k] = -ih_50[k]
                  + f_0 * lh_155[k];

        t_93[k] = -ih_51[k]
                  + f_0 * lh_156[k];

        t_94[k] = -ih_52[k]
                  + f_0 * lh_157[k];

        t_95[k] = -ih_53[k]
                  + f_0 * lh_158[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, ih_54, ih_55, ih_56, ih_57, ih_58, \
                         lh_159, lh_160, lh_161, lh_162, lh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = -ih_54[k]
                  + f_0 * lh_159[k];

        t_97[k] = -ih_55[k]
                  + f_0 * lh_160[k];

        t_98[k] = -ih_56[k]
                  + f_0 * lh_161[k];

        t_99[k] = -ih_57[k]
                  + f_0 * lh_162[k];

        t_100[k] = -ih_58[k]
                   + f_0 * lh_163[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, ih_59, ih_60, ih_61, ih_62, \
                         lh_164, lh_165, lh_166, lh_167, lh_168, \
                         lh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = -ih_59[k]
                   + f_0 * lh_164[k];

        t_102[k] = -ih_60[k]
                   + f_0 * lh_165[k];

        t_103[k] = -ih_61[k]
                   + f_0 * lh_166[k];

        t_104[k] = -ih_62[k]
                   + f_0 * lh_167[k];

        t_105[k] = f_0 * lh_168[k];

        t_106[k] = f_0 * lh_169[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, t_114, lh_170, \
                         lh_171, lh_172, lh_173, lh_174, lh_175, lh_176, \
                         lh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * lh_170[k];

        t_108[k] = f_0 * lh_171[k];

        t_109[k] = f_0 * lh_172[k];

        t_110[k] = f_0 * lh_173[k];

        t_111[k] = f_0 * lh_174[k];

        t_112[k] = f_0 * lh_175[k];

        t_113[k] = f_0 * lh_176[k];

        t_114[k] = f_0 * lh_177[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, t_121, t_122, lh_178, \
                         lh_179, lh_180, lh_181, lh_182, lh_183, lh_184, \
                         lh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_0 * lh_178[k];

        t_116[k] = f_0 * lh_179[k];

        t_117[k] = f_0 * lh_180[k];

        t_118[k] = f_0 * lh_181[k];

        t_119[k] = f_0 * lh_182[k];

        t_120[k] = f_0 * lh_183[k];

        t_121[k] = f_0 * lh_184[k];

        t_122[k] = f_0 * lh_185[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, t_127, t_128, ih_63, ih_64, ih_65, \
                         lh_186, lh_187, lh_188, lh_210, lh_211, \
                         lh_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_0 * lh_186[k];

        t_124[k] = f_0 * lh_187[k];

        t_125[k] = f_0 * lh_188[k];

        t_126[k] = -3.0 * ih_63[k]
                   + f_0 * lh_210[k];

        t_127[k] = -3.0 * ih_64[k]
                   + f_0 * lh_211[k];

        t_128[k] = -3.0 * ih_65[k]
                   + f_0 * lh_212[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, ih_66, ih_67, ih_68, ih_69, ih_70, \
                         lh_213, lh_214, lh_215, lh_216, lh_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -3.0 * ih_66[k]
                   + f_0 * lh_213[k];

        t_130[k] = -3.0 * ih_67[k]
                   + f_0 * lh_214[k];

        t_131[k] = -3.0 * ih_68[k]
                   + f_0 * lh_215[k];

        t_132[k] = -3.0 * ih_69[k]
                   + f_0 * lh_216[k];

        t_133[k] = -3.0 * ih_70[k]
                   + f_0 * lh_217[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ih_71, ih_72, ih_73, ih_74, ih_75, \
                         lh_218, lh_219, lh_220, lh_221, lh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -3.0 * ih_71[k]
                   + f_0 * lh_218[k];

        t_135[k] = -3.0 * ih_72[k]
                   + f_0 * lh_219[k];

        t_136[k] = -3.0 * ih_73[k]
                   + f_0 * lh_220[k];

        t_137[k] = -3.0 * ih_74[k]
                   + f_0 * lh_221[k];

        t_138[k] = -3.0 * ih_75[k]
                   + f_0 * lh_222[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ih_76, ih_77, ih_78, ih_79, ih_80, \
                         lh_223, lh_224, lh_225, lh_226, lh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -3.0 * ih_76[k]
                   + f_0 * lh_223[k];

        t_140[k] = -3.0 * ih_77[k]
                   + f_0 * lh_224[k];

        t_141[k] = -3.0 * ih_78[k]
                   + f_0 * lh_225[k];

        t_142[k] = -3.0 * ih_79[k]
                   + f_0 * lh_226[k];

        t_143[k] = -3.0 * ih_80[k]
                   + f_0 * lh_227[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, ih_81, ih_82, ih_83, ih_84, ih_85, \
                         lh_228, lh_229, lh_230, lh_231, lh_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -3.0 * ih_81[k]
                   + f_0 * lh_228[k];

        t_145[k] = -3.0 * ih_82[k]
                   + f_0 * lh_229[k];

        t_146[k] = -3.0 * ih_83[k]
                   + f_0 * lh_230[k];

        t_147[k] = -2.0 * ih_84[k]
                   + f_0 * lh_231[k];

        t_148[k] = -2.0 * ih_85[k]
                   + f_0 * lh_232[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, ih_86, ih_87, ih_88, ih_89, ih_90, \
                         lh_233, lh_234, lh_235, lh_236, lh_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -2.0 * ih_86[k]
                   + f_0 * lh_233[k];

        t_150[k] = -2.0 * ih_87[k]
                   + f_0 * lh_234[k];

        t_151[k] = -2.0 * ih_88[k]
                   + f_0 * lh_235[k];

        t_152[k] = -2.0 * ih_89[k]
                   + f_0 * lh_236[k];

        t_153[k] = -2.0 * ih_90[k]
                   + f_0 * lh_237[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, ih_91, ih_92, ih_93, ih_94, ih_95, \
                         lh_238, lh_239, lh_240, lh_241, lh_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -2.0 * ih_91[k]
                   + f_0 * lh_238[k];

        t_155[k] = -2.0 * ih_92[k]
                   + f_0 * lh_239[k];

        t_156[k] = -2.0 * ih_93[k]
                   + f_0 * lh_240[k];

        t_157[k] = -2.0 * ih_94[k]
                   + f_0 * lh_241[k];

        t_158[k] = -2.0 * ih_95[k]
                   + f_0 * lh_242[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, ih_96, ih_97, ih_98, ih_99, \
                         ih_100, lh_243, lh_244, lh_245, lh_246, \
                         lh_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -2.0 * ih_96[k]
                   + f_0 * lh_243[k];

        t_160[k] = -2.0 * ih_97[k]
                   + f_0 * lh_244[k];

        t_161[k] = -2.0 * ih_98[k]
                   + f_0 * lh_245[k];

        t_162[k] = -2.0 * ih_99[k]
                   + f_0 * lh_246[k];

        t_163[k] = -2.0 * ih_100[k]
                   + f_0 * lh_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, ih_101, ih_102, ih_103, ih_104, \
                         ih_105, lh_248, lh_249, lh_250, lh_251, \
                         lh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -2.0 * ih_101[k]
                   + f_0 * lh_248[k];

        t_165[k] = -2.0 * ih_102[k]
                   + f_0 * lh_249[k];

        t_166[k] = -2.0 * ih_103[k]
                   + f_0 * lh_250[k];

        t_167[k] = -2.0 * ih_104[k]
                   + f_0 * lh_251[k];

        t_168[k] = -ih_105[k]
                   + f_0 * lh_252[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
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

    const auto *ih_106 = buffer.data(ih + 106);
    const auto *ih_107 = buffer.data(ih + 107);
    const auto *ih_108 = buffer.data(ih + 108);
    const auto *ih_109 = buffer.data(ih + 109);
    const auto *ih_110 = buffer.data(ih + 110);
    const auto *ih_111 = buffer.data(ih + 111);
    const auto *ih_112 = buffer.data(ih + 112);
    const auto *ih_113 = buffer.data(ih + 113);
    const auto *ih_114 = buffer.data(ih + 114);
    const auto *ih_115 = buffer.data(ih + 115);
    const auto *ih_116 = buffer.data(ih + 116);
    const auto *ih_117 = buffer.data(ih + 117);
    const auto *ih_118 = buffer.data(ih + 118);
    const auto *ih_119 = buffer.data(ih + 119);
    const auto *ih_120 = buffer.data(ih + 120);
    const auto *ih_121 = buffer.data(ih + 121);
    const auto *ih_122 = buffer.data(ih + 122);
    const auto *ih_123 = buffer.data(ih + 123);
    const auto *ih_124 = buffer.data(ih + 124);
    const auto *ih_125 = buffer.data(ih + 125);
    const auto *ih_126 = buffer.data(ih + 126);
    const auto *ih_127 = buffer.data(ih + 127);
    const auto *ih_128 = buffer.data(ih + 128);
    const auto *ih_129 = buffer.data(ih + 129);
    const auto *ih_130 = buffer.data(ih + 130);
    const auto *ih_131 = buffer.data(ih + 131);
    const auto *ih_132 = buffer.data(ih + 132);
    const auto *ih_133 = buffer.data(ih + 133);
    const auto *ih_134 = buffer.data(ih + 134);
    const auto *ih_135 = buffer.data(ih + 135);
    const auto *ih_136 = buffer.data(ih + 136);
    const auto *ih_137 = buffer.data(ih + 137);
    const auto *ih_138 = buffer.data(ih + 138);
    const auto *ih_139 = buffer.data(ih + 139);
    const auto *ih_140 = buffer.data(ih + 140);
    const auto *ih_141 = buffer.data(ih + 141);
    const auto *ih_142 = buffer.data(ih + 142);
    const auto *ih_143 = buffer.data(ih + 143);
    const auto *ih_144 = buffer.data(ih + 144);
    const auto *ih_145 = buffer.data(ih + 145);
    const auto *ih_146 = buffer.data(ih + 146);
    const auto *ih_147 = buffer.data(ih + 147);
    const auto *ih_148 = buffer.data(ih + 148);
    const auto *ih_149 = buffer.data(ih + 149);
    const auto *ih_150 = buffer.data(ih + 150);
    const auto *ih_151 = buffer.data(ih + 151);
    const auto *ih_152 = buffer.data(ih + 152);
    const auto *ih_153 = buffer.data(ih + 153);
    const auto *ih_154 = buffer.data(ih + 154);
    const auto *ih_155 = buffer.data(ih + 155);
    const auto *ih_156 = buffer.data(ih + 156);
    const auto *ih_157 = buffer.data(ih + 157);
    const auto *ih_158 = buffer.data(ih + 158);
    const auto *ih_159 = buffer.data(ih + 159);
    const auto *ih_160 = buffer.data(ih + 160);
    const auto *ih_161 = buffer.data(ih + 161);
    const auto *ih_162 = buffer.data(ih + 162);
    const auto *ih_163 = buffer.data(ih + 163);
    const auto *ih_164 = buffer.data(ih + 164);
    const auto *ih_165 = buffer.data(ih + 165);
    const auto *ih_166 = buffer.data(ih + 166);
    const auto *ih_167 = buffer.data(ih + 167);
    const auto *ih_168 = buffer.data(ih + 168);
    const auto *ih_169 = buffer.data(ih + 169);
    const auto *ih_170 = buffer.data(ih + 170);
    const auto *ih_171 = buffer.data(ih + 171);
    const auto *ih_172 = buffer.data(ih + 172);
    const auto *ih_173 = buffer.data(ih + 173);
    const auto *ih_174 = buffer.data(ih + 174);
    const auto *ih_175 = buffer.data(ih + 175);
    const auto *ih_176 = buffer.data(ih + 176);
    const auto *ih_177 = buffer.data(ih + 177);
    const auto *ih_178 = buffer.data(ih + 178);
    const auto *ih_179 = buffer.data(ih + 179);
    const auto *ih_180 = buffer.data(ih + 180);
    const auto *ih_181 = buffer.data(ih + 181);
    const auto *ih_182 = buffer.data(ih + 182);
    const auto *ih_183 = buffer.data(ih + 183);
    const auto *ih_184 = buffer.data(ih + 184);
    const auto *ih_185 = buffer.data(ih + 185);
    const auto *ih_186 = buffer.data(ih + 186);
    const auto *ih_187 = buffer.data(ih + 187);
    const auto *ih_188 = buffer.data(ih + 188);
    const auto *ih_189 = buffer.data(ih + 189);
    const auto *ih_190 = buffer.data(ih + 190);
    const auto *ih_191 = buffer.data(ih + 191);
    const auto *ih_192 = buffer.data(ih + 192);
    const auto *ih_193 = buffer.data(ih + 193);
    const auto *ih_194 = buffer.data(ih + 194);
    const auto *ih_195 = buffer.data(ih + 195);
    const auto *ih_196 = buffer.data(ih + 196);
    const auto *ih_197 = buffer.data(ih + 197);
    const auto *ih_198 = buffer.data(ih + 198);
    const auto *ih_199 = buffer.data(ih + 199);
    const auto *ih_200 = buffer.data(ih + 200);
    const auto *ih_201 = buffer.data(ih + 201);
    const auto *ih_202 = buffer.data(ih + 202);
    const auto *ih_203 = buffer.data(ih + 203);
    const auto *ih_204 = buffer.data(ih + 204);
    const auto *ih_205 = buffer.data(ih + 205);
    const auto *ih_206 = buffer.data(ih + 206);
    const auto *ih_207 = buffer.data(ih + 207);
    const auto *ih_208 = buffer.data(ih + 208);
    const auto *ih_209 = buffer.data(ih + 209);
    const auto *ih_210 = buffer.data(ih + 210);
    const auto *ih_211 = buffer.data(ih + 211);
    const auto *ih_212 = buffer.data(ih + 212);
    const auto *ih_213 = buffer.data(ih + 213);
    const auto *ih_214 = buffer.data(ih + 214);
    const auto *ih_215 = buffer.data(ih + 215);
    const auto *ih_216 = buffer.data(ih + 216);
    const auto *ih_217 = buffer.data(ih + 217);
    const auto *ih_218 = buffer.data(ih + 218);
    const auto *ih_219 = buffer.data(ih + 219);
    const auto *ih_220 = buffer.data(ih + 220);
    const auto *ih_221 = buffer.data(ih + 221);
    const auto *ih_222 = buffer.data(ih + 222);
    const auto *ih_223 = buffer.data(ih + 223);
    const auto *ih_224 = buffer.data(ih + 224);

    const auto *lh_253 = buffer.data(lh + 253);
    const auto *lh_254 = buffer.data(lh + 254);
    const auto *lh_255 = buffer.data(lh + 255);
    const auto *lh_256 = buffer.data(lh + 256);
    const auto *lh_257 = buffer.data(lh + 257);
    const auto *lh_258 = buffer.data(lh + 258);
    const auto *lh_259 = buffer.data(lh + 259);
    const auto *lh_260 = buffer.data(lh + 260);
    const auto *lh_261 = buffer.data(lh + 261);
    const auto *lh_262 = buffer.data(lh + 262);
    const auto *lh_263 = buffer.data(lh + 263);
    const auto *lh_264 = buffer.data(lh + 264);
    const auto *lh_265 = buffer.data(lh + 265);
    const auto *lh_266 = buffer.data(lh + 266);
    const auto *lh_267 = buffer.data(lh + 267);
    const auto *lh_268 = buffer.data(lh + 268);
    const auto *lh_269 = buffer.data(lh + 269);
    const auto *lh_270 = buffer.data(lh + 270);
    const auto *lh_271 = buffer.data(lh + 271);
    const auto *lh_272 = buffer.data(lh + 272);
    const auto *lh_273 = buffer.data(lh + 273);
    const auto *lh_274 = buffer.data(lh + 274);
    const auto *lh_275 = buffer.data(lh + 275);
    const auto *lh_276 = buffer.data(lh + 276);
    const auto *lh_277 = buffer.data(lh + 277);
    const auto *lh_278 = buffer.data(lh + 278);
    const auto *lh_279 = buffer.data(lh + 279);
    const auto *lh_280 = buffer.data(lh + 280);
    const auto *lh_281 = buffer.data(lh + 281);
    const auto *lh_282 = buffer.data(lh + 282);
    const auto *lh_283 = buffer.data(lh + 283);
    const auto *lh_284 = buffer.data(lh + 284);
    const auto *lh_285 = buffer.data(lh + 285);
    const auto *lh_286 = buffer.data(lh + 286);
    const auto *lh_287 = buffer.data(lh + 287);
    const auto *lh_288 = buffer.data(lh + 288);
    const auto *lh_289 = buffer.data(lh + 289);
    const auto *lh_290 = buffer.data(lh + 290);
    const auto *lh_291 = buffer.data(lh + 291);
    const auto *lh_292 = buffer.data(lh + 292);
    const auto *lh_293 = buffer.data(lh + 293);
    const auto *lh_315 = buffer.data(lh + 315);
    const auto *lh_316 = buffer.data(lh + 316);
    const auto *lh_317 = buffer.data(lh + 317);
    const auto *lh_318 = buffer.data(lh + 318);
    const auto *lh_319 = buffer.data(lh + 319);
    const auto *lh_320 = buffer.data(lh + 320);
    const auto *lh_321 = buffer.data(lh + 321);
    const auto *lh_322 = buffer.data(lh + 322);
    const auto *lh_323 = buffer.data(lh + 323);
    const auto *lh_324 = buffer.data(lh + 324);
    const auto *lh_325 = buffer.data(lh + 325);
    const auto *lh_326 = buffer.data(lh + 326);
    const auto *lh_327 = buffer.data(lh + 327);
    const auto *lh_328 = buffer.data(lh + 328);
    const auto *lh_329 = buffer.data(lh + 329);
    const auto *lh_330 = buffer.data(lh + 330);
    const auto *lh_331 = buffer.data(lh + 331);
    const auto *lh_332 = buffer.data(lh + 332);
    const auto *lh_333 = buffer.data(lh + 333);
    const auto *lh_334 = buffer.data(lh + 334);
    const auto *lh_335 = buffer.data(lh + 335);
    const auto *lh_336 = buffer.data(lh + 336);
    const auto *lh_337 = buffer.data(lh + 337);
    const auto *lh_338 = buffer.data(lh + 338);
    const auto *lh_339 = buffer.data(lh + 339);
    const auto *lh_340 = buffer.data(lh + 340);
    const auto *lh_341 = buffer.data(lh + 341);
    const auto *lh_342 = buffer.data(lh + 342);
    const auto *lh_343 = buffer.data(lh + 343);
    const auto *lh_344 = buffer.data(lh + 344);
    const auto *lh_345 = buffer.data(lh + 345);
    const auto *lh_346 = buffer.data(lh + 346);
    const auto *lh_347 = buffer.data(lh + 347);
    const auto *lh_348 = buffer.data(lh + 348);
    const auto *lh_349 = buffer.data(lh + 349);
    const auto *lh_350 = buffer.data(lh + 350);
    const auto *lh_351 = buffer.data(lh + 351);
    const auto *lh_352 = buffer.data(lh + 352);
    const auto *lh_353 = buffer.data(lh + 353);
    const auto *lh_354 = buffer.data(lh + 354);
    const auto *lh_355 = buffer.data(lh + 355);
    const auto *lh_356 = buffer.data(lh + 356);
    const auto *lh_357 = buffer.data(lh + 357);
    const auto *lh_358 = buffer.data(lh + 358);
    const auto *lh_359 = buffer.data(lh + 359);
    const auto *lh_360 = buffer.data(lh + 360);
    const auto *lh_361 = buffer.data(lh + 361);
    const auto *lh_362 = buffer.data(lh + 362);
    const auto *lh_363 = buffer.data(lh + 363);
    const auto *lh_364 = buffer.data(lh + 364);
    const auto *lh_365 = buffer.data(lh + 365);
    const auto *lh_366 = buffer.data(lh + 366);
    const auto *lh_367 = buffer.data(lh + 367);
    const auto *lh_368 = buffer.data(lh + 368);
    const auto *lh_369 = buffer.data(lh + 369);
    const auto *lh_370 = buffer.data(lh + 370);
    const auto *lh_371 = buffer.data(lh + 371);
    const auto *lh_372 = buffer.data(lh + 372);
    const auto *lh_373 = buffer.data(lh + 373);
    const auto *lh_374 = buffer.data(lh + 374);
    const auto *lh_375 = buffer.data(lh + 375);
    const auto *lh_376 = buffer.data(lh + 376);
    const auto *lh_377 = buffer.data(lh + 377);
    const auto *lh_378 = buffer.data(lh + 378);
    const auto *lh_379 = buffer.data(lh + 379);
    const auto *lh_380 = buffer.data(lh + 380);
    const auto *lh_381 = buffer.data(lh + 381);
    const auto *lh_382 = buffer.data(lh + 382);
    const auto *lh_383 = buffer.data(lh + 383);
    const auto *lh_384 = buffer.data(lh + 384);
    const auto *lh_385 = buffer.data(lh + 385);
    const auto *lh_386 = buffer.data(lh + 386);
    const auto *lh_387 = buffer.data(lh + 387);
    const auto *lh_388 = buffer.data(lh + 388);
    const auto *lh_389 = buffer.data(lh + 389);
    const auto *lh_390 = buffer.data(lh + 390);
    const auto *lh_391 = buffer.data(lh + 391);
    const auto *lh_392 = buffer.data(lh + 392);
    const auto *lh_393 = buffer.data(lh + 393);
    const auto *lh_394 = buffer.data(lh + 394);
    const auto *lh_395 = buffer.data(lh + 395);
    const auto *lh_396 = buffer.data(lh + 396);
    const auto *lh_397 = buffer.data(lh + 397);
    const auto *lh_398 = buffer.data(lh + 398);
    const auto *lh_399 = buffer.data(lh + 399);
    const auto *lh_400 = buffer.data(lh + 400);
    const auto *lh_401 = buffer.data(lh + 401);
    const auto *lh_402 = buffer.data(lh + 402);
    const auto *lh_403 = buffer.data(lh + 403);
    const auto *lh_404 = buffer.data(lh + 404);
    const auto *lh_405 = buffer.data(lh + 405);
    const auto *lh_406 = buffer.data(lh + 406);
    const auto *lh_407 = buffer.data(lh + 407);
    const auto *lh_408 = buffer.data(lh + 408);
    const auto *lh_409 = buffer.data(lh + 409);
    const auto *lh_410 = buffer.data(lh + 410);
    const auto *lh_411 = buffer.data(lh + 411);
    const auto *lh_412 = buffer.data(lh + 412);
    const auto *lh_413 = buffer.data(lh + 413);
    const auto *lh_414 = buffer.data(lh + 414);
    const auto *lh_415 = buffer.data(lh + 415);
    const auto *lh_416 = buffer.data(lh + 416);
    const auto *lh_417 = buffer.data(lh + 417);
    const auto *lh_418 = buffer.data(lh + 418);
    const auto *lh_419 = buffer.data(lh + 419);
    const auto *lh_441 = buffer.data(lh + 441);
    const auto *lh_442 = buffer.data(lh + 442);
    const auto *lh_443 = buffer.data(lh + 443);
    const auto *lh_444 = buffer.data(lh + 444);
    const auto *lh_445 = buffer.data(lh + 445);
    const auto *lh_446 = buffer.data(lh + 446);
    const auto *lh_447 = buffer.data(lh + 447);
    const auto *lh_448 = buffer.data(lh + 448);
    const auto *lh_449 = buffer.data(lh + 449);
    const auto *lh_450 = buffer.data(lh + 450);
    const auto *lh_451 = buffer.data(lh + 451);
    const auto *lh_452 = buffer.data(lh + 452);
    const auto *lh_453 = buffer.data(lh + 453);
    const auto *lh_454 = buffer.data(lh + 454);
    const auto *lh_455 = buffer.data(lh + 455);

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, ih_106, ih_107, ih_108, ih_109, \
                         ih_110, lh_253, lh_254, lh_255, lh_256, \
                         lh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -ih_106[k]
                   + f_0 * lh_253[k];

        t_170[k] = -ih_107[k]
                   + f_0 * lh_254[k];

        t_171[k] = -ih_108[k]
                   + f_0 * lh_255[k];

        t_172[k] = -ih_109[k]
                   + f_0 * lh_256[k];

        t_173[k] = -ih_110[k]
                   + f_0 * lh_257[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, ih_111, ih_112, ih_113, ih_114, \
                         ih_115, lh_258, lh_259, lh_260, lh_261, \
                         lh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -ih_111[k]
                   + f_0 * lh_258[k];

        t_175[k] = -ih_112[k]
                   + f_0 * lh_259[k];

        t_176[k] = -ih_113[k]
                   + f_0 * lh_260[k];

        t_177[k] = -ih_114[k]
                   + f_0 * lh_261[k];

        t_178[k] = -ih_115[k]
                   + f_0 * lh_262[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, ih_116, ih_117, ih_118, ih_119, \
                         ih_120, lh_263, lh_264, lh_265, lh_266, \
                         lh_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -ih_116[k]
                   + f_0 * lh_263[k];

        t_180[k] = -ih_117[k]
                   + f_0 * lh_264[k];

        t_181[k] = -ih_118[k]
                   + f_0 * lh_265[k];

        t_182[k] = -ih_119[k]
                   + f_0 * lh_266[k];

        t_183[k] = -ih_120[k]
                   + f_0 * lh_267[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, ih_121, ih_122, ih_123, ih_124, \
                         ih_125, lh_268, lh_269, lh_270, lh_271, \
                         lh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -ih_121[k]
                   + f_0 * lh_268[k];

        t_185[k] = -ih_122[k]
                   + f_0 * lh_269[k];

        t_186[k] = -ih_123[k]
                   + f_0 * lh_270[k];

        t_187[k] = -ih_124[k]
                   + f_0 * lh_271[k];

        t_188[k] = -ih_125[k]
                   + f_0 * lh_272[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, t_196, lh_273, \
                         lh_274, lh_275, lh_276, lh_277, lh_278, lh_279, \
                         lh_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * lh_273[k];

        t_190[k] = f_0 * lh_274[k];

        t_191[k] = f_0 * lh_275[k];

        t_192[k] = f_0 * lh_276[k];

        t_193[k] = f_0 * lh_277[k];

        t_194[k] = f_0 * lh_278[k];

        t_195[k] = f_0 * lh_279[k];

        t_196[k] = f_0 * lh_280[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, t_202, t_203, t_204, lh_281, \
                         lh_282, lh_283, lh_284, lh_285, lh_286, lh_287, \
                         lh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_0 * lh_281[k];

        t_198[k] = f_0 * lh_282[k];

        t_199[k] = f_0 * lh_283[k];

        t_200[k] = f_0 * lh_284[k];

        t_201[k] = f_0 * lh_285[k];

        t_202[k] = f_0 * lh_286[k];

        t_203[k] = f_0 * lh_287[k];

        t_204[k] = f_0 * lh_288[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, t_210, t_211, ih_126, ih_127, \
                         lh_289, lh_290, lh_291, lh_292, lh_293, lh_315, \
                         lh_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_0 * lh_289[k];

        t_206[k] = f_0 * lh_290[k];

        t_207[k] = f_0 * lh_291[k];

        t_208[k] = f_0 * lh_292[k];

        t_209[k] = f_0 * lh_293[k];

        t_210[k] = -4.0 * ih_126[k]
                   + f_0 * lh_315[k];

        t_211[k] = -4.0 * ih_127[k]
                   + f_0 * lh_316[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, ih_128, ih_129, ih_130, ih_131, \
                         ih_132, lh_317, lh_318, lh_319, lh_320, \
                         lh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = -4.0 * ih_128[k]
                   + f_0 * lh_317[k];

        t_213[k] = -4.0 * ih_129[k]
                   + f_0 * lh_318[k];

        t_214[k] = -4.0 * ih_130[k]
                   + f_0 * lh_319[k];

        t_215[k] = -4.0 * ih_131[k]
                   + f_0 * lh_320[k];

        t_216[k] = -4.0 * ih_132[k]
                   + f_0 * lh_321[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, ih_133, ih_134, ih_135, ih_136, \
                         ih_137, lh_322, lh_323, lh_324, lh_325, \
                         lh_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = -4.0 * ih_133[k]
                   + f_0 * lh_322[k];

        t_218[k] = -4.0 * ih_134[k]
                   + f_0 * lh_323[k];

        t_219[k] = -4.0 * ih_135[k]
                   + f_0 * lh_324[k];

        t_220[k] = -4.0 * ih_136[k]
                   + f_0 * lh_325[k];

        t_221[k] = -4.0 * ih_137[k]
                   + f_0 * lh_326[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, ih_138, ih_139, ih_140, ih_141, \
                         ih_142, lh_327, lh_328, lh_329, lh_330, \
                         lh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -4.0 * ih_138[k]
                   + f_0 * lh_327[k];

        t_223[k] = -4.0 * ih_139[k]
                   + f_0 * lh_328[k];

        t_224[k] = -4.0 * ih_140[k]
                   + f_0 * lh_329[k];

        t_225[k] = -4.0 * ih_141[k]
                   + f_0 * lh_330[k];

        t_226[k] = -4.0 * ih_142[k]
                   + f_0 * lh_331[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, t_231, ih_143, ih_144, ih_145, ih_146, \
                         ih_147, lh_332, lh_333, lh_334, lh_335, \
                         lh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = -4.0 * ih_143[k]
                   + f_0 * lh_332[k];

        t_228[k] = -4.0 * ih_144[k]
                   + f_0 * lh_333[k];

        t_229[k] = -4.0 * ih_145[k]
                   + f_0 * lh_334[k];

        t_230[k] = -4.0 * ih_146[k]
                   + f_0 * lh_335[k];

        t_231[k] = -3.0 * ih_147[k]
                   + f_0 * lh_336[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, ih_148, ih_149, ih_150, ih_151, \
                         ih_152, lh_337, lh_338, lh_339, lh_340, \
                         lh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = -3.0 * ih_148[k]
                   + f_0 * lh_337[k];

        t_233[k] = -3.0 * ih_149[k]
                   + f_0 * lh_338[k];

        t_234[k] = -3.0 * ih_150[k]
                   + f_0 * lh_339[k];

        t_235[k] = -3.0 * ih_151[k]
                   + f_0 * lh_340[k];

        t_236[k] = -3.0 * ih_152[k]
                   + f_0 * lh_341[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, ih_153, ih_154, ih_155, ih_156, \
                         ih_157, lh_342, lh_343, lh_344, lh_345, \
                         lh_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = -3.0 * ih_153[k]
                   + f_0 * lh_342[k];

        t_238[k] = -3.0 * ih_154[k]
                   + f_0 * lh_343[k];

        t_239[k] = -3.0 * ih_155[k]
                   + f_0 * lh_344[k];

        t_240[k] = -3.0 * ih_156[k]
                   + f_0 * lh_345[k];

        t_241[k] = -3.0 * ih_157[k]
                   + f_0 * lh_346[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, ih_158, ih_159, ih_160, ih_161, \
                         ih_162, lh_347, lh_348, lh_349, lh_350, \
                         lh_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = -3.0 * ih_158[k]
                   + f_0 * lh_347[k];

        t_243[k] = -3.0 * ih_159[k]
                   + f_0 * lh_348[k];

        t_244[k] = -3.0 * ih_160[k]
                   + f_0 * lh_349[k];

        t_245[k] = -3.0 * ih_161[k]
                   + f_0 * lh_350[k];

        t_246[k] = -3.0 * ih_162[k]
                   + f_0 * lh_351[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, ih_163, ih_164, ih_165, ih_166, \
                         ih_167, lh_352, lh_353, lh_354, lh_355, \
                         lh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = -3.0 * ih_163[k]
                   + f_0 * lh_352[k];

        t_248[k] = -3.0 * ih_164[k]
                   + f_0 * lh_353[k];

        t_249[k] = -3.0 * ih_165[k]
                   + f_0 * lh_354[k];

        t_250[k] = -3.0 * ih_166[k]
                   + f_0 * lh_355[k];

        t_251[k] = -3.0 * ih_167[k]
                   + f_0 * lh_356[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, ih_168, ih_169, ih_170, ih_171, \
                         ih_172, lh_357, lh_358, lh_359, lh_360, \
                         lh_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -2.0 * ih_168[k]
                   + f_0 * lh_357[k];

        t_253[k] = -2.0 * ih_169[k]
                   + f_0 * lh_358[k];

        t_254[k] = -2.0 * ih_170[k]
                   + f_0 * lh_359[k];

        t_255[k] = -2.0 * ih_171[k]
                   + f_0 * lh_360[k];

        t_256[k] = -2.0 * ih_172[k]
                   + f_0 * lh_361[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, ih_173, ih_174, ih_175, ih_176, \
                         ih_177, lh_362, lh_363, lh_364, lh_365, \
                         lh_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = -2.0 * ih_173[k]
                   + f_0 * lh_362[k];

        t_258[k] = -2.0 * ih_174[k]
                   + f_0 * lh_363[k];

        t_259[k] = -2.0 * ih_175[k]
                   + f_0 * lh_364[k];

        t_260[k] = -2.0 * ih_176[k]
                   + f_0 * lh_365[k];

        t_261[k] = -2.0 * ih_177[k]
                   + f_0 * lh_366[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, ih_178, ih_179, ih_180, ih_181, \
                         ih_182, lh_367, lh_368, lh_369, lh_370, \
                         lh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = -2.0 * ih_178[k]
                   + f_0 * lh_367[k];

        t_263[k] = -2.0 * ih_179[k]
                   + f_0 * lh_368[k];

        t_264[k] = -2.0 * ih_180[k]
                   + f_0 * lh_369[k];

        t_265[k] = -2.0 * ih_181[k]
                   + f_0 * lh_370[k];

        t_266[k] = -2.0 * ih_182[k]
                   + f_0 * lh_371[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, ih_183, ih_184, ih_185, ih_186, \
                         ih_187, lh_372, lh_373, lh_374, lh_375, \
                         lh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = -2.0 * ih_183[k]
                   + f_0 * lh_372[k];

        t_268[k] = -2.0 * ih_184[k]
                   + f_0 * lh_373[k];

        t_269[k] = -2.0 * ih_185[k]
                   + f_0 * lh_374[k];

        t_270[k] = -2.0 * ih_186[k]
                   + f_0 * lh_375[k];

        t_271[k] = -2.0 * ih_187[k]
                   + f_0 * lh_376[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, ih_188, ih_189, ih_190, ih_191, \
                         ih_192, lh_377, lh_378, lh_379, lh_380, \
                         lh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -2.0 * ih_188[k]
                   + f_0 * lh_377[k];

        t_273[k] = -ih_189[k]
                   + f_0 * lh_378[k];

        t_274[k] = -ih_190[k]
                   + f_0 * lh_379[k];

        t_275[k] = -ih_191[k]
                   + f_0 * lh_380[k];

        t_276[k] = -ih_192[k]
                   + f_0 * lh_381[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, t_281, ih_193, ih_194, ih_195, ih_196, \
                         ih_197, lh_382, lh_383, lh_384, lh_385, \
                         lh_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = -ih_193[k]
                   + f_0 * lh_382[k];

        t_278[k] = -ih_194[k]
                   + f_0 * lh_383[k];

        t_279[k] = -ih_195[k]
                   + f_0 * lh_384[k];

        t_280[k] = -ih_196[k]
                   + f_0 * lh_385[k];

        t_281[k] = -ih_197[k]
                   + f_0 * lh_386[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, ih_198, ih_199, ih_200, ih_201, \
                         ih_202, lh_387, lh_388, lh_389, lh_390, \
                         lh_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = -ih_198[k]
                   + f_0 * lh_387[k];

        t_283[k] = -ih_199[k]
                   + f_0 * lh_388[k];

        t_284[k] = -ih_200[k]
                   + f_0 * lh_389[k];

        t_285[k] = -ih_201[k]
                   + f_0 * lh_390[k];

        t_286[k] = -ih_202[k]
                   + f_0 * lh_391[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, ih_203, ih_204, ih_205, ih_206, \
                         ih_207, lh_392, lh_393, lh_394, lh_395, \
                         lh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = -ih_203[k]
                   + f_0 * lh_392[k];

        t_288[k] = -ih_204[k]
                   + f_0 * lh_393[k];

        t_289[k] = -ih_205[k]
                   + f_0 * lh_394[k];

        t_290[k] = -ih_206[k]
                   + f_0 * lh_395[k];

        t_291[k] = -ih_207[k]
                   + f_0 * lh_396[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, t_297, t_298, ih_208, ih_209, \
                         lh_397, lh_398, lh_399, lh_400, lh_401, lh_402, \
                         lh_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = -ih_208[k]
                   + f_0 * lh_397[k];

        t_293[k] = -ih_209[k]
                   + f_0 * lh_398[k];

        t_294[k] = f_0 * lh_399[k];

        t_295[k] = f_0 * lh_400[k];

        t_296[k] = f_0 * lh_401[k];

        t_297[k] = f_0 * lh_402[k];

        t_298[k] = f_0 * lh_403[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, t_304, t_305, t_306, lh_404, \
                         lh_405, lh_406, lh_407, lh_408, lh_409, lh_410, \
                         lh_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_0 * lh_404[k];

        t_300[k] = f_0 * lh_405[k];

        t_301[k] = f_0 * lh_406[k];

        t_302[k] = f_0 * lh_407[k];

        t_303[k] = f_0 * lh_408[k];

        t_304[k] = f_0 * lh_409[k];

        t_305[k] = f_0 * lh_410[k];

        t_306[k] = f_0 * lh_411[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, t_312, t_313, t_314, lh_412, \
                         lh_413, lh_414, lh_415, lh_416, lh_417, lh_418, \
                         lh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_0 * lh_412[k];

        t_308[k] = f_0 * lh_413[k];

        t_309[k] = f_0 * lh_414[k];

        t_310[k] = f_0 * lh_415[k];

        t_311[k] = f_0 * lh_416[k];

        t_312[k] = f_0 * lh_417[k];

        t_313[k] = f_0 * lh_418[k];

        t_314[k] = f_0 * lh_419[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ih_210, ih_211, ih_212, ih_213, \
                         ih_214, lh_441, lh_442, lh_443, lh_444, \
                         lh_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -5.0 * ih_210[k]
                   + f_0 * lh_441[k];

        t_316[k] = -5.0 * ih_211[k]
                   + f_0 * lh_442[k];

        t_317[k] = -5.0 * ih_212[k]
                   + f_0 * lh_443[k];

        t_318[k] = -5.0 * ih_213[k]
                   + f_0 * lh_444[k];

        t_319[k] = -5.0 * ih_214[k]
                   + f_0 * lh_445[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ih_215, ih_216, ih_217, ih_218, \
                         ih_219, lh_446, lh_447, lh_448, lh_449, \
                         lh_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -5.0 * ih_215[k]
                   + f_0 * lh_446[k];

        t_321[k] = -5.0 * ih_216[k]
                   + f_0 * lh_447[k];

        t_322[k] = -5.0 * ih_217[k]
                   + f_0 * lh_448[k];

        t_323[k] = -5.0 * ih_218[k]
                   + f_0 * lh_449[k];

        t_324[k] = -5.0 * ih_219[k]
                   + f_0 * lh_450[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ih_220, ih_221, ih_222, ih_223, \
                         ih_224, lh_451, lh_452, lh_453, lh_454, \
                         lh_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -5.0 * ih_220[k]
                   + f_0 * lh_451[k];

        t_326[k] = -5.0 * ih_221[k]
                   + f_0 * lh_452[k];

        t_327[k] = -5.0 * ih_222[k]
                   + f_0 * lh_453[k];

        t_328[k] = -5.0 * ih_223[k]
                   + f_0 * lh_454[k];

        t_329[k] = -5.0 * ih_224[k]
                   + f_0 * lh_455[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
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
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);

    const auto *ih_225 = buffer.data(ih + 225);
    const auto *ih_226 = buffer.data(ih + 226);
    const auto *ih_227 = buffer.data(ih + 227);
    const auto *ih_228 = buffer.data(ih + 228);
    const auto *ih_229 = buffer.data(ih + 229);
    const auto *ih_230 = buffer.data(ih + 230);
    const auto *ih_231 = buffer.data(ih + 231);
    const auto *ih_232 = buffer.data(ih + 232);
    const auto *ih_233 = buffer.data(ih + 233);
    const auto *ih_234 = buffer.data(ih + 234);
    const auto *ih_235 = buffer.data(ih + 235);
    const auto *ih_236 = buffer.data(ih + 236);
    const auto *ih_237 = buffer.data(ih + 237);
    const auto *ih_238 = buffer.data(ih + 238);
    const auto *ih_239 = buffer.data(ih + 239);
    const auto *ih_240 = buffer.data(ih + 240);
    const auto *ih_241 = buffer.data(ih + 241);
    const auto *ih_242 = buffer.data(ih + 242);
    const auto *ih_243 = buffer.data(ih + 243);
    const auto *ih_244 = buffer.data(ih + 244);
    const auto *ih_245 = buffer.data(ih + 245);
    const auto *ih_246 = buffer.data(ih + 246);
    const auto *ih_247 = buffer.data(ih + 247);
    const auto *ih_248 = buffer.data(ih + 248);
    const auto *ih_249 = buffer.data(ih + 249);
    const auto *ih_250 = buffer.data(ih + 250);
    const auto *ih_251 = buffer.data(ih + 251);
    const auto *ih_252 = buffer.data(ih + 252);
    const auto *ih_253 = buffer.data(ih + 253);
    const auto *ih_254 = buffer.data(ih + 254);
    const auto *ih_255 = buffer.data(ih + 255);
    const auto *ih_256 = buffer.data(ih + 256);
    const auto *ih_257 = buffer.data(ih + 257);
    const auto *ih_258 = buffer.data(ih + 258);
    const auto *ih_259 = buffer.data(ih + 259);
    const auto *ih_260 = buffer.data(ih + 260);
    const auto *ih_261 = buffer.data(ih + 261);
    const auto *ih_262 = buffer.data(ih + 262);
    const auto *ih_263 = buffer.data(ih + 263);
    const auto *ih_264 = buffer.data(ih + 264);
    const auto *ih_265 = buffer.data(ih + 265);
    const auto *ih_266 = buffer.data(ih + 266);
    const auto *ih_267 = buffer.data(ih + 267);
    const auto *ih_268 = buffer.data(ih + 268);
    const auto *ih_269 = buffer.data(ih + 269);
    const auto *ih_270 = buffer.data(ih + 270);
    const auto *ih_271 = buffer.data(ih + 271);
    const auto *ih_272 = buffer.data(ih + 272);
    const auto *ih_273 = buffer.data(ih + 273);
    const auto *ih_274 = buffer.data(ih + 274);
    const auto *ih_275 = buffer.data(ih + 275);
    const auto *ih_276 = buffer.data(ih + 276);
    const auto *ih_277 = buffer.data(ih + 277);
    const auto *ih_278 = buffer.data(ih + 278);
    const auto *ih_279 = buffer.data(ih + 279);
    const auto *ih_280 = buffer.data(ih + 280);
    const auto *ih_281 = buffer.data(ih + 281);
    const auto *ih_282 = buffer.data(ih + 282);
    const auto *ih_283 = buffer.data(ih + 283);
    const auto *ih_284 = buffer.data(ih + 284);
    const auto *ih_285 = buffer.data(ih + 285);
    const auto *ih_286 = buffer.data(ih + 286);
    const auto *ih_287 = buffer.data(ih + 287);
    const auto *ih_288 = buffer.data(ih + 288);
    const auto *ih_289 = buffer.data(ih + 289);
    const auto *ih_290 = buffer.data(ih + 290);
    const auto *ih_291 = buffer.data(ih + 291);
    const auto *ih_292 = buffer.data(ih + 292);
    const auto *ih_293 = buffer.data(ih + 293);
    const auto *ih_294 = buffer.data(ih + 294);
    const auto *ih_295 = buffer.data(ih + 295);
    const auto *ih_296 = buffer.data(ih + 296);
    const auto *ih_297 = buffer.data(ih + 297);
    const auto *ih_298 = buffer.data(ih + 298);
    const auto *ih_299 = buffer.data(ih + 299);
    const auto *ih_300 = buffer.data(ih + 300);
    const auto *ih_301 = buffer.data(ih + 301);
    const auto *ih_302 = buffer.data(ih + 302);
    const auto *ih_303 = buffer.data(ih + 303);
    const auto *ih_304 = buffer.data(ih + 304);
    const auto *ih_305 = buffer.data(ih + 305);
    const auto *ih_306 = buffer.data(ih + 306);
    const auto *ih_307 = buffer.data(ih + 307);
    const auto *ih_308 = buffer.data(ih + 308);
    const auto *ih_309 = buffer.data(ih + 309);
    const auto *ih_310 = buffer.data(ih + 310);
    const auto *ih_311 = buffer.data(ih + 311);
    const auto *ih_312 = buffer.data(ih + 312);
    const auto *ih_313 = buffer.data(ih + 313);
    const auto *ih_314 = buffer.data(ih + 314);
    const auto *ih_315 = buffer.data(ih + 315);
    const auto *ih_316 = buffer.data(ih + 316);
    const auto *ih_317 = buffer.data(ih + 317);
    const auto *ih_318 = buffer.data(ih + 318);
    const auto *ih_319 = buffer.data(ih + 319);
    const auto *ih_320 = buffer.data(ih + 320);
    const auto *ih_321 = buffer.data(ih + 321);
    const auto *ih_322 = buffer.data(ih + 322);
    const auto *ih_323 = buffer.data(ih + 323);
    const auto *ih_324 = buffer.data(ih + 324);
    const auto *ih_325 = buffer.data(ih + 325);
    const auto *ih_326 = buffer.data(ih + 326);
    const auto *ih_327 = buffer.data(ih + 327);
    const auto *ih_328 = buffer.data(ih + 328);
    const auto *ih_329 = buffer.data(ih + 329);
    const auto *ih_330 = buffer.data(ih + 330);
    const auto *ih_331 = buffer.data(ih + 331);
    const auto *ih_332 = buffer.data(ih + 332);
    const auto *ih_333 = buffer.data(ih + 333);
    const auto *ih_334 = buffer.data(ih + 334);
    const auto *ih_335 = buffer.data(ih + 335);
    const auto *ih_336 = buffer.data(ih + 336);
    const auto *ih_337 = buffer.data(ih + 337);
    const auto *ih_338 = buffer.data(ih + 338);
    const auto *ih_339 = buffer.data(ih + 339);
    const auto *ih_340 = buffer.data(ih + 340);
    const auto *ih_341 = buffer.data(ih + 341);
    const auto *ih_342 = buffer.data(ih + 342);
    const auto *ih_343 = buffer.data(ih + 343);
    const auto *ih_344 = buffer.data(ih + 344);
    const auto *ih_345 = buffer.data(ih + 345);
    const auto *ih_346 = buffer.data(ih + 346);
    const auto *ih_347 = buffer.data(ih + 347);
    const auto *ih_348 = buffer.data(ih + 348);
    const auto *ih_349 = buffer.data(ih + 349);
    const auto *ih_350 = buffer.data(ih + 350);
    const auto *ih_351 = buffer.data(ih + 351);
    const auto *ih_352 = buffer.data(ih + 352);
    const auto *ih_353 = buffer.data(ih + 353);
    const auto *ih_354 = buffer.data(ih + 354);
    const auto *ih_355 = buffer.data(ih + 355);
    const auto *ih_356 = buffer.data(ih + 356);
    const auto *ih_357 = buffer.data(ih + 357);
    const auto *ih_358 = buffer.data(ih + 358);
    const auto *ih_359 = buffer.data(ih + 359);
    const auto *ih_360 = buffer.data(ih + 360);
    const auto *ih_361 = buffer.data(ih + 361);

    const auto *lh_456 = buffer.data(lh + 456);
    const auto *lh_457 = buffer.data(lh + 457);
    const auto *lh_458 = buffer.data(lh + 458);
    const auto *lh_459 = buffer.data(lh + 459);
    const auto *lh_460 = buffer.data(lh + 460);
    const auto *lh_461 = buffer.data(lh + 461);
    const auto *lh_462 = buffer.data(lh + 462);
    const auto *lh_463 = buffer.data(lh + 463);
    const auto *lh_464 = buffer.data(lh + 464);
    const auto *lh_465 = buffer.data(lh + 465);
    const auto *lh_466 = buffer.data(lh + 466);
    const auto *lh_467 = buffer.data(lh + 467);
    const auto *lh_468 = buffer.data(lh + 468);
    const auto *lh_469 = buffer.data(lh + 469);
    const auto *lh_470 = buffer.data(lh + 470);
    const auto *lh_471 = buffer.data(lh + 471);
    const auto *lh_472 = buffer.data(lh + 472);
    const auto *lh_473 = buffer.data(lh + 473);
    const auto *lh_474 = buffer.data(lh + 474);
    const auto *lh_475 = buffer.data(lh + 475);
    const auto *lh_476 = buffer.data(lh + 476);
    const auto *lh_477 = buffer.data(lh + 477);
    const auto *lh_478 = buffer.data(lh + 478);
    const auto *lh_479 = buffer.data(lh + 479);
    const auto *lh_480 = buffer.data(lh + 480);
    const auto *lh_481 = buffer.data(lh + 481);
    const auto *lh_482 = buffer.data(lh + 482);
    const auto *lh_483 = buffer.data(lh + 483);
    const auto *lh_484 = buffer.data(lh + 484);
    const auto *lh_485 = buffer.data(lh + 485);
    const auto *lh_486 = buffer.data(lh + 486);
    const auto *lh_487 = buffer.data(lh + 487);
    const auto *lh_488 = buffer.data(lh + 488);
    const auto *lh_489 = buffer.data(lh + 489);
    const auto *lh_490 = buffer.data(lh + 490);
    const auto *lh_491 = buffer.data(lh + 491);
    const auto *lh_492 = buffer.data(lh + 492);
    const auto *lh_493 = buffer.data(lh + 493);
    const auto *lh_494 = buffer.data(lh + 494);
    const auto *lh_495 = buffer.data(lh + 495);
    const auto *lh_496 = buffer.data(lh + 496);
    const auto *lh_497 = buffer.data(lh + 497);
    const auto *lh_498 = buffer.data(lh + 498);
    const auto *lh_499 = buffer.data(lh + 499);
    const auto *lh_500 = buffer.data(lh + 500);
    const auto *lh_501 = buffer.data(lh + 501);
    const auto *lh_502 = buffer.data(lh + 502);
    const auto *lh_503 = buffer.data(lh + 503);
    const auto *lh_504 = buffer.data(lh + 504);
    const auto *lh_505 = buffer.data(lh + 505);
    const auto *lh_506 = buffer.data(lh + 506);
    const auto *lh_507 = buffer.data(lh + 507);
    const auto *lh_508 = buffer.data(lh + 508);
    const auto *lh_509 = buffer.data(lh + 509);
    const auto *lh_510 = buffer.data(lh + 510);
    const auto *lh_511 = buffer.data(lh + 511);
    const auto *lh_512 = buffer.data(lh + 512);
    const auto *lh_513 = buffer.data(lh + 513);
    const auto *lh_514 = buffer.data(lh + 514);
    const auto *lh_515 = buffer.data(lh + 515);
    const auto *lh_516 = buffer.data(lh + 516);
    const auto *lh_517 = buffer.data(lh + 517);
    const auto *lh_518 = buffer.data(lh + 518);
    const auto *lh_519 = buffer.data(lh + 519);
    const auto *lh_520 = buffer.data(lh + 520);
    const auto *lh_521 = buffer.data(lh + 521);
    const auto *lh_522 = buffer.data(lh + 522);
    const auto *lh_523 = buffer.data(lh + 523);
    const auto *lh_524 = buffer.data(lh + 524);
    const auto *lh_525 = buffer.data(lh + 525);
    const auto *lh_526 = buffer.data(lh + 526);
    const auto *lh_527 = buffer.data(lh + 527);
    const auto *lh_528 = buffer.data(lh + 528);
    const auto *lh_529 = buffer.data(lh + 529);
    const auto *lh_530 = buffer.data(lh + 530);
    const auto *lh_531 = buffer.data(lh + 531);
    const auto *lh_532 = buffer.data(lh + 532);
    const auto *lh_533 = buffer.data(lh + 533);
    const auto *lh_534 = buffer.data(lh + 534);
    const auto *lh_535 = buffer.data(lh + 535);
    const auto *lh_536 = buffer.data(lh + 536);
    const auto *lh_537 = buffer.data(lh + 537);
    const auto *lh_538 = buffer.data(lh + 538);
    const auto *lh_539 = buffer.data(lh + 539);
    const auto *lh_540 = buffer.data(lh + 540);
    const auto *lh_541 = buffer.data(lh + 541);
    const auto *lh_542 = buffer.data(lh + 542);
    const auto *lh_543 = buffer.data(lh + 543);
    const auto *lh_544 = buffer.data(lh + 544);
    const auto *lh_545 = buffer.data(lh + 545);
    const auto *lh_546 = buffer.data(lh + 546);
    const auto *lh_547 = buffer.data(lh + 547);
    const auto *lh_548 = buffer.data(lh + 548);
    const auto *lh_549 = buffer.data(lh + 549);
    const auto *lh_550 = buffer.data(lh + 550);
    const auto *lh_551 = buffer.data(lh + 551);
    const auto *lh_552 = buffer.data(lh + 552);
    const auto *lh_553 = buffer.data(lh + 553);
    const auto *lh_554 = buffer.data(lh + 554);
    const auto *lh_555 = buffer.data(lh + 555);
    const auto *lh_556 = buffer.data(lh + 556);
    const auto *lh_557 = buffer.data(lh + 557);
    const auto *lh_558 = buffer.data(lh + 558);
    const auto *lh_559 = buffer.data(lh + 559);
    const auto *lh_560 = buffer.data(lh + 560);
    const auto *lh_561 = buffer.data(lh + 561);
    const auto *lh_562 = buffer.data(lh + 562);
    const auto *lh_563 = buffer.data(lh + 563);
    const auto *lh_564 = buffer.data(lh + 564);
    const auto *lh_565 = buffer.data(lh + 565);
    const auto *lh_566 = buffer.data(lh + 566);
    const auto *lh_588 = buffer.data(lh + 588);
    const auto *lh_589 = buffer.data(lh + 589);
    const auto *lh_590 = buffer.data(lh + 590);
    const auto *lh_591 = buffer.data(lh + 591);
    const auto *lh_592 = buffer.data(lh + 592);
    const auto *lh_593 = buffer.data(lh + 593);
    const auto *lh_594 = buffer.data(lh + 594);
    const auto *lh_595 = buffer.data(lh + 595);
    const auto *lh_596 = buffer.data(lh + 596);
    const auto *lh_597 = buffer.data(lh + 597);
    const auto *lh_598 = buffer.data(lh + 598);
    const auto *lh_599 = buffer.data(lh + 599);
    const auto *lh_600 = buffer.data(lh + 600);
    const auto *lh_601 = buffer.data(lh + 601);
    const auto *lh_602 = buffer.data(lh + 602);
    const auto *lh_603 = buffer.data(lh + 603);
    const auto *lh_604 = buffer.data(lh + 604);
    const auto *lh_605 = buffer.data(lh + 605);
    const auto *lh_606 = buffer.data(lh + 606);
    const auto *lh_607 = buffer.data(lh + 607);
    const auto *lh_608 = buffer.data(lh + 608);
    const auto *lh_609 = buffer.data(lh + 609);
    const auto *lh_610 = buffer.data(lh + 610);
    const auto *lh_611 = buffer.data(lh + 611);
    const auto *lh_612 = buffer.data(lh + 612);
    const auto *lh_613 = buffer.data(lh + 613);
    const auto *lh_614 = buffer.data(lh + 614);
    const auto *lh_615 = buffer.data(lh + 615);
    const auto *lh_616 = buffer.data(lh + 616);
    const auto *lh_617 = buffer.data(lh + 617);
    const auto *lh_618 = buffer.data(lh + 618);
    const auto *lh_619 = buffer.data(lh + 619);
    const auto *lh_620 = buffer.data(lh + 620);
    const auto *lh_621 = buffer.data(lh + 621);
    const auto *lh_622 = buffer.data(lh + 622);
    const auto *lh_623 = buffer.data(lh + 623);
    const auto *lh_624 = buffer.data(lh + 624);
    const auto *lh_625 = buffer.data(lh + 625);
    const auto *lh_626 = buffer.data(lh + 626);
    const auto *lh_627 = buffer.data(lh + 627);
    const auto *lh_628 = buffer.data(lh + 628);
    const auto *lh_629 = buffer.data(lh + 629);
    const auto *lh_630 = buffer.data(lh + 630);
    const auto *lh_631 = buffer.data(lh + 631);
    const auto *lh_632 = buffer.data(lh + 632);
    const auto *lh_633 = buffer.data(lh + 633);
    const auto *lh_634 = buffer.data(lh + 634);

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ih_225, ih_226, ih_227, ih_228, \
                         ih_229, lh_456, lh_457, lh_458, lh_459, \
                         lh_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -5.0 * ih_225[k]
                   + f_0 * lh_456[k];

        t_331[k] = -5.0 * ih_226[k]
                   + f_0 * lh_457[k];

        t_332[k] = -5.0 * ih_227[k]
                   + f_0 * lh_458[k];

        t_333[k] = -5.0 * ih_228[k]
                   + f_0 * lh_459[k];

        t_334[k] = -5.0 * ih_229[k]
                   + f_0 * lh_460[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ih_230, ih_231, ih_232, ih_233, \
                         ih_234, lh_461, lh_462, lh_463, lh_464, \
                         lh_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -5.0 * ih_230[k]
                   + f_0 * lh_461[k];

        t_336[k] = -4.0 * ih_231[k]
                   + f_0 * lh_462[k];

        t_337[k] = -4.0 * ih_232[k]
                   + f_0 * lh_463[k];

        t_338[k] = -4.0 * ih_233[k]
                   + f_0 * lh_464[k];

        t_339[k] = -4.0 * ih_234[k]
                   + f_0 * lh_465[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ih_235, ih_236, ih_237, ih_238, \
                         ih_239, lh_466, lh_467, lh_468, lh_469, \
                         lh_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -4.0 * ih_235[k]
                   + f_0 * lh_466[k];

        t_341[k] = -4.0 * ih_236[k]
                   + f_0 * lh_467[k];

        t_342[k] = -4.0 * ih_237[k]
                   + f_0 * lh_468[k];

        t_343[k] = -4.0 * ih_238[k]
                   + f_0 * lh_469[k];

        t_344[k] = -4.0 * ih_239[k]
                   + f_0 * lh_470[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ih_240, ih_241, ih_242, ih_243, \
                         ih_244, lh_471, lh_472, lh_473, lh_474, \
                         lh_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -4.0 * ih_240[k]
                   + f_0 * lh_471[k];

        t_346[k] = -4.0 * ih_241[k]
                   + f_0 * lh_472[k];

        t_347[k] = -4.0 * ih_242[k]
                   + f_0 * lh_473[k];

        t_348[k] = -4.0 * ih_243[k]
                   + f_0 * lh_474[k];

        t_349[k] = -4.0 * ih_244[k]
                   + f_0 * lh_475[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ih_245, ih_246, ih_247, ih_248, \
                         ih_249, lh_476, lh_477, lh_478, lh_479, \
                         lh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -4.0 * ih_245[k]
                   + f_0 * lh_476[k];

        t_351[k] = -4.0 * ih_246[k]
                   + f_0 * lh_477[k];

        t_352[k] = -4.0 * ih_247[k]
                   + f_0 * lh_478[k];

        t_353[k] = -4.0 * ih_248[k]
                   + f_0 * lh_479[k];

        t_354[k] = -4.0 * ih_249[k]
                   + f_0 * lh_480[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ih_250, ih_251, ih_252, ih_253, \
                         ih_254, lh_481, lh_482, lh_483, lh_484, \
                         lh_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -4.0 * ih_250[k]
                   + f_0 * lh_481[k];

        t_356[k] = -4.0 * ih_251[k]
                   + f_0 * lh_482[k];

        t_357[k] = -3.0 * ih_252[k]
                   + f_0 * lh_483[k];

        t_358[k] = -3.0 * ih_253[k]
                   + f_0 * lh_484[k];

        t_359[k] = -3.0 * ih_254[k]
                   + f_0 * lh_485[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ih_255, ih_256, ih_257, ih_258, \
                         ih_259, lh_486, lh_487, lh_488, lh_489, \
                         lh_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -3.0 * ih_255[k]
                   + f_0 * lh_486[k];

        t_361[k] = -3.0 * ih_256[k]
                   + f_0 * lh_487[k];

        t_362[k] = -3.0 * ih_257[k]
                   + f_0 * lh_488[k];

        t_363[k] = -3.0 * ih_258[k]
                   + f_0 * lh_489[k];

        t_364[k] = -3.0 * ih_259[k]
                   + f_0 * lh_490[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ih_260, ih_261, ih_262, ih_263, \
                         ih_264, lh_491, lh_492, lh_493, lh_494, \
                         lh_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -3.0 * ih_260[k]
                   + f_0 * lh_491[k];

        t_366[k] = -3.0 * ih_261[k]
                   + f_0 * lh_492[k];

        t_367[k] = -3.0 * ih_262[k]
                   + f_0 * lh_493[k];

        t_368[k] = -3.0 * ih_263[k]
                   + f_0 * lh_494[k];

        t_369[k] = -3.0 * ih_264[k]
                   + f_0 * lh_495[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ih_265, ih_266, ih_267, ih_268, \
                         ih_269, lh_496, lh_497, lh_498, lh_499, \
                         lh_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -3.0 * ih_265[k]
                   + f_0 * lh_496[k];

        t_371[k] = -3.0 * ih_266[k]
                   + f_0 * lh_497[k];

        t_372[k] = -3.0 * ih_267[k]
                   + f_0 * lh_498[k];

        t_373[k] = -3.0 * ih_268[k]
                   + f_0 * lh_499[k];

        t_374[k] = -3.0 * ih_269[k]
                   + f_0 * lh_500[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ih_270, ih_271, ih_272, ih_273, \
                         ih_274, lh_501, lh_502, lh_503, lh_504, \
                         lh_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -3.0 * ih_270[k]
                   + f_0 * lh_501[k];

        t_376[k] = -3.0 * ih_271[k]
                   + f_0 * lh_502[k];

        t_377[k] = -3.0 * ih_272[k]
                   + f_0 * lh_503[k];

        t_378[k] = -2.0 * ih_273[k]
                   + f_0 * lh_504[k];

        t_379[k] = -2.0 * ih_274[k]
                   + f_0 * lh_505[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ih_275, ih_276, ih_277, ih_278, \
                         ih_279, lh_506, lh_507, lh_508, lh_509, \
                         lh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * ih_275[k]
                   + f_0 * lh_506[k];

        t_381[k] = -2.0 * ih_276[k]
                   + f_0 * lh_507[k];

        t_382[k] = -2.0 * ih_277[k]
                   + f_0 * lh_508[k];

        t_383[k] = -2.0 * ih_278[k]
                   + f_0 * lh_509[k];

        t_384[k] = -2.0 * ih_279[k]
                   + f_0 * lh_510[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ih_280, ih_281, ih_282, ih_283, \
                         ih_284, lh_511, lh_512, lh_513, lh_514, \
                         lh_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * ih_280[k]
                   + f_0 * lh_511[k];

        t_386[k] = -2.0 * ih_281[k]
                   + f_0 * lh_512[k];

        t_387[k] = -2.0 * ih_282[k]
                   + f_0 * lh_513[k];

        t_388[k] = -2.0 * ih_283[k]
                   + f_0 * lh_514[k];

        t_389[k] = -2.0 * ih_284[k]
                   + f_0 * lh_515[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ih_285, ih_286, ih_287, ih_288, \
                         ih_289, lh_516, lh_517, lh_518, lh_519, \
                         lh_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -2.0 * ih_285[k]
                   + f_0 * lh_516[k];

        t_391[k] = -2.0 * ih_286[k]
                   + f_0 * lh_517[k];

        t_392[k] = -2.0 * ih_287[k]
                   + f_0 * lh_518[k];

        t_393[k] = -2.0 * ih_288[k]
                   + f_0 * lh_519[k];

        t_394[k] = -2.0 * ih_289[k]
                   + f_0 * lh_520[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ih_290, ih_291, ih_292, ih_293, \
                         ih_294, lh_521, lh_522, lh_523, lh_524, \
                         lh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -2.0 * ih_290[k]
                   + f_0 * lh_521[k];

        t_396[k] = -2.0 * ih_291[k]
                   + f_0 * lh_522[k];

        t_397[k] = -2.0 * ih_292[k]
                   + f_0 * lh_523[k];

        t_398[k] = -2.0 * ih_293[k]
                   + f_0 * lh_524[k];

        t_399[k] = -ih_294[k]
                   + f_0 * lh_525[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ih_295, ih_296, ih_297, ih_298, \
                         ih_299, lh_526, lh_527, lh_528, lh_529, \
                         lh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -ih_295[k]
                   + f_0 * lh_526[k];

        t_401[k] = -ih_296[k]
                   + f_0 * lh_527[k];

        t_402[k] = -ih_297[k]
                   + f_0 * lh_528[k];

        t_403[k] = -ih_298[k]
                   + f_0 * lh_529[k];

        t_404[k] = -ih_299[k]
                   + f_0 * lh_530[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ih_300, ih_301, ih_302, ih_303, \
                         ih_304, lh_531, lh_532, lh_533, lh_534, \
                         lh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -ih_300[k]
                   + f_0 * lh_531[k];

        t_406[k] = -ih_301[k]
                   + f_0 * lh_532[k];

        t_407[k] = -ih_302[k]
                   + f_0 * lh_533[k];

        t_408[k] = -ih_303[k]
                   + f_0 * lh_534[k];

        t_409[k] = -ih_304[k]
                   + f_0 * lh_535[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ih_305, ih_306, ih_307, ih_308, \
                         ih_309, lh_536, lh_537, lh_538, lh_539, \
                         lh_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -ih_305[k]
                   + f_0 * lh_536[k];

        t_411[k] = -ih_306[k]
                   + f_0 * lh_537[k];

        t_412[k] = -ih_307[k]
                   + f_0 * lh_538[k];

        t_413[k] = -ih_308[k]
                   + f_0 * lh_539[k];

        t_414[k] = -ih_309[k]
                   + f_0 * lh_540[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ih_310, ih_311, ih_312, ih_313, \
                         ih_314, lh_541, lh_542, lh_543, lh_544, \
                         lh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -ih_310[k]
                   + f_0 * lh_541[k];

        t_416[k] = -ih_311[k]
                   + f_0 * lh_542[k];

        t_417[k] = -ih_312[k]
                   + f_0 * lh_543[k];

        t_418[k] = -ih_313[k]
                   + f_0 * lh_544[k];

        t_419[k] = -ih_314[k]
                   + f_0 * lh_545[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, t_425, t_426, t_427, lh_546, \
                         lh_547, lh_548, lh_549, lh_550, lh_551, lh_552, \
                         lh_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_0 * lh_546[k];

        t_421[k] = f_0 * lh_547[k];

        t_422[k] = f_0 * lh_548[k];

        t_423[k] = f_0 * lh_549[k];

        t_424[k] = f_0 * lh_550[k];

        t_425[k] = f_0 * lh_551[k];

        t_426[k] = f_0 * lh_552[k];

        t_427[k] = f_0 * lh_553[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, t_433, t_434, t_435, lh_554, \
                         lh_555, lh_556, lh_557, lh_558, lh_559, lh_560, \
                         lh_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_0 * lh_554[k];

        t_429[k] = f_0 * lh_555[k];

        t_430[k] = f_0 * lh_556[k];

        t_431[k] = f_0 * lh_557[k];

        t_432[k] = f_0 * lh_558[k];

        t_433[k] = f_0 * lh_559[k];

        t_434[k] = f_0 * lh_560[k];

        t_435[k] = f_0 * lh_561[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, t_441, t_442, ih_315, ih_316, \
                         lh_562, lh_563, lh_564, lh_565, lh_566, lh_588, \
                         lh_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_0 * lh_562[k];

        t_437[k] = f_0 * lh_563[k];

        t_438[k] = f_0 * lh_564[k];

        t_439[k] = f_0 * lh_565[k];

        t_440[k] = f_0 * lh_566[k];

        t_441[k] = -6.0 * ih_315[k]
                   + f_0 * lh_588[k];

        t_442[k] = -6.0 * ih_316[k]
                   + f_0 * lh_589[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, t_447, ih_317, ih_318, ih_319, ih_320, \
                         ih_321, lh_590, lh_591, lh_592, lh_593, \
                         lh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = -6.0 * ih_317[k]
                   + f_0 * lh_590[k];

        t_444[k] = -6.0 * ih_318[k]
                   + f_0 * lh_591[k];

        t_445[k] = -6.0 * ih_319[k]
                   + f_0 * lh_592[k];

        t_446[k] = -6.0 * ih_320[k]
                   + f_0 * lh_593[k];

        t_447[k] = -6.0 * ih_321[k]
                   + f_0 * lh_594[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, t_452, ih_322, ih_323, ih_324, ih_325, \
                         ih_326, lh_595, lh_596, lh_597, lh_598, \
                         lh_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = -6.0 * ih_322[k]
                   + f_0 * lh_595[k];

        t_449[k] = -6.0 * ih_323[k]
                   + f_0 * lh_596[k];

        t_450[k] = -6.0 * ih_324[k]
                   + f_0 * lh_597[k];

        t_451[k] = -6.0 * ih_325[k]
                   + f_0 * lh_598[k];

        t_452[k] = -6.0 * ih_326[k]
                   + f_0 * lh_599[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, t_457, ih_327, ih_328, ih_329, ih_330, \
                         ih_331, lh_600, lh_601, lh_602, lh_603, \
                         lh_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = -6.0 * ih_327[k]
                   + f_0 * lh_600[k];

        t_454[k] = -6.0 * ih_328[k]
                   + f_0 * lh_601[k];

        t_455[k] = -6.0 * ih_329[k]
                   + f_0 * lh_602[k];

        t_456[k] = -6.0 * ih_330[k]
                   + f_0 * lh_603[k];

        t_457[k] = -6.0 * ih_331[k]
                   + f_0 * lh_604[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, t_462, ih_332, ih_333, ih_334, ih_335, \
                         ih_336, lh_605, lh_606, lh_607, lh_608, \
                         lh_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = -6.0 * ih_332[k]
                   + f_0 * lh_605[k];

        t_459[k] = -6.0 * ih_333[k]
                   + f_0 * lh_606[k];

        t_460[k] = -6.0 * ih_334[k]
                   + f_0 * lh_607[k];

        t_461[k] = -6.0 * ih_335[k]
                   + f_0 * lh_608[k];

        t_462[k] = -5.0 * ih_336[k]
                   + f_0 * lh_609[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, t_467, ih_337, ih_338, ih_339, ih_340, \
                         ih_341, lh_610, lh_611, lh_612, lh_613, \
                         lh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = -5.0 * ih_337[k]
                   + f_0 * lh_610[k];

        t_464[k] = -5.0 * ih_338[k]
                   + f_0 * lh_611[k];

        t_465[k] = -5.0 * ih_339[k]
                   + f_0 * lh_612[k];

        t_466[k] = -5.0 * ih_340[k]
                   + f_0 * lh_613[k];

        t_467[k] = -5.0 * ih_341[k]
                   + f_0 * lh_614[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, t_471, t_472, ih_342, ih_343, ih_344, ih_345, \
                         ih_346, lh_615, lh_616, lh_617, lh_618, \
                         lh_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = -5.0 * ih_342[k]
                   + f_0 * lh_615[k];

        t_469[k] = -5.0 * ih_343[k]
                   + f_0 * lh_616[k];

        t_470[k] = -5.0 * ih_344[k]
                   + f_0 * lh_617[k];

        t_471[k] = -5.0 * ih_345[k]
                   + f_0 * lh_618[k];

        t_472[k] = -5.0 * ih_346[k]
                   + f_0 * lh_619[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, t_476, t_477, ih_347, ih_348, ih_349, ih_350, \
                         ih_351, lh_620, lh_621, lh_622, lh_623, \
                         lh_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = -5.0 * ih_347[k]
                   + f_0 * lh_620[k];

        t_474[k] = -5.0 * ih_348[k]
                   + f_0 * lh_621[k];

        t_475[k] = -5.0 * ih_349[k]
                   + f_0 * lh_622[k];

        t_476[k] = -5.0 * ih_350[k]
                   + f_0 * lh_623[k];

        t_477[k] = -5.0 * ih_351[k]
                   + f_0 * lh_624[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, t_482, ih_352, ih_353, ih_354, ih_355, \
                         ih_356, lh_625, lh_626, lh_627, lh_628, \
                         lh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = -5.0 * ih_352[k]
                   + f_0 * lh_625[k];

        t_479[k] = -5.0 * ih_353[k]
                   + f_0 * lh_626[k];

        t_480[k] = -5.0 * ih_354[k]
                   + f_0 * lh_627[k];

        t_481[k] = -5.0 * ih_355[k]
                   + f_0 * lh_628[k];

        t_482[k] = -5.0 * ih_356[k]
                   + f_0 * lh_629[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, t_487, ih_357, ih_358, ih_359, ih_360, \
                         ih_361, lh_630, lh_631, lh_632, lh_633, \
                         lh_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = -4.0 * ih_357[k]
                   + f_0 * lh_630[k];

        t_484[k] = -4.0 * ih_358[k]
                   + f_0 * lh_631[k];

        t_485[k] = -4.0 * ih_359[k]
                   + f_0 * lh_632[k];

        t_486[k] = -4.0 * ih_360[k]
                   + f_0 * lh_633[k];

        t_487[k] = -4.0 * ih_361[k]
                   + f_0 * lh_634[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ih_362 = buffer.data(ih + 362);
    const auto *ih_363 = buffer.data(ih + 363);
    const auto *ih_364 = buffer.data(ih + 364);
    const auto *ih_365 = buffer.data(ih + 365);
    const auto *ih_366 = buffer.data(ih + 366);
    const auto *ih_367 = buffer.data(ih + 367);
    const auto *ih_368 = buffer.data(ih + 368);
    const auto *ih_369 = buffer.data(ih + 369);
    const auto *ih_370 = buffer.data(ih + 370);
    const auto *ih_371 = buffer.data(ih + 371);
    const auto *ih_372 = buffer.data(ih + 372);
    const auto *ih_373 = buffer.data(ih + 373);
    const auto *ih_374 = buffer.data(ih + 374);
    const auto *ih_375 = buffer.data(ih + 375);
    const auto *ih_376 = buffer.data(ih + 376);
    const auto *ih_377 = buffer.data(ih + 377);
    const auto *ih_378 = buffer.data(ih + 378);
    const auto *ih_379 = buffer.data(ih + 379);
    const auto *ih_380 = buffer.data(ih + 380);
    const auto *ih_381 = buffer.data(ih + 381);
    const auto *ih_382 = buffer.data(ih + 382);
    const auto *ih_383 = buffer.data(ih + 383);
    const auto *ih_384 = buffer.data(ih + 384);
    const auto *ih_385 = buffer.data(ih + 385);
    const auto *ih_386 = buffer.data(ih + 386);
    const auto *ih_387 = buffer.data(ih + 387);
    const auto *ih_388 = buffer.data(ih + 388);
    const auto *ih_389 = buffer.data(ih + 389);
    const auto *ih_390 = buffer.data(ih + 390);
    const auto *ih_391 = buffer.data(ih + 391);
    const auto *ih_392 = buffer.data(ih + 392);
    const auto *ih_393 = buffer.data(ih + 393);
    const auto *ih_394 = buffer.data(ih + 394);
    const auto *ih_395 = buffer.data(ih + 395);
    const auto *ih_396 = buffer.data(ih + 396);
    const auto *ih_397 = buffer.data(ih + 397);
    const auto *ih_398 = buffer.data(ih + 398);
    const auto *ih_399 = buffer.data(ih + 399);
    const auto *ih_400 = buffer.data(ih + 400);
    const auto *ih_401 = buffer.data(ih + 401);
    const auto *ih_402 = buffer.data(ih + 402);
    const auto *ih_403 = buffer.data(ih + 403);
    const auto *ih_404 = buffer.data(ih + 404);
    const auto *ih_405 = buffer.data(ih + 405);
    const auto *ih_406 = buffer.data(ih + 406);
    const auto *ih_407 = buffer.data(ih + 407);
    const auto *ih_408 = buffer.data(ih + 408);
    const auto *ih_409 = buffer.data(ih + 409);
    const auto *ih_410 = buffer.data(ih + 410);
    const auto *ih_411 = buffer.data(ih + 411);
    const auto *ih_412 = buffer.data(ih + 412);
    const auto *ih_413 = buffer.data(ih + 413);
    const auto *ih_414 = buffer.data(ih + 414);
    const auto *ih_415 = buffer.data(ih + 415);
    const auto *ih_416 = buffer.data(ih + 416);
    const auto *ih_417 = buffer.data(ih + 417);
    const auto *ih_418 = buffer.data(ih + 418);
    const auto *ih_419 = buffer.data(ih + 419);
    const auto *ih_420 = buffer.data(ih + 420);
    const auto *ih_421 = buffer.data(ih + 421);
    const auto *ih_422 = buffer.data(ih + 422);
    const auto *ih_423 = buffer.data(ih + 423);
    const auto *ih_424 = buffer.data(ih + 424);
    const auto *ih_425 = buffer.data(ih + 425);
    const auto *ih_426 = buffer.data(ih + 426);
    const auto *ih_427 = buffer.data(ih + 427);
    const auto *ih_428 = buffer.data(ih + 428);
    const auto *ih_429 = buffer.data(ih + 429);
    const auto *ih_430 = buffer.data(ih + 430);
    const auto *ih_431 = buffer.data(ih + 431);
    const auto *ih_432 = buffer.data(ih + 432);
    const auto *ih_433 = buffer.data(ih + 433);
    const auto *ih_434 = buffer.data(ih + 434);
    const auto *ih_435 = buffer.data(ih + 435);
    const auto *ih_436 = buffer.data(ih + 436);
    const auto *ih_437 = buffer.data(ih + 437);
    const auto *ih_438 = buffer.data(ih + 438);
    const auto *ih_439 = buffer.data(ih + 439);
    const auto *ih_440 = buffer.data(ih + 440);
    const auto *ih_441 = buffer.data(ih + 441);
    const auto *ih_442 = buffer.data(ih + 442);
    const auto *ih_443 = buffer.data(ih + 443);
    const auto *ih_444 = buffer.data(ih + 444);
    const auto *ih_445 = buffer.data(ih + 445);
    const auto *ih_446 = buffer.data(ih + 446);
    const auto *ih_447 = buffer.data(ih + 447);
    const auto *ih_448 = buffer.data(ih + 448);
    const auto *ih_449 = buffer.data(ih + 449);
    const auto *ih_450 = buffer.data(ih + 450);
    const auto *ih_451 = buffer.data(ih + 451);
    const auto *ih_452 = buffer.data(ih + 452);
    const auto *ih_453 = buffer.data(ih + 453);
    const auto *ih_454 = buffer.data(ih + 454);
    const auto *ih_455 = buffer.data(ih + 455);
    const auto *ih_456 = buffer.data(ih + 456);
    const auto *ih_457 = buffer.data(ih + 457);
    const auto *ih_458 = buffer.data(ih + 458);
    const auto *ih_459 = buffer.data(ih + 459);
    const auto *ih_460 = buffer.data(ih + 460);
    const auto *ih_461 = buffer.data(ih + 461);
    const auto *ih_462 = buffer.data(ih + 462);
    const auto *ih_463 = buffer.data(ih + 463);
    const auto *ih_464 = buffer.data(ih + 464);
    const auto *ih_465 = buffer.data(ih + 465);
    const auto *ih_466 = buffer.data(ih + 466);
    const auto *ih_467 = buffer.data(ih + 467);
    const auto *ih_468 = buffer.data(ih + 468);
    const auto *ih_469 = buffer.data(ih + 469);
    const auto *ih_470 = buffer.data(ih + 470);
    const auto *ih_471 = buffer.data(ih + 471);
    const auto *ih_472 = buffer.data(ih + 472);
    const auto *ih_473 = buffer.data(ih + 473);
    const auto *ih_474 = buffer.data(ih + 474);
    const auto *ih_475 = buffer.data(ih + 475);
    const auto *ih_476 = buffer.data(ih + 476);
    const auto *ih_477 = buffer.data(ih + 477);
    const auto *ih_478 = buffer.data(ih + 478);
    const auto *ih_479 = buffer.data(ih + 479);
    const auto *ih_480 = buffer.data(ih + 480);
    const auto *ih_481 = buffer.data(ih + 481);
    const auto *ih_482 = buffer.data(ih + 482);
    const auto *ih_483 = buffer.data(ih + 483);
    const auto *ih_484 = buffer.data(ih + 484);
    const auto *ih_485 = buffer.data(ih + 485);
    const auto *ih_486 = buffer.data(ih + 486);
    const auto *ih_487 = buffer.data(ih + 487);
    const auto *ih_488 = buffer.data(ih + 488);
    const auto *ih_489 = buffer.data(ih + 489);
    const auto *ih_490 = buffer.data(ih + 490);
    const auto *ih_491 = buffer.data(ih + 491);
    const auto *ih_492 = buffer.data(ih + 492);
    const auto *ih_493 = buffer.data(ih + 493);
    const auto *ih_494 = buffer.data(ih + 494);
    const auto *ih_495 = buffer.data(ih + 495);
    const auto *ih_496 = buffer.data(ih + 496);
    const auto *ih_497 = buffer.data(ih + 497);
    const auto *ih_498 = buffer.data(ih + 498);

    const auto *lh_635 = buffer.data(lh + 635);
    const auto *lh_636 = buffer.data(lh + 636);
    const auto *lh_637 = buffer.data(lh + 637);
    const auto *lh_638 = buffer.data(lh + 638);
    const auto *lh_639 = buffer.data(lh + 639);
    const auto *lh_640 = buffer.data(lh + 640);
    const auto *lh_641 = buffer.data(lh + 641);
    const auto *lh_642 = buffer.data(lh + 642);
    const auto *lh_643 = buffer.data(lh + 643);
    const auto *lh_644 = buffer.data(lh + 644);
    const auto *lh_645 = buffer.data(lh + 645);
    const auto *lh_646 = buffer.data(lh + 646);
    const auto *lh_647 = buffer.data(lh + 647);
    const auto *lh_648 = buffer.data(lh + 648);
    const auto *lh_649 = buffer.data(lh + 649);
    const auto *lh_650 = buffer.data(lh + 650);
    const auto *lh_651 = buffer.data(lh + 651);
    const auto *lh_652 = buffer.data(lh + 652);
    const auto *lh_653 = buffer.data(lh + 653);
    const auto *lh_654 = buffer.data(lh + 654);
    const auto *lh_655 = buffer.data(lh + 655);
    const auto *lh_656 = buffer.data(lh + 656);
    const auto *lh_657 = buffer.data(lh + 657);
    const auto *lh_658 = buffer.data(lh + 658);
    const auto *lh_659 = buffer.data(lh + 659);
    const auto *lh_660 = buffer.data(lh + 660);
    const auto *lh_661 = buffer.data(lh + 661);
    const auto *lh_662 = buffer.data(lh + 662);
    const auto *lh_663 = buffer.data(lh + 663);
    const auto *lh_664 = buffer.data(lh + 664);
    const auto *lh_665 = buffer.data(lh + 665);
    const auto *lh_666 = buffer.data(lh + 666);
    const auto *lh_667 = buffer.data(lh + 667);
    const auto *lh_668 = buffer.data(lh + 668);
    const auto *lh_669 = buffer.data(lh + 669);
    const auto *lh_670 = buffer.data(lh + 670);
    const auto *lh_671 = buffer.data(lh + 671);
    const auto *lh_672 = buffer.data(lh + 672);
    const auto *lh_673 = buffer.data(lh + 673);
    const auto *lh_674 = buffer.data(lh + 674);
    const auto *lh_675 = buffer.data(lh + 675);
    const auto *lh_676 = buffer.data(lh + 676);
    const auto *lh_677 = buffer.data(lh + 677);
    const auto *lh_678 = buffer.data(lh + 678);
    const auto *lh_679 = buffer.data(lh + 679);
    const auto *lh_680 = buffer.data(lh + 680);
    const auto *lh_681 = buffer.data(lh + 681);
    const auto *lh_682 = buffer.data(lh + 682);
    const auto *lh_683 = buffer.data(lh + 683);
    const auto *lh_684 = buffer.data(lh + 684);
    const auto *lh_685 = buffer.data(lh + 685);
    const auto *lh_686 = buffer.data(lh + 686);
    const auto *lh_687 = buffer.data(lh + 687);
    const auto *lh_688 = buffer.data(lh + 688);
    const auto *lh_689 = buffer.data(lh + 689);
    const auto *lh_690 = buffer.data(lh + 690);
    const auto *lh_691 = buffer.data(lh + 691);
    const auto *lh_692 = buffer.data(lh + 692);
    const auto *lh_693 = buffer.data(lh + 693);
    const auto *lh_694 = buffer.data(lh + 694);
    const auto *lh_695 = buffer.data(lh + 695);
    const auto *lh_696 = buffer.data(lh + 696);
    const auto *lh_697 = buffer.data(lh + 697);
    const auto *lh_698 = buffer.data(lh + 698);
    const auto *lh_699 = buffer.data(lh + 699);
    const auto *lh_700 = buffer.data(lh + 700);
    const auto *lh_701 = buffer.data(lh + 701);
    const auto *lh_702 = buffer.data(lh + 702);
    const auto *lh_703 = buffer.data(lh + 703);
    const auto *lh_704 = buffer.data(lh + 704);
    const auto *lh_705 = buffer.data(lh + 705);
    const auto *lh_706 = buffer.data(lh + 706);
    const auto *lh_707 = buffer.data(lh + 707);
    const auto *lh_708 = buffer.data(lh + 708);
    const auto *lh_709 = buffer.data(lh + 709);
    const auto *lh_710 = buffer.data(lh + 710);
    const auto *lh_711 = buffer.data(lh + 711);
    const auto *lh_712 = buffer.data(lh + 712);
    const auto *lh_713 = buffer.data(lh + 713);
    const auto *lh_714 = buffer.data(lh + 714);
    const auto *lh_715 = buffer.data(lh + 715);
    const auto *lh_716 = buffer.data(lh + 716);
    const auto *lh_717 = buffer.data(lh + 717);
    const auto *lh_718 = buffer.data(lh + 718);
    const auto *lh_719 = buffer.data(lh + 719);
    const auto *lh_720 = buffer.data(lh + 720);
    const auto *lh_721 = buffer.data(lh + 721);
    const auto *lh_722 = buffer.data(lh + 722);
    const auto *lh_723 = buffer.data(lh + 723);
    const auto *lh_724 = buffer.data(lh + 724);
    const auto *lh_725 = buffer.data(lh + 725);
    const auto *lh_726 = buffer.data(lh + 726);
    const auto *lh_727 = buffer.data(lh + 727);
    const auto *lh_728 = buffer.data(lh + 728);
    const auto *lh_729 = buffer.data(lh + 729);
    const auto *lh_730 = buffer.data(lh + 730);
    const auto *lh_731 = buffer.data(lh + 731);
    const auto *lh_732 = buffer.data(lh + 732);
    const auto *lh_733 = buffer.data(lh + 733);
    const auto *lh_734 = buffer.data(lh + 734);
    const auto *lh_756 = buffer.data(lh + 756);
    const auto *lh_757 = buffer.data(lh + 757);
    const auto *lh_758 = buffer.data(lh + 758);
    const auto *lh_759 = buffer.data(lh + 759);
    const auto *lh_760 = buffer.data(lh + 760);
    const auto *lh_761 = buffer.data(lh + 761);
    const auto *lh_762 = buffer.data(lh + 762);
    const auto *lh_763 = buffer.data(lh + 763);
    const auto *lh_764 = buffer.data(lh + 764);
    const auto *lh_765 = buffer.data(lh + 765);
    const auto *lh_766 = buffer.data(lh + 766);
    const auto *lh_767 = buffer.data(lh + 767);
    const auto *lh_768 = buffer.data(lh + 768);
    const auto *lh_769 = buffer.data(lh + 769);
    const auto *lh_770 = buffer.data(lh + 770);
    const auto *lh_771 = buffer.data(lh + 771);
    const auto *lh_772 = buffer.data(lh + 772);
    const auto *lh_773 = buffer.data(lh + 773);
    const auto *lh_774 = buffer.data(lh + 774);
    const auto *lh_775 = buffer.data(lh + 775);
    const auto *lh_776 = buffer.data(lh + 776);
    const auto *lh_777 = buffer.data(lh + 777);
    const auto *lh_778 = buffer.data(lh + 778);
    const auto *lh_779 = buffer.data(lh + 779);
    const auto *lh_780 = buffer.data(lh + 780);
    const auto *lh_781 = buffer.data(lh + 781);
    const auto *lh_782 = buffer.data(lh + 782);
    const auto *lh_783 = buffer.data(lh + 783);
    const auto *lh_784 = buffer.data(lh + 784);
    const auto *lh_785 = buffer.data(lh + 785);
    const auto *lh_786 = buffer.data(lh + 786);
    const auto *lh_787 = buffer.data(lh + 787);
    const auto *lh_788 = buffer.data(lh + 788);
    const auto *lh_789 = buffer.data(lh + 789);
    const auto *lh_790 = buffer.data(lh + 790);
    const auto *lh_791 = buffer.data(lh + 791);
    const auto *lh_792 = buffer.data(lh + 792);
    const auto *lh_793 = buffer.data(lh + 793);
    const auto *lh_794 = buffer.data(lh + 794);
    const auto *lh_795 = buffer.data(lh + 795);
    const auto *lh_796 = buffer.data(lh + 796);
    const auto *lh_797 = buffer.data(lh + 797);
    const auto *lh_798 = buffer.data(lh + 798);
    const auto *lh_799 = buffer.data(lh + 799);
    const auto *lh_800 = buffer.data(lh + 800);
    const auto *lh_801 = buffer.data(lh + 801);
    const auto *lh_802 = buffer.data(lh + 802);
    const auto *lh_803 = buffer.data(lh + 803);
    const auto *lh_804 = buffer.data(lh + 804);
    const auto *lh_805 = buffer.data(lh + 805);
    const auto *lh_806 = buffer.data(lh + 806);
    const auto *lh_807 = buffer.data(lh + 807);
    const auto *lh_808 = buffer.data(lh + 808);
    const auto *lh_809 = buffer.data(lh + 809);
    const auto *lh_810 = buffer.data(lh + 810);
    const auto *lh_811 = buffer.data(lh + 811);
    const auto *lh_812 = buffer.data(lh + 812);
    const auto *lh_813 = buffer.data(lh + 813);

#pragma omp simd aligned(t_488, t_489, t_490, t_491, t_492, ih_362, ih_363, ih_364, ih_365, \
                         ih_366, lh_635, lh_636, lh_637, lh_638, \
                         lh_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = -4.0 * ih_362[k]
                   + f_0 * lh_635[k];

        t_489[k] = -4.0 * ih_363[k]
                   + f_0 * lh_636[k];

        t_490[k] = -4.0 * ih_364[k]
                   + f_0 * lh_637[k];

        t_491[k] = -4.0 * ih_365[k]
                   + f_0 * lh_638[k];

        t_492[k] = -4.0 * ih_366[k]
                   + f_0 * lh_639[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, t_496, t_497, ih_367, ih_368, ih_369, ih_370, \
                         ih_371, lh_640, lh_641, lh_642, lh_643, \
                         lh_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = -4.0 * ih_367[k]
                   + f_0 * lh_640[k];

        t_494[k] = -4.0 * ih_368[k]
                   + f_0 * lh_641[k];

        t_495[k] = -4.0 * ih_369[k]
                   + f_0 * lh_642[k];

        t_496[k] = -4.0 * ih_370[k]
                   + f_0 * lh_643[k];

        t_497[k] = -4.0 * ih_371[k]
                   + f_0 * lh_644[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, t_501, t_502, ih_372, ih_373, ih_374, ih_375, \
                         ih_376, lh_645, lh_646, lh_647, lh_648, \
                         lh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = -4.0 * ih_372[k]
                   + f_0 * lh_645[k];

        t_499[k] = -4.0 * ih_373[k]
                   + f_0 * lh_646[k];

        t_500[k] = -4.0 * ih_374[k]
                   + f_0 * lh_647[k];

        t_501[k] = -4.0 * ih_375[k]
                   + f_0 * lh_648[k];

        t_502[k] = -4.0 * ih_376[k]
                   + f_0 * lh_649[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, t_507, ih_377, ih_378, ih_379, ih_380, \
                         ih_381, lh_650, lh_651, lh_652, lh_653, \
                         lh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = -4.0 * ih_377[k]
                   + f_0 * lh_650[k];

        t_504[k] = -3.0 * ih_378[k]
                   + f_0 * lh_651[k];

        t_505[k] = -3.0 * ih_379[k]
                   + f_0 * lh_652[k];

        t_506[k] = -3.0 * ih_380[k]
                   + f_0 * lh_653[k];

        t_507[k] = -3.0 * ih_381[k]
                   + f_0 * lh_654[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, t_512, ih_382, ih_383, ih_384, ih_385, \
                         ih_386, lh_655, lh_656, lh_657, lh_658, \
                         lh_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = -3.0 * ih_382[k]
                   + f_0 * lh_655[k];

        t_509[k] = -3.0 * ih_383[k]
                   + f_0 * lh_656[k];

        t_510[k] = -3.0 * ih_384[k]
                   + f_0 * lh_657[k];

        t_511[k] = -3.0 * ih_385[k]
                   + f_0 * lh_658[k];

        t_512[k] = -3.0 * ih_386[k]
                   + f_0 * lh_659[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, t_516, t_517, ih_387, ih_388, ih_389, ih_390, \
                         ih_391, lh_660, lh_661, lh_662, lh_663, \
                         lh_664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = -3.0 * ih_387[k]
                   + f_0 * lh_660[k];

        t_514[k] = -3.0 * ih_388[k]
                   + f_0 * lh_661[k];

        t_515[k] = -3.0 * ih_389[k]
                   + f_0 * lh_662[k];

        t_516[k] = -3.0 * ih_390[k]
                   + f_0 * lh_663[k];

        t_517[k] = -3.0 * ih_391[k]
                   + f_0 * lh_664[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, t_522, ih_392, ih_393, ih_394, ih_395, \
                         ih_396, lh_665, lh_666, lh_667, lh_668, \
                         lh_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = -3.0 * ih_392[k]
                   + f_0 * lh_665[k];

        t_519[k] = -3.0 * ih_393[k]
                   + f_0 * lh_666[k];

        t_520[k] = -3.0 * ih_394[k]
                   + f_0 * lh_667[k];

        t_521[k] = -3.0 * ih_395[k]
                   + f_0 * lh_668[k];

        t_522[k] = -3.0 * ih_396[k]
                   + f_0 * lh_669[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, t_527, ih_397, ih_398, ih_399, ih_400, \
                         ih_401, lh_670, lh_671, lh_672, lh_673, \
                         lh_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = -3.0 * ih_397[k]
                   + f_0 * lh_670[k];

        t_524[k] = -3.0 * ih_398[k]
                   + f_0 * lh_671[k];

        t_525[k] = -2.0 * ih_399[k]
                   + f_0 * lh_672[k];

        t_526[k] = -2.0 * ih_400[k]
                   + f_0 * lh_673[k];

        t_527[k] = -2.0 * ih_401[k]
                   + f_0 * lh_674[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, ih_402, ih_403, ih_404, ih_405, \
                         ih_406, lh_675, lh_676, lh_677, lh_678, \
                         lh_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = -2.0 * ih_402[k]
                   + f_0 * lh_675[k];

        t_529[k] = -2.0 * ih_403[k]
                   + f_0 * lh_676[k];

        t_530[k] = -2.0 * ih_404[k]
                   + f_0 * lh_677[k];

        t_531[k] = -2.0 * ih_405[k]
                   + f_0 * lh_678[k];

        t_532[k] = -2.0 * ih_406[k]
                   + f_0 * lh_679[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, t_537, ih_407, ih_408, ih_409, ih_410, \
                         ih_411, lh_680, lh_681, lh_682, lh_683, \
                         lh_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = -2.0 * ih_407[k]
                   + f_0 * lh_680[k];

        t_534[k] = -2.0 * ih_408[k]
                   + f_0 * lh_681[k];

        t_535[k] = -2.0 * ih_409[k]
                   + f_0 * lh_682[k];

        t_536[k] = -2.0 * ih_410[k]
                   + f_0 * lh_683[k];

        t_537[k] = -2.0 * ih_411[k]
                   + f_0 * lh_684[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, t_542, ih_412, ih_413, ih_414, ih_415, \
                         ih_416, lh_685, lh_686, lh_687, lh_688, \
                         lh_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = -2.0 * ih_412[k]
                   + f_0 * lh_685[k];

        t_539[k] = -2.0 * ih_413[k]
                   + f_0 * lh_686[k];

        t_540[k] = -2.0 * ih_414[k]
                   + f_0 * lh_687[k];

        t_541[k] = -2.0 * ih_415[k]
                   + f_0 * lh_688[k];

        t_542[k] = -2.0 * ih_416[k]
                   + f_0 * lh_689[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, t_547, ih_417, ih_418, ih_419, ih_420, \
                         ih_421, lh_690, lh_691, lh_692, lh_693, \
                         lh_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = -2.0 * ih_417[k]
                   + f_0 * lh_690[k];

        t_544[k] = -2.0 * ih_418[k]
                   + f_0 * lh_691[k];

        t_545[k] = -2.0 * ih_419[k]
                   + f_0 * lh_692[k];

        t_546[k] = -ih_420[k]
                   + f_0 * lh_693[k];

        t_547[k] = -ih_421[k]
                   + f_0 * lh_694[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, ih_422, ih_423, ih_424, ih_425, \
                         ih_426, lh_695, lh_696, lh_697, lh_698, \
                         lh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = -ih_422[k]
                   + f_0 * lh_695[k];

        t_549[k] = -ih_423[k]
                   + f_0 * lh_696[k];

        t_550[k] = -ih_424[k]
                   + f_0 * lh_697[k];

        t_551[k] = -ih_425[k]
                   + f_0 * lh_698[k];

        t_552[k] = -ih_426[k]
                   + f_0 * lh_699[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, t_557, ih_427, ih_428, ih_429, ih_430, \
                         ih_431, lh_700, lh_701, lh_702, lh_703, \
                         lh_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = -ih_427[k]
                   + f_0 * lh_700[k];

        t_554[k] = -ih_428[k]
                   + f_0 * lh_701[k];

        t_555[k] = -ih_429[k]
                   + f_0 * lh_702[k];

        t_556[k] = -ih_430[k]
                   + f_0 * lh_703[k];

        t_557[k] = -ih_431[k]
                   + f_0 * lh_704[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, t_562, ih_432, ih_433, ih_434, ih_435, \
                         ih_436, lh_705, lh_706, lh_707, lh_708, \
                         lh_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = -ih_432[k]
                   + f_0 * lh_705[k];

        t_559[k] = -ih_433[k]
                   + f_0 * lh_706[k];

        t_560[k] = -ih_434[k]
                   + f_0 * lh_707[k];

        t_561[k] = -ih_435[k]
                   + f_0 * lh_708[k];

        t_562[k] = -ih_436[k]
                   + f_0 * lh_709[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, t_567, t_568, ih_437, ih_438, ih_439, \
                         ih_440, lh_710, lh_711, lh_712, lh_713, lh_714, \
                         lh_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = -ih_437[k]
                   + f_0 * lh_710[k];

        t_564[k] = -ih_438[k]
                   + f_0 * lh_711[k];

        t_565[k] = -ih_439[k]
                   + f_0 * lh_712[k];

        t_566[k] = -ih_440[k]
                   + f_0 * lh_713[k];

        t_567[k] = f_0 * lh_714[k];

        t_568[k] = f_0 * lh_715[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, t_573, t_574, t_575, t_576, lh_716, \
                         lh_717, lh_718, lh_719, lh_720, lh_721, lh_722, \
                         lh_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_0 * lh_716[k];

        t_570[k] = f_0 * lh_717[k];

        t_571[k] = f_0 * lh_718[k];

        t_572[k] = f_0 * lh_719[k];

        t_573[k] = f_0 * lh_720[k];

        t_574[k] = f_0 * lh_721[k];

        t_575[k] = f_0 * lh_722[k];

        t_576[k] = f_0 * lh_723[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, t_582, t_583, t_584, lh_724, \
                         lh_725, lh_726, lh_727, lh_728, lh_729, lh_730, \
                         lh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_0 * lh_724[k];

        t_578[k] = f_0 * lh_725[k];

        t_579[k] = f_0 * lh_726[k];

        t_580[k] = f_0 * lh_727[k];

        t_581[k] = f_0 * lh_728[k];

        t_582[k] = f_0 * lh_729[k];

        t_583[k] = f_0 * lh_730[k];

        t_584[k] = f_0 * lh_731[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, t_590, ih_441, ih_442, ih_443, \
                         lh_732, lh_733, lh_734, lh_756, lh_757, \
                         lh_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_0 * lh_732[k];

        t_586[k] = f_0 * lh_733[k];

        t_587[k] = f_0 * lh_734[k];

        t_588[k] = -7.0 * ih_441[k]
                   + f_0 * lh_756[k];

        t_589[k] = -7.0 * ih_442[k]
                   + f_0 * lh_757[k];

        t_590[k] = -7.0 * ih_443[k]
                   + f_0 * lh_758[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, t_595, ih_444, ih_445, ih_446, ih_447, \
                         ih_448, lh_759, lh_760, lh_761, lh_762, \
                         lh_763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = -7.0 * ih_444[k]
                   + f_0 * lh_759[k];

        t_592[k] = -7.0 * ih_445[k]
                   + f_0 * lh_760[k];

        t_593[k] = -7.0 * ih_446[k]
                   + f_0 * lh_761[k];

        t_594[k] = -7.0 * ih_447[k]
                   + f_0 * lh_762[k];

        t_595[k] = -7.0 * ih_448[k]
                   + f_0 * lh_763[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, t_600, ih_449, ih_450, ih_451, ih_452, \
                         ih_453, lh_764, lh_765, lh_766, lh_767, \
                         lh_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = -7.0 * ih_449[k]
                   + f_0 * lh_764[k];

        t_597[k] = -7.0 * ih_450[k]
                   + f_0 * lh_765[k];

        t_598[k] = -7.0 * ih_451[k]
                   + f_0 * lh_766[k];

        t_599[k] = -7.0 * ih_452[k]
                   + f_0 * lh_767[k];

        t_600[k] = -7.0 * ih_453[k]
                   + f_0 * lh_768[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, t_604, t_605, ih_454, ih_455, ih_456, ih_457, \
                         ih_458, lh_769, lh_770, lh_771, lh_772, \
                         lh_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = -7.0 * ih_454[k]
                   + f_0 * lh_769[k];

        t_602[k] = -7.0 * ih_455[k]
                   + f_0 * lh_770[k];

        t_603[k] = -7.0 * ih_456[k]
                   + f_0 * lh_771[k];

        t_604[k] = -7.0 * ih_457[k]
                   + f_0 * lh_772[k];

        t_605[k] = -7.0 * ih_458[k]
                   + f_0 * lh_773[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, t_609, t_610, ih_459, ih_460, ih_461, ih_462, \
                         ih_463, lh_774, lh_775, lh_776, lh_777, \
                         lh_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = -7.0 * ih_459[k]
                   + f_0 * lh_774[k];

        t_607[k] = -7.0 * ih_460[k]
                   + f_0 * lh_775[k];

        t_608[k] = -7.0 * ih_461[k]
                   + f_0 * lh_776[k];

        t_609[k] = -6.0 * ih_462[k]
                   + f_0 * lh_777[k];

        t_610[k] = -6.0 * ih_463[k]
                   + f_0 * lh_778[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, t_614, t_615, ih_464, ih_465, ih_466, ih_467, \
                         ih_468, lh_779, lh_780, lh_781, lh_782, \
                         lh_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = -6.0 * ih_464[k]
                   + f_0 * lh_779[k];

        t_612[k] = -6.0 * ih_465[k]
                   + f_0 * lh_780[k];

        t_613[k] = -6.0 * ih_466[k]
                   + f_0 * lh_781[k];

        t_614[k] = -6.0 * ih_467[k]
                   + f_0 * lh_782[k];

        t_615[k] = -6.0 * ih_468[k]
                   + f_0 * lh_783[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, t_620, ih_469, ih_470, ih_471, ih_472, \
                         ih_473, lh_784, lh_785, lh_786, lh_787, \
                         lh_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = -6.0 * ih_469[k]
                   + f_0 * lh_784[k];

        t_617[k] = -6.0 * ih_470[k]
                   + f_0 * lh_785[k];

        t_618[k] = -6.0 * ih_471[k]
                   + f_0 * lh_786[k];

        t_619[k] = -6.0 * ih_472[k]
                   + f_0 * lh_787[k];

        t_620[k] = -6.0 * ih_473[k]
                   + f_0 * lh_788[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, t_624, t_625, ih_474, ih_475, ih_476, ih_477, \
                         ih_478, lh_789, lh_790, lh_791, lh_792, \
                         lh_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = -6.0 * ih_474[k]
                   + f_0 * lh_789[k];

        t_622[k] = -6.0 * ih_475[k]
                   + f_0 * lh_790[k];

        t_623[k] = -6.0 * ih_476[k]
                   + f_0 * lh_791[k];

        t_624[k] = -6.0 * ih_477[k]
                   + f_0 * lh_792[k];

        t_625[k] = -6.0 * ih_478[k]
                   + f_0 * lh_793[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, t_630, ih_479, ih_480, ih_481, ih_482, \
                         ih_483, lh_794, lh_795, lh_796, lh_797, \
                         lh_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = -6.0 * ih_479[k]
                   + f_0 * lh_794[k];

        t_627[k] = -6.0 * ih_480[k]
                   + f_0 * lh_795[k];

        t_628[k] = -6.0 * ih_481[k]
                   + f_0 * lh_796[k];

        t_629[k] = -6.0 * ih_482[k]
                   + f_0 * lh_797[k];

        t_630[k] = -5.0 * ih_483[k]
                   + f_0 * lh_798[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, t_635, ih_484, ih_485, ih_486, ih_487, \
                         ih_488, lh_799, lh_800, lh_801, lh_802, \
                         lh_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = -5.0 * ih_484[k]
                   + f_0 * lh_799[k];

        t_632[k] = -5.0 * ih_485[k]
                   + f_0 * lh_800[k];

        t_633[k] = -5.0 * ih_486[k]
                   + f_0 * lh_801[k];

        t_634[k] = -5.0 * ih_487[k]
                   + f_0 * lh_802[k];

        t_635[k] = -5.0 * ih_488[k]
                   + f_0 * lh_803[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, t_640, ih_489, ih_490, ih_491, ih_492, \
                         ih_493, lh_804, lh_805, lh_806, lh_807, \
                         lh_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = -5.0 * ih_489[k]
                   + f_0 * lh_804[k];

        t_637[k] = -5.0 * ih_490[k]
                   + f_0 * lh_805[k];

        t_638[k] = -5.0 * ih_491[k]
                   + f_0 * lh_806[k];

        t_639[k] = -5.0 * ih_492[k]
                   + f_0 * lh_807[k];

        t_640[k] = -5.0 * ih_493[k]
                   + f_0 * lh_808[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, t_645, ih_494, ih_495, ih_496, ih_497, \
                         ih_498, lh_809, lh_810, lh_811, lh_812, \
                         lh_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = -5.0 * ih_494[k]
                   + f_0 * lh_809[k];

        t_642[k] = -5.0 * ih_495[k]
                   + f_0 * lh_810[k];

        t_643[k] = -5.0 * ih_496[k]
                   + f_0 * lh_811[k];

        t_644[k] = -5.0 * ih_497[k]
                   + f_0 * lh_812[k];

        t_645[k] = -5.0 * ih_498[k]
                   + f_0 * lh_813[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ih_499 = buffer.data(ih + 499);
    const auto *ih_500 = buffer.data(ih + 500);
    const auto *ih_501 = buffer.data(ih + 501);
    const auto *ih_502 = buffer.data(ih + 502);
    const auto *ih_503 = buffer.data(ih + 503);
    const auto *ih_504 = buffer.data(ih + 504);
    const auto *ih_505 = buffer.data(ih + 505);
    const auto *ih_506 = buffer.data(ih + 506);
    const auto *ih_507 = buffer.data(ih + 507);
    const auto *ih_508 = buffer.data(ih + 508);
    const auto *ih_509 = buffer.data(ih + 509);
    const auto *ih_510 = buffer.data(ih + 510);
    const auto *ih_511 = buffer.data(ih + 511);
    const auto *ih_512 = buffer.data(ih + 512);
    const auto *ih_513 = buffer.data(ih + 513);
    const auto *ih_514 = buffer.data(ih + 514);
    const auto *ih_515 = buffer.data(ih + 515);
    const auto *ih_516 = buffer.data(ih + 516);
    const auto *ih_517 = buffer.data(ih + 517);
    const auto *ih_518 = buffer.data(ih + 518);
    const auto *ih_519 = buffer.data(ih + 519);
    const auto *ih_520 = buffer.data(ih + 520);
    const auto *ih_521 = buffer.data(ih + 521);
    const auto *ih_522 = buffer.data(ih + 522);
    const auto *ih_523 = buffer.data(ih + 523);
    const auto *ih_524 = buffer.data(ih + 524);
    const auto *ih_525 = buffer.data(ih + 525);
    const auto *ih_526 = buffer.data(ih + 526);
    const auto *ih_527 = buffer.data(ih + 527);
    const auto *ih_528 = buffer.data(ih + 528);
    const auto *ih_529 = buffer.data(ih + 529);
    const auto *ih_530 = buffer.data(ih + 530);
    const auto *ih_531 = buffer.data(ih + 531);
    const auto *ih_532 = buffer.data(ih + 532);
    const auto *ih_533 = buffer.data(ih + 533);
    const auto *ih_534 = buffer.data(ih + 534);
    const auto *ih_535 = buffer.data(ih + 535);
    const auto *ih_536 = buffer.data(ih + 536);
    const auto *ih_537 = buffer.data(ih + 537);
    const auto *ih_538 = buffer.data(ih + 538);
    const auto *ih_539 = buffer.data(ih + 539);
    const auto *ih_540 = buffer.data(ih + 540);
    const auto *ih_541 = buffer.data(ih + 541);
    const auto *ih_542 = buffer.data(ih + 542);
    const auto *ih_543 = buffer.data(ih + 543);
    const auto *ih_544 = buffer.data(ih + 544);
    const auto *ih_545 = buffer.data(ih + 545);
    const auto *ih_546 = buffer.data(ih + 546);
    const auto *ih_547 = buffer.data(ih + 547);
    const auto *ih_548 = buffer.data(ih + 548);
    const auto *ih_549 = buffer.data(ih + 549);
    const auto *ih_550 = buffer.data(ih + 550);
    const auto *ih_551 = buffer.data(ih + 551);
    const auto *ih_552 = buffer.data(ih + 552);
    const auto *ih_553 = buffer.data(ih + 553);
    const auto *ih_554 = buffer.data(ih + 554);
    const auto *ih_555 = buffer.data(ih + 555);
    const auto *ih_556 = buffer.data(ih + 556);
    const auto *ih_557 = buffer.data(ih + 557);
    const auto *ih_558 = buffer.data(ih + 558);
    const auto *ih_559 = buffer.data(ih + 559);
    const auto *ih_560 = buffer.data(ih + 560);
    const auto *ih_561 = buffer.data(ih + 561);
    const auto *ih_562 = buffer.data(ih + 562);
    const auto *ih_563 = buffer.data(ih + 563);
    const auto *ih_564 = buffer.data(ih + 564);
    const auto *ih_565 = buffer.data(ih + 565);
    const auto *ih_566 = buffer.data(ih + 566);
    const auto *ih_567 = buffer.data(ih + 567);
    const auto *ih_568 = buffer.data(ih + 568);
    const auto *ih_569 = buffer.data(ih + 569);
    const auto *ih_570 = buffer.data(ih + 570);
    const auto *ih_571 = buffer.data(ih + 571);
    const auto *ih_572 = buffer.data(ih + 572);
    const auto *ih_573 = buffer.data(ih + 573);
    const auto *ih_574 = buffer.data(ih + 574);
    const auto *ih_575 = buffer.data(ih + 575);
    const auto *ih_576 = buffer.data(ih + 576);
    const auto *ih_577 = buffer.data(ih + 577);
    const auto *ih_578 = buffer.data(ih + 578);
    const auto *ih_579 = buffer.data(ih + 579);
    const auto *ih_580 = buffer.data(ih + 580);
    const auto *ih_581 = buffer.data(ih + 581);
    const auto *ih_582 = buffer.data(ih + 582);
    const auto *ih_583 = buffer.data(ih + 583);
    const auto *ih_584 = buffer.data(ih + 584);
    const auto *ih_585 = buffer.data(ih + 585);
    const auto *ih_586 = buffer.data(ih + 586);
    const auto *ih_587 = buffer.data(ih + 587);

    const auto *lh_814 = buffer.data(lh + 814);
    const auto *lh_815 = buffer.data(lh + 815);
    const auto *lh_816 = buffer.data(lh + 816);
    const auto *lh_817 = buffer.data(lh + 817);
    const auto *lh_818 = buffer.data(lh + 818);
    const auto *lh_819 = buffer.data(lh + 819);
    const auto *lh_820 = buffer.data(lh + 820);
    const auto *lh_821 = buffer.data(lh + 821);
    const auto *lh_822 = buffer.data(lh + 822);
    const auto *lh_823 = buffer.data(lh + 823);
    const auto *lh_824 = buffer.data(lh + 824);
    const auto *lh_825 = buffer.data(lh + 825);
    const auto *lh_826 = buffer.data(lh + 826);
    const auto *lh_827 = buffer.data(lh + 827);
    const auto *lh_828 = buffer.data(lh + 828);
    const auto *lh_829 = buffer.data(lh + 829);
    const auto *lh_830 = buffer.data(lh + 830);
    const auto *lh_831 = buffer.data(lh + 831);
    const auto *lh_832 = buffer.data(lh + 832);
    const auto *lh_833 = buffer.data(lh + 833);
    const auto *lh_834 = buffer.data(lh + 834);
    const auto *lh_835 = buffer.data(lh + 835);
    const auto *lh_836 = buffer.data(lh + 836);
    const auto *lh_837 = buffer.data(lh + 837);
    const auto *lh_838 = buffer.data(lh + 838);
    const auto *lh_839 = buffer.data(lh + 839);
    const auto *lh_840 = buffer.data(lh + 840);
    const auto *lh_841 = buffer.data(lh + 841);
    const auto *lh_842 = buffer.data(lh + 842);
    const auto *lh_843 = buffer.data(lh + 843);
    const auto *lh_844 = buffer.data(lh + 844);
    const auto *lh_845 = buffer.data(lh + 845);
    const auto *lh_846 = buffer.data(lh + 846);
    const auto *lh_847 = buffer.data(lh + 847);
    const auto *lh_848 = buffer.data(lh + 848);
    const auto *lh_849 = buffer.data(lh + 849);
    const auto *lh_850 = buffer.data(lh + 850);
    const auto *lh_851 = buffer.data(lh + 851);
    const auto *lh_852 = buffer.data(lh + 852);
    const auto *lh_853 = buffer.data(lh + 853);
    const auto *lh_854 = buffer.data(lh + 854);
    const auto *lh_855 = buffer.data(lh + 855);
    const auto *lh_856 = buffer.data(lh + 856);
    const auto *lh_857 = buffer.data(lh + 857);
    const auto *lh_858 = buffer.data(lh + 858);
    const auto *lh_859 = buffer.data(lh + 859);
    const auto *lh_860 = buffer.data(lh + 860);
    const auto *lh_861 = buffer.data(lh + 861);
    const auto *lh_862 = buffer.data(lh + 862);
    const auto *lh_863 = buffer.data(lh + 863);
    const auto *lh_864 = buffer.data(lh + 864);
    const auto *lh_865 = buffer.data(lh + 865);
    const auto *lh_866 = buffer.data(lh + 866);
    const auto *lh_867 = buffer.data(lh + 867);
    const auto *lh_868 = buffer.data(lh + 868);
    const auto *lh_869 = buffer.data(lh + 869);
    const auto *lh_870 = buffer.data(lh + 870);
    const auto *lh_871 = buffer.data(lh + 871);
    const auto *lh_872 = buffer.data(lh + 872);
    const auto *lh_873 = buffer.data(lh + 873);
    const auto *lh_874 = buffer.data(lh + 874);
    const auto *lh_875 = buffer.data(lh + 875);
    const auto *lh_876 = buffer.data(lh + 876);
    const auto *lh_877 = buffer.data(lh + 877);
    const auto *lh_878 = buffer.data(lh + 878);
    const auto *lh_879 = buffer.data(lh + 879);
    const auto *lh_880 = buffer.data(lh + 880);
    const auto *lh_881 = buffer.data(lh + 881);
    const auto *lh_882 = buffer.data(lh + 882);
    const auto *lh_883 = buffer.data(lh + 883);
    const auto *lh_884 = buffer.data(lh + 884);
    const auto *lh_885 = buffer.data(lh + 885);
    const auto *lh_886 = buffer.data(lh + 886);
    const auto *lh_887 = buffer.data(lh + 887);
    const auto *lh_888 = buffer.data(lh + 888);
    const auto *lh_889 = buffer.data(lh + 889);
    const auto *lh_890 = buffer.data(lh + 890);
    const auto *lh_891 = buffer.data(lh + 891);
    const auto *lh_892 = buffer.data(lh + 892);
    const auto *lh_893 = buffer.data(lh + 893);
    const auto *lh_894 = buffer.data(lh + 894);
    const auto *lh_895 = buffer.data(lh + 895);
    const auto *lh_896 = buffer.data(lh + 896);
    const auto *lh_897 = buffer.data(lh + 897);
    const auto *lh_898 = buffer.data(lh + 898);
    const auto *lh_899 = buffer.data(lh + 899);
    const auto *lh_900 = buffer.data(lh + 900);
    const auto *lh_901 = buffer.data(lh + 901);
    const auto *lh_902 = buffer.data(lh + 902);
    const auto *lh_903 = buffer.data(lh + 903);
    const auto *lh_904 = buffer.data(lh + 904);
    const auto *lh_905 = buffer.data(lh + 905);
    const auto *lh_906 = buffer.data(lh + 906);
    const auto *lh_907 = buffer.data(lh + 907);
    const auto *lh_908 = buffer.data(lh + 908);
    const auto *lh_909 = buffer.data(lh + 909);
    const auto *lh_910 = buffer.data(lh + 910);
    const auto *lh_911 = buffer.data(lh + 911);
    const auto *lh_912 = buffer.data(lh + 912);
    const auto *lh_913 = buffer.data(lh + 913);
    const auto *lh_914 = buffer.data(lh + 914);
    const auto *lh_915 = buffer.data(lh + 915);
    const auto *lh_916 = buffer.data(lh + 916);
    const auto *lh_917 = buffer.data(lh + 917);
    const auto *lh_918 = buffer.data(lh + 918);
    const auto *lh_919 = buffer.data(lh + 919);
    const auto *lh_920 = buffer.data(lh + 920);
    const auto *lh_921 = buffer.data(lh + 921);
    const auto *lh_922 = buffer.data(lh + 922);
    const auto *lh_923 = buffer.data(lh + 923);

#pragma omp simd aligned(t_646, t_647, t_648, t_649, t_650, ih_499, ih_500, ih_501, ih_502, \
                         ih_503, lh_814, lh_815, lh_816, lh_817, \
                         lh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = -5.0 * ih_499[k]
                   + f_0 * lh_814[k];

        t_647[k] = -5.0 * ih_500[k]
                   + f_0 * lh_815[k];

        t_648[k] = -5.0 * ih_501[k]
                   + f_0 * lh_816[k];

        t_649[k] = -5.0 * ih_502[k]
                   + f_0 * lh_817[k];

        t_650[k] = -5.0 * ih_503[k]
                   + f_0 * lh_818[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, t_654, t_655, ih_504, ih_505, ih_506, ih_507, \
                         ih_508, lh_819, lh_820, lh_821, lh_822, \
                         lh_823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = -4.0 * ih_504[k]
                   + f_0 * lh_819[k];

        t_652[k] = -4.0 * ih_505[k]
                   + f_0 * lh_820[k];

        t_653[k] = -4.0 * ih_506[k]
                   + f_0 * lh_821[k];

        t_654[k] = -4.0 * ih_507[k]
                   + f_0 * lh_822[k];

        t_655[k] = -4.0 * ih_508[k]
                   + f_0 * lh_823[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, t_659, t_660, ih_509, ih_510, ih_511, ih_512, \
                         ih_513, lh_824, lh_825, lh_826, lh_827, \
                         lh_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = -4.0 * ih_509[k]
                   + f_0 * lh_824[k];

        t_657[k] = -4.0 * ih_510[k]
                   + f_0 * lh_825[k];

        t_658[k] = -4.0 * ih_511[k]
                   + f_0 * lh_826[k];

        t_659[k] = -4.0 * ih_512[k]
                   + f_0 * lh_827[k];

        t_660[k] = -4.0 * ih_513[k]
                   + f_0 * lh_828[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, ih_514, ih_515, ih_516, ih_517, \
                         ih_518, lh_829, lh_830, lh_831, lh_832, \
                         lh_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = -4.0 * ih_514[k]
                   + f_0 * lh_829[k];

        t_662[k] = -4.0 * ih_515[k]
                   + f_0 * lh_830[k];

        t_663[k] = -4.0 * ih_516[k]
                   + f_0 * lh_831[k];

        t_664[k] = -4.0 * ih_517[k]
                   + f_0 * lh_832[k];

        t_665[k] = -4.0 * ih_518[k]
                   + f_0 * lh_833[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, t_670, ih_519, ih_520, ih_521, ih_522, \
                         ih_523, lh_834, lh_835, lh_836, lh_837, \
                         lh_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = -4.0 * ih_519[k]
                   + f_0 * lh_834[k];

        t_667[k] = -4.0 * ih_520[k]
                   + f_0 * lh_835[k];

        t_668[k] = -4.0 * ih_521[k]
                   + f_0 * lh_836[k];

        t_669[k] = -4.0 * ih_522[k]
                   + f_0 * lh_837[k];

        t_670[k] = -4.0 * ih_523[k]
                   + f_0 * lh_838[k];
    }

#pragma omp simd aligned(t_671, t_672, t_673, t_674, t_675, ih_524, ih_525, ih_526, ih_527, \
                         ih_528, lh_839, lh_840, lh_841, lh_842, \
                         lh_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_671[k] = -4.0 * ih_524[k]
                   + f_0 * lh_839[k];

        t_672[k] = -3.0 * ih_525[k]
                   + f_0 * lh_840[k];

        t_673[k] = -3.0 * ih_526[k]
                   + f_0 * lh_841[k];

        t_674[k] = -3.0 * ih_527[k]
                   + f_0 * lh_842[k];

        t_675[k] = -3.0 * ih_528[k]
                   + f_0 * lh_843[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, t_680, ih_529, ih_530, ih_531, ih_532, \
                         ih_533, lh_844, lh_845, lh_846, lh_847, \
                         lh_848 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = -3.0 * ih_529[k]
                   + f_0 * lh_844[k];

        t_677[k] = -3.0 * ih_530[k]
                   + f_0 * lh_845[k];

        t_678[k] = -3.0 * ih_531[k]
                   + f_0 * lh_846[k];

        t_679[k] = -3.0 * ih_532[k]
                   + f_0 * lh_847[k];

        t_680[k] = -3.0 * ih_533[k]
                   + f_0 * lh_848[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, t_684, t_685, ih_534, ih_535, ih_536, ih_537, \
                         ih_538, lh_849, lh_850, lh_851, lh_852, \
                         lh_853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = -3.0 * ih_534[k]
                   + f_0 * lh_849[k];

        t_682[k] = -3.0 * ih_535[k]
                   + f_0 * lh_850[k];

        t_683[k] = -3.0 * ih_536[k]
                   + f_0 * lh_851[k];

        t_684[k] = -3.0 * ih_537[k]
                   + f_0 * lh_852[k];

        t_685[k] = -3.0 * ih_538[k]
                   + f_0 * lh_853[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, t_690, ih_539, ih_540, ih_541, ih_542, \
                         ih_543, lh_854, lh_855, lh_856, lh_857, \
                         lh_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = -3.0 * ih_539[k]
                   + f_0 * lh_854[k];

        t_687[k] = -3.0 * ih_540[k]
                   + f_0 * lh_855[k];

        t_688[k] = -3.0 * ih_541[k]
                   + f_0 * lh_856[k];

        t_689[k] = -3.0 * ih_542[k]
                   + f_0 * lh_857[k];

        t_690[k] = -3.0 * ih_543[k]
                   + f_0 * lh_858[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, t_695, ih_544, ih_545, ih_546, ih_547, \
                         ih_548, lh_859, lh_860, lh_861, lh_862, \
                         lh_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = -3.0 * ih_544[k]
                   + f_0 * lh_859[k];

        t_692[k] = -3.0 * ih_545[k]
                   + f_0 * lh_860[k];

        t_693[k] = -2.0 * ih_546[k]
                   + f_0 * lh_861[k];

        t_694[k] = -2.0 * ih_547[k]
                   + f_0 * lh_862[k];

        t_695[k] = -2.0 * ih_548[k]
                   + f_0 * lh_863[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, t_699, t_700, ih_549, ih_550, ih_551, ih_552, \
                         ih_553, lh_864, lh_865, lh_866, lh_867, \
                         lh_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = -2.0 * ih_549[k]
                   + f_0 * lh_864[k];

        t_697[k] = -2.0 * ih_550[k]
                   + f_0 * lh_865[k];

        t_698[k] = -2.0 * ih_551[k]
                   + f_0 * lh_866[k];

        t_699[k] = -2.0 * ih_552[k]
                   + f_0 * lh_867[k];

        t_700[k] = -2.0 * ih_553[k]
                   + f_0 * lh_868[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, t_705, ih_554, ih_555, ih_556, ih_557, \
                         ih_558, lh_869, lh_870, lh_871, lh_872, \
                         lh_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = -2.0 * ih_554[k]
                   + f_0 * lh_869[k];

        t_702[k] = -2.0 * ih_555[k]
                   + f_0 * lh_870[k];

        t_703[k] = -2.0 * ih_556[k]
                   + f_0 * lh_871[k];

        t_704[k] = -2.0 * ih_557[k]
                   + f_0 * lh_872[k];

        t_705[k] = -2.0 * ih_558[k]
                   + f_0 * lh_873[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, t_710, ih_559, ih_560, ih_561, ih_562, \
                         ih_563, lh_874, lh_875, lh_876, lh_877, \
                         lh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = -2.0 * ih_559[k]
                   + f_0 * lh_874[k];

        t_707[k] = -2.0 * ih_560[k]
                   + f_0 * lh_875[k];

        t_708[k] = -2.0 * ih_561[k]
                   + f_0 * lh_876[k];

        t_709[k] = -2.0 * ih_562[k]
                   + f_0 * lh_877[k];

        t_710[k] = -2.0 * ih_563[k]
                   + f_0 * lh_878[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, t_714, t_715, ih_564, ih_565, ih_566, ih_567, \
                         ih_568, lh_879, lh_880, lh_881, lh_882, \
                         lh_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = -2.0 * ih_564[k]
                   + f_0 * lh_879[k];

        t_712[k] = -2.0 * ih_565[k]
                   + f_0 * lh_880[k];

        t_713[k] = -2.0 * ih_566[k]
                   + f_0 * lh_881[k];

        t_714[k] = -ih_567[k]
                   + f_0 * lh_882[k];

        t_715[k] = -ih_568[k]
                   + f_0 * lh_883[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, t_720, ih_569, ih_570, ih_571, ih_572, \
                         ih_573, lh_884, lh_885, lh_886, lh_887, \
                         lh_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = -ih_569[k]
                   + f_0 * lh_884[k];

        t_717[k] = -ih_570[k]
                   + f_0 * lh_885[k];

        t_718[k] = -ih_571[k]
                   + f_0 * lh_886[k];

        t_719[k] = -ih_572[k]
                   + f_0 * lh_887[k];

        t_720[k] = -ih_573[k]
                   + f_0 * lh_888[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, t_725, ih_574, ih_575, ih_576, ih_577, \
                         ih_578, lh_889, lh_890, lh_891, lh_892, \
                         lh_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = -ih_574[k]
                   + f_0 * lh_889[k];

        t_722[k] = -ih_575[k]
                   + f_0 * lh_890[k];

        t_723[k] = -ih_576[k]
                   + f_0 * lh_891[k];

        t_724[k] = -ih_577[k]
                   + f_0 * lh_892[k];

        t_725[k] = -ih_578[k]
                   + f_0 * lh_893[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, t_730, ih_579, ih_580, ih_581, ih_582, \
                         ih_583, lh_894, lh_895, lh_896, lh_897, \
                         lh_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = -ih_579[k]
                   + f_0 * lh_894[k];

        t_727[k] = -ih_580[k]
                   + f_0 * lh_895[k];

        t_728[k] = -ih_581[k]
                   + f_0 * lh_896[k];

        t_729[k] = -ih_582[k]
                   + f_0 * lh_897[k];

        t_730[k] = -ih_583[k]
                   + f_0 * lh_898[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, t_734, t_735, t_736, ih_584, ih_585, ih_586, \
                         ih_587, lh_899, lh_900, lh_901, lh_902, lh_903, \
                         lh_904 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = -ih_584[k]
                   + f_0 * lh_899[k];

        t_732[k] = -ih_585[k]
                   + f_0 * lh_900[k];

        t_733[k] = -ih_586[k]
                   + f_0 * lh_901[k];

        t_734[k] = -ih_587[k]
                   + f_0 * lh_902[k];

        t_735[k] = f_0 * lh_903[k];

        t_736[k] = f_0 * lh_904[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, t_740, t_741, t_742, t_743, t_744, lh_905, \
                         lh_906, lh_907, lh_908, lh_909, lh_910, lh_911, \
                         lh_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_0 * lh_905[k];

        t_738[k] = f_0 * lh_906[k];

        t_739[k] = f_0 * lh_907[k];

        t_740[k] = f_0 * lh_908[k];

        t_741[k] = f_0 * lh_909[k];

        t_742[k] = f_0 * lh_910[k];

        t_743[k] = f_0 * lh_911[k];

        t_744[k] = f_0 * lh_912[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, t_750, t_751, t_752, lh_913, \
                         lh_914, lh_915, lh_916, lh_917, lh_918, lh_919, \
                         lh_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_0 * lh_913[k];

        t_746[k] = f_0 * lh_914[k];

        t_747[k] = f_0 * lh_915[k];

        t_748[k] = f_0 * lh_916[k];

        t_749[k] = f_0 * lh_917[k];

        t_750[k] = f_0 * lh_918[k];

        t_751[k] = f_0 * lh_919[k];

        t_752[k] = f_0 * lh_920[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, lh_921, lh_922, lh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_0 * lh_921[k];

        t_754[k] = f_0 * lh_922[k];

        t_755[k] = f_0 * lh_923[k];
    }
}

auto
compute_prim_geom_10_kh_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t ih, const size_t lh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kh_electron_repulsion_1_piece0(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_1_piece1(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_1_piece2(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_1_piece3(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_1_piece4(buffer, target, ih, lh, ncols, alpha);
}

static auto
compute_prim_geom_10_kh_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
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

    const auto *ih_0 = buffer.data(ih + 0);
    const auto *ih_1 = buffer.data(ih + 1);
    const auto *ih_2 = buffer.data(ih + 2);
    const auto *ih_3 = buffer.data(ih + 3);
    const auto *ih_4 = buffer.data(ih + 4);
    const auto *ih_5 = buffer.data(ih + 5);
    const auto *ih_6 = buffer.data(ih + 6);
    const auto *ih_7 = buffer.data(ih + 7);
    const auto *ih_8 = buffer.data(ih + 8);
    const auto *ih_9 = buffer.data(ih + 9);
    const auto *ih_10 = buffer.data(ih + 10);
    const auto *ih_11 = buffer.data(ih + 11);
    const auto *ih_12 = buffer.data(ih + 12);
    const auto *ih_13 = buffer.data(ih + 13);
    const auto *ih_14 = buffer.data(ih + 14);
    const auto *ih_15 = buffer.data(ih + 15);
    const auto *ih_16 = buffer.data(ih + 16);
    const auto *ih_17 = buffer.data(ih + 17);
    const auto *ih_18 = buffer.data(ih + 18);
    const auto *ih_19 = buffer.data(ih + 19);
    const auto *ih_20 = buffer.data(ih + 20);
    const auto *ih_21 = buffer.data(ih + 21);
    const auto *ih_22 = buffer.data(ih + 22);
    const auto *ih_23 = buffer.data(ih + 23);
    const auto *ih_24 = buffer.data(ih + 24);
    const auto *ih_25 = buffer.data(ih + 25);
    const auto *ih_26 = buffer.data(ih + 26);
    const auto *ih_27 = buffer.data(ih + 27);
    const auto *ih_28 = buffer.data(ih + 28);
    const auto *ih_29 = buffer.data(ih + 29);
    const auto *ih_30 = buffer.data(ih + 30);
    const auto *ih_31 = buffer.data(ih + 31);
    const auto *ih_32 = buffer.data(ih + 32);
    const auto *ih_33 = buffer.data(ih + 33);
    const auto *ih_34 = buffer.data(ih + 34);
    const auto *ih_35 = buffer.data(ih + 35);
    const auto *ih_36 = buffer.data(ih + 36);
    const auto *ih_37 = buffer.data(ih + 37);
    const auto *ih_38 = buffer.data(ih + 38);
    const auto *ih_39 = buffer.data(ih + 39);
    const auto *ih_40 = buffer.data(ih + 40);
    const auto *ih_41 = buffer.data(ih + 41);
    const auto *ih_42 = buffer.data(ih + 42);
    const auto *ih_43 = buffer.data(ih + 43);
    const auto *ih_44 = buffer.data(ih + 44);
    const auto *ih_45 = buffer.data(ih + 45);
    const auto *ih_46 = buffer.data(ih + 46);
    const auto *ih_47 = buffer.data(ih + 47);
    const auto *ih_48 = buffer.data(ih + 48);
    const auto *ih_49 = buffer.data(ih + 49);
    const auto *ih_50 = buffer.data(ih + 50);
    const auto *ih_51 = buffer.data(ih + 51);
    const auto *ih_52 = buffer.data(ih + 52);
    const auto *ih_53 = buffer.data(ih + 53);
    const auto *ih_54 = buffer.data(ih + 54);
    const auto *ih_55 = buffer.data(ih + 55);
    const auto *ih_56 = buffer.data(ih + 56);
    const auto *ih_57 = buffer.data(ih + 57);
    const auto *ih_58 = buffer.data(ih + 58);
    const auto *ih_59 = buffer.data(ih + 59);
    const auto *ih_60 = buffer.data(ih + 60);
    const auto *ih_61 = buffer.data(ih + 61);
    const auto *ih_62 = buffer.data(ih + 62);
    const auto *ih_63 = buffer.data(ih + 63);
    const auto *ih_64 = buffer.data(ih + 64);
    const auto *ih_65 = buffer.data(ih + 65);
    const auto *ih_66 = buffer.data(ih + 66);
    const auto *ih_67 = buffer.data(ih + 67);
    const auto *ih_68 = buffer.data(ih + 68);
    const auto *ih_69 = buffer.data(ih + 69);
    const auto *ih_70 = buffer.data(ih + 70);
    const auto *ih_71 = buffer.data(ih + 71);
    const auto *ih_72 = buffer.data(ih + 72);
    const auto *ih_73 = buffer.data(ih + 73);
    const auto *ih_74 = buffer.data(ih + 74);
    const auto *ih_75 = buffer.data(ih + 75);
    const auto *ih_76 = buffer.data(ih + 76);
    const auto *ih_77 = buffer.data(ih + 77);
    const auto *ih_78 = buffer.data(ih + 78);
    const auto *ih_79 = buffer.data(ih + 79);
    const auto *ih_80 = buffer.data(ih + 80);
    const auto *ih_81 = buffer.data(ih + 81);
    const auto *ih_82 = buffer.data(ih + 82);
    const auto *ih_83 = buffer.data(ih + 83);
    const auto *ih_84 = buffer.data(ih + 84);
    const auto *ih_85 = buffer.data(ih + 85);
    const auto *ih_86 = buffer.data(ih + 86);
    const auto *ih_87 = buffer.data(ih + 87);
    const auto *ih_88 = buffer.data(ih + 88);
    const auto *ih_89 = buffer.data(ih + 89);
    const auto *ih_90 = buffer.data(ih + 90);
    const auto *ih_91 = buffer.data(ih + 91);
    const auto *ih_92 = buffer.data(ih + 92);

    const auto *lh_42 = buffer.data(lh + 42);
    const auto *lh_43 = buffer.data(lh + 43);
    const auto *lh_44 = buffer.data(lh + 44);
    const auto *lh_45 = buffer.data(lh + 45);
    const auto *lh_46 = buffer.data(lh + 46);
    const auto *lh_47 = buffer.data(lh + 47);
    const auto *lh_48 = buffer.data(lh + 48);
    const auto *lh_49 = buffer.data(lh + 49);
    const auto *lh_50 = buffer.data(lh + 50);
    const auto *lh_51 = buffer.data(lh + 51);
    const auto *lh_52 = buffer.data(lh + 52);
    const auto *lh_53 = buffer.data(lh + 53);
    const auto *lh_54 = buffer.data(lh + 54);
    const auto *lh_55 = buffer.data(lh + 55);
    const auto *lh_56 = buffer.data(lh + 56);
    const auto *lh_57 = buffer.data(lh + 57);
    const auto *lh_58 = buffer.data(lh + 58);
    const auto *lh_59 = buffer.data(lh + 59);
    const auto *lh_60 = buffer.data(lh + 60);
    const auto *lh_61 = buffer.data(lh + 61);
    const auto *lh_62 = buffer.data(lh + 62);
    const auto *lh_84 = buffer.data(lh + 84);
    const auto *lh_85 = buffer.data(lh + 85);
    const auto *lh_86 = buffer.data(lh + 86);
    const auto *lh_87 = buffer.data(lh + 87);
    const auto *lh_88 = buffer.data(lh + 88);
    const auto *lh_89 = buffer.data(lh + 89);
    const auto *lh_90 = buffer.data(lh + 90);
    const auto *lh_91 = buffer.data(lh + 91);
    const auto *lh_92 = buffer.data(lh + 92);
    const auto *lh_93 = buffer.data(lh + 93);
    const auto *lh_94 = buffer.data(lh + 94);
    const auto *lh_95 = buffer.data(lh + 95);
    const auto *lh_96 = buffer.data(lh + 96);
    const auto *lh_97 = buffer.data(lh + 97);
    const auto *lh_98 = buffer.data(lh + 98);
    const auto *lh_99 = buffer.data(lh + 99);
    const auto *lh_100 = buffer.data(lh + 100);
    const auto *lh_101 = buffer.data(lh + 101);
    const auto *lh_102 = buffer.data(lh + 102);
    const auto *lh_103 = buffer.data(lh + 103);
    const auto *lh_104 = buffer.data(lh + 104);
    const auto *lh_105 = buffer.data(lh + 105);
    const auto *lh_106 = buffer.data(lh + 106);
    const auto *lh_107 = buffer.data(lh + 107);
    const auto *lh_108 = buffer.data(lh + 108);
    const auto *lh_109 = buffer.data(lh + 109);
    const auto *lh_110 = buffer.data(lh + 110);
    const auto *lh_111 = buffer.data(lh + 111);
    const auto *lh_112 = buffer.data(lh + 112);
    const auto *lh_113 = buffer.data(lh + 113);
    const auto *lh_114 = buffer.data(lh + 114);
    const auto *lh_115 = buffer.data(lh + 115);
    const auto *lh_116 = buffer.data(lh + 116);
    const auto *lh_117 = buffer.data(lh + 117);
    const auto *lh_118 = buffer.data(lh + 118);
    const auto *lh_119 = buffer.data(lh + 119);
    const auto *lh_120 = buffer.data(lh + 120);
    const auto *lh_121 = buffer.data(lh + 121);
    const auto *lh_122 = buffer.data(lh + 122);
    const auto *lh_123 = buffer.data(lh + 123);
    const auto *lh_124 = buffer.data(lh + 124);
    const auto *lh_125 = buffer.data(lh + 125);
    const auto *lh_147 = buffer.data(lh + 147);
    const auto *lh_148 = buffer.data(lh + 148);
    const auto *lh_149 = buffer.data(lh + 149);
    const auto *lh_150 = buffer.data(lh + 150);
    const auto *lh_151 = buffer.data(lh + 151);
    const auto *lh_152 = buffer.data(lh + 152);
    const auto *lh_153 = buffer.data(lh + 153);
    const auto *lh_154 = buffer.data(lh + 154);
    const auto *lh_155 = buffer.data(lh + 155);
    const auto *lh_156 = buffer.data(lh + 156);
    const auto *lh_157 = buffer.data(lh + 157);
    const auto *lh_158 = buffer.data(lh + 158);
    const auto *lh_159 = buffer.data(lh + 159);
    const auto *lh_160 = buffer.data(lh + 160);
    const auto *lh_161 = buffer.data(lh + 161);
    const auto *lh_162 = buffer.data(lh + 162);
    const auto *lh_163 = buffer.data(lh + 163);
    const auto *lh_164 = buffer.data(lh + 164);
    const auto *lh_165 = buffer.data(lh + 165);
    const auto *lh_166 = buffer.data(lh + 166);
    const auto *lh_167 = buffer.data(lh + 167);
    const auto *lh_168 = buffer.data(lh + 168);
    const auto *lh_169 = buffer.data(lh + 169);
    const auto *lh_170 = buffer.data(lh + 170);
    const auto *lh_171 = buffer.data(lh + 171);
    const auto *lh_172 = buffer.data(lh + 172);
    const auto *lh_173 = buffer.data(lh + 173);
    const auto *lh_174 = buffer.data(lh + 174);
    const auto *lh_175 = buffer.data(lh + 175);
    const auto *lh_176 = buffer.data(lh + 176);
    const auto *lh_177 = buffer.data(lh + 177);
    const auto *lh_178 = buffer.data(lh + 178);
    const auto *lh_179 = buffer.data(lh + 179);
    const auto *lh_180 = buffer.data(lh + 180);
    const auto *lh_181 = buffer.data(lh + 181);
    const auto *lh_182 = buffer.data(lh + 182);
    const auto *lh_183 = buffer.data(lh + 183);
    const auto *lh_184 = buffer.data(lh + 184);
    const auto *lh_185 = buffer.data(lh + 185);
    const auto *lh_186 = buffer.data(lh + 186);
    const auto *lh_187 = buffer.data(lh + 187);
    const auto *lh_188 = buffer.data(lh + 188);
    const auto *lh_189 = buffer.data(lh + 189);
    const auto *lh_190 = buffer.data(lh + 190);
    const auto *lh_191 = buffer.data(lh + 191);
    const auto *lh_192 = buffer.data(lh + 192);
    const auto *lh_193 = buffer.data(lh + 193);
    const auto *lh_194 = buffer.data(lh + 194);
    const auto *lh_195 = buffer.data(lh + 195);
    const auto *lh_196 = buffer.data(lh + 196);
    const auto *lh_197 = buffer.data(lh + 197);
    const auto *lh_198 = buffer.data(lh + 198);
    const auto *lh_199 = buffer.data(lh + 199);
    const auto *lh_200 = buffer.data(lh + 200);
    const auto *lh_201 = buffer.data(lh + 201);
    const auto *lh_202 = buffer.data(lh + 202);
    const auto *lh_203 = buffer.data(lh + 203);
    const auto *lh_204 = buffer.data(lh + 204);
    const auto *lh_205 = buffer.data(lh + 205);
    const auto *lh_206 = buffer.data(lh + 206);
    const auto *lh_207 = buffer.data(lh + 207);
    const auto *lh_208 = buffer.data(lh + 208);
    const auto *lh_209 = buffer.data(lh + 209);
    const auto *lh_231 = buffer.data(lh + 231);
    const auto *lh_232 = buffer.data(lh + 232);
    const auto *lh_233 = buffer.data(lh + 233);
    const auto *lh_234 = buffer.data(lh + 234);
    const auto *lh_235 = buffer.data(lh + 235);
    const auto *lh_236 = buffer.data(lh + 236);
    const auto *lh_237 = buffer.data(lh + 237);
    const auto *lh_238 = buffer.data(lh + 238);
    const auto *lh_239 = buffer.data(lh + 239);
    const auto *lh_240 = buffer.data(lh + 240);
    const auto *lh_241 = buffer.data(lh + 241);
    const auto *lh_242 = buffer.data(lh + 242);
    const auto *lh_243 = buffer.data(lh + 243);
    const auto *lh_244 = buffer.data(lh + 244);
    const auto *lh_245 = buffer.data(lh + 245);
    const auto *lh_246 = buffer.data(lh + 246);
    const auto *lh_247 = buffer.data(lh + 247);
    const auto *lh_248 = buffer.data(lh + 248);
    const auto *lh_249 = buffer.data(lh + 249);
    const auto *lh_250 = buffer.data(lh + 250);
    const auto *lh_251 = buffer.data(lh + 251);
    const auto *lh_252 = buffer.data(lh + 252);
    const auto *lh_253 = buffer.data(lh + 253);
    const auto *lh_254 = buffer.data(lh + 254);
    const auto *lh_255 = buffer.data(lh + 255);
    const auto *lh_256 = buffer.data(lh + 256);
    const auto *lh_257 = buffer.data(lh + 257);
    const auto *lh_258 = buffer.data(lh + 258);
    const auto *lh_259 = buffer.data(lh + 259);
    const auto *lh_260 = buffer.data(lh + 260);
    const auto *lh_261 = buffer.data(lh + 261);
    const auto *lh_262 = buffer.data(lh + 262);
    const auto *lh_263 = buffer.data(lh + 263);
    const auto *lh_264 = buffer.data(lh + 264);
    const auto *lh_265 = buffer.data(lh + 265);
    const auto *lh_266 = buffer.data(lh + 266);
    const auto *lh_267 = buffer.data(lh + 267);
    const auto *lh_268 = buffer.data(lh + 268);
    const auto *lh_269 = buffer.data(lh + 269);
    const auto *lh_270 = buffer.data(lh + 270);
    const auto *lh_271 = buffer.data(lh + 271);
    const auto *lh_272 = buffer.data(lh + 272);
    const auto *lh_273 = buffer.data(lh + 273);
    const auto *lh_274 = buffer.data(lh + 274);
    const auto *lh_275 = buffer.data(lh + 275);
    const auto *lh_276 = buffer.data(lh + 276);
    const auto *lh_277 = buffer.data(lh + 277);
    const auto *lh_278 = buffer.data(lh + 278);
    const auto *lh_279 = buffer.data(lh + 279);
    const auto *lh_280 = buffer.data(lh + 280);
    const auto *lh_281 = buffer.data(lh + 281);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, lh_42, lh_43, lh_44, lh_45, \
                         lh_46, lh_47, lh_48, lh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lh_42[k];

        t_1[k] = f_0 * lh_43[k];

        t_2[k] = f_0 * lh_44[k];

        t_3[k] = f_0 * lh_45[k];

        t_4[k] = f_0 * lh_46[k];

        t_5[k] = f_0 * lh_47[k];

        t_6[k] = f_0 * lh_48[k];

        t_7[k] = f_0 * lh_49[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, lh_50, lh_51, lh_52, \
                         lh_53, lh_54, lh_55, lh_56, lh_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * lh_50[k];

        t_9[k] = f_0 * lh_51[k];

        t_10[k] = f_0 * lh_52[k];

        t_11[k] = f_0 * lh_53[k];

        t_12[k] = f_0 * lh_54[k];

        t_13[k] = f_0 * lh_55[k];

        t_14[k] = f_0 * lh_56[k];

        t_15[k] = f_0 * lh_57[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, lh_58, lh_59, lh_60, \
                         lh_61, lh_62, lh_84, lh_85, lh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * lh_58[k];

        t_17[k] = f_0 * lh_59[k];

        t_18[k] = f_0 * lh_60[k];

        t_19[k] = f_0 * lh_61[k];

        t_20[k] = f_0 * lh_62[k];

        t_21[k] = f_0 * lh_84[k];

        t_22[k] = f_0 * lh_85[k];

        t_23[k] = f_0 * lh_86[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, lh_87, lh_88, lh_89, \
                         lh_90, lh_91, lh_92, lh_93, lh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * lh_87[k];

        t_25[k] = f_0 * lh_88[k];

        t_26[k] = f_0 * lh_89[k];

        t_27[k] = f_0 * lh_90[k];

        t_28[k] = f_0 * lh_91[k];

        t_29[k] = f_0 * lh_92[k];

        t_30[k] = f_0 * lh_93[k];

        t_31[k] = f_0 * lh_94[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, lh_95, lh_96, lh_97, \
                         lh_98, lh_99, lh_100, lh_101, lh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * lh_95[k];

        t_33[k] = f_0 * lh_96[k];

        t_34[k] = f_0 * lh_97[k];

        t_35[k] = f_0 * lh_98[k];

        t_36[k] = f_0 * lh_99[k];

        t_37[k] = f_0 * lh_100[k];

        t_38[k] = f_0 * lh_101[k];

        t_39[k] = f_0 * lh_102[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, ih_0, ih_1, ih_2, ih_3, lh_103, \
                         lh_104, lh_105, lh_106, lh_107, lh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * lh_103[k];

        t_41[k] = f_0 * lh_104[k];

        t_42[k] = -ih_0[k]
                  + f_0 * lh_105[k];

        t_43[k] = -ih_1[k]
                  + f_0 * lh_106[k];

        t_44[k] = -ih_2[k]
                  + f_0 * lh_107[k];

        t_45[k] = -ih_3[k]
                  + f_0 * lh_108[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, ih_4, ih_5, ih_6, ih_7, ih_8, lh_109, \
                         lh_110, lh_111, lh_112, lh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -ih_4[k]
                  + f_0 * lh_109[k];

        t_47[k] = -ih_5[k]
                  + f_0 * lh_110[k];

        t_48[k] = -ih_6[k]
                  + f_0 * lh_111[k];

        t_49[k] = -ih_7[k]
                  + f_0 * lh_112[k];

        t_50[k] = -ih_8[k]
                  + f_0 * lh_113[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, ih_9, ih_10, ih_11, ih_12, ih_13, \
                         lh_114, lh_115, lh_116, lh_117, lh_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -ih_9[k]
                  + f_0 * lh_114[k];

        t_52[k] = -ih_10[k]
                  + f_0 * lh_115[k];

        t_53[k] = -ih_11[k]
                  + f_0 * lh_116[k];

        t_54[k] = -ih_12[k]
                  + f_0 * lh_117[k];

        t_55[k] = -ih_13[k]
                  + f_0 * lh_118[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, ih_14, ih_15, ih_16, ih_17, ih_18, \
                         lh_119, lh_120, lh_121, lh_122, lh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -ih_14[k]
                  + f_0 * lh_119[k];

        t_57[k] = -ih_15[k]
                  + f_0 * lh_120[k];

        t_58[k] = -ih_16[k]
                  + f_0 * lh_121[k];

        t_59[k] = -ih_17[k]
                  + f_0 * lh_122[k];

        t_60[k] = -ih_18[k]
                  + f_0 * lh_123[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, t_67, ih_19, ih_20, lh_124, \
                         lh_125, lh_147, lh_148, lh_149, lh_150, \
                         lh_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -ih_19[k]
                  + f_0 * lh_124[k];

        t_62[k] = -ih_20[k]
                  + f_0 * lh_125[k];

        t_63[k] = f_0 * lh_147[k];

        t_64[k] = f_0 * lh_148[k];

        t_65[k] = f_0 * lh_149[k];

        t_66[k] = f_0 * lh_150[k];

        t_67[k] = f_0 * lh_151[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, t_74, t_75, lh_152, lh_153, \
                         lh_154, lh_155, lh_156, lh_157, lh_158, \
                         lh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * lh_152[k];

        t_69[k] = f_0 * lh_153[k];

        t_70[k] = f_0 * lh_154[k];

        t_71[k] = f_0 * lh_155[k];

        t_72[k] = f_0 * lh_156[k];

        t_73[k] = f_0 * lh_157[k];

        t_74[k] = f_0 * lh_158[k];

        t_75[k] = f_0 * lh_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, t_82, t_83, lh_160, lh_161, \
                         lh_162, lh_163, lh_164, lh_165, lh_166, \
                         lh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_0 * lh_160[k];

        t_77[k] = f_0 * lh_161[k];

        t_78[k] = f_0 * lh_162[k];

        t_79[k] = f_0 * lh_163[k];

        t_80[k] = f_0 * lh_164[k];

        t_81[k] = f_0 * lh_165[k];

        t_82[k] = f_0 * lh_166[k];

        t_83[k] = f_0 * lh_167[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ih_21, ih_22, ih_23, ih_24, ih_25, \
                         lh_168, lh_169, lh_170, lh_171, lh_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -ih_21[k]
                  + f_0 * lh_168[k];

        t_85[k] = -ih_22[k]
                  + f_0 * lh_169[k];

        t_86[k] = -ih_23[k]
                  + f_0 * lh_170[k];

        t_87[k] = -ih_24[k]
                  + f_0 * lh_171[k];

        t_88[k] = -ih_25[k]
                  + f_0 * lh_172[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ih_26, ih_27, ih_28, ih_29, ih_30, \
                         lh_173, lh_174, lh_175, lh_176, lh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -ih_26[k]
                  + f_0 * lh_173[k];

        t_90[k] = -ih_27[k]
                  + f_0 * lh_174[k];

        t_91[k] = -ih_28[k]
                  + f_0 * lh_175[k];

        t_92[k] = -ih_29[k]
                  + f_0 * lh_176[k];

        t_93[k] = -ih_30[k]
                  + f_0 * lh_177[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ih_31, ih_32, ih_33, ih_34, ih_35, \
                         lh_178, lh_179, lh_180, lh_181, lh_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -ih_31[k]
                  + f_0 * lh_178[k];

        t_95[k] = -ih_32[k]
                  + f_0 * lh_179[k];

        t_96[k] = -ih_33[k]
                  + f_0 * lh_180[k];

        t_97[k] = -ih_34[k]
                  + f_0 * lh_181[k];

        t_98[k] = -ih_35[k]
                  + f_0 * lh_182[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ih_36, ih_37, ih_38, ih_39, ih_40, \
                         lh_183, lh_184, lh_185, lh_186, lh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -ih_36[k]
                  + f_0 * lh_183[k];

        t_100[k] = -ih_37[k]
                   + f_0 * lh_184[k];

        t_101[k] = -ih_38[k]
                   + f_0 * lh_185[k];

        t_102[k] = -ih_39[k]
                   + f_0 * lh_186[k];

        t_103[k] = -ih_40[k]
                   + f_0 * lh_187[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, ih_41, ih_42, ih_43, ih_44, ih_45, \
                         lh_188, lh_189, lh_190, lh_191, lh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -ih_41[k]
                   + f_0 * lh_188[k];

        t_105[k] = -2.0 * ih_42[k]
                   + f_0 * lh_189[k];

        t_106[k] = -2.0 * ih_43[k]
                   + f_0 * lh_190[k];

        t_107[k] = -2.0 * ih_44[k]
                   + f_0 * lh_191[k];

        t_108[k] = -2.0 * ih_45[k]
                   + f_0 * lh_192[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, ih_46, ih_47, ih_48, ih_49, ih_50, \
                         lh_193, lh_194, lh_195, lh_196, lh_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -2.0 * ih_46[k]
                   + f_0 * lh_193[k];

        t_110[k] = -2.0 * ih_47[k]
                   + f_0 * lh_194[k];

        t_111[k] = -2.0 * ih_48[k]
                   + f_0 * lh_195[k];

        t_112[k] = -2.0 * ih_49[k]
                   + f_0 * lh_196[k];

        t_113[k] = -2.0 * ih_50[k]
                   + f_0 * lh_197[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, ih_51, ih_52, ih_53, ih_54, ih_55, \
                         lh_198, lh_199, lh_200, lh_201, lh_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -2.0 * ih_51[k]
                   + f_0 * lh_198[k];

        t_115[k] = -2.0 * ih_52[k]
                   + f_0 * lh_199[k];

        t_116[k] = -2.0 * ih_53[k]
                   + f_0 * lh_200[k];

        t_117[k] = -2.0 * ih_54[k]
                   + f_0 * lh_201[k];

        t_118[k] = -2.0 * ih_55[k]
                   + f_0 * lh_202[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, ih_56, ih_57, ih_58, ih_59, ih_60, \
                         lh_203, lh_204, lh_205, lh_206, lh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -2.0 * ih_56[k]
                   + f_0 * lh_203[k];

        t_120[k] = -2.0 * ih_57[k]
                   + f_0 * lh_204[k];

        t_121[k] = -2.0 * ih_58[k]
                   + f_0 * lh_205[k];

        t_122[k] = -2.0 * ih_59[k]
                   + f_0 * lh_206[k];

        t_123[k] = -2.0 * ih_60[k]
                   + f_0 * lh_207[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, t_129, t_130, ih_61, ih_62, \
                         lh_208, lh_209, lh_231, lh_232, lh_233, lh_234, \
                         lh_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * ih_61[k]
                   + f_0 * lh_208[k];

        t_125[k] = -2.0 * ih_62[k]
                   + f_0 * lh_209[k];

        t_126[k] = f_0 * lh_231[k];

        t_127[k] = f_0 * lh_232[k];

        t_128[k] = f_0 * lh_233[k];

        t_129[k] = f_0 * lh_234[k];

        t_130[k] = f_0 * lh_235[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, t_136, t_137, t_138, lh_236, \
                         lh_237, lh_238, lh_239, lh_240, lh_241, lh_242, \
                         lh_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_0 * lh_236[k];

        t_132[k] = f_0 * lh_237[k];

        t_133[k] = f_0 * lh_238[k];

        t_134[k] = f_0 * lh_239[k];

        t_135[k] = f_0 * lh_240[k];

        t_136[k] = f_0 * lh_241[k];

        t_137[k] = f_0 * lh_242[k];

        t_138[k] = f_0 * lh_243[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, t_145, t_146, lh_244, \
                         lh_245, lh_246, lh_247, lh_248, lh_249, lh_250, \
                         lh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_0 * lh_244[k];

        t_140[k] = f_0 * lh_245[k];

        t_141[k] = f_0 * lh_246[k];

        t_142[k] = f_0 * lh_247[k];

        t_143[k] = f_0 * lh_248[k];

        t_144[k] = f_0 * lh_249[k];

        t_145[k] = f_0 * lh_250[k];

        t_146[k] = f_0 * lh_251[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, ih_63, ih_64, ih_65, ih_66, ih_67, \
                         lh_252, lh_253, lh_254, lh_255, lh_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -ih_63[k]
                   + f_0 * lh_252[k];

        t_148[k] = -ih_64[k]
                   + f_0 * lh_253[k];

        t_149[k] = -ih_65[k]
                   + f_0 * lh_254[k];

        t_150[k] = -ih_66[k]
                   + f_0 * lh_255[k];

        t_151[k] = -ih_67[k]
                   + f_0 * lh_256[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, ih_68, ih_69, ih_70, ih_71, ih_72, \
                         lh_257, lh_258, lh_259, lh_260, lh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -ih_68[k]
                   + f_0 * lh_257[k];

        t_153[k] = -ih_69[k]
                   + f_0 * lh_258[k];

        t_154[k] = -ih_70[k]
                   + f_0 * lh_259[k];

        t_155[k] = -ih_71[k]
                   + f_0 * lh_260[k];

        t_156[k] = -ih_72[k]
                   + f_0 * lh_261[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, ih_73, ih_74, ih_75, ih_76, ih_77, \
                         lh_262, lh_263, lh_264, lh_265, lh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -ih_73[k]
                   + f_0 * lh_262[k];

        t_158[k] = -ih_74[k]
                   + f_0 * lh_263[k];

        t_159[k] = -ih_75[k]
                   + f_0 * lh_264[k];

        t_160[k] = -ih_76[k]
                   + f_0 * lh_265[k];

        t_161[k] = -ih_77[k]
                   + f_0 * lh_266[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, ih_78, ih_79, ih_80, ih_81, ih_82, \
                         lh_267, lh_268, lh_269, lh_270, lh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -ih_78[k]
                   + f_0 * lh_267[k];

        t_163[k] = -ih_79[k]
                   + f_0 * lh_268[k];

        t_164[k] = -ih_80[k]
                   + f_0 * lh_269[k];

        t_165[k] = -ih_81[k]
                   + f_0 * lh_270[k];

        t_166[k] = -ih_82[k]
                   + f_0 * lh_271[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, ih_83, ih_84, ih_85, ih_86, ih_87, \
                         lh_272, lh_273, lh_274, lh_275, lh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -ih_83[k]
                   + f_0 * lh_272[k];

        t_168[k] = -2.0 * ih_84[k]
                   + f_0 * lh_273[k];

        t_169[k] = -2.0 * ih_85[k]
                   + f_0 * lh_274[k];

        t_170[k] = -2.0 * ih_86[k]
                   + f_0 * lh_275[k];

        t_171[k] = -2.0 * ih_87[k]
                   + f_0 * lh_276[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, ih_88, ih_89, ih_90, ih_91, ih_92, \
                         lh_277, lh_278, lh_279, lh_280, lh_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -2.0 * ih_88[k]
                   + f_0 * lh_277[k];

        t_173[k] = -2.0 * ih_89[k]
                   + f_0 * lh_278[k];

        t_174[k] = -2.0 * ih_90[k]
                   + f_0 * lh_279[k];

        t_175[k] = -2.0 * ih_91[k]
                   + f_0 * lh_280[k];

        t_176[k] = -2.0 * ih_92[k]
                   + f_0 * lh_281[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
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

    const auto *ih_93 = buffer.data(ih + 93);
    const auto *ih_94 = buffer.data(ih + 94);
    const auto *ih_95 = buffer.data(ih + 95);
    const auto *ih_96 = buffer.data(ih + 96);
    const auto *ih_97 = buffer.data(ih + 97);
    const auto *ih_98 = buffer.data(ih + 98);
    const auto *ih_99 = buffer.data(ih + 99);
    const auto *ih_100 = buffer.data(ih + 100);
    const auto *ih_101 = buffer.data(ih + 101);
    const auto *ih_102 = buffer.data(ih + 102);
    const auto *ih_103 = buffer.data(ih + 103);
    const auto *ih_104 = buffer.data(ih + 104);
    const auto *ih_105 = buffer.data(ih + 105);
    const auto *ih_106 = buffer.data(ih + 106);
    const auto *ih_107 = buffer.data(ih + 107);
    const auto *ih_108 = buffer.data(ih + 108);
    const auto *ih_109 = buffer.data(ih + 109);
    const auto *ih_110 = buffer.data(ih + 110);
    const auto *ih_111 = buffer.data(ih + 111);
    const auto *ih_112 = buffer.data(ih + 112);
    const auto *ih_113 = buffer.data(ih + 113);
    const auto *ih_114 = buffer.data(ih + 114);
    const auto *ih_115 = buffer.data(ih + 115);
    const auto *ih_116 = buffer.data(ih + 116);
    const auto *ih_117 = buffer.data(ih + 117);
    const auto *ih_118 = buffer.data(ih + 118);
    const auto *ih_119 = buffer.data(ih + 119);
    const auto *ih_120 = buffer.data(ih + 120);
    const auto *ih_121 = buffer.data(ih + 121);
    const auto *ih_122 = buffer.data(ih + 122);
    const auto *ih_123 = buffer.data(ih + 123);
    const auto *ih_124 = buffer.data(ih + 124);
    const auto *ih_125 = buffer.data(ih + 125);
    const auto *ih_126 = buffer.data(ih + 126);
    const auto *ih_127 = buffer.data(ih + 127);
    const auto *ih_128 = buffer.data(ih + 128);
    const auto *ih_129 = buffer.data(ih + 129);
    const auto *ih_130 = buffer.data(ih + 130);
    const auto *ih_131 = buffer.data(ih + 131);
    const auto *ih_132 = buffer.data(ih + 132);
    const auto *ih_133 = buffer.data(ih + 133);
    const auto *ih_134 = buffer.data(ih + 134);
    const auto *ih_135 = buffer.data(ih + 135);
    const auto *ih_136 = buffer.data(ih + 136);
    const auto *ih_137 = buffer.data(ih + 137);
    const auto *ih_138 = buffer.data(ih + 138);
    const auto *ih_139 = buffer.data(ih + 139);
    const auto *ih_140 = buffer.data(ih + 140);
    const auto *ih_141 = buffer.data(ih + 141);
    const auto *ih_142 = buffer.data(ih + 142);
    const auto *ih_143 = buffer.data(ih + 143);
    const auto *ih_144 = buffer.data(ih + 144);
    const auto *ih_145 = buffer.data(ih + 145);
    const auto *ih_146 = buffer.data(ih + 146);
    const auto *ih_147 = buffer.data(ih + 147);
    const auto *ih_148 = buffer.data(ih + 148);
    const auto *ih_149 = buffer.data(ih + 149);
    const auto *ih_150 = buffer.data(ih + 150);
    const auto *ih_151 = buffer.data(ih + 151);
    const auto *ih_152 = buffer.data(ih + 152);
    const auto *ih_153 = buffer.data(ih + 153);
    const auto *ih_154 = buffer.data(ih + 154);
    const auto *ih_155 = buffer.data(ih + 155);
    const auto *ih_156 = buffer.data(ih + 156);
    const auto *ih_157 = buffer.data(ih + 157);
    const auto *ih_158 = buffer.data(ih + 158);
    const auto *ih_159 = buffer.data(ih + 159);
    const auto *ih_160 = buffer.data(ih + 160);
    const auto *ih_161 = buffer.data(ih + 161);
    const auto *ih_162 = buffer.data(ih + 162);
    const auto *ih_163 = buffer.data(ih + 163);
    const auto *ih_164 = buffer.data(ih + 164);
    const auto *ih_165 = buffer.data(ih + 165);
    const auto *ih_166 = buffer.data(ih + 166);
    const auto *ih_167 = buffer.data(ih + 167);
    const auto *ih_168 = buffer.data(ih + 168);
    const auto *ih_169 = buffer.data(ih + 169);
    const auto *ih_170 = buffer.data(ih + 170);
    const auto *ih_171 = buffer.data(ih + 171);
    const auto *ih_172 = buffer.data(ih + 172);
    const auto *ih_173 = buffer.data(ih + 173);
    const auto *ih_174 = buffer.data(ih + 174);
    const auto *ih_175 = buffer.data(ih + 175);
    const auto *ih_176 = buffer.data(ih + 176);
    const auto *ih_177 = buffer.data(ih + 177);
    const auto *ih_178 = buffer.data(ih + 178);
    const auto *ih_179 = buffer.data(ih + 179);
    const auto *ih_180 = buffer.data(ih + 180);
    const auto *ih_181 = buffer.data(ih + 181);
    const auto *ih_182 = buffer.data(ih + 182);
    const auto *ih_183 = buffer.data(ih + 183);
    const auto *ih_184 = buffer.data(ih + 184);
    const auto *ih_185 = buffer.data(ih + 185);
    const auto *ih_186 = buffer.data(ih + 186);
    const auto *ih_187 = buffer.data(ih + 187);
    const auto *ih_188 = buffer.data(ih + 188);
    const auto *ih_189 = buffer.data(ih + 189);
    const auto *ih_190 = buffer.data(ih + 190);
    const auto *ih_191 = buffer.data(ih + 191);
    const auto *ih_192 = buffer.data(ih + 192);
    const auto *ih_193 = buffer.data(ih + 193);
    const auto *ih_194 = buffer.data(ih + 194);
    const auto *ih_195 = buffer.data(ih + 195);
    const auto *ih_196 = buffer.data(ih + 196);
    const auto *ih_197 = buffer.data(ih + 197);
    const auto *ih_198 = buffer.data(ih + 198);
    const auto *ih_199 = buffer.data(ih + 199);
    const auto *ih_200 = buffer.data(ih + 200);
    const auto *ih_201 = buffer.data(ih + 201);
    const auto *ih_202 = buffer.data(ih + 202);
    const auto *ih_203 = buffer.data(ih + 203);
    const auto *ih_204 = buffer.data(ih + 204);
    const auto *ih_205 = buffer.data(ih + 205);
    const auto *ih_206 = buffer.data(ih + 206);
    const auto *ih_207 = buffer.data(ih + 207);
    const auto *ih_208 = buffer.data(ih + 208);
    const auto *ih_209 = buffer.data(ih + 209);
    const auto *ih_210 = buffer.data(ih + 210);
    const auto *ih_211 = buffer.data(ih + 211);

    const auto *lh_282 = buffer.data(lh + 282);
    const auto *lh_283 = buffer.data(lh + 283);
    const auto *lh_284 = buffer.data(lh + 284);
    const auto *lh_285 = buffer.data(lh + 285);
    const auto *lh_286 = buffer.data(lh + 286);
    const auto *lh_287 = buffer.data(lh + 287);
    const auto *lh_288 = buffer.data(lh + 288);
    const auto *lh_289 = buffer.data(lh + 289);
    const auto *lh_290 = buffer.data(lh + 290);
    const auto *lh_291 = buffer.data(lh + 291);
    const auto *lh_292 = buffer.data(lh + 292);
    const auto *lh_293 = buffer.data(lh + 293);
    const auto *lh_294 = buffer.data(lh + 294);
    const auto *lh_295 = buffer.data(lh + 295);
    const auto *lh_296 = buffer.data(lh + 296);
    const auto *lh_297 = buffer.data(lh + 297);
    const auto *lh_298 = buffer.data(lh + 298);
    const auto *lh_299 = buffer.data(lh + 299);
    const auto *lh_300 = buffer.data(lh + 300);
    const auto *lh_301 = buffer.data(lh + 301);
    const auto *lh_302 = buffer.data(lh + 302);
    const auto *lh_303 = buffer.data(lh + 303);
    const auto *lh_304 = buffer.data(lh + 304);
    const auto *lh_305 = buffer.data(lh + 305);
    const auto *lh_306 = buffer.data(lh + 306);
    const auto *lh_307 = buffer.data(lh + 307);
    const auto *lh_308 = buffer.data(lh + 308);
    const auto *lh_309 = buffer.data(lh + 309);
    const auto *lh_310 = buffer.data(lh + 310);
    const auto *lh_311 = buffer.data(lh + 311);
    const auto *lh_312 = buffer.data(lh + 312);
    const auto *lh_313 = buffer.data(lh + 313);
    const auto *lh_314 = buffer.data(lh + 314);
    const auto *lh_336 = buffer.data(lh + 336);
    const auto *lh_337 = buffer.data(lh + 337);
    const auto *lh_338 = buffer.data(lh + 338);
    const auto *lh_339 = buffer.data(lh + 339);
    const auto *lh_340 = buffer.data(lh + 340);
    const auto *lh_341 = buffer.data(lh + 341);
    const auto *lh_342 = buffer.data(lh + 342);
    const auto *lh_343 = buffer.data(lh + 343);
    const auto *lh_344 = buffer.data(lh + 344);
    const auto *lh_345 = buffer.data(lh + 345);
    const auto *lh_346 = buffer.data(lh + 346);
    const auto *lh_347 = buffer.data(lh + 347);
    const auto *lh_348 = buffer.data(lh + 348);
    const auto *lh_349 = buffer.data(lh + 349);
    const auto *lh_350 = buffer.data(lh + 350);
    const auto *lh_351 = buffer.data(lh + 351);
    const auto *lh_352 = buffer.data(lh + 352);
    const auto *lh_353 = buffer.data(lh + 353);
    const auto *lh_354 = buffer.data(lh + 354);
    const auto *lh_355 = buffer.data(lh + 355);
    const auto *lh_356 = buffer.data(lh + 356);
    const auto *lh_357 = buffer.data(lh + 357);
    const auto *lh_358 = buffer.data(lh + 358);
    const auto *lh_359 = buffer.data(lh + 359);
    const auto *lh_360 = buffer.data(lh + 360);
    const auto *lh_361 = buffer.data(lh + 361);
    const auto *lh_362 = buffer.data(lh + 362);
    const auto *lh_363 = buffer.data(lh + 363);
    const auto *lh_364 = buffer.data(lh + 364);
    const auto *lh_365 = buffer.data(lh + 365);
    const auto *lh_366 = buffer.data(lh + 366);
    const auto *lh_367 = buffer.data(lh + 367);
    const auto *lh_368 = buffer.data(lh + 368);
    const auto *lh_369 = buffer.data(lh + 369);
    const auto *lh_370 = buffer.data(lh + 370);
    const auto *lh_371 = buffer.data(lh + 371);
    const auto *lh_372 = buffer.data(lh + 372);
    const auto *lh_373 = buffer.data(lh + 373);
    const auto *lh_374 = buffer.data(lh + 374);
    const auto *lh_375 = buffer.data(lh + 375);
    const auto *lh_376 = buffer.data(lh + 376);
    const auto *lh_377 = buffer.data(lh + 377);
    const auto *lh_378 = buffer.data(lh + 378);
    const auto *lh_379 = buffer.data(lh + 379);
    const auto *lh_380 = buffer.data(lh + 380);
    const auto *lh_381 = buffer.data(lh + 381);
    const auto *lh_382 = buffer.data(lh + 382);
    const auto *lh_383 = buffer.data(lh + 383);
    const auto *lh_384 = buffer.data(lh + 384);
    const auto *lh_385 = buffer.data(lh + 385);
    const auto *lh_386 = buffer.data(lh + 386);
    const auto *lh_387 = buffer.data(lh + 387);
    const auto *lh_388 = buffer.data(lh + 388);
    const auto *lh_389 = buffer.data(lh + 389);
    const auto *lh_390 = buffer.data(lh + 390);
    const auto *lh_391 = buffer.data(lh + 391);
    const auto *lh_392 = buffer.data(lh + 392);
    const auto *lh_393 = buffer.data(lh + 393);
    const auto *lh_394 = buffer.data(lh + 394);
    const auto *lh_395 = buffer.data(lh + 395);
    const auto *lh_396 = buffer.data(lh + 396);
    const auto *lh_397 = buffer.data(lh + 397);
    const auto *lh_398 = buffer.data(lh + 398);
    const auto *lh_399 = buffer.data(lh + 399);
    const auto *lh_400 = buffer.data(lh + 400);
    const auto *lh_401 = buffer.data(lh + 401);
    const auto *lh_402 = buffer.data(lh + 402);
    const auto *lh_403 = buffer.data(lh + 403);
    const auto *lh_404 = buffer.data(lh + 404);
    const auto *lh_405 = buffer.data(lh + 405);
    const auto *lh_406 = buffer.data(lh + 406);
    const auto *lh_407 = buffer.data(lh + 407);
    const auto *lh_408 = buffer.data(lh + 408);
    const auto *lh_409 = buffer.data(lh + 409);
    const auto *lh_410 = buffer.data(lh + 410);
    const auto *lh_411 = buffer.data(lh + 411);
    const auto *lh_412 = buffer.data(lh + 412);
    const auto *lh_413 = buffer.data(lh + 413);
    const auto *lh_414 = buffer.data(lh + 414);
    const auto *lh_415 = buffer.data(lh + 415);
    const auto *lh_416 = buffer.data(lh + 416);
    const auto *lh_417 = buffer.data(lh + 417);
    const auto *lh_418 = buffer.data(lh + 418);
    const auto *lh_419 = buffer.data(lh + 419);
    const auto *lh_420 = buffer.data(lh + 420);
    const auto *lh_421 = buffer.data(lh + 421);
    const auto *lh_422 = buffer.data(lh + 422);
    const auto *lh_423 = buffer.data(lh + 423);
    const auto *lh_424 = buffer.data(lh + 424);
    const auto *lh_425 = buffer.data(lh + 425);
    const auto *lh_426 = buffer.data(lh + 426);
    const auto *lh_427 = buffer.data(lh + 427);
    const auto *lh_428 = buffer.data(lh + 428);
    const auto *lh_429 = buffer.data(lh + 429);
    const auto *lh_430 = buffer.data(lh + 430);
    const auto *lh_431 = buffer.data(lh + 431);
    const auto *lh_432 = buffer.data(lh + 432);
    const auto *lh_433 = buffer.data(lh + 433);
    const auto *lh_434 = buffer.data(lh + 434);
    const auto *lh_435 = buffer.data(lh + 435);
    const auto *lh_436 = buffer.data(lh + 436);
    const auto *lh_437 = buffer.data(lh + 437);
    const auto *lh_438 = buffer.data(lh + 438);
    const auto *lh_439 = buffer.data(lh + 439);
    const auto *lh_440 = buffer.data(lh + 440);
    const auto *lh_462 = buffer.data(lh + 462);
    const auto *lh_463 = buffer.data(lh + 463);
    const auto *lh_464 = buffer.data(lh + 464);
    const auto *lh_465 = buffer.data(lh + 465);
    const auto *lh_466 = buffer.data(lh + 466);
    const auto *lh_467 = buffer.data(lh + 467);
    const auto *lh_468 = buffer.data(lh + 468);
    const auto *lh_469 = buffer.data(lh + 469);
    const auto *lh_470 = buffer.data(lh + 470);
    const auto *lh_471 = buffer.data(lh + 471);
    const auto *lh_472 = buffer.data(lh + 472);
    const auto *lh_473 = buffer.data(lh + 473);
    const auto *lh_474 = buffer.data(lh + 474);
    const auto *lh_475 = buffer.data(lh + 475);
    const auto *lh_476 = buffer.data(lh + 476);
    const auto *lh_477 = buffer.data(lh + 477);
    const auto *lh_478 = buffer.data(lh + 478);
    const auto *lh_479 = buffer.data(lh + 479);
    const auto *lh_480 = buffer.data(lh + 480);
    const auto *lh_481 = buffer.data(lh + 481);
    const auto *lh_482 = buffer.data(lh + 482);
    const auto *lh_483 = buffer.data(lh + 483);
    const auto *lh_484 = buffer.data(lh + 484);

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, ih_93, ih_94, ih_95, ih_96, ih_97, \
                         lh_282, lh_283, lh_284, lh_285, lh_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -2.0 * ih_93[k]
                   + f_0 * lh_282[k];

        t_178[k] = -2.0 * ih_94[k]
                   + f_0 * lh_283[k];

        t_179[k] = -2.0 * ih_95[k]
                   + f_0 * lh_284[k];

        t_180[k] = -2.0 * ih_96[k]
                   + f_0 * lh_285[k];

        t_181[k] = -2.0 * ih_97[k]
                   + f_0 * lh_286[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, ih_98, ih_99, ih_100, ih_101, \
                         ih_102, lh_287, lh_288, lh_289, lh_290, \
                         lh_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * ih_98[k]
                   + f_0 * lh_287[k];

        t_183[k] = -2.0 * ih_99[k]
                   + f_0 * lh_288[k];

        t_184[k] = -2.0 * ih_100[k]
                   + f_0 * lh_289[k];

        t_185[k] = -2.0 * ih_101[k]
                   + f_0 * lh_290[k];

        t_186[k] = -2.0 * ih_102[k]
                   + f_0 * lh_291[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, ih_103, ih_104, ih_105, ih_106, \
                         ih_107, lh_292, lh_293, lh_294, lh_295, \
                         lh_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * ih_103[k]
                   + f_0 * lh_292[k];

        t_188[k] = -2.0 * ih_104[k]
                   + f_0 * lh_293[k];

        t_189[k] = -3.0 * ih_105[k]
                   + f_0 * lh_294[k];

        t_190[k] = -3.0 * ih_106[k]
                   + f_0 * lh_295[k];

        t_191[k] = -3.0 * ih_107[k]
                   + f_0 * lh_296[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, ih_108, ih_109, ih_110, ih_111, \
                         ih_112, lh_297, lh_298, lh_299, lh_300, \
                         lh_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -3.0 * ih_108[k]
                   + f_0 * lh_297[k];

        t_193[k] = -3.0 * ih_109[k]
                   + f_0 * lh_298[k];

        t_194[k] = -3.0 * ih_110[k]
                   + f_0 * lh_299[k];

        t_195[k] = -3.0 * ih_111[k]
                   + f_0 * lh_300[k];

        t_196[k] = -3.0 * ih_112[k]
                   + f_0 * lh_301[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, ih_113, ih_114, ih_115, ih_116, \
                         ih_117, lh_302, lh_303, lh_304, lh_305, \
                         lh_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -3.0 * ih_113[k]
                   + f_0 * lh_302[k];

        t_198[k] = -3.0 * ih_114[k]
                   + f_0 * lh_303[k];

        t_199[k] = -3.0 * ih_115[k]
                   + f_0 * lh_304[k];

        t_200[k] = -3.0 * ih_116[k]
                   + f_0 * lh_305[k];

        t_201[k] = -3.0 * ih_117[k]
                   + f_0 * lh_306[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, ih_118, ih_119, ih_120, ih_121, \
                         ih_122, lh_307, lh_308, lh_309, lh_310, \
                         lh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -3.0 * ih_118[k]
                   + f_0 * lh_307[k];

        t_203[k] = -3.0 * ih_119[k]
                   + f_0 * lh_308[k];

        t_204[k] = -3.0 * ih_120[k]
                   + f_0 * lh_309[k];

        t_205[k] = -3.0 * ih_121[k]
                   + f_0 * lh_310[k];

        t_206[k] = -3.0 * ih_122[k]
                   + f_0 * lh_311[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, ih_123, ih_124, ih_125, \
                         lh_312, lh_313, lh_314, lh_336, lh_337, \
                         lh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -3.0 * ih_123[k]
                   + f_0 * lh_312[k];

        t_208[k] = -3.0 * ih_124[k]
                   + f_0 * lh_313[k];

        t_209[k] = -3.0 * ih_125[k]
                   + f_0 * lh_314[k];

        t_210[k] = f_0 * lh_336[k];

        t_211[k] = f_0 * lh_337[k];

        t_212[k] = f_0 * lh_338[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, t_218, t_219, t_220, lh_339, \
                         lh_340, lh_341, lh_342, lh_343, lh_344, lh_345, \
                         lh_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_0 * lh_339[k];

        t_214[k] = f_0 * lh_340[k];

        t_215[k] = f_0 * lh_341[k];

        t_216[k] = f_0 * lh_342[k];

        t_217[k] = f_0 * lh_343[k];

        t_218[k] = f_0 * lh_344[k];

        t_219[k] = f_0 * lh_345[k];

        t_220[k] = f_0 * lh_346[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, t_226, t_227, t_228, lh_347, \
                         lh_348, lh_349, lh_350, lh_351, lh_352, lh_353, \
                         lh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_0 * lh_347[k];

        t_222[k] = f_0 * lh_348[k];

        t_223[k] = f_0 * lh_349[k];

        t_224[k] = f_0 * lh_350[k];

        t_225[k] = f_0 * lh_351[k];

        t_226[k] = f_0 * lh_352[k];

        t_227[k] = f_0 * lh_353[k];

        t_228[k] = f_0 * lh_354[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, t_234, ih_126, ih_127, ih_128, \
                         ih_129, lh_355, lh_356, lh_357, lh_358, lh_359, \
                         lh_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_0 * lh_355[k];

        t_230[k] = f_0 * lh_356[k];

        t_231[k] = -ih_126[k]
                   + f_0 * lh_357[k];

        t_232[k] = -ih_127[k]
                   + f_0 * lh_358[k];

        t_233[k] = -ih_128[k]
                   + f_0 * lh_359[k];

        t_234[k] = -ih_129[k]
                   + f_0 * lh_360[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ih_130, ih_131, ih_132, ih_133, \
                         ih_134, lh_361, lh_362, lh_363, lh_364, \
                         lh_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -ih_130[k]
                   + f_0 * lh_361[k];

        t_236[k] = -ih_131[k]
                   + f_0 * lh_362[k];

        t_237[k] = -ih_132[k]
                   + f_0 * lh_363[k];

        t_238[k] = -ih_133[k]
                   + f_0 * lh_364[k];

        t_239[k] = -ih_134[k]
                   + f_0 * lh_365[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ih_135, ih_136, ih_137, ih_138, \
                         ih_139, lh_366, lh_367, lh_368, lh_369, \
                         lh_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -ih_135[k]
                   + f_0 * lh_366[k];

        t_241[k] = -ih_136[k]
                   + f_0 * lh_367[k];

        t_242[k] = -ih_137[k]
                   + f_0 * lh_368[k];

        t_243[k] = -ih_138[k]
                   + f_0 * lh_369[k];

        t_244[k] = -ih_139[k]
                   + f_0 * lh_370[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ih_140, ih_141, ih_142, ih_143, \
                         ih_144, lh_371, lh_372, lh_373, lh_374, \
                         lh_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -ih_140[k]
                   + f_0 * lh_371[k];

        t_246[k] = -ih_141[k]
                   + f_0 * lh_372[k];

        t_247[k] = -ih_142[k]
                   + f_0 * lh_373[k];

        t_248[k] = -ih_143[k]
                   + f_0 * lh_374[k];

        t_249[k] = -ih_144[k]
                   + f_0 * lh_375[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ih_145, ih_146, ih_147, ih_148, \
                         ih_149, lh_376, lh_377, lh_378, lh_379, \
                         lh_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -ih_145[k]
                   + f_0 * lh_376[k];

        t_251[k] = -ih_146[k]
                   + f_0 * lh_377[k];

        t_252[k] = -2.0 * ih_147[k]
                   + f_0 * lh_378[k];

        t_253[k] = -2.0 * ih_148[k]
                   + f_0 * lh_379[k];

        t_254[k] = -2.0 * ih_149[k]
                   + f_0 * lh_380[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ih_150, ih_151, ih_152, ih_153, \
                         ih_154, lh_381, lh_382, lh_383, lh_384, \
                         lh_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * ih_150[k]
                   + f_0 * lh_381[k];

        t_256[k] = -2.0 * ih_151[k]
                   + f_0 * lh_382[k];

        t_257[k] = -2.0 * ih_152[k]
                   + f_0 * lh_383[k];

        t_258[k] = -2.0 * ih_153[k]
                   + f_0 * lh_384[k];

        t_259[k] = -2.0 * ih_154[k]
                   + f_0 * lh_385[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ih_155, ih_156, ih_157, ih_158, \
                         ih_159, lh_386, lh_387, lh_388, lh_389, \
                         lh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * ih_155[k]
                   + f_0 * lh_386[k];

        t_261[k] = -2.0 * ih_156[k]
                   + f_0 * lh_387[k];

        t_262[k] = -2.0 * ih_157[k]
                   + f_0 * lh_388[k];

        t_263[k] = -2.0 * ih_158[k]
                   + f_0 * lh_389[k];

        t_264[k] = -2.0 * ih_159[k]
                   + f_0 * lh_390[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ih_160, ih_161, ih_162, ih_163, \
                         ih_164, lh_391, lh_392, lh_393, lh_394, \
                         lh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * ih_160[k]
                   + f_0 * lh_391[k];

        t_266[k] = -2.0 * ih_161[k]
                   + f_0 * lh_392[k];

        t_267[k] = -2.0 * ih_162[k]
                   + f_0 * lh_393[k];

        t_268[k] = -2.0 * ih_163[k]
                   + f_0 * lh_394[k];

        t_269[k] = -2.0 * ih_164[k]
                   + f_0 * lh_395[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ih_165, ih_166, ih_167, ih_168, \
                         ih_169, lh_396, lh_397, lh_398, lh_399, \
                         lh_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * ih_165[k]
                   + f_0 * lh_396[k];

        t_271[k] = -2.0 * ih_166[k]
                   + f_0 * lh_397[k];

        t_272[k] = -2.0 * ih_167[k]
                   + f_0 * lh_398[k];

        t_273[k] = -3.0 * ih_168[k]
                   + f_0 * lh_399[k];

        t_274[k] = -3.0 * ih_169[k]
                   + f_0 * lh_400[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ih_170, ih_171, ih_172, ih_173, \
                         ih_174, lh_401, lh_402, lh_403, lh_404, \
                         lh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * ih_170[k]
                   + f_0 * lh_401[k];

        t_276[k] = -3.0 * ih_171[k]
                   + f_0 * lh_402[k];

        t_277[k] = -3.0 * ih_172[k]
                   + f_0 * lh_403[k];

        t_278[k] = -3.0 * ih_173[k]
                   + f_0 * lh_404[k];

        t_279[k] = -3.0 * ih_174[k]
                   + f_0 * lh_405[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ih_175, ih_176, ih_177, ih_178, \
                         ih_179, lh_406, lh_407, lh_408, lh_409, \
                         lh_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * ih_175[k]
                   + f_0 * lh_406[k];

        t_281[k] = -3.0 * ih_176[k]
                   + f_0 * lh_407[k];

        t_282[k] = -3.0 * ih_177[k]
                   + f_0 * lh_408[k];

        t_283[k] = -3.0 * ih_178[k]
                   + f_0 * lh_409[k];

        t_284[k] = -3.0 * ih_179[k]
                   + f_0 * lh_410[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ih_180, ih_181, ih_182, ih_183, \
                         ih_184, lh_411, lh_412, lh_413, lh_414, \
                         lh_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * ih_180[k]
                   + f_0 * lh_411[k];

        t_286[k] = -3.0 * ih_181[k]
                   + f_0 * lh_412[k];

        t_287[k] = -3.0 * ih_182[k]
                   + f_0 * lh_413[k];

        t_288[k] = -3.0 * ih_183[k]
                   + f_0 * lh_414[k];

        t_289[k] = -3.0 * ih_184[k]
                   + f_0 * lh_415[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ih_185, ih_186, ih_187, ih_188, \
                         ih_189, lh_416, lh_417, lh_418, lh_419, \
                         lh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * ih_185[k]
                   + f_0 * lh_416[k];

        t_291[k] = -3.0 * ih_186[k]
                   + f_0 * lh_417[k];

        t_292[k] = -3.0 * ih_187[k]
                   + f_0 * lh_418[k];

        t_293[k] = -3.0 * ih_188[k]
                   + f_0 * lh_419[k];

        t_294[k] = -4.0 * ih_189[k]
                   + f_0 * lh_420[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ih_190, ih_191, ih_192, ih_193, \
                         ih_194, lh_421, lh_422, lh_423, lh_424, \
                         lh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * ih_190[k]
                   + f_0 * lh_421[k];

        t_296[k] = -4.0 * ih_191[k]
                   + f_0 * lh_422[k];

        t_297[k] = -4.0 * ih_192[k]
                   + f_0 * lh_423[k];

        t_298[k] = -4.0 * ih_193[k]
                   + f_0 * lh_424[k];

        t_299[k] = -4.0 * ih_194[k]
                   + f_0 * lh_425[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ih_195, ih_196, ih_197, ih_198, \
                         ih_199, lh_426, lh_427, lh_428, lh_429, \
                         lh_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * ih_195[k]
                   + f_0 * lh_426[k];

        t_301[k] = -4.0 * ih_196[k]
                   + f_0 * lh_427[k];

        t_302[k] = -4.0 * ih_197[k]
                   + f_0 * lh_428[k];

        t_303[k] = -4.0 * ih_198[k]
                   + f_0 * lh_429[k];

        t_304[k] = -4.0 * ih_199[k]
                   + f_0 * lh_430[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ih_200, ih_201, ih_202, ih_203, \
                         ih_204, lh_431, lh_432, lh_433, lh_434, \
                         lh_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * ih_200[k]
                   + f_0 * lh_431[k];

        t_306[k] = -4.0 * ih_201[k]
                   + f_0 * lh_432[k];

        t_307[k] = -4.0 * ih_202[k]
                   + f_0 * lh_433[k];

        t_308[k] = -4.0 * ih_203[k]
                   + f_0 * lh_434[k];

        t_309[k] = -4.0 * ih_204[k]
                   + f_0 * lh_435[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ih_205, ih_206, ih_207, ih_208, \
                         ih_209, lh_436, lh_437, lh_438, lh_439, \
                         lh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -4.0 * ih_205[k]
                   + f_0 * lh_436[k];

        t_311[k] = -4.0 * ih_206[k]
                   + f_0 * lh_437[k];

        t_312[k] = -4.0 * ih_207[k]
                   + f_0 * lh_438[k];

        t_313[k] = -4.0 * ih_208[k]
                   + f_0 * lh_439[k];

        t_314[k] = -4.0 * ih_209[k]
                   + f_0 * lh_440[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, t_320, t_321, t_322, lh_462, \
                         lh_463, lh_464, lh_465, lh_466, lh_467, lh_468, \
                         lh_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_0 * lh_462[k];

        t_316[k] = f_0 * lh_463[k];

        t_317[k] = f_0 * lh_464[k];

        t_318[k] = f_0 * lh_465[k];

        t_319[k] = f_0 * lh_466[k];

        t_320[k] = f_0 * lh_467[k];

        t_321[k] = f_0 * lh_468[k];

        t_322[k] = f_0 * lh_469[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, t_328, t_329, t_330, lh_470, \
                         lh_471, lh_472, lh_473, lh_474, lh_475, lh_476, \
                         lh_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_0 * lh_470[k];

        t_324[k] = f_0 * lh_471[k];

        t_325[k] = f_0 * lh_472[k];

        t_326[k] = f_0 * lh_473[k];

        t_327[k] = f_0 * lh_474[k];

        t_328[k] = f_0 * lh_475[k];

        t_329[k] = f_0 * lh_476[k];

        t_330[k] = f_0 * lh_477[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, t_335, t_336, t_337, ih_210, ih_211, \
                         lh_478, lh_479, lh_480, lh_481, lh_482, lh_483, \
                         lh_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_0 * lh_478[k];

        t_332[k] = f_0 * lh_479[k];

        t_333[k] = f_0 * lh_480[k];

        t_334[k] = f_0 * lh_481[k];

        t_335[k] = f_0 * lh_482[k];

        t_336[k] = -ih_210[k]
                   + f_0 * lh_483[k];

        t_337[k] = -ih_211[k]
                   + f_0 * lh_484[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_495 = buffer.data(target + 495);

    const auto *ih_212 = buffer.data(ih + 212);
    const auto *ih_213 = buffer.data(ih + 213);
    const auto *ih_214 = buffer.data(ih + 214);
    const auto *ih_215 = buffer.data(ih + 215);
    const auto *ih_216 = buffer.data(ih + 216);
    const auto *ih_217 = buffer.data(ih + 217);
    const auto *ih_218 = buffer.data(ih + 218);
    const auto *ih_219 = buffer.data(ih + 219);
    const auto *ih_220 = buffer.data(ih + 220);
    const auto *ih_221 = buffer.data(ih + 221);
    const auto *ih_222 = buffer.data(ih + 222);
    const auto *ih_223 = buffer.data(ih + 223);
    const auto *ih_224 = buffer.data(ih + 224);
    const auto *ih_225 = buffer.data(ih + 225);
    const auto *ih_226 = buffer.data(ih + 226);
    const auto *ih_227 = buffer.data(ih + 227);
    const auto *ih_228 = buffer.data(ih + 228);
    const auto *ih_229 = buffer.data(ih + 229);
    const auto *ih_230 = buffer.data(ih + 230);
    const auto *ih_231 = buffer.data(ih + 231);
    const auto *ih_232 = buffer.data(ih + 232);
    const auto *ih_233 = buffer.data(ih + 233);
    const auto *ih_234 = buffer.data(ih + 234);
    const auto *ih_235 = buffer.data(ih + 235);
    const auto *ih_236 = buffer.data(ih + 236);
    const auto *ih_237 = buffer.data(ih + 237);
    const auto *ih_238 = buffer.data(ih + 238);
    const auto *ih_239 = buffer.data(ih + 239);
    const auto *ih_240 = buffer.data(ih + 240);
    const auto *ih_241 = buffer.data(ih + 241);
    const auto *ih_242 = buffer.data(ih + 242);
    const auto *ih_243 = buffer.data(ih + 243);
    const auto *ih_244 = buffer.data(ih + 244);
    const auto *ih_245 = buffer.data(ih + 245);
    const auto *ih_246 = buffer.data(ih + 246);
    const auto *ih_247 = buffer.data(ih + 247);
    const auto *ih_248 = buffer.data(ih + 248);
    const auto *ih_249 = buffer.data(ih + 249);
    const auto *ih_250 = buffer.data(ih + 250);
    const auto *ih_251 = buffer.data(ih + 251);
    const auto *ih_252 = buffer.data(ih + 252);
    const auto *ih_253 = buffer.data(ih + 253);
    const auto *ih_254 = buffer.data(ih + 254);
    const auto *ih_255 = buffer.data(ih + 255);
    const auto *ih_256 = buffer.data(ih + 256);
    const auto *ih_257 = buffer.data(ih + 257);
    const auto *ih_258 = buffer.data(ih + 258);
    const auto *ih_259 = buffer.data(ih + 259);
    const auto *ih_260 = buffer.data(ih + 260);
    const auto *ih_261 = buffer.data(ih + 261);
    const auto *ih_262 = buffer.data(ih + 262);
    const auto *ih_263 = buffer.data(ih + 263);
    const auto *ih_264 = buffer.data(ih + 264);
    const auto *ih_265 = buffer.data(ih + 265);
    const auto *ih_266 = buffer.data(ih + 266);
    const auto *ih_267 = buffer.data(ih + 267);
    const auto *ih_268 = buffer.data(ih + 268);
    const auto *ih_269 = buffer.data(ih + 269);
    const auto *ih_270 = buffer.data(ih + 270);
    const auto *ih_271 = buffer.data(ih + 271);
    const auto *ih_272 = buffer.data(ih + 272);
    const auto *ih_273 = buffer.data(ih + 273);
    const auto *ih_274 = buffer.data(ih + 274);
    const auto *ih_275 = buffer.data(ih + 275);
    const auto *ih_276 = buffer.data(ih + 276);
    const auto *ih_277 = buffer.data(ih + 277);
    const auto *ih_278 = buffer.data(ih + 278);
    const auto *ih_279 = buffer.data(ih + 279);
    const auto *ih_280 = buffer.data(ih + 280);
    const auto *ih_281 = buffer.data(ih + 281);
    const auto *ih_282 = buffer.data(ih + 282);
    const auto *ih_283 = buffer.data(ih + 283);
    const auto *ih_284 = buffer.data(ih + 284);
    const auto *ih_285 = buffer.data(ih + 285);
    const auto *ih_286 = buffer.data(ih + 286);
    const auto *ih_287 = buffer.data(ih + 287);
    const auto *ih_288 = buffer.data(ih + 288);
    const auto *ih_289 = buffer.data(ih + 289);
    const auto *ih_290 = buffer.data(ih + 290);
    const auto *ih_291 = buffer.data(ih + 291);
    const auto *ih_292 = buffer.data(ih + 292);
    const auto *ih_293 = buffer.data(ih + 293);
    const auto *ih_294 = buffer.data(ih + 294);
    const auto *ih_295 = buffer.data(ih + 295);
    const auto *ih_296 = buffer.data(ih + 296);
    const auto *ih_297 = buffer.data(ih + 297);
    const auto *ih_298 = buffer.data(ih + 298);
    const auto *ih_299 = buffer.data(ih + 299);
    const auto *ih_300 = buffer.data(ih + 300);
    const auto *ih_301 = buffer.data(ih + 301);
    const auto *ih_302 = buffer.data(ih + 302);
    const auto *ih_303 = buffer.data(ih + 303);
    const auto *ih_304 = buffer.data(ih + 304);
    const auto *ih_305 = buffer.data(ih + 305);
    const auto *ih_306 = buffer.data(ih + 306);
    const auto *ih_307 = buffer.data(ih + 307);
    const auto *ih_308 = buffer.data(ih + 308);
    const auto *ih_309 = buffer.data(ih + 309);
    const auto *ih_310 = buffer.data(ih + 310);
    const auto *ih_311 = buffer.data(ih + 311);
    const auto *ih_312 = buffer.data(ih + 312);
    const auto *ih_313 = buffer.data(ih + 313);
    const auto *ih_314 = buffer.data(ih + 314);
    const auto *ih_315 = buffer.data(ih + 315);
    const auto *ih_316 = buffer.data(ih + 316);
    const auto *ih_317 = buffer.data(ih + 317);
    const auto *ih_318 = buffer.data(ih + 318);
    const auto *ih_319 = buffer.data(ih + 319);
    const auto *ih_320 = buffer.data(ih + 320);
    const auto *ih_321 = buffer.data(ih + 321);
    const auto *ih_322 = buffer.data(ih + 322);
    const auto *ih_323 = buffer.data(ih + 323);
    const auto *ih_324 = buffer.data(ih + 324);
    const auto *ih_325 = buffer.data(ih + 325);
    const auto *ih_326 = buffer.data(ih + 326);
    const auto *ih_327 = buffer.data(ih + 327);
    const auto *ih_328 = buffer.data(ih + 328);
    const auto *ih_329 = buffer.data(ih + 329);
    const auto *ih_330 = buffer.data(ih + 330);
    const auto *ih_331 = buffer.data(ih + 331);
    const auto *ih_332 = buffer.data(ih + 332);
    const auto *ih_333 = buffer.data(ih + 333);
    const auto *ih_334 = buffer.data(ih + 334);
    const auto *ih_335 = buffer.data(ih + 335);
    const auto *ih_336 = buffer.data(ih + 336);
    const auto *ih_337 = buffer.data(ih + 337);
    const auto *ih_338 = buffer.data(ih + 338);
    const auto *ih_339 = buffer.data(ih + 339);
    const auto *ih_340 = buffer.data(ih + 340);
    const auto *ih_341 = buffer.data(ih + 341);
    const auto *ih_342 = buffer.data(ih + 342);
    const auto *ih_343 = buffer.data(ih + 343);
    const auto *ih_344 = buffer.data(ih + 344);
    const auto *ih_345 = buffer.data(ih + 345);
    const auto *ih_346 = buffer.data(ih + 346);
    const auto *ih_347 = buffer.data(ih + 347);
    const auto *ih_348 = buffer.data(ih + 348);

    const auto *lh_485 = buffer.data(lh + 485);
    const auto *lh_486 = buffer.data(lh + 486);
    const auto *lh_487 = buffer.data(lh + 487);
    const auto *lh_488 = buffer.data(lh + 488);
    const auto *lh_489 = buffer.data(lh + 489);
    const auto *lh_490 = buffer.data(lh + 490);
    const auto *lh_491 = buffer.data(lh + 491);
    const auto *lh_492 = buffer.data(lh + 492);
    const auto *lh_493 = buffer.data(lh + 493);
    const auto *lh_494 = buffer.data(lh + 494);
    const auto *lh_495 = buffer.data(lh + 495);
    const auto *lh_496 = buffer.data(lh + 496);
    const auto *lh_497 = buffer.data(lh + 497);
    const auto *lh_498 = buffer.data(lh + 498);
    const auto *lh_499 = buffer.data(lh + 499);
    const auto *lh_500 = buffer.data(lh + 500);
    const auto *lh_501 = buffer.data(lh + 501);
    const auto *lh_502 = buffer.data(lh + 502);
    const auto *lh_503 = buffer.data(lh + 503);
    const auto *lh_504 = buffer.data(lh + 504);
    const auto *lh_505 = buffer.data(lh + 505);
    const auto *lh_506 = buffer.data(lh + 506);
    const auto *lh_507 = buffer.data(lh + 507);
    const auto *lh_508 = buffer.data(lh + 508);
    const auto *lh_509 = buffer.data(lh + 509);
    const auto *lh_510 = buffer.data(lh + 510);
    const auto *lh_511 = buffer.data(lh + 511);
    const auto *lh_512 = buffer.data(lh + 512);
    const auto *lh_513 = buffer.data(lh + 513);
    const auto *lh_514 = buffer.data(lh + 514);
    const auto *lh_515 = buffer.data(lh + 515);
    const auto *lh_516 = buffer.data(lh + 516);
    const auto *lh_517 = buffer.data(lh + 517);
    const auto *lh_518 = buffer.data(lh + 518);
    const auto *lh_519 = buffer.data(lh + 519);
    const auto *lh_520 = buffer.data(lh + 520);
    const auto *lh_521 = buffer.data(lh + 521);
    const auto *lh_522 = buffer.data(lh + 522);
    const auto *lh_523 = buffer.data(lh + 523);
    const auto *lh_524 = buffer.data(lh + 524);
    const auto *lh_525 = buffer.data(lh + 525);
    const auto *lh_526 = buffer.data(lh + 526);
    const auto *lh_527 = buffer.data(lh + 527);
    const auto *lh_528 = buffer.data(lh + 528);
    const auto *lh_529 = buffer.data(lh + 529);
    const auto *lh_530 = buffer.data(lh + 530);
    const auto *lh_531 = buffer.data(lh + 531);
    const auto *lh_532 = buffer.data(lh + 532);
    const auto *lh_533 = buffer.data(lh + 533);
    const auto *lh_534 = buffer.data(lh + 534);
    const auto *lh_535 = buffer.data(lh + 535);
    const auto *lh_536 = buffer.data(lh + 536);
    const auto *lh_537 = buffer.data(lh + 537);
    const auto *lh_538 = buffer.data(lh + 538);
    const auto *lh_539 = buffer.data(lh + 539);
    const auto *lh_540 = buffer.data(lh + 540);
    const auto *lh_541 = buffer.data(lh + 541);
    const auto *lh_542 = buffer.data(lh + 542);
    const auto *lh_543 = buffer.data(lh + 543);
    const auto *lh_544 = buffer.data(lh + 544);
    const auto *lh_545 = buffer.data(lh + 545);
    const auto *lh_546 = buffer.data(lh + 546);
    const auto *lh_547 = buffer.data(lh + 547);
    const auto *lh_548 = buffer.data(lh + 548);
    const auto *lh_549 = buffer.data(lh + 549);
    const auto *lh_550 = buffer.data(lh + 550);
    const auto *lh_551 = buffer.data(lh + 551);
    const auto *lh_552 = buffer.data(lh + 552);
    const auto *lh_553 = buffer.data(lh + 553);
    const auto *lh_554 = buffer.data(lh + 554);
    const auto *lh_555 = buffer.data(lh + 555);
    const auto *lh_556 = buffer.data(lh + 556);
    const auto *lh_557 = buffer.data(lh + 557);
    const auto *lh_558 = buffer.data(lh + 558);
    const auto *lh_559 = buffer.data(lh + 559);
    const auto *lh_560 = buffer.data(lh + 560);
    const auto *lh_561 = buffer.data(lh + 561);
    const auto *lh_562 = buffer.data(lh + 562);
    const auto *lh_563 = buffer.data(lh + 563);
    const auto *lh_564 = buffer.data(lh + 564);
    const auto *lh_565 = buffer.data(lh + 565);
    const auto *lh_566 = buffer.data(lh + 566);
    const auto *lh_567 = buffer.data(lh + 567);
    const auto *lh_568 = buffer.data(lh + 568);
    const auto *lh_569 = buffer.data(lh + 569);
    const auto *lh_570 = buffer.data(lh + 570);
    const auto *lh_571 = buffer.data(lh + 571);
    const auto *lh_572 = buffer.data(lh + 572);
    const auto *lh_573 = buffer.data(lh + 573);
    const auto *lh_574 = buffer.data(lh + 574);
    const auto *lh_575 = buffer.data(lh + 575);
    const auto *lh_576 = buffer.data(lh + 576);
    const auto *lh_577 = buffer.data(lh + 577);
    const auto *lh_578 = buffer.data(lh + 578);
    const auto *lh_579 = buffer.data(lh + 579);
    const auto *lh_580 = buffer.data(lh + 580);
    const auto *lh_581 = buffer.data(lh + 581);
    const auto *lh_582 = buffer.data(lh + 582);
    const auto *lh_583 = buffer.data(lh + 583);
    const auto *lh_584 = buffer.data(lh + 584);
    const auto *lh_585 = buffer.data(lh + 585);
    const auto *lh_586 = buffer.data(lh + 586);
    const auto *lh_587 = buffer.data(lh + 587);
    const auto *lh_609 = buffer.data(lh + 609);
    const auto *lh_610 = buffer.data(lh + 610);
    const auto *lh_611 = buffer.data(lh + 611);
    const auto *lh_612 = buffer.data(lh + 612);
    const auto *lh_613 = buffer.data(lh + 613);
    const auto *lh_614 = buffer.data(lh + 614);
    const auto *lh_615 = buffer.data(lh + 615);
    const auto *lh_616 = buffer.data(lh + 616);
    const auto *lh_617 = buffer.data(lh + 617);
    const auto *lh_618 = buffer.data(lh + 618);
    const auto *lh_619 = buffer.data(lh + 619);
    const auto *lh_620 = buffer.data(lh + 620);
    const auto *lh_621 = buffer.data(lh + 621);
    const auto *lh_622 = buffer.data(lh + 622);
    const auto *lh_623 = buffer.data(lh + 623);
    const auto *lh_624 = buffer.data(lh + 624);
    const auto *lh_625 = buffer.data(lh + 625);
    const auto *lh_626 = buffer.data(lh + 626);
    const auto *lh_627 = buffer.data(lh + 627);
    const auto *lh_628 = buffer.data(lh + 628);
    const auto *lh_629 = buffer.data(lh + 629);
    const auto *lh_630 = buffer.data(lh + 630);
    const auto *lh_631 = buffer.data(lh + 631);
    const auto *lh_632 = buffer.data(lh + 632);
    const auto *lh_633 = buffer.data(lh + 633);
    const auto *lh_634 = buffer.data(lh + 634);
    const auto *lh_635 = buffer.data(lh + 635);
    const auto *lh_636 = buffer.data(lh + 636);
    const auto *lh_637 = buffer.data(lh + 637);
    const auto *lh_638 = buffer.data(lh + 638);
    const auto *lh_639 = buffer.data(lh + 639);
    const auto *lh_640 = buffer.data(lh + 640);
    const auto *lh_641 = buffer.data(lh + 641);
    const auto *lh_642 = buffer.data(lh + 642);
    const auto *lh_643 = buffer.data(lh + 643);
    const auto *lh_644 = buffer.data(lh + 644);
    const auto *lh_645 = buffer.data(lh + 645);
    const auto *lh_646 = buffer.data(lh + 646);
    const auto *lh_647 = buffer.data(lh + 647);
    const auto *lh_648 = buffer.data(lh + 648);
    const auto *lh_649 = buffer.data(lh + 649);
    const auto *lh_650 = buffer.data(lh + 650);
    const auto *lh_651 = buffer.data(lh + 651);
    const auto *lh_652 = buffer.data(lh + 652);
    const auto *lh_653 = buffer.data(lh + 653);
    const auto *lh_654 = buffer.data(lh + 654);
    const auto *lh_655 = buffer.data(lh + 655);
    const auto *lh_656 = buffer.data(lh + 656);
    const auto *lh_657 = buffer.data(lh + 657);
    const auto *lh_658 = buffer.data(lh + 658);
    const auto *lh_659 = buffer.data(lh + 659);
    const auto *lh_660 = buffer.data(lh + 660);
    const auto *lh_661 = buffer.data(lh + 661);
    const auto *lh_662 = buffer.data(lh + 662);
    const auto *lh_663 = buffer.data(lh + 663);

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, ih_212, ih_213, ih_214, ih_215, \
                         ih_216, lh_485, lh_486, lh_487, lh_488, \
                         lh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = -ih_212[k]
                   + f_0 * lh_485[k];

        t_339[k] = -ih_213[k]
                   + f_0 * lh_486[k];

        t_340[k] = -ih_214[k]
                   + f_0 * lh_487[k];

        t_341[k] = -ih_215[k]
                   + f_0 * lh_488[k];

        t_342[k] = -ih_216[k]
                   + f_0 * lh_489[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, ih_217, ih_218, ih_219, ih_220, \
                         ih_221, lh_490, lh_491, lh_492, lh_493, \
                         lh_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = -ih_217[k]
                   + f_0 * lh_490[k];

        t_344[k] = -ih_218[k]
                   + f_0 * lh_491[k];

        t_345[k] = -ih_219[k]
                   + f_0 * lh_492[k];

        t_346[k] = -ih_220[k]
                   + f_0 * lh_493[k];

        t_347[k] = -ih_221[k]
                   + f_0 * lh_494[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, ih_222, ih_223, ih_224, ih_225, \
                         ih_226, lh_495, lh_496, lh_497, lh_498, \
                         lh_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = -ih_222[k]
                   + f_0 * lh_495[k];

        t_349[k] = -ih_223[k]
                   + f_0 * lh_496[k];

        t_350[k] = -ih_224[k]
                   + f_0 * lh_497[k];

        t_351[k] = -ih_225[k]
                   + f_0 * lh_498[k];

        t_352[k] = -ih_226[k]
                   + f_0 * lh_499[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, ih_227, ih_228, ih_229, ih_230, \
                         ih_231, lh_500, lh_501, lh_502, lh_503, \
                         lh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = -ih_227[k]
                   + f_0 * lh_500[k];

        t_354[k] = -ih_228[k]
                   + f_0 * lh_501[k];

        t_355[k] = -ih_229[k]
                   + f_0 * lh_502[k];

        t_356[k] = -ih_230[k]
                   + f_0 * lh_503[k];

        t_357[k] = -2.0 * ih_231[k]
                   + f_0 * lh_504[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, ih_232, ih_233, ih_234, ih_235, \
                         ih_236, lh_505, lh_506, lh_507, lh_508, \
                         lh_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = -2.0 * ih_232[k]
                   + f_0 * lh_505[k];

        t_359[k] = -2.0 * ih_233[k]
                   + f_0 * lh_506[k];

        t_360[k] = -2.0 * ih_234[k]
                   + f_0 * lh_507[k];

        t_361[k] = -2.0 * ih_235[k]
                   + f_0 * lh_508[k];

        t_362[k] = -2.0 * ih_236[k]
                   + f_0 * lh_509[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, ih_237, ih_238, ih_239, ih_240, \
                         ih_241, lh_510, lh_511, lh_512, lh_513, \
                         lh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = -2.0 * ih_237[k]
                   + f_0 * lh_510[k];

        t_364[k] = -2.0 * ih_238[k]
                   + f_0 * lh_511[k];

        t_365[k] = -2.0 * ih_239[k]
                   + f_0 * lh_512[k];

        t_366[k] = -2.0 * ih_240[k]
                   + f_0 * lh_513[k];

        t_367[k] = -2.0 * ih_241[k]
                   + f_0 * lh_514[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, ih_242, ih_243, ih_244, ih_245, \
                         ih_246, lh_515, lh_516, lh_517, lh_518, \
                         lh_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = -2.0 * ih_242[k]
                   + f_0 * lh_515[k];

        t_369[k] = -2.0 * ih_243[k]
                   + f_0 * lh_516[k];

        t_370[k] = -2.0 * ih_244[k]
                   + f_0 * lh_517[k];

        t_371[k] = -2.0 * ih_245[k]
                   + f_0 * lh_518[k];

        t_372[k] = -2.0 * ih_246[k]
                   + f_0 * lh_519[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, ih_247, ih_248, ih_249, ih_250, \
                         ih_251, lh_520, lh_521, lh_522, lh_523, \
                         lh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = -2.0 * ih_247[k]
                   + f_0 * lh_520[k];

        t_374[k] = -2.0 * ih_248[k]
                   + f_0 * lh_521[k];

        t_375[k] = -2.0 * ih_249[k]
                   + f_0 * lh_522[k];

        t_376[k] = -2.0 * ih_250[k]
                   + f_0 * lh_523[k];

        t_377[k] = -2.0 * ih_251[k]
                   + f_0 * lh_524[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, ih_252, ih_253, ih_254, ih_255, \
                         ih_256, lh_525, lh_526, lh_527, lh_528, \
                         lh_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = -3.0 * ih_252[k]
                   + f_0 * lh_525[k];

        t_379[k] = -3.0 * ih_253[k]
                   + f_0 * lh_526[k];

        t_380[k] = -3.0 * ih_254[k]
                   + f_0 * lh_527[k];

        t_381[k] = -3.0 * ih_255[k]
                   + f_0 * lh_528[k];

        t_382[k] = -3.0 * ih_256[k]
                   + f_0 * lh_529[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, ih_257, ih_258, ih_259, ih_260, \
                         ih_261, lh_530, lh_531, lh_532, lh_533, \
                         lh_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = -3.0 * ih_257[k]
                   + f_0 * lh_530[k];

        t_384[k] = -3.0 * ih_258[k]
                   + f_0 * lh_531[k];

        t_385[k] = -3.0 * ih_259[k]
                   + f_0 * lh_532[k];

        t_386[k] = -3.0 * ih_260[k]
                   + f_0 * lh_533[k];

        t_387[k] = -3.0 * ih_261[k]
                   + f_0 * lh_534[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, ih_262, ih_263, ih_264, ih_265, \
                         ih_266, lh_535, lh_536, lh_537, lh_538, \
                         lh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = -3.0 * ih_262[k]
                   + f_0 * lh_535[k];

        t_389[k] = -3.0 * ih_263[k]
                   + f_0 * lh_536[k];

        t_390[k] = -3.0 * ih_264[k]
                   + f_0 * lh_537[k];

        t_391[k] = -3.0 * ih_265[k]
                   + f_0 * lh_538[k];

        t_392[k] = -3.0 * ih_266[k]
                   + f_0 * lh_539[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, t_397, ih_267, ih_268, ih_269, ih_270, \
                         ih_271, lh_540, lh_541, lh_542, lh_543, \
                         lh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = -3.0 * ih_267[k]
                   + f_0 * lh_540[k];

        t_394[k] = -3.0 * ih_268[k]
                   + f_0 * lh_541[k];

        t_395[k] = -3.0 * ih_269[k]
                   + f_0 * lh_542[k];

        t_396[k] = -3.0 * ih_270[k]
                   + f_0 * lh_543[k];

        t_397[k] = -3.0 * ih_271[k]
                   + f_0 * lh_544[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, t_402, ih_272, ih_273, ih_274, ih_275, \
                         ih_276, lh_545, lh_546, lh_547, lh_548, \
                         lh_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = -3.0 * ih_272[k]
                   + f_0 * lh_545[k];

        t_399[k] = -4.0 * ih_273[k]
                   + f_0 * lh_546[k];

        t_400[k] = -4.0 * ih_274[k]
                   + f_0 * lh_547[k];

        t_401[k] = -4.0 * ih_275[k]
                   + f_0 * lh_548[k];

        t_402[k] = -4.0 * ih_276[k]
                   + f_0 * lh_549[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, t_407, ih_277, ih_278, ih_279, ih_280, \
                         ih_281, lh_550, lh_551, lh_552, lh_553, \
                         lh_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = -4.0 * ih_277[k]
                   + f_0 * lh_550[k];

        t_404[k] = -4.0 * ih_278[k]
                   + f_0 * lh_551[k];

        t_405[k] = -4.0 * ih_279[k]
                   + f_0 * lh_552[k];

        t_406[k] = -4.0 * ih_280[k]
                   + f_0 * lh_553[k];

        t_407[k] = -4.0 * ih_281[k]
                   + f_0 * lh_554[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, ih_282, ih_283, ih_284, ih_285, \
                         ih_286, lh_555, lh_556, lh_557, lh_558, \
                         lh_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = -4.0 * ih_282[k]
                   + f_0 * lh_555[k];

        t_409[k] = -4.0 * ih_283[k]
                   + f_0 * lh_556[k];

        t_410[k] = -4.0 * ih_284[k]
                   + f_0 * lh_557[k];

        t_411[k] = -4.0 * ih_285[k]
                   + f_0 * lh_558[k];

        t_412[k] = -4.0 * ih_286[k]
                   + f_0 * lh_559[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, ih_287, ih_288, ih_289, ih_290, \
                         ih_291, lh_560, lh_561, lh_562, lh_563, \
                         lh_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = -4.0 * ih_287[k]
                   + f_0 * lh_560[k];

        t_414[k] = -4.0 * ih_288[k]
                   + f_0 * lh_561[k];

        t_415[k] = -4.0 * ih_289[k]
                   + f_0 * lh_562[k];

        t_416[k] = -4.0 * ih_290[k]
                   + f_0 * lh_563[k];

        t_417[k] = -4.0 * ih_291[k]
                   + f_0 * lh_564[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, t_422, ih_292, ih_293, ih_294, ih_295, \
                         ih_296, lh_565, lh_566, lh_567, lh_568, \
                         lh_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = -4.0 * ih_292[k]
                   + f_0 * lh_565[k];

        t_419[k] = -4.0 * ih_293[k]
                   + f_0 * lh_566[k];

        t_420[k] = -5.0 * ih_294[k]
                   + f_0 * lh_567[k];

        t_421[k] = -5.0 * ih_295[k]
                   + f_0 * lh_568[k];

        t_422[k] = -5.0 * ih_296[k]
                   + f_0 * lh_569[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, ih_297, ih_298, ih_299, ih_300, \
                         ih_301, lh_570, lh_571, lh_572, lh_573, \
                         lh_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = -5.0 * ih_297[k]
                   + f_0 * lh_570[k];

        t_424[k] = -5.0 * ih_298[k]
                   + f_0 * lh_571[k];

        t_425[k] = -5.0 * ih_299[k]
                   + f_0 * lh_572[k];

        t_426[k] = -5.0 * ih_300[k]
                   + f_0 * lh_573[k];

        t_427[k] = -5.0 * ih_301[k]
                   + f_0 * lh_574[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, ih_302, ih_303, ih_304, ih_305, \
                         ih_306, lh_575, lh_576, lh_577, lh_578, \
                         lh_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = -5.0 * ih_302[k]
                   + f_0 * lh_575[k];

        t_429[k] = -5.0 * ih_303[k]
                   + f_0 * lh_576[k];

        t_430[k] = -5.0 * ih_304[k]
                   + f_0 * lh_577[k];

        t_431[k] = -5.0 * ih_305[k]
                   + f_0 * lh_578[k];

        t_432[k] = -5.0 * ih_306[k]
                   + f_0 * lh_579[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, ih_307, ih_308, ih_309, ih_310, \
                         ih_311, lh_580, lh_581, lh_582, lh_583, \
                         lh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = -5.0 * ih_307[k]
                   + f_0 * lh_580[k];

        t_434[k] = -5.0 * ih_308[k]
                   + f_0 * lh_581[k];

        t_435[k] = -5.0 * ih_309[k]
                   + f_0 * lh_582[k];

        t_436[k] = -5.0 * ih_310[k]
                   + f_0 * lh_583[k];

        t_437[k] = -5.0 * ih_311[k]
                   + f_0 * lh_584[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, t_442, t_443, ih_312, ih_313, ih_314, \
                         lh_585, lh_586, lh_587, lh_609, lh_610, \
                         lh_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = -5.0 * ih_312[k]
                   + f_0 * lh_585[k];

        t_439[k] = -5.0 * ih_313[k]
                   + f_0 * lh_586[k];

        t_440[k] = -5.0 * ih_314[k]
                   + f_0 * lh_587[k];

        t_441[k] = f_0 * lh_609[k];

        t_442[k] = f_0 * lh_610[k];

        t_443[k] = f_0 * lh_611[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, t_449, t_450, t_451, lh_612, \
                         lh_613, lh_614, lh_615, lh_616, lh_617, lh_618, \
                         lh_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_0 * lh_612[k];

        t_445[k] = f_0 * lh_613[k];

        t_446[k] = f_0 * lh_614[k];

        t_447[k] = f_0 * lh_615[k];

        t_448[k] = f_0 * lh_616[k];

        t_449[k] = f_0 * lh_617[k];

        t_450[k] = f_0 * lh_618[k];

        t_451[k] = f_0 * lh_619[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, t_457, t_458, t_459, lh_620, \
                         lh_621, lh_622, lh_623, lh_624, lh_625, lh_626, \
                         lh_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_0 * lh_620[k];

        t_453[k] = f_0 * lh_621[k];

        t_454[k] = f_0 * lh_622[k];

        t_455[k] = f_0 * lh_623[k];

        t_456[k] = f_0 * lh_624[k];

        t_457[k] = f_0 * lh_625[k];

        t_458[k] = f_0 * lh_626[k];

        t_459[k] = f_0 * lh_627[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, t_465, ih_315, ih_316, ih_317, \
                         ih_318, lh_628, lh_629, lh_630, lh_631, lh_632, \
                         lh_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_0 * lh_628[k];

        t_461[k] = f_0 * lh_629[k];

        t_462[k] = -ih_315[k]
                   + f_0 * lh_630[k];

        t_463[k] = -ih_316[k]
                   + f_0 * lh_631[k];

        t_464[k] = -ih_317[k]
                   + f_0 * lh_632[k];

        t_465[k] = -ih_318[k]
                   + f_0 * lh_633[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, ih_319, ih_320, ih_321, ih_322, \
                         ih_323, lh_634, lh_635, lh_636, lh_637, \
                         lh_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = -ih_319[k]
                   + f_0 * lh_634[k];

        t_467[k] = -ih_320[k]
                   + f_0 * lh_635[k];

        t_468[k] = -ih_321[k]
                   + f_0 * lh_636[k];

        t_469[k] = -ih_322[k]
                   + f_0 * lh_637[k];

        t_470[k] = -ih_323[k]
                   + f_0 * lh_638[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, t_475, ih_324, ih_325, ih_326, ih_327, \
                         ih_328, lh_639, lh_640, lh_641, lh_642, \
                         lh_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = -ih_324[k]
                   + f_0 * lh_639[k];

        t_472[k] = -ih_325[k]
                   + f_0 * lh_640[k];

        t_473[k] = -ih_326[k]
                   + f_0 * lh_641[k];

        t_474[k] = -ih_327[k]
                   + f_0 * lh_642[k];

        t_475[k] = -ih_328[k]
                   + f_0 * lh_643[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, ih_329, ih_330, ih_331, ih_332, \
                         ih_333, lh_644, lh_645, lh_646, lh_647, \
                         lh_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = -ih_329[k]
                   + f_0 * lh_644[k];

        t_477[k] = -ih_330[k]
                   + f_0 * lh_645[k];

        t_478[k] = -ih_331[k]
                   + f_0 * lh_646[k];

        t_479[k] = -ih_332[k]
                   + f_0 * lh_647[k];

        t_480[k] = -ih_333[k]
                   + f_0 * lh_648[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, ih_334, ih_335, ih_336, ih_337, \
                         ih_338, lh_649, lh_650, lh_651, lh_652, \
                         lh_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = -ih_334[k]
                   + f_0 * lh_649[k];

        t_482[k] = -ih_335[k]
                   + f_0 * lh_650[k];

        t_483[k] = -2.0 * ih_336[k]
                   + f_0 * lh_651[k];

        t_484[k] = -2.0 * ih_337[k]
                   + f_0 * lh_652[k];

        t_485[k] = -2.0 * ih_338[k]
                   + f_0 * lh_653[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, ih_339, ih_340, ih_341, ih_342, \
                         ih_343, lh_654, lh_655, lh_656, lh_657, \
                         lh_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = -2.0 * ih_339[k]
                   + f_0 * lh_654[k];

        t_487[k] = -2.0 * ih_340[k]
                   + f_0 * lh_655[k];

        t_488[k] = -2.0 * ih_341[k]
                   + f_0 * lh_656[k];

        t_489[k] = -2.0 * ih_342[k]
                   + f_0 * lh_657[k];

        t_490[k] = -2.0 * ih_343[k]
                   + f_0 * lh_658[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, ih_344, ih_345, ih_346, ih_347, \
                         ih_348, lh_659, lh_660, lh_661, lh_662, \
                         lh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = -2.0 * ih_344[k]
                   + f_0 * lh_659[k];

        t_492[k] = -2.0 * ih_345[k]
                   + f_0 * lh_660[k];

        t_493[k] = -2.0 * ih_346[k]
                   + f_0 * lh_661[k];

        t_494[k] = -2.0 * ih_347[k]
                   + f_0 * lh_662[k];

        t_495[k] = -2.0 * ih_348[k]
                   + f_0 * lh_663[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ih_349 = buffer.data(ih + 349);
    const auto *ih_350 = buffer.data(ih + 350);
    const auto *ih_351 = buffer.data(ih + 351);
    const auto *ih_352 = buffer.data(ih + 352);
    const auto *ih_353 = buffer.data(ih + 353);
    const auto *ih_354 = buffer.data(ih + 354);
    const auto *ih_355 = buffer.data(ih + 355);
    const auto *ih_356 = buffer.data(ih + 356);
    const auto *ih_357 = buffer.data(ih + 357);
    const auto *ih_358 = buffer.data(ih + 358);
    const auto *ih_359 = buffer.data(ih + 359);
    const auto *ih_360 = buffer.data(ih + 360);
    const auto *ih_361 = buffer.data(ih + 361);
    const auto *ih_362 = buffer.data(ih + 362);
    const auto *ih_363 = buffer.data(ih + 363);
    const auto *ih_364 = buffer.data(ih + 364);
    const auto *ih_365 = buffer.data(ih + 365);
    const auto *ih_366 = buffer.data(ih + 366);
    const auto *ih_367 = buffer.data(ih + 367);
    const auto *ih_368 = buffer.data(ih + 368);
    const auto *ih_369 = buffer.data(ih + 369);
    const auto *ih_370 = buffer.data(ih + 370);
    const auto *ih_371 = buffer.data(ih + 371);
    const auto *ih_372 = buffer.data(ih + 372);
    const auto *ih_373 = buffer.data(ih + 373);
    const auto *ih_374 = buffer.data(ih + 374);
    const auto *ih_375 = buffer.data(ih + 375);
    const auto *ih_376 = buffer.data(ih + 376);
    const auto *ih_377 = buffer.data(ih + 377);
    const auto *ih_378 = buffer.data(ih + 378);
    const auto *ih_379 = buffer.data(ih + 379);
    const auto *ih_380 = buffer.data(ih + 380);
    const auto *ih_381 = buffer.data(ih + 381);
    const auto *ih_382 = buffer.data(ih + 382);
    const auto *ih_383 = buffer.data(ih + 383);
    const auto *ih_384 = buffer.data(ih + 384);
    const auto *ih_385 = buffer.data(ih + 385);
    const auto *ih_386 = buffer.data(ih + 386);
    const auto *ih_387 = buffer.data(ih + 387);
    const auto *ih_388 = buffer.data(ih + 388);
    const auto *ih_389 = buffer.data(ih + 389);
    const auto *ih_390 = buffer.data(ih + 390);
    const auto *ih_391 = buffer.data(ih + 391);
    const auto *ih_392 = buffer.data(ih + 392);
    const auto *ih_393 = buffer.data(ih + 393);
    const auto *ih_394 = buffer.data(ih + 394);
    const auto *ih_395 = buffer.data(ih + 395);
    const auto *ih_396 = buffer.data(ih + 396);
    const auto *ih_397 = buffer.data(ih + 397);
    const auto *ih_398 = buffer.data(ih + 398);
    const auto *ih_399 = buffer.data(ih + 399);
    const auto *ih_400 = buffer.data(ih + 400);
    const auto *ih_401 = buffer.data(ih + 401);
    const auto *ih_402 = buffer.data(ih + 402);
    const auto *ih_403 = buffer.data(ih + 403);
    const auto *ih_404 = buffer.data(ih + 404);
    const auto *ih_405 = buffer.data(ih + 405);
    const auto *ih_406 = buffer.data(ih + 406);
    const auto *ih_407 = buffer.data(ih + 407);
    const auto *ih_408 = buffer.data(ih + 408);
    const auto *ih_409 = buffer.data(ih + 409);
    const auto *ih_410 = buffer.data(ih + 410);
    const auto *ih_411 = buffer.data(ih + 411);
    const auto *ih_412 = buffer.data(ih + 412);
    const auto *ih_413 = buffer.data(ih + 413);
    const auto *ih_414 = buffer.data(ih + 414);
    const auto *ih_415 = buffer.data(ih + 415);
    const auto *ih_416 = buffer.data(ih + 416);
    const auto *ih_417 = buffer.data(ih + 417);
    const auto *ih_418 = buffer.data(ih + 418);
    const auto *ih_419 = buffer.data(ih + 419);
    const auto *ih_420 = buffer.data(ih + 420);
    const auto *ih_421 = buffer.data(ih + 421);
    const auto *ih_422 = buffer.data(ih + 422);
    const auto *ih_423 = buffer.data(ih + 423);
    const auto *ih_424 = buffer.data(ih + 424);
    const auto *ih_425 = buffer.data(ih + 425);
    const auto *ih_426 = buffer.data(ih + 426);
    const auto *ih_427 = buffer.data(ih + 427);
    const auto *ih_428 = buffer.data(ih + 428);
    const auto *ih_429 = buffer.data(ih + 429);
    const auto *ih_430 = buffer.data(ih + 430);
    const auto *ih_431 = buffer.data(ih + 431);
    const auto *ih_432 = buffer.data(ih + 432);
    const auto *ih_433 = buffer.data(ih + 433);
    const auto *ih_434 = buffer.data(ih + 434);
    const auto *ih_435 = buffer.data(ih + 435);
    const auto *ih_436 = buffer.data(ih + 436);
    const auto *ih_437 = buffer.data(ih + 437);
    const auto *ih_438 = buffer.data(ih + 438);
    const auto *ih_439 = buffer.data(ih + 439);
    const auto *ih_440 = buffer.data(ih + 440);
    const auto *ih_441 = buffer.data(ih + 441);
    const auto *ih_442 = buffer.data(ih + 442);
    const auto *ih_443 = buffer.data(ih + 443);
    const auto *ih_444 = buffer.data(ih + 444);
    const auto *ih_445 = buffer.data(ih + 445);
    const auto *ih_446 = buffer.data(ih + 446);
    const auto *ih_447 = buffer.data(ih + 447);
    const auto *ih_448 = buffer.data(ih + 448);
    const auto *ih_449 = buffer.data(ih + 449);
    const auto *ih_450 = buffer.data(ih + 450);
    const auto *ih_451 = buffer.data(ih + 451);
    const auto *ih_452 = buffer.data(ih + 452);
    const auto *ih_453 = buffer.data(ih + 453);
    const auto *ih_454 = buffer.data(ih + 454);
    const auto *ih_455 = buffer.data(ih + 455);
    const auto *ih_456 = buffer.data(ih + 456);
    const auto *ih_457 = buffer.data(ih + 457);
    const auto *ih_458 = buffer.data(ih + 458);
    const auto *ih_459 = buffer.data(ih + 459);
    const auto *ih_460 = buffer.data(ih + 460);
    const auto *ih_461 = buffer.data(ih + 461);
    const auto *ih_462 = buffer.data(ih + 462);
    const auto *ih_463 = buffer.data(ih + 463);
    const auto *ih_464 = buffer.data(ih + 464);
    const auto *ih_465 = buffer.data(ih + 465);
    const auto *ih_466 = buffer.data(ih + 466);
    const auto *ih_467 = buffer.data(ih + 467);
    const auto *ih_468 = buffer.data(ih + 468);
    const auto *ih_469 = buffer.data(ih + 469);
    const auto *ih_470 = buffer.data(ih + 470);
    const auto *ih_471 = buffer.data(ih + 471);
    const auto *ih_472 = buffer.data(ih + 472);
    const auto *ih_473 = buffer.data(ih + 473);
    const auto *ih_474 = buffer.data(ih + 474);
    const auto *ih_475 = buffer.data(ih + 475);
    const auto *ih_476 = buffer.data(ih + 476);
    const auto *ih_477 = buffer.data(ih + 477);
    const auto *ih_478 = buffer.data(ih + 478);
    const auto *ih_479 = buffer.data(ih + 479);
    const auto *ih_480 = buffer.data(ih + 480);
    const auto *ih_481 = buffer.data(ih + 481);
    const auto *ih_482 = buffer.data(ih + 482);
    const auto *ih_483 = buffer.data(ih + 483);
    const auto *ih_484 = buffer.data(ih + 484);
    const auto *ih_485 = buffer.data(ih + 485);

    const auto *lh_664 = buffer.data(lh + 664);
    const auto *lh_665 = buffer.data(lh + 665);
    const auto *lh_666 = buffer.data(lh + 666);
    const auto *lh_667 = buffer.data(lh + 667);
    const auto *lh_668 = buffer.data(lh + 668);
    const auto *lh_669 = buffer.data(lh + 669);
    const auto *lh_670 = buffer.data(lh + 670);
    const auto *lh_671 = buffer.data(lh + 671);
    const auto *lh_672 = buffer.data(lh + 672);
    const auto *lh_673 = buffer.data(lh + 673);
    const auto *lh_674 = buffer.data(lh + 674);
    const auto *lh_675 = buffer.data(lh + 675);
    const auto *lh_676 = buffer.data(lh + 676);
    const auto *lh_677 = buffer.data(lh + 677);
    const auto *lh_678 = buffer.data(lh + 678);
    const auto *lh_679 = buffer.data(lh + 679);
    const auto *lh_680 = buffer.data(lh + 680);
    const auto *lh_681 = buffer.data(lh + 681);
    const auto *lh_682 = buffer.data(lh + 682);
    const auto *lh_683 = buffer.data(lh + 683);
    const auto *lh_684 = buffer.data(lh + 684);
    const auto *lh_685 = buffer.data(lh + 685);
    const auto *lh_686 = buffer.data(lh + 686);
    const auto *lh_687 = buffer.data(lh + 687);
    const auto *lh_688 = buffer.data(lh + 688);
    const auto *lh_689 = buffer.data(lh + 689);
    const auto *lh_690 = buffer.data(lh + 690);
    const auto *lh_691 = buffer.data(lh + 691);
    const auto *lh_692 = buffer.data(lh + 692);
    const auto *lh_693 = buffer.data(lh + 693);
    const auto *lh_694 = buffer.data(lh + 694);
    const auto *lh_695 = buffer.data(lh + 695);
    const auto *lh_696 = buffer.data(lh + 696);
    const auto *lh_697 = buffer.data(lh + 697);
    const auto *lh_698 = buffer.data(lh + 698);
    const auto *lh_699 = buffer.data(lh + 699);
    const auto *lh_700 = buffer.data(lh + 700);
    const auto *lh_701 = buffer.data(lh + 701);
    const auto *lh_702 = buffer.data(lh + 702);
    const auto *lh_703 = buffer.data(lh + 703);
    const auto *lh_704 = buffer.data(lh + 704);
    const auto *lh_705 = buffer.data(lh + 705);
    const auto *lh_706 = buffer.data(lh + 706);
    const auto *lh_707 = buffer.data(lh + 707);
    const auto *lh_708 = buffer.data(lh + 708);
    const auto *lh_709 = buffer.data(lh + 709);
    const auto *lh_710 = buffer.data(lh + 710);
    const auto *lh_711 = buffer.data(lh + 711);
    const auto *lh_712 = buffer.data(lh + 712);
    const auto *lh_713 = buffer.data(lh + 713);
    const auto *lh_714 = buffer.data(lh + 714);
    const auto *lh_715 = buffer.data(lh + 715);
    const auto *lh_716 = buffer.data(lh + 716);
    const auto *lh_717 = buffer.data(lh + 717);
    const auto *lh_718 = buffer.data(lh + 718);
    const auto *lh_719 = buffer.data(lh + 719);
    const auto *lh_720 = buffer.data(lh + 720);
    const auto *lh_721 = buffer.data(lh + 721);
    const auto *lh_722 = buffer.data(lh + 722);
    const auto *lh_723 = buffer.data(lh + 723);
    const auto *lh_724 = buffer.data(lh + 724);
    const auto *lh_725 = buffer.data(lh + 725);
    const auto *lh_726 = buffer.data(lh + 726);
    const auto *lh_727 = buffer.data(lh + 727);
    const auto *lh_728 = buffer.data(lh + 728);
    const auto *lh_729 = buffer.data(lh + 729);
    const auto *lh_730 = buffer.data(lh + 730);
    const auto *lh_731 = buffer.data(lh + 731);
    const auto *lh_732 = buffer.data(lh + 732);
    const auto *lh_733 = buffer.data(lh + 733);
    const auto *lh_734 = buffer.data(lh + 734);
    const auto *lh_735 = buffer.data(lh + 735);
    const auto *lh_736 = buffer.data(lh + 736);
    const auto *lh_737 = buffer.data(lh + 737);
    const auto *lh_738 = buffer.data(lh + 738);
    const auto *lh_739 = buffer.data(lh + 739);
    const auto *lh_740 = buffer.data(lh + 740);
    const auto *lh_741 = buffer.data(lh + 741);
    const auto *lh_742 = buffer.data(lh + 742);
    const auto *lh_743 = buffer.data(lh + 743);
    const auto *lh_744 = buffer.data(lh + 744);
    const auto *lh_745 = buffer.data(lh + 745);
    const auto *lh_746 = buffer.data(lh + 746);
    const auto *lh_747 = buffer.data(lh + 747);
    const auto *lh_748 = buffer.data(lh + 748);
    const auto *lh_749 = buffer.data(lh + 749);
    const auto *lh_750 = buffer.data(lh + 750);
    const auto *lh_751 = buffer.data(lh + 751);
    const auto *lh_752 = buffer.data(lh + 752);
    const auto *lh_753 = buffer.data(lh + 753);
    const auto *lh_754 = buffer.data(lh + 754);
    const auto *lh_755 = buffer.data(lh + 755);
    const auto *lh_777 = buffer.data(lh + 777);
    const auto *lh_778 = buffer.data(lh + 778);
    const auto *lh_779 = buffer.data(lh + 779);
    const auto *lh_780 = buffer.data(lh + 780);
    const auto *lh_781 = buffer.data(lh + 781);
    const auto *lh_782 = buffer.data(lh + 782);
    const auto *lh_783 = buffer.data(lh + 783);
    const auto *lh_784 = buffer.data(lh + 784);
    const auto *lh_785 = buffer.data(lh + 785);
    const auto *lh_786 = buffer.data(lh + 786);
    const auto *lh_787 = buffer.data(lh + 787);
    const auto *lh_788 = buffer.data(lh + 788);
    const auto *lh_789 = buffer.data(lh + 789);
    const auto *lh_790 = buffer.data(lh + 790);
    const auto *lh_791 = buffer.data(lh + 791);
    const auto *lh_792 = buffer.data(lh + 792);
    const auto *lh_793 = buffer.data(lh + 793);
    const auto *lh_794 = buffer.data(lh + 794);
    const auto *lh_795 = buffer.data(lh + 795);
    const auto *lh_796 = buffer.data(lh + 796);
    const auto *lh_797 = buffer.data(lh + 797);
    const auto *lh_798 = buffer.data(lh + 798);
    const auto *lh_799 = buffer.data(lh + 799);
    const auto *lh_800 = buffer.data(lh + 800);
    const auto *lh_801 = buffer.data(lh + 801);
    const auto *lh_802 = buffer.data(lh + 802);
    const auto *lh_803 = buffer.data(lh + 803);
    const auto *lh_804 = buffer.data(lh + 804);
    const auto *lh_805 = buffer.data(lh + 805);
    const auto *lh_806 = buffer.data(lh + 806);
    const auto *lh_807 = buffer.data(lh + 807);
    const auto *lh_808 = buffer.data(lh + 808);
    const auto *lh_809 = buffer.data(lh + 809);
    const auto *lh_810 = buffer.data(lh + 810);
    const auto *lh_811 = buffer.data(lh + 811);
    const auto *lh_812 = buffer.data(lh + 812);
    const auto *lh_813 = buffer.data(lh + 813);
    const auto *lh_814 = buffer.data(lh + 814);
    const auto *lh_815 = buffer.data(lh + 815);
    const auto *lh_816 = buffer.data(lh + 816);
    const auto *lh_817 = buffer.data(lh + 817);
    const auto *lh_818 = buffer.data(lh + 818);
    const auto *lh_819 = buffer.data(lh + 819);
    const auto *lh_820 = buffer.data(lh + 820);
    const auto *lh_821 = buffer.data(lh + 821);
    const auto *lh_822 = buffer.data(lh + 822);
    const auto *lh_823 = buffer.data(lh + 823);
    const auto *lh_824 = buffer.data(lh + 824);
    const auto *lh_825 = buffer.data(lh + 825);
    const auto *lh_826 = buffer.data(lh + 826);
    const auto *lh_827 = buffer.data(lh + 827);
    const auto *lh_828 = buffer.data(lh + 828);
    const auto *lh_829 = buffer.data(lh + 829);
    const auto *lh_830 = buffer.data(lh + 830);
    const auto *lh_831 = buffer.data(lh + 831);
    const auto *lh_832 = buffer.data(lh + 832);
    const auto *lh_833 = buffer.data(lh + 833);
    const auto *lh_834 = buffer.data(lh + 834);
    const auto *lh_835 = buffer.data(lh + 835);
    const auto *lh_836 = buffer.data(lh + 836);
    const auto *lh_837 = buffer.data(lh + 837);
    const auto *lh_838 = buffer.data(lh + 838);
    const auto *lh_839 = buffer.data(lh + 839);
    const auto *lh_840 = buffer.data(lh + 840);
    const auto *lh_841 = buffer.data(lh + 841);
    const auto *lh_842 = buffer.data(lh + 842);

#pragma omp simd aligned(t_496, t_497, t_498, t_499, t_500, ih_349, ih_350, ih_351, ih_352, \
                         ih_353, lh_664, lh_665, lh_666, lh_667, \
                         lh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = -2.0 * ih_349[k]
                   + f_0 * lh_664[k];

        t_497[k] = -2.0 * ih_350[k]
                   + f_0 * lh_665[k];

        t_498[k] = -2.0 * ih_351[k]
                   + f_0 * lh_666[k];

        t_499[k] = -2.0 * ih_352[k]
                   + f_0 * lh_667[k];

        t_500[k] = -2.0 * ih_353[k]
                   + f_0 * lh_668[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, t_505, ih_354, ih_355, ih_356, ih_357, \
                         ih_358, lh_669, lh_670, lh_671, lh_672, \
                         lh_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = -2.0 * ih_354[k]
                   + f_0 * lh_669[k];

        t_502[k] = -2.0 * ih_355[k]
                   + f_0 * lh_670[k];

        t_503[k] = -2.0 * ih_356[k]
                   + f_0 * lh_671[k];

        t_504[k] = -3.0 * ih_357[k]
                   + f_0 * lh_672[k];

        t_505[k] = -3.0 * ih_358[k]
                   + f_0 * lh_673[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, t_509, t_510, ih_359, ih_360, ih_361, ih_362, \
                         ih_363, lh_674, lh_675, lh_676, lh_677, \
                         lh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = -3.0 * ih_359[k]
                   + f_0 * lh_674[k];

        t_507[k] = -3.0 * ih_360[k]
                   + f_0 * lh_675[k];

        t_508[k] = -3.0 * ih_361[k]
                   + f_0 * lh_676[k];

        t_509[k] = -3.0 * ih_362[k]
                   + f_0 * lh_677[k];

        t_510[k] = -3.0 * ih_363[k]
                   + f_0 * lh_678[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, t_515, ih_364, ih_365, ih_366, ih_367, \
                         ih_368, lh_679, lh_680, lh_681, lh_682, \
                         lh_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = -3.0 * ih_364[k]
                   + f_0 * lh_679[k];

        t_512[k] = -3.0 * ih_365[k]
                   + f_0 * lh_680[k];

        t_513[k] = -3.0 * ih_366[k]
                   + f_0 * lh_681[k];

        t_514[k] = -3.0 * ih_367[k]
                   + f_0 * lh_682[k];

        t_515[k] = -3.0 * ih_368[k]
                   + f_0 * lh_683[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, t_520, ih_369, ih_370, ih_371, ih_372, \
                         ih_373, lh_684, lh_685, lh_686, lh_687, \
                         lh_688 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = -3.0 * ih_369[k]
                   + f_0 * lh_684[k];

        t_517[k] = -3.0 * ih_370[k]
                   + f_0 * lh_685[k];

        t_518[k] = -3.0 * ih_371[k]
                   + f_0 * lh_686[k];

        t_519[k] = -3.0 * ih_372[k]
                   + f_0 * lh_687[k];

        t_520[k] = -3.0 * ih_373[k]
                   + f_0 * lh_688[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, t_525, ih_374, ih_375, ih_376, ih_377, \
                         ih_378, lh_689, lh_690, lh_691, lh_692, \
                         lh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = -3.0 * ih_374[k]
                   + f_0 * lh_689[k];

        t_522[k] = -3.0 * ih_375[k]
                   + f_0 * lh_690[k];

        t_523[k] = -3.0 * ih_376[k]
                   + f_0 * lh_691[k];

        t_524[k] = -3.0 * ih_377[k]
                   + f_0 * lh_692[k];

        t_525[k] = -4.0 * ih_378[k]
                   + f_0 * lh_693[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, ih_379, ih_380, ih_381, ih_382, \
                         ih_383, lh_694, lh_695, lh_696, lh_697, \
                         lh_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = -4.0 * ih_379[k]
                   + f_0 * lh_694[k];

        t_527[k] = -4.0 * ih_380[k]
                   + f_0 * lh_695[k];

        t_528[k] = -4.0 * ih_381[k]
                   + f_0 * lh_696[k];

        t_529[k] = -4.0 * ih_382[k]
                   + f_0 * lh_697[k];

        t_530[k] = -4.0 * ih_383[k]
                   + f_0 * lh_698[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, t_535, ih_384, ih_385, ih_386, ih_387, \
                         ih_388, lh_699, lh_700, lh_701, lh_702, \
                         lh_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = -4.0 * ih_384[k]
                   + f_0 * lh_699[k];

        t_532[k] = -4.0 * ih_385[k]
                   + f_0 * lh_700[k];

        t_533[k] = -4.0 * ih_386[k]
                   + f_0 * lh_701[k];

        t_534[k] = -4.0 * ih_387[k]
                   + f_0 * lh_702[k];

        t_535[k] = -4.0 * ih_388[k]
                   + f_0 * lh_703[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, ih_389, ih_390, ih_391, ih_392, \
                         ih_393, lh_704, lh_705, lh_706, lh_707, \
                         lh_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = -4.0 * ih_389[k]
                   + f_0 * lh_704[k];

        t_537[k] = -4.0 * ih_390[k]
                   + f_0 * lh_705[k];

        t_538[k] = -4.0 * ih_391[k]
                   + f_0 * lh_706[k];

        t_539[k] = -4.0 * ih_392[k]
                   + f_0 * lh_707[k];

        t_540[k] = -4.0 * ih_393[k]
                   + f_0 * lh_708[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, t_545, ih_394, ih_395, ih_396, ih_397, \
                         ih_398, lh_709, lh_710, lh_711, lh_712, \
                         lh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = -4.0 * ih_394[k]
                   + f_0 * lh_709[k];

        t_542[k] = -4.0 * ih_395[k]
                   + f_0 * lh_710[k];

        t_543[k] = -4.0 * ih_396[k]
                   + f_0 * lh_711[k];

        t_544[k] = -4.0 * ih_397[k]
                   + f_0 * lh_712[k];

        t_545[k] = -4.0 * ih_398[k]
                   + f_0 * lh_713[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, t_550, ih_399, ih_400, ih_401, ih_402, \
                         ih_403, lh_714, lh_715, lh_716, lh_717, \
                         lh_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = -5.0 * ih_399[k]
                   + f_0 * lh_714[k];

        t_547[k] = -5.0 * ih_400[k]
                   + f_0 * lh_715[k];

        t_548[k] = -5.0 * ih_401[k]
                   + f_0 * lh_716[k];

        t_549[k] = -5.0 * ih_402[k]
                   + f_0 * lh_717[k];

        t_550[k] = -5.0 * ih_403[k]
                   + f_0 * lh_718[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, t_554, t_555, ih_404, ih_405, ih_406, ih_407, \
                         ih_408, lh_719, lh_720, lh_721, lh_722, \
                         lh_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = -5.0 * ih_404[k]
                   + f_0 * lh_719[k];

        t_552[k] = -5.0 * ih_405[k]
                   + f_0 * lh_720[k];

        t_553[k] = -5.0 * ih_406[k]
                   + f_0 * lh_721[k];

        t_554[k] = -5.0 * ih_407[k]
                   + f_0 * lh_722[k];

        t_555[k] = -5.0 * ih_408[k]
                   + f_0 * lh_723[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, ih_409, ih_410, ih_411, ih_412, \
                         ih_413, lh_724, lh_725, lh_726, lh_727, \
                         lh_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = -5.0 * ih_409[k]
                   + f_0 * lh_724[k];

        t_557[k] = -5.0 * ih_410[k]
                   + f_0 * lh_725[k];

        t_558[k] = -5.0 * ih_411[k]
                   + f_0 * lh_726[k];

        t_559[k] = -5.0 * ih_412[k]
                   + f_0 * lh_727[k];

        t_560[k] = -5.0 * ih_413[k]
                   + f_0 * lh_728[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, t_565, ih_414, ih_415, ih_416, ih_417, \
                         ih_418, lh_729, lh_730, lh_731, lh_732, \
                         lh_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = -5.0 * ih_414[k]
                   + f_0 * lh_729[k];

        t_562[k] = -5.0 * ih_415[k]
                   + f_0 * lh_730[k];

        t_563[k] = -5.0 * ih_416[k]
                   + f_0 * lh_731[k];

        t_564[k] = -5.0 * ih_417[k]
                   + f_0 * lh_732[k];

        t_565[k] = -5.0 * ih_418[k]
                   + f_0 * lh_733[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, t_570, ih_419, ih_420, ih_421, ih_422, \
                         ih_423, lh_734, lh_735, lh_736, lh_737, \
                         lh_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = -5.0 * ih_419[k]
                   + f_0 * lh_734[k];

        t_567[k] = -6.0 * ih_420[k]
                   + f_0 * lh_735[k];

        t_568[k] = -6.0 * ih_421[k]
                   + f_0 * lh_736[k];

        t_569[k] = -6.0 * ih_422[k]
                   + f_0 * lh_737[k];

        t_570[k] = -6.0 * ih_423[k]
                   + f_0 * lh_738[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, t_574, t_575, ih_424, ih_425, ih_426, ih_427, \
                         ih_428, lh_739, lh_740, lh_741, lh_742, \
                         lh_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = -6.0 * ih_424[k]
                   + f_0 * lh_739[k];

        t_572[k] = -6.0 * ih_425[k]
                   + f_0 * lh_740[k];

        t_573[k] = -6.0 * ih_426[k]
                   + f_0 * lh_741[k];

        t_574[k] = -6.0 * ih_427[k]
                   + f_0 * lh_742[k];

        t_575[k] = -6.0 * ih_428[k]
                   + f_0 * lh_743[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, ih_429, ih_430, ih_431, ih_432, \
                         ih_433, lh_744, lh_745, lh_746, lh_747, \
                         lh_748 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = -6.0 * ih_429[k]
                   + f_0 * lh_744[k];

        t_577[k] = -6.0 * ih_430[k]
                   + f_0 * lh_745[k];

        t_578[k] = -6.0 * ih_431[k]
                   + f_0 * lh_746[k];

        t_579[k] = -6.0 * ih_432[k]
                   + f_0 * lh_747[k];

        t_580[k] = -6.0 * ih_433[k]
                   + f_0 * lh_748[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, t_584, t_585, ih_434, ih_435, ih_436, ih_437, \
                         ih_438, lh_749, lh_750, lh_751, lh_752, \
                         lh_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = -6.0 * ih_434[k]
                   + f_0 * lh_749[k];

        t_582[k] = -6.0 * ih_435[k]
                   + f_0 * lh_750[k];

        t_583[k] = -6.0 * ih_436[k]
                   + f_0 * lh_751[k];

        t_584[k] = -6.0 * ih_437[k]
                   + f_0 * lh_752[k];

        t_585[k] = -6.0 * ih_438[k]
                   + f_0 * lh_753[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, t_590, t_591, t_592, ih_439, ih_440, \
                         lh_754, lh_755, lh_777, lh_778, lh_779, lh_780, \
                         lh_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = -6.0 * ih_439[k]
                   + f_0 * lh_754[k];

        t_587[k] = -6.0 * ih_440[k]
                   + f_0 * lh_755[k];

        t_588[k] = f_0 * lh_777[k];

        t_589[k] = f_0 * lh_778[k];

        t_590[k] = f_0 * lh_779[k];

        t_591[k] = f_0 * lh_780[k];

        t_592[k] = f_0 * lh_781[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, t_596, t_597, t_598, t_599, t_600, lh_782, \
                         lh_783, lh_784, lh_785, lh_786, lh_787, lh_788, \
                         lh_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = f_0 * lh_782[k];

        t_594[k] = f_0 * lh_783[k];

        t_595[k] = f_0 * lh_784[k];

        t_596[k] = f_0 * lh_785[k];

        t_597[k] = f_0 * lh_786[k];

        t_598[k] = f_0 * lh_787[k];

        t_599[k] = f_0 * lh_788[k];

        t_600[k] = f_0 * lh_789[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, t_604, t_605, t_606, t_607, t_608, lh_790, \
                         lh_791, lh_792, lh_793, lh_794, lh_795, lh_796, \
                         lh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_0 * lh_790[k];

        t_602[k] = f_0 * lh_791[k];

        t_603[k] = f_0 * lh_792[k];

        t_604[k] = f_0 * lh_793[k];

        t_605[k] = f_0 * lh_794[k];

        t_606[k] = f_0 * lh_795[k];

        t_607[k] = f_0 * lh_796[k];

        t_608[k] = f_0 * lh_797[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, t_613, ih_441, ih_442, ih_443, ih_444, \
                         ih_445, lh_798, lh_799, lh_800, lh_801, \
                         lh_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = -ih_441[k]
                   + f_0 * lh_798[k];

        t_610[k] = -ih_442[k]
                   + f_0 * lh_799[k];

        t_611[k] = -ih_443[k]
                   + f_0 * lh_800[k];

        t_612[k] = -ih_444[k]
                   + f_0 * lh_801[k];

        t_613[k] = -ih_445[k]
                   + f_0 * lh_802[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, ih_446, ih_447, ih_448, ih_449, \
                         ih_450, lh_803, lh_804, lh_805, lh_806, \
                         lh_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = -ih_446[k]
                   + f_0 * lh_803[k];

        t_615[k] = -ih_447[k]
                   + f_0 * lh_804[k];

        t_616[k] = -ih_448[k]
                   + f_0 * lh_805[k];

        t_617[k] = -ih_449[k]
                   + f_0 * lh_806[k];

        t_618[k] = -ih_450[k]
                   + f_0 * lh_807[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, ih_451, ih_452, ih_453, ih_454, \
                         ih_455, lh_808, lh_809, lh_810, lh_811, \
                         lh_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = -ih_451[k]
                   + f_0 * lh_808[k];

        t_620[k] = -ih_452[k]
                   + f_0 * lh_809[k];

        t_621[k] = -ih_453[k]
                   + f_0 * lh_810[k];

        t_622[k] = -ih_454[k]
                   + f_0 * lh_811[k];

        t_623[k] = -ih_455[k]
                   + f_0 * lh_812[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, t_628, ih_456, ih_457, ih_458, ih_459, \
                         ih_460, lh_813, lh_814, lh_815, lh_816, \
                         lh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = -ih_456[k]
                   + f_0 * lh_813[k];

        t_625[k] = -ih_457[k]
                   + f_0 * lh_814[k];

        t_626[k] = -ih_458[k]
                   + f_0 * lh_815[k];

        t_627[k] = -ih_459[k]
                   + f_0 * lh_816[k];

        t_628[k] = -ih_460[k]
                   + f_0 * lh_817[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, t_633, ih_461, ih_462, ih_463, ih_464, \
                         ih_465, lh_818, lh_819, lh_820, lh_821, \
                         lh_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = -ih_461[k]
                   + f_0 * lh_818[k];

        t_630[k] = -2.0 * ih_462[k]
                   + f_0 * lh_819[k];

        t_631[k] = -2.0 * ih_463[k]
                   + f_0 * lh_820[k];

        t_632[k] = -2.0 * ih_464[k]
                   + f_0 * lh_821[k];

        t_633[k] = -2.0 * ih_465[k]
                   + f_0 * lh_822[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, ih_466, ih_467, ih_468, ih_469, \
                         ih_470, lh_823, lh_824, lh_825, lh_826, \
                         lh_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = -2.0 * ih_466[k]
                   + f_0 * lh_823[k];

        t_635[k] = -2.0 * ih_467[k]
                   + f_0 * lh_824[k];

        t_636[k] = -2.0 * ih_468[k]
                   + f_0 * lh_825[k];

        t_637[k] = -2.0 * ih_469[k]
                   + f_0 * lh_826[k];

        t_638[k] = -2.0 * ih_470[k]
                   + f_0 * lh_827[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, t_643, ih_471, ih_472, ih_473, ih_474, \
                         ih_475, lh_828, lh_829, lh_830, lh_831, \
                         lh_832 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = -2.0 * ih_471[k]
                   + f_0 * lh_828[k];

        t_640[k] = -2.0 * ih_472[k]
                   + f_0 * lh_829[k];

        t_641[k] = -2.0 * ih_473[k]
                   + f_0 * lh_830[k];

        t_642[k] = -2.0 * ih_474[k]
                   + f_0 * lh_831[k];

        t_643[k] = -2.0 * ih_475[k]
                   + f_0 * lh_832[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, t_648, ih_476, ih_477, ih_478, ih_479, \
                         ih_480, lh_833, lh_834, lh_835, lh_836, \
                         lh_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = -2.0 * ih_476[k]
                   + f_0 * lh_833[k];

        t_645[k] = -2.0 * ih_477[k]
                   + f_0 * lh_834[k];

        t_646[k] = -2.0 * ih_478[k]
                   + f_0 * lh_835[k];

        t_647[k] = -2.0 * ih_479[k]
                   + f_0 * lh_836[k];

        t_648[k] = -2.0 * ih_480[k]
                   + f_0 * lh_837[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, ih_481, ih_482, ih_483, ih_484, \
                         ih_485, lh_838, lh_839, lh_840, lh_841, \
                         lh_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = -2.0 * ih_481[k]
                   + f_0 * lh_838[k];

        t_650[k] = -2.0 * ih_482[k]
                   + f_0 * lh_839[k];

        t_651[k] = -3.0 * ih_483[k]
                   + f_0 * lh_840[k];

        t_652[k] = -3.0 * ih_484[k]
                   + f_0 * lh_841[k];

        t_653[k] = -3.0 * ih_485[k]
                   + f_0 * lh_842[k];
    }
}

static auto
compute_prim_geom_10_kh_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ih, const size_t lh,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ih_486 = buffer.data(ih + 486);
    const auto *ih_487 = buffer.data(ih + 487);
    const auto *ih_488 = buffer.data(ih + 488);
    const auto *ih_489 = buffer.data(ih + 489);
    const auto *ih_490 = buffer.data(ih + 490);
    const auto *ih_491 = buffer.data(ih + 491);
    const auto *ih_492 = buffer.data(ih + 492);
    const auto *ih_493 = buffer.data(ih + 493);
    const auto *ih_494 = buffer.data(ih + 494);
    const auto *ih_495 = buffer.data(ih + 495);
    const auto *ih_496 = buffer.data(ih + 496);
    const auto *ih_497 = buffer.data(ih + 497);
    const auto *ih_498 = buffer.data(ih + 498);
    const auto *ih_499 = buffer.data(ih + 499);
    const auto *ih_500 = buffer.data(ih + 500);
    const auto *ih_501 = buffer.data(ih + 501);
    const auto *ih_502 = buffer.data(ih + 502);
    const auto *ih_503 = buffer.data(ih + 503);
    const auto *ih_504 = buffer.data(ih + 504);
    const auto *ih_505 = buffer.data(ih + 505);
    const auto *ih_506 = buffer.data(ih + 506);
    const auto *ih_507 = buffer.data(ih + 507);
    const auto *ih_508 = buffer.data(ih + 508);
    const auto *ih_509 = buffer.data(ih + 509);
    const auto *ih_510 = buffer.data(ih + 510);
    const auto *ih_511 = buffer.data(ih + 511);
    const auto *ih_512 = buffer.data(ih + 512);
    const auto *ih_513 = buffer.data(ih + 513);
    const auto *ih_514 = buffer.data(ih + 514);
    const auto *ih_515 = buffer.data(ih + 515);
    const auto *ih_516 = buffer.data(ih + 516);
    const auto *ih_517 = buffer.data(ih + 517);
    const auto *ih_518 = buffer.data(ih + 518);
    const auto *ih_519 = buffer.data(ih + 519);
    const auto *ih_520 = buffer.data(ih + 520);
    const auto *ih_521 = buffer.data(ih + 521);
    const auto *ih_522 = buffer.data(ih + 522);
    const auto *ih_523 = buffer.data(ih + 523);
    const auto *ih_524 = buffer.data(ih + 524);
    const auto *ih_525 = buffer.data(ih + 525);
    const auto *ih_526 = buffer.data(ih + 526);
    const auto *ih_527 = buffer.data(ih + 527);
    const auto *ih_528 = buffer.data(ih + 528);
    const auto *ih_529 = buffer.data(ih + 529);
    const auto *ih_530 = buffer.data(ih + 530);
    const auto *ih_531 = buffer.data(ih + 531);
    const auto *ih_532 = buffer.data(ih + 532);
    const auto *ih_533 = buffer.data(ih + 533);
    const auto *ih_534 = buffer.data(ih + 534);
    const auto *ih_535 = buffer.data(ih + 535);
    const auto *ih_536 = buffer.data(ih + 536);
    const auto *ih_537 = buffer.data(ih + 537);
    const auto *ih_538 = buffer.data(ih + 538);
    const auto *ih_539 = buffer.data(ih + 539);
    const auto *ih_540 = buffer.data(ih + 540);
    const auto *ih_541 = buffer.data(ih + 541);
    const auto *ih_542 = buffer.data(ih + 542);
    const auto *ih_543 = buffer.data(ih + 543);
    const auto *ih_544 = buffer.data(ih + 544);
    const auto *ih_545 = buffer.data(ih + 545);
    const auto *ih_546 = buffer.data(ih + 546);
    const auto *ih_547 = buffer.data(ih + 547);
    const auto *ih_548 = buffer.data(ih + 548);
    const auto *ih_549 = buffer.data(ih + 549);
    const auto *ih_550 = buffer.data(ih + 550);
    const auto *ih_551 = buffer.data(ih + 551);
    const auto *ih_552 = buffer.data(ih + 552);
    const auto *ih_553 = buffer.data(ih + 553);
    const auto *ih_554 = buffer.data(ih + 554);
    const auto *ih_555 = buffer.data(ih + 555);
    const auto *ih_556 = buffer.data(ih + 556);
    const auto *ih_557 = buffer.data(ih + 557);
    const auto *ih_558 = buffer.data(ih + 558);
    const auto *ih_559 = buffer.data(ih + 559);
    const auto *ih_560 = buffer.data(ih + 560);
    const auto *ih_561 = buffer.data(ih + 561);
    const auto *ih_562 = buffer.data(ih + 562);
    const auto *ih_563 = buffer.data(ih + 563);
    const auto *ih_564 = buffer.data(ih + 564);
    const auto *ih_565 = buffer.data(ih + 565);
    const auto *ih_566 = buffer.data(ih + 566);
    const auto *ih_567 = buffer.data(ih + 567);
    const auto *ih_568 = buffer.data(ih + 568);
    const auto *ih_569 = buffer.data(ih + 569);
    const auto *ih_570 = buffer.data(ih + 570);
    const auto *ih_571 = buffer.data(ih + 571);
    const auto *ih_572 = buffer.data(ih + 572);
    const auto *ih_573 = buffer.data(ih + 573);
    const auto *ih_574 = buffer.data(ih + 574);
    const auto *ih_575 = buffer.data(ih + 575);
    const auto *ih_576 = buffer.data(ih + 576);
    const auto *ih_577 = buffer.data(ih + 577);
    const auto *ih_578 = buffer.data(ih + 578);
    const auto *ih_579 = buffer.data(ih + 579);
    const auto *ih_580 = buffer.data(ih + 580);
    const auto *ih_581 = buffer.data(ih + 581);
    const auto *ih_582 = buffer.data(ih + 582);
    const auto *ih_583 = buffer.data(ih + 583);
    const auto *ih_584 = buffer.data(ih + 584);
    const auto *ih_585 = buffer.data(ih + 585);
    const auto *ih_586 = buffer.data(ih + 586);
    const auto *ih_587 = buffer.data(ih + 587);

    const auto *lh_843 = buffer.data(lh + 843);
    const auto *lh_844 = buffer.data(lh + 844);
    const auto *lh_845 = buffer.data(lh + 845);
    const auto *lh_846 = buffer.data(lh + 846);
    const auto *lh_847 = buffer.data(lh + 847);
    const auto *lh_848 = buffer.data(lh + 848);
    const auto *lh_849 = buffer.data(lh + 849);
    const auto *lh_850 = buffer.data(lh + 850);
    const auto *lh_851 = buffer.data(lh + 851);
    const auto *lh_852 = buffer.data(lh + 852);
    const auto *lh_853 = buffer.data(lh + 853);
    const auto *lh_854 = buffer.data(lh + 854);
    const auto *lh_855 = buffer.data(lh + 855);
    const auto *lh_856 = buffer.data(lh + 856);
    const auto *lh_857 = buffer.data(lh + 857);
    const auto *lh_858 = buffer.data(lh + 858);
    const auto *lh_859 = buffer.data(lh + 859);
    const auto *lh_860 = buffer.data(lh + 860);
    const auto *lh_861 = buffer.data(lh + 861);
    const auto *lh_862 = buffer.data(lh + 862);
    const auto *lh_863 = buffer.data(lh + 863);
    const auto *lh_864 = buffer.data(lh + 864);
    const auto *lh_865 = buffer.data(lh + 865);
    const auto *lh_866 = buffer.data(lh + 866);
    const auto *lh_867 = buffer.data(lh + 867);
    const auto *lh_868 = buffer.data(lh + 868);
    const auto *lh_869 = buffer.data(lh + 869);
    const auto *lh_870 = buffer.data(lh + 870);
    const auto *lh_871 = buffer.data(lh + 871);
    const auto *lh_872 = buffer.data(lh + 872);
    const auto *lh_873 = buffer.data(lh + 873);
    const auto *lh_874 = buffer.data(lh + 874);
    const auto *lh_875 = buffer.data(lh + 875);
    const auto *lh_876 = buffer.data(lh + 876);
    const auto *lh_877 = buffer.data(lh + 877);
    const auto *lh_878 = buffer.data(lh + 878);
    const auto *lh_879 = buffer.data(lh + 879);
    const auto *lh_880 = buffer.data(lh + 880);
    const auto *lh_881 = buffer.data(lh + 881);
    const auto *lh_882 = buffer.data(lh + 882);
    const auto *lh_883 = buffer.data(lh + 883);
    const auto *lh_884 = buffer.data(lh + 884);
    const auto *lh_885 = buffer.data(lh + 885);
    const auto *lh_886 = buffer.data(lh + 886);
    const auto *lh_887 = buffer.data(lh + 887);
    const auto *lh_888 = buffer.data(lh + 888);
    const auto *lh_889 = buffer.data(lh + 889);
    const auto *lh_890 = buffer.data(lh + 890);
    const auto *lh_891 = buffer.data(lh + 891);
    const auto *lh_892 = buffer.data(lh + 892);
    const auto *lh_893 = buffer.data(lh + 893);
    const auto *lh_894 = buffer.data(lh + 894);
    const auto *lh_895 = buffer.data(lh + 895);
    const auto *lh_896 = buffer.data(lh + 896);
    const auto *lh_897 = buffer.data(lh + 897);
    const auto *lh_898 = buffer.data(lh + 898);
    const auto *lh_899 = buffer.data(lh + 899);
    const auto *lh_900 = buffer.data(lh + 900);
    const auto *lh_901 = buffer.data(lh + 901);
    const auto *lh_902 = buffer.data(lh + 902);
    const auto *lh_903 = buffer.data(lh + 903);
    const auto *lh_904 = buffer.data(lh + 904);
    const auto *lh_905 = buffer.data(lh + 905);
    const auto *lh_906 = buffer.data(lh + 906);
    const auto *lh_907 = buffer.data(lh + 907);
    const auto *lh_908 = buffer.data(lh + 908);
    const auto *lh_909 = buffer.data(lh + 909);
    const auto *lh_910 = buffer.data(lh + 910);
    const auto *lh_911 = buffer.data(lh + 911);
    const auto *lh_912 = buffer.data(lh + 912);
    const auto *lh_913 = buffer.data(lh + 913);
    const auto *lh_914 = buffer.data(lh + 914);
    const auto *lh_915 = buffer.data(lh + 915);
    const auto *lh_916 = buffer.data(lh + 916);
    const auto *lh_917 = buffer.data(lh + 917);
    const auto *lh_918 = buffer.data(lh + 918);
    const auto *lh_919 = buffer.data(lh + 919);
    const auto *lh_920 = buffer.data(lh + 920);
    const auto *lh_921 = buffer.data(lh + 921);
    const auto *lh_922 = buffer.data(lh + 922);
    const auto *lh_923 = buffer.data(lh + 923);
    const auto *lh_924 = buffer.data(lh + 924);
    const auto *lh_925 = buffer.data(lh + 925);
    const auto *lh_926 = buffer.data(lh + 926);
    const auto *lh_927 = buffer.data(lh + 927);
    const auto *lh_928 = buffer.data(lh + 928);
    const auto *lh_929 = buffer.data(lh + 929);
    const auto *lh_930 = buffer.data(lh + 930);
    const auto *lh_931 = buffer.data(lh + 931);
    const auto *lh_932 = buffer.data(lh + 932);
    const auto *lh_933 = buffer.data(lh + 933);
    const auto *lh_934 = buffer.data(lh + 934);
    const auto *lh_935 = buffer.data(lh + 935);
    const auto *lh_936 = buffer.data(lh + 936);
    const auto *lh_937 = buffer.data(lh + 937);
    const auto *lh_938 = buffer.data(lh + 938);
    const auto *lh_939 = buffer.data(lh + 939);
    const auto *lh_940 = buffer.data(lh + 940);
    const auto *lh_941 = buffer.data(lh + 941);
    const auto *lh_942 = buffer.data(lh + 942);
    const auto *lh_943 = buffer.data(lh + 943);
    const auto *lh_944 = buffer.data(lh + 944);

#pragma omp simd aligned(t_654, t_655, t_656, t_657, t_658, ih_486, ih_487, ih_488, ih_489, \
                         ih_490, lh_843, lh_844, lh_845, lh_846, \
                         lh_847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = -3.0 * ih_486[k]
                   + f_0 * lh_843[k];

        t_655[k] = -3.0 * ih_487[k]
                   + f_0 * lh_844[k];

        t_656[k] = -3.0 * ih_488[k]
                   + f_0 * lh_845[k];

        t_657[k] = -3.0 * ih_489[k]
                   + f_0 * lh_846[k];

        t_658[k] = -3.0 * ih_490[k]
                   + f_0 * lh_847[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, ih_491, ih_492, ih_493, ih_494, \
                         ih_495, lh_848, lh_849, lh_850, lh_851, \
                         lh_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = -3.0 * ih_491[k]
                   + f_0 * lh_848[k];

        t_660[k] = -3.0 * ih_492[k]
                   + f_0 * lh_849[k];

        t_661[k] = -3.0 * ih_493[k]
                   + f_0 * lh_850[k];

        t_662[k] = -3.0 * ih_494[k]
                   + f_0 * lh_851[k];

        t_663[k] = -3.0 * ih_495[k]
                   + f_0 * lh_852[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, t_668, ih_496, ih_497, ih_498, ih_499, \
                         ih_500, lh_853, lh_854, lh_855, lh_856, \
                         lh_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = -3.0 * ih_496[k]
                   + f_0 * lh_853[k];

        t_665[k] = -3.0 * ih_497[k]
                   + f_0 * lh_854[k];

        t_666[k] = -3.0 * ih_498[k]
                   + f_0 * lh_855[k];

        t_667[k] = -3.0 * ih_499[k]
                   + f_0 * lh_856[k];

        t_668[k] = -3.0 * ih_500[k]
                   + f_0 * lh_857[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, t_672, t_673, ih_501, ih_502, ih_503, ih_504, \
                         ih_505, lh_858, lh_859, lh_860, lh_861, \
                         lh_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = -3.0 * ih_501[k]
                   + f_0 * lh_858[k];

        t_670[k] = -3.0 * ih_502[k]
                   + f_0 * lh_859[k];

        t_671[k] = -3.0 * ih_503[k]
                   + f_0 * lh_860[k];

        t_672[k] = -4.0 * ih_504[k]
                   + f_0 * lh_861[k];

        t_673[k] = -4.0 * ih_505[k]
                   + f_0 * lh_862[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, t_677, t_678, ih_506, ih_507, ih_508, ih_509, \
                         ih_510, lh_863, lh_864, lh_865, lh_866, \
                         lh_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = -4.0 * ih_506[k]
                   + f_0 * lh_863[k];

        t_675[k] = -4.0 * ih_507[k]
                   + f_0 * lh_864[k];

        t_676[k] = -4.0 * ih_508[k]
                   + f_0 * lh_865[k];

        t_677[k] = -4.0 * ih_509[k]
                   + f_0 * lh_866[k];

        t_678[k] = -4.0 * ih_510[k]
                   + f_0 * lh_867[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, t_682, t_683, ih_511, ih_512, ih_513, ih_514, \
                         ih_515, lh_868, lh_869, lh_870, lh_871, \
                         lh_872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = -4.0 * ih_511[k]
                   + f_0 * lh_868[k];

        t_680[k] = -4.0 * ih_512[k]
                   + f_0 * lh_869[k];

        t_681[k] = -4.0 * ih_513[k]
                   + f_0 * lh_870[k];

        t_682[k] = -4.0 * ih_514[k]
                   + f_0 * lh_871[k];

        t_683[k] = -4.0 * ih_515[k]
                   + f_0 * lh_872[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, t_688, ih_516, ih_517, ih_518, ih_519, \
                         ih_520, lh_873, lh_874, lh_875, lh_876, \
                         lh_877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = -4.0 * ih_516[k]
                   + f_0 * lh_873[k];

        t_685[k] = -4.0 * ih_517[k]
                   + f_0 * lh_874[k];

        t_686[k] = -4.0 * ih_518[k]
                   + f_0 * lh_875[k];

        t_687[k] = -4.0 * ih_519[k]
                   + f_0 * lh_876[k];

        t_688[k] = -4.0 * ih_520[k]
                   + f_0 * lh_877[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, t_693, ih_521, ih_522, ih_523, ih_524, \
                         ih_525, lh_878, lh_879, lh_880, lh_881, \
                         lh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = -4.0 * ih_521[k]
                   + f_0 * lh_878[k];

        t_690[k] = -4.0 * ih_522[k]
                   + f_0 * lh_879[k];

        t_691[k] = -4.0 * ih_523[k]
                   + f_0 * lh_880[k];

        t_692[k] = -4.0 * ih_524[k]
                   + f_0 * lh_881[k];

        t_693[k] = -5.0 * ih_525[k]
                   + f_0 * lh_882[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, t_698, ih_526, ih_527, ih_528, ih_529, \
                         ih_530, lh_883, lh_884, lh_885, lh_886, \
                         lh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = -5.0 * ih_526[k]
                   + f_0 * lh_883[k];

        t_695[k] = -5.0 * ih_527[k]
                   + f_0 * lh_884[k];

        t_696[k] = -5.0 * ih_528[k]
                   + f_0 * lh_885[k];

        t_697[k] = -5.0 * ih_529[k]
                   + f_0 * lh_886[k];

        t_698[k] = -5.0 * ih_530[k]
                   + f_0 * lh_887[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, t_703, ih_531, ih_532, ih_533, ih_534, \
                         ih_535, lh_888, lh_889, lh_890, lh_891, \
                         lh_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = -5.0 * ih_531[k]
                   + f_0 * lh_888[k];

        t_700[k] = -5.0 * ih_532[k]
                   + f_0 * lh_889[k];

        t_701[k] = -5.0 * ih_533[k]
                   + f_0 * lh_890[k];

        t_702[k] = -5.0 * ih_534[k]
                   + f_0 * lh_891[k];

        t_703[k] = -5.0 * ih_535[k]
                   + f_0 * lh_892[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, ih_536, ih_537, ih_538, ih_539, \
                         ih_540, lh_893, lh_894, lh_895, lh_896, \
                         lh_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = -5.0 * ih_536[k]
                   + f_0 * lh_893[k];

        t_705[k] = -5.0 * ih_537[k]
                   + f_0 * lh_894[k];

        t_706[k] = -5.0 * ih_538[k]
                   + f_0 * lh_895[k];

        t_707[k] = -5.0 * ih_539[k]
                   + f_0 * lh_896[k];

        t_708[k] = -5.0 * ih_540[k]
                   + f_0 * lh_897[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, t_712, t_713, ih_541, ih_542, ih_543, ih_544, \
                         ih_545, lh_898, lh_899, lh_900, lh_901, \
                         lh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = -5.0 * ih_541[k]
                   + f_0 * lh_898[k];

        t_710[k] = -5.0 * ih_542[k]
                   + f_0 * lh_899[k];

        t_711[k] = -5.0 * ih_543[k]
                   + f_0 * lh_900[k];

        t_712[k] = -5.0 * ih_544[k]
                   + f_0 * lh_901[k];

        t_713[k] = -5.0 * ih_545[k]
                   + f_0 * lh_902[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, t_718, ih_546, ih_547, ih_548, ih_549, \
                         ih_550, lh_903, lh_904, lh_905, lh_906, \
                         lh_907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = -6.0 * ih_546[k]
                   + f_0 * lh_903[k];

        t_715[k] = -6.0 * ih_547[k]
                   + f_0 * lh_904[k];

        t_716[k] = -6.0 * ih_548[k]
                   + f_0 * lh_905[k];

        t_717[k] = -6.0 * ih_549[k]
                   + f_0 * lh_906[k];

        t_718[k] = -6.0 * ih_550[k]
                   + f_0 * lh_907[k];
    }

#pragma omp simd aligned(t_719, t_720, t_721, t_722, t_723, ih_551, ih_552, ih_553, ih_554, \
                         ih_555, lh_908, lh_909, lh_910, lh_911, \
                         lh_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_719[k] = -6.0 * ih_551[k]
                   + f_0 * lh_908[k];

        t_720[k] = -6.0 * ih_552[k]
                   + f_0 * lh_909[k];

        t_721[k] = -6.0 * ih_553[k]
                   + f_0 * lh_910[k];

        t_722[k] = -6.0 * ih_554[k]
                   + f_0 * lh_911[k];

        t_723[k] = -6.0 * ih_555[k]
                   + f_0 * lh_912[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, ih_556, ih_557, ih_558, ih_559, \
                         ih_560, lh_913, lh_914, lh_915, lh_916, \
                         lh_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = -6.0 * ih_556[k]
                   + f_0 * lh_913[k];

        t_725[k] = -6.0 * ih_557[k]
                   + f_0 * lh_914[k];

        t_726[k] = -6.0 * ih_558[k]
                   + f_0 * lh_915[k];

        t_727[k] = -6.0 * ih_559[k]
                   + f_0 * lh_916[k];

        t_728[k] = -6.0 * ih_560[k]
                   + f_0 * lh_917[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, t_733, ih_561, ih_562, ih_563, ih_564, \
                         ih_565, lh_918, lh_919, lh_920, lh_921, \
                         lh_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = -6.0 * ih_561[k]
                   + f_0 * lh_918[k];

        t_730[k] = -6.0 * ih_562[k]
                   + f_0 * lh_919[k];

        t_731[k] = -6.0 * ih_563[k]
                   + f_0 * lh_920[k];

        t_732[k] = -6.0 * ih_564[k]
                   + f_0 * lh_921[k];

        t_733[k] = -6.0 * ih_565[k]
                   + f_0 * lh_922[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, t_738, ih_566, ih_567, ih_568, ih_569, \
                         ih_570, lh_923, lh_924, lh_925, lh_926, \
                         lh_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = -6.0 * ih_566[k]
                   + f_0 * lh_923[k];

        t_735[k] = -7.0 * ih_567[k]
                   + f_0 * lh_924[k];

        t_736[k] = -7.0 * ih_568[k]
                   + f_0 * lh_925[k];

        t_737[k] = -7.0 * ih_569[k]
                   + f_0 * lh_926[k];

        t_738[k] = -7.0 * ih_570[k]
                   + f_0 * lh_927[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, t_743, ih_571, ih_572, ih_573, ih_574, \
                         ih_575, lh_928, lh_929, lh_930, lh_931, \
                         lh_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = -7.0 * ih_571[k]
                   + f_0 * lh_928[k];

        t_740[k] = -7.0 * ih_572[k]
                   + f_0 * lh_929[k];

        t_741[k] = -7.0 * ih_573[k]
                   + f_0 * lh_930[k];

        t_742[k] = -7.0 * ih_574[k]
                   + f_0 * lh_931[k];

        t_743[k] = -7.0 * ih_575[k]
                   + f_0 * lh_932[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, t_748, ih_576, ih_577, ih_578, ih_579, \
                         ih_580, lh_933, lh_934, lh_935, lh_936, \
                         lh_937 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = -7.0 * ih_576[k]
                   + f_0 * lh_933[k];

        t_745[k] = -7.0 * ih_577[k]
                   + f_0 * lh_934[k];

        t_746[k] = -7.0 * ih_578[k]
                   + f_0 * lh_935[k];

        t_747[k] = -7.0 * ih_579[k]
                   + f_0 * lh_936[k];

        t_748[k] = -7.0 * ih_580[k]
                   + f_0 * lh_937[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, t_752, t_753, ih_581, ih_582, ih_583, ih_584, \
                         ih_585, lh_938, lh_939, lh_940, lh_941, \
                         lh_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = -7.0 * ih_581[k]
                   + f_0 * lh_938[k];

        t_750[k] = -7.0 * ih_582[k]
                   + f_0 * lh_939[k];

        t_751[k] = -7.0 * ih_583[k]
                   + f_0 * lh_940[k];

        t_752[k] = -7.0 * ih_584[k]
                   + f_0 * lh_941[k];

        t_753[k] = -7.0 * ih_585[k]
                   + f_0 * lh_942[k];
    }

#pragma omp simd aligned(t_754, t_755, ih_586, ih_587, lh_943, lh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = -7.0 * ih_586[k]
                   + f_0 * lh_943[k];

        t_755[k] = -7.0 * ih_587[k]
                   + f_0 * lh_944[k];
    }
}

auto
compute_prim_geom_10_kh_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t ih, const size_t lh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kh_electron_repulsion_2_piece0(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_2_piece1(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_2_piece2(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_2_piece3(buffer, target, ih, lh, ncols, alpha);

    compute_prim_geom_10_kh_electron_repulsion_2_piece4(buffer, target, ih, lh, ncols, alpha);
}

}  // namespace simdt2ceri
