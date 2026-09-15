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


#include "SimdElectronRepulsionGeom10VrrRecFG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_fg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t dg, const size_t gg,
                                             const size_t ncols, const double alpha) -> void
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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_52 = buffer.data(dg + 52);
    const auto *dg_53 = buffer.data(dg + 53);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_60 = buffer.data(dg + 60);
    const auto *dg_61 = buffer.data(dg + 61);
    const auto *dg_62 = buffer.data(dg + 62);
    const auto *dg_63 = buffer.data(dg + 63);
    const auto *dg_64 = buffer.data(dg + 64);
    const auto *dg_65 = buffer.data(dg + 65);
    const auto *dg_66 = buffer.data(dg + 66);
    const auto *dg_67 = buffer.data(dg + 67);
    const auto *dg_68 = buffer.data(dg + 68);
    const auto *dg_69 = buffer.data(dg + 69);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_76 = buffer.data(dg + 76);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_79 = buffer.data(dg + 79);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_81 = buffer.data(dg + 81);
    const auto *dg_82 = buffer.data(dg + 82);
    const auto *dg_83 = buffer.data(dg + 83);
    const auto *dg_84 = buffer.data(dg + 84);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_88 = buffer.data(dg + 88);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_94 = buffer.data(gg + 94);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_97 = buffer.data(gg + 97);
    const auto *gg_98 = buffer.data(gg + 98);
    const auto *gg_99 = buffer.data(gg + 99);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_101 = buffer.data(gg + 101);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_106 = buffer.data(gg + 106);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_109 = buffer.data(gg + 109);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_113 = buffer.data(gg + 113);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_124 = buffer.data(gg + 124);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_128 = buffer.data(gg + 128);
    const auto *gg_129 = buffer.data(gg + 129);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_134 = buffer.data(gg + 134);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_136 = buffer.data(gg + 136);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_138 = buffer.data(gg + 138);
    const auto *gg_139 = buffer.data(gg + 139);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_141 = buffer.data(gg + 141);
    const auto *gg_142 = buffer.data(gg + 142);
    const auto *gg_143 = buffer.data(gg + 143);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_145 = buffer.data(gg + 145);
    const auto *gg_146 = buffer.data(gg + 146);
    const auto *gg_147 = buffer.data(gg + 147);
    const auto *gg_148 = buffer.data(gg + 148);
    const auto *gg_149 = buffer.data(gg + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, dg_0, dg_1, dg_2, dg_3, dg_4, gg_0, gg_1, \
                         gg_2, gg_3, gg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -3.0 * dg_0[k]
                 + f_0 * gg_0[k];

        t_1[k] = -3.0 * dg_1[k]
                 + f_0 * gg_1[k];

        t_2[k] = -3.0 * dg_2[k]
                 + f_0 * gg_2[k];

        t_3[k] = -3.0 * dg_3[k]
                 + f_0 * gg_3[k];

        t_4[k] = -3.0 * dg_4[k]
                 + f_0 * gg_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, dg_5, dg_6, dg_7, dg_8, dg_9, gg_5, gg_6, \
                         gg_7, gg_8, gg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -3.0 * dg_5[k]
                 + f_0 * gg_5[k];

        t_6[k] = -3.0 * dg_6[k]
                 + f_0 * gg_6[k];

        t_7[k] = -3.0 * dg_7[k]
                 + f_0 * gg_7[k];

        t_8[k] = -3.0 * dg_8[k]
                 + f_0 * gg_8[k];

        t_9[k] = -3.0 * dg_9[k]
                 + f_0 * gg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, dg_10, dg_11, dg_12, dg_13, dg_14, \
                         gg_10, gg_11, gg_12, gg_13, gg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -3.0 * dg_10[k]
                  + f_0 * gg_10[k];

        t_11[k] = -3.0 * dg_11[k]
                  + f_0 * gg_11[k];

        t_12[k] = -3.0 * dg_12[k]
                  + f_0 * gg_12[k];

        t_13[k] = -3.0 * dg_13[k]
                  + f_0 * gg_13[k];

        t_14[k] = -3.0 * dg_14[k]
                  + f_0 * gg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, dg_15, dg_16, dg_17, dg_18, dg_19, \
                         gg_15, gg_16, gg_17, gg_18, gg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -2.0 * dg_15[k]
                  + f_0 * gg_15[k];

        t_16[k] = -2.0 * dg_16[k]
                  + f_0 * gg_16[k];

        t_17[k] = -2.0 * dg_17[k]
                  + f_0 * gg_17[k];

        t_18[k] = -2.0 * dg_18[k]
                  + f_0 * gg_18[k];

        t_19[k] = -2.0 * dg_19[k]
                  + f_0 * gg_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, dg_20, dg_21, dg_22, dg_23, dg_24, \
                         gg_20, gg_21, gg_22, gg_23, gg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -2.0 * dg_20[k]
                  + f_0 * gg_20[k];

        t_21[k] = -2.0 * dg_21[k]
                  + f_0 * gg_21[k];

        t_22[k] = -2.0 * dg_22[k]
                  + f_0 * gg_22[k];

        t_23[k] = -2.0 * dg_23[k]
                  + f_0 * gg_23[k];

        t_24[k] = -2.0 * dg_24[k]
                  + f_0 * gg_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, dg_25, dg_26, dg_27, dg_28, dg_29, \
                         gg_25, gg_26, gg_27, gg_28, gg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -2.0 * dg_25[k]
                  + f_0 * gg_25[k];

        t_26[k] = -2.0 * dg_26[k]
                  + f_0 * gg_26[k];

        t_27[k] = -2.0 * dg_27[k]
                  + f_0 * gg_27[k];

        t_28[k] = -2.0 * dg_28[k]
                  + f_0 * gg_28[k];

        t_29[k] = -2.0 * dg_29[k]
                  + f_0 * gg_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, dg_30, dg_31, dg_32, dg_33, dg_34, \
                         gg_30, gg_31, gg_32, gg_33, gg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -2.0 * dg_30[k]
                  + f_0 * gg_30[k];

        t_31[k] = -2.0 * dg_31[k]
                  + f_0 * gg_31[k];

        t_32[k] = -2.0 * dg_32[k]
                  + f_0 * gg_32[k];

        t_33[k] = -2.0 * dg_33[k]
                  + f_0 * gg_33[k];

        t_34[k] = -2.0 * dg_34[k]
                  + f_0 * gg_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, dg_35, dg_36, dg_37, dg_38, dg_39, \
                         gg_35, gg_36, gg_37, gg_38, gg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -2.0 * dg_35[k]
                  + f_0 * gg_35[k];

        t_36[k] = -2.0 * dg_36[k]
                  + f_0 * gg_36[k];

        t_37[k] = -2.0 * dg_37[k]
                  + f_0 * gg_37[k];

        t_38[k] = -2.0 * dg_38[k]
                  + f_0 * gg_38[k];

        t_39[k] = -2.0 * dg_39[k]
                  + f_0 * gg_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, dg_40, dg_41, dg_42, dg_43, dg_44, \
                         gg_40, gg_41, gg_42, gg_43, gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -2.0 * dg_40[k]
                  + f_0 * gg_40[k];

        t_41[k] = -2.0 * dg_41[k]
                  + f_0 * gg_41[k];

        t_42[k] = -2.0 * dg_42[k]
                  + f_0 * gg_42[k];

        t_43[k] = -2.0 * dg_43[k]
                  + f_0 * gg_43[k];

        t_44[k] = -2.0 * dg_44[k]
                  + f_0 * gg_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, dg_45, dg_46, dg_47, dg_48, dg_49, \
                         gg_45, gg_46, gg_47, gg_48, gg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -dg_45[k]
                  + f_0 * gg_45[k];

        t_46[k] = -dg_46[k]
                  + f_0 * gg_46[k];

        t_47[k] = -dg_47[k]
                  + f_0 * gg_47[k];

        t_48[k] = -dg_48[k]
                  + f_0 * gg_48[k];

        t_49[k] = -dg_49[k]
                  + f_0 * gg_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, dg_50, dg_51, dg_52, dg_53, dg_54, \
                         gg_50, gg_51, gg_52, gg_53, gg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -dg_50[k]
                  + f_0 * gg_50[k];

        t_51[k] = -dg_51[k]
                  + f_0 * gg_51[k];

        t_52[k] = -dg_52[k]
                  + f_0 * gg_52[k];

        t_53[k] = -dg_53[k]
                  + f_0 * gg_53[k];

        t_54[k] = -dg_54[k]
                  + f_0 * gg_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, dg_55, dg_56, dg_57, dg_58, dg_59, \
                         gg_55, gg_56, gg_57, gg_58, gg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -dg_55[k]
                  + f_0 * gg_55[k];

        t_56[k] = -dg_56[k]
                  + f_0 * gg_56[k];

        t_57[k] = -dg_57[k]
                  + f_0 * gg_57[k];

        t_58[k] = -dg_58[k]
                  + f_0 * gg_58[k];

        t_59[k] = -dg_59[k]
                  + f_0 * gg_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, dg_60, dg_61, dg_62, dg_63, dg_64, \
                         gg_60, gg_61, gg_62, gg_63, gg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -dg_60[k]
                  + f_0 * gg_60[k];

        t_61[k] = -dg_61[k]
                  + f_0 * gg_61[k];

        t_62[k] = -dg_62[k]
                  + f_0 * gg_62[k];

        t_63[k] = -dg_63[k]
                  + f_0 * gg_63[k];

        t_64[k] = -dg_64[k]
                  + f_0 * gg_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, dg_65, dg_66, dg_67, dg_68, dg_69, \
                         gg_65, gg_66, gg_67, gg_68, gg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -dg_65[k]
                  + f_0 * gg_65[k];

        t_66[k] = -dg_66[k]
                  + f_0 * gg_66[k];

        t_67[k] = -dg_67[k]
                  + f_0 * gg_67[k];

        t_68[k] = -dg_68[k]
                  + f_0 * gg_68[k];

        t_69[k] = -dg_69[k]
                  + f_0 * gg_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, dg_70, dg_71, dg_72, dg_73, dg_74, \
                         gg_70, gg_71, gg_72, gg_73, gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -dg_70[k]
                  + f_0 * gg_70[k];

        t_71[k] = -dg_71[k]
                  + f_0 * gg_71[k];

        t_72[k] = -dg_72[k]
                  + f_0 * gg_72[k];

        t_73[k] = -dg_73[k]
                  + f_0 * gg_73[k];

        t_74[k] = -dg_74[k]
                  + f_0 * gg_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, dg_75, dg_76, dg_77, dg_78, dg_79, \
                         gg_75, gg_76, gg_77, gg_78, gg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -dg_75[k]
                  + f_0 * gg_75[k];

        t_76[k] = -dg_76[k]
                  + f_0 * gg_76[k];

        t_77[k] = -dg_77[k]
                  + f_0 * gg_77[k];

        t_78[k] = -dg_78[k]
                  + f_0 * gg_78[k];

        t_79[k] = -dg_79[k]
                  + f_0 * gg_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, dg_80, dg_81, dg_82, dg_83, dg_84, \
                         gg_80, gg_81, gg_82, gg_83, gg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -dg_80[k]
                  + f_0 * gg_80[k];

        t_81[k] = -dg_81[k]
                  + f_0 * gg_81[k];

        t_82[k] = -dg_82[k]
                  + f_0 * gg_82[k];

        t_83[k] = -dg_83[k]
                  + f_0 * gg_83[k];

        t_84[k] = -dg_84[k]
                  + f_0 * gg_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, dg_85, dg_86, dg_87, dg_88, dg_89, \
                         gg_85, gg_86, gg_87, gg_88, gg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -dg_85[k]
                  + f_0 * gg_85[k];

        t_86[k] = -dg_86[k]
                  + f_0 * gg_86[k];

        t_87[k] = -dg_87[k]
                  + f_0 * gg_87[k];

        t_88[k] = -dg_88[k]
                  + f_0 * gg_88[k];

        t_89[k] = -dg_89[k]
                  + f_0 * gg_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, gg_90, gg_91, gg_92, \
                         gg_93, gg_94, gg_95, gg_96, gg_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * gg_90[k];

        t_91[k] = f_0 * gg_91[k];

        t_92[k] = f_0 * gg_92[k];

        t_93[k] = f_0 * gg_93[k];

        t_94[k] = f_0 * gg_94[k];

        t_95[k] = f_0 * gg_95[k];

        t_96[k] = f_0 * gg_96[k];

        t_97[k] = f_0 * gg_97[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, t_104, t_105, gg_98, gg_99, \
                         gg_100, gg_101, gg_102, gg_103, gg_104, \
                         gg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * gg_98[k];

        t_99[k] = f_0 * gg_99[k];

        t_100[k] = f_0 * gg_100[k];

        t_101[k] = f_0 * gg_101[k];

        t_102[k] = f_0 * gg_102[k];

        t_103[k] = f_0 * gg_103[k];

        t_104[k] = f_0 * gg_104[k];

        t_105[k] = f_0 * gg_105[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, t_112, t_113, gg_106, \
                         gg_107, gg_108, gg_109, gg_110, gg_111, gg_112, \
                         gg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * gg_106[k];

        t_107[k] = f_0 * gg_107[k];

        t_108[k] = f_0 * gg_108[k];

        t_109[k] = f_0 * gg_109[k];

        t_110[k] = f_0 * gg_110[k];

        t_111[k] = f_0 * gg_111[k];

        t_112[k] = f_0 * gg_112[k];

        t_113[k] = f_0 * gg_113[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, t_120, t_121, gg_114, \
                         gg_115, gg_116, gg_117, gg_118, gg_119, gg_120, \
                         gg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * gg_114[k];

        t_115[k] = f_0 * gg_115[k];

        t_116[k] = f_0 * gg_116[k];

        t_117[k] = f_0 * gg_117[k];

        t_118[k] = f_0 * gg_118[k];

        t_119[k] = f_0 * gg_119[k];

        t_120[k] = f_0 * gg_120[k];

        t_121[k] = f_0 * gg_121[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, t_128, t_129, gg_122, \
                         gg_123, gg_124, gg_125, gg_126, gg_127, gg_128, \
                         gg_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * gg_122[k];

        t_123[k] = f_0 * gg_123[k];

        t_124[k] = f_0 * gg_124[k];

        t_125[k] = f_0 * gg_125[k];

        t_126[k] = f_0 * gg_126[k];

        t_127[k] = f_0 * gg_127[k];

        t_128[k] = f_0 * gg_128[k];

        t_129[k] = f_0 * gg_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, t_135, t_136, t_137, gg_130, \
                         gg_131, gg_132, gg_133, gg_134, gg_135, gg_136, \
                         gg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * gg_130[k];

        t_131[k] = f_0 * gg_131[k];

        t_132[k] = f_0 * gg_132[k];

        t_133[k] = f_0 * gg_133[k];

        t_134[k] = f_0 * gg_134[k];

        t_135[k] = f_0 * gg_135[k];

        t_136[k] = f_0 * gg_136[k];

        t_137[k] = f_0 * gg_137[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, t_144, t_145, gg_138, \
                         gg_139, gg_140, gg_141, gg_142, gg_143, gg_144, \
                         gg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * gg_138[k];

        t_139[k] = f_0 * gg_139[k];

        t_140[k] = f_0 * gg_140[k];

        t_141[k] = f_0 * gg_141[k];

        t_142[k] = f_0 * gg_142[k];

        t_143[k] = f_0 * gg_143[k];

        t_144[k] = f_0 * gg_144[k];

        t_145[k] = f_0 * gg_145[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, gg_146, gg_147, gg_148, \
                         gg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_0 * gg_146[k];

        t_147[k] = f_0 * gg_147[k];

        t_148[k] = f_0 * gg_148[k];

        t_149[k] = f_0 * gg_149[k];
    }
}

auto
compute_prim_geom_10_fg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t dg, const size_t gg,
                                             const size_t ncols, const double alpha) -> void
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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_52 = buffer.data(dg + 52);
    const auto *dg_53 = buffer.data(dg + 53);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_60 = buffer.data(dg + 60);
    const auto *dg_61 = buffer.data(dg + 61);
    const auto *dg_62 = buffer.data(dg + 62);
    const auto *dg_63 = buffer.data(dg + 63);
    const auto *dg_64 = buffer.data(dg + 64);
    const auto *dg_65 = buffer.data(dg + 65);
    const auto *dg_66 = buffer.data(dg + 66);
    const auto *dg_67 = buffer.data(dg + 67);
    const auto *dg_68 = buffer.data(dg + 68);
    const auto *dg_69 = buffer.data(dg + 69);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_76 = buffer.data(dg + 76);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_79 = buffer.data(dg + 79);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_81 = buffer.data(dg + 81);
    const auto *dg_82 = buffer.data(dg + 82);
    const auto *dg_83 = buffer.data(dg + 83);
    const auto *dg_84 = buffer.data(dg + 84);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_88 = buffer.data(dg + 88);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_94 = buffer.data(gg + 94);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_97 = buffer.data(gg + 97);
    const auto *gg_98 = buffer.data(gg + 98);
    const auto *gg_99 = buffer.data(gg + 99);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_101 = buffer.data(gg + 101);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_106 = buffer.data(gg + 106);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_109 = buffer.data(gg + 109);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_113 = buffer.data(gg + 113);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_124 = buffer.data(gg + 124);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_128 = buffer.data(gg + 128);
    const auto *gg_129 = buffer.data(gg + 129);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_134 = buffer.data(gg + 134);
    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_151 = buffer.data(gg + 151);
    const auto *gg_152 = buffer.data(gg + 152);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_154 = buffer.data(gg + 154);
    const auto *gg_155 = buffer.data(gg + 155);
    const auto *gg_156 = buffer.data(gg + 156);
    const auto *gg_157 = buffer.data(gg + 157);
    const auto *gg_158 = buffer.data(gg + 158);
    const auto *gg_159 = buffer.data(gg + 159);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_161 = buffer.data(gg + 161);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_163 = buffer.data(gg + 163);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_165 = buffer.data(gg + 165);
    const auto *gg_166 = buffer.data(gg + 166);
    const auto *gg_167 = buffer.data(gg + 167);
    const auto *gg_168 = buffer.data(gg + 168);
    const auto *gg_169 = buffer.data(gg + 169);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_171 = buffer.data(gg + 171);
    const auto *gg_172 = buffer.data(gg + 172);
    const auto *gg_173 = buffer.data(gg + 173);
    const auto *gg_174 = buffer.data(gg + 174);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_181 = buffer.data(gg + 181);
    const auto *gg_182 = buffer.data(gg + 182);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_184 = buffer.data(gg + 184);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_186 = buffer.data(gg + 186);
    const auto *gg_187 = buffer.data(gg + 187);
    const auto *gg_188 = buffer.data(gg + 188);
    const auto *gg_189 = buffer.data(gg + 189);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_195 = buffer.data(gg + 195);
    const auto *gg_196 = buffer.data(gg + 196);
    const auto *gg_197 = buffer.data(gg + 197);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_199 = buffer.data(gg + 199);
    const auto *gg_200 = buffer.data(gg + 200);
    const auto *gg_201 = buffer.data(gg + 201);
    const auto *gg_202 = buffer.data(gg + 202);
    const auto *gg_203 = buffer.data(gg + 203);
    const auto *gg_204 = buffer.data(gg + 204);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gg_15, gg_16, gg_17, gg_18, \
                         gg_19, gg_20, gg_21, gg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gg_15[k];

        t_1[k] = f_0 * gg_16[k];

        t_2[k] = f_0 * gg_17[k];

        t_3[k] = f_0 * gg_18[k];

        t_4[k] = f_0 * gg_19[k];

        t_5[k] = f_0 * gg_20[k];

        t_6[k] = f_0 * gg_21[k];

        t_7[k] = f_0 * gg_22[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, gg_23, gg_24, gg_25, gg_26, \
                         gg_27, gg_28, gg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gg_23[k];

        t_9[k] = f_0 * gg_24[k];

        t_10[k] = f_0 * gg_25[k];

        t_11[k] = f_0 * gg_26[k];

        t_12[k] = f_0 * gg_27[k];

        t_13[k] = f_0 * gg_28[k];

        t_14[k] = f_0 * gg_29[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, dg_0, dg_1, dg_2, dg_3, dg_4, gg_45, \
                         gg_46, gg_47, gg_48, gg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -dg_0[k]
                  + f_0 * gg_45[k];

        t_16[k] = -dg_1[k]
                  + f_0 * gg_46[k];

        t_17[k] = -dg_2[k]
                  + f_0 * gg_47[k];

        t_18[k] = -dg_3[k]
                  + f_0 * gg_48[k];

        t_19[k] = -dg_4[k]
                  + f_0 * gg_49[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, dg_5, dg_6, dg_7, dg_8, dg_9, gg_50, \
                         gg_51, gg_52, gg_53, gg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -dg_5[k]
                  + f_0 * gg_50[k];

        t_21[k] = -dg_6[k]
                  + f_0 * gg_51[k];

        t_22[k] = -dg_7[k]
                  + f_0 * gg_52[k];

        t_23[k] = -dg_8[k]
                  + f_0 * gg_53[k];

        t_24[k] = -dg_9[k]
                  + f_0 * gg_54[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, dg_10, dg_11, dg_12, dg_13, dg_14, \
                         gg_55, gg_56, gg_57, gg_58, gg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -dg_10[k]
                  + f_0 * gg_55[k];

        t_26[k] = -dg_11[k]
                  + f_0 * gg_56[k];

        t_27[k] = -dg_12[k]
                  + f_0 * gg_57[k];

        t_28[k] = -dg_13[k]
                  + f_0 * gg_58[k];

        t_29[k] = -dg_14[k]
                  + f_0 * gg_59[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, t_36, t_37, gg_60, gg_61, gg_62, \
                         gg_63, gg_64, gg_65, gg_66, gg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * gg_60[k];

        t_31[k] = f_0 * gg_61[k];

        t_32[k] = f_0 * gg_62[k];

        t_33[k] = f_0 * gg_63[k];

        t_34[k] = f_0 * gg_64[k];

        t_35[k] = f_0 * gg_65[k];

        t_36[k] = f_0 * gg_66[k];

        t_37[k] = f_0 * gg_67[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, gg_68, gg_69, gg_70, gg_71, \
                         gg_72, gg_73, gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * gg_68[k];

        t_39[k] = f_0 * gg_69[k];

        t_40[k] = f_0 * gg_70[k];

        t_41[k] = f_0 * gg_71[k];

        t_42[k] = f_0 * gg_72[k];

        t_43[k] = f_0 * gg_73[k];

        t_44[k] = f_0 * gg_74[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, dg_15, dg_16, dg_17, dg_18, dg_19, \
                         gg_90, gg_91, gg_92, gg_93, gg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * dg_15[k]
                  + f_0 * gg_90[k];

        t_46[k] = -2.0 * dg_16[k]
                  + f_0 * gg_91[k];

        t_47[k] = -2.0 * dg_17[k]
                  + f_0 * gg_92[k];

        t_48[k] = -2.0 * dg_18[k]
                  + f_0 * gg_93[k];

        t_49[k] = -2.0 * dg_19[k]
                  + f_0 * gg_94[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, dg_20, dg_21, dg_22, dg_23, dg_24, \
                         gg_95, gg_96, gg_97, gg_98, gg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * dg_20[k]
                  + f_0 * gg_95[k];

        t_51[k] = -2.0 * dg_21[k]
                  + f_0 * gg_96[k];

        t_52[k] = -2.0 * dg_22[k]
                  + f_0 * gg_97[k];

        t_53[k] = -2.0 * dg_23[k]
                  + f_0 * gg_98[k];

        t_54[k] = -2.0 * dg_24[k]
                  + f_0 * gg_99[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, dg_25, dg_26, dg_27, dg_28, dg_29, \
                         gg_100, gg_101, gg_102, gg_103, gg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * dg_25[k]
                  + f_0 * gg_100[k];

        t_56[k] = -2.0 * dg_26[k]
                  + f_0 * gg_101[k];

        t_57[k] = -2.0 * dg_27[k]
                  + f_0 * gg_102[k];

        t_58[k] = -2.0 * dg_28[k]
                  + f_0 * gg_103[k];

        t_59[k] = -2.0 * dg_29[k]
                  + f_0 * gg_104[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, dg_30, dg_31, dg_32, dg_33, dg_34, \
                         gg_105, gg_106, gg_107, gg_108, gg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -dg_30[k]
                  + f_0 * gg_105[k];

        t_61[k] = -dg_31[k]
                  + f_0 * gg_106[k];

        t_62[k] = -dg_32[k]
                  + f_0 * gg_107[k];

        t_63[k] = -dg_33[k]
                  + f_0 * gg_108[k];

        t_64[k] = -dg_34[k]
                  + f_0 * gg_109[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, dg_35, dg_36, dg_37, dg_38, dg_39, \
                         gg_110, gg_111, gg_112, gg_113, gg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -dg_35[k]
                  + f_0 * gg_110[k];

        t_66[k] = -dg_36[k]
                  + f_0 * gg_111[k];

        t_67[k] = -dg_37[k]
                  + f_0 * gg_112[k];

        t_68[k] = -dg_38[k]
                  + f_0 * gg_113[k];

        t_69[k] = -dg_39[k]
                  + f_0 * gg_114[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, dg_40, dg_41, dg_42, dg_43, dg_44, \
                         gg_115, gg_116, gg_117, gg_118, gg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -dg_40[k]
                  + f_0 * gg_115[k];

        t_71[k] = -dg_41[k]
                  + f_0 * gg_116[k];

        t_72[k] = -dg_42[k]
                  + f_0 * gg_117[k];

        t_73[k] = -dg_43[k]
                  + f_0 * gg_118[k];

        t_74[k] = -dg_44[k]
                  + f_0 * gg_119[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, t_80, t_81, t_82, gg_120, gg_121, \
                         gg_122, gg_123, gg_124, gg_125, gg_126, \
                         gg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * gg_120[k];

        t_76[k] = f_0 * gg_121[k];

        t_77[k] = f_0 * gg_122[k];

        t_78[k] = f_0 * gg_123[k];

        t_79[k] = f_0 * gg_124[k];

        t_80[k] = f_0 * gg_125[k];

        t_81[k] = f_0 * gg_126[k];

        t_82[k] = f_0 * gg_127[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, gg_128, gg_129, gg_130, \
                         gg_131, gg_132, gg_133, gg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * gg_128[k];

        t_84[k] = f_0 * gg_129[k];

        t_85[k] = f_0 * gg_130[k];

        t_86[k] = f_0 * gg_131[k];

        t_87[k] = f_0 * gg_132[k];

        t_88[k] = f_0 * gg_133[k];

        t_89[k] = f_0 * gg_134[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, dg_45, dg_46, dg_47, dg_48, dg_49, \
                         gg_150, gg_151, gg_152, gg_153, gg_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * dg_45[k]
                  + f_0 * gg_150[k];

        t_91[k] = -3.0 * dg_46[k]
                  + f_0 * gg_151[k];

        t_92[k] = -3.0 * dg_47[k]
                  + f_0 * gg_152[k];

        t_93[k] = -3.0 * dg_48[k]
                  + f_0 * gg_153[k];

        t_94[k] = -3.0 * dg_49[k]
                  + f_0 * gg_154[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, dg_50, dg_51, dg_52, dg_53, dg_54, \
                         gg_155, gg_156, gg_157, gg_158, gg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * dg_50[k]
                  + f_0 * gg_155[k];

        t_96[k] = -3.0 * dg_51[k]
                  + f_0 * gg_156[k];

        t_97[k] = -3.0 * dg_52[k]
                  + f_0 * gg_157[k];

        t_98[k] = -3.0 * dg_53[k]
                  + f_0 * gg_158[k];

        t_99[k] = -3.0 * dg_54[k]
                  + f_0 * gg_159[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, dg_55, dg_56, dg_57, dg_58, dg_59, \
                         gg_160, gg_161, gg_162, gg_163, gg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * dg_55[k]
                   + f_0 * gg_160[k];

        t_101[k] = -3.0 * dg_56[k]
                   + f_0 * gg_161[k];

        t_102[k] = -3.0 * dg_57[k]
                   + f_0 * gg_162[k];

        t_103[k] = -3.0 * dg_58[k]
                   + f_0 * gg_163[k];

        t_104[k] = -3.0 * dg_59[k]
                   + f_0 * gg_164[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, dg_60, dg_61, dg_62, dg_63, dg_64, \
                         gg_165, gg_166, gg_167, gg_168, gg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * dg_60[k]
                   + f_0 * gg_165[k];

        t_106[k] = -2.0 * dg_61[k]
                   + f_0 * gg_166[k];

        t_107[k] = -2.0 * dg_62[k]
                   + f_0 * gg_167[k];

        t_108[k] = -2.0 * dg_63[k]
                   + f_0 * gg_168[k];

        t_109[k] = -2.0 * dg_64[k]
                   + f_0 * gg_169[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, dg_65, dg_66, dg_67, dg_68, dg_69, \
                         gg_170, gg_171, gg_172, gg_173, gg_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * dg_65[k]
                   + f_0 * gg_170[k];

        t_111[k] = -2.0 * dg_66[k]
                   + f_0 * gg_171[k];

        t_112[k] = -2.0 * dg_67[k]
                   + f_0 * gg_172[k];

        t_113[k] = -2.0 * dg_68[k]
                   + f_0 * gg_173[k];

        t_114[k] = -2.0 * dg_69[k]
                   + f_0 * gg_174[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, dg_70, dg_71, dg_72, dg_73, dg_74, \
                         gg_175, gg_176, gg_177, gg_178, gg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * dg_70[k]
                   + f_0 * gg_175[k];

        t_116[k] = -2.0 * dg_71[k]
                   + f_0 * gg_176[k];

        t_117[k] = -2.0 * dg_72[k]
                   + f_0 * gg_177[k];

        t_118[k] = -2.0 * dg_73[k]
                   + f_0 * gg_178[k];

        t_119[k] = -2.0 * dg_74[k]
                   + f_0 * gg_179[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, dg_75, dg_76, dg_77, dg_78, dg_79, \
                         gg_180, gg_181, gg_182, gg_183, gg_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -dg_75[k]
                   + f_0 * gg_180[k];

        t_121[k] = -dg_76[k]
                   + f_0 * gg_181[k];

        t_122[k] = -dg_77[k]
                   + f_0 * gg_182[k];

        t_123[k] = -dg_78[k]
                   + f_0 * gg_183[k];

        t_124[k] = -dg_79[k]
                   + f_0 * gg_184[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, dg_80, dg_81, dg_82, dg_83, dg_84, \
                         gg_185, gg_186, gg_187, gg_188, gg_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -dg_80[k]
                   + f_0 * gg_185[k];

        t_126[k] = -dg_81[k]
                   + f_0 * gg_186[k];

        t_127[k] = -dg_82[k]
                   + f_0 * gg_187[k];

        t_128[k] = -dg_83[k]
                   + f_0 * gg_188[k];

        t_129[k] = -dg_84[k]
                   + f_0 * gg_189[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, dg_85, dg_86, dg_87, dg_88, dg_89, \
                         gg_190, gg_191, gg_192, gg_193, gg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -dg_85[k]
                   + f_0 * gg_190[k];

        t_131[k] = -dg_86[k]
                   + f_0 * gg_191[k];

        t_132[k] = -dg_87[k]
                   + f_0 * gg_192[k];

        t_133[k] = -dg_88[k]
                   + f_0 * gg_193[k];

        t_134[k] = -dg_89[k]
                   + f_0 * gg_194[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, t_140, t_141, t_142, gg_195, \
                         gg_196, gg_197, gg_198, gg_199, gg_200, gg_201, \
                         gg_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * gg_195[k];

        t_136[k] = f_0 * gg_196[k];

        t_137[k] = f_0 * gg_197[k];

        t_138[k] = f_0 * gg_198[k];

        t_139[k] = f_0 * gg_199[k];

        t_140[k] = f_0 * gg_200[k];

        t_141[k] = f_0 * gg_201[k];

        t_142[k] = f_0 * gg_202[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, t_148, t_149, gg_203, gg_204, \
                         gg_205, gg_206, gg_207, gg_208, gg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_0 * gg_203[k];

        t_144[k] = f_0 * gg_204[k];

        t_145[k] = f_0 * gg_205[k];

        t_146[k] = f_0 * gg_206[k];

        t_147[k] = f_0 * gg_207[k];

        t_148[k] = f_0 * gg_208[k];

        t_149[k] = f_0 * gg_209[k];
    }
}

auto
compute_prim_geom_10_fg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t dg, const size_t gg,
                                             const size_t ncols, const double alpha) -> void
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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_52 = buffer.data(dg + 52);
    const auto *dg_53 = buffer.data(dg + 53);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_60 = buffer.data(dg + 60);
    const auto *dg_61 = buffer.data(dg + 61);
    const auto *dg_62 = buffer.data(dg + 62);
    const auto *dg_63 = buffer.data(dg + 63);
    const auto *dg_64 = buffer.data(dg + 64);
    const auto *dg_65 = buffer.data(dg + 65);
    const auto *dg_66 = buffer.data(dg + 66);
    const auto *dg_67 = buffer.data(dg + 67);
    const auto *dg_68 = buffer.data(dg + 68);
    const auto *dg_69 = buffer.data(dg + 69);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_76 = buffer.data(dg + 76);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_79 = buffer.data(dg + 79);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_81 = buffer.data(dg + 81);
    const auto *dg_82 = buffer.data(dg + 82);
    const auto *dg_83 = buffer.data(dg + 83);
    const auto *dg_84 = buffer.data(dg + 84);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_88 = buffer.data(dg + 88);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_106 = buffer.data(gg + 106);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_109 = buffer.data(gg + 109);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_113 = buffer.data(gg + 113);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_124 = buffer.data(gg + 124);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_128 = buffer.data(gg + 128);
    const auto *gg_129 = buffer.data(gg + 129);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_134 = buffer.data(gg + 134);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_136 = buffer.data(gg + 136);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_138 = buffer.data(gg + 138);
    const auto *gg_139 = buffer.data(gg + 139);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_141 = buffer.data(gg + 141);
    const auto *gg_142 = buffer.data(gg + 142);
    const auto *gg_143 = buffer.data(gg + 143);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_145 = buffer.data(gg + 145);
    const auto *gg_146 = buffer.data(gg + 146);
    const auto *gg_147 = buffer.data(gg + 147);
    const auto *gg_148 = buffer.data(gg + 148);
    const auto *gg_149 = buffer.data(gg + 149);
    const auto *gg_165 = buffer.data(gg + 165);
    const auto *gg_166 = buffer.data(gg + 166);
    const auto *gg_167 = buffer.data(gg + 167);
    const auto *gg_168 = buffer.data(gg + 168);
    const auto *gg_169 = buffer.data(gg + 169);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_171 = buffer.data(gg + 171);
    const auto *gg_172 = buffer.data(gg + 172);
    const auto *gg_173 = buffer.data(gg + 173);
    const auto *gg_174 = buffer.data(gg + 174);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_181 = buffer.data(gg + 181);
    const auto *gg_182 = buffer.data(gg + 182);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_184 = buffer.data(gg + 184);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_186 = buffer.data(gg + 186);
    const auto *gg_187 = buffer.data(gg + 187);
    const auto *gg_188 = buffer.data(gg + 188);
    const auto *gg_189 = buffer.data(gg + 189);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_195 = buffer.data(gg + 195);
    const auto *gg_196 = buffer.data(gg + 196);
    const auto *gg_197 = buffer.data(gg + 197);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_199 = buffer.data(gg + 199);
    const auto *gg_200 = buffer.data(gg + 200);
    const auto *gg_201 = buffer.data(gg + 201);
    const auto *gg_202 = buffer.data(gg + 202);
    const auto *gg_203 = buffer.data(gg + 203);
    const auto *gg_204 = buffer.data(gg + 204);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_211 = buffer.data(gg + 211);
    const auto *gg_212 = buffer.data(gg + 212);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_214 = buffer.data(gg + 214);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_216 = buffer.data(gg + 216);
    const auto *gg_217 = buffer.data(gg + 217);
    const auto *gg_218 = buffer.data(gg + 218);
    const auto *gg_219 = buffer.data(gg + 219);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_221 = buffer.data(gg + 221);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_223 = buffer.data(gg + 223);
    const auto *gg_224 = buffer.data(gg + 224);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gg_30, gg_31, gg_32, gg_33, \
                         gg_34, gg_35, gg_36, gg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gg_30[k];

        t_1[k] = f_0 * gg_31[k];

        t_2[k] = f_0 * gg_32[k];

        t_3[k] = f_0 * gg_33[k];

        t_4[k] = f_0 * gg_34[k];

        t_5[k] = f_0 * gg_35[k];

        t_6[k] = f_0 * gg_36[k];

        t_7[k] = f_0 * gg_37[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, gg_38, gg_39, gg_40, \
                         gg_41, gg_42, gg_43, gg_44, gg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gg_38[k];

        t_9[k] = f_0 * gg_39[k];

        t_10[k] = f_0 * gg_40[k];

        t_11[k] = f_0 * gg_41[k];

        t_12[k] = f_0 * gg_42[k];

        t_13[k] = f_0 * gg_43[k];

        t_14[k] = f_0 * gg_44[k];

        t_15[k] = f_0 * gg_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, gg_61, gg_62, gg_63, \
                         gg_64, gg_65, gg_66, gg_67, gg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * gg_61[k];

        t_17[k] = f_0 * gg_62[k];

        t_18[k] = f_0 * gg_63[k];

        t_19[k] = f_0 * gg_64[k];

        t_20[k] = f_0 * gg_65[k];

        t_21[k] = f_0 * gg_66[k];

        t_22[k] = f_0 * gg_67[k];

        t_23[k] = f_0 * gg_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, dg_0, gg_69, gg_70, gg_71, \
                         gg_72, gg_73, gg_74, gg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * gg_69[k];

        t_25[k] = f_0 * gg_70[k];

        t_26[k] = f_0 * gg_71[k];

        t_27[k] = f_0 * gg_72[k];

        t_28[k] = f_0 * gg_73[k];

        t_29[k] = f_0 * gg_74[k];

        t_30[k] = -dg_0[k]
                  + f_0 * gg_75[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, dg_1, dg_2, dg_3, dg_4, dg_5, gg_76, \
                         gg_77, gg_78, gg_79, gg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -dg_1[k]
                  + f_0 * gg_76[k];

        t_32[k] = -dg_2[k]
                  + f_0 * gg_77[k];

        t_33[k] = -dg_3[k]
                  + f_0 * gg_78[k];

        t_34[k] = -dg_4[k]
                  + f_0 * gg_79[k];

        t_35[k] = -dg_5[k]
                  + f_0 * gg_80[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, dg_6, dg_7, dg_8, dg_9, dg_10, gg_81, \
                         gg_82, gg_83, gg_84, gg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -dg_6[k]
                  + f_0 * gg_81[k];

        t_37[k] = -dg_7[k]
                  + f_0 * gg_82[k];

        t_38[k] = -dg_8[k]
                  + f_0 * gg_83[k];

        t_39[k] = -dg_9[k]
                  + f_0 * gg_84[k];

        t_40[k] = -dg_10[k]
                  + f_0 * gg_85[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, t_46, dg_11, dg_12, dg_13, dg_14, \
                         gg_86, gg_87, gg_88, gg_89, gg_105, gg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -dg_11[k]
                  + f_0 * gg_86[k];

        t_42[k] = -dg_12[k]
                  + f_0 * gg_87[k];

        t_43[k] = -dg_13[k]
                  + f_0 * gg_88[k];

        t_44[k] = -dg_14[k]
                  + f_0 * gg_89[k];

        t_45[k] = f_0 * gg_105[k];

        t_46[k] = f_0 * gg_106[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, t_53, t_54, gg_107, gg_108, \
                         gg_109, gg_110, gg_111, gg_112, gg_113, \
                         gg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * gg_107[k];

        t_48[k] = f_0 * gg_108[k];

        t_49[k] = f_0 * gg_109[k];

        t_50[k] = f_0 * gg_110[k];

        t_51[k] = f_0 * gg_111[k];

        t_52[k] = f_0 * gg_112[k];

        t_53[k] = f_0 * gg_113[k];

        t_54[k] = f_0 * gg_114[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, dg_15, dg_16, gg_115, \
                         gg_116, gg_117, gg_118, gg_119, gg_120, \
                         gg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * gg_115[k];

        t_56[k] = f_0 * gg_116[k];

        t_57[k] = f_0 * gg_117[k];

        t_58[k] = f_0 * gg_118[k];

        t_59[k] = f_0 * gg_119[k];

        t_60[k] = -dg_15[k]
                  + f_0 * gg_120[k];

        t_61[k] = -dg_16[k]
                  + f_0 * gg_121[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, dg_17, dg_18, dg_19, dg_20, dg_21, \
                         gg_122, gg_123, gg_124, gg_125, gg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -dg_17[k]
                  + f_0 * gg_122[k];

        t_63[k] = -dg_18[k]
                  + f_0 * gg_123[k];

        t_64[k] = -dg_19[k]
                  + f_0 * gg_124[k];

        t_65[k] = -dg_20[k]
                  + f_0 * gg_125[k];

        t_66[k] = -dg_21[k]
                  + f_0 * gg_126[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, dg_22, dg_23, dg_24, dg_25, dg_26, \
                         gg_127, gg_128, gg_129, gg_130, gg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -dg_22[k]
                  + f_0 * gg_127[k];

        t_68[k] = -dg_23[k]
                  + f_0 * gg_128[k];

        t_69[k] = -dg_24[k]
                  + f_0 * gg_129[k];

        t_70[k] = -dg_25[k]
                  + f_0 * gg_130[k];

        t_71[k] = -dg_26[k]
                  + f_0 * gg_131[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, dg_27, dg_28, dg_29, dg_30, dg_31, \
                         gg_132, gg_133, gg_134, gg_135, gg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -dg_27[k]
                  + f_0 * gg_132[k];

        t_73[k] = -dg_28[k]
                  + f_0 * gg_133[k];

        t_74[k] = -dg_29[k]
                  + f_0 * gg_134[k];

        t_75[k] = -2.0 * dg_30[k]
                  + f_0 * gg_135[k];

        t_76[k] = -2.0 * dg_31[k]
                  + f_0 * gg_136[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, dg_32, dg_33, dg_34, dg_35, dg_36, \
                         gg_137, gg_138, gg_139, gg_140, gg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -2.0 * dg_32[k]
                  + f_0 * gg_137[k];

        t_78[k] = -2.0 * dg_33[k]
                  + f_0 * gg_138[k];

        t_79[k] = -2.0 * dg_34[k]
                  + f_0 * gg_139[k];

        t_80[k] = -2.0 * dg_35[k]
                  + f_0 * gg_140[k];

        t_81[k] = -2.0 * dg_36[k]
                  + f_0 * gg_141[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, dg_37, dg_38, dg_39, dg_40, dg_41, \
                         gg_142, gg_143, gg_144, gg_145, gg_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -2.0 * dg_37[k]
                  + f_0 * gg_142[k];

        t_83[k] = -2.0 * dg_38[k]
                  + f_0 * gg_143[k];

        t_84[k] = -2.0 * dg_39[k]
                  + f_0 * gg_144[k];

        t_85[k] = -2.0 * dg_40[k]
                  + f_0 * gg_145[k];

        t_86[k] = -2.0 * dg_41[k]
                  + f_0 * gg_146[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, dg_42, dg_43, dg_44, gg_147, \
                         gg_148, gg_149, gg_165, gg_166, gg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -2.0 * dg_42[k]
                  + f_0 * gg_147[k];

        t_88[k] = -2.0 * dg_43[k]
                  + f_0 * gg_148[k];

        t_89[k] = -2.0 * dg_44[k]
                  + f_0 * gg_149[k];

        t_90[k] = f_0 * gg_165[k];

        t_91[k] = f_0 * gg_166[k];

        t_92[k] = f_0 * gg_167[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, gg_168, gg_169, \
                         gg_170, gg_171, gg_172, gg_173, gg_174, \
                         gg_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * gg_168[k];

        t_94[k] = f_0 * gg_169[k];

        t_95[k] = f_0 * gg_170[k];

        t_96[k] = f_0 * gg_171[k];

        t_97[k] = f_0 * gg_172[k];

        t_98[k] = f_0 * gg_173[k];

        t_99[k] = f_0 * gg_174[k];

        t_100[k] = f_0 * gg_175[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, dg_45, dg_46, gg_176, \
                         gg_177, gg_178, gg_179, gg_180, gg_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * gg_176[k];

        t_102[k] = f_0 * gg_177[k];

        t_103[k] = f_0 * gg_178[k];

        t_104[k] = f_0 * gg_179[k];

        t_105[k] = -dg_45[k]
                   + f_0 * gg_180[k];

        t_106[k] = -dg_46[k]
                   + f_0 * gg_181[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, dg_47, dg_48, dg_49, dg_50, dg_51, \
                         gg_182, gg_183, gg_184, gg_185, gg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -dg_47[k]
                   + f_0 * gg_182[k];

        t_108[k] = -dg_48[k]
                   + f_0 * gg_183[k];

        t_109[k] = -dg_49[k]
                   + f_0 * gg_184[k];

        t_110[k] = -dg_50[k]
                   + f_0 * gg_185[k];

        t_111[k] = -dg_51[k]
                   + f_0 * gg_186[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, dg_52, dg_53, dg_54, dg_55, dg_56, \
                         gg_187, gg_188, gg_189, gg_190, gg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -dg_52[k]
                   + f_0 * gg_187[k];

        t_113[k] = -dg_53[k]
                   + f_0 * gg_188[k];

        t_114[k] = -dg_54[k]
                   + f_0 * gg_189[k];

        t_115[k] = -dg_55[k]
                   + f_0 * gg_190[k];

        t_116[k] = -dg_56[k]
                   + f_0 * gg_191[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, dg_57, dg_58, dg_59, dg_60, dg_61, \
                         gg_192, gg_193, gg_194, gg_195, gg_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -dg_57[k]
                   + f_0 * gg_192[k];

        t_118[k] = -dg_58[k]
                   + f_0 * gg_193[k];

        t_119[k] = -dg_59[k]
                   + f_0 * gg_194[k];

        t_120[k] = -2.0 * dg_60[k]
                   + f_0 * gg_195[k];

        t_121[k] = -2.0 * dg_61[k]
                   + f_0 * gg_196[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, dg_62, dg_63, dg_64, dg_65, dg_66, \
                         gg_197, gg_198, gg_199, gg_200, gg_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * dg_62[k]
                   + f_0 * gg_197[k];

        t_123[k] = -2.0 * dg_63[k]
                   + f_0 * gg_198[k];

        t_124[k] = -2.0 * dg_64[k]
                   + f_0 * gg_199[k];

        t_125[k] = -2.0 * dg_65[k]
                   + f_0 * gg_200[k];

        t_126[k] = -2.0 * dg_66[k]
                   + f_0 * gg_201[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, dg_67, dg_68, dg_69, dg_70, dg_71, \
                         gg_202, gg_203, gg_204, gg_205, gg_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * dg_67[k]
                   + f_0 * gg_202[k];

        t_128[k] = -2.0 * dg_68[k]
                   + f_0 * gg_203[k];

        t_129[k] = -2.0 * dg_69[k]
                   + f_0 * gg_204[k];

        t_130[k] = -2.0 * dg_70[k]
                   + f_0 * gg_205[k];

        t_131[k] = -2.0 * dg_71[k]
                   + f_0 * gg_206[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, dg_72, dg_73, dg_74, dg_75, dg_76, \
                         gg_207, gg_208, gg_209, gg_210, gg_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * dg_72[k]
                   + f_0 * gg_207[k];

        t_133[k] = -2.0 * dg_73[k]
                   + f_0 * gg_208[k];

        t_134[k] = -2.0 * dg_74[k]
                   + f_0 * gg_209[k];

        t_135[k] = -3.0 * dg_75[k]
                   + f_0 * gg_210[k];

        t_136[k] = -3.0 * dg_76[k]
                   + f_0 * gg_211[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, dg_77, dg_78, dg_79, dg_80, dg_81, \
                         gg_212, gg_213, gg_214, gg_215, gg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -3.0 * dg_77[k]
                   + f_0 * gg_212[k];

        t_138[k] = -3.0 * dg_78[k]
                   + f_0 * gg_213[k];

        t_139[k] = -3.0 * dg_79[k]
                   + f_0 * gg_214[k];

        t_140[k] = -3.0 * dg_80[k]
                   + f_0 * gg_215[k];

        t_141[k] = -3.0 * dg_81[k]
                   + f_0 * gg_216[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, dg_82, dg_83, dg_84, dg_85, dg_86, \
                         gg_217, gg_218, gg_219, gg_220, gg_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -3.0 * dg_82[k]
                   + f_0 * gg_217[k];

        t_143[k] = -3.0 * dg_83[k]
                   + f_0 * gg_218[k];

        t_144[k] = -3.0 * dg_84[k]
                   + f_0 * gg_219[k];

        t_145[k] = -3.0 * dg_85[k]
                   + f_0 * gg_220[k];

        t_146[k] = -3.0 * dg_86[k]
                   + f_0 * gg_221[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, dg_87, dg_88, dg_89, gg_222, gg_223, \
                         gg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -3.0 * dg_87[k]
                   + f_0 * gg_222[k];

        t_148[k] = -3.0 * dg_88[k]
                   + f_0 * gg_223[k];

        t_149[k] = -3.0 * dg_89[k]
                   + f_0 * gg_224[k];
    }
}

}  // namespace simdt2ceri
