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


#include "SimdElectronRepulsionGeom10VrrRecKK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_kk_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_0 = buffer.data(ik + 0);
    const auto *ik_1 = buffer.data(ik + 1);
    const auto *ik_2 = buffer.data(ik + 2);
    const auto *ik_3 = buffer.data(ik + 3);
    const auto *ik_4 = buffer.data(ik + 4);
    const auto *ik_5 = buffer.data(ik + 5);
    const auto *ik_6 = buffer.data(ik + 6);
    const auto *ik_7 = buffer.data(ik + 7);
    const auto *ik_8 = buffer.data(ik + 8);
    const auto *ik_9 = buffer.data(ik + 9);
    const auto *ik_10 = buffer.data(ik + 10);
    const auto *ik_11 = buffer.data(ik + 11);
    const auto *ik_12 = buffer.data(ik + 12);
    const auto *ik_13 = buffer.data(ik + 13);
    const auto *ik_14 = buffer.data(ik + 14);
    const auto *ik_15 = buffer.data(ik + 15);
    const auto *ik_16 = buffer.data(ik + 16);
    const auto *ik_17 = buffer.data(ik + 17);
    const auto *ik_18 = buffer.data(ik + 18);
    const auto *ik_19 = buffer.data(ik + 19);
    const auto *ik_20 = buffer.data(ik + 20);
    const auto *ik_21 = buffer.data(ik + 21);
    const auto *ik_22 = buffer.data(ik + 22);
    const auto *ik_23 = buffer.data(ik + 23);
    const auto *ik_24 = buffer.data(ik + 24);
    const auto *ik_25 = buffer.data(ik + 25);
    const auto *ik_26 = buffer.data(ik + 26);
    const auto *ik_27 = buffer.data(ik + 27);
    const auto *ik_28 = buffer.data(ik + 28);
    const auto *ik_29 = buffer.data(ik + 29);
    const auto *ik_30 = buffer.data(ik + 30);
    const auto *ik_31 = buffer.data(ik + 31);
    const auto *ik_32 = buffer.data(ik + 32);
    const auto *ik_33 = buffer.data(ik + 33);
    const auto *ik_34 = buffer.data(ik + 34);
    const auto *ik_35 = buffer.data(ik + 35);
    const auto *ik_36 = buffer.data(ik + 36);
    const auto *ik_37 = buffer.data(ik + 37);
    const auto *ik_38 = buffer.data(ik + 38);
    const auto *ik_39 = buffer.data(ik + 39);
    const auto *ik_40 = buffer.data(ik + 40);
    const auto *ik_41 = buffer.data(ik + 41);
    const auto *ik_42 = buffer.data(ik + 42);
    const auto *ik_43 = buffer.data(ik + 43);
    const auto *ik_44 = buffer.data(ik + 44);
    const auto *ik_45 = buffer.data(ik + 45);
    const auto *ik_46 = buffer.data(ik + 46);
    const auto *ik_47 = buffer.data(ik + 47);
    const auto *ik_48 = buffer.data(ik + 48);
    const auto *ik_49 = buffer.data(ik + 49);
    const auto *ik_50 = buffer.data(ik + 50);
    const auto *ik_51 = buffer.data(ik + 51);
    const auto *ik_52 = buffer.data(ik + 52);
    const auto *ik_53 = buffer.data(ik + 53);
    const auto *ik_54 = buffer.data(ik + 54);
    const auto *ik_55 = buffer.data(ik + 55);
    const auto *ik_56 = buffer.data(ik + 56);
    const auto *ik_57 = buffer.data(ik + 57);
    const auto *ik_58 = buffer.data(ik + 58);
    const auto *ik_59 = buffer.data(ik + 59);
    const auto *ik_60 = buffer.data(ik + 60);
    const auto *ik_61 = buffer.data(ik + 61);
    const auto *ik_62 = buffer.data(ik + 62);
    const auto *ik_63 = buffer.data(ik + 63);
    const auto *ik_64 = buffer.data(ik + 64);
    const auto *ik_65 = buffer.data(ik + 65);
    const auto *ik_66 = buffer.data(ik + 66);
    const auto *ik_67 = buffer.data(ik + 67);
    const auto *ik_68 = buffer.data(ik + 68);
    const auto *ik_69 = buffer.data(ik + 69);
    const auto *ik_70 = buffer.data(ik + 70);
    const auto *ik_71 = buffer.data(ik + 71);
    const auto *ik_72 = buffer.data(ik + 72);
    const auto *ik_73 = buffer.data(ik + 73);
    const auto *ik_74 = buffer.data(ik + 74);
    const auto *ik_75 = buffer.data(ik + 75);
    const auto *ik_76 = buffer.data(ik + 76);
    const auto *ik_77 = buffer.data(ik + 77);
    const auto *ik_78 = buffer.data(ik + 78);
    const auto *ik_79 = buffer.data(ik + 79);
    const auto *ik_80 = buffer.data(ik + 80);
    const auto *ik_81 = buffer.data(ik + 81);
    const auto *ik_82 = buffer.data(ik + 82);
    const auto *ik_83 = buffer.data(ik + 83);
    const auto *ik_84 = buffer.data(ik + 84);
    const auto *ik_85 = buffer.data(ik + 85);
    const auto *ik_86 = buffer.data(ik + 86);
    const auto *ik_87 = buffer.data(ik + 87);
    const auto *ik_88 = buffer.data(ik + 88);
    const auto *ik_89 = buffer.data(ik + 89);
    const auto *ik_90 = buffer.data(ik + 90);
    const auto *ik_91 = buffer.data(ik + 91);
    const auto *ik_92 = buffer.data(ik + 92);
    const auto *ik_93 = buffer.data(ik + 93);
    const auto *ik_94 = buffer.data(ik + 94);
    const auto *ik_95 = buffer.data(ik + 95);
    const auto *ik_96 = buffer.data(ik + 96);
    const auto *ik_97 = buffer.data(ik + 97);
    const auto *ik_98 = buffer.data(ik + 98);
    const auto *ik_99 = buffer.data(ik + 99);
    const auto *ik_100 = buffer.data(ik + 100);
    const auto *ik_101 = buffer.data(ik + 101);
    const auto *ik_102 = buffer.data(ik + 102);
    const auto *ik_103 = buffer.data(ik + 103);
    const auto *ik_104 = buffer.data(ik + 104);
    const auto *ik_105 = buffer.data(ik + 105);
    const auto *ik_106 = buffer.data(ik + 106);
    const auto *ik_107 = buffer.data(ik + 107);
    const auto *ik_108 = buffer.data(ik + 108);
    const auto *ik_109 = buffer.data(ik + 109);
    const auto *ik_110 = buffer.data(ik + 110);
    const auto *ik_111 = buffer.data(ik + 111);
    const auto *ik_112 = buffer.data(ik + 112);
    const auto *ik_113 = buffer.data(ik + 113);
    const auto *ik_114 = buffer.data(ik + 114);
    const auto *ik_115 = buffer.data(ik + 115);
    const auto *ik_116 = buffer.data(ik + 116);
    const auto *ik_117 = buffer.data(ik + 117);
    const auto *ik_118 = buffer.data(ik + 118);
    const auto *ik_119 = buffer.data(ik + 119);
    const auto *ik_120 = buffer.data(ik + 120);
    const auto *ik_121 = buffer.data(ik + 121);
    const auto *ik_122 = buffer.data(ik + 122);
    const auto *ik_123 = buffer.data(ik + 123);
    const auto *ik_124 = buffer.data(ik + 124);
    const auto *ik_125 = buffer.data(ik + 125);
    const auto *ik_126 = buffer.data(ik + 126);
    const auto *ik_127 = buffer.data(ik + 127);
    const auto *ik_128 = buffer.data(ik + 128);
    const auto *ik_129 = buffer.data(ik + 129);
    const auto *ik_130 = buffer.data(ik + 130);
    const auto *ik_131 = buffer.data(ik + 131);
    const auto *ik_132 = buffer.data(ik + 132);
    const auto *ik_133 = buffer.data(ik + 133);
    const auto *ik_134 = buffer.data(ik + 134);
    const auto *ik_135 = buffer.data(ik + 135);
    const auto *ik_136 = buffer.data(ik + 136);
    const auto *ik_137 = buffer.data(ik + 137);
    const auto *ik_138 = buffer.data(ik + 138);
    const auto *ik_139 = buffer.data(ik + 139);
    const auto *ik_140 = buffer.data(ik + 140);
    const auto *ik_141 = buffer.data(ik + 141);
    const auto *ik_142 = buffer.data(ik + 142);
    const auto *ik_143 = buffer.data(ik + 143);
    const auto *ik_144 = buffer.data(ik + 144);
    const auto *ik_145 = buffer.data(ik + 145);
    const auto *ik_146 = buffer.data(ik + 146);
    const auto *ik_147 = buffer.data(ik + 147);
    const auto *ik_148 = buffer.data(ik + 148);
    const auto *ik_149 = buffer.data(ik + 149);

    const auto *lk_0 = buffer.data(lk + 0);
    const auto *lk_1 = buffer.data(lk + 1);
    const auto *lk_2 = buffer.data(lk + 2);
    const auto *lk_3 = buffer.data(lk + 3);
    const auto *lk_4 = buffer.data(lk + 4);
    const auto *lk_5 = buffer.data(lk + 5);
    const auto *lk_6 = buffer.data(lk + 6);
    const auto *lk_7 = buffer.data(lk + 7);
    const auto *lk_8 = buffer.data(lk + 8);
    const auto *lk_9 = buffer.data(lk + 9);
    const auto *lk_10 = buffer.data(lk + 10);
    const auto *lk_11 = buffer.data(lk + 11);
    const auto *lk_12 = buffer.data(lk + 12);
    const auto *lk_13 = buffer.data(lk + 13);
    const auto *lk_14 = buffer.data(lk + 14);
    const auto *lk_15 = buffer.data(lk + 15);
    const auto *lk_16 = buffer.data(lk + 16);
    const auto *lk_17 = buffer.data(lk + 17);
    const auto *lk_18 = buffer.data(lk + 18);
    const auto *lk_19 = buffer.data(lk + 19);
    const auto *lk_20 = buffer.data(lk + 20);
    const auto *lk_21 = buffer.data(lk + 21);
    const auto *lk_22 = buffer.data(lk + 22);
    const auto *lk_23 = buffer.data(lk + 23);
    const auto *lk_24 = buffer.data(lk + 24);
    const auto *lk_25 = buffer.data(lk + 25);
    const auto *lk_26 = buffer.data(lk + 26);
    const auto *lk_27 = buffer.data(lk + 27);
    const auto *lk_28 = buffer.data(lk + 28);
    const auto *lk_29 = buffer.data(lk + 29);
    const auto *lk_30 = buffer.data(lk + 30);
    const auto *lk_31 = buffer.data(lk + 31);
    const auto *lk_32 = buffer.data(lk + 32);
    const auto *lk_33 = buffer.data(lk + 33);
    const auto *lk_34 = buffer.data(lk + 34);
    const auto *lk_35 = buffer.data(lk + 35);
    const auto *lk_36 = buffer.data(lk + 36);
    const auto *lk_37 = buffer.data(lk + 37);
    const auto *lk_38 = buffer.data(lk + 38);
    const auto *lk_39 = buffer.data(lk + 39);
    const auto *lk_40 = buffer.data(lk + 40);
    const auto *lk_41 = buffer.data(lk + 41);
    const auto *lk_42 = buffer.data(lk + 42);
    const auto *lk_43 = buffer.data(lk + 43);
    const auto *lk_44 = buffer.data(lk + 44);
    const auto *lk_45 = buffer.data(lk + 45);
    const auto *lk_46 = buffer.data(lk + 46);
    const auto *lk_47 = buffer.data(lk + 47);
    const auto *lk_48 = buffer.data(lk + 48);
    const auto *lk_49 = buffer.data(lk + 49);
    const auto *lk_50 = buffer.data(lk + 50);
    const auto *lk_51 = buffer.data(lk + 51);
    const auto *lk_52 = buffer.data(lk + 52);
    const auto *lk_53 = buffer.data(lk + 53);
    const auto *lk_54 = buffer.data(lk + 54);
    const auto *lk_55 = buffer.data(lk + 55);
    const auto *lk_56 = buffer.data(lk + 56);
    const auto *lk_57 = buffer.data(lk + 57);
    const auto *lk_58 = buffer.data(lk + 58);
    const auto *lk_59 = buffer.data(lk + 59);
    const auto *lk_60 = buffer.data(lk + 60);
    const auto *lk_61 = buffer.data(lk + 61);
    const auto *lk_62 = buffer.data(lk + 62);
    const auto *lk_63 = buffer.data(lk + 63);
    const auto *lk_64 = buffer.data(lk + 64);
    const auto *lk_65 = buffer.data(lk + 65);
    const auto *lk_66 = buffer.data(lk + 66);
    const auto *lk_67 = buffer.data(lk + 67);
    const auto *lk_68 = buffer.data(lk + 68);
    const auto *lk_69 = buffer.data(lk + 69);
    const auto *lk_70 = buffer.data(lk + 70);
    const auto *lk_71 = buffer.data(lk + 71);
    const auto *lk_72 = buffer.data(lk + 72);
    const auto *lk_73 = buffer.data(lk + 73);
    const auto *lk_74 = buffer.data(lk + 74);
    const auto *lk_75 = buffer.data(lk + 75);
    const auto *lk_76 = buffer.data(lk + 76);
    const auto *lk_77 = buffer.data(lk + 77);
    const auto *lk_78 = buffer.data(lk + 78);
    const auto *lk_79 = buffer.data(lk + 79);
    const auto *lk_80 = buffer.data(lk + 80);
    const auto *lk_81 = buffer.data(lk + 81);
    const auto *lk_82 = buffer.data(lk + 82);
    const auto *lk_83 = buffer.data(lk + 83);
    const auto *lk_84 = buffer.data(lk + 84);
    const auto *lk_85 = buffer.data(lk + 85);
    const auto *lk_86 = buffer.data(lk + 86);
    const auto *lk_87 = buffer.data(lk + 87);
    const auto *lk_88 = buffer.data(lk + 88);
    const auto *lk_89 = buffer.data(lk + 89);
    const auto *lk_90 = buffer.data(lk + 90);
    const auto *lk_91 = buffer.data(lk + 91);
    const auto *lk_92 = buffer.data(lk + 92);
    const auto *lk_93 = buffer.data(lk + 93);
    const auto *lk_94 = buffer.data(lk + 94);
    const auto *lk_95 = buffer.data(lk + 95);
    const auto *lk_96 = buffer.data(lk + 96);
    const auto *lk_97 = buffer.data(lk + 97);
    const auto *lk_98 = buffer.data(lk + 98);
    const auto *lk_99 = buffer.data(lk + 99);
    const auto *lk_100 = buffer.data(lk + 100);
    const auto *lk_101 = buffer.data(lk + 101);
    const auto *lk_102 = buffer.data(lk + 102);
    const auto *lk_103 = buffer.data(lk + 103);
    const auto *lk_104 = buffer.data(lk + 104);
    const auto *lk_105 = buffer.data(lk + 105);
    const auto *lk_106 = buffer.data(lk + 106);
    const auto *lk_107 = buffer.data(lk + 107);
    const auto *lk_108 = buffer.data(lk + 108);
    const auto *lk_109 = buffer.data(lk + 109);
    const auto *lk_110 = buffer.data(lk + 110);
    const auto *lk_111 = buffer.data(lk + 111);
    const auto *lk_112 = buffer.data(lk + 112);
    const auto *lk_113 = buffer.data(lk + 113);
    const auto *lk_114 = buffer.data(lk + 114);
    const auto *lk_115 = buffer.data(lk + 115);
    const auto *lk_116 = buffer.data(lk + 116);
    const auto *lk_117 = buffer.data(lk + 117);
    const auto *lk_118 = buffer.data(lk + 118);
    const auto *lk_119 = buffer.data(lk + 119);
    const auto *lk_120 = buffer.data(lk + 120);
    const auto *lk_121 = buffer.data(lk + 121);
    const auto *lk_122 = buffer.data(lk + 122);
    const auto *lk_123 = buffer.data(lk + 123);
    const auto *lk_124 = buffer.data(lk + 124);
    const auto *lk_125 = buffer.data(lk + 125);
    const auto *lk_126 = buffer.data(lk + 126);
    const auto *lk_127 = buffer.data(lk + 127);
    const auto *lk_128 = buffer.data(lk + 128);
    const auto *lk_129 = buffer.data(lk + 129);
    const auto *lk_130 = buffer.data(lk + 130);
    const auto *lk_131 = buffer.data(lk + 131);
    const auto *lk_132 = buffer.data(lk + 132);
    const auto *lk_133 = buffer.data(lk + 133);
    const auto *lk_134 = buffer.data(lk + 134);
    const auto *lk_135 = buffer.data(lk + 135);
    const auto *lk_136 = buffer.data(lk + 136);
    const auto *lk_137 = buffer.data(lk + 137);
    const auto *lk_138 = buffer.data(lk + 138);
    const auto *lk_139 = buffer.data(lk + 139);
    const auto *lk_140 = buffer.data(lk + 140);
    const auto *lk_141 = buffer.data(lk + 141);
    const auto *lk_142 = buffer.data(lk + 142);
    const auto *lk_143 = buffer.data(lk + 143);
    const auto *lk_144 = buffer.data(lk + 144);
    const auto *lk_145 = buffer.data(lk + 145);
    const auto *lk_146 = buffer.data(lk + 146);
    const auto *lk_147 = buffer.data(lk + 147);
    const auto *lk_148 = buffer.data(lk + 148);
    const auto *lk_149 = buffer.data(lk + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ik_0, ik_1, ik_2, ik_3, ik_4, lk_0, lk_1, \
                         lk_2, lk_3, lk_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -7.0 * ik_0[k]
                 + f_0 * lk_0[k];

        t_1[k] = -7.0 * ik_1[k]
                 + f_0 * lk_1[k];

        t_2[k] = -7.0 * ik_2[k]
                 + f_0 * lk_2[k];

        t_3[k] = -7.0 * ik_3[k]
                 + f_0 * lk_3[k];

        t_4[k] = -7.0 * ik_4[k]
                 + f_0 * lk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ik_5, ik_6, ik_7, ik_8, ik_9, lk_5, lk_6, \
                         lk_7, lk_8, lk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -7.0 * ik_5[k]
                 + f_0 * lk_5[k];

        t_6[k] = -7.0 * ik_6[k]
                 + f_0 * lk_6[k];

        t_7[k] = -7.0 * ik_7[k]
                 + f_0 * lk_7[k];

        t_8[k] = -7.0 * ik_8[k]
                 + f_0 * lk_8[k];

        t_9[k] = -7.0 * ik_9[k]
                 + f_0 * lk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ik_10, ik_11, ik_12, ik_13, ik_14, \
                         lk_10, lk_11, lk_12, lk_13, lk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -7.0 * ik_10[k]
                  + f_0 * lk_10[k];

        t_11[k] = -7.0 * ik_11[k]
                  + f_0 * lk_11[k];

        t_12[k] = -7.0 * ik_12[k]
                  + f_0 * lk_12[k];

        t_13[k] = -7.0 * ik_13[k]
                  + f_0 * lk_13[k];

        t_14[k] = -7.0 * ik_14[k]
                  + f_0 * lk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ik_15, ik_16, ik_17, ik_18, ik_19, \
                         lk_15, lk_16, lk_17, lk_18, lk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -7.0 * ik_15[k]
                  + f_0 * lk_15[k];

        t_16[k] = -7.0 * ik_16[k]
                  + f_0 * lk_16[k];

        t_17[k] = -7.0 * ik_17[k]
                  + f_0 * lk_17[k];

        t_18[k] = -7.0 * ik_18[k]
                  + f_0 * lk_18[k];

        t_19[k] = -7.0 * ik_19[k]
                  + f_0 * lk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ik_20, ik_21, ik_22, ik_23, ik_24, \
                         lk_20, lk_21, lk_22, lk_23, lk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -7.0 * ik_20[k]
                  + f_0 * lk_20[k];

        t_21[k] = -7.0 * ik_21[k]
                  + f_0 * lk_21[k];

        t_22[k] = -7.0 * ik_22[k]
                  + f_0 * lk_22[k];

        t_23[k] = -7.0 * ik_23[k]
                  + f_0 * lk_23[k];

        t_24[k] = -7.0 * ik_24[k]
                  + f_0 * lk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ik_25, ik_26, ik_27, ik_28, ik_29, \
                         lk_25, lk_26, lk_27, lk_28, lk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -7.0 * ik_25[k]
                  + f_0 * lk_25[k];

        t_26[k] = -7.0 * ik_26[k]
                  + f_0 * lk_26[k];

        t_27[k] = -7.0 * ik_27[k]
                  + f_0 * lk_27[k];

        t_28[k] = -7.0 * ik_28[k]
                  + f_0 * lk_28[k];

        t_29[k] = -7.0 * ik_29[k]
                  + f_0 * lk_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ik_30, ik_31, ik_32, ik_33, ik_34, \
                         lk_30, lk_31, lk_32, lk_33, lk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -7.0 * ik_30[k]
                  + f_0 * lk_30[k];

        t_31[k] = -7.0 * ik_31[k]
                  + f_0 * lk_31[k];

        t_32[k] = -7.0 * ik_32[k]
                  + f_0 * lk_32[k];

        t_33[k] = -7.0 * ik_33[k]
                  + f_0 * lk_33[k];

        t_34[k] = -7.0 * ik_34[k]
                  + f_0 * lk_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ik_35, ik_36, ik_37, ik_38, ik_39, \
                         lk_35, lk_36, lk_37, lk_38, lk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -7.0 * ik_35[k]
                  + f_0 * lk_35[k];

        t_36[k] = -6.0 * ik_36[k]
                  + f_0 * lk_36[k];

        t_37[k] = -6.0 * ik_37[k]
                  + f_0 * lk_37[k];

        t_38[k] = -6.0 * ik_38[k]
                  + f_0 * lk_38[k];

        t_39[k] = -6.0 * ik_39[k]
                  + f_0 * lk_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ik_40, ik_41, ik_42, ik_43, ik_44, \
                         lk_40, lk_41, lk_42, lk_43, lk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -6.0 * ik_40[k]
                  + f_0 * lk_40[k];

        t_41[k] = -6.0 * ik_41[k]
                  + f_0 * lk_41[k];

        t_42[k] = -6.0 * ik_42[k]
                  + f_0 * lk_42[k];

        t_43[k] = -6.0 * ik_43[k]
                  + f_0 * lk_43[k];

        t_44[k] = -6.0 * ik_44[k]
                  + f_0 * lk_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ik_45, ik_46, ik_47, ik_48, ik_49, \
                         lk_45, lk_46, lk_47, lk_48, lk_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -6.0 * ik_45[k]
                  + f_0 * lk_45[k];

        t_46[k] = -6.0 * ik_46[k]
                  + f_0 * lk_46[k];

        t_47[k] = -6.0 * ik_47[k]
                  + f_0 * lk_47[k];

        t_48[k] = -6.0 * ik_48[k]
                  + f_0 * lk_48[k];

        t_49[k] = -6.0 * ik_49[k]
                  + f_0 * lk_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ik_50, ik_51, ik_52, ik_53, ik_54, \
                         lk_50, lk_51, lk_52, lk_53, lk_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -6.0 * ik_50[k]
                  + f_0 * lk_50[k];

        t_51[k] = -6.0 * ik_51[k]
                  + f_0 * lk_51[k];

        t_52[k] = -6.0 * ik_52[k]
                  + f_0 * lk_52[k];

        t_53[k] = -6.0 * ik_53[k]
                  + f_0 * lk_53[k];

        t_54[k] = -6.0 * ik_54[k]
                  + f_0 * lk_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ik_55, ik_56, ik_57, ik_58, ik_59, \
                         lk_55, lk_56, lk_57, lk_58, lk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -6.0 * ik_55[k]
                  + f_0 * lk_55[k];

        t_56[k] = -6.0 * ik_56[k]
                  + f_0 * lk_56[k];

        t_57[k] = -6.0 * ik_57[k]
                  + f_0 * lk_57[k];

        t_58[k] = -6.0 * ik_58[k]
                  + f_0 * lk_58[k];

        t_59[k] = -6.0 * ik_59[k]
                  + f_0 * lk_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ik_60, ik_61, ik_62, ik_63, ik_64, \
                         lk_60, lk_61, lk_62, lk_63, lk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -6.0 * ik_60[k]
                  + f_0 * lk_60[k];

        t_61[k] = -6.0 * ik_61[k]
                  + f_0 * lk_61[k];

        t_62[k] = -6.0 * ik_62[k]
                  + f_0 * lk_62[k];

        t_63[k] = -6.0 * ik_63[k]
                  + f_0 * lk_63[k];

        t_64[k] = -6.0 * ik_64[k]
                  + f_0 * lk_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ik_65, ik_66, ik_67, ik_68, ik_69, \
                         lk_65, lk_66, lk_67, lk_68, lk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -6.0 * ik_65[k]
                  + f_0 * lk_65[k];

        t_66[k] = -6.0 * ik_66[k]
                  + f_0 * lk_66[k];

        t_67[k] = -6.0 * ik_67[k]
                  + f_0 * lk_67[k];

        t_68[k] = -6.0 * ik_68[k]
                  + f_0 * lk_68[k];

        t_69[k] = -6.0 * ik_69[k]
                  + f_0 * lk_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ik_70, ik_71, ik_72, ik_73, ik_74, \
                         lk_70, lk_71, lk_72, lk_73, lk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -6.0 * ik_70[k]
                  + f_0 * lk_70[k];

        t_71[k] = -6.0 * ik_71[k]
                  + f_0 * lk_71[k];

        t_72[k] = -6.0 * ik_72[k]
                  + f_0 * lk_72[k];

        t_73[k] = -6.0 * ik_73[k]
                  + f_0 * lk_73[k];

        t_74[k] = -6.0 * ik_74[k]
                  + f_0 * lk_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ik_75, ik_76, ik_77, ik_78, ik_79, \
                         lk_75, lk_76, lk_77, lk_78, lk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -6.0 * ik_75[k]
                  + f_0 * lk_75[k];

        t_76[k] = -6.0 * ik_76[k]
                  + f_0 * lk_76[k];

        t_77[k] = -6.0 * ik_77[k]
                  + f_0 * lk_77[k];

        t_78[k] = -6.0 * ik_78[k]
                  + f_0 * lk_78[k];

        t_79[k] = -6.0 * ik_79[k]
                  + f_0 * lk_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ik_80, ik_81, ik_82, ik_83, ik_84, \
                         lk_80, lk_81, lk_82, lk_83, lk_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -6.0 * ik_80[k]
                  + f_0 * lk_80[k];

        t_81[k] = -6.0 * ik_81[k]
                  + f_0 * lk_81[k];

        t_82[k] = -6.0 * ik_82[k]
                  + f_0 * lk_82[k];

        t_83[k] = -6.0 * ik_83[k]
                  + f_0 * lk_83[k];

        t_84[k] = -6.0 * ik_84[k]
                  + f_0 * lk_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ik_85, ik_86, ik_87, ik_88, ik_89, \
                         lk_85, lk_86, lk_87, lk_88, lk_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -6.0 * ik_85[k]
                  + f_0 * lk_85[k];

        t_86[k] = -6.0 * ik_86[k]
                  + f_0 * lk_86[k];

        t_87[k] = -6.0 * ik_87[k]
                  + f_0 * lk_87[k];

        t_88[k] = -6.0 * ik_88[k]
                  + f_0 * lk_88[k];

        t_89[k] = -6.0 * ik_89[k]
                  + f_0 * lk_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ik_90, ik_91, ik_92, ik_93, ik_94, \
                         lk_90, lk_91, lk_92, lk_93, lk_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -6.0 * ik_90[k]
                  + f_0 * lk_90[k];

        t_91[k] = -6.0 * ik_91[k]
                  + f_0 * lk_91[k];

        t_92[k] = -6.0 * ik_92[k]
                  + f_0 * lk_92[k];

        t_93[k] = -6.0 * ik_93[k]
                  + f_0 * lk_93[k];

        t_94[k] = -6.0 * ik_94[k]
                  + f_0 * lk_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ik_95, ik_96, ik_97, ik_98, ik_99, \
                         lk_95, lk_96, lk_97, lk_98, lk_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -6.0 * ik_95[k]
                  + f_0 * lk_95[k];

        t_96[k] = -6.0 * ik_96[k]
                  + f_0 * lk_96[k];

        t_97[k] = -6.0 * ik_97[k]
                  + f_0 * lk_97[k];

        t_98[k] = -6.0 * ik_98[k]
                  + f_0 * lk_98[k];

        t_99[k] = -6.0 * ik_99[k]
                  + f_0 * lk_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ik_100, ik_101, ik_102, ik_103, \
                         ik_104, lk_100, lk_101, lk_102, lk_103, \
                         lk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -6.0 * ik_100[k]
                   + f_0 * lk_100[k];

        t_101[k] = -6.0 * ik_101[k]
                   + f_0 * lk_101[k];

        t_102[k] = -6.0 * ik_102[k]
                   + f_0 * lk_102[k];

        t_103[k] = -6.0 * ik_103[k]
                   + f_0 * lk_103[k];

        t_104[k] = -6.0 * ik_104[k]
                   + f_0 * lk_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ik_105, ik_106, ik_107, ik_108, \
                         ik_109, lk_105, lk_106, lk_107, lk_108, \
                         lk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -6.0 * ik_105[k]
                   + f_0 * lk_105[k];

        t_106[k] = -6.0 * ik_106[k]
                   + f_0 * lk_106[k];

        t_107[k] = -6.0 * ik_107[k]
                   + f_0 * lk_107[k];

        t_108[k] = -5.0 * ik_108[k]
                   + f_0 * lk_108[k];

        t_109[k] = -5.0 * ik_109[k]
                   + f_0 * lk_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ik_110, ik_111, ik_112, ik_113, \
                         ik_114, lk_110, lk_111, lk_112, lk_113, \
                         lk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -5.0 * ik_110[k]
                   + f_0 * lk_110[k];

        t_111[k] = -5.0 * ik_111[k]
                   + f_0 * lk_111[k];

        t_112[k] = -5.0 * ik_112[k]
                   + f_0 * lk_112[k];

        t_113[k] = -5.0 * ik_113[k]
                   + f_0 * lk_113[k];

        t_114[k] = -5.0 * ik_114[k]
                   + f_0 * lk_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ik_115, ik_116, ik_117, ik_118, \
                         ik_119, lk_115, lk_116, lk_117, lk_118, \
                         lk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -5.0 * ik_115[k]
                   + f_0 * lk_115[k];

        t_116[k] = -5.0 * ik_116[k]
                   + f_0 * lk_116[k];

        t_117[k] = -5.0 * ik_117[k]
                   + f_0 * lk_117[k];

        t_118[k] = -5.0 * ik_118[k]
                   + f_0 * lk_118[k];

        t_119[k] = -5.0 * ik_119[k]
                   + f_0 * lk_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ik_120, ik_121, ik_122, ik_123, \
                         ik_124, lk_120, lk_121, lk_122, lk_123, \
                         lk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -5.0 * ik_120[k]
                   + f_0 * lk_120[k];

        t_121[k] = -5.0 * ik_121[k]
                   + f_0 * lk_121[k];

        t_122[k] = -5.0 * ik_122[k]
                   + f_0 * lk_122[k];

        t_123[k] = -5.0 * ik_123[k]
                   + f_0 * lk_123[k];

        t_124[k] = -5.0 * ik_124[k]
                   + f_0 * lk_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ik_125, ik_126, ik_127, ik_128, \
                         ik_129, lk_125, lk_126, lk_127, lk_128, \
                         lk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -5.0 * ik_125[k]
                   + f_0 * lk_125[k];

        t_126[k] = -5.0 * ik_126[k]
                   + f_0 * lk_126[k];

        t_127[k] = -5.0 * ik_127[k]
                   + f_0 * lk_127[k];

        t_128[k] = -5.0 * ik_128[k]
                   + f_0 * lk_128[k];

        t_129[k] = -5.0 * ik_129[k]
                   + f_0 * lk_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ik_130, ik_131, ik_132, ik_133, \
                         ik_134, lk_130, lk_131, lk_132, lk_133, \
                         lk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -5.0 * ik_130[k]
                   + f_0 * lk_130[k];

        t_131[k] = -5.0 * ik_131[k]
                   + f_0 * lk_131[k];

        t_132[k] = -5.0 * ik_132[k]
                   + f_0 * lk_132[k];

        t_133[k] = -5.0 * ik_133[k]
                   + f_0 * lk_133[k];

        t_134[k] = -5.0 * ik_134[k]
                   + f_0 * lk_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ik_135, ik_136, ik_137, ik_138, \
                         ik_139, lk_135, lk_136, lk_137, lk_138, \
                         lk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -5.0 * ik_135[k]
                   + f_0 * lk_135[k];

        t_136[k] = -5.0 * ik_136[k]
                   + f_0 * lk_136[k];

        t_137[k] = -5.0 * ik_137[k]
                   + f_0 * lk_137[k];

        t_138[k] = -5.0 * ik_138[k]
                   + f_0 * lk_138[k];

        t_139[k] = -5.0 * ik_139[k]
                   + f_0 * lk_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ik_140, ik_141, ik_142, ik_143, \
                         ik_144, lk_140, lk_141, lk_142, lk_143, \
                         lk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -5.0 * ik_140[k]
                   + f_0 * lk_140[k];

        t_141[k] = -5.0 * ik_141[k]
                   + f_0 * lk_141[k];

        t_142[k] = -5.0 * ik_142[k]
                   + f_0 * lk_142[k];

        t_143[k] = -5.0 * ik_143[k]
                   + f_0 * lk_143[k];

        t_144[k] = -5.0 * ik_144[k]
                   + f_0 * lk_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ik_145, ik_146, ik_147, ik_148, \
                         ik_149, lk_145, lk_146, lk_147, lk_148, \
                         lk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -5.0 * ik_145[k]
                   + f_0 * lk_145[k];

        t_146[k] = -5.0 * ik_146[k]
                   + f_0 * lk_146[k];

        t_147[k] = -5.0 * ik_147[k]
                   + f_0 * lk_147[k];

        t_148[k] = -5.0 * ik_148[k]
                   + f_0 * lk_148[k];

        t_149[k] = -5.0 * ik_149[k]
                   + f_0 * lk_149[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_150 = buffer.data(ik + 150);
    const auto *ik_151 = buffer.data(ik + 151);
    const auto *ik_152 = buffer.data(ik + 152);
    const auto *ik_153 = buffer.data(ik + 153);
    const auto *ik_154 = buffer.data(ik + 154);
    const auto *ik_155 = buffer.data(ik + 155);
    const auto *ik_156 = buffer.data(ik + 156);
    const auto *ik_157 = buffer.data(ik + 157);
    const auto *ik_158 = buffer.data(ik + 158);
    const auto *ik_159 = buffer.data(ik + 159);
    const auto *ik_160 = buffer.data(ik + 160);
    const auto *ik_161 = buffer.data(ik + 161);
    const auto *ik_162 = buffer.data(ik + 162);
    const auto *ik_163 = buffer.data(ik + 163);
    const auto *ik_164 = buffer.data(ik + 164);
    const auto *ik_165 = buffer.data(ik + 165);
    const auto *ik_166 = buffer.data(ik + 166);
    const auto *ik_167 = buffer.data(ik + 167);
    const auto *ik_168 = buffer.data(ik + 168);
    const auto *ik_169 = buffer.data(ik + 169);
    const auto *ik_170 = buffer.data(ik + 170);
    const auto *ik_171 = buffer.data(ik + 171);
    const auto *ik_172 = buffer.data(ik + 172);
    const auto *ik_173 = buffer.data(ik + 173);
    const auto *ik_174 = buffer.data(ik + 174);
    const auto *ik_175 = buffer.data(ik + 175);
    const auto *ik_176 = buffer.data(ik + 176);
    const auto *ik_177 = buffer.data(ik + 177);
    const auto *ik_178 = buffer.data(ik + 178);
    const auto *ik_179 = buffer.data(ik + 179);
    const auto *ik_180 = buffer.data(ik + 180);
    const auto *ik_181 = buffer.data(ik + 181);
    const auto *ik_182 = buffer.data(ik + 182);
    const auto *ik_183 = buffer.data(ik + 183);
    const auto *ik_184 = buffer.data(ik + 184);
    const auto *ik_185 = buffer.data(ik + 185);
    const auto *ik_186 = buffer.data(ik + 186);
    const auto *ik_187 = buffer.data(ik + 187);
    const auto *ik_188 = buffer.data(ik + 188);
    const auto *ik_189 = buffer.data(ik + 189);
    const auto *ik_190 = buffer.data(ik + 190);
    const auto *ik_191 = buffer.data(ik + 191);
    const auto *ik_192 = buffer.data(ik + 192);
    const auto *ik_193 = buffer.data(ik + 193);
    const auto *ik_194 = buffer.data(ik + 194);
    const auto *ik_195 = buffer.data(ik + 195);
    const auto *ik_196 = buffer.data(ik + 196);
    const auto *ik_197 = buffer.data(ik + 197);
    const auto *ik_198 = buffer.data(ik + 198);
    const auto *ik_199 = buffer.data(ik + 199);
    const auto *ik_200 = buffer.data(ik + 200);
    const auto *ik_201 = buffer.data(ik + 201);
    const auto *ik_202 = buffer.data(ik + 202);
    const auto *ik_203 = buffer.data(ik + 203);
    const auto *ik_204 = buffer.data(ik + 204);
    const auto *ik_205 = buffer.data(ik + 205);
    const auto *ik_206 = buffer.data(ik + 206);
    const auto *ik_207 = buffer.data(ik + 207);
    const auto *ik_208 = buffer.data(ik + 208);
    const auto *ik_209 = buffer.data(ik + 209);
    const auto *ik_210 = buffer.data(ik + 210);
    const auto *ik_211 = buffer.data(ik + 211);
    const auto *ik_212 = buffer.data(ik + 212);
    const auto *ik_213 = buffer.data(ik + 213);
    const auto *ik_214 = buffer.data(ik + 214);
    const auto *ik_215 = buffer.data(ik + 215);
    const auto *ik_216 = buffer.data(ik + 216);
    const auto *ik_217 = buffer.data(ik + 217);
    const auto *ik_218 = buffer.data(ik + 218);
    const auto *ik_219 = buffer.data(ik + 219);
    const auto *ik_220 = buffer.data(ik + 220);
    const auto *ik_221 = buffer.data(ik + 221);
    const auto *ik_222 = buffer.data(ik + 222);
    const auto *ik_223 = buffer.data(ik + 223);
    const auto *ik_224 = buffer.data(ik + 224);
    const auto *ik_225 = buffer.data(ik + 225);
    const auto *ik_226 = buffer.data(ik + 226);
    const auto *ik_227 = buffer.data(ik + 227);
    const auto *ik_228 = buffer.data(ik + 228);
    const auto *ik_229 = buffer.data(ik + 229);
    const auto *ik_230 = buffer.data(ik + 230);
    const auto *ik_231 = buffer.data(ik + 231);
    const auto *ik_232 = buffer.data(ik + 232);
    const auto *ik_233 = buffer.data(ik + 233);
    const auto *ik_234 = buffer.data(ik + 234);
    const auto *ik_235 = buffer.data(ik + 235);
    const auto *ik_236 = buffer.data(ik + 236);
    const auto *ik_237 = buffer.data(ik + 237);
    const auto *ik_238 = buffer.data(ik + 238);
    const auto *ik_239 = buffer.data(ik + 239);
    const auto *ik_240 = buffer.data(ik + 240);
    const auto *ik_241 = buffer.data(ik + 241);
    const auto *ik_242 = buffer.data(ik + 242);
    const auto *ik_243 = buffer.data(ik + 243);
    const auto *ik_244 = buffer.data(ik + 244);
    const auto *ik_245 = buffer.data(ik + 245);
    const auto *ik_246 = buffer.data(ik + 246);
    const auto *ik_247 = buffer.data(ik + 247);
    const auto *ik_248 = buffer.data(ik + 248);
    const auto *ik_249 = buffer.data(ik + 249);
    const auto *ik_250 = buffer.data(ik + 250);
    const auto *ik_251 = buffer.data(ik + 251);
    const auto *ik_252 = buffer.data(ik + 252);
    const auto *ik_253 = buffer.data(ik + 253);
    const auto *ik_254 = buffer.data(ik + 254);
    const auto *ik_255 = buffer.data(ik + 255);
    const auto *ik_256 = buffer.data(ik + 256);
    const auto *ik_257 = buffer.data(ik + 257);
    const auto *ik_258 = buffer.data(ik + 258);
    const auto *ik_259 = buffer.data(ik + 259);
    const auto *ik_260 = buffer.data(ik + 260);
    const auto *ik_261 = buffer.data(ik + 261);
    const auto *ik_262 = buffer.data(ik + 262);
    const auto *ik_263 = buffer.data(ik + 263);
    const auto *ik_264 = buffer.data(ik + 264);
    const auto *ik_265 = buffer.data(ik + 265);
    const auto *ik_266 = buffer.data(ik + 266);
    const auto *ik_267 = buffer.data(ik + 267);
    const auto *ik_268 = buffer.data(ik + 268);
    const auto *ik_269 = buffer.data(ik + 269);
    const auto *ik_270 = buffer.data(ik + 270);
    const auto *ik_271 = buffer.data(ik + 271);
    const auto *ik_272 = buffer.data(ik + 272);
    const auto *ik_273 = buffer.data(ik + 273);
    const auto *ik_274 = buffer.data(ik + 274);
    const auto *ik_275 = buffer.data(ik + 275);
    const auto *ik_276 = buffer.data(ik + 276);
    const auto *ik_277 = buffer.data(ik + 277);
    const auto *ik_278 = buffer.data(ik + 278);
    const auto *ik_279 = buffer.data(ik + 279);
    const auto *ik_280 = buffer.data(ik + 280);
    const auto *ik_281 = buffer.data(ik + 281);
    const auto *ik_282 = buffer.data(ik + 282);
    const auto *ik_283 = buffer.data(ik + 283);
    const auto *ik_284 = buffer.data(ik + 284);
    const auto *ik_285 = buffer.data(ik + 285);
    const auto *ik_286 = buffer.data(ik + 286);
    const auto *ik_287 = buffer.data(ik + 287);
    const auto *ik_288 = buffer.data(ik + 288);
    const auto *ik_289 = buffer.data(ik + 289);
    const auto *ik_290 = buffer.data(ik + 290);
    const auto *ik_291 = buffer.data(ik + 291);
    const auto *ik_292 = buffer.data(ik + 292);
    const auto *ik_293 = buffer.data(ik + 293);
    const auto *ik_294 = buffer.data(ik + 294);
    const auto *ik_295 = buffer.data(ik + 295);
    const auto *ik_296 = buffer.data(ik + 296);
    const auto *ik_297 = buffer.data(ik + 297);
    const auto *ik_298 = buffer.data(ik + 298);
    const auto *ik_299 = buffer.data(ik + 299);

    const auto *lk_150 = buffer.data(lk + 150);
    const auto *lk_151 = buffer.data(lk + 151);
    const auto *lk_152 = buffer.data(lk + 152);
    const auto *lk_153 = buffer.data(lk + 153);
    const auto *lk_154 = buffer.data(lk + 154);
    const auto *lk_155 = buffer.data(lk + 155);
    const auto *lk_156 = buffer.data(lk + 156);
    const auto *lk_157 = buffer.data(lk + 157);
    const auto *lk_158 = buffer.data(lk + 158);
    const auto *lk_159 = buffer.data(lk + 159);
    const auto *lk_160 = buffer.data(lk + 160);
    const auto *lk_161 = buffer.data(lk + 161);
    const auto *lk_162 = buffer.data(lk + 162);
    const auto *lk_163 = buffer.data(lk + 163);
    const auto *lk_164 = buffer.data(lk + 164);
    const auto *lk_165 = buffer.data(lk + 165);
    const auto *lk_166 = buffer.data(lk + 166);
    const auto *lk_167 = buffer.data(lk + 167);
    const auto *lk_168 = buffer.data(lk + 168);
    const auto *lk_169 = buffer.data(lk + 169);
    const auto *lk_170 = buffer.data(lk + 170);
    const auto *lk_171 = buffer.data(lk + 171);
    const auto *lk_172 = buffer.data(lk + 172);
    const auto *lk_173 = buffer.data(lk + 173);
    const auto *lk_174 = buffer.data(lk + 174);
    const auto *lk_175 = buffer.data(lk + 175);
    const auto *lk_176 = buffer.data(lk + 176);
    const auto *lk_177 = buffer.data(lk + 177);
    const auto *lk_178 = buffer.data(lk + 178);
    const auto *lk_179 = buffer.data(lk + 179);
    const auto *lk_180 = buffer.data(lk + 180);
    const auto *lk_181 = buffer.data(lk + 181);
    const auto *lk_182 = buffer.data(lk + 182);
    const auto *lk_183 = buffer.data(lk + 183);
    const auto *lk_184 = buffer.data(lk + 184);
    const auto *lk_185 = buffer.data(lk + 185);
    const auto *lk_186 = buffer.data(lk + 186);
    const auto *lk_187 = buffer.data(lk + 187);
    const auto *lk_188 = buffer.data(lk + 188);
    const auto *lk_189 = buffer.data(lk + 189);
    const auto *lk_190 = buffer.data(lk + 190);
    const auto *lk_191 = buffer.data(lk + 191);
    const auto *lk_192 = buffer.data(lk + 192);
    const auto *lk_193 = buffer.data(lk + 193);
    const auto *lk_194 = buffer.data(lk + 194);
    const auto *lk_195 = buffer.data(lk + 195);
    const auto *lk_196 = buffer.data(lk + 196);
    const auto *lk_197 = buffer.data(lk + 197);
    const auto *lk_198 = buffer.data(lk + 198);
    const auto *lk_199 = buffer.data(lk + 199);
    const auto *lk_200 = buffer.data(lk + 200);
    const auto *lk_201 = buffer.data(lk + 201);
    const auto *lk_202 = buffer.data(lk + 202);
    const auto *lk_203 = buffer.data(lk + 203);
    const auto *lk_204 = buffer.data(lk + 204);
    const auto *lk_205 = buffer.data(lk + 205);
    const auto *lk_206 = buffer.data(lk + 206);
    const auto *lk_207 = buffer.data(lk + 207);
    const auto *lk_208 = buffer.data(lk + 208);
    const auto *lk_209 = buffer.data(lk + 209);
    const auto *lk_210 = buffer.data(lk + 210);
    const auto *lk_211 = buffer.data(lk + 211);
    const auto *lk_212 = buffer.data(lk + 212);
    const auto *lk_213 = buffer.data(lk + 213);
    const auto *lk_214 = buffer.data(lk + 214);
    const auto *lk_215 = buffer.data(lk + 215);
    const auto *lk_216 = buffer.data(lk + 216);
    const auto *lk_217 = buffer.data(lk + 217);
    const auto *lk_218 = buffer.data(lk + 218);
    const auto *lk_219 = buffer.data(lk + 219);
    const auto *lk_220 = buffer.data(lk + 220);
    const auto *lk_221 = buffer.data(lk + 221);
    const auto *lk_222 = buffer.data(lk + 222);
    const auto *lk_223 = buffer.data(lk + 223);
    const auto *lk_224 = buffer.data(lk + 224);
    const auto *lk_225 = buffer.data(lk + 225);
    const auto *lk_226 = buffer.data(lk + 226);
    const auto *lk_227 = buffer.data(lk + 227);
    const auto *lk_228 = buffer.data(lk + 228);
    const auto *lk_229 = buffer.data(lk + 229);
    const auto *lk_230 = buffer.data(lk + 230);
    const auto *lk_231 = buffer.data(lk + 231);
    const auto *lk_232 = buffer.data(lk + 232);
    const auto *lk_233 = buffer.data(lk + 233);
    const auto *lk_234 = buffer.data(lk + 234);
    const auto *lk_235 = buffer.data(lk + 235);
    const auto *lk_236 = buffer.data(lk + 236);
    const auto *lk_237 = buffer.data(lk + 237);
    const auto *lk_238 = buffer.data(lk + 238);
    const auto *lk_239 = buffer.data(lk + 239);
    const auto *lk_240 = buffer.data(lk + 240);
    const auto *lk_241 = buffer.data(lk + 241);
    const auto *lk_242 = buffer.data(lk + 242);
    const auto *lk_243 = buffer.data(lk + 243);
    const auto *lk_244 = buffer.data(lk + 244);
    const auto *lk_245 = buffer.data(lk + 245);
    const auto *lk_246 = buffer.data(lk + 246);
    const auto *lk_247 = buffer.data(lk + 247);
    const auto *lk_248 = buffer.data(lk + 248);
    const auto *lk_249 = buffer.data(lk + 249);
    const auto *lk_250 = buffer.data(lk + 250);
    const auto *lk_251 = buffer.data(lk + 251);
    const auto *lk_252 = buffer.data(lk + 252);
    const auto *lk_253 = buffer.data(lk + 253);
    const auto *lk_254 = buffer.data(lk + 254);
    const auto *lk_255 = buffer.data(lk + 255);
    const auto *lk_256 = buffer.data(lk + 256);
    const auto *lk_257 = buffer.data(lk + 257);
    const auto *lk_258 = buffer.data(lk + 258);
    const auto *lk_259 = buffer.data(lk + 259);
    const auto *lk_260 = buffer.data(lk + 260);
    const auto *lk_261 = buffer.data(lk + 261);
    const auto *lk_262 = buffer.data(lk + 262);
    const auto *lk_263 = buffer.data(lk + 263);
    const auto *lk_264 = buffer.data(lk + 264);
    const auto *lk_265 = buffer.data(lk + 265);
    const auto *lk_266 = buffer.data(lk + 266);
    const auto *lk_267 = buffer.data(lk + 267);
    const auto *lk_268 = buffer.data(lk + 268);
    const auto *lk_269 = buffer.data(lk + 269);
    const auto *lk_270 = buffer.data(lk + 270);
    const auto *lk_271 = buffer.data(lk + 271);
    const auto *lk_272 = buffer.data(lk + 272);
    const auto *lk_273 = buffer.data(lk + 273);
    const auto *lk_274 = buffer.data(lk + 274);
    const auto *lk_275 = buffer.data(lk + 275);
    const auto *lk_276 = buffer.data(lk + 276);
    const auto *lk_277 = buffer.data(lk + 277);
    const auto *lk_278 = buffer.data(lk + 278);
    const auto *lk_279 = buffer.data(lk + 279);
    const auto *lk_280 = buffer.data(lk + 280);
    const auto *lk_281 = buffer.data(lk + 281);
    const auto *lk_282 = buffer.data(lk + 282);
    const auto *lk_283 = buffer.data(lk + 283);
    const auto *lk_284 = buffer.data(lk + 284);
    const auto *lk_285 = buffer.data(lk + 285);
    const auto *lk_286 = buffer.data(lk + 286);
    const auto *lk_287 = buffer.data(lk + 287);
    const auto *lk_288 = buffer.data(lk + 288);
    const auto *lk_289 = buffer.data(lk + 289);
    const auto *lk_290 = buffer.data(lk + 290);
    const auto *lk_291 = buffer.data(lk + 291);
    const auto *lk_292 = buffer.data(lk + 292);
    const auto *lk_293 = buffer.data(lk + 293);
    const auto *lk_294 = buffer.data(lk + 294);
    const auto *lk_295 = buffer.data(lk + 295);
    const auto *lk_296 = buffer.data(lk + 296);
    const auto *lk_297 = buffer.data(lk + 297);
    const auto *lk_298 = buffer.data(lk + 298);
    const auto *lk_299 = buffer.data(lk + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ik_150, ik_151, ik_152, ik_153, \
                         ik_154, lk_150, lk_151, lk_152, lk_153, \
                         lk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -5.0 * ik_150[k]
                   + f_0 * lk_150[k];

        t_151[k] = -5.0 * ik_151[k]
                   + f_0 * lk_151[k];

        t_152[k] = -5.0 * ik_152[k]
                   + f_0 * lk_152[k];

        t_153[k] = -5.0 * ik_153[k]
                   + f_0 * lk_153[k];

        t_154[k] = -5.0 * ik_154[k]
                   + f_0 * lk_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ik_155, ik_156, ik_157, ik_158, \
                         ik_159, lk_155, lk_156, lk_157, lk_158, \
                         lk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -5.0 * ik_155[k]
                   + f_0 * lk_155[k];

        t_156[k] = -5.0 * ik_156[k]
                   + f_0 * lk_156[k];

        t_157[k] = -5.0 * ik_157[k]
                   + f_0 * lk_157[k];

        t_158[k] = -5.0 * ik_158[k]
                   + f_0 * lk_158[k];

        t_159[k] = -5.0 * ik_159[k]
                   + f_0 * lk_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ik_160, ik_161, ik_162, ik_163, \
                         ik_164, lk_160, lk_161, lk_162, lk_163, \
                         lk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -5.0 * ik_160[k]
                   + f_0 * lk_160[k];

        t_161[k] = -5.0 * ik_161[k]
                   + f_0 * lk_161[k];

        t_162[k] = -5.0 * ik_162[k]
                   + f_0 * lk_162[k];

        t_163[k] = -5.0 * ik_163[k]
                   + f_0 * lk_163[k];

        t_164[k] = -5.0 * ik_164[k]
                   + f_0 * lk_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ik_165, ik_166, ik_167, ik_168, \
                         ik_169, lk_165, lk_166, lk_167, lk_168, \
                         lk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -5.0 * ik_165[k]
                   + f_0 * lk_165[k];

        t_166[k] = -5.0 * ik_166[k]
                   + f_0 * lk_166[k];

        t_167[k] = -5.0 * ik_167[k]
                   + f_0 * lk_167[k];

        t_168[k] = -5.0 * ik_168[k]
                   + f_0 * lk_168[k];

        t_169[k] = -5.0 * ik_169[k]
                   + f_0 * lk_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ik_170, ik_171, ik_172, ik_173, \
                         ik_174, lk_170, lk_171, lk_172, lk_173, \
                         lk_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -5.0 * ik_170[k]
                   + f_0 * lk_170[k];

        t_171[k] = -5.0 * ik_171[k]
                   + f_0 * lk_171[k];

        t_172[k] = -5.0 * ik_172[k]
                   + f_0 * lk_172[k];

        t_173[k] = -5.0 * ik_173[k]
                   + f_0 * lk_173[k];

        t_174[k] = -5.0 * ik_174[k]
                   + f_0 * lk_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ik_175, ik_176, ik_177, ik_178, \
                         ik_179, lk_175, lk_176, lk_177, lk_178, \
                         lk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -5.0 * ik_175[k]
                   + f_0 * lk_175[k];

        t_176[k] = -5.0 * ik_176[k]
                   + f_0 * lk_176[k];

        t_177[k] = -5.0 * ik_177[k]
                   + f_0 * lk_177[k];

        t_178[k] = -5.0 * ik_178[k]
                   + f_0 * lk_178[k];

        t_179[k] = -5.0 * ik_179[k]
                   + f_0 * lk_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ik_180, ik_181, ik_182, ik_183, \
                         ik_184, lk_180, lk_181, lk_182, lk_183, \
                         lk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -5.0 * ik_180[k]
                   + f_0 * lk_180[k];

        t_181[k] = -5.0 * ik_181[k]
                   + f_0 * lk_181[k];

        t_182[k] = -5.0 * ik_182[k]
                   + f_0 * lk_182[k];

        t_183[k] = -5.0 * ik_183[k]
                   + f_0 * lk_183[k];

        t_184[k] = -5.0 * ik_184[k]
                   + f_0 * lk_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ik_185, ik_186, ik_187, ik_188, \
                         ik_189, lk_185, lk_186, lk_187, lk_188, \
                         lk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -5.0 * ik_185[k]
                   + f_0 * lk_185[k];

        t_186[k] = -5.0 * ik_186[k]
                   + f_0 * lk_186[k];

        t_187[k] = -5.0 * ik_187[k]
                   + f_0 * lk_187[k];

        t_188[k] = -5.0 * ik_188[k]
                   + f_0 * lk_188[k];

        t_189[k] = -5.0 * ik_189[k]
                   + f_0 * lk_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ik_190, ik_191, ik_192, ik_193, \
                         ik_194, lk_190, lk_191, lk_192, lk_193, \
                         lk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -5.0 * ik_190[k]
                   + f_0 * lk_190[k];

        t_191[k] = -5.0 * ik_191[k]
                   + f_0 * lk_191[k];

        t_192[k] = -5.0 * ik_192[k]
                   + f_0 * lk_192[k];

        t_193[k] = -5.0 * ik_193[k]
                   + f_0 * lk_193[k];

        t_194[k] = -5.0 * ik_194[k]
                   + f_0 * lk_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ik_195, ik_196, ik_197, ik_198, \
                         ik_199, lk_195, lk_196, lk_197, lk_198, \
                         lk_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -5.0 * ik_195[k]
                   + f_0 * lk_195[k];

        t_196[k] = -5.0 * ik_196[k]
                   + f_0 * lk_196[k];

        t_197[k] = -5.0 * ik_197[k]
                   + f_0 * lk_197[k];

        t_198[k] = -5.0 * ik_198[k]
                   + f_0 * lk_198[k];

        t_199[k] = -5.0 * ik_199[k]
                   + f_0 * lk_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ik_200, ik_201, ik_202, ik_203, \
                         ik_204, lk_200, lk_201, lk_202, lk_203, \
                         lk_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -5.0 * ik_200[k]
                   + f_0 * lk_200[k];

        t_201[k] = -5.0 * ik_201[k]
                   + f_0 * lk_201[k];

        t_202[k] = -5.0 * ik_202[k]
                   + f_0 * lk_202[k];

        t_203[k] = -5.0 * ik_203[k]
                   + f_0 * lk_203[k];

        t_204[k] = -5.0 * ik_204[k]
                   + f_0 * lk_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ik_205, ik_206, ik_207, ik_208, \
                         ik_209, lk_205, lk_206, lk_207, lk_208, \
                         lk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -5.0 * ik_205[k]
                   + f_0 * lk_205[k];

        t_206[k] = -5.0 * ik_206[k]
                   + f_0 * lk_206[k];

        t_207[k] = -5.0 * ik_207[k]
                   + f_0 * lk_207[k];

        t_208[k] = -5.0 * ik_208[k]
                   + f_0 * lk_208[k];

        t_209[k] = -5.0 * ik_209[k]
                   + f_0 * lk_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ik_210, ik_211, ik_212, ik_213, \
                         ik_214, lk_210, lk_211, lk_212, lk_213, \
                         lk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -5.0 * ik_210[k]
                   + f_0 * lk_210[k];

        t_211[k] = -5.0 * ik_211[k]
                   + f_0 * lk_211[k];

        t_212[k] = -5.0 * ik_212[k]
                   + f_0 * lk_212[k];

        t_213[k] = -5.0 * ik_213[k]
                   + f_0 * lk_213[k];

        t_214[k] = -5.0 * ik_214[k]
                   + f_0 * lk_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ik_215, ik_216, ik_217, ik_218, \
                         ik_219, lk_215, lk_216, lk_217, lk_218, \
                         lk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -5.0 * ik_215[k]
                   + f_0 * lk_215[k];

        t_216[k] = -4.0 * ik_216[k]
                   + f_0 * lk_216[k];

        t_217[k] = -4.0 * ik_217[k]
                   + f_0 * lk_217[k];

        t_218[k] = -4.0 * ik_218[k]
                   + f_0 * lk_218[k];

        t_219[k] = -4.0 * ik_219[k]
                   + f_0 * lk_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ik_220, ik_221, ik_222, ik_223, \
                         ik_224, lk_220, lk_221, lk_222, lk_223, \
                         lk_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -4.0 * ik_220[k]
                   + f_0 * lk_220[k];

        t_221[k] = -4.0 * ik_221[k]
                   + f_0 * lk_221[k];

        t_222[k] = -4.0 * ik_222[k]
                   + f_0 * lk_222[k];

        t_223[k] = -4.0 * ik_223[k]
                   + f_0 * lk_223[k];

        t_224[k] = -4.0 * ik_224[k]
                   + f_0 * lk_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ik_225, ik_226, ik_227, ik_228, \
                         ik_229, lk_225, lk_226, lk_227, lk_228, \
                         lk_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -4.0 * ik_225[k]
                   + f_0 * lk_225[k];

        t_226[k] = -4.0 * ik_226[k]
                   + f_0 * lk_226[k];

        t_227[k] = -4.0 * ik_227[k]
                   + f_0 * lk_227[k];

        t_228[k] = -4.0 * ik_228[k]
                   + f_0 * lk_228[k];

        t_229[k] = -4.0 * ik_229[k]
                   + f_0 * lk_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ik_230, ik_231, ik_232, ik_233, \
                         ik_234, lk_230, lk_231, lk_232, lk_233, \
                         lk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -4.0 * ik_230[k]
                   + f_0 * lk_230[k];

        t_231[k] = -4.0 * ik_231[k]
                   + f_0 * lk_231[k];

        t_232[k] = -4.0 * ik_232[k]
                   + f_0 * lk_232[k];

        t_233[k] = -4.0 * ik_233[k]
                   + f_0 * lk_233[k];

        t_234[k] = -4.0 * ik_234[k]
                   + f_0 * lk_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ik_235, ik_236, ik_237, ik_238, \
                         ik_239, lk_235, lk_236, lk_237, lk_238, \
                         lk_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -4.0 * ik_235[k]
                   + f_0 * lk_235[k];

        t_236[k] = -4.0 * ik_236[k]
                   + f_0 * lk_236[k];

        t_237[k] = -4.0 * ik_237[k]
                   + f_0 * lk_237[k];

        t_238[k] = -4.0 * ik_238[k]
                   + f_0 * lk_238[k];

        t_239[k] = -4.0 * ik_239[k]
                   + f_0 * lk_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ik_240, ik_241, ik_242, ik_243, \
                         ik_244, lk_240, lk_241, lk_242, lk_243, \
                         lk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -4.0 * ik_240[k]
                   + f_0 * lk_240[k];

        t_241[k] = -4.0 * ik_241[k]
                   + f_0 * lk_241[k];

        t_242[k] = -4.0 * ik_242[k]
                   + f_0 * lk_242[k];

        t_243[k] = -4.0 * ik_243[k]
                   + f_0 * lk_243[k];

        t_244[k] = -4.0 * ik_244[k]
                   + f_0 * lk_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ik_245, ik_246, ik_247, ik_248, \
                         ik_249, lk_245, lk_246, lk_247, lk_248, \
                         lk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -4.0 * ik_245[k]
                   + f_0 * lk_245[k];

        t_246[k] = -4.0 * ik_246[k]
                   + f_0 * lk_246[k];

        t_247[k] = -4.0 * ik_247[k]
                   + f_0 * lk_247[k];

        t_248[k] = -4.0 * ik_248[k]
                   + f_0 * lk_248[k];

        t_249[k] = -4.0 * ik_249[k]
                   + f_0 * lk_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ik_250, ik_251, ik_252, ik_253, \
                         ik_254, lk_250, lk_251, lk_252, lk_253, \
                         lk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -4.0 * ik_250[k]
                   + f_0 * lk_250[k];

        t_251[k] = -4.0 * ik_251[k]
                   + f_0 * lk_251[k];

        t_252[k] = -4.0 * ik_252[k]
                   + f_0 * lk_252[k];

        t_253[k] = -4.0 * ik_253[k]
                   + f_0 * lk_253[k];

        t_254[k] = -4.0 * ik_254[k]
                   + f_0 * lk_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ik_255, ik_256, ik_257, ik_258, \
                         ik_259, lk_255, lk_256, lk_257, lk_258, \
                         lk_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -4.0 * ik_255[k]
                   + f_0 * lk_255[k];

        t_256[k] = -4.0 * ik_256[k]
                   + f_0 * lk_256[k];

        t_257[k] = -4.0 * ik_257[k]
                   + f_0 * lk_257[k];

        t_258[k] = -4.0 * ik_258[k]
                   + f_0 * lk_258[k];

        t_259[k] = -4.0 * ik_259[k]
                   + f_0 * lk_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ik_260, ik_261, ik_262, ik_263, \
                         ik_264, lk_260, lk_261, lk_262, lk_263, \
                         lk_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -4.0 * ik_260[k]
                   + f_0 * lk_260[k];

        t_261[k] = -4.0 * ik_261[k]
                   + f_0 * lk_261[k];

        t_262[k] = -4.0 * ik_262[k]
                   + f_0 * lk_262[k];

        t_263[k] = -4.0 * ik_263[k]
                   + f_0 * lk_263[k];

        t_264[k] = -4.0 * ik_264[k]
                   + f_0 * lk_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ik_265, ik_266, ik_267, ik_268, \
                         ik_269, lk_265, lk_266, lk_267, lk_268, \
                         lk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -4.0 * ik_265[k]
                   + f_0 * lk_265[k];

        t_266[k] = -4.0 * ik_266[k]
                   + f_0 * lk_266[k];

        t_267[k] = -4.0 * ik_267[k]
                   + f_0 * lk_267[k];

        t_268[k] = -4.0 * ik_268[k]
                   + f_0 * lk_268[k];

        t_269[k] = -4.0 * ik_269[k]
                   + f_0 * lk_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ik_270, ik_271, ik_272, ik_273, \
                         ik_274, lk_270, lk_271, lk_272, lk_273, \
                         lk_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -4.0 * ik_270[k]
                   + f_0 * lk_270[k];

        t_271[k] = -4.0 * ik_271[k]
                   + f_0 * lk_271[k];

        t_272[k] = -4.0 * ik_272[k]
                   + f_0 * lk_272[k];

        t_273[k] = -4.0 * ik_273[k]
                   + f_0 * lk_273[k];

        t_274[k] = -4.0 * ik_274[k]
                   + f_0 * lk_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ik_275, ik_276, ik_277, ik_278, \
                         ik_279, lk_275, lk_276, lk_277, lk_278, \
                         lk_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -4.0 * ik_275[k]
                   + f_0 * lk_275[k];

        t_276[k] = -4.0 * ik_276[k]
                   + f_0 * lk_276[k];

        t_277[k] = -4.0 * ik_277[k]
                   + f_0 * lk_277[k];

        t_278[k] = -4.0 * ik_278[k]
                   + f_0 * lk_278[k];

        t_279[k] = -4.0 * ik_279[k]
                   + f_0 * lk_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ik_280, ik_281, ik_282, ik_283, \
                         ik_284, lk_280, lk_281, lk_282, lk_283, \
                         lk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -4.0 * ik_280[k]
                   + f_0 * lk_280[k];

        t_281[k] = -4.0 * ik_281[k]
                   + f_0 * lk_281[k];

        t_282[k] = -4.0 * ik_282[k]
                   + f_0 * lk_282[k];

        t_283[k] = -4.0 * ik_283[k]
                   + f_0 * lk_283[k];

        t_284[k] = -4.0 * ik_284[k]
                   + f_0 * lk_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ik_285, ik_286, ik_287, ik_288, \
                         ik_289, lk_285, lk_286, lk_287, lk_288, \
                         lk_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -4.0 * ik_285[k]
                   + f_0 * lk_285[k];

        t_286[k] = -4.0 * ik_286[k]
                   + f_0 * lk_286[k];

        t_287[k] = -4.0 * ik_287[k]
                   + f_0 * lk_287[k];

        t_288[k] = -4.0 * ik_288[k]
                   + f_0 * lk_288[k];

        t_289[k] = -4.0 * ik_289[k]
                   + f_0 * lk_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ik_290, ik_291, ik_292, ik_293, \
                         ik_294, lk_290, lk_291, lk_292, lk_293, \
                         lk_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -4.0 * ik_290[k]
                   + f_0 * lk_290[k];

        t_291[k] = -4.0 * ik_291[k]
                   + f_0 * lk_291[k];

        t_292[k] = -4.0 * ik_292[k]
                   + f_0 * lk_292[k];

        t_293[k] = -4.0 * ik_293[k]
                   + f_0 * lk_293[k];

        t_294[k] = -4.0 * ik_294[k]
                   + f_0 * lk_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ik_295, ik_296, ik_297, ik_298, \
                         ik_299, lk_295, lk_296, lk_297, lk_298, \
                         lk_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * ik_295[k]
                   + f_0 * lk_295[k];

        t_296[k] = -4.0 * ik_296[k]
                   + f_0 * lk_296[k];

        t_297[k] = -4.0 * ik_297[k]
                   + f_0 * lk_297[k];

        t_298[k] = -4.0 * ik_298[k]
                   + f_0 * lk_298[k];

        t_299[k] = -4.0 * ik_299[k]
                   + f_0 * lk_299[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_300 = buffer.data(ik + 300);
    const auto *ik_301 = buffer.data(ik + 301);
    const auto *ik_302 = buffer.data(ik + 302);
    const auto *ik_303 = buffer.data(ik + 303);
    const auto *ik_304 = buffer.data(ik + 304);
    const auto *ik_305 = buffer.data(ik + 305);
    const auto *ik_306 = buffer.data(ik + 306);
    const auto *ik_307 = buffer.data(ik + 307);
    const auto *ik_308 = buffer.data(ik + 308);
    const auto *ik_309 = buffer.data(ik + 309);
    const auto *ik_310 = buffer.data(ik + 310);
    const auto *ik_311 = buffer.data(ik + 311);
    const auto *ik_312 = buffer.data(ik + 312);
    const auto *ik_313 = buffer.data(ik + 313);
    const auto *ik_314 = buffer.data(ik + 314);
    const auto *ik_315 = buffer.data(ik + 315);
    const auto *ik_316 = buffer.data(ik + 316);
    const auto *ik_317 = buffer.data(ik + 317);
    const auto *ik_318 = buffer.data(ik + 318);
    const auto *ik_319 = buffer.data(ik + 319);
    const auto *ik_320 = buffer.data(ik + 320);
    const auto *ik_321 = buffer.data(ik + 321);
    const auto *ik_322 = buffer.data(ik + 322);
    const auto *ik_323 = buffer.data(ik + 323);
    const auto *ik_324 = buffer.data(ik + 324);
    const auto *ik_325 = buffer.data(ik + 325);
    const auto *ik_326 = buffer.data(ik + 326);
    const auto *ik_327 = buffer.data(ik + 327);
    const auto *ik_328 = buffer.data(ik + 328);
    const auto *ik_329 = buffer.data(ik + 329);
    const auto *ik_330 = buffer.data(ik + 330);
    const auto *ik_331 = buffer.data(ik + 331);
    const auto *ik_332 = buffer.data(ik + 332);
    const auto *ik_333 = buffer.data(ik + 333);
    const auto *ik_334 = buffer.data(ik + 334);
    const auto *ik_335 = buffer.data(ik + 335);
    const auto *ik_336 = buffer.data(ik + 336);
    const auto *ik_337 = buffer.data(ik + 337);
    const auto *ik_338 = buffer.data(ik + 338);
    const auto *ik_339 = buffer.data(ik + 339);
    const auto *ik_340 = buffer.data(ik + 340);
    const auto *ik_341 = buffer.data(ik + 341);
    const auto *ik_342 = buffer.data(ik + 342);
    const auto *ik_343 = buffer.data(ik + 343);
    const auto *ik_344 = buffer.data(ik + 344);
    const auto *ik_345 = buffer.data(ik + 345);
    const auto *ik_346 = buffer.data(ik + 346);
    const auto *ik_347 = buffer.data(ik + 347);
    const auto *ik_348 = buffer.data(ik + 348);
    const auto *ik_349 = buffer.data(ik + 349);
    const auto *ik_350 = buffer.data(ik + 350);
    const auto *ik_351 = buffer.data(ik + 351);
    const auto *ik_352 = buffer.data(ik + 352);
    const auto *ik_353 = buffer.data(ik + 353);
    const auto *ik_354 = buffer.data(ik + 354);
    const auto *ik_355 = buffer.data(ik + 355);
    const auto *ik_356 = buffer.data(ik + 356);
    const auto *ik_357 = buffer.data(ik + 357);
    const auto *ik_358 = buffer.data(ik + 358);
    const auto *ik_359 = buffer.data(ik + 359);
    const auto *ik_360 = buffer.data(ik + 360);
    const auto *ik_361 = buffer.data(ik + 361);
    const auto *ik_362 = buffer.data(ik + 362);
    const auto *ik_363 = buffer.data(ik + 363);
    const auto *ik_364 = buffer.data(ik + 364);
    const auto *ik_365 = buffer.data(ik + 365);
    const auto *ik_366 = buffer.data(ik + 366);
    const auto *ik_367 = buffer.data(ik + 367);
    const auto *ik_368 = buffer.data(ik + 368);
    const auto *ik_369 = buffer.data(ik + 369);
    const auto *ik_370 = buffer.data(ik + 370);
    const auto *ik_371 = buffer.data(ik + 371);
    const auto *ik_372 = buffer.data(ik + 372);
    const auto *ik_373 = buffer.data(ik + 373);
    const auto *ik_374 = buffer.data(ik + 374);
    const auto *ik_375 = buffer.data(ik + 375);
    const auto *ik_376 = buffer.data(ik + 376);
    const auto *ik_377 = buffer.data(ik + 377);
    const auto *ik_378 = buffer.data(ik + 378);
    const auto *ik_379 = buffer.data(ik + 379);
    const auto *ik_380 = buffer.data(ik + 380);
    const auto *ik_381 = buffer.data(ik + 381);
    const auto *ik_382 = buffer.data(ik + 382);
    const auto *ik_383 = buffer.data(ik + 383);
    const auto *ik_384 = buffer.data(ik + 384);
    const auto *ik_385 = buffer.data(ik + 385);
    const auto *ik_386 = buffer.data(ik + 386);
    const auto *ik_387 = buffer.data(ik + 387);
    const auto *ik_388 = buffer.data(ik + 388);
    const auto *ik_389 = buffer.data(ik + 389);
    const auto *ik_390 = buffer.data(ik + 390);
    const auto *ik_391 = buffer.data(ik + 391);
    const auto *ik_392 = buffer.data(ik + 392);
    const auto *ik_393 = buffer.data(ik + 393);
    const auto *ik_394 = buffer.data(ik + 394);
    const auto *ik_395 = buffer.data(ik + 395);
    const auto *ik_396 = buffer.data(ik + 396);
    const auto *ik_397 = buffer.data(ik + 397);
    const auto *ik_398 = buffer.data(ik + 398);
    const auto *ik_399 = buffer.data(ik + 399);
    const auto *ik_400 = buffer.data(ik + 400);
    const auto *ik_401 = buffer.data(ik + 401);
    const auto *ik_402 = buffer.data(ik + 402);
    const auto *ik_403 = buffer.data(ik + 403);
    const auto *ik_404 = buffer.data(ik + 404);
    const auto *ik_405 = buffer.data(ik + 405);
    const auto *ik_406 = buffer.data(ik + 406);
    const auto *ik_407 = buffer.data(ik + 407);
    const auto *ik_408 = buffer.data(ik + 408);
    const auto *ik_409 = buffer.data(ik + 409);
    const auto *ik_410 = buffer.data(ik + 410);
    const auto *ik_411 = buffer.data(ik + 411);
    const auto *ik_412 = buffer.data(ik + 412);
    const auto *ik_413 = buffer.data(ik + 413);
    const auto *ik_414 = buffer.data(ik + 414);
    const auto *ik_415 = buffer.data(ik + 415);
    const auto *ik_416 = buffer.data(ik + 416);
    const auto *ik_417 = buffer.data(ik + 417);
    const auto *ik_418 = buffer.data(ik + 418);
    const auto *ik_419 = buffer.data(ik + 419);
    const auto *ik_420 = buffer.data(ik + 420);
    const auto *ik_421 = buffer.data(ik + 421);
    const auto *ik_422 = buffer.data(ik + 422);
    const auto *ik_423 = buffer.data(ik + 423);
    const auto *ik_424 = buffer.data(ik + 424);
    const auto *ik_425 = buffer.data(ik + 425);
    const auto *ik_426 = buffer.data(ik + 426);
    const auto *ik_427 = buffer.data(ik + 427);
    const auto *ik_428 = buffer.data(ik + 428);
    const auto *ik_429 = buffer.data(ik + 429);
    const auto *ik_430 = buffer.data(ik + 430);
    const auto *ik_431 = buffer.data(ik + 431);
    const auto *ik_432 = buffer.data(ik + 432);
    const auto *ik_433 = buffer.data(ik + 433);
    const auto *ik_434 = buffer.data(ik + 434);
    const auto *ik_435 = buffer.data(ik + 435);
    const auto *ik_436 = buffer.data(ik + 436);
    const auto *ik_437 = buffer.data(ik + 437);
    const auto *ik_438 = buffer.data(ik + 438);
    const auto *ik_439 = buffer.data(ik + 439);
    const auto *ik_440 = buffer.data(ik + 440);
    const auto *ik_441 = buffer.data(ik + 441);
    const auto *ik_442 = buffer.data(ik + 442);
    const auto *ik_443 = buffer.data(ik + 443);
    const auto *ik_444 = buffer.data(ik + 444);
    const auto *ik_445 = buffer.data(ik + 445);
    const auto *ik_446 = buffer.data(ik + 446);
    const auto *ik_447 = buffer.data(ik + 447);
    const auto *ik_448 = buffer.data(ik + 448);
    const auto *ik_449 = buffer.data(ik + 449);

    const auto *lk_300 = buffer.data(lk + 300);
    const auto *lk_301 = buffer.data(lk + 301);
    const auto *lk_302 = buffer.data(lk + 302);
    const auto *lk_303 = buffer.data(lk + 303);
    const auto *lk_304 = buffer.data(lk + 304);
    const auto *lk_305 = buffer.data(lk + 305);
    const auto *lk_306 = buffer.data(lk + 306);
    const auto *lk_307 = buffer.data(lk + 307);
    const auto *lk_308 = buffer.data(lk + 308);
    const auto *lk_309 = buffer.data(lk + 309);
    const auto *lk_310 = buffer.data(lk + 310);
    const auto *lk_311 = buffer.data(lk + 311);
    const auto *lk_312 = buffer.data(lk + 312);
    const auto *lk_313 = buffer.data(lk + 313);
    const auto *lk_314 = buffer.data(lk + 314);
    const auto *lk_315 = buffer.data(lk + 315);
    const auto *lk_316 = buffer.data(lk + 316);
    const auto *lk_317 = buffer.data(lk + 317);
    const auto *lk_318 = buffer.data(lk + 318);
    const auto *lk_319 = buffer.data(lk + 319);
    const auto *lk_320 = buffer.data(lk + 320);
    const auto *lk_321 = buffer.data(lk + 321);
    const auto *lk_322 = buffer.data(lk + 322);
    const auto *lk_323 = buffer.data(lk + 323);
    const auto *lk_324 = buffer.data(lk + 324);
    const auto *lk_325 = buffer.data(lk + 325);
    const auto *lk_326 = buffer.data(lk + 326);
    const auto *lk_327 = buffer.data(lk + 327);
    const auto *lk_328 = buffer.data(lk + 328);
    const auto *lk_329 = buffer.data(lk + 329);
    const auto *lk_330 = buffer.data(lk + 330);
    const auto *lk_331 = buffer.data(lk + 331);
    const auto *lk_332 = buffer.data(lk + 332);
    const auto *lk_333 = buffer.data(lk + 333);
    const auto *lk_334 = buffer.data(lk + 334);
    const auto *lk_335 = buffer.data(lk + 335);
    const auto *lk_336 = buffer.data(lk + 336);
    const auto *lk_337 = buffer.data(lk + 337);
    const auto *lk_338 = buffer.data(lk + 338);
    const auto *lk_339 = buffer.data(lk + 339);
    const auto *lk_340 = buffer.data(lk + 340);
    const auto *lk_341 = buffer.data(lk + 341);
    const auto *lk_342 = buffer.data(lk + 342);
    const auto *lk_343 = buffer.data(lk + 343);
    const auto *lk_344 = buffer.data(lk + 344);
    const auto *lk_345 = buffer.data(lk + 345);
    const auto *lk_346 = buffer.data(lk + 346);
    const auto *lk_347 = buffer.data(lk + 347);
    const auto *lk_348 = buffer.data(lk + 348);
    const auto *lk_349 = buffer.data(lk + 349);
    const auto *lk_350 = buffer.data(lk + 350);
    const auto *lk_351 = buffer.data(lk + 351);
    const auto *lk_352 = buffer.data(lk + 352);
    const auto *lk_353 = buffer.data(lk + 353);
    const auto *lk_354 = buffer.data(lk + 354);
    const auto *lk_355 = buffer.data(lk + 355);
    const auto *lk_356 = buffer.data(lk + 356);
    const auto *lk_357 = buffer.data(lk + 357);
    const auto *lk_358 = buffer.data(lk + 358);
    const auto *lk_359 = buffer.data(lk + 359);
    const auto *lk_360 = buffer.data(lk + 360);
    const auto *lk_361 = buffer.data(lk + 361);
    const auto *lk_362 = buffer.data(lk + 362);
    const auto *lk_363 = buffer.data(lk + 363);
    const auto *lk_364 = buffer.data(lk + 364);
    const auto *lk_365 = buffer.data(lk + 365);
    const auto *lk_366 = buffer.data(lk + 366);
    const auto *lk_367 = buffer.data(lk + 367);
    const auto *lk_368 = buffer.data(lk + 368);
    const auto *lk_369 = buffer.data(lk + 369);
    const auto *lk_370 = buffer.data(lk + 370);
    const auto *lk_371 = buffer.data(lk + 371);
    const auto *lk_372 = buffer.data(lk + 372);
    const auto *lk_373 = buffer.data(lk + 373);
    const auto *lk_374 = buffer.data(lk + 374);
    const auto *lk_375 = buffer.data(lk + 375);
    const auto *lk_376 = buffer.data(lk + 376);
    const auto *lk_377 = buffer.data(lk + 377);
    const auto *lk_378 = buffer.data(lk + 378);
    const auto *lk_379 = buffer.data(lk + 379);
    const auto *lk_380 = buffer.data(lk + 380);
    const auto *lk_381 = buffer.data(lk + 381);
    const auto *lk_382 = buffer.data(lk + 382);
    const auto *lk_383 = buffer.data(lk + 383);
    const auto *lk_384 = buffer.data(lk + 384);
    const auto *lk_385 = buffer.data(lk + 385);
    const auto *lk_386 = buffer.data(lk + 386);
    const auto *lk_387 = buffer.data(lk + 387);
    const auto *lk_388 = buffer.data(lk + 388);
    const auto *lk_389 = buffer.data(lk + 389);
    const auto *lk_390 = buffer.data(lk + 390);
    const auto *lk_391 = buffer.data(lk + 391);
    const auto *lk_392 = buffer.data(lk + 392);
    const auto *lk_393 = buffer.data(lk + 393);
    const auto *lk_394 = buffer.data(lk + 394);
    const auto *lk_395 = buffer.data(lk + 395);
    const auto *lk_396 = buffer.data(lk + 396);
    const auto *lk_397 = buffer.data(lk + 397);
    const auto *lk_398 = buffer.data(lk + 398);
    const auto *lk_399 = buffer.data(lk + 399);
    const auto *lk_400 = buffer.data(lk + 400);
    const auto *lk_401 = buffer.data(lk + 401);
    const auto *lk_402 = buffer.data(lk + 402);
    const auto *lk_403 = buffer.data(lk + 403);
    const auto *lk_404 = buffer.data(lk + 404);
    const auto *lk_405 = buffer.data(lk + 405);
    const auto *lk_406 = buffer.data(lk + 406);
    const auto *lk_407 = buffer.data(lk + 407);
    const auto *lk_408 = buffer.data(lk + 408);
    const auto *lk_409 = buffer.data(lk + 409);
    const auto *lk_410 = buffer.data(lk + 410);
    const auto *lk_411 = buffer.data(lk + 411);
    const auto *lk_412 = buffer.data(lk + 412);
    const auto *lk_413 = buffer.data(lk + 413);
    const auto *lk_414 = buffer.data(lk + 414);
    const auto *lk_415 = buffer.data(lk + 415);
    const auto *lk_416 = buffer.data(lk + 416);
    const auto *lk_417 = buffer.data(lk + 417);
    const auto *lk_418 = buffer.data(lk + 418);
    const auto *lk_419 = buffer.data(lk + 419);
    const auto *lk_420 = buffer.data(lk + 420);
    const auto *lk_421 = buffer.data(lk + 421);
    const auto *lk_422 = buffer.data(lk + 422);
    const auto *lk_423 = buffer.data(lk + 423);
    const auto *lk_424 = buffer.data(lk + 424);
    const auto *lk_425 = buffer.data(lk + 425);
    const auto *lk_426 = buffer.data(lk + 426);
    const auto *lk_427 = buffer.data(lk + 427);
    const auto *lk_428 = buffer.data(lk + 428);
    const auto *lk_429 = buffer.data(lk + 429);
    const auto *lk_430 = buffer.data(lk + 430);
    const auto *lk_431 = buffer.data(lk + 431);
    const auto *lk_432 = buffer.data(lk + 432);
    const auto *lk_433 = buffer.data(lk + 433);
    const auto *lk_434 = buffer.data(lk + 434);
    const auto *lk_435 = buffer.data(lk + 435);
    const auto *lk_436 = buffer.data(lk + 436);
    const auto *lk_437 = buffer.data(lk + 437);
    const auto *lk_438 = buffer.data(lk + 438);
    const auto *lk_439 = buffer.data(lk + 439);
    const auto *lk_440 = buffer.data(lk + 440);
    const auto *lk_441 = buffer.data(lk + 441);
    const auto *lk_442 = buffer.data(lk + 442);
    const auto *lk_443 = buffer.data(lk + 443);
    const auto *lk_444 = buffer.data(lk + 444);
    const auto *lk_445 = buffer.data(lk + 445);
    const auto *lk_446 = buffer.data(lk + 446);
    const auto *lk_447 = buffer.data(lk + 447);
    const auto *lk_448 = buffer.data(lk + 448);
    const auto *lk_449 = buffer.data(lk + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ik_300, ik_301, ik_302, ik_303, \
                         ik_304, lk_300, lk_301, lk_302, lk_303, \
                         lk_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * ik_300[k]
                   + f_0 * lk_300[k];

        t_301[k] = -4.0 * ik_301[k]
                   + f_0 * lk_301[k];

        t_302[k] = -4.0 * ik_302[k]
                   + f_0 * lk_302[k];

        t_303[k] = -4.0 * ik_303[k]
                   + f_0 * lk_303[k];

        t_304[k] = -4.0 * ik_304[k]
                   + f_0 * lk_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ik_305, ik_306, ik_307, ik_308, \
                         ik_309, lk_305, lk_306, lk_307, lk_308, \
                         lk_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * ik_305[k]
                   + f_0 * lk_305[k];

        t_306[k] = -4.0 * ik_306[k]
                   + f_0 * lk_306[k];

        t_307[k] = -4.0 * ik_307[k]
                   + f_0 * lk_307[k];

        t_308[k] = -4.0 * ik_308[k]
                   + f_0 * lk_308[k];

        t_309[k] = -4.0 * ik_309[k]
                   + f_0 * lk_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ik_310, ik_311, ik_312, ik_313, \
                         ik_314, lk_310, lk_311, lk_312, lk_313, \
                         lk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -4.0 * ik_310[k]
                   + f_0 * lk_310[k];

        t_311[k] = -4.0 * ik_311[k]
                   + f_0 * lk_311[k];

        t_312[k] = -4.0 * ik_312[k]
                   + f_0 * lk_312[k];

        t_313[k] = -4.0 * ik_313[k]
                   + f_0 * lk_313[k];

        t_314[k] = -4.0 * ik_314[k]
                   + f_0 * lk_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ik_315, ik_316, ik_317, ik_318, \
                         ik_319, lk_315, lk_316, lk_317, lk_318, \
                         lk_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -4.0 * ik_315[k]
                   + f_0 * lk_315[k];

        t_316[k] = -4.0 * ik_316[k]
                   + f_0 * lk_316[k];

        t_317[k] = -4.0 * ik_317[k]
                   + f_0 * lk_317[k];

        t_318[k] = -4.0 * ik_318[k]
                   + f_0 * lk_318[k];

        t_319[k] = -4.0 * ik_319[k]
                   + f_0 * lk_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ik_320, ik_321, ik_322, ik_323, \
                         ik_324, lk_320, lk_321, lk_322, lk_323, \
                         lk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -4.0 * ik_320[k]
                   + f_0 * lk_320[k];

        t_321[k] = -4.0 * ik_321[k]
                   + f_0 * lk_321[k];

        t_322[k] = -4.0 * ik_322[k]
                   + f_0 * lk_322[k];

        t_323[k] = -4.0 * ik_323[k]
                   + f_0 * lk_323[k];

        t_324[k] = -4.0 * ik_324[k]
                   + f_0 * lk_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ik_325, ik_326, ik_327, ik_328, \
                         ik_329, lk_325, lk_326, lk_327, lk_328, \
                         lk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -4.0 * ik_325[k]
                   + f_0 * lk_325[k];

        t_326[k] = -4.0 * ik_326[k]
                   + f_0 * lk_326[k];

        t_327[k] = -4.0 * ik_327[k]
                   + f_0 * lk_327[k];

        t_328[k] = -4.0 * ik_328[k]
                   + f_0 * lk_328[k];

        t_329[k] = -4.0 * ik_329[k]
                   + f_0 * lk_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ik_330, ik_331, ik_332, ik_333, \
                         ik_334, lk_330, lk_331, lk_332, lk_333, \
                         lk_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -4.0 * ik_330[k]
                   + f_0 * lk_330[k];

        t_331[k] = -4.0 * ik_331[k]
                   + f_0 * lk_331[k];

        t_332[k] = -4.0 * ik_332[k]
                   + f_0 * lk_332[k];

        t_333[k] = -4.0 * ik_333[k]
                   + f_0 * lk_333[k];

        t_334[k] = -4.0 * ik_334[k]
                   + f_0 * lk_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ik_335, ik_336, ik_337, ik_338, \
                         ik_339, lk_335, lk_336, lk_337, lk_338, \
                         lk_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -4.0 * ik_335[k]
                   + f_0 * lk_335[k];

        t_336[k] = -4.0 * ik_336[k]
                   + f_0 * lk_336[k];

        t_337[k] = -4.0 * ik_337[k]
                   + f_0 * lk_337[k];

        t_338[k] = -4.0 * ik_338[k]
                   + f_0 * lk_338[k];

        t_339[k] = -4.0 * ik_339[k]
                   + f_0 * lk_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ik_340, ik_341, ik_342, ik_343, \
                         ik_344, lk_340, lk_341, lk_342, lk_343, \
                         lk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -4.0 * ik_340[k]
                   + f_0 * lk_340[k];

        t_341[k] = -4.0 * ik_341[k]
                   + f_0 * lk_341[k];

        t_342[k] = -4.0 * ik_342[k]
                   + f_0 * lk_342[k];

        t_343[k] = -4.0 * ik_343[k]
                   + f_0 * lk_343[k];

        t_344[k] = -4.0 * ik_344[k]
                   + f_0 * lk_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ik_345, ik_346, ik_347, ik_348, \
                         ik_349, lk_345, lk_346, lk_347, lk_348, \
                         lk_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -4.0 * ik_345[k]
                   + f_0 * lk_345[k];

        t_346[k] = -4.0 * ik_346[k]
                   + f_0 * lk_346[k];

        t_347[k] = -4.0 * ik_347[k]
                   + f_0 * lk_347[k];

        t_348[k] = -4.0 * ik_348[k]
                   + f_0 * lk_348[k];

        t_349[k] = -4.0 * ik_349[k]
                   + f_0 * lk_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ik_350, ik_351, ik_352, ik_353, \
                         ik_354, lk_350, lk_351, lk_352, lk_353, \
                         lk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -4.0 * ik_350[k]
                   + f_0 * lk_350[k];

        t_351[k] = -4.0 * ik_351[k]
                   + f_0 * lk_351[k];

        t_352[k] = -4.0 * ik_352[k]
                   + f_0 * lk_352[k];

        t_353[k] = -4.0 * ik_353[k]
                   + f_0 * lk_353[k];

        t_354[k] = -4.0 * ik_354[k]
                   + f_0 * lk_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ik_355, ik_356, ik_357, ik_358, \
                         ik_359, lk_355, lk_356, lk_357, lk_358, \
                         lk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -4.0 * ik_355[k]
                   + f_0 * lk_355[k];

        t_356[k] = -4.0 * ik_356[k]
                   + f_0 * lk_356[k];

        t_357[k] = -4.0 * ik_357[k]
                   + f_0 * lk_357[k];

        t_358[k] = -4.0 * ik_358[k]
                   + f_0 * lk_358[k];

        t_359[k] = -4.0 * ik_359[k]
                   + f_0 * lk_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ik_360, ik_361, ik_362, ik_363, \
                         ik_364, lk_360, lk_361, lk_362, lk_363, \
                         lk_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -3.0 * ik_360[k]
                   + f_0 * lk_360[k];

        t_361[k] = -3.0 * ik_361[k]
                   + f_0 * lk_361[k];

        t_362[k] = -3.0 * ik_362[k]
                   + f_0 * lk_362[k];

        t_363[k] = -3.0 * ik_363[k]
                   + f_0 * lk_363[k];

        t_364[k] = -3.0 * ik_364[k]
                   + f_0 * lk_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ik_365, ik_366, ik_367, ik_368, \
                         ik_369, lk_365, lk_366, lk_367, lk_368, \
                         lk_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -3.0 * ik_365[k]
                   + f_0 * lk_365[k];

        t_366[k] = -3.0 * ik_366[k]
                   + f_0 * lk_366[k];

        t_367[k] = -3.0 * ik_367[k]
                   + f_0 * lk_367[k];

        t_368[k] = -3.0 * ik_368[k]
                   + f_0 * lk_368[k];

        t_369[k] = -3.0 * ik_369[k]
                   + f_0 * lk_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ik_370, ik_371, ik_372, ik_373, \
                         ik_374, lk_370, lk_371, lk_372, lk_373, \
                         lk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -3.0 * ik_370[k]
                   + f_0 * lk_370[k];

        t_371[k] = -3.0 * ik_371[k]
                   + f_0 * lk_371[k];

        t_372[k] = -3.0 * ik_372[k]
                   + f_0 * lk_372[k];

        t_373[k] = -3.0 * ik_373[k]
                   + f_0 * lk_373[k];

        t_374[k] = -3.0 * ik_374[k]
                   + f_0 * lk_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ik_375, ik_376, ik_377, ik_378, \
                         ik_379, lk_375, lk_376, lk_377, lk_378, \
                         lk_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -3.0 * ik_375[k]
                   + f_0 * lk_375[k];

        t_376[k] = -3.0 * ik_376[k]
                   + f_0 * lk_376[k];

        t_377[k] = -3.0 * ik_377[k]
                   + f_0 * lk_377[k];

        t_378[k] = -3.0 * ik_378[k]
                   + f_0 * lk_378[k];

        t_379[k] = -3.0 * ik_379[k]
                   + f_0 * lk_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ik_380, ik_381, ik_382, ik_383, \
                         ik_384, lk_380, lk_381, lk_382, lk_383, \
                         lk_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -3.0 * ik_380[k]
                   + f_0 * lk_380[k];

        t_381[k] = -3.0 * ik_381[k]
                   + f_0 * lk_381[k];

        t_382[k] = -3.0 * ik_382[k]
                   + f_0 * lk_382[k];

        t_383[k] = -3.0 * ik_383[k]
                   + f_0 * lk_383[k];

        t_384[k] = -3.0 * ik_384[k]
                   + f_0 * lk_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ik_385, ik_386, ik_387, ik_388, \
                         ik_389, lk_385, lk_386, lk_387, lk_388, \
                         lk_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -3.0 * ik_385[k]
                   + f_0 * lk_385[k];

        t_386[k] = -3.0 * ik_386[k]
                   + f_0 * lk_386[k];

        t_387[k] = -3.0 * ik_387[k]
                   + f_0 * lk_387[k];

        t_388[k] = -3.0 * ik_388[k]
                   + f_0 * lk_388[k];

        t_389[k] = -3.0 * ik_389[k]
                   + f_0 * lk_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ik_390, ik_391, ik_392, ik_393, \
                         ik_394, lk_390, lk_391, lk_392, lk_393, \
                         lk_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -3.0 * ik_390[k]
                   + f_0 * lk_390[k];

        t_391[k] = -3.0 * ik_391[k]
                   + f_0 * lk_391[k];

        t_392[k] = -3.0 * ik_392[k]
                   + f_0 * lk_392[k];

        t_393[k] = -3.0 * ik_393[k]
                   + f_0 * lk_393[k];

        t_394[k] = -3.0 * ik_394[k]
                   + f_0 * lk_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ik_395, ik_396, ik_397, ik_398, \
                         ik_399, lk_395, lk_396, lk_397, lk_398, \
                         lk_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -3.0 * ik_395[k]
                   + f_0 * lk_395[k];

        t_396[k] = -3.0 * ik_396[k]
                   + f_0 * lk_396[k];

        t_397[k] = -3.0 * ik_397[k]
                   + f_0 * lk_397[k];

        t_398[k] = -3.0 * ik_398[k]
                   + f_0 * lk_398[k];

        t_399[k] = -3.0 * ik_399[k]
                   + f_0 * lk_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ik_400, ik_401, ik_402, ik_403, \
                         ik_404, lk_400, lk_401, lk_402, lk_403, \
                         lk_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -3.0 * ik_400[k]
                   + f_0 * lk_400[k];

        t_401[k] = -3.0 * ik_401[k]
                   + f_0 * lk_401[k];

        t_402[k] = -3.0 * ik_402[k]
                   + f_0 * lk_402[k];

        t_403[k] = -3.0 * ik_403[k]
                   + f_0 * lk_403[k];

        t_404[k] = -3.0 * ik_404[k]
                   + f_0 * lk_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ik_405, ik_406, ik_407, ik_408, \
                         ik_409, lk_405, lk_406, lk_407, lk_408, \
                         lk_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -3.0 * ik_405[k]
                   + f_0 * lk_405[k];

        t_406[k] = -3.0 * ik_406[k]
                   + f_0 * lk_406[k];

        t_407[k] = -3.0 * ik_407[k]
                   + f_0 * lk_407[k];

        t_408[k] = -3.0 * ik_408[k]
                   + f_0 * lk_408[k];

        t_409[k] = -3.0 * ik_409[k]
                   + f_0 * lk_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ik_410, ik_411, ik_412, ik_413, \
                         ik_414, lk_410, lk_411, lk_412, lk_413, \
                         lk_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -3.0 * ik_410[k]
                   + f_0 * lk_410[k];

        t_411[k] = -3.0 * ik_411[k]
                   + f_0 * lk_411[k];

        t_412[k] = -3.0 * ik_412[k]
                   + f_0 * lk_412[k];

        t_413[k] = -3.0 * ik_413[k]
                   + f_0 * lk_413[k];

        t_414[k] = -3.0 * ik_414[k]
                   + f_0 * lk_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ik_415, ik_416, ik_417, ik_418, \
                         ik_419, lk_415, lk_416, lk_417, lk_418, \
                         lk_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -3.0 * ik_415[k]
                   + f_0 * lk_415[k];

        t_416[k] = -3.0 * ik_416[k]
                   + f_0 * lk_416[k];

        t_417[k] = -3.0 * ik_417[k]
                   + f_0 * lk_417[k];

        t_418[k] = -3.0 * ik_418[k]
                   + f_0 * lk_418[k];

        t_419[k] = -3.0 * ik_419[k]
                   + f_0 * lk_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ik_420, ik_421, ik_422, ik_423, \
                         ik_424, lk_420, lk_421, lk_422, lk_423, \
                         lk_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -3.0 * ik_420[k]
                   + f_0 * lk_420[k];

        t_421[k] = -3.0 * ik_421[k]
                   + f_0 * lk_421[k];

        t_422[k] = -3.0 * ik_422[k]
                   + f_0 * lk_422[k];

        t_423[k] = -3.0 * ik_423[k]
                   + f_0 * lk_423[k];

        t_424[k] = -3.0 * ik_424[k]
                   + f_0 * lk_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ik_425, ik_426, ik_427, ik_428, \
                         ik_429, lk_425, lk_426, lk_427, lk_428, \
                         lk_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -3.0 * ik_425[k]
                   + f_0 * lk_425[k];

        t_426[k] = -3.0 * ik_426[k]
                   + f_0 * lk_426[k];

        t_427[k] = -3.0 * ik_427[k]
                   + f_0 * lk_427[k];

        t_428[k] = -3.0 * ik_428[k]
                   + f_0 * lk_428[k];

        t_429[k] = -3.0 * ik_429[k]
                   + f_0 * lk_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ik_430, ik_431, ik_432, ik_433, \
                         ik_434, lk_430, lk_431, lk_432, lk_433, \
                         lk_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -3.0 * ik_430[k]
                   + f_0 * lk_430[k];

        t_431[k] = -3.0 * ik_431[k]
                   + f_0 * lk_431[k];

        t_432[k] = -3.0 * ik_432[k]
                   + f_0 * lk_432[k];

        t_433[k] = -3.0 * ik_433[k]
                   + f_0 * lk_433[k];

        t_434[k] = -3.0 * ik_434[k]
                   + f_0 * lk_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ik_435, ik_436, ik_437, ik_438, \
                         ik_439, lk_435, lk_436, lk_437, lk_438, \
                         lk_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -3.0 * ik_435[k]
                   + f_0 * lk_435[k];

        t_436[k] = -3.0 * ik_436[k]
                   + f_0 * lk_436[k];

        t_437[k] = -3.0 * ik_437[k]
                   + f_0 * lk_437[k];

        t_438[k] = -3.0 * ik_438[k]
                   + f_0 * lk_438[k];

        t_439[k] = -3.0 * ik_439[k]
                   + f_0 * lk_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ik_440, ik_441, ik_442, ik_443, \
                         ik_444, lk_440, lk_441, lk_442, lk_443, \
                         lk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -3.0 * ik_440[k]
                   + f_0 * lk_440[k];

        t_441[k] = -3.0 * ik_441[k]
                   + f_0 * lk_441[k];

        t_442[k] = -3.0 * ik_442[k]
                   + f_0 * lk_442[k];

        t_443[k] = -3.0 * ik_443[k]
                   + f_0 * lk_443[k];

        t_444[k] = -3.0 * ik_444[k]
                   + f_0 * lk_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ik_445, ik_446, ik_447, ik_448, \
                         ik_449, lk_445, lk_446, lk_447, lk_448, \
                         lk_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -3.0 * ik_445[k]
                   + f_0 * lk_445[k];

        t_446[k] = -3.0 * ik_446[k]
                   + f_0 * lk_446[k];

        t_447[k] = -3.0 * ik_447[k]
                   + f_0 * lk_447[k];

        t_448[k] = -3.0 * ik_448[k]
                   + f_0 * lk_448[k];

        t_449[k] = -3.0 * ik_449[k]
                   + f_0 * lk_449[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_450 = buffer.data(ik + 450);
    const auto *ik_451 = buffer.data(ik + 451);
    const auto *ik_452 = buffer.data(ik + 452);
    const auto *ik_453 = buffer.data(ik + 453);
    const auto *ik_454 = buffer.data(ik + 454);
    const auto *ik_455 = buffer.data(ik + 455);
    const auto *ik_456 = buffer.data(ik + 456);
    const auto *ik_457 = buffer.data(ik + 457);
    const auto *ik_458 = buffer.data(ik + 458);
    const auto *ik_459 = buffer.data(ik + 459);
    const auto *ik_460 = buffer.data(ik + 460);
    const auto *ik_461 = buffer.data(ik + 461);
    const auto *ik_462 = buffer.data(ik + 462);
    const auto *ik_463 = buffer.data(ik + 463);
    const auto *ik_464 = buffer.data(ik + 464);
    const auto *ik_465 = buffer.data(ik + 465);
    const auto *ik_466 = buffer.data(ik + 466);
    const auto *ik_467 = buffer.data(ik + 467);
    const auto *ik_468 = buffer.data(ik + 468);
    const auto *ik_469 = buffer.data(ik + 469);
    const auto *ik_470 = buffer.data(ik + 470);
    const auto *ik_471 = buffer.data(ik + 471);
    const auto *ik_472 = buffer.data(ik + 472);
    const auto *ik_473 = buffer.data(ik + 473);
    const auto *ik_474 = buffer.data(ik + 474);
    const auto *ik_475 = buffer.data(ik + 475);
    const auto *ik_476 = buffer.data(ik + 476);
    const auto *ik_477 = buffer.data(ik + 477);
    const auto *ik_478 = buffer.data(ik + 478);
    const auto *ik_479 = buffer.data(ik + 479);
    const auto *ik_480 = buffer.data(ik + 480);
    const auto *ik_481 = buffer.data(ik + 481);
    const auto *ik_482 = buffer.data(ik + 482);
    const auto *ik_483 = buffer.data(ik + 483);
    const auto *ik_484 = buffer.data(ik + 484);
    const auto *ik_485 = buffer.data(ik + 485);
    const auto *ik_486 = buffer.data(ik + 486);
    const auto *ik_487 = buffer.data(ik + 487);
    const auto *ik_488 = buffer.data(ik + 488);
    const auto *ik_489 = buffer.data(ik + 489);
    const auto *ik_490 = buffer.data(ik + 490);
    const auto *ik_491 = buffer.data(ik + 491);
    const auto *ik_492 = buffer.data(ik + 492);
    const auto *ik_493 = buffer.data(ik + 493);
    const auto *ik_494 = buffer.data(ik + 494);
    const auto *ik_495 = buffer.data(ik + 495);
    const auto *ik_496 = buffer.data(ik + 496);
    const auto *ik_497 = buffer.data(ik + 497);
    const auto *ik_498 = buffer.data(ik + 498);
    const auto *ik_499 = buffer.data(ik + 499);
    const auto *ik_500 = buffer.data(ik + 500);
    const auto *ik_501 = buffer.data(ik + 501);
    const auto *ik_502 = buffer.data(ik + 502);
    const auto *ik_503 = buffer.data(ik + 503);
    const auto *ik_504 = buffer.data(ik + 504);
    const auto *ik_505 = buffer.data(ik + 505);
    const auto *ik_506 = buffer.data(ik + 506);
    const auto *ik_507 = buffer.data(ik + 507);
    const auto *ik_508 = buffer.data(ik + 508);
    const auto *ik_509 = buffer.data(ik + 509);
    const auto *ik_510 = buffer.data(ik + 510);
    const auto *ik_511 = buffer.data(ik + 511);
    const auto *ik_512 = buffer.data(ik + 512);
    const auto *ik_513 = buffer.data(ik + 513);
    const auto *ik_514 = buffer.data(ik + 514);
    const auto *ik_515 = buffer.data(ik + 515);
    const auto *ik_516 = buffer.data(ik + 516);
    const auto *ik_517 = buffer.data(ik + 517);
    const auto *ik_518 = buffer.data(ik + 518);
    const auto *ik_519 = buffer.data(ik + 519);
    const auto *ik_520 = buffer.data(ik + 520);
    const auto *ik_521 = buffer.data(ik + 521);
    const auto *ik_522 = buffer.data(ik + 522);
    const auto *ik_523 = buffer.data(ik + 523);
    const auto *ik_524 = buffer.data(ik + 524);
    const auto *ik_525 = buffer.data(ik + 525);
    const auto *ik_526 = buffer.data(ik + 526);
    const auto *ik_527 = buffer.data(ik + 527);
    const auto *ik_528 = buffer.data(ik + 528);
    const auto *ik_529 = buffer.data(ik + 529);
    const auto *ik_530 = buffer.data(ik + 530);
    const auto *ik_531 = buffer.data(ik + 531);
    const auto *ik_532 = buffer.data(ik + 532);
    const auto *ik_533 = buffer.data(ik + 533);
    const auto *ik_534 = buffer.data(ik + 534);
    const auto *ik_535 = buffer.data(ik + 535);
    const auto *ik_536 = buffer.data(ik + 536);
    const auto *ik_537 = buffer.data(ik + 537);
    const auto *ik_538 = buffer.data(ik + 538);
    const auto *ik_539 = buffer.data(ik + 539);
    const auto *ik_540 = buffer.data(ik + 540);
    const auto *ik_541 = buffer.data(ik + 541);
    const auto *ik_542 = buffer.data(ik + 542);
    const auto *ik_543 = buffer.data(ik + 543);
    const auto *ik_544 = buffer.data(ik + 544);
    const auto *ik_545 = buffer.data(ik + 545);
    const auto *ik_546 = buffer.data(ik + 546);
    const auto *ik_547 = buffer.data(ik + 547);
    const auto *ik_548 = buffer.data(ik + 548);
    const auto *ik_549 = buffer.data(ik + 549);
    const auto *ik_550 = buffer.data(ik + 550);
    const auto *ik_551 = buffer.data(ik + 551);
    const auto *ik_552 = buffer.data(ik + 552);
    const auto *ik_553 = buffer.data(ik + 553);
    const auto *ik_554 = buffer.data(ik + 554);
    const auto *ik_555 = buffer.data(ik + 555);
    const auto *ik_556 = buffer.data(ik + 556);
    const auto *ik_557 = buffer.data(ik + 557);
    const auto *ik_558 = buffer.data(ik + 558);
    const auto *ik_559 = buffer.data(ik + 559);
    const auto *ik_560 = buffer.data(ik + 560);
    const auto *ik_561 = buffer.data(ik + 561);
    const auto *ik_562 = buffer.data(ik + 562);
    const auto *ik_563 = buffer.data(ik + 563);
    const auto *ik_564 = buffer.data(ik + 564);
    const auto *ik_565 = buffer.data(ik + 565);
    const auto *ik_566 = buffer.data(ik + 566);
    const auto *ik_567 = buffer.data(ik + 567);
    const auto *ik_568 = buffer.data(ik + 568);
    const auto *ik_569 = buffer.data(ik + 569);
    const auto *ik_570 = buffer.data(ik + 570);
    const auto *ik_571 = buffer.data(ik + 571);
    const auto *ik_572 = buffer.data(ik + 572);
    const auto *ik_573 = buffer.data(ik + 573);
    const auto *ik_574 = buffer.data(ik + 574);
    const auto *ik_575 = buffer.data(ik + 575);
    const auto *ik_576 = buffer.data(ik + 576);
    const auto *ik_577 = buffer.data(ik + 577);
    const auto *ik_578 = buffer.data(ik + 578);
    const auto *ik_579 = buffer.data(ik + 579);
    const auto *ik_580 = buffer.data(ik + 580);
    const auto *ik_581 = buffer.data(ik + 581);
    const auto *ik_582 = buffer.data(ik + 582);
    const auto *ik_583 = buffer.data(ik + 583);
    const auto *ik_584 = buffer.data(ik + 584);
    const auto *ik_585 = buffer.data(ik + 585);
    const auto *ik_586 = buffer.data(ik + 586);
    const auto *ik_587 = buffer.data(ik + 587);
    const auto *ik_588 = buffer.data(ik + 588);
    const auto *ik_589 = buffer.data(ik + 589);
    const auto *ik_590 = buffer.data(ik + 590);
    const auto *ik_591 = buffer.data(ik + 591);
    const auto *ik_592 = buffer.data(ik + 592);
    const auto *ik_593 = buffer.data(ik + 593);
    const auto *ik_594 = buffer.data(ik + 594);
    const auto *ik_595 = buffer.data(ik + 595);
    const auto *ik_596 = buffer.data(ik + 596);
    const auto *ik_597 = buffer.data(ik + 597);
    const auto *ik_598 = buffer.data(ik + 598);
    const auto *ik_599 = buffer.data(ik + 599);

    const auto *lk_450 = buffer.data(lk + 450);
    const auto *lk_451 = buffer.data(lk + 451);
    const auto *lk_452 = buffer.data(lk + 452);
    const auto *lk_453 = buffer.data(lk + 453);
    const auto *lk_454 = buffer.data(lk + 454);
    const auto *lk_455 = buffer.data(lk + 455);
    const auto *lk_456 = buffer.data(lk + 456);
    const auto *lk_457 = buffer.data(lk + 457);
    const auto *lk_458 = buffer.data(lk + 458);
    const auto *lk_459 = buffer.data(lk + 459);
    const auto *lk_460 = buffer.data(lk + 460);
    const auto *lk_461 = buffer.data(lk + 461);
    const auto *lk_462 = buffer.data(lk + 462);
    const auto *lk_463 = buffer.data(lk + 463);
    const auto *lk_464 = buffer.data(lk + 464);
    const auto *lk_465 = buffer.data(lk + 465);
    const auto *lk_466 = buffer.data(lk + 466);
    const auto *lk_467 = buffer.data(lk + 467);
    const auto *lk_468 = buffer.data(lk + 468);
    const auto *lk_469 = buffer.data(lk + 469);
    const auto *lk_470 = buffer.data(lk + 470);
    const auto *lk_471 = buffer.data(lk + 471);
    const auto *lk_472 = buffer.data(lk + 472);
    const auto *lk_473 = buffer.data(lk + 473);
    const auto *lk_474 = buffer.data(lk + 474);
    const auto *lk_475 = buffer.data(lk + 475);
    const auto *lk_476 = buffer.data(lk + 476);
    const auto *lk_477 = buffer.data(lk + 477);
    const auto *lk_478 = buffer.data(lk + 478);
    const auto *lk_479 = buffer.data(lk + 479);
    const auto *lk_480 = buffer.data(lk + 480);
    const auto *lk_481 = buffer.data(lk + 481);
    const auto *lk_482 = buffer.data(lk + 482);
    const auto *lk_483 = buffer.data(lk + 483);
    const auto *lk_484 = buffer.data(lk + 484);
    const auto *lk_485 = buffer.data(lk + 485);
    const auto *lk_486 = buffer.data(lk + 486);
    const auto *lk_487 = buffer.data(lk + 487);
    const auto *lk_488 = buffer.data(lk + 488);
    const auto *lk_489 = buffer.data(lk + 489);
    const auto *lk_490 = buffer.data(lk + 490);
    const auto *lk_491 = buffer.data(lk + 491);
    const auto *lk_492 = buffer.data(lk + 492);
    const auto *lk_493 = buffer.data(lk + 493);
    const auto *lk_494 = buffer.data(lk + 494);
    const auto *lk_495 = buffer.data(lk + 495);
    const auto *lk_496 = buffer.data(lk + 496);
    const auto *lk_497 = buffer.data(lk + 497);
    const auto *lk_498 = buffer.data(lk + 498);
    const auto *lk_499 = buffer.data(lk + 499);
    const auto *lk_500 = buffer.data(lk + 500);
    const auto *lk_501 = buffer.data(lk + 501);
    const auto *lk_502 = buffer.data(lk + 502);
    const auto *lk_503 = buffer.data(lk + 503);
    const auto *lk_504 = buffer.data(lk + 504);
    const auto *lk_505 = buffer.data(lk + 505);
    const auto *lk_506 = buffer.data(lk + 506);
    const auto *lk_507 = buffer.data(lk + 507);
    const auto *lk_508 = buffer.data(lk + 508);
    const auto *lk_509 = buffer.data(lk + 509);
    const auto *lk_510 = buffer.data(lk + 510);
    const auto *lk_511 = buffer.data(lk + 511);
    const auto *lk_512 = buffer.data(lk + 512);
    const auto *lk_513 = buffer.data(lk + 513);
    const auto *lk_514 = buffer.data(lk + 514);
    const auto *lk_515 = buffer.data(lk + 515);
    const auto *lk_516 = buffer.data(lk + 516);
    const auto *lk_517 = buffer.data(lk + 517);
    const auto *lk_518 = buffer.data(lk + 518);
    const auto *lk_519 = buffer.data(lk + 519);
    const auto *lk_520 = buffer.data(lk + 520);
    const auto *lk_521 = buffer.data(lk + 521);
    const auto *lk_522 = buffer.data(lk + 522);
    const auto *lk_523 = buffer.data(lk + 523);
    const auto *lk_524 = buffer.data(lk + 524);
    const auto *lk_525 = buffer.data(lk + 525);
    const auto *lk_526 = buffer.data(lk + 526);
    const auto *lk_527 = buffer.data(lk + 527);
    const auto *lk_528 = buffer.data(lk + 528);
    const auto *lk_529 = buffer.data(lk + 529);
    const auto *lk_530 = buffer.data(lk + 530);
    const auto *lk_531 = buffer.data(lk + 531);
    const auto *lk_532 = buffer.data(lk + 532);
    const auto *lk_533 = buffer.data(lk + 533);
    const auto *lk_534 = buffer.data(lk + 534);
    const auto *lk_535 = buffer.data(lk + 535);
    const auto *lk_536 = buffer.data(lk + 536);
    const auto *lk_537 = buffer.data(lk + 537);
    const auto *lk_538 = buffer.data(lk + 538);
    const auto *lk_539 = buffer.data(lk + 539);
    const auto *lk_540 = buffer.data(lk + 540);
    const auto *lk_541 = buffer.data(lk + 541);
    const auto *lk_542 = buffer.data(lk + 542);
    const auto *lk_543 = buffer.data(lk + 543);
    const auto *lk_544 = buffer.data(lk + 544);
    const auto *lk_545 = buffer.data(lk + 545);
    const auto *lk_546 = buffer.data(lk + 546);
    const auto *lk_547 = buffer.data(lk + 547);
    const auto *lk_548 = buffer.data(lk + 548);
    const auto *lk_549 = buffer.data(lk + 549);
    const auto *lk_550 = buffer.data(lk + 550);
    const auto *lk_551 = buffer.data(lk + 551);
    const auto *lk_552 = buffer.data(lk + 552);
    const auto *lk_553 = buffer.data(lk + 553);
    const auto *lk_554 = buffer.data(lk + 554);
    const auto *lk_555 = buffer.data(lk + 555);
    const auto *lk_556 = buffer.data(lk + 556);
    const auto *lk_557 = buffer.data(lk + 557);
    const auto *lk_558 = buffer.data(lk + 558);
    const auto *lk_559 = buffer.data(lk + 559);
    const auto *lk_560 = buffer.data(lk + 560);
    const auto *lk_561 = buffer.data(lk + 561);
    const auto *lk_562 = buffer.data(lk + 562);
    const auto *lk_563 = buffer.data(lk + 563);
    const auto *lk_564 = buffer.data(lk + 564);
    const auto *lk_565 = buffer.data(lk + 565);
    const auto *lk_566 = buffer.data(lk + 566);
    const auto *lk_567 = buffer.data(lk + 567);
    const auto *lk_568 = buffer.data(lk + 568);
    const auto *lk_569 = buffer.data(lk + 569);
    const auto *lk_570 = buffer.data(lk + 570);
    const auto *lk_571 = buffer.data(lk + 571);
    const auto *lk_572 = buffer.data(lk + 572);
    const auto *lk_573 = buffer.data(lk + 573);
    const auto *lk_574 = buffer.data(lk + 574);
    const auto *lk_575 = buffer.data(lk + 575);
    const auto *lk_576 = buffer.data(lk + 576);
    const auto *lk_577 = buffer.data(lk + 577);
    const auto *lk_578 = buffer.data(lk + 578);
    const auto *lk_579 = buffer.data(lk + 579);
    const auto *lk_580 = buffer.data(lk + 580);
    const auto *lk_581 = buffer.data(lk + 581);
    const auto *lk_582 = buffer.data(lk + 582);
    const auto *lk_583 = buffer.data(lk + 583);
    const auto *lk_584 = buffer.data(lk + 584);
    const auto *lk_585 = buffer.data(lk + 585);
    const auto *lk_586 = buffer.data(lk + 586);
    const auto *lk_587 = buffer.data(lk + 587);
    const auto *lk_588 = buffer.data(lk + 588);
    const auto *lk_589 = buffer.data(lk + 589);
    const auto *lk_590 = buffer.data(lk + 590);
    const auto *lk_591 = buffer.data(lk + 591);
    const auto *lk_592 = buffer.data(lk + 592);
    const auto *lk_593 = buffer.data(lk + 593);
    const auto *lk_594 = buffer.data(lk + 594);
    const auto *lk_595 = buffer.data(lk + 595);
    const auto *lk_596 = buffer.data(lk + 596);
    const auto *lk_597 = buffer.data(lk + 597);
    const auto *lk_598 = buffer.data(lk + 598);
    const auto *lk_599 = buffer.data(lk + 599);

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ik_450, ik_451, ik_452, ik_453, \
                         ik_454, lk_450, lk_451, lk_452, lk_453, \
                         lk_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -3.0 * ik_450[k]
                   + f_0 * lk_450[k];

        t_451[k] = -3.0 * ik_451[k]
                   + f_0 * lk_451[k];

        t_452[k] = -3.0 * ik_452[k]
                   + f_0 * lk_452[k];

        t_453[k] = -3.0 * ik_453[k]
                   + f_0 * lk_453[k];

        t_454[k] = -3.0 * ik_454[k]
                   + f_0 * lk_454[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ik_455, ik_456, ik_457, ik_458, \
                         ik_459, lk_455, lk_456, lk_457, lk_458, \
                         lk_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -3.0 * ik_455[k]
                   + f_0 * lk_455[k];

        t_456[k] = -3.0 * ik_456[k]
                   + f_0 * lk_456[k];

        t_457[k] = -3.0 * ik_457[k]
                   + f_0 * lk_457[k];

        t_458[k] = -3.0 * ik_458[k]
                   + f_0 * lk_458[k];

        t_459[k] = -3.0 * ik_459[k]
                   + f_0 * lk_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ik_460, ik_461, ik_462, ik_463, \
                         ik_464, lk_460, lk_461, lk_462, lk_463, \
                         lk_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -3.0 * ik_460[k]
                   + f_0 * lk_460[k];

        t_461[k] = -3.0 * ik_461[k]
                   + f_0 * lk_461[k];

        t_462[k] = -3.0 * ik_462[k]
                   + f_0 * lk_462[k];

        t_463[k] = -3.0 * ik_463[k]
                   + f_0 * lk_463[k];

        t_464[k] = -3.0 * ik_464[k]
                   + f_0 * lk_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ik_465, ik_466, ik_467, ik_468, \
                         ik_469, lk_465, lk_466, lk_467, lk_468, \
                         lk_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -3.0 * ik_465[k]
                   + f_0 * lk_465[k];

        t_466[k] = -3.0 * ik_466[k]
                   + f_0 * lk_466[k];

        t_467[k] = -3.0 * ik_467[k]
                   + f_0 * lk_467[k];

        t_468[k] = -3.0 * ik_468[k]
                   + f_0 * lk_468[k];

        t_469[k] = -3.0 * ik_469[k]
                   + f_0 * lk_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ik_470, ik_471, ik_472, ik_473, \
                         ik_474, lk_470, lk_471, lk_472, lk_473, \
                         lk_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -3.0 * ik_470[k]
                   + f_0 * lk_470[k];

        t_471[k] = -3.0 * ik_471[k]
                   + f_0 * lk_471[k];

        t_472[k] = -3.0 * ik_472[k]
                   + f_0 * lk_472[k];

        t_473[k] = -3.0 * ik_473[k]
                   + f_0 * lk_473[k];

        t_474[k] = -3.0 * ik_474[k]
                   + f_0 * lk_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ik_475, ik_476, ik_477, ik_478, \
                         ik_479, lk_475, lk_476, lk_477, lk_478, \
                         lk_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -3.0 * ik_475[k]
                   + f_0 * lk_475[k];

        t_476[k] = -3.0 * ik_476[k]
                   + f_0 * lk_476[k];

        t_477[k] = -3.0 * ik_477[k]
                   + f_0 * lk_477[k];

        t_478[k] = -3.0 * ik_478[k]
                   + f_0 * lk_478[k];

        t_479[k] = -3.0 * ik_479[k]
                   + f_0 * lk_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ik_480, ik_481, ik_482, ik_483, \
                         ik_484, lk_480, lk_481, lk_482, lk_483, \
                         lk_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -3.0 * ik_480[k]
                   + f_0 * lk_480[k];

        t_481[k] = -3.0 * ik_481[k]
                   + f_0 * lk_481[k];

        t_482[k] = -3.0 * ik_482[k]
                   + f_0 * lk_482[k];

        t_483[k] = -3.0 * ik_483[k]
                   + f_0 * lk_483[k];

        t_484[k] = -3.0 * ik_484[k]
                   + f_0 * lk_484[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ik_485, ik_486, ik_487, ik_488, \
                         ik_489, lk_485, lk_486, lk_487, lk_488, \
                         lk_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -3.0 * ik_485[k]
                   + f_0 * lk_485[k];

        t_486[k] = -3.0 * ik_486[k]
                   + f_0 * lk_486[k];

        t_487[k] = -3.0 * ik_487[k]
                   + f_0 * lk_487[k];

        t_488[k] = -3.0 * ik_488[k]
                   + f_0 * lk_488[k];

        t_489[k] = -3.0 * ik_489[k]
                   + f_0 * lk_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ik_490, ik_491, ik_492, ik_493, \
                         ik_494, lk_490, lk_491, lk_492, lk_493, \
                         lk_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -3.0 * ik_490[k]
                   + f_0 * lk_490[k];

        t_491[k] = -3.0 * ik_491[k]
                   + f_0 * lk_491[k];

        t_492[k] = -3.0 * ik_492[k]
                   + f_0 * lk_492[k];

        t_493[k] = -3.0 * ik_493[k]
                   + f_0 * lk_493[k];

        t_494[k] = -3.0 * ik_494[k]
                   + f_0 * lk_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ik_495, ik_496, ik_497, ik_498, \
                         ik_499, lk_495, lk_496, lk_497, lk_498, \
                         lk_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -3.0 * ik_495[k]
                   + f_0 * lk_495[k];

        t_496[k] = -3.0 * ik_496[k]
                   + f_0 * lk_496[k];

        t_497[k] = -3.0 * ik_497[k]
                   + f_0 * lk_497[k];

        t_498[k] = -3.0 * ik_498[k]
                   + f_0 * lk_498[k];

        t_499[k] = -3.0 * ik_499[k]
                   + f_0 * lk_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ik_500, ik_501, ik_502, ik_503, \
                         ik_504, lk_500, lk_501, lk_502, lk_503, \
                         lk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -3.0 * ik_500[k]
                   + f_0 * lk_500[k];

        t_501[k] = -3.0 * ik_501[k]
                   + f_0 * lk_501[k];

        t_502[k] = -3.0 * ik_502[k]
                   + f_0 * lk_502[k];

        t_503[k] = -3.0 * ik_503[k]
                   + f_0 * lk_503[k];

        t_504[k] = -3.0 * ik_504[k]
                   + f_0 * lk_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ik_505, ik_506, ik_507, ik_508, \
                         ik_509, lk_505, lk_506, lk_507, lk_508, \
                         lk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -3.0 * ik_505[k]
                   + f_0 * lk_505[k];

        t_506[k] = -3.0 * ik_506[k]
                   + f_0 * lk_506[k];

        t_507[k] = -3.0 * ik_507[k]
                   + f_0 * lk_507[k];

        t_508[k] = -3.0 * ik_508[k]
                   + f_0 * lk_508[k];

        t_509[k] = -3.0 * ik_509[k]
                   + f_0 * lk_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ik_510, ik_511, ik_512, ik_513, \
                         ik_514, lk_510, lk_511, lk_512, lk_513, \
                         lk_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -3.0 * ik_510[k]
                   + f_0 * lk_510[k];

        t_511[k] = -3.0 * ik_511[k]
                   + f_0 * lk_511[k];

        t_512[k] = -3.0 * ik_512[k]
                   + f_0 * lk_512[k];

        t_513[k] = -3.0 * ik_513[k]
                   + f_0 * lk_513[k];

        t_514[k] = -3.0 * ik_514[k]
                   + f_0 * lk_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ik_515, ik_516, ik_517, ik_518, \
                         ik_519, lk_515, lk_516, lk_517, lk_518, \
                         lk_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -3.0 * ik_515[k]
                   + f_0 * lk_515[k];

        t_516[k] = -3.0 * ik_516[k]
                   + f_0 * lk_516[k];

        t_517[k] = -3.0 * ik_517[k]
                   + f_0 * lk_517[k];

        t_518[k] = -3.0 * ik_518[k]
                   + f_0 * lk_518[k];

        t_519[k] = -3.0 * ik_519[k]
                   + f_0 * lk_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ik_520, ik_521, ik_522, ik_523, \
                         ik_524, lk_520, lk_521, lk_522, lk_523, \
                         lk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -3.0 * ik_520[k]
                   + f_0 * lk_520[k];

        t_521[k] = -3.0 * ik_521[k]
                   + f_0 * lk_521[k];

        t_522[k] = -3.0 * ik_522[k]
                   + f_0 * lk_522[k];

        t_523[k] = -3.0 * ik_523[k]
                   + f_0 * lk_523[k];

        t_524[k] = -3.0 * ik_524[k]
                   + f_0 * lk_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ik_525, ik_526, ik_527, ik_528, \
                         ik_529, lk_525, lk_526, lk_527, lk_528, \
                         lk_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -3.0 * ik_525[k]
                   + f_0 * lk_525[k];

        t_526[k] = -3.0 * ik_526[k]
                   + f_0 * lk_526[k];

        t_527[k] = -3.0 * ik_527[k]
                   + f_0 * lk_527[k];

        t_528[k] = -3.0 * ik_528[k]
                   + f_0 * lk_528[k];

        t_529[k] = -3.0 * ik_529[k]
                   + f_0 * lk_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ik_530, ik_531, ik_532, ik_533, \
                         ik_534, lk_530, lk_531, lk_532, lk_533, \
                         lk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -3.0 * ik_530[k]
                   + f_0 * lk_530[k];

        t_531[k] = -3.0 * ik_531[k]
                   + f_0 * lk_531[k];

        t_532[k] = -3.0 * ik_532[k]
                   + f_0 * lk_532[k];

        t_533[k] = -3.0 * ik_533[k]
                   + f_0 * lk_533[k];

        t_534[k] = -3.0 * ik_534[k]
                   + f_0 * lk_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ik_535, ik_536, ik_537, ik_538, \
                         ik_539, lk_535, lk_536, lk_537, lk_538, \
                         lk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -3.0 * ik_535[k]
                   + f_0 * lk_535[k];

        t_536[k] = -3.0 * ik_536[k]
                   + f_0 * lk_536[k];

        t_537[k] = -3.0 * ik_537[k]
                   + f_0 * lk_537[k];

        t_538[k] = -3.0 * ik_538[k]
                   + f_0 * lk_538[k];

        t_539[k] = -3.0 * ik_539[k]
                   + f_0 * lk_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ik_540, ik_541, ik_542, ik_543, \
                         ik_544, lk_540, lk_541, lk_542, lk_543, \
                         lk_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -2.0 * ik_540[k]
                   + f_0 * lk_540[k];

        t_541[k] = -2.0 * ik_541[k]
                   + f_0 * lk_541[k];

        t_542[k] = -2.0 * ik_542[k]
                   + f_0 * lk_542[k];

        t_543[k] = -2.0 * ik_543[k]
                   + f_0 * lk_543[k];

        t_544[k] = -2.0 * ik_544[k]
                   + f_0 * lk_544[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ik_545, ik_546, ik_547, ik_548, \
                         ik_549, lk_545, lk_546, lk_547, lk_548, \
                         lk_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -2.0 * ik_545[k]
                   + f_0 * lk_545[k];

        t_546[k] = -2.0 * ik_546[k]
                   + f_0 * lk_546[k];

        t_547[k] = -2.0 * ik_547[k]
                   + f_0 * lk_547[k];

        t_548[k] = -2.0 * ik_548[k]
                   + f_0 * lk_548[k];

        t_549[k] = -2.0 * ik_549[k]
                   + f_0 * lk_549[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ik_550, ik_551, ik_552, ik_553, \
                         ik_554, lk_550, lk_551, lk_552, lk_553, \
                         lk_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -2.0 * ik_550[k]
                   + f_0 * lk_550[k];

        t_551[k] = -2.0 * ik_551[k]
                   + f_0 * lk_551[k];

        t_552[k] = -2.0 * ik_552[k]
                   + f_0 * lk_552[k];

        t_553[k] = -2.0 * ik_553[k]
                   + f_0 * lk_553[k];

        t_554[k] = -2.0 * ik_554[k]
                   + f_0 * lk_554[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ik_555, ik_556, ik_557, ik_558, \
                         ik_559, lk_555, lk_556, lk_557, lk_558, \
                         lk_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -2.0 * ik_555[k]
                   + f_0 * lk_555[k];

        t_556[k] = -2.0 * ik_556[k]
                   + f_0 * lk_556[k];

        t_557[k] = -2.0 * ik_557[k]
                   + f_0 * lk_557[k];

        t_558[k] = -2.0 * ik_558[k]
                   + f_0 * lk_558[k];

        t_559[k] = -2.0 * ik_559[k]
                   + f_0 * lk_559[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, ik_560, ik_561, ik_562, ik_563, \
                         ik_564, lk_560, lk_561, lk_562, lk_563, \
                         lk_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = -2.0 * ik_560[k]
                   + f_0 * lk_560[k];

        t_561[k] = -2.0 * ik_561[k]
                   + f_0 * lk_561[k];

        t_562[k] = -2.0 * ik_562[k]
                   + f_0 * lk_562[k];

        t_563[k] = -2.0 * ik_563[k]
                   + f_0 * lk_563[k];

        t_564[k] = -2.0 * ik_564[k]
                   + f_0 * lk_564[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, ik_565, ik_566, ik_567, ik_568, \
                         ik_569, lk_565, lk_566, lk_567, lk_568, \
                         lk_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = -2.0 * ik_565[k]
                   + f_0 * lk_565[k];

        t_566[k] = -2.0 * ik_566[k]
                   + f_0 * lk_566[k];

        t_567[k] = -2.0 * ik_567[k]
                   + f_0 * lk_567[k];

        t_568[k] = -2.0 * ik_568[k]
                   + f_0 * lk_568[k];

        t_569[k] = -2.0 * ik_569[k]
                   + f_0 * lk_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ik_570, ik_571, ik_572, ik_573, \
                         ik_574, lk_570, lk_571, lk_572, lk_573, \
                         lk_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -2.0 * ik_570[k]
                   + f_0 * lk_570[k];

        t_571[k] = -2.0 * ik_571[k]
                   + f_0 * lk_571[k];

        t_572[k] = -2.0 * ik_572[k]
                   + f_0 * lk_572[k];

        t_573[k] = -2.0 * ik_573[k]
                   + f_0 * lk_573[k];

        t_574[k] = -2.0 * ik_574[k]
                   + f_0 * lk_574[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ik_575, ik_576, ik_577, ik_578, \
                         ik_579, lk_575, lk_576, lk_577, lk_578, \
                         lk_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = -2.0 * ik_575[k]
                   + f_0 * lk_575[k];

        t_576[k] = -2.0 * ik_576[k]
                   + f_0 * lk_576[k];

        t_577[k] = -2.0 * ik_577[k]
                   + f_0 * lk_577[k];

        t_578[k] = -2.0 * ik_578[k]
                   + f_0 * lk_578[k];

        t_579[k] = -2.0 * ik_579[k]
                   + f_0 * lk_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, ik_580, ik_581, ik_582, ik_583, \
                         ik_584, lk_580, lk_581, lk_582, lk_583, \
                         lk_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -2.0 * ik_580[k]
                   + f_0 * lk_580[k];

        t_581[k] = -2.0 * ik_581[k]
                   + f_0 * lk_581[k];

        t_582[k] = -2.0 * ik_582[k]
                   + f_0 * lk_582[k];

        t_583[k] = -2.0 * ik_583[k]
                   + f_0 * lk_583[k];

        t_584[k] = -2.0 * ik_584[k]
                   + f_0 * lk_584[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, ik_585, ik_586, ik_587, ik_588, \
                         ik_589, lk_585, lk_586, lk_587, lk_588, \
                         lk_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -2.0 * ik_585[k]
                   + f_0 * lk_585[k];

        t_586[k] = -2.0 * ik_586[k]
                   + f_0 * lk_586[k];

        t_587[k] = -2.0 * ik_587[k]
                   + f_0 * lk_587[k];

        t_588[k] = -2.0 * ik_588[k]
                   + f_0 * lk_588[k];

        t_589[k] = -2.0 * ik_589[k]
                   + f_0 * lk_589[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, ik_590, ik_591, ik_592, ik_593, \
                         ik_594, lk_590, lk_591, lk_592, lk_593, \
                         lk_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -2.0 * ik_590[k]
                   + f_0 * lk_590[k];

        t_591[k] = -2.0 * ik_591[k]
                   + f_0 * lk_591[k];

        t_592[k] = -2.0 * ik_592[k]
                   + f_0 * lk_592[k];

        t_593[k] = -2.0 * ik_593[k]
                   + f_0 * lk_593[k];

        t_594[k] = -2.0 * ik_594[k]
                   + f_0 * lk_594[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, ik_595, ik_596, ik_597, ik_598, \
                         ik_599, lk_595, lk_596, lk_597, lk_598, \
                         lk_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -2.0 * ik_595[k]
                   + f_0 * lk_595[k];

        t_596[k] = -2.0 * ik_596[k]
                   + f_0 * lk_596[k];

        t_597[k] = -2.0 * ik_597[k]
                   + f_0 * lk_597[k];

        t_598[k] = -2.0 * ik_598[k]
                   + f_0 * lk_598[k];

        t_599[k] = -2.0 * ik_599[k]
                   + f_0 * lk_599[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_600 = buffer.data(ik + 600);
    const auto *ik_601 = buffer.data(ik + 601);
    const auto *ik_602 = buffer.data(ik + 602);
    const auto *ik_603 = buffer.data(ik + 603);
    const auto *ik_604 = buffer.data(ik + 604);
    const auto *ik_605 = buffer.data(ik + 605);
    const auto *ik_606 = buffer.data(ik + 606);
    const auto *ik_607 = buffer.data(ik + 607);
    const auto *ik_608 = buffer.data(ik + 608);
    const auto *ik_609 = buffer.data(ik + 609);
    const auto *ik_610 = buffer.data(ik + 610);
    const auto *ik_611 = buffer.data(ik + 611);
    const auto *ik_612 = buffer.data(ik + 612);
    const auto *ik_613 = buffer.data(ik + 613);
    const auto *ik_614 = buffer.data(ik + 614);
    const auto *ik_615 = buffer.data(ik + 615);
    const auto *ik_616 = buffer.data(ik + 616);
    const auto *ik_617 = buffer.data(ik + 617);
    const auto *ik_618 = buffer.data(ik + 618);
    const auto *ik_619 = buffer.data(ik + 619);
    const auto *ik_620 = buffer.data(ik + 620);
    const auto *ik_621 = buffer.data(ik + 621);
    const auto *ik_622 = buffer.data(ik + 622);
    const auto *ik_623 = buffer.data(ik + 623);
    const auto *ik_624 = buffer.data(ik + 624);
    const auto *ik_625 = buffer.data(ik + 625);
    const auto *ik_626 = buffer.data(ik + 626);
    const auto *ik_627 = buffer.data(ik + 627);
    const auto *ik_628 = buffer.data(ik + 628);
    const auto *ik_629 = buffer.data(ik + 629);
    const auto *ik_630 = buffer.data(ik + 630);
    const auto *ik_631 = buffer.data(ik + 631);
    const auto *ik_632 = buffer.data(ik + 632);
    const auto *ik_633 = buffer.data(ik + 633);
    const auto *ik_634 = buffer.data(ik + 634);
    const auto *ik_635 = buffer.data(ik + 635);
    const auto *ik_636 = buffer.data(ik + 636);
    const auto *ik_637 = buffer.data(ik + 637);
    const auto *ik_638 = buffer.data(ik + 638);
    const auto *ik_639 = buffer.data(ik + 639);
    const auto *ik_640 = buffer.data(ik + 640);
    const auto *ik_641 = buffer.data(ik + 641);
    const auto *ik_642 = buffer.data(ik + 642);
    const auto *ik_643 = buffer.data(ik + 643);
    const auto *ik_644 = buffer.data(ik + 644);
    const auto *ik_645 = buffer.data(ik + 645);
    const auto *ik_646 = buffer.data(ik + 646);
    const auto *ik_647 = buffer.data(ik + 647);
    const auto *ik_648 = buffer.data(ik + 648);
    const auto *ik_649 = buffer.data(ik + 649);
    const auto *ik_650 = buffer.data(ik + 650);
    const auto *ik_651 = buffer.data(ik + 651);
    const auto *ik_652 = buffer.data(ik + 652);
    const auto *ik_653 = buffer.data(ik + 653);
    const auto *ik_654 = buffer.data(ik + 654);
    const auto *ik_655 = buffer.data(ik + 655);
    const auto *ik_656 = buffer.data(ik + 656);
    const auto *ik_657 = buffer.data(ik + 657);
    const auto *ik_658 = buffer.data(ik + 658);
    const auto *ik_659 = buffer.data(ik + 659);
    const auto *ik_660 = buffer.data(ik + 660);
    const auto *ik_661 = buffer.data(ik + 661);
    const auto *ik_662 = buffer.data(ik + 662);
    const auto *ik_663 = buffer.data(ik + 663);
    const auto *ik_664 = buffer.data(ik + 664);
    const auto *ik_665 = buffer.data(ik + 665);
    const auto *ik_666 = buffer.data(ik + 666);
    const auto *ik_667 = buffer.data(ik + 667);
    const auto *ik_668 = buffer.data(ik + 668);
    const auto *ik_669 = buffer.data(ik + 669);
    const auto *ik_670 = buffer.data(ik + 670);
    const auto *ik_671 = buffer.data(ik + 671);
    const auto *ik_672 = buffer.data(ik + 672);
    const auto *ik_673 = buffer.data(ik + 673);
    const auto *ik_674 = buffer.data(ik + 674);
    const auto *ik_675 = buffer.data(ik + 675);
    const auto *ik_676 = buffer.data(ik + 676);
    const auto *ik_677 = buffer.data(ik + 677);
    const auto *ik_678 = buffer.data(ik + 678);
    const auto *ik_679 = buffer.data(ik + 679);
    const auto *ik_680 = buffer.data(ik + 680);
    const auto *ik_681 = buffer.data(ik + 681);
    const auto *ik_682 = buffer.data(ik + 682);
    const auto *ik_683 = buffer.data(ik + 683);
    const auto *ik_684 = buffer.data(ik + 684);
    const auto *ik_685 = buffer.data(ik + 685);
    const auto *ik_686 = buffer.data(ik + 686);
    const auto *ik_687 = buffer.data(ik + 687);
    const auto *ik_688 = buffer.data(ik + 688);
    const auto *ik_689 = buffer.data(ik + 689);
    const auto *ik_690 = buffer.data(ik + 690);
    const auto *ik_691 = buffer.data(ik + 691);
    const auto *ik_692 = buffer.data(ik + 692);
    const auto *ik_693 = buffer.data(ik + 693);
    const auto *ik_694 = buffer.data(ik + 694);
    const auto *ik_695 = buffer.data(ik + 695);
    const auto *ik_696 = buffer.data(ik + 696);
    const auto *ik_697 = buffer.data(ik + 697);
    const auto *ik_698 = buffer.data(ik + 698);
    const auto *ik_699 = buffer.data(ik + 699);
    const auto *ik_700 = buffer.data(ik + 700);
    const auto *ik_701 = buffer.data(ik + 701);
    const auto *ik_702 = buffer.data(ik + 702);
    const auto *ik_703 = buffer.data(ik + 703);
    const auto *ik_704 = buffer.data(ik + 704);
    const auto *ik_705 = buffer.data(ik + 705);
    const auto *ik_706 = buffer.data(ik + 706);
    const auto *ik_707 = buffer.data(ik + 707);
    const auto *ik_708 = buffer.data(ik + 708);
    const auto *ik_709 = buffer.data(ik + 709);
    const auto *ik_710 = buffer.data(ik + 710);
    const auto *ik_711 = buffer.data(ik + 711);
    const auto *ik_712 = buffer.data(ik + 712);
    const auto *ik_713 = buffer.data(ik + 713);
    const auto *ik_714 = buffer.data(ik + 714);
    const auto *ik_715 = buffer.data(ik + 715);
    const auto *ik_716 = buffer.data(ik + 716);
    const auto *ik_717 = buffer.data(ik + 717);
    const auto *ik_718 = buffer.data(ik + 718);
    const auto *ik_719 = buffer.data(ik + 719);
    const auto *ik_720 = buffer.data(ik + 720);
    const auto *ik_721 = buffer.data(ik + 721);
    const auto *ik_722 = buffer.data(ik + 722);
    const auto *ik_723 = buffer.data(ik + 723);
    const auto *ik_724 = buffer.data(ik + 724);
    const auto *ik_725 = buffer.data(ik + 725);
    const auto *ik_726 = buffer.data(ik + 726);
    const auto *ik_727 = buffer.data(ik + 727);
    const auto *ik_728 = buffer.data(ik + 728);
    const auto *ik_729 = buffer.data(ik + 729);
    const auto *ik_730 = buffer.data(ik + 730);
    const auto *ik_731 = buffer.data(ik + 731);
    const auto *ik_732 = buffer.data(ik + 732);
    const auto *ik_733 = buffer.data(ik + 733);
    const auto *ik_734 = buffer.data(ik + 734);
    const auto *ik_735 = buffer.data(ik + 735);
    const auto *ik_736 = buffer.data(ik + 736);
    const auto *ik_737 = buffer.data(ik + 737);
    const auto *ik_738 = buffer.data(ik + 738);
    const auto *ik_739 = buffer.data(ik + 739);
    const auto *ik_740 = buffer.data(ik + 740);
    const auto *ik_741 = buffer.data(ik + 741);
    const auto *ik_742 = buffer.data(ik + 742);
    const auto *ik_743 = buffer.data(ik + 743);
    const auto *ik_744 = buffer.data(ik + 744);
    const auto *ik_745 = buffer.data(ik + 745);
    const auto *ik_746 = buffer.data(ik + 746);
    const auto *ik_747 = buffer.data(ik + 747);
    const auto *ik_748 = buffer.data(ik + 748);
    const auto *ik_749 = buffer.data(ik + 749);

    const auto *lk_600 = buffer.data(lk + 600);
    const auto *lk_601 = buffer.data(lk + 601);
    const auto *lk_602 = buffer.data(lk + 602);
    const auto *lk_603 = buffer.data(lk + 603);
    const auto *lk_604 = buffer.data(lk + 604);
    const auto *lk_605 = buffer.data(lk + 605);
    const auto *lk_606 = buffer.data(lk + 606);
    const auto *lk_607 = buffer.data(lk + 607);
    const auto *lk_608 = buffer.data(lk + 608);
    const auto *lk_609 = buffer.data(lk + 609);
    const auto *lk_610 = buffer.data(lk + 610);
    const auto *lk_611 = buffer.data(lk + 611);
    const auto *lk_612 = buffer.data(lk + 612);
    const auto *lk_613 = buffer.data(lk + 613);
    const auto *lk_614 = buffer.data(lk + 614);
    const auto *lk_615 = buffer.data(lk + 615);
    const auto *lk_616 = buffer.data(lk + 616);
    const auto *lk_617 = buffer.data(lk + 617);
    const auto *lk_618 = buffer.data(lk + 618);
    const auto *lk_619 = buffer.data(lk + 619);
    const auto *lk_620 = buffer.data(lk + 620);
    const auto *lk_621 = buffer.data(lk + 621);
    const auto *lk_622 = buffer.data(lk + 622);
    const auto *lk_623 = buffer.data(lk + 623);
    const auto *lk_624 = buffer.data(lk + 624);
    const auto *lk_625 = buffer.data(lk + 625);
    const auto *lk_626 = buffer.data(lk + 626);
    const auto *lk_627 = buffer.data(lk + 627);
    const auto *lk_628 = buffer.data(lk + 628);
    const auto *lk_629 = buffer.data(lk + 629);
    const auto *lk_630 = buffer.data(lk + 630);
    const auto *lk_631 = buffer.data(lk + 631);
    const auto *lk_632 = buffer.data(lk + 632);
    const auto *lk_633 = buffer.data(lk + 633);
    const auto *lk_634 = buffer.data(lk + 634);
    const auto *lk_635 = buffer.data(lk + 635);
    const auto *lk_636 = buffer.data(lk + 636);
    const auto *lk_637 = buffer.data(lk + 637);
    const auto *lk_638 = buffer.data(lk + 638);
    const auto *lk_639 = buffer.data(lk + 639);
    const auto *lk_640 = buffer.data(lk + 640);
    const auto *lk_641 = buffer.data(lk + 641);
    const auto *lk_642 = buffer.data(lk + 642);
    const auto *lk_643 = buffer.data(lk + 643);
    const auto *lk_644 = buffer.data(lk + 644);
    const auto *lk_645 = buffer.data(lk + 645);
    const auto *lk_646 = buffer.data(lk + 646);
    const auto *lk_647 = buffer.data(lk + 647);
    const auto *lk_648 = buffer.data(lk + 648);
    const auto *lk_649 = buffer.data(lk + 649);
    const auto *lk_650 = buffer.data(lk + 650);
    const auto *lk_651 = buffer.data(lk + 651);
    const auto *lk_652 = buffer.data(lk + 652);
    const auto *lk_653 = buffer.data(lk + 653);
    const auto *lk_654 = buffer.data(lk + 654);
    const auto *lk_655 = buffer.data(lk + 655);
    const auto *lk_656 = buffer.data(lk + 656);
    const auto *lk_657 = buffer.data(lk + 657);
    const auto *lk_658 = buffer.data(lk + 658);
    const auto *lk_659 = buffer.data(lk + 659);
    const auto *lk_660 = buffer.data(lk + 660);
    const auto *lk_661 = buffer.data(lk + 661);
    const auto *lk_662 = buffer.data(lk + 662);
    const auto *lk_663 = buffer.data(lk + 663);
    const auto *lk_664 = buffer.data(lk + 664);
    const auto *lk_665 = buffer.data(lk + 665);
    const auto *lk_666 = buffer.data(lk + 666);
    const auto *lk_667 = buffer.data(lk + 667);
    const auto *lk_668 = buffer.data(lk + 668);
    const auto *lk_669 = buffer.data(lk + 669);
    const auto *lk_670 = buffer.data(lk + 670);
    const auto *lk_671 = buffer.data(lk + 671);
    const auto *lk_672 = buffer.data(lk + 672);
    const auto *lk_673 = buffer.data(lk + 673);
    const auto *lk_674 = buffer.data(lk + 674);
    const auto *lk_675 = buffer.data(lk + 675);
    const auto *lk_676 = buffer.data(lk + 676);
    const auto *lk_677 = buffer.data(lk + 677);
    const auto *lk_678 = buffer.data(lk + 678);
    const auto *lk_679 = buffer.data(lk + 679);
    const auto *lk_680 = buffer.data(lk + 680);
    const auto *lk_681 = buffer.data(lk + 681);
    const auto *lk_682 = buffer.data(lk + 682);
    const auto *lk_683 = buffer.data(lk + 683);
    const auto *lk_684 = buffer.data(lk + 684);
    const auto *lk_685 = buffer.data(lk + 685);
    const auto *lk_686 = buffer.data(lk + 686);
    const auto *lk_687 = buffer.data(lk + 687);
    const auto *lk_688 = buffer.data(lk + 688);
    const auto *lk_689 = buffer.data(lk + 689);
    const auto *lk_690 = buffer.data(lk + 690);
    const auto *lk_691 = buffer.data(lk + 691);
    const auto *lk_692 = buffer.data(lk + 692);
    const auto *lk_693 = buffer.data(lk + 693);
    const auto *lk_694 = buffer.data(lk + 694);
    const auto *lk_695 = buffer.data(lk + 695);
    const auto *lk_696 = buffer.data(lk + 696);
    const auto *lk_697 = buffer.data(lk + 697);
    const auto *lk_698 = buffer.data(lk + 698);
    const auto *lk_699 = buffer.data(lk + 699);
    const auto *lk_700 = buffer.data(lk + 700);
    const auto *lk_701 = buffer.data(lk + 701);
    const auto *lk_702 = buffer.data(lk + 702);
    const auto *lk_703 = buffer.data(lk + 703);
    const auto *lk_704 = buffer.data(lk + 704);
    const auto *lk_705 = buffer.data(lk + 705);
    const auto *lk_706 = buffer.data(lk + 706);
    const auto *lk_707 = buffer.data(lk + 707);
    const auto *lk_708 = buffer.data(lk + 708);
    const auto *lk_709 = buffer.data(lk + 709);
    const auto *lk_710 = buffer.data(lk + 710);
    const auto *lk_711 = buffer.data(lk + 711);
    const auto *lk_712 = buffer.data(lk + 712);
    const auto *lk_713 = buffer.data(lk + 713);
    const auto *lk_714 = buffer.data(lk + 714);
    const auto *lk_715 = buffer.data(lk + 715);
    const auto *lk_716 = buffer.data(lk + 716);
    const auto *lk_717 = buffer.data(lk + 717);
    const auto *lk_718 = buffer.data(lk + 718);
    const auto *lk_719 = buffer.data(lk + 719);
    const auto *lk_720 = buffer.data(lk + 720);
    const auto *lk_721 = buffer.data(lk + 721);
    const auto *lk_722 = buffer.data(lk + 722);
    const auto *lk_723 = buffer.data(lk + 723);
    const auto *lk_724 = buffer.data(lk + 724);
    const auto *lk_725 = buffer.data(lk + 725);
    const auto *lk_726 = buffer.data(lk + 726);
    const auto *lk_727 = buffer.data(lk + 727);
    const auto *lk_728 = buffer.data(lk + 728);
    const auto *lk_729 = buffer.data(lk + 729);
    const auto *lk_730 = buffer.data(lk + 730);
    const auto *lk_731 = buffer.data(lk + 731);
    const auto *lk_732 = buffer.data(lk + 732);
    const auto *lk_733 = buffer.data(lk + 733);
    const auto *lk_734 = buffer.data(lk + 734);
    const auto *lk_735 = buffer.data(lk + 735);
    const auto *lk_736 = buffer.data(lk + 736);
    const auto *lk_737 = buffer.data(lk + 737);
    const auto *lk_738 = buffer.data(lk + 738);
    const auto *lk_739 = buffer.data(lk + 739);
    const auto *lk_740 = buffer.data(lk + 740);
    const auto *lk_741 = buffer.data(lk + 741);
    const auto *lk_742 = buffer.data(lk + 742);
    const auto *lk_743 = buffer.data(lk + 743);
    const auto *lk_744 = buffer.data(lk + 744);
    const auto *lk_745 = buffer.data(lk + 745);
    const auto *lk_746 = buffer.data(lk + 746);
    const auto *lk_747 = buffer.data(lk + 747);
    const auto *lk_748 = buffer.data(lk + 748);
    const auto *lk_749 = buffer.data(lk + 749);

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, ik_600, ik_601, ik_602, ik_603, \
                         ik_604, lk_600, lk_601, lk_602, lk_603, \
                         lk_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -2.0 * ik_600[k]
                   + f_0 * lk_600[k];

        t_601[k] = -2.0 * ik_601[k]
                   + f_0 * lk_601[k];

        t_602[k] = -2.0 * ik_602[k]
                   + f_0 * lk_602[k];

        t_603[k] = -2.0 * ik_603[k]
                   + f_0 * lk_603[k];

        t_604[k] = -2.0 * ik_604[k]
                   + f_0 * lk_604[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, ik_605, ik_606, ik_607, ik_608, \
                         ik_609, lk_605, lk_606, lk_607, lk_608, \
                         lk_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -2.0 * ik_605[k]
                   + f_0 * lk_605[k];

        t_606[k] = -2.0 * ik_606[k]
                   + f_0 * lk_606[k];

        t_607[k] = -2.0 * ik_607[k]
                   + f_0 * lk_607[k];

        t_608[k] = -2.0 * ik_608[k]
                   + f_0 * lk_608[k];

        t_609[k] = -2.0 * ik_609[k]
                   + f_0 * lk_609[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, ik_610, ik_611, ik_612, ik_613, \
                         ik_614, lk_610, lk_611, lk_612, lk_613, \
                         lk_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -2.0 * ik_610[k]
                   + f_0 * lk_610[k];

        t_611[k] = -2.0 * ik_611[k]
                   + f_0 * lk_611[k];

        t_612[k] = -2.0 * ik_612[k]
                   + f_0 * lk_612[k];

        t_613[k] = -2.0 * ik_613[k]
                   + f_0 * lk_613[k];

        t_614[k] = -2.0 * ik_614[k]
                   + f_0 * lk_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, ik_615, ik_616, ik_617, ik_618, \
                         ik_619, lk_615, lk_616, lk_617, lk_618, \
                         lk_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -2.0 * ik_615[k]
                   + f_0 * lk_615[k];

        t_616[k] = -2.0 * ik_616[k]
                   + f_0 * lk_616[k];

        t_617[k] = -2.0 * ik_617[k]
                   + f_0 * lk_617[k];

        t_618[k] = -2.0 * ik_618[k]
                   + f_0 * lk_618[k];

        t_619[k] = -2.0 * ik_619[k]
                   + f_0 * lk_619[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, ik_620, ik_621, ik_622, ik_623, \
                         ik_624, lk_620, lk_621, lk_622, lk_623, \
                         lk_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -2.0 * ik_620[k]
                   + f_0 * lk_620[k];

        t_621[k] = -2.0 * ik_621[k]
                   + f_0 * lk_621[k];

        t_622[k] = -2.0 * ik_622[k]
                   + f_0 * lk_622[k];

        t_623[k] = -2.0 * ik_623[k]
                   + f_0 * lk_623[k];

        t_624[k] = -2.0 * ik_624[k]
                   + f_0 * lk_624[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, ik_625, ik_626, ik_627, ik_628, \
                         ik_629, lk_625, lk_626, lk_627, lk_628, \
                         lk_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -2.0 * ik_625[k]
                   + f_0 * lk_625[k];

        t_626[k] = -2.0 * ik_626[k]
                   + f_0 * lk_626[k];

        t_627[k] = -2.0 * ik_627[k]
                   + f_0 * lk_627[k];

        t_628[k] = -2.0 * ik_628[k]
                   + f_0 * lk_628[k];

        t_629[k] = -2.0 * ik_629[k]
                   + f_0 * lk_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, ik_630, ik_631, ik_632, ik_633, \
                         ik_634, lk_630, lk_631, lk_632, lk_633, \
                         lk_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -2.0 * ik_630[k]
                   + f_0 * lk_630[k];

        t_631[k] = -2.0 * ik_631[k]
                   + f_0 * lk_631[k];

        t_632[k] = -2.0 * ik_632[k]
                   + f_0 * lk_632[k];

        t_633[k] = -2.0 * ik_633[k]
                   + f_0 * lk_633[k];

        t_634[k] = -2.0 * ik_634[k]
                   + f_0 * lk_634[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, ik_635, ik_636, ik_637, ik_638, \
                         ik_639, lk_635, lk_636, lk_637, lk_638, \
                         lk_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -2.0 * ik_635[k]
                   + f_0 * lk_635[k];

        t_636[k] = -2.0 * ik_636[k]
                   + f_0 * lk_636[k];

        t_637[k] = -2.0 * ik_637[k]
                   + f_0 * lk_637[k];

        t_638[k] = -2.0 * ik_638[k]
                   + f_0 * lk_638[k];

        t_639[k] = -2.0 * ik_639[k]
                   + f_0 * lk_639[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, ik_640, ik_641, ik_642, ik_643, \
                         ik_644, lk_640, lk_641, lk_642, lk_643, \
                         lk_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -2.0 * ik_640[k]
                   + f_0 * lk_640[k];

        t_641[k] = -2.0 * ik_641[k]
                   + f_0 * lk_641[k];

        t_642[k] = -2.0 * ik_642[k]
                   + f_0 * lk_642[k];

        t_643[k] = -2.0 * ik_643[k]
                   + f_0 * lk_643[k];

        t_644[k] = -2.0 * ik_644[k]
                   + f_0 * lk_644[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, ik_645, ik_646, ik_647, ik_648, \
                         ik_649, lk_645, lk_646, lk_647, lk_648, \
                         lk_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -2.0 * ik_645[k]
                   + f_0 * lk_645[k];

        t_646[k] = -2.0 * ik_646[k]
                   + f_0 * lk_646[k];

        t_647[k] = -2.0 * ik_647[k]
                   + f_0 * lk_647[k];

        t_648[k] = -2.0 * ik_648[k]
                   + f_0 * lk_648[k];

        t_649[k] = -2.0 * ik_649[k]
                   + f_0 * lk_649[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, ik_650, ik_651, ik_652, ik_653, \
                         ik_654, lk_650, lk_651, lk_652, lk_653, \
                         lk_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -2.0 * ik_650[k]
                   + f_0 * lk_650[k];

        t_651[k] = -2.0 * ik_651[k]
                   + f_0 * lk_651[k];

        t_652[k] = -2.0 * ik_652[k]
                   + f_0 * lk_652[k];

        t_653[k] = -2.0 * ik_653[k]
                   + f_0 * lk_653[k];

        t_654[k] = -2.0 * ik_654[k]
                   + f_0 * lk_654[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, ik_655, ik_656, ik_657, ik_658, \
                         ik_659, lk_655, lk_656, lk_657, lk_658, \
                         lk_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -2.0 * ik_655[k]
                   + f_0 * lk_655[k];

        t_656[k] = -2.0 * ik_656[k]
                   + f_0 * lk_656[k];

        t_657[k] = -2.0 * ik_657[k]
                   + f_0 * lk_657[k];

        t_658[k] = -2.0 * ik_658[k]
                   + f_0 * lk_658[k];

        t_659[k] = -2.0 * ik_659[k]
                   + f_0 * lk_659[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, ik_660, ik_661, ik_662, ik_663, \
                         ik_664, lk_660, lk_661, lk_662, lk_663, \
                         lk_664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -2.0 * ik_660[k]
                   + f_0 * lk_660[k];

        t_661[k] = -2.0 * ik_661[k]
                   + f_0 * lk_661[k];

        t_662[k] = -2.0 * ik_662[k]
                   + f_0 * lk_662[k];

        t_663[k] = -2.0 * ik_663[k]
                   + f_0 * lk_663[k];

        t_664[k] = -2.0 * ik_664[k]
                   + f_0 * lk_664[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, ik_665, ik_666, ik_667, ik_668, \
                         ik_669, lk_665, lk_666, lk_667, lk_668, \
                         lk_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -2.0 * ik_665[k]
                   + f_0 * lk_665[k];

        t_666[k] = -2.0 * ik_666[k]
                   + f_0 * lk_666[k];

        t_667[k] = -2.0 * ik_667[k]
                   + f_0 * lk_667[k];

        t_668[k] = -2.0 * ik_668[k]
                   + f_0 * lk_668[k];

        t_669[k] = -2.0 * ik_669[k]
                   + f_0 * lk_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, ik_670, ik_671, ik_672, ik_673, \
                         ik_674, lk_670, lk_671, lk_672, lk_673, \
                         lk_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -2.0 * ik_670[k]
                   + f_0 * lk_670[k];

        t_671[k] = -2.0 * ik_671[k]
                   + f_0 * lk_671[k];

        t_672[k] = -2.0 * ik_672[k]
                   + f_0 * lk_672[k];

        t_673[k] = -2.0 * ik_673[k]
                   + f_0 * lk_673[k];

        t_674[k] = -2.0 * ik_674[k]
                   + f_0 * lk_674[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, ik_675, ik_676, ik_677, ik_678, \
                         ik_679, lk_675, lk_676, lk_677, lk_678, \
                         lk_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -2.0 * ik_675[k]
                   + f_0 * lk_675[k];

        t_676[k] = -2.0 * ik_676[k]
                   + f_0 * lk_676[k];

        t_677[k] = -2.0 * ik_677[k]
                   + f_0 * lk_677[k];

        t_678[k] = -2.0 * ik_678[k]
                   + f_0 * lk_678[k];

        t_679[k] = -2.0 * ik_679[k]
                   + f_0 * lk_679[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, ik_680, ik_681, ik_682, ik_683, \
                         ik_684, lk_680, lk_681, lk_682, lk_683, \
                         lk_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -2.0 * ik_680[k]
                   + f_0 * lk_680[k];

        t_681[k] = -2.0 * ik_681[k]
                   + f_0 * lk_681[k];

        t_682[k] = -2.0 * ik_682[k]
                   + f_0 * lk_682[k];

        t_683[k] = -2.0 * ik_683[k]
                   + f_0 * lk_683[k];

        t_684[k] = -2.0 * ik_684[k]
                   + f_0 * lk_684[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, ik_685, ik_686, ik_687, ik_688, \
                         ik_689, lk_685, lk_686, lk_687, lk_688, \
                         lk_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -2.0 * ik_685[k]
                   + f_0 * lk_685[k];

        t_686[k] = -2.0 * ik_686[k]
                   + f_0 * lk_686[k];

        t_687[k] = -2.0 * ik_687[k]
                   + f_0 * lk_687[k];

        t_688[k] = -2.0 * ik_688[k]
                   + f_0 * lk_688[k];

        t_689[k] = -2.0 * ik_689[k]
                   + f_0 * lk_689[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, ik_690, ik_691, ik_692, ik_693, \
                         ik_694, lk_690, lk_691, lk_692, lk_693, \
                         lk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -2.0 * ik_690[k]
                   + f_0 * lk_690[k];

        t_691[k] = -2.0 * ik_691[k]
                   + f_0 * lk_691[k];

        t_692[k] = -2.0 * ik_692[k]
                   + f_0 * lk_692[k];

        t_693[k] = -2.0 * ik_693[k]
                   + f_0 * lk_693[k];

        t_694[k] = -2.0 * ik_694[k]
                   + f_0 * lk_694[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, ik_695, ik_696, ik_697, ik_698, \
                         ik_699, lk_695, lk_696, lk_697, lk_698, \
                         lk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -2.0 * ik_695[k]
                   + f_0 * lk_695[k];

        t_696[k] = -2.0 * ik_696[k]
                   + f_0 * lk_696[k];

        t_697[k] = -2.0 * ik_697[k]
                   + f_0 * lk_697[k];

        t_698[k] = -2.0 * ik_698[k]
                   + f_0 * lk_698[k];

        t_699[k] = -2.0 * ik_699[k]
                   + f_0 * lk_699[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, ik_700, ik_701, ik_702, ik_703, \
                         ik_704, lk_700, lk_701, lk_702, lk_703, \
                         lk_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -2.0 * ik_700[k]
                   + f_0 * lk_700[k];

        t_701[k] = -2.0 * ik_701[k]
                   + f_0 * lk_701[k];

        t_702[k] = -2.0 * ik_702[k]
                   + f_0 * lk_702[k];

        t_703[k] = -2.0 * ik_703[k]
                   + f_0 * lk_703[k];

        t_704[k] = -2.0 * ik_704[k]
                   + f_0 * lk_704[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, ik_705, ik_706, ik_707, ik_708, \
                         ik_709, lk_705, lk_706, lk_707, lk_708, \
                         lk_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -2.0 * ik_705[k]
                   + f_0 * lk_705[k];

        t_706[k] = -2.0 * ik_706[k]
                   + f_0 * lk_706[k];

        t_707[k] = -2.0 * ik_707[k]
                   + f_0 * lk_707[k];

        t_708[k] = -2.0 * ik_708[k]
                   + f_0 * lk_708[k];

        t_709[k] = -2.0 * ik_709[k]
                   + f_0 * lk_709[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, ik_710, ik_711, ik_712, ik_713, \
                         ik_714, lk_710, lk_711, lk_712, lk_713, \
                         lk_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -2.0 * ik_710[k]
                   + f_0 * lk_710[k];

        t_711[k] = -2.0 * ik_711[k]
                   + f_0 * lk_711[k];

        t_712[k] = -2.0 * ik_712[k]
                   + f_0 * lk_712[k];

        t_713[k] = -2.0 * ik_713[k]
                   + f_0 * lk_713[k];

        t_714[k] = -2.0 * ik_714[k]
                   + f_0 * lk_714[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, ik_715, ik_716, ik_717, ik_718, \
                         ik_719, lk_715, lk_716, lk_717, lk_718, \
                         lk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -2.0 * ik_715[k]
                   + f_0 * lk_715[k];

        t_716[k] = -2.0 * ik_716[k]
                   + f_0 * lk_716[k];

        t_717[k] = -2.0 * ik_717[k]
                   + f_0 * lk_717[k];

        t_718[k] = -2.0 * ik_718[k]
                   + f_0 * lk_718[k];

        t_719[k] = -2.0 * ik_719[k]
                   + f_0 * lk_719[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, ik_720, ik_721, ik_722, ik_723, \
                         ik_724, lk_720, lk_721, lk_722, lk_723, \
                         lk_724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -2.0 * ik_720[k]
                   + f_0 * lk_720[k];

        t_721[k] = -2.0 * ik_721[k]
                   + f_0 * lk_721[k];

        t_722[k] = -2.0 * ik_722[k]
                   + f_0 * lk_722[k];

        t_723[k] = -2.0 * ik_723[k]
                   + f_0 * lk_723[k];

        t_724[k] = -2.0 * ik_724[k]
                   + f_0 * lk_724[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, ik_725, ik_726, ik_727, ik_728, \
                         ik_729, lk_725, lk_726, lk_727, lk_728, \
                         lk_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -2.0 * ik_725[k]
                   + f_0 * lk_725[k];

        t_726[k] = -2.0 * ik_726[k]
                   + f_0 * lk_726[k];

        t_727[k] = -2.0 * ik_727[k]
                   + f_0 * lk_727[k];

        t_728[k] = -2.0 * ik_728[k]
                   + f_0 * lk_728[k];

        t_729[k] = -2.0 * ik_729[k]
                   + f_0 * lk_729[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, ik_730, ik_731, ik_732, ik_733, \
                         ik_734, lk_730, lk_731, lk_732, lk_733, \
                         lk_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -2.0 * ik_730[k]
                   + f_0 * lk_730[k];

        t_731[k] = -2.0 * ik_731[k]
                   + f_0 * lk_731[k];

        t_732[k] = -2.0 * ik_732[k]
                   + f_0 * lk_732[k];

        t_733[k] = -2.0 * ik_733[k]
                   + f_0 * lk_733[k];

        t_734[k] = -2.0 * ik_734[k]
                   + f_0 * lk_734[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, ik_735, ik_736, ik_737, ik_738, \
                         ik_739, lk_735, lk_736, lk_737, lk_738, \
                         lk_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -2.0 * ik_735[k]
                   + f_0 * lk_735[k];

        t_736[k] = -2.0 * ik_736[k]
                   + f_0 * lk_736[k];

        t_737[k] = -2.0 * ik_737[k]
                   + f_0 * lk_737[k];

        t_738[k] = -2.0 * ik_738[k]
                   + f_0 * lk_738[k];

        t_739[k] = -2.0 * ik_739[k]
                   + f_0 * lk_739[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, ik_740, ik_741, ik_742, ik_743, \
                         ik_744, lk_740, lk_741, lk_742, lk_743, \
                         lk_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -2.0 * ik_740[k]
                   + f_0 * lk_740[k];

        t_741[k] = -2.0 * ik_741[k]
                   + f_0 * lk_741[k];

        t_742[k] = -2.0 * ik_742[k]
                   + f_0 * lk_742[k];

        t_743[k] = -2.0 * ik_743[k]
                   + f_0 * lk_743[k];

        t_744[k] = -2.0 * ik_744[k]
                   + f_0 * lk_744[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, ik_745, ik_746, ik_747, ik_748, \
                         ik_749, lk_745, lk_746, lk_747, lk_748, \
                         lk_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -2.0 * ik_745[k]
                   + f_0 * lk_745[k];

        t_746[k] = -2.0 * ik_746[k]
                   + f_0 * lk_746[k];

        t_747[k] = -2.0 * ik_747[k]
                   + f_0 * lk_747[k];

        t_748[k] = -2.0 * ik_748[k]
                   + f_0 * lk_748[k];

        t_749[k] = -2.0 * ik_749[k]
                   + f_0 * lk_749[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_750 = buffer.data(ik + 750);
    const auto *ik_751 = buffer.data(ik + 751);
    const auto *ik_752 = buffer.data(ik + 752);
    const auto *ik_753 = buffer.data(ik + 753);
    const auto *ik_754 = buffer.data(ik + 754);
    const auto *ik_755 = buffer.data(ik + 755);
    const auto *ik_756 = buffer.data(ik + 756);
    const auto *ik_757 = buffer.data(ik + 757);
    const auto *ik_758 = buffer.data(ik + 758);
    const auto *ik_759 = buffer.data(ik + 759);
    const auto *ik_760 = buffer.data(ik + 760);
    const auto *ik_761 = buffer.data(ik + 761);
    const auto *ik_762 = buffer.data(ik + 762);
    const auto *ik_763 = buffer.data(ik + 763);
    const auto *ik_764 = buffer.data(ik + 764);
    const auto *ik_765 = buffer.data(ik + 765);
    const auto *ik_766 = buffer.data(ik + 766);
    const auto *ik_767 = buffer.data(ik + 767);
    const auto *ik_768 = buffer.data(ik + 768);
    const auto *ik_769 = buffer.data(ik + 769);
    const auto *ik_770 = buffer.data(ik + 770);
    const auto *ik_771 = buffer.data(ik + 771);
    const auto *ik_772 = buffer.data(ik + 772);
    const auto *ik_773 = buffer.data(ik + 773);
    const auto *ik_774 = buffer.data(ik + 774);
    const auto *ik_775 = buffer.data(ik + 775);
    const auto *ik_776 = buffer.data(ik + 776);
    const auto *ik_777 = buffer.data(ik + 777);
    const auto *ik_778 = buffer.data(ik + 778);
    const auto *ik_779 = buffer.data(ik + 779);
    const auto *ik_780 = buffer.data(ik + 780);
    const auto *ik_781 = buffer.data(ik + 781);
    const auto *ik_782 = buffer.data(ik + 782);
    const auto *ik_783 = buffer.data(ik + 783);
    const auto *ik_784 = buffer.data(ik + 784);
    const auto *ik_785 = buffer.data(ik + 785);
    const auto *ik_786 = buffer.data(ik + 786);
    const auto *ik_787 = buffer.data(ik + 787);
    const auto *ik_788 = buffer.data(ik + 788);
    const auto *ik_789 = buffer.data(ik + 789);
    const auto *ik_790 = buffer.data(ik + 790);
    const auto *ik_791 = buffer.data(ik + 791);
    const auto *ik_792 = buffer.data(ik + 792);
    const auto *ik_793 = buffer.data(ik + 793);
    const auto *ik_794 = buffer.data(ik + 794);
    const auto *ik_795 = buffer.data(ik + 795);
    const auto *ik_796 = buffer.data(ik + 796);
    const auto *ik_797 = buffer.data(ik + 797);
    const auto *ik_798 = buffer.data(ik + 798);
    const auto *ik_799 = buffer.data(ik + 799);
    const auto *ik_800 = buffer.data(ik + 800);
    const auto *ik_801 = buffer.data(ik + 801);
    const auto *ik_802 = buffer.data(ik + 802);
    const auto *ik_803 = buffer.data(ik + 803);
    const auto *ik_804 = buffer.data(ik + 804);
    const auto *ik_805 = buffer.data(ik + 805);
    const auto *ik_806 = buffer.data(ik + 806);
    const auto *ik_807 = buffer.data(ik + 807);
    const auto *ik_808 = buffer.data(ik + 808);
    const auto *ik_809 = buffer.data(ik + 809);
    const auto *ik_810 = buffer.data(ik + 810);
    const auto *ik_811 = buffer.data(ik + 811);
    const auto *ik_812 = buffer.data(ik + 812);
    const auto *ik_813 = buffer.data(ik + 813);
    const auto *ik_814 = buffer.data(ik + 814);
    const auto *ik_815 = buffer.data(ik + 815);
    const auto *ik_816 = buffer.data(ik + 816);
    const auto *ik_817 = buffer.data(ik + 817);
    const auto *ik_818 = buffer.data(ik + 818);
    const auto *ik_819 = buffer.data(ik + 819);
    const auto *ik_820 = buffer.data(ik + 820);
    const auto *ik_821 = buffer.data(ik + 821);
    const auto *ik_822 = buffer.data(ik + 822);
    const auto *ik_823 = buffer.data(ik + 823);
    const auto *ik_824 = buffer.data(ik + 824);
    const auto *ik_825 = buffer.data(ik + 825);
    const auto *ik_826 = buffer.data(ik + 826);
    const auto *ik_827 = buffer.data(ik + 827);
    const auto *ik_828 = buffer.data(ik + 828);
    const auto *ik_829 = buffer.data(ik + 829);
    const auto *ik_830 = buffer.data(ik + 830);
    const auto *ik_831 = buffer.data(ik + 831);
    const auto *ik_832 = buffer.data(ik + 832);
    const auto *ik_833 = buffer.data(ik + 833);
    const auto *ik_834 = buffer.data(ik + 834);
    const auto *ik_835 = buffer.data(ik + 835);
    const auto *ik_836 = buffer.data(ik + 836);
    const auto *ik_837 = buffer.data(ik + 837);
    const auto *ik_838 = buffer.data(ik + 838);
    const auto *ik_839 = buffer.data(ik + 839);
    const auto *ik_840 = buffer.data(ik + 840);
    const auto *ik_841 = buffer.data(ik + 841);
    const auto *ik_842 = buffer.data(ik + 842);
    const auto *ik_843 = buffer.data(ik + 843);
    const auto *ik_844 = buffer.data(ik + 844);
    const auto *ik_845 = buffer.data(ik + 845);
    const auto *ik_846 = buffer.data(ik + 846);
    const auto *ik_847 = buffer.data(ik + 847);
    const auto *ik_848 = buffer.data(ik + 848);
    const auto *ik_849 = buffer.data(ik + 849);
    const auto *ik_850 = buffer.data(ik + 850);
    const auto *ik_851 = buffer.data(ik + 851);
    const auto *ik_852 = buffer.data(ik + 852);
    const auto *ik_853 = buffer.data(ik + 853);
    const auto *ik_854 = buffer.data(ik + 854);
    const auto *ik_855 = buffer.data(ik + 855);
    const auto *ik_856 = buffer.data(ik + 856);
    const auto *ik_857 = buffer.data(ik + 857);
    const auto *ik_858 = buffer.data(ik + 858);
    const auto *ik_859 = buffer.data(ik + 859);
    const auto *ik_860 = buffer.data(ik + 860);
    const auto *ik_861 = buffer.data(ik + 861);
    const auto *ik_862 = buffer.data(ik + 862);
    const auto *ik_863 = buffer.data(ik + 863);
    const auto *ik_864 = buffer.data(ik + 864);
    const auto *ik_865 = buffer.data(ik + 865);
    const auto *ik_866 = buffer.data(ik + 866);
    const auto *ik_867 = buffer.data(ik + 867);
    const auto *ik_868 = buffer.data(ik + 868);
    const auto *ik_869 = buffer.data(ik + 869);
    const auto *ik_870 = buffer.data(ik + 870);
    const auto *ik_871 = buffer.data(ik + 871);
    const auto *ik_872 = buffer.data(ik + 872);
    const auto *ik_873 = buffer.data(ik + 873);
    const auto *ik_874 = buffer.data(ik + 874);
    const auto *ik_875 = buffer.data(ik + 875);
    const auto *ik_876 = buffer.data(ik + 876);
    const auto *ik_877 = buffer.data(ik + 877);
    const auto *ik_878 = buffer.data(ik + 878);
    const auto *ik_879 = buffer.data(ik + 879);
    const auto *ik_880 = buffer.data(ik + 880);
    const auto *ik_881 = buffer.data(ik + 881);
    const auto *ik_882 = buffer.data(ik + 882);
    const auto *ik_883 = buffer.data(ik + 883);
    const auto *ik_884 = buffer.data(ik + 884);
    const auto *ik_885 = buffer.data(ik + 885);
    const auto *ik_886 = buffer.data(ik + 886);
    const auto *ik_887 = buffer.data(ik + 887);
    const auto *ik_888 = buffer.data(ik + 888);
    const auto *ik_889 = buffer.data(ik + 889);
    const auto *ik_890 = buffer.data(ik + 890);
    const auto *ik_891 = buffer.data(ik + 891);
    const auto *ik_892 = buffer.data(ik + 892);
    const auto *ik_893 = buffer.data(ik + 893);
    const auto *ik_894 = buffer.data(ik + 894);
    const auto *ik_895 = buffer.data(ik + 895);
    const auto *ik_896 = buffer.data(ik + 896);
    const auto *ik_897 = buffer.data(ik + 897);
    const auto *ik_898 = buffer.data(ik + 898);
    const auto *ik_899 = buffer.data(ik + 899);

    const auto *lk_750 = buffer.data(lk + 750);
    const auto *lk_751 = buffer.data(lk + 751);
    const auto *lk_752 = buffer.data(lk + 752);
    const auto *lk_753 = buffer.data(lk + 753);
    const auto *lk_754 = buffer.data(lk + 754);
    const auto *lk_755 = buffer.data(lk + 755);
    const auto *lk_756 = buffer.data(lk + 756);
    const auto *lk_757 = buffer.data(lk + 757);
    const auto *lk_758 = buffer.data(lk + 758);
    const auto *lk_759 = buffer.data(lk + 759);
    const auto *lk_760 = buffer.data(lk + 760);
    const auto *lk_761 = buffer.data(lk + 761);
    const auto *lk_762 = buffer.data(lk + 762);
    const auto *lk_763 = buffer.data(lk + 763);
    const auto *lk_764 = buffer.data(lk + 764);
    const auto *lk_765 = buffer.data(lk + 765);
    const auto *lk_766 = buffer.data(lk + 766);
    const auto *lk_767 = buffer.data(lk + 767);
    const auto *lk_768 = buffer.data(lk + 768);
    const auto *lk_769 = buffer.data(lk + 769);
    const auto *lk_770 = buffer.data(lk + 770);
    const auto *lk_771 = buffer.data(lk + 771);
    const auto *lk_772 = buffer.data(lk + 772);
    const auto *lk_773 = buffer.data(lk + 773);
    const auto *lk_774 = buffer.data(lk + 774);
    const auto *lk_775 = buffer.data(lk + 775);
    const auto *lk_776 = buffer.data(lk + 776);
    const auto *lk_777 = buffer.data(lk + 777);
    const auto *lk_778 = buffer.data(lk + 778);
    const auto *lk_779 = buffer.data(lk + 779);
    const auto *lk_780 = buffer.data(lk + 780);
    const auto *lk_781 = buffer.data(lk + 781);
    const auto *lk_782 = buffer.data(lk + 782);
    const auto *lk_783 = buffer.data(lk + 783);
    const auto *lk_784 = buffer.data(lk + 784);
    const auto *lk_785 = buffer.data(lk + 785);
    const auto *lk_786 = buffer.data(lk + 786);
    const auto *lk_787 = buffer.data(lk + 787);
    const auto *lk_788 = buffer.data(lk + 788);
    const auto *lk_789 = buffer.data(lk + 789);
    const auto *lk_790 = buffer.data(lk + 790);
    const auto *lk_791 = buffer.data(lk + 791);
    const auto *lk_792 = buffer.data(lk + 792);
    const auto *lk_793 = buffer.data(lk + 793);
    const auto *lk_794 = buffer.data(lk + 794);
    const auto *lk_795 = buffer.data(lk + 795);
    const auto *lk_796 = buffer.data(lk + 796);
    const auto *lk_797 = buffer.data(lk + 797);
    const auto *lk_798 = buffer.data(lk + 798);
    const auto *lk_799 = buffer.data(lk + 799);
    const auto *lk_800 = buffer.data(lk + 800);
    const auto *lk_801 = buffer.data(lk + 801);
    const auto *lk_802 = buffer.data(lk + 802);
    const auto *lk_803 = buffer.data(lk + 803);
    const auto *lk_804 = buffer.data(lk + 804);
    const auto *lk_805 = buffer.data(lk + 805);
    const auto *lk_806 = buffer.data(lk + 806);
    const auto *lk_807 = buffer.data(lk + 807);
    const auto *lk_808 = buffer.data(lk + 808);
    const auto *lk_809 = buffer.data(lk + 809);
    const auto *lk_810 = buffer.data(lk + 810);
    const auto *lk_811 = buffer.data(lk + 811);
    const auto *lk_812 = buffer.data(lk + 812);
    const auto *lk_813 = buffer.data(lk + 813);
    const auto *lk_814 = buffer.data(lk + 814);
    const auto *lk_815 = buffer.data(lk + 815);
    const auto *lk_816 = buffer.data(lk + 816);
    const auto *lk_817 = buffer.data(lk + 817);
    const auto *lk_818 = buffer.data(lk + 818);
    const auto *lk_819 = buffer.data(lk + 819);
    const auto *lk_820 = buffer.data(lk + 820);
    const auto *lk_821 = buffer.data(lk + 821);
    const auto *lk_822 = buffer.data(lk + 822);
    const auto *lk_823 = buffer.data(lk + 823);
    const auto *lk_824 = buffer.data(lk + 824);
    const auto *lk_825 = buffer.data(lk + 825);
    const auto *lk_826 = buffer.data(lk + 826);
    const auto *lk_827 = buffer.data(lk + 827);
    const auto *lk_828 = buffer.data(lk + 828);
    const auto *lk_829 = buffer.data(lk + 829);
    const auto *lk_830 = buffer.data(lk + 830);
    const auto *lk_831 = buffer.data(lk + 831);
    const auto *lk_832 = buffer.data(lk + 832);
    const auto *lk_833 = buffer.data(lk + 833);
    const auto *lk_834 = buffer.data(lk + 834);
    const auto *lk_835 = buffer.data(lk + 835);
    const auto *lk_836 = buffer.data(lk + 836);
    const auto *lk_837 = buffer.data(lk + 837);
    const auto *lk_838 = buffer.data(lk + 838);
    const auto *lk_839 = buffer.data(lk + 839);
    const auto *lk_840 = buffer.data(lk + 840);
    const auto *lk_841 = buffer.data(lk + 841);
    const auto *lk_842 = buffer.data(lk + 842);
    const auto *lk_843 = buffer.data(lk + 843);
    const auto *lk_844 = buffer.data(lk + 844);
    const auto *lk_845 = buffer.data(lk + 845);
    const auto *lk_846 = buffer.data(lk + 846);
    const auto *lk_847 = buffer.data(lk + 847);
    const auto *lk_848 = buffer.data(lk + 848);
    const auto *lk_849 = buffer.data(lk + 849);
    const auto *lk_850 = buffer.data(lk + 850);
    const auto *lk_851 = buffer.data(lk + 851);
    const auto *lk_852 = buffer.data(lk + 852);
    const auto *lk_853 = buffer.data(lk + 853);
    const auto *lk_854 = buffer.data(lk + 854);
    const auto *lk_855 = buffer.data(lk + 855);
    const auto *lk_856 = buffer.data(lk + 856);
    const auto *lk_857 = buffer.data(lk + 857);
    const auto *lk_858 = buffer.data(lk + 858);
    const auto *lk_859 = buffer.data(lk + 859);
    const auto *lk_860 = buffer.data(lk + 860);
    const auto *lk_861 = buffer.data(lk + 861);
    const auto *lk_862 = buffer.data(lk + 862);
    const auto *lk_863 = buffer.data(lk + 863);
    const auto *lk_864 = buffer.data(lk + 864);
    const auto *lk_865 = buffer.data(lk + 865);
    const auto *lk_866 = buffer.data(lk + 866);
    const auto *lk_867 = buffer.data(lk + 867);
    const auto *lk_868 = buffer.data(lk + 868);
    const auto *lk_869 = buffer.data(lk + 869);
    const auto *lk_870 = buffer.data(lk + 870);
    const auto *lk_871 = buffer.data(lk + 871);
    const auto *lk_872 = buffer.data(lk + 872);
    const auto *lk_873 = buffer.data(lk + 873);
    const auto *lk_874 = buffer.data(lk + 874);
    const auto *lk_875 = buffer.data(lk + 875);
    const auto *lk_876 = buffer.data(lk + 876);
    const auto *lk_877 = buffer.data(lk + 877);
    const auto *lk_878 = buffer.data(lk + 878);
    const auto *lk_879 = buffer.data(lk + 879);
    const auto *lk_880 = buffer.data(lk + 880);
    const auto *lk_881 = buffer.data(lk + 881);
    const auto *lk_882 = buffer.data(lk + 882);
    const auto *lk_883 = buffer.data(lk + 883);
    const auto *lk_884 = buffer.data(lk + 884);
    const auto *lk_885 = buffer.data(lk + 885);
    const auto *lk_886 = buffer.data(lk + 886);
    const auto *lk_887 = buffer.data(lk + 887);
    const auto *lk_888 = buffer.data(lk + 888);
    const auto *lk_889 = buffer.data(lk + 889);
    const auto *lk_890 = buffer.data(lk + 890);
    const auto *lk_891 = buffer.data(lk + 891);
    const auto *lk_892 = buffer.data(lk + 892);
    const auto *lk_893 = buffer.data(lk + 893);
    const auto *lk_894 = buffer.data(lk + 894);
    const auto *lk_895 = buffer.data(lk + 895);
    const auto *lk_896 = buffer.data(lk + 896);
    const auto *lk_897 = buffer.data(lk + 897);
    const auto *lk_898 = buffer.data(lk + 898);
    const auto *lk_899 = buffer.data(lk + 899);

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, ik_750, ik_751, ik_752, ik_753, \
                         ik_754, lk_750, lk_751, lk_752, lk_753, \
                         lk_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -2.0 * ik_750[k]
                   + f_0 * lk_750[k];

        t_751[k] = -2.0 * ik_751[k]
                   + f_0 * lk_751[k];

        t_752[k] = -2.0 * ik_752[k]
                   + f_0 * lk_752[k];

        t_753[k] = -2.0 * ik_753[k]
                   + f_0 * lk_753[k];

        t_754[k] = -2.0 * ik_754[k]
                   + f_0 * lk_754[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, ik_755, ik_756, ik_757, ik_758, \
                         ik_759, lk_755, lk_756, lk_757, lk_758, \
                         lk_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -2.0 * ik_755[k]
                   + f_0 * lk_755[k];

        t_756[k] = -ik_756[k]
                   + f_0 * lk_756[k];

        t_757[k] = -ik_757[k]
                   + f_0 * lk_757[k];

        t_758[k] = -ik_758[k]
                   + f_0 * lk_758[k];

        t_759[k] = -ik_759[k]
                   + f_0 * lk_759[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, t_764, ik_760, ik_761, ik_762, ik_763, \
                         ik_764, lk_760, lk_761, lk_762, lk_763, \
                         lk_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = -ik_760[k]
                   + f_0 * lk_760[k];

        t_761[k] = -ik_761[k]
                   + f_0 * lk_761[k];

        t_762[k] = -ik_762[k]
                   + f_0 * lk_762[k];

        t_763[k] = -ik_763[k]
                   + f_0 * lk_763[k];

        t_764[k] = -ik_764[k]
                   + f_0 * lk_764[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, t_769, ik_765, ik_766, ik_767, ik_768, \
                         ik_769, lk_765, lk_766, lk_767, lk_768, \
                         lk_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = -ik_765[k]
                   + f_0 * lk_765[k];

        t_766[k] = -ik_766[k]
                   + f_0 * lk_766[k];

        t_767[k] = -ik_767[k]
                   + f_0 * lk_767[k];

        t_768[k] = -ik_768[k]
                   + f_0 * lk_768[k];

        t_769[k] = -ik_769[k]
                   + f_0 * lk_769[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, ik_770, ik_771, ik_772, ik_773, \
                         ik_774, lk_770, lk_771, lk_772, lk_773, \
                         lk_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = -ik_770[k]
                   + f_0 * lk_770[k];

        t_771[k] = -ik_771[k]
                   + f_0 * lk_771[k];

        t_772[k] = -ik_772[k]
                   + f_0 * lk_772[k];

        t_773[k] = -ik_773[k]
                   + f_0 * lk_773[k];

        t_774[k] = -ik_774[k]
                   + f_0 * lk_774[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, ik_775, ik_776, ik_777, ik_778, \
                         ik_779, lk_775, lk_776, lk_777, lk_778, \
                         lk_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = -ik_775[k]
                   + f_0 * lk_775[k];

        t_776[k] = -ik_776[k]
                   + f_0 * lk_776[k];

        t_777[k] = -ik_777[k]
                   + f_0 * lk_777[k];

        t_778[k] = -ik_778[k]
                   + f_0 * lk_778[k];

        t_779[k] = -ik_779[k]
                   + f_0 * lk_779[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, t_784, ik_780, ik_781, ik_782, ik_783, \
                         ik_784, lk_780, lk_781, lk_782, lk_783, \
                         lk_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = -ik_780[k]
                   + f_0 * lk_780[k];

        t_781[k] = -ik_781[k]
                   + f_0 * lk_781[k];

        t_782[k] = -ik_782[k]
                   + f_0 * lk_782[k];

        t_783[k] = -ik_783[k]
                   + f_0 * lk_783[k];

        t_784[k] = -ik_784[k]
                   + f_0 * lk_784[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, t_789, ik_785, ik_786, ik_787, ik_788, \
                         ik_789, lk_785, lk_786, lk_787, lk_788, \
                         lk_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = -ik_785[k]
                   + f_0 * lk_785[k];

        t_786[k] = -ik_786[k]
                   + f_0 * lk_786[k];

        t_787[k] = -ik_787[k]
                   + f_0 * lk_787[k];

        t_788[k] = -ik_788[k]
                   + f_0 * lk_788[k];

        t_789[k] = -ik_789[k]
                   + f_0 * lk_789[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, t_794, ik_790, ik_791, ik_792, ik_793, \
                         ik_794, lk_790, lk_791, lk_792, lk_793, \
                         lk_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = -ik_790[k]
                   + f_0 * lk_790[k];

        t_791[k] = -ik_791[k]
                   + f_0 * lk_791[k];

        t_792[k] = -ik_792[k]
                   + f_0 * lk_792[k];

        t_793[k] = -ik_793[k]
                   + f_0 * lk_793[k];

        t_794[k] = -ik_794[k]
                   + f_0 * lk_794[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, ik_795, ik_796, ik_797, ik_798, \
                         ik_799, lk_795, lk_796, lk_797, lk_798, \
                         lk_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = -ik_795[k]
                   + f_0 * lk_795[k];

        t_796[k] = -ik_796[k]
                   + f_0 * lk_796[k];

        t_797[k] = -ik_797[k]
                   + f_0 * lk_797[k];

        t_798[k] = -ik_798[k]
                   + f_0 * lk_798[k];

        t_799[k] = -ik_799[k]
                   + f_0 * lk_799[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, ik_800, ik_801, ik_802, ik_803, \
                         ik_804, lk_800, lk_801, lk_802, lk_803, \
                         lk_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = -ik_800[k]
                   + f_0 * lk_800[k];

        t_801[k] = -ik_801[k]
                   + f_0 * lk_801[k];

        t_802[k] = -ik_802[k]
                   + f_0 * lk_802[k];

        t_803[k] = -ik_803[k]
                   + f_0 * lk_803[k];

        t_804[k] = -ik_804[k]
                   + f_0 * lk_804[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, ik_805, ik_806, ik_807, ik_808, \
                         ik_809, lk_805, lk_806, lk_807, lk_808, \
                         lk_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = -ik_805[k]
                   + f_0 * lk_805[k];

        t_806[k] = -ik_806[k]
                   + f_0 * lk_806[k];

        t_807[k] = -ik_807[k]
                   + f_0 * lk_807[k];

        t_808[k] = -ik_808[k]
                   + f_0 * lk_808[k];

        t_809[k] = -ik_809[k]
                   + f_0 * lk_809[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, ik_810, ik_811, ik_812, ik_813, \
                         ik_814, lk_810, lk_811, lk_812, lk_813, \
                         lk_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = -ik_810[k]
                   + f_0 * lk_810[k];

        t_811[k] = -ik_811[k]
                   + f_0 * lk_811[k];

        t_812[k] = -ik_812[k]
                   + f_0 * lk_812[k];

        t_813[k] = -ik_813[k]
                   + f_0 * lk_813[k];

        t_814[k] = -ik_814[k]
                   + f_0 * lk_814[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, ik_815, ik_816, ik_817, ik_818, \
                         ik_819, lk_815, lk_816, lk_817, lk_818, \
                         lk_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = -ik_815[k]
                   + f_0 * lk_815[k];

        t_816[k] = -ik_816[k]
                   + f_0 * lk_816[k];

        t_817[k] = -ik_817[k]
                   + f_0 * lk_817[k];

        t_818[k] = -ik_818[k]
                   + f_0 * lk_818[k];

        t_819[k] = -ik_819[k]
                   + f_0 * lk_819[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, t_824, ik_820, ik_821, ik_822, ik_823, \
                         ik_824, lk_820, lk_821, lk_822, lk_823, \
                         lk_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = -ik_820[k]
                   + f_0 * lk_820[k];

        t_821[k] = -ik_821[k]
                   + f_0 * lk_821[k];

        t_822[k] = -ik_822[k]
                   + f_0 * lk_822[k];

        t_823[k] = -ik_823[k]
                   + f_0 * lk_823[k];

        t_824[k] = -ik_824[k]
                   + f_0 * lk_824[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, ik_825, ik_826, ik_827, ik_828, \
                         ik_829, lk_825, lk_826, lk_827, lk_828, \
                         lk_829 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = -ik_825[k]
                   + f_0 * lk_825[k];

        t_826[k] = -ik_826[k]
                   + f_0 * lk_826[k];

        t_827[k] = -ik_827[k]
                   + f_0 * lk_827[k];

        t_828[k] = -ik_828[k]
                   + f_0 * lk_828[k];

        t_829[k] = -ik_829[k]
                   + f_0 * lk_829[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, ik_830, ik_831, ik_832, ik_833, \
                         ik_834, lk_830, lk_831, lk_832, lk_833, \
                         lk_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = -ik_830[k]
                   + f_0 * lk_830[k];

        t_831[k] = -ik_831[k]
                   + f_0 * lk_831[k];

        t_832[k] = -ik_832[k]
                   + f_0 * lk_832[k];

        t_833[k] = -ik_833[k]
                   + f_0 * lk_833[k];

        t_834[k] = -ik_834[k]
                   + f_0 * lk_834[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, t_838, t_839, ik_835, ik_836, ik_837, ik_838, \
                         ik_839, lk_835, lk_836, lk_837, lk_838, \
                         lk_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = -ik_835[k]
                   + f_0 * lk_835[k];

        t_836[k] = -ik_836[k]
                   + f_0 * lk_836[k];

        t_837[k] = -ik_837[k]
                   + f_0 * lk_837[k];

        t_838[k] = -ik_838[k]
                   + f_0 * lk_838[k];

        t_839[k] = -ik_839[k]
                   + f_0 * lk_839[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, ik_840, ik_841, ik_842, ik_843, \
                         ik_844, lk_840, lk_841, lk_842, lk_843, \
                         lk_844 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = -ik_840[k]
                   + f_0 * lk_840[k];

        t_841[k] = -ik_841[k]
                   + f_0 * lk_841[k];

        t_842[k] = -ik_842[k]
                   + f_0 * lk_842[k];

        t_843[k] = -ik_843[k]
                   + f_0 * lk_843[k];

        t_844[k] = -ik_844[k]
                   + f_0 * lk_844[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, t_849, ik_845, ik_846, ik_847, ik_848, \
                         ik_849, lk_845, lk_846, lk_847, lk_848, \
                         lk_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = -ik_845[k]
                   + f_0 * lk_845[k];

        t_846[k] = -ik_846[k]
                   + f_0 * lk_846[k];

        t_847[k] = -ik_847[k]
                   + f_0 * lk_847[k];

        t_848[k] = -ik_848[k]
                   + f_0 * lk_848[k];

        t_849[k] = -ik_849[k]
                   + f_0 * lk_849[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, ik_850, ik_851, ik_852, ik_853, \
                         ik_854, lk_850, lk_851, lk_852, lk_853, \
                         lk_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = -ik_850[k]
                   + f_0 * lk_850[k];

        t_851[k] = -ik_851[k]
                   + f_0 * lk_851[k];

        t_852[k] = -ik_852[k]
                   + f_0 * lk_852[k];

        t_853[k] = -ik_853[k]
                   + f_0 * lk_853[k];

        t_854[k] = -ik_854[k]
                   + f_0 * lk_854[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, ik_855, ik_856, ik_857, ik_858, \
                         ik_859, lk_855, lk_856, lk_857, lk_858, \
                         lk_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = -ik_855[k]
                   + f_0 * lk_855[k];

        t_856[k] = -ik_856[k]
                   + f_0 * lk_856[k];

        t_857[k] = -ik_857[k]
                   + f_0 * lk_857[k];

        t_858[k] = -ik_858[k]
                   + f_0 * lk_858[k];

        t_859[k] = -ik_859[k]
                   + f_0 * lk_859[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, t_864, ik_860, ik_861, ik_862, ik_863, \
                         ik_864, lk_860, lk_861, lk_862, lk_863, \
                         lk_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = -ik_860[k]
                   + f_0 * lk_860[k];

        t_861[k] = -ik_861[k]
                   + f_0 * lk_861[k];

        t_862[k] = -ik_862[k]
                   + f_0 * lk_862[k];

        t_863[k] = -ik_863[k]
                   + f_0 * lk_863[k];

        t_864[k] = -ik_864[k]
                   + f_0 * lk_864[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, t_869, ik_865, ik_866, ik_867, ik_868, \
                         ik_869, lk_865, lk_866, lk_867, lk_868, \
                         lk_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = -ik_865[k]
                   + f_0 * lk_865[k];

        t_866[k] = -ik_866[k]
                   + f_0 * lk_866[k];

        t_867[k] = -ik_867[k]
                   + f_0 * lk_867[k];

        t_868[k] = -ik_868[k]
                   + f_0 * lk_868[k];

        t_869[k] = -ik_869[k]
                   + f_0 * lk_869[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, t_874, ik_870, ik_871, ik_872, ik_873, \
                         ik_874, lk_870, lk_871, lk_872, lk_873, \
                         lk_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = -ik_870[k]
                   + f_0 * lk_870[k];

        t_871[k] = -ik_871[k]
                   + f_0 * lk_871[k];

        t_872[k] = -ik_872[k]
                   + f_0 * lk_872[k];

        t_873[k] = -ik_873[k]
                   + f_0 * lk_873[k];

        t_874[k] = -ik_874[k]
                   + f_0 * lk_874[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, t_879, ik_875, ik_876, ik_877, ik_878, \
                         ik_879, lk_875, lk_876, lk_877, lk_878, \
                         lk_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = -ik_875[k]
                   + f_0 * lk_875[k];

        t_876[k] = -ik_876[k]
                   + f_0 * lk_876[k];

        t_877[k] = -ik_877[k]
                   + f_0 * lk_877[k];

        t_878[k] = -ik_878[k]
                   + f_0 * lk_878[k];

        t_879[k] = -ik_879[k]
                   + f_0 * lk_879[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, t_884, ik_880, ik_881, ik_882, ik_883, \
                         ik_884, lk_880, lk_881, lk_882, lk_883, \
                         lk_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = -ik_880[k]
                   + f_0 * lk_880[k];

        t_881[k] = -ik_881[k]
                   + f_0 * lk_881[k];

        t_882[k] = -ik_882[k]
                   + f_0 * lk_882[k];

        t_883[k] = -ik_883[k]
                   + f_0 * lk_883[k];

        t_884[k] = -ik_884[k]
                   + f_0 * lk_884[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, ik_885, ik_886, ik_887, ik_888, \
                         ik_889, lk_885, lk_886, lk_887, lk_888, \
                         lk_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = -ik_885[k]
                   + f_0 * lk_885[k];

        t_886[k] = -ik_886[k]
                   + f_0 * lk_886[k];

        t_887[k] = -ik_887[k]
                   + f_0 * lk_887[k];

        t_888[k] = -ik_888[k]
                   + f_0 * lk_888[k];

        t_889[k] = -ik_889[k]
                   + f_0 * lk_889[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, ik_890, ik_891, ik_892, ik_893, \
                         ik_894, lk_890, lk_891, lk_892, lk_893, \
                         lk_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = -ik_890[k]
                   + f_0 * lk_890[k];

        t_891[k] = -ik_891[k]
                   + f_0 * lk_891[k];

        t_892[k] = -ik_892[k]
                   + f_0 * lk_892[k];

        t_893[k] = -ik_893[k]
                   + f_0 * lk_893[k];

        t_894[k] = -ik_894[k]
                   + f_0 * lk_894[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, t_899, ik_895, ik_896, ik_897, ik_898, \
                         ik_899, lk_895, lk_896, lk_897, lk_898, \
                         lk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = -ik_895[k]
                   + f_0 * lk_895[k];

        t_896[k] = -ik_896[k]
                   + f_0 * lk_896[k];

        t_897[k] = -ik_897[k]
                   + f_0 * lk_897[k];

        t_898[k] = -ik_898[k]
                   + f_0 * lk_898[k];

        t_899[k] = -ik_899[k]
                   + f_0 * lk_899[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_0_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_900 = buffer.data(ik + 900);
    const auto *ik_901 = buffer.data(ik + 901);
    const auto *ik_902 = buffer.data(ik + 902);
    const auto *ik_903 = buffer.data(ik + 903);
    const auto *ik_904 = buffer.data(ik + 904);
    const auto *ik_905 = buffer.data(ik + 905);
    const auto *ik_906 = buffer.data(ik + 906);
    const auto *ik_907 = buffer.data(ik + 907);
    const auto *ik_908 = buffer.data(ik + 908);
    const auto *ik_909 = buffer.data(ik + 909);
    const auto *ik_910 = buffer.data(ik + 910);
    const auto *ik_911 = buffer.data(ik + 911);
    const auto *ik_912 = buffer.data(ik + 912);
    const auto *ik_913 = buffer.data(ik + 913);
    const auto *ik_914 = buffer.data(ik + 914);
    const auto *ik_915 = buffer.data(ik + 915);
    const auto *ik_916 = buffer.data(ik + 916);
    const auto *ik_917 = buffer.data(ik + 917);
    const auto *ik_918 = buffer.data(ik + 918);
    const auto *ik_919 = buffer.data(ik + 919);
    const auto *ik_920 = buffer.data(ik + 920);
    const auto *ik_921 = buffer.data(ik + 921);
    const auto *ik_922 = buffer.data(ik + 922);
    const auto *ik_923 = buffer.data(ik + 923);
    const auto *ik_924 = buffer.data(ik + 924);
    const auto *ik_925 = buffer.data(ik + 925);
    const auto *ik_926 = buffer.data(ik + 926);
    const auto *ik_927 = buffer.data(ik + 927);
    const auto *ik_928 = buffer.data(ik + 928);
    const auto *ik_929 = buffer.data(ik + 929);
    const auto *ik_930 = buffer.data(ik + 930);
    const auto *ik_931 = buffer.data(ik + 931);
    const auto *ik_932 = buffer.data(ik + 932);
    const auto *ik_933 = buffer.data(ik + 933);
    const auto *ik_934 = buffer.data(ik + 934);
    const auto *ik_935 = buffer.data(ik + 935);
    const auto *ik_936 = buffer.data(ik + 936);
    const auto *ik_937 = buffer.data(ik + 937);
    const auto *ik_938 = buffer.data(ik + 938);
    const auto *ik_939 = buffer.data(ik + 939);
    const auto *ik_940 = buffer.data(ik + 940);
    const auto *ik_941 = buffer.data(ik + 941);
    const auto *ik_942 = buffer.data(ik + 942);
    const auto *ik_943 = buffer.data(ik + 943);
    const auto *ik_944 = buffer.data(ik + 944);
    const auto *ik_945 = buffer.data(ik + 945);
    const auto *ik_946 = buffer.data(ik + 946);
    const auto *ik_947 = buffer.data(ik + 947);
    const auto *ik_948 = buffer.data(ik + 948);
    const auto *ik_949 = buffer.data(ik + 949);
    const auto *ik_950 = buffer.data(ik + 950);
    const auto *ik_951 = buffer.data(ik + 951);
    const auto *ik_952 = buffer.data(ik + 952);
    const auto *ik_953 = buffer.data(ik + 953);
    const auto *ik_954 = buffer.data(ik + 954);
    const auto *ik_955 = buffer.data(ik + 955);
    const auto *ik_956 = buffer.data(ik + 956);
    const auto *ik_957 = buffer.data(ik + 957);
    const auto *ik_958 = buffer.data(ik + 958);
    const auto *ik_959 = buffer.data(ik + 959);
    const auto *ik_960 = buffer.data(ik + 960);
    const auto *ik_961 = buffer.data(ik + 961);
    const auto *ik_962 = buffer.data(ik + 962);
    const auto *ik_963 = buffer.data(ik + 963);
    const auto *ik_964 = buffer.data(ik + 964);
    const auto *ik_965 = buffer.data(ik + 965);
    const auto *ik_966 = buffer.data(ik + 966);
    const auto *ik_967 = buffer.data(ik + 967);
    const auto *ik_968 = buffer.data(ik + 968);
    const auto *ik_969 = buffer.data(ik + 969);
    const auto *ik_970 = buffer.data(ik + 970);
    const auto *ik_971 = buffer.data(ik + 971);
    const auto *ik_972 = buffer.data(ik + 972);
    const auto *ik_973 = buffer.data(ik + 973);
    const auto *ik_974 = buffer.data(ik + 974);
    const auto *ik_975 = buffer.data(ik + 975);
    const auto *ik_976 = buffer.data(ik + 976);
    const auto *ik_977 = buffer.data(ik + 977);
    const auto *ik_978 = buffer.data(ik + 978);
    const auto *ik_979 = buffer.data(ik + 979);
    const auto *ik_980 = buffer.data(ik + 980);
    const auto *ik_981 = buffer.data(ik + 981);
    const auto *ik_982 = buffer.data(ik + 982);
    const auto *ik_983 = buffer.data(ik + 983);
    const auto *ik_984 = buffer.data(ik + 984);
    const auto *ik_985 = buffer.data(ik + 985);
    const auto *ik_986 = buffer.data(ik + 986);
    const auto *ik_987 = buffer.data(ik + 987);
    const auto *ik_988 = buffer.data(ik + 988);
    const auto *ik_989 = buffer.data(ik + 989);
    const auto *ik_990 = buffer.data(ik + 990);
    const auto *ik_991 = buffer.data(ik + 991);
    const auto *ik_992 = buffer.data(ik + 992);
    const auto *ik_993 = buffer.data(ik + 993);
    const auto *ik_994 = buffer.data(ik + 994);
    const auto *ik_995 = buffer.data(ik + 995);
    const auto *ik_996 = buffer.data(ik + 996);
    const auto *ik_997 = buffer.data(ik + 997);
    const auto *ik_998 = buffer.data(ik + 998);
    const auto *ik_999 = buffer.data(ik + 999);
    const auto *ik_1000 = buffer.data(ik + 1000);
    const auto *ik_1001 = buffer.data(ik + 1001);
    const auto *ik_1002 = buffer.data(ik + 1002);
    const auto *ik_1003 = buffer.data(ik + 1003);
    const auto *ik_1004 = buffer.data(ik + 1004);
    const auto *ik_1005 = buffer.data(ik + 1005);
    const auto *ik_1006 = buffer.data(ik + 1006);
    const auto *ik_1007 = buffer.data(ik + 1007);

    const auto *lk_900 = buffer.data(lk + 900);
    const auto *lk_901 = buffer.data(lk + 901);
    const auto *lk_902 = buffer.data(lk + 902);
    const auto *lk_903 = buffer.data(lk + 903);
    const auto *lk_904 = buffer.data(lk + 904);
    const auto *lk_905 = buffer.data(lk + 905);
    const auto *lk_906 = buffer.data(lk + 906);
    const auto *lk_907 = buffer.data(lk + 907);
    const auto *lk_908 = buffer.data(lk + 908);
    const auto *lk_909 = buffer.data(lk + 909);
    const auto *lk_910 = buffer.data(lk + 910);
    const auto *lk_911 = buffer.data(lk + 911);
    const auto *lk_912 = buffer.data(lk + 912);
    const auto *lk_913 = buffer.data(lk + 913);
    const auto *lk_914 = buffer.data(lk + 914);
    const auto *lk_915 = buffer.data(lk + 915);
    const auto *lk_916 = buffer.data(lk + 916);
    const auto *lk_917 = buffer.data(lk + 917);
    const auto *lk_918 = buffer.data(lk + 918);
    const auto *lk_919 = buffer.data(lk + 919);
    const auto *lk_920 = buffer.data(lk + 920);
    const auto *lk_921 = buffer.data(lk + 921);
    const auto *lk_922 = buffer.data(lk + 922);
    const auto *lk_923 = buffer.data(lk + 923);
    const auto *lk_924 = buffer.data(lk + 924);
    const auto *lk_925 = buffer.data(lk + 925);
    const auto *lk_926 = buffer.data(lk + 926);
    const auto *lk_927 = buffer.data(lk + 927);
    const auto *lk_928 = buffer.data(lk + 928);
    const auto *lk_929 = buffer.data(lk + 929);
    const auto *lk_930 = buffer.data(lk + 930);
    const auto *lk_931 = buffer.data(lk + 931);
    const auto *lk_932 = buffer.data(lk + 932);
    const auto *lk_933 = buffer.data(lk + 933);
    const auto *lk_934 = buffer.data(lk + 934);
    const auto *lk_935 = buffer.data(lk + 935);
    const auto *lk_936 = buffer.data(lk + 936);
    const auto *lk_937 = buffer.data(lk + 937);
    const auto *lk_938 = buffer.data(lk + 938);
    const auto *lk_939 = buffer.data(lk + 939);
    const auto *lk_940 = buffer.data(lk + 940);
    const auto *lk_941 = buffer.data(lk + 941);
    const auto *lk_942 = buffer.data(lk + 942);
    const auto *lk_943 = buffer.data(lk + 943);
    const auto *lk_944 = buffer.data(lk + 944);
    const auto *lk_945 = buffer.data(lk + 945);
    const auto *lk_946 = buffer.data(lk + 946);
    const auto *lk_947 = buffer.data(lk + 947);
    const auto *lk_948 = buffer.data(lk + 948);
    const auto *lk_949 = buffer.data(lk + 949);
    const auto *lk_950 = buffer.data(lk + 950);
    const auto *lk_951 = buffer.data(lk + 951);
    const auto *lk_952 = buffer.data(lk + 952);
    const auto *lk_953 = buffer.data(lk + 953);
    const auto *lk_954 = buffer.data(lk + 954);
    const auto *lk_955 = buffer.data(lk + 955);
    const auto *lk_956 = buffer.data(lk + 956);
    const auto *lk_957 = buffer.data(lk + 957);
    const auto *lk_958 = buffer.data(lk + 958);
    const auto *lk_959 = buffer.data(lk + 959);
    const auto *lk_960 = buffer.data(lk + 960);
    const auto *lk_961 = buffer.data(lk + 961);
    const auto *lk_962 = buffer.data(lk + 962);
    const auto *lk_963 = buffer.data(lk + 963);
    const auto *lk_964 = buffer.data(lk + 964);
    const auto *lk_965 = buffer.data(lk + 965);
    const auto *lk_966 = buffer.data(lk + 966);
    const auto *lk_967 = buffer.data(lk + 967);
    const auto *lk_968 = buffer.data(lk + 968);
    const auto *lk_969 = buffer.data(lk + 969);
    const auto *lk_970 = buffer.data(lk + 970);
    const auto *lk_971 = buffer.data(lk + 971);
    const auto *lk_972 = buffer.data(lk + 972);
    const auto *lk_973 = buffer.data(lk + 973);
    const auto *lk_974 = buffer.data(lk + 974);
    const auto *lk_975 = buffer.data(lk + 975);
    const auto *lk_976 = buffer.data(lk + 976);
    const auto *lk_977 = buffer.data(lk + 977);
    const auto *lk_978 = buffer.data(lk + 978);
    const auto *lk_979 = buffer.data(lk + 979);
    const auto *lk_980 = buffer.data(lk + 980);
    const auto *lk_981 = buffer.data(lk + 981);
    const auto *lk_982 = buffer.data(lk + 982);
    const auto *lk_983 = buffer.data(lk + 983);
    const auto *lk_984 = buffer.data(lk + 984);
    const auto *lk_985 = buffer.data(lk + 985);
    const auto *lk_986 = buffer.data(lk + 986);
    const auto *lk_987 = buffer.data(lk + 987);
    const auto *lk_988 = buffer.data(lk + 988);
    const auto *lk_989 = buffer.data(lk + 989);
    const auto *lk_990 = buffer.data(lk + 990);
    const auto *lk_991 = buffer.data(lk + 991);
    const auto *lk_992 = buffer.data(lk + 992);
    const auto *lk_993 = buffer.data(lk + 993);
    const auto *lk_994 = buffer.data(lk + 994);
    const auto *lk_995 = buffer.data(lk + 995);
    const auto *lk_996 = buffer.data(lk + 996);
    const auto *lk_997 = buffer.data(lk + 997);
    const auto *lk_998 = buffer.data(lk + 998);
    const auto *lk_999 = buffer.data(lk + 999);
    const auto *lk_1000 = buffer.data(lk + 1000);
    const auto *lk_1001 = buffer.data(lk + 1001);
    const auto *lk_1002 = buffer.data(lk + 1002);
    const auto *lk_1003 = buffer.data(lk + 1003);
    const auto *lk_1004 = buffer.data(lk + 1004);
    const auto *lk_1005 = buffer.data(lk + 1005);
    const auto *lk_1006 = buffer.data(lk + 1006);
    const auto *lk_1007 = buffer.data(lk + 1007);
    const auto *lk_1008 = buffer.data(lk + 1008);
    const auto *lk_1009 = buffer.data(lk + 1009);
    const auto *lk_1010 = buffer.data(lk + 1010);
    const auto *lk_1011 = buffer.data(lk + 1011);
    const auto *lk_1012 = buffer.data(lk + 1012);
    const auto *lk_1013 = buffer.data(lk + 1013);
    const auto *lk_1014 = buffer.data(lk + 1014);
    const auto *lk_1015 = buffer.data(lk + 1015);
    const auto *lk_1016 = buffer.data(lk + 1016);
    const auto *lk_1017 = buffer.data(lk + 1017);
    const auto *lk_1018 = buffer.data(lk + 1018);
    const auto *lk_1019 = buffer.data(lk + 1019);
    const auto *lk_1020 = buffer.data(lk + 1020);
    const auto *lk_1021 = buffer.data(lk + 1021);
    const auto *lk_1022 = buffer.data(lk + 1022);
    const auto *lk_1023 = buffer.data(lk + 1023);
    const auto *lk_1024 = buffer.data(lk + 1024);
    const auto *lk_1025 = buffer.data(lk + 1025);
    const auto *lk_1026 = buffer.data(lk + 1026);
    const auto *lk_1027 = buffer.data(lk + 1027);
    const auto *lk_1028 = buffer.data(lk + 1028);
    const auto *lk_1029 = buffer.data(lk + 1029);
    const auto *lk_1030 = buffer.data(lk + 1030);
    const auto *lk_1031 = buffer.data(lk + 1031);
    const auto *lk_1032 = buffer.data(lk + 1032);
    const auto *lk_1033 = buffer.data(lk + 1033);
    const auto *lk_1034 = buffer.data(lk + 1034);
    const auto *lk_1035 = buffer.data(lk + 1035);
    const auto *lk_1036 = buffer.data(lk + 1036);
    const auto *lk_1037 = buffer.data(lk + 1037);
    const auto *lk_1038 = buffer.data(lk + 1038);
    const auto *lk_1039 = buffer.data(lk + 1039);
    const auto *lk_1040 = buffer.data(lk + 1040);
    const auto *lk_1041 = buffer.data(lk + 1041);
    const auto *lk_1042 = buffer.data(lk + 1042);
    const auto *lk_1043 = buffer.data(lk + 1043);
    const auto *lk_1044 = buffer.data(lk + 1044);
    const auto *lk_1045 = buffer.data(lk + 1045);
    const auto *lk_1046 = buffer.data(lk + 1046);
    const auto *lk_1047 = buffer.data(lk + 1047);
    const auto *lk_1048 = buffer.data(lk + 1048);
    const auto *lk_1049 = buffer.data(lk + 1049);
    const auto *lk_1050 = buffer.data(lk + 1050);
    const auto *lk_1051 = buffer.data(lk + 1051);
    const auto *lk_1052 = buffer.data(lk + 1052);
    const auto *lk_1053 = buffer.data(lk + 1053);
    const auto *lk_1054 = buffer.data(lk + 1054);
    const auto *lk_1055 = buffer.data(lk + 1055);
    const auto *lk_1056 = buffer.data(lk + 1056);
    const auto *lk_1057 = buffer.data(lk + 1057);
    const auto *lk_1058 = buffer.data(lk + 1058);
    const auto *lk_1059 = buffer.data(lk + 1059);
    const auto *lk_1060 = buffer.data(lk + 1060);
    const auto *lk_1061 = buffer.data(lk + 1061);
    const auto *lk_1062 = buffer.data(lk + 1062);
    const auto *lk_1063 = buffer.data(lk + 1063);
    const auto *lk_1064 = buffer.data(lk + 1064);
    const auto *lk_1065 = buffer.data(lk + 1065);
    const auto *lk_1066 = buffer.data(lk + 1066);

#pragma omp simd aligned(t_900, t_901, t_902, t_903, t_904, ik_900, ik_901, ik_902, ik_903, \
                         ik_904, lk_900, lk_901, lk_902, lk_903, \
                         lk_904 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = -ik_900[k]
                   + f_0 * lk_900[k];

        t_901[k] = -ik_901[k]
                   + f_0 * lk_901[k];

        t_902[k] = -ik_902[k]
                   + f_0 * lk_902[k];

        t_903[k] = -ik_903[k]
                   + f_0 * lk_903[k];

        t_904[k] = -ik_904[k]
                   + f_0 * lk_904[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, t_909, ik_905, ik_906, ik_907, ik_908, \
                         ik_909, lk_905, lk_906, lk_907, lk_908, \
                         lk_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = -ik_905[k]
                   + f_0 * lk_905[k];

        t_906[k] = -ik_906[k]
                   + f_0 * lk_906[k];

        t_907[k] = -ik_907[k]
                   + f_0 * lk_907[k];

        t_908[k] = -ik_908[k]
                   + f_0 * lk_908[k];

        t_909[k] = -ik_909[k]
                   + f_0 * lk_909[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, t_914, ik_910, ik_911, ik_912, ik_913, \
                         ik_914, lk_910, lk_911, lk_912, lk_913, \
                         lk_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = -ik_910[k]
                   + f_0 * lk_910[k];

        t_911[k] = -ik_911[k]
                   + f_0 * lk_911[k];

        t_912[k] = -ik_912[k]
                   + f_0 * lk_912[k];

        t_913[k] = -ik_913[k]
                   + f_0 * lk_913[k];

        t_914[k] = -ik_914[k]
                   + f_0 * lk_914[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, t_919, ik_915, ik_916, ik_917, ik_918, \
                         ik_919, lk_915, lk_916, lk_917, lk_918, \
                         lk_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = -ik_915[k]
                   + f_0 * lk_915[k];

        t_916[k] = -ik_916[k]
                   + f_0 * lk_916[k];

        t_917[k] = -ik_917[k]
                   + f_0 * lk_917[k];

        t_918[k] = -ik_918[k]
                   + f_0 * lk_918[k];

        t_919[k] = -ik_919[k]
                   + f_0 * lk_919[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, t_924, ik_920, ik_921, ik_922, ik_923, \
                         ik_924, lk_920, lk_921, lk_922, lk_923, \
                         lk_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = -ik_920[k]
                   + f_0 * lk_920[k];

        t_921[k] = -ik_921[k]
                   + f_0 * lk_921[k];

        t_922[k] = -ik_922[k]
                   + f_0 * lk_922[k];

        t_923[k] = -ik_923[k]
                   + f_0 * lk_923[k];

        t_924[k] = -ik_924[k]
                   + f_0 * lk_924[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, t_929, ik_925, ik_926, ik_927, ik_928, \
                         ik_929, lk_925, lk_926, lk_927, lk_928, \
                         lk_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = -ik_925[k]
                   + f_0 * lk_925[k];

        t_926[k] = -ik_926[k]
                   + f_0 * lk_926[k];

        t_927[k] = -ik_927[k]
                   + f_0 * lk_927[k];

        t_928[k] = -ik_928[k]
                   + f_0 * lk_928[k];

        t_929[k] = -ik_929[k]
                   + f_0 * lk_929[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, ik_930, ik_931, ik_932, ik_933, \
                         ik_934, lk_930, lk_931, lk_932, lk_933, \
                         lk_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = -ik_930[k]
                   + f_0 * lk_930[k];

        t_931[k] = -ik_931[k]
                   + f_0 * lk_931[k];

        t_932[k] = -ik_932[k]
                   + f_0 * lk_932[k];

        t_933[k] = -ik_933[k]
                   + f_0 * lk_933[k];

        t_934[k] = -ik_934[k]
                   + f_0 * lk_934[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, t_939, ik_935, ik_936, ik_937, ik_938, \
                         ik_939, lk_935, lk_936, lk_937, lk_938, \
                         lk_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = -ik_935[k]
                   + f_0 * lk_935[k];

        t_936[k] = -ik_936[k]
                   + f_0 * lk_936[k];

        t_937[k] = -ik_937[k]
                   + f_0 * lk_937[k];

        t_938[k] = -ik_938[k]
                   + f_0 * lk_938[k];

        t_939[k] = -ik_939[k]
                   + f_0 * lk_939[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, t_944, ik_940, ik_941, ik_942, ik_943, \
                         ik_944, lk_940, lk_941, lk_942, lk_943, \
                         lk_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = -ik_940[k]
                   + f_0 * lk_940[k];

        t_941[k] = -ik_941[k]
                   + f_0 * lk_941[k];

        t_942[k] = -ik_942[k]
                   + f_0 * lk_942[k];

        t_943[k] = -ik_943[k]
                   + f_0 * lk_943[k];

        t_944[k] = -ik_944[k]
                   + f_0 * lk_944[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, t_949, ik_945, ik_946, ik_947, ik_948, \
                         ik_949, lk_945, lk_946, lk_947, lk_948, \
                         lk_949 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = -ik_945[k]
                   + f_0 * lk_945[k];

        t_946[k] = -ik_946[k]
                   + f_0 * lk_946[k];

        t_947[k] = -ik_947[k]
                   + f_0 * lk_947[k];

        t_948[k] = -ik_948[k]
                   + f_0 * lk_948[k];

        t_949[k] = -ik_949[k]
                   + f_0 * lk_949[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, t_954, ik_950, ik_951, ik_952, ik_953, \
                         ik_954, lk_950, lk_951, lk_952, lk_953, \
                         lk_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = -ik_950[k]
                   + f_0 * lk_950[k];

        t_951[k] = -ik_951[k]
                   + f_0 * lk_951[k];

        t_952[k] = -ik_952[k]
                   + f_0 * lk_952[k];

        t_953[k] = -ik_953[k]
                   + f_0 * lk_953[k];

        t_954[k] = -ik_954[k]
                   + f_0 * lk_954[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, t_959, ik_955, ik_956, ik_957, ik_958, \
                         ik_959, lk_955, lk_956, lk_957, lk_958, \
                         lk_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = -ik_955[k]
                   + f_0 * lk_955[k];

        t_956[k] = -ik_956[k]
                   + f_0 * lk_956[k];

        t_957[k] = -ik_957[k]
                   + f_0 * lk_957[k];

        t_958[k] = -ik_958[k]
                   + f_0 * lk_958[k];

        t_959[k] = -ik_959[k]
                   + f_0 * lk_959[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, t_963, t_964, ik_960, ik_961, ik_962, ik_963, \
                         ik_964, lk_960, lk_961, lk_962, lk_963, \
                         lk_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = -ik_960[k]
                   + f_0 * lk_960[k];

        t_961[k] = -ik_961[k]
                   + f_0 * lk_961[k];

        t_962[k] = -ik_962[k]
                   + f_0 * lk_962[k];

        t_963[k] = -ik_963[k]
                   + f_0 * lk_963[k];

        t_964[k] = -ik_964[k]
                   + f_0 * lk_964[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, t_969, ik_965, ik_966, ik_967, ik_968, \
                         ik_969, lk_965, lk_966, lk_967, lk_968, \
                         lk_969 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = -ik_965[k]
                   + f_0 * lk_965[k];

        t_966[k] = -ik_966[k]
                   + f_0 * lk_966[k];

        t_967[k] = -ik_967[k]
                   + f_0 * lk_967[k];

        t_968[k] = -ik_968[k]
                   + f_0 * lk_968[k];

        t_969[k] = -ik_969[k]
                   + f_0 * lk_969[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, t_974, ik_970, ik_971, ik_972, ik_973, \
                         ik_974, lk_970, lk_971, lk_972, lk_973, \
                         lk_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = -ik_970[k]
                   + f_0 * lk_970[k];

        t_971[k] = -ik_971[k]
                   + f_0 * lk_971[k];

        t_972[k] = -ik_972[k]
                   + f_0 * lk_972[k];

        t_973[k] = -ik_973[k]
                   + f_0 * lk_973[k];

        t_974[k] = -ik_974[k]
                   + f_0 * lk_974[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, t_979, ik_975, ik_976, ik_977, ik_978, \
                         ik_979, lk_975, lk_976, lk_977, lk_978, \
                         lk_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = -ik_975[k]
                   + f_0 * lk_975[k];

        t_976[k] = -ik_976[k]
                   + f_0 * lk_976[k];

        t_977[k] = -ik_977[k]
                   + f_0 * lk_977[k];

        t_978[k] = -ik_978[k]
                   + f_0 * lk_978[k];

        t_979[k] = -ik_979[k]
                   + f_0 * lk_979[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, t_984, ik_980, ik_981, ik_982, ik_983, \
                         ik_984, lk_980, lk_981, lk_982, lk_983, \
                         lk_984 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = -ik_980[k]
                   + f_0 * lk_980[k];

        t_981[k] = -ik_981[k]
                   + f_0 * lk_981[k];

        t_982[k] = -ik_982[k]
                   + f_0 * lk_982[k];

        t_983[k] = -ik_983[k]
                   + f_0 * lk_983[k];

        t_984[k] = -ik_984[k]
                   + f_0 * lk_984[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, t_989, ik_985, ik_986, ik_987, ik_988, \
                         ik_989, lk_985, lk_986, lk_987, lk_988, \
                         lk_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = -ik_985[k]
                   + f_0 * lk_985[k];

        t_986[k] = -ik_986[k]
                   + f_0 * lk_986[k];

        t_987[k] = -ik_987[k]
                   + f_0 * lk_987[k];

        t_988[k] = -ik_988[k]
                   + f_0 * lk_988[k];

        t_989[k] = -ik_989[k]
                   + f_0 * lk_989[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, t_993, t_994, ik_990, ik_991, ik_992, ik_993, \
                         ik_994, lk_990, lk_991, lk_992, lk_993, \
                         lk_994 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = -ik_990[k]
                   + f_0 * lk_990[k];

        t_991[k] = -ik_991[k]
                   + f_0 * lk_991[k];

        t_992[k] = -ik_992[k]
                   + f_0 * lk_992[k];

        t_993[k] = -ik_993[k]
                   + f_0 * lk_993[k];

        t_994[k] = -ik_994[k]
                   + f_0 * lk_994[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, ik_995, ik_996, ik_997, ik_998, \
                         ik_999, lk_995, lk_996, lk_997, lk_998, \
                         lk_999 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = -ik_995[k]
                   + f_0 * lk_995[k];

        t_996[k] = -ik_996[k]
                   + f_0 * lk_996[k];

        t_997[k] = -ik_997[k]
                   + f_0 * lk_997[k];

        t_998[k] = -ik_998[k]
                   + f_0 * lk_998[k];

        t_999[k] = -ik_999[k]
                   + f_0 * lk_999[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, t_1004, ik_1000, ik_1001, ik_1002, \
                         ik_1003, ik_1004, lk_1000, lk_1001, lk_1002, lk_1003, \
                         lk_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = -ik_1000[k]
                    + f_0 * lk_1000[k];

        t_1001[k] = -ik_1001[k]
                    + f_0 * lk_1001[k];

        t_1002[k] = -ik_1002[k]
                    + f_0 * lk_1002[k];

        t_1003[k] = -ik_1003[k]
                    + f_0 * lk_1003[k];

        t_1004[k] = -ik_1004[k]
                    + f_0 * lk_1004[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, t_1008, t_1009, t_1010, ik_1005, ik_1006, \
                         ik_1007, lk_1005, lk_1006, lk_1007, lk_1008, lk_1009, \
                         lk_1010 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = -ik_1005[k]
                    + f_0 * lk_1005[k];

        t_1006[k] = -ik_1006[k]
                    + f_0 * lk_1006[k];

        t_1007[k] = -ik_1007[k]
                    + f_0 * lk_1007[k];

        t_1008[k] = f_0 * lk_1008[k];

        t_1009[k] = f_0 * lk_1009[k];

        t_1010[k] = f_0 * lk_1010[k];
    }

#pragma omp simd aligned(t_1011, t_1012, t_1013, t_1014, t_1015, t_1016, t_1017, t_1018, \
                         lk_1011, lk_1012, lk_1013, lk_1014, lk_1015, lk_1016, lk_1017, \
                         lk_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1011[k] = f_0 * lk_1011[k];

        t_1012[k] = f_0 * lk_1012[k];

        t_1013[k] = f_0 * lk_1013[k];

        t_1014[k] = f_0 * lk_1014[k];

        t_1015[k] = f_0 * lk_1015[k];

        t_1016[k] = f_0 * lk_1016[k];

        t_1017[k] = f_0 * lk_1017[k];

        t_1018[k] = f_0 * lk_1018[k];
    }

#pragma omp simd aligned(t_1019, t_1020, t_1021, t_1022, t_1023, t_1024, t_1025, t_1026, \
                         lk_1019, lk_1020, lk_1021, lk_1022, lk_1023, lk_1024, lk_1025, \
                         lk_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1019[k] = f_0 * lk_1019[k];

        t_1020[k] = f_0 * lk_1020[k];

        t_1021[k] = f_0 * lk_1021[k];

        t_1022[k] = f_0 * lk_1022[k];

        t_1023[k] = f_0 * lk_1023[k];

        t_1024[k] = f_0 * lk_1024[k];

        t_1025[k] = f_0 * lk_1025[k];

        t_1026[k] = f_0 * lk_1026[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, t_1031, t_1032, t_1033, t_1034, \
                         lk_1027, lk_1028, lk_1029, lk_1030, lk_1031, lk_1032, lk_1033, \
                         lk_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_0 * lk_1027[k];

        t_1028[k] = f_0 * lk_1028[k];

        t_1029[k] = f_0 * lk_1029[k];

        t_1030[k] = f_0 * lk_1030[k];

        t_1031[k] = f_0 * lk_1031[k];

        t_1032[k] = f_0 * lk_1032[k];

        t_1033[k] = f_0 * lk_1033[k];

        t_1034[k] = f_0 * lk_1034[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, t_1038, t_1039, t_1040, t_1041, t_1042, \
                         lk_1035, lk_1036, lk_1037, lk_1038, lk_1039, lk_1040, lk_1041, \
                         lk_1042 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = f_0 * lk_1035[k];

        t_1036[k] = f_0 * lk_1036[k];

        t_1037[k] = f_0 * lk_1037[k];

        t_1038[k] = f_0 * lk_1038[k];

        t_1039[k] = f_0 * lk_1039[k];

        t_1040[k] = f_0 * lk_1040[k];

        t_1041[k] = f_0 * lk_1041[k];

        t_1042[k] = f_0 * lk_1042[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, t_1046, t_1047, t_1048, t_1049, t_1050, \
                         lk_1043, lk_1044, lk_1045, lk_1046, lk_1047, lk_1048, lk_1049, \
                         lk_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_0 * lk_1043[k];

        t_1044[k] = f_0 * lk_1044[k];

        t_1045[k] = f_0 * lk_1045[k];

        t_1046[k] = f_0 * lk_1046[k];

        t_1047[k] = f_0 * lk_1047[k];

        t_1048[k] = f_0 * lk_1048[k];

        t_1049[k] = f_0 * lk_1049[k];

        t_1050[k] = f_0 * lk_1050[k];
    }

#pragma omp simd aligned(t_1051, t_1052, t_1053, t_1054, t_1055, t_1056, t_1057, t_1058, \
                         lk_1051, lk_1052, lk_1053, lk_1054, lk_1055, lk_1056, lk_1057, \
                         lk_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1051[k] = f_0 * lk_1051[k];

        t_1052[k] = f_0 * lk_1052[k];

        t_1053[k] = f_0 * lk_1053[k];

        t_1054[k] = f_0 * lk_1054[k];

        t_1055[k] = f_0 * lk_1055[k];

        t_1056[k] = f_0 * lk_1056[k];

        t_1057[k] = f_0 * lk_1057[k];

        t_1058[k] = f_0 * lk_1058[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, t_1062, t_1063, t_1064, t_1065, t_1066, \
                         lk_1059, lk_1060, lk_1061, lk_1062, lk_1063, lk_1064, lk_1065, \
                         lk_1066 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = f_0 * lk_1059[k];

        t_1060[k] = f_0 * lk_1060[k];

        t_1061[k] = f_0 * lk_1061[k];

        t_1062[k] = f_0 * lk_1062[k];

        t_1063[k] = f_0 * lk_1063[k];

        t_1064[k] = f_0 * lk_1064[k];

        t_1065[k] = f_0 * lk_1065[k];

        t_1066[k] = f_0 * lk_1066[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_0_piece7(CSimdMatrix &buffer, const size_t target,
                                                    const size_t lk, const size_t ncols,
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
    auto *t_1260 = buffer.data(target + 1260);
    auto *t_1261 = buffer.data(target + 1261);
    auto *t_1262 = buffer.data(target + 1262);
    auto *t_1263 = buffer.data(target + 1263);
    auto *t_1264 = buffer.data(target + 1264);
    auto *t_1265 = buffer.data(target + 1265);
    auto *t_1266 = buffer.data(target + 1266);
    auto *t_1267 = buffer.data(target + 1267);
    auto *t_1268 = buffer.data(target + 1268);
    auto *t_1269 = buffer.data(target + 1269);
    auto *t_1270 = buffer.data(target + 1270);
    auto *t_1271 = buffer.data(target + 1271);
    auto *t_1272 = buffer.data(target + 1272);
    auto *t_1273 = buffer.data(target + 1273);
    auto *t_1274 = buffer.data(target + 1274);
    auto *t_1275 = buffer.data(target + 1275);
    auto *t_1276 = buffer.data(target + 1276);
    auto *t_1277 = buffer.data(target + 1277);
    auto *t_1278 = buffer.data(target + 1278);
    auto *t_1279 = buffer.data(target + 1279);
    auto *t_1280 = buffer.data(target + 1280);
    auto *t_1281 = buffer.data(target + 1281);
    auto *t_1282 = buffer.data(target + 1282);
    auto *t_1283 = buffer.data(target + 1283);
    auto *t_1284 = buffer.data(target + 1284);
    auto *t_1285 = buffer.data(target + 1285);
    auto *t_1286 = buffer.data(target + 1286);
    auto *t_1287 = buffer.data(target + 1287);
    auto *t_1288 = buffer.data(target + 1288);
    auto *t_1289 = buffer.data(target + 1289);
    auto *t_1290 = buffer.data(target + 1290);

    const auto *lk_1067 = buffer.data(lk + 1067);
    const auto *lk_1068 = buffer.data(lk + 1068);
    const auto *lk_1069 = buffer.data(lk + 1069);
    const auto *lk_1070 = buffer.data(lk + 1070);
    const auto *lk_1071 = buffer.data(lk + 1071);
    const auto *lk_1072 = buffer.data(lk + 1072);
    const auto *lk_1073 = buffer.data(lk + 1073);
    const auto *lk_1074 = buffer.data(lk + 1074);
    const auto *lk_1075 = buffer.data(lk + 1075);
    const auto *lk_1076 = buffer.data(lk + 1076);
    const auto *lk_1077 = buffer.data(lk + 1077);
    const auto *lk_1078 = buffer.data(lk + 1078);
    const auto *lk_1079 = buffer.data(lk + 1079);
    const auto *lk_1080 = buffer.data(lk + 1080);
    const auto *lk_1081 = buffer.data(lk + 1081);
    const auto *lk_1082 = buffer.data(lk + 1082);
    const auto *lk_1083 = buffer.data(lk + 1083);
    const auto *lk_1084 = buffer.data(lk + 1084);
    const auto *lk_1085 = buffer.data(lk + 1085);
    const auto *lk_1086 = buffer.data(lk + 1086);
    const auto *lk_1087 = buffer.data(lk + 1087);
    const auto *lk_1088 = buffer.data(lk + 1088);
    const auto *lk_1089 = buffer.data(lk + 1089);
    const auto *lk_1090 = buffer.data(lk + 1090);
    const auto *lk_1091 = buffer.data(lk + 1091);
    const auto *lk_1092 = buffer.data(lk + 1092);
    const auto *lk_1093 = buffer.data(lk + 1093);
    const auto *lk_1094 = buffer.data(lk + 1094);
    const auto *lk_1095 = buffer.data(lk + 1095);
    const auto *lk_1096 = buffer.data(lk + 1096);
    const auto *lk_1097 = buffer.data(lk + 1097);
    const auto *lk_1098 = buffer.data(lk + 1098);
    const auto *lk_1099 = buffer.data(lk + 1099);
    const auto *lk_1100 = buffer.data(lk + 1100);
    const auto *lk_1101 = buffer.data(lk + 1101);
    const auto *lk_1102 = buffer.data(lk + 1102);
    const auto *lk_1103 = buffer.data(lk + 1103);
    const auto *lk_1104 = buffer.data(lk + 1104);
    const auto *lk_1105 = buffer.data(lk + 1105);
    const auto *lk_1106 = buffer.data(lk + 1106);
    const auto *lk_1107 = buffer.data(lk + 1107);
    const auto *lk_1108 = buffer.data(lk + 1108);
    const auto *lk_1109 = buffer.data(lk + 1109);
    const auto *lk_1110 = buffer.data(lk + 1110);
    const auto *lk_1111 = buffer.data(lk + 1111);
    const auto *lk_1112 = buffer.data(lk + 1112);
    const auto *lk_1113 = buffer.data(lk + 1113);
    const auto *lk_1114 = buffer.data(lk + 1114);
    const auto *lk_1115 = buffer.data(lk + 1115);
    const auto *lk_1116 = buffer.data(lk + 1116);
    const auto *lk_1117 = buffer.data(lk + 1117);
    const auto *lk_1118 = buffer.data(lk + 1118);
    const auto *lk_1119 = buffer.data(lk + 1119);
    const auto *lk_1120 = buffer.data(lk + 1120);
    const auto *lk_1121 = buffer.data(lk + 1121);
    const auto *lk_1122 = buffer.data(lk + 1122);
    const auto *lk_1123 = buffer.data(lk + 1123);
    const auto *lk_1124 = buffer.data(lk + 1124);
    const auto *lk_1125 = buffer.data(lk + 1125);
    const auto *lk_1126 = buffer.data(lk + 1126);
    const auto *lk_1127 = buffer.data(lk + 1127);
    const auto *lk_1128 = buffer.data(lk + 1128);
    const auto *lk_1129 = buffer.data(lk + 1129);
    const auto *lk_1130 = buffer.data(lk + 1130);
    const auto *lk_1131 = buffer.data(lk + 1131);
    const auto *lk_1132 = buffer.data(lk + 1132);
    const auto *lk_1133 = buffer.data(lk + 1133);
    const auto *lk_1134 = buffer.data(lk + 1134);
    const auto *lk_1135 = buffer.data(lk + 1135);
    const auto *lk_1136 = buffer.data(lk + 1136);
    const auto *lk_1137 = buffer.data(lk + 1137);
    const auto *lk_1138 = buffer.data(lk + 1138);
    const auto *lk_1139 = buffer.data(lk + 1139);
    const auto *lk_1140 = buffer.data(lk + 1140);
    const auto *lk_1141 = buffer.data(lk + 1141);
    const auto *lk_1142 = buffer.data(lk + 1142);
    const auto *lk_1143 = buffer.data(lk + 1143);
    const auto *lk_1144 = buffer.data(lk + 1144);
    const auto *lk_1145 = buffer.data(lk + 1145);
    const auto *lk_1146 = buffer.data(lk + 1146);
    const auto *lk_1147 = buffer.data(lk + 1147);
    const auto *lk_1148 = buffer.data(lk + 1148);
    const auto *lk_1149 = buffer.data(lk + 1149);
    const auto *lk_1150 = buffer.data(lk + 1150);
    const auto *lk_1151 = buffer.data(lk + 1151);
    const auto *lk_1152 = buffer.data(lk + 1152);
    const auto *lk_1153 = buffer.data(lk + 1153);
    const auto *lk_1154 = buffer.data(lk + 1154);
    const auto *lk_1155 = buffer.data(lk + 1155);
    const auto *lk_1156 = buffer.data(lk + 1156);
    const auto *lk_1157 = buffer.data(lk + 1157);
    const auto *lk_1158 = buffer.data(lk + 1158);
    const auto *lk_1159 = buffer.data(lk + 1159);
    const auto *lk_1160 = buffer.data(lk + 1160);
    const auto *lk_1161 = buffer.data(lk + 1161);
    const auto *lk_1162 = buffer.data(lk + 1162);
    const auto *lk_1163 = buffer.data(lk + 1163);
    const auto *lk_1164 = buffer.data(lk + 1164);
    const auto *lk_1165 = buffer.data(lk + 1165);
    const auto *lk_1166 = buffer.data(lk + 1166);
    const auto *lk_1167 = buffer.data(lk + 1167);
    const auto *lk_1168 = buffer.data(lk + 1168);
    const auto *lk_1169 = buffer.data(lk + 1169);
    const auto *lk_1170 = buffer.data(lk + 1170);
    const auto *lk_1171 = buffer.data(lk + 1171);
    const auto *lk_1172 = buffer.data(lk + 1172);
    const auto *lk_1173 = buffer.data(lk + 1173);
    const auto *lk_1174 = buffer.data(lk + 1174);
    const auto *lk_1175 = buffer.data(lk + 1175);
    const auto *lk_1176 = buffer.data(lk + 1176);
    const auto *lk_1177 = buffer.data(lk + 1177);
    const auto *lk_1178 = buffer.data(lk + 1178);
    const auto *lk_1179 = buffer.data(lk + 1179);
    const auto *lk_1180 = buffer.data(lk + 1180);
    const auto *lk_1181 = buffer.data(lk + 1181);
    const auto *lk_1182 = buffer.data(lk + 1182);
    const auto *lk_1183 = buffer.data(lk + 1183);
    const auto *lk_1184 = buffer.data(lk + 1184);
    const auto *lk_1185 = buffer.data(lk + 1185);
    const auto *lk_1186 = buffer.data(lk + 1186);
    const auto *lk_1187 = buffer.data(lk + 1187);
    const auto *lk_1188 = buffer.data(lk + 1188);
    const auto *lk_1189 = buffer.data(lk + 1189);
    const auto *lk_1190 = buffer.data(lk + 1190);
    const auto *lk_1191 = buffer.data(lk + 1191);
    const auto *lk_1192 = buffer.data(lk + 1192);
    const auto *lk_1193 = buffer.data(lk + 1193);
    const auto *lk_1194 = buffer.data(lk + 1194);
    const auto *lk_1195 = buffer.data(lk + 1195);
    const auto *lk_1196 = buffer.data(lk + 1196);
    const auto *lk_1197 = buffer.data(lk + 1197);
    const auto *lk_1198 = buffer.data(lk + 1198);
    const auto *lk_1199 = buffer.data(lk + 1199);
    const auto *lk_1200 = buffer.data(lk + 1200);
    const auto *lk_1201 = buffer.data(lk + 1201);
    const auto *lk_1202 = buffer.data(lk + 1202);
    const auto *lk_1203 = buffer.data(lk + 1203);
    const auto *lk_1204 = buffer.data(lk + 1204);
    const auto *lk_1205 = buffer.data(lk + 1205);
    const auto *lk_1206 = buffer.data(lk + 1206);
    const auto *lk_1207 = buffer.data(lk + 1207);
    const auto *lk_1208 = buffer.data(lk + 1208);
    const auto *lk_1209 = buffer.data(lk + 1209);
    const auto *lk_1210 = buffer.data(lk + 1210);
    const auto *lk_1211 = buffer.data(lk + 1211);
    const auto *lk_1212 = buffer.data(lk + 1212);
    const auto *lk_1213 = buffer.data(lk + 1213);
    const auto *lk_1214 = buffer.data(lk + 1214);
    const auto *lk_1215 = buffer.data(lk + 1215);
    const auto *lk_1216 = buffer.data(lk + 1216);
    const auto *lk_1217 = buffer.data(lk + 1217);
    const auto *lk_1218 = buffer.data(lk + 1218);
    const auto *lk_1219 = buffer.data(lk + 1219);
    const auto *lk_1220 = buffer.data(lk + 1220);
    const auto *lk_1221 = buffer.data(lk + 1221);
    const auto *lk_1222 = buffer.data(lk + 1222);
    const auto *lk_1223 = buffer.data(lk + 1223);
    const auto *lk_1224 = buffer.data(lk + 1224);
    const auto *lk_1225 = buffer.data(lk + 1225);
    const auto *lk_1226 = buffer.data(lk + 1226);
    const auto *lk_1227 = buffer.data(lk + 1227);
    const auto *lk_1228 = buffer.data(lk + 1228);
    const auto *lk_1229 = buffer.data(lk + 1229);
    const auto *lk_1230 = buffer.data(lk + 1230);
    const auto *lk_1231 = buffer.data(lk + 1231);
    const auto *lk_1232 = buffer.data(lk + 1232);
    const auto *lk_1233 = buffer.data(lk + 1233);
    const auto *lk_1234 = buffer.data(lk + 1234);
    const auto *lk_1235 = buffer.data(lk + 1235);
    const auto *lk_1236 = buffer.data(lk + 1236);
    const auto *lk_1237 = buffer.data(lk + 1237);
    const auto *lk_1238 = buffer.data(lk + 1238);
    const auto *lk_1239 = buffer.data(lk + 1239);
    const auto *lk_1240 = buffer.data(lk + 1240);
    const auto *lk_1241 = buffer.data(lk + 1241);
    const auto *lk_1242 = buffer.data(lk + 1242);
    const auto *lk_1243 = buffer.data(lk + 1243);
    const auto *lk_1244 = buffer.data(lk + 1244);
    const auto *lk_1245 = buffer.data(lk + 1245);
    const auto *lk_1246 = buffer.data(lk + 1246);
    const auto *lk_1247 = buffer.data(lk + 1247);
    const auto *lk_1248 = buffer.data(lk + 1248);
    const auto *lk_1249 = buffer.data(lk + 1249);
    const auto *lk_1250 = buffer.data(lk + 1250);
    const auto *lk_1251 = buffer.data(lk + 1251);
    const auto *lk_1252 = buffer.data(lk + 1252);
    const auto *lk_1253 = buffer.data(lk + 1253);
    const auto *lk_1254 = buffer.data(lk + 1254);
    const auto *lk_1255 = buffer.data(lk + 1255);
    const auto *lk_1256 = buffer.data(lk + 1256);
    const auto *lk_1257 = buffer.data(lk + 1257);
    const auto *lk_1258 = buffer.data(lk + 1258);
    const auto *lk_1259 = buffer.data(lk + 1259);
    const auto *lk_1260 = buffer.data(lk + 1260);
    const auto *lk_1261 = buffer.data(lk + 1261);
    const auto *lk_1262 = buffer.data(lk + 1262);
    const auto *lk_1263 = buffer.data(lk + 1263);
    const auto *lk_1264 = buffer.data(lk + 1264);
    const auto *lk_1265 = buffer.data(lk + 1265);
    const auto *lk_1266 = buffer.data(lk + 1266);
    const auto *lk_1267 = buffer.data(lk + 1267);
    const auto *lk_1268 = buffer.data(lk + 1268);
    const auto *lk_1269 = buffer.data(lk + 1269);
    const auto *lk_1270 = buffer.data(lk + 1270);
    const auto *lk_1271 = buffer.data(lk + 1271);
    const auto *lk_1272 = buffer.data(lk + 1272);
    const auto *lk_1273 = buffer.data(lk + 1273);
    const auto *lk_1274 = buffer.data(lk + 1274);
    const auto *lk_1275 = buffer.data(lk + 1275);
    const auto *lk_1276 = buffer.data(lk + 1276);
    const auto *lk_1277 = buffer.data(lk + 1277);
    const auto *lk_1278 = buffer.data(lk + 1278);
    const auto *lk_1279 = buffer.data(lk + 1279);
    const auto *lk_1280 = buffer.data(lk + 1280);
    const auto *lk_1281 = buffer.data(lk + 1281);
    const auto *lk_1282 = buffer.data(lk + 1282);
    const auto *lk_1283 = buffer.data(lk + 1283);
    const auto *lk_1284 = buffer.data(lk + 1284);
    const auto *lk_1285 = buffer.data(lk + 1285);
    const auto *lk_1286 = buffer.data(lk + 1286);
    const auto *lk_1287 = buffer.data(lk + 1287);
    const auto *lk_1288 = buffer.data(lk + 1288);
    const auto *lk_1289 = buffer.data(lk + 1289);
    const auto *lk_1290 = buffer.data(lk + 1290);

#pragma omp simd aligned(t_1067, t_1068, t_1069, t_1070, t_1071, t_1072, t_1073, t_1074, \
                         lk_1067, lk_1068, lk_1069, lk_1070, lk_1071, lk_1072, lk_1073, \
                         lk_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_0 * lk_1067[k];

        t_1068[k] = f_0 * lk_1068[k];

        t_1069[k] = f_0 * lk_1069[k];

        t_1070[k] = f_0 * lk_1070[k];

        t_1071[k] = f_0 * lk_1071[k];

        t_1072[k] = f_0 * lk_1072[k];

        t_1073[k] = f_0 * lk_1073[k];

        t_1074[k] = f_0 * lk_1074[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, t_1078, t_1079, t_1080, t_1081, t_1082, \
                         lk_1075, lk_1076, lk_1077, lk_1078, lk_1079, lk_1080, lk_1081, \
                         lk_1082 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = f_0 * lk_1075[k];

        t_1076[k] = f_0 * lk_1076[k];

        t_1077[k] = f_0 * lk_1077[k];

        t_1078[k] = f_0 * lk_1078[k];

        t_1079[k] = f_0 * lk_1079[k];

        t_1080[k] = f_0 * lk_1080[k];

        t_1081[k] = f_0 * lk_1081[k];

        t_1082[k] = f_0 * lk_1082[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, t_1086, t_1087, t_1088, t_1089, t_1090, \
                         lk_1083, lk_1084, lk_1085, lk_1086, lk_1087, lk_1088, lk_1089, \
                         lk_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = f_0 * lk_1083[k];

        t_1084[k] = f_0 * lk_1084[k];

        t_1085[k] = f_0 * lk_1085[k];

        t_1086[k] = f_0 * lk_1086[k];

        t_1087[k] = f_0 * lk_1087[k];

        t_1088[k] = f_0 * lk_1088[k];

        t_1089[k] = f_0 * lk_1089[k];

        t_1090[k] = f_0 * lk_1090[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, t_1094, t_1095, t_1096, t_1097, t_1098, \
                         lk_1091, lk_1092, lk_1093, lk_1094, lk_1095, lk_1096, lk_1097, \
                         lk_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_0 * lk_1091[k];

        t_1092[k] = f_0 * lk_1092[k];

        t_1093[k] = f_0 * lk_1093[k];

        t_1094[k] = f_0 * lk_1094[k];

        t_1095[k] = f_0 * lk_1095[k];

        t_1096[k] = f_0 * lk_1096[k];

        t_1097[k] = f_0 * lk_1097[k];

        t_1098[k] = f_0 * lk_1098[k];
    }

#pragma omp simd aligned(t_1099, t_1100, t_1101, t_1102, t_1103, t_1104, t_1105, t_1106, \
                         lk_1099, lk_1100, lk_1101, lk_1102, lk_1103, lk_1104, lk_1105, \
                         lk_1106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1099[k] = f_0 * lk_1099[k];

        t_1100[k] = f_0 * lk_1100[k];

        t_1101[k] = f_0 * lk_1101[k];

        t_1102[k] = f_0 * lk_1102[k];

        t_1103[k] = f_0 * lk_1103[k];

        t_1104[k] = f_0 * lk_1104[k];

        t_1105[k] = f_0 * lk_1105[k];

        t_1106[k] = f_0 * lk_1106[k];
    }

#pragma omp simd aligned(t_1107, t_1108, t_1109, t_1110, t_1111, t_1112, t_1113, t_1114, \
                         lk_1107, lk_1108, lk_1109, lk_1110, lk_1111, lk_1112, lk_1113, \
                         lk_1114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1107[k] = f_0 * lk_1107[k];

        t_1108[k] = f_0 * lk_1108[k];

        t_1109[k] = f_0 * lk_1109[k];

        t_1110[k] = f_0 * lk_1110[k];

        t_1111[k] = f_0 * lk_1111[k];

        t_1112[k] = f_0 * lk_1112[k];

        t_1113[k] = f_0 * lk_1113[k];

        t_1114[k] = f_0 * lk_1114[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, t_1118, t_1119, t_1120, t_1121, t_1122, \
                         lk_1115, lk_1116, lk_1117, lk_1118, lk_1119, lk_1120, lk_1121, \
                         lk_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_0 * lk_1115[k];

        t_1116[k] = f_0 * lk_1116[k];

        t_1117[k] = f_0 * lk_1117[k];

        t_1118[k] = f_0 * lk_1118[k];

        t_1119[k] = f_0 * lk_1119[k];

        t_1120[k] = f_0 * lk_1120[k];

        t_1121[k] = f_0 * lk_1121[k];

        t_1122[k] = f_0 * lk_1122[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, t_1126, t_1127, t_1128, t_1129, t_1130, \
                         lk_1123, lk_1124, lk_1125, lk_1126, lk_1127, lk_1128, lk_1129, \
                         lk_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_0 * lk_1123[k];

        t_1124[k] = f_0 * lk_1124[k];

        t_1125[k] = f_0 * lk_1125[k];

        t_1126[k] = f_0 * lk_1126[k];

        t_1127[k] = f_0 * lk_1127[k];

        t_1128[k] = f_0 * lk_1128[k];

        t_1129[k] = f_0 * lk_1129[k];

        t_1130[k] = f_0 * lk_1130[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, t_1134, t_1135, t_1136, t_1137, t_1138, \
                         lk_1131, lk_1132, lk_1133, lk_1134, lk_1135, lk_1136, lk_1137, \
                         lk_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = f_0 * lk_1131[k];

        t_1132[k] = f_0 * lk_1132[k];

        t_1133[k] = f_0 * lk_1133[k];

        t_1134[k] = f_0 * lk_1134[k];

        t_1135[k] = f_0 * lk_1135[k];

        t_1136[k] = f_0 * lk_1136[k];

        t_1137[k] = f_0 * lk_1137[k];

        t_1138[k] = f_0 * lk_1138[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, t_1142, t_1143, t_1144, t_1145, t_1146, \
                         lk_1139, lk_1140, lk_1141, lk_1142, lk_1143, lk_1144, lk_1145, \
                         lk_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_0 * lk_1139[k];

        t_1140[k] = f_0 * lk_1140[k];

        t_1141[k] = f_0 * lk_1141[k];

        t_1142[k] = f_0 * lk_1142[k];

        t_1143[k] = f_0 * lk_1143[k];

        t_1144[k] = f_0 * lk_1144[k];

        t_1145[k] = f_0 * lk_1145[k];

        t_1146[k] = f_0 * lk_1146[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, t_1150, t_1151, t_1152, t_1153, t_1154, \
                         lk_1147, lk_1148, lk_1149, lk_1150, lk_1151, lk_1152, lk_1153, \
                         lk_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_0 * lk_1147[k];

        t_1148[k] = f_0 * lk_1148[k];

        t_1149[k] = f_0 * lk_1149[k];

        t_1150[k] = f_0 * lk_1150[k];

        t_1151[k] = f_0 * lk_1151[k];

        t_1152[k] = f_0 * lk_1152[k];

        t_1153[k] = f_0 * lk_1153[k];

        t_1154[k] = f_0 * lk_1154[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, t_1159, t_1160, t_1161, t_1162, \
                         lk_1155, lk_1156, lk_1157, lk_1158, lk_1159, lk_1160, lk_1161, \
                         lk_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_0 * lk_1155[k];

        t_1156[k] = f_0 * lk_1156[k];

        t_1157[k] = f_0 * lk_1157[k];

        t_1158[k] = f_0 * lk_1158[k];

        t_1159[k] = f_0 * lk_1159[k];

        t_1160[k] = f_0 * lk_1160[k];

        t_1161[k] = f_0 * lk_1161[k];

        t_1162[k] = f_0 * lk_1162[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, t_1166, t_1167, t_1168, t_1169, t_1170, \
                         lk_1163, lk_1164, lk_1165, lk_1166, lk_1167, lk_1168, lk_1169, \
                         lk_1170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_0 * lk_1163[k];

        t_1164[k] = f_0 * lk_1164[k];

        t_1165[k] = f_0 * lk_1165[k];

        t_1166[k] = f_0 * lk_1166[k];

        t_1167[k] = f_0 * lk_1167[k];

        t_1168[k] = f_0 * lk_1168[k];

        t_1169[k] = f_0 * lk_1169[k];

        t_1170[k] = f_0 * lk_1170[k];
    }

#pragma omp simd aligned(t_1171, t_1172, t_1173, t_1174, t_1175, t_1176, t_1177, t_1178, \
                         lk_1171, lk_1172, lk_1173, lk_1174, lk_1175, lk_1176, lk_1177, \
                         lk_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1171[k] = f_0 * lk_1171[k];

        t_1172[k] = f_0 * lk_1172[k];

        t_1173[k] = f_0 * lk_1173[k];

        t_1174[k] = f_0 * lk_1174[k];

        t_1175[k] = f_0 * lk_1175[k];

        t_1176[k] = f_0 * lk_1176[k];

        t_1177[k] = f_0 * lk_1177[k];

        t_1178[k] = f_0 * lk_1178[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, t_1182, t_1183, t_1184, t_1185, t_1186, \
                         lk_1179, lk_1180, lk_1181, lk_1182, lk_1183, lk_1184, lk_1185, \
                         lk_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = f_0 * lk_1179[k];

        t_1180[k] = f_0 * lk_1180[k];

        t_1181[k] = f_0 * lk_1181[k];

        t_1182[k] = f_0 * lk_1182[k];

        t_1183[k] = f_0 * lk_1183[k];

        t_1184[k] = f_0 * lk_1184[k];

        t_1185[k] = f_0 * lk_1185[k];

        t_1186[k] = f_0 * lk_1186[k];
    }

#pragma omp simd aligned(t_1187, t_1188, t_1189, t_1190, t_1191, t_1192, t_1193, t_1194, \
                         lk_1187, lk_1188, lk_1189, lk_1190, lk_1191, lk_1192, lk_1193, \
                         lk_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1187[k] = f_0 * lk_1187[k];

        t_1188[k] = f_0 * lk_1188[k];

        t_1189[k] = f_0 * lk_1189[k];

        t_1190[k] = f_0 * lk_1190[k];

        t_1191[k] = f_0 * lk_1191[k];

        t_1192[k] = f_0 * lk_1192[k];

        t_1193[k] = f_0 * lk_1193[k];

        t_1194[k] = f_0 * lk_1194[k];
    }

#pragma omp simd aligned(t_1195, t_1196, t_1197, t_1198, t_1199, t_1200, t_1201, t_1202, \
                         lk_1195, lk_1196, lk_1197, lk_1198, lk_1199, lk_1200, lk_1201, \
                         lk_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1195[k] = f_0 * lk_1195[k];

        t_1196[k] = f_0 * lk_1196[k];

        t_1197[k] = f_0 * lk_1197[k];

        t_1198[k] = f_0 * lk_1198[k];

        t_1199[k] = f_0 * lk_1199[k];

        t_1200[k] = f_0 * lk_1200[k];

        t_1201[k] = f_0 * lk_1201[k];

        t_1202[k] = f_0 * lk_1202[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, t_1206, t_1207, t_1208, t_1209, t_1210, \
                         lk_1203, lk_1204, lk_1205, lk_1206, lk_1207, lk_1208, lk_1209, \
                         lk_1210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = f_0 * lk_1203[k];

        t_1204[k] = f_0 * lk_1204[k];

        t_1205[k] = f_0 * lk_1205[k];

        t_1206[k] = f_0 * lk_1206[k];

        t_1207[k] = f_0 * lk_1207[k];

        t_1208[k] = f_0 * lk_1208[k];

        t_1209[k] = f_0 * lk_1209[k];

        t_1210[k] = f_0 * lk_1210[k];
    }

#pragma omp simd aligned(t_1211, t_1212, t_1213, t_1214, t_1215, t_1216, t_1217, t_1218, \
                         lk_1211, lk_1212, lk_1213, lk_1214, lk_1215, lk_1216, lk_1217, \
                         lk_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1211[k] = f_0 * lk_1211[k];

        t_1212[k] = f_0 * lk_1212[k];

        t_1213[k] = f_0 * lk_1213[k];

        t_1214[k] = f_0 * lk_1214[k];

        t_1215[k] = f_0 * lk_1215[k];

        t_1216[k] = f_0 * lk_1216[k];

        t_1217[k] = f_0 * lk_1217[k];

        t_1218[k] = f_0 * lk_1218[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, t_1222, t_1223, t_1224, t_1225, t_1226, \
                         lk_1219, lk_1220, lk_1221, lk_1222, lk_1223, lk_1224, lk_1225, \
                         lk_1226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_0 * lk_1219[k];

        t_1220[k] = f_0 * lk_1220[k];

        t_1221[k] = f_0 * lk_1221[k];

        t_1222[k] = f_0 * lk_1222[k];

        t_1223[k] = f_0 * lk_1223[k];

        t_1224[k] = f_0 * lk_1224[k];

        t_1225[k] = f_0 * lk_1225[k];

        t_1226[k] = f_0 * lk_1226[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, t_1230, t_1231, t_1232, t_1233, t_1234, \
                         lk_1227, lk_1228, lk_1229, lk_1230, lk_1231, lk_1232, lk_1233, \
                         lk_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = f_0 * lk_1227[k];

        t_1228[k] = f_0 * lk_1228[k];

        t_1229[k] = f_0 * lk_1229[k];

        t_1230[k] = f_0 * lk_1230[k];

        t_1231[k] = f_0 * lk_1231[k];

        t_1232[k] = f_0 * lk_1232[k];

        t_1233[k] = f_0 * lk_1233[k];

        t_1234[k] = f_0 * lk_1234[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, t_1238, t_1239, t_1240, t_1241, t_1242, \
                         lk_1235, lk_1236, lk_1237, lk_1238, lk_1239, lk_1240, lk_1241, \
                         lk_1242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_0 * lk_1235[k];

        t_1236[k] = f_0 * lk_1236[k];

        t_1237[k] = f_0 * lk_1237[k];

        t_1238[k] = f_0 * lk_1238[k];

        t_1239[k] = f_0 * lk_1239[k];

        t_1240[k] = f_0 * lk_1240[k];

        t_1241[k] = f_0 * lk_1241[k];

        t_1242[k] = f_0 * lk_1242[k];
    }

#pragma omp simd aligned(t_1243, t_1244, t_1245, t_1246, t_1247, t_1248, t_1249, t_1250, \
                         lk_1243, lk_1244, lk_1245, lk_1246, lk_1247, lk_1248, lk_1249, \
                         lk_1250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1243[k] = f_0 * lk_1243[k];

        t_1244[k] = f_0 * lk_1244[k];

        t_1245[k] = f_0 * lk_1245[k];

        t_1246[k] = f_0 * lk_1246[k];

        t_1247[k] = f_0 * lk_1247[k];

        t_1248[k] = f_0 * lk_1248[k];

        t_1249[k] = f_0 * lk_1249[k];

        t_1250[k] = f_0 * lk_1250[k];
    }

#pragma omp simd aligned(t_1251, t_1252, t_1253, t_1254, t_1255, t_1256, t_1257, t_1258, \
                         lk_1251, lk_1252, lk_1253, lk_1254, lk_1255, lk_1256, lk_1257, \
                         lk_1258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1251[k] = f_0 * lk_1251[k];

        t_1252[k] = f_0 * lk_1252[k];

        t_1253[k] = f_0 * lk_1253[k];

        t_1254[k] = f_0 * lk_1254[k];

        t_1255[k] = f_0 * lk_1255[k];

        t_1256[k] = f_0 * lk_1256[k];

        t_1257[k] = f_0 * lk_1257[k];

        t_1258[k] = f_0 * lk_1258[k];
    }

#pragma omp simd aligned(t_1259, t_1260, t_1261, t_1262, t_1263, t_1264, t_1265, t_1266, \
                         lk_1259, lk_1260, lk_1261, lk_1262, lk_1263, lk_1264, lk_1265, \
                         lk_1266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1259[k] = f_0 * lk_1259[k];

        t_1260[k] = f_0 * lk_1260[k];

        t_1261[k] = f_0 * lk_1261[k];

        t_1262[k] = f_0 * lk_1262[k];

        t_1263[k] = f_0 * lk_1263[k];

        t_1264[k] = f_0 * lk_1264[k];

        t_1265[k] = f_0 * lk_1265[k];

        t_1266[k] = f_0 * lk_1266[k];
    }

#pragma omp simd aligned(t_1267, t_1268, t_1269, t_1270, t_1271, t_1272, t_1273, t_1274, \
                         lk_1267, lk_1268, lk_1269, lk_1270, lk_1271, lk_1272, lk_1273, \
                         lk_1274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1267[k] = f_0 * lk_1267[k];

        t_1268[k] = f_0 * lk_1268[k];

        t_1269[k] = f_0 * lk_1269[k];

        t_1270[k] = f_0 * lk_1270[k];

        t_1271[k] = f_0 * lk_1271[k];

        t_1272[k] = f_0 * lk_1272[k];

        t_1273[k] = f_0 * lk_1273[k];

        t_1274[k] = f_0 * lk_1274[k];
    }

#pragma omp simd aligned(t_1275, t_1276, t_1277, t_1278, t_1279, t_1280, t_1281, t_1282, \
                         lk_1275, lk_1276, lk_1277, lk_1278, lk_1279, lk_1280, lk_1281, \
                         lk_1282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1275[k] = f_0 * lk_1275[k];

        t_1276[k] = f_0 * lk_1276[k];

        t_1277[k] = f_0 * lk_1277[k];

        t_1278[k] = f_0 * lk_1278[k];

        t_1279[k] = f_0 * lk_1279[k];

        t_1280[k] = f_0 * lk_1280[k];

        t_1281[k] = f_0 * lk_1281[k];

        t_1282[k] = f_0 * lk_1282[k];
    }

#pragma omp simd aligned(t_1283, t_1284, t_1285, t_1286, t_1287, t_1288, t_1289, t_1290, \
                         lk_1283, lk_1284, lk_1285, lk_1286, lk_1287, lk_1288, lk_1289, \
                         lk_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1283[k] = f_0 * lk_1283[k];

        t_1284[k] = f_0 * lk_1284[k];

        t_1285[k] = f_0 * lk_1285[k];

        t_1286[k] = f_0 * lk_1286[k];

        t_1287[k] = f_0 * lk_1287[k];

        t_1288[k] = f_0 * lk_1288[k];

        t_1289[k] = f_0 * lk_1289[k];

        t_1290[k] = f_0 * lk_1290[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_0_piece8(CSimdMatrix &buffer, const size_t target,
                                                    const size_t lk, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

    auto *t_1291 = buffer.data(target + 1291);
    auto *t_1292 = buffer.data(target + 1292);
    auto *t_1293 = buffer.data(target + 1293);
    auto *t_1294 = buffer.data(target + 1294);
    auto *t_1295 = buffer.data(target + 1295);

    const auto *lk_1291 = buffer.data(lk + 1291);
    const auto *lk_1292 = buffer.data(lk + 1292);
    const auto *lk_1293 = buffer.data(lk + 1293);
    const auto *lk_1294 = buffer.data(lk + 1294);
    const auto *lk_1295 = buffer.data(lk + 1295);

#pragma omp simd aligned(t_1291, t_1292, t_1293, t_1294, t_1295, lk_1291, lk_1292, lk_1293, \
                         lk_1294, lk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1291[k] = f_0 * lk_1291[k];

        t_1292[k] = f_0 * lk_1292[k];

        t_1293[k] = f_0 * lk_1293[k];

        t_1294[k] = f_0 * lk_1294[k];

        t_1295[k] = f_0 * lk_1295[k];
    }
}

auto
compute_prim_geom_10_kk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t ik, const size_t lk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kk_electron_repulsion_0_piece0(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_0_piece1(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_0_piece2(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_0_piece3(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_0_piece4(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_0_piece5(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_0_piece6(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_0_piece7(buffer, target, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_0_piece8(buffer, target, lk, ncols, alpha);
}

static auto
compute_prim_geom_10_kk_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_0 = buffer.data(ik + 0);
    const auto *ik_1 = buffer.data(ik + 1);
    const auto *ik_2 = buffer.data(ik + 2);
    const auto *ik_3 = buffer.data(ik + 3);
    const auto *ik_4 = buffer.data(ik + 4);
    const auto *ik_5 = buffer.data(ik + 5);
    const auto *ik_6 = buffer.data(ik + 6);
    const auto *ik_7 = buffer.data(ik + 7);
    const auto *ik_8 = buffer.data(ik + 8);
    const auto *ik_9 = buffer.data(ik + 9);
    const auto *ik_10 = buffer.data(ik + 10);
    const auto *ik_11 = buffer.data(ik + 11);
    const auto *ik_12 = buffer.data(ik + 12);
    const auto *ik_13 = buffer.data(ik + 13);
    const auto *ik_14 = buffer.data(ik + 14);
    const auto *ik_15 = buffer.data(ik + 15);
    const auto *ik_16 = buffer.data(ik + 16);
    const auto *ik_17 = buffer.data(ik + 17);
    const auto *ik_18 = buffer.data(ik + 18);
    const auto *ik_19 = buffer.data(ik + 19);
    const auto *ik_20 = buffer.data(ik + 20);
    const auto *ik_21 = buffer.data(ik + 21);
    const auto *ik_22 = buffer.data(ik + 22);
    const auto *ik_23 = buffer.data(ik + 23);
    const auto *ik_24 = buffer.data(ik + 24);
    const auto *ik_25 = buffer.data(ik + 25);
    const auto *ik_26 = buffer.data(ik + 26);
    const auto *ik_27 = buffer.data(ik + 27);
    const auto *ik_28 = buffer.data(ik + 28);
    const auto *ik_29 = buffer.data(ik + 29);
    const auto *ik_30 = buffer.data(ik + 30);
    const auto *ik_31 = buffer.data(ik + 31);
    const auto *ik_32 = buffer.data(ik + 32);
    const auto *ik_33 = buffer.data(ik + 33);
    const auto *ik_34 = buffer.data(ik + 34);
    const auto *ik_35 = buffer.data(ik + 35);
    const auto *ik_36 = buffer.data(ik + 36);
    const auto *ik_37 = buffer.data(ik + 37);
    const auto *ik_38 = buffer.data(ik + 38);
    const auto *ik_39 = buffer.data(ik + 39);
    const auto *ik_40 = buffer.data(ik + 40);
    const auto *ik_41 = buffer.data(ik + 41);
    const auto *ik_42 = buffer.data(ik + 42);
    const auto *ik_43 = buffer.data(ik + 43);
    const auto *ik_44 = buffer.data(ik + 44);
    const auto *ik_45 = buffer.data(ik + 45);
    const auto *ik_46 = buffer.data(ik + 46);
    const auto *ik_47 = buffer.data(ik + 47);
    const auto *ik_48 = buffer.data(ik + 48);
    const auto *ik_49 = buffer.data(ik + 49);
    const auto *ik_50 = buffer.data(ik + 50);
    const auto *ik_51 = buffer.data(ik + 51);
    const auto *ik_52 = buffer.data(ik + 52);
    const auto *ik_53 = buffer.data(ik + 53);
    const auto *ik_54 = buffer.data(ik + 54);
    const auto *ik_55 = buffer.data(ik + 55);
    const auto *ik_56 = buffer.data(ik + 56);
    const auto *ik_57 = buffer.data(ik + 57);
    const auto *ik_58 = buffer.data(ik + 58);
    const auto *ik_59 = buffer.data(ik + 59);
    const auto *ik_60 = buffer.data(ik + 60);
    const auto *ik_61 = buffer.data(ik + 61);
    const auto *ik_62 = buffer.data(ik + 62);
    const auto *ik_63 = buffer.data(ik + 63);
    const auto *ik_64 = buffer.data(ik + 64);
    const auto *ik_65 = buffer.data(ik + 65);
    const auto *ik_66 = buffer.data(ik + 66);
    const auto *ik_67 = buffer.data(ik + 67);
    const auto *ik_68 = buffer.data(ik + 68);
    const auto *ik_69 = buffer.data(ik + 69);
    const auto *ik_70 = buffer.data(ik + 70);
    const auto *ik_71 = buffer.data(ik + 71);
    const auto *ik_72 = buffer.data(ik + 72);
    const auto *ik_73 = buffer.data(ik + 73);
    const auto *ik_74 = buffer.data(ik + 74);
    const auto *ik_75 = buffer.data(ik + 75);
    const auto *ik_76 = buffer.data(ik + 76);
    const auto *ik_77 = buffer.data(ik + 77);
    const auto *ik_78 = buffer.data(ik + 78);
    const auto *ik_79 = buffer.data(ik + 79);
    const auto *ik_80 = buffer.data(ik + 80);
    const auto *ik_81 = buffer.data(ik + 81);
    const auto *ik_82 = buffer.data(ik + 82);
    const auto *ik_83 = buffer.data(ik + 83);
    const auto *ik_84 = buffer.data(ik + 84);
    const auto *ik_85 = buffer.data(ik + 85);
    const auto *ik_86 = buffer.data(ik + 86);
    const auto *ik_87 = buffer.data(ik + 87);
    const auto *ik_88 = buffer.data(ik + 88);
    const auto *ik_89 = buffer.data(ik + 89);
    const auto *ik_90 = buffer.data(ik + 90);
    const auto *ik_91 = buffer.data(ik + 91);
    const auto *ik_92 = buffer.data(ik + 92);
    const auto *ik_93 = buffer.data(ik + 93);
    const auto *ik_94 = buffer.data(ik + 94);
    const auto *ik_95 = buffer.data(ik + 95);
    const auto *ik_96 = buffer.data(ik + 96);
    const auto *ik_97 = buffer.data(ik + 97);
    const auto *ik_98 = buffer.data(ik + 98);
    const auto *ik_99 = buffer.data(ik + 99);

    const auto *lk_36 = buffer.data(lk + 36);
    const auto *lk_37 = buffer.data(lk + 37);
    const auto *lk_38 = buffer.data(lk + 38);
    const auto *lk_39 = buffer.data(lk + 39);
    const auto *lk_40 = buffer.data(lk + 40);
    const auto *lk_41 = buffer.data(lk + 41);
    const auto *lk_42 = buffer.data(lk + 42);
    const auto *lk_43 = buffer.data(lk + 43);
    const auto *lk_44 = buffer.data(lk + 44);
    const auto *lk_45 = buffer.data(lk + 45);
    const auto *lk_46 = buffer.data(lk + 46);
    const auto *lk_47 = buffer.data(lk + 47);
    const auto *lk_48 = buffer.data(lk + 48);
    const auto *lk_49 = buffer.data(lk + 49);
    const auto *lk_50 = buffer.data(lk + 50);
    const auto *lk_51 = buffer.data(lk + 51);
    const auto *lk_52 = buffer.data(lk + 52);
    const auto *lk_53 = buffer.data(lk + 53);
    const auto *lk_54 = buffer.data(lk + 54);
    const auto *lk_55 = buffer.data(lk + 55);
    const auto *lk_56 = buffer.data(lk + 56);
    const auto *lk_57 = buffer.data(lk + 57);
    const auto *lk_58 = buffer.data(lk + 58);
    const auto *lk_59 = buffer.data(lk + 59);
    const auto *lk_60 = buffer.data(lk + 60);
    const auto *lk_61 = buffer.data(lk + 61);
    const auto *lk_62 = buffer.data(lk + 62);
    const auto *lk_63 = buffer.data(lk + 63);
    const auto *lk_64 = buffer.data(lk + 64);
    const auto *lk_65 = buffer.data(lk + 65);
    const auto *lk_66 = buffer.data(lk + 66);
    const auto *lk_67 = buffer.data(lk + 67);
    const auto *lk_68 = buffer.data(lk + 68);
    const auto *lk_69 = buffer.data(lk + 69);
    const auto *lk_70 = buffer.data(lk + 70);
    const auto *lk_71 = buffer.data(lk + 71);
    const auto *lk_108 = buffer.data(lk + 108);
    const auto *lk_109 = buffer.data(lk + 109);
    const auto *lk_110 = buffer.data(lk + 110);
    const auto *lk_111 = buffer.data(lk + 111);
    const auto *lk_112 = buffer.data(lk + 112);
    const auto *lk_113 = buffer.data(lk + 113);
    const auto *lk_114 = buffer.data(lk + 114);
    const auto *lk_115 = buffer.data(lk + 115);
    const auto *lk_116 = buffer.data(lk + 116);
    const auto *lk_117 = buffer.data(lk + 117);
    const auto *lk_118 = buffer.data(lk + 118);
    const auto *lk_119 = buffer.data(lk + 119);
    const auto *lk_120 = buffer.data(lk + 120);
    const auto *lk_121 = buffer.data(lk + 121);
    const auto *lk_122 = buffer.data(lk + 122);
    const auto *lk_123 = buffer.data(lk + 123);
    const auto *lk_124 = buffer.data(lk + 124);
    const auto *lk_125 = buffer.data(lk + 125);
    const auto *lk_126 = buffer.data(lk + 126);
    const auto *lk_127 = buffer.data(lk + 127);
    const auto *lk_128 = buffer.data(lk + 128);
    const auto *lk_129 = buffer.data(lk + 129);
    const auto *lk_130 = buffer.data(lk + 130);
    const auto *lk_131 = buffer.data(lk + 131);
    const auto *lk_132 = buffer.data(lk + 132);
    const auto *lk_133 = buffer.data(lk + 133);
    const auto *lk_134 = buffer.data(lk + 134);
    const auto *lk_135 = buffer.data(lk + 135);
    const auto *lk_136 = buffer.data(lk + 136);
    const auto *lk_137 = buffer.data(lk + 137);
    const auto *lk_138 = buffer.data(lk + 138);
    const auto *lk_139 = buffer.data(lk + 139);
    const auto *lk_140 = buffer.data(lk + 140);
    const auto *lk_141 = buffer.data(lk + 141);
    const auto *lk_142 = buffer.data(lk + 142);
    const auto *lk_143 = buffer.data(lk + 143);
    const auto *lk_144 = buffer.data(lk + 144);
    const auto *lk_145 = buffer.data(lk + 145);
    const auto *lk_146 = buffer.data(lk + 146);
    const auto *lk_147 = buffer.data(lk + 147);
    const auto *lk_148 = buffer.data(lk + 148);
    const auto *lk_149 = buffer.data(lk + 149);
    const auto *lk_150 = buffer.data(lk + 150);
    const auto *lk_151 = buffer.data(lk + 151);
    const auto *lk_152 = buffer.data(lk + 152);
    const auto *lk_153 = buffer.data(lk + 153);
    const auto *lk_154 = buffer.data(lk + 154);
    const auto *lk_155 = buffer.data(lk + 155);
    const auto *lk_156 = buffer.data(lk + 156);
    const auto *lk_157 = buffer.data(lk + 157);
    const auto *lk_158 = buffer.data(lk + 158);
    const auto *lk_159 = buffer.data(lk + 159);
    const auto *lk_160 = buffer.data(lk + 160);
    const auto *lk_161 = buffer.data(lk + 161);
    const auto *lk_162 = buffer.data(lk + 162);
    const auto *lk_163 = buffer.data(lk + 163);
    const auto *lk_164 = buffer.data(lk + 164);
    const auto *lk_165 = buffer.data(lk + 165);
    const auto *lk_166 = buffer.data(lk + 166);
    const auto *lk_167 = buffer.data(lk + 167);
    const auto *lk_168 = buffer.data(lk + 168);
    const auto *lk_169 = buffer.data(lk + 169);
    const auto *lk_170 = buffer.data(lk + 170);
    const auto *lk_171 = buffer.data(lk + 171);
    const auto *lk_172 = buffer.data(lk + 172);
    const auto *lk_173 = buffer.data(lk + 173);
    const auto *lk_174 = buffer.data(lk + 174);
    const auto *lk_175 = buffer.data(lk + 175);
    const auto *lk_176 = buffer.data(lk + 176);
    const auto *lk_177 = buffer.data(lk + 177);
    const auto *lk_178 = buffer.data(lk + 178);
    const auto *lk_179 = buffer.data(lk + 179);
    const auto *lk_216 = buffer.data(lk + 216);
    const auto *lk_217 = buffer.data(lk + 217);
    const auto *lk_218 = buffer.data(lk + 218);
    const auto *lk_219 = buffer.data(lk + 219);
    const auto *lk_220 = buffer.data(lk + 220);
    const auto *lk_221 = buffer.data(lk + 221);
    const auto *lk_222 = buffer.data(lk + 222);
    const auto *lk_223 = buffer.data(lk + 223);
    const auto *lk_224 = buffer.data(lk + 224);
    const auto *lk_225 = buffer.data(lk + 225);
    const auto *lk_226 = buffer.data(lk + 226);
    const auto *lk_227 = buffer.data(lk + 227);
    const auto *lk_228 = buffer.data(lk + 228);
    const auto *lk_229 = buffer.data(lk + 229);
    const auto *lk_230 = buffer.data(lk + 230);
    const auto *lk_231 = buffer.data(lk + 231);
    const auto *lk_232 = buffer.data(lk + 232);
    const auto *lk_233 = buffer.data(lk + 233);
    const auto *lk_234 = buffer.data(lk + 234);
    const auto *lk_235 = buffer.data(lk + 235);
    const auto *lk_236 = buffer.data(lk + 236);
    const auto *lk_237 = buffer.data(lk + 237);
    const auto *lk_238 = buffer.data(lk + 238);
    const auto *lk_239 = buffer.data(lk + 239);
    const auto *lk_240 = buffer.data(lk + 240);
    const auto *lk_241 = buffer.data(lk + 241);
    const auto *lk_242 = buffer.data(lk + 242);
    const auto *lk_243 = buffer.data(lk + 243);
    const auto *lk_244 = buffer.data(lk + 244);
    const auto *lk_245 = buffer.data(lk + 245);
    const auto *lk_246 = buffer.data(lk + 246);
    const auto *lk_247 = buffer.data(lk + 247);
    const auto *lk_248 = buffer.data(lk + 248);
    const auto *lk_249 = buffer.data(lk + 249);
    const auto *lk_250 = buffer.data(lk + 250);
    const auto *lk_251 = buffer.data(lk + 251);
    const auto *lk_252 = buffer.data(lk + 252);
    const auto *lk_253 = buffer.data(lk + 253);
    const auto *lk_254 = buffer.data(lk + 254);
    const auto *lk_255 = buffer.data(lk + 255);
    const auto *lk_256 = buffer.data(lk + 256);
    const auto *lk_257 = buffer.data(lk + 257);
    const auto *lk_258 = buffer.data(lk + 258);
    const auto *lk_259 = buffer.data(lk + 259);
    const auto *lk_260 = buffer.data(lk + 260);
    const auto *lk_261 = buffer.data(lk + 261);
    const auto *lk_262 = buffer.data(lk + 262);
    const auto *lk_263 = buffer.data(lk + 263);
    const auto *lk_264 = buffer.data(lk + 264);
    const auto *lk_265 = buffer.data(lk + 265);
    const auto *lk_266 = buffer.data(lk + 266);
    const auto *lk_267 = buffer.data(lk + 267);
    const auto *lk_268 = buffer.data(lk + 268);
    const auto *lk_269 = buffer.data(lk + 269);
    const auto *lk_270 = buffer.data(lk + 270);
    const auto *lk_271 = buffer.data(lk + 271);
    const auto *lk_272 = buffer.data(lk + 272);
    const auto *lk_273 = buffer.data(lk + 273);
    const auto *lk_274 = buffer.data(lk + 274);
    const auto *lk_275 = buffer.data(lk + 275);
    const auto *lk_276 = buffer.data(lk + 276);
    const auto *lk_277 = buffer.data(lk + 277);
    const auto *lk_278 = buffer.data(lk + 278);
    const auto *lk_279 = buffer.data(lk + 279);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, lk_36, lk_37, lk_38, lk_39, \
                         lk_40, lk_41, lk_42, lk_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lk_36[k];

        t_1[k] = f_0 * lk_37[k];

        t_2[k] = f_0 * lk_38[k];

        t_3[k] = f_0 * lk_39[k];

        t_4[k] = f_0 * lk_40[k];

        t_5[k] = f_0 * lk_41[k];

        t_6[k] = f_0 * lk_42[k];

        t_7[k] = f_0 * lk_43[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, lk_44, lk_45, lk_46, \
                         lk_47, lk_48, lk_49, lk_50, lk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * lk_44[k];

        t_9[k] = f_0 * lk_45[k];

        t_10[k] = f_0 * lk_46[k];

        t_11[k] = f_0 * lk_47[k];

        t_12[k] = f_0 * lk_48[k];

        t_13[k] = f_0 * lk_49[k];

        t_14[k] = f_0 * lk_50[k];

        t_15[k] = f_0 * lk_51[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, lk_52, lk_53, lk_54, \
                         lk_55, lk_56, lk_57, lk_58, lk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * lk_52[k];

        t_17[k] = f_0 * lk_53[k];

        t_18[k] = f_0 * lk_54[k];

        t_19[k] = f_0 * lk_55[k];

        t_20[k] = f_0 * lk_56[k];

        t_21[k] = f_0 * lk_57[k];

        t_22[k] = f_0 * lk_58[k];

        t_23[k] = f_0 * lk_59[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, lk_60, lk_61, lk_62, \
                         lk_63, lk_64, lk_65, lk_66, lk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * lk_60[k];

        t_25[k] = f_0 * lk_61[k];

        t_26[k] = f_0 * lk_62[k];

        t_27[k] = f_0 * lk_63[k];

        t_28[k] = f_0 * lk_64[k];

        t_29[k] = f_0 * lk_65[k];

        t_30[k] = f_0 * lk_66[k];

        t_31[k] = f_0 * lk_67[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ik_0, ik_1, lk_68, lk_69, lk_70, \
                         lk_71, lk_108, lk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * lk_68[k];

        t_33[k] = f_0 * lk_69[k];

        t_34[k] = f_0 * lk_70[k];

        t_35[k] = f_0 * lk_71[k];

        t_36[k] = -ik_0[k]
                  + f_0 * lk_108[k];

        t_37[k] = -ik_1[k]
                  + f_0 * lk_109[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, ik_2, ik_3, ik_4, ik_5, ik_6, lk_110, \
                         lk_111, lk_112, lk_113, lk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -ik_2[k]
                  + f_0 * lk_110[k];

        t_39[k] = -ik_3[k]
                  + f_0 * lk_111[k];

        t_40[k] = -ik_4[k]
                  + f_0 * lk_112[k];

        t_41[k] = -ik_5[k]
                  + f_0 * lk_113[k];

        t_42[k] = -ik_6[k]
                  + f_0 * lk_114[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, ik_7, ik_8, ik_9, ik_10, ik_11, lk_115, \
                         lk_116, lk_117, lk_118, lk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -ik_7[k]
                  + f_0 * lk_115[k];

        t_44[k] = -ik_8[k]
                  + f_0 * lk_116[k];

        t_45[k] = -ik_9[k]
                  + f_0 * lk_117[k];

        t_46[k] = -ik_10[k]
                  + f_0 * lk_118[k];

        t_47[k] = -ik_11[k]
                  + f_0 * lk_119[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, ik_12, ik_13, ik_14, ik_15, ik_16, \
                         lk_120, lk_121, lk_122, lk_123, lk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -ik_12[k]
                  + f_0 * lk_120[k];

        t_49[k] = -ik_13[k]
                  + f_0 * lk_121[k];

        t_50[k] = -ik_14[k]
                  + f_0 * lk_122[k];

        t_51[k] = -ik_15[k]
                  + f_0 * lk_123[k];

        t_52[k] = -ik_16[k]
                  + f_0 * lk_124[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, ik_17, ik_18, ik_19, ik_20, ik_21, \
                         lk_125, lk_126, lk_127, lk_128, lk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -ik_17[k]
                  + f_0 * lk_125[k];

        t_54[k] = -ik_18[k]
                  + f_0 * lk_126[k];

        t_55[k] = -ik_19[k]
                  + f_0 * lk_127[k];

        t_56[k] = -ik_20[k]
                  + f_0 * lk_128[k];

        t_57[k] = -ik_21[k]
                  + f_0 * lk_129[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, ik_22, ik_23, ik_24, ik_25, ik_26, \
                         lk_130, lk_131, lk_132, lk_133, lk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -ik_22[k]
                  + f_0 * lk_130[k];

        t_59[k] = -ik_23[k]
                  + f_0 * lk_131[k];

        t_60[k] = -ik_24[k]
                  + f_0 * lk_132[k];

        t_61[k] = -ik_25[k]
                  + f_0 * lk_133[k];

        t_62[k] = -ik_26[k]
                  + f_0 * lk_134[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, ik_27, ik_28, ik_29, ik_30, ik_31, \
                         lk_135, lk_136, lk_137, lk_138, lk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -ik_27[k]
                  + f_0 * lk_135[k];

        t_64[k] = -ik_28[k]
                  + f_0 * lk_136[k];

        t_65[k] = -ik_29[k]
                  + f_0 * lk_137[k];

        t_66[k] = -ik_30[k]
                  + f_0 * lk_138[k];

        t_67[k] = -ik_31[k]
                  + f_0 * lk_139[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, ik_32, ik_33, ik_34, ik_35, \
                         lk_140, lk_141, lk_142, lk_143, lk_144, \
                         lk_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -ik_32[k]
                  + f_0 * lk_140[k];

        t_69[k] = -ik_33[k]
                  + f_0 * lk_141[k];

        t_70[k] = -ik_34[k]
                  + f_0 * lk_142[k];

        t_71[k] = -ik_35[k]
                  + f_0 * lk_143[k];

        t_72[k] = f_0 * lk_144[k];

        t_73[k] = f_0 * lk_145[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, t_80, t_81, lk_146, lk_147, \
                         lk_148, lk_149, lk_150, lk_151, lk_152, \
                         lk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * lk_146[k];

        t_75[k] = f_0 * lk_147[k];

        t_76[k] = f_0 * lk_148[k];

        t_77[k] = f_0 * lk_149[k];

        t_78[k] = f_0 * lk_150[k];

        t_79[k] = f_0 * lk_151[k];

        t_80[k] = f_0 * lk_152[k];

        t_81[k] = f_0 * lk_153[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, t_87, t_88, t_89, lk_154, lk_155, \
                         lk_156, lk_157, lk_158, lk_159, lk_160, \
                         lk_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * lk_154[k];

        t_83[k] = f_0 * lk_155[k];

        t_84[k] = f_0 * lk_156[k];

        t_85[k] = f_0 * lk_157[k];

        t_86[k] = f_0 * lk_158[k];

        t_87[k] = f_0 * lk_159[k];

        t_88[k] = f_0 * lk_160[k];

        t_89[k] = f_0 * lk_161[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, lk_162, lk_163, \
                         lk_164, lk_165, lk_166, lk_167, lk_168, \
                         lk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * lk_162[k];

        t_91[k] = f_0 * lk_163[k];

        t_92[k] = f_0 * lk_164[k];

        t_93[k] = f_0 * lk_165[k];

        t_94[k] = f_0 * lk_166[k];

        t_95[k] = f_0 * lk_167[k];

        t_96[k] = f_0 * lk_168[k];

        t_97[k] = f_0 * lk_169[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, t_104, t_105, lk_170, lk_171, \
                         lk_172, lk_173, lk_174, lk_175, lk_176, \
                         lk_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * lk_170[k];

        t_99[k] = f_0 * lk_171[k];

        t_100[k] = f_0 * lk_172[k];

        t_101[k] = f_0 * lk_173[k];

        t_102[k] = f_0 * lk_174[k];

        t_103[k] = f_0 * lk_175[k];

        t_104[k] = f_0 * lk_176[k];

        t_105[k] = f_0 * lk_177[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, ik_36, ik_37, ik_38, ik_39, \
                         lk_178, lk_179, lk_216, lk_217, lk_218, \
                         lk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * lk_178[k];

        t_107[k] = f_0 * lk_179[k];

        t_108[k] = -2.0 * ik_36[k]
                   + f_0 * lk_216[k];

        t_109[k] = -2.0 * ik_37[k]
                   + f_0 * lk_217[k];

        t_110[k] = -2.0 * ik_38[k]
                   + f_0 * lk_218[k];

        t_111[k] = -2.0 * ik_39[k]
                   + f_0 * lk_219[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, ik_40, ik_41, ik_42, ik_43, ik_44, \
                         lk_220, lk_221, lk_222, lk_223, lk_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -2.0 * ik_40[k]
                   + f_0 * lk_220[k];

        t_113[k] = -2.0 * ik_41[k]
                   + f_0 * lk_221[k];

        t_114[k] = -2.0 * ik_42[k]
                   + f_0 * lk_222[k];

        t_115[k] = -2.0 * ik_43[k]
                   + f_0 * lk_223[k];

        t_116[k] = -2.0 * ik_44[k]
                   + f_0 * lk_224[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, ik_45, ik_46, ik_47, ik_48, ik_49, \
                         lk_225, lk_226, lk_227, lk_228, lk_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -2.0 * ik_45[k]
                   + f_0 * lk_225[k];

        t_118[k] = -2.0 * ik_46[k]
                   + f_0 * lk_226[k];

        t_119[k] = -2.0 * ik_47[k]
                   + f_0 * lk_227[k];

        t_120[k] = -2.0 * ik_48[k]
                   + f_0 * lk_228[k];

        t_121[k] = -2.0 * ik_49[k]
                   + f_0 * lk_229[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, ik_50, ik_51, ik_52, ik_53, ik_54, \
                         lk_230, lk_231, lk_232, lk_233, lk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * ik_50[k]
                   + f_0 * lk_230[k];

        t_123[k] = -2.0 * ik_51[k]
                   + f_0 * lk_231[k];

        t_124[k] = -2.0 * ik_52[k]
                   + f_0 * lk_232[k];

        t_125[k] = -2.0 * ik_53[k]
                   + f_0 * lk_233[k];

        t_126[k] = -2.0 * ik_54[k]
                   + f_0 * lk_234[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, ik_55, ik_56, ik_57, ik_58, ik_59, \
                         lk_235, lk_236, lk_237, lk_238, lk_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * ik_55[k]
                   + f_0 * lk_235[k];

        t_128[k] = -2.0 * ik_56[k]
                   + f_0 * lk_236[k];

        t_129[k] = -2.0 * ik_57[k]
                   + f_0 * lk_237[k];

        t_130[k] = -2.0 * ik_58[k]
                   + f_0 * lk_238[k];

        t_131[k] = -2.0 * ik_59[k]
                   + f_0 * lk_239[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, ik_60, ik_61, ik_62, ik_63, ik_64, \
                         lk_240, lk_241, lk_242, lk_243, lk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * ik_60[k]
                   + f_0 * lk_240[k];

        t_133[k] = -2.0 * ik_61[k]
                   + f_0 * lk_241[k];

        t_134[k] = -2.0 * ik_62[k]
                   + f_0 * lk_242[k];

        t_135[k] = -2.0 * ik_63[k]
                   + f_0 * lk_243[k];

        t_136[k] = -2.0 * ik_64[k]
                   + f_0 * lk_244[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, ik_65, ik_66, ik_67, ik_68, ik_69, \
                         lk_245, lk_246, lk_247, lk_248, lk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -2.0 * ik_65[k]
                   + f_0 * lk_245[k];

        t_138[k] = -2.0 * ik_66[k]
                   + f_0 * lk_246[k];

        t_139[k] = -2.0 * ik_67[k]
                   + f_0 * lk_247[k];

        t_140[k] = -2.0 * ik_68[k]
                   + f_0 * lk_248[k];

        t_141[k] = -2.0 * ik_69[k]
                   + f_0 * lk_249[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, ik_70, ik_71, ik_72, ik_73, ik_74, \
                         lk_250, lk_251, lk_252, lk_253, lk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -2.0 * ik_70[k]
                   + f_0 * lk_250[k];

        t_143[k] = -2.0 * ik_71[k]
                   + f_0 * lk_251[k];

        t_144[k] = -ik_72[k]
                   + f_0 * lk_252[k];

        t_145[k] = -ik_73[k]
                   + f_0 * lk_253[k];

        t_146[k] = -ik_74[k]
                   + f_0 * lk_254[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, ik_75, ik_76, ik_77, ik_78, ik_79, \
                         lk_255, lk_256, lk_257, lk_258, lk_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -ik_75[k]
                   + f_0 * lk_255[k];

        t_148[k] = -ik_76[k]
                   + f_0 * lk_256[k];

        t_149[k] = -ik_77[k]
                   + f_0 * lk_257[k];

        t_150[k] = -ik_78[k]
                   + f_0 * lk_258[k];

        t_151[k] = -ik_79[k]
                   + f_0 * lk_259[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, ik_80, ik_81, ik_82, ik_83, ik_84, \
                         lk_260, lk_261, lk_262, lk_263, lk_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -ik_80[k]
                   + f_0 * lk_260[k];

        t_153[k] = -ik_81[k]
                   + f_0 * lk_261[k];

        t_154[k] = -ik_82[k]
                   + f_0 * lk_262[k];

        t_155[k] = -ik_83[k]
                   + f_0 * lk_263[k];

        t_156[k] = -ik_84[k]
                   + f_0 * lk_264[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, ik_85, ik_86, ik_87, ik_88, ik_89, \
                         lk_265, lk_266, lk_267, lk_268, lk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -ik_85[k]
                   + f_0 * lk_265[k];

        t_158[k] = -ik_86[k]
                   + f_0 * lk_266[k];

        t_159[k] = -ik_87[k]
                   + f_0 * lk_267[k];

        t_160[k] = -ik_88[k]
                   + f_0 * lk_268[k];

        t_161[k] = -ik_89[k]
                   + f_0 * lk_269[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, ik_90, ik_91, ik_92, ik_93, ik_94, \
                         lk_270, lk_271, lk_272, lk_273, lk_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -ik_90[k]
                   + f_0 * lk_270[k];

        t_163[k] = -ik_91[k]
                   + f_0 * lk_271[k];

        t_164[k] = -ik_92[k]
                   + f_0 * lk_272[k];

        t_165[k] = -ik_93[k]
                   + f_0 * lk_273[k];

        t_166[k] = -ik_94[k]
                   + f_0 * lk_274[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, ik_95, ik_96, ik_97, ik_98, ik_99, \
                         lk_275, lk_276, lk_277, lk_278, lk_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -ik_95[k]
                   + f_0 * lk_275[k];

        t_168[k] = -ik_96[k]
                   + f_0 * lk_276[k];

        t_169[k] = -ik_97[k]
                   + f_0 * lk_277[k];

        t_170[k] = -ik_98[k]
                   + f_0 * lk_278[k];

        t_171[k] = -ik_99[k]
                   + f_0 * lk_279[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_100 = buffer.data(ik + 100);
    const auto *ik_101 = buffer.data(ik + 101);
    const auto *ik_102 = buffer.data(ik + 102);
    const auto *ik_103 = buffer.data(ik + 103);
    const auto *ik_104 = buffer.data(ik + 104);
    const auto *ik_105 = buffer.data(ik + 105);
    const auto *ik_106 = buffer.data(ik + 106);
    const auto *ik_107 = buffer.data(ik + 107);
    const auto *ik_108 = buffer.data(ik + 108);
    const auto *ik_109 = buffer.data(ik + 109);
    const auto *ik_110 = buffer.data(ik + 110);
    const auto *ik_111 = buffer.data(ik + 111);
    const auto *ik_112 = buffer.data(ik + 112);
    const auto *ik_113 = buffer.data(ik + 113);
    const auto *ik_114 = buffer.data(ik + 114);
    const auto *ik_115 = buffer.data(ik + 115);
    const auto *ik_116 = buffer.data(ik + 116);
    const auto *ik_117 = buffer.data(ik + 117);
    const auto *ik_118 = buffer.data(ik + 118);
    const auto *ik_119 = buffer.data(ik + 119);
    const auto *ik_120 = buffer.data(ik + 120);
    const auto *ik_121 = buffer.data(ik + 121);
    const auto *ik_122 = buffer.data(ik + 122);
    const auto *ik_123 = buffer.data(ik + 123);
    const auto *ik_124 = buffer.data(ik + 124);
    const auto *ik_125 = buffer.data(ik + 125);
    const auto *ik_126 = buffer.data(ik + 126);
    const auto *ik_127 = buffer.data(ik + 127);
    const auto *ik_128 = buffer.data(ik + 128);
    const auto *ik_129 = buffer.data(ik + 129);
    const auto *ik_130 = buffer.data(ik + 130);
    const auto *ik_131 = buffer.data(ik + 131);
    const auto *ik_132 = buffer.data(ik + 132);
    const auto *ik_133 = buffer.data(ik + 133);
    const auto *ik_134 = buffer.data(ik + 134);
    const auto *ik_135 = buffer.data(ik + 135);
    const auto *ik_136 = buffer.data(ik + 136);
    const auto *ik_137 = buffer.data(ik + 137);
    const auto *ik_138 = buffer.data(ik + 138);
    const auto *ik_139 = buffer.data(ik + 139);
    const auto *ik_140 = buffer.data(ik + 140);
    const auto *ik_141 = buffer.data(ik + 141);
    const auto *ik_142 = buffer.data(ik + 142);
    const auto *ik_143 = buffer.data(ik + 143);
    const auto *ik_144 = buffer.data(ik + 144);
    const auto *ik_145 = buffer.data(ik + 145);
    const auto *ik_146 = buffer.data(ik + 146);
    const auto *ik_147 = buffer.data(ik + 147);
    const auto *ik_148 = buffer.data(ik + 148);
    const auto *ik_149 = buffer.data(ik + 149);
    const auto *ik_150 = buffer.data(ik + 150);
    const auto *ik_151 = buffer.data(ik + 151);
    const auto *ik_152 = buffer.data(ik + 152);
    const auto *ik_153 = buffer.data(ik + 153);
    const auto *ik_154 = buffer.data(ik + 154);
    const auto *ik_155 = buffer.data(ik + 155);
    const auto *ik_156 = buffer.data(ik + 156);
    const auto *ik_157 = buffer.data(ik + 157);
    const auto *ik_158 = buffer.data(ik + 158);
    const auto *ik_159 = buffer.data(ik + 159);
    const auto *ik_160 = buffer.data(ik + 160);
    const auto *ik_161 = buffer.data(ik + 161);
    const auto *ik_162 = buffer.data(ik + 162);
    const auto *ik_163 = buffer.data(ik + 163);
    const auto *ik_164 = buffer.data(ik + 164);
    const auto *ik_165 = buffer.data(ik + 165);
    const auto *ik_166 = buffer.data(ik + 166);
    const auto *ik_167 = buffer.data(ik + 167);
    const auto *ik_168 = buffer.data(ik + 168);
    const auto *ik_169 = buffer.data(ik + 169);
    const auto *ik_170 = buffer.data(ik + 170);
    const auto *ik_171 = buffer.data(ik + 171);
    const auto *ik_172 = buffer.data(ik + 172);
    const auto *ik_173 = buffer.data(ik + 173);
    const auto *ik_174 = buffer.data(ik + 174);
    const auto *ik_175 = buffer.data(ik + 175);
    const auto *ik_176 = buffer.data(ik + 176);
    const auto *ik_177 = buffer.data(ik + 177);
    const auto *ik_178 = buffer.data(ik + 178);
    const auto *ik_179 = buffer.data(ik + 179);
    const auto *ik_180 = buffer.data(ik + 180);
    const auto *ik_181 = buffer.data(ik + 181);
    const auto *ik_182 = buffer.data(ik + 182);
    const auto *ik_183 = buffer.data(ik + 183);
    const auto *ik_184 = buffer.data(ik + 184);
    const auto *ik_185 = buffer.data(ik + 185);
    const auto *ik_186 = buffer.data(ik + 186);
    const auto *ik_187 = buffer.data(ik + 187);
    const auto *ik_188 = buffer.data(ik + 188);
    const auto *ik_189 = buffer.data(ik + 189);
    const auto *ik_190 = buffer.data(ik + 190);
    const auto *ik_191 = buffer.data(ik + 191);
    const auto *ik_192 = buffer.data(ik + 192);
    const auto *ik_193 = buffer.data(ik + 193);
    const auto *ik_194 = buffer.data(ik + 194);
    const auto *ik_195 = buffer.data(ik + 195);
    const auto *ik_196 = buffer.data(ik + 196);
    const auto *ik_197 = buffer.data(ik + 197);
    const auto *ik_198 = buffer.data(ik + 198);
    const auto *ik_199 = buffer.data(ik + 199);
    const auto *ik_200 = buffer.data(ik + 200);
    const auto *ik_201 = buffer.data(ik + 201);
    const auto *ik_202 = buffer.data(ik + 202);
    const auto *ik_203 = buffer.data(ik + 203);
    const auto *ik_204 = buffer.data(ik + 204);
    const auto *ik_205 = buffer.data(ik + 205);
    const auto *ik_206 = buffer.data(ik + 206);
    const auto *ik_207 = buffer.data(ik + 207);
    const auto *ik_208 = buffer.data(ik + 208);
    const auto *ik_209 = buffer.data(ik + 209);
    const auto *ik_210 = buffer.data(ik + 210);
    const auto *ik_211 = buffer.data(ik + 211);
    const auto *ik_212 = buffer.data(ik + 212);
    const auto *ik_213 = buffer.data(ik + 213);
    const auto *ik_214 = buffer.data(ik + 214);
    const auto *ik_215 = buffer.data(ik + 215);

    const auto *lk_280 = buffer.data(lk + 280);
    const auto *lk_281 = buffer.data(lk + 281);
    const auto *lk_282 = buffer.data(lk + 282);
    const auto *lk_283 = buffer.data(lk + 283);
    const auto *lk_284 = buffer.data(lk + 284);
    const auto *lk_285 = buffer.data(lk + 285);
    const auto *lk_286 = buffer.data(lk + 286);
    const auto *lk_287 = buffer.data(lk + 287);
    const auto *lk_288 = buffer.data(lk + 288);
    const auto *lk_289 = buffer.data(lk + 289);
    const auto *lk_290 = buffer.data(lk + 290);
    const auto *lk_291 = buffer.data(lk + 291);
    const auto *lk_292 = buffer.data(lk + 292);
    const auto *lk_293 = buffer.data(lk + 293);
    const auto *lk_294 = buffer.data(lk + 294);
    const auto *lk_295 = buffer.data(lk + 295);
    const auto *lk_296 = buffer.data(lk + 296);
    const auto *lk_297 = buffer.data(lk + 297);
    const auto *lk_298 = buffer.data(lk + 298);
    const auto *lk_299 = buffer.data(lk + 299);
    const auto *lk_300 = buffer.data(lk + 300);
    const auto *lk_301 = buffer.data(lk + 301);
    const auto *lk_302 = buffer.data(lk + 302);
    const auto *lk_303 = buffer.data(lk + 303);
    const auto *lk_304 = buffer.data(lk + 304);
    const auto *lk_305 = buffer.data(lk + 305);
    const auto *lk_306 = buffer.data(lk + 306);
    const auto *lk_307 = buffer.data(lk + 307);
    const auto *lk_308 = buffer.data(lk + 308);
    const auto *lk_309 = buffer.data(lk + 309);
    const auto *lk_310 = buffer.data(lk + 310);
    const auto *lk_311 = buffer.data(lk + 311);
    const auto *lk_312 = buffer.data(lk + 312);
    const auto *lk_313 = buffer.data(lk + 313);
    const auto *lk_314 = buffer.data(lk + 314);
    const auto *lk_315 = buffer.data(lk + 315);
    const auto *lk_316 = buffer.data(lk + 316);
    const auto *lk_317 = buffer.data(lk + 317);
    const auto *lk_318 = buffer.data(lk + 318);
    const auto *lk_319 = buffer.data(lk + 319);
    const auto *lk_320 = buffer.data(lk + 320);
    const auto *lk_321 = buffer.data(lk + 321);
    const auto *lk_322 = buffer.data(lk + 322);
    const auto *lk_323 = buffer.data(lk + 323);
    const auto *lk_360 = buffer.data(lk + 360);
    const auto *lk_361 = buffer.data(lk + 361);
    const auto *lk_362 = buffer.data(lk + 362);
    const auto *lk_363 = buffer.data(lk + 363);
    const auto *lk_364 = buffer.data(lk + 364);
    const auto *lk_365 = buffer.data(lk + 365);
    const auto *lk_366 = buffer.data(lk + 366);
    const auto *lk_367 = buffer.data(lk + 367);
    const auto *lk_368 = buffer.data(lk + 368);
    const auto *lk_369 = buffer.data(lk + 369);
    const auto *lk_370 = buffer.data(lk + 370);
    const auto *lk_371 = buffer.data(lk + 371);
    const auto *lk_372 = buffer.data(lk + 372);
    const auto *lk_373 = buffer.data(lk + 373);
    const auto *lk_374 = buffer.data(lk + 374);
    const auto *lk_375 = buffer.data(lk + 375);
    const auto *lk_376 = buffer.data(lk + 376);
    const auto *lk_377 = buffer.data(lk + 377);
    const auto *lk_378 = buffer.data(lk + 378);
    const auto *lk_379 = buffer.data(lk + 379);
    const auto *lk_380 = buffer.data(lk + 380);
    const auto *lk_381 = buffer.data(lk + 381);
    const auto *lk_382 = buffer.data(lk + 382);
    const auto *lk_383 = buffer.data(lk + 383);
    const auto *lk_384 = buffer.data(lk + 384);
    const auto *lk_385 = buffer.data(lk + 385);
    const auto *lk_386 = buffer.data(lk + 386);
    const auto *lk_387 = buffer.data(lk + 387);
    const auto *lk_388 = buffer.data(lk + 388);
    const auto *lk_389 = buffer.data(lk + 389);
    const auto *lk_390 = buffer.data(lk + 390);
    const auto *lk_391 = buffer.data(lk + 391);
    const auto *lk_392 = buffer.data(lk + 392);
    const auto *lk_393 = buffer.data(lk + 393);
    const auto *lk_394 = buffer.data(lk + 394);
    const auto *lk_395 = buffer.data(lk + 395);
    const auto *lk_396 = buffer.data(lk + 396);
    const auto *lk_397 = buffer.data(lk + 397);
    const auto *lk_398 = buffer.data(lk + 398);
    const auto *lk_399 = buffer.data(lk + 399);
    const auto *lk_400 = buffer.data(lk + 400);
    const auto *lk_401 = buffer.data(lk + 401);
    const auto *lk_402 = buffer.data(lk + 402);
    const auto *lk_403 = buffer.data(lk + 403);
    const auto *lk_404 = buffer.data(lk + 404);
    const auto *lk_405 = buffer.data(lk + 405);
    const auto *lk_406 = buffer.data(lk + 406);
    const auto *lk_407 = buffer.data(lk + 407);
    const auto *lk_408 = buffer.data(lk + 408);
    const auto *lk_409 = buffer.data(lk + 409);
    const auto *lk_410 = buffer.data(lk + 410);
    const auto *lk_411 = buffer.data(lk + 411);
    const auto *lk_412 = buffer.data(lk + 412);
    const auto *lk_413 = buffer.data(lk + 413);
    const auto *lk_414 = buffer.data(lk + 414);
    const auto *lk_415 = buffer.data(lk + 415);
    const auto *lk_416 = buffer.data(lk + 416);
    const auto *lk_417 = buffer.data(lk + 417);
    const auto *lk_418 = buffer.data(lk + 418);
    const auto *lk_419 = buffer.data(lk + 419);
    const auto *lk_420 = buffer.data(lk + 420);
    const auto *lk_421 = buffer.data(lk + 421);
    const auto *lk_422 = buffer.data(lk + 422);
    const auto *lk_423 = buffer.data(lk + 423);
    const auto *lk_424 = buffer.data(lk + 424);
    const auto *lk_425 = buffer.data(lk + 425);
    const auto *lk_426 = buffer.data(lk + 426);
    const auto *lk_427 = buffer.data(lk + 427);
    const auto *lk_428 = buffer.data(lk + 428);
    const auto *lk_429 = buffer.data(lk + 429);
    const auto *lk_430 = buffer.data(lk + 430);
    const auto *lk_431 = buffer.data(lk + 431);
    const auto *lk_432 = buffer.data(lk + 432);
    const auto *lk_433 = buffer.data(lk + 433);
    const auto *lk_434 = buffer.data(lk + 434);
    const auto *lk_435 = buffer.data(lk + 435);
    const auto *lk_436 = buffer.data(lk + 436);
    const auto *lk_437 = buffer.data(lk + 437);
    const auto *lk_438 = buffer.data(lk + 438);
    const auto *lk_439 = buffer.data(lk + 439);
    const auto *lk_440 = buffer.data(lk + 440);
    const auto *lk_441 = buffer.data(lk + 441);
    const auto *lk_442 = buffer.data(lk + 442);
    const auto *lk_443 = buffer.data(lk + 443);
    const auto *lk_444 = buffer.data(lk + 444);
    const auto *lk_445 = buffer.data(lk + 445);
    const auto *lk_446 = buffer.data(lk + 446);
    const auto *lk_447 = buffer.data(lk + 447);
    const auto *lk_448 = buffer.data(lk + 448);
    const auto *lk_449 = buffer.data(lk + 449);
    const auto *lk_450 = buffer.data(lk + 450);
    const auto *lk_451 = buffer.data(lk + 451);
    const auto *lk_452 = buffer.data(lk + 452);
    const auto *lk_453 = buffer.data(lk + 453);
    const auto *lk_454 = buffer.data(lk + 454);
    const auto *lk_455 = buffer.data(lk + 455);
    const auto *lk_456 = buffer.data(lk + 456);
    const auto *lk_457 = buffer.data(lk + 457);
    const auto *lk_458 = buffer.data(lk + 458);
    const auto *lk_459 = buffer.data(lk + 459);
    const auto *lk_460 = buffer.data(lk + 460);
    const auto *lk_461 = buffer.data(lk + 461);
    const auto *lk_462 = buffer.data(lk + 462);
    const auto *lk_463 = buffer.data(lk + 463);
    const auto *lk_464 = buffer.data(lk + 464);
    const auto *lk_465 = buffer.data(lk + 465);
    const auto *lk_466 = buffer.data(lk + 466);
    const auto *lk_467 = buffer.data(lk + 467);
    const auto *lk_468 = buffer.data(lk + 468);
    const auto *lk_469 = buffer.data(lk + 469);
    const auto *lk_470 = buffer.data(lk + 470);
    const auto *lk_471 = buffer.data(lk + 471);
    const auto *lk_472 = buffer.data(lk + 472);
    const auto *lk_473 = buffer.data(lk + 473);
    const auto *lk_474 = buffer.data(lk + 474);
    const auto *lk_475 = buffer.data(lk + 475);
    const auto *lk_476 = buffer.data(lk + 476);
    const auto *lk_477 = buffer.data(lk + 477);

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, ik_100, ik_101, ik_102, ik_103, \
                         ik_104, lk_280, lk_281, lk_282, lk_283, \
                         lk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -ik_100[k]
                   + f_0 * lk_280[k];

        t_173[k] = -ik_101[k]
                   + f_0 * lk_281[k];

        t_174[k] = -ik_102[k]
                   + f_0 * lk_282[k];

        t_175[k] = -ik_103[k]
                   + f_0 * lk_283[k];

        t_176[k] = -ik_104[k]
                   + f_0 * lk_284[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, t_182, ik_105, ik_106, ik_107, \
                         lk_285, lk_286, lk_287, lk_288, lk_289, \
                         lk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -ik_105[k]
                   + f_0 * lk_285[k];

        t_178[k] = -ik_106[k]
                   + f_0 * lk_286[k];

        t_179[k] = -ik_107[k]
                   + f_0 * lk_287[k];

        t_180[k] = f_0 * lk_288[k];

        t_181[k] = f_0 * lk_289[k];

        t_182[k] = f_0 * lk_290[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, t_188, t_189, t_190, lk_291, \
                         lk_292, lk_293, lk_294, lk_295, lk_296, lk_297, \
                         lk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_0 * lk_291[k];

        t_184[k] = f_0 * lk_292[k];

        t_185[k] = f_0 * lk_293[k];

        t_186[k] = f_0 * lk_294[k];

        t_187[k] = f_0 * lk_295[k];

        t_188[k] = f_0 * lk_296[k];

        t_189[k] = f_0 * lk_297[k];

        t_190[k] = f_0 * lk_298[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, t_197, t_198, lk_299, \
                         lk_300, lk_301, lk_302, lk_303, lk_304, lk_305, \
                         lk_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_0 * lk_299[k];

        t_192[k] = f_0 * lk_300[k];

        t_193[k] = f_0 * lk_301[k];

        t_194[k] = f_0 * lk_302[k];

        t_195[k] = f_0 * lk_303[k];

        t_196[k] = f_0 * lk_304[k];

        t_197[k] = f_0 * lk_305[k];

        t_198[k] = f_0 * lk_306[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, t_204, t_205, t_206, lk_307, \
                         lk_308, lk_309, lk_310, lk_311, lk_312, lk_313, \
                         lk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_0 * lk_307[k];

        t_200[k] = f_0 * lk_308[k];

        t_201[k] = f_0 * lk_309[k];

        t_202[k] = f_0 * lk_310[k];

        t_203[k] = f_0 * lk_311[k];

        t_204[k] = f_0 * lk_312[k];

        t_205[k] = f_0 * lk_313[k];

        t_206[k] = f_0 * lk_314[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, t_213, t_214, lk_315, \
                         lk_316, lk_317, lk_318, lk_319, lk_320, lk_321, \
                         lk_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_0 * lk_315[k];

        t_208[k] = f_0 * lk_316[k];

        t_209[k] = f_0 * lk_317[k];

        t_210[k] = f_0 * lk_318[k];

        t_211[k] = f_0 * lk_319[k];

        t_212[k] = f_0 * lk_320[k];

        t_213[k] = f_0 * lk_321[k];

        t_214[k] = f_0 * lk_322[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ik_108, ik_109, ik_110, ik_111, \
                         lk_323, lk_360, lk_361, lk_362, lk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_0 * lk_323[k];

        t_216[k] = -3.0 * ik_108[k]
                   + f_0 * lk_360[k];

        t_217[k] = -3.0 * ik_109[k]
                   + f_0 * lk_361[k];

        t_218[k] = -3.0 * ik_110[k]
                   + f_0 * lk_362[k];

        t_219[k] = -3.0 * ik_111[k]
                   + f_0 * lk_363[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ik_112, ik_113, ik_114, ik_115, \
                         ik_116, lk_364, lk_365, lk_366, lk_367, \
                         lk_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -3.0 * ik_112[k]
                   + f_0 * lk_364[k];

        t_221[k] = -3.0 * ik_113[k]
                   + f_0 * lk_365[k];

        t_222[k] = -3.0 * ik_114[k]
                   + f_0 * lk_366[k];

        t_223[k] = -3.0 * ik_115[k]
                   + f_0 * lk_367[k];

        t_224[k] = -3.0 * ik_116[k]
                   + f_0 * lk_368[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ik_117, ik_118, ik_119, ik_120, \
                         ik_121, lk_369, lk_370, lk_371, lk_372, \
                         lk_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -3.0 * ik_117[k]
                   + f_0 * lk_369[k];

        t_226[k] = -3.0 * ik_118[k]
                   + f_0 * lk_370[k];

        t_227[k] = -3.0 * ik_119[k]
                   + f_0 * lk_371[k];

        t_228[k] = -3.0 * ik_120[k]
                   + f_0 * lk_372[k];

        t_229[k] = -3.0 * ik_121[k]
                   + f_0 * lk_373[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ik_122, ik_123, ik_124, ik_125, \
                         ik_126, lk_374, lk_375, lk_376, lk_377, \
                         lk_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -3.0 * ik_122[k]
                   + f_0 * lk_374[k];

        t_231[k] = -3.0 * ik_123[k]
                   + f_0 * lk_375[k];

        t_232[k] = -3.0 * ik_124[k]
                   + f_0 * lk_376[k];

        t_233[k] = -3.0 * ik_125[k]
                   + f_0 * lk_377[k];

        t_234[k] = -3.0 * ik_126[k]
                   + f_0 * lk_378[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ik_127, ik_128, ik_129, ik_130, \
                         ik_131, lk_379, lk_380, lk_381, lk_382, \
                         lk_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -3.0 * ik_127[k]
                   + f_0 * lk_379[k];

        t_236[k] = -3.0 * ik_128[k]
                   + f_0 * lk_380[k];

        t_237[k] = -3.0 * ik_129[k]
                   + f_0 * lk_381[k];

        t_238[k] = -3.0 * ik_130[k]
                   + f_0 * lk_382[k];

        t_239[k] = -3.0 * ik_131[k]
                   + f_0 * lk_383[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ik_132, ik_133, ik_134, ik_135, \
                         ik_136, lk_384, lk_385, lk_386, lk_387, \
                         lk_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * ik_132[k]
                   + f_0 * lk_384[k];

        t_241[k] = -3.0 * ik_133[k]
                   + f_0 * lk_385[k];

        t_242[k] = -3.0 * ik_134[k]
                   + f_0 * lk_386[k];

        t_243[k] = -3.0 * ik_135[k]
                   + f_0 * lk_387[k];

        t_244[k] = -3.0 * ik_136[k]
                   + f_0 * lk_388[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ik_137, ik_138, ik_139, ik_140, \
                         ik_141, lk_389, lk_390, lk_391, lk_392, \
                         lk_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * ik_137[k]
                   + f_0 * lk_389[k];

        t_246[k] = -3.0 * ik_138[k]
                   + f_0 * lk_390[k];

        t_247[k] = -3.0 * ik_139[k]
                   + f_0 * lk_391[k];

        t_248[k] = -3.0 * ik_140[k]
                   + f_0 * lk_392[k];

        t_249[k] = -3.0 * ik_141[k]
                   + f_0 * lk_393[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ik_142, ik_143, ik_144, ik_145, \
                         ik_146, lk_394, lk_395, lk_396, lk_397, \
                         lk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -3.0 * ik_142[k]
                   + f_0 * lk_394[k];

        t_251[k] = -3.0 * ik_143[k]
                   + f_0 * lk_395[k];

        t_252[k] = -2.0 * ik_144[k]
                   + f_0 * lk_396[k];

        t_253[k] = -2.0 * ik_145[k]
                   + f_0 * lk_397[k];

        t_254[k] = -2.0 * ik_146[k]
                   + f_0 * lk_398[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ik_147, ik_148, ik_149, ik_150, \
                         ik_151, lk_399, lk_400, lk_401, lk_402, \
                         lk_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * ik_147[k]
                   + f_0 * lk_399[k];

        t_256[k] = -2.0 * ik_148[k]
                   + f_0 * lk_400[k];

        t_257[k] = -2.0 * ik_149[k]
                   + f_0 * lk_401[k];

        t_258[k] = -2.0 * ik_150[k]
                   + f_0 * lk_402[k];

        t_259[k] = -2.0 * ik_151[k]
                   + f_0 * lk_403[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ik_152, ik_153, ik_154, ik_155, \
                         ik_156, lk_404, lk_405, lk_406, lk_407, \
                         lk_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * ik_152[k]
                   + f_0 * lk_404[k];

        t_261[k] = -2.0 * ik_153[k]
                   + f_0 * lk_405[k];

        t_262[k] = -2.0 * ik_154[k]
                   + f_0 * lk_406[k];

        t_263[k] = -2.0 * ik_155[k]
                   + f_0 * lk_407[k];

        t_264[k] = -2.0 * ik_156[k]
                   + f_0 * lk_408[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ik_157, ik_158, ik_159, ik_160, \
                         ik_161, lk_409, lk_410, lk_411, lk_412, \
                         lk_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * ik_157[k]
                   + f_0 * lk_409[k];

        t_266[k] = -2.0 * ik_158[k]
                   + f_0 * lk_410[k];

        t_267[k] = -2.0 * ik_159[k]
                   + f_0 * lk_411[k];

        t_268[k] = -2.0 * ik_160[k]
                   + f_0 * lk_412[k];

        t_269[k] = -2.0 * ik_161[k]
                   + f_0 * lk_413[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ik_162, ik_163, ik_164, ik_165, \
                         ik_166, lk_414, lk_415, lk_416, lk_417, \
                         lk_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * ik_162[k]
                   + f_0 * lk_414[k];

        t_271[k] = -2.0 * ik_163[k]
                   + f_0 * lk_415[k];

        t_272[k] = -2.0 * ik_164[k]
                   + f_0 * lk_416[k];

        t_273[k] = -2.0 * ik_165[k]
                   + f_0 * lk_417[k];

        t_274[k] = -2.0 * ik_166[k]
                   + f_0 * lk_418[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ik_167, ik_168, ik_169, ik_170, \
                         ik_171, lk_419, lk_420, lk_421, lk_422, \
                         lk_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * ik_167[k]
                   + f_0 * lk_419[k];

        t_276[k] = -2.0 * ik_168[k]
                   + f_0 * lk_420[k];

        t_277[k] = -2.0 * ik_169[k]
                   + f_0 * lk_421[k];

        t_278[k] = -2.0 * ik_170[k]
                   + f_0 * lk_422[k];

        t_279[k] = -2.0 * ik_171[k]
                   + f_0 * lk_423[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ik_172, ik_173, ik_174, ik_175, \
                         ik_176, lk_424, lk_425, lk_426, lk_427, \
                         lk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * ik_172[k]
                   + f_0 * lk_424[k];

        t_281[k] = -2.0 * ik_173[k]
                   + f_0 * lk_425[k];

        t_282[k] = -2.0 * ik_174[k]
                   + f_0 * lk_426[k];

        t_283[k] = -2.0 * ik_175[k]
                   + f_0 * lk_427[k];

        t_284[k] = -2.0 * ik_176[k]
                   + f_0 * lk_428[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ik_177, ik_178, ik_179, ik_180, \
                         ik_181, lk_429, lk_430, lk_431, lk_432, \
                         lk_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -2.0 * ik_177[k]
                   + f_0 * lk_429[k];

        t_286[k] = -2.0 * ik_178[k]
                   + f_0 * lk_430[k];

        t_287[k] = -2.0 * ik_179[k]
                   + f_0 * lk_431[k];

        t_288[k] = -ik_180[k]
                   + f_0 * lk_432[k];

        t_289[k] = -ik_181[k]
                   + f_0 * lk_433[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ik_182, ik_183, ik_184, ik_185, \
                         ik_186, lk_434, lk_435, lk_436, lk_437, \
                         lk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -ik_182[k]
                   + f_0 * lk_434[k];

        t_291[k] = -ik_183[k]
                   + f_0 * lk_435[k];

        t_292[k] = -ik_184[k]
                   + f_0 * lk_436[k];

        t_293[k] = -ik_185[k]
                   + f_0 * lk_437[k];

        t_294[k] = -ik_186[k]
                   + f_0 * lk_438[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ik_187, ik_188, ik_189, ik_190, \
                         ik_191, lk_439, lk_440, lk_441, lk_442, \
                         lk_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -ik_187[k]
                   + f_0 * lk_439[k];

        t_296[k] = -ik_188[k]
                   + f_0 * lk_440[k];

        t_297[k] = -ik_189[k]
                   + f_0 * lk_441[k];

        t_298[k] = -ik_190[k]
                   + f_0 * lk_442[k];

        t_299[k] = -ik_191[k]
                   + f_0 * lk_443[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ik_192, ik_193, ik_194, ik_195, \
                         ik_196, lk_444, lk_445, lk_446, lk_447, \
                         lk_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -ik_192[k]
                   + f_0 * lk_444[k];

        t_301[k] = -ik_193[k]
                   + f_0 * lk_445[k];

        t_302[k] = -ik_194[k]
                   + f_0 * lk_446[k];

        t_303[k] = -ik_195[k]
                   + f_0 * lk_447[k];

        t_304[k] = -ik_196[k]
                   + f_0 * lk_448[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ik_197, ik_198, ik_199, ik_200, \
                         ik_201, lk_449, lk_450, lk_451, lk_452, \
                         lk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -ik_197[k]
                   + f_0 * lk_449[k];

        t_306[k] = -ik_198[k]
                   + f_0 * lk_450[k];

        t_307[k] = -ik_199[k]
                   + f_0 * lk_451[k];

        t_308[k] = -ik_200[k]
                   + f_0 * lk_452[k];

        t_309[k] = -ik_201[k]
                   + f_0 * lk_453[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ik_202, ik_203, ik_204, ik_205, \
                         ik_206, lk_454, lk_455, lk_456, lk_457, \
                         lk_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -ik_202[k]
                   + f_0 * lk_454[k];

        t_311[k] = -ik_203[k]
                   + f_0 * lk_455[k];

        t_312[k] = -ik_204[k]
                   + f_0 * lk_456[k];

        t_313[k] = -ik_205[k]
                   + f_0 * lk_457[k];

        t_314[k] = -ik_206[k]
                   + f_0 * lk_458[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ik_207, ik_208, ik_209, ik_210, \
                         ik_211, lk_459, lk_460, lk_461, lk_462, \
                         lk_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -ik_207[k]
                   + f_0 * lk_459[k];

        t_316[k] = -ik_208[k]
                   + f_0 * lk_460[k];

        t_317[k] = -ik_209[k]
                   + f_0 * lk_461[k];

        t_318[k] = -ik_210[k]
                   + f_0 * lk_462[k];

        t_319[k] = -ik_211[k]
                   + f_0 * lk_463[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, t_325, ik_212, ik_213, ik_214, \
                         ik_215, lk_464, lk_465, lk_466, lk_467, lk_468, \
                         lk_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -ik_212[k]
                   + f_0 * lk_464[k];

        t_321[k] = -ik_213[k]
                   + f_0 * lk_465[k];

        t_322[k] = -ik_214[k]
                   + f_0 * lk_466[k];

        t_323[k] = -ik_215[k]
                   + f_0 * lk_467[k];

        t_324[k] = f_0 * lk_468[k];

        t_325[k] = f_0 * lk_469[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, t_332, t_333, lk_470, \
                         lk_471, lk_472, lk_473, lk_474, lk_475, lk_476, \
                         lk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_0 * lk_470[k];

        t_327[k] = f_0 * lk_471[k];

        t_328[k] = f_0 * lk_472[k];

        t_329[k] = f_0 * lk_473[k];

        t_330[k] = f_0 * lk_474[k];

        t_331[k] = f_0 * lk_475[k];

        t_332[k] = f_0 * lk_476[k];

        t_333[k] = f_0 * lk_477[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
                                                    const size_t ncols,
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

    const auto *ik_216 = buffer.data(ik + 216);
    const auto *ik_217 = buffer.data(ik + 217);
    const auto *ik_218 = buffer.data(ik + 218);
    const auto *ik_219 = buffer.data(ik + 219);
    const auto *ik_220 = buffer.data(ik + 220);
    const auto *ik_221 = buffer.data(ik + 221);
    const auto *ik_222 = buffer.data(ik + 222);
    const auto *ik_223 = buffer.data(ik + 223);
    const auto *ik_224 = buffer.data(ik + 224);
    const auto *ik_225 = buffer.data(ik + 225);
    const auto *ik_226 = buffer.data(ik + 226);
    const auto *ik_227 = buffer.data(ik + 227);
    const auto *ik_228 = buffer.data(ik + 228);
    const auto *ik_229 = buffer.data(ik + 229);
    const auto *ik_230 = buffer.data(ik + 230);
    const auto *ik_231 = buffer.data(ik + 231);
    const auto *ik_232 = buffer.data(ik + 232);
    const auto *ik_233 = buffer.data(ik + 233);
    const auto *ik_234 = buffer.data(ik + 234);
    const auto *ik_235 = buffer.data(ik + 235);
    const auto *ik_236 = buffer.data(ik + 236);
    const auto *ik_237 = buffer.data(ik + 237);
    const auto *ik_238 = buffer.data(ik + 238);
    const auto *ik_239 = buffer.data(ik + 239);
    const auto *ik_240 = buffer.data(ik + 240);
    const auto *ik_241 = buffer.data(ik + 241);
    const auto *ik_242 = buffer.data(ik + 242);
    const auto *ik_243 = buffer.data(ik + 243);
    const auto *ik_244 = buffer.data(ik + 244);
    const auto *ik_245 = buffer.data(ik + 245);
    const auto *ik_246 = buffer.data(ik + 246);
    const auto *ik_247 = buffer.data(ik + 247);
    const auto *ik_248 = buffer.data(ik + 248);
    const auto *ik_249 = buffer.data(ik + 249);
    const auto *ik_250 = buffer.data(ik + 250);
    const auto *ik_251 = buffer.data(ik + 251);
    const auto *ik_252 = buffer.data(ik + 252);
    const auto *ik_253 = buffer.data(ik + 253);
    const auto *ik_254 = buffer.data(ik + 254);
    const auto *ik_255 = buffer.data(ik + 255);
    const auto *ik_256 = buffer.data(ik + 256);
    const auto *ik_257 = buffer.data(ik + 257);
    const auto *ik_258 = buffer.data(ik + 258);
    const auto *ik_259 = buffer.data(ik + 259);
    const auto *ik_260 = buffer.data(ik + 260);
    const auto *ik_261 = buffer.data(ik + 261);
    const auto *ik_262 = buffer.data(ik + 262);
    const auto *ik_263 = buffer.data(ik + 263);
    const auto *ik_264 = buffer.data(ik + 264);
    const auto *ik_265 = buffer.data(ik + 265);
    const auto *ik_266 = buffer.data(ik + 266);
    const auto *ik_267 = buffer.data(ik + 267);
    const auto *ik_268 = buffer.data(ik + 268);
    const auto *ik_269 = buffer.data(ik + 269);
    const auto *ik_270 = buffer.data(ik + 270);
    const auto *ik_271 = buffer.data(ik + 271);
    const auto *ik_272 = buffer.data(ik + 272);
    const auto *ik_273 = buffer.data(ik + 273);
    const auto *ik_274 = buffer.data(ik + 274);
    const auto *ik_275 = buffer.data(ik + 275);
    const auto *ik_276 = buffer.data(ik + 276);
    const auto *ik_277 = buffer.data(ik + 277);
    const auto *ik_278 = buffer.data(ik + 278);
    const auto *ik_279 = buffer.data(ik + 279);
    const auto *ik_280 = buffer.data(ik + 280);
    const auto *ik_281 = buffer.data(ik + 281);
    const auto *ik_282 = buffer.data(ik + 282);
    const auto *ik_283 = buffer.data(ik + 283);
    const auto *ik_284 = buffer.data(ik + 284);
    const auto *ik_285 = buffer.data(ik + 285);
    const auto *ik_286 = buffer.data(ik + 286);
    const auto *ik_287 = buffer.data(ik + 287);
    const auto *ik_288 = buffer.data(ik + 288);
    const auto *ik_289 = buffer.data(ik + 289);
    const auto *ik_290 = buffer.data(ik + 290);
    const auto *ik_291 = buffer.data(ik + 291);
    const auto *ik_292 = buffer.data(ik + 292);
    const auto *ik_293 = buffer.data(ik + 293);
    const auto *ik_294 = buffer.data(ik + 294);
    const auto *ik_295 = buffer.data(ik + 295);
    const auto *ik_296 = buffer.data(ik + 296);
    const auto *ik_297 = buffer.data(ik + 297);
    const auto *ik_298 = buffer.data(ik + 298);
    const auto *ik_299 = buffer.data(ik + 299);
    const auto *ik_300 = buffer.data(ik + 300);
    const auto *ik_301 = buffer.data(ik + 301);
    const auto *ik_302 = buffer.data(ik + 302);
    const auto *ik_303 = buffer.data(ik + 303);
    const auto *ik_304 = buffer.data(ik + 304);
    const auto *ik_305 = buffer.data(ik + 305);
    const auto *ik_306 = buffer.data(ik + 306);
    const auto *ik_307 = buffer.data(ik + 307);
    const auto *ik_308 = buffer.data(ik + 308);
    const auto *ik_309 = buffer.data(ik + 309);
    const auto *ik_310 = buffer.data(ik + 310);
    const auto *ik_311 = buffer.data(ik + 311);
    const auto *ik_312 = buffer.data(ik + 312);
    const auto *ik_313 = buffer.data(ik + 313);
    const auto *ik_314 = buffer.data(ik + 314);
    const auto *ik_315 = buffer.data(ik + 315);
    const auto *ik_316 = buffer.data(ik + 316);
    const auto *ik_317 = buffer.data(ik + 317);
    const auto *ik_318 = buffer.data(ik + 318);
    const auto *ik_319 = buffer.data(ik + 319);
    const auto *ik_320 = buffer.data(ik + 320);
    const auto *ik_321 = buffer.data(ik + 321);
    const auto *ik_322 = buffer.data(ik + 322);
    const auto *ik_323 = buffer.data(ik + 323);
    const auto *ik_324 = buffer.data(ik + 324);
    const auto *ik_325 = buffer.data(ik + 325);
    const auto *ik_326 = buffer.data(ik + 326);
    const auto *ik_327 = buffer.data(ik + 327);
    const auto *ik_328 = buffer.data(ik + 328);
    const auto *ik_329 = buffer.data(ik + 329);
    const auto *ik_330 = buffer.data(ik + 330);
    const auto *ik_331 = buffer.data(ik + 331);
    const auto *ik_332 = buffer.data(ik + 332);
    const auto *ik_333 = buffer.data(ik + 333);
    const auto *ik_334 = buffer.data(ik + 334);
    const auto *ik_335 = buffer.data(ik + 335);
    const auto *ik_336 = buffer.data(ik + 336);
    const auto *ik_337 = buffer.data(ik + 337);
    const auto *ik_338 = buffer.data(ik + 338);
    const auto *ik_339 = buffer.data(ik + 339);
    const auto *ik_340 = buffer.data(ik + 340);
    const auto *ik_341 = buffer.data(ik + 341);
    const auto *ik_342 = buffer.data(ik + 342);
    const auto *ik_343 = buffer.data(ik + 343);
    const auto *ik_344 = buffer.data(ik + 344);
    const auto *ik_345 = buffer.data(ik + 345);
    const auto *ik_346 = buffer.data(ik + 346);
    const auto *ik_347 = buffer.data(ik + 347);
    const auto *ik_348 = buffer.data(ik + 348);
    const auto *ik_349 = buffer.data(ik + 349);

    const auto *lk_478 = buffer.data(lk + 478);
    const auto *lk_479 = buffer.data(lk + 479);
    const auto *lk_480 = buffer.data(lk + 480);
    const auto *lk_481 = buffer.data(lk + 481);
    const auto *lk_482 = buffer.data(lk + 482);
    const auto *lk_483 = buffer.data(lk + 483);
    const auto *lk_484 = buffer.data(lk + 484);
    const auto *lk_485 = buffer.data(lk + 485);
    const auto *lk_486 = buffer.data(lk + 486);
    const auto *lk_487 = buffer.data(lk + 487);
    const auto *lk_488 = buffer.data(lk + 488);
    const auto *lk_489 = buffer.data(lk + 489);
    const auto *lk_490 = buffer.data(lk + 490);
    const auto *lk_491 = buffer.data(lk + 491);
    const auto *lk_492 = buffer.data(lk + 492);
    const auto *lk_493 = buffer.data(lk + 493);
    const auto *lk_494 = buffer.data(lk + 494);
    const auto *lk_495 = buffer.data(lk + 495);
    const auto *lk_496 = buffer.data(lk + 496);
    const auto *lk_497 = buffer.data(lk + 497);
    const auto *lk_498 = buffer.data(lk + 498);
    const auto *lk_499 = buffer.data(lk + 499);
    const auto *lk_500 = buffer.data(lk + 500);
    const auto *lk_501 = buffer.data(lk + 501);
    const auto *lk_502 = buffer.data(lk + 502);
    const auto *lk_503 = buffer.data(lk + 503);
    const auto *lk_540 = buffer.data(lk + 540);
    const auto *lk_541 = buffer.data(lk + 541);
    const auto *lk_542 = buffer.data(lk + 542);
    const auto *lk_543 = buffer.data(lk + 543);
    const auto *lk_544 = buffer.data(lk + 544);
    const auto *lk_545 = buffer.data(lk + 545);
    const auto *lk_546 = buffer.data(lk + 546);
    const auto *lk_547 = buffer.data(lk + 547);
    const auto *lk_548 = buffer.data(lk + 548);
    const auto *lk_549 = buffer.data(lk + 549);
    const auto *lk_550 = buffer.data(lk + 550);
    const auto *lk_551 = buffer.data(lk + 551);
    const auto *lk_552 = buffer.data(lk + 552);
    const auto *lk_553 = buffer.data(lk + 553);
    const auto *lk_554 = buffer.data(lk + 554);
    const auto *lk_555 = buffer.data(lk + 555);
    const auto *lk_556 = buffer.data(lk + 556);
    const auto *lk_557 = buffer.data(lk + 557);
    const auto *lk_558 = buffer.data(lk + 558);
    const auto *lk_559 = buffer.data(lk + 559);
    const auto *lk_560 = buffer.data(lk + 560);
    const auto *lk_561 = buffer.data(lk + 561);
    const auto *lk_562 = buffer.data(lk + 562);
    const auto *lk_563 = buffer.data(lk + 563);
    const auto *lk_564 = buffer.data(lk + 564);
    const auto *lk_565 = buffer.data(lk + 565);
    const auto *lk_566 = buffer.data(lk + 566);
    const auto *lk_567 = buffer.data(lk + 567);
    const auto *lk_568 = buffer.data(lk + 568);
    const auto *lk_569 = buffer.data(lk + 569);
    const auto *lk_570 = buffer.data(lk + 570);
    const auto *lk_571 = buffer.data(lk + 571);
    const auto *lk_572 = buffer.data(lk + 572);
    const auto *lk_573 = buffer.data(lk + 573);
    const auto *lk_574 = buffer.data(lk + 574);
    const auto *lk_575 = buffer.data(lk + 575);
    const auto *lk_576 = buffer.data(lk + 576);
    const auto *lk_577 = buffer.data(lk + 577);
    const auto *lk_578 = buffer.data(lk + 578);
    const auto *lk_579 = buffer.data(lk + 579);
    const auto *lk_580 = buffer.data(lk + 580);
    const auto *lk_581 = buffer.data(lk + 581);
    const auto *lk_582 = buffer.data(lk + 582);
    const auto *lk_583 = buffer.data(lk + 583);
    const auto *lk_584 = buffer.data(lk + 584);
    const auto *lk_585 = buffer.data(lk + 585);
    const auto *lk_586 = buffer.data(lk + 586);
    const auto *lk_587 = buffer.data(lk + 587);
    const auto *lk_588 = buffer.data(lk + 588);
    const auto *lk_589 = buffer.data(lk + 589);
    const auto *lk_590 = buffer.data(lk + 590);
    const auto *lk_591 = buffer.data(lk + 591);
    const auto *lk_592 = buffer.data(lk + 592);
    const auto *lk_593 = buffer.data(lk + 593);
    const auto *lk_594 = buffer.data(lk + 594);
    const auto *lk_595 = buffer.data(lk + 595);
    const auto *lk_596 = buffer.data(lk + 596);
    const auto *lk_597 = buffer.data(lk + 597);
    const auto *lk_598 = buffer.data(lk + 598);
    const auto *lk_599 = buffer.data(lk + 599);
    const auto *lk_600 = buffer.data(lk + 600);
    const auto *lk_601 = buffer.data(lk + 601);
    const auto *lk_602 = buffer.data(lk + 602);
    const auto *lk_603 = buffer.data(lk + 603);
    const auto *lk_604 = buffer.data(lk + 604);
    const auto *lk_605 = buffer.data(lk + 605);
    const auto *lk_606 = buffer.data(lk + 606);
    const auto *lk_607 = buffer.data(lk + 607);
    const auto *lk_608 = buffer.data(lk + 608);
    const auto *lk_609 = buffer.data(lk + 609);
    const auto *lk_610 = buffer.data(lk + 610);
    const auto *lk_611 = buffer.data(lk + 611);
    const auto *lk_612 = buffer.data(lk + 612);
    const auto *lk_613 = buffer.data(lk + 613);
    const auto *lk_614 = buffer.data(lk + 614);
    const auto *lk_615 = buffer.data(lk + 615);
    const auto *lk_616 = buffer.data(lk + 616);
    const auto *lk_617 = buffer.data(lk + 617);
    const auto *lk_618 = buffer.data(lk + 618);
    const auto *lk_619 = buffer.data(lk + 619);
    const auto *lk_620 = buffer.data(lk + 620);
    const auto *lk_621 = buffer.data(lk + 621);
    const auto *lk_622 = buffer.data(lk + 622);
    const auto *lk_623 = buffer.data(lk + 623);
    const auto *lk_624 = buffer.data(lk + 624);
    const auto *lk_625 = buffer.data(lk + 625);
    const auto *lk_626 = buffer.data(lk + 626);
    const auto *lk_627 = buffer.data(lk + 627);
    const auto *lk_628 = buffer.data(lk + 628);
    const auto *lk_629 = buffer.data(lk + 629);
    const auto *lk_630 = buffer.data(lk + 630);
    const auto *lk_631 = buffer.data(lk + 631);
    const auto *lk_632 = buffer.data(lk + 632);
    const auto *lk_633 = buffer.data(lk + 633);
    const auto *lk_634 = buffer.data(lk + 634);
    const auto *lk_635 = buffer.data(lk + 635);
    const auto *lk_636 = buffer.data(lk + 636);
    const auto *lk_637 = buffer.data(lk + 637);
    const auto *lk_638 = buffer.data(lk + 638);
    const auto *lk_639 = buffer.data(lk + 639);
    const auto *lk_640 = buffer.data(lk + 640);
    const auto *lk_641 = buffer.data(lk + 641);
    const auto *lk_642 = buffer.data(lk + 642);
    const auto *lk_643 = buffer.data(lk + 643);
    const auto *lk_644 = buffer.data(lk + 644);
    const auto *lk_645 = buffer.data(lk + 645);
    const auto *lk_646 = buffer.data(lk + 646);
    const auto *lk_647 = buffer.data(lk + 647);
    const auto *lk_648 = buffer.data(lk + 648);
    const auto *lk_649 = buffer.data(lk + 649);
    const auto *lk_650 = buffer.data(lk + 650);
    const auto *lk_651 = buffer.data(lk + 651);
    const auto *lk_652 = buffer.data(lk + 652);
    const auto *lk_653 = buffer.data(lk + 653);
    const auto *lk_654 = buffer.data(lk + 654);
    const auto *lk_655 = buffer.data(lk + 655);
    const auto *lk_656 = buffer.data(lk + 656);
    const auto *lk_657 = buffer.data(lk + 657);
    const auto *lk_658 = buffer.data(lk + 658);
    const auto *lk_659 = buffer.data(lk + 659);
    const auto *lk_660 = buffer.data(lk + 660);
    const auto *lk_661 = buffer.data(lk + 661);
    const auto *lk_662 = buffer.data(lk + 662);
    const auto *lk_663 = buffer.data(lk + 663);
    const auto *lk_664 = buffer.data(lk + 664);
    const auto *lk_665 = buffer.data(lk + 665);
    const auto *lk_666 = buffer.data(lk + 666);
    const auto *lk_667 = buffer.data(lk + 667);
    const auto *lk_668 = buffer.data(lk + 668);
    const auto *lk_669 = buffer.data(lk + 669);
    const auto *lk_670 = buffer.data(lk + 670);
    const auto *lk_671 = buffer.data(lk + 671);
    const auto *lk_672 = buffer.data(lk + 672);
    const auto *lk_673 = buffer.data(lk + 673);

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, t_339, t_340, t_341, lk_478, \
                         lk_479, lk_480, lk_481, lk_482, lk_483, lk_484, \
                         lk_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_0 * lk_478[k];

        t_335[k] = f_0 * lk_479[k];

        t_336[k] = f_0 * lk_480[k];

        t_337[k] = f_0 * lk_481[k];

        t_338[k] = f_0 * lk_482[k];

        t_339[k] = f_0 * lk_483[k];

        t_340[k] = f_0 * lk_484[k];

        t_341[k] = f_0 * lk_485[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, t_347, t_348, t_349, lk_486, \
                         lk_487, lk_488, lk_489, lk_490, lk_491, lk_492, \
                         lk_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_0 * lk_486[k];

        t_343[k] = f_0 * lk_487[k];

        t_344[k] = f_0 * lk_488[k];

        t_345[k] = f_0 * lk_489[k];

        t_346[k] = f_0 * lk_490[k];

        t_347[k] = f_0 * lk_491[k];

        t_348[k] = f_0 * lk_492[k];

        t_349[k] = f_0 * lk_493[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, t_355, t_356, t_357, lk_494, \
                         lk_495, lk_496, lk_497, lk_498, lk_499, lk_500, \
                         lk_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_0 * lk_494[k];

        t_351[k] = f_0 * lk_495[k];

        t_352[k] = f_0 * lk_496[k];

        t_353[k] = f_0 * lk_497[k];

        t_354[k] = f_0 * lk_498[k];

        t_355[k] = f_0 * lk_499[k];

        t_356[k] = f_0 * lk_500[k];

        t_357[k] = f_0 * lk_501[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, t_363, ik_216, ik_217, ik_218, \
                         ik_219, lk_502, lk_503, lk_540, lk_541, lk_542, \
                         lk_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_0 * lk_502[k];

        t_359[k] = f_0 * lk_503[k];

        t_360[k] = -4.0 * ik_216[k]
                   + f_0 * lk_540[k];

        t_361[k] = -4.0 * ik_217[k]
                   + f_0 * lk_541[k];

        t_362[k] = -4.0 * ik_218[k]
                   + f_0 * lk_542[k];

        t_363[k] = -4.0 * ik_219[k]
                   + f_0 * lk_543[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, t_368, ik_220, ik_221, ik_222, ik_223, \
                         ik_224, lk_544, lk_545, lk_546, lk_547, \
                         lk_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = -4.0 * ik_220[k]
                   + f_0 * lk_544[k];

        t_365[k] = -4.0 * ik_221[k]
                   + f_0 * lk_545[k];

        t_366[k] = -4.0 * ik_222[k]
                   + f_0 * lk_546[k];

        t_367[k] = -4.0 * ik_223[k]
                   + f_0 * lk_547[k];

        t_368[k] = -4.0 * ik_224[k]
                   + f_0 * lk_548[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, t_373, ik_225, ik_226, ik_227, ik_228, \
                         ik_229, lk_549, lk_550, lk_551, lk_552, \
                         lk_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = -4.0 * ik_225[k]
                   + f_0 * lk_549[k];

        t_370[k] = -4.0 * ik_226[k]
                   + f_0 * lk_550[k];

        t_371[k] = -4.0 * ik_227[k]
                   + f_0 * lk_551[k];

        t_372[k] = -4.0 * ik_228[k]
                   + f_0 * lk_552[k];

        t_373[k] = -4.0 * ik_229[k]
                   + f_0 * lk_553[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, t_378, ik_230, ik_231, ik_232, ik_233, \
                         ik_234, lk_554, lk_555, lk_556, lk_557, \
                         lk_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = -4.0 * ik_230[k]
                   + f_0 * lk_554[k];

        t_375[k] = -4.0 * ik_231[k]
                   + f_0 * lk_555[k];

        t_376[k] = -4.0 * ik_232[k]
                   + f_0 * lk_556[k];

        t_377[k] = -4.0 * ik_233[k]
                   + f_0 * lk_557[k];

        t_378[k] = -4.0 * ik_234[k]
                   + f_0 * lk_558[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, ik_235, ik_236, ik_237, ik_238, \
                         ik_239, lk_559, lk_560, lk_561, lk_562, \
                         lk_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = -4.0 * ik_235[k]
                   + f_0 * lk_559[k];

        t_380[k] = -4.0 * ik_236[k]
                   + f_0 * lk_560[k];

        t_381[k] = -4.0 * ik_237[k]
                   + f_0 * lk_561[k];

        t_382[k] = -4.0 * ik_238[k]
                   + f_0 * lk_562[k];

        t_383[k] = -4.0 * ik_239[k]
                   + f_0 * lk_563[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, t_388, ik_240, ik_241, ik_242, ik_243, \
                         ik_244, lk_564, lk_565, lk_566, lk_567, \
                         lk_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = -4.0 * ik_240[k]
                   + f_0 * lk_564[k];

        t_385[k] = -4.0 * ik_241[k]
                   + f_0 * lk_565[k];

        t_386[k] = -4.0 * ik_242[k]
                   + f_0 * lk_566[k];

        t_387[k] = -4.0 * ik_243[k]
                   + f_0 * lk_567[k];

        t_388[k] = -4.0 * ik_244[k]
                   + f_0 * lk_568[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, ik_245, ik_246, ik_247, ik_248, \
                         ik_249, lk_569, lk_570, lk_571, lk_572, \
                         lk_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = -4.0 * ik_245[k]
                   + f_0 * lk_569[k];

        t_390[k] = -4.0 * ik_246[k]
                   + f_0 * lk_570[k];

        t_391[k] = -4.0 * ik_247[k]
                   + f_0 * lk_571[k];

        t_392[k] = -4.0 * ik_248[k]
                   + f_0 * lk_572[k];

        t_393[k] = -4.0 * ik_249[k]
                   + f_0 * lk_573[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, ik_250, ik_251, ik_252, ik_253, \
                         ik_254, lk_574, lk_575, lk_576, lk_577, \
                         lk_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = -4.0 * ik_250[k]
                   + f_0 * lk_574[k];

        t_395[k] = -4.0 * ik_251[k]
                   + f_0 * lk_575[k];

        t_396[k] = -3.0 * ik_252[k]
                   + f_0 * lk_576[k];

        t_397[k] = -3.0 * ik_253[k]
                   + f_0 * lk_577[k];

        t_398[k] = -3.0 * ik_254[k]
                   + f_0 * lk_578[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, t_403, ik_255, ik_256, ik_257, ik_258, \
                         ik_259, lk_579, lk_580, lk_581, lk_582, \
                         lk_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = -3.0 * ik_255[k]
                   + f_0 * lk_579[k];

        t_400[k] = -3.0 * ik_256[k]
                   + f_0 * lk_580[k];

        t_401[k] = -3.0 * ik_257[k]
                   + f_0 * lk_581[k];

        t_402[k] = -3.0 * ik_258[k]
                   + f_0 * lk_582[k];

        t_403[k] = -3.0 * ik_259[k]
                   + f_0 * lk_583[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, t_408, ik_260, ik_261, ik_262, ik_263, \
                         ik_264, lk_584, lk_585, lk_586, lk_587, \
                         lk_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = -3.0 * ik_260[k]
                   + f_0 * lk_584[k];

        t_405[k] = -3.0 * ik_261[k]
                   + f_0 * lk_585[k];

        t_406[k] = -3.0 * ik_262[k]
                   + f_0 * lk_586[k];

        t_407[k] = -3.0 * ik_263[k]
                   + f_0 * lk_587[k];

        t_408[k] = -3.0 * ik_264[k]
                   + f_0 * lk_588[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, ik_265, ik_266, ik_267, ik_268, \
                         ik_269, lk_589, lk_590, lk_591, lk_592, \
                         lk_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = -3.0 * ik_265[k]
                   + f_0 * lk_589[k];

        t_410[k] = -3.0 * ik_266[k]
                   + f_0 * lk_590[k];

        t_411[k] = -3.0 * ik_267[k]
                   + f_0 * lk_591[k];

        t_412[k] = -3.0 * ik_268[k]
                   + f_0 * lk_592[k];

        t_413[k] = -3.0 * ik_269[k]
                   + f_0 * lk_593[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, ik_270, ik_271, ik_272, ik_273, \
                         ik_274, lk_594, lk_595, lk_596, lk_597, \
                         lk_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = -3.0 * ik_270[k]
                   + f_0 * lk_594[k];

        t_415[k] = -3.0 * ik_271[k]
                   + f_0 * lk_595[k];

        t_416[k] = -3.0 * ik_272[k]
                   + f_0 * lk_596[k];

        t_417[k] = -3.0 * ik_273[k]
                   + f_0 * lk_597[k];

        t_418[k] = -3.0 * ik_274[k]
                   + f_0 * lk_598[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, t_423, ik_275, ik_276, ik_277, ik_278, \
                         ik_279, lk_599, lk_600, lk_601, lk_602, \
                         lk_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = -3.0 * ik_275[k]
                   + f_0 * lk_599[k];

        t_420[k] = -3.0 * ik_276[k]
                   + f_0 * lk_600[k];

        t_421[k] = -3.0 * ik_277[k]
                   + f_0 * lk_601[k];

        t_422[k] = -3.0 * ik_278[k]
                   + f_0 * lk_602[k];

        t_423[k] = -3.0 * ik_279[k]
                   + f_0 * lk_603[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, t_428, ik_280, ik_281, ik_282, ik_283, \
                         ik_284, lk_604, lk_605, lk_606, lk_607, \
                         lk_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = -3.0 * ik_280[k]
                   + f_0 * lk_604[k];

        t_425[k] = -3.0 * ik_281[k]
                   + f_0 * lk_605[k];

        t_426[k] = -3.0 * ik_282[k]
                   + f_0 * lk_606[k];

        t_427[k] = -3.0 * ik_283[k]
                   + f_0 * lk_607[k];

        t_428[k] = -3.0 * ik_284[k]
                   + f_0 * lk_608[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, ik_285, ik_286, ik_287, ik_288, \
                         ik_289, lk_609, lk_610, lk_611, lk_612, \
                         lk_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = -3.0 * ik_285[k]
                   + f_0 * lk_609[k];

        t_430[k] = -3.0 * ik_286[k]
                   + f_0 * lk_610[k];

        t_431[k] = -3.0 * ik_287[k]
                   + f_0 * lk_611[k];

        t_432[k] = -2.0 * ik_288[k]
                   + f_0 * lk_612[k];

        t_433[k] = -2.0 * ik_289[k]
                   + f_0 * lk_613[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, t_438, ik_290, ik_291, ik_292, ik_293, \
                         ik_294, lk_614, lk_615, lk_616, lk_617, \
                         lk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = -2.0 * ik_290[k]
                   + f_0 * lk_614[k];

        t_435[k] = -2.0 * ik_291[k]
                   + f_0 * lk_615[k];

        t_436[k] = -2.0 * ik_292[k]
                   + f_0 * lk_616[k];

        t_437[k] = -2.0 * ik_293[k]
                   + f_0 * lk_617[k];

        t_438[k] = -2.0 * ik_294[k]
                   + f_0 * lk_618[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, t_443, ik_295, ik_296, ik_297, ik_298, \
                         ik_299, lk_619, lk_620, lk_621, lk_622, \
                         lk_623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = -2.0 * ik_295[k]
                   + f_0 * lk_619[k];

        t_440[k] = -2.0 * ik_296[k]
                   + f_0 * lk_620[k];

        t_441[k] = -2.0 * ik_297[k]
                   + f_0 * lk_621[k];

        t_442[k] = -2.0 * ik_298[k]
                   + f_0 * lk_622[k];

        t_443[k] = -2.0 * ik_299[k]
                   + f_0 * lk_623[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, ik_300, ik_301, ik_302, ik_303, \
                         ik_304, lk_624, lk_625, lk_626, lk_627, \
                         lk_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = -2.0 * ik_300[k]
                   + f_0 * lk_624[k];

        t_445[k] = -2.0 * ik_301[k]
                   + f_0 * lk_625[k];

        t_446[k] = -2.0 * ik_302[k]
                   + f_0 * lk_626[k];

        t_447[k] = -2.0 * ik_303[k]
                   + f_0 * lk_627[k];

        t_448[k] = -2.0 * ik_304[k]
                   + f_0 * lk_628[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, ik_305, ik_306, ik_307, ik_308, \
                         ik_309, lk_629, lk_630, lk_631, lk_632, \
                         lk_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = -2.0 * ik_305[k]
                   + f_0 * lk_629[k];

        t_450[k] = -2.0 * ik_306[k]
                   + f_0 * lk_630[k];

        t_451[k] = -2.0 * ik_307[k]
                   + f_0 * lk_631[k];

        t_452[k] = -2.0 * ik_308[k]
                   + f_0 * lk_632[k];

        t_453[k] = -2.0 * ik_309[k]
                   + f_0 * lk_633[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, ik_310, ik_311, ik_312, ik_313, \
                         ik_314, lk_634, lk_635, lk_636, lk_637, \
                         lk_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = -2.0 * ik_310[k]
                   + f_0 * lk_634[k];

        t_455[k] = -2.0 * ik_311[k]
                   + f_0 * lk_635[k];

        t_456[k] = -2.0 * ik_312[k]
                   + f_0 * lk_636[k];

        t_457[k] = -2.0 * ik_313[k]
                   + f_0 * lk_637[k];

        t_458[k] = -2.0 * ik_314[k]
                   + f_0 * lk_638[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, t_463, ik_315, ik_316, ik_317, ik_318, \
                         ik_319, lk_639, lk_640, lk_641, lk_642, \
                         lk_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = -2.0 * ik_315[k]
                   + f_0 * lk_639[k];

        t_460[k] = -2.0 * ik_316[k]
                   + f_0 * lk_640[k];

        t_461[k] = -2.0 * ik_317[k]
                   + f_0 * lk_641[k];

        t_462[k] = -2.0 * ik_318[k]
                   + f_0 * lk_642[k];

        t_463[k] = -2.0 * ik_319[k]
                   + f_0 * lk_643[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, t_468, ik_320, ik_321, ik_322, ik_323, \
                         ik_324, lk_644, lk_645, lk_646, lk_647, \
                         lk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = -2.0 * ik_320[k]
                   + f_0 * lk_644[k];

        t_465[k] = -2.0 * ik_321[k]
                   + f_0 * lk_645[k];

        t_466[k] = -2.0 * ik_322[k]
                   + f_0 * lk_646[k];

        t_467[k] = -2.0 * ik_323[k]
                   + f_0 * lk_647[k];

        t_468[k] = -ik_324[k]
                   + f_0 * lk_648[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, t_473, ik_325, ik_326, ik_327, ik_328, \
                         ik_329, lk_649, lk_650, lk_651, lk_652, \
                         lk_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = -ik_325[k]
                   + f_0 * lk_649[k];

        t_470[k] = -ik_326[k]
                   + f_0 * lk_650[k];

        t_471[k] = -ik_327[k]
                   + f_0 * lk_651[k];

        t_472[k] = -ik_328[k]
                   + f_0 * lk_652[k];

        t_473[k] = -ik_329[k]
                   + f_0 * lk_653[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, t_478, ik_330, ik_331, ik_332, ik_333, \
                         ik_334, lk_654, lk_655, lk_656, lk_657, \
                         lk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = -ik_330[k]
                   + f_0 * lk_654[k];

        t_475[k] = -ik_331[k]
                   + f_0 * lk_655[k];

        t_476[k] = -ik_332[k]
                   + f_0 * lk_656[k];

        t_477[k] = -ik_333[k]
                   + f_0 * lk_657[k];

        t_478[k] = -ik_334[k]
                   + f_0 * lk_658[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, t_483, ik_335, ik_336, ik_337, ik_338, \
                         ik_339, lk_659, lk_660, lk_661, lk_662, \
                         lk_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = -ik_335[k]
                   + f_0 * lk_659[k];

        t_480[k] = -ik_336[k]
                   + f_0 * lk_660[k];

        t_481[k] = -ik_337[k]
                   + f_0 * lk_661[k];

        t_482[k] = -ik_338[k]
                   + f_0 * lk_662[k];

        t_483[k] = -ik_339[k]
                   + f_0 * lk_663[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, t_488, ik_340, ik_341, ik_342, ik_343, \
                         ik_344, lk_664, lk_665, lk_666, lk_667, \
                         lk_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = -ik_340[k]
                   + f_0 * lk_664[k];

        t_485[k] = -ik_341[k]
                   + f_0 * lk_665[k];

        t_486[k] = -ik_342[k]
                   + f_0 * lk_666[k];

        t_487[k] = -ik_343[k]
                   + f_0 * lk_667[k];

        t_488[k] = -ik_344[k]
                   + f_0 * lk_668[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, ik_345, ik_346, ik_347, ik_348, \
                         ik_349, lk_669, lk_670, lk_671, lk_672, \
                         lk_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = -ik_345[k]
                   + f_0 * lk_669[k];

        t_490[k] = -ik_346[k]
                   + f_0 * lk_670[k];

        t_491[k] = -ik_347[k]
                   + f_0 * lk_671[k];

        t_492[k] = -ik_348[k]
                   + f_0 * lk_672[k];

        t_493[k] = -ik_349[k]
                   + f_0 * lk_673[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ik_350 = buffer.data(ik + 350);
    const auto *ik_351 = buffer.data(ik + 351);
    const auto *ik_352 = buffer.data(ik + 352);
    const auto *ik_353 = buffer.data(ik + 353);
    const auto *ik_354 = buffer.data(ik + 354);
    const auto *ik_355 = buffer.data(ik + 355);
    const auto *ik_356 = buffer.data(ik + 356);
    const auto *ik_357 = buffer.data(ik + 357);
    const auto *ik_358 = buffer.data(ik + 358);
    const auto *ik_359 = buffer.data(ik + 359);
    const auto *ik_360 = buffer.data(ik + 360);
    const auto *ik_361 = buffer.data(ik + 361);
    const auto *ik_362 = buffer.data(ik + 362);
    const auto *ik_363 = buffer.data(ik + 363);
    const auto *ik_364 = buffer.data(ik + 364);
    const auto *ik_365 = buffer.data(ik + 365);
    const auto *ik_366 = buffer.data(ik + 366);
    const auto *ik_367 = buffer.data(ik + 367);
    const auto *ik_368 = buffer.data(ik + 368);
    const auto *ik_369 = buffer.data(ik + 369);
    const auto *ik_370 = buffer.data(ik + 370);
    const auto *ik_371 = buffer.data(ik + 371);
    const auto *ik_372 = buffer.data(ik + 372);
    const auto *ik_373 = buffer.data(ik + 373);
    const auto *ik_374 = buffer.data(ik + 374);
    const auto *ik_375 = buffer.data(ik + 375);
    const auto *ik_376 = buffer.data(ik + 376);
    const auto *ik_377 = buffer.data(ik + 377);
    const auto *ik_378 = buffer.data(ik + 378);
    const auto *ik_379 = buffer.data(ik + 379);
    const auto *ik_380 = buffer.data(ik + 380);
    const auto *ik_381 = buffer.data(ik + 381);
    const auto *ik_382 = buffer.data(ik + 382);
    const auto *ik_383 = buffer.data(ik + 383);
    const auto *ik_384 = buffer.data(ik + 384);
    const auto *ik_385 = buffer.data(ik + 385);
    const auto *ik_386 = buffer.data(ik + 386);
    const auto *ik_387 = buffer.data(ik + 387);
    const auto *ik_388 = buffer.data(ik + 388);
    const auto *ik_389 = buffer.data(ik + 389);
    const auto *ik_390 = buffer.data(ik + 390);
    const auto *ik_391 = buffer.data(ik + 391);
    const auto *ik_392 = buffer.data(ik + 392);
    const auto *ik_393 = buffer.data(ik + 393);
    const auto *ik_394 = buffer.data(ik + 394);
    const auto *ik_395 = buffer.data(ik + 395);
    const auto *ik_396 = buffer.data(ik + 396);
    const auto *ik_397 = buffer.data(ik + 397);
    const auto *ik_398 = buffer.data(ik + 398);
    const auto *ik_399 = buffer.data(ik + 399);
    const auto *ik_400 = buffer.data(ik + 400);
    const auto *ik_401 = buffer.data(ik + 401);
    const auto *ik_402 = buffer.data(ik + 402);
    const auto *ik_403 = buffer.data(ik + 403);
    const auto *ik_404 = buffer.data(ik + 404);
    const auto *ik_405 = buffer.data(ik + 405);
    const auto *ik_406 = buffer.data(ik + 406);
    const auto *ik_407 = buffer.data(ik + 407);
    const auto *ik_408 = buffer.data(ik + 408);
    const auto *ik_409 = buffer.data(ik + 409);
    const auto *ik_410 = buffer.data(ik + 410);
    const auto *ik_411 = buffer.data(ik + 411);
    const auto *ik_412 = buffer.data(ik + 412);
    const auto *ik_413 = buffer.data(ik + 413);
    const auto *ik_414 = buffer.data(ik + 414);
    const auto *ik_415 = buffer.data(ik + 415);
    const auto *ik_416 = buffer.data(ik + 416);
    const auto *ik_417 = buffer.data(ik + 417);
    const auto *ik_418 = buffer.data(ik + 418);
    const auto *ik_419 = buffer.data(ik + 419);
    const auto *ik_420 = buffer.data(ik + 420);
    const auto *ik_421 = buffer.data(ik + 421);
    const auto *ik_422 = buffer.data(ik + 422);
    const auto *ik_423 = buffer.data(ik + 423);
    const auto *ik_424 = buffer.data(ik + 424);
    const auto *ik_425 = buffer.data(ik + 425);
    const auto *ik_426 = buffer.data(ik + 426);
    const auto *ik_427 = buffer.data(ik + 427);
    const auto *ik_428 = buffer.data(ik + 428);
    const auto *ik_429 = buffer.data(ik + 429);
    const auto *ik_430 = buffer.data(ik + 430);
    const auto *ik_431 = buffer.data(ik + 431);
    const auto *ik_432 = buffer.data(ik + 432);
    const auto *ik_433 = buffer.data(ik + 433);
    const auto *ik_434 = buffer.data(ik + 434);
    const auto *ik_435 = buffer.data(ik + 435);
    const auto *ik_436 = buffer.data(ik + 436);
    const auto *ik_437 = buffer.data(ik + 437);
    const auto *ik_438 = buffer.data(ik + 438);
    const auto *ik_439 = buffer.data(ik + 439);
    const auto *ik_440 = buffer.data(ik + 440);
    const auto *ik_441 = buffer.data(ik + 441);
    const auto *ik_442 = buffer.data(ik + 442);
    const auto *ik_443 = buffer.data(ik + 443);
    const auto *ik_444 = buffer.data(ik + 444);
    const auto *ik_445 = buffer.data(ik + 445);
    const auto *ik_446 = buffer.data(ik + 446);
    const auto *ik_447 = buffer.data(ik + 447);
    const auto *ik_448 = buffer.data(ik + 448);
    const auto *ik_449 = buffer.data(ik + 449);
    const auto *ik_450 = buffer.data(ik + 450);
    const auto *ik_451 = buffer.data(ik + 451);
    const auto *ik_452 = buffer.data(ik + 452);
    const auto *ik_453 = buffer.data(ik + 453);
    const auto *ik_454 = buffer.data(ik + 454);
    const auto *ik_455 = buffer.data(ik + 455);
    const auto *ik_456 = buffer.data(ik + 456);
    const auto *ik_457 = buffer.data(ik + 457);
    const auto *ik_458 = buffer.data(ik + 458);
    const auto *ik_459 = buffer.data(ik + 459);
    const auto *ik_460 = buffer.data(ik + 460);
    const auto *ik_461 = buffer.data(ik + 461);
    const auto *ik_462 = buffer.data(ik + 462);
    const auto *ik_463 = buffer.data(ik + 463);
    const auto *ik_464 = buffer.data(ik + 464);
    const auto *ik_465 = buffer.data(ik + 465);
    const auto *ik_466 = buffer.data(ik + 466);
    const auto *ik_467 = buffer.data(ik + 467);
    const auto *ik_468 = buffer.data(ik + 468);
    const auto *ik_469 = buffer.data(ik + 469);
    const auto *ik_470 = buffer.data(ik + 470);
    const auto *ik_471 = buffer.data(ik + 471);
    const auto *ik_472 = buffer.data(ik + 472);
    const auto *ik_473 = buffer.data(ik + 473);
    const auto *ik_474 = buffer.data(ik + 474);
    const auto *ik_475 = buffer.data(ik + 475);
    const auto *ik_476 = buffer.data(ik + 476);

    const auto *lk_674 = buffer.data(lk + 674);
    const auto *lk_675 = buffer.data(lk + 675);
    const auto *lk_676 = buffer.data(lk + 676);
    const auto *lk_677 = buffer.data(lk + 677);
    const auto *lk_678 = buffer.data(lk + 678);
    const auto *lk_679 = buffer.data(lk + 679);
    const auto *lk_680 = buffer.data(lk + 680);
    const auto *lk_681 = buffer.data(lk + 681);
    const auto *lk_682 = buffer.data(lk + 682);
    const auto *lk_683 = buffer.data(lk + 683);
    const auto *lk_684 = buffer.data(lk + 684);
    const auto *lk_685 = buffer.data(lk + 685);
    const auto *lk_686 = buffer.data(lk + 686);
    const auto *lk_687 = buffer.data(lk + 687);
    const auto *lk_688 = buffer.data(lk + 688);
    const auto *lk_689 = buffer.data(lk + 689);
    const auto *lk_690 = buffer.data(lk + 690);
    const auto *lk_691 = buffer.data(lk + 691);
    const auto *lk_692 = buffer.data(lk + 692);
    const auto *lk_693 = buffer.data(lk + 693);
    const auto *lk_694 = buffer.data(lk + 694);
    const auto *lk_695 = buffer.data(lk + 695);
    const auto *lk_696 = buffer.data(lk + 696);
    const auto *lk_697 = buffer.data(lk + 697);
    const auto *lk_698 = buffer.data(lk + 698);
    const auto *lk_699 = buffer.data(lk + 699);
    const auto *lk_700 = buffer.data(lk + 700);
    const auto *lk_701 = buffer.data(lk + 701);
    const auto *lk_702 = buffer.data(lk + 702);
    const auto *lk_703 = buffer.data(lk + 703);
    const auto *lk_704 = buffer.data(lk + 704);
    const auto *lk_705 = buffer.data(lk + 705);
    const auto *lk_706 = buffer.data(lk + 706);
    const auto *lk_707 = buffer.data(lk + 707);
    const auto *lk_708 = buffer.data(lk + 708);
    const auto *lk_709 = buffer.data(lk + 709);
    const auto *lk_710 = buffer.data(lk + 710);
    const auto *lk_711 = buffer.data(lk + 711);
    const auto *lk_712 = buffer.data(lk + 712);
    const auto *lk_713 = buffer.data(lk + 713);
    const auto *lk_714 = buffer.data(lk + 714);
    const auto *lk_715 = buffer.data(lk + 715);
    const auto *lk_716 = buffer.data(lk + 716);
    const auto *lk_717 = buffer.data(lk + 717);
    const auto *lk_718 = buffer.data(lk + 718);
    const auto *lk_719 = buffer.data(lk + 719);
    const auto *lk_756 = buffer.data(lk + 756);
    const auto *lk_757 = buffer.data(lk + 757);
    const auto *lk_758 = buffer.data(lk + 758);
    const auto *lk_759 = buffer.data(lk + 759);
    const auto *lk_760 = buffer.data(lk + 760);
    const auto *lk_761 = buffer.data(lk + 761);
    const auto *lk_762 = buffer.data(lk + 762);
    const auto *lk_763 = buffer.data(lk + 763);
    const auto *lk_764 = buffer.data(lk + 764);
    const auto *lk_765 = buffer.data(lk + 765);
    const auto *lk_766 = buffer.data(lk + 766);
    const auto *lk_767 = buffer.data(lk + 767);
    const auto *lk_768 = buffer.data(lk + 768);
    const auto *lk_769 = buffer.data(lk + 769);
    const auto *lk_770 = buffer.data(lk + 770);
    const auto *lk_771 = buffer.data(lk + 771);
    const auto *lk_772 = buffer.data(lk + 772);
    const auto *lk_773 = buffer.data(lk + 773);
    const auto *lk_774 = buffer.data(lk + 774);
    const auto *lk_775 = buffer.data(lk + 775);
    const auto *lk_776 = buffer.data(lk + 776);
    const auto *lk_777 = buffer.data(lk + 777);
    const auto *lk_778 = buffer.data(lk + 778);
    const auto *lk_779 = buffer.data(lk + 779);
    const auto *lk_780 = buffer.data(lk + 780);
    const auto *lk_781 = buffer.data(lk + 781);
    const auto *lk_782 = buffer.data(lk + 782);
    const auto *lk_783 = buffer.data(lk + 783);
    const auto *lk_784 = buffer.data(lk + 784);
    const auto *lk_785 = buffer.data(lk + 785);
    const auto *lk_786 = buffer.data(lk + 786);
    const auto *lk_787 = buffer.data(lk + 787);
    const auto *lk_788 = buffer.data(lk + 788);
    const auto *lk_789 = buffer.data(lk + 789);
    const auto *lk_790 = buffer.data(lk + 790);
    const auto *lk_791 = buffer.data(lk + 791);
    const auto *lk_792 = buffer.data(lk + 792);
    const auto *lk_793 = buffer.data(lk + 793);
    const auto *lk_794 = buffer.data(lk + 794);
    const auto *lk_795 = buffer.data(lk + 795);
    const auto *lk_796 = buffer.data(lk + 796);
    const auto *lk_797 = buffer.data(lk + 797);
    const auto *lk_798 = buffer.data(lk + 798);
    const auto *lk_799 = buffer.data(lk + 799);
    const auto *lk_800 = buffer.data(lk + 800);
    const auto *lk_801 = buffer.data(lk + 801);
    const auto *lk_802 = buffer.data(lk + 802);
    const auto *lk_803 = buffer.data(lk + 803);
    const auto *lk_804 = buffer.data(lk + 804);
    const auto *lk_805 = buffer.data(lk + 805);
    const auto *lk_806 = buffer.data(lk + 806);
    const auto *lk_807 = buffer.data(lk + 807);
    const auto *lk_808 = buffer.data(lk + 808);
    const auto *lk_809 = buffer.data(lk + 809);
    const auto *lk_810 = buffer.data(lk + 810);
    const auto *lk_811 = buffer.data(lk + 811);
    const auto *lk_812 = buffer.data(lk + 812);
    const auto *lk_813 = buffer.data(lk + 813);
    const auto *lk_814 = buffer.data(lk + 814);
    const auto *lk_815 = buffer.data(lk + 815);
    const auto *lk_816 = buffer.data(lk + 816);
    const auto *lk_817 = buffer.data(lk + 817);
    const auto *lk_818 = buffer.data(lk + 818);
    const auto *lk_819 = buffer.data(lk + 819);
    const auto *lk_820 = buffer.data(lk + 820);
    const auto *lk_821 = buffer.data(lk + 821);
    const auto *lk_822 = buffer.data(lk + 822);
    const auto *lk_823 = buffer.data(lk + 823);
    const auto *lk_824 = buffer.data(lk + 824);
    const auto *lk_825 = buffer.data(lk + 825);
    const auto *lk_826 = buffer.data(lk + 826);
    const auto *lk_827 = buffer.data(lk + 827);
    const auto *lk_828 = buffer.data(lk + 828);
    const auto *lk_829 = buffer.data(lk + 829);
    const auto *lk_830 = buffer.data(lk + 830);
    const auto *lk_831 = buffer.data(lk + 831);
    const auto *lk_832 = buffer.data(lk + 832);
    const auto *lk_833 = buffer.data(lk + 833);
    const auto *lk_834 = buffer.data(lk + 834);
    const auto *lk_835 = buffer.data(lk + 835);
    const auto *lk_836 = buffer.data(lk + 836);
    const auto *lk_837 = buffer.data(lk + 837);
    const auto *lk_838 = buffer.data(lk + 838);
    const auto *lk_839 = buffer.data(lk + 839);
    const auto *lk_840 = buffer.data(lk + 840);
    const auto *lk_841 = buffer.data(lk + 841);
    const auto *lk_842 = buffer.data(lk + 842);
    const auto *lk_843 = buffer.data(lk + 843);
    const auto *lk_844 = buffer.data(lk + 844);
    const auto *lk_845 = buffer.data(lk + 845);
    const auto *lk_846 = buffer.data(lk + 846);
    const auto *lk_847 = buffer.data(lk + 847);
    const auto *lk_848 = buffer.data(lk + 848);
    const auto *lk_849 = buffer.data(lk + 849);
    const auto *lk_850 = buffer.data(lk + 850);
    const auto *lk_851 = buffer.data(lk + 851);
    const auto *lk_852 = buffer.data(lk + 852);
    const auto *lk_853 = buffer.data(lk + 853);
    const auto *lk_854 = buffer.data(lk + 854);
    const auto *lk_855 = buffer.data(lk + 855);
    const auto *lk_856 = buffer.data(lk + 856);
    const auto *lk_857 = buffer.data(lk + 857);
    const auto *lk_858 = buffer.data(lk + 858);
    const auto *lk_859 = buffer.data(lk + 859);
    const auto *lk_860 = buffer.data(lk + 860);
    const auto *lk_861 = buffer.data(lk + 861);
    const auto *lk_862 = buffer.data(lk + 862);
    const auto *lk_863 = buffer.data(lk + 863);
    const auto *lk_864 = buffer.data(lk + 864);
    const auto *lk_865 = buffer.data(lk + 865);
    const auto *lk_866 = buffer.data(lk + 866);
    const auto *lk_867 = buffer.data(lk + 867);
    const auto *lk_868 = buffer.data(lk + 868);
    const auto *lk_869 = buffer.data(lk + 869);
    const auto *lk_870 = buffer.data(lk + 870);
    const auto *lk_871 = buffer.data(lk + 871);
    const auto *lk_872 = buffer.data(lk + 872);

#pragma omp simd aligned(t_494, t_495, t_496, t_497, t_498, ik_350, ik_351, ik_352, ik_353, \
                         ik_354, lk_674, lk_675, lk_676, lk_677, \
                         lk_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = -ik_350[k]
                   + f_0 * lk_674[k];

        t_495[k] = -ik_351[k]
                   + f_0 * lk_675[k];

        t_496[k] = -ik_352[k]
                   + f_0 * lk_676[k];

        t_497[k] = -ik_353[k]
                   + f_0 * lk_677[k];

        t_498[k] = -ik_354[k]
                   + f_0 * lk_678[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, t_503, ik_355, ik_356, ik_357, ik_358, \
                         ik_359, lk_679, lk_680, lk_681, lk_682, \
                         lk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = -ik_355[k]
                   + f_0 * lk_679[k];

        t_500[k] = -ik_356[k]
                   + f_0 * lk_680[k];

        t_501[k] = -ik_357[k]
                   + f_0 * lk_681[k];

        t_502[k] = -ik_358[k]
                   + f_0 * lk_682[k];

        t_503[k] = -ik_359[k]
                   + f_0 * lk_683[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, t_509, t_510, t_511, lk_684, \
                         lk_685, lk_686, lk_687, lk_688, lk_689, lk_690, \
                         lk_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_0 * lk_684[k];

        t_505[k] = f_0 * lk_685[k];

        t_506[k] = f_0 * lk_686[k];

        t_507[k] = f_0 * lk_687[k];

        t_508[k] = f_0 * lk_688[k];

        t_509[k] = f_0 * lk_689[k];

        t_510[k] = f_0 * lk_690[k];

        t_511[k] = f_0 * lk_691[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, t_517, t_518, t_519, lk_692, \
                         lk_693, lk_694, lk_695, lk_696, lk_697, lk_698, \
                         lk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_0 * lk_692[k];

        t_513[k] = f_0 * lk_693[k];

        t_514[k] = f_0 * lk_694[k];

        t_515[k] = f_0 * lk_695[k];

        t_516[k] = f_0 * lk_696[k];

        t_517[k] = f_0 * lk_697[k];

        t_518[k] = f_0 * lk_698[k];

        t_519[k] = f_0 * lk_699[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, t_525, t_526, t_527, lk_700, \
                         lk_701, lk_702, lk_703, lk_704, lk_705, lk_706, \
                         lk_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_0 * lk_700[k];

        t_521[k] = f_0 * lk_701[k];

        t_522[k] = f_0 * lk_702[k];

        t_523[k] = f_0 * lk_703[k];

        t_524[k] = f_0 * lk_704[k];

        t_525[k] = f_0 * lk_705[k];

        t_526[k] = f_0 * lk_706[k];

        t_527[k] = f_0 * lk_707[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, t_533, t_534, t_535, lk_708, \
                         lk_709, lk_710, lk_711, lk_712, lk_713, lk_714, \
                         lk_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_0 * lk_708[k];

        t_529[k] = f_0 * lk_709[k];

        t_530[k] = f_0 * lk_710[k];

        t_531[k] = f_0 * lk_711[k];

        t_532[k] = f_0 * lk_712[k];

        t_533[k] = f_0 * lk_713[k];

        t_534[k] = f_0 * lk_714[k];

        t_535[k] = f_0 * lk_715[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, t_541, ik_360, ik_361, lk_716, \
                         lk_717, lk_718, lk_719, lk_756, lk_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_0 * lk_716[k];

        t_537[k] = f_0 * lk_717[k];

        t_538[k] = f_0 * lk_718[k];

        t_539[k] = f_0 * lk_719[k];

        t_540[k] = -5.0 * ik_360[k]
                   + f_0 * lk_756[k];

        t_541[k] = -5.0 * ik_361[k]
                   + f_0 * lk_757[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, ik_362, ik_363, ik_364, ik_365, \
                         ik_366, lk_758, lk_759, lk_760, lk_761, \
                         lk_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = -5.0 * ik_362[k]
                   + f_0 * lk_758[k];

        t_543[k] = -5.0 * ik_363[k]
                   + f_0 * lk_759[k];

        t_544[k] = -5.0 * ik_364[k]
                   + f_0 * lk_760[k];

        t_545[k] = -5.0 * ik_365[k]
                   + f_0 * lk_761[k];

        t_546[k] = -5.0 * ik_366[k]
                   + f_0 * lk_762[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, ik_367, ik_368, ik_369, ik_370, \
                         ik_371, lk_763, lk_764, lk_765, lk_766, \
                         lk_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = -5.0 * ik_367[k]
                   + f_0 * lk_763[k];

        t_548[k] = -5.0 * ik_368[k]
                   + f_0 * lk_764[k];

        t_549[k] = -5.0 * ik_369[k]
                   + f_0 * lk_765[k];

        t_550[k] = -5.0 * ik_370[k]
                   + f_0 * lk_766[k];

        t_551[k] = -5.0 * ik_371[k]
                   + f_0 * lk_767[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, t_556, ik_372, ik_373, ik_374, ik_375, \
                         ik_376, lk_768, lk_769, lk_770, lk_771, \
                         lk_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = -5.0 * ik_372[k]
                   + f_0 * lk_768[k];

        t_553[k] = -5.0 * ik_373[k]
                   + f_0 * lk_769[k];

        t_554[k] = -5.0 * ik_374[k]
                   + f_0 * lk_770[k];

        t_555[k] = -5.0 * ik_375[k]
                   + f_0 * lk_771[k];

        t_556[k] = -5.0 * ik_376[k]
                   + f_0 * lk_772[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, ik_377, ik_378, ik_379, ik_380, \
                         ik_381, lk_773, lk_774, lk_775, lk_776, \
                         lk_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = -5.0 * ik_377[k]
                   + f_0 * lk_773[k];

        t_558[k] = -5.0 * ik_378[k]
                   + f_0 * lk_774[k];

        t_559[k] = -5.0 * ik_379[k]
                   + f_0 * lk_775[k];

        t_560[k] = -5.0 * ik_380[k]
                   + f_0 * lk_776[k];

        t_561[k] = -5.0 * ik_381[k]
                   + f_0 * lk_777[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, ik_382, ik_383, ik_384, ik_385, \
                         ik_386, lk_778, lk_779, lk_780, lk_781, \
                         lk_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = -5.0 * ik_382[k]
                   + f_0 * lk_778[k];

        t_563[k] = -5.0 * ik_383[k]
                   + f_0 * lk_779[k];

        t_564[k] = -5.0 * ik_384[k]
                   + f_0 * lk_780[k];

        t_565[k] = -5.0 * ik_385[k]
                   + f_0 * lk_781[k];

        t_566[k] = -5.0 * ik_386[k]
                   + f_0 * lk_782[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, ik_387, ik_388, ik_389, ik_390, \
                         ik_391, lk_783, lk_784, lk_785, lk_786, \
                         lk_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -5.0 * ik_387[k]
                   + f_0 * lk_783[k];

        t_568[k] = -5.0 * ik_388[k]
                   + f_0 * lk_784[k];

        t_569[k] = -5.0 * ik_389[k]
                   + f_0 * lk_785[k];

        t_570[k] = -5.0 * ik_390[k]
                   + f_0 * lk_786[k];

        t_571[k] = -5.0 * ik_391[k]
                   + f_0 * lk_787[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, ik_392, ik_393, ik_394, ik_395, \
                         ik_396, lk_788, lk_789, lk_790, lk_791, \
                         lk_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = -5.0 * ik_392[k]
                   + f_0 * lk_788[k];

        t_573[k] = -5.0 * ik_393[k]
                   + f_0 * lk_789[k];

        t_574[k] = -5.0 * ik_394[k]
                   + f_0 * lk_790[k];

        t_575[k] = -5.0 * ik_395[k]
                   + f_0 * lk_791[k];

        t_576[k] = -4.0 * ik_396[k]
                   + f_0 * lk_792[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, ik_397, ik_398, ik_399, ik_400, \
                         ik_401, lk_793, lk_794, lk_795, lk_796, \
                         lk_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = -4.0 * ik_397[k]
                   + f_0 * lk_793[k];

        t_578[k] = -4.0 * ik_398[k]
                   + f_0 * lk_794[k];

        t_579[k] = -4.0 * ik_399[k]
                   + f_0 * lk_795[k];

        t_580[k] = -4.0 * ik_400[k]
                   + f_0 * lk_796[k];

        t_581[k] = -4.0 * ik_401[k]
                   + f_0 * lk_797[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, ik_402, ik_403, ik_404, ik_405, \
                         ik_406, lk_798, lk_799, lk_800, lk_801, \
                         lk_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = -4.0 * ik_402[k]
                   + f_0 * lk_798[k];

        t_583[k] = -4.0 * ik_403[k]
                   + f_0 * lk_799[k];

        t_584[k] = -4.0 * ik_404[k]
                   + f_0 * lk_800[k];

        t_585[k] = -4.0 * ik_405[k]
                   + f_0 * lk_801[k];

        t_586[k] = -4.0 * ik_406[k]
                   + f_0 * lk_802[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, t_591, ik_407, ik_408, ik_409, ik_410, \
                         ik_411, lk_803, lk_804, lk_805, lk_806, \
                         lk_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = -4.0 * ik_407[k]
                   + f_0 * lk_803[k];

        t_588[k] = -4.0 * ik_408[k]
                   + f_0 * lk_804[k];

        t_589[k] = -4.0 * ik_409[k]
                   + f_0 * lk_805[k];

        t_590[k] = -4.0 * ik_410[k]
                   + f_0 * lk_806[k];

        t_591[k] = -4.0 * ik_411[k]
                   + f_0 * lk_807[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, t_596, ik_412, ik_413, ik_414, ik_415, \
                         ik_416, lk_808, lk_809, lk_810, lk_811, \
                         lk_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = -4.0 * ik_412[k]
                   + f_0 * lk_808[k];

        t_593[k] = -4.0 * ik_413[k]
                   + f_0 * lk_809[k];

        t_594[k] = -4.0 * ik_414[k]
                   + f_0 * lk_810[k];

        t_595[k] = -4.0 * ik_415[k]
                   + f_0 * lk_811[k];

        t_596[k] = -4.0 * ik_416[k]
                   + f_0 * lk_812[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, t_600, t_601, ik_417, ik_418, ik_419, ik_420, \
                         ik_421, lk_813, lk_814, lk_815, lk_816, \
                         lk_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = -4.0 * ik_417[k]
                   + f_0 * lk_813[k];

        t_598[k] = -4.0 * ik_418[k]
                   + f_0 * lk_814[k];

        t_599[k] = -4.0 * ik_419[k]
                   + f_0 * lk_815[k];

        t_600[k] = -4.0 * ik_420[k]
                   + f_0 * lk_816[k];

        t_601[k] = -4.0 * ik_421[k]
                   + f_0 * lk_817[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, ik_422, ik_423, ik_424, ik_425, \
                         ik_426, lk_818, lk_819, lk_820, lk_821, \
                         lk_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = -4.0 * ik_422[k]
                   + f_0 * lk_818[k];

        t_603[k] = -4.0 * ik_423[k]
                   + f_0 * lk_819[k];

        t_604[k] = -4.0 * ik_424[k]
                   + f_0 * lk_820[k];

        t_605[k] = -4.0 * ik_425[k]
                   + f_0 * lk_821[k];

        t_606[k] = -4.0 * ik_426[k]
                   + f_0 * lk_822[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, t_611, ik_427, ik_428, ik_429, ik_430, \
                         ik_431, lk_823, lk_824, lk_825, lk_826, \
                         lk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = -4.0 * ik_427[k]
                   + f_0 * lk_823[k];

        t_608[k] = -4.0 * ik_428[k]
                   + f_0 * lk_824[k];

        t_609[k] = -4.0 * ik_429[k]
                   + f_0 * lk_825[k];

        t_610[k] = -4.0 * ik_430[k]
                   + f_0 * lk_826[k];

        t_611[k] = -4.0 * ik_431[k]
                   + f_0 * lk_827[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, ik_432, ik_433, ik_434, ik_435, \
                         ik_436, lk_828, lk_829, lk_830, lk_831, \
                         lk_832 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = -3.0 * ik_432[k]
                   + f_0 * lk_828[k];

        t_613[k] = -3.0 * ik_433[k]
                   + f_0 * lk_829[k];

        t_614[k] = -3.0 * ik_434[k]
                   + f_0 * lk_830[k];

        t_615[k] = -3.0 * ik_435[k]
                   + f_0 * lk_831[k];

        t_616[k] = -3.0 * ik_436[k]
                   + f_0 * lk_832[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, ik_437, ik_438, ik_439, ik_440, \
                         ik_441, lk_833, lk_834, lk_835, lk_836, \
                         lk_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = -3.0 * ik_437[k]
                   + f_0 * lk_833[k];

        t_618[k] = -3.0 * ik_438[k]
                   + f_0 * lk_834[k];

        t_619[k] = -3.0 * ik_439[k]
                   + f_0 * lk_835[k];

        t_620[k] = -3.0 * ik_440[k]
                   + f_0 * lk_836[k];

        t_621[k] = -3.0 * ik_441[k]
                   + f_0 * lk_837[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, ik_442, ik_443, ik_444, ik_445, \
                         ik_446, lk_838, lk_839, lk_840, lk_841, \
                         lk_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = -3.0 * ik_442[k]
                   + f_0 * lk_838[k];

        t_623[k] = -3.0 * ik_443[k]
                   + f_0 * lk_839[k];

        t_624[k] = -3.0 * ik_444[k]
                   + f_0 * lk_840[k];

        t_625[k] = -3.0 * ik_445[k]
                   + f_0 * lk_841[k];

        t_626[k] = -3.0 * ik_446[k]
                   + f_0 * lk_842[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, t_631, ik_447, ik_448, ik_449, ik_450, \
                         ik_451, lk_843, lk_844, lk_845, lk_846, \
                         lk_847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = -3.0 * ik_447[k]
                   + f_0 * lk_843[k];

        t_628[k] = -3.0 * ik_448[k]
                   + f_0 * lk_844[k];

        t_629[k] = -3.0 * ik_449[k]
                   + f_0 * lk_845[k];

        t_630[k] = -3.0 * ik_450[k]
                   + f_0 * lk_846[k];

        t_631[k] = -3.0 * ik_451[k]
                   + f_0 * lk_847[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, t_636, ik_452, ik_453, ik_454, ik_455, \
                         ik_456, lk_848, lk_849, lk_850, lk_851, \
                         lk_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = -3.0 * ik_452[k]
                   + f_0 * lk_848[k];

        t_633[k] = -3.0 * ik_453[k]
                   + f_0 * lk_849[k];

        t_634[k] = -3.0 * ik_454[k]
                   + f_0 * lk_850[k];

        t_635[k] = -3.0 * ik_455[k]
                   + f_0 * lk_851[k];

        t_636[k] = -3.0 * ik_456[k]
                   + f_0 * lk_852[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, t_641, ik_457, ik_458, ik_459, ik_460, \
                         ik_461, lk_853, lk_854, lk_855, lk_856, \
                         lk_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = -3.0 * ik_457[k]
                   + f_0 * lk_853[k];

        t_638[k] = -3.0 * ik_458[k]
                   + f_0 * lk_854[k];

        t_639[k] = -3.0 * ik_459[k]
                   + f_0 * lk_855[k];

        t_640[k] = -3.0 * ik_460[k]
                   + f_0 * lk_856[k];

        t_641[k] = -3.0 * ik_461[k]
                   + f_0 * lk_857[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, t_646, ik_462, ik_463, ik_464, ik_465, \
                         ik_466, lk_858, lk_859, lk_860, lk_861, \
                         lk_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = -3.0 * ik_462[k]
                   + f_0 * lk_858[k];

        t_643[k] = -3.0 * ik_463[k]
                   + f_0 * lk_859[k];

        t_644[k] = -3.0 * ik_464[k]
                   + f_0 * lk_860[k];

        t_645[k] = -3.0 * ik_465[k]
                   + f_0 * lk_861[k];

        t_646[k] = -3.0 * ik_466[k]
                   + f_0 * lk_862[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, t_651, ik_467, ik_468, ik_469, ik_470, \
                         ik_471, lk_863, lk_864, lk_865, lk_866, \
                         lk_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = -3.0 * ik_467[k]
                   + f_0 * lk_863[k];

        t_648[k] = -2.0 * ik_468[k]
                   + f_0 * lk_864[k];

        t_649[k] = -2.0 * ik_469[k]
                   + f_0 * lk_865[k];

        t_650[k] = -2.0 * ik_470[k]
                   + f_0 * lk_866[k];

        t_651[k] = -2.0 * ik_471[k]
                   + f_0 * lk_867[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, t_656, ik_472, ik_473, ik_474, ik_475, \
                         ik_476, lk_868, lk_869, lk_870, lk_871, \
                         lk_872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = -2.0 * ik_472[k]
                   + f_0 * lk_868[k];

        t_653[k] = -2.0 * ik_473[k]
                   + f_0 * lk_869[k];

        t_654[k] = -2.0 * ik_474[k]
                   + f_0 * lk_870[k];

        t_655[k] = -2.0 * ik_475[k]
                   + f_0 * lk_871[k];

        t_656[k] = -2.0 * ik_476[k]
                   + f_0 * lk_872[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_815 = buffer.data(target + 815);
    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);

    const auto *ik_477 = buffer.data(ik + 477);
    const auto *ik_478 = buffer.data(ik + 478);
    const auto *ik_479 = buffer.data(ik + 479);
    const auto *ik_480 = buffer.data(ik + 480);
    const auto *ik_481 = buffer.data(ik + 481);
    const auto *ik_482 = buffer.data(ik + 482);
    const auto *ik_483 = buffer.data(ik + 483);
    const auto *ik_484 = buffer.data(ik + 484);
    const auto *ik_485 = buffer.data(ik + 485);
    const auto *ik_486 = buffer.data(ik + 486);
    const auto *ik_487 = buffer.data(ik + 487);
    const auto *ik_488 = buffer.data(ik + 488);
    const auto *ik_489 = buffer.data(ik + 489);
    const auto *ik_490 = buffer.data(ik + 490);
    const auto *ik_491 = buffer.data(ik + 491);
    const auto *ik_492 = buffer.data(ik + 492);
    const auto *ik_493 = buffer.data(ik + 493);
    const auto *ik_494 = buffer.data(ik + 494);
    const auto *ik_495 = buffer.data(ik + 495);
    const auto *ik_496 = buffer.data(ik + 496);
    const auto *ik_497 = buffer.data(ik + 497);
    const auto *ik_498 = buffer.data(ik + 498);
    const auto *ik_499 = buffer.data(ik + 499);
    const auto *ik_500 = buffer.data(ik + 500);
    const auto *ik_501 = buffer.data(ik + 501);
    const auto *ik_502 = buffer.data(ik + 502);
    const auto *ik_503 = buffer.data(ik + 503);
    const auto *ik_504 = buffer.data(ik + 504);
    const auto *ik_505 = buffer.data(ik + 505);
    const auto *ik_506 = buffer.data(ik + 506);
    const auto *ik_507 = buffer.data(ik + 507);
    const auto *ik_508 = buffer.data(ik + 508);
    const auto *ik_509 = buffer.data(ik + 509);
    const auto *ik_510 = buffer.data(ik + 510);
    const auto *ik_511 = buffer.data(ik + 511);
    const auto *ik_512 = buffer.data(ik + 512);
    const auto *ik_513 = buffer.data(ik + 513);
    const auto *ik_514 = buffer.data(ik + 514);
    const auto *ik_515 = buffer.data(ik + 515);
    const auto *ik_516 = buffer.data(ik + 516);
    const auto *ik_517 = buffer.data(ik + 517);
    const auto *ik_518 = buffer.data(ik + 518);
    const auto *ik_519 = buffer.data(ik + 519);
    const auto *ik_520 = buffer.data(ik + 520);
    const auto *ik_521 = buffer.data(ik + 521);
    const auto *ik_522 = buffer.data(ik + 522);
    const auto *ik_523 = buffer.data(ik + 523);
    const auto *ik_524 = buffer.data(ik + 524);
    const auto *ik_525 = buffer.data(ik + 525);
    const auto *ik_526 = buffer.data(ik + 526);
    const auto *ik_527 = buffer.data(ik + 527);
    const auto *ik_528 = buffer.data(ik + 528);
    const auto *ik_529 = buffer.data(ik + 529);
    const auto *ik_530 = buffer.data(ik + 530);
    const auto *ik_531 = buffer.data(ik + 531);
    const auto *ik_532 = buffer.data(ik + 532);
    const auto *ik_533 = buffer.data(ik + 533);
    const auto *ik_534 = buffer.data(ik + 534);
    const auto *ik_535 = buffer.data(ik + 535);
    const auto *ik_536 = buffer.data(ik + 536);
    const auto *ik_537 = buffer.data(ik + 537);
    const auto *ik_538 = buffer.data(ik + 538);
    const auto *ik_539 = buffer.data(ik + 539);
    const auto *ik_540 = buffer.data(ik + 540);
    const auto *ik_541 = buffer.data(ik + 541);
    const auto *ik_542 = buffer.data(ik + 542);
    const auto *ik_543 = buffer.data(ik + 543);
    const auto *ik_544 = buffer.data(ik + 544);
    const auto *ik_545 = buffer.data(ik + 545);
    const auto *ik_546 = buffer.data(ik + 546);
    const auto *ik_547 = buffer.data(ik + 547);
    const auto *ik_548 = buffer.data(ik + 548);
    const auto *ik_549 = buffer.data(ik + 549);
    const auto *ik_550 = buffer.data(ik + 550);
    const auto *ik_551 = buffer.data(ik + 551);
    const auto *ik_552 = buffer.data(ik + 552);
    const auto *ik_553 = buffer.data(ik + 553);
    const auto *ik_554 = buffer.data(ik + 554);
    const auto *ik_555 = buffer.data(ik + 555);
    const auto *ik_556 = buffer.data(ik + 556);
    const auto *ik_557 = buffer.data(ik + 557);
    const auto *ik_558 = buffer.data(ik + 558);
    const auto *ik_559 = buffer.data(ik + 559);
    const auto *ik_560 = buffer.data(ik + 560);
    const auto *ik_561 = buffer.data(ik + 561);
    const auto *ik_562 = buffer.data(ik + 562);
    const auto *ik_563 = buffer.data(ik + 563);
    const auto *ik_564 = buffer.data(ik + 564);
    const auto *ik_565 = buffer.data(ik + 565);
    const auto *ik_566 = buffer.data(ik + 566);
    const auto *ik_567 = buffer.data(ik + 567);
    const auto *ik_568 = buffer.data(ik + 568);
    const auto *ik_569 = buffer.data(ik + 569);
    const auto *ik_570 = buffer.data(ik + 570);
    const auto *ik_571 = buffer.data(ik + 571);
    const auto *ik_572 = buffer.data(ik + 572);
    const auto *ik_573 = buffer.data(ik + 573);
    const auto *ik_574 = buffer.data(ik + 574);
    const auto *ik_575 = buffer.data(ik + 575);
    const auto *ik_576 = buffer.data(ik + 576);
    const auto *ik_577 = buffer.data(ik + 577);
    const auto *ik_578 = buffer.data(ik + 578);
    const auto *ik_579 = buffer.data(ik + 579);
    const auto *ik_580 = buffer.data(ik + 580);
    const auto *ik_581 = buffer.data(ik + 581);
    const auto *ik_582 = buffer.data(ik + 582);
    const auto *ik_583 = buffer.data(ik + 583);
    const auto *ik_584 = buffer.data(ik + 584);
    const auto *ik_585 = buffer.data(ik + 585);
    const auto *ik_586 = buffer.data(ik + 586);
    const auto *ik_587 = buffer.data(ik + 587);
    const auto *ik_588 = buffer.data(ik + 588);
    const auto *ik_589 = buffer.data(ik + 589);
    const auto *ik_590 = buffer.data(ik + 590);
    const auto *ik_591 = buffer.data(ik + 591);
    const auto *ik_592 = buffer.data(ik + 592);
    const auto *ik_593 = buffer.data(ik + 593);
    const auto *ik_594 = buffer.data(ik + 594);
    const auto *ik_595 = buffer.data(ik + 595);
    const auto *ik_596 = buffer.data(ik + 596);
    const auto *ik_597 = buffer.data(ik + 597);
    const auto *ik_598 = buffer.data(ik + 598);
    const auto *ik_599 = buffer.data(ik + 599);
    const auto *ik_600 = buffer.data(ik + 600);
    const auto *ik_601 = buffer.data(ik + 601);
    const auto *ik_602 = buffer.data(ik + 602);
    const auto *ik_603 = buffer.data(ik + 603);

    const auto *lk_873 = buffer.data(lk + 873);
    const auto *lk_874 = buffer.data(lk + 874);
    const auto *lk_875 = buffer.data(lk + 875);
    const auto *lk_876 = buffer.data(lk + 876);
    const auto *lk_877 = buffer.data(lk + 877);
    const auto *lk_878 = buffer.data(lk + 878);
    const auto *lk_879 = buffer.data(lk + 879);
    const auto *lk_880 = buffer.data(lk + 880);
    const auto *lk_881 = buffer.data(lk + 881);
    const auto *lk_882 = buffer.data(lk + 882);
    const auto *lk_883 = buffer.data(lk + 883);
    const auto *lk_884 = buffer.data(lk + 884);
    const auto *lk_885 = buffer.data(lk + 885);
    const auto *lk_886 = buffer.data(lk + 886);
    const auto *lk_887 = buffer.data(lk + 887);
    const auto *lk_888 = buffer.data(lk + 888);
    const auto *lk_889 = buffer.data(lk + 889);
    const auto *lk_890 = buffer.data(lk + 890);
    const auto *lk_891 = buffer.data(lk + 891);
    const auto *lk_892 = buffer.data(lk + 892);
    const auto *lk_893 = buffer.data(lk + 893);
    const auto *lk_894 = buffer.data(lk + 894);
    const auto *lk_895 = buffer.data(lk + 895);
    const auto *lk_896 = buffer.data(lk + 896);
    const auto *lk_897 = buffer.data(lk + 897);
    const auto *lk_898 = buffer.data(lk + 898);
    const auto *lk_899 = buffer.data(lk + 899);
    const auto *lk_900 = buffer.data(lk + 900);
    const auto *lk_901 = buffer.data(lk + 901);
    const auto *lk_902 = buffer.data(lk + 902);
    const auto *lk_903 = buffer.data(lk + 903);
    const auto *lk_904 = buffer.data(lk + 904);
    const auto *lk_905 = buffer.data(lk + 905);
    const auto *lk_906 = buffer.data(lk + 906);
    const auto *lk_907 = buffer.data(lk + 907);
    const auto *lk_908 = buffer.data(lk + 908);
    const auto *lk_909 = buffer.data(lk + 909);
    const auto *lk_910 = buffer.data(lk + 910);
    const auto *lk_911 = buffer.data(lk + 911);
    const auto *lk_912 = buffer.data(lk + 912);
    const auto *lk_913 = buffer.data(lk + 913);
    const auto *lk_914 = buffer.data(lk + 914);
    const auto *lk_915 = buffer.data(lk + 915);
    const auto *lk_916 = buffer.data(lk + 916);
    const auto *lk_917 = buffer.data(lk + 917);
    const auto *lk_918 = buffer.data(lk + 918);
    const auto *lk_919 = buffer.data(lk + 919);
    const auto *lk_920 = buffer.data(lk + 920);
    const auto *lk_921 = buffer.data(lk + 921);
    const auto *lk_922 = buffer.data(lk + 922);
    const auto *lk_923 = buffer.data(lk + 923);
    const auto *lk_924 = buffer.data(lk + 924);
    const auto *lk_925 = buffer.data(lk + 925);
    const auto *lk_926 = buffer.data(lk + 926);
    const auto *lk_927 = buffer.data(lk + 927);
    const auto *lk_928 = buffer.data(lk + 928);
    const auto *lk_929 = buffer.data(lk + 929);
    const auto *lk_930 = buffer.data(lk + 930);
    const auto *lk_931 = buffer.data(lk + 931);
    const auto *lk_932 = buffer.data(lk + 932);
    const auto *lk_933 = buffer.data(lk + 933);
    const auto *lk_934 = buffer.data(lk + 934);
    const auto *lk_935 = buffer.data(lk + 935);
    const auto *lk_936 = buffer.data(lk + 936);
    const auto *lk_937 = buffer.data(lk + 937);
    const auto *lk_938 = buffer.data(lk + 938);
    const auto *lk_939 = buffer.data(lk + 939);
    const auto *lk_940 = buffer.data(lk + 940);
    const auto *lk_941 = buffer.data(lk + 941);
    const auto *lk_942 = buffer.data(lk + 942);
    const auto *lk_943 = buffer.data(lk + 943);
    const auto *lk_944 = buffer.data(lk + 944);
    const auto *lk_945 = buffer.data(lk + 945);
    const auto *lk_946 = buffer.data(lk + 946);
    const auto *lk_947 = buffer.data(lk + 947);
    const auto *lk_948 = buffer.data(lk + 948);
    const auto *lk_949 = buffer.data(lk + 949);
    const auto *lk_950 = buffer.data(lk + 950);
    const auto *lk_951 = buffer.data(lk + 951);
    const auto *lk_952 = buffer.data(lk + 952);
    const auto *lk_953 = buffer.data(lk + 953);
    const auto *lk_954 = buffer.data(lk + 954);
    const auto *lk_955 = buffer.data(lk + 955);
    const auto *lk_956 = buffer.data(lk + 956);
    const auto *lk_957 = buffer.data(lk + 957);
    const auto *lk_958 = buffer.data(lk + 958);
    const auto *lk_959 = buffer.data(lk + 959);
    const auto *lk_960 = buffer.data(lk + 960);
    const auto *lk_961 = buffer.data(lk + 961);
    const auto *lk_962 = buffer.data(lk + 962);
    const auto *lk_963 = buffer.data(lk + 963);
    const auto *lk_964 = buffer.data(lk + 964);
    const auto *lk_965 = buffer.data(lk + 965);
    const auto *lk_966 = buffer.data(lk + 966);
    const auto *lk_967 = buffer.data(lk + 967);
    const auto *lk_968 = buffer.data(lk + 968);
    const auto *lk_969 = buffer.data(lk + 969);
    const auto *lk_970 = buffer.data(lk + 970);
    const auto *lk_971 = buffer.data(lk + 971);
    const auto *lk_1008 = buffer.data(lk + 1008);
    const auto *lk_1009 = buffer.data(lk + 1009);
    const auto *lk_1010 = buffer.data(lk + 1010);
    const auto *lk_1011 = buffer.data(lk + 1011);
    const auto *lk_1012 = buffer.data(lk + 1012);
    const auto *lk_1013 = buffer.data(lk + 1013);
    const auto *lk_1014 = buffer.data(lk + 1014);
    const auto *lk_1015 = buffer.data(lk + 1015);
    const auto *lk_1016 = buffer.data(lk + 1016);
    const auto *lk_1017 = buffer.data(lk + 1017);
    const auto *lk_1018 = buffer.data(lk + 1018);
    const auto *lk_1019 = buffer.data(lk + 1019);
    const auto *lk_1020 = buffer.data(lk + 1020);
    const auto *lk_1021 = buffer.data(lk + 1021);
    const auto *lk_1022 = buffer.data(lk + 1022);
    const auto *lk_1023 = buffer.data(lk + 1023);
    const auto *lk_1024 = buffer.data(lk + 1024);
    const auto *lk_1025 = buffer.data(lk + 1025);
    const auto *lk_1026 = buffer.data(lk + 1026);
    const auto *lk_1027 = buffer.data(lk + 1027);
    const auto *lk_1028 = buffer.data(lk + 1028);
    const auto *lk_1029 = buffer.data(lk + 1029);
    const auto *lk_1030 = buffer.data(lk + 1030);
    const auto *lk_1031 = buffer.data(lk + 1031);
    const auto *lk_1032 = buffer.data(lk + 1032);
    const auto *lk_1033 = buffer.data(lk + 1033);
    const auto *lk_1034 = buffer.data(lk + 1034);
    const auto *lk_1035 = buffer.data(lk + 1035);
    const auto *lk_1036 = buffer.data(lk + 1036);
    const auto *lk_1037 = buffer.data(lk + 1037);
    const auto *lk_1038 = buffer.data(lk + 1038);
    const auto *lk_1039 = buffer.data(lk + 1039);
    const auto *lk_1040 = buffer.data(lk + 1040);
    const auto *lk_1041 = buffer.data(lk + 1041);
    const auto *lk_1042 = buffer.data(lk + 1042);
    const auto *lk_1043 = buffer.data(lk + 1043);
    const auto *lk_1044 = buffer.data(lk + 1044);
    const auto *lk_1045 = buffer.data(lk + 1045);
    const auto *lk_1046 = buffer.data(lk + 1046);
    const auto *lk_1047 = buffer.data(lk + 1047);
    const auto *lk_1048 = buffer.data(lk + 1048);
    const auto *lk_1049 = buffer.data(lk + 1049);
    const auto *lk_1050 = buffer.data(lk + 1050);
    const auto *lk_1051 = buffer.data(lk + 1051);
    const auto *lk_1052 = buffer.data(lk + 1052);
    const auto *lk_1053 = buffer.data(lk + 1053);
    const auto *lk_1054 = buffer.data(lk + 1054);
    const auto *lk_1055 = buffer.data(lk + 1055);
    const auto *lk_1056 = buffer.data(lk + 1056);
    const auto *lk_1057 = buffer.data(lk + 1057);
    const auto *lk_1058 = buffer.data(lk + 1058);
    const auto *lk_1059 = buffer.data(lk + 1059);
    const auto *lk_1060 = buffer.data(lk + 1060);
    const auto *lk_1061 = buffer.data(lk + 1061);
    const auto *lk_1062 = buffer.data(lk + 1062);
    const auto *lk_1063 = buffer.data(lk + 1063);
    const auto *lk_1064 = buffer.data(lk + 1064);
    const auto *lk_1065 = buffer.data(lk + 1065);
    const auto *lk_1066 = buffer.data(lk + 1066);
    const auto *lk_1067 = buffer.data(lk + 1067);
    const auto *lk_1068 = buffer.data(lk + 1068);
    const auto *lk_1069 = buffer.data(lk + 1069);
    const auto *lk_1070 = buffer.data(lk + 1070);
    const auto *lk_1071 = buffer.data(lk + 1071);

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, ik_477, ik_478, ik_479, ik_480, \
                         ik_481, lk_873, lk_874, lk_875, lk_876, \
                         lk_877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = -2.0 * ik_477[k]
                   + f_0 * lk_873[k];

        t_658[k] = -2.0 * ik_478[k]
                   + f_0 * lk_874[k];

        t_659[k] = -2.0 * ik_479[k]
                   + f_0 * lk_875[k];

        t_660[k] = -2.0 * ik_480[k]
                   + f_0 * lk_876[k];

        t_661[k] = -2.0 * ik_481[k]
                   + f_0 * lk_877[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, t_666, ik_482, ik_483, ik_484, ik_485, \
                         ik_486, lk_878, lk_879, lk_880, lk_881, \
                         lk_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = -2.0 * ik_482[k]
                   + f_0 * lk_878[k];

        t_663[k] = -2.0 * ik_483[k]
                   + f_0 * lk_879[k];

        t_664[k] = -2.0 * ik_484[k]
                   + f_0 * lk_880[k];

        t_665[k] = -2.0 * ik_485[k]
                   + f_0 * lk_881[k];

        t_666[k] = -2.0 * ik_486[k]
                   + f_0 * lk_882[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, t_670, t_671, ik_487, ik_488, ik_489, ik_490, \
                         ik_491, lk_883, lk_884, lk_885, lk_886, \
                         lk_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = -2.0 * ik_487[k]
                   + f_0 * lk_883[k];

        t_668[k] = -2.0 * ik_488[k]
                   + f_0 * lk_884[k];

        t_669[k] = -2.0 * ik_489[k]
                   + f_0 * lk_885[k];

        t_670[k] = -2.0 * ik_490[k]
                   + f_0 * lk_886[k];

        t_671[k] = -2.0 * ik_491[k]
                   + f_0 * lk_887[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, t_676, ik_492, ik_493, ik_494, ik_495, \
                         ik_496, lk_888, lk_889, lk_890, lk_891, \
                         lk_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = -2.0 * ik_492[k]
                   + f_0 * lk_888[k];

        t_673[k] = -2.0 * ik_493[k]
                   + f_0 * lk_889[k];

        t_674[k] = -2.0 * ik_494[k]
                   + f_0 * lk_890[k];

        t_675[k] = -2.0 * ik_495[k]
                   + f_0 * lk_891[k];

        t_676[k] = -2.0 * ik_496[k]
                   + f_0 * lk_892[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, t_680, t_681, ik_497, ik_498, ik_499, ik_500, \
                         ik_501, lk_893, lk_894, lk_895, lk_896, \
                         lk_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = -2.0 * ik_497[k]
                   + f_0 * lk_893[k];

        t_678[k] = -2.0 * ik_498[k]
                   + f_0 * lk_894[k];

        t_679[k] = -2.0 * ik_499[k]
                   + f_0 * lk_895[k];

        t_680[k] = -2.0 * ik_500[k]
                   + f_0 * lk_896[k];

        t_681[k] = -2.0 * ik_501[k]
                   + f_0 * lk_897[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, t_686, ik_502, ik_503, ik_504, ik_505, \
                         ik_506, lk_898, lk_899, lk_900, lk_901, \
                         lk_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = -2.0 * ik_502[k]
                   + f_0 * lk_898[k];

        t_683[k] = -2.0 * ik_503[k]
                   + f_0 * lk_899[k];

        t_684[k] = -ik_504[k]
                   + f_0 * lk_900[k];

        t_685[k] = -ik_505[k]
                   + f_0 * lk_901[k];

        t_686[k] = -ik_506[k]
                   + f_0 * lk_902[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, t_691, ik_507, ik_508, ik_509, ik_510, \
                         ik_511, lk_903, lk_904, lk_905, lk_906, \
                         lk_907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = -ik_507[k]
                   + f_0 * lk_903[k];

        t_688[k] = -ik_508[k]
                   + f_0 * lk_904[k];

        t_689[k] = -ik_509[k]
                   + f_0 * lk_905[k];

        t_690[k] = -ik_510[k]
                   + f_0 * lk_906[k];

        t_691[k] = -ik_511[k]
                   + f_0 * lk_907[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, t_696, ik_512, ik_513, ik_514, ik_515, \
                         ik_516, lk_908, lk_909, lk_910, lk_911, \
                         lk_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = -ik_512[k]
                   + f_0 * lk_908[k];

        t_693[k] = -ik_513[k]
                   + f_0 * lk_909[k];

        t_694[k] = -ik_514[k]
                   + f_0 * lk_910[k];

        t_695[k] = -ik_515[k]
                   + f_0 * lk_911[k];

        t_696[k] = -ik_516[k]
                   + f_0 * lk_912[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, t_701, ik_517, ik_518, ik_519, ik_520, \
                         ik_521, lk_913, lk_914, lk_915, lk_916, \
                         lk_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = -ik_517[k]
                   + f_0 * lk_913[k];

        t_698[k] = -ik_518[k]
                   + f_0 * lk_914[k];

        t_699[k] = -ik_519[k]
                   + f_0 * lk_915[k];

        t_700[k] = -ik_520[k]
                   + f_0 * lk_916[k];

        t_701[k] = -ik_521[k]
                   + f_0 * lk_917[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, t_706, ik_522, ik_523, ik_524, ik_525, \
                         ik_526, lk_918, lk_919, lk_920, lk_921, \
                         lk_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = -ik_522[k]
                   + f_0 * lk_918[k];

        t_703[k] = -ik_523[k]
                   + f_0 * lk_919[k];

        t_704[k] = -ik_524[k]
                   + f_0 * lk_920[k];

        t_705[k] = -ik_525[k]
                   + f_0 * lk_921[k];

        t_706[k] = -ik_526[k]
                   + f_0 * lk_922[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, t_710, t_711, ik_527, ik_528, ik_529, ik_530, \
                         ik_531, lk_923, lk_924, lk_925, lk_926, \
                         lk_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = -ik_527[k]
                   + f_0 * lk_923[k];

        t_708[k] = -ik_528[k]
                   + f_0 * lk_924[k];

        t_709[k] = -ik_529[k]
                   + f_0 * lk_925[k];

        t_710[k] = -ik_530[k]
                   + f_0 * lk_926[k];

        t_711[k] = -ik_531[k]
                   + f_0 * lk_927[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, t_716, ik_532, ik_533, ik_534, ik_535, \
                         ik_536, lk_928, lk_929, lk_930, lk_931, \
                         lk_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = -ik_532[k]
                   + f_0 * lk_928[k];

        t_713[k] = -ik_533[k]
                   + f_0 * lk_929[k];

        t_714[k] = -ik_534[k]
                   + f_0 * lk_930[k];

        t_715[k] = -ik_535[k]
                   + f_0 * lk_931[k];

        t_716[k] = -ik_536[k]
                   + f_0 * lk_932[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, t_720, t_721, t_722, ik_537, ik_538, ik_539, \
                         lk_933, lk_934, lk_935, lk_936, lk_937, \
                         lk_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = -ik_537[k]
                   + f_0 * lk_933[k];

        t_718[k] = -ik_538[k]
                   + f_0 * lk_934[k];

        t_719[k] = -ik_539[k]
                   + f_0 * lk_935[k];

        t_720[k] = f_0 * lk_936[k];

        t_721[k] = f_0 * lk_937[k];

        t_722[k] = f_0 * lk_938[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, t_726, t_727, t_728, t_729, t_730, lk_939, \
                         lk_940, lk_941, lk_942, lk_943, lk_944, lk_945, \
                         lk_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_0 * lk_939[k];

        t_724[k] = f_0 * lk_940[k];

        t_725[k] = f_0 * lk_941[k];

        t_726[k] = f_0 * lk_942[k];

        t_727[k] = f_0 * lk_943[k];

        t_728[k] = f_0 * lk_944[k];

        t_729[k] = f_0 * lk_945[k];

        t_730[k] = f_0 * lk_946[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, t_734, t_735, t_736, t_737, t_738, lk_947, \
                         lk_948, lk_949, lk_950, lk_951, lk_952, lk_953, \
                         lk_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_0 * lk_947[k];

        t_732[k] = f_0 * lk_948[k];

        t_733[k] = f_0 * lk_949[k];

        t_734[k] = f_0 * lk_950[k];

        t_735[k] = f_0 * lk_951[k];

        t_736[k] = f_0 * lk_952[k];

        t_737[k] = f_0 * lk_953[k];

        t_738[k] = f_0 * lk_954[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, t_743, t_744, t_745, t_746, lk_955, \
                         lk_956, lk_957, lk_958, lk_959, lk_960, lk_961, \
                         lk_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_0 * lk_955[k];

        t_740[k] = f_0 * lk_956[k];

        t_741[k] = f_0 * lk_957[k];

        t_742[k] = f_0 * lk_958[k];

        t_743[k] = f_0 * lk_959[k];

        t_744[k] = f_0 * lk_960[k];

        t_745[k] = f_0 * lk_961[k];

        t_746[k] = f_0 * lk_962[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, t_750, t_751, t_752, t_753, t_754, lk_963, \
                         lk_964, lk_965, lk_966, lk_967, lk_968, lk_969, \
                         lk_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_0 * lk_963[k];

        t_748[k] = f_0 * lk_964[k];

        t_749[k] = f_0 * lk_965[k];

        t_750[k] = f_0 * lk_966[k];

        t_751[k] = f_0 * lk_967[k];

        t_752[k] = f_0 * lk_968[k];

        t_753[k] = f_0 * lk_969[k];

        t_754[k] = f_0 * lk_970[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, ik_540, ik_541, ik_542, ik_543, \
                         lk_971, lk_1008, lk_1009, lk_1010, lk_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_0 * lk_971[k];

        t_756[k] = -6.0 * ik_540[k]
                   + f_0 * lk_1008[k];

        t_757[k] = -6.0 * ik_541[k]
                   + f_0 * lk_1009[k];

        t_758[k] = -6.0 * ik_542[k]
                   + f_0 * lk_1010[k];

        t_759[k] = -6.0 * ik_543[k]
                   + f_0 * lk_1011[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, t_764, ik_544, ik_545, ik_546, ik_547, \
                         ik_548, lk_1012, lk_1013, lk_1014, lk_1015, \
                         lk_1016 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = -6.0 * ik_544[k]
                   + f_0 * lk_1012[k];

        t_761[k] = -6.0 * ik_545[k]
                   + f_0 * lk_1013[k];

        t_762[k] = -6.0 * ik_546[k]
                   + f_0 * lk_1014[k];

        t_763[k] = -6.0 * ik_547[k]
                   + f_0 * lk_1015[k];

        t_764[k] = -6.0 * ik_548[k]
                   + f_0 * lk_1016[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, t_769, ik_549, ik_550, ik_551, ik_552, \
                         ik_553, lk_1017, lk_1018, lk_1019, lk_1020, \
                         lk_1021 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = -6.0 * ik_549[k]
                   + f_0 * lk_1017[k];

        t_766[k] = -6.0 * ik_550[k]
                   + f_0 * lk_1018[k];

        t_767[k] = -6.0 * ik_551[k]
                   + f_0 * lk_1019[k];

        t_768[k] = -6.0 * ik_552[k]
                   + f_0 * lk_1020[k];

        t_769[k] = -6.0 * ik_553[k]
                   + f_0 * lk_1021[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, ik_554, ik_555, ik_556, ik_557, \
                         ik_558, lk_1022, lk_1023, lk_1024, lk_1025, \
                         lk_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = -6.0 * ik_554[k]
                   + f_0 * lk_1022[k];

        t_771[k] = -6.0 * ik_555[k]
                   + f_0 * lk_1023[k];

        t_772[k] = -6.0 * ik_556[k]
                   + f_0 * lk_1024[k];

        t_773[k] = -6.0 * ik_557[k]
                   + f_0 * lk_1025[k];

        t_774[k] = -6.0 * ik_558[k]
                   + f_0 * lk_1026[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, ik_559, ik_560, ik_561, ik_562, \
                         ik_563, lk_1027, lk_1028, lk_1029, lk_1030, \
                         lk_1031 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = -6.0 * ik_559[k]
                   + f_0 * lk_1027[k];

        t_776[k] = -6.0 * ik_560[k]
                   + f_0 * lk_1028[k];

        t_777[k] = -6.0 * ik_561[k]
                   + f_0 * lk_1029[k];

        t_778[k] = -6.0 * ik_562[k]
                   + f_0 * lk_1030[k];

        t_779[k] = -6.0 * ik_563[k]
                   + f_0 * lk_1031[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, t_784, ik_564, ik_565, ik_566, ik_567, \
                         ik_568, lk_1032, lk_1033, lk_1034, lk_1035, \
                         lk_1036 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = -6.0 * ik_564[k]
                   + f_0 * lk_1032[k];

        t_781[k] = -6.0 * ik_565[k]
                   + f_0 * lk_1033[k];

        t_782[k] = -6.0 * ik_566[k]
                   + f_0 * lk_1034[k];

        t_783[k] = -6.0 * ik_567[k]
                   + f_0 * lk_1035[k];

        t_784[k] = -6.0 * ik_568[k]
                   + f_0 * lk_1036[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, t_789, ik_569, ik_570, ik_571, ik_572, \
                         ik_573, lk_1037, lk_1038, lk_1039, lk_1040, \
                         lk_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = -6.0 * ik_569[k]
                   + f_0 * lk_1037[k];

        t_786[k] = -6.0 * ik_570[k]
                   + f_0 * lk_1038[k];

        t_787[k] = -6.0 * ik_571[k]
                   + f_0 * lk_1039[k];

        t_788[k] = -6.0 * ik_572[k]
                   + f_0 * lk_1040[k];

        t_789[k] = -6.0 * ik_573[k]
                   + f_0 * lk_1041[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, t_794, ik_574, ik_575, ik_576, ik_577, \
                         ik_578, lk_1042, lk_1043, lk_1044, lk_1045, \
                         lk_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = -6.0 * ik_574[k]
                   + f_0 * lk_1042[k];

        t_791[k] = -6.0 * ik_575[k]
                   + f_0 * lk_1043[k];

        t_792[k] = -5.0 * ik_576[k]
                   + f_0 * lk_1044[k];

        t_793[k] = -5.0 * ik_577[k]
                   + f_0 * lk_1045[k];

        t_794[k] = -5.0 * ik_578[k]
                   + f_0 * lk_1046[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, ik_579, ik_580, ik_581, ik_582, \
                         ik_583, lk_1047, lk_1048, lk_1049, lk_1050, \
                         lk_1051 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = -5.0 * ik_579[k]
                   + f_0 * lk_1047[k];

        t_796[k] = -5.0 * ik_580[k]
                   + f_0 * lk_1048[k];

        t_797[k] = -5.0 * ik_581[k]
                   + f_0 * lk_1049[k];

        t_798[k] = -5.0 * ik_582[k]
                   + f_0 * lk_1050[k];

        t_799[k] = -5.0 * ik_583[k]
                   + f_0 * lk_1051[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, ik_584, ik_585, ik_586, ik_587, \
                         ik_588, lk_1052, lk_1053, lk_1054, lk_1055, \
                         lk_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = -5.0 * ik_584[k]
                   + f_0 * lk_1052[k];

        t_801[k] = -5.0 * ik_585[k]
                   + f_0 * lk_1053[k];

        t_802[k] = -5.0 * ik_586[k]
                   + f_0 * lk_1054[k];

        t_803[k] = -5.0 * ik_587[k]
                   + f_0 * lk_1055[k];

        t_804[k] = -5.0 * ik_588[k]
                   + f_0 * lk_1056[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, ik_589, ik_590, ik_591, ik_592, \
                         ik_593, lk_1057, lk_1058, lk_1059, lk_1060, \
                         lk_1061 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = -5.0 * ik_589[k]
                   + f_0 * lk_1057[k];

        t_806[k] = -5.0 * ik_590[k]
                   + f_0 * lk_1058[k];

        t_807[k] = -5.0 * ik_591[k]
                   + f_0 * lk_1059[k];

        t_808[k] = -5.0 * ik_592[k]
                   + f_0 * lk_1060[k];

        t_809[k] = -5.0 * ik_593[k]
                   + f_0 * lk_1061[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, ik_594, ik_595, ik_596, ik_597, \
                         ik_598, lk_1062, lk_1063, lk_1064, lk_1065, \
                         lk_1066 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = -5.0 * ik_594[k]
                   + f_0 * lk_1062[k];

        t_811[k] = -5.0 * ik_595[k]
                   + f_0 * lk_1063[k];

        t_812[k] = -5.0 * ik_596[k]
                   + f_0 * lk_1064[k];

        t_813[k] = -5.0 * ik_597[k]
                   + f_0 * lk_1065[k];

        t_814[k] = -5.0 * ik_598[k]
                   + f_0 * lk_1066[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, ik_599, ik_600, ik_601, ik_602, \
                         ik_603, lk_1067, lk_1068, lk_1069, lk_1070, \
                         lk_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = -5.0 * ik_599[k]
                   + f_0 * lk_1067[k];

        t_816[k] = -5.0 * ik_600[k]
                   + f_0 * lk_1068[k];

        t_817[k] = -5.0 * ik_601[k]
                   + f_0 * lk_1069[k];

        t_818[k] = -5.0 * ik_602[k]
                   + f_0 * lk_1070[k];

        t_819[k] = -5.0 * ik_603[k]
                   + f_0 * lk_1071[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_1_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);

    const auto *ik_604 = buffer.data(ik + 604);
    const auto *ik_605 = buffer.data(ik + 605);
    const auto *ik_606 = buffer.data(ik + 606);
    const auto *ik_607 = buffer.data(ik + 607);
    const auto *ik_608 = buffer.data(ik + 608);
    const auto *ik_609 = buffer.data(ik + 609);
    const auto *ik_610 = buffer.data(ik + 610);
    const auto *ik_611 = buffer.data(ik + 611);
    const auto *ik_612 = buffer.data(ik + 612);
    const auto *ik_613 = buffer.data(ik + 613);
    const auto *ik_614 = buffer.data(ik + 614);
    const auto *ik_615 = buffer.data(ik + 615);
    const auto *ik_616 = buffer.data(ik + 616);
    const auto *ik_617 = buffer.data(ik + 617);
    const auto *ik_618 = buffer.data(ik + 618);
    const auto *ik_619 = buffer.data(ik + 619);
    const auto *ik_620 = buffer.data(ik + 620);
    const auto *ik_621 = buffer.data(ik + 621);
    const auto *ik_622 = buffer.data(ik + 622);
    const auto *ik_623 = buffer.data(ik + 623);
    const auto *ik_624 = buffer.data(ik + 624);
    const auto *ik_625 = buffer.data(ik + 625);
    const auto *ik_626 = buffer.data(ik + 626);
    const auto *ik_627 = buffer.data(ik + 627);
    const auto *ik_628 = buffer.data(ik + 628);
    const auto *ik_629 = buffer.data(ik + 629);
    const auto *ik_630 = buffer.data(ik + 630);
    const auto *ik_631 = buffer.data(ik + 631);
    const auto *ik_632 = buffer.data(ik + 632);
    const auto *ik_633 = buffer.data(ik + 633);
    const auto *ik_634 = buffer.data(ik + 634);
    const auto *ik_635 = buffer.data(ik + 635);
    const auto *ik_636 = buffer.data(ik + 636);
    const auto *ik_637 = buffer.data(ik + 637);
    const auto *ik_638 = buffer.data(ik + 638);
    const auto *ik_639 = buffer.data(ik + 639);
    const auto *ik_640 = buffer.data(ik + 640);
    const auto *ik_641 = buffer.data(ik + 641);
    const auto *ik_642 = buffer.data(ik + 642);
    const auto *ik_643 = buffer.data(ik + 643);
    const auto *ik_644 = buffer.data(ik + 644);
    const auto *ik_645 = buffer.data(ik + 645);
    const auto *ik_646 = buffer.data(ik + 646);
    const auto *ik_647 = buffer.data(ik + 647);
    const auto *ik_648 = buffer.data(ik + 648);
    const auto *ik_649 = buffer.data(ik + 649);
    const auto *ik_650 = buffer.data(ik + 650);
    const auto *ik_651 = buffer.data(ik + 651);
    const auto *ik_652 = buffer.data(ik + 652);
    const auto *ik_653 = buffer.data(ik + 653);
    const auto *ik_654 = buffer.data(ik + 654);
    const auto *ik_655 = buffer.data(ik + 655);
    const auto *ik_656 = buffer.data(ik + 656);
    const auto *ik_657 = buffer.data(ik + 657);
    const auto *ik_658 = buffer.data(ik + 658);
    const auto *ik_659 = buffer.data(ik + 659);
    const auto *ik_660 = buffer.data(ik + 660);
    const auto *ik_661 = buffer.data(ik + 661);
    const auto *ik_662 = buffer.data(ik + 662);
    const auto *ik_663 = buffer.data(ik + 663);
    const auto *ik_664 = buffer.data(ik + 664);
    const auto *ik_665 = buffer.data(ik + 665);
    const auto *ik_666 = buffer.data(ik + 666);
    const auto *ik_667 = buffer.data(ik + 667);
    const auto *ik_668 = buffer.data(ik + 668);
    const auto *ik_669 = buffer.data(ik + 669);
    const auto *ik_670 = buffer.data(ik + 670);
    const auto *ik_671 = buffer.data(ik + 671);
    const auto *ik_672 = buffer.data(ik + 672);
    const auto *ik_673 = buffer.data(ik + 673);
    const auto *ik_674 = buffer.data(ik + 674);
    const auto *ik_675 = buffer.data(ik + 675);
    const auto *ik_676 = buffer.data(ik + 676);
    const auto *ik_677 = buffer.data(ik + 677);
    const auto *ik_678 = buffer.data(ik + 678);
    const auto *ik_679 = buffer.data(ik + 679);
    const auto *ik_680 = buffer.data(ik + 680);
    const auto *ik_681 = buffer.data(ik + 681);
    const auto *ik_682 = buffer.data(ik + 682);
    const auto *ik_683 = buffer.data(ik + 683);
    const auto *ik_684 = buffer.data(ik + 684);
    const auto *ik_685 = buffer.data(ik + 685);
    const auto *ik_686 = buffer.data(ik + 686);
    const auto *ik_687 = buffer.data(ik + 687);
    const auto *ik_688 = buffer.data(ik + 688);
    const auto *ik_689 = buffer.data(ik + 689);
    const auto *ik_690 = buffer.data(ik + 690);
    const auto *ik_691 = buffer.data(ik + 691);
    const auto *ik_692 = buffer.data(ik + 692);
    const auto *ik_693 = buffer.data(ik + 693);
    const auto *ik_694 = buffer.data(ik + 694);
    const auto *ik_695 = buffer.data(ik + 695);
    const auto *ik_696 = buffer.data(ik + 696);
    const auto *ik_697 = buffer.data(ik + 697);
    const auto *ik_698 = buffer.data(ik + 698);
    const auto *ik_699 = buffer.data(ik + 699);
    const auto *ik_700 = buffer.data(ik + 700);
    const auto *ik_701 = buffer.data(ik + 701);
    const auto *ik_702 = buffer.data(ik + 702);
    const auto *ik_703 = buffer.data(ik + 703);
    const auto *ik_704 = buffer.data(ik + 704);
    const auto *ik_705 = buffer.data(ik + 705);
    const auto *ik_706 = buffer.data(ik + 706);
    const auto *ik_707 = buffer.data(ik + 707);
    const auto *ik_708 = buffer.data(ik + 708);
    const auto *ik_709 = buffer.data(ik + 709);
    const auto *ik_710 = buffer.data(ik + 710);
    const auto *ik_711 = buffer.data(ik + 711);
    const auto *ik_712 = buffer.data(ik + 712);
    const auto *ik_713 = buffer.data(ik + 713);
    const auto *ik_714 = buffer.data(ik + 714);
    const auto *ik_715 = buffer.data(ik + 715);
    const auto *ik_716 = buffer.data(ik + 716);
    const auto *ik_717 = buffer.data(ik + 717);
    const auto *ik_718 = buffer.data(ik + 718);
    const auto *ik_719 = buffer.data(ik + 719);
    const auto *ik_720 = buffer.data(ik + 720);
    const auto *ik_721 = buffer.data(ik + 721);
    const auto *ik_722 = buffer.data(ik + 722);
    const auto *ik_723 = buffer.data(ik + 723);
    const auto *ik_724 = buffer.data(ik + 724);
    const auto *ik_725 = buffer.data(ik + 725);
    const auto *ik_726 = buffer.data(ik + 726);
    const auto *ik_727 = buffer.data(ik + 727);
    const auto *ik_728 = buffer.data(ik + 728);
    const auto *ik_729 = buffer.data(ik + 729);
    const auto *ik_730 = buffer.data(ik + 730);
    const auto *ik_731 = buffer.data(ik + 731);
    const auto *ik_732 = buffer.data(ik + 732);
    const auto *ik_733 = buffer.data(ik + 733);
    const auto *ik_734 = buffer.data(ik + 734);
    const auto *ik_735 = buffer.data(ik + 735);
    const auto *ik_736 = buffer.data(ik + 736);
    const auto *ik_737 = buffer.data(ik + 737);
    const auto *ik_738 = buffer.data(ik + 738);
    const auto *ik_739 = buffer.data(ik + 739);
    const auto *ik_740 = buffer.data(ik + 740);
    const auto *ik_741 = buffer.data(ik + 741);
    const auto *ik_742 = buffer.data(ik + 742);
    const auto *ik_743 = buffer.data(ik + 743);
    const auto *ik_744 = buffer.data(ik + 744);
    const auto *ik_745 = buffer.data(ik + 745);
    const auto *ik_746 = buffer.data(ik + 746);
    const auto *ik_747 = buffer.data(ik + 747);
    const auto *ik_748 = buffer.data(ik + 748);
    const auto *ik_749 = buffer.data(ik + 749);
    const auto *ik_750 = buffer.data(ik + 750);
    const auto *ik_751 = buffer.data(ik + 751);
    const auto *ik_752 = buffer.data(ik + 752);
    const auto *ik_753 = buffer.data(ik + 753);

    const auto *lk_1072 = buffer.data(lk + 1072);
    const auto *lk_1073 = buffer.data(lk + 1073);
    const auto *lk_1074 = buffer.data(lk + 1074);
    const auto *lk_1075 = buffer.data(lk + 1075);
    const auto *lk_1076 = buffer.data(lk + 1076);
    const auto *lk_1077 = buffer.data(lk + 1077);
    const auto *lk_1078 = buffer.data(lk + 1078);
    const auto *lk_1079 = buffer.data(lk + 1079);
    const auto *lk_1080 = buffer.data(lk + 1080);
    const auto *lk_1081 = buffer.data(lk + 1081);
    const auto *lk_1082 = buffer.data(lk + 1082);
    const auto *lk_1083 = buffer.data(lk + 1083);
    const auto *lk_1084 = buffer.data(lk + 1084);
    const auto *lk_1085 = buffer.data(lk + 1085);
    const auto *lk_1086 = buffer.data(lk + 1086);
    const auto *lk_1087 = buffer.data(lk + 1087);
    const auto *lk_1088 = buffer.data(lk + 1088);
    const auto *lk_1089 = buffer.data(lk + 1089);
    const auto *lk_1090 = buffer.data(lk + 1090);
    const auto *lk_1091 = buffer.data(lk + 1091);
    const auto *lk_1092 = buffer.data(lk + 1092);
    const auto *lk_1093 = buffer.data(lk + 1093);
    const auto *lk_1094 = buffer.data(lk + 1094);
    const auto *lk_1095 = buffer.data(lk + 1095);
    const auto *lk_1096 = buffer.data(lk + 1096);
    const auto *lk_1097 = buffer.data(lk + 1097);
    const auto *lk_1098 = buffer.data(lk + 1098);
    const auto *lk_1099 = buffer.data(lk + 1099);
    const auto *lk_1100 = buffer.data(lk + 1100);
    const auto *lk_1101 = buffer.data(lk + 1101);
    const auto *lk_1102 = buffer.data(lk + 1102);
    const auto *lk_1103 = buffer.data(lk + 1103);
    const auto *lk_1104 = buffer.data(lk + 1104);
    const auto *lk_1105 = buffer.data(lk + 1105);
    const auto *lk_1106 = buffer.data(lk + 1106);
    const auto *lk_1107 = buffer.data(lk + 1107);
    const auto *lk_1108 = buffer.data(lk + 1108);
    const auto *lk_1109 = buffer.data(lk + 1109);
    const auto *lk_1110 = buffer.data(lk + 1110);
    const auto *lk_1111 = buffer.data(lk + 1111);
    const auto *lk_1112 = buffer.data(lk + 1112);
    const auto *lk_1113 = buffer.data(lk + 1113);
    const auto *lk_1114 = buffer.data(lk + 1114);
    const auto *lk_1115 = buffer.data(lk + 1115);
    const auto *lk_1116 = buffer.data(lk + 1116);
    const auto *lk_1117 = buffer.data(lk + 1117);
    const auto *lk_1118 = buffer.data(lk + 1118);
    const auto *lk_1119 = buffer.data(lk + 1119);
    const auto *lk_1120 = buffer.data(lk + 1120);
    const auto *lk_1121 = buffer.data(lk + 1121);
    const auto *lk_1122 = buffer.data(lk + 1122);
    const auto *lk_1123 = buffer.data(lk + 1123);
    const auto *lk_1124 = buffer.data(lk + 1124);
    const auto *lk_1125 = buffer.data(lk + 1125);
    const auto *lk_1126 = buffer.data(lk + 1126);
    const auto *lk_1127 = buffer.data(lk + 1127);
    const auto *lk_1128 = buffer.data(lk + 1128);
    const auto *lk_1129 = buffer.data(lk + 1129);
    const auto *lk_1130 = buffer.data(lk + 1130);
    const auto *lk_1131 = buffer.data(lk + 1131);
    const auto *lk_1132 = buffer.data(lk + 1132);
    const auto *lk_1133 = buffer.data(lk + 1133);
    const auto *lk_1134 = buffer.data(lk + 1134);
    const auto *lk_1135 = buffer.data(lk + 1135);
    const auto *lk_1136 = buffer.data(lk + 1136);
    const auto *lk_1137 = buffer.data(lk + 1137);
    const auto *lk_1138 = buffer.data(lk + 1138);
    const auto *lk_1139 = buffer.data(lk + 1139);
    const auto *lk_1140 = buffer.data(lk + 1140);
    const auto *lk_1141 = buffer.data(lk + 1141);
    const auto *lk_1142 = buffer.data(lk + 1142);
    const auto *lk_1143 = buffer.data(lk + 1143);
    const auto *lk_1144 = buffer.data(lk + 1144);
    const auto *lk_1145 = buffer.data(lk + 1145);
    const auto *lk_1146 = buffer.data(lk + 1146);
    const auto *lk_1147 = buffer.data(lk + 1147);
    const auto *lk_1148 = buffer.data(lk + 1148);
    const auto *lk_1149 = buffer.data(lk + 1149);
    const auto *lk_1150 = buffer.data(lk + 1150);
    const auto *lk_1151 = buffer.data(lk + 1151);
    const auto *lk_1152 = buffer.data(lk + 1152);
    const auto *lk_1153 = buffer.data(lk + 1153);
    const auto *lk_1154 = buffer.data(lk + 1154);
    const auto *lk_1155 = buffer.data(lk + 1155);
    const auto *lk_1156 = buffer.data(lk + 1156);
    const auto *lk_1157 = buffer.data(lk + 1157);
    const auto *lk_1158 = buffer.data(lk + 1158);
    const auto *lk_1159 = buffer.data(lk + 1159);
    const auto *lk_1160 = buffer.data(lk + 1160);
    const auto *lk_1161 = buffer.data(lk + 1161);
    const auto *lk_1162 = buffer.data(lk + 1162);
    const auto *lk_1163 = buffer.data(lk + 1163);
    const auto *lk_1164 = buffer.data(lk + 1164);
    const auto *lk_1165 = buffer.data(lk + 1165);
    const auto *lk_1166 = buffer.data(lk + 1166);
    const auto *lk_1167 = buffer.data(lk + 1167);
    const auto *lk_1168 = buffer.data(lk + 1168);
    const auto *lk_1169 = buffer.data(lk + 1169);
    const auto *lk_1170 = buffer.data(lk + 1170);
    const auto *lk_1171 = buffer.data(lk + 1171);
    const auto *lk_1172 = buffer.data(lk + 1172);
    const auto *lk_1173 = buffer.data(lk + 1173);
    const auto *lk_1174 = buffer.data(lk + 1174);
    const auto *lk_1175 = buffer.data(lk + 1175);
    const auto *lk_1176 = buffer.data(lk + 1176);
    const auto *lk_1177 = buffer.data(lk + 1177);
    const auto *lk_1178 = buffer.data(lk + 1178);
    const auto *lk_1179 = buffer.data(lk + 1179);
    const auto *lk_1180 = buffer.data(lk + 1180);
    const auto *lk_1181 = buffer.data(lk + 1181);
    const auto *lk_1182 = buffer.data(lk + 1182);
    const auto *lk_1183 = buffer.data(lk + 1183);
    const auto *lk_1184 = buffer.data(lk + 1184);
    const auto *lk_1185 = buffer.data(lk + 1185);
    const auto *lk_1186 = buffer.data(lk + 1186);
    const auto *lk_1187 = buffer.data(lk + 1187);
    const auto *lk_1188 = buffer.data(lk + 1188);
    const auto *lk_1189 = buffer.data(lk + 1189);
    const auto *lk_1190 = buffer.data(lk + 1190);
    const auto *lk_1191 = buffer.data(lk + 1191);
    const auto *lk_1192 = buffer.data(lk + 1192);
    const auto *lk_1193 = buffer.data(lk + 1193);
    const auto *lk_1194 = buffer.data(lk + 1194);
    const auto *lk_1195 = buffer.data(lk + 1195);
    const auto *lk_1196 = buffer.data(lk + 1196);
    const auto *lk_1197 = buffer.data(lk + 1197);
    const auto *lk_1198 = buffer.data(lk + 1198);
    const auto *lk_1199 = buffer.data(lk + 1199);
    const auto *lk_1200 = buffer.data(lk + 1200);
    const auto *lk_1201 = buffer.data(lk + 1201);
    const auto *lk_1202 = buffer.data(lk + 1202);
    const auto *lk_1203 = buffer.data(lk + 1203);
    const auto *lk_1204 = buffer.data(lk + 1204);
    const auto *lk_1205 = buffer.data(lk + 1205);
    const auto *lk_1206 = buffer.data(lk + 1206);
    const auto *lk_1207 = buffer.data(lk + 1207);
    const auto *lk_1208 = buffer.data(lk + 1208);
    const auto *lk_1209 = buffer.data(lk + 1209);
    const auto *lk_1210 = buffer.data(lk + 1210);
    const auto *lk_1211 = buffer.data(lk + 1211);
    const auto *lk_1212 = buffer.data(lk + 1212);
    const auto *lk_1213 = buffer.data(lk + 1213);
    const auto *lk_1214 = buffer.data(lk + 1214);
    const auto *lk_1215 = buffer.data(lk + 1215);
    const auto *lk_1216 = buffer.data(lk + 1216);
    const auto *lk_1217 = buffer.data(lk + 1217);
    const auto *lk_1218 = buffer.data(lk + 1218);
    const auto *lk_1219 = buffer.data(lk + 1219);
    const auto *lk_1220 = buffer.data(lk + 1220);
    const auto *lk_1221 = buffer.data(lk + 1221);

#pragma omp simd aligned(t_820, t_821, t_822, t_823, t_824, ik_604, ik_605, ik_606, ik_607, \
                         ik_608, lk_1072, lk_1073, lk_1074, lk_1075, \
                         lk_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = -5.0 * ik_604[k]
                   + f_0 * lk_1072[k];

        t_821[k] = -5.0 * ik_605[k]
                   + f_0 * lk_1073[k];

        t_822[k] = -5.0 * ik_606[k]
                   + f_0 * lk_1074[k];

        t_823[k] = -5.0 * ik_607[k]
                   + f_0 * lk_1075[k];

        t_824[k] = -5.0 * ik_608[k]
                   + f_0 * lk_1076[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, ik_609, ik_610, ik_611, ik_612, \
                         ik_613, lk_1077, lk_1078, lk_1079, lk_1080, \
                         lk_1081 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = -5.0 * ik_609[k]
                   + f_0 * lk_1077[k];

        t_826[k] = -5.0 * ik_610[k]
                   + f_0 * lk_1078[k];

        t_827[k] = -5.0 * ik_611[k]
                   + f_0 * lk_1079[k];

        t_828[k] = -4.0 * ik_612[k]
                   + f_0 * lk_1080[k];

        t_829[k] = -4.0 * ik_613[k]
                   + f_0 * lk_1081[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, ik_614, ik_615, ik_616, ik_617, \
                         ik_618, lk_1082, lk_1083, lk_1084, lk_1085, \
                         lk_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = -4.0 * ik_614[k]
                   + f_0 * lk_1082[k];

        t_831[k] = -4.0 * ik_615[k]
                   + f_0 * lk_1083[k];

        t_832[k] = -4.0 * ik_616[k]
                   + f_0 * lk_1084[k];

        t_833[k] = -4.0 * ik_617[k]
                   + f_0 * lk_1085[k];

        t_834[k] = -4.0 * ik_618[k]
                   + f_0 * lk_1086[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, t_838, t_839, ik_619, ik_620, ik_621, ik_622, \
                         ik_623, lk_1087, lk_1088, lk_1089, lk_1090, \
                         lk_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = -4.0 * ik_619[k]
                   + f_0 * lk_1087[k];

        t_836[k] = -4.0 * ik_620[k]
                   + f_0 * lk_1088[k];

        t_837[k] = -4.0 * ik_621[k]
                   + f_0 * lk_1089[k];

        t_838[k] = -4.0 * ik_622[k]
                   + f_0 * lk_1090[k];

        t_839[k] = -4.0 * ik_623[k]
                   + f_0 * lk_1091[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, ik_624, ik_625, ik_626, ik_627, \
                         ik_628, lk_1092, lk_1093, lk_1094, lk_1095, \
                         lk_1096 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = -4.0 * ik_624[k]
                   + f_0 * lk_1092[k];

        t_841[k] = -4.0 * ik_625[k]
                   + f_0 * lk_1093[k];

        t_842[k] = -4.0 * ik_626[k]
                   + f_0 * lk_1094[k];

        t_843[k] = -4.0 * ik_627[k]
                   + f_0 * lk_1095[k];

        t_844[k] = -4.0 * ik_628[k]
                   + f_0 * lk_1096[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, t_849, ik_629, ik_630, ik_631, ik_632, \
                         ik_633, lk_1097, lk_1098, lk_1099, lk_1100, \
                         lk_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = -4.0 * ik_629[k]
                   + f_0 * lk_1097[k];

        t_846[k] = -4.0 * ik_630[k]
                   + f_0 * lk_1098[k];

        t_847[k] = -4.0 * ik_631[k]
                   + f_0 * lk_1099[k];

        t_848[k] = -4.0 * ik_632[k]
                   + f_0 * lk_1100[k];

        t_849[k] = -4.0 * ik_633[k]
                   + f_0 * lk_1101[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, ik_634, ik_635, ik_636, ik_637, \
                         ik_638, lk_1102, lk_1103, lk_1104, lk_1105, \
                         lk_1106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = -4.0 * ik_634[k]
                   + f_0 * lk_1102[k];

        t_851[k] = -4.0 * ik_635[k]
                   + f_0 * lk_1103[k];

        t_852[k] = -4.0 * ik_636[k]
                   + f_0 * lk_1104[k];

        t_853[k] = -4.0 * ik_637[k]
                   + f_0 * lk_1105[k];

        t_854[k] = -4.0 * ik_638[k]
                   + f_0 * lk_1106[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, ik_639, ik_640, ik_641, ik_642, \
                         ik_643, lk_1107, lk_1108, lk_1109, lk_1110, \
                         lk_1111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = -4.0 * ik_639[k]
                   + f_0 * lk_1107[k];

        t_856[k] = -4.0 * ik_640[k]
                   + f_0 * lk_1108[k];

        t_857[k] = -4.0 * ik_641[k]
                   + f_0 * lk_1109[k];

        t_858[k] = -4.0 * ik_642[k]
                   + f_0 * lk_1110[k];

        t_859[k] = -4.0 * ik_643[k]
                   + f_0 * lk_1111[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, t_864, ik_644, ik_645, ik_646, ik_647, \
                         ik_648, lk_1112, lk_1113, lk_1114, lk_1115, \
                         lk_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = -4.0 * ik_644[k]
                   + f_0 * lk_1112[k];

        t_861[k] = -4.0 * ik_645[k]
                   + f_0 * lk_1113[k];

        t_862[k] = -4.0 * ik_646[k]
                   + f_0 * lk_1114[k];

        t_863[k] = -4.0 * ik_647[k]
                   + f_0 * lk_1115[k];

        t_864[k] = -3.0 * ik_648[k]
                   + f_0 * lk_1116[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, t_869, ik_649, ik_650, ik_651, ik_652, \
                         ik_653, lk_1117, lk_1118, lk_1119, lk_1120, \
                         lk_1121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = -3.0 * ik_649[k]
                   + f_0 * lk_1117[k];

        t_866[k] = -3.0 * ik_650[k]
                   + f_0 * lk_1118[k];

        t_867[k] = -3.0 * ik_651[k]
                   + f_0 * lk_1119[k];

        t_868[k] = -3.0 * ik_652[k]
                   + f_0 * lk_1120[k];

        t_869[k] = -3.0 * ik_653[k]
                   + f_0 * lk_1121[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, t_874, ik_654, ik_655, ik_656, ik_657, \
                         ik_658, lk_1122, lk_1123, lk_1124, lk_1125, \
                         lk_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = -3.0 * ik_654[k]
                   + f_0 * lk_1122[k];

        t_871[k] = -3.0 * ik_655[k]
                   + f_0 * lk_1123[k];

        t_872[k] = -3.0 * ik_656[k]
                   + f_0 * lk_1124[k];

        t_873[k] = -3.0 * ik_657[k]
                   + f_0 * lk_1125[k];

        t_874[k] = -3.0 * ik_658[k]
                   + f_0 * lk_1126[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, t_879, ik_659, ik_660, ik_661, ik_662, \
                         ik_663, lk_1127, lk_1128, lk_1129, lk_1130, \
                         lk_1131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = -3.0 * ik_659[k]
                   + f_0 * lk_1127[k];

        t_876[k] = -3.0 * ik_660[k]
                   + f_0 * lk_1128[k];

        t_877[k] = -3.0 * ik_661[k]
                   + f_0 * lk_1129[k];

        t_878[k] = -3.0 * ik_662[k]
                   + f_0 * lk_1130[k];

        t_879[k] = -3.0 * ik_663[k]
                   + f_0 * lk_1131[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, t_884, ik_664, ik_665, ik_666, ik_667, \
                         ik_668, lk_1132, lk_1133, lk_1134, lk_1135, \
                         lk_1136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = -3.0 * ik_664[k]
                   + f_0 * lk_1132[k];

        t_881[k] = -3.0 * ik_665[k]
                   + f_0 * lk_1133[k];

        t_882[k] = -3.0 * ik_666[k]
                   + f_0 * lk_1134[k];

        t_883[k] = -3.0 * ik_667[k]
                   + f_0 * lk_1135[k];

        t_884[k] = -3.0 * ik_668[k]
                   + f_0 * lk_1136[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, ik_669, ik_670, ik_671, ik_672, \
                         ik_673, lk_1137, lk_1138, lk_1139, lk_1140, \
                         lk_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = -3.0 * ik_669[k]
                   + f_0 * lk_1137[k];

        t_886[k] = -3.0 * ik_670[k]
                   + f_0 * lk_1138[k];

        t_887[k] = -3.0 * ik_671[k]
                   + f_0 * lk_1139[k];

        t_888[k] = -3.0 * ik_672[k]
                   + f_0 * lk_1140[k];

        t_889[k] = -3.0 * ik_673[k]
                   + f_0 * lk_1141[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, ik_674, ik_675, ik_676, ik_677, \
                         ik_678, lk_1142, lk_1143, lk_1144, lk_1145, \
                         lk_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = -3.0 * ik_674[k]
                   + f_0 * lk_1142[k];

        t_891[k] = -3.0 * ik_675[k]
                   + f_0 * lk_1143[k];

        t_892[k] = -3.0 * ik_676[k]
                   + f_0 * lk_1144[k];

        t_893[k] = -3.0 * ik_677[k]
                   + f_0 * lk_1145[k];

        t_894[k] = -3.0 * ik_678[k]
                   + f_0 * lk_1146[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, t_899, ik_679, ik_680, ik_681, ik_682, \
                         ik_683, lk_1147, lk_1148, lk_1149, lk_1150, \
                         lk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = -3.0 * ik_679[k]
                   + f_0 * lk_1147[k];

        t_896[k] = -3.0 * ik_680[k]
                   + f_0 * lk_1148[k];

        t_897[k] = -3.0 * ik_681[k]
                   + f_0 * lk_1149[k];

        t_898[k] = -3.0 * ik_682[k]
                   + f_0 * lk_1150[k];

        t_899[k] = -3.0 * ik_683[k]
                   + f_0 * lk_1151[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, t_904, ik_684, ik_685, ik_686, ik_687, \
                         ik_688, lk_1152, lk_1153, lk_1154, lk_1155, \
                         lk_1156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = -2.0 * ik_684[k]
                   + f_0 * lk_1152[k];

        t_901[k] = -2.0 * ik_685[k]
                   + f_0 * lk_1153[k];

        t_902[k] = -2.0 * ik_686[k]
                   + f_0 * lk_1154[k];

        t_903[k] = -2.0 * ik_687[k]
                   + f_0 * lk_1155[k];

        t_904[k] = -2.0 * ik_688[k]
                   + f_0 * lk_1156[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, t_909, ik_689, ik_690, ik_691, ik_692, \
                         ik_693, lk_1157, lk_1158, lk_1159, lk_1160, \
                         lk_1161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = -2.0 * ik_689[k]
                   + f_0 * lk_1157[k];

        t_906[k] = -2.0 * ik_690[k]
                   + f_0 * lk_1158[k];

        t_907[k] = -2.0 * ik_691[k]
                   + f_0 * lk_1159[k];

        t_908[k] = -2.0 * ik_692[k]
                   + f_0 * lk_1160[k];

        t_909[k] = -2.0 * ik_693[k]
                   + f_0 * lk_1161[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, t_914, ik_694, ik_695, ik_696, ik_697, \
                         ik_698, lk_1162, lk_1163, lk_1164, lk_1165, \
                         lk_1166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = -2.0 * ik_694[k]
                   + f_0 * lk_1162[k];

        t_911[k] = -2.0 * ik_695[k]
                   + f_0 * lk_1163[k];

        t_912[k] = -2.0 * ik_696[k]
                   + f_0 * lk_1164[k];

        t_913[k] = -2.0 * ik_697[k]
                   + f_0 * lk_1165[k];

        t_914[k] = -2.0 * ik_698[k]
                   + f_0 * lk_1166[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, t_919, ik_699, ik_700, ik_701, ik_702, \
                         ik_703, lk_1167, lk_1168, lk_1169, lk_1170, \
                         lk_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = -2.0 * ik_699[k]
                   + f_0 * lk_1167[k];

        t_916[k] = -2.0 * ik_700[k]
                   + f_0 * lk_1168[k];

        t_917[k] = -2.0 * ik_701[k]
                   + f_0 * lk_1169[k];

        t_918[k] = -2.0 * ik_702[k]
                   + f_0 * lk_1170[k];

        t_919[k] = -2.0 * ik_703[k]
                   + f_0 * lk_1171[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, t_924, ik_704, ik_705, ik_706, ik_707, \
                         ik_708, lk_1172, lk_1173, lk_1174, lk_1175, \
                         lk_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = -2.0 * ik_704[k]
                   + f_0 * lk_1172[k];

        t_921[k] = -2.0 * ik_705[k]
                   + f_0 * lk_1173[k];

        t_922[k] = -2.0 * ik_706[k]
                   + f_0 * lk_1174[k];

        t_923[k] = -2.0 * ik_707[k]
                   + f_0 * lk_1175[k];

        t_924[k] = -2.0 * ik_708[k]
                   + f_0 * lk_1176[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, t_929, ik_709, ik_710, ik_711, ik_712, \
                         ik_713, lk_1177, lk_1178, lk_1179, lk_1180, \
                         lk_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = -2.0 * ik_709[k]
                   + f_0 * lk_1177[k];

        t_926[k] = -2.0 * ik_710[k]
                   + f_0 * lk_1178[k];

        t_927[k] = -2.0 * ik_711[k]
                   + f_0 * lk_1179[k];

        t_928[k] = -2.0 * ik_712[k]
                   + f_0 * lk_1180[k];

        t_929[k] = -2.0 * ik_713[k]
                   + f_0 * lk_1181[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, ik_714, ik_715, ik_716, ik_717, \
                         ik_718, lk_1182, lk_1183, lk_1184, lk_1185, \
                         lk_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = -2.0 * ik_714[k]
                   + f_0 * lk_1182[k];

        t_931[k] = -2.0 * ik_715[k]
                   + f_0 * lk_1183[k];

        t_932[k] = -2.0 * ik_716[k]
                   + f_0 * lk_1184[k];

        t_933[k] = -2.0 * ik_717[k]
                   + f_0 * lk_1185[k];

        t_934[k] = -2.0 * ik_718[k]
                   + f_0 * lk_1186[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, t_939, ik_719, ik_720, ik_721, ik_722, \
                         ik_723, lk_1187, lk_1188, lk_1189, lk_1190, \
                         lk_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = -2.0 * ik_719[k]
                   + f_0 * lk_1187[k];

        t_936[k] = -ik_720[k]
                   + f_0 * lk_1188[k];

        t_937[k] = -ik_721[k]
                   + f_0 * lk_1189[k];

        t_938[k] = -ik_722[k]
                   + f_0 * lk_1190[k];

        t_939[k] = -ik_723[k]
                   + f_0 * lk_1191[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, t_944, ik_724, ik_725, ik_726, ik_727, \
                         ik_728, lk_1192, lk_1193, lk_1194, lk_1195, \
                         lk_1196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = -ik_724[k]
                   + f_0 * lk_1192[k];

        t_941[k] = -ik_725[k]
                   + f_0 * lk_1193[k];

        t_942[k] = -ik_726[k]
                   + f_0 * lk_1194[k];

        t_943[k] = -ik_727[k]
                   + f_0 * lk_1195[k];

        t_944[k] = -ik_728[k]
                   + f_0 * lk_1196[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, t_949, ik_729, ik_730, ik_731, ik_732, \
                         ik_733, lk_1197, lk_1198, lk_1199, lk_1200, \
                         lk_1201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = -ik_729[k]
                   + f_0 * lk_1197[k];

        t_946[k] = -ik_730[k]
                   + f_0 * lk_1198[k];

        t_947[k] = -ik_731[k]
                   + f_0 * lk_1199[k];

        t_948[k] = -ik_732[k]
                   + f_0 * lk_1200[k];

        t_949[k] = -ik_733[k]
                   + f_0 * lk_1201[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, t_954, ik_734, ik_735, ik_736, ik_737, \
                         ik_738, lk_1202, lk_1203, lk_1204, lk_1205, \
                         lk_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = -ik_734[k]
                   + f_0 * lk_1202[k];

        t_951[k] = -ik_735[k]
                   + f_0 * lk_1203[k];

        t_952[k] = -ik_736[k]
                   + f_0 * lk_1204[k];

        t_953[k] = -ik_737[k]
                   + f_0 * lk_1205[k];

        t_954[k] = -ik_738[k]
                   + f_0 * lk_1206[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, t_959, ik_739, ik_740, ik_741, ik_742, \
                         ik_743, lk_1207, lk_1208, lk_1209, lk_1210, \
                         lk_1211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = -ik_739[k]
                   + f_0 * lk_1207[k];

        t_956[k] = -ik_740[k]
                   + f_0 * lk_1208[k];

        t_957[k] = -ik_741[k]
                   + f_0 * lk_1209[k];

        t_958[k] = -ik_742[k]
                   + f_0 * lk_1210[k];

        t_959[k] = -ik_743[k]
                   + f_0 * lk_1211[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, t_963, t_964, ik_744, ik_745, ik_746, ik_747, \
                         ik_748, lk_1212, lk_1213, lk_1214, lk_1215, \
                         lk_1216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = -ik_744[k]
                   + f_0 * lk_1212[k];

        t_961[k] = -ik_745[k]
                   + f_0 * lk_1213[k];

        t_962[k] = -ik_746[k]
                   + f_0 * lk_1214[k];

        t_963[k] = -ik_747[k]
                   + f_0 * lk_1215[k];

        t_964[k] = -ik_748[k]
                   + f_0 * lk_1216[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, t_969, ik_749, ik_750, ik_751, ik_752, \
                         ik_753, lk_1217, lk_1218, lk_1219, lk_1220, \
                         lk_1221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = -ik_749[k]
                   + f_0 * lk_1217[k];

        t_966[k] = -ik_750[k]
                   + f_0 * lk_1218[k];

        t_967[k] = -ik_751[k]
                   + f_0 * lk_1219[k];

        t_968[k] = -ik_752[k]
                   + f_0 * lk_1220[k];

        t_969[k] = -ik_753[k]
                   + f_0 * lk_1221[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_1_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);

    const auto *ik_754 = buffer.data(ik + 754);
    const auto *ik_755 = buffer.data(ik + 755);
    const auto *ik_756 = buffer.data(ik + 756);
    const auto *ik_757 = buffer.data(ik + 757);
    const auto *ik_758 = buffer.data(ik + 758);
    const auto *ik_759 = buffer.data(ik + 759);
    const auto *ik_760 = buffer.data(ik + 760);
    const auto *ik_761 = buffer.data(ik + 761);
    const auto *ik_762 = buffer.data(ik + 762);
    const auto *ik_763 = buffer.data(ik + 763);
    const auto *ik_764 = buffer.data(ik + 764);
    const auto *ik_765 = buffer.data(ik + 765);
    const auto *ik_766 = buffer.data(ik + 766);
    const auto *ik_767 = buffer.data(ik + 767);
    const auto *ik_768 = buffer.data(ik + 768);
    const auto *ik_769 = buffer.data(ik + 769);
    const auto *ik_770 = buffer.data(ik + 770);
    const auto *ik_771 = buffer.data(ik + 771);
    const auto *ik_772 = buffer.data(ik + 772);
    const auto *ik_773 = buffer.data(ik + 773);
    const auto *ik_774 = buffer.data(ik + 774);
    const auto *ik_775 = buffer.data(ik + 775);
    const auto *ik_776 = buffer.data(ik + 776);
    const auto *ik_777 = buffer.data(ik + 777);
    const auto *ik_778 = buffer.data(ik + 778);
    const auto *ik_779 = buffer.data(ik + 779);
    const auto *ik_780 = buffer.data(ik + 780);
    const auto *ik_781 = buffer.data(ik + 781);
    const auto *ik_782 = buffer.data(ik + 782);
    const auto *ik_783 = buffer.data(ik + 783);
    const auto *ik_784 = buffer.data(ik + 784);
    const auto *ik_785 = buffer.data(ik + 785);
    const auto *ik_786 = buffer.data(ik + 786);
    const auto *ik_787 = buffer.data(ik + 787);
    const auto *ik_788 = buffer.data(ik + 788);
    const auto *ik_789 = buffer.data(ik + 789);
    const auto *ik_790 = buffer.data(ik + 790);
    const auto *ik_791 = buffer.data(ik + 791);
    const auto *ik_792 = buffer.data(ik + 792);
    const auto *ik_793 = buffer.data(ik + 793);
    const auto *ik_794 = buffer.data(ik + 794);
    const auto *ik_795 = buffer.data(ik + 795);
    const auto *ik_796 = buffer.data(ik + 796);
    const auto *ik_797 = buffer.data(ik + 797);
    const auto *ik_798 = buffer.data(ik + 798);
    const auto *ik_799 = buffer.data(ik + 799);
    const auto *ik_800 = buffer.data(ik + 800);
    const auto *ik_801 = buffer.data(ik + 801);
    const auto *ik_802 = buffer.data(ik + 802);
    const auto *ik_803 = buffer.data(ik + 803);
    const auto *ik_804 = buffer.data(ik + 804);
    const auto *ik_805 = buffer.data(ik + 805);
    const auto *ik_806 = buffer.data(ik + 806);
    const auto *ik_807 = buffer.data(ik + 807);
    const auto *ik_808 = buffer.data(ik + 808);
    const auto *ik_809 = buffer.data(ik + 809);
    const auto *ik_810 = buffer.data(ik + 810);
    const auto *ik_811 = buffer.data(ik + 811);
    const auto *ik_812 = buffer.data(ik + 812);
    const auto *ik_813 = buffer.data(ik + 813);
    const auto *ik_814 = buffer.data(ik + 814);
    const auto *ik_815 = buffer.data(ik + 815);
    const auto *ik_816 = buffer.data(ik + 816);
    const auto *ik_817 = buffer.data(ik + 817);
    const auto *ik_818 = buffer.data(ik + 818);
    const auto *ik_819 = buffer.data(ik + 819);
    const auto *ik_820 = buffer.data(ik + 820);
    const auto *ik_821 = buffer.data(ik + 821);
    const auto *ik_822 = buffer.data(ik + 822);
    const auto *ik_823 = buffer.data(ik + 823);
    const auto *ik_824 = buffer.data(ik + 824);
    const auto *ik_825 = buffer.data(ik + 825);
    const auto *ik_826 = buffer.data(ik + 826);
    const auto *ik_827 = buffer.data(ik + 827);
    const auto *ik_828 = buffer.data(ik + 828);
    const auto *ik_829 = buffer.data(ik + 829);
    const auto *ik_830 = buffer.data(ik + 830);
    const auto *ik_831 = buffer.data(ik + 831);
    const auto *ik_832 = buffer.data(ik + 832);
    const auto *ik_833 = buffer.data(ik + 833);
    const auto *ik_834 = buffer.data(ik + 834);
    const auto *ik_835 = buffer.data(ik + 835);
    const auto *ik_836 = buffer.data(ik + 836);
    const auto *ik_837 = buffer.data(ik + 837);
    const auto *ik_838 = buffer.data(ik + 838);
    const auto *ik_839 = buffer.data(ik + 839);
    const auto *ik_840 = buffer.data(ik + 840);
    const auto *ik_841 = buffer.data(ik + 841);
    const auto *ik_842 = buffer.data(ik + 842);
    const auto *ik_843 = buffer.data(ik + 843);
    const auto *ik_844 = buffer.data(ik + 844);
    const auto *ik_845 = buffer.data(ik + 845);
    const auto *ik_846 = buffer.data(ik + 846);
    const auto *ik_847 = buffer.data(ik + 847);
    const auto *ik_848 = buffer.data(ik + 848);
    const auto *ik_849 = buffer.data(ik + 849);
    const auto *ik_850 = buffer.data(ik + 850);
    const auto *ik_851 = buffer.data(ik + 851);
    const auto *ik_852 = buffer.data(ik + 852);
    const auto *ik_853 = buffer.data(ik + 853);
    const auto *ik_854 = buffer.data(ik + 854);
    const auto *ik_855 = buffer.data(ik + 855);
    const auto *ik_856 = buffer.data(ik + 856);
    const auto *ik_857 = buffer.data(ik + 857);
    const auto *ik_858 = buffer.data(ik + 858);
    const auto *ik_859 = buffer.data(ik + 859);
    const auto *ik_860 = buffer.data(ik + 860);
    const auto *ik_861 = buffer.data(ik + 861);
    const auto *ik_862 = buffer.data(ik + 862);
    const auto *ik_863 = buffer.data(ik + 863);
    const auto *ik_864 = buffer.data(ik + 864);
    const auto *ik_865 = buffer.data(ik + 865);
    const auto *ik_866 = buffer.data(ik + 866);
    const auto *ik_867 = buffer.data(ik + 867);
    const auto *ik_868 = buffer.data(ik + 868);
    const auto *ik_869 = buffer.data(ik + 869);
    const auto *ik_870 = buffer.data(ik + 870);
    const auto *ik_871 = buffer.data(ik + 871);
    const auto *ik_872 = buffer.data(ik + 872);
    const auto *ik_873 = buffer.data(ik + 873);
    const auto *ik_874 = buffer.data(ik + 874);
    const auto *ik_875 = buffer.data(ik + 875);
    const auto *ik_876 = buffer.data(ik + 876);
    const auto *ik_877 = buffer.data(ik + 877);
    const auto *ik_878 = buffer.data(ik + 878);
    const auto *ik_879 = buffer.data(ik + 879);
    const auto *ik_880 = buffer.data(ik + 880);

    const auto *lk_1222 = buffer.data(lk + 1222);
    const auto *lk_1223 = buffer.data(lk + 1223);
    const auto *lk_1224 = buffer.data(lk + 1224);
    const auto *lk_1225 = buffer.data(lk + 1225);
    const auto *lk_1226 = buffer.data(lk + 1226);
    const auto *lk_1227 = buffer.data(lk + 1227);
    const auto *lk_1228 = buffer.data(lk + 1228);
    const auto *lk_1229 = buffer.data(lk + 1229);
    const auto *lk_1230 = buffer.data(lk + 1230);
    const auto *lk_1231 = buffer.data(lk + 1231);
    const auto *lk_1232 = buffer.data(lk + 1232);
    const auto *lk_1233 = buffer.data(lk + 1233);
    const auto *lk_1234 = buffer.data(lk + 1234);
    const auto *lk_1235 = buffer.data(lk + 1235);
    const auto *lk_1236 = buffer.data(lk + 1236);
    const auto *lk_1237 = buffer.data(lk + 1237);
    const auto *lk_1238 = buffer.data(lk + 1238);
    const auto *lk_1239 = buffer.data(lk + 1239);
    const auto *lk_1240 = buffer.data(lk + 1240);
    const auto *lk_1241 = buffer.data(lk + 1241);
    const auto *lk_1242 = buffer.data(lk + 1242);
    const auto *lk_1243 = buffer.data(lk + 1243);
    const auto *lk_1244 = buffer.data(lk + 1244);
    const auto *lk_1245 = buffer.data(lk + 1245);
    const auto *lk_1246 = buffer.data(lk + 1246);
    const auto *lk_1247 = buffer.data(lk + 1247);
    const auto *lk_1248 = buffer.data(lk + 1248);
    const auto *lk_1249 = buffer.data(lk + 1249);
    const auto *lk_1250 = buffer.data(lk + 1250);
    const auto *lk_1251 = buffer.data(lk + 1251);
    const auto *lk_1252 = buffer.data(lk + 1252);
    const auto *lk_1253 = buffer.data(lk + 1253);
    const auto *lk_1254 = buffer.data(lk + 1254);
    const auto *lk_1255 = buffer.data(lk + 1255);
    const auto *lk_1256 = buffer.data(lk + 1256);
    const auto *lk_1257 = buffer.data(lk + 1257);
    const auto *lk_1258 = buffer.data(lk + 1258);
    const auto *lk_1259 = buffer.data(lk + 1259);
    const auto *lk_1296 = buffer.data(lk + 1296);
    const auto *lk_1297 = buffer.data(lk + 1297);
    const auto *lk_1298 = buffer.data(lk + 1298);
    const auto *lk_1299 = buffer.data(lk + 1299);
    const auto *lk_1300 = buffer.data(lk + 1300);
    const auto *lk_1301 = buffer.data(lk + 1301);
    const auto *lk_1302 = buffer.data(lk + 1302);
    const auto *lk_1303 = buffer.data(lk + 1303);
    const auto *lk_1304 = buffer.data(lk + 1304);
    const auto *lk_1305 = buffer.data(lk + 1305);
    const auto *lk_1306 = buffer.data(lk + 1306);
    const auto *lk_1307 = buffer.data(lk + 1307);
    const auto *lk_1308 = buffer.data(lk + 1308);
    const auto *lk_1309 = buffer.data(lk + 1309);
    const auto *lk_1310 = buffer.data(lk + 1310);
    const auto *lk_1311 = buffer.data(lk + 1311);
    const auto *lk_1312 = buffer.data(lk + 1312);
    const auto *lk_1313 = buffer.data(lk + 1313);
    const auto *lk_1314 = buffer.data(lk + 1314);
    const auto *lk_1315 = buffer.data(lk + 1315);
    const auto *lk_1316 = buffer.data(lk + 1316);
    const auto *lk_1317 = buffer.data(lk + 1317);
    const auto *lk_1318 = buffer.data(lk + 1318);
    const auto *lk_1319 = buffer.data(lk + 1319);
    const auto *lk_1320 = buffer.data(lk + 1320);
    const auto *lk_1321 = buffer.data(lk + 1321);
    const auto *lk_1322 = buffer.data(lk + 1322);
    const auto *lk_1323 = buffer.data(lk + 1323);
    const auto *lk_1324 = buffer.data(lk + 1324);
    const auto *lk_1325 = buffer.data(lk + 1325);
    const auto *lk_1326 = buffer.data(lk + 1326);
    const auto *lk_1327 = buffer.data(lk + 1327);
    const auto *lk_1328 = buffer.data(lk + 1328);
    const auto *lk_1329 = buffer.data(lk + 1329);
    const auto *lk_1330 = buffer.data(lk + 1330);
    const auto *lk_1331 = buffer.data(lk + 1331);
    const auto *lk_1332 = buffer.data(lk + 1332);
    const auto *lk_1333 = buffer.data(lk + 1333);
    const auto *lk_1334 = buffer.data(lk + 1334);
    const auto *lk_1335 = buffer.data(lk + 1335);
    const auto *lk_1336 = buffer.data(lk + 1336);
    const auto *lk_1337 = buffer.data(lk + 1337);
    const auto *lk_1338 = buffer.data(lk + 1338);
    const auto *lk_1339 = buffer.data(lk + 1339);
    const auto *lk_1340 = buffer.data(lk + 1340);
    const auto *lk_1341 = buffer.data(lk + 1341);
    const auto *lk_1342 = buffer.data(lk + 1342);
    const auto *lk_1343 = buffer.data(lk + 1343);
    const auto *lk_1344 = buffer.data(lk + 1344);
    const auto *lk_1345 = buffer.data(lk + 1345);
    const auto *lk_1346 = buffer.data(lk + 1346);
    const auto *lk_1347 = buffer.data(lk + 1347);
    const auto *lk_1348 = buffer.data(lk + 1348);
    const auto *lk_1349 = buffer.data(lk + 1349);
    const auto *lk_1350 = buffer.data(lk + 1350);
    const auto *lk_1351 = buffer.data(lk + 1351);
    const auto *lk_1352 = buffer.data(lk + 1352);
    const auto *lk_1353 = buffer.data(lk + 1353);
    const auto *lk_1354 = buffer.data(lk + 1354);
    const auto *lk_1355 = buffer.data(lk + 1355);
    const auto *lk_1356 = buffer.data(lk + 1356);
    const auto *lk_1357 = buffer.data(lk + 1357);
    const auto *lk_1358 = buffer.data(lk + 1358);
    const auto *lk_1359 = buffer.data(lk + 1359);
    const auto *lk_1360 = buffer.data(lk + 1360);
    const auto *lk_1361 = buffer.data(lk + 1361);
    const auto *lk_1362 = buffer.data(lk + 1362);
    const auto *lk_1363 = buffer.data(lk + 1363);
    const auto *lk_1364 = buffer.data(lk + 1364);
    const auto *lk_1365 = buffer.data(lk + 1365);
    const auto *lk_1366 = buffer.data(lk + 1366);
    const auto *lk_1367 = buffer.data(lk + 1367);
    const auto *lk_1368 = buffer.data(lk + 1368);
    const auto *lk_1369 = buffer.data(lk + 1369);
    const auto *lk_1370 = buffer.data(lk + 1370);
    const auto *lk_1371 = buffer.data(lk + 1371);
    const auto *lk_1372 = buffer.data(lk + 1372);
    const auto *lk_1373 = buffer.data(lk + 1373);
    const auto *lk_1374 = buffer.data(lk + 1374);
    const auto *lk_1375 = buffer.data(lk + 1375);
    const auto *lk_1376 = buffer.data(lk + 1376);
    const auto *lk_1377 = buffer.data(lk + 1377);
    const auto *lk_1378 = buffer.data(lk + 1378);
    const auto *lk_1379 = buffer.data(lk + 1379);
    const auto *lk_1380 = buffer.data(lk + 1380);
    const auto *lk_1381 = buffer.data(lk + 1381);
    const auto *lk_1382 = buffer.data(lk + 1382);
    const auto *lk_1383 = buffer.data(lk + 1383);
    const auto *lk_1384 = buffer.data(lk + 1384);
    const auto *lk_1385 = buffer.data(lk + 1385);
    const auto *lk_1386 = buffer.data(lk + 1386);
    const auto *lk_1387 = buffer.data(lk + 1387);
    const auto *lk_1388 = buffer.data(lk + 1388);
    const auto *lk_1389 = buffer.data(lk + 1389);
    const auto *lk_1390 = buffer.data(lk + 1390);
    const auto *lk_1391 = buffer.data(lk + 1391);
    const auto *lk_1392 = buffer.data(lk + 1392);
    const auto *lk_1393 = buffer.data(lk + 1393);
    const auto *lk_1394 = buffer.data(lk + 1394);
    const auto *lk_1395 = buffer.data(lk + 1395);
    const auto *lk_1396 = buffer.data(lk + 1396);
    const auto *lk_1397 = buffer.data(lk + 1397);
    const auto *lk_1398 = buffer.data(lk + 1398);
    const auto *lk_1399 = buffer.data(lk + 1399);
    const auto *lk_1400 = buffer.data(lk + 1400);
    const auto *lk_1401 = buffer.data(lk + 1401);
    const auto *lk_1402 = buffer.data(lk + 1402);
    const auto *lk_1403 = buffer.data(lk + 1403);
    const auto *lk_1404 = buffer.data(lk + 1404);
    const auto *lk_1405 = buffer.data(lk + 1405);
    const auto *lk_1406 = buffer.data(lk + 1406);
    const auto *lk_1407 = buffer.data(lk + 1407);
    const auto *lk_1408 = buffer.data(lk + 1408);
    const auto *lk_1409 = buffer.data(lk + 1409);
    const auto *lk_1410 = buffer.data(lk + 1410);
    const auto *lk_1411 = buffer.data(lk + 1411);
    const auto *lk_1412 = buffer.data(lk + 1412);
    const auto *lk_1413 = buffer.data(lk + 1413);
    const auto *lk_1414 = buffer.data(lk + 1414);
    const auto *lk_1415 = buffer.data(lk + 1415);
    const auto *lk_1416 = buffer.data(lk + 1416);
    const auto *lk_1417 = buffer.data(lk + 1417);
    const auto *lk_1418 = buffer.data(lk + 1418);
    const auto *lk_1419 = buffer.data(lk + 1419);
    const auto *lk_1420 = buffer.data(lk + 1420);

#pragma omp simd aligned(t_970, t_971, t_972, t_973, t_974, t_975, t_976, ik_754, ik_755, \
                         lk_1222, lk_1223, lk_1224, lk_1225, lk_1226, lk_1227, \
                         lk_1228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = -ik_754[k]
                   + f_0 * lk_1222[k];

        t_971[k] = -ik_755[k]
                   + f_0 * lk_1223[k];

        t_972[k] = f_0 * lk_1224[k];

        t_973[k] = f_0 * lk_1225[k];

        t_974[k] = f_0 * lk_1226[k];

        t_975[k] = f_0 * lk_1227[k];

        t_976[k] = f_0 * lk_1228[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, t_980, t_981, t_982, t_983, t_984, lk_1229, \
                         lk_1230, lk_1231, lk_1232, lk_1233, lk_1234, lk_1235, \
                         lk_1236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_0 * lk_1229[k];

        t_978[k] = f_0 * lk_1230[k];

        t_979[k] = f_0 * lk_1231[k];

        t_980[k] = f_0 * lk_1232[k];

        t_981[k] = f_0 * lk_1233[k];

        t_982[k] = f_0 * lk_1234[k];

        t_983[k] = f_0 * lk_1235[k];

        t_984[k] = f_0 * lk_1236[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, t_989, t_990, t_991, t_992, lk_1237, \
                         lk_1238, lk_1239, lk_1240, lk_1241, lk_1242, lk_1243, \
                         lk_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_0 * lk_1237[k];

        t_986[k] = f_0 * lk_1238[k];

        t_987[k] = f_0 * lk_1239[k];

        t_988[k] = f_0 * lk_1240[k];

        t_989[k] = f_0 * lk_1241[k];

        t_990[k] = f_0 * lk_1242[k];

        t_991[k] = f_0 * lk_1243[k];

        t_992[k] = f_0 * lk_1244[k];
    }

#pragma omp simd aligned(t_993, t_994, t_995, t_996, t_997, t_998, t_999, t_1000, lk_1245, \
                         lk_1246, lk_1247, lk_1248, lk_1249, lk_1250, lk_1251, \
                         lk_1252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_993[k] = f_0 * lk_1245[k];

        t_994[k] = f_0 * lk_1246[k];

        t_995[k] = f_0 * lk_1247[k];

        t_996[k] = f_0 * lk_1248[k];

        t_997[k] = f_0 * lk_1249[k];

        t_998[k] = f_0 * lk_1250[k];

        t_999[k] = f_0 * lk_1251[k];

        t_1000[k] = f_0 * lk_1252[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, t_1004, t_1005, t_1006, t_1007, lk_1253, \
                         lk_1254, lk_1255, lk_1256, lk_1257, lk_1258, \
                         lk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_0 * lk_1253[k];

        t_1002[k] = f_0 * lk_1254[k];

        t_1003[k] = f_0 * lk_1255[k];

        t_1004[k] = f_0 * lk_1256[k];

        t_1005[k] = f_0 * lk_1257[k];

        t_1006[k] = f_0 * lk_1258[k];

        t_1007[k] = f_0 * lk_1259[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, t_1012, ik_756, ik_757, ik_758, \
                         ik_759, ik_760, lk_1296, lk_1297, lk_1298, lk_1299, \
                         lk_1300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = -7.0 * ik_756[k]
                    + f_0 * lk_1296[k];

        t_1009[k] = -7.0 * ik_757[k]
                    + f_0 * lk_1297[k];

        t_1010[k] = -7.0 * ik_758[k]
                    + f_0 * lk_1298[k];

        t_1011[k] = -7.0 * ik_759[k]
                    + f_0 * lk_1299[k];

        t_1012[k] = -7.0 * ik_760[k]
                    + f_0 * lk_1300[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, t_1016, t_1017, ik_761, ik_762, ik_763, \
                         ik_764, ik_765, lk_1301, lk_1302, lk_1303, lk_1304, \
                         lk_1305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = -7.0 * ik_761[k]
                    + f_0 * lk_1301[k];

        t_1014[k] = -7.0 * ik_762[k]
                    + f_0 * lk_1302[k];

        t_1015[k] = -7.0 * ik_763[k]
                    + f_0 * lk_1303[k];

        t_1016[k] = -7.0 * ik_764[k]
                    + f_0 * lk_1304[k];

        t_1017[k] = -7.0 * ik_765[k]
                    + f_0 * lk_1305[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, t_1021, t_1022, ik_766, ik_767, ik_768, \
                         ik_769, ik_770, lk_1306, lk_1307, lk_1308, lk_1309, \
                         lk_1310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = -7.0 * ik_766[k]
                    + f_0 * lk_1306[k];

        t_1019[k] = -7.0 * ik_767[k]
                    + f_0 * lk_1307[k];

        t_1020[k] = -7.0 * ik_768[k]
                    + f_0 * lk_1308[k];

        t_1021[k] = -7.0 * ik_769[k]
                    + f_0 * lk_1309[k];

        t_1022[k] = -7.0 * ik_770[k]
                    + f_0 * lk_1310[k];
    }

#pragma omp simd aligned(t_1023, t_1024, t_1025, t_1026, t_1027, ik_771, ik_772, ik_773, \
                         ik_774, ik_775, lk_1311, lk_1312, lk_1313, lk_1314, \
                         lk_1315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1023[k] = -7.0 * ik_771[k]
                    + f_0 * lk_1311[k];

        t_1024[k] = -7.0 * ik_772[k]
                    + f_0 * lk_1312[k];

        t_1025[k] = -7.0 * ik_773[k]
                    + f_0 * lk_1313[k];

        t_1026[k] = -7.0 * ik_774[k]
                    + f_0 * lk_1314[k];

        t_1027[k] = -7.0 * ik_775[k]
                    + f_0 * lk_1315[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, t_1031, t_1032, ik_776, ik_777, ik_778, \
                         ik_779, ik_780, lk_1316, lk_1317, lk_1318, lk_1319, \
                         lk_1320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = -7.0 * ik_776[k]
                    + f_0 * lk_1316[k];

        t_1029[k] = -7.0 * ik_777[k]
                    + f_0 * lk_1317[k];

        t_1030[k] = -7.0 * ik_778[k]
                    + f_0 * lk_1318[k];

        t_1031[k] = -7.0 * ik_779[k]
                    + f_0 * lk_1319[k];

        t_1032[k] = -7.0 * ik_780[k]
                    + f_0 * lk_1320[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, t_1036, t_1037, ik_781, ik_782, ik_783, \
                         ik_784, ik_785, lk_1321, lk_1322, lk_1323, lk_1324, \
                         lk_1325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = -7.0 * ik_781[k]
                    + f_0 * lk_1321[k];

        t_1034[k] = -7.0 * ik_782[k]
                    + f_0 * lk_1322[k];

        t_1035[k] = -7.0 * ik_783[k]
                    + f_0 * lk_1323[k];

        t_1036[k] = -7.0 * ik_784[k]
                    + f_0 * lk_1324[k];

        t_1037[k] = -7.0 * ik_785[k]
                    + f_0 * lk_1325[k];
    }

#pragma omp simd aligned(t_1038, t_1039, t_1040, t_1041, t_1042, ik_786, ik_787, ik_788, \
                         ik_789, ik_790, lk_1326, lk_1327, lk_1328, lk_1329, \
                         lk_1330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1038[k] = -7.0 * ik_786[k]
                    + f_0 * lk_1326[k];

        t_1039[k] = -7.0 * ik_787[k]
                    + f_0 * lk_1327[k];

        t_1040[k] = -7.0 * ik_788[k]
                    + f_0 * lk_1328[k];

        t_1041[k] = -7.0 * ik_789[k]
                    + f_0 * lk_1329[k];

        t_1042[k] = -7.0 * ik_790[k]
                    + f_0 * lk_1330[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, t_1046, t_1047, ik_791, ik_792, ik_793, \
                         ik_794, ik_795, lk_1331, lk_1332, lk_1333, lk_1334, \
                         lk_1335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = -7.0 * ik_791[k]
                    + f_0 * lk_1331[k];

        t_1044[k] = -6.0 * ik_792[k]
                    + f_0 * lk_1332[k];

        t_1045[k] = -6.0 * ik_793[k]
                    + f_0 * lk_1333[k];

        t_1046[k] = -6.0 * ik_794[k]
                    + f_0 * lk_1334[k];

        t_1047[k] = -6.0 * ik_795[k]
                    + f_0 * lk_1335[k];
    }

#pragma omp simd aligned(t_1048, t_1049, t_1050, t_1051, t_1052, ik_796, ik_797, ik_798, \
                         ik_799, ik_800, lk_1336, lk_1337, lk_1338, lk_1339, \
                         lk_1340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1048[k] = -6.0 * ik_796[k]
                    + f_0 * lk_1336[k];

        t_1049[k] = -6.0 * ik_797[k]
                    + f_0 * lk_1337[k];

        t_1050[k] = -6.0 * ik_798[k]
                    + f_0 * lk_1338[k];

        t_1051[k] = -6.0 * ik_799[k]
                    + f_0 * lk_1339[k];

        t_1052[k] = -6.0 * ik_800[k]
                    + f_0 * lk_1340[k];
    }

#pragma omp simd aligned(t_1053, t_1054, t_1055, t_1056, t_1057, ik_801, ik_802, ik_803, \
                         ik_804, ik_805, lk_1341, lk_1342, lk_1343, lk_1344, \
                         lk_1345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1053[k] = -6.0 * ik_801[k]
                    + f_0 * lk_1341[k];

        t_1054[k] = -6.0 * ik_802[k]
                    + f_0 * lk_1342[k];

        t_1055[k] = -6.0 * ik_803[k]
                    + f_0 * lk_1343[k];

        t_1056[k] = -6.0 * ik_804[k]
                    + f_0 * lk_1344[k];

        t_1057[k] = -6.0 * ik_805[k]
                    + f_0 * lk_1345[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, t_1061, t_1062, ik_806, ik_807, ik_808, \
                         ik_809, ik_810, lk_1346, lk_1347, lk_1348, lk_1349, \
                         lk_1350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = -6.0 * ik_806[k]
                    + f_0 * lk_1346[k];

        t_1059[k] = -6.0 * ik_807[k]
                    + f_0 * lk_1347[k];

        t_1060[k] = -6.0 * ik_808[k]
                    + f_0 * lk_1348[k];

        t_1061[k] = -6.0 * ik_809[k]
                    + f_0 * lk_1349[k];

        t_1062[k] = -6.0 * ik_810[k]
                    + f_0 * lk_1350[k];
    }

#pragma omp simd aligned(t_1063, t_1064, t_1065, t_1066, t_1067, ik_811, ik_812, ik_813, \
                         ik_814, ik_815, lk_1351, lk_1352, lk_1353, lk_1354, \
                         lk_1355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1063[k] = -6.0 * ik_811[k]
                    + f_0 * lk_1351[k];

        t_1064[k] = -6.0 * ik_812[k]
                    + f_0 * lk_1352[k];

        t_1065[k] = -6.0 * ik_813[k]
                    + f_0 * lk_1353[k];

        t_1066[k] = -6.0 * ik_814[k]
                    + f_0 * lk_1354[k];

        t_1067[k] = -6.0 * ik_815[k]
                    + f_0 * lk_1355[k];
    }

#pragma omp simd aligned(t_1068, t_1069, t_1070, t_1071, t_1072, ik_816, ik_817, ik_818, \
                         ik_819, ik_820, lk_1356, lk_1357, lk_1358, lk_1359, \
                         lk_1360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1068[k] = -6.0 * ik_816[k]
                    + f_0 * lk_1356[k];

        t_1069[k] = -6.0 * ik_817[k]
                    + f_0 * lk_1357[k];

        t_1070[k] = -6.0 * ik_818[k]
                    + f_0 * lk_1358[k];

        t_1071[k] = -6.0 * ik_819[k]
                    + f_0 * lk_1359[k];

        t_1072[k] = -6.0 * ik_820[k]
                    + f_0 * lk_1360[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, t_1076, t_1077, ik_821, ik_822, ik_823, \
                         ik_824, ik_825, lk_1361, lk_1362, lk_1363, lk_1364, \
                         lk_1365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = -6.0 * ik_821[k]
                    + f_0 * lk_1361[k];

        t_1074[k] = -6.0 * ik_822[k]
                    + f_0 * lk_1362[k];

        t_1075[k] = -6.0 * ik_823[k]
                    + f_0 * lk_1363[k];

        t_1076[k] = -6.0 * ik_824[k]
                    + f_0 * lk_1364[k];

        t_1077[k] = -6.0 * ik_825[k]
                    + f_0 * lk_1365[k];
    }

#pragma omp simd aligned(t_1078, t_1079, t_1080, t_1081, t_1082, ik_826, ik_827, ik_828, \
                         ik_829, ik_830, lk_1366, lk_1367, lk_1368, lk_1369, \
                         lk_1370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1078[k] = -6.0 * ik_826[k]
                    + f_0 * lk_1366[k];

        t_1079[k] = -6.0 * ik_827[k]
                    + f_0 * lk_1367[k];

        t_1080[k] = -5.0 * ik_828[k]
                    + f_0 * lk_1368[k];

        t_1081[k] = -5.0 * ik_829[k]
                    + f_0 * lk_1369[k];

        t_1082[k] = -5.0 * ik_830[k]
                    + f_0 * lk_1370[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, t_1086, t_1087, ik_831, ik_832, ik_833, \
                         ik_834, ik_835, lk_1371, lk_1372, lk_1373, lk_1374, \
                         lk_1375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = -5.0 * ik_831[k]
                    + f_0 * lk_1371[k];

        t_1084[k] = -5.0 * ik_832[k]
                    + f_0 * lk_1372[k];

        t_1085[k] = -5.0 * ik_833[k]
                    + f_0 * lk_1373[k];

        t_1086[k] = -5.0 * ik_834[k]
                    + f_0 * lk_1374[k];

        t_1087[k] = -5.0 * ik_835[k]
                    + f_0 * lk_1375[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, t_1091, t_1092, ik_836, ik_837, ik_838, \
                         ik_839, ik_840, lk_1376, lk_1377, lk_1378, lk_1379, \
                         lk_1380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = -5.0 * ik_836[k]
                    + f_0 * lk_1376[k];

        t_1089[k] = -5.0 * ik_837[k]
                    + f_0 * lk_1377[k];

        t_1090[k] = -5.0 * ik_838[k]
                    + f_0 * lk_1378[k];

        t_1091[k] = -5.0 * ik_839[k]
                    + f_0 * lk_1379[k];

        t_1092[k] = -5.0 * ik_840[k]
                    + f_0 * lk_1380[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, t_1096, t_1097, ik_841, ik_842, ik_843, \
                         ik_844, ik_845, lk_1381, lk_1382, lk_1383, lk_1384, \
                         lk_1385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = -5.0 * ik_841[k]
                    + f_0 * lk_1381[k];

        t_1094[k] = -5.0 * ik_842[k]
                    + f_0 * lk_1382[k];

        t_1095[k] = -5.0 * ik_843[k]
                    + f_0 * lk_1383[k];

        t_1096[k] = -5.0 * ik_844[k]
                    + f_0 * lk_1384[k];

        t_1097[k] = -5.0 * ik_845[k]
                    + f_0 * lk_1385[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, t_1101, t_1102, ik_846, ik_847, ik_848, \
                         ik_849, ik_850, lk_1386, lk_1387, lk_1388, lk_1389, \
                         lk_1390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = -5.0 * ik_846[k]
                    + f_0 * lk_1386[k];

        t_1099[k] = -5.0 * ik_847[k]
                    + f_0 * lk_1387[k];

        t_1100[k] = -5.0 * ik_848[k]
                    + f_0 * lk_1388[k];

        t_1101[k] = -5.0 * ik_849[k]
                    + f_0 * lk_1389[k];

        t_1102[k] = -5.0 * ik_850[k]
                    + f_0 * lk_1390[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, t_1106, t_1107, ik_851, ik_852, ik_853, \
                         ik_854, ik_855, lk_1391, lk_1392, lk_1393, lk_1394, \
                         lk_1395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = -5.0 * ik_851[k]
                    + f_0 * lk_1391[k];

        t_1104[k] = -5.0 * ik_852[k]
                    + f_0 * lk_1392[k];

        t_1105[k] = -5.0 * ik_853[k]
                    + f_0 * lk_1393[k];

        t_1106[k] = -5.0 * ik_854[k]
                    + f_0 * lk_1394[k];

        t_1107[k] = -5.0 * ik_855[k]
                    + f_0 * lk_1395[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, t_1111, t_1112, ik_856, ik_857, ik_858, \
                         ik_859, ik_860, lk_1396, lk_1397, lk_1398, lk_1399, \
                         lk_1400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = -5.0 * ik_856[k]
                    + f_0 * lk_1396[k];

        t_1109[k] = -5.0 * ik_857[k]
                    + f_0 * lk_1397[k];

        t_1110[k] = -5.0 * ik_858[k]
                    + f_0 * lk_1398[k];

        t_1111[k] = -5.0 * ik_859[k]
                    + f_0 * lk_1399[k];

        t_1112[k] = -5.0 * ik_860[k]
                    + f_0 * lk_1400[k];
    }

#pragma omp simd aligned(t_1113, t_1114, t_1115, t_1116, t_1117, ik_861, ik_862, ik_863, \
                         ik_864, ik_865, lk_1401, lk_1402, lk_1403, lk_1404, \
                         lk_1405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1113[k] = -5.0 * ik_861[k]
                    + f_0 * lk_1401[k];

        t_1114[k] = -5.0 * ik_862[k]
                    + f_0 * lk_1402[k];

        t_1115[k] = -5.0 * ik_863[k]
                    + f_0 * lk_1403[k];

        t_1116[k] = -4.0 * ik_864[k]
                    + f_0 * lk_1404[k];

        t_1117[k] = -4.0 * ik_865[k]
                    + f_0 * lk_1405[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, t_1121, t_1122, ik_866, ik_867, ik_868, \
                         ik_869, ik_870, lk_1406, lk_1407, lk_1408, lk_1409, \
                         lk_1410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = -4.0 * ik_866[k]
                    + f_0 * lk_1406[k];

        t_1119[k] = -4.0 * ik_867[k]
                    + f_0 * lk_1407[k];

        t_1120[k] = -4.0 * ik_868[k]
                    + f_0 * lk_1408[k];

        t_1121[k] = -4.0 * ik_869[k]
                    + f_0 * lk_1409[k];

        t_1122[k] = -4.0 * ik_870[k]
                    + f_0 * lk_1410[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, t_1126, t_1127, ik_871, ik_872, ik_873, \
                         ik_874, ik_875, lk_1411, lk_1412, lk_1413, lk_1414, \
                         lk_1415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = -4.0 * ik_871[k]
                    + f_0 * lk_1411[k];

        t_1124[k] = -4.0 * ik_872[k]
                    + f_0 * lk_1412[k];

        t_1125[k] = -4.0 * ik_873[k]
                    + f_0 * lk_1413[k];

        t_1126[k] = -4.0 * ik_874[k]
                    + f_0 * lk_1414[k];

        t_1127[k] = -4.0 * ik_875[k]
                    + f_0 * lk_1415[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, t_1131, t_1132, ik_876, ik_877, ik_878, \
                         ik_879, ik_880, lk_1416, lk_1417, lk_1418, lk_1419, \
                         lk_1420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = -4.0 * ik_876[k]
                    + f_0 * lk_1416[k];

        t_1129[k] = -4.0 * ik_877[k]
                    + f_0 * lk_1417[k];

        t_1130[k] = -4.0 * ik_878[k]
                    + f_0 * lk_1418[k];

        t_1131[k] = -4.0 * ik_879[k]
                    + f_0 * lk_1419[k];

        t_1132[k] = -4.0 * ik_880[k]
                    + f_0 * lk_1420[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_1_piece7(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_1260 = buffer.data(target + 1260);
    auto *t_1261 = buffer.data(target + 1261);
    auto *t_1262 = buffer.data(target + 1262);
    auto *t_1263 = buffer.data(target + 1263);
    auto *t_1264 = buffer.data(target + 1264);
    auto *t_1265 = buffer.data(target + 1265);
    auto *t_1266 = buffer.data(target + 1266);
    auto *t_1267 = buffer.data(target + 1267);
    auto *t_1268 = buffer.data(target + 1268);
    auto *t_1269 = buffer.data(target + 1269);
    auto *t_1270 = buffer.data(target + 1270);
    auto *t_1271 = buffer.data(target + 1271);
    auto *t_1272 = buffer.data(target + 1272);
    auto *t_1273 = buffer.data(target + 1273);
    auto *t_1274 = buffer.data(target + 1274);
    auto *t_1275 = buffer.data(target + 1275);
    auto *t_1276 = buffer.data(target + 1276);
    auto *t_1277 = buffer.data(target + 1277);
    auto *t_1278 = buffer.data(target + 1278);
    auto *t_1279 = buffer.data(target + 1279);
    auto *t_1280 = buffer.data(target + 1280);
    auto *t_1281 = buffer.data(target + 1281);
    auto *t_1282 = buffer.data(target + 1282);
    auto *t_1283 = buffer.data(target + 1283);
    auto *t_1284 = buffer.data(target + 1284);
    auto *t_1285 = buffer.data(target + 1285);
    auto *t_1286 = buffer.data(target + 1286);
    auto *t_1287 = buffer.data(target + 1287);
    auto *t_1288 = buffer.data(target + 1288);
    auto *t_1289 = buffer.data(target + 1289);
    auto *t_1290 = buffer.data(target + 1290);
    auto *t_1291 = buffer.data(target + 1291);
    auto *t_1292 = buffer.data(target + 1292);
    auto *t_1293 = buffer.data(target + 1293);
    auto *t_1294 = buffer.data(target + 1294);
    auto *t_1295 = buffer.data(target + 1295);

    const auto *ik_881 = buffer.data(ik + 881);
    const auto *ik_882 = buffer.data(ik + 882);
    const auto *ik_883 = buffer.data(ik + 883);
    const auto *ik_884 = buffer.data(ik + 884);
    const auto *ik_885 = buffer.data(ik + 885);
    const auto *ik_886 = buffer.data(ik + 886);
    const auto *ik_887 = buffer.data(ik + 887);
    const auto *ik_888 = buffer.data(ik + 888);
    const auto *ik_889 = buffer.data(ik + 889);
    const auto *ik_890 = buffer.data(ik + 890);
    const auto *ik_891 = buffer.data(ik + 891);
    const auto *ik_892 = buffer.data(ik + 892);
    const auto *ik_893 = buffer.data(ik + 893);
    const auto *ik_894 = buffer.data(ik + 894);
    const auto *ik_895 = buffer.data(ik + 895);
    const auto *ik_896 = buffer.data(ik + 896);
    const auto *ik_897 = buffer.data(ik + 897);
    const auto *ik_898 = buffer.data(ik + 898);
    const auto *ik_899 = buffer.data(ik + 899);
    const auto *ik_900 = buffer.data(ik + 900);
    const auto *ik_901 = buffer.data(ik + 901);
    const auto *ik_902 = buffer.data(ik + 902);
    const auto *ik_903 = buffer.data(ik + 903);
    const auto *ik_904 = buffer.data(ik + 904);
    const auto *ik_905 = buffer.data(ik + 905);
    const auto *ik_906 = buffer.data(ik + 906);
    const auto *ik_907 = buffer.data(ik + 907);
    const auto *ik_908 = buffer.data(ik + 908);
    const auto *ik_909 = buffer.data(ik + 909);
    const auto *ik_910 = buffer.data(ik + 910);
    const auto *ik_911 = buffer.data(ik + 911);
    const auto *ik_912 = buffer.data(ik + 912);
    const auto *ik_913 = buffer.data(ik + 913);
    const auto *ik_914 = buffer.data(ik + 914);
    const auto *ik_915 = buffer.data(ik + 915);
    const auto *ik_916 = buffer.data(ik + 916);
    const auto *ik_917 = buffer.data(ik + 917);
    const auto *ik_918 = buffer.data(ik + 918);
    const auto *ik_919 = buffer.data(ik + 919);
    const auto *ik_920 = buffer.data(ik + 920);
    const auto *ik_921 = buffer.data(ik + 921);
    const auto *ik_922 = buffer.data(ik + 922);
    const auto *ik_923 = buffer.data(ik + 923);
    const auto *ik_924 = buffer.data(ik + 924);
    const auto *ik_925 = buffer.data(ik + 925);
    const auto *ik_926 = buffer.data(ik + 926);
    const auto *ik_927 = buffer.data(ik + 927);
    const auto *ik_928 = buffer.data(ik + 928);
    const auto *ik_929 = buffer.data(ik + 929);
    const auto *ik_930 = buffer.data(ik + 930);
    const auto *ik_931 = buffer.data(ik + 931);
    const auto *ik_932 = buffer.data(ik + 932);
    const auto *ik_933 = buffer.data(ik + 933);
    const auto *ik_934 = buffer.data(ik + 934);
    const auto *ik_935 = buffer.data(ik + 935);
    const auto *ik_936 = buffer.data(ik + 936);
    const auto *ik_937 = buffer.data(ik + 937);
    const auto *ik_938 = buffer.data(ik + 938);
    const auto *ik_939 = buffer.data(ik + 939);
    const auto *ik_940 = buffer.data(ik + 940);
    const auto *ik_941 = buffer.data(ik + 941);
    const auto *ik_942 = buffer.data(ik + 942);
    const auto *ik_943 = buffer.data(ik + 943);
    const auto *ik_944 = buffer.data(ik + 944);
    const auto *ik_945 = buffer.data(ik + 945);
    const auto *ik_946 = buffer.data(ik + 946);
    const auto *ik_947 = buffer.data(ik + 947);
    const auto *ik_948 = buffer.data(ik + 948);
    const auto *ik_949 = buffer.data(ik + 949);
    const auto *ik_950 = buffer.data(ik + 950);
    const auto *ik_951 = buffer.data(ik + 951);
    const auto *ik_952 = buffer.data(ik + 952);
    const auto *ik_953 = buffer.data(ik + 953);
    const auto *ik_954 = buffer.data(ik + 954);
    const auto *ik_955 = buffer.data(ik + 955);
    const auto *ik_956 = buffer.data(ik + 956);
    const auto *ik_957 = buffer.data(ik + 957);
    const auto *ik_958 = buffer.data(ik + 958);
    const auto *ik_959 = buffer.data(ik + 959);
    const auto *ik_960 = buffer.data(ik + 960);
    const auto *ik_961 = buffer.data(ik + 961);
    const auto *ik_962 = buffer.data(ik + 962);
    const auto *ik_963 = buffer.data(ik + 963);
    const auto *ik_964 = buffer.data(ik + 964);
    const auto *ik_965 = buffer.data(ik + 965);
    const auto *ik_966 = buffer.data(ik + 966);
    const auto *ik_967 = buffer.data(ik + 967);
    const auto *ik_968 = buffer.data(ik + 968);
    const auto *ik_969 = buffer.data(ik + 969);
    const auto *ik_970 = buffer.data(ik + 970);
    const auto *ik_971 = buffer.data(ik + 971);
    const auto *ik_972 = buffer.data(ik + 972);
    const auto *ik_973 = buffer.data(ik + 973);
    const auto *ik_974 = buffer.data(ik + 974);
    const auto *ik_975 = buffer.data(ik + 975);
    const auto *ik_976 = buffer.data(ik + 976);
    const auto *ik_977 = buffer.data(ik + 977);
    const auto *ik_978 = buffer.data(ik + 978);
    const auto *ik_979 = buffer.data(ik + 979);
    const auto *ik_980 = buffer.data(ik + 980);
    const auto *ik_981 = buffer.data(ik + 981);
    const auto *ik_982 = buffer.data(ik + 982);
    const auto *ik_983 = buffer.data(ik + 983);
    const auto *ik_984 = buffer.data(ik + 984);
    const auto *ik_985 = buffer.data(ik + 985);
    const auto *ik_986 = buffer.data(ik + 986);
    const auto *ik_987 = buffer.data(ik + 987);
    const auto *ik_988 = buffer.data(ik + 988);
    const auto *ik_989 = buffer.data(ik + 989);
    const auto *ik_990 = buffer.data(ik + 990);
    const auto *ik_991 = buffer.data(ik + 991);
    const auto *ik_992 = buffer.data(ik + 992);
    const auto *ik_993 = buffer.data(ik + 993);
    const auto *ik_994 = buffer.data(ik + 994);
    const auto *ik_995 = buffer.data(ik + 995);
    const auto *ik_996 = buffer.data(ik + 996);
    const auto *ik_997 = buffer.data(ik + 997);
    const auto *ik_998 = buffer.data(ik + 998);
    const auto *ik_999 = buffer.data(ik + 999);
    const auto *ik_1000 = buffer.data(ik + 1000);
    const auto *ik_1001 = buffer.data(ik + 1001);
    const auto *ik_1002 = buffer.data(ik + 1002);
    const auto *ik_1003 = buffer.data(ik + 1003);
    const auto *ik_1004 = buffer.data(ik + 1004);
    const auto *ik_1005 = buffer.data(ik + 1005);
    const auto *ik_1006 = buffer.data(ik + 1006);
    const auto *ik_1007 = buffer.data(ik + 1007);

    const auto *lk_1421 = buffer.data(lk + 1421);
    const auto *lk_1422 = buffer.data(lk + 1422);
    const auto *lk_1423 = buffer.data(lk + 1423);
    const auto *lk_1424 = buffer.data(lk + 1424);
    const auto *lk_1425 = buffer.data(lk + 1425);
    const auto *lk_1426 = buffer.data(lk + 1426);
    const auto *lk_1427 = buffer.data(lk + 1427);
    const auto *lk_1428 = buffer.data(lk + 1428);
    const auto *lk_1429 = buffer.data(lk + 1429);
    const auto *lk_1430 = buffer.data(lk + 1430);
    const auto *lk_1431 = buffer.data(lk + 1431);
    const auto *lk_1432 = buffer.data(lk + 1432);
    const auto *lk_1433 = buffer.data(lk + 1433);
    const auto *lk_1434 = buffer.data(lk + 1434);
    const auto *lk_1435 = buffer.data(lk + 1435);
    const auto *lk_1436 = buffer.data(lk + 1436);
    const auto *lk_1437 = buffer.data(lk + 1437);
    const auto *lk_1438 = buffer.data(lk + 1438);
    const auto *lk_1439 = buffer.data(lk + 1439);
    const auto *lk_1440 = buffer.data(lk + 1440);
    const auto *lk_1441 = buffer.data(lk + 1441);
    const auto *lk_1442 = buffer.data(lk + 1442);
    const auto *lk_1443 = buffer.data(lk + 1443);
    const auto *lk_1444 = buffer.data(lk + 1444);
    const auto *lk_1445 = buffer.data(lk + 1445);
    const auto *lk_1446 = buffer.data(lk + 1446);
    const auto *lk_1447 = buffer.data(lk + 1447);
    const auto *lk_1448 = buffer.data(lk + 1448);
    const auto *lk_1449 = buffer.data(lk + 1449);
    const auto *lk_1450 = buffer.data(lk + 1450);
    const auto *lk_1451 = buffer.data(lk + 1451);
    const auto *lk_1452 = buffer.data(lk + 1452);
    const auto *lk_1453 = buffer.data(lk + 1453);
    const auto *lk_1454 = buffer.data(lk + 1454);
    const auto *lk_1455 = buffer.data(lk + 1455);
    const auto *lk_1456 = buffer.data(lk + 1456);
    const auto *lk_1457 = buffer.data(lk + 1457);
    const auto *lk_1458 = buffer.data(lk + 1458);
    const auto *lk_1459 = buffer.data(lk + 1459);
    const auto *lk_1460 = buffer.data(lk + 1460);
    const auto *lk_1461 = buffer.data(lk + 1461);
    const auto *lk_1462 = buffer.data(lk + 1462);
    const auto *lk_1463 = buffer.data(lk + 1463);
    const auto *lk_1464 = buffer.data(lk + 1464);
    const auto *lk_1465 = buffer.data(lk + 1465);
    const auto *lk_1466 = buffer.data(lk + 1466);
    const auto *lk_1467 = buffer.data(lk + 1467);
    const auto *lk_1468 = buffer.data(lk + 1468);
    const auto *lk_1469 = buffer.data(lk + 1469);
    const auto *lk_1470 = buffer.data(lk + 1470);
    const auto *lk_1471 = buffer.data(lk + 1471);
    const auto *lk_1472 = buffer.data(lk + 1472);
    const auto *lk_1473 = buffer.data(lk + 1473);
    const auto *lk_1474 = buffer.data(lk + 1474);
    const auto *lk_1475 = buffer.data(lk + 1475);
    const auto *lk_1476 = buffer.data(lk + 1476);
    const auto *lk_1477 = buffer.data(lk + 1477);
    const auto *lk_1478 = buffer.data(lk + 1478);
    const auto *lk_1479 = buffer.data(lk + 1479);
    const auto *lk_1480 = buffer.data(lk + 1480);
    const auto *lk_1481 = buffer.data(lk + 1481);
    const auto *lk_1482 = buffer.data(lk + 1482);
    const auto *lk_1483 = buffer.data(lk + 1483);
    const auto *lk_1484 = buffer.data(lk + 1484);
    const auto *lk_1485 = buffer.data(lk + 1485);
    const auto *lk_1486 = buffer.data(lk + 1486);
    const auto *lk_1487 = buffer.data(lk + 1487);
    const auto *lk_1488 = buffer.data(lk + 1488);
    const auto *lk_1489 = buffer.data(lk + 1489);
    const auto *lk_1490 = buffer.data(lk + 1490);
    const auto *lk_1491 = buffer.data(lk + 1491);
    const auto *lk_1492 = buffer.data(lk + 1492);
    const auto *lk_1493 = buffer.data(lk + 1493);
    const auto *lk_1494 = buffer.data(lk + 1494);
    const auto *lk_1495 = buffer.data(lk + 1495);
    const auto *lk_1496 = buffer.data(lk + 1496);
    const auto *lk_1497 = buffer.data(lk + 1497);
    const auto *lk_1498 = buffer.data(lk + 1498);
    const auto *lk_1499 = buffer.data(lk + 1499);
    const auto *lk_1500 = buffer.data(lk + 1500);
    const auto *lk_1501 = buffer.data(lk + 1501);
    const auto *lk_1502 = buffer.data(lk + 1502);
    const auto *lk_1503 = buffer.data(lk + 1503);
    const auto *lk_1504 = buffer.data(lk + 1504);
    const auto *lk_1505 = buffer.data(lk + 1505);
    const auto *lk_1506 = buffer.data(lk + 1506);
    const auto *lk_1507 = buffer.data(lk + 1507);
    const auto *lk_1508 = buffer.data(lk + 1508);
    const auto *lk_1509 = buffer.data(lk + 1509);
    const auto *lk_1510 = buffer.data(lk + 1510);
    const auto *lk_1511 = buffer.data(lk + 1511);
    const auto *lk_1512 = buffer.data(lk + 1512);
    const auto *lk_1513 = buffer.data(lk + 1513);
    const auto *lk_1514 = buffer.data(lk + 1514);
    const auto *lk_1515 = buffer.data(lk + 1515);
    const auto *lk_1516 = buffer.data(lk + 1516);
    const auto *lk_1517 = buffer.data(lk + 1517);
    const auto *lk_1518 = buffer.data(lk + 1518);
    const auto *lk_1519 = buffer.data(lk + 1519);
    const auto *lk_1520 = buffer.data(lk + 1520);
    const auto *lk_1521 = buffer.data(lk + 1521);
    const auto *lk_1522 = buffer.data(lk + 1522);
    const auto *lk_1523 = buffer.data(lk + 1523);
    const auto *lk_1524 = buffer.data(lk + 1524);
    const auto *lk_1525 = buffer.data(lk + 1525);
    const auto *lk_1526 = buffer.data(lk + 1526);
    const auto *lk_1527 = buffer.data(lk + 1527);
    const auto *lk_1528 = buffer.data(lk + 1528);
    const auto *lk_1529 = buffer.data(lk + 1529);
    const auto *lk_1530 = buffer.data(lk + 1530);
    const auto *lk_1531 = buffer.data(lk + 1531);
    const auto *lk_1532 = buffer.data(lk + 1532);
    const auto *lk_1533 = buffer.data(lk + 1533);
    const auto *lk_1534 = buffer.data(lk + 1534);
    const auto *lk_1535 = buffer.data(lk + 1535);
    const auto *lk_1536 = buffer.data(lk + 1536);
    const auto *lk_1537 = buffer.data(lk + 1537);
    const auto *lk_1538 = buffer.data(lk + 1538);
    const auto *lk_1539 = buffer.data(lk + 1539);
    const auto *lk_1540 = buffer.data(lk + 1540);
    const auto *lk_1541 = buffer.data(lk + 1541);
    const auto *lk_1542 = buffer.data(lk + 1542);
    const auto *lk_1543 = buffer.data(lk + 1543);
    const auto *lk_1544 = buffer.data(lk + 1544);
    const auto *lk_1545 = buffer.data(lk + 1545);
    const auto *lk_1546 = buffer.data(lk + 1546);
    const auto *lk_1547 = buffer.data(lk + 1547);
    const auto *lk_1548 = buffer.data(lk + 1548);
    const auto *lk_1549 = buffer.data(lk + 1549);
    const auto *lk_1550 = buffer.data(lk + 1550);
    const auto *lk_1551 = buffer.data(lk + 1551);
    const auto *lk_1552 = buffer.data(lk + 1552);
    const auto *lk_1553 = buffer.data(lk + 1553);
    const auto *lk_1554 = buffer.data(lk + 1554);
    const auto *lk_1555 = buffer.data(lk + 1555);
    const auto *lk_1556 = buffer.data(lk + 1556);
    const auto *lk_1557 = buffer.data(lk + 1557);
    const auto *lk_1558 = buffer.data(lk + 1558);
    const auto *lk_1559 = buffer.data(lk + 1559);
    const auto *lk_1560 = buffer.data(lk + 1560);
    const auto *lk_1561 = buffer.data(lk + 1561);
    const auto *lk_1562 = buffer.data(lk + 1562);
    const auto *lk_1563 = buffer.data(lk + 1563);
    const auto *lk_1564 = buffer.data(lk + 1564);
    const auto *lk_1565 = buffer.data(lk + 1565);
    const auto *lk_1566 = buffer.data(lk + 1566);
    const auto *lk_1567 = buffer.data(lk + 1567);
    const auto *lk_1568 = buffer.data(lk + 1568);
    const auto *lk_1569 = buffer.data(lk + 1569);
    const auto *lk_1570 = buffer.data(lk + 1570);
    const auto *lk_1571 = buffer.data(lk + 1571);
    const auto *lk_1572 = buffer.data(lk + 1572);
    const auto *lk_1573 = buffer.data(lk + 1573);
    const auto *lk_1574 = buffer.data(lk + 1574);
    const auto *lk_1575 = buffer.data(lk + 1575);
    const auto *lk_1576 = buffer.data(lk + 1576);
    const auto *lk_1577 = buffer.data(lk + 1577);
    const auto *lk_1578 = buffer.data(lk + 1578);
    const auto *lk_1579 = buffer.data(lk + 1579);
    const auto *lk_1580 = buffer.data(lk + 1580);
    const auto *lk_1581 = buffer.data(lk + 1581);
    const auto *lk_1582 = buffer.data(lk + 1582);
    const auto *lk_1583 = buffer.data(lk + 1583);

#pragma omp simd aligned(t_1133, t_1134, t_1135, t_1136, t_1137, ik_881, ik_882, ik_883, \
                         ik_884, ik_885, lk_1421, lk_1422, lk_1423, lk_1424, \
                         lk_1425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = -4.0 * ik_881[k]
                    + f_0 * lk_1421[k];

        t_1134[k] = -4.0 * ik_882[k]
                    + f_0 * lk_1422[k];

        t_1135[k] = -4.0 * ik_883[k]
                    + f_0 * lk_1423[k];

        t_1136[k] = -4.0 * ik_884[k]
                    + f_0 * lk_1424[k];

        t_1137[k] = -4.0 * ik_885[k]
                    + f_0 * lk_1425[k];
    }

#pragma omp simd aligned(t_1138, t_1139, t_1140, t_1141, t_1142, ik_886, ik_887, ik_888, \
                         ik_889, ik_890, lk_1426, lk_1427, lk_1428, lk_1429, \
                         lk_1430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1138[k] = -4.0 * ik_886[k]
                    + f_0 * lk_1426[k];

        t_1139[k] = -4.0 * ik_887[k]
                    + f_0 * lk_1427[k];

        t_1140[k] = -4.0 * ik_888[k]
                    + f_0 * lk_1428[k];

        t_1141[k] = -4.0 * ik_889[k]
                    + f_0 * lk_1429[k];

        t_1142[k] = -4.0 * ik_890[k]
                    + f_0 * lk_1430[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, t_1146, t_1147, ik_891, ik_892, ik_893, \
                         ik_894, ik_895, lk_1431, lk_1432, lk_1433, lk_1434, \
                         lk_1435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = -4.0 * ik_891[k]
                    + f_0 * lk_1431[k];

        t_1144[k] = -4.0 * ik_892[k]
                    + f_0 * lk_1432[k];

        t_1145[k] = -4.0 * ik_893[k]
                    + f_0 * lk_1433[k];

        t_1146[k] = -4.0 * ik_894[k]
                    + f_0 * lk_1434[k];

        t_1147[k] = -4.0 * ik_895[k]
                    + f_0 * lk_1435[k];
    }

#pragma omp simd aligned(t_1148, t_1149, t_1150, t_1151, t_1152, ik_896, ik_897, ik_898, \
                         ik_899, ik_900, lk_1436, lk_1437, lk_1438, lk_1439, \
                         lk_1440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1148[k] = -4.0 * ik_896[k]
                    + f_0 * lk_1436[k];

        t_1149[k] = -4.0 * ik_897[k]
                    + f_0 * lk_1437[k];

        t_1150[k] = -4.0 * ik_898[k]
                    + f_0 * lk_1438[k];

        t_1151[k] = -4.0 * ik_899[k]
                    + f_0 * lk_1439[k];

        t_1152[k] = -3.0 * ik_900[k]
                    + f_0 * lk_1440[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, t_1157, ik_901, ik_902, ik_903, \
                         ik_904, ik_905, lk_1441, lk_1442, lk_1443, lk_1444, \
                         lk_1445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = -3.0 * ik_901[k]
                    + f_0 * lk_1441[k];

        t_1154[k] = -3.0 * ik_902[k]
                    + f_0 * lk_1442[k];

        t_1155[k] = -3.0 * ik_903[k]
                    + f_0 * lk_1443[k];

        t_1156[k] = -3.0 * ik_904[k]
                    + f_0 * lk_1444[k];

        t_1157[k] = -3.0 * ik_905[k]
                    + f_0 * lk_1445[k];
    }

#pragma omp simd aligned(t_1158, t_1159, t_1160, t_1161, t_1162, ik_906, ik_907, ik_908, \
                         ik_909, ik_910, lk_1446, lk_1447, lk_1448, lk_1449, \
                         lk_1450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1158[k] = -3.0 * ik_906[k]
                    + f_0 * lk_1446[k];

        t_1159[k] = -3.0 * ik_907[k]
                    + f_0 * lk_1447[k];

        t_1160[k] = -3.0 * ik_908[k]
                    + f_0 * lk_1448[k];

        t_1161[k] = -3.0 * ik_909[k]
                    + f_0 * lk_1449[k];

        t_1162[k] = -3.0 * ik_910[k]
                    + f_0 * lk_1450[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, t_1166, t_1167, ik_911, ik_912, ik_913, \
                         ik_914, ik_915, lk_1451, lk_1452, lk_1453, lk_1454, \
                         lk_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = -3.0 * ik_911[k]
                    + f_0 * lk_1451[k];

        t_1164[k] = -3.0 * ik_912[k]
                    + f_0 * lk_1452[k];

        t_1165[k] = -3.0 * ik_913[k]
                    + f_0 * lk_1453[k];

        t_1166[k] = -3.0 * ik_914[k]
                    + f_0 * lk_1454[k];

        t_1167[k] = -3.0 * ik_915[k]
                    + f_0 * lk_1455[k];
    }

#pragma omp simd aligned(t_1168, t_1169, t_1170, t_1171, t_1172, ik_916, ik_917, ik_918, \
                         ik_919, ik_920, lk_1456, lk_1457, lk_1458, lk_1459, \
                         lk_1460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1168[k] = -3.0 * ik_916[k]
                    + f_0 * lk_1456[k];

        t_1169[k] = -3.0 * ik_917[k]
                    + f_0 * lk_1457[k];

        t_1170[k] = -3.0 * ik_918[k]
                    + f_0 * lk_1458[k];

        t_1171[k] = -3.0 * ik_919[k]
                    + f_0 * lk_1459[k];

        t_1172[k] = -3.0 * ik_920[k]
                    + f_0 * lk_1460[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, t_1176, t_1177, ik_921, ik_922, ik_923, \
                         ik_924, ik_925, lk_1461, lk_1462, lk_1463, lk_1464, \
                         lk_1465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = -3.0 * ik_921[k]
                    + f_0 * lk_1461[k];

        t_1174[k] = -3.0 * ik_922[k]
                    + f_0 * lk_1462[k];

        t_1175[k] = -3.0 * ik_923[k]
                    + f_0 * lk_1463[k];

        t_1176[k] = -3.0 * ik_924[k]
                    + f_0 * lk_1464[k];

        t_1177[k] = -3.0 * ik_925[k]
                    + f_0 * lk_1465[k];
    }

#pragma omp simd aligned(t_1178, t_1179, t_1180, t_1181, t_1182, ik_926, ik_927, ik_928, \
                         ik_929, ik_930, lk_1466, lk_1467, lk_1468, lk_1469, \
                         lk_1470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1178[k] = -3.0 * ik_926[k]
                    + f_0 * lk_1466[k];

        t_1179[k] = -3.0 * ik_927[k]
                    + f_0 * lk_1467[k];

        t_1180[k] = -3.0 * ik_928[k]
                    + f_0 * lk_1468[k];

        t_1181[k] = -3.0 * ik_929[k]
                    + f_0 * lk_1469[k];

        t_1182[k] = -3.0 * ik_930[k]
                    + f_0 * lk_1470[k];
    }

#pragma omp simd aligned(t_1183, t_1184, t_1185, t_1186, t_1187, ik_931, ik_932, ik_933, \
                         ik_934, ik_935, lk_1471, lk_1472, lk_1473, lk_1474, \
                         lk_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1183[k] = -3.0 * ik_931[k]
                    + f_0 * lk_1471[k];

        t_1184[k] = -3.0 * ik_932[k]
                    + f_0 * lk_1472[k];

        t_1185[k] = -3.0 * ik_933[k]
                    + f_0 * lk_1473[k];

        t_1186[k] = -3.0 * ik_934[k]
                    + f_0 * lk_1474[k];

        t_1187[k] = -3.0 * ik_935[k]
                    + f_0 * lk_1475[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, t_1191, t_1192, ik_936, ik_937, ik_938, \
                         ik_939, ik_940, lk_1476, lk_1477, lk_1478, lk_1479, \
                         lk_1480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = -2.0 * ik_936[k]
                    + f_0 * lk_1476[k];

        t_1189[k] = -2.0 * ik_937[k]
                    + f_0 * lk_1477[k];

        t_1190[k] = -2.0 * ik_938[k]
                    + f_0 * lk_1478[k];

        t_1191[k] = -2.0 * ik_939[k]
                    + f_0 * lk_1479[k];

        t_1192[k] = -2.0 * ik_940[k]
                    + f_0 * lk_1480[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, t_1196, t_1197, ik_941, ik_942, ik_943, \
                         ik_944, ik_945, lk_1481, lk_1482, lk_1483, lk_1484, \
                         lk_1485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = -2.0 * ik_941[k]
                    + f_0 * lk_1481[k];

        t_1194[k] = -2.0 * ik_942[k]
                    + f_0 * lk_1482[k];

        t_1195[k] = -2.0 * ik_943[k]
                    + f_0 * lk_1483[k];

        t_1196[k] = -2.0 * ik_944[k]
                    + f_0 * lk_1484[k];

        t_1197[k] = -2.0 * ik_945[k]
                    + f_0 * lk_1485[k];
    }

#pragma omp simd aligned(t_1198, t_1199, t_1200, t_1201, t_1202, ik_946, ik_947, ik_948, \
                         ik_949, ik_950, lk_1486, lk_1487, lk_1488, lk_1489, \
                         lk_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1198[k] = -2.0 * ik_946[k]
                    + f_0 * lk_1486[k];

        t_1199[k] = -2.0 * ik_947[k]
                    + f_0 * lk_1487[k];

        t_1200[k] = -2.0 * ik_948[k]
                    + f_0 * lk_1488[k];

        t_1201[k] = -2.0 * ik_949[k]
                    + f_0 * lk_1489[k];

        t_1202[k] = -2.0 * ik_950[k]
                    + f_0 * lk_1490[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, t_1206, t_1207, ik_951, ik_952, ik_953, \
                         ik_954, ik_955, lk_1491, lk_1492, lk_1493, lk_1494, \
                         lk_1495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = -2.0 * ik_951[k]
                    + f_0 * lk_1491[k];

        t_1204[k] = -2.0 * ik_952[k]
                    + f_0 * lk_1492[k];

        t_1205[k] = -2.0 * ik_953[k]
                    + f_0 * lk_1493[k];

        t_1206[k] = -2.0 * ik_954[k]
                    + f_0 * lk_1494[k];

        t_1207[k] = -2.0 * ik_955[k]
                    + f_0 * lk_1495[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, t_1211, t_1212, ik_956, ik_957, ik_958, \
                         ik_959, ik_960, lk_1496, lk_1497, lk_1498, lk_1499, \
                         lk_1500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = -2.0 * ik_956[k]
                    + f_0 * lk_1496[k];

        t_1209[k] = -2.0 * ik_957[k]
                    + f_0 * lk_1497[k];

        t_1210[k] = -2.0 * ik_958[k]
                    + f_0 * lk_1498[k];

        t_1211[k] = -2.0 * ik_959[k]
                    + f_0 * lk_1499[k];

        t_1212[k] = -2.0 * ik_960[k]
                    + f_0 * lk_1500[k];
    }

#pragma omp simd aligned(t_1213, t_1214, t_1215, t_1216, t_1217, ik_961, ik_962, ik_963, \
                         ik_964, ik_965, lk_1501, lk_1502, lk_1503, lk_1504, \
                         lk_1505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1213[k] = -2.0 * ik_961[k]
                    + f_0 * lk_1501[k];

        t_1214[k] = -2.0 * ik_962[k]
                    + f_0 * lk_1502[k];

        t_1215[k] = -2.0 * ik_963[k]
                    + f_0 * lk_1503[k];

        t_1216[k] = -2.0 * ik_964[k]
                    + f_0 * lk_1504[k];

        t_1217[k] = -2.0 * ik_965[k]
                    + f_0 * lk_1505[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, t_1221, t_1222, ik_966, ik_967, ik_968, \
                         ik_969, ik_970, lk_1506, lk_1507, lk_1508, lk_1509, \
                         lk_1510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = -2.0 * ik_966[k]
                    + f_0 * lk_1506[k];

        t_1219[k] = -2.0 * ik_967[k]
                    + f_0 * lk_1507[k];

        t_1220[k] = -2.0 * ik_968[k]
                    + f_0 * lk_1508[k];

        t_1221[k] = -2.0 * ik_969[k]
                    + f_0 * lk_1509[k];

        t_1222[k] = -2.0 * ik_970[k]
                    + f_0 * lk_1510[k];
    }

#pragma omp simd aligned(t_1223, t_1224, t_1225, t_1226, t_1227, ik_971, ik_972, ik_973, \
                         ik_974, ik_975, lk_1511, lk_1512, lk_1513, lk_1514, \
                         lk_1515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1223[k] = -2.0 * ik_971[k]
                    + f_0 * lk_1511[k];

        t_1224[k] = -ik_972[k]
                    + f_0 * lk_1512[k];

        t_1225[k] = -ik_973[k]
                    + f_0 * lk_1513[k];

        t_1226[k] = -ik_974[k]
                    + f_0 * lk_1514[k];

        t_1227[k] = -ik_975[k]
                    + f_0 * lk_1515[k];
    }

#pragma omp simd aligned(t_1228, t_1229, t_1230, t_1231, t_1232, ik_976, ik_977, ik_978, \
                         ik_979, ik_980, lk_1516, lk_1517, lk_1518, lk_1519, \
                         lk_1520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1228[k] = -ik_976[k]
                    + f_0 * lk_1516[k];

        t_1229[k] = -ik_977[k]
                    + f_0 * lk_1517[k];

        t_1230[k] = -ik_978[k]
                    + f_0 * lk_1518[k];

        t_1231[k] = -ik_979[k]
                    + f_0 * lk_1519[k];

        t_1232[k] = -ik_980[k]
                    + f_0 * lk_1520[k];
    }

#pragma omp simd aligned(t_1233, t_1234, t_1235, t_1236, t_1237, ik_981, ik_982, ik_983, \
                         ik_984, ik_985, lk_1521, lk_1522, lk_1523, lk_1524, \
                         lk_1525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1233[k] = -ik_981[k]
                    + f_0 * lk_1521[k];

        t_1234[k] = -ik_982[k]
                    + f_0 * lk_1522[k];

        t_1235[k] = -ik_983[k]
                    + f_0 * lk_1523[k];

        t_1236[k] = -ik_984[k]
                    + f_0 * lk_1524[k];

        t_1237[k] = -ik_985[k]
                    + f_0 * lk_1525[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, t_1241, t_1242, ik_986, ik_987, ik_988, \
                         ik_989, ik_990, lk_1526, lk_1527, lk_1528, lk_1529, \
                         lk_1530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = -ik_986[k]
                    + f_0 * lk_1526[k];

        t_1239[k] = -ik_987[k]
                    + f_0 * lk_1527[k];

        t_1240[k] = -ik_988[k]
                    + f_0 * lk_1528[k];

        t_1241[k] = -ik_989[k]
                    + f_0 * lk_1529[k];

        t_1242[k] = -ik_990[k]
                    + f_0 * lk_1530[k];
    }

#pragma omp simd aligned(t_1243, t_1244, t_1245, t_1246, t_1247, ik_991, ik_992, ik_993, \
                         ik_994, ik_995, lk_1531, lk_1532, lk_1533, lk_1534, \
                         lk_1535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1243[k] = -ik_991[k]
                    + f_0 * lk_1531[k];

        t_1244[k] = -ik_992[k]
                    + f_0 * lk_1532[k];

        t_1245[k] = -ik_993[k]
                    + f_0 * lk_1533[k];

        t_1246[k] = -ik_994[k]
                    + f_0 * lk_1534[k];

        t_1247[k] = -ik_995[k]
                    + f_0 * lk_1535[k];
    }

#pragma omp simd aligned(t_1248, t_1249, t_1250, t_1251, t_1252, ik_996, ik_997, ik_998, \
                         ik_999, ik_1000, lk_1536, lk_1537, lk_1538, lk_1539, \
                         lk_1540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1248[k] = -ik_996[k]
                    + f_0 * lk_1536[k];

        t_1249[k] = -ik_997[k]
                    + f_0 * lk_1537[k];

        t_1250[k] = -ik_998[k]
                    + f_0 * lk_1538[k];

        t_1251[k] = -ik_999[k]
                    + f_0 * lk_1539[k];

        t_1252[k] = -ik_1000[k]
                    + f_0 * lk_1540[k];
    }

#pragma omp simd aligned(t_1253, t_1254, t_1255, t_1256, t_1257, ik_1001, ik_1002, ik_1003, \
                         ik_1004, ik_1005, lk_1541, lk_1542, lk_1543, lk_1544, \
                         lk_1545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1253[k] = -ik_1001[k]
                    + f_0 * lk_1541[k];

        t_1254[k] = -ik_1002[k]
                    + f_0 * lk_1542[k];

        t_1255[k] = -ik_1003[k]
                    + f_0 * lk_1543[k];

        t_1256[k] = -ik_1004[k]
                    + f_0 * lk_1544[k];

        t_1257[k] = -ik_1005[k]
                    + f_0 * lk_1545[k];
    }

#pragma omp simd aligned(t_1258, t_1259, t_1260, t_1261, t_1262, t_1263, t_1264, ik_1006, \
                         ik_1007, lk_1546, lk_1547, lk_1548, lk_1549, lk_1550, lk_1551, \
                         lk_1552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1258[k] = -ik_1006[k]
                    + f_0 * lk_1546[k];

        t_1259[k] = -ik_1007[k]
                    + f_0 * lk_1547[k];

        t_1260[k] = f_0 * lk_1548[k];

        t_1261[k] = f_0 * lk_1549[k];

        t_1262[k] = f_0 * lk_1550[k];

        t_1263[k] = f_0 * lk_1551[k];

        t_1264[k] = f_0 * lk_1552[k];
    }

#pragma omp simd aligned(t_1265, t_1266, t_1267, t_1268, t_1269, t_1270, t_1271, t_1272, \
                         lk_1553, lk_1554, lk_1555, lk_1556, lk_1557, lk_1558, lk_1559, \
                         lk_1560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1265[k] = f_0 * lk_1553[k];

        t_1266[k] = f_0 * lk_1554[k];

        t_1267[k] = f_0 * lk_1555[k];

        t_1268[k] = f_0 * lk_1556[k];

        t_1269[k] = f_0 * lk_1557[k];

        t_1270[k] = f_0 * lk_1558[k];

        t_1271[k] = f_0 * lk_1559[k];

        t_1272[k] = f_0 * lk_1560[k];
    }

#pragma omp simd aligned(t_1273, t_1274, t_1275, t_1276, t_1277, t_1278, t_1279, t_1280, \
                         lk_1561, lk_1562, lk_1563, lk_1564, lk_1565, lk_1566, lk_1567, \
                         lk_1568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1273[k] = f_0 * lk_1561[k];

        t_1274[k] = f_0 * lk_1562[k];

        t_1275[k] = f_0 * lk_1563[k];

        t_1276[k] = f_0 * lk_1564[k];

        t_1277[k] = f_0 * lk_1565[k];

        t_1278[k] = f_0 * lk_1566[k];

        t_1279[k] = f_0 * lk_1567[k];

        t_1280[k] = f_0 * lk_1568[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, t_1284, t_1285, t_1286, t_1287, t_1288, \
                         lk_1569, lk_1570, lk_1571, lk_1572, lk_1573, lk_1574, lk_1575, \
                         lk_1576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = f_0 * lk_1569[k];

        t_1282[k] = f_0 * lk_1570[k];

        t_1283[k] = f_0 * lk_1571[k];

        t_1284[k] = f_0 * lk_1572[k];

        t_1285[k] = f_0 * lk_1573[k];

        t_1286[k] = f_0 * lk_1574[k];

        t_1287[k] = f_0 * lk_1575[k];

        t_1288[k] = f_0 * lk_1576[k];
    }

#pragma omp simd aligned(t_1289, t_1290, t_1291, t_1292, t_1293, t_1294, t_1295, lk_1577, \
                         lk_1578, lk_1579, lk_1580, lk_1581, lk_1582, \
                         lk_1583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1289[k] = f_0 * lk_1577[k];

        t_1290[k] = f_0 * lk_1578[k];

        t_1291[k] = f_0 * lk_1579[k];

        t_1292[k] = f_0 * lk_1580[k];

        t_1293[k] = f_0 * lk_1581[k];

        t_1294[k] = f_0 * lk_1582[k];

        t_1295[k] = f_0 * lk_1583[k];
    }
}

auto
compute_prim_geom_10_kk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t ik, const size_t lk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kk_electron_repulsion_1_piece0(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_1_piece1(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_1_piece2(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_1_piece3(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_1_piece4(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_1_piece5(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_1_piece6(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_1_piece7(buffer, target, ik, lk, ncols, alpha);
}

static auto
compute_prim_geom_10_kk_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_0 = buffer.data(ik + 0);
    const auto *ik_1 = buffer.data(ik + 1);
    const auto *ik_2 = buffer.data(ik + 2);
    const auto *ik_3 = buffer.data(ik + 3);
    const auto *ik_4 = buffer.data(ik + 4);
    const auto *ik_5 = buffer.data(ik + 5);
    const auto *ik_6 = buffer.data(ik + 6);
    const auto *ik_7 = buffer.data(ik + 7);
    const auto *ik_8 = buffer.data(ik + 8);
    const auto *ik_9 = buffer.data(ik + 9);
    const auto *ik_10 = buffer.data(ik + 10);
    const auto *ik_11 = buffer.data(ik + 11);
    const auto *ik_12 = buffer.data(ik + 12);
    const auto *ik_13 = buffer.data(ik + 13);
    const auto *ik_14 = buffer.data(ik + 14);
    const auto *ik_15 = buffer.data(ik + 15);
    const auto *ik_16 = buffer.data(ik + 16);
    const auto *ik_17 = buffer.data(ik + 17);
    const auto *ik_18 = buffer.data(ik + 18);
    const auto *ik_19 = buffer.data(ik + 19);
    const auto *ik_20 = buffer.data(ik + 20);
    const auto *ik_21 = buffer.data(ik + 21);
    const auto *ik_22 = buffer.data(ik + 22);
    const auto *ik_23 = buffer.data(ik + 23);
    const auto *ik_24 = buffer.data(ik + 24);
    const auto *ik_25 = buffer.data(ik + 25);
    const auto *ik_26 = buffer.data(ik + 26);
    const auto *ik_27 = buffer.data(ik + 27);
    const auto *ik_28 = buffer.data(ik + 28);
    const auto *ik_29 = buffer.data(ik + 29);
    const auto *ik_30 = buffer.data(ik + 30);
    const auto *ik_31 = buffer.data(ik + 31);
    const auto *ik_32 = buffer.data(ik + 32);
    const auto *ik_33 = buffer.data(ik + 33);
    const auto *ik_34 = buffer.data(ik + 34);
    const auto *ik_35 = buffer.data(ik + 35);
    const auto *ik_36 = buffer.data(ik + 36);
    const auto *ik_37 = buffer.data(ik + 37);
    const auto *ik_38 = buffer.data(ik + 38);
    const auto *ik_39 = buffer.data(ik + 39);
    const auto *ik_40 = buffer.data(ik + 40);
    const auto *ik_41 = buffer.data(ik + 41);
    const auto *ik_42 = buffer.data(ik + 42);
    const auto *ik_43 = buffer.data(ik + 43);
    const auto *ik_44 = buffer.data(ik + 44);
    const auto *ik_45 = buffer.data(ik + 45);
    const auto *ik_46 = buffer.data(ik + 46);
    const auto *ik_47 = buffer.data(ik + 47);
    const auto *ik_48 = buffer.data(ik + 48);
    const auto *ik_49 = buffer.data(ik + 49);
    const auto *ik_50 = buffer.data(ik + 50);
    const auto *ik_51 = buffer.data(ik + 51);
    const auto *ik_52 = buffer.data(ik + 52);
    const auto *ik_53 = buffer.data(ik + 53);
    const auto *ik_54 = buffer.data(ik + 54);
    const auto *ik_55 = buffer.data(ik + 55);
    const auto *ik_56 = buffer.data(ik + 56);
    const auto *ik_57 = buffer.data(ik + 57);
    const auto *ik_58 = buffer.data(ik + 58);
    const auto *ik_59 = buffer.data(ik + 59);
    const auto *ik_60 = buffer.data(ik + 60);
    const auto *ik_61 = buffer.data(ik + 61);
    const auto *ik_62 = buffer.data(ik + 62);
    const auto *ik_63 = buffer.data(ik + 63);
    const auto *ik_64 = buffer.data(ik + 64);
    const auto *ik_65 = buffer.data(ik + 65);
    const auto *ik_66 = buffer.data(ik + 66);
    const auto *ik_67 = buffer.data(ik + 67);
    const auto *ik_68 = buffer.data(ik + 68);
    const auto *ik_69 = buffer.data(ik + 69);
    const auto *ik_70 = buffer.data(ik + 70);
    const auto *ik_71 = buffer.data(ik + 71);
    const auto *ik_72 = buffer.data(ik + 72);
    const auto *ik_73 = buffer.data(ik + 73);
    const auto *ik_74 = buffer.data(ik + 74);
    const auto *ik_75 = buffer.data(ik + 75);
    const auto *ik_76 = buffer.data(ik + 76);

    const auto *lk_72 = buffer.data(lk + 72);
    const auto *lk_73 = buffer.data(lk + 73);
    const auto *lk_74 = buffer.data(lk + 74);
    const auto *lk_75 = buffer.data(lk + 75);
    const auto *lk_76 = buffer.data(lk + 76);
    const auto *lk_77 = buffer.data(lk + 77);
    const auto *lk_78 = buffer.data(lk + 78);
    const auto *lk_79 = buffer.data(lk + 79);
    const auto *lk_80 = buffer.data(lk + 80);
    const auto *lk_81 = buffer.data(lk + 81);
    const auto *lk_82 = buffer.data(lk + 82);
    const auto *lk_83 = buffer.data(lk + 83);
    const auto *lk_84 = buffer.data(lk + 84);
    const auto *lk_85 = buffer.data(lk + 85);
    const auto *lk_86 = buffer.data(lk + 86);
    const auto *lk_87 = buffer.data(lk + 87);
    const auto *lk_88 = buffer.data(lk + 88);
    const auto *lk_89 = buffer.data(lk + 89);
    const auto *lk_90 = buffer.data(lk + 90);
    const auto *lk_91 = buffer.data(lk + 91);
    const auto *lk_92 = buffer.data(lk + 92);
    const auto *lk_93 = buffer.data(lk + 93);
    const auto *lk_94 = buffer.data(lk + 94);
    const auto *lk_95 = buffer.data(lk + 95);
    const auto *lk_96 = buffer.data(lk + 96);
    const auto *lk_97 = buffer.data(lk + 97);
    const auto *lk_98 = buffer.data(lk + 98);
    const auto *lk_99 = buffer.data(lk + 99);
    const auto *lk_100 = buffer.data(lk + 100);
    const auto *lk_101 = buffer.data(lk + 101);
    const auto *lk_102 = buffer.data(lk + 102);
    const auto *lk_103 = buffer.data(lk + 103);
    const auto *lk_104 = buffer.data(lk + 104);
    const auto *lk_105 = buffer.data(lk + 105);
    const auto *lk_106 = buffer.data(lk + 106);
    const auto *lk_107 = buffer.data(lk + 107);
    const auto *lk_144 = buffer.data(lk + 144);
    const auto *lk_145 = buffer.data(lk + 145);
    const auto *lk_146 = buffer.data(lk + 146);
    const auto *lk_147 = buffer.data(lk + 147);
    const auto *lk_148 = buffer.data(lk + 148);
    const auto *lk_149 = buffer.data(lk + 149);
    const auto *lk_150 = buffer.data(lk + 150);
    const auto *lk_151 = buffer.data(lk + 151);
    const auto *lk_152 = buffer.data(lk + 152);
    const auto *lk_153 = buffer.data(lk + 153);
    const auto *lk_154 = buffer.data(lk + 154);
    const auto *lk_155 = buffer.data(lk + 155);
    const auto *lk_156 = buffer.data(lk + 156);
    const auto *lk_157 = buffer.data(lk + 157);
    const auto *lk_158 = buffer.data(lk + 158);
    const auto *lk_159 = buffer.data(lk + 159);
    const auto *lk_160 = buffer.data(lk + 160);
    const auto *lk_161 = buffer.data(lk + 161);
    const auto *lk_162 = buffer.data(lk + 162);
    const auto *lk_163 = buffer.data(lk + 163);
    const auto *lk_164 = buffer.data(lk + 164);
    const auto *lk_165 = buffer.data(lk + 165);
    const auto *lk_166 = buffer.data(lk + 166);
    const auto *lk_167 = buffer.data(lk + 167);
    const auto *lk_168 = buffer.data(lk + 168);
    const auto *lk_169 = buffer.data(lk + 169);
    const auto *lk_170 = buffer.data(lk + 170);
    const auto *lk_171 = buffer.data(lk + 171);
    const auto *lk_172 = buffer.data(lk + 172);
    const auto *lk_173 = buffer.data(lk + 173);
    const auto *lk_174 = buffer.data(lk + 174);
    const auto *lk_175 = buffer.data(lk + 175);
    const auto *lk_176 = buffer.data(lk + 176);
    const auto *lk_177 = buffer.data(lk + 177);
    const auto *lk_178 = buffer.data(lk + 178);
    const auto *lk_179 = buffer.data(lk + 179);
    const auto *lk_180 = buffer.data(lk + 180);
    const auto *lk_181 = buffer.data(lk + 181);
    const auto *lk_182 = buffer.data(lk + 182);
    const auto *lk_183 = buffer.data(lk + 183);
    const auto *lk_184 = buffer.data(lk + 184);
    const auto *lk_185 = buffer.data(lk + 185);
    const auto *lk_186 = buffer.data(lk + 186);
    const auto *lk_187 = buffer.data(lk + 187);
    const auto *lk_188 = buffer.data(lk + 188);
    const auto *lk_189 = buffer.data(lk + 189);
    const auto *lk_190 = buffer.data(lk + 190);
    const auto *lk_191 = buffer.data(lk + 191);
    const auto *lk_192 = buffer.data(lk + 192);
    const auto *lk_193 = buffer.data(lk + 193);
    const auto *lk_194 = buffer.data(lk + 194);
    const auto *lk_195 = buffer.data(lk + 195);
    const auto *lk_196 = buffer.data(lk + 196);
    const auto *lk_197 = buffer.data(lk + 197);
    const auto *lk_198 = buffer.data(lk + 198);
    const auto *lk_199 = buffer.data(lk + 199);
    const auto *lk_200 = buffer.data(lk + 200);
    const auto *lk_201 = buffer.data(lk + 201);
    const auto *lk_202 = buffer.data(lk + 202);
    const auto *lk_203 = buffer.data(lk + 203);
    const auto *lk_204 = buffer.data(lk + 204);
    const auto *lk_205 = buffer.data(lk + 205);
    const auto *lk_206 = buffer.data(lk + 206);
    const auto *lk_207 = buffer.data(lk + 207);
    const auto *lk_208 = buffer.data(lk + 208);
    const auto *lk_209 = buffer.data(lk + 209);
    const auto *lk_210 = buffer.data(lk + 210);
    const auto *lk_211 = buffer.data(lk + 211);
    const auto *lk_212 = buffer.data(lk + 212);
    const auto *lk_213 = buffer.data(lk + 213);
    const auto *lk_214 = buffer.data(lk + 214);
    const auto *lk_215 = buffer.data(lk + 215);
    const auto *lk_252 = buffer.data(lk + 252);
    const auto *lk_253 = buffer.data(lk + 253);
    const auto *lk_254 = buffer.data(lk + 254);
    const auto *lk_255 = buffer.data(lk + 255);
    const auto *lk_256 = buffer.data(lk + 256);
    const auto *lk_257 = buffer.data(lk + 257);
    const auto *lk_258 = buffer.data(lk + 258);
    const auto *lk_259 = buffer.data(lk + 259);
    const auto *lk_260 = buffer.data(lk + 260);
    const auto *lk_261 = buffer.data(lk + 261);
    const auto *lk_262 = buffer.data(lk + 262);
    const auto *lk_263 = buffer.data(lk + 263);
    const auto *lk_264 = buffer.data(lk + 264);
    const auto *lk_265 = buffer.data(lk + 265);
    const auto *lk_266 = buffer.data(lk + 266);
    const auto *lk_267 = buffer.data(lk + 267);
    const auto *lk_268 = buffer.data(lk + 268);
    const auto *lk_269 = buffer.data(lk + 269);
    const auto *lk_270 = buffer.data(lk + 270);
    const auto *lk_271 = buffer.data(lk + 271);
    const auto *lk_272 = buffer.data(lk + 272);
    const auto *lk_273 = buffer.data(lk + 273);
    const auto *lk_274 = buffer.data(lk + 274);
    const auto *lk_275 = buffer.data(lk + 275);
    const auto *lk_276 = buffer.data(lk + 276);
    const auto *lk_277 = buffer.data(lk + 277);
    const auto *lk_278 = buffer.data(lk + 278);
    const auto *lk_279 = buffer.data(lk + 279);
    const auto *lk_280 = buffer.data(lk + 280);
    const auto *lk_281 = buffer.data(lk + 281);
    const auto *lk_282 = buffer.data(lk + 282);
    const auto *lk_283 = buffer.data(lk + 283);
    const auto *lk_284 = buffer.data(lk + 284);
    const auto *lk_285 = buffer.data(lk + 285);
    const auto *lk_286 = buffer.data(lk + 286);
    const auto *lk_287 = buffer.data(lk + 287);
    const auto *lk_288 = buffer.data(lk + 288);
    const auto *lk_289 = buffer.data(lk + 289);
    const auto *lk_290 = buffer.data(lk + 290);
    const auto *lk_291 = buffer.data(lk + 291);
    const auto *lk_292 = buffer.data(lk + 292);
    const auto *lk_293 = buffer.data(lk + 293);
    const auto *lk_294 = buffer.data(lk + 294);
    const auto *lk_295 = buffer.data(lk + 295);
    const auto *lk_296 = buffer.data(lk + 296);
    const auto *lk_297 = buffer.data(lk + 297);
    const auto *lk_298 = buffer.data(lk + 298);
    const auto *lk_299 = buffer.data(lk + 299);
    const auto *lk_300 = buffer.data(lk + 300);
    const auto *lk_301 = buffer.data(lk + 301);
    const auto *lk_302 = buffer.data(lk + 302);
    const auto *lk_303 = buffer.data(lk + 303);
    const auto *lk_304 = buffer.data(lk + 304);
    const auto *lk_305 = buffer.data(lk + 305);
    const auto *lk_306 = buffer.data(lk + 306);
    const auto *lk_307 = buffer.data(lk + 307);
    const auto *lk_308 = buffer.data(lk + 308);
    const auto *lk_309 = buffer.data(lk + 309);
    const auto *lk_310 = buffer.data(lk + 310);
    const auto *lk_311 = buffer.data(lk + 311);
    const auto *lk_312 = buffer.data(lk + 312);
    const auto *lk_313 = buffer.data(lk + 313);
    const auto *lk_314 = buffer.data(lk + 314);
    const auto *lk_315 = buffer.data(lk + 315);
    const auto *lk_316 = buffer.data(lk + 316);
    const auto *lk_317 = buffer.data(lk + 317);
    const auto *lk_318 = buffer.data(lk + 318);
    const auto *lk_319 = buffer.data(lk + 319);
    const auto *lk_320 = buffer.data(lk + 320);
    const auto *lk_321 = buffer.data(lk + 321);
    const auto *lk_322 = buffer.data(lk + 322);
    const auto *lk_323 = buffer.data(lk + 323);
    const auto *lk_324 = buffer.data(lk + 324);
    const auto *lk_325 = buffer.data(lk + 325);
    const auto *lk_326 = buffer.data(lk + 326);
    const auto *lk_327 = buffer.data(lk + 327);
    const auto *lk_328 = buffer.data(lk + 328);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, lk_72, lk_73, lk_74, lk_75, \
                         lk_76, lk_77, lk_78, lk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lk_72[k];

        t_1[k] = f_0 * lk_73[k];

        t_2[k] = f_0 * lk_74[k];

        t_3[k] = f_0 * lk_75[k];

        t_4[k] = f_0 * lk_76[k];

        t_5[k] = f_0 * lk_77[k];

        t_6[k] = f_0 * lk_78[k];

        t_7[k] = f_0 * lk_79[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, lk_80, lk_81, lk_82, \
                         lk_83, lk_84, lk_85, lk_86, lk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * lk_80[k];

        t_9[k] = f_0 * lk_81[k];

        t_10[k] = f_0 * lk_82[k];

        t_11[k] = f_0 * lk_83[k];

        t_12[k] = f_0 * lk_84[k];

        t_13[k] = f_0 * lk_85[k];

        t_14[k] = f_0 * lk_86[k];

        t_15[k] = f_0 * lk_87[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, lk_88, lk_89, lk_90, \
                         lk_91, lk_92, lk_93, lk_94, lk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * lk_88[k];

        t_17[k] = f_0 * lk_89[k];

        t_18[k] = f_0 * lk_90[k];

        t_19[k] = f_0 * lk_91[k];

        t_20[k] = f_0 * lk_92[k];

        t_21[k] = f_0 * lk_93[k];

        t_22[k] = f_0 * lk_94[k];

        t_23[k] = f_0 * lk_95[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, lk_96, lk_97, lk_98, \
                         lk_99, lk_100, lk_101, lk_102, lk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * lk_96[k];

        t_25[k] = f_0 * lk_97[k];

        t_26[k] = f_0 * lk_98[k];

        t_27[k] = f_0 * lk_99[k];

        t_28[k] = f_0 * lk_100[k];

        t_29[k] = f_0 * lk_101[k];

        t_30[k] = f_0 * lk_102[k];

        t_31[k] = f_0 * lk_103[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, lk_104, lk_105, \
                         lk_106, lk_107, lk_144, lk_145, lk_146, \
                         lk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * lk_104[k];

        t_33[k] = f_0 * lk_105[k];

        t_34[k] = f_0 * lk_106[k];

        t_35[k] = f_0 * lk_107[k];

        t_36[k] = f_0 * lk_144[k];

        t_37[k] = f_0 * lk_145[k];

        t_38[k] = f_0 * lk_146[k];

        t_39[k] = f_0 * lk_147[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, lk_148, lk_149, \
                         lk_150, lk_151, lk_152, lk_153, lk_154, \
                         lk_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * lk_148[k];

        t_41[k] = f_0 * lk_149[k];

        t_42[k] = f_0 * lk_150[k];

        t_43[k] = f_0 * lk_151[k];

        t_44[k] = f_0 * lk_152[k];

        t_45[k] = f_0 * lk_153[k];

        t_46[k] = f_0 * lk_154[k];

        t_47[k] = f_0 * lk_155[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, lk_156, lk_157, \
                         lk_158, lk_159, lk_160, lk_161, lk_162, \
                         lk_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * lk_156[k];

        t_49[k] = f_0 * lk_157[k];

        t_50[k] = f_0 * lk_158[k];

        t_51[k] = f_0 * lk_159[k];

        t_52[k] = f_0 * lk_160[k];

        t_53[k] = f_0 * lk_161[k];

        t_54[k] = f_0 * lk_162[k];

        t_55[k] = f_0 * lk_163[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, lk_164, lk_165, \
                         lk_166, lk_167, lk_168, lk_169, lk_170, \
                         lk_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * lk_164[k];

        t_57[k] = f_0 * lk_165[k];

        t_58[k] = f_0 * lk_166[k];

        t_59[k] = f_0 * lk_167[k];

        t_60[k] = f_0 * lk_168[k];

        t_61[k] = f_0 * lk_169[k];

        t_62[k] = f_0 * lk_170[k];

        t_63[k] = f_0 * lk_171[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, lk_172, lk_173, \
                         lk_174, lk_175, lk_176, lk_177, lk_178, \
                         lk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * lk_172[k];

        t_65[k] = f_0 * lk_173[k];

        t_66[k] = f_0 * lk_174[k];

        t_67[k] = f_0 * lk_175[k];

        t_68[k] = f_0 * lk_176[k];

        t_69[k] = f_0 * lk_177[k];

        t_70[k] = f_0 * lk_178[k];

        t_71[k] = f_0 * lk_179[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ik_0, ik_1, ik_2, ik_3, ik_4, lk_180, \
                         lk_181, lk_182, lk_183, lk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -ik_0[k]
                  + f_0 * lk_180[k];

        t_73[k] = -ik_1[k]
                  + f_0 * lk_181[k];

        t_74[k] = -ik_2[k]
                  + f_0 * lk_182[k];

        t_75[k] = -ik_3[k]
                  + f_0 * lk_183[k];

        t_76[k] = -ik_4[k]
                  + f_0 * lk_184[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, ik_5, ik_6, ik_7, ik_8, ik_9, lk_185, \
                         lk_186, lk_187, lk_188, lk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -ik_5[k]
                  + f_0 * lk_185[k];

        t_78[k] = -ik_6[k]
                  + f_0 * lk_186[k];

        t_79[k] = -ik_7[k]
                  + f_0 * lk_187[k];

        t_80[k] = -ik_8[k]
                  + f_0 * lk_188[k];

        t_81[k] = -ik_9[k]
                  + f_0 * lk_189[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, ik_10, ik_11, ik_12, ik_13, ik_14, \
                         lk_190, lk_191, lk_192, lk_193, lk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -ik_10[k]
                  + f_0 * lk_190[k];

        t_83[k] = -ik_11[k]
                  + f_0 * lk_191[k];

        t_84[k] = -ik_12[k]
                  + f_0 * lk_192[k];

        t_85[k] = -ik_13[k]
                  + f_0 * lk_193[k];

        t_86[k] = -ik_14[k]
                  + f_0 * lk_194[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, ik_15, ik_16, ik_17, ik_18, ik_19, \
                         lk_195, lk_196, lk_197, lk_198, lk_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -ik_15[k]
                  + f_0 * lk_195[k];

        t_88[k] = -ik_16[k]
                  + f_0 * lk_196[k];

        t_89[k] = -ik_17[k]
                  + f_0 * lk_197[k];

        t_90[k] = -ik_18[k]
                  + f_0 * lk_198[k];

        t_91[k] = -ik_19[k]
                  + f_0 * lk_199[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, ik_20, ik_21, ik_22, ik_23, ik_24, \
                         lk_200, lk_201, lk_202, lk_203, lk_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = -ik_20[k]
                  + f_0 * lk_200[k];

        t_93[k] = -ik_21[k]
                  + f_0 * lk_201[k];

        t_94[k] = -ik_22[k]
                  + f_0 * lk_202[k];

        t_95[k] = -ik_23[k]
                  + f_0 * lk_203[k];

        t_96[k] = -ik_24[k]
                  + f_0 * lk_204[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, ik_25, ik_26, ik_27, ik_28, ik_29, \
                         lk_205, lk_206, lk_207, lk_208, lk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -ik_25[k]
                  + f_0 * lk_205[k];

        t_98[k] = -ik_26[k]
                  + f_0 * lk_206[k];

        t_99[k] = -ik_27[k]
                  + f_0 * lk_207[k];

        t_100[k] = -ik_28[k]
                   + f_0 * lk_208[k];

        t_101[k] = -ik_29[k]
                   + f_0 * lk_209[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, ik_30, ik_31, ik_32, ik_33, ik_34, \
                         lk_210, lk_211, lk_212, lk_213, lk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -ik_30[k]
                   + f_0 * lk_210[k];

        t_103[k] = -ik_31[k]
                   + f_0 * lk_211[k];

        t_104[k] = -ik_32[k]
                   + f_0 * lk_212[k];

        t_105[k] = -ik_33[k]
                   + f_0 * lk_213[k];

        t_106[k] = -ik_34[k]
                   + f_0 * lk_214[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, ik_35, lk_215, \
                         lk_252, lk_253, lk_254, lk_255, lk_256, \
                         lk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -ik_35[k]
                   + f_0 * lk_215[k];

        t_108[k] = f_0 * lk_252[k];

        t_109[k] = f_0 * lk_253[k];

        t_110[k] = f_0 * lk_254[k];

        t_111[k] = f_0 * lk_255[k];

        t_112[k] = f_0 * lk_256[k];

        t_113[k] = f_0 * lk_257[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, t_120, t_121, lk_258, \
                         lk_259, lk_260, lk_261, lk_262, lk_263, lk_264, \
                         lk_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * lk_258[k];

        t_115[k] = f_0 * lk_259[k];

        t_116[k] = f_0 * lk_260[k];

        t_117[k] = f_0 * lk_261[k];

        t_118[k] = f_0 * lk_262[k];

        t_119[k] = f_0 * lk_263[k];

        t_120[k] = f_0 * lk_264[k];

        t_121[k] = f_0 * lk_265[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, t_128, t_129, lk_266, \
                         lk_267, lk_268, lk_269, lk_270, lk_271, lk_272, \
                         lk_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * lk_266[k];

        t_123[k] = f_0 * lk_267[k];

        t_124[k] = f_0 * lk_268[k];

        t_125[k] = f_0 * lk_269[k];

        t_126[k] = f_0 * lk_270[k];

        t_127[k] = f_0 * lk_271[k];

        t_128[k] = f_0 * lk_272[k];

        t_129[k] = f_0 * lk_273[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, t_135, t_136, t_137, lk_274, \
                         lk_275, lk_276, lk_277, lk_278, lk_279, lk_280, \
                         lk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * lk_274[k];

        t_131[k] = f_0 * lk_275[k];

        t_132[k] = f_0 * lk_276[k];

        t_133[k] = f_0 * lk_277[k];

        t_134[k] = f_0 * lk_278[k];

        t_135[k] = f_0 * lk_279[k];

        t_136[k] = f_0 * lk_280[k];

        t_137[k] = f_0 * lk_281[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, t_144, ik_36, lk_282, \
                         lk_283, lk_284, lk_285, lk_286, lk_287, \
                         lk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * lk_282[k];

        t_139[k] = f_0 * lk_283[k];

        t_140[k] = f_0 * lk_284[k];

        t_141[k] = f_0 * lk_285[k];

        t_142[k] = f_0 * lk_286[k];

        t_143[k] = f_0 * lk_287[k];

        t_144[k] = -ik_36[k]
                   + f_0 * lk_288[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ik_37, ik_38, ik_39, ik_40, ik_41, \
                         lk_289, lk_290, lk_291, lk_292, lk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -ik_37[k]
                   + f_0 * lk_289[k];

        t_146[k] = -ik_38[k]
                   + f_0 * lk_290[k];

        t_147[k] = -ik_39[k]
                   + f_0 * lk_291[k];

        t_148[k] = -ik_40[k]
                   + f_0 * lk_292[k];

        t_149[k] = -ik_41[k]
                   + f_0 * lk_293[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ik_42, ik_43, ik_44, ik_45, ik_46, \
                         lk_294, lk_295, lk_296, lk_297, lk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -ik_42[k]
                   + f_0 * lk_294[k];

        t_151[k] = -ik_43[k]
                   + f_0 * lk_295[k];

        t_152[k] = -ik_44[k]
                   + f_0 * lk_296[k];

        t_153[k] = -ik_45[k]
                   + f_0 * lk_297[k];

        t_154[k] = -ik_46[k]
                   + f_0 * lk_298[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ik_47, ik_48, ik_49, ik_50, ik_51, \
                         lk_299, lk_300, lk_301, lk_302, lk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -ik_47[k]
                   + f_0 * lk_299[k];

        t_156[k] = -ik_48[k]
                   + f_0 * lk_300[k];

        t_157[k] = -ik_49[k]
                   + f_0 * lk_301[k];

        t_158[k] = -ik_50[k]
                   + f_0 * lk_302[k];

        t_159[k] = -ik_51[k]
                   + f_0 * lk_303[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ik_52, ik_53, ik_54, ik_55, ik_56, \
                         lk_304, lk_305, lk_306, lk_307, lk_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -ik_52[k]
                   + f_0 * lk_304[k];

        t_161[k] = -ik_53[k]
                   + f_0 * lk_305[k];

        t_162[k] = -ik_54[k]
                   + f_0 * lk_306[k];

        t_163[k] = -ik_55[k]
                   + f_0 * lk_307[k];

        t_164[k] = -ik_56[k]
                   + f_0 * lk_308[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ik_57, ik_58, ik_59, ik_60, ik_61, \
                         lk_309, lk_310, lk_311, lk_312, lk_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -ik_57[k]
                   + f_0 * lk_309[k];

        t_166[k] = -ik_58[k]
                   + f_0 * lk_310[k];

        t_167[k] = -ik_59[k]
                   + f_0 * lk_311[k];

        t_168[k] = -ik_60[k]
                   + f_0 * lk_312[k];

        t_169[k] = -ik_61[k]
                   + f_0 * lk_313[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ik_62, ik_63, ik_64, ik_65, ik_66, \
                         lk_314, lk_315, lk_316, lk_317, lk_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -ik_62[k]
                   + f_0 * lk_314[k];

        t_171[k] = -ik_63[k]
                   + f_0 * lk_315[k];

        t_172[k] = -ik_64[k]
                   + f_0 * lk_316[k];

        t_173[k] = -ik_65[k]
                   + f_0 * lk_317[k];

        t_174[k] = -ik_66[k]
                   + f_0 * lk_318[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ik_67, ik_68, ik_69, ik_70, ik_71, \
                         lk_319, lk_320, lk_321, lk_322, lk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -ik_67[k]
                   + f_0 * lk_319[k];

        t_176[k] = -ik_68[k]
                   + f_0 * lk_320[k];

        t_177[k] = -ik_69[k]
                   + f_0 * lk_321[k];

        t_178[k] = -ik_70[k]
                   + f_0 * lk_322[k];

        t_179[k] = -ik_71[k]
                   + f_0 * lk_323[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ik_72, ik_73, ik_74, ik_75, ik_76, \
                         lk_324, lk_325, lk_326, lk_327, lk_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * ik_72[k]
                   + f_0 * lk_324[k];

        t_181[k] = -2.0 * ik_73[k]
                   + f_0 * lk_325[k];

        t_182[k] = -2.0 * ik_74[k]
                   + f_0 * lk_326[k];

        t_183[k] = -2.0 * ik_75[k]
                   + f_0 * lk_327[k];

        t_184[k] = -2.0 * ik_76[k]
                   + f_0 * lk_328[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_77 = buffer.data(ik + 77);
    const auto *ik_78 = buffer.data(ik + 78);
    const auto *ik_79 = buffer.data(ik + 79);
    const auto *ik_80 = buffer.data(ik + 80);
    const auto *ik_81 = buffer.data(ik + 81);
    const auto *ik_82 = buffer.data(ik + 82);
    const auto *ik_83 = buffer.data(ik + 83);
    const auto *ik_84 = buffer.data(ik + 84);
    const auto *ik_85 = buffer.data(ik + 85);
    const auto *ik_86 = buffer.data(ik + 86);
    const auto *ik_87 = buffer.data(ik + 87);
    const auto *ik_88 = buffer.data(ik + 88);
    const auto *ik_89 = buffer.data(ik + 89);
    const auto *ik_90 = buffer.data(ik + 90);
    const auto *ik_91 = buffer.data(ik + 91);
    const auto *ik_92 = buffer.data(ik + 92);
    const auto *ik_93 = buffer.data(ik + 93);
    const auto *ik_94 = buffer.data(ik + 94);
    const auto *ik_95 = buffer.data(ik + 95);
    const auto *ik_96 = buffer.data(ik + 96);
    const auto *ik_97 = buffer.data(ik + 97);
    const auto *ik_98 = buffer.data(ik + 98);
    const auto *ik_99 = buffer.data(ik + 99);
    const auto *ik_100 = buffer.data(ik + 100);
    const auto *ik_101 = buffer.data(ik + 101);
    const auto *ik_102 = buffer.data(ik + 102);
    const auto *ik_103 = buffer.data(ik + 103);
    const auto *ik_104 = buffer.data(ik + 104);
    const auto *ik_105 = buffer.data(ik + 105);
    const auto *ik_106 = buffer.data(ik + 106);
    const auto *ik_107 = buffer.data(ik + 107);
    const auto *ik_108 = buffer.data(ik + 108);
    const auto *ik_109 = buffer.data(ik + 109);
    const auto *ik_110 = buffer.data(ik + 110);
    const auto *ik_111 = buffer.data(ik + 111);
    const auto *ik_112 = buffer.data(ik + 112);
    const auto *ik_113 = buffer.data(ik + 113);
    const auto *ik_114 = buffer.data(ik + 114);
    const auto *ik_115 = buffer.data(ik + 115);
    const auto *ik_116 = buffer.data(ik + 116);
    const auto *ik_117 = buffer.data(ik + 117);
    const auto *ik_118 = buffer.data(ik + 118);
    const auto *ik_119 = buffer.data(ik + 119);
    const auto *ik_120 = buffer.data(ik + 120);
    const auto *ik_121 = buffer.data(ik + 121);
    const auto *ik_122 = buffer.data(ik + 122);
    const auto *ik_123 = buffer.data(ik + 123);
    const auto *ik_124 = buffer.data(ik + 124);
    const auto *ik_125 = buffer.data(ik + 125);
    const auto *ik_126 = buffer.data(ik + 126);
    const auto *ik_127 = buffer.data(ik + 127);
    const auto *ik_128 = buffer.data(ik + 128);
    const auto *ik_129 = buffer.data(ik + 129);
    const auto *ik_130 = buffer.data(ik + 130);
    const auto *ik_131 = buffer.data(ik + 131);
    const auto *ik_132 = buffer.data(ik + 132);
    const auto *ik_133 = buffer.data(ik + 133);
    const auto *ik_134 = buffer.data(ik + 134);
    const auto *ik_135 = buffer.data(ik + 135);
    const auto *ik_136 = buffer.data(ik + 136);
    const auto *ik_137 = buffer.data(ik + 137);
    const auto *ik_138 = buffer.data(ik + 138);
    const auto *ik_139 = buffer.data(ik + 139);
    const auto *ik_140 = buffer.data(ik + 140);
    const auto *ik_141 = buffer.data(ik + 141);
    const auto *ik_142 = buffer.data(ik + 142);
    const auto *ik_143 = buffer.data(ik + 143);
    const auto *ik_144 = buffer.data(ik + 144);
    const auto *ik_145 = buffer.data(ik + 145);
    const auto *ik_146 = buffer.data(ik + 146);
    const auto *ik_147 = buffer.data(ik + 147);
    const auto *ik_148 = buffer.data(ik + 148);
    const auto *ik_149 = buffer.data(ik + 149);
    const auto *ik_150 = buffer.data(ik + 150);
    const auto *ik_151 = buffer.data(ik + 151);
    const auto *ik_152 = buffer.data(ik + 152);
    const auto *ik_153 = buffer.data(ik + 153);
    const auto *ik_154 = buffer.data(ik + 154);
    const auto *ik_155 = buffer.data(ik + 155);
    const auto *ik_156 = buffer.data(ik + 156);
    const auto *ik_157 = buffer.data(ik + 157);
    const auto *ik_158 = buffer.data(ik + 158);
    const auto *ik_159 = buffer.data(ik + 159);
    const auto *ik_160 = buffer.data(ik + 160);
    const auto *ik_161 = buffer.data(ik + 161);
    const auto *ik_162 = buffer.data(ik + 162);
    const auto *ik_163 = buffer.data(ik + 163);
    const auto *ik_164 = buffer.data(ik + 164);
    const auto *ik_165 = buffer.data(ik + 165);
    const auto *ik_166 = buffer.data(ik + 166);
    const auto *ik_167 = buffer.data(ik + 167);
    const auto *ik_168 = buffer.data(ik + 168);
    const auto *ik_169 = buffer.data(ik + 169);
    const auto *ik_170 = buffer.data(ik + 170);
    const auto *ik_171 = buffer.data(ik + 171);
    const auto *ik_172 = buffer.data(ik + 172);
    const auto *ik_173 = buffer.data(ik + 173);
    const auto *ik_174 = buffer.data(ik + 174);
    const auto *ik_175 = buffer.data(ik + 175);
    const auto *ik_176 = buffer.data(ik + 176);
    const auto *ik_177 = buffer.data(ik + 177);
    const auto *ik_178 = buffer.data(ik + 178);
    const auto *ik_179 = buffer.data(ik + 179);
    const auto *ik_180 = buffer.data(ik + 180);
    const auto *ik_181 = buffer.data(ik + 181);
    const auto *ik_182 = buffer.data(ik + 182);
    const auto *ik_183 = buffer.data(ik + 183);
    const auto *ik_184 = buffer.data(ik + 184);
    const auto *ik_185 = buffer.data(ik + 185);
    const auto *ik_186 = buffer.data(ik + 186);
    const auto *ik_187 = buffer.data(ik + 187);
    const auto *ik_188 = buffer.data(ik + 188);
    const auto *ik_189 = buffer.data(ik + 189);
    const auto *ik_190 = buffer.data(ik + 190);
    const auto *ik_191 = buffer.data(ik + 191);
    const auto *ik_192 = buffer.data(ik + 192);
    const auto *ik_193 = buffer.data(ik + 193);
    const auto *ik_194 = buffer.data(ik + 194);
    const auto *ik_195 = buffer.data(ik + 195);
    const auto *ik_196 = buffer.data(ik + 196);
    const auto *ik_197 = buffer.data(ik + 197);
    const auto *ik_198 = buffer.data(ik + 198);
    const auto *ik_199 = buffer.data(ik + 199);
    const auto *ik_200 = buffer.data(ik + 200);
    const auto *ik_201 = buffer.data(ik + 201);
    const auto *ik_202 = buffer.data(ik + 202);
    const auto *ik_203 = buffer.data(ik + 203);

    const auto *lk_329 = buffer.data(lk + 329);
    const auto *lk_330 = buffer.data(lk + 330);
    const auto *lk_331 = buffer.data(lk + 331);
    const auto *lk_332 = buffer.data(lk + 332);
    const auto *lk_333 = buffer.data(lk + 333);
    const auto *lk_334 = buffer.data(lk + 334);
    const auto *lk_335 = buffer.data(lk + 335);
    const auto *lk_336 = buffer.data(lk + 336);
    const auto *lk_337 = buffer.data(lk + 337);
    const auto *lk_338 = buffer.data(lk + 338);
    const auto *lk_339 = buffer.data(lk + 339);
    const auto *lk_340 = buffer.data(lk + 340);
    const auto *lk_341 = buffer.data(lk + 341);
    const auto *lk_342 = buffer.data(lk + 342);
    const auto *lk_343 = buffer.data(lk + 343);
    const auto *lk_344 = buffer.data(lk + 344);
    const auto *lk_345 = buffer.data(lk + 345);
    const auto *lk_346 = buffer.data(lk + 346);
    const auto *lk_347 = buffer.data(lk + 347);
    const auto *lk_348 = buffer.data(lk + 348);
    const auto *lk_349 = buffer.data(lk + 349);
    const auto *lk_350 = buffer.data(lk + 350);
    const auto *lk_351 = buffer.data(lk + 351);
    const auto *lk_352 = buffer.data(lk + 352);
    const auto *lk_353 = buffer.data(lk + 353);
    const auto *lk_354 = buffer.data(lk + 354);
    const auto *lk_355 = buffer.data(lk + 355);
    const auto *lk_356 = buffer.data(lk + 356);
    const auto *lk_357 = buffer.data(lk + 357);
    const auto *lk_358 = buffer.data(lk + 358);
    const auto *lk_359 = buffer.data(lk + 359);
    const auto *lk_396 = buffer.data(lk + 396);
    const auto *lk_397 = buffer.data(lk + 397);
    const auto *lk_398 = buffer.data(lk + 398);
    const auto *lk_399 = buffer.data(lk + 399);
    const auto *lk_400 = buffer.data(lk + 400);
    const auto *lk_401 = buffer.data(lk + 401);
    const auto *lk_402 = buffer.data(lk + 402);
    const auto *lk_403 = buffer.data(lk + 403);
    const auto *lk_404 = buffer.data(lk + 404);
    const auto *lk_405 = buffer.data(lk + 405);
    const auto *lk_406 = buffer.data(lk + 406);
    const auto *lk_407 = buffer.data(lk + 407);
    const auto *lk_408 = buffer.data(lk + 408);
    const auto *lk_409 = buffer.data(lk + 409);
    const auto *lk_410 = buffer.data(lk + 410);
    const auto *lk_411 = buffer.data(lk + 411);
    const auto *lk_412 = buffer.data(lk + 412);
    const auto *lk_413 = buffer.data(lk + 413);
    const auto *lk_414 = buffer.data(lk + 414);
    const auto *lk_415 = buffer.data(lk + 415);
    const auto *lk_416 = buffer.data(lk + 416);
    const auto *lk_417 = buffer.data(lk + 417);
    const auto *lk_418 = buffer.data(lk + 418);
    const auto *lk_419 = buffer.data(lk + 419);
    const auto *lk_420 = buffer.data(lk + 420);
    const auto *lk_421 = buffer.data(lk + 421);
    const auto *lk_422 = buffer.data(lk + 422);
    const auto *lk_423 = buffer.data(lk + 423);
    const auto *lk_424 = buffer.data(lk + 424);
    const auto *lk_425 = buffer.data(lk + 425);
    const auto *lk_426 = buffer.data(lk + 426);
    const auto *lk_427 = buffer.data(lk + 427);
    const auto *lk_428 = buffer.data(lk + 428);
    const auto *lk_429 = buffer.data(lk + 429);
    const auto *lk_430 = buffer.data(lk + 430);
    const auto *lk_431 = buffer.data(lk + 431);
    const auto *lk_432 = buffer.data(lk + 432);
    const auto *lk_433 = buffer.data(lk + 433);
    const auto *lk_434 = buffer.data(lk + 434);
    const auto *lk_435 = buffer.data(lk + 435);
    const auto *lk_436 = buffer.data(lk + 436);
    const auto *lk_437 = buffer.data(lk + 437);
    const auto *lk_438 = buffer.data(lk + 438);
    const auto *lk_439 = buffer.data(lk + 439);
    const auto *lk_440 = buffer.data(lk + 440);
    const auto *lk_441 = buffer.data(lk + 441);
    const auto *lk_442 = buffer.data(lk + 442);
    const auto *lk_443 = buffer.data(lk + 443);
    const auto *lk_444 = buffer.data(lk + 444);
    const auto *lk_445 = buffer.data(lk + 445);
    const auto *lk_446 = buffer.data(lk + 446);
    const auto *lk_447 = buffer.data(lk + 447);
    const auto *lk_448 = buffer.data(lk + 448);
    const auto *lk_449 = buffer.data(lk + 449);
    const auto *lk_450 = buffer.data(lk + 450);
    const auto *lk_451 = buffer.data(lk + 451);
    const auto *lk_452 = buffer.data(lk + 452);
    const auto *lk_453 = buffer.data(lk + 453);
    const auto *lk_454 = buffer.data(lk + 454);
    const auto *lk_455 = buffer.data(lk + 455);
    const auto *lk_456 = buffer.data(lk + 456);
    const auto *lk_457 = buffer.data(lk + 457);
    const auto *lk_458 = buffer.data(lk + 458);
    const auto *lk_459 = buffer.data(lk + 459);
    const auto *lk_460 = buffer.data(lk + 460);
    const auto *lk_461 = buffer.data(lk + 461);
    const auto *lk_462 = buffer.data(lk + 462);
    const auto *lk_463 = buffer.data(lk + 463);
    const auto *lk_464 = buffer.data(lk + 464);
    const auto *lk_465 = buffer.data(lk + 465);
    const auto *lk_466 = buffer.data(lk + 466);
    const auto *lk_467 = buffer.data(lk + 467);
    const auto *lk_468 = buffer.data(lk + 468);
    const auto *lk_469 = buffer.data(lk + 469);
    const auto *lk_470 = buffer.data(lk + 470);
    const auto *lk_471 = buffer.data(lk + 471);
    const auto *lk_472 = buffer.data(lk + 472);
    const auto *lk_473 = buffer.data(lk + 473);
    const auto *lk_474 = buffer.data(lk + 474);
    const auto *lk_475 = buffer.data(lk + 475);
    const auto *lk_476 = buffer.data(lk + 476);
    const auto *lk_477 = buffer.data(lk + 477);
    const auto *lk_478 = buffer.data(lk + 478);
    const auto *lk_479 = buffer.data(lk + 479);
    const auto *lk_480 = buffer.data(lk + 480);
    const auto *lk_481 = buffer.data(lk + 481);
    const auto *lk_482 = buffer.data(lk + 482);
    const auto *lk_483 = buffer.data(lk + 483);
    const auto *lk_484 = buffer.data(lk + 484);
    const auto *lk_485 = buffer.data(lk + 485);
    const auto *lk_486 = buffer.data(lk + 486);
    const auto *lk_487 = buffer.data(lk + 487);
    const auto *lk_488 = buffer.data(lk + 488);
    const auto *lk_489 = buffer.data(lk + 489);
    const auto *lk_490 = buffer.data(lk + 490);
    const auto *lk_491 = buffer.data(lk + 491);
    const auto *lk_492 = buffer.data(lk + 492);
    const auto *lk_493 = buffer.data(lk + 493);
    const auto *lk_494 = buffer.data(lk + 494);
    const auto *lk_495 = buffer.data(lk + 495);
    const auto *lk_496 = buffer.data(lk + 496);
    const auto *lk_497 = buffer.data(lk + 497);
    const auto *lk_498 = buffer.data(lk + 498);
    const auto *lk_499 = buffer.data(lk + 499);
    const auto *lk_500 = buffer.data(lk + 500);
    const auto *lk_501 = buffer.data(lk + 501);
    const auto *lk_502 = buffer.data(lk + 502);
    const auto *lk_503 = buffer.data(lk + 503);
    const auto *lk_504 = buffer.data(lk + 504);
    const auto *lk_505 = buffer.data(lk + 505);
    const auto *lk_506 = buffer.data(lk + 506);
    const auto *lk_507 = buffer.data(lk + 507);
    const auto *lk_508 = buffer.data(lk + 508);
    const auto *lk_509 = buffer.data(lk + 509);
    const auto *lk_510 = buffer.data(lk + 510);
    const auto *lk_511 = buffer.data(lk + 511);
    const auto *lk_512 = buffer.data(lk + 512);
    const auto *lk_513 = buffer.data(lk + 513);
    const auto *lk_514 = buffer.data(lk + 514);
    const auto *lk_515 = buffer.data(lk + 515);
    const auto *lk_516 = buffer.data(lk + 516);
    const auto *lk_517 = buffer.data(lk + 517);
    const auto *lk_518 = buffer.data(lk + 518);
    const auto *lk_519 = buffer.data(lk + 519);
    const auto *lk_520 = buffer.data(lk + 520);
    const auto *lk_521 = buffer.data(lk + 521);
    const auto *lk_522 = buffer.data(lk + 522);
    const auto *lk_523 = buffer.data(lk + 523);
    const auto *lk_524 = buffer.data(lk + 524);
    const auto *lk_525 = buffer.data(lk + 525);
    const auto *lk_526 = buffer.data(lk + 526);
    const auto *lk_527 = buffer.data(lk + 527);

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ik_77, ik_78, ik_79, ik_80, ik_81, \
                         lk_329, lk_330, lk_331, lk_332, lk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * ik_77[k]
                   + f_0 * lk_329[k];

        t_186[k] = -2.0 * ik_78[k]
                   + f_0 * lk_330[k];

        t_187[k] = -2.0 * ik_79[k]
                   + f_0 * lk_331[k];

        t_188[k] = -2.0 * ik_80[k]
                   + f_0 * lk_332[k];

        t_189[k] = -2.0 * ik_81[k]
                   + f_0 * lk_333[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ik_82, ik_83, ik_84, ik_85, ik_86, \
                         lk_334, lk_335, lk_336, lk_337, lk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * ik_82[k]
                   + f_0 * lk_334[k];

        t_191[k] = -2.0 * ik_83[k]
                   + f_0 * lk_335[k];

        t_192[k] = -2.0 * ik_84[k]
                   + f_0 * lk_336[k];

        t_193[k] = -2.0 * ik_85[k]
                   + f_0 * lk_337[k];

        t_194[k] = -2.0 * ik_86[k]
                   + f_0 * lk_338[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ik_87, ik_88, ik_89, ik_90, ik_91, \
                         lk_339, lk_340, lk_341, lk_342, lk_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * ik_87[k]
                   + f_0 * lk_339[k];

        t_196[k] = -2.0 * ik_88[k]
                   + f_0 * lk_340[k];

        t_197[k] = -2.0 * ik_89[k]
                   + f_0 * lk_341[k];

        t_198[k] = -2.0 * ik_90[k]
                   + f_0 * lk_342[k];

        t_199[k] = -2.0 * ik_91[k]
                   + f_0 * lk_343[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ik_92, ik_93, ik_94, ik_95, ik_96, \
                         lk_344, lk_345, lk_346, lk_347, lk_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * ik_92[k]
                   + f_0 * lk_344[k];

        t_201[k] = -2.0 * ik_93[k]
                   + f_0 * lk_345[k];

        t_202[k] = -2.0 * ik_94[k]
                   + f_0 * lk_346[k];

        t_203[k] = -2.0 * ik_95[k]
                   + f_0 * lk_347[k];

        t_204[k] = -2.0 * ik_96[k]
                   + f_0 * lk_348[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ik_97, ik_98, ik_99, ik_100, \
                         ik_101, lk_349, lk_350, lk_351, lk_352, \
                         lk_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * ik_97[k]
                   + f_0 * lk_349[k];

        t_206[k] = -2.0 * ik_98[k]
                   + f_0 * lk_350[k];

        t_207[k] = -2.0 * ik_99[k]
                   + f_0 * lk_351[k];

        t_208[k] = -2.0 * ik_100[k]
                   + f_0 * lk_352[k];

        t_209[k] = -2.0 * ik_101[k]
                   + f_0 * lk_353[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ik_102, ik_103, ik_104, ik_105, \
                         ik_106, lk_354, lk_355, lk_356, lk_357, \
                         lk_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * ik_102[k]
                   + f_0 * lk_354[k];

        t_211[k] = -2.0 * ik_103[k]
                   + f_0 * lk_355[k];

        t_212[k] = -2.0 * ik_104[k]
                   + f_0 * lk_356[k];

        t_213[k] = -2.0 * ik_105[k]
                   + f_0 * lk_357[k];

        t_214[k] = -2.0 * ik_106[k]
                   + f_0 * lk_358[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, t_220, t_221, ik_107, lk_359, \
                         lk_396, lk_397, lk_398, lk_399, lk_400, \
                         lk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * ik_107[k]
                   + f_0 * lk_359[k];

        t_216[k] = f_0 * lk_396[k];

        t_217[k] = f_0 * lk_397[k];

        t_218[k] = f_0 * lk_398[k];

        t_219[k] = f_0 * lk_399[k];

        t_220[k] = f_0 * lk_400[k];

        t_221[k] = f_0 * lk_401[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, t_228, t_229, lk_402, \
                         lk_403, lk_404, lk_405, lk_406, lk_407, lk_408, \
                         lk_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_0 * lk_402[k];

        t_223[k] = f_0 * lk_403[k];

        t_224[k] = f_0 * lk_404[k];

        t_225[k] = f_0 * lk_405[k];

        t_226[k] = f_0 * lk_406[k];

        t_227[k] = f_0 * lk_407[k];

        t_228[k] = f_0 * lk_408[k];

        t_229[k] = f_0 * lk_409[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, t_235, t_236, t_237, lk_410, \
                         lk_411, lk_412, lk_413, lk_414, lk_415, lk_416, \
                         lk_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_0 * lk_410[k];

        t_231[k] = f_0 * lk_411[k];

        t_232[k] = f_0 * lk_412[k];

        t_233[k] = f_0 * lk_413[k];

        t_234[k] = f_0 * lk_414[k];

        t_235[k] = f_0 * lk_415[k];

        t_236[k] = f_0 * lk_416[k];

        t_237[k] = f_0 * lk_417[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, t_243, t_244, t_245, lk_418, \
                         lk_419, lk_420, lk_421, lk_422, lk_423, lk_424, \
                         lk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_0 * lk_418[k];

        t_239[k] = f_0 * lk_419[k];

        t_240[k] = f_0 * lk_420[k];

        t_241[k] = f_0 * lk_421[k];

        t_242[k] = f_0 * lk_422[k];

        t_243[k] = f_0 * lk_423[k];

        t_244[k] = f_0 * lk_424[k];

        t_245[k] = f_0 * lk_425[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, t_251, t_252, ik_108, lk_426, \
                         lk_427, lk_428, lk_429, lk_430, lk_431, \
                         lk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_0 * lk_426[k];

        t_247[k] = f_0 * lk_427[k];

        t_248[k] = f_0 * lk_428[k];

        t_249[k] = f_0 * lk_429[k];

        t_250[k] = f_0 * lk_430[k];

        t_251[k] = f_0 * lk_431[k];

        t_252[k] = -ik_108[k]
                   + f_0 * lk_432[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, ik_109, ik_110, ik_111, ik_112, \
                         ik_113, lk_433, lk_434, lk_435, lk_436, \
                         lk_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = -ik_109[k]
                   + f_0 * lk_433[k];

        t_254[k] = -ik_110[k]
                   + f_0 * lk_434[k];

        t_255[k] = -ik_111[k]
                   + f_0 * lk_435[k];

        t_256[k] = -ik_112[k]
                   + f_0 * lk_436[k];

        t_257[k] = -ik_113[k]
                   + f_0 * lk_437[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, ik_114, ik_115, ik_116, ik_117, \
                         ik_118, lk_438, lk_439, lk_440, lk_441, \
                         lk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = -ik_114[k]
                   + f_0 * lk_438[k];

        t_259[k] = -ik_115[k]
                   + f_0 * lk_439[k];

        t_260[k] = -ik_116[k]
                   + f_0 * lk_440[k];

        t_261[k] = -ik_117[k]
                   + f_0 * lk_441[k];

        t_262[k] = -ik_118[k]
                   + f_0 * lk_442[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, ik_119, ik_120, ik_121, ik_122, \
                         ik_123, lk_443, lk_444, lk_445, lk_446, \
                         lk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = -ik_119[k]
                   + f_0 * lk_443[k];

        t_264[k] = -ik_120[k]
                   + f_0 * lk_444[k];

        t_265[k] = -ik_121[k]
                   + f_0 * lk_445[k];

        t_266[k] = -ik_122[k]
                   + f_0 * lk_446[k];

        t_267[k] = -ik_123[k]
                   + f_0 * lk_447[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, ik_124, ik_125, ik_126, ik_127, \
                         ik_128, lk_448, lk_449, lk_450, lk_451, \
                         lk_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = -ik_124[k]
                   + f_0 * lk_448[k];

        t_269[k] = -ik_125[k]
                   + f_0 * lk_449[k];

        t_270[k] = -ik_126[k]
                   + f_0 * lk_450[k];

        t_271[k] = -ik_127[k]
                   + f_0 * lk_451[k];

        t_272[k] = -ik_128[k]
                   + f_0 * lk_452[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, ik_129, ik_130, ik_131, ik_132, \
                         ik_133, lk_453, lk_454, lk_455, lk_456, \
                         lk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = -ik_129[k]
                   + f_0 * lk_453[k];

        t_274[k] = -ik_130[k]
                   + f_0 * lk_454[k];

        t_275[k] = -ik_131[k]
                   + f_0 * lk_455[k];

        t_276[k] = -ik_132[k]
                   + f_0 * lk_456[k];

        t_277[k] = -ik_133[k]
                   + f_0 * lk_457[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, ik_134, ik_135, ik_136, ik_137, \
                         ik_138, lk_458, lk_459, lk_460, lk_461, \
                         lk_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = -ik_134[k]
                   + f_0 * lk_458[k];

        t_279[k] = -ik_135[k]
                   + f_0 * lk_459[k];

        t_280[k] = -ik_136[k]
                   + f_0 * lk_460[k];

        t_281[k] = -ik_137[k]
                   + f_0 * lk_461[k];

        t_282[k] = -ik_138[k]
                   + f_0 * lk_462[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, ik_139, ik_140, ik_141, ik_142, \
                         ik_143, lk_463, lk_464, lk_465, lk_466, \
                         lk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = -ik_139[k]
                   + f_0 * lk_463[k];

        t_284[k] = -ik_140[k]
                   + f_0 * lk_464[k];

        t_285[k] = -ik_141[k]
                   + f_0 * lk_465[k];

        t_286[k] = -ik_142[k]
                   + f_0 * lk_466[k];

        t_287[k] = -ik_143[k]
                   + f_0 * lk_467[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, ik_144, ik_145, ik_146, ik_147, \
                         ik_148, lk_468, lk_469, lk_470, lk_471, \
                         lk_472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = -2.0 * ik_144[k]
                   + f_0 * lk_468[k];

        t_289[k] = -2.0 * ik_145[k]
                   + f_0 * lk_469[k];

        t_290[k] = -2.0 * ik_146[k]
                   + f_0 * lk_470[k];

        t_291[k] = -2.0 * ik_147[k]
                   + f_0 * lk_471[k];

        t_292[k] = -2.0 * ik_148[k]
                   + f_0 * lk_472[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, t_297, ik_149, ik_150, ik_151, ik_152, \
                         ik_153, lk_473, lk_474, lk_475, lk_476, \
                         lk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = -2.0 * ik_149[k]
                   + f_0 * lk_473[k];

        t_294[k] = -2.0 * ik_150[k]
                   + f_0 * lk_474[k];

        t_295[k] = -2.0 * ik_151[k]
                   + f_0 * lk_475[k];

        t_296[k] = -2.0 * ik_152[k]
                   + f_0 * lk_476[k];

        t_297[k] = -2.0 * ik_153[k]
                   + f_0 * lk_477[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, ik_154, ik_155, ik_156, ik_157, \
                         ik_158, lk_478, lk_479, lk_480, lk_481, \
                         lk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = -2.0 * ik_154[k]
                   + f_0 * lk_478[k];

        t_299[k] = -2.0 * ik_155[k]
                   + f_0 * lk_479[k];

        t_300[k] = -2.0 * ik_156[k]
                   + f_0 * lk_480[k];

        t_301[k] = -2.0 * ik_157[k]
                   + f_0 * lk_481[k];

        t_302[k] = -2.0 * ik_158[k]
                   + f_0 * lk_482[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, ik_159, ik_160, ik_161, ik_162, \
                         ik_163, lk_483, lk_484, lk_485, lk_486, \
                         lk_487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = -2.0 * ik_159[k]
                   + f_0 * lk_483[k];

        t_304[k] = -2.0 * ik_160[k]
                   + f_0 * lk_484[k];

        t_305[k] = -2.0 * ik_161[k]
                   + f_0 * lk_485[k];

        t_306[k] = -2.0 * ik_162[k]
                   + f_0 * lk_486[k];

        t_307[k] = -2.0 * ik_163[k]
                   + f_0 * lk_487[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, ik_164, ik_165, ik_166, ik_167, \
                         ik_168, lk_488, lk_489, lk_490, lk_491, \
                         lk_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = -2.0 * ik_164[k]
                   + f_0 * lk_488[k];

        t_309[k] = -2.0 * ik_165[k]
                   + f_0 * lk_489[k];

        t_310[k] = -2.0 * ik_166[k]
                   + f_0 * lk_490[k];

        t_311[k] = -2.0 * ik_167[k]
                   + f_0 * lk_491[k];

        t_312[k] = -2.0 * ik_168[k]
                   + f_0 * lk_492[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, ik_169, ik_170, ik_171, ik_172, \
                         ik_173, lk_493, lk_494, lk_495, lk_496, \
                         lk_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = -2.0 * ik_169[k]
                   + f_0 * lk_493[k];

        t_314[k] = -2.0 * ik_170[k]
                   + f_0 * lk_494[k];

        t_315[k] = -2.0 * ik_171[k]
                   + f_0 * lk_495[k];

        t_316[k] = -2.0 * ik_172[k]
                   + f_0 * lk_496[k];

        t_317[k] = -2.0 * ik_173[k]
                   + f_0 * lk_497[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, ik_174, ik_175, ik_176, ik_177, \
                         ik_178, lk_498, lk_499, lk_500, lk_501, \
                         lk_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = -2.0 * ik_174[k]
                   + f_0 * lk_498[k];

        t_319[k] = -2.0 * ik_175[k]
                   + f_0 * lk_499[k];

        t_320[k] = -2.0 * ik_176[k]
                   + f_0 * lk_500[k];

        t_321[k] = -2.0 * ik_177[k]
                   + f_0 * lk_501[k];

        t_322[k] = -2.0 * ik_178[k]
                   + f_0 * lk_502[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, ik_179, ik_180, ik_181, ik_182, \
                         ik_183, lk_503, lk_504, lk_505, lk_506, \
                         lk_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = -2.0 * ik_179[k]
                   + f_0 * lk_503[k];

        t_324[k] = -3.0 * ik_180[k]
                   + f_0 * lk_504[k];

        t_325[k] = -3.0 * ik_181[k]
                   + f_0 * lk_505[k];

        t_326[k] = -3.0 * ik_182[k]
                   + f_0 * lk_506[k];

        t_327[k] = -3.0 * ik_183[k]
                   + f_0 * lk_507[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, ik_184, ik_185, ik_186, ik_187, \
                         ik_188, lk_508, lk_509, lk_510, lk_511, \
                         lk_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = -3.0 * ik_184[k]
                   + f_0 * lk_508[k];

        t_329[k] = -3.0 * ik_185[k]
                   + f_0 * lk_509[k];

        t_330[k] = -3.0 * ik_186[k]
                   + f_0 * lk_510[k];

        t_331[k] = -3.0 * ik_187[k]
                   + f_0 * lk_511[k];

        t_332[k] = -3.0 * ik_188[k]
                   + f_0 * lk_512[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, t_337, ik_189, ik_190, ik_191, ik_192, \
                         ik_193, lk_513, lk_514, lk_515, lk_516, \
                         lk_517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = -3.0 * ik_189[k]
                   + f_0 * lk_513[k];

        t_334[k] = -3.0 * ik_190[k]
                   + f_0 * lk_514[k];

        t_335[k] = -3.0 * ik_191[k]
                   + f_0 * lk_515[k];

        t_336[k] = -3.0 * ik_192[k]
                   + f_0 * lk_516[k];

        t_337[k] = -3.0 * ik_193[k]
                   + f_0 * lk_517[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, ik_194, ik_195, ik_196, ik_197, \
                         ik_198, lk_518, lk_519, lk_520, lk_521, \
                         lk_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = -3.0 * ik_194[k]
                   + f_0 * lk_518[k];

        t_339[k] = -3.0 * ik_195[k]
                   + f_0 * lk_519[k];

        t_340[k] = -3.0 * ik_196[k]
                   + f_0 * lk_520[k];

        t_341[k] = -3.0 * ik_197[k]
                   + f_0 * lk_521[k];

        t_342[k] = -3.0 * ik_198[k]
                   + f_0 * lk_522[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, ik_199, ik_200, ik_201, ik_202, \
                         ik_203, lk_523, lk_524, lk_525, lk_526, \
                         lk_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = -3.0 * ik_199[k]
                   + f_0 * lk_523[k];

        t_344[k] = -3.0 * ik_200[k]
                   + f_0 * lk_524[k];

        t_345[k] = -3.0 * ik_201[k]
                   + f_0 * lk_525[k];

        t_346[k] = -3.0 * ik_202[k]
                   + f_0 * lk_526[k];

        t_347[k] = -3.0 * ik_203[k]
                   + f_0 * lk_527[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
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

    const auto *ik_204 = buffer.data(ik + 204);
    const auto *ik_205 = buffer.data(ik + 205);
    const auto *ik_206 = buffer.data(ik + 206);
    const auto *ik_207 = buffer.data(ik + 207);
    const auto *ik_208 = buffer.data(ik + 208);
    const auto *ik_209 = buffer.data(ik + 209);
    const auto *ik_210 = buffer.data(ik + 210);
    const auto *ik_211 = buffer.data(ik + 211);
    const auto *ik_212 = buffer.data(ik + 212);
    const auto *ik_213 = buffer.data(ik + 213);
    const auto *ik_214 = buffer.data(ik + 214);
    const auto *ik_215 = buffer.data(ik + 215);
    const auto *ik_216 = buffer.data(ik + 216);
    const auto *ik_217 = buffer.data(ik + 217);
    const auto *ik_218 = buffer.data(ik + 218);
    const auto *ik_219 = buffer.data(ik + 219);
    const auto *ik_220 = buffer.data(ik + 220);
    const auto *ik_221 = buffer.data(ik + 221);
    const auto *ik_222 = buffer.data(ik + 222);
    const auto *ik_223 = buffer.data(ik + 223);
    const auto *ik_224 = buffer.data(ik + 224);
    const auto *ik_225 = buffer.data(ik + 225);
    const auto *ik_226 = buffer.data(ik + 226);
    const auto *ik_227 = buffer.data(ik + 227);
    const auto *ik_228 = buffer.data(ik + 228);
    const auto *ik_229 = buffer.data(ik + 229);
    const auto *ik_230 = buffer.data(ik + 230);
    const auto *ik_231 = buffer.data(ik + 231);
    const auto *ik_232 = buffer.data(ik + 232);
    const auto *ik_233 = buffer.data(ik + 233);
    const auto *ik_234 = buffer.data(ik + 234);
    const auto *ik_235 = buffer.data(ik + 235);
    const auto *ik_236 = buffer.data(ik + 236);
    const auto *ik_237 = buffer.data(ik + 237);
    const auto *ik_238 = buffer.data(ik + 238);
    const auto *ik_239 = buffer.data(ik + 239);
    const auto *ik_240 = buffer.data(ik + 240);
    const auto *ik_241 = buffer.data(ik + 241);
    const auto *ik_242 = buffer.data(ik + 242);
    const auto *ik_243 = buffer.data(ik + 243);
    const auto *ik_244 = buffer.data(ik + 244);
    const auto *ik_245 = buffer.data(ik + 245);
    const auto *ik_246 = buffer.data(ik + 246);
    const auto *ik_247 = buffer.data(ik + 247);
    const auto *ik_248 = buffer.data(ik + 248);
    const auto *ik_249 = buffer.data(ik + 249);
    const auto *ik_250 = buffer.data(ik + 250);
    const auto *ik_251 = buffer.data(ik + 251);
    const auto *ik_252 = buffer.data(ik + 252);
    const auto *ik_253 = buffer.data(ik + 253);
    const auto *ik_254 = buffer.data(ik + 254);
    const auto *ik_255 = buffer.data(ik + 255);
    const auto *ik_256 = buffer.data(ik + 256);
    const auto *ik_257 = buffer.data(ik + 257);
    const auto *ik_258 = buffer.data(ik + 258);
    const auto *ik_259 = buffer.data(ik + 259);
    const auto *ik_260 = buffer.data(ik + 260);
    const auto *ik_261 = buffer.data(ik + 261);
    const auto *ik_262 = buffer.data(ik + 262);
    const auto *ik_263 = buffer.data(ik + 263);
    const auto *ik_264 = buffer.data(ik + 264);
    const auto *ik_265 = buffer.data(ik + 265);
    const auto *ik_266 = buffer.data(ik + 266);
    const auto *ik_267 = buffer.data(ik + 267);
    const auto *ik_268 = buffer.data(ik + 268);
    const auto *ik_269 = buffer.data(ik + 269);
    const auto *ik_270 = buffer.data(ik + 270);
    const auto *ik_271 = buffer.data(ik + 271);
    const auto *ik_272 = buffer.data(ik + 272);
    const auto *ik_273 = buffer.data(ik + 273);
    const auto *ik_274 = buffer.data(ik + 274);
    const auto *ik_275 = buffer.data(ik + 275);
    const auto *ik_276 = buffer.data(ik + 276);
    const auto *ik_277 = buffer.data(ik + 277);
    const auto *ik_278 = buffer.data(ik + 278);
    const auto *ik_279 = buffer.data(ik + 279);
    const auto *ik_280 = buffer.data(ik + 280);
    const auto *ik_281 = buffer.data(ik + 281);
    const auto *ik_282 = buffer.data(ik + 282);
    const auto *ik_283 = buffer.data(ik + 283);
    const auto *ik_284 = buffer.data(ik + 284);
    const auto *ik_285 = buffer.data(ik + 285);
    const auto *ik_286 = buffer.data(ik + 286);
    const auto *ik_287 = buffer.data(ik + 287);
    const auto *ik_288 = buffer.data(ik + 288);
    const auto *ik_289 = buffer.data(ik + 289);
    const auto *ik_290 = buffer.data(ik + 290);
    const auto *ik_291 = buffer.data(ik + 291);
    const auto *ik_292 = buffer.data(ik + 292);
    const auto *ik_293 = buffer.data(ik + 293);
    const auto *ik_294 = buffer.data(ik + 294);
    const auto *ik_295 = buffer.data(ik + 295);
    const auto *ik_296 = buffer.data(ik + 296);
    const auto *ik_297 = buffer.data(ik + 297);
    const auto *ik_298 = buffer.data(ik + 298);
    const auto *ik_299 = buffer.data(ik + 299);
    const auto *ik_300 = buffer.data(ik + 300);
    const auto *ik_301 = buffer.data(ik + 301);
    const auto *ik_302 = buffer.data(ik + 302);
    const auto *ik_303 = buffer.data(ik + 303);
    const auto *ik_304 = buffer.data(ik + 304);
    const auto *ik_305 = buffer.data(ik + 305);
    const auto *ik_306 = buffer.data(ik + 306);
    const auto *ik_307 = buffer.data(ik + 307);
    const auto *ik_308 = buffer.data(ik + 308);
    const auto *ik_309 = buffer.data(ik + 309);
    const auto *ik_310 = buffer.data(ik + 310);
    const auto *ik_311 = buffer.data(ik + 311);
    const auto *ik_312 = buffer.data(ik + 312);
    const auto *ik_313 = buffer.data(ik + 313);
    const auto *ik_314 = buffer.data(ik + 314);
    const auto *ik_315 = buffer.data(ik + 315);
    const auto *ik_316 = buffer.data(ik + 316);
    const auto *ik_317 = buffer.data(ik + 317);
    const auto *ik_318 = buffer.data(ik + 318);
    const auto *ik_319 = buffer.data(ik + 319);
    const auto *ik_320 = buffer.data(ik + 320);
    const auto *ik_321 = buffer.data(ik + 321);
    const auto *ik_322 = buffer.data(ik + 322);
    const auto *ik_323 = buffer.data(ik + 323);
    const auto *ik_324 = buffer.data(ik + 324);
    const auto *ik_325 = buffer.data(ik + 325);
    const auto *ik_326 = buffer.data(ik + 326);
    const auto *ik_327 = buffer.data(ik + 327);
    const auto *ik_328 = buffer.data(ik + 328);
    const auto *ik_329 = buffer.data(ik + 329);
    const auto *ik_330 = buffer.data(ik + 330);

    const auto *lk_528 = buffer.data(lk + 528);
    const auto *lk_529 = buffer.data(lk + 529);
    const auto *lk_530 = buffer.data(lk + 530);
    const auto *lk_531 = buffer.data(lk + 531);
    const auto *lk_532 = buffer.data(lk + 532);
    const auto *lk_533 = buffer.data(lk + 533);
    const auto *lk_534 = buffer.data(lk + 534);
    const auto *lk_535 = buffer.data(lk + 535);
    const auto *lk_536 = buffer.data(lk + 536);
    const auto *lk_537 = buffer.data(lk + 537);
    const auto *lk_538 = buffer.data(lk + 538);
    const auto *lk_539 = buffer.data(lk + 539);
    const auto *lk_576 = buffer.data(lk + 576);
    const auto *lk_577 = buffer.data(lk + 577);
    const auto *lk_578 = buffer.data(lk + 578);
    const auto *lk_579 = buffer.data(lk + 579);
    const auto *lk_580 = buffer.data(lk + 580);
    const auto *lk_581 = buffer.data(lk + 581);
    const auto *lk_582 = buffer.data(lk + 582);
    const auto *lk_583 = buffer.data(lk + 583);
    const auto *lk_584 = buffer.data(lk + 584);
    const auto *lk_585 = buffer.data(lk + 585);
    const auto *lk_586 = buffer.data(lk + 586);
    const auto *lk_587 = buffer.data(lk + 587);
    const auto *lk_588 = buffer.data(lk + 588);
    const auto *lk_589 = buffer.data(lk + 589);
    const auto *lk_590 = buffer.data(lk + 590);
    const auto *lk_591 = buffer.data(lk + 591);
    const auto *lk_592 = buffer.data(lk + 592);
    const auto *lk_593 = buffer.data(lk + 593);
    const auto *lk_594 = buffer.data(lk + 594);
    const auto *lk_595 = buffer.data(lk + 595);
    const auto *lk_596 = buffer.data(lk + 596);
    const auto *lk_597 = buffer.data(lk + 597);
    const auto *lk_598 = buffer.data(lk + 598);
    const auto *lk_599 = buffer.data(lk + 599);
    const auto *lk_600 = buffer.data(lk + 600);
    const auto *lk_601 = buffer.data(lk + 601);
    const auto *lk_602 = buffer.data(lk + 602);
    const auto *lk_603 = buffer.data(lk + 603);
    const auto *lk_604 = buffer.data(lk + 604);
    const auto *lk_605 = buffer.data(lk + 605);
    const auto *lk_606 = buffer.data(lk + 606);
    const auto *lk_607 = buffer.data(lk + 607);
    const auto *lk_608 = buffer.data(lk + 608);
    const auto *lk_609 = buffer.data(lk + 609);
    const auto *lk_610 = buffer.data(lk + 610);
    const auto *lk_611 = buffer.data(lk + 611);
    const auto *lk_612 = buffer.data(lk + 612);
    const auto *lk_613 = buffer.data(lk + 613);
    const auto *lk_614 = buffer.data(lk + 614);
    const auto *lk_615 = buffer.data(lk + 615);
    const auto *lk_616 = buffer.data(lk + 616);
    const auto *lk_617 = buffer.data(lk + 617);
    const auto *lk_618 = buffer.data(lk + 618);
    const auto *lk_619 = buffer.data(lk + 619);
    const auto *lk_620 = buffer.data(lk + 620);
    const auto *lk_621 = buffer.data(lk + 621);
    const auto *lk_622 = buffer.data(lk + 622);
    const auto *lk_623 = buffer.data(lk + 623);
    const auto *lk_624 = buffer.data(lk + 624);
    const auto *lk_625 = buffer.data(lk + 625);
    const auto *lk_626 = buffer.data(lk + 626);
    const auto *lk_627 = buffer.data(lk + 627);
    const auto *lk_628 = buffer.data(lk + 628);
    const auto *lk_629 = buffer.data(lk + 629);
    const auto *lk_630 = buffer.data(lk + 630);
    const auto *lk_631 = buffer.data(lk + 631);
    const auto *lk_632 = buffer.data(lk + 632);
    const auto *lk_633 = buffer.data(lk + 633);
    const auto *lk_634 = buffer.data(lk + 634);
    const auto *lk_635 = buffer.data(lk + 635);
    const auto *lk_636 = buffer.data(lk + 636);
    const auto *lk_637 = buffer.data(lk + 637);
    const auto *lk_638 = buffer.data(lk + 638);
    const auto *lk_639 = buffer.data(lk + 639);
    const auto *lk_640 = buffer.data(lk + 640);
    const auto *lk_641 = buffer.data(lk + 641);
    const auto *lk_642 = buffer.data(lk + 642);
    const auto *lk_643 = buffer.data(lk + 643);
    const auto *lk_644 = buffer.data(lk + 644);
    const auto *lk_645 = buffer.data(lk + 645);
    const auto *lk_646 = buffer.data(lk + 646);
    const auto *lk_647 = buffer.data(lk + 647);
    const auto *lk_648 = buffer.data(lk + 648);
    const auto *lk_649 = buffer.data(lk + 649);
    const auto *lk_650 = buffer.data(lk + 650);
    const auto *lk_651 = buffer.data(lk + 651);
    const auto *lk_652 = buffer.data(lk + 652);
    const auto *lk_653 = buffer.data(lk + 653);
    const auto *lk_654 = buffer.data(lk + 654);
    const auto *lk_655 = buffer.data(lk + 655);
    const auto *lk_656 = buffer.data(lk + 656);
    const auto *lk_657 = buffer.data(lk + 657);
    const auto *lk_658 = buffer.data(lk + 658);
    const auto *lk_659 = buffer.data(lk + 659);
    const auto *lk_660 = buffer.data(lk + 660);
    const auto *lk_661 = buffer.data(lk + 661);
    const auto *lk_662 = buffer.data(lk + 662);
    const auto *lk_663 = buffer.data(lk + 663);
    const auto *lk_664 = buffer.data(lk + 664);
    const auto *lk_665 = buffer.data(lk + 665);
    const auto *lk_666 = buffer.data(lk + 666);
    const auto *lk_667 = buffer.data(lk + 667);
    const auto *lk_668 = buffer.data(lk + 668);
    const auto *lk_669 = buffer.data(lk + 669);
    const auto *lk_670 = buffer.data(lk + 670);
    const auto *lk_671 = buffer.data(lk + 671);
    const auto *lk_672 = buffer.data(lk + 672);
    const auto *lk_673 = buffer.data(lk + 673);
    const auto *lk_674 = buffer.data(lk + 674);
    const auto *lk_675 = buffer.data(lk + 675);
    const auto *lk_676 = buffer.data(lk + 676);
    const auto *lk_677 = buffer.data(lk + 677);
    const auto *lk_678 = buffer.data(lk + 678);
    const auto *lk_679 = buffer.data(lk + 679);
    const auto *lk_680 = buffer.data(lk + 680);
    const auto *lk_681 = buffer.data(lk + 681);
    const auto *lk_682 = buffer.data(lk + 682);
    const auto *lk_683 = buffer.data(lk + 683);
    const auto *lk_684 = buffer.data(lk + 684);
    const auto *lk_685 = buffer.data(lk + 685);
    const auto *lk_686 = buffer.data(lk + 686);
    const auto *lk_687 = buffer.data(lk + 687);
    const auto *lk_688 = buffer.data(lk + 688);
    const auto *lk_689 = buffer.data(lk + 689);
    const auto *lk_690 = buffer.data(lk + 690);
    const auto *lk_691 = buffer.data(lk + 691);
    const auto *lk_692 = buffer.data(lk + 692);
    const auto *lk_693 = buffer.data(lk + 693);
    const auto *lk_694 = buffer.data(lk + 694);
    const auto *lk_695 = buffer.data(lk + 695);
    const auto *lk_696 = buffer.data(lk + 696);
    const auto *lk_697 = buffer.data(lk + 697);
    const auto *lk_698 = buffer.data(lk + 698);
    const auto *lk_699 = buffer.data(lk + 699);
    const auto *lk_700 = buffer.data(lk + 700);
    const auto *lk_701 = buffer.data(lk + 701);
    const auto *lk_702 = buffer.data(lk + 702);
    const auto *lk_703 = buffer.data(lk + 703);
    const auto *lk_704 = buffer.data(lk + 704);
    const auto *lk_705 = buffer.data(lk + 705);
    const auto *lk_706 = buffer.data(lk + 706);
    const auto *lk_707 = buffer.data(lk + 707);
    const auto *lk_708 = buffer.data(lk + 708);
    const auto *lk_709 = buffer.data(lk + 709);
    const auto *lk_710 = buffer.data(lk + 710);
    const auto *lk_711 = buffer.data(lk + 711);
    const auto *lk_712 = buffer.data(lk + 712);
    const auto *lk_713 = buffer.data(lk + 713);
    const auto *lk_714 = buffer.data(lk + 714);
    const auto *lk_715 = buffer.data(lk + 715);
    const auto *lk_716 = buffer.data(lk + 716);
    const auto *lk_717 = buffer.data(lk + 717);
    const auto *lk_718 = buffer.data(lk + 718);
    const auto *lk_719 = buffer.data(lk + 719);
    const auto *lk_720 = buffer.data(lk + 720);
    const auto *lk_721 = buffer.data(lk + 721);
    const auto *lk_722 = buffer.data(lk + 722);
    const auto *lk_723 = buffer.data(lk + 723);
    const auto *lk_724 = buffer.data(lk + 724);
    const auto *lk_725 = buffer.data(lk + 725);
    const auto *lk_726 = buffer.data(lk + 726);

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, ik_204, ik_205, ik_206, ik_207, \
                         ik_208, lk_528, lk_529, lk_530, lk_531, \
                         lk_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = -3.0 * ik_204[k]
                   + f_0 * lk_528[k];

        t_349[k] = -3.0 * ik_205[k]
                   + f_0 * lk_529[k];

        t_350[k] = -3.0 * ik_206[k]
                   + f_0 * lk_530[k];

        t_351[k] = -3.0 * ik_207[k]
                   + f_0 * lk_531[k];

        t_352[k] = -3.0 * ik_208[k]
                   + f_0 * lk_532[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, ik_209, ik_210, ik_211, ik_212, \
                         ik_213, lk_533, lk_534, lk_535, lk_536, \
                         lk_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = -3.0 * ik_209[k]
                   + f_0 * lk_533[k];

        t_354[k] = -3.0 * ik_210[k]
                   + f_0 * lk_534[k];

        t_355[k] = -3.0 * ik_211[k]
                   + f_0 * lk_535[k];

        t_356[k] = -3.0 * ik_212[k]
                   + f_0 * lk_536[k];

        t_357[k] = -3.0 * ik_213[k]
                   + f_0 * lk_537[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, t_363, t_364, ik_214, ik_215, \
                         lk_538, lk_539, lk_576, lk_577, lk_578, lk_579, \
                         lk_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = -3.0 * ik_214[k]
                   + f_0 * lk_538[k];

        t_359[k] = -3.0 * ik_215[k]
                   + f_0 * lk_539[k];

        t_360[k] = f_0 * lk_576[k];

        t_361[k] = f_0 * lk_577[k];

        t_362[k] = f_0 * lk_578[k];

        t_363[k] = f_0 * lk_579[k];

        t_364[k] = f_0 * lk_580[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, t_370, t_371, t_372, lk_581, \
                         lk_582, lk_583, lk_584, lk_585, lk_586, lk_587, \
                         lk_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_0 * lk_581[k];

        t_366[k] = f_0 * lk_582[k];

        t_367[k] = f_0 * lk_583[k];

        t_368[k] = f_0 * lk_584[k];

        t_369[k] = f_0 * lk_585[k];

        t_370[k] = f_0 * lk_586[k];

        t_371[k] = f_0 * lk_587[k];

        t_372[k] = f_0 * lk_588[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, t_378, t_379, t_380, lk_589, \
                         lk_590, lk_591, lk_592, lk_593, lk_594, lk_595, \
                         lk_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_0 * lk_589[k];

        t_374[k] = f_0 * lk_590[k];

        t_375[k] = f_0 * lk_591[k];

        t_376[k] = f_0 * lk_592[k];

        t_377[k] = f_0 * lk_593[k];

        t_378[k] = f_0 * lk_594[k];

        t_379[k] = f_0 * lk_595[k];

        t_380[k] = f_0 * lk_596[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, t_385, t_386, t_387, t_388, lk_597, \
                         lk_598, lk_599, lk_600, lk_601, lk_602, lk_603, \
                         lk_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_0 * lk_597[k];

        t_382[k] = f_0 * lk_598[k];

        t_383[k] = f_0 * lk_599[k];

        t_384[k] = f_0 * lk_600[k];

        t_385[k] = f_0 * lk_601[k];

        t_386[k] = f_0 * lk_602[k];

        t_387[k] = f_0 * lk_603[k];

        t_388[k] = f_0 * lk_604[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, t_394, t_395, lk_605, lk_606, \
                         lk_607, lk_608, lk_609, lk_610, lk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_0 * lk_605[k];

        t_390[k] = f_0 * lk_606[k];

        t_391[k] = f_0 * lk_607[k];

        t_392[k] = f_0 * lk_608[k];

        t_393[k] = f_0 * lk_609[k];

        t_394[k] = f_0 * lk_610[k];

        t_395[k] = f_0 * lk_611[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, t_400, ik_216, ik_217, ik_218, ik_219, \
                         ik_220, lk_612, lk_613, lk_614, lk_615, \
                         lk_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = -ik_216[k]
                   + f_0 * lk_612[k];

        t_397[k] = -ik_217[k]
                   + f_0 * lk_613[k];

        t_398[k] = -ik_218[k]
                   + f_0 * lk_614[k];

        t_399[k] = -ik_219[k]
                   + f_0 * lk_615[k];

        t_400[k] = -ik_220[k]
                   + f_0 * lk_616[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, t_405, ik_221, ik_222, ik_223, ik_224, \
                         ik_225, lk_617, lk_618, lk_619, lk_620, \
                         lk_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = -ik_221[k]
                   + f_0 * lk_617[k];

        t_402[k] = -ik_222[k]
                   + f_0 * lk_618[k];

        t_403[k] = -ik_223[k]
                   + f_0 * lk_619[k];

        t_404[k] = -ik_224[k]
                   + f_0 * lk_620[k];

        t_405[k] = -ik_225[k]
                   + f_0 * lk_621[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, t_410, ik_226, ik_227, ik_228, ik_229, \
                         ik_230, lk_622, lk_623, lk_624, lk_625, \
                         lk_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = -ik_226[k]
                   + f_0 * lk_622[k];

        t_407[k] = -ik_227[k]
                   + f_0 * lk_623[k];

        t_408[k] = -ik_228[k]
                   + f_0 * lk_624[k];

        t_409[k] = -ik_229[k]
                   + f_0 * lk_625[k];

        t_410[k] = -ik_230[k]
                   + f_0 * lk_626[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, t_415, ik_231, ik_232, ik_233, ik_234, \
                         ik_235, lk_627, lk_628, lk_629, lk_630, \
                         lk_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = -ik_231[k]
                   + f_0 * lk_627[k];

        t_412[k] = -ik_232[k]
                   + f_0 * lk_628[k];

        t_413[k] = -ik_233[k]
                   + f_0 * lk_629[k];

        t_414[k] = -ik_234[k]
                   + f_0 * lk_630[k];

        t_415[k] = -ik_235[k]
                   + f_0 * lk_631[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, t_420, ik_236, ik_237, ik_238, ik_239, \
                         ik_240, lk_632, lk_633, lk_634, lk_635, \
                         lk_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = -ik_236[k]
                   + f_0 * lk_632[k];

        t_417[k] = -ik_237[k]
                   + f_0 * lk_633[k];

        t_418[k] = -ik_238[k]
                   + f_0 * lk_634[k];

        t_419[k] = -ik_239[k]
                   + f_0 * lk_635[k];

        t_420[k] = -ik_240[k]
                   + f_0 * lk_636[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, t_425, ik_241, ik_242, ik_243, ik_244, \
                         ik_245, lk_637, lk_638, lk_639, lk_640, \
                         lk_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = -ik_241[k]
                   + f_0 * lk_637[k];

        t_422[k] = -ik_242[k]
                   + f_0 * lk_638[k];

        t_423[k] = -ik_243[k]
                   + f_0 * lk_639[k];

        t_424[k] = -ik_244[k]
                   + f_0 * lk_640[k];

        t_425[k] = -ik_245[k]
                   + f_0 * lk_641[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, t_430, ik_246, ik_247, ik_248, ik_249, \
                         ik_250, lk_642, lk_643, lk_644, lk_645, \
                         lk_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = -ik_246[k]
                   + f_0 * lk_642[k];

        t_427[k] = -ik_247[k]
                   + f_0 * lk_643[k];

        t_428[k] = -ik_248[k]
                   + f_0 * lk_644[k];

        t_429[k] = -ik_249[k]
                   + f_0 * lk_645[k];

        t_430[k] = -ik_250[k]
                   + f_0 * lk_646[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, ik_251, ik_252, ik_253, ik_254, \
                         ik_255, lk_647, lk_648, lk_649, lk_650, \
                         lk_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = -ik_251[k]
                   + f_0 * lk_647[k];

        t_432[k] = -2.0 * ik_252[k]
                   + f_0 * lk_648[k];

        t_433[k] = -2.0 * ik_253[k]
                   + f_0 * lk_649[k];

        t_434[k] = -2.0 * ik_254[k]
                   + f_0 * lk_650[k];

        t_435[k] = -2.0 * ik_255[k]
                   + f_0 * lk_651[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, ik_256, ik_257, ik_258, ik_259, \
                         ik_260, lk_652, lk_653, lk_654, lk_655, \
                         lk_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = -2.0 * ik_256[k]
                   + f_0 * lk_652[k];

        t_437[k] = -2.0 * ik_257[k]
                   + f_0 * lk_653[k];

        t_438[k] = -2.0 * ik_258[k]
                   + f_0 * lk_654[k];

        t_439[k] = -2.0 * ik_259[k]
                   + f_0 * lk_655[k];

        t_440[k] = -2.0 * ik_260[k]
                   + f_0 * lk_656[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, t_445, ik_261, ik_262, ik_263, ik_264, \
                         ik_265, lk_657, lk_658, lk_659, lk_660, \
                         lk_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = -2.0 * ik_261[k]
                   + f_0 * lk_657[k];

        t_442[k] = -2.0 * ik_262[k]
                   + f_0 * lk_658[k];

        t_443[k] = -2.0 * ik_263[k]
                   + f_0 * lk_659[k];

        t_444[k] = -2.0 * ik_264[k]
                   + f_0 * lk_660[k];

        t_445[k] = -2.0 * ik_265[k]
                   + f_0 * lk_661[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, t_450, ik_266, ik_267, ik_268, ik_269, \
                         ik_270, lk_662, lk_663, lk_664, lk_665, \
                         lk_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = -2.0 * ik_266[k]
                   + f_0 * lk_662[k];

        t_447[k] = -2.0 * ik_267[k]
                   + f_0 * lk_663[k];

        t_448[k] = -2.0 * ik_268[k]
                   + f_0 * lk_664[k];

        t_449[k] = -2.0 * ik_269[k]
                   + f_0 * lk_665[k];

        t_450[k] = -2.0 * ik_270[k]
                   + f_0 * lk_666[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, t_455, ik_271, ik_272, ik_273, ik_274, \
                         ik_275, lk_667, lk_668, lk_669, lk_670, \
                         lk_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = -2.0 * ik_271[k]
                   + f_0 * lk_667[k];

        t_452[k] = -2.0 * ik_272[k]
                   + f_0 * lk_668[k];

        t_453[k] = -2.0 * ik_273[k]
                   + f_0 * lk_669[k];

        t_454[k] = -2.0 * ik_274[k]
                   + f_0 * lk_670[k];

        t_455[k] = -2.0 * ik_275[k]
                   + f_0 * lk_671[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, t_460, ik_276, ik_277, ik_278, ik_279, \
                         ik_280, lk_672, lk_673, lk_674, lk_675, \
                         lk_676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = -2.0 * ik_276[k]
                   + f_0 * lk_672[k];

        t_457[k] = -2.0 * ik_277[k]
                   + f_0 * lk_673[k];

        t_458[k] = -2.0 * ik_278[k]
                   + f_0 * lk_674[k];

        t_459[k] = -2.0 * ik_279[k]
                   + f_0 * lk_675[k];

        t_460[k] = -2.0 * ik_280[k]
                   + f_0 * lk_676[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, t_465, ik_281, ik_282, ik_283, ik_284, \
                         ik_285, lk_677, lk_678, lk_679, lk_680, \
                         lk_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = -2.0 * ik_281[k]
                   + f_0 * lk_677[k];

        t_462[k] = -2.0 * ik_282[k]
                   + f_0 * lk_678[k];

        t_463[k] = -2.0 * ik_283[k]
                   + f_0 * lk_679[k];

        t_464[k] = -2.0 * ik_284[k]
                   + f_0 * lk_680[k];

        t_465[k] = -2.0 * ik_285[k]
                   + f_0 * lk_681[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, ik_286, ik_287, ik_288, ik_289, \
                         ik_290, lk_682, lk_683, lk_684, lk_685, \
                         lk_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = -2.0 * ik_286[k]
                   + f_0 * lk_682[k];

        t_467[k] = -2.0 * ik_287[k]
                   + f_0 * lk_683[k];

        t_468[k] = -3.0 * ik_288[k]
                   + f_0 * lk_684[k];

        t_469[k] = -3.0 * ik_289[k]
                   + f_0 * lk_685[k];

        t_470[k] = -3.0 * ik_290[k]
                   + f_0 * lk_686[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, t_475, ik_291, ik_292, ik_293, ik_294, \
                         ik_295, lk_687, lk_688, lk_689, lk_690, \
                         lk_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = -3.0 * ik_291[k]
                   + f_0 * lk_687[k];

        t_472[k] = -3.0 * ik_292[k]
                   + f_0 * lk_688[k];

        t_473[k] = -3.0 * ik_293[k]
                   + f_0 * lk_689[k];

        t_474[k] = -3.0 * ik_294[k]
                   + f_0 * lk_690[k];

        t_475[k] = -3.0 * ik_295[k]
                   + f_0 * lk_691[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, ik_296, ik_297, ik_298, ik_299, \
                         ik_300, lk_692, lk_693, lk_694, lk_695, \
                         lk_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = -3.0 * ik_296[k]
                   + f_0 * lk_692[k];

        t_477[k] = -3.0 * ik_297[k]
                   + f_0 * lk_693[k];

        t_478[k] = -3.0 * ik_298[k]
                   + f_0 * lk_694[k];

        t_479[k] = -3.0 * ik_299[k]
                   + f_0 * lk_695[k];

        t_480[k] = -3.0 * ik_300[k]
                   + f_0 * lk_696[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, ik_301, ik_302, ik_303, ik_304, \
                         ik_305, lk_697, lk_698, lk_699, lk_700, \
                         lk_701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = -3.0 * ik_301[k]
                   + f_0 * lk_697[k];

        t_482[k] = -3.0 * ik_302[k]
                   + f_0 * lk_698[k];

        t_483[k] = -3.0 * ik_303[k]
                   + f_0 * lk_699[k];

        t_484[k] = -3.0 * ik_304[k]
                   + f_0 * lk_700[k];

        t_485[k] = -3.0 * ik_305[k]
                   + f_0 * lk_701[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, ik_306, ik_307, ik_308, ik_309, \
                         ik_310, lk_702, lk_703, lk_704, lk_705, \
                         lk_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = -3.0 * ik_306[k]
                   + f_0 * lk_702[k];

        t_487[k] = -3.0 * ik_307[k]
                   + f_0 * lk_703[k];

        t_488[k] = -3.0 * ik_308[k]
                   + f_0 * lk_704[k];

        t_489[k] = -3.0 * ik_309[k]
                   + f_0 * lk_705[k];

        t_490[k] = -3.0 * ik_310[k]
                   + f_0 * lk_706[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, ik_311, ik_312, ik_313, ik_314, \
                         ik_315, lk_707, lk_708, lk_709, lk_710, \
                         lk_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = -3.0 * ik_311[k]
                   + f_0 * lk_707[k];

        t_492[k] = -3.0 * ik_312[k]
                   + f_0 * lk_708[k];

        t_493[k] = -3.0 * ik_313[k]
                   + f_0 * lk_709[k];

        t_494[k] = -3.0 * ik_314[k]
                   + f_0 * lk_710[k];

        t_495[k] = -3.0 * ik_315[k]
                   + f_0 * lk_711[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, t_500, ik_316, ik_317, ik_318, ik_319, \
                         ik_320, lk_712, lk_713, lk_714, lk_715, \
                         lk_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = -3.0 * ik_316[k]
                   + f_0 * lk_712[k];

        t_497[k] = -3.0 * ik_317[k]
                   + f_0 * lk_713[k];

        t_498[k] = -3.0 * ik_318[k]
                   + f_0 * lk_714[k];

        t_499[k] = -3.0 * ik_319[k]
                   + f_0 * lk_715[k];

        t_500[k] = -3.0 * ik_320[k]
                   + f_0 * lk_716[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, t_505, ik_321, ik_322, ik_323, ik_324, \
                         ik_325, lk_717, lk_718, lk_719, lk_720, \
                         lk_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = -3.0 * ik_321[k]
                   + f_0 * lk_717[k];

        t_502[k] = -3.0 * ik_322[k]
                   + f_0 * lk_718[k];

        t_503[k] = -3.0 * ik_323[k]
                   + f_0 * lk_719[k];

        t_504[k] = -4.0 * ik_324[k]
                   + f_0 * lk_720[k];

        t_505[k] = -4.0 * ik_325[k]
                   + f_0 * lk_721[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, t_509, t_510, ik_326, ik_327, ik_328, ik_329, \
                         ik_330, lk_722, lk_723, lk_724, lk_725, \
                         lk_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = -4.0 * ik_326[k]
                   + f_0 * lk_722[k];

        t_507[k] = -4.0 * ik_327[k]
                   + f_0 * lk_723[k];

        t_508[k] = -4.0 * ik_328[k]
                   + f_0 * lk_724[k];

        t_509[k] = -4.0 * ik_329[k]
                   + f_0 * lk_725[k];

        t_510[k] = -4.0 * ik_330[k]
                   + f_0 * lk_726[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ik_331 = buffer.data(ik + 331);
    const auto *ik_332 = buffer.data(ik + 332);
    const auto *ik_333 = buffer.data(ik + 333);
    const auto *ik_334 = buffer.data(ik + 334);
    const auto *ik_335 = buffer.data(ik + 335);
    const auto *ik_336 = buffer.data(ik + 336);
    const auto *ik_337 = buffer.data(ik + 337);
    const auto *ik_338 = buffer.data(ik + 338);
    const auto *ik_339 = buffer.data(ik + 339);
    const auto *ik_340 = buffer.data(ik + 340);
    const auto *ik_341 = buffer.data(ik + 341);
    const auto *ik_342 = buffer.data(ik + 342);
    const auto *ik_343 = buffer.data(ik + 343);
    const auto *ik_344 = buffer.data(ik + 344);
    const auto *ik_345 = buffer.data(ik + 345);
    const auto *ik_346 = buffer.data(ik + 346);
    const auto *ik_347 = buffer.data(ik + 347);
    const auto *ik_348 = buffer.data(ik + 348);
    const auto *ik_349 = buffer.data(ik + 349);
    const auto *ik_350 = buffer.data(ik + 350);
    const auto *ik_351 = buffer.data(ik + 351);
    const auto *ik_352 = buffer.data(ik + 352);
    const auto *ik_353 = buffer.data(ik + 353);
    const auto *ik_354 = buffer.data(ik + 354);
    const auto *ik_355 = buffer.data(ik + 355);
    const auto *ik_356 = buffer.data(ik + 356);
    const auto *ik_357 = buffer.data(ik + 357);
    const auto *ik_358 = buffer.data(ik + 358);
    const auto *ik_359 = buffer.data(ik + 359);
    const auto *ik_360 = buffer.data(ik + 360);
    const auto *ik_361 = buffer.data(ik + 361);
    const auto *ik_362 = buffer.data(ik + 362);
    const auto *ik_363 = buffer.data(ik + 363);
    const auto *ik_364 = buffer.data(ik + 364);
    const auto *ik_365 = buffer.data(ik + 365);
    const auto *ik_366 = buffer.data(ik + 366);
    const auto *ik_367 = buffer.data(ik + 367);
    const auto *ik_368 = buffer.data(ik + 368);
    const auto *ik_369 = buffer.data(ik + 369);
    const auto *ik_370 = buffer.data(ik + 370);
    const auto *ik_371 = buffer.data(ik + 371);
    const auto *ik_372 = buffer.data(ik + 372);
    const auto *ik_373 = buffer.data(ik + 373);
    const auto *ik_374 = buffer.data(ik + 374);
    const auto *ik_375 = buffer.data(ik + 375);
    const auto *ik_376 = buffer.data(ik + 376);
    const auto *ik_377 = buffer.data(ik + 377);
    const auto *ik_378 = buffer.data(ik + 378);
    const auto *ik_379 = buffer.data(ik + 379);
    const auto *ik_380 = buffer.data(ik + 380);
    const auto *ik_381 = buffer.data(ik + 381);
    const auto *ik_382 = buffer.data(ik + 382);
    const auto *ik_383 = buffer.data(ik + 383);
    const auto *ik_384 = buffer.data(ik + 384);
    const auto *ik_385 = buffer.data(ik + 385);
    const auto *ik_386 = buffer.data(ik + 386);
    const auto *ik_387 = buffer.data(ik + 387);
    const auto *ik_388 = buffer.data(ik + 388);
    const auto *ik_389 = buffer.data(ik + 389);
    const auto *ik_390 = buffer.data(ik + 390);
    const auto *ik_391 = buffer.data(ik + 391);
    const auto *ik_392 = buffer.data(ik + 392);
    const auto *ik_393 = buffer.data(ik + 393);
    const auto *ik_394 = buffer.data(ik + 394);
    const auto *ik_395 = buffer.data(ik + 395);
    const auto *ik_396 = buffer.data(ik + 396);
    const auto *ik_397 = buffer.data(ik + 397);
    const auto *ik_398 = buffer.data(ik + 398);
    const auto *ik_399 = buffer.data(ik + 399);
    const auto *ik_400 = buffer.data(ik + 400);
    const auto *ik_401 = buffer.data(ik + 401);
    const auto *ik_402 = buffer.data(ik + 402);
    const auto *ik_403 = buffer.data(ik + 403);
    const auto *ik_404 = buffer.data(ik + 404);
    const auto *ik_405 = buffer.data(ik + 405);
    const auto *ik_406 = buffer.data(ik + 406);
    const auto *ik_407 = buffer.data(ik + 407);
    const auto *ik_408 = buffer.data(ik + 408);
    const auto *ik_409 = buffer.data(ik + 409);
    const auto *ik_410 = buffer.data(ik + 410);
    const auto *ik_411 = buffer.data(ik + 411);
    const auto *ik_412 = buffer.data(ik + 412);
    const auto *ik_413 = buffer.data(ik + 413);
    const auto *ik_414 = buffer.data(ik + 414);
    const auto *ik_415 = buffer.data(ik + 415);
    const auto *ik_416 = buffer.data(ik + 416);
    const auto *ik_417 = buffer.data(ik + 417);
    const auto *ik_418 = buffer.data(ik + 418);
    const auto *ik_419 = buffer.data(ik + 419);
    const auto *ik_420 = buffer.data(ik + 420);
    const auto *ik_421 = buffer.data(ik + 421);
    const auto *ik_422 = buffer.data(ik + 422);
    const auto *ik_423 = buffer.data(ik + 423);
    const auto *ik_424 = buffer.data(ik + 424);
    const auto *ik_425 = buffer.data(ik + 425);
    const auto *ik_426 = buffer.data(ik + 426);
    const auto *ik_427 = buffer.data(ik + 427);
    const auto *ik_428 = buffer.data(ik + 428);
    const auto *ik_429 = buffer.data(ik + 429);
    const auto *ik_430 = buffer.data(ik + 430);
    const auto *ik_431 = buffer.data(ik + 431);
    const auto *ik_432 = buffer.data(ik + 432);
    const auto *ik_433 = buffer.data(ik + 433);
    const auto *ik_434 = buffer.data(ik + 434);
    const auto *ik_435 = buffer.data(ik + 435);
    const auto *ik_436 = buffer.data(ik + 436);
    const auto *ik_437 = buffer.data(ik + 437);
    const auto *ik_438 = buffer.data(ik + 438);
    const auto *ik_439 = buffer.data(ik + 439);
    const auto *ik_440 = buffer.data(ik + 440);
    const auto *ik_441 = buffer.data(ik + 441);
    const auto *ik_442 = buffer.data(ik + 442);
    const auto *ik_443 = buffer.data(ik + 443);
    const auto *ik_444 = buffer.data(ik + 444);
    const auto *ik_445 = buffer.data(ik + 445);
    const auto *ik_446 = buffer.data(ik + 446);
    const auto *ik_447 = buffer.data(ik + 447);
    const auto *ik_448 = buffer.data(ik + 448);
    const auto *ik_449 = buffer.data(ik + 449);
    const auto *ik_450 = buffer.data(ik + 450);
    const auto *ik_451 = buffer.data(ik + 451);
    const auto *ik_452 = buffer.data(ik + 452);
    const auto *ik_453 = buffer.data(ik + 453);

    const auto *lk_727 = buffer.data(lk + 727);
    const auto *lk_728 = buffer.data(lk + 728);
    const auto *lk_729 = buffer.data(lk + 729);
    const auto *lk_730 = buffer.data(lk + 730);
    const auto *lk_731 = buffer.data(lk + 731);
    const auto *lk_732 = buffer.data(lk + 732);
    const auto *lk_733 = buffer.data(lk + 733);
    const auto *lk_734 = buffer.data(lk + 734);
    const auto *lk_735 = buffer.data(lk + 735);
    const auto *lk_736 = buffer.data(lk + 736);
    const auto *lk_737 = buffer.data(lk + 737);
    const auto *lk_738 = buffer.data(lk + 738);
    const auto *lk_739 = buffer.data(lk + 739);
    const auto *lk_740 = buffer.data(lk + 740);
    const auto *lk_741 = buffer.data(lk + 741);
    const auto *lk_742 = buffer.data(lk + 742);
    const auto *lk_743 = buffer.data(lk + 743);
    const auto *lk_744 = buffer.data(lk + 744);
    const auto *lk_745 = buffer.data(lk + 745);
    const auto *lk_746 = buffer.data(lk + 746);
    const auto *lk_747 = buffer.data(lk + 747);
    const auto *lk_748 = buffer.data(lk + 748);
    const auto *lk_749 = buffer.data(lk + 749);
    const auto *lk_750 = buffer.data(lk + 750);
    const auto *lk_751 = buffer.data(lk + 751);
    const auto *lk_752 = buffer.data(lk + 752);
    const auto *lk_753 = buffer.data(lk + 753);
    const auto *lk_754 = buffer.data(lk + 754);
    const auto *lk_755 = buffer.data(lk + 755);
    const auto *lk_792 = buffer.data(lk + 792);
    const auto *lk_793 = buffer.data(lk + 793);
    const auto *lk_794 = buffer.data(lk + 794);
    const auto *lk_795 = buffer.data(lk + 795);
    const auto *lk_796 = buffer.data(lk + 796);
    const auto *lk_797 = buffer.data(lk + 797);
    const auto *lk_798 = buffer.data(lk + 798);
    const auto *lk_799 = buffer.data(lk + 799);
    const auto *lk_800 = buffer.data(lk + 800);
    const auto *lk_801 = buffer.data(lk + 801);
    const auto *lk_802 = buffer.data(lk + 802);
    const auto *lk_803 = buffer.data(lk + 803);
    const auto *lk_804 = buffer.data(lk + 804);
    const auto *lk_805 = buffer.data(lk + 805);
    const auto *lk_806 = buffer.data(lk + 806);
    const auto *lk_807 = buffer.data(lk + 807);
    const auto *lk_808 = buffer.data(lk + 808);
    const auto *lk_809 = buffer.data(lk + 809);
    const auto *lk_810 = buffer.data(lk + 810);
    const auto *lk_811 = buffer.data(lk + 811);
    const auto *lk_812 = buffer.data(lk + 812);
    const auto *lk_813 = buffer.data(lk + 813);
    const auto *lk_814 = buffer.data(lk + 814);
    const auto *lk_815 = buffer.data(lk + 815);
    const auto *lk_816 = buffer.data(lk + 816);
    const auto *lk_817 = buffer.data(lk + 817);
    const auto *lk_818 = buffer.data(lk + 818);
    const auto *lk_819 = buffer.data(lk + 819);
    const auto *lk_820 = buffer.data(lk + 820);
    const auto *lk_821 = buffer.data(lk + 821);
    const auto *lk_822 = buffer.data(lk + 822);
    const auto *lk_823 = buffer.data(lk + 823);
    const auto *lk_824 = buffer.data(lk + 824);
    const auto *lk_825 = buffer.data(lk + 825);
    const auto *lk_826 = buffer.data(lk + 826);
    const auto *lk_827 = buffer.data(lk + 827);
    const auto *lk_828 = buffer.data(lk + 828);
    const auto *lk_829 = buffer.data(lk + 829);
    const auto *lk_830 = buffer.data(lk + 830);
    const auto *lk_831 = buffer.data(lk + 831);
    const auto *lk_832 = buffer.data(lk + 832);
    const auto *lk_833 = buffer.data(lk + 833);
    const auto *lk_834 = buffer.data(lk + 834);
    const auto *lk_835 = buffer.data(lk + 835);
    const auto *lk_836 = buffer.data(lk + 836);
    const auto *lk_837 = buffer.data(lk + 837);
    const auto *lk_838 = buffer.data(lk + 838);
    const auto *lk_839 = buffer.data(lk + 839);
    const auto *lk_840 = buffer.data(lk + 840);
    const auto *lk_841 = buffer.data(lk + 841);
    const auto *lk_842 = buffer.data(lk + 842);
    const auto *lk_843 = buffer.data(lk + 843);
    const auto *lk_844 = buffer.data(lk + 844);
    const auto *lk_845 = buffer.data(lk + 845);
    const auto *lk_846 = buffer.data(lk + 846);
    const auto *lk_847 = buffer.data(lk + 847);
    const auto *lk_848 = buffer.data(lk + 848);
    const auto *lk_849 = buffer.data(lk + 849);
    const auto *lk_850 = buffer.data(lk + 850);
    const auto *lk_851 = buffer.data(lk + 851);
    const auto *lk_852 = buffer.data(lk + 852);
    const auto *lk_853 = buffer.data(lk + 853);
    const auto *lk_854 = buffer.data(lk + 854);
    const auto *lk_855 = buffer.data(lk + 855);
    const auto *lk_856 = buffer.data(lk + 856);
    const auto *lk_857 = buffer.data(lk + 857);
    const auto *lk_858 = buffer.data(lk + 858);
    const auto *lk_859 = buffer.data(lk + 859);
    const auto *lk_860 = buffer.data(lk + 860);
    const auto *lk_861 = buffer.data(lk + 861);
    const auto *lk_862 = buffer.data(lk + 862);
    const auto *lk_863 = buffer.data(lk + 863);
    const auto *lk_864 = buffer.data(lk + 864);
    const auto *lk_865 = buffer.data(lk + 865);
    const auto *lk_866 = buffer.data(lk + 866);
    const auto *lk_867 = buffer.data(lk + 867);
    const auto *lk_868 = buffer.data(lk + 868);
    const auto *lk_869 = buffer.data(lk + 869);
    const auto *lk_870 = buffer.data(lk + 870);
    const auto *lk_871 = buffer.data(lk + 871);
    const auto *lk_872 = buffer.data(lk + 872);
    const auto *lk_873 = buffer.data(lk + 873);
    const auto *lk_874 = buffer.data(lk + 874);
    const auto *lk_875 = buffer.data(lk + 875);
    const auto *lk_876 = buffer.data(lk + 876);
    const auto *lk_877 = buffer.data(lk + 877);
    const auto *lk_878 = buffer.data(lk + 878);
    const auto *lk_879 = buffer.data(lk + 879);
    const auto *lk_880 = buffer.data(lk + 880);
    const auto *lk_881 = buffer.data(lk + 881);
    const auto *lk_882 = buffer.data(lk + 882);
    const auto *lk_883 = buffer.data(lk + 883);
    const auto *lk_884 = buffer.data(lk + 884);
    const auto *lk_885 = buffer.data(lk + 885);
    const auto *lk_886 = buffer.data(lk + 886);
    const auto *lk_887 = buffer.data(lk + 887);
    const auto *lk_888 = buffer.data(lk + 888);
    const auto *lk_889 = buffer.data(lk + 889);
    const auto *lk_890 = buffer.data(lk + 890);
    const auto *lk_891 = buffer.data(lk + 891);
    const auto *lk_892 = buffer.data(lk + 892);
    const auto *lk_893 = buffer.data(lk + 893);
    const auto *lk_894 = buffer.data(lk + 894);
    const auto *lk_895 = buffer.data(lk + 895);
    const auto *lk_896 = buffer.data(lk + 896);
    const auto *lk_897 = buffer.data(lk + 897);
    const auto *lk_898 = buffer.data(lk + 898);
    const auto *lk_899 = buffer.data(lk + 899);
    const auto *lk_900 = buffer.data(lk + 900);
    const auto *lk_901 = buffer.data(lk + 901);
    const auto *lk_902 = buffer.data(lk + 902);
    const auto *lk_903 = buffer.data(lk + 903);
    const auto *lk_904 = buffer.data(lk + 904);
    const auto *lk_905 = buffer.data(lk + 905);
    const auto *lk_906 = buffer.data(lk + 906);
    const auto *lk_907 = buffer.data(lk + 907);
    const auto *lk_908 = buffer.data(lk + 908);
    const auto *lk_909 = buffer.data(lk + 909);
    const auto *lk_910 = buffer.data(lk + 910);
    const auto *lk_911 = buffer.data(lk + 911);
    const auto *lk_912 = buffer.data(lk + 912);
    const auto *lk_913 = buffer.data(lk + 913);
    const auto *lk_914 = buffer.data(lk + 914);
    const auto *lk_915 = buffer.data(lk + 915);
    const auto *lk_916 = buffer.data(lk + 916);
    const auto *lk_917 = buffer.data(lk + 917);
    const auto *lk_918 = buffer.data(lk + 918);
    const auto *lk_919 = buffer.data(lk + 919);
    const auto *lk_920 = buffer.data(lk + 920);
    const auto *lk_921 = buffer.data(lk + 921);

#pragma omp simd aligned(t_511, t_512, t_513, t_514, t_515, ik_331, ik_332, ik_333, ik_334, \
                         ik_335, lk_727, lk_728, lk_729, lk_730, \
                         lk_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = -4.0 * ik_331[k]
                   + f_0 * lk_727[k];

        t_512[k] = -4.0 * ik_332[k]
                   + f_0 * lk_728[k];

        t_513[k] = -4.0 * ik_333[k]
                   + f_0 * lk_729[k];

        t_514[k] = -4.0 * ik_334[k]
                   + f_0 * lk_730[k];

        t_515[k] = -4.0 * ik_335[k]
                   + f_0 * lk_731[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, t_520, ik_336, ik_337, ik_338, ik_339, \
                         ik_340, lk_732, lk_733, lk_734, lk_735, \
                         lk_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = -4.0 * ik_336[k]
                   + f_0 * lk_732[k];

        t_517[k] = -4.0 * ik_337[k]
                   + f_0 * lk_733[k];

        t_518[k] = -4.0 * ik_338[k]
                   + f_0 * lk_734[k];

        t_519[k] = -4.0 * ik_339[k]
                   + f_0 * lk_735[k];

        t_520[k] = -4.0 * ik_340[k]
                   + f_0 * lk_736[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, t_525, ik_341, ik_342, ik_343, ik_344, \
                         ik_345, lk_737, lk_738, lk_739, lk_740, \
                         lk_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = -4.0 * ik_341[k]
                   + f_0 * lk_737[k];

        t_522[k] = -4.0 * ik_342[k]
                   + f_0 * lk_738[k];

        t_523[k] = -4.0 * ik_343[k]
                   + f_0 * lk_739[k];

        t_524[k] = -4.0 * ik_344[k]
                   + f_0 * lk_740[k];

        t_525[k] = -4.0 * ik_345[k]
                   + f_0 * lk_741[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, ik_346, ik_347, ik_348, ik_349, \
                         ik_350, lk_742, lk_743, lk_744, lk_745, \
                         lk_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = -4.0 * ik_346[k]
                   + f_0 * lk_742[k];

        t_527[k] = -4.0 * ik_347[k]
                   + f_0 * lk_743[k];

        t_528[k] = -4.0 * ik_348[k]
                   + f_0 * lk_744[k];

        t_529[k] = -4.0 * ik_349[k]
                   + f_0 * lk_745[k];

        t_530[k] = -4.0 * ik_350[k]
                   + f_0 * lk_746[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, t_535, ik_351, ik_352, ik_353, ik_354, \
                         ik_355, lk_747, lk_748, lk_749, lk_750, \
                         lk_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = -4.0 * ik_351[k]
                   + f_0 * lk_747[k];

        t_532[k] = -4.0 * ik_352[k]
                   + f_0 * lk_748[k];

        t_533[k] = -4.0 * ik_353[k]
                   + f_0 * lk_749[k];

        t_534[k] = -4.0 * ik_354[k]
                   + f_0 * lk_750[k];

        t_535[k] = -4.0 * ik_355[k]
                   + f_0 * lk_751[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, t_541, ik_356, ik_357, ik_358, \
                         ik_359, lk_752, lk_753, lk_754, lk_755, lk_792, \
                         lk_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = -4.0 * ik_356[k]
                   + f_0 * lk_752[k];

        t_537[k] = -4.0 * ik_357[k]
                   + f_0 * lk_753[k];

        t_538[k] = -4.0 * ik_358[k]
                   + f_0 * lk_754[k];

        t_539[k] = -4.0 * ik_359[k]
                   + f_0 * lk_755[k];

        t_540[k] = f_0 * lk_792[k];

        t_541[k] = f_0 * lk_793[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, t_547, t_548, t_549, lk_794, \
                         lk_795, lk_796, lk_797, lk_798, lk_799, lk_800, \
                         lk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_0 * lk_794[k];

        t_543[k] = f_0 * lk_795[k];

        t_544[k] = f_0 * lk_796[k];

        t_545[k] = f_0 * lk_797[k];

        t_546[k] = f_0 * lk_798[k];

        t_547[k] = f_0 * lk_799[k];

        t_548[k] = f_0 * lk_800[k];

        t_549[k] = f_0 * lk_801[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, t_555, t_556, t_557, lk_802, \
                         lk_803, lk_804, lk_805, lk_806, lk_807, lk_808, \
                         lk_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_0 * lk_802[k];

        t_551[k] = f_0 * lk_803[k];

        t_552[k] = f_0 * lk_804[k];

        t_553[k] = f_0 * lk_805[k];

        t_554[k] = f_0 * lk_806[k];

        t_555[k] = f_0 * lk_807[k];

        t_556[k] = f_0 * lk_808[k];

        t_557[k] = f_0 * lk_809[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, t_562, t_563, t_564, t_565, lk_810, \
                         lk_811, lk_812, lk_813, lk_814, lk_815, lk_816, \
                         lk_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_0 * lk_810[k];

        t_559[k] = f_0 * lk_811[k];

        t_560[k] = f_0 * lk_812[k];

        t_561[k] = f_0 * lk_813[k];

        t_562[k] = f_0 * lk_814[k];

        t_563[k] = f_0 * lk_815[k];

        t_564[k] = f_0 * lk_816[k];

        t_565[k] = f_0 * lk_817[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, t_570, t_571, t_572, t_573, lk_818, \
                         lk_819, lk_820, lk_821, lk_822, lk_823, lk_824, \
                         lk_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_0 * lk_818[k];

        t_567[k] = f_0 * lk_819[k];

        t_568[k] = f_0 * lk_820[k];

        t_569[k] = f_0 * lk_821[k];

        t_570[k] = f_0 * lk_822[k];

        t_571[k] = f_0 * lk_823[k];

        t_572[k] = f_0 * lk_824[k];

        t_573[k] = f_0 * lk_825[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, t_579, ik_360, ik_361, ik_362, \
                         ik_363, lk_826, lk_827, lk_828, lk_829, lk_830, \
                         lk_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_0 * lk_826[k];

        t_575[k] = f_0 * lk_827[k];

        t_576[k] = -ik_360[k]
                   + f_0 * lk_828[k];

        t_577[k] = -ik_361[k]
                   + f_0 * lk_829[k];

        t_578[k] = -ik_362[k]
                   + f_0 * lk_830[k];

        t_579[k] = -ik_363[k]
                   + f_0 * lk_831[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, ik_364, ik_365, ik_366, ik_367, \
                         ik_368, lk_832, lk_833, lk_834, lk_835, \
                         lk_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -ik_364[k]
                   + f_0 * lk_832[k];

        t_581[k] = -ik_365[k]
                   + f_0 * lk_833[k];

        t_582[k] = -ik_366[k]
                   + f_0 * lk_834[k];

        t_583[k] = -ik_367[k]
                   + f_0 * lk_835[k];

        t_584[k] = -ik_368[k]
                   + f_0 * lk_836[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, ik_369, ik_370, ik_371, ik_372, \
                         ik_373, lk_837, lk_838, lk_839, lk_840, \
                         lk_841 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -ik_369[k]
                   + f_0 * lk_837[k];

        t_586[k] = -ik_370[k]
                   + f_0 * lk_838[k];

        t_587[k] = -ik_371[k]
                   + f_0 * lk_839[k];

        t_588[k] = -ik_372[k]
                   + f_0 * lk_840[k];

        t_589[k] = -ik_373[k]
                   + f_0 * lk_841[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, ik_374, ik_375, ik_376, ik_377, \
                         ik_378, lk_842, lk_843, lk_844, lk_845, \
                         lk_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -ik_374[k]
                   + f_0 * lk_842[k];

        t_591[k] = -ik_375[k]
                   + f_0 * lk_843[k];

        t_592[k] = -ik_376[k]
                   + f_0 * lk_844[k];

        t_593[k] = -ik_377[k]
                   + f_0 * lk_845[k];

        t_594[k] = -ik_378[k]
                   + f_0 * lk_846[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, ik_379, ik_380, ik_381, ik_382, \
                         ik_383, lk_847, lk_848, lk_849, lk_850, \
                         lk_851 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -ik_379[k]
                   + f_0 * lk_847[k];

        t_596[k] = -ik_380[k]
                   + f_0 * lk_848[k];

        t_597[k] = -ik_381[k]
                   + f_0 * lk_849[k];

        t_598[k] = -ik_382[k]
                   + f_0 * lk_850[k];

        t_599[k] = -ik_383[k]
                   + f_0 * lk_851[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, ik_384, ik_385, ik_386, ik_387, \
                         ik_388, lk_852, lk_853, lk_854, lk_855, \
                         lk_856 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -ik_384[k]
                   + f_0 * lk_852[k];

        t_601[k] = -ik_385[k]
                   + f_0 * lk_853[k];

        t_602[k] = -ik_386[k]
                   + f_0 * lk_854[k];

        t_603[k] = -ik_387[k]
                   + f_0 * lk_855[k];

        t_604[k] = -ik_388[k]
                   + f_0 * lk_856[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, ik_389, ik_390, ik_391, ik_392, \
                         ik_393, lk_857, lk_858, lk_859, lk_860, \
                         lk_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -ik_389[k]
                   + f_0 * lk_857[k];

        t_606[k] = -ik_390[k]
                   + f_0 * lk_858[k];

        t_607[k] = -ik_391[k]
                   + f_0 * lk_859[k];

        t_608[k] = -ik_392[k]
                   + f_0 * lk_860[k];

        t_609[k] = -ik_393[k]
                   + f_0 * lk_861[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, ik_394, ik_395, ik_396, ik_397, \
                         ik_398, lk_862, lk_863, lk_864, lk_865, \
                         lk_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -ik_394[k]
                   + f_0 * lk_862[k];

        t_611[k] = -ik_395[k]
                   + f_0 * lk_863[k];

        t_612[k] = -2.0 * ik_396[k]
                   + f_0 * lk_864[k];

        t_613[k] = -2.0 * ik_397[k]
                   + f_0 * lk_865[k];

        t_614[k] = -2.0 * ik_398[k]
                   + f_0 * lk_866[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, ik_399, ik_400, ik_401, ik_402, \
                         ik_403, lk_867, lk_868, lk_869, lk_870, \
                         lk_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -2.0 * ik_399[k]
                   + f_0 * lk_867[k];

        t_616[k] = -2.0 * ik_400[k]
                   + f_0 * lk_868[k];

        t_617[k] = -2.0 * ik_401[k]
                   + f_0 * lk_869[k];

        t_618[k] = -2.0 * ik_402[k]
                   + f_0 * lk_870[k];

        t_619[k] = -2.0 * ik_403[k]
                   + f_0 * lk_871[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, ik_404, ik_405, ik_406, ik_407, \
                         ik_408, lk_872, lk_873, lk_874, lk_875, \
                         lk_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -2.0 * ik_404[k]
                   + f_0 * lk_872[k];

        t_621[k] = -2.0 * ik_405[k]
                   + f_0 * lk_873[k];

        t_622[k] = -2.0 * ik_406[k]
                   + f_0 * lk_874[k];

        t_623[k] = -2.0 * ik_407[k]
                   + f_0 * lk_875[k];

        t_624[k] = -2.0 * ik_408[k]
                   + f_0 * lk_876[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, ik_409, ik_410, ik_411, ik_412, \
                         ik_413, lk_877, lk_878, lk_879, lk_880, \
                         lk_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -2.0 * ik_409[k]
                   + f_0 * lk_877[k];

        t_626[k] = -2.0 * ik_410[k]
                   + f_0 * lk_878[k];

        t_627[k] = -2.0 * ik_411[k]
                   + f_0 * lk_879[k];

        t_628[k] = -2.0 * ik_412[k]
                   + f_0 * lk_880[k];

        t_629[k] = -2.0 * ik_413[k]
                   + f_0 * lk_881[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, ik_414, ik_415, ik_416, ik_417, \
                         ik_418, lk_882, lk_883, lk_884, lk_885, \
                         lk_886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -2.0 * ik_414[k]
                   + f_0 * lk_882[k];

        t_631[k] = -2.0 * ik_415[k]
                   + f_0 * lk_883[k];

        t_632[k] = -2.0 * ik_416[k]
                   + f_0 * lk_884[k];

        t_633[k] = -2.0 * ik_417[k]
                   + f_0 * lk_885[k];

        t_634[k] = -2.0 * ik_418[k]
                   + f_0 * lk_886[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, ik_419, ik_420, ik_421, ik_422, \
                         ik_423, lk_887, lk_888, lk_889, lk_890, \
                         lk_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -2.0 * ik_419[k]
                   + f_0 * lk_887[k];

        t_636[k] = -2.0 * ik_420[k]
                   + f_0 * lk_888[k];

        t_637[k] = -2.0 * ik_421[k]
                   + f_0 * lk_889[k];

        t_638[k] = -2.0 * ik_422[k]
                   + f_0 * lk_890[k];

        t_639[k] = -2.0 * ik_423[k]
                   + f_0 * lk_891[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, ik_424, ik_425, ik_426, ik_427, \
                         ik_428, lk_892, lk_893, lk_894, lk_895, \
                         lk_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -2.0 * ik_424[k]
                   + f_0 * lk_892[k];

        t_641[k] = -2.0 * ik_425[k]
                   + f_0 * lk_893[k];

        t_642[k] = -2.0 * ik_426[k]
                   + f_0 * lk_894[k];

        t_643[k] = -2.0 * ik_427[k]
                   + f_0 * lk_895[k];

        t_644[k] = -2.0 * ik_428[k]
                   + f_0 * lk_896[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, ik_429, ik_430, ik_431, ik_432, \
                         ik_433, lk_897, lk_898, lk_899, lk_900, \
                         lk_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -2.0 * ik_429[k]
                   + f_0 * lk_897[k];

        t_646[k] = -2.0 * ik_430[k]
                   + f_0 * lk_898[k];

        t_647[k] = -2.0 * ik_431[k]
                   + f_0 * lk_899[k];

        t_648[k] = -3.0 * ik_432[k]
                   + f_0 * lk_900[k];

        t_649[k] = -3.0 * ik_433[k]
                   + f_0 * lk_901[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, ik_434, ik_435, ik_436, ik_437, \
                         ik_438, lk_902, lk_903, lk_904, lk_905, \
                         lk_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -3.0 * ik_434[k]
                   + f_0 * lk_902[k];

        t_651[k] = -3.0 * ik_435[k]
                   + f_0 * lk_903[k];

        t_652[k] = -3.0 * ik_436[k]
                   + f_0 * lk_904[k];

        t_653[k] = -3.0 * ik_437[k]
                   + f_0 * lk_905[k];

        t_654[k] = -3.0 * ik_438[k]
                   + f_0 * lk_906[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, ik_439, ik_440, ik_441, ik_442, \
                         ik_443, lk_907, lk_908, lk_909, lk_910, \
                         lk_911 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -3.0 * ik_439[k]
                   + f_0 * lk_907[k];

        t_656[k] = -3.0 * ik_440[k]
                   + f_0 * lk_908[k];

        t_657[k] = -3.0 * ik_441[k]
                   + f_0 * lk_909[k];

        t_658[k] = -3.0 * ik_442[k]
                   + f_0 * lk_910[k];

        t_659[k] = -3.0 * ik_443[k]
                   + f_0 * lk_911[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, ik_444, ik_445, ik_446, ik_447, \
                         ik_448, lk_912, lk_913, lk_914, lk_915, \
                         lk_916 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -3.0 * ik_444[k]
                   + f_0 * lk_912[k];

        t_661[k] = -3.0 * ik_445[k]
                   + f_0 * lk_913[k];

        t_662[k] = -3.0 * ik_446[k]
                   + f_0 * lk_914[k];

        t_663[k] = -3.0 * ik_447[k]
                   + f_0 * lk_915[k];

        t_664[k] = -3.0 * ik_448[k]
                   + f_0 * lk_916[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, ik_449, ik_450, ik_451, ik_452, \
                         ik_453, lk_917, lk_918, lk_919, lk_920, \
                         lk_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -3.0 * ik_449[k]
                   + f_0 * lk_917[k];

        t_666[k] = -3.0 * ik_450[k]
                   + f_0 * lk_918[k];

        t_667[k] = -3.0 * ik_451[k]
                   + f_0 * lk_919[k];

        t_668[k] = -3.0 * ik_452[k]
                   + f_0 * lk_920[k];

        t_669[k] = -3.0 * ik_453[k]
                   + f_0 * lk_921[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ik_454 = buffer.data(ik + 454);
    const auto *ik_455 = buffer.data(ik + 455);
    const auto *ik_456 = buffer.data(ik + 456);
    const auto *ik_457 = buffer.data(ik + 457);
    const auto *ik_458 = buffer.data(ik + 458);
    const auto *ik_459 = buffer.data(ik + 459);
    const auto *ik_460 = buffer.data(ik + 460);
    const auto *ik_461 = buffer.data(ik + 461);
    const auto *ik_462 = buffer.data(ik + 462);
    const auto *ik_463 = buffer.data(ik + 463);
    const auto *ik_464 = buffer.data(ik + 464);
    const auto *ik_465 = buffer.data(ik + 465);
    const auto *ik_466 = buffer.data(ik + 466);
    const auto *ik_467 = buffer.data(ik + 467);
    const auto *ik_468 = buffer.data(ik + 468);
    const auto *ik_469 = buffer.data(ik + 469);
    const auto *ik_470 = buffer.data(ik + 470);
    const auto *ik_471 = buffer.data(ik + 471);
    const auto *ik_472 = buffer.data(ik + 472);
    const auto *ik_473 = buffer.data(ik + 473);
    const auto *ik_474 = buffer.data(ik + 474);
    const auto *ik_475 = buffer.data(ik + 475);
    const auto *ik_476 = buffer.data(ik + 476);
    const auto *ik_477 = buffer.data(ik + 477);
    const auto *ik_478 = buffer.data(ik + 478);
    const auto *ik_479 = buffer.data(ik + 479);
    const auto *ik_480 = buffer.data(ik + 480);
    const auto *ik_481 = buffer.data(ik + 481);
    const auto *ik_482 = buffer.data(ik + 482);
    const auto *ik_483 = buffer.data(ik + 483);
    const auto *ik_484 = buffer.data(ik + 484);
    const auto *ik_485 = buffer.data(ik + 485);
    const auto *ik_486 = buffer.data(ik + 486);
    const auto *ik_487 = buffer.data(ik + 487);
    const auto *ik_488 = buffer.data(ik + 488);
    const auto *ik_489 = buffer.data(ik + 489);
    const auto *ik_490 = buffer.data(ik + 490);
    const auto *ik_491 = buffer.data(ik + 491);
    const auto *ik_492 = buffer.data(ik + 492);
    const auto *ik_493 = buffer.data(ik + 493);
    const auto *ik_494 = buffer.data(ik + 494);
    const auto *ik_495 = buffer.data(ik + 495);
    const auto *ik_496 = buffer.data(ik + 496);
    const auto *ik_497 = buffer.data(ik + 497);
    const auto *ik_498 = buffer.data(ik + 498);
    const auto *ik_499 = buffer.data(ik + 499);
    const auto *ik_500 = buffer.data(ik + 500);
    const auto *ik_501 = buffer.data(ik + 501);
    const auto *ik_502 = buffer.data(ik + 502);
    const auto *ik_503 = buffer.data(ik + 503);
    const auto *ik_504 = buffer.data(ik + 504);
    const auto *ik_505 = buffer.data(ik + 505);
    const auto *ik_506 = buffer.data(ik + 506);
    const auto *ik_507 = buffer.data(ik + 507);
    const auto *ik_508 = buffer.data(ik + 508);
    const auto *ik_509 = buffer.data(ik + 509);
    const auto *ik_510 = buffer.data(ik + 510);
    const auto *ik_511 = buffer.data(ik + 511);
    const auto *ik_512 = buffer.data(ik + 512);
    const auto *ik_513 = buffer.data(ik + 513);
    const auto *ik_514 = buffer.data(ik + 514);
    const auto *ik_515 = buffer.data(ik + 515);
    const auto *ik_516 = buffer.data(ik + 516);
    const auto *ik_517 = buffer.data(ik + 517);
    const auto *ik_518 = buffer.data(ik + 518);
    const auto *ik_519 = buffer.data(ik + 519);
    const auto *ik_520 = buffer.data(ik + 520);
    const auto *ik_521 = buffer.data(ik + 521);
    const auto *ik_522 = buffer.data(ik + 522);
    const auto *ik_523 = buffer.data(ik + 523);
    const auto *ik_524 = buffer.data(ik + 524);
    const auto *ik_525 = buffer.data(ik + 525);
    const auto *ik_526 = buffer.data(ik + 526);
    const auto *ik_527 = buffer.data(ik + 527);
    const auto *ik_528 = buffer.data(ik + 528);
    const auto *ik_529 = buffer.data(ik + 529);
    const auto *ik_530 = buffer.data(ik + 530);
    const auto *ik_531 = buffer.data(ik + 531);
    const auto *ik_532 = buffer.data(ik + 532);
    const auto *ik_533 = buffer.data(ik + 533);
    const auto *ik_534 = buffer.data(ik + 534);
    const auto *ik_535 = buffer.data(ik + 535);
    const auto *ik_536 = buffer.data(ik + 536);
    const auto *ik_537 = buffer.data(ik + 537);
    const auto *ik_538 = buffer.data(ik + 538);
    const auto *ik_539 = buffer.data(ik + 539);
    const auto *ik_540 = buffer.data(ik + 540);
    const auto *ik_541 = buffer.data(ik + 541);
    const auto *ik_542 = buffer.data(ik + 542);
    const auto *ik_543 = buffer.data(ik + 543);
    const auto *ik_544 = buffer.data(ik + 544);
    const auto *ik_545 = buffer.data(ik + 545);
    const auto *ik_546 = buffer.data(ik + 546);
    const auto *ik_547 = buffer.data(ik + 547);
    const auto *ik_548 = buffer.data(ik + 548);
    const auto *ik_549 = buffer.data(ik + 549);
    const auto *ik_550 = buffer.data(ik + 550);
    const auto *ik_551 = buffer.data(ik + 551);
    const auto *ik_552 = buffer.data(ik + 552);
    const auto *ik_553 = buffer.data(ik + 553);
    const auto *ik_554 = buffer.data(ik + 554);
    const auto *ik_555 = buffer.data(ik + 555);
    const auto *ik_556 = buffer.data(ik + 556);
    const auto *ik_557 = buffer.data(ik + 557);
    const auto *ik_558 = buffer.data(ik + 558);
    const auto *ik_559 = buffer.data(ik + 559);
    const auto *ik_560 = buffer.data(ik + 560);
    const auto *ik_561 = buffer.data(ik + 561);
    const auto *ik_562 = buffer.data(ik + 562);
    const auto *ik_563 = buffer.data(ik + 563);
    const auto *ik_564 = buffer.data(ik + 564);
    const auto *ik_565 = buffer.data(ik + 565);
    const auto *ik_566 = buffer.data(ik + 566);
    const auto *ik_567 = buffer.data(ik + 567);
    const auto *ik_568 = buffer.data(ik + 568);
    const auto *ik_569 = buffer.data(ik + 569);
    const auto *ik_570 = buffer.data(ik + 570);
    const auto *ik_571 = buffer.data(ik + 571);
    const auto *ik_572 = buffer.data(ik + 572);
    const auto *ik_573 = buffer.data(ik + 573);
    const auto *ik_574 = buffer.data(ik + 574);
    const auto *ik_575 = buffer.data(ik + 575);
    const auto *ik_576 = buffer.data(ik + 576);
    const auto *ik_577 = buffer.data(ik + 577);
    const auto *ik_578 = buffer.data(ik + 578);
    const auto *ik_579 = buffer.data(ik + 579);
    const auto *ik_580 = buffer.data(ik + 580);

    const auto *lk_922 = buffer.data(lk + 922);
    const auto *lk_923 = buffer.data(lk + 923);
    const auto *lk_924 = buffer.data(lk + 924);
    const auto *lk_925 = buffer.data(lk + 925);
    const auto *lk_926 = buffer.data(lk + 926);
    const auto *lk_927 = buffer.data(lk + 927);
    const auto *lk_928 = buffer.data(lk + 928);
    const auto *lk_929 = buffer.data(lk + 929);
    const auto *lk_930 = buffer.data(lk + 930);
    const auto *lk_931 = buffer.data(lk + 931);
    const auto *lk_932 = buffer.data(lk + 932);
    const auto *lk_933 = buffer.data(lk + 933);
    const auto *lk_934 = buffer.data(lk + 934);
    const auto *lk_935 = buffer.data(lk + 935);
    const auto *lk_936 = buffer.data(lk + 936);
    const auto *lk_937 = buffer.data(lk + 937);
    const auto *lk_938 = buffer.data(lk + 938);
    const auto *lk_939 = buffer.data(lk + 939);
    const auto *lk_940 = buffer.data(lk + 940);
    const auto *lk_941 = buffer.data(lk + 941);
    const auto *lk_942 = buffer.data(lk + 942);
    const auto *lk_943 = buffer.data(lk + 943);
    const auto *lk_944 = buffer.data(lk + 944);
    const auto *lk_945 = buffer.data(lk + 945);
    const auto *lk_946 = buffer.data(lk + 946);
    const auto *lk_947 = buffer.data(lk + 947);
    const auto *lk_948 = buffer.data(lk + 948);
    const auto *lk_949 = buffer.data(lk + 949);
    const auto *lk_950 = buffer.data(lk + 950);
    const auto *lk_951 = buffer.data(lk + 951);
    const auto *lk_952 = buffer.data(lk + 952);
    const auto *lk_953 = buffer.data(lk + 953);
    const auto *lk_954 = buffer.data(lk + 954);
    const auto *lk_955 = buffer.data(lk + 955);
    const auto *lk_956 = buffer.data(lk + 956);
    const auto *lk_957 = buffer.data(lk + 957);
    const auto *lk_958 = buffer.data(lk + 958);
    const auto *lk_959 = buffer.data(lk + 959);
    const auto *lk_960 = buffer.data(lk + 960);
    const auto *lk_961 = buffer.data(lk + 961);
    const auto *lk_962 = buffer.data(lk + 962);
    const auto *lk_963 = buffer.data(lk + 963);
    const auto *lk_964 = buffer.data(lk + 964);
    const auto *lk_965 = buffer.data(lk + 965);
    const auto *lk_966 = buffer.data(lk + 966);
    const auto *lk_967 = buffer.data(lk + 967);
    const auto *lk_968 = buffer.data(lk + 968);
    const auto *lk_969 = buffer.data(lk + 969);
    const auto *lk_970 = buffer.data(lk + 970);
    const auto *lk_971 = buffer.data(lk + 971);
    const auto *lk_972 = buffer.data(lk + 972);
    const auto *lk_973 = buffer.data(lk + 973);
    const auto *lk_974 = buffer.data(lk + 974);
    const auto *lk_975 = buffer.data(lk + 975);
    const auto *lk_976 = buffer.data(lk + 976);
    const auto *lk_977 = buffer.data(lk + 977);
    const auto *lk_978 = buffer.data(lk + 978);
    const auto *lk_979 = buffer.data(lk + 979);
    const auto *lk_980 = buffer.data(lk + 980);
    const auto *lk_981 = buffer.data(lk + 981);
    const auto *lk_982 = buffer.data(lk + 982);
    const auto *lk_983 = buffer.data(lk + 983);
    const auto *lk_984 = buffer.data(lk + 984);
    const auto *lk_985 = buffer.data(lk + 985);
    const auto *lk_986 = buffer.data(lk + 986);
    const auto *lk_987 = buffer.data(lk + 987);
    const auto *lk_988 = buffer.data(lk + 988);
    const auto *lk_989 = buffer.data(lk + 989);
    const auto *lk_990 = buffer.data(lk + 990);
    const auto *lk_991 = buffer.data(lk + 991);
    const auto *lk_992 = buffer.data(lk + 992);
    const auto *lk_993 = buffer.data(lk + 993);
    const auto *lk_994 = buffer.data(lk + 994);
    const auto *lk_995 = buffer.data(lk + 995);
    const auto *lk_996 = buffer.data(lk + 996);
    const auto *lk_997 = buffer.data(lk + 997);
    const auto *lk_998 = buffer.data(lk + 998);
    const auto *lk_999 = buffer.data(lk + 999);
    const auto *lk_1000 = buffer.data(lk + 1000);
    const auto *lk_1001 = buffer.data(lk + 1001);
    const auto *lk_1002 = buffer.data(lk + 1002);
    const auto *lk_1003 = buffer.data(lk + 1003);
    const auto *lk_1004 = buffer.data(lk + 1004);
    const auto *lk_1005 = buffer.data(lk + 1005);
    const auto *lk_1006 = buffer.data(lk + 1006);
    const auto *lk_1007 = buffer.data(lk + 1007);
    const auto *lk_1044 = buffer.data(lk + 1044);
    const auto *lk_1045 = buffer.data(lk + 1045);
    const auto *lk_1046 = buffer.data(lk + 1046);
    const auto *lk_1047 = buffer.data(lk + 1047);
    const auto *lk_1048 = buffer.data(lk + 1048);
    const auto *lk_1049 = buffer.data(lk + 1049);
    const auto *lk_1050 = buffer.data(lk + 1050);
    const auto *lk_1051 = buffer.data(lk + 1051);
    const auto *lk_1052 = buffer.data(lk + 1052);
    const auto *lk_1053 = buffer.data(lk + 1053);
    const auto *lk_1054 = buffer.data(lk + 1054);
    const auto *lk_1055 = buffer.data(lk + 1055);
    const auto *lk_1056 = buffer.data(lk + 1056);
    const auto *lk_1057 = buffer.data(lk + 1057);
    const auto *lk_1058 = buffer.data(lk + 1058);
    const auto *lk_1059 = buffer.data(lk + 1059);
    const auto *lk_1060 = buffer.data(lk + 1060);
    const auto *lk_1061 = buffer.data(lk + 1061);
    const auto *lk_1062 = buffer.data(lk + 1062);
    const auto *lk_1063 = buffer.data(lk + 1063);
    const auto *lk_1064 = buffer.data(lk + 1064);
    const auto *lk_1065 = buffer.data(lk + 1065);
    const auto *lk_1066 = buffer.data(lk + 1066);
    const auto *lk_1067 = buffer.data(lk + 1067);
    const auto *lk_1068 = buffer.data(lk + 1068);
    const auto *lk_1069 = buffer.data(lk + 1069);
    const auto *lk_1070 = buffer.data(lk + 1070);
    const auto *lk_1071 = buffer.data(lk + 1071);
    const auto *lk_1072 = buffer.data(lk + 1072);
    const auto *lk_1073 = buffer.data(lk + 1073);
    const auto *lk_1074 = buffer.data(lk + 1074);
    const auto *lk_1075 = buffer.data(lk + 1075);
    const auto *lk_1076 = buffer.data(lk + 1076);
    const auto *lk_1077 = buffer.data(lk + 1077);
    const auto *lk_1078 = buffer.data(lk + 1078);
    const auto *lk_1079 = buffer.data(lk + 1079);
    const auto *lk_1080 = buffer.data(lk + 1080);
    const auto *lk_1081 = buffer.data(lk + 1081);
    const auto *lk_1082 = buffer.data(lk + 1082);
    const auto *lk_1083 = buffer.data(lk + 1083);
    const auto *lk_1084 = buffer.data(lk + 1084);
    const auto *lk_1085 = buffer.data(lk + 1085);
    const auto *lk_1086 = buffer.data(lk + 1086);
    const auto *lk_1087 = buffer.data(lk + 1087);
    const auto *lk_1088 = buffer.data(lk + 1088);
    const auto *lk_1089 = buffer.data(lk + 1089);
    const auto *lk_1090 = buffer.data(lk + 1090);
    const auto *lk_1091 = buffer.data(lk + 1091);
    const auto *lk_1092 = buffer.data(lk + 1092);
    const auto *lk_1093 = buffer.data(lk + 1093);
    const auto *lk_1094 = buffer.data(lk + 1094);
    const auto *lk_1095 = buffer.data(lk + 1095);
    const auto *lk_1096 = buffer.data(lk + 1096);
    const auto *lk_1097 = buffer.data(lk + 1097);
    const auto *lk_1098 = buffer.data(lk + 1098);
    const auto *lk_1099 = buffer.data(lk + 1099);
    const auto *lk_1100 = buffer.data(lk + 1100);
    const auto *lk_1101 = buffer.data(lk + 1101);
    const auto *lk_1102 = buffer.data(lk + 1102);
    const auto *lk_1103 = buffer.data(lk + 1103);
    const auto *lk_1104 = buffer.data(lk + 1104);
    const auto *lk_1105 = buffer.data(lk + 1105);
    const auto *lk_1106 = buffer.data(lk + 1106);
    const auto *lk_1107 = buffer.data(lk + 1107);
    const auto *lk_1108 = buffer.data(lk + 1108);
    const auto *lk_1109 = buffer.data(lk + 1109);
    const auto *lk_1110 = buffer.data(lk + 1110);
    const auto *lk_1111 = buffer.data(lk + 1111);
    const auto *lk_1112 = buffer.data(lk + 1112);
    const auto *lk_1113 = buffer.data(lk + 1113);
    const auto *lk_1114 = buffer.data(lk + 1114);
    const auto *lk_1115 = buffer.data(lk + 1115);
    const auto *lk_1116 = buffer.data(lk + 1116);
    const auto *lk_1117 = buffer.data(lk + 1117);
    const auto *lk_1118 = buffer.data(lk + 1118);
    const auto *lk_1119 = buffer.data(lk + 1119);
    const auto *lk_1120 = buffer.data(lk + 1120);

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, ik_454, ik_455, ik_456, ik_457, \
                         ik_458, lk_922, lk_923, lk_924, lk_925, \
                         lk_926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -3.0 * ik_454[k]
                   + f_0 * lk_922[k];

        t_671[k] = -3.0 * ik_455[k]
                   + f_0 * lk_923[k];

        t_672[k] = -3.0 * ik_456[k]
                   + f_0 * lk_924[k];

        t_673[k] = -3.0 * ik_457[k]
                   + f_0 * lk_925[k];

        t_674[k] = -3.0 * ik_458[k]
                   + f_0 * lk_926[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, ik_459, ik_460, ik_461, ik_462, \
                         ik_463, lk_927, lk_928, lk_929, lk_930, \
                         lk_931 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -3.0 * ik_459[k]
                   + f_0 * lk_927[k];

        t_676[k] = -3.0 * ik_460[k]
                   + f_0 * lk_928[k];

        t_677[k] = -3.0 * ik_461[k]
                   + f_0 * lk_929[k];

        t_678[k] = -3.0 * ik_462[k]
                   + f_0 * lk_930[k];

        t_679[k] = -3.0 * ik_463[k]
                   + f_0 * lk_931[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, ik_464, ik_465, ik_466, ik_467, \
                         ik_468, lk_932, lk_933, lk_934, lk_935, \
                         lk_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -3.0 * ik_464[k]
                   + f_0 * lk_932[k];

        t_681[k] = -3.0 * ik_465[k]
                   + f_0 * lk_933[k];

        t_682[k] = -3.0 * ik_466[k]
                   + f_0 * lk_934[k];

        t_683[k] = -3.0 * ik_467[k]
                   + f_0 * lk_935[k];

        t_684[k] = -4.0 * ik_468[k]
                   + f_0 * lk_936[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, ik_469, ik_470, ik_471, ik_472, \
                         ik_473, lk_937, lk_938, lk_939, lk_940, \
                         lk_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -4.0 * ik_469[k]
                   + f_0 * lk_937[k];

        t_686[k] = -4.0 * ik_470[k]
                   + f_0 * lk_938[k];

        t_687[k] = -4.0 * ik_471[k]
                   + f_0 * lk_939[k];

        t_688[k] = -4.0 * ik_472[k]
                   + f_0 * lk_940[k];

        t_689[k] = -4.0 * ik_473[k]
                   + f_0 * lk_941[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, ik_474, ik_475, ik_476, ik_477, \
                         ik_478, lk_942, lk_943, lk_944, lk_945, \
                         lk_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -4.0 * ik_474[k]
                   + f_0 * lk_942[k];

        t_691[k] = -4.0 * ik_475[k]
                   + f_0 * lk_943[k];

        t_692[k] = -4.0 * ik_476[k]
                   + f_0 * lk_944[k];

        t_693[k] = -4.0 * ik_477[k]
                   + f_0 * lk_945[k];

        t_694[k] = -4.0 * ik_478[k]
                   + f_0 * lk_946[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, ik_479, ik_480, ik_481, ik_482, \
                         ik_483, lk_947, lk_948, lk_949, lk_950, \
                         lk_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -4.0 * ik_479[k]
                   + f_0 * lk_947[k];

        t_696[k] = -4.0 * ik_480[k]
                   + f_0 * lk_948[k];

        t_697[k] = -4.0 * ik_481[k]
                   + f_0 * lk_949[k];

        t_698[k] = -4.0 * ik_482[k]
                   + f_0 * lk_950[k];

        t_699[k] = -4.0 * ik_483[k]
                   + f_0 * lk_951[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, ik_484, ik_485, ik_486, ik_487, \
                         ik_488, lk_952, lk_953, lk_954, lk_955, \
                         lk_956 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -4.0 * ik_484[k]
                   + f_0 * lk_952[k];

        t_701[k] = -4.0 * ik_485[k]
                   + f_0 * lk_953[k];

        t_702[k] = -4.0 * ik_486[k]
                   + f_0 * lk_954[k];

        t_703[k] = -4.0 * ik_487[k]
                   + f_0 * lk_955[k];

        t_704[k] = -4.0 * ik_488[k]
                   + f_0 * lk_956[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, ik_489, ik_490, ik_491, ik_492, \
                         ik_493, lk_957, lk_958, lk_959, lk_960, \
                         lk_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -4.0 * ik_489[k]
                   + f_0 * lk_957[k];

        t_706[k] = -4.0 * ik_490[k]
                   + f_0 * lk_958[k];

        t_707[k] = -4.0 * ik_491[k]
                   + f_0 * lk_959[k];

        t_708[k] = -4.0 * ik_492[k]
                   + f_0 * lk_960[k];

        t_709[k] = -4.0 * ik_493[k]
                   + f_0 * lk_961[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, ik_494, ik_495, ik_496, ik_497, \
                         ik_498, lk_962, lk_963, lk_964, lk_965, \
                         lk_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -4.0 * ik_494[k]
                   + f_0 * lk_962[k];

        t_711[k] = -4.0 * ik_495[k]
                   + f_0 * lk_963[k];

        t_712[k] = -4.0 * ik_496[k]
                   + f_0 * lk_964[k];

        t_713[k] = -4.0 * ik_497[k]
                   + f_0 * lk_965[k];

        t_714[k] = -4.0 * ik_498[k]
                   + f_0 * lk_966[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, ik_499, ik_500, ik_501, ik_502, \
                         ik_503, lk_967, lk_968, lk_969, lk_970, \
                         lk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -4.0 * ik_499[k]
                   + f_0 * lk_967[k];

        t_716[k] = -4.0 * ik_500[k]
                   + f_0 * lk_968[k];

        t_717[k] = -4.0 * ik_501[k]
                   + f_0 * lk_969[k];

        t_718[k] = -4.0 * ik_502[k]
                   + f_0 * lk_970[k];

        t_719[k] = -4.0 * ik_503[k]
                   + f_0 * lk_971[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, ik_504, ik_505, ik_506, ik_507, \
                         ik_508, lk_972, lk_973, lk_974, lk_975, \
                         lk_976 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -5.0 * ik_504[k]
                   + f_0 * lk_972[k];

        t_721[k] = -5.0 * ik_505[k]
                   + f_0 * lk_973[k];

        t_722[k] = -5.0 * ik_506[k]
                   + f_0 * lk_974[k];

        t_723[k] = -5.0 * ik_507[k]
                   + f_0 * lk_975[k];

        t_724[k] = -5.0 * ik_508[k]
                   + f_0 * lk_976[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, ik_509, ik_510, ik_511, ik_512, \
                         ik_513, lk_977, lk_978, lk_979, lk_980, \
                         lk_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -5.0 * ik_509[k]
                   + f_0 * lk_977[k];

        t_726[k] = -5.0 * ik_510[k]
                   + f_0 * lk_978[k];

        t_727[k] = -5.0 * ik_511[k]
                   + f_0 * lk_979[k];

        t_728[k] = -5.0 * ik_512[k]
                   + f_0 * lk_980[k];

        t_729[k] = -5.0 * ik_513[k]
                   + f_0 * lk_981[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, ik_514, ik_515, ik_516, ik_517, \
                         ik_518, lk_982, lk_983, lk_984, lk_985, \
                         lk_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -5.0 * ik_514[k]
                   + f_0 * lk_982[k];

        t_731[k] = -5.0 * ik_515[k]
                   + f_0 * lk_983[k];

        t_732[k] = -5.0 * ik_516[k]
                   + f_0 * lk_984[k];

        t_733[k] = -5.0 * ik_517[k]
                   + f_0 * lk_985[k];

        t_734[k] = -5.0 * ik_518[k]
                   + f_0 * lk_986[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, ik_519, ik_520, ik_521, ik_522, \
                         ik_523, lk_987, lk_988, lk_989, lk_990, \
                         lk_991 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -5.0 * ik_519[k]
                   + f_0 * lk_987[k];

        t_736[k] = -5.0 * ik_520[k]
                   + f_0 * lk_988[k];

        t_737[k] = -5.0 * ik_521[k]
                   + f_0 * lk_989[k];

        t_738[k] = -5.0 * ik_522[k]
                   + f_0 * lk_990[k];

        t_739[k] = -5.0 * ik_523[k]
                   + f_0 * lk_991[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, ik_524, ik_525, ik_526, ik_527, \
                         ik_528, lk_992, lk_993, lk_994, lk_995, \
                         lk_996 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -5.0 * ik_524[k]
                   + f_0 * lk_992[k];

        t_741[k] = -5.0 * ik_525[k]
                   + f_0 * lk_993[k];

        t_742[k] = -5.0 * ik_526[k]
                   + f_0 * lk_994[k];

        t_743[k] = -5.0 * ik_527[k]
                   + f_0 * lk_995[k];

        t_744[k] = -5.0 * ik_528[k]
                   + f_0 * lk_996[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, ik_529, ik_530, ik_531, ik_532, \
                         ik_533, lk_997, lk_998, lk_999, lk_1000, \
                         lk_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -5.0 * ik_529[k]
                   + f_0 * lk_997[k];

        t_746[k] = -5.0 * ik_530[k]
                   + f_0 * lk_998[k];

        t_747[k] = -5.0 * ik_531[k]
                   + f_0 * lk_999[k];

        t_748[k] = -5.0 * ik_532[k]
                   + f_0 * lk_1000[k];

        t_749[k] = -5.0 * ik_533[k]
                   + f_0 * lk_1001[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, ik_534, ik_535, ik_536, ik_537, \
                         ik_538, lk_1002, lk_1003, lk_1004, lk_1005, \
                         lk_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -5.0 * ik_534[k]
                   + f_0 * lk_1002[k];

        t_751[k] = -5.0 * ik_535[k]
                   + f_0 * lk_1003[k];

        t_752[k] = -5.0 * ik_536[k]
                   + f_0 * lk_1004[k];

        t_753[k] = -5.0 * ik_537[k]
                   + f_0 * lk_1005[k];

        t_754[k] = -5.0 * ik_538[k]
                   + f_0 * lk_1006[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, t_760, t_761, ik_539, lk_1007, \
                         lk_1044, lk_1045, lk_1046, lk_1047, lk_1048, \
                         lk_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -5.0 * ik_539[k]
                   + f_0 * lk_1007[k];

        t_756[k] = f_0 * lk_1044[k];

        t_757[k] = f_0 * lk_1045[k];

        t_758[k] = f_0 * lk_1046[k];

        t_759[k] = f_0 * lk_1047[k];

        t_760[k] = f_0 * lk_1048[k];

        t_761[k] = f_0 * lk_1049[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, t_766, t_767, t_768, t_769, lk_1050, \
                         lk_1051, lk_1052, lk_1053, lk_1054, lk_1055, lk_1056, \
                         lk_1057 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_0 * lk_1050[k];

        t_763[k] = f_0 * lk_1051[k];

        t_764[k] = f_0 * lk_1052[k];

        t_765[k] = f_0 * lk_1053[k];

        t_766[k] = f_0 * lk_1054[k];

        t_767[k] = f_0 * lk_1055[k];

        t_768[k] = f_0 * lk_1056[k];

        t_769[k] = f_0 * lk_1057[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, t_775, t_776, t_777, lk_1058, \
                         lk_1059, lk_1060, lk_1061, lk_1062, lk_1063, lk_1064, \
                         lk_1065 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_0 * lk_1058[k];

        t_771[k] = f_0 * lk_1059[k];

        t_772[k] = f_0 * lk_1060[k];

        t_773[k] = f_0 * lk_1061[k];

        t_774[k] = f_0 * lk_1062[k];

        t_775[k] = f_0 * lk_1063[k];

        t_776[k] = f_0 * lk_1064[k];

        t_777[k] = f_0 * lk_1065[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, t_783, t_784, t_785, lk_1066, \
                         lk_1067, lk_1068, lk_1069, lk_1070, lk_1071, lk_1072, \
                         lk_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_0 * lk_1066[k];

        t_779[k] = f_0 * lk_1067[k];

        t_780[k] = f_0 * lk_1068[k];

        t_781[k] = f_0 * lk_1069[k];

        t_782[k] = f_0 * lk_1070[k];

        t_783[k] = f_0 * lk_1071[k];

        t_784[k] = f_0 * lk_1072[k];

        t_785[k] = f_0 * lk_1073[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, t_790, t_791, t_792, ik_540, lk_1074, \
                         lk_1075, lk_1076, lk_1077, lk_1078, lk_1079, \
                         lk_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_0 * lk_1074[k];

        t_787[k] = f_0 * lk_1075[k];

        t_788[k] = f_0 * lk_1076[k];

        t_789[k] = f_0 * lk_1077[k];

        t_790[k] = f_0 * lk_1078[k];

        t_791[k] = f_0 * lk_1079[k];

        t_792[k] = -ik_540[k]
                   + f_0 * lk_1080[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, t_796, t_797, ik_541, ik_542, ik_543, ik_544, \
                         ik_545, lk_1081, lk_1082, lk_1083, lk_1084, \
                         lk_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = -ik_541[k]
                   + f_0 * lk_1081[k];

        t_794[k] = -ik_542[k]
                   + f_0 * lk_1082[k];

        t_795[k] = -ik_543[k]
                   + f_0 * lk_1083[k];

        t_796[k] = -ik_544[k]
                   + f_0 * lk_1084[k];

        t_797[k] = -ik_545[k]
                   + f_0 * lk_1085[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, t_801, t_802, ik_546, ik_547, ik_548, ik_549, \
                         ik_550, lk_1086, lk_1087, lk_1088, lk_1089, \
                         lk_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = -ik_546[k]
                   + f_0 * lk_1086[k];

        t_799[k] = -ik_547[k]
                   + f_0 * lk_1087[k];

        t_800[k] = -ik_548[k]
                   + f_0 * lk_1088[k];

        t_801[k] = -ik_549[k]
                   + f_0 * lk_1089[k];

        t_802[k] = -ik_550[k]
                   + f_0 * lk_1090[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, t_806, t_807, ik_551, ik_552, ik_553, ik_554, \
                         ik_555, lk_1091, lk_1092, lk_1093, lk_1094, \
                         lk_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = -ik_551[k]
                   + f_0 * lk_1091[k];

        t_804[k] = -ik_552[k]
                   + f_0 * lk_1092[k];

        t_805[k] = -ik_553[k]
                   + f_0 * lk_1093[k];

        t_806[k] = -ik_554[k]
                   + f_0 * lk_1094[k];

        t_807[k] = -ik_555[k]
                   + f_0 * lk_1095[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, t_812, ik_556, ik_557, ik_558, ik_559, \
                         ik_560, lk_1096, lk_1097, lk_1098, lk_1099, \
                         lk_1100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = -ik_556[k]
                   + f_0 * lk_1096[k];

        t_809[k] = -ik_557[k]
                   + f_0 * lk_1097[k];

        t_810[k] = -ik_558[k]
                   + f_0 * lk_1098[k];

        t_811[k] = -ik_559[k]
                   + f_0 * lk_1099[k];

        t_812[k] = -ik_560[k]
                   + f_0 * lk_1100[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, t_816, t_817, ik_561, ik_562, ik_563, ik_564, \
                         ik_565, lk_1101, lk_1102, lk_1103, lk_1104, \
                         lk_1105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = -ik_561[k]
                   + f_0 * lk_1101[k];

        t_814[k] = -ik_562[k]
                   + f_0 * lk_1102[k];

        t_815[k] = -ik_563[k]
                   + f_0 * lk_1103[k];

        t_816[k] = -ik_564[k]
                   + f_0 * lk_1104[k];

        t_817[k] = -ik_565[k]
                   + f_0 * lk_1105[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, t_821, t_822, ik_566, ik_567, ik_568, ik_569, \
                         ik_570, lk_1106, lk_1107, lk_1108, lk_1109, \
                         lk_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = -ik_566[k]
                   + f_0 * lk_1106[k];

        t_819[k] = -ik_567[k]
                   + f_0 * lk_1107[k];

        t_820[k] = -ik_568[k]
                   + f_0 * lk_1108[k];

        t_821[k] = -ik_569[k]
                   + f_0 * lk_1109[k];

        t_822[k] = -ik_570[k]
                   + f_0 * lk_1110[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, t_826, t_827, ik_571, ik_572, ik_573, ik_574, \
                         ik_575, lk_1111, lk_1112, lk_1113, lk_1114, \
                         lk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = -ik_571[k]
                   + f_0 * lk_1111[k];

        t_824[k] = -ik_572[k]
                   + f_0 * lk_1112[k];

        t_825[k] = -ik_573[k]
                   + f_0 * lk_1113[k];

        t_826[k] = -ik_574[k]
                   + f_0 * lk_1114[k];

        t_827[k] = -ik_575[k]
                   + f_0 * lk_1115[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, t_831, t_832, ik_576, ik_577, ik_578, ik_579, \
                         ik_580, lk_1116, lk_1117, lk_1118, lk_1119, \
                         lk_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = -2.0 * ik_576[k]
                   + f_0 * lk_1116[k];

        t_829[k] = -2.0 * ik_577[k]
                   + f_0 * lk_1117[k];

        t_830[k] = -2.0 * ik_578[k]
                   + f_0 * lk_1118[k];

        t_831[k] = -2.0 * ik_579[k]
                   + f_0 * lk_1119[k];

        t_832[k] = -2.0 * ik_580[k]
                   + f_0 * lk_1120[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_2_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);

    const auto *ik_581 = buffer.data(ik + 581);
    const auto *ik_582 = buffer.data(ik + 582);
    const auto *ik_583 = buffer.data(ik + 583);
    const auto *ik_584 = buffer.data(ik + 584);
    const auto *ik_585 = buffer.data(ik + 585);
    const auto *ik_586 = buffer.data(ik + 586);
    const auto *ik_587 = buffer.data(ik + 587);
    const auto *ik_588 = buffer.data(ik + 588);
    const auto *ik_589 = buffer.data(ik + 589);
    const auto *ik_590 = buffer.data(ik + 590);
    const auto *ik_591 = buffer.data(ik + 591);
    const auto *ik_592 = buffer.data(ik + 592);
    const auto *ik_593 = buffer.data(ik + 593);
    const auto *ik_594 = buffer.data(ik + 594);
    const auto *ik_595 = buffer.data(ik + 595);
    const auto *ik_596 = buffer.data(ik + 596);
    const auto *ik_597 = buffer.data(ik + 597);
    const auto *ik_598 = buffer.data(ik + 598);
    const auto *ik_599 = buffer.data(ik + 599);
    const auto *ik_600 = buffer.data(ik + 600);
    const auto *ik_601 = buffer.data(ik + 601);
    const auto *ik_602 = buffer.data(ik + 602);
    const auto *ik_603 = buffer.data(ik + 603);
    const auto *ik_604 = buffer.data(ik + 604);
    const auto *ik_605 = buffer.data(ik + 605);
    const auto *ik_606 = buffer.data(ik + 606);
    const auto *ik_607 = buffer.data(ik + 607);
    const auto *ik_608 = buffer.data(ik + 608);
    const auto *ik_609 = buffer.data(ik + 609);
    const auto *ik_610 = buffer.data(ik + 610);
    const auto *ik_611 = buffer.data(ik + 611);
    const auto *ik_612 = buffer.data(ik + 612);
    const auto *ik_613 = buffer.data(ik + 613);
    const auto *ik_614 = buffer.data(ik + 614);
    const auto *ik_615 = buffer.data(ik + 615);
    const auto *ik_616 = buffer.data(ik + 616);
    const auto *ik_617 = buffer.data(ik + 617);
    const auto *ik_618 = buffer.data(ik + 618);
    const auto *ik_619 = buffer.data(ik + 619);
    const auto *ik_620 = buffer.data(ik + 620);
    const auto *ik_621 = buffer.data(ik + 621);
    const auto *ik_622 = buffer.data(ik + 622);
    const auto *ik_623 = buffer.data(ik + 623);
    const auto *ik_624 = buffer.data(ik + 624);
    const auto *ik_625 = buffer.data(ik + 625);
    const auto *ik_626 = buffer.data(ik + 626);
    const auto *ik_627 = buffer.data(ik + 627);
    const auto *ik_628 = buffer.data(ik + 628);
    const auto *ik_629 = buffer.data(ik + 629);
    const auto *ik_630 = buffer.data(ik + 630);
    const auto *ik_631 = buffer.data(ik + 631);
    const auto *ik_632 = buffer.data(ik + 632);
    const auto *ik_633 = buffer.data(ik + 633);
    const auto *ik_634 = buffer.data(ik + 634);
    const auto *ik_635 = buffer.data(ik + 635);
    const auto *ik_636 = buffer.data(ik + 636);
    const auto *ik_637 = buffer.data(ik + 637);
    const auto *ik_638 = buffer.data(ik + 638);
    const auto *ik_639 = buffer.data(ik + 639);
    const auto *ik_640 = buffer.data(ik + 640);
    const auto *ik_641 = buffer.data(ik + 641);
    const auto *ik_642 = buffer.data(ik + 642);
    const auto *ik_643 = buffer.data(ik + 643);
    const auto *ik_644 = buffer.data(ik + 644);
    const auto *ik_645 = buffer.data(ik + 645);
    const auto *ik_646 = buffer.data(ik + 646);
    const auto *ik_647 = buffer.data(ik + 647);
    const auto *ik_648 = buffer.data(ik + 648);
    const auto *ik_649 = buffer.data(ik + 649);
    const auto *ik_650 = buffer.data(ik + 650);
    const auto *ik_651 = buffer.data(ik + 651);
    const auto *ik_652 = buffer.data(ik + 652);
    const auto *ik_653 = buffer.data(ik + 653);
    const auto *ik_654 = buffer.data(ik + 654);
    const auto *ik_655 = buffer.data(ik + 655);
    const auto *ik_656 = buffer.data(ik + 656);
    const auto *ik_657 = buffer.data(ik + 657);
    const auto *ik_658 = buffer.data(ik + 658);
    const auto *ik_659 = buffer.data(ik + 659);
    const auto *ik_660 = buffer.data(ik + 660);
    const auto *ik_661 = buffer.data(ik + 661);
    const auto *ik_662 = buffer.data(ik + 662);
    const auto *ik_663 = buffer.data(ik + 663);
    const auto *ik_664 = buffer.data(ik + 664);
    const auto *ik_665 = buffer.data(ik + 665);
    const auto *ik_666 = buffer.data(ik + 666);
    const auto *ik_667 = buffer.data(ik + 667);
    const auto *ik_668 = buffer.data(ik + 668);
    const auto *ik_669 = buffer.data(ik + 669);
    const auto *ik_670 = buffer.data(ik + 670);
    const auto *ik_671 = buffer.data(ik + 671);
    const auto *ik_672 = buffer.data(ik + 672);
    const auto *ik_673 = buffer.data(ik + 673);
    const auto *ik_674 = buffer.data(ik + 674);
    const auto *ik_675 = buffer.data(ik + 675);
    const auto *ik_676 = buffer.data(ik + 676);
    const auto *ik_677 = buffer.data(ik + 677);
    const auto *ik_678 = buffer.data(ik + 678);
    const auto *ik_679 = buffer.data(ik + 679);
    const auto *ik_680 = buffer.data(ik + 680);
    const auto *ik_681 = buffer.data(ik + 681);
    const auto *ik_682 = buffer.data(ik + 682);
    const auto *ik_683 = buffer.data(ik + 683);
    const auto *ik_684 = buffer.data(ik + 684);
    const auto *ik_685 = buffer.data(ik + 685);
    const auto *ik_686 = buffer.data(ik + 686);
    const auto *ik_687 = buffer.data(ik + 687);
    const auto *ik_688 = buffer.data(ik + 688);
    const auto *ik_689 = buffer.data(ik + 689);
    const auto *ik_690 = buffer.data(ik + 690);
    const auto *ik_691 = buffer.data(ik + 691);
    const auto *ik_692 = buffer.data(ik + 692);
    const auto *ik_693 = buffer.data(ik + 693);
    const auto *ik_694 = buffer.data(ik + 694);
    const auto *ik_695 = buffer.data(ik + 695);
    const auto *ik_696 = buffer.data(ik + 696);
    const auto *ik_697 = buffer.data(ik + 697);
    const auto *ik_698 = buffer.data(ik + 698);
    const auto *ik_699 = buffer.data(ik + 699);
    const auto *ik_700 = buffer.data(ik + 700);
    const auto *ik_701 = buffer.data(ik + 701);
    const auto *ik_702 = buffer.data(ik + 702);
    const auto *ik_703 = buffer.data(ik + 703);
    const auto *ik_704 = buffer.data(ik + 704);
    const auto *ik_705 = buffer.data(ik + 705);
    const auto *ik_706 = buffer.data(ik + 706);
    const auto *ik_707 = buffer.data(ik + 707);
    const auto *ik_708 = buffer.data(ik + 708);
    const auto *ik_709 = buffer.data(ik + 709);
    const auto *ik_710 = buffer.data(ik + 710);
    const auto *ik_711 = buffer.data(ik + 711);
    const auto *ik_712 = buffer.data(ik + 712);
    const auto *ik_713 = buffer.data(ik + 713);
    const auto *ik_714 = buffer.data(ik + 714);
    const auto *ik_715 = buffer.data(ik + 715);
    const auto *ik_716 = buffer.data(ik + 716);
    const auto *ik_717 = buffer.data(ik + 717);
    const auto *ik_718 = buffer.data(ik + 718);
    const auto *ik_719 = buffer.data(ik + 719);
    const auto *ik_720 = buffer.data(ik + 720);
    const auto *ik_721 = buffer.data(ik + 721);
    const auto *ik_722 = buffer.data(ik + 722);
    const auto *ik_723 = buffer.data(ik + 723);
    const auto *ik_724 = buffer.data(ik + 724);
    const auto *ik_725 = buffer.data(ik + 725);
    const auto *ik_726 = buffer.data(ik + 726);
    const auto *ik_727 = buffer.data(ik + 727);
    const auto *ik_728 = buffer.data(ik + 728);
    const auto *ik_729 = buffer.data(ik + 729);
    const auto *ik_730 = buffer.data(ik + 730);

    const auto *lk_1121 = buffer.data(lk + 1121);
    const auto *lk_1122 = buffer.data(lk + 1122);
    const auto *lk_1123 = buffer.data(lk + 1123);
    const auto *lk_1124 = buffer.data(lk + 1124);
    const auto *lk_1125 = buffer.data(lk + 1125);
    const auto *lk_1126 = buffer.data(lk + 1126);
    const auto *lk_1127 = buffer.data(lk + 1127);
    const auto *lk_1128 = buffer.data(lk + 1128);
    const auto *lk_1129 = buffer.data(lk + 1129);
    const auto *lk_1130 = buffer.data(lk + 1130);
    const auto *lk_1131 = buffer.data(lk + 1131);
    const auto *lk_1132 = buffer.data(lk + 1132);
    const auto *lk_1133 = buffer.data(lk + 1133);
    const auto *lk_1134 = buffer.data(lk + 1134);
    const auto *lk_1135 = buffer.data(lk + 1135);
    const auto *lk_1136 = buffer.data(lk + 1136);
    const auto *lk_1137 = buffer.data(lk + 1137);
    const auto *lk_1138 = buffer.data(lk + 1138);
    const auto *lk_1139 = buffer.data(lk + 1139);
    const auto *lk_1140 = buffer.data(lk + 1140);
    const auto *lk_1141 = buffer.data(lk + 1141);
    const auto *lk_1142 = buffer.data(lk + 1142);
    const auto *lk_1143 = buffer.data(lk + 1143);
    const auto *lk_1144 = buffer.data(lk + 1144);
    const auto *lk_1145 = buffer.data(lk + 1145);
    const auto *lk_1146 = buffer.data(lk + 1146);
    const auto *lk_1147 = buffer.data(lk + 1147);
    const auto *lk_1148 = buffer.data(lk + 1148);
    const auto *lk_1149 = buffer.data(lk + 1149);
    const auto *lk_1150 = buffer.data(lk + 1150);
    const auto *lk_1151 = buffer.data(lk + 1151);
    const auto *lk_1152 = buffer.data(lk + 1152);
    const auto *lk_1153 = buffer.data(lk + 1153);
    const auto *lk_1154 = buffer.data(lk + 1154);
    const auto *lk_1155 = buffer.data(lk + 1155);
    const auto *lk_1156 = buffer.data(lk + 1156);
    const auto *lk_1157 = buffer.data(lk + 1157);
    const auto *lk_1158 = buffer.data(lk + 1158);
    const auto *lk_1159 = buffer.data(lk + 1159);
    const auto *lk_1160 = buffer.data(lk + 1160);
    const auto *lk_1161 = buffer.data(lk + 1161);
    const auto *lk_1162 = buffer.data(lk + 1162);
    const auto *lk_1163 = buffer.data(lk + 1163);
    const auto *lk_1164 = buffer.data(lk + 1164);
    const auto *lk_1165 = buffer.data(lk + 1165);
    const auto *lk_1166 = buffer.data(lk + 1166);
    const auto *lk_1167 = buffer.data(lk + 1167);
    const auto *lk_1168 = buffer.data(lk + 1168);
    const auto *lk_1169 = buffer.data(lk + 1169);
    const auto *lk_1170 = buffer.data(lk + 1170);
    const auto *lk_1171 = buffer.data(lk + 1171);
    const auto *lk_1172 = buffer.data(lk + 1172);
    const auto *lk_1173 = buffer.data(lk + 1173);
    const auto *lk_1174 = buffer.data(lk + 1174);
    const auto *lk_1175 = buffer.data(lk + 1175);
    const auto *lk_1176 = buffer.data(lk + 1176);
    const auto *lk_1177 = buffer.data(lk + 1177);
    const auto *lk_1178 = buffer.data(lk + 1178);
    const auto *lk_1179 = buffer.data(lk + 1179);
    const auto *lk_1180 = buffer.data(lk + 1180);
    const auto *lk_1181 = buffer.data(lk + 1181);
    const auto *lk_1182 = buffer.data(lk + 1182);
    const auto *lk_1183 = buffer.data(lk + 1183);
    const auto *lk_1184 = buffer.data(lk + 1184);
    const auto *lk_1185 = buffer.data(lk + 1185);
    const auto *lk_1186 = buffer.data(lk + 1186);
    const auto *lk_1187 = buffer.data(lk + 1187);
    const auto *lk_1188 = buffer.data(lk + 1188);
    const auto *lk_1189 = buffer.data(lk + 1189);
    const auto *lk_1190 = buffer.data(lk + 1190);
    const auto *lk_1191 = buffer.data(lk + 1191);
    const auto *lk_1192 = buffer.data(lk + 1192);
    const auto *lk_1193 = buffer.data(lk + 1193);
    const auto *lk_1194 = buffer.data(lk + 1194);
    const auto *lk_1195 = buffer.data(lk + 1195);
    const auto *lk_1196 = buffer.data(lk + 1196);
    const auto *lk_1197 = buffer.data(lk + 1197);
    const auto *lk_1198 = buffer.data(lk + 1198);
    const auto *lk_1199 = buffer.data(lk + 1199);
    const auto *lk_1200 = buffer.data(lk + 1200);
    const auto *lk_1201 = buffer.data(lk + 1201);
    const auto *lk_1202 = buffer.data(lk + 1202);
    const auto *lk_1203 = buffer.data(lk + 1203);
    const auto *lk_1204 = buffer.data(lk + 1204);
    const auto *lk_1205 = buffer.data(lk + 1205);
    const auto *lk_1206 = buffer.data(lk + 1206);
    const auto *lk_1207 = buffer.data(lk + 1207);
    const auto *lk_1208 = buffer.data(lk + 1208);
    const auto *lk_1209 = buffer.data(lk + 1209);
    const auto *lk_1210 = buffer.data(lk + 1210);
    const auto *lk_1211 = buffer.data(lk + 1211);
    const auto *lk_1212 = buffer.data(lk + 1212);
    const auto *lk_1213 = buffer.data(lk + 1213);
    const auto *lk_1214 = buffer.data(lk + 1214);
    const auto *lk_1215 = buffer.data(lk + 1215);
    const auto *lk_1216 = buffer.data(lk + 1216);
    const auto *lk_1217 = buffer.data(lk + 1217);
    const auto *lk_1218 = buffer.data(lk + 1218);
    const auto *lk_1219 = buffer.data(lk + 1219);
    const auto *lk_1220 = buffer.data(lk + 1220);
    const auto *lk_1221 = buffer.data(lk + 1221);
    const auto *lk_1222 = buffer.data(lk + 1222);
    const auto *lk_1223 = buffer.data(lk + 1223);
    const auto *lk_1224 = buffer.data(lk + 1224);
    const auto *lk_1225 = buffer.data(lk + 1225);
    const auto *lk_1226 = buffer.data(lk + 1226);
    const auto *lk_1227 = buffer.data(lk + 1227);
    const auto *lk_1228 = buffer.data(lk + 1228);
    const auto *lk_1229 = buffer.data(lk + 1229);
    const auto *lk_1230 = buffer.data(lk + 1230);
    const auto *lk_1231 = buffer.data(lk + 1231);
    const auto *lk_1232 = buffer.data(lk + 1232);
    const auto *lk_1233 = buffer.data(lk + 1233);
    const auto *lk_1234 = buffer.data(lk + 1234);
    const auto *lk_1235 = buffer.data(lk + 1235);
    const auto *lk_1236 = buffer.data(lk + 1236);
    const auto *lk_1237 = buffer.data(lk + 1237);
    const auto *lk_1238 = buffer.data(lk + 1238);
    const auto *lk_1239 = buffer.data(lk + 1239);
    const auto *lk_1240 = buffer.data(lk + 1240);
    const auto *lk_1241 = buffer.data(lk + 1241);
    const auto *lk_1242 = buffer.data(lk + 1242);
    const auto *lk_1243 = buffer.data(lk + 1243);
    const auto *lk_1244 = buffer.data(lk + 1244);
    const auto *lk_1245 = buffer.data(lk + 1245);
    const auto *lk_1246 = buffer.data(lk + 1246);
    const auto *lk_1247 = buffer.data(lk + 1247);
    const auto *lk_1248 = buffer.data(lk + 1248);
    const auto *lk_1249 = buffer.data(lk + 1249);
    const auto *lk_1250 = buffer.data(lk + 1250);
    const auto *lk_1251 = buffer.data(lk + 1251);
    const auto *lk_1252 = buffer.data(lk + 1252);
    const auto *lk_1253 = buffer.data(lk + 1253);
    const auto *lk_1254 = buffer.data(lk + 1254);
    const auto *lk_1255 = buffer.data(lk + 1255);
    const auto *lk_1256 = buffer.data(lk + 1256);
    const auto *lk_1257 = buffer.data(lk + 1257);
    const auto *lk_1258 = buffer.data(lk + 1258);
    const auto *lk_1259 = buffer.data(lk + 1259);
    const auto *lk_1260 = buffer.data(lk + 1260);
    const auto *lk_1261 = buffer.data(lk + 1261);
    const auto *lk_1262 = buffer.data(lk + 1262);
    const auto *lk_1263 = buffer.data(lk + 1263);
    const auto *lk_1264 = buffer.data(lk + 1264);
    const auto *lk_1265 = buffer.data(lk + 1265);
    const auto *lk_1266 = buffer.data(lk + 1266);
    const auto *lk_1267 = buffer.data(lk + 1267);
    const auto *lk_1268 = buffer.data(lk + 1268);
    const auto *lk_1269 = buffer.data(lk + 1269);
    const auto *lk_1270 = buffer.data(lk + 1270);

#pragma omp simd aligned(t_833, t_834, t_835, t_836, t_837, ik_581, ik_582, ik_583, ik_584, \
                         ik_585, lk_1121, lk_1122, lk_1123, lk_1124, \
                         lk_1125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = -2.0 * ik_581[k]
                   + f_0 * lk_1121[k];

        t_834[k] = -2.0 * ik_582[k]
                   + f_0 * lk_1122[k];

        t_835[k] = -2.0 * ik_583[k]
                   + f_0 * lk_1123[k];

        t_836[k] = -2.0 * ik_584[k]
                   + f_0 * lk_1124[k];

        t_837[k] = -2.0 * ik_585[k]
                   + f_0 * lk_1125[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, t_841, t_842, ik_586, ik_587, ik_588, ik_589, \
                         ik_590, lk_1126, lk_1127, lk_1128, lk_1129, \
                         lk_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = -2.0 * ik_586[k]
                   + f_0 * lk_1126[k];

        t_839[k] = -2.0 * ik_587[k]
                   + f_0 * lk_1127[k];

        t_840[k] = -2.0 * ik_588[k]
                   + f_0 * lk_1128[k];

        t_841[k] = -2.0 * ik_589[k]
                   + f_0 * lk_1129[k];

        t_842[k] = -2.0 * ik_590[k]
                   + f_0 * lk_1130[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, t_846, t_847, ik_591, ik_592, ik_593, ik_594, \
                         ik_595, lk_1131, lk_1132, lk_1133, lk_1134, \
                         lk_1135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = -2.0 * ik_591[k]
                   + f_0 * lk_1131[k];

        t_844[k] = -2.0 * ik_592[k]
                   + f_0 * lk_1132[k];

        t_845[k] = -2.0 * ik_593[k]
                   + f_0 * lk_1133[k];

        t_846[k] = -2.0 * ik_594[k]
                   + f_0 * lk_1134[k];

        t_847[k] = -2.0 * ik_595[k]
                   + f_0 * lk_1135[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, t_852, ik_596, ik_597, ik_598, ik_599, \
                         ik_600, lk_1136, lk_1137, lk_1138, lk_1139, \
                         lk_1140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = -2.0 * ik_596[k]
                   + f_0 * lk_1136[k];

        t_849[k] = -2.0 * ik_597[k]
                   + f_0 * lk_1137[k];

        t_850[k] = -2.0 * ik_598[k]
                   + f_0 * lk_1138[k];

        t_851[k] = -2.0 * ik_599[k]
                   + f_0 * lk_1139[k];

        t_852[k] = -2.0 * ik_600[k]
                   + f_0 * lk_1140[k];
    }

#pragma omp simd aligned(t_853, t_854, t_855, t_856, t_857, ik_601, ik_602, ik_603, ik_604, \
                         ik_605, lk_1141, lk_1142, lk_1143, lk_1144, \
                         lk_1145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_853[k] = -2.0 * ik_601[k]
                   + f_0 * lk_1141[k];

        t_854[k] = -2.0 * ik_602[k]
                   + f_0 * lk_1142[k];

        t_855[k] = -2.0 * ik_603[k]
                   + f_0 * lk_1143[k];

        t_856[k] = -2.0 * ik_604[k]
                   + f_0 * lk_1144[k];

        t_857[k] = -2.0 * ik_605[k]
                   + f_0 * lk_1145[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, t_862, ik_606, ik_607, ik_608, ik_609, \
                         ik_610, lk_1146, lk_1147, lk_1148, lk_1149, \
                         lk_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = -2.0 * ik_606[k]
                   + f_0 * lk_1146[k];

        t_859[k] = -2.0 * ik_607[k]
                   + f_0 * lk_1147[k];

        t_860[k] = -2.0 * ik_608[k]
                   + f_0 * lk_1148[k];

        t_861[k] = -2.0 * ik_609[k]
                   + f_0 * lk_1149[k];

        t_862[k] = -2.0 * ik_610[k]
                   + f_0 * lk_1150[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, t_866, t_867, ik_611, ik_612, ik_613, ik_614, \
                         ik_615, lk_1151, lk_1152, lk_1153, lk_1154, \
                         lk_1155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = -2.0 * ik_611[k]
                   + f_0 * lk_1151[k];

        t_864[k] = -3.0 * ik_612[k]
                   + f_0 * lk_1152[k];

        t_865[k] = -3.0 * ik_613[k]
                   + f_0 * lk_1153[k];

        t_866[k] = -3.0 * ik_614[k]
                   + f_0 * lk_1154[k];

        t_867[k] = -3.0 * ik_615[k]
                   + f_0 * lk_1155[k];
    }

#pragma omp simd aligned(t_868, t_869, t_870, t_871, t_872, ik_616, ik_617, ik_618, ik_619, \
                         ik_620, lk_1156, lk_1157, lk_1158, lk_1159, \
                         lk_1160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_868[k] = -3.0 * ik_616[k]
                   + f_0 * lk_1156[k];

        t_869[k] = -3.0 * ik_617[k]
                   + f_0 * lk_1157[k];

        t_870[k] = -3.0 * ik_618[k]
                   + f_0 * lk_1158[k];

        t_871[k] = -3.0 * ik_619[k]
                   + f_0 * lk_1159[k];

        t_872[k] = -3.0 * ik_620[k]
                   + f_0 * lk_1160[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, t_877, ik_621, ik_622, ik_623, ik_624, \
                         ik_625, lk_1161, lk_1162, lk_1163, lk_1164, \
                         lk_1165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = -3.0 * ik_621[k]
                   + f_0 * lk_1161[k];

        t_874[k] = -3.0 * ik_622[k]
                   + f_0 * lk_1162[k];

        t_875[k] = -3.0 * ik_623[k]
                   + f_0 * lk_1163[k];

        t_876[k] = -3.0 * ik_624[k]
                   + f_0 * lk_1164[k];

        t_877[k] = -3.0 * ik_625[k]
                   + f_0 * lk_1165[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, t_881, t_882, ik_626, ik_627, ik_628, ik_629, \
                         ik_630, lk_1166, lk_1167, lk_1168, lk_1169, \
                         lk_1170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = -3.0 * ik_626[k]
                   + f_0 * lk_1166[k];

        t_879[k] = -3.0 * ik_627[k]
                   + f_0 * lk_1167[k];

        t_880[k] = -3.0 * ik_628[k]
                   + f_0 * lk_1168[k];

        t_881[k] = -3.0 * ik_629[k]
                   + f_0 * lk_1169[k];

        t_882[k] = -3.0 * ik_630[k]
                   + f_0 * lk_1170[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, t_887, ik_631, ik_632, ik_633, ik_634, \
                         ik_635, lk_1171, lk_1172, lk_1173, lk_1174, \
                         lk_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = -3.0 * ik_631[k]
                   + f_0 * lk_1171[k];

        t_884[k] = -3.0 * ik_632[k]
                   + f_0 * lk_1172[k];

        t_885[k] = -3.0 * ik_633[k]
                   + f_0 * lk_1173[k];

        t_886[k] = -3.0 * ik_634[k]
                   + f_0 * lk_1174[k];

        t_887[k] = -3.0 * ik_635[k]
                   + f_0 * lk_1175[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, t_892, ik_636, ik_637, ik_638, ik_639, \
                         ik_640, lk_1176, lk_1177, lk_1178, lk_1179, \
                         lk_1180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = -3.0 * ik_636[k]
                   + f_0 * lk_1176[k];

        t_889[k] = -3.0 * ik_637[k]
                   + f_0 * lk_1177[k];

        t_890[k] = -3.0 * ik_638[k]
                   + f_0 * lk_1178[k];

        t_891[k] = -3.0 * ik_639[k]
                   + f_0 * lk_1179[k];

        t_892[k] = -3.0 * ik_640[k]
                   + f_0 * lk_1180[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, t_896, t_897, ik_641, ik_642, ik_643, ik_644, \
                         ik_645, lk_1181, lk_1182, lk_1183, lk_1184, \
                         lk_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = -3.0 * ik_641[k]
                   + f_0 * lk_1181[k];

        t_894[k] = -3.0 * ik_642[k]
                   + f_0 * lk_1182[k];

        t_895[k] = -3.0 * ik_643[k]
                   + f_0 * lk_1183[k];

        t_896[k] = -3.0 * ik_644[k]
                   + f_0 * lk_1184[k];

        t_897[k] = -3.0 * ik_645[k]
                   + f_0 * lk_1185[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, t_902, ik_646, ik_647, ik_648, ik_649, \
                         ik_650, lk_1186, lk_1187, lk_1188, lk_1189, \
                         lk_1190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = -3.0 * ik_646[k]
                   + f_0 * lk_1186[k];

        t_899[k] = -3.0 * ik_647[k]
                   + f_0 * lk_1187[k];

        t_900[k] = -4.0 * ik_648[k]
                   + f_0 * lk_1188[k];

        t_901[k] = -4.0 * ik_649[k]
                   + f_0 * lk_1189[k];

        t_902[k] = -4.0 * ik_650[k]
                   + f_0 * lk_1190[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, t_906, t_907, ik_651, ik_652, ik_653, ik_654, \
                         ik_655, lk_1191, lk_1192, lk_1193, lk_1194, \
                         lk_1195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = -4.0 * ik_651[k]
                   + f_0 * lk_1191[k];

        t_904[k] = -4.0 * ik_652[k]
                   + f_0 * lk_1192[k];

        t_905[k] = -4.0 * ik_653[k]
                   + f_0 * lk_1193[k];

        t_906[k] = -4.0 * ik_654[k]
                   + f_0 * lk_1194[k];

        t_907[k] = -4.0 * ik_655[k]
                   + f_0 * lk_1195[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, t_911, t_912, ik_656, ik_657, ik_658, ik_659, \
                         ik_660, lk_1196, lk_1197, lk_1198, lk_1199, \
                         lk_1200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = -4.0 * ik_656[k]
                   + f_0 * lk_1196[k];

        t_909[k] = -4.0 * ik_657[k]
                   + f_0 * lk_1197[k];

        t_910[k] = -4.0 * ik_658[k]
                   + f_0 * lk_1198[k];

        t_911[k] = -4.0 * ik_659[k]
                   + f_0 * lk_1199[k];

        t_912[k] = -4.0 * ik_660[k]
                   + f_0 * lk_1200[k];
    }

#pragma omp simd aligned(t_913, t_914, t_915, t_916, t_917, ik_661, ik_662, ik_663, ik_664, \
                         ik_665, lk_1201, lk_1202, lk_1203, lk_1204, \
                         lk_1205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_913[k] = -4.0 * ik_661[k]
                   + f_0 * lk_1201[k];

        t_914[k] = -4.0 * ik_662[k]
                   + f_0 * lk_1202[k];

        t_915[k] = -4.0 * ik_663[k]
                   + f_0 * lk_1203[k];

        t_916[k] = -4.0 * ik_664[k]
                   + f_0 * lk_1204[k];

        t_917[k] = -4.0 * ik_665[k]
                   + f_0 * lk_1205[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, t_921, t_922, ik_666, ik_667, ik_668, ik_669, \
                         ik_670, lk_1206, lk_1207, lk_1208, lk_1209, \
                         lk_1210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = -4.0 * ik_666[k]
                   + f_0 * lk_1206[k];

        t_919[k] = -4.0 * ik_667[k]
                   + f_0 * lk_1207[k];

        t_920[k] = -4.0 * ik_668[k]
                   + f_0 * lk_1208[k];

        t_921[k] = -4.0 * ik_669[k]
                   + f_0 * lk_1209[k];

        t_922[k] = -4.0 * ik_670[k]
                   + f_0 * lk_1210[k];
    }

#pragma omp simd aligned(t_923, t_924, t_925, t_926, t_927, ik_671, ik_672, ik_673, ik_674, \
                         ik_675, lk_1211, lk_1212, lk_1213, lk_1214, \
                         lk_1215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_923[k] = -4.0 * ik_671[k]
                   + f_0 * lk_1211[k];

        t_924[k] = -4.0 * ik_672[k]
                   + f_0 * lk_1212[k];

        t_925[k] = -4.0 * ik_673[k]
                   + f_0 * lk_1213[k];

        t_926[k] = -4.0 * ik_674[k]
                   + f_0 * lk_1214[k];

        t_927[k] = -4.0 * ik_675[k]
                   + f_0 * lk_1215[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, t_931, t_932, ik_676, ik_677, ik_678, ik_679, \
                         ik_680, lk_1216, lk_1217, lk_1218, lk_1219, \
                         lk_1220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = -4.0 * ik_676[k]
                   + f_0 * lk_1216[k];

        t_929[k] = -4.0 * ik_677[k]
                   + f_0 * lk_1217[k];

        t_930[k] = -4.0 * ik_678[k]
                   + f_0 * lk_1218[k];

        t_931[k] = -4.0 * ik_679[k]
                   + f_0 * lk_1219[k];

        t_932[k] = -4.0 * ik_680[k]
                   + f_0 * lk_1220[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, t_937, ik_681, ik_682, ik_683, ik_684, \
                         ik_685, lk_1221, lk_1222, lk_1223, lk_1224, \
                         lk_1225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = -4.0 * ik_681[k]
                   + f_0 * lk_1221[k];

        t_934[k] = -4.0 * ik_682[k]
                   + f_0 * lk_1222[k];

        t_935[k] = -4.0 * ik_683[k]
                   + f_0 * lk_1223[k];

        t_936[k] = -5.0 * ik_684[k]
                   + f_0 * lk_1224[k];

        t_937[k] = -5.0 * ik_685[k]
                   + f_0 * lk_1225[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, t_942, ik_686, ik_687, ik_688, ik_689, \
                         ik_690, lk_1226, lk_1227, lk_1228, lk_1229, \
                         lk_1230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = -5.0 * ik_686[k]
                   + f_0 * lk_1226[k];

        t_939[k] = -5.0 * ik_687[k]
                   + f_0 * lk_1227[k];

        t_940[k] = -5.0 * ik_688[k]
                   + f_0 * lk_1228[k];

        t_941[k] = -5.0 * ik_689[k]
                   + f_0 * lk_1229[k];

        t_942[k] = -5.0 * ik_690[k]
                   + f_0 * lk_1230[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, t_946, t_947, ik_691, ik_692, ik_693, ik_694, \
                         ik_695, lk_1231, lk_1232, lk_1233, lk_1234, \
                         lk_1235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = -5.0 * ik_691[k]
                   + f_0 * lk_1231[k];

        t_944[k] = -5.0 * ik_692[k]
                   + f_0 * lk_1232[k];

        t_945[k] = -5.0 * ik_693[k]
                   + f_0 * lk_1233[k];

        t_946[k] = -5.0 * ik_694[k]
                   + f_0 * lk_1234[k];

        t_947[k] = -5.0 * ik_695[k]
                   + f_0 * lk_1235[k];
    }

#pragma omp simd aligned(t_948, t_949, t_950, t_951, t_952, ik_696, ik_697, ik_698, ik_699, \
                         ik_700, lk_1236, lk_1237, lk_1238, lk_1239, \
                         lk_1240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = -5.0 * ik_696[k]
                   + f_0 * lk_1236[k];

        t_949[k] = -5.0 * ik_697[k]
                   + f_0 * lk_1237[k];

        t_950[k] = -5.0 * ik_698[k]
                   + f_0 * lk_1238[k];

        t_951[k] = -5.0 * ik_699[k]
                   + f_0 * lk_1239[k];

        t_952[k] = -5.0 * ik_700[k]
                   + f_0 * lk_1240[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, t_957, ik_701, ik_702, ik_703, ik_704, \
                         ik_705, lk_1241, lk_1242, lk_1243, lk_1244, \
                         lk_1245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = -5.0 * ik_701[k]
                   + f_0 * lk_1241[k];

        t_954[k] = -5.0 * ik_702[k]
                   + f_0 * lk_1242[k];

        t_955[k] = -5.0 * ik_703[k]
                   + f_0 * lk_1243[k];

        t_956[k] = -5.0 * ik_704[k]
                   + f_0 * lk_1244[k];

        t_957[k] = -5.0 * ik_705[k]
                   + f_0 * lk_1245[k];
    }

#pragma omp simd aligned(t_958, t_959, t_960, t_961, t_962, ik_706, ik_707, ik_708, ik_709, \
                         ik_710, lk_1246, lk_1247, lk_1248, lk_1249, \
                         lk_1250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_958[k] = -5.0 * ik_706[k]
                   + f_0 * lk_1246[k];

        t_959[k] = -5.0 * ik_707[k]
                   + f_0 * lk_1247[k];

        t_960[k] = -5.0 * ik_708[k]
                   + f_0 * lk_1248[k];

        t_961[k] = -5.0 * ik_709[k]
                   + f_0 * lk_1249[k];

        t_962[k] = -5.0 * ik_710[k]
                   + f_0 * lk_1250[k];
    }

#pragma omp simd aligned(t_963, t_964, t_965, t_966, t_967, ik_711, ik_712, ik_713, ik_714, \
                         ik_715, lk_1251, lk_1252, lk_1253, lk_1254, \
                         lk_1255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_963[k] = -5.0 * ik_711[k]
                   + f_0 * lk_1251[k];

        t_964[k] = -5.0 * ik_712[k]
                   + f_0 * lk_1252[k];

        t_965[k] = -5.0 * ik_713[k]
                   + f_0 * lk_1253[k];

        t_966[k] = -5.0 * ik_714[k]
                   + f_0 * lk_1254[k];

        t_967[k] = -5.0 * ik_715[k]
                   + f_0 * lk_1255[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, t_972, ik_716, ik_717, ik_718, ik_719, \
                         ik_720, lk_1256, lk_1257, lk_1258, lk_1259, \
                         lk_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = -5.0 * ik_716[k]
                   + f_0 * lk_1256[k];

        t_969[k] = -5.0 * ik_717[k]
                   + f_0 * lk_1257[k];

        t_970[k] = -5.0 * ik_718[k]
                   + f_0 * lk_1258[k];

        t_971[k] = -5.0 * ik_719[k]
                   + f_0 * lk_1259[k];

        t_972[k] = -6.0 * ik_720[k]
                   + f_0 * lk_1260[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, t_976, t_977, ik_721, ik_722, ik_723, ik_724, \
                         ik_725, lk_1261, lk_1262, lk_1263, lk_1264, \
                         lk_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = -6.0 * ik_721[k]
                   + f_0 * lk_1261[k];

        t_974[k] = -6.0 * ik_722[k]
                   + f_0 * lk_1262[k];

        t_975[k] = -6.0 * ik_723[k]
                   + f_0 * lk_1263[k];

        t_976[k] = -6.0 * ik_724[k]
                   + f_0 * lk_1264[k];

        t_977[k] = -6.0 * ik_725[k]
                   + f_0 * lk_1265[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, t_982, ik_726, ik_727, ik_728, ik_729, \
                         ik_730, lk_1266, lk_1267, lk_1268, lk_1269, \
                         lk_1270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = -6.0 * ik_726[k]
                   + f_0 * lk_1266[k];

        t_979[k] = -6.0 * ik_727[k]
                   + f_0 * lk_1267[k];

        t_980[k] = -6.0 * ik_728[k]
                   + f_0 * lk_1268[k];

        t_981[k] = -6.0 * ik_729[k]
                   + f_0 * lk_1269[k];

        t_982[k] = -6.0 * ik_730[k]
                   + f_0 * lk_1270[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_2_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);

    const auto *ik_731 = buffer.data(ik + 731);
    const auto *ik_732 = buffer.data(ik + 732);
    const auto *ik_733 = buffer.data(ik + 733);
    const auto *ik_734 = buffer.data(ik + 734);
    const auto *ik_735 = buffer.data(ik + 735);
    const auto *ik_736 = buffer.data(ik + 736);
    const auto *ik_737 = buffer.data(ik + 737);
    const auto *ik_738 = buffer.data(ik + 738);
    const auto *ik_739 = buffer.data(ik + 739);
    const auto *ik_740 = buffer.data(ik + 740);
    const auto *ik_741 = buffer.data(ik + 741);
    const auto *ik_742 = buffer.data(ik + 742);
    const auto *ik_743 = buffer.data(ik + 743);
    const auto *ik_744 = buffer.data(ik + 744);
    const auto *ik_745 = buffer.data(ik + 745);
    const auto *ik_746 = buffer.data(ik + 746);
    const auto *ik_747 = buffer.data(ik + 747);
    const auto *ik_748 = buffer.data(ik + 748);
    const auto *ik_749 = buffer.data(ik + 749);
    const auto *ik_750 = buffer.data(ik + 750);
    const auto *ik_751 = buffer.data(ik + 751);
    const auto *ik_752 = buffer.data(ik + 752);
    const auto *ik_753 = buffer.data(ik + 753);
    const auto *ik_754 = buffer.data(ik + 754);
    const auto *ik_755 = buffer.data(ik + 755);
    const auto *ik_756 = buffer.data(ik + 756);
    const auto *ik_757 = buffer.data(ik + 757);
    const auto *ik_758 = buffer.data(ik + 758);
    const auto *ik_759 = buffer.data(ik + 759);
    const auto *ik_760 = buffer.data(ik + 760);
    const auto *ik_761 = buffer.data(ik + 761);
    const auto *ik_762 = buffer.data(ik + 762);
    const auto *ik_763 = buffer.data(ik + 763);
    const auto *ik_764 = buffer.data(ik + 764);
    const auto *ik_765 = buffer.data(ik + 765);
    const auto *ik_766 = buffer.data(ik + 766);
    const auto *ik_767 = buffer.data(ik + 767);
    const auto *ik_768 = buffer.data(ik + 768);
    const auto *ik_769 = buffer.data(ik + 769);
    const auto *ik_770 = buffer.data(ik + 770);
    const auto *ik_771 = buffer.data(ik + 771);
    const auto *ik_772 = buffer.data(ik + 772);
    const auto *ik_773 = buffer.data(ik + 773);
    const auto *ik_774 = buffer.data(ik + 774);
    const auto *ik_775 = buffer.data(ik + 775);
    const auto *ik_776 = buffer.data(ik + 776);
    const auto *ik_777 = buffer.data(ik + 777);
    const auto *ik_778 = buffer.data(ik + 778);
    const auto *ik_779 = buffer.data(ik + 779);
    const auto *ik_780 = buffer.data(ik + 780);
    const auto *ik_781 = buffer.data(ik + 781);
    const auto *ik_782 = buffer.data(ik + 782);
    const auto *ik_783 = buffer.data(ik + 783);
    const auto *ik_784 = buffer.data(ik + 784);
    const auto *ik_785 = buffer.data(ik + 785);
    const auto *ik_786 = buffer.data(ik + 786);
    const auto *ik_787 = buffer.data(ik + 787);
    const auto *ik_788 = buffer.data(ik + 788);
    const auto *ik_789 = buffer.data(ik + 789);
    const auto *ik_790 = buffer.data(ik + 790);
    const auto *ik_791 = buffer.data(ik + 791);
    const auto *ik_792 = buffer.data(ik + 792);
    const auto *ik_793 = buffer.data(ik + 793);
    const auto *ik_794 = buffer.data(ik + 794);
    const auto *ik_795 = buffer.data(ik + 795);
    const auto *ik_796 = buffer.data(ik + 796);
    const auto *ik_797 = buffer.data(ik + 797);
    const auto *ik_798 = buffer.data(ik + 798);
    const auto *ik_799 = buffer.data(ik + 799);
    const auto *ik_800 = buffer.data(ik + 800);
    const auto *ik_801 = buffer.data(ik + 801);
    const auto *ik_802 = buffer.data(ik + 802);
    const auto *ik_803 = buffer.data(ik + 803);
    const auto *ik_804 = buffer.data(ik + 804);
    const auto *ik_805 = buffer.data(ik + 805);
    const auto *ik_806 = buffer.data(ik + 806);
    const auto *ik_807 = buffer.data(ik + 807);
    const auto *ik_808 = buffer.data(ik + 808);
    const auto *ik_809 = buffer.data(ik + 809);
    const auto *ik_810 = buffer.data(ik + 810);
    const auto *ik_811 = buffer.data(ik + 811);
    const auto *ik_812 = buffer.data(ik + 812);
    const auto *ik_813 = buffer.data(ik + 813);
    const auto *ik_814 = buffer.data(ik + 814);
    const auto *ik_815 = buffer.data(ik + 815);
    const auto *ik_816 = buffer.data(ik + 816);
    const auto *ik_817 = buffer.data(ik + 817);
    const auto *ik_818 = buffer.data(ik + 818);
    const auto *ik_819 = buffer.data(ik + 819);
    const auto *ik_820 = buffer.data(ik + 820);
    const auto *ik_821 = buffer.data(ik + 821);
    const auto *ik_822 = buffer.data(ik + 822);
    const auto *ik_823 = buffer.data(ik + 823);
    const auto *ik_824 = buffer.data(ik + 824);
    const auto *ik_825 = buffer.data(ik + 825);
    const auto *ik_826 = buffer.data(ik + 826);
    const auto *ik_827 = buffer.data(ik + 827);
    const auto *ik_828 = buffer.data(ik + 828);
    const auto *ik_829 = buffer.data(ik + 829);
    const auto *ik_830 = buffer.data(ik + 830);
    const auto *ik_831 = buffer.data(ik + 831);
    const auto *ik_832 = buffer.data(ik + 832);
    const auto *ik_833 = buffer.data(ik + 833);
    const auto *ik_834 = buffer.data(ik + 834);
    const auto *ik_835 = buffer.data(ik + 835);
    const auto *ik_836 = buffer.data(ik + 836);
    const auto *ik_837 = buffer.data(ik + 837);
    const auto *ik_838 = buffer.data(ik + 838);
    const auto *ik_839 = buffer.data(ik + 839);
    const auto *ik_840 = buffer.data(ik + 840);
    const auto *ik_841 = buffer.data(ik + 841);
    const auto *ik_842 = buffer.data(ik + 842);
    const auto *ik_843 = buffer.data(ik + 843);
    const auto *ik_844 = buffer.data(ik + 844);
    const auto *ik_845 = buffer.data(ik + 845);
    const auto *ik_846 = buffer.data(ik + 846);
    const auto *ik_847 = buffer.data(ik + 847);
    const auto *ik_848 = buffer.data(ik + 848);
    const auto *ik_849 = buffer.data(ik + 849);
    const auto *ik_850 = buffer.data(ik + 850);
    const auto *ik_851 = buffer.data(ik + 851);
    const auto *ik_852 = buffer.data(ik + 852);
    const auto *ik_853 = buffer.data(ik + 853);
    const auto *ik_854 = buffer.data(ik + 854);
    const auto *ik_855 = buffer.data(ik + 855);
    const auto *ik_856 = buffer.data(ik + 856);
    const auto *ik_857 = buffer.data(ik + 857);

    const auto *lk_1271 = buffer.data(lk + 1271);
    const auto *lk_1272 = buffer.data(lk + 1272);
    const auto *lk_1273 = buffer.data(lk + 1273);
    const auto *lk_1274 = buffer.data(lk + 1274);
    const auto *lk_1275 = buffer.data(lk + 1275);
    const auto *lk_1276 = buffer.data(lk + 1276);
    const auto *lk_1277 = buffer.data(lk + 1277);
    const auto *lk_1278 = buffer.data(lk + 1278);
    const auto *lk_1279 = buffer.data(lk + 1279);
    const auto *lk_1280 = buffer.data(lk + 1280);
    const auto *lk_1281 = buffer.data(lk + 1281);
    const auto *lk_1282 = buffer.data(lk + 1282);
    const auto *lk_1283 = buffer.data(lk + 1283);
    const auto *lk_1284 = buffer.data(lk + 1284);
    const auto *lk_1285 = buffer.data(lk + 1285);
    const auto *lk_1286 = buffer.data(lk + 1286);
    const auto *lk_1287 = buffer.data(lk + 1287);
    const auto *lk_1288 = buffer.data(lk + 1288);
    const auto *lk_1289 = buffer.data(lk + 1289);
    const auto *lk_1290 = buffer.data(lk + 1290);
    const auto *lk_1291 = buffer.data(lk + 1291);
    const auto *lk_1292 = buffer.data(lk + 1292);
    const auto *lk_1293 = buffer.data(lk + 1293);
    const auto *lk_1294 = buffer.data(lk + 1294);
    const auto *lk_1295 = buffer.data(lk + 1295);
    const auto *lk_1332 = buffer.data(lk + 1332);
    const auto *lk_1333 = buffer.data(lk + 1333);
    const auto *lk_1334 = buffer.data(lk + 1334);
    const auto *lk_1335 = buffer.data(lk + 1335);
    const auto *lk_1336 = buffer.data(lk + 1336);
    const auto *lk_1337 = buffer.data(lk + 1337);
    const auto *lk_1338 = buffer.data(lk + 1338);
    const auto *lk_1339 = buffer.data(lk + 1339);
    const auto *lk_1340 = buffer.data(lk + 1340);
    const auto *lk_1341 = buffer.data(lk + 1341);
    const auto *lk_1342 = buffer.data(lk + 1342);
    const auto *lk_1343 = buffer.data(lk + 1343);
    const auto *lk_1344 = buffer.data(lk + 1344);
    const auto *lk_1345 = buffer.data(lk + 1345);
    const auto *lk_1346 = buffer.data(lk + 1346);
    const auto *lk_1347 = buffer.data(lk + 1347);
    const auto *lk_1348 = buffer.data(lk + 1348);
    const auto *lk_1349 = buffer.data(lk + 1349);
    const auto *lk_1350 = buffer.data(lk + 1350);
    const auto *lk_1351 = buffer.data(lk + 1351);
    const auto *lk_1352 = buffer.data(lk + 1352);
    const auto *lk_1353 = buffer.data(lk + 1353);
    const auto *lk_1354 = buffer.data(lk + 1354);
    const auto *lk_1355 = buffer.data(lk + 1355);
    const auto *lk_1356 = buffer.data(lk + 1356);
    const auto *lk_1357 = buffer.data(lk + 1357);
    const auto *lk_1358 = buffer.data(lk + 1358);
    const auto *lk_1359 = buffer.data(lk + 1359);
    const auto *lk_1360 = buffer.data(lk + 1360);
    const auto *lk_1361 = buffer.data(lk + 1361);
    const auto *lk_1362 = buffer.data(lk + 1362);
    const auto *lk_1363 = buffer.data(lk + 1363);
    const auto *lk_1364 = buffer.data(lk + 1364);
    const auto *lk_1365 = buffer.data(lk + 1365);
    const auto *lk_1366 = buffer.data(lk + 1366);
    const auto *lk_1367 = buffer.data(lk + 1367);
    const auto *lk_1368 = buffer.data(lk + 1368);
    const auto *lk_1369 = buffer.data(lk + 1369);
    const auto *lk_1370 = buffer.data(lk + 1370);
    const auto *lk_1371 = buffer.data(lk + 1371);
    const auto *lk_1372 = buffer.data(lk + 1372);
    const auto *lk_1373 = buffer.data(lk + 1373);
    const auto *lk_1374 = buffer.data(lk + 1374);
    const auto *lk_1375 = buffer.data(lk + 1375);
    const auto *lk_1376 = buffer.data(lk + 1376);
    const auto *lk_1377 = buffer.data(lk + 1377);
    const auto *lk_1378 = buffer.data(lk + 1378);
    const auto *lk_1379 = buffer.data(lk + 1379);
    const auto *lk_1380 = buffer.data(lk + 1380);
    const auto *lk_1381 = buffer.data(lk + 1381);
    const auto *lk_1382 = buffer.data(lk + 1382);
    const auto *lk_1383 = buffer.data(lk + 1383);
    const auto *lk_1384 = buffer.data(lk + 1384);
    const auto *lk_1385 = buffer.data(lk + 1385);
    const auto *lk_1386 = buffer.data(lk + 1386);
    const auto *lk_1387 = buffer.data(lk + 1387);
    const auto *lk_1388 = buffer.data(lk + 1388);
    const auto *lk_1389 = buffer.data(lk + 1389);
    const auto *lk_1390 = buffer.data(lk + 1390);
    const auto *lk_1391 = buffer.data(lk + 1391);
    const auto *lk_1392 = buffer.data(lk + 1392);
    const auto *lk_1393 = buffer.data(lk + 1393);
    const auto *lk_1394 = buffer.data(lk + 1394);
    const auto *lk_1395 = buffer.data(lk + 1395);
    const auto *lk_1396 = buffer.data(lk + 1396);
    const auto *lk_1397 = buffer.data(lk + 1397);
    const auto *lk_1398 = buffer.data(lk + 1398);
    const auto *lk_1399 = buffer.data(lk + 1399);
    const auto *lk_1400 = buffer.data(lk + 1400);
    const auto *lk_1401 = buffer.data(lk + 1401);
    const auto *lk_1402 = buffer.data(lk + 1402);
    const auto *lk_1403 = buffer.data(lk + 1403);
    const auto *lk_1404 = buffer.data(lk + 1404);
    const auto *lk_1405 = buffer.data(lk + 1405);
    const auto *lk_1406 = buffer.data(lk + 1406);
    const auto *lk_1407 = buffer.data(lk + 1407);
    const auto *lk_1408 = buffer.data(lk + 1408);
    const auto *lk_1409 = buffer.data(lk + 1409);
    const auto *lk_1410 = buffer.data(lk + 1410);
    const auto *lk_1411 = buffer.data(lk + 1411);
    const auto *lk_1412 = buffer.data(lk + 1412);
    const auto *lk_1413 = buffer.data(lk + 1413);
    const auto *lk_1414 = buffer.data(lk + 1414);
    const auto *lk_1415 = buffer.data(lk + 1415);
    const auto *lk_1416 = buffer.data(lk + 1416);
    const auto *lk_1417 = buffer.data(lk + 1417);
    const auto *lk_1418 = buffer.data(lk + 1418);
    const auto *lk_1419 = buffer.data(lk + 1419);
    const auto *lk_1420 = buffer.data(lk + 1420);
    const auto *lk_1421 = buffer.data(lk + 1421);
    const auto *lk_1422 = buffer.data(lk + 1422);
    const auto *lk_1423 = buffer.data(lk + 1423);
    const auto *lk_1424 = buffer.data(lk + 1424);
    const auto *lk_1425 = buffer.data(lk + 1425);
    const auto *lk_1426 = buffer.data(lk + 1426);
    const auto *lk_1427 = buffer.data(lk + 1427);
    const auto *lk_1428 = buffer.data(lk + 1428);
    const auto *lk_1429 = buffer.data(lk + 1429);
    const auto *lk_1430 = buffer.data(lk + 1430);
    const auto *lk_1431 = buffer.data(lk + 1431);
    const auto *lk_1432 = buffer.data(lk + 1432);
    const auto *lk_1433 = buffer.data(lk + 1433);
    const auto *lk_1434 = buffer.data(lk + 1434);
    const auto *lk_1435 = buffer.data(lk + 1435);
    const auto *lk_1436 = buffer.data(lk + 1436);
    const auto *lk_1437 = buffer.data(lk + 1437);
    const auto *lk_1438 = buffer.data(lk + 1438);
    const auto *lk_1439 = buffer.data(lk + 1439);
    const auto *lk_1440 = buffer.data(lk + 1440);
    const auto *lk_1441 = buffer.data(lk + 1441);
    const auto *lk_1442 = buffer.data(lk + 1442);
    const auto *lk_1443 = buffer.data(lk + 1443);
    const auto *lk_1444 = buffer.data(lk + 1444);
    const auto *lk_1445 = buffer.data(lk + 1445);
    const auto *lk_1446 = buffer.data(lk + 1446);
    const auto *lk_1447 = buffer.data(lk + 1447);
    const auto *lk_1448 = buffer.data(lk + 1448);
    const auto *lk_1449 = buffer.data(lk + 1449);
    const auto *lk_1450 = buffer.data(lk + 1450);
    const auto *lk_1451 = buffer.data(lk + 1451);
    const auto *lk_1452 = buffer.data(lk + 1452);
    const auto *lk_1453 = buffer.data(lk + 1453);
    const auto *lk_1454 = buffer.data(lk + 1454);
    const auto *lk_1455 = buffer.data(lk + 1455);
    const auto *lk_1456 = buffer.data(lk + 1456);
    const auto *lk_1457 = buffer.data(lk + 1457);
    const auto *lk_1458 = buffer.data(lk + 1458);
    const auto *lk_1459 = buffer.data(lk + 1459);
    const auto *lk_1460 = buffer.data(lk + 1460);
    const auto *lk_1461 = buffer.data(lk + 1461);
    const auto *lk_1462 = buffer.data(lk + 1462);
    const auto *lk_1463 = buffer.data(lk + 1463);
    const auto *lk_1464 = buffer.data(lk + 1464);
    const auto *lk_1465 = buffer.data(lk + 1465);
    const auto *lk_1466 = buffer.data(lk + 1466);
    const auto *lk_1467 = buffer.data(lk + 1467);
    const auto *lk_1468 = buffer.data(lk + 1468);
    const auto *lk_1469 = buffer.data(lk + 1469);

#pragma omp simd aligned(t_983, t_984, t_985, t_986, t_987, ik_731, ik_732, ik_733, ik_734, \
                         ik_735, lk_1271, lk_1272, lk_1273, lk_1274, \
                         lk_1275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = -6.0 * ik_731[k]
                   + f_0 * lk_1271[k];

        t_984[k] = -6.0 * ik_732[k]
                   + f_0 * lk_1272[k];

        t_985[k] = -6.0 * ik_733[k]
                   + f_0 * lk_1273[k];

        t_986[k] = -6.0 * ik_734[k]
                   + f_0 * lk_1274[k];

        t_987[k] = -6.0 * ik_735[k]
                   + f_0 * lk_1275[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, t_991, t_992, ik_736, ik_737, ik_738, ik_739, \
                         ik_740, lk_1276, lk_1277, lk_1278, lk_1279, \
                         lk_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = -6.0 * ik_736[k]
                   + f_0 * lk_1276[k];

        t_989[k] = -6.0 * ik_737[k]
                   + f_0 * lk_1277[k];

        t_990[k] = -6.0 * ik_738[k]
                   + f_0 * lk_1278[k];

        t_991[k] = -6.0 * ik_739[k]
                   + f_0 * lk_1279[k];

        t_992[k] = -6.0 * ik_740[k]
                   + f_0 * lk_1280[k];
    }

#pragma omp simd aligned(t_993, t_994, t_995, t_996, t_997, ik_741, ik_742, ik_743, ik_744, \
                         ik_745, lk_1281, lk_1282, lk_1283, lk_1284, \
                         lk_1285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_993[k] = -6.0 * ik_741[k]
                   + f_0 * lk_1281[k];

        t_994[k] = -6.0 * ik_742[k]
                   + f_0 * lk_1282[k];

        t_995[k] = -6.0 * ik_743[k]
                   + f_0 * lk_1283[k];

        t_996[k] = -6.0 * ik_744[k]
                   + f_0 * lk_1284[k];

        t_997[k] = -6.0 * ik_745[k]
                   + f_0 * lk_1285[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, t_1001, t_1002, ik_746, ik_747, ik_748, ik_749, \
                         ik_750, lk_1286, lk_1287, lk_1288, lk_1289, \
                         lk_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = -6.0 * ik_746[k]
                   + f_0 * lk_1286[k];

        t_999[k] = -6.0 * ik_747[k]
                   + f_0 * lk_1287[k];

        t_1000[k] = -6.0 * ik_748[k]
                    + f_0 * lk_1288[k];

        t_1001[k] = -6.0 * ik_749[k]
                    + f_0 * lk_1289[k];

        t_1002[k] = -6.0 * ik_750[k]
                    + f_0 * lk_1290[k];
    }

#pragma omp simd aligned(t_1003, t_1004, t_1005, t_1006, t_1007, ik_751, ik_752, ik_753, \
                         ik_754, ik_755, lk_1291, lk_1292, lk_1293, lk_1294, \
                         lk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1003[k] = -6.0 * ik_751[k]
                    + f_0 * lk_1291[k];

        t_1004[k] = -6.0 * ik_752[k]
                    + f_0 * lk_1292[k];

        t_1005[k] = -6.0 * ik_753[k]
                    + f_0 * lk_1293[k];

        t_1006[k] = -6.0 * ik_754[k]
                    + f_0 * lk_1294[k];

        t_1007[k] = -6.0 * ik_755[k]
                    + f_0 * lk_1295[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, t_1012, t_1013, t_1014, t_1015, \
                         lk_1332, lk_1333, lk_1334, lk_1335, lk_1336, lk_1337, lk_1338, \
                         lk_1339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_0 * lk_1332[k];

        t_1009[k] = f_0 * lk_1333[k];

        t_1010[k] = f_0 * lk_1334[k];

        t_1011[k] = f_0 * lk_1335[k];

        t_1012[k] = f_0 * lk_1336[k];

        t_1013[k] = f_0 * lk_1337[k];

        t_1014[k] = f_0 * lk_1338[k];

        t_1015[k] = f_0 * lk_1339[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, t_1020, t_1021, t_1022, t_1023, \
                         lk_1340, lk_1341, lk_1342, lk_1343, lk_1344, lk_1345, lk_1346, \
                         lk_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_0 * lk_1340[k];

        t_1017[k] = f_0 * lk_1341[k];

        t_1018[k] = f_0 * lk_1342[k];

        t_1019[k] = f_0 * lk_1343[k];

        t_1020[k] = f_0 * lk_1344[k];

        t_1021[k] = f_0 * lk_1345[k];

        t_1022[k] = f_0 * lk_1346[k];

        t_1023[k] = f_0 * lk_1347[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, t_1027, t_1028, t_1029, t_1030, t_1031, \
                         lk_1348, lk_1349, lk_1350, lk_1351, lk_1352, lk_1353, lk_1354, \
                         lk_1355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_0 * lk_1348[k];

        t_1025[k] = f_0 * lk_1349[k];

        t_1026[k] = f_0 * lk_1350[k];

        t_1027[k] = f_0 * lk_1351[k];

        t_1028[k] = f_0 * lk_1352[k];

        t_1029[k] = f_0 * lk_1353[k];

        t_1030[k] = f_0 * lk_1354[k];

        t_1031[k] = f_0 * lk_1355[k];
    }

#pragma omp simd aligned(t_1032, t_1033, t_1034, t_1035, t_1036, t_1037, t_1038, t_1039, \
                         lk_1356, lk_1357, lk_1358, lk_1359, lk_1360, lk_1361, lk_1362, \
                         lk_1363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1032[k] = f_0 * lk_1356[k];

        t_1033[k] = f_0 * lk_1357[k];

        t_1034[k] = f_0 * lk_1358[k];

        t_1035[k] = f_0 * lk_1359[k];

        t_1036[k] = f_0 * lk_1360[k];

        t_1037[k] = f_0 * lk_1361[k];

        t_1038[k] = f_0 * lk_1362[k];

        t_1039[k] = f_0 * lk_1363[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, t_1043, t_1044, t_1045, ik_756, ik_757, \
                         lk_1364, lk_1365, lk_1366, lk_1367, lk_1368, \
                         lk_1369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_0 * lk_1364[k];

        t_1041[k] = f_0 * lk_1365[k];

        t_1042[k] = f_0 * lk_1366[k];

        t_1043[k] = f_0 * lk_1367[k];

        t_1044[k] = -ik_756[k]
                    + f_0 * lk_1368[k];

        t_1045[k] = -ik_757[k]
                    + f_0 * lk_1369[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, t_1049, t_1050, ik_758, ik_759, ik_760, \
                         ik_761, ik_762, lk_1370, lk_1371, lk_1372, lk_1373, \
                         lk_1374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = -ik_758[k]
                    + f_0 * lk_1370[k];

        t_1047[k] = -ik_759[k]
                    + f_0 * lk_1371[k];

        t_1048[k] = -ik_760[k]
                    + f_0 * lk_1372[k];

        t_1049[k] = -ik_761[k]
                    + f_0 * lk_1373[k];

        t_1050[k] = -ik_762[k]
                    + f_0 * lk_1374[k];
    }

#pragma omp simd aligned(t_1051, t_1052, t_1053, t_1054, t_1055, ik_763, ik_764, ik_765, \
                         ik_766, ik_767, lk_1375, lk_1376, lk_1377, lk_1378, \
                         lk_1379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1051[k] = -ik_763[k]
                    + f_0 * lk_1375[k];

        t_1052[k] = -ik_764[k]
                    + f_0 * lk_1376[k];

        t_1053[k] = -ik_765[k]
                    + f_0 * lk_1377[k];

        t_1054[k] = -ik_766[k]
                    + f_0 * lk_1378[k];

        t_1055[k] = -ik_767[k]
                    + f_0 * lk_1379[k];
    }

#pragma omp simd aligned(t_1056, t_1057, t_1058, t_1059, t_1060, ik_768, ik_769, ik_770, \
                         ik_771, ik_772, lk_1380, lk_1381, lk_1382, lk_1383, \
                         lk_1384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1056[k] = -ik_768[k]
                    + f_0 * lk_1380[k];

        t_1057[k] = -ik_769[k]
                    + f_0 * lk_1381[k];

        t_1058[k] = -ik_770[k]
                    + f_0 * lk_1382[k];

        t_1059[k] = -ik_771[k]
                    + f_0 * lk_1383[k];

        t_1060[k] = -ik_772[k]
                    + f_0 * lk_1384[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, t_1064, t_1065, ik_773, ik_774, ik_775, \
                         ik_776, ik_777, lk_1385, lk_1386, lk_1387, lk_1388, \
                         lk_1389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = -ik_773[k]
                    + f_0 * lk_1385[k];

        t_1062[k] = -ik_774[k]
                    + f_0 * lk_1386[k];

        t_1063[k] = -ik_775[k]
                    + f_0 * lk_1387[k];

        t_1064[k] = -ik_776[k]
                    + f_0 * lk_1388[k];

        t_1065[k] = -ik_777[k]
                    + f_0 * lk_1389[k];
    }

#pragma omp simd aligned(t_1066, t_1067, t_1068, t_1069, t_1070, ik_778, ik_779, ik_780, \
                         ik_781, ik_782, lk_1390, lk_1391, lk_1392, lk_1393, \
                         lk_1394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1066[k] = -ik_778[k]
                    + f_0 * lk_1390[k];

        t_1067[k] = -ik_779[k]
                    + f_0 * lk_1391[k];

        t_1068[k] = -ik_780[k]
                    + f_0 * lk_1392[k];

        t_1069[k] = -ik_781[k]
                    + f_0 * lk_1393[k];

        t_1070[k] = -ik_782[k]
                    + f_0 * lk_1394[k];
    }

#pragma omp simd aligned(t_1071, t_1072, t_1073, t_1074, t_1075, ik_783, ik_784, ik_785, \
                         ik_786, ik_787, lk_1395, lk_1396, lk_1397, lk_1398, \
                         lk_1399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1071[k] = -ik_783[k]
                    + f_0 * lk_1395[k];

        t_1072[k] = -ik_784[k]
                    + f_0 * lk_1396[k];

        t_1073[k] = -ik_785[k]
                    + f_0 * lk_1397[k];

        t_1074[k] = -ik_786[k]
                    + f_0 * lk_1398[k];

        t_1075[k] = -ik_787[k]
                    + f_0 * lk_1399[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, t_1079, t_1080, ik_788, ik_789, ik_790, \
                         ik_791, ik_792, lk_1400, lk_1401, lk_1402, lk_1403, \
                         lk_1404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = -ik_788[k]
                    + f_0 * lk_1400[k];

        t_1077[k] = -ik_789[k]
                    + f_0 * lk_1401[k];

        t_1078[k] = -ik_790[k]
                    + f_0 * lk_1402[k];

        t_1079[k] = -ik_791[k]
                    + f_0 * lk_1403[k];

        t_1080[k] = -2.0 * ik_792[k]
                    + f_0 * lk_1404[k];
    }

#pragma omp simd aligned(t_1081, t_1082, t_1083, t_1084, t_1085, ik_793, ik_794, ik_795, \
                         ik_796, ik_797, lk_1405, lk_1406, lk_1407, lk_1408, \
                         lk_1409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1081[k] = -2.0 * ik_793[k]
                    + f_0 * lk_1405[k];

        t_1082[k] = -2.0 * ik_794[k]
                    + f_0 * lk_1406[k];

        t_1083[k] = -2.0 * ik_795[k]
                    + f_0 * lk_1407[k];

        t_1084[k] = -2.0 * ik_796[k]
                    + f_0 * lk_1408[k];

        t_1085[k] = -2.0 * ik_797[k]
                    + f_0 * lk_1409[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, t_1089, t_1090, ik_798, ik_799, ik_800, \
                         ik_801, ik_802, lk_1410, lk_1411, lk_1412, lk_1413, \
                         lk_1414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = -2.0 * ik_798[k]
                    + f_0 * lk_1410[k];

        t_1087[k] = -2.0 * ik_799[k]
                    + f_0 * lk_1411[k];

        t_1088[k] = -2.0 * ik_800[k]
                    + f_0 * lk_1412[k];

        t_1089[k] = -2.0 * ik_801[k]
                    + f_0 * lk_1413[k];

        t_1090[k] = -2.0 * ik_802[k]
                    + f_0 * lk_1414[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, t_1094, t_1095, ik_803, ik_804, ik_805, \
                         ik_806, ik_807, lk_1415, lk_1416, lk_1417, lk_1418, \
                         lk_1419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = -2.0 * ik_803[k]
                    + f_0 * lk_1415[k];

        t_1092[k] = -2.0 * ik_804[k]
                    + f_0 * lk_1416[k];

        t_1093[k] = -2.0 * ik_805[k]
                    + f_0 * lk_1417[k];

        t_1094[k] = -2.0 * ik_806[k]
                    + f_0 * lk_1418[k];

        t_1095[k] = -2.0 * ik_807[k]
                    + f_0 * lk_1419[k];
    }

#pragma omp simd aligned(t_1096, t_1097, t_1098, t_1099, t_1100, ik_808, ik_809, ik_810, \
                         ik_811, ik_812, lk_1420, lk_1421, lk_1422, lk_1423, \
                         lk_1424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1096[k] = -2.0 * ik_808[k]
                    + f_0 * lk_1420[k];

        t_1097[k] = -2.0 * ik_809[k]
                    + f_0 * lk_1421[k];

        t_1098[k] = -2.0 * ik_810[k]
                    + f_0 * lk_1422[k];

        t_1099[k] = -2.0 * ik_811[k]
                    + f_0 * lk_1423[k];

        t_1100[k] = -2.0 * ik_812[k]
                    + f_0 * lk_1424[k];
    }

#pragma omp simd aligned(t_1101, t_1102, t_1103, t_1104, t_1105, ik_813, ik_814, ik_815, \
                         ik_816, ik_817, lk_1425, lk_1426, lk_1427, lk_1428, \
                         lk_1429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1101[k] = -2.0 * ik_813[k]
                    + f_0 * lk_1425[k];

        t_1102[k] = -2.0 * ik_814[k]
                    + f_0 * lk_1426[k];

        t_1103[k] = -2.0 * ik_815[k]
                    + f_0 * lk_1427[k];

        t_1104[k] = -2.0 * ik_816[k]
                    + f_0 * lk_1428[k];

        t_1105[k] = -2.0 * ik_817[k]
                    + f_0 * lk_1429[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, t_1110, ik_818, ik_819, ik_820, \
                         ik_821, ik_822, lk_1430, lk_1431, lk_1432, lk_1433, \
                         lk_1434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = -2.0 * ik_818[k]
                    + f_0 * lk_1430[k];

        t_1107[k] = -2.0 * ik_819[k]
                    + f_0 * lk_1431[k];

        t_1108[k] = -2.0 * ik_820[k]
                    + f_0 * lk_1432[k];

        t_1109[k] = -2.0 * ik_821[k]
                    + f_0 * lk_1433[k];

        t_1110[k] = -2.0 * ik_822[k]
                    + f_0 * lk_1434[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, t_1114, t_1115, ik_823, ik_824, ik_825, \
                         ik_826, ik_827, lk_1435, lk_1436, lk_1437, lk_1438, \
                         lk_1439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = -2.0 * ik_823[k]
                    + f_0 * lk_1435[k];

        t_1112[k] = -2.0 * ik_824[k]
                    + f_0 * lk_1436[k];

        t_1113[k] = -2.0 * ik_825[k]
                    + f_0 * lk_1437[k];

        t_1114[k] = -2.0 * ik_826[k]
                    + f_0 * lk_1438[k];

        t_1115[k] = -2.0 * ik_827[k]
                    + f_0 * lk_1439[k];
    }

#pragma omp simd aligned(t_1116, t_1117, t_1118, t_1119, t_1120, ik_828, ik_829, ik_830, \
                         ik_831, ik_832, lk_1440, lk_1441, lk_1442, lk_1443, \
                         lk_1444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1116[k] = -3.0 * ik_828[k]
                    + f_0 * lk_1440[k];

        t_1117[k] = -3.0 * ik_829[k]
                    + f_0 * lk_1441[k];

        t_1118[k] = -3.0 * ik_830[k]
                    + f_0 * lk_1442[k];

        t_1119[k] = -3.0 * ik_831[k]
                    + f_0 * lk_1443[k];

        t_1120[k] = -3.0 * ik_832[k]
                    + f_0 * lk_1444[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, t_1124, t_1125, ik_833, ik_834, ik_835, \
                         ik_836, ik_837, lk_1445, lk_1446, lk_1447, lk_1448, \
                         lk_1449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = -3.0 * ik_833[k]
                    + f_0 * lk_1445[k];

        t_1122[k] = -3.0 * ik_834[k]
                    + f_0 * lk_1446[k];

        t_1123[k] = -3.0 * ik_835[k]
                    + f_0 * lk_1447[k];

        t_1124[k] = -3.0 * ik_836[k]
                    + f_0 * lk_1448[k];

        t_1125[k] = -3.0 * ik_837[k]
                    + f_0 * lk_1449[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, t_1129, t_1130, ik_838, ik_839, ik_840, \
                         ik_841, ik_842, lk_1450, lk_1451, lk_1452, lk_1453, \
                         lk_1454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = -3.0 * ik_838[k]
                    + f_0 * lk_1450[k];

        t_1127[k] = -3.0 * ik_839[k]
                    + f_0 * lk_1451[k];

        t_1128[k] = -3.0 * ik_840[k]
                    + f_0 * lk_1452[k];

        t_1129[k] = -3.0 * ik_841[k]
                    + f_0 * lk_1453[k];

        t_1130[k] = -3.0 * ik_842[k]
                    + f_0 * lk_1454[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, t_1134, t_1135, ik_843, ik_844, ik_845, \
                         ik_846, ik_847, lk_1455, lk_1456, lk_1457, lk_1458, \
                         lk_1459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = -3.0 * ik_843[k]
                    + f_0 * lk_1455[k];

        t_1132[k] = -3.0 * ik_844[k]
                    + f_0 * lk_1456[k];

        t_1133[k] = -3.0 * ik_845[k]
                    + f_0 * lk_1457[k];

        t_1134[k] = -3.0 * ik_846[k]
                    + f_0 * lk_1458[k];

        t_1135[k] = -3.0 * ik_847[k]
                    + f_0 * lk_1459[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, t_1139, t_1140, ik_848, ik_849, ik_850, \
                         ik_851, ik_852, lk_1460, lk_1461, lk_1462, lk_1463, \
                         lk_1464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = -3.0 * ik_848[k]
                    + f_0 * lk_1460[k];

        t_1137[k] = -3.0 * ik_849[k]
                    + f_0 * lk_1461[k];

        t_1138[k] = -3.0 * ik_850[k]
                    + f_0 * lk_1462[k];

        t_1139[k] = -3.0 * ik_851[k]
                    + f_0 * lk_1463[k];

        t_1140[k] = -3.0 * ik_852[k]
                    + f_0 * lk_1464[k];
    }

#pragma omp simd aligned(t_1141, t_1142, t_1143, t_1144, t_1145, ik_853, ik_854, ik_855, \
                         ik_856, ik_857, lk_1465, lk_1466, lk_1467, lk_1468, \
                         lk_1469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1141[k] = -3.0 * ik_853[k]
                    + f_0 * lk_1465[k];

        t_1142[k] = -3.0 * ik_854[k]
                    + f_0 * lk_1466[k];

        t_1143[k] = -3.0 * ik_855[k]
                    + f_0 * lk_1467[k];

        t_1144[k] = -3.0 * ik_856[k]
                    + f_0 * lk_1468[k];

        t_1145[k] = -3.0 * ik_857[k]
                    + f_0 * lk_1469[k];
    }
}

static auto
compute_prim_geom_10_kk_electron_repulsion_2_piece7(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t lk,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_1260 = buffer.data(target + 1260);
    auto *t_1261 = buffer.data(target + 1261);
    auto *t_1262 = buffer.data(target + 1262);
    auto *t_1263 = buffer.data(target + 1263);
    auto *t_1264 = buffer.data(target + 1264);
    auto *t_1265 = buffer.data(target + 1265);
    auto *t_1266 = buffer.data(target + 1266);
    auto *t_1267 = buffer.data(target + 1267);
    auto *t_1268 = buffer.data(target + 1268);
    auto *t_1269 = buffer.data(target + 1269);
    auto *t_1270 = buffer.data(target + 1270);
    auto *t_1271 = buffer.data(target + 1271);
    auto *t_1272 = buffer.data(target + 1272);
    auto *t_1273 = buffer.data(target + 1273);
    auto *t_1274 = buffer.data(target + 1274);
    auto *t_1275 = buffer.data(target + 1275);
    auto *t_1276 = buffer.data(target + 1276);
    auto *t_1277 = buffer.data(target + 1277);
    auto *t_1278 = buffer.data(target + 1278);
    auto *t_1279 = buffer.data(target + 1279);
    auto *t_1280 = buffer.data(target + 1280);
    auto *t_1281 = buffer.data(target + 1281);
    auto *t_1282 = buffer.data(target + 1282);
    auto *t_1283 = buffer.data(target + 1283);
    auto *t_1284 = buffer.data(target + 1284);
    auto *t_1285 = buffer.data(target + 1285);
    auto *t_1286 = buffer.data(target + 1286);
    auto *t_1287 = buffer.data(target + 1287);
    auto *t_1288 = buffer.data(target + 1288);
    auto *t_1289 = buffer.data(target + 1289);
    auto *t_1290 = buffer.data(target + 1290);
    auto *t_1291 = buffer.data(target + 1291);
    auto *t_1292 = buffer.data(target + 1292);
    auto *t_1293 = buffer.data(target + 1293);
    auto *t_1294 = buffer.data(target + 1294);
    auto *t_1295 = buffer.data(target + 1295);

    const auto *ik_858 = buffer.data(ik + 858);
    const auto *ik_859 = buffer.data(ik + 859);
    const auto *ik_860 = buffer.data(ik + 860);
    const auto *ik_861 = buffer.data(ik + 861);
    const auto *ik_862 = buffer.data(ik + 862);
    const auto *ik_863 = buffer.data(ik + 863);
    const auto *ik_864 = buffer.data(ik + 864);
    const auto *ik_865 = buffer.data(ik + 865);
    const auto *ik_866 = buffer.data(ik + 866);
    const auto *ik_867 = buffer.data(ik + 867);
    const auto *ik_868 = buffer.data(ik + 868);
    const auto *ik_869 = buffer.data(ik + 869);
    const auto *ik_870 = buffer.data(ik + 870);
    const auto *ik_871 = buffer.data(ik + 871);
    const auto *ik_872 = buffer.data(ik + 872);
    const auto *ik_873 = buffer.data(ik + 873);
    const auto *ik_874 = buffer.data(ik + 874);
    const auto *ik_875 = buffer.data(ik + 875);
    const auto *ik_876 = buffer.data(ik + 876);
    const auto *ik_877 = buffer.data(ik + 877);
    const auto *ik_878 = buffer.data(ik + 878);
    const auto *ik_879 = buffer.data(ik + 879);
    const auto *ik_880 = buffer.data(ik + 880);
    const auto *ik_881 = buffer.data(ik + 881);
    const auto *ik_882 = buffer.data(ik + 882);
    const auto *ik_883 = buffer.data(ik + 883);
    const auto *ik_884 = buffer.data(ik + 884);
    const auto *ik_885 = buffer.data(ik + 885);
    const auto *ik_886 = buffer.data(ik + 886);
    const auto *ik_887 = buffer.data(ik + 887);
    const auto *ik_888 = buffer.data(ik + 888);
    const auto *ik_889 = buffer.data(ik + 889);
    const auto *ik_890 = buffer.data(ik + 890);
    const auto *ik_891 = buffer.data(ik + 891);
    const auto *ik_892 = buffer.data(ik + 892);
    const auto *ik_893 = buffer.data(ik + 893);
    const auto *ik_894 = buffer.data(ik + 894);
    const auto *ik_895 = buffer.data(ik + 895);
    const auto *ik_896 = buffer.data(ik + 896);
    const auto *ik_897 = buffer.data(ik + 897);
    const auto *ik_898 = buffer.data(ik + 898);
    const auto *ik_899 = buffer.data(ik + 899);
    const auto *ik_900 = buffer.data(ik + 900);
    const auto *ik_901 = buffer.data(ik + 901);
    const auto *ik_902 = buffer.data(ik + 902);
    const auto *ik_903 = buffer.data(ik + 903);
    const auto *ik_904 = buffer.data(ik + 904);
    const auto *ik_905 = buffer.data(ik + 905);
    const auto *ik_906 = buffer.data(ik + 906);
    const auto *ik_907 = buffer.data(ik + 907);
    const auto *ik_908 = buffer.data(ik + 908);
    const auto *ik_909 = buffer.data(ik + 909);
    const auto *ik_910 = buffer.data(ik + 910);
    const auto *ik_911 = buffer.data(ik + 911);
    const auto *ik_912 = buffer.data(ik + 912);
    const auto *ik_913 = buffer.data(ik + 913);
    const auto *ik_914 = buffer.data(ik + 914);
    const auto *ik_915 = buffer.data(ik + 915);
    const auto *ik_916 = buffer.data(ik + 916);
    const auto *ik_917 = buffer.data(ik + 917);
    const auto *ik_918 = buffer.data(ik + 918);
    const auto *ik_919 = buffer.data(ik + 919);
    const auto *ik_920 = buffer.data(ik + 920);
    const auto *ik_921 = buffer.data(ik + 921);
    const auto *ik_922 = buffer.data(ik + 922);
    const auto *ik_923 = buffer.data(ik + 923);
    const auto *ik_924 = buffer.data(ik + 924);
    const auto *ik_925 = buffer.data(ik + 925);
    const auto *ik_926 = buffer.data(ik + 926);
    const auto *ik_927 = buffer.data(ik + 927);
    const auto *ik_928 = buffer.data(ik + 928);
    const auto *ik_929 = buffer.data(ik + 929);
    const auto *ik_930 = buffer.data(ik + 930);
    const auto *ik_931 = buffer.data(ik + 931);
    const auto *ik_932 = buffer.data(ik + 932);
    const auto *ik_933 = buffer.data(ik + 933);
    const auto *ik_934 = buffer.data(ik + 934);
    const auto *ik_935 = buffer.data(ik + 935);
    const auto *ik_936 = buffer.data(ik + 936);
    const auto *ik_937 = buffer.data(ik + 937);
    const auto *ik_938 = buffer.data(ik + 938);
    const auto *ik_939 = buffer.data(ik + 939);
    const auto *ik_940 = buffer.data(ik + 940);
    const auto *ik_941 = buffer.data(ik + 941);
    const auto *ik_942 = buffer.data(ik + 942);
    const auto *ik_943 = buffer.data(ik + 943);
    const auto *ik_944 = buffer.data(ik + 944);
    const auto *ik_945 = buffer.data(ik + 945);
    const auto *ik_946 = buffer.data(ik + 946);
    const auto *ik_947 = buffer.data(ik + 947);
    const auto *ik_948 = buffer.data(ik + 948);
    const auto *ik_949 = buffer.data(ik + 949);
    const auto *ik_950 = buffer.data(ik + 950);
    const auto *ik_951 = buffer.data(ik + 951);
    const auto *ik_952 = buffer.data(ik + 952);
    const auto *ik_953 = buffer.data(ik + 953);
    const auto *ik_954 = buffer.data(ik + 954);
    const auto *ik_955 = buffer.data(ik + 955);
    const auto *ik_956 = buffer.data(ik + 956);
    const auto *ik_957 = buffer.data(ik + 957);
    const auto *ik_958 = buffer.data(ik + 958);
    const auto *ik_959 = buffer.data(ik + 959);
    const auto *ik_960 = buffer.data(ik + 960);
    const auto *ik_961 = buffer.data(ik + 961);
    const auto *ik_962 = buffer.data(ik + 962);
    const auto *ik_963 = buffer.data(ik + 963);
    const auto *ik_964 = buffer.data(ik + 964);
    const auto *ik_965 = buffer.data(ik + 965);
    const auto *ik_966 = buffer.data(ik + 966);
    const auto *ik_967 = buffer.data(ik + 967);
    const auto *ik_968 = buffer.data(ik + 968);
    const auto *ik_969 = buffer.data(ik + 969);
    const auto *ik_970 = buffer.data(ik + 970);
    const auto *ik_971 = buffer.data(ik + 971);
    const auto *ik_972 = buffer.data(ik + 972);
    const auto *ik_973 = buffer.data(ik + 973);
    const auto *ik_974 = buffer.data(ik + 974);
    const auto *ik_975 = buffer.data(ik + 975);
    const auto *ik_976 = buffer.data(ik + 976);
    const auto *ik_977 = buffer.data(ik + 977);
    const auto *ik_978 = buffer.data(ik + 978);
    const auto *ik_979 = buffer.data(ik + 979);
    const auto *ik_980 = buffer.data(ik + 980);
    const auto *ik_981 = buffer.data(ik + 981);
    const auto *ik_982 = buffer.data(ik + 982);
    const auto *ik_983 = buffer.data(ik + 983);
    const auto *ik_984 = buffer.data(ik + 984);
    const auto *ik_985 = buffer.data(ik + 985);
    const auto *ik_986 = buffer.data(ik + 986);
    const auto *ik_987 = buffer.data(ik + 987);
    const auto *ik_988 = buffer.data(ik + 988);
    const auto *ik_989 = buffer.data(ik + 989);
    const auto *ik_990 = buffer.data(ik + 990);
    const auto *ik_991 = buffer.data(ik + 991);
    const auto *ik_992 = buffer.data(ik + 992);
    const auto *ik_993 = buffer.data(ik + 993);
    const auto *ik_994 = buffer.data(ik + 994);
    const auto *ik_995 = buffer.data(ik + 995);
    const auto *ik_996 = buffer.data(ik + 996);
    const auto *ik_997 = buffer.data(ik + 997);
    const auto *ik_998 = buffer.data(ik + 998);
    const auto *ik_999 = buffer.data(ik + 999);
    const auto *ik_1000 = buffer.data(ik + 1000);
    const auto *ik_1001 = buffer.data(ik + 1001);
    const auto *ik_1002 = buffer.data(ik + 1002);
    const auto *ik_1003 = buffer.data(ik + 1003);
    const auto *ik_1004 = buffer.data(ik + 1004);
    const auto *ik_1005 = buffer.data(ik + 1005);
    const auto *ik_1006 = buffer.data(ik + 1006);
    const auto *ik_1007 = buffer.data(ik + 1007);

    const auto *lk_1470 = buffer.data(lk + 1470);
    const auto *lk_1471 = buffer.data(lk + 1471);
    const auto *lk_1472 = buffer.data(lk + 1472);
    const auto *lk_1473 = buffer.data(lk + 1473);
    const auto *lk_1474 = buffer.data(lk + 1474);
    const auto *lk_1475 = buffer.data(lk + 1475);
    const auto *lk_1476 = buffer.data(lk + 1476);
    const auto *lk_1477 = buffer.data(lk + 1477);
    const auto *lk_1478 = buffer.data(lk + 1478);
    const auto *lk_1479 = buffer.data(lk + 1479);
    const auto *lk_1480 = buffer.data(lk + 1480);
    const auto *lk_1481 = buffer.data(lk + 1481);
    const auto *lk_1482 = buffer.data(lk + 1482);
    const auto *lk_1483 = buffer.data(lk + 1483);
    const auto *lk_1484 = buffer.data(lk + 1484);
    const auto *lk_1485 = buffer.data(lk + 1485);
    const auto *lk_1486 = buffer.data(lk + 1486);
    const auto *lk_1487 = buffer.data(lk + 1487);
    const auto *lk_1488 = buffer.data(lk + 1488);
    const auto *lk_1489 = buffer.data(lk + 1489);
    const auto *lk_1490 = buffer.data(lk + 1490);
    const auto *lk_1491 = buffer.data(lk + 1491);
    const auto *lk_1492 = buffer.data(lk + 1492);
    const auto *lk_1493 = buffer.data(lk + 1493);
    const auto *lk_1494 = buffer.data(lk + 1494);
    const auto *lk_1495 = buffer.data(lk + 1495);
    const auto *lk_1496 = buffer.data(lk + 1496);
    const auto *lk_1497 = buffer.data(lk + 1497);
    const auto *lk_1498 = buffer.data(lk + 1498);
    const auto *lk_1499 = buffer.data(lk + 1499);
    const auto *lk_1500 = buffer.data(lk + 1500);
    const auto *lk_1501 = buffer.data(lk + 1501);
    const auto *lk_1502 = buffer.data(lk + 1502);
    const auto *lk_1503 = buffer.data(lk + 1503);
    const auto *lk_1504 = buffer.data(lk + 1504);
    const auto *lk_1505 = buffer.data(lk + 1505);
    const auto *lk_1506 = buffer.data(lk + 1506);
    const auto *lk_1507 = buffer.data(lk + 1507);
    const auto *lk_1508 = buffer.data(lk + 1508);
    const auto *lk_1509 = buffer.data(lk + 1509);
    const auto *lk_1510 = buffer.data(lk + 1510);
    const auto *lk_1511 = buffer.data(lk + 1511);
    const auto *lk_1512 = buffer.data(lk + 1512);
    const auto *lk_1513 = buffer.data(lk + 1513);
    const auto *lk_1514 = buffer.data(lk + 1514);
    const auto *lk_1515 = buffer.data(lk + 1515);
    const auto *lk_1516 = buffer.data(lk + 1516);
    const auto *lk_1517 = buffer.data(lk + 1517);
    const auto *lk_1518 = buffer.data(lk + 1518);
    const auto *lk_1519 = buffer.data(lk + 1519);
    const auto *lk_1520 = buffer.data(lk + 1520);
    const auto *lk_1521 = buffer.data(lk + 1521);
    const auto *lk_1522 = buffer.data(lk + 1522);
    const auto *lk_1523 = buffer.data(lk + 1523);
    const auto *lk_1524 = buffer.data(lk + 1524);
    const auto *lk_1525 = buffer.data(lk + 1525);
    const auto *lk_1526 = buffer.data(lk + 1526);
    const auto *lk_1527 = buffer.data(lk + 1527);
    const auto *lk_1528 = buffer.data(lk + 1528);
    const auto *lk_1529 = buffer.data(lk + 1529);
    const auto *lk_1530 = buffer.data(lk + 1530);
    const auto *lk_1531 = buffer.data(lk + 1531);
    const auto *lk_1532 = buffer.data(lk + 1532);
    const auto *lk_1533 = buffer.data(lk + 1533);
    const auto *lk_1534 = buffer.data(lk + 1534);
    const auto *lk_1535 = buffer.data(lk + 1535);
    const auto *lk_1536 = buffer.data(lk + 1536);
    const auto *lk_1537 = buffer.data(lk + 1537);
    const auto *lk_1538 = buffer.data(lk + 1538);
    const auto *lk_1539 = buffer.data(lk + 1539);
    const auto *lk_1540 = buffer.data(lk + 1540);
    const auto *lk_1541 = buffer.data(lk + 1541);
    const auto *lk_1542 = buffer.data(lk + 1542);
    const auto *lk_1543 = buffer.data(lk + 1543);
    const auto *lk_1544 = buffer.data(lk + 1544);
    const auto *lk_1545 = buffer.data(lk + 1545);
    const auto *lk_1546 = buffer.data(lk + 1546);
    const auto *lk_1547 = buffer.data(lk + 1547);
    const auto *lk_1548 = buffer.data(lk + 1548);
    const auto *lk_1549 = buffer.data(lk + 1549);
    const auto *lk_1550 = buffer.data(lk + 1550);
    const auto *lk_1551 = buffer.data(lk + 1551);
    const auto *lk_1552 = buffer.data(lk + 1552);
    const auto *lk_1553 = buffer.data(lk + 1553);
    const auto *lk_1554 = buffer.data(lk + 1554);
    const auto *lk_1555 = buffer.data(lk + 1555);
    const auto *lk_1556 = buffer.data(lk + 1556);
    const auto *lk_1557 = buffer.data(lk + 1557);
    const auto *lk_1558 = buffer.data(lk + 1558);
    const auto *lk_1559 = buffer.data(lk + 1559);
    const auto *lk_1560 = buffer.data(lk + 1560);
    const auto *lk_1561 = buffer.data(lk + 1561);
    const auto *lk_1562 = buffer.data(lk + 1562);
    const auto *lk_1563 = buffer.data(lk + 1563);
    const auto *lk_1564 = buffer.data(lk + 1564);
    const auto *lk_1565 = buffer.data(lk + 1565);
    const auto *lk_1566 = buffer.data(lk + 1566);
    const auto *lk_1567 = buffer.data(lk + 1567);
    const auto *lk_1568 = buffer.data(lk + 1568);
    const auto *lk_1569 = buffer.data(lk + 1569);
    const auto *lk_1570 = buffer.data(lk + 1570);
    const auto *lk_1571 = buffer.data(lk + 1571);
    const auto *lk_1572 = buffer.data(lk + 1572);
    const auto *lk_1573 = buffer.data(lk + 1573);
    const auto *lk_1574 = buffer.data(lk + 1574);
    const auto *lk_1575 = buffer.data(lk + 1575);
    const auto *lk_1576 = buffer.data(lk + 1576);
    const auto *lk_1577 = buffer.data(lk + 1577);
    const auto *lk_1578 = buffer.data(lk + 1578);
    const auto *lk_1579 = buffer.data(lk + 1579);
    const auto *lk_1580 = buffer.data(lk + 1580);
    const auto *lk_1581 = buffer.data(lk + 1581);
    const auto *lk_1582 = buffer.data(lk + 1582);
    const auto *lk_1583 = buffer.data(lk + 1583);
    const auto *lk_1584 = buffer.data(lk + 1584);
    const auto *lk_1585 = buffer.data(lk + 1585);
    const auto *lk_1586 = buffer.data(lk + 1586);
    const auto *lk_1587 = buffer.data(lk + 1587);
    const auto *lk_1588 = buffer.data(lk + 1588);
    const auto *lk_1589 = buffer.data(lk + 1589);
    const auto *lk_1590 = buffer.data(lk + 1590);
    const auto *lk_1591 = buffer.data(lk + 1591);
    const auto *lk_1592 = buffer.data(lk + 1592);
    const auto *lk_1593 = buffer.data(lk + 1593);
    const auto *lk_1594 = buffer.data(lk + 1594);
    const auto *lk_1595 = buffer.data(lk + 1595);
    const auto *lk_1596 = buffer.data(lk + 1596);
    const auto *lk_1597 = buffer.data(lk + 1597);
    const auto *lk_1598 = buffer.data(lk + 1598);
    const auto *lk_1599 = buffer.data(lk + 1599);
    const auto *lk_1600 = buffer.data(lk + 1600);
    const auto *lk_1601 = buffer.data(lk + 1601);
    const auto *lk_1602 = buffer.data(lk + 1602);
    const auto *lk_1603 = buffer.data(lk + 1603);
    const auto *lk_1604 = buffer.data(lk + 1604);
    const auto *lk_1605 = buffer.data(lk + 1605);
    const auto *lk_1606 = buffer.data(lk + 1606);
    const auto *lk_1607 = buffer.data(lk + 1607);
    const auto *lk_1608 = buffer.data(lk + 1608);
    const auto *lk_1609 = buffer.data(lk + 1609);
    const auto *lk_1610 = buffer.data(lk + 1610);
    const auto *lk_1611 = buffer.data(lk + 1611);
    const auto *lk_1612 = buffer.data(lk + 1612);
    const auto *lk_1613 = buffer.data(lk + 1613);
    const auto *lk_1614 = buffer.data(lk + 1614);
    const auto *lk_1615 = buffer.data(lk + 1615);
    const auto *lk_1616 = buffer.data(lk + 1616);
    const auto *lk_1617 = buffer.data(lk + 1617);
    const auto *lk_1618 = buffer.data(lk + 1618);
    const auto *lk_1619 = buffer.data(lk + 1619);

#pragma omp simd aligned(t_1146, t_1147, t_1148, t_1149, t_1150, ik_858, ik_859, ik_860, \
                         ik_861, ik_862, lk_1470, lk_1471, lk_1472, lk_1473, \
                         lk_1474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = -3.0 * ik_858[k]
                    + f_0 * lk_1470[k];

        t_1147[k] = -3.0 * ik_859[k]
                    + f_0 * lk_1471[k];

        t_1148[k] = -3.0 * ik_860[k]
                    + f_0 * lk_1472[k];

        t_1149[k] = -3.0 * ik_861[k]
                    + f_0 * lk_1473[k];

        t_1150[k] = -3.0 * ik_862[k]
                    + f_0 * lk_1474[k];
    }

#pragma omp simd aligned(t_1151, t_1152, t_1153, t_1154, t_1155, ik_863, ik_864, ik_865, \
                         ik_866, ik_867, lk_1475, lk_1476, lk_1477, lk_1478, \
                         lk_1479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = -3.0 * ik_863[k]
                    + f_0 * lk_1475[k];

        t_1152[k] = -4.0 * ik_864[k]
                    + f_0 * lk_1476[k];

        t_1153[k] = -4.0 * ik_865[k]
                    + f_0 * lk_1477[k];

        t_1154[k] = -4.0 * ik_866[k]
                    + f_0 * lk_1478[k];

        t_1155[k] = -4.0 * ik_867[k]
                    + f_0 * lk_1479[k];
    }

#pragma omp simd aligned(t_1156, t_1157, t_1158, t_1159, t_1160, ik_868, ik_869, ik_870, \
                         ik_871, ik_872, lk_1480, lk_1481, lk_1482, lk_1483, \
                         lk_1484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1156[k] = -4.0 * ik_868[k]
                    + f_0 * lk_1480[k];

        t_1157[k] = -4.0 * ik_869[k]
                    + f_0 * lk_1481[k];

        t_1158[k] = -4.0 * ik_870[k]
                    + f_0 * lk_1482[k];

        t_1159[k] = -4.0 * ik_871[k]
                    + f_0 * lk_1483[k];

        t_1160[k] = -4.0 * ik_872[k]
                    + f_0 * lk_1484[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, t_1164, t_1165, ik_873, ik_874, ik_875, \
                         ik_876, ik_877, lk_1485, lk_1486, lk_1487, lk_1488, \
                         lk_1489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = -4.0 * ik_873[k]
                    + f_0 * lk_1485[k];

        t_1162[k] = -4.0 * ik_874[k]
                    + f_0 * lk_1486[k];

        t_1163[k] = -4.0 * ik_875[k]
                    + f_0 * lk_1487[k];

        t_1164[k] = -4.0 * ik_876[k]
                    + f_0 * lk_1488[k];

        t_1165[k] = -4.0 * ik_877[k]
                    + f_0 * lk_1489[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, t_1169, t_1170, ik_878, ik_879, ik_880, \
                         ik_881, ik_882, lk_1490, lk_1491, lk_1492, lk_1493, \
                         lk_1494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = -4.0 * ik_878[k]
                    + f_0 * lk_1490[k];

        t_1167[k] = -4.0 * ik_879[k]
                    + f_0 * lk_1491[k];

        t_1168[k] = -4.0 * ik_880[k]
                    + f_0 * lk_1492[k];

        t_1169[k] = -4.0 * ik_881[k]
                    + f_0 * lk_1493[k];

        t_1170[k] = -4.0 * ik_882[k]
                    + f_0 * lk_1494[k];
    }

#pragma omp simd aligned(t_1171, t_1172, t_1173, t_1174, t_1175, ik_883, ik_884, ik_885, \
                         ik_886, ik_887, lk_1495, lk_1496, lk_1497, lk_1498, \
                         lk_1499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1171[k] = -4.0 * ik_883[k]
                    + f_0 * lk_1495[k];

        t_1172[k] = -4.0 * ik_884[k]
                    + f_0 * lk_1496[k];

        t_1173[k] = -4.0 * ik_885[k]
                    + f_0 * lk_1497[k];

        t_1174[k] = -4.0 * ik_886[k]
                    + f_0 * lk_1498[k];

        t_1175[k] = -4.0 * ik_887[k]
                    + f_0 * lk_1499[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, t_1179, t_1180, ik_888, ik_889, ik_890, \
                         ik_891, ik_892, lk_1500, lk_1501, lk_1502, lk_1503, \
                         lk_1504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = -4.0 * ik_888[k]
                    + f_0 * lk_1500[k];

        t_1177[k] = -4.0 * ik_889[k]
                    + f_0 * lk_1501[k];

        t_1178[k] = -4.0 * ik_890[k]
                    + f_0 * lk_1502[k];

        t_1179[k] = -4.0 * ik_891[k]
                    + f_0 * lk_1503[k];

        t_1180[k] = -4.0 * ik_892[k]
                    + f_0 * lk_1504[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, t_1184, t_1185, ik_893, ik_894, ik_895, \
                         ik_896, ik_897, lk_1505, lk_1506, lk_1507, lk_1508, \
                         lk_1509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = -4.0 * ik_893[k]
                    + f_0 * lk_1505[k];

        t_1182[k] = -4.0 * ik_894[k]
                    + f_0 * lk_1506[k];

        t_1183[k] = -4.0 * ik_895[k]
                    + f_0 * lk_1507[k];

        t_1184[k] = -4.0 * ik_896[k]
                    + f_0 * lk_1508[k];

        t_1185[k] = -4.0 * ik_897[k]
                    + f_0 * lk_1509[k];
    }

#pragma omp simd aligned(t_1186, t_1187, t_1188, t_1189, t_1190, ik_898, ik_899, ik_900, \
                         ik_901, ik_902, lk_1510, lk_1511, lk_1512, lk_1513, \
                         lk_1514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1186[k] = -4.0 * ik_898[k]
                    + f_0 * lk_1510[k];

        t_1187[k] = -4.0 * ik_899[k]
                    + f_0 * lk_1511[k];

        t_1188[k] = -5.0 * ik_900[k]
                    + f_0 * lk_1512[k];

        t_1189[k] = -5.0 * ik_901[k]
                    + f_0 * lk_1513[k];

        t_1190[k] = -5.0 * ik_902[k]
                    + f_0 * lk_1514[k];
    }

#pragma omp simd aligned(t_1191, t_1192, t_1193, t_1194, t_1195, ik_903, ik_904, ik_905, \
                         ik_906, ik_907, lk_1515, lk_1516, lk_1517, lk_1518, \
                         lk_1519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1191[k] = -5.0 * ik_903[k]
                    + f_0 * lk_1515[k];

        t_1192[k] = -5.0 * ik_904[k]
                    + f_0 * lk_1516[k];

        t_1193[k] = -5.0 * ik_905[k]
                    + f_0 * lk_1517[k];

        t_1194[k] = -5.0 * ik_906[k]
                    + f_0 * lk_1518[k];

        t_1195[k] = -5.0 * ik_907[k]
                    + f_0 * lk_1519[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, t_1199, t_1200, ik_908, ik_909, ik_910, \
                         ik_911, ik_912, lk_1520, lk_1521, lk_1522, lk_1523, \
                         lk_1524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = -5.0 * ik_908[k]
                    + f_0 * lk_1520[k];

        t_1197[k] = -5.0 * ik_909[k]
                    + f_0 * lk_1521[k];

        t_1198[k] = -5.0 * ik_910[k]
                    + f_0 * lk_1522[k];

        t_1199[k] = -5.0 * ik_911[k]
                    + f_0 * lk_1523[k];

        t_1200[k] = -5.0 * ik_912[k]
                    + f_0 * lk_1524[k];
    }

#pragma omp simd aligned(t_1201, t_1202, t_1203, t_1204, t_1205, ik_913, ik_914, ik_915, \
                         ik_916, ik_917, lk_1525, lk_1526, lk_1527, lk_1528, \
                         lk_1529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1201[k] = -5.0 * ik_913[k]
                    + f_0 * lk_1525[k];

        t_1202[k] = -5.0 * ik_914[k]
                    + f_0 * lk_1526[k];

        t_1203[k] = -5.0 * ik_915[k]
                    + f_0 * lk_1527[k];

        t_1204[k] = -5.0 * ik_916[k]
                    + f_0 * lk_1528[k];

        t_1205[k] = -5.0 * ik_917[k]
                    + f_0 * lk_1529[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, t_1210, ik_918, ik_919, ik_920, \
                         ik_921, ik_922, lk_1530, lk_1531, lk_1532, lk_1533, \
                         lk_1534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = -5.0 * ik_918[k]
                    + f_0 * lk_1530[k];

        t_1207[k] = -5.0 * ik_919[k]
                    + f_0 * lk_1531[k];

        t_1208[k] = -5.0 * ik_920[k]
                    + f_0 * lk_1532[k];

        t_1209[k] = -5.0 * ik_921[k]
                    + f_0 * lk_1533[k];

        t_1210[k] = -5.0 * ik_922[k]
                    + f_0 * lk_1534[k];
    }

#pragma omp simd aligned(t_1211, t_1212, t_1213, t_1214, t_1215, ik_923, ik_924, ik_925, \
                         ik_926, ik_927, lk_1535, lk_1536, lk_1537, lk_1538, \
                         lk_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1211[k] = -5.0 * ik_923[k]
                    + f_0 * lk_1535[k];

        t_1212[k] = -5.0 * ik_924[k]
                    + f_0 * lk_1536[k];

        t_1213[k] = -5.0 * ik_925[k]
                    + f_0 * lk_1537[k];

        t_1214[k] = -5.0 * ik_926[k]
                    + f_0 * lk_1538[k];

        t_1215[k] = -5.0 * ik_927[k]
                    + f_0 * lk_1539[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, t_1219, t_1220, ik_928, ik_929, ik_930, \
                         ik_931, ik_932, lk_1540, lk_1541, lk_1542, lk_1543, \
                         lk_1544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = -5.0 * ik_928[k]
                    + f_0 * lk_1540[k];

        t_1217[k] = -5.0 * ik_929[k]
                    + f_0 * lk_1541[k];

        t_1218[k] = -5.0 * ik_930[k]
                    + f_0 * lk_1542[k];

        t_1219[k] = -5.0 * ik_931[k]
                    + f_0 * lk_1543[k];

        t_1220[k] = -5.0 * ik_932[k]
                    + f_0 * lk_1544[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, t_1224, t_1225, ik_933, ik_934, ik_935, \
                         ik_936, ik_937, lk_1545, lk_1546, lk_1547, lk_1548, \
                         lk_1549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = -5.0 * ik_933[k]
                    + f_0 * lk_1545[k];

        t_1222[k] = -5.0 * ik_934[k]
                    + f_0 * lk_1546[k];

        t_1223[k] = -5.0 * ik_935[k]
                    + f_0 * lk_1547[k];

        t_1224[k] = -6.0 * ik_936[k]
                    + f_0 * lk_1548[k];

        t_1225[k] = -6.0 * ik_937[k]
                    + f_0 * lk_1549[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, t_1229, t_1230, ik_938, ik_939, ik_940, \
                         ik_941, ik_942, lk_1550, lk_1551, lk_1552, lk_1553, \
                         lk_1554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = -6.0 * ik_938[k]
                    + f_0 * lk_1550[k];

        t_1227[k] = -6.0 * ik_939[k]
                    + f_0 * lk_1551[k];

        t_1228[k] = -6.0 * ik_940[k]
                    + f_0 * lk_1552[k];

        t_1229[k] = -6.0 * ik_941[k]
                    + f_0 * lk_1553[k];

        t_1230[k] = -6.0 * ik_942[k]
                    + f_0 * lk_1554[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, t_1235, ik_943, ik_944, ik_945, \
                         ik_946, ik_947, lk_1555, lk_1556, lk_1557, lk_1558, \
                         lk_1559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = -6.0 * ik_943[k]
                    + f_0 * lk_1555[k];

        t_1232[k] = -6.0 * ik_944[k]
                    + f_0 * lk_1556[k];

        t_1233[k] = -6.0 * ik_945[k]
                    + f_0 * lk_1557[k];

        t_1234[k] = -6.0 * ik_946[k]
                    + f_0 * lk_1558[k];

        t_1235[k] = -6.0 * ik_947[k]
                    + f_0 * lk_1559[k];
    }

#pragma omp simd aligned(t_1236, t_1237, t_1238, t_1239, t_1240, ik_948, ik_949, ik_950, \
                         ik_951, ik_952, lk_1560, lk_1561, lk_1562, lk_1563, \
                         lk_1564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1236[k] = -6.0 * ik_948[k]
                    + f_0 * lk_1560[k];

        t_1237[k] = -6.0 * ik_949[k]
                    + f_0 * lk_1561[k];

        t_1238[k] = -6.0 * ik_950[k]
                    + f_0 * lk_1562[k];

        t_1239[k] = -6.0 * ik_951[k]
                    + f_0 * lk_1563[k];

        t_1240[k] = -6.0 * ik_952[k]
                    + f_0 * lk_1564[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, t_1244, t_1245, ik_953, ik_954, ik_955, \
                         ik_956, ik_957, lk_1565, lk_1566, lk_1567, lk_1568, \
                         lk_1569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = -6.0 * ik_953[k]
                    + f_0 * lk_1565[k];

        t_1242[k] = -6.0 * ik_954[k]
                    + f_0 * lk_1566[k];

        t_1243[k] = -6.0 * ik_955[k]
                    + f_0 * lk_1567[k];

        t_1244[k] = -6.0 * ik_956[k]
                    + f_0 * lk_1568[k];

        t_1245[k] = -6.0 * ik_957[k]
                    + f_0 * lk_1569[k];
    }

#pragma omp simd aligned(t_1246, t_1247, t_1248, t_1249, t_1250, ik_958, ik_959, ik_960, \
                         ik_961, ik_962, lk_1570, lk_1571, lk_1572, lk_1573, \
                         lk_1574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1246[k] = -6.0 * ik_958[k]
                    + f_0 * lk_1570[k];

        t_1247[k] = -6.0 * ik_959[k]
                    + f_0 * lk_1571[k];

        t_1248[k] = -6.0 * ik_960[k]
                    + f_0 * lk_1572[k];

        t_1249[k] = -6.0 * ik_961[k]
                    + f_0 * lk_1573[k];

        t_1250[k] = -6.0 * ik_962[k]
                    + f_0 * lk_1574[k];
    }

#pragma omp simd aligned(t_1251, t_1252, t_1253, t_1254, t_1255, ik_963, ik_964, ik_965, \
                         ik_966, ik_967, lk_1575, lk_1576, lk_1577, lk_1578, \
                         lk_1579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1251[k] = -6.0 * ik_963[k]
                    + f_0 * lk_1575[k];

        t_1252[k] = -6.0 * ik_964[k]
                    + f_0 * lk_1576[k];

        t_1253[k] = -6.0 * ik_965[k]
                    + f_0 * lk_1577[k];

        t_1254[k] = -6.0 * ik_966[k]
                    + f_0 * lk_1578[k];

        t_1255[k] = -6.0 * ik_967[k]
                    + f_0 * lk_1579[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, t_1260, ik_968, ik_969, ik_970, \
                         ik_971, ik_972, lk_1580, lk_1581, lk_1582, lk_1583, \
                         lk_1584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = -6.0 * ik_968[k]
                    + f_0 * lk_1580[k];

        t_1257[k] = -6.0 * ik_969[k]
                    + f_0 * lk_1581[k];

        t_1258[k] = -6.0 * ik_970[k]
                    + f_0 * lk_1582[k];

        t_1259[k] = -6.0 * ik_971[k]
                    + f_0 * lk_1583[k];

        t_1260[k] = -7.0 * ik_972[k]
                    + f_0 * lk_1584[k];
    }

#pragma omp simd aligned(t_1261, t_1262, t_1263, t_1264, t_1265, ik_973, ik_974, ik_975, \
                         ik_976, ik_977, lk_1585, lk_1586, lk_1587, lk_1588, \
                         lk_1589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1261[k] = -7.0 * ik_973[k]
                    + f_0 * lk_1585[k];

        t_1262[k] = -7.0 * ik_974[k]
                    + f_0 * lk_1586[k];

        t_1263[k] = -7.0 * ik_975[k]
                    + f_0 * lk_1587[k];

        t_1264[k] = -7.0 * ik_976[k]
                    + f_0 * lk_1588[k];

        t_1265[k] = -7.0 * ik_977[k]
                    + f_0 * lk_1589[k];
    }

#pragma omp simd aligned(t_1266, t_1267, t_1268, t_1269, t_1270, ik_978, ik_979, ik_980, \
                         ik_981, ik_982, lk_1590, lk_1591, lk_1592, lk_1593, \
                         lk_1594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1266[k] = -7.0 * ik_978[k]
                    + f_0 * lk_1590[k];

        t_1267[k] = -7.0 * ik_979[k]
                    + f_0 * lk_1591[k];

        t_1268[k] = -7.0 * ik_980[k]
                    + f_0 * lk_1592[k];

        t_1269[k] = -7.0 * ik_981[k]
                    + f_0 * lk_1593[k];

        t_1270[k] = -7.0 * ik_982[k]
                    + f_0 * lk_1594[k];
    }

#pragma omp simd aligned(t_1271, t_1272, t_1273, t_1274, t_1275, ik_983, ik_984, ik_985, \
                         ik_986, ik_987, lk_1595, lk_1596, lk_1597, lk_1598, \
                         lk_1599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1271[k] = -7.0 * ik_983[k]
                    + f_0 * lk_1595[k];

        t_1272[k] = -7.0 * ik_984[k]
                    + f_0 * lk_1596[k];

        t_1273[k] = -7.0 * ik_985[k]
                    + f_0 * lk_1597[k];

        t_1274[k] = -7.0 * ik_986[k]
                    + f_0 * lk_1598[k];

        t_1275[k] = -7.0 * ik_987[k]
                    + f_0 * lk_1599[k];
    }

#pragma omp simd aligned(t_1276, t_1277, t_1278, t_1279, t_1280, ik_988, ik_989, ik_990, \
                         ik_991, ik_992, lk_1600, lk_1601, lk_1602, lk_1603, \
                         lk_1604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1276[k] = -7.0 * ik_988[k]
                    + f_0 * lk_1600[k];

        t_1277[k] = -7.0 * ik_989[k]
                    + f_0 * lk_1601[k];

        t_1278[k] = -7.0 * ik_990[k]
                    + f_0 * lk_1602[k];

        t_1279[k] = -7.0 * ik_991[k]
                    + f_0 * lk_1603[k];

        t_1280[k] = -7.0 * ik_992[k]
                    + f_0 * lk_1604[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, t_1284, t_1285, ik_993, ik_994, ik_995, \
                         ik_996, ik_997, lk_1605, lk_1606, lk_1607, lk_1608, \
                         lk_1609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = -7.0 * ik_993[k]
                    + f_0 * lk_1605[k];

        t_1282[k] = -7.0 * ik_994[k]
                    + f_0 * lk_1606[k];

        t_1283[k] = -7.0 * ik_995[k]
                    + f_0 * lk_1607[k];

        t_1284[k] = -7.0 * ik_996[k]
                    + f_0 * lk_1608[k];

        t_1285[k] = -7.0 * ik_997[k]
                    + f_0 * lk_1609[k];
    }

#pragma omp simd aligned(t_1286, t_1287, t_1288, t_1289, t_1290, ik_998, ik_999, ik_1000, \
                         ik_1001, ik_1002, lk_1610, lk_1611, lk_1612, lk_1613, \
                         lk_1614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1286[k] = -7.0 * ik_998[k]
                    + f_0 * lk_1610[k];

        t_1287[k] = -7.0 * ik_999[k]
                    + f_0 * lk_1611[k];

        t_1288[k] = -7.0 * ik_1000[k]
                    + f_0 * lk_1612[k];

        t_1289[k] = -7.0 * ik_1001[k]
                    + f_0 * lk_1613[k];

        t_1290[k] = -7.0 * ik_1002[k]
                    + f_0 * lk_1614[k];
    }

#pragma omp simd aligned(t_1291, t_1292, t_1293, t_1294, t_1295, ik_1003, ik_1004, ik_1005, \
                         ik_1006, ik_1007, lk_1615, lk_1616, lk_1617, lk_1618, \
                         lk_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1291[k] = -7.0 * ik_1003[k]
                    + f_0 * lk_1615[k];

        t_1292[k] = -7.0 * ik_1004[k]
                    + f_0 * lk_1616[k];

        t_1293[k] = -7.0 * ik_1005[k]
                    + f_0 * lk_1617[k];

        t_1294[k] = -7.0 * ik_1006[k]
                    + f_0 * lk_1618[k];

        t_1295[k] = -7.0 * ik_1007[k]
                    + f_0 * lk_1619[k];
    }
}

auto
compute_prim_geom_10_kk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t ik, const size_t lk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kk_electron_repulsion_2_piece0(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_2_piece1(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_2_piece2(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_2_piece3(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_2_piece4(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_2_piece5(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_2_piece6(buffer, target, ik, lk, ncols, alpha);

    compute_prim_geom_10_kk_electron_repulsion_2_piece7(buffer, target, ik, lk, ncols, alpha);
}

}  // namespace simdt2ceri
