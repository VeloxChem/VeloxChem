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


#include "SimdElectronRepulsionGeom10VrrRecHK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_hk_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, gk_0, gk_1, gk_2, gk_3, gk_4, ik_0, ik_1, \
                         ik_2, ik_3, ik_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -5.0 * gk_0[k]
                 + f_0 * ik_0[k];

        t_1[k] = -5.0 * gk_1[k]
                 + f_0 * ik_1[k];

        t_2[k] = -5.0 * gk_2[k]
                 + f_0 * ik_2[k];

        t_3[k] = -5.0 * gk_3[k]
                 + f_0 * ik_3[k];

        t_4[k] = -5.0 * gk_4[k]
                 + f_0 * ik_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, gk_5, gk_6, gk_7, gk_8, gk_9, ik_5, ik_6, \
                         ik_7, ik_8, ik_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -5.0 * gk_5[k]
                 + f_0 * ik_5[k];

        t_6[k] = -5.0 * gk_6[k]
                 + f_0 * ik_6[k];

        t_7[k] = -5.0 * gk_7[k]
                 + f_0 * ik_7[k];

        t_8[k] = -5.0 * gk_8[k]
                 + f_0 * ik_8[k];

        t_9[k] = -5.0 * gk_9[k]
                 + f_0 * ik_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, gk_10, gk_11, gk_12, gk_13, gk_14, \
                         ik_10, ik_11, ik_12, ik_13, ik_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -5.0 * gk_10[k]
                  + f_0 * ik_10[k];

        t_11[k] = -5.0 * gk_11[k]
                  + f_0 * ik_11[k];

        t_12[k] = -5.0 * gk_12[k]
                  + f_0 * ik_12[k];

        t_13[k] = -5.0 * gk_13[k]
                  + f_0 * ik_13[k];

        t_14[k] = -5.0 * gk_14[k]
                  + f_0 * ik_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, gk_15, gk_16, gk_17, gk_18, gk_19, \
                         ik_15, ik_16, ik_17, ik_18, ik_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -5.0 * gk_15[k]
                  + f_0 * ik_15[k];

        t_16[k] = -5.0 * gk_16[k]
                  + f_0 * ik_16[k];

        t_17[k] = -5.0 * gk_17[k]
                  + f_0 * ik_17[k];

        t_18[k] = -5.0 * gk_18[k]
                  + f_0 * ik_18[k];

        t_19[k] = -5.0 * gk_19[k]
                  + f_0 * ik_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, gk_20, gk_21, gk_22, gk_23, gk_24, \
                         ik_20, ik_21, ik_22, ik_23, ik_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -5.0 * gk_20[k]
                  + f_0 * ik_20[k];

        t_21[k] = -5.0 * gk_21[k]
                  + f_0 * ik_21[k];

        t_22[k] = -5.0 * gk_22[k]
                  + f_0 * ik_22[k];

        t_23[k] = -5.0 * gk_23[k]
                  + f_0 * ik_23[k];

        t_24[k] = -5.0 * gk_24[k]
                  + f_0 * ik_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, gk_25, gk_26, gk_27, gk_28, gk_29, \
                         ik_25, ik_26, ik_27, ik_28, ik_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -5.0 * gk_25[k]
                  + f_0 * ik_25[k];

        t_26[k] = -5.0 * gk_26[k]
                  + f_0 * ik_26[k];

        t_27[k] = -5.0 * gk_27[k]
                  + f_0 * ik_27[k];

        t_28[k] = -5.0 * gk_28[k]
                  + f_0 * ik_28[k];

        t_29[k] = -5.0 * gk_29[k]
                  + f_0 * ik_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, gk_30, gk_31, gk_32, gk_33, gk_34, \
                         ik_30, ik_31, ik_32, ik_33, ik_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -5.0 * gk_30[k]
                  + f_0 * ik_30[k];

        t_31[k] = -5.0 * gk_31[k]
                  + f_0 * ik_31[k];

        t_32[k] = -5.0 * gk_32[k]
                  + f_0 * ik_32[k];

        t_33[k] = -5.0 * gk_33[k]
                  + f_0 * ik_33[k];

        t_34[k] = -5.0 * gk_34[k]
                  + f_0 * ik_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, gk_35, gk_36, gk_37, gk_38, gk_39, \
                         ik_35, ik_36, ik_37, ik_38, ik_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -5.0 * gk_35[k]
                  + f_0 * ik_35[k];

        t_36[k] = -4.0 * gk_36[k]
                  + f_0 * ik_36[k];

        t_37[k] = -4.0 * gk_37[k]
                  + f_0 * ik_37[k];

        t_38[k] = -4.0 * gk_38[k]
                  + f_0 * ik_38[k];

        t_39[k] = -4.0 * gk_39[k]
                  + f_0 * ik_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, gk_40, gk_41, gk_42, gk_43, gk_44, \
                         ik_40, ik_41, ik_42, ik_43, ik_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -4.0 * gk_40[k]
                  + f_0 * ik_40[k];

        t_41[k] = -4.0 * gk_41[k]
                  + f_0 * ik_41[k];

        t_42[k] = -4.0 * gk_42[k]
                  + f_0 * ik_42[k];

        t_43[k] = -4.0 * gk_43[k]
                  + f_0 * ik_43[k];

        t_44[k] = -4.0 * gk_44[k]
                  + f_0 * ik_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, gk_45, gk_46, gk_47, gk_48, gk_49, \
                         ik_45, ik_46, ik_47, ik_48, ik_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -4.0 * gk_45[k]
                  + f_0 * ik_45[k];

        t_46[k] = -4.0 * gk_46[k]
                  + f_0 * ik_46[k];

        t_47[k] = -4.0 * gk_47[k]
                  + f_0 * ik_47[k];

        t_48[k] = -4.0 * gk_48[k]
                  + f_0 * ik_48[k];

        t_49[k] = -4.0 * gk_49[k]
                  + f_0 * ik_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, gk_50, gk_51, gk_52, gk_53, gk_54, \
                         ik_50, ik_51, ik_52, ik_53, ik_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -4.0 * gk_50[k]
                  + f_0 * ik_50[k];

        t_51[k] = -4.0 * gk_51[k]
                  + f_0 * ik_51[k];

        t_52[k] = -4.0 * gk_52[k]
                  + f_0 * ik_52[k];

        t_53[k] = -4.0 * gk_53[k]
                  + f_0 * ik_53[k];

        t_54[k] = -4.0 * gk_54[k]
                  + f_0 * ik_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, gk_55, gk_56, gk_57, gk_58, gk_59, \
                         ik_55, ik_56, ik_57, ik_58, ik_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -4.0 * gk_55[k]
                  + f_0 * ik_55[k];

        t_56[k] = -4.0 * gk_56[k]
                  + f_0 * ik_56[k];

        t_57[k] = -4.0 * gk_57[k]
                  + f_0 * ik_57[k];

        t_58[k] = -4.0 * gk_58[k]
                  + f_0 * ik_58[k];

        t_59[k] = -4.0 * gk_59[k]
                  + f_0 * ik_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, gk_60, gk_61, gk_62, gk_63, gk_64, \
                         ik_60, ik_61, ik_62, ik_63, ik_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -4.0 * gk_60[k]
                  + f_0 * ik_60[k];

        t_61[k] = -4.0 * gk_61[k]
                  + f_0 * ik_61[k];

        t_62[k] = -4.0 * gk_62[k]
                  + f_0 * ik_62[k];

        t_63[k] = -4.0 * gk_63[k]
                  + f_0 * ik_63[k];

        t_64[k] = -4.0 * gk_64[k]
                  + f_0 * ik_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, gk_65, gk_66, gk_67, gk_68, gk_69, \
                         ik_65, ik_66, ik_67, ik_68, ik_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -4.0 * gk_65[k]
                  + f_0 * ik_65[k];

        t_66[k] = -4.0 * gk_66[k]
                  + f_0 * ik_66[k];

        t_67[k] = -4.0 * gk_67[k]
                  + f_0 * ik_67[k];

        t_68[k] = -4.0 * gk_68[k]
                  + f_0 * ik_68[k];

        t_69[k] = -4.0 * gk_69[k]
                  + f_0 * ik_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, gk_70, gk_71, gk_72, gk_73, gk_74, \
                         ik_70, ik_71, ik_72, ik_73, ik_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -4.0 * gk_70[k]
                  + f_0 * ik_70[k];

        t_71[k] = -4.0 * gk_71[k]
                  + f_0 * ik_71[k];

        t_72[k] = -4.0 * gk_72[k]
                  + f_0 * ik_72[k];

        t_73[k] = -4.0 * gk_73[k]
                  + f_0 * ik_73[k];

        t_74[k] = -4.0 * gk_74[k]
                  + f_0 * ik_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, gk_75, gk_76, gk_77, gk_78, gk_79, \
                         ik_75, ik_76, ik_77, ik_78, ik_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -4.0 * gk_75[k]
                  + f_0 * ik_75[k];

        t_76[k] = -4.0 * gk_76[k]
                  + f_0 * ik_76[k];

        t_77[k] = -4.0 * gk_77[k]
                  + f_0 * ik_77[k];

        t_78[k] = -4.0 * gk_78[k]
                  + f_0 * ik_78[k];

        t_79[k] = -4.0 * gk_79[k]
                  + f_0 * ik_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, gk_80, gk_81, gk_82, gk_83, gk_84, \
                         ik_80, ik_81, ik_82, ik_83, ik_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -4.0 * gk_80[k]
                  + f_0 * ik_80[k];

        t_81[k] = -4.0 * gk_81[k]
                  + f_0 * ik_81[k];

        t_82[k] = -4.0 * gk_82[k]
                  + f_0 * ik_82[k];

        t_83[k] = -4.0 * gk_83[k]
                  + f_0 * ik_83[k];

        t_84[k] = -4.0 * gk_84[k]
                  + f_0 * ik_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, gk_85, gk_86, gk_87, gk_88, gk_89, \
                         ik_85, ik_86, ik_87, ik_88, ik_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * gk_85[k]
                  + f_0 * ik_85[k];

        t_86[k] = -4.0 * gk_86[k]
                  + f_0 * ik_86[k];

        t_87[k] = -4.0 * gk_87[k]
                  + f_0 * ik_87[k];

        t_88[k] = -4.0 * gk_88[k]
                  + f_0 * ik_88[k];

        t_89[k] = -4.0 * gk_89[k]
                  + f_0 * ik_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, gk_90, gk_91, gk_92, gk_93, gk_94, \
                         ik_90, ik_91, ik_92, ik_93, ik_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -4.0 * gk_90[k]
                  + f_0 * ik_90[k];

        t_91[k] = -4.0 * gk_91[k]
                  + f_0 * ik_91[k];

        t_92[k] = -4.0 * gk_92[k]
                  + f_0 * ik_92[k];

        t_93[k] = -4.0 * gk_93[k]
                  + f_0 * ik_93[k];

        t_94[k] = -4.0 * gk_94[k]
                  + f_0 * ik_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, gk_95, gk_96, gk_97, gk_98, gk_99, \
                         ik_95, ik_96, ik_97, ik_98, ik_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -4.0 * gk_95[k]
                  + f_0 * ik_95[k];

        t_96[k] = -4.0 * gk_96[k]
                  + f_0 * ik_96[k];

        t_97[k] = -4.0 * gk_97[k]
                  + f_0 * ik_97[k];

        t_98[k] = -4.0 * gk_98[k]
                  + f_0 * ik_98[k];

        t_99[k] = -4.0 * gk_99[k]
                  + f_0 * ik_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, gk_100, gk_101, gk_102, gk_103, \
                         gk_104, ik_100, ik_101, ik_102, ik_103, \
                         ik_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -4.0 * gk_100[k]
                   + f_0 * ik_100[k];

        t_101[k] = -4.0 * gk_101[k]
                   + f_0 * ik_101[k];

        t_102[k] = -4.0 * gk_102[k]
                   + f_0 * ik_102[k];

        t_103[k] = -4.0 * gk_103[k]
                   + f_0 * ik_103[k];

        t_104[k] = -4.0 * gk_104[k]
                   + f_0 * ik_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, gk_105, gk_106, gk_107, gk_108, \
                         gk_109, ik_105, ik_106, ik_107, ik_108, \
                         ik_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -4.0 * gk_105[k]
                   + f_0 * ik_105[k];

        t_106[k] = -4.0 * gk_106[k]
                   + f_0 * ik_106[k];

        t_107[k] = -4.0 * gk_107[k]
                   + f_0 * ik_107[k];

        t_108[k] = -3.0 * gk_108[k]
                   + f_0 * ik_108[k];

        t_109[k] = -3.0 * gk_109[k]
                   + f_0 * ik_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, gk_110, gk_111, gk_112, gk_113, \
                         gk_114, ik_110, ik_111, ik_112, ik_113, \
                         ik_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -3.0 * gk_110[k]
                   + f_0 * ik_110[k];

        t_111[k] = -3.0 * gk_111[k]
                   + f_0 * ik_111[k];

        t_112[k] = -3.0 * gk_112[k]
                   + f_0 * ik_112[k];

        t_113[k] = -3.0 * gk_113[k]
                   + f_0 * ik_113[k];

        t_114[k] = -3.0 * gk_114[k]
                   + f_0 * ik_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, gk_115, gk_116, gk_117, gk_118, \
                         gk_119, ik_115, ik_116, ik_117, ik_118, \
                         ik_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -3.0 * gk_115[k]
                   + f_0 * ik_115[k];

        t_116[k] = -3.0 * gk_116[k]
                   + f_0 * ik_116[k];

        t_117[k] = -3.0 * gk_117[k]
                   + f_0 * ik_117[k];

        t_118[k] = -3.0 * gk_118[k]
                   + f_0 * ik_118[k];

        t_119[k] = -3.0 * gk_119[k]
                   + f_0 * ik_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, gk_120, gk_121, gk_122, gk_123, \
                         gk_124, ik_120, ik_121, ik_122, ik_123, \
                         ik_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -3.0 * gk_120[k]
                   + f_0 * ik_120[k];

        t_121[k] = -3.0 * gk_121[k]
                   + f_0 * ik_121[k];

        t_122[k] = -3.0 * gk_122[k]
                   + f_0 * ik_122[k];

        t_123[k] = -3.0 * gk_123[k]
                   + f_0 * ik_123[k];

        t_124[k] = -3.0 * gk_124[k]
                   + f_0 * ik_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, gk_125, gk_126, gk_127, gk_128, \
                         gk_129, ik_125, ik_126, ik_127, ik_128, \
                         ik_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -3.0 * gk_125[k]
                   + f_0 * ik_125[k];

        t_126[k] = -3.0 * gk_126[k]
                   + f_0 * ik_126[k];

        t_127[k] = -3.0 * gk_127[k]
                   + f_0 * ik_127[k];

        t_128[k] = -3.0 * gk_128[k]
                   + f_0 * ik_128[k];

        t_129[k] = -3.0 * gk_129[k]
                   + f_0 * ik_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, gk_130, gk_131, gk_132, gk_133, \
                         gk_134, ik_130, ik_131, ik_132, ik_133, \
                         ik_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -3.0 * gk_130[k]
                   + f_0 * ik_130[k];

        t_131[k] = -3.0 * gk_131[k]
                   + f_0 * ik_131[k];

        t_132[k] = -3.0 * gk_132[k]
                   + f_0 * ik_132[k];

        t_133[k] = -3.0 * gk_133[k]
                   + f_0 * ik_133[k];

        t_134[k] = -3.0 * gk_134[k]
                   + f_0 * ik_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, gk_135, gk_136, gk_137, gk_138, \
                         gk_139, ik_135, ik_136, ik_137, ik_138, \
                         ik_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -3.0 * gk_135[k]
                   + f_0 * ik_135[k];

        t_136[k] = -3.0 * gk_136[k]
                   + f_0 * ik_136[k];

        t_137[k] = -3.0 * gk_137[k]
                   + f_0 * ik_137[k];

        t_138[k] = -3.0 * gk_138[k]
                   + f_0 * ik_138[k];

        t_139[k] = -3.0 * gk_139[k]
                   + f_0 * ik_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, gk_140, gk_141, gk_142, gk_143, \
                         gk_144, ik_140, ik_141, ik_142, ik_143, \
                         ik_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -3.0 * gk_140[k]
                   + f_0 * ik_140[k];

        t_141[k] = -3.0 * gk_141[k]
                   + f_0 * ik_141[k];

        t_142[k] = -3.0 * gk_142[k]
                   + f_0 * ik_142[k];

        t_143[k] = -3.0 * gk_143[k]
                   + f_0 * ik_143[k];

        t_144[k] = -3.0 * gk_144[k]
                   + f_0 * ik_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, gk_145, gk_146, gk_147, gk_148, \
                         gk_149, ik_145, ik_146, ik_147, ik_148, \
                         ik_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -3.0 * gk_145[k]
                   + f_0 * ik_145[k];

        t_146[k] = -3.0 * gk_146[k]
                   + f_0 * ik_146[k];

        t_147[k] = -3.0 * gk_147[k]
                   + f_0 * ik_147[k];

        t_148[k] = -3.0 * gk_148[k]
                   + f_0 * ik_148[k];

        t_149[k] = -3.0 * gk_149[k]
                   + f_0 * ik_149[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, gk_150, gk_151, gk_152, gk_153, \
                         gk_154, ik_150, ik_151, ik_152, ik_153, \
                         ik_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -3.0 * gk_150[k]
                   + f_0 * ik_150[k];

        t_151[k] = -3.0 * gk_151[k]
                   + f_0 * ik_151[k];

        t_152[k] = -3.0 * gk_152[k]
                   + f_0 * ik_152[k];

        t_153[k] = -3.0 * gk_153[k]
                   + f_0 * ik_153[k];

        t_154[k] = -3.0 * gk_154[k]
                   + f_0 * ik_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, gk_155, gk_156, gk_157, gk_158, \
                         gk_159, ik_155, ik_156, ik_157, ik_158, \
                         ik_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -3.0 * gk_155[k]
                   + f_0 * ik_155[k];

        t_156[k] = -3.0 * gk_156[k]
                   + f_0 * ik_156[k];

        t_157[k] = -3.0 * gk_157[k]
                   + f_0 * ik_157[k];

        t_158[k] = -3.0 * gk_158[k]
                   + f_0 * ik_158[k];

        t_159[k] = -3.0 * gk_159[k]
                   + f_0 * ik_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, gk_160, gk_161, gk_162, gk_163, \
                         gk_164, ik_160, ik_161, ik_162, ik_163, \
                         ik_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -3.0 * gk_160[k]
                   + f_0 * ik_160[k];

        t_161[k] = -3.0 * gk_161[k]
                   + f_0 * ik_161[k];

        t_162[k] = -3.0 * gk_162[k]
                   + f_0 * ik_162[k];

        t_163[k] = -3.0 * gk_163[k]
                   + f_0 * ik_163[k];

        t_164[k] = -3.0 * gk_164[k]
                   + f_0 * ik_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, gk_165, gk_166, gk_167, gk_168, \
                         gk_169, ik_165, ik_166, ik_167, ik_168, \
                         ik_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -3.0 * gk_165[k]
                   + f_0 * ik_165[k];

        t_166[k] = -3.0 * gk_166[k]
                   + f_0 * ik_166[k];

        t_167[k] = -3.0 * gk_167[k]
                   + f_0 * ik_167[k];

        t_168[k] = -3.0 * gk_168[k]
                   + f_0 * ik_168[k];

        t_169[k] = -3.0 * gk_169[k]
                   + f_0 * ik_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, gk_170, gk_171, gk_172, gk_173, \
                         gk_174, ik_170, ik_171, ik_172, ik_173, \
                         ik_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * gk_170[k]
                   + f_0 * ik_170[k];

        t_171[k] = -3.0 * gk_171[k]
                   + f_0 * ik_171[k];

        t_172[k] = -3.0 * gk_172[k]
                   + f_0 * ik_172[k];

        t_173[k] = -3.0 * gk_173[k]
                   + f_0 * ik_173[k];

        t_174[k] = -3.0 * gk_174[k]
                   + f_0 * ik_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, gk_175, gk_176, gk_177, gk_178, \
                         gk_179, ik_175, ik_176, ik_177, ik_178, \
                         ik_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * gk_175[k]
                   + f_0 * ik_175[k];

        t_176[k] = -3.0 * gk_176[k]
                   + f_0 * ik_176[k];

        t_177[k] = -3.0 * gk_177[k]
                   + f_0 * ik_177[k];

        t_178[k] = -3.0 * gk_178[k]
                   + f_0 * ik_178[k];

        t_179[k] = -3.0 * gk_179[k]
                   + f_0 * ik_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, gk_180, gk_181, gk_182, gk_183, \
                         gk_184, ik_180, ik_181, ik_182, ik_183, \
                         ik_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * gk_180[k]
                   + f_0 * ik_180[k];

        t_181[k] = -3.0 * gk_181[k]
                   + f_0 * ik_181[k];

        t_182[k] = -3.0 * gk_182[k]
                   + f_0 * ik_182[k];

        t_183[k] = -3.0 * gk_183[k]
                   + f_0 * ik_183[k];

        t_184[k] = -3.0 * gk_184[k]
                   + f_0 * ik_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, gk_185, gk_186, gk_187, gk_188, \
                         gk_189, ik_185, ik_186, ik_187, ik_188, \
                         ik_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * gk_185[k]
                   + f_0 * ik_185[k];

        t_186[k] = -3.0 * gk_186[k]
                   + f_0 * ik_186[k];

        t_187[k] = -3.0 * gk_187[k]
                   + f_0 * ik_187[k];

        t_188[k] = -3.0 * gk_188[k]
                   + f_0 * ik_188[k];

        t_189[k] = -3.0 * gk_189[k]
                   + f_0 * ik_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, gk_190, gk_191, gk_192, gk_193, \
                         gk_194, ik_190, ik_191, ik_192, ik_193, \
                         ik_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * gk_190[k]
                   + f_0 * ik_190[k];

        t_191[k] = -3.0 * gk_191[k]
                   + f_0 * ik_191[k];

        t_192[k] = -3.0 * gk_192[k]
                   + f_0 * ik_192[k];

        t_193[k] = -3.0 * gk_193[k]
                   + f_0 * ik_193[k];

        t_194[k] = -3.0 * gk_194[k]
                   + f_0 * ik_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, gk_195, gk_196, gk_197, gk_198, \
                         gk_199, ik_195, ik_196, ik_197, ik_198, \
                         ik_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * gk_195[k]
                   + f_0 * ik_195[k];

        t_196[k] = -3.0 * gk_196[k]
                   + f_0 * ik_196[k];

        t_197[k] = -3.0 * gk_197[k]
                   + f_0 * ik_197[k];

        t_198[k] = -3.0 * gk_198[k]
                   + f_0 * ik_198[k];

        t_199[k] = -3.0 * gk_199[k]
                   + f_0 * ik_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, gk_200, gk_201, gk_202, gk_203, \
                         gk_204, ik_200, ik_201, ik_202, ik_203, \
                         ik_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -3.0 * gk_200[k]
                   + f_0 * ik_200[k];

        t_201[k] = -3.0 * gk_201[k]
                   + f_0 * ik_201[k];

        t_202[k] = -3.0 * gk_202[k]
                   + f_0 * ik_202[k];

        t_203[k] = -3.0 * gk_203[k]
                   + f_0 * ik_203[k];

        t_204[k] = -3.0 * gk_204[k]
                   + f_0 * ik_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, gk_205, gk_206, gk_207, gk_208, \
                         gk_209, ik_205, ik_206, ik_207, ik_208, \
                         ik_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -3.0 * gk_205[k]
                   + f_0 * ik_205[k];

        t_206[k] = -3.0 * gk_206[k]
                   + f_0 * ik_206[k];

        t_207[k] = -3.0 * gk_207[k]
                   + f_0 * ik_207[k];

        t_208[k] = -3.0 * gk_208[k]
                   + f_0 * ik_208[k];

        t_209[k] = -3.0 * gk_209[k]
                   + f_0 * ik_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, gk_210, gk_211, gk_212, gk_213, \
                         gk_214, ik_210, ik_211, ik_212, ik_213, \
                         ik_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -3.0 * gk_210[k]
                   + f_0 * ik_210[k];

        t_211[k] = -3.0 * gk_211[k]
                   + f_0 * ik_211[k];

        t_212[k] = -3.0 * gk_212[k]
                   + f_0 * ik_212[k];

        t_213[k] = -3.0 * gk_213[k]
                   + f_0 * ik_213[k];

        t_214[k] = -3.0 * gk_214[k]
                   + f_0 * ik_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, gk_215, gk_216, gk_217, gk_218, \
                         gk_219, ik_215, ik_216, ik_217, ik_218, \
                         ik_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -3.0 * gk_215[k]
                   + f_0 * ik_215[k];

        t_216[k] = -2.0 * gk_216[k]
                   + f_0 * ik_216[k];

        t_217[k] = -2.0 * gk_217[k]
                   + f_0 * ik_217[k];

        t_218[k] = -2.0 * gk_218[k]
                   + f_0 * ik_218[k];

        t_219[k] = -2.0 * gk_219[k]
                   + f_0 * ik_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, gk_220, gk_221, gk_222, gk_223, \
                         gk_224, ik_220, ik_221, ik_222, ik_223, \
                         ik_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * gk_220[k]
                   + f_0 * ik_220[k];

        t_221[k] = -2.0 * gk_221[k]
                   + f_0 * ik_221[k];

        t_222[k] = -2.0 * gk_222[k]
                   + f_0 * ik_222[k];

        t_223[k] = -2.0 * gk_223[k]
                   + f_0 * ik_223[k];

        t_224[k] = -2.0 * gk_224[k]
                   + f_0 * ik_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, gk_225, gk_226, gk_227, gk_228, \
                         gk_229, ik_225, ik_226, ik_227, ik_228, \
                         ik_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -2.0 * gk_225[k]
                   + f_0 * ik_225[k];

        t_226[k] = -2.0 * gk_226[k]
                   + f_0 * ik_226[k];

        t_227[k] = -2.0 * gk_227[k]
                   + f_0 * ik_227[k];

        t_228[k] = -2.0 * gk_228[k]
                   + f_0 * ik_228[k];

        t_229[k] = -2.0 * gk_229[k]
                   + f_0 * ik_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, gk_230, gk_231, gk_232, gk_233, \
                         gk_234, ik_230, ik_231, ik_232, ik_233, \
                         ik_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -2.0 * gk_230[k]
                   + f_0 * ik_230[k];

        t_231[k] = -2.0 * gk_231[k]
                   + f_0 * ik_231[k];

        t_232[k] = -2.0 * gk_232[k]
                   + f_0 * ik_232[k];

        t_233[k] = -2.0 * gk_233[k]
                   + f_0 * ik_233[k];

        t_234[k] = -2.0 * gk_234[k]
                   + f_0 * ik_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, gk_235, gk_236, gk_237, gk_238, \
                         gk_239, ik_235, ik_236, ik_237, ik_238, \
                         ik_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -2.0 * gk_235[k]
                   + f_0 * ik_235[k];

        t_236[k] = -2.0 * gk_236[k]
                   + f_0 * ik_236[k];

        t_237[k] = -2.0 * gk_237[k]
                   + f_0 * ik_237[k];

        t_238[k] = -2.0 * gk_238[k]
                   + f_0 * ik_238[k];

        t_239[k] = -2.0 * gk_239[k]
                   + f_0 * ik_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, gk_240, gk_241, gk_242, gk_243, \
                         gk_244, ik_240, ik_241, ik_242, ik_243, \
                         ik_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -2.0 * gk_240[k]
                   + f_0 * ik_240[k];

        t_241[k] = -2.0 * gk_241[k]
                   + f_0 * ik_241[k];

        t_242[k] = -2.0 * gk_242[k]
                   + f_0 * ik_242[k];

        t_243[k] = -2.0 * gk_243[k]
                   + f_0 * ik_243[k];

        t_244[k] = -2.0 * gk_244[k]
                   + f_0 * ik_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, gk_245, gk_246, gk_247, gk_248, \
                         gk_249, ik_245, ik_246, ik_247, ik_248, \
                         ik_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -2.0 * gk_245[k]
                   + f_0 * ik_245[k];

        t_246[k] = -2.0 * gk_246[k]
                   + f_0 * ik_246[k];

        t_247[k] = -2.0 * gk_247[k]
                   + f_0 * ik_247[k];

        t_248[k] = -2.0 * gk_248[k]
                   + f_0 * ik_248[k];

        t_249[k] = -2.0 * gk_249[k]
                   + f_0 * ik_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, gk_250, gk_251, gk_252, gk_253, \
                         gk_254, ik_250, ik_251, ik_252, ik_253, \
                         ik_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * gk_250[k]
                   + f_0 * ik_250[k];

        t_251[k] = -2.0 * gk_251[k]
                   + f_0 * ik_251[k];

        t_252[k] = -2.0 * gk_252[k]
                   + f_0 * ik_252[k];

        t_253[k] = -2.0 * gk_253[k]
                   + f_0 * ik_253[k];

        t_254[k] = -2.0 * gk_254[k]
                   + f_0 * ik_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, gk_255, gk_256, gk_257, gk_258, \
                         gk_259, ik_255, ik_256, ik_257, ik_258, \
                         ik_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * gk_255[k]
                   + f_0 * ik_255[k];

        t_256[k] = -2.0 * gk_256[k]
                   + f_0 * ik_256[k];

        t_257[k] = -2.0 * gk_257[k]
                   + f_0 * ik_257[k];

        t_258[k] = -2.0 * gk_258[k]
                   + f_0 * ik_258[k];

        t_259[k] = -2.0 * gk_259[k]
                   + f_0 * ik_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, gk_260, gk_261, gk_262, gk_263, \
                         gk_264, ik_260, ik_261, ik_262, ik_263, \
                         ik_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * gk_260[k]
                   + f_0 * ik_260[k];

        t_261[k] = -2.0 * gk_261[k]
                   + f_0 * ik_261[k];

        t_262[k] = -2.0 * gk_262[k]
                   + f_0 * ik_262[k];

        t_263[k] = -2.0 * gk_263[k]
                   + f_0 * ik_263[k];

        t_264[k] = -2.0 * gk_264[k]
                   + f_0 * ik_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, gk_265, gk_266, gk_267, gk_268, \
                         gk_269, ik_265, ik_266, ik_267, ik_268, \
                         ik_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * gk_265[k]
                   + f_0 * ik_265[k];

        t_266[k] = -2.0 * gk_266[k]
                   + f_0 * ik_266[k];

        t_267[k] = -2.0 * gk_267[k]
                   + f_0 * ik_267[k];

        t_268[k] = -2.0 * gk_268[k]
                   + f_0 * ik_268[k];

        t_269[k] = -2.0 * gk_269[k]
                   + f_0 * ik_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, gk_270, gk_271, gk_272, gk_273, \
                         gk_274, ik_270, ik_271, ik_272, ik_273, \
                         ik_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * gk_270[k]
                   + f_0 * ik_270[k];

        t_271[k] = -2.0 * gk_271[k]
                   + f_0 * ik_271[k];

        t_272[k] = -2.0 * gk_272[k]
                   + f_0 * ik_272[k];

        t_273[k] = -2.0 * gk_273[k]
                   + f_0 * ik_273[k];

        t_274[k] = -2.0 * gk_274[k]
                   + f_0 * ik_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, gk_275, gk_276, gk_277, gk_278, \
                         gk_279, ik_275, ik_276, ik_277, ik_278, \
                         ik_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * gk_275[k]
                   + f_0 * ik_275[k];

        t_276[k] = -2.0 * gk_276[k]
                   + f_0 * ik_276[k];

        t_277[k] = -2.0 * gk_277[k]
                   + f_0 * ik_277[k];

        t_278[k] = -2.0 * gk_278[k]
                   + f_0 * ik_278[k];

        t_279[k] = -2.0 * gk_279[k]
                   + f_0 * ik_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, gk_280, gk_281, gk_282, gk_283, \
                         gk_284, ik_280, ik_281, ik_282, ik_283, \
                         ik_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * gk_280[k]
                   + f_0 * ik_280[k];

        t_281[k] = -2.0 * gk_281[k]
                   + f_0 * ik_281[k];

        t_282[k] = -2.0 * gk_282[k]
                   + f_0 * ik_282[k];

        t_283[k] = -2.0 * gk_283[k]
                   + f_0 * ik_283[k];

        t_284[k] = -2.0 * gk_284[k]
                   + f_0 * ik_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, gk_285, gk_286, gk_287, gk_288, \
                         gk_289, ik_285, ik_286, ik_287, ik_288, \
                         ik_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -2.0 * gk_285[k]
                   + f_0 * ik_285[k];

        t_286[k] = -2.0 * gk_286[k]
                   + f_0 * ik_286[k];

        t_287[k] = -2.0 * gk_287[k]
                   + f_0 * ik_287[k];

        t_288[k] = -2.0 * gk_288[k]
                   + f_0 * ik_288[k];

        t_289[k] = -2.0 * gk_289[k]
                   + f_0 * ik_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, gk_290, gk_291, gk_292, gk_293, \
                         gk_294, ik_290, ik_291, ik_292, ik_293, \
                         ik_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -2.0 * gk_290[k]
                   + f_0 * ik_290[k];

        t_291[k] = -2.0 * gk_291[k]
                   + f_0 * ik_291[k];

        t_292[k] = -2.0 * gk_292[k]
                   + f_0 * ik_292[k];

        t_293[k] = -2.0 * gk_293[k]
                   + f_0 * ik_293[k];

        t_294[k] = -2.0 * gk_294[k]
                   + f_0 * ik_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, gk_295, gk_296, gk_297, gk_298, \
                         gk_299, ik_295, ik_296, ik_297, ik_298, \
                         ik_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -2.0 * gk_295[k]
                   + f_0 * ik_295[k];

        t_296[k] = -2.0 * gk_296[k]
                   + f_0 * ik_296[k];

        t_297[k] = -2.0 * gk_297[k]
                   + f_0 * ik_297[k];

        t_298[k] = -2.0 * gk_298[k]
                   + f_0 * ik_298[k];

        t_299[k] = -2.0 * gk_299[k]
                   + f_0 * ik_299[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, gk_300, gk_301, gk_302, gk_303, \
                         gk_304, ik_300, ik_301, ik_302, ik_303, \
                         ik_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -2.0 * gk_300[k]
                   + f_0 * ik_300[k];

        t_301[k] = -2.0 * gk_301[k]
                   + f_0 * ik_301[k];

        t_302[k] = -2.0 * gk_302[k]
                   + f_0 * ik_302[k];

        t_303[k] = -2.0 * gk_303[k]
                   + f_0 * ik_303[k];

        t_304[k] = -2.0 * gk_304[k]
                   + f_0 * ik_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, gk_305, gk_306, gk_307, gk_308, \
                         gk_309, ik_305, ik_306, ik_307, ik_308, \
                         ik_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -2.0 * gk_305[k]
                   + f_0 * ik_305[k];

        t_306[k] = -2.0 * gk_306[k]
                   + f_0 * ik_306[k];

        t_307[k] = -2.0 * gk_307[k]
                   + f_0 * ik_307[k];

        t_308[k] = -2.0 * gk_308[k]
                   + f_0 * ik_308[k];

        t_309[k] = -2.0 * gk_309[k]
                   + f_0 * ik_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, gk_310, gk_311, gk_312, gk_313, \
                         gk_314, ik_310, ik_311, ik_312, ik_313, \
                         ik_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -2.0 * gk_310[k]
                   + f_0 * ik_310[k];

        t_311[k] = -2.0 * gk_311[k]
                   + f_0 * ik_311[k];

        t_312[k] = -2.0 * gk_312[k]
                   + f_0 * ik_312[k];

        t_313[k] = -2.0 * gk_313[k]
                   + f_0 * ik_313[k];

        t_314[k] = -2.0 * gk_314[k]
                   + f_0 * ik_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, gk_315, gk_316, gk_317, gk_318, \
                         gk_319, ik_315, ik_316, ik_317, ik_318, \
                         ik_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -2.0 * gk_315[k]
                   + f_0 * ik_315[k];

        t_316[k] = -2.0 * gk_316[k]
                   + f_0 * ik_316[k];

        t_317[k] = -2.0 * gk_317[k]
                   + f_0 * ik_317[k];

        t_318[k] = -2.0 * gk_318[k]
                   + f_0 * ik_318[k];

        t_319[k] = -2.0 * gk_319[k]
                   + f_0 * ik_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, gk_320, gk_321, gk_322, gk_323, \
                         gk_324, ik_320, ik_321, ik_322, ik_323, \
                         ik_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -2.0 * gk_320[k]
                   + f_0 * ik_320[k];

        t_321[k] = -2.0 * gk_321[k]
                   + f_0 * ik_321[k];

        t_322[k] = -2.0 * gk_322[k]
                   + f_0 * ik_322[k];

        t_323[k] = -2.0 * gk_323[k]
                   + f_0 * ik_323[k];

        t_324[k] = -2.0 * gk_324[k]
                   + f_0 * ik_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, gk_325, gk_326, gk_327, gk_328, \
                         gk_329, ik_325, ik_326, ik_327, ik_328, \
                         ik_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -2.0 * gk_325[k]
                   + f_0 * ik_325[k];

        t_326[k] = -2.0 * gk_326[k]
                   + f_0 * ik_326[k];

        t_327[k] = -2.0 * gk_327[k]
                   + f_0 * ik_327[k];

        t_328[k] = -2.0 * gk_328[k]
                   + f_0 * ik_328[k];

        t_329[k] = -2.0 * gk_329[k]
                   + f_0 * ik_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, gk_330, gk_331, gk_332, gk_333, \
                         gk_334, ik_330, ik_331, ik_332, ik_333, \
                         ik_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -2.0 * gk_330[k]
                   + f_0 * ik_330[k];

        t_331[k] = -2.0 * gk_331[k]
                   + f_0 * ik_331[k];

        t_332[k] = -2.0 * gk_332[k]
                   + f_0 * ik_332[k];

        t_333[k] = -2.0 * gk_333[k]
                   + f_0 * ik_333[k];

        t_334[k] = -2.0 * gk_334[k]
                   + f_0 * ik_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, gk_335, gk_336, gk_337, gk_338, \
                         gk_339, ik_335, ik_336, ik_337, ik_338, \
                         ik_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -2.0 * gk_335[k]
                   + f_0 * ik_335[k];

        t_336[k] = -2.0 * gk_336[k]
                   + f_0 * ik_336[k];

        t_337[k] = -2.0 * gk_337[k]
                   + f_0 * ik_337[k];

        t_338[k] = -2.0 * gk_338[k]
                   + f_0 * ik_338[k];

        t_339[k] = -2.0 * gk_339[k]
                   + f_0 * ik_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, gk_340, gk_341, gk_342, gk_343, \
                         gk_344, ik_340, ik_341, ik_342, ik_343, \
                         ik_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * gk_340[k]
                   + f_0 * ik_340[k];

        t_341[k] = -2.0 * gk_341[k]
                   + f_0 * ik_341[k];

        t_342[k] = -2.0 * gk_342[k]
                   + f_0 * ik_342[k];

        t_343[k] = -2.0 * gk_343[k]
                   + f_0 * ik_343[k];

        t_344[k] = -2.0 * gk_344[k]
                   + f_0 * ik_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, gk_345, gk_346, gk_347, gk_348, \
                         gk_349, ik_345, ik_346, ik_347, ik_348, \
                         ik_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * gk_345[k]
                   + f_0 * ik_345[k];

        t_346[k] = -2.0 * gk_346[k]
                   + f_0 * ik_346[k];

        t_347[k] = -2.0 * gk_347[k]
                   + f_0 * ik_347[k];

        t_348[k] = -2.0 * gk_348[k]
                   + f_0 * ik_348[k];

        t_349[k] = -2.0 * gk_349[k]
                   + f_0 * ik_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, gk_350, gk_351, gk_352, gk_353, \
                         gk_354, ik_350, ik_351, ik_352, ik_353, \
                         ik_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * gk_350[k]
                   + f_0 * ik_350[k];

        t_351[k] = -2.0 * gk_351[k]
                   + f_0 * ik_351[k];

        t_352[k] = -2.0 * gk_352[k]
                   + f_0 * ik_352[k];

        t_353[k] = -2.0 * gk_353[k]
                   + f_0 * ik_353[k];

        t_354[k] = -2.0 * gk_354[k]
                   + f_0 * ik_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, gk_355, gk_356, gk_357, gk_358, \
                         gk_359, ik_355, ik_356, ik_357, ik_358, \
                         ik_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * gk_355[k]
                   + f_0 * ik_355[k];

        t_356[k] = -2.0 * gk_356[k]
                   + f_0 * ik_356[k];

        t_357[k] = -2.0 * gk_357[k]
                   + f_0 * ik_357[k];

        t_358[k] = -2.0 * gk_358[k]
                   + f_0 * ik_358[k];

        t_359[k] = -2.0 * gk_359[k]
                   + f_0 * ik_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, gk_360, gk_361, gk_362, gk_363, \
                         gk_364, ik_360, ik_361, ik_362, ik_363, \
                         ik_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -gk_360[k]
                   + f_0 * ik_360[k];

        t_361[k] = -gk_361[k]
                   + f_0 * ik_361[k];

        t_362[k] = -gk_362[k]
                   + f_0 * ik_362[k];

        t_363[k] = -gk_363[k]
                   + f_0 * ik_363[k];

        t_364[k] = -gk_364[k]
                   + f_0 * ik_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, gk_365, gk_366, gk_367, gk_368, \
                         gk_369, ik_365, ik_366, ik_367, ik_368, \
                         ik_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -gk_365[k]
                   + f_0 * ik_365[k];

        t_366[k] = -gk_366[k]
                   + f_0 * ik_366[k];

        t_367[k] = -gk_367[k]
                   + f_0 * ik_367[k];

        t_368[k] = -gk_368[k]
                   + f_0 * ik_368[k];

        t_369[k] = -gk_369[k]
                   + f_0 * ik_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, gk_370, gk_371, gk_372, gk_373, \
                         gk_374, ik_370, ik_371, ik_372, ik_373, \
                         ik_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -gk_370[k]
                   + f_0 * ik_370[k];

        t_371[k] = -gk_371[k]
                   + f_0 * ik_371[k];

        t_372[k] = -gk_372[k]
                   + f_0 * ik_372[k];

        t_373[k] = -gk_373[k]
                   + f_0 * ik_373[k];

        t_374[k] = -gk_374[k]
                   + f_0 * ik_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, gk_375, gk_376, gk_377, gk_378, \
                         gk_379, ik_375, ik_376, ik_377, ik_378, \
                         ik_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -gk_375[k]
                   + f_0 * ik_375[k];

        t_376[k] = -gk_376[k]
                   + f_0 * ik_376[k];

        t_377[k] = -gk_377[k]
                   + f_0 * ik_377[k];

        t_378[k] = -gk_378[k]
                   + f_0 * ik_378[k];

        t_379[k] = -gk_379[k]
                   + f_0 * ik_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, gk_380, gk_381, gk_382, gk_383, \
                         gk_384, ik_380, ik_381, ik_382, ik_383, \
                         ik_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -gk_380[k]
                   + f_0 * ik_380[k];

        t_381[k] = -gk_381[k]
                   + f_0 * ik_381[k];

        t_382[k] = -gk_382[k]
                   + f_0 * ik_382[k];

        t_383[k] = -gk_383[k]
                   + f_0 * ik_383[k];

        t_384[k] = -gk_384[k]
                   + f_0 * ik_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, gk_385, gk_386, gk_387, gk_388, \
                         gk_389, ik_385, ik_386, ik_387, ik_388, \
                         ik_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -gk_385[k]
                   + f_0 * ik_385[k];

        t_386[k] = -gk_386[k]
                   + f_0 * ik_386[k];

        t_387[k] = -gk_387[k]
                   + f_0 * ik_387[k];

        t_388[k] = -gk_388[k]
                   + f_0 * ik_388[k];

        t_389[k] = -gk_389[k]
                   + f_0 * ik_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, gk_390, gk_391, gk_392, gk_393, \
                         gk_394, ik_390, ik_391, ik_392, ik_393, \
                         ik_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -gk_390[k]
                   + f_0 * ik_390[k];

        t_391[k] = -gk_391[k]
                   + f_0 * ik_391[k];

        t_392[k] = -gk_392[k]
                   + f_0 * ik_392[k];

        t_393[k] = -gk_393[k]
                   + f_0 * ik_393[k];

        t_394[k] = -gk_394[k]
                   + f_0 * ik_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, gk_395, gk_396, gk_397, gk_398, \
                         gk_399, ik_395, ik_396, ik_397, ik_398, \
                         ik_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -gk_395[k]
                   + f_0 * ik_395[k];

        t_396[k] = -gk_396[k]
                   + f_0 * ik_396[k];

        t_397[k] = -gk_397[k]
                   + f_0 * ik_397[k];

        t_398[k] = -gk_398[k]
                   + f_0 * ik_398[k];

        t_399[k] = -gk_399[k]
                   + f_0 * ik_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, gk_400, gk_401, gk_402, gk_403, \
                         gk_404, ik_400, ik_401, ik_402, ik_403, \
                         ik_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -gk_400[k]
                   + f_0 * ik_400[k];

        t_401[k] = -gk_401[k]
                   + f_0 * ik_401[k];

        t_402[k] = -gk_402[k]
                   + f_0 * ik_402[k];

        t_403[k] = -gk_403[k]
                   + f_0 * ik_403[k];

        t_404[k] = -gk_404[k]
                   + f_0 * ik_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, gk_405, gk_406, gk_407, gk_408, \
                         gk_409, ik_405, ik_406, ik_407, ik_408, \
                         ik_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -gk_405[k]
                   + f_0 * ik_405[k];

        t_406[k] = -gk_406[k]
                   + f_0 * ik_406[k];

        t_407[k] = -gk_407[k]
                   + f_0 * ik_407[k];

        t_408[k] = -gk_408[k]
                   + f_0 * ik_408[k];

        t_409[k] = -gk_409[k]
                   + f_0 * ik_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, gk_410, gk_411, gk_412, gk_413, \
                         gk_414, ik_410, ik_411, ik_412, ik_413, \
                         ik_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -gk_410[k]
                   + f_0 * ik_410[k];

        t_411[k] = -gk_411[k]
                   + f_0 * ik_411[k];

        t_412[k] = -gk_412[k]
                   + f_0 * ik_412[k];

        t_413[k] = -gk_413[k]
                   + f_0 * ik_413[k];

        t_414[k] = -gk_414[k]
                   + f_0 * ik_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, gk_415, gk_416, gk_417, gk_418, \
                         gk_419, ik_415, ik_416, ik_417, ik_418, \
                         ik_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -gk_415[k]
                   + f_0 * ik_415[k];

        t_416[k] = -gk_416[k]
                   + f_0 * ik_416[k];

        t_417[k] = -gk_417[k]
                   + f_0 * ik_417[k];

        t_418[k] = -gk_418[k]
                   + f_0 * ik_418[k];

        t_419[k] = -gk_419[k]
                   + f_0 * ik_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, gk_420, gk_421, gk_422, gk_423, \
                         gk_424, ik_420, ik_421, ik_422, ik_423, \
                         ik_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -gk_420[k]
                   + f_0 * ik_420[k];

        t_421[k] = -gk_421[k]
                   + f_0 * ik_421[k];

        t_422[k] = -gk_422[k]
                   + f_0 * ik_422[k];

        t_423[k] = -gk_423[k]
                   + f_0 * ik_423[k];

        t_424[k] = -gk_424[k]
                   + f_0 * ik_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, gk_425, gk_426, gk_427, gk_428, \
                         gk_429, ik_425, ik_426, ik_427, ik_428, \
                         ik_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -gk_425[k]
                   + f_0 * ik_425[k];

        t_426[k] = -gk_426[k]
                   + f_0 * ik_426[k];

        t_427[k] = -gk_427[k]
                   + f_0 * ik_427[k];

        t_428[k] = -gk_428[k]
                   + f_0 * ik_428[k];

        t_429[k] = -gk_429[k]
                   + f_0 * ik_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, gk_430, gk_431, gk_432, gk_433, \
                         gk_434, ik_430, ik_431, ik_432, ik_433, \
                         ik_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -gk_430[k]
                   + f_0 * ik_430[k];

        t_431[k] = -gk_431[k]
                   + f_0 * ik_431[k];

        t_432[k] = -gk_432[k]
                   + f_0 * ik_432[k];

        t_433[k] = -gk_433[k]
                   + f_0 * ik_433[k];

        t_434[k] = -gk_434[k]
                   + f_0 * ik_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, gk_435, gk_436, gk_437, gk_438, \
                         gk_439, ik_435, ik_436, ik_437, ik_438, \
                         ik_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -gk_435[k]
                   + f_0 * ik_435[k];

        t_436[k] = -gk_436[k]
                   + f_0 * ik_436[k];

        t_437[k] = -gk_437[k]
                   + f_0 * ik_437[k];

        t_438[k] = -gk_438[k]
                   + f_0 * ik_438[k];

        t_439[k] = -gk_439[k]
                   + f_0 * ik_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, gk_440, gk_441, gk_442, gk_443, \
                         gk_444, ik_440, ik_441, ik_442, ik_443, \
                         ik_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -gk_440[k]
                   + f_0 * ik_440[k];

        t_441[k] = -gk_441[k]
                   + f_0 * ik_441[k];

        t_442[k] = -gk_442[k]
                   + f_0 * ik_442[k];

        t_443[k] = -gk_443[k]
                   + f_0 * ik_443[k];

        t_444[k] = -gk_444[k]
                   + f_0 * ik_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, gk_445, gk_446, gk_447, gk_448, \
                         gk_449, ik_445, ik_446, ik_447, ik_448, \
                         ik_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -gk_445[k]
                   + f_0 * ik_445[k];

        t_446[k] = -gk_446[k]
                   + f_0 * ik_446[k];

        t_447[k] = -gk_447[k]
                   + f_0 * ik_447[k];

        t_448[k] = -gk_448[k]
                   + f_0 * ik_448[k];

        t_449[k] = -gk_449[k]
                   + f_0 * ik_449[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, gk_450, gk_451, gk_452, gk_453, \
                         gk_454, ik_450, ik_451, ik_452, ik_453, \
                         ik_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -gk_450[k]
                   + f_0 * ik_450[k];

        t_451[k] = -gk_451[k]
                   + f_0 * ik_451[k];

        t_452[k] = -gk_452[k]
                   + f_0 * ik_452[k];

        t_453[k] = -gk_453[k]
                   + f_0 * ik_453[k];

        t_454[k] = -gk_454[k]
                   + f_0 * ik_454[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, gk_455, gk_456, gk_457, gk_458, \
                         gk_459, ik_455, ik_456, ik_457, ik_458, \
                         ik_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -gk_455[k]
                   + f_0 * ik_455[k];

        t_456[k] = -gk_456[k]
                   + f_0 * ik_456[k];

        t_457[k] = -gk_457[k]
                   + f_0 * ik_457[k];

        t_458[k] = -gk_458[k]
                   + f_0 * ik_458[k];

        t_459[k] = -gk_459[k]
                   + f_0 * ik_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, gk_460, gk_461, gk_462, gk_463, \
                         gk_464, ik_460, ik_461, ik_462, ik_463, \
                         ik_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -gk_460[k]
                   + f_0 * ik_460[k];

        t_461[k] = -gk_461[k]
                   + f_0 * ik_461[k];

        t_462[k] = -gk_462[k]
                   + f_0 * ik_462[k];

        t_463[k] = -gk_463[k]
                   + f_0 * ik_463[k];

        t_464[k] = -gk_464[k]
                   + f_0 * ik_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, gk_465, gk_466, gk_467, gk_468, \
                         gk_469, ik_465, ik_466, ik_467, ik_468, \
                         ik_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -gk_465[k]
                   + f_0 * ik_465[k];

        t_466[k] = -gk_466[k]
                   + f_0 * ik_466[k];

        t_467[k] = -gk_467[k]
                   + f_0 * ik_467[k];

        t_468[k] = -gk_468[k]
                   + f_0 * ik_468[k];

        t_469[k] = -gk_469[k]
                   + f_0 * ik_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, gk_470, gk_471, gk_472, gk_473, \
                         gk_474, ik_470, ik_471, ik_472, ik_473, \
                         ik_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -gk_470[k]
                   + f_0 * ik_470[k];

        t_471[k] = -gk_471[k]
                   + f_0 * ik_471[k];

        t_472[k] = -gk_472[k]
                   + f_0 * ik_472[k];

        t_473[k] = -gk_473[k]
                   + f_0 * ik_473[k];

        t_474[k] = -gk_474[k]
                   + f_0 * ik_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, gk_475, gk_476, gk_477, gk_478, \
                         gk_479, ik_475, ik_476, ik_477, ik_478, \
                         ik_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -gk_475[k]
                   + f_0 * ik_475[k];

        t_476[k] = -gk_476[k]
                   + f_0 * ik_476[k];

        t_477[k] = -gk_477[k]
                   + f_0 * ik_477[k];

        t_478[k] = -gk_478[k]
                   + f_0 * ik_478[k];

        t_479[k] = -gk_479[k]
                   + f_0 * ik_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, gk_480, gk_481, gk_482, gk_483, \
                         gk_484, ik_480, ik_481, ik_482, ik_483, \
                         ik_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -gk_480[k]
                   + f_0 * ik_480[k];

        t_481[k] = -gk_481[k]
                   + f_0 * ik_481[k];

        t_482[k] = -gk_482[k]
                   + f_0 * ik_482[k];

        t_483[k] = -gk_483[k]
                   + f_0 * ik_483[k];

        t_484[k] = -gk_484[k]
                   + f_0 * ik_484[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, gk_485, gk_486, gk_487, gk_488, \
                         gk_489, ik_485, ik_486, ik_487, ik_488, \
                         ik_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -gk_485[k]
                   + f_0 * ik_485[k];

        t_486[k] = -gk_486[k]
                   + f_0 * ik_486[k];

        t_487[k] = -gk_487[k]
                   + f_0 * ik_487[k];

        t_488[k] = -gk_488[k]
                   + f_0 * ik_488[k];

        t_489[k] = -gk_489[k]
                   + f_0 * ik_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, gk_490, gk_491, gk_492, gk_493, \
                         gk_494, ik_490, ik_491, ik_492, ik_493, \
                         ik_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -gk_490[k]
                   + f_0 * ik_490[k];

        t_491[k] = -gk_491[k]
                   + f_0 * ik_491[k];

        t_492[k] = -gk_492[k]
                   + f_0 * ik_492[k];

        t_493[k] = -gk_493[k]
                   + f_0 * ik_493[k];

        t_494[k] = -gk_494[k]
                   + f_0 * ik_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, gk_495, gk_496, gk_497, gk_498, \
                         gk_499, ik_495, ik_496, ik_497, ik_498, \
                         ik_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -gk_495[k]
                   + f_0 * ik_495[k];

        t_496[k] = -gk_496[k]
                   + f_0 * ik_496[k];

        t_497[k] = -gk_497[k]
                   + f_0 * ik_497[k];

        t_498[k] = -gk_498[k]
                   + f_0 * ik_498[k];

        t_499[k] = -gk_499[k]
                   + f_0 * ik_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, gk_500, gk_501, gk_502, gk_503, \
                         gk_504, ik_500, ik_501, ik_502, ik_503, \
                         ik_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -gk_500[k]
                   + f_0 * ik_500[k];

        t_501[k] = -gk_501[k]
                   + f_0 * ik_501[k];

        t_502[k] = -gk_502[k]
                   + f_0 * ik_502[k];

        t_503[k] = -gk_503[k]
                   + f_0 * ik_503[k];

        t_504[k] = -gk_504[k]
                   + f_0 * ik_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, gk_505, gk_506, gk_507, gk_508, \
                         gk_509, ik_505, ik_506, ik_507, ik_508, \
                         ik_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -gk_505[k]
                   + f_0 * ik_505[k];

        t_506[k] = -gk_506[k]
                   + f_0 * ik_506[k];

        t_507[k] = -gk_507[k]
                   + f_0 * ik_507[k];

        t_508[k] = -gk_508[k]
                   + f_0 * ik_508[k];

        t_509[k] = -gk_509[k]
                   + f_0 * ik_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, gk_510, gk_511, gk_512, gk_513, \
                         gk_514, ik_510, ik_511, ik_512, ik_513, \
                         ik_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -gk_510[k]
                   + f_0 * ik_510[k];

        t_511[k] = -gk_511[k]
                   + f_0 * ik_511[k];

        t_512[k] = -gk_512[k]
                   + f_0 * ik_512[k];

        t_513[k] = -gk_513[k]
                   + f_0 * ik_513[k];

        t_514[k] = -gk_514[k]
                   + f_0 * ik_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, gk_515, gk_516, gk_517, gk_518, \
                         gk_519, ik_515, ik_516, ik_517, ik_518, \
                         ik_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -gk_515[k]
                   + f_0 * ik_515[k];

        t_516[k] = -gk_516[k]
                   + f_0 * ik_516[k];

        t_517[k] = -gk_517[k]
                   + f_0 * ik_517[k];

        t_518[k] = -gk_518[k]
                   + f_0 * ik_518[k];

        t_519[k] = -gk_519[k]
                   + f_0 * ik_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, gk_520, gk_521, gk_522, gk_523, \
                         gk_524, ik_520, ik_521, ik_522, ik_523, \
                         ik_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -gk_520[k]
                   + f_0 * ik_520[k];

        t_521[k] = -gk_521[k]
                   + f_0 * ik_521[k];

        t_522[k] = -gk_522[k]
                   + f_0 * ik_522[k];

        t_523[k] = -gk_523[k]
                   + f_0 * ik_523[k];

        t_524[k] = -gk_524[k]
                   + f_0 * ik_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, gk_525, gk_526, gk_527, gk_528, \
                         gk_529, ik_525, ik_526, ik_527, ik_528, \
                         ik_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -gk_525[k]
                   + f_0 * ik_525[k];

        t_526[k] = -gk_526[k]
                   + f_0 * ik_526[k];

        t_527[k] = -gk_527[k]
                   + f_0 * ik_527[k];

        t_528[k] = -gk_528[k]
                   + f_0 * ik_528[k];

        t_529[k] = -gk_529[k]
                   + f_0 * ik_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, gk_530, gk_531, gk_532, gk_533, \
                         gk_534, ik_530, ik_531, ik_532, ik_533, \
                         ik_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -gk_530[k]
                   + f_0 * ik_530[k];

        t_531[k] = -gk_531[k]
                   + f_0 * ik_531[k];

        t_532[k] = -gk_532[k]
                   + f_0 * ik_532[k];

        t_533[k] = -gk_533[k]
                   + f_0 * ik_533[k];

        t_534[k] = -gk_534[k]
                   + f_0 * ik_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, gk_535, gk_536, gk_537, gk_538, \
                         gk_539, ik_535, ik_536, ik_537, ik_538, \
                         ik_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -gk_535[k]
                   + f_0 * ik_535[k];

        t_536[k] = -gk_536[k]
                   + f_0 * ik_536[k];

        t_537[k] = -gk_537[k]
                   + f_0 * ik_537[k];

        t_538[k] = -gk_538[k]
                   + f_0 * ik_538[k];

        t_539[k] = -gk_539[k]
                   + f_0 * ik_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, t_545, t_546, t_547, ik_540, \
                         ik_541, ik_542, ik_543, ik_544, ik_545, ik_546, \
                         ik_547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_0 * ik_540[k];

        t_541[k] = f_0 * ik_541[k];

        t_542[k] = f_0 * ik_542[k];

        t_543[k] = f_0 * ik_543[k];

        t_544[k] = f_0 * ik_544[k];

        t_545[k] = f_0 * ik_545[k];

        t_546[k] = f_0 * ik_546[k];

        t_547[k] = f_0 * ik_547[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, t_553, t_554, t_555, ik_548, \
                         ik_549, ik_550, ik_551, ik_552, ik_553, ik_554, \
                         ik_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_0 * ik_548[k];

        t_549[k] = f_0 * ik_549[k];

        t_550[k] = f_0 * ik_550[k];

        t_551[k] = f_0 * ik_551[k];

        t_552[k] = f_0 * ik_552[k];

        t_553[k] = f_0 * ik_553[k];

        t_554[k] = f_0 * ik_554[k];

        t_555[k] = f_0 * ik_555[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, t_561, t_562, t_563, ik_556, \
                         ik_557, ik_558, ik_559, ik_560, ik_561, ik_562, \
                         ik_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_0 * ik_556[k];

        t_557[k] = f_0 * ik_557[k];

        t_558[k] = f_0 * ik_558[k];

        t_559[k] = f_0 * ik_559[k];

        t_560[k] = f_0 * ik_560[k];

        t_561[k] = f_0 * ik_561[k];

        t_562[k] = f_0 * ik_562[k];

        t_563[k] = f_0 * ik_563[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, t_568, t_569, t_570, t_571, ik_564, \
                         ik_565, ik_566, ik_567, ik_568, ik_569, ik_570, \
                         ik_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_0 * ik_564[k];

        t_565[k] = f_0 * ik_565[k];

        t_566[k] = f_0 * ik_566[k];

        t_567[k] = f_0 * ik_567[k];

        t_568[k] = f_0 * ik_568[k];

        t_569[k] = f_0 * ik_569[k];

        t_570[k] = f_0 * ik_570[k];

        t_571[k] = f_0 * ik_571[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, t_577, t_578, t_579, ik_572, \
                         ik_573, ik_574, ik_575, ik_576, ik_577, ik_578, \
                         ik_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_0 * ik_572[k];

        t_573[k] = f_0 * ik_573[k];

        t_574[k] = f_0 * ik_574[k];

        t_575[k] = f_0 * ik_575[k];

        t_576[k] = f_0 * ik_576[k];

        t_577[k] = f_0 * ik_577[k];

        t_578[k] = f_0 * ik_578[k];

        t_579[k] = f_0 * ik_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, t_585, t_586, t_587, ik_580, \
                         ik_581, ik_582, ik_583, ik_584, ik_585, ik_586, \
                         ik_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_0 * ik_580[k];

        t_581[k] = f_0 * ik_581[k];

        t_582[k] = f_0 * ik_582[k];

        t_583[k] = f_0 * ik_583[k];

        t_584[k] = f_0 * ik_584[k];

        t_585[k] = f_0 * ik_585[k];

        t_586[k] = f_0 * ik_586[k];

        t_587[k] = f_0 * ik_587[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, t_592, t_593, t_594, t_595, ik_588, \
                         ik_589, ik_590, ik_591, ik_592, ik_593, ik_594, \
                         ik_595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_0 * ik_588[k];

        t_589[k] = f_0 * ik_589[k];

        t_590[k] = f_0 * ik_590[k];

        t_591[k] = f_0 * ik_591[k];

        t_592[k] = f_0 * ik_592[k];

        t_593[k] = f_0 * ik_593[k];

        t_594[k] = f_0 * ik_594[k];

        t_595[k] = f_0 * ik_595[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, t_600, t_601, t_602, t_603, ik_596, \
                         ik_597, ik_598, ik_599, ik_600, ik_601, ik_602, \
                         ik_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_0 * ik_596[k];

        t_597[k] = f_0 * ik_597[k];

        t_598[k] = f_0 * ik_598[k];

        t_599[k] = f_0 * ik_599[k];

        t_600[k] = f_0 * ik_600[k];

        t_601[k] = f_0 * ik_601[k];

        t_602[k] = f_0 * ik_602[k];

        t_603[k] = f_0 * ik_603[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, t_609, t_610, t_611, ik_604, \
                         ik_605, ik_606, ik_607, ik_608, ik_609, ik_610, \
                         ik_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_0 * ik_604[k];

        t_605[k] = f_0 * ik_605[k];

        t_606[k] = f_0 * ik_606[k];

        t_607[k] = f_0 * ik_607[k];

        t_608[k] = f_0 * ik_608[k];

        t_609[k] = f_0 * ik_609[k];

        t_610[k] = f_0 * ik_610[k];

        t_611[k] = f_0 * ik_611[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, t_617, t_618, t_619, ik_612, \
                         ik_613, ik_614, ik_615, ik_616, ik_617, ik_618, \
                         ik_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_0 * ik_612[k];

        t_613[k] = f_0 * ik_613[k];

        t_614[k] = f_0 * ik_614[k];

        t_615[k] = f_0 * ik_615[k];

        t_616[k] = f_0 * ik_616[k];

        t_617[k] = f_0 * ik_617[k];

        t_618[k] = f_0 * ik_618[k];

        t_619[k] = f_0 * ik_619[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, t_625, t_626, t_627, ik_620, \
                         ik_621, ik_622, ik_623, ik_624, ik_625, ik_626, \
                         ik_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_0 * ik_620[k];

        t_621[k] = f_0 * ik_621[k];

        t_622[k] = f_0 * ik_622[k];

        t_623[k] = f_0 * ik_623[k];

        t_624[k] = f_0 * ik_624[k];

        t_625[k] = f_0 * ik_625[k];

        t_626[k] = f_0 * ik_626[k];

        t_627[k] = f_0 * ik_627[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ik, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    const auto *ik_754 = buffer.data(ik + 754);
    const auto *ik_755 = buffer.data(ik + 755);

#pragma omp simd aligned(t_628, t_629, t_630, t_631, t_632, t_633, t_634, t_635, ik_628, \
                         ik_629, ik_630, ik_631, ik_632, ik_633, ik_634, \
                         ik_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_0 * ik_628[k];

        t_629[k] = f_0 * ik_629[k];

        t_630[k] = f_0 * ik_630[k];

        t_631[k] = f_0 * ik_631[k];

        t_632[k] = f_0 * ik_632[k];

        t_633[k] = f_0 * ik_633[k];

        t_634[k] = f_0 * ik_634[k];

        t_635[k] = f_0 * ik_635[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, t_640, t_641, t_642, t_643, ik_636, \
                         ik_637, ik_638, ik_639, ik_640, ik_641, ik_642, \
                         ik_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_0 * ik_636[k];

        t_637[k] = f_0 * ik_637[k];

        t_638[k] = f_0 * ik_638[k];

        t_639[k] = f_0 * ik_639[k];

        t_640[k] = f_0 * ik_640[k];

        t_641[k] = f_0 * ik_641[k];

        t_642[k] = f_0 * ik_642[k];

        t_643[k] = f_0 * ik_643[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, t_648, t_649, t_650, t_651, ik_644, \
                         ik_645, ik_646, ik_647, ik_648, ik_649, ik_650, \
                         ik_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_0 * ik_644[k];

        t_645[k] = f_0 * ik_645[k];

        t_646[k] = f_0 * ik_646[k];

        t_647[k] = f_0 * ik_647[k];

        t_648[k] = f_0 * ik_648[k];

        t_649[k] = f_0 * ik_649[k];

        t_650[k] = f_0 * ik_650[k];

        t_651[k] = f_0 * ik_651[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, t_656, t_657, t_658, t_659, ik_652, \
                         ik_653, ik_654, ik_655, ik_656, ik_657, ik_658, \
                         ik_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_0 * ik_652[k];

        t_653[k] = f_0 * ik_653[k];

        t_654[k] = f_0 * ik_654[k];

        t_655[k] = f_0 * ik_655[k];

        t_656[k] = f_0 * ik_656[k];

        t_657[k] = f_0 * ik_657[k];

        t_658[k] = f_0 * ik_658[k];

        t_659[k] = f_0 * ik_659[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, t_665, t_666, t_667, ik_660, \
                         ik_661, ik_662, ik_663, ik_664, ik_665, ik_666, \
                         ik_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_0 * ik_660[k];

        t_661[k] = f_0 * ik_661[k];

        t_662[k] = f_0 * ik_662[k];

        t_663[k] = f_0 * ik_663[k];

        t_664[k] = f_0 * ik_664[k];

        t_665[k] = f_0 * ik_665[k];

        t_666[k] = f_0 * ik_666[k];

        t_667[k] = f_0 * ik_667[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, t_672, t_673, t_674, t_675, ik_668, \
                         ik_669, ik_670, ik_671, ik_672, ik_673, ik_674, \
                         ik_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_0 * ik_668[k];

        t_669[k] = f_0 * ik_669[k];

        t_670[k] = f_0 * ik_670[k];

        t_671[k] = f_0 * ik_671[k];

        t_672[k] = f_0 * ik_672[k];

        t_673[k] = f_0 * ik_673[k];

        t_674[k] = f_0 * ik_674[k];

        t_675[k] = f_0 * ik_675[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, t_680, t_681, t_682, t_683, ik_676, \
                         ik_677, ik_678, ik_679, ik_680, ik_681, ik_682, \
                         ik_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_0 * ik_676[k];

        t_677[k] = f_0 * ik_677[k];

        t_678[k] = f_0 * ik_678[k];

        t_679[k] = f_0 * ik_679[k];

        t_680[k] = f_0 * ik_680[k];

        t_681[k] = f_0 * ik_681[k];

        t_682[k] = f_0 * ik_682[k];

        t_683[k] = f_0 * ik_683[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, t_688, t_689, t_690, t_691, ik_684, \
                         ik_685, ik_686, ik_687, ik_688, ik_689, ik_690, \
                         ik_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_0 * ik_684[k];

        t_685[k] = f_0 * ik_685[k];

        t_686[k] = f_0 * ik_686[k];

        t_687[k] = f_0 * ik_687[k];

        t_688[k] = f_0 * ik_688[k];

        t_689[k] = f_0 * ik_689[k];

        t_690[k] = f_0 * ik_690[k];

        t_691[k] = f_0 * ik_691[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, t_696, t_697, t_698, t_699, ik_692, \
                         ik_693, ik_694, ik_695, ik_696, ik_697, ik_698, \
                         ik_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = f_0 * ik_692[k];

        t_693[k] = f_0 * ik_693[k];

        t_694[k] = f_0 * ik_694[k];

        t_695[k] = f_0 * ik_695[k];

        t_696[k] = f_0 * ik_696[k];

        t_697[k] = f_0 * ik_697[k];

        t_698[k] = f_0 * ik_698[k];

        t_699[k] = f_0 * ik_699[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, t_705, t_706, t_707, ik_700, \
                         ik_701, ik_702, ik_703, ik_704, ik_705, ik_706, \
                         ik_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_0 * ik_700[k];

        t_701[k] = f_0 * ik_701[k];

        t_702[k] = f_0 * ik_702[k];

        t_703[k] = f_0 * ik_703[k];

        t_704[k] = f_0 * ik_704[k];

        t_705[k] = f_0 * ik_705[k];

        t_706[k] = f_0 * ik_706[k];

        t_707[k] = f_0 * ik_707[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, t_711, t_712, t_713, t_714, t_715, ik_708, \
                         ik_709, ik_710, ik_711, ik_712, ik_713, ik_714, \
                         ik_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_0 * ik_708[k];

        t_709[k] = f_0 * ik_709[k];

        t_710[k] = f_0 * ik_710[k];

        t_711[k] = f_0 * ik_711[k];

        t_712[k] = f_0 * ik_712[k];

        t_713[k] = f_0 * ik_713[k];

        t_714[k] = f_0 * ik_714[k];

        t_715[k] = f_0 * ik_715[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, t_720, t_721, t_722, t_723, ik_716, \
                         ik_717, ik_718, ik_719, ik_720, ik_721, ik_722, \
                         ik_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_0 * ik_716[k];

        t_717[k] = f_0 * ik_717[k];

        t_718[k] = f_0 * ik_718[k];

        t_719[k] = f_0 * ik_719[k];

        t_720[k] = f_0 * ik_720[k];

        t_721[k] = f_0 * ik_721[k];

        t_722[k] = f_0 * ik_722[k];

        t_723[k] = f_0 * ik_723[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, t_729, t_730, t_731, ik_724, \
                         ik_725, ik_726, ik_727, ik_728, ik_729, ik_730, \
                         ik_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_0 * ik_724[k];

        t_725[k] = f_0 * ik_725[k];

        t_726[k] = f_0 * ik_726[k];

        t_727[k] = f_0 * ik_727[k];

        t_728[k] = f_0 * ik_728[k];

        t_729[k] = f_0 * ik_729[k];

        t_730[k] = f_0 * ik_730[k];

        t_731[k] = f_0 * ik_731[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, t_735, t_736, t_737, t_738, t_739, ik_732, \
                         ik_733, ik_734, ik_735, ik_736, ik_737, ik_738, \
                         ik_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_0 * ik_732[k];

        t_733[k] = f_0 * ik_733[k];

        t_734[k] = f_0 * ik_734[k];

        t_735[k] = f_0 * ik_735[k];

        t_736[k] = f_0 * ik_736[k];

        t_737[k] = f_0 * ik_737[k];

        t_738[k] = f_0 * ik_738[k];

        t_739[k] = f_0 * ik_739[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, t_745, t_746, t_747, ik_740, \
                         ik_741, ik_742, ik_743, ik_744, ik_745, ik_746, \
                         ik_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_0 * ik_740[k];

        t_741[k] = f_0 * ik_741[k];

        t_742[k] = f_0 * ik_742[k];

        t_743[k] = f_0 * ik_743[k];

        t_744[k] = f_0 * ik_744[k];

        t_745[k] = f_0 * ik_745[k];

        t_746[k] = f_0 * ik_746[k];

        t_747[k] = f_0 * ik_747[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, t_752, t_753, t_754, t_755, ik_748, \
                         ik_749, ik_750, ik_751, ik_752, ik_753, ik_754, \
                         ik_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_0 * ik_748[k];

        t_749[k] = f_0 * ik_749[k];

        t_750[k] = f_0 * ik_750[k];

        t_751[k] = f_0 * ik_751[k];

        t_752[k] = f_0 * ik_752[k];

        t_753[k] = f_0 * ik_753[k];

        t_754[k] = f_0 * ik_754[k];

        t_755[k] = f_0 * ik_755[k];
    }
}

auto
compute_prim_geom_10_hk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t gk, const size_t ik,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hk_electron_repulsion_0_piece0(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_0_piece1(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_0_piece2(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_0_piece3(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_0_piece4(buffer, target, ik, ncols, alpha);
}

static auto
compute_prim_geom_10_hk_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ik_36, ik_37, ik_38, ik_39, \
                         ik_40, ik_41, ik_42, ik_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ik_36[k];

        t_1[k] = f_0 * ik_37[k];

        t_2[k] = f_0 * ik_38[k];

        t_3[k] = f_0 * ik_39[k];

        t_4[k] = f_0 * ik_40[k];

        t_5[k] = f_0 * ik_41[k];

        t_6[k] = f_0 * ik_42[k];

        t_7[k] = f_0 * ik_43[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, ik_44, ik_45, ik_46, \
                         ik_47, ik_48, ik_49, ik_50, ik_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ik_44[k];

        t_9[k] = f_0 * ik_45[k];

        t_10[k] = f_0 * ik_46[k];

        t_11[k] = f_0 * ik_47[k];

        t_12[k] = f_0 * ik_48[k];

        t_13[k] = f_0 * ik_49[k];

        t_14[k] = f_0 * ik_50[k];

        t_15[k] = f_0 * ik_51[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, ik_52, ik_53, ik_54, \
                         ik_55, ik_56, ik_57, ik_58, ik_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * ik_52[k];

        t_17[k] = f_0 * ik_53[k];

        t_18[k] = f_0 * ik_54[k];

        t_19[k] = f_0 * ik_55[k];

        t_20[k] = f_0 * ik_56[k];

        t_21[k] = f_0 * ik_57[k];

        t_22[k] = f_0 * ik_58[k];

        t_23[k] = f_0 * ik_59[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, ik_60, ik_61, ik_62, \
                         ik_63, ik_64, ik_65, ik_66, ik_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * ik_60[k];

        t_25[k] = f_0 * ik_61[k];

        t_26[k] = f_0 * ik_62[k];

        t_27[k] = f_0 * ik_63[k];

        t_28[k] = f_0 * ik_64[k];

        t_29[k] = f_0 * ik_65[k];

        t_30[k] = f_0 * ik_66[k];

        t_31[k] = f_0 * ik_67[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, gk_0, gk_1, ik_68, ik_69, ik_70, \
                         ik_71, ik_108, ik_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ik_68[k];

        t_33[k] = f_0 * ik_69[k];

        t_34[k] = f_0 * ik_70[k];

        t_35[k] = f_0 * ik_71[k];

        t_36[k] = -gk_0[k]
                  + f_0 * ik_108[k];

        t_37[k] = -gk_1[k]
                  + f_0 * ik_109[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, gk_2, gk_3, gk_4, gk_5, gk_6, ik_110, \
                         ik_111, ik_112, ik_113, ik_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -gk_2[k]
                  + f_0 * ik_110[k];

        t_39[k] = -gk_3[k]
                  + f_0 * ik_111[k];

        t_40[k] = -gk_4[k]
                  + f_0 * ik_112[k];

        t_41[k] = -gk_5[k]
                  + f_0 * ik_113[k];

        t_42[k] = -gk_6[k]
                  + f_0 * ik_114[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, gk_7, gk_8, gk_9, gk_10, gk_11, ik_115, \
                         ik_116, ik_117, ik_118, ik_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -gk_7[k]
                  + f_0 * ik_115[k];

        t_44[k] = -gk_8[k]
                  + f_0 * ik_116[k];

        t_45[k] = -gk_9[k]
                  + f_0 * ik_117[k];

        t_46[k] = -gk_10[k]
                  + f_0 * ik_118[k];

        t_47[k] = -gk_11[k]
                  + f_0 * ik_119[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, gk_12, gk_13, gk_14, gk_15, gk_16, \
                         ik_120, ik_121, ik_122, ik_123, ik_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -gk_12[k]
                  + f_0 * ik_120[k];

        t_49[k] = -gk_13[k]
                  + f_0 * ik_121[k];

        t_50[k] = -gk_14[k]
                  + f_0 * ik_122[k];

        t_51[k] = -gk_15[k]
                  + f_0 * ik_123[k];

        t_52[k] = -gk_16[k]
                  + f_0 * ik_124[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, gk_17, gk_18, gk_19, gk_20, gk_21, \
                         ik_125, ik_126, ik_127, ik_128, ik_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -gk_17[k]
                  + f_0 * ik_125[k];

        t_54[k] = -gk_18[k]
                  + f_0 * ik_126[k];

        t_55[k] = -gk_19[k]
                  + f_0 * ik_127[k];

        t_56[k] = -gk_20[k]
                  + f_0 * ik_128[k];

        t_57[k] = -gk_21[k]
                  + f_0 * ik_129[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, gk_22, gk_23, gk_24, gk_25, gk_26, \
                         ik_130, ik_131, ik_132, ik_133, ik_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -gk_22[k]
                  + f_0 * ik_130[k];

        t_59[k] = -gk_23[k]
                  + f_0 * ik_131[k];

        t_60[k] = -gk_24[k]
                  + f_0 * ik_132[k];

        t_61[k] = -gk_25[k]
                  + f_0 * ik_133[k];

        t_62[k] = -gk_26[k]
                  + f_0 * ik_134[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, gk_27, gk_28, gk_29, gk_30, gk_31, \
                         ik_135, ik_136, ik_137, ik_138, ik_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -gk_27[k]
                  + f_0 * ik_135[k];

        t_64[k] = -gk_28[k]
                  + f_0 * ik_136[k];

        t_65[k] = -gk_29[k]
                  + f_0 * ik_137[k];

        t_66[k] = -gk_30[k]
                  + f_0 * ik_138[k];

        t_67[k] = -gk_31[k]
                  + f_0 * ik_139[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, gk_32, gk_33, gk_34, gk_35, \
                         ik_140, ik_141, ik_142, ik_143, ik_144, \
                         ik_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -gk_32[k]
                  + f_0 * ik_140[k];

        t_69[k] = -gk_33[k]
                  + f_0 * ik_141[k];

        t_70[k] = -gk_34[k]
                  + f_0 * ik_142[k];

        t_71[k] = -gk_35[k]
                  + f_0 * ik_143[k];

        t_72[k] = f_0 * ik_144[k];

        t_73[k] = f_0 * ik_145[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, t_80, t_81, ik_146, ik_147, \
                         ik_148, ik_149, ik_150, ik_151, ik_152, \
                         ik_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * ik_146[k];

        t_75[k] = f_0 * ik_147[k];

        t_76[k] = f_0 * ik_148[k];

        t_77[k] = f_0 * ik_149[k];

        t_78[k] = f_0 * ik_150[k];

        t_79[k] = f_0 * ik_151[k];

        t_80[k] = f_0 * ik_152[k];

        t_81[k] = f_0 * ik_153[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, t_87, t_88, t_89, ik_154, ik_155, \
                         ik_156, ik_157, ik_158, ik_159, ik_160, \
                         ik_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * ik_154[k];

        t_83[k] = f_0 * ik_155[k];

        t_84[k] = f_0 * ik_156[k];

        t_85[k] = f_0 * ik_157[k];

        t_86[k] = f_0 * ik_158[k];

        t_87[k] = f_0 * ik_159[k];

        t_88[k] = f_0 * ik_160[k];

        t_89[k] = f_0 * ik_161[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, ik_162, ik_163, \
                         ik_164, ik_165, ik_166, ik_167, ik_168, \
                         ik_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * ik_162[k];

        t_91[k] = f_0 * ik_163[k];

        t_92[k] = f_0 * ik_164[k];

        t_93[k] = f_0 * ik_165[k];

        t_94[k] = f_0 * ik_166[k];

        t_95[k] = f_0 * ik_167[k];

        t_96[k] = f_0 * ik_168[k];

        t_97[k] = f_0 * ik_169[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, t_104, t_105, ik_170, ik_171, \
                         ik_172, ik_173, ik_174, ik_175, ik_176, \
                         ik_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * ik_170[k];

        t_99[k] = f_0 * ik_171[k];

        t_100[k] = f_0 * ik_172[k];

        t_101[k] = f_0 * ik_173[k];

        t_102[k] = f_0 * ik_174[k];

        t_103[k] = f_0 * ik_175[k];

        t_104[k] = f_0 * ik_176[k];

        t_105[k] = f_0 * ik_177[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, gk_36, gk_37, gk_38, gk_39, \
                         ik_178, ik_179, ik_216, ik_217, ik_218, \
                         ik_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * ik_178[k];

        t_107[k] = f_0 * ik_179[k];

        t_108[k] = -2.0 * gk_36[k]
                   + f_0 * ik_216[k];

        t_109[k] = -2.0 * gk_37[k]
                   + f_0 * ik_217[k];

        t_110[k] = -2.0 * gk_38[k]
                   + f_0 * ik_218[k];

        t_111[k] = -2.0 * gk_39[k]
                   + f_0 * ik_219[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, gk_40, gk_41, gk_42, gk_43, gk_44, \
                         ik_220, ik_221, ik_222, ik_223, ik_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -2.0 * gk_40[k]
                   + f_0 * ik_220[k];

        t_113[k] = -2.0 * gk_41[k]
                   + f_0 * ik_221[k];

        t_114[k] = -2.0 * gk_42[k]
                   + f_0 * ik_222[k];

        t_115[k] = -2.0 * gk_43[k]
                   + f_0 * ik_223[k];

        t_116[k] = -2.0 * gk_44[k]
                   + f_0 * ik_224[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, gk_45, gk_46, gk_47, gk_48, gk_49, \
                         ik_225, ik_226, ik_227, ik_228, ik_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -2.0 * gk_45[k]
                   + f_0 * ik_225[k];

        t_118[k] = -2.0 * gk_46[k]
                   + f_0 * ik_226[k];

        t_119[k] = -2.0 * gk_47[k]
                   + f_0 * ik_227[k];

        t_120[k] = -2.0 * gk_48[k]
                   + f_0 * ik_228[k];

        t_121[k] = -2.0 * gk_49[k]
                   + f_0 * ik_229[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, gk_50, gk_51, gk_52, gk_53, gk_54, \
                         ik_230, ik_231, ik_232, ik_233, ik_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * gk_50[k]
                   + f_0 * ik_230[k];

        t_123[k] = -2.0 * gk_51[k]
                   + f_0 * ik_231[k];

        t_124[k] = -2.0 * gk_52[k]
                   + f_0 * ik_232[k];

        t_125[k] = -2.0 * gk_53[k]
                   + f_0 * ik_233[k];

        t_126[k] = -2.0 * gk_54[k]
                   + f_0 * ik_234[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, gk_55, gk_56, gk_57, gk_58, gk_59, \
                         ik_235, ik_236, ik_237, ik_238, ik_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * gk_55[k]
                   + f_0 * ik_235[k];

        t_128[k] = -2.0 * gk_56[k]
                   + f_0 * ik_236[k];

        t_129[k] = -2.0 * gk_57[k]
                   + f_0 * ik_237[k];

        t_130[k] = -2.0 * gk_58[k]
                   + f_0 * ik_238[k];

        t_131[k] = -2.0 * gk_59[k]
                   + f_0 * ik_239[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, gk_60, gk_61, gk_62, gk_63, gk_64, \
                         ik_240, ik_241, ik_242, ik_243, ik_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * gk_60[k]
                   + f_0 * ik_240[k];

        t_133[k] = -2.0 * gk_61[k]
                   + f_0 * ik_241[k];

        t_134[k] = -2.0 * gk_62[k]
                   + f_0 * ik_242[k];

        t_135[k] = -2.0 * gk_63[k]
                   + f_0 * ik_243[k];

        t_136[k] = -2.0 * gk_64[k]
                   + f_0 * ik_244[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, gk_65, gk_66, gk_67, gk_68, gk_69, \
                         ik_245, ik_246, ik_247, ik_248, ik_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -2.0 * gk_65[k]
                   + f_0 * ik_245[k];

        t_138[k] = -2.0 * gk_66[k]
                   + f_0 * ik_246[k];

        t_139[k] = -2.0 * gk_67[k]
                   + f_0 * ik_247[k];

        t_140[k] = -2.0 * gk_68[k]
                   + f_0 * ik_248[k];

        t_141[k] = -2.0 * gk_69[k]
                   + f_0 * ik_249[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, gk_70, gk_71, gk_72, gk_73, gk_74, \
                         ik_250, ik_251, ik_252, ik_253, ik_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -2.0 * gk_70[k]
                   + f_0 * ik_250[k];

        t_143[k] = -2.0 * gk_71[k]
                   + f_0 * ik_251[k];

        t_144[k] = -gk_72[k]
                   + f_0 * ik_252[k];

        t_145[k] = -gk_73[k]
                   + f_0 * ik_253[k];

        t_146[k] = -gk_74[k]
                   + f_0 * ik_254[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, gk_75, gk_76, gk_77, gk_78, gk_79, \
                         ik_255, ik_256, ik_257, ik_258, ik_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -gk_75[k]
                   + f_0 * ik_255[k];

        t_148[k] = -gk_76[k]
                   + f_0 * ik_256[k];

        t_149[k] = -gk_77[k]
                   + f_0 * ik_257[k];

        t_150[k] = -gk_78[k]
                   + f_0 * ik_258[k];

        t_151[k] = -gk_79[k]
                   + f_0 * ik_259[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, gk_80, gk_81, gk_82, gk_83, gk_84, \
                         ik_260, ik_261, ik_262, ik_263, ik_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -gk_80[k]
                   + f_0 * ik_260[k];

        t_153[k] = -gk_81[k]
                   + f_0 * ik_261[k];

        t_154[k] = -gk_82[k]
                   + f_0 * ik_262[k];

        t_155[k] = -gk_83[k]
                   + f_0 * ik_263[k];

        t_156[k] = -gk_84[k]
                   + f_0 * ik_264[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, gk_85, gk_86, gk_87, gk_88, gk_89, \
                         ik_265, ik_266, ik_267, ik_268, ik_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -gk_85[k]
                   + f_0 * ik_265[k];

        t_158[k] = -gk_86[k]
                   + f_0 * ik_266[k];

        t_159[k] = -gk_87[k]
                   + f_0 * ik_267[k];

        t_160[k] = -gk_88[k]
                   + f_0 * ik_268[k];

        t_161[k] = -gk_89[k]
                   + f_0 * ik_269[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, gk_90, gk_91, gk_92, gk_93, gk_94, \
                         ik_270, ik_271, ik_272, ik_273, ik_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -gk_90[k]
                   + f_0 * ik_270[k];

        t_163[k] = -gk_91[k]
                   + f_0 * ik_271[k];

        t_164[k] = -gk_92[k]
                   + f_0 * ik_272[k];

        t_165[k] = -gk_93[k]
                   + f_0 * ik_273[k];

        t_166[k] = -gk_94[k]
                   + f_0 * ik_274[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, gk_95, gk_96, gk_97, gk_98, gk_99, \
                         ik_275, ik_276, ik_277, ik_278, ik_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -gk_95[k]
                   + f_0 * ik_275[k];

        t_168[k] = -gk_96[k]
                   + f_0 * ik_276[k];

        t_169[k] = -gk_97[k]
                   + f_0 * ik_277[k];

        t_170[k] = -gk_98[k]
                   + f_0 * ik_278[k];

        t_171[k] = -gk_99[k]
                   + f_0 * ik_279[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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
    const auto *ik_477 = buffer.data(ik + 477);

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, gk_100, gk_101, gk_102, gk_103, \
                         gk_104, ik_280, ik_281, ik_282, ik_283, \
                         ik_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -gk_100[k]
                   + f_0 * ik_280[k];

        t_173[k] = -gk_101[k]
                   + f_0 * ik_281[k];

        t_174[k] = -gk_102[k]
                   + f_0 * ik_282[k];

        t_175[k] = -gk_103[k]
                   + f_0 * ik_283[k];

        t_176[k] = -gk_104[k]
                   + f_0 * ik_284[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, t_182, gk_105, gk_106, gk_107, \
                         ik_285, ik_286, ik_287, ik_288, ik_289, \
                         ik_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -gk_105[k]
                   + f_0 * ik_285[k];

        t_178[k] = -gk_106[k]
                   + f_0 * ik_286[k];

        t_179[k] = -gk_107[k]
                   + f_0 * ik_287[k];

        t_180[k] = f_0 * ik_288[k];

        t_181[k] = f_0 * ik_289[k];

        t_182[k] = f_0 * ik_290[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, t_188, t_189, t_190, ik_291, \
                         ik_292, ik_293, ik_294, ik_295, ik_296, ik_297, \
                         ik_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_0 * ik_291[k];

        t_184[k] = f_0 * ik_292[k];

        t_185[k] = f_0 * ik_293[k];

        t_186[k] = f_0 * ik_294[k];

        t_187[k] = f_0 * ik_295[k];

        t_188[k] = f_0 * ik_296[k];

        t_189[k] = f_0 * ik_297[k];

        t_190[k] = f_0 * ik_298[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, t_197, t_198, ik_299, \
                         ik_300, ik_301, ik_302, ik_303, ik_304, ik_305, \
                         ik_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_0 * ik_299[k];

        t_192[k] = f_0 * ik_300[k];

        t_193[k] = f_0 * ik_301[k];

        t_194[k] = f_0 * ik_302[k];

        t_195[k] = f_0 * ik_303[k];

        t_196[k] = f_0 * ik_304[k];

        t_197[k] = f_0 * ik_305[k];

        t_198[k] = f_0 * ik_306[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, t_204, t_205, t_206, ik_307, \
                         ik_308, ik_309, ik_310, ik_311, ik_312, ik_313, \
                         ik_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_0 * ik_307[k];

        t_200[k] = f_0 * ik_308[k];

        t_201[k] = f_0 * ik_309[k];

        t_202[k] = f_0 * ik_310[k];

        t_203[k] = f_0 * ik_311[k];

        t_204[k] = f_0 * ik_312[k];

        t_205[k] = f_0 * ik_313[k];

        t_206[k] = f_0 * ik_314[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, t_213, t_214, ik_315, \
                         ik_316, ik_317, ik_318, ik_319, ik_320, ik_321, \
                         ik_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_0 * ik_315[k];

        t_208[k] = f_0 * ik_316[k];

        t_209[k] = f_0 * ik_317[k];

        t_210[k] = f_0 * ik_318[k];

        t_211[k] = f_0 * ik_319[k];

        t_212[k] = f_0 * ik_320[k];

        t_213[k] = f_0 * ik_321[k];

        t_214[k] = f_0 * ik_322[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, gk_108, gk_109, gk_110, gk_111, \
                         ik_323, ik_360, ik_361, ik_362, ik_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_0 * ik_323[k];

        t_216[k] = -3.0 * gk_108[k]
                   + f_0 * ik_360[k];

        t_217[k] = -3.0 * gk_109[k]
                   + f_0 * ik_361[k];

        t_218[k] = -3.0 * gk_110[k]
                   + f_0 * ik_362[k];

        t_219[k] = -3.0 * gk_111[k]
                   + f_0 * ik_363[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, gk_112, gk_113, gk_114, gk_115, \
                         gk_116, ik_364, ik_365, ik_366, ik_367, \
                         ik_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -3.0 * gk_112[k]
                   + f_0 * ik_364[k];

        t_221[k] = -3.0 * gk_113[k]
                   + f_0 * ik_365[k];

        t_222[k] = -3.0 * gk_114[k]
                   + f_0 * ik_366[k];

        t_223[k] = -3.0 * gk_115[k]
                   + f_0 * ik_367[k];

        t_224[k] = -3.0 * gk_116[k]
                   + f_0 * ik_368[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, gk_117, gk_118, gk_119, gk_120, \
                         gk_121, ik_369, ik_370, ik_371, ik_372, \
                         ik_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -3.0 * gk_117[k]
                   + f_0 * ik_369[k];

        t_226[k] = -3.0 * gk_118[k]
                   + f_0 * ik_370[k];

        t_227[k] = -3.0 * gk_119[k]
                   + f_0 * ik_371[k];

        t_228[k] = -3.0 * gk_120[k]
                   + f_0 * ik_372[k];

        t_229[k] = -3.0 * gk_121[k]
                   + f_0 * ik_373[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, gk_122, gk_123, gk_124, gk_125, \
                         gk_126, ik_374, ik_375, ik_376, ik_377, \
                         ik_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -3.0 * gk_122[k]
                   + f_0 * ik_374[k];

        t_231[k] = -3.0 * gk_123[k]
                   + f_0 * ik_375[k];

        t_232[k] = -3.0 * gk_124[k]
                   + f_0 * ik_376[k];

        t_233[k] = -3.0 * gk_125[k]
                   + f_0 * ik_377[k];

        t_234[k] = -3.0 * gk_126[k]
                   + f_0 * ik_378[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, gk_127, gk_128, gk_129, gk_130, \
                         gk_131, ik_379, ik_380, ik_381, ik_382, \
                         ik_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -3.0 * gk_127[k]
                   + f_0 * ik_379[k];

        t_236[k] = -3.0 * gk_128[k]
                   + f_0 * ik_380[k];

        t_237[k] = -3.0 * gk_129[k]
                   + f_0 * ik_381[k];

        t_238[k] = -3.0 * gk_130[k]
                   + f_0 * ik_382[k];

        t_239[k] = -3.0 * gk_131[k]
                   + f_0 * ik_383[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, gk_132, gk_133, gk_134, gk_135, \
                         gk_136, ik_384, ik_385, ik_386, ik_387, \
                         ik_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * gk_132[k]
                   + f_0 * ik_384[k];

        t_241[k] = -3.0 * gk_133[k]
                   + f_0 * ik_385[k];

        t_242[k] = -3.0 * gk_134[k]
                   + f_0 * ik_386[k];

        t_243[k] = -3.0 * gk_135[k]
                   + f_0 * ik_387[k];

        t_244[k] = -3.0 * gk_136[k]
                   + f_0 * ik_388[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, gk_137, gk_138, gk_139, gk_140, \
                         gk_141, ik_389, ik_390, ik_391, ik_392, \
                         ik_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * gk_137[k]
                   + f_0 * ik_389[k];

        t_246[k] = -3.0 * gk_138[k]
                   + f_0 * ik_390[k];

        t_247[k] = -3.0 * gk_139[k]
                   + f_0 * ik_391[k];

        t_248[k] = -3.0 * gk_140[k]
                   + f_0 * ik_392[k];

        t_249[k] = -3.0 * gk_141[k]
                   + f_0 * ik_393[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, gk_142, gk_143, gk_144, gk_145, \
                         gk_146, ik_394, ik_395, ik_396, ik_397, \
                         ik_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -3.0 * gk_142[k]
                   + f_0 * ik_394[k];

        t_251[k] = -3.0 * gk_143[k]
                   + f_0 * ik_395[k];

        t_252[k] = -2.0 * gk_144[k]
                   + f_0 * ik_396[k];

        t_253[k] = -2.0 * gk_145[k]
                   + f_0 * ik_397[k];

        t_254[k] = -2.0 * gk_146[k]
                   + f_0 * ik_398[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, gk_147, gk_148, gk_149, gk_150, \
                         gk_151, ik_399, ik_400, ik_401, ik_402, \
                         ik_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * gk_147[k]
                   + f_0 * ik_399[k];

        t_256[k] = -2.0 * gk_148[k]
                   + f_0 * ik_400[k];

        t_257[k] = -2.0 * gk_149[k]
                   + f_0 * ik_401[k];

        t_258[k] = -2.0 * gk_150[k]
                   + f_0 * ik_402[k];

        t_259[k] = -2.0 * gk_151[k]
                   + f_0 * ik_403[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, gk_152, gk_153, gk_154, gk_155, \
                         gk_156, ik_404, ik_405, ik_406, ik_407, \
                         ik_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * gk_152[k]
                   + f_0 * ik_404[k];

        t_261[k] = -2.0 * gk_153[k]
                   + f_0 * ik_405[k];

        t_262[k] = -2.0 * gk_154[k]
                   + f_0 * ik_406[k];

        t_263[k] = -2.0 * gk_155[k]
                   + f_0 * ik_407[k];

        t_264[k] = -2.0 * gk_156[k]
                   + f_0 * ik_408[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, gk_157, gk_158, gk_159, gk_160, \
                         gk_161, ik_409, ik_410, ik_411, ik_412, \
                         ik_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * gk_157[k]
                   + f_0 * ik_409[k];

        t_266[k] = -2.0 * gk_158[k]
                   + f_0 * ik_410[k];

        t_267[k] = -2.0 * gk_159[k]
                   + f_0 * ik_411[k];

        t_268[k] = -2.0 * gk_160[k]
                   + f_0 * ik_412[k];

        t_269[k] = -2.0 * gk_161[k]
                   + f_0 * ik_413[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, gk_162, gk_163, gk_164, gk_165, \
                         gk_166, ik_414, ik_415, ik_416, ik_417, \
                         ik_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * gk_162[k]
                   + f_0 * ik_414[k];

        t_271[k] = -2.0 * gk_163[k]
                   + f_0 * ik_415[k];

        t_272[k] = -2.0 * gk_164[k]
                   + f_0 * ik_416[k];

        t_273[k] = -2.0 * gk_165[k]
                   + f_0 * ik_417[k];

        t_274[k] = -2.0 * gk_166[k]
                   + f_0 * ik_418[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, gk_167, gk_168, gk_169, gk_170, \
                         gk_171, ik_419, ik_420, ik_421, ik_422, \
                         ik_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * gk_167[k]
                   + f_0 * ik_419[k];

        t_276[k] = -2.0 * gk_168[k]
                   + f_0 * ik_420[k];

        t_277[k] = -2.0 * gk_169[k]
                   + f_0 * ik_421[k];

        t_278[k] = -2.0 * gk_170[k]
                   + f_0 * ik_422[k];

        t_279[k] = -2.0 * gk_171[k]
                   + f_0 * ik_423[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, gk_172, gk_173, gk_174, gk_175, \
                         gk_176, ik_424, ik_425, ik_426, ik_427, \
                         ik_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * gk_172[k]
                   + f_0 * ik_424[k];

        t_281[k] = -2.0 * gk_173[k]
                   + f_0 * ik_425[k];

        t_282[k] = -2.0 * gk_174[k]
                   + f_0 * ik_426[k];

        t_283[k] = -2.0 * gk_175[k]
                   + f_0 * ik_427[k];

        t_284[k] = -2.0 * gk_176[k]
                   + f_0 * ik_428[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, gk_177, gk_178, gk_179, gk_180, \
                         gk_181, ik_429, ik_430, ik_431, ik_432, \
                         ik_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -2.0 * gk_177[k]
                   + f_0 * ik_429[k];

        t_286[k] = -2.0 * gk_178[k]
                   + f_0 * ik_430[k];

        t_287[k] = -2.0 * gk_179[k]
                   + f_0 * ik_431[k];

        t_288[k] = -gk_180[k]
                   + f_0 * ik_432[k];

        t_289[k] = -gk_181[k]
                   + f_0 * ik_433[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, gk_182, gk_183, gk_184, gk_185, \
                         gk_186, ik_434, ik_435, ik_436, ik_437, \
                         ik_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -gk_182[k]
                   + f_0 * ik_434[k];

        t_291[k] = -gk_183[k]
                   + f_0 * ik_435[k];

        t_292[k] = -gk_184[k]
                   + f_0 * ik_436[k];

        t_293[k] = -gk_185[k]
                   + f_0 * ik_437[k];

        t_294[k] = -gk_186[k]
                   + f_0 * ik_438[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, gk_187, gk_188, gk_189, gk_190, \
                         gk_191, ik_439, ik_440, ik_441, ik_442, \
                         ik_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -gk_187[k]
                   + f_0 * ik_439[k];

        t_296[k] = -gk_188[k]
                   + f_0 * ik_440[k];

        t_297[k] = -gk_189[k]
                   + f_0 * ik_441[k];

        t_298[k] = -gk_190[k]
                   + f_0 * ik_442[k];

        t_299[k] = -gk_191[k]
                   + f_0 * ik_443[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, gk_192, gk_193, gk_194, gk_195, \
                         gk_196, ik_444, ik_445, ik_446, ik_447, \
                         ik_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -gk_192[k]
                   + f_0 * ik_444[k];

        t_301[k] = -gk_193[k]
                   + f_0 * ik_445[k];

        t_302[k] = -gk_194[k]
                   + f_0 * ik_446[k];

        t_303[k] = -gk_195[k]
                   + f_0 * ik_447[k];

        t_304[k] = -gk_196[k]
                   + f_0 * ik_448[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, gk_197, gk_198, gk_199, gk_200, \
                         gk_201, ik_449, ik_450, ik_451, ik_452, \
                         ik_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -gk_197[k]
                   + f_0 * ik_449[k];

        t_306[k] = -gk_198[k]
                   + f_0 * ik_450[k];

        t_307[k] = -gk_199[k]
                   + f_0 * ik_451[k];

        t_308[k] = -gk_200[k]
                   + f_0 * ik_452[k];

        t_309[k] = -gk_201[k]
                   + f_0 * ik_453[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, gk_202, gk_203, gk_204, gk_205, \
                         gk_206, ik_454, ik_455, ik_456, ik_457, \
                         ik_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -gk_202[k]
                   + f_0 * ik_454[k];

        t_311[k] = -gk_203[k]
                   + f_0 * ik_455[k];

        t_312[k] = -gk_204[k]
                   + f_0 * ik_456[k];

        t_313[k] = -gk_205[k]
                   + f_0 * ik_457[k];

        t_314[k] = -gk_206[k]
                   + f_0 * ik_458[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, gk_207, gk_208, gk_209, gk_210, \
                         gk_211, ik_459, ik_460, ik_461, ik_462, \
                         ik_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -gk_207[k]
                   + f_0 * ik_459[k];

        t_316[k] = -gk_208[k]
                   + f_0 * ik_460[k];

        t_317[k] = -gk_209[k]
                   + f_0 * ik_461[k];

        t_318[k] = -gk_210[k]
                   + f_0 * ik_462[k];

        t_319[k] = -gk_211[k]
                   + f_0 * ik_463[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, t_325, gk_212, gk_213, gk_214, \
                         gk_215, ik_464, ik_465, ik_466, ik_467, ik_468, \
                         ik_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -gk_212[k]
                   + f_0 * ik_464[k];

        t_321[k] = -gk_213[k]
                   + f_0 * ik_465[k];

        t_322[k] = -gk_214[k]
                   + f_0 * ik_466[k];

        t_323[k] = -gk_215[k]
                   + f_0 * ik_467[k];

        t_324[k] = f_0 * ik_468[k];

        t_325[k] = f_0 * ik_469[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, t_332, t_333, ik_470, \
                         ik_471, ik_472, ik_473, ik_474, ik_475, ik_476, \
                         ik_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_0 * ik_470[k];

        t_327[k] = f_0 * ik_471[k];

        t_328[k] = f_0 * ik_472[k];

        t_329[k] = f_0 * ik_473[k];

        t_330[k] = f_0 * ik_474[k];

        t_331[k] = f_0 * ik_475[k];

        t_332[k] = f_0 * ik_476[k];

        t_333[k] = f_0 * ik_477[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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
    const auto *gk_342 = buffer.data(gk + 342);
    const auto *gk_343 = buffer.data(gk + 343);
    const auto *gk_344 = buffer.data(gk + 344);
    const auto *gk_345 = buffer.data(gk + 345);
    const auto *gk_346 = buffer.data(gk + 346);
    const auto *gk_347 = buffer.data(gk + 347);
    const auto *gk_348 = buffer.data(gk + 348);
    const auto *gk_349 = buffer.data(gk + 349);

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

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, t_339, t_340, t_341, ik_478, \
                         ik_479, ik_480, ik_481, ik_482, ik_483, ik_484, \
                         ik_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_0 * ik_478[k];

        t_335[k] = f_0 * ik_479[k];

        t_336[k] = f_0 * ik_480[k];

        t_337[k] = f_0 * ik_481[k];

        t_338[k] = f_0 * ik_482[k];

        t_339[k] = f_0 * ik_483[k];

        t_340[k] = f_0 * ik_484[k];

        t_341[k] = f_0 * ik_485[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, t_347, t_348, t_349, ik_486, \
                         ik_487, ik_488, ik_489, ik_490, ik_491, ik_492, \
                         ik_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_0 * ik_486[k];

        t_343[k] = f_0 * ik_487[k];

        t_344[k] = f_0 * ik_488[k];

        t_345[k] = f_0 * ik_489[k];

        t_346[k] = f_0 * ik_490[k];

        t_347[k] = f_0 * ik_491[k];

        t_348[k] = f_0 * ik_492[k];

        t_349[k] = f_0 * ik_493[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, t_355, t_356, t_357, ik_494, \
                         ik_495, ik_496, ik_497, ik_498, ik_499, ik_500, \
                         ik_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_0 * ik_494[k];

        t_351[k] = f_0 * ik_495[k];

        t_352[k] = f_0 * ik_496[k];

        t_353[k] = f_0 * ik_497[k];

        t_354[k] = f_0 * ik_498[k];

        t_355[k] = f_0 * ik_499[k];

        t_356[k] = f_0 * ik_500[k];

        t_357[k] = f_0 * ik_501[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, t_363, gk_216, gk_217, gk_218, \
                         gk_219, ik_502, ik_503, ik_540, ik_541, ik_542, \
                         ik_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_0 * ik_502[k];

        t_359[k] = f_0 * ik_503[k];

        t_360[k] = -4.0 * gk_216[k]
                   + f_0 * ik_540[k];

        t_361[k] = -4.0 * gk_217[k]
                   + f_0 * ik_541[k];

        t_362[k] = -4.0 * gk_218[k]
                   + f_0 * ik_542[k];

        t_363[k] = -4.0 * gk_219[k]
                   + f_0 * ik_543[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, t_368, gk_220, gk_221, gk_222, gk_223, \
                         gk_224, ik_544, ik_545, ik_546, ik_547, \
                         ik_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = -4.0 * gk_220[k]
                   + f_0 * ik_544[k];

        t_365[k] = -4.0 * gk_221[k]
                   + f_0 * ik_545[k];

        t_366[k] = -4.0 * gk_222[k]
                   + f_0 * ik_546[k];

        t_367[k] = -4.0 * gk_223[k]
                   + f_0 * ik_547[k];

        t_368[k] = -4.0 * gk_224[k]
                   + f_0 * ik_548[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, t_373, gk_225, gk_226, gk_227, gk_228, \
                         gk_229, ik_549, ik_550, ik_551, ik_552, \
                         ik_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = -4.0 * gk_225[k]
                   + f_0 * ik_549[k];

        t_370[k] = -4.0 * gk_226[k]
                   + f_0 * ik_550[k];

        t_371[k] = -4.0 * gk_227[k]
                   + f_0 * ik_551[k];

        t_372[k] = -4.0 * gk_228[k]
                   + f_0 * ik_552[k];

        t_373[k] = -4.0 * gk_229[k]
                   + f_0 * ik_553[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, t_378, gk_230, gk_231, gk_232, gk_233, \
                         gk_234, ik_554, ik_555, ik_556, ik_557, \
                         ik_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = -4.0 * gk_230[k]
                   + f_0 * ik_554[k];

        t_375[k] = -4.0 * gk_231[k]
                   + f_0 * ik_555[k];

        t_376[k] = -4.0 * gk_232[k]
                   + f_0 * ik_556[k];

        t_377[k] = -4.0 * gk_233[k]
                   + f_0 * ik_557[k];

        t_378[k] = -4.0 * gk_234[k]
                   + f_0 * ik_558[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, gk_235, gk_236, gk_237, gk_238, \
                         gk_239, ik_559, ik_560, ik_561, ik_562, \
                         ik_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = -4.0 * gk_235[k]
                   + f_0 * ik_559[k];

        t_380[k] = -4.0 * gk_236[k]
                   + f_0 * ik_560[k];

        t_381[k] = -4.0 * gk_237[k]
                   + f_0 * ik_561[k];

        t_382[k] = -4.0 * gk_238[k]
                   + f_0 * ik_562[k];

        t_383[k] = -4.0 * gk_239[k]
                   + f_0 * ik_563[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, t_388, gk_240, gk_241, gk_242, gk_243, \
                         gk_244, ik_564, ik_565, ik_566, ik_567, \
                         ik_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = -4.0 * gk_240[k]
                   + f_0 * ik_564[k];

        t_385[k] = -4.0 * gk_241[k]
                   + f_0 * ik_565[k];

        t_386[k] = -4.0 * gk_242[k]
                   + f_0 * ik_566[k];

        t_387[k] = -4.0 * gk_243[k]
                   + f_0 * ik_567[k];

        t_388[k] = -4.0 * gk_244[k]
                   + f_0 * ik_568[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, gk_245, gk_246, gk_247, gk_248, \
                         gk_249, ik_569, ik_570, ik_571, ik_572, \
                         ik_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = -4.0 * gk_245[k]
                   + f_0 * ik_569[k];

        t_390[k] = -4.0 * gk_246[k]
                   + f_0 * ik_570[k];

        t_391[k] = -4.0 * gk_247[k]
                   + f_0 * ik_571[k];

        t_392[k] = -4.0 * gk_248[k]
                   + f_0 * ik_572[k];

        t_393[k] = -4.0 * gk_249[k]
                   + f_0 * ik_573[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, gk_250, gk_251, gk_252, gk_253, \
                         gk_254, ik_574, ik_575, ik_576, ik_577, \
                         ik_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = -4.0 * gk_250[k]
                   + f_0 * ik_574[k];

        t_395[k] = -4.0 * gk_251[k]
                   + f_0 * ik_575[k];

        t_396[k] = -3.0 * gk_252[k]
                   + f_0 * ik_576[k];

        t_397[k] = -3.0 * gk_253[k]
                   + f_0 * ik_577[k];

        t_398[k] = -3.0 * gk_254[k]
                   + f_0 * ik_578[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, t_403, gk_255, gk_256, gk_257, gk_258, \
                         gk_259, ik_579, ik_580, ik_581, ik_582, \
                         ik_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = -3.0 * gk_255[k]
                   + f_0 * ik_579[k];

        t_400[k] = -3.0 * gk_256[k]
                   + f_0 * ik_580[k];

        t_401[k] = -3.0 * gk_257[k]
                   + f_0 * ik_581[k];

        t_402[k] = -3.0 * gk_258[k]
                   + f_0 * ik_582[k];

        t_403[k] = -3.0 * gk_259[k]
                   + f_0 * ik_583[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, t_408, gk_260, gk_261, gk_262, gk_263, \
                         gk_264, ik_584, ik_585, ik_586, ik_587, \
                         ik_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = -3.0 * gk_260[k]
                   + f_0 * ik_584[k];

        t_405[k] = -3.0 * gk_261[k]
                   + f_0 * ik_585[k];

        t_406[k] = -3.0 * gk_262[k]
                   + f_0 * ik_586[k];

        t_407[k] = -3.0 * gk_263[k]
                   + f_0 * ik_587[k];

        t_408[k] = -3.0 * gk_264[k]
                   + f_0 * ik_588[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, gk_265, gk_266, gk_267, gk_268, \
                         gk_269, ik_589, ik_590, ik_591, ik_592, \
                         ik_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = -3.0 * gk_265[k]
                   + f_0 * ik_589[k];

        t_410[k] = -3.0 * gk_266[k]
                   + f_0 * ik_590[k];

        t_411[k] = -3.0 * gk_267[k]
                   + f_0 * ik_591[k];

        t_412[k] = -3.0 * gk_268[k]
                   + f_0 * ik_592[k];

        t_413[k] = -3.0 * gk_269[k]
                   + f_0 * ik_593[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, gk_270, gk_271, gk_272, gk_273, \
                         gk_274, ik_594, ik_595, ik_596, ik_597, \
                         ik_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = -3.0 * gk_270[k]
                   + f_0 * ik_594[k];

        t_415[k] = -3.0 * gk_271[k]
                   + f_0 * ik_595[k];

        t_416[k] = -3.0 * gk_272[k]
                   + f_0 * ik_596[k];

        t_417[k] = -3.0 * gk_273[k]
                   + f_0 * ik_597[k];

        t_418[k] = -3.0 * gk_274[k]
                   + f_0 * ik_598[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, t_423, gk_275, gk_276, gk_277, gk_278, \
                         gk_279, ik_599, ik_600, ik_601, ik_602, \
                         ik_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = -3.0 * gk_275[k]
                   + f_0 * ik_599[k];

        t_420[k] = -3.0 * gk_276[k]
                   + f_0 * ik_600[k];

        t_421[k] = -3.0 * gk_277[k]
                   + f_0 * ik_601[k];

        t_422[k] = -3.0 * gk_278[k]
                   + f_0 * ik_602[k];

        t_423[k] = -3.0 * gk_279[k]
                   + f_0 * ik_603[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, t_428, gk_280, gk_281, gk_282, gk_283, \
                         gk_284, ik_604, ik_605, ik_606, ik_607, \
                         ik_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = -3.0 * gk_280[k]
                   + f_0 * ik_604[k];

        t_425[k] = -3.0 * gk_281[k]
                   + f_0 * ik_605[k];

        t_426[k] = -3.0 * gk_282[k]
                   + f_0 * ik_606[k];

        t_427[k] = -3.0 * gk_283[k]
                   + f_0 * ik_607[k];

        t_428[k] = -3.0 * gk_284[k]
                   + f_0 * ik_608[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, gk_285, gk_286, gk_287, gk_288, \
                         gk_289, ik_609, ik_610, ik_611, ik_612, \
                         ik_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = -3.0 * gk_285[k]
                   + f_0 * ik_609[k];

        t_430[k] = -3.0 * gk_286[k]
                   + f_0 * ik_610[k];

        t_431[k] = -3.0 * gk_287[k]
                   + f_0 * ik_611[k];

        t_432[k] = -2.0 * gk_288[k]
                   + f_0 * ik_612[k];

        t_433[k] = -2.0 * gk_289[k]
                   + f_0 * ik_613[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, t_438, gk_290, gk_291, gk_292, gk_293, \
                         gk_294, ik_614, ik_615, ik_616, ik_617, \
                         ik_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = -2.0 * gk_290[k]
                   + f_0 * ik_614[k];

        t_435[k] = -2.0 * gk_291[k]
                   + f_0 * ik_615[k];

        t_436[k] = -2.0 * gk_292[k]
                   + f_0 * ik_616[k];

        t_437[k] = -2.0 * gk_293[k]
                   + f_0 * ik_617[k];

        t_438[k] = -2.0 * gk_294[k]
                   + f_0 * ik_618[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, t_443, gk_295, gk_296, gk_297, gk_298, \
                         gk_299, ik_619, ik_620, ik_621, ik_622, \
                         ik_623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = -2.0 * gk_295[k]
                   + f_0 * ik_619[k];

        t_440[k] = -2.0 * gk_296[k]
                   + f_0 * ik_620[k];

        t_441[k] = -2.0 * gk_297[k]
                   + f_0 * ik_621[k];

        t_442[k] = -2.0 * gk_298[k]
                   + f_0 * ik_622[k];

        t_443[k] = -2.0 * gk_299[k]
                   + f_0 * ik_623[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, gk_300, gk_301, gk_302, gk_303, \
                         gk_304, ik_624, ik_625, ik_626, ik_627, \
                         ik_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = -2.0 * gk_300[k]
                   + f_0 * ik_624[k];

        t_445[k] = -2.0 * gk_301[k]
                   + f_0 * ik_625[k];

        t_446[k] = -2.0 * gk_302[k]
                   + f_0 * ik_626[k];

        t_447[k] = -2.0 * gk_303[k]
                   + f_0 * ik_627[k];

        t_448[k] = -2.0 * gk_304[k]
                   + f_0 * ik_628[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, gk_305, gk_306, gk_307, gk_308, \
                         gk_309, ik_629, ik_630, ik_631, ik_632, \
                         ik_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = -2.0 * gk_305[k]
                   + f_0 * ik_629[k];

        t_450[k] = -2.0 * gk_306[k]
                   + f_0 * ik_630[k];

        t_451[k] = -2.0 * gk_307[k]
                   + f_0 * ik_631[k];

        t_452[k] = -2.0 * gk_308[k]
                   + f_0 * ik_632[k];

        t_453[k] = -2.0 * gk_309[k]
                   + f_0 * ik_633[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, gk_310, gk_311, gk_312, gk_313, \
                         gk_314, ik_634, ik_635, ik_636, ik_637, \
                         ik_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = -2.0 * gk_310[k]
                   + f_0 * ik_634[k];

        t_455[k] = -2.0 * gk_311[k]
                   + f_0 * ik_635[k];

        t_456[k] = -2.0 * gk_312[k]
                   + f_0 * ik_636[k];

        t_457[k] = -2.0 * gk_313[k]
                   + f_0 * ik_637[k];

        t_458[k] = -2.0 * gk_314[k]
                   + f_0 * ik_638[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, t_463, gk_315, gk_316, gk_317, gk_318, \
                         gk_319, ik_639, ik_640, ik_641, ik_642, \
                         ik_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = -2.0 * gk_315[k]
                   + f_0 * ik_639[k];

        t_460[k] = -2.0 * gk_316[k]
                   + f_0 * ik_640[k];

        t_461[k] = -2.0 * gk_317[k]
                   + f_0 * ik_641[k];

        t_462[k] = -2.0 * gk_318[k]
                   + f_0 * ik_642[k];

        t_463[k] = -2.0 * gk_319[k]
                   + f_0 * ik_643[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, t_468, gk_320, gk_321, gk_322, gk_323, \
                         gk_324, ik_644, ik_645, ik_646, ik_647, \
                         ik_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = -2.0 * gk_320[k]
                   + f_0 * ik_644[k];

        t_465[k] = -2.0 * gk_321[k]
                   + f_0 * ik_645[k];

        t_466[k] = -2.0 * gk_322[k]
                   + f_0 * ik_646[k];

        t_467[k] = -2.0 * gk_323[k]
                   + f_0 * ik_647[k];

        t_468[k] = -gk_324[k]
                   + f_0 * ik_648[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, t_473, gk_325, gk_326, gk_327, gk_328, \
                         gk_329, ik_649, ik_650, ik_651, ik_652, \
                         ik_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = -gk_325[k]
                   + f_0 * ik_649[k];

        t_470[k] = -gk_326[k]
                   + f_0 * ik_650[k];

        t_471[k] = -gk_327[k]
                   + f_0 * ik_651[k];

        t_472[k] = -gk_328[k]
                   + f_0 * ik_652[k];

        t_473[k] = -gk_329[k]
                   + f_0 * ik_653[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, t_478, gk_330, gk_331, gk_332, gk_333, \
                         gk_334, ik_654, ik_655, ik_656, ik_657, \
                         ik_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = -gk_330[k]
                   + f_0 * ik_654[k];

        t_475[k] = -gk_331[k]
                   + f_0 * ik_655[k];

        t_476[k] = -gk_332[k]
                   + f_0 * ik_656[k];

        t_477[k] = -gk_333[k]
                   + f_0 * ik_657[k];

        t_478[k] = -gk_334[k]
                   + f_0 * ik_658[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, t_483, gk_335, gk_336, gk_337, gk_338, \
                         gk_339, ik_659, ik_660, ik_661, ik_662, \
                         ik_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = -gk_335[k]
                   + f_0 * ik_659[k];

        t_480[k] = -gk_336[k]
                   + f_0 * ik_660[k];

        t_481[k] = -gk_337[k]
                   + f_0 * ik_661[k];

        t_482[k] = -gk_338[k]
                   + f_0 * ik_662[k];

        t_483[k] = -gk_339[k]
                   + f_0 * ik_663[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, t_488, gk_340, gk_341, gk_342, gk_343, \
                         gk_344, ik_664, ik_665, ik_666, ik_667, \
                         ik_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = -gk_340[k]
                   + f_0 * ik_664[k];

        t_485[k] = -gk_341[k]
                   + f_0 * ik_665[k];

        t_486[k] = -gk_342[k]
                   + f_0 * ik_666[k];

        t_487[k] = -gk_343[k]
                   + f_0 * ik_667[k];

        t_488[k] = -gk_344[k]
                   + f_0 * ik_668[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, gk_345, gk_346, gk_347, gk_348, \
                         gk_349, ik_669, ik_670, ik_671, ik_672, \
                         ik_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = -gk_345[k]
                   + f_0 * ik_669[k];

        t_490[k] = -gk_346[k]
                   + f_0 * ik_670[k];

        t_491[k] = -gk_347[k]
                   + f_0 * ik_671[k];

        t_492[k] = -gk_348[k]
                   + f_0 * ik_672[k];

        t_493[k] = -gk_349[k]
                   + f_0 * ik_673[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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

#pragma omp simd aligned(t_494, t_495, t_496, t_497, t_498, gk_350, gk_351, gk_352, gk_353, \
                         gk_354, ik_674, ik_675, ik_676, ik_677, \
                         ik_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = -gk_350[k]
                   + f_0 * ik_674[k];

        t_495[k] = -gk_351[k]
                   + f_0 * ik_675[k];

        t_496[k] = -gk_352[k]
                   + f_0 * ik_676[k];

        t_497[k] = -gk_353[k]
                   + f_0 * ik_677[k];

        t_498[k] = -gk_354[k]
                   + f_0 * ik_678[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, t_503, gk_355, gk_356, gk_357, gk_358, \
                         gk_359, ik_679, ik_680, ik_681, ik_682, \
                         ik_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = -gk_355[k]
                   + f_0 * ik_679[k];

        t_500[k] = -gk_356[k]
                   + f_0 * ik_680[k];

        t_501[k] = -gk_357[k]
                   + f_0 * ik_681[k];

        t_502[k] = -gk_358[k]
                   + f_0 * ik_682[k];

        t_503[k] = -gk_359[k]
                   + f_0 * ik_683[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, t_509, t_510, t_511, ik_684, \
                         ik_685, ik_686, ik_687, ik_688, ik_689, ik_690, \
                         ik_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_0 * ik_684[k];

        t_505[k] = f_0 * ik_685[k];

        t_506[k] = f_0 * ik_686[k];

        t_507[k] = f_0 * ik_687[k];

        t_508[k] = f_0 * ik_688[k];

        t_509[k] = f_0 * ik_689[k];

        t_510[k] = f_0 * ik_690[k];

        t_511[k] = f_0 * ik_691[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, t_517, t_518, t_519, ik_692, \
                         ik_693, ik_694, ik_695, ik_696, ik_697, ik_698, \
                         ik_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_0 * ik_692[k];

        t_513[k] = f_0 * ik_693[k];

        t_514[k] = f_0 * ik_694[k];

        t_515[k] = f_0 * ik_695[k];

        t_516[k] = f_0 * ik_696[k];

        t_517[k] = f_0 * ik_697[k];

        t_518[k] = f_0 * ik_698[k];

        t_519[k] = f_0 * ik_699[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, t_525, t_526, t_527, ik_700, \
                         ik_701, ik_702, ik_703, ik_704, ik_705, ik_706, \
                         ik_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_0 * ik_700[k];

        t_521[k] = f_0 * ik_701[k];

        t_522[k] = f_0 * ik_702[k];

        t_523[k] = f_0 * ik_703[k];

        t_524[k] = f_0 * ik_704[k];

        t_525[k] = f_0 * ik_705[k];

        t_526[k] = f_0 * ik_706[k];

        t_527[k] = f_0 * ik_707[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, t_533, t_534, t_535, ik_708, \
                         ik_709, ik_710, ik_711, ik_712, ik_713, ik_714, \
                         ik_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_0 * ik_708[k];

        t_529[k] = f_0 * ik_709[k];

        t_530[k] = f_0 * ik_710[k];

        t_531[k] = f_0 * ik_711[k];

        t_532[k] = f_0 * ik_712[k];

        t_533[k] = f_0 * ik_713[k];

        t_534[k] = f_0 * ik_714[k];

        t_535[k] = f_0 * ik_715[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, t_541, gk_360, gk_361, ik_716, \
                         ik_717, ik_718, ik_719, ik_756, ik_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_0 * ik_716[k];

        t_537[k] = f_0 * ik_717[k];

        t_538[k] = f_0 * ik_718[k];

        t_539[k] = f_0 * ik_719[k];

        t_540[k] = -5.0 * gk_360[k]
                   + f_0 * ik_756[k];

        t_541[k] = -5.0 * gk_361[k]
                   + f_0 * ik_757[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, gk_362, gk_363, gk_364, gk_365, \
                         gk_366, ik_758, ik_759, ik_760, ik_761, \
                         ik_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = -5.0 * gk_362[k]
                   + f_0 * ik_758[k];

        t_543[k] = -5.0 * gk_363[k]
                   + f_0 * ik_759[k];

        t_544[k] = -5.0 * gk_364[k]
                   + f_0 * ik_760[k];

        t_545[k] = -5.0 * gk_365[k]
                   + f_0 * ik_761[k];

        t_546[k] = -5.0 * gk_366[k]
                   + f_0 * ik_762[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, gk_367, gk_368, gk_369, gk_370, \
                         gk_371, ik_763, ik_764, ik_765, ik_766, \
                         ik_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = -5.0 * gk_367[k]
                   + f_0 * ik_763[k];

        t_548[k] = -5.0 * gk_368[k]
                   + f_0 * ik_764[k];

        t_549[k] = -5.0 * gk_369[k]
                   + f_0 * ik_765[k];

        t_550[k] = -5.0 * gk_370[k]
                   + f_0 * ik_766[k];

        t_551[k] = -5.0 * gk_371[k]
                   + f_0 * ik_767[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, t_556, gk_372, gk_373, gk_374, gk_375, \
                         gk_376, ik_768, ik_769, ik_770, ik_771, \
                         ik_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = -5.0 * gk_372[k]
                   + f_0 * ik_768[k];

        t_553[k] = -5.0 * gk_373[k]
                   + f_0 * ik_769[k];

        t_554[k] = -5.0 * gk_374[k]
                   + f_0 * ik_770[k];

        t_555[k] = -5.0 * gk_375[k]
                   + f_0 * ik_771[k];

        t_556[k] = -5.0 * gk_376[k]
                   + f_0 * ik_772[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, gk_377, gk_378, gk_379, gk_380, \
                         gk_381, ik_773, ik_774, ik_775, ik_776, \
                         ik_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = -5.0 * gk_377[k]
                   + f_0 * ik_773[k];

        t_558[k] = -5.0 * gk_378[k]
                   + f_0 * ik_774[k];

        t_559[k] = -5.0 * gk_379[k]
                   + f_0 * ik_775[k];

        t_560[k] = -5.0 * gk_380[k]
                   + f_0 * ik_776[k];

        t_561[k] = -5.0 * gk_381[k]
                   + f_0 * ik_777[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, gk_382, gk_383, gk_384, gk_385, \
                         gk_386, ik_778, ik_779, ik_780, ik_781, \
                         ik_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = -5.0 * gk_382[k]
                   + f_0 * ik_778[k];

        t_563[k] = -5.0 * gk_383[k]
                   + f_0 * ik_779[k];

        t_564[k] = -5.0 * gk_384[k]
                   + f_0 * ik_780[k];

        t_565[k] = -5.0 * gk_385[k]
                   + f_0 * ik_781[k];

        t_566[k] = -5.0 * gk_386[k]
                   + f_0 * ik_782[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, gk_387, gk_388, gk_389, gk_390, \
                         gk_391, ik_783, ik_784, ik_785, ik_786, \
                         ik_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -5.0 * gk_387[k]
                   + f_0 * ik_783[k];

        t_568[k] = -5.0 * gk_388[k]
                   + f_0 * ik_784[k];

        t_569[k] = -5.0 * gk_389[k]
                   + f_0 * ik_785[k];

        t_570[k] = -5.0 * gk_390[k]
                   + f_0 * ik_786[k];

        t_571[k] = -5.0 * gk_391[k]
                   + f_0 * ik_787[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, gk_392, gk_393, gk_394, gk_395, \
                         gk_396, ik_788, ik_789, ik_790, ik_791, \
                         ik_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = -5.0 * gk_392[k]
                   + f_0 * ik_788[k];

        t_573[k] = -5.0 * gk_393[k]
                   + f_0 * ik_789[k];

        t_574[k] = -5.0 * gk_394[k]
                   + f_0 * ik_790[k];

        t_575[k] = -5.0 * gk_395[k]
                   + f_0 * ik_791[k];

        t_576[k] = -4.0 * gk_396[k]
                   + f_0 * ik_792[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, gk_397, gk_398, gk_399, gk_400, \
                         gk_401, ik_793, ik_794, ik_795, ik_796, \
                         ik_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = -4.0 * gk_397[k]
                   + f_0 * ik_793[k];

        t_578[k] = -4.0 * gk_398[k]
                   + f_0 * ik_794[k];

        t_579[k] = -4.0 * gk_399[k]
                   + f_0 * ik_795[k];

        t_580[k] = -4.0 * gk_400[k]
                   + f_0 * ik_796[k];

        t_581[k] = -4.0 * gk_401[k]
                   + f_0 * ik_797[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, gk_402, gk_403, gk_404, gk_405, \
                         gk_406, ik_798, ik_799, ik_800, ik_801, \
                         ik_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = -4.0 * gk_402[k]
                   + f_0 * ik_798[k];

        t_583[k] = -4.0 * gk_403[k]
                   + f_0 * ik_799[k];

        t_584[k] = -4.0 * gk_404[k]
                   + f_0 * ik_800[k];

        t_585[k] = -4.0 * gk_405[k]
                   + f_0 * ik_801[k];

        t_586[k] = -4.0 * gk_406[k]
                   + f_0 * ik_802[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, t_591, gk_407, gk_408, gk_409, gk_410, \
                         gk_411, ik_803, ik_804, ik_805, ik_806, \
                         ik_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = -4.0 * gk_407[k]
                   + f_0 * ik_803[k];

        t_588[k] = -4.0 * gk_408[k]
                   + f_0 * ik_804[k];

        t_589[k] = -4.0 * gk_409[k]
                   + f_0 * ik_805[k];

        t_590[k] = -4.0 * gk_410[k]
                   + f_0 * ik_806[k];

        t_591[k] = -4.0 * gk_411[k]
                   + f_0 * ik_807[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, t_596, gk_412, gk_413, gk_414, gk_415, \
                         gk_416, ik_808, ik_809, ik_810, ik_811, \
                         ik_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = -4.0 * gk_412[k]
                   + f_0 * ik_808[k];

        t_593[k] = -4.0 * gk_413[k]
                   + f_0 * ik_809[k];

        t_594[k] = -4.0 * gk_414[k]
                   + f_0 * ik_810[k];

        t_595[k] = -4.0 * gk_415[k]
                   + f_0 * ik_811[k];

        t_596[k] = -4.0 * gk_416[k]
                   + f_0 * ik_812[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, t_600, t_601, gk_417, gk_418, gk_419, gk_420, \
                         gk_421, ik_813, ik_814, ik_815, ik_816, \
                         ik_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = -4.0 * gk_417[k]
                   + f_0 * ik_813[k];

        t_598[k] = -4.0 * gk_418[k]
                   + f_0 * ik_814[k];

        t_599[k] = -4.0 * gk_419[k]
                   + f_0 * ik_815[k];

        t_600[k] = -4.0 * gk_420[k]
                   + f_0 * ik_816[k];

        t_601[k] = -4.0 * gk_421[k]
                   + f_0 * ik_817[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, gk_422, gk_423, gk_424, gk_425, \
                         gk_426, ik_818, ik_819, ik_820, ik_821, \
                         ik_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = -4.0 * gk_422[k]
                   + f_0 * ik_818[k];

        t_603[k] = -4.0 * gk_423[k]
                   + f_0 * ik_819[k];

        t_604[k] = -4.0 * gk_424[k]
                   + f_0 * ik_820[k];

        t_605[k] = -4.0 * gk_425[k]
                   + f_0 * ik_821[k];

        t_606[k] = -4.0 * gk_426[k]
                   + f_0 * ik_822[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, t_611, gk_427, gk_428, gk_429, gk_430, \
                         gk_431, ik_823, ik_824, ik_825, ik_826, \
                         ik_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = -4.0 * gk_427[k]
                   + f_0 * ik_823[k];

        t_608[k] = -4.0 * gk_428[k]
                   + f_0 * ik_824[k];

        t_609[k] = -4.0 * gk_429[k]
                   + f_0 * ik_825[k];

        t_610[k] = -4.0 * gk_430[k]
                   + f_0 * ik_826[k];

        t_611[k] = -4.0 * gk_431[k]
                   + f_0 * ik_827[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, gk_432, gk_433, gk_434, gk_435, \
                         gk_436, ik_828, ik_829, ik_830, ik_831, \
                         ik_832 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = -3.0 * gk_432[k]
                   + f_0 * ik_828[k];

        t_613[k] = -3.0 * gk_433[k]
                   + f_0 * ik_829[k];

        t_614[k] = -3.0 * gk_434[k]
                   + f_0 * ik_830[k];

        t_615[k] = -3.0 * gk_435[k]
                   + f_0 * ik_831[k];

        t_616[k] = -3.0 * gk_436[k]
                   + f_0 * ik_832[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, gk_437, gk_438, gk_439, gk_440, \
                         gk_441, ik_833, ik_834, ik_835, ik_836, \
                         ik_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = -3.0 * gk_437[k]
                   + f_0 * ik_833[k];

        t_618[k] = -3.0 * gk_438[k]
                   + f_0 * ik_834[k];

        t_619[k] = -3.0 * gk_439[k]
                   + f_0 * ik_835[k];

        t_620[k] = -3.0 * gk_440[k]
                   + f_0 * ik_836[k];

        t_621[k] = -3.0 * gk_441[k]
                   + f_0 * ik_837[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, gk_442, gk_443, gk_444, gk_445, \
                         gk_446, ik_838, ik_839, ik_840, ik_841, \
                         ik_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = -3.0 * gk_442[k]
                   + f_0 * ik_838[k];

        t_623[k] = -3.0 * gk_443[k]
                   + f_0 * ik_839[k];

        t_624[k] = -3.0 * gk_444[k]
                   + f_0 * ik_840[k];

        t_625[k] = -3.0 * gk_445[k]
                   + f_0 * ik_841[k];

        t_626[k] = -3.0 * gk_446[k]
                   + f_0 * ik_842[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, t_631, gk_447, gk_448, gk_449, gk_450, \
                         gk_451, ik_843, ik_844, ik_845, ik_846, \
                         ik_847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = -3.0 * gk_447[k]
                   + f_0 * ik_843[k];

        t_628[k] = -3.0 * gk_448[k]
                   + f_0 * ik_844[k];

        t_629[k] = -3.0 * gk_449[k]
                   + f_0 * ik_845[k];

        t_630[k] = -3.0 * gk_450[k]
                   + f_0 * ik_846[k];

        t_631[k] = -3.0 * gk_451[k]
                   + f_0 * ik_847[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, t_636, gk_452, gk_453, gk_454, gk_455, \
                         gk_456, ik_848, ik_849, ik_850, ik_851, \
                         ik_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = -3.0 * gk_452[k]
                   + f_0 * ik_848[k];

        t_633[k] = -3.0 * gk_453[k]
                   + f_0 * ik_849[k];

        t_634[k] = -3.0 * gk_454[k]
                   + f_0 * ik_850[k];

        t_635[k] = -3.0 * gk_455[k]
                   + f_0 * ik_851[k];

        t_636[k] = -3.0 * gk_456[k]
                   + f_0 * ik_852[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, t_641, gk_457, gk_458, gk_459, gk_460, \
                         gk_461, ik_853, ik_854, ik_855, ik_856, \
                         ik_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = -3.0 * gk_457[k]
                   + f_0 * ik_853[k];

        t_638[k] = -3.0 * gk_458[k]
                   + f_0 * ik_854[k];

        t_639[k] = -3.0 * gk_459[k]
                   + f_0 * ik_855[k];

        t_640[k] = -3.0 * gk_460[k]
                   + f_0 * ik_856[k];

        t_641[k] = -3.0 * gk_461[k]
                   + f_0 * ik_857[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, t_646, gk_462, gk_463, gk_464, gk_465, \
                         gk_466, ik_858, ik_859, ik_860, ik_861, \
                         ik_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = -3.0 * gk_462[k]
                   + f_0 * ik_858[k];

        t_643[k] = -3.0 * gk_463[k]
                   + f_0 * ik_859[k];

        t_644[k] = -3.0 * gk_464[k]
                   + f_0 * ik_860[k];

        t_645[k] = -3.0 * gk_465[k]
                   + f_0 * ik_861[k];

        t_646[k] = -3.0 * gk_466[k]
                   + f_0 * ik_862[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, t_651, gk_467, gk_468, gk_469, gk_470, \
                         gk_471, ik_863, ik_864, ik_865, ik_866, \
                         ik_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = -3.0 * gk_467[k]
                   + f_0 * ik_863[k];

        t_648[k] = -2.0 * gk_468[k]
                   + f_0 * ik_864[k];

        t_649[k] = -2.0 * gk_469[k]
                   + f_0 * ik_865[k];

        t_650[k] = -2.0 * gk_470[k]
                   + f_0 * ik_866[k];

        t_651[k] = -2.0 * gk_471[k]
                   + f_0 * ik_867[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, t_656, gk_472, gk_473, gk_474, gk_475, \
                         gk_476, ik_868, ik_869, ik_870, ik_871, \
                         ik_872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = -2.0 * gk_472[k]
                   + f_0 * ik_868[k];

        t_653[k] = -2.0 * gk_473[k]
                   + f_0 * ik_869[k];

        t_654[k] = -2.0 * gk_474[k]
                   + f_0 * ik_870[k];

        t_655[k] = -2.0 * gk_475[k]
                   + f_0 * ik_871[k];

        t_656[k] = -2.0 * gk_476[k]
                   + f_0 * ik_872[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, gk_477, gk_478, gk_479, gk_480, \
                         gk_481, ik_873, ik_874, ik_875, ik_876, \
                         ik_877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = -2.0 * gk_477[k]
                   + f_0 * ik_873[k];

        t_658[k] = -2.0 * gk_478[k]
                   + f_0 * ik_874[k];

        t_659[k] = -2.0 * gk_479[k]
                   + f_0 * ik_875[k];

        t_660[k] = -2.0 * gk_480[k]
                   + f_0 * ik_876[k];

        t_661[k] = -2.0 * gk_481[k]
                   + f_0 * ik_877[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, t_666, gk_482, gk_483, gk_484, gk_485, \
                         gk_486, ik_878, ik_879, ik_880, ik_881, \
                         ik_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = -2.0 * gk_482[k]
                   + f_0 * ik_878[k];

        t_663[k] = -2.0 * gk_483[k]
                   + f_0 * ik_879[k];

        t_664[k] = -2.0 * gk_484[k]
                   + f_0 * ik_880[k];

        t_665[k] = -2.0 * gk_485[k]
                   + f_0 * ik_881[k];

        t_666[k] = -2.0 * gk_486[k]
                   + f_0 * ik_882[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, t_670, t_671, gk_487, gk_488, gk_489, gk_490, \
                         gk_491, ik_883, ik_884, ik_885, ik_886, \
                         ik_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = -2.0 * gk_487[k]
                   + f_0 * ik_883[k];

        t_668[k] = -2.0 * gk_488[k]
                   + f_0 * ik_884[k];

        t_669[k] = -2.0 * gk_489[k]
                   + f_0 * ik_885[k];

        t_670[k] = -2.0 * gk_490[k]
                   + f_0 * ik_886[k];

        t_671[k] = -2.0 * gk_491[k]
                   + f_0 * ik_887[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, t_676, gk_492, gk_493, gk_494, gk_495, \
                         gk_496, ik_888, ik_889, ik_890, ik_891, \
                         ik_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = -2.0 * gk_492[k]
                   + f_0 * ik_888[k];

        t_673[k] = -2.0 * gk_493[k]
                   + f_0 * ik_889[k];

        t_674[k] = -2.0 * gk_494[k]
                   + f_0 * ik_890[k];

        t_675[k] = -2.0 * gk_495[k]
                   + f_0 * ik_891[k];

        t_676[k] = -2.0 * gk_496[k]
                   + f_0 * ik_892[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, t_680, t_681, gk_497, gk_498, gk_499, gk_500, \
                         gk_501, ik_893, ik_894, ik_895, ik_896, \
                         ik_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = -2.0 * gk_497[k]
                   + f_0 * ik_893[k];

        t_678[k] = -2.0 * gk_498[k]
                   + f_0 * ik_894[k];

        t_679[k] = -2.0 * gk_499[k]
                   + f_0 * ik_895[k];

        t_680[k] = -2.0 * gk_500[k]
                   + f_0 * ik_896[k];

        t_681[k] = -2.0 * gk_501[k]
                   + f_0 * ik_897[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, t_686, gk_502, gk_503, gk_504, gk_505, \
                         gk_506, ik_898, ik_899, ik_900, ik_901, \
                         ik_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = -2.0 * gk_502[k]
                   + f_0 * ik_898[k];

        t_683[k] = -2.0 * gk_503[k]
                   + f_0 * ik_899[k];

        t_684[k] = -gk_504[k]
                   + f_0 * ik_900[k];

        t_685[k] = -gk_505[k]
                   + f_0 * ik_901[k];

        t_686[k] = -gk_506[k]
                   + f_0 * ik_902[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, t_691, gk_507, gk_508, gk_509, gk_510, \
                         gk_511, ik_903, ik_904, ik_905, ik_906, \
                         ik_907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = -gk_507[k]
                   + f_0 * ik_903[k];

        t_688[k] = -gk_508[k]
                   + f_0 * ik_904[k];

        t_689[k] = -gk_509[k]
                   + f_0 * ik_905[k];

        t_690[k] = -gk_510[k]
                   + f_0 * ik_906[k];

        t_691[k] = -gk_511[k]
                   + f_0 * ik_907[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, t_696, gk_512, gk_513, gk_514, gk_515, \
                         gk_516, ik_908, ik_909, ik_910, ik_911, \
                         ik_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = -gk_512[k]
                   + f_0 * ik_908[k];

        t_693[k] = -gk_513[k]
                   + f_0 * ik_909[k];

        t_694[k] = -gk_514[k]
                   + f_0 * ik_910[k];

        t_695[k] = -gk_515[k]
                   + f_0 * ik_911[k];

        t_696[k] = -gk_516[k]
                   + f_0 * ik_912[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, t_701, gk_517, gk_518, gk_519, gk_520, \
                         gk_521, ik_913, ik_914, ik_915, ik_916, \
                         ik_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = -gk_517[k]
                   + f_0 * ik_913[k];

        t_698[k] = -gk_518[k]
                   + f_0 * ik_914[k];

        t_699[k] = -gk_519[k]
                   + f_0 * ik_915[k];

        t_700[k] = -gk_520[k]
                   + f_0 * ik_916[k];

        t_701[k] = -gk_521[k]
                   + f_0 * ik_917[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, t_706, gk_522, gk_523, gk_524, gk_525, \
                         gk_526, ik_918, ik_919, ik_920, ik_921, \
                         ik_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = -gk_522[k]
                   + f_0 * ik_918[k];

        t_703[k] = -gk_523[k]
                   + f_0 * ik_919[k];

        t_704[k] = -gk_524[k]
                   + f_0 * ik_920[k];

        t_705[k] = -gk_525[k]
                   + f_0 * ik_921[k];

        t_706[k] = -gk_526[k]
                   + f_0 * ik_922[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, t_710, t_711, gk_527, gk_528, gk_529, gk_530, \
                         gk_531, ik_923, ik_924, ik_925, ik_926, \
                         ik_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = -gk_527[k]
                   + f_0 * ik_923[k];

        t_708[k] = -gk_528[k]
                   + f_0 * ik_924[k];

        t_709[k] = -gk_529[k]
                   + f_0 * ik_925[k];

        t_710[k] = -gk_530[k]
                   + f_0 * ik_926[k];

        t_711[k] = -gk_531[k]
                   + f_0 * ik_927[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, t_716, gk_532, gk_533, gk_534, gk_535, \
                         gk_536, ik_928, ik_929, ik_930, ik_931, \
                         ik_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = -gk_532[k]
                   + f_0 * ik_928[k];

        t_713[k] = -gk_533[k]
                   + f_0 * ik_929[k];

        t_714[k] = -gk_534[k]
                   + f_0 * ik_930[k];

        t_715[k] = -gk_535[k]
                   + f_0 * ik_931[k];

        t_716[k] = -gk_536[k]
                   + f_0 * ik_932[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, t_720, t_721, t_722, gk_537, gk_538, gk_539, \
                         ik_933, ik_934, ik_935, ik_936, ik_937, \
                         ik_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = -gk_537[k]
                   + f_0 * ik_933[k];

        t_718[k] = -gk_538[k]
                   + f_0 * ik_934[k];

        t_719[k] = -gk_539[k]
                   + f_0 * ik_935[k];

        t_720[k] = f_0 * ik_936[k];

        t_721[k] = f_0 * ik_937[k];

        t_722[k] = f_0 * ik_938[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, t_726, t_727, t_728, t_729, t_730, ik_939, \
                         ik_940, ik_941, ik_942, ik_943, ik_944, ik_945, \
                         ik_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_0 * ik_939[k];

        t_724[k] = f_0 * ik_940[k];

        t_725[k] = f_0 * ik_941[k];

        t_726[k] = f_0 * ik_942[k];

        t_727[k] = f_0 * ik_943[k];

        t_728[k] = f_0 * ik_944[k];

        t_729[k] = f_0 * ik_945[k];

        t_730[k] = f_0 * ik_946[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, t_734, t_735, t_736, t_737, t_738, ik_947, \
                         ik_948, ik_949, ik_950, ik_951, ik_952, ik_953, \
                         ik_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_0 * ik_947[k];

        t_732[k] = f_0 * ik_948[k];

        t_733[k] = f_0 * ik_949[k];

        t_734[k] = f_0 * ik_950[k];

        t_735[k] = f_0 * ik_951[k];

        t_736[k] = f_0 * ik_952[k];

        t_737[k] = f_0 * ik_953[k];

        t_738[k] = f_0 * ik_954[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, t_743, t_744, t_745, t_746, ik_955, \
                         ik_956, ik_957, ik_958, ik_959, ik_960, ik_961, \
                         ik_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_0 * ik_955[k];

        t_740[k] = f_0 * ik_956[k];

        t_741[k] = f_0 * ik_957[k];

        t_742[k] = f_0 * ik_958[k];

        t_743[k] = f_0 * ik_959[k];

        t_744[k] = f_0 * ik_960[k];

        t_745[k] = f_0 * ik_961[k];

        t_746[k] = f_0 * ik_962[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, t_750, t_751, t_752, t_753, t_754, ik_963, \
                         ik_964, ik_965, ik_966, ik_967, ik_968, ik_969, \
                         ik_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_0 * ik_963[k];

        t_748[k] = f_0 * ik_964[k];

        t_749[k] = f_0 * ik_965[k];

        t_750[k] = f_0 * ik_966[k];

        t_751[k] = f_0 * ik_967[k];

        t_752[k] = f_0 * ik_968[k];

        t_753[k] = f_0 * ik_969[k];

        t_754[k] = f_0 * ik_970[k];
    }

#pragma omp simd aligned(t_755, ik_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_0 * ik_971[k];
    }
}

auto
compute_prim_geom_10_hk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t gk, const size_t ik,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hk_electron_repulsion_1_piece0(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_1_piece1(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_1_piece2(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_1_piece3(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_1_piece4(buffer, target, gk, ik, ncols, alpha);
}

static auto
compute_prim_geom_10_hk_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ik_72, ik_73, ik_74, ik_75, \
                         ik_76, ik_77, ik_78, ik_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ik_72[k];

        t_1[k] = f_0 * ik_73[k];

        t_2[k] = f_0 * ik_74[k];

        t_3[k] = f_0 * ik_75[k];

        t_4[k] = f_0 * ik_76[k];

        t_5[k] = f_0 * ik_77[k];

        t_6[k] = f_0 * ik_78[k];

        t_7[k] = f_0 * ik_79[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, ik_80, ik_81, ik_82, \
                         ik_83, ik_84, ik_85, ik_86, ik_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ik_80[k];

        t_9[k] = f_0 * ik_81[k];

        t_10[k] = f_0 * ik_82[k];

        t_11[k] = f_0 * ik_83[k];

        t_12[k] = f_0 * ik_84[k];

        t_13[k] = f_0 * ik_85[k];

        t_14[k] = f_0 * ik_86[k];

        t_15[k] = f_0 * ik_87[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, ik_88, ik_89, ik_90, \
                         ik_91, ik_92, ik_93, ik_94, ik_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * ik_88[k];

        t_17[k] = f_0 * ik_89[k];

        t_18[k] = f_0 * ik_90[k];

        t_19[k] = f_0 * ik_91[k];

        t_20[k] = f_0 * ik_92[k];

        t_21[k] = f_0 * ik_93[k];

        t_22[k] = f_0 * ik_94[k];

        t_23[k] = f_0 * ik_95[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, ik_96, ik_97, ik_98, \
                         ik_99, ik_100, ik_101, ik_102, ik_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * ik_96[k];

        t_25[k] = f_0 * ik_97[k];

        t_26[k] = f_0 * ik_98[k];

        t_27[k] = f_0 * ik_99[k];

        t_28[k] = f_0 * ik_100[k];

        t_29[k] = f_0 * ik_101[k];

        t_30[k] = f_0 * ik_102[k];

        t_31[k] = f_0 * ik_103[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, ik_104, ik_105, \
                         ik_106, ik_107, ik_144, ik_145, ik_146, \
                         ik_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ik_104[k];

        t_33[k] = f_0 * ik_105[k];

        t_34[k] = f_0 * ik_106[k];

        t_35[k] = f_0 * ik_107[k];

        t_36[k] = f_0 * ik_144[k];

        t_37[k] = f_0 * ik_145[k];

        t_38[k] = f_0 * ik_146[k];

        t_39[k] = f_0 * ik_147[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, ik_148, ik_149, \
                         ik_150, ik_151, ik_152, ik_153, ik_154, \
                         ik_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * ik_148[k];

        t_41[k] = f_0 * ik_149[k];

        t_42[k] = f_0 * ik_150[k];

        t_43[k] = f_0 * ik_151[k];

        t_44[k] = f_0 * ik_152[k];

        t_45[k] = f_0 * ik_153[k];

        t_46[k] = f_0 * ik_154[k];

        t_47[k] = f_0 * ik_155[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, ik_156, ik_157, \
                         ik_158, ik_159, ik_160, ik_161, ik_162, \
                         ik_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * ik_156[k];

        t_49[k] = f_0 * ik_157[k];

        t_50[k] = f_0 * ik_158[k];

        t_51[k] = f_0 * ik_159[k];

        t_52[k] = f_0 * ik_160[k];

        t_53[k] = f_0 * ik_161[k];

        t_54[k] = f_0 * ik_162[k];

        t_55[k] = f_0 * ik_163[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, ik_164, ik_165, \
                         ik_166, ik_167, ik_168, ik_169, ik_170, \
                         ik_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * ik_164[k];

        t_57[k] = f_0 * ik_165[k];

        t_58[k] = f_0 * ik_166[k];

        t_59[k] = f_0 * ik_167[k];

        t_60[k] = f_0 * ik_168[k];

        t_61[k] = f_0 * ik_169[k];

        t_62[k] = f_0 * ik_170[k];

        t_63[k] = f_0 * ik_171[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, ik_172, ik_173, \
                         ik_174, ik_175, ik_176, ik_177, ik_178, \
                         ik_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * ik_172[k];

        t_65[k] = f_0 * ik_173[k];

        t_66[k] = f_0 * ik_174[k];

        t_67[k] = f_0 * ik_175[k];

        t_68[k] = f_0 * ik_176[k];

        t_69[k] = f_0 * ik_177[k];

        t_70[k] = f_0 * ik_178[k];

        t_71[k] = f_0 * ik_179[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, gk_0, gk_1, gk_2, gk_3, gk_4, ik_180, \
                         ik_181, ik_182, ik_183, ik_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -gk_0[k]
                  + f_0 * ik_180[k];

        t_73[k] = -gk_1[k]
                  + f_0 * ik_181[k];

        t_74[k] = -gk_2[k]
                  + f_0 * ik_182[k];

        t_75[k] = -gk_3[k]
                  + f_0 * ik_183[k];

        t_76[k] = -gk_4[k]
                  + f_0 * ik_184[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, gk_5, gk_6, gk_7, gk_8, gk_9, ik_185, \
                         ik_186, ik_187, ik_188, ik_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -gk_5[k]
                  + f_0 * ik_185[k];

        t_78[k] = -gk_6[k]
                  + f_0 * ik_186[k];

        t_79[k] = -gk_7[k]
                  + f_0 * ik_187[k];

        t_80[k] = -gk_8[k]
                  + f_0 * ik_188[k];

        t_81[k] = -gk_9[k]
                  + f_0 * ik_189[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, gk_10, gk_11, gk_12, gk_13, gk_14, \
                         ik_190, ik_191, ik_192, ik_193, ik_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -gk_10[k]
                  + f_0 * ik_190[k];

        t_83[k] = -gk_11[k]
                  + f_0 * ik_191[k];

        t_84[k] = -gk_12[k]
                  + f_0 * ik_192[k];

        t_85[k] = -gk_13[k]
                  + f_0 * ik_193[k];

        t_86[k] = -gk_14[k]
                  + f_0 * ik_194[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, gk_15, gk_16, gk_17, gk_18, gk_19, \
                         ik_195, ik_196, ik_197, ik_198, ik_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -gk_15[k]
                  + f_0 * ik_195[k];

        t_88[k] = -gk_16[k]
                  + f_0 * ik_196[k];

        t_89[k] = -gk_17[k]
                  + f_0 * ik_197[k];

        t_90[k] = -gk_18[k]
                  + f_0 * ik_198[k];

        t_91[k] = -gk_19[k]
                  + f_0 * ik_199[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, gk_20, gk_21, gk_22, gk_23, gk_24, \
                         ik_200, ik_201, ik_202, ik_203, ik_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = -gk_20[k]
                  + f_0 * ik_200[k];

        t_93[k] = -gk_21[k]
                  + f_0 * ik_201[k];

        t_94[k] = -gk_22[k]
                  + f_0 * ik_202[k];

        t_95[k] = -gk_23[k]
                  + f_0 * ik_203[k];

        t_96[k] = -gk_24[k]
                  + f_0 * ik_204[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, gk_25, gk_26, gk_27, gk_28, gk_29, \
                         ik_205, ik_206, ik_207, ik_208, ik_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -gk_25[k]
                  + f_0 * ik_205[k];

        t_98[k] = -gk_26[k]
                  + f_0 * ik_206[k];

        t_99[k] = -gk_27[k]
                  + f_0 * ik_207[k];

        t_100[k] = -gk_28[k]
                   + f_0 * ik_208[k];

        t_101[k] = -gk_29[k]
                   + f_0 * ik_209[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, gk_30, gk_31, gk_32, gk_33, gk_34, \
                         ik_210, ik_211, ik_212, ik_213, ik_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -gk_30[k]
                   + f_0 * ik_210[k];

        t_103[k] = -gk_31[k]
                   + f_0 * ik_211[k];

        t_104[k] = -gk_32[k]
                   + f_0 * ik_212[k];

        t_105[k] = -gk_33[k]
                   + f_0 * ik_213[k];

        t_106[k] = -gk_34[k]
                   + f_0 * ik_214[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, gk_35, ik_215, \
                         ik_252, ik_253, ik_254, ik_255, ik_256, \
                         ik_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -gk_35[k]
                   + f_0 * ik_215[k];

        t_108[k] = f_0 * ik_252[k];

        t_109[k] = f_0 * ik_253[k];

        t_110[k] = f_0 * ik_254[k];

        t_111[k] = f_0 * ik_255[k];

        t_112[k] = f_0 * ik_256[k];

        t_113[k] = f_0 * ik_257[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, t_120, t_121, ik_258, \
                         ik_259, ik_260, ik_261, ik_262, ik_263, ik_264, \
                         ik_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * ik_258[k];

        t_115[k] = f_0 * ik_259[k];

        t_116[k] = f_0 * ik_260[k];

        t_117[k] = f_0 * ik_261[k];

        t_118[k] = f_0 * ik_262[k];

        t_119[k] = f_0 * ik_263[k];

        t_120[k] = f_0 * ik_264[k];

        t_121[k] = f_0 * ik_265[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, t_128, t_129, ik_266, \
                         ik_267, ik_268, ik_269, ik_270, ik_271, ik_272, \
                         ik_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * ik_266[k];

        t_123[k] = f_0 * ik_267[k];

        t_124[k] = f_0 * ik_268[k];

        t_125[k] = f_0 * ik_269[k];

        t_126[k] = f_0 * ik_270[k];

        t_127[k] = f_0 * ik_271[k];

        t_128[k] = f_0 * ik_272[k];

        t_129[k] = f_0 * ik_273[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, t_135, t_136, t_137, ik_274, \
                         ik_275, ik_276, ik_277, ik_278, ik_279, ik_280, \
                         ik_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * ik_274[k];

        t_131[k] = f_0 * ik_275[k];

        t_132[k] = f_0 * ik_276[k];

        t_133[k] = f_0 * ik_277[k];

        t_134[k] = f_0 * ik_278[k];

        t_135[k] = f_0 * ik_279[k];

        t_136[k] = f_0 * ik_280[k];

        t_137[k] = f_0 * ik_281[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, t_144, gk_36, ik_282, \
                         ik_283, ik_284, ik_285, ik_286, ik_287, \
                         ik_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * ik_282[k];

        t_139[k] = f_0 * ik_283[k];

        t_140[k] = f_0 * ik_284[k];

        t_141[k] = f_0 * ik_285[k];

        t_142[k] = f_0 * ik_286[k];

        t_143[k] = f_0 * ik_287[k];

        t_144[k] = -gk_36[k]
                   + f_0 * ik_288[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, gk_37, gk_38, gk_39, gk_40, gk_41, \
                         ik_289, ik_290, ik_291, ik_292, ik_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -gk_37[k]
                   + f_0 * ik_289[k];

        t_146[k] = -gk_38[k]
                   + f_0 * ik_290[k];

        t_147[k] = -gk_39[k]
                   + f_0 * ik_291[k];

        t_148[k] = -gk_40[k]
                   + f_0 * ik_292[k];

        t_149[k] = -gk_41[k]
                   + f_0 * ik_293[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, gk_42, gk_43, gk_44, gk_45, gk_46, \
                         ik_294, ik_295, ik_296, ik_297, ik_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -gk_42[k]
                   + f_0 * ik_294[k];

        t_151[k] = -gk_43[k]
                   + f_0 * ik_295[k];

        t_152[k] = -gk_44[k]
                   + f_0 * ik_296[k];

        t_153[k] = -gk_45[k]
                   + f_0 * ik_297[k];

        t_154[k] = -gk_46[k]
                   + f_0 * ik_298[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, gk_47, gk_48, gk_49, gk_50, gk_51, \
                         ik_299, ik_300, ik_301, ik_302, ik_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -gk_47[k]
                   + f_0 * ik_299[k];

        t_156[k] = -gk_48[k]
                   + f_0 * ik_300[k];

        t_157[k] = -gk_49[k]
                   + f_0 * ik_301[k];

        t_158[k] = -gk_50[k]
                   + f_0 * ik_302[k];

        t_159[k] = -gk_51[k]
                   + f_0 * ik_303[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, gk_52, gk_53, gk_54, gk_55, gk_56, \
                         ik_304, ik_305, ik_306, ik_307, ik_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -gk_52[k]
                   + f_0 * ik_304[k];

        t_161[k] = -gk_53[k]
                   + f_0 * ik_305[k];

        t_162[k] = -gk_54[k]
                   + f_0 * ik_306[k];

        t_163[k] = -gk_55[k]
                   + f_0 * ik_307[k];

        t_164[k] = -gk_56[k]
                   + f_0 * ik_308[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, gk_57, gk_58, gk_59, gk_60, gk_61, \
                         ik_309, ik_310, ik_311, ik_312, ik_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -gk_57[k]
                   + f_0 * ik_309[k];

        t_166[k] = -gk_58[k]
                   + f_0 * ik_310[k];

        t_167[k] = -gk_59[k]
                   + f_0 * ik_311[k];

        t_168[k] = -gk_60[k]
                   + f_0 * ik_312[k];

        t_169[k] = -gk_61[k]
                   + f_0 * ik_313[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, gk_62, gk_63, gk_64, gk_65, gk_66, \
                         ik_314, ik_315, ik_316, ik_317, ik_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -gk_62[k]
                   + f_0 * ik_314[k];

        t_171[k] = -gk_63[k]
                   + f_0 * ik_315[k];

        t_172[k] = -gk_64[k]
                   + f_0 * ik_316[k];

        t_173[k] = -gk_65[k]
                   + f_0 * ik_317[k];

        t_174[k] = -gk_66[k]
                   + f_0 * ik_318[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, gk_67, gk_68, gk_69, gk_70, gk_71, \
                         ik_319, ik_320, ik_321, ik_322, ik_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -gk_67[k]
                   + f_0 * ik_319[k];

        t_176[k] = -gk_68[k]
                   + f_0 * ik_320[k];

        t_177[k] = -gk_69[k]
                   + f_0 * ik_321[k];

        t_178[k] = -gk_70[k]
                   + f_0 * ik_322[k];

        t_179[k] = -gk_71[k]
                   + f_0 * ik_323[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, gk_72, gk_73, gk_74, gk_75, gk_76, \
                         ik_324, ik_325, ik_326, ik_327, ik_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * gk_72[k]
                   + f_0 * ik_324[k];

        t_181[k] = -2.0 * gk_73[k]
                   + f_0 * ik_325[k];

        t_182[k] = -2.0 * gk_74[k]
                   + f_0 * ik_326[k];

        t_183[k] = -2.0 * gk_75[k]
                   + f_0 * ik_327[k];

        t_184[k] = -2.0 * gk_76[k]
                   + f_0 * ik_328[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, gk_77, gk_78, gk_79, gk_80, gk_81, \
                         ik_329, ik_330, ik_331, ik_332, ik_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * gk_77[k]
                   + f_0 * ik_329[k];

        t_186[k] = -2.0 * gk_78[k]
                   + f_0 * ik_330[k];

        t_187[k] = -2.0 * gk_79[k]
                   + f_0 * ik_331[k];

        t_188[k] = -2.0 * gk_80[k]
                   + f_0 * ik_332[k];

        t_189[k] = -2.0 * gk_81[k]
                   + f_0 * ik_333[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, gk_82, gk_83, gk_84, gk_85, gk_86, \
                         ik_334, ik_335, ik_336, ik_337, ik_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * gk_82[k]
                   + f_0 * ik_334[k];

        t_191[k] = -2.0 * gk_83[k]
                   + f_0 * ik_335[k];

        t_192[k] = -2.0 * gk_84[k]
                   + f_0 * ik_336[k];

        t_193[k] = -2.0 * gk_85[k]
                   + f_0 * ik_337[k];

        t_194[k] = -2.0 * gk_86[k]
                   + f_0 * ik_338[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, gk_87, gk_88, gk_89, gk_90, gk_91, \
                         ik_339, ik_340, ik_341, ik_342, ik_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * gk_87[k]
                   + f_0 * ik_339[k];

        t_196[k] = -2.0 * gk_88[k]
                   + f_0 * ik_340[k];

        t_197[k] = -2.0 * gk_89[k]
                   + f_0 * ik_341[k];

        t_198[k] = -2.0 * gk_90[k]
                   + f_0 * ik_342[k];

        t_199[k] = -2.0 * gk_91[k]
                   + f_0 * ik_343[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, gk_92, gk_93, gk_94, gk_95, gk_96, \
                         ik_344, ik_345, ik_346, ik_347, ik_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * gk_92[k]
                   + f_0 * ik_344[k];

        t_201[k] = -2.0 * gk_93[k]
                   + f_0 * ik_345[k];

        t_202[k] = -2.0 * gk_94[k]
                   + f_0 * ik_346[k];

        t_203[k] = -2.0 * gk_95[k]
                   + f_0 * ik_347[k];

        t_204[k] = -2.0 * gk_96[k]
                   + f_0 * ik_348[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, gk_97, gk_98, gk_99, gk_100, \
                         gk_101, ik_349, ik_350, ik_351, ik_352, \
                         ik_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * gk_97[k]
                   + f_0 * ik_349[k];

        t_206[k] = -2.0 * gk_98[k]
                   + f_0 * ik_350[k];

        t_207[k] = -2.0 * gk_99[k]
                   + f_0 * ik_351[k];

        t_208[k] = -2.0 * gk_100[k]
                   + f_0 * ik_352[k];

        t_209[k] = -2.0 * gk_101[k]
                   + f_0 * ik_353[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, gk_102, gk_103, gk_104, gk_105, \
                         gk_106, ik_354, ik_355, ik_356, ik_357, \
                         ik_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * gk_102[k]
                   + f_0 * ik_354[k];

        t_211[k] = -2.0 * gk_103[k]
                   + f_0 * ik_355[k];

        t_212[k] = -2.0 * gk_104[k]
                   + f_0 * ik_356[k];

        t_213[k] = -2.0 * gk_105[k]
                   + f_0 * ik_357[k];

        t_214[k] = -2.0 * gk_106[k]
                   + f_0 * ik_358[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, t_220, t_221, gk_107, ik_359, \
                         ik_396, ik_397, ik_398, ik_399, ik_400, \
                         ik_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * gk_107[k]
                   + f_0 * ik_359[k];

        t_216[k] = f_0 * ik_396[k];

        t_217[k] = f_0 * ik_397[k];

        t_218[k] = f_0 * ik_398[k];

        t_219[k] = f_0 * ik_399[k];

        t_220[k] = f_0 * ik_400[k];

        t_221[k] = f_0 * ik_401[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, t_228, t_229, ik_402, \
                         ik_403, ik_404, ik_405, ik_406, ik_407, ik_408, \
                         ik_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_0 * ik_402[k];

        t_223[k] = f_0 * ik_403[k];

        t_224[k] = f_0 * ik_404[k];

        t_225[k] = f_0 * ik_405[k];

        t_226[k] = f_0 * ik_406[k];

        t_227[k] = f_0 * ik_407[k];

        t_228[k] = f_0 * ik_408[k];

        t_229[k] = f_0 * ik_409[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, t_235, t_236, t_237, ik_410, \
                         ik_411, ik_412, ik_413, ik_414, ik_415, ik_416, \
                         ik_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_0 * ik_410[k];

        t_231[k] = f_0 * ik_411[k];

        t_232[k] = f_0 * ik_412[k];

        t_233[k] = f_0 * ik_413[k];

        t_234[k] = f_0 * ik_414[k];

        t_235[k] = f_0 * ik_415[k];

        t_236[k] = f_0 * ik_416[k];

        t_237[k] = f_0 * ik_417[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, t_243, t_244, t_245, ik_418, \
                         ik_419, ik_420, ik_421, ik_422, ik_423, ik_424, \
                         ik_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_0 * ik_418[k];

        t_239[k] = f_0 * ik_419[k];

        t_240[k] = f_0 * ik_420[k];

        t_241[k] = f_0 * ik_421[k];

        t_242[k] = f_0 * ik_422[k];

        t_243[k] = f_0 * ik_423[k];

        t_244[k] = f_0 * ik_424[k];

        t_245[k] = f_0 * ik_425[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, t_251, t_252, gk_108, ik_426, \
                         ik_427, ik_428, ik_429, ik_430, ik_431, \
                         ik_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_0 * ik_426[k];

        t_247[k] = f_0 * ik_427[k];

        t_248[k] = f_0 * ik_428[k];

        t_249[k] = f_0 * ik_429[k];

        t_250[k] = f_0 * ik_430[k];

        t_251[k] = f_0 * ik_431[k];

        t_252[k] = -gk_108[k]
                   + f_0 * ik_432[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, gk_109, gk_110, gk_111, gk_112, \
                         gk_113, ik_433, ik_434, ik_435, ik_436, \
                         ik_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = -gk_109[k]
                   + f_0 * ik_433[k];

        t_254[k] = -gk_110[k]
                   + f_0 * ik_434[k];

        t_255[k] = -gk_111[k]
                   + f_0 * ik_435[k];

        t_256[k] = -gk_112[k]
                   + f_0 * ik_436[k];

        t_257[k] = -gk_113[k]
                   + f_0 * ik_437[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, gk_114, gk_115, gk_116, gk_117, \
                         gk_118, ik_438, ik_439, ik_440, ik_441, \
                         ik_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = -gk_114[k]
                   + f_0 * ik_438[k];

        t_259[k] = -gk_115[k]
                   + f_0 * ik_439[k];

        t_260[k] = -gk_116[k]
                   + f_0 * ik_440[k];

        t_261[k] = -gk_117[k]
                   + f_0 * ik_441[k];

        t_262[k] = -gk_118[k]
                   + f_0 * ik_442[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, gk_119, gk_120, gk_121, gk_122, \
                         gk_123, ik_443, ik_444, ik_445, ik_446, \
                         ik_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = -gk_119[k]
                   + f_0 * ik_443[k];

        t_264[k] = -gk_120[k]
                   + f_0 * ik_444[k];

        t_265[k] = -gk_121[k]
                   + f_0 * ik_445[k];

        t_266[k] = -gk_122[k]
                   + f_0 * ik_446[k];

        t_267[k] = -gk_123[k]
                   + f_0 * ik_447[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, gk_124, gk_125, gk_126, gk_127, \
                         gk_128, ik_448, ik_449, ik_450, ik_451, \
                         ik_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = -gk_124[k]
                   + f_0 * ik_448[k];

        t_269[k] = -gk_125[k]
                   + f_0 * ik_449[k];

        t_270[k] = -gk_126[k]
                   + f_0 * ik_450[k];

        t_271[k] = -gk_127[k]
                   + f_0 * ik_451[k];

        t_272[k] = -gk_128[k]
                   + f_0 * ik_452[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, gk_129, gk_130, gk_131, gk_132, \
                         gk_133, ik_453, ik_454, ik_455, ik_456, \
                         ik_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = -gk_129[k]
                   + f_0 * ik_453[k];

        t_274[k] = -gk_130[k]
                   + f_0 * ik_454[k];

        t_275[k] = -gk_131[k]
                   + f_0 * ik_455[k];

        t_276[k] = -gk_132[k]
                   + f_0 * ik_456[k];

        t_277[k] = -gk_133[k]
                   + f_0 * ik_457[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, gk_134, gk_135, gk_136, gk_137, \
                         gk_138, ik_458, ik_459, ik_460, ik_461, \
                         ik_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = -gk_134[k]
                   + f_0 * ik_458[k];

        t_279[k] = -gk_135[k]
                   + f_0 * ik_459[k];

        t_280[k] = -gk_136[k]
                   + f_0 * ik_460[k];

        t_281[k] = -gk_137[k]
                   + f_0 * ik_461[k];

        t_282[k] = -gk_138[k]
                   + f_0 * ik_462[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, gk_139, gk_140, gk_141, gk_142, \
                         gk_143, ik_463, ik_464, ik_465, ik_466, \
                         ik_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = -gk_139[k]
                   + f_0 * ik_463[k];

        t_284[k] = -gk_140[k]
                   + f_0 * ik_464[k];

        t_285[k] = -gk_141[k]
                   + f_0 * ik_465[k];

        t_286[k] = -gk_142[k]
                   + f_0 * ik_466[k];

        t_287[k] = -gk_143[k]
                   + f_0 * ik_467[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, gk_144, gk_145, gk_146, gk_147, \
                         gk_148, ik_468, ik_469, ik_470, ik_471, \
                         ik_472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = -2.0 * gk_144[k]
                   + f_0 * ik_468[k];

        t_289[k] = -2.0 * gk_145[k]
                   + f_0 * ik_469[k];

        t_290[k] = -2.0 * gk_146[k]
                   + f_0 * ik_470[k];

        t_291[k] = -2.0 * gk_147[k]
                   + f_0 * ik_471[k];

        t_292[k] = -2.0 * gk_148[k]
                   + f_0 * ik_472[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, t_297, gk_149, gk_150, gk_151, gk_152, \
                         gk_153, ik_473, ik_474, ik_475, ik_476, \
                         ik_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = -2.0 * gk_149[k]
                   + f_0 * ik_473[k];

        t_294[k] = -2.0 * gk_150[k]
                   + f_0 * ik_474[k];

        t_295[k] = -2.0 * gk_151[k]
                   + f_0 * ik_475[k];

        t_296[k] = -2.0 * gk_152[k]
                   + f_0 * ik_476[k];

        t_297[k] = -2.0 * gk_153[k]
                   + f_0 * ik_477[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, gk_154, gk_155, gk_156, gk_157, \
                         gk_158, ik_478, ik_479, ik_480, ik_481, \
                         ik_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = -2.0 * gk_154[k]
                   + f_0 * ik_478[k];

        t_299[k] = -2.0 * gk_155[k]
                   + f_0 * ik_479[k];

        t_300[k] = -2.0 * gk_156[k]
                   + f_0 * ik_480[k];

        t_301[k] = -2.0 * gk_157[k]
                   + f_0 * ik_481[k];

        t_302[k] = -2.0 * gk_158[k]
                   + f_0 * ik_482[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, gk_159, gk_160, gk_161, gk_162, \
                         gk_163, ik_483, ik_484, ik_485, ik_486, \
                         ik_487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = -2.0 * gk_159[k]
                   + f_0 * ik_483[k];

        t_304[k] = -2.0 * gk_160[k]
                   + f_0 * ik_484[k];

        t_305[k] = -2.0 * gk_161[k]
                   + f_0 * ik_485[k];

        t_306[k] = -2.0 * gk_162[k]
                   + f_0 * ik_486[k];

        t_307[k] = -2.0 * gk_163[k]
                   + f_0 * ik_487[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, gk_164, gk_165, gk_166, gk_167, \
                         gk_168, ik_488, ik_489, ik_490, ik_491, \
                         ik_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = -2.0 * gk_164[k]
                   + f_0 * ik_488[k];

        t_309[k] = -2.0 * gk_165[k]
                   + f_0 * ik_489[k];

        t_310[k] = -2.0 * gk_166[k]
                   + f_0 * ik_490[k];

        t_311[k] = -2.0 * gk_167[k]
                   + f_0 * ik_491[k];

        t_312[k] = -2.0 * gk_168[k]
                   + f_0 * ik_492[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, gk_169, gk_170, gk_171, gk_172, \
                         gk_173, ik_493, ik_494, ik_495, ik_496, \
                         ik_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = -2.0 * gk_169[k]
                   + f_0 * ik_493[k];

        t_314[k] = -2.0 * gk_170[k]
                   + f_0 * ik_494[k];

        t_315[k] = -2.0 * gk_171[k]
                   + f_0 * ik_495[k];

        t_316[k] = -2.0 * gk_172[k]
                   + f_0 * ik_496[k];

        t_317[k] = -2.0 * gk_173[k]
                   + f_0 * ik_497[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, gk_174, gk_175, gk_176, gk_177, \
                         gk_178, ik_498, ik_499, ik_500, ik_501, \
                         ik_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = -2.0 * gk_174[k]
                   + f_0 * ik_498[k];

        t_319[k] = -2.0 * gk_175[k]
                   + f_0 * ik_499[k];

        t_320[k] = -2.0 * gk_176[k]
                   + f_0 * ik_500[k];

        t_321[k] = -2.0 * gk_177[k]
                   + f_0 * ik_501[k];

        t_322[k] = -2.0 * gk_178[k]
                   + f_0 * ik_502[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, gk_179, gk_180, gk_181, gk_182, \
                         gk_183, ik_503, ik_504, ik_505, ik_506, \
                         ik_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = -2.0 * gk_179[k]
                   + f_0 * ik_503[k];

        t_324[k] = -3.0 * gk_180[k]
                   + f_0 * ik_504[k];

        t_325[k] = -3.0 * gk_181[k]
                   + f_0 * ik_505[k];

        t_326[k] = -3.0 * gk_182[k]
                   + f_0 * ik_506[k];

        t_327[k] = -3.0 * gk_183[k]
                   + f_0 * ik_507[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, gk_184, gk_185, gk_186, gk_187, \
                         gk_188, ik_508, ik_509, ik_510, ik_511, \
                         ik_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = -3.0 * gk_184[k]
                   + f_0 * ik_508[k];

        t_329[k] = -3.0 * gk_185[k]
                   + f_0 * ik_509[k];

        t_330[k] = -3.0 * gk_186[k]
                   + f_0 * ik_510[k];

        t_331[k] = -3.0 * gk_187[k]
                   + f_0 * ik_511[k];

        t_332[k] = -3.0 * gk_188[k]
                   + f_0 * ik_512[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, t_337, gk_189, gk_190, gk_191, gk_192, \
                         gk_193, ik_513, ik_514, ik_515, ik_516, \
                         ik_517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = -3.0 * gk_189[k]
                   + f_0 * ik_513[k];

        t_334[k] = -3.0 * gk_190[k]
                   + f_0 * ik_514[k];

        t_335[k] = -3.0 * gk_191[k]
                   + f_0 * ik_515[k];

        t_336[k] = -3.0 * gk_192[k]
                   + f_0 * ik_516[k];

        t_337[k] = -3.0 * gk_193[k]
                   + f_0 * ik_517[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, gk_194, gk_195, gk_196, gk_197, \
                         gk_198, ik_518, ik_519, ik_520, ik_521, \
                         ik_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = -3.0 * gk_194[k]
                   + f_0 * ik_518[k];

        t_339[k] = -3.0 * gk_195[k]
                   + f_0 * ik_519[k];

        t_340[k] = -3.0 * gk_196[k]
                   + f_0 * ik_520[k];

        t_341[k] = -3.0 * gk_197[k]
                   + f_0 * ik_521[k];

        t_342[k] = -3.0 * gk_198[k]
                   + f_0 * ik_522[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, gk_199, gk_200, gk_201, gk_202, \
                         gk_203, ik_523, ik_524, ik_525, ik_526, \
                         ik_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = -3.0 * gk_199[k]
                   + f_0 * ik_523[k];

        t_344[k] = -3.0 * gk_200[k]
                   + f_0 * ik_524[k];

        t_345[k] = -3.0 * gk_201[k]
                   + f_0 * ik_525[k];

        t_346[k] = -3.0 * gk_202[k]
                   + f_0 * ik_526[k];

        t_347[k] = -3.0 * gk_203[k]
                   + f_0 * ik_527[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, gk_204, gk_205, gk_206, gk_207, \
                         gk_208, ik_528, ik_529, ik_530, ik_531, \
                         ik_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = -3.0 * gk_204[k]
                   + f_0 * ik_528[k];

        t_349[k] = -3.0 * gk_205[k]
                   + f_0 * ik_529[k];

        t_350[k] = -3.0 * gk_206[k]
                   + f_0 * ik_530[k];

        t_351[k] = -3.0 * gk_207[k]
                   + f_0 * ik_531[k];

        t_352[k] = -3.0 * gk_208[k]
                   + f_0 * ik_532[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, gk_209, gk_210, gk_211, gk_212, \
                         gk_213, ik_533, ik_534, ik_535, ik_536, \
                         ik_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = -3.0 * gk_209[k]
                   + f_0 * ik_533[k];

        t_354[k] = -3.0 * gk_210[k]
                   + f_0 * ik_534[k];

        t_355[k] = -3.0 * gk_211[k]
                   + f_0 * ik_535[k];

        t_356[k] = -3.0 * gk_212[k]
                   + f_0 * ik_536[k];

        t_357[k] = -3.0 * gk_213[k]
                   + f_0 * ik_537[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, t_363, t_364, gk_214, gk_215, \
                         ik_538, ik_539, ik_576, ik_577, ik_578, ik_579, \
                         ik_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = -3.0 * gk_214[k]
                   + f_0 * ik_538[k];

        t_359[k] = -3.0 * gk_215[k]
                   + f_0 * ik_539[k];

        t_360[k] = f_0 * ik_576[k];

        t_361[k] = f_0 * ik_577[k];

        t_362[k] = f_0 * ik_578[k];

        t_363[k] = f_0 * ik_579[k];

        t_364[k] = f_0 * ik_580[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, t_370, t_371, t_372, ik_581, \
                         ik_582, ik_583, ik_584, ik_585, ik_586, ik_587, \
                         ik_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_0 * ik_581[k];

        t_366[k] = f_0 * ik_582[k];

        t_367[k] = f_0 * ik_583[k];

        t_368[k] = f_0 * ik_584[k];

        t_369[k] = f_0 * ik_585[k];

        t_370[k] = f_0 * ik_586[k];

        t_371[k] = f_0 * ik_587[k];

        t_372[k] = f_0 * ik_588[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, t_378, t_379, t_380, ik_589, \
                         ik_590, ik_591, ik_592, ik_593, ik_594, ik_595, \
                         ik_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_0 * ik_589[k];

        t_374[k] = f_0 * ik_590[k];

        t_375[k] = f_0 * ik_591[k];

        t_376[k] = f_0 * ik_592[k];

        t_377[k] = f_0 * ik_593[k];

        t_378[k] = f_0 * ik_594[k];

        t_379[k] = f_0 * ik_595[k];

        t_380[k] = f_0 * ik_596[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, t_385, t_386, t_387, t_388, ik_597, \
                         ik_598, ik_599, ik_600, ik_601, ik_602, ik_603, \
                         ik_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_0 * ik_597[k];

        t_382[k] = f_0 * ik_598[k];

        t_383[k] = f_0 * ik_599[k];

        t_384[k] = f_0 * ik_600[k];

        t_385[k] = f_0 * ik_601[k];

        t_386[k] = f_0 * ik_602[k];

        t_387[k] = f_0 * ik_603[k];

        t_388[k] = f_0 * ik_604[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, t_394, t_395, ik_605, ik_606, \
                         ik_607, ik_608, ik_609, ik_610, ik_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_0 * ik_605[k];

        t_390[k] = f_0 * ik_606[k];

        t_391[k] = f_0 * ik_607[k];

        t_392[k] = f_0 * ik_608[k];

        t_393[k] = f_0 * ik_609[k];

        t_394[k] = f_0 * ik_610[k];

        t_395[k] = f_0 * ik_611[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, t_400, gk_216, gk_217, gk_218, gk_219, \
                         gk_220, ik_612, ik_613, ik_614, ik_615, \
                         ik_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = -gk_216[k]
                   + f_0 * ik_612[k];

        t_397[k] = -gk_217[k]
                   + f_0 * ik_613[k];

        t_398[k] = -gk_218[k]
                   + f_0 * ik_614[k];

        t_399[k] = -gk_219[k]
                   + f_0 * ik_615[k];

        t_400[k] = -gk_220[k]
                   + f_0 * ik_616[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, t_405, gk_221, gk_222, gk_223, gk_224, \
                         gk_225, ik_617, ik_618, ik_619, ik_620, \
                         ik_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = -gk_221[k]
                   + f_0 * ik_617[k];

        t_402[k] = -gk_222[k]
                   + f_0 * ik_618[k];

        t_403[k] = -gk_223[k]
                   + f_0 * ik_619[k];

        t_404[k] = -gk_224[k]
                   + f_0 * ik_620[k];

        t_405[k] = -gk_225[k]
                   + f_0 * ik_621[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, t_410, gk_226, gk_227, gk_228, gk_229, \
                         gk_230, ik_622, ik_623, ik_624, ik_625, \
                         ik_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = -gk_226[k]
                   + f_0 * ik_622[k];

        t_407[k] = -gk_227[k]
                   + f_0 * ik_623[k];

        t_408[k] = -gk_228[k]
                   + f_0 * ik_624[k];

        t_409[k] = -gk_229[k]
                   + f_0 * ik_625[k];

        t_410[k] = -gk_230[k]
                   + f_0 * ik_626[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, t_415, gk_231, gk_232, gk_233, gk_234, \
                         gk_235, ik_627, ik_628, ik_629, ik_630, \
                         ik_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = -gk_231[k]
                   + f_0 * ik_627[k];

        t_412[k] = -gk_232[k]
                   + f_0 * ik_628[k];

        t_413[k] = -gk_233[k]
                   + f_0 * ik_629[k];

        t_414[k] = -gk_234[k]
                   + f_0 * ik_630[k];

        t_415[k] = -gk_235[k]
                   + f_0 * ik_631[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, t_420, gk_236, gk_237, gk_238, gk_239, \
                         gk_240, ik_632, ik_633, ik_634, ik_635, \
                         ik_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = -gk_236[k]
                   + f_0 * ik_632[k];

        t_417[k] = -gk_237[k]
                   + f_0 * ik_633[k];

        t_418[k] = -gk_238[k]
                   + f_0 * ik_634[k];

        t_419[k] = -gk_239[k]
                   + f_0 * ik_635[k];

        t_420[k] = -gk_240[k]
                   + f_0 * ik_636[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, t_425, gk_241, gk_242, gk_243, gk_244, \
                         gk_245, ik_637, ik_638, ik_639, ik_640, \
                         ik_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = -gk_241[k]
                   + f_0 * ik_637[k];

        t_422[k] = -gk_242[k]
                   + f_0 * ik_638[k];

        t_423[k] = -gk_243[k]
                   + f_0 * ik_639[k];

        t_424[k] = -gk_244[k]
                   + f_0 * ik_640[k];

        t_425[k] = -gk_245[k]
                   + f_0 * ik_641[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, t_430, gk_246, gk_247, gk_248, gk_249, \
                         gk_250, ik_642, ik_643, ik_644, ik_645, \
                         ik_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = -gk_246[k]
                   + f_0 * ik_642[k];

        t_427[k] = -gk_247[k]
                   + f_0 * ik_643[k];

        t_428[k] = -gk_248[k]
                   + f_0 * ik_644[k];

        t_429[k] = -gk_249[k]
                   + f_0 * ik_645[k];

        t_430[k] = -gk_250[k]
                   + f_0 * ik_646[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, gk_251, gk_252, gk_253, gk_254, \
                         gk_255, ik_647, ik_648, ik_649, ik_650, \
                         ik_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = -gk_251[k]
                   + f_0 * ik_647[k];

        t_432[k] = -2.0 * gk_252[k]
                   + f_0 * ik_648[k];

        t_433[k] = -2.0 * gk_253[k]
                   + f_0 * ik_649[k];

        t_434[k] = -2.0 * gk_254[k]
                   + f_0 * ik_650[k];

        t_435[k] = -2.0 * gk_255[k]
                   + f_0 * ik_651[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, gk_256, gk_257, gk_258, gk_259, \
                         gk_260, ik_652, ik_653, ik_654, ik_655, \
                         ik_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = -2.0 * gk_256[k]
                   + f_0 * ik_652[k];

        t_437[k] = -2.0 * gk_257[k]
                   + f_0 * ik_653[k];

        t_438[k] = -2.0 * gk_258[k]
                   + f_0 * ik_654[k];

        t_439[k] = -2.0 * gk_259[k]
                   + f_0 * ik_655[k];

        t_440[k] = -2.0 * gk_260[k]
                   + f_0 * ik_656[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, t_445, gk_261, gk_262, gk_263, gk_264, \
                         gk_265, ik_657, ik_658, ik_659, ik_660, \
                         ik_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = -2.0 * gk_261[k]
                   + f_0 * ik_657[k];

        t_442[k] = -2.0 * gk_262[k]
                   + f_0 * ik_658[k];

        t_443[k] = -2.0 * gk_263[k]
                   + f_0 * ik_659[k];

        t_444[k] = -2.0 * gk_264[k]
                   + f_0 * ik_660[k];

        t_445[k] = -2.0 * gk_265[k]
                   + f_0 * ik_661[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, t_450, gk_266, gk_267, gk_268, gk_269, \
                         gk_270, ik_662, ik_663, ik_664, ik_665, \
                         ik_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = -2.0 * gk_266[k]
                   + f_0 * ik_662[k];

        t_447[k] = -2.0 * gk_267[k]
                   + f_0 * ik_663[k];

        t_448[k] = -2.0 * gk_268[k]
                   + f_0 * ik_664[k];

        t_449[k] = -2.0 * gk_269[k]
                   + f_0 * ik_665[k];

        t_450[k] = -2.0 * gk_270[k]
                   + f_0 * ik_666[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, t_455, gk_271, gk_272, gk_273, gk_274, \
                         gk_275, ik_667, ik_668, ik_669, ik_670, \
                         ik_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = -2.0 * gk_271[k]
                   + f_0 * ik_667[k];

        t_452[k] = -2.0 * gk_272[k]
                   + f_0 * ik_668[k];

        t_453[k] = -2.0 * gk_273[k]
                   + f_0 * ik_669[k];

        t_454[k] = -2.0 * gk_274[k]
                   + f_0 * ik_670[k];

        t_455[k] = -2.0 * gk_275[k]
                   + f_0 * ik_671[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, t_460, gk_276, gk_277, gk_278, gk_279, \
                         gk_280, ik_672, ik_673, ik_674, ik_675, \
                         ik_676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = -2.0 * gk_276[k]
                   + f_0 * ik_672[k];

        t_457[k] = -2.0 * gk_277[k]
                   + f_0 * ik_673[k];

        t_458[k] = -2.0 * gk_278[k]
                   + f_0 * ik_674[k];

        t_459[k] = -2.0 * gk_279[k]
                   + f_0 * ik_675[k];

        t_460[k] = -2.0 * gk_280[k]
                   + f_0 * ik_676[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, t_465, gk_281, gk_282, gk_283, gk_284, \
                         gk_285, ik_677, ik_678, ik_679, ik_680, \
                         ik_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = -2.0 * gk_281[k]
                   + f_0 * ik_677[k];

        t_462[k] = -2.0 * gk_282[k]
                   + f_0 * ik_678[k];

        t_463[k] = -2.0 * gk_283[k]
                   + f_0 * ik_679[k];

        t_464[k] = -2.0 * gk_284[k]
                   + f_0 * ik_680[k];

        t_465[k] = -2.0 * gk_285[k]
                   + f_0 * ik_681[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, gk_286, gk_287, gk_288, gk_289, \
                         gk_290, ik_682, ik_683, ik_684, ik_685, \
                         ik_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = -2.0 * gk_286[k]
                   + f_0 * ik_682[k];

        t_467[k] = -2.0 * gk_287[k]
                   + f_0 * ik_683[k];

        t_468[k] = -3.0 * gk_288[k]
                   + f_0 * ik_684[k];

        t_469[k] = -3.0 * gk_289[k]
                   + f_0 * ik_685[k];

        t_470[k] = -3.0 * gk_290[k]
                   + f_0 * ik_686[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, t_475, gk_291, gk_292, gk_293, gk_294, \
                         gk_295, ik_687, ik_688, ik_689, ik_690, \
                         ik_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = -3.0 * gk_291[k]
                   + f_0 * ik_687[k];

        t_472[k] = -3.0 * gk_292[k]
                   + f_0 * ik_688[k];

        t_473[k] = -3.0 * gk_293[k]
                   + f_0 * ik_689[k];

        t_474[k] = -3.0 * gk_294[k]
                   + f_0 * ik_690[k];

        t_475[k] = -3.0 * gk_295[k]
                   + f_0 * ik_691[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, gk_296, gk_297, gk_298, gk_299, \
                         gk_300, ik_692, ik_693, ik_694, ik_695, \
                         ik_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = -3.0 * gk_296[k]
                   + f_0 * ik_692[k];

        t_477[k] = -3.0 * gk_297[k]
                   + f_0 * ik_693[k];

        t_478[k] = -3.0 * gk_298[k]
                   + f_0 * ik_694[k];

        t_479[k] = -3.0 * gk_299[k]
                   + f_0 * ik_695[k];

        t_480[k] = -3.0 * gk_300[k]
                   + f_0 * ik_696[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, gk_301, gk_302, gk_303, gk_304, \
                         gk_305, ik_697, ik_698, ik_699, ik_700, \
                         ik_701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = -3.0 * gk_301[k]
                   + f_0 * ik_697[k];

        t_482[k] = -3.0 * gk_302[k]
                   + f_0 * ik_698[k];

        t_483[k] = -3.0 * gk_303[k]
                   + f_0 * ik_699[k];

        t_484[k] = -3.0 * gk_304[k]
                   + f_0 * ik_700[k];

        t_485[k] = -3.0 * gk_305[k]
                   + f_0 * ik_701[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, gk_306, gk_307, gk_308, gk_309, \
                         gk_310, ik_702, ik_703, ik_704, ik_705, \
                         ik_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = -3.0 * gk_306[k]
                   + f_0 * ik_702[k];

        t_487[k] = -3.0 * gk_307[k]
                   + f_0 * ik_703[k];

        t_488[k] = -3.0 * gk_308[k]
                   + f_0 * ik_704[k];

        t_489[k] = -3.0 * gk_309[k]
                   + f_0 * ik_705[k];

        t_490[k] = -3.0 * gk_310[k]
                   + f_0 * ik_706[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, gk_311, gk_312, gk_313, gk_314, \
                         gk_315, ik_707, ik_708, ik_709, ik_710, \
                         ik_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = -3.0 * gk_311[k]
                   + f_0 * ik_707[k];

        t_492[k] = -3.0 * gk_312[k]
                   + f_0 * ik_708[k];

        t_493[k] = -3.0 * gk_313[k]
                   + f_0 * ik_709[k];

        t_494[k] = -3.0 * gk_314[k]
                   + f_0 * ik_710[k];

        t_495[k] = -3.0 * gk_315[k]
                   + f_0 * ik_711[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, t_500, gk_316, gk_317, gk_318, gk_319, \
                         gk_320, ik_712, ik_713, ik_714, ik_715, \
                         ik_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = -3.0 * gk_316[k]
                   + f_0 * ik_712[k];

        t_497[k] = -3.0 * gk_317[k]
                   + f_0 * ik_713[k];

        t_498[k] = -3.0 * gk_318[k]
                   + f_0 * ik_714[k];

        t_499[k] = -3.0 * gk_319[k]
                   + f_0 * ik_715[k];

        t_500[k] = -3.0 * gk_320[k]
                   + f_0 * ik_716[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, t_505, gk_321, gk_322, gk_323, gk_324, \
                         gk_325, ik_717, ik_718, ik_719, ik_720, \
                         ik_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = -3.0 * gk_321[k]
                   + f_0 * ik_717[k];

        t_502[k] = -3.0 * gk_322[k]
                   + f_0 * ik_718[k];

        t_503[k] = -3.0 * gk_323[k]
                   + f_0 * ik_719[k];

        t_504[k] = -4.0 * gk_324[k]
                   + f_0 * ik_720[k];

        t_505[k] = -4.0 * gk_325[k]
                   + f_0 * ik_721[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, t_509, t_510, gk_326, gk_327, gk_328, gk_329, \
                         gk_330, ik_722, ik_723, ik_724, ik_725, \
                         ik_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = -4.0 * gk_326[k]
                   + f_0 * ik_722[k];

        t_507[k] = -4.0 * gk_327[k]
                   + f_0 * ik_723[k];

        t_508[k] = -4.0 * gk_328[k]
                   + f_0 * ik_724[k];

        t_509[k] = -4.0 * gk_329[k]
                   + f_0 * ik_725[k];

        t_510[k] = -4.0 * gk_330[k]
                   + f_0 * ik_726[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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
    const auto *ik_754 = buffer.data(ik + 754);
    const auto *ik_755 = buffer.data(ik + 755);
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

#pragma omp simd aligned(t_511, t_512, t_513, t_514, t_515, gk_331, gk_332, gk_333, gk_334, \
                         gk_335, ik_727, ik_728, ik_729, ik_730, \
                         ik_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = -4.0 * gk_331[k]
                   + f_0 * ik_727[k];

        t_512[k] = -4.0 * gk_332[k]
                   + f_0 * ik_728[k];

        t_513[k] = -4.0 * gk_333[k]
                   + f_0 * ik_729[k];

        t_514[k] = -4.0 * gk_334[k]
                   + f_0 * ik_730[k];

        t_515[k] = -4.0 * gk_335[k]
                   + f_0 * ik_731[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, t_520, gk_336, gk_337, gk_338, gk_339, \
                         gk_340, ik_732, ik_733, ik_734, ik_735, \
                         ik_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = -4.0 * gk_336[k]
                   + f_0 * ik_732[k];

        t_517[k] = -4.0 * gk_337[k]
                   + f_0 * ik_733[k];

        t_518[k] = -4.0 * gk_338[k]
                   + f_0 * ik_734[k];

        t_519[k] = -4.0 * gk_339[k]
                   + f_0 * ik_735[k];

        t_520[k] = -4.0 * gk_340[k]
                   + f_0 * ik_736[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, t_525, gk_341, gk_342, gk_343, gk_344, \
                         gk_345, ik_737, ik_738, ik_739, ik_740, \
                         ik_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = -4.0 * gk_341[k]
                   + f_0 * ik_737[k];

        t_522[k] = -4.0 * gk_342[k]
                   + f_0 * ik_738[k];

        t_523[k] = -4.0 * gk_343[k]
                   + f_0 * ik_739[k];

        t_524[k] = -4.0 * gk_344[k]
                   + f_0 * ik_740[k];

        t_525[k] = -4.0 * gk_345[k]
                   + f_0 * ik_741[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, gk_346, gk_347, gk_348, gk_349, \
                         gk_350, ik_742, ik_743, ik_744, ik_745, \
                         ik_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = -4.0 * gk_346[k]
                   + f_0 * ik_742[k];

        t_527[k] = -4.0 * gk_347[k]
                   + f_0 * ik_743[k];

        t_528[k] = -4.0 * gk_348[k]
                   + f_0 * ik_744[k];

        t_529[k] = -4.0 * gk_349[k]
                   + f_0 * ik_745[k];

        t_530[k] = -4.0 * gk_350[k]
                   + f_0 * ik_746[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, t_535, gk_351, gk_352, gk_353, gk_354, \
                         gk_355, ik_747, ik_748, ik_749, ik_750, \
                         ik_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = -4.0 * gk_351[k]
                   + f_0 * ik_747[k];

        t_532[k] = -4.0 * gk_352[k]
                   + f_0 * ik_748[k];

        t_533[k] = -4.0 * gk_353[k]
                   + f_0 * ik_749[k];

        t_534[k] = -4.0 * gk_354[k]
                   + f_0 * ik_750[k];

        t_535[k] = -4.0 * gk_355[k]
                   + f_0 * ik_751[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, t_541, gk_356, gk_357, gk_358, \
                         gk_359, ik_752, ik_753, ik_754, ik_755, ik_792, \
                         ik_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = -4.0 * gk_356[k]
                   + f_0 * ik_752[k];

        t_537[k] = -4.0 * gk_357[k]
                   + f_0 * ik_753[k];

        t_538[k] = -4.0 * gk_358[k]
                   + f_0 * ik_754[k];

        t_539[k] = -4.0 * gk_359[k]
                   + f_0 * ik_755[k];

        t_540[k] = f_0 * ik_792[k];

        t_541[k] = f_0 * ik_793[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, t_547, t_548, t_549, ik_794, \
                         ik_795, ik_796, ik_797, ik_798, ik_799, ik_800, \
                         ik_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_0 * ik_794[k];

        t_543[k] = f_0 * ik_795[k];

        t_544[k] = f_0 * ik_796[k];

        t_545[k] = f_0 * ik_797[k];

        t_546[k] = f_0 * ik_798[k];

        t_547[k] = f_0 * ik_799[k];

        t_548[k] = f_0 * ik_800[k];

        t_549[k] = f_0 * ik_801[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, t_555, t_556, t_557, ik_802, \
                         ik_803, ik_804, ik_805, ik_806, ik_807, ik_808, \
                         ik_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_0 * ik_802[k];

        t_551[k] = f_0 * ik_803[k];

        t_552[k] = f_0 * ik_804[k];

        t_553[k] = f_0 * ik_805[k];

        t_554[k] = f_0 * ik_806[k];

        t_555[k] = f_0 * ik_807[k];

        t_556[k] = f_0 * ik_808[k];

        t_557[k] = f_0 * ik_809[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, t_562, t_563, t_564, t_565, ik_810, \
                         ik_811, ik_812, ik_813, ik_814, ik_815, ik_816, \
                         ik_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_0 * ik_810[k];

        t_559[k] = f_0 * ik_811[k];

        t_560[k] = f_0 * ik_812[k];

        t_561[k] = f_0 * ik_813[k];

        t_562[k] = f_0 * ik_814[k];

        t_563[k] = f_0 * ik_815[k];

        t_564[k] = f_0 * ik_816[k];

        t_565[k] = f_0 * ik_817[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, t_570, t_571, t_572, t_573, ik_818, \
                         ik_819, ik_820, ik_821, ik_822, ik_823, ik_824, \
                         ik_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_0 * ik_818[k];

        t_567[k] = f_0 * ik_819[k];

        t_568[k] = f_0 * ik_820[k];

        t_569[k] = f_0 * ik_821[k];

        t_570[k] = f_0 * ik_822[k];

        t_571[k] = f_0 * ik_823[k];

        t_572[k] = f_0 * ik_824[k];

        t_573[k] = f_0 * ik_825[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, t_579, gk_360, gk_361, gk_362, \
                         gk_363, ik_826, ik_827, ik_828, ik_829, ik_830, \
                         ik_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_0 * ik_826[k];

        t_575[k] = f_0 * ik_827[k];

        t_576[k] = -gk_360[k]
                   + f_0 * ik_828[k];

        t_577[k] = -gk_361[k]
                   + f_0 * ik_829[k];

        t_578[k] = -gk_362[k]
                   + f_0 * ik_830[k];

        t_579[k] = -gk_363[k]
                   + f_0 * ik_831[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, gk_364, gk_365, gk_366, gk_367, \
                         gk_368, ik_832, ik_833, ik_834, ik_835, \
                         ik_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -gk_364[k]
                   + f_0 * ik_832[k];

        t_581[k] = -gk_365[k]
                   + f_0 * ik_833[k];

        t_582[k] = -gk_366[k]
                   + f_0 * ik_834[k];

        t_583[k] = -gk_367[k]
                   + f_0 * ik_835[k];

        t_584[k] = -gk_368[k]
                   + f_0 * ik_836[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, gk_369, gk_370, gk_371, gk_372, \
                         gk_373, ik_837, ik_838, ik_839, ik_840, \
                         ik_841 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -gk_369[k]
                   + f_0 * ik_837[k];

        t_586[k] = -gk_370[k]
                   + f_0 * ik_838[k];

        t_587[k] = -gk_371[k]
                   + f_0 * ik_839[k];

        t_588[k] = -gk_372[k]
                   + f_0 * ik_840[k];

        t_589[k] = -gk_373[k]
                   + f_0 * ik_841[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, gk_374, gk_375, gk_376, gk_377, \
                         gk_378, ik_842, ik_843, ik_844, ik_845, \
                         ik_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -gk_374[k]
                   + f_0 * ik_842[k];

        t_591[k] = -gk_375[k]
                   + f_0 * ik_843[k];

        t_592[k] = -gk_376[k]
                   + f_0 * ik_844[k];

        t_593[k] = -gk_377[k]
                   + f_0 * ik_845[k];

        t_594[k] = -gk_378[k]
                   + f_0 * ik_846[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, gk_379, gk_380, gk_381, gk_382, \
                         gk_383, ik_847, ik_848, ik_849, ik_850, \
                         ik_851 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -gk_379[k]
                   + f_0 * ik_847[k];

        t_596[k] = -gk_380[k]
                   + f_0 * ik_848[k];

        t_597[k] = -gk_381[k]
                   + f_0 * ik_849[k];

        t_598[k] = -gk_382[k]
                   + f_0 * ik_850[k];

        t_599[k] = -gk_383[k]
                   + f_0 * ik_851[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, gk_384, gk_385, gk_386, gk_387, \
                         gk_388, ik_852, ik_853, ik_854, ik_855, \
                         ik_856 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -gk_384[k]
                   + f_0 * ik_852[k];

        t_601[k] = -gk_385[k]
                   + f_0 * ik_853[k];

        t_602[k] = -gk_386[k]
                   + f_0 * ik_854[k];

        t_603[k] = -gk_387[k]
                   + f_0 * ik_855[k];

        t_604[k] = -gk_388[k]
                   + f_0 * ik_856[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, gk_389, gk_390, gk_391, gk_392, \
                         gk_393, ik_857, ik_858, ik_859, ik_860, \
                         ik_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -gk_389[k]
                   + f_0 * ik_857[k];

        t_606[k] = -gk_390[k]
                   + f_0 * ik_858[k];

        t_607[k] = -gk_391[k]
                   + f_0 * ik_859[k];

        t_608[k] = -gk_392[k]
                   + f_0 * ik_860[k];

        t_609[k] = -gk_393[k]
                   + f_0 * ik_861[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, gk_394, gk_395, gk_396, gk_397, \
                         gk_398, ik_862, ik_863, ik_864, ik_865, \
                         ik_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -gk_394[k]
                   + f_0 * ik_862[k];

        t_611[k] = -gk_395[k]
                   + f_0 * ik_863[k];

        t_612[k] = -2.0 * gk_396[k]
                   + f_0 * ik_864[k];

        t_613[k] = -2.0 * gk_397[k]
                   + f_0 * ik_865[k];

        t_614[k] = -2.0 * gk_398[k]
                   + f_0 * ik_866[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, gk_399, gk_400, gk_401, gk_402, \
                         gk_403, ik_867, ik_868, ik_869, ik_870, \
                         ik_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -2.0 * gk_399[k]
                   + f_0 * ik_867[k];

        t_616[k] = -2.0 * gk_400[k]
                   + f_0 * ik_868[k];

        t_617[k] = -2.0 * gk_401[k]
                   + f_0 * ik_869[k];

        t_618[k] = -2.0 * gk_402[k]
                   + f_0 * ik_870[k];

        t_619[k] = -2.0 * gk_403[k]
                   + f_0 * ik_871[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, gk_404, gk_405, gk_406, gk_407, \
                         gk_408, ik_872, ik_873, ik_874, ik_875, \
                         ik_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -2.0 * gk_404[k]
                   + f_0 * ik_872[k];

        t_621[k] = -2.0 * gk_405[k]
                   + f_0 * ik_873[k];

        t_622[k] = -2.0 * gk_406[k]
                   + f_0 * ik_874[k];

        t_623[k] = -2.0 * gk_407[k]
                   + f_0 * ik_875[k];

        t_624[k] = -2.0 * gk_408[k]
                   + f_0 * ik_876[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, gk_409, gk_410, gk_411, gk_412, \
                         gk_413, ik_877, ik_878, ik_879, ik_880, \
                         ik_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -2.0 * gk_409[k]
                   + f_0 * ik_877[k];

        t_626[k] = -2.0 * gk_410[k]
                   + f_0 * ik_878[k];

        t_627[k] = -2.0 * gk_411[k]
                   + f_0 * ik_879[k];

        t_628[k] = -2.0 * gk_412[k]
                   + f_0 * ik_880[k];

        t_629[k] = -2.0 * gk_413[k]
                   + f_0 * ik_881[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, gk_414, gk_415, gk_416, gk_417, \
                         gk_418, ik_882, ik_883, ik_884, ik_885, \
                         ik_886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -2.0 * gk_414[k]
                   + f_0 * ik_882[k];

        t_631[k] = -2.0 * gk_415[k]
                   + f_0 * ik_883[k];

        t_632[k] = -2.0 * gk_416[k]
                   + f_0 * ik_884[k];

        t_633[k] = -2.0 * gk_417[k]
                   + f_0 * ik_885[k];

        t_634[k] = -2.0 * gk_418[k]
                   + f_0 * ik_886[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, gk_419, gk_420, gk_421, gk_422, \
                         gk_423, ik_887, ik_888, ik_889, ik_890, \
                         ik_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -2.0 * gk_419[k]
                   + f_0 * ik_887[k];

        t_636[k] = -2.0 * gk_420[k]
                   + f_0 * ik_888[k];

        t_637[k] = -2.0 * gk_421[k]
                   + f_0 * ik_889[k];

        t_638[k] = -2.0 * gk_422[k]
                   + f_0 * ik_890[k];

        t_639[k] = -2.0 * gk_423[k]
                   + f_0 * ik_891[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, gk_424, gk_425, gk_426, gk_427, \
                         gk_428, ik_892, ik_893, ik_894, ik_895, \
                         ik_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -2.0 * gk_424[k]
                   + f_0 * ik_892[k];

        t_641[k] = -2.0 * gk_425[k]
                   + f_0 * ik_893[k];

        t_642[k] = -2.0 * gk_426[k]
                   + f_0 * ik_894[k];

        t_643[k] = -2.0 * gk_427[k]
                   + f_0 * ik_895[k];

        t_644[k] = -2.0 * gk_428[k]
                   + f_0 * ik_896[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, gk_429, gk_430, gk_431, gk_432, \
                         gk_433, ik_897, ik_898, ik_899, ik_900, \
                         ik_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -2.0 * gk_429[k]
                   + f_0 * ik_897[k];

        t_646[k] = -2.0 * gk_430[k]
                   + f_0 * ik_898[k];

        t_647[k] = -2.0 * gk_431[k]
                   + f_0 * ik_899[k];

        t_648[k] = -3.0 * gk_432[k]
                   + f_0 * ik_900[k];

        t_649[k] = -3.0 * gk_433[k]
                   + f_0 * ik_901[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, gk_434, gk_435, gk_436, gk_437, \
                         gk_438, ik_902, ik_903, ik_904, ik_905, \
                         ik_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -3.0 * gk_434[k]
                   + f_0 * ik_902[k];

        t_651[k] = -3.0 * gk_435[k]
                   + f_0 * ik_903[k];

        t_652[k] = -3.0 * gk_436[k]
                   + f_0 * ik_904[k];

        t_653[k] = -3.0 * gk_437[k]
                   + f_0 * ik_905[k];

        t_654[k] = -3.0 * gk_438[k]
                   + f_0 * ik_906[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, gk_439, gk_440, gk_441, gk_442, \
                         gk_443, ik_907, ik_908, ik_909, ik_910, \
                         ik_911 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -3.0 * gk_439[k]
                   + f_0 * ik_907[k];

        t_656[k] = -3.0 * gk_440[k]
                   + f_0 * ik_908[k];

        t_657[k] = -3.0 * gk_441[k]
                   + f_0 * ik_909[k];

        t_658[k] = -3.0 * gk_442[k]
                   + f_0 * ik_910[k];

        t_659[k] = -3.0 * gk_443[k]
                   + f_0 * ik_911[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, gk_444, gk_445, gk_446, gk_447, \
                         gk_448, ik_912, ik_913, ik_914, ik_915, \
                         ik_916 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -3.0 * gk_444[k]
                   + f_0 * ik_912[k];

        t_661[k] = -3.0 * gk_445[k]
                   + f_0 * ik_913[k];

        t_662[k] = -3.0 * gk_446[k]
                   + f_0 * ik_914[k];

        t_663[k] = -3.0 * gk_447[k]
                   + f_0 * ik_915[k];

        t_664[k] = -3.0 * gk_448[k]
                   + f_0 * ik_916[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, gk_449, gk_450, gk_451, gk_452, \
                         gk_453, ik_917, ik_918, ik_919, ik_920, \
                         ik_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -3.0 * gk_449[k]
                   + f_0 * ik_917[k];

        t_666[k] = -3.0 * gk_450[k]
                   + f_0 * ik_918[k];

        t_667[k] = -3.0 * gk_451[k]
                   + f_0 * ik_919[k];

        t_668[k] = -3.0 * gk_452[k]
                   + f_0 * ik_920[k];

        t_669[k] = -3.0 * gk_453[k]
                   + f_0 * ik_921[k];
    }
}

static auto
compute_prim_geom_10_hk_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gk, const size_t ik,
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

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, gk_454, gk_455, gk_456, gk_457, \
                         gk_458, ik_922, ik_923, ik_924, ik_925, \
                         ik_926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -3.0 * gk_454[k]
                   + f_0 * ik_922[k];

        t_671[k] = -3.0 * gk_455[k]
                   + f_0 * ik_923[k];

        t_672[k] = -3.0 * gk_456[k]
                   + f_0 * ik_924[k];

        t_673[k] = -3.0 * gk_457[k]
                   + f_0 * ik_925[k];

        t_674[k] = -3.0 * gk_458[k]
                   + f_0 * ik_926[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, gk_459, gk_460, gk_461, gk_462, \
                         gk_463, ik_927, ik_928, ik_929, ik_930, \
                         ik_931 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -3.0 * gk_459[k]
                   + f_0 * ik_927[k];

        t_676[k] = -3.0 * gk_460[k]
                   + f_0 * ik_928[k];

        t_677[k] = -3.0 * gk_461[k]
                   + f_0 * ik_929[k];

        t_678[k] = -3.0 * gk_462[k]
                   + f_0 * ik_930[k];

        t_679[k] = -3.0 * gk_463[k]
                   + f_0 * ik_931[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, gk_464, gk_465, gk_466, gk_467, \
                         gk_468, ik_932, ik_933, ik_934, ik_935, \
                         ik_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -3.0 * gk_464[k]
                   + f_0 * ik_932[k];

        t_681[k] = -3.0 * gk_465[k]
                   + f_0 * ik_933[k];

        t_682[k] = -3.0 * gk_466[k]
                   + f_0 * ik_934[k];

        t_683[k] = -3.0 * gk_467[k]
                   + f_0 * ik_935[k];

        t_684[k] = -4.0 * gk_468[k]
                   + f_0 * ik_936[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, gk_469, gk_470, gk_471, gk_472, \
                         gk_473, ik_937, ik_938, ik_939, ik_940, \
                         ik_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -4.0 * gk_469[k]
                   + f_0 * ik_937[k];

        t_686[k] = -4.0 * gk_470[k]
                   + f_0 * ik_938[k];

        t_687[k] = -4.0 * gk_471[k]
                   + f_0 * ik_939[k];

        t_688[k] = -4.0 * gk_472[k]
                   + f_0 * ik_940[k];

        t_689[k] = -4.0 * gk_473[k]
                   + f_0 * ik_941[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, gk_474, gk_475, gk_476, gk_477, \
                         gk_478, ik_942, ik_943, ik_944, ik_945, \
                         ik_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -4.0 * gk_474[k]
                   + f_0 * ik_942[k];

        t_691[k] = -4.0 * gk_475[k]
                   + f_0 * ik_943[k];

        t_692[k] = -4.0 * gk_476[k]
                   + f_0 * ik_944[k];

        t_693[k] = -4.0 * gk_477[k]
                   + f_0 * ik_945[k];

        t_694[k] = -4.0 * gk_478[k]
                   + f_0 * ik_946[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, gk_479, gk_480, gk_481, gk_482, \
                         gk_483, ik_947, ik_948, ik_949, ik_950, \
                         ik_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -4.0 * gk_479[k]
                   + f_0 * ik_947[k];

        t_696[k] = -4.0 * gk_480[k]
                   + f_0 * ik_948[k];

        t_697[k] = -4.0 * gk_481[k]
                   + f_0 * ik_949[k];

        t_698[k] = -4.0 * gk_482[k]
                   + f_0 * ik_950[k];

        t_699[k] = -4.0 * gk_483[k]
                   + f_0 * ik_951[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, gk_484, gk_485, gk_486, gk_487, \
                         gk_488, ik_952, ik_953, ik_954, ik_955, \
                         ik_956 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -4.0 * gk_484[k]
                   + f_0 * ik_952[k];

        t_701[k] = -4.0 * gk_485[k]
                   + f_0 * ik_953[k];

        t_702[k] = -4.0 * gk_486[k]
                   + f_0 * ik_954[k];

        t_703[k] = -4.0 * gk_487[k]
                   + f_0 * ik_955[k];

        t_704[k] = -4.0 * gk_488[k]
                   + f_0 * ik_956[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, gk_489, gk_490, gk_491, gk_492, \
                         gk_493, ik_957, ik_958, ik_959, ik_960, \
                         ik_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -4.0 * gk_489[k]
                   + f_0 * ik_957[k];

        t_706[k] = -4.0 * gk_490[k]
                   + f_0 * ik_958[k];

        t_707[k] = -4.0 * gk_491[k]
                   + f_0 * ik_959[k];

        t_708[k] = -4.0 * gk_492[k]
                   + f_0 * ik_960[k];

        t_709[k] = -4.0 * gk_493[k]
                   + f_0 * ik_961[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, gk_494, gk_495, gk_496, gk_497, \
                         gk_498, ik_962, ik_963, ik_964, ik_965, \
                         ik_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -4.0 * gk_494[k]
                   + f_0 * ik_962[k];

        t_711[k] = -4.0 * gk_495[k]
                   + f_0 * ik_963[k];

        t_712[k] = -4.0 * gk_496[k]
                   + f_0 * ik_964[k];

        t_713[k] = -4.0 * gk_497[k]
                   + f_0 * ik_965[k];

        t_714[k] = -4.0 * gk_498[k]
                   + f_0 * ik_966[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, gk_499, gk_500, gk_501, gk_502, \
                         gk_503, ik_967, ik_968, ik_969, ik_970, \
                         ik_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -4.0 * gk_499[k]
                   + f_0 * ik_967[k];

        t_716[k] = -4.0 * gk_500[k]
                   + f_0 * ik_968[k];

        t_717[k] = -4.0 * gk_501[k]
                   + f_0 * ik_969[k];

        t_718[k] = -4.0 * gk_502[k]
                   + f_0 * ik_970[k];

        t_719[k] = -4.0 * gk_503[k]
                   + f_0 * ik_971[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, gk_504, gk_505, gk_506, gk_507, \
                         gk_508, ik_972, ik_973, ik_974, ik_975, \
                         ik_976 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -5.0 * gk_504[k]
                   + f_0 * ik_972[k];

        t_721[k] = -5.0 * gk_505[k]
                   + f_0 * ik_973[k];

        t_722[k] = -5.0 * gk_506[k]
                   + f_0 * ik_974[k];

        t_723[k] = -5.0 * gk_507[k]
                   + f_0 * ik_975[k];

        t_724[k] = -5.0 * gk_508[k]
                   + f_0 * ik_976[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, gk_509, gk_510, gk_511, gk_512, \
                         gk_513, ik_977, ik_978, ik_979, ik_980, \
                         ik_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -5.0 * gk_509[k]
                   + f_0 * ik_977[k];

        t_726[k] = -5.0 * gk_510[k]
                   + f_0 * ik_978[k];

        t_727[k] = -5.0 * gk_511[k]
                   + f_0 * ik_979[k];

        t_728[k] = -5.0 * gk_512[k]
                   + f_0 * ik_980[k];

        t_729[k] = -5.0 * gk_513[k]
                   + f_0 * ik_981[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, gk_514, gk_515, gk_516, gk_517, \
                         gk_518, ik_982, ik_983, ik_984, ik_985, \
                         ik_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -5.0 * gk_514[k]
                   + f_0 * ik_982[k];

        t_731[k] = -5.0 * gk_515[k]
                   + f_0 * ik_983[k];

        t_732[k] = -5.0 * gk_516[k]
                   + f_0 * ik_984[k];

        t_733[k] = -5.0 * gk_517[k]
                   + f_0 * ik_985[k];

        t_734[k] = -5.0 * gk_518[k]
                   + f_0 * ik_986[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, gk_519, gk_520, gk_521, gk_522, \
                         gk_523, ik_987, ik_988, ik_989, ik_990, \
                         ik_991 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -5.0 * gk_519[k]
                   + f_0 * ik_987[k];

        t_736[k] = -5.0 * gk_520[k]
                   + f_0 * ik_988[k];

        t_737[k] = -5.0 * gk_521[k]
                   + f_0 * ik_989[k];

        t_738[k] = -5.0 * gk_522[k]
                   + f_0 * ik_990[k];

        t_739[k] = -5.0 * gk_523[k]
                   + f_0 * ik_991[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, gk_524, gk_525, gk_526, gk_527, \
                         gk_528, ik_992, ik_993, ik_994, ik_995, \
                         ik_996 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -5.0 * gk_524[k]
                   + f_0 * ik_992[k];

        t_741[k] = -5.0 * gk_525[k]
                   + f_0 * ik_993[k];

        t_742[k] = -5.0 * gk_526[k]
                   + f_0 * ik_994[k];

        t_743[k] = -5.0 * gk_527[k]
                   + f_0 * ik_995[k];

        t_744[k] = -5.0 * gk_528[k]
                   + f_0 * ik_996[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, gk_529, gk_530, gk_531, gk_532, \
                         gk_533, ik_997, ik_998, ik_999, ik_1000, \
                         ik_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -5.0 * gk_529[k]
                   + f_0 * ik_997[k];

        t_746[k] = -5.0 * gk_530[k]
                   + f_0 * ik_998[k];

        t_747[k] = -5.0 * gk_531[k]
                   + f_0 * ik_999[k];

        t_748[k] = -5.0 * gk_532[k]
                   + f_0 * ik_1000[k];

        t_749[k] = -5.0 * gk_533[k]
                   + f_0 * ik_1001[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, gk_534, gk_535, gk_536, gk_537, \
                         gk_538, ik_1002, ik_1003, ik_1004, ik_1005, \
                         ik_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -5.0 * gk_534[k]
                   + f_0 * ik_1002[k];

        t_751[k] = -5.0 * gk_535[k]
                   + f_0 * ik_1003[k];

        t_752[k] = -5.0 * gk_536[k]
                   + f_0 * ik_1004[k];

        t_753[k] = -5.0 * gk_537[k]
                   + f_0 * ik_1005[k];

        t_754[k] = -5.0 * gk_538[k]
                   + f_0 * ik_1006[k];
    }

#pragma omp simd aligned(t_755, gk_539, ik_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -5.0 * gk_539[k]
                   + f_0 * ik_1007[k];
    }
}

auto
compute_prim_geom_10_hk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t gk, const size_t ik,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hk_electron_repulsion_2_piece0(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_2_piece1(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_2_piece2(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_2_piece3(buffer, target, gk, ik, ncols, alpha);

    compute_prim_geom_10_hk_electron_repulsion_2_piece4(buffer, target, gk, ik, ncols, alpha);
}

}  // namespace simdt2ceri
