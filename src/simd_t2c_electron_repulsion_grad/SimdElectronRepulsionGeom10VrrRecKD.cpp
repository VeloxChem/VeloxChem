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


#include "SimdElectronRepulsionGeom10VrrRecKD.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_kd_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t id, const size_t ld,
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

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_73 = buffer.data(id + 73);
    const auto *id_74 = buffer.data(id + 74);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_79 = buffer.data(id + 79);
    const auto *id_80 = buffer.data(id + 80);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_85 = buffer.data(id + 85);
    const auto *id_86 = buffer.data(id + 86);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_88 = buffer.data(id + 88);
    const auto *id_89 = buffer.data(id + 89);
    const auto *id_90 = buffer.data(id + 90);
    const auto *id_91 = buffer.data(id + 91);
    const auto *id_92 = buffer.data(id + 92);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_94 = buffer.data(id + 94);
    const auto *id_95 = buffer.data(id + 95);
    const auto *id_96 = buffer.data(id + 96);
    const auto *id_97 = buffer.data(id + 97);
    const auto *id_98 = buffer.data(id + 98);
    const auto *id_99 = buffer.data(id + 99);
    const auto *id_100 = buffer.data(id + 100);
    const auto *id_101 = buffer.data(id + 101);
    const auto *id_102 = buffer.data(id + 102);
    const auto *id_103 = buffer.data(id + 103);
    const auto *id_104 = buffer.data(id + 104);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_106 = buffer.data(id + 106);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_108 = buffer.data(id + 108);
    const auto *id_109 = buffer.data(id + 109);
    const auto *id_110 = buffer.data(id + 110);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_115 = buffer.data(id + 115);
    const auto *id_116 = buffer.data(id + 116);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_119 = buffer.data(id + 119);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_121 = buffer.data(id + 121);
    const auto *id_122 = buffer.data(id + 122);
    const auto *id_123 = buffer.data(id + 123);
    const auto *id_124 = buffer.data(id + 124);
    const auto *id_125 = buffer.data(id + 125);
    const auto *id_126 = buffer.data(id + 126);
    const auto *id_127 = buffer.data(id + 127);
    const auto *id_128 = buffer.data(id + 128);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_130 = buffer.data(id + 130);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_132 = buffer.data(id + 132);
    const auto *id_133 = buffer.data(id + 133);
    const auto *id_134 = buffer.data(id + 134);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_138 = buffer.data(id + 138);
    const auto *id_139 = buffer.data(id + 139);
    const auto *id_140 = buffer.data(id + 140);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_144 = buffer.data(id + 144);
    const auto *id_145 = buffer.data(id + 145);
    const auto *id_146 = buffer.data(id + 146);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);
    const auto *ld_54 = buffer.data(ld + 54);
    const auto *ld_55 = buffer.data(ld + 55);
    const auto *ld_56 = buffer.data(ld + 56);
    const auto *ld_57 = buffer.data(ld + 57);
    const auto *ld_58 = buffer.data(ld + 58);
    const auto *ld_59 = buffer.data(ld + 59);
    const auto *ld_60 = buffer.data(ld + 60);
    const auto *ld_61 = buffer.data(ld + 61);
    const auto *ld_62 = buffer.data(ld + 62);
    const auto *ld_63 = buffer.data(ld + 63);
    const auto *ld_64 = buffer.data(ld + 64);
    const auto *ld_65 = buffer.data(ld + 65);
    const auto *ld_66 = buffer.data(ld + 66);
    const auto *ld_67 = buffer.data(ld + 67);
    const auto *ld_68 = buffer.data(ld + 68);
    const auto *ld_69 = buffer.data(ld + 69);
    const auto *ld_70 = buffer.data(ld + 70);
    const auto *ld_71 = buffer.data(ld + 71);
    const auto *ld_72 = buffer.data(ld + 72);
    const auto *ld_73 = buffer.data(ld + 73);
    const auto *ld_74 = buffer.data(ld + 74);
    const auto *ld_75 = buffer.data(ld + 75);
    const auto *ld_76 = buffer.data(ld + 76);
    const auto *ld_77 = buffer.data(ld + 77);
    const auto *ld_78 = buffer.data(ld + 78);
    const auto *ld_79 = buffer.data(ld + 79);
    const auto *ld_80 = buffer.data(ld + 80);
    const auto *ld_81 = buffer.data(ld + 81);
    const auto *ld_82 = buffer.data(ld + 82);
    const auto *ld_83 = buffer.data(ld + 83);
    const auto *ld_84 = buffer.data(ld + 84);
    const auto *ld_85 = buffer.data(ld + 85);
    const auto *ld_86 = buffer.data(ld + 86);
    const auto *ld_87 = buffer.data(ld + 87);
    const auto *ld_88 = buffer.data(ld + 88);
    const auto *ld_89 = buffer.data(ld + 89);
    const auto *ld_90 = buffer.data(ld + 90);
    const auto *ld_91 = buffer.data(ld + 91);
    const auto *ld_92 = buffer.data(ld + 92);
    const auto *ld_93 = buffer.data(ld + 93);
    const auto *ld_94 = buffer.data(ld + 94);
    const auto *ld_95 = buffer.data(ld + 95);
    const auto *ld_96 = buffer.data(ld + 96);
    const auto *ld_97 = buffer.data(ld + 97);
    const auto *ld_98 = buffer.data(ld + 98);
    const auto *ld_99 = buffer.data(ld + 99);
    const auto *ld_100 = buffer.data(ld + 100);
    const auto *ld_101 = buffer.data(ld + 101);
    const auto *ld_102 = buffer.data(ld + 102);
    const auto *ld_103 = buffer.data(ld + 103);
    const auto *ld_104 = buffer.data(ld + 104);
    const auto *ld_105 = buffer.data(ld + 105);
    const auto *ld_106 = buffer.data(ld + 106);
    const auto *ld_107 = buffer.data(ld + 107);
    const auto *ld_108 = buffer.data(ld + 108);
    const auto *ld_109 = buffer.data(ld + 109);
    const auto *ld_110 = buffer.data(ld + 110);
    const auto *ld_111 = buffer.data(ld + 111);
    const auto *ld_112 = buffer.data(ld + 112);
    const auto *ld_113 = buffer.data(ld + 113);
    const auto *ld_114 = buffer.data(ld + 114);
    const auto *ld_115 = buffer.data(ld + 115);
    const auto *ld_116 = buffer.data(ld + 116);
    const auto *ld_117 = buffer.data(ld + 117);
    const auto *ld_118 = buffer.data(ld + 118);
    const auto *ld_119 = buffer.data(ld + 119);
    const auto *ld_120 = buffer.data(ld + 120);
    const auto *ld_121 = buffer.data(ld + 121);
    const auto *ld_122 = buffer.data(ld + 122);
    const auto *ld_123 = buffer.data(ld + 123);
    const auto *ld_124 = buffer.data(ld + 124);
    const auto *ld_125 = buffer.data(ld + 125);
    const auto *ld_126 = buffer.data(ld + 126);
    const auto *ld_127 = buffer.data(ld + 127);
    const auto *ld_128 = buffer.data(ld + 128);
    const auto *ld_129 = buffer.data(ld + 129);
    const auto *ld_130 = buffer.data(ld + 130);
    const auto *ld_131 = buffer.data(ld + 131);
    const auto *ld_132 = buffer.data(ld + 132);
    const auto *ld_133 = buffer.data(ld + 133);
    const auto *ld_134 = buffer.data(ld + 134);
    const auto *ld_135 = buffer.data(ld + 135);
    const auto *ld_136 = buffer.data(ld + 136);
    const auto *ld_137 = buffer.data(ld + 137);
    const auto *ld_138 = buffer.data(ld + 138);
    const auto *ld_139 = buffer.data(ld + 139);
    const auto *ld_140 = buffer.data(ld + 140);
    const auto *ld_141 = buffer.data(ld + 141);
    const auto *ld_142 = buffer.data(ld + 142);
    const auto *ld_143 = buffer.data(ld + 143);
    const auto *ld_144 = buffer.data(ld + 144);
    const auto *ld_145 = buffer.data(ld + 145);
    const auto *ld_146 = buffer.data(ld + 146);
    const auto *ld_147 = buffer.data(ld + 147);
    const auto *ld_148 = buffer.data(ld + 148);
    const auto *ld_149 = buffer.data(ld + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, id_0, id_1, id_2, id_3, id_4, ld_0, ld_1, \
                         ld_2, ld_3, ld_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -7.0 * id_0[k]
                 + f_0 * ld_0[k];

        t_1[k] = -7.0 * id_1[k]
                 + f_0 * ld_1[k];

        t_2[k] = -7.0 * id_2[k]
                 + f_0 * ld_2[k];

        t_3[k] = -7.0 * id_3[k]
                 + f_0 * ld_3[k];

        t_4[k] = -7.0 * id_4[k]
                 + f_0 * ld_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, id_5, id_6, id_7, id_8, id_9, ld_5, ld_6, \
                         ld_7, ld_8, ld_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -7.0 * id_5[k]
                 + f_0 * ld_5[k];

        t_6[k] = -6.0 * id_6[k]
                 + f_0 * ld_6[k];

        t_7[k] = -6.0 * id_7[k]
                 + f_0 * ld_7[k];

        t_8[k] = -6.0 * id_8[k]
                 + f_0 * ld_8[k];

        t_9[k] = -6.0 * id_9[k]
                 + f_0 * ld_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, id_10, id_11, id_12, id_13, id_14, \
                         ld_10, ld_11, ld_12, ld_13, ld_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -6.0 * id_10[k]
                  + f_0 * ld_10[k];

        t_11[k] = -6.0 * id_11[k]
                  + f_0 * ld_11[k];

        t_12[k] = -6.0 * id_12[k]
                  + f_0 * ld_12[k];

        t_13[k] = -6.0 * id_13[k]
                  + f_0 * ld_13[k];

        t_14[k] = -6.0 * id_14[k]
                  + f_0 * ld_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, id_15, id_16, id_17, id_18, id_19, \
                         ld_15, ld_16, ld_17, ld_18, ld_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -6.0 * id_15[k]
                  + f_0 * ld_15[k];

        t_16[k] = -6.0 * id_16[k]
                  + f_0 * ld_16[k];

        t_17[k] = -6.0 * id_17[k]
                  + f_0 * ld_17[k];

        t_18[k] = -5.0 * id_18[k]
                  + f_0 * ld_18[k];

        t_19[k] = -5.0 * id_19[k]
                  + f_0 * ld_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, id_20, id_21, id_22, id_23, id_24, \
                         ld_20, ld_21, ld_22, ld_23, ld_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -5.0 * id_20[k]
                  + f_0 * ld_20[k];

        t_21[k] = -5.0 * id_21[k]
                  + f_0 * ld_21[k];

        t_22[k] = -5.0 * id_22[k]
                  + f_0 * ld_22[k];

        t_23[k] = -5.0 * id_23[k]
                  + f_0 * ld_23[k];

        t_24[k] = -5.0 * id_24[k]
                  + f_0 * ld_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, id_25, id_26, id_27, id_28, id_29, \
                         ld_25, ld_26, ld_27, ld_28, ld_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -5.0 * id_25[k]
                  + f_0 * ld_25[k];

        t_26[k] = -5.0 * id_26[k]
                  + f_0 * ld_26[k];

        t_27[k] = -5.0 * id_27[k]
                  + f_0 * ld_27[k];

        t_28[k] = -5.0 * id_28[k]
                  + f_0 * ld_28[k];

        t_29[k] = -5.0 * id_29[k]
                  + f_0 * ld_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, id_30, id_31, id_32, id_33, id_34, \
                         ld_30, ld_31, ld_32, ld_33, ld_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -5.0 * id_30[k]
                  + f_0 * ld_30[k];

        t_31[k] = -5.0 * id_31[k]
                  + f_0 * ld_31[k];

        t_32[k] = -5.0 * id_32[k]
                  + f_0 * ld_32[k];

        t_33[k] = -5.0 * id_33[k]
                  + f_0 * ld_33[k];

        t_34[k] = -5.0 * id_34[k]
                  + f_0 * ld_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, id_35, id_36, id_37, id_38, id_39, \
                         ld_35, ld_36, ld_37, ld_38, ld_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -5.0 * id_35[k]
                  + f_0 * ld_35[k];

        t_36[k] = -4.0 * id_36[k]
                  + f_0 * ld_36[k];

        t_37[k] = -4.0 * id_37[k]
                  + f_0 * ld_37[k];

        t_38[k] = -4.0 * id_38[k]
                  + f_0 * ld_38[k];

        t_39[k] = -4.0 * id_39[k]
                  + f_0 * ld_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, id_40, id_41, id_42, id_43, id_44, \
                         ld_40, ld_41, ld_42, ld_43, ld_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -4.0 * id_40[k]
                  + f_0 * ld_40[k];

        t_41[k] = -4.0 * id_41[k]
                  + f_0 * ld_41[k];

        t_42[k] = -4.0 * id_42[k]
                  + f_0 * ld_42[k];

        t_43[k] = -4.0 * id_43[k]
                  + f_0 * ld_43[k];

        t_44[k] = -4.0 * id_44[k]
                  + f_0 * ld_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, id_45, id_46, id_47, id_48, id_49, \
                         ld_45, ld_46, ld_47, ld_48, ld_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -4.0 * id_45[k]
                  + f_0 * ld_45[k];

        t_46[k] = -4.0 * id_46[k]
                  + f_0 * ld_46[k];

        t_47[k] = -4.0 * id_47[k]
                  + f_0 * ld_47[k];

        t_48[k] = -4.0 * id_48[k]
                  + f_0 * ld_48[k];

        t_49[k] = -4.0 * id_49[k]
                  + f_0 * ld_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, id_50, id_51, id_52, id_53, id_54, \
                         ld_50, ld_51, ld_52, ld_53, ld_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -4.0 * id_50[k]
                  + f_0 * ld_50[k];

        t_51[k] = -4.0 * id_51[k]
                  + f_0 * ld_51[k];

        t_52[k] = -4.0 * id_52[k]
                  + f_0 * ld_52[k];

        t_53[k] = -4.0 * id_53[k]
                  + f_0 * ld_53[k];

        t_54[k] = -4.0 * id_54[k]
                  + f_0 * ld_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, id_55, id_56, id_57, id_58, id_59, \
                         ld_55, ld_56, ld_57, ld_58, ld_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -4.0 * id_55[k]
                  + f_0 * ld_55[k];

        t_56[k] = -4.0 * id_56[k]
                  + f_0 * ld_56[k];

        t_57[k] = -4.0 * id_57[k]
                  + f_0 * ld_57[k];

        t_58[k] = -4.0 * id_58[k]
                  + f_0 * ld_58[k];

        t_59[k] = -4.0 * id_59[k]
                  + f_0 * ld_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, id_60, id_61, id_62, id_63, id_64, \
                         ld_60, ld_61, ld_62, ld_63, ld_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * id_60[k]
                  + f_0 * ld_60[k];

        t_61[k] = -3.0 * id_61[k]
                  + f_0 * ld_61[k];

        t_62[k] = -3.0 * id_62[k]
                  + f_0 * ld_62[k];

        t_63[k] = -3.0 * id_63[k]
                  + f_0 * ld_63[k];

        t_64[k] = -3.0 * id_64[k]
                  + f_0 * ld_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, id_65, id_66, id_67, id_68, id_69, \
                         ld_65, ld_66, ld_67, ld_68, ld_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * id_65[k]
                  + f_0 * ld_65[k];

        t_66[k] = -3.0 * id_66[k]
                  + f_0 * ld_66[k];

        t_67[k] = -3.0 * id_67[k]
                  + f_0 * ld_67[k];

        t_68[k] = -3.0 * id_68[k]
                  + f_0 * ld_68[k];

        t_69[k] = -3.0 * id_69[k]
                  + f_0 * ld_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, id_70, id_71, id_72, id_73, id_74, \
                         ld_70, ld_71, ld_72, ld_73, ld_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -3.0 * id_70[k]
                  + f_0 * ld_70[k];

        t_71[k] = -3.0 * id_71[k]
                  + f_0 * ld_71[k];

        t_72[k] = -3.0 * id_72[k]
                  + f_0 * ld_72[k];

        t_73[k] = -3.0 * id_73[k]
                  + f_0 * ld_73[k];

        t_74[k] = -3.0 * id_74[k]
                  + f_0 * ld_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, id_75, id_76, id_77, id_78, id_79, \
                         ld_75, ld_76, ld_77, ld_78, ld_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -3.0 * id_75[k]
                  + f_0 * ld_75[k];

        t_76[k] = -3.0 * id_76[k]
                  + f_0 * ld_76[k];

        t_77[k] = -3.0 * id_77[k]
                  + f_0 * ld_77[k];

        t_78[k] = -3.0 * id_78[k]
                  + f_0 * ld_78[k];

        t_79[k] = -3.0 * id_79[k]
                  + f_0 * ld_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, id_80, id_81, id_82, id_83, id_84, \
                         ld_80, ld_81, ld_82, ld_83, ld_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * id_80[k]
                  + f_0 * ld_80[k];

        t_81[k] = -3.0 * id_81[k]
                  + f_0 * ld_81[k];

        t_82[k] = -3.0 * id_82[k]
                  + f_0 * ld_82[k];

        t_83[k] = -3.0 * id_83[k]
                  + f_0 * ld_83[k];

        t_84[k] = -3.0 * id_84[k]
                  + f_0 * ld_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, id_85, id_86, id_87, id_88, id_89, \
                         ld_85, ld_86, ld_87, ld_88, ld_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -3.0 * id_85[k]
                  + f_0 * ld_85[k];

        t_86[k] = -3.0 * id_86[k]
                  + f_0 * ld_86[k];

        t_87[k] = -3.0 * id_87[k]
                  + f_0 * ld_87[k];

        t_88[k] = -3.0 * id_88[k]
                  + f_0 * ld_88[k];

        t_89[k] = -3.0 * id_89[k]
                  + f_0 * ld_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, id_90, id_91, id_92, id_93, id_94, \
                         ld_90, ld_91, ld_92, ld_93, ld_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * id_90[k]
                  + f_0 * ld_90[k];

        t_91[k] = -2.0 * id_91[k]
                  + f_0 * ld_91[k];

        t_92[k] = -2.0 * id_92[k]
                  + f_0 * ld_92[k];

        t_93[k] = -2.0 * id_93[k]
                  + f_0 * ld_93[k];

        t_94[k] = -2.0 * id_94[k]
                  + f_0 * ld_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, id_95, id_96, id_97, id_98, id_99, \
                         ld_95, ld_96, ld_97, ld_98, ld_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * id_95[k]
                  + f_0 * ld_95[k];

        t_96[k] = -2.0 * id_96[k]
                  + f_0 * ld_96[k];

        t_97[k] = -2.0 * id_97[k]
                  + f_0 * ld_97[k];

        t_98[k] = -2.0 * id_98[k]
                  + f_0 * ld_98[k];

        t_99[k] = -2.0 * id_99[k]
                  + f_0 * ld_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, id_100, id_101, id_102, id_103, \
                         id_104, ld_100, ld_101, ld_102, ld_103, \
                         ld_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * id_100[k]
                   + f_0 * ld_100[k];

        t_101[k] = -2.0 * id_101[k]
                   + f_0 * ld_101[k];

        t_102[k] = -2.0 * id_102[k]
                   + f_0 * ld_102[k];

        t_103[k] = -2.0 * id_103[k]
                   + f_0 * ld_103[k];

        t_104[k] = -2.0 * id_104[k]
                   + f_0 * ld_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, id_105, id_106, id_107, id_108, \
                         id_109, ld_105, ld_106, ld_107, ld_108, \
                         ld_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * id_105[k]
                   + f_0 * ld_105[k];

        t_106[k] = -2.0 * id_106[k]
                   + f_0 * ld_106[k];

        t_107[k] = -2.0 * id_107[k]
                   + f_0 * ld_107[k];

        t_108[k] = -2.0 * id_108[k]
                   + f_0 * ld_108[k];

        t_109[k] = -2.0 * id_109[k]
                   + f_0 * ld_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, id_110, id_111, id_112, id_113, \
                         id_114, ld_110, ld_111, ld_112, ld_113, \
                         ld_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * id_110[k]
                   + f_0 * ld_110[k];

        t_111[k] = -2.0 * id_111[k]
                   + f_0 * ld_111[k];

        t_112[k] = -2.0 * id_112[k]
                   + f_0 * ld_112[k];

        t_113[k] = -2.0 * id_113[k]
                   + f_0 * ld_113[k];

        t_114[k] = -2.0 * id_114[k]
                   + f_0 * ld_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, id_115, id_116, id_117, id_118, \
                         id_119, ld_115, ld_116, ld_117, ld_118, \
                         ld_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * id_115[k]
                   + f_0 * ld_115[k];

        t_116[k] = -2.0 * id_116[k]
                   + f_0 * ld_116[k];

        t_117[k] = -2.0 * id_117[k]
                   + f_0 * ld_117[k];

        t_118[k] = -2.0 * id_118[k]
                   + f_0 * ld_118[k];

        t_119[k] = -2.0 * id_119[k]
                   + f_0 * ld_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, id_120, id_121, id_122, id_123, \
                         id_124, ld_120, ld_121, ld_122, ld_123, \
                         ld_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -2.0 * id_120[k]
                   + f_0 * ld_120[k];

        t_121[k] = -2.0 * id_121[k]
                   + f_0 * ld_121[k];

        t_122[k] = -2.0 * id_122[k]
                   + f_0 * ld_122[k];

        t_123[k] = -2.0 * id_123[k]
                   + f_0 * ld_123[k];

        t_124[k] = -2.0 * id_124[k]
                   + f_0 * ld_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, id_125, id_126, id_127, id_128, \
                         id_129, ld_125, ld_126, ld_127, ld_128, \
                         ld_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -2.0 * id_125[k]
                   + f_0 * ld_125[k];

        t_126[k] = -id_126[k]
                   + f_0 * ld_126[k];

        t_127[k] = -id_127[k]
                   + f_0 * ld_127[k];

        t_128[k] = -id_128[k]
                   + f_0 * ld_128[k];

        t_129[k] = -id_129[k]
                   + f_0 * ld_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, id_130, id_131, id_132, id_133, \
                         id_134, ld_130, ld_131, ld_132, ld_133, \
                         ld_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -id_130[k]
                   + f_0 * ld_130[k];

        t_131[k] = -id_131[k]
                   + f_0 * ld_131[k];

        t_132[k] = -id_132[k]
                   + f_0 * ld_132[k];

        t_133[k] = -id_133[k]
                   + f_0 * ld_133[k];

        t_134[k] = -id_134[k]
                   + f_0 * ld_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, id_135, id_136, id_137, id_138, \
                         id_139, ld_135, ld_136, ld_137, ld_138, \
                         ld_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -id_135[k]
                   + f_0 * ld_135[k];

        t_136[k] = -id_136[k]
                   + f_0 * ld_136[k];

        t_137[k] = -id_137[k]
                   + f_0 * ld_137[k];

        t_138[k] = -id_138[k]
                   + f_0 * ld_138[k];

        t_139[k] = -id_139[k]
                   + f_0 * ld_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, id_140, id_141, id_142, id_143, \
                         id_144, ld_140, ld_141, ld_142, ld_143, \
                         ld_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -id_140[k]
                   + f_0 * ld_140[k];

        t_141[k] = -id_141[k]
                   + f_0 * ld_141[k];

        t_142[k] = -id_142[k]
                   + f_0 * ld_142[k];

        t_143[k] = -id_143[k]
                   + f_0 * ld_143[k];

        t_144[k] = -id_144[k]
                   + f_0 * ld_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, id_145, id_146, id_147, id_148, \
                         id_149, ld_145, ld_146, ld_147, ld_148, \
                         ld_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -id_145[k]
                   + f_0 * ld_145[k];

        t_146[k] = -id_146[k]
                   + f_0 * ld_146[k];

        t_147[k] = -id_147[k]
                   + f_0 * ld_147[k];

        t_148[k] = -id_148[k]
                   + f_0 * ld_148[k];

        t_149[k] = -id_149[k]
                   + f_0 * ld_149[k];
    }
}

static auto
compute_prim_geom_10_kd_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t id, const size_t ld,
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

    const auto *id_150 = buffer.data(id + 150);
    const auto *id_151 = buffer.data(id + 151);
    const auto *id_152 = buffer.data(id + 152);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_156 = buffer.data(id + 156);
    const auto *id_157 = buffer.data(id + 157);
    const auto *id_158 = buffer.data(id + 158);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_163 = buffer.data(id + 163);
    const auto *id_164 = buffer.data(id + 164);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_166 = buffer.data(id + 166);
    const auto *id_167 = buffer.data(id + 167);

    const auto *ld_150 = buffer.data(ld + 150);
    const auto *ld_151 = buffer.data(ld + 151);
    const auto *ld_152 = buffer.data(ld + 152);
    const auto *ld_153 = buffer.data(ld + 153);
    const auto *ld_154 = buffer.data(ld + 154);
    const auto *ld_155 = buffer.data(ld + 155);
    const auto *ld_156 = buffer.data(ld + 156);
    const auto *ld_157 = buffer.data(ld + 157);
    const auto *ld_158 = buffer.data(ld + 158);
    const auto *ld_159 = buffer.data(ld + 159);
    const auto *ld_160 = buffer.data(ld + 160);
    const auto *ld_161 = buffer.data(ld + 161);
    const auto *ld_162 = buffer.data(ld + 162);
    const auto *ld_163 = buffer.data(ld + 163);
    const auto *ld_164 = buffer.data(ld + 164);
    const auto *ld_165 = buffer.data(ld + 165);
    const auto *ld_166 = buffer.data(ld + 166);
    const auto *ld_167 = buffer.data(ld + 167);
    const auto *ld_168 = buffer.data(ld + 168);
    const auto *ld_169 = buffer.data(ld + 169);
    const auto *ld_170 = buffer.data(ld + 170);
    const auto *ld_171 = buffer.data(ld + 171);
    const auto *ld_172 = buffer.data(ld + 172);
    const auto *ld_173 = buffer.data(ld + 173);
    const auto *ld_174 = buffer.data(ld + 174);
    const auto *ld_175 = buffer.data(ld + 175);
    const auto *ld_176 = buffer.data(ld + 176);
    const auto *ld_177 = buffer.data(ld + 177);
    const auto *ld_178 = buffer.data(ld + 178);
    const auto *ld_179 = buffer.data(ld + 179);
    const auto *ld_180 = buffer.data(ld + 180);
    const auto *ld_181 = buffer.data(ld + 181);
    const auto *ld_182 = buffer.data(ld + 182);
    const auto *ld_183 = buffer.data(ld + 183);
    const auto *ld_184 = buffer.data(ld + 184);
    const auto *ld_185 = buffer.data(ld + 185);
    const auto *ld_186 = buffer.data(ld + 186);
    const auto *ld_187 = buffer.data(ld + 187);
    const auto *ld_188 = buffer.data(ld + 188);
    const auto *ld_189 = buffer.data(ld + 189);
    const auto *ld_190 = buffer.data(ld + 190);
    const auto *ld_191 = buffer.data(ld + 191);
    const auto *ld_192 = buffer.data(ld + 192);
    const auto *ld_193 = buffer.data(ld + 193);
    const auto *ld_194 = buffer.data(ld + 194);
    const auto *ld_195 = buffer.data(ld + 195);
    const auto *ld_196 = buffer.data(ld + 196);
    const auto *ld_197 = buffer.data(ld + 197);
    const auto *ld_198 = buffer.data(ld + 198);
    const auto *ld_199 = buffer.data(ld + 199);
    const auto *ld_200 = buffer.data(ld + 200);
    const auto *ld_201 = buffer.data(ld + 201);
    const auto *ld_202 = buffer.data(ld + 202);
    const auto *ld_203 = buffer.data(ld + 203);
    const auto *ld_204 = buffer.data(ld + 204);
    const auto *ld_205 = buffer.data(ld + 205);
    const auto *ld_206 = buffer.data(ld + 206);
    const auto *ld_207 = buffer.data(ld + 207);
    const auto *ld_208 = buffer.data(ld + 208);
    const auto *ld_209 = buffer.data(ld + 209);
    const auto *ld_210 = buffer.data(ld + 210);
    const auto *ld_211 = buffer.data(ld + 211);
    const auto *ld_212 = buffer.data(ld + 212);
    const auto *ld_213 = buffer.data(ld + 213);
    const auto *ld_214 = buffer.data(ld + 214);
    const auto *ld_215 = buffer.data(ld + 215);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, id_150, id_151, id_152, id_153, \
                         id_154, ld_150, ld_151, ld_152, ld_153, \
                         ld_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -id_150[k]
                   + f_0 * ld_150[k];

        t_151[k] = -id_151[k]
                   + f_0 * ld_151[k];

        t_152[k] = -id_152[k]
                   + f_0 * ld_152[k];

        t_153[k] = -id_153[k]
                   + f_0 * ld_153[k];

        t_154[k] = -id_154[k]
                   + f_0 * ld_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, id_155, id_156, id_157, id_158, \
                         id_159, ld_155, ld_156, ld_157, ld_158, \
                         ld_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -id_155[k]
                   + f_0 * ld_155[k];

        t_156[k] = -id_156[k]
                   + f_0 * ld_156[k];

        t_157[k] = -id_157[k]
                   + f_0 * ld_157[k];

        t_158[k] = -id_158[k]
                   + f_0 * ld_158[k];

        t_159[k] = -id_159[k]
                   + f_0 * ld_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, id_160, id_161, id_162, id_163, \
                         id_164, ld_160, ld_161, ld_162, ld_163, \
                         ld_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -id_160[k]
                   + f_0 * ld_160[k];

        t_161[k] = -id_161[k]
                   + f_0 * ld_161[k];

        t_162[k] = -id_162[k]
                   + f_0 * ld_162[k];

        t_163[k] = -id_163[k]
                   + f_0 * ld_163[k];

        t_164[k] = -id_164[k]
                   + f_0 * ld_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, t_170, id_165, id_166, id_167, \
                         ld_165, ld_166, ld_167, ld_168, ld_169, \
                         ld_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -id_165[k]
                   + f_0 * ld_165[k];

        t_166[k] = -id_166[k]
                   + f_0 * ld_166[k];

        t_167[k] = -id_167[k]
                   + f_0 * ld_167[k];

        t_168[k] = f_0 * ld_168[k];

        t_169[k] = f_0 * ld_169[k];

        t_170[k] = f_0 * ld_170[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, t_175, t_176, t_177, t_178, ld_171, \
                         ld_172, ld_173, ld_174, ld_175, ld_176, ld_177, \
                         ld_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_0 * ld_171[k];

        t_172[k] = f_0 * ld_172[k];

        t_173[k] = f_0 * ld_173[k];

        t_174[k] = f_0 * ld_174[k];

        t_175[k] = f_0 * ld_175[k];

        t_176[k] = f_0 * ld_176[k];

        t_177[k] = f_0 * ld_177[k];

        t_178[k] = f_0 * ld_178[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, t_184, t_185, t_186, ld_179, \
                         ld_180, ld_181, ld_182, ld_183, ld_184, ld_185, \
                         ld_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_0 * ld_179[k];

        t_180[k] = f_0 * ld_180[k];

        t_181[k] = f_0 * ld_181[k];

        t_182[k] = f_0 * ld_182[k];

        t_183[k] = f_0 * ld_183[k];

        t_184[k] = f_0 * ld_184[k];

        t_185[k] = f_0 * ld_185[k];

        t_186[k] = f_0 * ld_186[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, t_192, t_193, t_194, ld_187, \
                         ld_188, ld_189, ld_190, ld_191, ld_192, ld_193, \
                         ld_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_0 * ld_187[k];

        t_188[k] = f_0 * ld_188[k];

        t_189[k] = f_0 * ld_189[k];

        t_190[k] = f_0 * ld_190[k];

        t_191[k] = f_0 * ld_191[k];

        t_192[k] = f_0 * ld_192[k];

        t_193[k] = f_0 * ld_193[k];

        t_194[k] = f_0 * ld_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, t_200, t_201, t_202, ld_195, \
                         ld_196, ld_197, ld_198, ld_199, ld_200, ld_201, \
                         ld_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_0 * ld_195[k];

        t_196[k] = f_0 * ld_196[k];

        t_197[k] = f_0 * ld_197[k];

        t_198[k] = f_0 * ld_198[k];

        t_199[k] = f_0 * ld_199[k];

        t_200[k] = f_0 * ld_200[k];

        t_201[k] = f_0 * ld_201[k];

        t_202[k] = f_0 * ld_202[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, t_209, t_210, ld_203, \
                         ld_204, ld_205, ld_206, ld_207, ld_208, ld_209, \
                         ld_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_0 * ld_203[k];

        t_204[k] = f_0 * ld_204[k];

        t_205[k] = f_0 * ld_205[k];

        t_206[k] = f_0 * ld_206[k];

        t_207[k] = f_0 * ld_207[k];

        t_208[k] = f_0 * ld_208[k];

        t_209[k] = f_0 * ld_209[k];

        t_210[k] = f_0 * ld_210[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, ld_211, ld_212, ld_213, ld_214, \
                         ld_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_0 * ld_211[k];

        t_212[k] = f_0 * ld_212[k];

        t_213[k] = f_0 * ld_213[k];

        t_214[k] = f_0 * ld_214[k];

        t_215[k] = f_0 * ld_215[k];
    }
}

auto
compute_prim_geom_10_kd_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t id, const size_t ld,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kd_electron_repulsion_0_piece0(buffer, target, id, ld, ncols, alpha);

    compute_prim_geom_10_kd_electron_repulsion_0_piece1(buffer, target, id, ld, ncols, alpha);
}

static auto
compute_prim_geom_10_kd_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t id, const size_t ld,
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

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_73 = buffer.data(id + 73);
    const auto *id_74 = buffer.data(id + 74);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_79 = buffer.data(id + 79);
    const auto *id_80 = buffer.data(id + 80);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_85 = buffer.data(id + 85);
    const auto *id_86 = buffer.data(id + 86);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_88 = buffer.data(id + 88);
    const auto *id_89 = buffer.data(id + 89);
    const auto *id_90 = buffer.data(id + 90);
    const auto *id_91 = buffer.data(id + 91);
    const auto *id_92 = buffer.data(id + 92);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_94 = buffer.data(id + 94);
    const auto *id_95 = buffer.data(id + 95);
    const auto *id_96 = buffer.data(id + 96);
    const auto *id_97 = buffer.data(id + 97);
    const auto *id_98 = buffer.data(id + 98);
    const auto *id_99 = buffer.data(id + 99);
    const auto *id_100 = buffer.data(id + 100);
    const auto *id_101 = buffer.data(id + 101);
    const auto *id_102 = buffer.data(id + 102);
    const auto *id_103 = buffer.data(id + 103);
    const auto *id_104 = buffer.data(id + 104);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_106 = buffer.data(id + 106);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_108 = buffer.data(id + 108);
    const auto *id_109 = buffer.data(id + 109);
    const auto *id_110 = buffer.data(id + 110);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_115 = buffer.data(id + 115);
    const auto *id_116 = buffer.data(id + 116);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_119 = buffer.data(id + 119);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_121 = buffer.data(id + 121);
    const auto *id_122 = buffer.data(id + 122);
    const auto *id_123 = buffer.data(id + 123);
    const auto *id_124 = buffer.data(id + 124);
    const auto *id_125 = buffer.data(id + 125);

    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);
    const auto *ld_60 = buffer.data(ld + 60);
    const auto *ld_61 = buffer.data(ld + 61);
    const auto *ld_62 = buffer.data(ld + 62);
    const auto *ld_63 = buffer.data(ld + 63);
    const auto *ld_64 = buffer.data(ld + 64);
    const auto *ld_65 = buffer.data(ld + 65);
    const auto *ld_66 = buffer.data(ld + 66);
    const auto *ld_67 = buffer.data(ld + 67);
    const auto *ld_68 = buffer.data(ld + 68);
    const auto *ld_69 = buffer.data(ld + 69);
    const auto *ld_70 = buffer.data(ld + 70);
    const auto *ld_71 = buffer.data(ld + 71);
    const auto *ld_72 = buffer.data(ld + 72);
    const auto *ld_73 = buffer.data(ld + 73);
    const auto *ld_74 = buffer.data(ld + 74);
    const auto *ld_75 = buffer.data(ld + 75);
    const auto *ld_76 = buffer.data(ld + 76);
    const auto *ld_77 = buffer.data(ld + 77);
    const auto *ld_78 = buffer.data(ld + 78);
    const auto *ld_79 = buffer.data(ld + 79);
    const auto *ld_80 = buffer.data(ld + 80);
    const auto *ld_81 = buffer.data(ld + 81);
    const auto *ld_82 = buffer.data(ld + 82);
    const auto *ld_83 = buffer.data(ld + 83);
    const auto *ld_90 = buffer.data(ld + 90);
    const auto *ld_91 = buffer.data(ld + 91);
    const auto *ld_92 = buffer.data(ld + 92);
    const auto *ld_93 = buffer.data(ld + 93);
    const auto *ld_94 = buffer.data(ld + 94);
    const auto *ld_95 = buffer.data(ld + 95);
    const auto *ld_96 = buffer.data(ld + 96);
    const auto *ld_97 = buffer.data(ld + 97);
    const auto *ld_98 = buffer.data(ld + 98);
    const auto *ld_99 = buffer.data(ld + 99);
    const auto *ld_100 = buffer.data(ld + 100);
    const auto *ld_101 = buffer.data(ld + 101);
    const auto *ld_102 = buffer.data(ld + 102);
    const auto *ld_103 = buffer.data(ld + 103);
    const auto *ld_104 = buffer.data(ld + 104);
    const auto *ld_105 = buffer.data(ld + 105);
    const auto *ld_106 = buffer.data(ld + 106);
    const auto *ld_107 = buffer.data(ld + 107);
    const auto *ld_108 = buffer.data(ld + 108);
    const auto *ld_109 = buffer.data(ld + 109);
    const auto *ld_110 = buffer.data(ld + 110);
    const auto *ld_111 = buffer.data(ld + 111);
    const auto *ld_112 = buffer.data(ld + 112);
    const auto *ld_113 = buffer.data(ld + 113);
    const auto *ld_114 = buffer.data(ld + 114);
    const auto *ld_115 = buffer.data(ld + 115);
    const auto *ld_116 = buffer.data(ld + 116);
    const auto *ld_117 = buffer.data(ld + 117);
    const auto *ld_118 = buffer.data(ld + 118);
    const auto *ld_119 = buffer.data(ld + 119);
    const auto *ld_126 = buffer.data(ld + 126);
    const auto *ld_127 = buffer.data(ld + 127);
    const auto *ld_128 = buffer.data(ld + 128);
    const auto *ld_129 = buffer.data(ld + 129);
    const auto *ld_130 = buffer.data(ld + 130);
    const auto *ld_131 = buffer.data(ld + 131);
    const auto *ld_132 = buffer.data(ld + 132);
    const auto *ld_133 = buffer.data(ld + 133);
    const auto *ld_134 = buffer.data(ld + 134);
    const auto *ld_135 = buffer.data(ld + 135);
    const auto *ld_136 = buffer.data(ld + 136);
    const auto *ld_137 = buffer.data(ld + 137);
    const auto *ld_138 = buffer.data(ld + 138);
    const auto *ld_139 = buffer.data(ld + 139);
    const auto *ld_140 = buffer.data(ld + 140);
    const auto *ld_141 = buffer.data(ld + 141);
    const auto *ld_142 = buffer.data(ld + 142);
    const auto *ld_143 = buffer.data(ld + 143);
    const auto *ld_144 = buffer.data(ld + 144);
    const auto *ld_145 = buffer.data(ld + 145);
    const auto *ld_146 = buffer.data(ld + 146);
    const auto *ld_147 = buffer.data(ld + 147);
    const auto *ld_148 = buffer.data(ld + 148);
    const auto *ld_149 = buffer.data(ld + 149);
    const auto *ld_150 = buffer.data(ld + 150);
    const auto *ld_151 = buffer.data(ld + 151);
    const auto *ld_152 = buffer.data(ld + 152);
    const auto *ld_153 = buffer.data(ld + 153);
    const auto *ld_154 = buffer.data(ld + 154);
    const auto *ld_155 = buffer.data(ld + 155);
    const auto *ld_156 = buffer.data(ld + 156);
    const auto *ld_157 = buffer.data(ld + 157);
    const auto *ld_158 = buffer.data(ld + 158);
    const auto *ld_159 = buffer.data(ld + 159);
    const auto *ld_160 = buffer.data(ld + 160);
    const auto *ld_161 = buffer.data(ld + 161);
    const auto *ld_168 = buffer.data(ld + 168);
    const auto *ld_169 = buffer.data(ld + 169);
    const auto *ld_170 = buffer.data(ld + 170);
    const auto *ld_171 = buffer.data(ld + 171);
    const auto *ld_172 = buffer.data(ld + 172);
    const auto *ld_173 = buffer.data(ld + 173);
    const auto *ld_174 = buffer.data(ld + 174);
    const auto *ld_175 = buffer.data(ld + 175);
    const auto *ld_176 = buffer.data(ld + 176);
    const auto *ld_177 = buffer.data(ld + 177);
    const auto *ld_178 = buffer.data(ld + 178);
    const auto *ld_179 = buffer.data(ld + 179);
    const auto *ld_180 = buffer.data(ld + 180);
    const auto *ld_181 = buffer.data(ld + 181);
    const auto *ld_182 = buffer.data(ld + 182);
    const auto *ld_183 = buffer.data(ld + 183);
    const auto *ld_184 = buffer.data(ld + 184);
    const auto *ld_185 = buffer.data(ld + 185);
    const auto *ld_186 = buffer.data(ld + 186);
    const auto *ld_187 = buffer.data(ld + 187);
    const auto *ld_188 = buffer.data(ld + 188);
    const auto *ld_189 = buffer.data(ld + 189);
    const auto *ld_190 = buffer.data(ld + 190);
    const auto *ld_191 = buffer.data(ld + 191);
    const auto *ld_192 = buffer.data(ld + 192);
    const auto *ld_193 = buffer.data(ld + 193);
    const auto *ld_194 = buffer.data(ld + 194);
    const auto *ld_195 = buffer.data(ld + 195);
    const auto *ld_196 = buffer.data(ld + 196);
    const auto *ld_197 = buffer.data(ld + 197);
    const auto *ld_198 = buffer.data(ld + 198);
    const auto *ld_199 = buffer.data(ld + 199);
    const auto *ld_200 = buffer.data(ld + 200);
    const auto *ld_201 = buffer.data(ld + 201);
    const auto *ld_202 = buffer.data(ld + 202);
    const auto *ld_203 = buffer.data(ld + 203);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, id_0, ld_6, ld_7, ld_8, ld_9, \
                         ld_10, ld_11, ld_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ld_6[k];

        t_1[k] = f_0 * ld_7[k];

        t_2[k] = f_0 * ld_8[k];

        t_3[k] = f_0 * ld_9[k];

        t_4[k] = f_0 * ld_10[k];

        t_5[k] = f_0 * ld_11[k];

        t_6[k] = -id_0[k]
                 + f_0 * ld_18[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, id_1, id_2, id_3, id_4, id_5, ld_19, \
                         ld_20, ld_21, ld_22, ld_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -id_1[k]
                 + f_0 * ld_19[k];

        t_8[k] = -id_2[k]
                 + f_0 * ld_20[k];

        t_9[k] = -id_3[k]
                 + f_0 * ld_21[k];

        t_10[k] = -id_4[k]
                  + f_0 * ld_22[k];

        t_11[k] = -id_5[k]
                  + f_0 * ld_23[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, t_18, id_6, ld_24, ld_25, ld_26, \
                         ld_27, ld_28, ld_29, ld_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * ld_24[k];

        t_13[k] = f_0 * ld_25[k];

        t_14[k] = f_0 * ld_26[k];

        t_15[k] = f_0 * ld_27[k];

        t_16[k] = f_0 * ld_28[k];

        t_17[k] = f_0 * ld_29[k];

        t_18[k] = -2.0 * id_6[k]
                  + f_0 * ld_36[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, id_7, id_8, id_9, id_10, id_11, ld_37, \
                         ld_38, ld_39, ld_40, ld_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -2.0 * id_7[k]
                  + f_0 * ld_37[k];

        t_20[k] = -2.0 * id_8[k]
                  + f_0 * ld_38[k];

        t_21[k] = -2.0 * id_9[k]
                  + f_0 * ld_39[k];

        t_22[k] = -2.0 * id_10[k]
                  + f_0 * ld_40[k];

        t_23[k] = -2.0 * id_11[k]
                  + f_0 * ld_41[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, id_12, id_13, id_14, id_15, id_16, \
                         ld_42, ld_43, ld_44, ld_45, ld_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -id_12[k]
                  + f_0 * ld_42[k];

        t_25[k] = -id_13[k]
                  + f_0 * ld_43[k];

        t_26[k] = -id_14[k]
                  + f_0 * ld_44[k];

        t_27[k] = -id_15[k]
                  + f_0 * ld_45[k];

        t_28[k] = -id_16[k]
                  + f_0 * ld_46[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, t_35, id_17, ld_47, ld_48, ld_49, \
                         ld_50, ld_51, ld_52, ld_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -id_17[k]
                  + f_0 * ld_47[k];

        t_30[k] = f_0 * ld_48[k];

        t_31[k] = f_0 * ld_49[k];

        t_32[k] = f_0 * ld_50[k];

        t_33[k] = f_0 * ld_51[k];

        t_34[k] = f_0 * ld_52[k];

        t_35[k] = f_0 * ld_53[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, id_18, id_19, id_20, id_21, id_22, \
                         ld_60, ld_61, ld_62, ld_63, ld_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -3.0 * id_18[k]
                  + f_0 * ld_60[k];

        t_37[k] = -3.0 * id_19[k]
                  + f_0 * ld_61[k];

        t_38[k] = -3.0 * id_20[k]
                  + f_0 * ld_62[k];

        t_39[k] = -3.0 * id_21[k]
                  + f_0 * ld_63[k];

        t_40[k] = -3.0 * id_22[k]
                  + f_0 * ld_64[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, id_23, id_24, id_25, id_26, id_27, \
                         ld_65, ld_66, ld_67, ld_68, ld_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -3.0 * id_23[k]
                  + f_0 * ld_65[k];

        t_42[k] = -2.0 * id_24[k]
                  + f_0 * ld_66[k];

        t_43[k] = -2.0 * id_25[k]
                  + f_0 * ld_67[k];

        t_44[k] = -2.0 * id_26[k]
                  + f_0 * ld_68[k];

        t_45[k] = -2.0 * id_27[k]
                  + f_0 * ld_69[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, id_28, id_29, id_30, id_31, id_32, \
                         ld_70, ld_71, ld_72, ld_73, ld_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -2.0 * id_28[k]
                  + f_0 * ld_70[k];

        t_47[k] = -2.0 * id_29[k]
                  + f_0 * ld_71[k];

        t_48[k] = -id_30[k]
                  + f_0 * ld_72[k];

        t_49[k] = -id_31[k]
                  + f_0 * ld_73[k];

        t_50[k] = -id_32[k]
                  + f_0 * ld_74[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, id_33, id_34, id_35, ld_75, \
                         ld_76, ld_77, ld_78, ld_79, ld_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -id_33[k]
                  + f_0 * ld_75[k];

        t_52[k] = -id_34[k]
                  + f_0 * ld_76[k];

        t_53[k] = -id_35[k]
                  + f_0 * ld_77[k];

        t_54[k] = f_0 * ld_78[k];

        t_55[k] = f_0 * ld_79[k];

        t_56[k] = f_0 * ld_80[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, t_62, id_36, id_37, id_38, ld_81, \
                         ld_82, ld_83, ld_90, ld_91, ld_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_0 * ld_81[k];

        t_58[k] = f_0 * ld_82[k];

        t_59[k] = f_0 * ld_83[k];

        t_60[k] = -4.0 * id_36[k]
                  + f_0 * ld_90[k];

        t_61[k] = -4.0 * id_37[k]
                  + f_0 * ld_91[k];

        t_62[k] = -4.0 * id_38[k]
                  + f_0 * ld_92[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, id_39, id_40, id_41, id_42, id_43, \
                         ld_93, ld_94, ld_95, ld_96, ld_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -4.0 * id_39[k]
                  + f_0 * ld_93[k];

        t_64[k] = -4.0 * id_40[k]
                  + f_0 * ld_94[k];

        t_65[k] = -4.0 * id_41[k]
                  + f_0 * ld_95[k];

        t_66[k] = -3.0 * id_42[k]
                  + f_0 * ld_96[k];

        t_67[k] = -3.0 * id_43[k]
                  + f_0 * ld_97[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, id_44, id_45, id_46, id_47, id_48, \
                         ld_98, ld_99, ld_100, ld_101, ld_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -3.0 * id_44[k]
                  + f_0 * ld_98[k];

        t_69[k] = -3.0 * id_45[k]
                  + f_0 * ld_99[k];

        t_70[k] = -3.0 * id_46[k]
                  + f_0 * ld_100[k];

        t_71[k] = -3.0 * id_47[k]
                  + f_0 * ld_101[k];

        t_72[k] = -2.0 * id_48[k]
                  + f_0 * ld_102[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, id_49, id_50, id_51, id_52, id_53, \
                         ld_103, ld_104, ld_105, ld_106, ld_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = -2.0 * id_49[k]
                  + f_0 * ld_103[k];

        t_74[k] = -2.0 * id_50[k]
                  + f_0 * ld_104[k];

        t_75[k] = -2.0 * id_51[k]
                  + f_0 * ld_105[k];

        t_76[k] = -2.0 * id_52[k]
                  + f_0 * ld_106[k];

        t_77[k] = -2.0 * id_53[k]
                  + f_0 * ld_107[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, id_54, id_55, id_56, id_57, id_58, \
                         ld_108, ld_109, ld_110, ld_111, ld_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = -id_54[k]
                  + f_0 * ld_108[k];

        t_79[k] = -id_55[k]
                  + f_0 * ld_109[k];

        t_80[k] = -id_56[k]
                  + f_0 * ld_110[k];

        t_81[k] = -id_57[k]
                  + f_0 * ld_111[k];

        t_82[k] = -id_58[k]
                  + f_0 * ld_112[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, id_59, ld_113, ld_114, \
                         ld_115, ld_116, ld_117, ld_118, ld_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = -id_59[k]
                  + f_0 * ld_113[k];

        t_84[k] = f_0 * ld_114[k];

        t_85[k] = f_0 * ld_115[k];

        t_86[k] = f_0 * ld_116[k];

        t_87[k] = f_0 * ld_117[k];

        t_88[k] = f_0 * ld_118[k];

        t_89[k] = f_0 * ld_119[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, id_60, id_61, id_62, id_63, id_64, \
                         ld_126, ld_127, ld_128, ld_129, ld_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -5.0 * id_60[k]
                  + f_0 * ld_126[k];

        t_91[k] = -5.0 * id_61[k]
                  + f_0 * ld_127[k];

        t_92[k] = -5.0 * id_62[k]
                  + f_0 * ld_128[k];

        t_93[k] = -5.0 * id_63[k]
                  + f_0 * ld_129[k];

        t_94[k] = -5.0 * id_64[k]
                  + f_0 * ld_130[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, id_65, id_66, id_67, id_68, id_69, \
                         ld_131, ld_132, ld_133, ld_134, ld_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -5.0 * id_65[k]
                  + f_0 * ld_131[k];

        t_96[k] = -4.0 * id_66[k]
                  + f_0 * ld_132[k];

        t_97[k] = -4.0 * id_67[k]
                  + f_0 * ld_133[k];

        t_98[k] = -4.0 * id_68[k]
                  + f_0 * ld_134[k];

        t_99[k] = -4.0 * id_69[k]
                  + f_0 * ld_135[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, id_70, id_71, id_72, id_73, id_74, \
                         ld_136, ld_137, ld_138, ld_139, ld_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -4.0 * id_70[k]
                   + f_0 * ld_136[k];

        t_101[k] = -4.0 * id_71[k]
                   + f_0 * ld_137[k];

        t_102[k] = -3.0 * id_72[k]
                   + f_0 * ld_138[k];

        t_103[k] = -3.0 * id_73[k]
                   + f_0 * ld_139[k];

        t_104[k] = -3.0 * id_74[k]
                   + f_0 * ld_140[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, id_75, id_76, id_77, id_78, id_79, \
                         ld_141, ld_142, ld_143, ld_144, ld_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -3.0 * id_75[k]
                   + f_0 * ld_141[k];

        t_106[k] = -3.0 * id_76[k]
                   + f_0 * ld_142[k];

        t_107[k] = -3.0 * id_77[k]
                   + f_0 * ld_143[k];

        t_108[k] = -2.0 * id_78[k]
                   + f_0 * ld_144[k];

        t_109[k] = -2.0 * id_79[k]
                   + f_0 * ld_145[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, id_80, id_81, id_82, id_83, id_84, \
                         ld_146, ld_147, ld_148, ld_149, ld_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * id_80[k]
                   + f_0 * ld_146[k];

        t_111[k] = -2.0 * id_81[k]
                   + f_0 * ld_147[k];

        t_112[k] = -2.0 * id_82[k]
                   + f_0 * ld_148[k];

        t_113[k] = -2.0 * id_83[k]
                   + f_0 * ld_149[k];

        t_114[k] = -id_84[k]
                   + f_0 * ld_150[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, id_85, id_86, id_87, id_88, id_89, \
                         ld_151, ld_152, ld_153, ld_154, ld_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -id_85[k]
                   + f_0 * ld_151[k];

        t_116[k] = -id_86[k]
                   + f_0 * ld_152[k];

        t_117[k] = -id_87[k]
                   + f_0 * ld_153[k];

        t_118[k] = -id_88[k]
                   + f_0 * ld_154[k];

        t_119[k] = -id_89[k]
                   + f_0 * ld_155[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, t_126, id_90, ld_156, \
                         ld_157, ld_158, ld_159, ld_160, ld_161, \
                         ld_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_0 * ld_156[k];

        t_121[k] = f_0 * ld_157[k];

        t_122[k] = f_0 * ld_158[k];

        t_123[k] = f_0 * ld_159[k];

        t_124[k] = f_0 * ld_160[k];

        t_125[k] = f_0 * ld_161[k];

        t_126[k] = -6.0 * id_90[k]
                   + f_0 * ld_168[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, id_91, id_92, id_93, id_94, id_95, \
                         ld_169, ld_170, ld_171, ld_172, ld_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -6.0 * id_91[k]
                   + f_0 * ld_169[k];

        t_128[k] = -6.0 * id_92[k]
                   + f_0 * ld_170[k];

        t_129[k] = -6.0 * id_93[k]
                   + f_0 * ld_171[k];

        t_130[k] = -6.0 * id_94[k]
                   + f_0 * ld_172[k];

        t_131[k] = -6.0 * id_95[k]
                   + f_0 * ld_173[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, id_96, id_97, id_98, id_99, \
                         id_100, ld_174, ld_175, ld_176, ld_177, \
                         ld_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -5.0 * id_96[k]
                   + f_0 * ld_174[k];

        t_133[k] = -5.0 * id_97[k]
                   + f_0 * ld_175[k];

        t_134[k] = -5.0 * id_98[k]
                   + f_0 * ld_176[k];

        t_135[k] = -5.0 * id_99[k]
                   + f_0 * ld_177[k];

        t_136[k] = -5.0 * id_100[k]
                   + f_0 * ld_178[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, id_101, id_102, id_103, id_104, \
                         id_105, ld_179, ld_180, ld_181, ld_182, \
                         ld_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -5.0 * id_101[k]
                   + f_0 * ld_179[k];

        t_138[k] = -4.0 * id_102[k]
                   + f_0 * ld_180[k];

        t_139[k] = -4.0 * id_103[k]
                   + f_0 * ld_181[k];

        t_140[k] = -4.0 * id_104[k]
                   + f_0 * ld_182[k];

        t_141[k] = -4.0 * id_105[k]
                   + f_0 * ld_183[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, id_106, id_107, id_108, id_109, \
                         id_110, ld_184, ld_185, ld_186, ld_187, \
                         ld_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -4.0 * id_106[k]
                   + f_0 * ld_184[k];

        t_143[k] = -4.0 * id_107[k]
                   + f_0 * ld_185[k];

        t_144[k] = -3.0 * id_108[k]
                   + f_0 * ld_186[k];

        t_145[k] = -3.0 * id_109[k]
                   + f_0 * ld_187[k];

        t_146[k] = -3.0 * id_110[k]
                   + f_0 * ld_188[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, id_111, id_112, id_113, id_114, \
                         id_115, ld_189, ld_190, ld_191, ld_192, \
                         ld_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -3.0 * id_111[k]
                   + f_0 * ld_189[k];

        t_148[k] = -3.0 * id_112[k]
                   + f_0 * ld_190[k];

        t_149[k] = -3.0 * id_113[k]
                   + f_0 * ld_191[k];

        t_150[k] = -2.0 * id_114[k]
                   + f_0 * ld_192[k];

        t_151[k] = -2.0 * id_115[k]
                   + f_0 * ld_193[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, id_116, id_117, id_118, id_119, \
                         id_120, ld_194, ld_195, ld_196, ld_197, \
                         ld_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -2.0 * id_116[k]
                   + f_0 * ld_194[k];

        t_153[k] = -2.0 * id_117[k]
                   + f_0 * ld_195[k];

        t_154[k] = -2.0 * id_118[k]
                   + f_0 * ld_196[k];

        t_155[k] = -2.0 * id_119[k]
                   + f_0 * ld_197[k];

        t_156[k] = -id_120[k]
                   + f_0 * ld_198[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, id_121, id_122, id_123, id_124, \
                         id_125, ld_199, ld_200, ld_201, ld_202, \
                         ld_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -id_121[k]
                   + f_0 * ld_199[k];

        t_158[k] = -id_122[k]
                   + f_0 * ld_200[k];

        t_159[k] = -id_123[k]
                   + f_0 * ld_201[k];

        t_160[k] = -id_124[k]
                   + f_0 * ld_202[k];

        t_161[k] = -id_125[k]
                   + f_0 * ld_203[k];
    }
}

static auto
compute_prim_geom_10_kd_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t id, const size_t ld,
                                                    const size_t ncols,
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

    const auto *id_126 = buffer.data(id + 126);
    const auto *id_127 = buffer.data(id + 127);
    const auto *id_128 = buffer.data(id + 128);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_130 = buffer.data(id + 130);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_132 = buffer.data(id + 132);
    const auto *id_133 = buffer.data(id + 133);
    const auto *id_134 = buffer.data(id + 134);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_138 = buffer.data(id + 138);
    const auto *id_139 = buffer.data(id + 139);
    const auto *id_140 = buffer.data(id + 140);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_144 = buffer.data(id + 144);
    const auto *id_145 = buffer.data(id + 145);
    const auto *id_146 = buffer.data(id + 146);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);
    const auto *id_150 = buffer.data(id + 150);
    const auto *id_151 = buffer.data(id + 151);
    const auto *id_152 = buffer.data(id + 152);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_156 = buffer.data(id + 156);
    const auto *id_157 = buffer.data(id + 157);
    const auto *id_158 = buffer.data(id + 158);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_163 = buffer.data(id + 163);
    const auto *id_164 = buffer.data(id + 164);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_166 = buffer.data(id + 166);
    const auto *id_167 = buffer.data(id + 167);

    const auto *ld_204 = buffer.data(ld + 204);
    const auto *ld_205 = buffer.data(ld + 205);
    const auto *ld_206 = buffer.data(ld + 206);
    const auto *ld_207 = buffer.data(ld + 207);
    const auto *ld_208 = buffer.data(ld + 208);
    const auto *ld_209 = buffer.data(ld + 209);
    const auto *ld_216 = buffer.data(ld + 216);
    const auto *ld_217 = buffer.data(ld + 217);
    const auto *ld_218 = buffer.data(ld + 218);
    const auto *ld_219 = buffer.data(ld + 219);
    const auto *ld_220 = buffer.data(ld + 220);
    const auto *ld_221 = buffer.data(ld + 221);
    const auto *ld_222 = buffer.data(ld + 222);
    const auto *ld_223 = buffer.data(ld + 223);
    const auto *ld_224 = buffer.data(ld + 224);
    const auto *ld_225 = buffer.data(ld + 225);
    const auto *ld_226 = buffer.data(ld + 226);
    const auto *ld_227 = buffer.data(ld + 227);
    const auto *ld_228 = buffer.data(ld + 228);
    const auto *ld_229 = buffer.data(ld + 229);
    const auto *ld_230 = buffer.data(ld + 230);
    const auto *ld_231 = buffer.data(ld + 231);
    const auto *ld_232 = buffer.data(ld + 232);
    const auto *ld_233 = buffer.data(ld + 233);
    const auto *ld_234 = buffer.data(ld + 234);
    const auto *ld_235 = buffer.data(ld + 235);
    const auto *ld_236 = buffer.data(ld + 236);
    const auto *ld_237 = buffer.data(ld + 237);
    const auto *ld_238 = buffer.data(ld + 238);
    const auto *ld_239 = buffer.data(ld + 239);
    const auto *ld_240 = buffer.data(ld + 240);
    const auto *ld_241 = buffer.data(ld + 241);
    const auto *ld_242 = buffer.data(ld + 242);
    const auto *ld_243 = buffer.data(ld + 243);
    const auto *ld_244 = buffer.data(ld + 244);
    const auto *ld_245 = buffer.data(ld + 245);
    const auto *ld_246 = buffer.data(ld + 246);
    const auto *ld_247 = buffer.data(ld + 247);
    const auto *ld_248 = buffer.data(ld + 248);
    const auto *ld_249 = buffer.data(ld + 249);
    const auto *ld_250 = buffer.data(ld + 250);
    const auto *ld_251 = buffer.data(ld + 251);
    const auto *ld_252 = buffer.data(ld + 252);
    const auto *ld_253 = buffer.data(ld + 253);
    const auto *ld_254 = buffer.data(ld + 254);
    const auto *ld_255 = buffer.data(ld + 255);
    const auto *ld_256 = buffer.data(ld + 256);
    const auto *ld_257 = buffer.data(ld + 257);
    const auto *ld_258 = buffer.data(ld + 258);
    const auto *ld_259 = buffer.data(ld + 259);
    const auto *ld_260 = buffer.data(ld + 260);
    const auto *ld_261 = buffer.data(ld + 261);
    const auto *ld_262 = buffer.data(ld + 262);
    const auto *ld_263 = buffer.data(ld + 263);

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, t_167, t_168, id_126, ld_204, \
                         ld_205, ld_206, ld_207, ld_208, ld_209, \
                         ld_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_0 * ld_204[k];

        t_163[k] = f_0 * ld_205[k];

        t_164[k] = f_0 * ld_206[k];

        t_165[k] = f_0 * ld_207[k];

        t_166[k] = f_0 * ld_208[k];

        t_167[k] = f_0 * ld_209[k];

        t_168[k] = -7.0 * id_126[k]
                   + f_0 * ld_216[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, id_127, id_128, id_129, id_130, \
                         id_131, ld_217, ld_218, ld_219, ld_220, \
                         ld_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -7.0 * id_127[k]
                   + f_0 * ld_217[k];

        t_170[k] = -7.0 * id_128[k]
                   + f_0 * ld_218[k];

        t_171[k] = -7.0 * id_129[k]
                   + f_0 * ld_219[k];

        t_172[k] = -7.0 * id_130[k]
                   + f_0 * ld_220[k];

        t_173[k] = -7.0 * id_131[k]
                   + f_0 * ld_221[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, id_132, id_133, id_134, id_135, \
                         id_136, ld_222, ld_223, ld_224, ld_225, \
                         ld_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -6.0 * id_132[k]
                   + f_0 * ld_222[k];

        t_175[k] = -6.0 * id_133[k]
                   + f_0 * ld_223[k];

        t_176[k] = -6.0 * id_134[k]
                   + f_0 * ld_224[k];

        t_177[k] = -6.0 * id_135[k]
                   + f_0 * ld_225[k];

        t_178[k] = -6.0 * id_136[k]
                   + f_0 * ld_226[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, id_137, id_138, id_139, id_140, \
                         id_141, ld_227, ld_228, ld_229, ld_230, \
                         ld_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -6.0 * id_137[k]
                   + f_0 * ld_227[k];

        t_180[k] = -5.0 * id_138[k]
                   + f_0 * ld_228[k];

        t_181[k] = -5.0 * id_139[k]
                   + f_0 * ld_229[k];

        t_182[k] = -5.0 * id_140[k]
                   + f_0 * ld_230[k];

        t_183[k] = -5.0 * id_141[k]
                   + f_0 * ld_231[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, id_142, id_143, id_144, id_145, \
                         id_146, ld_232, ld_233, ld_234, ld_235, \
                         ld_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -5.0 * id_142[k]
                   + f_0 * ld_232[k];

        t_185[k] = -5.0 * id_143[k]
                   + f_0 * ld_233[k];

        t_186[k] = -4.0 * id_144[k]
                   + f_0 * ld_234[k];

        t_187[k] = -4.0 * id_145[k]
                   + f_0 * ld_235[k];

        t_188[k] = -4.0 * id_146[k]
                   + f_0 * ld_236[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, id_147, id_148, id_149, id_150, \
                         id_151, ld_237, ld_238, ld_239, ld_240, \
                         ld_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = -4.0 * id_147[k]
                   + f_0 * ld_237[k];

        t_190[k] = -4.0 * id_148[k]
                   + f_0 * ld_238[k];

        t_191[k] = -4.0 * id_149[k]
                   + f_0 * ld_239[k];

        t_192[k] = -3.0 * id_150[k]
                   + f_0 * ld_240[k];

        t_193[k] = -3.0 * id_151[k]
                   + f_0 * ld_241[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, id_152, id_153, id_154, id_155, \
                         id_156, ld_242, ld_243, ld_244, ld_245, \
                         ld_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = -3.0 * id_152[k]
                   + f_0 * ld_242[k];

        t_195[k] = -3.0 * id_153[k]
                   + f_0 * ld_243[k];

        t_196[k] = -3.0 * id_154[k]
                   + f_0 * ld_244[k];

        t_197[k] = -3.0 * id_155[k]
                   + f_0 * ld_245[k];

        t_198[k] = -2.0 * id_156[k]
                   + f_0 * ld_246[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, id_157, id_158, id_159, id_160, \
                         id_161, ld_247, ld_248, ld_249, ld_250, \
                         ld_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = -2.0 * id_157[k]
                   + f_0 * ld_247[k];

        t_200[k] = -2.0 * id_158[k]
                   + f_0 * ld_248[k];

        t_201[k] = -2.0 * id_159[k]
                   + f_0 * ld_249[k];

        t_202[k] = -2.0 * id_160[k]
                   + f_0 * ld_250[k];

        t_203[k] = -2.0 * id_161[k]
                   + f_0 * ld_251[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, id_162, id_163, id_164, id_165, \
                         id_166, ld_252, ld_253, ld_254, ld_255, \
                         ld_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = -id_162[k]
                   + f_0 * ld_252[k];

        t_205[k] = -id_163[k]
                   + f_0 * ld_253[k];

        t_206[k] = -id_164[k]
                   + f_0 * ld_254[k];

        t_207[k] = -id_165[k]
                   + f_0 * ld_255[k];

        t_208[k] = -id_166[k]
                   + f_0 * ld_256[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, t_214, t_215, id_167, ld_257, \
                         ld_258, ld_259, ld_260, ld_261, ld_262, \
                         ld_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = -id_167[k]
                   + f_0 * ld_257[k];

        t_210[k] = f_0 * ld_258[k];

        t_211[k] = f_0 * ld_259[k];

        t_212[k] = f_0 * ld_260[k];

        t_213[k] = f_0 * ld_261[k];

        t_214[k] = f_0 * ld_262[k];

        t_215[k] = f_0 * ld_263[k];
    }
}

auto
compute_prim_geom_10_kd_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t id, const size_t ld,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kd_electron_repulsion_1_piece0(buffer, target, id, ld, ncols, alpha);

    compute_prim_geom_10_kd_electron_repulsion_1_piece1(buffer, target, id, ld, ncols, alpha);
}

static auto
compute_prim_geom_10_kd_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t id, const size_t ld,
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

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_73 = buffer.data(id + 73);
    const auto *id_74 = buffer.data(id + 74);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_79 = buffer.data(id + 79);
    const auto *id_80 = buffer.data(id + 80);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_85 = buffer.data(id + 85);
    const auto *id_86 = buffer.data(id + 86);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_88 = buffer.data(id + 88);
    const auto *id_89 = buffer.data(id + 89);
    const auto *id_90 = buffer.data(id + 90);
    const auto *id_91 = buffer.data(id + 91);
    const auto *id_92 = buffer.data(id + 92);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_94 = buffer.data(id + 94);
    const auto *id_95 = buffer.data(id + 95);
    const auto *id_96 = buffer.data(id + 96);
    const auto *id_97 = buffer.data(id + 97);
    const auto *id_98 = buffer.data(id + 98);
    const auto *id_99 = buffer.data(id + 99);
    const auto *id_100 = buffer.data(id + 100);
    const auto *id_101 = buffer.data(id + 101);
    const auto *id_102 = buffer.data(id + 102);
    const auto *id_103 = buffer.data(id + 103);
    const auto *id_104 = buffer.data(id + 104);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_106 = buffer.data(id + 106);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_108 = buffer.data(id + 108);
    const auto *id_109 = buffer.data(id + 109);
    const auto *id_110 = buffer.data(id + 110);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_115 = buffer.data(id + 115);
    const auto *id_116 = buffer.data(id + 116);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_119 = buffer.data(id + 119);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_121 = buffer.data(id + 121);

    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);
    const auto *ld_54 = buffer.data(ld + 54);
    const auto *ld_55 = buffer.data(ld + 55);
    const auto *ld_56 = buffer.data(ld + 56);
    const auto *ld_57 = buffer.data(ld + 57);
    const auto *ld_58 = buffer.data(ld + 58);
    const auto *ld_59 = buffer.data(ld + 59);
    const auto *ld_66 = buffer.data(ld + 66);
    const auto *ld_67 = buffer.data(ld + 67);
    const auto *ld_68 = buffer.data(ld + 68);
    const auto *ld_69 = buffer.data(ld + 69);
    const auto *ld_70 = buffer.data(ld + 70);
    const auto *ld_71 = buffer.data(ld + 71);
    const auto *ld_72 = buffer.data(ld + 72);
    const auto *ld_73 = buffer.data(ld + 73);
    const auto *ld_74 = buffer.data(ld + 74);
    const auto *ld_75 = buffer.data(ld + 75);
    const auto *ld_76 = buffer.data(ld + 76);
    const auto *ld_77 = buffer.data(ld + 77);
    const auto *ld_78 = buffer.data(ld + 78);
    const auto *ld_79 = buffer.data(ld + 79);
    const auto *ld_80 = buffer.data(ld + 80);
    const auto *ld_81 = buffer.data(ld + 81);
    const auto *ld_82 = buffer.data(ld + 82);
    const auto *ld_83 = buffer.data(ld + 83);
    const auto *ld_84 = buffer.data(ld + 84);
    const auto *ld_85 = buffer.data(ld + 85);
    const auto *ld_86 = buffer.data(ld + 86);
    const auto *ld_87 = buffer.data(ld + 87);
    const auto *ld_88 = buffer.data(ld + 88);
    const auto *ld_89 = buffer.data(ld + 89);
    const auto *ld_96 = buffer.data(ld + 96);
    const auto *ld_97 = buffer.data(ld + 97);
    const auto *ld_98 = buffer.data(ld + 98);
    const auto *ld_99 = buffer.data(ld + 99);
    const auto *ld_100 = buffer.data(ld + 100);
    const auto *ld_101 = buffer.data(ld + 101);
    const auto *ld_102 = buffer.data(ld + 102);
    const auto *ld_103 = buffer.data(ld + 103);
    const auto *ld_104 = buffer.data(ld + 104);
    const auto *ld_105 = buffer.data(ld + 105);
    const auto *ld_106 = buffer.data(ld + 106);
    const auto *ld_107 = buffer.data(ld + 107);
    const auto *ld_108 = buffer.data(ld + 108);
    const auto *ld_109 = buffer.data(ld + 109);
    const auto *ld_110 = buffer.data(ld + 110);
    const auto *ld_111 = buffer.data(ld + 111);
    const auto *ld_112 = buffer.data(ld + 112);
    const auto *ld_113 = buffer.data(ld + 113);
    const auto *ld_114 = buffer.data(ld + 114);
    const auto *ld_115 = buffer.data(ld + 115);
    const auto *ld_116 = buffer.data(ld + 116);
    const auto *ld_117 = buffer.data(ld + 117);
    const auto *ld_118 = buffer.data(ld + 118);
    const auto *ld_119 = buffer.data(ld + 119);
    const auto *ld_120 = buffer.data(ld + 120);
    const auto *ld_121 = buffer.data(ld + 121);
    const auto *ld_122 = buffer.data(ld + 122);
    const auto *ld_123 = buffer.data(ld + 123);
    const auto *ld_124 = buffer.data(ld + 124);
    const auto *ld_125 = buffer.data(ld + 125);
    const auto *ld_132 = buffer.data(ld + 132);
    const auto *ld_133 = buffer.data(ld + 133);
    const auto *ld_134 = buffer.data(ld + 134);
    const auto *ld_135 = buffer.data(ld + 135);
    const auto *ld_136 = buffer.data(ld + 136);
    const auto *ld_137 = buffer.data(ld + 137);
    const auto *ld_138 = buffer.data(ld + 138);
    const auto *ld_139 = buffer.data(ld + 139);
    const auto *ld_140 = buffer.data(ld + 140);
    const auto *ld_141 = buffer.data(ld + 141);
    const auto *ld_142 = buffer.data(ld + 142);
    const auto *ld_143 = buffer.data(ld + 143);
    const auto *ld_144 = buffer.data(ld + 144);
    const auto *ld_145 = buffer.data(ld + 145);
    const auto *ld_146 = buffer.data(ld + 146);
    const auto *ld_147 = buffer.data(ld + 147);
    const auto *ld_148 = buffer.data(ld + 148);
    const auto *ld_149 = buffer.data(ld + 149);
    const auto *ld_150 = buffer.data(ld + 150);
    const auto *ld_151 = buffer.data(ld + 151);
    const auto *ld_152 = buffer.data(ld + 152);
    const auto *ld_153 = buffer.data(ld + 153);
    const auto *ld_154 = buffer.data(ld + 154);
    const auto *ld_155 = buffer.data(ld + 155);
    const auto *ld_156 = buffer.data(ld + 156);
    const auto *ld_157 = buffer.data(ld + 157);
    const auto *ld_158 = buffer.data(ld + 158);
    const auto *ld_159 = buffer.data(ld + 159);
    const auto *ld_160 = buffer.data(ld + 160);
    const auto *ld_161 = buffer.data(ld + 161);
    const auto *ld_162 = buffer.data(ld + 162);
    const auto *ld_163 = buffer.data(ld + 163);
    const auto *ld_164 = buffer.data(ld + 164);
    const auto *ld_165 = buffer.data(ld + 165);
    const auto *ld_166 = buffer.data(ld + 166);
    const auto *ld_167 = buffer.data(ld + 167);
    const auto *ld_174 = buffer.data(ld + 174);
    const auto *ld_175 = buffer.data(ld + 175);
    const auto *ld_176 = buffer.data(ld + 176);
    const auto *ld_177 = buffer.data(ld + 177);
    const auto *ld_178 = buffer.data(ld + 178);
    const auto *ld_179 = buffer.data(ld + 179);
    const auto *ld_180 = buffer.data(ld + 180);
    const auto *ld_181 = buffer.data(ld + 181);
    const auto *ld_182 = buffer.data(ld + 182);
    const auto *ld_183 = buffer.data(ld + 183);
    const auto *ld_184 = buffer.data(ld + 184);
    const auto *ld_185 = buffer.data(ld + 185);
    const auto *ld_186 = buffer.data(ld + 186);
    const auto *ld_187 = buffer.data(ld + 187);
    const auto *ld_188 = buffer.data(ld + 188);
    const auto *ld_189 = buffer.data(ld + 189);
    const auto *ld_190 = buffer.data(ld + 190);
    const auto *ld_191 = buffer.data(ld + 191);
    const auto *ld_192 = buffer.data(ld + 192);
    const auto *ld_193 = buffer.data(ld + 193);
    const auto *ld_194 = buffer.data(ld + 194);
    const auto *ld_195 = buffer.data(ld + 195);
    const auto *ld_196 = buffer.data(ld + 196);
    const auto *ld_197 = buffer.data(ld + 197);
    const auto *ld_198 = buffer.data(ld + 198);
    const auto *ld_199 = buffer.data(ld + 199);
    const auto *ld_200 = buffer.data(ld + 200);
    const auto *ld_201 = buffer.data(ld + 201);
    const auto *ld_202 = buffer.data(ld + 202);
    const auto *ld_203 = buffer.data(ld + 203);
    const auto *ld_204 = buffer.data(ld + 204);
    const auto *ld_205 = buffer.data(ld + 205);
    const auto *ld_206 = buffer.data(ld + 206);
    const auto *ld_207 = buffer.data(ld + 207);
    const auto *ld_208 = buffer.data(ld + 208);
    const auto *ld_209 = buffer.data(ld + 209);
    const auto *ld_210 = buffer.data(ld + 210);
    const auto *ld_211 = buffer.data(ld + 211);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ld_12, ld_13, ld_14, ld_15, \
                         ld_16, ld_17, ld_24, ld_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ld_12[k];

        t_1[k] = f_0 * ld_13[k];

        t_2[k] = f_0 * ld_14[k];

        t_3[k] = f_0 * ld_15[k];

        t_4[k] = f_0 * ld_16[k];

        t_5[k] = f_0 * ld_17[k];

        t_6[k] = f_0 * ld_24[k];

        t_7[k] = f_0 * ld_25[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, id_0, id_1, ld_26, ld_27, ld_28, \
                         ld_29, ld_30, ld_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ld_26[k];

        t_9[k] = f_0 * ld_27[k];

        t_10[k] = f_0 * ld_28[k];

        t_11[k] = f_0 * ld_29[k];

        t_12[k] = -id_0[k]
                  + f_0 * ld_30[k];

        t_13[k] = -id_1[k]
                  + f_0 * ld_31[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, id_2, id_3, id_4, id_5, ld_32, \
                         ld_33, ld_34, ld_35, ld_42, ld_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -id_2[k]
                  + f_0 * ld_32[k];

        t_15[k] = -id_3[k]
                  + f_0 * ld_33[k];

        t_16[k] = -id_4[k]
                  + f_0 * ld_34[k];

        t_17[k] = -id_5[k]
                  + f_0 * ld_35[k];

        t_18[k] = f_0 * ld_42[k];

        t_19[k] = f_0 * ld_43[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, id_6, id_7, ld_44, ld_45, ld_46, \
                         ld_47, ld_48, ld_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * ld_44[k];

        t_21[k] = f_0 * ld_45[k];

        t_22[k] = f_0 * ld_46[k];

        t_23[k] = f_0 * ld_47[k];

        t_24[k] = -id_6[k]
                  + f_0 * ld_48[k];

        t_25[k] = -id_7[k]
                  + f_0 * ld_49[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, id_8, id_9, id_10, id_11, id_12, ld_50, \
                         ld_51, ld_52, ld_53, ld_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -id_8[k]
                  + f_0 * ld_50[k];

        t_27[k] = -id_9[k]
                  + f_0 * ld_51[k];

        t_28[k] = -id_10[k]
                  + f_0 * ld_52[k];

        t_29[k] = -id_11[k]
                  + f_0 * ld_53[k];

        t_30[k] = -2.0 * id_12[k]
                  + f_0 * ld_54[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, id_13, id_14, id_15, id_16, id_17, \
                         ld_55, ld_56, ld_57, ld_58, ld_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -2.0 * id_13[k]
                  + f_0 * ld_55[k];

        t_32[k] = -2.0 * id_14[k]
                  + f_0 * ld_56[k];

        t_33[k] = -2.0 * id_15[k]
                  + f_0 * ld_57[k];

        t_34[k] = -2.0 * id_16[k]
                  + f_0 * ld_58[k];

        t_35[k] = -2.0 * id_17[k]
                  + f_0 * ld_59[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, t_42, id_18, ld_66, ld_67, ld_68, \
                         ld_69, ld_70, ld_71, ld_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * ld_66[k];

        t_37[k] = f_0 * ld_67[k];

        t_38[k] = f_0 * ld_68[k];

        t_39[k] = f_0 * ld_69[k];

        t_40[k] = f_0 * ld_70[k];

        t_41[k] = f_0 * ld_71[k];

        t_42[k] = -id_18[k]
                  + f_0 * ld_72[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, id_19, id_20, id_21, id_22, id_23, \
                         ld_73, ld_74, ld_75, ld_76, ld_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -id_19[k]
                  + f_0 * ld_73[k];

        t_44[k] = -id_20[k]
                  + f_0 * ld_74[k];

        t_45[k] = -id_21[k]
                  + f_0 * ld_75[k];

        t_46[k] = -id_22[k]
                  + f_0 * ld_76[k];

        t_47[k] = -id_23[k]
                  + f_0 * ld_77[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, id_24, id_25, id_26, id_27, id_28, \
                         ld_78, ld_79, ld_80, ld_81, ld_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -2.0 * id_24[k]
                  + f_0 * ld_78[k];

        t_49[k] = -2.0 * id_25[k]
                  + f_0 * ld_79[k];

        t_50[k] = -2.0 * id_26[k]
                  + f_0 * ld_80[k];

        t_51[k] = -2.0 * id_27[k]
                  + f_0 * ld_81[k];

        t_52[k] = -2.0 * id_28[k]
                  + f_0 * ld_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, id_29, id_30, id_31, id_32, id_33, \
                         ld_83, ld_84, ld_85, ld_86, ld_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -2.0 * id_29[k]
                  + f_0 * ld_83[k];

        t_54[k] = -3.0 * id_30[k]
                  + f_0 * ld_84[k];

        t_55[k] = -3.0 * id_31[k]
                  + f_0 * ld_85[k];

        t_56[k] = -3.0 * id_32[k]
                  + f_0 * ld_86[k];

        t_57[k] = -3.0 * id_33[k]
                  + f_0 * ld_87[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, t_63, t_64, id_34, id_35, ld_88, ld_89, \
                         ld_96, ld_97, ld_98, ld_99, ld_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -3.0 * id_34[k]
                  + f_0 * ld_88[k];

        t_59[k] = -3.0 * id_35[k]
                  + f_0 * ld_89[k];

        t_60[k] = f_0 * ld_96[k];

        t_61[k] = f_0 * ld_97[k];

        t_62[k] = f_0 * ld_98[k];

        t_63[k] = f_0 * ld_99[k];

        t_64[k] = f_0 * ld_100[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, id_36, id_37, id_38, id_39, ld_101, \
                         ld_102, ld_103, ld_104, ld_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_0 * ld_101[k];

        t_66[k] = -id_36[k]
                  + f_0 * ld_102[k];

        t_67[k] = -id_37[k]
                  + f_0 * ld_103[k];

        t_68[k] = -id_38[k]
                  + f_0 * ld_104[k];

        t_69[k] = -id_39[k]
                  + f_0 * ld_105[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, id_40, id_41, id_42, id_43, id_44, \
                         ld_106, ld_107, ld_108, ld_109, ld_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -id_40[k]
                  + f_0 * ld_106[k];

        t_71[k] = -id_41[k]
                  + f_0 * ld_107[k];

        t_72[k] = -2.0 * id_42[k]
                  + f_0 * ld_108[k];

        t_73[k] = -2.0 * id_43[k]
                  + f_0 * ld_109[k];

        t_74[k] = -2.0 * id_44[k]
                  + f_0 * ld_110[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, id_45, id_46, id_47, id_48, id_49, \
                         ld_111, ld_112, ld_113, ld_114, ld_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * id_45[k]
                  + f_0 * ld_111[k];

        t_76[k] = -2.0 * id_46[k]
                  + f_0 * ld_112[k];

        t_77[k] = -2.0 * id_47[k]
                  + f_0 * ld_113[k];

        t_78[k] = -3.0 * id_48[k]
                  + f_0 * ld_114[k];

        t_79[k] = -3.0 * id_49[k]
                  + f_0 * ld_115[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, id_50, id_51, id_52, id_53, id_54, \
                         ld_116, ld_117, ld_118, ld_119, ld_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * id_50[k]
                  + f_0 * ld_116[k];

        t_81[k] = -3.0 * id_51[k]
                  + f_0 * ld_117[k];

        t_82[k] = -3.0 * id_52[k]
                  + f_0 * ld_118[k];

        t_83[k] = -3.0 * id_53[k]
                  + f_0 * ld_119[k];

        t_84[k] = -4.0 * id_54[k]
                  + f_0 * ld_120[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, id_55, id_56, id_57, id_58, id_59, \
                         ld_121, ld_122, ld_123, ld_124, ld_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * id_55[k]
                  + f_0 * ld_121[k];

        t_86[k] = -4.0 * id_56[k]
                  + f_0 * ld_122[k];

        t_87[k] = -4.0 * id_57[k]
                  + f_0 * ld_123[k];

        t_88[k] = -4.0 * id_58[k]
                  + f_0 * ld_124[k];

        t_89[k] = -4.0 * id_59[k]
                  + f_0 * ld_125[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, id_60, ld_132, ld_133, \
                         ld_134, ld_135, ld_136, ld_137, ld_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * ld_132[k];

        t_91[k] = f_0 * ld_133[k];

        t_92[k] = f_0 * ld_134[k];

        t_93[k] = f_0 * ld_135[k];

        t_94[k] = f_0 * ld_136[k];

        t_95[k] = f_0 * ld_137[k];

        t_96[k] = -id_60[k]
                  + f_0 * ld_138[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, id_61, id_62, id_63, id_64, id_65, \
                         ld_139, ld_140, ld_141, ld_142, ld_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -id_61[k]
                  + f_0 * ld_139[k];

        t_98[k] = -id_62[k]
                  + f_0 * ld_140[k];

        t_99[k] = -id_63[k]
                  + f_0 * ld_141[k];

        t_100[k] = -id_64[k]
                   + f_0 * ld_142[k];

        t_101[k] = -id_65[k]
                   + f_0 * ld_143[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, id_66, id_67, id_68, id_69, id_70, \
                         ld_144, ld_145, ld_146, ld_147, ld_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -2.0 * id_66[k]
                   + f_0 * ld_144[k];

        t_103[k] = -2.0 * id_67[k]
                   + f_0 * ld_145[k];

        t_104[k] = -2.0 * id_68[k]
                   + f_0 * ld_146[k];

        t_105[k] = -2.0 * id_69[k]
                   + f_0 * ld_147[k];

        t_106[k] = -2.0 * id_70[k]
                   + f_0 * ld_148[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, id_71, id_72, id_73, id_74, id_75, \
                         ld_149, ld_150, ld_151, ld_152, ld_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -2.0 * id_71[k]
                   + f_0 * ld_149[k];

        t_108[k] = -3.0 * id_72[k]
                   + f_0 * ld_150[k];

        t_109[k] = -3.0 * id_73[k]
                   + f_0 * ld_151[k];

        t_110[k] = -3.0 * id_74[k]
                   + f_0 * ld_152[k];

        t_111[k] = -3.0 * id_75[k]
                   + f_0 * ld_153[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, id_76, id_77, id_78, id_79, id_80, \
                         ld_154, ld_155, ld_156, ld_157, ld_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -3.0 * id_76[k]
                   + f_0 * ld_154[k];

        t_113[k] = -3.0 * id_77[k]
                   + f_0 * ld_155[k];

        t_114[k] = -4.0 * id_78[k]
                   + f_0 * ld_156[k];

        t_115[k] = -4.0 * id_79[k]
                   + f_0 * ld_157[k];

        t_116[k] = -4.0 * id_80[k]
                   + f_0 * ld_158[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, id_81, id_82, id_83, id_84, id_85, \
                         ld_159, ld_160, ld_161, ld_162, ld_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -4.0 * id_81[k]
                   + f_0 * ld_159[k];

        t_118[k] = -4.0 * id_82[k]
                   + f_0 * ld_160[k];

        t_119[k] = -4.0 * id_83[k]
                   + f_0 * ld_161[k];

        t_120[k] = -5.0 * id_84[k]
                   + f_0 * ld_162[k];

        t_121[k] = -5.0 * id_85[k]
                   + f_0 * ld_163[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, id_86, id_87, id_88, id_89, \
                         ld_164, ld_165, ld_166, ld_167, ld_174, \
                         ld_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -5.0 * id_86[k]
                   + f_0 * ld_164[k];

        t_123[k] = -5.0 * id_87[k]
                   + f_0 * ld_165[k];

        t_124[k] = -5.0 * id_88[k]
                   + f_0 * ld_166[k];

        t_125[k] = -5.0 * id_89[k]
                   + f_0 * ld_167[k];

        t_126[k] = f_0 * ld_174[k];

        t_127[k] = f_0 * ld_175[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, id_90, id_91, ld_176, \
                         ld_177, ld_178, ld_179, ld_180, ld_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_0 * ld_176[k];

        t_129[k] = f_0 * ld_177[k];

        t_130[k] = f_0 * ld_178[k];

        t_131[k] = f_0 * ld_179[k];

        t_132[k] = -id_90[k]
                   + f_0 * ld_180[k];

        t_133[k] = -id_91[k]
                   + f_0 * ld_181[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, id_92, id_93, id_94, id_95, id_96, \
                         ld_182, ld_183, ld_184, ld_185, ld_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -id_92[k]
                   + f_0 * ld_182[k];

        t_135[k] = -id_93[k]
                   + f_0 * ld_183[k];

        t_136[k] = -id_94[k]
                   + f_0 * ld_184[k];

        t_137[k] = -id_95[k]
                   + f_0 * ld_185[k];

        t_138[k] = -2.0 * id_96[k]
                   + f_0 * ld_186[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, id_97, id_98, id_99, id_100, \
                         id_101, ld_187, ld_188, ld_189, ld_190, \
                         ld_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -2.0 * id_97[k]
                   + f_0 * ld_187[k];

        t_140[k] = -2.0 * id_98[k]
                   + f_0 * ld_188[k];

        t_141[k] = -2.0 * id_99[k]
                   + f_0 * ld_189[k];

        t_142[k] = -2.0 * id_100[k]
                   + f_0 * ld_190[k];

        t_143[k] = -2.0 * id_101[k]
                   + f_0 * ld_191[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, id_102, id_103, id_104, id_105, \
                         id_106, ld_192, ld_193, ld_194, ld_195, \
                         ld_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -3.0 * id_102[k]
                   + f_0 * ld_192[k];

        t_145[k] = -3.0 * id_103[k]
                   + f_0 * ld_193[k];

        t_146[k] = -3.0 * id_104[k]
                   + f_0 * ld_194[k];

        t_147[k] = -3.0 * id_105[k]
                   + f_0 * ld_195[k];

        t_148[k] = -3.0 * id_106[k]
                   + f_0 * ld_196[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, id_107, id_108, id_109, id_110, \
                         id_111, ld_197, ld_198, ld_199, ld_200, \
                         ld_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -3.0 * id_107[k]
                   + f_0 * ld_197[k];

        t_150[k] = -4.0 * id_108[k]
                   + f_0 * ld_198[k];

        t_151[k] = -4.0 * id_109[k]
                   + f_0 * ld_199[k];

        t_152[k] = -4.0 * id_110[k]
                   + f_0 * ld_200[k];

        t_153[k] = -4.0 * id_111[k]
                   + f_0 * ld_201[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, id_112, id_113, id_114, id_115, \
                         id_116, ld_202, ld_203, ld_204, ld_205, \
                         ld_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -4.0 * id_112[k]
                   + f_0 * ld_202[k];

        t_155[k] = -4.0 * id_113[k]
                   + f_0 * ld_203[k];

        t_156[k] = -5.0 * id_114[k]
                   + f_0 * ld_204[k];

        t_157[k] = -5.0 * id_115[k]
                   + f_0 * ld_205[k];

        t_158[k] = -5.0 * id_116[k]
                   + f_0 * ld_206[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, id_117, id_118, id_119, id_120, \
                         id_121, ld_207, ld_208, ld_209, ld_210, \
                         ld_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -5.0 * id_117[k]
                   + f_0 * ld_207[k];

        t_160[k] = -5.0 * id_118[k]
                   + f_0 * ld_208[k];

        t_161[k] = -5.0 * id_119[k]
                   + f_0 * ld_209[k];

        t_162[k] = -6.0 * id_120[k]
                   + f_0 * ld_210[k];

        t_163[k] = -6.0 * id_121[k]
                   + f_0 * ld_211[k];
    }
}

static auto
compute_prim_geom_10_kd_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t id, const size_t ld,
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

    const auto *id_122 = buffer.data(id + 122);
    const auto *id_123 = buffer.data(id + 123);
    const auto *id_124 = buffer.data(id + 124);
    const auto *id_125 = buffer.data(id + 125);
    const auto *id_126 = buffer.data(id + 126);
    const auto *id_127 = buffer.data(id + 127);
    const auto *id_128 = buffer.data(id + 128);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_130 = buffer.data(id + 130);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_132 = buffer.data(id + 132);
    const auto *id_133 = buffer.data(id + 133);
    const auto *id_134 = buffer.data(id + 134);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_138 = buffer.data(id + 138);
    const auto *id_139 = buffer.data(id + 139);
    const auto *id_140 = buffer.data(id + 140);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_144 = buffer.data(id + 144);
    const auto *id_145 = buffer.data(id + 145);
    const auto *id_146 = buffer.data(id + 146);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);
    const auto *id_150 = buffer.data(id + 150);
    const auto *id_151 = buffer.data(id + 151);
    const auto *id_152 = buffer.data(id + 152);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_156 = buffer.data(id + 156);
    const auto *id_157 = buffer.data(id + 157);
    const auto *id_158 = buffer.data(id + 158);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_163 = buffer.data(id + 163);
    const auto *id_164 = buffer.data(id + 164);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_166 = buffer.data(id + 166);
    const auto *id_167 = buffer.data(id + 167);

    const auto *ld_212 = buffer.data(ld + 212);
    const auto *ld_213 = buffer.data(ld + 213);
    const auto *ld_214 = buffer.data(ld + 214);
    const auto *ld_215 = buffer.data(ld + 215);
    const auto *ld_222 = buffer.data(ld + 222);
    const auto *ld_223 = buffer.data(ld + 223);
    const auto *ld_224 = buffer.data(ld + 224);
    const auto *ld_225 = buffer.data(ld + 225);
    const auto *ld_226 = buffer.data(ld + 226);
    const auto *ld_227 = buffer.data(ld + 227);
    const auto *ld_228 = buffer.data(ld + 228);
    const auto *ld_229 = buffer.data(ld + 229);
    const auto *ld_230 = buffer.data(ld + 230);
    const auto *ld_231 = buffer.data(ld + 231);
    const auto *ld_232 = buffer.data(ld + 232);
    const auto *ld_233 = buffer.data(ld + 233);
    const auto *ld_234 = buffer.data(ld + 234);
    const auto *ld_235 = buffer.data(ld + 235);
    const auto *ld_236 = buffer.data(ld + 236);
    const auto *ld_237 = buffer.data(ld + 237);
    const auto *ld_238 = buffer.data(ld + 238);
    const auto *ld_239 = buffer.data(ld + 239);
    const auto *ld_240 = buffer.data(ld + 240);
    const auto *ld_241 = buffer.data(ld + 241);
    const auto *ld_242 = buffer.data(ld + 242);
    const auto *ld_243 = buffer.data(ld + 243);
    const auto *ld_244 = buffer.data(ld + 244);
    const auto *ld_245 = buffer.data(ld + 245);
    const auto *ld_246 = buffer.data(ld + 246);
    const auto *ld_247 = buffer.data(ld + 247);
    const auto *ld_248 = buffer.data(ld + 248);
    const auto *ld_249 = buffer.data(ld + 249);
    const auto *ld_250 = buffer.data(ld + 250);
    const auto *ld_251 = buffer.data(ld + 251);
    const auto *ld_252 = buffer.data(ld + 252);
    const auto *ld_253 = buffer.data(ld + 253);
    const auto *ld_254 = buffer.data(ld + 254);
    const auto *ld_255 = buffer.data(ld + 255);
    const auto *ld_256 = buffer.data(ld + 256);
    const auto *ld_257 = buffer.data(ld + 257);
    const auto *ld_258 = buffer.data(ld + 258);
    const auto *ld_259 = buffer.data(ld + 259);
    const auto *ld_260 = buffer.data(ld + 260);
    const auto *ld_261 = buffer.data(ld + 261);
    const auto *ld_262 = buffer.data(ld + 262);
    const auto *ld_263 = buffer.data(ld + 263);
    const auto *ld_264 = buffer.data(ld + 264);
    const auto *ld_265 = buffer.data(ld + 265);
    const auto *ld_266 = buffer.data(ld + 266);
    const auto *ld_267 = buffer.data(ld + 267);
    const auto *ld_268 = buffer.data(ld + 268);
    const auto *ld_269 = buffer.data(ld + 269);

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, t_169, id_122, id_123, id_124, \
                         id_125, ld_212, ld_213, ld_214, ld_215, ld_222, \
                         ld_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -6.0 * id_122[k]
                   + f_0 * ld_212[k];

        t_165[k] = -6.0 * id_123[k]
                   + f_0 * ld_213[k];

        t_166[k] = -6.0 * id_124[k]
                   + f_0 * ld_214[k];

        t_167[k] = -6.0 * id_125[k]
                   + f_0 * ld_215[k];

        t_168[k] = f_0 * ld_222[k];

        t_169[k] = f_0 * ld_223[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, t_175, id_126, id_127, ld_224, \
                         ld_225, ld_226, ld_227, ld_228, ld_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_0 * ld_224[k];

        t_171[k] = f_0 * ld_225[k];

        t_172[k] = f_0 * ld_226[k];

        t_173[k] = f_0 * ld_227[k];

        t_174[k] = -id_126[k]
                   + f_0 * ld_228[k];

        t_175[k] = -id_127[k]
                   + f_0 * ld_229[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, t_180, id_128, id_129, id_130, id_131, \
                         id_132, ld_230, ld_231, ld_232, ld_233, \
                         ld_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = -id_128[k]
                   + f_0 * ld_230[k];

        t_177[k] = -id_129[k]
                   + f_0 * ld_231[k];

        t_178[k] = -id_130[k]
                   + f_0 * ld_232[k];

        t_179[k] = -id_131[k]
                   + f_0 * ld_233[k];

        t_180[k] = -2.0 * id_132[k]
                   + f_0 * ld_234[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, id_133, id_134, id_135, id_136, \
                         id_137, ld_235, ld_236, ld_237, ld_238, \
                         ld_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = -2.0 * id_133[k]
                   + f_0 * ld_235[k];

        t_182[k] = -2.0 * id_134[k]
                   + f_0 * ld_236[k];

        t_183[k] = -2.0 * id_135[k]
                   + f_0 * ld_237[k];

        t_184[k] = -2.0 * id_136[k]
                   + f_0 * ld_238[k];

        t_185[k] = -2.0 * id_137[k]
                   + f_0 * ld_239[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, id_138, id_139, id_140, id_141, \
                         id_142, ld_240, ld_241, ld_242, ld_243, \
                         ld_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = -3.0 * id_138[k]
                   + f_0 * ld_240[k];

        t_187[k] = -3.0 * id_139[k]
                   + f_0 * ld_241[k];

        t_188[k] = -3.0 * id_140[k]
                   + f_0 * ld_242[k];

        t_189[k] = -3.0 * id_141[k]
                   + f_0 * ld_243[k];

        t_190[k] = -3.0 * id_142[k]
                   + f_0 * ld_244[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, id_143, id_144, id_145, id_146, \
                         id_147, ld_245, ld_246, ld_247, ld_248, \
                         ld_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = -3.0 * id_143[k]
                   + f_0 * ld_245[k];

        t_192[k] = -4.0 * id_144[k]
                   + f_0 * ld_246[k];

        t_193[k] = -4.0 * id_145[k]
                   + f_0 * ld_247[k];

        t_194[k] = -4.0 * id_146[k]
                   + f_0 * ld_248[k];

        t_195[k] = -4.0 * id_147[k]
                   + f_0 * ld_249[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, id_148, id_149, id_150, id_151, \
                         id_152, ld_250, ld_251, ld_252, ld_253, \
                         ld_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -4.0 * id_148[k]
                   + f_0 * ld_250[k];

        t_197[k] = -4.0 * id_149[k]
                   + f_0 * ld_251[k];

        t_198[k] = -5.0 * id_150[k]
                   + f_0 * ld_252[k];

        t_199[k] = -5.0 * id_151[k]
                   + f_0 * ld_253[k];

        t_200[k] = -5.0 * id_152[k]
                   + f_0 * ld_254[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, id_153, id_154, id_155, id_156, \
                         id_157, ld_255, ld_256, ld_257, ld_258, \
                         ld_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = -5.0 * id_153[k]
                   + f_0 * ld_255[k];

        t_202[k] = -5.0 * id_154[k]
                   + f_0 * ld_256[k];

        t_203[k] = -5.0 * id_155[k]
                   + f_0 * ld_257[k];

        t_204[k] = -6.0 * id_156[k]
                   + f_0 * ld_258[k];

        t_205[k] = -6.0 * id_157[k]
                   + f_0 * ld_259[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, id_158, id_159, id_160, id_161, \
                         id_162, ld_260, ld_261, ld_262, ld_263, \
                         ld_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = -6.0 * id_158[k]
                   + f_0 * ld_260[k];

        t_207[k] = -6.0 * id_159[k]
                   + f_0 * ld_261[k];

        t_208[k] = -6.0 * id_160[k]
                   + f_0 * ld_262[k];

        t_209[k] = -6.0 * id_161[k]
                   + f_0 * ld_263[k];

        t_210[k] = -7.0 * id_162[k]
                   + f_0 * ld_264[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, id_163, id_164, id_165, id_166, \
                         id_167, ld_265, ld_266, ld_267, ld_268, \
                         ld_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = -7.0 * id_163[k]
                   + f_0 * ld_265[k];

        t_212[k] = -7.0 * id_164[k]
                   + f_0 * ld_266[k];

        t_213[k] = -7.0 * id_165[k]
                   + f_0 * ld_267[k];

        t_214[k] = -7.0 * id_166[k]
                   + f_0 * ld_268[k];

        t_215[k] = -7.0 * id_167[k]
                   + f_0 * ld_269[k];
    }
}

auto
compute_prim_geom_10_kd_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t id, const size_t ld,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kd_electron_repulsion_2_piece0(buffer, target, id, ld, ncols, alpha);

    compute_prim_geom_10_kd_electron_repulsion_2_piece1(buffer, target, id, ld, ncols, alpha);
}

}  // namespace simdt2ceri
