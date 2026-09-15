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


#include "SimdElectronRepulsionGeom10VrrRecHF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_hf_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gf, const size_t if_,
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

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_67 = buffer.data(gf + 67);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_71 = buffer.data(gf + 71);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_73 = buffer.data(gf + 73);
    const auto *gf_74 = buffer.data(gf + 74);
    const auto *gf_75 = buffer.data(gf + 75);
    const auto *gf_76 = buffer.data(gf + 76);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_81 = buffer.data(gf + 81);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_83 = buffer.data(gf + 83);
    const auto *gf_84 = buffer.data(gf + 84);
    const auto *gf_85 = buffer.data(gf + 85);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_89 = buffer.data(gf + 89);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_91 = buffer.data(gf + 91);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_93 = buffer.data(gf + 93);
    const auto *gf_94 = buffer.data(gf + 94);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_96 = buffer.data(gf + 96);
    const auto *gf_97 = buffer.data(gf + 97);
    const auto *gf_98 = buffer.data(gf + 98);
    const auto *gf_99 = buffer.data(gf + 99);
    const auto *gf_100 = buffer.data(gf + 100);
    const auto *gf_101 = buffer.data(gf + 101);
    const auto *gf_102 = buffer.data(gf + 102);
    const auto *gf_103 = buffer.data(gf + 103);
    const auto *gf_104 = buffer.data(gf + 104);
    const auto *gf_105 = buffer.data(gf + 105);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_107 = buffer.data(gf + 107);
    const auto *gf_108 = buffer.data(gf + 108);
    const auto *gf_109 = buffer.data(gf + 109);
    const auto *gf_110 = buffer.data(gf + 110);
    const auto *gf_111 = buffer.data(gf + 111);
    const auto *gf_112 = buffer.data(gf + 112);
    const auto *gf_113 = buffer.data(gf + 113);
    const auto *gf_114 = buffer.data(gf + 114);
    const auto *gf_115 = buffer.data(gf + 115);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_117 = buffer.data(gf + 117);
    const auto *gf_118 = buffer.data(gf + 118);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_120 = buffer.data(gf + 120);
    const auto *gf_121 = buffer.data(gf + 121);
    const auto *gf_122 = buffer.data(gf + 122);
    const auto *gf_123 = buffer.data(gf + 123);
    const auto *gf_124 = buffer.data(gf + 124);
    const auto *gf_125 = buffer.data(gf + 125);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_127 = buffer.data(gf + 127);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_130 = buffer.data(gf + 130);
    const auto *gf_131 = buffer.data(gf + 131);
    const auto *gf_132 = buffer.data(gf + 132);
    const auto *gf_133 = buffer.data(gf + 133);
    const auto *gf_134 = buffer.data(gf + 134);
    const auto *gf_135 = buffer.data(gf + 135);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_137 = buffer.data(gf + 137);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_140 = buffer.data(gf + 140);
    const auto *gf_141 = buffer.data(gf + 141);
    const auto *gf_142 = buffer.data(gf + 142);
    const auto *gf_143 = buffer.data(gf + 143);
    const auto *gf_144 = buffer.data(gf + 144);
    const auto *gf_145 = buffer.data(gf + 145);
    const auto *gf_146 = buffer.data(gf + 146);
    const auto *gf_147 = buffer.data(gf + 147);
    const auto *gf_148 = buffer.data(gf + 148);
    const auto *gf_149 = buffer.data(gf + 149);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__1 = buffer.data(if_ + 1);
    const auto *if__2 = buffer.data(if_ + 2);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__4 = buffer.data(if_ + 4);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__12 = buffer.data(if_ + 12);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__15 = buffer.data(if_ + 15);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__53 = buffer.data(if_ + 53);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__65 = buffer.data(if_ + 65);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__73 = buffer.data(if_ + 73);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__75 = buffer.data(if_ + 75);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__81 = buffer.data(if_ + 81);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__83 = buffer.data(if_ + 83);
    const auto *if__84 = buffer.data(if_ + 84);
    const auto *if__85 = buffer.data(if_ + 85);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__91 = buffer.data(if_ + 91);
    const auto *if__92 = buffer.data(if_ + 92);
    const auto *if__93 = buffer.data(if_ + 93);
    const auto *if__94 = buffer.data(if_ + 94);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__104 = buffer.data(if_ + 104);
    const auto *if__105 = buffer.data(if_ + 105);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__107 = buffer.data(if_ + 107);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__111 = buffer.data(if_ + 111);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__113 = buffer.data(if_ + 113);
    const auto *if__114 = buffer.data(if_ + 114);
    const auto *if__115 = buffer.data(if_ + 115);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__121 = buffer.data(if_ + 121);
    const auto *if__122 = buffer.data(if_ + 122);
    const auto *if__123 = buffer.data(if_ + 123);
    const auto *if__124 = buffer.data(if_ + 124);
    const auto *if__125 = buffer.data(if_ + 125);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__127 = buffer.data(if_ + 127);
    const auto *if__128 = buffer.data(if_ + 128);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__131 = buffer.data(if_ + 131);
    const auto *if__132 = buffer.data(if_ + 132);
    const auto *if__133 = buffer.data(if_ + 133);
    const auto *if__134 = buffer.data(if_ + 134);
    const auto *if__135 = buffer.data(if_ + 135);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__137 = buffer.data(if_ + 137);
    const auto *if__138 = buffer.data(if_ + 138);
    const auto *if__139 = buffer.data(if_ + 139);
    const auto *if__140 = buffer.data(if_ + 140);
    const auto *if__141 = buffer.data(if_ + 141);
    const auto *if__142 = buffer.data(if_ + 142);
    const auto *if__143 = buffer.data(if_ + 143);
    const auto *if__144 = buffer.data(if_ + 144);
    const auto *if__145 = buffer.data(if_ + 145);
    const auto *if__146 = buffer.data(if_ + 146);
    const auto *if__147 = buffer.data(if_ + 147);
    const auto *if__148 = buffer.data(if_ + 148);
    const auto *if__149 = buffer.data(if_ + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, gf_0, gf_1, gf_2, gf_3, gf_4, if__0, if__1, \
                         if__2, if__3, if__4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -5.0 * gf_0[k]
                 + f_0 * if__0[k];

        t_1[k] = -5.0 * gf_1[k]
                 + f_0 * if__1[k];

        t_2[k] = -5.0 * gf_2[k]
                 + f_0 * if__2[k];

        t_3[k] = -5.0 * gf_3[k]
                 + f_0 * if__3[k];

        t_4[k] = -5.0 * gf_4[k]
                 + f_0 * if__4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, gf_5, gf_6, gf_7, gf_8, gf_9, if__5, if__6, \
                         if__7, if__8, if__9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -5.0 * gf_5[k]
                 + f_0 * if__5[k];

        t_6[k] = -5.0 * gf_6[k]
                 + f_0 * if__6[k];

        t_7[k] = -5.0 * gf_7[k]
                 + f_0 * if__7[k];

        t_8[k] = -5.0 * gf_8[k]
                 + f_0 * if__8[k];

        t_9[k] = -5.0 * gf_9[k]
                 + f_0 * if__9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, gf_10, gf_11, gf_12, gf_13, gf_14, \
                         if__10, if__11, if__12, if__13, if__14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -4.0 * gf_10[k]
                  + f_0 * if__10[k];

        t_11[k] = -4.0 * gf_11[k]
                  + f_0 * if__11[k];

        t_12[k] = -4.0 * gf_12[k]
                  + f_0 * if__12[k];

        t_13[k] = -4.0 * gf_13[k]
                  + f_0 * if__13[k];

        t_14[k] = -4.0 * gf_14[k]
                  + f_0 * if__14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, gf_15, gf_16, gf_17, gf_18, gf_19, \
                         if__15, if__16, if__17, if__18, if__19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -4.0 * gf_15[k]
                  + f_0 * if__15[k];

        t_16[k] = -4.0 * gf_16[k]
                  + f_0 * if__16[k];

        t_17[k] = -4.0 * gf_17[k]
                  + f_0 * if__17[k];

        t_18[k] = -4.0 * gf_18[k]
                  + f_0 * if__18[k];

        t_19[k] = -4.0 * gf_19[k]
                  + f_0 * if__19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, gf_20, gf_21, gf_22, gf_23, gf_24, \
                         if__20, if__21, if__22, if__23, if__24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -4.0 * gf_20[k]
                  + f_0 * if__20[k];

        t_21[k] = -4.0 * gf_21[k]
                  + f_0 * if__21[k];

        t_22[k] = -4.0 * gf_22[k]
                  + f_0 * if__22[k];

        t_23[k] = -4.0 * gf_23[k]
                  + f_0 * if__23[k];

        t_24[k] = -4.0 * gf_24[k]
                  + f_0 * if__24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, gf_25, gf_26, gf_27, gf_28, gf_29, \
                         if__25, if__26, if__27, if__28, if__29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -4.0 * gf_25[k]
                  + f_0 * if__25[k];

        t_26[k] = -4.0 * gf_26[k]
                  + f_0 * if__26[k];

        t_27[k] = -4.0 * gf_27[k]
                  + f_0 * if__27[k];

        t_28[k] = -4.0 * gf_28[k]
                  + f_0 * if__28[k];

        t_29[k] = -4.0 * gf_29[k]
                  + f_0 * if__29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, gf_30, gf_31, gf_32, gf_33, gf_34, \
                         if__30, if__31, if__32, if__33, if__34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -3.0 * gf_30[k]
                  + f_0 * if__30[k];

        t_31[k] = -3.0 * gf_31[k]
                  + f_0 * if__31[k];

        t_32[k] = -3.0 * gf_32[k]
                  + f_0 * if__32[k];

        t_33[k] = -3.0 * gf_33[k]
                  + f_0 * if__33[k];

        t_34[k] = -3.0 * gf_34[k]
                  + f_0 * if__34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, gf_35, gf_36, gf_37, gf_38, gf_39, \
                         if__35, if__36, if__37, if__38, if__39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -3.0 * gf_35[k]
                  + f_0 * if__35[k];

        t_36[k] = -3.0 * gf_36[k]
                  + f_0 * if__36[k];

        t_37[k] = -3.0 * gf_37[k]
                  + f_0 * if__37[k];

        t_38[k] = -3.0 * gf_38[k]
                  + f_0 * if__38[k];

        t_39[k] = -3.0 * gf_39[k]
                  + f_0 * if__39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, gf_40, gf_41, gf_42, gf_43, gf_44, \
                         if__40, if__41, if__42, if__43, if__44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * gf_40[k]
                  + f_0 * if__40[k];

        t_41[k] = -3.0 * gf_41[k]
                  + f_0 * if__41[k];

        t_42[k] = -3.0 * gf_42[k]
                  + f_0 * if__42[k];

        t_43[k] = -3.0 * gf_43[k]
                  + f_0 * if__43[k];

        t_44[k] = -3.0 * gf_44[k]
                  + f_0 * if__44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, gf_45, gf_46, gf_47, gf_48, gf_49, \
                         if__45, if__46, if__47, if__48, if__49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -3.0 * gf_45[k]
                  + f_0 * if__45[k];

        t_46[k] = -3.0 * gf_46[k]
                  + f_0 * if__46[k];

        t_47[k] = -3.0 * gf_47[k]
                  + f_0 * if__47[k];

        t_48[k] = -3.0 * gf_48[k]
                  + f_0 * if__48[k];

        t_49[k] = -3.0 * gf_49[k]
                  + f_0 * if__49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, gf_50, gf_51, gf_52, gf_53, gf_54, \
                         if__50, if__51, if__52, if__53, if__54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -3.0 * gf_50[k]
                  + f_0 * if__50[k];

        t_51[k] = -3.0 * gf_51[k]
                  + f_0 * if__51[k];

        t_52[k] = -3.0 * gf_52[k]
                  + f_0 * if__52[k];

        t_53[k] = -3.0 * gf_53[k]
                  + f_0 * if__53[k];

        t_54[k] = -3.0 * gf_54[k]
                  + f_0 * if__54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, gf_55, gf_56, gf_57, gf_58, gf_59, \
                         if__55, if__56, if__57, if__58, if__59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -3.0 * gf_55[k]
                  + f_0 * if__55[k];

        t_56[k] = -3.0 * gf_56[k]
                  + f_0 * if__56[k];

        t_57[k] = -3.0 * gf_57[k]
                  + f_0 * if__57[k];

        t_58[k] = -3.0 * gf_58[k]
                  + f_0 * if__58[k];

        t_59[k] = -3.0 * gf_59[k]
                  + f_0 * if__59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, gf_60, gf_61, gf_62, gf_63, gf_64, \
                         if__60, if__61, if__62, if__63, if__64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -2.0 * gf_60[k]
                  + f_0 * if__60[k];

        t_61[k] = -2.0 * gf_61[k]
                  + f_0 * if__61[k];

        t_62[k] = -2.0 * gf_62[k]
                  + f_0 * if__62[k];

        t_63[k] = -2.0 * gf_63[k]
                  + f_0 * if__63[k];

        t_64[k] = -2.0 * gf_64[k]
                  + f_0 * if__64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, gf_65, gf_66, gf_67, gf_68, gf_69, \
                         if__65, if__66, if__67, if__68, if__69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -2.0 * gf_65[k]
                  + f_0 * if__65[k];

        t_66[k] = -2.0 * gf_66[k]
                  + f_0 * if__66[k];

        t_67[k] = -2.0 * gf_67[k]
                  + f_0 * if__67[k];

        t_68[k] = -2.0 * gf_68[k]
                  + f_0 * if__68[k];

        t_69[k] = -2.0 * gf_69[k]
                  + f_0 * if__69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, gf_70, gf_71, gf_72, gf_73, gf_74, \
                         if__70, if__71, if__72, if__73, if__74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * gf_70[k]
                  + f_0 * if__70[k];

        t_71[k] = -2.0 * gf_71[k]
                  + f_0 * if__71[k];

        t_72[k] = -2.0 * gf_72[k]
                  + f_0 * if__72[k];

        t_73[k] = -2.0 * gf_73[k]
                  + f_0 * if__73[k];

        t_74[k] = -2.0 * gf_74[k]
                  + f_0 * if__74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, gf_75, gf_76, gf_77, gf_78, gf_79, \
                         if__75, if__76, if__77, if__78, if__79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * gf_75[k]
                  + f_0 * if__75[k];

        t_76[k] = -2.0 * gf_76[k]
                  + f_0 * if__76[k];

        t_77[k] = -2.0 * gf_77[k]
                  + f_0 * if__77[k];

        t_78[k] = -2.0 * gf_78[k]
                  + f_0 * if__78[k];

        t_79[k] = -2.0 * gf_79[k]
                  + f_0 * if__79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, gf_80, gf_81, gf_82, gf_83, gf_84, \
                         if__80, if__81, if__82, if__83, if__84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -2.0 * gf_80[k]
                  + f_0 * if__80[k];

        t_81[k] = -2.0 * gf_81[k]
                  + f_0 * if__81[k];

        t_82[k] = -2.0 * gf_82[k]
                  + f_0 * if__82[k];

        t_83[k] = -2.0 * gf_83[k]
                  + f_0 * if__83[k];

        t_84[k] = -2.0 * gf_84[k]
                  + f_0 * if__84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, gf_85, gf_86, gf_87, gf_88, gf_89, \
                         if__85, if__86, if__87, if__88, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * gf_85[k]
                  + f_0 * if__85[k];

        t_86[k] = -2.0 * gf_86[k]
                  + f_0 * if__86[k];

        t_87[k] = -2.0 * gf_87[k]
                  + f_0 * if__87[k];

        t_88[k] = -2.0 * gf_88[k]
                  + f_0 * if__88[k];

        t_89[k] = -2.0 * gf_89[k]
                  + f_0 * if__89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, gf_90, gf_91, gf_92, gf_93, gf_94, \
                         if__90, if__91, if__92, if__93, if__94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * gf_90[k]
                  + f_0 * if__90[k];

        t_91[k] = -2.0 * gf_91[k]
                  + f_0 * if__91[k];

        t_92[k] = -2.0 * gf_92[k]
                  + f_0 * if__92[k];

        t_93[k] = -2.0 * gf_93[k]
                  + f_0 * if__93[k];

        t_94[k] = -2.0 * gf_94[k]
                  + f_0 * if__94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, gf_95, gf_96, gf_97, gf_98, gf_99, \
                         if__95, if__96, if__97, if__98, if__99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * gf_95[k]
                  + f_0 * if__95[k];

        t_96[k] = -2.0 * gf_96[k]
                  + f_0 * if__96[k];

        t_97[k] = -2.0 * gf_97[k]
                  + f_0 * if__97[k];

        t_98[k] = -2.0 * gf_98[k]
                  + f_0 * if__98[k];

        t_99[k] = -2.0 * gf_99[k]
                  + f_0 * if__99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, gf_100, gf_101, gf_102, gf_103, \
                         gf_104, if__100, if__101, if__102, if__103, \
                         if__104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -gf_100[k]
                   + f_0 * if__100[k];

        t_101[k] = -gf_101[k]
                   + f_0 * if__101[k];

        t_102[k] = -gf_102[k]
                   + f_0 * if__102[k];

        t_103[k] = -gf_103[k]
                   + f_0 * if__103[k];

        t_104[k] = -gf_104[k]
                   + f_0 * if__104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, gf_105, gf_106, gf_107, gf_108, \
                         gf_109, if__105, if__106, if__107, if__108, \
                         if__109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -gf_105[k]
                   + f_0 * if__105[k];

        t_106[k] = -gf_106[k]
                   + f_0 * if__106[k];

        t_107[k] = -gf_107[k]
                   + f_0 * if__107[k];

        t_108[k] = -gf_108[k]
                   + f_0 * if__108[k];

        t_109[k] = -gf_109[k]
                   + f_0 * if__109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, gf_110, gf_111, gf_112, gf_113, \
                         gf_114, if__110, if__111, if__112, if__113, \
                         if__114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -gf_110[k]
                   + f_0 * if__110[k];

        t_111[k] = -gf_111[k]
                   + f_0 * if__111[k];

        t_112[k] = -gf_112[k]
                   + f_0 * if__112[k];

        t_113[k] = -gf_113[k]
                   + f_0 * if__113[k];

        t_114[k] = -gf_114[k]
                   + f_0 * if__114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, gf_115, gf_116, gf_117, gf_118, \
                         gf_119, if__115, if__116, if__117, if__118, \
                         if__119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -gf_115[k]
                   + f_0 * if__115[k];

        t_116[k] = -gf_116[k]
                   + f_0 * if__116[k];

        t_117[k] = -gf_117[k]
                   + f_0 * if__117[k];

        t_118[k] = -gf_118[k]
                   + f_0 * if__118[k];

        t_119[k] = -gf_119[k]
                   + f_0 * if__119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, gf_120, gf_121, gf_122, gf_123, \
                         gf_124, if__120, if__121, if__122, if__123, \
                         if__124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -gf_120[k]
                   + f_0 * if__120[k];

        t_121[k] = -gf_121[k]
                   + f_0 * if__121[k];

        t_122[k] = -gf_122[k]
                   + f_0 * if__122[k];

        t_123[k] = -gf_123[k]
                   + f_0 * if__123[k];

        t_124[k] = -gf_124[k]
                   + f_0 * if__124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, gf_125, gf_126, gf_127, gf_128, \
                         gf_129, if__125, if__126, if__127, if__128, \
                         if__129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -gf_125[k]
                   + f_0 * if__125[k];

        t_126[k] = -gf_126[k]
                   + f_0 * if__126[k];

        t_127[k] = -gf_127[k]
                   + f_0 * if__127[k];

        t_128[k] = -gf_128[k]
                   + f_0 * if__128[k];

        t_129[k] = -gf_129[k]
                   + f_0 * if__129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, gf_130, gf_131, gf_132, gf_133, \
                         gf_134, if__130, if__131, if__132, if__133, \
                         if__134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -gf_130[k]
                   + f_0 * if__130[k];

        t_131[k] = -gf_131[k]
                   + f_0 * if__131[k];

        t_132[k] = -gf_132[k]
                   + f_0 * if__132[k];

        t_133[k] = -gf_133[k]
                   + f_0 * if__133[k];

        t_134[k] = -gf_134[k]
                   + f_0 * if__134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, gf_135, gf_136, gf_137, gf_138, \
                         gf_139, if__135, if__136, if__137, if__138, \
                         if__139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -gf_135[k]
                   + f_0 * if__135[k];

        t_136[k] = -gf_136[k]
                   + f_0 * if__136[k];

        t_137[k] = -gf_137[k]
                   + f_0 * if__137[k];

        t_138[k] = -gf_138[k]
                   + f_0 * if__138[k];

        t_139[k] = -gf_139[k]
                   + f_0 * if__139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, gf_140, gf_141, gf_142, gf_143, \
                         gf_144, if__140, if__141, if__142, if__143, \
                         if__144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -gf_140[k]
                   + f_0 * if__140[k];

        t_141[k] = -gf_141[k]
                   + f_0 * if__141[k];

        t_142[k] = -gf_142[k]
                   + f_0 * if__142[k];

        t_143[k] = -gf_143[k]
                   + f_0 * if__143[k];

        t_144[k] = -gf_144[k]
                   + f_0 * if__144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, gf_145, gf_146, gf_147, gf_148, \
                         gf_149, if__145, if__146, if__147, if__148, \
                         if__149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -gf_145[k]
                   + f_0 * if__145[k];

        t_146[k] = -gf_146[k]
                   + f_0 * if__146[k];

        t_147[k] = -gf_147[k]
                   + f_0 * if__147[k];

        t_148[k] = -gf_148[k]
                   + f_0 * if__148[k];

        t_149[k] = -gf_149[k]
                   + f_0 * if__149[k];
    }
}

static auto
compute_prim_geom_10_hf_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t if_, const size_t ncols,
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

    const auto *if__150 = buffer.data(if_ + 150);
    const auto *if__151 = buffer.data(if_ + 151);
    const auto *if__152 = buffer.data(if_ + 152);
    const auto *if__153 = buffer.data(if_ + 153);
    const auto *if__154 = buffer.data(if_ + 154);
    const auto *if__155 = buffer.data(if_ + 155);
    const auto *if__156 = buffer.data(if_ + 156);
    const auto *if__157 = buffer.data(if_ + 157);
    const auto *if__158 = buffer.data(if_ + 158);
    const auto *if__159 = buffer.data(if_ + 159);
    const auto *if__160 = buffer.data(if_ + 160);
    const auto *if__161 = buffer.data(if_ + 161);
    const auto *if__162 = buffer.data(if_ + 162);
    const auto *if__163 = buffer.data(if_ + 163);
    const auto *if__164 = buffer.data(if_ + 164);
    const auto *if__165 = buffer.data(if_ + 165);
    const auto *if__166 = buffer.data(if_ + 166);
    const auto *if__167 = buffer.data(if_ + 167);
    const auto *if__168 = buffer.data(if_ + 168);
    const auto *if__169 = buffer.data(if_ + 169);
    const auto *if__170 = buffer.data(if_ + 170);
    const auto *if__171 = buffer.data(if_ + 171);
    const auto *if__172 = buffer.data(if_ + 172);
    const auto *if__173 = buffer.data(if_ + 173);
    const auto *if__174 = buffer.data(if_ + 174);
    const auto *if__175 = buffer.data(if_ + 175);
    const auto *if__176 = buffer.data(if_ + 176);
    const auto *if__177 = buffer.data(if_ + 177);
    const auto *if__178 = buffer.data(if_ + 178);
    const auto *if__179 = buffer.data(if_ + 179);
    const auto *if__180 = buffer.data(if_ + 180);
    const auto *if__181 = buffer.data(if_ + 181);
    const auto *if__182 = buffer.data(if_ + 182);
    const auto *if__183 = buffer.data(if_ + 183);
    const auto *if__184 = buffer.data(if_ + 184);
    const auto *if__185 = buffer.data(if_ + 185);
    const auto *if__186 = buffer.data(if_ + 186);
    const auto *if__187 = buffer.data(if_ + 187);
    const auto *if__188 = buffer.data(if_ + 188);
    const auto *if__189 = buffer.data(if_ + 189);
    const auto *if__190 = buffer.data(if_ + 190);
    const auto *if__191 = buffer.data(if_ + 191);
    const auto *if__192 = buffer.data(if_ + 192);
    const auto *if__193 = buffer.data(if_ + 193);
    const auto *if__194 = buffer.data(if_ + 194);
    const auto *if__195 = buffer.data(if_ + 195);
    const auto *if__196 = buffer.data(if_ + 196);
    const auto *if__197 = buffer.data(if_ + 197);
    const auto *if__198 = buffer.data(if_ + 198);
    const auto *if__199 = buffer.data(if_ + 199);
    const auto *if__200 = buffer.data(if_ + 200);
    const auto *if__201 = buffer.data(if_ + 201);
    const auto *if__202 = buffer.data(if_ + 202);
    const auto *if__203 = buffer.data(if_ + 203);
    const auto *if__204 = buffer.data(if_ + 204);
    const auto *if__205 = buffer.data(if_ + 205);
    const auto *if__206 = buffer.data(if_ + 206);
    const auto *if__207 = buffer.data(if_ + 207);
    const auto *if__208 = buffer.data(if_ + 208);
    const auto *if__209 = buffer.data(if_ + 209);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, t_156, t_157, if__150, \
                         if__151, if__152, if__153, if__154, if__155, if__156, \
                         if__157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_0 * if__150[k];

        t_151[k] = f_0 * if__151[k];

        t_152[k] = f_0 * if__152[k];

        t_153[k] = f_0 * if__153[k];

        t_154[k] = f_0 * if__154[k];

        t_155[k] = f_0 * if__155[k];

        t_156[k] = f_0 * if__156[k];

        t_157[k] = f_0 * if__157[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, t_163, t_164, t_165, if__158, \
                         if__159, if__160, if__161, if__162, if__163, if__164, \
                         if__165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_0 * if__158[k];

        t_159[k] = f_0 * if__159[k];

        t_160[k] = f_0 * if__160[k];

        t_161[k] = f_0 * if__161[k];

        t_162[k] = f_0 * if__162[k];

        t_163[k] = f_0 * if__163[k];

        t_164[k] = f_0 * if__164[k];

        t_165[k] = f_0 * if__165[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, t_171, t_172, t_173, if__166, \
                         if__167, if__168, if__169, if__170, if__171, if__172, \
                         if__173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_0 * if__166[k];

        t_167[k] = f_0 * if__167[k];

        t_168[k] = f_0 * if__168[k];

        t_169[k] = f_0 * if__169[k];

        t_170[k] = f_0 * if__170[k];

        t_171[k] = f_0 * if__171[k];

        t_172[k] = f_0 * if__172[k];

        t_173[k] = f_0 * if__173[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, t_179, t_180, t_181, if__174, \
                         if__175, if__176, if__177, if__178, if__179, if__180, \
                         if__181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_0 * if__174[k];

        t_175[k] = f_0 * if__175[k];

        t_176[k] = f_0 * if__176[k];

        t_177[k] = f_0 * if__177[k];

        t_178[k] = f_0 * if__178[k];

        t_179[k] = f_0 * if__179[k];

        t_180[k] = f_0 * if__180[k];

        t_181[k] = f_0 * if__181[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, t_187, t_188, t_189, if__182, \
                         if__183, if__184, if__185, if__186, if__187, if__188, \
                         if__189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_0 * if__182[k];

        t_183[k] = f_0 * if__183[k];

        t_184[k] = f_0 * if__184[k];

        t_185[k] = f_0 * if__185[k];

        t_186[k] = f_0 * if__186[k];

        t_187[k] = f_0 * if__187[k];

        t_188[k] = f_0 * if__188[k];

        t_189[k] = f_0 * if__189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, t_195, t_196, t_197, if__190, \
                         if__191, if__192, if__193, if__194, if__195, if__196, \
                         if__197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_0 * if__190[k];

        t_191[k] = f_0 * if__191[k];

        t_192[k] = f_0 * if__192[k];

        t_193[k] = f_0 * if__193[k];

        t_194[k] = f_0 * if__194[k];

        t_195[k] = f_0 * if__195[k];

        t_196[k] = f_0 * if__196[k];

        t_197[k] = f_0 * if__197[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, t_203, t_204, t_205, if__198, \
                         if__199, if__200, if__201, if__202, if__203, if__204, \
                         if__205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_0 * if__198[k];

        t_199[k] = f_0 * if__199[k];

        t_200[k] = f_0 * if__200[k];

        t_201[k] = f_0 * if__201[k];

        t_202[k] = f_0 * if__202[k];

        t_203[k] = f_0 * if__203[k];

        t_204[k] = f_0 * if__204[k];

        t_205[k] = f_0 * if__205[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, if__206, if__207, if__208, \
                         if__209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_0 * if__206[k];

        t_207[k] = f_0 * if__207[k];

        t_208[k] = f_0 * if__208[k];

        t_209[k] = f_0 * if__209[k];
    }
}

auto
compute_prim_geom_10_hf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t gf, const size_t if_,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hf_electron_repulsion_0_piece0(buffer, target, gf, if_, ncols, alpha);

    compute_prim_geom_10_hf_electron_repulsion_0_piece1(buffer, target, if_, ncols, alpha);
}

static auto
compute_prim_geom_10_hf_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gf, const size_t if_,
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

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_67 = buffer.data(gf + 67);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_71 = buffer.data(gf + 71);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_73 = buffer.data(gf + 73);
    const auto *gf_74 = buffer.data(gf + 74);
    const auto *gf_75 = buffer.data(gf + 75);
    const auto *gf_76 = buffer.data(gf + 76);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_81 = buffer.data(gf + 81);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_83 = buffer.data(gf + 83);
    const auto *gf_84 = buffer.data(gf + 84);
    const auto *gf_85 = buffer.data(gf + 85);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_89 = buffer.data(gf + 89);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_91 = buffer.data(gf + 91);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_93 = buffer.data(gf + 93);
    const auto *gf_94 = buffer.data(gf + 94);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_96 = buffer.data(gf + 96);
    const auto *gf_97 = buffer.data(gf + 97);
    const auto *gf_98 = buffer.data(gf + 98);
    const auto *gf_99 = buffer.data(gf + 99);
    const auto *gf_100 = buffer.data(gf + 100);
    const auto *gf_101 = buffer.data(gf + 101);
    const auto *gf_102 = buffer.data(gf + 102);
    const auto *gf_103 = buffer.data(gf + 103);
    const auto *gf_104 = buffer.data(gf + 104);
    const auto *gf_105 = buffer.data(gf + 105);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_107 = buffer.data(gf + 107);
    const auto *gf_108 = buffer.data(gf + 108);
    const auto *gf_109 = buffer.data(gf + 109);
    const auto *gf_110 = buffer.data(gf + 110);
    const auto *gf_111 = buffer.data(gf + 111);
    const auto *gf_112 = buffer.data(gf + 112);
    const auto *gf_113 = buffer.data(gf + 113);
    const auto *gf_114 = buffer.data(gf + 114);
    const auto *gf_115 = buffer.data(gf + 115);
    const auto *gf_116 = buffer.data(gf + 116);

    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__12 = buffer.data(if_ + 12);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__15 = buffer.data(if_ + 15);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__65 = buffer.data(if_ + 65);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__73 = buffer.data(if_ + 73);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__75 = buffer.data(if_ + 75);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__81 = buffer.data(if_ + 81);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__83 = buffer.data(if_ + 83);
    const auto *if__84 = buffer.data(if_ + 84);
    const auto *if__85 = buffer.data(if_ + 85);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__104 = buffer.data(if_ + 104);
    const auto *if__105 = buffer.data(if_ + 105);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__107 = buffer.data(if_ + 107);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__111 = buffer.data(if_ + 111);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__113 = buffer.data(if_ + 113);
    const auto *if__114 = buffer.data(if_ + 114);
    const auto *if__115 = buffer.data(if_ + 115);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__121 = buffer.data(if_ + 121);
    const auto *if__122 = buffer.data(if_ + 122);
    const auto *if__123 = buffer.data(if_ + 123);
    const auto *if__124 = buffer.data(if_ + 124);
    const auto *if__125 = buffer.data(if_ + 125);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__127 = buffer.data(if_ + 127);
    const auto *if__128 = buffer.data(if_ + 128);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__131 = buffer.data(if_ + 131);
    const auto *if__132 = buffer.data(if_ + 132);
    const auto *if__133 = buffer.data(if_ + 133);
    const auto *if__134 = buffer.data(if_ + 134);
    const auto *if__135 = buffer.data(if_ + 135);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__137 = buffer.data(if_ + 137);
    const auto *if__138 = buffer.data(if_ + 138);
    const auto *if__139 = buffer.data(if_ + 139);
    const auto *if__150 = buffer.data(if_ + 150);
    const auto *if__151 = buffer.data(if_ + 151);
    const auto *if__152 = buffer.data(if_ + 152);
    const auto *if__153 = buffer.data(if_ + 153);
    const auto *if__154 = buffer.data(if_ + 154);
    const auto *if__155 = buffer.data(if_ + 155);
    const auto *if__156 = buffer.data(if_ + 156);
    const auto *if__157 = buffer.data(if_ + 157);
    const auto *if__158 = buffer.data(if_ + 158);
    const auto *if__159 = buffer.data(if_ + 159);
    const auto *if__160 = buffer.data(if_ + 160);
    const auto *if__161 = buffer.data(if_ + 161);
    const auto *if__162 = buffer.data(if_ + 162);
    const auto *if__163 = buffer.data(if_ + 163);
    const auto *if__164 = buffer.data(if_ + 164);
    const auto *if__165 = buffer.data(if_ + 165);
    const auto *if__166 = buffer.data(if_ + 166);
    const auto *if__167 = buffer.data(if_ + 167);
    const auto *if__168 = buffer.data(if_ + 168);
    const auto *if__169 = buffer.data(if_ + 169);
    const auto *if__170 = buffer.data(if_ + 170);
    const auto *if__171 = buffer.data(if_ + 171);
    const auto *if__172 = buffer.data(if_ + 172);
    const auto *if__173 = buffer.data(if_ + 173);
    const auto *if__174 = buffer.data(if_ + 174);
    const auto *if__175 = buffer.data(if_ + 175);
    const auto *if__176 = buffer.data(if_ + 176);
    const auto *if__177 = buffer.data(if_ + 177);
    const auto *if__178 = buffer.data(if_ + 178);
    const auto *if__179 = buffer.data(if_ + 179);
    const auto *if__180 = buffer.data(if_ + 180);
    const auto *if__181 = buffer.data(if_ + 181);
    const auto *if__182 = buffer.data(if_ + 182);
    const auto *if__183 = buffer.data(if_ + 183);
    const auto *if__184 = buffer.data(if_ + 184);
    const auto *if__185 = buffer.data(if_ + 185);
    const auto *if__186 = buffer.data(if_ + 186);
    const auto *if__187 = buffer.data(if_ + 187);
    const auto *if__188 = buffer.data(if_ + 188);
    const auto *if__189 = buffer.data(if_ + 189);
    const auto *if__190 = buffer.data(if_ + 190);
    const auto *if__191 = buffer.data(if_ + 191);
    const auto *if__192 = buffer.data(if_ + 192);
    const auto *if__193 = buffer.data(if_ + 193);
    const auto *if__194 = buffer.data(if_ + 194);
    const auto *if__195 = buffer.data(if_ + 195);
    const auto *if__196 = buffer.data(if_ + 196);
    const auto *if__197 = buffer.data(if_ + 197);
    const auto *if__198 = buffer.data(if_ + 198);
    const auto *if__199 = buffer.data(if_ + 199);
    const auto *if__210 = buffer.data(if_ + 210);
    const auto *if__211 = buffer.data(if_ + 211);
    const auto *if__212 = buffer.data(if_ + 212);
    const auto *if__213 = buffer.data(if_ + 213);
    const auto *if__214 = buffer.data(if_ + 214);
    const auto *if__215 = buffer.data(if_ + 215);
    const auto *if__216 = buffer.data(if_ + 216);
    const auto *if__217 = buffer.data(if_ + 217);
    const auto *if__218 = buffer.data(if_ + 218);
    const auto *if__219 = buffer.data(if_ + 219);
    const auto *if__220 = buffer.data(if_ + 220);
    const auto *if__221 = buffer.data(if_ + 221);
    const auto *if__222 = buffer.data(if_ + 222);
    const auto *if__223 = buffer.data(if_ + 223);
    const auto *if__224 = buffer.data(if_ + 224);
    const auto *if__225 = buffer.data(if_ + 225);
    const auto *if__226 = buffer.data(if_ + 226);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, if__10, if__11, if__12, \
                         if__13, if__14, if__15, if__16, if__17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * if__10[k];

        t_1[k] = f_0 * if__11[k];

        t_2[k] = f_0 * if__12[k];

        t_3[k] = f_0 * if__13[k];

        t_4[k] = f_0 * if__14[k];

        t_5[k] = f_0 * if__15[k];

        t_6[k] = f_0 * if__16[k];

        t_7[k] = f_0 * if__17[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, gf_0, gf_1, gf_2, gf_3, if__18, \
                         if__19, if__30, if__31, if__32, if__33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * if__18[k];

        t_9[k] = f_0 * if__19[k];

        t_10[k] = -gf_0[k]
                  + f_0 * if__30[k];

        t_11[k] = -gf_1[k]
                  + f_0 * if__31[k];

        t_12[k] = -gf_2[k]
                  + f_0 * if__32[k];

        t_13[k] = -gf_3[k]
                  + f_0 * if__33[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, gf_4, gf_5, gf_6, gf_7, gf_8, if__34, \
                         if__35, if__36, if__37, if__38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -gf_4[k]
                  + f_0 * if__34[k];

        t_15[k] = -gf_5[k]
                  + f_0 * if__35[k];

        t_16[k] = -gf_6[k]
                  + f_0 * if__36[k];

        t_17[k] = -gf_7[k]
                  + f_0 * if__37[k];

        t_18[k] = -gf_8[k]
                  + f_0 * if__38[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, t_25, gf_9, if__39, if__40, \
                         if__41, if__42, if__43, if__44, if__45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -gf_9[k]
                  + f_0 * if__39[k];

        t_20[k] = f_0 * if__40[k];

        t_21[k] = f_0 * if__41[k];

        t_22[k] = f_0 * if__42[k];

        t_23[k] = f_0 * if__43[k];

        t_24[k] = f_0 * if__44[k];

        t_25[k] = f_0 * if__45[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, gf_10, gf_11, if__46, if__47, \
                         if__48, if__49, if__60, if__61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * if__46[k];

        t_27[k] = f_0 * if__47[k];

        t_28[k] = f_0 * if__48[k];

        t_29[k] = f_0 * if__49[k];

        t_30[k] = -2.0 * gf_10[k]
                  + f_0 * if__60[k];

        t_31[k] = -2.0 * gf_11[k]
                  + f_0 * if__61[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, gf_12, gf_13, gf_14, gf_15, gf_16, \
                         if__62, if__63, if__64, if__65, if__66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -2.0 * gf_12[k]
                  + f_0 * if__62[k];

        t_33[k] = -2.0 * gf_13[k]
                  + f_0 * if__63[k];

        t_34[k] = -2.0 * gf_14[k]
                  + f_0 * if__64[k];

        t_35[k] = -2.0 * gf_15[k]
                  + f_0 * if__65[k];

        t_36[k] = -2.0 * gf_16[k]
                  + f_0 * if__66[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, gf_17, gf_18, gf_19, gf_20, gf_21, \
                         if__67, if__68, if__69, if__70, if__71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -2.0 * gf_17[k]
                  + f_0 * if__67[k];

        t_38[k] = -2.0 * gf_18[k]
                  + f_0 * if__68[k];

        t_39[k] = -2.0 * gf_19[k]
                  + f_0 * if__69[k];

        t_40[k] = -gf_20[k]
                  + f_0 * if__70[k];

        t_41[k] = -gf_21[k]
                  + f_0 * if__71[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, gf_22, gf_23, gf_24, gf_25, gf_26, \
                         if__72, if__73, if__74, if__75, if__76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -gf_22[k]
                  + f_0 * if__72[k];

        t_43[k] = -gf_23[k]
                  + f_0 * if__73[k];

        t_44[k] = -gf_24[k]
                  + f_0 * if__74[k];

        t_45[k] = -gf_25[k]
                  + f_0 * if__75[k];

        t_46[k] = -gf_26[k]
                  + f_0 * if__76[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, gf_27, gf_28, gf_29, if__77, \
                         if__78, if__79, if__80, if__81, if__82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -gf_27[k]
                  + f_0 * if__77[k];

        t_48[k] = -gf_28[k]
                  + f_0 * if__78[k];

        t_49[k] = -gf_29[k]
                  + f_0 * if__79[k];

        t_50[k] = f_0 * if__80[k];

        t_51[k] = f_0 * if__81[k];

        t_52[k] = f_0 * if__82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, t_59, if__83, if__84, if__85, \
                         if__86, if__87, if__88, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * if__83[k];

        t_54[k] = f_0 * if__84[k];

        t_55[k] = f_0 * if__85[k];

        t_56[k] = f_0 * if__86[k];

        t_57[k] = f_0 * if__87[k];

        t_58[k] = f_0 * if__88[k];

        t_59[k] = f_0 * if__89[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, gf_30, gf_31, gf_32, gf_33, gf_34, \
                         if__100, if__101, if__102, if__103, if__104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * gf_30[k]
                  + f_0 * if__100[k];

        t_61[k] = -3.0 * gf_31[k]
                  + f_0 * if__101[k];

        t_62[k] = -3.0 * gf_32[k]
                  + f_0 * if__102[k];

        t_63[k] = -3.0 * gf_33[k]
                  + f_0 * if__103[k];

        t_64[k] = -3.0 * gf_34[k]
                  + f_0 * if__104[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, gf_35, gf_36, gf_37, gf_38, gf_39, \
                         if__105, if__106, if__107, if__108, if__109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * gf_35[k]
                  + f_0 * if__105[k];

        t_66[k] = -3.0 * gf_36[k]
                  + f_0 * if__106[k];

        t_67[k] = -3.0 * gf_37[k]
                  + f_0 * if__107[k];

        t_68[k] = -3.0 * gf_38[k]
                  + f_0 * if__108[k];

        t_69[k] = -3.0 * gf_39[k]
                  + f_0 * if__109[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, gf_40, gf_41, gf_42, gf_43, gf_44, \
                         if__110, if__111, if__112, if__113, if__114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * gf_40[k]
                  + f_0 * if__110[k];

        t_71[k] = -2.0 * gf_41[k]
                  + f_0 * if__111[k];

        t_72[k] = -2.0 * gf_42[k]
                  + f_0 * if__112[k];

        t_73[k] = -2.0 * gf_43[k]
                  + f_0 * if__113[k];

        t_74[k] = -2.0 * gf_44[k]
                  + f_0 * if__114[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, gf_45, gf_46, gf_47, gf_48, gf_49, \
                         if__115, if__116, if__117, if__118, if__119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * gf_45[k]
                  + f_0 * if__115[k];

        t_76[k] = -2.0 * gf_46[k]
                  + f_0 * if__116[k];

        t_77[k] = -2.0 * gf_47[k]
                  + f_0 * if__117[k];

        t_78[k] = -2.0 * gf_48[k]
                  + f_0 * if__118[k];

        t_79[k] = -2.0 * gf_49[k]
                  + f_0 * if__119[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, gf_50, gf_51, gf_52, gf_53, gf_54, \
                         if__120, if__121, if__122, if__123, if__124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -gf_50[k]
                  + f_0 * if__120[k];

        t_81[k] = -gf_51[k]
                  + f_0 * if__121[k];

        t_82[k] = -gf_52[k]
                  + f_0 * if__122[k];

        t_83[k] = -gf_53[k]
                  + f_0 * if__123[k];

        t_84[k] = -gf_54[k]
                  + f_0 * if__124[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, gf_55, gf_56, gf_57, gf_58, gf_59, \
                         if__125, if__126, if__127, if__128, if__129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -gf_55[k]
                  + f_0 * if__125[k];

        t_86[k] = -gf_56[k]
                  + f_0 * if__126[k];

        t_87[k] = -gf_57[k]
                  + f_0 * if__127[k];

        t_88[k] = -gf_58[k]
                  + f_0 * if__128[k];

        t_89[k] = -gf_59[k]
                  + f_0 * if__129[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, if__130, if__131, \
                         if__132, if__133, if__134, if__135, if__136, \
                         if__137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * if__130[k];

        t_91[k] = f_0 * if__131[k];

        t_92[k] = f_0 * if__132[k];

        t_93[k] = f_0 * if__133[k];

        t_94[k] = f_0 * if__134[k];

        t_95[k] = f_0 * if__135[k];

        t_96[k] = f_0 * if__136[k];

        t_97[k] = f_0 * if__137[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, gf_60, gf_61, gf_62, gf_63, \
                         if__138, if__139, if__150, if__151, if__152, \
                         if__153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * if__138[k];

        t_99[k] = f_0 * if__139[k];

        t_100[k] = -4.0 * gf_60[k]
                   + f_0 * if__150[k];

        t_101[k] = -4.0 * gf_61[k]
                   + f_0 * if__151[k];

        t_102[k] = -4.0 * gf_62[k]
                   + f_0 * if__152[k];

        t_103[k] = -4.0 * gf_63[k]
                   + f_0 * if__153[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, gf_64, gf_65, gf_66, gf_67, gf_68, \
                         if__154, if__155, if__156, if__157, if__158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -4.0 * gf_64[k]
                   + f_0 * if__154[k];

        t_105[k] = -4.0 * gf_65[k]
                   + f_0 * if__155[k];

        t_106[k] = -4.0 * gf_66[k]
                   + f_0 * if__156[k];

        t_107[k] = -4.0 * gf_67[k]
                   + f_0 * if__157[k];

        t_108[k] = -4.0 * gf_68[k]
                   + f_0 * if__158[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, gf_69, gf_70, gf_71, gf_72, gf_73, \
                         if__159, if__160, if__161, if__162, if__163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -4.0 * gf_69[k]
                   + f_0 * if__159[k];

        t_110[k] = -3.0 * gf_70[k]
                   + f_0 * if__160[k];

        t_111[k] = -3.0 * gf_71[k]
                   + f_0 * if__161[k];

        t_112[k] = -3.0 * gf_72[k]
                   + f_0 * if__162[k];

        t_113[k] = -3.0 * gf_73[k]
                   + f_0 * if__163[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, gf_74, gf_75, gf_76, gf_77, gf_78, \
                         if__164, if__165, if__166, if__167, if__168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -3.0 * gf_74[k]
                   + f_0 * if__164[k];

        t_115[k] = -3.0 * gf_75[k]
                   + f_0 * if__165[k];

        t_116[k] = -3.0 * gf_76[k]
                   + f_0 * if__166[k];

        t_117[k] = -3.0 * gf_77[k]
                   + f_0 * if__167[k];

        t_118[k] = -3.0 * gf_78[k]
                   + f_0 * if__168[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, gf_79, gf_80, gf_81, gf_82, gf_83, \
                         if__169, if__170, if__171, if__172, if__173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -3.0 * gf_79[k]
                   + f_0 * if__169[k];

        t_120[k] = -2.0 * gf_80[k]
                   + f_0 * if__170[k];

        t_121[k] = -2.0 * gf_81[k]
                   + f_0 * if__171[k];

        t_122[k] = -2.0 * gf_82[k]
                   + f_0 * if__172[k];

        t_123[k] = -2.0 * gf_83[k]
                   + f_0 * if__173[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, gf_84, gf_85, gf_86, gf_87, gf_88, \
                         if__174, if__175, if__176, if__177, if__178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * gf_84[k]
                   + f_0 * if__174[k];

        t_125[k] = -2.0 * gf_85[k]
                   + f_0 * if__175[k];

        t_126[k] = -2.0 * gf_86[k]
                   + f_0 * if__176[k];

        t_127[k] = -2.0 * gf_87[k]
                   + f_0 * if__177[k];

        t_128[k] = -2.0 * gf_88[k]
                   + f_0 * if__178[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, gf_89, gf_90, gf_91, gf_92, gf_93, \
                         if__179, if__180, if__181, if__182, if__183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -2.0 * gf_89[k]
                   + f_0 * if__179[k];

        t_130[k] = -gf_90[k]
                   + f_0 * if__180[k];

        t_131[k] = -gf_91[k]
                   + f_0 * if__181[k];

        t_132[k] = -gf_92[k]
                   + f_0 * if__182[k];

        t_133[k] = -gf_93[k]
                   + f_0 * if__183[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, gf_94, gf_95, gf_96, gf_97, gf_98, \
                         if__184, if__185, if__186, if__187, if__188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -gf_94[k]
                   + f_0 * if__184[k];

        t_135[k] = -gf_95[k]
                   + f_0 * if__185[k];

        t_136[k] = -gf_96[k]
                   + f_0 * if__186[k];

        t_137[k] = -gf_97[k]
                   + f_0 * if__187[k];

        t_138[k] = -gf_98[k]
                   + f_0 * if__188[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, t_145, gf_99, if__189, \
                         if__190, if__191, if__192, if__193, if__194, \
                         if__195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -gf_99[k]
                   + f_0 * if__189[k];

        t_140[k] = f_0 * if__190[k];

        t_141[k] = f_0 * if__191[k];

        t_142[k] = f_0 * if__192[k];

        t_143[k] = f_0 * if__193[k];

        t_144[k] = f_0 * if__194[k];

        t_145[k] = f_0 * if__195[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, t_151, gf_100, gf_101, if__196, \
                         if__197, if__198, if__199, if__210, if__211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_0 * if__196[k];

        t_147[k] = f_0 * if__197[k];

        t_148[k] = f_0 * if__198[k];

        t_149[k] = f_0 * if__199[k];

        t_150[k] = -5.0 * gf_100[k]
                   + f_0 * if__210[k];

        t_151[k] = -5.0 * gf_101[k]
                   + f_0 * if__211[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, gf_102, gf_103, gf_104, gf_105, \
                         gf_106, if__212, if__213, if__214, if__215, \
                         if__216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -5.0 * gf_102[k]
                   + f_0 * if__212[k];

        t_153[k] = -5.0 * gf_103[k]
                   + f_0 * if__213[k];

        t_154[k] = -5.0 * gf_104[k]
                   + f_0 * if__214[k];

        t_155[k] = -5.0 * gf_105[k]
                   + f_0 * if__215[k];

        t_156[k] = -5.0 * gf_106[k]
                   + f_0 * if__216[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, gf_107, gf_108, gf_109, gf_110, \
                         gf_111, if__217, if__218, if__219, if__220, \
                         if__221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -5.0 * gf_107[k]
                   + f_0 * if__217[k];

        t_158[k] = -5.0 * gf_108[k]
                   + f_0 * if__218[k];

        t_159[k] = -5.0 * gf_109[k]
                   + f_0 * if__219[k];

        t_160[k] = -4.0 * gf_110[k]
                   + f_0 * if__220[k];

        t_161[k] = -4.0 * gf_111[k]
                   + f_0 * if__221[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, gf_112, gf_113, gf_114, gf_115, \
                         gf_116, if__222, if__223, if__224, if__225, \
                         if__226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -4.0 * gf_112[k]
                   + f_0 * if__222[k];

        t_163[k] = -4.0 * gf_113[k]
                   + f_0 * if__223[k];

        t_164[k] = -4.0 * gf_114[k]
                   + f_0 * if__224[k];

        t_165[k] = -4.0 * gf_115[k]
                   + f_0 * if__225[k];

        t_166[k] = -4.0 * gf_116[k]
                   + f_0 * if__226[k];
    }
}

static auto
compute_prim_geom_10_hf_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gf, const size_t if_,
                                                    const size_t ncols,
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

    const auto *gf_117 = buffer.data(gf + 117);
    const auto *gf_118 = buffer.data(gf + 118);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_120 = buffer.data(gf + 120);
    const auto *gf_121 = buffer.data(gf + 121);
    const auto *gf_122 = buffer.data(gf + 122);
    const auto *gf_123 = buffer.data(gf + 123);
    const auto *gf_124 = buffer.data(gf + 124);
    const auto *gf_125 = buffer.data(gf + 125);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_127 = buffer.data(gf + 127);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_130 = buffer.data(gf + 130);
    const auto *gf_131 = buffer.data(gf + 131);
    const auto *gf_132 = buffer.data(gf + 132);
    const auto *gf_133 = buffer.data(gf + 133);
    const auto *gf_134 = buffer.data(gf + 134);
    const auto *gf_135 = buffer.data(gf + 135);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_137 = buffer.data(gf + 137);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_140 = buffer.data(gf + 140);
    const auto *gf_141 = buffer.data(gf + 141);
    const auto *gf_142 = buffer.data(gf + 142);
    const auto *gf_143 = buffer.data(gf + 143);
    const auto *gf_144 = buffer.data(gf + 144);
    const auto *gf_145 = buffer.data(gf + 145);
    const auto *gf_146 = buffer.data(gf + 146);
    const auto *gf_147 = buffer.data(gf + 147);
    const auto *gf_148 = buffer.data(gf + 148);
    const auto *gf_149 = buffer.data(gf + 149);

    const auto *if__227 = buffer.data(if_ + 227);
    const auto *if__228 = buffer.data(if_ + 228);
    const auto *if__229 = buffer.data(if_ + 229);
    const auto *if__230 = buffer.data(if_ + 230);
    const auto *if__231 = buffer.data(if_ + 231);
    const auto *if__232 = buffer.data(if_ + 232);
    const auto *if__233 = buffer.data(if_ + 233);
    const auto *if__234 = buffer.data(if_ + 234);
    const auto *if__235 = buffer.data(if_ + 235);
    const auto *if__236 = buffer.data(if_ + 236);
    const auto *if__237 = buffer.data(if_ + 237);
    const auto *if__238 = buffer.data(if_ + 238);
    const auto *if__239 = buffer.data(if_ + 239);
    const auto *if__240 = buffer.data(if_ + 240);
    const auto *if__241 = buffer.data(if_ + 241);
    const auto *if__242 = buffer.data(if_ + 242);
    const auto *if__243 = buffer.data(if_ + 243);
    const auto *if__244 = buffer.data(if_ + 244);
    const auto *if__245 = buffer.data(if_ + 245);
    const auto *if__246 = buffer.data(if_ + 246);
    const auto *if__247 = buffer.data(if_ + 247);
    const auto *if__248 = buffer.data(if_ + 248);
    const auto *if__249 = buffer.data(if_ + 249);
    const auto *if__250 = buffer.data(if_ + 250);
    const auto *if__251 = buffer.data(if_ + 251);
    const auto *if__252 = buffer.data(if_ + 252);
    const auto *if__253 = buffer.data(if_ + 253);
    const auto *if__254 = buffer.data(if_ + 254);
    const auto *if__255 = buffer.data(if_ + 255);
    const auto *if__256 = buffer.data(if_ + 256);
    const auto *if__257 = buffer.data(if_ + 257);
    const auto *if__258 = buffer.data(if_ + 258);
    const auto *if__259 = buffer.data(if_ + 259);
    const auto *if__260 = buffer.data(if_ + 260);
    const auto *if__261 = buffer.data(if_ + 261);
    const auto *if__262 = buffer.data(if_ + 262);
    const auto *if__263 = buffer.data(if_ + 263);
    const auto *if__264 = buffer.data(if_ + 264);
    const auto *if__265 = buffer.data(if_ + 265);
    const auto *if__266 = buffer.data(if_ + 266);
    const auto *if__267 = buffer.data(if_ + 267);
    const auto *if__268 = buffer.data(if_ + 268);
    const auto *if__269 = buffer.data(if_ + 269);

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, gf_117, gf_118, gf_119, gf_120, \
                         gf_121, if__227, if__228, if__229, if__230, \
                         if__231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -4.0 * gf_117[k]
                   + f_0 * if__227[k];

        t_168[k] = -4.0 * gf_118[k]
                   + f_0 * if__228[k];

        t_169[k] = -4.0 * gf_119[k]
                   + f_0 * if__229[k];

        t_170[k] = -3.0 * gf_120[k]
                   + f_0 * if__230[k];

        t_171[k] = -3.0 * gf_121[k]
                   + f_0 * if__231[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, gf_122, gf_123, gf_124, gf_125, \
                         gf_126, if__232, if__233, if__234, if__235, \
                         if__236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -3.0 * gf_122[k]
                   + f_0 * if__232[k];

        t_173[k] = -3.0 * gf_123[k]
                   + f_0 * if__233[k];

        t_174[k] = -3.0 * gf_124[k]
                   + f_0 * if__234[k];

        t_175[k] = -3.0 * gf_125[k]
                   + f_0 * if__235[k];

        t_176[k] = -3.0 * gf_126[k]
                   + f_0 * if__236[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, gf_127, gf_128, gf_129, gf_130, \
                         gf_131, if__237, if__238, if__239, if__240, \
                         if__241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -3.0 * gf_127[k]
                   + f_0 * if__237[k];

        t_178[k] = -3.0 * gf_128[k]
                   + f_0 * if__238[k];

        t_179[k] = -3.0 * gf_129[k]
                   + f_0 * if__239[k];

        t_180[k] = -2.0 * gf_130[k]
                   + f_0 * if__240[k];

        t_181[k] = -2.0 * gf_131[k]
                   + f_0 * if__241[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, gf_132, gf_133, gf_134, gf_135, \
                         gf_136, if__242, if__243, if__244, if__245, \
                         if__246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * gf_132[k]
                   + f_0 * if__242[k];

        t_183[k] = -2.0 * gf_133[k]
                   + f_0 * if__243[k];

        t_184[k] = -2.0 * gf_134[k]
                   + f_0 * if__244[k];

        t_185[k] = -2.0 * gf_135[k]
                   + f_0 * if__245[k];

        t_186[k] = -2.0 * gf_136[k]
                   + f_0 * if__246[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, gf_137, gf_138, gf_139, gf_140, \
                         gf_141, if__247, if__248, if__249, if__250, \
                         if__251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * gf_137[k]
                   + f_0 * if__247[k];

        t_188[k] = -2.0 * gf_138[k]
                   + f_0 * if__248[k];

        t_189[k] = -2.0 * gf_139[k]
                   + f_0 * if__249[k];

        t_190[k] = -gf_140[k]
                   + f_0 * if__250[k];

        t_191[k] = -gf_141[k]
                   + f_0 * if__251[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, gf_142, gf_143, gf_144, gf_145, \
                         gf_146, if__252, if__253, if__254, if__255, \
                         if__256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -gf_142[k]
                   + f_0 * if__252[k];

        t_193[k] = -gf_143[k]
                   + f_0 * if__253[k];

        t_194[k] = -gf_144[k]
                   + f_0 * if__254[k];

        t_195[k] = -gf_145[k]
                   + f_0 * if__255[k];

        t_196[k] = -gf_146[k]
                   + f_0 * if__256[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, t_202, gf_147, gf_148, gf_149, \
                         if__257, if__258, if__259, if__260, if__261, \
                         if__262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -gf_147[k]
                   + f_0 * if__257[k];

        t_198[k] = -gf_148[k]
                   + f_0 * if__258[k];

        t_199[k] = -gf_149[k]
                   + f_0 * if__259[k];

        t_200[k] = f_0 * if__260[k];

        t_201[k] = f_0 * if__261[k];

        t_202[k] = f_0 * if__262[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, t_209, if__263, if__264, \
                         if__265, if__266, if__267, if__268, if__269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_0 * if__263[k];

        t_204[k] = f_0 * if__264[k];

        t_205[k] = f_0 * if__265[k];

        t_206[k] = f_0 * if__266[k];

        t_207[k] = f_0 * if__267[k];

        t_208[k] = f_0 * if__268[k];

        t_209[k] = f_0 * if__269[k];
    }
}

auto
compute_prim_geom_10_hf_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t gf, const size_t if_,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hf_electron_repulsion_1_piece0(buffer, target, gf, if_, ncols, alpha);

    compute_prim_geom_10_hf_electron_repulsion_1_piece1(buffer, target, gf, if_, ncols, alpha);
}

static auto
compute_prim_geom_10_hf_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gf, const size_t if_,
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

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_67 = buffer.data(gf + 67);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_71 = buffer.data(gf + 71);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_73 = buffer.data(gf + 73);
    const auto *gf_74 = buffer.data(gf + 74);
    const auto *gf_75 = buffer.data(gf + 75);
    const auto *gf_76 = buffer.data(gf + 76);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_81 = buffer.data(gf + 81);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_83 = buffer.data(gf + 83);
    const auto *gf_84 = buffer.data(gf + 84);
    const auto *gf_85 = buffer.data(gf + 85);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_89 = buffer.data(gf + 89);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_91 = buffer.data(gf + 91);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_93 = buffer.data(gf + 93);
    const auto *gf_94 = buffer.data(gf + 94);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_96 = buffer.data(gf + 96);
    const auto *gf_97 = buffer.data(gf + 97);
    const auto *gf_98 = buffer.data(gf + 98);
    const auto *gf_99 = buffer.data(gf + 99);
    const auto *gf_100 = buffer.data(gf + 100);
    const auto *gf_101 = buffer.data(gf + 101);
    const auto *gf_102 = buffer.data(gf + 102);
    const auto *gf_103 = buffer.data(gf + 103);
    const auto *gf_104 = buffer.data(gf + 104);
    const auto *gf_105 = buffer.data(gf + 105);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_107 = buffer.data(gf + 107);
    const auto *gf_108 = buffer.data(gf + 108);
    const auto *gf_109 = buffer.data(gf + 109);

    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__53 = buffer.data(if_ + 53);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__73 = buffer.data(if_ + 73);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__75 = buffer.data(if_ + 75);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__81 = buffer.data(if_ + 81);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__83 = buffer.data(if_ + 83);
    const auto *if__84 = buffer.data(if_ + 84);
    const auto *if__85 = buffer.data(if_ + 85);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__91 = buffer.data(if_ + 91);
    const auto *if__92 = buffer.data(if_ + 92);
    const auto *if__93 = buffer.data(if_ + 93);
    const auto *if__94 = buffer.data(if_ + 94);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__111 = buffer.data(if_ + 111);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__113 = buffer.data(if_ + 113);
    const auto *if__114 = buffer.data(if_ + 114);
    const auto *if__115 = buffer.data(if_ + 115);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__121 = buffer.data(if_ + 121);
    const auto *if__122 = buffer.data(if_ + 122);
    const auto *if__123 = buffer.data(if_ + 123);
    const auto *if__124 = buffer.data(if_ + 124);
    const auto *if__125 = buffer.data(if_ + 125);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__127 = buffer.data(if_ + 127);
    const auto *if__128 = buffer.data(if_ + 128);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__131 = buffer.data(if_ + 131);
    const auto *if__132 = buffer.data(if_ + 132);
    const auto *if__133 = buffer.data(if_ + 133);
    const auto *if__134 = buffer.data(if_ + 134);
    const auto *if__135 = buffer.data(if_ + 135);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__137 = buffer.data(if_ + 137);
    const auto *if__138 = buffer.data(if_ + 138);
    const auto *if__139 = buffer.data(if_ + 139);
    const auto *if__140 = buffer.data(if_ + 140);
    const auto *if__141 = buffer.data(if_ + 141);
    const auto *if__142 = buffer.data(if_ + 142);
    const auto *if__143 = buffer.data(if_ + 143);
    const auto *if__144 = buffer.data(if_ + 144);
    const auto *if__145 = buffer.data(if_ + 145);
    const auto *if__146 = buffer.data(if_ + 146);
    const auto *if__147 = buffer.data(if_ + 147);
    const auto *if__148 = buffer.data(if_ + 148);
    const auto *if__149 = buffer.data(if_ + 149);
    const auto *if__160 = buffer.data(if_ + 160);
    const auto *if__161 = buffer.data(if_ + 161);
    const auto *if__162 = buffer.data(if_ + 162);
    const auto *if__163 = buffer.data(if_ + 163);
    const auto *if__164 = buffer.data(if_ + 164);
    const auto *if__165 = buffer.data(if_ + 165);
    const auto *if__166 = buffer.data(if_ + 166);
    const auto *if__167 = buffer.data(if_ + 167);
    const auto *if__168 = buffer.data(if_ + 168);
    const auto *if__169 = buffer.data(if_ + 169);
    const auto *if__170 = buffer.data(if_ + 170);
    const auto *if__171 = buffer.data(if_ + 171);
    const auto *if__172 = buffer.data(if_ + 172);
    const auto *if__173 = buffer.data(if_ + 173);
    const auto *if__174 = buffer.data(if_ + 174);
    const auto *if__175 = buffer.data(if_ + 175);
    const auto *if__176 = buffer.data(if_ + 176);
    const auto *if__177 = buffer.data(if_ + 177);
    const auto *if__178 = buffer.data(if_ + 178);
    const auto *if__179 = buffer.data(if_ + 179);
    const auto *if__180 = buffer.data(if_ + 180);
    const auto *if__181 = buffer.data(if_ + 181);
    const auto *if__182 = buffer.data(if_ + 182);
    const auto *if__183 = buffer.data(if_ + 183);
    const auto *if__184 = buffer.data(if_ + 184);
    const auto *if__185 = buffer.data(if_ + 185);
    const auto *if__186 = buffer.data(if_ + 186);
    const auto *if__187 = buffer.data(if_ + 187);
    const auto *if__188 = buffer.data(if_ + 188);
    const auto *if__189 = buffer.data(if_ + 189);
    const auto *if__190 = buffer.data(if_ + 190);
    const auto *if__191 = buffer.data(if_ + 191);
    const auto *if__192 = buffer.data(if_ + 192);
    const auto *if__193 = buffer.data(if_ + 193);
    const auto *if__194 = buffer.data(if_ + 194);
    const auto *if__195 = buffer.data(if_ + 195);
    const auto *if__196 = buffer.data(if_ + 196);
    const auto *if__197 = buffer.data(if_ + 197);
    const auto *if__198 = buffer.data(if_ + 198);
    const auto *if__199 = buffer.data(if_ + 199);
    const auto *if__200 = buffer.data(if_ + 200);
    const auto *if__201 = buffer.data(if_ + 201);
    const auto *if__202 = buffer.data(if_ + 202);
    const auto *if__203 = buffer.data(if_ + 203);
    const auto *if__204 = buffer.data(if_ + 204);
    const auto *if__205 = buffer.data(if_ + 205);
    const auto *if__206 = buffer.data(if_ + 206);
    const auto *if__207 = buffer.data(if_ + 207);
    const auto *if__208 = buffer.data(if_ + 208);
    const auto *if__209 = buffer.data(if_ + 209);
    const auto *if__220 = buffer.data(if_ + 220);
    const auto *if__221 = buffer.data(if_ + 221);
    const auto *if__222 = buffer.data(if_ + 222);
    const auto *if__223 = buffer.data(if_ + 223);
    const auto *if__224 = buffer.data(if_ + 224);
    const auto *if__225 = buffer.data(if_ + 225);
    const auto *if__226 = buffer.data(if_ + 226);
    const auto *if__227 = buffer.data(if_ + 227);
    const auto *if__228 = buffer.data(if_ + 228);
    const auto *if__229 = buffer.data(if_ + 229);
    const auto *if__230 = buffer.data(if_ + 230);
    const auto *if__231 = buffer.data(if_ + 231);
    const auto *if__232 = buffer.data(if_ + 232);
    const auto *if__233 = buffer.data(if_ + 233);
    const auto *if__234 = buffer.data(if_ + 234);
    const auto *if__235 = buffer.data(if_ + 235);
    const auto *if__236 = buffer.data(if_ + 236);
    const auto *if__237 = buffer.data(if_ + 237);
    const auto *if__238 = buffer.data(if_ + 238);
    const auto *if__239 = buffer.data(if_ + 239);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, if__20, if__21, if__22, \
                         if__23, if__24, if__25, if__26, if__27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * if__20[k];

        t_1[k] = f_0 * if__21[k];

        t_2[k] = f_0 * if__22[k];

        t_3[k] = f_0 * if__23[k];

        t_4[k] = f_0 * if__24[k];

        t_5[k] = f_0 * if__25[k];

        t_6[k] = f_0 * if__26[k];

        t_7[k] = f_0 * if__27[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, if__28, if__29, if__40, \
                         if__41, if__42, if__43, if__44, if__45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * if__28[k];

        t_9[k] = f_0 * if__29[k];

        t_10[k] = f_0 * if__40[k];

        t_11[k] = f_0 * if__41[k];

        t_12[k] = f_0 * if__42[k];

        t_13[k] = f_0 * if__43[k];

        t_14[k] = f_0 * if__44[k];

        t_15[k] = f_0 * if__45[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, gf_0, gf_1, if__46, if__47, \
                         if__48, if__49, if__50, if__51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * if__46[k];

        t_17[k] = f_0 * if__47[k];

        t_18[k] = f_0 * if__48[k];

        t_19[k] = f_0 * if__49[k];

        t_20[k] = -gf_0[k]
                  + f_0 * if__50[k];

        t_21[k] = -gf_1[k]
                  + f_0 * if__51[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, gf_2, gf_3, gf_4, gf_5, gf_6, if__52, \
                         if__53, if__54, if__55, if__56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -gf_2[k]
                  + f_0 * if__52[k];

        t_23[k] = -gf_3[k]
                  + f_0 * if__53[k];

        t_24[k] = -gf_4[k]
                  + f_0 * if__54[k];

        t_25[k] = -gf_5[k]
                  + f_0 * if__55[k];

        t_26[k] = -gf_6[k]
                  + f_0 * if__56[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, gf_7, gf_8, gf_9, if__57, if__58, \
                         if__59, if__70, if__71, if__72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -gf_7[k]
                  + f_0 * if__57[k];

        t_28[k] = -gf_8[k]
                  + f_0 * if__58[k];

        t_29[k] = -gf_9[k]
                  + f_0 * if__59[k];

        t_30[k] = f_0 * if__70[k];

        t_31[k] = f_0 * if__71[k];

        t_32[k] = f_0 * if__72[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, t_39, if__73, if__74, if__75, \
                         if__76, if__77, if__78, if__79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * if__73[k];

        t_34[k] = f_0 * if__74[k];

        t_35[k] = f_0 * if__75[k];

        t_36[k] = f_0 * if__76[k];

        t_37[k] = f_0 * if__77[k];

        t_38[k] = f_0 * if__78[k];

        t_39[k] = f_0 * if__79[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, gf_10, gf_11, gf_12, gf_13, gf_14, \
                         if__80, if__81, if__82, if__83, if__84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -gf_10[k]
                  + f_0 * if__80[k];

        t_41[k] = -gf_11[k]
                  + f_0 * if__81[k];

        t_42[k] = -gf_12[k]
                  + f_0 * if__82[k];

        t_43[k] = -gf_13[k]
                  + f_0 * if__83[k];

        t_44[k] = -gf_14[k]
                  + f_0 * if__84[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, gf_15, gf_16, gf_17, gf_18, gf_19, \
                         if__85, if__86, if__87, if__88, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -gf_15[k]
                  + f_0 * if__85[k];

        t_46[k] = -gf_16[k]
                  + f_0 * if__86[k];

        t_47[k] = -gf_17[k]
                  + f_0 * if__87[k];

        t_48[k] = -gf_18[k]
                  + f_0 * if__88[k];

        t_49[k] = -gf_19[k]
                  + f_0 * if__89[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, gf_20, gf_21, gf_22, gf_23, gf_24, \
                         if__90, if__91, if__92, if__93, if__94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * gf_20[k]
                  + f_0 * if__90[k];

        t_51[k] = -2.0 * gf_21[k]
                  + f_0 * if__91[k];

        t_52[k] = -2.0 * gf_22[k]
                  + f_0 * if__92[k];

        t_53[k] = -2.0 * gf_23[k]
                  + f_0 * if__93[k];

        t_54[k] = -2.0 * gf_24[k]
                  + f_0 * if__94[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, gf_25, gf_26, gf_27, gf_28, gf_29, \
                         if__95, if__96, if__97, if__98, if__99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * gf_25[k]
                  + f_0 * if__95[k];

        t_56[k] = -2.0 * gf_26[k]
                  + f_0 * if__96[k];

        t_57[k] = -2.0 * gf_27[k]
                  + f_0 * if__97[k];

        t_58[k] = -2.0 * gf_28[k]
                  + f_0 * if__98[k];

        t_59[k] = -2.0 * gf_29[k]
                  + f_0 * if__99[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, t_66, t_67, if__110, if__111, \
                         if__112, if__113, if__114, if__115, if__116, \
                         if__117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * if__110[k];

        t_61[k] = f_0 * if__111[k];

        t_62[k] = f_0 * if__112[k];

        t_63[k] = f_0 * if__113[k];

        t_64[k] = f_0 * if__114[k];

        t_65[k] = f_0 * if__115[k];

        t_66[k] = f_0 * if__116[k];

        t_67[k] = f_0 * if__117[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, gf_30, gf_31, gf_32, gf_33, \
                         if__118, if__119, if__120, if__121, if__122, \
                         if__123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * if__118[k];

        t_69[k] = f_0 * if__119[k];

        t_70[k] = -gf_30[k]
                  + f_0 * if__120[k];

        t_71[k] = -gf_31[k]
                  + f_0 * if__121[k];

        t_72[k] = -gf_32[k]
                  + f_0 * if__122[k];

        t_73[k] = -gf_33[k]
                  + f_0 * if__123[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, gf_34, gf_35, gf_36, gf_37, gf_38, \
                         if__124, if__125, if__126, if__127, if__128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = -gf_34[k]
                  + f_0 * if__124[k];

        t_75[k] = -gf_35[k]
                  + f_0 * if__125[k];

        t_76[k] = -gf_36[k]
                  + f_0 * if__126[k];

        t_77[k] = -gf_37[k]
                  + f_0 * if__127[k];

        t_78[k] = -gf_38[k]
                  + f_0 * if__128[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, gf_39, gf_40, gf_41, gf_42, gf_43, \
                         if__129, if__130, if__131, if__132, if__133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -gf_39[k]
                  + f_0 * if__129[k];

        t_80[k] = -2.0 * gf_40[k]
                  + f_0 * if__130[k];

        t_81[k] = -2.0 * gf_41[k]
                  + f_0 * if__131[k];

        t_82[k] = -2.0 * gf_42[k]
                  + f_0 * if__132[k];

        t_83[k] = -2.0 * gf_43[k]
                  + f_0 * if__133[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, gf_44, gf_45, gf_46, gf_47, gf_48, \
                         if__134, if__135, if__136, if__137, if__138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -2.0 * gf_44[k]
                  + f_0 * if__134[k];

        t_85[k] = -2.0 * gf_45[k]
                  + f_0 * if__135[k];

        t_86[k] = -2.0 * gf_46[k]
                  + f_0 * if__136[k];

        t_87[k] = -2.0 * gf_47[k]
                  + f_0 * if__137[k];

        t_88[k] = -2.0 * gf_48[k]
                  + f_0 * if__138[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, gf_49, gf_50, gf_51, gf_52, gf_53, \
                         if__139, if__140, if__141, if__142, if__143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -2.0 * gf_49[k]
                  + f_0 * if__139[k];

        t_90[k] = -3.0 * gf_50[k]
                  + f_0 * if__140[k];

        t_91[k] = -3.0 * gf_51[k]
                  + f_0 * if__141[k];

        t_92[k] = -3.0 * gf_52[k]
                  + f_0 * if__142[k];

        t_93[k] = -3.0 * gf_53[k]
                  + f_0 * if__143[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, gf_54, gf_55, gf_56, gf_57, gf_58, \
                         if__144, if__145, if__146, if__147, if__148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -3.0 * gf_54[k]
                  + f_0 * if__144[k];

        t_95[k] = -3.0 * gf_55[k]
                  + f_0 * if__145[k];

        t_96[k] = -3.0 * gf_56[k]
                  + f_0 * if__146[k];

        t_97[k] = -3.0 * gf_57[k]
                  + f_0 * if__147[k];

        t_98[k] = -3.0 * gf_58[k]
                  + f_0 * if__148[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, t_105, gf_59, if__149, \
                         if__160, if__161, if__162, if__163, if__164, \
                         if__165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -3.0 * gf_59[k]
                  + f_0 * if__149[k];

        t_100[k] = f_0 * if__160[k];

        t_101[k] = f_0 * if__161[k];

        t_102[k] = f_0 * if__162[k];

        t_103[k] = f_0 * if__163[k];

        t_104[k] = f_0 * if__164[k];

        t_105[k] = f_0 * if__165[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, gf_60, gf_61, if__166, \
                         if__167, if__168, if__169, if__170, if__171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * if__166[k];

        t_107[k] = f_0 * if__167[k];

        t_108[k] = f_0 * if__168[k];

        t_109[k] = f_0 * if__169[k];

        t_110[k] = -gf_60[k]
                   + f_0 * if__170[k];

        t_111[k] = -gf_61[k]
                   + f_0 * if__171[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, gf_62, gf_63, gf_64, gf_65, gf_66, \
                         if__172, if__173, if__174, if__175, if__176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -gf_62[k]
                   + f_0 * if__172[k];

        t_113[k] = -gf_63[k]
                   + f_0 * if__173[k];

        t_114[k] = -gf_64[k]
                   + f_0 * if__174[k];

        t_115[k] = -gf_65[k]
                   + f_0 * if__175[k];

        t_116[k] = -gf_66[k]
                   + f_0 * if__176[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, gf_67, gf_68, gf_69, gf_70, gf_71, \
                         if__177, if__178, if__179, if__180, if__181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -gf_67[k]
                   + f_0 * if__177[k];

        t_118[k] = -gf_68[k]
                   + f_0 * if__178[k];

        t_119[k] = -gf_69[k]
                   + f_0 * if__179[k];

        t_120[k] = -2.0 * gf_70[k]
                   + f_0 * if__180[k];

        t_121[k] = -2.0 * gf_71[k]
                   + f_0 * if__181[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, gf_72, gf_73, gf_74, gf_75, gf_76, \
                         if__182, if__183, if__184, if__185, if__186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * gf_72[k]
                   + f_0 * if__182[k];

        t_123[k] = -2.0 * gf_73[k]
                   + f_0 * if__183[k];

        t_124[k] = -2.0 * gf_74[k]
                   + f_0 * if__184[k];

        t_125[k] = -2.0 * gf_75[k]
                   + f_0 * if__185[k];

        t_126[k] = -2.0 * gf_76[k]
                   + f_0 * if__186[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, gf_77, gf_78, gf_79, gf_80, gf_81, \
                         if__187, if__188, if__189, if__190, if__191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * gf_77[k]
                   + f_0 * if__187[k];

        t_128[k] = -2.0 * gf_78[k]
                   + f_0 * if__188[k];

        t_129[k] = -2.0 * gf_79[k]
                   + f_0 * if__189[k];

        t_130[k] = -3.0 * gf_80[k]
                   + f_0 * if__190[k];

        t_131[k] = -3.0 * gf_81[k]
                   + f_0 * if__191[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, gf_82, gf_83, gf_84, gf_85, gf_86, \
                         if__192, if__193, if__194, if__195, if__196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -3.0 * gf_82[k]
                   + f_0 * if__192[k];

        t_133[k] = -3.0 * gf_83[k]
                   + f_0 * if__193[k];

        t_134[k] = -3.0 * gf_84[k]
                   + f_0 * if__194[k];

        t_135[k] = -3.0 * gf_85[k]
                   + f_0 * if__195[k];

        t_136[k] = -3.0 * gf_86[k]
                   + f_0 * if__196[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, gf_87, gf_88, gf_89, gf_90, gf_91, \
                         if__197, if__198, if__199, if__200, if__201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -3.0 * gf_87[k]
                   + f_0 * if__197[k];

        t_138[k] = -3.0 * gf_88[k]
                   + f_0 * if__198[k];

        t_139[k] = -3.0 * gf_89[k]
                   + f_0 * if__199[k];

        t_140[k] = -4.0 * gf_90[k]
                   + f_0 * if__200[k];

        t_141[k] = -4.0 * gf_91[k]
                   + f_0 * if__201[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, gf_92, gf_93, gf_94, gf_95, gf_96, \
                         if__202, if__203, if__204, if__205, if__206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -4.0 * gf_92[k]
                   + f_0 * if__202[k];

        t_143[k] = -4.0 * gf_93[k]
                   + f_0 * if__203[k];

        t_144[k] = -4.0 * gf_94[k]
                   + f_0 * if__204[k];

        t_145[k] = -4.0 * gf_95[k]
                   + f_0 * if__205[k];

        t_146[k] = -4.0 * gf_96[k]
                   + f_0 * if__206[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, t_152, gf_97, gf_98, gf_99, \
                         if__207, if__208, if__209, if__220, if__221, \
                         if__222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -4.0 * gf_97[k]
                   + f_0 * if__207[k];

        t_148[k] = -4.0 * gf_98[k]
                   + f_0 * if__208[k];

        t_149[k] = -4.0 * gf_99[k]
                   + f_0 * if__209[k];

        t_150[k] = f_0 * if__220[k];

        t_151[k] = f_0 * if__221[k];

        t_152[k] = f_0 * if__222[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, t_158, t_159, if__223, if__224, \
                         if__225, if__226, if__227, if__228, if__229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_0 * if__223[k];

        t_154[k] = f_0 * if__224[k];

        t_155[k] = f_0 * if__225[k];

        t_156[k] = f_0 * if__226[k];

        t_157[k] = f_0 * if__227[k];

        t_158[k] = f_0 * if__228[k];

        t_159[k] = f_0 * if__229[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, gf_100, gf_101, gf_102, gf_103, \
                         gf_104, if__230, if__231, if__232, if__233, \
                         if__234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -gf_100[k]
                   + f_0 * if__230[k];

        t_161[k] = -gf_101[k]
                   + f_0 * if__231[k];

        t_162[k] = -gf_102[k]
                   + f_0 * if__232[k];

        t_163[k] = -gf_103[k]
                   + f_0 * if__233[k];

        t_164[k] = -gf_104[k]
                   + f_0 * if__234[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, gf_105, gf_106, gf_107, gf_108, \
                         gf_109, if__235, if__236, if__237, if__238, \
                         if__239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -gf_105[k]
                   + f_0 * if__235[k];

        t_166[k] = -gf_106[k]
                   + f_0 * if__236[k];

        t_167[k] = -gf_107[k]
                   + f_0 * if__237[k];

        t_168[k] = -gf_108[k]
                   + f_0 * if__238[k];

        t_169[k] = -gf_109[k]
                   + f_0 * if__239[k];
    }
}

static auto
compute_prim_geom_10_hf_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gf, const size_t if_,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gf_110 = buffer.data(gf + 110);
    const auto *gf_111 = buffer.data(gf + 111);
    const auto *gf_112 = buffer.data(gf + 112);
    const auto *gf_113 = buffer.data(gf + 113);
    const auto *gf_114 = buffer.data(gf + 114);
    const auto *gf_115 = buffer.data(gf + 115);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_117 = buffer.data(gf + 117);
    const auto *gf_118 = buffer.data(gf + 118);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_120 = buffer.data(gf + 120);
    const auto *gf_121 = buffer.data(gf + 121);
    const auto *gf_122 = buffer.data(gf + 122);
    const auto *gf_123 = buffer.data(gf + 123);
    const auto *gf_124 = buffer.data(gf + 124);
    const auto *gf_125 = buffer.data(gf + 125);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_127 = buffer.data(gf + 127);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_130 = buffer.data(gf + 130);
    const auto *gf_131 = buffer.data(gf + 131);
    const auto *gf_132 = buffer.data(gf + 132);
    const auto *gf_133 = buffer.data(gf + 133);
    const auto *gf_134 = buffer.data(gf + 134);
    const auto *gf_135 = buffer.data(gf + 135);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_137 = buffer.data(gf + 137);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_140 = buffer.data(gf + 140);
    const auto *gf_141 = buffer.data(gf + 141);
    const auto *gf_142 = buffer.data(gf + 142);
    const auto *gf_143 = buffer.data(gf + 143);
    const auto *gf_144 = buffer.data(gf + 144);
    const auto *gf_145 = buffer.data(gf + 145);
    const auto *gf_146 = buffer.data(gf + 146);
    const auto *gf_147 = buffer.data(gf + 147);
    const auto *gf_148 = buffer.data(gf + 148);
    const auto *gf_149 = buffer.data(gf + 149);

    const auto *if__240 = buffer.data(if_ + 240);
    const auto *if__241 = buffer.data(if_ + 241);
    const auto *if__242 = buffer.data(if_ + 242);
    const auto *if__243 = buffer.data(if_ + 243);
    const auto *if__244 = buffer.data(if_ + 244);
    const auto *if__245 = buffer.data(if_ + 245);
    const auto *if__246 = buffer.data(if_ + 246);
    const auto *if__247 = buffer.data(if_ + 247);
    const auto *if__248 = buffer.data(if_ + 248);
    const auto *if__249 = buffer.data(if_ + 249);
    const auto *if__250 = buffer.data(if_ + 250);
    const auto *if__251 = buffer.data(if_ + 251);
    const auto *if__252 = buffer.data(if_ + 252);
    const auto *if__253 = buffer.data(if_ + 253);
    const auto *if__254 = buffer.data(if_ + 254);
    const auto *if__255 = buffer.data(if_ + 255);
    const auto *if__256 = buffer.data(if_ + 256);
    const auto *if__257 = buffer.data(if_ + 257);
    const auto *if__258 = buffer.data(if_ + 258);
    const auto *if__259 = buffer.data(if_ + 259);
    const auto *if__260 = buffer.data(if_ + 260);
    const auto *if__261 = buffer.data(if_ + 261);
    const auto *if__262 = buffer.data(if_ + 262);
    const auto *if__263 = buffer.data(if_ + 263);
    const auto *if__264 = buffer.data(if_ + 264);
    const auto *if__265 = buffer.data(if_ + 265);
    const auto *if__266 = buffer.data(if_ + 266);
    const auto *if__267 = buffer.data(if_ + 267);
    const auto *if__268 = buffer.data(if_ + 268);
    const auto *if__269 = buffer.data(if_ + 269);
    const auto *if__270 = buffer.data(if_ + 270);
    const auto *if__271 = buffer.data(if_ + 271);
    const auto *if__272 = buffer.data(if_ + 272);
    const auto *if__273 = buffer.data(if_ + 273);
    const auto *if__274 = buffer.data(if_ + 274);
    const auto *if__275 = buffer.data(if_ + 275);
    const auto *if__276 = buffer.data(if_ + 276);
    const auto *if__277 = buffer.data(if_ + 277);
    const auto *if__278 = buffer.data(if_ + 278);
    const auto *if__279 = buffer.data(if_ + 279);

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, gf_110, gf_111, gf_112, gf_113, \
                         gf_114, if__240, if__241, if__242, if__243, \
                         if__244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -2.0 * gf_110[k]
                   + f_0 * if__240[k];

        t_171[k] = -2.0 * gf_111[k]
                   + f_0 * if__241[k];

        t_172[k] = -2.0 * gf_112[k]
                   + f_0 * if__242[k];

        t_173[k] = -2.0 * gf_113[k]
                   + f_0 * if__243[k];

        t_174[k] = -2.0 * gf_114[k]
                   + f_0 * if__244[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, gf_115, gf_116, gf_117, gf_118, \
                         gf_119, if__245, if__246, if__247, if__248, \
                         if__249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -2.0 * gf_115[k]
                   + f_0 * if__245[k];

        t_176[k] = -2.0 * gf_116[k]
                   + f_0 * if__246[k];

        t_177[k] = -2.0 * gf_117[k]
                   + f_0 * if__247[k];

        t_178[k] = -2.0 * gf_118[k]
                   + f_0 * if__248[k];

        t_179[k] = -2.0 * gf_119[k]
                   + f_0 * if__249[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, gf_120, gf_121, gf_122, gf_123, \
                         gf_124, if__250, if__251, if__252, if__253, \
                         if__254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * gf_120[k]
                   + f_0 * if__250[k];

        t_181[k] = -3.0 * gf_121[k]
                   + f_0 * if__251[k];

        t_182[k] = -3.0 * gf_122[k]
                   + f_0 * if__252[k];

        t_183[k] = -3.0 * gf_123[k]
                   + f_0 * if__253[k];

        t_184[k] = -3.0 * gf_124[k]
                   + f_0 * if__254[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, gf_125, gf_126, gf_127, gf_128, \
                         gf_129, if__255, if__256, if__257, if__258, \
                         if__259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * gf_125[k]
                   + f_0 * if__255[k];

        t_186[k] = -3.0 * gf_126[k]
                   + f_0 * if__256[k];

        t_187[k] = -3.0 * gf_127[k]
                   + f_0 * if__257[k];

        t_188[k] = -3.0 * gf_128[k]
                   + f_0 * if__258[k];

        t_189[k] = -3.0 * gf_129[k]
                   + f_0 * if__259[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, gf_130, gf_131, gf_132, gf_133, \
                         gf_134, if__260, if__261, if__262, if__263, \
                         if__264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -4.0 * gf_130[k]
                   + f_0 * if__260[k];

        t_191[k] = -4.0 * gf_131[k]
                   + f_0 * if__261[k];

        t_192[k] = -4.0 * gf_132[k]
                   + f_0 * if__262[k];

        t_193[k] = -4.0 * gf_133[k]
                   + f_0 * if__263[k];

        t_194[k] = -4.0 * gf_134[k]
                   + f_0 * if__264[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, gf_135, gf_136, gf_137, gf_138, \
                         gf_139, if__265, if__266, if__267, if__268, \
                         if__269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -4.0 * gf_135[k]
                   + f_0 * if__265[k];

        t_196[k] = -4.0 * gf_136[k]
                   + f_0 * if__266[k];

        t_197[k] = -4.0 * gf_137[k]
                   + f_0 * if__267[k];

        t_198[k] = -4.0 * gf_138[k]
                   + f_0 * if__268[k];

        t_199[k] = -4.0 * gf_139[k]
                   + f_0 * if__269[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, gf_140, gf_141, gf_142, gf_143, \
                         gf_144, if__270, if__271, if__272, if__273, \
                         if__274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -5.0 * gf_140[k]
                   + f_0 * if__270[k];

        t_201[k] = -5.0 * gf_141[k]
                   + f_0 * if__271[k];

        t_202[k] = -5.0 * gf_142[k]
                   + f_0 * if__272[k];

        t_203[k] = -5.0 * gf_143[k]
                   + f_0 * if__273[k];

        t_204[k] = -5.0 * gf_144[k]
                   + f_0 * if__274[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, gf_145, gf_146, gf_147, gf_148, \
                         gf_149, if__275, if__276, if__277, if__278, \
                         if__279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -5.0 * gf_145[k]
                   + f_0 * if__275[k];

        t_206[k] = -5.0 * gf_146[k]
                   + f_0 * if__276[k];

        t_207[k] = -5.0 * gf_147[k]
                   + f_0 * if__277[k];

        t_208[k] = -5.0 * gf_148[k]
                   + f_0 * if__278[k];

        t_209[k] = -5.0 * gf_149[k]
                   + f_0 * if__279[k];
    }
}

auto
compute_prim_geom_10_hf_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t gf, const size_t if_,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hf_electron_repulsion_2_piece0(buffer, target, gf, if_, ncols, alpha);

    compute_prim_geom_10_hf_electron_repulsion_2_piece1(buffer, target, gf, if_, ncols, alpha);
}

}  // namespace simdt2ceri
