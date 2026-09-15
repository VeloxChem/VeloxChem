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


#include "SimdElectronRepulsionGeom10VrrRecKF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_kf_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t if_, const size_t lf,
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

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);
    const auto *lf_108 = buffer.data(lf + 108);
    const auto *lf_109 = buffer.data(lf + 109);
    const auto *lf_110 = buffer.data(lf + 110);
    const auto *lf_111 = buffer.data(lf + 111);
    const auto *lf_112 = buffer.data(lf + 112);
    const auto *lf_113 = buffer.data(lf + 113);
    const auto *lf_114 = buffer.data(lf + 114);
    const auto *lf_115 = buffer.data(lf + 115);
    const auto *lf_116 = buffer.data(lf + 116);
    const auto *lf_117 = buffer.data(lf + 117);
    const auto *lf_118 = buffer.data(lf + 118);
    const auto *lf_119 = buffer.data(lf + 119);
    const auto *lf_120 = buffer.data(lf + 120);
    const auto *lf_121 = buffer.data(lf + 121);
    const auto *lf_122 = buffer.data(lf + 122);
    const auto *lf_123 = buffer.data(lf + 123);
    const auto *lf_124 = buffer.data(lf + 124);
    const auto *lf_125 = buffer.data(lf + 125);
    const auto *lf_126 = buffer.data(lf + 126);
    const auto *lf_127 = buffer.data(lf + 127);
    const auto *lf_128 = buffer.data(lf + 128);
    const auto *lf_129 = buffer.data(lf + 129);
    const auto *lf_130 = buffer.data(lf + 130);
    const auto *lf_131 = buffer.data(lf + 131);
    const auto *lf_132 = buffer.data(lf + 132);
    const auto *lf_133 = buffer.data(lf + 133);
    const auto *lf_134 = buffer.data(lf + 134);
    const auto *lf_135 = buffer.data(lf + 135);
    const auto *lf_136 = buffer.data(lf + 136);
    const auto *lf_137 = buffer.data(lf + 137);
    const auto *lf_138 = buffer.data(lf + 138);
    const auto *lf_139 = buffer.data(lf + 139);
    const auto *lf_140 = buffer.data(lf + 140);
    const auto *lf_141 = buffer.data(lf + 141);
    const auto *lf_142 = buffer.data(lf + 142);
    const auto *lf_143 = buffer.data(lf + 143);
    const auto *lf_144 = buffer.data(lf + 144);
    const auto *lf_145 = buffer.data(lf + 145);
    const auto *lf_146 = buffer.data(lf + 146);
    const auto *lf_147 = buffer.data(lf + 147);
    const auto *lf_148 = buffer.data(lf + 148);
    const auto *lf_149 = buffer.data(lf + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, if__0, if__1, if__2, if__3, if__4, lf_0, \
                         lf_1, lf_2, lf_3, lf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -7.0 * if__0[k]
                 + f_0 * lf_0[k];

        t_1[k] = -7.0 * if__1[k]
                 + f_0 * lf_1[k];

        t_2[k] = -7.0 * if__2[k]
                 + f_0 * lf_2[k];

        t_3[k] = -7.0 * if__3[k]
                 + f_0 * lf_3[k];

        t_4[k] = -7.0 * if__4[k]
                 + f_0 * lf_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, if__5, if__6, if__7, if__8, if__9, lf_5, \
                         lf_6, lf_7, lf_8, lf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -7.0 * if__5[k]
                 + f_0 * lf_5[k];

        t_6[k] = -7.0 * if__6[k]
                 + f_0 * lf_6[k];

        t_7[k] = -7.0 * if__7[k]
                 + f_0 * lf_7[k];

        t_8[k] = -7.0 * if__8[k]
                 + f_0 * lf_8[k];

        t_9[k] = -7.0 * if__9[k]
                 + f_0 * lf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, if__10, if__11, if__12, if__13, if__14, \
                         lf_10, lf_11, lf_12, lf_13, lf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -6.0 * if__10[k]
                  + f_0 * lf_10[k];

        t_11[k] = -6.0 * if__11[k]
                  + f_0 * lf_11[k];

        t_12[k] = -6.0 * if__12[k]
                  + f_0 * lf_12[k];

        t_13[k] = -6.0 * if__13[k]
                  + f_0 * lf_13[k];

        t_14[k] = -6.0 * if__14[k]
                  + f_0 * lf_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, if__15, if__16, if__17, if__18, if__19, \
                         lf_15, lf_16, lf_17, lf_18, lf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -6.0 * if__15[k]
                  + f_0 * lf_15[k];

        t_16[k] = -6.0 * if__16[k]
                  + f_0 * lf_16[k];

        t_17[k] = -6.0 * if__17[k]
                  + f_0 * lf_17[k];

        t_18[k] = -6.0 * if__18[k]
                  + f_0 * lf_18[k];

        t_19[k] = -6.0 * if__19[k]
                  + f_0 * lf_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, if__20, if__21, if__22, if__23, if__24, \
                         lf_20, lf_21, lf_22, lf_23, lf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -6.0 * if__20[k]
                  + f_0 * lf_20[k];

        t_21[k] = -6.0 * if__21[k]
                  + f_0 * lf_21[k];

        t_22[k] = -6.0 * if__22[k]
                  + f_0 * lf_22[k];

        t_23[k] = -6.0 * if__23[k]
                  + f_0 * lf_23[k];

        t_24[k] = -6.0 * if__24[k]
                  + f_0 * lf_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, if__25, if__26, if__27, if__28, if__29, \
                         lf_25, lf_26, lf_27, lf_28, lf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -6.0 * if__25[k]
                  + f_0 * lf_25[k];

        t_26[k] = -6.0 * if__26[k]
                  + f_0 * lf_26[k];

        t_27[k] = -6.0 * if__27[k]
                  + f_0 * lf_27[k];

        t_28[k] = -6.0 * if__28[k]
                  + f_0 * lf_28[k];

        t_29[k] = -6.0 * if__29[k]
                  + f_0 * lf_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, if__30, if__31, if__32, if__33, if__34, \
                         lf_30, lf_31, lf_32, lf_33, lf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -5.0 * if__30[k]
                  + f_0 * lf_30[k];

        t_31[k] = -5.0 * if__31[k]
                  + f_0 * lf_31[k];

        t_32[k] = -5.0 * if__32[k]
                  + f_0 * lf_32[k];

        t_33[k] = -5.0 * if__33[k]
                  + f_0 * lf_33[k];

        t_34[k] = -5.0 * if__34[k]
                  + f_0 * lf_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, if__35, if__36, if__37, if__38, if__39, \
                         lf_35, lf_36, lf_37, lf_38, lf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -5.0 * if__35[k]
                  + f_0 * lf_35[k];

        t_36[k] = -5.0 * if__36[k]
                  + f_0 * lf_36[k];

        t_37[k] = -5.0 * if__37[k]
                  + f_0 * lf_37[k];

        t_38[k] = -5.0 * if__38[k]
                  + f_0 * lf_38[k];

        t_39[k] = -5.0 * if__39[k]
                  + f_0 * lf_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, if__40, if__41, if__42, if__43, if__44, \
                         lf_40, lf_41, lf_42, lf_43, lf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -5.0 * if__40[k]
                  + f_0 * lf_40[k];

        t_41[k] = -5.0 * if__41[k]
                  + f_0 * lf_41[k];

        t_42[k] = -5.0 * if__42[k]
                  + f_0 * lf_42[k];

        t_43[k] = -5.0 * if__43[k]
                  + f_0 * lf_43[k];

        t_44[k] = -5.0 * if__44[k]
                  + f_0 * lf_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, if__45, if__46, if__47, if__48, if__49, \
                         lf_45, lf_46, lf_47, lf_48, lf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -5.0 * if__45[k]
                  + f_0 * lf_45[k];

        t_46[k] = -5.0 * if__46[k]
                  + f_0 * lf_46[k];

        t_47[k] = -5.0 * if__47[k]
                  + f_0 * lf_47[k];

        t_48[k] = -5.0 * if__48[k]
                  + f_0 * lf_48[k];

        t_49[k] = -5.0 * if__49[k]
                  + f_0 * lf_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, if__50, if__51, if__52, if__53, if__54, \
                         lf_50, lf_51, lf_52, lf_53, lf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -5.0 * if__50[k]
                  + f_0 * lf_50[k];

        t_51[k] = -5.0 * if__51[k]
                  + f_0 * lf_51[k];

        t_52[k] = -5.0 * if__52[k]
                  + f_0 * lf_52[k];

        t_53[k] = -5.0 * if__53[k]
                  + f_0 * lf_53[k];

        t_54[k] = -5.0 * if__54[k]
                  + f_0 * lf_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, if__55, if__56, if__57, if__58, if__59, \
                         lf_55, lf_56, lf_57, lf_58, lf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -5.0 * if__55[k]
                  + f_0 * lf_55[k];

        t_56[k] = -5.0 * if__56[k]
                  + f_0 * lf_56[k];

        t_57[k] = -5.0 * if__57[k]
                  + f_0 * lf_57[k];

        t_58[k] = -5.0 * if__58[k]
                  + f_0 * lf_58[k];

        t_59[k] = -5.0 * if__59[k]
                  + f_0 * lf_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, if__60, if__61, if__62, if__63, if__64, \
                         lf_60, lf_61, lf_62, lf_63, lf_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -4.0 * if__60[k]
                  + f_0 * lf_60[k];

        t_61[k] = -4.0 * if__61[k]
                  + f_0 * lf_61[k];

        t_62[k] = -4.0 * if__62[k]
                  + f_0 * lf_62[k];

        t_63[k] = -4.0 * if__63[k]
                  + f_0 * lf_63[k];

        t_64[k] = -4.0 * if__64[k]
                  + f_0 * lf_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, if__65, if__66, if__67, if__68, if__69, \
                         lf_65, lf_66, lf_67, lf_68, lf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -4.0 * if__65[k]
                  + f_0 * lf_65[k];

        t_66[k] = -4.0 * if__66[k]
                  + f_0 * lf_66[k];

        t_67[k] = -4.0 * if__67[k]
                  + f_0 * lf_67[k];

        t_68[k] = -4.0 * if__68[k]
                  + f_0 * lf_68[k];

        t_69[k] = -4.0 * if__69[k]
                  + f_0 * lf_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, if__70, if__71, if__72, if__73, if__74, \
                         lf_70, lf_71, lf_72, lf_73, lf_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -4.0 * if__70[k]
                  + f_0 * lf_70[k];

        t_71[k] = -4.0 * if__71[k]
                  + f_0 * lf_71[k];

        t_72[k] = -4.0 * if__72[k]
                  + f_0 * lf_72[k];

        t_73[k] = -4.0 * if__73[k]
                  + f_0 * lf_73[k];

        t_74[k] = -4.0 * if__74[k]
                  + f_0 * lf_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, if__75, if__76, if__77, if__78, if__79, \
                         lf_75, lf_76, lf_77, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -4.0 * if__75[k]
                  + f_0 * lf_75[k];

        t_76[k] = -4.0 * if__76[k]
                  + f_0 * lf_76[k];

        t_77[k] = -4.0 * if__77[k]
                  + f_0 * lf_77[k];

        t_78[k] = -4.0 * if__78[k]
                  + f_0 * lf_78[k];

        t_79[k] = -4.0 * if__79[k]
                  + f_0 * lf_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, if__80, if__81, if__82, if__83, if__84, \
                         lf_80, lf_81, lf_82, lf_83, lf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -4.0 * if__80[k]
                  + f_0 * lf_80[k];

        t_81[k] = -4.0 * if__81[k]
                  + f_0 * lf_81[k];

        t_82[k] = -4.0 * if__82[k]
                  + f_0 * lf_82[k];

        t_83[k] = -4.0 * if__83[k]
                  + f_0 * lf_83[k];

        t_84[k] = -4.0 * if__84[k]
                  + f_0 * lf_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, if__85, if__86, if__87, if__88, if__89, \
                         lf_85, lf_86, lf_87, lf_88, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * if__85[k]
                  + f_0 * lf_85[k];

        t_86[k] = -4.0 * if__86[k]
                  + f_0 * lf_86[k];

        t_87[k] = -4.0 * if__87[k]
                  + f_0 * lf_87[k];

        t_88[k] = -4.0 * if__88[k]
                  + f_0 * lf_88[k];

        t_89[k] = -4.0 * if__89[k]
                  + f_0 * lf_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, if__90, if__91, if__92, if__93, if__94, \
                         lf_90, lf_91, lf_92, lf_93, lf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -4.0 * if__90[k]
                  + f_0 * lf_90[k];

        t_91[k] = -4.0 * if__91[k]
                  + f_0 * lf_91[k];

        t_92[k] = -4.0 * if__92[k]
                  + f_0 * lf_92[k];

        t_93[k] = -4.0 * if__93[k]
                  + f_0 * lf_93[k];

        t_94[k] = -4.0 * if__94[k]
                  + f_0 * lf_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, if__95, if__96, if__97, if__98, if__99, \
                         lf_95, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -4.0 * if__95[k]
                  + f_0 * lf_95[k];

        t_96[k] = -4.0 * if__96[k]
                  + f_0 * lf_96[k];

        t_97[k] = -4.0 * if__97[k]
                  + f_0 * lf_97[k];

        t_98[k] = -4.0 * if__98[k]
                  + f_0 * lf_98[k];

        t_99[k] = -4.0 * if__99[k]
                  + f_0 * lf_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, if__100, if__101, if__102, \
                         if__103, if__104, lf_100, lf_101, lf_102, lf_103, \
                         lf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * if__100[k]
                   + f_0 * lf_100[k];

        t_101[k] = -3.0 * if__101[k]
                   + f_0 * lf_101[k];

        t_102[k] = -3.0 * if__102[k]
                   + f_0 * lf_102[k];

        t_103[k] = -3.0 * if__103[k]
                   + f_0 * lf_103[k];

        t_104[k] = -3.0 * if__104[k]
                   + f_0 * lf_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, if__105, if__106, if__107, \
                         if__108, if__109, lf_105, lf_106, lf_107, lf_108, \
                         lf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -3.0 * if__105[k]
                   + f_0 * lf_105[k];

        t_106[k] = -3.0 * if__106[k]
                   + f_0 * lf_106[k];

        t_107[k] = -3.0 * if__107[k]
                   + f_0 * lf_107[k];

        t_108[k] = -3.0 * if__108[k]
                   + f_0 * lf_108[k];

        t_109[k] = -3.0 * if__109[k]
                   + f_0 * lf_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, if__110, if__111, if__112, \
                         if__113, if__114, lf_110, lf_111, lf_112, lf_113, \
                         lf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -3.0 * if__110[k]
                   + f_0 * lf_110[k];

        t_111[k] = -3.0 * if__111[k]
                   + f_0 * lf_111[k];

        t_112[k] = -3.0 * if__112[k]
                   + f_0 * lf_112[k];

        t_113[k] = -3.0 * if__113[k]
                   + f_0 * lf_113[k];

        t_114[k] = -3.0 * if__114[k]
                   + f_0 * lf_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, if__115, if__116, if__117, \
                         if__118, if__119, lf_115, lf_116, lf_117, lf_118, \
                         lf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -3.0 * if__115[k]
                   + f_0 * lf_115[k];

        t_116[k] = -3.0 * if__116[k]
                   + f_0 * lf_116[k];

        t_117[k] = -3.0 * if__117[k]
                   + f_0 * lf_117[k];

        t_118[k] = -3.0 * if__118[k]
                   + f_0 * lf_118[k];

        t_119[k] = -3.0 * if__119[k]
                   + f_0 * lf_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, if__120, if__121, if__122, \
                         if__123, if__124, lf_120, lf_121, lf_122, lf_123, \
                         lf_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -3.0 * if__120[k]
                   + f_0 * lf_120[k];

        t_121[k] = -3.0 * if__121[k]
                   + f_0 * lf_121[k];

        t_122[k] = -3.0 * if__122[k]
                   + f_0 * lf_122[k];

        t_123[k] = -3.0 * if__123[k]
                   + f_0 * lf_123[k];

        t_124[k] = -3.0 * if__124[k]
                   + f_0 * lf_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, if__125, if__126, if__127, \
                         if__128, if__129, lf_125, lf_126, lf_127, lf_128, \
                         lf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -3.0 * if__125[k]
                   + f_0 * lf_125[k];

        t_126[k] = -3.0 * if__126[k]
                   + f_0 * lf_126[k];

        t_127[k] = -3.0 * if__127[k]
                   + f_0 * lf_127[k];

        t_128[k] = -3.0 * if__128[k]
                   + f_0 * lf_128[k];

        t_129[k] = -3.0 * if__129[k]
                   + f_0 * lf_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, if__130, if__131, if__132, \
                         if__133, if__134, lf_130, lf_131, lf_132, lf_133, \
                         lf_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -3.0 * if__130[k]
                   + f_0 * lf_130[k];

        t_131[k] = -3.0 * if__131[k]
                   + f_0 * lf_131[k];

        t_132[k] = -3.0 * if__132[k]
                   + f_0 * lf_132[k];

        t_133[k] = -3.0 * if__133[k]
                   + f_0 * lf_133[k];

        t_134[k] = -3.0 * if__134[k]
                   + f_0 * lf_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, if__135, if__136, if__137, \
                         if__138, if__139, lf_135, lf_136, lf_137, lf_138, \
                         lf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -3.0 * if__135[k]
                   + f_0 * lf_135[k];

        t_136[k] = -3.0 * if__136[k]
                   + f_0 * lf_136[k];

        t_137[k] = -3.0 * if__137[k]
                   + f_0 * lf_137[k];

        t_138[k] = -3.0 * if__138[k]
                   + f_0 * lf_138[k];

        t_139[k] = -3.0 * if__139[k]
                   + f_0 * lf_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, if__140, if__141, if__142, \
                         if__143, if__144, lf_140, lf_141, lf_142, lf_143, \
                         lf_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -3.0 * if__140[k]
                   + f_0 * lf_140[k];

        t_141[k] = -3.0 * if__141[k]
                   + f_0 * lf_141[k];

        t_142[k] = -3.0 * if__142[k]
                   + f_0 * lf_142[k];

        t_143[k] = -3.0 * if__143[k]
                   + f_0 * lf_143[k];

        t_144[k] = -3.0 * if__144[k]
                   + f_0 * lf_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, if__145, if__146, if__147, \
                         if__148, if__149, lf_145, lf_146, lf_147, lf_148, \
                         lf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -3.0 * if__145[k]
                   + f_0 * lf_145[k];

        t_146[k] = -3.0 * if__146[k]
                   + f_0 * lf_146[k];

        t_147[k] = -3.0 * if__147[k]
                   + f_0 * lf_147[k];

        t_148[k] = -3.0 * if__148[k]
                   + f_0 * lf_148[k];

        t_149[k] = -3.0 * if__149[k]
                   + f_0 * lf_149[k];
    }
}

static auto
compute_prim_geom_10_kf_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t if_, const size_t lf,
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

    const auto *lf_150 = buffer.data(lf + 150);
    const auto *lf_151 = buffer.data(lf + 151);
    const auto *lf_152 = buffer.data(lf + 152);
    const auto *lf_153 = buffer.data(lf + 153);
    const auto *lf_154 = buffer.data(lf + 154);
    const auto *lf_155 = buffer.data(lf + 155);
    const auto *lf_156 = buffer.data(lf + 156);
    const auto *lf_157 = buffer.data(lf + 157);
    const auto *lf_158 = buffer.data(lf + 158);
    const auto *lf_159 = buffer.data(lf + 159);
    const auto *lf_160 = buffer.data(lf + 160);
    const auto *lf_161 = buffer.data(lf + 161);
    const auto *lf_162 = buffer.data(lf + 162);
    const auto *lf_163 = buffer.data(lf + 163);
    const auto *lf_164 = buffer.data(lf + 164);
    const auto *lf_165 = buffer.data(lf + 165);
    const auto *lf_166 = buffer.data(lf + 166);
    const auto *lf_167 = buffer.data(lf + 167);
    const auto *lf_168 = buffer.data(lf + 168);
    const auto *lf_169 = buffer.data(lf + 169);
    const auto *lf_170 = buffer.data(lf + 170);
    const auto *lf_171 = buffer.data(lf + 171);
    const auto *lf_172 = buffer.data(lf + 172);
    const auto *lf_173 = buffer.data(lf + 173);
    const auto *lf_174 = buffer.data(lf + 174);
    const auto *lf_175 = buffer.data(lf + 175);
    const auto *lf_176 = buffer.data(lf + 176);
    const auto *lf_177 = buffer.data(lf + 177);
    const auto *lf_178 = buffer.data(lf + 178);
    const auto *lf_179 = buffer.data(lf + 179);
    const auto *lf_180 = buffer.data(lf + 180);
    const auto *lf_181 = buffer.data(lf + 181);
    const auto *lf_182 = buffer.data(lf + 182);
    const auto *lf_183 = buffer.data(lf + 183);
    const auto *lf_184 = buffer.data(lf + 184);
    const auto *lf_185 = buffer.data(lf + 185);
    const auto *lf_186 = buffer.data(lf + 186);
    const auto *lf_187 = buffer.data(lf + 187);
    const auto *lf_188 = buffer.data(lf + 188);
    const auto *lf_189 = buffer.data(lf + 189);
    const auto *lf_190 = buffer.data(lf + 190);
    const auto *lf_191 = buffer.data(lf + 191);
    const auto *lf_192 = buffer.data(lf + 192);
    const auto *lf_193 = buffer.data(lf + 193);
    const auto *lf_194 = buffer.data(lf + 194);
    const auto *lf_195 = buffer.data(lf + 195);
    const auto *lf_196 = buffer.data(lf + 196);
    const auto *lf_197 = buffer.data(lf + 197);
    const auto *lf_198 = buffer.data(lf + 198);
    const auto *lf_199 = buffer.data(lf + 199);
    const auto *lf_200 = buffer.data(lf + 200);
    const auto *lf_201 = buffer.data(lf + 201);
    const auto *lf_202 = buffer.data(lf + 202);
    const auto *lf_203 = buffer.data(lf + 203);
    const auto *lf_204 = buffer.data(lf + 204);
    const auto *lf_205 = buffer.data(lf + 205);
    const auto *lf_206 = buffer.data(lf + 206);
    const auto *lf_207 = buffer.data(lf + 207);
    const auto *lf_208 = buffer.data(lf + 208);
    const auto *lf_209 = buffer.data(lf + 209);
    const auto *lf_210 = buffer.data(lf + 210);
    const auto *lf_211 = buffer.data(lf + 211);
    const auto *lf_212 = buffer.data(lf + 212);
    const auto *lf_213 = buffer.data(lf + 213);
    const auto *lf_214 = buffer.data(lf + 214);
    const auto *lf_215 = buffer.data(lf + 215);
    const auto *lf_216 = buffer.data(lf + 216);
    const auto *lf_217 = buffer.data(lf + 217);
    const auto *lf_218 = buffer.data(lf + 218);
    const auto *lf_219 = buffer.data(lf + 219);
    const auto *lf_220 = buffer.data(lf + 220);
    const auto *lf_221 = buffer.data(lf + 221);
    const auto *lf_222 = buffer.data(lf + 222);
    const auto *lf_223 = buffer.data(lf + 223);
    const auto *lf_224 = buffer.data(lf + 224);
    const auto *lf_225 = buffer.data(lf + 225);
    const auto *lf_226 = buffer.data(lf + 226);
    const auto *lf_227 = buffer.data(lf + 227);
    const auto *lf_228 = buffer.data(lf + 228);
    const auto *lf_229 = buffer.data(lf + 229);
    const auto *lf_230 = buffer.data(lf + 230);
    const auto *lf_231 = buffer.data(lf + 231);
    const auto *lf_232 = buffer.data(lf + 232);
    const auto *lf_233 = buffer.data(lf + 233);
    const auto *lf_234 = buffer.data(lf + 234);
    const auto *lf_235 = buffer.data(lf + 235);
    const auto *lf_236 = buffer.data(lf + 236);
    const auto *lf_237 = buffer.data(lf + 237);
    const auto *lf_238 = buffer.data(lf + 238);
    const auto *lf_239 = buffer.data(lf + 239);
    const auto *lf_240 = buffer.data(lf + 240);
    const auto *lf_241 = buffer.data(lf + 241);
    const auto *lf_242 = buffer.data(lf + 242);
    const auto *lf_243 = buffer.data(lf + 243);
    const auto *lf_244 = buffer.data(lf + 244);
    const auto *lf_245 = buffer.data(lf + 245);
    const auto *lf_246 = buffer.data(lf + 246);
    const auto *lf_247 = buffer.data(lf + 247);
    const auto *lf_248 = buffer.data(lf + 248);
    const auto *lf_249 = buffer.data(lf + 249);
    const auto *lf_250 = buffer.data(lf + 250);
    const auto *lf_251 = buffer.data(lf + 251);
    const auto *lf_252 = buffer.data(lf + 252);
    const auto *lf_253 = buffer.data(lf + 253);
    const auto *lf_254 = buffer.data(lf + 254);
    const auto *lf_255 = buffer.data(lf + 255);
    const auto *lf_256 = buffer.data(lf + 256);
    const auto *lf_257 = buffer.data(lf + 257);
    const auto *lf_258 = buffer.data(lf + 258);
    const auto *lf_259 = buffer.data(lf + 259);
    const auto *lf_260 = buffer.data(lf + 260);
    const auto *lf_261 = buffer.data(lf + 261);
    const auto *lf_262 = buffer.data(lf + 262);
    const auto *lf_263 = buffer.data(lf + 263);
    const auto *lf_264 = buffer.data(lf + 264);
    const auto *lf_265 = buffer.data(lf + 265);
    const auto *lf_266 = buffer.data(lf + 266);
    const auto *lf_267 = buffer.data(lf + 267);
    const auto *lf_268 = buffer.data(lf + 268);
    const auto *lf_269 = buffer.data(lf + 269);
    const auto *lf_270 = buffer.data(lf + 270);
    const auto *lf_271 = buffer.data(lf + 271);
    const auto *lf_272 = buffer.data(lf + 272);
    const auto *lf_273 = buffer.data(lf + 273);
    const auto *lf_274 = buffer.data(lf + 274);
    const auto *lf_275 = buffer.data(lf + 275);
    const auto *lf_276 = buffer.data(lf + 276);
    const auto *lf_277 = buffer.data(lf + 277);
    const auto *lf_278 = buffer.data(lf + 278);
    const auto *lf_279 = buffer.data(lf + 279);
    const auto *lf_280 = buffer.data(lf + 280);
    const auto *lf_281 = buffer.data(lf + 281);
    const auto *lf_282 = buffer.data(lf + 282);
    const auto *lf_283 = buffer.data(lf + 283);
    const auto *lf_284 = buffer.data(lf + 284);
    const auto *lf_285 = buffer.data(lf + 285);
    const auto *lf_286 = buffer.data(lf + 286);
    const auto *lf_287 = buffer.data(lf + 287);
    const auto *lf_288 = buffer.data(lf + 288);
    const auto *lf_289 = buffer.data(lf + 289);
    const auto *lf_290 = buffer.data(lf + 290);
    const auto *lf_291 = buffer.data(lf + 291);
    const auto *lf_292 = buffer.data(lf + 292);
    const auto *lf_293 = buffer.data(lf + 293);
    const auto *lf_294 = buffer.data(lf + 294);
    const auto *lf_295 = buffer.data(lf + 295);
    const auto *lf_296 = buffer.data(lf + 296);
    const auto *lf_297 = buffer.data(lf + 297);
    const auto *lf_298 = buffer.data(lf + 298);
    const auto *lf_299 = buffer.data(lf + 299);
    const auto *lf_300 = buffer.data(lf + 300);
    const auto *lf_301 = buffer.data(lf + 301);
    const auto *lf_302 = buffer.data(lf + 302);
    const auto *lf_303 = buffer.data(lf + 303);
    const auto *lf_304 = buffer.data(lf + 304);
    const auto *lf_305 = buffer.data(lf + 305);
    const auto *lf_306 = buffer.data(lf + 306);
    const auto *lf_307 = buffer.data(lf + 307);
    const auto *lf_308 = buffer.data(lf + 308);
    const auto *lf_309 = buffer.data(lf + 309);
    const auto *lf_310 = buffer.data(lf + 310);
    const auto *lf_311 = buffer.data(lf + 311);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, if__150, if__151, if__152, \
                         if__153, if__154, lf_150, lf_151, lf_152, lf_153, \
                         lf_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -2.0 * if__150[k]
                   + f_0 * lf_150[k];

        t_151[k] = -2.0 * if__151[k]
                   + f_0 * lf_151[k];

        t_152[k] = -2.0 * if__152[k]
                   + f_0 * lf_152[k];

        t_153[k] = -2.0 * if__153[k]
                   + f_0 * lf_153[k];

        t_154[k] = -2.0 * if__154[k]
                   + f_0 * lf_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, if__155, if__156, if__157, \
                         if__158, if__159, lf_155, lf_156, lf_157, lf_158, \
                         lf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -2.0 * if__155[k]
                   + f_0 * lf_155[k];

        t_156[k] = -2.0 * if__156[k]
                   + f_0 * lf_156[k];

        t_157[k] = -2.0 * if__157[k]
                   + f_0 * lf_157[k];

        t_158[k] = -2.0 * if__158[k]
                   + f_0 * lf_158[k];

        t_159[k] = -2.0 * if__159[k]
                   + f_0 * lf_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, if__160, if__161, if__162, \
                         if__163, if__164, lf_160, lf_161, lf_162, lf_163, \
                         lf_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -2.0 * if__160[k]
                   + f_0 * lf_160[k];

        t_161[k] = -2.0 * if__161[k]
                   + f_0 * lf_161[k];

        t_162[k] = -2.0 * if__162[k]
                   + f_0 * lf_162[k];

        t_163[k] = -2.0 * if__163[k]
                   + f_0 * lf_163[k];

        t_164[k] = -2.0 * if__164[k]
                   + f_0 * lf_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, if__165, if__166, if__167, \
                         if__168, if__169, lf_165, lf_166, lf_167, lf_168, \
                         lf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -2.0 * if__165[k]
                   + f_0 * lf_165[k];

        t_166[k] = -2.0 * if__166[k]
                   + f_0 * lf_166[k];

        t_167[k] = -2.0 * if__167[k]
                   + f_0 * lf_167[k];

        t_168[k] = -2.0 * if__168[k]
                   + f_0 * lf_168[k];

        t_169[k] = -2.0 * if__169[k]
                   + f_0 * lf_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, if__170, if__171, if__172, \
                         if__173, if__174, lf_170, lf_171, lf_172, lf_173, \
                         lf_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -2.0 * if__170[k]
                   + f_0 * lf_170[k];

        t_171[k] = -2.0 * if__171[k]
                   + f_0 * lf_171[k];

        t_172[k] = -2.0 * if__172[k]
                   + f_0 * lf_172[k];

        t_173[k] = -2.0 * if__173[k]
                   + f_0 * lf_173[k];

        t_174[k] = -2.0 * if__174[k]
                   + f_0 * lf_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, if__175, if__176, if__177, \
                         if__178, if__179, lf_175, lf_176, lf_177, lf_178, \
                         lf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -2.0 * if__175[k]
                   + f_0 * lf_175[k];

        t_176[k] = -2.0 * if__176[k]
                   + f_0 * lf_176[k];

        t_177[k] = -2.0 * if__177[k]
                   + f_0 * lf_177[k];

        t_178[k] = -2.0 * if__178[k]
                   + f_0 * lf_178[k];

        t_179[k] = -2.0 * if__179[k]
                   + f_0 * lf_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, if__180, if__181, if__182, \
                         if__183, if__184, lf_180, lf_181, lf_182, lf_183, \
                         lf_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * if__180[k]
                   + f_0 * lf_180[k];

        t_181[k] = -2.0 * if__181[k]
                   + f_0 * lf_181[k];

        t_182[k] = -2.0 * if__182[k]
                   + f_0 * lf_182[k];

        t_183[k] = -2.0 * if__183[k]
                   + f_0 * lf_183[k];

        t_184[k] = -2.0 * if__184[k]
                   + f_0 * lf_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, if__185, if__186, if__187, \
                         if__188, if__189, lf_185, lf_186, lf_187, lf_188, \
                         lf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * if__185[k]
                   + f_0 * lf_185[k];

        t_186[k] = -2.0 * if__186[k]
                   + f_0 * lf_186[k];

        t_187[k] = -2.0 * if__187[k]
                   + f_0 * lf_187[k];

        t_188[k] = -2.0 * if__188[k]
                   + f_0 * lf_188[k];

        t_189[k] = -2.0 * if__189[k]
                   + f_0 * lf_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, if__190, if__191, if__192, \
                         if__193, if__194, lf_190, lf_191, lf_192, lf_193, \
                         lf_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * if__190[k]
                   + f_0 * lf_190[k];

        t_191[k] = -2.0 * if__191[k]
                   + f_0 * lf_191[k];

        t_192[k] = -2.0 * if__192[k]
                   + f_0 * lf_192[k];

        t_193[k] = -2.0 * if__193[k]
                   + f_0 * lf_193[k];

        t_194[k] = -2.0 * if__194[k]
                   + f_0 * lf_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, if__195, if__196, if__197, \
                         if__198, if__199, lf_195, lf_196, lf_197, lf_198, \
                         lf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * if__195[k]
                   + f_0 * lf_195[k];

        t_196[k] = -2.0 * if__196[k]
                   + f_0 * lf_196[k];

        t_197[k] = -2.0 * if__197[k]
                   + f_0 * lf_197[k];

        t_198[k] = -2.0 * if__198[k]
                   + f_0 * lf_198[k];

        t_199[k] = -2.0 * if__199[k]
                   + f_0 * lf_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, if__200, if__201, if__202, \
                         if__203, if__204, lf_200, lf_201, lf_202, lf_203, \
                         lf_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * if__200[k]
                   + f_0 * lf_200[k];

        t_201[k] = -2.0 * if__201[k]
                   + f_0 * lf_201[k];

        t_202[k] = -2.0 * if__202[k]
                   + f_0 * lf_202[k];

        t_203[k] = -2.0 * if__203[k]
                   + f_0 * lf_203[k];

        t_204[k] = -2.0 * if__204[k]
                   + f_0 * lf_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, if__205, if__206, if__207, \
                         if__208, if__209, lf_205, lf_206, lf_207, lf_208, \
                         lf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * if__205[k]
                   + f_0 * lf_205[k];

        t_206[k] = -2.0 * if__206[k]
                   + f_0 * lf_206[k];

        t_207[k] = -2.0 * if__207[k]
                   + f_0 * lf_207[k];

        t_208[k] = -2.0 * if__208[k]
                   + f_0 * lf_208[k];

        t_209[k] = -2.0 * if__209[k]
                   + f_0 * lf_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, if__210, if__211, if__212, \
                         if__213, if__214, lf_210, lf_211, lf_212, lf_213, \
                         lf_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -if__210[k]
                   + f_0 * lf_210[k];

        t_211[k] = -if__211[k]
                   + f_0 * lf_211[k];

        t_212[k] = -if__212[k]
                   + f_0 * lf_212[k];

        t_213[k] = -if__213[k]
                   + f_0 * lf_213[k];

        t_214[k] = -if__214[k]
                   + f_0 * lf_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, if__215, if__216, if__217, \
                         if__218, if__219, lf_215, lf_216, lf_217, lf_218, \
                         lf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -if__215[k]
                   + f_0 * lf_215[k];

        t_216[k] = -if__216[k]
                   + f_0 * lf_216[k];

        t_217[k] = -if__217[k]
                   + f_0 * lf_217[k];

        t_218[k] = -if__218[k]
                   + f_0 * lf_218[k];

        t_219[k] = -if__219[k]
                   + f_0 * lf_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, if__220, if__221, if__222, \
                         if__223, if__224, lf_220, lf_221, lf_222, lf_223, \
                         lf_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -if__220[k]
                   + f_0 * lf_220[k];

        t_221[k] = -if__221[k]
                   + f_0 * lf_221[k];

        t_222[k] = -if__222[k]
                   + f_0 * lf_222[k];

        t_223[k] = -if__223[k]
                   + f_0 * lf_223[k];

        t_224[k] = -if__224[k]
                   + f_0 * lf_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, if__225, if__226, if__227, \
                         if__228, if__229, lf_225, lf_226, lf_227, lf_228, \
                         lf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -if__225[k]
                   + f_0 * lf_225[k];

        t_226[k] = -if__226[k]
                   + f_0 * lf_226[k];

        t_227[k] = -if__227[k]
                   + f_0 * lf_227[k];

        t_228[k] = -if__228[k]
                   + f_0 * lf_228[k];

        t_229[k] = -if__229[k]
                   + f_0 * lf_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, if__230, if__231, if__232, \
                         if__233, if__234, lf_230, lf_231, lf_232, lf_233, \
                         lf_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -if__230[k]
                   + f_0 * lf_230[k];

        t_231[k] = -if__231[k]
                   + f_0 * lf_231[k];

        t_232[k] = -if__232[k]
                   + f_0 * lf_232[k];

        t_233[k] = -if__233[k]
                   + f_0 * lf_233[k];

        t_234[k] = -if__234[k]
                   + f_0 * lf_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, if__235, if__236, if__237, \
                         if__238, if__239, lf_235, lf_236, lf_237, lf_238, \
                         lf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -if__235[k]
                   + f_0 * lf_235[k];

        t_236[k] = -if__236[k]
                   + f_0 * lf_236[k];

        t_237[k] = -if__237[k]
                   + f_0 * lf_237[k];

        t_238[k] = -if__238[k]
                   + f_0 * lf_238[k];

        t_239[k] = -if__239[k]
                   + f_0 * lf_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, if__240, if__241, if__242, \
                         if__243, if__244, lf_240, lf_241, lf_242, lf_243, \
                         lf_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -if__240[k]
                   + f_0 * lf_240[k];

        t_241[k] = -if__241[k]
                   + f_0 * lf_241[k];

        t_242[k] = -if__242[k]
                   + f_0 * lf_242[k];

        t_243[k] = -if__243[k]
                   + f_0 * lf_243[k];

        t_244[k] = -if__244[k]
                   + f_0 * lf_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, if__245, if__246, if__247, \
                         if__248, if__249, lf_245, lf_246, lf_247, lf_248, \
                         lf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -if__245[k]
                   + f_0 * lf_245[k];

        t_246[k] = -if__246[k]
                   + f_0 * lf_246[k];

        t_247[k] = -if__247[k]
                   + f_0 * lf_247[k];

        t_248[k] = -if__248[k]
                   + f_0 * lf_248[k];

        t_249[k] = -if__249[k]
                   + f_0 * lf_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, if__250, if__251, if__252, \
                         if__253, if__254, lf_250, lf_251, lf_252, lf_253, \
                         lf_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -if__250[k]
                   + f_0 * lf_250[k];

        t_251[k] = -if__251[k]
                   + f_0 * lf_251[k];

        t_252[k] = -if__252[k]
                   + f_0 * lf_252[k];

        t_253[k] = -if__253[k]
                   + f_0 * lf_253[k];

        t_254[k] = -if__254[k]
                   + f_0 * lf_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, if__255, if__256, if__257, \
                         if__258, if__259, lf_255, lf_256, lf_257, lf_258, \
                         lf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -if__255[k]
                   + f_0 * lf_255[k];

        t_256[k] = -if__256[k]
                   + f_0 * lf_256[k];

        t_257[k] = -if__257[k]
                   + f_0 * lf_257[k];

        t_258[k] = -if__258[k]
                   + f_0 * lf_258[k];

        t_259[k] = -if__259[k]
                   + f_0 * lf_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, if__260, if__261, if__262, \
                         if__263, if__264, lf_260, lf_261, lf_262, lf_263, \
                         lf_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -if__260[k]
                   + f_0 * lf_260[k];

        t_261[k] = -if__261[k]
                   + f_0 * lf_261[k];

        t_262[k] = -if__262[k]
                   + f_0 * lf_262[k];

        t_263[k] = -if__263[k]
                   + f_0 * lf_263[k];

        t_264[k] = -if__264[k]
                   + f_0 * lf_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, if__265, if__266, if__267, \
                         if__268, if__269, lf_265, lf_266, lf_267, lf_268, \
                         lf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -if__265[k]
                   + f_0 * lf_265[k];

        t_266[k] = -if__266[k]
                   + f_0 * lf_266[k];

        t_267[k] = -if__267[k]
                   + f_0 * lf_267[k];

        t_268[k] = -if__268[k]
                   + f_0 * lf_268[k];

        t_269[k] = -if__269[k]
                   + f_0 * lf_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, if__270, if__271, if__272, \
                         if__273, if__274, lf_270, lf_271, lf_272, lf_273, \
                         lf_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -if__270[k]
                   + f_0 * lf_270[k];

        t_271[k] = -if__271[k]
                   + f_0 * lf_271[k];

        t_272[k] = -if__272[k]
                   + f_0 * lf_272[k];

        t_273[k] = -if__273[k]
                   + f_0 * lf_273[k];

        t_274[k] = -if__274[k]
                   + f_0 * lf_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, if__275, if__276, if__277, \
                         if__278, if__279, lf_275, lf_276, lf_277, lf_278, \
                         lf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -if__275[k]
                   + f_0 * lf_275[k];

        t_276[k] = -if__276[k]
                   + f_0 * lf_276[k];

        t_277[k] = -if__277[k]
                   + f_0 * lf_277[k];

        t_278[k] = -if__278[k]
                   + f_0 * lf_278[k];

        t_279[k] = -if__279[k]
                   + f_0 * lf_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, t_285, t_286, t_287, lf_280, \
                         lf_281, lf_282, lf_283, lf_284, lf_285, lf_286, \
                         lf_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_0 * lf_280[k];

        t_281[k] = f_0 * lf_281[k];

        t_282[k] = f_0 * lf_282[k];

        t_283[k] = f_0 * lf_283[k];

        t_284[k] = f_0 * lf_284[k];

        t_285[k] = f_0 * lf_285[k];

        t_286[k] = f_0 * lf_286[k];

        t_287[k] = f_0 * lf_287[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, t_293, t_294, t_295, lf_288, \
                         lf_289, lf_290, lf_291, lf_292, lf_293, lf_294, \
                         lf_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_0 * lf_288[k];

        t_289[k] = f_0 * lf_289[k];

        t_290[k] = f_0 * lf_290[k];

        t_291[k] = f_0 * lf_291[k];

        t_292[k] = f_0 * lf_292[k];

        t_293[k] = f_0 * lf_293[k];

        t_294[k] = f_0 * lf_294[k];

        t_295[k] = f_0 * lf_295[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, t_301, t_302, t_303, lf_296, \
                         lf_297, lf_298, lf_299, lf_300, lf_301, lf_302, \
                         lf_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_0 * lf_296[k];

        t_297[k] = f_0 * lf_297[k];

        t_298[k] = f_0 * lf_298[k];

        t_299[k] = f_0 * lf_299[k];

        t_300[k] = f_0 * lf_300[k];

        t_301[k] = f_0 * lf_301[k];

        t_302[k] = f_0 * lf_302[k];

        t_303[k] = f_0 * lf_303[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, t_309, t_310, t_311, lf_304, \
                         lf_305, lf_306, lf_307, lf_308, lf_309, lf_310, \
                         lf_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_0 * lf_304[k];

        t_305[k] = f_0 * lf_305[k];

        t_306[k] = f_0 * lf_306[k];

        t_307[k] = f_0 * lf_307[k];

        t_308[k] = f_0 * lf_308[k];

        t_309[k] = f_0 * lf_309[k];

        t_310[k] = f_0 * lf_310[k];

        t_311[k] = f_0 * lf_311[k];
    }
}

static auto
compute_prim_geom_10_kf_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t lf, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *lf_312 = buffer.data(lf + 312);
    const auto *lf_313 = buffer.data(lf + 313);
    const auto *lf_314 = buffer.data(lf + 314);
    const auto *lf_315 = buffer.data(lf + 315);
    const auto *lf_316 = buffer.data(lf + 316);
    const auto *lf_317 = buffer.data(lf + 317);
    const auto *lf_318 = buffer.data(lf + 318);
    const auto *lf_319 = buffer.data(lf + 319);
    const auto *lf_320 = buffer.data(lf + 320);
    const auto *lf_321 = buffer.data(lf + 321);
    const auto *lf_322 = buffer.data(lf + 322);
    const auto *lf_323 = buffer.data(lf + 323);
    const auto *lf_324 = buffer.data(lf + 324);
    const auto *lf_325 = buffer.data(lf + 325);
    const auto *lf_326 = buffer.data(lf + 326);
    const auto *lf_327 = buffer.data(lf + 327);
    const auto *lf_328 = buffer.data(lf + 328);
    const auto *lf_329 = buffer.data(lf + 329);
    const auto *lf_330 = buffer.data(lf + 330);
    const auto *lf_331 = buffer.data(lf + 331);
    const auto *lf_332 = buffer.data(lf + 332);
    const auto *lf_333 = buffer.data(lf + 333);
    const auto *lf_334 = buffer.data(lf + 334);
    const auto *lf_335 = buffer.data(lf + 335);
    const auto *lf_336 = buffer.data(lf + 336);
    const auto *lf_337 = buffer.data(lf + 337);
    const auto *lf_338 = buffer.data(lf + 338);
    const auto *lf_339 = buffer.data(lf + 339);
    const auto *lf_340 = buffer.data(lf + 340);
    const auto *lf_341 = buffer.data(lf + 341);
    const auto *lf_342 = buffer.data(lf + 342);
    const auto *lf_343 = buffer.data(lf + 343);
    const auto *lf_344 = buffer.data(lf + 344);
    const auto *lf_345 = buffer.data(lf + 345);
    const auto *lf_346 = buffer.data(lf + 346);
    const auto *lf_347 = buffer.data(lf + 347);
    const auto *lf_348 = buffer.data(lf + 348);
    const auto *lf_349 = buffer.data(lf + 349);
    const auto *lf_350 = buffer.data(lf + 350);
    const auto *lf_351 = buffer.data(lf + 351);
    const auto *lf_352 = buffer.data(lf + 352);
    const auto *lf_353 = buffer.data(lf + 353);
    const auto *lf_354 = buffer.data(lf + 354);
    const auto *lf_355 = buffer.data(lf + 355);
    const auto *lf_356 = buffer.data(lf + 356);
    const auto *lf_357 = buffer.data(lf + 357);
    const auto *lf_358 = buffer.data(lf + 358);
    const auto *lf_359 = buffer.data(lf + 359);

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, t_317, t_318, t_319, lf_312, \
                         lf_313, lf_314, lf_315, lf_316, lf_317, lf_318, \
                         lf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_0 * lf_312[k];

        t_313[k] = f_0 * lf_313[k];

        t_314[k] = f_0 * lf_314[k];

        t_315[k] = f_0 * lf_315[k];

        t_316[k] = f_0 * lf_316[k];

        t_317[k] = f_0 * lf_317[k];

        t_318[k] = f_0 * lf_318[k];

        t_319[k] = f_0 * lf_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, t_325, t_326, t_327, lf_320, \
                         lf_321, lf_322, lf_323, lf_324, lf_325, lf_326, \
                         lf_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_0 * lf_320[k];

        t_321[k] = f_0 * lf_321[k];

        t_322[k] = f_0 * lf_322[k];

        t_323[k] = f_0 * lf_323[k];

        t_324[k] = f_0 * lf_324[k];

        t_325[k] = f_0 * lf_325[k];

        t_326[k] = f_0 * lf_326[k];

        t_327[k] = f_0 * lf_327[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, t_333, t_334, t_335, lf_328, \
                         lf_329, lf_330, lf_331, lf_332, lf_333, lf_334, \
                         lf_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_0 * lf_328[k];

        t_329[k] = f_0 * lf_329[k];

        t_330[k] = f_0 * lf_330[k];

        t_331[k] = f_0 * lf_331[k];

        t_332[k] = f_0 * lf_332[k];

        t_333[k] = f_0 * lf_333[k];

        t_334[k] = f_0 * lf_334[k];

        t_335[k] = f_0 * lf_335[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, t_340, t_341, t_342, t_343, lf_336, \
                         lf_337, lf_338, lf_339, lf_340, lf_341, lf_342, \
                         lf_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_0 * lf_336[k];

        t_337[k] = f_0 * lf_337[k];

        t_338[k] = f_0 * lf_338[k];

        t_339[k] = f_0 * lf_339[k];

        t_340[k] = f_0 * lf_340[k];

        t_341[k] = f_0 * lf_341[k];

        t_342[k] = f_0 * lf_342[k];

        t_343[k] = f_0 * lf_343[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, t_349, t_350, t_351, lf_344, \
                         lf_345, lf_346, lf_347, lf_348, lf_349, lf_350, \
                         lf_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_0 * lf_344[k];

        t_345[k] = f_0 * lf_345[k];

        t_346[k] = f_0 * lf_346[k];

        t_347[k] = f_0 * lf_347[k];

        t_348[k] = f_0 * lf_348[k];

        t_349[k] = f_0 * lf_349[k];

        t_350[k] = f_0 * lf_350[k];

        t_351[k] = f_0 * lf_351[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, t_357, t_358, t_359, lf_352, \
                         lf_353, lf_354, lf_355, lf_356, lf_357, lf_358, \
                         lf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_0 * lf_352[k];

        t_353[k] = f_0 * lf_353[k];

        t_354[k] = f_0 * lf_354[k];

        t_355[k] = f_0 * lf_355[k];

        t_356[k] = f_0 * lf_356[k];

        t_357[k] = f_0 * lf_357[k];

        t_358[k] = f_0 * lf_358[k];

        t_359[k] = f_0 * lf_359[k];
    }
}

auto
compute_prim_geom_10_kf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t if_, const size_t lf,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kf_electron_repulsion_0_piece0(buffer, target, if_, lf, ncols, alpha);

    compute_prim_geom_10_kf_electron_repulsion_0_piece1(buffer, target, if_, lf, ncols, alpha);

    compute_prim_geom_10_kf_electron_repulsion_0_piece2(buffer, target, lf, ncols, alpha);
}

static auto
compute_prim_geom_10_kf_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t if_, const size_t lf,
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

    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);
    const auto *lf_108 = buffer.data(lf + 108);
    const auto *lf_109 = buffer.data(lf + 109);
    const auto *lf_110 = buffer.data(lf + 110);
    const auto *lf_111 = buffer.data(lf + 111);
    const auto *lf_112 = buffer.data(lf + 112);
    const auto *lf_113 = buffer.data(lf + 113);
    const auto *lf_114 = buffer.data(lf + 114);
    const auto *lf_115 = buffer.data(lf + 115);
    const auto *lf_116 = buffer.data(lf + 116);
    const auto *lf_117 = buffer.data(lf + 117);
    const auto *lf_118 = buffer.data(lf + 118);
    const auto *lf_119 = buffer.data(lf + 119);
    const auto *lf_120 = buffer.data(lf + 120);
    const auto *lf_121 = buffer.data(lf + 121);
    const auto *lf_122 = buffer.data(lf + 122);
    const auto *lf_123 = buffer.data(lf + 123);
    const auto *lf_124 = buffer.data(lf + 124);
    const auto *lf_125 = buffer.data(lf + 125);
    const auto *lf_126 = buffer.data(lf + 126);
    const auto *lf_127 = buffer.data(lf + 127);
    const auto *lf_128 = buffer.data(lf + 128);
    const auto *lf_129 = buffer.data(lf + 129);
    const auto *lf_130 = buffer.data(lf + 130);
    const auto *lf_131 = buffer.data(lf + 131);
    const auto *lf_132 = buffer.data(lf + 132);
    const auto *lf_133 = buffer.data(lf + 133);
    const auto *lf_134 = buffer.data(lf + 134);
    const auto *lf_135 = buffer.data(lf + 135);
    const auto *lf_136 = buffer.data(lf + 136);
    const auto *lf_137 = buffer.data(lf + 137);
    const auto *lf_138 = buffer.data(lf + 138);
    const auto *lf_139 = buffer.data(lf + 139);
    const auto *lf_150 = buffer.data(lf + 150);
    const auto *lf_151 = buffer.data(lf + 151);
    const auto *lf_152 = buffer.data(lf + 152);
    const auto *lf_153 = buffer.data(lf + 153);
    const auto *lf_154 = buffer.data(lf + 154);
    const auto *lf_155 = buffer.data(lf + 155);
    const auto *lf_156 = buffer.data(lf + 156);
    const auto *lf_157 = buffer.data(lf + 157);
    const auto *lf_158 = buffer.data(lf + 158);
    const auto *lf_159 = buffer.data(lf + 159);
    const auto *lf_160 = buffer.data(lf + 160);
    const auto *lf_161 = buffer.data(lf + 161);
    const auto *lf_162 = buffer.data(lf + 162);
    const auto *lf_163 = buffer.data(lf + 163);
    const auto *lf_164 = buffer.data(lf + 164);
    const auto *lf_165 = buffer.data(lf + 165);
    const auto *lf_166 = buffer.data(lf + 166);
    const auto *lf_167 = buffer.data(lf + 167);
    const auto *lf_168 = buffer.data(lf + 168);
    const auto *lf_169 = buffer.data(lf + 169);
    const auto *lf_170 = buffer.data(lf + 170);
    const auto *lf_171 = buffer.data(lf + 171);
    const auto *lf_172 = buffer.data(lf + 172);
    const auto *lf_173 = buffer.data(lf + 173);
    const auto *lf_174 = buffer.data(lf + 174);
    const auto *lf_175 = buffer.data(lf + 175);
    const auto *lf_176 = buffer.data(lf + 176);
    const auto *lf_177 = buffer.data(lf + 177);
    const auto *lf_178 = buffer.data(lf + 178);
    const auto *lf_179 = buffer.data(lf + 179);
    const auto *lf_180 = buffer.data(lf + 180);
    const auto *lf_181 = buffer.data(lf + 181);
    const auto *lf_182 = buffer.data(lf + 182);
    const auto *lf_183 = buffer.data(lf + 183);
    const auto *lf_184 = buffer.data(lf + 184);
    const auto *lf_185 = buffer.data(lf + 185);
    const auto *lf_186 = buffer.data(lf + 186);
    const auto *lf_187 = buffer.data(lf + 187);
    const auto *lf_188 = buffer.data(lf + 188);
    const auto *lf_189 = buffer.data(lf + 189);
    const auto *lf_190 = buffer.data(lf + 190);
    const auto *lf_191 = buffer.data(lf + 191);
    const auto *lf_192 = buffer.data(lf + 192);
    const auto *lf_193 = buffer.data(lf + 193);
    const auto *lf_194 = buffer.data(lf + 194);
    const auto *lf_195 = buffer.data(lf + 195);
    const auto *lf_196 = buffer.data(lf + 196);
    const auto *lf_197 = buffer.data(lf + 197);
    const auto *lf_198 = buffer.data(lf + 198);
    const auto *lf_199 = buffer.data(lf + 199);
    const auto *lf_210 = buffer.data(lf + 210);
    const auto *lf_211 = buffer.data(lf + 211);
    const auto *lf_212 = buffer.data(lf + 212);
    const auto *lf_213 = buffer.data(lf + 213);
    const auto *lf_214 = buffer.data(lf + 214);
    const auto *lf_215 = buffer.data(lf + 215);
    const auto *lf_216 = buffer.data(lf + 216);
    const auto *lf_217 = buffer.data(lf + 217);
    const auto *lf_218 = buffer.data(lf + 218);
    const auto *lf_219 = buffer.data(lf + 219);
    const auto *lf_220 = buffer.data(lf + 220);
    const auto *lf_221 = buffer.data(lf + 221);
    const auto *lf_222 = buffer.data(lf + 222);
    const auto *lf_223 = buffer.data(lf + 223);
    const auto *lf_224 = buffer.data(lf + 224);
    const auto *lf_225 = buffer.data(lf + 225);
    const auto *lf_226 = buffer.data(lf + 226);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, lf_10, lf_11, lf_12, lf_13, \
                         lf_14, lf_15, lf_16, lf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lf_10[k];

        t_1[k] = f_0 * lf_11[k];

        t_2[k] = f_0 * lf_12[k];

        t_3[k] = f_0 * lf_13[k];

        t_4[k] = f_0 * lf_14[k];

        t_5[k] = f_0 * lf_15[k];

        t_6[k] = f_0 * lf_16[k];

        t_7[k] = f_0 * lf_17[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, if__0, if__1, if__2, if__3, lf_18, \
                         lf_19, lf_30, lf_31, lf_32, lf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * lf_18[k];

        t_9[k] = f_0 * lf_19[k];

        t_10[k] = -if__0[k]
                  + f_0 * lf_30[k];

        t_11[k] = -if__1[k]
                  + f_0 * lf_31[k];

        t_12[k] = -if__2[k]
                  + f_0 * lf_32[k];

        t_13[k] = -if__3[k]
                  + f_0 * lf_33[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, if__4, if__5, if__6, if__7, if__8, \
                         lf_34, lf_35, lf_36, lf_37, lf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -if__4[k]
                  + f_0 * lf_34[k];

        t_15[k] = -if__5[k]
                  + f_0 * lf_35[k];

        t_16[k] = -if__6[k]
                  + f_0 * lf_36[k];

        t_17[k] = -if__7[k]
                  + f_0 * lf_37[k];

        t_18[k] = -if__8[k]
                  + f_0 * lf_38[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, t_25, if__9, lf_39, lf_40, lf_41, \
                         lf_42, lf_43, lf_44, lf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -if__9[k]
                  + f_0 * lf_39[k];

        t_20[k] = f_0 * lf_40[k];

        t_21[k] = f_0 * lf_41[k];

        t_22[k] = f_0 * lf_42[k];

        t_23[k] = f_0 * lf_43[k];

        t_24[k] = f_0 * lf_44[k];

        t_25[k] = f_0 * lf_45[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, if__10, if__11, lf_46, lf_47, \
                         lf_48, lf_49, lf_60, lf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * lf_46[k];

        t_27[k] = f_0 * lf_47[k];

        t_28[k] = f_0 * lf_48[k];

        t_29[k] = f_0 * lf_49[k];

        t_30[k] = -2.0 * if__10[k]
                  + f_0 * lf_60[k];

        t_31[k] = -2.0 * if__11[k]
                  + f_0 * lf_61[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, if__12, if__13, if__14, if__15, if__16, \
                         lf_62, lf_63, lf_64, lf_65, lf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -2.0 * if__12[k]
                  + f_0 * lf_62[k];

        t_33[k] = -2.0 * if__13[k]
                  + f_0 * lf_63[k];

        t_34[k] = -2.0 * if__14[k]
                  + f_0 * lf_64[k];

        t_35[k] = -2.0 * if__15[k]
                  + f_0 * lf_65[k];

        t_36[k] = -2.0 * if__16[k]
                  + f_0 * lf_66[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, if__17, if__18, if__19, if__20, if__21, \
                         lf_67, lf_68, lf_69, lf_70, lf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -2.0 * if__17[k]
                  + f_0 * lf_67[k];

        t_38[k] = -2.0 * if__18[k]
                  + f_0 * lf_68[k];

        t_39[k] = -2.0 * if__19[k]
                  + f_0 * lf_69[k];

        t_40[k] = -if__20[k]
                  + f_0 * lf_70[k];

        t_41[k] = -if__21[k]
                  + f_0 * lf_71[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, if__22, if__23, if__24, if__25, if__26, \
                         lf_72, lf_73, lf_74, lf_75, lf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -if__22[k]
                  + f_0 * lf_72[k];

        t_43[k] = -if__23[k]
                  + f_0 * lf_73[k];

        t_44[k] = -if__24[k]
                  + f_0 * lf_74[k];

        t_45[k] = -if__25[k]
                  + f_0 * lf_75[k];

        t_46[k] = -if__26[k]
                  + f_0 * lf_76[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, if__27, if__28, if__29, lf_77, \
                         lf_78, lf_79, lf_80, lf_81, lf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -if__27[k]
                  + f_0 * lf_77[k];

        t_48[k] = -if__28[k]
                  + f_0 * lf_78[k];

        t_49[k] = -if__29[k]
                  + f_0 * lf_79[k];

        t_50[k] = f_0 * lf_80[k];

        t_51[k] = f_0 * lf_81[k];

        t_52[k] = f_0 * lf_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, t_59, lf_83, lf_84, lf_85, lf_86, \
                         lf_87, lf_88, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * lf_83[k];

        t_54[k] = f_0 * lf_84[k];

        t_55[k] = f_0 * lf_85[k];

        t_56[k] = f_0 * lf_86[k];

        t_57[k] = f_0 * lf_87[k];

        t_58[k] = f_0 * lf_88[k];

        t_59[k] = f_0 * lf_89[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, if__30, if__31, if__32, if__33, if__34, \
                         lf_100, lf_101, lf_102, lf_103, lf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * if__30[k]
                  + f_0 * lf_100[k];

        t_61[k] = -3.0 * if__31[k]
                  + f_0 * lf_101[k];

        t_62[k] = -3.0 * if__32[k]
                  + f_0 * lf_102[k];

        t_63[k] = -3.0 * if__33[k]
                  + f_0 * lf_103[k];

        t_64[k] = -3.0 * if__34[k]
                  + f_0 * lf_104[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, if__35, if__36, if__37, if__38, if__39, \
                         lf_105, lf_106, lf_107, lf_108, lf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * if__35[k]
                  + f_0 * lf_105[k];

        t_66[k] = -3.0 * if__36[k]
                  + f_0 * lf_106[k];

        t_67[k] = -3.0 * if__37[k]
                  + f_0 * lf_107[k];

        t_68[k] = -3.0 * if__38[k]
                  + f_0 * lf_108[k];

        t_69[k] = -3.0 * if__39[k]
                  + f_0 * lf_109[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, if__40, if__41, if__42, if__43, if__44, \
                         lf_110, lf_111, lf_112, lf_113, lf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * if__40[k]
                  + f_0 * lf_110[k];

        t_71[k] = -2.0 * if__41[k]
                  + f_0 * lf_111[k];

        t_72[k] = -2.0 * if__42[k]
                  + f_0 * lf_112[k];

        t_73[k] = -2.0 * if__43[k]
                  + f_0 * lf_113[k];

        t_74[k] = -2.0 * if__44[k]
                  + f_0 * lf_114[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, if__45, if__46, if__47, if__48, if__49, \
                         lf_115, lf_116, lf_117, lf_118, lf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * if__45[k]
                  + f_0 * lf_115[k];

        t_76[k] = -2.0 * if__46[k]
                  + f_0 * lf_116[k];

        t_77[k] = -2.0 * if__47[k]
                  + f_0 * lf_117[k];

        t_78[k] = -2.0 * if__48[k]
                  + f_0 * lf_118[k];

        t_79[k] = -2.0 * if__49[k]
                  + f_0 * lf_119[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, if__50, if__51, if__52, if__53, if__54, \
                         lf_120, lf_121, lf_122, lf_123, lf_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -if__50[k]
                  + f_0 * lf_120[k];

        t_81[k] = -if__51[k]
                  + f_0 * lf_121[k];

        t_82[k] = -if__52[k]
                  + f_0 * lf_122[k];

        t_83[k] = -if__53[k]
                  + f_0 * lf_123[k];

        t_84[k] = -if__54[k]
                  + f_0 * lf_124[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, if__55, if__56, if__57, if__58, if__59, \
                         lf_125, lf_126, lf_127, lf_128, lf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -if__55[k]
                  + f_0 * lf_125[k];

        t_86[k] = -if__56[k]
                  + f_0 * lf_126[k];

        t_87[k] = -if__57[k]
                  + f_0 * lf_127[k];

        t_88[k] = -if__58[k]
                  + f_0 * lf_128[k];

        t_89[k] = -if__59[k]
                  + f_0 * lf_129[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, lf_130, lf_131, \
                         lf_132, lf_133, lf_134, lf_135, lf_136, \
                         lf_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * lf_130[k];

        t_91[k] = f_0 * lf_131[k];

        t_92[k] = f_0 * lf_132[k];

        t_93[k] = f_0 * lf_133[k];

        t_94[k] = f_0 * lf_134[k];

        t_95[k] = f_0 * lf_135[k];

        t_96[k] = f_0 * lf_136[k];

        t_97[k] = f_0 * lf_137[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, if__60, if__61, if__62, \
                         if__63, lf_138, lf_139, lf_150, lf_151, lf_152, \
                         lf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * lf_138[k];

        t_99[k] = f_0 * lf_139[k];

        t_100[k] = -4.0 * if__60[k]
                   + f_0 * lf_150[k];

        t_101[k] = -4.0 * if__61[k]
                   + f_0 * lf_151[k];

        t_102[k] = -4.0 * if__62[k]
                   + f_0 * lf_152[k];

        t_103[k] = -4.0 * if__63[k]
                   + f_0 * lf_153[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, if__64, if__65, if__66, if__67, \
                         if__68, lf_154, lf_155, lf_156, lf_157, \
                         lf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -4.0 * if__64[k]
                   + f_0 * lf_154[k];

        t_105[k] = -4.0 * if__65[k]
                   + f_0 * lf_155[k];

        t_106[k] = -4.0 * if__66[k]
                   + f_0 * lf_156[k];

        t_107[k] = -4.0 * if__67[k]
                   + f_0 * lf_157[k];

        t_108[k] = -4.0 * if__68[k]
                   + f_0 * lf_158[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, if__69, if__70, if__71, if__72, \
                         if__73, lf_159, lf_160, lf_161, lf_162, \
                         lf_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -4.0 * if__69[k]
                   + f_0 * lf_159[k];

        t_110[k] = -3.0 * if__70[k]
                   + f_0 * lf_160[k];

        t_111[k] = -3.0 * if__71[k]
                   + f_0 * lf_161[k];

        t_112[k] = -3.0 * if__72[k]
                   + f_0 * lf_162[k];

        t_113[k] = -3.0 * if__73[k]
                   + f_0 * lf_163[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, if__74, if__75, if__76, if__77, \
                         if__78, lf_164, lf_165, lf_166, lf_167, \
                         lf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -3.0 * if__74[k]
                   + f_0 * lf_164[k];

        t_115[k] = -3.0 * if__75[k]
                   + f_0 * lf_165[k];

        t_116[k] = -3.0 * if__76[k]
                   + f_0 * lf_166[k];

        t_117[k] = -3.0 * if__77[k]
                   + f_0 * lf_167[k];

        t_118[k] = -3.0 * if__78[k]
                   + f_0 * lf_168[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, if__79, if__80, if__81, if__82, \
                         if__83, lf_169, lf_170, lf_171, lf_172, \
                         lf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -3.0 * if__79[k]
                   + f_0 * lf_169[k];

        t_120[k] = -2.0 * if__80[k]
                   + f_0 * lf_170[k];

        t_121[k] = -2.0 * if__81[k]
                   + f_0 * lf_171[k];

        t_122[k] = -2.0 * if__82[k]
                   + f_0 * lf_172[k];

        t_123[k] = -2.0 * if__83[k]
                   + f_0 * lf_173[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, if__84, if__85, if__86, if__87, \
                         if__88, lf_174, lf_175, lf_176, lf_177, \
                         lf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * if__84[k]
                   + f_0 * lf_174[k];

        t_125[k] = -2.0 * if__85[k]
                   + f_0 * lf_175[k];

        t_126[k] = -2.0 * if__86[k]
                   + f_0 * lf_176[k];

        t_127[k] = -2.0 * if__87[k]
                   + f_0 * lf_177[k];

        t_128[k] = -2.0 * if__88[k]
                   + f_0 * lf_178[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, if__89, if__90, if__91, if__92, \
                         if__93, lf_179, lf_180, lf_181, lf_182, \
                         lf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -2.0 * if__89[k]
                   + f_0 * lf_179[k];

        t_130[k] = -if__90[k]
                   + f_0 * lf_180[k];

        t_131[k] = -if__91[k]
                   + f_0 * lf_181[k];

        t_132[k] = -if__92[k]
                   + f_0 * lf_182[k];

        t_133[k] = -if__93[k]
                   + f_0 * lf_183[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, if__94, if__95, if__96, if__97, \
                         if__98, lf_184, lf_185, lf_186, lf_187, \
                         lf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -if__94[k]
                   + f_0 * lf_184[k];

        t_135[k] = -if__95[k]
                   + f_0 * lf_185[k];

        t_136[k] = -if__96[k]
                   + f_0 * lf_186[k];

        t_137[k] = -if__97[k]
                   + f_0 * lf_187[k];

        t_138[k] = -if__98[k]
                   + f_0 * lf_188[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, t_145, if__99, lf_189, \
                         lf_190, lf_191, lf_192, lf_193, lf_194, \
                         lf_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -if__99[k]
                   + f_0 * lf_189[k];

        t_140[k] = f_0 * lf_190[k];

        t_141[k] = f_0 * lf_191[k];

        t_142[k] = f_0 * lf_192[k];

        t_143[k] = f_0 * lf_193[k];

        t_144[k] = f_0 * lf_194[k];

        t_145[k] = f_0 * lf_195[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, t_151, if__100, if__101, lf_196, \
                         lf_197, lf_198, lf_199, lf_210, lf_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_0 * lf_196[k];

        t_147[k] = f_0 * lf_197[k];

        t_148[k] = f_0 * lf_198[k];

        t_149[k] = f_0 * lf_199[k];

        t_150[k] = -5.0 * if__100[k]
                   + f_0 * lf_210[k];

        t_151[k] = -5.0 * if__101[k]
                   + f_0 * lf_211[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, if__102, if__103, if__104, \
                         if__105, if__106, lf_212, lf_213, lf_214, lf_215, \
                         lf_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -5.0 * if__102[k]
                   + f_0 * lf_212[k];

        t_153[k] = -5.0 * if__103[k]
                   + f_0 * lf_213[k];

        t_154[k] = -5.0 * if__104[k]
                   + f_0 * lf_214[k];

        t_155[k] = -5.0 * if__105[k]
                   + f_0 * lf_215[k];

        t_156[k] = -5.0 * if__106[k]
                   + f_0 * lf_216[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, if__107, if__108, if__109, \
                         if__110, if__111, lf_217, lf_218, lf_219, lf_220, \
                         lf_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -5.0 * if__107[k]
                   + f_0 * lf_217[k];

        t_158[k] = -5.0 * if__108[k]
                   + f_0 * lf_218[k];

        t_159[k] = -5.0 * if__109[k]
                   + f_0 * lf_219[k];

        t_160[k] = -4.0 * if__110[k]
                   + f_0 * lf_220[k];

        t_161[k] = -4.0 * if__111[k]
                   + f_0 * lf_221[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, if__112, if__113, if__114, \
                         if__115, if__116, lf_222, lf_223, lf_224, lf_225, \
                         lf_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -4.0 * if__112[k]
                   + f_0 * lf_222[k];

        t_163[k] = -4.0 * if__113[k]
                   + f_0 * lf_223[k];

        t_164[k] = -4.0 * if__114[k]
                   + f_0 * lf_224[k];

        t_165[k] = -4.0 * if__115[k]
                   + f_0 * lf_225[k];

        t_166[k] = -4.0 * if__116[k]
                   + f_0 * lf_226[k];
    }
}

static auto
compute_prim_geom_10_kf_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t if_, const size_t lf,
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

    const auto *lf_227 = buffer.data(lf + 227);
    const auto *lf_228 = buffer.data(lf + 228);
    const auto *lf_229 = buffer.data(lf + 229);
    const auto *lf_230 = buffer.data(lf + 230);
    const auto *lf_231 = buffer.data(lf + 231);
    const auto *lf_232 = buffer.data(lf + 232);
    const auto *lf_233 = buffer.data(lf + 233);
    const auto *lf_234 = buffer.data(lf + 234);
    const auto *lf_235 = buffer.data(lf + 235);
    const auto *lf_236 = buffer.data(lf + 236);
    const auto *lf_237 = buffer.data(lf + 237);
    const auto *lf_238 = buffer.data(lf + 238);
    const auto *lf_239 = buffer.data(lf + 239);
    const auto *lf_240 = buffer.data(lf + 240);
    const auto *lf_241 = buffer.data(lf + 241);
    const auto *lf_242 = buffer.data(lf + 242);
    const auto *lf_243 = buffer.data(lf + 243);
    const auto *lf_244 = buffer.data(lf + 244);
    const auto *lf_245 = buffer.data(lf + 245);
    const auto *lf_246 = buffer.data(lf + 246);
    const auto *lf_247 = buffer.data(lf + 247);
    const auto *lf_248 = buffer.data(lf + 248);
    const auto *lf_249 = buffer.data(lf + 249);
    const auto *lf_250 = buffer.data(lf + 250);
    const auto *lf_251 = buffer.data(lf + 251);
    const auto *lf_252 = buffer.data(lf + 252);
    const auto *lf_253 = buffer.data(lf + 253);
    const auto *lf_254 = buffer.data(lf + 254);
    const auto *lf_255 = buffer.data(lf + 255);
    const auto *lf_256 = buffer.data(lf + 256);
    const auto *lf_257 = buffer.data(lf + 257);
    const auto *lf_258 = buffer.data(lf + 258);
    const auto *lf_259 = buffer.data(lf + 259);
    const auto *lf_260 = buffer.data(lf + 260);
    const auto *lf_261 = buffer.data(lf + 261);
    const auto *lf_262 = buffer.data(lf + 262);
    const auto *lf_263 = buffer.data(lf + 263);
    const auto *lf_264 = buffer.data(lf + 264);
    const auto *lf_265 = buffer.data(lf + 265);
    const auto *lf_266 = buffer.data(lf + 266);
    const auto *lf_267 = buffer.data(lf + 267);
    const auto *lf_268 = buffer.data(lf + 268);
    const auto *lf_269 = buffer.data(lf + 269);
    const auto *lf_280 = buffer.data(lf + 280);
    const auto *lf_281 = buffer.data(lf + 281);
    const auto *lf_282 = buffer.data(lf + 282);
    const auto *lf_283 = buffer.data(lf + 283);
    const auto *lf_284 = buffer.data(lf + 284);
    const auto *lf_285 = buffer.data(lf + 285);
    const auto *lf_286 = buffer.data(lf + 286);
    const auto *lf_287 = buffer.data(lf + 287);
    const auto *lf_288 = buffer.data(lf + 288);
    const auto *lf_289 = buffer.data(lf + 289);
    const auto *lf_290 = buffer.data(lf + 290);
    const auto *lf_291 = buffer.data(lf + 291);
    const auto *lf_292 = buffer.data(lf + 292);
    const auto *lf_293 = buffer.data(lf + 293);
    const auto *lf_294 = buffer.data(lf + 294);
    const auto *lf_295 = buffer.data(lf + 295);
    const auto *lf_296 = buffer.data(lf + 296);
    const auto *lf_297 = buffer.data(lf + 297);
    const auto *lf_298 = buffer.data(lf + 298);
    const auto *lf_299 = buffer.data(lf + 299);
    const auto *lf_300 = buffer.data(lf + 300);
    const auto *lf_301 = buffer.data(lf + 301);
    const auto *lf_302 = buffer.data(lf + 302);
    const auto *lf_303 = buffer.data(lf + 303);
    const auto *lf_304 = buffer.data(lf + 304);
    const auto *lf_305 = buffer.data(lf + 305);
    const auto *lf_306 = buffer.data(lf + 306);
    const auto *lf_307 = buffer.data(lf + 307);
    const auto *lf_308 = buffer.data(lf + 308);
    const auto *lf_309 = buffer.data(lf + 309);
    const auto *lf_310 = buffer.data(lf + 310);
    const auto *lf_311 = buffer.data(lf + 311);
    const auto *lf_312 = buffer.data(lf + 312);
    const auto *lf_313 = buffer.data(lf + 313);
    const auto *lf_314 = buffer.data(lf + 314);
    const auto *lf_315 = buffer.data(lf + 315);
    const auto *lf_316 = buffer.data(lf + 316);
    const auto *lf_317 = buffer.data(lf + 317);
    const auto *lf_318 = buffer.data(lf + 318);
    const auto *lf_319 = buffer.data(lf + 319);
    const auto *lf_320 = buffer.data(lf + 320);
    const auto *lf_321 = buffer.data(lf + 321);
    const auto *lf_322 = buffer.data(lf + 322);
    const auto *lf_323 = buffer.data(lf + 323);
    const auto *lf_324 = buffer.data(lf + 324);
    const auto *lf_325 = buffer.data(lf + 325);
    const auto *lf_326 = buffer.data(lf + 326);
    const auto *lf_327 = buffer.data(lf + 327);
    const auto *lf_328 = buffer.data(lf + 328);
    const auto *lf_329 = buffer.data(lf + 329);
    const auto *lf_330 = buffer.data(lf + 330);
    const auto *lf_331 = buffer.data(lf + 331);
    const auto *lf_332 = buffer.data(lf + 332);
    const auto *lf_333 = buffer.data(lf + 333);
    const auto *lf_334 = buffer.data(lf + 334);
    const auto *lf_335 = buffer.data(lf + 335);
    const auto *lf_336 = buffer.data(lf + 336);
    const auto *lf_337 = buffer.data(lf + 337);
    const auto *lf_338 = buffer.data(lf + 338);
    const auto *lf_339 = buffer.data(lf + 339);
    const auto *lf_340 = buffer.data(lf + 340);
    const auto *lf_341 = buffer.data(lf + 341);
    const auto *lf_342 = buffer.data(lf + 342);
    const auto *lf_343 = buffer.data(lf + 343);
    const auto *lf_344 = buffer.data(lf + 344);
    const auto *lf_345 = buffer.data(lf + 345);
    const auto *lf_346 = buffer.data(lf + 346);
    const auto *lf_347 = buffer.data(lf + 347);
    const auto *lf_348 = buffer.data(lf + 348);
    const auto *lf_349 = buffer.data(lf + 349);
    const auto *lf_360 = buffer.data(lf + 360);
    const auto *lf_361 = buffer.data(lf + 361);
    const auto *lf_362 = buffer.data(lf + 362);
    const auto *lf_363 = buffer.data(lf + 363);
    const auto *lf_364 = buffer.data(lf + 364);
    const auto *lf_365 = buffer.data(lf + 365);
    const auto *lf_366 = buffer.data(lf + 366);
    const auto *lf_367 = buffer.data(lf + 367);
    const auto *lf_368 = buffer.data(lf + 368);
    const auto *lf_369 = buffer.data(lf + 369);
    const auto *lf_370 = buffer.data(lf + 370);
    const auto *lf_371 = buffer.data(lf + 371);
    const auto *lf_372 = buffer.data(lf + 372);
    const auto *lf_373 = buffer.data(lf + 373);
    const auto *lf_374 = buffer.data(lf + 374);
    const auto *lf_375 = buffer.data(lf + 375);
    const auto *lf_376 = buffer.data(lf + 376);
    const auto *lf_377 = buffer.data(lf + 377);
    const auto *lf_378 = buffer.data(lf + 378);
    const auto *lf_379 = buffer.data(lf + 379);
    const auto *lf_380 = buffer.data(lf + 380);
    const auto *lf_381 = buffer.data(lf + 381);
    const auto *lf_382 = buffer.data(lf + 382);
    const auto *lf_383 = buffer.data(lf + 383);
    const auto *lf_384 = buffer.data(lf + 384);
    const auto *lf_385 = buffer.data(lf + 385);
    const auto *lf_386 = buffer.data(lf + 386);
    const auto *lf_387 = buffer.data(lf + 387);
    const auto *lf_388 = buffer.data(lf + 388);
    const auto *lf_389 = buffer.data(lf + 389);
    const auto *lf_390 = buffer.data(lf + 390);
    const auto *lf_391 = buffer.data(lf + 391);
    const auto *lf_392 = buffer.data(lf + 392);
    const auto *lf_393 = buffer.data(lf + 393);
    const auto *lf_394 = buffer.data(lf + 394);
    const auto *lf_395 = buffer.data(lf + 395);
    const auto *lf_396 = buffer.data(lf + 396);
    const auto *lf_397 = buffer.data(lf + 397);
    const auto *lf_398 = buffer.data(lf + 398);
    const auto *lf_399 = buffer.data(lf + 399);
    const auto *lf_400 = buffer.data(lf + 400);
    const auto *lf_401 = buffer.data(lf + 401);
    const auto *lf_402 = buffer.data(lf + 402);
    const auto *lf_403 = buffer.data(lf + 403);

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, if__117, if__118, if__119, \
                         if__120, if__121, lf_227, lf_228, lf_229, lf_230, \
                         lf_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -4.0 * if__117[k]
                   + f_0 * lf_227[k];

        t_168[k] = -4.0 * if__118[k]
                   + f_0 * lf_228[k];

        t_169[k] = -4.0 * if__119[k]
                   + f_0 * lf_229[k];

        t_170[k] = -3.0 * if__120[k]
                   + f_0 * lf_230[k];

        t_171[k] = -3.0 * if__121[k]
                   + f_0 * lf_231[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, if__122, if__123, if__124, \
                         if__125, if__126, lf_232, lf_233, lf_234, lf_235, \
                         lf_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -3.0 * if__122[k]
                   + f_0 * lf_232[k];

        t_173[k] = -3.0 * if__123[k]
                   + f_0 * lf_233[k];

        t_174[k] = -3.0 * if__124[k]
                   + f_0 * lf_234[k];

        t_175[k] = -3.0 * if__125[k]
                   + f_0 * lf_235[k];

        t_176[k] = -3.0 * if__126[k]
                   + f_0 * lf_236[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, if__127, if__128, if__129, \
                         if__130, if__131, lf_237, lf_238, lf_239, lf_240, \
                         lf_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -3.0 * if__127[k]
                   + f_0 * lf_237[k];

        t_178[k] = -3.0 * if__128[k]
                   + f_0 * lf_238[k];

        t_179[k] = -3.0 * if__129[k]
                   + f_0 * lf_239[k];

        t_180[k] = -2.0 * if__130[k]
                   + f_0 * lf_240[k];

        t_181[k] = -2.0 * if__131[k]
                   + f_0 * lf_241[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, if__132, if__133, if__134, \
                         if__135, if__136, lf_242, lf_243, lf_244, lf_245, \
                         lf_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * if__132[k]
                   + f_0 * lf_242[k];

        t_183[k] = -2.0 * if__133[k]
                   + f_0 * lf_243[k];

        t_184[k] = -2.0 * if__134[k]
                   + f_0 * lf_244[k];

        t_185[k] = -2.0 * if__135[k]
                   + f_0 * lf_245[k];

        t_186[k] = -2.0 * if__136[k]
                   + f_0 * lf_246[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, if__137, if__138, if__139, \
                         if__140, if__141, lf_247, lf_248, lf_249, lf_250, \
                         lf_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * if__137[k]
                   + f_0 * lf_247[k];

        t_188[k] = -2.0 * if__138[k]
                   + f_0 * lf_248[k];

        t_189[k] = -2.0 * if__139[k]
                   + f_0 * lf_249[k];

        t_190[k] = -if__140[k]
                   + f_0 * lf_250[k];

        t_191[k] = -if__141[k]
                   + f_0 * lf_251[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, if__142, if__143, if__144, \
                         if__145, if__146, lf_252, lf_253, lf_254, lf_255, \
                         lf_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -if__142[k]
                   + f_0 * lf_252[k];

        t_193[k] = -if__143[k]
                   + f_0 * lf_253[k];

        t_194[k] = -if__144[k]
                   + f_0 * lf_254[k];

        t_195[k] = -if__145[k]
                   + f_0 * lf_255[k];

        t_196[k] = -if__146[k]
                   + f_0 * lf_256[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, t_202, if__147, if__148, if__149, \
                         lf_257, lf_258, lf_259, lf_260, lf_261, \
                         lf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -if__147[k]
                   + f_0 * lf_257[k];

        t_198[k] = -if__148[k]
                   + f_0 * lf_258[k];

        t_199[k] = -if__149[k]
                   + f_0 * lf_259[k];

        t_200[k] = f_0 * lf_260[k];

        t_201[k] = f_0 * lf_261[k];

        t_202[k] = f_0 * lf_262[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, t_209, lf_263, lf_264, \
                         lf_265, lf_266, lf_267, lf_268, lf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_0 * lf_263[k];

        t_204[k] = f_0 * lf_264[k];

        t_205[k] = f_0 * lf_265[k];

        t_206[k] = f_0 * lf_266[k];

        t_207[k] = f_0 * lf_267[k];

        t_208[k] = f_0 * lf_268[k];

        t_209[k] = f_0 * lf_269[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, if__150, if__151, if__152, \
                         if__153, if__154, lf_280, lf_281, lf_282, lf_283, \
                         lf_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -6.0 * if__150[k]
                   + f_0 * lf_280[k];

        t_211[k] = -6.0 * if__151[k]
                   + f_0 * lf_281[k];

        t_212[k] = -6.0 * if__152[k]
                   + f_0 * lf_282[k];

        t_213[k] = -6.0 * if__153[k]
                   + f_0 * lf_283[k];

        t_214[k] = -6.0 * if__154[k]
                   + f_0 * lf_284[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, if__155, if__156, if__157, \
                         if__158, if__159, lf_285, lf_286, lf_287, lf_288, \
                         lf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -6.0 * if__155[k]
                   + f_0 * lf_285[k];

        t_216[k] = -6.0 * if__156[k]
                   + f_0 * lf_286[k];

        t_217[k] = -6.0 * if__157[k]
                   + f_0 * lf_287[k];

        t_218[k] = -6.0 * if__158[k]
                   + f_0 * lf_288[k];

        t_219[k] = -6.0 * if__159[k]
                   + f_0 * lf_289[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, if__160, if__161, if__162, \
                         if__163, if__164, lf_290, lf_291, lf_292, lf_293, \
                         lf_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -5.0 * if__160[k]
                   + f_0 * lf_290[k];

        t_221[k] = -5.0 * if__161[k]
                   + f_0 * lf_291[k];

        t_222[k] = -5.0 * if__162[k]
                   + f_0 * lf_292[k];

        t_223[k] = -5.0 * if__163[k]
                   + f_0 * lf_293[k];

        t_224[k] = -5.0 * if__164[k]
                   + f_0 * lf_294[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, if__165, if__166, if__167, \
                         if__168, if__169, lf_295, lf_296, lf_297, lf_298, \
                         lf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -5.0 * if__165[k]
                   + f_0 * lf_295[k];

        t_226[k] = -5.0 * if__166[k]
                   + f_0 * lf_296[k];

        t_227[k] = -5.0 * if__167[k]
                   + f_0 * lf_297[k];

        t_228[k] = -5.0 * if__168[k]
                   + f_0 * lf_298[k];

        t_229[k] = -5.0 * if__169[k]
                   + f_0 * lf_299[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, if__170, if__171, if__172, \
                         if__173, if__174, lf_300, lf_301, lf_302, lf_303, \
                         lf_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -4.0 * if__170[k]
                   + f_0 * lf_300[k];

        t_231[k] = -4.0 * if__171[k]
                   + f_0 * lf_301[k];

        t_232[k] = -4.0 * if__172[k]
                   + f_0 * lf_302[k];

        t_233[k] = -4.0 * if__173[k]
                   + f_0 * lf_303[k];

        t_234[k] = -4.0 * if__174[k]
                   + f_0 * lf_304[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, if__175, if__176, if__177, \
                         if__178, if__179, lf_305, lf_306, lf_307, lf_308, \
                         lf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -4.0 * if__175[k]
                   + f_0 * lf_305[k];

        t_236[k] = -4.0 * if__176[k]
                   + f_0 * lf_306[k];

        t_237[k] = -4.0 * if__177[k]
                   + f_0 * lf_307[k];

        t_238[k] = -4.0 * if__178[k]
                   + f_0 * lf_308[k];

        t_239[k] = -4.0 * if__179[k]
                   + f_0 * lf_309[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, if__180, if__181, if__182, \
                         if__183, if__184, lf_310, lf_311, lf_312, lf_313, \
                         lf_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * if__180[k]
                   + f_0 * lf_310[k];

        t_241[k] = -3.0 * if__181[k]
                   + f_0 * lf_311[k];

        t_242[k] = -3.0 * if__182[k]
                   + f_0 * lf_312[k];

        t_243[k] = -3.0 * if__183[k]
                   + f_0 * lf_313[k];

        t_244[k] = -3.0 * if__184[k]
                   + f_0 * lf_314[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, if__185, if__186, if__187, \
                         if__188, if__189, lf_315, lf_316, lf_317, lf_318, \
                         lf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * if__185[k]
                   + f_0 * lf_315[k];

        t_246[k] = -3.0 * if__186[k]
                   + f_0 * lf_316[k];

        t_247[k] = -3.0 * if__187[k]
                   + f_0 * lf_317[k];

        t_248[k] = -3.0 * if__188[k]
                   + f_0 * lf_318[k];

        t_249[k] = -3.0 * if__189[k]
                   + f_0 * lf_319[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, if__190, if__191, if__192, \
                         if__193, if__194, lf_320, lf_321, lf_322, lf_323, \
                         lf_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * if__190[k]
                   + f_0 * lf_320[k];

        t_251[k] = -2.0 * if__191[k]
                   + f_0 * lf_321[k];

        t_252[k] = -2.0 * if__192[k]
                   + f_0 * lf_322[k];

        t_253[k] = -2.0 * if__193[k]
                   + f_0 * lf_323[k];

        t_254[k] = -2.0 * if__194[k]
                   + f_0 * lf_324[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, if__195, if__196, if__197, \
                         if__198, if__199, lf_325, lf_326, lf_327, lf_328, \
                         lf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * if__195[k]
                   + f_0 * lf_325[k];

        t_256[k] = -2.0 * if__196[k]
                   + f_0 * lf_326[k];

        t_257[k] = -2.0 * if__197[k]
                   + f_0 * lf_327[k];

        t_258[k] = -2.0 * if__198[k]
                   + f_0 * lf_328[k];

        t_259[k] = -2.0 * if__199[k]
                   + f_0 * lf_329[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, if__200, if__201, if__202, \
                         if__203, if__204, lf_330, lf_331, lf_332, lf_333, \
                         lf_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -if__200[k]
                   + f_0 * lf_330[k];

        t_261[k] = -if__201[k]
                   + f_0 * lf_331[k];

        t_262[k] = -if__202[k]
                   + f_0 * lf_332[k];

        t_263[k] = -if__203[k]
                   + f_0 * lf_333[k];

        t_264[k] = -if__204[k]
                   + f_0 * lf_334[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, if__205, if__206, if__207, \
                         if__208, if__209, lf_335, lf_336, lf_337, lf_338, \
                         lf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -if__205[k]
                   + f_0 * lf_335[k];

        t_266[k] = -if__206[k]
                   + f_0 * lf_336[k];

        t_267[k] = -if__207[k]
                   + f_0 * lf_337[k];

        t_268[k] = -if__208[k]
                   + f_0 * lf_338[k];

        t_269[k] = -if__209[k]
                   + f_0 * lf_339[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, t_276, t_277, lf_340, \
                         lf_341, lf_342, lf_343, lf_344, lf_345, lf_346, \
                         lf_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_0 * lf_340[k];

        t_271[k] = f_0 * lf_341[k];

        t_272[k] = f_0 * lf_342[k];

        t_273[k] = f_0 * lf_343[k];

        t_274[k] = f_0 * lf_344[k];

        t_275[k] = f_0 * lf_345[k];

        t_276[k] = f_0 * lf_346[k];

        t_277[k] = f_0 * lf_347[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, t_283, if__210, if__211, if__212, \
                         if__213, lf_348, lf_349, lf_360, lf_361, lf_362, \
                         lf_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_0 * lf_348[k];

        t_279[k] = f_0 * lf_349[k];

        t_280[k] = -7.0 * if__210[k]
                   + f_0 * lf_360[k];

        t_281[k] = -7.0 * if__211[k]
                   + f_0 * lf_361[k];

        t_282[k] = -7.0 * if__212[k]
                   + f_0 * lf_362[k];

        t_283[k] = -7.0 * if__213[k]
                   + f_0 * lf_363[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, if__214, if__215, if__216, \
                         if__217, if__218, lf_364, lf_365, lf_366, lf_367, \
                         lf_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = -7.0 * if__214[k]
                   + f_0 * lf_364[k];

        t_285[k] = -7.0 * if__215[k]
                   + f_0 * lf_365[k];

        t_286[k] = -7.0 * if__216[k]
                   + f_0 * lf_366[k];

        t_287[k] = -7.0 * if__217[k]
                   + f_0 * lf_367[k];

        t_288[k] = -7.0 * if__218[k]
                   + f_0 * lf_368[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, if__219, if__220, if__221, \
                         if__222, if__223, lf_369, lf_370, lf_371, lf_372, \
                         lf_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = -7.0 * if__219[k]
                   + f_0 * lf_369[k];

        t_290[k] = -6.0 * if__220[k]
                   + f_0 * lf_370[k];

        t_291[k] = -6.0 * if__221[k]
                   + f_0 * lf_371[k];

        t_292[k] = -6.0 * if__222[k]
                   + f_0 * lf_372[k];

        t_293[k] = -6.0 * if__223[k]
                   + f_0 * lf_373[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, if__224, if__225, if__226, \
                         if__227, if__228, lf_374, lf_375, lf_376, lf_377, \
                         lf_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = -6.0 * if__224[k]
                   + f_0 * lf_374[k];

        t_295[k] = -6.0 * if__225[k]
                   + f_0 * lf_375[k];

        t_296[k] = -6.0 * if__226[k]
                   + f_0 * lf_376[k];

        t_297[k] = -6.0 * if__227[k]
                   + f_0 * lf_377[k];

        t_298[k] = -6.0 * if__228[k]
                   + f_0 * lf_378[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, if__229, if__230, if__231, \
                         if__232, if__233, lf_379, lf_380, lf_381, lf_382, \
                         lf_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = -6.0 * if__229[k]
                   + f_0 * lf_379[k];

        t_300[k] = -5.0 * if__230[k]
                   + f_0 * lf_380[k];

        t_301[k] = -5.0 * if__231[k]
                   + f_0 * lf_381[k];

        t_302[k] = -5.0 * if__232[k]
                   + f_0 * lf_382[k];

        t_303[k] = -5.0 * if__233[k]
                   + f_0 * lf_383[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, if__234, if__235, if__236, \
                         if__237, if__238, lf_384, lf_385, lf_386, lf_387, \
                         lf_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = -5.0 * if__234[k]
                   + f_0 * lf_384[k];

        t_305[k] = -5.0 * if__235[k]
                   + f_0 * lf_385[k];

        t_306[k] = -5.0 * if__236[k]
                   + f_0 * lf_386[k];

        t_307[k] = -5.0 * if__237[k]
                   + f_0 * lf_387[k];

        t_308[k] = -5.0 * if__238[k]
                   + f_0 * lf_388[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, if__239, if__240, if__241, \
                         if__242, if__243, lf_389, lf_390, lf_391, lf_392, \
                         lf_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = -5.0 * if__239[k]
                   + f_0 * lf_389[k];

        t_310[k] = -4.0 * if__240[k]
                   + f_0 * lf_390[k];

        t_311[k] = -4.0 * if__241[k]
                   + f_0 * lf_391[k];

        t_312[k] = -4.0 * if__242[k]
                   + f_0 * lf_392[k];

        t_313[k] = -4.0 * if__243[k]
                   + f_0 * lf_393[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, t_318, if__244, if__245, if__246, \
                         if__247, if__248, lf_394, lf_395, lf_396, lf_397, \
                         lf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = -4.0 * if__244[k]
                   + f_0 * lf_394[k];

        t_315[k] = -4.0 * if__245[k]
                   + f_0 * lf_395[k];

        t_316[k] = -4.0 * if__246[k]
                   + f_0 * lf_396[k];

        t_317[k] = -4.0 * if__247[k]
                   + f_0 * lf_397[k];

        t_318[k] = -4.0 * if__248[k]
                   + f_0 * lf_398[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, t_323, if__249, if__250, if__251, \
                         if__252, if__253, lf_399, lf_400, lf_401, lf_402, \
                         lf_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = -4.0 * if__249[k]
                   + f_0 * lf_399[k];

        t_320[k] = -3.0 * if__250[k]
                   + f_0 * lf_400[k];

        t_321[k] = -3.0 * if__251[k]
                   + f_0 * lf_401[k];

        t_322[k] = -3.0 * if__252[k]
                   + f_0 * lf_402[k];

        t_323[k] = -3.0 * if__253[k]
                   + f_0 * lf_403[k];
    }
}

static auto
compute_prim_geom_10_kf_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t if_, const size_t lf,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *lf_404 = buffer.data(lf + 404);
    const auto *lf_405 = buffer.data(lf + 405);
    const auto *lf_406 = buffer.data(lf + 406);
    const auto *lf_407 = buffer.data(lf + 407);
    const auto *lf_408 = buffer.data(lf + 408);
    const auto *lf_409 = buffer.data(lf + 409);
    const auto *lf_410 = buffer.data(lf + 410);
    const auto *lf_411 = buffer.data(lf + 411);
    const auto *lf_412 = buffer.data(lf + 412);
    const auto *lf_413 = buffer.data(lf + 413);
    const auto *lf_414 = buffer.data(lf + 414);
    const auto *lf_415 = buffer.data(lf + 415);
    const auto *lf_416 = buffer.data(lf + 416);
    const auto *lf_417 = buffer.data(lf + 417);
    const auto *lf_418 = buffer.data(lf + 418);
    const auto *lf_419 = buffer.data(lf + 419);
    const auto *lf_420 = buffer.data(lf + 420);
    const auto *lf_421 = buffer.data(lf + 421);
    const auto *lf_422 = buffer.data(lf + 422);
    const auto *lf_423 = buffer.data(lf + 423);
    const auto *lf_424 = buffer.data(lf + 424);
    const auto *lf_425 = buffer.data(lf + 425);
    const auto *lf_426 = buffer.data(lf + 426);
    const auto *lf_427 = buffer.data(lf + 427);
    const auto *lf_428 = buffer.data(lf + 428);
    const auto *lf_429 = buffer.data(lf + 429);
    const auto *lf_430 = buffer.data(lf + 430);
    const auto *lf_431 = buffer.data(lf + 431);
    const auto *lf_432 = buffer.data(lf + 432);
    const auto *lf_433 = buffer.data(lf + 433);
    const auto *lf_434 = buffer.data(lf + 434);
    const auto *lf_435 = buffer.data(lf + 435);
    const auto *lf_436 = buffer.data(lf + 436);
    const auto *lf_437 = buffer.data(lf + 437);
    const auto *lf_438 = buffer.data(lf + 438);
    const auto *lf_439 = buffer.data(lf + 439);

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, if__254, if__255, if__256, \
                         if__257, if__258, lf_404, lf_405, lf_406, lf_407, \
                         lf_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = -3.0 * if__254[k]
                   + f_0 * lf_404[k];

        t_325[k] = -3.0 * if__255[k]
                   + f_0 * lf_405[k];

        t_326[k] = -3.0 * if__256[k]
                   + f_0 * lf_406[k];

        t_327[k] = -3.0 * if__257[k]
                   + f_0 * lf_407[k];

        t_328[k] = -3.0 * if__258[k]
                   + f_0 * lf_408[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, t_333, if__259, if__260, if__261, \
                         if__262, if__263, lf_409, lf_410, lf_411, lf_412, \
                         lf_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = -3.0 * if__259[k]
                   + f_0 * lf_409[k];

        t_330[k] = -2.0 * if__260[k]
                   + f_0 * lf_410[k];

        t_331[k] = -2.0 * if__261[k]
                   + f_0 * lf_411[k];

        t_332[k] = -2.0 * if__262[k]
                   + f_0 * lf_412[k];

        t_333[k] = -2.0 * if__263[k]
                   + f_0 * lf_413[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, if__264, if__265, if__266, \
                         if__267, if__268, lf_414, lf_415, lf_416, lf_417, \
                         lf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = -2.0 * if__264[k]
                   + f_0 * lf_414[k];

        t_335[k] = -2.0 * if__265[k]
                   + f_0 * lf_415[k];

        t_336[k] = -2.0 * if__266[k]
                   + f_0 * lf_416[k];

        t_337[k] = -2.0 * if__267[k]
                   + f_0 * lf_417[k];

        t_338[k] = -2.0 * if__268[k]
                   + f_0 * lf_418[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, t_343, if__269, if__270, if__271, \
                         if__272, if__273, lf_419, lf_420, lf_421, lf_422, \
                         lf_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = -2.0 * if__269[k]
                   + f_0 * lf_419[k];

        t_340[k] = -if__270[k]
                   + f_0 * lf_420[k];

        t_341[k] = -if__271[k]
                   + f_0 * lf_421[k];

        t_342[k] = -if__272[k]
                   + f_0 * lf_422[k];

        t_343[k] = -if__273[k]
                   + f_0 * lf_423[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, if__274, if__275, if__276, \
                         if__277, if__278, lf_424, lf_425, lf_426, lf_427, \
                         lf_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = -if__274[k]
                   + f_0 * lf_424[k];

        t_345[k] = -if__275[k]
                   + f_0 * lf_425[k];

        t_346[k] = -if__276[k]
                   + f_0 * lf_426[k];

        t_347[k] = -if__277[k]
                   + f_0 * lf_427[k];

        t_348[k] = -if__278[k]
                   + f_0 * lf_428[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, t_353, t_354, t_355, if__279, lf_429, \
                         lf_430, lf_431, lf_432, lf_433, lf_434, \
                         lf_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = -if__279[k]
                   + f_0 * lf_429[k];

        t_350[k] = f_0 * lf_430[k];

        t_351[k] = f_0 * lf_431[k];

        t_352[k] = f_0 * lf_432[k];

        t_353[k] = f_0 * lf_433[k];

        t_354[k] = f_0 * lf_434[k];

        t_355[k] = f_0 * lf_435[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, lf_436, lf_437, lf_438, \
                         lf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_0 * lf_436[k];

        t_357[k] = f_0 * lf_437[k];

        t_358[k] = f_0 * lf_438[k];

        t_359[k] = f_0 * lf_439[k];
    }
}

auto
compute_prim_geom_10_kf_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t if_, const size_t lf,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kf_electron_repulsion_1_piece0(buffer, target, if_, lf, ncols, alpha);

    compute_prim_geom_10_kf_electron_repulsion_1_piece1(buffer, target, if_, lf, ncols, alpha);

    compute_prim_geom_10_kf_electron_repulsion_1_piece2(buffer, target, if_, lf, ncols, alpha);
}

static auto
compute_prim_geom_10_kf_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t if_, const size_t lf,
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

    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_110 = buffer.data(lf + 110);
    const auto *lf_111 = buffer.data(lf + 111);
    const auto *lf_112 = buffer.data(lf + 112);
    const auto *lf_113 = buffer.data(lf + 113);
    const auto *lf_114 = buffer.data(lf + 114);
    const auto *lf_115 = buffer.data(lf + 115);
    const auto *lf_116 = buffer.data(lf + 116);
    const auto *lf_117 = buffer.data(lf + 117);
    const auto *lf_118 = buffer.data(lf + 118);
    const auto *lf_119 = buffer.data(lf + 119);
    const auto *lf_120 = buffer.data(lf + 120);
    const auto *lf_121 = buffer.data(lf + 121);
    const auto *lf_122 = buffer.data(lf + 122);
    const auto *lf_123 = buffer.data(lf + 123);
    const auto *lf_124 = buffer.data(lf + 124);
    const auto *lf_125 = buffer.data(lf + 125);
    const auto *lf_126 = buffer.data(lf + 126);
    const auto *lf_127 = buffer.data(lf + 127);
    const auto *lf_128 = buffer.data(lf + 128);
    const auto *lf_129 = buffer.data(lf + 129);
    const auto *lf_130 = buffer.data(lf + 130);
    const auto *lf_131 = buffer.data(lf + 131);
    const auto *lf_132 = buffer.data(lf + 132);
    const auto *lf_133 = buffer.data(lf + 133);
    const auto *lf_134 = buffer.data(lf + 134);
    const auto *lf_135 = buffer.data(lf + 135);
    const auto *lf_136 = buffer.data(lf + 136);
    const auto *lf_137 = buffer.data(lf + 137);
    const auto *lf_138 = buffer.data(lf + 138);
    const auto *lf_139 = buffer.data(lf + 139);
    const auto *lf_140 = buffer.data(lf + 140);
    const auto *lf_141 = buffer.data(lf + 141);
    const auto *lf_142 = buffer.data(lf + 142);
    const auto *lf_143 = buffer.data(lf + 143);
    const auto *lf_144 = buffer.data(lf + 144);
    const auto *lf_145 = buffer.data(lf + 145);
    const auto *lf_146 = buffer.data(lf + 146);
    const auto *lf_147 = buffer.data(lf + 147);
    const auto *lf_148 = buffer.data(lf + 148);
    const auto *lf_149 = buffer.data(lf + 149);
    const auto *lf_160 = buffer.data(lf + 160);
    const auto *lf_161 = buffer.data(lf + 161);
    const auto *lf_162 = buffer.data(lf + 162);
    const auto *lf_163 = buffer.data(lf + 163);
    const auto *lf_164 = buffer.data(lf + 164);
    const auto *lf_165 = buffer.data(lf + 165);
    const auto *lf_166 = buffer.data(lf + 166);
    const auto *lf_167 = buffer.data(lf + 167);
    const auto *lf_168 = buffer.data(lf + 168);
    const auto *lf_169 = buffer.data(lf + 169);
    const auto *lf_170 = buffer.data(lf + 170);
    const auto *lf_171 = buffer.data(lf + 171);
    const auto *lf_172 = buffer.data(lf + 172);
    const auto *lf_173 = buffer.data(lf + 173);
    const auto *lf_174 = buffer.data(lf + 174);
    const auto *lf_175 = buffer.data(lf + 175);
    const auto *lf_176 = buffer.data(lf + 176);
    const auto *lf_177 = buffer.data(lf + 177);
    const auto *lf_178 = buffer.data(lf + 178);
    const auto *lf_179 = buffer.data(lf + 179);
    const auto *lf_180 = buffer.data(lf + 180);
    const auto *lf_181 = buffer.data(lf + 181);
    const auto *lf_182 = buffer.data(lf + 182);
    const auto *lf_183 = buffer.data(lf + 183);
    const auto *lf_184 = buffer.data(lf + 184);
    const auto *lf_185 = buffer.data(lf + 185);
    const auto *lf_186 = buffer.data(lf + 186);
    const auto *lf_187 = buffer.data(lf + 187);
    const auto *lf_188 = buffer.data(lf + 188);
    const auto *lf_189 = buffer.data(lf + 189);
    const auto *lf_190 = buffer.data(lf + 190);
    const auto *lf_191 = buffer.data(lf + 191);
    const auto *lf_192 = buffer.data(lf + 192);
    const auto *lf_193 = buffer.data(lf + 193);
    const auto *lf_194 = buffer.data(lf + 194);
    const auto *lf_195 = buffer.data(lf + 195);
    const auto *lf_196 = buffer.data(lf + 196);
    const auto *lf_197 = buffer.data(lf + 197);
    const auto *lf_198 = buffer.data(lf + 198);
    const auto *lf_199 = buffer.data(lf + 199);
    const auto *lf_200 = buffer.data(lf + 200);
    const auto *lf_201 = buffer.data(lf + 201);
    const auto *lf_202 = buffer.data(lf + 202);
    const auto *lf_203 = buffer.data(lf + 203);
    const auto *lf_204 = buffer.data(lf + 204);
    const auto *lf_205 = buffer.data(lf + 205);
    const auto *lf_206 = buffer.data(lf + 206);
    const auto *lf_207 = buffer.data(lf + 207);
    const auto *lf_208 = buffer.data(lf + 208);
    const auto *lf_209 = buffer.data(lf + 209);
    const auto *lf_220 = buffer.data(lf + 220);
    const auto *lf_221 = buffer.data(lf + 221);
    const auto *lf_222 = buffer.data(lf + 222);
    const auto *lf_223 = buffer.data(lf + 223);
    const auto *lf_224 = buffer.data(lf + 224);
    const auto *lf_225 = buffer.data(lf + 225);
    const auto *lf_226 = buffer.data(lf + 226);
    const auto *lf_227 = buffer.data(lf + 227);
    const auto *lf_228 = buffer.data(lf + 228);
    const auto *lf_229 = buffer.data(lf + 229);
    const auto *lf_230 = buffer.data(lf + 230);
    const auto *lf_231 = buffer.data(lf + 231);
    const auto *lf_232 = buffer.data(lf + 232);
    const auto *lf_233 = buffer.data(lf + 233);
    const auto *lf_234 = buffer.data(lf + 234);
    const auto *lf_235 = buffer.data(lf + 235);
    const auto *lf_236 = buffer.data(lf + 236);
    const auto *lf_237 = buffer.data(lf + 237);
    const auto *lf_238 = buffer.data(lf + 238);
    const auto *lf_239 = buffer.data(lf + 239);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, lf_20, lf_21, lf_22, lf_23, \
                         lf_24, lf_25, lf_26, lf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lf_20[k];

        t_1[k] = f_0 * lf_21[k];

        t_2[k] = f_0 * lf_22[k];

        t_3[k] = f_0 * lf_23[k];

        t_4[k] = f_0 * lf_24[k];

        t_5[k] = f_0 * lf_25[k];

        t_6[k] = f_0 * lf_26[k];

        t_7[k] = f_0 * lf_27[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, lf_28, lf_29, lf_40, \
                         lf_41, lf_42, lf_43, lf_44, lf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * lf_28[k];

        t_9[k] = f_0 * lf_29[k];

        t_10[k] = f_0 * lf_40[k];

        t_11[k] = f_0 * lf_41[k];

        t_12[k] = f_0 * lf_42[k];

        t_13[k] = f_0 * lf_43[k];

        t_14[k] = f_0 * lf_44[k];

        t_15[k] = f_0 * lf_45[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, if__0, if__1, lf_46, lf_47, \
                         lf_48, lf_49, lf_50, lf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * lf_46[k];

        t_17[k] = f_0 * lf_47[k];

        t_18[k] = f_0 * lf_48[k];

        t_19[k] = f_0 * lf_49[k];

        t_20[k] = -if__0[k]
                  + f_0 * lf_50[k];

        t_21[k] = -if__1[k]
                  + f_0 * lf_51[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, if__2, if__3, if__4, if__5, if__6, \
                         lf_52, lf_53, lf_54, lf_55, lf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -if__2[k]
                  + f_0 * lf_52[k];

        t_23[k] = -if__3[k]
                  + f_0 * lf_53[k];

        t_24[k] = -if__4[k]
                  + f_0 * lf_54[k];

        t_25[k] = -if__5[k]
                  + f_0 * lf_55[k];

        t_26[k] = -if__6[k]
                  + f_0 * lf_56[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, if__7, if__8, if__9, lf_57, \
                         lf_58, lf_59, lf_70, lf_71, lf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -if__7[k]
                  + f_0 * lf_57[k];

        t_28[k] = -if__8[k]
                  + f_0 * lf_58[k];

        t_29[k] = -if__9[k]
                  + f_0 * lf_59[k];

        t_30[k] = f_0 * lf_70[k];

        t_31[k] = f_0 * lf_71[k];

        t_32[k] = f_0 * lf_72[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, t_39, lf_73, lf_74, lf_75, lf_76, \
                         lf_77, lf_78, lf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * lf_73[k];

        t_34[k] = f_0 * lf_74[k];

        t_35[k] = f_0 * lf_75[k];

        t_36[k] = f_0 * lf_76[k];

        t_37[k] = f_0 * lf_77[k];

        t_38[k] = f_0 * lf_78[k];

        t_39[k] = f_0 * lf_79[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, if__10, if__11, if__12, if__13, if__14, \
                         lf_80, lf_81, lf_82, lf_83, lf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -if__10[k]
                  + f_0 * lf_80[k];

        t_41[k] = -if__11[k]
                  + f_0 * lf_81[k];

        t_42[k] = -if__12[k]
                  + f_0 * lf_82[k];

        t_43[k] = -if__13[k]
                  + f_0 * lf_83[k];

        t_44[k] = -if__14[k]
                  + f_0 * lf_84[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, if__15, if__16, if__17, if__18, if__19, \
                         lf_85, lf_86, lf_87, lf_88, lf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -if__15[k]
                  + f_0 * lf_85[k];

        t_46[k] = -if__16[k]
                  + f_0 * lf_86[k];

        t_47[k] = -if__17[k]
                  + f_0 * lf_87[k];

        t_48[k] = -if__18[k]
                  + f_0 * lf_88[k];

        t_49[k] = -if__19[k]
                  + f_0 * lf_89[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, if__20, if__21, if__22, if__23, if__24, \
                         lf_90, lf_91, lf_92, lf_93, lf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * if__20[k]
                  + f_0 * lf_90[k];

        t_51[k] = -2.0 * if__21[k]
                  + f_0 * lf_91[k];

        t_52[k] = -2.0 * if__22[k]
                  + f_0 * lf_92[k];

        t_53[k] = -2.0 * if__23[k]
                  + f_0 * lf_93[k];

        t_54[k] = -2.0 * if__24[k]
                  + f_0 * lf_94[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, if__25, if__26, if__27, if__28, if__29, \
                         lf_95, lf_96, lf_97, lf_98, lf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * if__25[k]
                  + f_0 * lf_95[k];

        t_56[k] = -2.0 * if__26[k]
                  + f_0 * lf_96[k];

        t_57[k] = -2.0 * if__27[k]
                  + f_0 * lf_97[k];

        t_58[k] = -2.0 * if__28[k]
                  + f_0 * lf_98[k];

        t_59[k] = -2.0 * if__29[k]
                  + f_0 * lf_99[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, t_66, t_67, lf_110, lf_111, \
                         lf_112, lf_113, lf_114, lf_115, lf_116, \
                         lf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * lf_110[k];

        t_61[k] = f_0 * lf_111[k];

        t_62[k] = f_0 * lf_112[k];

        t_63[k] = f_0 * lf_113[k];

        t_64[k] = f_0 * lf_114[k];

        t_65[k] = f_0 * lf_115[k];

        t_66[k] = f_0 * lf_116[k];

        t_67[k] = f_0 * lf_117[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, if__30, if__31, if__32, if__33, \
                         lf_118, lf_119, lf_120, lf_121, lf_122, \
                         lf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * lf_118[k];

        t_69[k] = f_0 * lf_119[k];

        t_70[k] = -if__30[k]
                  + f_0 * lf_120[k];

        t_71[k] = -if__31[k]
                  + f_0 * lf_121[k];

        t_72[k] = -if__32[k]
                  + f_0 * lf_122[k];

        t_73[k] = -if__33[k]
                  + f_0 * lf_123[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, if__34, if__35, if__36, if__37, if__38, \
                         lf_124, lf_125, lf_126, lf_127, lf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = -if__34[k]
                  + f_0 * lf_124[k];

        t_75[k] = -if__35[k]
                  + f_0 * lf_125[k];

        t_76[k] = -if__36[k]
                  + f_0 * lf_126[k];

        t_77[k] = -if__37[k]
                  + f_0 * lf_127[k];

        t_78[k] = -if__38[k]
                  + f_0 * lf_128[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, if__39, if__40, if__41, if__42, if__43, \
                         lf_129, lf_130, lf_131, lf_132, lf_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -if__39[k]
                  + f_0 * lf_129[k];

        t_80[k] = -2.0 * if__40[k]
                  + f_0 * lf_130[k];

        t_81[k] = -2.0 * if__41[k]
                  + f_0 * lf_131[k];

        t_82[k] = -2.0 * if__42[k]
                  + f_0 * lf_132[k];

        t_83[k] = -2.0 * if__43[k]
                  + f_0 * lf_133[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, if__44, if__45, if__46, if__47, if__48, \
                         lf_134, lf_135, lf_136, lf_137, lf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -2.0 * if__44[k]
                  + f_0 * lf_134[k];

        t_85[k] = -2.0 * if__45[k]
                  + f_0 * lf_135[k];

        t_86[k] = -2.0 * if__46[k]
                  + f_0 * lf_136[k];

        t_87[k] = -2.0 * if__47[k]
                  + f_0 * lf_137[k];

        t_88[k] = -2.0 * if__48[k]
                  + f_0 * lf_138[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, if__49, if__50, if__51, if__52, if__53, \
                         lf_139, lf_140, lf_141, lf_142, lf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -2.0 * if__49[k]
                  + f_0 * lf_139[k];

        t_90[k] = -3.0 * if__50[k]
                  + f_0 * lf_140[k];

        t_91[k] = -3.0 * if__51[k]
                  + f_0 * lf_141[k];

        t_92[k] = -3.0 * if__52[k]
                  + f_0 * lf_142[k];

        t_93[k] = -3.0 * if__53[k]
                  + f_0 * lf_143[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, if__54, if__55, if__56, if__57, if__58, \
                         lf_144, lf_145, lf_146, lf_147, lf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -3.0 * if__54[k]
                  + f_0 * lf_144[k];

        t_95[k] = -3.0 * if__55[k]
                  + f_0 * lf_145[k];

        t_96[k] = -3.0 * if__56[k]
                  + f_0 * lf_146[k];

        t_97[k] = -3.0 * if__57[k]
                  + f_0 * lf_147[k];

        t_98[k] = -3.0 * if__58[k]
                  + f_0 * lf_148[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, t_105, if__59, lf_149, \
                         lf_160, lf_161, lf_162, lf_163, lf_164, \
                         lf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -3.0 * if__59[k]
                  + f_0 * lf_149[k];

        t_100[k] = f_0 * lf_160[k];

        t_101[k] = f_0 * lf_161[k];

        t_102[k] = f_0 * lf_162[k];

        t_103[k] = f_0 * lf_163[k];

        t_104[k] = f_0 * lf_164[k];

        t_105[k] = f_0 * lf_165[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, if__60, if__61, lf_166, \
                         lf_167, lf_168, lf_169, lf_170, lf_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * lf_166[k];

        t_107[k] = f_0 * lf_167[k];

        t_108[k] = f_0 * lf_168[k];

        t_109[k] = f_0 * lf_169[k];

        t_110[k] = -if__60[k]
                   + f_0 * lf_170[k];

        t_111[k] = -if__61[k]
                   + f_0 * lf_171[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, if__62, if__63, if__64, if__65, \
                         if__66, lf_172, lf_173, lf_174, lf_175, \
                         lf_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -if__62[k]
                   + f_0 * lf_172[k];

        t_113[k] = -if__63[k]
                   + f_0 * lf_173[k];

        t_114[k] = -if__64[k]
                   + f_0 * lf_174[k];

        t_115[k] = -if__65[k]
                   + f_0 * lf_175[k];

        t_116[k] = -if__66[k]
                   + f_0 * lf_176[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, if__67, if__68, if__69, if__70, \
                         if__71, lf_177, lf_178, lf_179, lf_180, \
                         lf_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -if__67[k]
                   + f_0 * lf_177[k];

        t_118[k] = -if__68[k]
                   + f_0 * lf_178[k];

        t_119[k] = -if__69[k]
                   + f_0 * lf_179[k];

        t_120[k] = -2.0 * if__70[k]
                   + f_0 * lf_180[k];

        t_121[k] = -2.0 * if__71[k]
                   + f_0 * lf_181[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, if__72, if__73, if__74, if__75, \
                         if__76, lf_182, lf_183, lf_184, lf_185, \
                         lf_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * if__72[k]
                   + f_0 * lf_182[k];

        t_123[k] = -2.0 * if__73[k]
                   + f_0 * lf_183[k];

        t_124[k] = -2.0 * if__74[k]
                   + f_0 * lf_184[k];

        t_125[k] = -2.0 * if__75[k]
                   + f_0 * lf_185[k];

        t_126[k] = -2.0 * if__76[k]
                   + f_0 * lf_186[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, if__77, if__78, if__79, if__80, \
                         if__81, lf_187, lf_188, lf_189, lf_190, \
                         lf_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * if__77[k]
                   + f_0 * lf_187[k];

        t_128[k] = -2.0 * if__78[k]
                   + f_0 * lf_188[k];

        t_129[k] = -2.0 * if__79[k]
                   + f_0 * lf_189[k];

        t_130[k] = -3.0 * if__80[k]
                   + f_0 * lf_190[k];

        t_131[k] = -3.0 * if__81[k]
                   + f_0 * lf_191[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, if__82, if__83, if__84, if__85, \
                         if__86, lf_192, lf_193, lf_194, lf_195, \
                         lf_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -3.0 * if__82[k]
                   + f_0 * lf_192[k];

        t_133[k] = -3.0 * if__83[k]
                   + f_0 * lf_193[k];

        t_134[k] = -3.0 * if__84[k]
                   + f_0 * lf_194[k];

        t_135[k] = -3.0 * if__85[k]
                   + f_0 * lf_195[k];

        t_136[k] = -3.0 * if__86[k]
                   + f_0 * lf_196[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, if__87, if__88, if__89, if__90, \
                         if__91, lf_197, lf_198, lf_199, lf_200, \
                         lf_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -3.0 * if__87[k]
                   + f_0 * lf_197[k];

        t_138[k] = -3.0 * if__88[k]
                   + f_0 * lf_198[k];

        t_139[k] = -3.0 * if__89[k]
                   + f_0 * lf_199[k];

        t_140[k] = -4.0 * if__90[k]
                   + f_0 * lf_200[k];

        t_141[k] = -4.0 * if__91[k]
                   + f_0 * lf_201[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, if__92, if__93, if__94, if__95, \
                         if__96, lf_202, lf_203, lf_204, lf_205, \
                         lf_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -4.0 * if__92[k]
                   + f_0 * lf_202[k];

        t_143[k] = -4.0 * if__93[k]
                   + f_0 * lf_203[k];

        t_144[k] = -4.0 * if__94[k]
                   + f_0 * lf_204[k];

        t_145[k] = -4.0 * if__95[k]
                   + f_0 * lf_205[k];

        t_146[k] = -4.0 * if__96[k]
                   + f_0 * lf_206[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, t_152, if__97, if__98, if__99, \
                         lf_207, lf_208, lf_209, lf_220, lf_221, \
                         lf_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -4.0 * if__97[k]
                   + f_0 * lf_207[k];

        t_148[k] = -4.0 * if__98[k]
                   + f_0 * lf_208[k];

        t_149[k] = -4.0 * if__99[k]
                   + f_0 * lf_209[k];

        t_150[k] = f_0 * lf_220[k];

        t_151[k] = f_0 * lf_221[k];

        t_152[k] = f_0 * lf_222[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, t_158, t_159, lf_223, lf_224, \
                         lf_225, lf_226, lf_227, lf_228, lf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_0 * lf_223[k];

        t_154[k] = f_0 * lf_224[k];

        t_155[k] = f_0 * lf_225[k];

        t_156[k] = f_0 * lf_226[k];

        t_157[k] = f_0 * lf_227[k];

        t_158[k] = f_0 * lf_228[k];

        t_159[k] = f_0 * lf_229[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, if__100, if__101, if__102, \
                         if__103, if__104, lf_230, lf_231, lf_232, lf_233, \
                         lf_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -if__100[k]
                   + f_0 * lf_230[k];

        t_161[k] = -if__101[k]
                   + f_0 * lf_231[k];

        t_162[k] = -if__102[k]
                   + f_0 * lf_232[k];

        t_163[k] = -if__103[k]
                   + f_0 * lf_233[k];

        t_164[k] = -if__104[k]
                   + f_0 * lf_234[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, if__105, if__106, if__107, \
                         if__108, if__109, lf_235, lf_236, lf_237, lf_238, \
                         lf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -if__105[k]
                   + f_0 * lf_235[k];

        t_166[k] = -if__106[k]
                   + f_0 * lf_236[k];

        t_167[k] = -if__107[k]
                   + f_0 * lf_237[k];

        t_168[k] = -if__108[k]
                   + f_0 * lf_238[k];

        t_169[k] = -if__109[k]
                   + f_0 * lf_239[k];
    }
}

static auto
compute_prim_geom_10_kf_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t if_, const size_t lf,
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

    const auto *lf_240 = buffer.data(lf + 240);
    const auto *lf_241 = buffer.data(lf + 241);
    const auto *lf_242 = buffer.data(lf + 242);
    const auto *lf_243 = buffer.data(lf + 243);
    const auto *lf_244 = buffer.data(lf + 244);
    const auto *lf_245 = buffer.data(lf + 245);
    const auto *lf_246 = buffer.data(lf + 246);
    const auto *lf_247 = buffer.data(lf + 247);
    const auto *lf_248 = buffer.data(lf + 248);
    const auto *lf_249 = buffer.data(lf + 249);
    const auto *lf_250 = buffer.data(lf + 250);
    const auto *lf_251 = buffer.data(lf + 251);
    const auto *lf_252 = buffer.data(lf + 252);
    const auto *lf_253 = buffer.data(lf + 253);
    const auto *lf_254 = buffer.data(lf + 254);
    const auto *lf_255 = buffer.data(lf + 255);
    const auto *lf_256 = buffer.data(lf + 256);
    const auto *lf_257 = buffer.data(lf + 257);
    const auto *lf_258 = buffer.data(lf + 258);
    const auto *lf_259 = buffer.data(lf + 259);
    const auto *lf_260 = buffer.data(lf + 260);
    const auto *lf_261 = buffer.data(lf + 261);
    const auto *lf_262 = buffer.data(lf + 262);
    const auto *lf_263 = buffer.data(lf + 263);
    const auto *lf_264 = buffer.data(lf + 264);
    const auto *lf_265 = buffer.data(lf + 265);
    const auto *lf_266 = buffer.data(lf + 266);
    const auto *lf_267 = buffer.data(lf + 267);
    const auto *lf_268 = buffer.data(lf + 268);
    const auto *lf_269 = buffer.data(lf + 269);
    const auto *lf_270 = buffer.data(lf + 270);
    const auto *lf_271 = buffer.data(lf + 271);
    const auto *lf_272 = buffer.data(lf + 272);
    const auto *lf_273 = buffer.data(lf + 273);
    const auto *lf_274 = buffer.data(lf + 274);
    const auto *lf_275 = buffer.data(lf + 275);
    const auto *lf_276 = buffer.data(lf + 276);
    const auto *lf_277 = buffer.data(lf + 277);
    const auto *lf_278 = buffer.data(lf + 278);
    const auto *lf_279 = buffer.data(lf + 279);
    const auto *lf_290 = buffer.data(lf + 290);
    const auto *lf_291 = buffer.data(lf + 291);
    const auto *lf_292 = buffer.data(lf + 292);
    const auto *lf_293 = buffer.data(lf + 293);
    const auto *lf_294 = buffer.data(lf + 294);
    const auto *lf_295 = buffer.data(lf + 295);
    const auto *lf_296 = buffer.data(lf + 296);
    const auto *lf_297 = buffer.data(lf + 297);
    const auto *lf_298 = buffer.data(lf + 298);
    const auto *lf_299 = buffer.data(lf + 299);
    const auto *lf_300 = buffer.data(lf + 300);
    const auto *lf_301 = buffer.data(lf + 301);
    const auto *lf_302 = buffer.data(lf + 302);
    const auto *lf_303 = buffer.data(lf + 303);
    const auto *lf_304 = buffer.data(lf + 304);
    const auto *lf_305 = buffer.data(lf + 305);
    const auto *lf_306 = buffer.data(lf + 306);
    const auto *lf_307 = buffer.data(lf + 307);
    const auto *lf_308 = buffer.data(lf + 308);
    const auto *lf_309 = buffer.data(lf + 309);
    const auto *lf_310 = buffer.data(lf + 310);
    const auto *lf_311 = buffer.data(lf + 311);
    const auto *lf_312 = buffer.data(lf + 312);
    const auto *lf_313 = buffer.data(lf + 313);
    const auto *lf_314 = buffer.data(lf + 314);
    const auto *lf_315 = buffer.data(lf + 315);
    const auto *lf_316 = buffer.data(lf + 316);
    const auto *lf_317 = buffer.data(lf + 317);
    const auto *lf_318 = buffer.data(lf + 318);
    const auto *lf_319 = buffer.data(lf + 319);
    const auto *lf_320 = buffer.data(lf + 320);
    const auto *lf_321 = buffer.data(lf + 321);
    const auto *lf_322 = buffer.data(lf + 322);
    const auto *lf_323 = buffer.data(lf + 323);
    const auto *lf_324 = buffer.data(lf + 324);
    const auto *lf_325 = buffer.data(lf + 325);
    const auto *lf_326 = buffer.data(lf + 326);
    const auto *lf_327 = buffer.data(lf + 327);
    const auto *lf_328 = buffer.data(lf + 328);
    const auto *lf_329 = buffer.data(lf + 329);
    const auto *lf_330 = buffer.data(lf + 330);
    const auto *lf_331 = buffer.data(lf + 331);
    const auto *lf_332 = buffer.data(lf + 332);
    const auto *lf_333 = buffer.data(lf + 333);
    const auto *lf_334 = buffer.data(lf + 334);
    const auto *lf_335 = buffer.data(lf + 335);
    const auto *lf_336 = buffer.data(lf + 336);
    const auto *lf_337 = buffer.data(lf + 337);
    const auto *lf_338 = buffer.data(lf + 338);
    const auto *lf_339 = buffer.data(lf + 339);
    const auto *lf_340 = buffer.data(lf + 340);
    const auto *lf_341 = buffer.data(lf + 341);
    const auto *lf_342 = buffer.data(lf + 342);
    const auto *lf_343 = buffer.data(lf + 343);
    const auto *lf_344 = buffer.data(lf + 344);
    const auto *lf_345 = buffer.data(lf + 345);
    const auto *lf_346 = buffer.data(lf + 346);
    const auto *lf_347 = buffer.data(lf + 347);
    const auto *lf_348 = buffer.data(lf + 348);
    const auto *lf_349 = buffer.data(lf + 349);
    const auto *lf_350 = buffer.data(lf + 350);
    const auto *lf_351 = buffer.data(lf + 351);
    const auto *lf_352 = buffer.data(lf + 352);
    const auto *lf_353 = buffer.data(lf + 353);
    const auto *lf_354 = buffer.data(lf + 354);
    const auto *lf_355 = buffer.data(lf + 355);
    const auto *lf_356 = buffer.data(lf + 356);
    const auto *lf_357 = buffer.data(lf + 357);
    const auto *lf_358 = buffer.data(lf + 358);
    const auto *lf_359 = buffer.data(lf + 359);
    const auto *lf_370 = buffer.data(lf + 370);
    const auto *lf_371 = buffer.data(lf + 371);
    const auto *lf_372 = buffer.data(lf + 372);
    const auto *lf_373 = buffer.data(lf + 373);
    const auto *lf_374 = buffer.data(lf + 374);
    const auto *lf_375 = buffer.data(lf + 375);
    const auto *lf_376 = buffer.data(lf + 376);
    const auto *lf_377 = buffer.data(lf + 377);
    const auto *lf_378 = buffer.data(lf + 378);
    const auto *lf_379 = buffer.data(lf + 379);
    const auto *lf_380 = buffer.data(lf + 380);
    const auto *lf_381 = buffer.data(lf + 381);
    const auto *lf_382 = buffer.data(lf + 382);
    const auto *lf_383 = buffer.data(lf + 383);
    const auto *lf_384 = buffer.data(lf + 384);
    const auto *lf_385 = buffer.data(lf + 385);
    const auto *lf_386 = buffer.data(lf + 386);
    const auto *lf_387 = buffer.data(lf + 387);
    const auto *lf_388 = buffer.data(lf + 388);
    const auto *lf_389 = buffer.data(lf + 389);
    const auto *lf_390 = buffer.data(lf + 390);
    const auto *lf_391 = buffer.data(lf + 391);
    const auto *lf_392 = buffer.data(lf + 392);
    const auto *lf_393 = buffer.data(lf + 393);
    const auto *lf_394 = buffer.data(lf + 394);
    const auto *lf_395 = buffer.data(lf + 395);
    const auto *lf_396 = buffer.data(lf + 396);
    const auto *lf_397 = buffer.data(lf + 397);
    const auto *lf_398 = buffer.data(lf + 398);
    const auto *lf_399 = buffer.data(lf + 399);
    const auto *lf_400 = buffer.data(lf + 400);
    const auto *lf_401 = buffer.data(lf + 401);
    const auto *lf_402 = buffer.data(lf + 402);
    const auto *lf_403 = buffer.data(lf + 403);
    const auto *lf_404 = buffer.data(lf + 404);
    const auto *lf_405 = buffer.data(lf + 405);
    const auto *lf_406 = buffer.data(lf + 406);
    const auto *lf_407 = buffer.data(lf + 407);
    const auto *lf_408 = buffer.data(lf + 408);
    const auto *lf_409 = buffer.data(lf + 409);
    const auto *lf_410 = buffer.data(lf + 410);
    const auto *lf_411 = buffer.data(lf + 411);
    const auto *lf_412 = buffer.data(lf + 412);
    const auto *lf_413 = buffer.data(lf + 413);
    const auto *lf_414 = buffer.data(lf + 414);
    const auto *lf_415 = buffer.data(lf + 415);
    const auto *lf_416 = buffer.data(lf + 416);

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, if__110, if__111, if__112, \
                         if__113, if__114, lf_240, lf_241, lf_242, lf_243, \
                         lf_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -2.0 * if__110[k]
                   + f_0 * lf_240[k];

        t_171[k] = -2.0 * if__111[k]
                   + f_0 * lf_241[k];

        t_172[k] = -2.0 * if__112[k]
                   + f_0 * lf_242[k];

        t_173[k] = -2.0 * if__113[k]
                   + f_0 * lf_243[k];

        t_174[k] = -2.0 * if__114[k]
                   + f_0 * lf_244[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, if__115, if__116, if__117, \
                         if__118, if__119, lf_245, lf_246, lf_247, lf_248, \
                         lf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -2.0 * if__115[k]
                   + f_0 * lf_245[k];

        t_176[k] = -2.0 * if__116[k]
                   + f_0 * lf_246[k];

        t_177[k] = -2.0 * if__117[k]
                   + f_0 * lf_247[k];

        t_178[k] = -2.0 * if__118[k]
                   + f_0 * lf_248[k];

        t_179[k] = -2.0 * if__119[k]
                   + f_0 * lf_249[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, if__120, if__121, if__122, \
                         if__123, if__124, lf_250, lf_251, lf_252, lf_253, \
                         lf_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * if__120[k]
                   + f_0 * lf_250[k];

        t_181[k] = -3.0 * if__121[k]
                   + f_0 * lf_251[k];

        t_182[k] = -3.0 * if__122[k]
                   + f_0 * lf_252[k];

        t_183[k] = -3.0 * if__123[k]
                   + f_0 * lf_253[k];

        t_184[k] = -3.0 * if__124[k]
                   + f_0 * lf_254[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, if__125, if__126, if__127, \
                         if__128, if__129, lf_255, lf_256, lf_257, lf_258, \
                         lf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * if__125[k]
                   + f_0 * lf_255[k];

        t_186[k] = -3.0 * if__126[k]
                   + f_0 * lf_256[k];

        t_187[k] = -3.0 * if__127[k]
                   + f_0 * lf_257[k];

        t_188[k] = -3.0 * if__128[k]
                   + f_0 * lf_258[k];

        t_189[k] = -3.0 * if__129[k]
                   + f_0 * lf_259[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, if__130, if__131, if__132, \
                         if__133, if__134, lf_260, lf_261, lf_262, lf_263, \
                         lf_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -4.0 * if__130[k]
                   + f_0 * lf_260[k];

        t_191[k] = -4.0 * if__131[k]
                   + f_0 * lf_261[k];

        t_192[k] = -4.0 * if__132[k]
                   + f_0 * lf_262[k];

        t_193[k] = -4.0 * if__133[k]
                   + f_0 * lf_263[k];

        t_194[k] = -4.0 * if__134[k]
                   + f_0 * lf_264[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, if__135, if__136, if__137, \
                         if__138, if__139, lf_265, lf_266, lf_267, lf_268, \
                         lf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -4.0 * if__135[k]
                   + f_0 * lf_265[k];

        t_196[k] = -4.0 * if__136[k]
                   + f_0 * lf_266[k];

        t_197[k] = -4.0 * if__137[k]
                   + f_0 * lf_267[k];

        t_198[k] = -4.0 * if__138[k]
                   + f_0 * lf_268[k];

        t_199[k] = -4.0 * if__139[k]
                   + f_0 * lf_269[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, if__140, if__141, if__142, \
                         if__143, if__144, lf_270, lf_271, lf_272, lf_273, \
                         lf_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -5.0 * if__140[k]
                   + f_0 * lf_270[k];

        t_201[k] = -5.0 * if__141[k]
                   + f_0 * lf_271[k];

        t_202[k] = -5.0 * if__142[k]
                   + f_0 * lf_272[k];

        t_203[k] = -5.0 * if__143[k]
                   + f_0 * lf_273[k];

        t_204[k] = -5.0 * if__144[k]
                   + f_0 * lf_274[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, if__145, if__146, if__147, \
                         if__148, if__149, lf_275, lf_276, lf_277, lf_278, \
                         lf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -5.0 * if__145[k]
                   + f_0 * lf_275[k];

        t_206[k] = -5.0 * if__146[k]
                   + f_0 * lf_276[k];

        t_207[k] = -5.0 * if__147[k]
                   + f_0 * lf_277[k];

        t_208[k] = -5.0 * if__148[k]
                   + f_0 * lf_278[k];

        t_209[k] = -5.0 * if__149[k]
                   + f_0 * lf_279[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, t_217, lf_290, \
                         lf_291, lf_292, lf_293, lf_294, lf_295, lf_296, \
                         lf_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * lf_290[k];

        t_211[k] = f_0 * lf_291[k];

        t_212[k] = f_0 * lf_292[k];

        t_213[k] = f_0 * lf_293[k];

        t_214[k] = f_0 * lf_294[k];

        t_215[k] = f_0 * lf_295[k];

        t_216[k] = f_0 * lf_296[k];

        t_217[k] = f_0 * lf_297[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, if__150, if__151, if__152, \
                         if__153, lf_298, lf_299, lf_300, lf_301, lf_302, \
                         lf_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * lf_298[k];

        t_219[k] = f_0 * lf_299[k];

        t_220[k] = -if__150[k]
                   + f_0 * lf_300[k];

        t_221[k] = -if__151[k]
                   + f_0 * lf_301[k];

        t_222[k] = -if__152[k]
                   + f_0 * lf_302[k];

        t_223[k] = -if__153[k]
                   + f_0 * lf_303[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, if__154, if__155, if__156, \
                         if__157, if__158, lf_304, lf_305, lf_306, lf_307, \
                         lf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = -if__154[k]
                   + f_0 * lf_304[k];

        t_225[k] = -if__155[k]
                   + f_0 * lf_305[k];

        t_226[k] = -if__156[k]
                   + f_0 * lf_306[k];

        t_227[k] = -if__157[k]
                   + f_0 * lf_307[k];

        t_228[k] = -if__158[k]
                   + f_0 * lf_308[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, if__159, if__160, if__161, \
                         if__162, if__163, lf_309, lf_310, lf_311, lf_312, \
                         lf_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = -if__159[k]
                   + f_0 * lf_309[k];

        t_230[k] = -2.0 * if__160[k]
                   + f_0 * lf_310[k];

        t_231[k] = -2.0 * if__161[k]
                   + f_0 * lf_311[k];

        t_232[k] = -2.0 * if__162[k]
                   + f_0 * lf_312[k];

        t_233[k] = -2.0 * if__163[k]
                   + f_0 * lf_313[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, if__164, if__165, if__166, \
                         if__167, if__168, lf_314, lf_315, lf_316, lf_317, \
                         lf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = -2.0 * if__164[k]
                   + f_0 * lf_314[k];

        t_235[k] = -2.0 * if__165[k]
                   + f_0 * lf_315[k];

        t_236[k] = -2.0 * if__166[k]
                   + f_0 * lf_316[k];

        t_237[k] = -2.0 * if__167[k]
                   + f_0 * lf_317[k];

        t_238[k] = -2.0 * if__168[k]
                   + f_0 * lf_318[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, if__169, if__170, if__171, \
                         if__172, if__173, lf_319, lf_320, lf_321, lf_322, \
                         lf_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = -2.0 * if__169[k]
                   + f_0 * lf_319[k];

        t_240[k] = -3.0 * if__170[k]
                   + f_0 * lf_320[k];

        t_241[k] = -3.0 * if__171[k]
                   + f_0 * lf_321[k];

        t_242[k] = -3.0 * if__172[k]
                   + f_0 * lf_322[k];

        t_243[k] = -3.0 * if__173[k]
                   + f_0 * lf_323[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, t_248, if__174, if__175, if__176, \
                         if__177, if__178, lf_324, lf_325, lf_326, lf_327, \
                         lf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = -3.0 * if__174[k]
                   + f_0 * lf_324[k];

        t_245[k] = -3.0 * if__175[k]
                   + f_0 * lf_325[k];

        t_246[k] = -3.0 * if__176[k]
                   + f_0 * lf_326[k];

        t_247[k] = -3.0 * if__177[k]
                   + f_0 * lf_327[k];

        t_248[k] = -3.0 * if__178[k]
                   + f_0 * lf_328[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, t_253, if__179, if__180, if__181, \
                         if__182, if__183, lf_329, lf_330, lf_331, lf_332, \
                         lf_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = -3.0 * if__179[k]
                   + f_0 * lf_329[k];

        t_250[k] = -4.0 * if__180[k]
                   + f_0 * lf_330[k];

        t_251[k] = -4.0 * if__181[k]
                   + f_0 * lf_331[k];

        t_252[k] = -4.0 * if__182[k]
                   + f_0 * lf_332[k];

        t_253[k] = -4.0 * if__183[k]
                   + f_0 * lf_333[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, t_258, if__184, if__185, if__186, \
                         if__187, if__188, lf_334, lf_335, lf_336, lf_337, \
                         lf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = -4.0 * if__184[k]
                   + f_0 * lf_334[k];

        t_255[k] = -4.0 * if__185[k]
                   + f_0 * lf_335[k];

        t_256[k] = -4.0 * if__186[k]
                   + f_0 * lf_336[k];

        t_257[k] = -4.0 * if__187[k]
                   + f_0 * lf_337[k];

        t_258[k] = -4.0 * if__188[k]
                   + f_0 * lf_338[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, if__189, if__190, if__191, \
                         if__192, if__193, lf_339, lf_340, lf_341, lf_342, \
                         lf_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = -4.0 * if__189[k]
                   + f_0 * lf_339[k];

        t_260[k] = -5.0 * if__190[k]
                   + f_0 * lf_340[k];

        t_261[k] = -5.0 * if__191[k]
                   + f_0 * lf_341[k];

        t_262[k] = -5.0 * if__192[k]
                   + f_0 * lf_342[k];

        t_263[k] = -5.0 * if__193[k]
                   + f_0 * lf_343[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, if__194, if__195, if__196, \
                         if__197, if__198, lf_344, lf_345, lf_346, lf_347, \
                         lf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = -5.0 * if__194[k]
                   + f_0 * lf_344[k];

        t_265[k] = -5.0 * if__195[k]
                   + f_0 * lf_345[k];

        t_266[k] = -5.0 * if__196[k]
                   + f_0 * lf_346[k];

        t_267[k] = -5.0 * if__197[k]
                   + f_0 * lf_347[k];

        t_268[k] = -5.0 * if__198[k]
                   + f_0 * lf_348[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, t_273, if__199, if__200, if__201, \
                         if__202, if__203, lf_349, lf_350, lf_351, lf_352, \
                         lf_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = -5.0 * if__199[k]
                   + f_0 * lf_349[k];

        t_270[k] = -6.0 * if__200[k]
                   + f_0 * lf_350[k];

        t_271[k] = -6.0 * if__201[k]
                   + f_0 * lf_351[k];

        t_272[k] = -6.0 * if__202[k]
                   + f_0 * lf_352[k];

        t_273[k] = -6.0 * if__203[k]
                   + f_0 * lf_353[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, if__204, if__205, if__206, \
                         if__207, if__208, lf_354, lf_355, lf_356, lf_357, \
                         lf_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = -6.0 * if__204[k]
                   + f_0 * lf_354[k];

        t_275[k] = -6.0 * if__205[k]
                   + f_0 * lf_355[k];

        t_276[k] = -6.0 * if__206[k]
                   + f_0 * lf_356[k];

        t_277[k] = -6.0 * if__207[k]
                   + f_0 * lf_357[k];

        t_278[k] = -6.0 * if__208[k]
                   + f_0 * lf_358[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, t_283, t_284, t_285, if__209, lf_359, \
                         lf_370, lf_371, lf_372, lf_373, lf_374, \
                         lf_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = -6.0 * if__209[k]
                   + f_0 * lf_359[k];

        t_280[k] = f_0 * lf_370[k];

        t_281[k] = f_0 * lf_371[k];

        t_282[k] = f_0 * lf_372[k];

        t_283[k] = f_0 * lf_373[k];

        t_284[k] = f_0 * lf_374[k];

        t_285[k] = f_0 * lf_375[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, t_291, if__210, if__211, lf_376, \
                         lf_377, lf_378, lf_379, lf_380, lf_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * lf_376[k];

        t_287[k] = f_0 * lf_377[k];

        t_288[k] = f_0 * lf_378[k];

        t_289[k] = f_0 * lf_379[k];

        t_290[k] = -if__210[k]
                   + f_0 * lf_380[k];

        t_291[k] = -if__211[k]
                   + f_0 * lf_381[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, if__212, if__213, if__214, \
                         if__215, if__216, lf_382, lf_383, lf_384, lf_385, \
                         lf_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = -if__212[k]
                   + f_0 * lf_382[k];

        t_293[k] = -if__213[k]
                   + f_0 * lf_383[k];

        t_294[k] = -if__214[k]
                   + f_0 * lf_384[k];

        t_295[k] = -if__215[k]
                   + f_0 * lf_385[k];

        t_296[k] = -if__216[k]
                   + f_0 * lf_386[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, if__217, if__218, if__219, \
                         if__220, if__221, lf_387, lf_388, lf_389, lf_390, \
                         lf_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = -if__217[k]
                   + f_0 * lf_387[k];

        t_298[k] = -if__218[k]
                   + f_0 * lf_388[k];

        t_299[k] = -if__219[k]
                   + f_0 * lf_389[k];

        t_300[k] = -2.0 * if__220[k]
                   + f_0 * lf_390[k];

        t_301[k] = -2.0 * if__221[k]
                   + f_0 * lf_391[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, if__222, if__223, if__224, \
                         if__225, if__226, lf_392, lf_393, lf_394, lf_395, \
                         lf_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = -2.0 * if__222[k]
                   + f_0 * lf_392[k];

        t_303[k] = -2.0 * if__223[k]
                   + f_0 * lf_393[k];

        t_304[k] = -2.0 * if__224[k]
                   + f_0 * lf_394[k];

        t_305[k] = -2.0 * if__225[k]
                   + f_0 * lf_395[k];

        t_306[k] = -2.0 * if__226[k]
                   + f_0 * lf_396[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, if__227, if__228, if__229, \
                         if__230, if__231, lf_397, lf_398, lf_399, lf_400, \
                         lf_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = -2.0 * if__227[k]
                   + f_0 * lf_397[k];

        t_308[k] = -2.0 * if__228[k]
                   + f_0 * lf_398[k];

        t_309[k] = -2.0 * if__229[k]
                   + f_0 * lf_399[k];

        t_310[k] = -3.0 * if__230[k]
                   + f_0 * lf_400[k];

        t_311[k] = -3.0 * if__231[k]
                   + f_0 * lf_401[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, if__232, if__233, if__234, \
                         if__235, if__236, lf_402, lf_403, lf_404, lf_405, \
                         lf_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -3.0 * if__232[k]
                   + f_0 * lf_402[k];

        t_313[k] = -3.0 * if__233[k]
                   + f_0 * lf_403[k];

        t_314[k] = -3.0 * if__234[k]
                   + f_0 * lf_404[k];

        t_315[k] = -3.0 * if__235[k]
                   + f_0 * lf_405[k];

        t_316[k] = -3.0 * if__236[k]
                   + f_0 * lf_406[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, if__237, if__238, if__239, \
                         if__240, if__241, lf_407, lf_408, lf_409, lf_410, \
                         lf_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = -3.0 * if__237[k]
                   + f_0 * lf_407[k];

        t_318[k] = -3.0 * if__238[k]
                   + f_0 * lf_408[k];

        t_319[k] = -3.0 * if__239[k]
                   + f_0 * lf_409[k];

        t_320[k] = -4.0 * if__240[k]
                   + f_0 * lf_410[k];

        t_321[k] = -4.0 * if__241[k]
                   + f_0 * lf_411[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, if__242, if__243, if__244, \
                         if__245, if__246, lf_412, lf_413, lf_414, lf_415, \
                         lf_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = -4.0 * if__242[k]
                   + f_0 * lf_412[k];

        t_323[k] = -4.0 * if__243[k]
                   + f_0 * lf_413[k];

        t_324[k] = -4.0 * if__244[k]
                   + f_0 * lf_414[k];

        t_325[k] = -4.0 * if__245[k]
                   + f_0 * lf_415[k];

        t_326[k] = -4.0 * if__246[k]
                   + f_0 * lf_416[k];
    }
}

static auto
compute_prim_geom_10_kf_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t if_, const size_t lf,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *lf_417 = buffer.data(lf + 417);
    const auto *lf_418 = buffer.data(lf + 418);
    const auto *lf_419 = buffer.data(lf + 419);
    const auto *lf_420 = buffer.data(lf + 420);
    const auto *lf_421 = buffer.data(lf + 421);
    const auto *lf_422 = buffer.data(lf + 422);
    const auto *lf_423 = buffer.data(lf + 423);
    const auto *lf_424 = buffer.data(lf + 424);
    const auto *lf_425 = buffer.data(lf + 425);
    const auto *lf_426 = buffer.data(lf + 426);
    const auto *lf_427 = buffer.data(lf + 427);
    const auto *lf_428 = buffer.data(lf + 428);
    const auto *lf_429 = buffer.data(lf + 429);
    const auto *lf_430 = buffer.data(lf + 430);
    const auto *lf_431 = buffer.data(lf + 431);
    const auto *lf_432 = buffer.data(lf + 432);
    const auto *lf_433 = buffer.data(lf + 433);
    const auto *lf_434 = buffer.data(lf + 434);
    const auto *lf_435 = buffer.data(lf + 435);
    const auto *lf_436 = buffer.data(lf + 436);
    const auto *lf_437 = buffer.data(lf + 437);
    const auto *lf_438 = buffer.data(lf + 438);
    const auto *lf_439 = buffer.data(lf + 439);
    const auto *lf_440 = buffer.data(lf + 440);
    const auto *lf_441 = buffer.data(lf + 441);
    const auto *lf_442 = buffer.data(lf + 442);
    const auto *lf_443 = buffer.data(lf + 443);
    const auto *lf_444 = buffer.data(lf + 444);
    const auto *lf_445 = buffer.data(lf + 445);
    const auto *lf_446 = buffer.data(lf + 446);
    const auto *lf_447 = buffer.data(lf + 447);
    const auto *lf_448 = buffer.data(lf + 448);
    const auto *lf_449 = buffer.data(lf + 449);

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, if__247, if__248, if__249, \
                         if__250, if__251, lf_417, lf_418, lf_419, lf_420, \
                         lf_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = -4.0 * if__247[k]
                   + f_0 * lf_417[k];

        t_328[k] = -4.0 * if__248[k]
                   + f_0 * lf_418[k];

        t_329[k] = -4.0 * if__249[k]
                   + f_0 * lf_419[k];

        t_330[k] = -5.0 * if__250[k]
                   + f_0 * lf_420[k];

        t_331[k] = -5.0 * if__251[k]
                   + f_0 * lf_421[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, if__252, if__253, if__254, \
                         if__255, if__256, lf_422, lf_423, lf_424, lf_425, \
                         lf_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -5.0 * if__252[k]
                   + f_0 * lf_422[k];

        t_333[k] = -5.0 * if__253[k]
                   + f_0 * lf_423[k];

        t_334[k] = -5.0 * if__254[k]
                   + f_0 * lf_424[k];

        t_335[k] = -5.0 * if__255[k]
                   + f_0 * lf_425[k];

        t_336[k] = -5.0 * if__256[k]
                   + f_0 * lf_426[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, if__257, if__258, if__259, \
                         if__260, if__261, lf_427, lf_428, lf_429, lf_430, \
                         lf_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -5.0 * if__257[k]
                   + f_0 * lf_427[k];

        t_338[k] = -5.0 * if__258[k]
                   + f_0 * lf_428[k];

        t_339[k] = -5.0 * if__259[k]
                   + f_0 * lf_429[k];

        t_340[k] = -6.0 * if__260[k]
                   + f_0 * lf_430[k];

        t_341[k] = -6.0 * if__261[k]
                   + f_0 * lf_431[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, if__262, if__263, if__264, \
                         if__265, if__266, lf_432, lf_433, lf_434, lf_435, \
                         lf_436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -6.0 * if__262[k]
                   + f_0 * lf_432[k];

        t_343[k] = -6.0 * if__263[k]
                   + f_0 * lf_433[k];

        t_344[k] = -6.0 * if__264[k]
                   + f_0 * lf_434[k];

        t_345[k] = -6.0 * if__265[k]
                   + f_0 * lf_435[k];

        t_346[k] = -6.0 * if__266[k]
                   + f_0 * lf_436[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, if__267, if__268, if__269, \
                         if__270, if__271, lf_437, lf_438, lf_439, lf_440, \
                         lf_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -6.0 * if__267[k]
                   + f_0 * lf_437[k];

        t_348[k] = -6.0 * if__268[k]
                   + f_0 * lf_438[k];

        t_349[k] = -6.0 * if__269[k]
                   + f_0 * lf_439[k];

        t_350[k] = -7.0 * if__270[k]
                   + f_0 * lf_440[k];

        t_351[k] = -7.0 * if__271[k]
                   + f_0 * lf_441[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, if__272, if__273, if__274, \
                         if__275, if__276, lf_442, lf_443, lf_444, lf_445, \
                         lf_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -7.0 * if__272[k]
                   + f_0 * lf_442[k];

        t_353[k] = -7.0 * if__273[k]
                   + f_0 * lf_443[k];

        t_354[k] = -7.0 * if__274[k]
                   + f_0 * lf_444[k];

        t_355[k] = -7.0 * if__275[k]
                   + f_0 * lf_445[k];

        t_356[k] = -7.0 * if__276[k]
                   + f_0 * lf_446[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, if__277, if__278, if__279, lf_447, lf_448, \
                         lf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -7.0 * if__277[k]
                   + f_0 * lf_447[k];

        t_358[k] = -7.0 * if__278[k]
                   + f_0 * lf_448[k];

        t_359[k] = -7.0 * if__279[k]
                   + f_0 * lf_449[k];
    }
}

auto
compute_prim_geom_10_kf_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t if_, const size_t lf,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_kf_electron_repulsion_2_piece0(buffer, target, if_, lf, ncols, alpha);

    compute_prim_geom_10_kf_electron_repulsion_2_piece1(buffer, target, if_, lf, ncols, alpha);

    compute_prim_geom_10_kf_electron_repulsion_2_piece2(buffer, target, if_, lf, ncols, alpha);
}

}  // namespace simdt2ceri
