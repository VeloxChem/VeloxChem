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


#include "SimdElectronRepulsionGeom10VrrRecLF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_lf_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kf, const size_t mf,
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

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_1 = buffer.data(kf + 1);
    const auto *kf_2 = buffer.data(kf + 2);
    const auto *kf_3 = buffer.data(kf + 3);
    const auto *kf_4 = buffer.data(kf + 4);
    const auto *kf_5 = buffer.data(kf + 5);
    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_12 = buffer.data(kf + 12);
    const auto *kf_13 = buffer.data(kf + 13);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_15 = buffer.data(kf + 15);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_17 = buffer.data(kf + 17);
    const auto *kf_18 = buffer.data(kf + 18);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_21 = buffer.data(kf + 21);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_24 = buffer.data(kf + 24);
    const auto *kf_25 = buffer.data(kf + 25);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_37 = buffer.data(kf + 37);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_49 = buffer.data(kf + 49);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_51 = buffer.data(kf + 51);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_53 = buffer.data(kf + 53);
    const auto *kf_54 = buffer.data(kf + 54);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_62 = buffer.data(kf + 62);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_65 = buffer.data(kf + 65);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_73 = buffer.data(kf + 73);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_81 = buffer.data(kf + 81);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_85 = buffer.data(kf + 85);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_93 = buffer.data(kf + 93);
    const auto *kf_94 = buffer.data(kf + 94);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_97 = buffer.data(kf + 97);
    const auto *kf_98 = buffer.data(kf + 98);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_101 = buffer.data(kf + 101);
    const auto *kf_102 = buffer.data(kf + 102);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_104 = buffer.data(kf + 104);
    const auto *kf_105 = buffer.data(kf + 105);
    const auto *kf_106 = buffer.data(kf + 106);
    const auto *kf_107 = buffer.data(kf + 107);
    const auto *kf_108 = buffer.data(kf + 108);
    const auto *kf_109 = buffer.data(kf + 109);
    const auto *kf_110 = buffer.data(kf + 110);
    const auto *kf_111 = buffer.data(kf + 111);
    const auto *kf_112 = buffer.data(kf + 112);
    const auto *kf_113 = buffer.data(kf + 113);
    const auto *kf_114 = buffer.data(kf + 114);
    const auto *kf_115 = buffer.data(kf + 115);
    const auto *kf_116 = buffer.data(kf + 116);
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_119 = buffer.data(kf + 119);
    const auto *kf_120 = buffer.data(kf + 120);
    const auto *kf_121 = buffer.data(kf + 121);
    const auto *kf_122 = buffer.data(kf + 122);
    const auto *kf_123 = buffer.data(kf + 123);
    const auto *kf_124 = buffer.data(kf + 124);
    const auto *kf_125 = buffer.data(kf + 125);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_127 = buffer.data(kf + 127);
    const auto *kf_128 = buffer.data(kf + 128);
    const auto *kf_129 = buffer.data(kf + 129);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_131 = buffer.data(kf + 131);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_133 = buffer.data(kf + 133);
    const auto *kf_134 = buffer.data(kf + 134);
    const auto *kf_135 = buffer.data(kf + 135);
    const auto *kf_136 = buffer.data(kf + 136);
    const auto *kf_137 = buffer.data(kf + 137);
    const auto *kf_138 = buffer.data(kf + 138);
    const auto *kf_139 = buffer.data(kf + 139);
    const auto *kf_140 = buffer.data(kf + 140);
    const auto *kf_141 = buffer.data(kf + 141);
    const auto *kf_142 = buffer.data(kf + 142);
    const auto *kf_143 = buffer.data(kf + 143);
    const auto *kf_144 = buffer.data(kf + 144);
    const auto *kf_145 = buffer.data(kf + 145);
    const auto *kf_146 = buffer.data(kf + 146);
    const auto *kf_147 = buffer.data(kf + 147);
    const auto *kf_148 = buffer.data(kf + 148);
    const auto *kf_149 = buffer.data(kf + 149);

    const auto *mf_0 = buffer.data(mf + 0);
    const auto *mf_1 = buffer.data(mf + 1);
    const auto *mf_2 = buffer.data(mf + 2);
    const auto *mf_3 = buffer.data(mf + 3);
    const auto *mf_4 = buffer.data(mf + 4);
    const auto *mf_5 = buffer.data(mf + 5);
    const auto *mf_6 = buffer.data(mf + 6);
    const auto *mf_7 = buffer.data(mf + 7);
    const auto *mf_8 = buffer.data(mf + 8);
    const auto *mf_9 = buffer.data(mf + 9);
    const auto *mf_10 = buffer.data(mf + 10);
    const auto *mf_11 = buffer.data(mf + 11);
    const auto *mf_12 = buffer.data(mf + 12);
    const auto *mf_13 = buffer.data(mf + 13);
    const auto *mf_14 = buffer.data(mf + 14);
    const auto *mf_15 = buffer.data(mf + 15);
    const auto *mf_16 = buffer.data(mf + 16);
    const auto *mf_17 = buffer.data(mf + 17);
    const auto *mf_18 = buffer.data(mf + 18);
    const auto *mf_19 = buffer.data(mf + 19);
    const auto *mf_20 = buffer.data(mf + 20);
    const auto *mf_21 = buffer.data(mf + 21);
    const auto *mf_22 = buffer.data(mf + 22);
    const auto *mf_23 = buffer.data(mf + 23);
    const auto *mf_24 = buffer.data(mf + 24);
    const auto *mf_25 = buffer.data(mf + 25);
    const auto *mf_26 = buffer.data(mf + 26);
    const auto *mf_27 = buffer.data(mf + 27);
    const auto *mf_28 = buffer.data(mf + 28);
    const auto *mf_29 = buffer.data(mf + 29);
    const auto *mf_30 = buffer.data(mf + 30);
    const auto *mf_31 = buffer.data(mf + 31);
    const auto *mf_32 = buffer.data(mf + 32);
    const auto *mf_33 = buffer.data(mf + 33);
    const auto *mf_34 = buffer.data(mf + 34);
    const auto *mf_35 = buffer.data(mf + 35);
    const auto *mf_36 = buffer.data(mf + 36);
    const auto *mf_37 = buffer.data(mf + 37);
    const auto *mf_38 = buffer.data(mf + 38);
    const auto *mf_39 = buffer.data(mf + 39);
    const auto *mf_40 = buffer.data(mf + 40);
    const auto *mf_41 = buffer.data(mf + 41);
    const auto *mf_42 = buffer.data(mf + 42);
    const auto *mf_43 = buffer.data(mf + 43);
    const auto *mf_44 = buffer.data(mf + 44);
    const auto *mf_45 = buffer.data(mf + 45);
    const auto *mf_46 = buffer.data(mf + 46);
    const auto *mf_47 = buffer.data(mf + 47);
    const auto *mf_48 = buffer.data(mf + 48);
    const auto *mf_49 = buffer.data(mf + 49);
    const auto *mf_50 = buffer.data(mf + 50);
    const auto *mf_51 = buffer.data(mf + 51);
    const auto *mf_52 = buffer.data(mf + 52);
    const auto *mf_53 = buffer.data(mf + 53);
    const auto *mf_54 = buffer.data(mf + 54);
    const auto *mf_55 = buffer.data(mf + 55);
    const auto *mf_56 = buffer.data(mf + 56);
    const auto *mf_57 = buffer.data(mf + 57);
    const auto *mf_58 = buffer.data(mf + 58);
    const auto *mf_59 = buffer.data(mf + 59);
    const auto *mf_60 = buffer.data(mf + 60);
    const auto *mf_61 = buffer.data(mf + 61);
    const auto *mf_62 = buffer.data(mf + 62);
    const auto *mf_63 = buffer.data(mf + 63);
    const auto *mf_64 = buffer.data(mf + 64);
    const auto *mf_65 = buffer.data(mf + 65);
    const auto *mf_66 = buffer.data(mf + 66);
    const auto *mf_67 = buffer.data(mf + 67);
    const auto *mf_68 = buffer.data(mf + 68);
    const auto *mf_69 = buffer.data(mf + 69);
    const auto *mf_70 = buffer.data(mf + 70);
    const auto *mf_71 = buffer.data(mf + 71);
    const auto *mf_72 = buffer.data(mf + 72);
    const auto *mf_73 = buffer.data(mf + 73);
    const auto *mf_74 = buffer.data(mf + 74);
    const auto *mf_75 = buffer.data(mf + 75);
    const auto *mf_76 = buffer.data(mf + 76);
    const auto *mf_77 = buffer.data(mf + 77);
    const auto *mf_78 = buffer.data(mf + 78);
    const auto *mf_79 = buffer.data(mf + 79);
    const auto *mf_80 = buffer.data(mf + 80);
    const auto *mf_81 = buffer.data(mf + 81);
    const auto *mf_82 = buffer.data(mf + 82);
    const auto *mf_83 = buffer.data(mf + 83);
    const auto *mf_84 = buffer.data(mf + 84);
    const auto *mf_85 = buffer.data(mf + 85);
    const auto *mf_86 = buffer.data(mf + 86);
    const auto *mf_87 = buffer.data(mf + 87);
    const auto *mf_88 = buffer.data(mf + 88);
    const auto *mf_89 = buffer.data(mf + 89);
    const auto *mf_90 = buffer.data(mf + 90);
    const auto *mf_91 = buffer.data(mf + 91);
    const auto *mf_92 = buffer.data(mf + 92);
    const auto *mf_93 = buffer.data(mf + 93);
    const auto *mf_94 = buffer.data(mf + 94);
    const auto *mf_95 = buffer.data(mf + 95);
    const auto *mf_96 = buffer.data(mf + 96);
    const auto *mf_97 = buffer.data(mf + 97);
    const auto *mf_98 = buffer.data(mf + 98);
    const auto *mf_99 = buffer.data(mf + 99);
    const auto *mf_100 = buffer.data(mf + 100);
    const auto *mf_101 = buffer.data(mf + 101);
    const auto *mf_102 = buffer.data(mf + 102);
    const auto *mf_103 = buffer.data(mf + 103);
    const auto *mf_104 = buffer.data(mf + 104);
    const auto *mf_105 = buffer.data(mf + 105);
    const auto *mf_106 = buffer.data(mf + 106);
    const auto *mf_107 = buffer.data(mf + 107);
    const auto *mf_108 = buffer.data(mf + 108);
    const auto *mf_109 = buffer.data(mf + 109);
    const auto *mf_110 = buffer.data(mf + 110);
    const auto *mf_111 = buffer.data(mf + 111);
    const auto *mf_112 = buffer.data(mf + 112);
    const auto *mf_113 = buffer.data(mf + 113);
    const auto *mf_114 = buffer.data(mf + 114);
    const auto *mf_115 = buffer.data(mf + 115);
    const auto *mf_116 = buffer.data(mf + 116);
    const auto *mf_117 = buffer.data(mf + 117);
    const auto *mf_118 = buffer.data(mf + 118);
    const auto *mf_119 = buffer.data(mf + 119);
    const auto *mf_120 = buffer.data(mf + 120);
    const auto *mf_121 = buffer.data(mf + 121);
    const auto *mf_122 = buffer.data(mf + 122);
    const auto *mf_123 = buffer.data(mf + 123);
    const auto *mf_124 = buffer.data(mf + 124);
    const auto *mf_125 = buffer.data(mf + 125);
    const auto *mf_126 = buffer.data(mf + 126);
    const auto *mf_127 = buffer.data(mf + 127);
    const auto *mf_128 = buffer.data(mf + 128);
    const auto *mf_129 = buffer.data(mf + 129);
    const auto *mf_130 = buffer.data(mf + 130);
    const auto *mf_131 = buffer.data(mf + 131);
    const auto *mf_132 = buffer.data(mf + 132);
    const auto *mf_133 = buffer.data(mf + 133);
    const auto *mf_134 = buffer.data(mf + 134);
    const auto *mf_135 = buffer.data(mf + 135);
    const auto *mf_136 = buffer.data(mf + 136);
    const auto *mf_137 = buffer.data(mf + 137);
    const auto *mf_138 = buffer.data(mf + 138);
    const auto *mf_139 = buffer.data(mf + 139);
    const auto *mf_140 = buffer.data(mf + 140);
    const auto *mf_141 = buffer.data(mf + 141);
    const auto *mf_142 = buffer.data(mf + 142);
    const auto *mf_143 = buffer.data(mf + 143);
    const auto *mf_144 = buffer.data(mf + 144);
    const auto *mf_145 = buffer.data(mf + 145);
    const auto *mf_146 = buffer.data(mf + 146);
    const auto *mf_147 = buffer.data(mf + 147);
    const auto *mf_148 = buffer.data(mf + 148);
    const auto *mf_149 = buffer.data(mf + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, kf_0, kf_1, kf_2, kf_3, kf_4, mf_0, mf_1, \
                         mf_2, mf_3, mf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -8.0 * kf_0[k]
                 + f_0 * mf_0[k];

        t_1[k] = -8.0 * kf_1[k]
                 + f_0 * mf_1[k];

        t_2[k] = -8.0 * kf_2[k]
                 + f_0 * mf_2[k];

        t_3[k] = -8.0 * kf_3[k]
                 + f_0 * mf_3[k];

        t_4[k] = -8.0 * kf_4[k]
                 + f_0 * mf_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, kf_5, kf_6, kf_7, kf_8, kf_9, mf_5, mf_6, \
                         mf_7, mf_8, mf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -8.0 * kf_5[k]
                 + f_0 * mf_5[k];

        t_6[k] = -8.0 * kf_6[k]
                 + f_0 * mf_6[k];

        t_7[k] = -8.0 * kf_7[k]
                 + f_0 * mf_7[k];

        t_8[k] = -8.0 * kf_8[k]
                 + f_0 * mf_8[k];

        t_9[k] = -8.0 * kf_9[k]
                 + f_0 * mf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, kf_10, kf_11, kf_12, kf_13, kf_14, \
                         mf_10, mf_11, mf_12, mf_13, mf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -7.0 * kf_10[k]
                  + f_0 * mf_10[k];

        t_11[k] = -7.0 * kf_11[k]
                  + f_0 * mf_11[k];

        t_12[k] = -7.0 * kf_12[k]
                  + f_0 * mf_12[k];

        t_13[k] = -7.0 * kf_13[k]
                  + f_0 * mf_13[k];

        t_14[k] = -7.0 * kf_14[k]
                  + f_0 * mf_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, kf_15, kf_16, kf_17, kf_18, kf_19, \
                         mf_15, mf_16, mf_17, mf_18, mf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -7.0 * kf_15[k]
                  + f_0 * mf_15[k];

        t_16[k] = -7.0 * kf_16[k]
                  + f_0 * mf_16[k];

        t_17[k] = -7.0 * kf_17[k]
                  + f_0 * mf_17[k];

        t_18[k] = -7.0 * kf_18[k]
                  + f_0 * mf_18[k];

        t_19[k] = -7.0 * kf_19[k]
                  + f_0 * mf_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, kf_20, kf_21, kf_22, kf_23, kf_24, \
                         mf_20, mf_21, mf_22, mf_23, mf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -7.0 * kf_20[k]
                  + f_0 * mf_20[k];

        t_21[k] = -7.0 * kf_21[k]
                  + f_0 * mf_21[k];

        t_22[k] = -7.0 * kf_22[k]
                  + f_0 * mf_22[k];

        t_23[k] = -7.0 * kf_23[k]
                  + f_0 * mf_23[k];

        t_24[k] = -7.0 * kf_24[k]
                  + f_0 * mf_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, kf_25, kf_26, kf_27, kf_28, kf_29, \
                         mf_25, mf_26, mf_27, mf_28, mf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -7.0 * kf_25[k]
                  + f_0 * mf_25[k];

        t_26[k] = -7.0 * kf_26[k]
                  + f_0 * mf_26[k];

        t_27[k] = -7.0 * kf_27[k]
                  + f_0 * mf_27[k];

        t_28[k] = -7.0 * kf_28[k]
                  + f_0 * mf_28[k];

        t_29[k] = -7.0 * kf_29[k]
                  + f_0 * mf_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, kf_30, kf_31, kf_32, kf_33, kf_34, \
                         mf_30, mf_31, mf_32, mf_33, mf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -6.0 * kf_30[k]
                  + f_0 * mf_30[k];

        t_31[k] = -6.0 * kf_31[k]
                  + f_0 * mf_31[k];

        t_32[k] = -6.0 * kf_32[k]
                  + f_0 * mf_32[k];

        t_33[k] = -6.0 * kf_33[k]
                  + f_0 * mf_33[k];

        t_34[k] = -6.0 * kf_34[k]
                  + f_0 * mf_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, kf_35, kf_36, kf_37, kf_38, kf_39, \
                         mf_35, mf_36, mf_37, mf_38, mf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -6.0 * kf_35[k]
                  + f_0 * mf_35[k];

        t_36[k] = -6.0 * kf_36[k]
                  + f_0 * mf_36[k];

        t_37[k] = -6.0 * kf_37[k]
                  + f_0 * mf_37[k];

        t_38[k] = -6.0 * kf_38[k]
                  + f_0 * mf_38[k];

        t_39[k] = -6.0 * kf_39[k]
                  + f_0 * mf_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, kf_40, kf_41, kf_42, kf_43, kf_44, \
                         mf_40, mf_41, mf_42, mf_43, mf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -6.0 * kf_40[k]
                  + f_0 * mf_40[k];

        t_41[k] = -6.0 * kf_41[k]
                  + f_0 * mf_41[k];

        t_42[k] = -6.0 * kf_42[k]
                  + f_0 * mf_42[k];

        t_43[k] = -6.0 * kf_43[k]
                  + f_0 * mf_43[k];

        t_44[k] = -6.0 * kf_44[k]
                  + f_0 * mf_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, kf_45, kf_46, kf_47, kf_48, kf_49, \
                         mf_45, mf_46, mf_47, mf_48, mf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -6.0 * kf_45[k]
                  + f_0 * mf_45[k];

        t_46[k] = -6.0 * kf_46[k]
                  + f_0 * mf_46[k];

        t_47[k] = -6.0 * kf_47[k]
                  + f_0 * mf_47[k];

        t_48[k] = -6.0 * kf_48[k]
                  + f_0 * mf_48[k];

        t_49[k] = -6.0 * kf_49[k]
                  + f_0 * mf_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, kf_50, kf_51, kf_52, kf_53, kf_54, \
                         mf_50, mf_51, mf_52, mf_53, mf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -6.0 * kf_50[k]
                  + f_0 * mf_50[k];

        t_51[k] = -6.0 * kf_51[k]
                  + f_0 * mf_51[k];

        t_52[k] = -6.0 * kf_52[k]
                  + f_0 * mf_52[k];

        t_53[k] = -6.0 * kf_53[k]
                  + f_0 * mf_53[k];

        t_54[k] = -6.0 * kf_54[k]
                  + f_0 * mf_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, kf_55, kf_56, kf_57, kf_58, kf_59, \
                         mf_55, mf_56, mf_57, mf_58, mf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -6.0 * kf_55[k]
                  + f_0 * mf_55[k];

        t_56[k] = -6.0 * kf_56[k]
                  + f_0 * mf_56[k];

        t_57[k] = -6.0 * kf_57[k]
                  + f_0 * mf_57[k];

        t_58[k] = -6.0 * kf_58[k]
                  + f_0 * mf_58[k];

        t_59[k] = -6.0 * kf_59[k]
                  + f_0 * mf_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, kf_60, kf_61, kf_62, kf_63, kf_64, \
                         mf_60, mf_61, mf_62, mf_63, mf_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -5.0 * kf_60[k]
                  + f_0 * mf_60[k];

        t_61[k] = -5.0 * kf_61[k]
                  + f_0 * mf_61[k];

        t_62[k] = -5.0 * kf_62[k]
                  + f_0 * mf_62[k];

        t_63[k] = -5.0 * kf_63[k]
                  + f_0 * mf_63[k];

        t_64[k] = -5.0 * kf_64[k]
                  + f_0 * mf_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, kf_65, kf_66, kf_67, kf_68, kf_69, \
                         mf_65, mf_66, mf_67, mf_68, mf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -5.0 * kf_65[k]
                  + f_0 * mf_65[k];

        t_66[k] = -5.0 * kf_66[k]
                  + f_0 * mf_66[k];

        t_67[k] = -5.0 * kf_67[k]
                  + f_0 * mf_67[k];

        t_68[k] = -5.0 * kf_68[k]
                  + f_0 * mf_68[k];

        t_69[k] = -5.0 * kf_69[k]
                  + f_0 * mf_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, kf_70, kf_71, kf_72, kf_73, kf_74, \
                         mf_70, mf_71, mf_72, mf_73, mf_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -5.0 * kf_70[k]
                  + f_0 * mf_70[k];

        t_71[k] = -5.0 * kf_71[k]
                  + f_0 * mf_71[k];

        t_72[k] = -5.0 * kf_72[k]
                  + f_0 * mf_72[k];

        t_73[k] = -5.0 * kf_73[k]
                  + f_0 * mf_73[k];

        t_74[k] = -5.0 * kf_74[k]
                  + f_0 * mf_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, kf_75, kf_76, kf_77, kf_78, kf_79, \
                         mf_75, mf_76, mf_77, mf_78, mf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -5.0 * kf_75[k]
                  + f_0 * mf_75[k];

        t_76[k] = -5.0 * kf_76[k]
                  + f_0 * mf_76[k];

        t_77[k] = -5.0 * kf_77[k]
                  + f_0 * mf_77[k];

        t_78[k] = -5.0 * kf_78[k]
                  + f_0 * mf_78[k];

        t_79[k] = -5.0 * kf_79[k]
                  + f_0 * mf_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, kf_80, kf_81, kf_82, kf_83, kf_84, \
                         mf_80, mf_81, mf_82, mf_83, mf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -5.0 * kf_80[k]
                  + f_0 * mf_80[k];

        t_81[k] = -5.0 * kf_81[k]
                  + f_0 * mf_81[k];

        t_82[k] = -5.0 * kf_82[k]
                  + f_0 * mf_82[k];

        t_83[k] = -5.0 * kf_83[k]
                  + f_0 * mf_83[k];

        t_84[k] = -5.0 * kf_84[k]
                  + f_0 * mf_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, kf_85, kf_86, kf_87, kf_88, kf_89, \
                         mf_85, mf_86, mf_87, mf_88, mf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -5.0 * kf_85[k]
                  + f_0 * mf_85[k];

        t_86[k] = -5.0 * kf_86[k]
                  + f_0 * mf_86[k];

        t_87[k] = -5.0 * kf_87[k]
                  + f_0 * mf_87[k];

        t_88[k] = -5.0 * kf_88[k]
                  + f_0 * mf_88[k];

        t_89[k] = -5.0 * kf_89[k]
                  + f_0 * mf_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, kf_90, kf_91, kf_92, kf_93, kf_94, \
                         mf_90, mf_91, mf_92, mf_93, mf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -5.0 * kf_90[k]
                  + f_0 * mf_90[k];

        t_91[k] = -5.0 * kf_91[k]
                  + f_0 * mf_91[k];

        t_92[k] = -5.0 * kf_92[k]
                  + f_0 * mf_92[k];

        t_93[k] = -5.0 * kf_93[k]
                  + f_0 * mf_93[k];

        t_94[k] = -5.0 * kf_94[k]
                  + f_0 * mf_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, kf_95, kf_96, kf_97, kf_98, kf_99, \
                         mf_95, mf_96, mf_97, mf_98, mf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -5.0 * kf_95[k]
                  + f_0 * mf_95[k];

        t_96[k] = -5.0 * kf_96[k]
                  + f_0 * mf_96[k];

        t_97[k] = -5.0 * kf_97[k]
                  + f_0 * mf_97[k];

        t_98[k] = -5.0 * kf_98[k]
                  + f_0 * mf_98[k];

        t_99[k] = -5.0 * kf_99[k]
                  + f_0 * mf_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, kf_100, kf_101, kf_102, kf_103, \
                         kf_104, mf_100, mf_101, mf_102, mf_103, \
                         mf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -4.0 * kf_100[k]
                   + f_0 * mf_100[k];

        t_101[k] = -4.0 * kf_101[k]
                   + f_0 * mf_101[k];

        t_102[k] = -4.0 * kf_102[k]
                   + f_0 * mf_102[k];

        t_103[k] = -4.0 * kf_103[k]
                   + f_0 * mf_103[k];

        t_104[k] = -4.0 * kf_104[k]
                   + f_0 * mf_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, kf_105, kf_106, kf_107, kf_108, \
                         kf_109, mf_105, mf_106, mf_107, mf_108, \
                         mf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -4.0 * kf_105[k]
                   + f_0 * mf_105[k];

        t_106[k] = -4.0 * kf_106[k]
                   + f_0 * mf_106[k];

        t_107[k] = -4.0 * kf_107[k]
                   + f_0 * mf_107[k];

        t_108[k] = -4.0 * kf_108[k]
                   + f_0 * mf_108[k];

        t_109[k] = -4.0 * kf_109[k]
                   + f_0 * mf_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, kf_110, kf_111, kf_112, kf_113, \
                         kf_114, mf_110, mf_111, mf_112, mf_113, \
                         mf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -4.0 * kf_110[k]
                   + f_0 * mf_110[k];

        t_111[k] = -4.0 * kf_111[k]
                   + f_0 * mf_111[k];

        t_112[k] = -4.0 * kf_112[k]
                   + f_0 * mf_112[k];

        t_113[k] = -4.0 * kf_113[k]
                   + f_0 * mf_113[k];

        t_114[k] = -4.0 * kf_114[k]
                   + f_0 * mf_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, kf_115, kf_116, kf_117, kf_118, \
                         kf_119, mf_115, mf_116, mf_117, mf_118, \
                         mf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -4.0 * kf_115[k]
                   + f_0 * mf_115[k];

        t_116[k] = -4.0 * kf_116[k]
                   + f_0 * mf_116[k];

        t_117[k] = -4.0 * kf_117[k]
                   + f_0 * mf_117[k];

        t_118[k] = -4.0 * kf_118[k]
                   + f_0 * mf_118[k];

        t_119[k] = -4.0 * kf_119[k]
                   + f_0 * mf_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, kf_120, kf_121, kf_122, kf_123, \
                         kf_124, mf_120, mf_121, mf_122, mf_123, \
                         mf_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -4.0 * kf_120[k]
                   + f_0 * mf_120[k];

        t_121[k] = -4.0 * kf_121[k]
                   + f_0 * mf_121[k];

        t_122[k] = -4.0 * kf_122[k]
                   + f_0 * mf_122[k];

        t_123[k] = -4.0 * kf_123[k]
                   + f_0 * mf_123[k];

        t_124[k] = -4.0 * kf_124[k]
                   + f_0 * mf_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, kf_125, kf_126, kf_127, kf_128, \
                         kf_129, mf_125, mf_126, mf_127, mf_128, \
                         mf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -4.0 * kf_125[k]
                   + f_0 * mf_125[k];

        t_126[k] = -4.0 * kf_126[k]
                   + f_0 * mf_126[k];

        t_127[k] = -4.0 * kf_127[k]
                   + f_0 * mf_127[k];

        t_128[k] = -4.0 * kf_128[k]
                   + f_0 * mf_128[k];

        t_129[k] = -4.0 * kf_129[k]
                   + f_0 * mf_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, kf_130, kf_131, kf_132, kf_133, \
                         kf_134, mf_130, mf_131, mf_132, mf_133, \
                         mf_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -4.0 * kf_130[k]
                   + f_0 * mf_130[k];

        t_131[k] = -4.0 * kf_131[k]
                   + f_0 * mf_131[k];

        t_132[k] = -4.0 * kf_132[k]
                   + f_0 * mf_132[k];

        t_133[k] = -4.0 * kf_133[k]
                   + f_0 * mf_133[k];

        t_134[k] = -4.0 * kf_134[k]
                   + f_0 * mf_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, kf_135, kf_136, kf_137, kf_138, \
                         kf_139, mf_135, mf_136, mf_137, mf_138, \
                         mf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -4.0 * kf_135[k]
                   + f_0 * mf_135[k];

        t_136[k] = -4.0 * kf_136[k]
                   + f_0 * mf_136[k];

        t_137[k] = -4.0 * kf_137[k]
                   + f_0 * mf_137[k];

        t_138[k] = -4.0 * kf_138[k]
                   + f_0 * mf_138[k];

        t_139[k] = -4.0 * kf_139[k]
                   + f_0 * mf_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, kf_140, kf_141, kf_142, kf_143, \
                         kf_144, mf_140, mf_141, mf_142, mf_143, \
                         mf_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -4.0 * kf_140[k]
                   + f_0 * mf_140[k];

        t_141[k] = -4.0 * kf_141[k]
                   + f_0 * mf_141[k];

        t_142[k] = -4.0 * kf_142[k]
                   + f_0 * mf_142[k];

        t_143[k] = -4.0 * kf_143[k]
                   + f_0 * mf_143[k];

        t_144[k] = -4.0 * kf_144[k]
                   + f_0 * mf_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, kf_145, kf_146, kf_147, kf_148, \
                         kf_149, mf_145, mf_146, mf_147, mf_148, \
                         mf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -4.0 * kf_145[k]
                   + f_0 * mf_145[k];

        t_146[k] = -4.0 * kf_146[k]
                   + f_0 * mf_146[k];

        t_147[k] = -4.0 * kf_147[k]
                   + f_0 * mf_147[k];

        t_148[k] = -4.0 * kf_148[k]
                   + f_0 * mf_148[k];

        t_149[k] = -4.0 * kf_149[k]
                   + f_0 * mf_149[k];
    }
}

static auto
compute_prim_geom_10_lf_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kf, const size_t mf,
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

    const auto *kf_150 = buffer.data(kf + 150);
    const auto *kf_151 = buffer.data(kf + 151);
    const auto *kf_152 = buffer.data(kf + 152);
    const auto *kf_153 = buffer.data(kf + 153);
    const auto *kf_154 = buffer.data(kf + 154);
    const auto *kf_155 = buffer.data(kf + 155);
    const auto *kf_156 = buffer.data(kf + 156);
    const auto *kf_157 = buffer.data(kf + 157);
    const auto *kf_158 = buffer.data(kf + 158);
    const auto *kf_159 = buffer.data(kf + 159);
    const auto *kf_160 = buffer.data(kf + 160);
    const auto *kf_161 = buffer.data(kf + 161);
    const auto *kf_162 = buffer.data(kf + 162);
    const auto *kf_163 = buffer.data(kf + 163);
    const auto *kf_164 = buffer.data(kf + 164);
    const auto *kf_165 = buffer.data(kf + 165);
    const auto *kf_166 = buffer.data(kf + 166);
    const auto *kf_167 = buffer.data(kf + 167);
    const auto *kf_168 = buffer.data(kf + 168);
    const auto *kf_169 = buffer.data(kf + 169);
    const auto *kf_170 = buffer.data(kf + 170);
    const auto *kf_171 = buffer.data(kf + 171);
    const auto *kf_172 = buffer.data(kf + 172);
    const auto *kf_173 = buffer.data(kf + 173);
    const auto *kf_174 = buffer.data(kf + 174);
    const auto *kf_175 = buffer.data(kf + 175);
    const auto *kf_176 = buffer.data(kf + 176);
    const auto *kf_177 = buffer.data(kf + 177);
    const auto *kf_178 = buffer.data(kf + 178);
    const auto *kf_179 = buffer.data(kf + 179);
    const auto *kf_180 = buffer.data(kf + 180);
    const auto *kf_181 = buffer.data(kf + 181);
    const auto *kf_182 = buffer.data(kf + 182);
    const auto *kf_183 = buffer.data(kf + 183);
    const auto *kf_184 = buffer.data(kf + 184);
    const auto *kf_185 = buffer.data(kf + 185);
    const auto *kf_186 = buffer.data(kf + 186);
    const auto *kf_187 = buffer.data(kf + 187);
    const auto *kf_188 = buffer.data(kf + 188);
    const auto *kf_189 = buffer.data(kf + 189);
    const auto *kf_190 = buffer.data(kf + 190);
    const auto *kf_191 = buffer.data(kf + 191);
    const auto *kf_192 = buffer.data(kf + 192);
    const auto *kf_193 = buffer.data(kf + 193);
    const auto *kf_194 = buffer.data(kf + 194);
    const auto *kf_195 = buffer.data(kf + 195);
    const auto *kf_196 = buffer.data(kf + 196);
    const auto *kf_197 = buffer.data(kf + 197);
    const auto *kf_198 = buffer.data(kf + 198);
    const auto *kf_199 = buffer.data(kf + 199);
    const auto *kf_200 = buffer.data(kf + 200);
    const auto *kf_201 = buffer.data(kf + 201);
    const auto *kf_202 = buffer.data(kf + 202);
    const auto *kf_203 = buffer.data(kf + 203);
    const auto *kf_204 = buffer.data(kf + 204);
    const auto *kf_205 = buffer.data(kf + 205);
    const auto *kf_206 = buffer.data(kf + 206);
    const auto *kf_207 = buffer.data(kf + 207);
    const auto *kf_208 = buffer.data(kf + 208);
    const auto *kf_209 = buffer.data(kf + 209);
    const auto *kf_210 = buffer.data(kf + 210);
    const auto *kf_211 = buffer.data(kf + 211);
    const auto *kf_212 = buffer.data(kf + 212);
    const auto *kf_213 = buffer.data(kf + 213);
    const auto *kf_214 = buffer.data(kf + 214);
    const auto *kf_215 = buffer.data(kf + 215);
    const auto *kf_216 = buffer.data(kf + 216);
    const auto *kf_217 = buffer.data(kf + 217);
    const auto *kf_218 = buffer.data(kf + 218);
    const auto *kf_219 = buffer.data(kf + 219);
    const auto *kf_220 = buffer.data(kf + 220);
    const auto *kf_221 = buffer.data(kf + 221);
    const auto *kf_222 = buffer.data(kf + 222);
    const auto *kf_223 = buffer.data(kf + 223);
    const auto *kf_224 = buffer.data(kf + 224);
    const auto *kf_225 = buffer.data(kf + 225);
    const auto *kf_226 = buffer.data(kf + 226);
    const auto *kf_227 = buffer.data(kf + 227);
    const auto *kf_228 = buffer.data(kf + 228);
    const auto *kf_229 = buffer.data(kf + 229);
    const auto *kf_230 = buffer.data(kf + 230);
    const auto *kf_231 = buffer.data(kf + 231);
    const auto *kf_232 = buffer.data(kf + 232);
    const auto *kf_233 = buffer.data(kf + 233);
    const auto *kf_234 = buffer.data(kf + 234);
    const auto *kf_235 = buffer.data(kf + 235);
    const auto *kf_236 = buffer.data(kf + 236);
    const auto *kf_237 = buffer.data(kf + 237);
    const auto *kf_238 = buffer.data(kf + 238);
    const auto *kf_239 = buffer.data(kf + 239);
    const auto *kf_240 = buffer.data(kf + 240);
    const auto *kf_241 = buffer.data(kf + 241);
    const auto *kf_242 = buffer.data(kf + 242);
    const auto *kf_243 = buffer.data(kf + 243);
    const auto *kf_244 = buffer.data(kf + 244);
    const auto *kf_245 = buffer.data(kf + 245);
    const auto *kf_246 = buffer.data(kf + 246);
    const auto *kf_247 = buffer.data(kf + 247);
    const auto *kf_248 = buffer.data(kf + 248);
    const auto *kf_249 = buffer.data(kf + 249);
    const auto *kf_250 = buffer.data(kf + 250);
    const auto *kf_251 = buffer.data(kf + 251);
    const auto *kf_252 = buffer.data(kf + 252);
    const auto *kf_253 = buffer.data(kf + 253);
    const auto *kf_254 = buffer.data(kf + 254);
    const auto *kf_255 = buffer.data(kf + 255);
    const auto *kf_256 = buffer.data(kf + 256);
    const auto *kf_257 = buffer.data(kf + 257);
    const auto *kf_258 = buffer.data(kf + 258);
    const auto *kf_259 = buffer.data(kf + 259);
    const auto *kf_260 = buffer.data(kf + 260);
    const auto *kf_261 = buffer.data(kf + 261);
    const auto *kf_262 = buffer.data(kf + 262);
    const auto *kf_263 = buffer.data(kf + 263);
    const auto *kf_264 = buffer.data(kf + 264);
    const auto *kf_265 = buffer.data(kf + 265);
    const auto *kf_266 = buffer.data(kf + 266);
    const auto *kf_267 = buffer.data(kf + 267);
    const auto *kf_268 = buffer.data(kf + 268);
    const auto *kf_269 = buffer.data(kf + 269);
    const auto *kf_270 = buffer.data(kf + 270);
    const auto *kf_271 = buffer.data(kf + 271);
    const auto *kf_272 = buffer.data(kf + 272);
    const auto *kf_273 = buffer.data(kf + 273);
    const auto *kf_274 = buffer.data(kf + 274);
    const auto *kf_275 = buffer.data(kf + 275);
    const auto *kf_276 = buffer.data(kf + 276);
    const auto *kf_277 = buffer.data(kf + 277);
    const auto *kf_278 = buffer.data(kf + 278);
    const auto *kf_279 = buffer.data(kf + 279);
    const auto *kf_280 = buffer.data(kf + 280);
    const auto *kf_281 = buffer.data(kf + 281);
    const auto *kf_282 = buffer.data(kf + 282);
    const auto *kf_283 = buffer.data(kf + 283);
    const auto *kf_284 = buffer.data(kf + 284);
    const auto *kf_285 = buffer.data(kf + 285);
    const auto *kf_286 = buffer.data(kf + 286);
    const auto *kf_287 = buffer.data(kf + 287);
    const auto *kf_288 = buffer.data(kf + 288);
    const auto *kf_289 = buffer.data(kf + 289);
    const auto *kf_290 = buffer.data(kf + 290);
    const auto *kf_291 = buffer.data(kf + 291);
    const auto *kf_292 = buffer.data(kf + 292);
    const auto *kf_293 = buffer.data(kf + 293);
    const auto *kf_294 = buffer.data(kf + 294);
    const auto *kf_295 = buffer.data(kf + 295);
    const auto *kf_296 = buffer.data(kf + 296);
    const auto *kf_297 = buffer.data(kf + 297);
    const auto *kf_298 = buffer.data(kf + 298);
    const auto *kf_299 = buffer.data(kf + 299);

    const auto *mf_150 = buffer.data(mf + 150);
    const auto *mf_151 = buffer.data(mf + 151);
    const auto *mf_152 = buffer.data(mf + 152);
    const auto *mf_153 = buffer.data(mf + 153);
    const auto *mf_154 = buffer.data(mf + 154);
    const auto *mf_155 = buffer.data(mf + 155);
    const auto *mf_156 = buffer.data(mf + 156);
    const auto *mf_157 = buffer.data(mf + 157);
    const auto *mf_158 = buffer.data(mf + 158);
    const auto *mf_159 = buffer.data(mf + 159);
    const auto *mf_160 = buffer.data(mf + 160);
    const auto *mf_161 = buffer.data(mf + 161);
    const auto *mf_162 = buffer.data(mf + 162);
    const auto *mf_163 = buffer.data(mf + 163);
    const auto *mf_164 = buffer.data(mf + 164);
    const auto *mf_165 = buffer.data(mf + 165);
    const auto *mf_166 = buffer.data(mf + 166);
    const auto *mf_167 = buffer.data(mf + 167);
    const auto *mf_168 = buffer.data(mf + 168);
    const auto *mf_169 = buffer.data(mf + 169);
    const auto *mf_170 = buffer.data(mf + 170);
    const auto *mf_171 = buffer.data(mf + 171);
    const auto *mf_172 = buffer.data(mf + 172);
    const auto *mf_173 = buffer.data(mf + 173);
    const auto *mf_174 = buffer.data(mf + 174);
    const auto *mf_175 = buffer.data(mf + 175);
    const auto *mf_176 = buffer.data(mf + 176);
    const auto *mf_177 = buffer.data(mf + 177);
    const auto *mf_178 = buffer.data(mf + 178);
    const auto *mf_179 = buffer.data(mf + 179);
    const auto *mf_180 = buffer.data(mf + 180);
    const auto *mf_181 = buffer.data(mf + 181);
    const auto *mf_182 = buffer.data(mf + 182);
    const auto *mf_183 = buffer.data(mf + 183);
    const auto *mf_184 = buffer.data(mf + 184);
    const auto *mf_185 = buffer.data(mf + 185);
    const auto *mf_186 = buffer.data(mf + 186);
    const auto *mf_187 = buffer.data(mf + 187);
    const auto *mf_188 = buffer.data(mf + 188);
    const auto *mf_189 = buffer.data(mf + 189);
    const auto *mf_190 = buffer.data(mf + 190);
    const auto *mf_191 = buffer.data(mf + 191);
    const auto *mf_192 = buffer.data(mf + 192);
    const auto *mf_193 = buffer.data(mf + 193);
    const auto *mf_194 = buffer.data(mf + 194);
    const auto *mf_195 = buffer.data(mf + 195);
    const auto *mf_196 = buffer.data(mf + 196);
    const auto *mf_197 = buffer.data(mf + 197);
    const auto *mf_198 = buffer.data(mf + 198);
    const auto *mf_199 = buffer.data(mf + 199);
    const auto *mf_200 = buffer.data(mf + 200);
    const auto *mf_201 = buffer.data(mf + 201);
    const auto *mf_202 = buffer.data(mf + 202);
    const auto *mf_203 = buffer.data(mf + 203);
    const auto *mf_204 = buffer.data(mf + 204);
    const auto *mf_205 = buffer.data(mf + 205);
    const auto *mf_206 = buffer.data(mf + 206);
    const auto *mf_207 = buffer.data(mf + 207);
    const auto *mf_208 = buffer.data(mf + 208);
    const auto *mf_209 = buffer.data(mf + 209);
    const auto *mf_210 = buffer.data(mf + 210);
    const auto *mf_211 = buffer.data(mf + 211);
    const auto *mf_212 = buffer.data(mf + 212);
    const auto *mf_213 = buffer.data(mf + 213);
    const auto *mf_214 = buffer.data(mf + 214);
    const auto *mf_215 = buffer.data(mf + 215);
    const auto *mf_216 = buffer.data(mf + 216);
    const auto *mf_217 = buffer.data(mf + 217);
    const auto *mf_218 = buffer.data(mf + 218);
    const auto *mf_219 = buffer.data(mf + 219);
    const auto *mf_220 = buffer.data(mf + 220);
    const auto *mf_221 = buffer.data(mf + 221);
    const auto *mf_222 = buffer.data(mf + 222);
    const auto *mf_223 = buffer.data(mf + 223);
    const auto *mf_224 = buffer.data(mf + 224);
    const auto *mf_225 = buffer.data(mf + 225);
    const auto *mf_226 = buffer.data(mf + 226);
    const auto *mf_227 = buffer.data(mf + 227);
    const auto *mf_228 = buffer.data(mf + 228);
    const auto *mf_229 = buffer.data(mf + 229);
    const auto *mf_230 = buffer.data(mf + 230);
    const auto *mf_231 = buffer.data(mf + 231);
    const auto *mf_232 = buffer.data(mf + 232);
    const auto *mf_233 = buffer.data(mf + 233);
    const auto *mf_234 = buffer.data(mf + 234);
    const auto *mf_235 = buffer.data(mf + 235);
    const auto *mf_236 = buffer.data(mf + 236);
    const auto *mf_237 = buffer.data(mf + 237);
    const auto *mf_238 = buffer.data(mf + 238);
    const auto *mf_239 = buffer.data(mf + 239);
    const auto *mf_240 = buffer.data(mf + 240);
    const auto *mf_241 = buffer.data(mf + 241);
    const auto *mf_242 = buffer.data(mf + 242);
    const auto *mf_243 = buffer.data(mf + 243);
    const auto *mf_244 = buffer.data(mf + 244);
    const auto *mf_245 = buffer.data(mf + 245);
    const auto *mf_246 = buffer.data(mf + 246);
    const auto *mf_247 = buffer.data(mf + 247);
    const auto *mf_248 = buffer.data(mf + 248);
    const auto *mf_249 = buffer.data(mf + 249);
    const auto *mf_250 = buffer.data(mf + 250);
    const auto *mf_251 = buffer.data(mf + 251);
    const auto *mf_252 = buffer.data(mf + 252);
    const auto *mf_253 = buffer.data(mf + 253);
    const auto *mf_254 = buffer.data(mf + 254);
    const auto *mf_255 = buffer.data(mf + 255);
    const auto *mf_256 = buffer.data(mf + 256);
    const auto *mf_257 = buffer.data(mf + 257);
    const auto *mf_258 = buffer.data(mf + 258);
    const auto *mf_259 = buffer.data(mf + 259);
    const auto *mf_260 = buffer.data(mf + 260);
    const auto *mf_261 = buffer.data(mf + 261);
    const auto *mf_262 = buffer.data(mf + 262);
    const auto *mf_263 = buffer.data(mf + 263);
    const auto *mf_264 = buffer.data(mf + 264);
    const auto *mf_265 = buffer.data(mf + 265);
    const auto *mf_266 = buffer.data(mf + 266);
    const auto *mf_267 = buffer.data(mf + 267);
    const auto *mf_268 = buffer.data(mf + 268);
    const auto *mf_269 = buffer.data(mf + 269);
    const auto *mf_270 = buffer.data(mf + 270);
    const auto *mf_271 = buffer.data(mf + 271);
    const auto *mf_272 = buffer.data(mf + 272);
    const auto *mf_273 = buffer.data(mf + 273);
    const auto *mf_274 = buffer.data(mf + 274);
    const auto *mf_275 = buffer.data(mf + 275);
    const auto *mf_276 = buffer.data(mf + 276);
    const auto *mf_277 = buffer.data(mf + 277);
    const auto *mf_278 = buffer.data(mf + 278);
    const auto *mf_279 = buffer.data(mf + 279);
    const auto *mf_280 = buffer.data(mf + 280);
    const auto *mf_281 = buffer.data(mf + 281);
    const auto *mf_282 = buffer.data(mf + 282);
    const auto *mf_283 = buffer.data(mf + 283);
    const auto *mf_284 = buffer.data(mf + 284);
    const auto *mf_285 = buffer.data(mf + 285);
    const auto *mf_286 = buffer.data(mf + 286);
    const auto *mf_287 = buffer.data(mf + 287);
    const auto *mf_288 = buffer.data(mf + 288);
    const auto *mf_289 = buffer.data(mf + 289);
    const auto *mf_290 = buffer.data(mf + 290);
    const auto *mf_291 = buffer.data(mf + 291);
    const auto *mf_292 = buffer.data(mf + 292);
    const auto *mf_293 = buffer.data(mf + 293);
    const auto *mf_294 = buffer.data(mf + 294);
    const auto *mf_295 = buffer.data(mf + 295);
    const auto *mf_296 = buffer.data(mf + 296);
    const auto *mf_297 = buffer.data(mf + 297);
    const auto *mf_298 = buffer.data(mf + 298);
    const auto *mf_299 = buffer.data(mf + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, kf_150, kf_151, kf_152, kf_153, \
                         kf_154, mf_150, mf_151, mf_152, mf_153, \
                         mf_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -3.0 * kf_150[k]
                   + f_0 * mf_150[k];

        t_151[k] = -3.0 * kf_151[k]
                   + f_0 * mf_151[k];

        t_152[k] = -3.0 * kf_152[k]
                   + f_0 * mf_152[k];

        t_153[k] = -3.0 * kf_153[k]
                   + f_0 * mf_153[k];

        t_154[k] = -3.0 * kf_154[k]
                   + f_0 * mf_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, kf_155, kf_156, kf_157, kf_158, \
                         kf_159, mf_155, mf_156, mf_157, mf_158, \
                         mf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -3.0 * kf_155[k]
                   + f_0 * mf_155[k];

        t_156[k] = -3.0 * kf_156[k]
                   + f_0 * mf_156[k];

        t_157[k] = -3.0 * kf_157[k]
                   + f_0 * mf_157[k];

        t_158[k] = -3.0 * kf_158[k]
                   + f_0 * mf_158[k];

        t_159[k] = -3.0 * kf_159[k]
                   + f_0 * mf_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, kf_160, kf_161, kf_162, kf_163, \
                         kf_164, mf_160, mf_161, mf_162, mf_163, \
                         mf_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -3.0 * kf_160[k]
                   + f_0 * mf_160[k];

        t_161[k] = -3.0 * kf_161[k]
                   + f_0 * mf_161[k];

        t_162[k] = -3.0 * kf_162[k]
                   + f_0 * mf_162[k];

        t_163[k] = -3.0 * kf_163[k]
                   + f_0 * mf_163[k];

        t_164[k] = -3.0 * kf_164[k]
                   + f_0 * mf_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, kf_165, kf_166, kf_167, kf_168, \
                         kf_169, mf_165, mf_166, mf_167, mf_168, \
                         mf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -3.0 * kf_165[k]
                   + f_0 * mf_165[k];

        t_166[k] = -3.0 * kf_166[k]
                   + f_0 * mf_166[k];

        t_167[k] = -3.0 * kf_167[k]
                   + f_0 * mf_167[k];

        t_168[k] = -3.0 * kf_168[k]
                   + f_0 * mf_168[k];

        t_169[k] = -3.0 * kf_169[k]
                   + f_0 * mf_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, kf_170, kf_171, kf_172, kf_173, \
                         kf_174, mf_170, mf_171, mf_172, mf_173, \
                         mf_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * kf_170[k]
                   + f_0 * mf_170[k];

        t_171[k] = -3.0 * kf_171[k]
                   + f_0 * mf_171[k];

        t_172[k] = -3.0 * kf_172[k]
                   + f_0 * mf_172[k];

        t_173[k] = -3.0 * kf_173[k]
                   + f_0 * mf_173[k];

        t_174[k] = -3.0 * kf_174[k]
                   + f_0 * mf_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, kf_175, kf_176, kf_177, kf_178, \
                         kf_179, mf_175, mf_176, mf_177, mf_178, \
                         mf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * kf_175[k]
                   + f_0 * mf_175[k];

        t_176[k] = -3.0 * kf_176[k]
                   + f_0 * mf_176[k];

        t_177[k] = -3.0 * kf_177[k]
                   + f_0 * mf_177[k];

        t_178[k] = -3.0 * kf_178[k]
                   + f_0 * mf_178[k];

        t_179[k] = -3.0 * kf_179[k]
                   + f_0 * mf_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, kf_180, kf_181, kf_182, kf_183, \
                         kf_184, mf_180, mf_181, mf_182, mf_183, \
                         mf_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * kf_180[k]
                   + f_0 * mf_180[k];

        t_181[k] = -3.0 * kf_181[k]
                   + f_0 * mf_181[k];

        t_182[k] = -3.0 * kf_182[k]
                   + f_0 * mf_182[k];

        t_183[k] = -3.0 * kf_183[k]
                   + f_0 * mf_183[k];

        t_184[k] = -3.0 * kf_184[k]
                   + f_0 * mf_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, kf_185, kf_186, kf_187, kf_188, \
                         kf_189, mf_185, mf_186, mf_187, mf_188, \
                         mf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * kf_185[k]
                   + f_0 * mf_185[k];

        t_186[k] = -3.0 * kf_186[k]
                   + f_0 * mf_186[k];

        t_187[k] = -3.0 * kf_187[k]
                   + f_0 * mf_187[k];

        t_188[k] = -3.0 * kf_188[k]
                   + f_0 * mf_188[k];

        t_189[k] = -3.0 * kf_189[k]
                   + f_0 * mf_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, kf_190, kf_191, kf_192, kf_193, \
                         kf_194, mf_190, mf_191, mf_192, mf_193, \
                         mf_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * kf_190[k]
                   + f_0 * mf_190[k];

        t_191[k] = -3.0 * kf_191[k]
                   + f_0 * mf_191[k];

        t_192[k] = -3.0 * kf_192[k]
                   + f_0 * mf_192[k];

        t_193[k] = -3.0 * kf_193[k]
                   + f_0 * mf_193[k];

        t_194[k] = -3.0 * kf_194[k]
                   + f_0 * mf_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, kf_195, kf_196, kf_197, kf_198, \
                         kf_199, mf_195, mf_196, mf_197, mf_198, \
                         mf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * kf_195[k]
                   + f_0 * mf_195[k];

        t_196[k] = -3.0 * kf_196[k]
                   + f_0 * mf_196[k];

        t_197[k] = -3.0 * kf_197[k]
                   + f_0 * mf_197[k];

        t_198[k] = -3.0 * kf_198[k]
                   + f_0 * mf_198[k];

        t_199[k] = -3.0 * kf_199[k]
                   + f_0 * mf_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, kf_200, kf_201, kf_202, kf_203, \
                         kf_204, mf_200, mf_201, mf_202, mf_203, \
                         mf_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -3.0 * kf_200[k]
                   + f_0 * mf_200[k];

        t_201[k] = -3.0 * kf_201[k]
                   + f_0 * mf_201[k];

        t_202[k] = -3.0 * kf_202[k]
                   + f_0 * mf_202[k];

        t_203[k] = -3.0 * kf_203[k]
                   + f_0 * mf_203[k];

        t_204[k] = -3.0 * kf_204[k]
                   + f_0 * mf_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, kf_205, kf_206, kf_207, kf_208, \
                         kf_209, mf_205, mf_206, mf_207, mf_208, \
                         mf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -3.0 * kf_205[k]
                   + f_0 * mf_205[k];

        t_206[k] = -3.0 * kf_206[k]
                   + f_0 * mf_206[k];

        t_207[k] = -3.0 * kf_207[k]
                   + f_0 * mf_207[k];

        t_208[k] = -3.0 * kf_208[k]
                   + f_0 * mf_208[k];

        t_209[k] = -3.0 * kf_209[k]
                   + f_0 * mf_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, kf_210, kf_211, kf_212, kf_213, \
                         kf_214, mf_210, mf_211, mf_212, mf_213, \
                         mf_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * kf_210[k]
                   + f_0 * mf_210[k];

        t_211[k] = -2.0 * kf_211[k]
                   + f_0 * mf_211[k];

        t_212[k] = -2.0 * kf_212[k]
                   + f_0 * mf_212[k];

        t_213[k] = -2.0 * kf_213[k]
                   + f_0 * mf_213[k];

        t_214[k] = -2.0 * kf_214[k]
                   + f_0 * mf_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, kf_215, kf_216, kf_217, kf_218, \
                         kf_219, mf_215, mf_216, mf_217, mf_218, \
                         mf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * kf_215[k]
                   + f_0 * mf_215[k];

        t_216[k] = -2.0 * kf_216[k]
                   + f_0 * mf_216[k];

        t_217[k] = -2.0 * kf_217[k]
                   + f_0 * mf_217[k];

        t_218[k] = -2.0 * kf_218[k]
                   + f_0 * mf_218[k];

        t_219[k] = -2.0 * kf_219[k]
                   + f_0 * mf_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, kf_220, kf_221, kf_222, kf_223, \
                         kf_224, mf_220, mf_221, mf_222, mf_223, \
                         mf_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * kf_220[k]
                   + f_0 * mf_220[k];

        t_221[k] = -2.0 * kf_221[k]
                   + f_0 * mf_221[k];

        t_222[k] = -2.0 * kf_222[k]
                   + f_0 * mf_222[k];

        t_223[k] = -2.0 * kf_223[k]
                   + f_0 * mf_223[k];

        t_224[k] = -2.0 * kf_224[k]
                   + f_0 * mf_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, kf_225, kf_226, kf_227, kf_228, \
                         kf_229, mf_225, mf_226, mf_227, mf_228, \
                         mf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -2.0 * kf_225[k]
                   + f_0 * mf_225[k];

        t_226[k] = -2.0 * kf_226[k]
                   + f_0 * mf_226[k];

        t_227[k] = -2.0 * kf_227[k]
                   + f_0 * mf_227[k];

        t_228[k] = -2.0 * kf_228[k]
                   + f_0 * mf_228[k];

        t_229[k] = -2.0 * kf_229[k]
                   + f_0 * mf_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, kf_230, kf_231, kf_232, kf_233, \
                         kf_234, mf_230, mf_231, mf_232, mf_233, \
                         mf_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -2.0 * kf_230[k]
                   + f_0 * mf_230[k];

        t_231[k] = -2.0 * kf_231[k]
                   + f_0 * mf_231[k];

        t_232[k] = -2.0 * kf_232[k]
                   + f_0 * mf_232[k];

        t_233[k] = -2.0 * kf_233[k]
                   + f_0 * mf_233[k];

        t_234[k] = -2.0 * kf_234[k]
                   + f_0 * mf_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, kf_235, kf_236, kf_237, kf_238, \
                         kf_239, mf_235, mf_236, mf_237, mf_238, \
                         mf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -2.0 * kf_235[k]
                   + f_0 * mf_235[k];

        t_236[k] = -2.0 * kf_236[k]
                   + f_0 * mf_236[k];

        t_237[k] = -2.0 * kf_237[k]
                   + f_0 * mf_237[k];

        t_238[k] = -2.0 * kf_238[k]
                   + f_0 * mf_238[k];

        t_239[k] = -2.0 * kf_239[k]
                   + f_0 * mf_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, kf_240, kf_241, kf_242, kf_243, \
                         kf_244, mf_240, mf_241, mf_242, mf_243, \
                         mf_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -2.0 * kf_240[k]
                   + f_0 * mf_240[k];

        t_241[k] = -2.0 * kf_241[k]
                   + f_0 * mf_241[k];

        t_242[k] = -2.0 * kf_242[k]
                   + f_0 * mf_242[k];

        t_243[k] = -2.0 * kf_243[k]
                   + f_0 * mf_243[k];

        t_244[k] = -2.0 * kf_244[k]
                   + f_0 * mf_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, kf_245, kf_246, kf_247, kf_248, \
                         kf_249, mf_245, mf_246, mf_247, mf_248, \
                         mf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -2.0 * kf_245[k]
                   + f_0 * mf_245[k];

        t_246[k] = -2.0 * kf_246[k]
                   + f_0 * mf_246[k];

        t_247[k] = -2.0 * kf_247[k]
                   + f_0 * mf_247[k];

        t_248[k] = -2.0 * kf_248[k]
                   + f_0 * mf_248[k];

        t_249[k] = -2.0 * kf_249[k]
                   + f_0 * mf_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, kf_250, kf_251, kf_252, kf_253, \
                         kf_254, mf_250, mf_251, mf_252, mf_253, \
                         mf_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * kf_250[k]
                   + f_0 * mf_250[k];

        t_251[k] = -2.0 * kf_251[k]
                   + f_0 * mf_251[k];

        t_252[k] = -2.0 * kf_252[k]
                   + f_0 * mf_252[k];

        t_253[k] = -2.0 * kf_253[k]
                   + f_0 * mf_253[k];

        t_254[k] = -2.0 * kf_254[k]
                   + f_0 * mf_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, kf_255, kf_256, kf_257, kf_258, \
                         kf_259, mf_255, mf_256, mf_257, mf_258, \
                         mf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * kf_255[k]
                   + f_0 * mf_255[k];

        t_256[k] = -2.0 * kf_256[k]
                   + f_0 * mf_256[k];

        t_257[k] = -2.0 * kf_257[k]
                   + f_0 * mf_257[k];

        t_258[k] = -2.0 * kf_258[k]
                   + f_0 * mf_258[k];

        t_259[k] = -2.0 * kf_259[k]
                   + f_0 * mf_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, kf_260, kf_261, kf_262, kf_263, \
                         kf_264, mf_260, mf_261, mf_262, mf_263, \
                         mf_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * kf_260[k]
                   + f_0 * mf_260[k];

        t_261[k] = -2.0 * kf_261[k]
                   + f_0 * mf_261[k];

        t_262[k] = -2.0 * kf_262[k]
                   + f_0 * mf_262[k];

        t_263[k] = -2.0 * kf_263[k]
                   + f_0 * mf_263[k];

        t_264[k] = -2.0 * kf_264[k]
                   + f_0 * mf_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, kf_265, kf_266, kf_267, kf_268, \
                         kf_269, mf_265, mf_266, mf_267, mf_268, \
                         mf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * kf_265[k]
                   + f_0 * mf_265[k];

        t_266[k] = -2.0 * kf_266[k]
                   + f_0 * mf_266[k];

        t_267[k] = -2.0 * kf_267[k]
                   + f_0 * mf_267[k];

        t_268[k] = -2.0 * kf_268[k]
                   + f_0 * mf_268[k];

        t_269[k] = -2.0 * kf_269[k]
                   + f_0 * mf_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, kf_270, kf_271, kf_272, kf_273, \
                         kf_274, mf_270, mf_271, mf_272, mf_273, \
                         mf_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * kf_270[k]
                   + f_0 * mf_270[k];

        t_271[k] = -2.0 * kf_271[k]
                   + f_0 * mf_271[k];

        t_272[k] = -2.0 * kf_272[k]
                   + f_0 * mf_272[k];

        t_273[k] = -2.0 * kf_273[k]
                   + f_0 * mf_273[k];

        t_274[k] = -2.0 * kf_274[k]
                   + f_0 * mf_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, kf_275, kf_276, kf_277, kf_278, \
                         kf_279, mf_275, mf_276, mf_277, mf_278, \
                         mf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * kf_275[k]
                   + f_0 * mf_275[k];

        t_276[k] = -2.0 * kf_276[k]
                   + f_0 * mf_276[k];

        t_277[k] = -2.0 * kf_277[k]
                   + f_0 * mf_277[k];

        t_278[k] = -2.0 * kf_278[k]
                   + f_0 * mf_278[k];

        t_279[k] = -2.0 * kf_279[k]
                   + f_0 * mf_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, kf_280, kf_281, kf_282, kf_283, \
                         kf_284, mf_280, mf_281, mf_282, mf_283, \
                         mf_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -kf_280[k]
                   + f_0 * mf_280[k];

        t_281[k] = -kf_281[k]
                   + f_0 * mf_281[k];

        t_282[k] = -kf_282[k]
                   + f_0 * mf_282[k];

        t_283[k] = -kf_283[k]
                   + f_0 * mf_283[k];

        t_284[k] = -kf_284[k]
                   + f_0 * mf_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, kf_285, kf_286, kf_287, kf_288, \
                         kf_289, mf_285, mf_286, mf_287, mf_288, \
                         mf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -kf_285[k]
                   + f_0 * mf_285[k];

        t_286[k] = -kf_286[k]
                   + f_0 * mf_286[k];

        t_287[k] = -kf_287[k]
                   + f_0 * mf_287[k];

        t_288[k] = -kf_288[k]
                   + f_0 * mf_288[k];

        t_289[k] = -kf_289[k]
                   + f_0 * mf_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, kf_290, kf_291, kf_292, kf_293, \
                         kf_294, mf_290, mf_291, mf_292, mf_293, \
                         mf_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -kf_290[k]
                   + f_0 * mf_290[k];

        t_291[k] = -kf_291[k]
                   + f_0 * mf_291[k];

        t_292[k] = -kf_292[k]
                   + f_0 * mf_292[k];

        t_293[k] = -kf_293[k]
                   + f_0 * mf_293[k];

        t_294[k] = -kf_294[k]
                   + f_0 * mf_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, kf_295, kf_296, kf_297, kf_298, \
                         kf_299, mf_295, mf_296, mf_297, mf_298, \
                         mf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -kf_295[k]
                   + f_0 * mf_295[k];

        t_296[k] = -kf_296[k]
                   + f_0 * mf_296[k];

        t_297[k] = -kf_297[k]
                   + f_0 * mf_297[k];

        t_298[k] = -kf_298[k]
                   + f_0 * mf_298[k];

        t_299[k] = -kf_299[k]
                   + f_0 * mf_299[k];
    }
}

static auto
compute_prim_geom_10_lf_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kf, const size_t mf,
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

    const auto *kf_300 = buffer.data(kf + 300);
    const auto *kf_301 = buffer.data(kf + 301);
    const auto *kf_302 = buffer.data(kf + 302);
    const auto *kf_303 = buffer.data(kf + 303);
    const auto *kf_304 = buffer.data(kf + 304);
    const auto *kf_305 = buffer.data(kf + 305);
    const auto *kf_306 = buffer.data(kf + 306);
    const auto *kf_307 = buffer.data(kf + 307);
    const auto *kf_308 = buffer.data(kf + 308);
    const auto *kf_309 = buffer.data(kf + 309);
    const auto *kf_310 = buffer.data(kf + 310);
    const auto *kf_311 = buffer.data(kf + 311);
    const auto *kf_312 = buffer.data(kf + 312);
    const auto *kf_313 = buffer.data(kf + 313);
    const auto *kf_314 = buffer.data(kf + 314);
    const auto *kf_315 = buffer.data(kf + 315);
    const auto *kf_316 = buffer.data(kf + 316);
    const auto *kf_317 = buffer.data(kf + 317);
    const auto *kf_318 = buffer.data(kf + 318);
    const auto *kf_319 = buffer.data(kf + 319);
    const auto *kf_320 = buffer.data(kf + 320);
    const auto *kf_321 = buffer.data(kf + 321);
    const auto *kf_322 = buffer.data(kf + 322);
    const auto *kf_323 = buffer.data(kf + 323);
    const auto *kf_324 = buffer.data(kf + 324);
    const auto *kf_325 = buffer.data(kf + 325);
    const auto *kf_326 = buffer.data(kf + 326);
    const auto *kf_327 = buffer.data(kf + 327);
    const auto *kf_328 = buffer.data(kf + 328);
    const auto *kf_329 = buffer.data(kf + 329);
    const auto *kf_330 = buffer.data(kf + 330);
    const auto *kf_331 = buffer.data(kf + 331);
    const auto *kf_332 = buffer.data(kf + 332);
    const auto *kf_333 = buffer.data(kf + 333);
    const auto *kf_334 = buffer.data(kf + 334);
    const auto *kf_335 = buffer.data(kf + 335);
    const auto *kf_336 = buffer.data(kf + 336);
    const auto *kf_337 = buffer.data(kf + 337);
    const auto *kf_338 = buffer.data(kf + 338);
    const auto *kf_339 = buffer.data(kf + 339);
    const auto *kf_340 = buffer.data(kf + 340);
    const auto *kf_341 = buffer.data(kf + 341);
    const auto *kf_342 = buffer.data(kf + 342);
    const auto *kf_343 = buffer.data(kf + 343);
    const auto *kf_344 = buffer.data(kf + 344);
    const auto *kf_345 = buffer.data(kf + 345);
    const auto *kf_346 = buffer.data(kf + 346);
    const auto *kf_347 = buffer.data(kf + 347);
    const auto *kf_348 = buffer.data(kf + 348);
    const auto *kf_349 = buffer.data(kf + 349);
    const auto *kf_350 = buffer.data(kf + 350);
    const auto *kf_351 = buffer.data(kf + 351);
    const auto *kf_352 = buffer.data(kf + 352);
    const auto *kf_353 = buffer.data(kf + 353);
    const auto *kf_354 = buffer.data(kf + 354);
    const auto *kf_355 = buffer.data(kf + 355);
    const auto *kf_356 = buffer.data(kf + 356);
    const auto *kf_357 = buffer.data(kf + 357);
    const auto *kf_358 = buffer.data(kf + 358);
    const auto *kf_359 = buffer.data(kf + 359);

    const auto *mf_300 = buffer.data(mf + 300);
    const auto *mf_301 = buffer.data(mf + 301);
    const auto *mf_302 = buffer.data(mf + 302);
    const auto *mf_303 = buffer.data(mf + 303);
    const auto *mf_304 = buffer.data(mf + 304);
    const auto *mf_305 = buffer.data(mf + 305);
    const auto *mf_306 = buffer.data(mf + 306);
    const auto *mf_307 = buffer.data(mf + 307);
    const auto *mf_308 = buffer.data(mf + 308);
    const auto *mf_309 = buffer.data(mf + 309);
    const auto *mf_310 = buffer.data(mf + 310);
    const auto *mf_311 = buffer.data(mf + 311);
    const auto *mf_312 = buffer.data(mf + 312);
    const auto *mf_313 = buffer.data(mf + 313);
    const auto *mf_314 = buffer.data(mf + 314);
    const auto *mf_315 = buffer.data(mf + 315);
    const auto *mf_316 = buffer.data(mf + 316);
    const auto *mf_317 = buffer.data(mf + 317);
    const auto *mf_318 = buffer.data(mf + 318);
    const auto *mf_319 = buffer.data(mf + 319);
    const auto *mf_320 = buffer.data(mf + 320);
    const auto *mf_321 = buffer.data(mf + 321);
    const auto *mf_322 = buffer.data(mf + 322);
    const auto *mf_323 = buffer.data(mf + 323);
    const auto *mf_324 = buffer.data(mf + 324);
    const auto *mf_325 = buffer.data(mf + 325);
    const auto *mf_326 = buffer.data(mf + 326);
    const auto *mf_327 = buffer.data(mf + 327);
    const auto *mf_328 = buffer.data(mf + 328);
    const auto *mf_329 = buffer.data(mf + 329);
    const auto *mf_330 = buffer.data(mf + 330);
    const auto *mf_331 = buffer.data(mf + 331);
    const auto *mf_332 = buffer.data(mf + 332);
    const auto *mf_333 = buffer.data(mf + 333);
    const auto *mf_334 = buffer.data(mf + 334);
    const auto *mf_335 = buffer.data(mf + 335);
    const auto *mf_336 = buffer.data(mf + 336);
    const auto *mf_337 = buffer.data(mf + 337);
    const auto *mf_338 = buffer.data(mf + 338);
    const auto *mf_339 = buffer.data(mf + 339);
    const auto *mf_340 = buffer.data(mf + 340);
    const auto *mf_341 = buffer.data(mf + 341);
    const auto *mf_342 = buffer.data(mf + 342);
    const auto *mf_343 = buffer.data(mf + 343);
    const auto *mf_344 = buffer.data(mf + 344);
    const auto *mf_345 = buffer.data(mf + 345);
    const auto *mf_346 = buffer.data(mf + 346);
    const auto *mf_347 = buffer.data(mf + 347);
    const auto *mf_348 = buffer.data(mf + 348);
    const auto *mf_349 = buffer.data(mf + 349);
    const auto *mf_350 = buffer.data(mf + 350);
    const auto *mf_351 = buffer.data(mf + 351);
    const auto *mf_352 = buffer.data(mf + 352);
    const auto *mf_353 = buffer.data(mf + 353);
    const auto *mf_354 = buffer.data(mf + 354);
    const auto *mf_355 = buffer.data(mf + 355);
    const auto *mf_356 = buffer.data(mf + 356);
    const auto *mf_357 = buffer.data(mf + 357);
    const auto *mf_358 = buffer.data(mf + 358);
    const auto *mf_359 = buffer.data(mf + 359);
    const auto *mf_360 = buffer.data(mf + 360);
    const auto *mf_361 = buffer.data(mf + 361);
    const auto *mf_362 = buffer.data(mf + 362);
    const auto *mf_363 = buffer.data(mf + 363);
    const auto *mf_364 = buffer.data(mf + 364);
    const auto *mf_365 = buffer.data(mf + 365);
    const auto *mf_366 = buffer.data(mf + 366);
    const auto *mf_367 = buffer.data(mf + 367);
    const auto *mf_368 = buffer.data(mf + 368);
    const auto *mf_369 = buffer.data(mf + 369);
    const auto *mf_370 = buffer.data(mf + 370);
    const auto *mf_371 = buffer.data(mf + 371);
    const auto *mf_372 = buffer.data(mf + 372);
    const auto *mf_373 = buffer.data(mf + 373);
    const auto *mf_374 = buffer.data(mf + 374);
    const auto *mf_375 = buffer.data(mf + 375);
    const auto *mf_376 = buffer.data(mf + 376);
    const auto *mf_377 = buffer.data(mf + 377);
    const auto *mf_378 = buffer.data(mf + 378);
    const auto *mf_379 = buffer.data(mf + 379);
    const auto *mf_380 = buffer.data(mf + 380);
    const auto *mf_381 = buffer.data(mf + 381);
    const auto *mf_382 = buffer.data(mf + 382);
    const auto *mf_383 = buffer.data(mf + 383);
    const auto *mf_384 = buffer.data(mf + 384);
    const auto *mf_385 = buffer.data(mf + 385);
    const auto *mf_386 = buffer.data(mf + 386);
    const auto *mf_387 = buffer.data(mf + 387);
    const auto *mf_388 = buffer.data(mf + 388);
    const auto *mf_389 = buffer.data(mf + 389);
    const auto *mf_390 = buffer.data(mf + 390);
    const auto *mf_391 = buffer.data(mf + 391);
    const auto *mf_392 = buffer.data(mf + 392);
    const auto *mf_393 = buffer.data(mf + 393);
    const auto *mf_394 = buffer.data(mf + 394);
    const auto *mf_395 = buffer.data(mf + 395);
    const auto *mf_396 = buffer.data(mf + 396);
    const auto *mf_397 = buffer.data(mf + 397);
    const auto *mf_398 = buffer.data(mf + 398);
    const auto *mf_399 = buffer.data(mf + 399);
    const auto *mf_400 = buffer.data(mf + 400);
    const auto *mf_401 = buffer.data(mf + 401);
    const auto *mf_402 = buffer.data(mf + 402);
    const auto *mf_403 = buffer.data(mf + 403);
    const auto *mf_404 = buffer.data(mf + 404);
    const auto *mf_405 = buffer.data(mf + 405);
    const auto *mf_406 = buffer.data(mf + 406);
    const auto *mf_407 = buffer.data(mf + 407);
    const auto *mf_408 = buffer.data(mf + 408);
    const auto *mf_409 = buffer.data(mf + 409);
    const auto *mf_410 = buffer.data(mf + 410);
    const auto *mf_411 = buffer.data(mf + 411);
    const auto *mf_412 = buffer.data(mf + 412);
    const auto *mf_413 = buffer.data(mf + 413);
    const auto *mf_414 = buffer.data(mf + 414);
    const auto *mf_415 = buffer.data(mf + 415);
    const auto *mf_416 = buffer.data(mf + 416);
    const auto *mf_417 = buffer.data(mf + 417);
    const auto *mf_418 = buffer.data(mf + 418);
    const auto *mf_419 = buffer.data(mf + 419);
    const auto *mf_420 = buffer.data(mf + 420);
    const auto *mf_421 = buffer.data(mf + 421);
    const auto *mf_422 = buffer.data(mf + 422);
    const auto *mf_423 = buffer.data(mf + 423);
    const auto *mf_424 = buffer.data(mf + 424);
    const auto *mf_425 = buffer.data(mf + 425);
    const auto *mf_426 = buffer.data(mf + 426);
    const auto *mf_427 = buffer.data(mf + 427);
    const auto *mf_428 = buffer.data(mf + 428);
    const auto *mf_429 = buffer.data(mf + 429);
    const auto *mf_430 = buffer.data(mf + 430);
    const auto *mf_431 = buffer.data(mf + 431);
    const auto *mf_432 = buffer.data(mf + 432);
    const auto *mf_433 = buffer.data(mf + 433);
    const auto *mf_434 = buffer.data(mf + 434);
    const auto *mf_435 = buffer.data(mf + 435);
    const auto *mf_436 = buffer.data(mf + 436);
    const auto *mf_437 = buffer.data(mf + 437);
    const auto *mf_438 = buffer.data(mf + 438);
    const auto *mf_439 = buffer.data(mf + 439);
    const auto *mf_440 = buffer.data(mf + 440);
    const auto *mf_441 = buffer.data(mf + 441);
    const auto *mf_442 = buffer.data(mf + 442);
    const auto *mf_443 = buffer.data(mf + 443);
    const auto *mf_444 = buffer.data(mf + 444);
    const auto *mf_445 = buffer.data(mf + 445);
    const auto *mf_446 = buffer.data(mf + 446);
    const auto *mf_447 = buffer.data(mf + 447);
    const auto *mf_448 = buffer.data(mf + 448);
    const auto *mf_449 = buffer.data(mf + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, kf_300, kf_301, kf_302, kf_303, \
                         kf_304, mf_300, mf_301, mf_302, mf_303, \
                         mf_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -kf_300[k]
                   + f_0 * mf_300[k];

        t_301[k] = -kf_301[k]
                   + f_0 * mf_301[k];

        t_302[k] = -kf_302[k]
                   + f_0 * mf_302[k];

        t_303[k] = -kf_303[k]
                   + f_0 * mf_303[k];

        t_304[k] = -kf_304[k]
                   + f_0 * mf_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, kf_305, kf_306, kf_307, kf_308, \
                         kf_309, mf_305, mf_306, mf_307, mf_308, \
                         mf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -kf_305[k]
                   + f_0 * mf_305[k];

        t_306[k] = -kf_306[k]
                   + f_0 * mf_306[k];

        t_307[k] = -kf_307[k]
                   + f_0 * mf_307[k];

        t_308[k] = -kf_308[k]
                   + f_0 * mf_308[k];

        t_309[k] = -kf_309[k]
                   + f_0 * mf_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, kf_310, kf_311, kf_312, kf_313, \
                         kf_314, mf_310, mf_311, mf_312, mf_313, \
                         mf_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -kf_310[k]
                   + f_0 * mf_310[k];

        t_311[k] = -kf_311[k]
                   + f_0 * mf_311[k];

        t_312[k] = -kf_312[k]
                   + f_0 * mf_312[k];

        t_313[k] = -kf_313[k]
                   + f_0 * mf_313[k];

        t_314[k] = -kf_314[k]
                   + f_0 * mf_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, kf_315, kf_316, kf_317, kf_318, \
                         kf_319, mf_315, mf_316, mf_317, mf_318, \
                         mf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -kf_315[k]
                   + f_0 * mf_315[k];

        t_316[k] = -kf_316[k]
                   + f_0 * mf_316[k];

        t_317[k] = -kf_317[k]
                   + f_0 * mf_317[k];

        t_318[k] = -kf_318[k]
                   + f_0 * mf_318[k];

        t_319[k] = -kf_319[k]
                   + f_0 * mf_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, kf_320, kf_321, kf_322, kf_323, \
                         kf_324, mf_320, mf_321, mf_322, mf_323, \
                         mf_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -kf_320[k]
                   + f_0 * mf_320[k];

        t_321[k] = -kf_321[k]
                   + f_0 * mf_321[k];

        t_322[k] = -kf_322[k]
                   + f_0 * mf_322[k];

        t_323[k] = -kf_323[k]
                   + f_0 * mf_323[k];

        t_324[k] = -kf_324[k]
                   + f_0 * mf_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, kf_325, kf_326, kf_327, kf_328, \
                         kf_329, mf_325, mf_326, mf_327, mf_328, \
                         mf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -kf_325[k]
                   + f_0 * mf_325[k];

        t_326[k] = -kf_326[k]
                   + f_0 * mf_326[k];

        t_327[k] = -kf_327[k]
                   + f_0 * mf_327[k];

        t_328[k] = -kf_328[k]
                   + f_0 * mf_328[k];

        t_329[k] = -kf_329[k]
                   + f_0 * mf_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, kf_330, kf_331, kf_332, kf_333, \
                         kf_334, mf_330, mf_331, mf_332, mf_333, \
                         mf_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -kf_330[k]
                   + f_0 * mf_330[k];

        t_331[k] = -kf_331[k]
                   + f_0 * mf_331[k];

        t_332[k] = -kf_332[k]
                   + f_0 * mf_332[k];

        t_333[k] = -kf_333[k]
                   + f_0 * mf_333[k];

        t_334[k] = -kf_334[k]
                   + f_0 * mf_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, kf_335, kf_336, kf_337, kf_338, \
                         kf_339, mf_335, mf_336, mf_337, mf_338, \
                         mf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -kf_335[k]
                   + f_0 * mf_335[k];

        t_336[k] = -kf_336[k]
                   + f_0 * mf_336[k];

        t_337[k] = -kf_337[k]
                   + f_0 * mf_337[k];

        t_338[k] = -kf_338[k]
                   + f_0 * mf_338[k];

        t_339[k] = -kf_339[k]
                   + f_0 * mf_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, kf_340, kf_341, kf_342, kf_343, \
                         kf_344, mf_340, mf_341, mf_342, mf_343, \
                         mf_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -kf_340[k]
                   + f_0 * mf_340[k];

        t_341[k] = -kf_341[k]
                   + f_0 * mf_341[k];

        t_342[k] = -kf_342[k]
                   + f_0 * mf_342[k];

        t_343[k] = -kf_343[k]
                   + f_0 * mf_343[k];

        t_344[k] = -kf_344[k]
                   + f_0 * mf_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, kf_345, kf_346, kf_347, kf_348, \
                         kf_349, mf_345, mf_346, mf_347, mf_348, \
                         mf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -kf_345[k]
                   + f_0 * mf_345[k];

        t_346[k] = -kf_346[k]
                   + f_0 * mf_346[k];

        t_347[k] = -kf_347[k]
                   + f_0 * mf_347[k];

        t_348[k] = -kf_348[k]
                   + f_0 * mf_348[k];

        t_349[k] = -kf_349[k]
                   + f_0 * mf_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, kf_350, kf_351, kf_352, kf_353, \
                         kf_354, mf_350, mf_351, mf_352, mf_353, \
                         mf_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -kf_350[k]
                   + f_0 * mf_350[k];

        t_351[k] = -kf_351[k]
                   + f_0 * mf_351[k];

        t_352[k] = -kf_352[k]
                   + f_0 * mf_352[k];

        t_353[k] = -kf_353[k]
                   + f_0 * mf_353[k];

        t_354[k] = -kf_354[k]
                   + f_0 * mf_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, kf_355, kf_356, kf_357, kf_358, \
                         kf_359, mf_355, mf_356, mf_357, mf_358, \
                         mf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -kf_355[k]
                   + f_0 * mf_355[k];

        t_356[k] = -kf_356[k]
                   + f_0 * mf_356[k];

        t_357[k] = -kf_357[k]
                   + f_0 * mf_357[k];

        t_358[k] = -kf_358[k]
                   + f_0 * mf_358[k];

        t_359[k] = -kf_359[k]
                   + f_0 * mf_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, t_365, t_366, t_367, mf_360, \
                         mf_361, mf_362, mf_363, mf_364, mf_365, mf_366, \
                         mf_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_0 * mf_360[k];

        t_361[k] = f_0 * mf_361[k];

        t_362[k] = f_0 * mf_362[k];

        t_363[k] = f_0 * mf_363[k];

        t_364[k] = f_0 * mf_364[k];

        t_365[k] = f_0 * mf_365[k];

        t_366[k] = f_0 * mf_366[k];

        t_367[k] = f_0 * mf_367[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, t_373, t_374, t_375, mf_368, \
                         mf_369, mf_370, mf_371, mf_372, mf_373, mf_374, \
                         mf_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_0 * mf_368[k];

        t_369[k] = f_0 * mf_369[k];

        t_370[k] = f_0 * mf_370[k];

        t_371[k] = f_0 * mf_371[k];

        t_372[k] = f_0 * mf_372[k];

        t_373[k] = f_0 * mf_373[k];

        t_374[k] = f_0 * mf_374[k];

        t_375[k] = f_0 * mf_375[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, t_380, t_381, t_382, t_383, mf_376, \
                         mf_377, mf_378, mf_379, mf_380, mf_381, mf_382, \
                         mf_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_0 * mf_376[k];

        t_377[k] = f_0 * mf_377[k];

        t_378[k] = f_0 * mf_378[k];

        t_379[k] = f_0 * mf_379[k];

        t_380[k] = f_0 * mf_380[k];

        t_381[k] = f_0 * mf_381[k];

        t_382[k] = f_0 * mf_382[k];

        t_383[k] = f_0 * mf_383[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, t_388, t_389, t_390, t_391, mf_384, \
                         mf_385, mf_386, mf_387, mf_388, mf_389, mf_390, \
                         mf_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_0 * mf_384[k];

        t_385[k] = f_0 * mf_385[k];

        t_386[k] = f_0 * mf_386[k];

        t_387[k] = f_0 * mf_387[k];

        t_388[k] = f_0 * mf_388[k];

        t_389[k] = f_0 * mf_389[k];

        t_390[k] = f_0 * mf_390[k];

        t_391[k] = f_0 * mf_391[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, t_397, t_398, t_399, mf_392, \
                         mf_393, mf_394, mf_395, mf_396, mf_397, mf_398, \
                         mf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_0 * mf_392[k];

        t_393[k] = f_0 * mf_393[k];

        t_394[k] = f_0 * mf_394[k];

        t_395[k] = f_0 * mf_395[k];

        t_396[k] = f_0 * mf_396[k];

        t_397[k] = f_0 * mf_397[k];

        t_398[k] = f_0 * mf_398[k];

        t_399[k] = f_0 * mf_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, t_405, t_406, t_407, mf_400, \
                         mf_401, mf_402, mf_403, mf_404, mf_405, mf_406, \
                         mf_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_0 * mf_400[k];

        t_401[k] = f_0 * mf_401[k];

        t_402[k] = f_0 * mf_402[k];

        t_403[k] = f_0 * mf_403[k];

        t_404[k] = f_0 * mf_404[k];

        t_405[k] = f_0 * mf_405[k];

        t_406[k] = f_0 * mf_406[k];

        t_407[k] = f_0 * mf_407[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, t_413, t_414, t_415, mf_408, \
                         mf_409, mf_410, mf_411, mf_412, mf_413, mf_414, \
                         mf_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_0 * mf_408[k];

        t_409[k] = f_0 * mf_409[k];

        t_410[k] = f_0 * mf_410[k];

        t_411[k] = f_0 * mf_411[k];

        t_412[k] = f_0 * mf_412[k];

        t_413[k] = f_0 * mf_413[k];

        t_414[k] = f_0 * mf_414[k];

        t_415[k] = f_0 * mf_415[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, t_420, t_421, t_422, t_423, mf_416, \
                         mf_417, mf_418, mf_419, mf_420, mf_421, mf_422, \
                         mf_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_0 * mf_416[k];

        t_417[k] = f_0 * mf_417[k];

        t_418[k] = f_0 * mf_418[k];

        t_419[k] = f_0 * mf_419[k];

        t_420[k] = f_0 * mf_420[k];

        t_421[k] = f_0 * mf_421[k];

        t_422[k] = f_0 * mf_422[k];

        t_423[k] = f_0 * mf_423[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, t_428, t_429, t_430, t_431, mf_424, \
                         mf_425, mf_426, mf_427, mf_428, mf_429, mf_430, \
                         mf_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_0 * mf_424[k];

        t_425[k] = f_0 * mf_425[k];

        t_426[k] = f_0 * mf_426[k];

        t_427[k] = f_0 * mf_427[k];

        t_428[k] = f_0 * mf_428[k];

        t_429[k] = f_0 * mf_429[k];

        t_430[k] = f_0 * mf_430[k];

        t_431[k] = f_0 * mf_431[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, t_436, t_437, t_438, t_439, mf_432, \
                         mf_433, mf_434, mf_435, mf_436, mf_437, mf_438, \
                         mf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_0 * mf_432[k];

        t_433[k] = f_0 * mf_433[k];

        t_434[k] = f_0 * mf_434[k];

        t_435[k] = f_0 * mf_435[k];

        t_436[k] = f_0 * mf_436[k];

        t_437[k] = f_0 * mf_437[k];

        t_438[k] = f_0 * mf_438[k];

        t_439[k] = f_0 * mf_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, t_445, t_446, t_447, mf_440, \
                         mf_441, mf_442, mf_443, mf_444, mf_445, mf_446, \
                         mf_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_0 * mf_440[k];

        t_441[k] = f_0 * mf_441[k];

        t_442[k] = f_0 * mf_442[k];

        t_443[k] = f_0 * mf_443[k];

        t_444[k] = f_0 * mf_444[k];

        t_445[k] = f_0 * mf_445[k];

        t_446[k] = f_0 * mf_446[k];

        t_447[k] = f_0 * mf_447[k];
    }

#pragma omp simd aligned(t_448, t_449, mf_448, mf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = f_0 * mf_448[k];

        t_449[k] = f_0 * mf_449[k];
    }
}

auto
compute_prim_geom_10_lf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t kf, const size_t mf,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_lf_electron_repulsion_0_piece0(buffer, target, kf, mf, ncols, alpha);

    compute_prim_geom_10_lf_electron_repulsion_0_piece1(buffer, target, kf, mf, ncols, alpha);

    compute_prim_geom_10_lf_electron_repulsion_0_piece2(buffer, target, kf, mf, ncols, alpha);
}

static auto
compute_prim_geom_10_lf_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kf, const size_t mf,
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

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_1 = buffer.data(kf + 1);
    const auto *kf_2 = buffer.data(kf + 2);
    const auto *kf_3 = buffer.data(kf + 3);
    const auto *kf_4 = buffer.data(kf + 4);
    const auto *kf_5 = buffer.data(kf + 5);
    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_12 = buffer.data(kf + 12);
    const auto *kf_13 = buffer.data(kf + 13);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_15 = buffer.data(kf + 15);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_17 = buffer.data(kf + 17);
    const auto *kf_18 = buffer.data(kf + 18);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_21 = buffer.data(kf + 21);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_24 = buffer.data(kf + 24);
    const auto *kf_25 = buffer.data(kf + 25);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_37 = buffer.data(kf + 37);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_49 = buffer.data(kf + 49);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_51 = buffer.data(kf + 51);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_53 = buffer.data(kf + 53);
    const auto *kf_54 = buffer.data(kf + 54);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_62 = buffer.data(kf + 62);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_65 = buffer.data(kf + 65);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_73 = buffer.data(kf + 73);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_81 = buffer.data(kf + 81);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_85 = buffer.data(kf + 85);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_93 = buffer.data(kf + 93);
    const auto *kf_94 = buffer.data(kf + 94);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_97 = buffer.data(kf + 97);
    const auto *kf_98 = buffer.data(kf + 98);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_101 = buffer.data(kf + 101);
    const auto *kf_102 = buffer.data(kf + 102);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_104 = buffer.data(kf + 104);
    const auto *kf_105 = buffer.data(kf + 105);
    const auto *kf_106 = buffer.data(kf + 106);
    const auto *kf_107 = buffer.data(kf + 107);
    const auto *kf_108 = buffer.data(kf + 108);
    const auto *kf_109 = buffer.data(kf + 109);
    const auto *kf_110 = buffer.data(kf + 110);
    const auto *kf_111 = buffer.data(kf + 111);
    const auto *kf_112 = buffer.data(kf + 112);
    const auto *kf_113 = buffer.data(kf + 113);
    const auto *kf_114 = buffer.data(kf + 114);
    const auto *kf_115 = buffer.data(kf + 115);
    const auto *kf_116 = buffer.data(kf + 116);

    const auto *mf_10 = buffer.data(mf + 10);
    const auto *mf_11 = buffer.data(mf + 11);
    const auto *mf_12 = buffer.data(mf + 12);
    const auto *mf_13 = buffer.data(mf + 13);
    const auto *mf_14 = buffer.data(mf + 14);
    const auto *mf_15 = buffer.data(mf + 15);
    const auto *mf_16 = buffer.data(mf + 16);
    const auto *mf_17 = buffer.data(mf + 17);
    const auto *mf_18 = buffer.data(mf + 18);
    const auto *mf_19 = buffer.data(mf + 19);
    const auto *mf_30 = buffer.data(mf + 30);
    const auto *mf_31 = buffer.data(mf + 31);
    const auto *mf_32 = buffer.data(mf + 32);
    const auto *mf_33 = buffer.data(mf + 33);
    const auto *mf_34 = buffer.data(mf + 34);
    const auto *mf_35 = buffer.data(mf + 35);
    const auto *mf_36 = buffer.data(mf + 36);
    const auto *mf_37 = buffer.data(mf + 37);
    const auto *mf_38 = buffer.data(mf + 38);
    const auto *mf_39 = buffer.data(mf + 39);
    const auto *mf_40 = buffer.data(mf + 40);
    const auto *mf_41 = buffer.data(mf + 41);
    const auto *mf_42 = buffer.data(mf + 42);
    const auto *mf_43 = buffer.data(mf + 43);
    const auto *mf_44 = buffer.data(mf + 44);
    const auto *mf_45 = buffer.data(mf + 45);
    const auto *mf_46 = buffer.data(mf + 46);
    const auto *mf_47 = buffer.data(mf + 47);
    const auto *mf_48 = buffer.data(mf + 48);
    const auto *mf_49 = buffer.data(mf + 49);
    const auto *mf_60 = buffer.data(mf + 60);
    const auto *mf_61 = buffer.data(mf + 61);
    const auto *mf_62 = buffer.data(mf + 62);
    const auto *mf_63 = buffer.data(mf + 63);
    const auto *mf_64 = buffer.data(mf + 64);
    const auto *mf_65 = buffer.data(mf + 65);
    const auto *mf_66 = buffer.data(mf + 66);
    const auto *mf_67 = buffer.data(mf + 67);
    const auto *mf_68 = buffer.data(mf + 68);
    const auto *mf_69 = buffer.data(mf + 69);
    const auto *mf_70 = buffer.data(mf + 70);
    const auto *mf_71 = buffer.data(mf + 71);
    const auto *mf_72 = buffer.data(mf + 72);
    const auto *mf_73 = buffer.data(mf + 73);
    const auto *mf_74 = buffer.data(mf + 74);
    const auto *mf_75 = buffer.data(mf + 75);
    const auto *mf_76 = buffer.data(mf + 76);
    const auto *mf_77 = buffer.data(mf + 77);
    const auto *mf_78 = buffer.data(mf + 78);
    const auto *mf_79 = buffer.data(mf + 79);
    const auto *mf_80 = buffer.data(mf + 80);
    const auto *mf_81 = buffer.data(mf + 81);
    const auto *mf_82 = buffer.data(mf + 82);
    const auto *mf_83 = buffer.data(mf + 83);
    const auto *mf_84 = buffer.data(mf + 84);
    const auto *mf_85 = buffer.data(mf + 85);
    const auto *mf_86 = buffer.data(mf + 86);
    const auto *mf_87 = buffer.data(mf + 87);
    const auto *mf_88 = buffer.data(mf + 88);
    const auto *mf_89 = buffer.data(mf + 89);
    const auto *mf_100 = buffer.data(mf + 100);
    const auto *mf_101 = buffer.data(mf + 101);
    const auto *mf_102 = buffer.data(mf + 102);
    const auto *mf_103 = buffer.data(mf + 103);
    const auto *mf_104 = buffer.data(mf + 104);
    const auto *mf_105 = buffer.data(mf + 105);
    const auto *mf_106 = buffer.data(mf + 106);
    const auto *mf_107 = buffer.data(mf + 107);
    const auto *mf_108 = buffer.data(mf + 108);
    const auto *mf_109 = buffer.data(mf + 109);
    const auto *mf_110 = buffer.data(mf + 110);
    const auto *mf_111 = buffer.data(mf + 111);
    const auto *mf_112 = buffer.data(mf + 112);
    const auto *mf_113 = buffer.data(mf + 113);
    const auto *mf_114 = buffer.data(mf + 114);
    const auto *mf_115 = buffer.data(mf + 115);
    const auto *mf_116 = buffer.data(mf + 116);
    const auto *mf_117 = buffer.data(mf + 117);
    const auto *mf_118 = buffer.data(mf + 118);
    const auto *mf_119 = buffer.data(mf + 119);
    const auto *mf_120 = buffer.data(mf + 120);
    const auto *mf_121 = buffer.data(mf + 121);
    const auto *mf_122 = buffer.data(mf + 122);
    const auto *mf_123 = buffer.data(mf + 123);
    const auto *mf_124 = buffer.data(mf + 124);
    const auto *mf_125 = buffer.data(mf + 125);
    const auto *mf_126 = buffer.data(mf + 126);
    const auto *mf_127 = buffer.data(mf + 127);
    const auto *mf_128 = buffer.data(mf + 128);
    const auto *mf_129 = buffer.data(mf + 129);
    const auto *mf_130 = buffer.data(mf + 130);
    const auto *mf_131 = buffer.data(mf + 131);
    const auto *mf_132 = buffer.data(mf + 132);
    const auto *mf_133 = buffer.data(mf + 133);
    const auto *mf_134 = buffer.data(mf + 134);
    const auto *mf_135 = buffer.data(mf + 135);
    const auto *mf_136 = buffer.data(mf + 136);
    const auto *mf_137 = buffer.data(mf + 137);
    const auto *mf_138 = buffer.data(mf + 138);
    const auto *mf_139 = buffer.data(mf + 139);
    const auto *mf_150 = buffer.data(mf + 150);
    const auto *mf_151 = buffer.data(mf + 151);
    const auto *mf_152 = buffer.data(mf + 152);
    const auto *mf_153 = buffer.data(mf + 153);
    const auto *mf_154 = buffer.data(mf + 154);
    const auto *mf_155 = buffer.data(mf + 155);
    const auto *mf_156 = buffer.data(mf + 156);
    const auto *mf_157 = buffer.data(mf + 157);
    const auto *mf_158 = buffer.data(mf + 158);
    const auto *mf_159 = buffer.data(mf + 159);
    const auto *mf_160 = buffer.data(mf + 160);
    const auto *mf_161 = buffer.data(mf + 161);
    const auto *mf_162 = buffer.data(mf + 162);
    const auto *mf_163 = buffer.data(mf + 163);
    const auto *mf_164 = buffer.data(mf + 164);
    const auto *mf_165 = buffer.data(mf + 165);
    const auto *mf_166 = buffer.data(mf + 166);
    const auto *mf_167 = buffer.data(mf + 167);
    const auto *mf_168 = buffer.data(mf + 168);
    const auto *mf_169 = buffer.data(mf + 169);
    const auto *mf_170 = buffer.data(mf + 170);
    const auto *mf_171 = buffer.data(mf + 171);
    const auto *mf_172 = buffer.data(mf + 172);
    const auto *mf_173 = buffer.data(mf + 173);
    const auto *mf_174 = buffer.data(mf + 174);
    const auto *mf_175 = buffer.data(mf + 175);
    const auto *mf_176 = buffer.data(mf + 176);
    const auto *mf_177 = buffer.data(mf + 177);
    const auto *mf_178 = buffer.data(mf + 178);
    const auto *mf_179 = buffer.data(mf + 179);
    const auto *mf_180 = buffer.data(mf + 180);
    const auto *mf_181 = buffer.data(mf + 181);
    const auto *mf_182 = buffer.data(mf + 182);
    const auto *mf_183 = buffer.data(mf + 183);
    const auto *mf_184 = buffer.data(mf + 184);
    const auto *mf_185 = buffer.data(mf + 185);
    const auto *mf_186 = buffer.data(mf + 186);
    const auto *mf_187 = buffer.data(mf + 187);
    const auto *mf_188 = buffer.data(mf + 188);
    const auto *mf_189 = buffer.data(mf + 189);
    const auto *mf_190 = buffer.data(mf + 190);
    const auto *mf_191 = buffer.data(mf + 191);
    const auto *mf_192 = buffer.data(mf + 192);
    const auto *mf_193 = buffer.data(mf + 193);
    const auto *mf_194 = buffer.data(mf + 194);
    const auto *mf_195 = buffer.data(mf + 195);
    const auto *mf_196 = buffer.data(mf + 196);
    const auto *mf_197 = buffer.data(mf + 197);
    const auto *mf_198 = buffer.data(mf + 198);
    const auto *mf_199 = buffer.data(mf + 199);
    const auto *mf_210 = buffer.data(mf + 210);
    const auto *mf_211 = buffer.data(mf + 211);
    const auto *mf_212 = buffer.data(mf + 212);
    const auto *mf_213 = buffer.data(mf + 213);
    const auto *mf_214 = buffer.data(mf + 214);
    const auto *mf_215 = buffer.data(mf + 215);
    const auto *mf_216 = buffer.data(mf + 216);
    const auto *mf_217 = buffer.data(mf + 217);
    const auto *mf_218 = buffer.data(mf + 218);
    const auto *mf_219 = buffer.data(mf + 219);
    const auto *mf_220 = buffer.data(mf + 220);
    const auto *mf_221 = buffer.data(mf + 221);
    const auto *mf_222 = buffer.data(mf + 222);
    const auto *mf_223 = buffer.data(mf + 223);
    const auto *mf_224 = buffer.data(mf + 224);
    const auto *mf_225 = buffer.data(mf + 225);
    const auto *mf_226 = buffer.data(mf + 226);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, mf_10, mf_11, mf_12, mf_13, \
                         mf_14, mf_15, mf_16, mf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * mf_10[k];

        t_1[k] = f_0 * mf_11[k];

        t_2[k] = f_0 * mf_12[k];

        t_3[k] = f_0 * mf_13[k];

        t_4[k] = f_0 * mf_14[k];

        t_5[k] = f_0 * mf_15[k];

        t_6[k] = f_0 * mf_16[k];

        t_7[k] = f_0 * mf_17[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, kf_0, kf_1, kf_2, kf_3, mf_18, \
                         mf_19, mf_30, mf_31, mf_32, mf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * mf_18[k];

        t_9[k] = f_0 * mf_19[k];

        t_10[k] = -kf_0[k]
                  + f_0 * mf_30[k];

        t_11[k] = -kf_1[k]
                  + f_0 * mf_31[k];

        t_12[k] = -kf_2[k]
                  + f_0 * mf_32[k];

        t_13[k] = -kf_3[k]
                  + f_0 * mf_33[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, kf_4, kf_5, kf_6, kf_7, kf_8, mf_34, \
                         mf_35, mf_36, mf_37, mf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -kf_4[k]
                  + f_0 * mf_34[k];

        t_15[k] = -kf_5[k]
                  + f_0 * mf_35[k];

        t_16[k] = -kf_6[k]
                  + f_0 * mf_36[k];

        t_17[k] = -kf_7[k]
                  + f_0 * mf_37[k];

        t_18[k] = -kf_8[k]
                  + f_0 * mf_38[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, t_25, kf_9, mf_39, mf_40, mf_41, \
                         mf_42, mf_43, mf_44, mf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -kf_9[k]
                  + f_0 * mf_39[k];

        t_20[k] = f_0 * mf_40[k];

        t_21[k] = f_0 * mf_41[k];

        t_22[k] = f_0 * mf_42[k];

        t_23[k] = f_0 * mf_43[k];

        t_24[k] = f_0 * mf_44[k];

        t_25[k] = f_0 * mf_45[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, kf_10, kf_11, mf_46, mf_47, \
                         mf_48, mf_49, mf_60, mf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * mf_46[k];

        t_27[k] = f_0 * mf_47[k];

        t_28[k] = f_0 * mf_48[k];

        t_29[k] = f_0 * mf_49[k];

        t_30[k] = -2.0 * kf_10[k]
                  + f_0 * mf_60[k];

        t_31[k] = -2.0 * kf_11[k]
                  + f_0 * mf_61[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, kf_12, kf_13, kf_14, kf_15, kf_16, \
                         mf_62, mf_63, mf_64, mf_65, mf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -2.0 * kf_12[k]
                  + f_0 * mf_62[k];

        t_33[k] = -2.0 * kf_13[k]
                  + f_0 * mf_63[k];

        t_34[k] = -2.0 * kf_14[k]
                  + f_0 * mf_64[k];

        t_35[k] = -2.0 * kf_15[k]
                  + f_0 * mf_65[k];

        t_36[k] = -2.0 * kf_16[k]
                  + f_0 * mf_66[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, kf_17, kf_18, kf_19, kf_20, kf_21, \
                         mf_67, mf_68, mf_69, mf_70, mf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -2.0 * kf_17[k]
                  + f_0 * mf_67[k];

        t_38[k] = -2.0 * kf_18[k]
                  + f_0 * mf_68[k];

        t_39[k] = -2.0 * kf_19[k]
                  + f_0 * mf_69[k];

        t_40[k] = -kf_20[k]
                  + f_0 * mf_70[k];

        t_41[k] = -kf_21[k]
                  + f_0 * mf_71[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, kf_22, kf_23, kf_24, kf_25, kf_26, \
                         mf_72, mf_73, mf_74, mf_75, mf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -kf_22[k]
                  + f_0 * mf_72[k];

        t_43[k] = -kf_23[k]
                  + f_0 * mf_73[k];

        t_44[k] = -kf_24[k]
                  + f_0 * mf_74[k];

        t_45[k] = -kf_25[k]
                  + f_0 * mf_75[k];

        t_46[k] = -kf_26[k]
                  + f_0 * mf_76[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, kf_27, kf_28, kf_29, mf_77, \
                         mf_78, mf_79, mf_80, mf_81, mf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -kf_27[k]
                  + f_0 * mf_77[k];

        t_48[k] = -kf_28[k]
                  + f_0 * mf_78[k];

        t_49[k] = -kf_29[k]
                  + f_0 * mf_79[k];

        t_50[k] = f_0 * mf_80[k];

        t_51[k] = f_0 * mf_81[k];

        t_52[k] = f_0 * mf_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, t_59, mf_83, mf_84, mf_85, mf_86, \
                         mf_87, mf_88, mf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * mf_83[k];

        t_54[k] = f_0 * mf_84[k];

        t_55[k] = f_0 * mf_85[k];

        t_56[k] = f_0 * mf_86[k];

        t_57[k] = f_0 * mf_87[k];

        t_58[k] = f_0 * mf_88[k];

        t_59[k] = f_0 * mf_89[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, kf_30, kf_31, kf_32, kf_33, kf_34, \
                         mf_100, mf_101, mf_102, mf_103, mf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * kf_30[k]
                  + f_0 * mf_100[k];

        t_61[k] = -3.0 * kf_31[k]
                  + f_0 * mf_101[k];

        t_62[k] = -3.0 * kf_32[k]
                  + f_0 * mf_102[k];

        t_63[k] = -3.0 * kf_33[k]
                  + f_0 * mf_103[k];

        t_64[k] = -3.0 * kf_34[k]
                  + f_0 * mf_104[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, kf_35, kf_36, kf_37, kf_38, kf_39, \
                         mf_105, mf_106, mf_107, mf_108, mf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * kf_35[k]
                  + f_0 * mf_105[k];

        t_66[k] = -3.0 * kf_36[k]
                  + f_0 * mf_106[k];

        t_67[k] = -3.0 * kf_37[k]
                  + f_0 * mf_107[k];

        t_68[k] = -3.0 * kf_38[k]
                  + f_0 * mf_108[k];

        t_69[k] = -3.0 * kf_39[k]
                  + f_0 * mf_109[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, kf_40, kf_41, kf_42, kf_43, kf_44, \
                         mf_110, mf_111, mf_112, mf_113, mf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * kf_40[k]
                  + f_0 * mf_110[k];

        t_71[k] = -2.0 * kf_41[k]
                  + f_0 * mf_111[k];

        t_72[k] = -2.0 * kf_42[k]
                  + f_0 * mf_112[k];

        t_73[k] = -2.0 * kf_43[k]
                  + f_0 * mf_113[k];

        t_74[k] = -2.0 * kf_44[k]
                  + f_0 * mf_114[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, kf_45, kf_46, kf_47, kf_48, kf_49, \
                         mf_115, mf_116, mf_117, mf_118, mf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * kf_45[k]
                  + f_0 * mf_115[k];

        t_76[k] = -2.0 * kf_46[k]
                  + f_0 * mf_116[k];

        t_77[k] = -2.0 * kf_47[k]
                  + f_0 * mf_117[k];

        t_78[k] = -2.0 * kf_48[k]
                  + f_0 * mf_118[k];

        t_79[k] = -2.0 * kf_49[k]
                  + f_0 * mf_119[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, kf_50, kf_51, kf_52, kf_53, kf_54, \
                         mf_120, mf_121, mf_122, mf_123, mf_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -kf_50[k]
                  + f_0 * mf_120[k];

        t_81[k] = -kf_51[k]
                  + f_0 * mf_121[k];

        t_82[k] = -kf_52[k]
                  + f_0 * mf_122[k];

        t_83[k] = -kf_53[k]
                  + f_0 * mf_123[k];

        t_84[k] = -kf_54[k]
                  + f_0 * mf_124[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, kf_55, kf_56, kf_57, kf_58, kf_59, \
                         mf_125, mf_126, mf_127, mf_128, mf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -kf_55[k]
                  + f_0 * mf_125[k];

        t_86[k] = -kf_56[k]
                  + f_0 * mf_126[k];

        t_87[k] = -kf_57[k]
                  + f_0 * mf_127[k];

        t_88[k] = -kf_58[k]
                  + f_0 * mf_128[k];

        t_89[k] = -kf_59[k]
                  + f_0 * mf_129[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, mf_130, mf_131, \
                         mf_132, mf_133, mf_134, mf_135, mf_136, \
                         mf_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * mf_130[k];

        t_91[k] = f_0 * mf_131[k];

        t_92[k] = f_0 * mf_132[k];

        t_93[k] = f_0 * mf_133[k];

        t_94[k] = f_0 * mf_134[k];

        t_95[k] = f_0 * mf_135[k];

        t_96[k] = f_0 * mf_136[k];

        t_97[k] = f_0 * mf_137[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, kf_60, kf_61, kf_62, kf_63, \
                         mf_138, mf_139, mf_150, mf_151, mf_152, \
                         mf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * mf_138[k];

        t_99[k] = f_0 * mf_139[k];

        t_100[k] = -4.0 * kf_60[k]
                   + f_0 * mf_150[k];

        t_101[k] = -4.0 * kf_61[k]
                   + f_0 * mf_151[k];

        t_102[k] = -4.0 * kf_62[k]
                   + f_0 * mf_152[k];

        t_103[k] = -4.0 * kf_63[k]
                   + f_0 * mf_153[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, kf_64, kf_65, kf_66, kf_67, kf_68, \
                         mf_154, mf_155, mf_156, mf_157, mf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -4.0 * kf_64[k]
                   + f_0 * mf_154[k];

        t_105[k] = -4.0 * kf_65[k]
                   + f_0 * mf_155[k];

        t_106[k] = -4.0 * kf_66[k]
                   + f_0 * mf_156[k];

        t_107[k] = -4.0 * kf_67[k]
                   + f_0 * mf_157[k];

        t_108[k] = -4.0 * kf_68[k]
                   + f_0 * mf_158[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, kf_69, kf_70, kf_71, kf_72, kf_73, \
                         mf_159, mf_160, mf_161, mf_162, mf_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -4.0 * kf_69[k]
                   + f_0 * mf_159[k];

        t_110[k] = -3.0 * kf_70[k]
                   + f_0 * mf_160[k];

        t_111[k] = -3.0 * kf_71[k]
                   + f_0 * mf_161[k];

        t_112[k] = -3.0 * kf_72[k]
                   + f_0 * mf_162[k];

        t_113[k] = -3.0 * kf_73[k]
                   + f_0 * mf_163[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, kf_74, kf_75, kf_76, kf_77, kf_78, \
                         mf_164, mf_165, mf_166, mf_167, mf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -3.0 * kf_74[k]
                   + f_0 * mf_164[k];

        t_115[k] = -3.0 * kf_75[k]
                   + f_0 * mf_165[k];

        t_116[k] = -3.0 * kf_76[k]
                   + f_0 * mf_166[k];

        t_117[k] = -3.0 * kf_77[k]
                   + f_0 * mf_167[k];

        t_118[k] = -3.0 * kf_78[k]
                   + f_0 * mf_168[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, kf_79, kf_80, kf_81, kf_82, kf_83, \
                         mf_169, mf_170, mf_171, mf_172, mf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -3.0 * kf_79[k]
                   + f_0 * mf_169[k];

        t_120[k] = -2.0 * kf_80[k]
                   + f_0 * mf_170[k];

        t_121[k] = -2.0 * kf_81[k]
                   + f_0 * mf_171[k];

        t_122[k] = -2.0 * kf_82[k]
                   + f_0 * mf_172[k];

        t_123[k] = -2.0 * kf_83[k]
                   + f_0 * mf_173[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, kf_84, kf_85, kf_86, kf_87, kf_88, \
                         mf_174, mf_175, mf_176, mf_177, mf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * kf_84[k]
                   + f_0 * mf_174[k];

        t_125[k] = -2.0 * kf_85[k]
                   + f_0 * mf_175[k];

        t_126[k] = -2.0 * kf_86[k]
                   + f_0 * mf_176[k];

        t_127[k] = -2.0 * kf_87[k]
                   + f_0 * mf_177[k];

        t_128[k] = -2.0 * kf_88[k]
                   + f_0 * mf_178[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, kf_89, kf_90, kf_91, kf_92, kf_93, \
                         mf_179, mf_180, mf_181, mf_182, mf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -2.0 * kf_89[k]
                   + f_0 * mf_179[k];

        t_130[k] = -kf_90[k]
                   + f_0 * mf_180[k];

        t_131[k] = -kf_91[k]
                   + f_0 * mf_181[k];

        t_132[k] = -kf_92[k]
                   + f_0 * mf_182[k];

        t_133[k] = -kf_93[k]
                   + f_0 * mf_183[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, kf_94, kf_95, kf_96, kf_97, kf_98, \
                         mf_184, mf_185, mf_186, mf_187, mf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -kf_94[k]
                   + f_0 * mf_184[k];

        t_135[k] = -kf_95[k]
                   + f_0 * mf_185[k];

        t_136[k] = -kf_96[k]
                   + f_0 * mf_186[k];

        t_137[k] = -kf_97[k]
                   + f_0 * mf_187[k];

        t_138[k] = -kf_98[k]
                   + f_0 * mf_188[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, t_145, kf_99, mf_189, \
                         mf_190, mf_191, mf_192, mf_193, mf_194, \
                         mf_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -kf_99[k]
                   + f_0 * mf_189[k];

        t_140[k] = f_0 * mf_190[k];

        t_141[k] = f_0 * mf_191[k];

        t_142[k] = f_0 * mf_192[k];

        t_143[k] = f_0 * mf_193[k];

        t_144[k] = f_0 * mf_194[k];

        t_145[k] = f_0 * mf_195[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, t_151, kf_100, kf_101, mf_196, \
                         mf_197, mf_198, mf_199, mf_210, mf_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_0 * mf_196[k];

        t_147[k] = f_0 * mf_197[k];

        t_148[k] = f_0 * mf_198[k];

        t_149[k] = f_0 * mf_199[k];

        t_150[k] = -5.0 * kf_100[k]
                   + f_0 * mf_210[k];

        t_151[k] = -5.0 * kf_101[k]
                   + f_0 * mf_211[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, kf_102, kf_103, kf_104, kf_105, \
                         kf_106, mf_212, mf_213, mf_214, mf_215, \
                         mf_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -5.0 * kf_102[k]
                   + f_0 * mf_212[k];

        t_153[k] = -5.0 * kf_103[k]
                   + f_0 * mf_213[k];

        t_154[k] = -5.0 * kf_104[k]
                   + f_0 * mf_214[k];

        t_155[k] = -5.0 * kf_105[k]
                   + f_0 * mf_215[k];

        t_156[k] = -5.0 * kf_106[k]
                   + f_0 * mf_216[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, kf_107, kf_108, kf_109, kf_110, \
                         kf_111, mf_217, mf_218, mf_219, mf_220, \
                         mf_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -5.0 * kf_107[k]
                   + f_0 * mf_217[k];

        t_158[k] = -5.0 * kf_108[k]
                   + f_0 * mf_218[k];

        t_159[k] = -5.0 * kf_109[k]
                   + f_0 * mf_219[k];

        t_160[k] = -4.0 * kf_110[k]
                   + f_0 * mf_220[k];

        t_161[k] = -4.0 * kf_111[k]
                   + f_0 * mf_221[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, kf_112, kf_113, kf_114, kf_115, \
                         kf_116, mf_222, mf_223, mf_224, mf_225, \
                         mf_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -4.0 * kf_112[k]
                   + f_0 * mf_222[k];

        t_163[k] = -4.0 * kf_113[k]
                   + f_0 * mf_223[k];

        t_164[k] = -4.0 * kf_114[k]
                   + f_0 * mf_224[k];

        t_165[k] = -4.0 * kf_115[k]
                   + f_0 * mf_225[k];

        t_166[k] = -4.0 * kf_116[k]
                   + f_0 * mf_226[k];
    }
}

static auto
compute_prim_geom_10_lf_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kf, const size_t mf,
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

    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_119 = buffer.data(kf + 119);
    const auto *kf_120 = buffer.data(kf + 120);
    const auto *kf_121 = buffer.data(kf + 121);
    const auto *kf_122 = buffer.data(kf + 122);
    const auto *kf_123 = buffer.data(kf + 123);
    const auto *kf_124 = buffer.data(kf + 124);
    const auto *kf_125 = buffer.data(kf + 125);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_127 = buffer.data(kf + 127);
    const auto *kf_128 = buffer.data(kf + 128);
    const auto *kf_129 = buffer.data(kf + 129);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_131 = buffer.data(kf + 131);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_133 = buffer.data(kf + 133);
    const auto *kf_134 = buffer.data(kf + 134);
    const auto *kf_135 = buffer.data(kf + 135);
    const auto *kf_136 = buffer.data(kf + 136);
    const auto *kf_137 = buffer.data(kf + 137);
    const auto *kf_138 = buffer.data(kf + 138);
    const auto *kf_139 = buffer.data(kf + 139);
    const auto *kf_140 = buffer.data(kf + 140);
    const auto *kf_141 = buffer.data(kf + 141);
    const auto *kf_142 = buffer.data(kf + 142);
    const auto *kf_143 = buffer.data(kf + 143);
    const auto *kf_144 = buffer.data(kf + 144);
    const auto *kf_145 = buffer.data(kf + 145);
    const auto *kf_146 = buffer.data(kf + 146);
    const auto *kf_147 = buffer.data(kf + 147);
    const auto *kf_148 = buffer.data(kf + 148);
    const auto *kf_149 = buffer.data(kf + 149);
    const auto *kf_150 = buffer.data(kf + 150);
    const auto *kf_151 = buffer.data(kf + 151);
    const auto *kf_152 = buffer.data(kf + 152);
    const auto *kf_153 = buffer.data(kf + 153);
    const auto *kf_154 = buffer.data(kf + 154);
    const auto *kf_155 = buffer.data(kf + 155);
    const auto *kf_156 = buffer.data(kf + 156);
    const auto *kf_157 = buffer.data(kf + 157);
    const auto *kf_158 = buffer.data(kf + 158);
    const auto *kf_159 = buffer.data(kf + 159);
    const auto *kf_160 = buffer.data(kf + 160);
    const auto *kf_161 = buffer.data(kf + 161);
    const auto *kf_162 = buffer.data(kf + 162);
    const auto *kf_163 = buffer.data(kf + 163);
    const auto *kf_164 = buffer.data(kf + 164);
    const auto *kf_165 = buffer.data(kf + 165);
    const auto *kf_166 = buffer.data(kf + 166);
    const auto *kf_167 = buffer.data(kf + 167);
    const auto *kf_168 = buffer.data(kf + 168);
    const auto *kf_169 = buffer.data(kf + 169);
    const auto *kf_170 = buffer.data(kf + 170);
    const auto *kf_171 = buffer.data(kf + 171);
    const auto *kf_172 = buffer.data(kf + 172);
    const auto *kf_173 = buffer.data(kf + 173);
    const auto *kf_174 = buffer.data(kf + 174);
    const auto *kf_175 = buffer.data(kf + 175);
    const auto *kf_176 = buffer.data(kf + 176);
    const auto *kf_177 = buffer.data(kf + 177);
    const auto *kf_178 = buffer.data(kf + 178);
    const auto *kf_179 = buffer.data(kf + 179);
    const auto *kf_180 = buffer.data(kf + 180);
    const auto *kf_181 = buffer.data(kf + 181);
    const auto *kf_182 = buffer.data(kf + 182);
    const auto *kf_183 = buffer.data(kf + 183);
    const auto *kf_184 = buffer.data(kf + 184);
    const auto *kf_185 = buffer.data(kf + 185);
    const auto *kf_186 = buffer.data(kf + 186);
    const auto *kf_187 = buffer.data(kf + 187);
    const auto *kf_188 = buffer.data(kf + 188);
    const auto *kf_189 = buffer.data(kf + 189);
    const auto *kf_190 = buffer.data(kf + 190);
    const auto *kf_191 = buffer.data(kf + 191);
    const auto *kf_192 = buffer.data(kf + 192);
    const auto *kf_193 = buffer.data(kf + 193);
    const auto *kf_194 = buffer.data(kf + 194);
    const auto *kf_195 = buffer.data(kf + 195);
    const auto *kf_196 = buffer.data(kf + 196);
    const auto *kf_197 = buffer.data(kf + 197);
    const auto *kf_198 = buffer.data(kf + 198);
    const auto *kf_199 = buffer.data(kf + 199);
    const auto *kf_200 = buffer.data(kf + 200);
    const auto *kf_201 = buffer.data(kf + 201);
    const auto *kf_202 = buffer.data(kf + 202);
    const auto *kf_203 = buffer.data(kf + 203);
    const auto *kf_204 = buffer.data(kf + 204);
    const auto *kf_205 = buffer.data(kf + 205);
    const auto *kf_206 = buffer.data(kf + 206);
    const auto *kf_207 = buffer.data(kf + 207);
    const auto *kf_208 = buffer.data(kf + 208);
    const auto *kf_209 = buffer.data(kf + 209);
    const auto *kf_210 = buffer.data(kf + 210);
    const auto *kf_211 = buffer.data(kf + 211);
    const auto *kf_212 = buffer.data(kf + 212);
    const auto *kf_213 = buffer.data(kf + 213);
    const auto *kf_214 = buffer.data(kf + 214);
    const auto *kf_215 = buffer.data(kf + 215);
    const auto *kf_216 = buffer.data(kf + 216);
    const auto *kf_217 = buffer.data(kf + 217);
    const auto *kf_218 = buffer.data(kf + 218);
    const auto *kf_219 = buffer.data(kf + 219);
    const auto *kf_220 = buffer.data(kf + 220);
    const auto *kf_221 = buffer.data(kf + 221);
    const auto *kf_222 = buffer.data(kf + 222);
    const auto *kf_223 = buffer.data(kf + 223);
    const auto *kf_224 = buffer.data(kf + 224);
    const auto *kf_225 = buffer.data(kf + 225);
    const auto *kf_226 = buffer.data(kf + 226);
    const auto *kf_227 = buffer.data(kf + 227);
    const auto *kf_228 = buffer.data(kf + 228);
    const auto *kf_229 = buffer.data(kf + 229);
    const auto *kf_230 = buffer.data(kf + 230);
    const auto *kf_231 = buffer.data(kf + 231);
    const auto *kf_232 = buffer.data(kf + 232);
    const auto *kf_233 = buffer.data(kf + 233);
    const auto *kf_234 = buffer.data(kf + 234);
    const auto *kf_235 = buffer.data(kf + 235);
    const auto *kf_236 = buffer.data(kf + 236);
    const auto *kf_237 = buffer.data(kf + 237);
    const auto *kf_238 = buffer.data(kf + 238);
    const auto *kf_239 = buffer.data(kf + 239);
    const auto *kf_240 = buffer.data(kf + 240);
    const auto *kf_241 = buffer.data(kf + 241);
    const auto *kf_242 = buffer.data(kf + 242);
    const auto *kf_243 = buffer.data(kf + 243);
    const auto *kf_244 = buffer.data(kf + 244);
    const auto *kf_245 = buffer.data(kf + 245);
    const auto *kf_246 = buffer.data(kf + 246);
    const auto *kf_247 = buffer.data(kf + 247);
    const auto *kf_248 = buffer.data(kf + 248);
    const auto *kf_249 = buffer.data(kf + 249);
    const auto *kf_250 = buffer.data(kf + 250);
    const auto *kf_251 = buffer.data(kf + 251);
    const auto *kf_252 = buffer.data(kf + 252);
    const auto *kf_253 = buffer.data(kf + 253);

    const auto *mf_227 = buffer.data(mf + 227);
    const auto *mf_228 = buffer.data(mf + 228);
    const auto *mf_229 = buffer.data(mf + 229);
    const auto *mf_230 = buffer.data(mf + 230);
    const auto *mf_231 = buffer.data(mf + 231);
    const auto *mf_232 = buffer.data(mf + 232);
    const auto *mf_233 = buffer.data(mf + 233);
    const auto *mf_234 = buffer.data(mf + 234);
    const auto *mf_235 = buffer.data(mf + 235);
    const auto *mf_236 = buffer.data(mf + 236);
    const auto *mf_237 = buffer.data(mf + 237);
    const auto *mf_238 = buffer.data(mf + 238);
    const auto *mf_239 = buffer.data(mf + 239);
    const auto *mf_240 = buffer.data(mf + 240);
    const auto *mf_241 = buffer.data(mf + 241);
    const auto *mf_242 = buffer.data(mf + 242);
    const auto *mf_243 = buffer.data(mf + 243);
    const auto *mf_244 = buffer.data(mf + 244);
    const auto *mf_245 = buffer.data(mf + 245);
    const auto *mf_246 = buffer.data(mf + 246);
    const auto *mf_247 = buffer.data(mf + 247);
    const auto *mf_248 = buffer.data(mf + 248);
    const auto *mf_249 = buffer.data(mf + 249);
    const auto *mf_250 = buffer.data(mf + 250);
    const auto *mf_251 = buffer.data(mf + 251);
    const auto *mf_252 = buffer.data(mf + 252);
    const auto *mf_253 = buffer.data(mf + 253);
    const auto *mf_254 = buffer.data(mf + 254);
    const auto *mf_255 = buffer.data(mf + 255);
    const auto *mf_256 = buffer.data(mf + 256);
    const auto *mf_257 = buffer.data(mf + 257);
    const auto *mf_258 = buffer.data(mf + 258);
    const auto *mf_259 = buffer.data(mf + 259);
    const auto *mf_260 = buffer.data(mf + 260);
    const auto *mf_261 = buffer.data(mf + 261);
    const auto *mf_262 = buffer.data(mf + 262);
    const auto *mf_263 = buffer.data(mf + 263);
    const auto *mf_264 = buffer.data(mf + 264);
    const auto *mf_265 = buffer.data(mf + 265);
    const auto *mf_266 = buffer.data(mf + 266);
    const auto *mf_267 = buffer.data(mf + 267);
    const auto *mf_268 = buffer.data(mf + 268);
    const auto *mf_269 = buffer.data(mf + 269);
    const auto *mf_280 = buffer.data(mf + 280);
    const auto *mf_281 = buffer.data(mf + 281);
    const auto *mf_282 = buffer.data(mf + 282);
    const auto *mf_283 = buffer.data(mf + 283);
    const auto *mf_284 = buffer.data(mf + 284);
    const auto *mf_285 = buffer.data(mf + 285);
    const auto *mf_286 = buffer.data(mf + 286);
    const auto *mf_287 = buffer.data(mf + 287);
    const auto *mf_288 = buffer.data(mf + 288);
    const auto *mf_289 = buffer.data(mf + 289);
    const auto *mf_290 = buffer.data(mf + 290);
    const auto *mf_291 = buffer.data(mf + 291);
    const auto *mf_292 = buffer.data(mf + 292);
    const auto *mf_293 = buffer.data(mf + 293);
    const auto *mf_294 = buffer.data(mf + 294);
    const auto *mf_295 = buffer.data(mf + 295);
    const auto *mf_296 = buffer.data(mf + 296);
    const auto *mf_297 = buffer.data(mf + 297);
    const auto *mf_298 = buffer.data(mf + 298);
    const auto *mf_299 = buffer.data(mf + 299);
    const auto *mf_300 = buffer.data(mf + 300);
    const auto *mf_301 = buffer.data(mf + 301);
    const auto *mf_302 = buffer.data(mf + 302);
    const auto *mf_303 = buffer.data(mf + 303);
    const auto *mf_304 = buffer.data(mf + 304);
    const auto *mf_305 = buffer.data(mf + 305);
    const auto *mf_306 = buffer.data(mf + 306);
    const auto *mf_307 = buffer.data(mf + 307);
    const auto *mf_308 = buffer.data(mf + 308);
    const auto *mf_309 = buffer.data(mf + 309);
    const auto *mf_310 = buffer.data(mf + 310);
    const auto *mf_311 = buffer.data(mf + 311);
    const auto *mf_312 = buffer.data(mf + 312);
    const auto *mf_313 = buffer.data(mf + 313);
    const auto *mf_314 = buffer.data(mf + 314);
    const auto *mf_315 = buffer.data(mf + 315);
    const auto *mf_316 = buffer.data(mf + 316);
    const auto *mf_317 = buffer.data(mf + 317);
    const auto *mf_318 = buffer.data(mf + 318);
    const auto *mf_319 = buffer.data(mf + 319);
    const auto *mf_320 = buffer.data(mf + 320);
    const auto *mf_321 = buffer.data(mf + 321);
    const auto *mf_322 = buffer.data(mf + 322);
    const auto *mf_323 = buffer.data(mf + 323);
    const auto *mf_324 = buffer.data(mf + 324);
    const auto *mf_325 = buffer.data(mf + 325);
    const auto *mf_326 = buffer.data(mf + 326);
    const auto *mf_327 = buffer.data(mf + 327);
    const auto *mf_328 = buffer.data(mf + 328);
    const auto *mf_329 = buffer.data(mf + 329);
    const auto *mf_330 = buffer.data(mf + 330);
    const auto *mf_331 = buffer.data(mf + 331);
    const auto *mf_332 = buffer.data(mf + 332);
    const auto *mf_333 = buffer.data(mf + 333);
    const auto *mf_334 = buffer.data(mf + 334);
    const auto *mf_335 = buffer.data(mf + 335);
    const auto *mf_336 = buffer.data(mf + 336);
    const auto *mf_337 = buffer.data(mf + 337);
    const auto *mf_338 = buffer.data(mf + 338);
    const auto *mf_339 = buffer.data(mf + 339);
    const auto *mf_340 = buffer.data(mf + 340);
    const auto *mf_341 = buffer.data(mf + 341);
    const auto *mf_342 = buffer.data(mf + 342);
    const auto *mf_343 = buffer.data(mf + 343);
    const auto *mf_344 = buffer.data(mf + 344);
    const auto *mf_345 = buffer.data(mf + 345);
    const auto *mf_346 = buffer.data(mf + 346);
    const auto *mf_347 = buffer.data(mf + 347);
    const auto *mf_348 = buffer.data(mf + 348);
    const auto *mf_349 = buffer.data(mf + 349);
    const auto *mf_360 = buffer.data(mf + 360);
    const auto *mf_361 = buffer.data(mf + 361);
    const auto *mf_362 = buffer.data(mf + 362);
    const auto *mf_363 = buffer.data(mf + 363);
    const auto *mf_364 = buffer.data(mf + 364);
    const auto *mf_365 = buffer.data(mf + 365);
    const auto *mf_366 = buffer.data(mf + 366);
    const auto *mf_367 = buffer.data(mf + 367);
    const auto *mf_368 = buffer.data(mf + 368);
    const auto *mf_369 = buffer.data(mf + 369);
    const auto *mf_370 = buffer.data(mf + 370);
    const auto *mf_371 = buffer.data(mf + 371);
    const auto *mf_372 = buffer.data(mf + 372);
    const auto *mf_373 = buffer.data(mf + 373);
    const auto *mf_374 = buffer.data(mf + 374);
    const auto *mf_375 = buffer.data(mf + 375);
    const auto *mf_376 = buffer.data(mf + 376);
    const auto *mf_377 = buffer.data(mf + 377);
    const auto *mf_378 = buffer.data(mf + 378);
    const auto *mf_379 = buffer.data(mf + 379);
    const auto *mf_380 = buffer.data(mf + 380);
    const auto *mf_381 = buffer.data(mf + 381);
    const auto *mf_382 = buffer.data(mf + 382);
    const auto *mf_383 = buffer.data(mf + 383);
    const auto *mf_384 = buffer.data(mf + 384);
    const auto *mf_385 = buffer.data(mf + 385);
    const auto *mf_386 = buffer.data(mf + 386);
    const auto *mf_387 = buffer.data(mf + 387);
    const auto *mf_388 = buffer.data(mf + 388);
    const auto *mf_389 = buffer.data(mf + 389);
    const auto *mf_390 = buffer.data(mf + 390);
    const auto *mf_391 = buffer.data(mf + 391);
    const auto *mf_392 = buffer.data(mf + 392);
    const auto *mf_393 = buffer.data(mf + 393);
    const auto *mf_394 = buffer.data(mf + 394);
    const auto *mf_395 = buffer.data(mf + 395);
    const auto *mf_396 = buffer.data(mf + 396);
    const auto *mf_397 = buffer.data(mf + 397);
    const auto *mf_398 = buffer.data(mf + 398);
    const auto *mf_399 = buffer.data(mf + 399);
    const auto *mf_400 = buffer.data(mf + 400);
    const auto *mf_401 = buffer.data(mf + 401);
    const auto *mf_402 = buffer.data(mf + 402);
    const auto *mf_403 = buffer.data(mf + 403);

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, kf_117, kf_118, kf_119, kf_120, \
                         kf_121, mf_227, mf_228, mf_229, mf_230, \
                         mf_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -4.0 * kf_117[k]
                   + f_0 * mf_227[k];

        t_168[k] = -4.0 * kf_118[k]
                   + f_0 * mf_228[k];

        t_169[k] = -4.0 * kf_119[k]
                   + f_0 * mf_229[k];

        t_170[k] = -3.0 * kf_120[k]
                   + f_0 * mf_230[k];

        t_171[k] = -3.0 * kf_121[k]
                   + f_0 * mf_231[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, kf_122, kf_123, kf_124, kf_125, \
                         kf_126, mf_232, mf_233, mf_234, mf_235, \
                         mf_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -3.0 * kf_122[k]
                   + f_0 * mf_232[k];

        t_173[k] = -3.0 * kf_123[k]
                   + f_0 * mf_233[k];

        t_174[k] = -3.0 * kf_124[k]
                   + f_0 * mf_234[k];

        t_175[k] = -3.0 * kf_125[k]
                   + f_0 * mf_235[k];

        t_176[k] = -3.0 * kf_126[k]
                   + f_0 * mf_236[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, kf_127, kf_128, kf_129, kf_130, \
                         kf_131, mf_237, mf_238, mf_239, mf_240, \
                         mf_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -3.0 * kf_127[k]
                   + f_0 * mf_237[k];

        t_178[k] = -3.0 * kf_128[k]
                   + f_0 * mf_238[k];

        t_179[k] = -3.0 * kf_129[k]
                   + f_0 * mf_239[k];

        t_180[k] = -2.0 * kf_130[k]
                   + f_0 * mf_240[k];

        t_181[k] = -2.0 * kf_131[k]
                   + f_0 * mf_241[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, kf_132, kf_133, kf_134, kf_135, \
                         kf_136, mf_242, mf_243, mf_244, mf_245, \
                         mf_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * kf_132[k]
                   + f_0 * mf_242[k];

        t_183[k] = -2.0 * kf_133[k]
                   + f_0 * mf_243[k];

        t_184[k] = -2.0 * kf_134[k]
                   + f_0 * mf_244[k];

        t_185[k] = -2.0 * kf_135[k]
                   + f_0 * mf_245[k];

        t_186[k] = -2.0 * kf_136[k]
                   + f_0 * mf_246[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, kf_137, kf_138, kf_139, kf_140, \
                         kf_141, mf_247, mf_248, mf_249, mf_250, \
                         mf_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * kf_137[k]
                   + f_0 * mf_247[k];

        t_188[k] = -2.0 * kf_138[k]
                   + f_0 * mf_248[k];

        t_189[k] = -2.0 * kf_139[k]
                   + f_0 * mf_249[k];

        t_190[k] = -kf_140[k]
                   + f_0 * mf_250[k];

        t_191[k] = -kf_141[k]
                   + f_0 * mf_251[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, kf_142, kf_143, kf_144, kf_145, \
                         kf_146, mf_252, mf_253, mf_254, mf_255, \
                         mf_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -kf_142[k]
                   + f_0 * mf_252[k];

        t_193[k] = -kf_143[k]
                   + f_0 * mf_253[k];

        t_194[k] = -kf_144[k]
                   + f_0 * mf_254[k];

        t_195[k] = -kf_145[k]
                   + f_0 * mf_255[k];

        t_196[k] = -kf_146[k]
                   + f_0 * mf_256[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, t_202, kf_147, kf_148, kf_149, \
                         mf_257, mf_258, mf_259, mf_260, mf_261, \
                         mf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -kf_147[k]
                   + f_0 * mf_257[k];

        t_198[k] = -kf_148[k]
                   + f_0 * mf_258[k];

        t_199[k] = -kf_149[k]
                   + f_0 * mf_259[k];

        t_200[k] = f_0 * mf_260[k];

        t_201[k] = f_0 * mf_261[k];

        t_202[k] = f_0 * mf_262[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, t_209, mf_263, mf_264, \
                         mf_265, mf_266, mf_267, mf_268, mf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_0 * mf_263[k];

        t_204[k] = f_0 * mf_264[k];

        t_205[k] = f_0 * mf_265[k];

        t_206[k] = f_0 * mf_266[k];

        t_207[k] = f_0 * mf_267[k];

        t_208[k] = f_0 * mf_268[k];

        t_209[k] = f_0 * mf_269[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, kf_150, kf_151, kf_152, kf_153, \
                         kf_154, mf_280, mf_281, mf_282, mf_283, \
                         mf_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -6.0 * kf_150[k]
                   + f_0 * mf_280[k];

        t_211[k] = -6.0 * kf_151[k]
                   + f_0 * mf_281[k];

        t_212[k] = -6.0 * kf_152[k]
                   + f_0 * mf_282[k];

        t_213[k] = -6.0 * kf_153[k]
                   + f_0 * mf_283[k];

        t_214[k] = -6.0 * kf_154[k]
                   + f_0 * mf_284[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, kf_155, kf_156, kf_157, kf_158, \
                         kf_159, mf_285, mf_286, mf_287, mf_288, \
                         mf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -6.0 * kf_155[k]
                   + f_0 * mf_285[k];

        t_216[k] = -6.0 * kf_156[k]
                   + f_0 * mf_286[k];

        t_217[k] = -6.0 * kf_157[k]
                   + f_0 * mf_287[k];

        t_218[k] = -6.0 * kf_158[k]
                   + f_0 * mf_288[k];

        t_219[k] = -6.0 * kf_159[k]
                   + f_0 * mf_289[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, kf_160, kf_161, kf_162, kf_163, \
                         kf_164, mf_290, mf_291, mf_292, mf_293, \
                         mf_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -5.0 * kf_160[k]
                   + f_0 * mf_290[k];

        t_221[k] = -5.0 * kf_161[k]
                   + f_0 * mf_291[k];

        t_222[k] = -5.0 * kf_162[k]
                   + f_0 * mf_292[k];

        t_223[k] = -5.0 * kf_163[k]
                   + f_0 * mf_293[k];

        t_224[k] = -5.0 * kf_164[k]
                   + f_0 * mf_294[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, kf_165, kf_166, kf_167, kf_168, \
                         kf_169, mf_295, mf_296, mf_297, mf_298, \
                         mf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -5.0 * kf_165[k]
                   + f_0 * mf_295[k];

        t_226[k] = -5.0 * kf_166[k]
                   + f_0 * mf_296[k];

        t_227[k] = -5.0 * kf_167[k]
                   + f_0 * mf_297[k];

        t_228[k] = -5.0 * kf_168[k]
                   + f_0 * mf_298[k];

        t_229[k] = -5.0 * kf_169[k]
                   + f_0 * mf_299[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, kf_170, kf_171, kf_172, kf_173, \
                         kf_174, mf_300, mf_301, mf_302, mf_303, \
                         mf_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -4.0 * kf_170[k]
                   + f_0 * mf_300[k];

        t_231[k] = -4.0 * kf_171[k]
                   + f_0 * mf_301[k];

        t_232[k] = -4.0 * kf_172[k]
                   + f_0 * mf_302[k];

        t_233[k] = -4.0 * kf_173[k]
                   + f_0 * mf_303[k];

        t_234[k] = -4.0 * kf_174[k]
                   + f_0 * mf_304[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, kf_175, kf_176, kf_177, kf_178, \
                         kf_179, mf_305, mf_306, mf_307, mf_308, \
                         mf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -4.0 * kf_175[k]
                   + f_0 * mf_305[k];

        t_236[k] = -4.0 * kf_176[k]
                   + f_0 * mf_306[k];

        t_237[k] = -4.0 * kf_177[k]
                   + f_0 * mf_307[k];

        t_238[k] = -4.0 * kf_178[k]
                   + f_0 * mf_308[k];

        t_239[k] = -4.0 * kf_179[k]
                   + f_0 * mf_309[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, kf_180, kf_181, kf_182, kf_183, \
                         kf_184, mf_310, mf_311, mf_312, mf_313, \
                         mf_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * kf_180[k]
                   + f_0 * mf_310[k];

        t_241[k] = -3.0 * kf_181[k]
                   + f_0 * mf_311[k];

        t_242[k] = -3.0 * kf_182[k]
                   + f_0 * mf_312[k];

        t_243[k] = -3.0 * kf_183[k]
                   + f_0 * mf_313[k];

        t_244[k] = -3.0 * kf_184[k]
                   + f_0 * mf_314[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, kf_185, kf_186, kf_187, kf_188, \
                         kf_189, mf_315, mf_316, mf_317, mf_318, \
                         mf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * kf_185[k]
                   + f_0 * mf_315[k];

        t_246[k] = -3.0 * kf_186[k]
                   + f_0 * mf_316[k];

        t_247[k] = -3.0 * kf_187[k]
                   + f_0 * mf_317[k];

        t_248[k] = -3.0 * kf_188[k]
                   + f_0 * mf_318[k];

        t_249[k] = -3.0 * kf_189[k]
                   + f_0 * mf_319[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, kf_190, kf_191, kf_192, kf_193, \
                         kf_194, mf_320, mf_321, mf_322, mf_323, \
                         mf_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * kf_190[k]
                   + f_0 * mf_320[k];

        t_251[k] = -2.0 * kf_191[k]
                   + f_0 * mf_321[k];

        t_252[k] = -2.0 * kf_192[k]
                   + f_0 * mf_322[k];

        t_253[k] = -2.0 * kf_193[k]
                   + f_0 * mf_323[k];

        t_254[k] = -2.0 * kf_194[k]
                   + f_0 * mf_324[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, kf_195, kf_196, kf_197, kf_198, \
                         kf_199, mf_325, mf_326, mf_327, mf_328, \
                         mf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * kf_195[k]
                   + f_0 * mf_325[k];

        t_256[k] = -2.0 * kf_196[k]
                   + f_0 * mf_326[k];

        t_257[k] = -2.0 * kf_197[k]
                   + f_0 * mf_327[k];

        t_258[k] = -2.0 * kf_198[k]
                   + f_0 * mf_328[k];

        t_259[k] = -2.0 * kf_199[k]
                   + f_0 * mf_329[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, kf_200, kf_201, kf_202, kf_203, \
                         kf_204, mf_330, mf_331, mf_332, mf_333, \
                         mf_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -kf_200[k]
                   + f_0 * mf_330[k];

        t_261[k] = -kf_201[k]
                   + f_0 * mf_331[k];

        t_262[k] = -kf_202[k]
                   + f_0 * mf_332[k];

        t_263[k] = -kf_203[k]
                   + f_0 * mf_333[k];

        t_264[k] = -kf_204[k]
                   + f_0 * mf_334[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, kf_205, kf_206, kf_207, kf_208, \
                         kf_209, mf_335, mf_336, mf_337, mf_338, \
                         mf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -kf_205[k]
                   + f_0 * mf_335[k];

        t_266[k] = -kf_206[k]
                   + f_0 * mf_336[k];

        t_267[k] = -kf_207[k]
                   + f_0 * mf_337[k];

        t_268[k] = -kf_208[k]
                   + f_0 * mf_338[k];

        t_269[k] = -kf_209[k]
                   + f_0 * mf_339[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, t_276, t_277, mf_340, \
                         mf_341, mf_342, mf_343, mf_344, mf_345, mf_346, \
                         mf_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_0 * mf_340[k];

        t_271[k] = f_0 * mf_341[k];

        t_272[k] = f_0 * mf_342[k];

        t_273[k] = f_0 * mf_343[k];

        t_274[k] = f_0 * mf_344[k];

        t_275[k] = f_0 * mf_345[k];

        t_276[k] = f_0 * mf_346[k];

        t_277[k] = f_0 * mf_347[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, t_283, kf_210, kf_211, kf_212, \
                         kf_213, mf_348, mf_349, mf_360, mf_361, mf_362, \
                         mf_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_0 * mf_348[k];

        t_279[k] = f_0 * mf_349[k];

        t_280[k] = -7.0 * kf_210[k]
                   + f_0 * mf_360[k];

        t_281[k] = -7.0 * kf_211[k]
                   + f_0 * mf_361[k];

        t_282[k] = -7.0 * kf_212[k]
                   + f_0 * mf_362[k];

        t_283[k] = -7.0 * kf_213[k]
                   + f_0 * mf_363[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, kf_214, kf_215, kf_216, kf_217, \
                         kf_218, mf_364, mf_365, mf_366, mf_367, \
                         mf_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = -7.0 * kf_214[k]
                   + f_0 * mf_364[k];

        t_285[k] = -7.0 * kf_215[k]
                   + f_0 * mf_365[k];

        t_286[k] = -7.0 * kf_216[k]
                   + f_0 * mf_366[k];

        t_287[k] = -7.0 * kf_217[k]
                   + f_0 * mf_367[k];

        t_288[k] = -7.0 * kf_218[k]
                   + f_0 * mf_368[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, kf_219, kf_220, kf_221, kf_222, \
                         kf_223, mf_369, mf_370, mf_371, mf_372, \
                         mf_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = -7.0 * kf_219[k]
                   + f_0 * mf_369[k];

        t_290[k] = -6.0 * kf_220[k]
                   + f_0 * mf_370[k];

        t_291[k] = -6.0 * kf_221[k]
                   + f_0 * mf_371[k];

        t_292[k] = -6.0 * kf_222[k]
                   + f_0 * mf_372[k];

        t_293[k] = -6.0 * kf_223[k]
                   + f_0 * mf_373[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, kf_224, kf_225, kf_226, kf_227, \
                         kf_228, mf_374, mf_375, mf_376, mf_377, \
                         mf_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = -6.0 * kf_224[k]
                   + f_0 * mf_374[k];

        t_295[k] = -6.0 * kf_225[k]
                   + f_0 * mf_375[k];

        t_296[k] = -6.0 * kf_226[k]
                   + f_0 * mf_376[k];

        t_297[k] = -6.0 * kf_227[k]
                   + f_0 * mf_377[k];

        t_298[k] = -6.0 * kf_228[k]
                   + f_0 * mf_378[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, kf_229, kf_230, kf_231, kf_232, \
                         kf_233, mf_379, mf_380, mf_381, mf_382, \
                         mf_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = -6.0 * kf_229[k]
                   + f_0 * mf_379[k];

        t_300[k] = -5.0 * kf_230[k]
                   + f_0 * mf_380[k];

        t_301[k] = -5.0 * kf_231[k]
                   + f_0 * mf_381[k];

        t_302[k] = -5.0 * kf_232[k]
                   + f_0 * mf_382[k];

        t_303[k] = -5.0 * kf_233[k]
                   + f_0 * mf_383[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, kf_234, kf_235, kf_236, kf_237, \
                         kf_238, mf_384, mf_385, mf_386, mf_387, \
                         mf_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = -5.0 * kf_234[k]
                   + f_0 * mf_384[k];

        t_305[k] = -5.0 * kf_235[k]
                   + f_0 * mf_385[k];

        t_306[k] = -5.0 * kf_236[k]
                   + f_0 * mf_386[k];

        t_307[k] = -5.0 * kf_237[k]
                   + f_0 * mf_387[k];

        t_308[k] = -5.0 * kf_238[k]
                   + f_0 * mf_388[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, kf_239, kf_240, kf_241, kf_242, \
                         kf_243, mf_389, mf_390, mf_391, mf_392, \
                         mf_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = -5.0 * kf_239[k]
                   + f_0 * mf_389[k];

        t_310[k] = -4.0 * kf_240[k]
                   + f_0 * mf_390[k];

        t_311[k] = -4.0 * kf_241[k]
                   + f_0 * mf_391[k];

        t_312[k] = -4.0 * kf_242[k]
                   + f_0 * mf_392[k];

        t_313[k] = -4.0 * kf_243[k]
                   + f_0 * mf_393[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, t_318, kf_244, kf_245, kf_246, kf_247, \
                         kf_248, mf_394, mf_395, mf_396, mf_397, \
                         mf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = -4.0 * kf_244[k]
                   + f_0 * mf_394[k];

        t_315[k] = -4.0 * kf_245[k]
                   + f_0 * mf_395[k];

        t_316[k] = -4.0 * kf_246[k]
                   + f_0 * mf_396[k];

        t_317[k] = -4.0 * kf_247[k]
                   + f_0 * mf_397[k];

        t_318[k] = -4.0 * kf_248[k]
                   + f_0 * mf_398[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, t_323, kf_249, kf_250, kf_251, kf_252, \
                         kf_253, mf_399, mf_400, mf_401, mf_402, \
                         mf_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = -4.0 * kf_249[k]
                   + f_0 * mf_399[k];

        t_320[k] = -3.0 * kf_250[k]
                   + f_0 * mf_400[k];

        t_321[k] = -3.0 * kf_251[k]
                   + f_0 * mf_401[k];

        t_322[k] = -3.0 * kf_252[k]
                   + f_0 * mf_402[k];

        t_323[k] = -3.0 * kf_253[k]
                   + f_0 * mf_403[k];
    }
}

static auto
compute_prim_geom_10_lf_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kf, const size_t mf,
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

    const auto *kf_254 = buffer.data(kf + 254);
    const auto *kf_255 = buffer.data(kf + 255);
    const auto *kf_256 = buffer.data(kf + 256);
    const auto *kf_257 = buffer.data(kf + 257);
    const auto *kf_258 = buffer.data(kf + 258);
    const auto *kf_259 = buffer.data(kf + 259);
    const auto *kf_260 = buffer.data(kf + 260);
    const auto *kf_261 = buffer.data(kf + 261);
    const auto *kf_262 = buffer.data(kf + 262);
    const auto *kf_263 = buffer.data(kf + 263);
    const auto *kf_264 = buffer.data(kf + 264);
    const auto *kf_265 = buffer.data(kf + 265);
    const auto *kf_266 = buffer.data(kf + 266);
    const auto *kf_267 = buffer.data(kf + 267);
    const auto *kf_268 = buffer.data(kf + 268);
    const auto *kf_269 = buffer.data(kf + 269);
    const auto *kf_270 = buffer.data(kf + 270);
    const auto *kf_271 = buffer.data(kf + 271);
    const auto *kf_272 = buffer.data(kf + 272);
    const auto *kf_273 = buffer.data(kf + 273);
    const auto *kf_274 = buffer.data(kf + 274);
    const auto *kf_275 = buffer.data(kf + 275);
    const auto *kf_276 = buffer.data(kf + 276);
    const auto *kf_277 = buffer.data(kf + 277);
    const auto *kf_278 = buffer.data(kf + 278);
    const auto *kf_279 = buffer.data(kf + 279);
    const auto *kf_280 = buffer.data(kf + 280);
    const auto *kf_281 = buffer.data(kf + 281);
    const auto *kf_282 = buffer.data(kf + 282);
    const auto *kf_283 = buffer.data(kf + 283);
    const auto *kf_284 = buffer.data(kf + 284);
    const auto *kf_285 = buffer.data(kf + 285);
    const auto *kf_286 = buffer.data(kf + 286);
    const auto *kf_287 = buffer.data(kf + 287);
    const auto *kf_288 = buffer.data(kf + 288);
    const auto *kf_289 = buffer.data(kf + 289);
    const auto *kf_290 = buffer.data(kf + 290);
    const auto *kf_291 = buffer.data(kf + 291);
    const auto *kf_292 = buffer.data(kf + 292);
    const auto *kf_293 = buffer.data(kf + 293);
    const auto *kf_294 = buffer.data(kf + 294);
    const auto *kf_295 = buffer.data(kf + 295);
    const auto *kf_296 = buffer.data(kf + 296);
    const auto *kf_297 = buffer.data(kf + 297);
    const auto *kf_298 = buffer.data(kf + 298);
    const auto *kf_299 = buffer.data(kf + 299);
    const auto *kf_300 = buffer.data(kf + 300);
    const auto *kf_301 = buffer.data(kf + 301);
    const auto *kf_302 = buffer.data(kf + 302);
    const auto *kf_303 = buffer.data(kf + 303);
    const auto *kf_304 = buffer.data(kf + 304);
    const auto *kf_305 = buffer.data(kf + 305);
    const auto *kf_306 = buffer.data(kf + 306);
    const auto *kf_307 = buffer.data(kf + 307);
    const auto *kf_308 = buffer.data(kf + 308);
    const auto *kf_309 = buffer.data(kf + 309);
    const auto *kf_310 = buffer.data(kf + 310);
    const auto *kf_311 = buffer.data(kf + 311);
    const auto *kf_312 = buffer.data(kf + 312);
    const auto *kf_313 = buffer.data(kf + 313);
    const auto *kf_314 = buffer.data(kf + 314);
    const auto *kf_315 = buffer.data(kf + 315);
    const auto *kf_316 = buffer.data(kf + 316);
    const auto *kf_317 = buffer.data(kf + 317);
    const auto *kf_318 = buffer.data(kf + 318);
    const auto *kf_319 = buffer.data(kf + 319);
    const auto *kf_320 = buffer.data(kf + 320);
    const auto *kf_321 = buffer.data(kf + 321);
    const auto *kf_322 = buffer.data(kf + 322);
    const auto *kf_323 = buffer.data(kf + 323);
    const auto *kf_324 = buffer.data(kf + 324);
    const auto *kf_325 = buffer.data(kf + 325);
    const auto *kf_326 = buffer.data(kf + 326);
    const auto *kf_327 = buffer.data(kf + 327);
    const auto *kf_328 = buffer.data(kf + 328);
    const auto *kf_329 = buffer.data(kf + 329);
    const auto *kf_330 = buffer.data(kf + 330);
    const auto *kf_331 = buffer.data(kf + 331);
    const auto *kf_332 = buffer.data(kf + 332);
    const auto *kf_333 = buffer.data(kf + 333);
    const auto *kf_334 = buffer.data(kf + 334);
    const auto *kf_335 = buffer.data(kf + 335);
    const auto *kf_336 = buffer.data(kf + 336);
    const auto *kf_337 = buffer.data(kf + 337);
    const auto *kf_338 = buffer.data(kf + 338);
    const auto *kf_339 = buffer.data(kf + 339);
    const auto *kf_340 = buffer.data(kf + 340);
    const auto *kf_341 = buffer.data(kf + 341);
    const auto *kf_342 = buffer.data(kf + 342);
    const auto *kf_343 = buffer.data(kf + 343);
    const auto *kf_344 = buffer.data(kf + 344);
    const auto *kf_345 = buffer.data(kf + 345);
    const auto *kf_346 = buffer.data(kf + 346);
    const auto *kf_347 = buffer.data(kf + 347);
    const auto *kf_348 = buffer.data(kf + 348);
    const auto *kf_349 = buffer.data(kf + 349);
    const auto *kf_350 = buffer.data(kf + 350);
    const auto *kf_351 = buffer.data(kf + 351);
    const auto *kf_352 = buffer.data(kf + 352);
    const auto *kf_353 = buffer.data(kf + 353);
    const auto *kf_354 = buffer.data(kf + 354);
    const auto *kf_355 = buffer.data(kf + 355);
    const auto *kf_356 = buffer.data(kf + 356);
    const auto *kf_357 = buffer.data(kf + 357);
    const auto *kf_358 = buffer.data(kf + 358);
    const auto *kf_359 = buffer.data(kf + 359);

    const auto *mf_404 = buffer.data(mf + 404);
    const auto *mf_405 = buffer.data(mf + 405);
    const auto *mf_406 = buffer.data(mf + 406);
    const auto *mf_407 = buffer.data(mf + 407);
    const auto *mf_408 = buffer.data(mf + 408);
    const auto *mf_409 = buffer.data(mf + 409);
    const auto *mf_410 = buffer.data(mf + 410);
    const auto *mf_411 = buffer.data(mf + 411);
    const auto *mf_412 = buffer.data(mf + 412);
    const auto *mf_413 = buffer.data(mf + 413);
    const auto *mf_414 = buffer.data(mf + 414);
    const auto *mf_415 = buffer.data(mf + 415);
    const auto *mf_416 = buffer.data(mf + 416);
    const auto *mf_417 = buffer.data(mf + 417);
    const auto *mf_418 = buffer.data(mf + 418);
    const auto *mf_419 = buffer.data(mf + 419);
    const auto *mf_420 = buffer.data(mf + 420);
    const auto *mf_421 = buffer.data(mf + 421);
    const auto *mf_422 = buffer.data(mf + 422);
    const auto *mf_423 = buffer.data(mf + 423);
    const auto *mf_424 = buffer.data(mf + 424);
    const auto *mf_425 = buffer.data(mf + 425);
    const auto *mf_426 = buffer.data(mf + 426);
    const auto *mf_427 = buffer.data(mf + 427);
    const auto *mf_428 = buffer.data(mf + 428);
    const auto *mf_429 = buffer.data(mf + 429);
    const auto *mf_430 = buffer.data(mf + 430);
    const auto *mf_431 = buffer.data(mf + 431);
    const auto *mf_432 = buffer.data(mf + 432);
    const auto *mf_433 = buffer.data(mf + 433);
    const auto *mf_434 = buffer.data(mf + 434);
    const auto *mf_435 = buffer.data(mf + 435);
    const auto *mf_436 = buffer.data(mf + 436);
    const auto *mf_437 = buffer.data(mf + 437);
    const auto *mf_438 = buffer.data(mf + 438);
    const auto *mf_439 = buffer.data(mf + 439);
    const auto *mf_450 = buffer.data(mf + 450);
    const auto *mf_451 = buffer.data(mf + 451);
    const auto *mf_452 = buffer.data(mf + 452);
    const auto *mf_453 = buffer.data(mf + 453);
    const auto *mf_454 = buffer.data(mf + 454);
    const auto *mf_455 = buffer.data(mf + 455);
    const auto *mf_456 = buffer.data(mf + 456);
    const auto *mf_457 = buffer.data(mf + 457);
    const auto *mf_458 = buffer.data(mf + 458);
    const auto *mf_459 = buffer.data(mf + 459);
    const auto *mf_460 = buffer.data(mf + 460);
    const auto *mf_461 = buffer.data(mf + 461);
    const auto *mf_462 = buffer.data(mf + 462);
    const auto *mf_463 = buffer.data(mf + 463);
    const auto *mf_464 = buffer.data(mf + 464);
    const auto *mf_465 = buffer.data(mf + 465);
    const auto *mf_466 = buffer.data(mf + 466);
    const auto *mf_467 = buffer.data(mf + 467);
    const auto *mf_468 = buffer.data(mf + 468);
    const auto *mf_469 = buffer.data(mf + 469);
    const auto *mf_470 = buffer.data(mf + 470);
    const auto *mf_471 = buffer.data(mf + 471);
    const auto *mf_472 = buffer.data(mf + 472);
    const auto *mf_473 = buffer.data(mf + 473);
    const auto *mf_474 = buffer.data(mf + 474);
    const auto *mf_475 = buffer.data(mf + 475);
    const auto *mf_476 = buffer.data(mf + 476);
    const auto *mf_477 = buffer.data(mf + 477);
    const auto *mf_478 = buffer.data(mf + 478);
    const auto *mf_479 = buffer.data(mf + 479);
    const auto *mf_480 = buffer.data(mf + 480);
    const auto *mf_481 = buffer.data(mf + 481);
    const auto *mf_482 = buffer.data(mf + 482);
    const auto *mf_483 = buffer.data(mf + 483);
    const auto *mf_484 = buffer.data(mf + 484);
    const auto *mf_485 = buffer.data(mf + 485);
    const auto *mf_486 = buffer.data(mf + 486);
    const auto *mf_487 = buffer.data(mf + 487);
    const auto *mf_488 = buffer.data(mf + 488);
    const auto *mf_489 = buffer.data(mf + 489);
    const auto *mf_490 = buffer.data(mf + 490);
    const auto *mf_491 = buffer.data(mf + 491);
    const auto *mf_492 = buffer.data(mf + 492);
    const auto *mf_493 = buffer.data(mf + 493);
    const auto *mf_494 = buffer.data(mf + 494);
    const auto *mf_495 = buffer.data(mf + 495);
    const auto *mf_496 = buffer.data(mf + 496);
    const auto *mf_497 = buffer.data(mf + 497);
    const auto *mf_498 = buffer.data(mf + 498);
    const auto *mf_499 = buffer.data(mf + 499);
    const auto *mf_500 = buffer.data(mf + 500);
    const auto *mf_501 = buffer.data(mf + 501);
    const auto *mf_502 = buffer.data(mf + 502);
    const auto *mf_503 = buffer.data(mf + 503);
    const auto *mf_504 = buffer.data(mf + 504);
    const auto *mf_505 = buffer.data(mf + 505);
    const auto *mf_506 = buffer.data(mf + 506);
    const auto *mf_507 = buffer.data(mf + 507);
    const auto *mf_508 = buffer.data(mf + 508);
    const auto *mf_509 = buffer.data(mf + 509);
    const auto *mf_510 = buffer.data(mf + 510);
    const auto *mf_511 = buffer.data(mf + 511);
    const auto *mf_512 = buffer.data(mf + 512);
    const auto *mf_513 = buffer.data(mf + 513);
    const auto *mf_514 = buffer.data(mf + 514);
    const auto *mf_515 = buffer.data(mf + 515);
    const auto *mf_516 = buffer.data(mf + 516);
    const auto *mf_517 = buffer.data(mf + 517);
    const auto *mf_518 = buffer.data(mf + 518);
    const auto *mf_519 = buffer.data(mf + 519);
    const auto *mf_520 = buffer.data(mf + 520);
    const auto *mf_521 = buffer.data(mf + 521);
    const auto *mf_522 = buffer.data(mf + 522);
    const auto *mf_523 = buffer.data(mf + 523);
    const auto *mf_524 = buffer.data(mf + 524);
    const auto *mf_525 = buffer.data(mf + 525);
    const auto *mf_526 = buffer.data(mf + 526);
    const auto *mf_527 = buffer.data(mf + 527);
    const auto *mf_528 = buffer.data(mf + 528);
    const auto *mf_529 = buffer.data(mf + 529);
    const auto *mf_530 = buffer.data(mf + 530);
    const auto *mf_531 = buffer.data(mf + 531);
    const auto *mf_532 = buffer.data(mf + 532);
    const auto *mf_533 = buffer.data(mf + 533);
    const auto *mf_534 = buffer.data(mf + 534);
    const auto *mf_535 = buffer.data(mf + 535);
    const auto *mf_536 = buffer.data(mf + 536);
    const auto *mf_537 = buffer.data(mf + 537);
    const auto *mf_538 = buffer.data(mf + 538);
    const auto *mf_539 = buffer.data(mf + 539);

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, kf_254, kf_255, kf_256, kf_257, \
                         kf_258, mf_404, mf_405, mf_406, mf_407, \
                         mf_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = -3.0 * kf_254[k]
                   + f_0 * mf_404[k];

        t_325[k] = -3.0 * kf_255[k]
                   + f_0 * mf_405[k];

        t_326[k] = -3.0 * kf_256[k]
                   + f_0 * mf_406[k];

        t_327[k] = -3.0 * kf_257[k]
                   + f_0 * mf_407[k];

        t_328[k] = -3.0 * kf_258[k]
                   + f_0 * mf_408[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, t_333, kf_259, kf_260, kf_261, kf_262, \
                         kf_263, mf_409, mf_410, mf_411, mf_412, \
                         mf_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = -3.0 * kf_259[k]
                   + f_0 * mf_409[k];

        t_330[k] = -2.0 * kf_260[k]
                   + f_0 * mf_410[k];

        t_331[k] = -2.0 * kf_261[k]
                   + f_0 * mf_411[k];

        t_332[k] = -2.0 * kf_262[k]
                   + f_0 * mf_412[k];

        t_333[k] = -2.0 * kf_263[k]
                   + f_0 * mf_413[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, kf_264, kf_265, kf_266, kf_267, \
                         kf_268, mf_414, mf_415, mf_416, mf_417, \
                         mf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = -2.0 * kf_264[k]
                   + f_0 * mf_414[k];

        t_335[k] = -2.0 * kf_265[k]
                   + f_0 * mf_415[k];

        t_336[k] = -2.0 * kf_266[k]
                   + f_0 * mf_416[k];

        t_337[k] = -2.0 * kf_267[k]
                   + f_0 * mf_417[k];

        t_338[k] = -2.0 * kf_268[k]
                   + f_0 * mf_418[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, t_343, kf_269, kf_270, kf_271, kf_272, \
                         kf_273, mf_419, mf_420, mf_421, mf_422, \
                         mf_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = -2.0 * kf_269[k]
                   + f_0 * mf_419[k];

        t_340[k] = -kf_270[k]
                   + f_0 * mf_420[k];

        t_341[k] = -kf_271[k]
                   + f_0 * mf_421[k];

        t_342[k] = -kf_272[k]
                   + f_0 * mf_422[k];

        t_343[k] = -kf_273[k]
                   + f_0 * mf_423[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, kf_274, kf_275, kf_276, kf_277, \
                         kf_278, mf_424, mf_425, mf_426, mf_427, \
                         mf_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = -kf_274[k]
                   + f_0 * mf_424[k];

        t_345[k] = -kf_275[k]
                   + f_0 * mf_425[k];

        t_346[k] = -kf_276[k]
                   + f_0 * mf_426[k];

        t_347[k] = -kf_277[k]
                   + f_0 * mf_427[k];

        t_348[k] = -kf_278[k]
                   + f_0 * mf_428[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, t_353, t_354, t_355, kf_279, mf_429, \
                         mf_430, mf_431, mf_432, mf_433, mf_434, \
                         mf_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = -kf_279[k]
                   + f_0 * mf_429[k];

        t_350[k] = f_0 * mf_430[k];

        t_351[k] = f_0 * mf_431[k];

        t_352[k] = f_0 * mf_432[k];

        t_353[k] = f_0 * mf_433[k];

        t_354[k] = f_0 * mf_434[k];

        t_355[k] = f_0 * mf_435[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, t_360, t_361, kf_280, kf_281, mf_436, \
                         mf_437, mf_438, mf_439, mf_450, mf_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_0 * mf_436[k];

        t_357[k] = f_0 * mf_437[k];

        t_358[k] = f_0 * mf_438[k];

        t_359[k] = f_0 * mf_439[k];

        t_360[k] = -8.0 * kf_280[k]
                   + f_0 * mf_450[k];

        t_361[k] = -8.0 * kf_281[k]
                   + f_0 * mf_451[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, kf_282, kf_283, kf_284, kf_285, \
                         kf_286, mf_452, mf_453, mf_454, mf_455, \
                         mf_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -8.0 * kf_282[k]
                   + f_0 * mf_452[k];

        t_363[k] = -8.0 * kf_283[k]
                   + f_0 * mf_453[k];

        t_364[k] = -8.0 * kf_284[k]
                   + f_0 * mf_454[k];

        t_365[k] = -8.0 * kf_285[k]
                   + f_0 * mf_455[k];

        t_366[k] = -8.0 * kf_286[k]
                   + f_0 * mf_456[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, kf_287, kf_288, kf_289, kf_290, \
                         kf_291, mf_457, mf_458, mf_459, mf_460, \
                         mf_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -8.0 * kf_287[k]
                   + f_0 * mf_457[k];

        t_368[k] = -8.0 * kf_288[k]
                   + f_0 * mf_458[k];

        t_369[k] = -8.0 * kf_289[k]
                   + f_0 * mf_459[k];

        t_370[k] = -7.0 * kf_290[k]
                   + f_0 * mf_460[k];

        t_371[k] = -7.0 * kf_291[k]
                   + f_0 * mf_461[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, kf_292, kf_293, kf_294, kf_295, \
                         kf_296, mf_462, mf_463, mf_464, mf_465, \
                         mf_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -7.0 * kf_292[k]
                   + f_0 * mf_462[k];

        t_373[k] = -7.0 * kf_293[k]
                   + f_0 * mf_463[k];

        t_374[k] = -7.0 * kf_294[k]
                   + f_0 * mf_464[k];

        t_375[k] = -7.0 * kf_295[k]
                   + f_0 * mf_465[k];

        t_376[k] = -7.0 * kf_296[k]
                   + f_0 * mf_466[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, kf_297, kf_298, kf_299, kf_300, \
                         kf_301, mf_467, mf_468, mf_469, mf_470, \
                         mf_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -7.0 * kf_297[k]
                   + f_0 * mf_467[k];

        t_378[k] = -7.0 * kf_298[k]
                   + f_0 * mf_468[k];

        t_379[k] = -7.0 * kf_299[k]
                   + f_0 * mf_469[k];

        t_380[k] = -6.0 * kf_300[k]
                   + f_0 * mf_470[k];

        t_381[k] = -6.0 * kf_301[k]
                   + f_0 * mf_471[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, kf_302, kf_303, kf_304, kf_305, \
                         kf_306, mf_472, mf_473, mf_474, mf_475, \
                         mf_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -6.0 * kf_302[k]
                   + f_0 * mf_472[k];

        t_383[k] = -6.0 * kf_303[k]
                   + f_0 * mf_473[k];

        t_384[k] = -6.0 * kf_304[k]
                   + f_0 * mf_474[k];

        t_385[k] = -6.0 * kf_305[k]
                   + f_0 * mf_475[k];

        t_386[k] = -6.0 * kf_306[k]
                   + f_0 * mf_476[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, kf_307, kf_308, kf_309, kf_310, \
                         kf_311, mf_477, mf_478, mf_479, mf_480, \
                         mf_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -6.0 * kf_307[k]
                   + f_0 * mf_477[k];

        t_388[k] = -6.0 * kf_308[k]
                   + f_0 * mf_478[k];

        t_389[k] = -6.0 * kf_309[k]
                   + f_0 * mf_479[k];

        t_390[k] = -5.0 * kf_310[k]
                   + f_0 * mf_480[k];

        t_391[k] = -5.0 * kf_311[k]
                   + f_0 * mf_481[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, kf_312, kf_313, kf_314, kf_315, \
                         kf_316, mf_482, mf_483, mf_484, mf_485, \
                         mf_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -5.0 * kf_312[k]
                   + f_0 * mf_482[k];

        t_393[k] = -5.0 * kf_313[k]
                   + f_0 * mf_483[k];

        t_394[k] = -5.0 * kf_314[k]
                   + f_0 * mf_484[k];

        t_395[k] = -5.0 * kf_315[k]
                   + f_0 * mf_485[k];

        t_396[k] = -5.0 * kf_316[k]
                   + f_0 * mf_486[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, kf_317, kf_318, kf_319, kf_320, \
                         kf_321, mf_487, mf_488, mf_489, mf_490, \
                         mf_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -5.0 * kf_317[k]
                   + f_0 * mf_487[k];

        t_398[k] = -5.0 * kf_318[k]
                   + f_0 * mf_488[k];

        t_399[k] = -5.0 * kf_319[k]
                   + f_0 * mf_489[k];

        t_400[k] = -4.0 * kf_320[k]
                   + f_0 * mf_490[k];

        t_401[k] = -4.0 * kf_321[k]
                   + f_0 * mf_491[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, kf_322, kf_323, kf_324, kf_325, \
                         kf_326, mf_492, mf_493, mf_494, mf_495, \
                         mf_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -4.0 * kf_322[k]
                   + f_0 * mf_492[k];

        t_403[k] = -4.0 * kf_323[k]
                   + f_0 * mf_493[k];

        t_404[k] = -4.0 * kf_324[k]
                   + f_0 * mf_494[k];

        t_405[k] = -4.0 * kf_325[k]
                   + f_0 * mf_495[k];

        t_406[k] = -4.0 * kf_326[k]
                   + f_0 * mf_496[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, kf_327, kf_328, kf_329, kf_330, \
                         kf_331, mf_497, mf_498, mf_499, mf_500, \
                         mf_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -4.0 * kf_327[k]
                   + f_0 * mf_497[k];

        t_408[k] = -4.0 * kf_328[k]
                   + f_0 * mf_498[k];

        t_409[k] = -4.0 * kf_329[k]
                   + f_0 * mf_499[k];

        t_410[k] = -3.0 * kf_330[k]
                   + f_0 * mf_500[k];

        t_411[k] = -3.0 * kf_331[k]
                   + f_0 * mf_501[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, kf_332, kf_333, kf_334, kf_335, \
                         kf_336, mf_502, mf_503, mf_504, mf_505, \
                         mf_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -3.0 * kf_332[k]
                   + f_0 * mf_502[k];

        t_413[k] = -3.0 * kf_333[k]
                   + f_0 * mf_503[k];

        t_414[k] = -3.0 * kf_334[k]
                   + f_0 * mf_504[k];

        t_415[k] = -3.0 * kf_335[k]
                   + f_0 * mf_505[k];

        t_416[k] = -3.0 * kf_336[k]
                   + f_0 * mf_506[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, kf_337, kf_338, kf_339, kf_340, \
                         kf_341, mf_507, mf_508, mf_509, mf_510, \
                         mf_511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -3.0 * kf_337[k]
                   + f_0 * mf_507[k];

        t_418[k] = -3.0 * kf_338[k]
                   + f_0 * mf_508[k];

        t_419[k] = -3.0 * kf_339[k]
                   + f_0 * mf_509[k];

        t_420[k] = -2.0 * kf_340[k]
                   + f_0 * mf_510[k];

        t_421[k] = -2.0 * kf_341[k]
                   + f_0 * mf_511[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, t_426, kf_342, kf_343, kf_344, kf_345, \
                         kf_346, mf_512, mf_513, mf_514, mf_515, \
                         mf_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = -2.0 * kf_342[k]
                   + f_0 * mf_512[k];

        t_423[k] = -2.0 * kf_343[k]
                   + f_0 * mf_513[k];

        t_424[k] = -2.0 * kf_344[k]
                   + f_0 * mf_514[k];

        t_425[k] = -2.0 * kf_345[k]
                   + f_0 * mf_515[k];

        t_426[k] = -2.0 * kf_346[k]
                   + f_0 * mf_516[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, t_431, kf_347, kf_348, kf_349, kf_350, \
                         kf_351, mf_517, mf_518, mf_519, mf_520, \
                         mf_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = -2.0 * kf_347[k]
                   + f_0 * mf_517[k];

        t_428[k] = -2.0 * kf_348[k]
                   + f_0 * mf_518[k];

        t_429[k] = -2.0 * kf_349[k]
                   + f_0 * mf_519[k];

        t_430[k] = -kf_350[k]
                   + f_0 * mf_520[k];

        t_431[k] = -kf_351[k]
                   + f_0 * mf_521[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, t_436, kf_352, kf_353, kf_354, kf_355, \
                         kf_356, mf_522, mf_523, mf_524, mf_525, \
                         mf_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = -kf_352[k]
                   + f_0 * mf_522[k];

        t_433[k] = -kf_353[k]
                   + f_0 * mf_523[k];

        t_434[k] = -kf_354[k]
                   + f_0 * mf_524[k];

        t_435[k] = -kf_355[k]
                   + f_0 * mf_525[k];

        t_436[k] = -kf_356[k]
                   + f_0 * mf_526[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, t_441, t_442, kf_357, kf_358, kf_359, \
                         mf_527, mf_528, mf_529, mf_530, mf_531, \
                         mf_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = -kf_357[k]
                   + f_0 * mf_527[k];

        t_438[k] = -kf_358[k]
                   + f_0 * mf_528[k];

        t_439[k] = -kf_359[k]
                   + f_0 * mf_529[k];

        t_440[k] = f_0 * mf_530[k];

        t_441[k] = f_0 * mf_531[k];

        t_442[k] = f_0 * mf_532[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, t_447, t_448, t_449, mf_533, mf_534, \
                         mf_535, mf_536, mf_537, mf_538, mf_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_0 * mf_533[k];

        t_444[k] = f_0 * mf_534[k];

        t_445[k] = f_0 * mf_535[k];

        t_446[k] = f_0 * mf_536[k];

        t_447[k] = f_0 * mf_537[k];

        t_448[k] = f_0 * mf_538[k];

        t_449[k] = f_0 * mf_539[k];
    }
}

auto
compute_prim_geom_10_lf_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t kf, const size_t mf,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_lf_electron_repulsion_1_piece0(buffer, target, kf, mf, ncols, alpha);

    compute_prim_geom_10_lf_electron_repulsion_1_piece1(buffer, target, kf, mf, ncols, alpha);

    compute_prim_geom_10_lf_electron_repulsion_1_piece2(buffer, target, kf, mf, ncols, alpha);
}

static auto
compute_prim_geom_10_lf_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kf, const size_t mf,
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

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_1 = buffer.data(kf + 1);
    const auto *kf_2 = buffer.data(kf + 2);
    const auto *kf_3 = buffer.data(kf + 3);
    const auto *kf_4 = buffer.data(kf + 4);
    const auto *kf_5 = buffer.data(kf + 5);
    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_12 = buffer.data(kf + 12);
    const auto *kf_13 = buffer.data(kf + 13);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_15 = buffer.data(kf + 15);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_17 = buffer.data(kf + 17);
    const auto *kf_18 = buffer.data(kf + 18);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_21 = buffer.data(kf + 21);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_24 = buffer.data(kf + 24);
    const auto *kf_25 = buffer.data(kf + 25);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_37 = buffer.data(kf + 37);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_49 = buffer.data(kf + 49);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_51 = buffer.data(kf + 51);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_53 = buffer.data(kf + 53);
    const auto *kf_54 = buffer.data(kf + 54);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_62 = buffer.data(kf + 62);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_65 = buffer.data(kf + 65);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_73 = buffer.data(kf + 73);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_81 = buffer.data(kf + 81);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_85 = buffer.data(kf + 85);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_93 = buffer.data(kf + 93);
    const auto *kf_94 = buffer.data(kf + 94);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_97 = buffer.data(kf + 97);
    const auto *kf_98 = buffer.data(kf + 98);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_101 = buffer.data(kf + 101);
    const auto *kf_102 = buffer.data(kf + 102);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_104 = buffer.data(kf + 104);
    const auto *kf_105 = buffer.data(kf + 105);
    const auto *kf_106 = buffer.data(kf + 106);
    const auto *kf_107 = buffer.data(kf + 107);
    const auto *kf_108 = buffer.data(kf + 108);
    const auto *kf_109 = buffer.data(kf + 109);

    const auto *mf_20 = buffer.data(mf + 20);
    const auto *mf_21 = buffer.data(mf + 21);
    const auto *mf_22 = buffer.data(mf + 22);
    const auto *mf_23 = buffer.data(mf + 23);
    const auto *mf_24 = buffer.data(mf + 24);
    const auto *mf_25 = buffer.data(mf + 25);
    const auto *mf_26 = buffer.data(mf + 26);
    const auto *mf_27 = buffer.data(mf + 27);
    const auto *mf_28 = buffer.data(mf + 28);
    const auto *mf_29 = buffer.data(mf + 29);
    const auto *mf_40 = buffer.data(mf + 40);
    const auto *mf_41 = buffer.data(mf + 41);
    const auto *mf_42 = buffer.data(mf + 42);
    const auto *mf_43 = buffer.data(mf + 43);
    const auto *mf_44 = buffer.data(mf + 44);
    const auto *mf_45 = buffer.data(mf + 45);
    const auto *mf_46 = buffer.data(mf + 46);
    const auto *mf_47 = buffer.data(mf + 47);
    const auto *mf_48 = buffer.data(mf + 48);
    const auto *mf_49 = buffer.data(mf + 49);
    const auto *mf_50 = buffer.data(mf + 50);
    const auto *mf_51 = buffer.data(mf + 51);
    const auto *mf_52 = buffer.data(mf + 52);
    const auto *mf_53 = buffer.data(mf + 53);
    const auto *mf_54 = buffer.data(mf + 54);
    const auto *mf_55 = buffer.data(mf + 55);
    const auto *mf_56 = buffer.data(mf + 56);
    const auto *mf_57 = buffer.data(mf + 57);
    const auto *mf_58 = buffer.data(mf + 58);
    const auto *mf_59 = buffer.data(mf + 59);
    const auto *mf_70 = buffer.data(mf + 70);
    const auto *mf_71 = buffer.data(mf + 71);
    const auto *mf_72 = buffer.data(mf + 72);
    const auto *mf_73 = buffer.data(mf + 73);
    const auto *mf_74 = buffer.data(mf + 74);
    const auto *mf_75 = buffer.data(mf + 75);
    const auto *mf_76 = buffer.data(mf + 76);
    const auto *mf_77 = buffer.data(mf + 77);
    const auto *mf_78 = buffer.data(mf + 78);
    const auto *mf_79 = buffer.data(mf + 79);
    const auto *mf_80 = buffer.data(mf + 80);
    const auto *mf_81 = buffer.data(mf + 81);
    const auto *mf_82 = buffer.data(mf + 82);
    const auto *mf_83 = buffer.data(mf + 83);
    const auto *mf_84 = buffer.data(mf + 84);
    const auto *mf_85 = buffer.data(mf + 85);
    const auto *mf_86 = buffer.data(mf + 86);
    const auto *mf_87 = buffer.data(mf + 87);
    const auto *mf_88 = buffer.data(mf + 88);
    const auto *mf_89 = buffer.data(mf + 89);
    const auto *mf_90 = buffer.data(mf + 90);
    const auto *mf_91 = buffer.data(mf + 91);
    const auto *mf_92 = buffer.data(mf + 92);
    const auto *mf_93 = buffer.data(mf + 93);
    const auto *mf_94 = buffer.data(mf + 94);
    const auto *mf_95 = buffer.data(mf + 95);
    const auto *mf_96 = buffer.data(mf + 96);
    const auto *mf_97 = buffer.data(mf + 97);
    const auto *mf_98 = buffer.data(mf + 98);
    const auto *mf_99 = buffer.data(mf + 99);
    const auto *mf_110 = buffer.data(mf + 110);
    const auto *mf_111 = buffer.data(mf + 111);
    const auto *mf_112 = buffer.data(mf + 112);
    const auto *mf_113 = buffer.data(mf + 113);
    const auto *mf_114 = buffer.data(mf + 114);
    const auto *mf_115 = buffer.data(mf + 115);
    const auto *mf_116 = buffer.data(mf + 116);
    const auto *mf_117 = buffer.data(mf + 117);
    const auto *mf_118 = buffer.data(mf + 118);
    const auto *mf_119 = buffer.data(mf + 119);
    const auto *mf_120 = buffer.data(mf + 120);
    const auto *mf_121 = buffer.data(mf + 121);
    const auto *mf_122 = buffer.data(mf + 122);
    const auto *mf_123 = buffer.data(mf + 123);
    const auto *mf_124 = buffer.data(mf + 124);
    const auto *mf_125 = buffer.data(mf + 125);
    const auto *mf_126 = buffer.data(mf + 126);
    const auto *mf_127 = buffer.data(mf + 127);
    const auto *mf_128 = buffer.data(mf + 128);
    const auto *mf_129 = buffer.data(mf + 129);
    const auto *mf_130 = buffer.data(mf + 130);
    const auto *mf_131 = buffer.data(mf + 131);
    const auto *mf_132 = buffer.data(mf + 132);
    const auto *mf_133 = buffer.data(mf + 133);
    const auto *mf_134 = buffer.data(mf + 134);
    const auto *mf_135 = buffer.data(mf + 135);
    const auto *mf_136 = buffer.data(mf + 136);
    const auto *mf_137 = buffer.data(mf + 137);
    const auto *mf_138 = buffer.data(mf + 138);
    const auto *mf_139 = buffer.data(mf + 139);
    const auto *mf_140 = buffer.data(mf + 140);
    const auto *mf_141 = buffer.data(mf + 141);
    const auto *mf_142 = buffer.data(mf + 142);
    const auto *mf_143 = buffer.data(mf + 143);
    const auto *mf_144 = buffer.data(mf + 144);
    const auto *mf_145 = buffer.data(mf + 145);
    const auto *mf_146 = buffer.data(mf + 146);
    const auto *mf_147 = buffer.data(mf + 147);
    const auto *mf_148 = buffer.data(mf + 148);
    const auto *mf_149 = buffer.data(mf + 149);
    const auto *mf_160 = buffer.data(mf + 160);
    const auto *mf_161 = buffer.data(mf + 161);
    const auto *mf_162 = buffer.data(mf + 162);
    const auto *mf_163 = buffer.data(mf + 163);
    const auto *mf_164 = buffer.data(mf + 164);
    const auto *mf_165 = buffer.data(mf + 165);
    const auto *mf_166 = buffer.data(mf + 166);
    const auto *mf_167 = buffer.data(mf + 167);
    const auto *mf_168 = buffer.data(mf + 168);
    const auto *mf_169 = buffer.data(mf + 169);
    const auto *mf_170 = buffer.data(mf + 170);
    const auto *mf_171 = buffer.data(mf + 171);
    const auto *mf_172 = buffer.data(mf + 172);
    const auto *mf_173 = buffer.data(mf + 173);
    const auto *mf_174 = buffer.data(mf + 174);
    const auto *mf_175 = buffer.data(mf + 175);
    const auto *mf_176 = buffer.data(mf + 176);
    const auto *mf_177 = buffer.data(mf + 177);
    const auto *mf_178 = buffer.data(mf + 178);
    const auto *mf_179 = buffer.data(mf + 179);
    const auto *mf_180 = buffer.data(mf + 180);
    const auto *mf_181 = buffer.data(mf + 181);
    const auto *mf_182 = buffer.data(mf + 182);
    const auto *mf_183 = buffer.data(mf + 183);
    const auto *mf_184 = buffer.data(mf + 184);
    const auto *mf_185 = buffer.data(mf + 185);
    const auto *mf_186 = buffer.data(mf + 186);
    const auto *mf_187 = buffer.data(mf + 187);
    const auto *mf_188 = buffer.data(mf + 188);
    const auto *mf_189 = buffer.data(mf + 189);
    const auto *mf_190 = buffer.data(mf + 190);
    const auto *mf_191 = buffer.data(mf + 191);
    const auto *mf_192 = buffer.data(mf + 192);
    const auto *mf_193 = buffer.data(mf + 193);
    const auto *mf_194 = buffer.data(mf + 194);
    const auto *mf_195 = buffer.data(mf + 195);
    const auto *mf_196 = buffer.data(mf + 196);
    const auto *mf_197 = buffer.data(mf + 197);
    const auto *mf_198 = buffer.data(mf + 198);
    const auto *mf_199 = buffer.data(mf + 199);
    const auto *mf_200 = buffer.data(mf + 200);
    const auto *mf_201 = buffer.data(mf + 201);
    const auto *mf_202 = buffer.data(mf + 202);
    const auto *mf_203 = buffer.data(mf + 203);
    const auto *mf_204 = buffer.data(mf + 204);
    const auto *mf_205 = buffer.data(mf + 205);
    const auto *mf_206 = buffer.data(mf + 206);
    const auto *mf_207 = buffer.data(mf + 207);
    const auto *mf_208 = buffer.data(mf + 208);
    const auto *mf_209 = buffer.data(mf + 209);
    const auto *mf_220 = buffer.data(mf + 220);
    const auto *mf_221 = buffer.data(mf + 221);
    const auto *mf_222 = buffer.data(mf + 222);
    const auto *mf_223 = buffer.data(mf + 223);
    const auto *mf_224 = buffer.data(mf + 224);
    const auto *mf_225 = buffer.data(mf + 225);
    const auto *mf_226 = buffer.data(mf + 226);
    const auto *mf_227 = buffer.data(mf + 227);
    const auto *mf_228 = buffer.data(mf + 228);
    const auto *mf_229 = buffer.data(mf + 229);
    const auto *mf_230 = buffer.data(mf + 230);
    const auto *mf_231 = buffer.data(mf + 231);
    const auto *mf_232 = buffer.data(mf + 232);
    const auto *mf_233 = buffer.data(mf + 233);
    const auto *mf_234 = buffer.data(mf + 234);
    const auto *mf_235 = buffer.data(mf + 235);
    const auto *mf_236 = buffer.data(mf + 236);
    const auto *mf_237 = buffer.data(mf + 237);
    const auto *mf_238 = buffer.data(mf + 238);
    const auto *mf_239 = buffer.data(mf + 239);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, mf_20, mf_21, mf_22, mf_23, \
                         mf_24, mf_25, mf_26, mf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * mf_20[k];

        t_1[k] = f_0 * mf_21[k];

        t_2[k] = f_0 * mf_22[k];

        t_3[k] = f_0 * mf_23[k];

        t_4[k] = f_0 * mf_24[k];

        t_5[k] = f_0 * mf_25[k];

        t_6[k] = f_0 * mf_26[k];

        t_7[k] = f_0 * mf_27[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, mf_28, mf_29, mf_40, \
                         mf_41, mf_42, mf_43, mf_44, mf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * mf_28[k];

        t_9[k] = f_0 * mf_29[k];

        t_10[k] = f_0 * mf_40[k];

        t_11[k] = f_0 * mf_41[k];

        t_12[k] = f_0 * mf_42[k];

        t_13[k] = f_0 * mf_43[k];

        t_14[k] = f_0 * mf_44[k];

        t_15[k] = f_0 * mf_45[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, kf_0, kf_1, mf_46, mf_47, mf_48, \
                         mf_49, mf_50, mf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * mf_46[k];

        t_17[k] = f_0 * mf_47[k];

        t_18[k] = f_0 * mf_48[k];

        t_19[k] = f_0 * mf_49[k];

        t_20[k] = -kf_0[k]
                  + f_0 * mf_50[k];

        t_21[k] = -kf_1[k]
                  + f_0 * mf_51[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, kf_2, kf_3, kf_4, kf_5, kf_6, mf_52, \
                         mf_53, mf_54, mf_55, mf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -kf_2[k]
                  + f_0 * mf_52[k];

        t_23[k] = -kf_3[k]
                  + f_0 * mf_53[k];

        t_24[k] = -kf_4[k]
                  + f_0 * mf_54[k];

        t_25[k] = -kf_5[k]
                  + f_0 * mf_55[k];

        t_26[k] = -kf_6[k]
                  + f_0 * mf_56[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, kf_7, kf_8, kf_9, mf_57, mf_58, \
                         mf_59, mf_70, mf_71, mf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -kf_7[k]
                  + f_0 * mf_57[k];

        t_28[k] = -kf_8[k]
                  + f_0 * mf_58[k];

        t_29[k] = -kf_9[k]
                  + f_0 * mf_59[k];

        t_30[k] = f_0 * mf_70[k];

        t_31[k] = f_0 * mf_71[k];

        t_32[k] = f_0 * mf_72[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, t_39, mf_73, mf_74, mf_75, mf_76, \
                         mf_77, mf_78, mf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * mf_73[k];

        t_34[k] = f_0 * mf_74[k];

        t_35[k] = f_0 * mf_75[k];

        t_36[k] = f_0 * mf_76[k];

        t_37[k] = f_0 * mf_77[k];

        t_38[k] = f_0 * mf_78[k];

        t_39[k] = f_0 * mf_79[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, kf_10, kf_11, kf_12, kf_13, kf_14, \
                         mf_80, mf_81, mf_82, mf_83, mf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -kf_10[k]
                  + f_0 * mf_80[k];

        t_41[k] = -kf_11[k]
                  + f_0 * mf_81[k];

        t_42[k] = -kf_12[k]
                  + f_0 * mf_82[k];

        t_43[k] = -kf_13[k]
                  + f_0 * mf_83[k];

        t_44[k] = -kf_14[k]
                  + f_0 * mf_84[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, kf_15, kf_16, kf_17, kf_18, kf_19, \
                         mf_85, mf_86, mf_87, mf_88, mf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -kf_15[k]
                  + f_0 * mf_85[k];

        t_46[k] = -kf_16[k]
                  + f_0 * mf_86[k];

        t_47[k] = -kf_17[k]
                  + f_0 * mf_87[k];

        t_48[k] = -kf_18[k]
                  + f_0 * mf_88[k];

        t_49[k] = -kf_19[k]
                  + f_0 * mf_89[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, kf_20, kf_21, kf_22, kf_23, kf_24, \
                         mf_90, mf_91, mf_92, mf_93, mf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * kf_20[k]
                  + f_0 * mf_90[k];

        t_51[k] = -2.0 * kf_21[k]
                  + f_0 * mf_91[k];

        t_52[k] = -2.0 * kf_22[k]
                  + f_0 * mf_92[k];

        t_53[k] = -2.0 * kf_23[k]
                  + f_0 * mf_93[k];

        t_54[k] = -2.0 * kf_24[k]
                  + f_0 * mf_94[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, kf_25, kf_26, kf_27, kf_28, kf_29, \
                         mf_95, mf_96, mf_97, mf_98, mf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * kf_25[k]
                  + f_0 * mf_95[k];

        t_56[k] = -2.0 * kf_26[k]
                  + f_0 * mf_96[k];

        t_57[k] = -2.0 * kf_27[k]
                  + f_0 * mf_97[k];

        t_58[k] = -2.0 * kf_28[k]
                  + f_0 * mf_98[k];

        t_59[k] = -2.0 * kf_29[k]
                  + f_0 * mf_99[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, t_66, t_67, mf_110, mf_111, \
                         mf_112, mf_113, mf_114, mf_115, mf_116, \
                         mf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * mf_110[k];

        t_61[k] = f_0 * mf_111[k];

        t_62[k] = f_0 * mf_112[k];

        t_63[k] = f_0 * mf_113[k];

        t_64[k] = f_0 * mf_114[k];

        t_65[k] = f_0 * mf_115[k];

        t_66[k] = f_0 * mf_116[k];

        t_67[k] = f_0 * mf_117[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, kf_30, kf_31, kf_32, kf_33, \
                         mf_118, mf_119, mf_120, mf_121, mf_122, \
                         mf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * mf_118[k];

        t_69[k] = f_0 * mf_119[k];

        t_70[k] = -kf_30[k]
                  + f_0 * mf_120[k];

        t_71[k] = -kf_31[k]
                  + f_0 * mf_121[k];

        t_72[k] = -kf_32[k]
                  + f_0 * mf_122[k];

        t_73[k] = -kf_33[k]
                  + f_0 * mf_123[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, kf_34, kf_35, kf_36, kf_37, kf_38, \
                         mf_124, mf_125, mf_126, mf_127, mf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = -kf_34[k]
                  + f_0 * mf_124[k];

        t_75[k] = -kf_35[k]
                  + f_0 * mf_125[k];

        t_76[k] = -kf_36[k]
                  + f_0 * mf_126[k];

        t_77[k] = -kf_37[k]
                  + f_0 * mf_127[k];

        t_78[k] = -kf_38[k]
                  + f_0 * mf_128[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, kf_39, kf_40, kf_41, kf_42, kf_43, \
                         mf_129, mf_130, mf_131, mf_132, mf_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -kf_39[k]
                  + f_0 * mf_129[k];

        t_80[k] = -2.0 * kf_40[k]
                  + f_0 * mf_130[k];

        t_81[k] = -2.0 * kf_41[k]
                  + f_0 * mf_131[k];

        t_82[k] = -2.0 * kf_42[k]
                  + f_0 * mf_132[k];

        t_83[k] = -2.0 * kf_43[k]
                  + f_0 * mf_133[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, kf_44, kf_45, kf_46, kf_47, kf_48, \
                         mf_134, mf_135, mf_136, mf_137, mf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -2.0 * kf_44[k]
                  + f_0 * mf_134[k];

        t_85[k] = -2.0 * kf_45[k]
                  + f_0 * mf_135[k];

        t_86[k] = -2.0 * kf_46[k]
                  + f_0 * mf_136[k];

        t_87[k] = -2.0 * kf_47[k]
                  + f_0 * mf_137[k];

        t_88[k] = -2.0 * kf_48[k]
                  + f_0 * mf_138[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, kf_49, kf_50, kf_51, kf_52, kf_53, \
                         mf_139, mf_140, mf_141, mf_142, mf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -2.0 * kf_49[k]
                  + f_0 * mf_139[k];

        t_90[k] = -3.0 * kf_50[k]
                  + f_0 * mf_140[k];

        t_91[k] = -3.0 * kf_51[k]
                  + f_0 * mf_141[k];

        t_92[k] = -3.0 * kf_52[k]
                  + f_0 * mf_142[k];

        t_93[k] = -3.0 * kf_53[k]
                  + f_0 * mf_143[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, kf_54, kf_55, kf_56, kf_57, kf_58, \
                         mf_144, mf_145, mf_146, mf_147, mf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -3.0 * kf_54[k]
                  + f_0 * mf_144[k];

        t_95[k] = -3.0 * kf_55[k]
                  + f_0 * mf_145[k];

        t_96[k] = -3.0 * kf_56[k]
                  + f_0 * mf_146[k];

        t_97[k] = -3.0 * kf_57[k]
                  + f_0 * mf_147[k];

        t_98[k] = -3.0 * kf_58[k]
                  + f_0 * mf_148[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, t_105, kf_59, mf_149, \
                         mf_160, mf_161, mf_162, mf_163, mf_164, \
                         mf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -3.0 * kf_59[k]
                  + f_0 * mf_149[k];

        t_100[k] = f_0 * mf_160[k];

        t_101[k] = f_0 * mf_161[k];

        t_102[k] = f_0 * mf_162[k];

        t_103[k] = f_0 * mf_163[k];

        t_104[k] = f_0 * mf_164[k];

        t_105[k] = f_0 * mf_165[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, kf_60, kf_61, mf_166, \
                         mf_167, mf_168, mf_169, mf_170, mf_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * mf_166[k];

        t_107[k] = f_0 * mf_167[k];

        t_108[k] = f_0 * mf_168[k];

        t_109[k] = f_0 * mf_169[k];

        t_110[k] = -kf_60[k]
                   + f_0 * mf_170[k];

        t_111[k] = -kf_61[k]
                   + f_0 * mf_171[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, kf_62, kf_63, kf_64, kf_65, kf_66, \
                         mf_172, mf_173, mf_174, mf_175, mf_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -kf_62[k]
                   + f_0 * mf_172[k];

        t_113[k] = -kf_63[k]
                   + f_0 * mf_173[k];

        t_114[k] = -kf_64[k]
                   + f_0 * mf_174[k];

        t_115[k] = -kf_65[k]
                   + f_0 * mf_175[k];

        t_116[k] = -kf_66[k]
                   + f_0 * mf_176[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, kf_67, kf_68, kf_69, kf_70, kf_71, \
                         mf_177, mf_178, mf_179, mf_180, mf_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -kf_67[k]
                   + f_0 * mf_177[k];

        t_118[k] = -kf_68[k]
                   + f_0 * mf_178[k];

        t_119[k] = -kf_69[k]
                   + f_0 * mf_179[k];

        t_120[k] = -2.0 * kf_70[k]
                   + f_0 * mf_180[k];

        t_121[k] = -2.0 * kf_71[k]
                   + f_0 * mf_181[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, kf_72, kf_73, kf_74, kf_75, kf_76, \
                         mf_182, mf_183, mf_184, mf_185, mf_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * kf_72[k]
                   + f_0 * mf_182[k];

        t_123[k] = -2.0 * kf_73[k]
                   + f_0 * mf_183[k];

        t_124[k] = -2.0 * kf_74[k]
                   + f_0 * mf_184[k];

        t_125[k] = -2.0 * kf_75[k]
                   + f_0 * mf_185[k];

        t_126[k] = -2.0 * kf_76[k]
                   + f_0 * mf_186[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, kf_77, kf_78, kf_79, kf_80, kf_81, \
                         mf_187, mf_188, mf_189, mf_190, mf_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * kf_77[k]
                   + f_0 * mf_187[k];

        t_128[k] = -2.0 * kf_78[k]
                   + f_0 * mf_188[k];

        t_129[k] = -2.0 * kf_79[k]
                   + f_0 * mf_189[k];

        t_130[k] = -3.0 * kf_80[k]
                   + f_0 * mf_190[k];

        t_131[k] = -3.0 * kf_81[k]
                   + f_0 * mf_191[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, kf_82, kf_83, kf_84, kf_85, kf_86, \
                         mf_192, mf_193, mf_194, mf_195, mf_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -3.0 * kf_82[k]
                   + f_0 * mf_192[k];

        t_133[k] = -3.0 * kf_83[k]
                   + f_0 * mf_193[k];

        t_134[k] = -3.0 * kf_84[k]
                   + f_0 * mf_194[k];

        t_135[k] = -3.0 * kf_85[k]
                   + f_0 * mf_195[k];

        t_136[k] = -3.0 * kf_86[k]
                   + f_0 * mf_196[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, kf_87, kf_88, kf_89, kf_90, kf_91, \
                         mf_197, mf_198, mf_199, mf_200, mf_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -3.0 * kf_87[k]
                   + f_0 * mf_197[k];

        t_138[k] = -3.0 * kf_88[k]
                   + f_0 * mf_198[k];

        t_139[k] = -3.0 * kf_89[k]
                   + f_0 * mf_199[k];

        t_140[k] = -4.0 * kf_90[k]
                   + f_0 * mf_200[k];

        t_141[k] = -4.0 * kf_91[k]
                   + f_0 * mf_201[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, kf_92, kf_93, kf_94, kf_95, kf_96, \
                         mf_202, mf_203, mf_204, mf_205, mf_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -4.0 * kf_92[k]
                   + f_0 * mf_202[k];

        t_143[k] = -4.0 * kf_93[k]
                   + f_0 * mf_203[k];

        t_144[k] = -4.0 * kf_94[k]
                   + f_0 * mf_204[k];

        t_145[k] = -4.0 * kf_95[k]
                   + f_0 * mf_205[k];

        t_146[k] = -4.0 * kf_96[k]
                   + f_0 * mf_206[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, t_152, kf_97, kf_98, kf_99, \
                         mf_207, mf_208, mf_209, mf_220, mf_221, \
                         mf_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -4.0 * kf_97[k]
                   + f_0 * mf_207[k];

        t_148[k] = -4.0 * kf_98[k]
                   + f_0 * mf_208[k];

        t_149[k] = -4.0 * kf_99[k]
                   + f_0 * mf_209[k];

        t_150[k] = f_0 * mf_220[k];

        t_151[k] = f_0 * mf_221[k];

        t_152[k] = f_0 * mf_222[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, t_158, t_159, mf_223, mf_224, \
                         mf_225, mf_226, mf_227, mf_228, mf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_0 * mf_223[k];

        t_154[k] = f_0 * mf_224[k];

        t_155[k] = f_0 * mf_225[k];

        t_156[k] = f_0 * mf_226[k];

        t_157[k] = f_0 * mf_227[k];

        t_158[k] = f_0 * mf_228[k];

        t_159[k] = f_0 * mf_229[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, kf_100, kf_101, kf_102, kf_103, \
                         kf_104, mf_230, mf_231, mf_232, mf_233, \
                         mf_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -kf_100[k]
                   + f_0 * mf_230[k];

        t_161[k] = -kf_101[k]
                   + f_0 * mf_231[k];

        t_162[k] = -kf_102[k]
                   + f_0 * mf_232[k];

        t_163[k] = -kf_103[k]
                   + f_0 * mf_233[k];

        t_164[k] = -kf_104[k]
                   + f_0 * mf_234[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, kf_105, kf_106, kf_107, kf_108, \
                         kf_109, mf_235, mf_236, mf_237, mf_238, \
                         mf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -kf_105[k]
                   + f_0 * mf_235[k];

        t_166[k] = -kf_106[k]
                   + f_0 * mf_236[k];

        t_167[k] = -kf_107[k]
                   + f_0 * mf_237[k];

        t_168[k] = -kf_108[k]
                   + f_0 * mf_238[k];

        t_169[k] = -kf_109[k]
                   + f_0 * mf_239[k];
    }
}

static auto
compute_prim_geom_10_lf_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kf, const size_t mf,
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

    const auto *kf_110 = buffer.data(kf + 110);
    const auto *kf_111 = buffer.data(kf + 111);
    const auto *kf_112 = buffer.data(kf + 112);
    const auto *kf_113 = buffer.data(kf + 113);
    const auto *kf_114 = buffer.data(kf + 114);
    const auto *kf_115 = buffer.data(kf + 115);
    const auto *kf_116 = buffer.data(kf + 116);
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_119 = buffer.data(kf + 119);
    const auto *kf_120 = buffer.data(kf + 120);
    const auto *kf_121 = buffer.data(kf + 121);
    const auto *kf_122 = buffer.data(kf + 122);
    const auto *kf_123 = buffer.data(kf + 123);
    const auto *kf_124 = buffer.data(kf + 124);
    const auto *kf_125 = buffer.data(kf + 125);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_127 = buffer.data(kf + 127);
    const auto *kf_128 = buffer.data(kf + 128);
    const auto *kf_129 = buffer.data(kf + 129);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_131 = buffer.data(kf + 131);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_133 = buffer.data(kf + 133);
    const auto *kf_134 = buffer.data(kf + 134);
    const auto *kf_135 = buffer.data(kf + 135);
    const auto *kf_136 = buffer.data(kf + 136);
    const auto *kf_137 = buffer.data(kf + 137);
    const auto *kf_138 = buffer.data(kf + 138);
    const auto *kf_139 = buffer.data(kf + 139);
    const auto *kf_140 = buffer.data(kf + 140);
    const auto *kf_141 = buffer.data(kf + 141);
    const auto *kf_142 = buffer.data(kf + 142);
    const auto *kf_143 = buffer.data(kf + 143);
    const auto *kf_144 = buffer.data(kf + 144);
    const auto *kf_145 = buffer.data(kf + 145);
    const auto *kf_146 = buffer.data(kf + 146);
    const auto *kf_147 = buffer.data(kf + 147);
    const auto *kf_148 = buffer.data(kf + 148);
    const auto *kf_149 = buffer.data(kf + 149);
    const auto *kf_150 = buffer.data(kf + 150);
    const auto *kf_151 = buffer.data(kf + 151);
    const auto *kf_152 = buffer.data(kf + 152);
    const auto *kf_153 = buffer.data(kf + 153);
    const auto *kf_154 = buffer.data(kf + 154);
    const auto *kf_155 = buffer.data(kf + 155);
    const auto *kf_156 = buffer.data(kf + 156);
    const auto *kf_157 = buffer.data(kf + 157);
    const auto *kf_158 = buffer.data(kf + 158);
    const auto *kf_159 = buffer.data(kf + 159);
    const auto *kf_160 = buffer.data(kf + 160);
    const auto *kf_161 = buffer.data(kf + 161);
    const auto *kf_162 = buffer.data(kf + 162);
    const auto *kf_163 = buffer.data(kf + 163);
    const auto *kf_164 = buffer.data(kf + 164);
    const auto *kf_165 = buffer.data(kf + 165);
    const auto *kf_166 = buffer.data(kf + 166);
    const auto *kf_167 = buffer.data(kf + 167);
    const auto *kf_168 = buffer.data(kf + 168);
    const auto *kf_169 = buffer.data(kf + 169);
    const auto *kf_170 = buffer.data(kf + 170);
    const auto *kf_171 = buffer.data(kf + 171);
    const auto *kf_172 = buffer.data(kf + 172);
    const auto *kf_173 = buffer.data(kf + 173);
    const auto *kf_174 = buffer.data(kf + 174);
    const auto *kf_175 = buffer.data(kf + 175);
    const auto *kf_176 = buffer.data(kf + 176);
    const auto *kf_177 = buffer.data(kf + 177);
    const auto *kf_178 = buffer.data(kf + 178);
    const auto *kf_179 = buffer.data(kf + 179);
    const auto *kf_180 = buffer.data(kf + 180);
    const auto *kf_181 = buffer.data(kf + 181);
    const auto *kf_182 = buffer.data(kf + 182);
    const auto *kf_183 = buffer.data(kf + 183);
    const auto *kf_184 = buffer.data(kf + 184);
    const auto *kf_185 = buffer.data(kf + 185);
    const auto *kf_186 = buffer.data(kf + 186);
    const auto *kf_187 = buffer.data(kf + 187);
    const auto *kf_188 = buffer.data(kf + 188);
    const auto *kf_189 = buffer.data(kf + 189);
    const auto *kf_190 = buffer.data(kf + 190);
    const auto *kf_191 = buffer.data(kf + 191);
    const auto *kf_192 = buffer.data(kf + 192);
    const auto *kf_193 = buffer.data(kf + 193);
    const auto *kf_194 = buffer.data(kf + 194);
    const auto *kf_195 = buffer.data(kf + 195);
    const auto *kf_196 = buffer.data(kf + 196);
    const auto *kf_197 = buffer.data(kf + 197);
    const auto *kf_198 = buffer.data(kf + 198);
    const auto *kf_199 = buffer.data(kf + 199);
    const auto *kf_200 = buffer.data(kf + 200);
    const auto *kf_201 = buffer.data(kf + 201);
    const auto *kf_202 = buffer.data(kf + 202);
    const auto *kf_203 = buffer.data(kf + 203);
    const auto *kf_204 = buffer.data(kf + 204);
    const auto *kf_205 = buffer.data(kf + 205);
    const auto *kf_206 = buffer.data(kf + 206);
    const auto *kf_207 = buffer.data(kf + 207);
    const auto *kf_208 = buffer.data(kf + 208);
    const auto *kf_209 = buffer.data(kf + 209);
    const auto *kf_210 = buffer.data(kf + 210);
    const auto *kf_211 = buffer.data(kf + 211);
    const auto *kf_212 = buffer.data(kf + 212);
    const auto *kf_213 = buffer.data(kf + 213);
    const auto *kf_214 = buffer.data(kf + 214);
    const auto *kf_215 = buffer.data(kf + 215);
    const auto *kf_216 = buffer.data(kf + 216);
    const auto *kf_217 = buffer.data(kf + 217);
    const auto *kf_218 = buffer.data(kf + 218);
    const auto *kf_219 = buffer.data(kf + 219);
    const auto *kf_220 = buffer.data(kf + 220);
    const auto *kf_221 = buffer.data(kf + 221);
    const auto *kf_222 = buffer.data(kf + 222);
    const auto *kf_223 = buffer.data(kf + 223);
    const auto *kf_224 = buffer.data(kf + 224);
    const auto *kf_225 = buffer.data(kf + 225);
    const auto *kf_226 = buffer.data(kf + 226);
    const auto *kf_227 = buffer.data(kf + 227);
    const auto *kf_228 = buffer.data(kf + 228);
    const auto *kf_229 = buffer.data(kf + 229);
    const auto *kf_230 = buffer.data(kf + 230);
    const auto *kf_231 = buffer.data(kf + 231);
    const auto *kf_232 = buffer.data(kf + 232);
    const auto *kf_233 = buffer.data(kf + 233);
    const auto *kf_234 = buffer.data(kf + 234);
    const auto *kf_235 = buffer.data(kf + 235);
    const auto *kf_236 = buffer.data(kf + 236);
    const auto *kf_237 = buffer.data(kf + 237);
    const auto *kf_238 = buffer.data(kf + 238);
    const auto *kf_239 = buffer.data(kf + 239);
    const auto *kf_240 = buffer.data(kf + 240);
    const auto *kf_241 = buffer.data(kf + 241);
    const auto *kf_242 = buffer.data(kf + 242);
    const auto *kf_243 = buffer.data(kf + 243);
    const auto *kf_244 = buffer.data(kf + 244);
    const auto *kf_245 = buffer.data(kf + 245);
    const auto *kf_246 = buffer.data(kf + 246);

    const auto *mf_240 = buffer.data(mf + 240);
    const auto *mf_241 = buffer.data(mf + 241);
    const auto *mf_242 = buffer.data(mf + 242);
    const auto *mf_243 = buffer.data(mf + 243);
    const auto *mf_244 = buffer.data(mf + 244);
    const auto *mf_245 = buffer.data(mf + 245);
    const auto *mf_246 = buffer.data(mf + 246);
    const auto *mf_247 = buffer.data(mf + 247);
    const auto *mf_248 = buffer.data(mf + 248);
    const auto *mf_249 = buffer.data(mf + 249);
    const auto *mf_250 = buffer.data(mf + 250);
    const auto *mf_251 = buffer.data(mf + 251);
    const auto *mf_252 = buffer.data(mf + 252);
    const auto *mf_253 = buffer.data(mf + 253);
    const auto *mf_254 = buffer.data(mf + 254);
    const auto *mf_255 = buffer.data(mf + 255);
    const auto *mf_256 = buffer.data(mf + 256);
    const auto *mf_257 = buffer.data(mf + 257);
    const auto *mf_258 = buffer.data(mf + 258);
    const auto *mf_259 = buffer.data(mf + 259);
    const auto *mf_260 = buffer.data(mf + 260);
    const auto *mf_261 = buffer.data(mf + 261);
    const auto *mf_262 = buffer.data(mf + 262);
    const auto *mf_263 = buffer.data(mf + 263);
    const auto *mf_264 = buffer.data(mf + 264);
    const auto *mf_265 = buffer.data(mf + 265);
    const auto *mf_266 = buffer.data(mf + 266);
    const auto *mf_267 = buffer.data(mf + 267);
    const auto *mf_268 = buffer.data(mf + 268);
    const auto *mf_269 = buffer.data(mf + 269);
    const auto *mf_270 = buffer.data(mf + 270);
    const auto *mf_271 = buffer.data(mf + 271);
    const auto *mf_272 = buffer.data(mf + 272);
    const auto *mf_273 = buffer.data(mf + 273);
    const auto *mf_274 = buffer.data(mf + 274);
    const auto *mf_275 = buffer.data(mf + 275);
    const auto *mf_276 = buffer.data(mf + 276);
    const auto *mf_277 = buffer.data(mf + 277);
    const auto *mf_278 = buffer.data(mf + 278);
    const auto *mf_279 = buffer.data(mf + 279);
    const auto *mf_290 = buffer.data(mf + 290);
    const auto *mf_291 = buffer.data(mf + 291);
    const auto *mf_292 = buffer.data(mf + 292);
    const auto *mf_293 = buffer.data(mf + 293);
    const auto *mf_294 = buffer.data(mf + 294);
    const auto *mf_295 = buffer.data(mf + 295);
    const auto *mf_296 = buffer.data(mf + 296);
    const auto *mf_297 = buffer.data(mf + 297);
    const auto *mf_298 = buffer.data(mf + 298);
    const auto *mf_299 = buffer.data(mf + 299);
    const auto *mf_300 = buffer.data(mf + 300);
    const auto *mf_301 = buffer.data(mf + 301);
    const auto *mf_302 = buffer.data(mf + 302);
    const auto *mf_303 = buffer.data(mf + 303);
    const auto *mf_304 = buffer.data(mf + 304);
    const auto *mf_305 = buffer.data(mf + 305);
    const auto *mf_306 = buffer.data(mf + 306);
    const auto *mf_307 = buffer.data(mf + 307);
    const auto *mf_308 = buffer.data(mf + 308);
    const auto *mf_309 = buffer.data(mf + 309);
    const auto *mf_310 = buffer.data(mf + 310);
    const auto *mf_311 = buffer.data(mf + 311);
    const auto *mf_312 = buffer.data(mf + 312);
    const auto *mf_313 = buffer.data(mf + 313);
    const auto *mf_314 = buffer.data(mf + 314);
    const auto *mf_315 = buffer.data(mf + 315);
    const auto *mf_316 = buffer.data(mf + 316);
    const auto *mf_317 = buffer.data(mf + 317);
    const auto *mf_318 = buffer.data(mf + 318);
    const auto *mf_319 = buffer.data(mf + 319);
    const auto *mf_320 = buffer.data(mf + 320);
    const auto *mf_321 = buffer.data(mf + 321);
    const auto *mf_322 = buffer.data(mf + 322);
    const auto *mf_323 = buffer.data(mf + 323);
    const auto *mf_324 = buffer.data(mf + 324);
    const auto *mf_325 = buffer.data(mf + 325);
    const auto *mf_326 = buffer.data(mf + 326);
    const auto *mf_327 = buffer.data(mf + 327);
    const auto *mf_328 = buffer.data(mf + 328);
    const auto *mf_329 = buffer.data(mf + 329);
    const auto *mf_330 = buffer.data(mf + 330);
    const auto *mf_331 = buffer.data(mf + 331);
    const auto *mf_332 = buffer.data(mf + 332);
    const auto *mf_333 = buffer.data(mf + 333);
    const auto *mf_334 = buffer.data(mf + 334);
    const auto *mf_335 = buffer.data(mf + 335);
    const auto *mf_336 = buffer.data(mf + 336);
    const auto *mf_337 = buffer.data(mf + 337);
    const auto *mf_338 = buffer.data(mf + 338);
    const auto *mf_339 = buffer.data(mf + 339);
    const auto *mf_340 = buffer.data(mf + 340);
    const auto *mf_341 = buffer.data(mf + 341);
    const auto *mf_342 = buffer.data(mf + 342);
    const auto *mf_343 = buffer.data(mf + 343);
    const auto *mf_344 = buffer.data(mf + 344);
    const auto *mf_345 = buffer.data(mf + 345);
    const auto *mf_346 = buffer.data(mf + 346);
    const auto *mf_347 = buffer.data(mf + 347);
    const auto *mf_348 = buffer.data(mf + 348);
    const auto *mf_349 = buffer.data(mf + 349);
    const auto *mf_350 = buffer.data(mf + 350);
    const auto *mf_351 = buffer.data(mf + 351);
    const auto *mf_352 = buffer.data(mf + 352);
    const auto *mf_353 = buffer.data(mf + 353);
    const auto *mf_354 = buffer.data(mf + 354);
    const auto *mf_355 = buffer.data(mf + 355);
    const auto *mf_356 = buffer.data(mf + 356);
    const auto *mf_357 = buffer.data(mf + 357);
    const auto *mf_358 = buffer.data(mf + 358);
    const auto *mf_359 = buffer.data(mf + 359);
    const auto *mf_370 = buffer.data(mf + 370);
    const auto *mf_371 = buffer.data(mf + 371);
    const auto *mf_372 = buffer.data(mf + 372);
    const auto *mf_373 = buffer.data(mf + 373);
    const auto *mf_374 = buffer.data(mf + 374);
    const auto *mf_375 = buffer.data(mf + 375);
    const auto *mf_376 = buffer.data(mf + 376);
    const auto *mf_377 = buffer.data(mf + 377);
    const auto *mf_378 = buffer.data(mf + 378);
    const auto *mf_379 = buffer.data(mf + 379);
    const auto *mf_380 = buffer.data(mf + 380);
    const auto *mf_381 = buffer.data(mf + 381);
    const auto *mf_382 = buffer.data(mf + 382);
    const auto *mf_383 = buffer.data(mf + 383);
    const auto *mf_384 = buffer.data(mf + 384);
    const auto *mf_385 = buffer.data(mf + 385);
    const auto *mf_386 = buffer.data(mf + 386);
    const auto *mf_387 = buffer.data(mf + 387);
    const auto *mf_388 = buffer.data(mf + 388);
    const auto *mf_389 = buffer.data(mf + 389);
    const auto *mf_390 = buffer.data(mf + 390);
    const auto *mf_391 = buffer.data(mf + 391);
    const auto *mf_392 = buffer.data(mf + 392);
    const auto *mf_393 = buffer.data(mf + 393);
    const auto *mf_394 = buffer.data(mf + 394);
    const auto *mf_395 = buffer.data(mf + 395);
    const auto *mf_396 = buffer.data(mf + 396);
    const auto *mf_397 = buffer.data(mf + 397);
    const auto *mf_398 = buffer.data(mf + 398);
    const auto *mf_399 = buffer.data(mf + 399);
    const auto *mf_400 = buffer.data(mf + 400);
    const auto *mf_401 = buffer.data(mf + 401);
    const auto *mf_402 = buffer.data(mf + 402);
    const auto *mf_403 = buffer.data(mf + 403);
    const auto *mf_404 = buffer.data(mf + 404);
    const auto *mf_405 = buffer.data(mf + 405);
    const auto *mf_406 = buffer.data(mf + 406);
    const auto *mf_407 = buffer.data(mf + 407);
    const auto *mf_408 = buffer.data(mf + 408);
    const auto *mf_409 = buffer.data(mf + 409);
    const auto *mf_410 = buffer.data(mf + 410);
    const auto *mf_411 = buffer.data(mf + 411);
    const auto *mf_412 = buffer.data(mf + 412);
    const auto *mf_413 = buffer.data(mf + 413);
    const auto *mf_414 = buffer.data(mf + 414);
    const auto *mf_415 = buffer.data(mf + 415);
    const auto *mf_416 = buffer.data(mf + 416);

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, kf_110, kf_111, kf_112, kf_113, \
                         kf_114, mf_240, mf_241, mf_242, mf_243, \
                         mf_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -2.0 * kf_110[k]
                   + f_0 * mf_240[k];

        t_171[k] = -2.0 * kf_111[k]
                   + f_0 * mf_241[k];

        t_172[k] = -2.0 * kf_112[k]
                   + f_0 * mf_242[k];

        t_173[k] = -2.0 * kf_113[k]
                   + f_0 * mf_243[k];

        t_174[k] = -2.0 * kf_114[k]
                   + f_0 * mf_244[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, kf_115, kf_116, kf_117, kf_118, \
                         kf_119, mf_245, mf_246, mf_247, mf_248, \
                         mf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -2.0 * kf_115[k]
                   + f_0 * mf_245[k];

        t_176[k] = -2.0 * kf_116[k]
                   + f_0 * mf_246[k];

        t_177[k] = -2.0 * kf_117[k]
                   + f_0 * mf_247[k];

        t_178[k] = -2.0 * kf_118[k]
                   + f_0 * mf_248[k];

        t_179[k] = -2.0 * kf_119[k]
                   + f_0 * mf_249[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, kf_120, kf_121, kf_122, kf_123, \
                         kf_124, mf_250, mf_251, mf_252, mf_253, \
                         mf_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * kf_120[k]
                   + f_0 * mf_250[k];

        t_181[k] = -3.0 * kf_121[k]
                   + f_0 * mf_251[k];

        t_182[k] = -3.0 * kf_122[k]
                   + f_0 * mf_252[k];

        t_183[k] = -3.0 * kf_123[k]
                   + f_0 * mf_253[k];

        t_184[k] = -3.0 * kf_124[k]
                   + f_0 * mf_254[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, kf_125, kf_126, kf_127, kf_128, \
                         kf_129, mf_255, mf_256, mf_257, mf_258, \
                         mf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * kf_125[k]
                   + f_0 * mf_255[k];

        t_186[k] = -3.0 * kf_126[k]
                   + f_0 * mf_256[k];

        t_187[k] = -3.0 * kf_127[k]
                   + f_0 * mf_257[k];

        t_188[k] = -3.0 * kf_128[k]
                   + f_0 * mf_258[k];

        t_189[k] = -3.0 * kf_129[k]
                   + f_0 * mf_259[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, kf_130, kf_131, kf_132, kf_133, \
                         kf_134, mf_260, mf_261, mf_262, mf_263, \
                         mf_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -4.0 * kf_130[k]
                   + f_0 * mf_260[k];

        t_191[k] = -4.0 * kf_131[k]
                   + f_0 * mf_261[k];

        t_192[k] = -4.0 * kf_132[k]
                   + f_0 * mf_262[k];

        t_193[k] = -4.0 * kf_133[k]
                   + f_0 * mf_263[k];

        t_194[k] = -4.0 * kf_134[k]
                   + f_0 * mf_264[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, kf_135, kf_136, kf_137, kf_138, \
                         kf_139, mf_265, mf_266, mf_267, mf_268, \
                         mf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -4.0 * kf_135[k]
                   + f_0 * mf_265[k];

        t_196[k] = -4.0 * kf_136[k]
                   + f_0 * mf_266[k];

        t_197[k] = -4.0 * kf_137[k]
                   + f_0 * mf_267[k];

        t_198[k] = -4.0 * kf_138[k]
                   + f_0 * mf_268[k];

        t_199[k] = -4.0 * kf_139[k]
                   + f_0 * mf_269[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, kf_140, kf_141, kf_142, kf_143, \
                         kf_144, mf_270, mf_271, mf_272, mf_273, \
                         mf_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -5.0 * kf_140[k]
                   + f_0 * mf_270[k];

        t_201[k] = -5.0 * kf_141[k]
                   + f_0 * mf_271[k];

        t_202[k] = -5.0 * kf_142[k]
                   + f_0 * mf_272[k];

        t_203[k] = -5.0 * kf_143[k]
                   + f_0 * mf_273[k];

        t_204[k] = -5.0 * kf_144[k]
                   + f_0 * mf_274[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, kf_145, kf_146, kf_147, kf_148, \
                         kf_149, mf_275, mf_276, mf_277, mf_278, \
                         mf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -5.0 * kf_145[k]
                   + f_0 * mf_275[k];

        t_206[k] = -5.0 * kf_146[k]
                   + f_0 * mf_276[k];

        t_207[k] = -5.0 * kf_147[k]
                   + f_0 * mf_277[k];

        t_208[k] = -5.0 * kf_148[k]
                   + f_0 * mf_278[k];

        t_209[k] = -5.0 * kf_149[k]
                   + f_0 * mf_279[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, t_217, mf_290, \
                         mf_291, mf_292, mf_293, mf_294, mf_295, mf_296, \
                         mf_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * mf_290[k];

        t_211[k] = f_0 * mf_291[k];

        t_212[k] = f_0 * mf_292[k];

        t_213[k] = f_0 * mf_293[k];

        t_214[k] = f_0 * mf_294[k];

        t_215[k] = f_0 * mf_295[k];

        t_216[k] = f_0 * mf_296[k];

        t_217[k] = f_0 * mf_297[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, kf_150, kf_151, kf_152, \
                         kf_153, mf_298, mf_299, mf_300, mf_301, mf_302, \
                         mf_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * mf_298[k];

        t_219[k] = f_0 * mf_299[k];

        t_220[k] = -kf_150[k]
                   + f_0 * mf_300[k];

        t_221[k] = -kf_151[k]
                   + f_0 * mf_301[k];

        t_222[k] = -kf_152[k]
                   + f_0 * mf_302[k];

        t_223[k] = -kf_153[k]
                   + f_0 * mf_303[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, kf_154, kf_155, kf_156, kf_157, \
                         kf_158, mf_304, mf_305, mf_306, mf_307, \
                         mf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = -kf_154[k]
                   + f_0 * mf_304[k];

        t_225[k] = -kf_155[k]
                   + f_0 * mf_305[k];

        t_226[k] = -kf_156[k]
                   + f_0 * mf_306[k];

        t_227[k] = -kf_157[k]
                   + f_0 * mf_307[k];

        t_228[k] = -kf_158[k]
                   + f_0 * mf_308[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, kf_159, kf_160, kf_161, kf_162, \
                         kf_163, mf_309, mf_310, mf_311, mf_312, \
                         mf_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = -kf_159[k]
                   + f_0 * mf_309[k];

        t_230[k] = -2.0 * kf_160[k]
                   + f_0 * mf_310[k];

        t_231[k] = -2.0 * kf_161[k]
                   + f_0 * mf_311[k];

        t_232[k] = -2.0 * kf_162[k]
                   + f_0 * mf_312[k];

        t_233[k] = -2.0 * kf_163[k]
                   + f_0 * mf_313[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, kf_164, kf_165, kf_166, kf_167, \
                         kf_168, mf_314, mf_315, mf_316, mf_317, \
                         mf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = -2.0 * kf_164[k]
                   + f_0 * mf_314[k];

        t_235[k] = -2.0 * kf_165[k]
                   + f_0 * mf_315[k];

        t_236[k] = -2.0 * kf_166[k]
                   + f_0 * mf_316[k];

        t_237[k] = -2.0 * kf_167[k]
                   + f_0 * mf_317[k];

        t_238[k] = -2.0 * kf_168[k]
                   + f_0 * mf_318[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, kf_169, kf_170, kf_171, kf_172, \
                         kf_173, mf_319, mf_320, mf_321, mf_322, \
                         mf_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = -2.0 * kf_169[k]
                   + f_0 * mf_319[k];

        t_240[k] = -3.0 * kf_170[k]
                   + f_0 * mf_320[k];

        t_241[k] = -3.0 * kf_171[k]
                   + f_0 * mf_321[k];

        t_242[k] = -3.0 * kf_172[k]
                   + f_0 * mf_322[k];

        t_243[k] = -3.0 * kf_173[k]
                   + f_0 * mf_323[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, t_248, kf_174, kf_175, kf_176, kf_177, \
                         kf_178, mf_324, mf_325, mf_326, mf_327, \
                         mf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = -3.0 * kf_174[k]
                   + f_0 * mf_324[k];

        t_245[k] = -3.0 * kf_175[k]
                   + f_0 * mf_325[k];

        t_246[k] = -3.0 * kf_176[k]
                   + f_0 * mf_326[k];

        t_247[k] = -3.0 * kf_177[k]
                   + f_0 * mf_327[k];

        t_248[k] = -3.0 * kf_178[k]
                   + f_0 * mf_328[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, t_253, kf_179, kf_180, kf_181, kf_182, \
                         kf_183, mf_329, mf_330, mf_331, mf_332, \
                         mf_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = -3.0 * kf_179[k]
                   + f_0 * mf_329[k];

        t_250[k] = -4.0 * kf_180[k]
                   + f_0 * mf_330[k];

        t_251[k] = -4.0 * kf_181[k]
                   + f_0 * mf_331[k];

        t_252[k] = -4.0 * kf_182[k]
                   + f_0 * mf_332[k];

        t_253[k] = -4.0 * kf_183[k]
                   + f_0 * mf_333[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, t_258, kf_184, kf_185, kf_186, kf_187, \
                         kf_188, mf_334, mf_335, mf_336, mf_337, \
                         mf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = -4.0 * kf_184[k]
                   + f_0 * mf_334[k];

        t_255[k] = -4.0 * kf_185[k]
                   + f_0 * mf_335[k];

        t_256[k] = -4.0 * kf_186[k]
                   + f_0 * mf_336[k];

        t_257[k] = -4.0 * kf_187[k]
                   + f_0 * mf_337[k];

        t_258[k] = -4.0 * kf_188[k]
                   + f_0 * mf_338[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, kf_189, kf_190, kf_191, kf_192, \
                         kf_193, mf_339, mf_340, mf_341, mf_342, \
                         mf_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = -4.0 * kf_189[k]
                   + f_0 * mf_339[k];

        t_260[k] = -5.0 * kf_190[k]
                   + f_0 * mf_340[k];

        t_261[k] = -5.0 * kf_191[k]
                   + f_0 * mf_341[k];

        t_262[k] = -5.0 * kf_192[k]
                   + f_0 * mf_342[k];

        t_263[k] = -5.0 * kf_193[k]
                   + f_0 * mf_343[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, kf_194, kf_195, kf_196, kf_197, \
                         kf_198, mf_344, mf_345, mf_346, mf_347, \
                         mf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = -5.0 * kf_194[k]
                   + f_0 * mf_344[k];

        t_265[k] = -5.0 * kf_195[k]
                   + f_0 * mf_345[k];

        t_266[k] = -5.0 * kf_196[k]
                   + f_0 * mf_346[k];

        t_267[k] = -5.0 * kf_197[k]
                   + f_0 * mf_347[k];

        t_268[k] = -5.0 * kf_198[k]
                   + f_0 * mf_348[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, t_273, kf_199, kf_200, kf_201, kf_202, \
                         kf_203, mf_349, mf_350, mf_351, mf_352, \
                         mf_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = -5.0 * kf_199[k]
                   + f_0 * mf_349[k];

        t_270[k] = -6.0 * kf_200[k]
                   + f_0 * mf_350[k];

        t_271[k] = -6.0 * kf_201[k]
                   + f_0 * mf_351[k];

        t_272[k] = -6.0 * kf_202[k]
                   + f_0 * mf_352[k];

        t_273[k] = -6.0 * kf_203[k]
                   + f_0 * mf_353[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, kf_204, kf_205, kf_206, kf_207, \
                         kf_208, mf_354, mf_355, mf_356, mf_357, \
                         mf_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = -6.0 * kf_204[k]
                   + f_0 * mf_354[k];

        t_275[k] = -6.0 * kf_205[k]
                   + f_0 * mf_355[k];

        t_276[k] = -6.0 * kf_206[k]
                   + f_0 * mf_356[k];

        t_277[k] = -6.0 * kf_207[k]
                   + f_0 * mf_357[k];

        t_278[k] = -6.0 * kf_208[k]
                   + f_0 * mf_358[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, t_283, t_284, t_285, kf_209, mf_359, \
                         mf_370, mf_371, mf_372, mf_373, mf_374, \
                         mf_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = -6.0 * kf_209[k]
                   + f_0 * mf_359[k];

        t_280[k] = f_0 * mf_370[k];

        t_281[k] = f_0 * mf_371[k];

        t_282[k] = f_0 * mf_372[k];

        t_283[k] = f_0 * mf_373[k];

        t_284[k] = f_0 * mf_374[k];

        t_285[k] = f_0 * mf_375[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, t_291, kf_210, kf_211, mf_376, \
                         mf_377, mf_378, mf_379, mf_380, mf_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * mf_376[k];

        t_287[k] = f_0 * mf_377[k];

        t_288[k] = f_0 * mf_378[k];

        t_289[k] = f_0 * mf_379[k];

        t_290[k] = -kf_210[k]
                   + f_0 * mf_380[k];

        t_291[k] = -kf_211[k]
                   + f_0 * mf_381[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, kf_212, kf_213, kf_214, kf_215, \
                         kf_216, mf_382, mf_383, mf_384, mf_385, \
                         mf_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = -kf_212[k]
                   + f_0 * mf_382[k];

        t_293[k] = -kf_213[k]
                   + f_0 * mf_383[k];

        t_294[k] = -kf_214[k]
                   + f_0 * mf_384[k];

        t_295[k] = -kf_215[k]
                   + f_0 * mf_385[k];

        t_296[k] = -kf_216[k]
                   + f_0 * mf_386[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, t_301, kf_217, kf_218, kf_219, kf_220, \
                         kf_221, mf_387, mf_388, mf_389, mf_390, \
                         mf_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = -kf_217[k]
                   + f_0 * mf_387[k];

        t_298[k] = -kf_218[k]
                   + f_0 * mf_388[k];

        t_299[k] = -kf_219[k]
                   + f_0 * mf_389[k];

        t_300[k] = -2.0 * kf_220[k]
                   + f_0 * mf_390[k];

        t_301[k] = -2.0 * kf_221[k]
                   + f_0 * mf_391[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, kf_222, kf_223, kf_224, kf_225, \
                         kf_226, mf_392, mf_393, mf_394, mf_395, \
                         mf_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = -2.0 * kf_222[k]
                   + f_0 * mf_392[k];

        t_303[k] = -2.0 * kf_223[k]
                   + f_0 * mf_393[k];

        t_304[k] = -2.0 * kf_224[k]
                   + f_0 * mf_394[k];

        t_305[k] = -2.0 * kf_225[k]
                   + f_0 * mf_395[k];

        t_306[k] = -2.0 * kf_226[k]
                   + f_0 * mf_396[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, kf_227, kf_228, kf_229, kf_230, \
                         kf_231, mf_397, mf_398, mf_399, mf_400, \
                         mf_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = -2.0 * kf_227[k]
                   + f_0 * mf_397[k];

        t_308[k] = -2.0 * kf_228[k]
                   + f_0 * mf_398[k];

        t_309[k] = -2.0 * kf_229[k]
                   + f_0 * mf_399[k];

        t_310[k] = -3.0 * kf_230[k]
                   + f_0 * mf_400[k];

        t_311[k] = -3.0 * kf_231[k]
                   + f_0 * mf_401[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, kf_232, kf_233, kf_234, kf_235, \
                         kf_236, mf_402, mf_403, mf_404, mf_405, \
                         mf_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -3.0 * kf_232[k]
                   + f_0 * mf_402[k];

        t_313[k] = -3.0 * kf_233[k]
                   + f_0 * mf_403[k];

        t_314[k] = -3.0 * kf_234[k]
                   + f_0 * mf_404[k];

        t_315[k] = -3.0 * kf_235[k]
                   + f_0 * mf_405[k];

        t_316[k] = -3.0 * kf_236[k]
                   + f_0 * mf_406[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, kf_237, kf_238, kf_239, kf_240, \
                         kf_241, mf_407, mf_408, mf_409, mf_410, \
                         mf_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = -3.0 * kf_237[k]
                   + f_0 * mf_407[k];

        t_318[k] = -3.0 * kf_238[k]
                   + f_0 * mf_408[k];

        t_319[k] = -3.0 * kf_239[k]
                   + f_0 * mf_409[k];

        t_320[k] = -4.0 * kf_240[k]
                   + f_0 * mf_410[k];

        t_321[k] = -4.0 * kf_241[k]
                   + f_0 * mf_411[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, kf_242, kf_243, kf_244, kf_245, \
                         kf_246, mf_412, mf_413, mf_414, mf_415, \
                         mf_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = -4.0 * kf_242[k]
                   + f_0 * mf_412[k];

        t_323[k] = -4.0 * kf_243[k]
                   + f_0 * mf_413[k];

        t_324[k] = -4.0 * kf_244[k]
                   + f_0 * mf_414[k];

        t_325[k] = -4.0 * kf_245[k]
                   + f_0 * mf_415[k];

        t_326[k] = -4.0 * kf_246[k]
                   + f_0 * mf_416[k];
    }
}

static auto
compute_prim_geom_10_lf_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kf, const size_t mf,
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

    const auto *kf_247 = buffer.data(kf + 247);
    const auto *kf_248 = buffer.data(kf + 248);
    const auto *kf_249 = buffer.data(kf + 249);
    const auto *kf_250 = buffer.data(kf + 250);
    const auto *kf_251 = buffer.data(kf + 251);
    const auto *kf_252 = buffer.data(kf + 252);
    const auto *kf_253 = buffer.data(kf + 253);
    const auto *kf_254 = buffer.data(kf + 254);
    const auto *kf_255 = buffer.data(kf + 255);
    const auto *kf_256 = buffer.data(kf + 256);
    const auto *kf_257 = buffer.data(kf + 257);
    const auto *kf_258 = buffer.data(kf + 258);
    const auto *kf_259 = buffer.data(kf + 259);
    const auto *kf_260 = buffer.data(kf + 260);
    const auto *kf_261 = buffer.data(kf + 261);
    const auto *kf_262 = buffer.data(kf + 262);
    const auto *kf_263 = buffer.data(kf + 263);
    const auto *kf_264 = buffer.data(kf + 264);
    const auto *kf_265 = buffer.data(kf + 265);
    const auto *kf_266 = buffer.data(kf + 266);
    const auto *kf_267 = buffer.data(kf + 267);
    const auto *kf_268 = buffer.data(kf + 268);
    const auto *kf_269 = buffer.data(kf + 269);
    const auto *kf_270 = buffer.data(kf + 270);
    const auto *kf_271 = buffer.data(kf + 271);
    const auto *kf_272 = buffer.data(kf + 272);
    const auto *kf_273 = buffer.data(kf + 273);
    const auto *kf_274 = buffer.data(kf + 274);
    const auto *kf_275 = buffer.data(kf + 275);
    const auto *kf_276 = buffer.data(kf + 276);
    const auto *kf_277 = buffer.data(kf + 277);
    const auto *kf_278 = buffer.data(kf + 278);
    const auto *kf_279 = buffer.data(kf + 279);
    const auto *kf_280 = buffer.data(kf + 280);
    const auto *kf_281 = buffer.data(kf + 281);
    const auto *kf_282 = buffer.data(kf + 282);
    const auto *kf_283 = buffer.data(kf + 283);
    const auto *kf_284 = buffer.data(kf + 284);
    const auto *kf_285 = buffer.data(kf + 285);
    const auto *kf_286 = buffer.data(kf + 286);
    const auto *kf_287 = buffer.data(kf + 287);
    const auto *kf_288 = buffer.data(kf + 288);
    const auto *kf_289 = buffer.data(kf + 289);
    const auto *kf_290 = buffer.data(kf + 290);
    const auto *kf_291 = buffer.data(kf + 291);
    const auto *kf_292 = buffer.data(kf + 292);
    const auto *kf_293 = buffer.data(kf + 293);
    const auto *kf_294 = buffer.data(kf + 294);
    const auto *kf_295 = buffer.data(kf + 295);
    const auto *kf_296 = buffer.data(kf + 296);
    const auto *kf_297 = buffer.data(kf + 297);
    const auto *kf_298 = buffer.data(kf + 298);
    const auto *kf_299 = buffer.data(kf + 299);
    const auto *kf_300 = buffer.data(kf + 300);
    const auto *kf_301 = buffer.data(kf + 301);
    const auto *kf_302 = buffer.data(kf + 302);
    const auto *kf_303 = buffer.data(kf + 303);
    const auto *kf_304 = buffer.data(kf + 304);
    const auto *kf_305 = buffer.data(kf + 305);
    const auto *kf_306 = buffer.data(kf + 306);
    const auto *kf_307 = buffer.data(kf + 307);
    const auto *kf_308 = buffer.data(kf + 308);
    const auto *kf_309 = buffer.data(kf + 309);
    const auto *kf_310 = buffer.data(kf + 310);
    const auto *kf_311 = buffer.data(kf + 311);
    const auto *kf_312 = buffer.data(kf + 312);
    const auto *kf_313 = buffer.data(kf + 313);
    const auto *kf_314 = buffer.data(kf + 314);
    const auto *kf_315 = buffer.data(kf + 315);
    const auto *kf_316 = buffer.data(kf + 316);
    const auto *kf_317 = buffer.data(kf + 317);
    const auto *kf_318 = buffer.data(kf + 318);
    const auto *kf_319 = buffer.data(kf + 319);
    const auto *kf_320 = buffer.data(kf + 320);
    const auto *kf_321 = buffer.data(kf + 321);
    const auto *kf_322 = buffer.data(kf + 322);
    const auto *kf_323 = buffer.data(kf + 323);
    const auto *kf_324 = buffer.data(kf + 324);
    const auto *kf_325 = buffer.data(kf + 325);
    const auto *kf_326 = buffer.data(kf + 326);
    const auto *kf_327 = buffer.data(kf + 327);
    const auto *kf_328 = buffer.data(kf + 328);
    const auto *kf_329 = buffer.data(kf + 329);
    const auto *kf_330 = buffer.data(kf + 330);
    const auto *kf_331 = buffer.data(kf + 331);
    const auto *kf_332 = buffer.data(kf + 332);
    const auto *kf_333 = buffer.data(kf + 333);
    const auto *kf_334 = buffer.data(kf + 334);
    const auto *kf_335 = buffer.data(kf + 335);
    const auto *kf_336 = buffer.data(kf + 336);
    const auto *kf_337 = buffer.data(kf + 337);
    const auto *kf_338 = buffer.data(kf + 338);
    const auto *kf_339 = buffer.data(kf + 339);
    const auto *kf_340 = buffer.data(kf + 340);
    const auto *kf_341 = buffer.data(kf + 341);
    const auto *kf_342 = buffer.data(kf + 342);
    const auto *kf_343 = buffer.data(kf + 343);
    const auto *kf_344 = buffer.data(kf + 344);
    const auto *kf_345 = buffer.data(kf + 345);
    const auto *kf_346 = buffer.data(kf + 346);
    const auto *kf_347 = buffer.data(kf + 347);
    const auto *kf_348 = buffer.data(kf + 348);
    const auto *kf_349 = buffer.data(kf + 349);
    const auto *kf_350 = buffer.data(kf + 350);
    const auto *kf_351 = buffer.data(kf + 351);
    const auto *kf_352 = buffer.data(kf + 352);
    const auto *kf_353 = buffer.data(kf + 353);
    const auto *kf_354 = buffer.data(kf + 354);
    const auto *kf_355 = buffer.data(kf + 355);
    const auto *kf_356 = buffer.data(kf + 356);
    const auto *kf_357 = buffer.data(kf + 357);
    const auto *kf_358 = buffer.data(kf + 358);
    const auto *kf_359 = buffer.data(kf + 359);

    const auto *mf_417 = buffer.data(mf + 417);
    const auto *mf_418 = buffer.data(mf + 418);
    const auto *mf_419 = buffer.data(mf + 419);
    const auto *mf_420 = buffer.data(mf + 420);
    const auto *mf_421 = buffer.data(mf + 421);
    const auto *mf_422 = buffer.data(mf + 422);
    const auto *mf_423 = buffer.data(mf + 423);
    const auto *mf_424 = buffer.data(mf + 424);
    const auto *mf_425 = buffer.data(mf + 425);
    const auto *mf_426 = buffer.data(mf + 426);
    const auto *mf_427 = buffer.data(mf + 427);
    const auto *mf_428 = buffer.data(mf + 428);
    const auto *mf_429 = buffer.data(mf + 429);
    const auto *mf_430 = buffer.data(mf + 430);
    const auto *mf_431 = buffer.data(mf + 431);
    const auto *mf_432 = buffer.data(mf + 432);
    const auto *mf_433 = buffer.data(mf + 433);
    const auto *mf_434 = buffer.data(mf + 434);
    const auto *mf_435 = buffer.data(mf + 435);
    const auto *mf_436 = buffer.data(mf + 436);
    const auto *mf_437 = buffer.data(mf + 437);
    const auto *mf_438 = buffer.data(mf + 438);
    const auto *mf_439 = buffer.data(mf + 439);
    const auto *mf_440 = buffer.data(mf + 440);
    const auto *mf_441 = buffer.data(mf + 441);
    const auto *mf_442 = buffer.data(mf + 442);
    const auto *mf_443 = buffer.data(mf + 443);
    const auto *mf_444 = buffer.data(mf + 444);
    const auto *mf_445 = buffer.data(mf + 445);
    const auto *mf_446 = buffer.data(mf + 446);
    const auto *mf_447 = buffer.data(mf + 447);
    const auto *mf_448 = buffer.data(mf + 448);
    const auto *mf_449 = buffer.data(mf + 449);
    const auto *mf_460 = buffer.data(mf + 460);
    const auto *mf_461 = buffer.data(mf + 461);
    const auto *mf_462 = buffer.data(mf + 462);
    const auto *mf_463 = buffer.data(mf + 463);
    const auto *mf_464 = buffer.data(mf + 464);
    const auto *mf_465 = buffer.data(mf + 465);
    const auto *mf_466 = buffer.data(mf + 466);
    const auto *mf_467 = buffer.data(mf + 467);
    const auto *mf_468 = buffer.data(mf + 468);
    const auto *mf_469 = buffer.data(mf + 469);
    const auto *mf_470 = buffer.data(mf + 470);
    const auto *mf_471 = buffer.data(mf + 471);
    const auto *mf_472 = buffer.data(mf + 472);
    const auto *mf_473 = buffer.data(mf + 473);
    const auto *mf_474 = buffer.data(mf + 474);
    const auto *mf_475 = buffer.data(mf + 475);
    const auto *mf_476 = buffer.data(mf + 476);
    const auto *mf_477 = buffer.data(mf + 477);
    const auto *mf_478 = buffer.data(mf + 478);
    const auto *mf_479 = buffer.data(mf + 479);
    const auto *mf_480 = buffer.data(mf + 480);
    const auto *mf_481 = buffer.data(mf + 481);
    const auto *mf_482 = buffer.data(mf + 482);
    const auto *mf_483 = buffer.data(mf + 483);
    const auto *mf_484 = buffer.data(mf + 484);
    const auto *mf_485 = buffer.data(mf + 485);
    const auto *mf_486 = buffer.data(mf + 486);
    const auto *mf_487 = buffer.data(mf + 487);
    const auto *mf_488 = buffer.data(mf + 488);
    const auto *mf_489 = buffer.data(mf + 489);
    const auto *mf_490 = buffer.data(mf + 490);
    const auto *mf_491 = buffer.data(mf + 491);
    const auto *mf_492 = buffer.data(mf + 492);
    const auto *mf_493 = buffer.data(mf + 493);
    const auto *mf_494 = buffer.data(mf + 494);
    const auto *mf_495 = buffer.data(mf + 495);
    const auto *mf_496 = buffer.data(mf + 496);
    const auto *mf_497 = buffer.data(mf + 497);
    const auto *mf_498 = buffer.data(mf + 498);
    const auto *mf_499 = buffer.data(mf + 499);
    const auto *mf_500 = buffer.data(mf + 500);
    const auto *mf_501 = buffer.data(mf + 501);
    const auto *mf_502 = buffer.data(mf + 502);
    const auto *mf_503 = buffer.data(mf + 503);
    const auto *mf_504 = buffer.data(mf + 504);
    const auto *mf_505 = buffer.data(mf + 505);
    const auto *mf_506 = buffer.data(mf + 506);
    const auto *mf_507 = buffer.data(mf + 507);
    const auto *mf_508 = buffer.data(mf + 508);
    const auto *mf_509 = buffer.data(mf + 509);
    const auto *mf_510 = buffer.data(mf + 510);
    const auto *mf_511 = buffer.data(mf + 511);
    const auto *mf_512 = buffer.data(mf + 512);
    const auto *mf_513 = buffer.data(mf + 513);
    const auto *mf_514 = buffer.data(mf + 514);
    const auto *mf_515 = buffer.data(mf + 515);
    const auto *mf_516 = buffer.data(mf + 516);
    const auto *mf_517 = buffer.data(mf + 517);
    const auto *mf_518 = buffer.data(mf + 518);
    const auto *mf_519 = buffer.data(mf + 519);
    const auto *mf_520 = buffer.data(mf + 520);
    const auto *mf_521 = buffer.data(mf + 521);
    const auto *mf_522 = buffer.data(mf + 522);
    const auto *mf_523 = buffer.data(mf + 523);
    const auto *mf_524 = buffer.data(mf + 524);
    const auto *mf_525 = buffer.data(mf + 525);
    const auto *mf_526 = buffer.data(mf + 526);
    const auto *mf_527 = buffer.data(mf + 527);
    const auto *mf_528 = buffer.data(mf + 528);
    const auto *mf_529 = buffer.data(mf + 529);
    const auto *mf_530 = buffer.data(mf + 530);
    const auto *mf_531 = buffer.data(mf + 531);
    const auto *mf_532 = buffer.data(mf + 532);
    const auto *mf_533 = buffer.data(mf + 533);
    const auto *mf_534 = buffer.data(mf + 534);
    const auto *mf_535 = buffer.data(mf + 535);
    const auto *mf_536 = buffer.data(mf + 536);
    const auto *mf_537 = buffer.data(mf + 537);
    const auto *mf_538 = buffer.data(mf + 538);
    const auto *mf_539 = buffer.data(mf + 539);
    const auto *mf_540 = buffer.data(mf + 540);
    const auto *mf_541 = buffer.data(mf + 541);
    const auto *mf_542 = buffer.data(mf + 542);
    const auto *mf_543 = buffer.data(mf + 543);
    const auto *mf_544 = buffer.data(mf + 544);
    const auto *mf_545 = buffer.data(mf + 545);
    const auto *mf_546 = buffer.data(mf + 546);
    const auto *mf_547 = buffer.data(mf + 547);
    const auto *mf_548 = buffer.data(mf + 548);
    const auto *mf_549 = buffer.data(mf + 549);

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, kf_247, kf_248, kf_249, kf_250, \
                         kf_251, mf_417, mf_418, mf_419, mf_420, \
                         mf_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = -4.0 * kf_247[k]
                   + f_0 * mf_417[k];

        t_328[k] = -4.0 * kf_248[k]
                   + f_0 * mf_418[k];

        t_329[k] = -4.0 * kf_249[k]
                   + f_0 * mf_419[k];

        t_330[k] = -5.0 * kf_250[k]
                   + f_0 * mf_420[k];

        t_331[k] = -5.0 * kf_251[k]
                   + f_0 * mf_421[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, kf_252, kf_253, kf_254, kf_255, \
                         kf_256, mf_422, mf_423, mf_424, mf_425, \
                         mf_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -5.0 * kf_252[k]
                   + f_0 * mf_422[k];

        t_333[k] = -5.0 * kf_253[k]
                   + f_0 * mf_423[k];

        t_334[k] = -5.0 * kf_254[k]
                   + f_0 * mf_424[k];

        t_335[k] = -5.0 * kf_255[k]
                   + f_0 * mf_425[k];

        t_336[k] = -5.0 * kf_256[k]
                   + f_0 * mf_426[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, kf_257, kf_258, kf_259, kf_260, \
                         kf_261, mf_427, mf_428, mf_429, mf_430, \
                         mf_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -5.0 * kf_257[k]
                   + f_0 * mf_427[k];

        t_338[k] = -5.0 * kf_258[k]
                   + f_0 * mf_428[k];

        t_339[k] = -5.0 * kf_259[k]
                   + f_0 * mf_429[k];

        t_340[k] = -6.0 * kf_260[k]
                   + f_0 * mf_430[k];

        t_341[k] = -6.0 * kf_261[k]
                   + f_0 * mf_431[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, kf_262, kf_263, kf_264, kf_265, \
                         kf_266, mf_432, mf_433, mf_434, mf_435, \
                         mf_436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -6.0 * kf_262[k]
                   + f_0 * mf_432[k];

        t_343[k] = -6.0 * kf_263[k]
                   + f_0 * mf_433[k];

        t_344[k] = -6.0 * kf_264[k]
                   + f_0 * mf_434[k];

        t_345[k] = -6.0 * kf_265[k]
                   + f_0 * mf_435[k];

        t_346[k] = -6.0 * kf_266[k]
                   + f_0 * mf_436[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, kf_267, kf_268, kf_269, kf_270, \
                         kf_271, mf_437, mf_438, mf_439, mf_440, \
                         mf_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -6.0 * kf_267[k]
                   + f_0 * mf_437[k];

        t_348[k] = -6.0 * kf_268[k]
                   + f_0 * mf_438[k];

        t_349[k] = -6.0 * kf_269[k]
                   + f_0 * mf_439[k];

        t_350[k] = -7.0 * kf_270[k]
                   + f_0 * mf_440[k];

        t_351[k] = -7.0 * kf_271[k]
                   + f_0 * mf_441[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, kf_272, kf_273, kf_274, kf_275, \
                         kf_276, mf_442, mf_443, mf_444, mf_445, \
                         mf_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -7.0 * kf_272[k]
                   + f_0 * mf_442[k];

        t_353[k] = -7.0 * kf_273[k]
                   + f_0 * mf_443[k];

        t_354[k] = -7.0 * kf_274[k]
                   + f_0 * mf_444[k];

        t_355[k] = -7.0 * kf_275[k]
                   + f_0 * mf_445[k];

        t_356[k] = -7.0 * kf_276[k]
                   + f_0 * mf_446[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, t_362, kf_277, kf_278, kf_279, \
                         mf_447, mf_448, mf_449, mf_460, mf_461, \
                         mf_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -7.0 * kf_277[k]
                   + f_0 * mf_447[k];

        t_358[k] = -7.0 * kf_278[k]
                   + f_0 * mf_448[k];

        t_359[k] = -7.0 * kf_279[k]
                   + f_0 * mf_449[k];

        t_360[k] = f_0 * mf_460[k];

        t_361[k] = f_0 * mf_461[k];

        t_362[k] = f_0 * mf_462[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, t_368, t_369, mf_463, mf_464, \
                         mf_465, mf_466, mf_467, mf_468, mf_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_0 * mf_463[k];

        t_364[k] = f_0 * mf_464[k];

        t_365[k] = f_0 * mf_465[k];

        t_366[k] = f_0 * mf_466[k];

        t_367[k] = f_0 * mf_467[k];

        t_368[k] = f_0 * mf_468[k];

        t_369[k] = f_0 * mf_469[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, kf_280, kf_281, kf_282, kf_283, \
                         kf_284, mf_470, mf_471, mf_472, mf_473, \
                         mf_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -kf_280[k]
                   + f_0 * mf_470[k];

        t_371[k] = -kf_281[k]
                   + f_0 * mf_471[k];

        t_372[k] = -kf_282[k]
                   + f_0 * mf_472[k];

        t_373[k] = -kf_283[k]
                   + f_0 * mf_473[k];

        t_374[k] = -kf_284[k]
                   + f_0 * mf_474[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, kf_285, kf_286, kf_287, kf_288, \
                         kf_289, mf_475, mf_476, mf_477, mf_478, \
                         mf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -kf_285[k]
                   + f_0 * mf_475[k];

        t_376[k] = -kf_286[k]
                   + f_0 * mf_476[k];

        t_377[k] = -kf_287[k]
                   + f_0 * mf_477[k];

        t_378[k] = -kf_288[k]
                   + f_0 * mf_478[k];

        t_379[k] = -kf_289[k]
                   + f_0 * mf_479[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, kf_290, kf_291, kf_292, kf_293, \
                         kf_294, mf_480, mf_481, mf_482, mf_483, \
                         mf_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * kf_290[k]
                   + f_0 * mf_480[k];

        t_381[k] = -2.0 * kf_291[k]
                   + f_0 * mf_481[k];

        t_382[k] = -2.0 * kf_292[k]
                   + f_0 * mf_482[k];

        t_383[k] = -2.0 * kf_293[k]
                   + f_0 * mf_483[k];

        t_384[k] = -2.0 * kf_294[k]
                   + f_0 * mf_484[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, kf_295, kf_296, kf_297, kf_298, \
                         kf_299, mf_485, mf_486, mf_487, mf_488, \
                         mf_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * kf_295[k]
                   + f_0 * mf_485[k];

        t_386[k] = -2.0 * kf_296[k]
                   + f_0 * mf_486[k];

        t_387[k] = -2.0 * kf_297[k]
                   + f_0 * mf_487[k];

        t_388[k] = -2.0 * kf_298[k]
                   + f_0 * mf_488[k];

        t_389[k] = -2.0 * kf_299[k]
                   + f_0 * mf_489[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, kf_300, kf_301, kf_302, kf_303, \
                         kf_304, mf_490, mf_491, mf_492, mf_493, \
                         mf_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -3.0 * kf_300[k]
                   + f_0 * mf_490[k];

        t_391[k] = -3.0 * kf_301[k]
                   + f_0 * mf_491[k];

        t_392[k] = -3.0 * kf_302[k]
                   + f_0 * mf_492[k];

        t_393[k] = -3.0 * kf_303[k]
                   + f_0 * mf_493[k];

        t_394[k] = -3.0 * kf_304[k]
                   + f_0 * mf_494[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, kf_305, kf_306, kf_307, kf_308, \
                         kf_309, mf_495, mf_496, mf_497, mf_498, \
                         mf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -3.0 * kf_305[k]
                   + f_0 * mf_495[k];

        t_396[k] = -3.0 * kf_306[k]
                   + f_0 * mf_496[k];

        t_397[k] = -3.0 * kf_307[k]
                   + f_0 * mf_497[k];

        t_398[k] = -3.0 * kf_308[k]
                   + f_0 * mf_498[k];

        t_399[k] = -3.0 * kf_309[k]
                   + f_0 * mf_499[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, kf_310, kf_311, kf_312, kf_313, \
                         kf_314, mf_500, mf_501, mf_502, mf_503, \
                         mf_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -4.0 * kf_310[k]
                   + f_0 * mf_500[k];

        t_401[k] = -4.0 * kf_311[k]
                   + f_0 * mf_501[k];

        t_402[k] = -4.0 * kf_312[k]
                   + f_0 * mf_502[k];

        t_403[k] = -4.0 * kf_313[k]
                   + f_0 * mf_503[k];

        t_404[k] = -4.0 * kf_314[k]
                   + f_0 * mf_504[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, kf_315, kf_316, kf_317, kf_318, \
                         kf_319, mf_505, mf_506, mf_507, mf_508, \
                         mf_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -4.0 * kf_315[k]
                   + f_0 * mf_505[k];

        t_406[k] = -4.0 * kf_316[k]
                   + f_0 * mf_506[k];

        t_407[k] = -4.0 * kf_317[k]
                   + f_0 * mf_507[k];

        t_408[k] = -4.0 * kf_318[k]
                   + f_0 * mf_508[k];

        t_409[k] = -4.0 * kf_319[k]
                   + f_0 * mf_509[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, kf_320, kf_321, kf_322, kf_323, \
                         kf_324, mf_510, mf_511, mf_512, mf_513, \
                         mf_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -5.0 * kf_320[k]
                   + f_0 * mf_510[k];

        t_411[k] = -5.0 * kf_321[k]
                   + f_0 * mf_511[k];

        t_412[k] = -5.0 * kf_322[k]
                   + f_0 * mf_512[k];

        t_413[k] = -5.0 * kf_323[k]
                   + f_0 * mf_513[k];

        t_414[k] = -5.0 * kf_324[k]
                   + f_0 * mf_514[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, kf_325, kf_326, kf_327, kf_328, \
                         kf_329, mf_515, mf_516, mf_517, mf_518, \
                         mf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -5.0 * kf_325[k]
                   + f_0 * mf_515[k];

        t_416[k] = -5.0 * kf_326[k]
                   + f_0 * mf_516[k];

        t_417[k] = -5.0 * kf_327[k]
                   + f_0 * mf_517[k];

        t_418[k] = -5.0 * kf_328[k]
                   + f_0 * mf_518[k];

        t_419[k] = -5.0 * kf_329[k]
                   + f_0 * mf_519[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, kf_330, kf_331, kf_332, kf_333, \
                         kf_334, mf_520, mf_521, mf_522, mf_523, \
                         mf_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -6.0 * kf_330[k]
                   + f_0 * mf_520[k];

        t_421[k] = -6.0 * kf_331[k]
                   + f_0 * mf_521[k];

        t_422[k] = -6.0 * kf_332[k]
                   + f_0 * mf_522[k];

        t_423[k] = -6.0 * kf_333[k]
                   + f_0 * mf_523[k];

        t_424[k] = -6.0 * kf_334[k]
                   + f_0 * mf_524[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, kf_335, kf_336, kf_337, kf_338, \
                         kf_339, mf_525, mf_526, mf_527, mf_528, \
                         mf_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -6.0 * kf_335[k]
                   + f_0 * mf_525[k];

        t_426[k] = -6.0 * kf_336[k]
                   + f_0 * mf_526[k];

        t_427[k] = -6.0 * kf_337[k]
                   + f_0 * mf_527[k];

        t_428[k] = -6.0 * kf_338[k]
                   + f_0 * mf_528[k];

        t_429[k] = -6.0 * kf_339[k]
                   + f_0 * mf_529[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, kf_340, kf_341, kf_342, kf_343, \
                         kf_344, mf_530, mf_531, mf_532, mf_533, \
                         mf_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -7.0 * kf_340[k]
                   + f_0 * mf_530[k];

        t_431[k] = -7.0 * kf_341[k]
                   + f_0 * mf_531[k];

        t_432[k] = -7.0 * kf_342[k]
                   + f_0 * mf_532[k];

        t_433[k] = -7.0 * kf_343[k]
                   + f_0 * mf_533[k];

        t_434[k] = -7.0 * kf_344[k]
                   + f_0 * mf_534[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, kf_345, kf_346, kf_347, kf_348, \
                         kf_349, mf_535, mf_536, mf_537, mf_538, \
                         mf_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -7.0 * kf_345[k]
                   + f_0 * mf_535[k];

        t_436[k] = -7.0 * kf_346[k]
                   + f_0 * mf_536[k];

        t_437[k] = -7.0 * kf_347[k]
                   + f_0 * mf_537[k];

        t_438[k] = -7.0 * kf_348[k]
                   + f_0 * mf_538[k];

        t_439[k] = -7.0 * kf_349[k]
                   + f_0 * mf_539[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, kf_350, kf_351, kf_352, kf_353, \
                         kf_354, mf_540, mf_541, mf_542, mf_543, \
                         mf_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -8.0 * kf_350[k]
                   + f_0 * mf_540[k];

        t_441[k] = -8.0 * kf_351[k]
                   + f_0 * mf_541[k];

        t_442[k] = -8.0 * kf_352[k]
                   + f_0 * mf_542[k];

        t_443[k] = -8.0 * kf_353[k]
                   + f_0 * mf_543[k];

        t_444[k] = -8.0 * kf_354[k]
                   + f_0 * mf_544[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, kf_355, kf_356, kf_357, kf_358, \
                         kf_359, mf_545, mf_546, mf_547, mf_548, \
                         mf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -8.0 * kf_355[k]
                   + f_0 * mf_545[k];

        t_446[k] = -8.0 * kf_356[k]
                   + f_0 * mf_546[k];

        t_447[k] = -8.0 * kf_357[k]
                   + f_0 * mf_547[k];

        t_448[k] = -8.0 * kf_358[k]
                   + f_0 * mf_548[k];

        t_449[k] = -8.0 * kf_359[k]
                   + f_0 * mf_549[k];
    }
}

auto
compute_prim_geom_10_lf_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t kf, const size_t mf,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_lf_electron_repulsion_2_piece0(buffer, target, kf, mf, ncols, alpha);

    compute_prim_geom_10_lf_electron_repulsion_2_piece1(buffer, target, kf, mf, ncols, alpha);

    compute_prim_geom_10_lf_electron_repulsion_2_piece2(buffer, target, kf, mf, ncols, alpha);
}

}  // namespace simdt2ceri
