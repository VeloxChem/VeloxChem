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


#include "SimdElectronRepulsionGeom10VrrRecIL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_il_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
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

    const auto *hl_0 = buffer.data(hl + 0);
    const auto *hl_1 = buffer.data(hl + 1);
    const auto *hl_2 = buffer.data(hl + 2);
    const auto *hl_3 = buffer.data(hl + 3);
    const auto *hl_4 = buffer.data(hl + 4);
    const auto *hl_5 = buffer.data(hl + 5);
    const auto *hl_6 = buffer.data(hl + 6);
    const auto *hl_7 = buffer.data(hl + 7);
    const auto *hl_8 = buffer.data(hl + 8);
    const auto *hl_9 = buffer.data(hl + 9);
    const auto *hl_10 = buffer.data(hl + 10);
    const auto *hl_11 = buffer.data(hl + 11);
    const auto *hl_12 = buffer.data(hl + 12);
    const auto *hl_13 = buffer.data(hl + 13);
    const auto *hl_14 = buffer.data(hl + 14);
    const auto *hl_15 = buffer.data(hl + 15);
    const auto *hl_16 = buffer.data(hl + 16);
    const auto *hl_17 = buffer.data(hl + 17);
    const auto *hl_18 = buffer.data(hl + 18);
    const auto *hl_19 = buffer.data(hl + 19);
    const auto *hl_20 = buffer.data(hl + 20);
    const auto *hl_21 = buffer.data(hl + 21);
    const auto *hl_22 = buffer.data(hl + 22);
    const auto *hl_23 = buffer.data(hl + 23);
    const auto *hl_24 = buffer.data(hl + 24);
    const auto *hl_25 = buffer.data(hl + 25);
    const auto *hl_26 = buffer.data(hl + 26);
    const auto *hl_27 = buffer.data(hl + 27);
    const auto *hl_28 = buffer.data(hl + 28);
    const auto *hl_29 = buffer.data(hl + 29);
    const auto *hl_30 = buffer.data(hl + 30);
    const auto *hl_31 = buffer.data(hl + 31);
    const auto *hl_32 = buffer.data(hl + 32);
    const auto *hl_33 = buffer.data(hl + 33);
    const auto *hl_34 = buffer.data(hl + 34);
    const auto *hl_35 = buffer.data(hl + 35);
    const auto *hl_36 = buffer.data(hl + 36);
    const auto *hl_37 = buffer.data(hl + 37);
    const auto *hl_38 = buffer.data(hl + 38);
    const auto *hl_39 = buffer.data(hl + 39);
    const auto *hl_40 = buffer.data(hl + 40);
    const auto *hl_41 = buffer.data(hl + 41);
    const auto *hl_42 = buffer.data(hl + 42);
    const auto *hl_43 = buffer.data(hl + 43);
    const auto *hl_44 = buffer.data(hl + 44);
    const auto *hl_45 = buffer.data(hl + 45);
    const auto *hl_46 = buffer.data(hl + 46);
    const auto *hl_47 = buffer.data(hl + 47);
    const auto *hl_48 = buffer.data(hl + 48);
    const auto *hl_49 = buffer.data(hl + 49);
    const auto *hl_50 = buffer.data(hl + 50);
    const auto *hl_51 = buffer.data(hl + 51);
    const auto *hl_52 = buffer.data(hl + 52);
    const auto *hl_53 = buffer.data(hl + 53);
    const auto *hl_54 = buffer.data(hl + 54);
    const auto *hl_55 = buffer.data(hl + 55);
    const auto *hl_56 = buffer.data(hl + 56);
    const auto *hl_57 = buffer.data(hl + 57);
    const auto *hl_58 = buffer.data(hl + 58);
    const auto *hl_59 = buffer.data(hl + 59);
    const auto *hl_60 = buffer.data(hl + 60);
    const auto *hl_61 = buffer.data(hl + 61);
    const auto *hl_62 = buffer.data(hl + 62);
    const auto *hl_63 = buffer.data(hl + 63);
    const auto *hl_64 = buffer.data(hl + 64);
    const auto *hl_65 = buffer.data(hl + 65);
    const auto *hl_66 = buffer.data(hl + 66);
    const auto *hl_67 = buffer.data(hl + 67);
    const auto *hl_68 = buffer.data(hl + 68);
    const auto *hl_69 = buffer.data(hl + 69);
    const auto *hl_70 = buffer.data(hl + 70);
    const auto *hl_71 = buffer.data(hl + 71);
    const auto *hl_72 = buffer.data(hl + 72);
    const auto *hl_73 = buffer.data(hl + 73);
    const auto *hl_74 = buffer.data(hl + 74);
    const auto *hl_75 = buffer.data(hl + 75);
    const auto *hl_76 = buffer.data(hl + 76);
    const auto *hl_77 = buffer.data(hl + 77);
    const auto *hl_78 = buffer.data(hl + 78);
    const auto *hl_79 = buffer.data(hl + 79);
    const auto *hl_80 = buffer.data(hl + 80);
    const auto *hl_81 = buffer.data(hl + 81);
    const auto *hl_82 = buffer.data(hl + 82);
    const auto *hl_83 = buffer.data(hl + 83);
    const auto *hl_84 = buffer.data(hl + 84);
    const auto *hl_85 = buffer.data(hl + 85);
    const auto *hl_86 = buffer.data(hl + 86);
    const auto *hl_87 = buffer.data(hl + 87);
    const auto *hl_88 = buffer.data(hl + 88);
    const auto *hl_89 = buffer.data(hl + 89);
    const auto *hl_90 = buffer.data(hl + 90);
    const auto *hl_91 = buffer.data(hl + 91);
    const auto *hl_92 = buffer.data(hl + 92);
    const auto *hl_93 = buffer.data(hl + 93);
    const auto *hl_94 = buffer.data(hl + 94);
    const auto *hl_95 = buffer.data(hl + 95);
    const auto *hl_96 = buffer.data(hl + 96);
    const auto *hl_97 = buffer.data(hl + 97);
    const auto *hl_98 = buffer.data(hl + 98);
    const auto *hl_99 = buffer.data(hl + 99);
    const auto *hl_100 = buffer.data(hl + 100);
    const auto *hl_101 = buffer.data(hl + 101);
    const auto *hl_102 = buffer.data(hl + 102);
    const auto *hl_103 = buffer.data(hl + 103);
    const auto *hl_104 = buffer.data(hl + 104);
    const auto *hl_105 = buffer.data(hl + 105);
    const auto *hl_106 = buffer.data(hl + 106);
    const auto *hl_107 = buffer.data(hl + 107);
    const auto *hl_108 = buffer.data(hl + 108);
    const auto *hl_109 = buffer.data(hl + 109);
    const auto *hl_110 = buffer.data(hl + 110);
    const auto *hl_111 = buffer.data(hl + 111);
    const auto *hl_112 = buffer.data(hl + 112);
    const auto *hl_113 = buffer.data(hl + 113);
    const auto *hl_114 = buffer.data(hl + 114);
    const auto *hl_115 = buffer.data(hl + 115);
    const auto *hl_116 = buffer.data(hl + 116);
    const auto *hl_117 = buffer.data(hl + 117);
    const auto *hl_118 = buffer.data(hl + 118);
    const auto *hl_119 = buffer.data(hl + 119);
    const auto *hl_120 = buffer.data(hl + 120);
    const auto *hl_121 = buffer.data(hl + 121);
    const auto *hl_122 = buffer.data(hl + 122);
    const auto *hl_123 = buffer.data(hl + 123);
    const auto *hl_124 = buffer.data(hl + 124);
    const auto *hl_125 = buffer.data(hl + 125);
    const auto *hl_126 = buffer.data(hl + 126);
    const auto *hl_127 = buffer.data(hl + 127);
    const auto *hl_128 = buffer.data(hl + 128);
    const auto *hl_129 = buffer.data(hl + 129);
    const auto *hl_130 = buffer.data(hl + 130);
    const auto *hl_131 = buffer.data(hl + 131);
    const auto *hl_132 = buffer.data(hl + 132);
    const auto *hl_133 = buffer.data(hl + 133);
    const auto *hl_134 = buffer.data(hl + 134);
    const auto *hl_135 = buffer.data(hl + 135);
    const auto *hl_136 = buffer.data(hl + 136);
    const auto *hl_137 = buffer.data(hl + 137);
    const auto *hl_138 = buffer.data(hl + 138);
    const auto *hl_139 = buffer.data(hl + 139);
    const auto *hl_140 = buffer.data(hl + 140);
    const auto *hl_141 = buffer.data(hl + 141);
    const auto *hl_142 = buffer.data(hl + 142);
    const auto *hl_143 = buffer.data(hl + 143);
    const auto *hl_144 = buffer.data(hl + 144);
    const auto *hl_145 = buffer.data(hl + 145);
    const auto *hl_146 = buffer.data(hl + 146);
    const auto *hl_147 = buffer.data(hl + 147);
    const auto *hl_148 = buffer.data(hl + 148);
    const auto *hl_149 = buffer.data(hl + 149);

    const auto *kl_0 = buffer.data(kl + 0);
    const auto *kl_1 = buffer.data(kl + 1);
    const auto *kl_2 = buffer.data(kl + 2);
    const auto *kl_3 = buffer.data(kl + 3);
    const auto *kl_4 = buffer.data(kl + 4);
    const auto *kl_5 = buffer.data(kl + 5);
    const auto *kl_6 = buffer.data(kl + 6);
    const auto *kl_7 = buffer.data(kl + 7);
    const auto *kl_8 = buffer.data(kl + 8);
    const auto *kl_9 = buffer.data(kl + 9);
    const auto *kl_10 = buffer.data(kl + 10);
    const auto *kl_11 = buffer.data(kl + 11);
    const auto *kl_12 = buffer.data(kl + 12);
    const auto *kl_13 = buffer.data(kl + 13);
    const auto *kl_14 = buffer.data(kl + 14);
    const auto *kl_15 = buffer.data(kl + 15);
    const auto *kl_16 = buffer.data(kl + 16);
    const auto *kl_17 = buffer.data(kl + 17);
    const auto *kl_18 = buffer.data(kl + 18);
    const auto *kl_19 = buffer.data(kl + 19);
    const auto *kl_20 = buffer.data(kl + 20);
    const auto *kl_21 = buffer.data(kl + 21);
    const auto *kl_22 = buffer.data(kl + 22);
    const auto *kl_23 = buffer.data(kl + 23);
    const auto *kl_24 = buffer.data(kl + 24);
    const auto *kl_25 = buffer.data(kl + 25);
    const auto *kl_26 = buffer.data(kl + 26);
    const auto *kl_27 = buffer.data(kl + 27);
    const auto *kl_28 = buffer.data(kl + 28);
    const auto *kl_29 = buffer.data(kl + 29);
    const auto *kl_30 = buffer.data(kl + 30);
    const auto *kl_31 = buffer.data(kl + 31);
    const auto *kl_32 = buffer.data(kl + 32);
    const auto *kl_33 = buffer.data(kl + 33);
    const auto *kl_34 = buffer.data(kl + 34);
    const auto *kl_35 = buffer.data(kl + 35);
    const auto *kl_36 = buffer.data(kl + 36);
    const auto *kl_37 = buffer.data(kl + 37);
    const auto *kl_38 = buffer.data(kl + 38);
    const auto *kl_39 = buffer.data(kl + 39);
    const auto *kl_40 = buffer.data(kl + 40);
    const auto *kl_41 = buffer.data(kl + 41);
    const auto *kl_42 = buffer.data(kl + 42);
    const auto *kl_43 = buffer.data(kl + 43);
    const auto *kl_44 = buffer.data(kl + 44);
    const auto *kl_45 = buffer.data(kl + 45);
    const auto *kl_46 = buffer.data(kl + 46);
    const auto *kl_47 = buffer.data(kl + 47);
    const auto *kl_48 = buffer.data(kl + 48);
    const auto *kl_49 = buffer.data(kl + 49);
    const auto *kl_50 = buffer.data(kl + 50);
    const auto *kl_51 = buffer.data(kl + 51);
    const auto *kl_52 = buffer.data(kl + 52);
    const auto *kl_53 = buffer.data(kl + 53);
    const auto *kl_54 = buffer.data(kl + 54);
    const auto *kl_55 = buffer.data(kl + 55);
    const auto *kl_56 = buffer.data(kl + 56);
    const auto *kl_57 = buffer.data(kl + 57);
    const auto *kl_58 = buffer.data(kl + 58);
    const auto *kl_59 = buffer.data(kl + 59);
    const auto *kl_60 = buffer.data(kl + 60);
    const auto *kl_61 = buffer.data(kl + 61);
    const auto *kl_62 = buffer.data(kl + 62);
    const auto *kl_63 = buffer.data(kl + 63);
    const auto *kl_64 = buffer.data(kl + 64);
    const auto *kl_65 = buffer.data(kl + 65);
    const auto *kl_66 = buffer.data(kl + 66);
    const auto *kl_67 = buffer.data(kl + 67);
    const auto *kl_68 = buffer.data(kl + 68);
    const auto *kl_69 = buffer.data(kl + 69);
    const auto *kl_70 = buffer.data(kl + 70);
    const auto *kl_71 = buffer.data(kl + 71);
    const auto *kl_72 = buffer.data(kl + 72);
    const auto *kl_73 = buffer.data(kl + 73);
    const auto *kl_74 = buffer.data(kl + 74);
    const auto *kl_75 = buffer.data(kl + 75);
    const auto *kl_76 = buffer.data(kl + 76);
    const auto *kl_77 = buffer.data(kl + 77);
    const auto *kl_78 = buffer.data(kl + 78);
    const auto *kl_79 = buffer.data(kl + 79);
    const auto *kl_80 = buffer.data(kl + 80);
    const auto *kl_81 = buffer.data(kl + 81);
    const auto *kl_82 = buffer.data(kl + 82);
    const auto *kl_83 = buffer.data(kl + 83);
    const auto *kl_84 = buffer.data(kl + 84);
    const auto *kl_85 = buffer.data(kl + 85);
    const auto *kl_86 = buffer.data(kl + 86);
    const auto *kl_87 = buffer.data(kl + 87);
    const auto *kl_88 = buffer.data(kl + 88);
    const auto *kl_89 = buffer.data(kl + 89);
    const auto *kl_90 = buffer.data(kl + 90);
    const auto *kl_91 = buffer.data(kl + 91);
    const auto *kl_92 = buffer.data(kl + 92);
    const auto *kl_93 = buffer.data(kl + 93);
    const auto *kl_94 = buffer.data(kl + 94);
    const auto *kl_95 = buffer.data(kl + 95);
    const auto *kl_96 = buffer.data(kl + 96);
    const auto *kl_97 = buffer.data(kl + 97);
    const auto *kl_98 = buffer.data(kl + 98);
    const auto *kl_99 = buffer.data(kl + 99);
    const auto *kl_100 = buffer.data(kl + 100);
    const auto *kl_101 = buffer.data(kl + 101);
    const auto *kl_102 = buffer.data(kl + 102);
    const auto *kl_103 = buffer.data(kl + 103);
    const auto *kl_104 = buffer.data(kl + 104);
    const auto *kl_105 = buffer.data(kl + 105);
    const auto *kl_106 = buffer.data(kl + 106);
    const auto *kl_107 = buffer.data(kl + 107);
    const auto *kl_108 = buffer.data(kl + 108);
    const auto *kl_109 = buffer.data(kl + 109);
    const auto *kl_110 = buffer.data(kl + 110);
    const auto *kl_111 = buffer.data(kl + 111);
    const auto *kl_112 = buffer.data(kl + 112);
    const auto *kl_113 = buffer.data(kl + 113);
    const auto *kl_114 = buffer.data(kl + 114);
    const auto *kl_115 = buffer.data(kl + 115);
    const auto *kl_116 = buffer.data(kl + 116);
    const auto *kl_117 = buffer.data(kl + 117);
    const auto *kl_118 = buffer.data(kl + 118);
    const auto *kl_119 = buffer.data(kl + 119);
    const auto *kl_120 = buffer.data(kl + 120);
    const auto *kl_121 = buffer.data(kl + 121);
    const auto *kl_122 = buffer.data(kl + 122);
    const auto *kl_123 = buffer.data(kl + 123);
    const auto *kl_124 = buffer.data(kl + 124);
    const auto *kl_125 = buffer.data(kl + 125);
    const auto *kl_126 = buffer.data(kl + 126);
    const auto *kl_127 = buffer.data(kl + 127);
    const auto *kl_128 = buffer.data(kl + 128);
    const auto *kl_129 = buffer.data(kl + 129);
    const auto *kl_130 = buffer.data(kl + 130);
    const auto *kl_131 = buffer.data(kl + 131);
    const auto *kl_132 = buffer.data(kl + 132);
    const auto *kl_133 = buffer.data(kl + 133);
    const auto *kl_134 = buffer.data(kl + 134);
    const auto *kl_135 = buffer.data(kl + 135);
    const auto *kl_136 = buffer.data(kl + 136);
    const auto *kl_137 = buffer.data(kl + 137);
    const auto *kl_138 = buffer.data(kl + 138);
    const auto *kl_139 = buffer.data(kl + 139);
    const auto *kl_140 = buffer.data(kl + 140);
    const auto *kl_141 = buffer.data(kl + 141);
    const auto *kl_142 = buffer.data(kl + 142);
    const auto *kl_143 = buffer.data(kl + 143);
    const auto *kl_144 = buffer.data(kl + 144);
    const auto *kl_145 = buffer.data(kl + 145);
    const auto *kl_146 = buffer.data(kl + 146);
    const auto *kl_147 = buffer.data(kl + 147);
    const auto *kl_148 = buffer.data(kl + 148);
    const auto *kl_149 = buffer.data(kl + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, hl_0, hl_1, hl_2, hl_3, hl_4, kl_0, kl_1, \
                         kl_2, kl_3, kl_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -6.0 * hl_0[k]
                 + f_0 * kl_0[k];

        t_1[k] = -6.0 * hl_1[k]
                 + f_0 * kl_1[k];

        t_2[k] = -6.0 * hl_2[k]
                 + f_0 * kl_2[k];

        t_3[k] = -6.0 * hl_3[k]
                 + f_0 * kl_3[k];

        t_4[k] = -6.0 * hl_4[k]
                 + f_0 * kl_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, hl_5, hl_6, hl_7, hl_8, hl_9, kl_5, kl_6, \
                         kl_7, kl_8, kl_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -6.0 * hl_5[k]
                 + f_0 * kl_5[k];

        t_6[k] = -6.0 * hl_6[k]
                 + f_0 * kl_6[k];

        t_7[k] = -6.0 * hl_7[k]
                 + f_0 * kl_7[k];

        t_8[k] = -6.0 * hl_8[k]
                 + f_0 * kl_8[k];

        t_9[k] = -6.0 * hl_9[k]
                 + f_0 * kl_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, hl_10, hl_11, hl_12, hl_13, hl_14, \
                         kl_10, kl_11, kl_12, kl_13, kl_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -6.0 * hl_10[k]
                  + f_0 * kl_10[k];

        t_11[k] = -6.0 * hl_11[k]
                  + f_0 * kl_11[k];

        t_12[k] = -6.0 * hl_12[k]
                  + f_0 * kl_12[k];

        t_13[k] = -6.0 * hl_13[k]
                  + f_0 * kl_13[k];

        t_14[k] = -6.0 * hl_14[k]
                  + f_0 * kl_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, hl_15, hl_16, hl_17, hl_18, hl_19, \
                         kl_15, kl_16, kl_17, kl_18, kl_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -6.0 * hl_15[k]
                  + f_0 * kl_15[k];

        t_16[k] = -6.0 * hl_16[k]
                  + f_0 * kl_16[k];

        t_17[k] = -6.0 * hl_17[k]
                  + f_0 * kl_17[k];

        t_18[k] = -6.0 * hl_18[k]
                  + f_0 * kl_18[k];

        t_19[k] = -6.0 * hl_19[k]
                  + f_0 * kl_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, hl_20, hl_21, hl_22, hl_23, hl_24, \
                         kl_20, kl_21, kl_22, kl_23, kl_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -6.0 * hl_20[k]
                  + f_0 * kl_20[k];

        t_21[k] = -6.0 * hl_21[k]
                  + f_0 * kl_21[k];

        t_22[k] = -6.0 * hl_22[k]
                  + f_0 * kl_22[k];

        t_23[k] = -6.0 * hl_23[k]
                  + f_0 * kl_23[k];

        t_24[k] = -6.0 * hl_24[k]
                  + f_0 * kl_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, hl_25, hl_26, hl_27, hl_28, hl_29, \
                         kl_25, kl_26, kl_27, kl_28, kl_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -6.0 * hl_25[k]
                  + f_0 * kl_25[k];

        t_26[k] = -6.0 * hl_26[k]
                  + f_0 * kl_26[k];

        t_27[k] = -6.0 * hl_27[k]
                  + f_0 * kl_27[k];

        t_28[k] = -6.0 * hl_28[k]
                  + f_0 * kl_28[k];

        t_29[k] = -6.0 * hl_29[k]
                  + f_0 * kl_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, hl_30, hl_31, hl_32, hl_33, hl_34, \
                         kl_30, kl_31, kl_32, kl_33, kl_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -6.0 * hl_30[k]
                  + f_0 * kl_30[k];

        t_31[k] = -6.0 * hl_31[k]
                  + f_0 * kl_31[k];

        t_32[k] = -6.0 * hl_32[k]
                  + f_0 * kl_32[k];

        t_33[k] = -6.0 * hl_33[k]
                  + f_0 * kl_33[k];

        t_34[k] = -6.0 * hl_34[k]
                  + f_0 * kl_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, hl_35, hl_36, hl_37, hl_38, hl_39, \
                         kl_35, kl_36, kl_37, kl_38, kl_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -6.0 * hl_35[k]
                  + f_0 * kl_35[k];

        t_36[k] = -6.0 * hl_36[k]
                  + f_0 * kl_36[k];

        t_37[k] = -6.0 * hl_37[k]
                  + f_0 * kl_37[k];

        t_38[k] = -6.0 * hl_38[k]
                  + f_0 * kl_38[k];

        t_39[k] = -6.0 * hl_39[k]
                  + f_0 * kl_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, hl_40, hl_41, hl_42, hl_43, hl_44, \
                         kl_40, kl_41, kl_42, kl_43, kl_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -6.0 * hl_40[k]
                  + f_0 * kl_40[k];

        t_41[k] = -6.0 * hl_41[k]
                  + f_0 * kl_41[k];

        t_42[k] = -6.0 * hl_42[k]
                  + f_0 * kl_42[k];

        t_43[k] = -6.0 * hl_43[k]
                  + f_0 * kl_43[k];

        t_44[k] = -6.0 * hl_44[k]
                  + f_0 * kl_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hl_45, hl_46, hl_47, hl_48, hl_49, \
                         kl_45, kl_46, kl_47, kl_48, kl_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -5.0 * hl_45[k]
                  + f_0 * kl_45[k];

        t_46[k] = -5.0 * hl_46[k]
                  + f_0 * kl_46[k];

        t_47[k] = -5.0 * hl_47[k]
                  + f_0 * kl_47[k];

        t_48[k] = -5.0 * hl_48[k]
                  + f_0 * kl_48[k];

        t_49[k] = -5.0 * hl_49[k]
                  + f_0 * kl_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hl_50, hl_51, hl_52, hl_53, hl_54, \
                         kl_50, kl_51, kl_52, kl_53, kl_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -5.0 * hl_50[k]
                  + f_0 * kl_50[k];

        t_51[k] = -5.0 * hl_51[k]
                  + f_0 * kl_51[k];

        t_52[k] = -5.0 * hl_52[k]
                  + f_0 * kl_52[k];

        t_53[k] = -5.0 * hl_53[k]
                  + f_0 * kl_53[k];

        t_54[k] = -5.0 * hl_54[k]
                  + f_0 * kl_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, hl_55, hl_56, hl_57, hl_58, hl_59, \
                         kl_55, kl_56, kl_57, kl_58, kl_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -5.0 * hl_55[k]
                  + f_0 * kl_55[k];

        t_56[k] = -5.0 * hl_56[k]
                  + f_0 * kl_56[k];

        t_57[k] = -5.0 * hl_57[k]
                  + f_0 * kl_57[k];

        t_58[k] = -5.0 * hl_58[k]
                  + f_0 * kl_58[k];

        t_59[k] = -5.0 * hl_59[k]
                  + f_0 * kl_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, hl_60, hl_61, hl_62, hl_63, hl_64, \
                         kl_60, kl_61, kl_62, kl_63, kl_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -5.0 * hl_60[k]
                  + f_0 * kl_60[k];

        t_61[k] = -5.0 * hl_61[k]
                  + f_0 * kl_61[k];

        t_62[k] = -5.0 * hl_62[k]
                  + f_0 * kl_62[k];

        t_63[k] = -5.0 * hl_63[k]
                  + f_0 * kl_63[k];

        t_64[k] = -5.0 * hl_64[k]
                  + f_0 * kl_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, hl_65, hl_66, hl_67, hl_68, hl_69, \
                         kl_65, kl_66, kl_67, kl_68, kl_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -5.0 * hl_65[k]
                  + f_0 * kl_65[k];

        t_66[k] = -5.0 * hl_66[k]
                  + f_0 * kl_66[k];

        t_67[k] = -5.0 * hl_67[k]
                  + f_0 * kl_67[k];

        t_68[k] = -5.0 * hl_68[k]
                  + f_0 * kl_68[k];

        t_69[k] = -5.0 * hl_69[k]
                  + f_0 * kl_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, hl_70, hl_71, hl_72, hl_73, hl_74, \
                         kl_70, kl_71, kl_72, kl_73, kl_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -5.0 * hl_70[k]
                  + f_0 * kl_70[k];

        t_71[k] = -5.0 * hl_71[k]
                  + f_0 * kl_71[k];

        t_72[k] = -5.0 * hl_72[k]
                  + f_0 * kl_72[k];

        t_73[k] = -5.0 * hl_73[k]
                  + f_0 * kl_73[k];

        t_74[k] = -5.0 * hl_74[k]
                  + f_0 * kl_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, hl_75, hl_76, hl_77, hl_78, hl_79, \
                         kl_75, kl_76, kl_77, kl_78, kl_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -5.0 * hl_75[k]
                  + f_0 * kl_75[k];

        t_76[k] = -5.0 * hl_76[k]
                  + f_0 * kl_76[k];

        t_77[k] = -5.0 * hl_77[k]
                  + f_0 * kl_77[k];

        t_78[k] = -5.0 * hl_78[k]
                  + f_0 * kl_78[k];

        t_79[k] = -5.0 * hl_79[k]
                  + f_0 * kl_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, hl_80, hl_81, hl_82, hl_83, hl_84, \
                         kl_80, kl_81, kl_82, kl_83, kl_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -5.0 * hl_80[k]
                  + f_0 * kl_80[k];

        t_81[k] = -5.0 * hl_81[k]
                  + f_0 * kl_81[k];

        t_82[k] = -5.0 * hl_82[k]
                  + f_0 * kl_82[k];

        t_83[k] = -5.0 * hl_83[k]
                  + f_0 * kl_83[k];

        t_84[k] = -5.0 * hl_84[k]
                  + f_0 * kl_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, hl_85, hl_86, hl_87, hl_88, hl_89, \
                         kl_85, kl_86, kl_87, kl_88, kl_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -5.0 * hl_85[k]
                  + f_0 * kl_85[k];

        t_86[k] = -5.0 * hl_86[k]
                  + f_0 * kl_86[k];

        t_87[k] = -5.0 * hl_87[k]
                  + f_0 * kl_87[k];

        t_88[k] = -5.0 * hl_88[k]
                  + f_0 * kl_88[k];

        t_89[k] = -5.0 * hl_89[k]
                  + f_0 * kl_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, hl_90, hl_91, hl_92, hl_93, hl_94, \
                         kl_90, kl_91, kl_92, kl_93, kl_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -5.0 * hl_90[k]
                  + f_0 * kl_90[k];

        t_91[k] = -5.0 * hl_91[k]
                  + f_0 * kl_91[k];

        t_92[k] = -5.0 * hl_92[k]
                  + f_0 * kl_92[k];

        t_93[k] = -5.0 * hl_93[k]
                  + f_0 * kl_93[k];

        t_94[k] = -5.0 * hl_94[k]
                  + f_0 * kl_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, hl_95, hl_96, hl_97, hl_98, hl_99, \
                         kl_95, kl_96, kl_97, kl_98, kl_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -5.0 * hl_95[k]
                  + f_0 * kl_95[k];

        t_96[k] = -5.0 * hl_96[k]
                  + f_0 * kl_96[k];

        t_97[k] = -5.0 * hl_97[k]
                  + f_0 * kl_97[k];

        t_98[k] = -5.0 * hl_98[k]
                  + f_0 * kl_98[k];

        t_99[k] = -5.0 * hl_99[k]
                  + f_0 * kl_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, hl_100, hl_101, hl_102, hl_103, \
                         hl_104, kl_100, kl_101, kl_102, kl_103, \
                         kl_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -5.0 * hl_100[k]
                   + f_0 * kl_100[k];

        t_101[k] = -5.0 * hl_101[k]
                   + f_0 * kl_101[k];

        t_102[k] = -5.0 * hl_102[k]
                   + f_0 * kl_102[k];

        t_103[k] = -5.0 * hl_103[k]
                   + f_0 * kl_103[k];

        t_104[k] = -5.0 * hl_104[k]
                   + f_0 * kl_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, hl_105, hl_106, hl_107, hl_108, \
                         hl_109, kl_105, kl_106, kl_107, kl_108, \
                         kl_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -5.0 * hl_105[k]
                   + f_0 * kl_105[k];

        t_106[k] = -5.0 * hl_106[k]
                   + f_0 * kl_106[k];

        t_107[k] = -5.0 * hl_107[k]
                   + f_0 * kl_107[k];

        t_108[k] = -5.0 * hl_108[k]
                   + f_0 * kl_108[k];

        t_109[k] = -5.0 * hl_109[k]
                   + f_0 * kl_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, hl_110, hl_111, hl_112, hl_113, \
                         hl_114, kl_110, kl_111, kl_112, kl_113, \
                         kl_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -5.0 * hl_110[k]
                   + f_0 * kl_110[k];

        t_111[k] = -5.0 * hl_111[k]
                   + f_0 * kl_111[k];

        t_112[k] = -5.0 * hl_112[k]
                   + f_0 * kl_112[k];

        t_113[k] = -5.0 * hl_113[k]
                   + f_0 * kl_113[k];

        t_114[k] = -5.0 * hl_114[k]
                   + f_0 * kl_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, hl_115, hl_116, hl_117, hl_118, \
                         hl_119, kl_115, kl_116, kl_117, kl_118, \
                         kl_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -5.0 * hl_115[k]
                   + f_0 * kl_115[k];

        t_116[k] = -5.0 * hl_116[k]
                   + f_0 * kl_116[k];

        t_117[k] = -5.0 * hl_117[k]
                   + f_0 * kl_117[k];

        t_118[k] = -5.0 * hl_118[k]
                   + f_0 * kl_118[k];

        t_119[k] = -5.0 * hl_119[k]
                   + f_0 * kl_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, hl_120, hl_121, hl_122, hl_123, \
                         hl_124, kl_120, kl_121, kl_122, kl_123, \
                         kl_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -5.0 * hl_120[k]
                   + f_0 * kl_120[k];

        t_121[k] = -5.0 * hl_121[k]
                   + f_0 * kl_121[k];

        t_122[k] = -5.0 * hl_122[k]
                   + f_0 * kl_122[k];

        t_123[k] = -5.0 * hl_123[k]
                   + f_0 * kl_123[k];

        t_124[k] = -5.0 * hl_124[k]
                   + f_0 * kl_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, hl_125, hl_126, hl_127, hl_128, \
                         hl_129, kl_125, kl_126, kl_127, kl_128, \
                         kl_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -5.0 * hl_125[k]
                   + f_0 * kl_125[k];

        t_126[k] = -5.0 * hl_126[k]
                   + f_0 * kl_126[k];

        t_127[k] = -5.0 * hl_127[k]
                   + f_0 * kl_127[k];

        t_128[k] = -5.0 * hl_128[k]
                   + f_0 * kl_128[k];

        t_129[k] = -5.0 * hl_129[k]
                   + f_0 * kl_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, hl_130, hl_131, hl_132, hl_133, \
                         hl_134, kl_130, kl_131, kl_132, kl_133, \
                         kl_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -5.0 * hl_130[k]
                   + f_0 * kl_130[k];

        t_131[k] = -5.0 * hl_131[k]
                   + f_0 * kl_131[k];

        t_132[k] = -5.0 * hl_132[k]
                   + f_0 * kl_132[k];

        t_133[k] = -5.0 * hl_133[k]
                   + f_0 * kl_133[k];

        t_134[k] = -5.0 * hl_134[k]
                   + f_0 * kl_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, hl_135, hl_136, hl_137, hl_138, \
                         hl_139, kl_135, kl_136, kl_137, kl_138, \
                         kl_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -4.0 * hl_135[k]
                   + f_0 * kl_135[k];

        t_136[k] = -4.0 * hl_136[k]
                   + f_0 * kl_136[k];

        t_137[k] = -4.0 * hl_137[k]
                   + f_0 * kl_137[k];

        t_138[k] = -4.0 * hl_138[k]
                   + f_0 * kl_138[k];

        t_139[k] = -4.0 * hl_139[k]
                   + f_0 * kl_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, hl_140, hl_141, hl_142, hl_143, \
                         hl_144, kl_140, kl_141, kl_142, kl_143, \
                         kl_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -4.0 * hl_140[k]
                   + f_0 * kl_140[k];

        t_141[k] = -4.0 * hl_141[k]
                   + f_0 * kl_141[k];

        t_142[k] = -4.0 * hl_142[k]
                   + f_0 * kl_142[k];

        t_143[k] = -4.0 * hl_143[k]
                   + f_0 * kl_143[k];

        t_144[k] = -4.0 * hl_144[k]
                   + f_0 * kl_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, hl_145, hl_146, hl_147, hl_148, \
                         hl_149, kl_145, kl_146, kl_147, kl_148, \
                         kl_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -4.0 * hl_145[k]
                   + f_0 * kl_145[k];

        t_146[k] = -4.0 * hl_146[k]
                   + f_0 * kl_146[k];

        t_147[k] = -4.0 * hl_147[k]
                   + f_0 * kl_147[k];

        t_148[k] = -4.0 * hl_148[k]
                   + f_0 * kl_148[k];

        t_149[k] = -4.0 * hl_149[k]
                   + f_0 * kl_149[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
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

    const auto *hl_150 = buffer.data(hl + 150);
    const auto *hl_151 = buffer.data(hl + 151);
    const auto *hl_152 = buffer.data(hl + 152);
    const auto *hl_153 = buffer.data(hl + 153);
    const auto *hl_154 = buffer.data(hl + 154);
    const auto *hl_155 = buffer.data(hl + 155);
    const auto *hl_156 = buffer.data(hl + 156);
    const auto *hl_157 = buffer.data(hl + 157);
    const auto *hl_158 = buffer.data(hl + 158);
    const auto *hl_159 = buffer.data(hl + 159);
    const auto *hl_160 = buffer.data(hl + 160);
    const auto *hl_161 = buffer.data(hl + 161);
    const auto *hl_162 = buffer.data(hl + 162);
    const auto *hl_163 = buffer.data(hl + 163);
    const auto *hl_164 = buffer.data(hl + 164);
    const auto *hl_165 = buffer.data(hl + 165);
    const auto *hl_166 = buffer.data(hl + 166);
    const auto *hl_167 = buffer.data(hl + 167);
    const auto *hl_168 = buffer.data(hl + 168);
    const auto *hl_169 = buffer.data(hl + 169);
    const auto *hl_170 = buffer.data(hl + 170);
    const auto *hl_171 = buffer.data(hl + 171);
    const auto *hl_172 = buffer.data(hl + 172);
    const auto *hl_173 = buffer.data(hl + 173);
    const auto *hl_174 = buffer.data(hl + 174);
    const auto *hl_175 = buffer.data(hl + 175);
    const auto *hl_176 = buffer.data(hl + 176);
    const auto *hl_177 = buffer.data(hl + 177);
    const auto *hl_178 = buffer.data(hl + 178);
    const auto *hl_179 = buffer.data(hl + 179);
    const auto *hl_180 = buffer.data(hl + 180);
    const auto *hl_181 = buffer.data(hl + 181);
    const auto *hl_182 = buffer.data(hl + 182);
    const auto *hl_183 = buffer.data(hl + 183);
    const auto *hl_184 = buffer.data(hl + 184);
    const auto *hl_185 = buffer.data(hl + 185);
    const auto *hl_186 = buffer.data(hl + 186);
    const auto *hl_187 = buffer.data(hl + 187);
    const auto *hl_188 = buffer.data(hl + 188);
    const auto *hl_189 = buffer.data(hl + 189);
    const auto *hl_190 = buffer.data(hl + 190);
    const auto *hl_191 = buffer.data(hl + 191);
    const auto *hl_192 = buffer.data(hl + 192);
    const auto *hl_193 = buffer.data(hl + 193);
    const auto *hl_194 = buffer.data(hl + 194);
    const auto *hl_195 = buffer.data(hl + 195);
    const auto *hl_196 = buffer.data(hl + 196);
    const auto *hl_197 = buffer.data(hl + 197);
    const auto *hl_198 = buffer.data(hl + 198);
    const auto *hl_199 = buffer.data(hl + 199);
    const auto *hl_200 = buffer.data(hl + 200);
    const auto *hl_201 = buffer.data(hl + 201);
    const auto *hl_202 = buffer.data(hl + 202);
    const auto *hl_203 = buffer.data(hl + 203);
    const auto *hl_204 = buffer.data(hl + 204);
    const auto *hl_205 = buffer.data(hl + 205);
    const auto *hl_206 = buffer.data(hl + 206);
    const auto *hl_207 = buffer.data(hl + 207);
    const auto *hl_208 = buffer.data(hl + 208);
    const auto *hl_209 = buffer.data(hl + 209);
    const auto *hl_210 = buffer.data(hl + 210);
    const auto *hl_211 = buffer.data(hl + 211);
    const auto *hl_212 = buffer.data(hl + 212);
    const auto *hl_213 = buffer.data(hl + 213);
    const auto *hl_214 = buffer.data(hl + 214);
    const auto *hl_215 = buffer.data(hl + 215);
    const auto *hl_216 = buffer.data(hl + 216);
    const auto *hl_217 = buffer.data(hl + 217);
    const auto *hl_218 = buffer.data(hl + 218);
    const auto *hl_219 = buffer.data(hl + 219);
    const auto *hl_220 = buffer.data(hl + 220);
    const auto *hl_221 = buffer.data(hl + 221);
    const auto *hl_222 = buffer.data(hl + 222);
    const auto *hl_223 = buffer.data(hl + 223);
    const auto *hl_224 = buffer.data(hl + 224);
    const auto *hl_225 = buffer.data(hl + 225);
    const auto *hl_226 = buffer.data(hl + 226);
    const auto *hl_227 = buffer.data(hl + 227);
    const auto *hl_228 = buffer.data(hl + 228);
    const auto *hl_229 = buffer.data(hl + 229);
    const auto *hl_230 = buffer.data(hl + 230);
    const auto *hl_231 = buffer.data(hl + 231);
    const auto *hl_232 = buffer.data(hl + 232);
    const auto *hl_233 = buffer.data(hl + 233);
    const auto *hl_234 = buffer.data(hl + 234);
    const auto *hl_235 = buffer.data(hl + 235);
    const auto *hl_236 = buffer.data(hl + 236);
    const auto *hl_237 = buffer.data(hl + 237);
    const auto *hl_238 = buffer.data(hl + 238);
    const auto *hl_239 = buffer.data(hl + 239);
    const auto *hl_240 = buffer.data(hl + 240);
    const auto *hl_241 = buffer.data(hl + 241);
    const auto *hl_242 = buffer.data(hl + 242);
    const auto *hl_243 = buffer.data(hl + 243);
    const auto *hl_244 = buffer.data(hl + 244);
    const auto *hl_245 = buffer.data(hl + 245);
    const auto *hl_246 = buffer.data(hl + 246);
    const auto *hl_247 = buffer.data(hl + 247);
    const auto *hl_248 = buffer.data(hl + 248);
    const auto *hl_249 = buffer.data(hl + 249);
    const auto *hl_250 = buffer.data(hl + 250);
    const auto *hl_251 = buffer.data(hl + 251);
    const auto *hl_252 = buffer.data(hl + 252);
    const auto *hl_253 = buffer.data(hl + 253);
    const auto *hl_254 = buffer.data(hl + 254);
    const auto *hl_255 = buffer.data(hl + 255);
    const auto *hl_256 = buffer.data(hl + 256);
    const auto *hl_257 = buffer.data(hl + 257);
    const auto *hl_258 = buffer.data(hl + 258);
    const auto *hl_259 = buffer.data(hl + 259);
    const auto *hl_260 = buffer.data(hl + 260);
    const auto *hl_261 = buffer.data(hl + 261);
    const auto *hl_262 = buffer.data(hl + 262);
    const auto *hl_263 = buffer.data(hl + 263);
    const auto *hl_264 = buffer.data(hl + 264);
    const auto *hl_265 = buffer.data(hl + 265);
    const auto *hl_266 = buffer.data(hl + 266);
    const auto *hl_267 = buffer.data(hl + 267);
    const auto *hl_268 = buffer.data(hl + 268);
    const auto *hl_269 = buffer.data(hl + 269);
    const auto *hl_270 = buffer.data(hl + 270);
    const auto *hl_271 = buffer.data(hl + 271);
    const auto *hl_272 = buffer.data(hl + 272);
    const auto *hl_273 = buffer.data(hl + 273);
    const auto *hl_274 = buffer.data(hl + 274);
    const auto *hl_275 = buffer.data(hl + 275);
    const auto *hl_276 = buffer.data(hl + 276);
    const auto *hl_277 = buffer.data(hl + 277);
    const auto *hl_278 = buffer.data(hl + 278);
    const auto *hl_279 = buffer.data(hl + 279);
    const auto *hl_280 = buffer.data(hl + 280);
    const auto *hl_281 = buffer.data(hl + 281);
    const auto *hl_282 = buffer.data(hl + 282);
    const auto *hl_283 = buffer.data(hl + 283);
    const auto *hl_284 = buffer.data(hl + 284);
    const auto *hl_285 = buffer.data(hl + 285);
    const auto *hl_286 = buffer.data(hl + 286);
    const auto *hl_287 = buffer.data(hl + 287);
    const auto *hl_288 = buffer.data(hl + 288);
    const auto *hl_289 = buffer.data(hl + 289);
    const auto *hl_290 = buffer.data(hl + 290);
    const auto *hl_291 = buffer.data(hl + 291);
    const auto *hl_292 = buffer.data(hl + 292);
    const auto *hl_293 = buffer.data(hl + 293);
    const auto *hl_294 = buffer.data(hl + 294);
    const auto *hl_295 = buffer.data(hl + 295);
    const auto *hl_296 = buffer.data(hl + 296);
    const auto *hl_297 = buffer.data(hl + 297);
    const auto *hl_298 = buffer.data(hl + 298);
    const auto *hl_299 = buffer.data(hl + 299);

    const auto *kl_150 = buffer.data(kl + 150);
    const auto *kl_151 = buffer.data(kl + 151);
    const auto *kl_152 = buffer.data(kl + 152);
    const auto *kl_153 = buffer.data(kl + 153);
    const auto *kl_154 = buffer.data(kl + 154);
    const auto *kl_155 = buffer.data(kl + 155);
    const auto *kl_156 = buffer.data(kl + 156);
    const auto *kl_157 = buffer.data(kl + 157);
    const auto *kl_158 = buffer.data(kl + 158);
    const auto *kl_159 = buffer.data(kl + 159);
    const auto *kl_160 = buffer.data(kl + 160);
    const auto *kl_161 = buffer.data(kl + 161);
    const auto *kl_162 = buffer.data(kl + 162);
    const auto *kl_163 = buffer.data(kl + 163);
    const auto *kl_164 = buffer.data(kl + 164);
    const auto *kl_165 = buffer.data(kl + 165);
    const auto *kl_166 = buffer.data(kl + 166);
    const auto *kl_167 = buffer.data(kl + 167);
    const auto *kl_168 = buffer.data(kl + 168);
    const auto *kl_169 = buffer.data(kl + 169);
    const auto *kl_170 = buffer.data(kl + 170);
    const auto *kl_171 = buffer.data(kl + 171);
    const auto *kl_172 = buffer.data(kl + 172);
    const auto *kl_173 = buffer.data(kl + 173);
    const auto *kl_174 = buffer.data(kl + 174);
    const auto *kl_175 = buffer.data(kl + 175);
    const auto *kl_176 = buffer.data(kl + 176);
    const auto *kl_177 = buffer.data(kl + 177);
    const auto *kl_178 = buffer.data(kl + 178);
    const auto *kl_179 = buffer.data(kl + 179);
    const auto *kl_180 = buffer.data(kl + 180);
    const auto *kl_181 = buffer.data(kl + 181);
    const auto *kl_182 = buffer.data(kl + 182);
    const auto *kl_183 = buffer.data(kl + 183);
    const auto *kl_184 = buffer.data(kl + 184);
    const auto *kl_185 = buffer.data(kl + 185);
    const auto *kl_186 = buffer.data(kl + 186);
    const auto *kl_187 = buffer.data(kl + 187);
    const auto *kl_188 = buffer.data(kl + 188);
    const auto *kl_189 = buffer.data(kl + 189);
    const auto *kl_190 = buffer.data(kl + 190);
    const auto *kl_191 = buffer.data(kl + 191);
    const auto *kl_192 = buffer.data(kl + 192);
    const auto *kl_193 = buffer.data(kl + 193);
    const auto *kl_194 = buffer.data(kl + 194);
    const auto *kl_195 = buffer.data(kl + 195);
    const auto *kl_196 = buffer.data(kl + 196);
    const auto *kl_197 = buffer.data(kl + 197);
    const auto *kl_198 = buffer.data(kl + 198);
    const auto *kl_199 = buffer.data(kl + 199);
    const auto *kl_200 = buffer.data(kl + 200);
    const auto *kl_201 = buffer.data(kl + 201);
    const auto *kl_202 = buffer.data(kl + 202);
    const auto *kl_203 = buffer.data(kl + 203);
    const auto *kl_204 = buffer.data(kl + 204);
    const auto *kl_205 = buffer.data(kl + 205);
    const auto *kl_206 = buffer.data(kl + 206);
    const auto *kl_207 = buffer.data(kl + 207);
    const auto *kl_208 = buffer.data(kl + 208);
    const auto *kl_209 = buffer.data(kl + 209);
    const auto *kl_210 = buffer.data(kl + 210);
    const auto *kl_211 = buffer.data(kl + 211);
    const auto *kl_212 = buffer.data(kl + 212);
    const auto *kl_213 = buffer.data(kl + 213);
    const auto *kl_214 = buffer.data(kl + 214);
    const auto *kl_215 = buffer.data(kl + 215);
    const auto *kl_216 = buffer.data(kl + 216);
    const auto *kl_217 = buffer.data(kl + 217);
    const auto *kl_218 = buffer.data(kl + 218);
    const auto *kl_219 = buffer.data(kl + 219);
    const auto *kl_220 = buffer.data(kl + 220);
    const auto *kl_221 = buffer.data(kl + 221);
    const auto *kl_222 = buffer.data(kl + 222);
    const auto *kl_223 = buffer.data(kl + 223);
    const auto *kl_224 = buffer.data(kl + 224);
    const auto *kl_225 = buffer.data(kl + 225);
    const auto *kl_226 = buffer.data(kl + 226);
    const auto *kl_227 = buffer.data(kl + 227);
    const auto *kl_228 = buffer.data(kl + 228);
    const auto *kl_229 = buffer.data(kl + 229);
    const auto *kl_230 = buffer.data(kl + 230);
    const auto *kl_231 = buffer.data(kl + 231);
    const auto *kl_232 = buffer.data(kl + 232);
    const auto *kl_233 = buffer.data(kl + 233);
    const auto *kl_234 = buffer.data(kl + 234);
    const auto *kl_235 = buffer.data(kl + 235);
    const auto *kl_236 = buffer.data(kl + 236);
    const auto *kl_237 = buffer.data(kl + 237);
    const auto *kl_238 = buffer.data(kl + 238);
    const auto *kl_239 = buffer.data(kl + 239);
    const auto *kl_240 = buffer.data(kl + 240);
    const auto *kl_241 = buffer.data(kl + 241);
    const auto *kl_242 = buffer.data(kl + 242);
    const auto *kl_243 = buffer.data(kl + 243);
    const auto *kl_244 = buffer.data(kl + 244);
    const auto *kl_245 = buffer.data(kl + 245);
    const auto *kl_246 = buffer.data(kl + 246);
    const auto *kl_247 = buffer.data(kl + 247);
    const auto *kl_248 = buffer.data(kl + 248);
    const auto *kl_249 = buffer.data(kl + 249);
    const auto *kl_250 = buffer.data(kl + 250);
    const auto *kl_251 = buffer.data(kl + 251);
    const auto *kl_252 = buffer.data(kl + 252);
    const auto *kl_253 = buffer.data(kl + 253);
    const auto *kl_254 = buffer.data(kl + 254);
    const auto *kl_255 = buffer.data(kl + 255);
    const auto *kl_256 = buffer.data(kl + 256);
    const auto *kl_257 = buffer.data(kl + 257);
    const auto *kl_258 = buffer.data(kl + 258);
    const auto *kl_259 = buffer.data(kl + 259);
    const auto *kl_260 = buffer.data(kl + 260);
    const auto *kl_261 = buffer.data(kl + 261);
    const auto *kl_262 = buffer.data(kl + 262);
    const auto *kl_263 = buffer.data(kl + 263);
    const auto *kl_264 = buffer.data(kl + 264);
    const auto *kl_265 = buffer.data(kl + 265);
    const auto *kl_266 = buffer.data(kl + 266);
    const auto *kl_267 = buffer.data(kl + 267);
    const auto *kl_268 = buffer.data(kl + 268);
    const auto *kl_269 = buffer.data(kl + 269);
    const auto *kl_270 = buffer.data(kl + 270);
    const auto *kl_271 = buffer.data(kl + 271);
    const auto *kl_272 = buffer.data(kl + 272);
    const auto *kl_273 = buffer.data(kl + 273);
    const auto *kl_274 = buffer.data(kl + 274);
    const auto *kl_275 = buffer.data(kl + 275);
    const auto *kl_276 = buffer.data(kl + 276);
    const auto *kl_277 = buffer.data(kl + 277);
    const auto *kl_278 = buffer.data(kl + 278);
    const auto *kl_279 = buffer.data(kl + 279);
    const auto *kl_280 = buffer.data(kl + 280);
    const auto *kl_281 = buffer.data(kl + 281);
    const auto *kl_282 = buffer.data(kl + 282);
    const auto *kl_283 = buffer.data(kl + 283);
    const auto *kl_284 = buffer.data(kl + 284);
    const auto *kl_285 = buffer.data(kl + 285);
    const auto *kl_286 = buffer.data(kl + 286);
    const auto *kl_287 = buffer.data(kl + 287);
    const auto *kl_288 = buffer.data(kl + 288);
    const auto *kl_289 = buffer.data(kl + 289);
    const auto *kl_290 = buffer.data(kl + 290);
    const auto *kl_291 = buffer.data(kl + 291);
    const auto *kl_292 = buffer.data(kl + 292);
    const auto *kl_293 = buffer.data(kl + 293);
    const auto *kl_294 = buffer.data(kl + 294);
    const auto *kl_295 = buffer.data(kl + 295);
    const auto *kl_296 = buffer.data(kl + 296);
    const auto *kl_297 = buffer.data(kl + 297);
    const auto *kl_298 = buffer.data(kl + 298);
    const auto *kl_299 = buffer.data(kl + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, hl_150, hl_151, hl_152, hl_153, \
                         hl_154, kl_150, kl_151, kl_152, kl_153, \
                         kl_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -4.0 * hl_150[k]
                   + f_0 * kl_150[k];

        t_151[k] = -4.0 * hl_151[k]
                   + f_0 * kl_151[k];

        t_152[k] = -4.0 * hl_152[k]
                   + f_0 * kl_152[k];

        t_153[k] = -4.0 * hl_153[k]
                   + f_0 * kl_153[k];

        t_154[k] = -4.0 * hl_154[k]
                   + f_0 * kl_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, hl_155, hl_156, hl_157, hl_158, \
                         hl_159, kl_155, kl_156, kl_157, kl_158, \
                         kl_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -4.0 * hl_155[k]
                   + f_0 * kl_155[k];

        t_156[k] = -4.0 * hl_156[k]
                   + f_0 * kl_156[k];

        t_157[k] = -4.0 * hl_157[k]
                   + f_0 * kl_157[k];

        t_158[k] = -4.0 * hl_158[k]
                   + f_0 * kl_158[k];

        t_159[k] = -4.0 * hl_159[k]
                   + f_0 * kl_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, hl_160, hl_161, hl_162, hl_163, \
                         hl_164, kl_160, kl_161, kl_162, kl_163, \
                         kl_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -4.0 * hl_160[k]
                   + f_0 * kl_160[k];

        t_161[k] = -4.0 * hl_161[k]
                   + f_0 * kl_161[k];

        t_162[k] = -4.0 * hl_162[k]
                   + f_0 * kl_162[k];

        t_163[k] = -4.0 * hl_163[k]
                   + f_0 * kl_163[k];

        t_164[k] = -4.0 * hl_164[k]
                   + f_0 * kl_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, hl_165, hl_166, hl_167, hl_168, \
                         hl_169, kl_165, kl_166, kl_167, kl_168, \
                         kl_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -4.0 * hl_165[k]
                   + f_0 * kl_165[k];

        t_166[k] = -4.0 * hl_166[k]
                   + f_0 * kl_166[k];

        t_167[k] = -4.0 * hl_167[k]
                   + f_0 * kl_167[k];

        t_168[k] = -4.0 * hl_168[k]
                   + f_0 * kl_168[k];

        t_169[k] = -4.0 * hl_169[k]
                   + f_0 * kl_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, hl_170, hl_171, hl_172, hl_173, \
                         hl_174, kl_170, kl_171, kl_172, kl_173, \
                         kl_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -4.0 * hl_170[k]
                   + f_0 * kl_170[k];

        t_171[k] = -4.0 * hl_171[k]
                   + f_0 * kl_171[k];

        t_172[k] = -4.0 * hl_172[k]
                   + f_0 * kl_172[k];

        t_173[k] = -4.0 * hl_173[k]
                   + f_0 * kl_173[k];

        t_174[k] = -4.0 * hl_174[k]
                   + f_0 * kl_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, hl_175, hl_176, hl_177, hl_178, \
                         hl_179, kl_175, kl_176, kl_177, kl_178, \
                         kl_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -4.0 * hl_175[k]
                   + f_0 * kl_175[k];

        t_176[k] = -4.0 * hl_176[k]
                   + f_0 * kl_176[k];

        t_177[k] = -4.0 * hl_177[k]
                   + f_0 * kl_177[k];

        t_178[k] = -4.0 * hl_178[k]
                   + f_0 * kl_178[k];

        t_179[k] = -4.0 * hl_179[k]
                   + f_0 * kl_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, hl_180, hl_181, hl_182, hl_183, \
                         hl_184, kl_180, kl_181, kl_182, kl_183, \
                         kl_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -4.0 * hl_180[k]
                   + f_0 * kl_180[k];

        t_181[k] = -4.0 * hl_181[k]
                   + f_0 * kl_181[k];

        t_182[k] = -4.0 * hl_182[k]
                   + f_0 * kl_182[k];

        t_183[k] = -4.0 * hl_183[k]
                   + f_0 * kl_183[k];

        t_184[k] = -4.0 * hl_184[k]
                   + f_0 * kl_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, hl_185, hl_186, hl_187, hl_188, \
                         hl_189, kl_185, kl_186, kl_187, kl_188, \
                         kl_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -4.0 * hl_185[k]
                   + f_0 * kl_185[k];

        t_186[k] = -4.0 * hl_186[k]
                   + f_0 * kl_186[k];

        t_187[k] = -4.0 * hl_187[k]
                   + f_0 * kl_187[k];

        t_188[k] = -4.0 * hl_188[k]
                   + f_0 * kl_188[k];

        t_189[k] = -4.0 * hl_189[k]
                   + f_0 * kl_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, hl_190, hl_191, hl_192, hl_193, \
                         hl_194, kl_190, kl_191, kl_192, kl_193, \
                         kl_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -4.0 * hl_190[k]
                   + f_0 * kl_190[k];

        t_191[k] = -4.0 * hl_191[k]
                   + f_0 * kl_191[k];

        t_192[k] = -4.0 * hl_192[k]
                   + f_0 * kl_192[k];

        t_193[k] = -4.0 * hl_193[k]
                   + f_0 * kl_193[k];

        t_194[k] = -4.0 * hl_194[k]
                   + f_0 * kl_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, hl_195, hl_196, hl_197, hl_198, \
                         hl_199, kl_195, kl_196, kl_197, kl_198, \
                         kl_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -4.0 * hl_195[k]
                   + f_0 * kl_195[k];

        t_196[k] = -4.0 * hl_196[k]
                   + f_0 * kl_196[k];

        t_197[k] = -4.0 * hl_197[k]
                   + f_0 * kl_197[k];

        t_198[k] = -4.0 * hl_198[k]
                   + f_0 * kl_198[k];

        t_199[k] = -4.0 * hl_199[k]
                   + f_0 * kl_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, hl_200, hl_201, hl_202, hl_203, \
                         hl_204, kl_200, kl_201, kl_202, kl_203, \
                         kl_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -4.0 * hl_200[k]
                   + f_0 * kl_200[k];

        t_201[k] = -4.0 * hl_201[k]
                   + f_0 * kl_201[k];

        t_202[k] = -4.0 * hl_202[k]
                   + f_0 * kl_202[k];

        t_203[k] = -4.0 * hl_203[k]
                   + f_0 * kl_203[k];

        t_204[k] = -4.0 * hl_204[k]
                   + f_0 * kl_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, hl_205, hl_206, hl_207, hl_208, \
                         hl_209, kl_205, kl_206, kl_207, kl_208, \
                         kl_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -4.0 * hl_205[k]
                   + f_0 * kl_205[k];

        t_206[k] = -4.0 * hl_206[k]
                   + f_0 * kl_206[k];

        t_207[k] = -4.0 * hl_207[k]
                   + f_0 * kl_207[k];

        t_208[k] = -4.0 * hl_208[k]
                   + f_0 * kl_208[k];

        t_209[k] = -4.0 * hl_209[k]
                   + f_0 * kl_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, hl_210, hl_211, hl_212, hl_213, \
                         hl_214, kl_210, kl_211, kl_212, kl_213, \
                         kl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -4.0 * hl_210[k]
                   + f_0 * kl_210[k];

        t_211[k] = -4.0 * hl_211[k]
                   + f_0 * kl_211[k];

        t_212[k] = -4.0 * hl_212[k]
                   + f_0 * kl_212[k];

        t_213[k] = -4.0 * hl_213[k]
                   + f_0 * kl_213[k];

        t_214[k] = -4.0 * hl_214[k]
                   + f_0 * kl_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, hl_215, hl_216, hl_217, hl_218, \
                         hl_219, kl_215, kl_216, kl_217, kl_218, \
                         kl_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -4.0 * hl_215[k]
                   + f_0 * kl_215[k];

        t_216[k] = -4.0 * hl_216[k]
                   + f_0 * kl_216[k];

        t_217[k] = -4.0 * hl_217[k]
                   + f_0 * kl_217[k];

        t_218[k] = -4.0 * hl_218[k]
                   + f_0 * kl_218[k];

        t_219[k] = -4.0 * hl_219[k]
                   + f_0 * kl_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, hl_220, hl_221, hl_222, hl_223, \
                         hl_224, kl_220, kl_221, kl_222, kl_223, \
                         kl_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -4.0 * hl_220[k]
                   + f_0 * kl_220[k];

        t_221[k] = -4.0 * hl_221[k]
                   + f_0 * kl_221[k];

        t_222[k] = -4.0 * hl_222[k]
                   + f_0 * kl_222[k];

        t_223[k] = -4.0 * hl_223[k]
                   + f_0 * kl_223[k];

        t_224[k] = -4.0 * hl_224[k]
                   + f_0 * kl_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, hl_225, hl_226, hl_227, hl_228, \
                         hl_229, kl_225, kl_226, kl_227, kl_228, \
                         kl_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -4.0 * hl_225[k]
                   + f_0 * kl_225[k];

        t_226[k] = -4.0 * hl_226[k]
                   + f_0 * kl_226[k];

        t_227[k] = -4.0 * hl_227[k]
                   + f_0 * kl_227[k];

        t_228[k] = -4.0 * hl_228[k]
                   + f_0 * kl_228[k];

        t_229[k] = -4.0 * hl_229[k]
                   + f_0 * kl_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, hl_230, hl_231, hl_232, hl_233, \
                         hl_234, kl_230, kl_231, kl_232, kl_233, \
                         kl_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -4.0 * hl_230[k]
                   + f_0 * kl_230[k];

        t_231[k] = -4.0 * hl_231[k]
                   + f_0 * kl_231[k];

        t_232[k] = -4.0 * hl_232[k]
                   + f_0 * kl_232[k];

        t_233[k] = -4.0 * hl_233[k]
                   + f_0 * kl_233[k];

        t_234[k] = -4.0 * hl_234[k]
                   + f_0 * kl_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, hl_235, hl_236, hl_237, hl_238, \
                         hl_239, kl_235, kl_236, kl_237, kl_238, \
                         kl_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -4.0 * hl_235[k]
                   + f_0 * kl_235[k];

        t_236[k] = -4.0 * hl_236[k]
                   + f_0 * kl_236[k];

        t_237[k] = -4.0 * hl_237[k]
                   + f_0 * kl_237[k];

        t_238[k] = -4.0 * hl_238[k]
                   + f_0 * kl_238[k];

        t_239[k] = -4.0 * hl_239[k]
                   + f_0 * kl_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, hl_240, hl_241, hl_242, hl_243, \
                         hl_244, kl_240, kl_241, kl_242, kl_243, \
                         kl_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -4.0 * hl_240[k]
                   + f_0 * kl_240[k];

        t_241[k] = -4.0 * hl_241[k]
                   + f_0 * kl_241[k];

        t_242[k] = -4.0 * hl_242[k]
                   + f_0 * kl_242[k];

        t_243[k] = -4.0 * hl_243[k]
                   + f_0 * kl_243[k];

        t_244[k] = -4.0 * hl_244[k]
                   + f_0 * kl_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, hl_245, hl_246, hl_247, hl_248, \
                         hl_249, kl_245, kl_246, kl_247, kl_248, \
                         kl_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -4.0 * hl_245[k]
                   + f_0 * kl_245[k];

        t_246[k] = -4.0 * hl_246[k]
                   + f_0 * kl_246[k];

        t_247[k] = -4.0 * hl_247[k]
                   + f_0 * kl_247[k];

        t_248[k] = -4.0 * hl_248[k]
                   + f_0 * kl_248[k];

        t_249[k] = -4.0 * hl_249[k]
                   + f_0 * kl_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, hl_250, hl_251, hl_252, hl_253, \
                         hl_254, kl_250, kl_251, kl_252, kl_253, \
                         kl_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -4.0 * hl_250[k]
                   + f_0 * kl_250[k];

        t_251[k] = -4.0 * hl_251[k]
                   + f_0 * kl_251[k];

        t_252[k] = -4.0 * hl_252[k]
                   + f_0 * kl_252[k];

        t_253[k] = -4.0 * hl_253[k]
                   + f_0 * kl_253[k];

        t_254[k] = -4.0 * hl_254[k]
                   + f_0 * kl_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, hl_255, hl_256, hl_257, hl_258, \
                         hl_259, kl_255, kl_256, kl_257, kl_258, \
                         kl_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -4.0 * hl_255[k]
                   + f_0 * kl_255[k];

        t_256[k] = -4.0 * hl_256[k]
                   + f_0 * kl_256[k];

        t_257[k] = -4.0 * hl_257[k]
                   + f_0 * kl_257[k];

        t_258[k] = -4.0 * hl_258[k]
                   + f_0 * kl_258[k];

        t_259[k] = -4.0 * hl_259[k]
                   + f_0 * kl_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, hl_260, hl_261, hl_262, hl_263, \
                         hl_264, kl_260, kl_261, kl_262, kl_263, \
                         kl_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -4.0 * hl_260[k]
                   + f_0 * kl_260[k];

        t_261[k] = -4.0 * hl_261[k]
                   + f_0 * kl_261[k];

        t_262[k] = -4.0 * hl_262[k]
                   + f_0 * kl_262[k];

        t_263[k] = -4.0 * hl_263[k]
                   + f_0 * kl_263[k];

        t_264[k] = -4.0 * hl_264[k]
                   + f_0 * kl_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, hl_265, hl_266, hl_267, hl_268, \
                         hl_269, kl_265, kl_266, kl_267, kl_268, \
                         kl_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -4.0 * hl_265[k]
                   + f_0 * kl_265[k];

        t_266[k] = -4.0 * hl_266[k]
                   + f_0 * kl_266[k];

        t_267[k] = -4.0 * hl_267[k]
                   + f_0 * kl_267[k];

        t_268[k] = -4.0 * hl_268[k]
                   + f_0 * kl_268[k];

        t_269[k] = -4.0 * hl_269[k]
                   + f_0 * kl_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, hl_270, hl_271, hl_272, hl_273, \
                         hl_274, kl_270, kl_271, kl_272, kl_273, \
                         kl_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -3.0 * hl_270[k]
                   + f_0 * kl_270[k];

        t_271[k] = -3.0 * hl_271[k]
                   + f_0 * kl_271[k];

        t_272[k] = -3.0 * hl_272[k]
                   + f_0 * kl_272[k];

        t_273[k] = -3.0 * hl_273[k]
                   + f_0 * kl_273[k];

        t_274[k] = -3.0 * hl_274[k]
                   + f_0 * kl_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, hl_275, hl_276, hl_277, hl_278, \
                         hl_279, kl_275, kl_276, kl_277, kl_278, \
                         kl_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * hl_275[k]
                   + f_0 * kl_275[k];

        t_276[k] = -3.0 * hl_276[k]
                   + f_0 * kl_276[k];

        t_277[k] = -3.0 * hl_277[k]
                   + f_0 * kl_277[k];

        t_278[k] = -3.0 * hl_278[k]
                   + f_0 * kl_278[k];

        t_279[k] = -3.0 * hl_279[k]
                   + f_0 * kl_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, hl_280, hl_281, hl_282, hl_283, \
                         hl_284, kl_280, kl_281, kl_282, kl_283, \
                         kl_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * hl_280[k]
                   + f_0 * kl_280[k];

        t_281[k] = -3.0 * hl_281[k]
                   + f_0 * kl_281[k];

        t_282[k] = -3.0 * hl_282[k]
                   + f_0 * kl_282[k];

        t_283[k] = -3.0 * hl_283[k]
                   + f_0 * kl_283[k];

        t_284[k] = -3.0 * hl_284[k]
                   + f_0 * kl_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, hl_285, hl_286, hl_287, hl_288, \
                         hl_289, kl_285, kl_286, kl_287, kl_288, \
                         kl_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * hl_285[k]
                   + f_0 * kl_285[k];

        t_286[k] = -3.0 * hl_286[k]
                   + f_0 * kl_286[k];

        t_287[k] = -3.0 * hl_287[k]
                   + f_0 * kl_287[k];

        t_288[k] = -3.0 * hl_288[k]
                   + f_0 * kl_288[k];

        t_289[k] = -3.0 * hl_289[k]
                   + f_0 * kl_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, hl_290, hl_291, hl_292, hl_293, \
                         hl_294, kl_290, kl_291, kl_292, kl_293, \
                         kl_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * hl_290[k]
                   + f_0 * kl_290[k];

        t_291[k] = -3.0 * hl_291[k]
                   + f_0 * kl_291[k];

        t_292[k] = -3.0 * hl_292[k]
                   + f_0 * kl_292[k];

        t_293[k] = -3.0 * hl_293[k]
                   + f_0 * kl_293[k];

        t_294[k] = -3.0 * hl_294[k]
                   + f_0 * kl_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, hl_295, hl_296, hl_297, hl_298, \
                         hl_299, kl_295, kl_296, kl_297, kl_298, \
                         kl_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -3.0 * hl_295[k]
                   + f_0 * kl_295[k];

        t_296[k] = -3.0 * hl_296[k]
                   + f_0 * kl_296[k];

        t_297[k] = -3.0 * hl_297[k]
                   + f_0 * kl_297[k];

        t_298[k] = -3.0 * hl_298[k]
                   + f_0 * kl_298[k];

        t_299[k] = -3.0 * hl_299[k]
                   + f_0 * kl_299[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
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

    const auto *hl_300 = buffer.data(hl + 300);
    const auto *hl_301 = buffer.data(hl + 301);
    const auto *hl_302 = buffer.data(hl + 302);
    const auto *hl_303 = buffer.data(hl + 303);
    const auto *hl_304 = buffer.data(hl + 304);
    const auto *hl_305 = buffer.data(hl + 305);
    const auto *hl_306 = buffer.data(hl + 306);
    const auto *hl_307 = buffer.data(hl + 307);
    const auto *hl_308 = buffer.data(hl + 308);
    const auto *hl_309 = buffer.data(hl + 309);
    const auto *hl_310 = buffer.data(hl + 310);
    const auto *hl_311 = buffer.data(hl + 311);
    const auto *hl_312 = buffer.data(hl + 312);
    const auto *hl_313 = buffer.data(hl + 313);
    const auto *hl_314 = buffer.data(hl + 314);
    const auto *hl_315 = buffer.data(hl + 315);
    const auto *hl_316 = buffer.data(hl + 316);
    const auto *hl_317 = buffer.data(hl + 317);
    const auto *hl_318 = buffer.data(hl + 318);
    const auto *hl_319 = buffer.data(hl + 319);
    const auto *hl_320 = buffer.data(hl + 320);
    const auto *hl_321 = buffer.data(hl + 321);
    const auto *hl_322 = buffer.data(hl + 322);
    const auto *hl_323 = buffer.data(hl + 323);
    const auto *hl_324 = buffer.data(hl + 324);
    const auto *hl_325 = buffer.data(hl + 325);
    const auto *hl_326 = buffer.data(hl + 326);
    const auto *hl_327 = buffer.data(hl + 327);
    const auto *hl_328 = buffer.data(hl + 328);
    const auto *hl_329 = buffer.data(hl + 329);
    const auto *hl_330 = buffer.data(hl + 330);
    const auto *hl_331 = buffer.data(hl + 331);
    const auto *hl_332 = buffer.data(hl + 332);
    const auto *hl_333 = buffer.data(hl + 333);
    const auto *hl_334 = buffer.data(hl + 334);
    const auto *hl_335 = buffer.data(hl + 335);
    const auto *hl_336 = buffer.data(hl + 336);
    const auto *hl_337 = buffer.data(hl + 337);
    const auto *hl_338 = buffer.data(hl + 338);
    const auto *hl_339 = buffer.data(hl + 339);
    const auto *hl_340 = buffer.data(hl + 340);
    const auto *hl_341 = buffer.data(hl + 341);
    const auto *hl_342 = buffer.data(hl + 342);
    const auto *hl_343 = buffer.data(hl + 343);
    const auto *hl_344 = buffer.data(hl + 344);
    const auto *hl_345 = buffer.data(hl + 345);
    const auto *hl_346 = buffer.data(hl + 346);
    const auto *hl_347 = buffer.data(hl + 347);
    const auto *hl_348 = buffer.data(hl + 348);
    const auto *hl_349 = buffer.data(hl + 349);
    const auto *hl_350 = buffer.data(hl + 350);
    const auto *hl_351 = buffer.data(hl + 351);
    const auto *hl_352 = buffer.data(hl + 352);
    const auto *hl_353 = buffer.data(hl + 353);
    const auto *hl_354 = buffer.data(hl + 354);
    const auto *hl_355 = buffer.data(hl + 355);
    const auto *hl_356 = buffer.data(hl + 356);
    const auto *hl_357 = buffer.data(hl + 357);
    const auto *hl_358 = buffer.data(hl + 358);
    const auto *hl_359 = buffer.data(hl + 359);
    const auto *hl_360 = buffer.data(hl + 360);
    const auto *hl_361 = buffer.data(hl + 361);
    const auto *hl_362 = buffer.data(hl + 362);
    const auto *hl_363 = buffer.data(hl + 363);
    const auto *hl_364 = buffer.data(hl + 364);
    const auto *hl_365 = buffer.data(hl + 365);
    const auto *hl_366 = buffer.data(hl + 366);
    const auto *hl_367 = buffer.data(hl + 367);
    const auto *hl_368 = buffer.data(hl + 368);
    const auto *hl_369 = buffer.data(hl + 369);
    const auto *hl_370 = buffer.data(hl + 370);
    const auto *hl_371 = buffer.data(hl + 371);
    const auto *hl_372 = buffer.data(hl + 372);
    const auto *hl_373 = buffer.data(hl + 373);
    const auto *hl_374 = buffer.data(hl + 374);
    const auto *hl_375 = buffer.data(hl + 375);
    const auto *hl_376 = buffer.data(hl + 376);
    const auto *hl_377 = buffer.data(hl + 377);
    const auto *hl_378 = buffer.data(hl + 378);
    const auto *hl_379 = buffer.data(hl + 379);
    const auto *hl_380 = buffer.data(hl + 380);
    const auto *hl_381 = buffer.data(hl + 381);
    const auto *hl_382 = buffer.data(hl + 382);
    const auto *hl_383 = buffer.data(hl + 383);
    const auto *hl_384 = buffer.data(hl + 384);
    const auto *hl_385 = buffer.data(hl + 385);
    const auto *hl_386 = buffer.data(hl + 386);
    const auto *hl_387 = buffer.data(hl + 387);
    const auto *hl_388 = buffer.data(hl + 388);
    const auto *hl_389 = buffer.data(hl + 389);
    const auto *hl_390 = buffer.data(hl + 390);
    const auto *hl_391 = buffer.data(hl + 391);
    const auto *hl_392 = buffer.data(hl + 392);
    const auto *hl_393 = buffer.data(hl + 393);
    const auto *hl_394 = buffer.data(hl + 394);
    const auto *hl_395 = buffer.data(hl + 395);
    const auto *hl_396 = buffer.data(hl + 396);
    const auto *hl_397 = buffer.data(hl + 397);
    const auto *hl_398 = buffer.data(hl + 398);
    const auto *hl_399 = buffer.data(hl + 399);
    const auto *hl_400 = buffer.data(hl + 400);
    const auto *hl_401 = buffer.data(hl + 401);
    const auto *hl_402 = buffer.data(hl + 402);
    const auto *hl_403 = buffer.data(hl + 403);
    const auto *hl_404 = buffer.data(hl + 404);
    const auto *hl_405 = buffer.data(hl + 405);
    const auto *hl_406 = buffer.data(hl + 406);
    const auto *hl_407 = buffer.data(hl + 407);
    const auto *hl_408 = buffer.data(hl + 408);
    const auto *hl_409 = buffer.data(hl + 409);
    const auto *hl_410 = buffer.data(hl + 410);
    const auto *hl_411 = buffer.data(hl + 411);
    const auto *hl_412 = buffer.data(hl + 412);
    const auto *hl_413 = buffer.data(hl + 413);
    const auto *hl_414 = buffer.data(hl + 414);
    const auto *hl_415 = buffer.data(hl + 415);
    const auto *hl_416 = buffer.data(hl + 416);
    const auto *hl_417 = buffer.data(hl + 417);
    const auto *hl_418 = buffer.data(hl + 418);
    const auto *hl_419 = buffer.data(hl + 419);
    const auto *hl_420 = buffer.data(hl + 420);
    const auto *hl_421 = buffer.data(hl + 421);
    const auto *hl_422 = buffer.data(hl + 422);
    const auto *hl_423 = buffer.data(hl + 423);
    const auto *hl_424 = buffer.data(hl + 424);
    const auto *hl_425 = buffer.data(hl + 425);
    const auto *hl_426 = buffer.data(hl + 426);
    const auto *hl_427 = buffer.data(hl + 427);
    const auto *hl_428 = buffer.data(hl + 428);
    const auto *hl_429 = buffer.data(hl + 429);
    const auto *hl_430 = buffer.data(hl + 430);
    const auto *hl_431 = buffer.data(hl + 431);
    const auto *hl_432 = buffer.data(hl + 432);
    const auto *hl_433 = buffer.data(hl + 433);
    const auto *hl_434 = buffer.data(hl + 434);
    const auto *hl_435 = buffer.data(hl + 435);
    const auto *hl_436 = buffer.data(hl + 436);
    const auto *hl_437 = buffer.data(hl + 437);
    const auto *hl_438 = buffer.data(hl + 438);
    const auto *hl_439 = buffer.data(hl + 439);
    const auto *hl_440 = buffer.data(hl + 440);
    const auto *hl_441 = buffer.data(hl + 441);
    const auto *hl_442 = buffer.data(hl + 442);
    const auto *hl_443 = buffer.data(hl + 443);
    const auto *hl_444 = buffer.data(hl + 444);
    const auto *hl_445 = buffer.data(hl + 445);
    const auto *hl_446 = buffer.data(hl + 446);
    const auto *hl_447 = buffer.data(hl + 447);
    const auto *hl_448 = buffer.data(hl + 448);
    const auto *hl_449 = buffer.data(hl + 449);

    const auto *kl_300 = buffer.data(kl + 300);
    const auto *kl_301 = buffer.data(kl + 301);
    const auto *kl_302 = buffer.data(kl + 302);
    const auto *kl_303 = buffer.data(kl + 303);
    const auto *kl_304 = buffer.data(kl + 304);
    const auto *kl_305 = buffer.data(kl + 305);
    const auto *kl_306 = buffer.data(kl + 306);
    const auto *kl_307 = buffer.data(kl + 307);
    const auto *kl_308 = buffer.data(kl + 308);
    const auto *kl_309 = buffer.data(kl + 309);
    const auto *kl_310 = buffer.data(kl + 310);
    const auto *kl_311 = buffer.data(kl + 311);
    const auto *kl_312 = buffer.data(kl + 312);
    const auto *kl_313 = buffer.data(kl + 313);
    const auto *kl_314 = buffer.data(kl + 314);
    const auto *kl_315 = buffer.data(kl + 315);
    const auto *kl_316 = buffer.data(kl + 316);
    const auto *kl_317 = buffer.data(kl + 317);
    const auto *kl_318 = buffer.data(kl + 318);
    const auto *kl_319 = buffer.data(kl + 319);
    const auto *kl_320 = buffer.data(kl + 320);
    const auto *kl_321 = buffer.data(kl + 321);
    const auto *kl_322 = buffer.data(kl + 322);
    const auto *kl_323 = buffer.data(kl + 323);
    const auto *kl_324 = buffer.data(kl + 324);
    const auto *kl_325 = buffer.data(kl + 325);
    const auto *kl_326 = buffer.data(kl + 326);
    const auto *kl_327 = buffer.data(kl + 327);
    const auto *kl_328 = buffer.data(kl + 328);
    const auto *kl_329 = buffer.data(kl + 329);
    const auto *kl_330 = buffer.data(kl + 330);
    const auto *kl_331 = buffer.data(kl + 331);
    const auto *kl_332 = buffer.data(kl + 332);
    const auto *kl_333 = buffer.data(kl + 333);
    const auto *kl_334 = buffer.data(kl + 334);
    const auto *kl_335 = buffer.data(kl + 335);
    const auto *kl_336 = buffer.data(kl + 336);
    const auto *kl_337 = buffer.data(kl + 337);
    const auto *kl_338 = buffer.data(kl + 338);
    const auto *kl_339 = buffer.data(kl + 339);
    const auto *kl_340 = buffer.data(kl + 340);
    const auto *kl_341 = buffer.data(kl + 341);
    const auto *kl_342 = buffer.data(kl + 342);
    const auto *kl_343 = buffer.data(kl + 343);
    const auto *kl_344 = buffer.data(kl + 344);
    const auto *kl_345 = buffer.data(kl + 345);
    const auto *kl_346 = buffer.data(kl + 346);
    const auto *kl_347 = buffer.data(kl + 347);
    const auto *kl_348 = buffer.data(kl + 348);
    const auto *kl_349 = buffer.data(kl + 349);
    const auto *kl_350 = buffer.data(kl + 350);
    const auto *kl_351 = buffer.data(kl + 351);
    const auto *kl_352 = buffer.data(kl + 352);
    const auto *kl_353 = buffer.data(kl + 353);
    const auto *kl_354 = buffer.data(kl + 354);
    const auto *kl_355 = buffer.data(kl + 355);
    const auto *kl_356 = buffer.data(kl + 356);
    const auto *kl_357 = buffer.data(kl + 357);
    const auto *kl_358 = buffer.data(kl + 358);
    const auto *kl_359 = buffer.data(kl + 359);
    const auto *kl_360 = buffer.data(kl + 360);
    const auto *kl_361 = buffer.data(kl + 361);
    const auto *kl_362 = buffer.data(kl + 362);
    const auto *kl_363 = buffer.data(kl + 363);
    const auto *kl_364 = buffer.data(kl + 364);
    const auto *kl_365 = buffer.data(kl + 365);
    const auto *kl_366 = buffer.data(kl + 366);
    const auto *kl_367 = buffer.data(kl + 367);
    const auto *kl_368 = buffer.data(kl + 368);
    const auto *kl_369 = buffer.data(kl + 369);
    const auto *kl_370 = buffer.data(kl + 370);
    const auto *kl_371 = buffer.data(kl + 371);
    const auto *kl_372 = buffer.data(kl + 372);
    const auto *kl_373 = buffer.data(kl + 373);
    const auto *kl_374 = buffer.data(kl + 374);
    const auto *kl_375 = buffer.data(kl + 375);
    const auto *kl_376 = buffer.data(kl + 376);
    const auto *kl_377 = buffer.data(kl + 377);
    const auto *kl_378 = buffer.data(kl + 378);
    const auto *kl_379 = buffer.data(kl + 379);
    const auto *kl_380 = buffer.data(kl + 380);
    const auto *kl_381 = buffer.data(kl + 381);
    const auto *kl_382 = buffer.data(kl + 382);
    const auto *kl_383 = buffer.data(kl + 383);
    const auto *kl_384 = buffer.data(kl + 384);
    const auto *kl_385 = buffer.data(kl + 385);
    const auto *kl_386 = buffer.data(kl + 386);
    const auto *kl_387 = buffer.data(kl + 387);
    const auto *kl_388 = buffer.data(kl + 388);
    const auto *kl_389 = buffer.data(kl + 389);
    const auto *kl_390 = buffer.data(kl + 390);
    const auto *kl_391 = buffer.data(kl + 391);
    const auto *kl_392 = buffer.data(kl + 392);
    const auto *kl_393 = buffer.data(kl + 393);
    const auto *kl_394 = buffer.data(kl + 394);
    const auto *kl_395 = buffer.data(kl + 395);
    const auto *kl_396 = buffer.data(kl + 396);
    const auto *kl_397 = buffer.data(kl + 397);
    const auto *kl_398 = buffer.data(kl + 398);
    const auto *kl_399 = buffer.data(kl + 399);
    const auto *kl_400 = buffer.data(kl + 400);
    const auto *kl_401 = buffer.data(kl + 401);
    const auto *kl_402 = buffer.data(kl + 402);
    const auto *kl_403 = buffer.data(kl + 403);
    const auto *kl_404 = buffer.data(kl + 404);
    const auto *kl_405 = buffer.data(kl + 405);
    const auto *kl_406 = buffer.data(kl + 406);
    const auto *kl_407 = buffer.data(kl + 407);
    const auto *kl_408 = buffer.data(kl + 408);
    const auto *kl_409 = buffer.data(kl + 409);
    const auto *kl_410 = buffer.data(kl + 410);
    const auto *kl_411 = buffer.data(kl + 411);
    const auto *kl_412 = buffer.data(kl + 412);
    const auto *kl_413 = buffer.data(kl + 413);
    const auto *kl_414 = buffer.data(kl + 414);
    const auto *kl_415 = buffer.data(kl + 415);
    const auto *kl_416 = buffer.data(kl + 416);
    const auto *kl_417 = buffer.data(kl + 417);
    const auto *kl_418 = buffer.data(kl + 418);
    const auto *kl_419 = buffer.data(kl + 419);
    const auto *kl_420 = buffer.data(kl + 420);
    const auto *kl_421 = buffer.data(kl + 421);
    const auto *kl_422 = buffer.data(kl + 422);
    const auto *kl_423 = buffer.data(kl + 423);
    const auto *kl_424 = buffer.data(kl + 424);
    const auto *kl_425 = buffer.data(kl + 425);
    const auto *kl_426 = buffer.data(kl + 426);
    const auto *kl_427 = buffer.data(kl + 427);
    const auto *kl_428 = buffer.data(kl + 428);
    const auto *kl_429 = buffer.data(kl + 429);
    const auto *kl_430 = buffer.data(kl + 430);
    const auto *kl_431 = buffer.data(kl + 431);
    const auto *kl_432 = buffer.data(kl + 432);
    const auto *kl_433 = buffer.data(kl + 433);
    const auto *kl_434 = buffer.data(kl + 434);
    const auto *kl_435 = buffer.data(kl + 435);
    const auto *kl_436 = buffer.data(kl + 436);
    const auto *kl_437 = buffer.data(kl + 437);
    const auto *kl_438 = buffer.data(kl + 438);
    const auto *kl_439 = buffer.data(kl + 439);
    const auto *kl_440 = buffer.data(kl + 440);
    const auto *kl_441 = buffer.data(kl + 441);
    const auto *kl_442 = buffer.data(kl + 442);
    const auto *kl_443 = buffer.data(kl + 443);
    const auto *kl_444 = buffer.data(kl + 444);
    const auto *kl_445 = buffer.data(kl + 445);
    const auto *kl_446 = buffer.data(kl + 446);
    const auto *kl_447 = buffer.data(kl + 447);
    const auto *kl_448 = buffer.data(kl + 448);
    const auto *kl_449 = buffer.data(kl + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, hl_300, hl_301, hl_302, hl_303, \
                         hl_304, kl_300, kl_301, kl_302, kl_303, \
                         kl_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -3.0 * hl_300[k]
                   + f_0 * kl_300[k];

        t_301[k] = -3.0 * hl_301[k]
                   + f_0 * kl_301[k];

        t_302[k] = -3.0 * hl_302[k]
                   + f_0 * kl_302[k];

        t_303[k] = -3.0 * hl_303[k]
                   + f_0 * kl_303[k];

        t_304[k] = -3.0 * hl_304[k]
                   + f_0 * kl_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, hl_305, hl_306, hl_307, hl_308, \
                         hl_309, kl_305, kl_306, kl_307, kl_308, \
                         kl_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -3.0 * hl_305[k]
                   + f_0 * kl_305[k];

        t_306[k] = -3.0 * hl_306[k]
                   + f_0 * kl_306[k];

        t_307[k] = -3.0 * hl_307[k]
                   + f_0 * kl_307[k];

        t_308[k] = -3.0 * hl_308[k]
                   + f_0 * kl_308[k];

        t_309[k] = -3.0 * hl_309[k]
                   + f_0 * kl_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, hl_310, hl_311, hl_312, hl_313, \
                         hl_314, kl_310, kl_311, kl_312, kl_313, \
                         kl_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * hl_310[k]
                   + f_0 * kl_310[k];

        t_311[k] = -3.0 * hl_311[k]
                   + f_0 * kl_311[k];

        t_312[k] = -3.0 * hl_312[k]
                   + f_0 * kl_312[k];

        t_313[k] = -3.0 * hl_313[k]
                   + f_0 * kl_313[k];

        t_314[k] = -3.0 * hl_314[k]
                   + f_0 * kl_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, hl_315, hl_316, hl_317, hl_318, \
                         hl_319, kl_315, kl_316, kl_317, kl_318, \
                         kl_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -3.0 * hl_315[k]
                   + f_0 * kl_315[k];

        t_316[k] = -3.0 * hl_316[k]
                   + f_0 * kl_316[k];

        t_317[k] = -3.0 * hl_317[k]
                   + f_0 * kl_317[k];

        t_318[k] = -3.0 * hl_318[k]
                   + f_0 * kl_318[k];

        t_319[k] = -3.0 * hl_319[k]
                   + f_0 * kl_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, hl_320, hl_321, hl_322, hl_323, \
                         hl_324, kl_320, kl_321, kl_322, kl_323, \
                         kl_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -3.0 * hl_320[k]
                   + f_0 * kl_320[k];

        t_321[k] = -3.0 * hl_321[k]
                   + f_0 * kl_321[k];

        t_322[k] = -3.0 * hl_322[k]
                   + f_0 * kl_322[k];

        t_323[k] = -3.0 * hl_323[k]
                   + f_0 * kl_323[k];

        t_324[k] = -3.0 * hl_324[k]
                   + f_0 * kl_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, hl_325, hl_326, hl_327, hl_328, \
                         hl_329, kl_325, kl_326, kl_327, kl_328, \
                         kl_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -3.0 * hl_325[k]
                   + f_0 * kl_325[k];

        t_326[k] = -3.0 * hl_326[k]
                   + f_0 * kl_326[k];

        t_327[k] = -3.0 * hl_327[k]
                   + f_0 * kl_327[k];

        t_328[k] = -3.0 * hl_328[k]
                   + f_0 * kl_328[k];

        t_329[k] = -3.0 * hl_329[k]
                   + f_0 * kl_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, hl_330, hl_331, hl_332, hl_333, \
                         hl_334, kl_330, kl_331, kl_332, kl_333, \
                         kl_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -3.0 * hl_330[k]
                   + f_0 * kl_330[k];

        t_331[k] = -3.0 * hl_331[k]
                   + f_0 * kl_331[k];

        t_332[k] = -3.0 * hl_332[k]
                   + f_0 * kl_332[k];

        t_333[k] = -3.0 * hl_333[k]
                   + f_0 * kl_333[k];

        t_334[k] = -3.0 * hl_334[k]
                   + f_0 * kl_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, hl_335, hl_336, hl_337, hl_338, \
                         hl_339, kl_335, kl_336, kl_337, kl_338, \
                         kl_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -3.0 * hl_335[k]
                   + f_0 * kl_335[k];

        t_336[k] = -3.0 * hl_336[k]
                   + f_0 * kl_336[k];

        t_337[k] = -3.0 * hl_337[k]
                   + f_0 * kl_337[k];

        t_338[k] = -3.0 * hl_338[k]
                   + f_0 * kl_338[k];

        t_339[k] = -3.0 * hl_339[k]
                   + f_0 * kl_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, hl_340, hl_341, hl_342, hl_343, \
                         hl_344, kl_340, kl_341, kl_342, kl_343, \
                         kl_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -3.0 * hl_340[k]
                   + f_0 * kl_340[k];

        t_341[k] = -3.0 * hl_341[k]
                   + f_0 * kl_341[k];

        t_342[k] = -3.0 * hl_342[k]
                   + f_0 * kl_342[k];

        t_343[k] = -3.0 * hl_343[k]
                   + f_0 * kl_343[k];

        t_344[k] = -3.0 * hl_344[k]
                   + f_0 * kl_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, hl_345, hl_346, hl_347, hl_348, \
                         hl_349, kl_345, kl_346, kl_347, kl_348, \
                         kl_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -3.0 * hl_345[k]
                   + f_0 * kl_345[k];

        t_346[k] = -3.0 * hl_346[k]
                   + f_0 * kl_346[k];

        t_347[k] = -3.0 * hl_347[k]
                   + f_0 * kl_347[k];

        t_348[k] = -3.0 * hl_348[k]
                   + f_0 * kl_348[k];

        t_349[k] = -3.0 * hl_349[k]
                   + f_0 * kl_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, hl_350, hl_351, hl_352, hl_353, \
                         hl_354, kl_350, kl_351, kl_352, kl_353, \
                         kl_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -3.0 * hl_350[k]
                   + f_0 * kl_350[k];

        t_351[k] = -3.0 * hl_351[k]
                   + f_0 * kl_351[k];

        t_352[k] = -3.0 * hl_352[k]
                   + f_0 * kl_352[k];

        t_353[k] = -3.0 * hl_353[k]
                   + f_0 * kl_353[k];

        t_354[k] = -3.0 * hl_354[k]
                   + f_0 * kl_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, hl_355, hl_356, hl_357, hl_358, \
                         hl_359, kl_355, kl_356, kl_357, kl_358, \
                         kl_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -3.0 * hl_355[k]
                   + f_0 * kl_355[k];

        t_356[k] = -3.0 * hl_356[k]
                   + f_0 * kl_356[k];

        t_357[k] = -3.0 * hl_357[k]
                   + f_0 * kl_357[k];

        t_358[k] = -3.0 * hl_358[k]
                   + f_0 * kl_358[k];

        t_359[k] = -3.0 * hl_359[k]
                   + f_0 * kl_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, hl_360, hl_361, hl_362, hl_363, \
                         hl_364, kl_360, kl_361, kl_362, kl_363, \
                         kl_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -3.0 * hl_360[k]
                   + f_0 * kl_360[k];

        t_361[k] = -3.0 * hl_361[k]
                   + f_0 * kl_361[k];

        t_362[k] = -3.0 * hl_362[k]
                   + f_0 * kl_362[k];

        t_363[k] = -3.0 * hl_363[k]
                   + f_0 * kl_363[k];

        t_364[k] = -3.0 * hl_364[k]
                   + f_0 * kl_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, hl_365, hl_366, hl_367, hl_368, \
                         hl_369, kl_365, kl_366, kl_367, kl_368, \
                         kl_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -3.0 * hl_365[k]
                   + f_0 * kl_365[k];

        t_366[k] = -3.0 * hl_366[k]
                   + f_0 * kl_366[k];

        t_367[k] = -3.0 * hl_367[k]
                   + f_0 * kl_367[k];

        t_368[k] = -3.0 * hl_368[k]
                   + f_0 * kl_368[k];

        t_369[k] = -3.0 * hl_369[k]
                   + f_0 * kl_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, hl_370, hl_371, hl_372, hl_373, \
                         hl_374, kl_370, kl_371, kl_372, kl_373, \
                         kl_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -3.0 * hl_370[k]
                   + f_0 * kl_370[k];

        t_371[k] = -3.0 * hl_371[k]
                   + f_0 * kl_371[k];

        t_372[k] = -3.0 * hl_372[k]
                   + f_0 * kl_372[k];

        t_373[k] = -3.0 * hl_373[k]
                   + f_0 * kl_373[k];

        t_374[k] = -3.0 * hl_374[k]
                   + f_0 * kl_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, hl_375, hl_376, hl_377, hl_378, \
                         hl_379, kl_375, kl_376, kl_377, kl_378, \
                         kl_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -3.0 * hl_375[k]
                   + f_0 * kl_375[k];

        t_376[k] = -3.0 * hl_376[k]
                   + f_0 * kl_376[k];

        t_377[k] = -3.0 * hl_377[k]
                   + f_0 * kl_377[k];

        t_378[k] = -3.0 * hl_378[k]
                   + f_0 * kl_378[k];

        t_379[k] = -3.0 * hl_379[k]
                   + f_0 * kl_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, hl_380, hl_381, hl_382, hl_383, \
                         hl_384, kl_380, kl_381, kl_382, kl_383, \
                         kl_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -3.0 * hl_380[k]
                   + f_0 * kl_380[k];

        t_381[k] = -3.0 * hl_381[k]
                   + f_0 * kl_381[k];

        t_382[k] = -3.0 * hl_382[k]
                   + f_0 * kl_382[k];

        t_383[k] = -3.0 * hl_383[k]
                   + f_0 * kl_383[k];

        t_384[k] = -3.0 * hl_384[k]
                   + f_0 * kl_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, hl_385, hl_386, hl_387, hl_388, \
                         hl_389, kl_385, kl_386, kl_387, kl_388, \
                         kl_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -3.0 * hl_385[k]
                   + f_0 * kl_385[k];

        t_386[k] = -3.0 * hl_386[k]
                   + f_0 * kl_386[k];

        t_387[k] = -3.0 * hl_387[k]
                   + f_0 * kl_387[k];

        t_388[k] = -3.0 * hl_388[k]
                   + f_0 * kl_388[k];

        t_389[k] = -3.0 * hl_389[k]
                   + f_0 * kl_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, hl_390, hl_391, hl_392, hl_393, \
                         hl_394, kl_390, kl_391, kl_392, kl_393, \
                         kl_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -3.0 * hl_390[k]
                   + f_0 * kl_390[k];

        t_391[k] = -3.0 * hl_391[k]
                   + f_0 * kl_391[k];

        t_392[k] = -3.0 * hl_392[k]
                   + f_0 * kl_392[k];

        t_393[k] = -3.0 * hl_393[k]
                   + f_0 * kl_393[k];

        t_394[k] = -3.0 * hl_394[k]
                   + f_0 * kl_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, hl_395, hl_396, hl_397, hl_398, \
                         hl_399, kl_395, kl_396, kl_397, kl_398, \
                         kl_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -3.0 * hl_395[k]
                   + f_0 * kl_395[k];

        t_396[k] = -3.0 * hl_396[k]
                   + f_0 * kl_396[k];

        t_397[k] = -3.0 * hl_397[k]
                   + f_0 * kl_397[k];

        t_398[k] = -3.0 * hl_398[k]
                   + f_0 * kl_398[k];

        t_399[k] = -3.0 * hl_399[k]
                   + f_0 * kl_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, hl_400, hl_401, hl_402, hl_403, \
                         hl_404, kl_400, kl_401, kl_402, kl_403, \
                         kl_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -3.0 * hl_400[k]
                   + f_0 * kl_400[k];

        t_401[k] = -3.0 * hl_401[k]
                   + f_0 * kl_401[k];

        t_402[k] = -3.0 * hl_402[k]
                   + f_0 * kl_402[k];

        t_403[k] = -3.0 * hl_403[k]
                   + f_0 * kl_403[k];

        t_404[k] = -3.0 * hl_404[k]
                   + f_0 * kl_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, hl_405, hl_406, hl_407, hl_408, \
                         hl_409, kl_405, kl_406, kl_407, kl_408, \
                         kl_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -3.0 * hl_405[k]
                   + f_0 * kl_405[k];

        t_406[k] = -3.0 * hl_406[k]
                   + f_0 * kl_406[k];

        t_407[k] = -3.0 * hl_407[k]
                   + f_0 * kl_407[k];

        t_408[k] = -3.0 * hl_408[k]
                   + f_0 * kl_408[k];

        t_409[k] = -3.0 * hl_409[k]
                   + f_0 * kl_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, hl_410, hl_411, hl_412, hl_413, \
                         hl_414, kl_410, kl_411, kl_412, kl_413, \
                         kl_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -3.0 * hl_410[k]
                   + f_0 * kl_410[k];

        t_411[k] = -3.0 * hl_411[k]
                   + f_0 * kl_411[k];

        t_412[k] = -3.0 * hl_412[k]
                   + f_0 * kl_412[k];

        t_413[k] = -3.0 * hl_413[k]
                   + f_0 * kl_413[k];

        t_414[k] = -3.0 * hl_414[k]
                   + f_0 * kl_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, hl_415, hl_416, hl_417, hl_418, \
                         hl_419, kl_415, kl_416, kl_417, kl_418, \
                         kl_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -3.0 * hl_415[k]
                   + f_0 * kl_415[k];

        t_416[k] = -3.0 * hl_416[k]
                   + f_0 * kl_416[k];

        t_417[k] = -3.0 * hl_417[k]
                   + f_0 * kl_417[k];

        t_418[k] = -3.0 * hl_418[k]
                   + f_0 * kl_418[k];

        t_419[k] = -3.0 * hl_419[k]
                   + f_0 * kl_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, hl_420, hl_421, hl_422, hl_423, \
                         hl_424, kl_420, kl_421, kl_422, kl_423, \
                         kl_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -3.0 * hl_420[k]
                   + f_0 * kl_420[k];

        t_421[k] = -3.0 * hl_421[k]
                   + f_0 * kl_421[k];

        t_422[k] = -3.0 * hl_422[k]
                   + f_0 * kl_422[k];

        t_423[k] = -3.0 * hl_423[k]
                   + f_0 * kl_423[k];

        t_424[k] = -3.0 * hl_424[k]
                   + f_0 * kl_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, hl_425, hl_426, hl_427, hl_428, \
                         hl_429, kl_425, kl_426, kl_427, kl_428, \
                         kl_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -3.0 * hl_425[k]
                   + f_0 * kl_425[k];

        t_426[k] = -3.0 * hl_426[k]
                   + f_0 * kl_426[k];

        t_427[k] = -3.0 * hl_427[k]
                   + f_0 * kl_427[k];

        t_428[k] = -3.0 * hl_428[k]
                   + f_0 * kl_428[k];

        t_429[k] = -3.0 * hl_429[k]
                   + f_0 * kl_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, hl_430, hl_431, hl_432, hl_433, \
                         hl_434, kl_430, kl_431, kl_432, kl_433, \
                         kl_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -3.0 * hl_430[k]
                   + f_0 * kl_430[k];

        t_431[k] = -3.0 * hl_431[k]
                   + f_0 * kl_431[k];

        t_432[k] = -3.0 * hl_432[k]
                   + f_0 * kl_432[k];

        t_433[k] = -3.0 * hl_433[k]
                   + f_0 * kl_433[k];

        t_434[k] = -3.0 * hl_434[k]
                   + f_0 * kl_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, hl_435, hl_436, hl_437, hl_438, \
                         hl_439, kl_435, kl_436, kl_437, kl_438, \
                         kl_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -3.0 * hl_435[k]
                   + f_0 * kl_435[k];

        t_436[k] = -3.0 * hl_436[k]
                   + f_0 * kl_436[k];

        t_437[k] = -3.0 * hl_437[k]
                   + f_0 * kl_437[k];

        t_438[k] = -3.0 * hl_438[k]
                   + f_0 * kl_438[k];

        t_439[k] = -3.0 * hl_439[k]
                   + f_0 * kl_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, hl_440, hl_441, hl_442, hl_443, \
                         hl_444, kl_440, kl_441, kl_442, kl_443, \
                         kl_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -3.0 * hl_440[k]
                   + f_0 * kl_440[k];

        t_441[k] = -3.0 * hl_441[k]
                   + f_0 * kl_441[k];

        t_442[k] = -3.0 * hl_442[k]
                   + f_0 * kl_442[k];

        t_443[k] = -3.0 * hl_443[k]
                   + f_0 * kl_443[k];

        t_444[k] = -3.0 * hl_444[k]
                   + f_0 * kl_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, hl_445, hl_446, hl_447, hl_448, \
                         hl_449, kl_445, kl_446, kl_447, kl_448, \
                         kl_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -3.0 * hl_445[k]
                   + f_0 * kl_445[k];

        t_446[k] = -3.0 * hl_446[k]
                   + f_0 * kl_446[k];

        t_447[k] = -3.0 * hl_447[k]
                   + f_0 * kl_447[k];

        t_448[k] = -3.0 * hl_448[k]
                   + f_0 * kl_448[k];

        t_449[k] = -3.0 * hl_449[k]
                   + f_0 * kl_449[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
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

    const auto *hl_450 = buffer.data(hl + 450);
    const auto *hl_451 = buffer.data(hl + 451);
    const auto *hl_452 = buffer.data(hl + 452);
    const auto *hl_453 = buffer.data(hl + 453);
    const auto *hl_454 = buffer.data(hl + 454);
    const auto *hl_455 = buffer.data(hl + 455);
    const auto *hl_456 = buffer.data(hl + 456);
    const auto *hl_457 = buffer.data(hl + 457);
    const auto *hl_458 = buffer.data(hl + 458);
    const auto *hl_459 = buffer.data(hl + 459);
    const auto *hl_460 = buffer.data(hl + 460);
    const auto *hl_461 = buffer.data(hl + 461);
    const auto *hl_462 = buffer.data(hl + 462);
    const auto *hl_463 = buffer.data(hl + 463);
    const auto *hl_464 = buffer.data(hl + 464);
    const auto *hl_465 = buffer.data(hl + 465);
    const auto *hl_466 = buffer.data(hl + 466);
    const auto *hl_467 = buffer.data(hl + 467);
    const auto *hl_468 = buffer.data(hl + 468);
    const auto *hl_469 = buffer.data(hl + 469);
    const auto *hl_470 = buffer.data(hl + 470);
    const auto *hl_471 = buffer.data(hl + 471);
    const auto *hl_472 = buffer.data(hl + 472);
    const auto *hl_473 = buffer.data(hl + 473);
    const auto *hl_474 = buffer.data(hl + 474);
    const auto *hl_475 = buffer.data(hl + 475);
    const auto *hl_476 = buffer.data(hl + 476);
    const auto *hl_477 = buffer.data(hl + 477);
    const auto *hl_478 = buffer.data(hl + 478);
    const auto *hl_479 = buffer.data(hl + 479);
    const auto *hl_480 = buffer.data(hl + 480);
    const auto *hl_481 = buffer.data(hl + 481);
    const auto *hl_482 = buffer.data(hl + 482);
    const auto *hl_483 = buffer.data(hl + 483);
    const auto *hl_484 = buffer.data(hl + 484);
    const auto *hl_485 = buffer.data(hl + 485);
    const auto *hl_486 = buffer.data(hl + 486);
    const auto *hl_487 = buffer.data(hl + 487);
    const auto *hl_488 = buffer.data(hl + 488);
    const auto *hl_489 = buffer.data(hl + 489);
    const auto *hl_490 = buffer.data(hl + 490);
    const auto *hl_491 = buffer.data(hl + 491);
    const auto *hl_492 = buffer.data(hl + 492);
    const auto *hl_493 = buffer.data(hl + 493);
    const auto *hl_494 = buffer.data(hl + 494);
    const auto *hl_495 = buffer.data(hl + 495);
    const auto *hl_496 = buffer.data(hl + 496);
    const auto *hl_497 = buffer.data(hl + 497);
    const auto *hl_498 = buffer.data(hl + 498);
    const auto *hl_499 = buffer.data(hl + 499);
    const auto *hl_500 = buffer.data(hl + 500);
    const auto *hl_501 = buffer.data(hl + 501);
    const auto *hl_502 = buffer.data(hl + 502);
    const auto *hl_503 = buffer.data(hl + 503);
    const auto *hl_504 = buffer.data(hl + 504);
    const auto *hl_505 = buffer.data(hl + 505);
    const auto *hl_506 = buffer.data(hl + 506);
    const auto *hl_507 = buffer.data(hl + 507);
    const auto *hl_508 = buffer.data(hl + 508);
    const auto *hl_509 = buffer.data(hl + 509);
    const auto *hl_510 = buffer.data(hl + 510);
    const auto *hl_511 = buffer.data(hl + 511);
    const auto *hl_512 = buffer.data(hl + 512);
    const auto *hl_513 = buffer.data(hl + 513);
    const auto *hl_514 = buffer.data(hl + 514);
    const auto *hl_515 = buffer.data(hl + 515);
    const auto *hl_516 = buffer.data(hl + 516);
    const auto *hl_517 = buffer.data(hl + 517);
    const auto *hl_518 = buffer.data(hl + 518);
    const auto *hl_519 = buffer.data(hl + 519);
    const auto *hl_520 = buffer.data(hl + 520);
    const auto *hl_521 = buffer.data(hl + 521);
    const auto *hl_522 = buffer.data(hl + 522);
    const auto *hl_523 = buffer.data(hl + 523);
    const auto *hl_524 = buffer.data(hl + 524);
    const auto *hl_525 = buffer.data(hl + 525);
    const auto *hl_526 = buffer.data(hl + 526);
    const auto *hl_527 = buffer.data(hl + 527);
    const auto *hl_528 = buffer.data(hl + 528);
    const auto *hl_529 = buffer.data(hl + 529);
    const auto *hl_530 = buffer.data(hl + 530);
    const auto *hl_531 = buffer.data(hl + 531);
    const auto *hl_532 = buffer.data(hl + 532);
    const auto *hl_533 = buffer.data(hl + 533);
    const auto *hl_534 = buffer.data(hl + 534);
    const auto *hl_535 = buffer.data(hl + 535);
    const auto *hl_536 = buffer.data(hl + 536);
    const auto *hl_537 = buffer.data(hl + 537);
    const auto *hl_538 = buffer.data(hl + 538);
    const auto *hl_539 = buffer.data(hl + 539);
    const auto *hl_540 = buffer.data(hl + 540);
    const auto *hl_541 = buffer.data(hl + 541);
    const auto *hl_542 = buffer.data(hl + 542);
    const auto *hl_543 = buffer.data(hl + 543);
    const auto *hl_544 = buffer.data(hl + 544);
    const auto *hl_545 = buffer.data(hl + 545);
    const auto *hl_546 = buffer.data(hl + 546);
    const auto *hl_547 = buffer.data(hl + 547);
    const auto *hl_548 = buffer.data(hl + 548);
    const auto *hl_549 = buffer.data(hl + 549);
    const auto *hl_550 = buffer.data(hl + 550);
    const auto *hl_551 = buffer.data(hl + 551);
    const auto *hl_552 = buffer.data(hl + 552);
    const auto *hl_553 = buffer.data(hl + 553);
    const auto *hl_554 = buffer.data(hl + 554);
    const auto *hl_555 = buffer.data(hl + 555);
    const auto *hl_556 = buffer.data(hl + 556);
    const auto *hl_557 = buffer.data(hl + 557);
    const auto *hl_558 = buffer.data(hl + 558);
    const auto *hl_559 = buffer.data(hl + 559);
    const auto *hl_560 = buffer.data(hl + 560);
    const auto *hl_561 = buffer.data(hl + 561);
    const auto *hl_562 = buffer.data(hl + 562);
    const auto *hl_563 = buffer.data(hl + 563);
    const auto *hl_564 = buffer.data(hl + 564);
    const auto *hl_565 = buffer.data(hl + 565);
    const auto *hl_566 = buffer.data(hl + 566);
    const auto *hl_567 = buffer.data(hl + 567);
    const auto *hl_568 = buffer.data(hl + 568);
    const auto *hl_569 = buffer.data(hl + 569);
    const auto *hl_570 = buffer.data(hl + 570);
    const auto *hl_571 = buffer.data(hl + 571);
    const auto *hl_572 = buffer.data(hl + 572);
    const auto *hl_573 = buffer.data(hl + 573);
    const auto *hl_574 = buffer.data(hl + 574);
    const auto *hl_575 = buffer.data(hl + 575);
    const auto *hl_576 = buffer.data(hl + 576);
    const auto *hl_577 = buffer.data(hl + 577);
    const auto *hl_578 = buffer.data(hl + 578);
    const auto *hl_579 = buffer.data(hl + 579);
    const auto *hl_580 = buffer.data(hl + 580);
    const auto *hl_581 = buffer.data(hl + 581);
    const auto *hl_582 = buffer.data(hl + 582);
    const auto *hl_583 = buffer.data(hl + 583);
    const auto *hl_584 = buffer.data(hl + 584);
    const auto *hl_585 = buffer.data(hl + 585);
    const auto *hl_586 = buffer.data(hl + 586);
    const auto *hl_587 = buffer.data(hl + 587);
    const auto *hl_588 = buffer.data(hl + 588);
    const auto *hl_589 = buffer.data(hl + 589);
    const auto *hl_590 = buffer.data(hl + 590);
    const auto *hl_591 = buffer.data(hl + 591);
    const auto *hl_592 = buffer.data(hl + 592);
    const auto *hl_593 = buffer.data(hl + 593);
    const auto *hl_594 = buffer.data(hl + 594);
    const auto *hl_595 = buffer.data(hl + 595);
    const auto *hl_596 = buffer.data(hl + 596);
    const auto *hl_597 = buffer.data(hl + 597);
    const auto *hl_598 = buffer.data(hl + 598);
    const auto *hl_599 = buffer.data(hl + 599);

    const auto *kl_450 = buffer.data(kl + 450);
    const auto *kl_451 = buffer.data(kl + 451);
    const auto *kl_452 = buffer.data(kl + 452);
    const auto *kl_453 = buffer.data(kl + 453);
    const auto *kl_454 = buffer.data(kl + 454);
    const auto *kl_455 = buffer.data(kl + 455);
    const auto *kl_456 = buffer.data(kl + 456);
    const auto *kl_457 = buffer.data(kl + 457);
    const auto *kl_458 = buffer.data(kl + 458);
    const auto *kl_459 = buffer.data(kl + 459);
    const auto *kl_460 = buffer.data(kl + 460);
    const auto *kl_461 = buffer.data(kl + 461);
    const auto *kl_462 = buffer.data(kl + 462);
    const auto *kl_463 = buffer.data(kl + 463);
    const auto *kl_464 = buffer.data(kl + 464);
    const auto *kl_465 = buffer.data(kl + 465);
    const auto *kl_466 = buffer.data(kl + 466);
    const auto *kl_467 = buffer.data(kl + 467);
    const auto *kl_468 = buffer.data(kl + 468);
    const auto *kl_469 = buffer.data(kl + 469);
    const auto *kl_470 = buffer.data(kl + 470);
    const auto *kl_471 = buffer.data(kl + 471);
    const auto *kl_472 = buffer.data(kl + 472);
    const auto *kl_473 = buffer.data(kl + 473);
    const auto *kl_474 = buffer.data(kl + 474);
    const auto *kl_475 = buffer.data(kl + 475);
    const auto *kl_476 = buffer.data(kl + 476);
    const auto *kl_477 = buffer.data(kl + 477);
    const auto *kl_478 = buffer.data(kl + 478);
    const auto *kl_479 = buffer.data(kl + 479);
    const auto *kl_480 = buffer.data(kl + 480);
    const auto *kl_481 = buffer.data(kl + 481);
    const auto *kl_482 = buffer.data(kl + 482);
    const auto *kl_483 = buffer.data(kl + 483);
    const auto *kl_484 = buffer.data(kl + 484);
    const auto *kl_485 = buffer.data(kl + 485);
    const auto *kl_486 = buffer.data(kl + 486);
    const auto *kl_487 = buffer.data(kl + 487);
    const auto *kl_488 = buffer.data(kl + 488);
    const auto *kl_489 = buffer.data(kl + 489);
    const auto *kl_490 = buffer.data(kl + 490);
    const auto *kl_491 = buffer.data(kl + 491);
    const auto *kl_492 = buffer.data(kl + 492);
    const auto *kl_493 = buffer.data(kl + 493);
    const auto *kl_494 = buffer.data(kl + 494);
    const auto *kl_495 = buffer.data(kl + 495);
    const auto *kl_496 = buffer.data(kl + 496);
    const auto *kl_497 = buffer.data(kl + 497);
    const auto *kl_498 = buffer.data(kl + 498);
    const auto *kl_499 = buffer.data(kl + 499);
    const auto *kl_500 = buffer.data(kl + 500);
    const auto *kl_501 = buffer.data(kl + 501);
    const auto *kl_502 = buffer.data(kl + 502);
    const auto *kl_503 = buffer.data(kl + 503);
    const auto *kl_504 = buffer.data(kl + 504);
    const auto *kl_505 = buffer.data(kl + 505);
    const auto *kl_506 = buffer.data(kl + 506);
    const auto *kl_507 = buffer.data(kl + 507);
    const auto *kl_508 = buffer.data(kl + 508);
    const auto *kl_509 = buffer.data(kl + 509);
    const auto *kl_510 = buffer.data(kl + 510);
    const auto *kl_511 = buffer.data(kl + 511);
    const auto *kl_512 = buffer.data(kl + 512);
    const auto *kl_513 = buffer.data(kl + 513);
    const auto *kl_514 = buffer.data(kl + 514);
    const auto *kl_515 = buffer.data(kl + 515);
    const auto *kl_516 = buffer.data(kl + 516);
    const auto *kl_517 = buffer.data(kl + 517);
    const auto *kl_518 = buffer.data(kl + 518);
    const auto *kl_519 = buffer.data(kl + 519);
    const auto *kl_520 = buffer.data(kl + 520);
    const auto *kl_521 = buffer.data(kl + 521);
    const auto *kl_522 = buffer.data(kl + 522);
    const auto *kl_523 = buffer.data(kl + 523);
    const auto *kl_524 = buffer.data(kl + 524);
    const auto *kl_525 = buffer.data(kl + 525);
    const auto *kl_526 = buffer.data(kl + 526);
    const auto *kl_527 = buffer.data(kl + 527);
    const auto *kl_528 = buffer.data(kl + 528);
    const auto *kl_529 = buffer.data(kl + 529);
    const auto *kl_530 = buffer.data(kl + 530);
    const auto *kl_531 = buffer.data(kl + 531);
    const auto *kl_532 = buffer.data(kl + 532);
    const auto *kl_533 = buffer.data(kl + 533);
    const auto *kl_534 = buffer.data(kl + 534);
    const auto *kl_535 = buffer.data(kl + 535);
    const auto *kl_536 = buffer.data(kl + 536);
    const auto *kl_537 = buffer.data(kl + 537);
    const auto *kl_538 = buffer.data(kl + 538);
    const auto *kl_539 = buffer.data(kl + 539);
    const auto *kl_540 = buffer.data(kl + 540);
    const auto *kl_541 = buffer.data(kl + 541);
    const auto *kl_542 = buffer.data(kl + 542);
    const auto *kl_543 = buffer.data(kl + 543);
    const auto *kl_544 = buffer.data(kl + 544);
    const auto *kl_545 = buffer.data(kl + 545);
    const auto *kl_546 = buffer.data(kl + 546);
    const auto *kl_547 = buffer.data(kl + 547);
    const auto *kl_548 = buffer.data(kl + 548);
    const auto *kl_549 = buffer.data(kl + 549);
    const auto *kl_550 = buffer.data(kl + 550);
    const auto *kl_551 = buffer.data(kl + 551);
    const auto *kl_552 = buffer.data(kl + 552);
    const auto *kl_553 = buffer.data(kl + 553);
    const auto *kl_554 = buffer.data(kl + 554);
    const auto *kl_555 = buffer.data(kl + 555);
    const auto *kl_556 = buffer.data(kl + 556);
    const auto *kl_557 = buffer.data(kl + 557);
    const auto *kl_558 = buffer.data(kl + 558);
    const auto *kl_559 = buffer.data(kl + 559);
    const auto *kl_560 = buffer.data(kl + 560);
    const auto *kl_561 = buffer.data(kl + 561);
    const auto *kl_562 = buffer.data(kl + 562);
    const auto *kl_563 = buffer.data(kl + 563);
    const auto *kl_564 = buffer.data(kl + 564);
    const auto *kl_565 = buffer.data(kl + 565);
    const auto *kl_566 = buffer.data(kl + 566);
    const auto *kl_567 = buffer.data(kl + 567);
    const auto *kl_568 = buffer.data(kl + 568);
    const auto *kl_569 = buffer.data(kl + 569);
    const auto *kl_570 = buffer.data(kl + 570);
    const auto *kl_571 = buffer.data(kl + 571);
    const auto *kl_572 = buffer.data(kl + 572);
    const auto *kl_573 = buffer.data(kl + 573);
    const auto *kl_574 = buffer.data(kl + 574);
    const auto *kl_575 = buffer.data(kl + 575);
    const auto *kl_576 = buffer.data(kl + 576);
    const auto *kl_577 = buffer.data(kl + 577);
    const auto *kl_578 = buffer.data(kl + 578);
    const auto *kl_579 = buffer.data(kl + 579);
    const auto *kl_580 = buffer.data(kl + 580);
    const auto *kl_581 = buffer.data(kl + 581);
    const auto *kl_582 = buffer.data(kl + 582);
    const auto *kl_583 = buffer.data(kl + 583);
    const auto *kl_584 = buffer.data(kl + 584);
    const auto *kl_585 = buffer.data(kl + 585);
    const auto *kl_586 = buffer.data(kl + 586);
    const auto *kl_587 = buffer.data(kl + 587);
    const auto *kl_588 = buffer.data(kl + 588);
    const auto *kl_589 = buffer.data(kl + 589);
    const auto *kl_590 = buffer.data(kl + 590);
    const auto *kl_591 = buffer.data(kl + 591);
    const auto *kl_592 = buffer.data(kl + 592);
    const auto *kl_593 = buffer.data(kl + 593);
    const auto *kl_594 = buffer.data(kl + 594);
    const auto *kl_595 = buffer.data(kl + 595);
    const auto *kl_596 = buffer.data(kl + 596);
    const auto *kl_597 = buffer.data(kl + 597);
    const auto *kl_598 = buffer.data(kl + 598);
    const auto *kl_599 = buffer.data(kl + 599);

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, hl_450, hl_451, hl_452, hl_453, \
                         hl_454, kl_450, kl_451, kl_452, kl_453, \
                         kl_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -2.0 * hl_450[k]
                   + f_0 * kl_450[k];

        t_451[k] = -2.0 * hl_451[k]
                   + f_0 * kl_451[k];

        t_452[k] = -2.0 * hl_452[k]
                   + f_0 * kl_452[k];

        t_453[k] = -2.0 * hl_453[k]
                   + f_0 * kl_453[k];

        t_454[k] = -2.0 * hl_454[k]
                   + f_0 * kl_454[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, hl_455, hl_456, hl_457, hl_458, \
                         hl_459, kl_455, kl_456, kl_457, kl_458, \
                         kl_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -2.0 * hl_455[k]
                   + f_0 * kl_455[k];

        t_456[k] = -2.0 * hl_456[k]
                   + f_0 * kl_456[k];

        t_457[k] = -2.0 * hl_457[k]
                   + f_0 * kl_457[k];

        t_458[k] = -2.0 * hl_458[k]
                   + f_0 * kl_458[k];

        t_459[k] = -2.0 * hl_459[k]
                   + f_0 * kl_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, hl_460, hl_461, hl_462, hl_463, \
                         hl_464, kl_460, kl_461, kl_462, kl_463, \
                         kl_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -2.0 * hl_460[k]
                   + f_0 * kl_460[k];

        t_461[k] = -2.0 * hl_461[k]
                   + f_0 * kl_461[k];

        t_462[k] = -2.0 * hl_462[k]
                   + f_0 * kl_462[k];

        t_463[k] = -2.0 * hl_463[k]
                   + f_0 * kl_463[k];

        t_464[k] = -2.0 * hl_464[k]
                   + f_0 * kl_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, hl_465, hl_466, hl_467, hl_468, \
                         hl_469, kl_465, kl_466, kl_467, kl_468, \
                         kl_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -2.0 * hl_465[k]
                   + f_0 * kl_465[k];

        t_466[k] = -2.0 * hl_466[k]
                   + f_0 * kl_466[k];

        t_467[k] = -2.0 * hl_467[k]
                   + f_0 * kl_467[k];

        t_468[k] = -2.0 * hl_468[k]
                   + f_0 * kl_468[k];

        t_469[k] = -2.0 * hl_469[k]
                   + f_0 * kl_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, hl_470, hl_471, hl_472, hl_473, \
                         hl_474, kl_470, kl_471, kl_472, kl_473, \
                         kl_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -2.0 * hl_470[k]
                   + f_0 * kl_470[k];

        t_471[k] = -2.0 * hl_471[k]
                   + f_0 * kl_471[k];

        t_472[k] = -2.0 * hl_472[k]
                   + f_0 * kl_472[k];

        t_473[k] = -2.0 * hl_473[k]
                   + f_0 * kl_473[k];

        t_474[k] = -2.0 * hl_474[k]
                   + f_0 * kl_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, hl_475, hl_476, hl_477, hl_478, \
                         hl_479, kl_475, kl_476, kl_477, kl_478, \
                         kl_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -2.0 * hl_475[k]
                   + f_0 * kl_475[k];

        t_476[k] = -2.0 * hl_476[k]
                   + f_0 * kl_476[k];

        t_477[k] = -2.0 * hl_477[k]
                   + f_0 * kl_477[k];

        t_478[k] = -2.0 * hl_478[k]
                   + f_0 * kl_478[k];

        t_479[k] = -2.0 * hl_479[k]
                   + f_0 * kl_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, hl_480, hl_481, hl_482, hl_483, \
                         hl_484, kl_480, kl_481, kl_482, kl_483, \
                         kl_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -2.0 * hl_480[k]
                   + f_0 * kl_480[k];

        t_481[k] = -2.0 * hl_481[k]
                   + f_0 * kl_481[k];

        t_482[k] = -2.0 * hl_482[k]
                   + f_0 * kl_482[k];

        t_483[k] = -2.0 * hl_483[k]
                   + f_0 * kl_483[k];

        t_484[k] = -2.0 * hl_484[k]
                   + f_0 * kl_484[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, hl_485, hl_486, hl_487, hl_488, \
                         hl_489, kl_485, kl_486, kl_487, kl_488, \
                         kl_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -2.0 * hl_485[k]
                   + f_0 * kl_485[k];

        t_486[k] = -2.0 * hl_486[k]
                   + f_0 * kl_486[k];

        t_487[k] = -2.0 * hl_487[k]
                   + f_0 * kl_487[k];

        t_488[k] = -2.0 * hl_488[k]
                   + f_0 * kl_488[k];

        t_489[k] = -2.0 * hl_489[k]
                   + f_0 * kl_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, hl_490, hl_491, hl_492, hl_493, \
                         hl_494, kl_490, kl_491, kl_492, kl_493, \
                         kl_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -2.0 * hl_490[k]
                   + f_0 * kl_490[k];

        t_491[k] = -2.0 * hl_491[k]
                   + f_0 * kl_491[k];

        t_492[k] = -2.0 * hl_492[k]
                   + f_0 * kl_492[k];

        t_493[k] = -2.0 * hl_493[k]
                   + f_0 * kl_493[k];

        t_494[k] = -2.0 * hl_494[k]
                   + f_0 * kl_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, hl_495, hl_496, hl_497, hl_498, \
                         hl_499, kl_495, kl_496, kl_497, kl_498, \
                         kl_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -2.0 * hl_495[k]
                   + f_0 * kl_495[k];

        t_496[k] = -2.0 * hl_496[k]
                   + f_0 * kl_496[k];

        t_497[k] = -2.0 * hl_497[k]
                   + f_0 * kl_497[k];

        t_498[k] = -2.0 * hl_498[k]
                   + f_0 * kl_498[k];

        t_499[k] = -2.0 * hl_499[k]
                   + f_0 * kl_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, hl_500, hl_501, hl_502, hl_503, \
                         hl_504, kl_500, kl_501, kl_502, kl_503, \
                         kl_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -2.0 * hl_500[k]
                   + f_0 * kl_500[k];

        t_501[k] = -2.0 * hl_501[k]
                   + f_0 * kl_501[k];

        t_502[k] = -2.0 * hl_502[k]
                   + f_0 * kl_502[k];

        t_503[k] = -2.0 * hl_503[k]
                   + f_0 * kl_503[k];

        t_504[k] = -2.0 * hl_504[k]
                   + f_0 * kl_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, hl_505, hl_506, hl_507, hl_508, \
                         hl_509, kl_505, kl_506, kl_507, kl_508, \
                         kl_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -2.0 * hl_505[k]
                   + f_0 * kl_505[k];

        t_506[k] = -2.0 * hl_506[k]
                   + f_0 * kl_506[k];

        t_507[k] = -2.0 * hl_507[k]
                   + f_0 * kl_507[k];

        t_508[k] = -2.0 * hl_508[k]
                   + f_0 * kl_508[k];

        t_509[k] = -2.0 * hl_509[k]
                   + f_0 * kl_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, hl_510, hl_511, hl_512, hl_513, \
                         hl_514, kl_510, kl_511, kl_512, kl_513, \
                         kl_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -2.0 * hl_510[k]
                   + f_0 * kl_510[k];

        t_511[k] = -2.0 * hl_511[k]
                   + f_0 * kl_511[k];

        t_512[k] = -2.0 * hl_512[k]
                   + f_0 * kl_512[k];

        t_513[k] = -2.0 * hl_513[k]
                   + f_0 * kl_513[k];

        t_514[k] = -2.0 * hl_514[k]
                   + f_0 * kl_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, hl_515, hl_516, hl_517, hl_518, \
                         hl_519, kl_515, kl_516, kl_517, kl_518, \
                         kl_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -2.0 * hl_515[k]
                   + f_0 * kl_515[k];

        t_516[k] = -2.0 * hl_516[k]
                   + f_0 * kl_516[k];

        t_517[k] = -2.0 * hl_517[k]
                   + f_0 * kl_517[k];

        t_518[k] = -2.0 * hl_518[k]
                   + f_0 * kl_518[k];

        t_519[k] = -2.0 * hl_519[k]
                   + f_0 * kl_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, hl_520, hl_521, hl_522, hl_523, \
                         hl_524, kl_520, kl_521, kl_522, kl_523, \
                         kl_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -2.0 * hl_520[k]
                   + f_0 * kl_520[k];

        t_521[k] = -2.0 * hl_521[k]
                   + f_0 * kl_521[k];

        t_522[k] = -2.0 * hl_522[k]
                   + f_0 * kl_522[k];

        t_523[k] = -2.0 * hl_523[k]
                   + f_0 * kl_523[k];

        t_524[k] = -2.0 * hl_524[k]
                   + f_0 * kl_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, hl_525, hl_526, hl_527, hl_528, \
                         hl_529, kl_525, kl_526, kl_527, kl_528, \
                         kl_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -2.0 * hl_525[k]
                   + f_0 * kl_525[k];

        t_526[k] = -2.0 * hl_526[k]
                   + f_0 * kl_526[k];

        t_527[k] = -2.0 * hl_527[k]
                   + f_0 * kl_527[k];

        t_528[k] = -2.0 * hl_528[k]
                   + f_0 * kl_528[k];

        t_529[k] = -2.0 * hl_529[k]
                   + f_0 * kl_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, hl_530, hl_531, hl_532, hl_533, \
                         hl_534, kl_530, kl_531, kl_532, kl_533, \
                         kl_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -2.0 * hl_530[k]
                   + f_0 * kl_530[k];

        t_531[k] = -2.0 * hl_531[k]
                   + f_0 * kl_531[k];

        t_532[k] = -2.0 * hl_532[k]
                   + f_0 * kl_532[k];

        t_533[k] = -2.0 * hl_533[k]
                   + f_0 * kl_533[k];

        t_534[k] = -2.0 * hl_534[k]
                   + f_0 * kl_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, hl_535, hl_536, hl_537, hl_538, \
                         hl_539, kl_535, kl_536, kl_537, kl_538, \
                         kl_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -2.0 * hl_535[k]
                   + f_0 * kl_535[k];

        t_536[k] = -2.0 * hl_536[k]
                   + f_0 * kl_536[k];

        t_537[k] = -2.0 * hl_537[k]
                   + f_0 * kl_537[k];

        t_538[k] = -2.0 * hl_538[k]
                   + f_0 * kl_538[k];

        t_539[k] = -2.0 * hl_539[k]
                   + f_0 * kl_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, hl_540, hl_541, hl_542, hl_543, \
                         hl_544, kl_540, kl_541, kl_542, kl_543, \
                         kl_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -2.0 * hl_540[k]
                   + f_0 * kl_540[k];

        t_541[k] = -2.0 * hl_541[k]
                   + f_0 * kl_541[k];

        t_542[k] = -2.0 * hl_542[k]
                   + f_0 * kl_542[k];

        t_543[k] = -2.0 * hl_543[k]
                   + f_0 * kl_543[k];

        t_544[k] = -2.0 * hl_544[k]
                   + f_0 * kl_544[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, hl_545, hl_546, hl_547, hl_548, \
                         hl_549, kl_545, kl_546, kl_547, kl_548, \
                         kl_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -2.0 * hl_545[k]
                   + f_0 * kl_545[k];

        t_546[k] = -2.0 * hl_546[k]
                   + f_0 * kl_546[k];

        t_547[k] = -2.0 * hl_547[k]
                   + f_0 * kl_547[k];

        t_548[k] = -2.0 * hl_548[k]
                   + f_0 * kl_548[k];

        t_549[k] = -2.0 * hl_549[k]
                   + f_0 * kl_549[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, hl_550, hl_551, hl_552, hl_553, \
                         hl_554, kl_550, kl_551, kl_552, kl_553, \
                         kl_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -2.0 * hl_550[k]
                   + f_0 * kl_550[k];

        t_551[k] = -2.0 * hl_551[k]
                   + f_0 * kl_551[k];

        t_552[k] = -2.0 * hl_552[k]
                   + f_0 * kl_552[k];

        t_553[k] = -2.0 * hl_553[k]
                   + f_0 * kl_553[k];

        t_554[k] = -2.0 * hl_554[k]
                   + f_0 * kl_554[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, hl_555, hl_556, hl_557, hl_558, \
                         hl_559, kl_555, kl_556, kl_557, kl_558, \
                         kl_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -2.0 * hl_555[k]
                   + f_0 * kl_555[k];

        t_556[k] = -2.0 * hl_556[k]
                   + f_0 * kl_556[k];

        t_557[k] = -2.0 * hl_557[k]
                   + f_0 * kl_557[k];

        t_558[k] = -2.0 * hl_558[k]
                   + f_0 * kl_558[k];

        t_559[k] = -2.0 * hl_559[k]
                   + f_0 * kl_559[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, hl_560, hl_561, hl_562, hl_563, \
                         hl_564, kl_560, kl_561, kl_562, kl_563, \
                         kl_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = -2.0 * hl_560[k]
                   + f_0 * kl_560[k];

        t_561[k] = -2.0 * hl_561[k]
                   + f_0 * kl_561[k];

        t_562[k] = -2.0 * hl_562[k]
                   + f_0 * kl_562[k];

        t_563[k] = -2.0 * hl_563[k]
                   + f_0 * kl_563[k];

        t_564[k] = -2.0 * hl_564[k]
                   + f_0 * kl_564[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, hl_565, hl_566, hl_567, hl_568, \
                         hl_569, kl_565, kl_566, kl_567, kl_568, \
                         kl_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = -2.0 * hl_565[k]
                   + f_0 * kl_565[k];

        t_566[k] = -2.0 * hl_566[k]
                   + f_0 * kl_566[k];

        t_567[k] = -2.0 * hl_567[k]
                   + f_0 * kl_567[k];

        t_568[k] = -2.0 * hl_568[k]
                   + f_0 * kl_568[k];

        t_569[k] = -2.0 * hl_569[k]
                   + f_0 * kl_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, hl_570, hl_571, hl_572, hl_573, \
                         hl_574, kl_570, kl_571, kl_572, kl_573, \
                         kl_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -2.0 * hl_570[k]
                   + f_0 * kl_570[k];

        t_571[k] = -2.0 * hl_571[k]
                   + f_0 * kl_571[k];

        t_572[k] = -2.0 * hl_572[k]
                   + f_0 * kl_572[k];

        t_573[k] = -2.0 * hl_573[k]
                   + f_0 * kl_573[k];

        t_574[k] = -2.0 * hl_574[k]
                   + f_0 * kl_574[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, hl_575, hl_576, hl_577, hl_578, \
                         hl_579, kl_575, kl_576, kl_577, kl_578, \
                         kl_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = -2.0 * hl_575[k]
                   + f_0 * kl_575[k];

        t_576[k] = -2.0 * hl_576[k]
                   + f_0 * kl_576[k];

        t_577[k] = -2.0 * hl_577[k]
                   + f_0 * kl_577[k];

        t_578[k] = -2.0 * hl_578[k]
                   + f_0 * kl_578[k];

        t_579[k] = -2.0 * hl_579[k]
                   + f_0 * kl_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, hl_580, hl_581, hl_582, hl_583, \
                         hl_584, kl_580, kl_581, kl_582, kl_583, \
                         kl_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -2.0 * hl_580[k]
                   + f_0 * kl_580[k];

        t_581[k] = -2.0 * hl_581[k]
                   + f_0 * kl_581[k];

        t_582[k] = -2.0 * hl_582[k]
                   + f_0 * kl_582[k];

        t_583[k] = -2.0 * hl_583[k]
                   + f_0 * kl_583[k];

        t_584[k] = -2.0 * hl_584[k]
                   + f_0 * kl_584[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, hl_585, hl_586, hl_587, hl_588, \
                         hl_589, kl_585, kl_586, kl_587, kl_588, \
                         kl_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -2.0 * hl_585[k]
                   + f_0 * kl_585[k];

        t_586[k] = -2.0 * hl_586[k]
                   + f_0 * kl_586[k];

        t_587[k] = -2.0 * hl_587[k]
                   + f_0 * kl_587[k];

        t_588[k] = -2.0 * hl_588[k]
                   + f_0 * kl_588[k];

        t_589[k] = -2.0 * hl_589[k]
                   + f_0 * kl_589[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, hl_590, hl_591, hl_592, hl_593, \
                         hl_594, kl_590, kl_591, kl_592, kl_593, \
                         kl_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -2.0 * hl_590[k]
                   + f_0 * kl_590[k];

        t_591[k] = -2.0 * hl_591[k]
                   + f_0 * kl_591[k];

        t_592[k] = -2.0 * hl_592[k]
                   + f_0 * kl_592[k];

        t_593[k] = -2.0 * hl_593[k]
                   + f_0 * kl_593[k];

        t_594[k] = -2.0 * hl_594[k]
                   + f_0 * kl_594[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, hl_595, hl_596, hl_597, hl_598, \
                         hl_599, kl_595, kl_596, kl_597, kl_598, \
                         kl_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -2.0 * hl_595[k]
                   + f_0 * kl_595[k];

        t_596[k] = -2.0 * hl_596[k]
                   + f_0 * kl_596[k];

        t_597[k] = -2.0 * hl_597[k]
                   + f_0 * kl_597[k];

        t_598[k] = -2.0 * hl_598[k]
                   + f_0 * kl_598[k];

        t_599[k] = -2.0 * hl_599[k]
                   + f_0 * kl_599[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
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

    const auto *hl_600 = buffer.data(hl + 600);
    const auto *hl_601 = buffer.data(hl + 601);
    const auto *hl_602 = buffer.data(hl + 602);
    const auto *hl_603 = buffer.data(hl + 603);
    const auto *hl_604 = buffer.data(hl + 604);
    const auto *hl_605 = buffer.data(hl + 605);
    const auto *hl_606 = buffer.data(hl + 606);
    const auto *hl_607 = buffer.data(hl + 607);
    const auto *hl_608 = buffer.data(hl + 608);
    const auto *hl_609 = buffer.data(hl + 609);
    const auto *hl_610 = buffer.data(hl + 610);
    const auto *hl_611 = buffer.data(hl + 611);
    const auto *hl_612 = buffer.data(hl + 612);
    const auto *hl_613 = buffer.data(hl + 613);
    const auto *hl_614 = buffer.data(hl + 614);
    const auto *hl_615 = buffer.data(hl + 615);
    const auto *hl_616 = buffer.data(hl + 616);
    const auto *hl_617 = buffer.data(hl + 617);
    const auto *hl_618 = buffer.data(hl + 618);
    const auto *hl_619 = buffer.data(hl + 619);
    const auto *hl_620 = buffer.data(hl + 620);
    const auto *hl_621 = buffer.data(hl + 621);
    const auto *hl_622 = buffer.data(hl + 622);
    const auto *hl_623 = buffer.data(hl + 623);
    const auto *hl_624 = buffer.data(hl + 624);
    const auto *hl_625 = buffer.data(hl + 625);
    const auto *hl_626 = buffer.data(hl + 626);
    const auto *hl_627 = buffer.data(hl + 627);
    const auto *hl_628 = buffer.data(hl + 628);
    const auto *hl_629 = buffer.data(hl + 629);
    const auto *hl_630 = buffer.data(hl + 630);
    const auto *hl_631 = buffer.data(hl + 631);
    const auto *hl_632 = buffer.data(hl + 632);
    const auto *hl_633 = buffer.data(hl + 633);
    const auto *hl_634 = buffer.data(hl + 634);
    const auto *hl_635 = buffer.data(hl + 635);
    const auto *hl_636 = buffer.data(hl + 636);
    const auto *hl_637 = buffer.data(hl + 637);
    const auto *hl_638 = buffer.data(hl + 638);
    const auto *hl_639 = buffer.data(hl + 639);
    const auto *hl_640 = buffer.data(hl + 640);
    const auto *hl_641 = buffer.data(hl + 641);
    const auto *hl_642 = buffer.data(hl + 642);
    const auto *hl_643 = buffer.data(hl + 643);
    const auto *hl_644 = buffer.data(hl + 644);
    const auto *hl_645 = buffer.data(hl + 645);
    const auto *hl_646 = buffer.data(hl + 646);
    const auto *hl_647 = buffer.data(hl + 647);
    const auto *hl_648 = buffer.data(hl + 648);
    const auto *hl_649 = buffer.data(hl + 649);
    const auto *hl_650 = buffer.data(hl + 650);
    const auto *hl_651 = buffer.data(hl + 651);
    const auto *hl_652 = buffer.data(hl + 652);
    const auto *hl_653 = buffer.data(hl + 653);
    const auto *hl_654 = buffer.data(hl + 654);
    const auto *hl_655 = buffer.data(hl + 655);
    const auto *hl_656 = buffer.data(hl + 656);
    const auto *hl_657 = buffer.data(hl + 657);
    const auto *hl_658 = buffer.data(hl + 658);
    const auto *hl_659 = buffer.data(hl + 659);
    const auto *hl_660 = buffer.data(hl + 660);
    const auto *hl_661 = buffer.data(hl + 661);
    const auto *hl_662 = buffer.data(hl + 662);
    const auto *hl_663 = buffer.data(hl + 663);
    const auto *hl_664 = buffer.data(hl + 664);
    const auto *hl_665 = buffer.data(hl + 665);
    const auto *hl_666 = buffer.data(hl + 666);
    const auto *hl_667 = buffer.data(hl + 667);
    const auto *hl_668 = buffer.data(hl + 668);
    const auto *hl_669 = buffer.data(hl + 669);
    const auto *hl_670 = buffer.data(hl + 670);
    const auto *hl_671 = buffer.data(hl + 671);
    const auto *hl_672 = buffer.data(hl + 672);
    const auto *hl_673 = buffer.data(hl + 673);
    const auto *hl_674 = buffer.data(hl + 674);
    const auto *hl_675 = buffer.data(hl + 675);
    const auto *hl_676 = buffer.data(hl + 676);
    const auto *hl_677 = buffer.data(hl + 677);
    const auto *hl_678 = buffer.data(hl + 678);
    const auto *hl_679 = buffer.data(hl + 679);
    const auto *hl_680 = buffer.data(hl + 680);
    const auto *hl_681 = buffer.data(hl + 681);
    const auto *hl_682 = buffer.data(hl + 682);
    const auto *hl_683 = buffer.data(hl + 683);
    const auto *hl_684 = buffer.data(hl + 684);
    const auto *hl_685 = buffer.data(hl + 685);
    const auto *hl_686 = buffer.data(hl + 686);
    const auto *hl_687 = buffer.data(hl + 687);
    const auto *hl_688 = buffer.data(hl + 688);
    const auto *hl_689 = buffer.data(hl + 689);
    const auto *hl_690 = buffer.data(hl + 690);
    const auto *hl_691 = buffer.data(hl + 691);
    const auto *hl_692 = buffer.data(hl + 692);
    const auto *hl_693 = buffer.data(hl + 693);
    const auto *hl_694 = buffer.data(hl + 694);
    const auto *hl_695 = buffer.data(hl + 695);
    const auto *hl_696 = buffer.data(hl + 696);
    const auto *hl_697 = buffer.data(hl + 697);
    const auto *hl_698 = buffer.data(hl + 698);
    const auto *hl_699 = buffer.data(hl + 699);
    const auto *hl_700 = buffer.data(hl + 700);
    const auto *hl_701 = buffer.data(hl + 701);
    const auto *hl_702 = buffer.data(hl + 702);
    const auto *hl_703 = buffer.data(hl + 703);
    const auto *hl_704 = buffer.data(hl + 704);
    const auto *hl_705 = buffer.data(hl + 705);
    const auto *hl_706 = buffer.data(hl + 706);
    const auto *hl_707 = buffer.data(hl + 707);
    const auto *hl_708 = buffer.data(hl + 708);
    const auto *hl_709 = buffer.data(hl + 709);
    const auto *hl_710 = buffer.data(hl + 710);
    const auto *hl_711 = buffer.data(hl + 711);
    const auto *hl_712 = buffer.data(hl + 712);
    const auto *hl_713 = buffer.data(hl + 713);
    const auto *hl_714 = buffer.data(hl + 714);
    const auto *hl_715 = buffer.data(hl + 715);
    const auto *hl_716 = buffer.data(hl + 716);
    const auto *hl_717 = buffer.data(hl + 717);
    const auto *hl_718 = buffer.data(hl + 718);
    const auto *hl_719 = buffer.data(hl + 719);
    const auto *hl_720 = buffer.data(hl + 720);
    const auto *hl_721 = buffer.data(hl + 721);
    const auto *hl_722 = buffer.data(hl + 722);
    const auto *hl_723 = buffer.data(hl + 723);
    const auto *hl_724 = buffer.data(hl + 724);
    const auto *hl_725 = buffer.data(hl + 725);
    const auto *hl_726 = buffer.data(hl + 726);
    const auto *hl_727 = buffer.data(hl + 727);
    const auto *hl_728 = buffer.data(hl + 728);
    const auto *hl_729 = buffer.data(hl + 729);
    const auto *hl_730 = buffer.data(hl + 730);
    const auto *hl_731 = buffer.data(hl + 731);
    const auto *hl_732 = buffer.data(hl + 732);
    const auto *hl_733 = buffer.data(hl + 733);
    const auto *hl_734 = buffer.data(hl + 734);
    const auto *hl_735 = buffer.data(hl + 735);
    const auto *hl_736 = buffer.data(hl + 736);
    const auto *hl_737 = buffer.data(hl + 737);
    const auto *hl_738 = buffer.data(hl + 738);
    const auto *hl_739 = buffer.data(hl + 739);
    const auto *hl_740 = buffer.data(hl + 740);
    const auto *hl_741 = buffer.data(hl + 741);
    const auto *hl_742 = buffer.data(hl + 742);
    const auto *hl_743 = buffer.data(hl + 743);
    const auto *hl_744 = buffer.data(hl + 744);
    const auto *hl_745 = buffer.data(hl + 745);
    const auto *hl_746 = buffer.data(hl + 746);
    const auto *hl_747 = buffer.data(hl + 747);
    const auto *hl_748 = buffer.data(hl + 748);
    const auto *hl_749 = buffer.data(hl + 749);

    const auto *kl_600 = buffer.data(kl + 600);
    const auto *kl_601 = buffer.data(kl + 601);
    const auto *kl_602 = buffer.data(kl + 602);
    const auto *kl_603 = buffer.data(kl + 603);
    const auto *kl_604 = buffer.data(kl + 604);
    const auto *kl_605 = buffer.data(kl + 605);
    const auto *kl_606 = buffer.data(kl + 606);
    const auto *kl_607 = buffer.data(kl + 607);
    const auto *kl_608 = buffer.data(kl + 608);
    const auto *kl_609 = buffer.data(kl + 609);
    const auto *kl_610 = buffer.data(kl + 610);
    const auto *kl_611 = buffer.data(kl + 611);
    const auto *kl_612 = buffer.data(kl + 612);
    const auto *kl_613 = buffer.data(kl + 613);
    const auto *kl_614 = buffer.data(kl + 614);
    const auto *kl_615 = buffer.data(kl + 615);
    const auto *kl_616 = buffer.data(kl + 616);
    const auto *kl_617 = buffer.data(kl + 617);
    const auto *kl_618 = buffer.data(kl + 618);
    const auto *kl_619 = buffer.data(kl + 619);
    const auto *kl_620 = buffer.data(kl + 620);
    const auto *kl_621 = buffer.data(kl + 621);
    const auto *kl_622 = buffer.data(kl + 622);
    const auto *kl_623 = buffer.data(kl + 623);
    const auto *kl_624 = buffer.data(kl + 624);
    const auto *kl_625 = buffer.data(kl + 625);
    const auto *kl_626 = buffer.data(kl + 626);
    const auto *kl_627 = buffer.data(kl + 627);
    const auto *kl_628 = buffer.data(kl + 628);
    const auto *kl_629 = buffer.data(kl + 629);
    const auto *kl_630 = buffer.data(kl + 630);
    const auto *kl_631 = buffer.data(kl + 631);
    const auto *kl_632 = buffer.data(kl + 632);
    const auto *kl_633 = buffer.data(kl + 633);
    const auto *kl_634 = buffer.data(kl + 634);
    const auto *kl_635 = buffer.data(kl + 635);
    const auto *kl_636 = buffer.data(kl + 636);
    const auto *kl_637 = buffer.data(kl + 637);
    const auto *kl_638 = buffer.data(kl + 638);
    const auto *kl_639 = buffer.data(kl + 639);
    const auto *kl_640 = buffer.data(kl + 640);
    const auto *kl_641 = buffer.data(kl + 641);
    const auto *kl_642 = buffer.data(kl + 642);
    const auto *kl_643 = buffer.data(kl + 643);
    const auto *kl_644 = buffer.data(kl + 644);
    const auto *kl_645 = buffer.data(kl + 645);
    const auto *kl_646 = buffer.data(kl + 646);
    const auto *kl_647 = buffer.data(kl + 647);
    const auto *kl_648 = buffer.data(kl + 648);
    const auto *kl_649 = buffer.data(kl + 649);
    const auto *kl_650 = buffer.data(kl + 650);
    const auto *kl_651 = buffer.data(kl + 651);
    const auto *kl_652 = buffer.data(kl + 652);
    const auto *kl_653 = buffer.data(kl + 653);
    const auto *kl_654 = buffer.data(kl + 654);
    const auto *kl_655 = buffer.data(kl + 655);
    const auto *kl_656 = buffer.data(kl + 656);
    const auto *kl_657 = buffer.data(kl + 657);
    const auto *kl_658 = buffer.data(kl + 658);
    const auto *kl_659 = buffer.data(kl + 659);
    const auto *kl_660 = buffer.data(kl + 660);
    const auto *kl_661 = buffer.data(kl + 661);
    const auto *kl_662 = buffer.data(kl + 662);
    const auto *kl_663 = buffer.data(kl + 663);
    const auto *kl_664 = buffer.data(kl + 664);
    const auto *kl_665 = buffer.data(kl + 665);
    const auto *kl_666 = buffer.data(kl + 666);
    const auto *kl_667 = buffer.data(kl + 667);
    const auto *kl_668 = buffer.data(kl + 668);
    const auto *kl_669 = buffer.data(kl + 669);
    const auto *kl_670 = buffer.data(kl + 670);
    const auto *kl_671 = buffer.data(kl + 671);
    const auto *kl_672 = buffer.data(kl + 672);
    const auto *kl_673 = buffer.data(kl + 673);
    const auto *kl_674 = buffer.data(kl + 674);
    const auto *kl_675 = buffer.data(kl + 675);
    const auto *kl_676 = buffer.data(kl + 676);
    const auto *kl_677 = buffer.data(kl + 677);
    const auto *kl_678 = buffer.data(kl + 678);
    const auto *kl_679 = buffer.data(kl + 679);
    const auto *kl_680 = buffer.data(kl + 680);
    const auto *kl_681 = buffer.data(kl + 681);
    const auto *kl_682 = buffer.data(kl + 682);
    const auto *kl_683 = buffer.data(kl + 683);
    const auto *kl_684 = buffer.data(kl + 684);
    const auto *kl_685 = buffer.data(kl + 685);
    const auto *kl_686 = buffer.data(kl + 686);
    const auto *kl_687 = buffer.data(kl + 687);
    const auto *kl_688 = buffer.data(kl + 688);
    const auto *kl_689 = buffer.data(kl + 689);
    const auto *kl_690 = buffer.data(kl + 690);
    const auto *kl_691 = buffer.data(kl + 691);
    const auto *kl_692 = buffer.data(kl + 692);
    const auto *kl_693 = buffer.data(kl + 693);
    const auto *kl_694 = buffer.data(kl + 694);
    const auto *kl_695 = buffer.data(kl + 695);
    const auto *kl_696 = buffer.data(kl + 696);
    const auto *kl_697 = buffer.data(kl + 697);
    const auto *kl_698 = buffer.data(kl + 698);
    const auto *kl_699 = buffer.data(kl + 699);
    const auto *kl_700 = buffer.data(kl + 700);
    const auto *kl_701 = buffer.data(kl + 701);
    const auto *kl_702 = buffer.data(kl + 702);
    const auto *kl_703 = buffer.data(kl + 703);
    const auto *kl_704 = buffer.data(kl + 704);
    const auto *kl_705 = buffer.data(kl + 705);
    const auto *kl_706 = buffer.data(kl + 706);
    const auto *kl_707 = buffer.data(kl + 707);
    const auto *kl_708 = buffer.data(kl + 708);
    const auto *kl_709 = buffer.data(kl + 709);
    const auto *kl_710 = buffer.data(kl + 710);
    const auto *kl_711 = buffer.data(kl + 711);
    const auto *kl_712 = buffer.data(kl + 712);
    const auto *kl_713 = buffer.data(kl + 713);
    const auto *kl_714 = buffer.data(kl + 714);
    const auto *kl_715 = buffer.data(kl + 715);
    const auto *kl_716 = buffer.data(kl + 716);
    const auto *kl_717 = buffer.data(kl + 717);
    const auto *kl_718 = buffer.data(kl + 718);
    const auto *kl_719 = buffer.data(kl + 719);
    const auto *kl_720 = buffer.data(kl + 720);
    const auto *kl_721 = buffer.data(kl + 721);
    const auto *kl_722 = buffer.data(kl + 722);
    const auto *kl_723 = buffer.data(kl + 723);
    const auto *kl_724 = buffer.data(kl + 724);
    const auto *kl_725 = buffer.data(kl + 725);
    const auto *kl_726 = buffer.data(kl + 726);
    const auto *kl_727 = buffer.data(kl + 727);
    const auto *kl_728 = buffer.data(kl + 728);
    const auto *kl_729 = buffer.data(kl + 729);
    const auto *kl_730 = buffer.data(kl + 730);
    const auto *kl_731 = buffer.data(kl + 731);
    const auto *kl_732 = buffer.data(kl + 732);
    const auto *kl_733 = buffer.data(kl + 733);
    const auto *kl_734 = buffer.data(kl + 734);
    const auto *kl_735 = buffer.data(kl + 735);
    const auto *kl_736 = buffer.data(kl + 736);
    const auto *kl_737 = buffer.data(kl + 737);
    const auto *kl_738 = buffer.data(kl + 738);
    const auto *kl_739 = buffer.data(kl + 739);
    const auto *kl_740 = buffer.data(kl + 740);
    const auto *kl_741 = buffer.data(kl + 741);
    const auto *kl_742 = buffer.data(kl + 742);
    const auto *kl_743 = buffer.data(kl + 743);
    const auto *kl_744 = buffer.data(kl + 744);
    const auto *kl_745 = buffer.data(kl + 745);
    const auto *kl_746 = buffer.data(kl + 746);
    const auto *kl_747 = buffer.data(kl + 747);
    const auto *kl_748 = buffer.data(kl + 748);
    const auto *kl_749 = buffer.data(kl + 749);

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, hl_600, hl_601, hl_602, hl_603, \
                         hl_604, kl_600, kl_601, kl_602, kl_603, \
                         kl_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -2.0 * hl_600[k]
                   + f_0 * kl_600[k];

        t_601[k] = -2.0 * hl_601[k]
                   + f_0 * kl_601[k];

        t_602[k] = -2.0 * hl_602[k]
                   + f_0 * kl_602[k];

        t_603[k] = -2.0 * hl_603[k]
                   + f_0 * kl_603[k];

        t_604[k] = -2.0 * hl_604[k]
                   + f_0 * kl_604[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, hl_605, hl_606, hl_607, hl_608, \
                         hl_609, kl_605, kl_606, kl_607, kl_608, \
                         kl_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -2.0 * hl_605[k]
                   + f_0 * kl_605[k];

        t_606[k] = -2.0 * hl_606[k]
                   + f_0 * kl_606[k];

        t_607[k] = -2.0 * hl_607[k]
                   + f_0 * kl_607[k];

        t_608[k] = -2.0 * hl_608[k]
                   + f_0 * kl_608[k];

        t_609[k] = -2.0 * hl_609[k]
                   + f_0 * kl_609[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, hl_610, hl_611, hl_612, hl_613, \
                         hl_614, kl_610, kl_611, kl_612, kl_613, \
                         kl_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -2.0 * hl_610[k]
                   + f_0 * kl_610[k];

        t_611[k] = -2.0 * hl_611[k]
                   + f_0 * kl_611[k];

        t_612[k] = -2.0 * hl_612[k]
                   + f_0 * kl_612[k];

        t_613[k] = -2.0 * hl_613[k]
                   + f_0 * kl_613[k];

        t_614[k] = -2.0 * hl_614[k]
                   + f_0 * kl_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, hl_615, hl_616, hl_617, hl_618, \
                         hl_619, kl_615, kl_616, kl_617, kl_618, \
                         kl_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -2.0 * hl_615[k]
                   + f_0 * kl_615[k];

        t_616[k] = -2.0 * hl_616[k]
                   + f_0 * kl_616[k];

        t_617[k] = -2.0 * hl_617[k]
                   + f_0 * kl_617[k];

        t_618[k] = -2.0 * hl_618[k]
                   + f_0 * kl_618[k];

        t_619[k] = -2.0 * hl_619[k]
                   + f_0 * kl_619[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, hl_620, hl_621, hl_622, hl_623, \
                         hl_624, kl_620, kl_621, kl_622, kl_623, \
                         kl_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -2.0 * hl_620[k]
                   + f_0 * kl_620[k];

        t_621[k] = -2.0 * hl_621[k]
                   + f_0 * kl_621[k];

        t_622[k] = -2.0 * hl_622[k]
                   + f_0 * kl_622[k];

        t_623[k] = -2.0 * hl_623[k]
                   + f_0 * kl_623[k];

        t_624[k] = -2.0 * hl_624[k]
                   + f_0 * kl_624[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, hl_625, hl_626, hl_627, hl_628, \
                         hl_629, kl_625, kl_626, kl_627, kl_628, \
                         kl_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -2.0 * hl_625[k]
                   + f_0 * kl_625[k];

        t_626[k] = -2.0 * hl_626[k]
                   + f_0 * kl_626[k];

        t_627[k] = -2.0 * hl_627[k]
                   + f_0 * kl_627[k];

        t_628[k] = -2.0 * hl_628[k]
                   + f_0 * kl_628[k];

        t_629[k] = -2.0 * hl_629[k]
                   + f_0 * kl_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, hl_630, hl_631, hl_632, hl_633, \
                         hl_634, kl_630, kl_631, kl_632, kl_633, \
                         kl_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -2.0 * hl_630[k]
                   + f_0 * kl_630[k];

        t_631[k] = -2.0 * hl_631[k]
                   + f_0 * kl_631[k];

        t_632[k] = -2.0 * hl_632[k]
                   + f_0 * kl_632[k];

        t_633[k] = -2.0 * hl_633[k]
                   + f_0 * kl_633[k];

        t_634[k] = -2.0 * hl_634[k]
                   + f_0 * kl_634[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, hl_635, hl_636, hl_637, hl_638, \
                         hl_639, kl_635, kl_636, kl_637, kl_638, \
                         kl_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -2.0 * hl_635[k]
                   + f_0 * kl_635[k];

        t_636[k] = -2.0 * hl_636[k]
                   + f_0 * kl_636[k];

        t_637[k] = -2.0 * hl_637[k]
                   + f_0 * kl_637[k];

        t_638[k] = -2.0 * hl_638[k]
                   + f_0 * kl_638[k];

        t_639[k] = -2.0 * hl_639[k]
                   + f_0 * kl_639[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, hl_640, hl_641, hl_642, hl_643, \
                         hl_644, kl_640, kl_641, kl_642, kl_643, \
                         kl_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -2.0 * hl_640[k]
                   + f_0 * kl_640[k];

        t_641[k] = -2.0 * hl_641[k]
                   + f_0 * kl_641[k];

        t_642[k] = -2.0 * hl_642[k]
                   + f_0 * kl_642[k];

        t_643[k] = -2.0 * hl_643[k]
                   + f_0 * kl_643[k];

        t_644[k] = -2.0 * hl_644[k]
                   + f_0 * kl_644[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, hl_645, hl_646, hl_647, hl_648, \
                         hl_649, kl_645, kl_646, kl_647, kl_648, \
                         kl_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -2.0 * hl_645[k]
                   + f_0 * kl_645[k];

        t_646[k] = -2.0 * hl_646[k]
                   + f_0 * kl_646[k];

        t_647[k] = -2.0 * hl_647[k]
                   + f_0 * kl_647[k];

        t_648[k] = -2.0 * hl_648[k]
                   + f_0 * kl_648[k];

        t_649[k] = -2.0 * hl_649[k]
                   + f_0 * kl_649[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, hl_650, hl_651, hl_652, hl_653, \
                         hl_654, kl_650, kl_651, kl_652, kl_653, \
                         kl_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -2.0 * hl_650[k]
                   + f_0 * kl_650[k];

        t_651[k] = -2.0 * hl_651[k]
                   + f_0 * kl_651[k];

        t_652[k] = -2.0 * hl_652[k]
                   + f_0 * kl_652[k];

        t_653[k] = -2.0 * hl_653[k]
                   + f_0 * kl_653[k];

        t_654[k] = -2.0 * hl_654[k]
                   + f_0 * kl_654[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, hl_655, hl_656, hl_657, hl_658, \
                         hl_659, kl_655, kl_656, kl_657, kl_658, \
                         kl_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -2.0 * hl_655[k]
                   + f_0 * kl_655[k];

        t_656[k] = -2.0 * hl_656[k]
                   + f_0 * kl_656[k];

        t_657[k] = -2.0 * hl_657[k]
                   + f_0 * kl_657[k];

        t_658[k] = -2.0 * hl_658[k]
                   + f_0 * kl_658[k];

        t_659[k] = -2.0 * hl_659[k]
                   + f_0 * kl_659[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, hl_660, hl_661, hl_662, hl_663, \
                         hl_664, kl_660, kl_661, kl_662, kl_663, \
                         kl_664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -2.0 * hl_660[k]
                   + f_0 * kl_660[k];

        t_661[k] = -2.0 * hl_661[k]
                   + f_0 * kl_661[k];

        t_662[k] = -2.0 * hl_662[k]
                   + f_0 * kl_662[k];

        t_663[k] = -2.0 * hl_663[k]
                   + f_0 * kl_663[k];

        t_664[k] = -2.0 * hl_664[k]
                   + f_0 * kl_664[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, hl_665, hl_666, hl_667, hl_668, \
                         hl_669, kl_665, kl_666, kl_667, kl_668, \
                         kl_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -2.0 * hl_665[k]
                   + f_0 * kl_665[k];

        t_666[k] = -2.0 * hl_666[k]
                   + f_0 * kl_666[k];

        t_667[k] = -2.0 * hl_667[k]
                   + f_0 * kl_667[k];

        t_668[k] = -2.0 * hl_668[k]
                   + f_0 * kl_668[k];

        t_669[k] = -2.0 * hl_669[k]
                   + f_0 * kl_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, hl_670, hl_671, hl_672, hl_673, \
                         hl_674, kl_670, kl_671, kl_672, kl_673, \
                         kl_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -2.0 * hl_670[k]
                   + f_0 * kl_670[k];

        t_671[k] = -2.0 * hl_671[k]
                   + f_0 * kl_671[k];

        t_672[k] = -2.0 * hl_672[k]
                   + f_0 * kl_672[k];

        t_673[k] = -2.0 * hl_673[k]
                   + f_0 * kl_673[k];

        t_674[k] = -2.0 * hl_674[k]
                   + f_0 * kl_674[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, hl_675, hl_676, hl_677, hl_678, \
                         hl_679, kl_675, kl_676, kl_677, kl_678, \
                         kl_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -hl_675[k]
                   + f_0 * kl_675[k];

        t_676[k] = -hl_676[k]
                   + f_0 * kl_676[k];

        t_677[k] = -hl_677[k]
                   + f_0 * kl_677[k];

        t_678[k] = -hl_678[k]
                   + f_0 * kl_678[k];

        t_679[k] = -hl_679[k]
                   + f_0 * kl_679[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, hl_680, hl_681, hl_682, hl_683, \
                         hl_684, kl_680, kl_681, kl_682, kl_683, \
                         kl_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -hl_680[k]
                   + f_0 * kl_680[k];

        t_681[k] = -hl_681[k]
                   + f_0 * kl_681[k];

        t_682[k] = -hl_682[k]
                   + f_0 * kl_682[k];

        t_683[k] = -hl_683[k]
                   + f_0 * kl_683[k];

        t_684[k] = -hl_684[k]
                   + f_0 * kl_684[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, hl_685, hl_686, hl_687, hl_688, \
                         hl_689, kl_685, kl_686, kl_687, kl_688, \
                         kl_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -hl_685[k]
                   + f_0 * kl_685[k];

        t_686[k] = -hl_686[k]
                   + f_0 * kl_686[k];

        t_687[k] = -hl_687[k]
                   + f_0 * kl_687[k];

        t_688[k] = -hl_688[k]
                   + f_0 * kl_688[k];

        t_689[k] = -hl_689[k]
                   + f_0 * kl_689[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, hl_690, hl_691, hl_692, hl_693, \
                         hl_694, kl_690, kl_691, kl_692, kl_693, \
                         kl_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -hl_690[k]
                   + f_0 * kl_690[k];

        t_691[k] = -hl_691[k]
                   + f_0 * kl_691[k];

        t_692[k] = -hl_692[k]
                   + f_0 * kl_692[k];

        t_693[k] = -hl_693[k]
                   + f_0 * kl_693[k];

        t_694[k] = -hl_694[k]
                   + f_0 * kl_694[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, hl_695, hl_696, hl_697, hl_698, \
                         hl_699, kl_695, kl_696, kl_697, kl_698, \
                         kl_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -hl_695[k]
                   + f_0 * kl_695[k];

        t_696[k] = -hl_696[k]
                   + f_0 * kl_696[k];

        t_697[k] = -hl_697[k]
                   + f_0 * kl_697[k];

        t_698[k] = -hl_698[k]
                   + f_0 * kl_698[k];

        t_699[k] = -hl_699[k]
                   + f_0 * kl_699[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, hl_700, hl_701, hl_702, hl_703, \
                         hl_704, kl_700, kl_701, kl_702, kl_703, \
                         kl_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -hl_700[k]
                   + f_0 * kl_700[k];

        t_701[k] = -hl_701[k]
                   + f_0 * kl_701[k];

        t_702[k] = -hl_702[k]
                   + f_0 * kl_702[k];

        t_703[k] = -hl_703[k]
                   + f_0 * kl_703[k];

        t_704[k] = -hl_704[k]
                   + f_0 * kl_704[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, hl_705, hl_706, hl_707, hl_708, \
                         hl_709, kl_705, kl_706, kl_707, kl_708, \
                         kl_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -hl_705[k]
                   + f_0 * kl_705[k];

        t_706[k] = -hl_706[k]
                   + f_0 * kl_706[k];

        t_707[k] = -hl_707[k]
                   + f_0 * kl_707[k];

        t_708[k] = -hl_708[k]
                   + f_0 * kl_708[k];

        t_709[k] = -hl_709[k]
                   + f_0 * kl_709[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, hl_710, hl_711, hl_712, hl_713, \
                         hl_714, kl_710, kl_711, kl_712, kl_713, \
                         kl_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -hl_710[k]
                   + f_0 * kl_710[k];

        t_711[k] = -hl_711[k]
                   + f_0 * kl_711[k];

        t_712[k] = -hl_712[k]
                   + f_0 * kl_712[k];

        t_713[k] = -hl_713[k]
                   + f_0 * kl_713[k];

        t_714[k] = -hl_714[k]
                   + f_0 * kl_714[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, hl_715, hl_716, hl_717, hl_718, \
                         hl_719, kl_715, kl_716, kl_717, kl_718, \
                         kl_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -hl_715[k]
                   + f_0 * kl_715[k];

        t_716[k] = -hl_716[k]
                   + f_0 * kl_716[k];

        t_717[k] = -hl_717[k]
                   + f_0 * kl_717[k];

        t_718[k] = -hl_718[k]
                   + f_0 * kl_718[k];

        t_719[k] = -hl_719[k]
                   + f_0 * kl_719[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, hl_720, hl_721, hl_722, hl_723, \
                         hl_724, kl_720, kl_721, kl_722, kl_723, \
                         kl_724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -hl_720[k]
                   + f_0 * kl_720[k];

        t_721[k] = -hl_721[k]
                   + f_0 * kl_721[k];

        t_722[k] = -hl_722[k]
                   + f_0 * kl_722[k];

        t_723[k] = -hl_723[k]
                   + f_0 * kl_723[k];

        t_724[k] = -hl_724[k]
                   + f_0 * kl_724[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, hl_725, hl_726, hl_727, hl_728, \
                         hl_729, kl_725, kl_726, kl_727, kl_728, \
                         kl_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -hl_725[k]
                   + f_0 * kl_725[k];

        t_726[k] = -hl_726[k]
                   + f_0 * kl_726[k];

        t_727[k] = -hl_727[k]
                   + f_0 * kl_727[k];

        t_728[k] = -hl_728[k]
                   + f_0 * kl_728[k];

        t_729[k] = -hl_729[k]
                   + f_0 * kl_729[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, hl_730, hl_731, hl_732, hl_733, \
                         hl_734, kl_730, kl_731, kl_732, kl_733, \
                         kl_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -hl_730[k]
                   + f_0 * kl_730[k];

        t_731[k] = -hl_731[k]
                   + f_0 * kl_731[k];

        t_732[k] = -hl_732[k]
                   + f_0 * kl_732[k];

        t_733[k] = -hl_733[k]
                   + f_0 * kl_733[k];

        t_734[k] = -hl_734[k]
                   + f_0 * kl_734[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, hl_735, hl_736, hl_737, hl_738, \
                         hl_739, kl_735, kl_736, kl_737, kl_738, \
                         kl_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -hl_735[k]
                   + f_0 * kl_735[k];

        t_736[k] = -hl_736[k]
                   + f_0 * kl_736[k];

        t_737[k] = -hl_737[k]
                   + f_0 * kl_737[k];

        t_738[k] = -hl_738[k]
                   + f_0 * kl_738[k];

        t_739[k] = -hl_739[k]
                   + f_0 * kl_739[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, hl_740, hl_741, hl_742, hl_743, \
                         hl_744, kl_740, kl_741, kl_742, kl_743, \
                         kl_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -hl_740[k]
                   + f_0 * kl_740[k];

        t_741[k] = -hl_741[k]
                   + f_0 * kl_741[k];

        t_742[k] = -hl_742[k]
                   + f_0 * kl_742[k];

        t_743[k] = -hl_743[k]
                   + f_0 * kl_743[k];

        t_744[k] = -hl_744[k]
                   + f_0 * kl_744[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, hl_745, hl_746, hl_747, hl_748, \
                         hl_749, kl_745, kl_746, kl_747, kl_748, \
                         kl_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -hl_745[k]
                   + f_0 * kl_745[k];

        t_746[k] = -hl_746[k]
                   + f_0 * kl_746[k];

        t_747[k] = -hl_747[k]
                   + f_0 * kl_747[k];

        t_748[k] = -hl_748[k]
                   + f_0 * kl_748[k];

        t_749[k] = -hl_749[k]
                   + f_0 * kl_749[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
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

    const auto *hl_750 = buffer.data(hl + 750);
    const auto *hl_751 = buffer.data(hl + 751);
    const auto *hl_752 = buffer.data(hl + 752);
    const auto *hl_753 = buffer.data(hl + 753);
    const auto *hl_754 = buffer.data(hl + 754);
    const auto *hl_755 = buffer.data(hl + 755);
    const auto *hl_756 = buffer.data(hl + 756);
    const auto *hl_757 = buffer.data(hl + 757);
    const auto *hl_758 = buffer.data(hl + 758);
    const auto *hl_759 = buffer.data(hl + 759);
    const auto *hl_760 = buffer.data(hl + 760);
    const auto *hl_761 = buffer.data(hl + 761);
    const auto *hl_762 = buffer.data(hl + 762);
    const auto *hl_763 = buffer.data(hl + 763);
    const auto *hl_764 = buffer.data(hl + 764);
    const auto *hl_765 = buffer.data(hl + 765);
    const auto *hl_766 = buffer.data(hl + 766);
    const auto *hl_767 = buffer.data(hl + 767);
    const auto *hl_768 = buffer.data(hl + 768);
    const auto *hl_769 = buffer.data(hl + 769);
    const auto *hl_770 = buffer.data(hl + 770);
    const auto *hl_771 = buffer.data(hl + 771);
    const auto *hl_772 = buffer.data(hl + 772);
    const auto *hl_773 = buffer.data(hl + 773);
    const auto *hl_774 = buffer.data(hl + 774);
    const auto *hl_775 = buffer.data(hl + 775);
    const auto *hl_776 = buffer.data(hl + 776);
    const auto *hl_777 = buffer.data(hl + 777);
    const auto *hl_778 = buffer.data(hl + 778);
    const auto *hl_779 = buffer.data(hl + 779);
    const auto *hl_780 = buffer.data(hl + 780);
    const auto *hl_781 = buffer.data(hl + 781);
    const auto *hl_782 = buffer.data(hl + 782);
    const auto *hl_783 = buffer.data(hl + 783);
    const auto *hl_784 = buffer.data(hl + 784);
    const auto *hl_785 = buffer.data(hl + 785);
    const auto *hl_786 = buffer.data(hl + 786);
    const auto *hl_787 = buffer.data(hl + 787);
    const auto *hl_788 = buffer.data(hl + 788);
    const auto *hl_789 = buffer.data(hl + 789);
    const auto *hl_790 = buffer.data(hl + 790);
    const auto *hl_791 = buffer.data(hl + 791);
    const auto *hl_792 = buffer.data(hl + 792);
    const auto *hl_793 = buffer.data(hl + 793);
    const auto *hl_794 = buffer.data(hl + 794);
    const auto *hl_795 = buffer.data(hl + 795);
    const auto *hl_796 = buffer.data(hl + 796);
    const auto *hl_797 = buffer.data(hl + 797);
    const auto *hl_798 = buffer.data(hl + 798);
    const auto *hl_799 = buffer.data(hl + 799);
    const auto *hl_800 = buffer.data(hl + 800);
    const auto *hl_801 = buffer.data(hl + 801);
    const auto *hl_802 = buffer.data(hl + 802);
    const auto *hl_803 = buffer.data(hl + 803);
    const auto *hl_804 = buffer.data(hl + 804);
    const auto *hl_805 = buffer.data(hl + 805);
    const auto *hl_806 = buffer.data(hl + 806);
    const auto *hl_807 = buffer.data(hl + 807);
    const auto *hl_808 = buffer.data(hl + 808);
    const auto *hl_809 = buffer.data(hl + 809);
    const auto *hl_810 = buffer.data(hl + 810);
    const auto *hl_811 = buffer.data(hl + 811);
    const auto *hl_812 = buffer.data(hl + 812);
    const auto *hl_813 = buffer.data(hl + 813);
    const auto *hl_814 = buffer.data(hl + 814);
    const auto *hl_815 = buffer.data(hl + 815);
    const auto *hl_816 = buffer.data(hl + 816);
    const auto *hl_817 = buffer.data(hl + 817);
    const auto *hl_818 = buffer.data(hl + 818);
    const auto *hl_819 = buffer.data(hl + 819);
    const auto *hl_820 = buffer.data(hl + 820);
    const auto *hl_821 = buffer.data(hl + 821);
    const auto *hl_822 = buffer.data(hl + 822);
    const auto *hl_823 = buffer.data(hl + 823);
    const auto *hl_824 = buffer.data(hl + 824);
    const auto *hl_825 = buffer.data(hl + 825);
    const auto *hl_826 = buffer.data(hl + 826);
    const auto *hl_827 = buffer.data(hl + 827);
    const auto *hl_828 = buffer.data(hl + 828);
    const auto *hl_829 = buffer.data(hl + 829);
    const auto *hl_830 = buffer.data(hl + 830);
    const auto *hl_831 = buffer.data(hl + 831);
    const auto *hl_832 = buffer.data(hl + 832);
    const auto *hl_833 = buffer.data(hl + 833);
    const auto *hl_834 = buffer.data(hl + 834);
    const auto *hl_835 = buffer.data(hl + 835);
    const auto *hl_836 = buffer.data(hl + 836);
    const auto *hl_837 = buffer.data(hl + 837);
    const auto *hl_838 = buffer.data(hl + 838);
    const auto *hl_839 = buffer.data(hl + 839);
    const auto *hl_840 = buffer.data(hl + 840);
    const auto *hl_841 = buffer.data(hl + 841);
    const auto *hl_842 = buffer.data(hl + 842);
    const auto *hl_843 = buffer.data(hl + 843);
    const auto *hl_844 = buffer.data(hl + 844);
    const auto *hl_845 = buffer.data(hl + 845);
    const auto *hl_846 = buffer.data(hl + 846);
    const auto *hl_847 = buffer.data(hl + 847);
    const auto *hl_848 = buffer.data(hl + 848);
    const auto *hl_849 = buffer.data(hl + 849);
    const auto *hl_850 = buffer.data(hl + 850);
    const auto *hl_851 = buffer.data(hl + 851);
    const auto *hl_852 = buffer.data(hl + 852);
    const auto *hl_853 = buffer.data(hl + 853);
    const auto *hl_854 = buffer.data(hl + 854);
    const auto *hl_855 = buffer.data(hl + 855);
    const auto *hl_856 = buffer.data(hl + 856);
    const auto *hl_857 = buffer.data(hl + 857);
    const auto *hl_858 = buffer.data(hl + 858);
    const auto *hl_859 = buffer.data(hl + 859);
    const auto *hl_860 = buffer.data(hl + 860);
    const auto *hl_861 = buffer.data(hl + 861);
    const auto *hl_862 = buffer.data(hl + 862);
    const auto *hl_863 = buffer.data(hl + 863);
    const auto *hl_864 = buffer.data(hl + 864);
    const auto *hl_865 = buffer.data(hl + 865);
    const auto *hl_866 = buffer.data(hl + 866);
    const auto *hl_867 = buffer.data(hl + 867);
    const auto *hl_868 = buffer.data(hl + 868);
    const auto *hl_869 = buffer.data(hl + 869);
    const auto *hl_870 = buffer.data(hl + 870);
    const auto *hl_871 = buffer.data(hl + 871);
    const auto *hl_872 = buffer.data(hl + 872);
    const auto *hl_873 = buffer.data(hl + 873);
    const auto *hl_874 = buffer.data(hl + 874);
    const auto *hl_875 = buffer.data(hl + 875);
    const auto *hl_876 = buffer.data(hl + 876);
    const auto *hl_877 = buffer.data(hl + 877);
    const auto *hl_878 = buffer.data(hl + 878);
    const auto *hl_879 = buffer.data(hl + 879);
    const auto *hl_880 = buffer.data(hl + 880);
    const auto *hl_881 = buffer.data(hl + 881);
    const auto *hl_882 = buffer.data(hl + 882);
    const auto *hl_883 = buffer.data(hl + 883);
    const auto *hl_884 = buffer.data(hl + 884);
    const auto *hl_885 = buffer.data(hl + 885);
    const auto *hl_886 = buffer.data(hl + 886);
    const auto *hl_887 = buffer.data(hl + 887);
    const auto *hl_888 = buffer.data(hl + 888);
    const auto *hl_889 = buffer.data(hl + 889);
    const auto *hl_890 = buffer.data(hl + 890);
    const auto *hl_891 = buffer.data(hl + 891);
    const auto *hl_892 = buffer.data(hl + 892);
    const auto *hl_893 = buffer.data(hl + 893);
    const auto *hl_894 = buffer.data(hl + 894);
    const auto *hl_895 = buffer.data(hl + 895);
    const auto *hl_896 = buffer.data(hl + 896);
    const auto *hl_897 = buffer.data(hl + 897);
    const auto *hl_898 = buffer.data(hl + 898);
    const auto *hl_899 = buffer.data(hl + 899);

    const auto *kl_750 = buffer.data(kl + 750);
    const auto *kl_751 = buffer.data(kl + 751);
    const auto *kl_752 = buffer.data(kl + 752);
    const auto *kl_753 = buffer.data(kl + 753);
    const auto *kl_754 = buffer.data(kl + 754);
    const auto *kl_755 = buffer.data(kl + 755);
    const auto *kl_756 = buffer.data(kl + 756);
    const auto *kl_757 = buffer.data(kl + 757);
    const auto *kl_758 = buffer.data(kl + 758);
    const auto *kl_759 = buffer.data(kl + 759);
    const auto *kl_760 = buffer.data(kl + 760);
    const auto *kl_761 = buffer.data(kl + 761);
    const auto *kl_762 = buffer.data(kl + 762);
    const auto *kl_763 = buffer.data(kl + 763);
    const auto *kl_764 = buffer.data(kl + 764);
    const auto *kl_765 = buffer.data(kl + 765);
    const auto *kl_766 = buffer.data(kl + 766);
    const auto *kl_767 = buffer.data(kl + 767);
    const auto *kl_768 = buffer.data(kl + 768);
    const auto *kl_769 = buffer.data(kl + 769);
    const auto *kl_770 = buffer.data(kl + 770);
    const auto *kl_771 = buffer.data(kl + 771);
    const auto *kl_772 = buffer.data(kl + 772);
    const auto *kl_773 = buffer.data(kl + 773);
    const auto *kl_774 = buffer.data(kl + 774);
    const auto *kl_775 = buffer.data(kl + 775);
    const auto *kl_776 = buffer.data(kl + 776);
    const auto *kl_777 = buffer.data(kl + 777);
    const auto *kl_778 = buffer.data(kl + 778);
    const auto *kl_779 = buffer.data(kl + 779);
    const auto *kl_780 = buffer.data(kl + 780);
    const auto *kl_781 = buffer.data(kl + 781);
    const auto *kl_782 = buffer.data(kl + 782);
    const auto *kl_783 = buffer.data(kl + 783);
    const auto *kl_784 = buffer.data(kl + 784);
    const auto *kl_785 = buffer.data(kl + 785);
    const auto *kl_786 = buffer.data(kl + 786);
    const auto *kl_787 = buffer.data(kl + 787);
    const auto *kl_788 = buffer.data(kl + 788);
    const auto *kl_789 = buffer.data(kl + 789);
    const auto *kl_790 = buffer.data(kl + 790);
    const auto *kl_791 = buffer.data(kl + 791);
    const auto *kl_792 = buffer.data(kl + 792);
    const auto *kl_793 = buffer.data(kl + 793);
    const auto *kl_794 = buffer.data(kl + 794);
    const auto *kl_795 = buffer.data(kl + 795);
    const auto *kl_796 = buffer.data(kl + 796);
    const auto *kl_797 = buffer.data(kl + 797);
    const auto *kl_798 = buffer.data(kl + 798);
    const auto *kl_799 = buffer.data(kl + 799);
    const auto *kl_800 = buffer.data(kl + 800);
    const auto *kl_801 = buffer.data(kl + 801);
    const auto *kl_802 = buffer.data(kl + 802);
    const auto *kl_803 = buffer.data(kl + 803);
    const auto *kl_804 = buffer.data(kl + 804);
    const auto *kl_805 = buffer.data(kl + 805);
    const auto *kl_806 = buffer.data(kl + 806);
    const auto *kl_807 = buffer.data(kl + 807);
    const auto *kl_808 = buffer.data(kl + 808);
    const auto *kl_809 = buffer.data(kl + 809);
    const auto *kl_810 = buffer.data(kl + 810);
    const auto *kl_811 = buffer.data(kl + 811);
    const auto *kl_812 = buffer.data(kl + 812);
    const auto *kl_813 = buffer.data(kl + 813);
    const auto *kl_814 = buffer.data(kl + 814);
    const auto *kl_815 = buffer.data(kl + 815);
    const auto *kl_816 = buffer.data(kl + 816);
    const auto *kl_817 = buffer.data(kl + 817);
    const auto *kl_818 = buffer.data(kl + 818);
    const auto *kl_819 = buffer.data(kl + 819);
    const auto *kl_820 = buffer.data(kl + 820);
    const auto *kl_821 = buffer.data(kl + 821);
    const auto *kl_822 = buffer.data(kl + 822);
    const auto *kl_823 = buffer.data(kl + 823);
    const auto *kl_824 = buffer.data(kl + 824);
    const auto *kl_825 = buffer.data(kl + 825);
    const auto *kl_826 = buffer.data(kl + 826);
    const auto *kl_827 = buffer.data(kl + 827);
    const auto *kl_828 = buffer.data(kl + 828);
    const auto *kl_829 = buffer.data(kl + 829);
    const auto *kl_830 = buffer.data(kl + 830);
    const auto *kl_831 = buffer.data(kl + 831);
    const auto *kl_832 = buffer.data(kl + 832);
    const auto *kl_833 = buffer.data(kl + 833);
    const auto *kl_834 = buffer.data(kl + 834);
    const auto *kl_835 = buffer.data(kl + 835);
    const auto *kl_836 = buffer.data(kl + 836);
    const auto *kl_837 = buffer.data(kl + 837);
    const auto *kl_838 = buffer.data(kl + 838);
    const auto *kl_839 = buffer.data(kl + 839);
    const auto *kl_840 = buffer.data(kl + 840);
    const auto *kl_841 = buffer.data(kl + 841);
    const auto *kl_842 = buffer.data(kl + 842);
    const auto *kl_843 = buffer.data(kl + 843);
    const auto *kl_844 = buffer.data(kl + 844);
    const auto *kl_845 = buffer.data(kl + 845);
    const auto *kl_846 = buffer.data(kl + 846);
    const auto *kl_847 = buffer.data(kl + 847);
    const auto *kl_848 = buffer.data(kl + 848);
    const auto *kl_849 = buffer.data(kl + 849);
    const auto *kl_850 = buffer.data(kl + 850);
    const auto *kl_851 = buffer.data(kl + 851);
    const auto *kl_852 = buffer.data(kl + 852);
    const auto *kl_853 = buffer.data(kl + 853);
    const auto *kl_854 = buffer.data(kl + 854);
    const auto *kl_855 = buffer.data(kl + 855);
    const auto *kl_856 = buffer.data(kl + 856);
    const auto *kl_857 = buffer.data(kl + 857);
    const auto *kl_858 = buffer.data(kl + 858);
    const auto *kl_859 = buffer.data(kl + 859);
    const auto *kl_860 = buffer.data(kl + 860);
    const auto *kl_861 = buffer.data(kl + 861);
    const auto *kl_862 = buffer.data(kl + 862);
    const auto *kl_863 = buffer.data(kl + 863);
    const auto *kl_864 = buffer.data(kl + 864);
    const auto *kl_865 = buffer.data(kl + 865);
    const auto *kl_866 = buffer.data(kl + 866);
    const auto *kl_867 = buffer.data(kl + 867);
    const auto *kl_868 = buffer.data(kl + 868);
    const auto *kl_869 = buffer.data(kl + 869);
    const auto *kl_870 = buffer.data(kl + 870);
    const auto *kl_871 = buffer.data(kl + 871);
    const auto *kl_872 = buffer.data(kl + 872);
    const auto *kl_873 = buffer.data(kl + 873);
    const auto *kl_874 = buffer.data(kl + 874);
    const auto *kl_875 = buffer.data(kl + 875);
    const auto *kl_876 = buffer.data(kl + 876);
    const auto *kl_877 = buffer.data(kl + 877);
    const auto *kl_878 = buffer.data(kl + 878);
    const auto *kl_879 = buffer.data(kl + 879);
    const auto *kl_880 = buffer.data(kl + 880);
    const auto *kl_881 = buffer.data(kl + 881);
    const auto *kl_882 = buffer.data(kl + 882);
    const auto *kl_883 = buffer.data(kl + 883);
    const auto *kl_884 = buffer.data(kl + 884);
    const auto *kl_885 = buffer.data(kl + 885);
    const auto *kl_886 = buffer.data(kl + 886);
    const auto *kl_887 = buffer.data(kl + 887);
    const auto *kl_888 = buffer.data(kl + 888);
    const auto *kl_889 = buffer.data(kl + 889);
    const auto *kl_890 = buffer.data(kl + 890);
    const auto *kl_891 = buffer.data(kl + 891);
    const auto *kl_892 = buffer.data(kl + 892);
    const auto *kl_893 = buffer.data(kl + 893);
    const auto *kl_894 = buffer.data(kl + 894);
    const auto *kl_895 = buffer.data(kl + 895);
    const auto *kl_896 = buffer.data(kl + 896);
    const auto *kl_897 = buffer.data(kl + 897);
    const auto *kl_898 = buffer.data(kl + 898);
    const auto *kl_899 = buffer.data(kl + 899);

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, hl_750, hl_751, hl_752, hl_753, \
                         hl_754, kl_750, kl_751, kl_752, kl_753, \
                         kl_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -hl_750[k]
                   + f_0 * kl_750[k];

        t_751[k] = -hl_751[k]
                   + f_0 * kl_751[k];

        t_752[k] = -hl_752[k]
                   + f_0 * kl_752[k];

        t_753[k] = -hl_753[k]
                   + f_0 * kl_753[k];

        t_754[k] = -hl_754[k]
                   + f_0 * kl_754[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, hl_755, hl_756, hl_757, hl_758, \
                         hl_759, kl_755, kl_756, kl_757, kl_758, \
                         kl_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -hl_755[k]
                   + f_0 * kl_755[k];

        t_756[k] = -hl_756[k]
                   + f_0 * kl_756[k];

        t_757[k] = -hl_757[k]
                   + f_0 * kl_757[k];

        t_758[k] = -hl_758[k]
                   + f_0 * kl_758[k];

        t_759[k] = -hl_759[k]
                   + f_0 * kl_759[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, t_764, hl_760, hl_761, hl_762, hl_763, \
                         hl_764, kl_760, kl_761, kl_762, kl_763, \
                         kl_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = -hl_760[k]
                   + f_0 * kl_760[k];

        t_761[k] = -hl_761[k]
                   + f_0 * kl_761[k];

        t_762[k] = -hl_762[k]
                   + f_0 * kl_762[k];

        t_763[k] = -hl_763[k]
                   + f_0 * kl_763[k];

        t_764[k] = -hl_764[k]
                   + f_0 * kl_764[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, t_769, hl_765, hl_766, hl_767, hl_768, \
                         hl_769, kl_765, kl_766, kl_767, kl_768, \
                         kl_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = -hl_765[k]
                   + f_0 * kl_765[k];

        t_766[k] = -hl_766[k]
                   + f_0 * kl_766[k];

        t_767[k] = -hl_767[k]
                   + f_0 * kl_767[k];

        t_768[k] = -hl_768[k]
                   + f_0 * kl_768[k];

        t_769[k] = -hl_769[k]
                   + f_0 * kl_769[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, hl_770, hl_771, hl_772, hl_773, \
                         hl_774, kl_770, kl_771, kl_772, kl_773, \
                         kl_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = -hl_770[k]
                   + f_0 * kl_770[k];

        t_771[k] = -hl_771[k]
                   + f_0 * kl_771[k];

        t_772[k] = -hl_772[k]
                   + f_0 * kl_772[k];

        t_773[k] = -hl_773[k]
                   + f_0 * kl_773[k];

        t_774[k] = -hl_774[k]
                   + f_0 * kl_774[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, hl_775, hl_776, hl_777, hl_778, \
                         hl_779, kl_775, kl_776, kl_777, kl_778, \
                         kl_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = -hl_775[k]
                   + f_0 * kl_775[k];

        t_776[k] = -hl_776[k]
                   + f_0 * kl_776[k];

        t_777[k] = -hl_777[k]
                   + f_0 * kl_777[k];

        t_778[k] = -hl_778[k]
                   + f_0 * kl_778[k];

        t_779[k] = -hl_779[k]
                   + f_0 * kl_779[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, t_784, hl_780, hl_781, hl_782, hl_783, \
                         hl_784, kl_780, kl_781, kl_782, kl_783, \
                         kl_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = -hl_780[k]
                   + f_0 * kl_780[k];

        t_781[k] = -hl_781[k]
                   + f_0 * kl_781[k];

        t_782[k] = -hl_782[k]
                   + f_0 * kl_782[k];

        t_783[k] = -hl_783[k]
                   + f_0 * kl_783[k];

        t_784[k] = -hl_784[k]
                   + f_0 * kl_784[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, t_789, hl_785, hl_786, hl_787, hl_788, \
                         hl_789, kl_785, kl_786, kl_787, kl_788, \
                         kl_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = -hl_785[k]
                   + f_0 * kl_785[k];

        t_786[k] = -hl_786[k]
                   + f_0 * kl_786[k];

        t_787[k] = -hl_787[k]
                   + f_0 * kl_787[k];

        t_788[k] = -hl_788[k]
                   + f_0 * kl_788[k];

        t_789[k] = -hl_789[k]
                   + f_0 * kl_789[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, t_794, hl_790, hl_791, hl_792, hl_793, \
                         hl_794, kl_790, kl_791, kl_792, kl_793, \
                         kl_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = -hl_790[k]
                   + f_0 * kl_790[k];

        t_791[k] = -hl_791[k]
                   + f_0 * kl_791[k];

        t_792[k] = -hl_792[k]
                   + f_0 * kl_792[k];

        t_793[k] = -hl_793[k]
                   + f_0 * kl_793[k];

        t_794[k] = -hl_794[k]
                   + f_0 * kl_794[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, hl_795, hl_796, hl_797, hl_798, \
                         hl_799, kl_795, kl_796, kl_797, kl_798, \
                         kl_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = -hl_795[k]
                   + f_0 * kl_795[k];

        t_796[k] = -hl_796[k]
                   + f_0 * kl_796[k];

        t_797[k] = -hl_797[k]
                   + f_0 * kl_797[k];

        t_798[k] = -hl_798[k]
                   + f_0 * kl_798[k];

        t_799[k] = -hl_799[k]
                   + f_0 * kl_799[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, hl_800, hl_801, hl_802, hl_803, \
                         hl_804, kl_800, kl_801, kl_802, kl_803, \
                         kl_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = -hl_800[k]
                   + f_0 * kl_800[k];

        t_801[k] = -hl_801[k]
                   + f_0 * kl_801[k];

        t_802[k] = -hl_802[k]
                   + f_0 * kl_802[k];

        t_803[k] = -hl_803[k]
                   + f_0 * kl_803[k];

        t_804[k] = -hl_804[k]
                   + f_0 * kl_804[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, hl_805, hl_806, hl_807, hl_808, \
                         hl_809, kl_805, kl_806, kl_807, kl_808, \
                         kl_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = -hl_805[k]
                   + f_0 * kl_805[k];

        t_806[k] = -hl_806[k]
                   + f_0 * kl_806[k];

        t_807[k] = -hl_807[k]
                   + f_0 * kl_807[k];

        t_808[k] = -hl_808[k]
                   + f_0 * kl_808[k];

        t_809[k] = -hl_809[k]
                   + f_0 * kl_809[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, hl_810, hl_811, hl_812, hl_813, \
                         hl_814, kl_810, kl_811, kl_812, kl_813, \
                         kl_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = -hl_810[k]
                   + f_0 * kl_810[k];

        t_811[k] = -hl_811[k]
                   + f_0 * kl_811[k];

        t_812[k] = -hl_812[k]
                   + f_0 * kl_812[k];

        t_813[k] = -hl_813[k]
                   + f_0 * kl_813[k];

        t_814[k] = -hl_814[k]
                   + f_0 * kl_814[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, hl_815, hl_816, hl_817, hl_818, \
                         hl_819, kl_815, kl_816, kl_817, kl_818, \
                         kl_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = -hl_815[k]
                   + f_0 * kl_815[k];

        t_816[k] = -hl_816[k]
                   + f_0 * kl_816[k];

        t_817[k] = -hl_817[k]
                   + f_0 * kl_817[k];

        t_818[k] = -hl_818[k]
                   + f_0 * kl_818[k];

        t_819[k] = -hl_819[k]
                   + f_0 * kl_819[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, t_824, hl_820, hl_821, hl_822, hl_823, \
                         hl_824, kl_820, kl_821, kl_822, kl_823, \
                         kl_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = -hl_820[k]
                   + f_0 * kl_820[k];

        t_821[k] = -hl_821[k]
                   + f_0 * kl_821[k];

        t_822[k] = -hl_822[k]
                   + f_0 * kl_822[k];

        t_823[k] = -hl_823[k]
                   + f_0 * kl_823[k];

        t_824[k] = -hl_824[k]
                   + f_0 * kl_824[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, hl_825, hl_826, hl_827, hl_828, \
                         hl_829, kl_825, kl_826, kl_827, kl_828, \
                         kl_829 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = -hl_825[k]
                   + f_0 * kl_825[k];

        t_826[k] = -hl_826[k]
                   + f_0 * kl_826[k];

        t_827[k] = -hl_827[k]
                   + f_0 * kl_827[k];

        t_828[k] = -hl_828[k]
                   + f_0 * kl_828[k];

        t_829[k] = -hl_829[k]
                   + f_0 * kl_829[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, hl_830, hl_831, hl_832, hl_833, \
                         hl_834, kl_830, kl_831, kl_832, kl_833, \
                         kl_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = -hl_830[k]
                   + f_0 * kl_830[k];

        t_831[k] = -hl_831[k]
                   + f_0 * kl_831[k];

        t_832[k] = -hl_832[k]
                   + f_0 * kl_832[k];

        t_833[k] = -hl_833[k]
                   + f_0 * kl_833[k];

        t_834[k] = -hl_834[k]
                   + f_0 * kl_834[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, t_838, t_839, hl_835, hl_836, hl_837, hl_838, \
                         hl_839, kl_835, kl_836, kl_837, kl_838, \
                         kl_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = -hl_835[k]
                   + f_0 * kl_835[k];

        t_836[k] = -hl_836[k]
                   + f_0 * kl_836[k];

        t_837[k] = -hl_837[k]
                   + f_0 * kl_837[k];

        t_838[k] = -hl_838[k]
                   + f_0 * kl_838[k];

        t_839[k] = -hl_839[k]
                   + f_0 * kl_839[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, hl_840, hl_841, hl_842, hl_843, \
                         hl_844, kl_840, kl_841, kl_842, kl_843, \
                         kl_844 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = -hl_840[k]
                   + f_0 * kl_840[k];

        t_841[k] = -hl_841[k]
                   + f_0 * kl_841[k];

        t_842[k] = -hl_842[k]
                   + f_0 * kl_842[k];

        t_843[k] = -hl_843[k]
                   + f_0 * kl_843[k];

        t_844[k] = -hl_844[k]
                   + f_0 * kl_844[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, t_849, hl_845, hl_846, hl_847, hl_848, \
                         hl_849, kl_845, kl_846, kl_847, kl_848, \
                         kl_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = -hl_845[k]
                   + f_0 * kl_845[k];

        t_846[k] = -hl_846[k]
                   + f_0 * kl_846[k];

        t_847[k] = -hl_847[k]
                   + f_0 * kl_847[k];

        t_848[k] = -hl_848[k]
                   + f_0 * kl_848[k];

        t_849[k] = -hl_849[k]
                   + f_0 * kl_849[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, hl_850, hl_851, hl_852, hl_853, \
                         hl_854, kl_850, kl_851, kl_852, kl_853, \
                         kl_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = -hl_850[k]
                   + f_0 * kl_850[k];

        t_851[k] = -hl_851[k]
                   + f_0 * kl_851[k];

        t_852[k] = -hl_852[k]
                   + f_0 * kl_852[k];

        t_853[k] = -hl_853[k]
                   + f_0 * kl_853[k];

        t_854[k] = -hl_854[k]
                   + f_0 * kl_854[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, hl_855, hl_856, hl_857, hl_858, \
                         hl_859, kl_855, kl_856, kl_857, kl_858, \
                         kl_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = -hl_855[k]
                   + f_0 * kl_855[k];

        t_856[k] = -hl_856[k]
                   + f_0 * kl_856[k];

        t_857[k] = -hl_857[k]
                   + f_0 * kl_857[k];

        t_858[k] = -hl_858[k]
                   + f_0 * kl_858[k];

        t_859[k] = -hl_859[k]
                   + f_0 * kl_859[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, t_864, hl_860, hl_861, hl_862, hl_863, \
                         hl_864, kl_860, kl_861, kl_862, kl_863, \
                         kl_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = -hl_860[k]
                   + f_0 * kl_860[k];

        t_861[k] = -hl_861[k]
                   + f_0 * kl_861[k];

        t_862[k] = -hl_862[k]
                   + f_0 * kl_862[k];

        t_863[k] = -hl_863[k]
                   + f_0 * kl_863[k];

        t_864[k] = -hl_864[k]
                   + f_0 * kl_864[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, t_869, hl_865, hl_866, hl_867, hl_868, \
                         hl_869, kl_865, kl_866, kl_867, kl_868, \
                         kl_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = -hl_865[k]
                   + f_0 * kl_865[k];

        t_866[k] = -hl_866[k]
                   + f_0 * kl_866[k];

        t_867[k] = -hl_867[k]
                   + f_0 * kl_867[k];

        t_868[k] = -hl_868[k]
                   + f_0 * kl_868[k];

        t_869[k] = -hl_869[k]
                   + f_0 * kl_869[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, t_874, hl_870, hl_871, hl_872, hl_873, \
                         hl_874, kl_870, kl_871, kl_872, kl_873, \
                         kl_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = -hl_870[k]
                   + f_0 * kl_870[k];

        t_871[k] = -hl_871[k]
                   + f_0 * kl_871[k];

        t_872[k] = -hl_872[k]
                   + f_0 * kl_872[k];

        t_873[k] = -hl_873[k]
                   + f_0 * kl_873[k];

        t_874[k] = -hl_874[k]
                   + f_0 * kl_874[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, t_879, hl_875, hl_876, hl_877, hl_878, \
                         hl_879, kl_875, kl_876, kl_877, kl_878, \
                         kl_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = -hl_875[k]
                   + f_0 * kl_875[k];

        t_876[k] = -hl_876[k]
                   + f_0 * kl_876[k];

        t_877[k] = -hl_877[k]
                   + f_0 * kl_877[k];

        t_878[k] = -hl_878[k]
                   + f_0 * kl_878[k];

        t_879[k] = -hl_879[k]
                   + f_0 * kl_879[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, t_884, hl_880, hl_881, hl_882, hl_883, \
                         hl_884, kl_880, kl_881, kl_882, kl_883, \
                         kl_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = -hl_880[k]
                   + f_0 * kl_880[k];

        t_881[k] = -hl_881[k]
                   + f_0 * kl_881[k];

        t_882[k] = -hl_882[k]
                   + f_0 * kl_882[k];

        t_883[k] = -hl_883[k]
                   + f_0 * kl_883[k];

        t_884[k] = -hl_884[k]
                   + f_0 * kl_884[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, hl_885, hl_886, hl_887, hl_888, \
                         hl_889, kl_885, kl_886, kl_887, kl_888, \
                         kl_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = -hl_885[k]
                   + f_0 * kl_885[k];

        t_886[k] = -hl_886[k]
                   + f_0 * kl_886[k];

        t_887[k] = -hl_887[k]
                   + f_0 * kl_887[k];

        t_888[k] = -hl_888[k]
                   + f_0 * kl_888[k];

        t_889[k] = -hl_889[k]
                   + f_0 * kl_889[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, hl_890, hl_891, hl_892, hl_893, \
                         hl_894, kl_890, kl_891, kl_892, kl_893, \
                         kl_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = -hl_890[k]
                   + f_0 * kl_890[k];

        t_891[k] = -hl_891[k]
                   + f_0 * kl_891[k];

        t_892[k] = -hl_892[k]
                   + f_0 * kl_892[k];

        t_893[k] = -hl_893[k]
                   + f_0 * kl_893[k];

        t_894[k] = -hl_894[k]
                   + f_0 * kl_894[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, t_899, hl_895, hl_896, hl_897, hl_898, \
                         hl_899, kl_895, kl_896, kl_897, kl_898, \
                         kl_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = -hl_895[k]
                   + f_0 * kl_895[k];

        t_896[k] = -hl_896[k]
                   + f_0 * kl_896[k];

        t_897[k] = -hl_897[k]
                   + f_0 * kl_897[k];

        t_898[k] = -hl_898[k]
                   + f_0 * kl_898[k];

        t_899[k] = -hl_899[k]
                   + f_0 * kl_899[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_0_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
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

    const auto *hl_900 = buffer.data(hl + 900);
    const auto *hl_901 = buffer.data(hl + 901);
    const auto *hl_902 = buffer.data(hl + 902);
    const auto *hl_903 = buffer.data(hl + 903);
    const auto *hl_904 = buffer.data(hl + 904);
    const auto *hl_905 = buffer.data(hl + 905);
    const auto *hl_906 = buffer.data(hl + 906);
    const auto *hl_907 = buffer.data(hl + 907);
    const auto *hl_908 = buffer.data(hl + 908);
    const auto *hl_909 = buffer.data(hl + 909);
    const auto *hl_910 = buffer.data(hl + 910);
    const auto *hl_911 = buffer.data(hl + 911);
    const auto *hl_912 = buffer.data(hl + 912);
    const auto *hl_913 = buffer.data(hl + 913);
    const auto *hl_914 = buffer.data(hl + 914);
    const auto *hl_915 = buffer.data(hl + 915);
    const auto *hl_916 = buffer.data(hl + 916);
    const auto *hl_917 = buffer.data(hl + 917);
    const auto *hl_918 = buffer.data(hl + 918);
    const auto *hl_919 = buffer.data(hl + 919);
    const auto *hl_920 = buffer.data(hl + 920);
    const auto *hl_921 = buffer.data(hl + 921);
    const auto *hl_922 = buffer.data(hl + 922);
    const auto *hl_923 = buffer.data(hl + 923);
    const auto *hl_924 = buffer.data(hl + 924);
    const auto *hl_925 = buffer.data(hl + 925);
    const auto *hl_926 = buffer.data(hl + 926);
    const auto *hl_927 = buffer.data(hl + 927);
    const auto *hl_928 = buffer.data(hl + 928);
    const auto *hl_929 = buffer.data(hl + 929);
    const auto *hl_930 = buffer.data(hl + 930);
    const auto *hl_931 = buffer.data(hl + 931);
    const auto *hl_932 = buffer.data(hl + 932);
    const auto *hl_933 = buffer.data(hl + 933);
    const auto *hl_934 = buffer.data(hl + 934);
    const auto *hl_935 = buffer.data(hl + 935);
    const auto *hl_936 = buffer.data(hl + 936);
    const auto *hl_937 = buffer.data(hl + 937);
    const auto *hl_938 = buffer.data(hl + 938);
    const auto *hl_939 = buffer.data(hl + 939);
    const auto *hl_940 = buffer.data(hl + 940);
    const auto *hl_941 = buffer.data(hl + 941);
    const auto *hl_942 = buffer.data(hl + 942);
    const auto *hl_943 = buffer.data(hl + 943);
    const auto *hl_944 = buffer.data(hl + 944);

    const auto *kl_900 = buffer.data(kl + 900);
    const auto *kl_901 = buffer.data(kl + 901);
    const auto *kl_902 = buffer.data(kl + 902);
    const auto *kl_903 = buffer.data(kl + 903);
    const auto *kl_904 = buffer.data(kl + 904);
    const auto *kl_905 = buffer.data(kl + 905);
    const auto *kl_906 = buffer.data(kl + 906);
    const auto *kl_907 = buffer.data(kl + 907);
    const auto *kl_908 = buffer.data(kl + 908);
    const auto *kl_909 = buffer.data(kl + 909);
    const auto *kl_910 = buffer.data(kl + 910);
    const auto *kl_911 = buffer.data(kl + 911);
    const auto *kl_912 = buffer.data(kl + 912);
    const auto *kl_913 = buffer.data(kl + 913);
    const auto *kl_914 = buffer.data(kl + 914);
    const auto *kl_915 = buffer.data(kl + 915);
    const auto *kl_916 = buffer.data(kl + 916);
    const auto *kl_917 = buffer.data(kl + 917);
    const auto *kl_918 = buffer.data(kl + 918);
    const auto *kl_919 = buffer.data(kl + 919);
    const auto *kl_920 = buffer.data(kl + 920);
    const auto *kl_921 = buffer.data(kl + 921);
    const auto *kl_922 = buffer.data(kl + 922);
    const auto *kl_923 = buffer.data(kl + 923);
    const auto *kl_924 = buffer.data(kl + 924);
    const auto *kl_925 = buffer.data(kl + 925);
    const auto *kl_926 = buffer.data(kl + 926);
    const auto *kl_927 = buffer.data(kl + 927);
    const auto *kl_928 = buffer.data(kl + 928);
    const auto *kl_929 = buffer.data(kl + 929);
    const auto *kl_930 = buffer.data(kl + 930);
    const auto *kl_931 = buffer.data(kl + 931);
    const auto *kl_932 = buffer.data(kl + 932);
    const auto *kl_933 = buffer.data(kl + 933);
    const auto *kl_934 = buffer.data(kl + 934);
    const auto *kl_935 = buffer.data(kl + 935);
    const auto *kl_936 = buffer.data(kl + 936);
    const auto *kl_937 = buffer.data(kl + 937);
    const auto *kl_938 = buffer.data(kl + 938);
    const auto *kl_939 = buffer.data(kl + 939);
    const auto *kl_940 = buffer.data(kl + 940);
    const auto *kl_941 = buffer.data(kl + 941);
    const auto *kl_942 = buffer.data(kl + 942);
    const auto *kl_943 = buffer.data(kl + 943);
    const auto *kl_944 = buffer.data(kl + 944);
    const auto *kl_945 = buffer.data(kl + 945);
    const auto *kl_946 = buffer.data(kl + 946);
    const auto *kl_947 = buffer.data(kl + 947);
    const auto *kl_948 = buffer.data(kl + 948);
    const auto *kl_949 = buffer.data(kl + 949);
    const auto *kl_950 = buffer.data(kl + 950);
    const auto *kl_951 = buffer.data(kl + 951);
    const auto *kl_952 = buffer.data(kl + 952);
    const auto *kl_953 = buffer.data(kl + 953);
    const auto *kl_954 = buffer.data(kl + 954);
    const auto *kl_955 = buffer.data(kl + 955);
    const auto *kl_956 = buffer.data(kl + 956);
    const auto *kl_957 = buffer.data(kl + 957);
    const auto *kl_958 = buffer.data(kl + 958);
    const auto *kl_959 = buffer.data(kl + 959);
    const auto *kl_960 = buffer.data(kl + 960);
    const auto *kl_961 = buffer.data(kl + 961);
    const auto *kl_962 = buffer.data(kl + 962);
    const auto *kl_963 = buffer.data(kl + 963);
    const auto *kl_964 = buffer.data(kl + 964);
    const auto *kl_965 = buffer.data(kl + 965);
    const auto *kl_966 = buffer.data(kl + 966);
    const auto *kl_967 = buffer.data(kl + 967);
    const auto *kl_968 = buffer.data(kl + 968);
    const auto *kl_969 = buffer.data(kl + 969);
    const auto *kl_970 = buffer.data(kl + 970);
    const auto *kl_971 = buffer.data(kl + 971);
    const auto *kl_972 = buffer.data(kl + 972);
    const auto *kl_973 = buffer.data(kl + 973);
    const auto *kl_974 = buffer.data(kl + 974);
    const auto *kl_975 = buffer.data(kl + 975);
    const auto *kl_976 = buffer.data(kl + 976);
    const auto *kl_977 = buffer.data(kl + 977);
    const auto *kl_978 = buffer.data(kl + 978);
    const auto *kl_979 = buffer.data(kl + 979);
    const auto *kl_980 = buffer.data(kl + 980);
    const auto *kl_981 = buffer.data(kl + 981);
    const auto *kl_982 = buffer.data(kl + 982);
    const auto *kl_983 = buffer.data(kl + 983);
    const auto *kl_984 = buffer.data(kl + 984);
    const auto *kl_985 = buffer.data(kl + 985);
    const auto *kl_986 = buffer.data(kl + 986);
    const auto *kl_987 = buffer.data(kl + 987);
    const auto *kl_988 = buffer.data(kl + 988);
    const auto *kl_989 = buffer.data(kl + 989);
    const auto *kl_990 = buffer.data(kl + 990);
    const auto *kl_991 = buffer.data(kl + 991);
    const auto *kl_992 = buffer.data(kl + 992);
    const auto *kl_993 = buffer.data(kl + 993);
    const auto *kl_994 = buffer.data(kl + 994);
    const auto *kl_995 = buffer.data(kl + 995);
    const auto *kl_996 = buffer.data(kl + 996);
    const auto *kl_997 = buffer.data(kl + 997);
    const auto *kl_998 = buffer.data(kl + 998);
    const auto *kl_999 = buffer.data(kl + 999);
    const auto *kl_1000 = buffer.data(kl + 1000);
    const auto *kl_1001 = buffer.data(kl + 1001);
    const auto *kl_1002 = buffer.data(kl + 1002);
    const auto *kl_1003 = buffer.data(kl + 1003);
    const auto *kl_1004 = buffer.data(kl + 1004);
    const auto *kl_1005 = buffer.data(kl + 1005);
    const auto *kl_1006 = buffer.data(kl + 1006);
    const auto *kl_1007 = buffer.data(kl + 1007);
    const auto *kl_1008 = buffer.data(kl + 1008);
    const auto *kl_1009 = buffer.data(kl + 1009);
    const auto *kl_1010 = buffer.data(kl + 1010);
    const auto *kl_1011 = buffer.data(kl + 1011);
    const auto *kl_1012 = buffer.data(kl + 1012);
    const auto *kl_1013 = buffer.data(kl + 1013);
    const auto *kl_1014 = buffer.data(kl + 1014);
    const auto *kl_1015 = buffer.data(kl + 1015);
    const auto *kl_1016 = buffer.data(kl + 1016);
    const auto *kl_1017 = buffer.data(kl + 1017);
    const auto *kl_1018 = buffer.data(kl + 1018);
    const auto *kl_1019 = buffer.data(kl + 1019);
    const auto *kl_1020 = buffer.data(kl + 1020);
    const auto *kl_1021 = buffer.data(kl + 1021);
    const auto *kl_1022 = buffer.data(kl + 1022);
    const auto *kl_1023 = buffer.data(kl + 1023);
    const auto *kl_1024 = buffer.data(kl + 1024);
    const auto *kl_1025 = buffer.data(kl + 1025);
    const auto *kl_1026 = buffer.data(kl + 1026);
    const auto *kl_1027 = buffer.data(kl + 1027);
    const auto *kl_1028 = buffer.data(kl + 1028);
    const auto *kl_1029 = buffer.data(kl + 1029);
    const auto *kl_1030 = buffer.data(kl + 1030);
    const auto *kl_1031 = buffer.data(kl + 1031);
    const auto *kl_1032 = buffer.data(kl + 1032);
    const auto *kl_1033 = buffer.data(kl + 1033);
    const auto *kl_1034 = buffer.data(kl + 1034);
    const auto *kl_1035 = buffer.data(kl + 1035);
    const auto *kl_1036 = buffer.data(kl + 1036);
    const auto *kl_1037 = buffer.data(kl + 1037);
    const auto *kl_1038 = buffer.data(kl + 1038);
    const auto *kl_1039 = buffer.data(kl + 1039);
    const auto *kl_1040 = buffer.data(kl + 1040);
    const auto *kl_1041 = buffer.data(kl + 1041);
    const auto *kl_1042 = buffer.data(kl + 1042);
    const auto *kl_1043 = buffer.data(kl + 1043);
    const auto *kl_1044 = buffer.data(kl + 1044);
    const auto *kl_1045 = buffer.data(kl + 1045);
    const auto *kl_1046 = buffer.data(kl + 1046);
    const auto *kl_1047 = buffer.data(kl + 1047);
    const auto *kl_1048 = buffer.data(kl + 1048);
    const auto *kl_1049 = buffer.data(kl + 1049);
    const auto *kl_1050 = buffer.data(kl + 1050);
    const auto *kl_1051 = buffer.data(kl + 1051);
    const auto *kl_1052 = buffer.data(kl + 1052);
    const auto *kl_1053 = buffer.data(kl + 1053);
    const auto *kl_1054 = buffer.data(kl + 1054);
    const auto *kl_1055 = buffer.data(kl + 1055);
    const auto *kl_1056 = buffer.data(kl + 1056);
    const auto *kl_1057 = buffer.data(kl + 1057);
    const auto *kl_1058 = buffer.data(kl + 1058);
    const auto *kl_1059 = buffer.data(kl + 1059);
    const auto *kl_1060 = buffer.data(kl + 1060);
    const auto *kl_1061 = buffer.data(kl + 1061);
    const auto *kl_1062 = buffer.data(kl + 1062);
    const auto *kl_1063 = buffer.data(kl + 1063);
    const auto *kl_1064 = buffer.data(kl + 1064);
    const auto *kl_1065 = buffer.data(kl + 1065);
    const auto *kl_1066 = buffer.data(kl + 1066);
    const auto *kl_1067 = buffer.data(kl + 1067);
    const auto *kl_1068 = buffer.data(kl + 1068);
    const auto *kl_1069 = buffer.data(kl + 1069);
    const auto *kl_1070 = buffer.data(kl + 1070);
    const auto *kl_1071 = buffer.data(kl + 1071);
    const auto *kl_1072 = buffer.data(kl + 1072);
    const auto *kl_1073 = buffer.data(kl + 1073);
    const auto *kl_1074 = buffer.data(kl + 1074);
    const auto *kl_1075 = buffer.data(kl + 1075);
    const auto *kl_1076 = buffer.data(kl + 1076);
    const auto *kl_1077 = buffer.data(kl + 1077);
    const auto *kl_1078 = buffer.data(kl + 1078);
    const auto *kl_1079 = buffer.data(kl + 1079);
    const auto *kl_1080 = buffer.data(kl + 1080);
    const auto *kl_1081 = buffer.data(kl + 1081);
    const auto *kl_1082 = buffer.data(kl + 1082);
    const auto *kl_1083 = buffer.data(kl + 1083);
    const auto *kl_1084 = buffer.data(kl + 1084);
    const auto *kl_1085 = buffer.data(kl + 1085);
    const auto *kl_1086 = buffer.data(kl + 1086);
    const auto *kl_1087 = buffer.data(kl + 1087);
    const auto *kl_1088 = buffer.data(kl + 1088);
    const auto *kl_1089 = buffer.data(kl + 1089);
    const auto *kl_1090 = buffer.data(kl + 1090);
    const auto *kl_1091 = buffer.data(kl + 1091);
    const auto *kl_1092 = buffer.data(kl + 1092);
    const auto *kl_1093 = buffer.data(kl + 1093);
    const auto *kl_1094 = buffer.data(kl + 1094);
    const auto *kl_1095 = buffer.data(kl + 1095);
    const auto *kl_1096 = buffer.data(kl + 1096);
    const auto *kl_1097 = buffer.data(kl + 1097);
    const auto *kl_1098 = buffer.data(kl + 1098);
    const auto *kl_1099 = buffer.data(kl + 1099);
    const auto *kl_1100 = buffer.data(kl + 1100);
    const auto *kl_1101 = buffer.data(kl + 1101);
    const auto *kl_1102 = buffer.data(kl + 1102);
    const auto *kl_1103 = buffer.data(kl + 1103);
    const auto *kl_1104 = buffer.data(kl + 1104);

#pragma omp simd aligned(t_900, t_901, t_902, t_903, t_904, hl_900, hl_901, hl_902, hl_903, \
                         hl_904, kl_900, kl_901, kl_902, kl_903, \
                         kl_904 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = -hl_900[k]
                   + f_0 * kl_900[k];

        t_901[k] = -hl_901[k]
                   + f_0 * kl_901[k];

        t_902[k] = -hl_902[k]
                   + f_0 * kl_902[k];

        t_903[k] = -hl_903[k]
                   + f_0 * kl_903[k];

        t_904[k] = -hl_904[k]
                   + f_0 * kl_904[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, t_909, hl_905, hl_906, hl_907, hl_908, \
                         hl_909, kl_905, kl_906, kl_907, kl_908, \
                         kl_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = -hl_905[k]
                   + f_0 * kl_905[k];

        t_906[k] = -hl_906[k]
                   + f_0 * kl_906[k];

        t_907[k] = -hl_907[k]
                   + f_0 * kl_907[k];

        t_908[k] = -hl_908[k]
                   + f_0 * kl_908[k];

        t_909[k] = -hl_909[k]
                   + f_0 * kl_909[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, t_914, hl_910, hl_911, hl_912, hl_913, \
                         hl_914, kl_910, kl_911, kl_912, kl_913, \
                         kl_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = -hl_910[k]
                   + f_0 * kl_910[k];

        t_911[k] = -hl_911[k]
                   + f_0 * kl_911[k];

        t_912[k] = -hl_912[k]
                   + f_0 * kl_912[k];

        t_913[k] = -hl_913[k]
                   + f_0 * kl_913[k];

        t_914[k] = -hl_914[k]
                   + f_0 * kl_914[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, t_919, hl_915, hl_916, hl_917, hl_918, \
                         hl_919, kl_915, kl_916, kl_917, kl_918, \
                         kl_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = -hl_915[k]
                   + f_0 * kl_915[k];

        t_916[k] = -hl_916[k]
                   + f_0 * kl_916[k];

        t_917[k] = -hl_917[k]
                   + f_0 * kl_917[k];

        t_918[k] = -hl_918[k]
                   + f_0 * kl_918[k];

        t_919[k] = -hl_919[k]
                   + f_0 * kl_919[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, t_924, hl_920, hl_921, hl_922, hl_923, \
                         hl_924, kl_920, kl_921, kl_922, kl_923, \
                         kl_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = -hl_920[k]
                   + f_0 * kl_920[k];

        t_921[k] = -hl_921[k]
                   + f_0 * kl_921[k];

        t_922[k] = -hl_922[k]
                   + f_0 * kl_922[k];

        t_923[k] = -hl_923[k]
                   + f_0 * kl_923[k];

        t_924[k] = -hl_924[k]
                   + f_0 * kl_924[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, t_929, hl_925, hl_926, hl_927, hl_928, \
                         hl_929, kl_925, kl_926, kl_927, kl_928, \
                         kl_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = -hl_925[k]
                   + f_0 * kl_925[k];

        t_926[k] = -hl_926[k]
                   + f_0 * kl_926[k];

        t_927[k] = -hl_927[k]
                   + f_0 * kl_927[k];

        t_928[k] = -hl_928[k]
                   + f_0 * kl_928[k];

        t_929[k] = -hl_929[k]
                   + f_0 * kl_929[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, hl_930, hl_931, hl_932, hl_933, \
                         hl_934, kl_930, kl_931, kl_932, kl_933, \
                         kl_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = -hl_930[k]
                   + f_0 * kl_930[k];

        t_931[k] = -hl_931[k]
                   + f_0 * kl_931[k];

        t_932[k] = -hl_932[k]
                   + f_0 * kl_932[k];

        t_933[k] = -hl_933[k]
                   + f_0 * kl_933[k];

        t_934[k] = -hl_934[k]
                   + f_0 * kl_934[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, t_939, hl_935, hl_936, hl_937, hl_938, \
                         hl_939, kl_935, kl_936, kl_937, kl_938, \
                         kl_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = -hl_935[k]
                   + f_0 * kl_935[k];

        t_936[k] = -hl_936[k]
                   + f_0 * kl_936[k];

        t_937[k] = -hl_937[k]
                   + f_0 * kl_937[k];

        t_938[k] = -hl_938[k]
                   + f_0 * kl_938[k];

        t_939[k] = -hl_939[k]
                   + f_0 * kl_939[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, t_944, hl_940, hl_941, hl_942, hl_943, \
                         hl_944, kl_940, kl_941, kl_942, kl_943, \
                         kl_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = -hl_940[k]
                   + f_0 * kl_940[k];

        t_941[k] = -hl_941[k]
                   + f_0 * kl_941[k];

        t_942[k] = -hl_942[k]
                   + f_0 * kl_942[k];

        t_943[k] = -hl_943[k]
                   + f_0 * kl_943[k];

        t_944[k] = -hl_944[k]
                   + f_0 * kl_944[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, t_949, t_950, t_951, t_952, kl_945, \
                         kl_946, kl_947, kl_948, kl_949, kl_950, kl_951, \
                         kl_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = f_0 * kl_945[k];

        t_946[k] = f_0 * kl_946[k];

        t_947[k] = f_0 * kl_947[k];

        t_948[k] = f_0 * kl_948[k];

        t_949[k] = f_0 * kl_949[k];

        t_950[k] = f_0 * kl_950[k];

        t_951[k] = f_0 * kl_951[k];

        t_952[k] = f_0 * kl_952[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, t_957, t_958, t_959, t_960, kl_953, \
                         kl_954, kl_955, kl_956, kl_957, kl_958, kl_959, \
                         kl_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_0 * kl_953[k];

        t_954[k] = f_0 * kl_954[k];

        t_955[k] = f_0 * kl_955[k];

        t_956[k] = f_0 * kl_956[k];

        t_957[k] = f_0 * kl_957[k];

        t_958[k] = f_0 * kl_958[k];

        t_959[k] = f_0 * kl_959[k];

        t_960[k] = f_0 * kl_960[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, t_964, t_965, t_966, t_967, t_968, kl_961, \
                         kl_962, kl_963, kl_964, kl_965, kl_966, kl_967, \
                         kl_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_0 * kl_961[k];

        t_962[k] = f_0 * kl_962[k];

        t_963[k] = f_0 * kl_963[k];

        t_964[k] = f_0 * kl_964[k];

        t_965[k] = f_0 * kl_965[k];

        t_966[k] = f_0 * kl_966[k];

        t_967[k] = f_0 * kl_967[k];

        t_968[k] = f_0 * kl_968[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, t_972, t_973, t_974, t_975, t_976, kl_969, \
                         kl_970, kl_971, kl_972, kl_973, kl_974, kl_975, \
                         kl_976 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_0 * kl_969[k];

        t_970[k] = f_0 * kl_970[k];

        t_971[k] = f_0 * kl_971[k];

        t_972[k] = f_0 * kl_972[k];

        t_973[k] = f_0 * kl_973[k];

        t_974[k] = f_0 * kl_974[k];

        t_975[k] = f_0 * kl_975[k];

        t_976[k] = f_0 * kl_976[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, t_980, t_981, t_982, t_983, t_984, kl_977, \
                         kl_978, kl_979, kl_980, kl_981, kl_982, kl_983, \
                         kl_984 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_0 * kl_977[k];

        t_978[k] = f_0 * kl_978[k];

        t_979[k] = f_0 * kl_979[k];

        t_980[k] = f_0 * kl_980[k];

        t_981[k] = f_0 * kl_981[k];

        t_982[k] = f_0 * kl_982[k];

        t_983[k] = f_0 * kl_983[k];

        t_984[k] = f_0 * kl_984[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, t_989, t_990, t_991, t_992, kl_985, \
                         kl_986, kl_987, kl_988, kl_989, kl_990, kl_991, \
                         kl_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_0 * kl_985[k];

        t_986[k] = f_0 * kl_986[k];

        t_987[k] = f_0 * kl_987[k];

        t_988[k] = f_0 * kl_988[k];

        t_989[k] = f_0 * kl_989[k];

        t_990[k] = f_0 * kl_990[k];

        t_991[k] = f_0 * kl_991[k];

        t_992[k] = f_0 * kl_992[k];
    }

#pragma omp simd aligned(t_993, t_994, t_995, t_996, t_997, t_998, t_999, t_1000, kl_993, \
                         kl_994, kl_995, kl_996, kl_997, kl_998, kl_999, \
                         kl_1000 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_993[k] = f_0 * kl_993[k];

        t_994[k] = f_0 * kl_994[k];

        t_995[k] = f_0 * kl_995[k];

        t_996[k] = f_0 * kl_996[k];

        t_997[k] = f_0 * kl_997[k];

        t_998[k] = f_0 * kl_998[k];

        t_999[k] = f_0 * kl_999[k];

        t_1000[k] = f_0 * kl_1000[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, t_1004, t_1005, t_1006, t_1007, t_1008, \
                         kl_1001, kl_1002, kl_1003, kl_1004, kl_1005, kl_1006, kl_1007, \
                         kl_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_0 * kl_1001[k];

        t_1002[k] = f_0 * kl_1002[k];

        t_1003[k] = f_0 * kl_1003[k];

        t_1004[k] = f_0 * kl_1004[k];

        t_1005[k] = f_0 * kl_1005[k];

        t_1006[k] = f_0 * kl_1006[k];

        t_1007[k] = f_0 * kl_1007[k];

        t_1008[k] = f_0 * kl_1008[k];
    }

#pragma omp simd aligned(t_1009, t_1010, t_1011, t_1012, t_1013, t_1014, t_1015, t_1016, \
                         kl_1009, kl_1010, kl_1011, kl_1012, kl_1013, kl_1014, kl_1015, \
                         kl_1016 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1009[k] = f_0 * kl_1009[k];

        t_1010[k] = f_0 * kl_1010[k];

        t_1011[k] = f_0 * kl_1011[k];

        t_1012[k] = f_0 * kl_1012[k];

        t_1013[k] = f_0 * kl_1013[k];

        t_1014[k] = f_0 * kl_1014[k];

        t_1015[k] = f_0 * kl_1015[k];

        t_1016[k] = f_0 * kl_1016[k];
    }

#pragma omp simd aligned(t_1017, t_1018, t_1019, t_1020, t_1021, t_1022, t_1023, t_1024, \
                         kl_1017, kl_1018, kl_1019, kl_1020, kl_1021, kl_1022, kl_1023, \
                         kl_1024 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1017[k] = f_0 * kl_1017[k];

        t_1018[k] = f_0 * kl_1018[k];

        t_1019[k] = f_0 * kl_1019[k];

        t_1020[k] = f_0 * kl_1020[k];

        t_1021[k] = f_0 * kl_1021[k];

        t_1022[k] = f_0 * kl_1022[k];

        t_1023[k] = f_0 * kl_1023[k];

        t_1024[k] = f_0 * kl_1024[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, t_1029, t_1030, t_1031, t_1032, \
                         kl_1025, kl_1026, kl_1027, kl_1028, kl_1029, kl_1030, kl_1031, \
                         kl_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_0 * kl_1025[k];

        t_1026[k] = f_0 * kl_1026[k];

        t_1027[k] = f_0 * kl_1027[k];

        t_1028[k] = f_0 * kl_1028[k];

        t_1029[k] = f_0 * kl_1029[k];

        t_1030[k] = f_0 * kl_1030[k];

        t_1031[k] = f_0 * kl_1031[k];

        t_1032[k] = f_0 * kl_1032[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, t_1036, t_1037, t_1038, t_1039, t_1040, \
                         kl_1033, kl_1034, kl_1035, kl_1036, kl_1037, kl_1038, kl_1039, \
                         kl_1040 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = f_0 * kl_1033[k];

        t_1034[k] = f_0 * kl_1034[k];

        t_1035[k] = f_0 * kl_1035[k];

        t_1036[k] = f_0 * kl_1036[k];

        t_1037[k] = f_0 * kl_1037[k];

        t_1038[k] = f_0 * kl_1038[k];

        t_1039[k] = f_0 * kl_1039[k];

        t_1040[k] = f_0 * kl_1040[k];
    }

#pragma omp simd aligned(t_1041, t_1042, t_1043, t_1044, t_1045, t_1046, t_1047, t_1048, \
                         kl_1041, kl_1042, kl_1043, kl_1044, kl_1045, kl_1046, kl_1047, \
                         kl_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1041[k] = f_0 * kl_1041[k];

        t_1042[k] = f_0 * kl_1042[k];

        t_1043[k] = f_0 * kl_1043[k];

        t_1044[k] = f_0 * kl_1044[k];

        t_1045[k] = f_0 * kl_1045[k];

        t_1046[k] = f_0 * kl_1046[k];

        t_1047[k] = f_0 * kl_1047[k];

        t_1048[k] = f_0 * kl_1048[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, t_1052, t_1053, t_1054, t_1055, t_1056, \
                         kl_1049, kl_1050, kl_1051, kl_1052, kl_1053, kl_1054, kl_1055, \
                         kl_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_0 * kl_1049[k];

        t_1050[k] = f_0 * kl_1050[k];

        t_1051[k] = f_0 * kl_1051[k];

        t_1052[k] = f_0 * kl_1052[k];

        t_1053[k] = f_0 * kl_1053[k];

        t_1054[k] = f_0 * kl_1054[k];

        t_1055[k] = f_0 * kl_1055[k];

        t_1056[k] = f_0 * kl_1056[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, t_1060, t_1061, t_1062, t_1063, t_1064, \
                         kl_1057, kl_1058, kl_1059, kl_1060, kl_1061, kl_1062, kl_1063, \
                         kl_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = f_0 * kl_1057[k];

        t_1058[k] = f_0 * kl_1058[k];

        t_1059[k] = f_0 * kl_1059[k];

        t_1060[k] = f_0 * kl_1060[k];

        t_1061[k] = f_0 * kl_1061[k];

        t_1062[k] = f_0 * kl_1062[k];

        t_1063[k] = f_0 * kl_1063[k];

        t_1064[k] = f_0 * kl_1064[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, t_1069, t_1070, t_1071, t_1072, \
                         kl_1065, kl_1066, kl_1067, kl_1068, kl_1069, kl_1070, kl_1071, \
                         kl_1072 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_0 * kl_1065[k];

        t_1066[k] = f_0 * kl_1066[k];

        t_1067[k] = f_0 * kl_1067[k];

        t_1068[k] = f_0 * kl_1068[k];

        t_1069[k] = f_0 * kl_1069[k];

        t_1070[k] = f_0 * kl_1070[k];

        t_1071[k] = f_0 * kl_1071[k];

        t_1072[k] = f_0 * kl_1072[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, t_1076, t_1077, t_1078, t_1079, t_1080, \
                         kl_1073, kl_1074, kl_1075, kl_1076, kl_1077, kl_1078, kl_1079, \
                         kl_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_0 * kl_1073[k];

        t_1074[k] = f_0 * kl_1074[k];

        t_1075[k] = f_0 * kl_1075[k];

        t_1076[k] = f_0 * kl_1076[k];

        t_1077[k] = f_0 * kl_1077[k];

        t_1078[k] = f_0 * kl_1078[k];

        t_1079[k] = f_0 * kl_1079[k];

        t_1080[k] = f_0 * kl_1080[k];
    }

#pragma omp simd aligned(t_1081, t_1082, t_1083, t_1084, t_1085, t_1086, t_1087, t_1088, \
                         kl_1081, kl_1082, kl_1083, kl_1084, kl_1085, kl_1086, kl_1087, \
                         kl_1088 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1081[k] = f_0 * kl_1081[k];

        t_1082[k] = f_0 * kl_1082[k];

        t_1083[k] = f_0 * kl_1083[k];

        t_1084[k] = f_0 * kl_1084[k];

        t_1085[k] = f_0 * kl_1085[k];

        t_1086[k] = f_0 * kl_1086[k];

        t_1087[k] = f_0 * kl_1087[k];

        t_1088[k] = f_0 * kl_1088[k];
    }

#pragma omp simd aligned(t_1089, t_1090, t_1091, t_1092, t_1093, t_1094, t_1095, t_1096, \
                         kl_1089, kl_1090, kl_1091, kl_1092, kl_1093, kl_1094, kl_1095, \
                         kl_1096 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = f_0 * kl_1089[k];

        t_1090[k] = f_0 * kl_1090[k];

        t_1091[k] = f_0 * kl_1091[k];

        t_1092[k] = f_0 * kl_1092[k];

        t_1093[k] = f_0 * kl_1093[k];

        t_1094[k] = f_0 * kl_1094[k];

        t_1095[k] = f_0 * kl_1095[k];

        t_1096[k] = f_0 * kl_1096[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, t_1100, t_1101, t_1102, t_1103, t_1104, \
                         kl_1097, kl_1098, kl_1099, kl_1100, kl_1101, kl_1102, kl_1103, \
                         kl_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_0 * kl_1097[k];

        t_1098[k] = f_0 * kl_1098[k];

        t_1099[k] = f_0 * kl_1099[k];

        t_1100[k] = f_0 * kl_1100[k];

        t_1101[k] = f_0 * kl_1101[k];

        t_1102[k] = f_0 * kl_1102[k];

        t_1103[k] = f_0 * kl_1103[k];

        t_1104[k] = f_0 * kl_1104[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_0_piece7(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kl, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *kl_1105 = buffer.data(kl + 1105);
    const auto *kl_1106 = buffer.data(kl + 1106);
    const auto *kl_1107 = buffer.data(kl + 1107);
    const auto *kl_1108 = buffer.data(kl + 1108);
    const auto *kl_1109 = buffer.data(kl + 1109);
    const auto *kl_1110 = buffer.data(kl + 1110);
    const auto *kl_1111 = buffer.data(kl + 1111);
    const auto *kl_1112 = buffer.data(kl + 1112);
    const auto *kl_1113 = buffer.data(kl + 1113);
    const auto *kl_1114 = buffer.data(kl + 1114);
    const auto *kl_1115 = buffer.data(kl + 1115);
    const auto *kl_1116 = buffer.data(kl + 1116);
    const auto *kl_1117 = buffer.data(kl + 1117);
    const auto *kl_1118 = buffer.data(kl + 1118);
    const auto *kl_1119 = buffer.data(kl + 1119);
    const auto *kl_1120 = buffer.data(kl + 1120);
    const auto *kl_1121 = buffer.data(kl + 1121);
    const auto *kl_1122 = buffer.data(kl + 1122);
    const auto *kl_1123 = buffer.data(kl + 1123);
    const auto *kl_1124 = buffer.data(kl + 1124);
    const auto *kl_1125 = buffer.data(kl + 1125);
    const auto *kl_1126 = buffer.data(kl + 1126);
    const auto *kl_1127 = buffer.data(kl + 1127);
    const auto *kl_1128 = buffer.data(kl + 1128);
    const auto *kl_1129 = buffer.data(kl + 1129);
    const auto *kl_1130 = buffer.data(kl + 1130);
    const auto *kl_1131 = buffer.data(kl + 1131);
    const auto *kl_1132 = buffer.data(kl + 1132);
    const auto *kl_1133 = buffer.data(kl + 1133);
    const auto *kl_1134 = buffer.data(kl + 1134);
    const auto *kl_1135 = buffer.data(kl + 1135);
    const auto *kl_1136 = buffer.data(kl + 1136);
    const auto *kl_1137 = buffer.data(kl + 1137);
    const auto *kl_1138 = buffer.data(kl + 1138);
    const auto *kl_1139 = buffer.data(kl + 1139);
    const auto *kl_1140 = buffer.data(kl + 1140);
    const auto *kl_1141 = buffer.data(kl + 1141);
    const auto *kl_1142 = buffer.data(kl + 1142);
    const auto *kl_1143 = buffer.data(kl + 1143);
    const auto *kl_1144 = buffer.data(kl + 1144);
    const auto *kl_1145 = buffer.data(kl + 1145);
    const auto *kl_1146 = buffer.data(kl + 1146);
    const auto *kl_1147 = buffer.data(kl + 1147);
    const auto *kl_1148 = buffer.data(kl + 1148);
    const auto *kl_1149 = buffer.data(kl + 1149);
    const auto *kl_1150 = buffer.data(kl + 1150);
    const auto *kl_1151 = buffer.data(kl + 1151);
    const auto *kl_1152 = buffer.data(kl + 1152);
    const auto *kl_1153 = buffer.data(kl + 1153);
    const auto *kl_1154 = buffer.data(kl + 1154);
    const auto *kl_1155 = buffer.data(kl + 1155);
    const auto *kl_1156 = buffer.data(kl + 1156);
    const auto *kl_1157 = buffer.data(kl + 1157);
    const auto *kl_1158 = buffer.data(kl + 1158);
    const auto *kl_1159 = buffer.data(kl + 1159);
    const auto *kl_1160 = buffer.data(kl + 1160);
    const auto *kl_1161 = buffer.data(kl + 1161);
    const auto *kl_1162 = buffer.data(kl + 1162);
    const auto *kl_1163 = buffer.data(kl + 1163);
    const auto *kl_1164 = buffer.data(kl + 1164);
    const auto *kl_1165 = buffer.data(kl + 1165);
    const auto *kl_1166 = buffer.data(kl + 1166);
    const auto *kl_1167 = buffer.data(kl + 1167);
    const auto *kl_1168 = buffer.data(kl + 1168);
    const auto *kl_1169 = buffer.data(kl + 1169);
    const auto *kl_1170 = buffer.data(kl + 1170);
    const auto *kl_1171 = buffer.data(kl + 1171);
    const auto *kl_1172 = buffer.data(kl + 1172);
    const auto *kl_1173 = buffer.data(kl + 1173);
    const auto *kl_1174 = buffer.data(kl + 1174);
    const auto *kl_1175 = buffer.data(kl + 1175);
    const auto *kl_1176 = buffer.data(kl + 1176);
    const auto *kl_1177 = buffer.data(kl + 1177);
    const auto *kl_1178 = buffer.data(kl + 1178);
    const auto *kl_1179 = buffer.data(kl + 1179);
    const auto *kl_1180 = buffer.data(kl + 1180);
    const auto *kl_1181 = buffer.data(kl + 1181);
    const auto *kl_1182 = buffer.data(kl + 1182);
    const auto *kl_1183 = buffer.data(kl + 1183);
    const auto *kl_1184 = buffer.data(kl + 1184);
    const auto *kl_1185 = buffer.data(kl + 1185);
    const auto *kl_1186 = buffer.data(kl + 1186);
    const auto *kl_1187 = buffer.data(kl + 1187);
    const auto *kl_1188 = buffer.data(kl + 1188);
    const auto *kl_1189 = buffer.data(kl + 1189);
    const auto *kl_1190 = buffer.data(kl + 1190);
    const auto *kl_1191 = buffer.data(kl + 1191);
    const auto *kl_1192 = buffer.data(kl + 1192);
    const auto *kl_1193 = buffer.data(kl + 1193);
    const auto *kl_1194 = buffer.data(kl + 1194);
    const auto *kl_1195 = buffer.data(kl + 1195);
    const auto *kl_1196 = buffer.data(kl + 1196);
    const auto *kl_1197 = buffer.data(kl + 1197);
    const auto *kl_1198 = buffer.data(kl + 1198);
    const auto *kl_1199 = buffer.data(kl + 1199);
    const auto *kl_1200 = buffer.data(kl + 1200);
    const auto *kl_1201 = buffer.data(kl + 1201);
    const auto *kl_1202 = buffer.data(kl + 1202);
    const auto *kl_1203 = buffer.data(kl + 1203);
    const auto *kl_1204 = buffer.data(kl + 1204);
    const auto *kl_1205 = buffer.data(kl + 1205);
    const auto *kl_1206 = buffer.data(kl + 1206);
    const auto *kl_1207 = buffer.data(kl + 1207);
    const auto *kl_1208 = buffer.data(kl + 1208);
    const auto *kl_1209 = buffer.data(kl + 1209);
    const auto *kl_1210 = buffer.data(kl + 1210);
    const auto *kl_1211 = buffer.data(kl + 1211);
    const auto *kl_1212 = buffer.data(kl + 1212);
    const auto *kl_1213 = buffer.data(kl + 1213);
    const auto *kl_1214 = buffer.data(kl + 1214);
    const auto *kl_1215 = buffer.data(kl + 1215);
    const auto *kl_1216 = buffer.data(kl + 1216);
    const auto *kl_1217 = buffer.data(kl + 1217);
    const auto *kl_1218 = buffer.data(kl + 1218);
    const auto *kl_1219 = buffer.data(kl + 1219);
    const auto *kl_1220 = buffer.data(kl + 1220);
    const auto *kl_1221 = buffer.data(kl + 1221);
    const auto *kl_1222 = buffer.data(kl + 1222);
    const auto *kl_1223 = buffer.data(kl + 1223);
    const auto *kl_1224 = buffer.data(kl + 1224);
    const auto *kl_1225 = buffer.data(kl + 1225);
    const auto *kl_1226 = buffer.data(kl + 1226);
    const auto *kl_1227 = buffer.data(kl + 1227);
    const auto *kl_1228 = buffer.data(kl + 1228);
    const auto *kl_1229 = buffer.data(kl + 1229);
    const auto *kl_1230 = buffer.data(kl + 1230);
    const auto *kl_1231 = buffer.data(kl + 1231);
    const auto *kl_1232 = buffer.data(kl + 1232);
    const auto *kl_1233 = buffer.data(kl + 1233);
    const auto *kl_1234 = buffer.data(kl + 1234);
    const auto *kl_1235 = buffer.data(kl + 1235);
    const auto *kl_1236 = buffer.data(kl + 1236);
    const auto *kl_1237 = buffer.data(kl + 1237);
    const auto *kl_1238 = buffer.data(kl + 1238);
    const auto *kl_1239 = buffer.data(kl + 1239);
    const auto *kl_1240 = buffer.data(kl + 1240);
    const auto *kl_1241 = buffer.data(kl + 1241);
    const auto *kl_1242 = buffer.data(kl + 1242);
    const auto *kl_1243 = buffer.data(kl + 1243);
    const auto *kl_1244 = buffer.data(kl + 1244);
    const auto *kl_1245 = buffer.data(kl + 1245);
    const auto *kl_1246 = buffer.data(kl + 1246);
    const auto *kl_1247 = buffer.data(kl + 1247);
    const auto *kl_1248 = buffer.data(kl + 1248);
    const auto *kl_1249 = buffer.data(kl + 1249);
    const auto *kl_1250 = buffer.data(kl + 1250);
    const auto *kl_1251 = buffer.data(kl + 1251);
    const auto *kl_1252 = buffer.data(kl + 1252);
    const auto *kl_1253 = buffer.data(kl + 1253);
    const auto *kl_1254 = buffer.data(kl + 1254);
    const auto *kl_1255 = buffer.data(kl + 1255);
    const auto *kl_1256 = buffer.data(kl + 1256);
    const auto *kl_1257 = buffer.data(kl + 1257);
    const auto *kl_1258 = buffer.data(kl + 1258);
    const auto *kl_1259 = buffer.data(kl + 1259);

#pragma omp simd aligned(t_1105, t_1106, t_1107, t_1108, t_1109, t_1110, t_1111, t_1112, \
                         kl_1105, kl_1106, kl_1107, kl_1108, kl_1109, kl_1110, kl_1111, \
                         kl_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1105[k] = f_0 * kl_1105[k];

        t_1106[k] = f_0 * kl_1106[k];

        t_1107[k] = f_0 * kl_1107[k];

        t_1108[k] = f_0 * kl_1108[k];

        t_1109[k] = f_0 * kl_1109[k];

        t_1110[k] = f_0 * kl_1110[k];

        t_1111[k] = f_0 * kl_1111[k];

        t_1112[k] = f_0 * kl_1112[k];
    }

#pragma omp simd aligned(t_1113, t_1114, t_1115, t_1116, t_1117, t_1118, t_1119, t_1120, \
                         kl_1113, kl_1114, kl_1115, kl_1116, kl_1117, kl_1118, kl_1119, \
                         kl_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1113[k] = f_0 * kl_1113[k];

        t_1114[k] = f_0 * kl_1114[k];

        t_1115[k] = f_0 * kl_1115[k];

        t_1116[k] = f_0 * kl_1116[k];

        t_1117[k] = f_0 * kl_1117[k];

        t_1118[k] = f_0 * kl_1118[k];

        t_1119[k] = f_0 * kl_1119[k];

        t_1120[k] = f_0 * kl_1120[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, t_1124, t_1125, t_1126, t_1127, t_1128, \
                         kl_1121, kl_1122, kl_1123, kl_1124, kl_1125, kl_1126, kl_1127, \
                         kl_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_0 * kl_1121[k];

        t_1122[k] = f_0 * kl_1122[k];

        t_1123[k] = f_0 * kl_1123[k];

        t_1124[k] = f_0 * kl_1124[k];

        t_1125[k] = f_0 * kl_1125[k];

        t_1126[k] = f_0 * kl_1126[k];

        t_1127[k] = f_0 * kl_1127[k];

        t_1128[k] = f_0 * kl_1128[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, t_1132, t_1133, t_1134, t_1135, t_1136, \
                         kl_1129, kl_1130, kl_1131, kl_1132, kl_1133, kl_1134, kl_1135, \
                         kl_1136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_0 * kl_1129[k];

        t_1130[k] = f_0 * kl_1130[k];

        t_1131[k] = f_0 * kl_1131[k];

        t_1132[k] = f_0 * kl_1132[k];

        t_1133[k] = f_0 * kl_1133[k];

        t_1134[k] = f_0 * kl_1134[k];

        t_1135[k] = f_0 * kl_1135[k];

        t_1136[k] = f_0 * kl_1136[k];
    }

#pragma omp simd aligned(t_1137, t_1138, t_1139, t_1140, t_1141, t_1142, t_1143, t_1144, \
                         kl_1137, kl_1138, kl_1139, kl_1140, kl_1141, kl_1142, kl_1143, \
                         kl_1144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1137[k] = f_0 * kl_1137[k];

        t_1138[k] = f_0 * kl_1138[k];

        t_1139[k] = f_0 * kl_1139[k];

        t_1140[k] = f_0 * kl_1140[k];

        t_1141[k] = f_0 * kl_1141[k];

        t_1142[k] = f_0 * kl_1142[k];

        t_1143[k] = f_0 * kl_1143[k];

        t_1144[k] = f_0 * kl_1144[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, t_1148, t_1149, t_1150, t_1151, t_1152, \
                         kl_1145, kl_1146, kl_1147, kl_1148, kl_1149, kl_1150, kl_1151, \
                         kl_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = f_0 * kl_1145[k];

        t_1146[k] = f_0 * kl_1146[k];

        t_1147[k] = f_0 * kl_1147[k];

        t_1148[k] = f_0 * kl_1148[k];

        t_1149[k] = f_0 * kl_1149[k];

        t_1150[k] = f_0 * kl_1150[k];

        t_1151[k] = f_0 * kl_1151[k];

        t_1152[k] = f_0 * kl_1152[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, t_1157, t_1158, t_1159, t_1160, \
                         kl_1153, kl_1154, kl_1155, kl_1156, kl_1157, kl_1158, kl_1159, \
                         kl_1160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_0 * kl_1153[k];

        t_1154[k] = f_0 * kl_1154[k];

        t_1155[k] = f_0 * kl_1155[k];

        t_1156[k] = f_0 * kl_1156[k];

        t_1157[k] = f_0 * kl_1157[k];

        t_1158[k] = f_0 * kl_1158[k];

        t_1159[k] = f_0 * kl_1159[k];

        t_1160[k] = f_0 * kl_1160[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, t_1164, t_1165, t_1166, t_1167, t_1168, \
                         kl_1161, kl_1162, kl_1163, kl_1164, kl_1165, kl_1166, kl_1167, \
                         kl_1168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = f_0 * kl_1161[k];

        t_1162[k] = f_0 * kl_1162[k];

        t_1163[k] = f_0 * kl_1163[k];

        t_1164[k] = f_0 * kl_1164[k];

        t_1165[k] = f_0 * kl_1165[k];

        t_1166[k] = f_0 * kl_1166[k];

        t_1167[k] = f_0 * kl_1167[k];

        t_1168[k] = f_0 * kl_1168[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, t_1172, t_1173, t_1174, t_1175, t_1176, \
                         kl_1169, kl_1170, kl_1171, kl_1172, kl_1173, kl_1174, kl_1175, \
                         kl_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_0 * kl_1169[k];

        t_1170[k] = f_0 * kl_1170[k];

        t_1171[k] = f_0 * kl_1171[k];

        t_1172[k] = f_0 * kl_1172[k];

        t_1173[k] = f_0 * kl_1173[k];

        t_1174[k] = f_0 * kl_1174[k];

        t_1175[k] = f_0 * kl_1175[k];

        t_1176[k] = f_0 * kl_1176[k];
    }

#pragma omp simd aligned(t_1177, t_1178, t_1179, t_1180, t_1181, t_1182, t_1183, t_1184, \
                         kl_1177, kl_1178, kl_1179, kl_1180, kl_1181, kl_1182, kl_1183, \
                         kl_1184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1177[k] = f_0 * kl_1177[k];

        t_1178[k] = f_0 * kl_1178[k];

        t_1179[k] = f_0 * kl_1179[k];

        t_1180[k] = f_0 * kl_1180[k];

        t_1181[k] = f_0 * kl_1181[k];

        t_1182[k] = f_0 * kl_1182[k];

        t_1183[k] = f_0 * kl_1183[k];

        t_1184[k] = f_0 * kl_1184[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, t_1188, t_1189, t_1190, t_1191, t_1192, \
                         kl_1185, kl_1186, kl_1187, kl_1188, kl_1189, kl_1190, kl_1191, \
                         kl_1192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_0 * kl_1185[k];

        t_1186[k] = f_0 * kl_1186[k];

        t_1187[k] = f_0 * kl_1187[k];

        t_1188[k] = f_0 * kl_1188[k];

        t_1189[k] = f_0 * kl_1189[k];

        t_1190[k] = f_0 * kl_1190[k];

        t_1191[k] = f_0 * kl_1191[k];

        t_1192[k] = f_0 * kl_1192[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, t_1196, t_1197, t_1198, t_1199, t_1200, \
                         kl_1193, kl_1194, kl_1195, kl_1196, kl_1197, kl_1198, kl_1199, \
                         kl_1200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_0 * kl_1193[k];

        t_1194[k] = f_0 * kl_1194[k];

        t_1195[k] = f_0 * kl_1195[k];

        t_1196[k] = f_0 * kl_1196[k];

        t_1197[k] = f_0 * kl_1197[k];

        t_1198[k] = f_0 * kl_1198[k];

        t_1199[k] = f_0 * kl_1199[k];

        t_1200[k] = f_0 * kl_1200[k];
    }

#pragma omp simd aligned(t_1201, t_1202, t_1203, t_1204, t_1205, t_1206, t_1207, t_1208, \
                         kl_1201, kl_1202, kl_1203, kl_1204, kl_1205, kl_1206, kl_1207, \
                         kl_1208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1201[k] = f_0 * kl_1201[k];

        t_1202[k] = f_0 * kl_1202[k];

        t_1203[k] = f_0 * kl_1203[k];

        t_1204[k] = f_0 * kl_1204[k];

        t_1205[k] = f_0 * kl_1205[k];

        t_1206[k] = f_0 * kl_1206[k];

        t_1207[k] = f_0 * kl_1207[k];

        t_1208[k] = f_0 * kl_1208[k];
    }

#pragma omp simd aligned(t_1209, t_1210, t_1211, t_1212, t_1213, t_1214, t_1215, t_1216, \
                         kl_1209, kl_1210, kl_1211, kl_1212, kl_1213, kl_1214, kl_1215, \
                         kl_1216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1209[k] = f_0 * kl_1209[k];

        t_1210[k] = f_0 * kl_1210[k];

        t_1211[k] = f_0 * kl_1211[k];

        t_1212[k] = f_0 * kl_1212[k];

        t_1213[k] = f_0 * kl_1213[k];

        t_1214[k] = f_0 * kl_1214[k];

        t_1215[k] = f_0 * kl_1215[k];

        t_1216[k] = f_0 * kl_1216[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, t_1220, t_1221, t_1222, t_1223, t_1224, \
                         kl_1217, kl_1218, kl_1219, kl_1220, kl_1221, kl_1222, kl_1223, \
                         kl_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = f_0 * kl_1217[k];

        t_1218[k] = f_0 * kl_1218[k];

        t_1219[k] = f_0 * kl_1219[k];

        t_1220[k] = f_0 * kl_1220[k];

        t_1221[k] = f_0 * kl_1221[k];

        t_1222[k] = f_0 * kl_1222[k];

        t_1223[k] = f_0 * kl_1223[k];

        t_1224[k] = f_0 * kl_1224[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, t_1228, t_1229, t_1230, t_1231, t_1232, \
                         kl_1225, kl_1226, kl_1227, kl_1228, kl_1229, kl_1230, kl_1231, \
                         kl_1232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = f_0 * kl_1225[k];

        t_1226[k] = f_0 * kl_1226[k];

        t_1227[k] = f_0 * kl_1227[k];

        t_1228[k] = f_0 * kl_1228[k];

        t_1229[k] = f_0 * kl_1229[k];

        t_1230[k] = f_0 * kl_1230[k];

        t_1231[k] = f_0 * kl_1231[k];

        t_1232[k] = f_0 * kl_1232[k];
    }

#pragma omp simd aligned(t_1233, t_1234, t_1235, t_1236, t_1237, t_1238, t_1239, t_1240, \
                         kl_1233, kl_1234, kl_1235, kl_1236, kl_1237, kl_1238, kl_1239, \
                         kl_1240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1233[k] = f_0 * kl_1233[k];

        t_1234[k] = f_0 * kl_1234[k];

        t_1235[k] = f_0 * kl_1235[k];

        t_1236[k] = f_0 * kl_1236[k];

        t_1237[k] = f_0 * kl_1237[k];

        t_1238[k] = f_0 * kl_1238[k];

        t_1239[k] = f_0 * kl_1239[k];

        t_1240[k] = f_0 * kl_1240[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, t_1244, t_1245, t_1246, t_1247, t_1248, \
                         kl_1241, kl_1242, kl_1243, kl_1244, kl_1245, kl_1246, kl_1247, \
                         kl_1248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_0 * kl_1241[k];

        t_1242[k] = f_0 * kl_1242[k];

        t_1243[k] = f_0 * kl_1243[k];

        t_1244[k] = f_0 * kl_1244[k];

        t_1245[k] = f_0 * kl_1245[k];

        t_1246[k] = f_0 * kl_1246[k];

        t_1247[k] = f_0 * kl_1247[k];

        t_1248[k] = f_0 * kl_1248[k];
    }

#pragma omp simd aligned(t_1249, t_1250, t_1251, t_1252, t_1253, t_1254, t_1255, t_1256, \
                         kl_1249, kl_1250, kl_1251, kl_1252, kl_1253, kl_1254, kl_1255, \
                         kl_1256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1249[k] = f_0 * kl_1249[k];

        t_1250[k] = f_0 * kl_1250[k];

        t_1251[k] = f_0 * kl_1251[k];

        t_1252[k] = f_0 * kl_1252[k];

        t_1253[k] = f_0 * kl_1253[k];

        t_1254[k] = f_0 * kl_1254[k];

        t_1255[k] = f_0 * kl_1255[k];

        t_1256[k] = f_0 * kl_1256[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, kl_1257, kl_1258, \
                         kl_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_0 * kl_1257[k];

        t_1258[k] = f_0 * kl_1258[k];

        t_1259[k] = f_0 * kl_1259[k];
    }
}

auto
compute_prim_geom_10_il_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t hl, const size_t kl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_il_electron_repulsion_0_piece0(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_0_piece1(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_0_piece2(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_0_piece3(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_0_piece4(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_0_piece5(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_0_piece6(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_0_piece7(buffer, target, kl, ncols, alpha);
}

static auto
compute_prim_geom_10_il_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
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

    const auto *hl_0 = buffer.data(hl + 0);
    const auto *hl_1 = buffer.data(hl + 1);
    const auto *hl_2 = buffer.data(hl + 2);
    const auto *hl_3 = buffer.data(hl + 3);
    const auto *hl_4 = buffer.data(hl + 4);
    const auto *hl_5 = buffer.data(hl + 5);
    const auto *hl_6 = buffer.data(hl + 6);
    const auto *hl_7 = buffer.data(hl + 7);
    const auto *hl_8 = buffer.data(hl + 8);
    const auto *hl_9 = buffer.data(hl + 9);
    const auto *hl_10 = buffer.data(hl + 10);
    const auto *hl_11 = buffer.data(hl + 11);
    const auto *hl_12 = buffer.data(hl + 12);
    const auto *hl_13 = buffer.data(hl + 13);
    const auto *hl_14 = buffer.data(hl + 14);
    const auto *hl_15 = buffer.data(hl + 15);
    const auto *hl_16 = buffer.data(hl + 16);
    const auto *hl_17 = buffer.data(hl + 17);
    const auto *hl_18 = buffer.data(hl + 18);
    const auto *hl_19 = buffer.data(hl + 19);
    const auto *hl_20 = buffer.data(hl + 20);
    const auto *hl_21 = buffer.data(hl + 21);
    const auto *hl_22 = buffer.data(hl + 22);
    const auto *hl_23 = buffer.data(hl + 23);
    const auto *hl_24 = buffer.data(hl + 24);
    const auto *hl_25 = buffer.data(hl + 25);
    const auto *hl_26 = buffer.data(hl + 26);
    const auto *hl_27 = buffer.data(hl + 27);
    const auto *hl_28 = buffer.data(hl + 28);
    const auto *hl_29 = buffer.data(hl + 29);
    const auto *hl_30 = buffer.data(hl + 30);
    const auto *hl_31 = buffer.data(hl + 31);
    const auto *hl_32 = buffer.data(hl + 32);
    const auto *hl_33 = buffer.data(hl + 33);
    const auto *hl_34 = buffer.data(hl + 34);
    const auto *hl_35 = buffer.data(hl + 35);
    const auto *hl_36 = buffer.data(hl + 36);
    const auto *hl_37 = buffer.data(hl + 37);
    const auto *hl_38 = buffer.data(hl + 38);
    const auto *hl_39 = buffer.data(hl + 39);
    const auto *hl_40 = buffer.data(hl + 40);
    const auto *hl_41 = buffer.data(hl + 41);
    const auto *hl_42 = buffer.data(hl + 42);
    const auto *hl_43 = buffer.data(hl + 43);
    const auto *hl_44 = buffer.data(hl + 44);
    const auto *hl_45 = buffer.data(hl + 45);
    const auto *hl_46 = buffer.data(hl + 46);
    const auto *hl_47 = buffer.data(hl + 47);
    const auto *hl_48 = buffer.data(hl + 48);
    const auto *hl_49 = buffer.data(hl + 49);
    const auto *hl_50 = buffer.data(hl + 50);
    const auto *hl_51 = buffer.data(hl + 51);
    const auto *hl_52 = buffer.data(hl + 52);
    const auto *hl_53 = buffer.data(hl + 53);
    const auto *hl_54 = buffer.data(hl + 54);
    const auto *hl_55 = buffer.data(hl + 55);
    const auto *hl_56 = buffer.data(hl + 56);
    const auto *hl_57 = buffer.data(hl + 57);
    const auto *hl_58 = buffer.data(hl + 58);
    const auto *hl_59 = buffer.data(hl + 59);
    const auto *hl_60 = buffer.data(hl + 60);
    const auto *hl_61 = buffer.data(hl + 61);
    const auto *hl_62 = buffer.data(hl + 62);
    const auto *hl_63 = buffer.data(hl + 63);
    const auto *hl_64 = buffer.data(hl + 64);
    const auto *hl_65 = buffer.data(hl + 65);
    const auto *hl_66 = buffer.data(hl + 66);
    const auto *hl_67 = buffer.data(hl + 67);
    const auto *hl_68 = buffer.data(hl + 68);
    const auto *hl_69 = buffer.data(hl + 69);
    const auto *hl_70 = buffer.data(hl + 70);
    const auto *hl_71 = buffer.data(hl + 71);
    const auto *hl_72 = buffer.data(hl + 72);
    const auto *hl_73 = buffer.data(hl + 73);
    const auto *hl_74 = buffer.data(hl + 74);
    const auto *hl_75 = buffer.data(hl + 75);
    const auto *hl_76 = buffer.data(hl + 76);
    const auto *hl_77 = buffer.data(hl + 77);
    const auto *hl_78 = buffer.data(hl + 78);
    const auto *hl_79 = buffer.data(hl + 79);
    const auto *hl_80 = buffer.data(hl + 80);
    const auto *hl_81 = buffer.data(hl + 81);
    const auto *hl_82 = buffer.data(hl + 82);
    const auto *hl_83 = buffer.data(hl + 83);
    const auto *hl_84 = buffer.data(hl + 84);
    const auto *hl_85 = buffer.data(hl + 85);
    const auto *hl_86 = buffer.data(hl + 86);
    const auto *hl_87 = buffer.data(hl + 87);
    const auto *hl_88 = buffer.data(hl + 88);

    const auto *kl_45 = buffer.data(kl + 45);
    const auto *kl_46 = buffer.data(kl + 46);
    const auto *kl_47 = buffer.data(kl + 47);
    const auto *kl_48 = buffer.data(kl + 48);
    const auto *kl_49 = buffer.data(kl + 49);
    const auto *kl_50 = buffer.data(kl + 50);
    const auto *kl_51 = buffer.data(kl + 51);
    const auto *kl_52 = buffer.data(kl + 52);
    const auto *kl_53 = buffer.data(kl + 53);
    const auto *kl_54 = buffer.data(kl + 54);
    const auto *kl_55 = buffer.data(kl + 55);
    const auto *kl_56 = buffer.data(kl + 56);
    const auto *kl_57 = buffer.data(kl + 57);
    const auto *kl_58 = buffer.data(kl + 58);
    const auto *kl_59 = buffer.data(kl + 59);
    const auto *kl_60 = buffer.data(kl + 60);
    const auto *kl_61 = buffer.data(kl + 61);
    const auto *kl_62 = buffer.data(kl + 62);
    const auto *kl_63 = buffer.data(kl + 63);
    const auto *kl_64 = buffer.data(kl + 64);
    const auto *kl_65 = buffer.data(kl + 65);
    const auto *kl_66 = buffer.data(kl + 66);
    const auto *kl_67 = buffer.data(kl + 67);
    const auto *kl_68 = buffer.data(kl + 68);
    const auto *kl_69 = buffer.data(kl + 69);
    const auto *kl_70 = buffer.data(kl + 70);
    const auto *kl_71 = buffer.data(kl + 71);
    const auto *kl_72 = buffer.data(kl + 72);
    const auto *kl_73 = buffer.data(kl + 73);
    const auto *kl_74 = buffer.data(kl + 74);
    const auto *kl_75 = buffer.data(kl + 75);
    const auto *kl_76 = buffer.data(kl + 76);
    const auto *kl_77 = buffer.data(kl + 77);
    const auto *kl_78 = buffer.data(kl + 78);
    const auto *kl_79 = buffer.data(kl + 79);
    const auto *kl_80 = buffer.data(kl + 80);
    const auto *kl_81 = buffer.data(kl + 81);
    const auto *kl_82 = buffer.data(kl + 82);
    const auto *kl_83 = buffer.data(kl + 83);
    const auto *kl_84 = buffer.data(kl + 84);
    const auto *kl_85 = buffer.data(kl + 85);
    const auto *kl_86 = buffer.data(kl + 86);
    const auto *kl_87 = buffer.data(kl + 87);
    const auto *kl_88 = buffer.data(kl + 88);
    const auto *kl_89 = buffer.data(kl + 89);
    const auto *kl_135 = buffer.data(kl + 135);
    const auto *kl_136 = buffer.data(kl + 136);
    const auto *kl_137 = buffer.data(kl + 137);
    const auto *kl_138 = buffer.data(kl + 138);
    const auto *kl_139 = buffer.data(kl + 139);
    const auto *kl_140 = buffer.data(kl + 140);
    const auto *kl_141 = buffer.data(kl + 141);
    const auto *kl_142 = buffer.data(kl + 142);
    const auto *kl_143 = buffer.data(kl + 143);
    const auto *kl_144 = buffer.data(kl + 144);
    const auto *kl_145 = buffer.data(kl + 145);
    const auto *kl_146 = buffer.data(kl + 146);
    const auto *kl_147 = buffer.data(kl + 147);
    const auto *kl_148 = buffer.data(kl + 148);
    const auto *kl_149 = buffer.data(kl + 149);
    const auto *kl_150 = buffer.data(kl + 150);
    const auto *kl_151 = buffer.data(kl + 151);
    const auto *kl_152 = buffer.data(kl + 152);
    const auto *kl_153 = buffer.data(kl + 153);
    const auto *kl_154 = buffer.data(kl + 154);
    const auto *kl_155 = buffer.data(kl + 155);
    const auto *kl_156 = buffer.data(kl + 156);
    const auto *kl_157 = buffer.data(kl + 157);
    const auto *kl_158 = buffer.data(kl + 158);
    const auto *kl_159 = buffer.data(kl + 159);
    const auto *kl_160 = buffer.data(kl + 160);
    const auto *kl_161 = buffer.data(kl + 161);
    const auto *kl_162 = buffer.data(kl + 162);
    const auto *kl_163 = buffer.data(kl + 163);
    const auto *kl_164 = buffer.data(kl + 164);
    const auto *kl_165 = buffer.data(kl + 165);
    const auto *kl_166 = buffer.data(kl + 166);
    const auto *kl_167 = buffer.data(kl + 167);
    const auto *kl_168 = buffer.data(kl + 168);
    const auto *kl_169 = buffer.data(kl + 169);
    const auto *kl_170 = buffer.data(kl + 170);
    const auto *kl_171 = buffer.data(kl + 171);
    const auto *kl_172 = buffer.data(kl + 172);
    const auto *kl_173 = buffer.data(kl + 173);
    const auto *kl_174 = buffer.data(kl + 174);
    const auto *kl_175 = buffer.data(kl + 175);
    const auto *kl_176 = buffer.data(kl + 176);
    const auto *kl_177 = buffer.data(kl + 177);
    const auto *kl_178 = buffer.data(kl + 178);
    const auto *kl_179 = buffer.data(kl + 179);
    const auto *kl_180 = buffer.data(kl + 180);
    const auto *kl_181 = buffer.data(kl + 181);
    const auto *kl_182 = buffer.data(kl + 182);
    const auto *kl_183 = buffer.data(kl + 183);
    const auto *kl_184 = buffer.data(kl + 184);
    const auto *kl_185 = buffer.data(kl + 185);
    const auto *kl_186 = buffer.data(kl + 186);
    const auto *kl_187 = buffer.data(kl + 187);
    const auto *kl_188 = buffer.data(kl + 188);
    const auto *kl_189 = buffer.data(kl + 189);
    const auto *kl_190 = buffer.data(kl + 190);
    const auto *kl_191 = buffer.data(kl + 191);
    const auto *kl_192 = buffer.data(kl + 192);
    const auto *kl_193 = buffer.data(kl + 193);
    const auto *kl_194 = buffer.data(kl + 194);
    const auto *kl_195 = buffer.data(kl + 195);
    const auto *kl_196 = buffer.data(kl + 196);
    const auto *kl_197 = buffer.data(kl + 197);
    const auto *kl_198 = buffer.data(kl + 198);
    const auto *kl_199 = buffer.data(kl + 199);
    const auto *kl_200 = buffer.data(kl + 200);
    const auto *kl_201 = buffer.data(kl + 201);
    const auto *kl_202 = buffer.data(kl + 202);
    const auto *kl_203 = buffer.data(kl + 203);
    const auto *kl_204 = buffer.data(kl + 204);
    const auto *kl_205 = buffer.data(kl + 205);
    const auto *kl_206 = buffer.data(kl + 206);
    const auto *kl_207 = buffer.data(kl + 207);
    const auto *kl_208 = buffer.data(kl + 208);
    const auto *kl_209 = buffer.data(kl + 209);
    const auto *kl_210 = buffer.data(kl + 210);
    const auto *kl_211 = buffer.data(kl + 211);
    const auto *kl_212 = buffer.data(kl + 212);
    const auto *kl_213 = buffer.data(kl + 213);
    const auto *kl_214 = buffer.data(kl + 214);
    const auto *kl_215 = buffer.data(kl + 215);
    const auto *kl_216 = buffer.data(kl + 216);
    const auto *kl_217 = buffer.data(kl + 217);
    const auto *kl_218 = buffer.data(kl + 218);
    const auto *kl_219 = buffer.data(kl + 219);
    const auto *kl_220 = buffer.data(kl + 220);
    const auto *kl_221 = buffer.data(kl + 221);
    const auto *kl_222 = buffer.data(kl + 222);
    const auto *kl_223 = buffer.data(kl + 223);
    const auto *kl_224 = buffer.data(kl + 224);
    const auto *kl_270 = buffer.data(kl + 270);
    const auto *kl_271 = buffer.data(kl + 271);
    const auto *kl_272 = buffer.data(kl + 272);
    const auto *kl_273 = buffer.data(kl + 273);
    const auto *kl_274 = buffer.data(kl + 274);
    const auto *kl_275 = buffer.data(kl + 275);
    const auto *kl_276 = buffer.data(kl + 276);
    const auto *kl_277 = buffer.data(kl + 277);
    const auto *kl_278 = buffer.data(kl + 278);
    const auto *kl_279 = buffer.data(kl + 279);
    const auto *kl_280 = buffer.data(kl + 280);
    const auto *kl_281 = buffer.data(kl + 281);
    const auto *kl_282 = buffer.data(kl + 282);
    const auto *kl_283 = buffer.data(kl + 283);
    const auto *kl_284 = buffer.data(kl + 284);
    const auto *kl_285 = buffer.data(kl + 285);
    const auto *kl_286 = buffer.data(kl + 286);
    const auto *kl_287 = buffer.data(kl + 287);
    const auto *kl_288 = buffer.data(kl + 288);
    const auto *kl_289 = buffer.data(kl + 289);
    const auto *kl_290 = buffer.data(kl + 290);
    const auto *kl_291 = buffer.data(kl + 291);
    const auto *kl_292 = buffer.data(kl + 292);
    const auto *kl_293 = buffer.data(kl + 293);
    const auto *kl_294 = buffer.data(kl + 294);
    const auto *kl_295 = buffer.data(kl + 295);
    const auto *kl_296 = buffer.data(kl + 296);
    const auto *kl_297 = buffer.data(kl + 297);
    const auto *kl_298 = buffer.data(kl + 298);
    const auto *kl_299 = buffer.data(kl + 299);
    const auto *kl_300 = buffer.data(kl + 300);
    const auto *kl_301 = buffer.data(kl + 301);
    const auto *kl_302 = buffer.data(kl + 302);
    const auto *kl_303 = buffer.data(kl + 303);
    const auto *kl_304 = buffer.data(kl + 304);
    const auto *kl_305 = buffer.data(kl + 305);
    const auto *kl_306 = buffer.data(kl + 306);
    const auto *kl_307 = buffer.data(kl + 307);
    const auto *kl_308 = buffer.data(kl + 308);
    const auto *kl_309 = buffer.data(kl + 309);
    const auto *kl_310 = buffer.data(kl + 310);
    const auto *kl_311 = buffer.data(kl + 311);
    const auto *kl_312 = buffer.data(kl + 312);
    const auto *kl_313 = buffer.data(kl + 313);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, kl_45, kl_46, kl_47, kl_48, \
                         kl_49, kl_50, kl_51, kl_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kl_45[k];

        t_1[k] = f_0 * kl_46[k];

        t_2[k] = f_0 * kl_47[k];

        t_3[k] = f_0 * kl_48[k];

        t_4[k] = f_0 * kl_49[k];

        t_5[k] = f_0 * kl_50[k];

        t_6[k] = f_0 * kl_51[k];

        t_7[k] = f_0 * kl_52[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, kl_53, kl_54, kl_55, \
                         kl_56, kl_57, kl_58, kl_59, kl_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * kl_53[k];

        t_9[k] = f_0 * kl_54[k];

        t_10[k] = f_0 * kl_55[k];

        t_11[k] = f_0 * kl_56[k];

        t_12[k] = f_0 * kl_57[k];

        t_13[k] = f_0 * kl_58[k];

        t_14[k] = f_0 * kl_59[k];

        t_15[k] = f_0 * kl_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, kl_61, kl_62, kl_63, \
                         kl_64, kl_65, kl_66, kl_67, kl_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * kl_61[k];

        t_17[k] = f_0 * kl_62[k];

        t_18[k] = f_0 * kl_63[k];

        t_19[k] = f_0 * kl_64[k];

        t_20[k] = f_0 * kl_65[k];

        t_21[k] = f_0 * kl_66[k];

        t_22[k] = f_0 * kl_67[k];

        t_23[k] = f_0 * kl_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, kl_69, kl_70, kl_71, \
                         kl_72, kl_73, kl_74, kl_75, kl_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * kl_69[k];

        t_25[k] = f_0 * kl_70[k];

        t_26[k] = f_0 * kl_71[k];

        t_27[k] = f_0 * kl_72[k];

        t_28[k] = f_0 * kl_73[k];

        t_29[k] = f_0 * kl_74[k];

        t_30[k] = f_0 * kl_75[k];

        t_31[k] = f_0 * kl_76[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, kl_77, kl_78, kl_79, \
                         kl_80, kl_81, kl_82, kl_83, kl_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * kl_77[k];

        t_33[k] = f_0 * kl_78[k];

        t_34[k] = f_0 * kl_79[k];

        t_35[k] = f_0 * kl_80[k];

        t_36[k] = f_0 * kl_81[k];

        t_37[k] = f_0 * kl_82[k];

        t_38[k] = f_0 * kl_83[k];

        t_39[k] = f_0 * kl_84[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, hl_0, hl_1, kl_85, kl_86, \
                         kl_87, kl_88, kl_89, kl_135, kl_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * kl_85[k];

        t_41[k] = f_0 * kl_86[k];

        t_42[k] = f_0 * kl_87[k];

        t_43[k] = f_0 * kl_88[k];

        t_44[k] = f_0 * kl_89[k];

        t_45[k] = -hl_0[k]
                  + f_0 * kl_135[k];

        t_46[k] = -hl_1[k]
                  + f_0 * kl_136[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, hl_2, hl_3, hl_4, hl_5, hl_6, kl_137, \
                         kl_138, kl_139, kl_140, kl_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -hl_2[k]
                  + f_0 * kl_137[k];

        t_48[k] = -hl_3[k]
                  + f_0 * kl_138[k];

        t_49[k] = -hl_4[k]
                  + f_0 * kl_139[k];

        t_50[k] = -hl_5[k]
                  + f_0 * kl_140[k];

        t_51[k] = -hl_6[k]
                  + f_0 * kl_141[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, hl_7, hl_8, hl_9, hl_10, hl_11, kl_142, \
                         kl_143, kl_144, kl_145, kl_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = -hl_7[k]
                  + f_0 * kl_142[k];

        t_53[k] = -hl_8[k]
                  + f_0 * kl_143[k];

        t_54[k] = -hl_9[k]
                  + f_0 * kl_144[k];

        t_55[k] = -hl_10[k]
                  + f_0 * kl_145[k];

        t_56[k] = -hl_11[k]
                  + f_0 * kl_146[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, hl_12, hl_13, hl_14, hl_15, hl_16, \
                         kl_147, kl_148, kl_149, kl_150, kl_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -hl_12[k]
                  + f_0 * kl_147[k];

        t_58[k] = -hl_13[k]
                  + f_0 * kl_148[k];

        t_59[k] = -hl_14[k]
                  + f_0 * kl_149[k];

        t_60[k] = -hl_15[k]
                  + f_0 * kl_150[k];

        t_61[k] = -hl_16[k]
                  + f_0 * kl_151[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, hl_17, hl_18, hl_19, hl_20, hl_21, \
                         kl_152, kl_153, kl_154, kl_155, kl_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -hl_17[k]
                  + f_0 * kl_152[k];

        t_63[k] = -hl_18[k]
                  + f_0 * kl_153[k];

        t_64[k] = -hl_19[k]
                  + f_0 * kl_154[k];

        t_65[k] = -hl_20[k]
                  + f_0 * kl_155[k];

        t_66[k] = -hl_21[k]
                  + f_0 * kl_156[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, hl_22, hl_23, hl_24, hl_25, hl_26, \
                         kl_157, kl_158, kl_159, kl_160, kl_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -hl_22[k]
                  + f_0 * kl_157[k];

        t_68[k] = -hl_23[k]
                  + f_0 * kl_158[k];

        t_69[k] = -hl_24[k]
                  + f_0 * kl_159[k];

        t_70[k] = -hl_25[k]
                  + f_0 * kl_160[k];

        t_71[k] = -hl_26[k]
                  + f_0 * kl_161[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, hl_27, hl_28, hl_29, hl_30, hl_31, \
                         kl_162, kl_163, kl_164, kl_165, kl_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -hl_27[k]
                  + f_0 * kl_162[k];

        t_73[k] = -hl_28[k]
                  + f_0 * kl_163[k];

        t_74[k] = -hl_29[k]
                  + f_0 * kl_164[k];

        t_75[k] = -hl_30[k]
                  + f_0 * kl_165[k];

        t_76[k] = -hl_31[k]
                  + f_0 * kl_166[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, hl_32, hl_33, hl_34, hl_35, hl_36, \
                         kl_167, kl_168, kl_169, kl_170, kl_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -hl_32[k]
                  + f_0 * kl_167[k];

        t_78[k] = -hl_33[k]
                  + f_0 * kl_168[k];

        t_79[k] = -hl_34[k]
                  + f_0 * kl_169[k];

        t_80[k] = -hl_35[k]
                  + f_0 * kl_170[k];

        t_81[k] = -hl_36[k]
                  + f_0 * kl_171[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, hl_37, hl_38, hl_39, hl_40, hl_41, \
                         kl_172, kl_173, kl_174, kl_175, kl_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -hl_37[k]
                  + f_0 * kl_172[k];

        t_83[k] = -hl_38[k]
                  + f_0 * kl_173[k];

        t_84[k] = -hl_39[k]
                  + f_0 * kl_174[k];

        t_85[k] = -hl_40[k]
                  + f_0 * kl_175[k];

        t_86[k] = -hl_41[k]
                  + f_0 * kl_176[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, hl_42, hl_43, hl_44, kl_177, \
                         kl_178, kl_179, kl_180, kl_181, kl_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -hl_42[k]
                  + f_0 * kl_177[k];

        t_88[k] = -hl_43[k]
                  + f_0 * kl_178[k];

        t_89[k] = -hl_44[k]
                  + f_0 * kl_179[k];

        t_90[k] = f_0 * kl_180[k];

        t_91[k] = f_0 * kl_181[k];

        t_92[k] = f_0 * kl_182[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, kl_183, kl_184, \
                         kl_185, kl_186, kl_187, kl_188, kl_189, \
                         kl_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * kl_183[k];

        t_94[k] = f_0 * kl_184[k];

        t_95[k] = f_0 * kl_185[k];

        t_96[k] = f_0 * kl_186[k];

        t_97[k] = f_0 * kl_187[k];

        t_98[k] = f_0 * kl_188[k];

        t_99[k] = f_0 * kl_189[k];

        t_100[k] = f_0 * kl_190[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, t_107, t_108, kl_191, \
                         kl_192, kl_193, kl_194, kl_195, kl_196, kl_197, \
                         kl_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * kl_191[k];

        t_102[k] = f_0 * kl_192[k];

        t_103[k] = f_0 * kl_193[k];

        t_104[k] = f_0 * kl_194[k];

        t_105[k] = f_0 * kl_195[k];

        t_106[k] = f_0 * kl_196[k];

        t_107[k] = f_0 * kl_197[k];

        t_108[k] = f_0 * kl_198[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, t_114, t_115, t_116, kl_199, \
                         kl_200, kl_201, kl_202, kl_203, kl_204, kl_205, \
                         kl_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_0 * kl_199[k];

        t_110[k] = f_0 * kl_200[k];

        t_111[k] = f_0 * kl_201[k];

        t_112[k] = f_0 * kl_202[k];

        t_113[k] = f_0 * kl_203[k];

        t_114[k] = f_0 * kl_204[k];

        t_115[k] = f_0 * kl_205[k];

        t_116[k] = f_0 * kl_206[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, t_122, t_123, t_124, kl_207, \
                         kl_208, kl_209, kl_210, kl_211, kl_212, kl_213, \
                         kl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_0 * kl_207[k];

        t_118[k] = f_0 * kl_208[k];

        t_119[k] = f_0 * kl_209[k];

        t_120[k] = f_0 * kl_210[k];

        t_121[k] = f_0 * kl_211[k];

        t_122[k] = f_0 * kl_212[k];

        t_123[k] = f_0 * kl_213[k];

        t_124[k] = f_0 * kl_214[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, t_131, t_132, kl_215, \
                         kl_216, kl_217, kl_218, kl_219, kl_220, kl_221, \
                         kl_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * kl_215[k];

        t_126[k] = f_0 * kl_216[k];

        t_127[k] = f_0 * kl_217[k];

        t_128[k] = f_0 * kl_218[k];

        t_129[k] = f_0 * kl_219[k];

        t_130[k] = f_0 * kl_220[k];

        t_131[k] = f_0 * kl_221[k];

        t_132[k] = f_0 * kl_222[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, t_138, hl_45, hl_46, hl_47, hl_48, \
                         kl_223, kl_224, kl_270, kl_271, kl_272, \
                         kl_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_0 * kl_223[k];

        t_134[k] = f_0 * kl_224[k];

        t_135[k] = -2.0 * hl_45[k]
                   + f_0 * kl_270[k];

        t_136[k] = -2.0 * hl_46[k]
                   + f_0 * kl_271[k];

        t_137[k] = -2.0 * hl_47[k]
                   + f_0 * kl_272[k];

        t_138[k] = -2.0 * hl_48[k]
                   + f_0 * kl_273[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, hl_49, hl_50, hl_51, hl_52, hl_53, \
                         kl_274, kl_275, kl_276, kl_277, kl_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -2.0 * hl_49[k]
                   + f_0 * kl_274[k];

        t_140[k] = -2.0 * hl_50[k]
                   + f_0 * kl_275[k];

        t_141[k] = -2.0 * hl_51[k]
                   + f_0 * kl_276[k];

        t_142[k] = -2.0 * hl_52[k]
                   + f_0 * kl_277[k];

        t_143[k] = -2.0 * hl_53[k]
                   + f_0 * kl_278[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, hl_54, hl_55, hl_56, hl_57, hl_58, \
                         kl_279, kl_280, kl_281, kl_282, kl_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -2.0 * hl_54[k]
                   + f_0 * kl_279[k];

        t_145[k] = -2.0 * hl_55[k]
                   + f_0 * kl_280[k];

        t_146[k] = -2.0 * hl_56[k]
                   + f_0 * kl_281[k];

        t_147[k] = -2.0 * hl_57[k]
                   + f_0 * kl_282[k];

        t_148[k] = -2.0 * hl_58[k]
                   + f_0 * kl_283[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, hl_59, hl_60, hl_61, hl_62, hl_63, \
                         kl_284, kl_285, kl_286, kl_287, kl_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -2.0 * hl_59[k]
                   + f_0 * kl_284[k];

        t_150[k] = -2.0 * hl_60[k]
                   + f_0 * kl_285[k];

        t_151[k] = -2.0 * hl_61[k]
                   + f_0 * kl_286[k];

        t_152[k] = -2.0 * hl_62[k]
                   + f_0 * kl_287[k];

        t_153[k] = -2.0 * hl_63[k]
                   + f_0 * kl_288[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, hl_64, hl_65, hl_66, hl_67, hl_68, \
                         kl_289, kl_290, kl_291, kl_292, kl_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -2.0 * hl_64[k]
                   + f_0 * kl_289[k];

        t_155[k] = -2.0 * hl_65[k]
                   + f_0 * kl_290[k];

        t_156[k] = -2.0 * hl_66[k]
                   + f_0 * kl_291[k];

        t_157[k] = -2.0 * hl_67[k]
                   + f_0 * kl_292[k];

        t_158[k] = -2.0 * hl_68[k]
                   + f_0 * kl_293[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, hl_69, hl_70, hl_71, hl_72, hl_73, \
                         kl_294, kl_295, kl_296, kl_297, kl_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -2.0 * hl_69[k]
                   + f_0 * kl_294[k];

        t_160[k] = -2.0 * hl_70[k]
                   + f_0 * kl_295[k];

        t_161[k] = -2.0 * hl_71[k]
                   + f_0 * kl_296[k];

        t_162[k] = -2.0 * hl_72[k]
                   + f_0 * kl_297[k];

        t_163[k] = -2.0 * hl_73[k]
                   + f_0 * kl_298[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, hl_74, hl_75, hl_76, hl_77, hl_78, \
                         kl_299, kl_300, kl_301, kl_302, kl_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -2.0 * hl_74[k]
                   + f_0 * kl_299[k];

        t_165[k] = -2.0 * hl_75[k]
                   + f_0 * kl_300[k];

        t_166[k] = -2.0 * hl_76[k]
                   + f_0 * kl_301[k];

        t_167[k] = -2.0 * hl_77[k]
                   + f_0 * kl_302[k];

        t_168[k] = -2.0 * hl_78[k]
                   + f_0 * kl_303[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, hl_79, hl_80, hl_81, hl_82, hl_83, \
                         kl_304, kl_305, kl_306, kl_307, kl_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -2.0 * hl_79[k]
                   + f_0 * kl_304[k];

        t_170[k] = -2.0 * hl_80[k]
                   + f_0 * kl_305[k];

        t_171[k] = -2.0 * hl_81[k]
                   + f_0 * kl_306[k];

        t_172[k] = -2.0 * hl_82[k]
                   + f_0 * kl_307[k];

        t_173[k] = -2.0 * hl_83[k]
                   + f_0 * kl_308[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, hl_84, hl_85, hl_86, hl_87, hl_88, \
                         kl_309, kl_310, kl_311, kl_312, kl_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -2.0 * hl_84[k]
                   + f_0 * kl_309[k];

        t_175[k] = -2.0 * hl_85[k]
                   + f_0 * kl_310[k];

        t_176[k] = -2.0 * hl_86[k]
                   + f_0 * kl_311[k];

        t_177[k] = -2.0 * hl_87[k]
                   + f_0 * kl_312[k];

        t_178[k] = -2.0 * hl_88[k]
                   + f_0 * kl_313[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
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

    const auto *hl_89 = buffer.data(hl + 89);
    const auto *hl_90 = buffer.data(hl + 90);
    const auto *hl_91 = buffer.data(hl + 91);
    const auto *hl_92 = buffer.data(hl + 92);
    const auto *hl_93 = buffer.data(hl + 93);
    const auto *hl_94 = buffer.data(hl + 94);
    const auto *hl_95 = buffer.data(hl + 95);
    const auto *hl_96 = buffer.data(hl + 96);
    const auto *hl_97 = buffer.data(hl + 97);
    const auto *hl_98 = buffer.data(hl + 98);
    const auto *hl_99 = buffer.data(hl + 99);
    const auto *hl_100 = buffer.data(hl + 100);
    const auto *hl_101 = buffer.data(hl + 101);
    const auto *hl_102 = buffer.data(hl + 102);
    const auto *hl_103 = buffer.data(hl + 103);
    const auto *hl_104 = buffer.data(hl + 104);
    const auto *hl_105 = buffer.data(hl + 105);
    const auto *hl_106 = buffer.data(hl + 106);
    const auto *hl_107 = buffer.data(hl + 107);
    const auto *hl_108 = buffer.data(hl + 108);
    const auto *hl_109 = buffer.data(hl + 109);
    const auto *hl_110 = buffer.data(hl + 110);
    const auto *hl_111 = buffer.data(hl + 111);
    const auto *hl_112 = buffer.data(hl + 112);
    const auto *hl_113 = buffer.data(hl + 113);
    const auto *hl_114 = buffer.data(hl + 114);
    const auto *hl_115 = buffer.data(hl + 115);
    const auto *hl_116 = buffer.data(hl + 116);
    const auto *hl_117 = buffer.data(hl + 117);
    const auto *hl_118 = buffer.data(hl + 118);
    const auto *hl_119 = buffer.data(hl + 119);
    const auto *hl_120 = buffer.data(hl + 120);
    const auto *hl_121 = buffer.data(hl + 121);
    const auto *hl_122 = buffer.data(hl + 122);
    const auto *hl_123 = buffer.data(hl + 123);
    const auto *hl_124 = buffer.data(hl + 124);
    const auto *hl_125 = buffer.data(hl + 125);
    const auto *hl_126 = buffer.data(hl + 126);
    const auto *hl_127 = buffer.data(hl + 127);
    const auto *hl_128 = buffer.data(hl + 128);
    const auto *hl_129 = buffer.data(hl + 129);
    const auto *hl_130 = buffer.data(hl + 130);
    const auto *hl_131 = buffer.data(hl + 131);
    const auto *hl_132 = buffer.data(hl + 132);
    const auto *hl_133 = buffer.data(hl + 133);
    const auto *hl_134 = buffer.data(hl + 134);
    const auto *hl_135 = buffer.data(hl + 135);
    const auto *hl_136 = buffer.data(hl + 136);
    const auto *hl_137 = buffer.data(hl + 137);
    const auto *hl_138 = buffer.data(hl + 138);
    const auto *hl_139 = buffer.data(hl + 139);
    const auto *hl_140 = buffer.data(hl + 140);
    const auto *hl_141 = buffer.data(hl + 141);
    const auto *hl_142 = buffer.data(hl + 142);
    const auto *hl_143 = buffer.data(hl + 143);
    const auto *hl_144 = buffer.data(hl + 144);
    const auto *hl_145 = buffer.data(hl + 145);
    const auto *hl_146 = buffer.data(hl + 146);
    const auto *hl_147 = buffer.data(hl + 147);
    const auto *hl_148 = buffer.data(hl + 148);
    const auto *hl_149 = buffer.data(hl + 149);
    const auto *hl_150 = buffer.data(hl + 150);
    const auto *hl_151 = buffer.data(hl + 151);
    const auto *hl_152 = buffer.data(hl + 152);
    const auto *hl_153 = buffer.data(hl + 153);
    const auto *hl_154 = buffer.data(hl + 154);
    const auto *hl_155 = buffer.data(hl + 155);
    const auto *hl_156 = buffer.data(hl + 156);
    const auto *hl_157 = buffer.data(hl + 157);
    const auto *hl_158 = buffer.data(hl + 158);
    const auto *hl_159 = buffer.data(hl + 159);
    const auto *hl_160 = buffer.data(hl + 160);
    const auto *hl_161 = buffer.data(hl + 161);
    const auto *hl_162 = buffer.data(hl + 162);
    const auto *hl_163 = buffer.data(hl + 163);
    const auto *hl_164 = buffer.data(hl + 164);
    const auto *hl_165 = buffer.data(hl + 165);
    const auto *hl_166 = buffer.data(hl + 166);
    const auto *hl_167 = buffer.data(hl + 167);
    const auto *hl_168 = buffer.data(hl + 168);
    const auto *hl_169 = buffer.data(hl + 169);
    const auto *hl_170 = buffer.data(hl + 170);
    const auto *hl_171 = buffer.data(hl + 171);
    const auto *hl_172 = buffer.data(hl + 172);
    const auto *hl_173 = buffer.data(hl + 173);
    const auto *hl_174 = buffer.data(hl + 174);
    const auto *hl_175 = buffer.data(hl + 175);
    const auto *hl_176 = buffer.data(hl + 176);
    const auto *hl_177 = buffer.data(hl + 177);
    const auto *hl_178 = buffer.data(hl + 178);
    const auto *hl_179 = buffer.data(hl + 179);
    const auto *hl_180 = buffer.data(hl + 180);
    const auto *hl_181 = buffer.data(hl + 181);
    const auto *hl_182 = buffer.data(hl + 182);
    const auto *hl_183 = buffer.data(hl + 183);
    const auto *hl_184 = buffer.data(hl + 184);
    const auto *hl_185 = buffer.data(hl + 185);
    const auto *hl_186 = buffer.data(hl + 186);
    const auto *hl_187 = buffer.data(hl + 187);
    const auto *hl_188 = buffer.data(hl + 188);
    const auto *hl_189 = buffer.data(hl + 189);
    const auto *hl_190 = buffer.data(hl + 190);
    const auto *hl_191 = buffer.data(hl + 191);
    const auto *hl_192 = buffer.data(hl + 192);
    const auto *hl_193 = buffer.data(hl + 193);
    const auto *hl_194 = buffer.data(hl + 194);
    const auto *hl_195 = buffer.data(hl + 195);
    const auto *hl_196 = buffer.data(hl + 196);
    const auto *hl_197 = buffer.data(hl + 197);
    const auto *hl_198 = buffer.data(hl + 198);
    const auto *hl_199 = buffer.data(hl + 199);
    const auto *hl_200 = buffer.data(hl + 200);
    const auto *hl_201 = buffer.data(hl + 201);
    const auto *hl_202 = buffer.data(hl + 202);
    const auto *hl_203 = buffer.data(hl + 203);
    const auto *hl_204 = buffer.data(hl + 204);
    const auto *hl_205 = buffer.data(hl + 205);
    const auto *hl_206 = buffer.data(hl + 206);
    const auto *hl_207 = buffer.data(hl + 207);
    const auto *hl_208 = buffer.data(hl + 208);
    const auto *hl_209 = buffer.data(hl + 209);

    const auto *kl_314 = buffer.data(kl + 314);
    const auto *kl_315 = buffer.data(kl + 315);
    const auto *kl_316 = buffer.data(kl + 316);
    const auto *kl_317 = buffer.data(kl + 317);
    const auto *kl_318 = buffer.data(kl + 318);
    const auto *kl_319 = buffer.data(kl + 319);
    const auto *kl_320 = buffer.data(kl + 320);
    const auto *kl_321 = buffer.data(kl + 321);
    const auto *kl_322 = buffer.data(kl + 322);
    const auto *kl_323 = buffer.data(kl + 323);
    const auto *kl_324 = buffer.data(kl + 324);
    const auto *kl_325 = buffer.data(kl + 325);
    const auto *kl_326 = buffer.data(kl + 326);
    const auto *kl_327 = buffer.data(kl + 327);
    const auto *kl_328 = buffer.data(kl + 328);
    const auto *kl_329 = buffer.data(kl + 329);
    const auto *kl_330 = buffer.data(kl + 330);
    const auto *kl_331 = buffer.data(kl + 331);
    const auto *kl_332 = buffer.data(kl + 332);
    const auto *kl_333 = buffer.data(kl + 333);
    const auto *kl_334 = buffer.data(kl + 334);
    const auto *kl_335 = buffer.data(kl + 335);
    const auto *kl_336 = buffer.data(kl + 336);
    const auto *kl_337 = buffer.data(kl + 337);
    const auto *kl_338 = buffer.data(kl + 338);
    const auto *kl_339 = buffer.data(kl + 339);
    const auto *kl_340 = buffer.data(kl + 340);
    const auto *kl_341 = buffer.data(kl + 341);
    const auto *kl_342 = buffer.data(kl + 342);
    const auto *kl_343 = buffer.data(kl + 343);
    const auto *kl_344 = buffer.data(kl + 344);
    const auto *kl_345 = buffer.data(kl + 345);
    const auto *kl_346 = buffer.data(kl + 346);
    const auto *kl_347 = buffer.data(kl + 347);
    const auto *kl_348 = buffer.data(kl + 348);
    const auto *kl_349 = buffer.data(kl + 349);
    const auto *kl_350 = buffer.data(kl + 350);
    const auto *kl_351 = buffer.data(kl + 351);
    const auto *kl_352 = buffer.data(kl + 352);
    const auto *kl_353 = buffer.data(kl + 353);
    const auto *kl_354 = buffer.data(kl + 354);
    const auto *kl_355 = buffer.data(kl + 355);
    const auto *kl_356 = buffer.data(kl + 356);
    const auto *kl_357 = buffer.data(kl + 357);
    const auto *kl_358 = buffer.data(kl + 358);
    const auto *kl_359 = buffer.data(kl + 359);
    const auto *kl_360 = buffer.data(kl + 360);
    const auto *kl_361 = buffer.data(kl + 361);
    const auto *kl_362 = buffer.data(kl + 362);
    const auto *kl_363 = buffer.data(kl + 363);
    const auto *kl_364 = buffer.data(kl + 364);
    const auto *kl_365 = buffer.data(kl + 365);
    const auto *kl_366 = buffer.data(kl + 366);
    const auto *kl_367 = buffer.data(kl + 367);
    const auto *kl_368 = buffer.data(kl + 368);
    const auto *kl_369 = buffer.data(kl + 369);
    const auto *kl_370 = buffer.data(kl + 370);
    const auto *kl_371 = buffer.data(kl + 371);
    const auto *kl_372 = buffer.data(kl + 372);
    const auto *kl_373 = buffer.data(kl + 373);
    const auto *kl_374 = buffer.data(kl + 374);
    const auto *kl_375 = buffer.data(kl + 375);
    const auto *kl_376 = buffer.data(kl + 376);
    const auto *kl_377 = buffer.data(kl + 377);
    const auto *kl_378 = buffer.data(kl + 378);
    const auto *kl_379 = buffer.data(kl + 379);
    const auto *kl_380 = buffer.data(kl + 380);
    const auto *kl_381 = buffer.data(kl + 381);
    const auto *kl_382 = buffer.data(kl + 382);
    const auto *kl_383 = buffer.data(kl + 383);
    const auto *kl_384 = buffer.data(kl + 384);
    const auto *kl_385 = buffer.data(kl + 385);
    const auto *kl_386 = buffer.data(kl + 386);
    const auto *kl_387 = buffer.data(kl + 387);
    const auto *kl_388 = buffer.data(kl + 388);
    const auto *kl_389 = buffer.data(kl + 389);
    const auto *kl_390 = buffer.data(kl + 390);
    const auto *kl_391 = buffer.data(kl + 391);
    const auto *kl_392 = buffer.data(kl + 392);
    const auto *kl_393 = buffer.data(kl + 393);
    const auto *kl_394 = buffer.data(kl + 394);
    const auto *kl_395 = buffer.data(kl + 395);
    const auto *kl_396 = buffer.data(kl + 396);
    const auto *kl_397 = buffer.data(kl + 397);
    const auto *kl_398 = buffer.data(kl + 398);
    const auto *kl_399 = buffer.data(kl + 399);
    const auto *kl_400 = buffer.data(kl + 400);
    const auto *kl_401 = buffer.data(kl + 401);
    const auto *kl_402 = buffer.data(kl + 402);
    const auto *kl_403 = buffer.data(kl + 403);
    const auto *kl_404 = buffer.data(kl + 404);
    const auto *kl_450 = buffer.data(kl + 450);
    const auto *kl_451 = buffer.data(kl + 451);
    const auto *kl_452 = buffer.data(kl + 452);
    const auto *kl_453 = buffer.data(kl + 453);
    const auto *kl_454 = buffer.data(kl + 454);
    const auto *kl_455 = buffer.data(kl + 455);
    const auto *kl_456 = buffer.data(kl + 456);
    const auto *kl_457 = buffer.data(kl + 457);
    const auto *kl_458 = buffer.data(kl + 458);
    const auto *kl_459 = buffer.data(kl + 459);
    const auto *kl_460 = buffer.data(kl + 460);
    const auto *kl_461 = buffer.data(kl + 461);
    const auto *kl_462 = buffer.data(kl + 462);
    const auto *kl_463 = buffer.data(kl + 463);
    const auto *kl_464 = buffer.data(kl + 464);
    const auto *kl_465 = buffer.data(kl + 465);
    const auto *kl_466 = buffer.data(kl + 466);
    const auto *kl_467 = buffer.data(kl + 467);
    const auto *kl_468 = buffer.data(kl + 468);
    const auto *kl_469 = buffer.data(kl + 469);
    const auto *kl_470 = buffer.data(kl + 470);
    const auto *kl_471 = buffer.data(kl + 471);
    const auto *kl_472 = buffer.data(kl + 472);
    const auto *kl_473 = buffer.data(kl + 473);
    const auto *kl_474 = buffer.data(kl + 474);
    const auto *kl_475 = buffer.data(kl + 475);
    const auto *kl_476 = buffer.data(kl + 476);
    const auto *kl_477 = buffer.data(kl + 477);
    const auto *kl_478 = buffer.data(kl + 478);
    const auto *kl_479 = buffer.data(kl + 479);
    const auto *kl_480 = buffer.data(kl + 480);
    const auto *kl_481 = buffer.data(kl + 481);
    const auto *kl_482 = buffer.data(kl + 482);
    const auto *kl_483 = buffer.data(kl + 483);
    const auto *kl_484 = buffer.data(kl + 484);
    const auto *kl_485 = buffer.data(kl + 485);
    const auto *kl_486 = buffer.data(kl + 486);
    const auto *kl_487 = buffer.data(kl + 487);
    const auto *kl_488 = buffer.data(kl + 488);
    const auto *kl_489 = buffer.data(kl + 489);
    const auto *kl_490 = buffer.data(kl + 490);
    const auto *kl_491 = buffer.data(kl + 491);
    const auto *kl_492 = buffer.data(kl + 492);
    const auto *kl_493 = buffer.data(kl + 493);
    const auto *kl_494 = buffer.data(kl + 494);
    const auto *kl_495 = buffer.data(kl + 495);
    const auto *kl_496 = buffer.data(kl + 496);
    const auto *kl_497 = buffer.data(kl + 497);
    const auto *kl_498 = buffer.data(kl + 498);
    const auto *kl_499 = buffer.data(kl + 499);
    const auto *kl_500 = buffer.data(kl + 500);
    const auto *kl_501 = buffer.data(kl + 501);
    const auto *kl_502 = buffer.data(kl + 502);
    const auto *kl_503 = buffer.data(kl + 503);
    const auto *kl_504 = buffer.data(kl + 504);
    const auto *kl_505 = buffer.data(kl + 505);
    const auto *kl_506 = buffer.data(kl + 506);
    const auto *kl_507 = buffer.data(kl + 507);
    const auto *kl_508 = buffer.data(kl + 508);
    const auto *kl_509 = buffer.data(kl + 509);
    const auto *kl_510 = buffer.data(kl + 510);
    const auto *kl_511 = buffer.data(kl + 511);
    const auto *kl_512 = buffer.data(kl + 512);
    const auto *kl_513 = buffer.data(kl + 513);
    const auto *kl_514 = buffer.data(kl + 514);
    const auto *kl_515 = buffer.data(kl + 515);
    const auto *kl_516 = buffer.data(kl + 516);
    const auto *kl_517 = buffer.data(kl + 517);
    const auto *kl_518 = buffer.data(kl + 518);
    const auto *kl_519 = buffer.data(kl + 519);
    const auto *kl_520 = buffer.data(kl + 520);
    const auto *kl_521 = buffer.data(kl + 521);
    const auto *kl_522 = buffer.data(kl + 522);
    const auto *kl_523 = buffer.data(kl + 523);
    const auto *kl_524 = buffer.data(kl + 524);

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, hl_89, hl_90, hl_91, hl_92, hl_93, \
                         kl_314, kl_315, kl_316, kl_317, kl_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -2.0 * hl_89[k]
                   + f_0 * kl_314[k];

        t_180[k] = -hl_90[k]
                   + f_0 * kl_315[k];

        t_181[k] = -hl_91[k]
                   + f_0 * kl_316[k];

        t_182[k] = -hl_92[k]
                   + f_0 * kl_317[k];

        t_183[k] = -hl_93[k]
                   + f_0 * kl_318[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, hl_94, hl_95, hl_96, hl_97, hl_98, \
                         kl_319, kl_320, kl_321, kl_322, kl_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -hl_94[k]
                   + f_0 * kl_319[k];

        t_185[k] = -hl_95[k]
                   + f_0 * kl_320[k];

        t_186[k] = -hl_96[k]
                   + f_0 * kl_321[k];

        t_187[k] = -hl_97[k]
                   + f_0 * kl_322[k];

        t_188[k] = -hl_98[k]
                   + f_0 * kl_323[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, hl_99, hl_100, hl_101, hl_102, \
                         hl_103, kl_324, kl_325, kl_326, kl_327, \
                         kl_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = -hl_99[k]
                   + f_0 * kl_324[k];

        t_190[k] = -hl_100[k]
                   + f_0 * kl_325[k];

        t_191[k] = -hl_101[k]
                   + f_0 * kl_326[k];

        t_192[k] = -hl_102[k]
                   + f_0 * kl_327[k];

        t_193[k] = -hl_103[k]
                   + f_0 * kl_328[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, hl_104, hl_105, hl_106, hl_107, \
                         hl_108, kl_329, kl_330, kl_331, kl_332, \
                         kl_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = -hl_104[k]
                   + f_0 * kl_329[k];

        t_195[k] = -hl_105[k]
                   + f_0 * kl_330[k];

        t_196[k] = -hl_106[k]
                   + f_0 * kl_331[k];

        t_197[k] = -hl_107[k]
                   + f_0 * kl_332[k];

        t_198[k] = -hl_108[k]
                   + f_0 * kl_333[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, hl_109, hl_110, hl_111, hl_112, \
                         hl_113, kl_334, kl_335, kl_336, kl_337, \
                         kl_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = -hl_109[k]
                   + f_0 * kl_334[k];

        t_200[k] = -hl_110[k]
                   + f_0 * kl_335[k];

        t_201[k] = -hl_111[k]
                   + f_0 * kl_336[k];

        t_202[k] = -hl_112[k]
                   + f_0 * kl_337[k];

        t_203[k] = -hl_113[k]
                   + f_0 * kl_338[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, hl_114, hl_115, hl_116, hl_117, \
                         hl_118, kl_339, kl_340, kl_341, kl_342, \
                         kl_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = -hl_114[k]
                   + f_0 * kl_339[k];

        t_205[k] = -hl_115[k]
                   + f_0 * kl_340[k];

        t_206[k] = -hl_116[k]
                   + f_0 * kl_341[k];

        t_207[k] = -hl_117[k]
                   + f_0 * kl_342[k];

        t_208[k] = -hl_118[k]
                   + f_0 * kl_343[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, hl_119, hl_120, hl_121, hl_122, \
                         hl_123, kl_344, kl_345, kl_346, kl_347, \
                         kl_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = -hl_119[k]
                   + f_0 * kl_344[k];

        t_210[k] = -hl_120[k]
                   + f_0 * kl_345[k];

        t_211[k] = -hl_121[k]
                   + f_0 * kl_346[k];

        t_212[k] = -hl_122[k]
                   + f_0 * kl_347[k];

        t_213[k] = -hl_123[k]
                   + f_0 * kl_348[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, hl_124, hl_125, hl_126, hl_127, \
                         hl_128, kl_349, kl_350, kl_351, kl_352, \
                         kl_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = -hl_124[k]
                   + f_0 * kl_349[k];

        t_215[k] = -hl_125[k]
                   + f_0 * kl_350[k];

        t_216[k] = -hl_126[k]
                   + f_0 * kl_351[k];

        t_217[k] = -hl_127[k]
                   + f_0 * kl_352[k];

        t_218[k] = -hl_128[k]
                   + f_0 * kl_353[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, hl_129, hl_130, hl_131, hl_132, \
                         hl_133, kl_354, kl_355, kl_356, kl_357, \
                         kl_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = -hl_129[k]
                   + f_0 * kl_354[k];

        t_220[k] = -hl_130[k]
                   + f_0 * kl_355[k];

        t_221[k] = -hl_131[k]
                   + f_0 * kl_356[k];

        t_222[k] = -hl_132[k]
                   + f_0 * kl_357[k];

        t_223[k] = -hl_133[k]
                   + f_0 * kl_358[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, t_229, t_230, hl_134, kl_359, \
                         kl_360, kl_361, kl_362, kl_363, kl_364, \
                         kl_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = -hl_134[k]
                   + f_0 * kl_359[k];

        t_225[k] = f_0 * kl_360[k];

        t_226[k] = f_0 * kl_361[k];

        t_227[k] = f_0 * kl_362[k];

        t_228[k] = f_0 * kl_363[k];

        t_229[k] = f_0 * kl_364[k];

        t_230[k] = f_0 * kl_365[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, t_236, t_237, t_238, kl_366, \
                         kl_367, kl_368, kl_369, kl_370, kl_371, kl_372, \
                         kl_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_0 * kl_366[k];

        t_232[k] = f_0 * kl_367[k];

        t_233[k] = f_0 * kl_368[k];

        t_234[k] = f_0 * kl_369[k];

        t_235[k] = f_0 * kl_370[k];

        t_236[k] = f_0 * kl_371[k];

        t_237[k] = f_0 * kl_372[k];

        t_238[k] = f_0 * kl_373[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, t_244, t_245, t_246, kl_374, \
                         kl_375, kl_376, kl_377, kl_378, kl_379, kl_380, \
                         kl_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_0 * kl_374[k];

        t_240[k] = f_0 * kl_375[k];

        t_241[k] = f_0 * kl_376[k];

        t_242[k] = f_0 * kl_377[k];

        t_243[k] = f_0 * kl_378[k];

        t_244[k] = f_0 * kl_379[k];

        t_245[k] = f_0 * kl_380[k];

        t_246[k] = f_0 * kl_381[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, t_252, t_253, t_254, kl_382, \
                         kl_383, kl_384, kl_385, kl_386, kl_387, kl_388, \
                         kl_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_0 * kl_382[k];

        t_248[k] = f_0 * kl_383[k];

        t_249[k] = f_0 * kl_384[k];

        t_250[k] = f_0 * kl_385[k];

        t_251[k] = f_0 * kl_386[k];

        t_252[k] = f_0 * kl_387[k];

        t_253[k] = f_0 * kl_388[k];

        t_254[k] = f_0 * kl_389[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, t_260, t_261, t_262, kl_390, \
                         kl_391, kl_392, kl_393, kl_394, kl_395, kl_396, \
                         kl_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_0 * kl_390[k];

        t_256[k] = f_0 * kl_391[k];

        t_257[k] = f_0 * kl_392[k];

        t_258[k] = f_0 * kl_393[k];

        t_259[k] = f_0 * kl_394[k];

        t_260[k] = f_0 * kl_395[k];

        t_261[k] = f_0 * kl_396[k];

        t_262[k] = f_0 * kl_397[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, t_268, t_269, kl_398, kl_399, \
                         kl_400, kl_401, kl_402, kl_403, kl_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_0 * kl_398[k];

        t_264[k] = f_0 * kl_399[k];

        t_265[k] = f_0 * kl_400[k];

        t_266[k] = f_0 * kl_401[k];

        t_267[k] = f_0 * kl_402[k];

        t_268[k] = f_0 * kl_403[k];

        t_269[k] = f_0 * kl_404[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, hl_135, hl_136, hl_137, hl_138, \
                         hl_139, kl_450, kl_451, kl_452, kl_453, \
                         kl_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -3.0 * hl_135[k]
                   + f_0 * kl_450[k];

        t_271[k] = -3.0 * hl_136[k]
                   + f_0 * kl_451[k];

        t_272[k] = -3.0 * hl_137[k]
                   + f_0 * kl_452[k];

        t_273[k] = -3.0 * hl_138[k]
                   + f_0 * kl_453[k];

        t_274[k] = -3.0 * hl_139[k]
                   + f_0 * kl_454[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, hl_140, hl_141, hl_142, hl_143, \
                         hl_144, kl_455, kl_456, kl_457, kl_458, \
                         kl_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * hl_140[k]
                   + f_0 * kl_455[k];

        t_276[k] = -3.0 * hl_141[k]
                   + f_0 * kl_456[k];

        t_277[k] = -3.0 * hl_142[k]
                   + f_0 * kl_457[k];

        t_278[k] = -3.0 * hl_143[k]
                   + f_0 * kl_458[k];

        t_279[k] = -3.0 * hl_144[k]
                   + f_0 * kl_459[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, hl_145, hl_146, hl_147, hl_148, \
                         hl_149, kl_460, kl_461, kl_462, kl_463, \
                         kl_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * hl_145[k]
                   + f_0 * kl_460[k];

        t_281[k] = -3.0 * hl_146[k]
                   + f_0 * kl_461[k];

        t_282[k] = -3.0 * hl_147[k]
                   + f_0 * kl_462[k];

        t_283[k] = -3.0 * hl_148[k]
                   + f_0 * kl_463[k];

        t_284[k] = -3.0 * hl_149[k]
                   + f_0 * kl_464[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, hl_150, hl_151, hl_152, hl_153, \
                         hl_154, kl_465, kl_466, kl_467, kl_468, \
                         kl_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * hl_150[k]
                   + f_0 * kl_465[k];

        t_286[k] = -3.0 * hl_151[k]
                   + f_0 * kl_466[k];

        t_287[k] = -3.0 * hl_152[k]
                   + f_0 * kl_467[k];

        t_288[k] = -3.0 * hl_153[k]
                   + f_0 * kl_468[k];

        t_289[k] = -3.0 * hl_154[k]
                   + f_0 * kl_469[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, hl_155, hl_156, hl_157, hl_158, \
                         hl_159, kl_470, kl_471, kl_472, kl_473, \
                         kl_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * hl_155[k]
                   + f_0 * kl_470[k];

        t_291[k] = -3.0 * hl_156[k]
                   + f_0 * kl_471[k];

        t_292[k] = -3.0 * hl_157[k]
                   + f_0 * kl_472[k];

        t_293[k] = -3.0 * hl_158[k]
                   + f_0 * kl_473[k];

        t_294[k] = -3.0 * hl_159[k]
                   + f_0 * kl_474[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, hl_160, hl_161, hl_162, hl_163, \
                         hl_164, kl_475, kl_476, kl_477, kl_478, \
                         kl_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -3.0 * hl_160[k]
                   + f_0 * kl_475[k];

        t_296[k] = -3.0 * hl_161[k]
                   + f_0 * kl_476[k];

        t_297[k] = -3.0 * hl_162[k]
                   + f_0 * kl_477[k];

        t_298[k] = -3.0 * hl_163[k]
                   + f_0 * kl_478[k];

        t_299[k] = -3.0 * hl_164[k]
                   + f_0 * kl_479[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, hl_165, hl_166, hl_167, hl_168, \
                         hl_169, kl_480, kl_481, kl_482, kl_483, \
                         kl_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -3.0 * hl_165[k]
                   + f_0 * kl_480[k];

        t_301[k] = -3.0 * hl_166[k]
                   + f_0 * kl_481[k];

        t_302[k] = -3.0 * hl_167[k]
                   + f_0 * kl_482[k];

        t_303[k] = -3.0 * hl_168[k]
                   + f_0 * kl_483[k];

        t_304[k] = -3.0 * hl_169[k]
                   + f_0 * kl_484[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, hl_170, hl_171, hl_172, hl_173, \
                         hl_174, kl_485, kl_486, kl_487, kl_488, \
                         kl_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -3.0 * hl_170[k]
                   + f_0 * kl_485[k];

        t_306[k] = -3.0 * hl_171[k]
                   + f_0 * kl_486[k];

        t_307[k] = -3.0 * hl_172[k]
                   + f_0 * kl_487[k];

        t_308[k] = -3.0 * hl_173[k]
                   + f_0 * kl_488[k];

        t_309[k] = -3.0 * hl_174[k]
                   + f_0 * kl_489[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, hl_175, hl_176, hl_177, hl_178, \
                         hl_179, kl_490, kl_491, kl_492, kl_493, \
                         kl_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * hl_175[k]
                   + f_0 * kl_490[k];

        t_311[k] = -3.0 * hl_176[k]
                   + f_0 * kl_491[k];

        t_312[k] = -3.0 * hl_177[k]
                   + f_0 * kl_492[k];

        t_313[k] = -3.0 * hl_178[k]
                   + f_0 * kl_493[k];

        t_314[k] = -3.0 * hl_179[k]
                   + f_0 * kl_494[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, hl_180, hl_181, hl_182, hl_183, \
                         hl_184, kl_495, kl_496, kl_497, kl_498, \
                         kl_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -2.0 * hl_180[k]
                   + f_0 * kl_495[k];

        t_316[k] = -2.0 * hl_181[k]
                   + f_0 * kl_496[k];

        t_317[k] = -2.0 * hl_182[k]
                   + f_0 * kl_497[k];

        t_318[k] = -2.0 * hl_183[k]
                   + f_0 * kl_498[k];

        t_319[k] = -2.0 * hl_184[k]
                   + f_0 * kl_499[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, hl_185, hl_186, hl_187, hl_188, \
                         hl_189, kl_500, kl_501, kl_502, kl_503, \
                         kl_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -2.0 * hl_185[k]
                   + f_0 * kl_500[k];

        t_321[k] = -2.0 * hl_186[k]
                   + f_0 * kl_501[k];

        t_322[k] = -2.0 * hl_187[k]
                   + f_0 * kl_502[k];

        t_323[k] = -2.0 * hl_188[k]
                   + f_0 * kl_503[k];

        t_324[k] = -2.0 * hl_189[k]
                   + f_0 * kl_504[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, hl_190, hl_191, hl_192, hl_193, \
                         hl_194, kl_505, kl_506, kl_507, kl_508, \
                         kl_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -2.0 * hl_190[k]
                   + f_0 * kl_505[k];

        t_326[k] = -2.0 * hl_191[k]
                   + f_0 * kl_506[k];

        t_327[k] = -2.0 * hl_192[k]
                   + f_0 * kl_507[k];

        t_328[k] = -2.0 * hl_193[k]
                   + f_0 * kl_508[k];

        t_329[k] = -2.0 * hl_194[k]
                   + f_0 * kl_509[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, hl_195, hl_196, hl_197, hl_198, \
                         hl_199, kl_510, kl_511, kl_512, kl_513, \
                         kl_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -2.0 * hl_195[k]
                   + f_0 * kl_510[k];

        t_331[k] = -2.0 * hl_196[k]
                   + f_0 * kl_511[k];

        t_332[k] = -2.0 * hl_197[k]
                   + f_0 * kl_512[k];

        t_333[k] = -2.0 * hl_198[k]
                   + f_0 * kl_513[k];

        t_334[k] = -2.0 * hl_199[k]
                   + f_0 * kl_514[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, hl_200, hl_201, hl_202, hl_203, \
                         hl_204, kl_515, kl_516, kl_517, kl_518, \
                         kl_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -2.0 * hl_200[k]
                   + f_0 * kl_515[k];

        t_336[k] = -2.0 * hl_201[k]
                   + f_0 * kl_516[k];

        t_337[k] = -2.0 * hl_202[k]
                   + f_0 * kl_517[k];

        t_338[k] = -2.0 * hl_203[k]
                   + f_0 * kl_518[k];

        t_339[k] = -2.0 * hl_204[k]
                   + f_0 * kl_519[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, hl_205, hl_206, hl_207, hl_208, \
                         hl_209, kl_520, kl_521, kl_522, kl_523, \
                         kl_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * hl_205[k]
                   + f_0 * kl_520[k];

        t_341[k] = -2.0 * hl_206[k]
                   + f_0 * kl_521[k];

        t_342[k] = -2.0 * hl_207[k]
                   + f_0 * kl_522[k];

        t_343[k] = -2.0 * hl_208[k]
                   + f_0 * kl_523[k];

        t_344[k] = -2.0 * hl_209[k]
                   + f_0 * kl_524[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hl_210 = buffer.data(hl + 210);
    const auto *hl_211 = buffer.data(hl + 211);
    const auto *hl_212 = buffer.data(hl + 212);
    const auto *hl_213 = buffer.data(hl + 213);
    const auto *hl_214 = buffer.data(hl + 214);
    const auto *hl_215 = buffer.data(hl + 215);
    const auto *hl_216 = buffer.data(hl + 216);
    const auto *hl_217 = buffer.data(hl + 217);
    const auto *hl_218 = buffer.data(hl + 218);
    const auto *hl_219 = buffer.data(hl + 219);
    const auto *hl_220 = buffer.data(hl + 220);
    const auto *hl_221 = buffer.data(hl + 221);
    const auto *hl_222 = buffer.data(hl + 222);
    const auto *hl_223 = buffer.data(hl + 223);
    const auto *hl_224 = buffer.data(hl + 224);
    const auto *hl_225 = buffer.data(hl + 225);
    const auto *hl_226 = buffer.data(hl + 226);
    const auto *hl_227 = buffer.data(hl + 227);
    const auto *hl_228 = buffer.data(hl + 228);
    const auto *hl_229 = buffer.data(hl + 229);
    const auto *hl_230 = buffer.data(hl + 230);
    const auto *hl_231 = buffer.data(hl + 231);
    const auto *hl_232 = buffer.data(hl + 232);
    const auto *hl_233 = buffer.data(hl + 233);
    const auto *hl_234 = buffer.data(hl + 234);
    const auto *hl_235 = buffer.data(hl + 235);
    const auto *hl_236 = buffer.data(hl + 236);
    const auto *hl_237 = buffer.data(hl + 237);
    const auto *hl_238 = buffer.data(hl + 238);
    const auto *hl_239 = buffer.data(hl + 239);
    const auto *hl_240 = buffer.data(hl + 240);
    const auto *hl_241 = buffer.data(hl + 241);
    const auto *hl_242 = buffer.data(hl + 242);
    const auto *hl_243 = buffer.data(hl + 243);
    const auto *hl_244 = buffer.data(hl + 244);
    const auto *hl_245 = buffer.data(hl + 245);
    const auto *hl_246 = buffer.data(hl + 246);
    const auto *hl_247 = buffer.data(hl + 247);
    const auto *hl_248 = buffer.data(hl + 248);
    const auto *hl_249 = buffer.data(hl + 249);
    const auto *hl_250 = buffer.data(hl + 250);
    const auto *hl_251 = buffer.data(hl + 251);
    const auto *hl_252 = buffer.data(hl + 252);
    const auto *hl_253 = buffer.data(hl + 253);
    const auto *hl_254 = buffer.data(hl + 254);
    const auto *hl_255 = buffer.data(hl + 255);
    const auto *hl_256 = buffer.data(hl + 256);
    const auto *hl_257 = buffer.data(hl + 257);
    const auto *hl_258 = buffer.data(hl + 258);
    const auto *hl_259 = buffer.data(hl + 259);
    const auto *hl_260 = buffer.data(hl + 260);
    const auto *hl_261 = buffer.data(hl + 261);
    const auto *hl_262 = buffer.data(hl + 262);
    const auto *hl_263 = buffer.data(hl + 263);
    const auto *hl_264 = buffer.data(hl + 264);
    const auto *hl_265 = buffer.data(hl + 265);
    const auto *hl_266 = buffer.data(hl + 266);
    const auto *hl_267 = buffer.data(hl + 267);
    const auto *hl_268 = buffer.data(hl + 268);
    const auto *hl_269 = buffer.data(hl + 269);
    const auto *hl_270 = buffer.data(hl + 270);
    const auto *hl_271 = buffer.data(hl + 271);
    const auto *hl_272 = buffer.data(hl + 272);
    const auto *hl_273 = buffer.data(hl + 273);
    const auto *hl_274 = buffer.data(hl + 274);
    const auto *hl_275 = buffer.data(hl + 275);
    const auto *hl_276 = buffer.data(hl + 276);
    const auto *hl_277 = buffer.data(hl + 277);
    const auto *hl_278 = buffer.data(hl + 278);
    const auto *hl_279 = buffer.data(hl + 279);
    const auto *hl_280 = buffer.data(hl + 280);
    const auto *hl_281 = buffer.data(hl + 281);
    const auto *hl_282 = buffer.data(hl + 282);
    const auto *hl_283 = buffer.data(hl + 283);
    const auto *hl_284 = buffer.data(hl + 284);
    const auto *hl_285 = buffer.data(hl + 285);
    const auto *hl_286 = buffer.data(hl + 286);
    const auto *hl_287 = buffer.data(hl + 287);
    const auto *hl_288 = buffer.data(hl + 288);
    const auto *hl_289 = buffer.data(hl + 289);
    const auto *hl_290 = buffer.data(hl + 290);
    const auto *hl_291 = buffer.data(hl + 291);
    const auto *hl_292 = buffer.data(hl + 292);
    const auto *hl_293 = buffer.data(hl + 293);
    const auto *hl_294 = buffer.data(hl + 294);
    const auto *hl_295 = buffer.data(hl + 295);
    const auto *hl_296 = buffer.data(hl + 296);
    const auto *hl_297 = buffer.data(hl + 297);
    const auto *hl_298 = buffer.data(hl + 298);
    const auto *hl_299 = buffer.data(hl + 299);
    const auto *hl_300 = buffer.data(hl + 300);
    const auto *hl_301 = buffer.data(hl + 301);
    const auto *hl_302 = buffer.data(hl + 302);
    const auto *hl_303 = buffer.data(hl + 303);
    const auto *hl_304 = buffer.data(hl + 304);
    const auto *hl_305 = buffer.data(hl + 305);
    const auto *hl_306 = buffer.data(hl + 306);
    const auto *hl_307 = buffer.data(hl + 307);
    const auto *hl_308 = buffer.data(hl + 308);
    const auto *hl_309 = buffer.data(hl + 309);
    const auto *hl_310 = buffer.data(hl + 310);
    const auto *hl_311 = buffer.data(hl + 311);
    const auto *hl_312 = buffer.data(hl + 312);
    const auto *hl_313 = buffer.data(hl + 313);
    const auto *hl_314 = buffer.data(hl + 314);
    const auto *hl_315 = buffer.data(hl + 315);
    const auto *hl_316 = buffer.data(hl + 316);
    const auto *hl_317 = buffer.data(hl + 317);
    const auto *hl_318 = buffer.data(hl + 318);
    const auto *hl_319 = buffer.data(hl + 319);
    const auto *hl_320 = buffer.data(hl + 320);
    const auto *hl_321 = buffer.data(hl + 321);
    const auto *hl_322 = buffer.data(hl + 322);
    const auto *hl_323 = buffer.data(hl + 323);
    const auto *hl_324 = buffer.data(hl + 324);
    const auto *hl_325 = buffer.data(hl + 325);
    const auto *hl_326 = buffer.data(hl + 326);

    const auto *kl_525 = buffer.data(kl + 525);
    const auto *kl_526 = buffer.data(kl + 526);
    const auto *kl_527 = buffer.data(kl + 527);
    const auto *kl_528 = buffer.data(kl + 528);
    const auto *kl_529 = buffer.data(kl + 529);
    const auto *kl_530 = buffer.data(kl + 530);
    const auto *kl_531 = buffer.data(kl + 531);
    const auto *kl_532 = buffer.data(kl + 532);
    const auto *kl_533 = buffer.data(kl + 533);
    const auto *kl_534 = buffer.data(kl + 534);
    const auto *kl_535 = buffer.data(kl + 535);
    const auto *kl_536 = buffer.data(kl + 536);
    const auto *kl_537 = buffer.data(kl + 537);
    const auto *kl_538 = buffer.data(kl + 538);
    const auto *kl_539 = buffer.data(kl + 539);
    const auto *kl_540 = buffer.data(kl + 540);
    const auto *kl_541 = buffer.data(kl + 541);
    const auto *kl_542 = buffer.data(kl + 542);
    const auto *kl_543 = buffer.data(kl + 543);
    const auto *kl_544 = buffer.data(kl + 544);
    const auto *kl_545 = buffer.data(kl + 545);
    const auto *kl_546 = buffer.data(kl + 546);
    const auto *kl_547 = buffer.data(kl + 547);
    const auto *kl_548 = buffer.data(kl + 548);
    const auto *kl_549 = buffer.data(kl + 549);
    const auto *kl_550 = buffer.data(kl + 550);
    const auto *kl_551 = buffer.data(kl + 551);
    const auto *kl_552 = buffer.data(kl + 552);
    const auto *kl_553 = buffer.data(kl + 553);
    const auto *kl_554 = buffer.data(kl + 554);
    const auto *kl_555 = buffer.data(kl + 555);
    const auto *kl_556 = buffer.data(kl + 556);
    const auto *kl_557 = buffer.data(kl + 557);
    const auto *kl_558 = buffer.data(kl + 558);
    const auto *kl_559 = buffer.data(kl + 559);
    const auto *kl_560 = buffer.data(kl + 560);
    const auto *kl_561 = buffer.data(kl + 561);
    const auto *kl_562 = buffer.data(kl + 562);
    const auto *kl_563 = buffer.data(kl + 563);
    const auto *kl_564 = buffer.data(kl + 564);
    const auto *kl_565 = buffer.data(kl + 565);
    const auto *kl_566 = buffer.data(kl + 566);
    const auto *kl_567 = buffer.data(kl + 567);
    const auto *kl_568 = buffer.data(kl + 568);
    const auto *kl_569 = buffer.data(kl + 569);
    const auto *kl_570 = buffer.data(kl + 570);
    const auto *kl_571 = buffer.data(kl + 571);
    const auto *kl_572 = buffer.data(kl + 572);
    const auto *kl_573 = buffer.data(kl + 573);
    const auto *kl_574 = buffer.data(kl + 574);
    const auto *kl_575 = buffer.data(kl + 575);
    const auto *kl_576 = buffer.data(kl + 576);
    const auto *kl_577 = buffer.data(kl + 577);
    const auto *kl_578 = buffer.data(kl + 578);
    const auto *kl_579 = buffer.data(kl + 579);
    const auto *kl_580 = buffer.data(kl + 580);
    const auto *kl_581 = buffer.data(kl + 581);
    const auto *kl_582 = buffer.data(kl + 582);
    const auto *kl_583 = buffer.data(kl + 583);
    const auto *kl_584 = buffer.data(kl + 584);
    const auto *kl_585 = buffer.data(kl + 585);
    const auto *kl_586 = buffer.data(kl + 586);
    const auto *kl_587 = buffer.data(kl + 587);
    const auto *kl_588 = buffer.data(kl + 588);
    const auto *kl_589 = buffer.data(kl + 589);
    const auto *kl_590 = buffer.data(kl + 590);
    const auto *kl_591 = buffer.data(kl + 591);
    const auto *kl_592 = buffer.data(kl + 592);
    const auto *kl_593 = buffer.data(kl + 593);
    const auto *kl_594 = buffer.data(kl + 594);
    const auto *kl_595 = buffer.data(kl + 595);
    const auto *kl_596 = buffer.data(kl + 596);
    const auto *kl_597 = buffer.data(kl + 597);
    const auto *kl_598 = buffer.data(kl + 598);
    const auto *kl_599 = buffer.data(kl + 599);
    const auto *kl_600 = buffer.data(kl + 600);
    const auto *kl_601 = buffer.data(kl + 601);
    const auto *kl_602 = buffer.data(kl + 602);
    const auto *kl_603 = buffer.data(kl + 603);
    const auto *kl_604 = buffer.data(kl + 604);
    const auto *kl_605 = buffer.data(kl + 605);
    const auto *kl_606 = buffer.data(kl + 606);
    const auto *kl_607 = buffer.data(kl + 607);
    const auto *kl_608 = buffer.data(kl + 608);
    const auto *kl_609 = buffer.data(kl + 609);
    const auto *kl_610 = buffer.data(kl + 610);
    const auto *kl_611 = buffer.data(kl + 611);
    const auto *kl_612 = buffer.data(kl + 612);
    const auto *kl_613 = buffer.data(kl + 613);
    const auto *kl_614 = buffer.data(kl + 614);
    const auto *kl_615 = buffer.data(kl + 615);
    const auto *kl_616 = buffer.data(kl + 616);
    const auto *kl_617 = buffer.data(kl + 617);
    const auto *kl_618 = buffer.data(kl + 618);
    const auto *kl_619 = buffer.data(kl + 619);
    const auto *kl_620 = buffer.data(kl + 620);
    const auto *kl_621 = buffer.data(kl + 621);
    const auto *kl_622 = buffer.data(kl + 622);
    const auto *kl_623 = buffer.data(kl + 623);
    const auto *kl_624 = buffer.data(kl + 624);
    const auto *kl_625 = buffer.data(kl + 625);
    const auto *kl_626 = buffer.data(kl + 626);
    const auto *kl_627 = buffer.data(kl + 627);
    const auto *kl_628 = buffer.data(kl + 628);
    const auto *kl_629 = buffer.data(kl + 629);
    const auto *kl_675 = buffer.data(kl + 675);
    const auto *kl_676 = buffer.data(kl + 676);
    const auto *kl_677 = buffer.data(kl + 677);
    const auto *kl_678 = buffer.data(kl + 678);
    const auto *kl_679 = buffer.data(kl + 679);
    const auto *kl_680 = buffer.data(kl + 680);
    const auto *kl_681 = buffer.data(kl + 681);
    const auto *kl_682 = buffer.data(kl + 682);
    const auto *kl_683 = buffer.data(kl + 683);
    const auto *kl_684 = buffer.data(kl + 684);
    const auto *kl_685 = buffer.data(kl + 685);
    const auto *kl_686 = buffer.data(kl + 686);
    const auto *kl_687 = buffer.data(kl + 687);
    const auto *kl_688 = buffer.data(kl + 688);
    const auto *kl_689 = buffer.data(kl + 689);
    const auto *kl_690 = buffer.data(kl + 690);
    const auto *kl_691 = buffer.data(kl + 691);
    const auto *kl_692 = buffer.data(kl + 692);
    const auto *kl_693 = buffer.data(kl + 693);
    const auto *kl_694 = buffer.data(kl + 694);
    const auto *kl_695 = buffer.data(kl + 695);
    const auto *kl_696 = buffer.data(kl + 696);
    const auto *kl_697 = buffer.data(kl + 697);
    const auto *kl_698 = buffer.data(kl + 698);
    const auto *kl_699 = buffer.data(kl + 699);
    const auto *kl_700 = buffer.data(kl + 700);
    const auto *kl_701 = buffer.data(kl + 701);
    const auto *kl_702 = buffer.data(kl + 702);
    const auto *kl_703 = buffer.data(kl + 703);
    const auto *kl_704 = buffer.data(kl + 704);
    const auto *kl_705 = buffer.data(kl + 705);
    const auto *kl_706 = buffer.data(kl + 706);
    const auto *kl_707 = buffer.data(kl + 707);
    const auto *kl_708 = buffer.data(kl + 708);
    const auto *kl_709 = buffer.data(kl + 709);
    const auto *kl_710 = buffer.data(kl + 710);
    const auto *kl_711 = buffer.data(kl + 711);
    const auto *kl_712 = buffer.data(kl + 712);
    const auto *kl_713 = buffer.data(kl + 713);
    const auto *kl_714 = buffer.data(kl + 714);
    const auto *kl_715 = buffer.data(kl + 715);
    const auto *kl_716 = buffer.data(kl + 716);
    const auto *kl_717 = buffer.data(kl + 717);
    const auto *kl_718 = buffer.data(kl + 718);
    const auto *kl_719 = buffer.data(kl + 719);
    const auto *kl_720 = buffer.data(kl + 720);
    const auto *kl_721 = buffer.data(kl + 721);
    const auto *kl_722 = buffer.data(kl + 722);
    const auto *kl_723 = buffer.data(kl + 723);
    const auto *kl_724 = buffer.data(kl + 724);
    const auto *kl_725 = buffer.data(kl + 725);
    const auto *kl_726 = buffer.data(kl + 726);
    const auto *kl_727 = buffer.data(kl + 727);
    const auto *kl_728 = buffer.data(kl + 728);
    const auto *kl_729 = buffer.data(kl + 729);
    const auto *kl_730 = buffer.data(kl + 730);
    const auto *kl_731 = buffer.data(kl + 731);

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, hl_210, hl_211, hl_212, hl_213, \
                         hl_214, kl_525, kl_526, kl_527, kl_528, \
                         kl_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * hl_210[k]
                   + f_0 * kl_525[k];

        t_346[k] = -2.0 * hl_211[k]
                   + f_0 * kl_526[k];

        t_347[k] = -2.0 * hl_212[k]
                   + f_0 * kl_527[k];

        t_348[k] = -2.0 * hl_213[k]
                   + f_0 * kl_528[k];

        t_349[k] = -2.0 * hl_214[k]
                   + f_0 * kl_529[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, hl_215, hl_216, hl_217, hl_218, \
                         hl_219, kl_530, kl_531, kl_532, kl_533, \
                         kl_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * hl_215[k]
                   + f_0 * kl_530[k];

        t_351[k] = -2.0 * hl_216[k]
                   + f_0 * kl_531[k];

        t_352[k] = -2.0 * hl_217[k]
                   + f_0 * kl_532[k];

        t_353[k] = -2.0 * hl_218[k]
                   + f_0 * kl_533[k];

        t_354[k] = -2.0 * hl_219[k]
                   + f_0 * kl_534[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, hl_220, hl_221, hl_222, hl_223, \
                         hl_224, kl_535, kl_536, kl_537, kl_538, \
                         kl_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * hl_220[k]
                   + f_0 * kl_535[k];

        t_356[k] = -2.0 * hl_221[k]
                   + f_0 * kl_536[k];

        t_357[k] = -2.0 * hl_222[k]
                   + f_0 * kl_537[k];

        t_358[k] = -2.0 * hl_223[k]
                   + f_0 * kl_538[k];

        t_359[k] = -2.0 * hl_224[k]
                   + f_0 * kl_539[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, hl_225, hl_226, hl_227, hl_228, \
                         hl_229, kl_540, kl_541, kl_542, kl_543, \
                         kl_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -hl_225[k]
                   + f_0 * kl_540[k];

        t_361[k] = -hl_226[k]
                   + f_0 * kl_541[k];

        t_362[k] = -hl_227[k]
                   + f_0 * kl_542[k];

        t_363[k] = -hl_228[k]
                   + f_0 * kl_543[k];

        t_364[k] = -hl_229[k]
                   + f_0 * kl_544[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, hl_230, hl_231, hl_232, hl_233, \
                         hl_234, kl_545, kl_546, kl_547, kl_548, \
                         kl_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -hl_230[k]
                   + f_0 * kl_545[k];

        t_366[k] = -hl_231[k]
                   + f_0 * kl_546[k];

        t_367[k] = -hl_232[k]
                   + f_0 * kl_547[k];

        t_368[k] = -hl_233[k]
                   + f_0 * kl_548[k];

        t_369[k] = -hl_234[k]
                   + f_0 * kl_549[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, hl_235, hl_236, hl_237, hl_238, \
                         hl_239, kl_550, kl_551, kl_552, kl_553, \
                         kl_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -hl_235[k]
                   + f_0 * kl_550[k];

        t_371[k] = -hl_236[k]
                   + f_0 * kl_551[k];

        t_372[k] = -hl_237[k]
                   + f_0 * kl_552[k];

        t_373[k] = -hl_238[k]
                   + f_0 * kl_553[k];

        t_374[k] = -hl_239[k]
                   + f_0 * kl_554[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, hl_240, hl_241, hl_242, hl_243, \
                         hl_244, kl_555, kl_556, kl_557, kl_558, \
                         kl_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -hl_240[k]
                   + f_0 * kl_555[k];

        t_376[k] = -hl_241[k]
                   + f_0 * kl_556[k];

        t_377[k] = -hl_242[k]
                   + f_0 * kl_557[k];

        t_378[k] = -hl_243[k]
                   + f_0 * kl_558[k];

        t_379[k] = -hl_244[k]
                   + f_0 * kl_559[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, hl_245, hl_246, hl_247, hl_248, \
                         hl_249, kl_560, kl_561, kl_562, kl_563, \
                         kl_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -hl_245[k]
                   + f_0 * kl_560[k];

        t_381[k] = -hl_246[k]
                   + f_0 * kl_561[k];

        t_382[k] = -hl_247[k]
                   + f_0 * kl_562[k];

        t_383[k] = -hl_248[k]
                   + f_0 * kl_563[k];

        t_384[k] = -hl_249[k]
                   + f_0 * kl_564[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, hl_250, hl_251, hl_252, hl_253, \
                         hl_254, kl_565, kl_566, kl_567, kl_568, \
                         kl_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -hl_250[k]
                   + f_0 * kl_565[k];

        t_386[k] = -hl_251[k]
                   + f_0 * kl_566[k];

        t_387[k] = -hl_252[k]
                   + f_0 * kl_567[k];

        t_388[k] = -hl_253[k]
                   + f_0 * kl_568[k];

        t_389[k] = -hl_254[k]
                   + f_0 * kl_569[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, hl_255, hl_256, hl_257, hl_258, \
                         hl_259, kl_570, kl_571, kl_572, kl_573, \
                         kl_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -hl_255[k]
                   + f_0 * kl_570[k];

        t_391[k] = -hl_256[k]
                   + f_0 * kl_571[k];

        t_392[k] = -hl_257[k]
                   + f_0 * kl_572[k];

        t_393[k] = -hl_258[k]
                   + f_0 * kl_573[k];

        t_394[k] = -hl_259[k]
                   + f_0 * kl_574[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, hl_260, hl_261, hl_262, hl_263, \
                         hl_264, kl_575, kl_576, kl_577, kl_578, \
                         kl_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -hl_260[k]
                   + f_0 * kl_575[k];

        t_396[k] = -hl_261[k]
                   + f_0 * kl_576[k];

        t_397[k] = -hl_262[k]
                   + f_0 * kl_577[k];

        t_398[k] = -hl_263[k]
                   + f_0 * kl_578[k];

        t_399[k] = -hl_264[k]
                   + f_0 * kl_579[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, hl_265, hl_266, hl_267, hl_268, \
                         hl_269, kl_580, kl_581, kl_582, kl_583, \
                         kl_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -hl_265[k]
                   + f_0 * kl_580[k];

        t_401[k] = -hl_266[k]
                   + f_0 * kl_581[k];

        t_402[k] = -hl_267[k]
                   + f_0 * kl_582[k];

        t_403[k] = -hl_268[k]
                   + f_0 * kl_583[k];

        t_404[k] = -hl_269[k]
                   + f_0 * kl_584[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, kl_585, \
                         kl_586, kl_587, kl_588, kl_589, kl_590, kl_591, \
                         kl_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * kl_585[k];

        t_406[k] = f_0 * kl_586[k];

        t_407[k] = f_0 * kl_587[k];

        t_408[k] = f_0 * kl_588[k];

        t_409[k] = f_0 * kl_589[k];

        t_410[k] = f_0 * kl_590[k];

        t_411[k] = f_0 * kl_591[k];

        t_412[k] = f_0 * kl_592[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, t_420, kl_593, \
                         kl_594, kl_595, kl_596, kl_597, kl_598, kl_599, \
                         kl_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * kl_593[k];

        t_414[k] = f_0 * kl_594[k];

        t_415[k] = f_0 * kl_595[k];

        t_416[k] = f_0 * kl_596[k];

        t_417[k] = f_0 * kl_597[k];

        t_418[k] = f_0 * kl_598[k];

        t_419[k] = f_0 * kl_599[k];

        t_420[k] = f_0 * kl_600[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, t_425, t_426, t_427, t_428, kl_601, \
                         kl_602, kl_603, kl_604, kl_605, kl_606, kl_607, \
                         kl_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_0 * kl_601[k];

        t_422[k] = f_0 * kl_602[k];

        t_423[k] = f_0 * kl_603[k];

        t_424[k] = f_0 * kl_604[k];

        t_425[k] = f_0 * kl_605[k];

        t_426[k] = f_0 * kl_606[k];

        t_427[k] = f_0 * kl_607[k];

        t_428[k] = f_0 * kl_608[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, t_434, t_435, t_436, kl_609, \
                         kl_610, kl_611, kl_612, kl_613, kl_614, kl_615, \
                         kl_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_0 * kl_609[k];

        t_430[k] = f_0 * kl_610[k];

        t_431[k] = f_0 * kl_611[k];

        t_432[k] = f_0 * kl_612[k];

        t_433[k] = f_0 * kl_613[k];

        t_434[k] = f_0 * kl_614[k];

        t_435[k] = f_0 * kl_615[k];

        t_436[k] = f_0 * kl_616[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, t_441, t_442, t_443, t_444, kl_617, \
                         kl_618, kl_619, kl_620, kl_621, kl_622, kl_623, \
                         kl_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_0 * kl_617[k];

        t_438[k] = f_0 * kl_618[k];

        t_439[k] = f_0 * kl_619[k];

        t_440[k] = f_0 * kl_620[k];

        t_441[k] = f_0 * kl_621[k];

        t_442[k] = f_0 * kl_622[k];

        t_443[k] = f_0 * kl_623[k];

        t_444[k] = f_0 * kl_624[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, t_450, t_451, hl_270, hl_271, \
                         kl_625, kl_626, kl_627, kl_628, kl_629, kl_675, \
                         kl_676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_0 * kl_625[k];

        t_446[k] = f_0 * kl_626[k];

        t_447[k] = f_0 * kl_627[k];

        t_448[k] = f_0 * kl_628[k];

        t_449[k] = f_0 * kl_629[k];

        t_450[k] = -4.0 * hl_270[k]
                   + f_0 * kl_675[k];

        t_451[k] = -4.0 * hl_271[k]
                   + f_0 * kl_676[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, hl_272, hl_273, hl_274, hl_275, \
                         hl_276, kl_677, kl_678, kl_679, kl_680, \
                         kl_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = -4.0 * hl_272[k]
                   + f_0 * kl_677[k];

        t_453[k] = -4.0 * hl_273[k]
                   + f_0 * kl_678[k];

        t_454[k] = -4.0 * hl_274[k]
                   + f_0 * kl_679[k];

        t_455[k] = -4.0 * hl_275[k]
                   + f_0 * kl_680[k];

        t_456[k] = -4.0 * hl_276[k]
                   + f_0 * kl_681[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, t_461, hl_277, hl_278, hl_279, hl_280, \
                         hl_281, kl_682, kl_683, kl_684, kl_685, \
                         kl_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = -4.0 * hl_277[k]
                   + f_0 * kl_682[k];

        t_458[k] = -4.0 * hl_278[k]
                   + f_0 * kl_683[k];

        t_459[k] = -4.0 * hl_279[k]
                   + f_0 * kl_684[k];

        t_460[k] = -4.0 * hl_280[k]
                   + f_0 * kl_685[k];

        t_461[k] = -4.0 * hl_281[k]
                   + f_0 * kl_686[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, hl_282, hl_283, hl_284, hl_285, \
                         hl_286, kl_687, kl_688, kl_689, kl_690, \
                         kl_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = -4.0 * hl_282[k]
                   + f_0 * kl_687[k];

        t_463[k] = -4.0 * hl_283[k]
                   + f_0 * kl_688[k];

        t_464[k] = -4.0 * hl_284[k]
                   + f_0 * kl_689[k];

        t_465[k] = -4.0 * hl_285[k]
                   + f_0 * kl_690[k];

        t_466[k] = -4.0 * hl_286[k]
                   + f_0 * kl_691[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, hl_287, hl_288, hl_289, hl_290, \
                         hl_291, kl_692, kl_693, kl_694, kl_695, \
                         kl_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = -4.0 * hl_287[k]
                   + f_0 * kl_692[k];

        t_468[k] = -4.0 * hl_288[k]
                   + f_0 * kl_693[k];

        t_469[k] = -4.0 * hl_289[k]
                   + f_0 * kl_694[k];

        t_470[k] = -4.0 * hl_290[k]
                   + f_0 * kl_695[k];

        t_471[k] = -4.0 * hl_291[k]
                   + f_0 * kl_696[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, hl_292, hl_293, hl_294, hl_295, \
                         hl_296, kl_697, kl_698, kl_699, kl_700, \
                         kl_701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = -4.0 * hl_292[k]
                   + f_0 * kl_697[k];

        t_473[k] = -4.0 * hl_293[k]
                   + f_0 * kl_698[k];

        t_474[k] = -4.0 * hl_294[k]
                   + f_0 * kl_699[k];

        t_475[k] = -4.0 * hl_295[k]
                   + f_0 * kl_700[k];

        t_476[k] = -4.0 * hl_296[k]
                   + f_0 * kl_701[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, hl_297, hl_298, hl_299, hl_300, \
                         hl_301, kl_702, kl_703, kl_704, kl_705, \
                         kl_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = -4.0 * hl_297[k]
                   + f_0 * kl_702[k];

        t_478[k] = -4.0 * hl_298[k]
                   + f_0 * kl_703[k];

        t_479[k] = -4.0 * hl_299[k]
                   + f_0 * kl_704[k];

        t_480[k] = -4.0 * hl_300[k]
                   + f_0 * kl_705[k];

        t_481[k] = -4.0 * hl_301[k]
                   + f_0 * kl_706[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, t_486, hl_302, hl_303, hl_304, hl_305, \
                         hl_306, kl_707, kl_708, kl_709, kl_710, \
                         kl_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = -4.0 * hl_302[k]
                   + f_0 * kl_707[k];

        t_483[k] = -4.0 * hl_303[k]
                   + f_0 * kl_708[k];

        t_484[k] = -4.0 * hl_304[k]
                   + f_0 * kl_709[k];

        t_485[k] = -4.0 * hl_305[k]
                   + f_0 * kl_710[k];

        t_486[k] = -4.0 * hl_306[k]
                   + f_0 * kl_711[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, hl_307, hl_308, hl_309, hl_310, \
                         hl_311, kl_712, kl_713, kl_714, kl_715, \
                         kl_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = -4.0 * hl_307[k]
                   + f_0 * kl_712[k];

        t_488[k] = -4.0 * hl_308[k]
                   + f_0 * kl_713[k];

        t_489[k] = -4.0 * hl_309[k]
                   + f_0 * kl_714[k];

        t_490[k] = -4.0 * hl_310[k]
                   + f_0 * kl_715[k];

        t_491[k] = -4.0 * hl_311[k]
                   + f_0 * kl_716[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, hl_312, hl_313, hl_314, hl_315, \
                         hl_316, kl_717, kl_718, kl_719, kl_720, \
                         kl_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = -4.0 * hl_312[k]
                   + f_0 * kl_717[k];

        t_493[k] = -4.0 * hl_313[k]
                   + f_0 * kl_718[k];

        t_494[k] = -4.0 * hl_314[k]
                   + f_0 * kl_719[k];

        t_495[k] = -3.0 * hl_315[k]
                   + f_0 * kl_720[k];

        t_496[k] = -3.0 * hl_316[k]
                   + f_0 * kl_721[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, t_500, t_501, hl_317, hl_318, hl_319, hl_320, \
                         hl_321, kl_722, kl_723, kl_724, kl_725, \
                         kl_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = -3.0 * hl_317[k]
                   + f_0 * kl_722[k];

        t_498[k] = -3.0 * hl_318[k]
                   + f_0 * kl_723[k];

        t_499[k] = -3.0 * hl_319[k]
                   + f_0 * kl_724[k];

        t_500[k] = -3.0 * hl_320[k]
                   + f_0 * kl_725[k];

        t_501[k] = -3.0 * hl_321[k]
                   + f_0 * kl_726[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, hl_322, hl_323, hl_324, hl_325, \
                         hl_326, kl_727, kl_728, kl_729, kl_730, \
                         kl_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = -3.0 * hl_322[k]
                   + f_0 * kl_727[k];

        t_503[k] = -3.0 * hl_323[k]
                   + f_0 * kl_728[k];

        t_504[k] = -3.0 * hl_324[k]
                   + f_0 * kl_729[k];

        t_505[k] = -3.0 * hl_325[k]
                   + f_0 * kl_730[k];

        t_506[k] = -3.0 * hl_326[k]
                   + f_0 * kl_731[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hl_327 = buffer.data(hl + 327);
    const auto *hl_328 = buffer.data(hl + 328);
    const auto *hl_329 = buffer.data(hl + 329);
    const auto *hl_330 = buffer.data(hl + 330);
    const auto *hl_331 = buffer.data(hl + 331);
    const auto *hl_332 = buffer.data(hl + 332);
    const auto *hl_333 = buffer.data(hl + 333);
    const auto *hl_334 = buffer.data(hl + 334);
    const auto *hl_335 = buffer.data(hl + 335);
    const auto *hl_336 = buffer.data(hl + 336);
    const auto *hl_337 = buffer.data(hl + 337);
    const auto *hl_338 = buffer.data(hl + 338);
    const auto *hl_339 = buffer.data(hl + 339);
    const auto *hl_340 = buffer.data(hl + 340);
    const auto *hl_341 = buffer.data(hl + 341);
    const auto *hl_342 = buffer.data(hl + 342);
    const auto *hl_343 = buffer.data(hl + 343);
    const auto *hl_344 = buffer.data(hl + 344);
    const auto *hl_345 = buffer.data(hl + 345);
    const auto *hl_346 = buffer.data(hl + 346);
    const auto *hl_347 = buffer.data(hl + 347);
    const auto *hl_348 = buffer.data(hl + 348);
    const auto *hl_349 = buffer.data(hl + 349);
    const auto *hl_350 = buffer.data(hl + 350);
    const auto *hl_351 = buffer.data(hl + 351);
    const auto *hl_352 = buffer.data(hl + 352);
    const auto *hl_353 = buffer.data(hl + 353);
    const auto *hl_354 = buffer.data(hl + 354);
    const auto *hl_355 = buffer.data(hl + 355);
    const auto *hl_356 = buffer.data(hl + 356);
    const auto *hl_357 = buffer.data(hl + 357);
    const auto *hl_358 = buffer.data(hl + 358);
    const auto *hl_359 = buffer.data(hl + 359);
    const auto *hl_360 = buffer.data(hl + 360);
    const auto *hl_361 = buffer.data(hl + 361);
    const auto *hl_362 = buffer.data(hl + 362);
    const auto *hl_363 = buffer.data(hl + 363);
    const auto *hl_364 = buffer.data(hl + 364);
    const auto *hl_365 = buffer.data(hl + 365);
    const auto *hl_366 = buffer.data(hl + 366);
    const auto *hl_367 = buffer.data(hl + 367);
    const auto *hl_368 = buffer.data(hl + 368);
    const auto *hl_369 = buffer.data(hl + 369);
    const auto *hl_370 = buffer.data(hl + 370);
    const auto *hl_371 = buffer.data(hl + 371);
    const auto *hl_372 = buffer.data(hl + 372);
    const auto *hl_373 = buffer.data(hl + 373);
    const auto *hl_374 = buffer.data(hl + 374);
    const auto *hl_375 = buffer.data(hl + 375);
    const auto *hl_376 = buffer.data(hl + 376);
    const auto *hl_377 = buffer.data(hl + 377);
    const auto *hl_378 = buffer.data(hl + 378);
    const auto *hl_379 = buffer.data(hl + 379);
    const auto *hl_380 = buffer.data(hl + 380);
    const auto *hl_381 = buffer.data(hl + 381);
    const auto *hl_382 = buffer.data(hl + 382);
    const auto *hl_383 = buffer.data(hl + 383);
    const auto *hl_384 = buffer.data(hl + 384);
    const auto *hl_385 = buffer.data(hl + 385);
    const auto *hl_386 = buffer.data(hl + 386);
    const auto *hl_387 = buffer.data(hl + 387);
    const auto *hl_388 = buffer.data(hl + 388);
    const auto *hl_389 = buffer.data(hl + 389);
    const auto *hl_390 = buffer.data(hl + 390);
    const auto *hl_391 = buffer.data(hl + 391);
    const auto *hl_392 = buffer.data(hl + 392);
    const auto *hl_393 = buffer.data(hl + 393);
    const auto *hl_394 = buffer.data(hl + 394);
    const auto *hl_395 = buffer.data(hl + 395);
    const auto *hl_396 = buffer.data(hl + 396);
    const auto *hl_397 = buffer.data(hl + 397);
    const auto *hl_398 = buffer.data(hl + 398);
    const auto *hl_399 = buffer.data(hl + 399);
    const auto *hl_400 = buffer.data(hl + 400);
    const auto *hl_401 = buffer.data(hl + 401);
    const auto *hl_402 = buffer.data(hl + 402);
    const auto *hl_403 = buffer.data(hl + 403);
    const auto *hl_404 = buffer.data(hl + 404);
    const auto *hl_405 = buffer.data(hl + 405);
    const auto *hl_406 = buffer.data(hl + 406);
    const auto *hl_407 = buffer.data(hl + 407);
    const auto *hl_408 = buffer.data(hl + 408);
    const auto *hl_409 = buffer.data(hl + 409);
    const auto *hl_410 = buffer.data(hl + 410);
    const auto *hl_411 = buffer.data(hl + 411);
    const auto *hl_412 = buffer.data(hl + 412);
    const auto *hl_413 = buffer.data(hl + 413);
    const auto *hl_414 = buffer.data(hl + 414);
    const auto *hl_415 = buffer.data(hl + 415);
    const auto *hl_416 = buffer.data(hl + 416);
    const auto *hl_417 = buffer.data(hl + 417);
    const auto *hl_418 = buffer.data(hl + 418);
    const auto *hl_419 = buffer.data(hl + 419);
    const auto *hl_420 = buffer.data(hl + 420);
    const auto *hl_421 = buffer.data(hl + 421);
    const auto *hl_422 = buffer.data(hl + 422);
    const auto *hl_423 = buffer.data(hl + 423);
    const auto *hl_424 = buffer.data(hl + 424);
    const auto *hl_425 = buffer.data(hl + 425);
    const auto *hl_426 = buffer.data(hl + 426);
    const auto *hl_427 = buffer.data(hl + 427);
    const auto *hl_428 = buffer.data(hl + 428);
    const auto *hl_429 = buffer.data(hl + 429);
    const auto *hl_430 = buffer.data(hl + 430);
    const auto *hl_431 = buffer.data(hl + 431);
    const auto *hl_432 = buffer.data(hl + 432);
    const auto *hl_433 = buffer.data(hl + 433);
    const auto *hl_434 = buffer.data(hl + 434);
    const auto *hl_435 = buffer.data(hl + 435);
    const auto *hl_436 = buffer.data(hl + 436);
    const auto *hl_437 = buffer.data(hl + 437);
    const auto *hl_438 = buffer.data(hl + 438);
    const auto *hl_439 = buffer.data(hl + 439);
    const auto *hl_440 = buffer.data(hl + 440);
    const auto *hl_441 = buffer.data(hl + 441);
    const auto *hl_442 = buffer.data(hl + 442);
    const auto *hl_443 = buffer.data(hl + 443);
    const auto *hl_444 = buffer.data(hl + 444);
    const auto *hl_445 = buffer.data(hl + 445);
    const auto *hl_446 = buffer.data(hl + 446);
    const auto *hl_447 = buffer.data(hl + 447);
    const auto *hl_448 = buffer.data(hl + 448);
    const auto *hl_449 = buffer.data(hl + 449);

    const auto *kl_732 = buffer.data(kl + 732);
    const auto *kl_733 = buffer.data(kl + 733);
    const auto *kl_734 = buffer.data(kl + 734);
    const auto *kl_735 = buffer.data(kl + 735);
    const auto *kl_736 = buffer.data(kl + 736);
    const auto *kl_737 = buffer.data(kl + 737);
    const auto *kl_738 = buffer.data(kl + 738);
    const auto *kl_739 = buffer.data(kl + 739);
    const auto *kl_740 = buffer.data(kl + 740);
    const auto *kl_741 = buffer.data(kl + 741);
    const auto *kl_742 = buffer.data(kl + 742);
    const auto *kl_743 = buffer.data(kl + 743);
    const auto *kl_744 = buffer.data(kl + 744);
    const auto *kl_745 = buffer.data(kl + 745);
    const auto *kl_746 = buffer.data(kl + 746);
    const auto *kl_747 = buffer.data(kl + 747);
    const auto *kl_748 = buffer.data(kl + 748);
    const auto *kl_749 = buffer.data(kl + 749);
    const auto *kl_750 = buffer.data(kl + 750);
    const auto *kl_751 = buffer.data(kl + 751);
    const auto *kl_752 = buffer.data(kl + 752);
    const auto *kl_753 = buffer.data(kl + 753);
    const auto *kl_754 = buffer.data(kl + 754);
    const auto *kl_755 = buffer.data(kl + 755);
    const auto *kl_756 = buffer.data(kl + 756);
    const auto *kl_757 = buffer.data(kl + 757);
    const auto *kl_758 = buffer.data(kl + 758);
    const auto *kl_759 = buffer.data(kl + 759);
    const auto *kl_760 = buffer.data(kl + 760);
    const auto *kl_761 = buffer.data(kl + 761);
    const auto *kl_762 = buffer.data(kl + 762);
    const auto *kl_763 = buffer.data(kl + 763);
    const auto *kl_764 = buffer.data(kl + 764);
    const auto *kl_765 = buffer.data(kl + 765);
    const auto *kl_766 = buffer.data(kl + 766);
    const auto *kl_767 = buffer.data(kl + 767);
    const auto *kl_768 = buffer.data(kl + 768);
    const auto *kl_769 = buffer.data(kl + 769);
    const auto *kl_770 = buffer.data(kl + 770);
    const auto *kl_771 = buffer.data(kl + 771);
    const auto *kl_772 = buffer.data(kl + 772);
    const auto *kl_773 = buffer.data(kl + 773);
    const auto *kl_774 = buffer.data(kl + 774);
    const auto *kl_775 = buffer.data(kl + 775);
    const auto *kl_776 = buffer.data(kl + 776);
    const auto *kl_777 = buffer.data(kl + 777);
    const auto *kl_778 = buffer.data(kl + 778);
    const auto *kl_779 = buffer.data(kl + 779);
    const auto *kl_780 = buffer.data(kl + 780);
    const auto *kl_781 = buffer.data(kl + 781);
    const auto *kl_782 = buffer.data(kl + 782);
    const auto *kl_783 = buffer.data(kl + 783);
    const auto *kl_784 = buffer.data(kl + 784);
    const auto *kl_785 = buffer.data(kl + 785);
    const auto *kl_786 = buffer.data(kl + 786);
    const auto *kl_787 = buffer.data(kl + 787);
    const auto *kl_788 = buffer.data(kl + 788);
    const auto *kl_789 = buffer.data(kl + 789);
    const auto *kl_790 = buffer.data(kl + 790);
    const auto *kl_791 = buffer.data(kl + 791);
    const auto *kl_792 = buffer.data(kl + 792);
    const auto *kl_793 = buffer.data(kl + 793);
    const auto *kl_794 = buffer.data(kl + 794);
    const auto *kl_795 = buffer.data(kl + 795);
    const auto *kl_796 = buffer.data(kl + 796);
    const auto *kl_797 = buffer.data(kl + 797);
    const auto *kl_798 = buffer.data(kl + 798);
    const auto *kl_799 = buffer.data(kl + 799);
    const auto *kl_800 = buffer.data(kl + 800);
    const auto *kl_801 = buffer.data(kl + 801);
    const auto *kl_802 = buffer.data(kl + 802);
    const auto *kl_803 = buffer.data(kl + 803);
    const auto *kl_804 = buffer.data(kl + 804);
    const auto *kl_805 = buffer.data(kl + 805);
    const auto *kl_806 = buffer.data(kl + 806);
    const auto *kl_807 = buffer.data(kl + 807);
    const auto *kl_808 = buffer.data(kl + 808);
    const auto *kl_809 = buffer.data(kl + 809);
    const auto *kl_810 = buffer.data(kl + 810);
    const auto *kl_811 = buffer.data(kl + 811);
    const auto *kl_812 = buffer.data(kl + 812);
    const auto *kl_813 = buffer.data(kl + 813);
    const auto *kl_814 = buffer.data(kl + 814);
    const auto *kl_815 = buffer.data(kl + 815);
    const auto *kl_816 = buffer.data(kl + 816);
    const auto *kl_817 = buffer.data(kl + 817);
    const auto *kl_818 = buffer.data(kl + 818);
    const auto *kl_819 = buffer.data(kl + 819);
    const auto *kl_820 = buffer.data(kl + 820);
    const auto *kl_821 = buffer.data(kl + 821);
    const auto *kl_822 = buffer.data(kl + 822);
    const auto *kl_823 = buffer.data(kl + 823);
    const auto *kl_824 = buffer.data(kl + 824);
    const auto *kl_825 = buffer.data(kl + 825);
    const auto *kl_826 = buffer.data(kl + 826);
    const auto *kl_827 = buffer.data(kl + 827);
    const auto *kl_828 = buffer.data(kl + 828);
    const auto *kl_829 = buffer.data(kl + 829);
    const auto *kl_830 = buffer.data(kl + 830);
    const auto *kl_831 = buffer.data(kl + 831);
    const auto *kl_832 = buffer.data(kl + 832);
    const auto *kl_833 = buffer.data(kl + 833);
    const auto *kl_834 = buffer.data(kl + 834);
    const auto *kl_835 = buffer.data(kl + 835);
    const auto *kl_836 = buffer.data(kl + 836);
    const auto *kl_837 = buffer.data(kl + 837);
    const auto *kl_838 = buffer.data(kl + 838);
    const auto *kl_839 = buffer.data(kl + 839);
    const auto *kl_840 = buffer.data(kl + 840);
    const auto *kl_841 = buffer.data(kl + 841);
    const auto *kl_842 = buffer.data(kl + 842);
    const auto *kl_843 = buffer.data(kl + 843);
    const auto *kl_844 = buffer.data(kl + 844);
    const auto *kl_845 = buffer.data(kl + 845);
    const auto *kl_846 = buffer.data(kl + 846);
    const auto *kl_847 = buffer.data(kl + 847);
    const auto *kl_848 = buffer.data(kl + 848);
    const auto *kl_849 = buffer.data(kl + 849);
    const auto *kl_850 = buffer.data(kl + 850);
    const auto *kl_851 = buffer.data(kl + 851);
    const auto *kl_852 = buffer.data(kl + 852);
    const auto *kl_853 = buffer.data(kl + 853);
    const auto *kl_854 = buffer.data(kl + 854);
    const auto *kl_855 = buffer.data(kl + 855);
    const auto *kl_856 = buffer.data(kl + 856);
    const auto *kl_857 = buffer.data(kl + 857);
    const auto *kl_858 = buffer.data(kl + 858);
    const auto *kl_859 = buffer.data(kl + 859);
    const auto *kl_860 = buffer.data(kl + 860);
    const auto *kl_861 = buffer.data(kl + 861);
    const auto *kl_862 = buffer.data(kl + 862);
    const auto *kl_863 = buffer.data(kl + 863);
    const auto *kl_864 = buffer.data(kl + 864);
    const auto *kl_865 = buffer.data(kl + 865);
    const auto *kl_866 = buffer.data(kl + 866);
    const auto *kl_867 = buffer.data(kl + 867);
    const auto *kl_868 = buffer.data(kl + 868);
    const auto *kl_869 = buffer.data(kl + 869);
    const auto *kl_870 = buffer.data(kl + 870);
    const auto *kl_871 = buffer.data(kl + 871);
    const auto *kl_872 = buffer.data(kl + 872);
    const auto *kl_873 = buffer.data(kl + 873);
    const auto *kl_874 = buffer.data(kl + 874);
    const auto *kl_875 = buffer.data(kl + 875);
    const auto *kl_876 = buffer.data(kl + 876);
    const auto *kl_877 = buffer.data(kl + 877);
    const auto *kl_878 = buffer.data(kl + 878);
    const auto *kl_879 = buffer.data(kl + 879);
    const auto *kl_880 = buffer.data(kl + 880);
    const auto *kl_881 = buffer.data(kl + 881);
    const auto *kl_882 = buffer.data(kl + 882);
    const auto *kl_883 = buffer.data(kl + 883);
    const auto *kl_884 = buffer.data(kl + 884);
    const auto *kl_885 = buffer.data(kl + 885);
    const auto *kl_886 = buffer.data(kl + 886);
    const auto *kl_887 = buffer.data(kl + 887);
    const auto *kl_888 = buffer.data(kl + 888);
    const auto *kl_889 = buffer.data(kl + 889);
    const auto *kl_890 = buffer.data(kl + 890);
    const auto *kl_891 = buffer.data(kl + 891);
    const auto *kl_892 = buffer.data(kl + 892);
    const auto *kl_893 = buffer.data(kl + 893);
    const auto *kl_894 = buffer.data(kl + 894);
    const auto *kl_895 = buffer.data(kl + 895);
    const auto *kl_896 = buffer.data(kl + 896);
    const auto *kl_897 = buffer.data(kl + 897);

#pragma omp simd aligned(t_507, t_508, t_509, t_510, t_511, hl_327, hl_328, hl_329, hl_330, \
                         hl_331, kl_732, kl_733, kl_734, kl_735, \
                         kl_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = -3.0 * hl_327[k]
                   + f_0 * kl_732[k];

        t_508[k] = -3.0 * hl_328[k]
                   + f_0 * kl_733[k];

        t_509[k] = -3.0 * hl_329[k]
                   + f_0 * kl_734[k];

        t_510[k] = -3.0 * hl_330[k]
                   + f_0 * kl_735[k];

        t_511[k] = -3.0 * hl_331[k]
                   + f_0 * kl_736[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, hl_332, hl_333, hl_334, hl_335, \
                         hl_336, kl_737, kl_738, kl_739, kl_740, \
                         kl_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = -3.0 * hl_332[k]
                   + f_0 * kl_737[k];

        t_513[k] = -3.0 * hl_333[k]
                   + f_0 * kl_738[k];

        t_514[k] = -3.0 * hl_334[k]
                   + f_0 * kl_739[k];

        t_515[k] = -3.0 * hl_335[k]
                   + f_0 * kl_740[k];

        t_516[k] = -3.0 * hl_336[k]
                   + f_0 * kl_741[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, hl_337, hl_338, hl_339, hl_340, \
                         hl_341, kl_742, kl_743, kl_744, kl_745, \
                         kl_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = -3.0 * hl_337[k]
                   + f_0 * kl_742[k];

        t_518[k] = -3.0 * hl_338[k]
                   + f_0 * kl_743[k];

        t_519[k] = -3.0 * hl_339[k]
                   + f_0 * kl_744[k];

        t_520[k] = -3.0 * hl_340[k]
                   + f_0 * kl_745[k];

        t_521[k] = -3.0 * hl_341[k]
                   + f_0 * kl_746[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, hl_342, hl_343, hl_344, hl_345, \
                         hl_346, kl_747, kl_748, kl_749, kl_750, \
                         kl_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = -3.0 * hl_342[k]
                   + f_0 * kl_747[k];

        t_523[k] = -3.0 * hl_343[k]
                   + f_0 * kl_748[k];

        t_524[k] = -3.0 * hl_344[k]
                   + f_0 * kl_749[k];

        t_525[k] = -3.0 * hl_345[k]
                   + f_0 * kl_750[k];

        t_526[k] = -3.0 * hl_346[k]
                   + f_0 * kl_751[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, hl_347, hl_348, hl_349, hl_350, \
                         hl_351, kl_752, kl_753, kl_754, kl_755, \
                         kl_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = -3.0 * hl_347[k]
                   + f_0 * kl_752[k];

        t_528[k] = -3.0 * hl_348[k]
                   + f_0 * kl_753[k];

        t_529[k] = -3.0 * hl_349[k]
                   + f_0 * kl_754[k];

        t_530[k] = -3.0 * hl_350[k]
                   + f_0 * kl_755[k];

        t_531[k] = -3.0 * hl_351[k]
                   + f_0 * kl_756[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, hl_352, hl_353, hl_354, hl_355, \
                         hl_356, kl_757, kl_758, kl_759, kl_760, \
                         kl_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = -3.0 * hl_352[k]
                   + f_0 * kl_757[k];

        t_533[k] = -3.0 * hl_353[k]
                   + f_0 * kl_758[k];

        t_534[k] = -3.0 * hl_354[k]
                   + f_0 * kl_759[k];

        t_535[k] = -3.0 * hl_355[k]
                   + f_0 * kl_760[k];

        t_536[k] = -3.0 * hl_356[k]
                   + f_0 * kl_761[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, t_541, hl_357, hl_358, hl_359, hl_360, \
                         hl_361, kl_762, kl_763, kl_764, kl_765, \
                         kl_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = -3.0 * hl_357[k]
                   + f_0 * kl_762[k];

        t_538[k] = -3.0 * hl_358[k]
                   + f_0 * kl_763[k];

        t_539[k] = -3.0 * hl_359[k]
                   + f_0 * kl_764[k];

        t_540[k] = -2.0 * hl_360[k]
                   + f_0 * kl_765[k];

        t_541[k] = -2.0 * hl_361[k]
                   + f_0 * kl_766[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, hl_362, hl_363, hl_364, hl_365, \
                         hl_366, kl_767, kl_768, kl_769, kl_770, \
                         kl_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = -2.0 * hl_362[k]
                   + f_0 * kl_767[k];

        t_543[k] = -2.0 * hl_363[k]
                   + f_0 * kl_768[k];

        t_544[k] = -2.0 * hl_364[k]
                   + f_0 * kl_769[k];

        t_545[k] = -2.0 * hl_365[k]
                   + f_0 * kl_770[k];

        t_546[k] = -2.0 * hl_366[k]
                   + f_0 * kl_771[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, hl_367, hl_368, hl_369, hl_370, \
                         hl_371, kl_772, kl_773, kl_774, kl_775, \
                         kl_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = -2.0 * hl_367[k]
                   + f_0 * kl_772[k];

        t_548[k] = -2.0 * hl_368[k]
                   + f_0 * kl_773[k];

        t_549[k] = -2.0 * hl_369[k]
                   + f_0 * kl_774[k];

        t_550[k] = -2.0 * hl_370[k]
                   + f_0 * kl_775[k];

        t_551[k] = -2.0 * hl_371[k]
                   + f_0 * kl_776[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, t_556, hl_372, hl_373, hl_374, hl_375, \
                         hl_376, kl_777, kl_778, kl_779, kl_780, \
                         kl_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = -2.0 * hl_372[k]
                   + f_0 * kl_777[k];

        t_553[k] = -2.0 * hl_373[k]
                   + f_0 * kl_778[k];

        t_554[k] = -2.0 * hl_374[k]
                   + f_0 * kl_779[k];

        t_555[k] = -2.0 * hl_375[k]
                   + f_0 * kl_780[k];

        t_556[k] = -2.0 * hl_376[k]
                   + f_0 * kl_781[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, hl_377, hl_378, hl_379, hl_380, \
                         hl_381, kl_782, kl_783, kl_784, kl_785, \
                         kl_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = -2.0 * hl_377[k]
                   + f_0 * kl_782[k];

        t_558[k] = -2.0 * hl_378[k]
                   + f_0 * kl_783[k];

        t_559[k] = -2.0 * hl_379[k]
                   + f_0 * kl_784[k];

        t_560[k] = -2.0 * hl_380[k]
                   + f_0 * kl_785[k];

        t_561[k] = -2.0 * hl_381[k]
                   + f_0 * kl_786[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, hl_382, hl_383, hl_384, hl_385, \
                         hl_386, kl_787, kl_788, kl_789, kl_790, \
                         kl_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = -2.0 * hl_382[k]
                   + f_0 * kl_787[k];

        t_563[k] = -2.0 * hl_383[k]
                   + f_0 * kl_788[k];

        t_564[k] = -2.0 * hl_384[k]
                   + f_0 * kl_789[k];

        t_565[k] = -2.0 * hl_385[k]
                   + f_0 * kl_790[k];

        t_566[k] = -2.0 * hl_386[k]
                   + f_0 * kl_791[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, hl_387, hl_388, hl_389, hl_390, \
                         hl_391, kl_792, kl_793, kl_794, kl_795, \
                         kl_796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -2.0 * hl_387[k]
                   + f_0 * kl_792[k];

        t_568[k] = -2.0 * hl_388[k]
                   + f_0 * kl_793[k];

        t_569[k] = -2.0 * hl_389[k]
                   + f_0 * kl_794[k];

        t_570[k] = -2.0 * hl_390[k]
                   + f_0 * kl_795[k];

        t_571[k] = -2.0 * hl_391[k]
                   + f_0 * kl_796[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, hl_392, hl_393, hl_394, hl_395, \
                         hl_396, kl_797, kl_798, kl_799, kl_800, \
                         kl_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = -2.0 * hl_392[k]
                   + f_0 * kl_797[k];

        t_573[k] = -2.0 * hl_393[k]
                   + f_0 * kl_798[k];

        t_574[k] = -2.0 * hl_394[k]
                   + f_0 * kl_799[k];

        t_575[k] = -2.0 * hl_395[k]
                   + f_0 * kl_800[k];

        t_576[k] = -2.0 * hl_396[k]
                   + f_0 * kl_801[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, hl_397, hl_398, hl_399, hl_400, \
                         hl_401, kl_802, kl_803, kl_804, kl_805, \
                         kl_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = -2.0 * hl_397[k]
                   + f_0 * kl_802[k];

        t_578[k] = -2.0 * hl_398[k]
                   + f_0 * kl_803[k];

        t_579[k] = -2.0 * hl_399[k]
                   + f_0 * kl_804[k];

        t_580[k] = -2.0 * hl_400[k]
                   + f_0 * kl_805[k];

        t_581[k] = -2.0 * hl_401[k]
                   + f_0 * kl_806[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, hl_402, hl_403, hl_404, hl_405, \
                         hl_406, kl_807, kl_808, kl_809, kl_810, \
                         kl_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = -2.0 * hl_402[k]
                   + f_0 * kl_807[k];

        t_583[k] = -2.0 * hl_403[k]
                   + f_0 * kl_808[k];

        t_584[k] = -2.0 * hl_404[k]
                   + f_0 * kl_809[k];

        t_585[k] = -hl_405[k]
                   + f_0 * kl_810[k];

        t_586[k] = -hl_406[k]
                   + f_0 * kl_811[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, t_591, hl_407, hl_408, hl_409, hl_410, \
                         hl_411, kl_812, kl_813, kl_814, kl_815, \
                         kl_816 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = -hl_407[k]
                   + f_0 * kl_812[k];

        t_588[k] = -hl_408[k]
                   + f_0 * kl_813[k];

        t_589[k] = -hl_409[k]
                   + f_0 * kl_814[k];

        t_590[k] = -hl_410[k]
                   + f_0 * kl_815[k];

        t_591[k] = -hl_411[k]
                   + f_0 * kl_816[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, t_596, hl_412, hl_413, hl_414, hl_415, \
                         hl_416, kl_817, kl_818, kl_819, kl_820, \
                         kl_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = -hl_412[k]
                   + f_0 * kl_817[k];

        t_593[k] = -hl_413[k]
                   + f_0 * kl_818[k];

        t_594[k] = -hl_414[k]
                   + f_0 * kl_819[k];

        t_595[k] = -hl_415[k]
                   + f_0 * kl_820[k];

        t_596[k] = -hl_416[k]
                   + f_0 * kl_821[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, t_600, t_601, hl_417, hl_418, hl_419, hl_420, \
                         hl_421, kl_822, kl_823, kl_824, kl_825, \
                         kl_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = -hl_417[k]
                   + f_0 * kl_822[k];

        t_598[k] = -hl_418[k]
                   + f_0 * kl_823[k];

        t_599[k] = -hl_419[k]
                   + f_0 * kl_824[k];

        t_600[k] = -hl_420[k]
                   + f_0 * kl_825[k];

        t_601[k] = -hl_421[k]
                   + f_0 * kl_826[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, hl_422, hl_423, hl_424, hl_425, \
                         hl_426, kl_827, kl_828, kl_829, kl_830, \
                         kl_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = -hl_422[k]
                   + f_0 * kl_827[k];

        t_603[k] = -hl_423[k]
                   + f_0 * kl_828[k];

        t_604[k] = -hl_424[k]
                   + f_0 * kl_829[k];

        t_605[k] = -hl_425[k]
                   + f_0 * kl_830[k];

        t_606[k] = -hl_426[k]
                   + f_0 * kl_831[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, t_611, hl_427, hl_428, hl_429, hl_430, \
                         hl_431, kl_832, kl_833, kl_834, kl_835, \
                         kl_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = -hl_427[k]
                   + f_0 * kl_832[k];

        t_608[k] = -hl_428[k]
                   + f_0 * kl_833[k];

        t_609[k] = -hl_429[k]
                   + f_0 * kl_834[k];

        t_610[k] = -hl_430[k]
                   + f_0 * kl_835[k];

        t_611[k] = -hl_431[k]
                   + f_0 * kl_836[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, hl_432, hl_433, hl_434, hl_435, \
                         hl_436, kl_837, kl_838, kl_839, kl_840, \
                         kl_841 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = -hl_432[k]
                   + f_0 * kl_837[k];

        t_613[k] = -hl_433[k]
                   + f_0 * kl_838[k];

        t_614[k] = -hl_434[k]
                   + f_0 * kl_839[k];

        t_615[k] = -hl_435[k]
                   + f_0 * kl_840[k];

        t_616[k] = -hl_436[k]
                   + f_0 * kl_841[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, hl_437, hl_438, hl_439, hl_440, \
                         hl_441, kl_842, kl_843, kl_844, kl_845, \
                         kl_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = -hl_437[k]
                   + f_0 * kl_842[k];

        t_618[k] = -hl_438[k]
                   + f_0 * kl_843[k];

        t_619[k] = -hl_439[k]
                   + f_0 * kl_844[k];

        t_620[k] = -hl_440[k]
                   + f_0 * kl_845[k];

        t_621[k] = -hl_441[k]
                   + f_0 * kl_846[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, hl_442, hl_443, hl_444, hl_445, \
                         hl_446, kl_847, kl_848, kl_849, kl_850, \
                         kl_851 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = -hl_442[k]
                   + f_0 * kl_847[k];

        t_623[k] = -hl_443[k]
                   + f_0 * kl_848[k];

        t_624[k] = -hl_444[k]
                   + f_0 * kl_849[k];

        t_625[k] = -hl_445[k]
                   + f_0 * kl_850[k];

        t_626[k] = -hl_446[k]
                   + f_0 * kl_851[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, t_631, t_632, hl_447, hl_448, hl_449, \
                         kl_852, kl_853, kl_854, kl_855, kl_856, \
                         kl_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = -hl_447[k]
                   + f_0 * kl_852[k];

        t_628[k] = -hl_448[k]
                   + f_0 * kl_853[k];

        t_629[k] = -hl_449[k]
                   + f_0 * kl_854[k];

        t_630[k] = f_0 * kl_855[k];

        t_631[k] = f_0 * kl_856[k];

        t_632[k] = f_0 * kl_857[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, t_636, t_637, t_638, t_639, t_640, kl_858, \
                         kl_859, kl_860, kl_861, kl_862, kl_863, kl_864, \
                         kl_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_0 * kl_858[k];

        t_634[k] = f_0 * kl_859[k];

        t_635[k] = f_0 * kl_860[k];

        t_636[k] = f_0 * kl_861[k];

        t_637[k] = f_0 * kl_862[k];

        t_638[k] = f_0 * kl_863[k];

        t_639[k] = f_0 * kl_864[k];

        t_640[k] = f_0 * kl_865[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, t_645, t_646, t_647, t_648, kl_866, \
                         kl_867, kl_868, kl_869, kl_870, kl_871, kl_872, \
                         kl_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_0 * kl_866[k];

        t_642[k] = f_0 * kl_867[k];

        t_643[k] = f_0 * kl_868[k];

        t_644[k] = f_0 * kl_869[k];

        t_645[k] = f_0 * kl_870[k];

        t_646[k] = f_0 * kl_871[k];

        t_647[k] = f_0 * kl_872[k];

        t_648[k] = f_0 * kl_873[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, t_654, t_655, t_656, kl_874, \
                         kl_875, kl_876, kl_877, kl_878, kl_879, kl_880, \
                         kl_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_0 * kl_874[k];

        t_650[k] = f_0 * kl_875[k];

        t_651[k] = f_0 * kl_876[k];

        t_652[k] = f_0 * kl_877[k];

        t_653[k] = f_0 * kl_878[k];

        t_654[k] = f_0 * kl_879[k];

        t_655[k] = f_0 * kl_880[k];

        t_656[k] = f_0 * kl_881[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, t_662, t_663, t_664, kl_882, \
                         kl_883, kl_884, kl_885, kl_886, kl_887, kl_888, \
                         kl_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_0 * kl_882[k];

        t_658[k] = f_0 * kl_883[k];

        t_659[k] = f_0 * kl_884[k];

        t_660[k] = f_0 * kl_885[k];

        t_661[k] = f_0 * kl_886[k];

        t_662[k] = f_0 * kl_887[k];

        t_663[k] = f_0 * kl_888[k];

        t_664[k] = f_0 * kl_889[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, t_670, t_671, t_672, kl_890, \
                         kl_891, kl_892, kl_893, kl_894, kl_895, kl_896, \
                         kl_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_0 * kl_890[k];

        t_666[k] = f_0 * kl_891[k];

        t_667[k] = f_0 * kl_892[k];

        t_668[k] = f_0 * kl_893[k];

        t_669[k] = f_0 * kl_894[k];

        t_670[k] = f_0 * kl_895[k];

        t_671[k] = f_0 * kl_896[k];

        t_672[k] = f_0 * kl_897[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hl_450 = buffer.data(hl + 450);
    const auto *hl_451 = buffer.data(hl + 451);
    const auto *hl_452 = buffer.data(hl + 452);
    const auto *hl_453 = buffer.data(hl + 453);
    const auto *hl_454 = buffer.data(hl + 454);
    const auto *hl_455 = buffer.data(hl + 455);
    const auto *hl_456 = buffer.data(hl + 456);
    const auto *hl_457 = buffer.data(hl + 457);
    const auto *hl_458 = buffer.data(hl + 458);
    const auto *hl_459 = buffer.data(hl + 459);
    const auto *hl_460 = buffer.data(hl + 460);
    const auto *hl_461 = buffer.data(hl + 461);
    const auto *hl_462 = buffer.data(hl + 462);
    const auto *hl_463 = buffer.data(hl + 463);
    const auto *hl_464 = buffer.data(hl + 464);
    const auto *hl_465 = buffer.data(hl + 465);
    const auto *hl_466 = buffer.data(hl + 466);
    const auto *hl_467 = buffer.data(hl + 467);
    const auto *hl_468 = buffer.data(hl + 468);
    const auto *hl_469 = buffer.data(hl + 469);
    const auto *hl_470 = buffer.data(hl + 470);
    const auto *hl_471 = buffer.data(hl + 471);
    const auto *hl_472 = buffer.data(hl + 472);
    const auto *hl_473 = buffer.data(hl + 473);
    const auto *hl_474 = buffer.data(hl + 474);
    const auto *hl_475 = buffer.data(hl + 475);
    const auto *hl_476 = buffer.data(hl + 476);
    const auto *hl_477 = buffer.data(hl + 477);
    const auto *hl_478 = buffer.data(hl + 478);
    const auto *hl_479 = buffer.data(hl + 479);
    const auto *hl_480 = buffer.data(hl + 480);
    const auto *hl_481 = buffer.data(hl + 481);
    const auto *hl_482 = buffer.data(hl + 482);
    const auto *hl_483 = buffer.data(hl + 483);
    const auto *hl_484 = buffer.data(hl + 484);
    const auto *hl_485 = buffer.data(hl + 485);
    const auto *hl_486 = buffer.data(hl + 486);
    const auto *hl_487 = buffer.data(hl + 487);
    const auto *hl_488 = buffer.data(hl + 488);
    const auto *hl_489 = buffer.data(hl + 489);
    const auto *hl_490 = buffer.data(hl + 490);
    const auto *hl_491 = buffer.data(hl + 491);
    const auto *hl_492 = buffer.data(hl + 492);
    const auto *hl_493 = buffer.data(hl + 493);
    const auto *hl_494 = buffer.data(hl + 494);
    const auto *hl_495 = buffer.data(hl + 495);
    const auto *hl_496 = buffer.data(hl + 496);
    const auto *hl_497 = buffer.data(hl + 497);
    const auto *hl_498 = buffer.data(hl + 498);
    const auto *hl_499 = buffer.data(hl + 499);
    const auto *hl_500 = buffer.data(hl + 500);
    const auto *hl_501 = buffer.data(hl + 501);
    const auto *hl_502 = buffer.data(hl + 502);
    const auto *hl_503 = buffer.data(hl + 503);
    const auto *hl_504 = buffer.data(hl + 504);
    const auto *hl_505 = buffer.data(hl + 505);
    const auto *hl_506 = buffer.data(hl + 506);
    const auto *hl_507 = buffer.data(hl + 507);
    const auto *hl_508 = buffer.data(hl + 508);
    const auto *hl_509 = buffer.data(hl + 509);
    const auto *hl_510 = buffer.data(hl + 510);
    const auto *hl_511 = buffer.data(hl + 511);
    const auto *hl_512 = buffer.data(hl + 512);
    const auto *hl_513 = buffer.data(hl + 513);
    const auto *hl_514 = buffer.data(hl + 514);
    const auto *hl_515 = buffer.data(hl + 515);
    const auto *hl_516 = buffer.data(hl + 516);
    const auto *hl_517 = buffer.data(hl + 517);
    const auto *hl_518 = buffer.data(hl + 518);
    const auto *hl_519 = buffer.data(hl + 519);
    const auto *hl_520 = buffer.data(hl + 520);
    const auto *hl_521 = buffer.data(hl + 521);
    const auto *hl_522 = buffer.data(hl + 522);
    const auto *hl_523 = buffer.data(hl + 523);
    const auto *hl_524 = buffer.data(hl + 524);
    const auto *hl_525 = buffer.data(hl + 525);
    const auto *hl_526 = buffer.data(hl + 526);
    const auto *hl_527 = buffer.data(hl + 527);
    const auto *hl_528 = buffer.data(hl + 528);
    const auto *hl_529 = buffer.data(hl + 529);
    const auto *hl_530 = buffer.data(hl + 530);
    const auto *hl_531 = buffer.data(hl + 531);
    const auto *hl_532 = buffer.data(hl + 532);
    const auto *hl_533 = buffer.data(hl + 533);
    const auto *hl_534 = buffer.data(hl + 534);
    const auto *hl_535 = buffer.data(hl + 535);
    const auto *hl_536 = buffer.data(hl + 536);
    const auto *hl_537 = buffer.data(hl + 537);
    const auto *hl_538 = buffer.data(hl + 538);
    const auto *hl_539 = buffer.data(hl + 539);
    const auto *hl_540 = buffer.data(hl + 540);
    const auto *hl_541 = buffer.data(hl + 541);
    const auto *hl_542 = buffer.data(hl + 542);
    const auto *hl_543 = buffer.data(hl + 543);
    const auto *hl_544 = buffer.data(hl + 544);
    const auto *hl_545 = buffer.data(hl + 545);
    const auto *hl_546 = buffer.data(hl + 546);
    const auto *hl_547 = buffer.data(hl + 547);
    const auto *hl_548 = buffer.data(hl + 548);
    const auto *hl_549 = buffer.data(hl + 549);
    const auto *hl_550 = buffer.data(hl + 550);
    const auto *hl_551 = buffer.data(hl + 551);
    const auto *hl_552 = buffer.data(hl + 552);
    const auto *hl_553 = buffer.data(hl + 553);
    const auto *hl_554 = buffer.data(hl + 554);
    const auto *hl_555 = buffer.data(hl + 555);
    const auto *hl_556 = buffer.data(hl + 556);
    const auto *hl_557 = buffer.data(hl + 557);
    const auto *hl_558 = buffer.data(hl + 558);
    const auto *hl_559 = buffer.data(hl + 559);
    const auto *hl_560 = buffer.data(hl + 560);
    const auto *hl_561 = buffer.data(hl + 561);
    const auto *hl_562 = buffer.data(hl + 562);
    const auto *hl_563 = buffer.data(hl + 563);
    const auto *hl_564 = buffer.data(hl + 564);
    const auto *hl_565 = buffer.data(hl + 565);
    const auto *hl_566 = buffer.data(hl + 566);
    const auto *hl_567 = buffer.data(hl + 567);
    const auto *hl_568 = buffer.data(hl + 568);
    const auto *hl_569 = buffer.data(hl + 569);
    const auto *hl_570 = buffer.data(hl + 570);
    const auto *hl_571 = buffer.data(hl + 571);
    const auto *hl_572 = buffer.data(hl + 572);
    const auto *hl_573 = buffer.data(hl + 573);
    const auto *hl_574 = buffer.data(hl + 574);
    const auto *hl_575 = buffer.data(hl + 575);
    const auto *hl_576 = buffer.data(hl + 576);
    const auto *hl_577 = buffer.data(hl + 577);
    const auto *hl_578 = buffer.data(hl + 578);
    const auto *hl_579 = buffer.data(hl + 579);
    const auto *hl_580 = buffer.data(hl + 580);
    const auto *hl_581 = buffer.data(hl + 581);
    const auto *hl_582 = buffer.data(hl + 582);
    const auto *hl_583 = buffer.data(hl + 583);
    const auto *hl_584 = buffer.data(hl + 584);
    const auto *hl_585 = buffer.data(hl + 585);
    const auto *hl_586 = buffer.data(hl + 586);
    const auto *hl_587 = buffer.data(hl + 587);
    const auto *hl_588 = buffer.data(hl + 588);
    const auto *hl_589 = buffer.data(hl + 589);
    const auto *hl_590 = buffer.data(hl + 590);
    const auto *hl_591 = buffer.data(hl + 591);
    const auto *hl_592 = buffer.data(hl + 592);
    const auto *hl_593 = buffer.data(hl + 593);
    const auto *hl_594 = buffer.data(hl + 594);
    const auto *hl_595 = buffer.data(hl + 595);
    const auto *hl_596 = buffer.data(hl + 596);
    const auto *hl_597 = buffer.data(hl + 597);
    const auto *hl_598 = buffer.data(hl + 598);

    const auto *kl_898 = buffer.data(kl + 898);
    const auto *kl_899 = buffer.data(kl + 899);
    const auto *kl_945 = buffer.data(kl + 945);
    const auto *kl_946 = buffer.data(kl + 946);
    const auto *kl_947 = buffer.data(kl + 947);
    const auto *kl_948 = buffer.data(kl + 948);
    const auto *kl_949 = buffer.data(kl + 949);
    const auto *kl_950 = buffer.data(kl + 950);
    const auto *kl_951 = buffer.data(kl + 951);
    const auto *kl_952 = buffer.data(kl + 952);
    const auto *kl_953 = buffer.data(kl + 953);
    const auto *kl_954 = buffer.data(kl + 954);
    const auto *kl_955 = buffer.data(kl + 955);
    const auto *kl_956 = buffer.data(kl + 956);
    const auto *kl_957 = buffer.data(kl + 957);
    const auto *kl_958 = buffer.data(kl + 958);
    const auto *kl_959 = buffer.data(kl + 959);
    const auto *kl_960 = buffer.data(kl + 960);
    const auto *kl_961 = buffer.data(kl + 961);
    const auto *kl_962 = buffer.data(kl + 962);
    const auto *kl_963 = buffer.data(kl + 963);
    const auto *kl_964 = buffer.data(kl + 964);
    const auto *kl_965 = buffer.data(kl + 965);
    const auto *kl_966 = buffer.data(kl + 966);
    const auto *kl_967 = buffer.data(kl + 967);
    const auto *kl_968 = buffer.data(kl + 968);
    const auto *kl_969 = buffer.data(kl + 969);
    const auto *kl_970 = buffer.data(kl + 970);
    const auto *kl_971 = buffer.data(kl + 971);
    const auto *kl_972 = buffer.data(kl + 972);
    const auto *kl_973 = buffer.data(kl + 973);
    const auto *kl_974 = buffer.data(kl + 974);
    const auto *kl_975 = buffer.data(kl + 975);
    const auto *kl_976 = buffer.data(kl + 976);
    const auto *kl_977 = buffer.data(kl + 977);
    const auto *kl_978 = buffer.data(kl + 978);
    const auto *kl_979 = buffer.data(kl + 979);
    const auto *kl_980 = buffer.data(kl + 980);
    const auto *kl_981 = buffer.data(kl + 981);
    const auto *kl_982 = buffer.data(kl + 982);
    const auto *kl_983 = buffer.data(kl + 983);
    const auto *kl_984 = buffer.data(kl + 984);
    const auto *kl_985 = buffer.data(kl + 985);
    const auto *kl_986 = buffer.data(kl + 986);
    const auto *kl_987 = buffer.data(kl + 987);
    const auto *kl_988 = buffer.data(kl + 988);
    const auto *kl_989 = buffer.data(kl + 989);
    const auto *kl_990 = buffer.data(kl + 990);
    const auto *kl_991 = buffer.data(kl + 991);
    const auto *kl_992 = buffer.data(kl + 992);
    const auto *kl_993 = buffer.data(kl + 993);
    const auto *kl_994 = buffer.data(kl + 994);
    const auto *kl_995 = buffer.data(kl + 995);
    const auto *kl_996 = buffer.data(kl + 996);
    const auto *kl_997 = buffer.data(kl + 997);
    const auto *kl_998 = buffer.data(kl + 998);
    const auto *kl_999 = buffer.data(kl + 999);
    const auto *kl_1000 = buffer.data(kl + 1000);
    const auto *kl_1001 = buffer.data(kl + 1001);
    const auto *kl_1002 = buffer.data(kl + 1002);
    const auto *kl_1003 = buffer.data(kl + 1003);
    const auto *kl_1004 = buffer.data(kl + 1004);
    const auto *kl_1005 = buffer.data(kl + 1005);
    const auto *kl_1006 = buffer.data(kl + 1006);
    const auto *kl_1007 = buffer.data(kl + 1007);
    const auto *kl_1008 = buffer.data(kl + 1008);
    const auto *kl_1009 = buffer.data(kl + 1009);
    const auto *kl_1010 = buffer.data(kl + 1010);
    const auto *kl_1011 = buffer.data(kl + 1011);
    const auto *kl_1012 = buffer.data(kl + 1012);
    const auto *kl_1013 = buffer.data(kl + 1013);
    const auto *kl_1014 = buffer.data(kl + 1014);
    const auto *kl_1015 = buffer.data(kl + 1015);
    const auto *kl_1016 = buffer.data(kl + 1016);
    const auto *kl_1017 = buffer.data(kl + 1017);
    const auto *kl_1018 = buffer.data(kl + 1018);
    const auto *kl_1019 = buffer.data(kl + 1019);
    const auto *kl_1020 = buffer.data(kl + 1020);
    const auto *kl_1021 = buffer.data(kl + 1021);
    const auto *kl_1022 = buffer.data(kl + 1022);
    const auto *kl_1023 = buffer.data(kl + 1023);
    const auto *kl_1024 = buffer.data(kl + 1024);
    const auto *kl_1025 = buffer.data(kl + 1025);
    const auto *kl_1026 = buffer.data(kl + 1026);
    const auto *kl_1027 = buffer.data(kl + 1027);
    const auto *kl_1028 = buffer.data(kl + 1028);
    const auto *kl_1029 = buffer.data(kl + 1029);
    const auto *kl_1030 = buffer.data(kl + 1030);
    const auto *kl_1031 = buffer.data(kl + 1031);
    const auto *kl_1032 = buffer.data(kl + 1032);
    const auto *kl_1033 = buffer.data(kl + 1033);
    const auto *kl_1034 = buffer.data(kl + 1034);
    const auto *kl_1035 = buffer.data(kl + 1035);
    const auto *kl_1036 = buffer.data(kl + 1036);
    const auto *kl_1037 = buffer.data(kl + 1037);
    const auto *kl_1038 = buffer.data(kl + 1038);
    const auto *kl_1039 = buffer.data(kl + 1039);
    const auto *kl_1040 = buffer.data(kl + 1040);
    const auto *kl_1041 = buffer.data(kl + 1041);
    const auto *kl_1042 = buffer.data(kl + 1042);
    const auto *kl_1043 = buffer.data(kl + 1043);
    const auto *kl_1044 = buffer.data(kl + 1044);
    const auto *kl_1045 = buffer.data(kl + 1045);
    const auto *kl_1046 = buffer.data(kl + 1046);
    const auto *kl_1047 = buffer.data(kl + 1047);
    const auto *kl_1048 = buffer.data(kl + 1048);
    const auto *kl_1049 = buffer.data(kl + 1049);
    const auto *kl_1050 = buffer.data(kl + 1050);
    const auto *kl_1051 = buffer.data(kl + 1051);
    const auto *kl_1052 = buffer.data(kl + 1052);
    const auto *kl_1053 = buffer.data(kl + 1053);
    const auto *kl_1054 = buffer.data(kl + 1054);
    const auto *kl_1055 = buffer.data(kl + 1055);
    const auto *kl_1056 = buffer.data(kl + 1056);
    const auto *kl_1057 = buffer.data(kl + 1057);
    const auto *kl_1058 = buffer.data(kl + 1058);
    const auto *kl_1059 = buffer.data(kl + 1059);
    const auto *kl_1060 = buffer.data(kl + 1060);
    const auto *kl_1061 = buffer.data(kl + 1061);
    const auto *kl_1062 = buffer.data(kl + 1062);
    const auto *kl_1063 = buffer.data(kl + 1063);
    const auto *kl_1064 = buffer.data(kl + 1064);
    const auto *kl_1065 = buffer.data(kl + 1065);
    const auto *kl_1066 = buffer.data(kl + 1066);
    const auto *kl_1067 = buffer.data(kl + 1067);
    const auto *kl_1068 = buffer.data(kl + 1068);
    const auto *kl_1069 = buffer.data(kl + 1069);
    const auto *kl_1070 = buffer.data(kl + 1070);
    const auto *kl_1071 = buffer.data(kl + 1071);
    const auto *kl_1072 = buffer.data(kl + 1072);
    const auto *kl_1073 = buffer.data(kl + 1073);
    const auto *kl_1074 = buffer.data(kl + 1074);
    const auto *kl_1075 = buffer.data(kl + 1075);
    const auto *kl_1076 = buffer.data(kl + 1076);
    const auto *kl_1077 = buffer.data(kl + 1077);
    const auto *kl_1078 = buffer.data(kl + 1078);
    const auto *kl_1079 = buffer.data(kl + 1079);
    const auto *kl_1080 = buffer.data(kl + 1080);
    const auto *kl_1081 = buffer.data(kl + 1081);
    const auto *kl_1082 = buffer.data(kl + 1082);
    const auto *kl_1083 = buffer.data(kl + 1083);
    const auto *kl_1084 = buffer.data(kl + 1084);
    const auto *kl_1085 = buffer.data(kl + 1085);
    const auto *kl_1086 = buffer.data(kl + 1086);
    const auto *kl_1087 = buffer.data(kl + 1087);
    const auto *kl_1088 = buffer.data(kl + 1088);
    const auto *kl_1089 = buffer.data(kl + 1089);
    const auto *kl_1090 = buffer.data(kl + 1090);
    const auto *kl_1091 = buffer.data(kl + 1091);
    const auto *kl_1092 = buffer.data(kl + 1092);
    const auto *kl_1093 = buffer.data(kl + 1093);

#pragma omp simd aligned(t_673, t_674, t_675, t_676, t_677, t_678, hl_450, hl_451, hl_452, \
                         hl_453, kl_898, kl_899, kl_945, kl_946, kl_947, \
                         kl_948 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_0 * kl_898[k];

        t_674[k] = f_0 * kl_899[k];

        t_675[k] = -5.0 * hl_450[k]
                   + f_0 * kl_945[k];

        t_676[k] = -5.0 * hl_451[k]
                   + f_0 * kl_946[k];

        t_677[k] = -5.0 * hl_452[k]
                   + f_0 * kl_947[k];

        t_678[k] = -5.0 * hl_453[k]
                   + f_0 * kl_948[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, t_682, t_683, hl_454, hl_455, hl_456, hl_457, \
                         hl_458, kl_949, kl_950, kl_951, kl_952, \
                         kl_953 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = -5.0 * hl_454[k]
                   + f_0 * kl_949[k];

        t_680[k] = -5.0 * hl_455[k]
                   + f_0 * kl_950[k];

        t_681[k] = -5.0 * hl_456[k]
                   + f_0 * kl_951[k];

        t_682[k] = -5.0 * hl_457[k]
                   + f_0 * kl_952[k];

        t_683[k] = -5.0 * hl_458[k]
                   + f_0 * kl_953[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, t_688, hl_459, hl_460, hl_461, hl_462, \
                         hl_463, kl_954, kl_955, kl_956, kl_957, \
                         kl_958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = -5.0 * hl_459[k]
                   + f_0 * kl_954[k];

        t_685[k] = -5.0 * hl_460[k]
                   + f_0 * kl_955[k];

        t_686[k] = -5.0 * hl_461[k]
                   + f_0 * kl_956[k];

        t_687[k] = -5.0 * hl_462[k]
                   + f_0 * kl_957[k];

        t_688[k] = -5.0 * hl_463[k]
                   + f_0 * kl_958[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, t_693, hl_464, hl_465, hl_466, hl_467, \
                         hl_468, kl_959, kl_960, kl_961, kl_962, \
                         kl_963 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = -5.0 * hl_464[k]
                   + f_0 * kl_959[k];

        t_690[k] = -5.0 * hl_465[k]
                   + f_0 * kl_960[k];

        t_691[k] = -5.0 * hl_466[k]
                   + f_0 * kl_961[k];

        t_692[k] = -5.0 * hl_467[k]
                   + f_0 * kl_962[k];

        t_693[k] = -5.0 * hl_468[k]
                   + f_0 * kl_963[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, t_698, hl_469, hl_470, hl_471, hl_472, \
                         hl_473, kl_964, kl_965, kl_966, kl_967, \
                         kl_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = -5.0 * hl_469[k]
                   + f_0 * kl_964[k];

        t_695[k] = -5.0 * hl_470[k]
                   + f_0 * kl_965[k];

        t_696[k] = -5.0 * hl_471[k]
                   + f_0 * kl_966[k];

        t_697[k] = -5.0 * hl_472[k]
                   + f_0 * kl_967[k];

        t_698[k] = -5.0 * hl_473[k]
                   + f_0 * kl_968[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, t_703, hl_474, hl_475, hl_476, hl_477, \
                         hl_478, kl_969, kl_970, kl_971, kl_972, \
                         kl_973 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = -5.0 * hl_474[k]
                   + f_0 * kl_969[k];

        t_700[k] = -5.0 * hl_475[k]
                   + f_0 * kl_970[k];

        t_701[k] = -5.0 * hl_476[k]
                   + f_0 * kl_971[k];

        t_702[k] = -5.0 * hl_477[k]
                   + f_0 * kl_972[k];

        t_703[k] = -5.0 * hl_478[k]
                   + f_0 * kl_973[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, hl_479, hl_480, hl_481, hl_482, \
                         hl_483, kl_974, kl_975, kl_976, kl_977, \
                         kl_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = -5.0 * hl_479[k]
                   + f_0 * kl_974[k];

        t_705[k] = -5.0 * hl_480[k]
                   + f_0 * kl_975[k];

        t_706[k] = -5.0 * hl_481[k]
                   + f_0 * kl_976[k];

        t_707[k] = -5.0 * hl_482[k]
                   + f_0 * kl_977[k];

        t_708[k] = -5.0 * hl_483[k]
                   + f_0 * kl_978[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, t_712, t_713, hl_484, hl_485, hl_486, hl_487, \
                         hl_488, kl_979, kl_980, kl_981, kl_982, \
                         kl_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = -5.0 * hl_484[k]
                   + f_0 * kl_979[k];

        t_710[k] = -5.0 * hl_485[k]
                   + f_0 * kl_980[k];

        t_711[k] = -5.0 * hl_486[k]
                   + f_0 * kl_981[k];

        t_712[k] = -5.0 * hl_487[k]
                   + f_0 * kl_982[k];

        t_713[k] = -5.0 * hl_488[k]
                   + f_0 * kl_983[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, t_718, hl_489, hl_490, hl_491, hl_492, \
                         hl_493, kl_984, kl_985, kl_986, kl_987, \
                         kl_988 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = -5.0 * hl_489[k]
                   + f_0 * kl_984[k];

        t_715[k] = -5.0 * hl_490[k]
                   + f_0 * kl_985[k];

        t_716[k] = -5.0 * hl_491[k]
                   + f_0 * kl_986[k];

        t_717[k] = -5.0 * hl_492[k]
                   + f_0 * kl_987[k];

        t_718[k] = -5.0 * hl_493[k]
                   + f_0 * kl_988[k];
    }

#pragma omp simd aligned(t_719, t_720, t_721, t_722, t_723, hl_494, hl_495, hl_496, hl_497, \
                         hl_498, kl_989, kl_990, kl_991, kl_992, \
                         kl_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_719[k] = -5.0 * hl_494[k]
                   + f_0 * kl_989[k];

        t_720[k] = -4.0 * hl_495[k]
                   + f_0 * kl_990[k];

        t_721[k] = -4.0 * hl_496[k]
                   + f_0 * kl_991[k];

        t_722[k] = -4.0 * hl_497[k]
                   + f_0 * kl_992[k];

        t_723[k] = -4.0 * hl_498[k]
                   + f_0 * kl_993[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, hl_499, hl_500, hl_501, hl_502, \
                         hl_503, kl_994, kl_995, kl_996, kl_997, \
                         kl_998 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = -4.0 * hl_499[k]
                   + f_0 * kl_994[k];

        t_725[k] = -4.0 * hl_500[k]
                   + f_0 * kl_995[k];

        t_726[k] = -4.0 * hl_501[k]
                   + f_0 * kl_996[k];

        t_727[k] = -4.0 * hl_502[k]
                   + f_0 * kl_997[k];

        t_728[k] = -4.0 * hl_503[k]
                   + f_0 * kl_998[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, t_733, hl_504, hl_505, hl_506, hl_507, \
                         hl_508, kl_999, kl_1000, kl_1001, kl_1002, \
                         kl_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = -4.0 * hl_504[k]
                   + f_0 * kl_999[k];

        t_730[k] = -4.0 * hl_505[k]
                   + f_0 * kl_1000[k];

        t_731[k] = -4.0 * hl_506[k]
                   + f_0 * kl_1001[k];

        t_732[k] = -4.0 * hl_507[k]
                   + f_0 * kl_1002[k];

        t_733[k] = -4.0 * hl_508[k]
                   + f_0 * kl_1003[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, t_738, hl_509, hl_510, hl_511, hl_512, \
                         hl_513, kl_1004, kl_1005, kl_1006, kl_1007, \
                         kl_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = -4.0 * hl_509[k]
                   + f_0 * kl_1004[k];

        t_735[k] = -4.0 * hl_510[k]
                   + f_0 * kl_1005[k];

        t_736[k] = -4.0 * hl_511[k]
                   + f_0 * kl_1006[k];

        t_737[k] = -4.0 * hl_512[k]
                   + f_0 * kl_1007[k];

        t_738[k] = -4.0 * hl_513[k]
                   + f_0 * kl_1008[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, t_743, hl_514, hl_515, hl_516, hl_517, \
                         hl_518, kl_1009, kl_1010, kl_1011, kl_1012, \
                         kl_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = -4.0 * hl_514[k]
                   + f_0 * kl_1009[k];

        t_740[k] = -4.0 * hl_515[k]
                   + f_0 * kl_1010[k];

        t_741[k] = -4.0 * hl_516[k]
                   + f_0 * kl_1011[k];

        t_742[k] = -4.0 * hl_517[k]
                   + f_0 * kl_1012[k];

        t_743[k] = -4.0 * hl_518[k]
                   + f_0 * kl_1013[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, t_748, hl_519, hl_520, hl_521, hl_522, \
                         hl_523, kl_1014, kl_1015, kl_1016, kl_1017, \
                         kl_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = -4.0 * hl_519[k]
                   + f_0 * kl_1014[k];

        t_745[k] = -4.0 * hl_520[k]
                   + f_0 * kl_1015[k];

        t_746[k] = -4.0 * hl_521[k]
                   + f_0 * kl_1016[k];

        t_747[k] = -4.0 * hl_522[k]
                   + f_0 * kl_1017[k];

        t_748[k] = -4.0 * hl_523[k]
                   + f_0 * kl_1018[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, t_752, t_753, hl_524, hl_525, hl_526, hl_527, \
                         hl_528, kl_1019, kl_1020, kl_1021, kl_1022, \
                         kl_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = -4.0 * hl_524[k]
                   + f_0 * kl_1019[k];

        t_750[k] = -4.0 * hl_525[k]
                   + f_0 * kl_1020[k];

        t_751[k] = -4.0 * hl_526[k]
                   + f_0 * kl_1021[k];

        t_752[k] = -4.0 * hl_527[k]
                   + f_0 * kl_1022[k];

        t_753[k] = -4.0 * hl_528[k]
                   + f_0 * kl_1023[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, t_758, hl_529, hl_530, hl_531, hl_532, \
                         hl_533, kl_1024, kl_1025, kl_1026, kl_1027, \
                         kl_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = -4.0 * hl_529[k]
                   + f_0 * kl_1024[k];

        t_755[k] = -4.0 * hl_530[k]
                   + f_0 * kl_1025[k];

        t_756[k] = -4.0 * hl_531[k]
                   + f_0 * kl_1026[k];

        t_757[k] = -4.0 * hl_532[k]
                   + f_0 * kl_1027[k];

        t_758[k] = -4.0 * hl_533[k]
                   + f_0 * kl_1028[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, t_762, t_763, hl_534, hl_535, hl_536, hl_537, \
                         hl_538, kl_1029, kl_1030, kl_1031, kl_1032, \
                         kl_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = -4.0 * hl_534[k]
                   + f_0 * kl_1029[k];

        t_760[k] = -4.0 * hl_535[k]
                   + f_0 * kl_1030[k];

        t_761[k] = -4.0 * hl_536[k]
                   + f_0 * kl_1031[k];

        t_762[k] = -4.0 * hl_537[k]
                   + f_0 * kl_1032[k];

        t_763[k] = -4.0 * hl_538[k]
                   + f_0 * kl_1033[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, t_768, hl_539, hl_540, hl_541, hl_542, \
                         hl_543, kl_1034, kl_1035, kl_1036, kl_1037, \
                         kl_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = -4.0 * hl_539[k]
                   + f_0 * kl_1034[k];

        t_765[k] = -3.0 * hl_540[k]
                   + f_0 * kl_1035[k];

        t_766[k] = -3.0 * hl_541[k]
                   + f_0 * kl_1036[k];

        t_767[k] = -3.0 * hl_542[k]
                   + f_0 * kl_1037[k];

        t_768[k] = -3.0 * hl_543[k]
                   + f_0 * kl_1038[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, t_772, t_773, hl_544, hl_545, hl_546, hl_547, \
                         hl_548, kl_1039, kl_1040, kl_1041, kl_1042, \
                         kl_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = -3.0 * hl_544[k]
                   + f_0 * kl_1039[k];

        t_770[k] = -3.0 * hl_545[k]
                   + f_0 * kl_1040[k];

        t_771[k] = -3.0 * hl_546[k]
                   + f_0 * kl_1041[k];

        t_772[k] = -3.0 * hl_547[k]
                   + f_0 * kl_1042[k];

        t_773[k] = -3.0 * hl_548[k]
                   + f_0 * kl_1043[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, t_778, hl_549, hl_550, hl_551, hl_552, \
                         hl_553, kl_1044, kl_1045, kl_1046, kl_1047, \
                         kl_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = -3.0 * hl_549[k]
                   + f_0 * kl_1044[k];

        t_775[k] = -3.0 * hl_550[k]
                   + f_0 * kl_1045[k];

        t_776[k] = -3.0 * hl_551[k]
                   + f_0 * kl_1046[k];

        t_777[k] = -3.0 * hl_552[k]
                   + f_0 * kl_1047[k];

        t_778[k] = -3.0 * hl_553[k]
                   + f_0 * kl_1048[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, t_782, t_783, hl_554, hl_555, hl_556, hl_557, \
                         hl_558, kl_1049, kl_1050, kl_1051, kl_1052, \
                         kl_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = -3.0 * hl_554[k]
                   + f_0 * kl_1049[k];

        t_780[k] = -3.0 * hl_555[k]
                   + f_0 * kl_1050[k];

        t_781[k] = -3.0 * hl_556[k]
                   + f_0 * kl_1051[k];

        t_782[k] = -3.0 * hl_557[k]
                   + f_0 * kl_1052[k];

        t_783[k] = -3.0 * hl_558[k]
                   + f_0 * kl_1053[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, t_788, hl_559, hl_560, hl_561, hl_562, \
                         hl_563, kl_1054, kl_1055, kl_1056, kl_1057, \
                         kl_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = -3.0 * hl_559[k]
                   + f_0 * kl_1054[k];

        t_785[k] = -3.0 * hl_560[k]
                   + f_0 * kl_1055[k];

        t_786[k] = -3.0 * hl_561[k]
                   + f_0 * kl_1056[k];

        t_787[k] = -3.0 * hl_562[k]
                   + f_0 * kl_1057[k];

        t_788[k] = -3.0 * hl_563[k]
                   + f_0 * kl_1058[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, t_792, t_793, hl_564, hl_565, hl_566, hl_567, \
                         hl_568, kl_1059, kl_1060, kl_1061, kl_1062, \
                         kl_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = -3.0 * hl_564[k]
                   + f_0 * kl_1059[k];

        t_790[k] = -3.0 * hl_565[k]
                   + f_0 * kl_1060[k];

        t_791[k] = -3.0 * hl_566[k]
                   + f_0 * kl_1061[k];

        t_792[k] = -3.0 * hl_567[k]
                   + f_0 * kl_1062[k];

        t_793[k] = -3.0 * hl_568[k]
                   + f_0 * kl_1063[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, t_798, hl_569, hl_570, hl_571, hl_572, \
                         hl_573, kl_1064, kl_1065, kl_1066, kl_1067, \
                         kl_1068 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = -3.0 * hl_569[k]
                   + f_0 * kl_1064[k];

        t_795[k] = -3.0 * hl_570[k]
                   + f_0 * kl_1065[k];

        t_796[k] = -3.0 * hl_571[k]
                   + f_0 * kl_1066[k];

        t_797[k] = -3.0 * hl_572[k]
                   + f_0 * kl_1067[k];

        t_798[k] = -3.0 * hl_573[k]
                   + f_0 * kl_1068[k];
    }

#pragma omp simd aligned(t_799, t_800, t_801, t_802, t_803, hl_574, hl_575, hl_576, hl_577, \
                         hl_578, kl_1069, kl_1070, kl_1071, kl_1072, \
                         kl_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = -3.0 * hl_574[k]
                   + f_0 * kl_1069[k];

        t_800[k] = -3.0 * hl_575[k]
                   + f_0 * kl_1070[k];

        t_801[k] = -3.0 * hl_576[k]
                   + f_0 * kl_1071[k];

        t_802[k] = -3.0 * hl_577[k]
                   + f_0 * kl_1072[k];

        t_803[k] = -3.0 * hl_578[k]
                   + f_0 * kl_1073[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, t_807, t_808, hl_579, hl_580, hl_581, hl_582, \
                         hl_583, kl_1074, kl_1075, kl_1076, kl_1077, \
                         kl_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = -3.0 * hl_579[k]
                   + f_0 * kl_1074[k];

        t_805[k] = -3.0 * hl_580[k]
                   + f_0 * kl_1075[k];

        t_806[k] = -3.0 * hl_581[k]
                   + f_0 * kl_1076[k];

        t_807[k] = -3.0 * hl_582[k]
                   + f_0 * kl_1077[k];

        t_808[k] = -3.0 * hl_583[k]
                   + f_0 * kl_1078[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, t_812, t_813, hl_584, hl_585, hl_586, hl_587, \
                         hl_588, kl_1079, kl_1080, kl_1081, kl_1082, \
                         kl_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = -3.0 * hl_584[k]
                   + f_0 * kl_1079[k];

        t_810[k] = -2.0 * hl_585[k]
                   + f_0 * kl_1080[k];

        t_811[k] = -2.0 * hl_586[k]
                   + f_0 * kl_1081[k];

        t_812[k] = -2.0 * hl_587[k]
                   + f_0 * kl_1082[k];

        t_813[k] = -2.0 * hl_588[k]
                   + f_0 * kl_1083[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, t_817, t_818, hl_589, hl_590, hl_591, hl_592, \
                         hl_593, kl_1084, kl_1085, kl_1086, kl_1087, \
                         kl_1088 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = -2.0 * hl_589[k]
                   + f_0 * kl_1084[k];

        t_815[k] = -2.0 * hl_590[k]
                   + f_0 * kl_1085[k];

        t_816[k] = -2.0 * hl_591[k]
                   + f_0 * kl_1086[k];

        t_817[k] = -2.0 * hl_592[k]
                   + f_0 * kl_1087[k];

        t_818[k] = -2.0 * hl_593[k]
                   + f_0 * kl_1088[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, t_822, t_823, hl_594, hl_595, hl_596, hl_597, \
                         hl_598, kl_1089, kl_1090, kl_1091, kl_1092, \
                         kl_1093 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = -2.0 * hl_594[k]
                   + f_0 * kl_1089[k];

        t_820[k] = -2.0 * hl_595[k]
                   + f_0 * kl_1090[k];

        t_821[k] = -2.0 * hl_596[k]
                   + f_0 * kl_1091[k];

        t_822[k] = -2.0 * hl_597[k]
                   + f_0 * kl_1092[k];

        t_823[k] = -2.0 * hl_598[k]
                   + f_0 * kl_1093[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_1_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hl_599 = buffer.data(hl + 599);
    const auto *hl_600 = buffer.data(hl + 600);
    const auto *hl_601 = buffer.data(hl + 601);
    const auto *hl_602 = buffer.data(hl + 602);
    const auto *hl_603 = buffer.data(hl + 603);
    const auto *hl_604 = buffer.data(hl + 604);
    const auto *hl_605 = buffer.data(hl + 605);
    const auto *hl_606 = buffer.data(hl + 606);
    const auto *hl_607 = buffer.data(hl + 607);
    const auto *hl_608 = buffer.data(hl + 608);
    const auto *hl_609 = buffer.data(hl + 609);
    const auto *hl_610 = buffer.data(hl + 610);
    const auto *hl_611 = buffer.data(hl + 611);
    const auto *hl_612 = buffer.data(hl + 612);
    const auto *hl_613 = buffer.data(hl + 613);
    const auto *hl_614 = buffer.data(hl + 614);
    const auto *hl_615 = buffer.data(hl + 615);
    const auto *hl_616 = buffer.data(hl + 616);
    const auto *hl_617 = buffer.data(hl + 617);
    const auto *hl_618 = buffer.data(hl + 618);
    const auto *hl_619 = buffer.data(hl + 619);
    const auto *hl_620 = buffer.data(hl + 620);
    const auto *hl_621 = buffer.data(hl + 621);
    const auto *hl_622 = buffer.data(hl + 622);
    const auto *hl_623 = buffer.data(hl + 623);
    const auto *hl_624 = buffer.data(hl + 624);
    const auto *hl_625 = buffer.data(hl + 625);
    const auto *hl_626 = buffer.data(hl + 626);
    const auto *hl_627 = buffer.data(hl + 627);
    const auto *hl_628 = buffer.data(hl + 628);
    const auto *hl_629 = buffer.data(hl + 629);
    const auto *hl_630 = buffer.data(hl + 630);
    const auto *hl_631 = buffer.data(hl + 631);
    const auto *hl_632 = buffer.data(hl + 632);
    const auto *hl_633 = buffer.data(hl + 633);
    const auto *hl_634 = buffer.data(hl + 634);
    const auto *hl_635 = buffer.data(hl + 635);
    const auto *hl_636 = buffer.data(hl + 636);
    const auto *hl_637 = buffer.data(hl + 637);
    const auto *hl_638 = buffer.data(hl + 638);
    const auto *hl_639 = buffer.data(hl + 639);
    const auto *hl_640 = buffer.data(hl + 640);
    const auto *hl_641 = buffer.data(hl + 641);
    const auto *hl_642 = buffer.data(hl + 642);
    const auto *hl_643 = buffer.data(hl + 643);
    const auto *hl_644 = buffer.data(hl + 644);
    const auto *hl_645 = buffer.data(hl + 645);
    const auto *hl_646 = buffer.data(hl + 646);
    const auto *hl_647 = buffer.data(hl + 647);
    const auto *hl_648 = buffer.data(hl + 648);
    const auto *hl_649 = buffer.data(hl + 649);
    const auto *hl_650 = buffer.data(hl + 650);
    const auto *hl_651 = buffer.data(hl + 651);
    const auto *hl_652 = buffer.data(hl + 652);
    const auto *hl_653 = buffer.data(hl + 653);
    const auto *hl_654 = buffer.data(hl + 654);
    const auto *hl_655 = buffer.data(hl + 655);
    const auto *hl_656 = buffer.data(hl + 656);
    const auto *hl_657 = buffer.data(hl + 657);
    const auto *hl_658 = buffer.data(hl + 658);
    const auto *hl_659 = buffer.data(hl + 659);
    const auto *hl_660 = buffer.data(hl + 660);
    const auto *hl_661 = buffer.data(hl + 661);
    const auto *hl_662 = buffer.data(hl + 662);
    const auto *hl_663 = buffer.data(hl + 663);
    const auto *hl_664 = buffer.data(hl + 664);
    const auto *hl_665 = buffer.data(hl + 665);
    const auto *hl_666 = buffer.data(hl + 666);
    const auto *hl_667 = buffer.data(hl + 667);
    const auto *hl_668 = buffer.data(hl + 668);
    const auto *hl_669 = buffer.data(hl + 669);
    const auto *hl_670 = buffer.data(hl + 670);
    const auto *hl_671 = buffer.data(hl + 671);
    const auto *hl_672 = buffer.data(hl + 672);
    const auto *hl_673 = buffer.data(hl + 673);
    const auto *hl_674 = buffer.data(hl + 674);
    const auto *hl_675 = buffer.data(hl + 675);
    const auto *hl_676 = buffer.data(hl + 676);
    const auto *hl_677 = buffer.data(hl + 677);
    const auto *hl_678 = buffer.data(hl + 678);
    const auto *hl_679 = buffer.data(hl + 679);
    const auto *hl_680 = buffer.data(hl + 680);
    const auto *hl_681 = buffer.data(hl + 681);
    const auto *hl_682 = buffer.data(hl + 682);
    const auto *hl_683 = buffer.data(hl + 683);
    const auto *hl_684 = buffer.data(hl + 684);
    const auto *hl_685 = buffer.data(hl + 685);
    const auto *hl_686 = buffer.data(hl + 686);
    const auto *hl_687 = buffer.data(hl + 687);
    const auto *hl_688 = buffer.data(hl + 688);
    const auto *hl_689 = buffer.data(hl + 689);
    const auto *hl_690 = buffer.data(hl + 690);
    const auto *hl_691 = buffer.data(hl + 691);
    const auto *hl_692 = buffer.data(hl + 692);
    const auto *hl_693 = buffer.data(hl + 693);
    const auto *hl_694 = buffer.data(hl + 694);
    const auto *hl_695 = buffer.data(hl + 695);
    const auto *hl_696 = buffer.data(hl + 696);
    const auto *hl_697 = buffer.data(hl + 697);
    const auto *hl_698 = buffer.data(hl + 698);
    const auto *hl_699 = buffer.data(hl + 699);
    const auto *hl_700 = buffer.data(hl + 700);
    const auto *hl_701 = buffer.data(hl + 701);
    const auto *hl_702 = buffer.data(hl + 702);
    const auto *hl_703 = buffer.data(hl + 703);
    const auto *hl_704 = buffer.data(hl + 704);
    const auto *hl_705 = buffer.data(hl + 705);
    const auto *hl_706 = buffer.data(hl + 706);
    const auto *hl_707 = buffer.data(hl + 707);
    const auto *hl_708 = buffer.data(hl + 708);
    const auto *hl_709 = buffer.data(hl + 709);
    const auto *hl_710 = buffer.data(hl + 710);
    const auto *hl_711 = buffer.data(hl + 711);
    const auto *hl_712 = buffer.data(hl + 712);
    const auto *hl_713 = buffer.data(hl + 713);
    const auto *hl_714 = buffer.data(hl + 714);
    const auto *hl_715 = buffer.data(hl + 715);
    const auto *hl_716 = buffer.data(hl + 716);
    const auto *hl_717 = buffer.data(hl + 717);
    const auto *hl_718 = buffer.data(hl + 718);
    const auto *hl_719 = buffer.data(hl + 719);

    const auto *kl_1094 = buffer.data(kl + 1094);
    const auto *kl_1095 = buffer.data(kl + 1095);
    const auto *kl_1096 = buffer.data(kl + 1096);
    const auto *kl_1097 = buffer.data(kl + 1097);
    const auto *kl_1098 = buffer.data(kl + 1098);
    const auto *kl_1099 = buffer.data(kl + 1099);
    const auto *kl_1100 = buffer.data(kl + 1100);
    const auto *kl_1101 = buffer.data(kl + 1101);
    const auto *kl_1102 = buffer.data(kl + 1102);
    const auto *kl_1103 = buffer.data(kl + 1103);
    const auto *kl_1104 = buffer.data(kl + 1104);
    const auto *kl_1105 = buffer.data(kl + 1105);
    const auto *kl_1106 = buffer.data(kl + 1106);
    const auto *kl_1107 = buffer.data(kl + 1107);
    const auto *kl_1108 = buffer.data(kl + 1108);
    const auto *kl_1109 = buffer.data(kl + 1109);
    const auto *kl_1110 = buffer.data(kl + 1110);
    const auto *kl_1111 = buffer.data(kl + 1111);
    const auto *kl_1112 = buffer.data(kl + 1112);
    const auto *kl_1113 = buffer.data(kl + 1113);
    const auto *kl_1114 = buffer.data(kl + 1114);
    const auto *kl_1115 = buffer.data(kl + 1115);
    const auto *kl_1116 = buffer.data(kl + 1116);
    const auto *kl_1117 = buffer.data(kl + 1117);
    const auto *kl_1118 = buffer.data(kl + 1118);
    const auto *kl_1119 = buffer.data(kl + 1119);
    const auto *kl_1120 = buffer.data(kl + 1120);
    const auto *kl_1121 = buffer.data(kl + 1121);
    const auto *kl_1122 = buffer.data(kl + 1122);
    const auto *kl_1123 = buffer.data(kl + 1123);
    const auto *kl_1124 = buffer.data(kl + 1124);
    const auto *kl_1125 = buffer.data(kl + 1125);
    const auto *kl_1126 = buffer.data(kl + 1126);
    const auto *kl_1127 = buffer.data(kl + 1127);
    const auto *kl_1128 = buffer.data(kl + 1128);
    const auto *kl_1129 = buffer.data(kl + 1129);
    const auto *kl_1130 = buffer.data(kl + 1130);
    const auto *kl_1131 = buffer.data(kl + 1131);
    const auto *kl_1132 = buffer.data(kl + 1132);
    const auto *kl_1133 = buffer.data(kl + 1133);
    const auto *kl_1134 = buffer.data(kl + 1134);
    const auto *kl_1135 = buffer.data(kl + 1135);
    const auto *kl_1136 = buffer.data(kl + 1136);
    const auto *kl_1137 = buffer.data(kl + 1137);
    const auto *kl_1138 = buffer.data(kl + 1138);
    const auto *kl_1139 = buffer.data(kl + 1139);
    const auto *kl_1140 = buffer.data(kl + 1140);
    const auto *kl_1141 = buffer.data(kl + 1141);
    const auto *kl_1142 = buffer.data(kl + 1142);
    const auto *kl_1143 = buffer.data(kl + 1143);
    const auto *kl_1144 = buffer.data(kl + 1144);
    const auto *kl_1145 = buffer.data(kl + 1145);
    const auto *kl_1146 = buffer.data(kl + 1146);
    const auto *kl_1147 = buffer.data(kl + 1147);
    const auto *kl_1148 = buffer.data(kl + 1148);
    const auto *kl_1149 = buffer.data(kl + 1149);
    const auto *kl_1150 = buffer.data(kl + 1150);
    const auto *kl_1151 = buffer.data(kl + 1151);
    const auto *kl_1152 = buffer.data(kl + 1152);
    const auto *kl_1153 = buffer.data(kl + 1153);
    const auto *kl_1154 = buffer.data(kl + 1154);
    const auto *kl_1155 = buffer.data(kl + 1155);
    const auto *kl_1156 = buffer.data(kl + 1156);
    const auto *kl_1157 = buffer.data(kl + 1157);
    const auto *kl_1158 = buffer.data(kl + 1158);
    const auto *kl_1159 = buffer.data(kl + 1159);
    const auto *kl_1160 = buffer.data(kl + 1160);
    const auto *kl_1161 = buffer.data(kl + 1161);
    const auto *kl_1162 = buffer.data(kl + 1162);
    const auto *kl_1163 = buffer.data(kl + 1163);
    const auto *kl_1164 = buffer.data(kl + 1164);
    const auto *kl_1165 = buffer.data(kl + 1165);
    const auto *kl_1166 = buffer.data(kl + 1166);
    const auto *kl_1167 = buffer.data(kl + 1167);
    const auto *kl_1168 = buffer.data(kl + 1168);
    const auto *kl_1169 = buffer.data(kl + 1169);
    const auto *kl_1170 = buffer.data(kl + 1170);
    const auto *kl_1171 = buffer.data(kl + 1171);
    const auto *kl_1172 = buffer.data(kl + 1172);
    const auto *kl_1173 = buffer.data(kl + 1173);
    const auto *kl_1174 = buffer.data(kl + 1174);
    const auto *kl_1175 = buffer.data(kl + 1175);
    const auto *kl_1176 = buffer.data(kl + 1176);
    const auto *kl_1177 = buffer.data(kl + 1177);
    const auto *kl_1178 = buffer.data(kl + 1178);
    const auto *kl_1179 = buffer.data(kl + 1179);
    const auto *kl_1180 = buffer.data(kl + 1180);
    const auto *kl_1181 = buffer.data(kl + 1181);
    const auto *kl_1182 = buffer.data(kl + 1182);
    const auto *kl_1183 = buffer.data(kl + 1183);
    const auto *kl_1184 = buffer.data(kl + 1184);
    const auto *kl_1185 = buffer.data(kl + 1185);
    const auto *kl_1186 = buffer.data(kl + 1186);
    const auto *kl_1187 = buffer.data(kl + 1187);
    const auto *kl_1188 = buffer.data(kl + 1188);
    const auto *kl_1189 = buffer.data(kl + 1189);
    const auto *kl_1190 = buffer.data(kl + 1190);
    const auto *kl_1191 = buffer.data(kl + 1191);
    const auto *kl_1192 = buffer.data(kl + 1192);
    const auto *kl_1193 = buffer.data(kl + 1193);
    const auto *kl_1194 = buffer.data(kl + 1194);
    const auto *kl_1195 = buffer.data(kl + 1195);
    const auto *kl_1196 = buffer.data(kl + 1196);
    const auto *kl_1197 = buffer.data(kl + 1197);
    const auto *kl_1198 = buffer.data(kl + 1198);
    const auto *kl_1199 = buffer.data(kl + 1199);
    const auto *kl_1200 = buffer.data(kl + 1200);
    const auto *kl_1201 = buffer.data(kl + 1201);
    const auto *kl_1202 = buffer.data(kl + 1202);
    const auto *kl_1203 = buffer.data(kl + 1203);
    const auto *kl_1204 = buffer.data(kl + 1204);
    const auto *kl_1205 = buffer.data(kl + 1205);
    const auto *kl_1206 = buffer.data(kl + 1206);
    const auto *kl_1207 = buffer.data(kl + 1207);
    const auto *kl_1208 = buffer.data(kl + 1208);
    const auto *kl_1209 = buffer.data(kl + 1209);
    const auto *kl_1210 = buffer.data(kl + 1210);
    const auto *kl_1211 = buffer.data(kl + 1211);
    const auto *kl_1212 = buffer.data(kl + 1212);
    const auto *kl_1213 = buffer.data(kl + 1213);
    const auto *kl_1214 = buffer.data(kl + 1214);
    const auto *kl_1260 = buffer.data(kl + 1260);
    const auto *kl_1261 = buffer.data(kl + 1261);
    const auto *kl_1262 = buffer.data(kl + 1262);
    const auto *kl_1263 = buffer.data(kl + 1263);
    const auto *kl_1264 = buffer.data(kl + 1264);
    const auto *kl_1265 = buffer.data(kl + 1265);
    const auto *kl_1266 = buffer.data(kl + 1266);
    const auto *kl_1267 = buffer.data(kl + 1267);
    const auto *kl_1268 = buffer.data(kl + 1268);
    const auto *kl_1269 = buffer.data(kl + 1269);
    const auto *kl_1270 = buffer.data(kl + 1270);
    const auto *kl_1271 = buffer.data(kl + 1271);
    const auto *kl_1272 = buffer.data(kl + 1272);
    const auto *kl_1273 = buffer.data(kl + 1273);
    const auto *kl_1274 = buffer.data(kl + 1274);
    const auto *kl_1275 = buffer.data(kl + 1275);
    const auto *kl_1276 = buffer.data(kl + 1276);
    const auto *kl_1277 = buffer.data(kl + 1277);
    const auto *kl_1278 = buffer.data(kl + 1278);
    const auto *kl_1279 = buffer.data(kl + 1279);
    const auto *kl_1280 = buffer.data(kl + 1280);
    const auto *kl_1281 = buffer.data(kl + 1281);
    const auto *kl_1282 = buffer.data(kl + 1282);
    const auto *kl_1283 = buffer.data(kl + 1283);
    const auto *kl_1284 = buffer.data(kl + 1284);
    const auto *kl_1285 = buffer.data(kl + 1285);
    const auto *kl_1286 = buffer.data(kl + 1286);
    const auto *kl_1287 = buffer.data(kl + 1287);
    const auto *kl_1288 = buffer.data(kl + 1288);
    const auto *kl_1289 = buffer.data(kl + 1289);
    const auto *kl_1290 = buffer.data(kl + 1290);
    const auto *kl_1291 = buffer.data(kl + 1291);
    const auto *kl_1292 = buffer.data(kl + 1292);
    const auto *kl_1293 = buffer.data(kl + 1293);
    const auto *kl_1294 = buffer.data(kl + 1294);
    const auto *kl_1295 = buffer.data(kl + 1295);
    const auto *kl_1296 = buffer.data(kl + 1296);
    const auto *kl_1297 = buffer.data(kl + 1297);
    const auto *kl_1298 = buffer.data(kl + 1298);
    const auto *kl_1299 = buffer.data(kl + 1299);
    const auto *kl_1300 = buffer.data(kl + 1300);
    const auto *kl_1301 = buffer.data(kl + 1301);
    const auto *kl_1302 = buffer.data(kl + 1302);
    const auto *kl_1303 = buffer.data(kl + 1303);
    const auto *kl_1304 = buffer.data(kl + 1304);

#pragma omp simd aligned(t_824, t_825, t_826, t_827, t_828, hl_599, hl_600, hl_601, hl_602, \
                         hl_603, kl_1094, kl_1095, kl_1096, kl_1097, \
                         kl_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = -2.0 * hl_599[k]
                   + f_0 * kl_1094[k];

        t_825[k] = -2.0 * hl_600[k]
                   + f_0 * kl_1095[k];

        t_826[k] = -2.0 * hl_601[k]
                   + f_0 * kl_1096[k];

        t_827[k] = -2.0 * hl_602[k]
                   + f_0 * kl_1097[k];

        t_828[k] = -2.0 * hl_603[k]
                   + f_0 * kl_1098[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, t_833, hl_604, hl_605, hl_606, hl_607, \
                         hl_608, kl_1099, kl_1100, kl_1101, kl_1102, \
                         kl_1103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = -2.0 * hl_604[k]
                   + f_0 * kl_1099[k];

        t_830[k] = -2.0 * hl_605[k]
                   + f_0 * kl_1100[k];

        t_831[k] = -2.0 * hl_606[k]
                   + f_0 * kl_1101[k];

        t_832[k] = -2.0 * hl_607[k]
                   + f_0 * kl_1102[k];

        t_833[k] = -2.0 * hl_608[k]
                   + f_0 * kl_1103[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, t_837, t_838, hl_609, hl_610, hl_611, hl_612, \
                         hl_613, kl_1104, kl_1105, kl_1106, kl_1107, \
                         kl_1108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = -2.0 * hl_609[k]
                   + f_0 * kl_1104[k];

        t_835[k] = -2.0 * hl_610[k]
                   + f_0 * kl_1105[k];

        t_836[k] = -2.0 * hl_611[k]
                   + f_0 * kl_1106[k];

        t_837[k] = -2.0 * hl_612[k]
                   + f_0 * kl_1107[k];

        t_838[k] = -2.0 * hl_613[k]
                   + f_0 * kl_1108[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, t_842, t_843, hl_614, hl_615, hl_616, hl_617, \
                         hl_618, kl_1109, kl_1110, kl_1111, kl_1112, \
                         kl_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = -2.0 * hl_614[k]
                   + f_0 * kl_1109[k];

        t_840[k] = -2.0 * hl_615[k]
                   + f_0 * kl_1110[k];

        t_841[k] = -2.0 * hl_616[k]
                   + f_0 * kl_1111[k];

        t_842[k] = -2.0 * hl_617[k]
                   + f_0 * kl_1112[k];

        t_843[k] = -2.0 * hl_618[k]
                   + f_0 * kl_1113[k];
    }

#pragma omp simd aligned(t_844, t_845, t_846, t_847, t_848, hl_619, hl_620, hl_621, hl_622, \
                         hl_623, kl_1114, kl_1115, kl_1116, kl_1117, \
                         kl_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_844[k] = -2.0 * hl_619[k]
                   + f_0 * kl_1114[k];

        t_845[k] = -2.0 * hl_620[k]
                   + f_0 * kl_1115[k];

        t_846[k] = -2.0 * hl_621[k]
                   + f_0 * kl_1116[k];

        t_847[k] = -2.0 * hl_622[k]
                   + f_0 * kl_1117[k];

        t_848[k] = -2.0 * hl_623[k]
                   + f_0 * kl_1118[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, t_853, hl_624, hl_625, hl_626, hl_627, \
                         hl_628, kl_1119, kl_1120, kl_1121, kl_1122, \
                         kl_1123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = -2.0 * hl_624[k]
                   + f_0 * kl_1119[k];

        t_850[k] = -2.0 * hl_625[k]
                   + f_0 * kl_1120[k];

        t_851[k] = -2.0 * hl_626[k]
                   + f_0 * kl_1121[k];

        t_852[k] = -2.0 * hl_627[k]
                   + f_0 * kl_1122[k];

        t_853[k] = -2.0 * hl_628[k]
                   + f_0 * kl_1123[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, t_858, hl_629, hl_630, hl_631, hl_632, \
                         hl_633, kl_1124, kl_1125, kl_1126, kl_1127, \
                         kl_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = -2.0 * hl_629[k]
                   + f_0 * kl_1124[k];

        t_855[k] = -hl_630[k]
                   + f_0 * kl_1125[k];

        t_856[k] = -hl_631[k]
                   + f_0 * kl_1126[k];

        t_857[k] = -hl_632[k]
                   + f_0 * kl_1127[k];

        t_858[k] = -hl_633[k]
                   + f_0 * kl_1128[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, t_862, t_863, hl_634, hl_635, hl_636, hl_637, \
                         hl_638, kl_1129, kl_1130, kl_1131, kl_1132, \
                         kl_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = -hl_634[k]
                   + f_0 * kl_1129[k];

        t_860[k] = -hl_635[k]
                   + f_0 * kl_1130[k];

        t_861[k] = -hl_636[k]
                   + f_0 * kl_1131[k];

        t_862[k] = -hl_637[k]
                   + f_0 * kl_1132[k];

        t_863[k] = -hl_638[k]
                   + f_0 * kl_1133[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, t_867, t_868, hl_639, hl_640, hl_641, hl_642, \
                         hl_643, kl_1134, kl_1135, kl_1136, kl_1137, \
                         kl_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = -hl_639[k]
                   + f_0 * kl_1134[k];

        t_865[k] = -hl_640[k]
                   + f_0 * kl_1135[k];

        t_866[k] = -hl_641[k]
                   + f_0 * kl_1136[k];

        t_867[k] = -hl_642[k]
                   + f_0 * kl_1137[k];

        t_868[k] = -hl_643[k]
                   + f_0 * kl_1138[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, t_873, hl_644, hl_645, hl_646, hl_647, \
                         hl_648, kl_1139, kl_1140, kl_1141, kl_1142, \
                         kl_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = -hl_644[k]
                   + f_0 * kl_1139[k];

        t_870[k] = -hl_645[k]
                   + f_0 * kl_1140[k];

        t_871[k] = -hl_646[k]
                   + f_0 * kl_1141[k];

        t_872[k] = -hl_647[k]
                   + f_0 * kl_1142[k];

        t_873[k] = -hl_648[k]
                   + f_0 * kl_1143[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, t_878, hl_649, hl_650, hl_651, hl_652, \
                         hl_653, kl_1144, kl_1145, kl_1146, kl_1147, \
                         kl_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = -hl_649[k]
                   + f_0 * kl_1144[k];

        t_875[k] = -hl_650[k]
                   + f_0 * kl_1145[k];

        t_876[k] = -hl_651[k]
                   + f_0 * kl_1146[k];

        t_877[k] = -hl_652[k]
                   + f_0 * kl_1147[k];

        t_878[k] = -hl_653[k]
                   + f_0 * kl_1148[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, t_882, t_883, hl_654, hl_655, hl_656, hl_657, \
                         hl_658, kl_1149, kl_1150, kl_1151, kl_1152, \
                         kl_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = -hl_654[k]
                   + f_0 * kl_1149[k];

        t_880[k] = -hl_655[k]
                   + f_0 * kl_1150[k];

        t_881[k] = -hl_656[k]
                   + f_0 * kl_1151[k];

        t_882[k] = -hl_657[k]
                   + f_0 * kl_1152[k];

        t_883[k] = -hl_658[k]
                   + f_0 * kl_1153[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, t_888, hl_659, hl_660, hl_661, hl_662, \
                         hl_663, kl_1154, kl_1155, kl_1156, kl_1157, \
                         kl_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = -hl_659[k]
                   + f_0 * kl_1154[k];

        t_885[k] = -hl_660[k]
                   + f_0 * kl_1155[k];

        t_886[k] = -hl_661[k]
                   + f_0 * kl_1156[k];

        t_887[k] = -hl_662[k]
                   + f_0 * kl_1157[k];

        t_888[k] = -hl_663[k]
                   + f_0 * kl_1158[k];
    }

#pragma omp simd aligned(t_889, t_890, t_891, t_892, t_893, hl_664, hl_665, hl_666, hl_667, \
                         hl_668, kl_1159, kl_1160, kl_1161, kl_1162, \
                         kl_1163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = -hl_664[k]
                   + f_0 * kl_1159[k];

        t_890[k] = -hl_665[k]
                   + f_0 * kl_1160[k];

        t_891[k] = -hl_666[k]
                   + f_0 * kl_1161[k];

        t_892[k] = -hl_667[k]
                   + f_0 * kl_1162[k];

        t_893[k] = -hl_668[k]
                   + f_0 * kl_1163[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, t_898, hl_669, hl_670, hl_671, hl_672, \
                         hl_673, kl_1164, kl_1165, kl_1166, kl_1167, \
                         kl_1168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = -hl_669[k]
                   + f_0 * kl_1164[k];

        t_895[k] = -hl_670[k]
                   + f_0 * kl_1165[k];

        t_896[k] = -hl_671[k]
                   + f_0 * kl_1166[k];

        t_897[k] = -hl_672[k]
                   + f_0 * kl_1167[k];

        t_898[k] = -hl_673[k]
                   + f_0 * kl_1168[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, t_903, t_904, t_905, hl_674, kl_1169, \
                         kl_1170, kl_1171, kl_1172, kl_1173, kl_1174, \
                         kl_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = -hl_674[k]
                   + f_0 * kl_1169[k];

        t_900[k] = f_0 * kl_1170[k];

        t_901[k] = f_0 * kl_1171[k];

        t_902[k] = f_0 * kl_1172[k];

        t_903[k] = f_0 * kl_1173[k];

        t_904[k] = f_0 * kl_1174[k];

        t_905[k] = f_0 * kl_1175[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, t_909, t_910, t_911, t_912, t_913, kl_1176, \
                         kl_1177, kl_1178, kl_1179, kl_1180, kl_1181, kl_1182, \
                         kl_1183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_0 * kl_1176[k];

        t_907[k] = f_0 * kl_1177[k];

        t_908[k] = f_0 * kl_1178[k];

        t_909[k] = f_0 * kl_1179[k];

        t_910[k] = f_0 * kl_1180[k];

        t_911[k] = f_0 * kl_1181[k];

        t_912[k] = f_0 * kl_1182[k];

        t_913[k] = f_0 * kl_1183[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, t_918, t_919, t_920, t_921, kl_1184, \
                         kl_1185, kl_1186, kl_1187, kl_1188, kl_1189, kl_1190, \
                         kl_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_0 * kl_1184[k];

        t_915[k] = f_0 * kl_1185[k];

        t_916[k] = f_0 * kl_1186[k];

        t_917[k] = f_0 * kl_1187[k];

        t_918[k] = f_0 * kl_1188[k];

        t_919[k] = f_0 * kl_1189[k];

        t_920[k] = f_0 * kl_1190[k];

        t_921[k] = f_0 * kl_1191[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, t_926, t_927, t_928, t_929, kl_1192, \
                         kl_1193, kl_1194, kl_1195, kl_1196, kl_1197, kl_1198, \
                         kl_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_0 * kl_1192[k];

        t_923[k] = f_0 * kl_1193[k];

        t_924[k] = f_0 * kl_1194[k];

        t_925[k] = f_0 * kl_1195[k];

        t_926[k] = f_0 * kl_1196[k];

        t_927[k] = f_0 * kl_1197[k];

        t_928[k] = f_0 * kl_1198[k];

        t_929[k] = f_0 * kl_1199[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, t_935, t_936, t_937, kl_1200, \
                         kl_1201, kl_1202, kl_1203, kl_1204, kl_1205, kl_1206, \
                         kl_1207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_0 * kl_1200[k];

        t_931[k] = f_0 * kl_1201[k];

        t_932[k] = f_0 * kl_1202[k];

        t_933[k] = f_0 * kl_1203[k];

        t_934[k] = f_0 * kl_1204[k];

        t_935[k] = f_0 * kl_1205[k];

        t_936[k] = f_0 * kl_1206[k];

        t_937[k] = f_0 * kl_1207[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, t_942, t_943, t_944, kl_1208, kl_1209, \
                         kl_1210, kl_1211, kl_1212, kl_1213, kl_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_0 * kl_1208[k];

        t_939[k] = f_0 * kl_1209[k];

        t_940[k] = f_0 * kl_1210[k];

        t_941[k] = f_0 * kl_1211[k];

        t_942[k] = f_0 * kl_1212[k];

        t_943[k] = f_0 * kl_1213[k];

        t_944[k] = f_0 * kl_1214[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, t_949, hl_675, hl_676, hl_677, hl_678, \
                         hl_679, kl_1260, kl_1261, kl_1262, kl_1263, \
                         kl_1264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = -6.0 * hl_675[k]
                   + f_0 * kl_1260[k];

        t_946[k] = -6.0 * hl_676[k]
                   + f_0 * kl_1261[k];

        t_947[k] = -6.0 * hl_677[k]
                   + f_0 * kl_1262[k];

        t_948[k] = -6.0 * hl_678[k]
                   + f_0 * kl_1263[k];

        t_949[k] = -6.0 * hl_679[k]
                   + f_0 * kl_1264[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, t_954, hl_680, hl_681, hl_682, hl_683, \
                         hl_684, kl_1265, kl_1266, kl_1267, kl_1268, \
                         kl_1269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = -6.0 * hl_680[k]
                   + f_0 * kl_1265[k];

        t_951[k] = -6.0 * hl_681[k]
                   + f_0 * kl_1266[k];

        t_952[k] = -6.0 * hl_682[k]
                   + f_0 * kl_1267[k];

        t_953[k] = -6.0 * hl_683[k]
                   + f_0 * kl_1268[k];

        t_954[k] = -6.0 * hl_684[k]
                   + f_0 * kl_1269[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, t_959, hl_685, hl_686, hl_687, hl_688, \
                         hl_689, kl_1270, kl_1271, kl_1272, kl_1273, \
                         kl_1274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = -6.0 * hl_685[k]
                   + f_0 * kl_1270[k];

        t_956[k] = -6.0 * hl_686[k]
                   + f_0 * kl_1271[k];

        t_957[k] = -6.0 * hl_687[k]
                   + f_0 * kl_1272[k];

        t_958[k] = -6.0 * hl_688[k]
                   + f_0 * kl_1273[k];

        t_959[k] = -6.0 * hl_689[k]
                   + f_0 * kl_1274[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, t_963, t_964, hl_690, hl_691, hl_692, hl_693, \
                         hl_694, kl_1275, kl_1276, kl_1277, kl_1278, \
                         kl_1279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = -6.0 * hl_690[k]
                   + f_0 * kl_1275[k];

        t_961[k] = -6.0 * hl_691[k]
                   + f_0 * kl_1276[k];

        t_962[k] = -6.0 * hl_692[k]
                   + f_0 * kl_1277[k];

        t_963[k] = -6.0 * hl_693[k]
                   + f_0 * kl_1278[k];

        t_964[k] = -6.0 * hl_694[k]
                   + f_0 * kl_1279[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, t_969, hl_695, hl_696, hl_697, hl_698, \
                         hl_699, kl_1280, kl_1281, kl_1282, kl_1283, \
                         kl_1284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = -6.0 * hl_695[k]
                   + f_0 * kl_1280[k];

        t_966[k] = -6.0 * hl_696[k]
                   + f_0 * kl_1281[k];

        t_967[k] = -6.0 * hl_697[k]
                   + f_0 * kl_1282[k];

        t_968[k] = -6.0 * hl_698[k]
                   + f_0 * kl_1283[k];

        t_969[k] = -6.0 * hl_699[k]
                   + f_0 * kl_1284[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, t_974, hl_700, hl_701, hl_702, hl_703, \
                         hl_704, kl_1285, kl_1286, kl_1287, kl_1288, \
                         kl_1289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = -6.0 * hl_700[k]
                   + f_0 * kl_1285[k];

        t_971[k] = -6.0 * hl_701[k]
                   + f_0 * kl_1286[k];

        t_972[k] = -6.0 * hl_702[k]
                   + f_0 * kl_1287[k];

        t_973[k] = -6.0 * hl_703[k]
                   + f_0 * kl_1288[k];

        t_974[k] = -6.0 * hl_704[k]
                   + f_0 * kl_1289[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, t_979, hl_705, hl_706, hl_707, hl_708, \
                         hl_709, kl_1290, kl_1291, kl_1292, kl_1293, \
                         kl_1294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = -6.0 * hl_705[k]
                   + f_0 * kl_1290[k];

        t_976[k] = -6.0 * hl_706[k]
                   + f_0 * kl_1291[k];

        t_977[k] = -6.0 * hl_707[k]
                   + f_0 * kl_1292[k];

        t_978[k] = -6.0 * hl_708[k]
                   + f_0 * kl_1293[k];

        t_979[k] = -6.0 * hl_709[k]
                   + f_0 * kl_1294[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, t_984, hl_710, hl_711, hl_712, hl_713, \
                         hl_714, kl_1295, kl_1296, kl_1297, kl_1298, \
                         kl_1299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = -6.0 * hl_710[k]
                   + f_0 * kl_1295[k];

        t_981[k] = -6.0 * hl_711[k]
                   + f_0 * kl_1296[k];

        t_982[k] = -6.0 * hl_712[k]
                   + f_0 * kl_1297[k];

        t_983[k] = -6.0 * hl_713[k]
                   + f_0 * kl_1298[k];

        t_984[k] = -6.0 * hl_714[k]
                   + f_0 * kl_1299[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, t_989, hl_715, hl_716, hl_717, hl_718, \
                         hl_719, kl_1300, kl_1301, kl_1302, kl_1303, \
                         kl_1304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = -6.0 * hl_715[k]
                   + f_0 * kl_1300[k];

        t_986[k] = -6.0 * hl_716[k]
                   + f_0 * kl_1301[k];

        t_987[k] = -6.0 * hl_717[k]
                   + f_0 * kl_1302[k];

        t_988[k] = -6.0 * hl_718[k]
                   + f_0 * kl_1303[k];

        t_989[k] = -6.0 * hl_719[k]
                   + f_0 * kl_1304[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_1_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hl_720 = buffer.data(hl + 720);
    const auto *hl_721 = buffer.data(hl + 721);
    const auto *hl_722 = buffer.data(hl + 722);
    const auto *hl_723 = buffer.data(hl + 723);
    const auto *hl_724 = buffer.data(hl + 724);
    const auto *hl_725 = buffer.data(hl + 725);
    const auto *hl_726 = buffer.data(hl + 726);
    const auto *hl_727 = buffer.data(hl + 727);
    const auto *hl_728 = buffer.data(hl + 728);
    const auto *hl_729 = buffer.data(hl + 729);
    const auto *hl_730 = buffer.data(hl + 730);
    const auto *hl_731 = buffer.data(hl + 731);
    const auto *hl_732 = buffer.data(hl + 732);
    const auto *hl_733 = buffer.data(hl + 733);
    const auto *hl_734 = buffer.data(hl + 734);
    const auto *hl_735 = buffer.data(hl + 735);
    const auto *hl_736 = buffer.data(hl + 736);
    const auto *hl_737 = buffer.data(hl + 737);
    const auto *hl_738 = buffer.data(hl + 738);
    const auto *hl_739 = buffer.data(hl + 739);
    const auto *hl_740 = buffer.data(hl + 740);
    const auto *hl_741 = buffer.data(hl + 741);
    const auto *hl_742 = buffer.data(hl + 742);
    const auto *hl_743 = buffer.data(hl + 743);
    const auto *hl_744 = buffer.data(hl + 744);
    const auto *hl_745 = buffer.data(hl + 745);
    const auto *hl_746 = buffer.data(hl + 746);
    const auto *hl_747 = buffer.data(hl + 747);
    const auto *hl_748 = buffer.data(hl + 748);
    const auto *hl_749 = buffer.data(hl + 749);
    const auto *hl_750 = buffer.data(hl + 750);
    const auto *hl_751 = buffer.data(hl + 751);
    const auto *hl_752 = buffer.data(hl + 752);
    const auto *hl_753 = buffer.data(hl + 753);
    const auto *hl_754 = buffer.data(hl + 754);
    const auto *hl_755 = buffer.data(hl + 755);
    const auto *hl_756 = buffer.data(hl + 756);
    const auto *hl_757 = buffer.data(hl + 757);
    const auto *hl_758 = buffer.data(hl + 758);
    const auto *hl_759 = buffer.data(hl + 759);
    const auto *hl_760 = buffer.data(hl + 760);
    const auto *hl_761 = buffer.data(hl + 761);
    const auto *hl_762 = buffer.data(hl + 762);
    const auto *hl_763 = buffer.data(hl + 763);
    const auto *hl_764 = buffer.data(hl + 764);
    const auto *hl_765 = buffer.data(hl + 765);
    const auto *hl_766 = buffer.data(hl + 766);
    const auto *hl_767 = buffer.data(hl + 767);
    const auto *hl_768 = buffer.data(hl + 768);
    const auto *hl_769 = buffer.data(hl + 769);
    const auto *hl_770 = buffer.data(hl + 770);
    const auto *hl_771 = buffer.data(hl + 771);
    const auto *hl_772 = buffer.data(hl + 772);
    const auto *hl_773 = buffer.data(hl + 773);
    const auto *hl_774 = buffer.data(hl + 774);
    const auto *hl_775 = buffer.data(hl + 775);
    const auto *hl_776 = buffer.data(hl + 776);
    const auto *hl_777 = buffer.data(hl + 777);
    const auto *hl_778 = buffer.data(hl + 778);
    const auto *hl_779 = buffer.data(hl + 779);
    const auto *hl_780 = buffer.data(hl + 780);
    const auto *hl_781 = buffer.data(hl + 781);
    const auto *hl_782 = buffer.data(hl + 782);
    const auto *hl_783 = buffer.data(hl + 783);
    const auto *hl_784 = buffer.data(hl + 784);
    const auto *hl_785 = buffer.data(hl + 785);
    const auto *hl_786 = buffer.data(hl + 786);
    const auto *hl_787 = buffer.data(hl + 787);
    const auto *hl_788 = buffer.data(hl + 788);
    const auto *hl_789 = buffer.data(hl + 789);
    const auto *hl_790 = buffer.data(hl + 790);
    const auto *hl_791 = buffer.data(hl + 791);
    const auto *hl_792 = buffer.data(hl + 792);
    const auto *hl_793 = buffer.data(hl + 793);
    const auto *hl_794 = buffer.data(hl + 794);
    const auto *hl_795 = buffer.data(hl + 795);
    const auto *hl_796 = buffer.data(hl + 796);
    const auto *hl_797 = buffer.data(hl + 797);
    const auto *hl_798 = buffer.data(hl + 798);
    const auto *hl_799 = buffer.data(hl + 799);
    const auto *hl_800 = buffer.data(hl + 800);
    const auto *hl_801 = buffer.data(hl + 801);
    const auto *hl_802 = buffer.data(hl + 802);
    const auto *hl_803 = buffer.data(hl + 803);
    const auto *hl_804 = buffer.data(hl + 804);
    const auto *hl_805 = buffer.data(hl + 805);
    const auto *hl_806 = buffer.data(hl + 806);
    const auto *hl_807 = buffer.data(hl + 807);
    const auto *hl_808 = buffer.data(hl + 808);
    const auto *hl_809 = buffer.data(hl + 809);
    const auto *hl_810 = buffer.data(hl + 810);
    const auto *hl_811 = buffer.data(hl + 811);
    const auto *hl_812 = buffer.data(hl + 812);
    const auto *hl_813 = buffer.data(hl + 813);
    const auto *hl_814 = buffer.data(hl + 814);
    const auto *hl_815 = buffer.data(hl + 815);
    const auto *hl_816 = buffer.data(hl + 816);
    const auto *hl_817 = buffer.data(hl + 817);
    const auto *hl_818 = buffer.data(hl + 818);
    const auto *hl_819 = buffer.data(hl + 819);
    const auto *hl_820 = buffer.data(hl + 820);
    const auto *hl_821 = buffer.data(hl + 821);
    const auto *hl_822 = buffer.data(hl + 822);
    const auto *hl_823 = buffer.data(hl + 823);
    const auto *hl_824 = buffer.data(hl + 824);
    const auto *hl_825 = buffer.data(hl + 825);
    const auto *hl_826 = buffer.data(hl + 826);
    const auto *hl_827 = buffer.data(hl + 827);
    const auto *hl_828 = buffer.data(hl + 828);
    const auto *hl_829 = buffer.data(hl + 829);
    const auto *hl_830 = buffer.data(hl + 830);
    const auto *hl_831 = buffer.data(hl + 831);
    const auto *hl_832 = buffer.data(hl + 832);
    const auto *hl_833 = buffer.data(hl + 833);
    const auto *hl_834 = buffer.data(hl + 834);
    const auto *hl_835 = buffer.data(hl + 835);
    const auto *hl_836 = buffer.data(hl + 836);
    const auto *hl_837 = buffer.data(hl + 837);
    const auto *hl_838 = buffer.data(hl + 838);
    const auto *hl_839 = buffer.data(hl + 839);
    const auto *hl_840 = buffer.data(hl + 840);
    const auto *hl_841 = buffer.data(hl + 841);
    const auto *hl_842 = buffer.data(hl + 842);
    const auto *hl_843 = buffer.data(hl + 843);
    const auto *hl_844 = buffer.data(hl + 844);
    const auto *hl_845 = buffer.data(hl + 845);
    const auto *hl_846 = buffer.data(hl + 846);
    const auto *hl_847 = buffer.data(hl + 847);
    const auto *hl_848 = buffer.data(hl + 848);
    const auto *hl_849 = buffer.data(hl + 849);
    const auto *hl_850 = buffer.data(hl + 850);
    const auto *hl_851 = buffer.data(hl + 851);
    const auto *hl_852 = buffer.data(hl + 852);
    const auto *hl_853 = buffer.data(hl + 853);
    const auto *hl_854 = buffer.data(hl + 854);
    const auto *hl_855 = buffer.data(hl + 855);
    const auto *hl_856 = buffer.data(hl + 856);
    const auto *hl_857 = buffer.data(hl + 857);
    const auto *hl_858 = buffer.data(hl + 858);
    const auto *hl_859 = buffer.data(hl + 859);
    const auto *hl_860 = buffer.data(hl + 860);
    const auto *hl_861 = buffer.data(hl + 861);
    const auto *hl_862 = buffer.data(hl + 862);
    const auto *hl_863 = buffer.data(hl + 863);
    const auto *hl_864 = buffer.data(hl + 864);
    const auto *hl_865 = buffer.data(hl + 865);
    const auto *hl_866 = buffer.data(hl + 866);
    const auto *hl_867 = buffer.data(hl + 867);
    const auto *hl_868 = buffer.data(hl + 868);
    const auto *hl_869 = buffer.data(hl + 869);

    const auto *kl_1305 = buffer.data(kl + 1305);
    const auto *kl_1306 = buffer.data(kl + 1306);
    const auto *kl_1307 = buffer.data(kl + 1307);
    const auto *kl_1308 = buffer.data(kl + 1308);
    const auto *kl_1309 = buffer.data(kl + 1309);
    const auto *kl_1310 = buffer.data(kl + 1310);
    const auto *kl_1311 = buffer.data(kl + 1311);
    const auto *kl_1312 = buffer.data(kl + 1312);
    const auto *kl_1313 = buffer.data(kl + 1313);
    const auto *kl_1314 = buffer.data(kl + 1314);
    const auto *kl_1315 = buffer.data(kl + 1315);
    const auto *kl_1316 = buffer.data(kl + 1316);
    const auto *kl_1317 = buffer.data(kl + 1317);
    const auto *kl_1318 = buffer.data(kl + 1318);
    const auto *kl_1319 = buffer.data(kl + 1319);
    const auto *kl_1320 = buffer.data(kl + 1320);
    const auto *kl_1321 = buffer.data(kl + 1321);
    const auto *kl_1322 = buffer.data(kl + 1322);
    const auto *kl_1323 = buffer.data(kl + 1323);
    const auto *kl_1324 = buffer.data(kl + 1324);
    const auto *kl_1325 = buffer.data(kl + 1325);
    const auto *kl_1326 = buffer.data(kl + 1326);
    const auto *kl_1327 = buffer.data(kl + 1327);
    const auto *kl_1328 = buffer.data(kl + 1328);
    const auto *kl_1329 = buffer.data(kl + 1329);
    const auto *kl_1330 = buffer.data(kl + 1330);
    const auto *kl_1331 = buffer.data(kl + 1331);
    const auto *kl_1332 = buffer.data(kl + 1332);
    const auto *kl_1333 = buffer.data(kl + 1333);
    const auto *kl_1334 = buffer.data(kl + 1334);
    const auto *kl_1335 = buffer.data(kl + 1335);
    const auto *kl_1336 = buffer.data(kl + 1336);
    const auto *kl_1337 = buffer.data(kl + 1337);
    const auto *kl_1338 = buffer.data(kl + 1338);
    const auto *kl_1339 = buffer.data(kl + 1339);
    const auto *kl_1340 = buffer.data(kl + 1340);
    const auto *kl_1341 = buffer.data(kl + 1341);
    const auto *kl_1342 = buffer.data(kl + 1342);
    const auto *kl_1343 = buffer.data(kl + 1343);
    const auto *kl_1344 = buffer.data(kl + 1344);
    const auto *kl_1345 = buffer.data(kl + 1345);
    const auto *kl_1346 = buffer.data(kl + 1346);
    const auto *kl_1347 = buffer.data(kl + 1347);
    const auto *kl_1348 = buffer.data(kl + 1348);
    const auto *kl_1349 = buffer.data(kl + 1349);
    const auto *kl_1350 = buffer.data(kl + 1350);
    const auto *kl_1351 = buffer.data(kl + 1351);
    const auto *kl_1352 = buffer.data(kl + 1352);
    const auto *kl_1353 = buffer.data(kl + 1353);
    const auto *kl_1354 = buffer.data(kl + 1354);
    const auto *kl_1355 = buffer.data(kl + 1355);
    const auto *kl_1356 = buffer.data(kl + 1356);
    const auto *kl_1357 = buffer.data(kl + 1357);
    const auto *kl_1358 = buffer.data(kl + 1358);
    const auto *kl_1359 = buffer.data(kl + 1359);
    const auto *kl_1360 = buffer.data(kl + 1360);
    const auto *kl_1361 = buffer.data(kl + 1361);
    const auto *kl_1362 = buffer.data(kl + 1362);
    const auto *kl_1363 = buffer.data(kl + 1363);
    const auto *kl_1364 = buffer.data(kl + 1364);
    const auto *kl_1365 = buffer.data(kl + 1365);
    const auto *kl_1366 = buffer.data(kl + 1366);
    const auto *kl_1367 = buffer.data(kl + 1367);
    const auto *kl_1368 = buffer.data(kl + 1368);
    const auto *kl_1369 = buffer.data(kl + 1369);
    const auto *kl_1370 = buffer.data(kl + 1370);
    const auto *kl_1371 = buffer.data(kl + 1371);
    const auto *kl_1372 = buffer.data(kl + 1372);
    const auto *kl_1373 = buffer.data(kl + 1373);
    const auto *kl_1374 = buffer.data(kl + 1374);
    const auto *kl_1375 = buffer.data(kl + 1375);
    const auto *kl_1376 = buffer.data(kl + 1376);
    const auto *kl_1377 = buffer.data(kl + 1377);
    const auto *kl_1378 = buffer.data(kl + 1378);
    const auto *kl_1379 = buffer.data(kl + 1379);
    const auto *kl_1380 = buffer.data(kl + 1380);
    const auto *kl_1381 = buffer.data(kl + 1381);
    const auto *kl_1382 = buffer.data(kl + 1382);
    const auto *kl_1383 = buffer.data(kl + 1383);
    const auto *kl_1384 = buffer.data(kl + 1384);
    const auto *kl_1385 = buffer.data(kl + 1385);
    const auto *kl_1386 = buffer.data(kl + 1386);
    const auto *kl_1387 = buffer.data(kl + 1387);
    const auto *kl_1388 = buffer.data(kl + 1388);
    const auto *kl_1389 = buffer.data(kl + 1389);
    const auto *kl_1390 = buffer.data(kl + 1390);
    const auto *kl_1391 = buffer.data(kl + 1391);
    const auto *kl_1392 = buffer.data(kl + 1392);
    const auto *kl_1393 = buffer.data(kl + 1393);
    const auto *kl_1394 = buffer.data(kl + 1394);
    const auto *kl_1395 = buffer.data(kl + 1395);
    const auto *kl_1396 = buffer.data(kl + 1396);
    const auto *kl_1397 = buffer.data(kl + 1397);
    const auto *kl_1398 = buffer.data(kl + 1398);
    const auto *kl_1399 = buffer.data(kl + 1399);
    const auto *kl_1400 = buffer.data(kl + 1400);
    const auto *kl_1401 = buffer.data(kl + 1401);
    const auto *kl_1402 = buffer.data(kl + 1402);
    const auto *kl_1403 = buffer.data(kl + 1403);
    const auto *kl_1404 = buffer.data(kl + 1404);
    const auto *kl_1405 = buffer.data(kl + 1405);
    const auto *kl_1406 = buffer.data(kl + 1406);
    const auto *kl_1407 = buffer.data(kl + 1407);
    const auto *kl_1408 = buffer.data(kl + 1408);
    const auto *kl_1409 = buffer.data(kl + 1409);
    const auto *kl_1410 = buffer.data(kl + 1410);
    const auto *kl_1411 = buffer.data(kl + 1411);
    const auto *kl_1412 = buffer.data(kl + 1412);
    const auto *kl_1413 = buffer.data(kl + 1413);
    const auto *kl_1414 = buffer.data(kl + 1414);
    const auto *kl_1415 = buffer.data(kl + 1415);
    const auto *kl_1416 = buffer.data(kl + 1416);
    const auto *kl_1417 = buffer.data(kl + 1417);
    const auto *kl_1418 = buffer.data(kl + 1418);
    const auto *kl_1419 = buffer.data(kl + 1419);
    const auto *kl_1420 = buffer.data(kl + 1420);
    const auto *kl_1421 = buffer.data(kl + 1421);
    const auto *kl_1422 = buffer.data(kl + 1422);
    const auto *kl_1423 = buffer.data(kl + 1423);
    const auto *kl_1424 = buffer.data(kl + 1424);
    const auto *kl_1425 = buffer.data(kl + 1425);
    const auto *kl_1426 = buffer.data(kl + 1426);
    const auto *kl_1427 = buffer.data(kl + 1427);
    const auto *kl_1428 = buffer.data(kl + 1428);
    const auto *kl_1429 = buffer.data(kl + 1429);
    const auto *kl_1430 = buffer.data(kl + 1430);
    const auto *kl_1431 = buffer.data(kl + 1431);
    const auto *kl_1432 = buffer.data(kl + 1432);
    const auto *kl_1433 = buffer.data(kl + 1433);
    const auto *kl_1434 = buffer.data(kl + 1434);
    const auto *kl_1435 = buffer.data(kl + 1435);
    const auto *kl_1436 = buffer.data(kl + 1436);
    const auto *kl_1437 = buffer.data(kl + 1437);
    const auto *kl_1438 = buffer.data(kl + 1438);
    const auto *kl_1439 = buffer.data(kl + 1439);
    const auto *kl_1440 = buffer.data(kl + 1440);
    const auto *kl_1441 = buffer.data(kl + 1441);
    const auto *kl_1442 = buffer.data(kl + 1442);
    const auto *kl_1443 = buffer.data(kl + 1443);
    const auto *kl_1444 = buffer.data(kl + 1444);
    const auto *kl_1445 = buffer.data(kl + 1445);
    const auto *kl_1446 = buffer.data(kl + 1446);
    const auto *kl_1447 = buffer.data(kl + 1447);
    const auto *kl_1448 = buffer.data(kl + 1448);
    const auto *kl_1449 = buffer.data(kl + 1449);
    const auto *kl_1450 = buffer.data(kl + 1450);
    const auto *kl_1451 = buffer.data(kl + 1451);
    const auto *kl_1452 = buffer.data(kl + 1452);
    const auto *kl_1453 = buffer.data(kl + 1453);
    const auto *kl_1454 = buffer.data(kl + 1454);

#pragma omp simd aligned(t_990, t_991, t_992, t_993, t_994, hl_720, hl_721, hl_722, hl_723, \
                         hl_724, kl_1305, kl_1306, kl_1307, kl_1308, \
                         kl_1309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = -5.0 * hl_720[k]
                   + f_0 * kl_1305[k];

        t_991[k] = -5.0 * hl_721[k]
                   + f_0 * kl_1306[k];

        t_992[k] = -5.0 * hl_722[k]
                   + f_0 * kl_1307[k];

        t_993[k] = -5.0 * hl_723[k]
                   + f_0 * kl_1308[k];

        t_994[k] = -5.0 * hl_724[k]
                   + f_0 * kl_1309[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, hl_725, hl_726, hl_727, hl_728, \
                         hl_729, kl_1310, kl_1311, kl_1312, kl_1313, \
                         kl_1314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = -5.0 * hl_725[k]
                   + f_0 * kl_1310[k];

        t_996[k] = -5.0 * hl_726[k]
                   + f_0 * kl_1311[k];

        t_997[k] = -5.0 * hl_727[k]
                   + f_0 * kl_1312[k];

        t_998[k] = -5.0 * hl_728[k]
                   + f_0 * kl_1313[k];

        t_999[k] = -5.0 * hl_729[k]
                   + f_0 * kl_1314[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, t_1004, hl_730, hl_731, hl_732, \
                         hl_733, hl_734, kl_1315, kl_1316, kl_1317, kl_1318, \
                         kl_1319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = -5.0 * hl_730[k]
                    + f_0 * kl_1315[k];

        t_1001[k] = -5.0 * hl_731[k]
                    + f_0 * kl_1316[k];

        t_1002[k] = -5.0 * hl_732[k]
                    + f_0 * kl_1317[k];

        t_1003[k] = -5.0 * hl_733[k]
                    + f_0 * kl_1318[k];

        t_1004[k] = -5.0 * hl_734[k]
                    + f_0 * kl_1319[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, t_1008, t_1009, hl_735, hl_736, hl_737, \
                         hl_738, hl_739, kl_1320, kl_1321, kl_1322, kl_1323, \
                         kl_1324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = -5.0 * hl_735[k]
                    + f_0 * kl_1320[k];

        t_1006[k] = -5.0 * hl_736[k]
                    + f_0 * kl_1321[k];

        t_1007[k] = -5.0 * hl_737[k]
                    + f_0 * kl_1322[k];

        t_1008[k] = -5.0 * hl_738[k]
                    + f_0 * kl_1323[k];

        t_1009[k] = -5.0 * hl_739[k]
                    + f_0 * kl_1324[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, t_1013, t_1014, hl_740, hl_741, hl_742, \
                         hl_743, hl_744, kl_1325, kl_1326, kl_1327, kl_1328, \
                         kl_1329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = -5.0 * hl_740[k]
                    + f_0 * kl_1325[k];

        t_1011[k] = -5.0 * hl_741[k]
                    + f_0 * kl_1326[k];

        t_1012[k] = -5.0 * hl_742[k]
                    + f_0 * kl_1327[k];

        t_1013[k] = -5.0 * hl_743[k]
                    + f_0 * kl_1328[k];

        t_1014[k] = -5.0 * hl_744[k]
                    + f_0 * kl_1329[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, t_1018, t_1019, hl_745, hl_746, hl_747, \
                         hl_748, hl_749, kl_1330, kl_1331, kl_1332, kl_1333, \
                         kl_1334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = -5.0 * hl_745[k]
                    + f_0 * kl_1330[k];

        t_1016[k] = -5.0 * hl_746[k]
                    + f_0 * kl_1331[k];

        t_1017[k] = -5.0 * hl_747[k]
                    + f_0 * kl_1332[k];

        t_1018[k] = -5.0 * hl_748[k]
                    + f_0 * kl_1333[k];

        t_1019[k] = -5.0 * hl_749[k]
                    + f_0 * kl_1334[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, t_1024, hl_750, hl_751, hl_752, \
                         hl_753, hl_754, kl_1335, kl_1336, kl_1337, kl_1338, \
                         kl_1339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = -5.0 * hl_750[k]
                    + f_0 * kl_1335[k];

        t_1021[k] = -5.0 * hl_751[k]
                    + f_0 * kl_1336[k];

        t_1022[k] = -5.0 * hl_752[k]
                    + f_0 * kl_1337[k];

        t_1023[k] = -5.0 * hl_753[k]
                    + f_0 * kl_1338[k];

        t_1024[k] = -5.0 * hl_754[k]
                    + f_0 * kl_1339[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, t_1029, hl_755, hl_756, hl_757, \
                         hl_758, hl_759, kl_1340, kl_1341, kl_1342, kl_1343, \
                         kl_1344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = -5.0 * hl_755[k]
                    + f_0 * kl_1340[k];

        t_1026[k] = -5.0 * hl_756[k]
                    + f_0 * kl_1341[k];

        t_1027[k] = -5.0 * hl_757[k]
                    + f_0 * kl_1342[k];

        t_1028[k] = -5.0 * hl_758[k]
                    + f_0 * kl_1343[k];

        t_1029[k] = -5.0 * hl_759[k]
                    + f_0 * kl_1344[k];
    }

#pragma omp simd aligned(t_1030, t_1031, t_1032, t_1033, t_1034, hl_760, hl_761, hl_762, \
                         hl_763, hl_764, kl_1345, kl_1346, kl_1347, kl_1348, \
                         kl_1349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1030[k] = -5.0 * hl_760[k]
                    + f_0 * kl_1345[k];

        t_1031[k] = -5.0 * hl_761[k]
                    + f_0 * kl_1346[k];

        t_1032[k] = -5.0 * hl_762[k]
                    + f_0 * kl_1347[k];

        t_1033[k] = -5.0 * hl_763[k]
                    + f_0 * kl_1348[k];

        t_1034[k] = -5.0 * hl_764[k]
                    + f_0 * kl_1349[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, t_1038, t_1039, hl_765, hl_766, hl_767, \
                         hl_768, hl_769, kl_1350, kl_1351, kl_1352, kl_1353, \
                         kl_1354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = -4.0 * hl_765[k]
                    + f_0 * kl_1350[k];

        t_1036[k] = -4.0 * hl_766[k]
                    + f_0 * kl_1351[k];

        t_1037[k] = -4.0 * hl_767[k]
                    + f_0 * kl_1352[k];

        t_1038[k] = -4.0 * hl_768[k]
                    + f_0 * kl_1353[k];

        t_1039[k] = -4.0 * hl_769[k]
                    + f_0 * kl_1354[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, t_1043, t_1044, hl_770, hl_771, hl_772, \
                         hl_773, hl_774, kl_1355, kl_1356, kl_1357, kl_1358, \
                         kl_1359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = -4.0 * hl_770[k]
                    + f_0 * kl_1355[k];

        t_1041[k] = -4.0 * hl_771[k]
                    + f_0 * kl_1356[k];

        t_1042[k] = -4.0 * hl_772[k]
                    + f_0 * kl_1357[k];

        t_1043[k] = -4.0 * hl_773[k]
                    + f_0 * kl_1358[k];

        t_1044[k] = -4.0 * hl_774[k]
                    + f_0 * kl_1359[k];
    }

#pragma omp simd aligned(t_1045, t_1046, t_1047, t_1048, t_1049, hl_775, hl_776, hl_777, \
                         hl_778, hl_779, kl_1360, kl_1361, kl_1362, kl_1363, \
                         kl_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1045[k] = -4.0 * hl_775[k]
                    + f_0 * kl_1360[k];

        t_1046[k] = -4.0 * hl_776[k]
                    + f_0 * kl_1361[k];

        t_1047[k] = -4.0 * hl_777[k]
                    + f_0 * kl_1362[k];

        t_1048[k] = -4.0 * hl_778[k]
                    + f_0 * kl_1363[k];

        t_1049[k] = -4.0 * hl_779[k]
                    + f_0 * kl_1364[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, t_1054, hl_780, hl_781, hl_782, \
                         hl_783, hl_784, kl_1365, kl_1366, kl_1367, kl_1368, \
                         kl_1369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = -4.0 * hl_780[k]
                    + f_0 * kl_1365[k];

        t_1051[k] = -4.0 * hl_781[k]
                    + f_0 * kl_1366[k];

        t_1052[k] = -4.0 * hl_782[k]
                    + f_0 * kl_1367[k];

        t_1053[k] = -4.0 * hl_783[k]
                    + f_0 * kl_1368[k];

        t_1054[k] = -4.0 * hl_784[k]
                    + f_0 * kl_1369[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, t_1058, t_1059, hl_785, hl_786, hl_787, \
                         hl_788, hl_789, kl_1370, kl_1371, kl_1372, kl_1373, \
                         kl_1374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = -4.0 * hl_785[k]
                    + f_0 * kl_1370[k];

        t_1056[k] = -4.0 * hl_786[k]
                    + f_0 * kl_1371[k];

        t_1057[k] = -4.0 * hl_787[k]
                    + f_0 * kl_1372[k];

        t_1058[k] = -4.0 * hl_788[k]
                    + f_0 * kl_1373[k];

        t_1059[k] = -4.0 * hl_789[k]
                    + f_0 * kl_1374[k];
    }

#pragma omp simd aligned(t_1060, t_1061, t_1062, t_1063, t_1064, hl_790, hl_791, hl_792, \
                         hl_793, hl_794, kl_1375, kl_1376, kl_1377, kl_1378, \
                         kl_1379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1060[k] = -4.0 * hl_790[k]
                    + f_0 * kl_1375[k];

        t_1061[k] = -4.0 * hl_791[k]
                    + f_0 * kl_1376[k];

        t_1062[k] = -4.0 * hl_792[k]
                    + f_0 * kl_1377[k];

        t_1063[k] = -4.0 * hl_793[k]
                    + f_0 * kl_1378[k];

        t_1064[k] = -4.0 * hl_794[k]
                    + f_0 * kl_1379[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, t_1069, hl_795, hl_796, hl_797, \
                         hl_798, hl_799, kl_1380, kl_1381, kl_1382, kl_1383, \
                         kl_1384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = -4.0 * hl_795[k]
                    + f_0 * kl_1380[k];

        t_1066[k] = -4.0 * hl_796[k]
                    + f_0 * kl_1381[k];

        t_1067[k] = -4.0 * hl_797[k]
                    + f_0 * kl_1382[k];

        t_1068[k] = -4.0 * hl_798[k]
                    + f_0 * kl_1383[k];

        t_1069[k] = -4.0 * hl_799[k]
                    + f_0 * kl_1384[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, t_1073, t_1074, hl_800, hl_801, hl_802, \
                         hl_803, hl_804, kl_1385, kl_1386, kl_1387, kl_1388, \
                         kl_1389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = -4.0 * hl_800[k]
                    + f_0 * kl_1385[k];

        t_1071[k] = -4.0 * hl_801[k]
                    + f_0 * kl_1386[k];

        t_1072[k] = -4.0 * hl_802[k]
                    + f_0 * kl_1387[k];

        t_1073[k] = -4.0 * hl_803[k]
                    + f_0 * kl_1388[k];

        t_1074[k] = -4.0 * hl_804[k]
                    + f_0 * kl_1389[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, t_1078, t_1079, hl_805, hl_806, hl_807, \
                         hl_808, hl_809, kl_1390, kl_1391, kl_1392, kl_1393, \
                         kl_1394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = -4.0 * hl_805[k]
                    + f_0 * kl_1390[k];

        t_1076[k] = -4.0 * hl_806[k]
                    + f_0 * kl_1391[k];

        t_1077[k] = -4.0 * hl_807[k]
                    + f_0 * kl_1392[k];

        t_1078[k] = -4.0 * hl_808[k]
                    + f_0 * kl_1393[k];

        t_1079[k] = -4.0 * hl_809[k]
                    + f_0 * kl_1394[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, t_1083, t_1084, hl_810, hl_811, hl_812, \
                         hl_813, hl_814, kl_1395, kl_1396, kl_1397, kl_1398, \
                         kl_1399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = -3.0 * hl_810[k]
                    + f_0 * kl_1395[k];

        t_1081[k] = -3.0 * hl_811[k]
                    + f_0 * kl_1396[k];

        t_1082[k] = -3.0 * hl_812[k]
                    + f_0 * kl_1397[k];

        t_1083[k] = -3.0 * hl_813[k]
                    + f_0 * kl_1398[k];

        t_1084[k] = -3.0 * hl_814[k]
                    + f_0 * kl_1399[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, t_1088, t_1089, hl_815, hl_816, hl_817, \
                         hl_818, hl_819, kl_1400, kl_1401, kl_1402, kl_1403, \
                         kl_1404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = -3.0 * hl_815[k]
                    + f_0 * kl_1400[k];

        t_1086[k] = -3.0 * hl_816[k]
                    + f_0 * kl_1401[k];

        t_1087[k] = -3.0 * hl_817[k]
                    + f_0 * kl_1402[k];

        t_1088[k] = -3.0 * hl_818[k]
                    + f_0 * kl_1403[k];

        t_1089[k] = -3.0 * hl_819[k]
                    + f_0 * kl_1404[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, t_1093, t_1094, hl_820, hl_821, hl_822, \
                         hl_823, hl_824, kl_1405, kl_1406, kl_1407, kl_1408, \
                         kl_1409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = -3.0 * hl_820[k]
                    + f_0 * kl_1405[k];

        t_1091[k] = -3.0 * hl_821[k]
                    + f_0 * kl_1406[k];

        t_1092[k] = -3.0 * hl_822[k]
                    + f_0 * kl_1407[k];

        t_1093[k] = -3.0 * hl_823[k]
                    + f_0 * kl_1408[k];

        t_1094[k] = -3.0 * hl_824[k]
                    + f_0 * kl_1409[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, t_1098, t_1099, hl_825, hl_826, hl_827, \
                         hl_828, hl_829, kl_1410, kl_1411, kl_1412, kl_1413, \
                         kl_1414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = -3.0 * hl_825[k]
                    + f_0 * kl_1410[k];

        t_1096[k] = -3.0 * hl_826[k]
                    + f_0 * kl_1411[k];

        t_1097[k] = -3.0 * hl_827[k]
                    + f_0 * kl_1412[k];

        t_1098[k] = -3.0 * hl_828[k]
                    + f_0 * kl_1413[k];

        t_1099[k] = -3.0 * hl_829[k]
                    + f_0 * kl_1414[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, t_1104, hl_830, hl_831, hl_832, \
                         hl_833, hl_834, kl_1415, kl_1416, kl_1417, kl_1418, \
                         kl_1419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = -3.0 * hl_830[k]
                    + f_0 * kl_1415[k];

        t_1101[k] = -3.0 * hl_831[k]
                    + f_0 * kl_1416[k];

        t_1102[k] = -3.0 * hl_832[k]
                    + f_0 * kl_1417[k];

        t_1103[k] = -3.0 * hl_833[k]
                    + f_0 * kl_1418[k];

        t_1104[k] = -3.0 * hl_834[k]
                    + f_0 * kl_1419[k];
    }

#pragma omp simd aligned(t_1105, t_1106, t_1107, t_1108, t_1109, hl_835, hl_836, hl_837, \
                         hl_838, hl_839, kl_1420, kl_1421, kl_1422, kl_1423, \
                         kl_1424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1105[k] = -3.0 * hl_835[k]
                    + f_0 * kl_1420[k];

        t_1106[k] = -3.0 * hl_836[k]
                    + f_0 * kl_1421[k];

        t_1107[k] = -3.0 * hl_837[k]
                    + f_0 * kl_1422[k];

        t_1108[k] = -3.0 * hl_838[k]
                    + f_0 * kl_1423[k];

        t_1109[k] = -3.0 * hl_839[k]
                    + f_0 * kl_1424[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, t_1114, hl_840, hl_841, hl_842, \
                         hl_843, hl_844, kl_1425, kl_1426, kl_1427, kl_1428, \
                         kl_1429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = -3.0 * hl_840[k]
                    + f_0 * kl_1425[k];

        t_1111[k] = -3.0 * hl_841[k]
                    + f_0 * kl_1426[k];

        t_1112[k] = -3.0 * hl_842[k]
                    + f_0 * kl_1427[k];

        t_1113[k] = -3.0 * hl_843[k]
                    + f_0 * kl_1428[k];

        t_1114[k] = -3.0 * hl_844[k]
                    + f_0 * kl_1429[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, t_1118, t_1119, hl_845, hl_846, hl_847, \
                         hl_848, hl_849, kl_1430, kl_1431, kl_1432, kl_1433, \
                         kl_1434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = -3.0 * hl_845[k]
                    + f_0 * kl_1430[k];

        t_1116[k] = -3.0 * hl_846[k]
                    + f_0 * kl_1431[k];

        t_1117[k] = -3.0 * hl_847[k]
                    + f_0 * kl_1432[k];

        t_1118[k] = -3.0 * hl_848[k]
                    + f_0 * kl_1433[k];

        t_1119[k] = -3.0 * hl_849[k]
                    + f_0 * kl_1434[k];
    }

#pragma omp simd aligned(t_1120, t_1121, t_1122, t_1123, t_1124, hl_850, hl_851, hl_852, \
                         hl_853, hl_854, kl_1435, kl_1436, kl_1437, kl_1438, \
                         kl_1439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1120[k] = -3.0 * hl_850[k]
                    + f_0 * kl_1435[k];

        t_1121[k] = -3.0 * hl_851[k]
                    + f_0 * kl_1436[k];

        t_1122[k] = -3.0 * hl_852[k]
                    + f_0 * kl_1437[k];

        t_1123[k] = -3.0 * hl_853[k]
                    + f_0 * kl_1438[k];

        t_1124[k] = -3.0 * hl_854[k]
                    + f_0 * kl_1439[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, t_1128, t_1129, hl_855, hl_856, hl_857, \
                         hl_858, hl_859, kl_1440, kl_1441, kl_1442, kl_1443, \
                         kl_1444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = -2.0 * hl_855[k]
                    + f_0 * kl_1440[k];

        t_1126[k] = -2.0 * hl_856[k]
                    + f_0 * kl_1441[k];

        t_1127[k] = -2.0 * hl_857[k]
                    + f_0 * kl_1442[k];

        t_1128[k] = -2.0 * hl_858[k]
                    + f_0 * kl_1443[k];

        t_1129[k] = -2.0 * hl_859[k]
                    + f_0 * kl_1444[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, t_1133, t_1134, hl_860, hl_861, hl_862, \
                         hl_863, hl_864, kl_1445, kl_1446, kl_1447, kl_1448, \
                         kl_1449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = -2.0 * hl_860[k]
                    + f_0 * kl_1445[k];

        t_1131[k] = -2.0 * hl_861[k]
                    + f_0 * kl_1446[k];

        t_1132[k] = -2.0 * hl_862[k]
                    + f_0 * kl_1447[k];

        t_1133[k] = -2.0 * hl_863[k]
                    + f_0 * kl_1448[k];

        t_1134[k] = -2.0 * hl_864[k]
                    + f_0 * kl_1449[k];
    }

#pragma omp simd aligned(t_1135, t_1136, t_1137, t_1138, t_1139, hl_865, hl_866, hl_867, \
                         hl_868, hl_869, kl_1450, kl_1451, kl_1452, kl_1453, \
                         kl_1454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1135[k] = -2.0 * hl_865[k]
                    + f_0 * kl_1450[k];

        t_1136[k] = -2.0 * hl_866[k]
                    + f_0 * kl_1451[k];

        t_1137[k] = -2.0 * hl_867[k]
                    + f_0 * kl_1452[k];

        t_1138[k] = -2.0 * hl_868[k]
                    + f_0 * kl_1453[k];

        t_1139[k] = -2.0 * hl_869[k]
                    + f_0 * kl_1454[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_1_piece7(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hl_870 = buffer.data(hl + 870);
    const auto *hl_871 = buffer.data(hl + 871);
    const auto *hl_872 = buffer.data(hl + 872);
    const auto *hl_873 = buffer.data(hl + 873);
    const auto *hl_874 = buffer.data(hl + 874);
    const auto *hl_875 = buffer.data(hl + 875);
    const auto *hl_876 = buffer.data(hl + 876);
    const auto *hl_877 = buffer.data(hl + 877);
    const auto *hl_878 = buffer.data(hl + 878);
    const auto *hl_879 = buffer.data(hl + 879);
    const auto *hl_880 = buffer.data(hl + 880);
    const auto *hl_881 = buffer.data(hl + 881);
    const auto *hl_882 = buffer.data(hl + 882);
    const auto *hl_883 = buffer.data(hl + 883);
    const auto *hl_884 = buffer.data(hl + 884);
    const auto *hl_885 = buffer.data(hl + 885);
    const auto *hl_886 = buffer.data(hl + 886);
    const auto *hl_887 = buffer.data(hl + 887);
    const auto *hl_888 = buffer.data(hl + 888);
    const auto *hl_889 = buffer.data(hl + 889);
    const auto *hl_890 = buffer.data(hl + 890);
    const auto *hl_891 = buffer.data(hl + 891);
    const auto *hl_892 = buffer.data(hl + 892);
    const auto *hl_893 = buffer.data(hl + 893);
    const auto *hl_894 = buffer.data(hl + 894);
    const auto *hl_895 = buffer.data(hl + 895);
    const auto *hl_896 = buffer.data(hl + 896);
    const auto *hl_897 = buffer.data(hl + 897);
    const auto *hl_898 = buffer.data(hl + 898);
    const auto *hl_899 = buffer.data(hl + 899);
    const auto *hl_900 = buffer.data(hl + 900);
    const auto *hl_901 = buffer.data(hl + 901);
    const auto *hl_902 = buffer.data(hl + 902);
    const auto *hl_903 = buffer.data(hl + 903);
    const auto *hl_904 = buffer.data(hl + 904);
    const auto *hl_905 = buffer.data(hl + 905);
    const auto *hl_906 = buffer.data(hl + 906);
    const auto *hl_907 = buffer.data(hl + 907);
    const auto *hl_908 = buffer.data(hl + 908);
    const auto *hl_909 = buffer.data(hl + 909);
    const auto *hl_910 = buffer.data(hl + 910);
    const auto *hl_911 = buffer.data(hl + 911);
    const auto *hl_912 = buffer.data(hl + 912);
    const auto *hl_913 = buffer.data(hl + 913);
    const auto *hl_914 = buffer.data(hl + 914);
    const auto *hl_915 = buffer.data(hl + 915);
    const auto *hl_916 = buffer.data(hl + 916);
    const auto *hl_917 = buffer.data(hl + 917);
    const auto *hl_918 = buffer.data(hl + 918);
    const auto *hl_919 = buffer.data(hl + 919);
    const auto *hl_920 = buffer.data(hl + 920);
    const auto *hl_921 = buffer.data(hl + 921);
    const auto *hl_922 = buffer.data(hl + 922);
    const auto *hl_923 = buffer.data(hl + 923);
    const auto *hl_924 = buffer.data(hl + 924);
    const auto *hl_925 = buffer.data(hl + 925);
    const auto *hl_926 = buffer.data(hl + 926);
    const auto *hl_927 = buffer.data(hl + 927);
    const auto *hl_928 = buffer.data(hl + 928);
    const auto *hl_929 = buffer.data(hl + 929);
    const auto *hl_930 = buffer.data(hl + 930);
    const auto *hl_931 = buffer.data(hl + 931);
    const auto *hl_932 = buffer.data(hl + 932);
    const auto *hl_933 = buffer.data(hl + 933);
    const auto *hl_934 = buffer.data(hl + 934);
    const auto *hl_935 = buffer.data(hl + 935);
    const auto *hl_936 = buffer.data(hl + 936);
    const auto *hl_937 = buffer.data(hl + 937);
    const auto *hl_938 = buffer.data(hl + 938);
    const auto *hl_939 = buffer.data(hl + 939);
    const auto *hl_940 = buffer.data(hl + 940);
    const auto *hl_941 = buffer.data(hl + 941);
    const auto *hl_942 = buffer.data(hl + 942);
    const auto *hl_943 = buffer.data(hl + 943);
    const auto *hl_944 = buffer.data(hl + 944);

    const auto *kl_1455 = buffer.data(kl + 1455);
    const auto *kl_1456 = buffer.data(kl + 1456);
    const auto *kl_1457 = buffer.data(kl + 1457);
    const auto *kl_1458 = buffer.data(kl + 1458);
    const auto *kl_1459 = buffer.data(kl + 1459);
    const auto *kl_1460 = buffer.data(kl + 1460);
    const auto *kl_1461 = buffer.data(kl + 1461);
    const auto *kl_1462 = buffer.data(kl + 1462);
    const auto *kl_1463 = buffer.data(kl + 1463);
    const auto *kl_1464 = buffer.data(kl + 1464);
    const auto *kl_1465 = buffer.data(kl + 1465);
    const auto *kl_1466 = buffer.data(kl + 1466);
    const auto *kl_1467 = buffer.data(kl + 1467);
    const auto *kl_1468 = buffer.data(kl + 1468);
    const auto *kl_1469 = buffer.data(kl + 1469);
    const auto *kl_1470 = buffer.data(kl + 1470);
    const auto *kl_1471 = buffer.data(kl + 1471);
    const auto *kl_1472 = buffer.data(kl + 1472);
    const auto *kl_1473 = buffer.data(kl + 1473);
    const auto *kl_1474 = buffer.data(kl + 1474);
    const auto *kl_1475 = buffer.data(kl + 1475);
    const auto *kl_1476 = buffer.data(kl + 1476);
    const auto *kl_1477 = buffer.data(kl + 1477);
    const auto *kl_1478 = buffer.data(kl + 1478);
    const auto *kl_1479 = buffer.data(kl + 1479);
    const auto *kl_1480 = buffer.data(kl + 1480);
    const auto *kl_1481 = buffer.data(kl + 1481);
    const auto *kl_1482 = buffer.data(kl + 1482);
    const auto *kl_1483 = buffer.data(kl + 1483);
    const auto *kl_1484 = buffer.data(kl + 1484);
    const auto *kl_1485 = buffer.data(kl + 1485);
    const auto *kl_1486 = buffer.data(kl + 1486);
    const auto *kl_1487 = buffer.data(kl + 1487);
    const auto *kl_1488 = buffer.data(kl + 1488);
    const auto *kl_1489 = buffer.data(kl + 1489);
    const auto *kl_1490 = buffer.data(kl + 1490);
    const auto *kl_1491 = buffer.data(kl + 1491);
    const auto *kl_1492 = buffer.data(kl + 1492);
    const auto *kl_1493 = buffer.data(kl + 1493);
    const auto *kl_1494 = buffer.data(kl + 1494);
    const auto *kl_1495 = buffer.data(kl + 1495);
    const auto *kl_1496 = buffer.data(kl + 1496);
    const auto *kl_1497 = buffer.data(kl + 1497);
    const auto *kl_1498 = buffer.data(kl + 1498);
    const auto *kl_1499 = buffer.data(kl + 1499);
    const auto *kl_1500 = buffer.data(kl + 1500);
    const auto *kl_1501 = buffer.data(kl + 1501);
    const auto *kl_1502 = buffer.data(kl + 1502);
    const auto *kl_1503 = buffer.data(kl + 1503);
    const auto *kl_1504 = buffer.data(kl + 1504);
    const auto *kl_1505 = buffer.data(kl + 1505);
    const auto *kl_1506 = buffer.data(kl + 1506);
    const auto *kl_1507 = buffer.data(kl + 1507);
    const auto *kl_1508 = buffer.data(kl + 1508);
    const auto *kl_1509 = buffer.data(kl + 1509);
    const auto *kl_1510 = buffer.data(kl + 1510);
    const auto *kl_1511 = buffer.data(kl + 1511);
    const auto *kl_1512 = buffer.data(kl + 1512);
    const auto *kl_1513 = buffer.data(kl + 1513);
    const auto *kl_1514 = buffer.data(kl + 1514);
    const auto *kl_1515 = buffer.data(kl + 1515);
    const auto *kl_1516 = buffer.data(kl + 1516);
    const auto *kl_1517 = buffer.data(kl + 1517);
    const auto *kl_1518 = buffer.data(kl + 1518);
    const auto *kl_1519 = buffer.data(kl + 1519);
    const auto *kl_1520 = buffer.data(kl + 1520);
    const auto *kl_1521 = buffer.data(kl + 1521);
    const auto *kl_1522 = buffer.data(kl + 1522);
    const auto *kl_1523 = buffer.data(kl + 1523);
    const auto *kl_1524 = buffer.data(kl + 1524);
    const auto *kl_1525 = buffer.data(kl + 1525);
    const auto *kl_1526 = buffer.data(kl + 1526);
    const auto *kl_1527 = buffer.data(kl + 1527);
    const auto *kl_1528 = buffer.data(kl + 1528);
    const auto *kl_1529 = buffer.data(kl + 1529);
    const auto *kl_1530 = buffer.data(kl + 1530);
    const auto *kl_1531 = buffer.data(kl + 1531);
    const auto *kl_1532 = buffer.data(kl + 1532);
    const auto *kl_1533 = buffer.data(kl + 1533);
    const auto *kl_1534 = buffer.data(kl + 1534);
    const auto *kl_1535 = buffer.data(kl + 1535);
    const auto *kl_1536 = buffer.data(kl + 1536);
    const auto *kl_1537 = buffer.data(kl + 1537);
    const auto *kl_1538 = buffer.data(kl + 1538);
    const auto *kl_1539 = buffer.data(kl + 1539);
    const auto *kl_1540 = buffer.data(kl + 1540);
    const auto *kl_1541 = buffer.data(kl + 1541);
    const auto *kl_1542 = buffer.data(kl + 1542);
    const auto *kl_1543 = buffer.data(kl + 1543);
    const auto *kl_1544 = buffer.data(kl + 1544);
    const auto *kl_1545 = buffer.data(kl + 1545);
    const auto *kl_1546 = buffer.data(kl + 1546);
    const auto *kl_1547 = buffer.data(kl + 1547);
    const auto *kl_1548 = buffer.data(kl + 1548);
    const auto *kl_1549 = buffer.data(kl + 1549);
    const auto *kl_1550 = buffer.data(kl + 1550);
    const auto *kl_1551 = buffer.data(kl + 1551);
    const auto *kl_1552 = buffer.data(kl + 1552);
    const auto *kl_1553 = buffer.data(kl + 1553);
    const auto *kl_1554 = buffer.data(kl + 1554);
    const auto *kl_1555 = buffer.data(kl + 1555);
    const auto *kl_1556 = buffer.data(kl + 1556);
    const auto *kl_1557 = buffer.data(kl + 1557);
    const auto *kl_1558 = buffer.data(kl + 1558);
    const auto *kl_1559 = buffer.data(kl + 1559);
    const auto *kl_1560 = buffer.data(kl + 1560);
    const auto *kl_1561 = buffer.data(kl + 1561);
    const auto *kl_1562 = buffer.data(kl + 1562);
    const auto *kl_1563 = buffer.data(kl + 1563);
    const auto *kl_1564 = buffer.data(kl + 1564);
    const auto *kl_1565 = buffer.data(kl + 1565);
    const auto *kl_1566 = buffer.data(kl + 1566);
    const auto *kl_1567 = buffer.data(kl + 1567);
    const auto *kl_1568 = buffer.data(kl + 1568);
    const auto *kl_1569 = buffer.data(kl + 1569);
    const auto *kl_1570 = buffer.data(kl + 1570);
    const auto *kl_1571 = buffer.data(kl + 1571);
    const auto *kl_1572 = buffer.data(kl + 1572);
    const auto *kl_1573 = buffer.data(kl + 1573);
    const auto *kl_1574 = buffer.data(kl + 1574);

#pragma omp simd aligned(t_1140, t_1141, t_1142, t_1143, t_1144, hl_870, hl_871, hl_872, \
                         hl_873, hl_874, kl_1455, kl_1456, kl_1457, kl_1458, \
                         kl_1459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = -2.0 * hl_870[k]
                    + f_0 * kl_1455[k];

        t_1141[k] = -2.0 * hl_871[k]
                    + f_0 * kl_1456[k];

        t_1142[k] = -2.0 * hl_872[k]
                    + f_0 * kl_1457[k];

        t_1143[k] = -2.0 * hl_873[k]
                    + f_0 * kl_1458[k];

        t_1144[k] = -2.0 * hl_874[k]
                    + f_0 * kl_1459[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, t_1148, t_1149, hl_875, hl_876, hl_877, \
                         hl_878, hl_879, kl_1460, kl_1461, kl_1462, kl_1463, \
                         kl_1464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = -2.0 * hl_875[k]
                    + f_0 * kl_1460[k];

        t_1146[k] = -2.0 * hl_876[k]
                    + f_0 * kl_1461[k];

        t_1147[k] = -2.0 * hl_877[k]
                    + f_0 * kl_1462[k];

        t_1148[k] = -2.0 * hl_878[k]
                    + f_0 * kl_1463[k];

        t_1149[k] = -2.0 * hl_879[k]
                    + f_0 * kl_1464[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, t_1153, t_1154, hl_880, hl_881, hl_882, \
                         hl_883, hl_884, kl_1465, kl_1466, kl_1467, kl_1468, \
                         kl_1469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = -2.0 * hl_880[k]
                    + f_0 * kl_1465[k];

        t_1151[k] = -2.0 * hl_881[k]
                    + f_0 * kl_1466[k];

        t_1152[k] = -2.0 * hl_882[k]
                    + f_0 * kl_1467[k];

        t_1153[k] = -2.0 * hl_883[k]
                    + f_0 * kl_1468[k];

        t_1154[k] = -2.0 * hl_884[k]
                    + f_0 * kl_1469[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, t_1159, hl_885, hl_886, hl_887, \
                         hl_888, hl_889, kl_1470, kl_1471, kl_1472, kl_1473, \
                         kl_1474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = -2.0 * hl_885[k]
                    + f_0 * kl_1470[k];

        t_1156[k] = -2.0 * hl_886[k]
                    + f_0 * kl_1471[k];

        t_1157[k] = -2.0 * hl_887[k]
                    + f_0 * kl_1472[k];

        t_1158[k] = -2.0 * hl_888[k]
                    + f_0 * kl_1473[k];

        t_1159[k] = -2.0 * hl_889[k]
                    + f_0 * kl_1474[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, t_1163, t_1164, hl_890, hl_891, hl_892, \
                         hl_893, hl_894, kl_1475, kl_1476, kl_1477, kl_1478, \
                         kl_1479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = -2.0 * hl_890[k]
                    + f_0 * kl_1475[k];

        t_1161[k] = -2.0 * hl_891[k]
                    + f_0 * kl_1476[k];

        t_1162[k] = -2.0 * hl_892[k]
                    + f_0 * kl_1477[k];

        t_1163[k] = -2.0 * hl_893[k]
                    + f_0 * kl_1478[k];

        t_1164[k] = -2.0 * hl_894[k]
                    + f_0 * kl_1479[k];
    }

#pragma omp simd aligned(t_1165, t_1166, t_1167, t_1168, t_1169, hl_895, hl_896, hl_897, \
                         hl_898, hl_899, kl_1480, kl_1481, kl_1482, kl_1483, \
                         kl_1484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1165[k] = -2.0 * hl_895[k]
                    + f_0 * kl_1480[k];

        t_1166[k] = -2.0 * hl_896[k]
                    + f_0 * kl_1481[k];

        t_1167[k] = -2.0 * hl_897[k]
                    + f_0 * kl_1482[k];

        t_1168[k] = -2.0 * hl_898[k]
                    + f_0 * kl_1483[k];

        t_1169[k] = -2.0 * hl_899[k]
                    + f_0 * kl_1484[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, t_1173, t_1174, hl_900, hl_901, hl_902, \
                         hl_903, hl_904, kl_1485, kl_1486, kl_1487, kl_1488, \
                         kl_1489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = -hl_900[k]
                    + f_0 * kl_1485[k];

        t_1171[k] = -hl_901[k]
                    + f_0 * kl_1486[k];

        t_1172[k] = -hl_902[k]
                    + f_0 * kl_1487[k];

        t_1173[k] = -hl_903[k]
                    + f_0 * kl_1488[k];

        t_1174[k] = -hl_904[k]
                    + f_0 * kl_1489[k];
    }

#pragma omp simd aligned(t_1175, t_1176, t_1177, t_1178, t_1179, hl_905, hl_906, hl_907, \
                         hl_908, hl_909, kl_1490, kl_1491, kl_1492, kl_1493, \
                         kl_1494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1175[k] = -hl_905[k]
                    + f_0 * kl_1490[k];

        t_1176[k] = -hl_906[k]
                    + f_0 * kl_1491[k];

        t_1177[k] = -hl_907[k]
                    + f_0 * kl_1492[k];

        t_1178[k] = -hl_908[k]
                    + f_0 * kl_1493[k];

        t_1179[k] = -hl_909[k]
                    + f_0 * kl_1494[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, t_1183, t_1184, hl_910, hl_911, hl_912, \
                         hl_913, hl_914, kl_1495, kl_1496, kl_1497, kl_1498, \
                         kl_1499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = -hl_910[k]
                    + f_0 * kl_1495[k];

        t_1181[k] = -hl_911[k]
                    + f_0 * kl_1496[k];

        t_1182[k] = -hl_912[k]
                    + f_0 * kl_1497[k];

        t_1183[k] = -hl_913[k]
                    + f_0 * kl_1498[k];

        t_1184[k] = -hl_914[k]
                    + f_0 * kl_1499[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, t_1188, t_1189, hl_915, hl_916, hl_917, \
                         hl_918, hl_919, kl_1500, kl_1501, kl_1502, kl_1503, \
                         kl_1504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = -hl_915[k]
                    + f_0 * kl_1500[k];

        t_1186[k] = -hl_916[k]
                    + f_0 * kl_1501[k];

        t_1187[k] = -hl_917[k]
                    + f_0 * kl_1502[k];

        t_1188[k] = -hl_918[k]
                    + f_0 * kl_1503[k];

        t_1189[k] = -hl_919[k]
                    + f_0 * kl_1504[k];
    }

#pragma omp simd aligned(t_1190, t_1191, t_1192, t_1193, t_1194, hl_920, hl_921, hl_922, \
                         hl_923, hl_924, kl_1505, kl_1506, kl_1507, kl_1508, \
                         kl_1509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1190[k] = -hl_920[k]
                    + f_0 * kl_1505[k];

        t_1191[k] = -hl_921[k]
                    + f_0 * kl_1506[k];

        t_1192[k] = -hl_922[k]
                    + f_0 * kl_1507[k];

        t_1193[k] = -hl_923[k]
                    + f_0 * kl_1508[k];

        t_1194[k] = -hl_924[k]
                    + f_0 * kl_1509[k];
    }

#pragma omp simd aligned(t_1195, t_1196, t_1197, t_1198, t_1199, hl_925, hl_926, hl_927, \
                         hl_928, hl_929, kl_1510, kl_1511, kl_1512, kl_1513, \
                         kl_1514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1195[k] = -hl_925[k]
                    + f_0 * kl_1510[k];

        t_1196[k] = -hl_926[k]
                    + f_0 * kl_1511[k];

        t_1197[k] = -hl_927[k]
                    + f_0 * kl_1512[k];

        t_1198[k] = -hl_928[k]
                    + f_0 * kl_1513[k];

        t_1199[k] = -hl_929[k]
                    + f_0 * kl_1514[k];
    }

#pragma omp simd aligned(t_1200, t_1201, t_1202, t_1203, t_1204, hl_930, hl_931, hl_932, \
                         hl_933, hl_934, kl_1515, kl_1516, kl_1517, kl_1518, \
                         kl_1519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1200[k] = -hl_930[k]
                    + f_0 * kl_1515[k];

        t_1201[k] = -hl_931[k]
                    + f_0 * kl_1516[k];

        t_1202[k] = -hl_932[k]
                    + f_0 * kl_1517[k];

        t_1203[k] = -hl_933[k]
                    + f_0 * kl_1518[k];

        t_1204[k] = -hl_934[k]
                    + f_0 * kl_1519[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, t_1208, t_1209, hl_935, hl_936, hl_937, \
                         hl_938, hl_939, kl_1520, kl_1521, kl_1522, kl_1523, \
                         kl_1524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = -hl_935[k]
                    + f_0 * kl_1520[k];

        t_1206[k] = -hl_936[k]
                    + f_0 * kl_1521[k];

        t_1207[k] = -hl_937[k]
                    + f_0 * kl_1522[k];

        t_1208[k] = -hl_938[k]
                    + f_0 * kl_1523[k];

        t_1209[k] = -hl_939[k]
                    + f_0 * kl_1524[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, t_1214, hl_940, hl_941, hl_942, \
                         hl_943, hl_944, kl_1525, kl_1526, kl_1527, kl_1528, \
                         kl_1529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = -hl_940[k]
                    + f_0 * kl_1525[k];

        t_1211[k] = -hl_941[k]
                    + f_0 * kl_1526[k];

        t_1212[k] = -hl_942[k]
                    + f_0 * kl_1527[k];

        t_1213[k] = -hl_943[k]
                    + f_0 * kl_1528[k];

        t_1214[k] = -hl_944[k]
                    + f_0 * kl_1529[k];
    }

#pragma omp simd aligned(t_1215, t_1216, t_1217, t_1218, t_1219, t_1220, t_1221, t_1222, \
                         kl_1530, kl_1531, kl_1532, kl_1533, kl_1534, kl_1535, kl_1536, \
                         kl_1537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1215[k] = f_0 * kl_1530[k];

        t_1216[k] = f_0 * kl_1531[k];

        t_1217[k] = f_0 * kl_1532[k];

        t_1218[k] = f_0 * kl_1533[k];

        t_1219[k] = f_0 * kl_1534[k];

        t_1220[k] = f_0 * kl_1535[k];

        t_1221[k] = f_0 * kl_1536[k];

        t_1222[k] = f_0 * kl_1537[k];
    }

#pragma omp simd aligned(t_1223, t_1224, t_1225, t_1226, t_1227, t_1228, t_1229, t_1230, \
                         kl_1538, kl_1539, kl_1540, kl_1541, kl_1542, kl_1543, kl_1544, \
                         kl_1545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1223[k] = f_0 * kl_1538[k];

        t_1224[k] = f_0 * kl_1539[k];

        t_1225[k] = f_0 * kl_1540[k];

        t_1226[k] = f_0 * kl_1541[k];

        t_1227[k] = f_0 * kl_1542[k];

        t_1228[k] = f_0 * kl_1543[k];

        t_1229[k] = f_0 * kl_1544[k];

        t_1230[k] = f_0 * kl_1545[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, t_1235, t_1236, t_1237, t_1238, \
                         kl_1546, kl_1547, kl_1548, kl_1549, kl_1550, kl_1551, kl_1552, \
                         kl_1553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = f_0 * kl_1546[k];

        t_1232[k] = f_0 * kl_1547[k];

        t_1233[k] = f_0 * kl_1548[k];

        t_1234[k] = f_0 * kl_1549[k];

        t_1235[k] = f_0 * kl_1550[k];

        t_1236[k] = f_0 * kl_1551[k];

        t_1237[k] = f_0 * kl_1552[k];

        t_1238[k] = f_0 * kl_1553[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, t_1242, t_1243, t_1244, t_1245, t_1246, \
                         kl_1554, kl_1555, kl_1556, kl_1557, kl_1558, kl_1559, kl_1560, \
                         kl_1561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = f_0 * kl_1554[k];

        t_1240[k] = f_0 * kl_1555[k];

        t_1241[k] = f_0 * kl_1556[k];

        t_1242[k] = f_0 * kl_1557[k];

        t_1243[k] = f_0 * kl_1558[k];

        t_1244[k] = f_0 * kl_1559[k];

        t_1245[k] = f_0 * kl_1560[k];

        t_1246[k] = f_0 * kl_1561[k];
    }

#pragma omp simd aligned(t_1247, t_1248, t_1249, t_1250, t_1251, t_1252, t_1253, t_1254, \
                         kl_1562, kl_1563, kl_1564, kl_1565, kl_1566, kl_1567, kl_1568, \
                         kl_1569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1247[k] = f_0 * kl_1562[k];

        t_1248[k] = f_0 * kl_1563[k];

        t_1249[k] = f_0 * kl_1564[k];

        t_1250[k] = f_0 * kl_1565[k];

        t_1251[k] = f_0 * kl_1566[k];

        t_1252[k] = f_0 * kl_1567[k];

        t_1253[k] = f_0 * kl_1568[k];

        t_1254[k] = f_0 * kl_1569[k];
    }

#pragma omp simd aligned(t_1255, t_1256, t_1257, t_1258, t_1259, kl_1570, kl_1571, kl_1572, \
                         kl_1573, kl_1574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1255[k] = f_0 * kl_1570[k];

        t_1256[k] = f_0 * kl_1571[k];

        t_1257[k] = f_0 * kl_1572[k];

        t_1258[k] = f_0 * kl_1573[k];

        t_1259[k] = f_0 * kl_1574[k];
    }
}

auto
compute_prim_geom_10_il_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t hl, const size_t kl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_il_electron_repulsion_1_piece0(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_1_piece1(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_1_piece2(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_1_piece3(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_1_piece4(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_1_piece5(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_1_piece6(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_1_piece7(buffer, target, hl, kl, ncols, alpha);
}

static auto
compute_prim_geom_10_il_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
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

    const auto *hl_0 = buffer.data(hl + 0);
    const auto *hl_1 = buffer.data(hl + 1);
    const auto *hl_2 = buffer.data(hl + 2);
    const auto *hl_3 = buffer.data(hl + 3);
    const auto *hl_4 = buffer.data(hl + 4);
    const auto *hl_5 = buffer.data(hl + 5);
    const auto *hl_6 = buffer.data(hl + 6);
    const auto *hl_7 = buffer.data(hl + 7);
    const auto *hl_8 = buffer.data(hl + 8);
    const auto *hl_9 = buffer.data(hl + 9);
    const auto *hl_10 = buffer.data(hl + 10);
    const auto *hl_11 = buffer.data(hl + 11);
    const auto *hl_12 = buffer.data(hl + 12);
    const auto *hl_13 = buffer.data(hl + 13);
    const auto *hl_14 = buffer.data(hl + 14);
    const auto *hl_15 = buffer.data(hl + 15);
    const auto *hl_16 = buffer.data(hl + 16);
    const auto *hl_17 = buffer.data(hl + 17);
    const auto *hl_18 = buffer.data(hl + 18);
    const auto *hl_19 = buffer.data(hl + 19);
    const auto *hl_20 = buffer.data(hl + 20);
    const auto *hl_21 = buffer.data(hl + 21);
    const auto *hl_22 = buffer.data(hl + 22);
    const auto *hl_23 = buffer.data(hl + 23);
    const auto *hl_24 = buffer.data(hl + 24);
    const auto *hl_25 = buffer.data(hl + 25);
    const auto *hl_26 = buffer.data(hl + 26);
    const auto *hl_27 = buffer.data(hl + 27);
    const auto *hl_28 = buffer.data(hl + 28);
    const auto *hl_29 = buffer.data(hl + 29);
    const auto *hl_30 = buffer.data(hl + 30);
    const auto *hl_31 = buffer.data(hl + 31);
    const auto *hl_32 = buffer.data(hl + 32);
    const auto *hl_33 = buffer.data(hl + 33);
    const auto *hl_34 = buffer.data(hl + 34);
    const auto *hl_35 = buffer.data(hl + 35);
    const auto *hl_36 = buffer.data(hl + 36);
    const auto *hl_37 = buffer.data(hl + 37);
    const auto *hl_38 = buffer.data(hl + 38);
    const auto *hl_39 = buffer.data(hl + 39);
    const auto *hl_40 = buffer.data(hl + 40);
    const auto *hl_41 = buffer.data(hl + 41);
    const auto *hl_42 = buffer.data(hl + 42);
    const auto *hl_43 = buffer.data(hl + 43);
    const auto *hl_44 = buffer.data(hl + 44);
    const auto *hl_45 = buffer.data(hl + 45);
    const auto *hl_46 = buffer.data(hl + 46);
    const auto *hl_47 = buffer.data(hl + 47);
    const auto *hl_48 = buffer.data(hl + 48);
    const auto *hl_49 = buffer.data(hl + 49);
    const auto *hl_50 = buffer.data(hl + 50);
    const auto *hl_51 = buffer.data(hl + 51);
    const auto *hl_52 = buffer.data(hl + 52);
    const auto *hl_53 = buffer.data(hl + 53);
    const auto *hl_54 = buffer.data(hl + 54);
    const auto *hl_55 = buffer.data(hl + 55);
    const auto *hl_56 = buffer.data(hl + 56);
    const auto *hl_57 = buffer.data(hl + 57);
    const auto *hl_58 = buffer.data(hl + 58);
    const auto *hl_59 = buffer.data(hl + 59);

    const auto *kl_90 = buffer.data(kl + 90);
    const auto *kl_91 = buffer.data(kl + 91);
    const auto *kl_92 = buffer.data(kl + 92);
    const auto *kl_93 = buffer.data(kl + 93);
    const auto *kl_94 = buffer.data(kl + 94);
    const auto *kl_95 = buffer.data(kl + 95);
    const auto *kl_96 = buffer.data(kl + 96);
    const auto *kl_97 = buffer.data(kl + 97);
    const auto *kl_98 = buffer.data(kl + 98);
    const auto *kl_99 = buffer.data(kl + 99);
    const auto *kl_100 = buffer.data(kl + 100);
    const auto *kl_101 = buffer.data(kl + 101);
    const auto *kl_102 = buffer.data(kl + 102);
    const auto *kl_103 = buffer.data(kl + 103);
    const auto *kl_104 = buffer.data(kl + 104);
    const auto *kl_105 = buffer.data(kl + 105);
    const auto *kl_106 = buffer.data(kl + 106);
    const auto *kl_107 = buffer.data(kl + 107);
    const auto *kl_108 = buffer.data(kl + 108);
    const auto *kl_109 = buffer.data(kl + 109);
    const auto *kl_110 = buffer.data(kl + 110);
    const auto *kl_111 = buffer.data(kl + 111);
    const auto *kl_112 = buffer.data(kl + 112);
    const auto *kl_113 = buffer.data(kl + 113);
    const auto *kl_114 = buffer.data(kl + 114);
    const auto *kl_115 = buffer.data(kl + 115);
    const auto *kl_116 = buffer.data(kl + 116);
    const auto *kl_117 = buffer.data(kl + 117);
    const auto *kl_118 = buffer.data(kl + 118);
    const auto *kl_119 = buffer.data(kl + 119);
    const auto *kl_120 = buffer.data(kl + 120);
    const auto *kl_121 = buffer.data(kl + 121);
    const auto *kl_122 = buffer.data(kl + 122);
    const auto *kl_123 = buffer.data(kl + 123);
    const auto *kl_124 = buffer.data(kl + 124);
    const auto *kl_125 = buffer.data(kl + 125);
    const auto *kl_126 = buffer.data(kl + 126);
    const auto *kl_127 = buffer.data(kl + 127);
    const auto *kl_128 = buffer.data(kl + 128);
    const auto *kl_129 = buffer.data(kl + 129);
    const auto *kl_130 = buffer.data(kl + 130);
    const auto *kl_131 = buffer.data(kl + 131);
    const auto *kl_132 = buffer.data(kl + 132);
    const auto *kl_133 = buffer.data(kl + 133);
    const auto *kl_134 = buffer.data(kl + 134);
    const auto *kl_180 = buffer.data(kl + 180);
    const auto *kl_181 = buffer.data(kl + 181);
    const auto *kl_182 = buffer.data(kl + 182);
    const auto *kl_183 = buffer.data(kl + 183);
    const auto *kl_184 = buffer.data(kl + 184);
    const auto *kl_185 = buffer.data(kl + 185);
    const auto *kl_186 = buffer.data(kl + 186);
    const auto *kl_187 = buffer.data(kl + 187);
    const auto *kl_188 = buffer.data(kl + 188);
    const auto *kl_189 = buffer.data(kl + 189);
    const auto *kl_190 = buffer.data(kl + 190);
    const auto *kl_191 = buffer.data(kl + 191);
    const auto *kl_192 = buffer.data(kl + 192);
    const auto *kl_193 = buffer.data(kl + 193);
    const auto *kl_194 = buffer.data(kl + 194);
    const auto *kl_195 = buffer.data(kl + 195);
    const auto *kl_196 = buffer.data(kl + 196);
    const auto *kl_197 = buffer.data(kl + 197);
    const auto *kl_198 = buffer.data(kl + 198);
    const auto *kl_199 = buffer.data(kl + 199);
    const auto *kl_200 = buffer.data(kl + 200);
    const auto *kl_201 = buffer.data(kl + 201);
    const auto *kl_202 = buffer.data(kl + 202);
    const auto *kl_203 = buffer.data(kl + 203);
    const auto *kl_204 = buffer.data(kl + 204);
    const auto *kl_205 = buffer.data(kl + 205);
    const auto *kl_206 = buffer.data(kl + 206);
    const auto *kl_207 = buffer.data(kl + 207);
    const auto *kl_208 = buffer.data(kl + 208);
    const auto *kl_209 = buffer.data(kl + 209);
    const auto *kl_210 = buffer.data(kl + 210);
    const auto *kl_211 = buffer.data(kl + 211);
    const auto *kl_212 = buffer.data(kl + 212);
    const auto *kl_213 = buffer.data(kl + 213);
    const auto *kl_214 = buffer.data(kl + 214);
    const auto *kl_215 = buffer.data(kl + 215);
    const auto *kl_216 = buffer.data(kl + 216);
    const auto *kl_217 = buffer.data(kl + 217);
    const auto *kl_218 = buffer.data(kl + 218);
    const auto *kl_219 = buffer.data(kl + 219);
    const auto *kl_220 = buffer.data(kl + 220);
    const auto *kl_221 = buffer.data(kl + 221);
    const auto *kl_222 = buffer.data(kl + 222);
    const auto *kl_223 = buffer.data(kl + 223);
    const auto *kl_224 = buffer.data(kl + 224);
    const auto *kl_225 = buffer.data(kl + 225);
    const auto *kl_226 = buffer.data(kl + 226);
    const auto *kl_227 = buffer.data(kl + 227);
    const auto *kl_228 = buffer.data(kl + 228);
    const auto *kl_229 = buffer.data(kl + 229);
    const auto *kl_230 = buffer.data(kl + 230);
    const auto *kl_231 = buffer.data(kl + 231);
    const auto *kl_232 = buffer.data(kl + 232);
    const auto *kl_233 = buffer.data(kl + 233);
    const auto *kl_234 = buffer.data(kl + 234);
    const auto *kl_235 = buffer.data(kl + 235);
    const auto *kl_236 = buffer.data(kl + 236);
    const auto *kl_237 = buffer.data(kl + 237);
    const auto *kl_238 = buffer.data(kl + 238);
    const auto *kl_239 = buffer.data(kl + 239);
    const auto *kl_240 = buffer.data(kl + 240);
    const auto *kl_241 = buffer.data(kl + 241);
    const auto *kl_242 = buffer.data(kl + 242);
    const auto *kl_243 = buffer.data(kl + 243);
    const auto *kl_244 = buffer.data(kl + 244);
    const auto *kl_245 = buffer.data(kl + 245);
    const auto *kl_246 = buffer.data(kl + 246);
    const auto *kl_247 = buffer.data(kl + 247);
    const auto *kl_248 = buffer.data(kl + 248);
    const auto *kl_249 = buffer.data(kl + 249);
    const auto *kl_250 = buffer.data(kl + 250);
    const auto *kl_251 = buffer.data(kl + 251);
    const auto *kl_252 = buffer.data(kl + 252);
    const auto *kl_253 = buffer.data(kl + 253);
    const auto *kl_254 = buffer.data(kl + 254);
    const auto *kl_255 = buffer.data(kl + 255);
    const auto *kl_256 = buffer.data(kl + 256);
    const auto *kl_257 = buffer.data(kl + 257);
    const auto *kl_258 = buffer.data(kl + 258);
    const auto *kl_259 = buffer.data(kl + 259);
    const auto *kl_260 = buffer.data(kl + 260);
    const auto *kl_261 = buffer.data(kl + 261);
    const auto *kl_262 = buffer.data(kl + 262);
    const auto *kl_263 = buffer.data(kl + 263);
    const auto *kl_264 = buffer.data(kl + 264);
    const auto *kl_265 = buffer.data(kl + 265);
    const auto *kl_266 = buffer.data(kl + 266);
    const auto *kl_267 = buffer.data(kl + 267);
    const auto *kl_268 = buffer.data(kl + 268);
    const auto *kl_269 = buffer.data(kl + 269);
    const auto *kl_315 = buffer.data(kl + 315);
    const auto *kl_316 = buffer.data(kl + 316);
    const auto *kl_317 = buffer.data(kl + 317);
    const auto *kl_318 = buffer.data(kl + 318);
    const auto *kl_319 = buffer.data(kl + 319);
    const auto *kl_320 = buffer.data(kl + 320);
    const auto *kl_321 = buffer.data(kl + 321);
    const auto *kl_322 = buffer.data(kl + 322);
    const auto *kl_323 = buffer.data(kl + 323);
    const auto *kl_324 = buffer.data(kl + 324);
    const auto *kl_325 = buffer.data(kl + 325);
    const auto *kl_326 = buffer.data(kl + 326);
    const auto *kl_327 = buffer.data(kl + 327);
    const auto *kl_328 = buffer.data(kl + 328);
    const auto *kl_329 = buffer.data(kl + 329);
    const auto *kl_330 = buffer.data(kl + 330);
    const auto *kl_331 = buffer.data(kl + 331);
    const auto *kl_332 = buffer.data(kl + 332);
    const auto *kl_333 = buffer.data(kl + 333);
    const auto *kl_334 = buffer.data(kl + 334);
    const auto *kl_335 = buffer.data(kl + 335);
    const auto *kl_336 = buffer.data(kl + 336);
    const auto *kl_337 = buffer.data(kl + 337);
    const auto *kl_338 = buffer.data(kl + 338);
    const auto *kl_339 = buffer.data(kl + 339);
    const auto *kl_340 = buffer.data(kl + 340);
    const auto *kl_341 = buffer.data(kl + 341);
    const auto *kl_342 = buffer.data(kl + 342);
    const auto *kl_343 = buffer.data(kl + 343);
    const auto *kl_344 = buffer.data(kl + 344);
    const auto *kl_345 = buffer.data(kl + 345);
    const auto *kl_346 = buffer.data(kl + 346);
    const auto *kl_347 = buffer.data(kl + 347);
    const auto *kl_348 = buffer.data(kl + 348);
    const auto *kl_349 = buffer.data(kl + 349);
    const auto *kl_350 = buffer.data(kl + 350);
    const auto *kl_351 = buffer.data(kl + 351);
    const auto *kl_352 = buffer.data(kl + 352);
    const auto *kl_353 = buffer.data(kl + 353);
    const auto *kl_354 = buffer.data(kl + 354);
    const auto *kl_355 = buffer.data(kl + 355);
    const auto *kl_356 = buffer.data(kl + 356);
    const auto *kl_357 = buffer.data(kl + 357);
    const auto *kl_358 = buffer.data(kl + 358);
    const auto *kl_359 = buffer.data(kl + 359);
    const auto *kl_360 = buffer.data(kl + 360);
    const auto *kl_361 = buffer.data(kl + 361);
    const auto *kl_362 = buffer.data(kl + 362);
    const auto *kl_363 = buffer.data(kl + 363);
    const auto *kl_364 = buffer.data(kl + 364);
    const auto *kl_365 = buffer.data(kl + 365);
    const auto *kl_366 = buffer.data(kl + 366);
    const auto *kl_367 = buffer.data(kl + 367);
    const auto *kl_368 = buffer.data(kl + 368);
    const auto *kl_369 = buffer.data(kl + 369);
    const auto *kl_370 = buffer.data(kl + 370);
    const auto *kl_371 = buffer.data(kl + 371);
    const auto *kl_372 = buffer.data(kl + 372);
    const auto *kl_373 = buffer.data(kl + 373);
    const auto *kl_374 = buffer.data(kl + 374);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, kl_90, kl_91, kl_92, kl_93, \
                         kl_94, kl_95, kl_96, kl_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kl_90[k];

        t_1[k] = f_0 * kl_91[k];

        t_2[k] = f_0 * kl_92[k];

        t_3[k] = f_0 * kl_93[k];

        t_4[k] = f_0 * kl_94[k];

        t_5[k] = f_0 * kl_95[k];

        t_6[k] = f_0 * kl_96[k];

        t_7[k] = f_0 * kl_97[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, kl_98, kl_99, kl_100, \
                         kl_101, kl_102, kl_103, kl_104, kl_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * kl_98[k];

        t_9[k] = f_0 * kl_99[k];

        t_10[k] = f_0 * kl_100[k];

        t_11[k] = f_0 * kl_101[k];

        t_12[k] = f_0 * kl_102[k];

        t_13[k] = f_0 * kl_103[k];

        t_14[k] = f_0 * kl_104[k];

        t_15[k] = f_0 * kl_105[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, kl_106, kl_107, \
                         kl_108, kl_109, kl_110, kl_111, kl_112, \
                         kl_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * kl_106[k];

        t_17[k] = f_0 * kl_107[k];

        t_18[k] = f_0 * kl_108[k];

        t_19[k] = f_0 * kl_109[k];

        t_20[k] = f_0 * kl_110[k];

        t_21[k] = f_0 * kl_111[k];

        t_22[k] = f_0 * kl_112[k];

        t_23[k] = f_0 * kl_113[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, kl_114, kl_115, \
                         kl_116, kl_117, kl_118, kl_119, kl_120, \
                         kl_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * kl_114[k];

        t_25[k] = f_0 * kl_115[k];

        t_26[k] = f_0 * kl_116[k];

        t_27[k] = f_0 * kl_117[k];

        t_28[k] = f_0 * kl_118[k];

        t_29[k] = f_0 * kl_119[k];

        t_30[k] = f_0 * kl_120[k];

        t_31[k] = f_0 * kl_121[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, kl_122, kl_123, \
                         kl_124, kl_125, kl_126, kl_127, kl_128, \
                         kl_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * kl_122[k];

        t_33[k] = f_0 * kl_123[k];

        t_34[k] = f_0 * kl_124[k];

        t_35[k] = f_0 * kl_125[k];

        t_36[k] = f_0 * kl_126[k];

        t_37[k] = f_0 * kl_127[k];

        t_38[k] = f_0 * kl_128[k];

        t_39[k] = f_0 * kl_129[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, kl_130, kl_131, \
                         kl_132, kl_133, kl_134, kl_180, kl_181, \
                         kl_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * kl_130[k];

        t_41[k] = f_0 * kl_131[k];

        t_42[k] = f_0 * kl_132[k];

        t_43[k] = f_0 * kl_133[k];

        t_44[k] = f_0 * kl_134[k];

        t_45[k] = f_0 * kl_180[k];

        t_46[k] = f_0 * kl_181[k];

        t_47[k] = f_0 * kl_182[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, kl_183, kl_184, \
                         kl_185, kl_186, kl_187, kl_188, kl_189, \
                         kl_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * kl_183[k];

        t_49[k] = f_0 * kl_184[k];

        t_50[k] = f_0 * kl_185[k];

        t_51[k] = f_0 * kl_186[k];

        t_52[k] = f_0 * kl_187[k];

        t_53[k] = f_0 * kl_188[k];

        t_54[k] = f_0 * kl_189[k];

        t_55[k] = f_0 * kl_190[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, kl_191, kl_192, \
                         kl_193, kl_194, kl_195, kl_196, kl_197, \
                         kl_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * kl_191[k];

        t_57[k] = f_0 * kl_192[k];

        t_58[k] = f_0 * kl_193[k];

        t_59[k] = f_0 * kl_194[k];

        t_60[k] = f_0 * kl_195[k];

        t_61[k] = f_0 * kl_196[k];

        t_62[k] = f_0 * kl_197[k];

        t_63[k] = f_0 * kl_198[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, kl_199, kl_200, \
                         kl_201, kl_202, kl_203, kl_204, kl_205, \
                         kl_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * kl_199[k];

        t_65[k] = f_0 * kl_200[k];

        t_66[k] = f_0 * kl_201[k];

        t_67[k] = f_0 * kl_202[k];

        t_68[k] = f_0 * kl_203[k];

        t_69[k] = f_0 * kl_204[k];

        t_70[k] = f_0 * kl_205[k];

        t_71[k] = f_0 * kl_206[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, t_78, t_79, kl_207, kl_208, \
                         kl_209, kl_210, kl_211, kl_212, kl_213, \
                         kl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * kl_207[k];

        t_73[k] = f_0 * kl_208[k];

        t_74[k] = f_0 * kl_209[k];

        t_75[k] = f_0 * kl_210[k];

        t_76[k] = f_0 * kl_211[k];

        t_77[k] = f_0 * kl_212[k];

        t_78[k] = f_0 * kl_213[k];

        t_79[k] = f_0 * kl_214[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, t_86, t_87, kl_215, kl_216, \
                         kl_217, kl_218, kl_219, kl_220, kl_221, \
                         kl_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * kl_215[k];

        t_81[k] = f_0 * kl_216[k];

        t_82[k] = f_0 * kl_217[k];

        t_83[k] = f_0 * kl_218[k];

        t_84[k] = f_0 * kl_219[k];

        t_85[k] = f_0 * kl_220[k];

        t_86[k] = f_0 * kl_221[k];

        t_87[k] = f_0 * kl_222[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, t_93, hl_0, hl_1, hl_2, hl_3, kl_223, \
                         kl_224, kl_225, kl_226, kl_227, kl_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_0 * kl_223[k];

        t_89[k] = f_0 * kl_224[k];

        t_90[k] = -hl_0[k]
                  + f_0 * kl_225[k];

        t_91[k] = -hl_1[k]
                  + f_0 * kl_226[k];

        t_92[k] = -hl_2[k]
                  + f_0 * kl_227[k];

        t_93[k] = -hl_3[k]
                  + f_0 * kl_228[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, hl_4, hl_5, hl_6, hl_7, hl_8, kl_229, \
                         kl_230, kl_231, kl_232, kl_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -hl_4[k]
                  + f_0 * kl_229[k];

        t_95[k] = -hl_5[k]
                  + f_0 * kl_230[k];

        t_96[k] = -hl_6[k]
                  + f_0 * kl_231[k];

        t_97[k] = -hl_7[k]
                  + f_0 * kl_232[k];

        t_98[k] = -hl_8[k]
                  + f_0 * kl_233[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, hl_9, hl_10, hl_11, hl_12, hl_13, \
                         kl_234, kl_235, kl_236, kl_237, kl_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -hl_9[k]
                  + f_0 * kl_234[k];

        t_100[k] = -hl_10[k]
                   + f_0 * kl_235[k];

        t_101[k] = -hl_11[k]
                   + f_0 * kl_236[k];

        t_102[k] = -hl_12[k]
                   + f_0 * kl_237[k];

        t_103[k] = -hl_13[k]
                   + f_0 * kl_238[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, hl_14, hl_15, hl_16, hl_17, hl_18, \
                         kl_239, kl_240, kl_241, kl_242, kl_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -hl_14[k]
                   + f_0 * kl_239[k];

        t_105[k] = -hl_15[k]
                   + f_0 * kl_240[k];

        t_106[k] = -hl_16[k]
                   + f_0 * kl_241[k];

        t_107[k] = -hl_17[k]
                   + f_0 * kl_242[k];

        t_108[k] = -hl_18[k]
                   + f_0 * kl_243[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, hl_19, hl_20, hl_21, hl_22, hl_23, \
                         kl_244, kl_245, kl_246, kl_247, kl_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -hl_19[k]
                   + f_0 * kl_244[k];

        t_110[k] = -hl_20[k]
                   + f_0 * kl_245[k];

        t_111[k] = -hl_21[k]
                   + f_0 * kl_246[k];

        t_112[k] = -hl_22[k]
                   + f_0 * kl_247[k];

        t_113[k] = -hl_23[k]
                   + f_0 * kl_248[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, hl_24, hl_25, hl_26, hl_27, hl_28, \
                         kl_249, kl_250, kl_251, kl_252, kl_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -hl_24[k]
                   + f_0 * kl_249[k];

        t_115[k] = -hl_25[k]
                   + f_0 * kl_250[k];

        t_116[k] = -hl_26[k]
                   + f_0 * kl_251[k];

        t_117[k] = -hl_27[k]
                   + f_0 * kl_252[k];

        t_118[k] = -hl_28[k]
                   + f_0 * kl_253[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, hl_29, hl_30, hl_31, hl_32, hl_33, \
                         kl_254, kl_255, kl_256, kl_257, kl_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -hl_29[k]
                   + f_0 * kl_254[k];

        t_120[k] = -hl_30[k]
                   + f_0 * kl_255[k];

        t_121[k] = -hl_31[k]
                   + f_0 * kl_256[k];

        t_122[k] = -hl_32[k]
                   + f_0 * kl_257[k];

        t_123[k] = -hl_33[k]
                   + f_0 * kl_258[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, hl_34, hl_35, hl_36, hl_37, hl_38, \
                         kl_259, kl_260, kl_261, kl_262, kl_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -hl_34[k]
                   + f_0 * kl_259[k];

        t_125[k] = -hl_35[k]
                   + f_0 * kl_260[k];

        t_126[k] = -hl_36[k]
                   + f_0 * kl_261[k];

        t_127[k] = -hl_37[k]
                   + f_0 * kl_262[k];

        t_128[k] = -hl_38[k]
                   + f_0 * kl_263[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, hl_39, hl_40, hl_41, hl_42, hl_43, \
                         kl_264, kl_265, kl_266, kl_267, kl_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -hl_39[k]
                   + f_0 * kl_264[k];

        t_130[k] = -hl_40[k]
                   + f_0 * kl_265[k];

        t_131[k] = -hl_41[k]
                   + f_0 * kl_266[k];

        t_132[k] = -hl_42[k]
                   + f_0 * kl_267[k];

        t_133[k] = -hl_43[k]
                   + f_0 * kl_268[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, t_139, t_140, hl_44, kl_269, \
                         kl_315, kl_316, kl_317, kl_318, kl_319, \
                         kl_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -hl_44[k]
                   + f_0 * kl_269[k];

        t_135[k] = f_0 * kl_315[k];

        t_136[k] = f_0 * kl_316[k];

        t_137[k] = f_0 * kl_317[k];

        t_138[k] = f_0 * kl_318[k];

        t_139[k] = f_0 * kl_319[k];

        t_140[k] = f_0 * kl_320[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, t_146, t_147, t_148, kl_321, \
                         kl_322, kl_323, kl_324, kl_325, kl_326, kl_327, \
                         kl_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_0 * kl_321[k];

        t_142[k] = f_0 * kl_322[k];

        t_143[k] = f_0 * kl_323[k];

        t_144[k] = f_0 * kl_324[k];

        t_145[k] = f_0 * kl_325[k];

        t_146[k] = f_0 * kl_326[k];

        t_147[k] = f_0 * kl_327[k];

        t_148[k] = f_0 * kl_328[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, t_154, t_155, t_156, kl_329, \
                         kl_330, kl_331, kl_332, kl_333, kl_334, kl_335, \
                         kl_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * kl_329[k];

        t_150[k] = f_0 * kl_330[k];

        t_151[k] = f_0 * kl_331[k];

        t_152[k] = f_0 * kl_332[k];

        t_153[k] = f_0 * kl_333[k];

        t_154[k] = f_0 * kl_334[k];

        t_155[k] = f_0 * kl_335[k];

        t_156[k] = f_0 * kl_336[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, t_162, t_163, t_164, kl_337, \
                         kl_338, kl_339, kl_340, kl_341, kl_342, kl_343, \
                         kl_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_0 * kl_337[k];

        t_158[k] = f_0 * kl_338[k];

        t_159[k] = f_0 * kl_339[k];

        t_160[k] = f_0 * kl_340[k];

        t_161[k] = f_0 * kl_341[k];

        t_162[k] = f_0 * kl_342[k];

        t_163[k] = f_0 * kl_343[k];

        t_164[k] = f_0 * kl_344[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, t_170, t_171, t_172, kl_345, \
                         kl_346, kl_347, kl_348, kl_349, kl_350, kl_351, \
                         kl_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_0 * kl_345[k];

        t_166[k] = f_0 * kl_346[k];

        t_167[k] = f_0 * kl_347[k];

        t_168[k] = f_0 * kl_348[k];

        t_169[k] = f_0 * kl_349[k];

        t_170[k] = f_0 * kl_350[k];

        t_171[k] = f_0 * kl_351[k];

        t_172[k] = f_0 * kl_352[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, t_179, kl_353, kl_354, \
                         kl_355, kl_356, kl_357, kl_358, kl_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_0 * kl_353[k];

        t_174[k] = f_0 * kl_354[k];

        t_175[k] = f_0 * kl_355[k];

        t_176[k] = f_0 * kl_356[k];

        t_177[k] = f_0 * kl_357[k];

        t_178[k] = f_0 * kl_358[k];

        t_179[k] = f_0 * kl_359[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, hl_45, hl_46, hl_47, hl_48, hl_49, \
                         kl_360, kl_361, kl_362, kl_363, kl_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -hl_45[k]
                   + f_0 * kl_360[k];

        t_181[k] = -hl_46[k]
                   + f_0 * kl_361[k];

        t_182[k] = -hl_47[k]
                   + f_0 * kl_362[k];

        t_183[k] = -hl_48[k]
                   + f_0 * kl_363[k];

        t_184[k] = -hl_49[k]
                   + f_0 * kl_364[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, hl_50, hl_51, hl_52, hl_53, hl_54, \
                         kl_365, kl_366, kl_367, kl_368, kl_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -hl_50[k]
                   + f_0 * kl_365[k];

        t_186[k] = -hl_51[k]
                   + f_0 * kl_366[k];

        t_187[k] = -hl_52[k]
                   + f_0 * kl_367[k];

        t_188[k] = -hl_53[k]
                   + f_0 * kl_368[k];

        t_189[k] = -hl_54[k]
                   + f_0 * kl_369[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, hl_55, hl_56, hl_57, hl_58, hl_59, \
                         kl_370, kl_371, kl_372, kl_373, kl_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -hl_55[k]
                   + f_0 * kl_370[k];

        t_191[k] = -hl_56[k]
                   + f_0 * kl_371[k];

        t_192[k] = -hl_57[k]
                   + f_0 * kl_372[k];

        t_193[k] = -hl_58[k]
                   + f_0 * kl_373[k];

        t_194[k] = -hl_59[k]
                   + f_0 * kl_374[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
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
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);

    const auto *hl_60 = buffer.data(hl + 60);
    const auto *hl_61 = buffer.data(hl + 61);
    const auto *hl_62 = buffer.data(hl + 62);
    const auto *hl_63 = buffer.data(hl + 63);
    const auto *hl_64 = buffer.data(hl + 64);
    const auto *hl_65 = buffer.data(hl + 65);
    const auto *hl_66 = buffer.data(hl + 66);
    const auto *hl_67 = buffer.data(hl + 67);
    const auto *hl_68 = buffer.data(hl + 68);
    const auto *hl_69 = buffer.data(hl + 69);
    const auto *hl_70 = buffer.data(hl + 70);
    const auto *hl_71 = buffer.data(hl + 71);
    const auto *hl_72 = buffer.data(hl + 72);
    const auto *hl_73 = buffer.data(hl + 73);
    const auto *hl_74 = buffer.data(hl + 74);
    const auto *hl_75 = buffer.data(hl + 75);
    const auto *hl_76 = buffer.data(hl + 76);
    const auto *hl_77 = buffer.data(hl + 77);
    const auto *hl_78 = buffer.data(hl + 78);
    const auto *hl_79 = buffer.data(hl + 79);
    const auto *hl_80 = buffer.data(hl + 80);
    const auto *hl_81 = buffer.data(hl + 81);
    const auto *hl_82 = buffer.data(hl + 82);
    const auto *hl_83 = buffer.data(hl + 83);
    const auto *hl_84 = buffer.data(hl + 84);
    const auto *hl_85 = buffer.data(hl + 85);
    const auto *hl_86 = buffer.data(hl + 86);
    const auto *hl_87 = buffer.data(hl + 87);
    const auto *hl_88 = buffer.data(hl + 88);
    const auto *hl_89 = buffer.data(hl + 89);
    const auto *hl_90 = buffer.data(hl + 90);
    const auto *hl_91 = buffer.data(hl + 91);
    const auto *hl_92 = buffer.data(hl + 92);
    const auto *hl_93 = buffer.data(hl + 93);
    const auto *hl_94 = buffer.data(hl + 94);
    const auto *hl_95 = buffer.data(hl + 95);
    const auto *hl_96 = buffer.data(hl + 96);
    const auto *hl_97 = buffer.data(hl + 97);
    const auto *hl_98 = buffer.data(hl + 98);
    const auto *hl_99 = buffer.data(hl + 99);
    const auto *hl_100 = buffer.data(hl + 100);
    const auto *hl_101 = buffer.data(hl + 101);
    const auto *hl_102 = buffer.data(hl + 102);
    const auto *hl_103 = buffer.data(hl + 103);
    const auto *hl_104 = buffer.data(hl + 104);
    const auto *hl_105 = buffer.data(hl + 105);
    const auto *hl_106 = buffer.data(hl + 106);
    const auto *hl_107 = buffer.data(hl + 107);
    const auto *hl_108 = buffer.data(hl + 108);
    const auto *hl_109 = buffer.data(hl + 109);
    const auto *hl_110 = buffer.data(hl + 110);
    const auto *hl_111 = buffer.data(hl + 111);
    const auto *hl_112 = buffer.data(hl + 112);
    const auto *hl_113 = buffer.data(hl + 113);
    const auto *hl_114 = buffer.data(hl + 114);
    const auto *hl_115 = buffer.data(hl + 115);
    const auto *hl_116 = buffer.data(hl + 116);
    const auto *hl_117 = buffer.data(hl + 117);
    const auto *hl_118 = buffer.data(hl + 118);
    const auto *hl_119 = buffer.data(hl + 119);
    const auto *hl_120 = buffer.data(hl + 120);
    const auto *hl_121 = buffer.data(hl + 121);
    const auto *hl_122 = buffer.data(hl + 122);
    const auto *hl_123 = buffer.data(hl + 123);
    const auto *hl_124 = buffer.data(hl + 124);
    const auto *hl_125 = buffer.data(hl + 125);
    const auto *hl_126 = buffer.data(hl + 126);
    const auto *hl_127 = buffer.data(hl + 127);
    const auto *hl_128 = buffer.data(hl + 128);
    const auto *hl_129 = buffer.data(hl + 129);
    const auto *hl_130 = buffer.data(hl + 130);
    const auto *hl_131 = buffer.data(hl + 131);
    const auto *hl_132 = buffer.data(hl + 132);
    const auto *hl_133 = buffer.data(hl + 133);
    const auto *hl_134 = buffer.data(hl + 134);
    const auto *hl_135 = buffer.data(hl + 135);
    const auto *hl_136 = buffer.data(hl + 136);
    const auto *hl_137 = buffer.data(hl + 137);
    const auto *hl_138 = buffer.data(hl + 138);
    const auto *hl_139 = buffer.data(hl + 139);
    const auto *hl_140 = buffer.data(hl + 140);
    const auto *hl_141 = buffer.data(hl + 141);
    const auto *hl_142 = buffer.data(hl + 142);
    const auto *hl_143 = buffer.data(hl + 143);
    const auto *hl_144 = buffer.data(hl + 144);
    const auto *hl_145 = buffer.data(hl + 145);
    const auto *hl_146 = buffer.data(hl + 146);
    const auto *hl_147 = buffer.data(hl + 147);
    const auto *hl_148 = buffer.data(hl + 148);
    const auto *hl_149 = buffer.data(hl + 149);
    const auto *hl_150 = buffer.data(hl + 150);
    const auto *hl_151 = buffer.data(hl + 151);
    const auto *hl_152 = buffer.data(hl + 152);
    const auto *hl_153 = buffer.data(hl + 153);
    const auto *hl_154 = buffer.data(hl + 154);
    const auto *hl_155 = buffer.data(hl + 155);
    const auto *hl_156 = buffer.data(hl + 156);
    const auto *hl_157 = buffer.data(hl + 157);
    const auto *hl_158 = buffer.data(hl + 158);
    const auto *hl_159 = buffer.data(hl + 159);
    const auto *hl_160 = buffer.data(hl + 160);
    const auto *hl_161 = buffer.data(hl + 161);
    const auto *hl_162 = buffer.data(hl + 162);
    const auto *hl_163 = buffer.data(hl + 163);
    const auto *hl_164 = buffer.data(hl + 164);
    const auto *hl_165 = buffer.data(hl + 165);
    const auto *hl_166 = buffer.data(hl + 166);
    const auto *hl_167 = buffer.data(hl + 167);
    const auto *hl_168 = buffer.data(hl + 168);
    const auto *hl_169 = buffer.data(hl + 169);
    const auto *hl_170 = buffer.data(hl + 170);
    const auto *hl_171 = buffer.data(hl + 171);
    const auto *hl_172 = buffer.data(hl + 172);
    const auto *hl_173 = buffer.data(hl + 173);
    const auto *hl_174 = buffer.data(hl + 174);
    const auto *hl_175 = buffer.data(hl + 175);
    const auto *hl_176 = buffer.data(hl + 176);

    const auto *kl_375 = buffer.data(kl + 375);
    const auto *kl_376 = buffer.data(kl + 376);
    const auto *kl_377 = buffer.data(kl + 377);
    const auto *kl_378 = buffer.data(kl + 378);
    const auto *kl_379 = buffer.data(kl + 379);
    const auto *kl_380 = buffer.data(kl + 380);
    const auto *kl_381 = buffer.data(kl + 381);
    const auto *kl_382 = buffer.data(kl + 382);
    const auto *kl_383 = buffer.data(kl + 383);
    const auto *kl_384 = buffer.data(kl + 384);
    const auto *kl_385 = buffer.data(kl + 385);
    const auto *kl_386 = buffer.data(kl + 386);
    const auto *kl_387 = buffer.data(kl + 387);
    const auto *kl_388 = buffer.data(kl + 388);
    const auto *kl_389 = buffer.data(kl + 389);
    const auto *kl_390 = buffer.data(kl + 390);
    const auto *kl_391 = buffer.data(kl + 391);
    const auto *kl_392 = buffer.data(kl + 392);
    const auto *kl_393 = buffer.data(kl + 393);
    const auto *kl_394 = buffer.data(kl + 394);
    const auto *kl_395 = buffer.data(kl + 395);
    const auto *kl_396 = buffer.data(kl + 396);
    const auto *kl_397 = buffer.data(kl + 397);
    const auto *kl_398 = buffer.data(kl + 398);
    const auto *kl_399 = buffer.data(kl + 399);
    const auto *kl_400 = buffer.data(kl + 400);
    const auto *kl_401 = buffer.data(kl + 401);
    const auto *kl_402 = buffer.data(kl + 402);
    const auto *kl_403 = buffer.data(kl + 403);
    const auto *kl_404 = buffer.data(kl + 404);
    const auto *kl_405 = buffer.data(kl + 405);
    const auto *kl_406 = buffer.data(kl + 406);
    const auto *kl_407 = buffer.data(kl + 407);
    const auto *kl_408 = buffer.data(kl + 408);
    const auto *kl_409 = buffer.data(kl + 409);
    const auto *kl_410 = buffer.data(kl + 410);
    const auto *kl_411 = buffer.data(kl + 411);
    const auto *kl_412 = buffer.data(kl + 412);
    const auto *kl_413 = buffer.data(kl + 413);
    const auto *kl_414 = buffer.data(kl + 414);
    const auto *kl_415 = buffer.data(kl + 415);
    const auto *kl_416 = buffer.data(kl + 416);
    const auto *kl_417 = buffer.data(kl + 417);
    const auto *kl_418 = buffer.data(kl + 418);
    const auto *kl_419 = buffer.data(kl + 419);
    const auto *kl_420 = buffer.data(kl + 420);
    const auto *kl_421 = buffer.data(kl + 421);
    const auto *kl_422 = buffer.data(kl + 422);
    const auto *kl_423 = buffer.data(kl + 423);
    const auto *kl_424 = buffer.data(kl + 424);
    const auto *kl_425 = buffer.data(kl + 425);
    const auto *kl_426 = buffer.data(kl + 426);
    const auto *kl_427 = buffer.data(kl + 427);
    const auto *kl_428 = buffer.data(kl + 428);
    const auto *kl_429 = buffer.data(kl + 429);
    const auto *kl_430 = buffer.data(kl + 430);
    const auto *kl_431 = buffer.data(kl + 431);
    const auto *kl_432 = buffer.data(kl + 432);
    const auto *kl_433 = buffer.data(kl + 433);
    const auto *kl_434 = buffer.data(kl + 434);
    const auto *kl_435 = buffer.data(kl + 435);
    const auto *kl_436 = buffer.data(kl + 436);
    const auto *kl_437 = buffer.data(kl + 437);
    const auto *kl_438 = buffer.data(kl + 438);
    const auto *kl_439 = buffer.data(kl + 439);
    const auto *kl_440 = buffer.data(kl + 440);
    const auto *kl_441 = buffer.data(kl + 441);
    const auto *kl_442 = buffer.data(kl + 442);
    const auto *kl_443 = buffer.data(kl + 443);
    const auto *kl_444 = buffer.data(kl + 444);
    const auto *kl_445 = buffer.data(kl + 445);
    const auto *kl_446 = buffer.data(kl + 446);
    const auto *kl_447 = buffer.data(kl + 447);
    const auto *kl_448 = buffer.data(kl + 448);
    const auto *kl_449 = buffer.data(kl + 449);
    const auto *kl_495 = buffer.data(kl + 495);
    const auto *kl_496 = buffer.data(kl + 496);
    const auto *kl_497 = buffer.data(kl + 497);
    const auto *kl_498 = buffer.data(kl + 498);
    const auto *kl_499 = buffer.data(kl + 499);
    const auto *kl_500 = buffer.data(kl + 500);
    const auto *kl_501 = buffer.data(kl + 501);
    const auto *kl_502 = buffer.data(kl + 502);
    const auto *kl_503 = buffer.data(kl + 503);
    const auto *kl_504 = buffer.data(kl + 504);
    const auto *kl_505 = buffer.data(kl + 505);
    const auto *kl_506 = buffer.data(kl + 506);
    const auto *kl_507 = buffer.data(kl + 507);
    const auto *kl_508 = buffer.data(kl + 508);
    const auto *kl_509 = buffer.data(kl + 509);
    const auto *kl_510 = buffer.data(kl + 510);
    const auto *kl_511 = buffer.data(kl + 511);
    const auto *kl_512 = buffer.data(kl + 512);
    const auto *kl_513 = buffer.data(kl + 513);
    const auto *kl_514 = buffer.data(kl + 514);
    const auto *kl_515 = buffer.data(kl + 515);
    const auto *kl_516 = buffer.data(kl + 516);
    const auto *kl_517 = buffer.data(kl + 517);
    const auto *kl_518 = buffer.data(kl + 518);
    const auto *kl_519 = buffer.data(kl + 519);
    const auto *kl_520 = buffer.data(kl + 520);
    const auto *kl_521 = buffer.data(kl + 521);
    const auto *kl_522 = buffer.data(kl + 522);
    const auto *kl_523 = buffer.data(kl + 523);
    const auto *kl_524 = buffer.data(kl + 524);
    const auto *kl_525 = buffer.data(kl + 525);
    const auto *kl_526 = buffer.data(kl + 526);
    const auto *kl_527 = buffer.data(kl + 527);
    const auto *kl_528 = buffer.data(kl + 528);
    const auto *kl_529 = buffer.data(kl + 529);
    const auto *kl_530 = buffer.data(kl + 530);
    const auto *kl_531 = buffer.data(kl + 531);
    const auto *kl_532 = buffer.data(kl + 532);
    const auto *kl_533 = buffer.data(kl + 533);
    const auto *kl_534 = buffer.data(kl + 534);
    const auto *kl_535 = buffer.data(kl + 535);
    const auto *kl_536 = buffer.data(kl + 536);
    const auto *kl_537 = buffer.data(kl + 537);
    const auto *kl_538 = buffer.data(kl + 538);
    const auto *kl_539 = buffer.data(kl + 539);
    const auto *kl_540 = buffer.data(kl + 540);
    const auto *kl_541 = buffer.data(kl + 541);
    const auto *kl_542 = buffer.data(kl + 542);
    const auto *kl_543 = buffer.data(kl + 543);
    const auto *kl_544 = buffer.data(kl + 544);
    const auto *kl_545 = buffer.data(kl + 545);
    const auto *kl_546 = buffer.data(kl + 546);
    const auto *kl_547 = buffer.data(kl + 547);
    const auto *kl_548 = buffer.data(kl + 548);
    const auto *kl_549 = buffer.data(kl + 549);
    const auto *kl_550 = buffer.data(kl + 550);
    const auto *kl_551 = buffer.data(kl + 551);
    const auto *kl_552 = buffer.data(kl + 552);
    const auto *kl_553 = buffer.data(kl + 553);
    const auto *kl_554 = buffer.data(kl + 554);
    const auto *kl_555 = buffer.data(kl + 555);
    const auto *kl_556 = buffer.data(kl + 556);
    const auto *kl_557 = buffer.data(kl + 557);
    const auto *kl_558 = buffer.data(kl + 558);
    const auto *kl_559 = buffer.data(kl + 559);
    const auto *kl_560 = buffer.data(kl + 560);
    const auto *kl_561 = buffer.data(kl + 561);
    const auto *kl_562 = buffer.data(kl + 562);
    const auto *kl_563 = buffer.data(kl + 563);
    const auto *kl_564 = buffer.data(kl + 564);
    const auto *kl_565 = buffer.data(kl + 565);
    const auto *kl_566 = buffer.data(kl + 566);
    const auto *kl_567 = buffer.data(kl + 567);
    const auto *kl_568 = buffer.data(kl + 568);
    const auto *kl_569 = buffer.data(kl + 569);
    const auto *kl_570 = buffer.data(kl + 570);
    const auto *kl_571 = buffer.data(kl + 571);
    const auto *kl_572 = buffer.data(kl + 572);
    const auto *kl_573 = buffer.data(kl + 573);
    const auto *kl_574 = buffer.data(kl + 574);
    const auto *kl_575 = buffer.data(kl + 575);
    const auto *kl_576 = buffer.data(kl + 576);
    const auto *kl_577 = buffer.data(kl + 577);
    const auto *kl_578 = buffer.data(kl + 578);
    const auto *kl_579 = buffer.data(kl + 579);
    const auto *kl_580 = buffer.data(kl + 580);
    const auto *kl_581 = buffer.data(kl + 581);

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, hl_60, hl_61, hl_62, hl_63, hl_64, \
                         kl_375, kl_376, kl_377, kl_378, kl_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -hl_60[k]
                   + f_0 * kl_375[k];

        t_196[k] = -hl_61[k]
                   + f_0 * kl_376[k];

        t_197[k] = -hl_62[k]
                   + f_0 * kl_377[k];

        t_198[k] = -hl_63[k]
                   + f_0 * kl_378[k];

        t_199[k] = -hl_64[k]
                   + f_0 * kl_379[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, hl_65, hl_66, hl_67, hl_68, hl_69, \
                         kl_380, kl_381, kl_382, kl_383, kl_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -hl_65[k]
                   + f_0 * kl_380[k];

        t_201[k] = -hl_66[k]
                   + f_0 * kl_381[k];

        t_202[k] = -hl_67[k]
                   + f_0 * kl_382[k];

        t_203[k] = -hl_68[k]
                   + f_0 * kl_383[k];

        t_204[k] = -hl_69[k]
                   + f_0 * kl_384[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, hl_70, hl_71, hl_72, hl_73, hl_74, \
                         kl_385, kl_386, kl_387, kl_388, kl_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -hl_70[k]
                   + f_0 * kl_385[k];

        t_206[k] = -hl_71[k]
                   + f_0 * kl_386[k];

        t_207[k] = -hl_72[k]
                   + f_0 * kl_387[k];

        t_208[k] = -hl_73[k]
                   + f_0 * kl_388[k];

        t_209[k] = -hl_74[k]
                   + f_0 * kl_389[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, hl_75, hl_76, hl_77, hl_78, hl_79, \
                         kl_390, kl_391, kl_392, kl_393, kl_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -hl_75[k]
                   + f_0 * kl_390[k];

        t_211[k] = -hl_76[k]
                   + f_0 * kl_391[k];

        t_212[k] = -hl_77[k]
                   + f_0 * kl_392[k];

        t_213[k] = -hl_78[k]
                   + f_0 * kl_393[k];

        t_214[k] = -hl_79[k]
                   + f_0 * kl_394[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, hl_80, hl_81, hl_82, hl_83, hl_84, \
                         kl_395, kl_396, kl_397, kl_398, kl_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -hl_80[k]
                   + f_0 * kl_395[k];

        t_216[k] = -hl_81[k]
                   + f_0 * kl_396[k];

        t_217[k] = -hl_82[k]
                   + f_0 * kl_397[k];

        t_218[k] = -hl_83[k]
                   + f_0 * kl_398[k];

        t_219[k] = -hl_84[k]
                   + f_0 * kl_399[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, hl_85, hl_86, hl_87, hl_88, hl_89, \
                         kl_400, kl_401, kl_402, kl_403, kl_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -hl_85[k]
                   + f_0 * kl_400[k];

        t_221[k] = -hl_86[k]
                   + f_0 * kl_401[k];

        t_222[k] = -hl_87[k]
                   + f_0 * kl_402[k];

        t_223[k] = -hl_88[k]
                   + f_0 * kl_403[k];

        t_224[k] = -hl_89[k]
                   + f_0 * kl_404[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, hl_90, hl_91, hl_92, hl_93, hl_94, \
                         kl_405, kl_406, kl_407, kl_408, kl_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -2.0 * hl_90[k]
                   + f_0 * kl_405[k];

        t_226[k] = -2.0 * hl_91[k]
                   + f_0 * kl_406[k];

        t_227[k] = -2.0 * hl_92[k]
                   + f_0 * kl_407[k];

        t_228[k] = -2.0 * hl_93[k]
                   + f_0 * kl_408[k];

        t_229[k] = -2.0 * hl_94[k]
                   + f_0 * kl_409[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, hl_95, hl_96, hl_97, hl_98, hl_99, \
                         kl_410, kl_411, kl_412, kl_413, kl_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -2.0 * hl_95[k]
                   + f_0 * kl_410[k];

        t_231[k] = -2.0 * hl_96[k]
                   + f_0 * kl_411[k];

        t_232[k] = -2.0 * hl_97[k]
                   + f_0 * kl_412[k];

        t_233[k] = -2.0 * hl_98[k]
                   + f_0 * kl_413[k];

        t_234[k] = -2.0 * hl_99[k]
                   + f_0 * kl_414[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, hl_100, hl_101, hl_102, hl_103, \
                         hl_104, kl_415, kl_416, kl_417, kl_418, \
                         kl_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -2.0 * hl_100[k]
                   + f_0 * kl_415[k];

        t_236[k] = -2.0 * hl_101[k]
                   + f_0 * kl_416[k];

        t_237[k] = -2.0 * hl_102[k]
                   + f_0 * kl_417[k];

        t_238[k] = -2.0 * hl_103[k]
                   + f_0 * kl_418[k];

        t_239[k] = -2.0 * hl_104[k]
                   + f_0 * kl_419[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, hl_105, hl_106, hl_107, hl_108, \
                         hl_109, kl_420, kl_421, kl_422, kl_423, \
                         kl_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -2.0 * hl_105[k]
                   + f_0 * kl_420[k];

        t_241[k] = -2.0 * hl_106[k]
                   + f_0 * kl_421[k];

        t_242[k] = -2.0 * hl_107[k]
                   + f_0 * kl_422[k];

        t_243[k] = -2.0 * hl_108[k]
                   + f_0 * kl_423[k];

        t_244[k] = -2.0 * hl_109[k]
                   + f_0 * kl_424[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, hl_110, hl_111, hl_112, hl_113, \
                         hl_114, kl_425, kl_426, kl_427, kl_428, \
                         kl_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -2.0 * hl_110[k]
                   + f_0 * kl_425[k];

        t_246[k] = -2.0 * hl_111[k]
                   + f_0 * kl_426[k];

        t_247[k] = -2.0 * hl_112[k]
                   + f_0 * kl_427[k];

        t_248[k] = -2.0 * hl_113[k]
                   + f_0 * kl_428[k];

        t_249[k] = -2.0 * hl_114[k]
                   + f_0 * kl_429[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, hl_115, hl_116, hl_117, hl_118, \
                         hl_119, kl_430, kl_431, kl_432, kl_433, \
                         kl_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * hl_115[k]
                   + f_0 * kl_430[k];

        t_251[k] = -2.0 * hl_116[k]
                   + f_0 * kl_431[k];

        t_252[k] = -2.0 * hl_117[k]
                   + f_0 * kl_432[k];

        t_253[k] = -2.0 * hl_118[k]
                   + f_0 * kl_433[k];

        t_254[k] = -2.0 * hl_119[k]
                   + f_0 * kl_434[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, hl_120, hl_121, hl_122, hl_123, \
                         hl_124, kl_435, kl_436, kl_437, kl_438, \
                         kl_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * hl_120[k]
                   + f_0 * kl_435[k];

        t_256[k] = -2.0 * hl_121[k]
                   + f_0 * kl_436[k];

        t_257[k] = -2.0 * hl_122[k]
                   + f_0 * kl_437[k];

        t_258[k] = -2.0 * hl_123[k]
                   + f_0 * kl_438[k];

        t_259[k] = -2.0 * hl_124[k]
                   + f_0 * kl_439[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, hl_125, hl_126, hl_127, hl_128, \
                         hl_129, kl_440, kl_441, kl_442, kl_443, \
                         kl_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * hl_125[k]
                   + f_0 * kl_440[k];

        t_261[k] = -2.0 * hl_126[k]
                   + f_0 * kl_441[k];

        t_262[k] = -2.0 * hl_127[k]
                   + f_0 * kl_442[k];

        t_263[k] = -2.0 * hl_128[k]
                   + f_0 * kl_443[k];

        t_264[k] = -2.0 * hl_129[k]
                   + f_0 * kl_444[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, hl_130, hl_131, hl_132, hl_133, \
                         hl_134, kl_445, kl_446, kl_447, kl_448, \
                         kl_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * hl_130[k]
                   + f_0 * kl_445[k];

        t_266[k] = -2.0 * hl_131[k]
                   + f_0 * kl_446[k];

        t_267[k] = -2.0 * hl_132[k]
                   + f_0 * kl_447[k];

        t_268[k] = -2.0 * hl_133[k]
                   + f_0 * kl_448[k];

        t_269[k] = -2.0 * hl_134[k]
                   + f_0 * kl_449[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, t_276, t_277, kl_495, \
                         kl_496, kl_497, kl_498, kl_499, kl_500, kl_501, \
                         kl_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_0 * kl_495[k];

        t_271[k] = f_0 * kl_496[k];

        t_272[k] = f_0 * kl_497[k];

        t_273[k] = f_0 * kl_498[k];

        t_274[k] = f_0 * kl_499[k];

        t_275[k] = f_0 * kl_500[k];

        t_276[k] = f_0 * kl_501[k];

        t_277[k] = f_0 * kl_502[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, t_283, t_284, t_285, kl_503, \
                         kl_504, kl_505, kl_506, kl_507, kl_508, kl_509, \
                         kl_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_0 * kl_503[k];

        t_279[k] = f_0 * kl_504[k];

        t_280[k] = f_0 * kl_505[k];

        t_281[k] = f_0 * kl_506[k];

        t_282[k] = f_0 * kl_507[k];

        t_283[k] = f_0 * kl_508[k];

        t_284[k] = f_0 * kl_509[k];

        t_285[k] = f_0 * kl_510[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, t_291, t_292, t_293, kl_511, \
                         kl_512, kl_513, kl_514, kl_515, kl_516, kl_517, \
                         kl_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * kl_511[k];

        t_287[k] = f_0 * kl_512[k];

        t_288[k] = f_0 * kl_513[k];

        t_289[k] = f_0 * kl_514[k];

        t_290[k] = f_0 * kl_515[k];

        t_291[k] = f_0 * kl_516[k];

        t_292[k] = f_0 * kl_517[k];

        t_293[k] = f_0 * kl_518[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, t_299, t_300, t_301, kl_519, \
                         kl_520, kl_521, kl_522, kl_523, kl_524, kl_525, \
                         kl_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_0 * kl_519[k];

        t_295[k] = f_0 * kl_520[k];

        t_296[k] = f_0 * kl_521[k];

        t_297[k] = f_0 * kl_522[k];

        t_298[k] = f_0 * kl_523[k];

        t_299[k] = f_0 * kl_524[k];

        t_300[k] = f_0 * kl_525[k];

        t_301[k] = f_0 * kl_526[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, t_307, t_308, t_309, kl_527, \
                         kl_528, kl_529, kl_530, kl_531, kl_532, kl_533, \
                         kl_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_0 * kl_527[k];

        t_303[k] = f_0 * kl_528[k];

        t_304[k] = f_0 * kl_529[k];

        t_305[k] = f_0 * kl_530[k];

        t_306[k] = f_0 * kl_531[k];

        t_307[k] = f_0 * kl_532[k];

        t_308[k] = f_0 * kl_533[k];

        t_309[k] = f_0 * kl_534[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, t_315, t_316, hl_135, hl_136, \
                         kl_535, kl_536, kl_537, kl_538, kl_539, kl_540, \
                         kl_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_0 * kl_535[k];

        t_311[k] = f_0 * kl_536[k];

        t_312[k] = f_0 * kl_537[k];

        t_313[k] = f_0 * kl_538[k];

        t_314[k] = f_0 * kl_539[k];

        t_315[k] = -hl_135[k]
                   + f_0 * kl_540[k];

        t_316[k] = -hl_136[k]
                   + f_0 * kl_541[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, hl_137, hl_138, hl_139, hl_140, \
                         hl_141, kl_542, kl_543, kl_544, kl_545, \
                         kl_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = -hl_137[k]
                   + f_0 * kl_542[k];

        t_318[k] = -hl_138[k]
                   + f_0 * kl_543[k];

        t_319[k] = -hl_139[k]
                   + f_0 * kl_544[k];

        t_320[k] = -hl_140[k]
                   + f_0 * kl_545[k];

        t_321[k] = -hl_141[k]
                   + f_0 * kl_546[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, hl_142, hl_143, hl_144, hl_145, \
                         hl_146, kl_547, kl_548, kl_549, kl_550, \
                         kl_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = -hl_142[k]
                   + f_0 * kl_547[k];

        t_323[k] = -hl_143[k]
                   + f_0 * kl_548[k];

        t_324[k] = -hl_144[k]
                   + f_0 * kl_549[k];

        t_325[k] = -hl_145[k]
                   + f_0 * kl_550[k];

        t_326[k] = -hl_146[k]
                   + f_0 * kl_551[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, hl_147, hl_148, hl_149, hl_150, \
                         hl_151, kl_552, kl_553, kl_554, kl_555, \
                         kl_556 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = -hl_147[k]
                   + f_0 * kl_552[k];

        t_328[k] = -hl_148[k]
                   + f_0 * kl_553[k];

        t_329[k] = -hl_149[k]
                   + f_0 * kl_554[k];

        t_330[k] = -hl_150[k]
                   + f_0 * kl_555[k];

        t_331[k] = -hl_151[k]
                   + f_0 * kl_556[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, hl_152, hl_153, hl_154, hl_155, \
                         hl_156, kl_557, kl_558, kl_559, kl_560, \
                         kl_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -hl_152[k]
                   + f_0 * kl_557[k];

        t_333[k] = -hl_153[k]
                   + f_0 * kl_558[k];

        t_334[k] = -hl_154[k]
                   + f_0 * kl_559[k];

        t_335[k] = -hl_155[k]
                   + f_0 * kl_560[k];

        t_336[k] = -hl_156[k]
                   + f_0 * kl_561[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, hl_157, hl_158, hl_159, hl_160, \
                         hl_161, kl_562, kl_563, kl_564, kl_565, \
                         kl_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -hl_157[k]
                   + f_0 * kl_562[k];

        t_338[k] = -hl_158[k]
                   + f_0 * kl_563[k];

        t_339[k] = -hl_159[k]
                   + f_0 * kl_564[k];

        t_340[k] = -hl_160[k]
                   + f_0 * kl_565[k];

        t_341[k] = -hl_161[k]
                   + f_0 * kl_566[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, hl_162, hl_163, hl_164, hl_165, \
                         hl_166, kl_567, kl_568, kl_569, kl_570, \
                         kl_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -hl_162[k]
                   + f_0 * kl_567[k];

        t_343[k] = -hl_163[k]
                   + f_0 * kl_568[k];

        t_344[k] = -hl_164[k]
                   + f_0 * kl_569[k];

        t_345[k] = -hl_165[k]
                   + f_0 * kl_570[k];

        t_346[k] = -hl_166[k]
                   + f_0 * kl_571[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, hl_167, hl_168, hl_169, hl_170, \
                         hl_171, kl_572, kl_573, kl_574, kl_575, \
                         kl_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -hl_167[k]
                   + f_0 * kl_572[k];

        t_348[k] = -hl_168[k]
                   + f_0 * kl_573[k];

        t_349[k] = -hl_169[k]
                   + f_0 * kl_574[k];

        t_350[k] = -hl_170[k]
                   + f_0 * kl_575[k];

        t_351[k] = -hl_171[k]
                   + f_0 * kl_576[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, hl_172, hl_173, hl_174, hl_175, \
                         hl_176, kl_577, kl_578, kl_579, kl_580, \
                         kl_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -hl_172[k]
                   + f_0 * kl_577[k];

        t_353[k] = -hl_173[k]
                   + f_0 * kl_578[k];

        t_354[k] = -hl_174[k]
                   + f_0 * kl_579[k];

        t_355[k] = -hl_175[k]
                   + f_0 * kl_580[k];

        t_356[k] = -hl_176[k]
                   + f_0 * kl_581[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_511 = buffer.data(target + 511);
    auto *t_512 = buffer.data(target + 512);
    auto *t_513 = buffer.data(target + 513);
    auto *t_514 = buffer.data(target + 514);
    auto *t_515 = buffer.data(target + 515);
    auto *t_516 = buffer.data(target + 516);
    auto *t_517 = buffer.data(target + 517);
    auto *t_518 = buffer.data(target + 518);

    const auto *hl_177 = buffer.data(hl + 177);
    const auto *hl_178 = buffer.data(hl + 178);
    const auto *hl_179 = buffer.data(hl + 179);
    const auto *hl_180 = buffer.data(hl + 180);
    const auto *hl_181 = buffer.data(hl + 181);
    const auto *hl_182 = buffer.data(hl + 182);
    const auto *hl_183 = buffer.data(hl + 183);
    const auto *hl_184 = buffer.data(hl + 184);
    const auto *hl_185 = buffer.data(hl + 185);
    const auto *hl_186 = buffer.data(hl + 186);
    const auto *hl_187 = buffer.data(hl + 187);
    const auto *hl_188 = buffer.data(hl + 188);
    const auto *hl_189 = buffer.data(hl + 189);
    const auto *hl_190 = buffer.data(hl + 190);
    const auto *hl_191 = buffer.data(hl + 191);
    const auto *hl_192 = buffer.data(hl + 192);
    const auto *hl_193 = buffer.data(hl + 193);
    const auto *hl_194 = buffer.data(hl + 194);
    const auto *hl_195 = buffer.data(hl + 195);
    const auto *hl_196 = buffer.data(hl + 196);
    const auto *hl_197 = buffer.data(hl + 197);
    const auto *hl_198 = buffer.data(hl + 198);
    const auto *hl_199 = buffer.data(hl + 199);
    const auto *hl_200 = buffer.data(hl + 200);
    const auto *hl_201 = buffer.data(hl + 201);
    const auto *hl_202 = buffer.data(hl + 202);
    const auto *hl_203 = buffer.data(hl + 203);
    const auto *hl_204 = buffer.data(hl + 204);
    const auto *hl_205 = buffer.data(hl + 205);
    const auto *hl_206 = buffer.data(hl + 206);
    const auto *hl_207 = buffer.data(hl + 207);
    const auto *hl_208 = buffer.data(hl + 208);
    const auto *hl_209 = buffer.data(hl + 209);
    const auto *hl_210 = buffer.data(hl + 210);
    const auto *hl_211 = buffer.data(hl + 211);
    const auto *hl_212 = buffer.data(hl + 212);
    const auto *hl_213 = buffer.data(hl + 213);
    const auto *hl_214 = buffer.data(hl + 214);
    const auto *hl_215 = buffer.data(hl + 215);
    const auto *hl_216 = buffer.data(hl + 216);
    const auto *hl_217 = buffer.data(hl + 217);
    const auto *hl_218 = buffer.data(hl + 218);
    const auto *hl_219 = buffer.data(hl + 219);
    const auto *hl_220 = buffer.data(hl + 220);
    const auto *hl_221 = buffer.data(hl + 221);
    const auto *hl_222 = buffer.data(hl + 222);
    const auto *hl_223 = buffer.data(hl + 223);
    const auto *hl_224 = buffer.data(hl + 224);
    const auto *hl_225 = buffer.data(hl + 225);
    const auto *hl_226 = buffer.data(hl + 226);
    const auto *hl_227 = buffer.data(hl + 227);
    const auto *hl_228 = buffer.data(hl + 228);
    const auto *hl_229 = buffer.data(hl + 229);
    const auto *hl_230 = buffer.data(hl + 230);
    const auto *hl_231 = buffer.data(hl + 231);
    const auto *hl_232 = buffer.data(hl + 232);
    const auto *hl_233 = buffer.data(hl + 233);
    const auto *hl_234 = buffer.data(hl + 234);
    const auto *hl_235 = buffer.data(hl + 235);
    const auto *hl_236 = buffer.data(hl + 236);
    const auto *hl_237 = buffer.data(hl + 237);
    const auto *hl_238 = buffer.data(hl + 238);
    const auto *hl_239 = buffer.data(hl + 239);
    const auto *hl_240 = buffer.data(hl + 240);
    const auto *hl_241 = buffer.data(hl + 241);
    const auto *hl_242 = buffer.data(hl + 242);
    const auto *hl_243 = buffer.data(hl + 243);
    const auto *hl_244 = buffer.data(hl + 244);
    const auto *hl_245 = buffer.data(hl + 245);
    const auto *hl_246 = buffer.data(hl + 246);
    const auto *hl_247 = buffer.data(hl + 247);
    const auto *hl_248 = buffer.data(hl + 248);
    const auto *hl_249 = buffer.data(hl + 249);
    const auto *hl_250 = buffer.data(hl + 250);
    const auto *hl_251 = buffer.data(hl + 251);
    const auto *hl_252 = buffer.data(hl + 252);
    const auto *hl_253 = buffer.data(hl + 253);
    const auto *hl_254 = buffer.data(hl + 254);
    const auto *hl_255 = buffer.data(hl + 255);
    const auto *hl_256 = buffer.data(hl + 256);
    const auto *hl_257 = buffer.data(hl + 257);
    const auto *hl_258 = buffer.data(hl + 258);
    const auto *hl_259 = buffer.data(hl + 259);
    const auto *hl_260 = buffer.data(hl + 260);
    const auto *hl_261 = buffer.data(hl + 261);
    const auto *hl_262 = buffer.data(hl + 262);
    const auto *hl_263 = buffer.data(hl + 263);
    const auto *hl_264 = buffer.data(hl + 264);
    const auto *hl_265 = buffer.data(hl + 265);
    const auto *hl_266 = buffer.data(hl + 266);
    const auto *hl_267 = buffer.data(hl + 267);
    const auto *hl_268 = buffer.data(hl + 268);
    const auto *hl_269 = buffer.data(hl + 269);
    const auto *hl_270 = buffer.data(hl + 270);
    const auto *hl_271 = buffer.data(hl + 271);
    const auto *hl_272 = buffer.data(hl + 272);
    const auto *hl_273 = buffer.data(hl + 273);
    const auto *hl_274 = buffer.data(hl + 274);
    const auto *hl_275 = buffer.data(hl + 275);
    const auto *hl_276 = buffer.data(hl + 276);
    const auto *hl_277 = buffer.data(hl + 277);
    const auto *hl_278 = buffer.data(hl + 278);
    const auto *hl_279 = buffer.data(hl + 279);
    const auto *hl_280 = buffer.data(hl + 280);
    const auto *hl_281 = buffer.data(hl + 281);
    const auto *hl_282 = buffer.data(hl + 282);
    const auto *hl_283 = buffer.data(hl + 283);
    const auto *hl_284 = buffer.data(hl + 284);
    const auto *hl_285 = buffer.data(hl + 285);
    const auto *hl_286 = buffer.data(hl + 286);
    const auto *hl_287 = buffer.data(hl + 287);
    const auto *hl_288 = buffer.data(hl + 288);
    const auto *hl_289 = buffer.data(hl + 289);
    const auto *hl_290 = buffer.data(hl + 290);
    const auto *hl_291 = buffer.data(hl + 291);
    const auto *hl_292 = buffer.data(hl + 292);
    const auto *hl_293 = buffer.data(hl + 293);

    const auto *kl_582 = buffer.data(kl + 582);
    const auto *kl_583 = buffer.data(kl + 583);
    const auto *kl_584 = buffer.data(kl + 584);
    const auto *kl_585 = buffer.data(kl + 585);
    const auto *kl_586 = buffer.data(kl + 586);
    const auto *kl_587 = buffer.data(kl + 587);
    const auto *kl_588 = buffer.data(kl + 588);
    const auto *kl_589 = buffer.data(kl + 589);
    const auto *kl_590 = buffer.data(kl + 590);
    const auto *kl_591 = buffer.data(kl + 591);
    const auto *kl_592 = buffer.data(kl + 592);
    const auto *kl_593 = buffer.data(kl + 593);
    const auto *kl_594 = buffer.data(kl + 594);
    const auto *kl_595 = buffer.data(kl + 595);
    const auto *kl_596 = buffer.data(kl + 596);
    const auto *kl_597 = buffer.data(kl + 597);
    const auto *kl_598 = buffer.data(kl + 598);
    const auto *kl_599 = buffer.data(kl + 599);
    const auto *kl_600 = buffer.data(kl + 600);
    const auto *kl_601 = buffer.data(kl + 601);
    const auto *kl_602 = buffer.data(kl + 602);
    const auto *kl_603 = buffer.data(kl + 603);
    const auto *kl_604 = buffer.data(kl + 604);
    const auto *kl_605 = buffer.data(kl + 605);
    const auto *kl_606 = buffer.data(kl + 606);
    const auto *kl_607 = buffer.data(kl + 607);
    const auto *kl_608 = buffer.data(kl + 608);
    const auto *kl_609 = buffer.data(kl + 609);
    const auto *kl_610 = buffer.data(kl + 610);
    const auto *kl_611 = buffer.data(kl + 611);
    const auto *kl_612 = buffer.data(kl + 612);
    const auto *kl_613 = buffer.data(kl + 613);
    const auto *kl_614 = buffer.data(kl + 614);
    const auto *kl_615 = buffer.data(kl + 615);
    const auto *kl_616 = buffer.data(kl + 616);
    const auto *kl_617 = buffer.data(kl + 617);
    const auto *kl_618 = buffer.data(kl + 618);
    const auto *kl_619 = buffer.data(kl + 619);
    const auto *kl_620 = buffer.data(kl + 620);
    const auto *kl_621 = buffer.data(kl + 621);
    const auto *kl_622 = buffer.data(kl + 622);
    const auto *kl_623 = buffer.data(kl + 623);
    const auto *kl_624 = buffer.data(kl + 624);
    const auto *kl_625 = buffer.data(kl + 625);
    const auto *kl_626 = buffer.data(kl + 626);
    const auto *kl_627 = buffer.data(kl + 627);
    const auto *kl_628 = buffer.data(kl + 628);
    const auto *kl_629 = buffer.data(kl + 629);
    const auto *kl_630 = buffer.data(kl + 630);
    const auto *kl_631 = buffer.data(kl + 631);
    const auto *kl_632 = buffer.data(kl + 632);
    const auto *kl_633 = buffer.data(kl + 633);
    const auto *kl_634 = buffer.data(kl + 634);
    const auto *kl_635 = buffer.data(kl + 635);
    const auto *kl_636 = buffer.data(kl + 636);
    const auto *kl_637 = buffer.data(kl + 637);
    const auto *kl_638 = buffer.data(kl + 638);
    const auto *kl_639 = buffer.data(kl + 639);
    const auto *kl_640 = buffer.data(kl + 640);
    const auto *kl_641 = buffer.data(kl + 641);
    const auto *kl_642 = buffer.data(kl + 642);
    const auto *kl_643 = buffer.data(kl + 643);
    const auto *kl_644 = buffer.data(kl + 644);
    const auto *kl_645 = buffer.data(kl + 645);
    const auto *kl_646 = buffer.data(kl + 646);
    const auto *kl_647 = buffer.data(kl + 647);
    const auto *kl_648 = buffer.data(kl + 648);
    const auto *kl_649 = buffer.data(kl + 649);
    const auto *kl_650 = buffer.data(kl + 650);
    const auto *kl_651 = buffer.data(kl + 651);
    const auto *kl_652 = buffer.data(kl + 652);
    const auto *kl_653 = buffer.data(kl + 653);
    const auto *kl_654 = buffer.data(kl + 654);
    const auto *kl_655 = buffer.data(kl + 655);
    const auto *kl_656 = buffer.data(kl + 656);
    const auto *kl_657 = buffer.data(kl + 657);
    const auto *kl_658 = buffer.data(kl + 658);
    const auto *kl_659 = buffer.data(kl + 659);
    const auto *kl_660 = buffer.data(kl + 660);
    const auto *kl_661 = buffer.data(kl + 661);
    const auto *kl_662 = buffer.data(kl + 662);
    const auto *kl_663 = buffer.data(kl + 663);
    const auto *kl_664 = buffer.data(kl + 664);
    const auto *kl_665 = buffer.data(kl + 665);
    const auto *kl_666 = buffer.data(kl + 666);
    const auto *kl_667 = buffer.data(kl + 667);
    const auto *kl_668 = buffer.data(kl + 668);
    const auto *kl_669 = buffer.data(kl + 669);
    const auto *kl_670 = buffer.data(kl + 670);
    const auto *kl_671 = buffer.data(kl + 671);
    const auto *kl_672 = buffer.data(kl + 672);
    const auto *kl_673 = buffer.data(kl + 673);
    const auto *kl_674 = buffer.data(kl + 674);
    const auto *kl_720 = buffer.data(kl + 720);
    const auto *kl_721 = buffer.data(kl + 721);
    const auto *kl_722 = buffer.data(kl + 722);
    const auto *kl_723 = buffer.data(kl + 723);
    const auto *kl_724 = buffer.data(kl + 724);
    const auto *kl_725 = buffer.data(kl + 725);
    const auto *kl_726 = buffer.data(kl + 726);
    const auto *kl_727 = buffer.data(kl + 727);
    const auto *kl_728 = buffer.data(kl + 728);
    const auto *kl_729 = buffer.data(kl + 729);
    const auto *kl_730 = buffer.data(kl + 730);
    const auto *kl_731 = buffer.data(kl + 731);
    const auto *kl_732 = buffer.data(kl + 732);
    const auto *kl_733 = buffer.data(kl + 733);
    const auto *kl_734 = buffer.data(kl + 734);
    const auto *kl_735 = buffer.data(kl + 735);
    const auto *kl_736 = buffer.data(kl + 736);
    const auto *kl_737 = buffer.data(kl + 737);
    const auto *kl_738 = buffer.data(kl + 738);
    const auto *kl_739 = buffer.data(kl + 739);
    const auto *kl_740 = buffer.data(kl + 740);
    const auto *kl_741 = buffer.data(kl + 741);
    const auto *kl_742 = buffer.data(kl + 742);
    const auto *kl_743 = buffer.data(kl + 743);
    const auto *kl_744 = buffer.data(kl + 744);
    const auto *kl_745 = buffer.data(kl + 745);
    const auto *kl_746 = buffer.data(kl + 746);
    const auto *kl_747 = buffer.data(kl + 747);
    const auto *kl_748 = buffer.data(kl + 748);
    const auto *kl_749 = buffer.data(kl + 749);
    const auto *kl_750 = buffer.data(kl + 750);
    const auto *kl_751 = buffer.data(kl + 751);
    const auto *kl_752 = buffer.data(kl + 752);
    const auto *kl_753 = buffer.data(kl + 753);
    const auto *kl_754 = buffer.data(kl + 754);
    const auto *kl_755 = buffer.data(kl + 755);
    const auto *kl_756 = buffer.data(kl + 756);
    const auto *kl_757 = buffer.data(kl + 757);
    const auto *kl_758 = buffer.data(kl + 758);
    const auto *kl_759 = buffer.data(kl + 759);
    const auto *kl_760 = buffer.data(kl + 760);
    const auto *kl_761 = buffer.data(kl + 761);
    const auto *kl_762 = buffer.data(kl + 762);
    const auto *kl_763 = buffer.data(kl + 763);
    const auto *kl_764 = buffer.data(kl + 764);
    const auto *kl_765 = buffer.data(kl + 765);
    const auto *kl_766 = buffer.data(kl + 766);
    const auto *kl_767 = buffer.data(kl + 767);
    const auto *kl_768 = buffer.data(kl + 768);
    const auto *kl_769 = buffer.data(kl + 769);
    const auto *kl_770 = buffer.data(kl + 770);
    const auto *kl_771 = buffer.data(kl + 771);
    const auto *kl_772 = buffer.data(kl + 772);
    const auto *kl_773 = buffer.data(kl + 773);
    const auto *kl_774 = buffer.data(kl + 774);
    const auto *kl_775 = buffer.data(kl + 775);
    const auto *kl_776 = buffer.data(kl + 776);
    const auto *kl_777 = buffer.data(kl + 777);
    const auto *kl_778 = buffer.data(kl + 778);
    const auto *kl_779 = buffer.data(kl + 779);
    const auto *kl_780 = buffer.data(kl + 780);
    const auto *kl_781 = buffer.data(kl + 781);
    const auto *kl_782 = buffer.data(kl + 782);
    const auto *kl_783 = buffer.data(kl + 783);
    const auto *kl_784 = buffer.data(kl + 784);
    const auto *kl_785 = buffer.data(kl + 785);
    const auto *kl_786 = buffer.data(kl + 786);
    const auto *kl_787 = buffer.data(kl + 787);
    const auto *kl_788 = buffer.data(kl + 788);

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, hl_177, hl_178, hl_179, hl_180, \
                         hl_181, kl_582, kl_583, kl_584, kl_585, \
                         kl_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -hl_177[k]
                   + f_0 * kl_582[k];

        t_358[k] = -hl_178[k]
                   + f_0 * kl_583[k];

        t_359[k] = -hl_179[k]
                   + f_0 * kl_584[k];

        t_360[k] = -2.0 * hl_180[k]
                   + f_0 * kl_585[k];

        t_361[k] = -2.0 * hl_181[k]
                   + f_0 * kl_586[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, hl_182, hl_183, hl_184, hl_185, \
                         hl_186, kl_587, kl_588, kl_589, kl_590, \
                         kl_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -2.0 * hl_182[k]
                   + f_0 * kl_587[k];

        t_363[k] = -2.0 * hl_183[k]
                   + f_0 * kl_588[k];

        t_364[k] = -2.0 * hl_184[k]
                   + f_0 * kl_589[k];

        t_365[k] = -2.0 * hl_185[k]
                   + f_0 * kl_590[k];

        t_366[k] = -2.0 * hl_186[k]
                   + f_0 * kl_591[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, hl_187, hl_188, hl_189, hl_190, \
                         hl_191, kl_592, kl_593, kl_594, kl_595, \
                         kl_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -2.0 * hl_187[k]
                   + f_0 * kl_592[k];

        t_368[k] = -2.0 * hl_188[k]
                   + f_0 * kl_593[k];

        t_369[k] = -2.0 * hl_189[k]
                   + f_0 * kl_594[k];

        t_370[k] = -2.0 * hl_190[k]
                   + f_0 * kl_595[k];

        t_371[k] = -2.0 * hl_191[k]
                   + f_0 * kl_596[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, hl_192, hl_193, hl_194, hl_195, \
                         hl_196, kl_597, kl_598, kl_599, kl_600, \
                         kl_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -2.0 * hl_192[k]
                   + f_0 * kl_597[k];

        t_373[k] = -2.0 * hl_193[k]
                   + f_0 * kl_598[k];

        t_374[k] = -2.0 * hl_194[k]
                   + f_0 * kl_599[k];

        t_375[k] = -2.0 * hl_195[k]
                   + f_0 * kl_600[k];

        t_376[k] = -2.0 * hl_196[k]
                   + f_0 * kl_601[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, hl_197, hl_198, hl_199, hl_200, \
                         hl_201, kl_602, kl_603, kl_604, kl_605, \
                         kl_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -2.0 * hl_197[k]
                   + f_0 * kl_602[k];

        t_378[k] = -2.0 * hl_198[k]
                   + f_0 * kl_603[k];

        t_379[k] = -2.0 * hl_199[k]
                   + f_0 * kl_604[k];

        t_380[k] = -2.0 * hl_200[k]
                   + f_0 * kl_605[k];

        t_381[k] = -2.0 * hl_201[k]
                   + f_0 * kl_606[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, hl_202, hl_203, hl_204, hl_205, \
                         hl_206, kl_607, kl_608, kl_609, kl_610, \
                         kl_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -2.0 * hl_202[k]
                   + f_0 * kl_607[k];

        t_383[k] = -2.0 * hl_203[k]
                   + f_0 * kl_608[k];

        t_384[k] = -2.0 * hl_204[k]
                   + f_0 * kl_609[k];

        t_385[k] = -2.0 * hl_205[k]
                   + f_0 * kl_610[k];

        t_386[k] = -2.0 * hl_206[k]
                   + f_0 * kl_611[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, hl_207, hl_208, hl_209, hl_210, \
                         hl_211, kl_612, kl_613, kl_614, kl_615, \
                         kl_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -2.0 * hl_207[k]
                   + f_0 * kl_612[k];

        t_388[k] = -2.0 * hl_208[k]
                   + f_0 * kl_613[k];

        t_389[k] = -2.0 * hl_209[k]
                   + f_0 * kl_614[k];

        t_390[k] = -2.0 * hl_210[k]
                   + f_0 * kl_615[k];

        t_391[k] = -2.0 * hl_211[k]
                   + f_0 * kl_616[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, hl_212, hl_213, hl_214, hl_215, \
                         hl_216, kl_617, kl_618, kl_619, kl_620, \
                         kl_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -2.0 * hl_212[k]
                   + f_0 * kl_617[k];

        t_393[k] = -2.0 * hl_213[k]
                   + f_0 * kl_618[k];

        t_394[k] = -2.0 * hl_214[k]
                   + f_0 * kl_619[k];

        t_395[k] = -2.0 * hl_215[k]
                   + f_0 * kl_620[k];

        t_396[k] = -2.0 * hl_216[k]
                   + f_0 * kl_621[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, hl_217, hl_218, hl_219, hl_220, \
                         hl_221, kl_622, kl_623, kl_624, kl_625, \
                         kl_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -2.0 * hl_217[k]
                   + f_0 * kl_622[k];

        t_398[k] = -2.0 * hl_218[k]
                   + f_0 * kl_623[k];

        t_399[k] = -2.0 * hl_219[k]
                   + f_0 * kl_624[k];

        t_400[k] = -2.0 * hl_220[k]
                   + f_0 * kl_625[k];

        t_401[k] = -2.0 * hl_221[k]
                   + f_0 * kl_626[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, hl_222, hl_223, hl_224, hl_225, \
                         hl_226, kl_627, kl_628, kl_629, kl_630, \
                         kl_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -2.0 * hl_222[k]
                   + f_0 * kl_627[k];

        t_403[k] = -2.0 * hl_223[k]
                   + f_0 * kl_628[k];

        t_404[k] = -2.0 * hl_224[k]
                   + f_0 * kl_629[k];

        t_405[k] = -3.0 * hl_225[k]
                   + f_0 * kl_630[k];

        t_406[k] = -3.0 * hl_226[k]
                   + f_0 * kl_631[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, hl_227, hl_228, hl_229, hl_230, \
                         hl_231, kl_632, kl_633, kl_634, kl_635, \
                         kl_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -3.0 * hl_227[k]
                   + f_0 * kl_632[k];

        t_408[k] = -3.0 * hl_228[k]
                   + f_0 * kl_633[k];

        t_409[k] = -3.0 * hl_229[k]
                   + f_0 * kl_634[k];

        t_410[k] = -3.0 * hl_230[k]
                   + f_0 * kl_635[k];

        t_411[k] = -3.0 * hl_231[k]
                   + f_0 * kl_636[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, hl_232, hl_233, hl_234, hl_235, \
                         hl_236, kl_637, kl_638, kl_639, kl_640, \
                         kl_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -3.0 * hl_232[k]
                   + f_0 * kl_637[k];

        t_413[k] = -3.0 * hl_233[k]
                   + f_0 * kl_638[k];

        t_414[k] = -3.0 * hl_234[k]
                   + f_0 * kl_639[k];

        t_415[k] = -3.0 * hl_235[k]
                   + f_0 * kl_640[k];

        t_416[k] = -3.0 * hl_236[k]
                   + f_0 * kl_641[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, hl_237, hl_238, hl_239, hl_240, \
                         hl_241, kl_642, kl_643, kl_644, kl_645, \
                         kl_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -3.0 * hl_237[k]
                   + f_0 * kl_642[k];

        t_418[k] = -3.0 * hl_238[k]
                   + f_0 * kl_643[k];

        t_419[k] = -3.0 * hl_239[k]
                   + f_0 * kl_644[k];

        t_420[k] = -3.0 * hl_240[k]
                   + f_0 * kl_645[k];

        t_421[k] = -3.0 * hl_241[k]
                   + f_0 * kl_646[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, t_426, hl_242, hl_243, hl_244, hl_245, \
                         hl_246, kl_647, kl_648, kl_649, kl_650, \
                         kl_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = -3.0 * hl_242[k]
                   + f_0 * kl_647[k];

        t_423[k] = -3.0 * hl_243[k]
                   + f_0 * kl_648[k];

        t_424[k] = -3.0 * hl_244[k]
                   + f_0 * kl_649[k];

        t_425[k] = -3.0 * hl_245[k]
                   + f_0 * kl_650[k];

        t_426[k] = -3.0 * hl_246[k]
                   + f_0 * kl_651[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, t_431, hl_247, hl_248, hl_249, hl_250, \
                         hl_251, kl_652, kl_653, kl_654, kl_655, \
                         kl_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = -3.0 * hl_247[k]
                   + f_0 * kl_652[k];

        t_428[k] = -3.0 * hl_248[k]
                   + f_0 * kl_653[k];

        t_429[k] = -3.0 * hl_249[k]
                   + f_0 * kl_654[k];

        t_430[k] = -3.0 * hl_250[k]
                   + f_0 * kl_655[k];

        t_431[k] = -3.0 * hl_251[k]
                   + f_0 * kl_656[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, t_436, hl_252, hl_253, hl_254, hl_255, \
                         hl_256, kl_657, kl_658, kl_659, kl_660, \
                         kl_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = -3.0 * hl_252[k]
                   + f_0 * kl_657[k];

        t_433[k] = -3.0 * hl_253[k]
                   + f_0 * kl_658[k];

        t_434[k] = -3.0 * hl_254[k]
                   + f_0 * kl_659[k];

        t_435[k] = -3.0 * hl_255[k]
                   + f_0 * kl_660[k];

        t_436[k] = -3.0 * hl_256[k]
                   + f_0 * kl_661[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, t_441, hl_257, hl_258, hl_259, hl_260, \
                         hl_261, kl_662, kl_663, kl_664, kl_665, \
                         kl_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = -3.0 * hl_257[k]
                   + f_0 * kl_662[k];

        t_438[k] = -3.0 * hl_258[k]
                   + f_0 * kl_663[k];

        t_439[k] = -3.0 * hl_259[k]
                   + f_0 * kl_664[k];

        t_440[k] = -3.0 * hl_260[k]
                   + f_0 * kl_665[k];

        t_441[k] = -3.0 * hl_261[k]
                   + f_0 * kl_666[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, t_446, hl_262, hl_263, hl_264, hl_265, \
                         hl_266, kl_667, kl_668, kl_669, kl_670, \
                         kl_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = -3.0 * hl_262[k]
                   + f_0 * kl_667[k];

        t_443[k] = -3.0 * hl_263[k]
                   + f_0 * kl_668[k];

        t_444[k] = -3.0 * hl_264[k]
                   + f_0 * kl_669[k];

        t_445[k] = -3.0 * hl_265[k]
                   + f_0 * kl_670[k];

        t_446[k] = -3.0 * hl_266[k]
                   + f_0 * kl_671[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, t_451, t_452, hl_267, hl_268, hl_269, \
                         kl_672, kl_673, kl_674, kl_720, kl_721, \
                         kl_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = -3.0 * hl_267[k]
                   + f_0 * kl_672[k];

        t_448[k] = -3.0 * hl_268[k]
                   + f_0 * kl_673[k];

        t_449[k] = -3.0 * hl_269[k]
                   + f_0 * kl_674[k];

        t_450[k] = f_0 * kl_720[k];

        t_451[k] = f_0 * kl_721[k];

        t_452[k] = f_0 * kl_722[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, t_457, t_458, t_459, t_460, kl_723, \
                         kl_724, kl_725, kl_726, kl_727, kl_728, kl_729, \
                         kl_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_0 * kl_723[k];

        t_454[k] = f_0 * kl_724[k];

        t_455[k] = f_0 * kl_725[k];

        t_456[k] = f_0 * kl_726[k];

        t_457[k] = f_0 * kl_727[k];

        t_458[k] = f_0 * kl_728[k];

        t_459[k] = f_0 * kl_729[k];

        t_460[k] = f_0 * kl_730[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, t_465, t_466, t_467, t_468, kl_731, \
                         kl_732, kl_733, kl_734, kl_735, kl_736, kl_737, \
                         kl_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_0 * kl_731[k];

        t_462[k] = f_0 * kl_732[k];

        t_463[k] = f_0 * kl_733[k];

        t_464[k] = f_0 * kl_734[k];

        t_465[k] = f_0 * kl_735[k];

        t_466[k] = f_0 * kl_736[k];

        t_467[k] = f_0 * kl_737[k];

        t_468[k] = f_0 * kl_738[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, t_473, t_474, t_475, t_476, kl_739, \
                         kl_740, kl_741, kl_742, kl_743, kl_744, kl_745, \
                         kl_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_0 * kl_739[k];

        t_470[k] = f_0 * kl_740[k];

        t_471[k] = f_0 * kl_741[k];

        t_472[k] = f_0 * kl_742[k];

        t_473[k] = f_0 * kl_743[k];

        t_474[k] = f_0 * kl_744[k];

        t_475[k] = f_0 * kl_745[k];

        t_476[k] = f_0 * kl_746[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, t_482, t_483, t_484, kl_747, \
                         kl_748, kl_749, kl_750, kl_751, kl_752, kl_753, \
                         kl_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_0 * kl_747[k];

        t_478[k] = f_0 * kl_748[k];

        t_479[k] = f_0 * kl_749[k];

        t_480[k] = f_0 * kl_750[k];

        t_481[k] = f_0 * kl_751[k];

        t_482[k] = f_0 * kl_752[k];

        t_483[k] = f_0 * kl_753[k];

        t_484[k] = f_0 * kl_754[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, t_490, t_491, t_492, kl_755, \
                         kl_756, kl_757, kl_758, kl_759, kl_760, kl_761, \
                         kl_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_0 * kl_755[k];

        t_486[k] = f_0 * kl_756[k];

        t_487[k] = f_0 * kl_757[k];

        t_488[k] = f_0 * kl_758[k];

        t_489[k] = f_0 * kl_759[k];

        t_490[k] = f_0 * kl_760[k];

        t_491[k] = f_0 * kl_761[k];

        t_492[k] = f_0 * kl_762[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, t_496, t_497, t_498, hl_270, hl_271, hl_272, \
                         hl_273, kl_763, kl_764, kl_765, kl_766, kl_767, \
                         kl_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = f_0 * kl_763[k];

        t_494[k] = f_0 * kl_764[k];

        t_495[k] = -hl_270[k]
                   + f_0 * kl_765[k];

        t_496[k] = -hl_271[k]
                   + f_0 * kl_766[k];

        t_497[k] = -hl_272[k]
                   + f_0 * kl_767[k];

        t_498[k] = -hl_273[k]
                   + f_0 * kl_768[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, t_503, hl_274, hl_275, hl_276, hl_277, \
                         hl_278, kl_769, kl_770, kl_771, kl_772, \
                         kl_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = -hl_274[k]
                   + f_0 * kl_769[k];

        t_500[k] = -hl_275[k]
                   + f_0 * kl_770[k];

        t_501[k] = -hl_276[k]
                   + f_0 * kl_771[k];

        t_502[k] = -hl_277[k]
                   + f_0 * kl_772[k];

        t_503[k] = -hl_278[k]
                   + f_0 * kl_773[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, hl_279, hl_280, hl_281, hl_282, \
                         hl_283, kl_774, kl_775, kl_776, kl_777, \
                         kl_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = -hl_279[k]
                   + f_0 * kl_774[k];

        t_505[k] = -hl_280[k]
                   + f_0 * kl_775[k];

        t_506[k] = -hl_281[k]
                   + f_0 * kl_776[k];

        t_507[k] = -hl_282[k]
                   + f_0 * kl_777[k];

        t_508[k] = -hl_283[k]
                   + f_0 * kl_778[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, t_513, hl_284, hl_285, hl_286, hl_287, \
                         hl_288, kl_779, kl_780, kl_781, kl_782, \
                         kl_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = -hl_284[k]
                   + f_0 * kl_779[k];

        t_510[k] = -hl_285[k]
                   + f_0 * kl_780[k];

        t_511[k] = -hl_286[k]
                   + f_0 * kl_781[k];

        t_512[k] = -hl_287[k]
                   + f_0 * kl_782[k];

        t_513[k] = -hl_288[k]
                   + f_0 * kl_783[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, hl_289, hl_290, hl_291, hl_292, \
                         hl_293, kl_784, kl_785, kl_786, kl_787, \
                         kl_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = -hl_289[k]
                   + f_0 * kl_784[k];

        t_515[k] = -hl_290[k]
                   + f_0 * kl_785[k];

        t_516[k] = -hl_291[k]
                   + f_0 * kl_786[k];

        t_517[k] = -hl_292[k]
                   + f_0 * kl_787[k];

        t_518[k] = -hl_293[k]
                   + f_0 * kl_788[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hl_294 = buffer.data(hl + 294);
    const auto *hl_295 = buffer.data(hl + 295);
    const auto *hl_296 = buffer.data(hl + 296);
    const auto *hl_297 = buffer.data(hl + 297);
    const auto *hl_298 = buffer.data(hl + 298);
    const auto *hl_299 = buffer.data(hl + 299);
    const auto *hl_300 = buffer.data(hl + 300);
    const auto *hl_301 = buffer.data(hl + 301);
    const auto *hl_302 = buffer.data(hl + 302);
    const auto *hl_303 = buffer.data(hl + 303);
    const auto *hl_304 = buffer.data(hl + 304);
    const auto *hl_305 = buffer.data(hl + 305);
    const auto *hl_306 = buffer.data(hl + 306);
    const auto *hl_307 = buffer.data(hl + 307);
    const auto *hl_308 = buffer.data(hl + 308);
    const auto *hl_309 = buffer.data(hl + 309);
    const auto *hl_310 = buffer.data(hl + 310);
    const auto *hl_311 = buffer.data(hl + 311);
    const auto *hl_312 = buffer.data(hl + 312);
    const auto *hl_313 = buffer.data(hl + 313);
    const auto *hl_314 = buffer.data(hl + 314);
    const auto *hl_315 = buffer.data(hl + 315);
    const auto *hl_316 = buffer.data(hl + 316);
    const auto *hl_317 = buffer.data(hl + 317);
    const auto *hl_318 = buffer.data(hl + 318);
    const auto *hl_319 = buffer.data(hl + 319);
    const auto *hl_320 = buffer.data(hl + 320);
    const auto *hl_321 = buffer.data(hl + 321);
    const auto *hl_322 = buffer.data(hl + 322);
    const auto *hl_323 = buffer.data(hl + 323);
    const auto *hl_324 = buffer.data(hl + 324);
    const auto *hl_325 = buffer.data(hl + 325);
    const auto *hl_326 = buffer.data(hl + 326);
    const auto *hl_327 = buffer.data(hl + 327);
    const auto *hl_328 = buffer.data(hl + 328);
    const auto *hl_329 = buffer.data(hl + 329);
    const auto *hl_330 = buffer.data(hl + 330);
    const auto *hl_331 = buffer.data(hl + 331);
    const auto *hl_332 = buffer.data(hl + 332);
    const auto *hl_333 = buffer.data(hl + 333);
    const auto *hl_334 = buffer.data(hl + 334);
    const auto *hl_335 = buffer.data(hl + 335);
    const auto *hl_336 = buffer.data(hl + 336);
    const auto *hl_337 = buffer.data(hl + 337);
    const auto *hl_338 = buffer.data(hl + 338);
    const auto *hl_339 = buffer.data(hl + 339);
    const auto *hl_340 = buffer.data(hl + 340);
    const auto *hl_341 = buffer.data(hl + 341);
    const auto *hl_342 = buffer.data(hl + 342);
    const auto *hl_343 = buffer.data(hl + 343);
    const auto *hl_344 = buffer.data(hl + 344);
    const auto *hl_345 = buffer.data(hl + 345);
    const auto *hl_346 = buffer.data(hl + 346);
    const auto *hl_347 = buffer.data(hl + 347);
    const auto *hl_348 = buffer.data(hl + 348);
    const auto *hl_349 = buffer.data(hl + 349);
    const auto *hl_350 = buffer.data(hl + 350);
    const auto *hl_351 = buffer.data(hl + 351);
    const auto *hl_352 = buffer.data(hl + 352);
    const auto *hl_353 = buffer.data(hl + 353);
    const auto *hl_354 = buffer.data(hl + 354);
    const auto *hl_355 = buffer.data(hl + 355);
    const auto *hl_356 = buffer.data(hl + 356);
    const auto *hl_357 = buffer.data(hl + 357);
    const auto *hl_358 = buffer.data(hl + 358);
    const auto *hl_359 = buffer.data(hl + 359);
    const auto *hl_360 = buffer.data(hl + 360);
    const auto *hl_361 = buffer.data(hl + 361);
    const auto *hl_362 = buffer.data(hl + 362);
    const auto *hl_363 = buffer.data(hl + 363);
    const auto *hl_364 = buffer.data(hl + 364);
    const auto *hl_365 = buffer.data(hl + 365);
    const auto *hl_366 = buffer.data(hl + 366);
    const auto *hl_367 = buffer.data(hl + 367);
    const auto *hl_368 = buffer.data(hl + 368);
    const auto *hl_369 = buffer.data(hl + 369);
    const auto *hl_370 = buffer.data(hl + 370);
    const auto *hl_371 = buffer.data(hl + 371);
    const auto *hl_372 = buffer.data(hl + 372);
    const auto *hl_373 = buffer.data(hl + 373);
    const auto *hl_374 = buffer.data(hl + 374);
    const auto *hl_375 = buffer.data(hl + 375);
    const auto *hl_376 = buffer.data(hl + 376);
    const auto *hl_377 = buffer.data(hl + 377);
    const auto *hl_378 = buffer.data(hl + 378);
    const auto *hl_379 = buffer.data(hl + 379);
    const auto *hl_380 = buffer.data(hl + 380);
    const auto *hl_381 = buffer.data(hl + 381);
    const auto *hl_382 = buffer.data(hl + 382);
    const auto *hl_383 = buffer.data(hl + 383);
    const auto *hl_384 = buffer.data(hl + 384);
    const auto *hl_385 = buffer.data(hl + 385);
    const auto *hl_386 = buffer.data(hl + 386);
    const auto *hl_387 = buffer.data(hl + 387);
    const auto *hl_388 = buffer.data(hl + 388);
    const auto *hl_389 = buffer.data(hl + 389);
    const auto *hl_390 = buffer.data(hl + 390);
    const auto *hl_391 = buffer.data(hl + 391);
    const auto *hl_392 = buffer.data(hl + 392);
    const auto *hl_393 = buffer.data(hl + 393);
    const auto *hl_394 = buffer.data(hl + 394);
    const auto *hl_395 = buffer.data(hl + 395);
    const auto *hl_396 = buffer.data(hl + 396);
    const auto *hl_397 = buffer.data(hl + 397);
    const auto *hl_398 = buffer.data(hl + 398);
    const auto *hl_399 = buffer.data(hl + 399);
    const auto *hl_400 = buffer.data(hl + 400);
    const auto *hl_401 = buffer.data(hl + 401);
    const auto *hl_402 = buffer.data(hl + 402);
    const auto *hl_403 = buffer.data(hl + 403);
    const auto *hl_404 = buffer.data(hl + 404);
    const auto *hl_405 = buffer.data(hl + 405);
    const auto *hl_406 = buffer.data(hl + 406);
    const auto *hl_407 = buffer.data(hl + 407);
    const auto *hl_408 = buffer.data(hl + 408);
    const auto *hl_409 = buffer.data(hl + 409);
    const auto *hl_410 = buffer.data(hl + 410);
    const auto *hl_411 = buffer.data(hl + 411);
    const auto *hl_412 = buffer.data(hl + 412);
    const auto *hl_413 = buffer.data(hl + 413);
    const auto *hl_414 = buffer.data(hl + 414);
    const auto *hl_415 = buffer.data(hl + 415);
    const auto *hl_416 = buffer.data(hl + 416);
    const auto *hl_417 = buffer.data(hl + 417);
    const auto *hl_418 = buffer.data(hl + 418);
    const auto *hl_419 = buffer.data(hl + 419);
    const auto *hl_420 = buffer.data(hl + 420);
    const auto *hl_421 = buffer.data(hl + 421);
    const auto *hl_422 = buffer.data(hl + 422);
    const auto *hl_423 = buffer.data(hl + 423);
    const auto *hl_424 = buffer.data(hl + 424);
    const auto *hl_425 = buffer.data(hl + 425);
    const auto *hl_426 = buffer.data(hl + 426);
    const auto *hl_427 = buffer.data(hl + 427);
    const auto *hl_428 = buffer.data(hl + 428);
    const auto *hl_429 = buffer.data(hl + 429);
    const auto *hl_430 = buffer.data(hl + 430);
    const auto *hl_431 = buffer.data(hl + 431);
    const auto *hl_432 = buffer.data(hl + 432);
    const auto *hl_433 = buffer.data(hl + 433);
    const auto *hl_434 = buffer.data(hl + 434);
    const auto *hl_435 = buffer.data(hl + 435);
    const auto *hl_436 = buffer.data(hl + 436);
    const auto *hl_437 = buffer.data(hl + 437);
    const auto *hl_438 = buffer.data(hl + 438);
    const auto *hl_439 = buffer.data(hl + 439);
    const auto *hl_440 = buffer.data(hl + 440);
    const auto *hl_441 = buffer.data(hl + 441);
    const auto *hl_442 = buffer.data(hl + 442);
    const auto *hl_443 = buffer.data(hl + 443);

    const auto *kl_789 = buffer.data(kl + 789);
    const auto *kl_790 = buffer.data(kl + 790);
    const auto *kl_791 = buffer.data(kl + 791);
    const auto *kl_792 = buffer.data(kl + 792);
    const auto *kl_793 = buffer.data(kl + 793);
    const auto *kl_794 = buffer.data(kl + 794);
    const auto *kl_795 = buffer.data(kl + 795);
    const auto *kl_796 = buffer.data(kl + 796);
    const auto *kl_797 = buffer.data(kl + 797);
    const auto *kl_798 = buffer.data(kl + 798);
    const auto *kl_799 = buffer.data(kl + 799);
    const auto *kl_800 = buffer.data(kl + 800);
    const auto *kl_801 = buffer.data(kl + 801);
    const auto *kl_802 = buffer.data(kl + 802);
    const auto *kl_803 = buffer.data(kl + 803);
    const auto *kl_804 = buffer.data(kl + 804);
    const auto *kl_805 = buffer.data(kl + 805);
    const auto *kl_806 = buffer.data(kl + 806);
    const auto *kl_807 = buffer.data(kl + 807);
    const auto *kl_808 = buffer.data(kl + 808);
    const auto *kl_809 = buffer.data(kl + 809);
    const auto *kl_810 = buffer.data(kl + 810);
    const auto *kl_811 = buffer.data(kl + 811);
    const auto *kl_812 = buffer.data(kl + 812);
    const auto *kl_813 = buffer.data(kl + 813);
    const auto *kl_814 = buffer.data(kl + 814);
    const auto *kl_815 = buffer.data(kl + 815);
    const auto *kl_816 = buffer.data(kl + 816);
    const auto *kl_817 = buffer.data(kl + 817);
    const auto *kl_818 = buffer.data(kl + 818);
    const auto *kl_819 = buffer.data(kl + 819);
    const auto *kl_820 = buffer.data(kl + 820);
    const auto *kl_821 = buffer.data(kl + 821);
    const auto *kl_822 = buffer.data(kl + 822);
    const auto *kl_823 = buffer.data(kl + 823);
    const auto *kl_824 = buffer.data(kl + 824);
    const auto *kl_825 = buffer.data(kl + 825);
    const auto *kl_826 = buffer.data(kl + 826);
    const auto *kl_827 = buffer.data(kl + 827);
    const auto *kl_828 = buffer.data(kl + 828);
    const auto *kl_829 = buffer.data(kl + 829);
    const auto *kl_830 = buffer.data(kl + 830);
    const auto *kl_831 = buffer.data(kl + 831);
    const auto *kl_832 = buffer.data(kl + 832);
    const auto *kl_833 = buffer.data(kl + 833);
    const auto *kl_834 = buffer.data(kl + 834);
    const auto *kl_835 = buffer.data(kl + 835);
    const auto *kl_836 = buffer.data(kl + 836);
    const auto *kl_837 = buffer.data(kl + 837);
    const auto *kl_838 = buffer.data(kl + 838);
    const auto *kl_839 = buffer.data(kl + 839);
    const auto *kl_840 = buffer.data(kl + 840);
    const auto *kl_841 = buffer.data(kl + 841);
    const auto *kl_842 = buffer.data(kl + 842);
    const auto *kl_843 = buffer.data(kl + 843);
    const auto *kl_844 = buffer.data(kl + 844);
    const auto *kl_845 = buffer.data(kl + 845);
    const auto *kl_846 = buffer.data(kl + 846);
    const auto *kl_847 = buffer.data(kl + 847);
    const auto *kl_848 = buffer.data(kl + 848);
    const auto *kl_849 = buffer.data(kl + 849);
    const auto *kl_850 = buffer.data(kl + 850);
    const auto *kl_851 = buffer.data(kl + 851);
    const auto *kl_852 = buffer.data(kl + 852);
    const auto *kl_853 = buffer.data(kl + 853);
    const auto *kl_854 = buffer.data(kl + 854);
    const auto *kl_855 = buffer.data(kl + 855);
    const auto *kl_856 = buffer.data(kl + 856);
    const auto *kl_857 = buffer.data(kl + 857);
    const auto *kl_858 = buffer.data(kl + 858);
    const auto *kl_859 = buffer.data(kl + 859);
    const auto *kl_860 = buffer.data(kl + 860);
    const auto *kl_861 = buffer.data(kl + 861);
    const auto *kl_862 = buffer.data(kl + 862);
    const auto *kl_863 = buffer.data(kl + 863);
    const auto *kl_864 = buffer.data(kl + 864);
    const auto *kl_865 = buffer.data(kl + 865);
    const auto *kl_866 = buffer.data(kl + 866);
    const auto *kl_867 = buffer.data(kl + 867);
    const auto *kl_868 = buffer.data(kl + 868);
    const auto *kl_869 = buffer.data(kl + 869);
    const auto *kl_870 = buffer.data(kl + 870);
    const auto *kl_871 = buffer.data(kl + 871);
    const auto *kl_872 = buffer.data(kl + 872);
    const auto *kl_873 = buffer.data(kl + 873);
    const auto *kl_874 = buffer.data(kl + 874);
    const auto *kl_875 = buffer.data(kl + 875);
    const auto *kl_876 = buffer.data(kl + 876);
    const auto *kl_877 = buffer.data(kl + 877);
    const auto *kl_878 = buffer.data(kl + 878);
    const auto *kl_879 = buffer.data(kl + 879);
    const auto *kl_880 = buffer.data(kl + 880);
    const auto *kl_881 = buffer.data(kl + 881);
    const auto *kl_882 = buffer.data(kl + 882);
    const auto *kl_883 = buffer.data(kl + 883);
    const auto *kl_884 = buffer.data(kl + 884);
    const auto *kl_885 = buffer.data(kl + 885);
    const auto *kl_886 = buffer.data(kl + 886);
    const auto *kl_887 = buffer.data(kl + 887);
    const auto *kl_888 = buffer.data(kl + 888);
    const auto *kl_889 = buffer.data(kl + 889);
    const auto *kl_890 = buffer.data(kl + 890);
    const auto *kl_891 = buffer.data(kl + 891);
    const auto *kl_892 = buffer.data(kl + 892);
    const auto *kl_893 = buffer.data(kl + 893);
    const auto *kl_894 = buffer.data(kl + 894);
    const auto *kl_895 = buffer.data(kl + 895);
    const auto *kl_896 = buffer.data(kl + 896);
    const auto *kl_897 = buffer.data(kl + 897);
    const auto *kl_898 = buffer.data(kl + 898);
    const auto *kl_899 = buffer.data(kl + 899);
    const auto *kl_900 = buffer.data(kl + 900);
    const auto *kl_901 = buffer.data(kl + 901);
    const auto *kl_902 = buffer.data(kl + 902);
    const auto *kl_903 = buffer.data(kl + 903);
    const auto *kl_904 = buffer.data(kl + 904);
    const auto *kl_905 = buffer.data(kl + 905);
    const auto *kl_906 = buffer.data(kl + 906);
    const auto *kl_907 = buffer.data(kl + 907);
    const auto *kl_908 = buffer.data(kl + 908);
    const auto *kl_909 = buffer.data(kl + 909);
    const auto *kl_910 = buffer.data(kl + 910);
    const auto *kl_911 = buffer.data(kl + 911);
    const auto *kl_912 = buffer.data(kl + 912);
    const auto *kl_913 = buffer.data(kl + 913);
    const auto *kl_914 = buffer.data(kl + 914);
    const auto *kl_915 = buffer.data(kl + 915);
    const auto *kl_916 = buffer.data(kl + 916);
    const auto *kl_917 = buffer.data(kl + 917);
    const auto *kl_918 = buffer.data(kl + 918);
    const auto *kl_919 = buffer.data(kl + 919);
    const auto *kl_920 = buffer.data(kl + 920);
    const auto *kl_921 = buffer.data(kl + 921);
    const auto *kl_922 = buffer.data(kl + 922);
    const auto *kl_923 = buffer.data(kl + 923);
    const auto *kl_924 = buffer.data(kl + 924);
    const auto *kl_925 = buffer.data(kl + 925);
    const auto *kl_926 = buffer.data(kl + 926);
    const auto *kl_927 = buffer.data(kl + 927);
    const auto *kl_928 = buffer.data(kl + 928);
    const auto *kl_929 = buffer.data(kl + 929);
    const auto *kl_930 = buffer.data(kl + 930);
    const auto *kl_931 = buffer.data(kl + 931);
    const auto *kl_932 = buffer.data(kl + 932);
    const auto *kl_933 = buffer.data(kl + 933);
    const auto *kl_934 = buffer.data(kl + 934);
    const auto *kl_935 = buffer.data(kl + 935);
    const auto *kl_936 = buffer.data(kl + 936);
    const auto *kl_937 = buffer.data(kl + 937);
    const auto *kl_938 = buffer.data(kl + 938);

#pragma omp simd aligned(t_519, t_520, t_521, t_522, t_523, hl_294, hl_295, hl_296, hl_297, \
                         hl_298, kl_789, kl_790, kl_791, kl_792, \
                         kl_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = -hl_294[k]
                   + f_0 * kl_789[k];

        t_520[k] = -hl_295[k]
                   + f_0 * kl_790[k];

        t_521[k] = -hl_296[k]
                   + f_0 * kl_791[k];

        t_522[k] = -hl_297[k]
                   + f_0 * kl_792[k];

        t_523[k] = -hl_298[k]
                   + f_0 * kl_793[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, hl_299, hl_300, hl_301, hl_302, \
                         hl_303, kl_794, kl_795, kl_796, kl_797, \
                         kl_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = -hl_299[k]
                   + f_0 * kl_794[k];

        t_525[k] = -hl_300[k]
                   + f_0 * kl_795[k];

        t_526[k] = -hl_301[k]
                   + f_0 * kl_796[k];

        t_527[k] = -hl_302[k]
                   + f_0 * kl_797[k];

        t_528[k] = -hl_303[k]
                   + f_0 * kl_798[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, t_533, hl_304, hl_305, hl_306, hl_307, \
                         hl_308, kl_799, kl_800, kl_801, kl_802, \
                         kl_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = -hl_304[k]
                   + f_0 * kl_799[k];

        t_530[k] = -hl_305[k]
                   + f_0 * kl_800[k];

        t_531[k] = -hl_306[k]
                   + f_0 * kl_801[k];

        t_532[k] = -hl_307[k]
                   + f_0 * kl_802[k];

        t_533[k] = -hl_308[k]
                   + f_0 * kl_803[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, t_538, hl_309, hl_310, hl_311, hl_312, \
                         hl_313, kl_804, kl_805, kl_806, kl_807, \
                         kl_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = -hl_309[k]
                   + f_0 * kl_804[k];

        t_535[k] = -hl_310[k]
                   + f_0 * kl_805[k];

        t_536[k] = -hl_311[k]
                   + f_0 * kl_806[k];

        t_537[k] = -hl_312[k]
                   + f_0 * kl_807[k];

        t_538[k] = -hl_313[k]
                   + f_0 * kl_808[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, t_542, t_543, hl_314, hl_315, hl_316, hl_317, \
                         hl_318, kl_809, kl_810, kl_811, kl_812, \
                         kl_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = -hl_314[k]
                   + f_0 * kl_809[k];

        t_540[k] = -2.0 * hl_315[k]
                   + f_0 * kl_810[k];

        t_541[k] = -2.0 * hl_316[k]
                   + f_0 * kl_811[k];

        t_542[k] = -2.0 * hl_317[k]
                   + f_0 * kl_812[k];

        t_543[k] = -2.0 * hl_318[k]
                   + f_0 * kl_813[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, t_548, hl_319, hl_320, hl_321, hl_322, \
                         hl_323, kl_814, kl_815, kl_816, kl_817, \
                         kl_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = -2.0 * hl_319[k]
                   + f_0 * kl_814[k];

        t_545[k] = -2.0 * hl_320[k]
                   + f_0 * kl_815[k];

        t_546[k] = -2.0 * hl_321[k]
                   + f_0 * kl_816[k];

        t_547[k] = -2.0 * hl_322[k]
                   + f_0 * kl_817[k];

        t_548[k] = -2.0 * hl_323[k]
                   + f_0 * kl_818[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, t_553, hl_324, hl_325, hl_326, hl_327, \
                         hl_328, kl_819, kl_820, kl_821, kl_822, \
                         kl_823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = -2.0 * hl_324[k]
                   + f_0 * kl_819[k];

        t_550[k] = -2.0 * hl_325[k]
                   + f_0 * kl_820[k];

        t_551[k] = -2.0 * hl_326[k]
                   + f_0 * kl_821[k];

        t_552[k] = -2.0 * hl_327[k]
                   + f_0 * kl_822[k];

        t_553[k] = -2.0 * hl_328[k]
                   + f_0 * kl_823[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, t_558, hl_329, hl_330, hl_331, hl_332, \
                         hl_333, kl_824, kl_825, kl_826, kl_827, \
                         kl_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = -2.0 * hl_329[k]
                   + f_0 * kl_824[k];

        t_555[k] = -2.0 * hl_330[k]
                   + f_0 * kl_825[k];

        t_556[k] = -2.0 * hl_331[k]
                   + f_0 * kl_826[k];

        t_557[k] = -2.0 * hl_332[k]
                   + f_0 * kl_827[k];

        t_558[k] = -2.0 * hl_333[k]
                   + f_0 * kl_828[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, t_563, hl_334, hl_335, hl_336, hl_337, \
                         hl_338, kl_829, kl_830, kl_831, kl_832, \
                         kl_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = -2.0 * hl_334[k]
                   + f_0 * kl_829[k];

        t_560[k] = -2.0 * hl_335[k]
                   + f_0 * kl_830[k];

        t_561[k] = -2.0 * hl_336[k]
                   + f_0 * kl_831[k];

        t_562[k] = -2.0 * hl_337[k]
                   + f_0 * kl_832[k];

        t_563[k] = -2.0 * hl_338[k]
                   + f_0 * kl_833[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, t_568, hl_339, hl_340, hl_341, hl_342, \
                         hl_343, kl_834, kl_835, kl_836, kl_837, \
                         kl_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = -2.0 * hl_339[k]
                   + f_0 * kl_834[k];

        t_565[k] = -2.0 * hl_340[k]
                   + f_0 * kl_835[k];

        t_566[k] = -2.0 * hl_341[k]
                   + f_0 * kl_836[k];

        t_567[k] = -2.0 * hl_342[k]
                   + f_0 * kl_837[k];

        t_568[k] = -2.0 * hl_343[k]
                   + f_0 * kl_838[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, t_572, t_573, hl_344, hl_345, hl_346, hl_347, \
                         hl_348, kl_839, kl_840, kl_841, kl_842, \
                         kl_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = -2.0 * hl_344[k]
                   + f_0 * kl_839[k];

        t_570[k] = -2.0 * hl_345[k]
                   + f_0 * kl_840[k];

        t_571[k] = -2.0 * hl_346[k]
                   + f_0 * kl_841[k];

        t_572[k] = -2.0 * hl_347[k]
                   + f_0 * kl_842[k];

        t_573[k] = -2.0 * hl_348[k]
                   + f_0 * kl_843[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, hl_349, hl_350, hl_351, hl_352, \
                         hl_353, kl_844, kl_845, kl_846, kl_847, \
                         kl_848 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = -2.0 * hl_349[k]
                   + f_0 * kl_844[k];

        t_575[k] = -2.0 * hl_350[k]
                   + f_0 * kl_845[k];

        t_576[k] = -2.0 * hl_351[k]
                   + f_0 * kl_846[k];

        t_577[k] = -2.0 * hl_352[k]
                   + f_0 * kl_847[k];

        t_578[k] = -2.0 * hl_353[k]
                   + f_0 * kl_848[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, t_582, t_583, hl_354, hl_355, hl_356, hl_357, \
                         hl_358, kl_849, kl_850, kl_851, kl_852, \
                         kl_853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = -2.0 * hl_354[k]
                   + f_0 * kl_849[k];

        t_580[k] = -2.0 * hl_355[k]
                   + f_0 * kl_850[k];

        t_581[k] = -2.0 * hl_356[k]
                   + f_0 * kl_851[k];

        t_582[k] = -2.0 * hl_357[k]
                   + f_0 * kl_852[k];

        t_583[k] = -2.0 * hl_358[k]
                   + f_0 * kl_853[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, t_588, hl_359, hl_360, hl_361, hl_362, \
                         hl_363, kl_854, kl_855, kl_856, kl_857, \
                         kl_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = -2.0 * hl_359[k]
                   + f_0 * kl_854[k];

        t_585[k] = -3.0 * hl_360[k]
                   + f_0 * kl_855[k];

        t_586[k] = -3.0 * hl_361[k]
                   + f_0 * kl_856[k];

        t_587[k] = -3.0 * hl_362[k]
                   + f_0 * kl_857[k];

        t_588[k] = -3.0 * hl_363[k]
                   + f_0 * kl_858[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, t_593, hl_364, hl_365, hl_366, hl_367, \
                         hl_368, kl_859, kl_860, kl_861, kl_862, \
                         kl_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = -3.0 * hl_364[k]
                   + f_0 * kl_859[k];

        t_590[k] = -3.0 * hl_365[k]
                   + f_0 * kl_860[k];

        t_591[k] = -3.0 * hl_366[k]
                   + f_0 * kl_861[k];

        t_592[k] = -3.0 * hl_367[k]
                   + f_0 * kl_862[k];

        t_593[k] = -3.0 * hl_368[k]
                   + f_0 * kl_863[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, t_598, hl_369, hl_370, hl_371, hl_372, \
                         hl_373, kl_864, kl_865, kl_866, kl_867, \
                         kl_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = -3.0 * hl_369[k]
                   + f_0 * kl_864[k];

        t_595[k] = -3.0 * hl_370[k]
                   + f_0 * kl_865[k];

        t_596[k] = -3.0 * hl_371[k]
                   + f_0 * kl_866[k];

        t_597[k] = -3.0 * hl_372[k]
                   + f_0 * kl_867[k];

        t_598[k] = -3.0 * hl_373[k]
                   + f_0 * kl_868[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, hl_374, hl_375, hl_376, hl_377, \
                         hl_378, kl_869, kl_870, kl_871, kl_872, \
                         kl_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = -3.0 * hl_374[k]
                   + f_0 * kl_869[k];

        t_600[k] = -3.0 * hl_375[k]
                   + f_0 * kl_870[k];

        t_601[k] = -3.0 * hl_376[k]
                   + f_0 * kl_871[k];

        t_602[k] = -3.0 * hl_377[k]
                   + f_0 * kl_872[k];

        t_603[k] = -3.0 * hl_378[k]
                   + f_0 * kl_873[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, hl_379, hl_380, hl_381, hl_382, \
                         hl_383, kl_874, kl_875, kl_876, kl_877, \
                         kl_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = -3.0 * hl_379[k]
                   + f_0 * kl_874[k];

        t_605[k] = -3.0 * hl_380[k]
                   + f_0 * kl_875[k];

        t_606[k] = -3.0 * hl_381[k]
                   + f_0 * kl_876[k];

        t_607[k] = -3.0 * hl_382[k]
                   + f_0 * kl_877[k];

        t_608[k] = -3.0 * hl_383[k]
                   + f_0 * kl_878[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, t_613, hl_384, hl_385, hl_386, hl_387, \
                         hl_388, kl_879, kl_880, kl_881, kl_882, \
                         kl_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = -3.0 * hl_384[k]
                   + f_0 * kl_879[k];

        t_610[k] = -3.0 * hl_385[k]
                   + f_0 * kl_880[k];

        t_611[k] = -3.0 * hl_386[k]
                   + f_0 * kl_881[k];

        t_612[k] = -3.0 * hl_387[k]
                   + f_0 * kl_882[k];

        t_613[k] = -3.0 * hl_388[k]
                   + f_0 * kl_883[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, hl_389, hl_390, hl_391, hl_392, \
                         hl_393, kl_884, kl_885, kl_886, kl_887, \
                         kl_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = -3.0 * hl_389[k]
                   + f_0 * kl_884[k];

        t_615[k] = -3.0 * hl_390[k]
                   + f_0 * kl_885[k];

        t_616[k] = -3.0 * hl_391[k]
                   + f_0 * kl_886[k];

        t_617[k] = -3.0 * hl_392[k]
                   + f_0 * kl_887[k];

        t_618[k] = -3.0 * hl_393[k]
                   + f_0 * kl_888[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, hl_394, hl_395, hl_396, hl_397, \
                         hl_398, kl_889, kl_890, kl_891, kl_892, \
                         kl_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = -3.0 * hl_394[k]
                   + f_0 * kl_889[k];

        t_620[k] = -3.0 * hl_395[k]
                   + f_0 * kl_890[k];

        t_621[k] = -3.0 * hl_396[k]
                   + f_0 * kl_891[k];

        t_622[k] = -3.0 * hl_397[k]
                   + f_0 * kl_892[k];

        t_623[k] = -3.0 * hl_398[k]
                   + f_0 * kl_893[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, t_628, hl_399, hl_400, hl_401, hl_402, \
                         hl_403, kl_894, kl_895, kl_896, kl_897, \
                         kl_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = -3.0 * hl_399[k]
                   + f_0 * kl_894[k];

        t_625[k] = -3.0 * hl_400[k]
                   + f_0 * kl_895[k];

        t_626[k] = -3.0 * hl_401[k]
                   + f_0 * kl_896[k];

        t_627[k] = -3.0 * hl_402[k]
                   + f_0 * kl_897[k];

        t_628[k] = -3.0 * hl_403[k]
                   + f_0 * kl_898[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, t_633, hl_404, hl_405, hl_406, hl_407, \
                         hl_408, kl_899, kl_900, kl_901, kl_902, \
                         kl_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = -3.0 * hl_404[k]
                   + f_0 * kl_899[k];

        t_630[k] = -4.0 * hl_405[k]
                   + f_0 * kl_900[k];

        t_631[k] = -4.0 * hl_406[k]
                   + f_0 * kl_901[k];

        t_632[k] = -4.0 * hl_407[k]
                   + f_0 * kl_902[k];

        t_633[k] = -4.0 * hl_408[k]
                   + f_0 * kl_903[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, hl_409, hl_410, hl_411, hl_412, \
                         hl_413, kl_904, kl_905, kl_906, kl_907, \
                         kl_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = -4.0 * hl_409[k]
                   + f_0 * kl_904[k];

        t_635[k] = -4.0 * hl_410[k]
                   + f_0 * kl_905[k];

        t_636[k] = -4.0 * hl_411[k]
                   + f_0 * kl_906[k];

        t_637[k] = -4.0 * hl_412[k]
                   + f_0 * kl_907[k];

        t_638[k] = -4.0 * hl_413[k]
                   + f_0 * kl_908[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, t_643, hl_414, hl_415, hl_416, hl_417, \
                         hl_418, kl_909, kl_910, kl_911, kl_912, \
                         kl_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = -4.0 * hl_414[k]
                   + f_0 * kl_909[k];

        t_640[k] = -4.0 * hl_415[k]
                   + f_0 * kl_910[k];

        t_641[k] = -4.0 * hl_416[k]
                   + f_0 * kl_911[k];

        t_642[k] = -4.0 * hl_417[k]
                   + f_0 * kl_912[k];

        t_643[k] = -4.0 * hl_418[k]
                   + f_0 * kl_913[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, t_648, hl_419, hl_420, hl_421, hl_422, \
                         hl_423, kl_914, kl_915, kl_916, kl_917, \
                         kl_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = -4.0 * hl_419[k]
                   + f_0 * kl_914[k];

        t_645[k] = -4.0 * hl_420[k]
                   + f_0 * kl_915[k];

        t_646[k] = -4.0 * hl_421[k]
                   + f_0 * kl_916[k];

        t_647[k] = -4.0 * hl_422[k]
                   + f_0 * kl_917[k];

        t_648[k] = -4.0 * hl_423[k]
                   + f_0 * kl_918[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, hl_424, hl_425, hl_426, hl_427, \
                         hl_428, kl_919, kl_920, kl_921, kl_922, \
                         kl_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = -4.0 * hl_424[k]
                   + f_0 * kl_919[k];

        t_650[k] = -4.0 * hl_425[k]
                   + f_0 * kl_920[k];

        t_651[k] = -4.0 * hl_426[k]
                   + f_0 * kl_921[k];

        t_652[k] = -4.0 * hl_427[k]
                   + f_0 * kl_922[k];

        t_653[k] = -4.0 * hl_428[k]
                   + f_0 * kl_923[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, t_657, t_658, hl_429, hl_430, hl_431, hl_432, \
                         hl_433, kl_924, kl_925, kl_926, kl_927, \
                         kl_928 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = -4.0 * hl_429[k]
                   + f_0 * kl_924[k];

        t_655[k] = -4.0 * hl_430[k]
                   + f_0 * kl_925[k];

        t_656[k] = -4.0 * hl_431[k]
                   + f_0 * kl_926[k];

        t_657[k] = -4.0 * hl_432[k]
                   + f_0 * kl_927[k];

        t_658[k] = -4.0 * hl_433[k]
                   + f_0 * kl_928[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, hl_434, hl_435, hl_436, hl_437, \
                         hl_438, kl_929, kl_930, kl_931, kl_932, \
                         kl_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = -4.0 * hl_434[k]
                   + f_0 * kl_929[k];

        t_660[k] = -4.0 * hl_435[k]
                   + f_0 * kl_930[k];

        t_661[k] = -4.0 * hl_436[k]
                   + f_0 * kl_931[k];

        t_662[k] = -4.0 * hl_437[k]
                   + f_0 * kl_932[k];

        t_663[k] = -4.0 * hl_438[k]
                   + f_0 * kl_933[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, t_668, hl_439, hl_440, hl_441, hl_442, \
                         hl_443, kl_934, kl_935, kl_936, kl_937, \
                         kl_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = -4.0 * hl_439[k]
                   + f_0 * kl_934[k];

        t_665[k] = -4.0 * hl_440[k]
                   + f_0 * kl_935[k];

        t_666[k] = -4.0 * hl_441[k]
                   + f_0 * kl_936[k];

        t_667[k] = -4.0 * hl_442[k]
                   + f_0 * kl_937[k];

        t_668[k] = -4.0 * hl_443[k]
                   + f_0 * kl_938[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hl_444 = buffer.data(hl + 444);
    const auto *hl_445 = buffer.data(hl + 445);
    const auto *hl_446 = buffer.data(hl + 446);
    const auto *hl_447 = buffer.data(hl + 447);
    const auto *hl_448 = buffer.data(hl + 448);
    const auto *hl_449 = buffer.data(hl + 449);
    const auto *hl_450 = buffer.data(hl + 450);
    const auto *hl_451 = buffer.data(hl + 451);
    const auto *hl_452 = buffer.data(hl + 452);
    const auto *hl_453 = buffer.data(hl + 453);
    const auto *hl_454 = buffer.data(hl + 454);
    const auto *hl_455 = buffer.data(hl + 455);
    const auto *hl_456 = buffer.data(hl + 456);
    const auto *hl_457 = buffer.data(hl + 457);
    const auto *hl_458 = buffer.data(hl + 458);
    const auto *hl_459 = buffer.data(hl + 459);
    const auto *hl_460 = buffer.data(hl + 460);
    const auto *hl_461 = buffer.data(hl + 461);
    const auto *hl_462 = buffer.data(hl + 462);
    const auto *hl_463 = buffer.data(hl + 463);
    const auto *hl_464 = buffer.data(hl + 464);
    const auto *hl_465 = buffer.data(hl + 465);
    const auto *hl_466 = buffer.data(hl + 466);
    const auto *hl_467 = buffer.data(hl + 467);
    const auto *hl_468 = buffer.data(hl + 468);
    const auto *hl_469 = buffer.data(hl + 469);
    const auto *hl_470 = buffer.data(hl + 470);
    const auto *hl_471 = buffer.data(hl + 471);
    const auto *hl_472 = buffer.data(hl + 472);
    const auto *hl_473 = buffer.data(hl + 473);
    const auto *hl_474 = buffer.data(hl + 474);
    const auto *hl_475 = buffer.data(hl + 475);
    const auto *hl_476 = buffer.data(hl + 476);
    const auto *hl_477 = buffer.data(hl + 477);
    const auto *hl_478 = buffer.data(hl + 478);
    const auto *hl_479 = buffer.data(hl + 479);
    const auto *hl_480 = buffer.data(hl + 480);
    const auto *hl_481 = buffer.data(hl + 481);
    const auto *hl_482 = buffer.data(hl + 482);
    const auto *hl_483 = buffer.data(hl + 483);
    const auto *hl_484 = buffer.data(hl + 484);
    const auto *hl_485 = buffer.data(hl + 485);
    const auto *hl_486 = buffer.data(hl + 486);
    const auto *hl_487 = buffer.data(hl + 487);
    const auto *hl_488 = buffer.data(hl + 488);
    const auto *hl_489 = buffer.data(hl + 489);
    const auto *hl_490 = buffer.data(hl + 490);
    const auto *hl_491 = buffer.data(hl + 491);
    const auto *hl_492 = buffer.data(hl + 492);
    const auto *hl_493 = buffer.data(hl + 493);
    const auto *hl_494 = buffer.data(hl + 494);
    const auto *hl_495 = buffer.data(hl + 495);
    const auto *hl_496 = buffer.data(hl + 496);
    const auto *hl_497 = buffer.data(hl + 497);
    const auto *hl_498 = buffer.data(hl + 498);
    const auto *hl_499 = buffer.data(hl + 499);
    const auto *hl_500 = buffer.data(hl + 500);
    const auto *hl_501 = buffer.data(hl + 501);
    const auto *hl_502 = buffer.data(hl + 502);
    const auto *hl_503 = buffer.data(hl + 503);
    const auto *hl_504 = buffer.data(hl + 504);
    const auto *hl_505 = buffer.data(hl + 505);
    const auto *hl_506 = buffer.data(hl + 506);
    const auto *hl_507 = buffer.data(hl + 507);
    const auto *hl_508 = buffer.data(hl + 508);
    const auto *hl_509 = buffer.data(hl + 509);
    const auto *hl_510 = buffer.data(hl + 510);
    const auto *hl_511 = buffer.data(hl + 511);
    const auto *hl_512 = buffer.data(hl + 512);
    const auto *hl_513 = buffer.data(hl + 513);
    const auto *hl_514 = buffer.data(hl + 514);
    const auto *hl_515 = buffer.data(hl + 515);
    const auto *hl_516 = buffer.data(hl + 516);
    const auto *hl_517 = buffer.data(hl + 517);
    const auto *hl_518 = buffer.data(hl + 518);
    const auto *hl_519 = buffer.data(hl + 519);
    const auto *hl_520 = buffer.data(hl + 520);
    const auto *hl_521 = buffer.data(hl + 521);
    const auto *hl_522 = buffer.data(hl + 522);
    const auto *hl_523 = buffer.data(hl + 523);
    const auto *hl_524 = buffer.data(hl + 524);
    const auto *hl_525 = buffer.data(hl + 525);
    const auto *hl_526 = buffer.data(hl + 526);
    const auto *hl_527 = buffer.data(hl + 527);
    const auto *hl_528 = buffer.data(hl + 528);
    const auto *hl_529 = buffer.data(hl + 529);
    const auto *hl_530 = buffer.data(hl + 530);
    const auto *hl_531 = buffer.data(hl + 531);
    const auto *hl_532 = buffer.data(hl + 532);
    const auto *hl_533 = buffer.data(hl + 533);
    const auto *hl_534 = buffer.data(hl + 534);
    const auto *hl_535 = buffer.data(hl + 535);
    const auto *hl_536 = buffer.data(hl + 536);
    const auto *hl_537 = buffer.data(hl + 537);
    const auto *hl_538 = buffer.data(hl + 538);
    const auto *hl_539 = buffer.data(hl + 539);
    const auto *hl_540 = buffer.data(hl + 540);
    const auto *hl_541 = buffer.data(hl + 541);
    const auto *hl_542 = buffer.data(hl + 542);
    const auto *hl_543 = buffer.data(hl + 543);
    const auto *hl_544 = buffer.data(hl + 544);
    const auto *hl_545 = buffer.data(hl + 545);
    const auto *hl_546 = buffer.data(hl + 546);
    const auto *hl_547 = buffer.data(hl + 547);
    const auto *hl_548 = buffer.data(hl + 548);
    const auto *hl_549 = buffer.data(hl + 549);
    const auto *hl_550 = buffer.data(hl + 550);
    const auto *hl_551 = buffer.data(hl + 551);
    const auto *hl_552 = buffer.data(hl + 552);
    const auto *hl_553 = buffer.data(hl + 553);
    const auto *hl_554 = buffer.data(hl + 554);
    const auto *hl_555 = buffer.data(hl + 555);
    const auto *hl_556 = buffer.data(hl + 556);
    const auto *hl_557 = buffer.data(hl + 557);
    const auto *hl_558 = buffer.data(hl + 558);
    const auto *hl_559 = buffer.data(hl + 559);
    const auto *hl_560 = buffer.data(hl + 560);
    const auto *hl_561 = buffer.data(hl + 561);
    const auto *hl_562 = buffer.data(hl + 562);
    const auto *hl_563 = buffer.data(hl + 563);
    const auto *hl_564 = buffer.data(hl + 564);

    const auto *kl_939 = buffer.data(kl + 939);
    const auto *kl_940 = buffer.data(kl + 940);
    const auto *kl_941 = buffer.data(kl + 941);
    const auto *kl_942 = buffer.data(kl + 942);
    const auto *kl_943 = buffer.data(kl + 943);
    const auto *kl_944 = buffer.data(kl + 944);
    const auto *kl_990 = buffer.data(kl + 990);
    const auto *kl_991 = buffer.data(kl + 991);
    const auto *kl_992 = buffer.data(kl + 992);
    const auto *kl_993 = buffer.data(kl + 993);
    const auto *kl_994 = buffer.data(kl + 994);
    const auto *kl_995 = buffer.data(kl + 995);
    const auto *kl_996 = buffer.data(kl + 996);
    const auto *kl_997 = buffer.data(kl + 997);
    const auto *kl_998 = buffer.data(kl + 998);
    const auto *kl_999 = buffer.data(kl + 999);
    const auto *kl_1000 = buffer.data(kl + 1000);
    const auto *kl_1001 = buffer.data(kl + 1001);
    const auto *kl_1002 = buffer.data(kl + 1002);
    const auto *kl_1003 = buffer.data(kl + 1003);
    const auto *kl_1004 = buffer.data(kl + 1004);
    const auto *kl_1005 = buffer.data(kl + 1005);
    const auto *kl_1006 = buffer.data(kl + 1006);
    const auto *kl_1007 = buffer.data(kl + 1007);
    const auto *kl_1008 = buffer.data(kl + 1008);
    const auto *kl_1009 = buffer.data(kl + 1009);
    const auto *kl_1010 = buffer.data(kl + 1010);
    const auto *kl_1011 = buffer.data(kl + 1011);
    const auto *kl_1012 = buffer.data(kl + 1012);
    const auto *kl_1013 = buffer.data(kl + 1013);
    const auto *kl_1014 = buffer.data(kl + 1014);
    const auto *kl_1015 = buffer.data(kl + 1015);
    const auto *kl_1016 = buffer.data(kl + 1016);
    const auto *kl_1017 = buffer.data(kl + 1017);
    const auto *kl_1018 = buffer.data(kl + 1018);
    const auto *kl_1019 = buffer.data(kl + 1019);
    const auto *kl_1020 = buffer.data(kl + 1020);
    const auto *kl_1021 = buffer.data(kl + 1021);
    const auto *kl_1022 = buffer.data(kl + 1022);
    const auto *kl_1023 = buffer.data(kl + 1023);
    const auto *kl_1024 = buffer.data(kl + 1024);
    const auto *kl_1025 = buffer.data(kl + 1025);
    const auto *kl_1026 = buffer.data(kl + 1026);
    const auto *kl_1027 = buffer.data(kl + 1027);
    const auto *kl_1028 = buffer.data(kl + 1028);
    const auto *kl_1029 = buffer.data(kl + 1029);
    const auto *kl_1030 = buffer.data(kl + 1030);
    const auto *kl_1031 = buffer.data(kl + 1031);
    const auto *kl_1032 = buffer.data(kl + 1032);
    const auto *kl_1033 = buffer.data(kl + 1033);
    const auto *kl_1034 = buffer.data(kl + 1034);
    const auto *kl_1035 = buffer.data(kl + 1035);
    const auto *kl_1036 = buffer.data(kl + 1036);
    const auto *kl_1037 = buffer.data(kl + 1037);
    const auto *kl_1038 = buffer.data(kl + 1038);
    const auto *kl_1039 = buffer.data(kl + 1039);
    const auto *kl_1040 = buffer.data(kl + 1040);
    const auto *kl_1041 = buffer.data(kl + 1041);
    const auto *kl_1042 = buffer.data(kl + 1042);
    const auto *kl_1043 = buffer.data(kl + 1043);
    const auto *kl_1044 = buffer.data(kl + 1044);
    const auto *kl_1045 = buffer.data(kl + 1045);
    const auto *kl_1046 = buffer.data(kl + 1046);
    const auto *kl_1047 = buffer.data(kl + 1047);
    const auto *kl_1048 = buffer.data(kl + 1048);
    const auto *kl_1049 = buffer.data(kl + 1049);
    const auto *kl_1050 = buffer.data(kl + 1050);
    const auto *kl_1051 = buffer.data(kl + 1051);
    const auto *kl_1052 = buffer.data(kl + 1052);
    const auto *kl_1053 = buffer.data(kl + 1053);
    const auto *kl_1054 = buffer.data(kl + 1054);
    const auto *kl_1055 = buffer.data(kl + 1055);
    const auto *kl_1056 = buffer.data(kl + 1056);
    const auto *kl_1057 = buffer.data(kl + 1057);
    const auto *kl_1058 = buffer.data(kl + 1058);
    const auto *kl_1059 = buffer.data(kl + 1059);
    const auto *kl_1060 = buffer.data(kl + 1060);
    const auto *kl_1061 = buffer.data(kl + 1061);
    const auto *kl_1062 = buffer.data(kl + 1062);
    const auto *kl_1063 = buffer.data(kl + 1063);
    const auto *kl_1064 = buffer.data(kl + 1064);
    const auto *kl_1065 = buffer.data(kl + 1065);
    const auto *kl_1066 = buffer.data(kl + 1066);
    const auto *kl_1067 = buffer.data(kl + 1067);
    const auto *kl_1068 = buffer.data(kl + 1068);
    const auto *kl_1069 = buffer.data(kl + 1069);
    const auto *kl_1070 = buffer.data(kl + 1070);
    const auto *kl_1071 = buffer.data(kl + 1071);
    const auto *kl_1072 = buffer.data(kl + 1072);
    const auto *kl_1073 = buffer.data(kl + 1073);
    const auto *kl_1074 = buffer.data(kl + 1074);
    const auto *kl_1075 = buffer.data(kl + 1075);
    const auto *kl_1076 = buffer.data(kl + 1076);
    const auto *kl_1077 = buffer.data(kl + 1077);
    const auto *kl_1078 = buffer.data(kl + 1078);
    const auto *kl_1079 = buffer.data(kl + 1079);
    const auto *kl_1080 = buffer.data(kl + 1080);
    const auto *kl_1081 = buffer.data(kl + 1081);
    const auto *kl_1082 = buffer.data(kl + 1082);
    const auto *kl_1083 = buffer.data(kl + 1083);
    const auto *kl_1084 = buffer.data(kl + 1084);
    const auto *kl_1085 = buffer.data(kl + 1085);
    const auto *kl_1086 = buffer.data(kl + 1086);
    const auto *kl_1087 = buffer.data(kl + 1087);
    const auto *kl_1088 = buffer.data(kl + 1088);
    const auto *kl_1089 = buffer.data(kl + 1089);
    const auto *kl_1090 = buffer.data(kl + 1090);
    const auto *kl_1091 = buffer.data(kl + 1091);
    const auto *kl_1092 = buffer.data(kl + 1092);
    const auto *kl_1093 = buffer.data(kl + 1093);
    const auto *kl_1094 = buffer.data(kl + 1094);
    const auto *kl_1095 = buffer.data(kl + 1095);
    const auto *kl_1096 = buffer.data(kl + 1096);
    const auto *kl_1097 = buffer.data(kl + 1097);
    const auto *kl_1098 = buffer.data(kl + 1098);
    const auto *kl_1099 = buffer.data(kl + 1099);
    const auto *kl_1100 = buffer.data(kl + 1100);
    const auto *kl_1101 = buffer.data(kl + 1101);
    const auto *kl_1102 = buffer.data(kl + 1102);
    const auto *kl_1103 = buffer.data(kl + 1103);
    const auto *kl_1104 = buffer.data(kl + 1104);
    const auto *kl_1105 = buffer.data(kl + 1105);
    const auto *kl_1106 = buffer.data(kl + 1106);
    const auto *kl_1107 = buffer.data(kl + 1107);
    const auto *kl_1108 = buffer.data(kl + 1108);
    const auto *kl_1109 = buffer.data(kl + 1109);
    const auto *kl_1110 = buffer.data(kl + 1110);
    const auto *kl_1111 = buffer.data(kl + 1111);
    const auto *kl_1112 = buffer.data(kl + 1112);
    const auto *kl_1113 = buffer.data(kl + 1113);
    const auto *kl_1114 = buffer.data(kl + 1114);
    const auto *kl_1115 = buffer.data(kl + 1115);
    const auto *kl_1116 = buffer.data(kl + 1116);
    const auto *kl_1117 = buffer.data(kl + 1117);
    const auto *kl_1118 = buffer.data(kl + 1118);
    const auto *kl_1119 = buffer.data(kl + 1119);
    const auto *kl_1120 = buffer.data(kl + 1120);
    const auto *kl_1121 = buffer.data(kl + 1121);
    const auto *kl_1122 = buffer.data(kl + 1122);
    const auto *kl_1123 = buffer.data(kl + 1123);
    const auto *kl_1124 = buffer.data(kl + 1124);
    const auto *kl_1125 = buffer.data(kl + 1125);
    const auto *kl_1126 = buffer.data(kl + 1126);
    const auto *kl_1127 = buffer.data(kl + 1127);
    const auto *kl_1128 = buffer.data(kl + 1128);
    const auto *kl_1129 = buffer.data(kl + 1129);
    const auto *kl_1130 = buffer.data(kl + 1130);
    const auto *kl_1131 = buffer.data(kl + 1131);
    const auto *kl_1132 = buffer.data(kl + 1132);
    const auto *kl_1133 = buffer.data(kl + 1133);
    const auto *kl_1134 = buffer.data(kl + 1134);
    const auto *kl_1135 = buffer.data(kl + 1135);
    const auto *kl_1136 = buffer.data(kl + 1136);
    const auto *kl_1137 = buffer.data(kl + 1137);
    const auto *kl_1138 = buffer.data(kl + 1138);
    const auto *kl_1139 = buffer.data(kl + 1139);
    const auto *kl_1140 = buffer.data(kl + 1140);
    const auto *kl_1141 = buffer.data(kl + 1141);
    const auto *kl_1142 = buffer.data(kl + 1142);
    const auto *kl_1143 = buffer.data(kl + 1143);
    const auto *kl_1144 = buffer.data(kl + 1144);
    const auto *kl_1145 = buffer.data(kl + 1145);
    const auto *kl_1146 = buffer.data(kl + 1146);
    const auto *kl_1147 = buffer.data(kl + 1147);
    const auto *kl_1148 = buffer.data(kl + 1148);
    const auto *kl_1149 = buffer.data(kl + 1149);

#pragma omp simd aligned(t_669, t_670, t_671, t_672, t_673, hl_444, hl_445, hl_446, hl_447, \
                         hl_448, kl_939, kl_940, kl_941, kl_942, \
                         kl_943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = -4.0 * hl_444[k]
                   + f_0 * kl_939[k];

        t_670[k] = -4.0 * hl_445[k]
                   + f_0 * kl_940[k];

        t_671[k] = -4.0 * hl_446[k]
                   + f_0 * kl_941[k];

        t_672[k] = -4.0 * hl_447[k]
                   + f_0 * kl_942[k];

        t_673[k] = -4.0 * hl_448[k]
                   + f_0 * kl_943[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, t_677, t_678, t_679, t_680, hl_449, kl_944, \
                         kl_990, kl_991, kl_992, kl_993, kl_994, \
                         kl_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = -4.0 * hl_449[k]
                   + f_0 * kl_944[k];

        t_675[k] = f_0 * kl_990[k];

        t_676[k] = f_0 * kl_991[k];

        t_677[k] = f_0 * kl_992[k];

        t_678[k] = f_0 * kl_993[k];

        t_679[k] = f_0 * kl_994[k];

        t_680[k] = f_0 * kl_995[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, t_684, t_685, t_686, t_687, t_688, kl_996, \
                         kl_997, kl_998, kl_999, kl_1000, kl_1001, kl_1002, \
                         kl_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = f_0 * kl_996[k];

        t_682[k] = f_0 * kl_997[k];

        t_683[k] = f_0 * kl_998[k];

        t_684[k] = f_0 * kl_999[k];

        t_685[k] = f_0 * kl_1000[k];

        t_686[k] = f_0 * kl_1001[k];

        t_687[k] = f_0 * kl_1002[k];

        t_688[k] = f_0 * kl_1003[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, t_693, t_694, t_695, t_696, kl_1004, \
                         kl_1005, kl_1006, kl_1007, kl_1008, kl_1009, kl_1010, \
                         kl_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = f_0 * kl_1004[k];

        t_690[k] = f_0 * kl_1005[k];

        t_691[k] = f_0 * kl_1006[k];

        t_692[k] = f_0 * kl_1007[k];

        t_693[k] = f_0 * kl_1008[k];

        t_694[k] = f_0 * kl_1009[k];

        t_695[k] = f_0 * kl_1010[k];

        t_696[k] = f_0 * kl_1011[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, t_701, t_702, t_703, t_704, kl_1012, \
                         kl_1013, kl_1014, kl_1015, kl_1016, kl_1017, kl_1018, \
                         kl_1019 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_0 * kl_1012[k];

        t_698[k] = f_0 * kl_1013[k];

        t_699[k] = f_0 * kl_1014[k];

        t_700[k] = f_0 * kl_1015[k];

        t_701[k] = f_0 * kl_1016[k];

        t_702[k] = f_0 * kl_1017[k];

        t_703[k] = f_0 * kl_1018[k];

        t_704[k] = f_0 * kl_1019[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, t_710, t_711, t_712, kl_1020, \
                         kl_1021, kl_1022, kl_1023, kl_1024, kl_1025, kl_1026, \
                         kl_1027 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_0 * kl_1020[k];

        t_706[k] = f_0 * kl_1021[k];

        t_707[k] = f_0 * kl_1022[k];

        t_708[k] = f_0 * kl_1023[k];

        t_709[k] = f_0 * kl_1024[k];

        t_710[k] = f_0 * kl_1025[k];

        t_711[k] = f_0 * kl_1026[k];

        t_712[k] = f_0 * kl_1027[k];
    }

#pragma omp simd aligned(t_713, t_714, t_715, t_716, t_717, t_718, t_719, kl_1028, kl_1029, \
                         kl_1030, kl_1031, kl_1032, kl_1033, kl_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_713[k] = f_0 * kl_1028[k];

        t_714[k] = f_0 * kl_1029[k];

        t_715[k] = f_0 * kl_1030[k];

        t_716[k] = f_0 * kl_1031[k];

        t_717[k] = f_0 * kl_1032[k];

        t_718[k] = f_0 * kl_1033[k];

        t_719[k] = f_0 * kl_1034[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, hl_450, hl_451, hl_452, hl_453, \
                         hl_454, kl_1035, kl_1036, kl_1037, kl_1038, \
                         kl_1039 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -hl_450[k]
                   + f_0 * kl_1035[k];

        t_721[k] = -hl_451[k]
                   + f_0 * kl_1036[k];

        t_722[k] = -hl_452[k]
                   + f_0 * kl_1037[k];

        t_723[k] = -hl_453[k]
                   + f_0 * kl_1038[k];

        t_724[k] = -hl_454[k]
                   + f_0 * kl_1039[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, hl_455, hl_456, hl_457, hl_458, \
                         hl_459, kl_1040, kl_1041, kl_1042, kl_1043, \
                         kl_1044 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -hl_455[k]
                   + f_0 * kl_1040[k];

        t_726[k] = -hl_456[k]
                   + f_0 * kl_1041[k];

        t_727[k] = -hl_457[k]
                   + f_0 * kl_1042[k];

        t_728[k] = -hl_458[k]
                   + f_0 * kl_1043[k];

        t_729[k] = -hl_459[k]
                   + f_0 * kl_1044[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, hl_460, hl_461, hl_462, hl_463, \
                         hl_464, kl_1045, kl_1046, kl_1047, kl_1048, \
                         kl_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -hl_460[k]
                   + f_0 * kl_1045[k];

        t_731[k] = -hl_461[k]
                   + f_0 * kl_1046[k];

        t_732[k] = -hl_462[k]
                   + f_0 * kl_1047[k];

        t_733[k] = -hl_463[k]
                   + f_0 * kl_1048[k];

        t_734[k] = -hl_464[k]
                   + f_0 * kl_1049[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, hl_465, hl_466, hl_467, hl_468, \
                         hl_469, kl_1050, kl_1051, kl_1052, kl_1053, \
                         kl_1054 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -hl_465[k]
                   + f_0 * kl_1050[k];

        t_736[k] = -hl_466[k]
                   + f_0 * kl_1051[k];

        t_737[k] = -hl_467[k]
                   + f_0 * kl_1052[k];

        t_738[k] = -hl_468[k]
                   + f_0 * kl_1053[k];

        t_739[k] = -hl_469[k]
                   + f_0 * kl_1054[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, hl_470, hl_471, hl_472, hl_473, \
                         hl_474, kl_1055, kl_1056, kl_1057, kl_1058, \
                         kl_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -hl_470[k]
                   + f_0 * kl_1055[k];

        t_741[k] = -hl_471[k]
                   + f_0 * kl_1056[k];

        t_742[k] = -hl_472[k]
                   + f_0 * kl_1057[k];

        t_743[k] = -hl_473[k]
                   + f_0 * kl_1058[k];

        t_744[k] = -hl_474[k]
                   + f_0 * kl_1059[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, hl_475, hl_476, hl_477, hl_478, \
                         hl_479, kl_1060, kl_1061, kl_1062, kl_1063, \
                         kl_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -hl_475[k]
                   + f_0 * kl_1060[k];

        t_746[k] = -hl_476[k]
                   + f_0 * kl_1061[k];

        t_747[k] = -hl_477[k]
                   + f_0 * kl_1062[k];

        t_748[k] = -hl_478[k]
                   + f_0 * kl_1063[k];

        t_749[k] = -hl_479[k]
                   + f_0 * kl_1064[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, hl_480, hl_481, hl_482, hl_483, \
                         hl_484, kl_1065, kl_1066, kl_1067, kl_1068, \
                         kl_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -hl_480[k]
                   + f_0 * kl_1065[k];

        t_751[k] = -hl_481[k]
                   + f_0 * kl_1066[k];

        t_752[k] = -hl_482[k]
                   + f_0 * kl_1067[k];

        t_753[k] = -hl_483[k]
                   + f_0 * kl_1068[k];

        t_754[k] = -hl_484[k]
                   + f_0 * kl_1069[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, hl_485, hl_486, hl_487, hl_488, \
                         hl_489, kl_1070, kl_1071, kl_1072, kl_1073, \
                         kl_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -hl_485[k]
                   + f_0 * kl_1070[k];

        t_756[k] = -hl_486[k]
                   + f_0 * kl_1071[k];

        t_757[k] = -hl_487[k]
                   + f_0 * kl_1072[k];

        t_758[k] = -hl_488[k]
                   + f_0 * kl_1073[k];

        t_759[k] = -hl_489[k]
                   + f_0 * kl_1074[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, t_764, hl_490, hl_491, hl_492, hl_493, \
                         hl_494, kl_1075, kl_1076, kl_1077, kl_1078, \
                         kl_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = -hl_490[k]
                   + f_0 * kl_1075[k];

        t_761[k] = -hl_491[k]
                   + f_0 * kl_1076[k];

        t_762[k] = -hl_492[k]
                   + f_0 * kl_1077[k];

        t_763[k] = -hl_493[k]
                   + f_0 * kl_1078[k];

        t_764[k] = -hl_494[k]
                   + f_0 * kl_1079[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, t_769, hl_495, hl_496, hl_497, hl_498, \
                         hl_499, kl_1080, kl_1081, kl_1082, kl_1083, \
                         kl_1084 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = -2.0 * hl_495[k]
                   + f_0 * kl_1080[k];

        t_766[k] = -2.0 * hl_496[k]
                   + f_0 * kl_1081[k];

        t_767[k] = -2.0 * hl_497[k]
                   + f_0 * kl_1082[k];

        t_768[k] = -2.0 * hl_498[k]
                   + f_0 * kl_1083[k];

        t_769[k] = -2.0 * hl_499[k]
                   + f_0 * kl_1084[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, hl_500, hl_501, hl_502, hl_503, \
                         hl_504, kl_1085, kl_1086, kl_1087, kl_1088, \
                         kl_1089 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = -2.0 * hl_500[k]
                   + f_0 * kl_1085[k];

        t_771[k] = -2.0 * hl_501[k]
                   + f_0 * kl_1086[k];

        t_772[k] = -2.0 * hl_502[k]
                   + f_0 * kl_1087[k];

        t_773[k] = -2.0 * hl_503[k]
                   + f_0 * kl_1088[k];

        t_774[k] = -2.0 * hl_504[k]
                   + f_0 * kl_1089[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, hl_505, hl_506, hl_507, hl_508, \
                         hl_509, kl_1090, kl_1091, kl_1092, kl_1093, \
                         kl_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = -2.0 * hl_505[k]
                   + f_0 * kl_1090[k];

        t_776[k] = -2.0 * hl_506[k]
                   + f_0 * kl_1091[k];

        t_777[k] = -2.0 * hl_507[k]
                   + f_0 * kl_1092[k];

        t_778[k] = -2.0 * hl_508[k]
                   + f_0 * kl_1093[k];

        t_779[k] = -2.0 * hl_509[k]
                   + f_0 * kl_1094[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, t_784, hl_510, hl_511, hl_512, hl_513, \
                         hl_514, kl_1095, kl_1096, kl_1097, kl_1098, \
                         kl_1099 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = -2.0 * hl_510[k]
                   + f_0 * kl_1095[k];

        t_781[k] = -2.0 * hl_511[k]
                   + f_0 * kl_1096[k];

        t_782[k] = -2.0 * hl_512[k]
                   + f_0 * kl_1097[k];

        t_783[k] = -2.0 * hl_513[k]
                   + f_0 * kl_1098[k];

        t_784[k] = -2.0 * hl_514[k]
                   + f_0 * kl_1099[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, t_789, hl_515, hl_516, hl_517, hl_518, \
                         hl_519, kl_1100, kl_1101, kl_1102, kl_1103, \
                         kl_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = -2.0 * hl_515[k]
                   + f_0 * kl_1100[k];

        t_786[k] = -2.0 * hl_516[k]
                   + f_0 * kl_1101[k];

        t_787[k] = -2.0 * hl_517[k]
                   + f_0 * kl_1102[k];

        t_788[k] = -2.0 * hl_518[k]
                   + f_0 * kl_1103[k];

        t_789[k] = -2.0 * hl_519[k]
                   + f_0 * kl_1104[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, t_794, hl_520, hl_521, hl_522, hl_523, \
                         hl_524, kl_1105, kl_1106, kl_1107, kl_1108, \
                         kl_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = -2.0 * hl_520[k]
                   + f_0 * kl_1105[k];

        t_791[k] = -2.0 * hl_521[k]
                   + f_0 * kl_1106[k];

        t_792[k] = -2.0 * hl_522[k]
                   + f_0 * kl_1107[k];

        t_793[k] = -2.0 * hl_523[k]
                   + f_0 * kl_1108[k];

        t_794[k] = -2.0 * hl_524[k]
                   + f_0 * kl_1109[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, hl_525, hl_526, hl_527, hl_528, \
                         hl_529, kl_1110, kl_1111, kl_1112, kl_1113, \
                         kl_1114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = -2.0 * hl_525[k]
                   + f_0 * kl_1110[k];

        t_796[k] = -2.0 * hl_526[k]
                   + f_0 * kl_1111[k];

        t_797[k] = -2.0 * hl_527[k]
                   + f_0 * kl_1112[k];

        t_798[k] = -2.0 * hl_528[k]
                   + f_0 * kl_1113[k];

        t_799[k] = -2.0 * hl_529[k]
                   + f_0 * kl_1114[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, hl_530, hl_531, hl_532, hl_533, \
                         hl_534, kl_1115, kl_1116, kl_1117, kl_1118, \
                         kl_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = -2.0 * hl_530[k]
                   + f_0 * kl_1115[k];

        t_801[k] = -2.0 * hl_531[k]
                   + f_0 * kl_1116[k];

        t_802[k] = -2.0 * hl_532[k]
                   + f_0 * kl_1117[k];

        t_803[k] = -2.0 * hl_533[k]
                   + f_0 * kl_1118[k];

        t_804[k] = -2.0 * hl_534[k]
                   + f_0 * kl_1119[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, hl_535, hl_536, hl_537, hl_538, \
                         hl_539, kl_1120, kl_1121, kl_1122, kl_1123, \
                         kl_1124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = -2.0 * hl_535[k]
                   + f_0 * kl_1120[k];

        t_806[k] = -2.0 * hl_536[k]
                   + f_0 * kl_1121[k];

        t_807[k] = -2.0 * hl_537[k]
                   + f_0 * kl_1122[k];

        t_808[k] = -2.0 * hl_538[k]
                   + f_0 * kl_1123[k];

        t_809[k] = -2.0 * hl_539[k]
                   + f_0 * kl_1124[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, hl_540, hl_541, hl_542, hl_543, \
                         hl_544, kl_1125, kl_1126, kl_1127, kl_1128, \
                         kl_1129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = -3.0 * hl_540[k]
                   + f_0 * kl_1125[k];

        t_811[k] = -3.0 * hl_541[k]
                   + f_0 * kl_1126[k];

        t_812[k] = -3.0 * hl_542[k]
                   + f_0 * kl_1127[k];

        t_813[k] = -3.0 * hl_543[k]
                   + f_0 * kl_1128[k];

        t_814[k] = -3.0 * hl_544[k]
                   + f_0 * kl_1129[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, hl_545, hl_546, hl_547, hl_548, \
                         hl_549, kl_1130, kl_1131, kl_1132, kl_1133, \
                         kl_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = -3.0 * hl_545[k]
                   + f_0 * kl_1130[k];

        t_816[k] = -3.0 * hl_546[k]
                   + f_0 * kl_1131[k];

        t_817[k] = -3.0 * hl_547[k]
                   + f_0 * kl_1132[k];

        t_818[k] = -3.0 * hl_548[k]
                   + f_0 * kl_1133[k];

        t_819[k] = -3.0 * hl_549[k]
                   + f_0 * kl_1134[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, t_824, hl_550, hl_551, hl_552, hl_553, \
                         hl_554, kl_1135, kl_1136, kl_1137, kl_1138, \
                         kl_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = -3.0 * hl_550[k]
                   + f_0 * kl_1135[k];

        t_821[k] = -3.0 * hl_551[k]
                   + f_0 * kl_1136[k];

        t_822[k] = -3.0 * hl_552[k]
                   + f_0 * kl_1137[k];

        t_823[k] = -3.0 * hl_553[k]
                   + f_0 * kl_1138[k];

        t_824[k] = -3.0 * hl_554[k]
                   + f_0 * kl_1139[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, hl_555, hl_556, hl_557, hl_558, \
                         hl_559, kl_1140, kl_1141, kl_1142, kl_1143, \
                         kl_1144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = -3.0 * hl_555[k]
                   + f_0 * kl_1140[k];

        t_826[k] = -3.0 * hl_556[k]
                   + f_0 * kl_1141[k];

        t_827[k] = -3.0 * hl_557[k]
                   + f_0 * kl_1142[k];

        t_828[k] = -3.0 * hl_558[k]
                   + f_0 * kl_1143[k];

        t_829[k] = -3.0 * hl_559[k]
                   + f_0 * kl_1144[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, hl_560, hl_561, hl_562, hl_563, \
                         hl_564, kl_1145, kl_1146, kl_1147, kl_1148, \
                         kl_1149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = -3.0 * hl_560[k]
                   + f_0 * kl_1145[k];

        t_831[k] = -3.0 * hl_561[k]
                   + f_0 * kl_1146[k];

        t_832[k] = -3.0 * hl_562[k]
                   + f_0 * kl_1147[k];

        t_833[k] = -3.0 * hl_563[k]
                   + f_0 * kl_1148[k];

        t_834[k] = -3.0 * hl_564[k]
                   + f_0 * kl_1149[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_2_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hl_565 = buffer.data(hl + 565);
    const auto *hl_566 = buffer.data(hl + 566);
    const auto *hl_567 = buffer.data(hl + 567);
    const auto *hl_568 = buffer.data(hl + 568);
    const auto *hl_569 = buffer.data(hl + 569);
    const auto *hl_570 = buffer.data(hl + 570);
    const auto *hl_571 = buffer.data(hl + 571);
    const auto *hl_572 = buffer.data(hl + 572);
    const auto *hl_573 = buffer.data(hl + 573);
    const auto *hl_574 = buffer.data(hl + 574);
    const auto *hl_575 = buffer.data(hl + 575);
    const auto *hl_576 = buffer.data(hl + 576);
    const auto *hl_577 = buffer.data(hl + 577);
    const auto *hl_578 = buffer.data(hl + 578);
    const auto *hl_579 = buffer.data(hl + 579);
    const auto *hl_580 = buffer.data(hl + 580);
    const auto *hl_581 = buffer.data(hl + 581);
    const auto *hl_582 = buffer.data(hl + 582);
    const auto *hl_583 = buffer.data(hl + 583);
    const auto *hl_584 = buffer.data(hl + 584);
    const auto *hl_585 = buffer.data(hl + 585);
    const auto *hl_586 = buffer.data(hl + 586);
    const auto *hl_587 = buffer.data(hl + 587);
    const auto *hl_588 = buffer.data(hl + 588);
    const auto *hl_589 = buffer.data(hl + 589);
    const auto *hl_590 = buffer.data(hl + 590);
    const auto *hl_591 = buffer.data(hl + 591);
    const auto *hl_592 = buffer.data(hl + 592);
    const auto *hl_593 = buffer.data(hl + 593);
    const auto *hl_594 = buffer.data(hl + 594);
    const auto *hl_595 = buffer.data(hl + 595);
    const auto *hl_596 = buffer.data(hl + 596);
    const auto *hl_597 = buffer.data(hl + 597);
    const auto *hl_598 = buffer.data(hl + 598);
    const auto *hl_599 = buffer.data(hl + 599);
    const auto *hl_600 = buffer.data(hl + 600);
    const auto *hl_601 = buffer.data(hl + 601);
    const auto *hl_602 = buffer.data(hl + 602);
    const auto *hl_603 = buffer.data(hl + 603);
    const auto *hl_604 = buffer.data(hl + 604);
    const auto *hl_605 = buffer.data(hl + 605);
    const auto *hl_606 = buffer.data(hl + 606);
    const auto *hl_607 = buffer.data(hl + 607);
    const auto *hl_608 = buffer.data(hl + 608);
    const auto *hl_609 = buffer.data(hl + 609);
    const auto *hl_610 = buffer.data(hl + 610);
    const auto *hl_611 = buffer.data(hl + 611);
    const auto *hl_612 = buffer.data(hl + 612);
    const auto *hl_613 = buffer.data(hl + 613);
    const auto *hl_614 = buffer.data(hl + 614);
    const auto *hl_615 = buffer.data(hl + 615);
    const auto *hl_616 = buffer.data(hl + 616);
    const auto *hl_617 = buffer.data(hl + 617);
    const auto *hl_618 = buffer.data(hl + 618);
    const auto *hl_619 = buffer.data(hl + 619);
    const auto *hl_620 = buffer.data(hl + 620);
    const auto *hl_621 = buffer.data(hl + 621);
    const auto *hl_622 = buffer.data(hl + 622);
    const auto *hl_623 = buffer.data(hl + 623);
    const auto *hl_624 = buffer.data(hl + 624);
    const auto *hl_625 = buffer.data(hl + 625);
    const auto *hl_626 = buffer.data(hl + 626);
    const auto *hl_627 = buffer.data(hl + 627);
    const auto *hl_628 = buffer.data(hl + 628);
    const auto *hl_629 = buffer.data(hl + 629);
    const auto *hl_630 = buffer.data(hl + 630);
    const auto *hl_631 = buffer.data(hl + 631);
    const auto *hl_632 = buffer.data(hl + 632);
    const auto *hl_633 = buffer.data(hl + 633);
    const auto *hl_634 = buffer.data(hl + 634);
    const auto *hl_635 = buffer.data(hl + 635);
    const auto *hl_636 = buffer.data(hl + 636);
    const auto *hl_637 = buffer.data(hl + 637);
    const auto *hl_638 = buffer.data(hl + 638);
    const auto *hl_639 = buffer.data(hl + 639);
    const auto *hl_640 = buffer.data(hl + 640);
    const auto *hl_641 = buffer.data(hl + 641);
    const auto *hl_642 = buffer.data(hl + 642);
    const auto *hl_643 = buffer.data(hl + 643);
    const auto *hl_644 = buffer.data(hl + 644);
    const auto *hl_645 = buffer.data(hl + 645);
    const auto *hl_646 = buffer.data(hl + 646);
    const auto *hl_647 = buffer.data(hl + 647);
    const auto *hl_648 = buffer.data(hl + 648);
    const auto *hl_649 = buffer.data(hl + 649);
    const auto *hl_650 = buffer.data(hl + 650);
    const auto *hl_651 = buffer.data(hl + 651);
    const auto *hl_652 = buffer.data(hl + 652);
    const auto *hl_653 = buffer.data(hl + 653);
    const auto *hl_654 = buffer.data(hl + 654);
    const auto *hl_655 = buffer.data(hl + 655);
    const auto *hl_656 = buffer.data(hl + 656);
    const auto *hl_657 = buffer.data(hl + 657);
    const auto *hl_658 = buffer.data(hl + 658);
    const auto *hl_659 = buffer.data(hl + 659);
    const auto *hl_660 = buffer.data(hl + 660);
    const auto *hl_661 = buffer.data(hl + 661);
    const auto *hl_662 = buffer.data(hl + 662);
    const auto *hl_663 = buffer.data(hl + 663);
    const auto *hl_664 = buffer.data(hl + 664);
    const auto *hl_665 = buffer.data(hl + 665);
    const auto *hl_666 = buffer.data(hl + 666);
    const auto *hl_667 = buffer.data(hl + 667);
    const auto *hl_668 = buffer.data(hl + 668);
    const auto *hl_669 = buffer.data(hl + 669);
    const auto *hl_670 = buffer.data(hl + 670);
    const auto *hl_671 = buffer.data(hl + 671);
    const auto *hl_672 = buffer.data(hl + 672);
    const auto *hl_673 = buffer.data(hl + 673);
    const auto *hl_674 = buffer.data(hl + 674);
    const auto *hl_675 = buffer.data(hl + 675);
    const auto *hl_676 = buffer.data(hl + 676);
    const auto *hl_677 = buffer.data(hl + 677);
    const auto *hl_678 = buffer.data(hl + 678);
    const auto *hl_679 = buffer.data(hl + 679);
    const auto *hl_680 = buffer.data(hl + 680);
    const auto *hl_681 = buffer.data(hl + 681);

    const auto *kl_1150 = buffer.data(kl + 1150);
    const auto *kl_1151 = buffer.data(kl + 1151);
    const auto *kl_1152 = buffer.data(kl + 1152);
    const auto *kl_1153 = buffer.data(kl + 1153);
    const auto *kl_1154 = buffer.data(kl + 1154);
    const auto *kl_1155 = buffer.data(kl + 1155);
    const auto *kl_1156 = buffer.data(kl + 1156);
    const auto *kl_1157 = buffer.data(kl + 1157);
    const auto *kl_1158 = buffer.data(kl + 1158);
    const auto *kl_1159 = buffer.data(kl + 1159);
    const auto *kl_1160 = buffer.data(kl + 1160);
    const auto *kl_1161 = buffer.data(kl + 1161);
    const auto *kl_1162 = buffer.data(kl + 1162);
    const auto *kl_1163 = buffer.data(kl + 1163);
    const auto *kl_1164 = buffer.data(kl + 1164);
    const auto *kl_1165 = buffer.data(kl + 1165);
    const auto *kl_1166 = buffer.data(kl + 1166);
    const auto *kl_1167 = buffer.data(kl + 1167);
    const auto *kl_1168 = buffer.data(kl + 1168);
    const auto *kl_1169 = buffer.data(kl + 1169);
    const auto *kl_1170 = buffer.data(kl + 1170);
    const auto *kl_1171 = buffer.data(kl + 1171);
    const auto *kl_1172 = buffer.data(kl + 1172);
    const auto *kl_1173 = buffer.data(kl + 1173);
    const auto *kl_1174 = buffer.data(kl + 1174);
    const auto *kl_1175 = buffer.data(kl + 1175);
    const auto *kl_1176 = buffer.data(kl + 1176);
    const auto *kl_1177 = buffer.data(kl + 1177);
    const auto *kl_1178 = buffer.data(kl + 1178);
    const auto *kl_1179 = buffer.data(kl + 1179);
    const auto *kl_1180 = buffer.data(kl + 1180);
    const auto *kl_1181 = buffer.data(kl + 1181);
    const auto *kl_1182 = buffer.data(kl + 1182);
    const auto *kl_1183 = buffer.data(kl + 1183);
    const auto *kl_1184 = buffer.data(kl + 1184);
    const auto *kl_1185 = buffer.data(kl + 1185);
    const auto *kl_1186 = buffer.data(kl + 1186);
    const auto *kl_1187 = buffer.data(kl + 1187);
    const auto *kl_1188 = buffer.data(kl + 1188);
    const auto *kl_1189 = buffer.data(kl + 1189);
    const auto *kl_1190 = buffer.data(kl + 1190);
    const auto *kl_1191 = buffer.data(kl + 1191);
    const auto *kl_1192 = buffer.data(kl + 1192);
    const auto *kl_1193 = buffer.data(kl + 1193);
    const auto *kl_1194 = buffer.data(kl + 1194);
    const auto *kl_1195 = buffer.data(kl + 1195);
    const auto *kl_1196 = buffer.data(kl + 1196);
    const auto *kl_1197 = buffer.data(kl + 1197);
    const auto *kl_1198 = buffer.data(kl + 1198);
    const auto *kl_1199 = buffer.data(kl + 1199);
    const auto *kl_1200 = buffer.data(kl + 1200);
    const auto *kl_1201 = buffer.data(kl + 1201);
    const auto *kl_1202 = buffer.data(kl + 1202);
    const auto *kl_1203 = buffer.data(kl + 1203);
    const auto *kl_1204 = buffer.data(kl + 1204);
    const auto *kl_1205 = buffer.data(kl + 1205);
    const auto *kl_1206 = buffer.data(kl + 1206);
    const auto *kl_1207 = buffer.data(kl + 1207);
    const auto *kl_1208 = buffer.data(kl + 1208);
    const auto *kl_1209 = buffer.data(kl + 1209);
    const auto *kl_1210 = buffer.data(kl + 1210);
    const auto *kl_1211 = buffer.data(kl + 1211);
    const auto *kl_1212 = buffer.data(kl + 1212);
    const auto *kl_1213 = buffer.data(kl + 1213);
    const auto *kl_1214 = buffer.data(kl + 1214);
    const auto *kl_1215 = buffer.data(kl + 1215);
    const auto *kl_1216 = buffer.data(kl + 1216);
    const auto *kl_1217 = buffer.data(kl + 1217);
    const auto *kl_1218 = buffer.data(kl + 1218);
    const auto *kl_1219 = buffer.data(kl + 1219);
    const auto *kl_1220 = buffer.data(kl + 1220);
    const auto *kl_1221 = buffer.data(kl + 1221);
    const auto *kl_1222 = buffer.data(kl + 1222);
    const auto *kl_1223 = buffer.data(kl + 1223);
    const auto *kl_1224 = buffer.data(kl + 1224);
    const auto *kl_1225 = buffer.data(kl + 1225);
    const auto *kl_1226 = buffer.data(kl + 1226);
    const auto *kl_1227 = buffer.data(kl + 1227);
    const auto *kl_1228 = buffer.data(kl + 1228);
    const auto *kl_1229 = buffer.data(kl + 1229);
    const auto *kl_1230 = buffer.data(kl + 1230);
    const auto *kl_1231 = buffer.data(kl + 1231);
    const auto *kl_1232 = buffer.data(kl + 1232);
    const auto *kl_1233 = buffer.data(kl + 1233);
    const auto *kl_1234 = buffer.data(kl + 1234);
    const auto *kl_1235 = buffer.data(kl + 1235);
    const auto *kl_1236 = buffer.data(kl + 1236);
    const auto *kl_1237 = buffer.data(kl + 1237);
    const auto *kl_1238 = buffer.data(kl + 1238);
    const auto *kl_1239 = buffer.data(kl + 1239);
    const auto *kl_1240 = buffer.data(kl + 1240);
    const auto *kl_1241 = buffer.data(kl + 1241);
    const auto *kl_1242 = buffer.data(kl + 1242);
    const auto *kl_1243 = buffer.data(kl + 1243);
    const auto *kl_1244 = buffer.data(kl + 1244);
    const auto *kl_1245 = buffer.data(kl + 1245);
    const auto *kl_1246 = buffer.data(kl + 1246);
    const auto *kl_1247 = buffer.data(kl + 1247);
    const auto *kl_1248 = buffer.data(kl + 1248);
    const auto *kl_1249 = buffer.data(kl + 1249);
    const auto *kl_1250 = buffer.data(kl + 1250);
    const auto *kl_1251 = buffer.data(kl + 1251);
    const auto *kl_1252 = buffer.data(kl + 1252);
    const auto *kl_1253 = buffer.data(kl + 1253);
    const auto *kl_1254 = buffer.data(kl + 1254);
    const auto *kl_1255 = buffer.data(kl + 1255);
    const auto *kl_1256 = buffer.data(kl + 1256);
    const auto *kl_1257 = buffer.data(kl + 1257);
    const auto *kl_1258 = buffer.data(kl + 1258);
    const auto *kl_1259 = buffer.data(kl + 1259);
    const auto *kl_1305 = buffer.data(kl + 1305);
    const auto *kl_1306 = buffer.data(kl + 1306);
    const auto *kl_1307 = buffer.data(kl + 1307);
    const auto *kl_1308 = buffer.data(kl + 1308);
    const auto *kl_1309 = buffer.data(kl + 1309);
    const auto *kl_1310 = buffer.data(kl + 1310);
    const auto *kl_1311 = buffer.data(kl + 1311);
    const auto *kl_1312 = buffer.data(kl + 1312);
    const auto *kl_1313 = buffer.data(kl + 1313);
    const auto *kl_1314 = buffer.data(kl + 1314);
    const auto *kl_1315 = buffer.data(kl + 1315);
    const auto *kl_1316 = buffer.data(kl + 1316);
    const auto *kl_1317 = buffer.data(kl + 1317);
    const auto *kl_1318 = buffer.data(kl + 1318);
    const auto *kl_1319 = buffer.data(kl + 1319);
    const auto *kl_1320 = buffer.data(kl + 1320);
    const auto *kl_1321 = buffer.data(kl + 1321);
    const auto *kl_1322 = buffer.data(kl + 1322);
    const auto *kl_1323 = buffer.data(kl + 1323);
    const auto *kl_1324 = buffer.data(kl + 1324);
    const auto *kl_1325 = buffer.data(kl + 1325);
    const auto *kl_1326 = buffer.data(kl + 1326);
    const auto *kl_1327 = buffer.data(kl + 1327);
    const auto *kl_1328 = buffer.data(kl + 1328);
    const auto *kl_1329 = buffer.data(kl + 1329);
    const auto *kl_1330 = buffer.data(kl + 1330);
    const auto *kl_1331 = buffer.data(kl + 1331);
    const auto *kl_1332 = buffer.data(kl + 1332);
    const auto *kl_1333 = buffer.data(kl + 1333);
    const auto *kl_1334 = buffer.data(kl + 1334);
    const auto *kl_1335 = buffer.data(kl + 1335);
    const auto *kl_1336 = buffer.data(kl + 1336);
    const auto *kl_1337 = buffer.data(kl + 1337);
    const auto *kl_1338 = buffer.data(kl + 1338);
    const auto *kl_1339 = buffer.data(kl + 1339);
    const auto *kl_1340 = buffer.data(kl + 1340);
    const auto *kl_1341 = buffer.data(kl + 1341);
    const auto *kl_1342 = buffer.data(kl + 1342);
    const auto *kl_1343 = buffer.data(kl + 1343);
    const auto *kl_1344 = buffer.data(kl + 1344);
    const auto *kl_1345 = buffer.data(kl + 1345);
    const auto *kl_1346 = buffer.data(kl + 1346);
    const auto *kl_1347 = buffer.data(kl + 1347);
    const auto *kl_1348 = buffer.data(kl + 1348);
    const auto *kl_1349 = buffer.data(kl + 1349);
    const auto *kl_1350 = buffer.data(kl + 1350);
    const auto *kl_1351 = buffer.data(kl + 1351);
    const auto *kl_1352 = buffer.data(kl + 1352);
    const auto *kl_1353 = buffer.data(kl + 1353);
    const auto *kl_1354 = buffer.data(kl + 1354);
    const auto *kl_1355 = buffer.data(kl + 1355);
    const auto *kl_1356 = buffer.data(kl + 1356);

#pragma omp simd aligned(t_835, t_836, t_837, t_838, t_839, hl_565, hl_566, hl_567, hl_568, \
                         hl_569, kl_1150, kl_1151, kl_1152, kl_1153, \
                         kl_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = -3.0 * hl_565[k]
                   + f_0 * kl_1150[k];

        t_836[k] = -3.0 * hl_566[k]
                   + f_0 * kl_1151[k];

        t_837[k] = -3.0 * hl_567[k]
                   + f_0 * kl_1152[k];

        t_838[k] = -3.0 * hl_568[k]
                   + f_0 * kl_1153[k];

        t_839[k] = -3.0 * hl_569[k]
                   + f_0 * kl_1154[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, hl_570, hl_571, hl_572, hl_573, \
                         hl_574, kl_1155, kl_1156, kl_1157, kl_1158, \
                         kl_1159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = -3.0 * hl_570[k]
                   + f_0 * kl_1155[k];

        t_841[k] = -3.0 * hl_571[k]
                   + f_0 * kl_1156[k];

        t_842[k] = -3.0 * hl_572[k]
                   + f_0 * kl_1157[k];

        t_843[k] = -3.0 * hl_573[k]
                   + f_0 * kl_1158[k];

        t_844[k] = -3.0 * hl_574[k]
                   + f_0 * kl_1159[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, t_849, hl_575, hl_576, hl_577, hl_578, \
                         hl_579, kl_1160, kl_1161, kl_1162, kl_1163, \
                         kl_1164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = -3.0 * hl_575[k]
                   + f_0 * kl_1160[k];

        t_846[k] = -3.0 * hl_576[k]
                   + f_0 * kl_1161[k];

        t_847[k] = -3.0 * hl_577[k]
                   + f_0 * kl_1162[k];

        t_848[k] = -3.0 * hl_578[k]
                   + f_0 * kl_1163[k];

        t_849[k] = -3.0 * hl_579[k]
                   + f_0 * kl_1164[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, hl_580, hl_581, hl_582, hl_583, \
                         hl_584, kl_1165, kl_1166, kl_1167, kl_1168, \
                         kl_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = -3.0 * hl_580[k]
                   + f_0 * kl_1165[k];

        t_851[k] = -3.0 * hl_581[k]
                   + f_0 * kl_1166[k];

        t_852[k] = -3.0 * hl_582[k]
                   + f_0 * kl_1167[k];

        t_853[k] = -3.0 * hl_583[k]
                   + f_0 * kl_1168[k];

        t_854[k] = -3.0 * hl_584[k]
                   + f_0 * kl_1169[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, hl_585, hl_586, hl_587, hl_588, \
                         hl_589, kl_1170, kl_1171, kl_1172, kl_1173, \
                         kl_1174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = -4.0 * hl_585[k]
                   + f_0 * kl_1170[k];

        t_856[k] = -4.0 * hl_586[k]
                   + f_0 * kl_1171[k];

        t_857[k] = -4.0 * hl_587[k]
                   + f_0 * kl_1172[k];

        t_858[k] = -4.0 * hl_588[k]
                   + f_0 * kl_1173[k];

        t_859[k] = -4.0 * hl_589[k]
                   + f_0 * kl_1174[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, t_864, hl_590, hl_591, hl_592, hl_593, \
                         hl_594, kl_1175, kl_1176, kl_1177, kl_1178, \
                         kl_1179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = -4.0 * hl_590[k]
                   + f_0 * kl_1175[k];

        t_861[k] = -4.0 * hl_591[k]
                   + f_0 * kl_1176[k];

        t_862[k] = -4.0 * hl_592[k]
                   + f_0 * kl_1177[k];

        t_863[k] = -4.0 * hl_593[k]
                   + f_0 * kl_1178[k];

        t_864[k] = -4.0 * hl_594[k]
                   + f_0 * kl_1179[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, t_869, hl_595, hl_596, hl_597, hl_598, \
                         hl_599, kl_1180, kl_1181, kl_1182, kl_1183, \
                         kl_1184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = -4.0 * hl_595[k]
                   + f_0 * kl_1180[k];

        t_866[k] = -4.0 * hl_596[k]
                   + f_0 * kl_1181[k];

        t_867[k] = -4.0 * hl_597[k]
                   + f_0 * kl_1182[k];

        t_868[k] = -4.0 * hl_598[k]
                   + f_0 * kl_1183[k];

        t_869[k] = -4.0 * hl_599[k]
                   + f_0 * kl_1184[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, t_874, hl_600, hl_601, hl_602, hl_603, \
                         hl_604, kl_1185, kl_1186, kl_1187, kl_1188, \
                         kl_1189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = -4.0 * hl_600[k]
                   + f_0 * kl_1185[k];

        t_871[k] = -4.0 * hl_601[k]
                   + f_0 * kl_1186[k];

        t_872[k] = -4.0 * hl_602[k]
                   + f_0 * kl_1187[k];

        t_873[k] = -4.0 * hl_603[k]
                   + f_0 * kl_1188[k];

        t_874[k] = -4.0 * hl_604[k]
                   + f_0 * kl_1189[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, t_879, hl_605, hl_606, hl_607, hl_608, \
                         hl_609, kl_1190, kl_1191, kl_1192, kl_1193, \
                         kl_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = -4.0 * hl_605[k]
                   + f_0 * kl_1190[k];

        t_876[k] = -4.0 * hl_606[k]
                   + f_0 * kl_1191[k];

        t_877[k] = -4.0 * hl_607[k]
                   + f_0 * kl_1192[k];

        t_878[k] = -4.0 * hl_608[k]
                   + f_0 * kl_1193[k];

        t_879[k] = -4.0 * hl_609[k]
                   + f_0 * kl_1194[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, t_884, hl_610, hl_611, hl_612, hl_613, \
                         hl_614, kl_1195, kl_1196, kl_1197, kl_1198, \
                         kl_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = -4.0 * hl_610[k]
                   + f_0 * kl_1195[k];

        t_881[k] = -4.0 * hl_611[k]
                   + f_0 * kl_1196[k];

        t_882[k] = -4.0 * hl_612[k]
                   + f_0 * kl_1197[k];

        t_883[k] = -4.0 * hl_613[k]
                   + f_0 * kl_1198[k];

        t_884[k] = -4.0 * hl_614[k]
                   + f_0 * kl_1199[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, hl_615, hl_616, hl_617, hl_618, \
                         hl_619, kl_1200, kl_1201, kl_1202, kl_1203, \
                         kl_1204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = -4.0 * hl_615[k]
                   + f_0 * kl_1200[k];

        t_886[k] = -4.0 * hl_616[k]
                   + f_0 * kl_1201[k];

        t_887[k] = -4.0 * hl_617[k]
                   + f_0 * kl_1202[k];

        t_888[k] = -4.0 * hl_618[k]
                   + f_0 * kl_1203[k];

        t_889[k] = -4.0 * hl_619[k]
                   + f_0 * kl_1204[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, hl_620, hl_621, hl_622, hl_623, \
                         hl_624, kl_1205, kl_1206, kl_1207, kl_1208, \
                         kl_1209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = -4.0 * hl_620[k]
                   + f_0 * kl_1205[k];

        t_891[k] = -4.0 * hl_621[k]
                   + f_0 * kl_1206[k];

        t_892[k] = -4.0 * hl_622[k]
                   + f_0 * kl_1207[k];

        t_893[k] = -4.0 * hl_623[k]
                   + f_0 * kl_1208[k];

        t_894[k] = -4.0 * hl_624[k]
                   + f_0 * kl_1209[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, t_899, hl_625, hl_626, hl_627, hl_628, \
                         hl_629, kl_1210, kl_1211, kl_1212, kl_1213, \
                         kl_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = -4.0 * hl_625[k]
                   + f_0 * kl_1210[k];

        t_896[k] = -4.0 * hl_626[k]
                   + f_0 * kl_1211[k];

        t_897[k] = -4.0 * hl_627[k]
                   + f_0 * kl_1212[k];

        t_898[k] = -4.0 * hl_628[k]
                   + f_0 * kl_1213[k];

        t_899[k] = -4.0 * hl_629[k]
                   + f_0 * kl_1214[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, t_904, hl_630, hl_631, hl_632, hl_633, \
                         hl_634, kl_1215, kl_1216, kl_1217, kl_1218, \
                         kl_1219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = -5.0 * hl_630[k]
                   + f_0 * kl_1215[k];

        t_901[k] = -5.0 * hl_631[k]
                   + f_0 * kl_1216[k];

        t_902[k] = -5.0 * hl_632[k]
                   + f_0 * kl_1217[k];

        t_903[k] = -5.0 * hl_633[k]
                   + f_0 * kl_1218[k];

        t_904[k] = -5.0 * hl_634[k]
                   + f_0 * kl_1219[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, t_909, hl_635, hl_636, hl_637, hl_638, \
                         hl_639, kl_1220, kl_1221, kl_1222, kl_1223, \
                         kl_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = -5.0 * hl_635[k]
                   + f_0 * kl_1220[k];

        t_906[k] = -5.0 * hl_636[k]
                   + f_0 * kl_1221[k];

        t_907[k] = -5.0 * hl_637[k]
                   + f_0 * kl_1222[k];

        t_908[k] = -5.0 * hl_638[k]
                   + f_0 * kl_1223[k];

        t_909[k] = -5.0 * hl_639[k]
                   + f_0 * kl_1224[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, t_914, hl_640, hl_641, hl_642, hl_643, \
                         hl_644, kl_1225, kl_1226, kl_1227, kl_1228, \
                         kl_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = -5.0 * hl_640[k]
                   + f_0 * kl_1225[k];

        t_911[k] = -5.0 * hl_641[k]
                   + f_0 * kl_1226[k];

        t_912[k] = -5.0 * hl_642[k]
                   + f_0 * kl_1227[k];

        t_913[k] = -5.0 * hl_643[k]
                   + f_0 * kl_1228[k];

        t_914[k] = -5.0 * hl_644[k]
                   + f_0 * kl_1229[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, t_919, hl_645, hl_646, hl_647, hl_648, \
                         hl_649, kl_1230, kl_1231, kl_1232, kl_1233, \
                         kl_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = -5.0 * hl_645[k]
                   + f_0 * kl_1230[k];

        t_916[k] = -5.0 * hl_646[k]
                   + f_0 * kl_1231[k];

        t_917[k] = -5.0 * hl_647[k]
                   + f_0 * kl_1232[k];

        t_918[k] = -5.0 * hl_648[k]
                   + f_0 * kl_1233[k];

        t_919[k] = -5.0 * hl_649[k]
                   + f_0 * kl_1234[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, t_924, hl_650, hl_651, hl_652, hl_653, \
                         hl_654, kl_1235, kl_1236, kl_1237, kl_1238, \
                         kl_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = -5.0 * hl_650[k]
                   + f_0 * kl_1235[k];

        t_921[k] = -5.0 * hl_651[k]
                   + f_0 * kl_1236[k];

        t_922[k] = -5.0 * hl_652[k]
                   + f_0 * kl_1237[k];

        t_923[k] = -5.0 * hl_653[k]
                   + f_0 * kl_1238[k];

        t_924[k] = -5.0 * hl_654[k]
                   + f_0 * kl_1239[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, t_929, hl_655, hl_656, hl_657, hl_658, \
                         hl_659, kl_1240, kl_1241, kl_1242, kl_1243, \
                         kl_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = -5.0 * hl_655[k]
                   + f_0 * kl_1240[k];

        t_926[k] = -5.0 * hl_656[k]
                   + f_0 * kl_1241[k];

        t_927[k] = -5.0 * hl_657[k]
                   + f_0 * kl_1242[k];

        t_928[k] = -5.0 * hl_658[k]
                   + f_0 * kl_1243[k];

        t_929[k] = -5.0 * hl_659[k]
                   + f_0 * kl_1244[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, hl_660, hl_661, hl_662, hl_663, \
                         hl_664, kl_1245, kl_1246, kl_1247, kl_1248, \
                         kl_1249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = -5.0 * hl_660[k]
                   + f_0 * kl_1245[k];

        t_931[k] = -5.0 * hl_661[k]
                   + f_0 * kl_1246[k];

        t_932[k] = -5.0 * hl_662[k]
                   + f_0 * kl_1247[k];

        t_933[k] = -5.0 * hl_663[k]
                   + f_0 * kl_1248[k];

        t_934[k] = -5.0 * hl_664[k]
                   + f_0 * kl_1249[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, t_939, hl_665, hl_666, hl_667, hl_668, \
                         hl_669, kl_1250, kl_1251, kl_1252, kl_1253, \
                         kl_1254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = -5.0 * hl_665[k]
                   + f_0 * kl_1250[k];

        t_936[k] = -5.0 * hl_666[k]
                   + f_0 * kl_1251[k];

        t_937[k] = -5.0 * hl_667[k]
                   + f_0 * kl_1252[k];

        t_938[k] = -5.0 * hl_668[k]
                   + f_0 * kl_1253[k];

        t_939[k] = -5.0 * hl_669[k]
                   + f_0 * kl_1254[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, t_944, hl_670, hl_671, hl_672, hl_673, \
                         hl_674, kl_1255, kl_1256, kl_1257, kl_1258, \
                         kl_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = -5.0 * hl_670[k]
                   + f_0 * kl_1255[k];

        t_941[k] = -5.0 * hl_671[k]
                   + f_0 * kl_1256[k];

        t_942[k] = -5.0 * hl_672[k]
                   + f_0 * kl_1257[k];

        t_943[k] = -5.0 * hl_673[k]
                   + f_0 * kl_1258[k];

        t_944[k] = -5.0 * hl_674[k]
                   + f_0 * kl_1259[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, t_949, t_950, t_951, t_952, kl_1305, \
                         kl_1306, kl_1307, kl_1308, kl_1309, kl_1310, kl_1311, \
                         kl_1312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = f_0 * kl_1305[k];

        t_946[k] = f_0 * kl_1306[k];

        t_947[k] = f_0 * kl_1307[k];

        t_948[k] = f_0 * kl_1308[k];

        t_949[k] = f_0 * kl_1309[k];

        t_950[k] = f_0 * kl_1310[k];

        t_951[k] = f_0 * kl_1311[k];

        t_952[k] = f_0 * kl_1312[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, t_957, t_958, t_959, t_960, kl_1313, \
                         kl_1314, kl_1315, kl_1316, kl_1317, kl_1318, kl_1319, \
                         kl_1320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_0 * kl_1313[k];

        t_954[k] = f_0 * kl_1314[k];

        t_955[k] = f_0 * kl_1315[k];

        t_956[k] = f_0 * kl_1316[k];

        t_957[k] = f_0 * kl_1317[k];

        t_958[k] = f_0 * kl_1318[k];

        t_959[k] = f_0 * kl_1319[k];

        t_960[k] = f_0 * kl_1320[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, t_964, t_965, t_966, t_967, t_968, kl_1321, \
                         kl_1322, kl_1323, kl_1324, kl_1325, kl_1326, kl_1327, \
                         kl_1328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_0 * kl_1321[k];

        t_962[k] = f_0 * kl_1322[k];

        t_963[k] = f_0 * kl_1323[k];

        t_964[k] = f_0 * kl_1324[k];

        t_965[k] = f_0 * kl_1325[k];

        t_966[k] = f_0 * kl_1326[k];

        t_967[k] = f_0 * kl_1327[k];

        t_968[k] = f_0 * kl_1328[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, t_972, t_973, t_974, t_975, t_976, kl_1329, \
                         kl_1330, kl_1331, kl_1332, kl_1333, kl_1334, kl_1335, \
                         kl_1336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_0 * kl_1329[k];

        t_970[k] = f_0 * kl_1330[k];

        t_971[k] = f_0 * kl_1331[k];

        t_972[k] = f_0 * kl_1332[k];

        t_973[k] = f_0 * kl_1333[k];

        t_974[k] = f_0 * kl_1334[k];

        t_975[k] = f_0 * kl_1335[k];

        t_976[k] = f_0 * kl_1336[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, t_980, t_981, t_982, t_983, t_984, kl_1337, \
                         kl_1338, kl_1339, kl_1340, kl_1341, kl_1342, kl_1343, \
                         kl_1344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_0 * kl_1337[k];

        t_978[k] = f_0 * kl_1338[k];

        t_979[k] = f_0 * kl_1339[k];

        t_980[k] = f_0 * kl_1340[k];

        t_981[k] = f_0 * kl_1341[k];

        t_982[k] = f_0 * kl_1342[k];

        t_983[k] = f_0 * kl_1343[k];

        t_984[k] = f_0 * kl_1344[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, t_989, t_990, t_991, hl_675, hl_676, \
                         kl_1345, kl_1346, kl_1347, kl_1348, kl_1349, kl_1350, \
                         kl_1351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_0 * kl_1345[k];

        t_986[k] = f_0 * kl_1346[k];

        t_987[k] = f_0 * kl_1347[k];

        t_988[k] = f_0 * kl_1348[k];

        t_989[k] = f_0 * kl_1349[k];

        t_990[k] = -hl_675[k]
                   + f_0 * kl_1350[k];

        t_991[k] = -hl_676[k]
                   + f_0 * kl_1351[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, t_996, hl_677, hl_678, hl_679, hl_680, \
                         hl_681, kl_1352, kl_1353, kl_1354, kl_1355, \
                         kl_1356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = -hl_677[k]
                   + f_0 * kl_1352[k];

        t_993[k] = -hl_678[k]
                   + f_0 * kl_1353[k];

        t_994[k] = -hl_679[k]
                   + f_0 * kl_1354[k];

        t_995[k] = -hl_680[k]
                   + f_0 * kl_1355[k];

        t_996[k] = -hl_681[k]
                   + f_0 * kl_1356[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_2_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_1146 = buffer.data(target + 1146);

    const auto *hl_682 = buffer.data(hl + 682);
    const auto *hl_683 = buffer.data(hl + 683);
    const auto *hl_684 = buffer.data(hl + 684);
    const auto *hl_685 = buffer.data(hl + 685);
    const auto *hl_686 = buffer.data(hl + 686);
    const auto *hl_687 = buffer.data(hl + 687);
    const auto *hl_688 = buffer.data(hl + 688);
    const auto *hl_689 = buffer.data(hl + 689);
    const auto *hl_690 = buffer.data(hl + 690);
    const auto *hl_691 = buffer.data(hl + 691);
    const auto *hl_692 = buffer.data(hl + 692);
    const auto *hl_693 = buffer.data(hl + 693);
    const auto *hl_694 = buffer.data(hl + 694);
    const auto *hl_695 = buffer.data(hl + 695);
    const auto *hl_696 = buffer.data(hl + 696);
    const auto *hl_697 = buffer.data(hl + 697);
    const auto *hl_698 = buffer.data(hl + 698);
    const auto *hl_699 = buffer.data(hl + 699);
    const auto *hl_700 = buffer.data(hl + 700);
    const auto *hl_701 = buffer.data(hl + 701);
    const auto *hl_702 = buffer.data(hl + 702);
    const auto *hl_703 = buffer.data(hl + 703);
    const auto *hl_704 = buffer.data(hl + 704);
    const auto *hl_705 = buffer.data(hl + 705);
    const auto *hl_706 = buffer.data(hl + 706);
    const auto *hl_707 = buffer.data(hl + 707);
    const auto *hl_708 = buffer.data(hl + 708);
    const auto *hl_709 = buffer.data(hl + 709);
    const auto *hl_710 = buffer.data(hl + 710);
    const auto *hl_711 = buffer.data(hl + 711);
    const auto *hl_712 = buffer.data(hl + 712);
    const auto *hl_713 = buffer.data(hl + 713);
    const auto *hl_714 = buffer.data(hl + 714);
    const auto *hl_715 = buffer.data(hl + 715);
    const auto *hl_716 = buffer.data(hl + 716);
    const auto *hl_717 = buffer.data(hl + 717);
    const auto *hl_718 = buffer.data(hl + 718);
    const auto *hl_719 = buffer.data(hl + 719);
    const auto *hl_720 = buffer.data(hl + 720);
    const auto *hl_721 = buffer.data(hl + 721);
    const auto *hl_722 = buffer.data(hl + 722);
    const auto *hl_723 = buffer.data(hl + 723);
    const auto *hl_724 = buffer.data(hl + 724);
    const auto *hl_725 = buffer.data(hl + 725);
    const auto *hl_726 = buffer.data(hl + 726);
    const auto *hl_727 = buffer.data(hl + 727);
    const auto *hl_728 = buffer.data(hl + 728);
    const auto *hl_729 = buffer.data(hl + 729);
    const auto *hl_730 = buffer.data(hl + 730);
    const auto *hl_731 = buffer.data(hl + 731);
    const auto *hl_732 = buffer.data(hl + 732);
    const auto *hl_733 = buffer.data(hl + 733);
    const auto *hl_734 = buffer.data(hl + 734);
    const auto *hl_735 = buffer.data(hl + 735);
    const auto *hl_736 = buffer.data(hl + 736);
    const auto *hl_737 = buffer.data(hl + 737);
    const auto *hl_738 = buffer.data(hl + 738);
    const auto *hl_739 = buffer.data(hl + 739);
    const auto *hl_740 = buffer.data(hl + 740);
    const auto *hl_741 = buffer.data(hl + 741);
    const auto *hl_742 = buffer.data(hl + 742);
    const auto *hl_743 = buffer.data(hl + 743);
    const auto *hl_744 = buffer.data(hl + 744);
    const auto *hl_745 = buffer.data(hl + 745);
    const auto *hl_746 = buffer.data(hl + 746);
    const auto *hl_747 = buffer.data(hl + 747);
    const auto *hl_748 = buffer.data(hl + 748);
    const auto *hl_749 = buffer.data(hl + 749);
    const auto *hl_750 = buffer.data(hl + 750);
    const auto *hl_751 = buffer.data(hl + 751);
    const auto *hl_752 = buffer.data(hl + 752);
    const auto *hl_753 = buffer.data(hl + 753);
    const auto *hl_754 = buffer.data(hl + 754);
    const auto *hl_755 = buffer.data(hl + 755);
    const auto *hl_756 = buffer.data(hl + 756);
    const auto *hl_757 = buffer.data(hl + 757);
    const auto *hl_758 = buffer.data(hl + 758);
    const auto *hl_759 = buffer.data(hl + 759);
    const auto *hl_760 = buffer.data(hl + 760);
    const auto *hl_761 = buffer.data(hl + 761);
    const auto *hl_762 = buffer.data(hl + 762);
    const auto *hl_763 = buffer.data(hl + 763);
    const auto *hl_764 = buffer.data(hl + 764);
    const auto *hl_765 = buffer.data(hl + 765);
    const auto *hl_766 = buffer.data(hl + 766);
    const auto *hl_767 = buffer.data(hl + 767);
    const auto *hl_768 = buffer.data(hl + 768);
    const auto *hl_769 = buffer.data(hl + 769);
    const auto *hl_770 = buffer.data(hl + 770);
    const auto *hl_771 = buffer.data(hl + 771);
    const auto *hl_772 = buffer.data(hl + 772);
    const auto *hl_773 = buffer.data(hl + 773);
    const auto *hl_774 = buffer.data(hl + 774);
    const auto *hl_775 = buffer.data(hl + 775);
    const auto *hl_776 = buffer.data(hl + 776);
    const auto *hl_777 = buffer.data(hl + 777);
    const auto *hl_778 = buffer.data(hl + 778);
    const auto *hl_779 = buffer.data(hl + 779);
    const auto *hl_780 = buffer.data(hl + 780);
    const auto *hl_781 = buffer.data(hl + 781);
    const auto *hl_782 = buffer.data(hl + 782);
    const auto *hl_783 = buffer.data(hl + 783);
    const auto *hl_784 = buffer.data(hl + 784);
    const auto *hl_785 = buffer.data(hl + 785);
    const auto *hl_786 = buffer.data(hl + 786);
    const auto *hl_787 = buffer.data(hl + 787);
    const auto *hl_788 = buffer.data(hl + 788);
    const auto *hl_789 = buffer.data(hl + 789);
    const auto *hl_790 = buffer.data(hl + 790);
    const auto *hl_791 = buffer.data(hl + 791);
    const auto *hl_792 = buffer.data(hl + 792);
    const auto *hl_793 = buffer.data(hl + 793);
    const auto *hl_794 = buffer.data(hl + 794);
    const auto *hl_795 = buffer.data(hl + 795);
    const auto *hl_796 = buffer.data(hl + 796);
    const auto *hl_797 = buffer.data(hl + 797);
    const auto *hl_798 = buffer.data(hl + 798);
    const auto *hl_799 = buffer.data(hl + 799);
    const auto *hl_800 = buffer.data(hl + 800);
    const auto *hl_801 = buffer.data(hl + 801);
    const auto *hl_802 = buffer.data(hl + 802);
    const auto *hl_803 = buffer.data(hl + 803);
    const auto *hl_804 = buffer.data(hl + 804);
    const auto *hl_805 = buffer.data(hl + 805);
    const auto *hl_806 = buffer.data(hl + 806);
    const auto *hl_807 = buffer.data(hl + 807);
    const auto *hl_808 = buffer.data(hl + 808);
    const auto *hl_809 = buffer.data(hl + 809);
    const auto *hl_810 = buffer.data(hl + 810);
    const auto *hl_811 = buffer.data(hl + 811);
    const auto *hl_812 = buffer.data(hl + 812);
    const auto *hl_813 = buffer.data(hl + 813);
    const auto *hl_814 = buffer.data(hl + 814);
    const auto *hl_815 = buffer.data(hl + 815);
    const auto *hl_816 = buffer.data(hl + 816);
    const auto *hl_817 = buffer.data(hl + 817);
    const auto *hl_818 = buffer.data(hl + 818);
    const auto *hl_819 = buffer.data(hl + 819);
    const auto *hl_820 = buffer.data(hl + 820);
    const auto *hl_821 = buffer.data(hl + 821);
    const auto *hl_822 = buffer.data(hl + 822);
    const auto *hl_823 = buffer.data(hl + 823);
    const auto *hl_824 = buffer.data(hl + 824);
    const auto *hl_825 = buffer.data(hl + 825);
    const auto *hl_826 = buffer.data(hl + 826);
    const auto *hl_827 = buffer.data(hl + 827);
    const auto *hl_828 = buffer.data(hl + 828);
    const auto *hl_829 = buffer.data(hl + 829);
    const auto *hl_830 = buffer.data(hl + 830);
    const auto *hl_831 = buffer.data(hl + 831);

    const auto *kl_1357 = buffer.data(kl + 1357);
    const auto *kl_1358 = buffer.data(kl + 1358);
    const auto *kl_1359 = buffer.data(kl + 1359);
    const auto *kl_1360 = buffer.data(kl + 1360);
    const auto *kl_1361 = buffer.data(kl + 1361);
    const auto *kl_1362 = buffer.data(kl + 1362);
    const auto *kl_1363 = buffer.data(kl + 1363);
    const auto *kl_1364 = buffer.data(kl + 1364);
    const auto *kl_1365 = buffer.data(kl + 1365);
    const auto *kl_1366 = buffer.data(kl + 1366);
    const auto *kl_1367 = buffer.data(kl + 1367);
    const auto *kl_1368 = buffer.data(kl + 1368);
    const auto *kl_1369 = buffer.data(kl + 1369);
    const auto *kl_1370 = buffer.data(kl + 1370);
    const auto *kl_1371 = buffer.data(kl + 1371);
    const auto *kl_1372 = buffer.data(kl + 1372);
    const auto *kl_1373 = buffer.data(kl + 1373);
    const auto *kl_1374 = buffer.data(kl + 1374);
    const auto *kl_1375 = buffer.data(kl + 1375);
    const auto *kl_1376 = buffer.data(kl + 1376);
    const auto *kl_1377 = buffer.data(kl + 1377);
    const auto *kl_1378 = buffer.data(kl + 1378);
    const auto *kl_1379 = buffer.data(kl + 1379);
    const auto *kl_1380 = buffer.data(kl + 1380);
    const auto *kl_1381 = buffer.data(kl + 1381);
    const auto *kl_1382 = buffer.data(kl + 1382);
    const auto *kl_1383 = buffer.data(kl + 1383);
    const auto *kl_1384 = buffer.data(kl + 1384);
    const auto *kl_1385 = buffer.data(kl + 1385);
    const auto *kl_1386 = buffer.data(kl + 1386);
    const auto *kl_1387 = buffer.data(kl + 1387);
    const auto *kl_1388 = buffer.data(kl + 1388);
    const auto *kl_1389 = buffer.data(kl + 1389);
    const auto *kl_1390 = buffer.data(kl + 1390);
    const auto *kl_1391 = buffer.data(kl + 1391);
    const auto *kl_1392 = buffer.data(kl + 1392);
    const auto *kl_1393 = buffer.data(kl + 1393);
    const auto *kl_1394 = buffer.data(kl + 1394);
    const auto *kl_1395 = buffer.data(kl + 1395);
    const auto *kl_1396 = buffer.data(kl + 1396);
    const auto *kl_1397 = buffer.data(kl + 1397);
    const auto *kl_1398 = buffer.data(kl + 1398);
    const auto *kl_1399 = buffer.data(kl + 1399);
    const auto *kl_1400 = buffer.data(kl + 1400);
    const auto *kl_1401 = buffer.data(kl + 1401);
    const auto *kl_1402 = buffer.data(kl + 1402);
    const auto *kl_1403 = buffer.data(kl + 1403);
    const auto *kl_1404 = buffer.data(kl + 1404);
    const auto *kl_1405 = buffer.data(kl + 1405);
    const auto *kl_1406 = buffer.data(kl + 1406);
    const auto *kl_1407 = buffer.data(kl + 1407);
    const auto *kl_1408 = buffer.data(kl + 1408);
    const auto *kl_1409 = buffer.data(kl + 1409);
    const auto *kl_1410 = buffer.data(kl + 1410);
    const auto *kl_1411 = buffer.data(kl + 1411);
    const auto *kl_1412 = buffer.data(kl + 1412);
    const auto *kl_1413 = buffer.data(kl + 1413);
    const auto *kl_1414 = buffer.data(kl + 1414);
    const auto *kl_1415 = buffer.data(kl + 1415);
    const auto *kl_1416 = buffer.data(kl + 1416);
    const auto *kl_1417 = buffer.data(kl + 1417);
    const auto *kl_1418 = buffer.data(kl + 1418);
    const auto *kl_1419 = buffer.data(kl + 1419);
    const auto *kl_1420 = buffer.data(kl + 1420);
    const auto *kl_1421 = buffer.data(kl + 1421);
    const auto *kl_1422 = buffer.data(kl + 1422);
    const auto *kl_1423 = buffer.data(kl + 1423);
    const auto *kl_1424 = buffer.data(kl + 1424);
    const auto *kl_1425 = buffer.data(kl + 1425);
    const auto *kl_1426 = buffer.data(kl + 1426);
    const auto *kl_1427 = buffer.data(kl + 1427);
    const auto *kl_1428 = buffer.data(kl + 1428);
    const auto *kl_1429 = buffer.data(kl + 1429);
    const auto *kl_1430 = buffer.data(kl + 1430);
    const auto *kl_1431 = buffer.data(kl + 1431);
    const auto *kl_1432 = buffer.data(kl + 1432);
    const auto *kl_1433 = buffer.data(kl + 1433);
    const auto *kl_1434 = buffer.data(kl + 1434);
    const auto *kl_1435 = buffer.data(kl + 1435);
    const auto *kl_1436 = buffer.data(kl + 1436);
    const auto *kl_1437 = buffer.data(kl + 1437);
    const auto *kl_1438 = buffer.data(kl + 1438);
    const auto *kl_1439 = buffer.data(kl + 1439);
    const auto *kl_1440 = buffer.data(kl + 1440);
    const auto *kl_1441 = buffer.data(kl + 1441);
    const auto *kl_1442 = buffer.data(kl + 1442);
    const auto *kl_1443 = buffer.data(kl + 1443);
    const auto *kl_1444 = buffer.data(kl + 1444);
    const auto *kl_1445 = buffer.data(kl + 1445);
    const auto *kl_1446 = buffer.data(kl + 1446);
    const auto *kl_1447 = buffer.data(kl + 1447);
    const auto *kl_1448 = buffer.data(kl + 1448);
    const auto *kl_1449 = buffer.data(kl + 1449);
    const auto *kl_1450 = buffer.data(kl + 1450);
    const auto *kl_1451 = buffer.data(kl + 1451);
    const auto *kl_1452 = buffer.data(kl + 1452);
    const auto *kl_1453 = buffer.data(kl + 1453);
    const auto *kl_1454 = buffer.data(kl + 1454);
    const auto *kl_1455 = buffer.data(kl + 1455);
    const auto *kl_1456 = buffer.data(kl + 1456);
    const auto *kl_1457 = buffer.data(kl + 1457);
    const auto *kl_1458 = buffer.data(kl + 1458);
    const auto *kl_1459 = buffer.data(kl + 1459);
    const auto *kl_1460 = buffer.data(kl + 1460);
    const auto *kl_1461 = buffer.data(kl + 1461);
    const auto *kl_1462 = buffer.data(kl + 1462);
    const auto *kl_1463 = buffer.data(kl + 1463);
    const auto *kl_1464 = buffer.data(kl + 1464);
    const auto *kl_1465 = buffer.data(kl + 1465);
    const auto *kl_1466 = buffer.data(kl + 1466);
    const auto *kl_1467 = buffer.data(kl + 1467);
    const auto *kl_1468 = buffer.data(kl + 1468);
    const auto *kl_1469 = buffer.data(kl + 1469);
    const auto *kl_1470 = buffer.data(kl + 1470);
    const auto *kl_1471 = buffer.data(kl + 1471);
    const auto *kl_1472 = buffer.data(kl + 1472);
    const auto *kl_1473 = buffer.data(kl + 1473);
    const auto *kl_1474 = buffer.data(kl + 1474);
    const auto *kl_1475 = buffer.data(kl + 1475);
    const auto *kl_1476 = buffer.data(kl + 1476);
    const auto *kl_1477 = buffer.data(kl + 1477);
    const auto *kl_1478 = buffer.data(kl + 1478);
    const auto *kl_1479 = buffer.data(kl + 1479);
    const auto *kl_1480 = buffer.data(kl + 1480);
    const auto *kl_1481 = buffer.data(kl + 1481);
    const auto *kl_1482 = buffer.data(kl + 1482);
    const auto *kl_1483 = buffer.data(kl + 1483);
    const auto *kl_1484 = buffer.data(kl + 1484);
    const auto *kl_1485 = buffer.data(kl + 1485);
    const auto *kl_1486 = buffer.data(kl + 1486);
    const auto *kl_1487 = buffer.data(kl + 1487);
    const auto *kl_1488 = buffer.data(kl + 1488);
    const auto *kl_1489 = buffer.data(kl + 1489);
    const auto *kl_1490 = buffer.data(kl + 1490);
    const auto *kl_1491 = buffer.data(kl + 1491);
    const auto *kl_1492 = buffer.data(kl + 1492);
    const auto *kl_1493 = buffer.data(kl + 1493);
    const auto *kl_1494 = buffer.data(kl + 1494);
    const auto *kl_1495 = buffer.data(kl + 1495);
    const auto *kl_1496 = buffer.data(kl + 1496);
    const auto *kl_1497 = buffer.data(kl + 1497);
    const auto *kl_1498 = buffer.data(kl + 1498);
    const auto *kl_1499 = buffer.data(kl + 1499);
    const auto *kl_1500 = buffer.data(kl + 1500);
    const auto *kl_1501 = buffer.data(kl + 1501);
    const auto *kl_1502 = buffer.data(kl + 1502);
    const auto *kl_1503 = buffer.data(kl + 1503);
    const auto *kl_1504 = buffer.data(kl + 1504);
    const auto *kl_1505 = buffer.data(kl + 1505);
    const auto *kl_1506 = buffer.data(kl + 1506);

#pragma omp simd aligned(t_997, t_998, t_999, t_1000, t_1001, hl_682, hl_683, hl_684, hl_685, \
                         hl_686, kl_1357, kl_1358, kl_1359, kl_1360, \
                         kl_1361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_997[k] = -hl_682[k]
                   + f_0 * kl_1357[k];

        t_998[k] = -hl_683[k]
                   + f_0 * kl_1358[k];

        t_999[k] = -hl_684[k]
                   + f_0 * kl_1359[k];

        t_1000[k] = -hl_685[k]
                    + f_0 * kl_1360[k];

        t_1001[k] = -hl_686[k]
                    + f_0 * kl_1361[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, t_1005, t_1006, hl_687, hl_688, hl_689, \
                         hl_690, hl_691, kl_1362, kl_1363, kl_1364, kl_1365, \
                         kl_1366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = -hl_687[k]
                    + f_0 * kl_1362[k];

        t_1003[k] = -hl_688[k]
                    + f_0 * kl_1363[k];

        t_1004[k] = -hl_689[k]
                    + f_0 * kl_1364[k];

        t_1005[k] = -hl_690[k]
                    + f_0 * kl_1365[k];

        t_1006[k] = -hl_691[k]
                    + f_0 * kl_1366[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, t_1010, t_1011, hl_692, hl_693, hl_694, \
                         hl_695, hl_696, kl_1367, kl_1368, kl_1369, kl_1370, \
                         kl_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = -hl_692[k]
                    + f_0 * kl_1367[k];

        t_1008[k] = -hl_693[k]
                    + f_0 * kl_1368[k];

        t_1009[k] = -hl_694[k]
                    + f_0 * kl_1369[k];

        t_1010[k] = -hl_695[k]
                    + f_0 * kl_1370[k];

        t_1011[k] = -hl_696[k]
                    + f_0 * kl_1371[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, t_1015, t_1016, hl_697, hl_698, hl_699, \
                         hl_700, hl_701, kl_1372, kl_1373, kl_1374, kl_1375, \
                         kl_1376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = -hl_697[k]
                    + f_0 * kl_1372[k];

        t_1013[k] = -hl_698[k]
                    + f_0 * kl_1373[k];

        t_1014[k] = -hl_699[k]
                    + f_0 * kl_1374[k];

        t_1015[k] = -hl_700[k]
                    + f_0 * kl_1375[k];

        t_1016[k] = -hl_701[k]
                    + f_0 * kl_1376[k];
    }

#pragma omp simd aligned(t_1017, t_1018, t_1019, t_1020, t_1021, hl_702, hl_703, hl_704, \
                         hl_705, hl_706, kl_1377, kl_1378, kl_1379, kl_1380, \
                         kl_1381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1017[k] = -hl_702[k]
                    + f_0 * kl_1377[k];

        t_1018[k] = -hl_703[k]
                    + f_0 * kl_1378[k];

        t_1019[k] = -hl_704[k]
                    + f_0 * kl_1379[k];

        t_1020[k] = -hl_705[k]
                    + f_0 * kl_1380[k];

        t_1021[k] = -hl_706[k]
                    + f_0 * kl_1381[k];
    }

#pragma omp simd aligned(t_1022, t_1023, t_1024, t_1025, t_1026, hl_707, hl_708, hl_709, \
                         hl_710, hl_711, kl_1382, kl_1383, kl_1384, kl_1385, \
                         kl_1386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1022[k] = -hl_707[k]
                    + f_0 * kl_1382[k];

        t_1023[k] = -hl_708[k]
                    + f_0 * kl_1383[k];

        t_1024[k] = -hl_709[k]
                    + f_0 * kl_1384[k];

        t_1025[k] = -hl_710[k]
                    + f_0 * kl_1385[k];

        t_1026[k] = -hl_711[k]
                    + f_0 * kl_1386[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, t_1031, hl_712, hl_713, hl_714, \
                         hl_715, hl_716, kl_1387, kl_1388, kl_1389, kl_1390, \
                         kl_1391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = -hl_712[k]
                    + f_0 * kl_1387[k];

        t_1028[k] = -hl_713[k]
                    + f_0 * kl_1388[k];

        t_1029[k] = -hl_714[k]
                    + f_0 * kl_1389[k];

        t_1030[k] = -hl_715[k]
                    + f_0 * kl_1390[k];

        t_1031[k] = -hl_716[k]
                    + f_0 * kl_1391[k];
    }

#pragma omp simd aligned(t_1032, t_1033, t_1034, t_1035, t_1036, hl_717, hl_718, hl_719, \
                         hl_720, hl_721, kl_1392, kl_1393, kl_1394, kl_1395, \
                         kl_1396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1032[k] = -hl_717[k]
                    + f_0 * kl_1392[k];

        t_1033[k] = -hl_718[k]
                    + f_0 * kl_1393[k];

        t_1034[k] = -hl_719[k]
                    + f_0 * kl_1394[k];

        t_1035[k] = -2.0 * hl_720[k]
                    + f_0 * kl_1395[k];

        t_1036[k] = -2.0 * hl_721[k]
                    + f_0 * kl_1396[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, t_1040, t_1041, hl_722, hl_723, hl_724, \
                         hl_725, hl_726, kl_1397, kl_1398, kl_1399, kl_1400, \
                         kl_1401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = -2.0 * hl_722[k]
                    + f_0 * kl_1397[k];

        t_1038[k] = -2.0 * hl_723[k]
                    + f_0 * kl_1398[k];

        t_1039[k] = -2.0 * hl_724[k]
                    + f_0 * kl_1399[k];

        t_1040[k] = -2.0 * hl_725[k]
                    + f_0 * kl_1400[k];

        t_1041[k] = -2.0 * hl_726[k]
                    + f_0 * kl_1401[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, t_1046, hl_727, hl_728, hl_729, \
                         hl_730, hl_731, kl_1402, kl_1403, kl_1404, kl_1405, \
                         kl_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = -2.0 * hl_727[k]
                    + f_0 * kl_1402[k];

        t_1043[k] = -2.0 * hl_728[k]
                    + f_0 * kl_1403[k];

        t_1044[k] = -2.0 * hl_729[k]
                    + f_0 * kl_1404[k];

        t_1045[k] = -2.0 * hl_730[k]
                    + f_0 * kl_1405[k];

        t_1046[k] = -2.0 * hl_731[k]
                    + f_0 * kl_1406[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, t_1050, t_1051, hl_732, hl_733, hl_734, \
                         hl_735, hl_736, kl_1407, kl_1408, kl_1409, kl_1410, \
                         kl_1411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = -2.0 * hl_732[k]
                    + f_0 * kl_1407[k];

        t_1048[k] = -2.0 * hl_733[k]
                    + f_0 * kl_1408[k];

        t_1049[k] = -2.0 * hl_734[k]
                    + f_0 * kl_1409[k];

        t_1050[k] = -2.0 * hl_735[k]
                    + f_0 * kl_1410[k];

        t_1051[k] = -2.0 * hl_736[k]
                    + f_0 * kl_1411[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, t_1055, t_1056, hl_737, hl_738, hl_739, \
                         hl_740, hl_741, kl_1412, kl_1413, kl_1414, kl_1415, \
                         kl_1416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = -2.0 * hl_737[k]
                    + f_0 * kl_1412[k];

        t_1053[k] = -2.0 * hl_738[k]
                    + f_0 * kl_1413[k];

        t_1054[k] = -2.0 * hl_739[k]
                    + f_0 * kl_1414[k];

        t_1055[k] = -2.0 * hl_740[k]
                    + f_0 * kl_1415[k];

        t_1056[k] = -2.0 * hl_741[k]
                    + f_0 * kl_1416[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, t_1060, t_1061, hl_742, hl_743, hl_744, \
                         hl_745, hl_746, kl_1417, kl_1418, kl_1419, kl_1420, \
                         kl_1421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = -2.0 * hl_742[k]
                    + f_0 * kl_1417[k];

        t_1058[k] = -2.0 * hl_743[k]
                    + f_0 * kl_1418[k];

        t_1059[k] = -2.0 * hl_744[k]
                    + f_0 * kl_1419[k];

        t_1060[k] = -2.0 * hl_745[k]
                    + f_0 * kl_1420[k];

        t_1061[k] = -2.0 * hl_746[k]
                    + f_0 * kl_1421[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, t_1065, t_1066, hl_747, hl_748, hl_749, \
                         hl_750, hl_751, kl_1422, kl_1423, kl_1424, kl_1425, \
                         kl_1426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = -2.0 * hl_747[k]
                    + f_0 * kl_1422[k];

        t_1063[k] = -2.0 * hl_748[k]
                    + f_0 * kl_1423[k];

        t_1064[k] = -2.0 * hl_749[k]
                    + f_0 * kl_1424[k];

        t_1065[k] = -2.0 * hl_750[k]
                    + f_0 * kl_1425[k];

        t_1066[k] = -2.0 * hl_751[k]
                    + f_0 * kl_1426[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, t_1070, t_1071, hl_752, hl_753, hl_754, \
                         hl_755, hl_756, kl_1427, kl_1428, kl_1429, kl_1430, \
                         kl_1431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = -2.0 * hl_752[k]
                    + f_0 * kl_1427[k];

        t_1068[k] = -2.0 * hl_753[k]
                    + f_0 * kl_1428[k];

        t_1069[k] = -2.0 * hl_754[k]
                    + f_0 * kl_1429[k];

        t_1070[k] = -2.0 * hl_755[k]
                    + f_0 * kl_1430[k];

        t_1071[k] = -2.0 * hl_756[k]
                    + f_0 * kl_1431[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, t_1075, t_1076, hl_757, hl_758, hl_759, \
                         hl_760, hl_761, kl_1432, kl_1433, kl_1434, kl_1435, \
                         kl_1436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = -2.0 * hl_757[k]
                    + f_0 * kl_1432[k];

        t_1073[k] = -2.0 * hl_758[k]
                    + f_0 * kl_1433[k];

        t_1074[k] = -2.0 * hl_759[k]
                    + f_0 * kl_1434[k];

        t_1075[k] = -2.0 * hl_760[k]
                    + f_0 * kl_1435[k];

        t_1076[k] = -2.0 * hl_761[k]
                    + f_0 * kl_1436[k];
    }

#pragma omp simd aligned(t_1077, t_1078, t_1079, t_1080, t_1081, hl_762, hl_763, hl_764, \
                         hl_765, hl_766, kl_1437, kl_1438, kl_1439, kl_1440, \
                         kl_1441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = -2.0 * hl_762[k]
                    + f_0 * kl_1437[k];

        t_1078[k] = -2.0 * hl_763[k]
                    + f_0 * kl_1438[k];

        t_1079[k] = -2.0 * hl_764[k]
                    + f_0 * kl_1439[k];

        t_1080[k] = -3.0 * hl_765[k]
                    + f_0 * kl_1440[k];

        t_1081[k] = -3.0 * hl_766[k]
                    + f_0 * kl_1441[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, t_1085, t_1086, hl_767, hl_768, hl_769, \
                         hl_770, hl_771, kl_1442, kl_1443, kl_1444, kl_1445, \
                         kl_1446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = -3.0 * hl_767[k]
                    + f_0 * kl_1442[k];

        t_1083[k] = -3.0 * hl_768[k]
                    + f_0 * kl_1443[k];

        t_1084[k] = -3.0 * hl_769[k]
                    + f_0 * kl_1444[k];

        t_1085[k] = -3.0 * hl_770[k]
                    + f_0 * kl_1445[k];

        t_1086[k] = -3.0 * hl_771[k]
                    + f_0 * kl_1446[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, t_1090, t_1091, hl_772, hl_773, hl_774, \
                         hl_775, hl_776, kl_1447, kl_1448, kl_1449, kl_1450, \
                         kl_1451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = -3.0 * hl_772[k]
                    + f_0 * kl_1447[k];

        t_1088[k] = -3.0 * hl_773[k]
                    + f_0 * kl_1448[k];

        t_1089[k] = -3.0 * hl_774[k]
                    + f_0 * kl_1449[k];

        t_1090[k] = -3.0 * hl_775[k]
                    + f_0 * kl_1450[k];

        t_1091[k] = -3.0 * hl_776[k]
                    + f_0 * kl_1451[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, t_1095, t_1096, hl_777, hl_778, hl_779, \
                         hl_780, hl_781, kl_1452, kl_1453, kl_1454, kl_1455, \
                         kl_1456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = -3.0 * hl_777[k]
                    + f_0 * kl_1452[k];

        t_1093[k] = -3.0 * hl_778[k]
                    + f_0 * kl_1453[k];

        t_1094[k] = -3.0 * hl_779[k]
                    + f_0 * kl_1454[k];

        t_1095[k] = -3.0 * hl_780[k]
                    + f_0 * kl_1455[k];

        t_1096[k] = -3.0 * hl_781[k]
                    + f_0 * kl_1456[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, t_1100, t_1101, hl_782, hl_783, hl_784, \
                         hl_785, hl_786, kl_1457, kl_1458, kl_1459, kl_1460, \
                         kl_1461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = -3.0 * hl_782[k]
                    + f_0 * kl_1457[k];

        t_1098[k] = -3.0 * hl_783[k]
                    + f_0 * kl_1458[k];

        t_1099[k] = -3.0 * hl_784[k]
                    + f_0 * kl_1459[k];

        t_1100[k] = -3.0 * hl_785[k]
                    + f_0 * kl_1460[k];

        t_1101[k] = -3.0 * hl_786[k]
                    + f_0 * kl_1461[k];
    }

#pragma omp simd aligned(t_1102, t_1103, t_1104, t_1105, t_1106, hl_787, hl_788, hl_789, \
                         hl_790, hl_791, kl_1462, kl_1463, kl_1464, kl_1465, \
                         kl_1466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1102[k] = -3.0 * hl_787[k]
                    + f_0 * kl_1462[k];

        t_1103[k] = -3.0 * hl_788[k]
                    + f_0 * kl_1463[k];

        t_1104[k] = -3.0 * hl_789[k]
                    + f_0 * kl_1464[k];

        t_1105[k] = -3.0 * hl_790[k]
                    + f_0 * kl_1465[k];

        t_1106[k] = -3.0 * hl_791[k]
                    + f_0 * kl_1466[k];
    }

#pragma omp simd aligned(t_1107, t_1108, t_1109, t_1110, t_1111, hl_792, hl_793, hl_794, \
                         hl_795, hl_796, kl_1467, kl_1468, kl_1469, kl_1470, \
                         kl_1471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1107[k] = -3.0 * hl_792[k]
                    + f_0 * kl_1467[k];

        t_1108[k] = -3.0 * hl_793[k]
                    + f_0 * kl_1468[k];

        t_1109[k] = -3.0 * hl_794[k]
                    + f_0 * kl_1469[k];

        t_1110[k] = -3.0 * hl_795[k]
                    + f_0 * kl_1470[k];

        t_1111[k] = -3.0 * hl_796[k]
                    + f_0 * kl_1471[k];
    }

#pragma omp simd aligned(t_1112, t_1113, t_1114, t_1115, t_1116, hl_797, hl_798, hl_799, \
                         hl_800, hl_801, kl_1472, kl_1473, kl_1474, kl_1475, \
                         kl_1476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1112[k] = -3.0 * hl_797[k]
                    + f_0 * kl_1472[k];

        t_1113[k] = -3.0 * hl_798[k]
                    + f_0 * kl_1473[k];

        t_1114[k] = -3.0 * hl_799[k]
                    + f_0 * kl_1474[k];

        t_1115[k] = -3.0 * hl_800[k]
                    + f_0 * kl_1475[k];

        t_1116[k] = -3.0 * hl_801[k]
                    + f_0 * kl_1476[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, t_1120, t_1121, hl_802, hl_803, hl_804, \
                         hl_805, hl_806, kl_1477, kl_1478, kl_1479, kl_1480, \
                         kl_1481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = -3.0 * hl_802[k]
                    + f_0 * kl_1477[k];

        t_1118[k] = -3.0 * hl_803[k]
                    + f_0 * kl_1478[k];

        t_1119[k] = -3.0 * hl_804[k]
                    + f_0 * kl_1479[k];

        t_1120[k] = -3.0 * hl_805[k]
                    + f_0 * kl_1480[k];

        t_1121[k] = -3.0 * hl_806[k]
                    + f_0 * kl_1481[k];
    }

#pragma omp simd aligned(t_1122, t_1123, t_1124, t_1125, t_1126, hl_807, hl_808, hl_809, \
                         hl_810, hl_811, kl_1482, kl_1483, kl_1484, kl_1485, \
                         kl_1486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1122[k] = -3.0 * hl_807[k]
                    + f_0 * kl_1482[k];

        t_1123[k] = -3.0 * hl_808[k]
                    + f_0 * kl_1483[k];

        t_1124[k] = -3.0 * hl_809[k]
                    + f_0 * kl_1484[k];

        t_1125[k] = -4.0 * hl_810[k]
                    + f_0 * kl_1485[k];

        t_1126[k] = -4.0 * hl_811[k]
                    + f_0 * kl_1486[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, t_1130, t_1131, hl_812, hl_813, hl_814, \
                         hl_815, hl_816, kl_1487, kl_1488, kl_1489, kl_1490, \
                         kl_1491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = -4.0 * hl_812[k]
                    + f_0 * kl_1487[k];

        t_1128[k] = -4.0 * hl_813[k]
                    + f_0 * kl_1488[k];

        t_1129[k] = -4.0 * hl_814[k]
                    + f_0 * kl_1489[k];

        t_1130[k] = -4.0 * hl_815[k]
                    + f_0 * kl_1490[k];

        t_1131[k] = -4.0 * hl_816[k]
                    + f_0 * kl_1491[k];
    }

#pragma omp simd aligned(t_1132, t_1133, t_1134, t_1135, t_1136, hl_817, hl_818, hl_819, \
                         hl_820, hl_821, kl_1492, kl_1493, kl_1494, kl_1495, \
                         kl_1496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = -4.0 * hl_817[k]
                    + f_0 * kl_1492[k];

        t_1133[k] = -4.0 * hl_818[k]
                    + f_0 * kl_1493[k];

        t_1134[k] = -4.0 * hl_819[k]
                    + f_0 * kl_1494[k];

        t_1135[k] = -4.0 * hl_820[k]
                    + f_0 * kl_1495[k];

        t_1136[k] = -4.0 * hl_821[k]
                    + f_0 * kl_1496[k];
    }

#pragma omp simd aligned(t_1137, t_1138, t_1139, t_1140, t_1141, hl_822, hl_823, hl_824, \
                         hl_825, hl_826, kl_1497, kl_1498, kl_1499, kl_1500, \
                         kl_1501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1137[k] = -4.0 * hl_822[k]
                    + f_0 * kl_1497[k];

        t_1138[k] = -4.0 * hl_823[k]
                    + f_0 * kl_1498[k];

        t_1139[k] = -4.0 * hl_824[k]
                    + f_0 * kl_1499[k];

        t_1140[k] = -4.0 * hl_825[k]
                    + f_0 * kl_1500[k];

        t_1141[k] = -4.0 * hl_826[k]
                    + f_0 * kl_1501[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, t_1145, t_1146, hl_827, hl_828, hl_829, \
                         hl_830, hl_831, kl_1502, kl_1503, kl_1504, kl_1505, \
                         kl_1506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = -4.0 * hl_827[k]
                    + f_0 * kl_1502[k];

        t_1143[k] = -4.0 * hl_828[k]
                    + f_0 * kl_1503[k];

        t_1144[k] = -4.0 * hl_829[k]
                    + f_0 * kl_1504[k];

        t_1145[k] = -4.0 * hl_830[k]
                    + f_0 * kl_1505[k];

        t_1146[k] = -4.0 * hl_831[k]
                    + f_0 * kl_1506[k];
    }
}

static auto
compute_prim_geom_10_il_electron_repulsion_2_piece7(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hl, const size_t kl,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hl_832 = buffer.data(hl + 832);
    const auto *hl_833 = buffer.data(hl + 833);
    const auto *hl_834 = buffer.data(hl + 834);
    const auto *hl_835 = buffer.data(hl + 835);
    const auto *hl_836 = buffer.data(hl + 836);
    const auto *hl_837 = buffer.data(hl + 837);
    const auto *hl_838 = buffer.data(hl + 838);
    const auto *hl_839 = buffer.data(hl + 839);
    const auto *hl_840 = buffer.data(hl + 840);
    const auto *hl_841 = buffer.data(hl + 841);
    const auto *hl_842 = buffer.data(hl + 842);
    const auto *hl_843 = buffer.data(hl + 843);
    const auto *hl_844 = buffer.data(hl + 844);
    const auto *hl_845 = buffer.data(hl + 845);
    const auto *hl_846 = buffer.data(hl + 846);
    const auto *hl_847 = buffer.data(hl + 847);
    const auto *hl_848 = buffer.data(hl + 848);
    const auto *hl_849 = buffer.data(hl + 849);
    const auto *hl_850 = buffer.data(hl + 850);
    const auto *hl_851 = buffer.data(hl + 851);
    const auto *hl_852 = buffer.data(hl + 852);
    const auto *hl_853 = buffer.data(hl + 853);
    const auto *hl_854 = buffer.data(hl + 854);
    const auto *hl_855 = buffer.data(hl + 855);
    const auto *hl_856 = buffer.data(hl + 856);
    const auto *hl_857 = buffer.data(hl + 857);
    const auto *hl_858 = buffer.data(hl + 858);
    const auto *hl_859 = buffer.data(hl + 859);
    const auto *hl_860 = buffer.data(hl + 860);
    const auto *hl_861 = buffer.data(hl + 861);
    const auto *hl_862 = buffer.data(hl + 862);
    const auto *hl_863 = buffer.data(hl + 863);
    const auto *hl_864 = buffer.data(hl + 864);
    const auto *hl_865 = buffer.data(hl + 865);
    const auto *hl_866 = buffer.data(hl + 866);
    const auto *hl_867 = buffer.data(hl + 867);
    const auto *hl_868 = buffer.data(hl + 868);
    const auto *hl_869 = buffer.data(hl + 869);
    const auto *hl_870 = buffer.data(hl + 870);
    const auto *hl_871 = buffer.data(hl + 871);
    const auto *hl_872 = buffer.data(hl + 872);
    const auto *hl_873 = buffer.data(hl + 873);
    const auto *hl_874 = buffer.data(hl + 874);
    const auto *hl_875 = buffer.data(hl + 875);
    const auto *hl_876 = buffer.data(hl + 876);
    const auto *hl_877 = buffer.data(hl + 877);
    const auto *hl_878 = buffer.data(hl + 878);
    const auto *hl_879 = buffer.data(hl + 879);
    const auto *hl_880 = buffer.data(hl + 880);
    const auto *hl_881 = buffer.data(hl + 881);
    const auto *hl_882 = buffer.data(hl + 882);
    const auto *hl_883 = buffer.data(hl + 883);
    const auto *hl_884 = buffer.data(hl + 884);
    const auto *hl_885 = buffer.data(hl + 885);
    const auto *hl_886 = buffer.data(hl + 886);
    const auto *hl_887 = buffer.data(hl + 887);
    const auto *hl_888 = buffer.data(hl + 888);
    const auto *hl_889 = buffer.data(hl + 889);
    const auto *hl_890 = buffer.data(hl + 890);
    const auto *hl_891 = buffer.data(hl + 891);
    const auto *hl_892 = buffer.data(hl + 892);
    const auto *hl_893 = buffer.data(hl + 893);
    const auto *hl_894 = buffer.data(hl + 894);
    const auto *hl_895 = buffer.data(hl + 895);
    const auto *hl_896 = buffer.data(hl + 896);
    const auto *hl_897 = buffer.data(hl + 897);
    const auto *hl_898 = buffer.data(hl + 898);
    const auto *hl_899 = buffer.data(hl + 899);
    const auto *hl_900 = buffer.data(hl + 900);
    const auto *hl_901 = buffer.data(hl + 901);
    const auto *hl_902 = buffer.data(hl + 902);
    const auto *hl_903 = buffer.data(hl + 903);
    const auto *hl_904 = buffer.data(hl + 904);
    const auto *hl_905 = buffer.data(hl + 905);
    const auto *hl_906 = buffer.data(hl + 906);
    const auto *hl_907 = buffer.data(hl + 907);
    const auto *hl_908 = buffer.data(hl + 908);
    const auto *hl_909 = buffer.data(hl + 909);
    const auto *hl_910 = buffer.data(hl + 910);
    const auto *hl_911 = buffer.data(hl + 911);
    const auto *hl_912 = buffer.data(hl + 912);
    const auto *hl_913 = buffer.data(hl + 913);
    const auto *hl_914 = buffer.data(hl + 914);
    const auto *hl_915 = buffer.data(hl + 915);
    const auto *hl_916 = buffer.data(hl + 916);
    const auto *hl_917 = buffer.data(hl + 917);
    const auto *hl_918 = buffer.data(hl + 918);
    const auto *hl_919 = buffer.data(hl + 919);
    const auto *hl_920 = buffer.data(hl + 920);
    const auto *hl_921 = buffer.data(hl + 921);
    const auto *hl_922 = buffer.data(hl + 922);
    const auto *hl_923 = buffer.data(hl + 923);
    const auto *hl_924 = buffer.data(hl + 924);
    const auto *hl_925 = buffer.data(hl + 925);
    const auto *hl_926 = buffer.data(hl + 926);
    const auto *hl_927 = buffer.data(hl + 927);
    const auto *hl_928 = buffer.data(hl + 928);
    const auto *hl_929 = buffer.data(hl + 929);
    const auto *hl_930 = buffer.data(hl + 930);
    const auto *hl_931 = buffer.data(hl + 931);
    const auto *hl_932 = buffer.data(hl + 932);
    const auto *hl_933 = buffer.data(hl + 933);
    const auto *hl_934 = buffer.data(hl + 934);
    const auto *hl_935 = buffer.data(hl + 935);
    const auto *hl_936 = buffer.data(hl + 936);
    const auto *hl_937 = buffer.data(hl + 937);
    const auto *hl_938 = buffer.data(hl + 938);
    const auto *hl_939 = buffer.data(hl + 939);
    const auto *hl_940 = buffer.data(hl + 940);
    const auto *hl_941 = buffer.data(hl + 941);
    const auto *hl_942 = buffer.data(hl + 942);
    const auto *hl_943 = buffer.data(hl + 943);
    const auto *hl_944 = buffer.data(hl + 944);

    const auto *kl_1507 = buffer.data(kl + 1507);
    const auto *kl_1508 = buffer.data(kl + 1508);
    const auto *kl_1509 = buffer.data(kl + 1509);
    const auto *kl_1510 = buffer.data(kl + 1510);
    const auto *kl_1511 = buffer.data(kl + 1511);
    const auto *kl_1512 = buffer.data(kl + 1512);
    const auto *kl_1513 = buffer.data(kl + 1513);
    const auto *kl_1514 = buffer.data(kl + 1514);
    const auto *kl_1515 = buffer.data(kl + 1515);
    const auto *kl_1516 = buffer.data(kl + 1516);
    const auto *kl_1517 = buffer.data(kl + 1517);
    const auto *kl_1518 = buffer.data(kl + 1518);
    const auto *kl_1519 = buffer.data(kl + 1519);
    const auto *kl_1520 = buffer.data(kl + 1520);
    const auto *kl_1521 = buffer.data(kl + 1521);
    const auto *kl_1522 = buffer.data(kl + 1522);
    const auto *kl_1523 = buffer.data(kl + 1523);
    const auto *kl_1524 = buffer.data(kl + 1524);
    const auto *kl_1525 = buffer.data(kl + 1525);
    const auto *kl_1526 = buffer.data(kl + 1526);
    const auto *kl_1527 = buffer.data(kl + 1527);
    const auto *kl_1528 = buffer.data(kl + 1528);
    const auto *kl_1529 = buffer.data(kl + 1529);
    const auto *kl_1530 = buffer.data(kl + 1530);
    const auto *kl_1531 = buffer.data(kl + 1531);
    const auto *kl_1532 = buffer.data(kl + 1532);
    const auto *kl_1533 = buffer.data(kl + 1533);
    const auto *kl_1534 = buffer.data(kl + 1534);
    const auto *kl_1535 = buffer.data(kl + 1535);
    const auto *kl_1536 = buffer.data(kl + 1536);
    const auto *kl_1537 = buffer.data(kl + 1537);
    const auto *kl_1538 = buffer.data(kl + 1538);
    const auto *kl_1539 = buffer.data(kl + 1539);
    const auto *kl_1540 = buffer.data(kl + 1540);
    const auto *kl_1541 = buffer.data(kl + 1541);
    const auto *kl_1542 = buffer.data(kl + 1542);
    const auto *kl_1543 = buffer.data(kl + 1543);
    const auto *kl_1544 = buffer.data(kl + 1544);
    const auto *kl_1545 = buffer.data(kl + 1545);
    const auto *kl_1546 = buffer.data(kl + 1546);
    const auto *kl_1547 = buffer.data(kl + 1547);
    const auto *kl_1548 = buffer.data(kl + 1548);
    const auto *kl_1549 = buffer.data(kl + 1549);
    const auto *kl_1550 = buffer.data(kl + 1550);
    const auto *kl_1551 = buffer.data(kl + 1551);
    const auto *kl_1552 = buffer.data(kl + 1552);
    const auto *kl_1553 = buffer.data(kl + 1553);
    const auto *kl_1554 = buffer.data(kl + 1554);
    const auto *kl_1555 = buffer.data(kl + 1555);
    const auto *kl_1556 = buffer.data(kl + 1556);
    const auto *kl_1557 = buffer.data(kl + 1557);
    const auto *kl_1558 = buffer.data(kl + 1558);
    const auto *kl_1559 = buffer.data(kl + 1559);
    const auto *kl_1560 = buffer.data(kl + 1560);
    const auto *kl_1561 = buffer.data(kl + 1561);
    const auto *kl_1562 = buffer.data(kl + 1562);
    const auto *kl_1563 = buffer.data(kl + 1563);
    const auto *kl_1564 = buffer.data(kl + 1564);
    const auto *kl_1565 = buffer.data(kl + 1565);
    const auto *kl_1566 = buffer.data(kl + 1566);
    const auto *kl_1567 = buffer.data(kl + 1567);
    const auto *kl_1568 = buffer.data(kl + 1568);
    const auto *kl_1569 = buffer.data(kl + 1569);
    const auto *kl_1570 = buffer.data(kl + 1570);
    const auto *kl_1571 = buffer.data(kl + 1571);
    const auto *kl_1572 = buffer.data(kl + 1572);
    const auto *kl_1573 = buffer.data(kl + 1573);
    const auto *kl_1574 = buffer.data(kl + 1574);
    const auto *kl_1575 = buffer.data(kl + 1575);
    const auto *kl_1576 = buffer.data(kl + 1576);
    const auto *kl_1577 = buffer.data(kl + 1577);
    const auto *kl_1578 = buffer.data(kl + 1578);
    const auto *kl_1579 = buffer.data(kl + 1579);
    const auto *kl_1580 = buffer.data(kl + 1580);
    const auto *kl_1581 = buffer.data(kl + 1581);
    const auto *kl_1582 = buffer.data(kl + 1582);
    const auto *kl_1583 = buffer.data(kl + 1583);
    const auto *kl_1584 = buffer.data(kl + 1584);
    const auto *kl_1585 = buffer.data(kl + 1585);
    const auto *kl_1586 = buffer.data(kl + 1586);
    const auto *kl_1587 = buffer.data(kl + 1587);
    const auto *kl_1588 = buffer.data(kl + 1588);
    const auto *kl_1589 = buffer.data(kl + 1589);
    const auto *kl_1590 = buffer.data(kl + 1590);
    const auto *kl_1591 = buffer.data(kl + 1591);
    const auto *kl_1592 = buffer.data(kl + 1592);
    const auto *kl_1593 = buffer.data(kl + 1593);
    const auto *kl_1594 = buffer.data(kl + 1594);
    const auto *kl_1595 = buffer.data(kl + 1595);
    const auto *kl_1596 = buffer.data(kl + 1596);
    const auto *kl_1597 = buffer.data(kl + 1597);
    const auto *kl_1598 = buffer.data(kl + 1598);
    const auto *kl_1599 = buffer.data(kl + 1599);
    const auto *kl_1600 = buffer.data(kl + 1600);
    const auto *kl_1601 = buffer.data(kl + 1601);
    const auto *kl_1602 = buffer.data(kl + 1602);
    const auto *kl_1603 = buffer.data(kl + 1603);
    const auto *kl_1604 = buffer.data(kl + 1604);
    const auto *kl_1605 = buffer.data(kl + 1605);
    const auto *kl_1606 = buffer.data(kl + 1606);
    const auto *kl_1607 = buffer.data(kl + 1607);
    const auto *kl_1608 = buffer.data(kl + 1608);
    const auto *kl_1609 = buffer.data(kl + 1609);
    const auto *kl_1610 = buffer.data(kl + 1610);
    const auto *kl_1611 = buffer.data(kl + 1611);
    const auto *kl_1612 = buffer.data(kl + 1612);
    const auto *kl_1613 = buffer.data(kl + 1613);
    const auto *kl_1614 = buffer.data(kl + 1614);
    const auto *kl_1615 = buffer.data(kl + 1615);
    const auto *kl_1616 = buffer.data(kl + 1616);
    const auto *kl_1617 = buffer.data(kl + 1617);
    const auto *kl_1618 = buffer.data(kl + 1618);
    const auto *kl_1619 = buffer.data(kl + 1619);

#pragma omp simd aligned(t_1147, t_1148, t_1149, t_1150, t_1151, hl_832, hl_833, hl_834, \
                         hl_835, hl_836, kl_1507, kl_1508, kl_1509, kl_1510, \
                         kl_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = -4.0 * hl_832[k]
                    + f_0 * kl_1507[k];

        t_1148[k] = -4.0 * hl_833[k]
                    + f_0 * kl_1508[k];

        t_1149[k] = -4.0 * hl_834[k]
                    + f_0 * kl_1509[k];

        t_1150[k] = -4.0 * hl_835[k]
                    + f_0 * kl_1510[k];

        t_1151[k] = -4.0 * hl_836[k]
                    + f_0 * kl_1511[k];
    }

#pragma omp simd aligned(t_1152, t_1153, t_1154, t_1155, t_1156, hl_837, hl_838, hl_839, \
                         hl_840, hl_841, kl_1512, kl_1513, kl_1514, kl_1515, \
                         kl_1516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1152[k] = -4.0 * hl_837[k]
                    + f_0 * kl_1512[k];

        t_1153[k] = -4.0 * hl_838[k]
                    + f_0 * kl_1513[k];

        t_1154[k] = -4.0 * hl_839[k]
                    + f_0 * kl_1514[k];

        t_1155[k] = -4.0 * hl_840[k]
                    + f_0 * kl_1515[k];

        t_1156[k] = -4.0 * hl_841[k]
                    + f_0 * kl_1516[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, t_1160, t_1161, hl_842, hl_843, hl_844, \
                         hl_845, hl_846, kl_1517, kl_1518, kl_1519, kl_1520, \
                         kl_1521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = -4.0 * hl_842[k]
                    + f_0 * kl_1517[k];

        t_1158[k] = -4.0 * hl_843[k]
                    + f_0 * kl_1518[k];

        t_1159[k] = -4.0 * hl_844[k]
                    + f_0 * kl_1519[k];

        t_1160[k] = -4.0 * hl_845[k]
                    + f_0 * kl_1520[k];

        t_1161[k] = -4.0 * hl_846[k]
                    + f_0 * kl_1521[k];
    }

#pragma omp simd aligned(t_1162, t_1163, t_1164, t_1165, t_1166, hl_847, hl_848, hl_849, \
                         hl_850, hl_851, kl_1522, kl_1523, kl_1524, kl_1525, \
                         kl_1526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1162[k] = -4.0 * hl_847[k]
                    + f_0 * kl_1522[k];

        t_1163[k] = -4.0 * hl_848[k]
                    + f_0 * kl_1523[k];

        t_1164[k] = -4.0 * hl_849[k]
                    + f_0 * kl_1524[k];

        t_1165[k] = -4.0 * hl_850[k]
                    + f_0 * kl_1525[k];

        t_1166[k] = -4.0 * hl_851[k]
                    + f_0 * kl_1526[k];
    }

#pragma omp simd aligned(t_1167, t_1168, t_1169, t_1170, t_1171, hl_852, hl_853, hl_854, \
                         hl_855, hl_856, kl_1527, kl_1528, kl_1529, kl_1530, \
                         kl_1531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1167[k] = -4.0 * hl_852[k]
                    + f_0 * kl_1527[k];

        t_1168[k] = -4.0 * hl_853[k]
                    + f_0 * kl_1528[k];

        t_1169[k] = -4.0 * hl_854[k]
                    + f_0 * kl_1529[k];

        t_1170[k] = -5.0 * hl_855[k]
                    + f_0 * kl_1530[k];

        t_1171[k] = -5.0 * hl_856[k]
                    + f_0 * kl_1531[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, t_1175, t_1176, hl_857, hl_858, hl_859, \
                         hl_860, hl_861, kl_1532, kl_1533, kl_1534, kl_1535, \
                         kl_1536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = -5.0 * hl_857[k]
                    + f_0 * kl_1532[k];

        t_1173[k] = -5.0 * hl_858[k]
                    + f_0 * kl_1533[k];

        t_1174[k] = -5.0 * hl_859[k]
                    + f_0 * kl_1534[k];

        t_1175[k] = -5.0 * hl_860[k]
                    + f_0 * kl_1535[k];

        t_1176[k] = -5.0 * hl_861[k]
                    + f_0 * kl_1536[k];
    }

#pragma omp simd aligned(t_1177, t_1178, t_1179, t_1180, t_1181, hl_862, hl_863, hl_864, \
                         hl_865, hl_866, kl_1537, kl_1538, kl_1539, kl_1540, \
                         kl_1541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1177[k] = -5.0 * hl_862[k]
                    + f_0 * kl_1537[k];

        t_1178[k] = -5.0 * hl_863[k]
                    + f_0 * kl_1538[k];

        t_1179[k] = -5.0 * hl_864[k]
                    + f_0 * kl_1539[k];

        t_1180[k] = -5.0 * hl_865[k]
                    + f_0 * kl_1540[k];

        t_1181[k] = -5.0 * hl_866[k]
                    + f_0 * kl_1541[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, t_1185, t_1186, hl_867, hl_868, hl_869, \
                         hl_870, hl_871, kl_1542, kl_1543, kl_1544, kl_1545, \
                         kl_1546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = -5.0 * hl_867[k]
                    + f_0 * kl_1542[k];

        t_1183[k] = -5.0 * hl_868[k]
                    + f_0 * kl_1543[k];

        t_1184[k] = -5.0 * hl_869[k]
                    + f_0 * kl_1544[k];

        t_1185[k] = -5.0 * hl_870[k]
                    + f_0 * kl_1545[k];

        t_1186[k] = -5.0 * hl_871[k]
                    + f_0 * kl_1546[k];
    }

#pragma omp simd aligned(t_1187, t_1188, t_1189, t_1190, t_1191, hl_872, hl_873, hl_874, \
                         hl_875, hl_876, kl_1547, kl_1548, kl_1549, kl_1550, \
                         kl_1551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1187[k] = -5.0 * hl_872[k]
                    + f_0 * kl_1547[k];

        t_1188[k] = -5.0 * hl_873[k]
                    + f_0 * kl_1548[k];

        t_1189[k] = -5.0 * hl_874[k]
                    + f_0 * kl_1549[k];

        t_1190[k] = -5.0 * hl_875[k]
                    + f_0 * kl_1550[k];

        t_1191[k] = -5.0 * hl_876[k]
                    + f_0 * kl_1551[k];
    }

#pragma omp simd aligned(t_1192, t_1193, t_1194, t_1195, t_1196, hl_877, hl_878, hl_879, \
                         hl_880, hl_881, kl_1552, kl_1553, kl_1554, kl_1555, \
                         kl_1556 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1192[k] = -5.0 * hl_877[k]
                    + f_0 * kl_1552[k];

        t_1193[k] = -5.0 * hl_878[k]
                    + f_0 * kl_1553[k];

        t_1194[k] = -5.0 * hl_879[k]
                    + f_0 * kl_1554[k];

        t_1195[k] = -5.0 * hl_880[k]
                    + f_0 * kl_1555[k];

        t_1196[k] = -5.0 * hl_881[k]
                    + f_0 * kl_1556[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, t_1200, t_1201, hl_882, hl_883, hl_884, \
                         hl_885, hl_886, kl_1557, kl_1558, kl_1559, kl_1560, \
                         kl_1561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = -5.0 * hl_882[k]
                    + f_0 * kl_1557[k];

        t_1198[k] = -5.0 * hl_883[k]
                    + f_0 * kl_1558[k];

        t_1199[k] = -5.0 * hl_884[k]
                    + f_0 * kl_1559[k];

        t_1200[k] = -5.0 * hl_885[k]
                    + f_0 * kl_1560[k];

        t_1201[k] = -5.0 * hl_886[k]
                    + f_0 * kl_1561[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, t_1205, t_1206, hl_887, hl_888, hl_889, \
                         hl_890, hl_891, kl_1562, kl_1563, kl_1564, kl_1565, \
                         kl_1566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = -5.0 * hl_887[k]
                    + f_0 * kl_1562[k];

        t_1203[k] = -5.0 * hl_888[k]
                    + f_0 * kl_1563[k];

        t_1204[k] = -5.0 * hl_889[k]
                    + f_0 * kl_1564[k];

        t_1205[k] = -5.0 * hl_890[k]
                    + f_0 * kl_1565[k];

        t_1206[k] = -5.0 * hl_891[k]
                    + f_0 * kl_1566[k];
    }

#pragma omp simd aligned(t_1207, t_1208, t_1209, t_1210, t_1211, hl_892, hl_893, hl_894, \
                         hl_895, hl_896, kl_1567, kl_1568, kl_1569, kl_1570, \
                         kl_1571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1207[k] = -5.0 * hl_892[k]
                    + f_0 * kl_1567[k];

        t_1208[k] = -5.0 * hl_893[k]
                    + f_0 * kl_1568[k];

        t_1209[k] = -5.0 * hl_894[k]
                    + f_0 * kl_1569[k];

        t_1210[k] = -5.0 * hl_895[k]
                    + f_0 * kl_1570[k];

        t_1211[k] = -5.0 * hl_896[k]
                    + f_0 * kl_1571[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, t_1215, t_1216, hl_897, hl_898, hl_899, \
                         hl_900, hl_901, kl_1572, kl_1573, kl_1574, kl_1575, \
                         kl_1576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = -5.0 * hl_897[k]
                    + f_0 * kl_1572[k];

        t_1213[k] = -5.0 * hl_898[k]
                    + f_0 * kl_1573[k];

        t_1214[k] = -5.0 * hl_899[k]
                    + f_0 * kl_1574[k];

        t_1215[k] = -6.0 * hl_900[k]
                    + f_0 * kl_1575[k];

        t_1216[k] = -6.0 * hl_901[k]
                    + f_0 * kl_1576[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, t_1220, t_1221, hl_902, hl_903, hl_904, \
                         hl_905, hl_906, kl_1577, kl_1578, kl_1579, kl_1580, \
                         kl_1581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = -6.0 * hl_902[k]
                    + f_0 * kl_1577[k];

        t_1218[k] = -6.0 * hl_903[k]
                    + f_0 * kl_1578[k];

        t_1219[k] = -6.0 * hl_904[k]
                    + f_0 * kl_1579[k];

        t_1220[k] = -6.0 * hl_905[k]
                    + f_0 * kl_1580[k];

        t_1221[k] = -6.0 * hl_906[k]
                    + f_0 * kl_1581[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, t_1226, hl_907, hl_908, hl_909, \
                         hl_910, hl_911, kl_1582, kl_1583, kl_1584, kl_1585, \
                         kl_1586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = -6.0 * hl_907[k]
                    + f_0 * kl_1582[k];

        t_1223[k] = -6.0 * hl_908[k]
                    + f_0 * kl_1583[k];

        t_1224[k] = -6.0 * hl_909[k]
                    + f_0 * kl_1584[k];

        t_1225[k] = -6.0 * hl_910[k]
                    + f_0 * kl_1585[k];

        t_1226[k] = -6.0 * hl_911[k]
                    + f_0 * kl_1586[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, t_1230, t_1231, hl_912, hl_913, hl_914, \
                         hl_915, hl_916, kl_1587, kl_1588, kl_1589, kl_1590, \
                         kl_1591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = -6.0 * hl_912[k]
                    + f_0 * kl_1587[k];

        t_1228[k] = -6.0 * hl_913[k]
                    + f_0 * kl_1588[k];

        t_1229[k] = -6.0 * hl_914[k]
                    + f_0 * kl_1589[k];

        t_1230[k] = -6.0 * hl_915[k]
                    + f_0 * kl_1590[k];

        t_1231[k] = -6.0 * hl_916[k]
                    + f_0 * kl_1591[k];
    }

#pragma omp simd aligned(t_1232, t_1233, t_1234, t_1235, t_1236, hl_917, hl_918, hl_919, \
                         hl_920, hl_921, kl_1592, kl_1593, kl_1594, kl_1595, \
                         kl_1596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1232[k] = -6.0 * hl_917[k]
                    + f_0 * kl_1592[k];

        t_1233[k] = -6.0 * hl_918[k]
                    + f_0 * kl_1593[k];

        t_1234[k] = -6.0 * hl_919[k]
                    + f_0 * kl_1594[k];

        t_1235[k] = -6.0 * hl_920[k]
                    + f_0 * kl_1595[k];

        t_1236[k] = -6.0 * hl_921[k]
                    + f_0 * kl_1596[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, t_1241, hl_922, hl_923, hl_924, \
                         hl_925, hl_926, kl_1597, kl_1598, kl_1599, kl_1600, \
                         kl_1601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = -6.0 * hl_922[k]
                    + f_0 * kl_1597[k];

        t_1238[k] = -6.0 * hl_923[k]
                    + f_0 * kl_1598[k];

        t_1239[k] = -6.0 * hl_924[k]
                    + f_0 * kl_1599[k];

        t_1240[k] = -6.0 * hl_925[k]
                    + f_0 * kl_1600[k];

        t_1241[k] = -6.0 * hl_926[k]
                    + f_0 * kl_1601[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, t_1245, t_1246, hl_927, hl_928, hl_929, \
                         hl_930, hl_931, kl_1602, kl_1603, kl_1604, kl_1605, \
                         kl_1606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = -6.0 * hl_927[k]
                    + f_0 * kl_1602[k];

        t_1243[k] = -6.0 * hl_928[k]
                    + f_0 * kl_1603[k];

        t_1244[k] = -6.0 * hl_929[k]
                    + f_0 * kl_1604[k];

        t_1245[k] = -6.0 * hl_930[k]
                    + f_0 * kl_1605[k];

        t_1246[k] = -6.0 * hl_931[k]
                    + f_0 * kl_1606[k];
    }

#pragma omp simd aligned(t_1247, t_1248, t_1249, t_1250, t_1251, hl_932, hl_933, hl_934, \
                         hl_935, hl_936, kl_1607, kl_1608, kl_1609, kl_1610, \
                         kl_1611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1247[k] = -6.0 * hl_932[k]
                    + f_0 * kl_1607[k];

        t_1248[k] = -6.0 * hl_933[k]
                    + f_0 * kl_1608[k];

        t_1249[k] = -6.0 * hl_934[k]
                    + f_0 * kl_1609[k];

        t_1250[k] = -6.0 * hl_935[k]
                    + f_0 * kl_1610[k];

        t_1251[k] = -6.0 * hl_936[k]
                    + f_0 * kl_1611[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, t_1256, hl_937, hl_938, hl_939, \
                         hl_940, hl_941, kl_1612, kl_1613, kl_1614, kl_1615, \
                         kl_1616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = -6.0 * hl_937[k]
                    + f_0 * kl_1612[k];

        t_1253[k] = -6.0 * hl_938[k]
                    + f_0 * kl_1613[k];

        t_1254[k] = -6.0 * hl_939[k]
                    + f_0 * kl_1614[k];

        t_1255[k] = -6.0 * hl_940[k]
                    + f_0 * kl_1615[k];

        t_1256[k] = -6.0 * hl_941[k]
                    + f_0 * kl_1616[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, hl_942, hl_943, hl_944, kl_1617, kl_1618, \
                         kl_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = -6.0 * hl_942[k]
                    + f_0 * kl_1617[k];

        t_1258[k] = -6.0 * hl_943[k]
                    + f_0 * kl_1618[k];

        t_1259[k] = -6.0 * hl_944[k]
                    + f_0 * kl_1619[k];
    }
}

auto
compute_prim_geom_10_il_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t hl, const size_t kl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_il_electron_repulsion_2_piece0(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_2_piece1(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_2_piece2(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_2_piece3(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_2_piece4(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_2_piece5(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_2_piece6(buffer, target, hl, kl, ncols, alpha);

    compute_prim_geom_10_il_electron_repulsion_2_piece7(buffer, target, hl, kl, ncols, alpha);
}

}  // namespace simdt2ceri
