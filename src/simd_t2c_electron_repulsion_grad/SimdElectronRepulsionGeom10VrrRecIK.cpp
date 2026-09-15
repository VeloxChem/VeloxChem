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


#include "SimdElectronRepulsionGeom10VrrRecIK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_ik_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_1 = buffer.data(hk + 1);
    const auto *hk_2 = buffer.data(hk + 2);
    const auto *hk_3 = buffer.data(hk + 3);
    const auto *hk_4 = buffer.data(hk + 4);
    const auto *hk_5 = buffer.data(hk + 5);
    const auto *hk_6 = buffer.data(hk + 6);
    const auto *hk_7 = buffer.data(hk + 7);
    const auto *hk_8 = buffer.data(hk + 8);
    const auto *hk_9 = buffer.data(hk + 9);
    const auto *hk_10 = buffer.data(hk + 10);
    const auto *hk_11 = buffer.data(hk + 11);
    const auto *hk_12 = buffer.data(hk + 12);
    const auto *hk_13 = buffer.data(hk + 13);
    const auto *hk_14 = buffer.data(hk + 14);
    const auto *hk_15 = buffer.data(hk + 15);
    const auto *hk_16 = buffer.data(hk + 16);
    const auto *hk_17 = buffer.data(hk + 17);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_21 = buffer.data(hk + 21);
    const auto *hk_22 = buffer.data(hk + 22);
    const auto *hk_23 = buffer.data(hk + 23);
    const auto *hk_24 = buffer.data(hk + 24);
    const auto *hk_25 = buffer.data(hk + 25);
    const auto *hk_26 = buffer.data(hk + 26);
    const auto *hk_27 = buffer.data(hk + 27);
    const auto *hk_28 = buffer.data(hk + 28);
    const auto *hk_29 = buffer.data(hk + 29);
    const auto *hk_30 = buffer.data(hk + 30);
    const auto *hk_31 = buffer.data(hk + 31);
    const auto *hk_32 = buffer.data(hk + 32);
    const auto *hk_33 = buffer.data(hk + 33);
    const auto *hk_34 = buffer.data(hk + 34);
    const auto *hk_35 = buffer.data(hk + 35);
    const auto *hk_36 = buffer.data(hk + 36);
    const auto *hk_37 = buffer.data(hk + 37);
    const auto *hk_38 = buffer.data(hk + 38);
    const auto *hk_39 = buffer.data(hk + 39);
    const auto *hk_40 = buffer.data(hk + 40);
    const auto *hk_41 = buffer.data(hk + 41);
    const auto *hk_42 = buffer.data(hk + 42);
    const auto *hk_43 = buffer.data(hk + 43);
    const auto *hk_44 = buffer.data(hk + 44);
    const auto *hk_45 = buffer.data(hk + 45);
    const auto *hk_46 = buffer.data(hk + 46);
    const auto *hk_47 = buffer.data(hk + 47);
    const auto *hk_48 = buffer.data(hk + 48);
    const auto *hk_49 = buffer.data(hk + 49);
    const auto *hk_50 = buffer.data(hk + 50);
    const auto *hk_51 = buffer.data(hk + 51);
    const auto *hk_52 = buffer.data(hk + 52);
    const auto *hk_53 = buffer.data(hk + 53);
    const auto *hk_54 = buffer.data(hk + 54);
    const auto *hk_55 = buffer.data(hk + 55);
    const auto *hk_56 = buffer.data(hk + 56);
    const auto *hk_57 = buffer.data(hk + 57);
    const auto *hk_58 = buffer.data(hk + 58);
    const auto *hk_59 = buffer.data(hk + 59);
    const auto *hk_60 = buffer.data(hk + 60);
    const auto *hk_61 = buffer.data(hk + 61);
    const auto *hk_62 = buffer.data(hk + 62);
    const auto *hk_63 = buffer.data(hk + 63);
    const auto *hk_64 = buffer.data(hk + 64);
    const auto *hk_65 = buffer.data(hk + 65);
    const auto *hk_66 = buffer.data(hk + 66);
    const auto *hk_67 = buffer.data(hk + 67);
    const auto *hk_68 = buffer.data(hk + 68);
    const auto *hk_69 = buffer.data(hk + 69);
    const auto *hk_70 = buffer.data(hk + 70);
    const auto *hk_71 = buffer.data(hk + 71);
    const auto *hk_72 = buffer.data(hk + 72);
    const auto *hk_73 = buffer.data(hk + 73);
    const auto *hk_74 = buffer.data(hk + 74);
    const auto *hk_75 = buffer.data(hk + 75);
    const auto *hk_76 = buffer.data(hk + 76);
    const auto *hk_77 = buffer.data(hk + 77);
    const auto *hk_78 = buffer.data(hk + 78);
    const auto *hk_79 = buffer.data(hk + 79);
    const auto *hk_80 = buffer.data(hk + 80);
    const auto *hk_81 = buffer.data(hk + 81);
    const auto *hk_82 = buffer.data(hk + 82);
    const auto *hk_83 = buffer.data(hk + 83);
    const auto *hk_84 = buffer.data(hk + 84);
    const auto *hk_85 = buffer.data(hk + 85);
    const auto *hk_86 = buffer.data(hk + 86);
    const auto *hk_87 = buffer.data(hk + 87);
    const auto *hk_88 = buffer.data(hk + 88);
    const auto *hk_89 = buffer.data(hk + 89);
    const auto *hk_90 = buffer.data(hk + 90);
    const auto *hk_91 = buffer.data(hk + 91);
    const auto *hk_92 = buffer.data(hk + 92);
    const auto *hk_93 = buffer.data(hk + 93);
    const auto *hk_94 = buffer.data(hk + 94);
    const auto *hk_95 = buffer.data(hk + 95);
    const auto *hk_96 = buffer.data(hk + 96);
    const auto *hk_97 = buffer.data(hk + 97);
    const auto *hk_98 = buffer.data(hk + 98);
    const auto *hk_99 = buffer.data(hk + 99);
    const auto *hk_100 = buffer.data(hk + 100);
    const auto *hk_101 = buffer.data(hk + 101);
    const auto *hk_102 = buffer.data(hk + 102);
    const auto *hk_103 = buffer.data(hk + 103);
    const auto *hk_104 = buffer.data(hk + 104);
    const auto *hk_105 = buffer.data(hk + 105);
    const auto *hk_106 = buffer.data(hk + 106);
    const auto *hk_107 = buffer.data(hk + 107);
    const auto *hk_108 = buffer.data(hk + 108);
    const auto *hk_109 = buffer.data(hk + 109);
    const auto *hk_110 = buffer.data(hk + 110);
    const auto *hk_111 = buffer.data(hk + 111);
    const auto *hk_112 = buffer.data(hk + 112);
    const auto *hk_113 = buffer.data(hk + 113);
    const auto *hk_114 = buffer.data(hk + 114);
    const auto *hk_115 = buffer.data(hk + 115);
    const auto *hk_116 = buffer.data(hk + 116);
    const auto *hk_117 = buffer.data(hk + 117);
    const auto *hk_118 = buffer.data(hk + 118);
    const auto *hk_119 = buffer.data(hk + 119);
    const auto *hk_120 = buffer.data(hk + 120);
    const auto *hk_121 = buffer.data(hk + 121);
    const auto *hk_122 = buffer.data(hk + 122);
    const auto *hk_123 = buffer.data(hk + 123);
    const auto *hk_124 = buffer.data(hk + 124);
    const auto *hk_125 = buffer.data(hk + 125);
    const auto *hk_126 = buffer.data(hk + 126);
    const auto *hk_127 = buffer.data(hk + 127);
    const auto *hk_128 = buffer.data(hk + 128);
    const auto *hk_129 = buffer.data(hk + 129);
    const auto *hk_130 = buffer.data(hk + 130);
    const auto *hk_131 = buffer.data(hk + 131);
    const auto *hk_132 = buffer.data(hk + 132);
    const auto *hk_133 = buffer.data(hk + 133);
    const auto *hk_134 = buffer.data(hk + 134);
    const auto *hk_135 = buffer.data(hk + 135);
    const auto *hk_136 = buffer.data(hk + 136);
    const auto *hk_137 = buffer.data(hk + 137);
    const auto *hk_138 = buffer.data(hk + 138);
    const auto *hk_139 = buffer.data(hk + 139);
    const auto *hk_140 = buffer.data(hk + 140);
    const auto *hk_141 = buffer.data(hk + 141);
    const auto *hk_142 = buffer.data(hk + 142);
    const auto *hk_143 = buffer.data(hk + 143);
    const auto *hk_144 = buffer.data(hk + 144);
    const auto *hk_145 = buffer.data(hk + 145);
    const auto *hk_146 = buffer.data(hk + 146);
    const auto *hk_147 = buffer.data(hk + 147);
    const auto *hk_148 = buffer.data(hk + 148);
    const auto *hk_149 = buffer.data(hk + 149);

    const auto *kk_0 = buffer.data(kk + 0);
    const auto *kk_1 = buffer.data(kk + 1);
    const auto *kk_2 = buffer.data(kk + 2);
    const auto *kk_3 = buffer.data(kk + 3);
    const auto *kk_4 = buffer.data(kk + 4);
    const auto *kk_5 = buffer.data(kk + 5);
    const auto *kk_6 = buffer.data(kk + 6);
    const auto *kk_7 = buffer.data(kk + 7);
    const auto *kk_8 = buffer.data(kk + 8);
    const auto *kk_9 = buffer.data(kk + 9);
    const auto *kk_10 = buffer.data(kk + 10);
    const auto *kk_11 = buffer.data(kk + 11);
    const auto *kk_12 = buffer.data(kk + 12);
    const auto *kk_13 = buffer.data(kk + 13);
    const auto *kk_14 = buffer.data(kk + 14);
    const auto *kk_15 = buffer.data(kk + 15);
    const auto *kk_16 = buffer.data(kk + 16);
    const auto *kk_17 = buffer.data(kk + 17);
    const auto *kk_18 = buffer.data(kk + 18);
    const auto *kk_19 = buffer.data(kk + 19);
    const auto *kk_20 = buffer.data(kk + 20);
    const auto *kk_21 = buffer.data(kk + 21);
    const auto *kk_22 = buffer.data(kk + 22);
    const auto *kk_23 = buffer.data(kk + 23);
    const auto *kk_24 = buffer.data(kk + 24);
    const auto *kk_25 = buffer.data(kk + 25);
    const auto *kk_26 = buffer.data(kk + 26);
    const auto *kk_27 = buffer.data(kk + 27);
    const auto *kk_28 = buffer.data(kk + 28);
    const auto *kk_29 = buffer.data(kk + 29);
    const auto *kk_30 = buffer.data(kk + 30);
    const auto *kk_31 = buffer.data(kk + 31);
    const auto *kk_32 = buffer.data(kk + 32);
    const auto *kk_33 = buffer.data(kk + 33);
    const auto *kk_34 = buffer.data(kk + 34);
    const auto *kk_35 = buffer.data(kk + 35);
    const auto *kk_36 = buffer.data(kk + 36);
    const auto *kk_37 = buffer.data(kk + 37);
    const auto *kk_38 = buffer.data(kk + 38);
    const auto *kk_39 = buffer.data(kk + 39);
    const auto *kk_40 = buffer.data(kk + 40);
    const auto *kk_41 = buffer.data(kk + 41);
    const auto *kk_42 = buffer.data(kk + 42);
    const auto *kk_43 = buffer.data(kk + 43);
    const auto *kk_44 = buffer.data(kk + 44);
    const auto *kk_45 = buffer.data(kk + 45);
    const auto *kk_46 = buffer.data(kk + 46);
    const auto *kk_47 = buffer.data(kk + 47);
    const auto *kk_48 = buffer.data(kk + 48);
    const auto *kk_49 = buffer.data(kk + 49);
    const auto *kk_50 = buffer.data(kk + 50);
    const auto *kk_51 = buffer.data(kk + 51);
    const auto *kk_52 = buffer.data(kk + 52);
    const auto *kk_53 = buffer.data(kk + 53);
    const auto *kk_54 = buffer.data(kk + 54);
    const auto *kk_55 = buffer.data(kk + 55);
    const auto *kk_56 = buffer.data(kk + 56);
    const auto *kk_57 = buffer.data(kk + 57);
    const auto *kk_58 = buffer.data(kk + 58);
    const auto *kk_59 = buffer.data(kk + 59);
    const auto *kk_60 = buffer.data(kk + 60);
    const auto *kk_61 = buffer.data(kk + 61);
    const auto *kk_62 = buffer.data(kk + 62);
    const auto *kk_63 = buffer.data(kk + 63);
    const auto *kk_64 = buffer.data(kk + 64);
    const auto *kk_65 = buffer.data(kk + 65);
    const auto *kk_66 = buffer.data(kk + 66);
    const auto *kk_67 = buffer.data(kk + 67);
    const auto *kk_68 = buffer.data(kk + 68);
    const auto *kk_69 = buffer.data(kk + 69);
    const auto *kk_70 = buffer.data(kk + 70);
    const auto *kk_71 = buffer.data(kk + 71);
    const auto *kk_72 = buffer.data(kk + 72);
    const auto *kk_73 = buffer.data(kk + 73);
    const auto *kk_74 = buffer.data(kk + 74);
    const auto *kk_75 = buffer.data(kk + 75);
    const auto *kk_76 = buffer.data(kk + 76);
    const auto *kk_77 = buffer.data(kk + 77);
    const auto *kk_78 = buffer.data(kk + 78);
    const auto *kk_79 = buffer.data(kk + 79);
    const auto *kk_80 = buffer.data(kk + 80);
    const auto *kk_81 = buffer.data(kk + 81);
    const auto *kk_82 = buffer.data(kk + 82);
    const auto *kk_83 = buffer.data(kk + 83);
    const auto *kk_84 = buffer.data(kk + 84);
    const auto *kk_85 = buffer.data(kk + 85);
    const auto *kk_86 = buffer.data(kk + 86);
    const auto *kk_87 = buffer.data(kk + 87);
    const auto *kk_88 = buffer.data(kk + 88);
    const auto *kk_89 = buffer.data(kk + 89);
    const auto *kk_90 = buffer.data(kk + 90);
    const auto *kk_91 = buffer.data(kk + 91);
    const auto *kk_92 = buffer.data(kk + 92);
    const auto *kk_93 = buffer.data(kk + 93);
    const auto *kk_94 = buffer.data(kk + 94);
    const auto *kk_95 = buffer.data(kk + 95);
    const auto *kk_96 = buffer.data(kk + 96);
    const auto *kk_97 = buffer.data(kk + 97);
    const auto *kk_98 = buffer.data(kk + 98);
    const auto *kk_99 = buffer.data(kk + 99);
    const auto *kk_100 = buffer.data(kk + 100);
    const auto *kk_101 = buffer.data(kk + 101);
    const auto *kk_102 = buffer.data(kk + 102);
    const auto *kk_103 = buffer.data(kk + 103);
    const auto *kk_104 = buffer.data(kk + 104);
    const auto *kk_105 = buffer.data(kk + 105);
    const auto *kk_106 = buffer.data(kk + 106);
    const auto *kk_107 = buffer.data(kk + 107);
    const auto *kk_108 = buffer.data(kk + 108);
    const auto *kk_109 = buffer.data(kk + 109);
    const auto *kk_110 = buffer.data(kk + 110);
    const auto *kk_111 = buffer.data(kk + 111);
    const auto *kk_112 = buffer.data(kk + 112);
    const auto *kk_113 = buffer.data(kk + 113);
    const auto *kk_114 = buffer.data(kk + 114);
    const auto *kk_115 = buffer.data(kk + 115);
    const auto *kk_116 = buffer.data(kk + 116);
    const auto *kk_117 = buffer.data(kk + 117);
    const auto *kk_118 = buffer.data(kk + 118);
    const auto *kk_119 = buffer.data(kk + 119);
    const auto *kk_120 = buffer.data(kk + 120);
    const auto *kk_121 = buffer.data(kk + 121);
    const auto *kk_122 = buffer.data(kk + 122);
    const auto *kk_123 = buffer.data(kk + 123);
    const auto *kk_124 = buffer.data(kk + 124);
    const auto *kk_125 = buffer.data(kk + 125);
    const auto *kk_126 = buffer.data(kk + 126);
    const auto *kk_127 = buffer.data(kk + 127);
    const auto *kk_128 = buffer.data(kk + 128);
    const auto *kk_129 = buffer.data(kk + 129);
    const auto *kk_130 = buffer.data(kk + 130);
    const auto *kk_131 = buffer.data(kk + 131);
    const auto *kk_132 = buffer.data(kk + 132);
    const auto *kk_133 = buffer.data(kk + 133);
    const auto *kk_134 = buffer.data(kk + 134);
    const auto *kk_135 = buffer.data(kk + 135);
    const auto *kk_136 = buffer.data(kk + 136);
    const auto *kk_137 = buffer.data(kk + 137);
    const auto *kk_138 = buffer.data(kk + 138);
    const auto *kk_139 = buffer.data(kk + 139);
    const auto *kk_140 = buffer.data(kk + 140);
    const auto *kk_141 = buffer.data(kk + 141);
    const auto *kk_142 = buffer.data(kk + 142);
    const auto *kk_143 = buffer.data(kk + 143);
    const auto *kk_144 = buffer.data(kk + 144);
    const auto *kk_145 = buffer.data(kk + 145);
    const auto *kk_146 = buffer.data(kk + 146);
    const auto *kk_147 = buffer.data(kk + 147);
    const auto *kk_148 = buffer.data(kk + 148);
    const auto *kk_149 = buffer.data(kk + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, hk_0, hk_1, hk_2, hk_3, hk_4, kk_0, kk_1, \
                         kk_2, kk_3, kk_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -6.0 * hk_0[k]
                 + f_0 * kk_0[k];

        t_1[k] = -6.0 * hk_1[k]
                 + f_0 * kk_1[k];

        t_2[k] = -6.0 * hk_2[k]
                 + f_0 * kk_2[k];

        t_3[k] = -6.0 * hk_3[k]
                 + f_0 * kk_3[k];

        t_4[k] = -6.0 * hk_4[k]
                 + f_0 * kk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, hk_5, hk_6, hk_7, hk_8, hk_9, kk_5, kk_6, \
                         kk_7, kk_8, kk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -6.0 * hk_5[k]
                 + f_0 * kk_5[k];

        t_6[k] = -6.0 * hk_6[k]
                 + f_0 * kk_6[k];

        t_7[k] = -6.0 * hk_7[k]
                 + f_0 * kk_7[k];

        t_8[k] = -6.0 * hk_8[k]
                 + f_0 * kk_8[k];

        t_9[k] = -6.0 * hk_9[k]
                 + f_0 * kk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, hk_10, hk_11, hk_12, hk_13, hk_14, \
                         kk_10, kk_11, kk_12, kk_13, kk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -6.0 * hk_10[k]
                  + f_0 * kk_10[k];

        t_11[k] = -6.0 * hk_11[k]
                  + f_0 * kk_11[k];

        t_12[k] = -6.0 * hk_12[k]
                  + f_0 * kk_12[k];

        t_13[k] = -6.0 * hk_13[k]
                  + f_0 * kk_13[k];

        t_14[k] = -6.0 * hk_14[k]
                  + f_0 * kk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, hk_15, hk_16, hk_17, hk_18, hk_19, \
                         kk_15, kk_16, kk_17, kk_18, kk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -6.0 * hk_15[k]
                  + f_0 * kk_15[k];

        t_16[k] = -6.0 * hk_16[k]
                  + f_0 * kk_16[k];

        t_17[k] = -6.0 * hk_17[k]
                  + f_0 * kk_17[k];

        t_18[k] = -6.0 * hk_18[k]
                  + f_0 * kk_18[k];

        t_19[k] = -6.0 * hk_19[k]
                  + f_0 * kk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, hk_20, hk_21, hk_22, hk_23, hk_24, \
                         kk_20, kk_21, kk_22, kk_23, kk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -6.0 * hk_20[k]
                  + f_0 * kk_20[k];

        t_21[k] = -6.0 * hk_21[k]
                  + f_0 * kk_21[k];

        t_22[k] = -6.0 * hk_22[k]
                  + f_0 * kk_22[k];

        t_23[k] = -6.0 * hk_23[k]
                  + f_0 * kk_23[k];

        t_24[k] = -6.0 * hk_24[k]
                  + f_0 * kk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, hk_25, hk_26, hk_27, hk_28, hk_29, \
                         kk_25, kk_26, kk_27, kk_28, kk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -6.0 * hk_25[k]
                  + f_0 * kk_25[k];

        t_26[k] = -6.0 * hk_26[k]
                  + f_0 * kk_26[k];

        t_27[k] = -6.0 * hk_27[k]
                  + f_0 * kk_27[k];

        t_28[k] = -6.0 * hk_28[k]
                  + f_0 * kk_28[k];

        t_29[k] = -6.0 * hk_29[k]
                  + f_0 * kk_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, hk_30, hk_31, hk_32, hk_33, hk_34, \
                         kk_30, kk_31, kk_32, kk_33, kk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -6.0 * hk_30[k]
                  + f_0 * kk_30[k];

        t_31[k] = -6.0 * hk_31[k]
                  + f_0 * kk_31[k];

        t_32[k] = -6.0 * hk_32[k]
                  + f_0 * kk_32[k];

        t_33[k] = -6.0 * hk_33[k]
                  + f_0 * kk_33[k];

        t_34[k] = -6.0 * hk_34[k]
                  + f_0 * kk_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, hk_35, hk_36, hk_37, hk_38, hk_39, \
                         kk_35, kk_36, kk_37, kk_38, kk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -6.0 * hk_35[k]
                  + f_0 * kk_35[k];

        t_36[k] = -5.0 * hk_36[k]
                  + f_0 * kk_36[k];

        t_37[k] = -5.0 * hk_37[k]
                  + f_0 * kk_37[k];

        t_38[k] = -5.0 * hk_38[k]
                  + f_0 * kk_38[k];

        t_39[k] = -5.0 * hk_39[k]
                  + f_0 * kk_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, hk_40, hk_41, hk_42, hk_43, hk_44, \
                         kk_40, kk_41, kk_42, kk_43, kk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -5.0 * hk_40[k]
                  + f_0 * kk_40[k];

        t_41[k] = -5.0 * hk_41[k]
                  + f_0 * kk_41[k];

        t_42[k] = -5.0 * hk_42[k]
                  + f_0 * kk_42[k];

        t_43[k] = -5.0 * hk_43[k]
                  + f_0 * kk_43[k];

        t_44[k] = -5.0 * hk_44[k]
                  + f_0 * kk_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hk_45, hk_46, hk_47, hk_48, hk_49, \
                         kk_45, kk_46, kk_47, kk_48, kk_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -5.0 * hk_45[k]
                  + f_0 * kk_45[k];

        t_46[k] = -5.0 * hk_46[k]
                  + f_0 * kk_46[k];

        t_47[k] = -5.0 * hk_47[k]
                  + f_0 * kk_47[k];

        t_48[k] = -5.0 * hk_48[k]
                  + f_0 * kk_48[k];

        t_49[k] = -5.0 * hk_49[k]
                  + f_0 * kk_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hk_50, hk_51, hk_52, hk_53, hk_54, \
                         kk_50, kk_51, kk_52, kk_53, kk_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -5.0 * hk_50[k]
                  + f_0 * kk_50[k];

        t_51[k] = -5.0 * hk_51[k]
                  + f_0 * kk_51[k];

        t_52[k] = -5.0 * hk_52[k]
                  + f_0 * kk_52[k];

        t_53[k] = -5.0 * hk_53[k]
                  + f_0 * kk_53[k];

        t_54[k] = -5.0 * hk_54[k]
                  + f_0 * kk_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, hk_55, hk_56, hk_57, hk_58, hk_59, \
                         kk_55, kk_56, kk_57, kk_58, kk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -5.0 * hk_55[k]
                  + f_0 * kk_55[k];

        t_56[k] = -5.0 * hk_56[k]
                  + f_0 * kk_56[k];

        t_57[k] = -5.0 * hk_57[k]
                  + f_0 * kk_57[k];

        t_58[k] = -5.0 * hk_58[k]
                  + f_0 * kk_58[k];

        t_59[k] = -5.0 * hk_59[k]
                  + f_0 * kk_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, hk_60, hk_61, hk_62, hk_63, hk_64, \
                         kk_60, kk_61, kk_62, kk_63, kk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -5.0 * hk_60[k]
                  + f_0 * kk_60[k];

        t_61[k] = -5.0 * hk_61[k]
                  + f_0 * kk_61[k];

        t_62[k] = -5.0 * hk_62[k]
                  + f_0 * kk_62[k];

        t_63[k] = -5.0 * hk_63[k]
                  + f_0 * kk_63[k];

        t_64[k] = -5.0 * hk_64[k]
                  + f_0 * kk_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, hk_65, hk_66, hk_67, hk_68, hk_69, \
                         kk_65, kk_66, kk_67, kk_68, kk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -5.0 * hk_65[k]
                  + f_0 * kk_65[k];

        t_66[k] = -5.0 * hk_66[k]
                  + f_0 * kk_66[k];

        t_67[k] = -5.0 * hk_67[k]
                  + f_0 * kk_67[k];

        t_68[k] = -5.0 * hk_68[k]
                  + f_0 * kk_68[k];

        t_69[k] = -5.0 * hk_69[k]
                  + f_0 * kk_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, hk_70, hk_71, hk_72, hk_73, hk_74, \
                         kk_70, kk_71, kk_72, kk_73, kk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -5.0 * hk_70[k]
                  + f_0 * kk_70[k];

        t_71[k] = -5.0 * hk_71[k]
                  + f_0 * kk_71[k];

        t_72[k] = -5.0 * hk_72[k]
                  + f_0 * kk_72[k];

        t_73[k] = -5.0 * hk_73[k]
                  + f_0 * kk_73[k];

        t_74[k] = -5.0 * hk_74[k]
                  + f_0 * kk_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, hk_75, hk_76, hk_77, hk_78, hk_79, \
                         kk_75, kk_76, kk_77, kk_78, kk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -5.0 * hk_75[k]
                  + f_0 * kk_75[k];

        t_76[k] = -5.0 * hk_76[k]
                  + f_0 * kk_76[k];

        t_77[k] = -5.0 * hk_77[k]
                  + f_0 * kk_77[k];

        t_78[k] = -5.0 * hk_78[k]
                  + f_0 * kk_78[k];

        t_79[k] = -5.0 * hk_79[k]
                  + f_0 * kk_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, hk_80, hk_81, hk_82, hk_83, hk_84, \
                         kk_80, kk_81, kk_82, kk_83, kk_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -5.0 * hk_80[k]
                  + f_0 * kk_80[k];

        t_81[k] = -5.0 * hk_81[k]
                  + f_0 * kk_81[k];

        t_82[k] = -5.0 * hk_82[k]
                  + f_0 * kk_82[k];

        t_83[k] = -5.0 * hk_83[k]
                  + f_0 * kk_83[k];

        t_84[k] = -5.0 * hk_84[k]
                  + f_0 * kk_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, hk_85, hk_86, hk_87, hk_88, hk_89, \
                         kk_85, kk_86, kk_87, kk_88, kk_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -5.0 * hk_85[k]
                  + f_0 * kk_85[k];

        t_86[k] = -5.0 * hk_86[k]
                  + f_0 * kk_86[k];

        t_87[k] = -5.0 * hk_87[k]
                  + f_0 * kk_87[k];

        t_88[k] = -5.0 * hk_88[k]
                  + f_0 * kk_88[k];

        t_89[k] = -5.0 * hk_89[k]
                  + f_0 * kk_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, hk_90, hk_91, hk_92, hk_93, hk_94, \
                         kk_90, kk_91, kk_92, kk_93, kk_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -5.0 * hk_90[k]
                  + f_0 * kk_90[k];

        t_91[k] = -5.0 * hk_91[k]
                  + f_0 * kk_91[k];

        t_92[k] = -5.0 * hk_92[k]
                  + f_0 * kk_92[k];

        t_93[k] = -5.0 * hk_93[k]
                  + f_0 * kk_93[k];

        t_94[k] = -5.0 * hk_94[k]
                  + f_0 * kk_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, hk_95, hk_96, hk_97, hk_98, hk_99, \
                         kk_95, kk_96, kk_97, kk_98, kk_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -5.0 * hk_95[k]
                  + f_0 * kk_95[k];

        t_96[k] = -5.0 * hk_96[k]
                  + f_0 * kk_96[k];

        t_97[k] = -5.0 * hk_97[k]
                  + f_0 * kk_97[k];

        t_98[k] = -5.0 * hk_98[k]
                  + f_0 * kk_98[k];

        t_99[k] = -5.0 * hk_99[k]
                  + f_0 * kk_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, hk_100, hk_101, hk_102, hk_103, \
                         hk_104, kk_100, kk_101, kk_102, kk_103, \
                         kk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -5.0 * hk_100[k]
                   + f_0 * kk_100[k];

        t_101[k] = -5.0 * hk_101[k]
                   + f_0 * kk_101[k];

        t_102[k] = -5.0 * hk_102[k]
                   + f_0 * kk_102[k];

        t_103[k] = -5.0 * hk_103[k]
                   + f_0 * kk_103[k];

        t_104[k] = -5.0 * hk_104[k]
                   + f_0 * kk_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, hk_105, hk_106, hk_107, hk_108, \
                         hk_109, kk_105, kk_106, kk_107, kk_108, \
                         kk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -5.0 * hk_105[k]
                   + f_0 * kk_105[k];

        t_106[k] = -5.0 * hk_106[k]
                   + f_0 * kk_106[k];

        t_107[k] = -5.0 * hk_107[k]
                   + f_0 * kk_107[k];

        t_108[k] = -4.0 * hk_108[k]
                   + f_0 * kk_108[k];

        t_109[k] = -4.0 * hk_109[k]
                   + f_0 * kk_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, hk_110, hk_111, hk_112, hk_113, \
                         hk_114, kk_110, kk_111, kk_112, kk_113, \
                         kk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -4.0 * hk_110[k]
                   + f_0 * kk_110[k];

        t_111[k] = -4.0 * hk_111[k]
                   + f_0 * kk_111[k];

        t_112[k] = -4.0 * hk_112[k]
                   + f_0 * kk_112[k];

        t_113[k] = -4.0 * hk_113[k]
                   + f_0 * kk_113[k];

        t_114[k] = -4.0 * hk_114[k]
                   + f_0 * kk_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, hk_115, hk_116, hk_117, hk_118, \
                         hk_119, kk_115, kk_116, kk_117, kk_118, \
                         kk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -4.0 * hk_115[k]
                   + f_0 * kk_115[k];

        t_116[k] = -4.0 * hk_116[k]
                   + f_0 * kk_116[k];

        t_117[k] = -4.0 * hk_117[k]
                   + f_0 * kk_117[k];

        t_118[k] = -4.0 * hk_118[k]
                   + f_0 * kk_118[k];

        t_119[k] = -4.0 * hk_119[k]
                   + f_0 * kk_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, hk_120, hk_121, hk_122, hk_123, \
                         hk_124, kk_120, kk_121, kk_122, kk_123, \
                         kk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -4.0 * hk_120[k]
                   + f_0 * kk_120[k];

        t_121[k] = -4.0 * hk_121[k]
                   + f_0 * kk_121[k];

        t_122[k] = -4.0 * hk_122[k]
                   + f_0 * kk_122[k];

        t_123[k] = -4.0 * hk_123[k]
                   + f_0 * kk_123[k];

        t_124[k] = -4.0 * hk_124[k]
                   + f_0 * kk_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, hk_125, hk_126, hk_127, hk_128, \
                         hk_129, kk_125, kk_126, kk_127, kk_128, \
                         kk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -4.0 * hk_125[k]
                   + f_0 * kk_125[k];

        t_126[k] = -4.0 * hk_126[k]
                   + f_0 * kk_126[k];

        t_127[k] = -4.0 * hk_127[k]
                   + f_0 * kk_127[k];

        t_128[k] = -4.0 * hk_128[k]
                   + f_0 * kk_128[k];

        t_129[k] = -4.0 * hk_129[k]
                   + f_0 * kk_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, hk_130, hk_131, hk_132, hk_133, \
                         hk_134, kk_130, kk_131, kk_132, kk_133, \
                         kk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -4.0 * hk_130[k]
                   + f_0 * kk_130[k];

        t_131[k] = -4.0 * hk_131[k]
                   + f_0 * kk_131[k];

        t_132[k] = -4.0 * hk_132[k]
                   + f_0 * kk_132[k];

        t_133[k] = -4.0 * hk_133[k]
                   + f_0 * kk_133[k];

        t_134[k] = -4.0 * hk_134[k]
                   + f_0 * kk_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, hk_135, hk_136, hk_137, hk_138, \
                         hk_139, kk_135, kk_136, kk_137, kk_138, \
                         kk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -4.0 * hk_135[k]
                   + f_0 * kk_135[k];

        t_136[k] = -4.0 * hk_136[k]
                   + f_0 * kk_136[k];

        t_137[k] = -4.0 * hk_137[k]
                   + f_0 * kk_137[k];

        t_138[k] = -4.0 * hk_138[k]
                   + f_0 * kk_138[k];

        t_139[k] = -4.0 * hk_139[k]
                   + f_0 * kk_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, hk_140, hk_141, hk_142, hk_143, \
                         hk_144, kk_140, kk_141, kk_142, kk_143, \
                         kk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -4.0 * hk_140[k]
                   + f_0 * kk_140[k];

        t_141[k] = -4.0 * hk_141[k]
                   + f_0 * kk_141[k];

        t_142[k] = -4.0 * hk_142[k]
                   + f_0 * kk_142[k];

        t_143[k] = -4.0 * hk_143[k]
                   + f_0 * kk_143[k];

        t_144[k] = -4.0 * hk_144[k]
                   + f_0 * kk_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, hk_145, hk_146, hk_147, hk_148, \
                         hk_149, kk_145, kk_146, kk_147, kk_148, \
                         kk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -4.0 * hk_145[k]
                   + f_0 * kk_145[k];

        t_146[k] = -4.0 * hk_146[k]
                   + f_0 * kk_146[k];

        t_147[k] = -4.0 * hk_147[k]
                   + f_0 * kk_147[k];

        t_148[k] = -4.0 * hk_148[k]
                   + f_0 * kk_148[k];

        t_149[k] = -4.0 * hk_149[k]
                   + f_0 * kk_149[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_150 = buffer.data(hk + 150);
    const auto *hk_151 = buffer.data(hk + 151);
    const auto *hk_152 = buffer.data(hk + 152);
    const auto *hk_153 = buffer.data(hk + 153);
    const auto *hk_154 = buffer.data(hk + 154);
    const auto *hk_155 = buffer.data(hk + 155);
    const auto *hk_156 = buffer.data(hk + 156);
    const auto *hk_157 = buffer.data(hk + 157);
    const auto *hk_158 = buffer.data(hk + 158);
    const auto *hk_159 = buffer.data(hk + 159);
    const auto *hk_160 = buffer.data(hk + 160);
    const auto *hk_161 = buffer.data(hk + 161);
    const auto *hk_162 = buffer.data(hk + 162);
    const auto *hk_163 = buffer.data(hk + 163);
    const auto *hk_164 = buffer.data(hk + 164);
    const auto *hk_165 = buffer.data(hk + 165);
    const auto *hk_166 = buffer.data(hk + 166);
    const auto *hk_167 = buffer.data(hk + 167);
    const auto *hk_168 = buffer.data(hk + 168);
    const auto *hk_169 = buffer.data(hk + 169);
    const auto *hk_170 = buffer.data(hk + 170);
    const auto *hk_171 = buffer.data(hk + 171);
    const auto *hk_172 = buffer.data(hk + 172);
    const auto *hk_173 = buffer.data(hk + 173);
    const auto *hk_174 = buffer.data(hk + 174);
    const auto *hk_175 = buffer.data(hk + 175);
    const auto *hk_176 = buffer.data(hk + 176);
    const auto *hk_177 = buffer.data(hk + 177);
    const auto *hk_178 = buffer.data(hk + 178);
    const auto *hk_179 = buffer.data(hk + 179);
    const auto *hk_180 = buffer.data(hk + 180);
    const auto *hk_181 = buffer.data(hk + 181);
    const auto *hk_182 = buffer.data(hk + 182);
    const auto *hk_183 = buffer.data(hk + 183);
    const auto *hk_184 = buffer.data(hk + 184);
    const auto *hk_185 = buffer.data(hk + 185);
    const auto *hk_186 = buffer.data(hk + 186);
    const auto *hk_187 = buffer.data(hk + 187);
    const auto *hk_188 = buffer.data(hk + 188);
    const auto *hk_189 = buffer.data(hk + 189);
    const auto *hk_190 = buffer.data(hk + 190);
    const auto *hk_191 = buffer.data(hk + 191);
    const auto *hk_192 = buffer.data(hk + 192);
    const auto *hk_193 = buffer.data(hk + 193);
    const auto *hk_194 = buffer.data(hk + 194);
    const auto *hk_195 = buffer.data(hk + 195);
    const auto *hk_196 = buffer.data(hk + 196);
    const auto *hk_197 = buffer.data(hk + 197);
    const auto *hk_198 = buffer.data(hk + 198);
    const auto *hk_199 = buffer.data(hk + 199);
    const auto *hk_200 = buffer.data(hk + 200);
    const auto *hk_201 = buffer.data(hk + 201);
    const auto *hk_202 = buffer.data(hk + 202);
    const auto *hk_203 = buffer.data(hk + 203);
    const auto *hk_204 = buffer.data(hk + 204);
    const auto *hk_205 = buffer.data(hk + 205);
    const auto *hk_206 = buffer.data(hk + 206);
    const auto *hk_207 = buffer.data(hk + 207);
    const auto *hk_208 = buffer.data(hk + 208);
    const auto *hk_209 = buffer.data(hk + 209);
    const auto *hk_210 = buffer.data(hk + 210);
    const auto *hk_211 = buffer.data(hk + 211);
    const auto *hk_212 = buffer.data(hk + 212);
    const auto *hk_213 = buffer.data(hk + 213);
    const auto *hk_214 = buffer.data(hk + 214);
    const auto *hk_215 = buffer.data(hk + 215);
    const auto *hk_216 = buffer.data(hk + 216);
    const auto *hk_217 = buffer.data(hk + 217);
    const auto *hk_218 = buffer.data(hk + 218);
    const auto *hk_219 = buffer.data(hk + 219);
    const auto *hk_220 = buffer.data(hk + 220);
    const auto *hk_221 = buffer.data(hk + 221);
    const auto *hk_222 = buffer.data(hk + 222);
    const auto *hk_223 = buffer.data(hk + 223);
    const auto *hk_224 = buffer.data(hk + 224);
    const auto *hk_225 = buffer.data(hk + 225);
    const auto *hk_226 = buffer.data(hk + 226);
    const auto *hk_227 = buffer.data(hk + 227);
    const auto *hk_228 = buffer.data(hk + 228);
    const auto *hk_229 = buffer.data(hk + 229);
    const auto *hk_230 = buffer.data(hk + 230);
    const auto *hk_231 = buffer.data(hk + 231);
    const auto *hk_232 = buffer.data(hk + 232);
    const auto *hk_233 = buffer.data(hk + 233);
    const auto *hk_234 = buffer.data(hk + 234);
    const auto *hk_235 = buffer.data(hk + 235);
    const auto *hk_236 = buffer.data(hk + 236);
    const auto *hk_237 = buffer.data(hk + 237);
    const auto *hk_238 = buffer.data(hk + 238);
    const auto *hk_239 = buffer.data(hk + 239);
    const auto *hk_240 = buffer.data(hk + 240);
    const auto *hk_241 = buffer.data(hk + 241);
    const auto *hk_242 = buffer.data(hk + 242);
    const auto *hk_243 = buffer.data(hk + 243);
    const auto *hk_244 = buffer.data(hk + 244);
    const auto *hk_245 = buffer.data(hk + 245);
    const auto *hk_246 = buffer.data(hk + 246);
    const auto *hk_247 = buffer.data(hk + 247);
    const auto *hk_248 = buffer.data(hk + 248);
    const auto *hk_249 = buffer.data(hk + 249);
    const auto *hk_250 = buffer.data(hk + 250);
    const auto *hk_251 = buffer.data(hk + 251);
    const auto *hk_252 = buffer.data(hk + 252);
    const auto *hk_253 = buffer.data(hk + 253);
    const auto *hk_254 = buffer.data(hk + 254);
    const auto *hk_255 = buffer.data(hk + 255);
    const auto *hk_256 = buffer.data(hk + 256);
    const auto *hk_257 = buffer.data(hk + 257);
    const auto *hk_258 = buffer.data(hk + 258);
    const auto *hk_259 = buffer.data(hk + 259);
    const auto *hk_260 = buffer.data(hk + 260);
    const auto *hk_261 = buffer.data(hk + 261);
    const auto *hk_262 = buffer.data(hk + 262);
    const auto *hk_263 = buffer.data(hk + 263);
    const auto *hk_264 = buffer.data(hk + 264);
    const auto *hk_265 = buffer.data(hk + 265);
    const auto *hk_266 = buffer.data(hk + 266);
    const auto *hk_267 = buffer.data(hk + 267);
    const auto *hk_268 = buffer.data(hk + 268);
    const auto *hk_269 = buffer.data(hk + 269);
    const auto *hk_270 = buffer.data(hk + 270);
    const auto *hk_271 = buffer.data(hk + 271);
    const auto *hk_272 = buffer.data(hk + 272);
    const auto *hk_273 = buffer.data(hk + 273);
    const auto *hk_274 = buffer.data(hk + 274);
    const auto *hk_275 = buffer.data(hk + 275);
    const auto *hk_276 = buffer.data(hk + 276);
    const auto *hk_277 = buffer.data(hk + 277);
    const auto *hk_278 = buffer.data(hk + 278);
    const auto *hk_279 = buffer.data(hk + 279);
    const auto *hk_280 = buffer.data(hk + 280);
    const auto *hk_281 = buffer.data(hk + 281);
    const auto *hk_282 = buffer.data(hk + 282);
    const auto *hk_283 = buffer.data(hk + 283);
    const auto *hk_284 = buffer.data(hk + 284);
    const auto *hk_285 = buffer.data(hk + 285);
    const auto *hk_286 = buffer.data(hk + 286);
    const auto *hk_287 = buffer.data(hk + 287);
    const auto *hk_288 = buffer.data(hk + 288);
    const auto *hk_289 = buffer.data(hk + 289);
    const auto *hk_290 = buffer.data(hk + 290);
    const auto *hk_291 = buffer.data(hk + 291);
    const auto *hk_292 = buffer.data(hk + 292);
    const auto *hk_293 = buffer.data(hk + 293);
    const auto *hk_294 = buffer.data(hk + 294);
    const auto *hk_295 = buffer.data(hk + 295);
    const auto *hk_296 = buffer.data(hk + 296);
    const auto *hk_297 = buffer.data(hk + 297);
    const auto *hk_298 = buffer.data(hk + 298);
    const auto *hk_299 = buffer.data(hk + 299);

    const auto *kk_150 = buffer.data(kk + 150);
    const auto *kk_151 = buffer.data(kk + 151);
    const auto *kk_152 = buffer.data(kk + 152);
    const auto *kk_153 = buffer.data(kk + 153);
    const auto *kk_154 = buffer.data(kk + 154);
    const auto *kk_155 = buffer.data(kk + 155);
    const auto *kk_156 = buffer.data(kk + 156);
    const auto *kk_157 = buffer.data(kk + 157);
    const auto *kk_158 = buffer.data(kk + 158);
    const auto *kk_159 = buffer.data(kk + 159);
    const auto *kk_160 = buffer.data(kk + 160);
    const auto *kk_161 = buffer.data(kk + 161);
    const auto *kk_162 = buffer.data(kk + 162);
    const auto *kk_163 = buffer.data(kk + 163);
    const auto *kk_164 = buffer.data(kk + 164);
    const auto *kk_165 = buffer.data(kk + 165);
    const auto *kk_166 = buffer.data(kk + 166);
    const auto *kk_167 = buffer.data(kk + 167);
    const auto *kk_168 = buffer.data(kk + 168);
    const auto *kk_169 = buffer.data(kk + 169);
    const auto *kk_170 = buffer.data(kk + 170);
    const auto *kk_171 = buffer.data(kk + 171);
    const auto *kk_172 = buffer.data(kk + 172);
    const auto *kk_173 = buffer.data(kk + 173);
    const auto *kk_174 = buffer.data(kk + 174);
    const auto *kk_175 = buffer.data(kk + 175);
    const auto *kk_176 = buffer.data(kk + 176);
    const auto *kk_177 = buffer.data(kk + 177);
    const auto *kk_178 = buffer.data(kk + 178);
    const auto *kk_179 = buffer.data(kk + 179);
    const auto *kk_180 = buffer.data(kk + 180);
    const auto *kk_181 = buffer.data(kk + 181);
    const auto *kk_182 = buffer.data(kk + 182);
    const auto *kk_183 = buffer.data(kk + 183);
    const auto *kk_184 = buffer.data(kk + 184);
    const auto *kk_185 = buffer.data(kk + 185);
    const auto *kk_186 = buffer.data(kk + 186);
    const auto *kk_187 = buffer.data(kk + 187);
    const auto *kk_188 = buffer.data(kk + 188);
    const auto *kk_189 = buffer.data(kk + 189);
    const auto *kk_190 = buffer.data(kk + 190);
    const auto *kk_191 = buffer.data(kk + 191);
    const auto *kk_192 = buffer.data(kk + 192);
    const auto *kk_193 = buffer.data(kk + 193);
    const auto *kk_194 = buffer.data(kk + 194);
    const auto *kk_195 = buffer.data(kk + 195);
    const auto *kk_196 = buffer.data(kk + 196);
    const auto *kk_197 = buffer.data(kk + 197);
    const auto *kk_198 = buffer.data(kk + 198);
    const auto *kk_199 = buffer.data(kk + 199);
    const auto *kk_200 = buffer.data(kk + 200);
    const auto *kk_201 = buffer.data(kk + 201);
    const auto *kk_202 = buffer.data(kk + 202);
    const auto *kk_203 = buffer.data(kk + 203);
    const auto *kk_204 = buffer.data(kk + 204);
    const auto *kk_205 = buffer.data(kk + 205);
    const auto *kk_206 = buffer.data(kk + 206);
    const auto *kk_207 = buffer.data(kk + 207);
    const auto *kk_208 = buffer.data(kk + 208);
    const auto *kk_209 = buffer.data(kk + 209);
    const auto *kk_210 = buffer.data(kk + 210);
    const auto *kk_211 = buffer.data(kk + 211);
    const auto *kk_212 = buffer.data(kk + 212);
    const auto *kk_213 = buffer.data(kk + 213);
    const auto *kk_214 = buffer.data(kk + 214);
    const auto *kk_215 = buffer.data(kk + 215);
    const auto *kk_216 = buffer.data(kk + 216);
    const auto *kk_217 = buffer.data(kk + 217);
    const auto *kk_218 = buffer.data(kk + 218);
    const auto *kk_219 = buffer.data(kk + 219);
    const auto *kk_220 = buffer.data(kk + 220);
    const auto *kk_221 = buffer.data(kk + 221);
    const auto *kk_222 = buffer.data(kk + 222);
    const auto *kk_223 = buffer.data(kk + 223);
    const auto *kk_224 = buffer.data(kk + 224);
    const auto *kk_225 = buffer.data(kk + 225);
    const auto *kk_226 = buffer.data(kk + 226);
    const auto *kk_227 = buffer.data(kk + 227);
    const auto *kk_228 = buffer.data(kk + 228);
    const auto *kk_229 = buffer.data(kk + 229);
    const auto *kk_230 = buffer.data(kk + 230);
    const auto *kk_231 = buffer.data(kk + 231);
    const auto *kk_232 = buffer.data(kk + 232);
    const auto *kk_233 = buffer.data(kk + 233);
    const auto *kk_234 = buffer.data(kk + 234);
    const auto *kk_235 = buffer.data(kk + 235);
    const auto *kk_236 = buffer.data(kk + 236);
    const auto *kk_237 = buffer.data(kk + 237);
    const auto *kk_238 = buffer.data(kk + 238);
    const auto *kk_239 = buffer.data(kk + 239);
    const auto *kk_240 = buffer.data(kk + 240);
    const auto *kk_241 = buffer.data(kk + 241);
    const auto *kk_242 = buffer.data(kk + 242);
    const auto *kk_243 = buffer.data(kk + 243);
    const auto *kk_244 = buffer.data(kk + 244);
    const auto *kk_245 = buffer.data(kk + 245);
    const auto *kk_246 = buffer.data(kk + 246);
    const auto *kk_247 = buffer.data(kk + 247);
    const auto *kk_248 = buffer.data(kk + 248);
    const auto *kk_249 = buffer.data(kk + 249);
    const auto *kk_250 = buffer.data(kk + 250);
    const auto *kk_251 = buffer.data(kk + 251);
    const auto *kk_252 = buffer.data(kk + 252);
    const auto *kk_253 = buffer.data(kk + 253);
    const auto *kk_254 = buffer.data(kk + 254);
    const auto *kk_255 = buffer.data(kk + 255);
    const auto *kk_256 = buffer.data(kk + 256);
    const auto *kk_257 = buffer.data(kk + 257);
    const auto *kk_258 = buffer.data(kk + 258);
    const auto *kk_259 = buffer.data(kk + 259);
    const auto *kk_260 = buffer.data(kk + 260);
    const auto *kk_261 = buffer.data(kk + 261);
    const auto *kk_262 = buffer.data(kk + 262);
    const auto *kk_263 = buffer.data(kk + 263);
    const auto *kk_264 = buffer.data(kk + 264);
    const auto *kk_265 = buffer.data(kk + 265);
    const auto *kk_266 = buffer.data(kk + 266);
    const auto *kk_267 = buffer.data(kk + 267);
    const auto *kk_268 = buffer.data(kk + 268);
    const auto *kk_269 = buffer.data(kk + 269);
    const auto *kk_270 = buffer.data(kk + 270);
    const auto *kk_271 = buffer.data(kk + 271);
    const auto *kk_272 = buffer.data(kk + 272);
    const auto *kk_273 = buffer.data(kk + 273);
    const auto *kk_274 = buffer.data(kk + 274);
    const auto *kk_275 = buffer.data(kk + 275);
    const auto *kk_276 = buffer.data(kk + 276);
    const auto *kk_277 = buffer.data(kk + 277);
    const auto *kk_278 = buffer.data(kk + 278);
    const auto *kk_279 = buffer.data(kk + 279);
    const auto *kk_280 = buffer.data(kk + 280);
    const auto *kk_281 = buffer.data(kk + 281);
    const auto *kk_282 = buffer.data(kk + 282);
    const auto *kk_283 = buffer.data(kk + 283);
    const auto *kk_284 = buffer.data(kk + 284);
    const auto *kk_285 = buffer.data(kk + 285);
    const auto *kk_286 = buffer.data(kk + 286);
    const auto *kk_287 = buffer.data(kk + 287);
    const auto *kk_288 = buffer.data(kk + 288);
    const auto *kk_289 = buffer.data(kk + 289);
    const auto *kk_290 = buffer.data(kk + 290);
    const auto *kk_291 = buffer.data(kk + 291);
    const auto *kk_292 = buffer.data(kk + 292);
    const auto *kk_293 = buffer.data(kk + 293);
    const auto *kk_294 = buffer.data(kk + 294);
    const auto *kk_295 = buffer.data(kk + 295);
    const auto *kk_296 = buffer.data(kk + 296);
    const auto *kk_297 = buffer.data(kk + 297);
    const auto *kk_298 = buffer.data(kk + 298);
    const auto *kk_299 = buffer.data(kk + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, hk_150, hk_151, hk_152, hk_153, \
                         hk_154, kk_150, kk_151, kk_152, kk_153, \
                         kk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -4.0 * hk_150[k]
                   + f_0 * kk_150[k];

        t_151[k] = -4.0 * hk_151[k]
                   + f_0 * kk_151[k];

        t_152[k] = -4.0 * hk_152[k]
                   + f_0 * kk_152[k];

        t_153[k] = -4.0 * hk_153[k]
                   + f_0 * kk_153[k];

        t_154[k] = -4.0 * hk_154[k]
                   + f_0 * kk_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, hk_155, hk_156, hk_157, hk_158, \
                         hk_159, kk_155, kk_156, kk_157, kk_158, \
                         kk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -4.0 * hk_155[k]
                   + f_0 * kk_155[k];

        t_156[k] = -4.0 * hk_156[k]
                   + f_0 * kk_156[k];

        t_157[k] = -4.0 * hk_157[k]
                   + f_0 * kk_157[k];

        t_158[k] = -4.0 * hk_158[k]
                   + f_0 * kk_158[k];

        t_159[k] = -4.0 * hk_159[k]
                   + f_0 * kk_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, hk_160, hk_161, hk_162, hk_163, \
                         hk_164, kk_160, kk_161, kk_162, kk_163, \
                         kk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -4.0 * hk_160[k]
                   + f_0 * kk_160[k];

        t_161[k] = -4.0 * hk_161[k]
                   + f_0 * kk_161[k];

        t_162[k] = -4.0 * hk_162[k]
                   + f_0 * kk_162[k];

        t_163[k] = -4.0 * hk_163[k]
                   + f_0 * kk_163[k];

        t_164[k] = -4.0 * hk_164[k]
                   + f_0 * kk_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, hk_165, hk_166, hk_167, hk_168, \
                         hk_169, kk_165, kk_166, kk_167, kk_168, \
                         kk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -4.0 * hk_165[k]
                   + f_0 * kk_165[k];

        t_166[k] = -4.0 * hk_166[k]
                   + f_0 * kk_166[k];

        t_167[k] = -4.0 * hk_167[k]
                   + f_0 * kk_167[k];

        t_168[k] = -4.0 * hk_168[k]
                   + f_0 * kk_168[k];

        t_169[k] = -4.0 * hk_169[k]
                   + f_0 * kk_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, hk_170, hk_171, hk_172, hk_173, \
                         hk_174, kk_170, kk_171, kk_172, kk_173, \
                         kk_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -4.0 * hk_170[k]
                   + f_0 * kk_170[k];

        t_171[k] = -4.0 * hk_171[k]
                   + f_0 * kk_171[k];

        t_172[k] = -4.0 * hk_172[k]
                   + f_0 * kk_172[k];

        t_173[k] = -4.0 * hk_173[k]
                   + f_0 * kk_173[k];

        t_174[k] = -4.0 * hk_174[k]
                   + f_0 * kk_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, hk_175, hk_176, hk_177, hk_178, \
                         hk_179, kk_175, kk_176, kk_177, kk_178, \
                         kk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -4.0 * hk_175[k]
                   + f_0 * kk_175[k];

        t_176[k] = -4.0 * hk_176[k]
                   + f_0 * kk_176[k];

        t_177[k] = -4.0 * hk_177[k]
                   + f_0 * kk_177[k];

        t_178[k] = -4.0 * hk_178[k]
                   + f_0 * kk_178[k];

        t_179[k] = -4.0 * hk_179[k]
                   + f_0 * kk_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, hk_180, hk_181, hk_182, hk_183, \
                         hk_184, kk_180, kk_181, kk_182, kk_183, \
                         kk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -4.0 * hk_180[k]
                   + f_0 * kk_180[k];

        t_181[k] = -4.0 * hk_181[k]
                   + f_0 * kk_181[k];

        t_182[k] = -4.0 * hk_182[k]
                   + f_0 * kk_182[k];

        t_183[k] = -4.0 * hk_183[k]
                   + f_0 * kk_183[k];

        t_184[k] = -4.0 * hk_184[k]
                   + f_0 * kk_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, hk_185, hk_186, hk_187, hk_188, \
                         hk_189, kk_185, kk_186, kk_187, kk_188, \
                         kk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -4.0 * hk_185[k]
                   + f_0 * kk_185[k];

        t_186[k] = -4.0 * hk_186[k]
                   + f_0 * kk_186[k];

        t_187[k] = -4.0 * hk_187[k]
                   + f_0 * kk_187[k];

        t_188[k] = -4.0 * hk_188[k]
                   + f_0 * kk_188[k];

        t_189[k] = -4.0 * hk_189[k]
                   + f_0 * kk_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, hk_190, hk_191, hk_192, hk_193, \
                         hk_194, kk_190, kk_191, kk_192, kk_193, \
                         kk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -4.0 * hk_190[k]
                   + f_0 * kk_190[k];

        t_191[k] = -4.0 * hk_191[k]
                   + f_0 * kk_191[k];

        t_192[k] = -4.0 * hk_192[k]
                   + f_0 * kk_192[k];

        t_193[k] = -4.0 * hk_193[k]
                   + f_0 * kk_193[k];

        t_194[k] = -4.0 * hk_194[k]
                   + f_0 * kk_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, hk_195, hk_196, hk_197, hk_198, \
                         hk_199, kk_195, kk_196, kk_197, kk_198, \
                         kk_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -4.0 * hk_195[k]
                   + f_0 * kk_195[k];

        t_196[k] = -4.0 * hk_196[k]
                   + f_0 * kk_196[k];

        t_197[k] = -4.0 * hk_197[k]
                   + f_0 * kk_197[k];

        t_198[k] = -4.0 * hk_198[k]
                   + f_0 * kk_198[k];

        t_199[k] = -4.0 * hk_199[k]
                   + f_0 * kk_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, hk_200, hk_201, hk_202, hk_203, \
                         hk_204, kk_200, kk_201, kk_202, kk_203, \
                         kk_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -4.0 * hk_200[k]
                   + f_0 * kk_200[k];

        t_201[k] = -4.0 * hk_201[k]
                   + f_0 * kk_201[k];

        t_202[k] = -4.0 * hk_202[k]
                   + f_0 * kk_202[k];

        t_203[k] = -4.0 * hk_203[k]
                   + f_0 * kk_203[k];

        t_204[k] = -4.0 * hk_204[k]
                   + f_0 * kk_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, hk_205, hk_206, hk_207, hk_208, \
                         hk_209, kk_205, kk_206, kk_207, kk_208, \
                         kk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -4.0 * hk_205[k]
                   + f_0 * kk_205[k];

        t_206[k] = -4.0 * hk_206[k]
                   + f_0 * kk_206[k];

        t_207[k] = -4.0 * hk_207[k]
                   + f_0 * kk_207[k];

        t_208[k] = -4.0 * hk_208[k]
                   + f_0 * kk_208[k];

        t_209[k] = -4.0 * hk_209[k]
                   + f_0 * kk_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, hk_210, hk_211, hk_212, hk_213, \
                         hk_214, kk_210, kk_211, kk_212, kk_213, \
                         kk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -4.0 * hk_210[k]
                   + f_0 * kk_210[k];

        t_211[k] = -4.0 * hk_211[k]
                   + f_0 * kk_211[k];

        t_212[k] = -4.0 * hk_212[k]
                   + f_0 * kk_212[k];

        t_213[k] = -4.0 * hk_213[k]
                   + f_0 * kk_213[k];

        t_214[k] = -4.0 * hk_214[k]
                   + f_0 * kk_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, hk_215, hk_216, hk_217, hk_218, \
                         hk_219, kk_215, kk_216, kk_217, kk_218, \
                         kk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -4.0 * hk_215[k]
                   + f_0 * kk_215[k];

        t_216[k] = -3.0 * hk_216[k]
                   + f_0 * kk_216[k];

        t_217[k] = -3.0 * hk_217[k]
                   + f_0 * kk_217[k];

        t_218[k] = -3.0 * hk_218[k]
                   + f_0 * kk_218[k];

        t_219[k] = -3.0 * hk_219[k]
                   + f_0 * kk_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, hk_220, hk_221, hk_222, hk_223, \
                         hk_224, kk_220, kk_221, kk_222, kk_223, \
                         kk_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -3.0 * hk_220[k]
                   + f_0 * kk_220[k];

        t_221[k] = -3.0 * hk_221[k]
                   + f_0 * kk_221[k];

        t_222[k] = -3.0 * hk_222[k]
                   + f_0 * kk_222[k];

        t_223[k] = -3.0 * hk_223[k]
                   + f_0 * kk_223[k];

        t_224[k] = -3.0 * hk_224[k]
                   + f_0 * kk_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, hk_225, hk_226, hk_227, hk_228, \
                         hk_229, kk_225, kk_226, kk_227, kk_228, \
                         kk_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -3.0 * hk_225[k]
                   + f_0 * kk_225[k];

        t_226[k] = -3.0 * hk_226[k]
                   + f_0 * kk_226[k];

        t_227[k] = -3.0 * hk_227[k]
                   + f_0 * kk_227[k];

        t_228[k] = -3.0 * hk_228[k]
                   + f_0 * kk_228[k];

        t_229[k] = -3.0 * hk_229[k]
                   + f_0 * kk_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, hk_230, hk_231, hk_232, hk_233, \
                         hk_234, kk_230, kk_231, kk_232, kk_233, \
                         kk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -3.0 * hk_230[k]
                   + f_0 * kk_230[k];

        t_231[k] = -3.0 * hk_231[k]
                   + f_0 * kk_231[k];

        t_232[k] = -3.0 * hk_232[k]
                   + f_0 * kk_232[k];

        t_233[k] = -3.0 * hk_233[k]
                   + f_0 * kk_233[k];

        t_234[k] = -3.0 * hk_234[k]
                   + f_0 * kk_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, hk_235, hk_236, hk_237, hk_238, \
                         hk_239, kk_235, kk_236, kk_237, kk_238, \
                         kk_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -3.0 * hk_235[k]
                   + f_0 * kk_235[k];

        t_236[k] = -3.0 * hk_236[k]
                   + f_0 * kk_236[k];

        t_237[k] = -3.0 * hk_237[k]
                   + f_0 * kk_237[k];

        t_238[k] = -3.0 * hk_238[k]
                   + f_0 * kk_238[k];

        t_239[k] = -3.0 * hk_239[k]
                   + f_0 * kk_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, hk_240, hk_241, hk_242, hk_243, \
                         hk_244, kk_240, kk_241, kk_242, kk_243, \
                         kk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * hk_240[k]
                   + f_0 * kk_240[k];

        t_241[k] = -3.0 * hk_241[k]
                   + f_0 * kk_241[k];

        t_242[k] = -3.0 * hk_242[k]
                   + f_0 * kk_242[k];

        t_243[k] = -3.0 * hk_243[k]
                   + f_0 * kk_243[k];

        t_244[k] = -3.0 * hk_244[k]
                   + f_0 * kk_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, hk_245, hk_246, hk_247, hk_248, \
                         hk_249, kk_245, kk_246, kk_247, kk_248, \
                         kk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * hk_245[k]
                   + f_0 * kk_245[k];

        t_246[k] = -3.0 * hk_246[k]
                   + f_0 * kk_246[k];

        t_247[k] = -3.0 * hk_247[k]
                   + f_0 * kk_247[k];

        t_248[k] = -3.0 * hk_248[k]
                   + f_0 * kk_248[k];

        t_249[k] = -3.0 * hk_249[k]
                   + f_0 * kk_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, hk_250, hk_251, hk_252, hk_253, \
                         hk_254, kk_250, kk_251, kk_252, kk_253, \
                         kk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -3.0 * hk_250[k]
                   + f_0 * kk_250[k];

        t_251[k] = -3.0 * hk_251[k]
                   + f_0 * kk_251[k];

        t_252[k] = -3.0 * hk_252[k]
                   + f_0 * kk_252[k];

        t_253[k] = -3.0 * hk_253[k]
                   + f_0 * kk_253[k];

        t_254[k] = -3.0 * hk_254[k]
                   + f_0 * kk_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, hk_255, hk_256, hk_257, hk_258, \
                         hk_259, kk_255, kk_256, kk_257, kk_258, \
                         kk_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -3.0 * hk_255[k]
                   + f_0 * kk_255[k];

        t_256[k] = -3.0 * hk_256[k]
                   + f_0 * kk_256[k];

        t_257[k] = -3.0 * hk_257[k]
                   + f_0 * kk_257[k];

        t_258[k] = -3.0 * hk_258[k]
                   + f_0 * kk_258[k];

        t_259[k] = -3.0 * hk_259[k]
                   + f_0 * kk_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, hk_260, hk_261, hk_262, hk_263, \
                         hk_264, kk_260, kk_261, kk_262, kk_263, \
                         kk_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -3.0 * hk_260[k]
                   + f_0 * kk_260[k];

        t_261[k] = -3.0 * hk_261[k]
                   + f_0 * kk_261[k];

        t_262[k] = -3.0 * hk_262[k]
                   + f_0 * kk_262[k];

        t_263[k] = -3.0 * hk_263[k]
                   + f_0 * kk_263[k];

        t_264[k] = -3.0 * hk_264[k]
                   + f_0 * kk_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, hk_265, hk_266, hk_267, hk_268, \
                         hk_269, kk_265, kk_266, kk_267, kk_268, \
                         kk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -3.0 * hk_265[k]
                   + f_0 * kk_265[k];

        t_266[k] = -3.0 * hk_266[k]
                   + f_0 * kk_266[k];

        t_267[k] = -3.0 * hk_267[k]
                   + f_0 * kk_267[k];

        t_268[k] = -3.0 * hk_268[k]
                   + f_0 * kk_268[k];

        t_269[k] = -3.0 * hk_269[k]
                   + f_0 * kk_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, hk_270, hk_271, hk_272, hk_273, \
                         hk_274, kk_270, kk_271, kk_272, kk_273, \
                         kk_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -3.0 * hk_270[k]
                   + f_0 * kk_270[k];

        t_271[k] = -3.0 * hk_271[k]
                   + f_0 * kk_271[k];

        t_272[k] = -3.0 * hk_272[k]
                   + f_0 * kk_272[k];

        t_273[k] = -3.0 * hk_273[k]
                   + f_0 * kk_273[k];

        t_274[k] = -3.0 * hk_274[k]
                   + f_0 * kk_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, hk_275, hk_276, hk_277, hk_278, \
                         hk_279, kk_275, kk_276, kk_277, kk_278, \
                         kk_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * hk_275[k]
                   + f_0 * kk_275[k];

        t_276[k] = -3.0 * hk_276[k]
                   + f_0 * kk_276[k];

        t_277[k] = -3.0 * hk_277[k]
                   + f_0 * kk_277[k];

        t_278[k] = -3.0 * hk_278[k]
                   + f_0 * kk_278[k];

        t_279[k] = -3.0 * hk_279[k]
                   + f_0 * kk_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, hk_280, hk_281, hk_282, hk_283, \
                         hk_284, kk_280, kk_281, kk_282, kk_283, \
                         kk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * hk_280[k]
                   + f_0 * kk_280[k];

        t_281[k] = -3.0 * hk_281[k]
                   + f_0 * kk_281[k];

        t_282[k] = -3.0 * hk_282[k]
                   + f_0 * kk_282[k];

        t_283[k] = -3.0 * hk_283[k]
                   + f_0 * kk_283[k];

        t_284[k] = -3.0 * hk_284[k]
                   + f_0 * kk_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, hk_285, hk_286, hk_287, hk_288, \
                         hk_289, kk_285, kk_286, kk_287, kk_288, \
                         kk_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * hk_285[k]
                   + f_0 * kk_285[k];

        t_286[k] = -3.0 * hk_286[k]
                   + f_0 * kk_286[k];

        t_287[k] = -3.0 * hk_287[k]
                   + f_0 * kk_287[k];

        t_288[k] = -3.0 * hk_288[k]
                   + f_0 * kk_288[k];

        t_289[k] = -3.0 * hk_289[k]
                   + f_0 * kk_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, hk_290, hk_291, hk_292, hk_293, \
                         hk_294, kk_290, kk_291, kk_292, kk_293, \
                         kk_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * hk_290[k]
                   + f_0 * kk_290[k];

        t_291[k] = -3.0 * hk_291[k]
                   + f_0 * kk_291[k];

        t_292[k] = -3.0 * hk_292[k]
                   + f_0 * kk_292[k];

        t_293[k] = -3.0 * hk_293[k]
                   + f_0 * kk_293[k];

        t_294[k] = -3.0 * hk_294[k]
                   + f_0 * kk_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, hk_295, hk_296, hk_297, hk_298, \
                         hk_299, kk_295, kk_296, kk_297, kk_298, \
                         kk_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -3.0 * hk_295[k]
                   + f_0 * kk_295[k];

        t_296[k] = -3.0 * hk_296[k]
                   + f_0 * kk_296[k];

        t_297[k] = -3.0 * hk_297[k]
                   + f_0 * kk_297[k];

        t_298[k] = -3.0 * hk_298[k]
                   + f_0 * kk_298[k];

        t_299[k] = -3.0 * hk_299[k]
                   + f_0 * kk_299[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_300 = buffer.data(hk + 300);
    const auto *hk_301 = buffer.data(hk + 301);
    const auto *hk_302 = buffer.data(hk + 302);
    const auto *hk_303 = buffer.data(hk + 303);
    const auto *hk_304 = buffer.data(hk + 304);
    const auto *hk_305 = buffer.data(hk + 305);
    const auto *hk_306 = buffer.data(hk + 306);
    const auto *hk_307 = buffer.data(hk + 307);
    const auto *hk_308 = buffer.data(hk + 308);
    const auto *hk_309 = buffer.data(hk + 309);
    const auto *hk_310 = buffer.data(hk + 310);
    const auto *hk_311 = buffer.data(hk + 311);
    const auto *hk_312 = buffer.data(hk + 312);
    const auto *hk_313 = buffer.data(hk + 313);
    const auto *hk_314 = buffer.data(hk + 314);
    const auto *hk_315 = buffer.data(hk + 315);
    const auto *hk_316 = buffer.data(hk + 316);
    const auto *hk_317 = buffer.data(hk + 317);
    const auto *hk_318 = buffer.data(hk + 318);
    const auto *hk_319 = buffer.data(hk + 319);
    const auto *hk_320 = buffer.data(hk + 320);
    const auto *hk_321 = buffer.data(hk + 321);
    const auto *hk_322 = buffer.data(hk + 322);
    const auto *hk_323 = buffer.data(hk + 323);
    const auto *hk_324 = buffer.data(hk + 324);
    const auto *hk_325 = buffer.data(hk + 325);
    const auto *hk_326 = buffer.data(hk + 326);
    const auto *hk_327 = buffer.data(hk + 327);
    const auto *hk_328 = buffer.data(hk + 328);
    const auto *hk_329 = buffer.data(hk + 329);
    const auto *hk_330 = buffer.data(hk + 330);
    const auto *hk_331 = buffer.data(hk + 331);
    const auto *hk_332 = buffer.data(hk + 332);
    const auto *hk_333 = buffer.data(hk + 333);
    const auto *hk_334 = buffer.data(hk + 334);
    const auto *hk_335 = buffer.data(hk + 335);
    const auto *hk_336 = buffer.data(hk + 336);
    const auto *hk_337 = buffer.data(hk + 337);
    const auto *hk_338 = buffer.data(hk + 338);
    const auto *hk_339 = buffer.data(hk + 339);
    const auto *hk_340 = buffer.data(hk + 340);
    const auto *hk_341 = buffer.data(hk + 341);
    const auto *hk_342 = buffer.data(hk + 342);
    const auto *hk_343 = buffer.data(hk + 343);
    const auto *hk_344 = buffer.data(hk + 344);
    const auto *hk_345 = buffer.data(hk + 345);
    const auto *hk_346 = buffer.data(hk + 346);
    const auto *hk_347 = buffer.data(hk + 347);
    const auto *hk_348 = buffer.data(hk + 348);
    const auto *hk_349 = buffer.data(hk + 349);
    const auto *hk_350 = buffer.data(hk + 350);
    const auto *hk_351 = buffer.data(hk + 351);
    const auto *hk_352 = buffer.data(hk + 352);
    const auto *hk_353 = buffer.data(hk + 353);
    const auto *hk_354 = buffer.data(hk + 354);
    const auto *hk_355 = buffer.data(hk + 355);
    const auto *hk_356 = buffer.data(hk + 356);
    const auto *hk_357 = buffer.data(hk + 357);
    const auto *hk_358 = buffer.data(hk + 358);
    const auto *hk_359 = buffer.data(hk + 359);
    const auto *hk_360 = buffer.data(hk + 360);
    const auto *hk_361 = buffer.data(hk + 361);
    const auto *hk_362 = buffer.data(hk + 362);
    const auto *hk_363 = buffer.data(hk + 363);
    const auto *hk_364 = buffer.data(hk + 364);
    const auto *hk_365 = buffer.data(hk + 365);
    const auto *hk_366 = buffer.data(hk + 366);
    const auto *hk_367 = buffer.data(hk + 367);
    const auto *hk_368 = buffer.data(hk + 368);
    const auto *hk_369 = buffer.data(hk + 369);
    const auto *hk_370 = buffer.data(hk + 370);
    const auto *hk_371 = buffer.data(hk + 371);
    const auto *hk_372 = buffer.data(hk + 372);
    const auto *hk_373 = buffer.data(hk + 373);
    const auto *hk_374 = buffer.data(hk + 374);
    const auto *hk_375 = buffer.data(hk + 375);
    const auto *hk_376 = buffer.data(hk + 376);
    const auto *hk_377 = buffer.data(hk + 377);
    const auto *hk_378 = buffer.data(hk + 378);
    const auto *hk_379 = buffer.data(hk + 379);
    const auto *hk_380 = buffer.data(hk + 380);
    const auto *hk_381 = buffer.data(hk + 381);
    const auto *hk_382 = buffer.data(hk + 382);
    const auto *hk_383 = buffer.data(hk + 383);
    const auto *hk_384 = buffer.data(hk + 384);
    const auto *hk_385 = buffer.data(hk + 385);
    const auto *hk_386 = buffer.data(hk + 386);
    const auto *hk_387 = buffer.data(hk + 387);
    const auto *hk_388 = buffer.data(hk + 388);
    const auto *hk_389 = buffer.data(hk + 389);
    const auto *hk_390 = buffer.data(hk + 390);
    const auto *hk_391 = buffer.data(hk + 391);
    const auto *hk_392 = buffer.data(hk + 392);
    const auto *hk_393 = buffer.data(hk + 393);
    const auto *hk_394 = buffer.data(hk + 394);
    const auto *hk_395 = buffer.data(hk + 395);
    const auto *hk_396 = buffer.data(hk + 396);
    const auto *hk_397 = buffer.data(hk + 397);
    const auto *hk_398 = buffer.data(hk + 398);
    const auto *hk_399 = buffer.data(hk + 399);
    const auto *hk_400 = buffer.data(hk + 400);
    const auto *hk_401 = buffer.data(hk + 401);
    const auto *hk_402 = buffer.data(hk + 402);
    const auto *hk_403 = buffer.data(hk + 403);
    const auto *hk_404 = buffer.data(hk + 404);
    const auto *hk_405 = buffer.data(hk + 405);
    const auto *hk_406 = buffer.data(hk + 406);
    const auto *hk_407 = buffer.data(hk + 407);
    const auto *hk_408 = buffer.data(hk + 408);
    const auto *hk_409 = buffer.data(hk + 409);
    const auto *hk_410 = buffer.data(hk + 410);
    const auto *hk_411 = buffer.data(hk + 411);
    const auto *hk_412 = buffer.data(hk + 412);
    const auto *hk_413 = buffer.data(hk + 413);
    const auto *hk_414 = buffer.data(hk + 414);
    const auto *hk_415 = buffer.data(hk + 415);
    const auto *hk_416 = buffer.data(hk + 416);
    const auto *hk_417 = buffer.data(hk + 417);
    const auto *hk_418 = buffer.data(hk + 418);
    const auto *hk_419 = buffer.data(hk + 419);
    const auto *hk_420 = buffer.data(hk + 420);
    const auto *hk_421 = buffer.data(hk + 421);
    const auto *hk_422 = buffer.data(hk + 422);
    const auto *hk_423 = buffer.data(hk + 423);
    const auto *hk_424 = buffer.data(hk + 424);
    const auto *hk_425 = buffer.data(hk + 425);
    const auto *hk_426 = buffer.data(hk + 426);
    const auto *hk_427 = buffer.data(hk + 427);
    const auto *hk_428 = buffer.data(hk + 428);
    const auto *hk_429 = buffer.data(hk + 429);
    const auto *hk_430 = buffer.data(hk + 430);
    const auto *hk_431 = buffer.data(hk + 431);
    const auto *hk_432 = buffer.data(hk + 432);
    const auto *hk_433 = buffer.data(hk + 433);
    const auto *hk_434 = buffer.data(hk + 434);
    const auto *hk_435 = buffer.data(hk + 435);
    const auto *hk_436 = buffer.data(hk + 436);
    const auto *hk_437 = buffer.data(hk + 437);
    const auto *hk_438 = buffer.data(hk + 438);
    const auto *hk_439 = buffer.data(hk + 439);
    const auto *hk_440 = buffer.data(hk + 440);
    const auto *hk_441 = buffer.data(hk + 441);
    const auto *hk_442 = buffer.data(hk + 442);
    const auto *hk_443 = buffer.data(hk + 443);
    const auto *hk_444 = buffer.data(hk + 444);
    const auto *hk_445 = buffer.data(hk + 445);
    const auto *hk_446 = buffer.data(hk + 446);
    const auto *hk_447 = buffer.data(hk + 447);
    const auto *hk_448 = buffer.data(hk + 448);
    const auto *hk_449 = buffer.data(hk + 449);

    const auto *kk_300 = buffer.data(kk + 300);
    const auto *kk_301 = buffer.data(kk + 301);
    const auto *kk_302 = buffer.data(kk + 302);
    const auto *kk_303 = buffer.data(kk + 303);
    const auto *kk_304 = buffer.data(kk + 304);
    const auto *kk_305 = buffer.data(kk + 305);
    const auto *kk_306 = buffer.data(kk + 306);
    const auto *kk_307 = buffer.data(kk + 307);
    const auto *kk_308 = buffer.data(kk + 308);
    const auto *kk_309 = buffer.data(kk + 309);
    const auto *kk_310 = buffer.data(kk + 310);
    const auto *kk_311 = buffer.data(kk + 311);
    const auto *kk_312 = buffer.data(kk + 312);
    const auto *kk_313 = buffer.data(kk + 313);
    const auto *kk_314 = buffer.data(kk + 314);
    const auto *kk_315 = buffer.data(kk + 315);
    const auto *kk_316 = buffer.data(kk + 316);
    const auto *kk_317 = buffer.data(kk + 317);
    const auto *kk_318 = buffer.data(kk + 318);
    const auto *kk_319 = buffer.data(kk + 319);
    const auto *kk_320 = buffer.data(kk + 320);
    const auto *kk_321 = buffer.data(kk + 321);
    const auto *kk_322 = buffer.data(kk + 322);
    const auto *kk_323 = buffer.data(kk + 323);
    const auto *kk_324 = buffer.data(kk + 324);
    const auto *kk_325 = buffer.data(kk + 325);
    const auto *kk_326 = buffer.data(kk + 326);
    const auto *kk_327 = buffer.data(kk + 327);
    const auto *kk_328 = buffer.data(kk + 328);
    const auto *kk_329 = buffer.data(kk + 329);
    const auto *kk_330 = buffer.data(kk + 330);
    const auto *kk_331 = buffer.data(kk + 331);
    const auto *kk_332 = buffer.data(kk + 332);
    const auto *kk_333 = buffer.data(kk + 333);
    const auto *kk_334 = buffer.data(kk + 334);
    const auto *kk_335 = buffer.data(kk + 335);
    const auto *kk_336 = buffer.data(kk + 336);
    const auto *kk_337 = buffer.data(kk + 337);
    const auto *kk_338 = buffer.data(kk + 338);
    const auto *kk_339 = buffer.data(kk + 339);
    const auto *kk_340 = buffer.data(kk + 340);
    const auto *kk_341 = buffer.data(kk + 341);
    const auto *kk_342 = buffer.data(kk + 342);
    const auto *kk_343 = buffer.data(kk + 343);
    const auto *kk_344 = buffer.data(kk + 344);
    const auto *kk_345 = buffer.data(kk + 345);
    const auto *kk_346 = buffer.data(kk + 346);
    const auto *kk_347 = buffer.data(kk + 347);
    const auto *kk_348 = buffer.data(kk + 348);
    const auto *kk_349 = buffer.data(kk + 349);
    const auto *kk_350 = buffer.data(kk + 350);
    const auto *kk_351 = buffer.data(kk + 351);
    const auto *kk_352 = buffer.data(kk + 352);
    const auto *kk_353 = buffer.data(kk + 353);
    const auto *kk_354 = buffer.data(kk + 354);
    const auto *kk_355 = buffer.data(kk + 355);
    const auto *kk_356 = buffer.data(kk + 356);
    const auto *kk_357 = buffer.data(kk + 357);
    const auto *kk_358 = buffer.data(kk + 358);
    const auto *kk_359 = buffer.data(kk + 359);
    const auto *kk_360 = buffer.data(kk + 360);
    const auto *kk_361 = buffer.data(kk + 361);
    const auto *kk_362 = buffer.data(kk + 362);
    const auto *kk_363 = buffer.data(kk + 363);
    const auto *kk_364 = buffer.data(kk + 364);
    const auto *kk_365 = buffer.data(kk + 365);
    const auto *kk_366 = buffer.data(kk + 366);
    const auto *kk_367 = buffer.data(kk + 367);
    const auto *kk_368 = buffer.data(kk + 368);
    const auto *kk_369 = buffer.data(kk + 369);
    const auto *kk_370 = buffer.data(kk + 370);
    const auto *kk_371 = buffer.data(kk + 371);
    const auto *kk_372 = buffer.data(kk + 372);
    const auto *kk_373 = buffer.data(kk + 373);
    const auto *kk_374 = buffer.data(kk + 374);
    const auto *kk_375 = buffer.data(kk + 375);
    const auto *kk_376 = buffer.data(kk + 376);
    const auto *kk_377 = buffer.data(kk + 377);
    const auto *kk_378 = buffer.data(kk + 378);
    const auto *kk_379 = buffer.data(kk + 379);
    const auto *kk_380 = buffer.data(kk + 380);
    const auto *kk_381 = buffer.data(kk + 381);
    const auto *kk_382 = buffer.data(kk + 382);
    const auto *kk_383 = buffer.data(kk + 383);
    const auto *kk_384 = buffer.data(kk + 384);
    const auto *kk_385 = buffer.data(kk + 385);
    const auto *kk_386 = buffer.data(kk + 386);
    const auto *kk_387 = buffer.data(kk + 387);
    const auto *kk_388 = buffer.data(kk + 388);
    const auto *kk_389 = buffer.data(kk + 389);
    const auto *kk_390 = buffer.data(kk + 390);
    const auto *kk_391 = buffer.data(kk + 391);
    const auto *kk_392 = buffer.data(kk + 392);
    const auto *kk_393 = buffer.data(kk + 393);
    const auto *kk_394 = buffer.data(kk + 394);
    const auto *kk_395 = buffer.data(kk + 395);
    const auto *kk_396 = buffer.data(kk + 396);
    const auto *kk_397 = buffer.data(kk + 397);
    const auto *kk_398 = buffer.data(kk + 398);
    const auto *kk_399 = buffer.data(kk + 399);
    const auto *kk_400 = buffer.data(kk + 400);
    const auto *kk_401 = buffer.data(kk + 401);
    const auto *kk_402 = buffer.data(kk + 402);
    const auto *kk_403 = buffer.data(kk + 403);
    const auto *kk_404 = buffer.data(kk + 404);
    const auto *kk_405 = buffer.data(kk + 405);
    const auto *kk_406 = buffer.data(kk + 406);
    const auto *kk_407 = buffer.data(kk + 407);
    const auto *kk_408 = buffer.data(kk + 408);
    const auto *kk_409 = buffer.data(kk + 409);
    const auto *kk_410 = buffer.data(kk + 410);
    const auto *kk_411 = buffer.data(kk + 411);
    const auto *kk_412 = buffer.data(kk + 412);
    const auto *kk_413 = buffer.data(kk + 413);
    const auto *kk_414 = buffer.data(kk + 414);
    const auto *kk_415 = buffer.data(kk + 415);
    const auto *kk_416 = buffer.data(kk + 416);
    const auto *kk_417 = buffer.data(kk + 417);
    const auto *kk_418 = buffer.data(kk + 418);
    const auto *kk_419 = buffer.data(kk + 419);
    const auto *kk_420 = buffer.data(kk + 420);
    const auto *kk_421 = buffer.data(kk + 421);
    const auto *kk_422 = buffer.data(kk + 422);
    const auto *kk_423 = buffer.data(kk + 423);
    const auto *kk_424 = buffer.data(kk + 424);
    const auto *kk_425 = buffer.data(kk + 425);
    const auto *kk_426 = buffer.data(kk + 426);
    const auto *kk_427 = buffer.data(kk + 427);
    const auto *kk_428 = buffer.data(kk + 428);
    const auto *kk_429 = buffer.data(kk + 429);
    const auto *kk_430 = buffer.data(kk + 430);
    const auto *kk_431 = buffer.data(kk + 431);
    const auto *kk_432 = buffer.data(kk + 432);
    const auto *kk_433 = buffer.data(kk + 433);
    const auto *kk_434 = buffer.data(kk + 434);
    const auto *kk_435 = buffer.data(kk + 435);
    const auto *kk_436 = buffer.data(kk + 436);
    const auto *kk_437 = buffer.data(kk + 437);
    const auto *kk_438 = buffer.data(kk + 438);
    const auto *kk_439 = buffer.data(kk + 439);
    const auto *kk_440 = buffer.data(kk + 440);
    const auto *kk_441 = buffer.data(kk + 441);
    const auto *kk_442 = buffer.data(kk + 442);
    const auto *kk_443 = buffer.data(kk + 443);
    const auto *kk_444 = buffer.data(kk + 444);
    const auto *kk_445 = buffer.data(kk + 445);
    const auto *kk_446 = buffer.data(kk + 446);
    const auto *kk_447 = buffer.data(kk + 447);
    const auto *kk_448 = buffer.data(kk + 448);
    const auto *kk_449 = buffer.data(kk + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, hk_300, hk_301, hk_302, hk_303, \
                         hk_304, kk_300, kk_301, kk_302, kk_303, \
                         kk_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -3.0 * hk_300[k]
                   + f_0 * kk_300[k];

        t_301[k] = -3.0 * hk_301[k]
                   + f_0 * kk_301[k];

        t_302[k] = -3.0 * hk_302[k]
                   + f_0 * kk_302[k];

        t_303[k] = -3.0 * hk_303[k]
                   + f_0 * kk_303[k];

        t_304[k] = -3.0 * hk_304[k]
                   + f_0 * kk_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, hk_305, hk_306, hk_307, hk_308, \
                         hk_309, kk_305, kk_306, kk_307, kk_308, \
                         kk_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -3.0 * hk_305[k]
                   + f_0 * kk_305[k];

        t_306[k] = -3.0 * hk_306[k]
                   + f_0 * kk_306[k];

        t_307[k] = -3.0 * hk_307[k]
                   + f_0 * kk_307[k];

        t_308[k] = -3.0 * hk_308[k]
                   + f_0 * kk_308[k];

        t_309[k] = -3.0 * hk_309[k]
                   + f_0 * kk_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, hk_310, hk_311, hk_312, hk_313, \
                         hk_314, kk_310, kk_311, kk_312, kk_313, \
                         kk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * hk_310[k]
                   + f_0 * kk_310[k];

        t_311[k] = -3.0 * hk_311[k]
                   + f_0 * kk_311[k];

        t_312[k] = -3.0 * hk_312[k]
                   + f_0 * kk_312[k];

        t_313[k] = -3.0 * hk_313[k]
                   + f_0 * kk_313[k];

        t_314[k] = -3.0 * hk_314[k]
                   + f_0 * kk_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, hk_315, hk_316, hk_317, hk_318, \
                         hk_319, kk_315, kk_316, kk_317, kk_318, \
                         kk_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -3.0 * hk_315[k]
                   + f_0 * kk_315[k];

        t_316[k] = -3.0 * hk_316[k]
                   + f_0 * kk_316[k];

        t_317[k] = -3.0 * hk_317[k]
                   + f_0 * kk_317[k];

        t_318[k] = -3.0 * hk_318[k]
                   + f_0 * kk_318[k];

        t_319[k] = -3.0 * hk_319[k]
                   + f_0 * kk_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, hk_320, hk_321, hk_322, hk_323, \
                         hk_324, kk_320, kk_321, kk_322, kk_323, \
                         kk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -3.0 * hk_320[k]
                   + f_0 * kk_320[k];

        t_321[k] = -3.0 * hk_321[k]
                   + f_0 * kk_321[k];

        t_322[k] = -3.0 * hk_322[k]
                   + f_0 * kk_322[k];

        t_323[k] = -3.0 * hk_323[k]
                   + f_0 * kk_323[k];

        t_324[k] = -3.0 * hk_324[k]
                   + f_0 * kk_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, hk_325, hk_326, hk_327, hk_328, \
                         hk_329, kk_325, kk_326, kk_327, kk_328, \
                         kk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -3.0 * hk_325[k]
                   + f_0 * kk_325[k];

        t_326[k] = -3.0 * hk_326[k]
                   + f_0 * kk_326[k];

        t_327[k] = -3.0 * hk_327[k]
                   + f_0 * kk_327[k];

        t_328[k] = -3.0 * hk_328[k]
                   + f_0 * kk_328[k];

        t_329[k] = -3.0 * hk_329[k]
                   + f_0 * kk_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, hk_330, hk_331, hk_332, hk_333, \
                         hk_334, kk_330, kk_331, kk_332, kk_333, \
                         kk_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -3.0 * hk_330[k]
                   + f_0 * kk_330[k];

        t_331[k] = -3.0 * hk_331[k]
                   + f_0 * kk_331[k];

        t_332[k] = -3.0 * hk_332[k]
                   + f_0 * kk_332[k];

        t_333[k] = -3.0 * hk_333[k]
                   + f_0 * kk_333[k];

        t_334[k] = -3.0 * hk_334[k]
                   + f_0 * kk_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, hk_335, hk_336, hk_337, hk_338, \
                         hk_339, kk_335, kk_336, kk_337, kk_338, \
                         kk_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -3.0 * hk_335[k]
                   + f_0 * kk_335[k];

        t_336[k] = -3.0 * hk_336[k]
                   + f_0 * kk_336[k];

        t_337[k] = -3.0 * hk_337[k]
                   + f_0 * kk_337[k];

        t_338[k] = -3.0 * hk_338[k]
                   + f_0 * kk_338[k];

        t_339[k] = -3.0 * hk_339[k]
                   + f_0 * kk_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, hk_340, hk_341, hk_342, hk_343, \
                         hk_344, kk_340, kk_341, kk_342, kk_343, \
                         kk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -3.0 * hk_340[k]
                   + f_0 * kk_340[k];

        t_341[k] = -3.0 * hk_341[k]
                   + f_0 * kk_341[k];

        t_342[k] = -3.0 * hk_342[k]
                   + f_0 * kk_342[k];

        t_343[k] = -3.0 * hk_343[k]
                   + f_0 * kk_343[k];

        t_344[k] = -3.0 * hk_344[k]
                   + f_0 * kk_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, hk_345, hk_346, hk_347, hk_348, \
                         hk_349, kk_345, kk_346, kk_347, kk_348, \
                         kk_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -3.0 * hk_345[k]
                   + f_0 * kk_345[k];

        t_346[k] = -3.0 * hk_346[k]
                   + f_0 * kk_346[k];

        t_347[k] = -3.0 * hk_347[k]
                   + f_0 * kk_347[k];

        t_348[k] = -3.0 * hk_348[k]
                   + f_0 * kk_348[k];

        t_349[k] = -3.0 * hk_349[k]
                   + f_0 * kk_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, hk_350, hk_351, hk_352, hk_353, \
                         hk_354, kk_350, kk_351, kk_352, kk_353, \
                         kk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -3.0 * hk_350[k]
                   + f_0 * kk_350[k];

        t_351[k] = -3.0 * hk_351[k]
                   + f_0 * kk_351[k];

        t_352[k] = -3.0 * hk_352[k]
                   + f_0 * kk_352[k];

        t_353[k] = -3.0 * hk_353[k]
                   + f_0 * kk_353[k];

        t_354[k] = -3.0 * hk_354[k]
                   + f_0 * kk_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, hk_355, hk_356, hk_357, hk_358, \
                         hk_359, kk_355, kk_356, kk_357, kk_358, \
                         kk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -3.0 * hk_355[k]
                   + f_0 * kk_355[k];

        t_356[k] = -3.0 * hk_356[k]
                   + f_0 * kk_356[k];

        t_357[k] = -3.0 * hk_357[k]
                   + f_0 * kk_357[k];

        t_358[k] = -3.0 * hk_358[k]
                   + f_0 * kk_358[k];

        t_359[k] = -3.0 * hk_359[k]
                   + f_0 * kk_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, hk_360, hk_361, hk_362, hk_363, \
                         hk_364, kk_360, kk_361, kk_362, kk_363, \
                         kk_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -2.0 * hk_360[k]
                   + f_0 * kk_360[k];

        t_361[k] = -2.0 * hk_361[k]
                   + f_0 * kk_361[k];

        t_362[k] = -2.0 * hk_362[k]
                   + f_0 * kk_362[k];

        t_363[k] = -2.0 * hk_363[k]
                   + f_0 * kk_363[k];

        t_364[k] = -2.0 * hk_364[k]
                   + f_0 * kk_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, hk_365, hk_366, hk_367, hk_368, \
                         hk_369, kk_365, kk_366, kk_367, kk_368, \
                         kk_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -2.0 * hk_365[k]
                   + f_0 * kk_365[k];

        t_366[k] = -2.0 * hk_366[k]
                   + f_0 * kk_366[k];

        t_367[k] = -2.0 * hk_367[k]
                   + f_0 * kk_367[k];

        t_368[k] = -2.0 * hk_368[k]
                   + f_0 * kk_368[k];

        t_369[k] = -2.0 * hk_369[k]
                   + f_0 * kk_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, hk_370, hk_371, hk_372, hk_373, \
                         hk_374, kk_370, kk_371, kk_372, kk_373, \
                         kk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -2.0 * hk_370[k]
                   + f_0 * kk_370[k];

        t_371[k] = -2.0 * hk_371[k]
                   + f_0 * kk_371[k];

        t_372[k] = -2.0 * hk_372[k]
                   + f_0 * kk_372[k];

        t_373[k] = -2.0 * hk_373[k]
                   + f_0 * kk_373[k];

        t_374[k] = -2.0 * hk_374[k]
                   + f_0 * kk_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, hk_375, hk_376, hk_377, hk_378, \
                         hk_379, kk_375, kk_376, kk_377, kk_378, \
                         kk_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -2.0 * hk_375[k]
                   + f_0 * kk_375[k];

        t_376[k] = -2.0 * hk_376[k]
                   + f_0 * kk_376[k];

        t_377[k] = -2.0 * hk_377[k]
                   + f_0 * kk_377[k];

        t_378[k] = -2.0 * hk_378[k]
                   + f_0 * kk_378[k];

        t_379[k] = -2.0 * hk_379[k]
                   + f_0 * kk_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, hk_380, hk_381, hk_382, hk_383, \
                         hk_384, kk_380, kk_381, kk_382, kk_383, \
                         kk_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * hk_380[k]
                   + f_0 * kk_380[k];

        t_381[k] = -2.0 * hk_381[k]
                   + f_0 * kk_381[k];

        t_382[k] = -2.0 * hk_382[k]
                   + f_0 * kk_382[k];

        t_383[k] = -2.0 * hk_383[k]
                   + f_0 * kk_383[k];

        t_384[k] = -2.0 * hk_384[k]
                   + f_0 * kk_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, hk_385, hk_386, hk_387, hk_388, \
                         hk_389, kk_385, kk_386, kk_387, kk_388, \
                         kk_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * hk_385[k]
                   + f_0 * kk_385[k];

        t_386[k] = -2.0 * hk_386[k]
                   + f_0 * kk_386[k];

        t_387[k] = -2.0 * hk_387[k]
                   + f_0 * kk_387[k];

        t_388[k] = -2.0 * hk_388[k]
                   + f_0 * kk_388[k];

        t_389[k] = -2.0 * hk_389[k]
                   + f_0 * kk_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, hk_390, hk_391, hk_392, hk_393, \
                         hk_394, kk_390, kk_391, kk_392, kk_393, \
                         kk_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -2.0 * hk_390[k]
                   + f_0 * kk_390[k];

        t_391[k] = -2.0 * hk_391[k]
                   + f_0 * kk_391[k];

        t_392[k] = -2.0 * hk_392[k]
                   + f_0 * kk_392[k];

        t_393[k] = -2.0 * hk_393[k]
                   + f_0 * kk_393[k];

        t_394[k] = -2.0 * hk_394[k]
                   + f_0 * kk_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, hk_395, hk_396, hk_397, hk_398, \
                         hk_399, kk_395, kk_396, kk_397, kk_398, \
                         kk_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -2.0 * hk_395[k]
                   + f_0 * kk_395[k];

        t_396[k] = -2.0 * hk_396[k]
                   + f_0 * kk_396[k];

        t_397[k] = -2.0 * hk_397[k]
                   + f_0 * kk_397[k];

        t_398[k] = -2.0 * hk_398[k]
                   + f_0 * kk_398[k];

        t_399[k] = -2.0 * hk_399[k]
                   + f_0 * kk_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, hk_400, hk_401, hk_402, hk_403, \
                         hk_404, kk_400, kk_401, kk_402, kk_403, \
                         kk_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -2.0 * hk_400[k]
                   + f_0 * kk_400[k];

        t_401[k] = -2.0 * hk_401[k]
                   + f_0 * kk_401[k];

        t_402[k] = -2.0 * hk_402[k]
                   + f_0 * kk_402[k];

        t_403[k] = -2.0 * hk_403[k]
                   + f_0 * kk_403[k];

        t_404[k] = -2.0 * hk_404[k]
                   + f_0 * kk_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, hk_405, hk_406, hk_407, hk_408, \
                         hk_409, kk_405, kk_406, kk_407, kk_408, \
                         kk_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -2.0 * hk_405[k]
                   + f_0 * kk_405[k];

        t_406[k] = -2.0 * hk_406[k]
                   + f_0 * kk_406[k];

        t_407[k] = -2.0 * hk_407[k]
                   + f_0 * kk_407[k];

        t_408[k] = -2.0 * hk_408[k]
                   + f_0 * kk_408[k];

        t_409[k] = -2.0 * hk_409[k]
                   + f_0 * kk_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, hk_410, hk_411, hk_412, hk_413, \
                         hk_414, kk_410, kk_411, kk_412, kk_413, \
                         kk_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -2.0 * hk_410[k]
                   + f_0 * kk_410[k];

        t_411[k] = -2.0 * hk_411[k]
                   + f_0 * kk_411[k];

        t_412[k] = -2.0 * hk_412[k]
                   + f_0 * kk_412[k];

        t_413[k] = -2.0 * hk_413[k]
                   + f_0 * kk_413[k];

        t_414[k] = -2.0 * hk_414[k]
                   + f_0 * kk_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, hk_415, hk_416, hk_417, hk_418, \
                         hk_419, kk_415, kk_416, kk_417, kk_418, \
                         kk_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -2.0 * hk_415[k]
                   + f_0 * kk_415[k];

        t_416[k] = -2.0 * hk_416[k]
                   + f_0 * kk_416[k];

        t_417[k] = -2.0 * hk_417[k]
                   + f_0 * kk_417[k];

        t_418[k] = -2.0 * hk_418[k]
                   + f_0 * kk_418[k];

        t_419[k] = -2.0 * hk_419[k]
                   + f_0 * kk_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, hk_420, hk_421, hk_422, hk_423, \
                         hk_424, kk_420, kk_421, kk_422, kk_423, \
                         kk_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -2.0 * hk_420[k]
                   + f_0 * kk_420[k];

        t_421[k] = -2.0 * hk_421[k]
                   + f_0 * kk_421[k];

        t_422[k] = -2.0 * hk_422[k]
                   + f_0 * kk_422[k];

        t_423[k] = -2.0 * hk_423[k]
                   + f_0 * kk_423[k];

        t_424[k] = -2.0 * hk_424[k]
                   + f_0 * kk_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, hk_425, hk_426, hk_427, hk_428, \
                         hk_429, kk_425, kk_426, kk_427, kk_428, \
                         kk_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -2.0 * hk_425[k]
                   + f_0 * kk_425[k];

        t_426[k] = -2.0 * hk_426[k]
                   + f_0 * kk_426[k];

        t_427[k] = -2.0 * hk_427[k]
                   + f_0 * kk_427[k];

        t_428[k] = -2.0 * hk_428[k]
                   + f_0 * kk_428[k];

        t_429[k] = -2.0 * hk_429[k]
                   + f_0 * kk_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, hk_430, hk_431, hk_432, hk_433, \
                         hk_434, kk_430, kk_431, kk_432, kk_433, \
                         kk_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -2.0 * hk_430[k]
                   + f_0 * kk_430[k];

        t_431[k] = -2.0 * hk_431[k]
                   + f_0 * kk_431[k];

        t_432[k] = -2.0 * hk_432[k]
                   + f_0 * kk_432[k];

        t_433[k] = -2.0 * hk_433[k]
                   + f_0 * kk_433[k];

        t_434[k] = -2.0 * hk_434[k]
                   + f_0 * kk_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, hk_435, hk_436, hk_437, hk_438, \
                         hk_439, kk_435, kk_436, kk_437, kk_438, \
                         kk_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -2.0 * hk_435[k]
                   + f_0 * kk_435[k];

        t_436[k] = -2.0 * hk_436[k]
                   + f_0 * kk_436[k];

        t_437[k] = -2.0 * hk_437[k]
                   + f_0 * kk_437[k];

        t_438[k] = -2.0 * hk_438[k]
                   + f_0 * kk_438[k];

        t_439[k] = -2.0 * hk_439[k]
                   + f_0 * kk_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, hk_440, hk_441, hk_442, hk_443, \
                         hk_444, kk_440, kk_441, kk_442, kk_443, \
                         kk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -2.0 * hk_440[k]
                   + f_0 * kk_440[k];

        t_441[k] = -2.0 * hk_441[k]
                   + f_0 * kk_441[k];

        t_442[k] = -2.0 * hk_442[k]
                   + f_0 * kk_442[k];

        t_443[k] = -2.0 * hk_443[k]
                   + f_0 * kk_443[k];

        t_444[k] = -2.0 * hk_444[k]
                   + f_0 * kk_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, hk_445, hk_446, hk_447, hk_448, \
                         hk_449, kk_445, kk_446, kk_447, kk_448, \
                         kk_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -2.0 * hk_445[k]
                   + f_0 * kk_445[k];

        t_446[k] = -2.0 * hk_446[k]
                   + f_0 * kk_446[k];

        t_447[k] = -2.0 * hk_447[k]
                   + f_0 * kk_447[k];

        t_448[k] = -2.0 * hk_448[k]
                   + f_0 * kk_448[k];

        t_449[k] = -2.0 * hk_449[k]
                   + f_0 * kk_449[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_450 = buffer.data(hk + 450);
    const auto *hk_451 = buffer.data(hk + 451);
    const auto *hk_452 = buffer.data(hk + 452);
    const auto *hk_453 = buffer.data(hk + 453);
    const auto *hk_454 = buffer.data(hk + 454);
    const auto *hk_455 = buffer.data(hk + 455);
    const auto *hk_456 = buffer.data(hk + 456);
    const auto *hk_457 = buffer.data(hk + 457);
    const auto *hk_458 = buffer.data(hk + 458);
    const auto *hk_459 = buffer.data(hk + 459);
    const auto *hk_460 = buffer.data(hk + 460);
    const auto *hk_461 = buffer.data(hk + 461);
    const auto *hk_462 = buffer.data(hk + 462);
    const auto *hk_463 = buffer.data(hk + 463);
    const auto *hk_464 = buffer.data(hk + 464);
    const auto *hk_465 = buffer.data(hk + 465);
    const auto *hk_466 = buffer.data(hk + 466);
    const auto *hk_467 = buffer.data(hk + 467);
    const auto *hk_468 = buffer.data(hk + 468);
    const auto *hk_469 = buffer.data(hk + 469);
    const auto *hk_470 = buffer.data(hk + 470);
    const auto *hk_471 = buffer.data(hk + 471);
    const auto *hk_472 = buffer.data(hk + 472);
    const auto *hk_473 = buffer.data(hk + 473);
    const auto *hk_474 = buffer.data(hk + 474);
    const auto *hk_475 = buffer.data(hk + 475);
    const auto *hk_476 = buffer.data(hk + 476);
    const auto *hk_477 = buffer.data(hk + 477);
    const auto *hk_478 = buffer.data(hk + 478);
    const auto *hk_479 = buffer.data(hk + 479);
    const auto *hk_480 = buffer.data(hk + 480);
    const auto *hk_481 = buffer.data(hk + 481);
    const auto *hk_482 = buffer.data(hk + 482);
    const auto *hk_483 = buffer.data(hk + 483);
    const auto *hk_484 = buffer.data(hk + 484);
    const auto *hk_485 = buffer.data(hk + 485);
    const auto *hk_486 = buffer.data(hk + 486);
    const auto *hk_487 = buffer.data(hk + 487);
    const auto *hk_488 = buffer.data(hk + 488);
    const auto *hk_489 = buffer.data(hk + 489);
    const auto *hk_490 = buffer.data(hk + 490);
    const auto *hk_491 = buffer.data(hk + 491);
    const auto *hk_492 = buffer.data(hk + 492);
    const auto *hk_493 = buffer.data(hk + 493);
    const auto *hk_494 = buffer.data(hk + 494);
    const auto *hk_495 = buffer.data(hk + 495);
    const auto *hk_496 = buffer.data(hk + 496);
    const auto *hk_497 = buffer.data(hk + 497);
    const auto *hk_498 = buffer.data(hk + 498);
    const auto *hk_499 = buffer.data(hk + 499);
    const auto *hk_500 = buffer.data(hk + 500);
    const auto *hk_501 = buffer.data(hk + 501);
    const auto *hk_502 = buffer.data(hk + 502);
    const auto *hk_503 = buffer.data(hk + 503);
    const auto *hk_504 = buffer.data(hk + 504);
    const auto *hk_505 = buffer.data(hk + 505);
    const auto *hk_506 = buffer.data(hk + 506);
    const auto *hk_507 = buffer.data(hk + 507);
    const auto *hk_508 = buffer.data(hk + 508);
    const auto *hk_509 = buffer.data(hk + 509);
    const auto *hk_510 = buffer.data(hk + 510);
    const auto *hk_511 = buffer.data(hk + 511);
    const auto *hk_512 = buffer.data(hk + 512);
    const auto *hk_513 = buffer.data(hk + 513);
    const auto *hk_514 = buffer.data(hk + 514);
    const auto *hk_515 = buffer.data(hk + 515);
    const auto *hk_516 = buffer.data(hk + 516);
    const auto *hk_517 = buffer.data(hk + 517);
    const auto *hk_518 = buffer.data(hk + 518);
    const auto *hk_519 = buffer.data(hk + 519);
    const auto *hk_520 = buffer.data(hk + 520);
    const auto *hk_521 = buffer.data(hk + 521);
    const auto *hk_522 = buffer.data(hk + 522);
    const auto *hk_523 = buffer.data(hk + 523);
    const auto *hk_524 = buffer.data(hk + 524);
    const auto *hk_525 = buffer.data(hk + 525);
    const auto *hk_526 = buffer.data(hk + 526);
    const auto *hk_527 = buffer.data(hk + 527);
    const auto *hk_528 = buffer.data(hk + 528);
    const auto *hk_529 = buffer.data(hk + 529);
    const auto *hk_530 = buffer.data(hk + 530);
    const auto *hk_531 = buffer.data(hk + 531);
    const auto *hk_532 = buffer.data(hk + 532);
    const auto *hk_533 = buffer.data(hk + 533);
    const auto *hk_534 = buffer.data(hk + 534);
    const auto *hk_535 = buffer.data(hk + 535);
    const auto *hk_536 = buffer.data(hk + 536);
    const auto *hk_537 = buffer.data(hk + 537);
    const auto *hk_538 = buffer.data(hk + 538);
    const auto *hk_539 = buffer.data(hk + 539);
    const auto *hk_540 = buffer.data(hk + 540);
    const auto *hk_541 = buffer.data(hk + 541);
    const auto *hk_542 = buffer.data(hk + 542);
    const auto *hk_543 = buffer.data(hk + 543);
    const auto *hk_544 = buffer.data(hk + 544);
    const auto *hk_545 = buffer.data(hk + 545);
    const auto *hk_546 = buffer.data(hk + 546);
    const auto *hk_547 = buffer.data(hk + 547);
    const auto *hk_548 = buffer.data(hk + 548);
    const auto *hk_549 = buffer.data(hk + 549);
    const auto *hk_550 = buffer.data(hk + 550);
    const auto *hk_551 = buffer.data(hk + 551);
    const auto *hk_552 = buffer.data(hk + 552);
    const auto *hk_553 = buffer.data(hk + 553);
    const auto *hk_554 = buffer.data(hk + 554);
    const auto *hk_555 = buffer.data(hk + 555);
    const auto *hk_556 = buffer.data(hk + 556);
    const auto *hk_557 = buffer.data(hk + 557);
    const auto *hk_558 = buffer.data(hk + 558);
    const auto *hk_559 = buffer.data(hk + 559);
    const auto *hk_560 = buffer.data(hk + 560);
    const auto *hk_561 = buffer.data(hk + 561);
    const auto *hk_562 = buffer.data(hk + 562);
    const auto *hk_563 = buffer.data(hk + 563);
    const auto *hk_564 = buffer.data(hk + 564);
    const auto *hk_565 = buffer.data(hk + 565);
    const auto *hk_566 = buffer.data(hk + 566);
    const auto *hk_567 = buffer.data(hk + 567);
    const auto *hk_568 = buffer.data(hk + 568);
    const auto *hk_569 = buffer.data(hk + 569);
    const auto *hk_570 = buffer.data(hk + 570);
    const auto *hk_571 = buffer.data(hk + 571);
    const auto *hk_572 = buffer.data(hk + 572);
    const auto *hk_573 = buffer.data(hk + 573);
    const auto *hk_574 = buffer.data(hk + 574);
    const auto *hk_575 = buffer.data(hk + 575);
    const auto *hk_576 = buffer.data(hk + 576);
    const auto *hk_577 = buffer.data(hk + 577);
    const auto *hk_578 = buffer.data(hk + 578);
    const auto *hk_579 = buffer.data(hk + 579);
    const auto *hk_580 = buffer.data(hk + 580);
    const auto *hk_581 = buffer.data(hk + 581);
    const auto *hk_582 = buffer.data(hk + 582);
    const auto *hk_583 = buffer.data(hk + 583);
    const auto *hk_584 = buffer.data(hk + 584);
    const auto *hk_585 = buffer.data(hk + 585);
    const auto *hk_586 = buffer.data(hk + 586);
    const auto *hk_587 = buffer.data(hk + 587);
    const auto *hk_588 = buffer.data(hk + 588);
    const auto *hk_589 = buffer.data(hk + 589);
    const auto *hk_590 = buffer.data(hk + 590);
    const auto *hk_591 = buffer.data(hk + 591);
    const auto *hk_592 = buffer.data(hk + 592);
    const auto *hk_593 = buffer.data(hk + 593);
    const auto *hk_594 = buffer.data(hk + 594);
    const auto *hk_595 = buffer.data(hk + 595);
    const auto *hk_596 = buffer.data(hk + 596);
    const auto *hk_597 = buffer.data(hk + 597);
    const auto *hk_598 = buffer.data(hk + 598);
    const auto *hk_599 = buffer.data(hk + 599);

    const auto *kk_450 = buffer.data(kk + 450);
    const auto *kk_451 = buffer.data(kk + 451);
    const auto *kk_452 = buffer.data(kk + 452);
    const auto *kk_453 = buffer.data(kk + 453);
    const auto *kk_454 = buffer.data(kk + 454);
    const auto *kk_455 = buffer.data(kk + 455);
    const auto *kk_456 = buffer.data(kk + 456);
    const auto *kk_457 = buffer.data(kk + 457);
    const auto *kk_458 = buffer.data(kk + 458);
    const auto *kk_459 = buffer.data(kk + 459);
    const auto *kk_460 = buffer.data(kk + 460);
    const auto *kk_461 = buffer.data(kk + 461);
    const auto *kk_462 = buffer.data(kk + 462);
    const auto *kk_463 = buffer.data(kk + 463);
    const auto *kk_464 = buffer.data(kk + 464);
    const auto *kk_465 = buffer.data(kk + 465);
    const auto *kk_466 = buffer.data(kk + 466);
    const auto *kk_467 = buffer.data(kk + 467);
    const auto *kk_468 = buffer.data(kk + 468);
    const auto *kk_469 = buffer.data(kk + 469);
    const auto *kk_470 = buffer.data(kk + 470);
    const auto *kk_471 = buffer.data(kk + 471);
    const auto *kk_472 = buffer.data(kk + 472);
    const auto *kk_473 = buffer.data(kk + 473);
    const auto *kk_474 = buffer.data(kk + 474);
    const auto *kk_475 = buffer.data(kk + 475);
    const auto *kk_476 = buffer.data(kk + 476);
    const auto *kk_477 = buffer.data(kk + 477);
    const auto *kk_478 = buffer.data(kk + 478);
    const auto *kk_479 = buffer.data(kk + 479);
    const auto *kk_480 = buffer.data(kk + 480);
    const auto *kk_481 = buffer.data(kk + 481);
    const auto *kk_482 = buffer.data(kk + 482);
    const auto *kk_483 = buffer.data(kk + 483);
    const auto *kk_484 = buffer.data(kk + 484);
    const auto *kk_485 = buffer.data(kk + 485);
    const auto *kk_486 = buffer.data(kk + 486);
    const auto *kk_487 = buffer.data(kk + 487);
    const auto *kk_488 = buffer.data(kk + 488);
    const auto *kk_489 = buffer.data(kk + 489);
    const auto *kk_490 = buffer.data(kk + 490);
    const auto *kk_491 = buffer.data(kk + 491);
    const auto *kk_492 = buffer.data(kk + 492);
    const auto *kk_493 = buffer.data(kk + 493);
    const auto *kk_494 = buffer.data(kk + 494);
    const auto *kk_495 = buffer.data(kk + 495);
    const auto *kk_496 = buffer.data(kk + 496);
    const auto *kk_497 = buffer.data(kk + 497);
    const auto *kk_498 = buffer.data(kk + 498);
    const auto *kk_499 = buffer.data(kk + 499);
    const auto *kk_500 = buffer.data(kk + 500);
    const auto *kk_501 = buffer.data(kk + 501);
    const auto *kk_502 = buffer.data(kk + 502);
    const auto *kk_503 = buffer.data(kk + 503);
    const auto *kk_504 = buffer.data(kk + 504);
    const auto *kk_505 = buffer.data(kk + 505);
    const auto *kk_506 = buffer.data(kk + 506);
    const auto *kk_507 = buffer.data(kk + 507);
    const auto *kk_508 = buffer.data(kk + 508);
    const auto *kk_509 = buffer.data(kk + 509);
    const auto *kk_510 = buffer.data(kk + 510);
    const auto *kk_511 = buffer.data(kk + 511);
    const auto *kk_512 = buffer.data(kk + 512);
    const auto *kk_513 = buffer.data(kk + 513);
    const auto *kk_514 = buffer.data(kk + 514);
    const auto *kk_515 = buffer.data(kk + 515);
    const auto *kk_516 = buffer.data(kk + 516);
    const auto *kk_517 = buffer.data(kk + 517);
    const auto *kk_518 = buffer.data(kk + 518);
    const auto *kk_519 = buffer.data(kk + 519);
    const auto *kk_520 = buffer.data(kk + 520);
    const auto *kk_521 = buffer.data(kk + 521);
    const auto *kk_522 = buffer.data(kk + 522);
    const auto *kk_523 = buffer.data(kk + 523);
    const auto *kk_524 = buffer.data(kk + 524);
    const auto *kk_525 = buffer.data(kk + 525);
    const auto *kk_526 = buffer.data(kk + 526);
    const auto *kk_527 = buffer.data(kk + 527);
    const auto *kk_528 = buffer.data(kk + 528);
    const auto *kk_529 = buffer.data(kk + 529);
    const auto *kk_530 = buffer.data(kk + 530);
    const auto *kk_531 = buffer.data(kk + 531);
    const auto *kk_532 = buffer.data(kk + 532);
    const auto *kk_533 = buffer.data(kk + 533);
    const auto *kk_534 = buffer.data(kk + 534);
    const auto *kk_535 = buffer.data(kk + 535);
    const auto *kk_536 = buffer.data(kk + 536);
    const auto *kk_537 = buffer.data(kk + 537);
    const auto *kk_538 = buffer.data(kk + 538);
    const auto *kk_539 = buffer.data(kk + 539);
    const auto *kk_540 = buffer.data(kk + 540);
    const auto *kk_541 = buffer.data(kk + 541);
    const auto *kk_542 = buffer.data(kk + 542);
    const auto *kk_543 = buffer.data(kk + 543);
    const auto *kk_544 = buffer.data(kk + 544);
    const auto *kk_545 = buffer.data(kk + 545);
    const auto *kk_546 = buffer.data(kk + 546);
    const auto *kk_547 = buffer.data(kk + 547);
    const auto *kk_548 = buffer.data(kk + 548);
    const auto *kk_549 = buffer.data(kk + 549);
    const auto *kk_550 = buffer.data(kk + 550);
    const auto *kk_551 = buffer.data(kk + 551);
    const auto *kk_552 = buffer.data(kk + 552);
    const auto *kk_553 = buffer.data(kk + 553);
    const auto *kk_554 = buffer.data(kk + 554);
    const auto *kk_555 = buffer.data(kk + 555);
    const auto *kk_556 = buffer.data(kk + 556);
    const auto *kk_557 = buffer.data(kk + 557);
    const auto *kk_558 = buffer.data(kk + 558);
    const auto *kk_559 = buffer.data(kk + 559);
    const auto *kk_560 = buffer.data(kk + 560);
    const auto *kk_561 = buffer.data(kk + 561);
    const auto *kk_562 = buffer.data(kk + 562);
    const auto *kk_563 = buffer.data(kk + 563);
    const auto *kk_564 = buffer.data(kk + 564);
    const auto *kk_565 = buffer.data(kk + 565);
    const auto *kk_566 = buffer.data(kk + 566);
    const auto *kk_567 = buffer.data(kk + 567);
    const auto *kk_568 = buffer.data(kk + 568);
    const auto *kk_569 = buffer.data(kk + 569);
    const auto *kk_570 = buffer.data(kk + 570);
    const auto *kk_571 = buffer.data(kk + 571);
    const auto *kk_572 = buffer.data(kk + 572);
    const auto *kk_573 = buffer.data(kk + 573);
    const auto *kk_574 = buffer.data(kk + 574);
    const auto *kk_575 = buffer.data(kk + 575);
    const auto *kk_576 = buffer.data(kk + 576);
    const auto *kk_577 = buffer.data(kk + 577);
    const auto *kk_578 = buffer.data(kk + 578);
    const auto *kk_579 = buffer.data(kk + 579);
    const auto *kk_580 = buffer.data(kk + 580);
    const auto *kk_581 = buffer.data(kk + 581);
    const auto *kk_582 = buffer.data(kk + 582);
    const auto *kk_583 = buffer.data(kk + 583);
    const auto *kk_584 = buffer.data(kk + 584);
    const auto *kk_585 = buffer.data(kk + 585);
    const auto *kk_586 = buffer.data(kk + 586);
    const auto *kk_587 = buffer.data(kk + 587);
    const auto *kk_588 = buffer.data(kk + 588);
    const auto *kk_589 = buffer.data(kk + 589);
    const auto *kk_590 = buffer.data(kk + 590);
    const auto *kk_591 = buffer.data(kk + 591);
    const auto *kk_592 = buffer.data(kk + 592);
    const auto *kk_593 = buffer.data(kk + 593);
    const auto *kk_594 = buffer.data(kk + 594);
    const auto *kk_595 = buffer.data(kk + 595);
    const auto *kk_596 = buffer.data(kk + 596);
    const auto *kk_597 = buffer.data(kk + 597);
    const auto *kk_598 = buffer.data(kk + 598);
    const auto *kk_599 = buffer.data(kk + 599);

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, hk_450, hk_451, hk_452, hk_453, \
                         hk_454, kk_450, kk_451, kk_452, kk_453, \
                         kk_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -2.0 * hk_450[k]
                   + f_0 * kk_450[k];

        t_451[k] = -2.0 * hk_451[k]
                   + f_0 * kk_451[k];

        t_452[k] = -2.0 * hk_452[k]
                   + f_0 * kk_452[k];

        t_453[k] = -2.0 * hk_453[k]
                   + f_0 * kk_453[k];

        t_454[k] = -2.0 * hk_454[k]
                   + f_0 * kk_454[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, hk_455, hk_456, hk_457, hk_458, \
                         hk_459, kk_455, kk_456, kk_457, kk_458, \
                         kk_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -2.0 * hk_455[k]
                   + f_0 * kk_455[k];

        t_456[k] = -2.0 * hk_456[k]
                   + f_0 * kk_456[k];

        t_457[k] = -2.0 * hk_457[k]
                   + f_0 * kk_457[k];

        t_458[k] = -2.0 * hk_458[k]
                   + f_0 * kk_458[k];

        t_459[k] = -2.0 * hk_459[k]
                   + f_0 * kk_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, hk_460, hk_461, hk_462, hk_463, \
                         hk_464, kk_460, kk_461, kk_462, kk_463, \
                         kk_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -2.0 * hk_460[k]
                   + f_0 * kk_460[k];

        t_461[k] = -2.0 * hk_461[k]
                   + f_0 * kk_461[k];

        t_462[k] = -2.0 * hk_462[k]
                   + f_0 * kk_462[k];

        t_463[k] = -2.0 * hk_463[k]
                   + f_0 * kk_463[k];

        t_464[k] = -2.0 * hk_464[k]
                   + f_0 * kk_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, hk_465, hk_466, hk_467, hk_468, \
                         hk_469, kk_465, kk_466, kk_467, kk_468, \
                         kk_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -2.0 * hk_465[k]
                   + f_0 * kk_465[k];

        t_466[k] = -2.0 * hk_466[k]
                   + f_0 * kk_466[k];

        t_467[k] = -2.0 * hk_467[k]
                   + f_0 * kk_467[k];

        t_468[k] = -2.0 * hk_468[k]
                   + f_0 * kk_468[k];

        t_469[k] = -2.0 * hk_469[k]
                   + f_0 * kk_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, hk_470, hk_471, hk_472, hk_473, \
                         hk_474, kk_470, kk_471, kk_472, kk_473, \
                         kk_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -2.0 * hk_470[k]
                   + f_0 * kk_470[k];

        t_471[k] = -2.0 * hk_471[k]
                   + f_0 * kk_471[k];

        t_472[k] = -2.0 * hk_472[k]
                   + f_0 * kk_472[k];

        t_473[k] = -2.0 * hk_473[k]
                   + f_0 * kk_473[k];

        t_474[k] = -2.0 * hk_474[k]
                   + f_0 * kk_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, hk_475, hk_476, hk_477, hk_478, \
                         hk_479, kk_475, kk_476, kk_477, kk_478, \
                         kk_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -2.0 * hk_475[k]
                   + f_0 * kk_475[k];

        t_476[k] = -2.0 * hk_476[k]
                   + f_0 * kk_476[k];

        t_477[k] = -2.0 * hk_477[k]
                   + f_0 * kk_477[k];

        t_478[k] = -2.0 * hk_478[k]
                   + f_0 * kk_478[k];

        t_479[k] = -2.0 * hk_479[k]
                   + f_0 * kk_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, hk_480, hk_481, hk_482, hk_483, \
                         hk_484, kk_480, kk_481, kk_482, kk_483, \
                         kk_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -2.0 * hk_480[k]
                   + f_0 * kk_480[k];

        t_481[k] = -2.0 * hk_481[k]
                   + f_0 * kk_481[k];

        t_482[k] = -2.0 * hk_482[k]
                   + f_0 * kk_482[k];

        t_483[k] = -2.0 * hk_483[k]
                   + f_0 * kk_483[k];

        t_484[k] = -2.0 * hk_484[k]
                   + f_0 * kk_484[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, hk_485, hk_486, hk_487, hk_488, \
                         hk_489, kk_485, kk_486, kk_487, kk_488, \
                         kk_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -2.0 * hk_485[k]
                   + f_0 * kk_485[k];

        t_486[k] = -2.0 * hk_486[k]
                   + f_0 * kk_486[k];

        t_487[k] = -2.0 * hk_487[k]
                   + f_0 * kk_487[k];

        t_488[k] = -2.0 * hk_488[k]
                   + f_0 * kk_488[k];

        t_489[k] = -2.0 * hk_489[k]
                   + f_0 * kk_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, hk_490, hk_491, hk_492, hk_493, \
                         hk_494, kk_490, kk_491, kk_492, kk_493, \
                         kk_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -2.0 * hk_490[k]
                   + f_0 * kk_490[k];

        t_491[k] = -2.0 * hk_491[k]
                   + f_0 * kk_491[k];

        t_492[k] = -2.0 * hk_492[k]
                   + f_0 * kk_492[k];

        t_493[k] = -2.0 * hk_493[k]
                   + f_0 * kk_493[k];

        t_494[k] = -2.0 * hk_494[k]
                   + f_0 * kk_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, hk_495, hk_496, hk_497, hk_498, \
                         hk_499, kk_495, kk_496, kk_497, kk_498, \
                         kk_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -2.0 * hk_495[k]
                   + f_0 * kk_495[k];

        t_496[k] = -2.0 * hk_496[k]
                   + f_0 * kk_496[k];

        t_497[k] = -2.0 * hk_497[k]
                   + f_0 * kk_497[k];

        t_498[k] = -2.0 * hk_498[k]
                   + f_0 * kk_498[k];

        t_499[k] = -2.0 * hk_499[k]
                   + f_0 * kk_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, hk_500, hk_501, hk_502, hk_503, \
                         hk_504, kk_500, kk_501, kk_502, kk_503, \
                         kk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -2.0 * hk_500[k]
                   + f_0 * kk_500[k];

        t_501[k] = -2.0 * hk_501[k]
                   + f_0 * kk_501[k];

        t_502[k] = -2.0 * hk_502[k]
                   + f_0 * kk_502[k];

        t_503[k] = -2.0 * hk_503[k]
                   + f_0 * kk_503[k];

        t_504[k] = -2.0 * hk_504[k]
                   + f_0 * kk_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, hk_505, hk_506, hk_507, hk_508, \
                         hk_509, kk_505, kk_506, kk_507, kk_508, \
                         kk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -2.0 * hk_505[k]
                   + f_0 * kk_505[k];

        t_506[k] = -2.0 * hk_506[k]
                   + f_0 * kk_506[k];

        t_507[k] = -2.0 * hk_507[k]
                   + f_0 * kk_507[k];

        t_508[k] = -2.0 * hk_508[k]
                   + f_0 * kk_508[k];

        t_509[k] = -2.0 * hk_509[k]
                   + f_0 * kk_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, hk_510, hk_511, hk_512, hk_513, \
                         hk_514, kk_510, kk_511, kk_512, kk_513, \
                         kk_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -2.0 * hk_510[k]
                   + f_0 * kk_510[k];

        t_511[k] = -2.0 * hk_511[k]
                   + f_0 * kk_511[k];

        t_512[k] = -2.0 * hk_512[k]
                   + f_0 * kk_512[k];

        t_513[k] = -2.0 * hk_513[k]
                   + f_0 * kk_513[k];

        t_514[k] = -2.0 * hk_514[k]
                   + f_0 * kk_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, hk_515, hk_516, hk_517, hk_518, \
                         hk_519, kk_515, kk_516, kk_517, kk_518, \
                         kk_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -2.0 * hk_515[k]
                   + f_0 * kk_515[k];

        t_516[k] = -2.0 * hk_516[k]
                   + f_0 * kk_516[k];

        t_517[k] = -2.0 * hk_517[k]
                   + f_0 * kk_517[k];

        t_518[k] = -2.0 * hk_518[k]
                   + f_0 * kk_518[k];

        t_519[k] = -2.0 * hk_519[k]
                   + f_0 * kk_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, hk_520, hk_521, hk_522, hk_523, \
                         hk_524, kk_520, kk_521, kk_522, kk_523, \
                         kk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -2.0 * hk_520[k]
                   + f_0 * kk_520[k];

        t_521[k] = -2.0 * hk_521[k]
                   + f_0 * kk_521[k];

        t_522[k] = -2.0 * hk_522[k]
                   + f_0 * kk_522[k];

        t_523[k] = -2.0 * hk_523[k]
                   + f_0 * kk_523[k];

        t_524[k] = -2.0 * hk_524[k]
                   + f_0 * kk_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, hk_525, hk_526, hk_527, hk_528, \
                         hk_529, kk_525, kk_526, kk_527, kk_528, \
                         kk_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -2.0 * hk_525[k]
                   + f_0 * kk_525[k];

        t_526[k] = -2.0 * hk_526[k]
                   + f_0 * kk_526[k];

        t_527[k] = -2.0 * hk_527[k]
                   + f_0 * kk_527[k];

        t_528[k] = -2.0 * hk_528[k]
                   + f_0 * kk_528[k];

        t_529[k] = -2.0 * hk_529[k]
                   + f_0 * kk_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, hk_530, hk_531, hk_532, hk_533, \
                         hk_534, kk_530, kk_531, kk_532, kk_533, \
                         kk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -2.0 * hk_530[k]
                   + f_0 * kk_530[k];

        t_531[k] = -2.0 * hk_531[k]
                   + f_0 * kk_531[k];

        t_532[k] = -2.0 * hk_532[k]
                   + f_0 * kk_532[k];

        t_533[k] = -2.0 * hk_533[k]
                   + f_0 * kk_533[k];

        t_534[k] = -2.0 * hk_534[k]
                   + f_0 * kk_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, hk_535, hk_536, hk_537, hk_538, \
                         hk_539, kk_535, kk_536, kk_537, kk_538, \
                         kk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -2.0 * hk_535[k]
                   + f_0 * kk_535[k];

        t_536[k] = -2.0 * hk_536[k]
                   + f_0 * kk_536[k];

        t_537[k] = -2.0 * hk_537[k]
                   + f_0 * kk_537[k];

        t_538[k] = -2.0 * hk_538[k]
                   + f_0 * kk_538[k];

        t_539[k] = -2.0 * hk_539[k]
                   + f_0 * kk_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, hk_540, hk_541, hk_542, hk_543, \
                         hk_544, kk_540, kk_541, kk_542, kk_543, \
                         kk_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -hk_540[k]
                   + f_0 * kk_540[k];

        t_541[k] = -hk_541[k]
                   + f_0 * kk_541[k];

        t_542[k] = -hk_542[k]
                   + f_0 * kk_542[k];

        t_543[k] = -hk_543[k]
                   + f_0 * kk_543[k];

        t_544[k] = -hk_544[k]
                   + f_0 * kk_544[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, hk_545, hk_546, hk_547, hk_548, \
                         hk_549, kk_545, kk_546, kk_547, kk_548, \
                         kk_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -hk_545[k]
                   + f_0 * kk_545[k];

        t_546[k] = -hk_546[k]
                   + f_0 * kk_546[k];

        t_547[k] = -hk_547[k]
                   + f_0 * kk_547[k];

        t_548[k] = -hk_548[k]
                   + f_0 * kk_548[k];

        t_549[k] = -hk_549[k]
                   + f_0 * kk_549[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, hk_550, hk_551, hk_552, hk_553, \
                         hk_554, kk_550, kk_551, kk_552, kk_553, \
                         kk_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -hk_550[k]
                   + f_0 * kk_550[k];

        t_551[k] = -hk_551[k]
                   + f_0 * kk_551[k];

        t_552[k] = -hk_552[k]
                   + f_0 * kk_552[k];

        t_553[k] = -hk_553[k]
                   + f_0 * kk_553[k];

        t_554[k] = -hk_554[k]
                   + f_0 * kk_554[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, hk_555, hk_556, hk_557, hk_558, \
                         hk_559, kk_555, kk_556, kk_557, kk_558, \
                         kk_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -hk_555[k]
                   + f_0 * kk_555[k];

        t_556[k] = -hk_556[k]
                   + f_0 * kk_556[k];

        t_557[k] = -hk_557[k]
                   + f_0 * kk_557[k];

        t_558[k] = -hk_558[k]
                   + f_0 * kk_558[k];

        t_559[k] = -hk_559[k]
                   + f_0 * kk_559[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, hk_560, hk_561, hk_562, hk_563, \
                         hk_564, kk_560, kk_561, kk_562, kk_563, \
                         kk_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = -hk_560[k]
                   + f_0 * kk_560[k];

        t_561[k] = -hk_561[k]
                   + f_0 * kk_561[k];

        t_562[k] = -hk_562[k]
                   + f_0 * kk_562[k];

        t_563[k] = -hk_563[k]
                   + f_0 * kk_563[k];

        t_564[k] = -hk_564[k]
                   + f_0 * kk_564[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, hk_565, hk_566, hk_567, hk_568, \
                         hk_569, kk_565, kk_566, kk_567, kk_568, \
                         kk_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = -hk_565[k]
                   + f_0 * kk_565[k];

        t_566[k] = -hk_566[k]
                   + f_0 * kk_566[k];

        t_567[k] = -hk_567[k]
                   + f_0 * kk_567[k];

        t_568[k] = -hk_568[k]
                   + f_0 * kk_568[k];

        t_569[k] = -hk_569[k]
                   + f_0 * kk_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, hk_570, hk_571, hk_572, hk_573, \
                         hk_574, kk_570, kk_571, kk_572, kk_573, \
                         kk_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -hk_570[k]
                   + f_0 * kk_570[k];

        t_571[k] = -hk_571[k]
                   + f_0 * kk_571[k];

        t_572[k] = -hk_572[k]
                   + f_0 * kk_572[k];

        t_573[k] = -hk_573[k]
                   + f_0 * kk_573[k];

        t_574[k] = -hk_574[k]
                   + f_0 * kk_574[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, hk_575, hk_576, hk_577, hk_578, \
                         hk_579, kk_575, kk_576, kk_577, kk_578, \
                         kk_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = -hk_575[k]
                   + f_0 * kk_575[k];

        t_576[k] = -hk_576[k]
                   + f_0 * kk_576[k];

        t_577[k] = -hk_577[k]
                   + f_0 * kk_577[k];

        t_578[k] = -hk_578[k]
                   + f_0 * kk_578[k];

        t_579[k] = -hk_579[k]
                   + f_0 * kk_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, hk_580, hk_581, hk_582, hk_583, \
                         hk_584, kk_580, kk_581, kk_582, kk_583, \
                         kk_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -hk_580[k]
                   + f_0 * kk_580[k];

        t_581[k] = -hk_581[k]
                   + f_0 * kk_581[k];

        t_582[k] = -hk_582[k]
                   + f_0 * kk_582[k];

        t_583[k] = -hk_583[k]
                   + f_0 * kk_583[k];

        t_584[k] = -hk_584[k]
                   + f_0 * kk_584[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, hk_585, hk_586, hk_587, hk_588, \
                         hk_589, kk_585, kk_586, kk_587, kk_588, \
                         kk_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -hk_585[k]
                   + f_0 * kk_585[k];

        t_586[k] = -hk_586[k]
                   + f_0 * kk_586[k];

        t_587[k] = -hk_587[k]
                   + f_0 * kk_587[k];

        t_588[k] = -hk_588[k]
                   + f_0 * kk_588[k];

        t_589[k] = -hk_589[k]
                   + f_0 * kk_589[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, hk_590, hk_591, hk_592, hk_593, \
                         hk_594, kk_590, kk_591, kk_592, kk_593, \
                         kk_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -hk_590[k]
                   + f_0 * kk_590[k];

        t_591[k] = -hk_591[k]
                   + f_0 * kk_591[k];

        t_592[k] = -hk_592[k]
                   + f_0 * kk_592[k];

        t_593[k] = -hk_593[k]
                   + f_0 * kk_593[k];

        t_594[k] = -hk_594[k]
                   + f_0 * kk_594[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, hk_595, hk_596, hk_597, hk_598, \
                         hk_599, kk_595, kk_596, kk_597, kk_598, \
                         kk_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -hk_595[k]
                   + f_0 * kk_595[k];

        t_596[k] = -hk_596[k]
                   + f_0 * kk_596[k];

        t_597[k] = -hk_597[k]
                   + f_0 * kk_597[k];

        t_598[k] = -hk_598[k]
                   + f_0 * kk_598[k];

        t_599[k] = -hk_599[k]
                   + f_0 * kk_599[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_600 = buffer.data(hk + 600);
    const auto *hk_601 = buffer.data(hk + 601);
    const auto *hk_602 = buffer.data(hk + 602);
    const auto *hk_603 = buffer.data(hk + 603);
    const auto *hk_604 = buffer.data(hk + 604);
    const auto *hk_605 = buffer.data(hk + 605);
    const auto *hk_606 = buffer.data(hk + 606);
    const auto *hk_607 = buffer.data(hk + 607);
    const auto *hk_608 = buffer.data(hk + 608);
    const auto *hk_609 = buffer.data(hk + 609);
    const auto *hk_610 = buffer.data(hk + 610);
    const auto *hk_611 = buffer.data(hk + 611);
    const auto *hk_612 = buffer.data(hk + 612);
    const auto *hk_613 = buffer.data(hk + 613);
    const auto *hk_614 = buffer.data(hk + 614);
    const auto *hk_615 = buffer.data(hk + 615);
    const auto *hk_616 = buffer.data(hk + 616);
    const auto *hk_617 = buffer.data(hk + 617);
    const auto *hk_618 = buffer.data(hk + 618);
    const auto *hk_619 = buffer.data(hk + 619);
    const auto *hk_620 = buffer.data(hk + 620);
    const auto *hk_621 = buffer.data(hk + 621);
    const auto *hk_622 = buffer.data(hk + 622);
    const auto *hk_623 = buffer.data(hk + 623);
    const auto *hk_624 = buffer.data(hk + 624);
    const auto *hk_625 = buffer.data(hk + 625);
    const auto *hk_626 = buffer.data(hk + 626);
    const auto *hk_627 = buffer.data(hk + 627);
    const auto *hk_628 = buffer.data(hk + 628);
    const auto *hk_629 = buffer.data(hk + 629);
    const auto *hk_630 = buffer.data(hk + 630);
    const auto *hk_631 = buffer.data(hk + 631);
    const auto *hk_632 = buffer.data(hk + 632);
    const auto *hk_633 = buffer.data(hk + 633);
    const auto *hk_634 = buffer.data(hk + 634);
    const auto *hk_635 = buffer.data(hk + 635);
    const auto *hk_636 = buffer.data(hk + 636);
    const auto *hk_637 = buffer.data(hk + 637);
    const auto *hk_638 = buffer.data(hk + 638);
    const auto *hk_639 = buffer.data(hk + 639);
    const auto *hk_640 = buffer.data(hk + 640);
    const auto *hk_641 = buffer.data(hk + 641);
    const auto *hk_642 = buffer.data(hk + 642);
    const auto *hk_643 = buffer.data(hk + 643);
    const auto *hk_644 = buffer.data(hk + 644);
    const auto *hk_645 = buffer.data(hk + 645);
    const auto *hk_646 = buffer.data(hk + 646);
    const auto *hk_647 = buffer.data(hk + 647);
    const auto *hk_648 = buffer.data(hk + 648);
    const auto *hk_649 = buffer.data(hk + 649);
    const auto *hk_650 = buffer.data(hk + 650);
    const auto *hk_651 = buffer.data(hk + 651);
    const auto *hk_652 = buffer.data(hk + 652);
    const auto *hk_653 = buffer.data(hk + 653);
    const auto *hk_654 = buffer.data(hk + 654);
    const auto *hk_655 = buffer.data(hk + 655);
    const auto *hk_656 = buffer.data(hk + 656);
    const auto *hk_657 = buffer.data(hk + 657);
    const auto *hk_658 = buffer.data(hk + 658);
    const auto *hk_659 = buffer.data(hk + 659);
    const auto *hk_660 = buffer.data(hk + 660);
    const auto *hk_661 = buffer.data(hk + 661);
    const auto *hk_662 = buffer.data(hk + 662);
    const auto *hk_663 = buffer.data(hk + 663);
    const auto *hk_664 = buffer.data(hk + 664);
    const auto *hk_665 = buffer.data(hk + 665);
    const auto *hk_666 = buffer.data(hk + 666);
    const auto *hk_667 = buffer.data(hk + 667);
    const auto *hk_668 = buffer.data(hk + 668);
    const auto *hk_669 = buffer.data(hk + 669);
    const auto *hk_670 = buffer.data(hk + 670);
    const auto *hk_671 = buffer.data(hk + 671);
    const auto *hk_672 = buffer.data(hk + 672);
    const auto *hk_673 = buffer.data(hk + 673);
    const auto *hk_674 = buffer.data(hk + 674);
    const auto *hk_675 = buffer.data(hk + 675);
    const auto *hk_676 = buffer.data(hk + 676);
    const auto *hk_677 = buffer.data(hk + 677);
    const auto *hk_678 = buffer.data(hk + 678);
    const auto *hk_679 = buffer.data(hk + 679);
    const auto *hk_680 = buffer.data(hk + 680);
    const auto *hk_681 = buffer.data(hk + 681);
    const auto *hk_682 = buffer.data(hk + 682);
    const auto *hk_683 = buffer.data(hk + 683);
    const auto *hk_684 = buffer.data(hk + 684);
    const auto *hk_685 = buffer.data(hk + 685);
    const auto *hk_686 = buffer.data(hk + 686);
    const auto *hk_687 = buffer.data(hk + 687);
    const auto *hk_688 = buffer.data(hk + 688);
    const auto *hk_689 = buffer.data(hk + 689);
    const auto *hk_690 = buffer.data(hk + 690);
    const auto *hk_691 = buffer.data(hk + 691);
    const auto *hk_692 = buffer.data(hk + 692);
    const auto *hk_693 = buffer.data(hk + 693);
    const auto *hk_694 = buffer.data(hk + 694);
    const auto *hk_695 = buffer.data(hk + 695);
    const auto *hk_696 = buffer.data(hk + 696);
    const auto *hk_697 = buffer.data(hk + 697);
    const auto *hk_698 = buffer.data(hk + 698);
    const auto *hk_699 = buffer.data(hk + 699);
    const auto *hk_700 = buffer.data(hk + 700);
    const auto *hk_701 = buffer.data(hk + 701);
    const auto *hk_702 = buffer.data(hk + 702);
    const auto *hk_703 = buffer.data(hk + 703);
    const auto *hk_704 = buffer.data(hk + 704);
    const auto *hk_705 = buffer.data(hk + 705);
    const auto *hk_706 = buffer.data(hk + 706);
    const auto *hk_707 = buffer.data(hk + 707);
    const auto *hk_708 = buffer.data(hk + 708);
    const auto *hk_709 = buffer.data(hk + 709);
    const auto *hk_710 = buffer.data(hk + 710);
    const auto *hk_711 = buffer.data(hk + 711);
    const auto *hk_712 = buffer.data(hk + 712);
    const auto *hk_713 = buffer.data(hk + 713);
    const auto *hk_714 = buffer.data(hk + 714);
    const auto *hk_715 = buffer.data(hk + 715);
    const auto *hk_716 = buffer.data(hk + 716);
    const auto *hk_717 = buffer.data(hk + 717);
    const auto *hk_718 = buffer.data(hk + 718);
    const auto *hk_719 = buffer.data(hk + 719);
    const auto *hk_720 = buffer.data(hk + 720);
    const auto *hk_721 = buffer.data(hk + 721);
    const auto *hk_722 = buffer.data(hk + 722);
    const auto *hk_723 = buffer.data(hk + 723);
    const auto *hk_724 = buffer.data(hk + 724);
    const auto *hk_725 = buffer.data(hk + 725);
    const auto *hk_726 = buffer.data(hk + 726);
    const auto *hk_727 = buffer.data(hk + 727);
    const auto *hk_728 = buffer.data(hk + 728);
    const auto *hk_729 = buffer.data(hk + 729);
    const auto *hk_730 = buffer.data(hk + 730);
    const auto *hk_731 = buffer.data(hk + 731);
    const auto *hk_732 = buffer.data(hk + 732);
    const auto *hk_733 = buffer.data(hk + 733);
    const auto *hk_734 = buffer.data(hk + 734);
    const auto *hk_735 = buffer.data(hk + 735);
    const auto *hk_736 = buffer.data(hk + 736);
    const auto *hk_737 = buffer.data(hk + 737);
    const auto *hk_738 = buffer.data(hk + 738);
    const auto *hk_739 = buffer.data(hk + 739);
    const auto *hk_740 = buffer.data(hk + 740);
    const auto *hk_741 = buffer.data(hk + 741);
    const auto *hk_742 = buffer.data(hk + 742);
    const auto *hk_743 = buffer.data(hk + 743);
    const auto *hk_744 = buffer.data(hk + 744);
    const auto *hk_745 = buffer.data(hk + 745);
    const auto *hk_746 = buffer.data(hk + 746);
    const auto *hk_747 = buffer.data(hk + 747);
    const auto *hk_748 = buffer.data(hk + 748);
    const auto *hk_749 = buffer.data(hk + 749);

    const auto *kk_600 = buffer.data(kk + 600);
    const auto *kk_601 = buffer.data(kk + 601);
    const auto *kk_602 = buffer.data(kk + 602);
    const auto *kk_603 = buffer.data(kk + 603);
    const auto *kk_604 = buffer.data(kk + 604);
    const auto *kk_605 = buffer.data(kk + 605);
    const auto *kk_606 = buffer.data(kk + 606);
    const auto *kk_607 = buffer.data(kk + 607);
    const auto *kk_608 = buffer.data(kk + 608);
    const auto *kk_609 = buffer.data(kk + 609);
    const auto *kk_610 = buffer.data(kk + 610);
    const auto *kk_611 = buffer.data(kk + 611);
    const auto *kk_612 = buffer.data(kk + 612);
    const auto *kk_613 = buffer.data(kk + 613);
    const auto *kk_614 = buffer.data(kk + 614);
    const auto *kk_615 = buffer.data(kk + 615);
    const auto *kk_616 = buffer.data(kk + 616);
    const auto *kk_617 = buffer.data(kk + 617);
    const auto *kk_618 = buffer.data(kk + 618);
    const auto *kk_619 = buffer.data(kk + 619);
    const auto *kk_620 = buffer.data(kk + 620);
    const auto *kk_621 = buffer.data(kk + 621);
    const auto *kk_622 = buffer.data(kk + 622);
    const auto *kk_623 = buffer.data(kk + 623);
    const auto *kk_624 = buffer.data(kk + 624);
    const auto *kk_625 = buffer.data(kk + 625);
    const auto *kk_626 = buffer.data(kk + 626);
    const auto *kk_627 = buffer.data(kk + 627);
    const auto *kk_628 = buffer.data(kk + 628);
    const auto *kk_629 = buffer.data(kk + 629);
    const auto *kk_630 = buffer.data(kk + 630);
    const auto *kk_631 = buffer.data(kk + 631);
    const auto *kk_632 = buffer.data(kk + 632);
    const auto *kk_633 = buffer.data(kk + 633);
    const auto *kk_634 = buffer.data(kk + 634);
    const auto *kk_635 = buffer.data(kk + 635);
    const auto *kk_636 = buffer.data(kk + 636);
    const auto *kk_637 = buffer.data(kk + 637);
    const auto *kk_638 = buffer.data(kk + 638);
    const auto *kk_639 = buffer.data(kk + 639);
    const auto *kk_640 = buffer.data(kk + 640);
    const auto *kk_641 = buffer.data(kk + 641);
    const auto *kk_642 = buffer.data(kk + 642);
    const auto *kk_643 = buffer.data(kk + 643);
    const auto *kk_644 = buffer.data(kk + 644);
    const auto *kk_645 = buffer.data(kk + 645);
    const auto *kk_646 = buffer.data(kk + 646);
    const auto *kk_647 = buffer.data(kk + 647);
    const auto *kk_648 = buffer.data(kk + 648);
    const auto *kk_649 = buffer.data(kk + 649);
    const auto *kk_650 = buffer.data(kk + 650);
    const auto *kk_651 = buffer.data(kk + 651);
    const auto *kk_652 = buffer.data(kk + 652);
    const auto *kk_653 = buffer.data(kk + 653);
    const auto *kk_654 = buffer.data(kk + 654);
    const auto *kk_655 = buffer.data(kk + 655);
    const auto *kk_656 = buffer.data(kk + 656);
    const auto *kk_657 = buffer.data(kk + 657);
    const auto *kk_658 = buffer.data(kk + 658);
    const auto *kk_659 = buffer.data(kk + 659);
    const auto *kk_660 = buffer.data(kk + 660);
    const auto *kk_661 = buffer.data(kk + 661);
    const auto *kk_662 = buffer.data(kk + 662);
    const auto *kk_663 = buffer.data(kk + 663);
    const auto *kk_664 = buffer.data(kk + 664);
    const auto *kk_665 = buffer.data(kk + 665);
    const auto *kk_666 = buffer.data(kk + 666);
    const auto *kk_667 = buffer.data(kk + 667);
    const auto *kk_668 = buffer.data(kk + 668);
    const auto *kk_669 = buffer.data(kk + 669);
    const auto *kk_670 = buffer.data(kk + 670);
    const auto *kk_671 = buffer.data(kk + 671);
    const auto *kk_672 = buffer.data(kk + 672);
    const auto *kk_673 = buffer.data(kk + 673);
    const auto *kk_674 = buffer.data(kk + 674);
    const auto *kk_675 = buffer.data(kk + 675);
    const auto *kk_676 = buffer.data(kk + 676);
    const auto *kk_677 = buffer.data(kk + 677);
    const auto *kk_678 = buffer.data(kk + 678);
    const auto *kk_679 = buffer.data(kk + 679);
    const auto *kk_680 = buffer.data(kk + 680);
    const auto *kk_681 = buffer.data(kk + 681);
    const auto *kk_682 = buffer.data(kk + 682);
    const auto *kk_683 = buffer.data(kk + 683);
    const auto *kk_684 = buffer.data(kk + 684);
    const auto *kk_685 = buffer.data(kk + 685);
    const auto *kk_686 = buffer.data(kk + 686);
    const auto *kk_687 = buffer.data(kk + 687);
    const auto *kk_688 = buffer.data(kk + 688);
    const auto *kk_689 = buffer.data(kk + 689);
    const auto *kk_690 = buffer.data(kk + 690);
    const auto *kk_691 = buffer.data(kk + 691);
    const auto *kk_692 = buffer.data(kk + 692);
    const auto *kk_693 = buffer.data(kk + 693);
    const auto *kk_694 = buffer.data(kk + 694);
    const auto *kk_695 = buffer.data(kk + 695);
    const auto *kk_696 = buffer.data(kk + 696);
    const auto *kk_697 = buffer.data(kk + 697);
    const auto *kk_698 = buffer.data(kk + 698);
    const auto *kk_699 = buffer.data(kk + 699);
    const auto *kk_700 = buffer.data(kk + 700);
    const auto *kk_701 = buffer.data(kk + 701);
    const auto *kk_702 = buffer.data(kk + 702);
    const auto *kk_703 = buffer.data(kk + 703);
    const auto *kk_704 = buffer.data(kk + 704);
    const auto *kk_705 = buffer.data(kk + 705);
    const auto *kk_706 = buffer.data(kk + 706);
    const auto *kk_707 = buffer.data(kk + 707);
    const auto *kk_708 = buffer.data(kk + 708);
    const auto *kk_709 = buffer.data(kk + 709);
    const auto *kk_710 = buffer.data(kk + 710);
    const auto *kk_711 = buffer.data(kk + 711);
    const auto *kk_712 = buffer.data(kk + 712);
    const auto *kk_713 = buffer.data(kk + 713);
    const auto *kk_714 = buffer.data(kk + 714);
    const auto *kk_715 = buffer.data(kk + 715);
    const auto *kk_716 = buffer.data(kk + 716);
    const auto *kk_717 = buffer.data(kk + 717);
    const auto *kk_718 = buffer.data(kk + 718);
    const auto *kk_719 = buffer.data(kk + 719);
    const auto *kk_720 = buffer.data(kk + 720);
    const auto *kk_721 = buffer.data(kk + 721);
    const auto *kk_722 = buffer.data(kk + 722);
    const auto *kk_723 = buffer.data(kk + 723);
    const auto *kk_724 = buffer.data(kk + 724);
    const auto *kk_725 = buffer.data(kk + 725);
    const auto *kk_726 = buffer.data(kk + 726);
    const auto *kk_727 = buffer.data(kk + 727);
    const auto *kk_728 = buffer.data(kk + 728);
    const auto *kk_729 = buffer.data(kk + 729);
    const auto *kk_730 = buffer.data(kk + 730);
    const auto *kk_731 = buffer.data(kk + 731);
    const auto *kk_732 = buffer.data(kk + 732);
    const auto *kk_733 = buffer.data(kk + 733);
    const auto *kk_734 = buffer.data(kk + 734);
    const auto *kk_735 = buffer.data(kk + 735);
    const auto *kk_736 = buffer.data(kk + 736);
    const auto *kk_737 = buffer.data(kk + 737);
    const auto *kk_738 = buffer.data(kk + 738);
    const auto *kk_739 = buffer.data(kk + 739);
    const auto *kk_740 = buffer.data(kk + 740);
    const auto *kk_741 = buffer.data(kk + 741);
    const auto *kk_742 = buffer.data(kk + 742);
    const auto *kk_743 = buffer.data(kk + 743);
    const auto *kk_744 = buffer.data(kk + 744);
    const auto *kk_745 = buffer.data(kk + 745);
    const auto *kk_746 = buffer.data(kk + 746);
    const auto *kk_747 = buffer.data(kk + 747);
    const auto *kk_748 = buffer.data(kk + 748);
    const auto *kk_749 = buffer.data(kk + 749);

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, hk_600, hk_601, hk_602, hk_603, \
                         hk_604, kk_600, kk_601, kk_602, kk_603, \
                         kk_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -hk_600[k]
                   + f_0 * kk_600[k];

        t_601[k] = -hk_601[k]
                   + f_0 * kk_601[k];

        t_602[k] = -hk_602[k]
                   + f_0 * kk_602[k];

        t_603[k] = -hk_603[k]
                   + f_0 * kk_603[k];

        t_604[k] = -hk_604[k]
                   + f_0 * kk_604[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, hk_605, hk_606, hk_607, hk_608, \
                         hk_609, kk_605, kk_606, kk_607, kk_608, \
                         kk_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -hk_605[k]
                   + f_0 * kk_605[k];

        t_606[k] = -hk_606[k]
                   + f_0 * kk_606[k];

        t_607[k] = -hk_607[k]
                   + f_0 * kk_607[k];

        t_608[k] = -hk_608[k]
                   + f_0 * kk_608[k];

        t_609[k] = -hk_609[k]
                   + f_0 * kk_609[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, hk_610, hk_611, hk_612, hk_613, \
                         hk_614, kk_610, kk_611, kk_612, kk_613, \
                         kk_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -hk_610[k]
                   + f_0 * kk_610[k];

        t_611[k] = -hk_611[k]
                   + f_0 * kk_611[k];

        t_612[k] = -hk_612[k]
                   + f_0 * kk_612[k];

        t_613[k] = -hk_613[k]
                   + f_0 * kk_613[k];

        t_614[k] = -hk_614[k]
                   + f_0 * kk_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, hk_615, hk_616, hk_617, hk_618, \
                         hk_619, kk_615, kk_616, kk_617, kk_618, \
                         kk_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -hk_615[k]
                   + f_0 * kk_615[k];

        t_616[k] = -hk_616[k]
                   + f_0 * kk_616[k];

        t_617[k] = -hk_617[k]
                   + f_0 * kk_617[k];

        t_618[k] = -hk_618[k]
                   + f_0 * kk_618[k];

        t_619[k] = -hk_619[k]
                   + f_0 * kk_619[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, hk_620, hk_621, hk_622, hk_623, \
                         hk_624, kk_620, kk_621, kk_622, kk_623, \
                         kk_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -hk_620[k]
                   + f_0 * kk_620[k];

        t_621[k] = -hk_621[k]
                   + f_0 * kk_621[k];

        t_622[k] = -hk_622[k]
                   + f_0 * kk_622[k];

        t_623[k] = -hk_623[k]
                   + f_0 * kk_623[k];

        t_624[k] = -hk_624[k]
                   + f_0 * kk_624[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, hk_625, hk_626, hk_627, hk_628, \
                         hk_629, kk_625, kk_626, kk_627, kk_628, \
                         kk_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -hk_625[k]
                   + f_0 * kk_625[k];

        t_626[k] = -hk_626[k]
                   + f_0 * kk_626[k];

        t_627[k] = -hk_627[k]
                   + f_0 * kk_627[k];

        t_628[k] = -hk_628[k]
                   + f_0 * kk_628[k];

        t_629[k] = -hk_629[k]
                   + f_0 * kk_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, hk_630, hk_631, hk_632, hk_633, \
                         hk_634, kk_630, kk_631, kk_632, kk_633, \
                         kk_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -hk_630[k]
                   + f_0 * kk_630[k];

        t_631[k] = -hk_631[k]
                   + f_0 * kk_631[k];

        t_632[k] = -hk_632[k]
                   + f_0 * kk_632[k];

        t_633[k] = -hk_633[k]
                   + f_0 * kk_633[k];

        t_634[k] = -hk_634[k]
                   + f_0 * kk_634[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, hk_635, hk_636, hk_637, hk_638, \
                         hk_639, kk_635, kk_636, kk_637, kk_638, \
                         kk_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -hk_635[k]
                   + f_0 * kk_635[k];

        t_636[k] = -hk_636[k]
                   + f_0 * kk_636[k];

        t_637[k] = -hk_637[k]
                   + f_0 * kk_637[k];

        t_638[k] = -hk_638[k]
                   + f_0 * kk_638[k];

        t_639[k] = -hk_639[k]
                   + f_0 * kk_639[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, hk_640, hk_641, hk_642, hk_643, \
                         hk_644, kk_640, kk_641, kk_642, kk_643, \
                         kk_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -hk_640[k]
                   + f_0 * kk_640[k];

        t_641[k] = -hk_641[k]
                   + f_0 * kk_641[k];

        t_642[k] = -hk_642[k]
                   + f_0 * kk_642[k];

        t_643[k] = -hk_643[k]
                   + f_0 * kk_643[k];

        t_644[k] = -hk_644[k]
                   + f_0 * kk_644[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, hk_645, hk_646, hk_647, hk_648, \
                         hk_649, kk_645, kk_646, kk_647, kk_648, \
                         kk_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -hk_645[k]
                   + f_0 * kk_645[k];

        t_646[k] = -hk_646[k]
                   + f_0 * kk_646[k];

        t_647[k] = -hk_647[k]
                   + f_0 * kk_647[k];

        t_648[k] = -hk_648[k]
                   + f_0 * kk_648[k];

        t_649[k] = -hk_649[k]
                   + f_0 * kk_649[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, hk_650, hk_651, hk_652, hk_653, \
                         hk_654, kk_650, kk_651, kk_652, kk_653, \
                         kk_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -hk_650[k]
                   + f_0 * kk_650[k];

        t_651[k] = -hk_651[k]
                   + f_0 * kk_651[k];

        t_652[k] = -hk_652[k]
                   + f_0 * kk_652[k];

        t_653[k] = -hk_653[k]
                   + f_0 * kk_653[k];

        t_654[k] = -hk_654[k]
                   + f_0 * kk_654[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, hk_655, hk_656, hk_657, hk_658, \
                         hk_659, kk_655, kk_656, kk_657, kk_658, \
                         kk_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -hk_655[k]
                   + f_0 * kk_655[k];

        t_656[k] = -hk_656[k]
                   + f_0 * kk_656[k];

        t_657[k] = -hk_657[k]
                   + f_0 * kk_657[k];

        t_658[k] = -hk_658[k]
                   + f_0 * kk_658[k];

        t_659[k] = -hk_659[k]
                   + f_0 * kk_659[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, hk_660, hk_661, hk_662, hk_663, \
                         hk_664, kk_660, kk_661, kk_662, kk_663, \
                         kk_664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -hk_660[k]
                   + f_0 * kk_660[k];

        t_661[k] = -hk_661[k]
                   + f_0 * kk_661[k];

        t_662[k] = -hk_662[k]
                   + f_0 * kk_662[k];

        t_663[k] = -hk_663[k]
                   + f_0 * kk_663[k];

        t_664[k] = -hk_664[k]
                   + f_0 * kk_664[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, hk_665, hk_666, hk_667, hk_668, \
                         hk_669, kk_665, kk_666, kk_667, kk_668, \
                         kk_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -hk_665[k]
                   + f_0 * kk_665[k];

        t_666[k] = -hk_666[k]
                   + f_0 * kk_666[k];

        t_667[k] = -hk_667[k]
                   + f_0 * kk_667[k];

        t_668[k] = -hk_668[k]
                   + f_0 * kk_668[k];

        t_669[k] = -hk_669[k]
                   + f_0 * kk_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, hk_670, hk_671, hk_672, hk_673, \
                         hk_674, kk_670, kk_671, kk_672, kk_673, \
                         kk_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -hk_670[k]
                   + f_0 * kk_670[k];

        t_671[k] = -hk_671[k]
                   + f_0 * kk_671[k];

        t_672[k] = -hk_672[k]
                   + f_0 * kk_672[k];

        t_673[k] = -hk_673[k]
                   + f_0 * kk_673[k];

        t_674[k] = -hk_674[k]
                   + f_0 * kk_674[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, hk_675, hk_676, hk_677, hk_678, \
                         hk_679, kk_675, kk_676, kk_677, kk_678, \
                         kk_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -hk_675[k]
                   + f_0 * kk_675[k];

        t_676[k] = -hk_676[k]
                   + f_0 * kk_676[k];

        t_677[k] = -hk_677[k]
                   + f_0 * kk_677[k];

        t_678[k] = -hk_678[k]
                   + f_0 * kk_678[k];

        t_679[k] = -hk_679[k]
                   + f_0 * kk_679[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, hk_680, hk_681, hk_682, hk_683, \
                         hk_684, kk_680, kk_681, kk_682, kk_683, \
                         kk_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -hk_680[k]
                   + f_0 * kk_680[k];

        t_681[k] = -hk_681[k]
                   + f_0 * kk_681[k];

        t_682[k] = -hk_682[k]
                   + f_0 * kk_682[k];

        t_683[k] = -hk_683[k]
                   + f_0 * kk_683[k];

        t_684[k] = -hk_684[k]
                   + f_0 * kk_684[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, hk_685, hk_686, hk_687, hk_688, \
                         hk_689, kk_685, kk_686, kk_687, kk_688, \
                         kk_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -hk_685[k]
                   + f_0 * kk_685[k];

        t_686[k] = -hk_686[k]
                   + f_0 * kk_686[k];

        t_687[k] = -hk_687[k]
                   + f_0 * kk_687[k];

        t_688[k] = -hk_688[k]
                   + f_0 * kk_688[k];

        t_689[k] = -hk_689[k]
                   + f_0 * kk_689[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, hk_690, hk_691, hk_692, hk_693, \
                         hk_694, kk_690, kk_691, kk_692, kk_693, \
                         kk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -hk_690[k]
                   + f_0 * kk_690[k];

        t_691[k] = -hk_691[k]
                   + f_0 * kk_691[k];

        t_692[k] = -hk_692[k]
                   + f_0 * kk_692[k];

        t_693[k] = -hk_693[k]
                   + f_0 * kk_693[k];

        t_694[k] = -hk_694[k]
                   + f_0 * kk_694[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, hk_695, hk_696, hk_697, hk_698, \
                         hk_699, kk_695, kk_696, kk_697, kk_698, \
                         kk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -hk_695[k]
                   + f_0 * kk_695[k];

        t_696[k] = -hk_696[k]
                   + f_0 * kk_696[k];

        t_697[k] = -hk_697[k]
                   + f_0 * kk_697[k];

        t_698[k] = -hk_698[k]
                   + f_0 * kk_698[k];

        t_699[k] = -hk_699[k]
                   + f_0 * kk_699[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, hk_700, hk_701, hk_702, hk_703, \
                         hk_704, kk_700, kk_701, kk_702, kk_703, \
                         kk_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -hk_700[k]
                   + f_0 * kk_700[k];

        t_701[k] = -hk_701[k]
                   + f_0 * kk_701[k];

        t_702[k] = -hk_702[k]
                   + f_0 * kk_702[k];

        t_703[k] = -hk_703[k]
                   + f_0 * kk_703[k];

        t_704[k] = -hk_704[k]
                   + f_0 * kk_704[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, hk_705, hk_706, hk_707, hk_708, \
                         hk_709, kk_705, kk_706, kk_707, kk_708, \
                         kk_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -hk_705[k]
                   + f_0 * kk_705[k];

        t_706[k] = -hk_706[k]
                   + f_0 * kk_706[k];

        t_707[k] = -hk_707[k]
                   + f_0 * kk_707[k];

        t_708[k] = -hk_708[k]
                   + f_0 * kk_708[k];

        t_709[k] = -hk_709[k]
                   + f_0 * kk_709[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, hk_710, hk_711, hk_712, hk_713, \
                         hk_714, kk_710, kk_711, kk_712, kk_713, \
                         kk_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -hk_710[k]
                   + f_0 * kk_710[k];

        t_711[k] = -hk_711[k]
                   + f_0 * kk_711[k];

        t_712[k] = -hk_712[k]
                   + f_0 * kk_712[k];

        t_713[k] = -hk_713[k]
                   + f_0 * kk_713[k];

        t_714[k] = -hk_714[k]
                   + f_0 * kk_714[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, hk_715, hk_716, hk_717, hk_718, \
                         hk_719, kk_715, kk_716, kk_717, kk_718, \
                         kk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -hk_715[k]
                   + f_0 * kk_715[k];

        t_716[k] = -hk_716[k]
                   + f_0 * kk_716[k];

        t_717[k] = -hk_717[k]
                   + f_0 * kk_717[k];

        t_718[k] = -hk_718[k]
                   + f_0 * kk_718[k];

        t_719[k] = -hk_719[k]
                   + f_0 * kk_719[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, hk_720, hk_721, hk_722, hk_723, \
                         hk_724, kk_720, kk_721, kk_722, kk_723, \
                         kk_724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -hk_720[k]
                   + f_0 * kk_720[k];

        t_721[k] = -hk_721[k]
                   + f_0 * kk_721[k];

        t_722[k] = -hk_722[k]
                   + f_0 * kk_722[k];

        t_723[k] = -hk_723[k]
                   + f_0 * kk_723[k];

        t_724[k] = -hk_724[k]
                   + f_0 * kk_724[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, hk_725, hk_726, hk_727, hk_728, \
                         hk_729, kk_725, kk_726, kk_727, kk_728, \
                         kk_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -hk_725[k]
                   + f_0 * kk_725[k];

        t_726[k] = -hk_726[k]
                   + f_0 * kk_726[k];

        t_727[k] = -hk_727[k]
                   + f_0 * kk_727[k];

        t_728[k] = -hk_728[k]
                   + f_0 * kk_728[k];

        t_729[k] = -hk_729[k]
                   + f_0 * kk_729[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, hk_730, hk_731, hk_732, hk_733, \
                         hk_734, kk_730, kk_731, kk_732, kk_733, \
                         kk_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -hk_730[k]
                   + f_0 * kk_730[k];

        t_731[k] = -hk_731[k]
                   + f_0 * kk_731[k];

        t_732[k] = -hk_732[k]
                   + f_0 * kk_732[k];

        t_733[k] = -hk_733[k]
                   + f_0 * kk_733[k];

        t_734[k] = -hk_734[k]
                   + f_0 * kk_734[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, hk_735, hk_736, hk_737, hk_738, \
                         hk_739, kk_735, kk_736, kk_737, kk_738, \
                         kk_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -hk_735[k]
                   + f_0 * kk_735[k];

        t_736[k] = -hk_736[k]
                   + f_0 * kk_736[k];

        t_737[k] = -hk_737[k]
                   + f_0 * kk_737[k];

        t_738[k] = -hk_738[k]
                   + f_0 * kk_738[k];

        t_739[k] = -hk_739[k]
                   + f_0 * kk_739[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, hk_740, hk_741, hk_742, hk_743, \
                         hk_744, kk_740, kk_741, kk_742, kk_743, \
                         kk_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -hk_740[k]
                   + f_0 * kk_740[k];

        t_741[k] = -hk_741[k]
                   + f_0 * kk_741[k];

        t_742[k] = -hk_742[k]
                   + f_0 * kk_742[k];

        t_743[k] = -hk_743[k]
                   + f_0 * kk_743[k];

        t_744[k] = -hk_744[k]
                   + f_0 * kk_744[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, hk_745, hk_746, hk_747, hk_748, \
                         hk_749, kk_745, kk_746, kk_747, kk_748, \
                         kk_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -hk_745[k]
                   + f_0 * kk_745[k];

        t_746[k] = -hk_746[k]
                   + f_0 * kk_746[k];

        t_747[k] = -hk_747[k]
                   + f_0 * kk_747[k];

        t_748[k] = -hk_748[k]
                   + f_0 * kk_748[k];

        t_749[k] = -hk_749[k]
                   + f_0 * kk_749[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_750 = buffer.data(hk + 750);
    const auto *hk_751 = buffer.data(hk + 751);
    const auto *hk_752 = buffer.data(hk + 752);
    const auto *hk_753 = buffer.data(hk + 753);
    const auto *hk_754 = buffer.data(hk + 754);
    const auto *hk_755 = buffer.data(hk + 755);

    const auto *kk_750 = buffer.data(kk + 750);
    const auto *kk_751 = buffer.data(kk + 751);
    const auto *kk_752 = buffer.data(kk + 752);
    const auto *kk_753 = buffer.data(kk + 753);
    const auto *kk_754 = buffer.data(kk + 754);
    const auto *kk_755 = buffer.data(kk + 755);
    const auto *kk_756 = buffer.data(kk + 756);
    const auto *kk_757 = buffer.data(kk + 757);
    const auto *kk_758 = buffer.data(kk + 758);
    const auto *kk_759 = buffer.data(kk + 759);
    const auto *kk_760 = buffer.data(kk + 760);
    const auto *kk_761 = buffer.data(kk + 761);
    const auto *kk_762 = buffer.data(kk + 762);
    const auto *kk_763 = buffer.data(kk + 763);
    const auto *kk_764 = buffer.data(kk + 764);
    const auto *kk_765 = buffer.data(kk + 765);
    const auto *kk_766 = buffer.data(kk + 766);
    const auto *kk_767 = buffer.data(kk + 767);
    const auto *kk_768 = buffer.data(kk + 768);
    const auto *kk_769 = buffer.data(kk + 769);
    const auto *kk_770 = buffer.data(kk + 770);
    const auto *kk_771 = buffer.data(kk + 771);
    const auto *kk_772 = buffer.data(kk + 772);
    const auto *kk_773 = buffer.data(kk + 773);
    const auto *kk_774 = buffer.data(kk + 774);
    const auto *kk_775 = buffer.data(kk + 775);
    const auto *kk_776 = buffer.data(kk + 776);
    const auto *kk_777 = buffer.data(kk + 777);
    const auto *kk_778 = buffer.data(kk + 778);
    const auto *kk_779 = buffer.data(kk + 779);
    const auto *kk_780 = buffer.data(kk + 780);
    const auto *kk_781 = buffer.data(kk + 781);
    const auto *kk_782 = buffer.data(kk + 782);
    const auto *kk_783 = buffer.data(kk + 783);
    const auto *kk_784 = buffer.data(kk + 784);
    const auto *kk_785 = buffer.data(kk + 785);
    const auto *kk_786 = buffer.data(kk + 786);
    const auto *kk_787 = buffer.data(kk + 787);
    const auto *kk_788 = buffer.data(kk + 788);
    const auto *kk_789 = buffer.data(kk + 789);
    const auto *kk_790 = buffer.data(kk + 790);
    const auto *kk_791 = buffer.data(kk + 791);
    const auto *kk_792 = buffer.data(kk + 792);
    const auto *kk_793 = buffer.data(kk + 793);
    const auto *kk_794 = buffer.data(kk + 794);
    const auto *kk_795 = buffer.data(kk + 795);
    const auto *kk_796 = buffer.data(kk + 796);
    const auto *kk_797 = buffer.data(kk + 797);
    const auto *kk_798 = buffer.data(kk + 798);
    const auto *kk_799 = buffer.data(kk + 799);
    const auto *kk_800 = buffer.data(kk + 800);
    const auto *kk_801 = buffer.data(kk + 801);
    const auto *kk_802 = buffer.data(kk + 802);
    const auto *kk_803 = buffer.data(kk + 803);
    const auto *kk_804 = buffer.data(kk + 804);
    const auto *kk_805 = buffer.data(kk + 805);
    const auto *kk_806 = buffer.data(kk + 806);
    const auto *kk_807 = buffer.data(kk + 807);
    const auto *kk_808 = buffer.data(kk + 808);
    const auto *kk_809 = buffer.data(kk + 809);
    const auto *kk_810 = buffer.data(kk + 810);
    const auto *kk_811 = buffer.data(kk + 811);
    const auto *kk_812 = buffer.data(kk + 812);
    const auto *kk_813 = buffer.data(kk + 813);
    const auto *kk_814 = buffer.data(kk + 814);
    const auto *kk_815 = buffer.data(kk + 815);
    const auto *kk_816 = buffer.data(kk + 816);
    const auto *kk_817 = buffer.data(kk + 817);
    const auto *kk_818 = buffer.data(kk + 818);
    const auto *kk_819 = buffer.data(kk + 819);
    const auto *kk_820 = buffer.data(kk + 820);
    const auto *kk_821 = buffer.data(kk + 821);
    const auto *kk_822 = buffer.data(kk + 822);
    const auto *kk_823 = buffer.data(kk + 823);
    const auto *kk_824 = buffer.data(kk + 824);
    const auto *kk_825 = buffer.data(kk + 825);
    const auto *kk_826 = buffer.data(kk + 826);
    const auto *kk_827 = buffer.data(kk + 827);
    const auto *kk_828 = buffer.data(kk + 828);
    const auto *kk_829 = buffer.data(kk + 829);
    const auto *kk_830 = buffer.data(kk + 830);
    const auto *kk_831 = buffer.data(kk + 831);
    const auto *kk_832 = buffer.data(kk + 832);
    const auto *kk_833 = buffer.data(kk + 833);
    const auto *kk_834 = buffer.data(kk + 834);
    const auto *kk_835 = buffer.data(kk + 835);
    const auto *kk_836 = buffer.data(kk + 836);
    const auto *kk_837 = buffer.data(kk + 837);
    const auto *kk_838 = buffer.data(kk + 838);
    const auto *kk_839 = buffer.data(kk + 839);
    const auto *kk_840 = buffer.data(kk + 840);
    const auto *kk_841 = buffer.data(kk + 841);
    const auto *kk_842 = buffer.data(kk + 842);
    const auto *kk_843 = buffer.data(kk + 843);
    const auto *kk_844 = buffer.data(kk + 844);
    const auto *kk_845 = buffer.data(kk + 845);
    const auto *kk_846 = buffer.data(kk + 846);
    const auto *kk_847 = buffer.data(kk + 847);
    const auto *kk_848 = buffer.data(kk + 848);
    const auto *kk_849 = buffer.data(kk + 849);
    const auto *kk_850 = buffer.data(kk + 850);
    const auto *kk_851 = buffer.data(kk + 851);
    const auto *kk_852 = buffer.data(kk + 852);
    const auto *kk_853 = buffer.data(kk + 853);
    const auto *kk_854 = buffer.data(kk + 854);
    const auto *kk_855 = buffer.data(kk + 855);
    const auto *kk_856 = buffer.data(kk + 856);
    const auto *kk_857 = buffer.data(kk + 857);
    const auto *kk_858 = buffer.data(kk + 858);
    const auto *kk_859 = buffer.data(kk + 859);
    const auto *kk_860 = buffer.data(kk + 860);
    const auto *kk_861 = buffer.data(kk + 861);
    const auto *kk_862 = buffer.data(kk + 862);
    const auto *kk_863 = buffer.data(kk + 863);
    const auto *kk_864 = buffer.data(kk + 864);
    const auto *kk_865 = buffer.data(kk + 865);
    const auto *kk_866 = buffer.data(kk + 866);
    const auto *kk_867 = buffer.data(kk + 867);
    const auto *kk_868 = buffer.data(kk + 868);
    const auto *kk_869 = buffer.data(kk + 869);
    const auto *kk_870 = buffer.data(kk + 870);
    const auto *kk_871 = buffer.data(kk + 871);
    const auto *kk_872 = buffer.data(kk + 872);
    const auto *kk_873 = buffer.data(kk + 873);
    const auto *kk_874 = buffer.data(kk + 874);
    const auto *kk_875 = buffer.data(kk + 875);
    const auto *kk_876 = buffer.data(kk + 876);
    const auto *kk_877 = buffer.data(kk + 877);
    const auto *kk_878 = buffer.data(kk + 878);
    const auto *kk_879 = buffer.data(kk + 879);
    const auto *kk_880 = buffer.data(kk + 880);
    const auto *kk_881 = buffer.data(kk + 881);
    const auto *kk_882 = buffer.data(kk + 882);
    const auto *kk_883 = buffer.data(kk + 883);
    const auto *kk_884 = buffer.data(kk + 884);
    const auto *kk_885 = buffer.data(kk + 885);
    const auto *kk_886 = buffer.data(kk + 886);
    const auto *kk_887 = buffer.data(kk + 887);
    const auto *kk_888 = buffer.data(kk + 888);
    const auto *kk_889 = buffer.data(kk + 889);
    const auto *kk_890 = buffer.data(kk + 890);
    const auto *kk_891 = buffer.data(kk + 891);
    const auto *kk_892 = buffer.data(kk + 892);
    const auto *kk_893 = buffer.data(kk + 893);
    const auto *kk_894 = buffer.data(kk + 894);
    const auto *kk_895 = buffer.data(kk + 895);
    const auto *kk_896 = buffer.data(kk + 896);
    const auto *kk_897 = buffer.data(kk + 897);
    const auto *kk_898 = buffer.data(kk + 898);
    const auto *kk_899 = buffer.data(kk + 899);
    const auto *kk_900 = buffer.data(kk + 900);
    const auto *kk_901 = buffer.data(kk + 901);
    const auto *kk_902 = buffer.data(kk + 902);
    const auto *kk_903 = buffer.data(kk + 903);
    const auto *kk_904 = buffer.data(kk + 904);
    const auto *kk_905 = buffer.data(kk + 905);
    const auto *kk_906 = buffer.data(kk + 906);
    const auto *kk_907 = buffer.data(kk + 907);
    const auto *kk_908 = buffer.data(kk + 908);
    const auto *kk_909 = buffer.data(kk + 909);
    const auto *kk_910 = buffer.data(kk + 910);
    const auto *kk_911 = buffer.data(kk + 911);
    const auto *kk_912 = buffer.data(kk + 912);
    const auto *kk_913 = buffer.data(kk + 913);
    const auto *kk_914 = buffer.data(kk + 914);
    const auto *kk_915 = buffer.data(kk + 915);
    const auto *kk_916 = buffer.data(kk + 916);
    const auto *kk_917 = buffer.data(kk + 917);
    const auto *kk_918 = buffer.data(kk + 918);
    const auto *kk_919 = buffer.data(kk + 919);
    const auto *kk_920 = buffer.data(kk + 920);
    const auto *kk_921 = buffer.data(kk + 921);
    const auto *kk_922 = buffer.data(kk + 922);
    const auto *kk_923 = buffer.data(kk + 923);
    const auto *kk_924 = buffer.data(kk + 924);
    const auto *kk_925 = buffer.data(kk + 925);
    const auto *kk_926 = buffer.data(kk + 926);
    const auto *kk_927 = buffer.data(kk + 927);
    const auto *kk_928 = buffer.data(kk + 928);
    const auto *kk_929 = buffer.data(kk + 929);
    const auto *kk_930 = buffer.data(kk + 930);
    const auto *kk_931 = buffer.data(kk + 931);
    const auto *kk_932 = buffer.data(kk + 932);
    const auto *kk_933 = buffer.data(kk + 933);
    const auto *kk_934 = buffer.data(kk + 934);
    const auto *kk_935 = buffer.data(kk + 935);
    const auto *kk_936 = buffer.data(kk + 936);
    const auto *kk_937 = buffer.data(kk + 937);
    const auto *kk_938 = buffer.data(kk + 938);
    const auto *kk_939 = buffer.data(kk + 939);
    const auto *kk_940 = buffer.data(kk + 940);
    const auto *kk_941 = buffer.data(kk + 941);
    const auto *kk_942 = buffer.data(kk + 942);
    const auto *kk_943 = buffer.data(kk + 943);
    const auto *kk_944 = buffer.data(kk + 944);
    const auto *kk_945 = buffer.data(kk + 945);
    const auto *kk_946 = buffer.data(kk + 946);
    const auto *kk_947 = buffer.data(kk + 947);
    const auto *kk_948 = buffer.data(kk + 948);
    const auto *kk_949 = buffer.data(kk + 949);
    const auto *kk_950 = buffer.data(kk + 950);
    const auto *kk_951 = buffer.data(kk + 951);
    const auto *kk_952 = buffer.data(kk + 952);
    const auto *kk_953 = buffer.data(kk + 953);
    const auto *kk_954 = buffer.data(kk + 954);
    const auto *kk_955 = buffer.data(kk + 955);
    const auto *kk_956 = buffer.data(kk + 956);
    const auto *kk_957 = buffer.data(kk + 957);
    const auto *kk_958 = buffer.data(kk + 958);
    const auto *kk_959 = buffer.data(kk + 959);
    const auto *kk_960 = buffer.data(kk + 960);
    const auto *kk_961 = buffer.data(kk + 961);
    const auto *kk_962 = buffer.data(kk + 962);
    const auto *kk_963 = buffer.data(kk + 963);
    const auto *kk_964 = buffer.data(kk + 964);
    const auto *kk_965 = buffer.data(kk + 965);
    const auto *kk_966 = buffer.data(kk + 966);
    const auto *kk_967 = buffer.data(kk + 967);
    const auto *kk_968 = buffer.data(kk + 968);
    const auto *kk_969 = buffer.data(kk + 969);

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, hk_750, hk_751, hk_752, hk_753, \
                         hk_754, kk_750, kk_751, kk_752, kk_753, \
                         kk_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -hk_750[k]
                   + f_0 * kk_750[k];

        t_751[k] = -hk_751[k]
                   + f_0 * kk_751[k];

        t_752[k] = -hk_752[k]
                   + f_0 * kk_752[k];

        t_753[k] = -hk_753[k]
                   + f_0 * kk_753[k];

        t_754[k] = -hk_754[k]
                   + f_0 * kk_754[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, t_760, t_761, hk_755, kk_755, \
                         kk_756, kk_757, kk_758, kk_759, kk_760, \
                         kk_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -hk_755[k]
                   + f_0 * kk_755[k];

        t_756[k] = f_0 * kk_756[k];

        t_757[k] = f_0 * kk_757[k];

        t_758[k] = f_0 * kk_758[k];

        t_759[k] = f_0 * kk_759[k];

        t_760[k] = f_0 * kk_760[k];

        t_761[k] = f_0 * kk_761[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, t_766, t_767, t_768, t_769, kk_762, \
                         kk_763, kk_764, kk_765, kk_766, kk_767, kk_768, \
                         kk_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_0 * kk_762[k];

        t_763[k] = f_0 * kk_763[k];

        t_764[k] = f_0 * kk_764[k];

        t_765[k] = f_0 * kk_765[k];

        t_766[k] = f_0 * kk_766[k];

        t_767[k] = f_0 * kk_767[k];

        t_768[k] = f_0 * kk_768[k];

        t_769[k] = f_0 * kk_769[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, t_775, t_776, t_777, kk_770, \
                         kk_771, kk_772, kk_773, kk_774, kk_775, kk_776, \
                         kk_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_0 * kk_770[k];

        t_771[k] = f_0 * kk_771[k];

        t_772[k] = f_0 * kk_772[k];

        t_773[k] = f_0 * kk_773[k];

        t_774[k] = f_0 * kk_774[k];

        t_775[k] = f_0 * kk_775[k];

        t_776[k] = f_0 * kk_776[k];

        t_777[k] = f_0 * kk_777[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, t_783, t_784, t_785, kk_778, \
                         kk_779, kk_780, kk_781, kk_782, kk_783, kk_784, \
                         kk_785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_0 * kk_778[k];

        t_779[k] = f_0 * kk_779[k];

        t_780[k] = f_0 * kk_780[k];

        t_781[k] = f_0 * kk_781[k];

        t_782[k] = f_0 * kk_782[k];

        t_783[k] = f_0 * kk_783[k];

        t_784[k] = f_0 * kk_784[k];

        t_785[k] = f_0 * kk_785[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, t_790, t_791, t_792, t_793, kk_786, \
                         kk_787, kk_788, kk_789, kk_790, kk_791, kk_792, \
                         kk_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_0 * kk_786[k];

        t_787[k] = f_0 * kk_787[k];

        t_788[k] = f_0 * kk_788[k];

        t_789[k] = f_0 * kk_789[k];

        t_790[k] = f_0 * kk_790[k];

        t_791[k] = f_0 * kk_791[k];

        t_792[k] = f_0 * kk_792[k];

        t_793[k] = f_0 * kk_793[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, t_798, t_799, t_800, t_801, kk_794, \
                         kk_795, kk_796, kk_797, kk_798, kk_799, kk_800, \
                         kk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_0 * kk_794[k];

        t_795[k] = f_0 * kk_795[k];

        t_796[k] = f_0 * kk_796[k];

        t_797[k] = f_0 * kk_797[k];

        t_798[k] = f_0 * kk_798[k];

        t_799[k] = f_0 * kk_799[k];

        t_800[k] = f_0 * kk_800[k];

        t_801[k] = f_0 * kk_801[k];
    }

#pragma omp simd aligned(t_802, t_803, t_804, t_805, t_806, t_807, t_808, t_809, kk_802, \
                         kk_803, kk_804, kk_805, kk_806, kk_807, kk_808, \
                         kk_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_802[k] = f_0 * kk_802[k];

        t_803[k] = f_0 * kk_803[k];

        t_804[k] = f_0 * kk_804[k];

        t_805[k] = f_0 * kk_805[k];

        t_806[k] = f_0 * kk_806[k];

        t_807[k] = f_0 * kk_807[k];

        t_808[k] = f_0 * kk_808[k];

        t_809[k] = f_0 * kk_809[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, t_815, t_816, t_817, kk_810, \
                         kk_811, kk_812, kk_813, kk_814, kk_815, kk_816, \
                         kk_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_0 * kk_810[k];

        t_811[k] = f_0 * kk_811[k];

        t_812[k] = f_0 * kk_812[k];

        t_813[k] = f_0 * kk_813[k];

        t_814[k] = f_0 * kk_814[k];

        t_815[k] = f_0 * kk_815[k];

        t_816[k] = f_0 * kk_816[k];

        t_817[k] = f_0 * kk_817[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, t_821, t_822, t_823, t_824, t_825, kk_818, \
                         kk_819, kk_820, kk_821, kk_822, kk_823, kk_824, \
                         kk_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_0 * kk_818[k];

        t_819[k] = f_0 * kk_819[k];

        t_820[k] = f_0 * kk_820[k];

        t_821[k] = f_0 * kk_821[k];

        t_822[k] = f_0 * kk_822[k];

        t_823[k] = f_0 * kk_823[k];

        t_824[k] = f_0 * kk_824[k];

        t_825[k] = f_0 * kk_825[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, t_830, t_831, t_832, t_833, kk_826, \
                         kk_827, kk_828, kk_829, kk_830, kk_831, kk_832, \
                         kk_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_0 * kk_826[k];

        t_827[k] = f_0 * kk_827[k];

        t_828[k] = f_0 * kk_828[k];

        t_829[k] = f_0 * kk_829[k];

        t_830[k] = f_0 * kk_830[k];

        t_831[k] = f_0 * kk_831[k];

        t_832[k] = f_0 * kk_832[k];

        t_833[k] = f_0 * kk_833[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, t_837, t_838, t_839, t_840, t_841, kk_834, \
                         kk_835, kk_836, kk_837, kk_838, kk_839, kk_840, \
                         kk_841 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_0 * kk_834[k];

        t_835[k] = f_0 * kk_835[k];

        t_836[k] = f_0 * kk_836[k];

        t_837[k] = f_0 * kk_837[k];

        t_838[k] = f_0 * kk_838[k];

        t_839[k] = f_0 * kk_839[k];

        t_840[k] = f_0 * kk_840[k];

        t_841[k] = f_0 * kk_841[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, t_845, t_846, t_847, t_848, t_849, kk_842, \
                         kk_843, kk_844, kk_845, kk_846, kk_847, kk_848, \
                         kk_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_0 * kk_842[k];

        t_843[k] = f_0 * kk_843[k];

        t_844[k] = f_0 * kk_844[k];

        t_845[k] = f_0 * kk_845[k];

        t_846[k] = f_0 * kk_846[k];

        t_847[k] = f_0 * kk_847[k];

        t_848[k] = f_0 * kk_848[k];

        t_849[k] = f_0 * kk_849[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, t_855, t_856, t_857, kk_850, \
                         kk_851, kk_852, kk_853, kk_854, kk_855, kk_856, \
                         kk_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = f_0 * kk_850[k];

        t_851[k] = f_0 * kk_851[k];

        t_852[k] = f_0 * kk_852[k];

        t_853[k] = f_0 * kk_853[k];

        t_854[k] = f_0 * kk_854[k];

        t_855[k] = f_0 * kk_855[k];

        t_856[k] = f_0 * kk_856[k];

        t_857[k] = f_0 * kk_857[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, t_862, t_863, t_864, t_865, kk_858, \
                         kk_859, kk_860, kk_861, kk_862, kk_863, kk_864, \
                         kk_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_0 * kk_858[k];

        t_859[k] = f_0 * kk_859[k];

        t_860[k] = f_0 * kk_860[k];

        t_861[k] = f_0 * kk_861[k];

        t_862[k] = f_0 * kk_862[k];

        t_863[k] = f_0 * kk_863[k];

        t_864[k] = f_0 * kk_864[k];

        t_865[k] = f_0 * kk_865[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, t_870, t_871, t_872, t_873, kk_866, \
                         kk_867, kk_868, kk_869, kk_870, kk_871, kk_872, \
                         kk_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_0 * kk_866[k];

        t_867[k] = f_0 * kk_867[k];

        t_868[k] = f_0 * kk_868[k];

        t_869[k] = f_0 * kk_869[k];

        t_870[k] = f_0 * kk_870[k];

        t_871[k] = f_0 * kk_871[k];

        t_872[k] = f_0 * kk_872[k];

        t_873[k] = f_0 * kk_873[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, t_878, t_879, t_880, t_881, kk_874, \
                         kk_875, kk_876, kk_877, kk_878, kk_879, kk_880, \
                         kk_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_0 * kk_874[k];

        t_875[k] = f_0 * kk_875[k];

        t_876[k] = f_0 * kk_876[k];

        t_877[k] = f_0 * kk_877[k];

        t_878[k] = f_0 * kk_878[k];

        t_879[k] = f_0 * kk_879[k];

        t_880[k] = f_0 * kk_880[k];

        t_881[k] = f_0 * kk_881[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, t_886, t_887, t_888, t_889, kk_882, \
                         kk_883, kk_884, kk_885, kk_886, kk_887, kk_888, \
                         kk_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_0 * kk_882[k];

        t_883[k] = f_0 * kk_883[k];

        t_884[k] = f_0 * kk_884[k];

        t_885[k] = f_0 * kk_885[k];

        t_886[k] = f_0 * kk_886[k];

        t_887[k] = f_0 * kk_887[k];

        t_888[k] = f_0 * kk_888[k];

        t_889[k] = f_0 * kk_889[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, t_895, t_896, t_897, kk_890, \
                         kk_891, kk_892, kk_893, kk_894, kk_895, kk_896, \
                         kk_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_0 * kk_890[k];

        t_891[k] = f_0 * kk_891[k];

        t_892[k] = f_0 * kk_892[k];

        t_893[k] = f_0 * kk_893[k];

        t_894[k] = f_0 * kk_894[k];

        t_895[k] = f_0 * kk_895[k];

        t_896[k] = f_0 * kk_896[k];

        t_897[k] = f_0 * kk_897[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, t_902, t_903, t_904, t_905, kk_898, \
                         kk_899, kk_900, kk_901, kk_902, kk_903, kk_904, \
                         kk_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_0 * kk_898[k];

        t_899[k] = f_0 * kk_899[k];

        t_900[k] = f_0 * kk_900[k];

        t_901[k] = f_0 * kk_901[k];

        t_902[k] = f_0 * kk_902[k];

        t_903[k] = f_0 * kk_903[k];

        t_904[k] = f_0 * kk_904[k];

        t_905[k] = f_0 * kk_905[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, t_909, t_910, t_911, t_912, t_913, kk_906, \
                         kk_907, kk_908, kk_909, kk_910, kk_911, kk_912, \
                         kk_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_0 * kk_906[k];

        t_907[k] = f_0 * kk_907[k];

        t_908[k] = f_0 * kk_908[k];

        t_909[k] = f_0 * kk_909[k];

        t_910[k] = f_0 * kk_910[k];

        t_911[k] = f_0 * kk_911[k];

        t_912[k] = f_0 * kk_912[k];

        t_913[k] = f_0 * kk_913[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, t_918, t_919, t_920, t_921, kk_914, \
                         kk_915, kk_916, kk_917, kk_918, kk_919, kk_920, \
                         kk_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_0 * kk_914[k];

        t_915[k] = f_0 * kk_915[k];

        t_916[k] = f_0 * kk_916[k];

        t_917[k] = f_0 * kk_917[k];

        t_918[k] = f_0 * kk_918[k];

        t_919[k] = f_0 * kk_919[k];

        t_920[k] = f_0 * kk_920[k];

        t_921[k] = f_0 * kk_921[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, t_926, t_927, t_928, t_929, kk_922, \
                         kk_923, kk_924, kk_925, kk_926, kk_927, kk_928, \
                         kk_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_0 * kk_922[k];

        t_923[k] = f_0 * kk_923[k];

        t_924[k] = f_0 * kk_924[k];

        t_925[k] = f_0 * kk_925[k];

        t_926[k] = f_0 * kk_926[k];

        t_927[k] = f_0 * kk_927[k];

        t_928[k] = f_0 * kk_928[k];

        t_929[k] = f_0 * kk_929[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, t_935, t_936, t_937, kk_930, \
                         kk_931, kk_932, kk_933, kk_934, kk_935, kk_936, \
                         kk_937 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_0 * kk_930[k];

        t_931[k] = f_0 * kk_931[k];

        t_932[k] = f_0 * kk_932[k];

        t_933[k] = f_0 * kk_933[k];

        t_934[k] = f_0 * kk_934[k];

        t_935[k] = f_0 * kk_935[k];

        t_936[k] = f_0 * kk_936[k];

        t_937[k] = f_0 * kk_937[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, t_942, t_943, t_944, t_945, kk_938, \
                         kk_939, kk_940, kk_941, kk_942, kk_943, kk_944, \
                         kk_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_0 * kk_938[k];

        t_939[k] = f_0 * kk_939[k];

        t_940[k] = f_0 * kk_940[k];

        t_941[k] = f_0 * kk_941[k];

        t_942[k] = f_0 * kk_942[k];

        t_943[k] = f_0 * kk_943[k];

        t_944[k] = f_0 * kk_944[k];

        t_945[k] = f_0 * kk_945[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, t_949, t_950, t_951, t_952, t_953, kk_946, \
                         kk_947, kk_948, kk_949, kk_950, kk_951, kk_952, \
                         kk_953 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_0 * kk_946[k];

        t_947[k] = f_0 * kk_947[k];

        t_948[k] = f_0 * kk_948[k];

        t_949[k] = f_0 * kk_949[k];

        t_950[k] = f_0 * kk_950[k];

        t_951[k] = f_0 * kk_951[k];

        t_952[k] = f_0 * kk_952[k];

        t_953[k] = f_0 * kk_953[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, t_957, t_958, t_959, t_960, t_961, kk_954, \
                         kk_955, kk_956, kk_957, kk_958, kk_959, kk_960, \
                         kk_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_0 * kk_954[k];

        t_955[k] = f_0 * kk_955[k];

        t_956[k] = f_0 * kk_956[k];

        t_957[k] = f_0 * kk_957[k];

        t_958[k] = f_0 * kk_958[k];

        t_959[k] = f_0 * kk_959[k];

        t_960[k] = f_0 * kk_960[k];

        t_961[k] = f_0 * kk_961[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, t_966, t_967, t_968, t_969, kk_962, \
                         kk_963, kk_964, kk_965, kk_966, kk_967, kk_968, \
                         kk_969 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_0 * kk_962[k];

        t_963[k] = f_0 * kk_963[k];

        t_964[k] = f_0 * kk_964[k];

        t_965[k] = f_0 * kk_965[k];

        t_966[k] = f_0 * kk_966[k];

        t_967[k] = f_0 * kk_967[k];

        t_968[k] = f_0 * kk_968[k];

        t_969[k] = f_0 * kk_969[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_0_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kk, const size_t ncols,
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

    const auto *kk_970 = buffer.data(kk + 970);
    const auto *kk_971 = buffer.data(kk + 971);
    const auto *kk_972 = buffer.data(kk + 972);
    const auto *kk_973 = buffer.data(kk + 973);
    const auto *kk_974 = buffer.data(kk + 974);
    const auto *kk_975 = buffer.data(kk + 975);
    const auto *kk_976 = buffer.data(kk + 976);
    const auto *kk_977 = buffer.data(kk + 977);
    const auto *kk_978 = buffer.data(kk + 978);
    const auto *kk_979 = buffer.data(kk + 979);
    const auto *kk_980 = buffer.data(kk + 980);
    const auto *kk_981 = buffer.data(kk + 981);
    const auto *kk_982 = buffer.data(kk + 982);
    const auto *kk_983 = buffer.data(kk + 983);
    const auto *kk_984 = buffer.data(kk + 984);
    const auto *kk_985 = buffer.data(kk + 985);
    const auto *kk_986 = buffer.data(kk + 986);
    const auto *kk_987 = buffer.data(kk + 987);
    const auto *kk_988 = buffer.data(kk + 988);
    const auto *kk_989 = buffer.data(kk + 989);
    const auto *kk_990 = buffer.data(kk + 990);
    const auto *kk_991 = buffer.data(kk + 991);
    const auto *kk_992 = buffer.data(kk + 992);
    const auto *kk_993 = buffer.data(kk + 993);
    const auto *kk_994 = buffer.data(kk + 994);
    const auto *kk_995 = buffer.data(kk + 995);
    const auto *kk_996 = buffer.data(kk + 996);
    const auto *kk_997 = buffer.data(kk + 997);
    const auto *kk_998 = buffer.data(kk + 998);
    const auto *kk_999 = buffer.data(kk + 999);
    const auto *kk_1000 = buffer.data(kk + 1000);
    const auto *kk_1001 = buffer.data(kk + 1001);
    const auto *kk_1002 = buffer.data(kk + 1002);
    const auto *kk_1003 = buffer.data(kk + 1003);
    const auto *kk_1004 = buffer.data(kk + 1004);
    const auto *kk_1005 = buffer.data(kk + 1005);
    const auto *kk_1006 = buffer.data(kk + 1006);
    const auto *kk_1007 = buffer.data(kk + 1007);

#pragma omp simd aligned(t_970, t_971, t_972, t_973, t_974, t_975, t_976, t_977, kk_970, \
                         kk_971, kk_972, kk_973, kk_974, kk_975, kk_976, \
                         kk_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = f_0 * kk_970[k];

        t_971[k] = f_0 * kk_971[k];

        t_972[k] = f_0 * kk_972[k];

        t_973[k] = f_0 * kk_973[k];

        t_974[k] = f_0 * kk_974[k];

        t_975[k] = f_0 * kk_975[k];

        t_976[k] = f_0 * kk_976[k];

        t_977[k] = f_0 * kk_977[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, t_982, t_983, t_984, t_985, kk_978, \
                         kk_979, kk_980, kk_981, kk_982, kk_983, kk_984, \
                         kk_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = f_0 * kk_978[k];

        t_979[k] = f_0 * kk_979[k];

        t_980[k] = f_0 * kk_980[k];

        t_981[k] = f_0 * kk_981[k];

        t_982[k] = f_0 * kk_982[k];

        t_983[k] = f_0 * kk_983[k];

        t_984[k] = f_0 * kk_984[k];

        t_985[k] = f_0 * kk_985[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, t_989, t_990, t_991, t_992, t_993, kk_986, \
                         kk_987, kk_988, kk_989, kk_990, kk_991, kk_992, \
                         kk_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_0 * kk_986[k];

        t_987[k] = f_0 * kk_987[k];

        t_988[k] = f_0 * kk_988[k];

        t_989[k] = f_0 * kk_989[k];

        t_990[k] = f_0 * kk_990[k];

        t_991[k] = f_0 * kk_991[k];

        t_992[k] = f_0 * kk_992[k];

        t_993[k] = f_0 * kk_993[k];
    }

#pragma omp simd aligned(t_994, t_995, t_996, t_997, t_998, t_999, t_1000, t_1001, kk_994, \
                         kk_995, kk_996, kk_997, kk_998, kk_999, kk_1000, \
                         kk_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_994[k] = f_0 * kk_994[k];

        t_995[k] = f_0 * kk_995[k];

        t_996[k] = f_0 * kk_996[k];

        t_997[k] = f_0 * kk_997[k];

        t_998[k] = f_0 * kk_998[k];

        t_999[k] = f_0 * kk_999[k];

        t_1000[k] = f_0 * kk_1000[k];

        t_1001[k] = f_0 * kk_1001[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, t_1005, t_1006, t_1007, kk_1002, kk_1003, \
                         kk_1004, kk_1005, kk_1006, kk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_0 * kk_1002[k];

        t_1003[k] = f_0 * kk_1003[k];

        t_1004[k] = f_0 * kk_1004[k];

        t_1005[k] = f_0 * kk_1005[k];

        t_1006[k] = f_0 * kk_1006[k];

        t_1007[k] = f_0 * kk_1007[k];
    }
}

auto
compute_prim_geom_10_ik_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t hk, const size_t kk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ik_electron_repulsion_0_piece0(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_0_piece1(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_0_piece2(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_0_piece3(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_0_piece4(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_0_piece5(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_0_piece6(buffer, target, kk, ncols, alpha);
}

static auto
compute_prim_geom_10_ik_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_1 = buffer.data(hk + 1);
    const auto *hk_2 = buffer.data(hk + 2);
    const auto *hk_3 = buffer.data(hk + 3);
    const auto *hk_4 = buffer.data(hk + 4);
    const auto *hk_5 = buffer.data(hk + 5);
    const auto *hk_6 = buffer.data(hk + 6);
    const auto *hk_7 = buffer.data(hk + 7);
    const auto *hk_8 = buffer.data(hk + 8);
    const auto *hk_9 = buffer.data(hk + 9);
    const auto *hk_10 = buffer.data(hk + 10);
    const auto *hk_11 = buffer.data(hk + 11);
    const auto *hk_12 = buffer.data(hk + 12);
    const auto *hk_13 = buffer.data(hk + 13);
    const auto *hk_14 = buffer.data(hk + 14);
    const auto *hk_15 = buffer.data(hk + 15);
    const auto *hk_16 = buffer.data(hk + 16);
    const auto *hk_17 = buffer.data(hk + 17);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_21 = buffer.data(hk + 21);
    const auto *hk_22 = buffer.data(hk + 22);
    const auto *hk_23 = buffer.data(hk + 23);
    const auto *hk_24 = buffer.data(hk + 24);
    const auto *hk_25 = buffer.data(hk + 25);
    const auto *hk_26 = buffer.data(hk + 26);
    const auto *hk_27 = buffer.data(hk + 27);
    const auto *hk_28 = buffer.data(hk + 28);
    const auto *hk_29 = buffer.data(hk + 29);
    const auto *hk_30 = buffer.data(hk + 30);
    const auto *hk_31 = buffer.data(hk + 31);
    const auto *hk_32 = buffer.data(hk + 32);
    const auto *hk_33 = buffer.data(hk + 33);
    const auto *hk_34 = buffer.data(hk + 34);
    const auto *hk_35 = buffer.data(hk + 35);
    const auto *hk_36 = buffer.data(hk + 36);
    const auto *hk_37 = buffer.data(hk + 37);
    const auto *hk_38 = buffer.data(hk + 38);
    const auto *hk_39 = buffer.data(hk + 39);
    const auto *hk_40 = buffer.data(hk + 40);
    const auto *hk_41 = buffer.data(hk + 41);
    const auto *hk_42 = buffer.data(hk + 42);
    const auto *hk_43 = buffer.data(hk + 43);
    const auto *hk_44 = buffer.data(hk + 44);
    const auto *hk_45 = buffer.data(hk + 45);
    const auto *hk_46 = buffer.data(hk + 46);
    const auto *hk_47 = buffer.data(hk + 47);
    const auto *hk_48 = buffer.data(hk + 48);
    const auto *hk_49 = buffer.data(hk + 49);
    const auto *hk_50 = buffer.data(hk + 50);
    const auto *hk_51 = buffer.data(hk + 51);
    const auto *hk_52 = buffer.data(hk + 52);
    const auto *hk_53 = buffer.data(hk + 53);
    const auto *hk_54 = buffer.data(hk + 54);
    const auto *hk_55 = buffer.data(hk + 55);
    const auto *hk_56 = buffer.data(hk + 56);
    const auto *hk_57 = buffer.data(hk + 57);
    const auto *hk_58 = buffer.data(hk + 58);
    const auto *hk_59 = buffer.data(hk + 59);
    const auto *hk_60 = buffer.data(hk + 60);
    const auto *hk_61 = buffer.data(hk + 61);
    const auto *hk_62 = buffer.data(hk + 62);
    const auto *hk_63 = buffer.data(hk + 63);
    const auto *hk_64 = buffer.data(hk + 64);
    const auto *hk_65 = buffer.data(hk + 65);
    const auto *hk_66 = buffer.data(hk + 66);
    const auto *hk_67 = buffer.data(hk + 67);
    const auto *hk_68 = buffer.data(hk + 68);
    const auto *hk_69 = buffer.data(hk + 69);
    const auto *hk_70 = buffer.data(hk + 70);
    const auto *hk_71 = buffer.data(hk + 71);
    const auto *hk_72 = buffer.data(hk + 72);
    const auto *hk_73 = buffer.data(hk + 73);
    const auto *hk_74 = buffer.data(hk + 74);
    const auto *hk_75 = buffer.data(hk + 75);
    const auto *hk_76 = buffer.data(hk + 76);
    const auto *hk_77 = buffer.data(hk + 77);
    const auto *hk_78 = buffer.data(hk + 78);
    const auto *hk_79 = buffer.data(hk + 79);
    const auto *hk_80 = buffer.data(hk + 80);
    const auto *hk_81 = buffer.data(hk + 81);
    const auto *hk_82 = buffer.data(hk + 82);
    const auto *hk_83 = buffer.data(hk + 83);
    const auto *hk_84 = buffer.data(hk + 84);
    const auto *hk_85 = buffer.data(hk + 85);
    const auto *hk_86 = buffer.data(hk + 86);
    const auto *hk_87 = buffer.data(hk + 87);
    const auto *hk_88 = buffer.data(hk + 88);
    const auto *hk_89 = buffer.data(hk + 89);
    const auto *hk_90 = buffer.data(hk + 90);
    const auto *hk_91 = buffer.data(hk + 91);
    const auto *hk_92 = buffer.data(hk + 92);
    const auto *hk_93 = buffer.data(hk + 93);
    const auto *hk_94 = buffer.data(hk + 94);
    const auto *hk_95 = buffer.data(hk + 95);
    const auto *hk_96 = buffer.data(hk + 96);
    const auto *hk_97 = buffer.data(hk + 97);
    const auto *hk_98 = buffer.data(hk + 98);
    const auto *hk_99 = buffer.data(hk + 99);

    const auto *kk_36 = buffer.data(kk + 36);
    const auto *kk_37 = buffer.data(kk + 37);
    const auto *kk_38 = buffer.data(kk + 38);
    const auto *kk_39 = buffer.data(kk + 39);
    const auto *kk_40 = buffer.data(kk + 40);
    const auto *kk_41 = buffer.data(kk + 41);
    const auto *kk_42 = buffer.data(kk + 42);
    const auto *kk_43 = buffer.data(kk + 43);
    const auto *kk_44 = buffer.data(kk + 44);
    const auto *kk_45 = buffer.data(kk + 45);
    const auto *kk_46 = buffer.data(kk + 46);
    const auto *kk_47 = buffer.data(kk + 47);
    const auto *kk_48 = buffer.data(kk + 48);
    const auto *kk_49 = buffer.data(kk + 49);
    const auto *kk_50 = buffer.data(kk + 50);
    const auto *kk_51 = buffer.data(kk + 51);
    const auto *kk_52 = buffer.data(kk + 52);
    const auto *kk_53 = buffer.data(kk + 53);
    const auto *kk_54 = buffer.data(kk + 54);
    const auto *kk_55 = buffer.data(kk + 55);
    const auto *kk_56 = buffer.data(kk + 56);
    const auto *kk_57 = buffer.data(kk + 57);
    const auto *kk_58 = buffer.data(kk + 58);
    const auto *kk_59 = buffer.data(kk + 59);
    const auto *kk_60 = buffer.data(kk + 60);
    const auto *kk_61 = buffer.data(kk + 61);
    const auto *kk_62 = buffer.data(kk + 62);
    const auto *kk_63 = buffer.data(kk + 63);
    const auto *kk_64 = buffer.data(kk + 64);
    const auto *kk_65 = buffer.data(kk + 65);
    const auto *kk_66 = buffer.data(kk + 66);
    const auto *kk_67 = buffer.data(kk + 67);
    const auto *kk_68 = buffer.data(kk + 68);
    const auto *kk_69 = buffer.data(kk + 69);
    const auto *kk_70 = buffer.data(kk + 70);
    const auto *kk_71 = buffer.data(kk + 71);
    const auto *kk_108 = buffer.data(kk + 108);
    const auto *kk_109 = buffer.data(kk + 109);
    const auto *kk_110 = buffer.data(kk + 110);
    const auto *kk_111 = buffer.data(kk + 111);
    const auto *kk_112 = buffer.data(kk + 112);
    const auto *kk_113 = buffer.data(kk + 113);
    const auto *kk_114 = buffer.data(kk + 114);
    const auto *kk_115 = buffer.data(kk + 115);
    const auto *kk_116 = buffer.data(kk + 116);
    const auto *kk_117 = buffer.data(kk + 117);
    const auto *kk_118 = buffer.data(kk + 118);
    const auto *kk_119 = buffer.data(kk + 119);
    const auto *kk_120 = buffer.data(kk + 120);
    const auto *kk_121 = buffer.data(kk + 121);
    const auto *kk_122 = buffer.data(kk + 122);
    const auto *kk_123 = buffer.data(kk + 123);
    const auto *kk_124 = buffer.data(kk + 124);
    const auto *kk_125 = buffer.data(kk + 125);
    const auto *kk_126 = buffer.data(kk + 126);
    const auto *kk_127 = buffer.data(kk + 127);
    const auto *kk_128 = buffer.data(kk + 128);
    const auto *kk_129 = buffer.data(kk + 129);
    const auto *kk_130 = buffer.data(kk + 130);
    const auto *kk_131 = buffer.data(kk + 131);
    const auto *kk_132 = buffer.data(kk + 132);
    const auto *kk_133 = buffer.data(kk + 133);
    const auto *kk_134 = buffer.data(kk + 134);
    const auto *kk_135 = buffer.data(kk + 135);
    const auto *kk_136 = buffer.data(kk + 136);
    const auto *kk_137 = buffer.data(kk + 137);
    const auto *kk_138 = buffer.data(kk + 138);
    const auto *kk_139 = buffer.data(kk + 139);
    const auto *kk_140 = buffer.data(kk + 140);
    const auto *kk_141 = buffer.data(kk + 141);
    const auto *kk_142 = buffer.data(kk + 142);
    const auto *kk_143 = buffer.data(kk + 143);
    const auto *kk_144 = buffer.data(kk + 144);
    const auto *kk_145 = buffer.data(kk + 145);
    const auto *kk_146 = buffer.data(kk + 146);
    const auto *kk_147 = buffer.data(kk + 147);
    const auto *kk_148 = buffer.data(kk + 148);
    const auto *kk_149 = buffer.data(kk + 149);
    const auto *kk_150 = buffer.data(kk + 150);
    const auto *kk_151 = buffer.data(kk + 151);
    const auto *kk_152 = buffer.data(kk + 152);
    const auto *kk_153 = buffer.data(kk + 153);
    const auto *kk_154 = buffer.data(kk + 154);
    const auto *kk_155 = buffer.data(kk + 155);
    const auto *kk_156 = buffer.data(kk + 156);
    const auto *kk_157 = buffer.data(kk + 157);
    const auto *kk_158 = buffer.data(kk + 158);
    const auto *kk_159 = buffer.data(kk + 159);
    const auto *kk_160 = buffer.data(kk + 160);
    const auto *kk_161 = buffer.data(kk + 161);
    const auto *kk_162 = buffer.data(kk + 162);
    const auto *kk_163 = buffer.data(kk + 163);
    const auto *kk_164 = buffer.data(kk + 164);
    const auto *kk_165 = buffer.data(kk + 165);
    const auto *kk_166 = buffer.data(kk + 166);
    const auto *kk_167 = buffer.data(kk + 167);
    const auto *kk_168 = buffer.data(kk + 168);
    const auto *kk_169 = buffer.data(kk + 169);
    const auto *kk_170 = buffer.data(kk + 170);
    const auto *kk_171 = buffer.data(kk + 171);
    const auto *kk_172 = buffer.data(kk + 172);
    const auto *kk_173 = buffer.data(kk + 173);
    const auto *kk_174 = buffer.data(kk + 174);
    const auto *kk_175 = buffer.data(kk + 175);
    const auto *kk_176 = buffer.data(kk + 176);
    const auto *kk_177 = buffer.data(kk + 177);
    const auto *kk_178 = buffer.data(kk + 178);
    const auto *kk_179 = buffer.data(kk + 179);
    const auto *kk_216 = buffer.data(kk + 216);
    const auto *kk_217 = buffer.data(kk + 217);
    const auto *kk_218 = buffer.data(kk + 218);
    const auto *kk_219 = buffer.data(kk + 219);
    const auto *kk_220 = buffer.data(kk + 220);
    const auto *kk_221 = buffer.data(kk + 221);
    const auto *kk_222 = buffer.data(kk + 222);
    const auto *kk_223 = buffer.data(kk + 223);
    const auto *kk_224 = buffer.data(kk + 224);
    const auto *kk_225 = buffer.data(kk + 225);
    const auto *kk_226 = buffer.data(kk + 226);
    const auto *kk_227 = buffer.data(kk + 227);
    const auto *kk_228 = buffer.data(kk + 228);
    const auto *kk_229 = buffer.data(kk + 229);
    const auto *kk_230 = buffer.data(kk + 230);
    const auto *kk_231 = buffer.data(kk + 231);
    const auto *kk_232 = buffer.data(kk + 232);
    const auto *kk_233 = buffer.data(kk + 233);
    const auto *kk_234 = buffer.data(kk + 234);
    const auto *kk_235 = buffer.data(kk + 235);
    const auto *kk_236 = buffer.data(kk + 236);
    const auto *kk_237 = buffer.data(kk + 237);
    const auto *kk_238 = buffer.data(kk + 238);
    const auto *kk_239 = buffer.data(kk + 239);
    const auto *kk_240 = buffer.data(kk + 240);
    const auto *kk_241 = buffer.data(kk + 241);
    const auto *kk_242 = buffer.data(kk + 242);
    const auto *kk_243 = buffer.data(kk + 243);
    const auto *kk_244 = buffer.data(kk + 244);
    const auto *kk_245 = buffer.data(kk + 245);
    const auto *kk_246 = buffer.data(kk + 246);
    const auto *kk_247 = buffer.data(kk + 247);
    const auto *kk_248 = buffer.data(kk + 248);
    const auto *kk_249 = buffer.data(kk + 249);
    const auto *kk_250 = buffer.data(kk + 250);
    const auto *kk_251 = buffer.data(kk + 251);
    const auto *kk_252 = buffer.data(kk + 252);
    const auto *kk_253 = buffer.data(kk + 253);
    const auto *kk_254 = buffer.data(kk + 254);
    const auto *kk_255 = buffer.data(kk + 255);
    const auto *kk_256 = buffer.data(kk + 256);
    const auto *kk_257 = buffer.data(kk + 257);
    const auto *kk_258 = buffer.data(kk + 258);
    const auto *kk_259 = buffer.data(kk + 259);
    const auto *kk_260 = buffer.data(kk + 260);
    const auto *kk_261 = buffer.data(kk + 261);
    const auto *kk_262 = buffer.data(kk + 262);
    const auto *kk_263 = buffer.data(kk + 263);
    const auto *kk_264 = buffer.data(kk + 264);
    const auto *kk_265 = buffer.data(kk + 265);
    const auto *kk_266 = buffer.data(kk + 266);
    const auto *kk_267 = buffer.data(kk + 267);
    const auto *kk_268 = buffer.data(kk + 268);
    const auto *kk_269 = buffer.data(kk + 269);
    const auto *kk_270 = buffer.data(kk + 270);
    const auto *kk_271 = buffer.data(kk + 271);
    const auto *kk_272 = buffer.data(kk + 272);
    const auto *kk_273 = buffer.data(kk + 273);
    const auto *kk_274 = buffer.data(kk + 274);
    const auto *kk_275 = buffer.data(kk + 275);
    const auto *kk_276 = buffer.data(kk + 276);
    const auto *kk_277 = buffer.data(kk + 277);
    const auto *kk_278 = buffer.data(kk + 278);
    const auto *kk_279 = buffer.data(kk + 279);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, kk_36, kk_37, kk_38, kk_39, \
                         kk_40, kk_41, kk_42, kk_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kk_36[k];

        t_1[k] = f_0 * kk_37[k];

        t_2[k] = f_0 * kk_38[k];

        t_3[k] = f_0 * kk_39[k];

        t_4[k] = f_0 * kk_40[k];

        t_5[k] = f_0 * kk_41[k];

        t_6[k] = f_0 * kk_42[k];

        t_7[k] = f_0 * kk_43[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, kk_44, kk_45, kk_46, \
                         kk_47, kk_48, kk_49, kk_50, kk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * kk_44[k];

        t_9[k] = f_0 * kk_45[k];

        t_10[k] = f_0 * kk_46[k];

        t_11[k] = f_0 * kk_47[k];

        t_12[k] = f_0 * kk_48[k];

        t_13[k] = f_0 * kk_49[k];

        t_14[k] = f_0 * kk_50[k];

        t_15[k] = f_0 * kk_51[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, kk_52, kk_53, kk_54, \
                         kk_55, kk_56, kk_57, kk_58, kk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * kk_52[k];

        t_17[k] = f_0 * kk_53[k];

        t_18[k] = f_0 * kk_54[k];

        t_19[k] = f_0 * kk_55[k];

        t_20[k] = f_0 * kk_56[k];

        t_21[k] = f_0 * kk_57[k];

        t_22[k] = f_0 * kk_58[k];

        t_23[k] = f_0 * kk_59[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, kk_60, kk_61, kk_62, \
                         kk_63, kk_64, kk_65, kk_66, kk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * kk_60[k];

        t_25[k] = f_0 * kk_61[k];

        t_26[k] = f_0 * kk_62[k];

        t_27[k] = f_0 * kk_63[k];

        t_28[k] = f_0 * kk_64[k];

        t_29[k] = f_0 * kk_65[k];

        t_30[k] = f_0 * kk_66[k];

        t_31[k] = f_0 * kk_67[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, hk_0, hk_1, kk_68, kk_69, kk_70, \
                         kk_71, kk_108, kk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * kk_68[k];

        t_33[k] = f_0 * kk_69[k];

        t_34[k] = f_0 * kk_70[k];

        t_35[k] = f_0 * kk_71[k];

        t_36[k] = -hk_0[k]
                  + f_0 * kk_108[k];

        t_37[k] = -hk_1[k]
                  + f_0 * kk_109[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, hk_2, hk_3, hk_4, hk_5, hk_6, kk_110, \
                         kk_111, kk_112, kk_113, kk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -hk_2[k]
                  + f_0 * kk_110[k];

        t_39[k] = -hk_3[k]
                  + f_0 * kk_111[k];

        t_40[k] = -hk_4[k]
                  + f_0 * kk_112[k];

        t_41[k] = -hk_5[k]
                  + f_0 * kk_113[k];

        t_42[k] = -hk_6[k]
                  + f_0 * kk_114[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, hk_7, hk_8, hk_9, hk_10, hk_11, kk_115, \
                         kk_116, kk_117, kk_118, kk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -hk_7[k]
                  + f_0 * kk_115[k];

        t_44[k] = -hk_8[k]
                  + f_0 * kk_116[k];

        t_45[k] = -hk_9[k]
                  + f_0 * kk_117[k];

        t_46[k] = -hk_10[k]
                  + f_0 * kk_118[k];

        t_47[k] = -hk_11[k]
                  + f_0 * kk_119[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, hk_12, hk_13, hk_14, hk_15, hk_16, \
                         kk_120, kk_121, kk_122, kk_123, kk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -hk_12[k]
                  + f_0 * kk_120[k];

        t_49[k] = -hk_13[k]
                  + f_0 * kk_121[k];

        t_50[k] = -hk_14[k]
                  + f_0 * kk_122[k];

        t_51[k] = -hk_15[k]
                  + f_0 * kk_123[k];

        t_52[k] = -hk_16[k]
                  + f_0 * kk_124[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, hk_17, hk_18, hk_19, hk_20, hk_21, \
                         kk_125, kk_126, kk_127, kk_128, kk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -hk_17[k]
                  + f_0 * kk_125[k];

        t_54[k] = -hk_18[k]
                  + f_0 * kk_126[k];

        t_55[k] = -hk_19[k]
                  + f_0 * kk_127[k];

        t_56[k] = -hk_20[k]
                  + f_0 * kk_128[k];

        t_57[k] = -hk_21[k]
                  + f_0 * kk_129[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, hk_22, hk_23, hk_24, hk_25, hk_26, \
                         kk_130, kk_131, kk_132, kk_133, kk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -hk_22[k]
                  + f_0 * kk_130[k];

        t_59[k] = -hk_23[k]
                  + f_0 * kk_131[k];

        t_60[k] = -hk_24[k]
                  + f_0 * kk_132[k];

        t_61[k] = -hk_25[k]
                  + f_0 * kk_133[k];

        t_62[k] = -hk_26[k]
                  + f_0 * kk_134[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, hk_27, hk_28, hk_29, hk_30, hk_31, \
                         kk_135, kk_136, kk_137, kk_138, kk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -hk_27[k]
                  + f_0 * kk_135[k];

        t_64[k] = -hk_28[k]
                  + f_0 * kk_136[k];

        t_65[k] = -hk_29[k]
                  + f_0 * kk_137[k];

        t_66[k] = -hk_30[k]
                  + f_0 * kk_138[k];

        t_67[k] = -hk_31[k]
                  + f_0 * kk_139[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, hk_32, hk_33, hk_34, hk_35, \
                         kk_140, kk_141, kk_142, kk_143, kk_144, \
                         kk_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -hk_32[k]
                  + f_0 * kk_140[k];

        t_69[k] = -hk_33[k]
                  + f_0 * kk_141[k];

        t_70[k] = -hk_34[k]
                  + f_0 * kk_142[k];

        t_71[k] = -hk_35[k]
                  + f_0 * kk_143[k];

        t_72[k] = f_0 * kk_144[k];

        t_73[k] = f_0 * kk_145[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, t_80, t_81, kk_146, kk_147, \
                         kk_148, kk_149, kk_150, kk_151, kk_152, \
                         kk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * kk_146[k];

        t_75[k] = f_0 * kk_147[k];

        t_76[k] = f_0 * kk_148[k];

        t_77[k] = f_0 * kk_149[k];

        t_78[k] = f_0 * kk_150[k];

        t_79[k] = f_0 * kk_151[k];

        t_80[k] = f_0 * kk_152[k];

        t_81[k] = f_0 * kk_153[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, t_87, t_88, t_89, kk_154, kk_155, \
                         kk_156, kk_157, kk_158, kk_159, kk_160, \
                         kk_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * kk_154[k];

        t_83[k] = f_0 * kk_155[k];

        t_84[k] = f_0 * kk_156[k];

        t_85[k] = f_0 * kk_157[k];

        t_86[k] = f_0 * kk_158[k];

        t_87[k] = f_0 * kk_159[k];

        t_88[k] = f_0 * kk_160[k];

        t_89[k] = f_0 * kk_161[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, kk_162, kk_163, \
                         kk_164, kk_165, kk_166, kk_167, kk_168, \
                         kk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * kk_162[k];

        t_91[k] = f_0 * kk_163[k];

        t_92[k] = f_0 * kk_164[k];

        t_93[k] = f_0 * kk_165[k];

        t_94[k] = f_0 * kk_166[k];

        t_95[k] = f_0 * kk_167[k];

        t_96[k] = f_0 * kk_168[k];

        t_97[k] = f_0 * kk_169[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, t_104, t_105, kk_170, kk_171, \
                         kk_172, kk_173, kk_174, kk_175, kk_176, \
                         kk_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * kk_170[k];

        t_99[k] = f_0 * kk_171[k];

        t_100[k] = f_0 * kk_172[k];

        t_101[k] = f_0 * kk_173[k];

        t_102[k] = f_0 * kk_174[k];

        t_103[k] = f_0 * kk_175[k];

        t_104[k] = f_0 * kk_176[k];

        t_105[k] = f_0 * kk_177[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, hk_36, hk_37, hk_38, hk_39, \
                         kk_178, kk_179, kk_216, kk_217, kk_218, \
                         kk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * kk_178[k];

        t_107[k] = f_0 * kk_179[k];

        t_108[k] = -2.0 * hk_36[k]
                   + f_0 * kk_216[k];

        t_109[k] = -2.0 * hk_37[k]
                   + f_0 * kk_217[k];

        t_110[k] = -2.0 * hk_38[k]
                   + f_0 * kk_218[k];

        t_111[k] = -2.0 * hk_39[k]
                   + f_0 * kk_219[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, hk_40, hk_41, hk_42, hk_43, hk_44, \
                         kk_220, kk_221, kk_222, kk_223, kk_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -2.0 * hk_40[k]
                   + f_0 * kk_220[k];

        t_113[k] = -2.0 * hk_41[k]
                   + f_0 * kk_221[k];

        t_114[k] = -2.0 * hk_42[k]
                   + f_0 * kk_222[k];

        t_115[k] = -2.0 * hk_43[k]
                   + f_0 * kk_223[k];

        t_116[k] = -2.0 * hk_44[k]
                   + f_0 * kk_224[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, hk_45, hk_46, hk_47, hk_48, hk_49, \
                         kk_225, kk_226, kk_227, kk_228, kk_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -2.0 * hk_45[k]
                   + f_0 * kk_225[k];

        t_118[k] = -2.0 * hk_46[k]
                   + f_0 * kk_226[k];

        t_119[k] = -2.0 * hk_47[k]
                   + f_0 * kk_227[k];

        t_120[k] = -2.0 * hk_48[k]
                   + f_0 * kk_228[k];

        t_121[k] = -2.0 * hk_49[k]
                   + f_0 * kk_229[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, hk_50, hk_51, hk_52, hk_53, hk_54, \
                         kk_230, kk_231, kk_232, kk_233, kk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * hk_50[k]
                   + f_0 * kk_230[k];

        t_123[k] = -2.0 * hk_51[k]
                   + f_0 * kk_231[k];

        t_124[k] = -2.0 * hk_52[k]
                   + f_0 * kk_232[k];

        t_125[k] = -2.0 * hk_53[k]
                   + f_0 * kk_233[k];

        t_126[k] = -2.0 * hk_54[k]
                   + f_0 * kk_234[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, hk_55, hk_56, hk_57, hk_58, hk_59, \
                         kk_235, kk_236, kk_237, kk_238, kk_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * hk_55[k]
                   + f_0 * kk_235[k];

        t_128[k] = -2.0 * hk_56[k]
                   + f_0 * kk_236[k];

        t_129[k] = -2.0 * hk_57[k]
                   + f_0 * kk_237[k];

        t_130[k] = -2.0 * hk_58[k]
                   + f_0 * kk_238[k];

        t_131[k] = -2.0 * hk_59[k]
                   + f_0 * kk_239[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, hk_60, hk_61, hk_62, hk_63, hk_64, \
                         kk_240, kk_241, kk_242, kk_243, kk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * hk_60[k]
                   + f_0 * kk_240[k];

        t_133[k] = -2.0 * hk_61[k]
                   + f_0 * kk_241[k];

        t_134[k] = -2.0 * hk_62[k]
                   + f_0 * kk_242[k];

        t_135[k] = -2.0 * hk_63[k]
                   + f_0 * kk_243[k];

        t_136[k] = -2.0 * hk_64[k]
                   + f_0 * kk_244[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, hk_65, hk_66, hk_67, hk_68, hk_69, \
                         kk_245, kk_246, kk_247, kk_248, kk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -2.0 * hk_65[k]
                   + f_0 * kk_245[k];

        t_138[k] = -2.0 * hk_66[k]
                   + f_0 * kk_246[k];

        t_139[k] = -2.0 * hk_67[k]
                   + f_0 * kk_247[k];

        t_140[k] = -2.0 * hk_68[k]
                   + f_0 * kk_248[k];

        t_141[k] = -2.0 * hk_69[k]
                   + f_0 * kk_249[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, hk_70, hk_71, hk_72, hk_73, hk_74, \
                         kk_250, kk_251, kk_252, kk_253, kk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -2.0 * hk_70[k]
                   + f_0 * kk_250[k];

        t_143[k] = -2.0 * hk_71[k]
                   + f_0 * kk_251[k];

        t_144[k] = -hk_72[k]
                   + f_0 * kk_252[k];

        t_145[k] = -hk_73[k]
                   + f_0 * kk_253[k];

        t_146[k] = -hk_74[k]
                   + f_0 * kk_254[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, hk_75, hk_76, hk_77, hk_78, hk_79, \
                         kk_255, kk_256, kk_257, kk_258, kk_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -hk_75[k]
                   + f_0 * kk_255[k];

        t_148[k] = -hk_76[k]
                   + f_0 * kk_256[k];

        t_149[k] = -hk_77[k]
                   + f_0 * kk_257[k];

        t_150[k] = -hk_78[k]
                   + f_0 * kk_258[k];

        t_151[k] = -hk_79[k]
                   + f_0 * kk_259[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, hk_80, hk_81, hk_82, hk_83, hk_84, \
                         kk_260, kk_261, kk_262, kk_263, kk_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -hk_80[k]
                   + f_0 * kk_260[k];

        t_153[k] = -hk_81[k]
                   + f_0 * kk_261[k];

        t_154[k] = -hk_82[k]
                   + f_0 * kk_262[k];

        t_155[k] = -hk_83[k]
                   + f_0 * kk_263[k];

        t_156[k] = -hk_84[k]
                   + f_0 * kk_264[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, hk_85, hk_86, hk_87, hk_88, hk_89, \
                         kk_265, kk_266, kk_267, kk_268, kk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -hk_85[k]
                   + f_0 * kk_265[k];

        t_158[k] = -hk_86[k]
                   + f_0 * kk_266[k];

        t_159[k] = -hk_87[k]
                   + f_0 * kk_267[k];

        t_160[k] = -hk_88[k]
                   + f_0 * kk_268[k];

        t_161[k] = -hk_89[k]
                   + f_0 * kk_269[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, hk_90, hk_91, hk_92, hk_93, hk_94, \
                         kk_270, kk_271, kk_272, kk_273, kk_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -hk_90[k]
                   + f_0 * kk_270[k];

        t_163[k] = -hk_91[k]
                   + f_0 * kk_271[k];

        t_164[k] = -hk_92[k]
                   + f_0 * kk_272[k];

        t_165[k] = -hk_93[k]
                   + f_0 * kk_273[k];

        t_166[k] = -hk_94[k]
                   + f_0 * kk_274[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, hk_95, hk_96, hk_97, hk_98, hk_99, \
                         kk_275, kk_276, kk_277, kk_278, kk_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -hk_95[k]
                   + f_0 * kk_275[k];

        t_168[k] = -hk_96[k]
                   + f_0 * kk_276[k];

        t_169[k] = -hk_97[k]
                   + f_0 * kk_277[k];

        t_170[k] = -hk_98[k]
                   + f_0 * kk_278[k];

        t_171[k] = -hk_99[k]
                   + f_0 * kk_279[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_100 = buffer.data(hk + 100);
    const auto *hk_101 = buffer.data(hk + 101);
    const auto *hk_102 = buffer.data(hk + 102);
    const auto *hk_103 = buffer.data(hk + 103);
    const auto *hk_104 = buffer.data(hk + 104);
    const auto *hk_105 = buffer.data(hk + 105);
    const auto *hk_106 = buffer.data(hk + 106);
    const auto *hk_107 = buffer.data(hk + 107);
    const auto *hk_108 = buffer.data(hk + 108);
    const auto *hk_109 = buffer.data(hk + 109);
    const auto *hk_110 = buffer.data(hk + 110);
    const auto *hk_111 = buffer.data(hk + 111);
    const auto *hk_112 = buffer.data(hk + 112);
    const auto *hk_113 = buffer.data(hk + 113);
    const auto *hk_114 = buffer.data(hk + 114);
    const auto *hk_115 = buffer.data(hk + 115);
    const auto *hk_116 = buffer.data(hk + 116);
    const auto *hk_117 = buffer.data(hk + 117);
    const auto *hk_118 = buffer.data(hk + 118);
    const auto *hk_119 = buffer.data(hk + 119);
    const auto *hk_120 = buffer.data(hk + 120);
    const auto *hk_121 = buffer.data(hk + 121);
    const auto *hk_122 = buffer.data(hk + 122);
    const auto *hk_123 = buffer.data(hk + 123);
    const auto *hk_124 = buffer.data(hk + 124);
    const auto *hk_125 = buffer.data(hk + 125);
    const auto *hk_126 = buffer.data(hk + 126);
    const auto *hk_127 = buffer.data(hk + 127);
    const auto *hk_128 = buffer.data(hk + 128);
    const auto *hk_129 = buffer.data(hk + 129);
    const auto *hk_130 = buffer.data(hk + 130);
    const auto *hk_131 = buffer.data(hk + 131);
    const auto *hk_132 = buffer.data(hk + 132);
    const auto *hk_133 = buffer.data(hk + 133);
    const auto *hk_134 = buffer.data(hk + 134);
    const auto *hk_135 = buffer.data(hk + 135);
    const auto *hk_136 = buffer.data(hk + 136);
    const auto *hk_137 = buffer.data(hk + 137);
    const auto *hk_138 = buffer.data(hk + 138);
    const auto *hk_139 = buffer.data(hk + 139);
    const auto *hk_140 = buffer.data(hk + 140);
    const auto *hk_141 = buffer.data(hk + 141);
    const auto *hk_142 = buffer.data(hk + 142);
    const auto *hk_143 = buffer.data(hk + 143);
    const auto *hk_144 = buffer.data(hk + 144);
    const auto *hk_145 = buffer.data(hk + 145);
    const auto *hk_146 = buffer.data(hk + 146);
    const auto *hk_147 = buffer.data(hk + 147);
    const auto *hk_148 = buffer.data(hk + 148);
    const auto *hk_149 = buffer.data(hk + 149);
    const auto *hk_150 = buffer.data(hk + 150);
    const auto *hk_151 = buffer.data(hk + 151);
    const auto *hk_152 = buffer.data(hk + 152);
    const auto *hk_153 = buffer.data(hk + 153);
    const auto *hk_154 = buffer.data(hk + 154);
    const auto *hk_155 = buffer.data(hk + 155);
    const auto *hk_156 = buffer.data(hk + 156);
    const auto *hk_157 = buffer.data(hk + 157);
    const auto *hk_158 = buffer.data(hk + 158);
    const auto *hk_159 = buffer.data(hk + 159);
    const auto *hk_160 = buffer.data(hk + 160);
    const auto *hk_161 = buffer.data(hk + 161);
    const auto *hk_162 = buffer.data(hk + 162);
    const auto *hk_163 = buffer.data(hk + 163);
    const auto *hk_164 = buffer.data(hk + 164);
    const auto *hk_165 = buffer.data(hk + 165);
    const auto *hk_166 = buffer.data(hk + 166);
    const auto *hk_167 = buffer.data(hk + 167);
    const auto *hk_168 = buffer.data(hk + 168);
    const auto *hk_169 = buffer.data(hk + 169);
    const auto *hk_170 = buffer.data(hk + 170);
    const auto *hk_171 = buffer.data(hk + 171);
    const auto *hk_172 = buffer.data(hk + 172);
    const auto *hk_173 = buffer.data(hk + 173);
    const auto *hk_174 = buffer.data(hk + 174);
    const auto *hk_175 = buffer.data(hk + 175);
    const auto *hk_176 = buffer.data(hk + 176);
    const auto *hk_177 = buffer.data(hk + 177);
    const auto *hk_178 = buffer.data(hk + 178);
    const auto *hk_179 = buffer.data(hk + 179);
    const auto *hk_180 = buffer.data(hk + 180);
    const auto *hk_181 = buffer.data(hk + 181);
    const auto *hk_182 = buffer.data(hk + 182);
    const auto *hk_183 = buffer.data(hk + 183);
    const auto *hk_184 = buffer.data(hk + 184);
    const auto *hk_185 = buffer.data(hk + 185);
    const auto *hk_186 = buffer.data(hk + 186);
    const auto *hk_187 = buffer.data(hk + 187);
    const auto *hk_188 = buffer.data(hk + 188);
    const auto *hk_189 = buffer.data(hk + 189);
    const auto *hk_190 = buffer.data(hk + 190);
    const auto *hk_191 = buffer.data(hk + 191);
    const auto *hk_192 = buffer.data(hk + 192);
    const auto *hk_193 = buffer.data(hk + 193);
    const auto *hk_194 = buffer.data(hk + 194);
    const auto *hk_195 = buffer.data(hk + 195);
    const auto *hk_196 = buffer.data(hk + 196);
    const auto *hk_197 = buffer.data(hk + 197);
    const auto *hk_198 = buffer.data(hk + 198);
    const auto *hk_199 = buffer.data(hk + 199);
    const auto *hk_200 = buffer.data(hk + 200);
    const auto *hk_201 = buffer.data(hk + 201);
    const auto *hk_202 = buffer.data(hk + 202);
    const auto *hk_203 = buffer.data(hk + 203);
    const auto *hk_204 = buffer.data(hk + 204);
    const auto *hk_205 = buffer.data(hk + 205);
    const auto *hk_206 = buffer.data(hk + 206);
    const auto *hk_207 = buffer.data(hk + 207);
    const auto *hk_208 = buffer.data(hk + 208);
    const auto *hk_209 = buffer.data(hk + 209);
    const auto *hk_210 = buffer.data(hk + 210);
    const auto *hk_211 = buffer.data(hk + 211);
    const auto *hk_212 = buffer.data(hk + 212);
    const auto *hk_213 = buffer.data(hk + 213);
    const auto *hk_214 = buffer.data(hk + 214);
    const auto *hk_215 = buffer.data(hk + 215);

    const auto *kk_280 = buffer.data(kk + 280);
    const auto *kk_281 = buffer.data(kk + 281);
    const auto *kk_282 = buffer.data(kk + 282);
    const auto *kk_283 = buffer.data(kk + 283);
    const auto *kk_284 = buffer.data(kk + 284);
    const auto *kk_285 = buffer.data(kk + 285);
    const auto *kk_286 = buffer.data(kk + 286);
    const auto *kk_287 = buffer.data(kk + 287);
    const auto *kk_288 = buffer.data(kk + 288);
    const auto *kk_289 = buffer.data(kk + 289);
    const auto *kk_290 = buffer.data(kk + 290);
    const auto *kk_291 = buffer.data(kk + 291);
    const auto *kk_292 = buffer.data(kk + 292);
    const auto *kk_293 = buffer.data(kk + 293);
    const auto *kk_294 = buffer.data(kk + 294);
    const auto *kk_295 = buffer.data(kk + 295);
    const auto *kk_296 = buffer.data(kk + 296);
    const auto *kk_297 = buffer.data(kk + 297);
    const auto *kk_298 = buffer.data(kk + 298);
    const auto *kk_299 = buffer.data(kk + 299);
    const auto *kk_300 = buffer.data(kk + 300);
    const auto *kk_301 = buffer.data(kk + 301);
    const auto *kk_302 = buffer.data(kk + 302);
    const auto *kk_303 = buffer.data(kk + 303);
    const auto *kk_304 = buffer.data(kk + 304);
    const auto *kk_305 = buffer.data(kk + 305);
    const auto *kk_306 = buffer.data(kk + 306);
    const auto *kk_307 = buffer.data(kk + 307);
    const auto *kk_308 = buffer.data(kk + 308);
    const auto *kk_309 = buffer.data(kk + 309);
    const auto *kk_310 = buffer.data(kk + 310);
    const auto *kk_311 = buffer.data(kk + 311);
    const auto *kk_312 = buffer.data(kk + 312);
    const auto *kk_313 = buffer.data(kk + 313);
    const auto *kk_314 = buffer.data(kk + 314);
    const auto *kk_315 = buffer.data(kk + 315);
    const auto *kk_316 = buffer.data(kk + 316);
    const auto *kk_317 = buffer.data(kk + 317);
    const auto *kk_318 = buffer.data(kk + 318);
    const auto *kk_319 = buffer.data(kk + 319);
    const auto *kk_320 = buffer.data(kk + 320);
    const auto *kk_321 = buffer.data(kk + 321);
    const auto *kk_322 = buffer.data(kk + 322);
    const auto *kk_323 = buffer.data(kk + 323);
    const auto *kk_360 = buffer.data(kk + 360);
    const auto *kk_361 = buffer.data(kk + 361);
    const auto *kk_362 = buffer.data(kk + 362);
    const auto *kk_363 = buffer.data(kk + 363);
    const auto *kk_364 = buffer.data(kk + 364);
    const auto *kk_365 = buffer.data(kk + 365);
    const auto *kk_366 = buffer.data(kk + 366);
    const auto *kk_367 = buffer.data(kk + 367);
    const auto *kk_368 = buffer.data(kk + 368);
    const auto *kk_369 = buffer.data(kk + 369);
    const auto *kk_370 = buffer.data(kk + 370);
    const auto *kk_371 = buffer.data(kk + 371);
    const auto *kk_372 = buffer.data(kk + 372);
    const auto *kk_373 = buffer.data(kk + 373);
    const auto *kk_374 = buffer.data(kk + 374);
    const auto *kk_375 = buffer.data(kk + 375);
    const auto *kk_376 = buffer.data(kk + 376);
    const auto *kk_377 = buffer.data(kk + 377);
    const auto *kk_378 = buffer.data(kk + 378);
    const auto *kk_379 = buffer.data(kk + 379);
    const auto *kk_380 = buffer.data(kk + 380);
    const auto *kk_381 = buffer.data(kk + 381);
    const auto *kk_382 = buffer.data(kk + 382);
    const auto *kk_383 = buffer.data(kk + 383);
    const auto *kk_384 = buffer.data(kk + 384);
    const auto *kk_385 = buffer.data(kk + 385);
    const auto *kk_386 = buffer.data(kk + 386);
    const auto *kk_387 = buffer.data(kk + 387);
    const auto *kk_388 = buffer.data(kk + 388);
    const auto *kk_389 = buffer.data(kk + 389);
    const auto *kk_390 = buffer.data(kk + 390);
    const auto *kk_391 = buffer.data(kk + 391);
    const auto *kk_392 = buffer.data(kk + 392);
    const auto *kk_393 = buffer.data(kk + 393);
    const auto *kk_394 = buffer.data(kk + 394);
    const auto *kk_395 = buffer.data(kk + 395);
    const auto *kk_396 = buffer.data(kk + 396);
    const auto *kk_397 = buffer.data(kk + 397);
    const auto *kk_398 = buffer.data(kk + 398);
    const auto *kk_399 = buffer.data(kk + 399);
    const auto *kk_400 = buffer.data(kk + 400);
    const auto *kk_401 = buffer.data(kk + 401);
    const auto *kk_402 = buffer.data(kk + 402);
    const auto *kk_403 = buffer.data(kk + 403);
    const auto *kk_404 = buffer.data(kk + 404);
    const auto *kk_405 = buffer.data(kk + 405);
    const auto *kk_406 = buffer.data(kk + 406);
    const auto *kk_407 = buffer.data(kk + 407);
    const auto *kk_408 = buffer.data(kk + 408);
    const auto *kk_409 = buffer.data(kk + 409);
    const auto *kk_410 = buffer.data(kk + 410);
    const auto *kk_411 = buffer.data(kk + 411);
    const auto *kk_412 = buffer.data(kk + 412);
    const auto *kk_413 = buffer.data(kk + 413);
    const auto *kk_414 = buffer.data(kk + 414);
    const auto *kk_415 = buffer.data(kk + 415);
    const auto *kk_416 = buffer.data(kk + 416);
    const auto *kk_417 = buffer.data(kk + 417);
    const auto *kk_418 = buffer.data(kk + 418);
    const auto *kk_419 = buffer.data(kk + 419);
    const auto *kk_420 = buffer.data(kk + 420);
    const auto *kk_421 = buffer.data(kk + 421);
    const auto *kk_422 = buffer.data(kk + 422);
    const auto *kk_423 = buffer.data(kk + 423);
    const auto *kk_424 = buffer.data(kk + 424);
    const auto *kk_425 = buffer.data(kk + 425);
    const auto *kk_426 = buffer.data(kk + 426);
    const auto *kk_427 = buffer.data(kk + 427);
    const auto *kk_428 = buffer.data(kk + 428);
    const auto *kk_429 = buffer.data(kk + 429);
    const auto *kk_430 = buffer.data(kk + 430);
    const auto *kk_431 = buffer.data(kk + 431);
    const auto *kk_432 = buffer.data(kk + 432);
    const auto *kk_433 = buffer.data(kk + 433);
    const auto *kk_434 = buffer.data(kk + 434);
    const auto *kk_435 = buffer.data(kk + 435);
    const auto *kk_436 = buffer.data(kk + 436);
    const auto *kk_437 = buffer.data(kk + 437);
    const auto *kk_438 = buffer.data(kk + 438);
    const auto *kk_439 = buffer.data(kk + 439);
    const auto *kk_440 = buffer.data(kk + 440);
    const auto *kk_441 = buffer.data(kk + 441);
    const auto *kk_442 = buffer.data(kk + 442);
    const auto *kk_443 = buffer.data(kk + 443);
    const auto *kk_444 = buffer.data(kk + 444);
    const auto *kk_445 = buffer.data(kk + 445);
    const auto *kk_446 = buffer.data(kk + 446);
    const auto *kk_447 = buffer.data(kk + 447);
    const auto *kk_448 = buffer.data(kk + 448);
    const auto *kk_449 = buffer.data(kk + 449);
    const auto *kk_450 = buffer.data(kk + 450);
    const auto *kk_451 = buffer.data(kk + 451);
    const auto *kk_452 = buffer.data(kk + 452);
    const auto *kk_453 = buffer.data(kk + 453);
    const auto *kk_454 = buffer.data(kk + 454);
    const auto *kk_455 = buffer.data(kk + 455);
    const auto *kk_456 = buffer.data(kk + 456);
    const auto *kk_457 = buffer.data(kk + 457);
    const auto *kk_458 = buffer.data(kk + 458);
    const auto *kk_459 = buffer.data(kk + 459);
    const auto *kk_460 = buffer.data(kk + 460);
    const auto *kk_461 = buffer.data(kk + 461);
    const auto *kk_462 = buffer.data(kk + 462);
    const auto *kk_463 = buffer.data(kk + 463);
    const auto *kk_464 = buffer.data(kk + 464);
    const auto *kk_465 = buffer.data(kk + 465);
    const auto *kk_466 = buffer.data(kk + 466);
    const auto *kk_467 = buffer.data(kk + 467);
    const auto *kk_468 = buffer.data(kk + 468);
    const auto *kk_469 = buffer.data(kk + 469);
    const auto *kk_470 = buffer.data(kk + 470);
    const auto *kk_471 = buffer.data(kk + 471);
    const auto *kk_472 = buffer.data(kk + 472);
    const auto *kk_473 = buffer.data(kk + 473);
    const auto *kk_474 = buffer.data(kk + 474);
    const auto *kk_475 = buffer.data(kk + 475);
    const auto *kk_476 = buffer.data(kk + 476);
    const auto *kk_477 = buffer.data(kk + 477);

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, hk_100, hk_101, hk_102, hk_103, \
                         hk_104, kk_280, kk_281, kk_282, kk_283, \
                         kk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -hk_100[k]
                   + f_0 * kk_280[k];

        t_173[k] = -hk_101[k]
                   + f_0 * kk_281[k];

        t_174[k] = -hk_102[k]
                   + f_0 * kk_282[k];

        t_175[k] = -hk_103[k]
                   + f_0 * kk_283[k];

        t_176[k] = -hk_104[k]
                   + f_0 * kk_284[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, t_182, hk_105, hk_106, hk_107, \
                         kk_285, kk_286, kk_287, kk_288, kk_289, \
                         kk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -hk_105[k]
                   + f_0 * kk_285[k];

        t_178[k] = -hk_106[k]
                   + f_0 * kk_286[k];

        t_179[k] = -hk_107[k]
                   + f_0 * kk_287[k];

        t_180[k] = f_0 * kk_288[k];

        t_181[k] = f_0 * kk_289[k];

        t_182[k] = f_0 * kk_290[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, t_188, t_189, t_190, kk_291, \
                         kk_292, kk_293, kk_294, kk_295, kk_296, kk_297, \
                         kk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_0 * kk_291[k];

        t_184[k] = f_0 * kk_292[k];

        t_185[k] = f_0 * kk_293[k];

        t_186[k] = f_0 * kk_294[k];

        t_187[k] = f_0 * kk_295[k];

        t_188[k] = f_0 * kk_296[k];

        t_189[k] = f_0 * kk_297[k];

        t_190[k] = f_0 * kk_298[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, t_197, t_198, kk_299, \
                         kk_300, kk_301, kk_302, kk_303, kk_304, kk_305, \
                         kk_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_0 * kk_299[k];

        t_192[k] = f_0 * kk_300[k];

        t_193[k] = f_0 * kk_301[k];

        t_194[k] = f_0 * kk_302[k];

        t_195[k] = f_0 * kk_303[k];

        t_196[k] = f_0 * kk_304[k];

        t_197[k] = f_0 * kk_305[k];

        t_198[k] = f_0 * kk_306[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, t_204, t_205, t_206, kk_307, \
                         kk_308, kk_309, kk_310, kk_311, kk_312, kk_313, \
                         kk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_0 * kk_307[k];

        t_200[k] = f_0 * kk_308[k];

        t_201[k] = f_0 * kk_309[k];

        t_202[k] = f_0 * kk_310[k];

        t_203[k] = f_0 * kk_311[k];

        t_204[k] = f_0 * kk_312[k];

        t_205[k] = f_0 * kk_313[k];

        t_206[k] = f_0 * kk_314[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, t_213, t_214, kk_315, \
                         kk_316, kk_317, kk_318, kk_319, kk_320, kk_321, \
                         kk_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_0 * kk_315[k];

        t_208[k] = f_0 * kk_316[k];

        t_209[k] = f_0 * kk_317[k];

        t_210[k] = f_0 * kk_318[k];

        t_211[k] = f_0 * kk_319[k];

        t_212[k] = f_0 * kk_320[k];

        t_213[k] = f_0 * kk_321[k];

        t_214[k] = f_0 * kk_322[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, hk_108, hk_109, hk_110, hk_111, \
                         kk_323, kk_360, kk_361, kk_362, kk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_0 * kk_323[k];

        t_216[k] = -3.0 * hk_108[k]
                   + f_0 * kk_360[k];

        t_217[k] = -3.0 * hk_109[k]
                   + f_0 * kk_361[k];

        t_218[k] = -3.0 * hk_110[k]
                   + f_0 * kk_362[k];

        t_219[k] = -3.0 * hk_111[k]
                   + f_0 * kk_363[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, hk_112, hk_113, hk_114, hk_115, \
                         hk_116, kk_364, kk_365, kk_366, kk_367, \
                         kk_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -3.0 * hk_112[k]
                   + f_0 * kk_364[k];

        t_221[k] = -3.0 * hk_113[k]
                   + f_0 * kk_365[k];

        t_222[k] = -3.0 * hk_114[k]
                   + f_0 * kk_366[k];

        t_223[k] = -3.0 * hk_115[k]
                   + f_0 * kk_367[k];

        t_224[k] = -3.0 * hk_116[k]
                   + f_0 * kk_368[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, hk_117, hk_118, hk_119, hk_120, \
                         hk_121, kk_369, kk_370, kk_371, kk_372, \
                         kk_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -3.0 * hk_117[k]
                   + f_0 * kk_369[k];

        t_226[k] = -3.0 * hk_118[k]
                   + f_0 * kk_370[k];

        t_227[k] = -3.0 * hk_119[k]
                   + f_0 * kk_371[k];

        t_228[k] = -3.0 * hk_120[k]
                   + f_0 * kk_372[k];

        t_229[k] = -3.0 * hk_121[k]
                   + f_0 * kk_373[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, hk_122, hk_123, hk_124, hk_125, \
                         hk_126, kk_374, kk_375, kk_376, kk_377, \
                         kk_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -3.0 * hk_122[k]
                   + f_0 * kk_374[k];

        t_231[k] = -3.0 * hk_123[k]
                   + f_0 * kk_375[k];

        t_232[k] = -3.0 * hk_124[k]
                   + f_0 * kk_376[k];

        t_233[k] = -3.0 * hk_125[k]
                   + f_0 * kk_377[k];

        t_234[k] = -3.0 * hk_126[k]
                   + f_0 * kk_378[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, hk_127, hk_128, hk_129, hk_130, \
                         hk_131, kk_379, kk_380, kk_381, kk_382, \
                         kk_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -3.0 * hk_127[k]
                   + f_0 * kk_379[k];

        t_236[k] = -3.0 * hk_128[k]
                   + f_0 * kk_380[k];

        t_237[k] = -3.0 * hk_129[k]
                   + f_0 * kk_381[k];

        t_238[k] = -3.0 * hk_130[k]
                   + f_0 * kk_382[k];

        t_239[k] = -3.0 * hk_131[k]
                   + f_0 * kk_383[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, hk_132, hk_133, hk_134, hk_135, \
                         hk_136, kk_384, kk_385, kk_386, kk_387, \
                         kk_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * hk_132[k]
                   + f_0 * kk_384[k];

        t_241[k] = -3.0 * hk_133[k]
                   + f_0 * kk_385[k];

        t_242[k] = -3.0 * hk_134[k]
                   + f_0 * kk_386[k];

        t_243[k] = -3.0 * hk_135[k]
                   + f_0 * kk_387[k];

        t_244[k] = -3.0 * hk_136[k]
                   + f_0 * kk_388[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, hk_137, hk_138, hk_139, hk_140, \
                         hk_141, kk_389, kk_390, kk_391, kk_392, \
                         kk_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * hk_137[k]
                   + f_0 * kk_389[k];

        t_246[k] = -3.0 * hk_138[k]
                   + f_0 * kk_390[k];

        t_247[k] = -3.0 * hk_139[k]
                   + f_0 * kk_391[k];

        t_248[k] = -3.0 * hk_140[k]
                   + f_0 * kk_392[k];

        t_249[k] = -3.0 * hk_141[k]
                   + f_0 * kk_393[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, hk_142, hk_143, hk_144, hk_145, \
                         hk_146, kk_394, kk_395, kk_396, kk_397, \
                         kk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -3.0 * hk_142[k]
                   + f_0 * kk_394[k];

        t_251[k] = -3.0 * hk_143[k]
                   + f_0 * kk_395[k];

        t_252[k] = -2.0 * hk_144[k]
                   + f_0 * kk_396[k];

        t_253[k] = -2.0 * hk_145[k]
                   + f_0 * kk_397[k];

        t_254[k] = -2.0 * hk_146[k]
                   + f_0 * kk_398[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, hk_147, hk_148, hk_149, hk_150, \
                         hk_151, kk_399, kk_400, kk_401, kk_402, \
                         kk_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * hk_147[k]
                   + f_0 * kk_399[k];

        t_256[k] = -2.0 * hk_148[k]
                   + f_0 * kk_400[k];

        t_257[k] = -2.0 * hk_149[k]
                   + f_0 * kk_401[k];

        t_258[k] = -2.0 * hk_150[k]
                   + f_0 * kk_402[k];

        t_259[k] = -2.0 * hk_151[k]
                   + f_0 * kk_403[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, hk_152, hk_153, hk_154, hk_155, \
                         hk_156, kk_404, kk_405, kk_406, kk_407, \
                         kk_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * hk_152[k]
                   + f_0 * kk_404[k];

        t_261[k] = -2.0 * hk_153[k]
                   + f_0 * kk_405[k];

        t_262[k] = -2.0 * hk_154[k]
                   + f_0 * kk_406[k];

        t_263[k] = -2.0 * hk_155[k]
                   + f_0 * kk_407[k];

        t_264[k] = -2.0 * hk_156[k]
                   + f_0 * kk_408[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, hk_157, hk_158, hk_159, hk_160, \
                         hk_161, kk_409, kk_410, kk_411, kk_412, \
                         kk_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * hk_157[k]
                   + f_0 * kk_409[k];

        t_266[k] = -2.0 * hk_158[k]
                   + f_0 * kk_410[k];

        t_267[k] = -2.0 * hk_159[k]
                   + f_0 * kk_411[k];

        t_268[k] = -2.0 * hk_160[k]
                   + f_0 * kk_412[k];

        t_269[k] = -2.0 * hk_161[k]
                   + f_0 * kk_413[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, hk_162, hk_163, hk_164, hk_165, \
                         hk_166, kk_414, kk_415, kk_416, kk_417, \
                         kk_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * hk_162[k]
                   + f_0 * kk_414[k];

        t_271[k] = -2.0 * hk_163[k]
                   + f_0 * kk_415[k];

        t_272[k] = -2.0 * hk_164[k]
                   + f_0 * kk_416[k];

        t_273[k] = -2.0 * hk_165[k]
                   + f_0 * kk_417[k];

        t_274[k] = -2.0 * hk_166[k]
                   + f_0 * kk_418[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, hk_167, hk_168, hk_169, hk_170, \
                         hk_171, kk_419, kk_420, kk_421, kk_422, \
                         kk_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * hk_167[k]
                   + f_0 * kk_419[k];

        t_276[k] = -2.0 * hk_168[k]
                   + f_0 * kk_420[k];

        t_277[k] = -2.0 * hk_169[k]
                   + f_0 * kk_421[k];

        t_278[k] = -2.0 * hk_170[k]
                   + f_0 * kk_422[k];

        t_279[k] = -2.0 * hk_171[k]
                   + f_0 * kk_423[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, hk_172, hk_173, hk_174, hk_175, \
                         hk_176, kk_424, kk_425, kk_426, kk_427, \
                         kk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * hk_172[k]
                   + f_0 * kk_424[k];

        t_281[k] = -2.0 * hk_173[k]
                   + f_0 * kk_425[k];

        t_282[k] = -2.0 * hk_174[k]
                   + f_0 * kk_426[k];

        t_283[k] = -2.0 * hk_175[k]
                   + f_0 * kk_427[k];

        t_284[k] = -2.0 * hk_176[k]
                   + f_0 * kk_428[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, hk_177, hk_178, hk_179, hk_180, \
                         hk_181, kk_429, kk_430, kk_431, kk_432, \
                         kk_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -2.0 * hk_177[k]
                   + f_0 * kk_429[k];

        t_286[k] = -2.0 * hk_178[k]
                   + f_0 * kk_430[k];

        t_287[k] = -2.0 * hk_179[k]
                   + f_0 * kk_431[k];

        t_288[k] = -hk_180[k]
                   + f_0 * kk_432[k];

        t_289[k] = -hk_181[k]
                   + f_0 * kk_433[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, hk_182, hk_183, hk_184, hk_185, \
                         hk_186, kk_434, kk_435, kk_436, kk_437, \
                         kk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -hk_182[k]
                   + f_0 * kk_434[k];

        t_291[k] = -hk_183[k]
                   + f_0 * kk_435[k];

        t_292[k] = -hk_184[k]
                   + f_0 * kk_436[k];

        t_293[k] = -hk_185[k]
                   + f_0 * kk_437[k];

        t_294[k] = -hk_186[k]
                   + f_0 * kk_438[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, hk_187, hk_188, hk_189, hk_190, \
                         hk_191, kk_439, kk_440, kk_441, kk_442, \
                         kk_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -hk_187[k]
                   + f_0 * kk_439[k];

        t_296[k] = -hk_188[k]
                   + f_0 * kk_440[k];

        t_297[k] = -hk_189[k]
                   + f_0 * kk_441[k];

        t_298[k] = -hk_190[k]
                   + f_0 * kk_442[k];

        t_299[k] = -hk_191[k]
                   + f_0 * kk_443[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, hk_192, hk_193, hk_194, hk_195, \
                         hk_196, kk_444, kk_445, kk_446, kk_447, \
                         kk_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -hk_192[k]
                   + f_0 * kk_444[k];

        t_301[k] = -hk_193[k]
                   + f_0 * kk_445[k];

        t_302[k] = -hk_194[k]
                   + f_0 * kk_446[k];

        t_303[k] = -hk_195[k]
                   + f_0 * kk_447[k];

        t_304[k] = -hk_196[k]
                   + f_0 * kk_448[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, hk_197, hk_198, hk_199, hk_200, \
                         hk_201, kk_449, kk_450, kk_451, kk_452, \
                         kk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -hk_197[k]
                   + f_0 * kk_449[k];

        t_306[k] = -hk_198[k]
                   + f_0 * kk_450[k];

        t_307[k] = -hk_199[k]
                   + f_0 * kk_451[k];

        t_308[k] = -hk_200[k]
                   + f_0 * kk_452[k];

        t_309[k] = -hk_201[k]
                   + f_0 * kk_453[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, hk_202, hk_203, hk_204, hk_205, \
                         hk_206, kk_454, kk_455, kk_456, kk_457, \
                         kk_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -hk_202[k]
                   + f_0 * kk_454[k];

        t_311[k] = -hk_203[k]
                   + f_0 * kk_455[k];

        t_312[k] = -hk_204[k]
                   + f_0 * kk_456[k];

        t_313[k] = -hk_205[k]
                   + f_0 * kk_457[k];

        t_314[k] = -hk_206[k]
                   + f_0 * kk_458[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, hk_207, hk_208, hk_209, hk_210, \
                         hk_211, kk_459, kk_460, kk_461, kk_462, \
                         kk_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -hk_207[k]
                   + f_0 * kk_459[k];

        t_316[k] = -hk_208[k]
                   + f_0 * kk_460[k];

        t_317[k] = -hk_209[k]
                   + f_0 * kk_461[k];

        t_318[k] = -hk_210[k]
                   + f_0 * kk_462[k];

        t_319[k] = -hk_211[k]
                   + f_0 * kk_463[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, t_325, hk_212, hk_213, hk_214, \
                         hk_215, kk_464, kk_465, kk_466, kk_467, kk_468, \
                         kk_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -hk_212[k]
                   + f_0 * kk_464[k];

        t_321[k] = -hk_213[k]
                   + f_0 * kk_465[k];

        t_322[k] = -hk_214[k]
                   + f_0 * kk_466[k];

        t_323[k] = -hk_215[k]
                   + f_0 * kk_467[k];

        t_324[k] = f_0 * kk_468[k];

        t_325[k] = f_0 * kk_469[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, t_332, t_333, kk_470, \
                         kk_471, kk_472, kk_473, kk_474, kk_475, kk_476, \
                         kk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_0 * kk_470[k];

        t_327[k] = f_0 * kk_471[k];

        t_328[k] = f_0 * kk_472[k];

        t_329[k] = f_0 * kk_473[k];

        t_330[k] = f_0 * kk_474[k];

        t_331[k] = f_0 * kk_475[k];

        t_332[k] = f_0 * kk_476[k];

        t_333[k] = f_0 * kk_477[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_216 = buffer.data(hk + 216);
    const auto *hk_217 = buffer.data(hk + 217);
    const auto *hk_218 = buffer.data(hk + 218);
    const auto *hk_219 = buffer.data(hk + 219);
    const auto *hk_220 = buffer.data(hk + 220);
    const auto *hk_221 = buffer.data(hk + 221);
    const auto *hk_222 = buffer.data(hk + 222);
    const auto *hk_223 = buffer.data(hk + 223);
    const auto *hk_224 = buffer.data(hk + 224);
    const auto *hk_225 = buffer.data(hk + 225);
    const auto *hk_226 = buffer.data(hk + 226);
    const auto *hk_227 = buffer.data(hk + 227);
    const auto *hk_228 = buffer.data(hk + 228);
    const auto *hk_229 = buffer.data(hk + 229);
    const auto *hk_230 = buffer.data(hk + 230);
    const auto *hk_231 = buffer.data(hk + 231);
    const auto *hk_232 = buffer.data(hk + 232);
    const auto *hk_233 = buffer.data(hk + 233);
    const auto *hk_234 = buffer.data(hk + 234);
    const auto *hk_235 = buffer.data(hk + 235);
    const auto *hk_236 = buffer.data(hk + 236);
    const auto *hk_237 = buffer.data(hk + 237);
    const auto *hk_238 = buffer.data(hk + 238);
    const auto *hk_239 = buffer.data(hk + 239);
    const auto *hk_240 = buffer.data(hk + 240);
    const auto *hk_241 = buffer.data(hk + 241);
    const auto *hk_242 = buffer.data(hk + 242);
    const auto *hk_243 = buffer.data(hk + 243);
    const auto *hk_244 = buffer.data(hk + 244);
    const auto *hk_245 = buffer.data(hk + 245);
    const auto *hk_246 = buffer.data(hk + 246);
    const auto *hk_247 = buffer.data(hk + 247);
    const auto *hk_248 = buffer.data(hk + 248);
    const auto *hk_249 = buffer.data(hk + 249);
    const auto *hk_250 = buffer.data(hk + 250);
    const auto *hk_251 = buffer.data(hk + 251);
    const auto *hk_252 = buffer.data(hk + 252);
    const auto *hk_253 = buffer.data(hk + 253);
    const auto *hk_254 = buffer.data(hk + 254);
    const auto *hk_255 = buffer.data(hk + 255);
    const auto *hk_256 = buffer.data(hk + 256);
    const auto *hk_257 = buffer.data(hk + 257);
    const auto *hk_258 = buffer.data(hk + 258);
    const auto *hk_259 = buffer.data(hk + 259);
    const auto *hk_260 = buffer.data(hk + 260);
    const auto *hk_261 = buffer.data(hk + 261);
    const auto *hk_262 = buffer.data(hk + 262);
    const auto *hk_263 = buffer.data(hk + 263);
    const auto *hk_264 = buffer.data(hk + 264);
    const auto *hk_265 = buffer.data(hk + 265);
    const auto *hk_266 = buffer.data(hk + 266);
    const auto *hk_267 = buffer.data(hk + 267);
    const auto *hk_268 = buffer.data(hk + 268);
    const auto *hk_269 = buffer.data(hk + 269);
    const auto *hk_270 = buffer.data(hk + 270);
    const auto *hk_271 = buffer.data(hk + 271);
    const auto *hk_272 = buffer.data(hk + 272);
    const auto *hk_273 = buffer.data(hk + 273);
    const auto *hk_274 = buffer.data(hk + 274);
    const auto *hk_275 = buffer.data(hk + 275);
    const auto *hk_276 = buffer.data(hk + 276);
    const auto *hk_277 = buffer.data(hk + 277);
    const auto *hk_278 = buffer.data(hk + 278);
    const auto *hk_279 = buffer.data(hk + 279);
    const auto *hk_280 = buffer.data(hk + 280);
    const auto *hk_281 = buffer.data(hk + 281);
    const auto *hk_282 = buffer.data(hk + 282);
    const auto *hk_283 = buffer.data(hk + 283);
    const auto *hk_284 = buffer.data(hk + 284);
    const auto *hk_285 = buffer.data(hk + 285);
    const auto *hk_286 = buffer.data(hk + 286);
    const auto *hk_287 = buffer.data(hk + 287);
    const auto *hk_288 = buffer.data(hk + 288);
    const auto *hk_289 = buffer.data(hk + 289);
    const auto *hk_290 = buffer.data(hk + 290);
    const auto *hk_291 = buffer.data(hk + 291);
    const auto *hk_292 = buffer.data(hk + 292);
    const auto *hk_293 = buffer.data(hk + 293);
    const auto *hk_294 = buffer.data(hk + 294);
    const auto *hk_295 = buffer.data(hk + 295);
    const auto *hk_296 = buffer.data(hk + 296);
    const auto *hk_297 = buffer.data(hk + 297);
    const auto *hk_298 = buffer.data(hk + 298);
    const auto *hk_299 = buffer.data(hk + 299);
    const auto *hk_300 = buffer.data(hk + 300);
    const auto *hk_301 = buffer.data(hk + 301);
    const auto *hk_302 = buffer.data(hk + 302);
    const auto *hk_303 = buffer.data(hk + 303);
    const auto *hk_304 = buffer.data(hk + 304);
    const auto *hk_305 = buffer.data(hk + 305);
    const auto *hk_306 = buffer.data(hk + 306);
    const auto *hk_307 = buffer.data(hk + 307);
    const auto *hk_308 = buffer.data(hk + 308);
    const auto *hk_309 = buffer.data(hk + 309);
    const auto *hk_310 = buffer.data(hk + 310);
    const auto *hk_311 = buffer.data(hk + 311);
    const auto *hk_312 = buffer.data(hk + 312);
    const auto *hk_313 = buffer.data(hk + 313);
    const auto *hk_314 = buffer.data(hk + 314);
    const auto *hk_315 = buffer.data(hk + 315);
    const auto *hk_316 = buffer.data(hk + 316);
    const auto *hk_317 = buffer.data(hk + 317);
    const auto *hk_318 = buffer.data(hk + 318);
    const auto *hk_319 = buffer.data(hk + 319);
    const auto *hk_320 = buffer.data(hk + 320);
    const auto *hk_321 = buffer.data(hk + 321);
    const auto *hk_322 = buffer.data(hk + 322);
    const auto *hk_323 = buffer.data(hk + 323);
    const auto *hk_324 = buffer.data(hk + 324);
    const auto *hk_325 = buffer.data(hk + 325);
    const auto *hk_326 = buffer.data(hk + 326);
    const auto *hk_327 = buffer.data(hk + 327);
    const auto *hk_328 = buffer.data(hk + 328);
    const auto *hk_329 = buffer.data(hk + 329);
    const auto *hk_330 = buffer.data(hk + 330);
    const auto *hk_331 = buffer.data(hk + 331);
    const auto *hk_332 = buffer.data(hk + 332);
    const auto *hk_333 = buffer.data(hk + 333);
    const auto *hk_334 = buffer.data(hk + 334);
    const auto *hk_335 = buffer.data(hk + 335);
    const auto *hk_336 = buffer.data(hk + 336);
    const auto *hk_337 = buffer.data(hk + 337);
    const auto *hk_338 = buffer.data(hk + 338);
    const auto *hk_339 = buffer.data(hk + 339);
    const auto *hk_340 = buffer.data(hk + 340);
    const auto *hk_341 = buffer.data(hk + 341);
    const auto *hk_342 = buffer.data(hk + 342);
    const auto *hk_343 = buffer.data(hk + 343);
    const auto *hk_344 = buffer.data(hk + 344);
    const auto *hk_345 = buffer.data(hk + 345);
    const auto *hk_346 = buffer.data(hk + 346);
    const auto *hk_347 = buffer.data(hk + 347);
    const auto *hk_348 = buffer.data(hk + 348);
    const auto *hk_349 = buffer.data(hk + 349);

    const auto *kk_478 = buffer.data(kk + 478);
    const auto *kk_479 = buffer.data(kk + 479);
    const auto *kk_480 = buffer.data(kk + 480);
    const auto *kk_481 = buffer.data(kk + 481);
    const auto *kk_482 = buffer.data(kk + 482);
    const auto *kk_483 = buffer.data(kk + 483);
    const auto *kk_484 = buffer.data(kk + 484);
    const auto *kk_485 = buffer.data(kk + 485);
    const auto *kk_486 = buffer.data(kk + 486);
    const auto *kk_487 = buffer.data(kk + 487);
    const auto *kk_488 = buffer.data(kk + 488);
    const auto *kk_489 = buffer.data(kk + 489);
    const auto *kk_490 = buffer.data(kk + 490);
    const auto *kk_491 = buffer.data(kk + 491);
    const auto *kk_492 = buffer.data(kk + 492);
    const auto *kk_493 = buffer.data(kk + 493);
    const auto *kk_494 = buffer.data(kk + 494);
    const auto *kk_495 = buffer.data(kk + 495);
    const auto *kk_496 = buffer.data(kk + 496);
    const auto *kk_497 = buffer.data(kk + 497);
    const auto *kk_498 = buffer.data(kk + 498);
    const auto *kk_499 = buffer.data(kk + 499);
    const auto *kk_500 = buffer.data(kk + 500);
    const auto *kk_501 = buffer.data(kk + 501);
    const auto *kk_502 = buffer.data(kk + 502);
    const auto *kk_503 = buffer.data(kk + 503);
    const auto *kk_540 = buffer.data(kk + 540);
    const auto *kk_541 = buffer.data(kk + 541);
    const auto *kk_542 = buffer.data(kk + 542);
    const auto *kk_543 = buffer.data(kk + 543);
    const auto *kk_544 = buffer.data(kk + 544);
    const auto *kk_545 = buffer.data(kk + 545);
    const auto *kk_546 = buffer.data(kk + 546);
    const auto *kk_547 = buffer.data(kk + 547);
    const auto *kk_548 = buffer.data(kk + 548);
    const auto *kk_549 = buffer.data(kk + 549);
    const auto *kk_550 = buffer.data(kk + 550);
    const auto *kk_551 = buffer.data(kk + 551);
    const auto *kk_552 = buffer.data(kk + 552);
    const auto *kk_553 = buffer.data(kk + 553);
    const auto *kk_554 = buffer.data(kk + 554);
    const auto *kk_555 = buffer.data(kk + 555);
    const auto *kk_556 = buffer.data(kk + 556);
    const auto *kk_557 = buffer.data(kk + 557);
    const auto *kk_558 = buffer.data(kk + 558);
    const auto *kk_559 = buffer.data(kk + 559);
    const auto *kk_560 = buffer.data(kk + 560);
    const auto *kk_561 = buffer.data(kk + 561);
    const auto *kk_562 = buffer.data(kk + 562);
    const auto *kk_563 = buffer.data(kk + 563);
    const auto *kk_564 = buffer.data(kk + 564);
    const auto *kk_565 = buffer.data(kk + 565);
    const auto *kk_566 = buffer.data(kk + 566);
    const auto *kk_567 = buffer.data(kk + 567);
    const auto *kk_568 = buffer.data(kk + 568);
    const auto *kk_569 = buffer.data(kk + 569);
    const auto *kk_570 = buffer.data(kk + 570);
    const auto *kk_571 = buffer.data(kk + 571);
    const auto *kk_572 = buffer.data(kk + 572);
    const auto *kk_573 = buffer.data(kk + 573);
    const auto *kk_574 = buffer.data(kk + 574);
    const auto *kk_575 = buffer.data(kk + 575);
    const auto *kk_576 = buffer.data(kk + 576);
    const auto *kk_577 = buffer.data(kk + 577);
    const auto *kk_578 = buffer.data(kk + 578);
    const auto *kk_579 = buffer.data(kk + 579);
    const auto *kk_580 = buffer.data(kk + 580);
    const auto *kk_581 = buffer.data(kk + 581);
    const auto *kk_582 = buffer.data(kk + 582);
    const auto *kk_583 = buffer.data(kk + 583);
    const auto *kk_584 = buffer.data(kk + 584);
    const auto *kk_585 = buffer.data(kk + 585);
    const auto *kk_586 = buffer.data(kk + 586);
    const auto *kk_587 = buffer.data(kk + 587);
    const auto *kk_588 = buffer.data(kk + 588);
    const auto *kk_589 = buffer.data(kk + 589);
    const auto *kk_590 = buffer.data(kk + 590);
    const auto *kk_591 = buffer.data(kk + 591);
    const auto *kk_592 = buffer.data(kk + 592);
    const auto *kk_593 = buffer.data(kk + 593);
    const auto *kk_594 = buffer.data(kk + 594);
    const auto *kk_595 = buffer.data(kk + 595);
    const auto *kk_596 = buffer.data(kk + 596);
    const auto *kk_597 = buffer.data(kk + 597);
    const auto *kk_598 = buffer.data(kk + 598);
    const auto *kk_599 = buffer.data(kk + 599);
    const auto *kk_600 = buffer.data(kk + 600);
    const auto *kk_601 = buffer.data(kk + 601);
    const auto *kk_602 = buffer.data(kk + 602);
    const auto *kk_603 = buffer.data(kk + 603);
    const auto *kk_604 = buffer.data(kk + 604);
    const auto *kk_605 = buffer.data(kk + 605);
    const auto *kk_606 = buffer.data(kk + 606);
    const auto *kk_607 = buffer.data(kk + 607);
    const auto *kk_608 = buffer.data(kk + 608);
    const auto *kk_609 = buffer.data(kk + 609);
    const auto *kk_610 = buffer.data(kk + 610);
    const auto *kk_611 = buffer.data(kk + 611);
    const auto *kk_612 = buffer.data(kk + 612);
    const auto *kk_613 = buffer.data(kk + 613);
    const auto *kk_614 = buffer.data(kk + 614);
    const auto *kk_615 = buffer.data(kk + 615);
    const auto *kk_616 = buffer.data(kk + 616);
    const auto *kk_617 = buffer.data(kk + 617);
    const auto *kk_618 = buffer.data(kk + 618);
    const auto *kk_619 = buffer.data(kk + 619);
    const auto *kk_620 = buffer.data(kk + 620);
    const auto *kk_621 = buffer.data(kk + 621);
    const auto *kk_622 = buffer.data(kk + 622);
    const auto *kk_623 = buffer.data(kk + 623);
    const auto *kk_624 = buffer.data(kk + 624);
    const auto *kk_625 = buffer.data(kk + 625);
    const auto *kk_626 = buffer.data(kk + 626);
    const auto *kk_627 = buffer.data(kk + 627);
    const auto *kk_628 = buffer.data(kk + 628);
    const auto *kk_629 = buffer.data(kk + 629);
    const auto *kk_630 = buffer.data(kk + 630);
    const auto *kk_631 = buffer.data(kk + 631);
    const auto *kk_632 = buffer.data(kk + 632);
    const auto *kk_633 = buffer.data(kk + 633);
    const auto *kk_634 = buffer.data(kk + 634);
    const auto *kk_635 = buffer.data(kk + 635);
    const auto *kk_636 = buffer.data(kk + 636);
    const auto *kk_637 = buffer.data(kk + 637);
    const auto *kk_638 = buffer.data(kk + 638);
    const auto *kk_639 = buffer.data(kk + 639);
    const auto *kk_640 = buffer.data(kk + 640);
    const auto *kk_641 = buffer.data(kk + 641);
    const auto *kk_642 = buffer.data(kk + 642);
    const auto *kk_643 = buffer.data(kk + 643);
    const auto *kk_644 = buffer.data(kk + 644);
    const auto *kk_645 = buffer.data(kk + 645);
    const auto *kk_646 = buffer.data(kk + 646);
    const auto *kk_647 = buffer.data(kk + 647);
    const auto *kk_648 = buffer.data(kk + 648);
    const auto *kk_649 = buffer.data(kk + 649);
    const auto *kk_650 = buffer.data(kk + 650);
    const auto *kk_651 = buffer.data(kk + 651);
    const auto *kk_652 = buffer.data(kk + 652);
    const auto *kk_653 = buffer.data(kk + 653);
    const auto *kk_654 = buffer.data(kk + 654);
    const auto *kk_655 = buffer.data(kk + 655);
    const auto *kk_656 = buffer.data(kk + 656);
    const auto *kk_657 = buffer.data(kk + 657);
    const auto *kk_658 = buffer.data(kk + 658);
    const auto *kk_659 = buffer.data(kk + 659);
    const auto *kk_660 = buffer.data(kk + 660);
    const auto *kk_661 = buffer.data(kk + 661);
    const auto *kk_662 = buffer.data(kk + 662);
    const auto *kk_663 = buffer.data(kk + 663);
    const auto *kk_664 = buffer.data(kk + 664);
    const auto *kk_665 = buffer.data(kk + 665);
    const auto *kk_666 = buffer.data(kk + 666);
    const auto *kk_667 = buffer.data(kk + 667);
    const auto *kk_668 = buffer.data(kk + 668);
    const auto *kk_669 = buffer.data(kk + 669);
    const auto *kk_670 = buffer.data(kk + 670);
    const auto *kk_671 = buffer.data(kk + 671);
    const auto *kk_672 = buffer.data(kk + 672);
    const auto *kk_673 = buffer.data(kk + 673);

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, t_339, t_340, t_341, kk_478, \
                         kk_479, kk_480, kk_481, kk_482, kk_483, kk_484, \
                         kk_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_0 * kk_478[k];

        t_335[k] = f_0 * kk_479[k];

        t_336[k] = f_0 * kk_480[k];

        t_337[k] = f_0 * kk_481[k];

        t_338[k] = f_0 * kk_482[k];

        t_339[k] = f_0 * kk_483[k];

        t_340[k] = f_0 * kk_484[k];

        t_341[k] = f_0 * kk_485[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, t_347, t_348, t_349, kk_486, \
                         kk_487, kk_488, kk_489, kk_490, kk_491, kk_492, \
                         kk_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_0 * kk_486[k];

        t_343[k] = f_0 * kk_487[k];

        t_344[k] = f_0 * kk_488[k];

        t_345[k] = f_0 * kk_489[k];

        t_346[k] = f_0 * kk_490[k];

        t_347[k] = f_0 * kk_491[k];

        t_348[k] = f_0 * kk_492[k];

        t_349[k] = f_0 * kk_493[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, t_355, t_356, t_357, kk_494, \
                         kk_495, kk_496, kk_497, kk_498, kk_499, kk_500, \
                         kk_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_0 * kk_494[k];

        t_351[k] = f_0 * kk_495[k];

        t_352[k] = f_0 * kk_496[k];

        t_353[k] = f_0 * kk_497[k];

        t_354[k] = f_0 * kk_498[k];

        t_355[k] = f_0 * kk_499[k];

        t_356[k] = f_0 * kk_500[k];

        t_357[k] = f_0 * kk_501[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, t_363, hk_216, hk_217, hk_218, \
                         hk_219, kk_502, kk_503, kk_540, kk_541, kk_542, \
                         kk_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_0 * kk_502[k];

        t_359[k] = f_0 * kk_503[k];

        t_360[k] = -4.0 * hk_216[k]
                   + f_0 * kk_540[k];

        t_361[k] = -4.0 * hk_217[k]
                   + f_0 * kk_541[k];

        t_362[k] = -4.0 * hk_218[k]
                   + f_0 * kk_542[k];

        t_363[k] = -4.0 * hk_219[k]
                   + f_0 * kk_543[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, t_368, hk_220, hk_221, hk_222, hk_223, \
                         hk_224, kk_544, kk_545, kk_546, kk_547, \
                         kk_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = -4.0 * hk_220[k]
                   + f_0 * kk_544[k];

        t_365[k] = -4.0 * hk_221[k]
                   + f_0 * kk_545[k];

        t_366[k] = -4.0 * hk_222[k]
                   + f_0 * kk_546[k];

        t_367[k] = -4.0 * hk_223[k]
                   + f_0 * kk_547[k];

        t_368[k] = -4.0 * hk_224[k]
                   + f_0 * kk_548[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, t_373, hk_225, hk_226, hk_227, hk_228, \
                         hk_229, kk_549, kk_550, kk_551, kk_552, \
                         kk_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = -4.0 * hk_225[k]
                   + f_0 * kk_549[k];

        t_370[k] = -4.0 * hk_226[k]
                   + f_0 * kk_550[k];

        t_371[k] = -4.0 * hk_227[k]
                   + f_0 * kk_551[k];

        t_372[k] = -4.0 * hk_228[k]
                   + f_0 * kk_552[k];

        t_373[k] = -4.0 * hk_229[k]
                   + f_0 * kk_553[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, t_378, hk_230, hk_231, hk_232, hk_233, \
                         hk_234, kk_554, kk_555, kk_556, kk_557, \
                         kk_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = -4.0 * hk_230[k]
                   + f_0 * kk_554[k];

        t_375[k] = -4.0 * hk_231[k]
                   + f_0 * kk_555[k];

        t_376[k] = -4.0 * hk_232[k]
                   + f_0 * kk_556[k];

        t_377[k] = -4.0 * hk_233[k]
                   + f_0 * kk_557[k];

        t_378[k] = -4.0 * hk_234[k]
                   + f_0 * kk_558[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, hk_235, hk_236, hk_237, hk_238, \
                         hk_239, kk_559, kk_560, kk_561, kk_562, \
                         kk_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = -4.0 * hk_235[k]
                   + f_0 * kk_559[k];

        t_380[k] = -4.0 * hk_236[k]
                   + f_0 * kk_560[k];

        t_381[k] = -4.0 * hk_237[k]
                   + f_0 * kk_561[k];

        t_382[k] = -4.0 * hk_238[k]
                   + f_0 * kk_562[k];

        t_383[k] = -4.0 * hk_239[k]
                   + f_0 * kk_563[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, t_388, hk_240, hk_241, hk_242, hk_243, \
                         hk_244, kk_564, kk_565, kk_566, kk_567, \
                         kk_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = -4.0 * hk_240[k]
                   + f_0 * kk_564[k];

        t_385[k] = -4.0 * hk_241[k]
                   + f_0 * kk_565[k];

        t_386[k] = -4.0 * hk_242[k]
                   + f_0 * kk_566[k];

        t_387[k] = -4.0 * hk_243[k]
                   + f_0 * kk_567[k];

        t_388[k] = -4.0 * hk_244[k]
                   + f_0 * kk_568[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, hk_245, hk_246, hk_247, hk_248, \
                         hk_249, kk_569, kk_570, kk_571, kk_572, \
                         kk_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = -4.0 * hk_245[k]
                   + f_0 * kk_569[k];

        t_390[k] = -4.0 * hk_246[k]
                   + f_0 * kk_570[k];

        t_391[k] = -4.0 * hk_247[k]
                   + f_0 * kk_571[k];

        t_392[k] = -4.0 * hk_248[k]
                   + f_0 * kk_572[k];

        t_393[k] = -4.0 * hk_249[k]
                   + f_0 * kk_573[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, hk_250, hk_251, hk_252, hk_253, \
                         hk_254, kk_574, kk_575, kk_576, kk_577, \
                         kk_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = -4.0 * hk_250[k]
                   + f_0 * kk_574[k];

        t_395[k] = -4.0 * hk_251[k]
                   + f_0 * kk_575[k];

        t_396[k] = -3.0 * hk_252[k]
                   + f_0 * kk_576[k];

        t_397[k] = -3.0 * hk_253[k]
                   + f_0 * kk_577[k];

        t_398[k] = -3.0 * hk_254[k]
                   + f_0 * kk_578[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, t_403, hk_255, hk_256, hk_257, hk_258, \
                         hk_259, kk_579, kk_580, kk_581, kk_582, \
                         kk_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = -3.0 * hk_255[k]
                   + f_0 * kk_579[k];

        t_400[k] = -3.0 * hk_256[k]
                   + f_0 * kk_580[k];

        t_401[k] = -3.0 * hk_257[k]
                   + f_0 * kk_581[k];

        t_402[k] = -3.0 * hk_258[k]
                   + f_0 * kk_582[k];

        t_403[k] = -3.0 * hk_259[k]
                   + f_0 * kk_583[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, t_408, hk_260, hk_261, hk_262, hk_263, \
                         hk_264, kk_584, kk_585, kk_586, kk_587, \
                         kk_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = -3.0 * hk_260[k]
                   + f_0 * kk_584[k];

        t_405[k] = -3.0 * hk_261[k]
                   + f_0 * kk_585[k];

        t_406[k] = -3.0 * hk_262[k]
                   + f_0 * kk_586[k];

        t_407[k] = -3.0 * hk_263[k]
                   + f_0 * kk_587[k];

        t_408[k] = -3.0 * hk_264[k]
                   + f_0 * kk_588[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, hk_265, hk_266, hk_267, hk_268, \
                         hk_269, kk_589, kk_590, kk_591, kk_592, \
                         kk_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = -3.0 * hk_265[k]
                   + f_0 * kk_589[k];

        t_410[k] = -3.0 * hk_266[k]
                   + f_0 * kk_590[k];

        t_411[k] = -3.0 * hk_267[k]
                   + f_0 * kk_591[k];

        t_412[k] = -3.0 * hk_268[k]
                   + f_0 * kk_592[k];

        t_413[k] = -3.0 * hk_269[k]
                   + f_0 * kk_593[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, hk_270, hk_271, hk_272, hk_273, \
                         hk_274, kk_594, kk_595, kk_596, kk_597, \
                         kk_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = -3.0 * hk_270[k]
                   + f_0 * kk_594[k];

        t_415[k] = -3.0 * hk_271[k]
                   + f_0 * kk_595[k];

        t_416[k] = -3.0 * hk_272[k]
                   + f_0 * kk_596[k];

        t_417[k] = -3.0 * hk_273[k]
                   + f_0 * kk_597[k];

        t_418[k] = -3.0 * hk_274[k]
                   + f_0 * kk_598[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, t_423, hk_275, hk_276, hk_277, hk_278, \
                         hk_279, kk_599, kk_600, kk_601, kk_602, \
                         kk_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = -3.0 * hk_275[k]
                   + f_0 * kk_599[k];

        t_420[k] = -3.0 * hk_276[k]
                   + f_0 * kk_600[k];

        t_421[k] = -3.0 * hk_277[k]
                   + f_0 * kk_601[k];

        t_422[k] = -3.0 * hk_278[k]
                   + f_0 * kk_602[k];

        t_423[k] = -3.0 * hk_279[k]
                   + f_0 * kk_603[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, t_428, hk_280, hk_281, hk_282, hk_283, \
                         hk_284, kk_604, kk_605, kk_606, kk_607, \
                         kk_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = -3.0 * hk_280[k]
                   + f_0 * kk_604[k];

        t_425[k] = -3.0 * hk_281[k]
                   + f_0 * kk_605[k];

        t_426[k] = -3.0 * hk_282[k]
                   + f_0 * kk_606[k];

        t_427[k] = -3.0 * hk_283[k]
                   + f_0 * kk_607[k];

        t_428[k] = -3.0 * hk_284[k]
                   + f_0 * kk_608[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, hk_285, hk_286, hk_287, hk_288, \
                         hk_289, kk_609, kk_610, kk_611, kk_612, \
                         kk_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = -3.0 * hk_285[k]
                   + f_0 * kk_609[k];

        t_430[k] = -3.0 * hk_286[k]
                   + f_0 * kk_610[k];

        t_431[k] = -3.0 * hk_287[k]
                   + f_0 * kk_611[k];

        t_432[k] = -2.0 * hk_288[k]
                   + f_0 * kk_612[k];

        t_433[k] = -2.0 * hk_289[k]
                   + f_0 * kk_613[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, t_438, hk_290, hk_291, hk_292, hk_293, \
                         hk_294, kk_614, kk_615, kk_616, kk_617, \
                         kk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = -2.0 * hk_290[k]
                   + f_0 * kk_614[k];

        t_435[k] = -2.0 * hk_291[k]
                   + f_0 * kk_615[k];

        t_436[k] = -2.0 * hk_292[k]
                   + f_0 * kk_616[k];

        t_437[k] = -2.0 * hk_293[k]
                   + f_0 * kk_617[k];

        t_438[k] = -2.0 * hk_294[k]
                   + f_0 * kk_618[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, t_443, hk_295, hk_296, hk_297, hk_298, \
                         hk_299, kk_619, kk_620, kk_621, kk_622, \
                         kk_623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = -2.0 * hk_295[k]
                   + f_0 * kk_619[k];

        t_440[k] = -2.0 * hk_296[k]
                   + f_0 * kk_620[k];

        t_441[k] = -2.0 * hk_297[k]
                   + f_0 * kk_621[k];

        t_442[k] = -2.0 * hk_298[k]
                   + f_0 * kk_622[k];

        t_443[k] = -2.0 * hk_299[k]
                   + f_0 * kk_623[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, hk_300, hk_301, hk_302, hk_303, \
                         hk_304, kk_624, kk_625, kk_626, kk_627, \
                         kk_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = -2.0 * hk_300[k]
                   + f_0 * kk_624[k];

        t_445[k] = -2.0 * hk_301[k]
                   + f_0 * kk_625[k];

        t_446[k] = -2.0 * hk_302[k]
                   + f_0 * kk_626[k];

        t_447[k] = -2.0 * hk_303[k]
                   + f_0 * kk_627[k];

        t_448[k] = -2.0 * hk_304[k]
                   + f_0 * kk_628[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, hk_305, hk_306, hk_307, hk_308, \
                         hk_309, kk_629, kk_630, kk_631, kk_632, \
                         kk_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = -2.0 * hk_305[k]
                   + f_0 * kk_629[k];

        t_450[k] = -2.0 * hk_306[k]
                   + f_0 * kk_630[k];

        t_451[k] = -2.0 * hk_307[k]
                   + f_0 * kk_631[k];

        t_452[k] = -2.0 * hk_308[k]
                   + f_0 * kk_632[k];

        t_453[k] = -2.0 * hk_309[k]
                   + f_0 * kk_633[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, hk_310, hk_311, hk_312, hk_313, \
                         hk_314, kk_634, kk_635, kk_636, kk_637, \
                         kk_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = -2.0 * hk_310[k]
                   + f_0 * kk_634[k];

        t_455[k] = -2.0 * hk_311[k]
                   + f_0 * kk_635[k];

        t_456[k] = -2.0 * hk_312[k]
                   + f_0 * kk_636[k];

        t_457[k] = -2.0 * hk_313[k]
                   + f_0 * kk_637[k];

        t_458[k] = -2.0 * hk_314[k]
                   + f_0 * kk_638[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, t_463, hk_315, hk_316, hk_317, hk_318, \
                         hk_319, kk_639, kk_640, kk_641, kk_642, \
                         kk_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = -2.0 * hk_315[k]
                   + f_0 * kk_639[k];

        t_460[k] = -2.0 * hk_316[k]
                   + f_0 * kk_640[k];

        t_461[k] = -2.0 * hk_317[k]
                   + f_0 * kk_641[k];

        t_462[k] = -2.0 * hk_318[k]
                   + f_0 * kk_642[k];

        t_463[k] = -2.0 * hk_319[k]
                   + f_0 * kk_643[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, t_468, hk_320, hk_321, hk_322, hk_323, \
                         hk_324, kk_644, kk_645, kk_646, kk_647, \
                         kk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = -2.0 * hk_320[k]
                   + f_0 * kk_644[k];

        t_465[k] = -2.0 * hk_321[k]
                   + f_0 * kk_645[k];

        t_466[k] = -2.0 * hk_322[k]
                   + f_0 * kk_646[k];

        t_467[k] = -2.0 * hk_323[k]
                   + f_0 * kk_647[k];

        t_468[k] = -hk_324[k]
                   + f_0 * kk_648[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, t_473, hk_325, hk_326, hk_327, hk_328, \
                         hk_329, kk_649, kk_650, kk_651, kk_652, \
                         kk_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = -hk_325[k]
                   + f_0 * kk_649[k];

        t_470[k] = -hk_326[k]
                   + f_0 * kk_650[k];

        t_471[k] = -hk_327[k]
                   + f_0 * kk_651[k];

        t_472[k] = -hk_328[k]
                   + f_0 * kk_652[k];

        t_473[k] = -hk_329[k]
                   + f_0 * kk_653[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, t_478, hk_330, hk_331, hk_332, hk_333, \
                         hk_334, kk_654, kk_655, kk_656, kk_657, \
                         kk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = -hk_330[k]
                   + f_0 * kk_654[k];

        t_475[k] = -hk_331[k]
                   + f_0 * kk_655[k];

        t_476[k] = -hk_332[k]
                   + f_0 * kk_656[k];

        t_477[k] = -hk_333[k]
                   + f_0 * kk_657[k];

        t_478[k] = -hk_334[k]
                   + f_0 * kk_658[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, t_483, hk_335, hk_336, hk_337, hk_338, \
                         hk_339, kk_659, kk_660, kk_661, kk_662, \
                         kk_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = -hk_335[k]
                   + f_0 * kk_659[k];

        t_480[k] = -hk_336[k]
                   + f_0 * kk_660[k];

        t_481[k] = -hk_337[k]
                   + f_0 * kk_661[k];

        t_482[k] = -hk_338[k]
                   + f_0 * kk_662[k];

        t_483[k] = -hk_339[k]
                   + f_0 * kk_663[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, t_488, hk_340, hk_341, hk_342, hk_343, \
                         hk_344, kk_664, kk_665, kk_666, kk_667, \
                         kk_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = -hk_340[k]
                   + f_0 * kk_664[k];

        t_485[k] = -hk_341[k]
                   + f_0 * kk_665[k];

        t_486[k] = -hk_342[k]
                   + f_0 * kk_666[k];

        t_487[k] = -hk_343[k]
                   + f_0 * kk_667[k];

        t_488[k] = -hk_344[k]
                   + f_0 * kk_668[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, hk_345, hk_346, hk_347, hk_348, \
                         hk_349, kk_669, kk_670, kk_671, kk_672, \
                         kk_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = -hk_345[k]
                   + f_0 * kk_669[k];

        t_490[k] = -hk_346[k]
                   + f_0 * kk_670[k];

        t_491[k] = -hk_347[k]
                   + f_0 * kk_671[k];

        t_492[k] = -hk_348[k]
                   + f_0 * kk_672[k];

        t_493[k] = -hk_349[k]
                   + f_0 * kk_673[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_350 = buffer.data(hk + 350);
    const auto *hk_351 = buffer.data(hk + 351);
    const auto *hk_352 = buffer.data(hk + 352);
    const auto *hk_353 = buffer.data(hk + 353);
    const auto *hk_354 = buffer.data(hk + 354);
    const auto *hk_355 = buffer.data(hk + 355);
    const auto *hk_356 = buffer.data(hk + 356);
    const auto *hk_357 = buffer.data(hk + 357);
    const auto *hk_358 = buffer.data(hk + 358);
    const auto *hk_359 = buffer.data(hk + 359);
    const auto *hk_360 = buffer.data(hk + 360);
    const auto *hk_361 = buffer.data(hk + 361);
    const auto *hk_362 = buffer.data(hk + 362);
    const auto *hk_363 = buffer.data(hk + 363);
    const auto *hk_364 = buffer.data(hk + 364);
    const auto *hk_365 = buffer.data(hk + 365);
    const auto *hk_366 = buffer.data(hk + 366);
    const auto *hk_367 = buffer.data(hk + 367);
    const auto *hk_368 = buffer.data(hk + 368);
    const auto *hk_369 = buffer.data(hk + 369);
    const auto *hk_370 = buffer.data(hk + 370);
    const auto *hk_371 = buffer.data(hk + 371);
    const auto *hk_372 = buffer.data(hk + 372);
    const auto *hk_373 = buffer.data(hk + 373);
    const auto *hk_374 = buffer.data(hk + 374);
    const auto *hk_375 = buffer.data(hk + 375);
    const auto *hk_376 = buffer.data(hk + 376);
    const auto *hk_377 = buffer.data(hk + 377);
    const auto *hk_378 = buffer.data(hk + 378);
    const auto *hk_379 = buffer.data(hk + 379);
    const auto *hk_380 = buffer.data(hk + 380);
    const auto *hk_381 = buffer.data(hk + 381);
    const auto *hk_382 = buffer.data(hk + 382);
    const auto *hk_383 = buffer.data(hk + 383);
    const auto *hk_384 = buffer.data(hk + 384);
    const auto *hk_385 = buffer.data(hk + 385);
    const auto *hk_386 = buffer.data(hk + 386);
    const auto *hk_387 = buffer.data(hk + 387);
    const auto *hk_388 = buffer.data(hk + 388);
    const auto *hk_389 = buffer.data(hk + 389);
    const auto *hk_390 = buffer.data(hk + 390);
    const auto *hk_391 = buffer.data(hk + 391);
    const auto *hk_392 = buffer.data(hk + 392);
    const auto *hk_393 = buffer.data(hk + 393);
    const auto *hk_394 = buffer.data(hk + 394);
    const auto *hk_395 = buffer.data(hk + 395);
    const auto *hk_396 = buffer.data(hk + 396);
    const auto *hk_397 = buffer.data(hk + 397);
    const auto *hk_398 = buffer.data(hk + 398);
    const auto *hk_399 = buffer.data(hk + 399);
    const auto *hk_400 = buffer.data(hk + 400);
    const auto *hk_401 = buffer.data(hk + 401);
    const auto *hk_402 = buffer.data(hk + 402);
    const auto *hk_403 = buffer.data(hk + 403);
    const auto *hk_404 = buffer.data(hk + 404);
    const auto *hk_405 = buffer.data(hk + 405);
    const auto *hk_406 = buffer.data(hk + 406);
    const auto *hk_407 = buffer.data(hk + 407);
    const auto *hk_408 = buffer.data(hk + 408);
    const auto *hk_409 = buffer.data(hk + 409);
    const auto *hk_410 = buffer.data(hk + 410);
    const auto *hk_411 = buffer.data(hk + 411);
    const auto *hk_412 = buffer.data(hk + 412);
    const auto *hk_413 = buffer.data(hk + 413);
    const auto *hk_414 = buffer.data(hk + 414);
    const auto *hk_415 = buffer.data(hk + 415);
    const auto *hk_416 = buffer.data(hk + 416);
    const auto *hk_417 = buffer.data(hk + 417);
    const auto *hk_418 = buffer.data(hk + 418);
    const auto *hk_419 = buffer.data(hk + 419);
    const auto *hk_420 = buffer.data(hk + 420);
    const auto *hk_421 = buffer.data(hk + 421);
    const auto *hk_422 = buffer.data(hk + 422);
    const auto *hk_423 = buffer.data(hk + 423);
    const auto *hk_424 = buffer.data(hk + 424);
    const auto *hk_425 = buffer.data(hk + 425);
    const auto *hk_426 = buffer.data(hk + 426);
    const auto *hk_427 = buffer.data(hk + 427);
    const auto *hk_428 = buffer.data(hk + 428);
    const auto *hk_429 = buffer.data(hk + 429);
    const auto *hk_430 = buffer.data(hk + 430);
    const auto *hk_431 = buffer.data(hk + 431);
    const auto *hk_432 = buffer.data(hk + 432);
    const auto *hk_433 = buffer.data(hk + 433);
    const auto *hk_434 = buffer.data(hk + 434);
    const auto *hk_435 = buffer.data(hk + 435);
    const auto *hk_436 = buffer.data(hk + 436);
    const auto *hk_437 = buffer.data(hk + 437);
    const auto *hk_438 = buffer.data(hk + 438);
    const auto *hk_439 = buffer.data(hk + 439);
    const auto *hk_440 = buffer.data(hk + 440);
    const auto *hk_441 = buffer.data(hk + 441);
    const auto *hk_442 = buffer.data(hk + 442);
    const auto *hk_443 = buffer.data(hk + 443);
    const auto *hk_444 = buffer.data(hk + 444);
    const auto *hk_445 = buffer.data(hk + 445);
    const auto *hk_446 = buffer.data(hk + 446);
    const auto *hk_447 = buffer.data(hk + 447);
    const auto *hk_448 = buffer.data(hk + 448);
    const auto *hk_449 = buffer.data(hk + 449);
    const auto *hk_450 = buffer.data(hk + 450);
    const auto *hk_451 = buffer.data(hk + 451);
    const auto *hk_452 = buffer.data(hk + 452);
    const auto *hk_453 = buffer.data(hk + 453);
    const auto *hk_454 = buffer.data(hk + 454);
    const auto *hk_455 = buffer.data(hk + 455);
    const auto *hk_456 = buffer.data(hk + 456);
    const auto *hk_457 = buffer.data(hk + 457);
    const auto *hk_458 = buffer.data(hk + 458);
    const auto *hk_459 = buffer.data(hk + 459);
    const auto *hk_460 = buffer.data(hk + 460);
    const auto *hk_461 = buffer.data(hk + 461);
    const auto *hk_462 = buffer.data(hk + 462);
    const auto *hk_463 = buffer.data(hk + 463);
    const auto *hk_464 = buffer.data(hk + 464);
    const auto *hk_465 = buffer.data(hk + 465);
    const auto *hk_466 = buffer.data(hk + 466);
    const auto *hk_467 = buffer.data(hk + 467);
    const auto *hk_468 = buffer.data(hk + 468);
    const auto *hk_469 = buffer.data(hk + 469);
    const auto *hk_470 = buffer.data(hk + 470);
    const auto *hk_471 = buffer.data(hk + 471);
    const auto *hk_472 = buffer.data(hk + 472);
    const auto *hk_473 = buffer.data(hk + 473);
    const auto *hk_474 = buffer.data(hk + 474);
    const auto *hk_475 = buffer.data(hk + 475);
    const auto *hk_476 = buffer.data(hk + 476);

    const auto *kk_674 = buffer.data(kk + 674);
    const auto *kk_675 = buffer.data(kk + 675);
    const auto *kk_676 = buffer.data(kk + 676);
    const auto *kk_677 = buffer.data(kk + 677);
    const auto *kk_678 = buffer.data(kk + 678);
    const auto *kk_679 = buffer.data(kk + 679);
    const auto *kk_680 = buffer.data(kk + 680);
    const auto *kk_681 = buffer.data(kk + 681);
    const auto *kk_682 = buffer.data(kk + 682);
    const auto *kk_683 = buffer.data(kk + 683);
    const auto *kk_684 = buffer.data(kk + 684);
    const auto *kk_685 = buffer.data(kk + 685);
    const auto *kk_686 = buffer.data(kk + 686);
    const auto *kk_687 = buffer.data(kk + 687);
    const auto *kk_688 = buffer.data(kk + 688);
    const auto *kk_689 = buffer.data(kk + 689);
    const auto *kk_690 = buffer.data(kk + 690);
    const auto *kk_691 = buffer.data(kk + 691);
    const auto *kk_692 = buffer.data(kk + 692);
    const auto *kk_693 = buffer.data(kk + 693);
    const auto *kk_694 = buffer.data(kk + 694);
    const auto *kk_695 = buffer.data(kk + 695);
    const auto *kk_696 = buffer.data(kk + 696);
    const auto *kk_697 = buffer.data(kk + 697);
    const auto *kk_698 = buffer.data(kk + 698);
    const auto *kk_699 = buffer.data(kk + 699);
    const auto *kk_700 = buffer.data(kk + 700);
    const auto *kk_701 = buffer.data(kk + 701);
    const auto *kk_702 = buffer.data(kk + 702);
    const auto *kk_703 = buffer.data(kk + 703);
    const auto *kk_704 = buffer.data(kk + 704);
    const auto *kk_705 = buffer.data(kk + 705);
    const auto *kk_706 = buffer.data(kk + 706);
    const auto *kk_707 = buffer.data(kk + 707);
    const auto *kk_708 = buffer.data(kk + 708);
    const auto *kk_709 = buffer.data(kk + 709);
    const auto *kk_710 = buffer.data(kk + 710);
    const auto *kk_711 = buffer.data(kk + 711);
    const auto *kk_712 = buffer.data(kk + 712);
    const auto *kk_713 = buffer.data(kk + 713);
    const auto *kk_714 = buffer.data(kk + 714);
    const auto *kk_715 = buffer.data(kk + 715);
    const auto *kk_716 = buffer.data(kk + 716);
    const auto *kk_717 = buffer.data(kk + 717);
    const auto *kk_718 = buffer.data(kk + 718);
    const auto *kk_719 = buffer.data(kk + 719);
    const auto *kk_756 = buffer.data(kk + 756);
    const auto *kk_757 = buffer.data(kk + 757);
    const auto *kk_758 = buffer.data(kk + 758);
    const auto *kk_759 = buffer.data(kk + 759);
    const auto *kk_760 = buffer.data(kk + 760);
    const auto *kk_761 = buffer.data(kk + 761);
    const auto *kk_762 = buffer.data(kk + 762);
    const auto *kk_763 = buffer.data(kk + 763);
    const auto *kk_764 = buffer.data(kk + 764);
    const auto *kk_765 = buffer.data(kk + 765);
    const auto *kk_766 = buffer.data(kk + 766);
    const auto *kk_767 = buffer.data(kk + 767);
    const auto *kk_768 = buffer.data(kk + 768);
    const auto *kk_769 = buffer.data(kk + 769);
    const auto *kk_770 = buffer.data(kk + 770);
    const auto *kk_771 = buffer.data(kk + 771);
    const auto *kk_772 = buffer.data(kk + 772);
    const auto *kk_773 = buffer.data(kk + 773);
    const auto *kk_774 = buffer.data(kk + 774);
    const auto *kk_775 = buffer.data(kk + 775);
    const auto *kk_776 = buffer.data(kk + 776);
    const auto *kk_777 = buffer.data(kk + 777);
    const auto *kk_778 = buffer.data(kk + 778);
    const auto *kk_779 = buffer.data(kk + 779);
    const auto *kk_780 = buffer.data(kk + 780);
    const auto *kk_781 = buffer.data(kk + 781);
    const auto *kk_782 = buffer.data(kk + 782);
    const auto *kk_783 = buffer.data(kk + 783);
    const auto *kk_784 = buffer.data(kk + 784);
    const auto *kk_785 = buffer.data(kk + 785);
    const auto *kk_786 = buffer.data(kk + 786);
    const auto *kk_787 = buffer.data(kk + 787);
    const auto *kk_788 = buffer.data(kk + 788);
    const auto *kk_789 = buffer.data(kk + 789);
    const auto *kk_790 = buffer.data(kk + 790);
    const auto *kk_791 = buffer.data(kk + 791);
    const auto *kk_792 = buffer.data(kk + 792);
    const auto *kk_793 = buffer.data(kk + 793);
    const auto *kk_794 = buffer.data(kk + 794);
    const auto *kk_795 = buffer.data(kk + 795);
    const auto *kk_796 = buffer.data(kk + 796);
    const auto *kk_797 = buffer.data(kk + 797);
    const auto *kk_798 = buffer.data(kk + 798);
    const auto *kk_799 = buffer.data(kk + 799);
    const auto *kk_800 = buffer.data(kk + 800);
    const auto *kk_801 = buffer.data(kk + 801);
    const auto *kk_802 = buffer.data(kk + 802);
    const auto *kk_803 = buffer.data(kk + 803);
    const auto *kk_804 = buffer.data(kk + 804);
    const auto *kk_805 = buffer.data(kk + 805);
    const auto *kk_806 = buffer.data(kk + 806);
    const auto *kk_807 = buffer.data(kk + 807);
    const auto *kk_808 = buffer.data(kk + 808);
    const auto *kk_809 = buffer.data(kk + 809);
    const auto *kk_810 = buffer.data(kk + 810);
    const auto *kk_811 = buffer.data(kk + 811);
    const auto *kk_812 = buffer.data(kk + 812);
    const auto *kk_813 = buffer.data(kk + 813);
    const auto *kk_814 = buffer.data(kk + 814);
    const auto *kk_815 = buffer.data(kk + 815);
    const auto *kk_816 = buffer.data(kk + 816);
    const auto *kk_817 = buffer.data(kk + 817);
    const auto *kk_818 = buffer.data(kk + 818);
    const auto *kk_819 = buffer.data(kk + 819);
    const auto *kk_820 = buffer.data(kk + 820);
    const auto *kk_821 = buffer.data(kk + 821);
    const auto *kk_822 = buffer.data(kk + 822);
    const auto *kk_823 = buffer.data(kk + 823);
    const auto *kk_824 = buffer.data(kk + 824);
    const auto *kk_825 = buffer.data(kk + 825);
    const auto *kk_826 = buffer.data(kk + 826);
    const auto *kk_827 = buffer.data(kk + 827);
    const auto *kk_828 = buffer.data(kk + 828);
    const auto *kk_829 = buffer.data(kk + 829);
    const auto *kk_830 = buffer.data(kk + 830);
    const auto *kk_831 = buffer.data(kk + 831);
    const auto *kk_832 = buffer.data(kk + 832);
    const auto *kk_833 = buffer.data(kk + 833);
    const auto *kk_834 = buffer.data(kk + 834);
    const auto *kk_835 = buffer.data(kk + 835);
    const auto *kk_836 = buffer.data(kk + 836);
    const auto *kk_837 = buffer.data(kk + 837);
    const auto *kk_838 = buffer.data(kk + 838);
    const auto *kk_839 = buffer.data(kk + 839);
    const auto *kk_840 = buffer.data(kk + 840);
    const auto *kk_841 = buffer.data(kk + 841);
    const auto *kk_842 = buffer.data(kk + 842);
    const auto *kk_843 = buffer.data(kk + 843);
    const auto *kk_844 = buffer.data(kk + 844);
    const auto *kk_845 = buffer.data(kk + 845);
    const auto *kk_846 = buffer.data(kk + 846);
    const auto *kk_847 = buffer.data(kk + 847);
    const auto *kk_848 = buffer.data(kk + 848);
    const auto *kk_849 = buffer.data(kk + 849);
    const auto *kk_850 = buffer.data(kk + 850);
    const auto *kk_851 = buffer.data(kk + 851);
    const auto *kk_852 = buffer.data(kk + 852);
    const auto *kk_853 = buffer.data(kk + 853);
    const auto *kk_854 = buffer.data(kk + 854);
    const auto *kk_855 = buffer.data(kk + 855);
    const auto *kk_856 = buffer.data(kk + 856);
    const auto *kk_857 = buffer.data(kk + 857);
    const auto *kk_858 = buffer.data(kk + 858);
    const auto *kk_859 = buffer.data(kk + 859);
    const auto *kk_860 = buffer.data(kk + 860);
    const auto *kk_861 = buffer.data(kk + 861);
    const auto *kk_862 = buffer.data(kk + 862);
    const auto *kk_863 = buffer.data(kk + 863);
    const auto *kk_864 = buffer.data(kk + 864);
    const auto *kk_865 = buffer.data(kk + 865);
    const auto *kk_866 = buffer.data(kk + 866);
    const auto *kk_867 = buffer.data(kk + 867);
    const auto *kk_868 = buffer.data(kk + 868);
    const auto *kk_869 = buffer.data(kk + 869);
    const auto *kk_870 = buffer.data(kk + 870);
    const auto *kk_871 = buffer.data(kk + 871);
    const auto *kk_872 = buffer.data(kk + 872);

#pragma omp simd aligned(t_494, t_495, t_496, t_497, t_498, hk_350, hk_351, hk_352, hk_353, \
                         hk_354, kk_674, kk_675, kk_676, kk_677, \
                         kk_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = -hk_350[k]
                   + f_0 * kk_674[k];

        t_495[k] = -hk_351[k]
                   + f_0 * kk_675[k];

        t_496[k] = -hk_352[k]
                   + f_0 * kk_676[k];

        t_497[k] = -hk_353[k]
                   + f_0 * kk_677[k];

        t_498[k] = -hk_354[k]
                   + f_0 * kk_678[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, t_503, hk_355, hk_356, hk_357, hk_358, \
                         hk_359, kk_679, kk_680, kk_681, kk_682, \
                         kk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = -hk_355[k]
                   + f_0 * kk_679[k];

        t_500[k] = -hk_356[k]
                   + f_0 * kk_680[k];

        t_501[k] = -hk_357[k]
                   + f_0 * kk_681[k];

        t_502[k] = -hk_358[k]
                   + f_0 * kk_682[k];

        t_503[k] = -hk_359[k]
                   + f_0 * kk_683[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, t_509, t_510, t_511, kk_684, \
                         kk_685, kk_686, kk_687, kk_688, kk_689, kk_690, \
                         kk_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_0 * kk_684[k];

        t_505[k] = f_0 * kk_685[k];

        t_506[k] = f_0 * kk_686[k];

        t_507[k] = f_0 * kk_687[k];

        t_508[k] = f_0 * kk_688[k];

        t_509[k] = f_0 * kk_689[k];

        t_510[k] = f_0 * kk_690[k];

        t_511[k] = f_0 * kk_691[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, t_517, t_518, t_519, kk_692, \
                         kk_693, kk_694, kk_695, kk_696, kk_697, kk_698, \
                         kk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_0 * kk_692[k];

        t_513[k] = f_0 * kk_693[k];

        t_514[k] = f_0 * kk_694[k];

        t_515[k] = f_0 * kk_695[k];

        t_516[k] = f_0 * kk_696[k];

        t_517[k] = f_0 * kk_697[k];

        t_518[k] = f_0 * kk_698[k];

        t_519[k] = f_0 * kk_699[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, t_525, t_526, t_527, kk_700, \
                         kk_701, kk_702, kk_703, kk_704, kk_705, kk_706, \
                         kk_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_0 * kk_700[k];

        t_521[k] = f_0 * kk_701[k];

        t_522[k] = f_0 * kk_702[k];

        t_523[k] = f_0 * kk_703[k];

        t_524[k] = f_0 * kk_704[k];

        t_525[k] = f_0 * kk_705[k];

        t_526[k] = f_0 * kk_706[k];

        t_527[k] = f_0 * kk_707[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, t_533, t_534, t_535, kk_708, \
                         kk_709, kk_710, kk_711, kk_712, kk_713, kk_714, \
                         kk_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_0 * kk_708[k];

        t_529[k] = f_0 * kk_709[k];

        t_530[k] = f_0 * kk_710[k];

        t_531[k] = f_0 * kk_711[k];

        t_532[k] = f_0 * kk_712[k];

        t_533[k] = f_0 * kk_713[k];

        t_534[k] = f_0 * kk_714[k];

        t_535[k] = f_0 * kk_715[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, t_541, hk_360, hk_361, kk_716, \
                         kk_717, kk_718, kk_719, kk_756, kk_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_0 * kk_716[k];

        t_537[k] = f_0 * kk_717[k];

        t_538[k] = f_0 * kk_718[k];

        t_539[k] = f_0 * kk_719[k];

        t_540[k] = -5.0 * hk_360[k]
                   + f_0 * kk_756[k];

        t_541[k] = -5.0 * hk_361[k]
                   + f_0 * kk_757[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, hk_362, hk_363, hk_364, hk_365, \
                         hk_366, kk_758, kk_759, kk_760, kk_761, \
                         kk_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = -5.0 * hk_362[k]
                   + f_0 * kk_758[k];

        t_543[k] = -5.0 * hk_363[k]
                   + f_0 * kk_759[k];

        t_544[k] = -5.0 * hk_364[k]
                   + f_0 * kk_760[k];

        t_545[k] = -5.0 * hk_365[k]
                   + f_0 * kk_761[k];

        t_546[k] = -5.0 * hk_366[k]
                   + f_0 * kk_762[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, hk_367, hk_368, hk_369, hk_370, \
                         hk_371, kk_763, kk_764, kk_765, kk_766, \
                         kk_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = -5.0 * hk_367[k]
                   + f_0 * kk_763[k];

        t_548[k] = -5.0 * hk_368[k]
                   + f_0 * kk_764[k];

        t_549[k] = -5.0 * hk_369[k]
                   + f_0 * kk_765[k];

        t_550[k] = -5.0 * hk_370[k]
                   + f_0 * kk_766[k];

        t_551[k] = -5.0 * hk_371[k]
                   + f_0 * kk_767[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, t_556, hk_372, hk_373, hk_374, hk_375, \
                         hk_376, kk_768, kk_769, kk_770, kk_771, \
                         kk_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = -5.0 * hk_372[k]
                   + f_0 * kk_768[k];

        t_553[k] = -5.0 * hk_373[k]
                   + f_0 * kk_769[k];

        t_554[k] = -5.0 * hk_374[k]
                   + f_0 * kk_770[k];

        t_555[k] = -5.0 * hk_375[k]
                   + f_0 * kk_771[k];

        t_556[k] = -5.0 * hk_376[k]
                   + f_0 * kk_772[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, hk_377, hk_378, hk_379, hk_380, \
                         hk_381, kk_773, kk_774, kk_775, kk_776, \
                         kk_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = -5.0 * hk_377[k]
                   + f_0 * kk_773[k];

        t_558[k] = -5.0 * hk_378[k]
                   + f_0 * kk_774[k];

        t_559[k] = -5.0 * hk_379[k]
                   + f_0 * kk_775[k];

        t_560[k] = -5.0 * hk_380[k]
                   + f_0 * kk_776[k];

        t_561[k] = -5.0 * hk_381[k]
                   + f_0 * kk_777[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, hk_382, hk_383, hk_384, hk_385, \
                         hk_386, kk_778, kk_779, kk_780, kk_781, \
                         kk_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = -5.0 * hk_382[k]
                   + f_0 * kk_778[k];

        t_563[k] = -5.0 * hk_383[k]
                   + f_0 * kk_779[k];

        t_564[k] = -5.0 * hk_384[k]
                   + f_0 * kk_780[k];

        t_565[k] = -5.0 * hk_385[k]
                   + f_0 * kk_781[k];

        t_566[k] = -5.0 * hk_386[k]
                   + f_0 * kk_782[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, hk_387, hk_388, hk_389, hk_390, \
                         hk_391, kk_783, kk_784, kk_785, kk_786, \
                         kk_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -5.0 * hk_387[k]
                   + f_0 * kk_783[k];

        t_568[k] = -5.0 * hk_388[k]
                   + f_0 * kk_784[k];

        t_569[k] = -5.0 * hk_389[k]
                   + f_0 * kk_785[k];

        t_570[k] = -5.0 * hk_390[k]
                   + f_0 * kk_786[k];

        t_571[k] = -5.0 * hk_391[k]
                   + f_0 * kk_787[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, hk_392, hk_393, hk_394, hk_395, \
                         hk_396, kk_788, kk_789, kk_790, kk_791, \
                         kk_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = -5.0 * hk_392[k]
                   + f_0 * kk_788[k];

        t_573[k] = -5.0 * hk_393[k]
                   + f_0 * kk_789[k];

        t_574[k] = -5.0 * hk_394[k]
                   + f_0 * kk_790[k];

        t_575[k] = -5.0 * hk_395[k]
                   + f_0 * kk_791[k];

        t_576[k] = -4.0 * hk_396[k]
                   + f_0 * kk_792[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, hk_397, hk_398, hk_399, hk_400, \
                         hk_401, kk_793, kk_794, kk_795, kk_796, \
                         kk_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = -4.0 * hk_397[k]
                   + f_0 * kk_793[k];

        t_578[k] = -4.0 * hk_398[k]
                   + f_0 * kk_794[k];

        t_579[k] = -4.0 * hk_399[k]
                   + f_0 * kk_795[k];

        t_580[k] = -4.0 * hk_400[k]
                   + f_0 * kk_796[k];

        t_581[k] = -4.0 * hk_401[k]
                   + f_0 * kk_797[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, hk_402, hk_403, hk_404, hk_405, \
                         hk_406, kk_798, kk_799, kk_800, kk_801, \
                         kk_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = -4.0 * hk_402[k]
                   + f_0 * kk_798[k];

        t_583[k] = -4.0 * hk_403[k]
                   + f_0 * kk_799[k];

        t_584[k] = -4.0 * hk_404[k]
                   + f_0 * kk_800[k];

        t_585[k] = -4.0 * hk_405[k]
                   + f_0 * kk_801[k];

        t_586[k] = -4.0 * hk_406[k]
                   + f_0 * kk_802[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, t_591, hk_407, hk_408, hk_409, hk_410, \
                         hk_411, kk_803, kk_804, kk_805, kk_806, \
                         kk_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = -4.0 * hk_407[k]
                   + f_0 * kk_803[k];

        t_588[k] = -4.0 * hk_408[k]
                   + f_0 * kk_804[k];

        t_589[k] = -4.0 * hk_409[k]
                   + f_0 * kk_805[k];

        t_590[k] = -4.0 * hk_410[k]
                   + f_0 * kk_806[k];

        t_591[k] = -4.0 * hk_411[k]
                   + f_0 * kk_807[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, t_596, hk_412, hk_413, hk_414, hk_415, \
                         hk_416, kk_808, kk_809, kk_810, kk_811, \
                         kk_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = -4.0 * hk_412[k]
                   + f_0 * kk_808[k];

        t_593[k] = -4.0 * hk_413[k]
                   + f_0 * kk_809[k];

        t_594[k] = -4.0 * hk_414[k]
                   + f_0 * kk_810[k];

        t_595[k] = -4.0 * hk_415[k]
                   + f_0 * kk_811[k];

        t_596[k] = -4.0 * hk_416[k]
                   + f_0 * kk_812[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, t_600, t_601, hk_417, hk_418, hk_419, hk_420, \
                         hk_421, kk_813, kk_814, kk_815, kk_816, \
                         kk_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = -4.0 * hk_417[k]
                   + f_0 * kk_813[k];

        t_598[k] = -4.0 * hk_418[k]
                   + f_0 * kk_814[k];

        t_599[k] = -4.0 * hk_419[k]
                   + f_0 * kk_815[k];

        t_600[k] = -4.0 * hk_420[k]
                   + f_0 * kk_816[k];

        t_601[k] = -4.0 * hk_421[k]
                   + f_0 * kk_817[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, hk_422, hk_423, hk_424, hk_425, \
                         hk_426, kk_818, kk_819, kk_820, kk_821, \
                         kk_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = -4.0 * hk_422[k]
                   + f_0 * kk_818[k];

        t_603[k] = -4.0 * hk_423[k]
                   + f_0 * kk_819[k];

        t_604[k] = -4.0 * hk_424[k]
                   + f_0 * kk_820[k];

        t_605[k] = -4.0 * hk_425[k]
                   + f_0 * kk_821[k];

        t_606[k] = -4.0 * hk_426[k]
                   + f_0 * kk_822[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, t_611, hk_427, hk_428, hk_429, hk_430, \
                         hk_431, kk_823, kk_824, kk_825, kk_826, \
                         kk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = -4.0 * hk_427[k]
                   + f_0 * kk_823[k];

        t_608[k] = -4.0 * hk_428[k]
                   + f_0 * kk_824[k];

        t_609[k] = -4.0 * hk_429[k]
                   + f_0 * kk_825[k];

        t_610[k] = -4.0 * hk_430[k]
                   + f_0 * kk_826[k];

        t_611[k] = -4.0 * hk_431[k]
                   + f_0 * kk_827[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, hk_432, hk_433, hk_434, hk_435, \
                         hk_436, kk_828, kk_829, kk_830, kk_831, \
                         kk_832 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = -3.0 * hk_432[k]
                   + f_0 * kk_828[k];

        t_613[k] = -3.0 * hk_433[k]
                   + f_0 * kk_829[k];

        t_614[k] = -3.0 * hk_434[k]
                   + f_0 * kk_830[k];

        t_615[k] = -3.0 * hk_435[k]
                   + f_0 * kk_831[k];

        t_616[k] = -3.0 * hk_436[k]
                   + f_0 * kk_832[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, hk_437, hk_438, hk_439, hk_440, \
                         hk_441, kk_833, kk_834, kk_835, kk_836, \
                         kk_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = -3.0 * hk_437[k]
                   + f_0 * kk_833[k];

        t_618[k] = -3.0 * hk_438[k]
                   + f_0 * kk_834[k];

        t_619[k] = -3.0 * hk_439[k]
                   + f_0 * kk_835[k];

        t_620[k] = -3.0 * hk_440[k]
                   + f_0 * kk_836[k];

        t_621[k] = -3.0 * hk_441[k]
                   + f_0 * kk_837[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, hk_442, hk_443, hk_444, hk_445, \
                         hk_446, kk_838, kk_839, kk_840, kk_841, \
                         kk_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = -3.0 * hk_442[k]
                   + f_0 * kk_838[k];

        t_623[k] = -3.0 * hk_443[k]
                   + f_0 * kk_839[k];

        t_624[k] = -3.0 * hk_444[k]
                   + f_0 * kk_840[k];

        t_625[k] = -3.0 * hk_445[k]
                   + f_0 * kk_841[k];

        t_626[k] = -3.0 * hk_446[k]
                   + f_0 * kk_842[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, t_631, hk_447, hk_448, hk_449, hk_450, \
                         hk_451, kk_843, kk_844, kk_845, kk_846, \
                         kk_847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = -3.0 * hk_447[k]
                   + f_0 * kk_843[k];

        t_628[k] = -3.0 * hk_448[k]
                   + f_0 * kk_844[k];

        t_629[k] = -3.0 * hk_449[k]
                   + f_0 * kk_845[k];

        t_630[k] = -3.0 * hk_450[k]
                   + f_0 * kk_846[k];

        t_631[k] = -3.0 * hk_451[k]
                   + f_0 * kk_847[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, t_636, hk_452, hk_453, hk_454, hk_455, \
                         hk_456, kk_848, kk_849, kk_850, kk_851, \
                         kk_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = -3.0 * hk_452[k]
                   + f_0 * kk_848[k];

        t_633[k] = -3.0 * hk_453[k]
                   + f_0 * kk_849[k];

        t_634[k] = -3.0 * hk_454[k]
                   + f_0 * kk_850[k];

        t_635[k] = -3.0 * hk_455[k]
                   + f_0 * kk_851[k];

        t_636[k] = -3.0 * hk_456[k]
                   + f_0 * kk_852[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, t_641, hk_457, hk_458, hk_459, hk_460, \
                         hk_461, kk_853, kk_854, kk_855, kk_856, \
                         kk_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = -3.0 * hk_457[k]
                   + f_0 * kk_853[k];

        t_638[k] = -3.0 * hk_458[k]
                   + f_0 * kk_854[k];

        t_639[k] = -3.0 * hk_459[k]
                   + f_0 * kk_855[k];

        t_640[k] = -3.0 * hk_460[k]
                   + f_0 * kk_856[k];

        t_641[k] = -3.0 * hk_461[k]
                   + f_0 * kk_857[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, t_646, hk_462, hk_463, hk_464, hk_465, \
                         hk_466, kk_858, kk_859, kk_860, kk_861, \
                         kk_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = -3.0 * hk_462[k]
                   + f_0 * kk_858[k];

        t_643[k] = -3.0 * hk_463[k]
                   + f_0 * kk_859[k];

        t_644[k] = -3.0 * hk_464[k]
                   + f_0 * kk_860[k];

        t_645[k] = -3.0 * hk_465[k]
                   + f_0 * kk_861[k];

        t_646[k] = -3.0 * hk_466[k]
                   + f_0 * kk_862[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, t_651, hk_467, hk_468, hk_469, hk_470, \
                         hk_471, kk_863, kk_864, kk_865, kk_866, \
                         kk_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = -3.0 * hk_467[k]
                   + f_0 * kk_863[k];

        t_648[k] = -2.0 * hk_468[k]
                   + f_0 * kk_864[k];

        t_649[k] = -2.0 * hk_469[k]
                   + f_0 * kk_865[k];

        t_650[k] = -2.0 * hk_470[k]
                   + f_0 * kk_866[k];

        t_651[k] = -2.0 * hk_471[k]
                   + f_0 * kk_867[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, t_656, hk_472, hk_473, hk_474, hk_475, \
                         hk_476, kk_868, kk_869, kk_870, kk_871, \
                         kk_872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = -2.0 * hk_472[k]
                   + f_0 * kk_868[k];

        t_653[k] = -2.0 * hk_473[k]
                   + f_0 * kk_869[k];

        t_654[k] = -2.0 * hk_474[k]
                   + f_0 * kk_870[k];

        t_655[k] = -2.0 * hk_475[k]
                   + f_0 * kk_871[k];

        t_656[k] = -2.0 * hk_476[k]
                   + f_0 * kk_872[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_477 = buffer.data(hk + 477);
    const auto *hk_478 = buffer.data(hk + 478);
    const auto *hk_479 = buffer.data(hk + 479);
    const auto *hk_480 = buffer.data(hk + 480);
    const auto *hk_481 = buffer.data(hk + 481);
    const auto *hk_482 = buffer.data(hk + 482);
    const auto *hk_483 = buffer.data(hk + 483);
    const auto *hk_484 = buffer.data(hk + 484);
    const auto *hk_485 = buffer.data(hk + 485);
    const auto *hk_486 = buffer.data(hk + 486);
    const auto *hk_487 = buffer.data(hk + 487);
    const auto *hk_488 = buffer.data(hk + 488);
    const auto *hk_489 = buffer.data(hk + 489);
    const auto *hk_490 = buffer.data(hk + 490);
    const auto *hk_491 = buffer.data(hk + 491);
    const auto *hk_492 = buffer.data(hk + 492);
    const auto *hk_493 = buffer.data(hk + 493);
    const auto *hk_494 = buffer.data(hk + 494);
    const auto *hk_495 = buffer.data(hk + 495);
    const auto *hk_496 = buffer.data(hk + 496);
    const auto *hk_497 = buffer.data(hk + 497);
    const auto *hk_498 = buffer.data(hk + 498);
    const auto *hk_499 = buffer.data(hk + 499);
    const auto *hk_500 = buffer.data(hk + 500);
    const auto *hk_501 = buffer.data(hk + 501);
    const auto *hk_502 = buffer.data(hk + 502);
    const auto *hk_503 = buffer.data(hk + 503);
    const auto *hk_504 = buffer.data(hk + 504);
    const auto *hk_505 = buffer.data(hk + 505);
    const auto *hk_506 = buffer.data(hk + 506);
    const auto *hk_507 = buffer.data(hk + 507);
    const auto *hk_508 = buffer.data(hk + 508);
    const auto *hk_509 = buffer.data(hk + 509);
    const auto *hk_510 = buffer.data(hk + 510);
    const auto *hk_511 = buffer.data(hk + 511);
    const auto *hk_512 = buffer.data(hk + 512);
    const auto *hk_513 = buffer.data(hk + 513);
    const auto *hk_514 = buffer.data(hk + 514);
    const auto *hk_515 = buffer.data(hk + 515);
    const auto *hk_516 = buffer.data(hk + 516);
    const auto *hk_517 = buffer.data(hk + 517);
    const auto *hk_518 = buffer.data(hk + 518);
    const auto *hk_519 = buffer.data(hk + 519);
    const auto *hk_520 = buffer.data(hk + 520);
    const auto *hk_521 = buffer.data(hk + 521);
    const auto *hk_522 = buffer.data(hk + 522);
    const auto *hk_523 = buffer.data(hk + 523);
    const auto *hk_524 = buffer.data(hk + 524);
    const auto *hk_525 = buffer.data(hk + 525);
    const auto *hk_526 = buffer.data(hk + 526);
    const auto *hk_527 = buffer.data(hk + 527);
    const auto *hk_528 = buffer.data(hk + 528);
    const auto *hk_529 = buffer.data(hk + 529);
    const auto *hk_530 = buffer.data(hk + 530);
    const auto *hk_531 = buffer.data(hk + 531);
    const auto *hk_532 = buffer.data(hk + 532);
    const auto *hk_533 = buffer.data(hk + 533);
    const auto *hk_534 = buffer.data(hk + 534);
    const auto *hk_535 = buffer.data(hk + 535);
    const auto *hk_536 = buffer.data(hk + 536);
    const auto *hk_537 = buffer.data(hk + 537);
    const auto *hk_538 = buffer.data(hk + 538);
    const auto *hk_539 = buffer.data(hk + 539);
    const auto *hk_540 = buffer.data(hk + 540);
    const auto *hk_541 = buffer.data(hk + 541);
    const auto *hk_542 = buffer.data(hk + 542);
    const auto *hk_543 = buffer.data(hk + 543);
    const auto *hk_544 = buffer.data(hk + 544);
    const auto *hk_545 = buffer.data(hk + 545);
    const auto *hk_546 = buffer.data(hk + 546);
    const auto *hk_547 = buffer.data(hk + 547);
    const auto *hk_548 = buffer.data(hk + 548);
    const auto *hk_549 = buffer.data(hk + 549);
    const auto *hk_550 = buffer.data(hk + 550);
    const auto *hk_551 = buffer.data(hk + 551);
    const auto *hk_552 = buffer.data(hk + 552);
    const auto *hk_553 = buffer.data(hk + 553);
    const auto *hk_554 = buffer.data(hk + 554);
    const auto *hk_555 = buffer.data(hk + 555);
    const auto *hk_556 = buffer.data(hk + 556);
    const auto *hk_557 = buffer.data(hk + 557);
    const auto *hk_558 = buffer.data(hk + 558);
    const auto *hk_559 = buffer.data(hk + 559);
    const auto *hk_560 = buffer.data(hk + 560);
    const auto *hk_561 = buffer.data(hk + 561);
    const auto *hk_562 = buffer.data(hk + 562);
    const auto *hk_563 = buffer.data(hk + 563);
    const auto *hk_564 = buffer.data(hk + 564);
    const auto *hk_565 = buffer.data(hk + 565);
    const auto *hk_566 = buffer.data(hk + 566);
    const auto *hk_567 = buffer.data(hk + 567);
    const auto *hk_568 = buffer.data(hk + 568);
    const auto *hk_569 = buffer.data(hk + 569);
    const auto *hk_570 = buffer.data(hk + 570);
    const auto *hk_571 = buffer.data(hk + 571);
    const auto *hk_572 = buffer.data(hk + 572);
    const auto *hk_573 = buffer.data(hk + 573);
    const auto *hk_574 = buffer.data(hk + 574);
    const auto *hk_575 = buffer.data(hk + 575);
    const auto *hk_576 = buffer.data(hk + 576);
    const auto *hk_577 = buffer.data(hk + 577);
    const auto *hk_578 = buffer.data(hk + 578);
    const auto *hk_579 = buffer.data(hk + 579);
    const auto *hk_580 = buffer.data(hk + 580);
    const auto *hk_581 = buffer.data(hk + 581);
    const auto *hk_582 = buffer.data(hk + 582);
    const auto *hk_583 = buffer.data(hk + 583);
    const auto *hk_584 = buffer.data(hk + 584);
    const auto *hk_585 = buffer.data(hk + 585);
    const auto *hk_586 = buffer.data(hk + 586);
    const auto *hk_587 = buffer.data(hk + 587);
    const auto *hk_588 = buffer.data(hk + 588);
    const auto *hk_589 = buffer.data(hk + 589);
    const auto *hk_590 = buffer.data(hk + 590);
    const auto *hk_591 = buffer.data(hk + 591);
    const auto *hk_592 = buffer.data(hk + 592);
    const auto *hk_593 = buffer.data(hk + 593);
    const auto *hk_594 = buffer.data(hk + 594);
    const auto *hk_595 = buffer.data(hk + 595);
    const auto *hk_596 = buffer.data(hk + 596);
    const auto *hk_597 = buffer.data(hk + 597);
    const auto *hk_598 = buffer.data(hk + 598);
    const auto *hk_599 = buffer.data(hk + 599);
    const auto *hk_600 = buffer.data(hk + 600);
    const auto *hk_601 = buffer.data(hk + 601);
    const auto *hk_602 = buffer.data(hk + 602);
    const auto *hk_603 = buffer.data(hk + 603);

    const auto *kk_873 = buffer.data(kk + 873);
    const auto *kk_874 = buffer.data(kk + 874);
    const auto *kk_875 = buffer.data(kk + 875);
    const auto *kk_876 = buffer.data(kk + 876);
    const auto *kk_877 = buffer.data(kk + 877);
    const auto *kk_878 = buffer.data(kk + 878);
    const auto *kk_879 = buffer.data(kk + 879);
    const auto *kk_880 = buffer.data(kk + 880);
    const auto *kk_881 = buffer.data(kk + 881);
    const auto *kk_882 = buffer.data(kk + 882);
    const auto *kk_883 = buffer.data(kk + 883);
    const auto *kk_884 = buffer.data(kk + 884);
    const auto *kk_885 = buffer.data(kk + 885);
    const auto *kk_886 = buffer.data(kk + 886);
    const auto *kk_887 = buffer.data(kk + 887);
    const auto *kk_888 = buffer.data(kk + 888);
    const auto *kk_889 = buffer.data(kk + 889);
    const auto *kk_890 = buffer.data(kk + 890);
    const auto *kk_891 = buffer.data(kk + 891);
    const auto *kk_892 = buffer.data(kk + 892);
    const auto *kk_893 = buffer.data(kk + 893);
    const auto *kk_894 = buffer.data(kk + 894);
    const auto *kk_895 = buffer.data(kk + 895);
    const auto *kk_896 = buffer.data(kk + 896);
    const auto *kk_897 = buffer.data(kk + 897);
    const auto *kk_898 = buffer.data(kk + 898);
    const auto *kk_899 = buffer.data(kk + 899);
    const auto *kk_900 = buffer.data(kk + 900);
    const auto *kk_901 = buffer.data(kk + 901);
    const auto *kk_902 = buffer.data(kk + 902);
    const auto *kk_903 = buffer.data(kk + 903);
    const auto *kk_904 = buffer.data(kk + 904);
    const auto *kk_905 = buffer.data(kk + 905);
    const auto *kk_906 = buffer.data(kk + 906);
    const auto *kk_907 = buffer.data(kk + 907);
    const auto *kk_908 = buffer.data(kk + 908);
    const auto *kk_909 = buffer.data(kk + 909);
    const auto *kk_910 = buffer.data(kk + 910);
    const auto *kk_911 = buffer.data(kk + 911);
    const auto *kk_912 = buffer.data(kk + 912);
    const auto *kk_913 = buffer.data(kk + 913);
    const auto *kk_914 = buffer.data(kk + 914);
    const auto *kk_915 = buffer.data(kk + 915);
    const auto *kk_916 = buffer.data(kk + 916);
    const auto *kk_917 = buffer.data(kk + 917);
    const auto *kk_918 = buffer.data(kk + 918);
    const auto *kk_919 = buffer.data(kk + 919);
    const auto *kk_920 = buffer.data(kk + 920);
    const auto *kk_921 = buffer.data(kk + 921);
    const auto *kk_922 = buffer.data(kk + 922);
    const auto *kk_923 = buffer.data(kk + 923);
    const auto *kk_924 = buffer.data(kk + 924);
    const auto *kk_925 = buffer.data(kk + 925);
    const auto *kk_926 = buffer.data(kk + 926);
    const auto *kk_927 = buffer.data(kk + 927);
    const auto *kk_928 = buffer.data(kk + 928);
    const auto *kk_929 = buffer.data(kk + 929);
    const auto *kk_930 = buffer.data(kk + 930);
    const auto *kk_931 = buffer.data(kk + 931);
    const auto *kk_932 = buffer.data(kk + 932);
    const auto *kk_933 = buffer.data(kk + 933);
    const auto *kk_934 = buffer.data(kk + 934);
    const auto *kk_935 = buffer.data(kk + 935);
    const auto *kk_936 = buffer.data(kk + 936);
    const auto *kk_937 = buffer.data(kk + 937);
    const auto *kk_938 = buffer.data(kk + 938);
    const auto *kk_939 = buffer.data(kk + 939);
    const auto *kk_940 = buffer.data(kk + 940);
    const auto *kk_941 = buffer.data(kk + 941);
    const auto *kk_942 = buffer.data(kk + 942);
    const auto *kk_943 = buffer.data(kk + 943);
    const auto *kk_944 = buffer.data(kk + 944);
    const auto *kk_945 = buffer.data(kk + 945);
    const auto *kk_946 = buffer.data(kk + 946);
    const auto *kk_947 = buffer.data(kk + 947);
    const auto *kk_948 = buffer.data(kk + 948);
    const auto *kk_949 = buffer.data(kk + 949);
    const auto *kk_950 = buffer.data(kk + 950);
    const auto *kk_951 = buffer.data(kk + 951);
    const auto *kk_952 = buffer.data(kk + 952);
    const auto *kk_953 = buffer.data(kk + 953);
    const auto *kk_954 = buffer.data(kk + 954);
    const auto *kk_955 = buffer.data(kk + 955);
    const auto *kk_956 = buffer.data(kk + 956);
    const auto *kk_957 = buffer.data(kk + 957);
    const auto *kk_958 = buffer.data(kk + 958);
    const auto *kk_959 = buffer.data(kk + 959);
    const auto *kk_960 = buffer.data(kk + 960);
    const auto *kk_961 = buffer.data(kk + 961);
    const auto *kk_962 = buffer.data(kk + 962);
    const auto *kk_963 = buffer.data(kk + 963);
    const auto *kk_964 = buffer.data(kk + 964);
    const auto *kk_965 = buffer.data(kk + 965);
    const auto *kk_966 = buffer.data(kk + 966);
    const auto *kk_967 = buffer.data(kk + 967);
    const auto *kk_968 = buffer.data(kk + 968);
    const auto *kk_969 = buffer.data(kk + 969);
    const auto *kk_970 = buffer.data(kk + 970);
    const auto *kk_971 = buffer.data(kk + 971);
    const auto *kk_1008 = buffer.data(kk + 1008);
    const auto *kk_1009 = buffer.data(kk + 1009);
    const auto *kk_1010 = buffer.data(kk + 1010);
    const auto *kk_1011 = buffer.data(kk + 1011);
    const auto *kk_1012 = buffer.data(kk + 1012);
    const auto *kk_1013 = buffer.data(kk + 1013);
    const auto *kk_1014 = buffer.data(kk + 1014);
    const auto *kk_1015 = buffer.data(kk + 1015);
    const auto *kk_1016 = buffer.data(kk + 1016);
    const auto *kk_1017 = buffer.data(kk + 1017);
    const auto *kk_1018 = buffer.data(kk + 1018);
    const auto *kk_1019 = buffer.data(kk + 1019);
    const auto *kk_1020 = buffer.data(kk + 1020);
    const auto *kk_1021 = buffer.data(kk + 1021);
    const auto *kk_1022 = buffer.data(kk + 1022);
    const auto *kk_1023 = buffer.data(kk + 1023);
    const auto *kk_1024 = buffer.data(kk + 1024);
    const auto *kk_1025 = buffer.data(kk + 1025);
    const auto *kk_1026 = buffer.data(kk + 1026);
    const auto *kk_1027 = buffer.data(kk + 1027);
    const auto *kk_1028 = buffer.data(kk + 1028);
    const auto *kk_1029 = buffer.data(kk + 1029);
    const auto *kk_1030 = buffer.data(kk + 1030);
    const auto *kk_1031 = buffer.data(kk + 1031);
    const auto *kk_1032 = buffer.data(kk + 1032);
    const auto *kk_1033 = buffer.data(kk + 1033);
    const auto *kk_1034 = buffer.data(kk + 1034);
    const auto *kk_1035 = buffer.data(kk + 1035);
    const auto *kk_1036 = buffer.data(kk + 1036);
    const auto *kk_1037 = buffer.data(kk + 1037);
    const auto *kk_1038 = buffer.data(kk + 1038);
    const auto *kk_1039 = buffer.data(kk + 1039);
    const auto *kk_1040 = buffer.data(kk + 1040);
    const auto *kk_1041 = buffer.data(kk + 1041);
    const auto *kk_1042 = buffer.data(kk + 1042);
    const auto *kk_1043 = buffer.data(kk + 1043);
    const auto *kk_1044 = buffer.data(kk + 1044);
    const auto *kk_1045 = buffer.data(kk + 1045);
    const auto *kk_1046 = buffer.data(kk + 1046);
    const auto *kk_1047 = buffer.data(kk + 1047);
    const auto *kk_1048 = buffer.data(kk + 1048);
    const auto *kk_1049 = buffer.data(kk + 1049);
    const auto *kk_1050 = buffer.data(kk + 1050);
    const auto *kk_1051 = buffer.data(kk + 1051);
    const auto *kk_1052 = buffer.data(kk + 1052);
    const auto *kk_1053 = buffer.data(kk + 1053);
    const auto *kk_1054 = buffer.data(kk + 1054);
    const auto *kk_1055 = buffer.data(kk + 1055);
    const auto *kk_1056 = buffer.data(kk + 1056);
    const auto *kk_1057 = buffer.data(kk + 1057);
    const auto *kk_1058 = buffer.data(kk + 1058);
    const auto *kk_1059 = buffer.data(kk + 1059);
    const auto *kk_1060 = buffer.data(kk + 1060);
    const auto *kk_1061 = buffer.data(kk + 1061);
    const auto *kk_1062 = buffer.data(kk + 1062);
    const auto *kk_1063 = buffer.data(kk + 1063);
    const auto *kk_1064 = buffer.data(kk + 1064);
    const auto *kk_1065 = buffer.data(kk + 1065);
    const auto *kk_1066 = buffer.data(kk + 1066);
    const auto *kk_1067 = buffer.data(kk + 1067);
    const auto *kk_1068 = buffer.data(kk + 1068);
    const auto *kk_1069 = buffer.data(kk + 1069);
    const auto *kk_1070 = buffer.data(kk + 1070);
    const auto *kk_1071 = buffer.data(kk + 1071);

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, hk_477, hk_478, hk_479, hk_480, \
                         hk_481, kk_873, kk_874, kk_875, kk_876, \
                         kk_877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = -2.0 * hk_477[k]
                   + f_0 * kk_873[k];

        t_658[k] = -2.0 * hk_478[k]
                   + f_0 * kk_874[k];

        t_659[k] = -2.0 * hk_479[k]
                   + f_0 * kk_875[k];

        t_660[k] = -2.0 * hk_480[k]
                   + f_0 * kk_876[k];

        t_661[k] = -2.0 * hk_481[k]
                   + f_0 * kk_877[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, t_666, hk_482, hk_483, hk_484, hk_485, \
                         hk_486, kk_878, kk_879, kk_880, kk_881, \
                         kk_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = -2.0 * hk_482[k]
                   + f_0 * kk_878[k];

        t_663[k] = -2.0 * hk_483[k]
                   + f_0 * kk_879[k];

        t_664[k] = -2.0 * hk_484[k]
                   + f_0 * kk_880[k];

        t_665[k] = -2.0 * hk_485[k]
                   + f_0 * kk_881[k];

        t_666[k] = -2.0 * hk_486[k]
                   + f_0 * kk_882[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, t_670, t_671, hk_487, hk_488, hk_489, hk_490, \
                         hk_491, kk_883, kk_884, kk_885, kk_886, \
                         kk_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = -2.0 * hk_487[k]
                   + f_0 * kk_883[k];

        t_668[k] = -2.0 * hk_488[k]
                   + f_0 * kk_884[k];

        t_669[k] = -2.0 * hk_489[k]
                   + f_0 * kk_885[k];

        t_670[k] = -2.0 * hk_490[k]
                   + f_0 * kk_886[k];

        t_671[k] = -2.0 * hk_491[k]
                   + f_0 * kk_887[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, t_676, hk_492, hk_493, hk_494, hk_495, \
                         hk_496, kk_888, kk_889, kk_890, kk_891, \
                         kk_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = -2.0 * hk_492[k]
                   + f_0 * kk_888[k];

        t_673[k] = -2.0 * hk_493[k]
                   + f_0 * kk_889[k];

        t_674[k] = -2.0 * hk_494[k]
                   + f_0 * kk_890[k];

        t_675[k] = -2.0 * hk_495[k]
                   + f_0 * kk_891[k];

        t_676[k] = -2.0 * hk_496[k]
                   + f_0 * kk_892[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, t_680, t_681, hk_497, hk_498, hk_499, hk_500, \
                         hk_501, kk_893, kk_894, kk_895, kk_896, \
                         kk_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = -2.0 * hk_497[k]
                   + f_0 * kk_893[k];

        t_678[k] = -2.0 * hk_498[k]
                   + f_0 * kk_894[k];

        t_679[k] = -2.0 * hk_499[k]
                   + f_0 * kk_895[k];

        t_680[k] = -2.0 * hk_500[k]
                   + f_0 * kk_896[k];

        t_681[k] = -2.0 * hk_501[k]
                   + f_0 * kk_897[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, t_686, hk_502, hk_503, hk_504, hk_505, \
                         hk_506, kk_898, kk_899, kk_900, kk_901, \
                         kk_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = -2.0 * hk_502[k]
                   + f_0 * kk_898[k];

        t_683[k] = -2.0 * hk_503[k]
                   + f_0 * kk_899[k];

        t_684[k] = -hk_504[k]
                   + f_0 * kk_900[k];

        t_685[k] = -hk_505[k]
                   + f_0 * kk_901[k];

        t_686[k] = -hk_506[k]
                   + f_0 * kk_902[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, t_691, hk_507, hk_508, hk_509, hk_510, \
                         hk_511, kk_903, kk_904, kk_905, kk_906, \
                         kk_907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = -hk_507[k]
                   + f_0 * kk_903[k];

        t_688[k] = -hk_508[k]
                   + f_0 * kk_904[k];

        t_689[k] = -hk_509[k]
                   + f_0 * kk_905[k];

        t_690[k] = -hk_510[k]
                   + f_0 * kk_906[k];

        t_691[k] = -hk_511[k]
                   + f_0 * kk_907[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, t_696, hk_512, hk_513, hk_514, hk_515, \
                         hk_516, kk_908, kk_909, kk_910, kk_911, \
                         kk_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = -hk_512[k]
                   + f_0 * kk_908[k];

        t_693[k] = -hk_513[k]
                   + f_0 * kk_909[k];

        t_694[k] = -hk_514[k]
                   + f_0 * kk_910[k];

        t_695[k] = -hk_515[k]
                   + f_0 * kk_911[k];

        t_696[k] = -hk_516[k]
                   + f_0 * kk_912[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, t_701, hk_517, hk_518, hk_519, hk_520, \
                         hk_521, kk_913, kk_914, kk_915, kk_916, \
                         kk_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = -hk_517[k]
                   + f_0 * kk_913[k];

        t_698[k] = -hk_518[k]
                   + f_0 * kk_914[k];

        t_699[k] = -hk_519[k]
                   + f_0 * kk_915[k];

        t_700[k] = -hk_520[k]
                   + f_0 * kk_916[k];

        t_701[k] = -hk_521[k]
                   + f_0 * kk_917[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, t_706, hk_522, hk_523, hk_524, hk_525, \
                         hk_526, kk_918, kk_919, kk_920, kk_921, \
                         kk_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = -hk_522[k]
                   + f_0 * kk_918[k];

        t_703[k] = -hk_523[k]
                   + f_0 * kk_919[k];

        t_704[k] = -hk_524[k]
                   + f_0 * kk_920[k];

        t_705[k] = -hk_525[k]
                   + f_0 * kk_921[k];

        t_706[k] = -hk_526[k]
                   + f_0 * kk_922[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, t_710, t_711, hk_527, hk_528, hk_529, hk_530, \
                         hk_531, kk_923, kk_924, kk_925, kk_926, \
                         kk_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = -hk_527[k]
                   + f_0 * kk_923[k];

        t_708[k] = -hk_528[k]
                   + f_0 * kk_924[k];

        t_709[k] = -hk_529[k]
                   + f_0 * kk_925[k];

        t_710[k] = -hk_530[k]
                   + f_0 * kk_926[k];

        t_711[k] = -hk_531[k]
                   + f_0 * kk_927[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, t_716, hk_532, hk_533, hk_534, hk_535, \
                         hk_536, kk_928, kk_929, kk_930, kk_931, \
                         kk_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = -hk_532[k]
                   + f_0 * kk_928[k];

        t_713[k] = -hk_533[k]
                   + f_0 * kk_929[k];

        t_714[k] = -hk_534[k]
                   + f_0 * kk_930[k];

        t_715[k] = -hk_535[k]
                   + f_0 * kk_931[k];

        t_716[k] = -hk_536[k]
                   + f_0 * kk_932[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, t_720, t_721, t_722, hk_537, hk_538, hk_539, \
                         kk_933, kk_934, kk_935, kk_936, kk_937, \
                         kk_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = -hk_537[k]
                   + f_0 * kk_933[k];

        t_718[k] = -hk_538[k]
                   + f_0 * kk_934[k];

        t_719[k] = -hk_539[k]
                   + f_0 * kk_935[k];

        t_720[k] = f_0 * kk_936[k];

        t_721[k] = f_0 * kk_937[k];

        t_722[k] = f_0 * kk_938[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, t_726, t_727, t_728, t_729, t_730, kk_939, \
                         kk_940, kk_941, kk_942, kk_943, kk_944, kk_945, \
                         kk_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_0 * kk_939[k];

        t_724[k] = f_0 * kk_940[k];

        t_725[k] = f_0 * kk_941[k];

        t_726[k] = f_0 * kk_942[k];

        t_727[k] = f_0 * kk_943[k];

        t_728[k] = f_0 * kk_944[k];

        t_729[k] = f_0 * kk_945[k];

        t_730[k] = f_0 * kk_946[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, t_734, t_735, t_736, t_737, t_738, kk_947, \
                         kk_948, kk_949, kk_950, kk_951, kk_952, kk_953, \
                         kk_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_0 * kk_947[k];

        t_732[k] = f_0 * kk_948[k];

        t_733[k] = f_0 * kk_949[k];

        t_734[k] = f_0 * kk_950[k];

        t_735[k] = f_0 * kk_951[k];

        t_736[k] = f_0 * kk_952[k];

        t_737[k] = f_0 * kk_953[k];

        t_738[k] = f_0 * kk_954[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, t_743, t_744, t_745, t_746, kk_955, \
                         kk_956, kk_957, kk_958, kk_959, kk_960, kk_961, \
                         kk_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_0 * kk_955[k];

        t_740[k] = f_0 * kk_956[k];

        t_741[k] = f_0 * kk_957[k];

        t_742[k] = f_0 * kk_958[k];

        t_743[k] = f_0 * kk_959[k];

        t_744[k] = f_0 * kk_960[k];

        t_745[k] = f_0 * kk_961[k];

        t_746[k] = f_0 * kk_962[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, t_750, t_751, t_752, t_753, t_754, kk_963, \
                         kk_964, kk_965, kk_966, kk_967, kk_968, kk_969, \
                         kk_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_0 * kk_963[k];

        t_748[k] = f_0 * kk_964[k];

        t_749[k] = f_0 * kk_965[k];

        t_750[k] = f_0 * kk_966[k];

        t_751[k] = f_0 * kk_967[k];

        t_752[k] = f_0 * kk_968[k];

        t_753[k] = f_0 * kk_969[k];

        t_754[k] = f_0 * kk_970[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, hk_540, hk_541, hk_542, hk_543, \
                         kk_971, kk_1008, kk_1009, kk_1010, kk_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_0 * kk_971[k];

        t_756[k] = -6.0 * hk_540[k]
                   + f_0 * kk_1008[k];

        t_757[k] = -6.0 * hk_541[k]
                   + f_0 * kk_1009[k];

        t_758[k] = -6.0 * hk_542[k]
                   + f_0 * kk_1010[k];

        t_759[k] = -6.0 * hk_543[k]
                   + f_0 * kk_1011[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, t_764, hk_544, hk_545, hk_546, hk_547, \
                         hk_548, kk_1012, kk_1013, kk_1014, kk_1015, \
                         kk_1016 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = -6.0 * hk_544[k]
                   + f_0 * kk_1012[k];

        t_761[k] = -6.0 * hk_545[k]
                   + f_0 * kk_1013[k];

        t_762[k] = -6.0 * hk_546[k]
                   + f_0 * kk_1014[k];

        t_763[k] = -6.0 * hk_547[k]
                   + f_0 * kk_1015[k];

        t_764[k] = -6.0 * hk_548[k]
                   + f_0 * kk_1016[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, t_769, hk_549, hk_550, hk_551, hk_552, \
                         hk_553, kk_1017, kk_1018, kk_1019, kk_1020, \
                         kk_1021 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = -6.0 * hk_549[k]
                   + f_0 * kk_1017[k];

        t_766[k] = -6.0 * hk_550[k]
                   + f_0 * kk_1018[k];

        t_767[k] = -6.0 * hk_551[k]
                   + f_0 * kk_1019[k];

        t_768[k] = -6.0 * hk_552[k]
                   + f_0 * kk_1020[k];

        t_769[k] = -6.0 * hk_553[k]
                   + f_0 * kk_1021[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, hk_554, hk_555, hk_556, hk_557, \
                         hk_558, kk_1022, kk_1023, kk_1024, kk_1025, \
                         kk_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = -6.0 * hk_554[k]
                   + f_0 * kk_1022[k];

        t_771[k] = -6.0 * hk_555[k]
                   + f_0 * kk_1023[k];

        t_772[k] = -6.0 * hk_556[k]
                   + f_0 * kk_1024[k];

        t_773[k] = -6.0 * hk_557[k]
                   + f_0 * kk_1025[k];

        t_774[k] = -6.0 * hk_558[k]
                   + f_0 * kk_1026[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, hk_559, hk_560, hk_561, hk_562, \
                         hk_563, kk_1027, kk_1028, kk_1029, kk_1030, \
                         kk_1031 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = -6.0 * hk_559[k]
                   + f_0 * kk_1027[k];

        t_776[k] = -6.0 * hk_560[k]
                   + f_0 * kk_1028[k];

        t_777[k] = -6.0 * hk_561[k]
                   + f_0 * kk_1029[k];

        t_778[k] = -6.0 * hk_562[k]
                   + f_0 * kk_1030[k];

        t_779[k] = -6.0 * hk_563[k]
                   + f_0 * kk_1031[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, t_784, hk_564, hk_565, hk_566, hk_567, \
                         hk_568, kk_1032, kk_1033, kk_1034, kk_1035, \
                         kk_1036 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = -6.0 * hk_564[k]
                   + f_0 * kk_1032[k];

        t_781[k] = -6.0 * hk_565[k]
                   + f_0 * kk_1033[k];

        t_782[k] = -6.0 * hk_566[k]
                   + f_0 * kk_1034[k];

        t_783[k] = -6.0 * hk_567[k]
                   + f_0 * kk_1035[k];

        t_784[k] = -6.0 * hk_568[k]
                   + f_0 * kk_1036[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, t_789, hk_569, hk_570, hk_571, hk_572, \
                         hk_573, kk_1037, kk_1038, kk_1039, kk_1040, \
                         kk_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = -6.0 * hk_569[k]
                   + f_0 * kk_1037[k];

        t_786[k] = -6.0 * hk_570[k]
                   + f_0 * kk_1038[k];

        t_787[k] = -6.0 * hk_571[k]
                   + f_0 * kk_1039[k];

        t_788[k] = -6.0 * hk_572[k]
                   + f_0 * kk_1040[k];

        t_789[k] = -6.0 * hk_573[k]
                   + f_0 * kk_1041[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, t_794, hk_574, hk_575, hk_576, hk_577, \
                         hk_578, kk_1042, kk_1043, kk_1044, kk_1045, \
                         kk_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = -6.0 * hk_574[k]
                   + f_0 * kk_1042[k];

        t_791[k] = -6.0 * hk_575[k]
                   + f_0 * kk_1043[k];

        t_792[k] = -5.0 * hk_576[k]
                   + f_0 * kk_1044[k];

        t_793[k] = -5.0 * hk_577[k]
                   + f_0 * kk_1045[k];

        t_794[k] = -5.0 * hk_578[k]
                   + f_0 * kk_1046[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, hk_579, hk_580, hk_581, hk_582, \
                         hk_583, kk_1047, kk_1048, kk_1049, kk_1050, \
                         kk_1051 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = -5.0 * hk_579[k]
                   + f_0 * kk_1047[k];

        t_796[k] = -5.0 * hk_580[k]
                   + f_0 * kk_1048[k];

        t_797[k] = -5.0 * hk_581[k]
                   + f_0 * kk_1049[k];

        t_798[k] = -5.0 * hk_582[k]
                   + f_0 * kk_1050[k];

        t_799[k] = -5.0 * hk_583[k]
                   + f_0 * kk_1051[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, hk_584, hk_585, hk_586, hk_587, \
                         hk_588, kk_1052, kk_1053, kk_1054, kk_1055, \
                         kk_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = -5.0 * hk_584[k]
                   + f_0 * kk_1052[k];

        t_801[k] = -5.0 * hk_585[k]
                   + f_0 * kk_1053[k];

        t_802[k] = -5.0 * hk_586[k]
                   + f_0 * kk_1054[k];

        t_803[k] = -5.0 * hk_587[k]
                   + f_0 * kk_1055[k];

        t_804[k] = -5.0 * hk_588[k]
                   + f_0 * kk_1056[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, hk_589, hk_590, hk_591, hk_592, \
                         hk_593, kk_1057, kk_1058, kk_1059, kk_1060, \
                         kk_1061 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = -5.0 * hk_589[k]
                   + f_0 * kk_1057[k];

        t_806[k] = -5.0 * hk_590[k]
                   + f_0 * kk_1058[k];

        t_807[k] = -5.0 * hk_591[k]
                   + f_0 * kk_1059[k];

        t_808[k] = -5.0 * hk_592[k]
                   + f_0 * kk_1060[k];

        t_809[k] = -5.0 * hk_593[k]
                   + f_0 * kk_1061[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, hk_594, hk_595, hk_596, hk_597, \
                         hk_598, kk_1062, kk_1063, kk_1064, kk_1065, \
                         kk_1066 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = -5.0 * hk_594[k]
                   + f_0 * kk_1062[k];

        t_811[k] = -5.0 * hk_595[k]
                   + f_0 * kk_1063[k];

        t_812[k] = -5.0 * hk_596[k]
                   + f_0 * kk_1064[k];

        t_813[k] = -5.0 * hk_597[k]
                   + f_0 * kk_1065[k];

        t_814[k] = -5.0 * hk_598[k]
                   + f_0 * kk_1066[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, hk_599, hk_600, hk_601, hk_602, \
                         hk_603, kk_1067, kk_1068, kk_1069, kk_1070, \
                         kk_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = -5.0 * hk_599[k]
                   + f_0 * kk_1067[k];

        t_816[k] = -5.0 * hk_600[k]
                   + f_0 * kk_1068[k];

        t_817[k] = -5.0 * hk_601[k]
                   + f_0 * kk_1069[k];

        t_818[k] = -5.0 * hk_602[k]
                   + f_0 * kk_1070[k];

        t_819[k] = -5.0 * hk_603[k]
                   + f_0 * kk_1071[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_1_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_604 = buffer.data(hk + 604);
    const auto *hk_605 = buffer.data(hk + 605);
    const auto *hk_606 = buffer.data(hk + 606);
    const auto *hk_607 = buffer.data(hk + 607);
    const auto *hk_608 = buffer.data(hk + 608);
    const auto *hk_609 = buffer.data(hk + 609);
    const auto *hk_610 = buffer.data(hk + 610);
    const auto *hk_611 = buffer.data(hk + 611);
    const auto *hk_612 = buffer.data(hk + 612);
    const auto *hk_613 = buffer.data(hk + 613);
    const auto *hk_614 = buffer.data(hk + 614);
    const auto *hk_615 = buffer.data(hk + 615);
    const auto *hk_616 = buffer.data(hk + 616);
    const auto *hk_617 = buffer.data(hk + 617);
    const auto *hk_618 = buffer.data(hk + 618);
    const auto *hk_619 = buffer.data(hk + 619);
    const auto *hk_620 = buffer.data(hk + 620);
    const auto *hk_621 = buffer.data(hk + 621);
    const auto *hk_622 = buffer.data(hk + 622);
    const auto *hk_623 = buffer.data(hk + 623);
    const auto *hk_624 = buffer.data(hk + 624);
    const auto *hk_625 = buffer.data(hk + 625);
    const auto *hk_626 = buffer.data(hk + 626);
    const auto *hk_627 = buffer.data(hk + 627);
    const auto *hk_628 = buffer.data(hk + 628);
    const auto *hk_629 = buffer.data(hk + 629);
    const auto *hk_630 = buffer.data(hk + 630);
    const auto *hk_631 = buffer.data(hk + 631);
    const auto *hk_632 = buffer.data(hk + 632);
    const auto *hk_633 = buffer.data(hk + 633);
    const auto *hk_634 = buffer.data(hk + 634);
    const auto *hk_635 = buffer.data(hk + 635);
    const auto *hk_636 = buffer.data(hk + 636);
    const auto *hk_637 = buffer.data(hk + 637);
    const auto *hk_638 = buffer.data(hk + 638);
    const auto *hk_639 = buffer.data(hk + 639);
    const auto *hk_640 = buffer.data(hk + 640);
    const auto *hk_641 = buffer.data(hk + 641);
    const auto *hk_642 = buffer.data(hk + 642);
    const auto *hk_643 = buffer.data(hk + 643);
    const auto *hk_644 = buffer.data(hk + 644);
    const auto *hk_645 = buffer.data(hk + 645);
    const auto *hk_646 = buffer.data(hk + 646);
    const auto *hk_647 = buffer.data(hk + 647);
    const auto *hk_648 = buffer.data(hk + 648);
    const auto *hk_649 = buffer.data(hk + 649);
    const auto *hk_650 = buffer.data(hk + 650);
    const auto *hk_651 = buffer.data(hk + 651);
    const auto *hk_652 = buffer.data(hk + 652);
    const auto *hk_653 = buffer.data(hk + 653);
    const auto *hk_654 = buffer.data(hk + 654);
    const auto *hk_655 = buffer.data(hk + 655);
    const auto *hk_656 = buffer.data(hk + 656);
    const auto *hk_657 = buffer.data(hk + 657);
    const auto *hk_658 = buffer.data(hk + 658);
    const auto *hk_659 = buffer.data(hk + 659);
    const auto *hk_660 = buffer.data(hk + 660);
    const auto *hk_661 = buffer.data(hk + 661);
    const auto *hk_662 = buffer.data(hk + 662);
    const auto *hk_663 = buffer.data(hk + 663);
    const auto *hk_664 = buffer.data(hk + 664);
    const auto *hk_665 = buffer.data(hk + 665);
    const auto *hk_666 = buffer.data(hk + 666);
    const auto *hk_667 = buffer.data(hk + 667);
    const auto *hk_668 = buffer.data(hk + 668);
    const auto *hk_669 = buffer.data(hk + 669);
    const auto *hk_670 = buffer.data(hk + 670);
    const auto *hk_671 = buffer.data(hk + 671);
    const auto *hk_672 = buffer.data(hk + 672);
    const auto *hk_673 = buffer.data(hk + 673);
    const auto *hk_674 = buffer.data(hk + 674);
    const auto *hk_675 = buffer.data(hk + 675);
    const auto *hk_676 = buffer.data(hk + 676);
    const auto *hk_677 = buffer.data(hk + 677);
    const auto *hk_678 = buffer.data(hk + 678);
    const auto *hk_679 = buffer.data(hk + 679);
    const auto *hk_680 = buffer.data(hk + 680);
    const auto *hk_681 = buffer.data(hk + 681);
    const auto *hk_682 = buffer.data(hk + 682);
    const auto *hk_683 = buffer.data(hk + 683);
    const auto *hk_684 = buffer.data(hk + 684);
    const auto *hk_685 = buffer.data(hk + 685);
    const auto *hk_686 = buffer.data(hk + 686);
    const auto *hk_687 = buffer.data(hk + 687);
    const auto *hk_688 = buffer.data(hk + 688);
    const auto *hk_689 = buffer.data(hk + 689);
    const auto *hk_690 = buffer.data(hk + 690);
    const auto *hk_691 = buffer.data(hk + 691);
    const auto *hk_692 = buffer.data(hk + 692);
    const auto *hk_693 = buffer.data(hk + 693);
    const auto *hk_694 = buffer.data(hk + 694);
    const auto *hk_695 = buffer.data(hk + 695);
    const auto *hk_696 = buffer.data(hk + 696);
    const auto *hk_697 = buffer.data(hk + 697);
    const auto *hk_698 = buffer.data(hk + 698);
    const auto *hk_699 = buffer.data(hk + 699);
    const auto *hk_700 = buffer.data(hk + 700);
    const auto *hk_701 = buffer.data(hk + 701);
    const auto *hk_702 = buffer.data(hk + 702);
    const auto *hk_703 = buffer.data(hk + 703);
    const auto *hk_704 = buffer.data(hk + 704);
    const auto *hk_705 = buffer.data(hk + 705);
    const auto *hk_706 = buffer.data(hk + 706);
    const auto *hk_707 = buffer.data(hk + 707);
    const auto *hk_708 = buffer.data(hk + 708);
    const auto *hk_709 = buffer.data(hk + 709);
    const auto *hk_710 = buffer.data(hk + 710);
    const auto *hk_711 = buffer.data(hk + 711);
    const auto *hk_712 = buffer.data(hk + 712);
    const auto *hk_713 = buffer.data(hk + 713);
    const auto *hk_714 = buffer.data(hk + 714);
    const auto *hk_715 = buffer.data(hk + 715);
    const auto *hk_716 = buffer.data(hk + 716);
    const auto *hk_717 = buffer.data(hk + 717);
    const auto *hk_718 = buffer.data(hk + 718);
    const auto *hk_719 = buffer.data(hk + 719);
    const auto *hk_720 = buffer.data(hk + 720);
    const auto *hk_721 = buffer.data(hk + 721);
    const auto *hk_722 = buffer.data(hk + 722);
    const auto *hk_723 = buffer.data(hk + 723);
    const auto *hk_724 = buffer.data(hk + 724);
    const auto *hk_725 = buffer.data(hk + 725);
    const auto *hk_726 = buffer.data(hk + 726);
    const auto *hk_727 = buffer.data(hk + 727);
    const auto *hk_728 = buffer.data(hk + 728);
    const auto *hk_729 = buffer.data(hk + 729);
    const auto *hk_730 = buffer.data(hk + 730);
    const auto *hk_731 = buffer.data(hk + 731);
    const auto *hk_732 = buffer.data(hk + 732);
    const auto *hk_733 = buffer.data(hk + 733);
    const auto *hk_734 = buffer.data(hk + 734);
    const auto *hk_735 = buffer.data(hk + 735);
    const auto *hk_736 = buffer.data(hk + 736);
    const auto *hk_737 = buffer.data(hk + 737);
    const auto *hk_738 = buffer.data(hk + 738);
    const auto *hk_739 = buffer.data(hk + 739);
    const auto *hk_740 = buffer.data(hk + 740);
    const auto *hk_741 = buffer.data(hk + 741);
    const auto *hk_742 = buffer.data(hk + 742);
    const auto *hk_743 = buffer.data(hk + 743);
    const auto *hk_744 = buffer.data(hk + 744);
    const auto *hk_745 = buffer.data(hk + 745);
    const auto *hk_746 = buffer.data(hk + 746);
    const auto *hk_747 = buffer.data(hk + 747);
    const auto *hk_748 = buffer.data(hk + 748);
    const auto *hk_749 = buffer.data(hk + 749);
    const auto *hk_750 = buffer.data(hk + 750);
    const auto *hk_751 = buffer.data(hk + 751);
    const auto *hk_752 = buffer.data(hk + 752);
    const auto *hk_753 = buffer.data(hk + 753);

    const auto *kk_1072 = buffer.data(kk + 1072);
    const auto *kk_1073 = buffer.data(kk + 1073);
    const auto *kk_1074 = buffer.data(kk + 1074);
    const auto *kk_1075 = buffer.data(kk + 1075);
    const auto *kk_1076 = buffer.data(kk + 1076);
    const auto *kk_1077 = buffer.data(kk + 1077);
    const auto *kk_1078 = buffer.data(kk + 1078);
    const auto *kk_1079 = buffer.data(kk + 1079);
    const auto *kk_1080 = buffer.data(kk + 1080);
    const auto *kk_1081 = buffer.data(kk + 1081);
    const auto *kk_1082 = buffer.data(kk + 1082);
    const auto *kk_1083 = buffer.data(kk + 1083);
    const auto *kk_1084 = buffer.data(kk + 1084);
    const auto *kk_1085 = buffer.data(kk + 1085);
    const auto *kk_1086 = buffer.data(kk + 1086);
    const auto *kk_1087 = buffer.data(kk + 1087);
    const auto *kk_1088 = buffer.data(kk + 1088);
    const auto *kk_1089 = buffer.data(kk + 1089);
    const auto *kk_1090 = buffer.data(kk + 1090);
    const auto *kk_1091 = buffer.data(kk + 1091);
    const auto *kk_1092 = buffer.data(kk + 1092);
    const auto *kk_1093 = buffer.data(kk + 1093);
    const auto *kk_1094 = buffer.data(kk + 1094);
    const auto *kk_1095 = buffer.data(kk + 1095);
    const auto *kk_1096 = buffer.data(kk + 1096);
    const auto *kk_1097 = buffer.data(kk + 1097);
    const auto *kk_1098 = buffer.data(kk + 1098);
    const auto *kk_1099 = buffer.data(kk + 1099);
    const auto *kk_1100 = buffer.data(kk + 1100);
    const auto *kk_1101 = buffer.data(kk + 1101);
    const auto *kk_1102 = buffer.data(kk + 1102);
    const auto *kk_1103 = buffer.data(kk + 1103);
    const auto *kk_1104 = buffer.data(kk + 1104);
    const auto *kk_1105 = buffer.data(kk + 1105);
    const auto *kk_1106 = buffer.data(kk + 1106);
    const auto *kk_1107 = buffer.data(kk + 1107);
    const auto *kk_1108 = buffer.data(kk + 1108);
    const auto *kk_1109 = buffer.data(kk + 1109);
    const auto *kk_1110 = buffer.data(kk + 1110);
    const auto *kk_1111 = buffer.data(kk + 1111);
    const auto *kk_1112 = buffer.data(kk + 1112);
    const auto *kk_1113 = buffer.data(kk + 1113);
    const auto *kk_1114 = buffer.data(kk + 1114);
    const auto *kk_1115 = buffer.data(kk + 1115);
    const auto *kk_1116 = buffer.data(kk + 1116);
    const auto *kk_1117 = buffer.data(kk + 1117);
    const auto *kk_1118 = buffer.data(kk + 1118);
    const auto *kk_1119 = buffer.data(kk + 1119);
    const auto *kk_1120 = buffer.data(kk + 1120);
    const auto *kk_1121 = buffer.data(kk + 1121);
    const auto *kk_1122 = buffer.data(kk + 1122);
    const auto *kk_1123 = buffer.data(kk + 1123);
    const auto *kk_1124 = buffer.data(kk + 1124);
    const auto *kk_1125 = buffer.data(kk + 1125);
    const auto *kk_1126 = buffer.data(kk + 1126);
    const auto *kk_1127 = buffer.data(kk + 1127);
    const auto *kk_1128 = buffer.data(kk + 1128);
    const auto *kk_1129 = buffer.data(kk + 1129);
    const auto *kk_1130 = buffer.data(kk + 1130);
    const auto *kk_1131 = buffer.data(kk + 1131);
    const auto *kk_1132 = buffer.data(kk + 1132);
    const auto *kk_1133 = buffer.data(kk + 1133);
    const auto *kk_1134 = buffer.data(kk + 1134);
    const auto *kk_1135 = buffer.data(kk + 1135);
    const auto *kk_1136 = buffer.data(kk + 1136);
    const auto *kk_1137 = buffer.data(kk + 1137);
    const auto *kk_1138 = buffer.data(kk + 1138);
    const auto *kk_1139 = buffer.data(kk + 1139);
    const auto *kk_1140 = buffer.data(kk + 1140);
    const auto *kk_1141 = buffer.data(kk + 1141);
    const auto *kk_1142 = buffer.data(kk + 1142);
    const auto *kk_1143 = buffer.data(kk + 1143);
    const auto *kk_1144 = buffer.data(kk + 1144);
    const auto *kk_1145 = buffer.data(kk + 1145);
    const auto *kk_1146 = buffer.data(kk + 1146);
    const auto *kk_1147 = buffer.data(kk + 1147);
    const auto *kk_1148 = buffer.data(kk + 1148);
    const auto *kk_1149 = buffer.data(kk + 1149);
    const auto *kk_1150 = buffer.data(kk + 1150);
    const auto *kk_1151 = buffer.data(kk + 1151);
    const auto *kk_1152 = buffer.data(kk + 1152);
    const auto *kk_1153 = buffer.data(kk + 1153);
    const auto *kk_1154 = buffer.data(kk + 1154);
    const auto *kk_1155 = buffer.data(kk + 1155);
    const auto *kk_1156 = buffer.data(kk + 1156);
    const auto *kk_1157 = buffer.data(kk + 1157);
    const auto *kk_1158 = buffer.data(kk + 1158);
    const auto *kk_1159 = buffer.data(kk + 1159);
    const auto *kk_1160 = buffer.data(kk + 1160);
    const auto *kk_1161 = buffer.data(kk + 1161);
    const auto *kk_1162 = buffer.data(kk + 1162);
    const auto *kk_1163 = buffer.data(kk + 1163);
    const auto *kk_1164 = buffer.data(kk + 1164);
    const auto *kk_1165 = buffer.data(kk + 1165);
    const auto *kk_1166 = buffer.data(kk + 1166);
    const auto *kk_1167 = buffer.data(kk + 1167);
    const auto *kk_1168 = buffer.data(kk + 1168);
    const auto *kk_1169 = buffer.data(kk + 1169);
    const auto *kk_1170 = buffer.data(kk + 1170);
    const auto *kk_1171 = buffer.data(kk + 1171);
    const auto *kk_1172 = buffer.data(kk + 1172);
    const auto *kk_1173 = buffer.data(kk + 1173);
    const auto *kk_1174 = buffer.data(kk + 1174);
    const auto *kk_1175 = buffer.data(kk + 1175);
    const auto *kk_1176 = buffer.data(kk + 1176);
    const auto *kk_1177 = buffer.data(kk + 1177);
    const auto *kk_1178 = buffer.data(kk + 1178);
    const auto *kk_1179 = buffer.data(kk + 1179);
    const auto *kk_1180 = buffer.data(kk + 1180);
    const auto *kk_1181 = buffer.data(kk + 1181);
    const auto *kk_1182 = buffer.data(kk + 1182);
    const auto *kk_1183 = buffer.data(kk + 1183);
    const auto *kk_1184 = buffer.data(kk + 1184);
    const auto *kk_1185 = buffer.data(kk + 1185);
    const auto *kk_1186 = buffer.data(kk + 1186);
    const auto *kk_1187 = buffer.data(kk + 1187);
    const auto *kk_1188 = buffer.data(kk + 1188);
    const auto *kk_1189 = buffer.data(kk + 1189);
    const auto *kk_1190 = buffer.data(kk + 1190);
    const auto *kk_1191 = buffer.data(kk + 1191);
    const auto *kk_1192 = buffer.data(kk + 1192);
    const auto *kk_1193 = buffer.data(kk + 1193);
    const auto *kk_1194 = buffer.data(kk + 1194);
    const auto *kk_1195 = buffer.data(kk + 1195);
    const auto *kk_1196 = buffer.data(kk + 1196);
    const auto *kk_1197 = buffer.data(kk + 1197);
    const auto *kk_1198 = buffer.data(kk + 1198);
    const auto *kk_1199 = buffer.data(kk + 1199);
    const auto *kk_1200 = buffer.data(kk + 1200);
    const auto *kk_1201 = buffer.data(kk + 1201);
    const auto *kk_1202 = buffer.data(kk + 1202);
    const auto *kk_1203 = buffer.data(kk + 1203);
    const auto *kk_1204 = buffer.data(kk + 1204);
    const auto *kk_1205 = buffer.data(kk + 1205);
    const auto *kk_1206 = buffer.data(kk + 1206);
    const auto *kk_1207 = buffer.data(kk + 1207);
    const auto *kk_1208 = buffer.data(kk + 1208);
    const auto *kk_1209 = buffer.data(kk + 1209);
    const auto *kk_1210 = buffer.data(kk + 1210);
    const auto *kk_1211 = buffer.data(kk + 1211);
    const auto *kk_1212 = buffer.data(kk + 1212);
    const auto *kk_1213 = buffer.data(kk + 1213);
    const auto *kk_1214 = buffer.data(kk + 1214);
    const auto *kk_1215 = buffer.data(kk + 1215);
    const auto *kk_1216 = buffer.data(kk + 1216);
    const auto *kk_1217 = buffer.data(kk + 1217);
    const auto *kk_1218 = buffer.data(kk + 1218);
    const auto *kk_1219 = buffer.data(kk + 1219);
    const auto *kk_1220 = buffer.data(kk + 1220);
    const auto *kk_1221 = buffer.data(kk + 1221);

#pragma omp simd aligned(t_820, t_821, t_822, t_823, t_824, hk_604, hk_605, hk_606, hk_607, \
                         hk_608, kk_1072, kk_1073, kk_1074, kk_1075, \
                         kk_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = -5.0 * hk_604[k]
                   + f_0 * kk_1072[k];

        t_821[k] = -5.0 * hk_605[k]
                   + f_0 * kk_1073[k];

        t_822[k] = -5.0 * hk_606[k]
                   + f_0 * kk_1074[k];

        t_823[k] = -5.0 * hk_607[k]
                   + f_0 * kk_1075[k];

        t_824[k] = -5.0 * hk_608[k]
                   + f_0 * kk_1076[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, hk_609, hk_610, hk_611, hk_612, \
                         hk_613, kk_1077, kk_1078, kk_1079, kk_1080, \
                         kk_1081 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = -5.0 * hk_609[k]
                   + f_0 * kk_1077[k];

        t_826[k] = -5.0 * hk_610[k]
                   + f_0 * kk_1078[k];

        t_827[k] = -5.0 * hk_611[k]
                   + f_0 * kk_1079[k];

        t_828[k] = -4.0 * hk_612[k]
                   + f_0 * kk_1080[k];

        t_829[k] = -4.0 * hk_613[k]
                   + f_0 * kk_1081[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, hk_614, hk_615, hk_616, hk_617, \
                         hk_618, kk_1082, kk_1083, kk_1084, kk_1085, \
                         kk_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = -4.0 * hk_614[k]
                   + f_0 * kk_1082[k];

        t_831[k] = -4.0 * hk_615[k]
                   + f_0 * kk_1083[k];

        t_832[k] = -4.0 * hk_616[k]
                   + f_0 * kk_1084[k];

        t_833[k] = -4.0 * hk_617[k]
                   + f_0 * kk_1085[k];

        t_834[k] = -4.0 * hk_618[k]
                   + f_0 * kk_1086[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, t_838, t_839, hk_619, hk_620, hk_621, hk_622, \
                         hk_623, kk_1087, kk_1088, kk_1089, kk_1090, \
                         kk_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = -4.0 * hk_619[k]
                   + f_0 * kk_1087[k];

        t_836[k] = -4.0 * hk_620[k]
                   + f_0 * kk_1088[k];

        t_837[k] = -4.0 * hk_621[k]
                   + f_0 * kk_1089[k];

        t_838[k] = -4.0 * hk_622[k]
                   + f_0 * kk_1090[k];

        t_839[k] = -4.0 * hk_623[k]
                   + f_0 * kk_1091[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, hk_624, hk_625, hk_626, hk_627, \
                         hk_628, kk_1092, kk_1093, kk_1094, kk_1095, \
                         kk_1096 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = -4.0 * hk_624[k]
                   + f_0 * kk_1092[k];

        t_841[k] = -4.0 * hk_625[k]
                   + f_0 * kk_1093[k];

        t_842[k] = -4.0 * hk_626[k]
                   + f_0 * kk_1094[k];

        t_843[k] = -4.0 * hk_627[k]
                   + f_0 * kk_1095[k];

        t_844[k] = -4.0 * hk_628[k]
                   + f_0 * kk_1096[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, t_849, hk_629, hk_630, hk_631, hk_632, \
                         hk_633, kk_1097, kk_1098, kk_1099, kk_1100, \
                         kk_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = -4.0 * hk_629[k]
                   + f_0 * kk_1097[k];

        t_846[k] = -4.0 * hk_630[k]
                   + f_0 * kk_1098[k];

        t_847[k] = -4.0 * hk_631[k]
                   + f_0 * kk_1099[k];

        t_848[k] = -4.0 * hk_632[k]
                   + f_0 * kk_1100[k];

        t_849[k] = -4.0 * hk_633[k]
                   + f_0 * kk_1101[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, hk_634, hk_635, hk_636, hk_637, \
                         hk_638, kk_1102, kk_1103, kk_1104, kk_1105, \
                         kk_1106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = -4.0 * hk_634[k]
                   + f_0 * kk_1102[k];

        t_851[k] = -4.0 * hk_635[k]
                   + f_0 * kk_1103[k];

        t_852[k] = -4.0 * hk_636[k]
                   + f_0 * kk_1104[k];

        t_853[k] = -4.0 * hk_637[k]
                   + f_0 * kk_1105[k];

        t_854[k] = -4.0 * hk_638[k]
                   + f_0 * kk_1106[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, hk_639, hk_640, hk_641, hk_642, \
                         hk_643, kk_1107, kk_1108, kk_1109, kk_1110, \
                         kk_1111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = -4.0 * hk_639[k]
                   + f_0 * kk_1107[k];

        t_856[k] = -4.0 * hk_640[k]
                   + f_0 * kk_1108[k];

        t_857[k] = -4.0 * hk_641[k]
                   + f_0 * kk_1109[k];

        t_858[k] = -4.0 * hk_642[k]
                   + f_0 * kk_1110[k];

        t_859[k] = -4.0 * hk_643[k]
                   + f_0 * kk_1111[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, t_864, hk_644, hk_645, hk_646, hk_647, \
                         hk_648, kk_1112, kk_1113, kk_1114, kk_1115, \
                         kk_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = -4.0 * hk_644[k]
                   + f_0 * kk_1112[k];

        t_861[k] = -4.0 * hk_645[k]
                   + f_0 * kk_1113[k];

        t_862[k] = -4.0 * hk_646[k]
                   + f_0 * kk_1114[k];

        t_863[k] = -4.0 * hk_647[k]
                   + f_0 * kk_1115[k];

        t_864[k] = -3.0 * hk_648[k]
                   + f_0 * kk_1116[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, t_869, hk_649, hk_650, hk_651, hk_652, \
                         hk_653, kk_1117, kk_1118, kk_1119, kk_1120, \
                         kk_1121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = -3.0 * hk_649[k]
                   + f_0 * kk_1117[k];

        t_866[k] = -3.0 * hk_650[k]
                   + f_0 * kk_1118[k];

        t_867[k] = -3.0 * hk_651[k]
                   + f_0 * kk_1119[k];

        t_868[k] = -3.0 * hk_652[k]
                   + f_0 * kk_1120[k];

        t_869[k] = -3.0 * hk_653[k]
                   + f_0 * kk_1121[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, t_874, hk_654, hk_655, hk_656, hk_657, \
                         hk_658, kk_1122, kk_1123, kk_1124, kk_1125, \
                         kk_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = -3.0 * hk_654[k]
                   + f_0 * kk_1122[k];

        t_871[k] = -3.0 * hk_655[k]
                   + f_0 * kk_1123[k];

        t_872[k] = -3.0 * hk_656[k]
                   + f_0 * kk_1124[k];

        t_873[k] = -3.0 * hk_657[k]
                   + f_0 * kk_1125[k];

        t_874[k] = -3.0 * hk_658[k]
                   + f_0 * kk_1126[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, t_879, hk_659, hk_660, hk_661, hk_662, \
                         hk_663, kk_1127, kk_1128, kk_1129, kk_1130, \
                         kk_1131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = -3.0 * hk_659[k]
                   + f_0 * kk_1127[k];

        t_876[k] = -3.0 * hk_660[k]
                   + f_0 * kk_1128[k];

        t_877[k] = -3.0 * hk_661[k]
                   + f_0 * kk_1129[k];

        t_878[k] = -3.0 * hk_662[k]
                   + f_0 * kk_1130[k];

        t_879[k] = -3.0 * hk_663[k]
                   + f_0 * kk_1131[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, t_884, hk_664, hk_665, hk_666, hk_667, \
                         hk_668, kk_1132, kk_1133, kk_1134, kk_1135, \
                         kk_1136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = -3.0 * hk_664[k]
                   + f_0 * kk_1132[k];

        t_881[k] = -3.0 * hk_665[k]
                   + f_0 * kk_1133[k];

        t_882[k] = -3.0 * hk_666[k]
                   + f_0 * kk_1134[k];

        t_883[k] = -3.0 * hk_667[k]
                   + f_0 * kk_1135[k];

        t_884[k] = -3.0 * hk_668[k]
                   + f_0 * kk_1136[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, hk_669, hk_670, hk_671, hk_672, \
                         hk_673, kk_1137, kk_1138, kk_1139, kk_1140, \
                         kk_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = -3.0 * hk_669[k]
                   + f_0 * kk_1137[k];

        t_886[k] = -3.0 * hk_670[k]
                   + f_0 * kk_1138[k];

        t_887[k] = -3.0 * hk_671[k]
                   + f_0 * kk_1139[k];

        t_888[k] = -3.0 * hk_672[k]
                   + f_0 * kk_1140[k];

        t_889[k] = -3.0 * hk_673[k]
                   + f_0 * kk_1141[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, hk_674, hk_675, hk_676, hk_677, \
                         hk_678, kk_1142, kk_1143, kk_1144, kk_1145, \
                         kk_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = -3.0 * hk_674[k]
                   + f_0 * kk_1142[k];

        t_891[k] = -3.0 * hk_675[k]
                   + f_0 * kk_1143[k];

        t_892[k] = -3.0 * hk_676[k]
                   + f_0 * kk_1144[k];

        t_893[k] = -3.0 * hk_677[k]
                   + f_0 * kk_1145[k];

        t_894[k] = -3.0 * hk_678[k]
                   + f_0 * kk_1146[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, t_899, hk_679, hk_680, hk_681, hk_682, \
                         hk_683, kk_1147, kk_1148, kk_1149, kk_1150, \
                         kk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = -3.0 * hk_679[k]
                   + f_0 * kk_1147[k];

        t_896[k] = -3.0 * hk_680[k]
                   + f_0 * kk_1148[k];

        t_897[k] = -3.0 * hk_681[k]
                   + f_0 * kk_1149[k];

        t_898[k] = -3.0 * hk_682[k]
                   + f_0 * kk_1150[k];

        t_899[k] = -3.0 * hk_683[k]
                   + f_0 * kk_1151[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, t_904, hk_684, hk_685, hk_686, hk_687, \
                         hk_688, kk_1152, kk_1153, kk_1154, kk_1155, \
                         kk_1156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = -2.0 * hk_684[k]
                   + f_0 * kk_1152[k];

        t_901[k] = -2.0 * hk_685[k]
                   + f_0 * kk_1153[k];

        t_902[k] = -2.0 * hk_686[k]
                   + f_0 * kk_1154[k];

        t_903[k] = -2.0 * hk_687[k]
                   + f_0 * kk_1155[k];

        t_904[k] = -2.0 * hk_688[k]
                   + f_0 * kk_1156[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, t_909, hk_689, hk_690, hk_691, hk_692, \
                         hk_693, kk_1157, kk_1158, kk_1159, kk_1160, \
                         kk_1161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = -2.0 * hk_689[k]
                   + f_0 * kk_1157[k];

        t_906[k] = -2.0 * hk_690[k]
                   + f_0 * kk_1158[k];

        t_907[k] = -2.0 * hk_691[k]
                   + f_0 * kk_1159[k];

        t_908[k] = -2.0 * hk_692[k]
                   + f_0 * kk_1160[k];

        t_909[k] = -2.0 * hk_693[k]
                   + f_0 * kk_1161[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, t_914, hk_694, hk_695, hk_696, hk_697, \
                         hk_698, kk_1162, kk_1163, kk_1164, kk_1165, \
                         kk_1166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = -2.0 * hk_694[k]
                   + f_0 * kk_1162[k];

        t_911[k] = -2.0 * hk_695[k]
                   + f_0 * kk_1163[k];

        t_912[k] = -2.0 * hk_696[k]
                   + f_0 * kk_1164[k];

        t_913[k] = -2.0 * hk_697[k]
                   + f_0 * kk_1165[k];

        t_914[k] = -2.0 * hk_698[k]
                   + f_0 * kk_1166[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, t_919, hk_699, hk_700, hk_701, hk_702, \
                         hk_703, kk_1167, kk_1168, kk_1169, kk_1170, \
                         kk_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = -2.0 * hk_699[k]
                   + f_0 * kk_1167[k];

        t_916[k] = -2.0 * hk_700[k]
                   + f_0 * kk_1168[k];

        t_917[k] = -2.0 * hk_701[k]
                   + f_0 * kk_1169[k];

        t_918[k] = -2.0 * hk_702[k]
                   + f_0 * kk_1170[k];

        t_919[k] = -2.0 * hk_703[k]
                   + f_0 * kk_1171[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, t_924, hk_704, hk_705, hk_706, hk_707, \
                         hk_708, kk_1172, kk_1173, kk_1174, kk_1175, \
                         kk_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = -2.0 * hk_704[k]
                   + f_0 * kk_1172[k];

        t_921[k] = -2.0 * hk_705[k]
                   + f_0 * kk_1173[k];

        t_922[k] = -2.0 * hk_706[k]
                   + f_0 * kk_1174[k];

        t_923[k] = -2.0 * hk_707[k]
                   + f_0 * kk_1175[k];

        t_924[k] = -2.0 * hk_708[k]
                   + f_0 * kk_1176[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, t_929, hk_709, hk_710, hk_711, hk_712, \
                         hk_713, kk_1177, kk_1178, kk_1179, kk_1180, \
                         kk_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = -2.0 * hk_709[k]
                   + f_0 * kk_1177[k];

        t_926[k] = -2.0 * hk_710[k]
                   + f_0 * kk_1178[k];

        t_927[k] = -2.0 * hk_711[k]
                   + f_0 * kk_1179[k];

        t_928[k] = -2.0 * hk_712[k]
                   + f_0 * kk_1180[k];

        t_929[k] = -2.0 * hk_713[k]
                   + f_0 * kk_1181[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, hk_714, hk_715, hk_716, hk_717, \
                         hk_718, kk_1182, kk_1183, kk_1184, kk_1185, \
                         kk_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = -2.0 * hk_714[k]
                   + f_0 * kk_1182[k];

        t_931[k] = -2.0 * hk_715[k]
                   + f_0 * kk_1183[k];

        t_932[k] = -2.0 * hk_716[k]
                   + f_0 * kk_1184[k];

        t_933[k] = -2.0 * hk_717[k]
                   + f_0 * kk_1185[k];

        t_934[k] = -2.0 * hk_718[k]
                   + f_0 * kk_1186[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, t_939, hk_719, hk_720, hk_721, hk_722, \
                         hk_723, kk_1187, kk_1188, kk_1189, kk_1190, \
                         kk_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = -2.0 * hk_719[k]
                   + f_0 * kk_1187[k];

        t_936[k] = -hk_720[k]
                   + f_0 * kk_1188[k];

        t_937[k] = -hk_721[k]
                   + f_0 * kk_1189[k];

        t_938[k] = -hk_722[k]
                   + f_0 * kk_1190[k];

        t_939[k] = -hk_723[k]
                   + f_0 * kk_1191[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, t_944, hk_724, hk_725, hk_726, hk_727, \
                         hk_728, kk_1192, kk_1193, kk_1194, kk_1195, \
                         kk_1196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = -hk_724[k]
                   + f_0 * kk_1192[k];

        t_941[k] = -hk_725[k]
                   + f_0 * kk_1193[k];

        t_942[k] = -hk_726[k]
                   + f_0 * kk_1194[k];

        t_943[k] = -hk_727[k]
                   + f_0 * kk_1195[k];

        t_944[k] = -hk_728[k]
                   + f_0 * kk_1196[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, t_949, hk_729, hk_730, hk_731, hk_732, \
                         hk_733, kk_1197, kk_1198, kk_1199, kk_1200, \
                         kk_1201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = -hk_729[k]
                   + f_0 * kk_1197[k];

        t_946[k] = -hk_730[k]
                   + f_0 * kk_1198[k];

        t_947[k] = -hk_731[k]
                   + f_0 * kk_1199[k];

        t_948[k] = -hk_732[k]
                   + f_0 * kk_1200[k];

        t_949[k] = -hk_733[k]
                   + f_0 * kk_1201[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, t_954, hk_734, hk_735, hk_736, hk_737, \
                         hk_738, kk_1202, kk_1203, kk_1204, kk_1205, \
                         kk_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = -hk_734[k]
                   + f_0 * kk_1202[k];

        t_951[k] = -hk_735[k]
                   + f_0 * kk_1203[k];

        t_952[k] = -hk_736[k]
                   + f_0 * kk_1204[k];

        t_953[k] = -hk_737[k]
                   + f_0 * kk_1205[k];

        t_954[k] = -hk_738[k]
                   + f_0 * kk_1206[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, t_959, hk_739, hk_740, hk_741, hk_742, \
                         hk_743, kk_1207, kk_1208, kk_1209, kk_1210, \
                         kk_1211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = -hk_739[k]
                   + f_0 * kk_1207[k];

        t_956[k] = -hk_740[k]
                   + f_0 * kk_1208[k];

        t_957[k] = -hk_741[k]
                   + f_0 * kk_1209[k];

        t_958[k] = -hk_742[k]
                   + f_0 * kk_1210[k];

        t_959[k] = -hk_743[k]
                   + f_0 * kk_1211[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, t_963, t_964, hk_744, hk_745, hk_746, hk_747, \
                         hk_748, kk_1212, kk_1213, kk_1214, kk_1215, \
                         kk_1216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = -hk_744[k]
                   + f_0 * kk_1212[k];

        t_961[k] = -hk_745[k]
                   + f_0 * kk_1213[k];

        t_962[k] = -hk_746[k]
                   + f_0 * kk_1214[k];

        t_963[k] = -hk_747[k]
                   + f_0 * kk_1215[k];

        t_964[k] = -hk_748[k]
                   + f_0 * kk_1216[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, t_969, hk_749, hk_750, hk_751, hk_752, \
                         hk_753, kk_1217, kk_1218, kk_1219, kk_1220, \
                         kk_1221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = -hk_749[k]
                   + f_0 * kk_1217[k];

        t_966[k] = -hk_750[k]
                   + f_0 * kk_1218[k];

        t_967[k] = -hk_751[k]
                   + f_0 * kk_1219[k];

        t_968[k] = -hk_752[k]
                   + f_0 * kk_1220[k];

        t_969[k] = -hk_753[k]
                   + f_0 * kk_1221[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_1_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_754 = buffer.data(hk + 754);
    const auto *hk_755 = buffer.data(hk + 755);

    const auto *kk_1222 = buffer.data(kk + 1222);
    const auto *kk_1223 = buffer.data(kk + 1223);
    const auto *kk_1224 = buffer.data(kk + 1224);
    const auto *kk_1225 = buffer.data(kk + 1225);
    const auto *kk_1226 = buffer.data(kk + 1226);
    const auto *kk_1227 = buffer.data(kk + 1227);
    const auto *kk_1228 = buffer.data(kk + 1228);
    const auto *kk_1229 = buffer.data(kk + 1229);
    const auto *kk_1230 = buffer.data(kk + 1230);
    const auto *kk_1231 = buffer.data(kk + 1231);
    const auto *kk_1232 = buffer.data(kk + 1232);
    const auto *kk_1233 = buffer.data(kk + 1233);
    const auto *kk_1234 = buffer.data(kk + 1234);
    const auto *kk_1235 = buffer.data(kk + 1235);
    const auto *kk_1236 = buffer.data(kk + 1236);
    const auto *kk_1237 = buffer.data(kk + 1237);
    const auto *kk_1238 = buffer.data(kk + 1238);
    const auto *kk_1239 = buffer.data(kk + 1239);
    const auto *kk_1240 = buffer.data(kk + 1240);
    const auto *kk_1241 = buffer.data(kk + 1241);
    const auto *kk_1242 = buffer.data(kk + 1242);
    const auto *kk_1243 = buffer.data(kk + 1243);
    const auto *kk_1244 = buffer.data(kk + 1244);
    const auto *kk_1245 = buffer.data(kk + 1245);
    const auto *kk_1246 = buffer.data(kk + 1246);
    const auto *kk_1247 = buffer.data(kk + 1247);
    const auto *kk_1248 = buffer.data(kk + 1248);
    const auto *kk_1249 = buffer.data(kk + 1249);
    const auto *kk_1250 = buffer.data(kk + 1250);
    const auto *kk_1251 = buffer.data(kk + 1251);
    const auto *kk_1252 = buffer.data(kk + 1252);
    const auto *kk_1253 = buffer.data(kk + 1253);
    const auto *kk_1254 = buffer.data(kk + 1254);
    const auto *kk_1255 = buffer.data(kk + 1255);
    const auto *kk_1256 = buffer.data(kk + 1256);
    const auto *kk_1257 = buffer.data(kk + 1257);
    const auto *kk_1258 = buffer.data(kk + 1258);
    const auto *kk_1259 = buffer.data(kk + 1259);

#pragma omp simd aligned(t_970, t_971, t_972, t_973, t_974, t_975, t_976, hk_754, hk_755, \
                         kk_1222, kk_1223, kk_1224, kk_1225, kk_1226, kk_1227, \
                         kk_1228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = -hk_754[k]
                   + f_0 * kk_1222[k];

        t_971[k] = -hk_755[k]
                   + f_0 * kk_1223[k];

        t_972[k] = f_0 * kk_1224[k];

        t_973[k] = f_0 * kk_1225[k];

        t_974[k] = f_0 * kk_1226[k];

        t_975[k] = f_0 * kk_1227[k];

        t_976[k] = f_0 * kk_1228[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, t_980, t_981, t_982, t_983, t_984, kk_1229, \
                         kk_1230, kk_1231, kk_1232, kk_1233, kk_1234, kk_1235, \
                         kk_1236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_0 * kk_1229[k];

        t_978[k] = f_0 * kk_1230[k];

        t_979[k] = f_0 * kk_1231[k];

        t_980[k] = f_0 * kk_1232[k];

        t_981[k] = f_0 * kk_1233[k];

        t_982[k] = f_0 * kk_1234[k];

        t_983[k] = f_0 * kk_1235[k];

        t_984[k] = f_0 * kk_1236[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, t_989, t_990, t_991, t_992, kk_1237, \
                         kk_1238, kk_1239, kk_1240, kk_1241, kk_1242, kk_1243, \
                         kk_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_0 * kk_1237[k];

        t_986[k] = f_0 * kk_1238[k];

        t_987[k] = f_0 * kk_1239[k];

        t_988[k] = f_0 * kk_1240[k];

        t_989[k] = f_0 * kk_1241[k];

        t_990[k] = f_0 * kk_1242[k];

        t_991[k] = f_0 * kk_1243[k];

        t_992[k] = f_0 * kk_1244[k];
    }

#pragma omp simd aligned(t_993, t_994, t_995, t_996, t_997, t_998, t_999, t_1000, kk_1245, \
                         kk_1246, kk_1247, kk_1248, kk_1249, kk_1250, kk_1251, \
                         kk_1252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_993[k] = f_0 * kk_1245[k];

        t_994[k] = f_0 * kk_1246[k];

        t_995[k] = f_0 * kk_1247[k];

        t_996[k] = f_0 * kk_1248[k];

        t_997[k] = f_0 * kk_1249[k];

        t_998[k] = f_0 * kk_1250[k];

        t_999[k] = f_0 * kk_1251[k];

        t_1000[k] = f_0 * kk_1252[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, t_1004, t_1005, t_1006, t_1007, kk_1253, \
                         kk_1254, kk_1255, kk_1256, kk_1257, kk_1258, \
                         kk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_0 * kk_1253[k];

        t_1002[k] = f_0 * kk_1254[k];

        t_1003[k] = f_0 * kk_1255[k];

        t_1004[k] = f_0 * kk_1256[k];

        t_1005[k] = f_0 * kk_1257[k];

        t_1006[k] = f_0 * kk_1258[k];

        t_1007[k] = f_0 * kk_1259[k];
    }
}

auto
compute_prim_geom_10_ik_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t hk, const size_t kk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ik_electron_repulsion_1_piece0(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_1_piece1(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_1_piece2(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_1_piece3(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_1_piece4(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_1_piece5(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_1_piece6(buffer, target, hk, kk, ncols, alpha);
}

static auto
compute_prim_geom_10_ik_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_1 = buffer.data(hk + 1);
    const auto *hk_2 = buffer.data(hk + 2);
    const auto *hk_3 = buffer.data(hk + 3);
    const auto *hk_4 = buffer.data(hk + 4);
    const auto *hk_5 = buffer.data(hk + 5);
    const auto *hk_6 = buffer.data(hk + 6);
    const auto *hk_7 = buffer.data(hk + 7);
    const auto *hk_8 = buffer.data(hk + 8);
    const auto *hk_9 = buffer.data(hk + 9);
    const auto *hk_10 = buffer.data(hk + 10);
    const auto *hk_11 = buffer.data(hk + 11);
    const auto *hk_12 = buffer.data(hk + 12);
    const auto *hk_13 = buffer.data(hk + 13);
    const auto *hk_14 = buffer.data(hk + 14);
    const auto *hk_15 = buffer.data(hk + 15);
    const auto *hk_16 = buffer.data(hk + 16);
    const auto *hk_17 = buffer.data(hk + 17);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_21 = buffer.data(hk + 21);
    const auto *hk_22 = buffer.data(hk + 22);
    const auto *hk_23 = buffer.data(hk + 23);
    const auto *hk_24 = buffer.data(hk + 24);
    const auto *hk_25 = buffer.data(hk + 25);
    const auto *hk_26 = buffer.data(hk + 26);
    const auto *hk_27 = buffer.data(hk + 27);
    const auto *hk_28 = buffer.data(hk + 28);
    const auto *hk_29 = buffer.data(hk + 29);
    const auto *hk_30 = buffer.data(hk + 30);
    const auto *hk_31 = buffer.data(hk + 31);
    const auto *hk_32 = buffer.data(hk + 32);
    const auto *hk_33 = buffer.data(hk + 33);
    const auto *hk_34 = buffer.data(hk + 34);
    const auto *hk_35 = buffer.data(hk + 35);
    const auto *hk_36 = buffer.data(hk + 36);
    const auto *hk_37 = buffer.data(hk + 37);
    const auto *hk_38 = buffer.data(hk + 38);
    const auto *hk_39 = buffer.data(hk + 39);
    const auto *hk_40 = buffer.data(hk + 40);
    const auto *hk_41 = buffer.data(hk + 41);
    const auto *hk_42 = buffer.data(hk + 42);
    const auto *hk_43 = buffer.data(hk + 43);
    const auto *hk_44 = buffer.data(hk + 44);
    const auto *hk_45 = buffer.data(hk + 45);
    const auto *hk_46 = buffer.data(hk + 46);
    const auto *hk_47 = buffer.data(hk + 47);
    const auto *hk_48 = buffer.data(hk + 48);
    const auto *hk_49 = buffer.data(hk + 49);
    const auto *hk_50 = buffer.data(hk + 50);
    const auto *hk_51 = buffer.data(hk + 51);
    const auto *hk_52 = buffer.data(hk + 52);
    const auto *hk_53 = buffer.data(hk + 53);
    const auto *hk_54 = buffer.data(hk + 54);
    const auto *hk_55 = buffer.data(hk + 55);
    const auto *hk_56 = buffer.data(hk + 56);
    const auto *hk_57 = buffer.data(hk + 57);
    const auto *hk_58 = buffer.data(hk + 58);
    const auto *hk_59 = buffer.data(hk + 59);
    const auto *hk_60 = buffer.data(hk + 60);
    const auto *hk_61 = buffer.data(hk + 61);
    const auto *hk_62 = buffer.data(hk + 62);
    const auto *hk_63 = buffer.data(hk + 63);
    const auto *hk_64 = buffer.data(hk + 64);
    const auto *hk_65 = buffer.data(hk + 65);
    const auto *hk_66 = buffer.data(hk + 66);
    const auto *hk_67 = buffer.data(hk + 67);
    const auto *hk_68 = buffer.data(hk + 68);
    const auto *hk_69 = buffer.data(hk + 69);
    const auto *hk_70 = buffer.data(hk + 70);
    const auto *hk_71 = buffer.data(hk + 71);
    const auto *hk_72 = buffer.data(hk + 72);
    const auto *hk_73 = buffer.data(hk + 73);
    const auto *hk_74 = buffer.data(hk + 74);
    const auto *hk_75 = buffer.data(hk + 75);
    const auto *hk_76 = buffer.data(hk + 76);

    const auto *kk_72 = buffer.data(kk + 72);
    const auto *kk_73 = buffer.data(kk + 73);
    const auto *kk_74 = buffer.data(kk + 74);
    const auto *kk_75 = buffer.data(kk + 75);
    const auto *kk_76 = buffer.data(kk + 76);
    const auto *kk_77 = buffer.data(kk + 77);
    const auto *kk_78 = buffer.data(kk + 78);
    const auto *kk_79 = buffer.data(kk + 79);
    const auto *kk_80 = buffer.data(kk + 80);
    const auto *kk_81 = buffer.data(kk + 81);
    const auto *kk_82 = buffer.data(kk + 82);
    const auto *kk_83 = buffer.data(kk + 83);
    const auto *kk_84 = buffer.data(kk + 84);
    const auto *kk_85 = buffer.data(kk + 85);
    const auto *kk_86 = buffer.data(kk + 86);
    const auto *kk_87 = buffer.data(kk + 87);
    const auto *kk_88 = buffer.data(kk + 88);
    const auto *kk_89 = buffer.data(kk + 89);
    const auto *kk_90 = buffer.data(kk + 90);
    const auto *kk_91 = buffer.data(kk + 91);
    const auto *kk_92 = buffer.data(kk + 92);
    const auto *kk_93 = buffer.data(kk + 93);
    const auto *kk_94 = buffer.data(kk + 94);
    const auto *kk_95 = buffer.data(kk + 95);
    const auto *kk_96 = buffer.data(kk + 96);
    const auto *kk_97 = buffer.data(kk + 97);
    const auto *kk_98 = buffer.data(kk + 98);
    const auto *kk_99 = buffer.data(kk + 99);
    const auto *kk_100 = buffer.data(kk + 100);
    const auto *kk_101 = buffer.data(kk + 101);
    const auto *kk_102 = buffer.data(kk + 102);
    const auto *kk_103 = buffer.data(kk + 103);
    const auto *kk_104 = buffer.data(kk + 104);
    const auto *kk_105 = buffer.data(kk + 105);
    const auto *kk_106 = buffer.data(kk + 106);
    const auto *kk_107 = buffer.data(kk + 107);
    const auto *kk_144 = buffer.data(kk + 144);
    const auto *kk_145 = buffer.data(kk + 145);
    const auto *kk_146 = buffer.data(kk + 146);
    const auto *kk_147 = buffer.data(kk + 147);
    const auto *kk_148 = buffer.data(kk + 148);
    const auto *kk_149 = buffer.data(kk + 149);
    const auto *kk_150 = buffer.data(kk + 150);
    const auto *kk_151 = buffer.data(kk + 151);
    const auto *kk_152 = buffer.data(kk + 152);
    const auto *kk_153 = buffer.data(kk + 153);
    const auto *kk_154 = buffer.data(kk + 154);
    const auto *kk_155 = buffer.data(kk + 155);
    const auto *kk_156 = buffer.data(kk + 156);
    const auto *kk_157 = buffer.data(kk + 157);
    const auto *kk_158 = buffer.data(kk + 158);
    const auto *kk_159 = buffer.data(kk + 159);
    const auto *kk_160 = buffer.data(kk + 160);
    const auto *kk_161 = buffer.data(kk + 161);
    const auto *kk_162 = buffer.data(kk + 162);
    const auto *kk_163 = buffer.data(kk + 163);
    const auto *kk_164 = buffer.data(kk + 164);
    const auto *kk_165 = buffer.data(kk + 165);
    const auto *kk_166 = buffer.data(kk + 166);
    const auto *kk_167 = buffer.data(kk + 167);
    const auto *kk_168 = buffer.data(kk + 168);
    const auto *kk_169 = buffer.data(kk + 169);
    const auto *kk_170 = buffer.data(kk + 170);
    const auto *kk_171 = buffer.data(kk + 171);
    const auto *kk_172 = buffer.data(kk + 172);
    const auto *kk_173 = buffer.data(kk + 173);
    const auto *kk_174 = buffer.data(kk + 174);
    const auto *kk_175 = buffer.data(kk + 175);
    const auto *kk_176 = buffer.data(kk + 176);
    const auto *kk_177 = buffer.data(kk + 177);
    const auto *kk_178 = buffer.data(kk + 178);
    const auto *kk_179 = buffer.data(kk + 179);
    const auto *kk_180 = buffer.data(kk + 180);
    const auto *kk_181 = buffer.data(kk + 181);
    const auto *kk_182 = buffer.data(kk + 182);
    const auto *kk_183 = buffer.data(kk + 183);
    const auto *kk_184 = buffer.data(kk + 184);
    const auto *kk_185 = buffer.data(kk + 185);
    const auto *kk_186 = buffer.data(kk + 186);
    const auto *kk_187 = buffer.data(kk + 187);
    const auto *kk_188 = buffer.data(kk + 188);
    const auto *kk_189 = buffer.data(kk + 189);
    const auto *kk_190 = buffer.data(kk + 190);
    const auto *kk_191 = buffer.data(kk + 191);
    const auto *kk_192 = buffer.data(kk + 192);
    const auto *kk_193 = buffer.data(kk + 193);
    const auto *kk_194 = buffer.data(kk + 194);
    const auto *kk_195 = buffer.data(kk + 195);
    const auto *kk_196 = buffer.data(kk + 196);
    const auto *kk_197 = buffer.data(kk + 197);
    const auto *kk_198 = buffer.data(kk + 198);
    const auto *kk_199 = buffer.data(kk + 199);
    const auto *kk_200 = buffer.data(kk + 200);
    const auto *kk_201 = buffer.data(kk + 201);
    const auto *kk_202 = buffer.data(kk + 202);
    const auto *kk_203 = buffer.data(kk + 203);
    const auto *kk_204 = buffer.data(kk + 204);
    const auto *kk_205 = buffer.data(kk + 205);
    const auto *kk_206 = buffer.data(kk + 206);
    const auto *kk_207 = buffer.data(kk + 207);
    const auto *kk_208 = buffer.data(kk + 208);
    const auto *kk_209 = buffer.data(kk + 209);
    const auto *kk_210 = buffer.data(kk + 210);
    const auto *kk_211 = buffer.data(kk + 211);
    const auto *kk_212 = buffer.data(kk + 212);
    const auto *kk_213 = buffer.data(kk + 213);
    const auto *kk_214 = buffer.data(kk + 214);
    const auto *kk_215 = buffer.data(kk + 215);
    const auto *kk_252 = buffer.data(kk + 252);
    const auto *kk_253 = buffer.data(kk + 253);
    const auto *kk_254 = buffer.data(kk + 254);
    const auto *kk_255 = buffer.data(kk + 255);
    const auto *kk_256 = buffer.data(kk + 256);
    const auto *kk_257 = buffer.data(kk + 257);
    const auto *kk_258 = buffer.data(kk + 258);
    const auto *kk_259 = buffer.data(kk + 259);
    const auto *kk_260 = buffer.data(kk + 260);
    const auto *kk_261 = buffer.data(kk + 261);
    const auto *kk_262 = buffer.data(kk + 262);
    const auto *kk_263 = buffer.data(kk + 263);
    const auto *kk_264 = buffer.data(kk + 264);
    const auto *kk_265 = buffer.data(kk + 265);
    const auto *kk_266 = buffer.data(kk + 266);
    const auto *kk_267 = buffer.data(kk + 267);
    const auto *kk_268 = buffer.data(kk + 268);
    const auto *kk_269 = buffer.data(kk + 269);
    const auto *kk_270 = buffer.data(kk + 270);
    const auto *kk_271 = buffer.data(kk + 271);
    const auto *kk_272 = buffer.data(kk + 272);
    const auto *kk_273 = buffer.data(kk + 273);
    const auto *kk_274 = buffer.data(kk + 274);
    const auto *kk_275 = buffer.data(kk + 275);
    const auto *kk_276 = buffer.data(kk + 276);
    const auto *kk_277 = buffer.data(kk + 277);
    const auto *kk_278 = buffer.data(kk + 278);
    const auto *kk_279 = buffer.data(kk + 279);
    const auto *kk_280 = buffer.data(kk + 280);
    const auto *kk_281 = buffer.data(kk + 281);
    const auto *kk_282 = buffer.data(kk + 282);
    const auto *kk_283 = buffer.data(kk + 283);
    const auto *kk_284 = buffer.data(kk + 284);
    const auto *kk_285 = buffer.data(kk + 285);
    const auto *kk_286 = buffer.data(kk + 286);
    const auto *kk_287 = buffer.data(kk + 287);
    const auto *kk_288 = buffer.data(kk + 288);
    const auto *kk_289 = buffer.data(kk + 289);
    const auto *kk_290 = buffer.data(kk + 290);
    const auto *kk_291 = buffer.data(kk + 291);
    const auto *kk_292 = buffer.data(kk + 292);
    const auto *kk_293 = buffer.data(kk + 293);
    const auto *kk_294 = buffer.data(kk + 294);
    const auto *kk_295 = buffer.data(kk + 295);
    const auto *kk_296 = buffer.data(kk + 296);
    const auto *kk_297 = buffer.data(kk + 297);
    const auto *kk_298 = buffer.data(kk + 298);
    const auto *kk_299 = buffer.data(kk + 299);
    const auto *kk_300 = buffer.data(kk + 300);
    const auto *kk_301 = buffer.data(kk + 301);
    const auto *kk_302 = buffer.data(kk + 302);
    const auto *kk_303 = buffer.data(kk + 303);
    const auto *kk_304 = buffer.data(kk + 304);
    const auto *kk_305 = buffer.data(kk + 305);
    const auto *kk_306 = buffer.data(kk + 306);
    const auto *kk_307 = buffer.data(kk + 307);
    const auto *kk_308 = buffer.data(kk + 308);
    const auto *kk_309 = buffer.data(kk + 309);
    const auto *kk_310 = buffer.data(kk + 310);
    const auto *kk_311 = buffer.data(kk + 311);
    const auto *kk_312 = buffer.data(kk + 312);
    const auto *kk_313 = buffer.data(kk + 313);
    const auto *kk_314 = buffer.data(kk + 314);
    const auto *kk_315 = buffer.data(kk + 315);
    const auto *kk_316 = buffer.data(kk + 316);
    const auto *kk_317 = buffer.data(kk + 317);
    const auto *kk_318 = buffer.data(kk + 318);
    const auto *kk_319 = buffer.data(kk + 319);
    const auto *kk_320 = buffer.data(kk + 320);
    const auto *kk_321 = buffer.data(kk + 321);
    const auto *kk_322 = buffer.data(kk + 322);
    const auto *kk_323 = buffer.data(kk + 323);
    const auto *kk_324 = buffer.data(kk + 324);
    const auto *kk_325 = buffer.data(kk + 325);
    const auto *kk_326 = buffer.data(kk + 326);
    const auto *kk_327 = buffer.data(kk + 327);
    const auto *kk_328 = buffer.data(kk + 328);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, kk_72, kk_73, kk_74, kk_75, \
                         kk_76, kk_77, kk_78, kk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kk_72[k];

        t_1[k] = f_0 * kk_73[k];

        t_2[k] = f_0 * kk_74[k];

        t_3[k] = f_0 * kk_75[k];

        t_4[k] = f_0 * kk_76[k];

        t_5[k] = f_0 * kk_77[k];

        t_6[k] = f_0 * kk_78[k];

        t_7[k] = f_0 * kk_79[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, kk_80, kk_81, kk_82, \
                         kk_83, kk_84, kk_85, kk_86, kk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * kk_80[k];

        t_9[k] = f_0 * kk_81[k];

        t_10[k] = f_0 * kk_82[k];

        t_11[k] = f_0 * kk_83[k];

        t_12[k] = f_0 * kk_84[k];

        t_13[k] = f_0 * kk_85[k];

        t_14[k] = f_0 * kk_86[k];

        t_15[k] = f_0 * kk_87[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, kk_88, kk_89, kk_90, \
                         kk_91, kk_92, kk_93, kk_94, kk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * kk_88[k];

        t_17[k] = f_0 * kk_89[k];

        t_18[k] = f_0 * kk_90[k];

        t_19[k] = f_0 * kk_91[k];

        t_20[k] = f_0 * kk_92[k];

        t_21[k] = f_0 * kk_93[k];

        t_22[k] = f_0 * kk_94[k];

        t_23[k] = f_0 * kk_95[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, kk_96, kk_97, kk_98, \
                         kk_99, kk_100, kk_101, kk_102, kk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * kk_96[k];

        t_25[k] = f_0 * kk_97[k];

        t_26[k] = f_0 * kk_98[k];

        t_27[k] = f_0 * kk_99[k];

        t_28[k] = f_0 * kk_100[k];

        t_29[k] = f_0 * kk_101[k];

        t_30[k] = f_0 * kk_102[k];

        t_31[k] = f_0 * kk_103[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, kk_104, kk_105, \
                         kk_106, kk_107, kk_144, kk_145, kk_146, \
                         kk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * kk_104[k];

        t_33[k] = f_0 * kk_105[k];

        t_34[k] = f_0 * kk_106[k];

        t_35[k] = f_0 * kk_107[k];

        t_36[k] = f_0 * kk_144[k];

        t_37[k] = f_0 * kk_145[k];

        t_38[k] = f_0 * kk_146[k];

        t_39[k] = f_0 * kk_147[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, kk_148, kk_149, \
                         kk_150, kk_151, kk_152, kk_153, kk_154, \
                         kk_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * kk_148[k];

        t_41[k] = f_0 * kk_149[k];

        t_42[k] = f_0 * kk_150[k];

        t_43[k] = f_0 * kk_151[k];

        t_44[k] = f_0 * kk_152[k];

        t_45[k] = f_0 * kk_153[k];

        t_46[k] = f_0 * kk_154[k];

        t_47[k] = f_0 * kk_155[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, kk_156, kk_157, \
                         kk_158, kk_159, kk_160, kk_161, kk_162, \
                         kk_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * kk_156[k];

        t_49[k] = f_0 * kk_157[k];

        t_50[k] = f_0 * kk_158[k];

        t_51[k] = f_0 * kk_159[k];

        t_52[k] = f_0 * kk_160[k];

        t_53[k] = f_0 * kk_161[k];

        t_54[k] = f_0 * kk_162[k];

        t_55[k] = f_0 * kk_163[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, kk_164, kk_165, \
                         kk_166, kk_167, kk_168, kk_169, kk_170, \
                         kk_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * kk_164[k];

        t_57[k] = f_0 * kk_165[k];

        t_58[k] = f_0 * kk_166[k];

        t_59[k] = f_0 * kk_167[k];

        t_60[k] = f_0 * kk_168[k];

        t_61[k] = f_0 * kk_169[k];

        t_62[k] = f_0 * kk_170[k];

        t_63[k] = f_0 * kk_171[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, kk_172, kk_173, \
                         kk_174, kk_175, kk_176, kk_177, kk_178, \
                         kk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * kk_172[k];

        t_65[k] = f_0 * kk_173[k];

        t_66[k] = f_0 * kk_174[k];

        t_67[k] = f_0 * kk_175[k];

        t_68[k] = f_0 * kk_176[k];

        t_69[k] = f_0 * kk_177[k];

        t_70[k] = f_0 * kk_178[k];

        t_71[k] = f_0 * kk_179[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, hk_0, hk_1, hk_2, hk_3, hk_4, kk_180, \
                         kk_181, kk_182, kk_183, kk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -hk_0[k]
                  + f_0 * kk_180[k];

        t_73[k] = -hk_1[k]
                  + f_0 * kk_181[k];

        t_74[k] = -hk_2[k]
                  + f_0 * kk_182[k];

        t_75[k] = -hk_3[k]
                  + f_0 * kk_183[k];

        t_76[k] = -hk_4[k]
                  + f_0 * kk_184[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, hk_5, hk_6, hk_7, hk_8, hk_9, kk_185, \
                         kk_186, kk_187, kk_188, kk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -hk_5[k]
                  + f_0 * kk_185[k];

        t_78[k] = -hk_6[k]
                  + f_0 * kk_186[k];

        t_79[k] = -hk_7[k]
                  + f_0 * kk_187[k];

        t_80[k] = -hk_8[k]
                  + f_0 * kk_188[k];

        t_81[k] = -hk_9[k]
                  + f_0 * kk_189[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, hk_10, hk_11, hk_12, hk_13, hk_14, \
                         kk_190, kk_191, kk_192, kk_193, kk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -hk_10[k]
                  + f_0 * kk_190[k];

        t_83[k] = -hk_11[k]
                  + f_0 * kk_191[k];

        t_84[k] = -hk_12[k]
                  + f_0 * kk_192[k];

        t_85[k] = -hk_13[k]
                  + f_0 * kk_193[k];

        t_86[k] = -hk_14[k]
                  + f_0 * kk_194[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, hk_15, hk_16, hk_17, hk_18, hk_19, \
                         kk_195, kk_196, kk_197, kk_198, kk_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -hk_15[k]
                  + f_0 * kk_195[k];

        t_88[k] = -hk_16[k]
                  + f_0 * kk_196[k];

        t_89[k] = -hk_17[k]
                  + f_0 * kk_197[k];

        t_90[k] = -hk_18[k]
                  + f_0 * kk_198[k];

        t_91[k] = -hk_19[k]
                  + f_0 * kk_199[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, hk_20, hk_21, hk_22, hk_23, hk_24, \
                         kk_200, kk_201, kk_202, kk_203, kk_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = -hk_20[k]
                  + f_0 * kk_200[k];

        t_93[k] = -hk_21[k]
                  + f_0 * kk_201[k];

        t_94[k] = -hk_22[k]
                  + f_0 * kk_202[k];

        t_95[k] = -hk_23[k]
                  + f_0 * kk_203[k];

        t_96[k] = -hk_24[k]
                  + f_0 * kk_204[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, hk_25, hk_26, hk_27, hk_28, hk_29, \
                         kk_205, kk_206, kk_207, kk_208, kk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -hk_25[k]
                  + f_0 * kk_205[k];

        t_98[k] = -hk_26[k]
                  + f_0 * kk_206[k];

        t_99[k] = -hk_27[k]
                  + f_0 * kk_207[k];

        t_100[k] = -hk_28[k]
                   + f_0 * kk_208[k];

        t_101[k] = -hk_29[k]
                   + f_0 * kk_209[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, hk_30, hk_31, hk_32, hk_33, hk_34, \
                         kk_210, kk_211, kk_212, kk_213, kk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -hk_30[k]
                   + f_0 * kk_210[k];

        t_103[k] = -hk_31[k]
                   + f_0 * kk_211[k];

        t_104[k] = -hk_32[k]
                   + f_0 * kk_212[k];

        t_105[k] = -hk_33[k]
                   + f_0 * kk_213[k];

        t_106[k] = -hk_34[k]
                   + f_0 * kk_214[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, hk_35, kk_215, \
                         kk_252, kk_253, kk_254, kk_255, kk_256, \
                         kk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -hk_35[k]
                   + f_0 * kk_215[k];

        t_108[k] = f_0 * kk_252[k];

        t_109[k] = f_0 * kk_253[k];

        t_110[k] = f_0 * kk_254[k];

        t_111[k] = f_0 * kk_255[k];

        t_112[k] = f_0 * kk_256[k];

        t_113[k] = f_0 * kk_257[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, t_120, t_121, kk_258, \
                         kk_259, kk_260, kk_261, kk_262, kk_263, kk_264, \
                         kk_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * kk_258[k];

        t_115[k] = f_0 * kk_259[k];

        t_116[k] = f_0 * kk_260[k];

        t_117[k] = f_0 * kk_261[k];

        t_118[k] = f_0 * kk_262[k];

        t_119[k] = f_0 * kk_263[k];

        t_120[k] = f_0 * kk_264[k];

        t_121[k] = f_0 * kk_265[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, t_128, t_129, kk_266, \
                         kk_267, kk_268, kk_269, kk_270, kk_271, kk_272, \
                         kk_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * kk_266[k];

        t_123[k] = f_0 * kk_267[k];

        t_124[k] = f_0 * kk_268[k];

        t_125[k] = f_0 * kk_269[k];

        t_126[k] = f_0 * kk_270[k];

        t_127[k] = f_0 * kk_271[k];

        t_128[k] = f_0 * kk_272[k];

        t_129[k] = f_0 * kk_273[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, t_135, t_136, t_137, kk_274, \
                         kk_275, kk_276, kk_277, kk_278, kk_279, kk_280, \
                         kk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * kk_274[k];

        t_131[k] = f_0 * kk_275[k];

        t_132[k] = f_0 * kk_276[k];

        t_133[k] = f_0 * kk_277[k];

        t_134[k] = f_0 * kk_278[k];

        t_135[k] = f_0 * kk_279[k];

        t_136[k] = f_0 * kk_280[k];

        t_137[k] = f_0 * kk_281[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, t_144, hk_36, kk_282, \
                         kk_283, kk_284, kk_285, kk_286, kk_287, \
                         kk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * kk_282[k];

        t_139[k] = f_0 * kk_283[k];

        t_140[k] = f_0 * kk_284[k];

        t_141[k] = f_0 * kk_285[k];

        t_142[k] = f_0 * kk_286[k];

        t_143[k] = f_0 * kk_287[k];

        t_144[k] = -hk_36[k]
                   + f_0 * kk_288[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, hk_37, hk_38, hk_39, hk_40, hk_41, \
                         kk_289, kk_290, kk_291, kk_292, kk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -hk_37[k]
                   + f_0 * kk_289[k];

        t_146[k] = -hk_38[k]
                   + f_0 * kk_290[k];

        t_147[k] = -hk_39[k]
                   + f_0 * kk_291[k];

        t_148[k] = -hk_40[k]
                   + f_0 * kk_292[k];

        t_149[k] = -hk_41[k]
                   + f_0 * kk_293[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, hk_42, hk_43, hk_44, hk_45, hk_46, \
                         kk_294, kk_295, kk_296, kk_297, kk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -hk_42[k]
                   + f_0 * kk_294[k];

        t_151[k] = -hk_43[k]
                   + f_0 * kk_295[k];

        t_152[k] = -hk_44[k]
                   + f_0 * kk_296[k];

        t_153[k] = -hk_45[k]
                   + f_0 * kk_297[k];

        t_154[k] = -hk_46[k]
                   + f_0 * kk_298[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, hk_47, hk_48, hk_49, hk_50, hk_51, \
                         kk_299, kk_300, kk_301, kk_302, kk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -hk_47[k]
                   + f_0 * kk_299[k];

        t_156[k] = -hk_48[k]
                   + f_0 * kk_300[k];

        t_157[k] = -hk_49[k]
                   + f_0 * kk_301[k];

        t_158[k] = -hk_50[k]
                   + f_0 * kk_302[k];

        t_159[k] = -hk_51[k]
                   + f_0 * kk_303[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, hk_52, hk_53, hk_54, hk_55, hk_56, \
                         kk_304, kk_305, kk_306, kk_307, kk_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -hk_52[k]
                   + f_0 * kk_304[k];

        t_161[k] = -hk_53[k]
                   + f_0 * kk_305[k];

        t_162[k] = -hk_54[k]
                   + f_0 * kk_306[k];

        t_163[k] = -hk_55[k]
                   + f_0 * kk_307[k];

        t_164[k] = -hk_56[k]
                   + f_0 * kk_308[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, hk_57, hk_58, hk_59, hk_60, hk_61, \
                         kk_309, kk_310, kk_311, kk_312, kk_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -hk_57[k]
                   + f_0 * kk_309[k];

        t_166[k] = -hk_58[k]
                   + f_0 * kk_310[k];

        t_167[k] = -hk_59[k]
                   + f_0 * kk_311[k];

        t_168[k] = -hk_60[k]
                   + f_0 * kk_312[k];

        t_169[k] = -hk_61[k]
                   + f_0 * kk_313[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, hk_62, hk_63, hk_64, hk_65, hk_66, \
                         kk_314, kk_315, kk_316, kk_317, kk_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -hk_62[k]
                   + f_0 * kk_314[k];

        t_171[k] = -hk_63[k]
                   + f_0 * kk_315[k];

        t_172[k] = -hk_64[k]
                   + f_0 * kk_316[k];

        t_173[k] = -hk_65[k]
                   + f_0 * kk_317[k];

        t_174[k] = -hk_66[k]
                   + f_0 * kk_318[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, hk_67, hk_68, hk_69, hk_70, hk_71, \
                         kk_319, kk_320, kk_321, kk_322, kk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -hk_67[k]
                   + f_0 * kk_319[k];

        t_176[k] = -hk_68[k]
                   + f_0 * kk_320[k];

        t_177[k] = -hk_69[k]
                   + f_0 * kk_321[k];

        t_178[k] = -hk_70[k]
                   + f_0 * kk_322[k];

        t_179[k] = -hk_71[k]
                   + f_0 * kk_323[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, hk_72, hk_73, hk_74, hk_75, hk_76, \
                         kk_324, kk_325, kk_326, kk_327, kk_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * hk_72[k]
                   + f_0 * kk_324[k];

        t_181[k] = -2.0 * hk_73[k]
                   + f_0 * kk_325[k];

        t_182[k] = -2.0 * hk_74[k]
                   + f_0 * kk_326[k];

        t_183[k] = -2.0 * hk_75[k]
                   + f_0 * kk_327[k];

        t_184[k] = -2.0 * hk_76[k]
                   + f_0 * kk_328[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_77 = buffer.data(hk + 77);
    const auto *hk_78 = buffer.data(hk + 78);
    const auto *hk_79 = buffer.data(hk + 79);
    const auto *hk_80 = buffer.data(hk + 80);
    const auto *hk_81 = buffer.data(hk + 81);
    const auto *hk_82 = buffer.data(hk + 82);
    const auto *hk_83 = buffer.data(hk + 83);
    const auto *hk_84 = buffer.data(hk + 84);
    const auto *hk_85 = buffer.data(hk + 85);
    const auto *hk_86 = buffer.data(hk + 86);
    const auto *hk_87 = buffer.data(hk + 87);
    const auto *hk_88 = buffer.data(hk + 88);
    const auto *hk_89 = buffer.data(hk + 89);
    const auto *hk_90 = buffer.data(hk + 90);
    const auto *hk_91 = buffer.data(hk + 91);
    const auto *hk_92 = buffer.data(hk + 92);
    const auto *hk_93 = buffer.data(hk + 93);
    const auto *hk_94 = buffer.data(hk + 94);
    const auto *hk_95 = buffer.data(hk + 95);
    const auto *hk_96 = buffer.data(hk + 96);
    const auto *hk_97 = buffer.data(hk + 97);
    const auto *hk_98 = buffer.data(hk + 98);
    const auto *hk_99 = buffer.data(hk + 99);
    const auto *hk_100 = buffer.data(hk + 100);
    const auto *hk_101 = buffer.data(hk + 101);
    const auto *hk_102 = buffer.data(hk + 102);
    const auto *hk_103 = buffer.data(hk + 103);
    const auto *hk_104 = buffer.data(hk + 104);
    const auto *hk_105 = buffer.data(hk + 105);
    const auto *hk_106 = buffer.data(hk + 106);
    const auto *hk_107 = buffer.data(hk + 107);
    const auto *hk_108 = buffer.data(hk + 108);
    const auto *hk_109 = buffer.data(hk + 109);
    const auto *hk_110 = buffer.data(hk + 110);
    const auto *hk_111 = buffer.data(hk + 111);
    const auto *hk_112 = buffer.data(hk + 112);
    const auto *hk_113 = buffer.data(hk + 113);
    const auto *hk_114 = buffer.data(hk + 114);
    const auto *hk_115 = buffer.data(hk + 115);
    const auto *hk_116 = buffer.data(hk + 116);
    const auto *hk_117 = buffer.data(hk + 117);
    const auto *hk_118 = buffer.data(hk + 118);
    const auto *hk_119 = buffer.data(hk + 119);
    const auto *hk_120 = buffer.data(hk + 120);
    const auto *hk_121 = buffer.data(hk + 121);
    const auto *hk_122 = buffer.data(hk + 122);
    const auto *hk_123 = buffer.data(hk + 123);
    const auto *hk_124 = buffer.data(hk + 124);
    const auto *hk_125 = buffer.data(hk + 125);
    const auto *hk_126 = buffer.data(hk + 126);
    const auto *hk_127 = buffer.data(hk + 127);
    const auto *hk_128 = buffer.data(hk + 128);
    const auto *hk_129 = buffer.data(hk + 129);
    const auto *hk_130 = buffer.data(hk + 130);
    const auto *hk_131 = buffer.data(hk + 131);
    const auto *hk_132 = buffer.data(hk + 132);
    const auto *hk_133 = buffer.data(hk + 133);
    const auto *hk_134 = buffer.data(hk + 134);
    const auto *hk_135 = buffer.data(hk + 135);
    const auto *hk_136 = buffer.data(hk + 136);
    const auto *hk_137 = buffer.data(hk + 137);
    const auto *hk_138 = buffer.data(hk + 138);
    const auto *hk_139 = buffer.data(hk + 139);
    const auto *hk_140 = buffer.data(hk + 140);
    const auto *hk_141 = buffer.data(hk + 141);
    const auto *hk_142 = buffer.data(hk + 142);
    const auto *hk_143 = buffer.data(hk + 143);
    const auto *hk_144 = buffer.data(hk + 144);
    const auto *hk_145 = buffer.data(hk + 145);
    const auto *hk_146 = buffer.data(hk + 146);
    const auto *hk_147 = buffer.data(hk + 147);
    const auto *hk_148 = buffer.data(hk + 148);
    const auto *hk_149 = buffer.data(hk + 149);
    const auto *hk_150 = buffer.data(hk + 150);
    const auto *hk_151 = buffer.data(hk + 151);
    const auto *hk_152 = buffer.data(hk + 152);
    const auto *hk_153 = buffer.data(hk + 153);
    const auto *hk_154 = buffer.data(hk + 154);
    const auto *hk_155 = buffer.data(hk + 155);
    const auto *hk_156 = buffer.data(hk + 156);
    const auto *hk_157 = buffer.data(hk + 157);
    const auto *hk_158 = buffer.data(hk + 158);
    const auto *hk_159 = buffer.data(hk + 159);
    const auto *hk_160 = buffer.data(hk + 160);
    const auto *hk_161 = buffer.data(hk + 161);
    const auto *hk_162 = buffer.data(hk + 162);
    const auto *hk_163 = buffer.data(hk + 163);
    const auto *hk_164 = buffer.data(hk + 164);
    const auto *hk_165 = buffer.data(hk + 165);
    const auto *hk_166 = buffer.data(hk + 166);
    const auto *hk_167 = buffer.data(hk + 167);
    const auto *hk_168 = buffer.data(hk + 168);
    const auto *hk_169 = buffer.data(hk + 169);
    const auto *hk_170 = buffer.data(hk + 170);
    const auto *hk_171 = buffer.data(hk + 171);
    const auto *hk_172 = buffer.data(hk + 172);
    const auto *hk_173 = buffer.data(hk + 173);
    const auto *hk_174 = buffer.data(hk + 174);
    const auto *hk_175 = buffer.data(hk + 175);
    const auto *hk_176 = buffer.data(hk + 176);
    const auto *hk_177 = buffer.data(hk + 177);
    const auto *hk_178 = buffer.data(hk + 178);
    const auto *hk_179 = buffer.data(hk + 179);
    const auto *hk_180 = buffer.data(hk + 180);
    const auto *hk_181 = buffer.data(hk + 181);
    const auto *hk_182 = buffer.data(hk + 182);
    const auto *hk_183 = buffer.data(hk + 183);
    const auto *hk_184 = buffer.data(hk + 184);
    const auto *hk_185 = buffer.data(hk + 185);
    const auto *hk_186 = buffer.data(hk + 186);
    const auto *hk_187 = buffer.data(hk + 187);
    const auto *hk_188 = buffer.data(hk + 188);
    const auto *hk_189 = buffer.data(hk + 189);
    const auto *hk_190 = buffer.data(hk + 190);
    const auto *hk_191 = buffer.data(hk + 191);
    const auto *hk_192 = buffer.data(hk + 192);
    const auto *hk_193 = buffer.data(hk + 193);
    const auto *hk_194 = buffer.data(hk + 194);
    const auto *hk_195 = buffer.data(hk + 195);
    const auto *hk_196 = buffer.data(hk + 196);
    const auto *hk_197 = buffer.data(hk + 197);
    const auto *hk_198 = buffer.data(hk + 198);
    const auto *hk_199 = buffer.data(hk + 199);
    const auto *hk_200 = buffer.data(hk + 200);
    const auto *hk_201 = buffer.data(hk + 201);
    const auto *hk_202 = buffer.data(hk + 202);
    const auto *hk_203 = buffer.data(hk + 203);

    const auto *kk_329 = buffer.data(kk + 329);
    const auto *kk_330 = buffer.data(kk + 330);
    const auto *kk_331 = buffer.data(kk + 331);
    const auto *kk_332 = buffer.data(kk + 332);
    const auto *kk_333 = buffer.data(kk + 333);
    const auto *kk_334 = buffer.data(kk + 334);
    const auto *kk_335 = buffer.data(kk + 335);
    const auto *kk_336 = buffer.data(kk + 336);
    const auto *kk_337 = buffer.data(kk + 337);
    const auto *kk_338 = buffer.data(kk + 338);
    const auto *kk_339 = buffer.data(kk + 339);
    const auto *kk_340 = buffer.data(kk + 340);
    const auto *kk_341 = buffer.data(kk + 341);
    const auto *kk_342 = buffer.data(kk + 342);
    const auto *kk_343 = buffer.data(kk + 343);
    const auto *kk_344 = buffer.data(kk + 344);
    const auto *kk_345 = buffer.data(kk + 345);
    const auto *kk_346 = buffer.data(kk + 346);
    const auto *kk_347 = buffer.data(kk + 347);
    const auto *kk_348 = buffer.data(kk + 348);
    const auto *kk_349 = buffer.data(kk + 349);
    const auto *kk_350 = buffer.data(kk + 350);
    const auto *kk_351 = buffer.data(kk + 351);
    const auto *kk_352 = buffer.data(kk + 352);
    const auto *kk_353 = buffer.data(kk + 353);
    const auto *kk_354 = buffer.data(kk + 354);
    const auto *kk_355 = buffer.data(kk + 355);
    const auto *kk_356 = buffer.data(kk + 356);
    const auto *kk_357 = buffer.data(kk + 357);
    const auto *kk_358 = buffer.data(kk + 358);
    const auto *kk_359 = buffer.data(kk + 359);
    const auto *kk_396 = buffer.data(kk + 396);
    const auto *kk_397 = buffer.data(kk + 397);
    const auto *kk_398 = buffer.data(kk + 398);
    const auto *kk_399 = buffer.data(kk + 399);
    const auto *kk_400 = buffer.data(kk + 400);
    const auto *kk_401 = buffer.data(kk + 401);
    const auto *kk_402 = buffer.data(kk + 402);
    const auto *kk_403 = buffer.data(kk + 403);
    const auto *kk_404 = buffer.data(kk + 404);
    const auto *kk_405 = buffer.data(kk + 405);
    const auto *kk_406 = buffer.data(kk + 406);
    const auto *kk_407 = buffer.data(kk + 407);
    const auto *kk_408 = buffer.data(kk + 408);
    const auto *kk_409 = buffer.data(kk + 409);
    const auto *kk_410 = buffer.data(kk + 410);
    const auto *kk_411 = buffer.data(kk + 411);
    const auto *kk_412 = buffer.data(kk + 412);
    const auto *kk_413 = buffer.data(kk + 413);
    const auto *kk_414 = buffer.data(kk + 414);
    const auto *kk_415 = buffer.data(kk + 415);
    const auto *kk_416 = buffer.data(kk + 416);
    const auto *kk_417 = buffer.data(kk + 417);
    const auto *kk_418 = buffer.data(kk + 418);
    const auto *kk_419 = buffer.data(kk + 419);
    const auto *kk_420 = buffer.data(kk + 420);
    const auto *kk_421 = buffer.data(kk + 421);
    const auto *kk_422 = buffer.data(kk + 422);
    const auto *kk_423 = buffer.data(kk + 423);
    const auto *kk_424 = buffer.data(kk + 424);
    const auto *kk_425 = buffer.data(kk + 425);
    const auto *kk_426 = buffer.data(kk + 426);
    const auto *kk_427 = buffer.data(kk + 427);
    const auto *kk_428 = buffer.data(kk + 428);
    const auto *kk_429 = buffer.data(kk + 429);
    const auto *kk_430 = buffer.data(kk + 430);
    const auto *kk_431 = buffer.data(kk + 431);
    const auto *kk_432 = buffer.data(kk + 432);
    const auto *kk_433 = buffer.data(kk + 433);
    const auto *kk_434 = buffer.data(kk + 434);
    const auto *kk_435 = buffer.data(kk + 435);
    const auto *kk_436 = buffer.data(kk + 436);
    const auto *kk_437 = buffer.data(kk + 437);
    const auto *kk_438 = buffer.data(kk + 438);
    const auto *kk_439 = buffer.data(kk + 439);
    const auto *kk_440 = buffer.data(kk + 440);
    const auto *kk_441 = buffer.data(kk + 441);
    const auto *kk_442 = buffer.data(kk + 442);
    const auto *kk_443 = buffer.data(kk + 443);
    const auto *kk_444 = buffer.data(kk + 444);
    const auto *kk_445 = buffer.data(kk + 445);
    const auto *kk_446 = buffer.data(kk + 446);
    const auto *kk_447 = buffer.data(kk + 447);
    const auto *kk_448 = buffer.data(kk + 448);
    const auto *kk_449 = buffer.data(kk + 449);
    const auto *kk_450 = buffer.data(kk + 450);
    const auto *kk_451 = buffer.data(kk + 451);
    const auto *kk_452 = buffer.data(kk + 452);
    const auto *kk_453 = buffer.data(kk + 453);
    const auto *kk_454 = buffer.data(kk + 454);
    const auto *kk_455 = buffer.data(kk + 455);
    const auto *kk_456 = buffer.data(kk + 456);
    const auto *kk_457 = buffer.data(kk + 457);
    const auto *kk_458 = buffer.data(kk + 458);
    const auto *kk_459 = buffer.data(kk + 459);
    const auto *kk_460 = buffer.data(kk + 460);
    const auto *kk_461 = buffer.data(kk + 461);
    const auto *kk_462 = buffer.data(kk + 462);
    const auto *kk_463 = buffer.data(kk + 463);
    const auto *kk_464 = buffer.data(kk + 464);
    const auto *kk_465 = buffer.data(kk + 465);
    const auto *kk_466 = buffer.data(kk + 466);
    const auto *kk_467 = buffer.data(kk + 467);
    const auto *kk_468 = buffer.data(kk + 468);
    const auto *kk_469 = buffer.data(kk + 469);
    const auto *kk_470 = buffer.data(kk + 470);
    const auto *kk_471 = buffer.data(kk + 471);
    const auto *kk_472 = buffer.data(kk + 472);
    const auto *kk_473 = buffer.data(kk + 473);
    const auto *kk_474 = buffer.data(kk + 474);
    const auto *kk_475 = buffer.data(kk + 475);
    const auto *kk_476 = buffer.data(kk + 476);
    const auto *kk_477 = buffer.data(kk + 477);
    const auto *kk_478 = buffer.data(kk + 478);
    const auto *kk_479 = buffer.data(kk + 479);
    const auto *kk_480 = buffer.data(kk + 480);
    const auto *kk_481 = buffer.data(kk + 481);
    const auto *kk_482 = buffer.data(kk + 482);
    const auto *kk_483 = buffer.data(kk + 483);
    const auto *kk_484 = buffer.data(kk + 484);
    const auto *kk_485 = buffer.data(kk + 485);
    const auto *kk_486 = buffer.data(kk + 486);
    const auto *kk_487 = buffer.data(kk + 487);
    const auto *kk_488 = buffer.data(kk + 488);
    const auto *kk_489 = buffer.data(kk + 489);
    const auto *kk_490 = buffer.data(kk + 490);
    const auto *kk_491 = buffer.data(kk + 491);
    const auto *kk_492 = buffer.data(kk + 492);
    const auto *kk_493 = buffer.data(kk + 493);
    const auto *kk_494 = buffer.data(kk + 494);
    const auto *kk_495 = buffer.data(kk + 495);
    const auto *kk_496 = buffer.data(kk + 496);
    const auto *kk_497 = buffer.data(kk + 497);
    const auto *kk_498 = buffer.data(kk + 498);
    const auto *kk_499 = buffer.data(kk + 499);
    const auto *kk_500 = buffer.data(kk + 500);
    const auto *kk_501 = buffer.data(kk + 501);
    const auto *kk_502 = buffer.data(kk + 502);
    const auto *kk_503 = buffer.data(kk + 503);
    const auto *kk_504 = buffer.data(kk + 504);
    const auto *kk_505 = buffer.data(kk + 505);
    const auto *kk_506 = buffer.data(kk + 506);
    const auto *kk_507 = buffer.data(kk + 507);
    const auto *kk_508 = buffer.data(kk + 508);
    const auto *kk_509 = buffer.data(kk + 509);
    const auto *kk_510 = buffer.data(kk + 510);
    const auto *kk_511 = buffer.data(kk + 511);
    const auto *kk_512 = buffer.data(kk + 512);
    const auto *kk_513 = buffer.data(kk + 513);
    const auto *kk_514 = buffer.data(kk + 514);
    const auto *kk_515 = buffer.data(kk + 515);
    const auto *kk_516 = buffer.data(kk + 516);
    const auto *kk_517 = buffer.data(kk + 517);
    const auto *kk_518 = buffer.data(kk + 518);
    const auto *kk_519 = buffer.data(kk + 519);
    const auto *kk_520 = buffer.data(kk + 520);
    const auto *kk_521 = buffer.data(kk + 521);
    const auto *kk_522 = buffer.data(kk + 522);
    const auto *kk_523 = buffer.data(kk + 523);
    const auto *kk_524 = buffer.data(kk + 524);
    const auto *kk_525 = buffer.data(kk + 525);
    const auto *kk_526 = buffer.data(kk + 526);
    const auto *kk_527 = buffer.data(kk + 527);

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, hk_77, hk_78, hk_79, hk_80, hk_81, \
                         kk_329, kk_330, kk_331, kk_332, kk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * hk_77[k]
                   + f_0 * kk_329[k];

        t_186[k] = -2.0 * hk_78[k]
                   + f_0 * kk_330[k];

        t_187[k] = -2.0 * hk_79[k]
                   + f_0 * kk_331[k];

        t_188[k] = -2.0 * hk_80[k]
                   + f_0 * kk_332[k];

        t_189[k] = -2.0 * hk_81[k]
                   + f_0 * kk_333[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, hk_82, hk_83, hk_84, hk_85, hk_86, \
                         kk_334, kk_335, kk_336, kk_337, kk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * hk_82[k]
                   + f_0 * kk_334[k];

        t_191[k] = -2.0 * hk_83[k]
                   + f_0 * kk_335[k];

        t_192[k] = -2.0 * hk_84[k]
                   + f_0 * kk_336[k];

        t_193[k] = -2.0 * hk_85[k]
                   + f_0 * kk_337[k];

        t_194[k] = -2.0 * hk_86[k]
                   + f_0 * kk_338[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, hk_87, hk_88, hk_89, hk_90, hk_91, \
                         kk_339, kk_340, kk_341, kk_342, kk_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * hk_87[k]
                   + f_0 * kk_339[k];

        t_196[k] = -2.0 * hk_88[k]
                   + f_0 * kk_340[k];

        t_197[k] = -2.0 * hk_89[k]
                   + f_0 * kk_341[k];

        t_198[k] = -2.0 * hk_90[k]
                   + f_0 * kk_342[k];

        t_199[k] = -2.0 * hk_91[k]
                   + f_0 * kk_343[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, hk_92, hk_93, hk_94, hk_95, hk_96, \
                         kk_344, kk_345, kk_346, kk_347, kk_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * hk_92[k]
                   + f_0 * kk_344[k];

        t_201[k] = -2.0 * hk_93[k]
                   + f_0 * kk_345[k];

        t_202[k] = -2.0 * hk_94[k]
                   + f_0 * kk_346[k];

        t_203[k] = -2.0 * hk_95[k]
                   + f_0 * kk_347[k];

        t_204[k] = -2.0 * hk_96[k]
                   + f_0 * kk_348[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, hk_97, hk_98, hk_99, hk_100, \
                         hk_101, kk_349, kk_350, kk_351, kk_352, \
                         kk_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * hk_97[k]
                   + f_0 * kk_349[k];

        t_206[k] = -2.0 * hk_98[k]
                   + f_0 * kk_350[k];

        t_207[k] = -2.0 * hk_99[k]
                   + f_0 * kk_351[k];

        t_208[k] = -2.0 * hk_100[k]
                   + f_0 * kk_352[k];

        t_209[k] = -2.0 * hk_101[k]
                   + f_0 * kk_353[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, hk_102, hk_103, hk_104, hk_105, \
                         hk_106, kk_354, kk_355, kk_356, kk_357, \
                         kk_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * hk_102[k]
                   + f_0 * kk_354[k];

        t_211[k] = -2.0 * hk_103[k]
                   + f_0 * kk_355[k];

        t_212[k] = -2.0 * hk_104[k]
                   + f_0 * kk_356[k];

        t_213[k] = -2.0 * hk_105[k]
                   + f_0 * kk_357[k];

        t_214[k] = -2.0 * hk_106[k]
                   + f_0 * kk_358[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, t_220, t_221, hk_107, kk_359, \
                         kk_396, kk_397, kk_398, kk_399, kk_400, \
                         kk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * hk_107[k]
                   + f_0 * kk_359[k];

        t_216[k] = f_0 * kk_396[k];

        t_217[k] = f_0 * kk_397[k];

        t_218[k] = f_0 * kk_398[k];

        t_219[k] = f_0 * kk_399[k];

        t_220[k] = f_0 * kk_400[k];

        t_221[k] = f_0 * kk_401[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, t_228, t_229, kk_402, \
                         kk_403, kk_404, kk_405, kk_406, kk_407, kk_408, \
                         kk_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_0 * kk_402[k];

        t_223[k] = f_0 * kk_403[k];

        t_224[k] = f_0 * kk_404[k];

        t_225[k] = f_0 * kk_405[k];

        t_226[k] = f_0 * kk_406[k];

        t_227[k] = f_0 * kk_407[k];

        t_228[k] = f_0 * kk_408[k];

        t_229[k] = f_0 * kk_409[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, t_235, t_236, t_237, kk_410, \
                         kk_411, kk_412, kk_413, kk_414, kk_415, kk_416, \
                         kk_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_0 * kk_410[k];

        t_231[k] = f_0 * kk_411[k];

        t_232[k] = f_0 * kk_412[k];

        t_233[k] = f_0 * kk_413[k];

        t_234[k] = f_0 * kk_414[k];

        t_235[k] = f_0 * kk_415[k];

        t_236[k] = f_0 * kk_416[k];

        t_237[k] = f_0 * kk_417[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, t_243, t_244, t_245, kk_418, \
                         kk_419, kk_420, kk_421, kk_422, kk_423, kk_424, \
                         kk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_0 * kk_418[k];

        t_239[k] = f_0 * kk_419[k];

        t_240[k] = f_0 * kk_420[k];

        t_241[k] = f_0 * kk_421[k];

        t_242[k] = f_0 * kk_422[k];

        t_243[k] = f_0 * kk_423[k];

        t_244[k] = f_0 * kk_424[k];

        t_245[k] = f_0 * kk_425[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, t_251, t_252, hk_108, kk_426, \
                         kk_427, kk_428, kk_429, kk_430, kk_431, \
                         kk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_0 * kk_426[k];

        t_247[k] = f_0 * kk_427[k];

        t_248[k] = f_0 * kk_428[k];

        t_249[k] = f_0 * kk_429[k];

        t_250[k] = f_0 * kk_430[k];

        t_251[k] = f_0 * kk_431[k];

        t_252[k] = -hk_108[k]
                   + f_0 * kk_432[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, hk_109, hk_110, hk_111, hk_112, \
                         hk_113, kk_433, kk_434, kk_435, kk_436, \
                         kk_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = -hk_109[k]
                   + f_0 * kk_433[k];

        t_254[k] = -hk_110[k]
                   + f_0 * kk_434[k];

        t_255[k] = -hk_111[k]
                   + f_0 * kk_435[k];

        t_256[k] = -hk_112[k]
                   + f_0 * kk_436[k];

        t_257[k] = -hk_113[k]
                   + f_0 * kk_437[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, hk_114, hk_115, hk_116, hk_117, \
                         hk_118, kk_438, kk_439, kk_440, kk_441, \
                         kk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = -hk_114[k]
                   + f_0 * kk_438[k];

        t_259[k] = -hk_115[k]
                   + f_0 * kk_439[k];

        t_260[k] = -hk_116[k]
                   + f_0 * kk_440[k];

        t_261[k] = -hk_117[k]
                   + f_0 * kk_441[k];

        t_262[k] = -hk_118[k]
                   + f_0 * kk_442[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, hk_119, hk_120, hk_121, hk_122, \
                         hk_123, kk_443, kk_444, kk_445, kk_446, \
                         kk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = -hk_119[k]
                   + f_0 * kk_443[k];

        t_264[k] = -hk_120[k]
                   + f_0 * kk_444[k];

        t_265[k] = -hk_121[k]
                   + f_0 * kk_445[k];

        t_266[k] = -hk_122[k]
                   + f_0 * kk_446[k];

        t_267[k] = -hk_123[k]
                   + f_0 * kk_447[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, hk_124, hk_125, hk_126, hk_127, \
                         hk_128, kk_448, kk_449, kk_450, kk_451, \
                         kk_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = -hk_124[k]
                   + f_0 * kk_448[k];

        t_269[k] = -hk_125[k]
                   + f_0 * kk_449[k];

        t_270[k] = -hk_126[k]
                   + f_0 * kk_450[k];

        t_271[k] = -hk_127[k]
                   + f_0 * kk_451[k];

        t_272[k] = -hk_128[k]
                   + f_0 * kk_452[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, hk_129, hk_130, hk_131, hk_132, \
                         hk_133, kk_453, kk_454, kk_455, kk_456, \
                         kk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = -hk_129[k]
                   + f_0 * kk_453[k];

        t_274[k] = -hk_130[k]
                   + f_0 * kk_454[k];

        t_275[k] = -hk_131[k]
                   + f_0 * kk_455[k];

        t_276[k] = -hk_132[k]
                   + f_0 * kk_456[k];

        t_277[k] = -hk_133[k]
                   + f_0 * kk_457[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, hk_134, hk_135, hk_136, hk_137, \
                         hk_138, kk_458, kk_459, kk_460, kk_461, \
                         kk_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = -hk_134[k]
                   + f_0 * kk_458[k];

        t_279[k] = -hk_135[k]
                   + f_0 * kk_459[k];

        t_280[k] = -hk_136[k]
                   + f_0 * kk_460[k];

        t_281[k] = -hk_137[k]
                   + f_0 * kk_461[k];

        t_282[k] = -hk_138[k]
                   + f_0 * kk_462[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, hk_139, hk_140, hk_141, hk_142, \
                         hk_143, kk_463, kk_464, kk_465, kk_466, \
                         kk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = -hk_139[k]
                   + f_0 * kk_463[k];

        t_284[k] = -hk_140[k]
                   + f_0 * kk_464[k];

        t_285[k] = -hk_141[k]
                   + f_0 * kk_465[k];

        t_286[k] = -hk_142[k]
                   + f_0 * kk_466[k];

        t_287[k] = -hk_143[k]
                   + f_0 * kk_467[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, hk_144, hk_145, hk_146, hk_147, \
                         hk_148, kk_468, kk_469, kk_470, kk_471, \
                         kk_472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = -2.0 * hk_144[k]
                   + f_0 * kk_468[k];

        t_289[k] = -2.0 * hk_145[k]
                   + f_0 * kk_469[k];

        t_290[k] = -2.0 * hk_146[k]
                   + f_0 * kk_470[k];

        t_291[k] = -2.0 * hk_147[k]
                   + f_0 * kk_471[k];

        t_292[k] = -2.0 * hk_148[k]
                   + f_0 * kk_472[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, t_297, hk_149, hk_150, hk_151, hk_152, \
                         hk_153, kk_473, kk_474, kk_475, kk_476, \
                         kk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = -2.0 * hk_149[k]
                   + f_0 * kk_473[k];

        t_294[k] = -2.0 * hk_150[k]
                   + f_0 * kk_474[k];

        t_295[k] = -2.0 * hk_151[k]
                   + f_0 * kk_475[k];

        t_296[k] = -2.0 * hk_152[k]
                   + f_0 * kk_476[k];

        t_297[k] = -2.0 * hk_153[k]
                   + f_0 * kk_477[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, hk_154, hk_155, hk_156, hk_157, \
                         hk_158, kk_478, kk_479, kk_480, kk_481, \
                         kk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = -2.0 * hk_154[k]
                   + f_0 * kk_478[k];

        t_299[k] = -2.0 * hk_155[k]
                   + f_0 * kk_479[k];

        t_300[k] = -2.0 * hk_156[k]
                   + f_0 * kk_480[k];

        t_301[k] = -2.0 * hk_157[k]
                   + f_0 * kk_481[k];

        t_302[k] = -2.0 * hk_158[k]
                   + f_0 * kk_482[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, hk_159, hk_160, hk_161, hk_162, \
                         hk_163, kk_483, kk_484, kk_485, kk_486, \
                         kk_487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = -2.0 * hk_159[k]
                   + f_0 * kk_483[k];

        t_304[k] = -2.0 * hk_160[k]
                   + f_0 * kk_484[k];

        t_305[k] = -2.0 * hk_161[k]
                   + f_0 * kk_485[k];

        t_306[k] = -2.0 * hk_162[k]
                   + f_0 * kk_486[k];

        t_307[k] = -2.0 * hk_163[k]
                   + f_0 * kk_487[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, hk_164, hk_165, hk_166, hk_167, \
                         hk_168, kk_488, kk_489, kk_490, kk_491, \
                         kk_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = -2.0 * hk_164[k]
                   + f_0 * kk_488[k];

        t_309[k] = -2.0 * hk_165[k]
                   + f_0 * kk_489[k];

        t_310[k] = -2.0 * hk_166[k]
                   + f_0 * kk_490[k];

        t_311[k] = -2.0 * hk_167[k]
                   + f_0 * kk_491[k];

        t_312[k] = -2.0 * hk_168[k]
                   + f_0 * kk_492[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, hk_169, hk_170, hk_171, hk_172, \
                         hk_173, kk_493, kk_494, kk_495, kk_496, \
                         kk_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = -2.0 * hk_169[k]
                   + f_0 * kk_493[k];

        t_314[k] = -2.0 * hk_170[k]
                   + f_0 * kk_494[k];

        t_315[k] = -2.0 * hk_171[k]
                   + f_0 * kk_495[k];

        t_316[k] = -2.0 * hk_172[k]
                   + f_0 * kk_496[k];

        t_317[k] = -2.0 * hk_173[k]
                   + f_0 * kk_497[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, hk_174, hk_175, hk_176, hk_177, \
                         hk_178, kk_498, kk_499, kk_500, kk_501, \
                         kk_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = -2.0 * hk_174[k]
                   + f_0 * kk_498[k];

        t_319[k] = -2.0 * hk_175[k]
                   + f_0 * kk_499[k];

        t_320[k] = -2.0 * hk_176[k]
                   + f_0 * kk_500[k];

        t_321[k] = -2.0 * hk_177[k]
                   + f_0 * kk_501[k];

        t_322[k] = -2.0 * hk_178[k]
                   + f_0 * kk_502[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, hk_179, hk_180, hk_181, hk_182, \
                         hk_183, kk_503, kk_504, kk_505, kk_506, \
                         kk_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = -2.0 * hk_179[k]
                   + f_0 * kk_503[k];

        t_324[k] = -3.0 * hk_180[k]
                   + f_0 * kk_504[k];

        t_325[k] = -3.0 * hk_181[k]
                   + f_0 * kk_505[k];

        t_326[k] = -3.0 * hk_182[k]
                   + f_0 * kk_506[k];

        t_327[k] = -3.0 * hk_183[k]
                   + f_0 * kk_507[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, hk_184, hk_185, hk_186, hk_187, \
                         hk_188, kk_508, kk_509, kk_510, kk_511, \
                         kk_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = -3.0 * hk_184[k]
                   + f_0 * kk_508[k];

        t_329[k] = -3.0 * hk_185[k]
                   + f_0 * kk_509[k];

        t_330[k] = -3.0 * hk_186[k]
                   + f_0 * kk_510[k];

        t_331[k] = -3.0 * hk_187[k]
                   + f_0 * kk_511[k];

        t_332[k] = -3.0 * hk_188[k]
                   + f_0 * kk_512[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, t_337, hk_189, hk_190, hk_191, hk_192, \
                         hk_193, kk_513, kk_514, kk_515, kk_516, \
                         kk_517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = -3.0 * hk_189[k]
                   + f_0 * kk_513[k];

        t_334[k] = -3.0 * hk_190[k]
                   + f_0 * kk_514[k];

        t_335[k] = -3.0 * hk_191[k]
                   + f_0 * kk_515[k];

        t_336[k] = -3.0 * hk_192[k]
                   + f_0 * kk_516[k];

        t_337[k] = -3.0 * hk_193[k]
                   + f_0 * kk_517[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, hk_194, hk_195, hk_196, hk_197, \
                         hk_198, kk_518, kk_519, kk_520, kk_521, \
                         kk_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = -3.0 * hk_194[k]
                   + f_0 * kk_518[k];

        t_339[k] = -3.0 * hk_195[k]
                   + f_0 * kk_519[k];

        t_340[k] = -3.0 * hk_196[k]
                   + f_0 * kk_520[k];

        t_341[k] = -3.0 * hk_197[k]
                   + f_0 * kk_521[k];

        t_342[k] = -3.0 * hk_198[k]
                   + f_0 * kk_522[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, hk_199, hk_200, hk_201, hk_202, \
                         hk_203, kk_523, kk_524, kk_525, kk_526, \
                         kk_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = -3.0 * hk_199[k]
                   + f_0 * kk_523[k];

        t_344[k] = -3.0 * hk_200[k]
                   + f_0 * kk_524[k];

        t_345[k] = -3.0 * hk_201[k]
                   + f_0 * kk_525[k];

        t_346[k] = -3.0 * hk_202[k]
                   + f_0 * kk_526[k];

        t_347[k] = -3.0 * hk_203[k]
                   + f_0 * kk_527[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_204 = buffer.data(hk + 204);
    const auto *hk_205 = buffer.data(hk + 205);
    const auto *hk_206 = buffer.data(hk + 206);
    const auto *hk_207 = buffer.data(hk + 207);
    const auto *hk_208 = buffer.data(hk + 208);
    const auto *hk_209 = buffer.data(hk + 209);
    const auto *hk_210 = buffer.data(hk + 210);
    const auto *hk_211 = buffer.data(hk + 211);
    const auto *hk_212 = buffer.data(hk + 212);
    const auto *hk_213 = buffer.data(hk + 213);
    const auto *hk_214 = buffer.data(hk + 214);
    const auto *hk_215 = buffer.data(hk + 215);
    const auto *hk_216 = buffer.data(hk + 216);
    const auto *hk_217 = buffer.data(hk + 217);
    const auto *hk_218 = buffer.data(hk + 218);
    const auto *hk_219 = buffer.data(hk + 219);
    const auto *hk_220 = buffer.data(hk + 220);
    const auto *hk_221 = buffer.data(hk + 221);
    const auto *hk_222 = buffer.data(hk + 222);
    const auto *hk_223 = buffer.data(hk + 223);
    const auto *hk_224 = buffer.data(hk + 224);
    const auto *hk_225 = buffer.data(hk + 225);
    const auto *hk_226 = buffer.data(hk + 226);
    const auto *hk_227 = buffer.data(hk + 227);
    const auto *hk_228 = buffer.data(hk + 228);
    const auto *hk_229 = buffer.data(hk + 229);
    const auto *hk_230 = buffer.data(hk + 230);
    const auto *hk_231 = buffer.data(hk + 231);
    const auto *hk_232 = buffer.data(hk + 232);
    const auto *hk_233 = buffer.data(hk + 233);
    const auto *hk_234 = buffer.data(hk + 234);
    const auto *hk_235 = buffer.data(hk + 235);
    const auto *hk_236 = buffer.data(hk + 236);
    const auto *hk_237 = buffer.data(hk + 237);
    const auto *hk_238 = buffer.data(hk + 238);
    const auto *hk_239 = buffer.data(hk + 239);
    const auto *hk_240 = buffer.data(hk + 240);
    const auto *hk_241 = buffer.data(hk + 241);
    const auto *hk_242 = buffer.data(hk + 242);
    const auto *hk_243 = buffer.data(hk + 243);
    const auto *hk_244 = buffer.data(hk + 244);
    const auto *hk_245 = buffer.data(hk + 245);
    const auto *hk_246 = buffer.data(hk + 246);
    const auto *hk_247 = buffer.data(hk + 247);
    const auto *hk_248 = buffer.data(hk + 248);
    const auto *hk_249 = buffer.data(hk + 249);
    const auto *hk_250 = buffer.data(hk + 250);
    const auto *hk_251 = buffer.data(hk + 251);
    const auto *hk_252 = buffer.data(hk + 252);
    const auto *hk_253 = buffer.data(hk + 253);
    const auto *hk_254 = buffer.data(hk + 254);
    const auto *hk_255 = buffer.data(hk + 255);
    const auto *hk_256 = buffer.data(hk + 256);
    const auto *hk_257 = buffer.data(hk + 257);
    const auto *hk_258 = buffer.data(hk + 258);
    const auto *hk_259 = buffer.data(hk + 259);
    const auto *hk_260 = buffer.data(hk + 260);
    const auto *hk_261 = buffer.data(hk + 261);
    const auto *hk_262 = buffer.data(hk + 262);
    const auto *hk_263 = buffer.data(hk + 263);
    const auto *hk_264 = buffer.data(hk + 264);
    const auto *hk_265 = buffer.data(hk + 265);
    const auto *hk_266 = buffer.data(hk + 266);
    const auto *hk_267 = buffer.data(hk + 267);
    const auto *hk_268 = buffer.data(hk + 268);
    const auto *hk_269 = buffer.data(hk + 269);
    const auto *hk_270 = buffer.data(hk + 270);
    const auto *hk_271 = buffer.data(hk + 271);
    const auto *hk_272 = buffer.data(hk + 272);
    const auto *hk_273 = buffer.data(hk + 273);
    const auto *hk_274 = buffer.data(hk + 274);
    const auto *hk_275 = buffer.data(hk + 275);
    const auto *hk_276 = buffer.data(hk + 276);
    const auto *hk_277 = buffer.data(hk + 277);
    const auto *hk_278 = buffer.data(hk + 278);
    const auto *hk_279 = buffer.data(hk + 279);
    const auto *hk_280 = buffer.data(hk + 280);
    const auto *hk_281 = buffer.data(hk + 281);
    const auto *hk_282 = buffer.data(hk + 282);
    const auto *hk_283 = buffer.data(hk + 283);
    const auto *hk_284 = buffer.data(hk + 284);
    const auto *hk_285 = buffer.data(hk + 285);
    const auto *hk_286 = buffer.data(hk + 286);
    const auto *hk_287 = buffer.data(hk + 287);
    const auto *hk_288 = buffer.data(hk + 288);
    const auto *hk_289 = buffer.data(hk + 289);
    const auto *hk_290 = buffer.data(hk + 290);
    const auto *hk_291 = buffer.data(hk + 291);
    const auto *hk_292 = buffer.data(hk + 292);
    const auto *hk_293 = buffer.data(hk + 293);
    const auto *hk_294 = buffer.data(hk + 294);
    const auto *hk_295 = buffer.data(hk + 295);
    const auto *hk_296 = buffer.data(hk + 296);
    const auto *hk_297 = buffer.data(hk + 297);
    const auto *hk_298 = buffer.data(hk + 298);
    const auto *hk_299 = buffer.data(hk + 299);
    const auto *hk_300 = buffer.data(hk + 300);
    const auto *hk_301 = buffer.data(hk + 301);
    const auto *hk_302 = buffer.data(hk + 302);
    const auto *hk_303 = buffer.data(hk + 303);
    const auto *hk_304 = buffer.data(hk + 304);
    const auto *hk_305 = buffer.data(hk + 305);
    const auto *hk_306 = buffer.data(hk + 306);
    const auto *hk_307 = buffer.data(hk + 307);
    const auto *hk_308 = buffer.data(hk + 308);
    const auto *hk_309 = buffer.data(hk + 309);
    const auto *hk_310 = buffer.data(hk + 310);
    const auto *hk_311 = buffer.data(hk + 311);
    const auto *hk_312 = buffer.data(hk + 312);
    const auto *hk_313 = buffer.data(hk + 313);
    const auto *hk_314 = buffer.data(hk + 314);
    const auto *hk_315 = buffer.data(hk + 315);
    const auto *hk_316 = buffer.data(hk + 316);
    const auto *hk_317 = buffer.data(hk + 317);
    const auto *hk_318 = buffer.data(hk + 318);
    const auto *hk_319 = buffer.data(hk + 319);
    const auto *hk_320 = buffer.data(hk + 320);
    const auto *hk_321 = buffer.data(hk + 321);
    const auto *hk_322 = buffer.data(hk + 322);
    const auto *hk_323 = buffer.data(hk + 323);
    const auto *hk_324 = buffer.data(hk + 324);
    const auto *hk_325 = buffer.data(hk + 325);
    const auto *hk_326 = buffer.data(hk + 326);
    const auto *hk_327 = buffer.data(hk + 327);
    const auto *hk_328 = buffer.data(hk + 328);
    const auto *hk_329 = buffer.data(hk + 329);
    const auto *hk_330 = buffer.data(hk + 330);

    const auto *kk_528 = buffer.data(kk + 528);
    const auto *kk_529 = buffer.data(kk + 529);
    const auto *kk_530 = buffer.data(kk + 530);
    const auto *kk_531 = buffer.data(kk + 531);
    const auto *kk_532 = buffer.data(kk + 532);
    const auto *kk_533 = buffer.data(kk + 533);
    const auto *kk_534 = buffer.data(kk + 534);
    const auto *kk_535 = buffer.data(kk + 535);
    const auto *kk_536 = buffer.data(kk + 536);
    const auto *kk_537 = buffer.data(kk + 537);
    const auto *kk_538 = buffer.data(kk + 538);
    const auto *kk_539 = buffer.data(kk + 539);
    const auto *kk_576 = buffer.data(kk + 576);
    const auto *kk_577 = buffer.data(kk + 577);
    const auto *kk_578 = buffer.data(kk + 578);
    const auto *kk_579 = buffer.data(kk + 579);
    const auto *kk_580 = buffer.data(kk + 580);
    const auto *kk_581 = buffer.data(kk + 581);
    const auto *kk_582 = buffer.data(kk + 582);
    const auto *kk_583 = buffer.data(kk + 583);
    const auto *kk_584 = buffer.data(kk + 584);
    const auto *kk_585 = buffer.data(kk + 585);
    const auto *kk_586 = buffer.data(kk + 586);
    const auto *kk_587 = buffer.data(kk + 587);
    const auto *kk_588 = buffer.data(kk + 588);
    const auto *kk_589 = buffer.data(kk + 589);
    const auto *kk_590 = buffer.data(kk + 590);
    const auto *kk_591 = buffer.data(kk + 591);
    const auto *kk_592 = buffer.data(kk + 592);
    const auto *kk_593 = buffer.data(kk + 593);
    const auto *kk_594 = buffer.data(kk + 594);
    const auto *kk_595 = buffer.data(kk + 595);
    const auto *kk_596 = buffer.data(kk + 596);
    const auto *kk_597 = buffer.data(kk + 597);
    const auto *kk_598 = buffer.data(kk + 598);
    const auto *kk_599 = buffer.data(kk + 599);
    const auto *kk_600 = buffer.data(kk + 600);
    const auto *kk_601 = buffer.data(kk + 601);
    const auto *kk_602 = buffer.data(kk + 602);
    const auto *kk_603 = buffer.data(kk + 603);
    const auto *kk_604 = buffer.data(kk + 604);
    const auto *kk_605 = buffer.data(kk + 605);
    const auto *kk_606 = buffer.data(kk + 606);
    const auto *kk_607 = buffer.data(kk + 607);
    const auto *kk_608 = buffer.data(kk + 608);
    const auto *kk_609 = buffer.data(kk + 609);
    const auto *kk_610 = buffer.data(kk + 610);
    const auto *kk_611 = buffer.data(kk + 611);
    const auto *kk_612 = buffer.data(kk + 612);
    const auto *kk_613 = buffer.data(kk + 613);
    const auto *kk_614 = buffer.data(kk + 614);
    const auto *kk_615 = buffer.data(kk + 615);
    const auto *kk_616 = buffer.data(kk + 616);
    const auto *kk_617 = buffer.data(kk + 617);
    const auto *kk_618 = buffer.data(kk + 618);
    const auto *kk_619 = buffer.data(kk + 619);
    const auto *kk_620 = buffer.data(kk + 620);
    const auto *kk_621 = buffer.data(kk + 621);
    const auto *kk_622 = buffer.data(kk + 622);
    const auto *kk_623 = buffer.data(kk + 623);
    const auto *kk_624 = buffer.data(kk + 624);
    const auto *kk_625 = buffer.data(kk + 625);
    const auto *kk_626 = buffer.data(kk + 626);
    const auto *kk_627 = buffer.data(kk + 627);
    const auto *kk_628 = buffer.data(kk + 628);
    const auto *kk_629 = buffer.data(kk + 629);
    const auto *kk_630 = buffer.data(kk + 630);
    const auto *kk_631 = buffer.data(kk + 631);
    const auto *kk_632 = buffer.data(kk + 632);
    const auto *kk_633 = buffer.data(kk + 633);
    const auto *kk_634 = buffer.data(kk + 634);
    const auto *kk_635 = buffer.data(kk + 635);
    const auto *kk_636 = buffer.data(kk + 636);
    const auto *kk_637 = buffer.data(kk + 637);
    const auto *kk_638 = buffer.data(kk + 638);
    const auto *kk_639 = buffer.data(kk + 639);
    const auto *kk_640 = buffer.data(kk + 640);
    const auto *kk_641 = buffer.data(kk + 641);
    const auto *kk_642 = buffer.data(kk + 642);
    const auto *kk_643 = buffer.data(kk + 643);
    const auto *kk_644 = buffer.data(kk + 644);
    const auto *kk_645 = buffer.data(kk + 645);
    const auto *kk_646 = buffer.data(kk + 646);
    const auto *kk_647 = buffer.data(kk + 647);
    const auto *kk_648 = buffer.data(kk + 648);
    const auto *kk_649 = buffer.data(kk + 649);
    const auto *kk_650 = buffer.data(kk + 650);
    const auto *kk_651 = buffer.data(kk + 651);
    const auto *kk_652 = buffer.data(kk + 652);
    const auto *kk_653 = buffer.data(kk + 653);
    const auto *kk_654 = buffer.data(kk + 654);
    const auto *kk_655 = buffer.data(kk + 655);
    const auto *kk_656 = buffer.data(kk + 656);
    const auto *kk_657 = buffer.data(kk + 657);
    const auto *kk_658 = buffer.data(kk + 658);
    const auto *kk_659 = buffer.data(kk + 659);
    const auto *kk_660 = buffer.data(kk + 660);
    const auto *kk_661 = buffer.data(kk + 661);
    const auto *kk_662 = buffer.data(kk + 662);
    const auto *kk_663 = buffer.data(kk + 663);
    const auto *kk_664 = buffer.data(kk + 664);
    const auto *kk_665 = buffer.data(kk + 665);
    const auto *kk_666 = buffer.data(kk + 666);
    const auto *kk_667 = buffer.data(kk + 667);
    const auto *kk_668 = buffer.data(kk + 668);
    const auto *kk_669 = buffer.data(kk + 669);
    const auto *kk_670 = buffer.data(kk + 670);
    const auto *kk_671 = buffer.data(kk + 671);
    const auto *kk_672 = buffer.data(kk + 672);
    const auto *kk_673 = buffer.data(kk + 673);
    const auto *kk_674 = buffer.data(kk + 674);
    const auto *kk_675 = buffer.data(kk + 675);
    const auto *kk_676 = buffer.data(kk + 676);
    const auto *kk_677 = buffer.data(kk + 677);
    const auto *kk_678 = buffer.data(kk + 678);
    const auto *kk_679 = buffer.data(kk + 679);
    const auto *kk_680 = buffer.data(kk + 680);
    const auto *kk_681 = buffer.data(kk + 681);
    const auto *kk_682 = buffer.data(kk + 682);
    const auto *kk_683 = buffer.data(kk + 683);
    const auto *kk_684 = buffer.data(kk + 684);
    const auto *kk_685 = buffer.data(kk + 685);
    const auto *kk_686 = buffer.data(kk + 686);
    const auto *kk_687 = buffer.data(kk + 687);
    const auto *kk_688 = buffer.data(kk + 688);
    const auto *kk_689 = buffer.data(kk + 689);
    const auto *kk_690 = buffer.data(kk + 690);
    const auto *kk_691 = buffer.data(kk + 691);
    const auto *kk_692 = buffer.data(kk + 692);
    const auto *kk_693 = buffer.data(kk + 693);
    const auto *kk_694 = buffer.data(kk + 694);
    const auto *kk_695 = buffer.data(kk + 695);
    const auto *kk_696 = buffer.data(kk + 696);
    const auto *kk_697 = buffer.data(kk + 697);
    const auto *kk_698 = buffer.data(kk + 698);
    const auto *kk_699 = buffer.data(kk + 699);
    const auto *kk_700 = buffer.data(kk + 700);
    const auto *kk_701 = buffer.data(kk + 701);
    const auto *kk_702 = buffer.data(kk + 702);
    const auto *kk_703 = buffer.data(kk + 703);
    const auto *kk_704 = buffer.data(kk + 704);
    const auto *kk_705 = buffer.data(kk + 705);
    const auto *kk_706 = buffer.data(kk + 706);
    const auto *kk_707 = buffer.data(kk + 707);
    const auto *kk_708 = buffer.data(kk + 708);
    const auto *kk_709 = buffer.data(kk + 709);
    const auto *kk_710 = buffer.data(kk + 710);
    const auto *kk_711 = buffer.data(kk + 711);
    const auto *kk_712 = buffer.data(kk + 712);
    const auto *kk_713 = buffer.data(kk + 713);
    const auto *kk_714 = buffer.data(kk + 714);
    const auto *kk_715 = buffer.data(kk + 715);
    const auto *kk_716 = buffer.data(kk + 716);
    const auto *kk_717 = buffer.data(kk + 717);
    const auto *kk_718 = buffer.data(kk + 718);
    const auto *kk_719 = buffer.data(kk + 719);
    const auto *kk_720 = buffer.data(kk + 720);
    const auto *kk_721 = buffer.data(kk + 721);
    const auto *kk_722 = buffer.data(kk + 722);
    const auto *kk_723 = buffer.data(kk + 723);
    const auto *kk_724 = buffer.data(kk + 724);
    const auto *kk_725 = buffer.data(kk + 725);
    const auto *kk_726 = buffer.data(kk + 726);

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, hk_204, hk_205, hk_206, hk_207, \
                         hk_208, kk_528, kk_529, kk_530, kk_531, \
                         kk_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = -3.0 * hk_204[k]
                   + f_0 * kk_528[k];

        t_349[k] = -3.0 * hk_205[k]
                   + f_0 * kk_529[k];

        t_350[k] = -3.0 * hk_206[k]
                   + f_0 * kk_530[k];

        t_351[k] = -3.0 * hk_207[k]
                   + f_0 * kk_531[k];

        t_352[k] = -3.0 * hk_208[k]
                   + f_0 * kk_532[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, hk_209, hk_210, hk_211, hk_212, \
                         hk_213, kk_533, kk_534, kk_535, kk_536, \
                         kk_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = -3.0 * hk_209[k]
                   + f_0 * kk_533[k];

        t_354[k] = -3.0 * hk_210[k]
                   + f_0 * kk_534[k];

        t_355[k] = -3.0 * hk_211[k]
                   + f_0 * kk_535[k];

        t_356[k] = -3.0 * hk_212[k]
                   + f_0 * kk_536[k];

        t_357[k] = -3.0 * hk_213[k]
                   + f_0 * kk_537[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, t_363, t_364, hk_214, hk_215, \
                         kk_538, kk_539, kk_576, kk_577, kk_578, kk_579, \
                         kk_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = -3.0 * hk_214[k]
                   + f_0 * kk_538[k];

        t_359[k] = -3.0 * hk_215[k]
                   + f_0 * kk_539[k];

        t_360[k] = f_0 * kk_576[k];

        t_361[k] = f_0 * kk_577[k];

        t_362[k] = f_0 * kk_578[k];

        t_363[k] = f_0 * kk_579[k];

        t_364[k] = f_0 * kk_580[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, t_370, t_371, t_372, kk_581, \
                         kk_582, kk_583, kk_584, kk_585, kk_586, kk_587, \
                         kk_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_0 * kk_581[k];

        t_366[k] = f_0 * kk_582[k];

        t_367[k] = f_0 * kk_583[k];

        t_368[k] = f_0 * kk_584[k];

        t_369[k] = f_0 * kk_585[k];

        t_370[k] = f_0 * kk_586[k];

        t_371[k] = f_0 * kk_587[k];

        t_372[k] = f_0 * kk_588[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, t_378, t_379, t_380, kk_589, \
                         kk_590, kk_591, kk_592, kk_593, kk_594, kk_595, \
                         kk_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_0 * kk_589[k];

        t_374[k] = f_0 * kk_590[k];

        t_375[k] = f_0 * kk_591[k];

        t_376[k] = f_0 * kk_592[k];

        t_377[k] = f_0 * kk_593[k];

        t_378[k] = f_0 * kk_594[k];

        t_379[k] = f_0 * kk_595[k];

        t_380[k] = f_0 * kk_596[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, t_385, t_386, t_387, t_388, kk_597, \
                         kk_598, kk_599, kk_600, kk_601, kk_602, kk_603, \
                         kk_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_0 * kk_597[k];

        t_382[k] = f_0 * kk_598[k];

        t_383[k] = f_0 * kk_599[k];

        t_384[k] = f_0 * kk_600[k];

        t_385[k] = f_0 * kk_601[k];

        t_386[k] = f_0 * kk_602[k];

        t_387[k] = f_0 * kk_603[k];

        t_388[k] = f_0 * kk_604[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, t_394, t_395, kk_605, kk_606, \
                         kk_607, kk_608, kk_609, kk_610, kk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_0 * kk_605[k];

        t_390[k] = f_0 * kk_606[k];

        t_391[k] = f_0 * kk_607[k];

        t_392[k] = f_0 * kk_608[k];

        t_393[k] = f_0 * kk_609[k];

        t_394[k] = f_0 * kk_610[k];

        t_395[k] = f_0 * kk_611[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, t_400, hk_216, hk_217, hk_218, hk_219, \
                         hk_220, kk_612, kk_613, kk_614, kk_615, \
                         kk_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = -hk_216[k]
                   + f_0 * kk_612[k];

        t_397[k] = -hk_217[k]
                   + f_0 * kk_613[k];

        t_398[k] = -hk_218[k]
                   + f_0 * kk_614[k];

        t_399[k] = -hk_219[k]
                   + f_0 * kk_615[k];

        t_400[k] = -hk_220[k]
                   + f_0 * kk_616[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, t_405, hk_221, hk_222, hk_223, hk_224, \
                         hk_225, kk_617, kk_618, kk_619, kk_620, \
                         kk_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = -hk_221[k]
                   + f_0 * kk_617[k];

        t_402[k] = -hk_222[k]
                   + f_0 * kk_618[k];

        t_403[k] = -hk_223[k]
                   + f_0 * kk_619[k];

        t_404[k] = -hk_224[k]
                   + f_0 * kk_620[k];

        t_405[k] = -hk_225[k]
                   + f_0 * kk_621[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, t_410, hk_226, hk_227, hk_228, hk_229, \
                         hk_230, kk_622, kk_623, kk_624, kk_625, \
                         kk_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = -hk_226[k]
                   + f_0 * kk_622[k];

        t_407[k] = -hk_227[k]
                   + f_0 * kk_623[k];

        t_408[k] = -hk_228[k]
                   + f_0 * kk_624[k];

        t_409[k] = -hk_229[k]
                   + f_0 * kk_625[k];

        t_410[k] = -hk_230[k]
                   + f_0 * kk_626[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, t_415, hk_231, hk_232, hk_233, hk_234, \
                         hk_235, kk_627, kk_628, kk_629, kk_630, \
                         kk_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = -hk_231[k]
                   + f_0 * kk_627[k];

        t_412[k] = -hk_232[k]
                   + f_0 * kk_628[k];

        t_413[k] = -hk_233[k]
                   + f_0 * kk_629[k];

        t_414[k] = -hk_234[k]
                   + f_0 * kk_630[k];

        t_415[k] = -hk_235[k]
                   + f_0 * kk_631[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, t_420, hk_236, hk_237, hk_238, hk_239, \
                         hk_240, kk_632, kk_633, kk_634, kk_635, \
                         kk_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = -hk_236[k]
                   + f_0 * kk_632[k];

        t_417[k] = -hk_237[k]
                   + f_0 * kk_633[k];

        t_418[k] = -hk_238[k]
                   + f_0 * kk_634[k];

        t_419[k] = -hk_239[k]
                   + f_0 * kk_635[k];

        t_420[k] = -hk_240[k]
                   + f_0 * kk_636[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, t_425, hk_241, hk_242, hk_243, hk_244, \
                         hk_245, kk_637, kk_638, kk_639, kk_640, \
                         kk_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = -hk_241[k]
                   + f_0 * kk_637[k];

        t_422[k] = -hk_242[k]
                   + f_0 * kk_638[k];

        t_423[k] = -hk_243[k]
                   + f_0 * kk_639[k];

        t_424[k] = -hk_244[k]
                   + f_0 * kk_640[k];

        t_425[k] = -hk_245[k]
                   + f_0 * kk_641[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, t_430, hk_246, hk_247, hk_248, hk_249, \
                         hk_250, kk_642, kk_643, kk_644, kk_645, \
                         kk_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = -hk_246[k]
                   + f_0 * kk_642[k];

        t_427[k] = -hk_247[k]
                   + f_0 * kk_643[k];

        t_428[k] = -hk_248[k]
                   + f_0 * kk_644[k];

        t_429[k] = -hk_249[k]
                   + f_0 * kk_645[k];

        t_430[k] = -hk_250[k]
                   + f_0 * kk_646[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, hk_251, hk_252, hk_253, hk_254, \
                         hk_255, kk_647, kk_648, kk_649, kk_650, \
                         kk_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = -hk_251[k]
                   + f_0 * kk_647[k];

        t_432[k] = -2.0 * hk_252[k]
                   + f_0 * kk_648[k];

        t_433[k] = -2.0 * hk_253[k]
                   + f_0 * kk_649[k];

        t_434[k] = -2.0 * hk_254[k]
                   + f_0 * kk_650[k];

        t_435[k] = -2.0 * hk_255[k]
                   + f_0 * kk_651[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, hk_256, hk_257, hk_258, hk_259, \
                         hk_260, kk_652, kk_653, kk_654, kk_655, \
                         kk_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = -2.0 * hk_256[k]
                   + f_0 * kk_652[k];

        t_437[k] = -2.0 * hk_257[k]
                   + f_0 * kk_653[k];

        t_438[k] = -2.0 * hk_258[k]
                   + f_0 * kk_654[k];

        t_439[k] = -2.0 * hk_259[k]
                   + f_0 * kk_655[k];

        t_440[k] = -2.0 * hk_260[k]
                   + f_0 * kk_656[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, t_445, hk_261, hk_262, hk_263, hk_264, \
                         hk_265, kk_657, kk_658, kk_659, kk_660, \
                         kk_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = -2.0 * hk_261[k]
                   + f_0 * kk_657[k];

        t_442[k] = -2.0 * hk_262[k]
                   + f_0 * kk_658[k];

        t_443[k] = -2.0 * hk_263[k]
                   + f_0 * kk_659[k];

        t_444[k] = -2.0 * hk_264[k]
                   + f_0 * kk_660[k];

        t_445[k] = -2.0 * hk_265[k]
                   + f_0 * kk_661[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, t_450, hk_266, hk_267, hk_268, hk_269, \
                         hk_270, kk_662, kk_663, kk_664, kk_665, \
                         kk_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = -2.0 * hk_266[k]
                   + f_0 * kk_662[k];

        t_447[k] = -2.0 * hk_267[k]
                   + f_0 * kk_663[k];

        t_448[k] = -2.0 * hk_268[k]
                   + f_0 * kk_664[k];

        t_449[k] = -2.0 * hk_269[k]
                   + f_0 * kk_665[k];

        t_450[k] = -2.0 * hk_270[k]
                   + f_0 * kk_666[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, t_455, hk_271, hk_272, hk_273, hk_274, \
                         hk_275, kk_667, kk_668, kk_669, kk_670, \
                         kk_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = -2.0 * hk_271[k]
                   + f_0 * kk_667[k];

        t_452[k] = -2.0 * hk_272[k]
                   + f_0 * kk_668[k];

        t_453[k] = -2.0 * hk_273[k]
                   + f_0 * kk_669[k];

        t_454[k] = -2.0 * hk_274[k]
                   + f_0 * kk_670[k];

        t_455[k] = -2.0 * hk_275[k]
                   + f_0 * kk_671[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, t_460, hk_276, hk_277, hk_278, hk_279, \
                         hk_280, kk_672, kk_673, kk_674, kk_675, \
                         kk_676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = -2.0 * hk_276[k]
                   + f_0 * kk_672[k];

        t_457[k] = -2.0 * hk_277[k]
                   + f_0 * kk_673[k];

        t_458[k] = -2.0 * hk_278[k]
                   + f_0 * kk_674[k];

        t_459[k] = -2.0 * hk_279[k]
                   + f_0 * kk_675[k];

        t_460[k] = -2.0 * hk_280[k]
                   + f_0 * kk_676[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, t_465, hk_281, hk_282, hk_283, hk_284, \
                         hk_285, kk_677, kk_678, kk_679, kk_680, \
                         kk_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = -2.0 * hk_281[k]
                   + f_0 * kk_677[k];

        t_462[k] = -2.0 * hk_282[k]
                   + f_0 * kk_678[k];

        t_463[k] = -2.0 * hk_283[k]
                   + f_0 * kk_679[k];

        t_464[k] = -2.0 * hk_284[k]
                   + f_0 * kk_680[k];

        t_465[k] = -2.0 * hk_285[k]
                   + f_0 * kk_681[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, hk_286, hk_287, hk_288, hk_289, \
                         hk_290, kk_682, kk_683, kk_684, kk_685, \
                         kk_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = -2.0 * hk_286[k]
                   + f_0 * kk_682[k];

        t_467[k] = -2.0 * hk_287[k]
                   + f_0 * kk_683[k];

        t_468[k] = -3.0 * hk_288[k]
                   + f_0 * kk_684[k];

        t_469[k] = -3.0 * hk_289[k]
                   + f_0 * kk_685[k];

        t_470[k] = -3.0 * hk_290[k]
                   + f_0 * kk_686[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, t_475, hk_291, hk_292, hk_293, hk_294, \
                         hk_295, kk_687, kk_688, kk_689, kk_690, \
                         kk_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = -3.0 * hk_291[k]
                   + f_0 * kk_687[k];

        t_472[k] = -3.0 * hk_292[k]
                   + f_0 * kk_688[k];

        t_473[k] = -3.0 * hk_293[k]
                   + f_0 * kk_689[k];

        t_474[k] = -3.0 * hk_294[k]
                   + f_0 * kk_690[k];

        t_475[k] = -3.0 * hk_295[k]
                   + f_0 * kk_691[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, hk_296, hk_297, hk_298, hk_299, \
                         hk_300, kk_692, kk_693, kk_694, kk_695, \
                         kk_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = -3.0 * hk_296[k]
                   + f_0 * kk_692[k];

        t_477[k] = -3.0 * hk_297[k]
                   + f_0 * kk_693[k];

        t_478[k] = -3.0 * hk_298[k]
                   + f_0 * kk_694[k];

        t_479[k] = -3.0 * hk_299[k]
                   + f_0 * kk_695[k];

        t_480[k] = -3.0 * hk_300[k]
                   + f_0 * kk_696[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, hk_301, hk_302, hk_303, hk_304, \
                         hk_305, kk_697, kk_698, kk_699, kk_700, \
                         kk_701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = -3.0 * hk_301[k]
                   + f_0 * kk_697[k];

        t_482[k] = -3.0 * hk_302[k]
                   + f_0 * kk_698[k];

        t_483[k] = -3.0 * hk_303[k]
                   + f_0 * kk_699[k];

        t_484[k] = -3.0 * hk_304[k]
                   + f_0 * kk_700[k];

        t_485[k] = -3.0 * hk_305[k]
                   + f_0 * kk_701[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, hk_306, hk_307, hk_308, hk_309, \
                         hk_310, kk_702, kk_703, kk_704, kk_705, \
                         kk_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = -3.0 * hk_306[k]
                   + f_0 * kk_702[k];

        t_487[k] = -3.0 * hk_307[k]
                   + f_0 * kk_703[k];

        t_488[k] = -3.0 * hk_308[k]
                   + f_0 * kk_704[k];

        t_489[k] = -3.0 * hk_309[k]
                   + f_0 * kk_705[k];

        t_490[k] = -3.0 * hk_310[k]
                   + f_0 * kk_706[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, hk_311, hk_312, hk_313, hk_314, \
                         hk_315, kk_707, kk_708, kk_709, kk_710, \
                         kk_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = -3.0 * hk_311[k]
                   + f_0 * kk_707[k];

        t_492[k] = -3.0 * hk_312[k]
                   + f_0 * kk_708[k];

        t_493[k] = -3.0 * hk_313[k]
                   + f_0 * kk_709[k];

        t_494[k] = -3.0 * hk_314[k]
                   + f_0 * kk_710[k];

        t_495[k] = -3.0 * hk_315[k]
                   + f_0 * kk_711[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, t_500, hk_316, hk_317, hk_318, hk_319, \
                         hk_320, kk_712, kk_713, kk_714, kk_715, \
                         kk_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = -3.0 * hk_316[k]
                   + f_0 * kk_712[k];

        t_497[k] = -3.0 * hk_317[k]
                   + f_0 * kk_713[k];

        t_498[k] = -3.0 * hk_318[k]
                   + f_0 * kk_714[k];

        t_499[k] = -3.0 * hk_319[k]
                   + f_0 * kk_715[k];

        t_500[k] = -3.0 * hk_320[k]
                   + f_0 * kk_716[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, t_505, hk_321, hk_322, hk_323, hk_324, \
                         hk_325, kk_717, kk_718, kk_719, kk_720, \
                         kk_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = -3.0 * hk_321[k]
                   + f_0 * kk_717[k];

        t_502[k] = -3.0 * hk_322[k]
                   + f_0 * kk_718[k];

        t_503[k] = -3.0 * hk_323[k]
                   + f_0 * kk_719[k];

        t_504[k] = -4.0 * hk_324[k]
                   + f_0 * kk_720[k];

        t_505[k] = -4.0 * hk_325[k]
                   + f_0 * kk_721[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, t_509, t_510, hk_326, hk_327, hk_328, hk_329, \
                         hk_330, kk_722, kk_723, kk_724, kk_725, \
                         kk_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = -4.0 * hk_326[k]
                   + f_0 * kk_722[k];

        t_507[k] = -4.0 * hk_327[k]
                   + f_0 * kk_723[k];

        t_508[k] = -4.0 * hk_328[k]
                   + f_0 * kk_724[k];

        t_509[k] = -4.0 * hk_329[k]
                   + f_0 * kk_725[k];

        t_510[k] = -4.0 * hk_330[k]
                   + f_0 * kk_726[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_331 = buffer.data(hk + 331);
    const auto *hk_332 = buffer.data(hk + 332);
    const auto *hk_333 = buffer.data(hk + 333);
    const auto *hk_334 = buffer.data(hk + 334);
    const auto *hk_335 = buffer.data(hk + 335);
    const auto *hk_336 = buffer.data(hk + 336);
    const auto *hk_337 = buffer.data(hk + 337);
    const auto *hk_338 = buffer.data(hk + 338);
    const auto *hk_339 = buffer.data(hk + 339);
    const auto *hk_340 = buffer.data(hk + 340);
    const auto *hk_341 = buffer.data(hk + 341);
    const auto *hk_342 = buffer.data(hk + 342);
    const auto *hk_343 = buffer.data(hk + 343);
    const auto *hk_344 = buffer.data(hk + 344);
    const auto *hk_345 = buffer.data(hk + 345);
    const auto *hk_346 = buffer.data(hk + 346);
    const auto *hk_347 = buffer.data(hk + 347);
    const auto *hk_348 = buffer.data(hk + 348);
    const auto *hk_349 = buffer.data(hk + 349);
    const auto *hk_350 = buffer.data(hk + 350);
    const auto *hk_351 = buffer.data(hk + 351);
    const auto *hk_352 = buffer.data(hk + 352);
    const auto *hk_353 = buffer.data(hk + 353);
    const auto *hk_354 = buffer.data(hk + 354);
    const auto *hk_355 = buffer.data(hk + 355);
    const auto *hk_356 = buffer.data(hk + 356);
    const auto *hk_357 = buffer.data(hk + 357);
    const auto *hk_358 = buffer.data(hk + 358);
    const auto *hk_359 = buffer.data(hk + 359);
    const auto *hk_360 = buffer.data(hk + 360);
    const auto *hk_361 = buffer.data(hk + 361);
    const auto *hk_362 = buffer.data(hk + 362);
    const auto *hk_363 = buffer.data(hk + 363);
    const auto *hk_364 = buffer.data(hk + 364);
    const auto *hk_365 = buffer.data(hk + 365);
    const auto *hk_366 = buffer.data(hk + 366);
    const auto *hk_367 = buffer.data(hk + 367);
    const auto *hk_368 = buffer.data(hk + 368);
    const auto *hk_369 = buffer.data(hk + 369);
    const auto *hk_370 = buffer.data(hk + 370);
    const auto *hk_371 = buffer.data(hk + 371);
    const auto *hk_372 = buffer.data(hk + 372);
    const auto *hk_373 = buffer.data(hk + 373);
    const auto *hk_374 = buffer.data(hk + 374);
    const auto *hk_375 = buffer.data(hk + 375);
    const auto *hk_376 = buffer.data(hk + 376);
    const auto *hk_377 = buffer.data(hk + 377);
    const auto *hk_378 = buffer.data(hk + 378);
    const auto *hk_379 = buffer.data(hk + 379);
    const auto *hk_380 = buffer.data(hk + 380);
    const auto *hk_381 = buffer.data(hk + 381);
    const auto *hk_382 = buffer.data(hk + 382);
    const auto *hk_383 = buffer.data(hk + 383);
    const auto *hk_384 = buffer.data(hk + 384);
    const auto *hk_385 = buffer.data(hk + 385);
    const auto *hk_386 = buffer.data(hk + 386);
    const auto *hk_387 = buffer.data(hk + 387);
    const auto *hk_388 = buffer.data(hk + 388);
    const auto *hk_389 = buffer.data(hk + 389);
    const auto *hk_390 = buffer.data(hk + 390);
    const auto *hk_391 = buffer.data(hk + 391);
    const auto *hk_392 = buffer.data(hk + 392);
    const auto *hk_393 = buffer.data(hk + 393);
    const auto *hk_394 = buffer.data(hk + 394);
    const auto *hk_395 = buffer.data(hk + 395);
    const auto *hk_396 = buffer.data(hk + 396);
    const auto *hk_397 = buffer.data(hk + 397);
    const auto *hk_398 = buffer.data(hk + 398);
    const auto *hk_399 = buffer.data(hk + 399);
    const auto *hk_400 = buffer.data(hk + 400);
    const auto *hk_401 = buffer.data(hk + 401);
    const auto *hk_402 = buffer.data(hk + 402);
    const auto *hk_403 = buffer.data(hk + 403);
    const auto *hk_404 = buffer.data(hk + 404);
    const auto *hk_405 = buffer.data(hk + 405);
    const auto *hk_406 = buffer.data(hk + 406);
    const auto *hk_407 = buffer.data(hk + 407);
    const auto *hk_408 = buffer.data(hk + 408);
    const auto *hk_409 = buffer.data(hk + 409);
    const auto *hk_410 = buffer.data(hk + 410);
    const auto *hk_411 = buffer.data(hk + 411);
    const auto *hk_412 = buffer.data(hk + 412);
    const auto *hk_413 = buffer.data(hk + 413);
    const auto *hk_414 = buffer.data(hk + 414);
    const auto *hk_415 = buffer.data(hk + 415);
    const auto *hk_416 = buffer.data(hk + 416);
    const auto *hk_417 = buffer.data(hk + 417);
    const auto *hk_418 = buffer.data(hk + 418);
    const auto *hk_419 = buffer.data(hk + 419);
    const auto *hk_420 = buffer.data(hk + 420);
    const auto *hk_421 = buffer.data(hk + 421);
    const auto *hk_422 = buffer.data(hk + 422);
    const auto *hk_423 = buffer.data(hk + 423);
    const auto *hk_424 = buffer.data(hk + 424);
    const auto *hk_425 = buffer.data(hk + 425);
    const auto *hk_426 = buffer.data(hk + 426);
    const auto *hk_427 = buffer.data(hk + 427);
    const auto *hk_428 = buffer.data(hk + 428);
    const auto *hk_429 = buffer.data(hk + 429);
    const auto *hk_430 = buffer.data(hk + 430);
    const auto *hk_431 = buffer.data(hk + 431);
    const auto *hk_432 = buffer.data(hk + 432);
    const auto *hk_433 = buffer.data(hk + 433);
    const auto *hk_434 = buffer.data(hk + 434);
    const auto *hk_435 = buffer.data(hk + 435);
    const auto *hk_436 = buffer.data(hk + 436);
    const auto *hk_437 = buffer.data(hk + 437);
    const auto *hk_438 = buffer.data(hk + 438);
    const auto *hk_439 = buffer.data(hk + 439);
    const auto *hk_440 = buffer.data(hk + 440);
    const auto *hk_441 = buffer.data(hk + 441);
    const auto *hk_442 = buffer.data(hk + 442);
    const auto *hk_443 = buffer.data(hk + 443);
    const auto *hk_444 = buffer.data(hk + 444);
    const auto *hk_445 = buffer.data(hk + 445);
    const auto *hk_446 = buffer.data(hk + 446);
    const auto *hk_447 = buffer.data(hk + 447);
    const auto *hk_448 = buffer.data(hk + 448);
    const auto *hk_449 = buffer.data(hk + 449);
    const auto *hk_450 = buffer.data(hk + 450);
    const auto *hk_451 = buffer.data(hk + 451);
    const auto *hk_452 = buffer.data(hk + 452);
    const auto *hk_453 = buffer.data(hk + 453);

    const auto *kk_727 = buffer.data(kk + 727);
    const auto *kk_728 = buffer.data(kk + 728);
    const auto *kk_729 = buffer.data(kk + 729);
    const auto *kk_730 = buffer.data(kk + 730);
    const auto *kk_731 = buffer.data(kk + 731);
    const auto *kk_732 = buffer.data(kk + 732);
    const auto *kk_733 = buffer.data(kk + 733);
    const auto *kk_734 = buffer.data(kk + 734);
    const auto *kk_735 = buffer.data(kk + 735);
    const auto *kk_736 = buffer.data(kk + 736);
    const auto *kk_737 = buffer.data(kk + 737);
    const auto *kk_738 = buffer.data(kk + 738);
    const auto *kk_739 = buffer.data(kk + 739);
    const auto *kk_740 = buffer.data(kk + 740);
    const auto *kk_741 = buffer.data(kk + 741);
    const auto *kk_742 = buffer.data(kk + 742);
    const auto *kk_743 = buffer.data(kk + 743);
    const auto *kk_744 = buffer.data(kk + 744);
    const auto *kk_745 = buffer.data(kk + 745);
    const auto *kk_746 = buffer.data(kk + 746);
    const auto *kk_747 = buffer.data(kk + 747);
    const auto *kk_748 = buffer.data(kk + 748);
    const auto *kk_749 = buffer.data(kk + 749);
    const auto *kk_750 = buffer.data(kk + 750);
    const auto *kk_751 = buffer.data(kk + 751);
    const auto *kk_752 = buffer.data(kk + 752);
    const auto *kk_753 = buffer.data(kk + 753);
    const auto *kk_754 = buffer.data(kk + 754);
    const auto *kk_755 = buffer.data(kk + 755);
    const auto *kk_792 = buffer.data(kk + 792);
    const auto *kk_793 = buffer.data(kk + 793);
    const auto *kk_794 = buffer.data(kk + 794);
    const auto *kk_795 = buffer.data(kk + 795);
    const auto *kk_796 = buffer.data(kk + 796);
    const auto *kk_797 = buffer.data(kk + 797);
    const auto *kk_798 = buffer.data(kk + 798);
    const auto *kk_799 = buffer.data(kk + 799);
    const auto *kk_800 = buffer.data(kk + 800);
    const auto *kk_801 = buffer.data(kk + 801);
    const auto *kk_802 = buffer.data(kk + 802);
    const auto *kk_803 = buffer.data(kk + 803);
    const auto *kk_804 = buffer.data(kk + 804);
    const auto *kk_805 = buffer.data(kk + 805);
    const auto *kk_806 = buffer.data(kk + 806);
    const auto *kk_807 = buffer.data(kk + 807);
    const auto *kk_808 = buffer.data(kk + 808);
    const auto *kk_809 = buffer.data(kk + 809);
    const auto *kk_810 = buffer.data(kk + 810);
    const auto *kk_811 = buffer.data(kk + 811);
    const auto *kk_812 = buffer.data(kk + 812);
    const auto *kk_813 = buffer.data(kk + 813);
    const auto *kk_814 = buffer.data(kk + 814);
    const auto *kk_815 = buffer.data(kk + 815);
    const auto *kk_816 = buffer.data(kk + 816);
    const auto *kk_817 = buffer.data(kk + 817);
    const auto *kk_818 = buffer.data(kk + 818);
    const auto *kk_819 = buffer.data(kk + 819);
    const auto *kk_820 = buffer.data(kk + 820);
    const auto *kk_821 = buffer.data(kk + 821);
    const auto *kk_822 = buffer.data(kk + 822);
    const auto *kk_823 = buffer.data(kk + 823);
    const auto *kk_824 = buffer.data(kk + 824);
    const auto *kk_825 = buffer.data(kk + 825);
    const auto *kk_826 = buffer.data(kk + 826);
    const auto *kk_827 = buffer.data(kk + 827);
    const auto *kk_828 = buffer.data(kk + 828);
    const auto *kk_829 = buffer.data(kk + 829);
    const auto *kk_830 = buffer.data(kk + 830);
    const auto *kk_831 = buffer.data(kk + 831);
    const auto *kk_832 = buffer.data(kk + 832);
    const auto *kk_833 = buffer.data(kk + 833);
    const auto *kk_834 = buffer.data(kk + 834);
    const auto *kk_835 = buffer.data(kk + 835);
    const auto *kk_836 = buffer.data(kk + 836);
    const auto *kk_837 = buffer.data(kk + 837);
    const auto *kk_838 = buffer.data(kk + 838);
    const auto *kk_839 = buffer.data(kk + 839);
    const auto *kk_840 = buffer.data(kk + 840);
    const auto *kk_841 = buffer.data(kk + 841);
    const auto *kk_842 = buffer.data(kk + 842);
    const auto *kk_843 = buffer.data(kk + 843);
    const auto *kk_844 = buffer.data(kk + 844);
    const auto *kk_845 = buffer.data(kk + 845);
    const auto *kk_846 = buffer.data(kk + 846);
    const auto *kk_847 = buffer.data(kk + 847);
    const auto *kk_848 = buffer.data(kk + 848);
    const auto *kk_849 = buffer.data(kk + 849);
    const auto *kk_850 = buffer.data(kk + 850);
    const auto *kk_851 = buffer.data(kk + 851);
    const auto *kk_852 = buffer.data(kk + 852);
    const auto *kk_853 = buffer.data(kk + 853);
    const auto *kk_854 = buffer.data(kk + 854);
    const auto *kk_855 = buffer.data(kk + 855);
    const auto *kk_856 = buffer.data(kk + 856);
    const auto *kk_857 = buffer.data(kk + 857);
    const auto *kk_858 = buffer.data(kk + 858);
    const auto *kk_859 = buffer.data(kk + 859);
    const auto *kk_860 = buffer.data(kk + 860);
    const auto *kk_861 = buffer.data(kk + 861);
    const auto *kk_862 = buffer.data(kk + 862);
    const auto *kk_863 = buffer.data(kk + 863);
    const auto *kk_864 = buffer.data(kk + 864);
    const auto *kk_865 = buffer.data(kk + 865);
    const auto *kk_866 = buffer.data(kk + 866);
    const auto *kk_867 = buffer.data(kk + 867);
    const auto *kk_868 = buffer.data(kk + 868);
    const auto *kk_869 = buffer.data(kk + 869);
    const auto *kk_870 = buffer.data(kk + 870);
    const auto *kk_871 = buffer.data(kk + 871);
    const auto *kk_872 = buffer.data(kk + 872);
    const auto *kk_873 = buffer.data(kk + 873);
    const auto *kk_874 = buffer.data(kk + 874);
    const auto *kk_875 = buffer.data(kk + 875);
    const auto *kk_876 = buffer.data(kk + 876);
    const auto *kk_877 = buffer.data(kk + 877);
    const auto *kk_878 = buffer.data(kk + 878);
    const auto *kk_879 = buffer.data(kk + 879);
    const auto *kk_880 = buffer.data(kk + 880);
    const auto *kk_881 = buffer.data(kk + 881);
    const auto *kk_882 = buffer.data(kk + 882);
    const auto *kk_883 = buffer.data(kk + 883);
    const auto *kk_884 = buffer.data(kk + 884);
    const auto *kk_885 = buffer.data(kk + 885);
    const auto *kk_886 = buffer.data(kk + 886);
    const auto *kk_887 = buffer.data(kk + 887);
    const auto *kk_888 = buffer.data(kk + 888);
    const auto *kk_889 = buffer.data(kk + 889);
    const auto *kk_890 = buffer.data(kk + 890);
    const auto *kk_891 = buffer.data(kk + 891);
    const auto *kk_892 = buffer.data(kk + 892);
    const auto *kk_893 = buffer.data(kk + 893);
    const auto *kk_894 = buffer.data(kk + 894);
    const auto *kk_895 = buffer.data(kk + 895);
    const auto *kk_896 = buffer.data(kk + 896);
    const auto *kk_897 = buffer.data(kk + 897);
    const auto *kk_898 = buffer.data(kk + 898);
    const auto *kk_899 = buffer.data(kk + 899);
    const auto *kk_900 = buffer.data(kk + 900);
    const auto *kk_901 = buffer.data(kk + 901);
    const auto *kk_902 = buffer.data(kk + 902);
    const auto *kk_903 = buffer.data(kk + 903);
    const auto *kk_904 = buffer.data(kk + 904);
    const auto *kk_905 = buffer.data(kk + 905);
    const auto *kk_906 = buffer.data(kk + 906);
    const auto *kk_907 = buffer.data(kk + 907);
    const auto *kk_908 = buffer.data(kk + 908);
    const auto *kk_909 = buffer.data(kk + 909);
    const auto *kk_910 = buffer.data(kk + 910);
    const auto *kk_911 = buffer.data(kk + 911);
    const auto *kk_912 = buffer.data(kk + 912);
    const auto *kk_913 = buffer.data(kk + 913);
    const auto *kk_914 = buffer.data(kk + 914);
    const auto *kk_915 = buffer.data(kk + 915);
    const auto *kk_916 = buffer.data(kk + 916);
    const auto *kk_917 = buffer.data(kk + 917);
    const auto *kk_918 = buffer.data(kk + 918);
    const auto *kk_919 = buffer.data(kk + 919);
    const auto *kk_920 = buffer.data(kk + 920);
    const auto *kk_921 = buffer.data(kk + 921);

#pragma omp simd aligned(t_511, t_512, t_513, t_514, t_515, hk_331, hk_332, hk_333, hk_334, \
                         hk_335, kk_727, kk_728, kk_729, kk_730, \
                         kk_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = -4.0 * hk_331[k]
                   + f_0 * kk_727[k];

        t_512[k] = -4.0 * hk_332[k]
                   + f_0 * kk_728[k];

        t_513[k] = -4.0 * hk_333[k]
                   + f_0 * kk_729[k];

        t_514[k] = -4.0 * hk_334[k]
                   + f_0 * kk_730[k];

        t_515[k] = -4.0 * hk_335[k]
                   + f_0 * kk_731[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, t_520, hk_336, hk_337, hk_338, hk_339, \
                         hk_340, kk_732, kk_733, kk_734, kk_735, \
                         kk_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = -4.0 * hk_336[k]
                   + f_0 * kk_732[k];

        t_517[k] = -4.0 * hk_337[k]
                   + f_0 * kk_733[k];

        t_518[k] = -4.0 * hk_338[k]
                   + f_0 * kk_734[k];

        t_519[k] = -4.0 * hk_339[k]
                   + f_0 * kk_735[k];

        t_520[k] = -4.0 * hk_340[k]
                   + f_0 * kk_736[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, t_525, hk_341, hk_342, hk_343, hk_344, \
                         hk_345, kk_737, kk_738, kk_739, kk_740, \
                         kk_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = -4.0 * hk_341[k]
                   + f_0 * kk_737[k];

        t_522[k] = -4.0 * hk_342[k]
                   + f_0 * kk_738[k];

        t_523[k] = -4.0 * hk_343[k]
                   + f_0 * kk_739[k];

        t_524[k] = -4.0 * hk_344[k]
                   + f_0 * kk_740[k];

        t_525[k] = -4.0 * hk_345[k]
                   + f_0 * kk_741[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, hk_346, hk_347, hk_348, hk_349, \
                         hk_350, kk_742, kk_743, kk_744, kk_745, \
                         kk_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = -4.0 * hk_346[k]
                   + f_0 * kk_742[k];

        t_527[k] = -4.0 * hk_347[k]
                   + f_0 * kk_743[k];

        t_528[k] = -4.0 * hk_348[k]
                   + f_0 * kk_744[k];

        t_529[k] = -4.0 * hk_349[k]
                   + f_0 * kk_745[k];

        t_530[k] = -4.0 * hk_350[k]
                   + f_0 * kk_746[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, t_535, hk_351, hk_352, hk_353, hk_354, \
                         hk_355, kk_747, kk_748, kk_749, kk_750, \
                         kk_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = -4.0 * hk_351[k]
                   + f_0 * kk_747[k];

        t_532[k] = -4.0 * hk_352[k]
                   + f_0 * kk_748[k];

        t_533[k] = -4.0 * hk_353[k]
                   + f_0 * kk_749[k];

        t_534[k] = -4.0 * hk_354[k]
                   + f_0 * kk_750[k];

        t_535[k] = -4.0 * hk_355[k]
                   + f_0 * kk_751[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, t_541, hk_356, hk_357, hk_358, \
                         hk_359, kk_752, kk_753, kk_754, kk_755, kk_792, \
                         kk_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = -4.0 * hk_356[k]
                   + f_0 * kk_752[k];

        t_537[k] = -4.0 * hk_357[k]
                   + f_0 * kk_753[k];

        t_538[k] = -4.0 * hk_358[k]
                   + f_0 * kk_754[k];

        t_539[k] = -4.0 * hk_359[k]
                   + f_0 * kk_755[k];

        t_540[k] = f_0 * kk_792[k];

        t_541[k] = f_0 * kk_793[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, t_547, t_548, t_549, kk_794, \
                         kk_795, kk_796, kk_797, kk_798, kk_799, kk_800, \
                         kk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_0 * kk_794[k];

        t_543[k] = f_0 * kk_795[k];

        t_544[k] = f_0 * kk_796[k];

        t_545[k] = f_0 * kk_797[k];

        t_546[k] = f_0 * kk_798[k];

        t_547[k] = f_0 * kk_799[k];

        t_548[k] = f_0 * kk_800[k];

        t_549[k] = f_0 * kk_801[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, t_555, t_556, t_557, kk_802, \
                         kk_803, kk_804, kk_805, kk_806, kk_807, kk_808, \
                         kk_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_0 * kk_802[k];

        t_551[k] = f_0 * kk_803[k];

        t_552[k] = f_0 * kk_804[k];

        t_553[k] = f_0 * kk_805[k];

        t_554[k] = f_0 * kk_806[k];

        t_555[k] = f_0 * kk_807[k];

        t_556[k] = f_0 * kk_808[k];

        t_557[k] = f_0 * kk_809[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, t_562, t_563, t_564, t_565, kk_810, \
                         kk_811, kk_812, kk_813, kk_814, kk_815, kk_816, \
                         kk_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_0 * kk_810[k];

        t_559[k] = f_0 * kk_811[k];

        t_560[k] = f_0 * kk_812[k];

        t_561[k] = f_0 * kk_813[k];

        t_562[k] = f_0 * kk_814[k];

        t_563[k] = f_0 * kk_815[k];

        t_564[k] = f_0 * kk_816[k];

        t_565[k] = f_0 * kk_817[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, t_570, t_571, t_572, t_573, kk_818, \
                         kk_819, kk_820, kk_821, kk_822, kk_823, kk_824, \
                         kk_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_0 * kk_818[k];

        t_567[k] = f_0 * kk_819[k];

        t_568[k] = f_0 * kk_820[k];

        t_569[k] = f_0 * kk_821[k];

        t_570[k] = f_0 * kk_822[k];

        t_571[k] = f_0 * kk_823[k];

        t_572[k] = f_0 * kk_824[k];

        t_573[k] = f_0 * kk_825[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, t_579, hk_360, hk_361, hk_362, \
                         hk_363, kk_826, kk_827, kk_828, kk_829, kk_830, \
                         kk_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_0 * kk_826[k];

        t_575[k] = f_0 * kk_827[k];

        t_576[k] = -hk_360[k]
                   + f_0 * kk_828[k];

        t_577[k] = -hk_361[k]
                   + f_0 * kk_829[k];

        t_578[k] = -hk_362[k]
                   + f_0 * kk_830[k];

        t_579[k] = -hk_363[k]
                   + f_0 * kk_831[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, hk_364, hk_365, hk_366, hk_367, \
                         hk_368, kk_832, kk_833, kk_834, kk_835, \
                         kk_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -hk_364[k]
                   + f_0 * kk_832[k];

        t_581[k] = -hk_365[k]
                   + f_0 * kk_833[k];

        t_582[k] = -hk_366[k]
                   + f_0 * kk_834[k];

        t_583[k] = -hk_367[k]
                   + f_0 * kk_835[k];

        t_584[k] = -hk_368[k]
                   + f_0 * kk_836[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, hk_369, hk_370, hk_371, hk_372, \
                         hk_373, kk_837, kk_838, kk_839, kk_840, \
                         kk_841 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -hk_369[k]
                   + f_0 * kk_837[k];

        t_586[k] = -hk_370[k]
                   + f_0 * kk_838[k];

        t_587[k] = -hk_371[k]
                   + f_0 * kk_839[k];

        t_588[k] = -hk_372[k]
                   + f_0 * kk_840[k];

        t_589[k] = -hk_373[k]
                   + f_0 * kk_841[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, hk_374, hk_375, hk_376, hk_377, \
                         hk_378, kk_842, kk_843, kk_844, kk_845, \
                         kk_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -hk_374[k]
                   + f_0 * kk_842[k];

        t_591[k] = -hk_375[k]
                   + f_0 * kk_843[k];

        t_592[k] = -hk_376[k]
                   + f_0 * kk_844[k];

        t_593[k] = -hk_377[k]
                   + f_0 * kk_845[k];

        t_594[k] = -hk_378[k]
                   + f_0 * kk_846[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, hk_379, hk_380, hk_381, hk_382, \
                         hk_383, kk_847, kk_848, kk_849, kk_850, \
                         kk_851 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -hk_379[k]
                   + f_0 * kk_847[k];

        t_596[k] = -hk_380[k]
                   + f_0 * kk_848[k];

        t_597[k] = -hk_381[k]
                   + f_0 * kk_849[k];

        t_598[k] = -hk_382[k]
                   + f_0 * kk_850[k];

        t_599[k] = -hk_383[k]
                   + f_0 * kk_851[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, hk_384, hk_385, hk_386, hk_387, \
                         hk_388, kk_852, kk_853, kk_854, kk_855, \
                         kk_856 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -hk_384[k]
                   + f_0 * kk_852[k];

        t_601[k] = -hk_385[k]
                   + f_0 * kk_853[k];

        t_602[k] = -hk_386[k]
                   + f_0 * kk_854[k];

        t_603[k] = -hk_387[k]
                   + f_0 * kk_855[k];

        t_604[k] = -hk_388[k]
                   + f_0 * kk_856[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, hk_389, hk_390, hk_391, hk_392, \
                         hk_393, kk_857, kk_858, kk_859, kk_860, \
                         kk_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -hk_389[k]
                   + f_0 * kk_857[k];

        t_606[k] = -hk_390[k]
                   + f_0 * kk_858[k];

        t_607[k] = -hk_391[k]
                   + f_0 * kk_859[k];

        t_608[k] = -hk_392[k]
                   + f_0 * kk_860[k];

        t_609[k] = -hk_393[k]
                   + f_0 * kk_861[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, hk_394, hk_395, hk_396, hk_397, \
                         hk_398, kk_862, kk_863, kk_864, kk_865, \
                         kk_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -hk_394[k]
                   + f_0 * kk_862[k];

        t_611[k] = -hk_395[k]
                   + f_0 * kk_863[k];

        t_612[k] = -2.0 * hk_396[k]
                   + f_0 * kk_864[k];

        t_613[k] = -2.0 * hk_397[k]
                   + f_0 * kk_865[k];

        t_614[k] = -2.0 * hk_398[k]
                   + f_0 * kk_866[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, hk_399, hk_400, hk_401, hk_402, \
                         hk_403, kk_867, kk_868, kk_869, kk_870, \
                         kk_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -2.0 * hk_399[k]
                   + f_0 * kk_867[k];

        t_616[k] = -2.0 * hk_400[k]
                   + f_0 * kk_868[k];

        t_617[k] = -2.0 * hk_401[k]
                   + f_0 * kk_869[k];

        t_618[k] = -2.0 * hk_402[k]
                   + f_0 * kk_870[k];

        t_619[k] = -2.0 * hk_403[k]
                   + f_0 * kk_871[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, hk_404, hk_405, hk_406, hk_407, \
                         hk_408, kk_872, kk_873, kk_874, kk_875, \
                         kk_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -2.0 * hk_404[k]
                   + f_0 * kk_872[k];

        t_621[k] = -2.0 * hk_405[k]
                   + f_0 * kk_873[k];

        t_622[k] = -2.0 * hk_406[k]
                   + f_0 * kk_874[k];

        t_623[k] = -2.0 * hk_407[k]
                   + f_0 * kk_875[k];

        t_624[k] = -2.0 * hk_408[k]
                   + f_0 * kk_876[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, hk_409, hk_410, hk_411, hk_412, \
                         hk_413, kk_877, kk_878, kk_879, kk_880, \
                         kk_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -2.0 * hk_409[k]
                   + f_0 * kk_877[k];

        t_626[k] = -2.0 * hk_410[k]
                   + f_0 * kk_878[k];

        t_627[k] = -2.0 * hk_411[k]
                   + f_0 * kk_879[k];

        t_628[k] = -2.0 * hk_412[k]
                   + f_0 * kk_880[k];

        t_629[k] = -2.0 * hk_413[k]
                   + f_0 * kk_881[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, hk_414, hk_415, hk_416, hk_417, \
                         hk_418, kk_882, kk_883, kk_884, kk_885, \
                         kk_886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -2.0 * hk_414[k]
                   + f_0 * kk_882[k];

        t_631[k] = -2.0 * hk_415[k]
                   + f_0 * kk_883[k];

        t_632[k] = -2.0 * hk_416[k]
                   + f_0 * kk_884[k];

        t_633[k] = -2.0 * hk_417[k]
                   + f_0 * kk_885[k];

        t_634[k] = -2.0 * hk_418[k]
                   + f_0 * kk_886[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, hk_419, hk_420, hk_421, hk_422, \
                         hk_423, kk_887, kk_888, kk_889, kk_890, \
                         kk_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -2.0 * hk_419[k]
                   + f_0 * kk_887[k];

        t_636[k] = -2.0 * hk_420[k]
                   + f_0 * kk_888[k];

        t_637[k] = -2.0 * hk_421[k]
                   + f_0 * kk_889[k];

        t_638[k] = -2.0 * hk_422[k]
                   + f_0 * kk_890[k];

        t_639[k] = -2.0 * hk_423[k]
                   + f_0 * kk_891[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, hk_424, hk_425, hk_426, hk_427, \
                         hk_428, kk_892, kk_893, kk_894, kk_895, \
                         kk_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -2.0 * hk_424[k]
                   + f_0 * kk_892[k];

        t_641[k] = -2.0 * hk_425[k]
                   + f_0 * kk_893[k];

        t_642[k] = -2.0 * hk_426[k]
                   + f_0 * kk_894[k];

        t_643[k] = -2.0 * hk_427[k]
                   + f_0 * kk_895[k];

        t_644[k] = -2.0 * hk_428[k]
                   + f_0 * kk_896[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, hk_429, hk_430, hk_431, hk_432, \
                         hk_433, kk_897, kk_898, kk_899, kk_900, \
                         kk_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -2.0 * hk_429[k]
                   + f_0 * kk_897[k];

        t_646[k] = -2.0 * hk_430[k]
                   + f_0 * kk_898[k];

        t_647[k] = -2.0 * hk_431[k]
                   + f_0 * kk_899[k];

        t_648[k] = -3.0 * hk_432[k]
                   + f_0 * kk_900[k];

        t_649[k] = -3.0 * hk_433[k]
                   + f_0 * kk_901[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, hk_434, hk_435, hk_436, hk_437, \
                         hk_438, kk_902, kk_903, kk_904, kk_905, \
                         kk_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -3.0 * hk_434[k]
                   + f_0 * kk_902[k];

        t_651[k] = -3.0 * hk_435[k]
                   + f_0 * kk_903[k];

        t_652[k] = -3.0 * hk_436[k]
                   + f_0 * kk_904[k];

        t_653[k] = -3.0 * hk_437[k]
                   + f_0 * kk_905[k];

        t_654[k] = -3.0 * hk_438[k]
                   + f_0 * kk_906[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, hk_439, hk_440, hk_441, hk_442, \
                         hk_443, kk_907, kk_908, kk_909, kk_910, \
                         kk_911 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -3.0 * hk_439[k]
                   + f_0 * kk_907[k];

        t_656[k] = -3.0 * hk_440[k]
                   + f_0 * kk_908[k];

        t_657[k] = -3.0 * hk_441[k]
                   + f_0 * kk_909[k];

        t_658[k] = -3.0 * hk_442[k]
                   + f_0 * kk_910[k];

        t_659[k] = -3.0 * hk_443[k]
                   + f_0 * kk_911[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, hk_444, hk_445, hk_446, hk_447, \
                         hk_448, kk_912, kk_913, kk_914, kk_915, \
                         kk_916 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -3.0 * hk_444[k]
                   + f_0 * kk_912[k];

        t_661[k] = -3.0 * hk_445[k]
                   + f_0 * kk_913[k];

        t_662[k] = -3.0 * hk_446[k]
                   + f_0 * kk_914[k];

        t_663[k] = -3.0 * hk_447[k]
                   + f_0 * kk_915[k];

        t_664[k] = -3.0 * hk_448[k]
                   + f_0 * kk_916[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, hk_449, hk_450, hk_451, hk_452, \
                         hk_453, kk_917, kk_918, kk_919, kk_920, \
                         kk_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -3.0 * hk_449[k]
                   + f_0 * kk_917[k];

        t_666[k] = -3.0 * hk_450[k]
                   + f_0 * kk_918[k];

        t_667[k] = -3.0 * hk_451[k]
                   + f_0 * kk_919[k];

        t_668[k] = -3.0 * hk_452[k]
                   + f_0 * kk_920[k];

        t_669[k] = -3.0 * hk_453[k]
                   + f_0 * kk_921[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_454 = buffer.data(hk + 454);
    const auto *hk_455 = buffer.data(hk + 455);
    const auto *hk_456 = buffer.data(hk + 456);
    const auto *hk_457 = buffer.data(hk + 457);
    const auto *hk_458 = buffer.data(hk + 458);
    const auto *hk_459 = buffer.data(hk + 459);
    const auto *hk_460 = buffer.data(hk + 460);
    const auto *hk_461 = buffer.data(hk + 461);
    const auto *hk_462 = buffer.data(hk + 462);
    const auto *hk_463 = buffer.data(hk + 463);
    const auto *hk_464 = buffer.data(hk + 464);
    const auto *hk_465 = buffer.data(hk + 465);
    const auto *hk_466 = buffer.data(hk + 466);
    const auto *hk_467 = buffer.data(hk + 467);
    const auto *hk_468 = buffer.data(hk + 468);
    const auto *hk_469 = buffer.data(hk + 469);
    const auto *hk_470 = buffer.data(hk + 470);
    const auto *hk_471 = buffer.data(hk + 471);
    const auto *hk_472 = buffer.data(hk + 472);
    const auto *hk_473 = buffer.data(hk + 473);
    const auto *hk_474 = buffer.data(hk + 474);
    const auto *hk_475 = buffer.data(hk + 475);
    const auto *hk_476 = buffer.data(hk + 476);
    const auto *hk_477 = buffer.data(hk + 477);
    const auto *hk_478 = buffer.data(hk + 478);
    const auto *hk_479 = buffer.data(hk + 479);
    const auto *hk_480 = buffer.data(hk + 480);
    const auto *hk_481 = buffer.data(hk + 481);
    const auto *hk_482 = buffer.data(hk + 482);
    const auto *hk_483 = buffer.data(hk + 483);
    const auto *hk_484 = buffer.data(hk + 484);
    const auto *hk_485 = buffer.data(hk + 485);
    const auto *hk_486 = buffer.data(hk + 486);
    const auto *hk_487 = buffer.data(hk + 487);
    const auto *hk_488 = buffer.data(hk + 488);
    const auto *hk_489 = buffer.data(hk + 489);
    const auto *hk_490 = buffer.data(hk + 490);
    const auto *hk_491 = buffer.data(hk + 491);
    const auto *hk_492 = buffer.data(hk + 492);
    const auto *hk_493 = buffer.data(hk + 493);
    const auto *hk_494 = buffer.data(hk + 494);
    const auto *hk_495 = buffer.data(hk + 495);
    const auto *hk_496 = buffer.data(hk + 496);
    const auto *hk_497 = buffer.data(hk + 497);
    const auto *hk_498 = buffer.data(hk + 498);
    const auto *hk_499 = buffer.data(hk + 499);
    const auto *hk_500 = buffer.data(hk + 500);
    const auto *hk_501 = buffer.data(hk + 501);
    const auto *hk_502 = buffer.data(hk + 502);
    const auto *hk_503 = buffer.data(hk + 503);
    const auto *hk_504 = buffer.data(hk + 504);
    const auto *hk_505 = buffer.data(hk + 505);
    const auto *hk_506 = buffer.data(hk + 506);
    const auto *hk_507 = buffer.data(hk + 507);
    const auto *hk_508 = buffer.data(hk + 508);
    const auto *hk_509 = buffer.data(hk + 509);
    const auto *hk_510 = buffer.data(hk + 510);
    const auto *hk_511 = buffer.data(hk + 511);
    const auto *hk_512 = buffer.data(hk + 512);
    const auto *hk_513 = buffer.data(hk + 513);
    const auto *hk_514 = buffer.data(hk + 514);
    const auto *hk_515 = buffer.data(hk + 515);
    const auto *hk_516 = buffer.data(hk + 516);
    const auto *hk_517 = buffer.data(hk + 517);
    const auto *hk_518 = buffer.data(hk + 518);
    const auto *hk_519 = buffer.data(hk + 519);
    const auto *hk_520 = buffer.data(hk + 520);
    const auto *hk_521 = buffer.data(hk + 521);
    const auto *hk_522 = buffer.data(hk + 522);
    const auto *hk_523 = buffer.data(hk + 523);
    const auto *hk_524 = buffer.data(hk + 524);
    const auto *hk_525 = buffer.data(hk + 525);
    const auto *hk_526 = buffer.data(hk + 526);
    const auto *hk_527 = buffer.data(hk + 527);
    const auto *hk_528 = buffer.data(hk + 528);
    const auto *hk_529 = buffer.data(hk + 529);
    const auto *hk_530 = buffer.data(hk + 530);
    const auto *hk_531 = buffer.data(hk + 531);
    const auto *hk_532 = buffer.data(hk + 532);
    const auto *hk_533 = buffer.data(hk + 533);
    const auto *hk_534 = buffer.data(hk + 534);
    const auto *hk_535 = buffer.data(hk + 535);
    const auto *hk_536 = buffer.data(hk + 536);
    const auto *hk_537 = buffer.data(hk + 537);
    const auto *hk_538 = buffer.data(hk + 538);
    const auto *hk_539 = buffer.data(hk + 539);
    const auto *hk_540 = buffer.data(hk + 540);
    const auto *hk_541 = buffer.data(hk + 541);
    const auto *hk_542 = buffer.data(hk + 542);
    const auto *hk_543 = buffer.data(hk + 543);
    const auto *hk_544 = buffer.data(hk + 544);
    const auto *hk_545 = buffer.data(hk + 545);
    const auto *hk_546 = buffer.data(hk + 546);
    const auto *hk_547 = buffer.data(hk + 547);
    const auto *hk_548 = buffer.data(hk + 548);
    const auto *hk_549 = buffer.data(hk + 549);
    const auto *hk_550 = buffer.data(hk + 550);
    const auto *hk_551 = buffer.data(hk + 551);
    const auto *hk_552 = buffer.data(hk + 552);
    const auto *hk_553 = buffer.data(hk + 553);
    const auto *hk_554 = buffer.data(hk + 554);
    const auto *hk_555 = buffer.data(hk + 555);
    const auto *hk_556 = buffer.data(hk + 556);
    const auto *hk_557 = buffer.data(hk + 557);
    const auto *hk_558 = buffer.data(hk + 558);
    const auto *hk_559 = buffer.data(hk + 559);
    const auto *hk_560 = buffer.data(hk + 560);
    const auto *hk_561 = buffer.data(hk + 561);
    const auto *hk_562 = buffer.data(hk + 562);
    const auto *hk_563 = buffer.data(hk + 563);
    const auto *hk_564 = buffer.data(hk + 564);
    const auto *hk_565 = buffer.data(hk + 565);
    const auto *hk_566 = buffer.data(hk + 566);
    const auto *hk_567 = buffer.data(hk + 567);
    const auto *hk_568 = buffer.data(hk + 568);
    const auto *hk_569 = buffer.data(hk + 569);
    const auto *hk_570 = buffer.data(hk + 570);
    const auto *hk_571 = buffer.data(hk + 571);
    const auto *hk_572 = buffer.data(hk + 572);
    const auto *hk_573 = buffer.data(hk + 573);
    const auto *hk_574 = buffer.data(hk + 574);
    const auto *hk_575 = buffer.data(hk + 575);
    const auto *hk_576 = buffer.data(hk + 576);
    const auto *hk_577 = buffer.data(hk + 577);
    const auto *hk_578 = buffer.data(hk + 578);
    const auto *hk_579 = buffer.data(hk + 579);
    const auto *hk_580 = buffer.data(hk + 580);

    const auto *kk_922 = buffer.data(kk + 922);
    const auto *kk_923 = buffer.data(kk + 923);
    const auto *kk_924 = buffer.data(kk + 924);
    const auto *kk_925 = buffer.data(kk + 925);
    const auto *kk_926 = buffer.data(kk + 926);
    const auto *kk_927 = buffer.data(kk + 927);
    const auto *kk_928 = buffer.data(kk + 928);
    const auto *kk_929 = buffer.data(kk + 929);
    const auto *kk_930 = buffer.data(kk + 930);
    const auto *kk_931 = buffer.data(kk + 931);
    const auto *kk_932 = buffer.data(kk + 932);
    const auto *kk_933 = buffer.data(kk + 933);
    const auto *kk_934 = buffer.data(kk + 934);
    const auto *kk_935 = buffer.data(kk + 935);
    const auto *kk_936 = buffer.data(kk + 936);
    const auto *kk_937 = buffer.data(kk + 937);
    const auto *kk_938 = buffer.data(kk + 938);
    const auto *kk_939 = buffer.data(kk + 939);
    const auto *kk_940 = buffer.data(kk + 940);
    const auto *kk_941 = buffer.data(kk + 941);
    const auto *kk_942 = buffer.data(kk + 942);
    const auto *kk_943 = buffer.data(kk + 943);
    const auto *kk_944 = buffer.data(kk + 944);
    const auto *kk_945 = buffer.data(kk + 945);
    const auto *kk_946 = buffer.data(kk + 946);
    const auto *kk_947 = buffer.data(kk + 947);
    const auto *kk_948 = buffer.data(kk + 948);
    const auto *kk_949 = buffer.data(kk + 949);
    const auto *kk_950 = buffer.data(kk + 950);
    const auto *kk_951 = buffer.data(kk + 951);
    const auto *kk_952 = buffer.data(kk + 952);
    const auto *kk_953 = buffer.data(kk + 953);
    const auto *kk_954 = buffer.data(kk + 954);
    const auto *kk_955 = buffer.data(kk + 955);
    const auto *kk_956 = buffer.data(kk + 956);
    const auto *kk_957 = buffer.data(kk + 957);
    const auto *kk_958 = buffer.data(kk + 958);
    const auto *kk_959 = buffer.data(kk + 959);
    const auto *kk_960 = buffer.data(kk + 960);
    const auto *kk_961 = buffer.data(kk + 961);
    const auto *kk_962 = buffer.data(kk + 962);
    const auto *kk_963 = buffer.data(kk + 963);
    const auto *kk_964 = buffer.data(kk + 964);
    const auto *kk_965 = buffer.data(kk + 965);
    const auto *kk_966 = buffer.data(kk + 966);
    const auto *kk_967 = buffer.data(kk + 967);
    const auto *kk_968 = buffer.data(kk + 968);
    const auto *kk_969 = buffer.data(kk + 969);
    const auto *kk_970 = buffer.data(kk + 970);
    const auto *kk_971 = buffer.data(kk + 971);
    const auto *kk_972 = buffer.data(kk + 972);
    const auto *kk_973 = buffer.data(kk + 973);
    const auto *kk_974 = buffer.data(kk + 974);
    const auto *kk_975 = buffer.data(kk + 975);
    const auto *kk_976 = buffer.data(kk + 976);
    const auto *kk_977 = buffer.data(kk + 977);
    const auto *kk_978 = buffer.data(kk + 978);
    const auto *kk_979 = buffer.data(kk + 979);
    const auto *kk_980 = buffer.data(kk + 980);
    const auto *kk_981 = buffer.data(kk + 981);
    const auto *kk_982 = buffer.data(kk + 982);
    const auto *kk_983 = buffer.data(kk + 983);
    const auto *kk_984 = buffer.data(kk + 984);
    const auto *kk_985 = buffer.data(kk + 985);
    const auto *kk_986 = buffer.data(kk + 986);
    const auto *kk_987 = buffer.data(kk + 987);
    const auto *kk_988 = buffer.data(kk + 988);
    const auto *kk_989 = buffer.data(kk + 989);
    const auto *kk_990 = buffer.data(kk + 990);
    const auto *kk_991 = buffer.data(kk + 991);
    const auto *kk_992 = buffer.data(kk + 992);
    const auto *kk_993 = buffer.data(kk + 993);
    const auto *kk_994 = buffer.data(kk + 994);
    const auto *kk_995 = buffer.data(kk + 995);
    const auto *kk_996 = buffer.data(kk + 996);
    const auto *kk_997 = buffer.data(kk + 997);
    const auto *kk_998 = buffer.data(kk + 998);
    const auto *kk_999 = buffer.data(kk + 999);
    const auto *kk_1000 = buffer.data(kk + 1000);
    const auto *kk_1001 = buffer.data(kk + 1001);
    const auto *kk_1002 = buffer.data(kk + 1002);
    const auto *kk_1003 = buffer.data(kk + 1003);
    const auto *kk_1004 = buffer.data(kk + 1004);
    const auto *kk_1005 = buffer.data(kk + 1005);
    const auto *kk_1006 = buffer.data(kk + 1006);
    const auto *kk_1007 = buffer.data(kk + 1007);
    const auto *kk_1044 = buffer.data(kk + 1044);
    const auto *kk_1045 = buffer.data(kk + 1045);
    const auto *kk_1046 = buffer.data(kk + 1046);
    const auto *kk_1047 = buffer.data(kk + 1047);
    const auto *kk_1048 = buffer.data(kk + 1048);
    const auto *kk_1049 = buffer.data(kk + 1049);
    const auto *kk_1050 = buffer.data(kk + 1050);
    const auto *kk_1051 = buffer.data(kk + 1051);
    const auto *kk_1052 = buffer.data(kk + 1052);
    const auto *kk_1053 = buffer.data(kk + 1053);
    const auto *kk_1054 = buffer.data(kk + 1054);
    const auto *kk_1055 = buffer.data(kk + 1055);
    const auto *kk_1056 = buffer.data(kk + 1056);
    const auto *kk_1057 = buffer.data(kk + 1057);
    const auto *kk_1058 = buffer.data(kk + 1058);
    const auto *kk_1059 = buffer.data(kk + 1059);
    const auto *kk_1060 = buffer.data(kk + 1060);
    const auto *kk_1061 = buffer.data(kk + 1061);
    const auto *kk_1062 = buffer.data(kk + 1062);
    const auto *kk_1063 = buffer.data(kk + 1063);
    const auto *kk_1064 = buffer.data(kk + 1064);
    const auto *kk_1065 = buffer.data(kk + 1065);
    const auto *kk_1066 = buffer.data(kk + 1066);
    const auto *kk_1067 = buffer.data(kk + 1067);
    const auto *kk_1068 = buffer.data(kk + 1068);
    const auto *kk_1069 = buffer.data(kk + 1069);
    const auto *kk_1070 = buffer.data(kk + 1070);
    const auto *kk_1071 = buffer.data(kk + 1071);
    const auto *kk_1072 = buffer.data(kk + 1072);
    const auto *kk_1073 = buffer.data(kk + 1073);
    const auto *kk_1074 = buffer.data(kk + 1074);
    const auto *kk_1075 = buffer.data(kk + 1075);
    const auto *kk_1076 = buffer.data(kk + 1076);
    const auto *kk_1077 = buffer.data(kk + 1077);
    const auto *kk_1078 = buffer.data(kk + 1078);
    const auto *kk_1079 = buffer.data(kk + 1079);
    const auto *kk_1080 = buffer.data(kk + 1080);
    const auto *kk_1081 = buffer.data(kk + 1081);
    const auto *kk_1082 = buffer.data(kk + 1082);
    const auto *kk_1083 = buffer.data(kk + 1083);
    const auto *kk_1084 = buffer.data(kk + 1084);
    const auto *kk_1085 = buffer.data(kk + 1085);
    const auto *kk_1086 = buffer.data(kk + 1086);
    const auto *kk_1087 = buffer.data(kk + 1087);
    const auto *kk_1088 = buffer.data(kk + 1088);
    const auto *kk_1089 = buffer.data(kk + 1089);
    const auto *kk_1090 = buffer.data(kk + 1090);
    const auto *kk_1091 = buffer.data(kk + 1091);
    const auto *kk_1092 = buffer.data(kk + 1092);
    const auto *kk_1093 = buffer.data(kk + 1093);
    const auto *kk_1094 = buffer.data(kk + 1094);
    const auto *kk_1095 = buffer.data(kk + 1095);
    const auto *kk_1096 = buffer.data(kk + 1096);
    const auto *kk_1097 = buffer.data(kk + 1097);
    const auto *kk_1098 = buffer.data(kk + 1098);
    const auto *kk_1099 = buffer.data(kk + 1099);
    const auto *kk_1100 = buffer.data(kk + 1100);
    const auto *kk_1101 = buffer.data(kk + 1101);
    const auto *kk_1102 = buffer.data(kk + 1102);
    const auto *kk_1103 = buffer.data(kk + 1103);
    const auto *kk_1104 = buffer.data(kk + 1104);
    const auto *kk_1105 = buffer.data(kk + 1105);
    const auto *kk_1106 = buffer.data(kk + 1106);
    const auto *kk_1107 = buffer.data(kk + 1107);
    const auto *kk_1108 = buffer.data(kk + 1108);
    const auto *kk_1109 = buffer.data(kk + 1109);
    const auto *kk_1110 = buffer.data(kk + 1110);
    const auto *kk_1111 = buffer.data(kk + 1111);
    const auto *kk_1112 = buffer.data(kk + 1112);
    const auto *kk_1113 = buffer.data(kk + 1113);
    const auto *kk_1114 = buffer.data(kk + 1114);
    const auto *kk_1115 = buffer.data(kk + 1115);
    const auto *kk_1116 = buffer.data(kk + 1116);
    const auto *kk_1117 = buffer.data(kk + 1117);
    const auto *kk_1118 = buffer.data(kk + 1118);
    const auto *kk_1119 = buffer.data(kk + 1119);
    const auto *kk_1120 = buffer.data(kk + 1120);

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, hk_454, hk_455, hk_456, hk_457, \
                         hk_458, kk_922, kk_923, kk_924, kk_925, \
                         kk_926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -3.0 * hk_454[k]
                   + f_0 * kk_922[k];

        t_671[k] = -3.0 * hk_455[k]
                   + f_0 * kk_923[k];

        t_672[k] = -3.0 * hk_456[k]
                   + f_0 * kk_924[k];

        t_673[k] = -3.0 * hk_457[k]
                   + f_0 * kk_925[k];

        t_674[k] = -3.0 * hk_458[k]
                   + f_0 * kk_926[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, hk_459, hk_460, hk_461, hk_462, \
                         hk_463, kk_927, kk_928, kk_929, kk_930, \
                         kk_931 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -3.0 * hk_459[k]
                   + f_0 * kk_927[k];

        t_676[k] = -3.0 * hk_460[k]
                   + f_0 * kk_928[k];

        t_677[k] = -3.0 * hk_461[k]
                   + f_0 * kk_929[k];

        t_678[k] = -3.0 * hk_462[k]
                   + f_0 * kk_930[k];

        t_679[k] = -3.0 * hk_463[k]
                   + f_0 * kk_931[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, hk_464, hk_465, hk_466, hk_467, \
                         hk_468, kk_932, kk_933, kk_934, kk_935, \
                         kk_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -3.0 * hk_464[k]
                   + f_0 * kk_932[k];

        t_681[k] = -3.0 * hk_465[k]
                   + f_0 * kk_933[k];

        t_682[k] = -3.0 * hk_466[k]
                   + f_0 * kk_934[k];

        t_683[k] = -3.0 * hk_467[k]
                   + f_0 * kk_935[k];

        t_684[k] = -4.0 * hk_468[k]
                   + f_0 * kk_936[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, hk_469, hk_470, hk_471, hk_472, \
                         hk_473, kk_937, kk_938, kk_939, kk_940, \
                         kk_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -4.0 * hk_469[k]
                   + f_0 * kk_937[k];

        t_686[k] = -4.0 * hk_470[k]
                   + f_0 * kk_938[k];

        t_687[k] = -4.0 * hk_471[k]
                   + f_0 * kk_939[k];

        t_688[k] = -4.0 * hk_472[k]
                   + f_0 * kk_940[k];

        t_689[k] = -4.0 * hk_473[k]
                   + f_0 * kk_941[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, hk_474, hk_475, hk_476, hk_477, \
                         hk_478, kk_942, kk_943, kk_944, kk_945, \
                         kk_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -4.0 * hk_474[k]
                   + f_0 * kk_942[k];

        t_691[k] = -4.0 * hk_475[k]
                   + f_0 * kk_943[k];

        t_692[k] = -4.0 * hk_476[k]
                   + f_0 * kk_944[k];

        t_693[k] = -4.0 * hk_477[k]
                   + f_0 * kk_945[k];

        t_694[k] = -4.0 * hk_478[k]
                   + f_0 * kk_946[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, hk_479, hk_480, hk_481, hk_482, \
                         hk_483, kk_947, kk_948, kk_949, kk_950, \
                         kk_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -4.0 * hk_479[k]
                   + f_0 * kk_947[k];

        t_696[k] = -4.0 * hk_480[k]
                   + f_0 * kk_948[k];

        t_697[k] = -4.0 * hk_481[k]
                   + f_0 * kk_949[k];

        t_698[k] = -4.0 * hk_482[k]
                   + f_0 * kk_950[k];

        t_699[k] = -4.0 * hk_483[k]
                   + f_0 * kk_951[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, hk_484, hk_485, hk_486, hk_487, \
                         hk_488, kk_952, kk_953, kk_954, kk_955, \
                         kk_956 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -4.0 * hk_484[k]
                   + f_0 * kk_952[k];

        t_701[k] = -4.0 * hk_485[k]
                   + f_0 * kk_953[k];

        t_702[k] = -4.0 * hk_486[k]
                   + f_0 * kk_954[k];

        t_703[k] = -4.0 * hk_487[k]
                   + f_0 * kk_955[k];

        t_704[k] = -4.0 * hk_488[k]
                   + f_0 * kk_956[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, hk_489, hk_490, hk_491, hk_492, \
                         hk_493, kk_957, kk_958, kk_959, kk_960, \
                         kk_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -4.0 * hk_489[k]
                   + f_0 * kk_957[k];

        t_706[k] = -4.0 * hk_490[k]
                   + f_0 * kk_958[k];

        t_707[k] = -4.0 * hk_491[k]
                   + f_0 * kk_959[k];

        t_708[k] = -4.0 * hk_492[k]
                   + f_0 * kk_960[k];

        t_709[k] = -4.0 * hk_493[k]
                   + f_0 * kk_961[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, hk_494, hk_495, hk_496, hk_497, \
                         hk_498, kk_962, kk_963, kk_964, kk_965, \
                         kk_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -4.0 * hk_494[k]
                   + f_0 * kk_962[k];

        t_711[k] = -4.0 * hk_495[k]
                   + f_0 * kk_963[k];

        t_712[k] = -4.0 * hk_496[k]
                   + f_0 * kk_964[k];

        t_713[k] = -4.0 * hk_497[k]
                   + f_0 * kk_965[k];

        t_714[k] = -4.0 * hk_498[k]
                   + f_0 * kk_966[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, hk_499, hk_500, hk_501, hk_502, \
                         hk_503, kk_967, kk_968, kk_969, kk_970, \
                         kk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -4.0 * hk_499[k]
                   + f_0 * kk_967[k];

        t_716[k] = -4.0 * hk_500[k]
                   + f_0 * kk_968[k];

        t_717[k] = -4.0 * hk_501[k]
                   + f_0 * kk_969[k];

        t_718[k] = -4.0 * hk_502[k]
                   + f_0 * kk_970[k];

        t_719[k] = -4.0 * hk_503[k]
                   + f_0 * kk_971[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, hk_504, hk_505, hk_506, hk_507, \
                         hk_508, kk_972, kk_973, kk_974, kk_975, \
                         kk_976 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -5.0 * hk_504[k]
                   + f_0 * kk_972[k];

        t_721[k] = -5.0 * hk_505[k]
                   + f_0 * kk_973[k];

        t_722[k] = -5.0 * hk_506[k]
                   + f_0 * kk_974[k];

        t_723[k] = -5.0 * hk_507[k]
                   + f_0 * kk_975[k];

        t_724[k] = -5.0 * hk_508[k]
                   + f_0 * kk_976[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, hk_509, hk_510, hk_511, hk_512, \
                         hk_513, kk_977, kk_978, kk_979, kk_980, \
                         kk_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -5.0 * hk_509[k]
                   + f_0 * kk_977[k];

        t_726[k] = -5.0 * hk_510[k]
                   + f_0 * kk_978[k];

        t_727[k] = -5.0 * hk_511[k]
                   + f_0 * kk_979[k];

        t_728[k] = -5.0 * hk_512[k]
                   + f_0 * kk_980[k];

        t_729[k] = -5.0 * hk_513[k]
                   + f_0 * kk_981[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, hk_514, hk_515, hk_516, hk_517, \
                         hk_518, kk_982, kk_983, kk_984, kk_985, \
                         kk_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -5.0 * hk_514[k]
                   + f_0 * kk_982[k];

        t_731[k] = -5.0 * hk_515[k]
                   + f_0 * kk_983[k];

        t_732[k] = -5.0 * hk_516[k]
                   + f_0 * kk_984[k];

        t_733[k] = -5.0 * hk_517[k]
                   + f_0 * kk_985[k];

        t_734[k] = -5.0 * hk_518[k]
                   + f_0 * kk_986[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, hk_519, hk_520, hk_521, hk_522, \
                         hk_523, kk_987, kk_988, kk_989, kk_990, \
                         kk_991 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -5.0 * hk_519[k]
                   + f_0 * kk_987[k];

        t_736[k] = -5.0 * hk_520[k]
                   + f_0 * kk_988[k];

        t_737[k] = -5.0 * hk_521[k]
                   + f_0 * kk_989[k];

        t_738[k] = -5.0 * hk_522[k]
                   + f_0 * kk_990[k];

        t_739[k] = -5.0 * hk_523[k]
                   + f_0 * kk_991[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, hk_524, hk_525, hk_526, hk_527, \
                         hk_528, kk_992, kk_993, kk_994, kk_995, \
                         kk_996 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -5.0 * hk_524[k]
                   + f_0 * kk_992[k];

        t_741[k] = -5.0 * hk_525[k]
                   + f_0 * kk_993[k];

        t_742[k] = -5.0 * hk_526[k]
                   + f_0 * kk_994[k];

        t_743[k] = -5.0 * hk_527[k]
                   + f_0 * kk_995[k];

        t_744[k] = -5.0 * hk_528[k]
                   + f_0 * kk_996[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, hk_529, hk_530, hk_531, hk_532, \
                         hk_533, kk_997, kk_998, kk_999, kk_1000, \
                         kk_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -5.0 * hk_529[k]
                   + f_0 * kk_997[k];

        t_746[k] = -5.0 * hk_530[k]
                   + f_0 * kk_998[k];

        t_747[k] = -5.0 * hk_531[k]
                   + f_0 * kk_999[k];

        t_748[k] = -5.0 * hk_532[k]
                   + f_0 * kk_1000[k];

        t_749[k] = -5.0 * hk_533[k]
                   + f_0 * kk_1001[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, hk_534, hk_535, hk_536, hk_537, \
                         hk_538, kk_1002, kk_1003, kk_1004, kk_1005, \
                         kk_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -5.0 * hk_534[k]
                   + f_0 * kk_1002[k];

        t_751[k] = -5.0 * hk_535[k]
                   + f_0 * kk_1003[k];

        t_752[k] = -5.0 * hk_536[k]
                   + f_0 * kk_1004[k];

        t_753[k] = -5.0 * hk_537[k]
                   + f_0 * kk_1005[k];

        t_754[k] = -5.0 * hk_538[k]
                   + f_0 * kk_1006[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, t_760, t_761, hk_539, kk_1007, \
                         kk_1044, kk_1045, kk_1046, kk_1047, kk_1048, \
                         kk_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -5.0 * hk_539[k]
                   + f_0 * kk_1007[k];

        t_756[k] = f_0 * kk_1044[k];

        t_757[k] = f_0 * kk_1045[k];

        t_758[k] = f_0 * kk_1046[k];

        t_759[k] = f_0 * kk_1047[k];

        t_760[k] = f_0 * kk_1048[k];

        t_761[k] = f_0 * kk_1049[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, t_766, t_767, t_768, t_769, kk_1050, \
                         kk_1051, kk_1052, kk_1053, kk_1054, kk_1055, kk_1056, \
                         kk_1057 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_0 * kk_1050[k];

        t_763[k] = f_0 * kk_1051[k];

        t_764[k] = f_0 * kk_1052[k];

        t_765[k] = f_0 * kk_1053[k];

        t_766[k] = f_0 * kk_1054[k];

        t_767[k] = f_0 * kk_1055[k];

        t_768[k] = f_0 * kk_1056[k];

        t_769[k] = f_0 * kk_1057[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, t_775, t_776, t_777, kk_1058, \
                         kk_1059, kk_1060, kk_1061, kk_1062, kk_1063, kk_1064, \
                         kk_1065 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_0 * kk_1058[k];

        t_771[k] = f_0 * kk_1059[k];

        t_772[k] = f_0 * kk_1060[k];

        t_773[k] = f_0 * kk_1061[k];

        t_774[k] = f_0 * kk_1062[k];

        t_775[k] = f_0 * kk_1063[k];

        t_776[k] = f_0 * kk_1064[k];

        t_777[k] = f_0 * kk_1065[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, t_783, t_784, t_785, kk_1066, \
                         kk_1067, kk_1068, kk_1069, kk_1070, kk_1071, kk_1072, \
                         kk_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_0 * kk_1066[k];

        t_779[k] = f_0 * kk_1067[k];

        t_780[k] = f_0 * kk_1068[k];

        t_781[k] = f_0 * kk_1069[k];

        t_782[k] = f_0 * kk_1070[k];

        t_783[k] = f_0 * kk_1071[k];

        t_784[k] = f_0 * kk_1072[k];

        t_785[k] = f_0 * kk_1073[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, t_790, t_791, t_792, hk_540, kk_1074, \
                         kk_1075, kk_1076, kk_1077, kk_1078, kk_1079, \
                         kk_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_0 * kk_1074[k];

        t_787[k] = f_0 * kk_1075[k];

        t_788[k] = f_0 * kk_1076[k];

        t_789[k] = f_0 * kk_1077[k];

        t_790[k] = f_0 * kk_1078[k];

        t_791[k] = f_0 * kk_1079[k];

        t_792[k] = -hk_540[k]
                   + f_0 * kk_1080[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, t_796, t_797, hk_541, hk_542, hk_543, hk_544, \
                         hk_545, kk_1081, kk_1082, kk_1083, kk_1084, \
                         kk_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = -hk_541[k]
                   + f_0 * kk_1081[k];

        t_794[k] = -hk_542[k]
                   + f_0 * kk_1082[k];

        t_795[k] = -hk_543[k]
                   + f_0 * kk_1083[k];

        t_796[k] = -hk_544[k]
                   + f_0 * kk_1084[k];

        t_797[k] = -hk_545[k]
                   + f_0 * kk_1085[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, t_801, t_802, hk_546, hk_547, hk_548, hk_549, \
                         hk_550, kk_1086, kk_1087, kk_1088, kk_1089, \
                         kk_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = -hk_546[k]
                   + f_0 * kk_1086[k];

        t_799[k] = -hk_547[k]
                   + f_0 * kk_1087[k];

        t_800[k] = -hk_548[k]
                   + f_0 * kk_1088[k];

        t_801[k] = -hk_549[k]
                   + f_0 * kk_1089[k];

        t_802[k] = -hk_550[k]
                   + f_0 * kk_1090[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, t_806, t_807, hk_551, hk_552, hk_553, hk_554, \
                         hk_555, kk_1091, kk_1092, kk_1093, kk_1094, \
                         kk_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = -hk_551[k]
                   + f_0 * kk_1091[k];

        t_804[k] = -hk_552[k]
                   + f_0 * kk_1092[k];

        t_805[k] = -hk_553[k]
                   + f_0 * kk_1093[k];

        t_806[k] = -hk_554[k]
                   + f_0 * kk_1094[k];

        t_807[k] = -hk_555[k]
                   + f_0 * kk_1095[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, t_812, hk_556, hk_557, hk_558, hk_559, \
                         hk_560, kk_1096, kk_1097, kk_1098, kk_1099, \
                         kk_1100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = -hk_556[k]
                   + f_0 * kk_1096[k];

        t_809[k] = -hk_557[k]
                   + f_0 * kk_1097[k];

        t_810[k] = -hk_558[k]
                   + f_0 * kk_1098[k];

        t_811[k] = -hk_559[k]
                   + f_0 * kk_1099[k];

        t_812[k] = -hk_560[k]
                   + f_0 * kk_1100[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, t_816, t_817, hk_561, hk_562, hk_563, hk_564, \
                         hk_565, kk_1101, kk_1102, kk_1103, kk_1104, \
                         kk_1105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = -hk_561[k]
                   + f_0 * kk_1101[k];

        t_814[k] = -hk_562[k]
                   + f_0 * kk_1102[k];

        t_815[k] = -hk_563[k]
                   + f_0 * kk_1103[k];

        t_816[k] = -hk_564[k]
                   + f_0 * kk_1104[k];

        t_817[k] = -hk_565[k]
                   + f_0 * kk_1105[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, t_821, t_822, hk_566, hk_567, hk_568, hk_569, \
                         hk_570, kk_1106, kk_1107, kk_1108, kk_1109, \
                         kk_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = -hk_566[k]
                   + f_0 * kk_1106[k];

        t_819[k] = -hk_567[k]
                   + f_0 * kk_1107[k];

        t_820[k] = -hk_568[k]
                   + f_0 * kk_1108[k];

        t_821[k] = -hk_569[k]
                   + f_0 * kk_1109[k];

        t_822[k] = -hk_570[k]
                   + f_0 * kk_1110[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, t_826, t_827, hk_571, hk_572, hk_573, hk_574, \
                         hk_575, kk_1111, kk_1112, kk_1113, kk_1114, \
                         kk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = -hk_571[k]
                   + f_0 * kk_1111[k];

        t_824[k] = -hk_572[k]
                   + f_0 * kk_1112[k];

        t_825[k] = -hk_573[k]
                   + f_0 * kk_1113[k];

        t_826[k] = -hk_574[k]
                   + f_0 * kk_1114[k];

        t_827[k] = -hk_575[k]
                   + f_0 * kk_1115[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, t_831, t_832, hk_576, hk_577, hk_578, hk_579, \
                         hk_580, kk_1116, kk_1117, kk_1118, kk_1119, \
                         kk_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = -2.0 * hk_576[k]
                   + f_0 * kk_1116[k];

        t_829[k] = -2.0 * hk_577[k]
                   + f_0 * kk_1117[k];

        t_830[k] = -2.0 * hk_578[k]
                   + f_0 * kk_1118[k];

        t_831[k] = -2.0 * hk_579[k]
                   + f_0 * kk_1119[k];

        t_832[k] = -2.0 * hk_580[k]
                   + f_0 * kk_1120[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_2_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_581 = buffer.data(hk + 581);
    const auto *hk_582 = buffer.data(hk + 582);
    const auto *hk_583 = buffer.data(hk + 583);
    const auto *hk_584 = buffer.data(hk + 584);
    const auto *hk_585 = buffer.data(hk + 585);
    const auto *hk_586 = buffer.data(hk + 586);
    const auto *hk_587 = buffer.data(hk + 587);
    const auto *hk_588 = buffer.data(hk + 588);
    const auto *hk_589 = buffer.data(hk + 589);
    const auto *hk_590 = buffer.data(hk + 590);
    const auto *hk_591 = buffer.data(hk + 591);
    const auto *hk_592 = buffer.data(hk + 592);
    const auto *hk_593 = buffer.data(hk + 593);
    const auto *hk_594 = buffer.data(hk + 594);
    const auto *hk_595 = buffer.data(hk + 595);
    const auto *hk_596 = buffer.data(hk + 596);
    const auto *hk_597 = buffer.data(hk + 597);
    const auto *hk_598 = buffer.data(hk + 598);
    const auto *hk_599 = buffer.data(hk + 599);
    const auto *hk_600 = buffer.data(hk + 600);
    const auto *hk_601 = buffer.data(hk + 601);
    const auto *hk_602 = buffer.data(hk + 602);
    const auto *hk_603 = buffer.data(hk + 603);
    const auto *hk_604 = buffer.data(hk + 604);
    const auto *hk_605 = buffer.data(hk + 605);
    const auto *hk_606 = buffer.data(hk + 606);
    const auto *hk_607 = buffer.data(hk + 607);
    const auto *hk_608 = buffer.data(hk + 608);
    const auto *hk_609 = buffer.data(hk + 609);
    const auto *hk_610 = buffer.data(hk + 610);
    const auto *hk_611 = buffer.data(hk + 611);
    const auto *hk_612 = buffer.data(hk + 612);
    const auto *hk_613 = buffer.data(hk + 613);
    const auto *hk_614 = buffer.data(hk + 614);
    const auto *hk_615 = buffer.data(hk + 615);
    const auto *hk_616 = buffer.data(hk + 616);
    const auto *hk_617 = buffer.data(hk + 617);
    const auto *hk_618 = buffer.data(hk + 618);
    const auto *hk_619 = buffer.data(hk + 619);
    const auto *hk_620 = buffer.data(hk + 620);
    const auto *hk_621 = buffer.data(hk + 621);
    const auto *hk_622 = buffer.data(hk + 622);
    const auto *hk_623 = buffer.data(hk + 623);
    const auto *hk_624 = buffer.data(hk + 624);
    const auto *hk_625 = buffer.data(hk + 625);
    const auto *hk_626 = buffer.data(hk + 626);
    const auto *hk_627 = buffer.data(hk + 627);
    const auto *hk_628 = buffer.data(hk + 628);
    const auto *hk_629 = buffer.data(hk + 629);
    const auto *hk_630 = buffer.data(hk + 630);
    const auto *hk_631 = buffer.data(hk + 631);
    const auto *hk_632 = buffer.data(hk + 632);
    const auto *hk_633 = buffer.data(hk + 633);
    const auto *hk_634 = buffer.data(hk + 634);
    const auto *hk_635 = buffer.data(hk + 635);
    const auto *hk_636 = buffer.data(hk + 636);
    const auto *hk_637 = buffer.data(hk + 637);
    const auto *hk_638 = buffer.data(hk + 638);
    const auto *hk_639 = buffer.data(hk + 639);
    const auto *hk_640 = buffer.data(hk + 640);
    const auto *hk_641 = buffer.data(hk + 641);
    const auto *hk_642 = buffer.data(hk + 642);
    const auto *hk_643 = buffer.data(hk + 643);
    const auto *hk_644 = buffer.data(hk + 644);
    const auto *hk_645 = buffer.data(hk + 645);
    const auto *hk_646 = buffer.data(hk + 646);
    const auto *hk_647 = buffer.data(hk + 647);
    const auto *hk_648 = buffer.data(hk + 648);
    const auto *hk_649 = buffer.data(hk + 649);
    const auto *hk_650 = buffer.data(hk + 650);
    const auto *hk_651 = buffer.data(hk + 651);
    const auto *hk_652 = buffer.data(hk + 652);
    const auto *hk_653 = buffer.data(hk + 653);
    const auto *hk_654 = buffer.data(hk + 654);
    const auto *hk_655 = buffer.data(hk + 655);
    const auto *hk_656 = buffer.data(hk + 656);
    const auto *hk_657 = buffer.data(hk + 657);
    const auto *hk_658 = buffer.data(hk + 658);
    const auto *hk_659 = buffer.data(hk + 659);
    const auto *hk_660 = buffer.data(hk + 660);
    const auto *hk_661 = buffer.data(hk + 661);
    const auto *hk_662 = buffer.data(hk + 662);
    const auto *hk_663 = buffer.data(hk + 663);
    const auto *hk_664 = buffer.data(hk + 664);
    const auto *hk_665 = buffer.data(hk + 665);
    const auto *hk_666 = buffer.data(hk + 666);
    const auto *hk_667 = buffer.data(hk + 667);
    const auto *hk_668 = buffer.data(hk + 668);
    const auto *hk_669 = buffer.data(hk + 669);
    const auto *hk_670 = buffer.data(hk + 670);
    const auto *hk_671 = buffer.data(hk + 671);
    const auto *hk_672 = buffer.data(hk + 672);
    const auto *hk_673 = buffer.data(hk + 673);
    const auto *hk_674 = buffer.data(hk + 674);
    const auto *hk_675 = buffer.data(hk + 675);
    const auto *hk_676 = buffer.data(hk + 676);
    const auto *hk_677 = buffer.data(hk + 677);
    const auto *hk_678 = buffer.data(hk + 678);
    const auto *hk_679 = buffer.data(hk + 679);
    const auto *hk_680 = buffer.data(hk + 680);
    const auto *hk_681 = buffer.data(hk + 681);
    const auto *hk_682 = buffer.data(hk + 682);
    const auto *hk_683 = buffer.data(hk + 683);
    const auto *hk_684 = buffer.data(hk + 684);
    const auto *hk_685 = buffer.data(hk + 685);
    const auto *hk_686 = buffer.data(hk + 686);
    const auto *hk_687 = buffer.data(hk + 687);
    const auto *hk_688 = buffer.data(hk + 688);
    const auto *hk_689 = buffer.data(hk + 689);
    const auto *hk_690 = buffer.data(hk + 690);
    const auto *hk_691 = buffer.data(hk + 691);
    const auto *hk_692 = buffer.data(hk + 692);
    const auto *hk_693 = buffer.data(hk + 693);
    const auto *hk_694 = buffer.data(hk + 694);
    const auto *hk_695 = buffer.data(hk + 695);
    const auto *hk_696 = buffer.data(hk + 696);
    const auto *hk_697 = buffer.data(hk + 697);
    const auto *hk_698 = buffer.data(hk + 698);
    const auto *hk_699 = buffer.data(hk + 699);
    const auto *hk_700 = buffer.data(hk + 700);
    const auto *hk_701 = buffer.data(hk + 701);
    const auto *hk_702 = buffer.data(hk + 702);
    const auto *hk_703 = buffer.data(hk + 703);
    const auto *hk_704 = buffer.data(hk + 704);
    const auto *hk_705 = buffer.data(hk + 705);
    const auto *hk_706 = buffer.data(hk + 706);
    const auto *hk_707 = buffer.data(hk + 707);
    const auto *hk_708 = buffer.data(hk + 708);
    const auto *hk_709 = buffer.data(hk + 709);
    const auto *hk_710 = buffer.data(hk + 710);
    const auto *hk_711 = buffer.data(hk + 711);
    const auto *hk_712 = buffer.data(hk + 712);
    const auto *hk_713 = buffer.data(hk + 713);
    const auto *hk_714 = buffer.data(hk + 714);
    const auto *hk_715 = buffer.data(hk + 715);
    const auto *hk_716 = buffer.data(hk + 716);
    const auto *hk_717 = buffer.data(hk + 717);
    const auto *hk_718 = buffer.data(hk + 718);
    const auto *hk_719 = buffer.data(hk + 719);
    const auto *hk_720 = buffer.data(hk + 720);
    const auto *hk_721 = buffer.data(hk + 721);
    const auto *hk_722 = buffer.data(hk + 722);
    const auto *hk_723 = buffer.data(hk + 723);
    const auto *hk_724 = buffer.data(hk + 724);
    const auto *hk_725 = buffer.data(hk + 725);
    const auto *hk_726 = buffer.data(hk + 726);
    const auto *hk_727 = buffer.data(hk + 727);
    const auto *hk_728 = buffer.data(hk + 728);
    const auto *hk_729 = buffer.data(hk + 729);
    const auto *hk_730 = buffer.data(hk + 730);

    const auto *kk_1121 = buffer.data(kk + 1121);
    const auto *kk_1122 = buffer.data(kk + 1122);
    const auto *kk_1123 = buffer.data(kk + 1123);
    const auto *kk_1124 = buffer.data(kk + 1124);
    const auto *kk_1125 = buffer.data(kk + 1125);
    const auto *kk_1126 = buffer.data(kk + 1126);
    const auto *kk_1127 = buffer.data(kk + 1127);
    const auto *kk_1128 = buffer.data(kk + 1128);
    const auto *kk_1129 = buffer.data(kk + 1129);
    const auto *kk_1130 = buffer.data(kk + 1130);
    const auto *kk_1131 = buffer.data(kk + 1131);
    const auto *kk_1132 = buffer.data(kk + 1132);
    const auto *kk_1133 = buffer.data(kk + 1133);
    const auto *kk_1134 = buffer.data(kk + 1134);
    const auto *kk_1135 = buffer.data(kk + 1135);
    const auto *kk_1136 = buffer.data(kk + 1136);
    const auto *kk_1137 = buffer.data(kk + 1137);
    const auto *kk_1138 = buffer.data(kk + 1138);
    const auto *kk_1139 = buffer.data(kk + 1139);
    const auto *kk_1140 = buffer.data(kk + 1140);
    const auto *kk_1141 = buffer.data(kk + 1141);
    const auto *kk_1142 = buffer.data(kk + 1142);
    const auto *kk_1143 = buffer.data(kk + 1143);
    const auto *kk_1144 = buffer.data(kk + 1144);
    const auto *kk_1145 = buffer.data(kk + 1145);
    const auto *kk_1146 = buffer.data(kk + 1146);
    const auto *kk_1147 = buffer.data(kk + 1147);
    const auto *kk_1148 = buffer.data(kk + 1148);
    const auto *kk_1149 = buffer.data(kk + 1149);
    const auto *kk_1150 = buffer.data(kk + 1150);
    const auto *kk_1151 = buffer.data(kk + 1151);
    const auto *kk_1152 = buffer.data(kk + 1152);
    const auto *kk_1153 = buffer.data(kk + 1153);
    const auto *kk_1154 = buffer.data(kk + 1154);
    const auto *kk_1155 = buffer.data(kk + 1155);
    const auto *kk_1156 = buffer.data(kk + 1156);
    const auto *kk_1157 = buffer.data(kk + 1157);
    const auto *kk_1158 = buffer.data(kk + 1158);
    const auto *kk_1159 = buffer.data(kk + 1159);
    const auto *kk_1160 = buffer.data(kk + 1160);
    const auto *kk_1161 = buffer.data(kk + 1161);
    const auto *kk_1162 = buffer.data(kk + 1162);
    const auto *kk_1163 = buffer.data(kk + 1163);
    const auto *kk_1164 = buffer.data(kk + 1164);
    const auto *kk_1165 = buffer.data(kk + 1165);
    const auto *kk_1166 = buffer.data(kk + 1166);
    const auto *kk_1167 = buffer.data(kk + 1167);
    const auto *kk_1168 = buffer.data(kk + 1168);
    const auto *kk_1169 = buffer.data(kk + 1169);
    const auto *kk_1170 = buffer.data(kk + 1170);
    const auto *kk_1171 = buffer.data(kk + 1171);
    const auto *kk_1172 = buffer.data(kk + 1172);
    const auto *kk_1173 = buffer.data(kk + 1173);
    const auto *kk_1174 = buffer.data(kk + 1174);
    const auto *kk_1175 = buffer.data(kk + 1175);
    const auto *kk_1176 = buffer.data(kk + 1176);
    const auto *kk_1177 = buffer.data(kk + 1177);
    const auto *kk_1178 = buffer.data(kk + 1178);
    const auto *kk_1179 = buffer.data(kk + 1179);
    const auto *kk_1180 = buffer.data(kk + 1180);
    const auto *kk_1181 = buffer.data(kk + 1181);
    const auto *kk_1182 = buffer.data(kk + 1182);
    const auto *kk_1183 = buffer.data(kk + 1183);
    const auto *kk_1184 = buffer.data(kk + 1184);
    const auto *kk_1185 = buffer.data(kk + 1185);
    const auto *kk_1186 = buffer.data(kk + 1186);
    const auto *kk_1187 = buffer.data(kk + 1187);
    const auto *kk_1188 = buffer.data(kk + 1188);
    const auto *kk_1189 = buffer.data(kk + 1189);
    const auto *kk_1190 = buffer.data(kk + 1190);
    const auto *kk_1191 = buffer.data(kk + 1191);
    const auto *kk_1192 = buffer.data(kk + 1192);
    const auto *kk_1193 = buffer.data(kk + 1193);
    const auto *kk_1194 = buffer.data(kk + 1194);
    const auto *kk_1195 = buffer.data(kk + 1195);
    const auto *kk_1196 = buffer.data(kk + 1196);
    const auto *kk_1197 = buffer.data(kk + 1197);
    const auto *kk_1198 = buffer.data(kk + 1198);
    const auto *kk_1199 = buffer.data(kk + 1199);
    const auto *kk_1200 = buffer.data(kk + 1200);
    const auto *kk_1201 = buffer.data(kk + 1201);
    const auto *kk_1202 = buffer.data(kk + 1202);
    const auto *kk_1203 = buffer.data(kk + 1203);
    const auto *kk_1204 = buffer.data(kk + 1204);
    const auto *kk_1205 = buffer.data(kk + 1205);
    const auto *kk_1206 = buffer.data(kk + 1206);
    const auto *kk_1207 = buffer.data(kk + 1207);
    const auto *kk_1208 = buffer.data(kk + 1208);
    const auto *kk_1209 = buffer.data(kk + 1209);
    const auto *kk_1210 = buffer.data(kk + 1210);
    const auto *kk_1211 = buffer.data(kk + 1211);
    const auto *kk_1212 = buffer.data(kk + 1212);
    const auto *kk_1213 = buffer.data(kk + 1213);
    const auto *kk_1214 = buffer.data(kk + 1214);
    const auto *kk_1215 = buffer.data(kk + 1215);
    const auto *kk_1216 = buffer.data(kk + 1216);
    const auto *kk_1217 = buffer.data(kk + 1217);
    const auto *kk_1218 = buffer.data(kk + 1218);
    const auto *kk_1219 = buffer.data(kk + 1219);
    const auto *kk_1220 = buffer.data(kk + 1220);
    const auto *kk_1221 = buffer.data(kk + 1221);
    const auto *kk_1222 = buffer.data(kk + 1222);
    const auto *kk_1223 = buffer.data(kk + 1223);
    const auto *kk_1224 = buffer.data(kk + 1224);
    const auto *kk_1225 = buffer.data(kk + 1225);
    const auto *kk_1226 = buffer.data(kk + 1226);
    const auto *kk_1227 = buffer.data(kk + 1227);
    const auto *kk_1228 = buffer.data(kk + 1228);
    const auto *kk_1229 = buffer.data(kk + 1229);
    const auto *kk_1230 = buffer.data(kk + 1230);
    const auto *kk_1231 = buffer.data(kk + 1231);
    const auto *kk_1232 = buffer.data(kk + 1232);
    const auto *kk_1233 = buffer.data(kk + 1233);
    const auto *kk_1234 = buffer.data(kk + 1234);
    const auto *kk_1235 = buffer.data(kk + 1235);
    const auto *kk_1236 = buffer.data(kk + 1236);
    const auto *kk_1237 = buffer.data(kk + 1237);
    const auto *kk_1238 = buffer.data(kk + 1238);
    const auto *kk_1239 = buffer.data(kk + 1239);
    const auto *kk_1240 = buffer.data(kk + 1240);
    const auto *kk_1241 = buffer.data(kk + 1241);
    const auto *kk_1242 = buffer.data(kk + 1242);
    const auto *kk_1243 = buffer.data(kk + 1243);
    const auto *kk_1244 = buffer.data(kk + 1244);
    const auto *kk_1245 = buffer.data(kk + 1245);
    const auto *kk_1246 = buffer.data(kk + 1246);
    const auto *kk_1247 = buffer.data(kk + 1247);
    const auto *kk_1248 = buffer.data(kk + 1248);
    const auto *kk_1249 = buffer.data(kk + 1249);
    const auto *kk_1250 = buffer.data(kk + 1250);
    const auto *kk_1251 = buffer.data(kk + 1251);
    const auto *kk_1252 = buffer.data(kk + 1252);
    const auto *kk_1253 = buffer.data(kk + 1253);
    const auto *kk_1254 = buffer.data(kk + 1254);
    const auto *kk_1255 = buffer.data(kk + 1255);
    const auto *kk_1256 = buffer.data(kk + 1256);
    const auto *kk_1257 = buffer.data(kk + 1257);
    const auto *kk_1258 = buffer.data(kk + 1258);
    const auto *kk_1259 = buffer.data(kk + 1259);
    const auto *kk_1260 = buffer.data(kk + 1260);
    const auto *kk_1261 = buffer.data(kk + 1261);
    const auto *kk_1262 = buffer.data(kk + 1262);
    const auto *kk_1263 = buffer.data(kk + 1263);
    const auto *kk_1264 = buffer.data(kk + 1264);
    const auto *kk_1265 = buffer.data(kk + 1265);
    const auto *kk_1266 = buffer.data(kk + 1266);
    const auto *kk_1267 = buffer.data(kk + 1267);
    const auto *kk_1268 = buffer.data(kk + 1268);
    const auto *kk_1269 = buffer.data(kk + 1269);
    const auto *kk_1270 = buffer.data(kk + 1270);

#pragma omp simd aligned(t_833, t_834, t_835, t_836, t_837, hk_581, hk_582, hk_583, hk_584, \
                         hk_585, kk_1121, kk_1122, kk_1123, kk_1124, \
                         kk_1125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = -2.0 * hk_581[k]
                   + f_0 * kk_1121[k];

        t_834[k] = -2.0 * hk_582[k]
                   + f_0 * kk_1122[k];

        t_835[k] = -2.0 * hk_583[k]
                   + f_0 * kk_1123[k];

        t_836[k] = -2.0 * hk_584[k]
                   + f_0 * kk_1124[k];

        t_837[k] = -2.0 * hk_585[k]
                   + f_0 * kk_1125[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, t_841, t_842, hk_586, hk_587, hk_588, hk_589, \
                         hk_590, kk_1126, kk_1127, kk_1128, kk_1129, \
                         kk_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = -2.0 * hk_586[k]
                   + f_0 * kk_1126[k];

        t_839[k] = -2.0 * hk_587[k]
                   + f_0 * kk_1127[k];

        t_840[k] = -2.0 * hk_588[k]
                   + f_0 * kk_1128[k];

        t_841[k] = -2.0 * hk_589[k]
                   + f_0 * kk_1129[k];

        t_842[k] = -2.0 * hk_590[k]
                   + f_0 * kk_1130[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, t_846, t_847, hk_591, hk_592, hk_593, hk_594, \
                         hk_595, kk_1131, kk_1132, kk_1133, kk_1134, \
                         kk_1135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = -2.0 * hk_591[k]
                   + f_0 * kk_1131[k];

        t_844[k] = -2.0 * hk_592[k]
                   + f_0 * kk_1132[k];

        t_845[k] = -2.0 * hk_593[k]
                   + f_0 * kk_1133[k];

        t_846[k] = -2.0 * hk_594[k]
                   + f_0 * kk_1134[k];

        t_847[k] = -2.0 * hk_595[k]
                   + f_0 * kk_1135[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, t_852, hk_596, hk_597, hk_598, hk_599, \
                         hk_600, kk_1136, kk_1137, kk_1138, kk_1139, \
                         kk_1140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = -2.0 * hk_596[k]
                   + f_0 * kk_1136[k];

        t_849[k] = -2.0 * hk_597[k]
                   + f_0 * kk_1137[k];

        t_850[k] = -2.0 * hk_598[k]
                   + f_0 * kk_1138[k];

        t_851[k] = -2.0 * hk_599[k]
                   + f_0 * kk_1139[k];

        t_852[k] = -2.0 * hk_600[k]
                   + f_0 * kk_1140[k];
    }

#pragma omp simd aligned(t_853, t_854, t_855, t_856, t_857, hk_601, hk_602, hk_603, hk_604, \
                         hk_605, kk_1141, kk_1142, kk_1143, kk_1144, \
                         kk_1145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_853[k] = -2.0 * hk_601[k]
                   + f_0 * kk_1141[k];

        t_854[k] = -2.0 * hk_602[k]
                   + f_0 * kk_1142[k];

        t_855[k] = -2.0 * hk_603[k]
                   + f_0 * kk_1143[k];

        t_856[k] = -2.0 * hk_604[k]
                   + f_0 * kk_1144[k];

        t_857[k] = -2.0 * hk_605[k]
                   + f_0 * kk_1145[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, t_862, hk_606, hk_607, hk_608, hk_609, \
                         hk_610, kk_1146, kk_1147, kk_1148, kk_1149, \
                         kk_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = -2.0 * hk_606[k]
                   + f_0 * kk_1146[k];

        t_859[k] = -2.0 * hk_607[k]
                   + f_0 * kk_1147[k];

        t_860[k] = -2.0 * hk_608[k]
                   + f_0 * kk_1148[k];

        t_861[k] = -2.0 * hk_609[k]
                   + f_0 * kk_1149[k];

        t_862[k] = -2.0 * hk_610[k]
                   + f_0 * kk_1150[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, t_866, t_867, hk_611, hk_612, hk_613, hk_614, \
                         hk_615, kk_1151, kk_1152, kk_1153, kk_1154, \
                         kk_1155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = -2.0 * hk_611[k]
                   + f_0 * kk_1151[k];

        t_864[k] = -3.0 * hk_612[k]
                   + f_0 * kk_1152[k];

        t_865[k] = -3.0 * hk_613[k]
                   + f_0 * kk_1153[k];

        t_866[k] = -3.0 * hk_614[k]
                   + f_0 * kk_1154[k];

        t_867[k] = -3.0 * hk_615[k]
                   + f_0 * kk_1155[k];
    }

#pragma omp simd aligned(t_868, t_869, t_870, t_871, t_872, hk_616, hk_617, hk_618, hk_619, \
                         hk_620, kk_1156, kk_1157, kk_1158, kk_1159, \
                         kk_1160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_868[k] = -3.0 * hk_616[k]
                   + f_0 * kk_1156[k];

        t_869[k] = -3.0 * hk_617[k]
                   + f_0 * kk_1157[k];

        t_870[k] = -3.0 * hk_618[k]
                   + f_0 * kk_1158[k];

        t_871[k] = -3.0 * hk_619[k]
                   + f_0 * kk_1159[k];

        t_872[k] = -3.0 * hk_620[k]
                   + f_0 * kk_1160[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, t_877, hk_621, hk_622, hk_623, hk_624, \
                         hk_625, kk_1161, kk_1162, kk_1163, kk_1164, \
                         kk_1165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = -3.0 * hk_621[k]
                   + f_0 * kk_1161[k];

        t_874[k] = -3.0 * hk_622[k]
                   + f_0 * kk_1162[k];

        t_875[k] = -3.0 * hk_623[k]
                   + f_0 * kk_1163[k];

        t_876[k] = -3.0 * hk_624[k]
                   + f_0 * kk_1164[k];

        t_877[k] = -3.0 * hk_625[k]
                   + f_0 * kk_1165[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, t_881, t_882, hk_626, hk_627, hk_628, hk_629, \
                         hk_630, kk_1166, kk_1167, kk_1168, kk_1169, \
                         kk_1170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = -3.0 * hk_626[k]
                   + f_0 * kk_1166[k];

        t_879[k] = -3.0 * hk_627[k]
                   + f_0 * kk_1167[k];

        t_880[k] = -3.0 * hk_628[k]
                   + f_0 * kk_1168[k];

        t_881[k] = -3.0 * hk_629[k]
                   + f_0 * kk_1169[k];

        t_882[k] = -3.0 * hk_630[k]
                   + f_0 * kk_1170[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, t_887, hk_631, hk_632, hk_633, hk_634, \
                         hk_635, kk_1171, kk_1172, kk_1173, kk_1174, \
                         kk_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = -3.0 * hk_631[k]
                   + f_0 * kk_1171[k];

        t_884[k] = -3.0 * hk_632[k]
                   + f_0 * kk_1172[k];

        t_885[k] = -3.0 * hk_633[k]
                   + f_0 * kk_1173[k];

        t_886[k] = -3.0 * hk_634[k]
                   + f_0 * kk_1174[k];

        t_887[k] = -3.0 * hk_635[k]
                   + f_0 * kk_1175[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, t_892, hk_636, hk_637, hk_638, hk_639, \
                         hk_640, kk_1176, kk_1177, kk_1178, kk_1179, \
                         kk_1180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = -3.0 * hk_636[k]
                   + f_0 * kk_1176[k];

        t_889[k] = -3.0 * hk_637[k]
                   + f_0 * kk_1177[k];

        t_890[k] = -3.0 * hk_638[k]
                   + f_0 * kk_1178[k];

        t_891[k] = -3.0 * hk_639[k]
                   + f_0 * kk_1179[k];

        t_892[k] = -3.0 * hk_640[k]
                   + f_0 * kk_1180[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, t_896, t_897, hk_641, hk_642, hk_643, hk_644, \
                         hk_645, kk_1181, kk_1182, kk_1183, kk_1184, \
                         kk_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = -3.0 * hk_641[k]
                   + f_0 * kk_1181[k];

        t_894[k] = -3.0 * hk_642[k]
                   + f_0 * kk_1182[k];

        t_895[k] = -3.0 * hk_643[k]
                   + f_0 * kk_1183[k];

        t_896[k] = -3.0 * hk_644[k]
                   + f_0 * kk_1184[k];

        t_897[k] = -3.0 * hk_645[k]
                   + f_0 * kk_1185[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, t_902, hk_646, hk_647, hk_648, hk_649, \
                         hk_650, kk_1186, kk_1187, kk_1188, kk_1189, \
                         kk_1190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = -3.0 * hk_646[k]
                   + f_0 * kk_1186[k];

        t_899[k] = -3.0 * hk_647[k]
                   + f_0 * kk_1187[k];

        t_900[k] = -4.0 * hk_648[k]
                   + f_0 * kk_1188[k];

        t_901[k] = -4.0 * hk_649[k]
                   + f_0 * kk_1189[k];

        t_902[k] = -4.0 * hk_650[k]
                   + f_0 * kk_1190[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, t_906, t_907, hk_651, hk_652, hk_653, hk_654, \
                         hk_655, kk_1191, kk_1192, kk_1193, kk_1194, \
                         kk_1195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = -4.0 * hk_651[k]
                   + f_0 * kk_1191[k];

        t_904[k] = -4.0 * hk_652[k]
                   + f_0 * kk_1192[k];

        t_905[k] = -4.0 * hk_653[k]
                   + f_0 * kk_1193[k];

        t_906[k] = -4.0 * hk_654[k]
                   + f_0 * kk_1194[k];

        t_907[k] = -4.0 * hk_655[k]
                   + f_0 * kk_1195[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, t_911, t_912, hk_656, hk_657, hk_658, hk_659, \
                         hk_660, kk_1196, kk_1197, kk_1198, kk_1199, \
                         kk_1200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = -4.0 * hk_656[k]
                   + f_0 * kk_1196[k];

        t_909[k] = -4.0 * hk_657[k]
                   + f_0 * kk_1197[k];

        t_910[k] = -4.0 * hk_658[k]
                   + f_0 * kk_1198[k];

        t_911[k] = -4.0 * hk_659[k]
                   + f_0 * kk_1199[k];

        t_912[k] = -4.0 * hk_660[k]
                   + f_0 * kk_1200[k];
    }

#pragma omp simd aligned(t_913, t_914, t_915, t_916, t_917, hk_661, hk_662, hk_663, hk_664, \
                         hk_665, kk_1201, kk_1202, kk_1203, kk_1204, \
                         kk_1205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_913[k] = -4.0 * hk_661[k]
                   + f_0 * kk_1201[k];

        t_914[k] = -4.0 * hk_662[k]
                   + f_0 * kk_1202[k];

        t_915[k] = -4.0 * hk_663[k]
                   + f_0 * kk_1203[k];

        t_916[k] = -4.0 * hk_664[k]
                   + f_0 * kk_1204[k];

        t_917[k] = -4.0 * hk_665[k]
                   + f_0 * kk_1205[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, t_921, t_922, hk_666, hk_667, hk_668, hk_669, \
                         hk_670, kk_1206, kk_1207, kk_1208, kk_1209, \
                         kk_1210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = -4.0 * hk_666[k]
                   + f_0 * kk_1206[k];

        t_919[k] = -4.0 * hk_667[k]
                   + f_0 * kk_1207[k];

        t_920[k] = -4.0 * hk_668[k]
                   + f_0 * kk_1208[k];

        t_921[k] = -4.0 * hk_669[k]
                   + f_0 * kk_1209[k];

        t_922[k] = -4.0 * hk_670[k]
                   + f_0 * kk_1210[k];
    }

#pragma omp simd aligned(t_923, t_924, t_925, t_926, t_927, hk_671, hk_672, hk_673, hk_674, \
                         hk_675, kk_1211, kk_1212, kk_1213, kk_1214, \
                         kk_1215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_923[k] = -4.0 * hk_671[k]
                   + f_0 * kk_1211[k];

        t_924[k] = -4.0 * hk_672[k]
                   + f_0 * kk_1212[k];

        t_925[k] = -4.0 * hk_673[k]
                   + f_0 * kk_1213[k];

        t_926[k] = -4.0 * hk_674[k]
                   + f_0 * kk_1214[k];

        t_927[k] = -4.0 * hk_675[k]
                   + f_0 * kk_1215[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, t_931, t_932, hk_676, hk_677, hk_678, hk_679, \
                         hk_680, kk_1216, kk_1217, kk_1218, kk_1219, \
                         kk_1220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = -4.0 * hk_676[k]
                   + f_0 * kk_1216[k];

        t_929[k] = -4.0 * hk_677[k]
                   + f_0 * kk_1217[k];

        t_930[k] = -4.0 * hk_678[k]
                   + f_0 * kk_1218[k];

        t_931[k] = -4.0 * hk_679[k]
                   + f_0 * kk_1219[k];

        t_932[k] = -4.0 * hk_680[k]
                   + f_0 * kk_1220[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, t_937, hk_681, hk_682, hk_683, hk_684, \
                         hk_685, kk_1221, kk_1222, kk_1223, kk_1224, \
                         kk_1225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = -4.0 * hk_681[k]
                   + f_0 * kk_1221[k];

        t_934[k] = -4.0 * hk_682[k]
                   + f_0 * kk_1222[k];

        t_935[k] = -4.0 * hk_683[k]
                   + f_0 * kk_1223[k];

        t_936[k] = -5.0 * hk_684[k]
                   + f_0 * kk_1224[k];

        t_937[k] = -5.0 * hk_685[k]
                   + f_0 * kk_1225[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, t_942, hk_686, hk_687, hk_688, hk_689, \
                         hk_690, kk_1226, kk_1227, kk_1228, kk_1229, \
                         kk_1230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = -5.0 * hk_686[k]
                   + f_0 * kk_1226[k];

        t_939[k] = -5.0 * hk_687[k]
                   + f_0 * kk_1227[k];

        t_940[k] = -5.0 * hk_688[k]
                   + f_0 * kk_1228[k];

        t_941[k] = -5.0 * hk_689[k]
                   + f_0 * kk_1229[k];

        t_942[k] = -5.0 * hk_690[k]
                   + f_0 * kk_1230[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, t_946, t_947, hk_691, hk_692, hk_693, hk_694, \
                         hk_695, kk_1231, kk_1232, kk_1233, kk_1234, \
                         kk_1235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = -5.0 * hk_691[k]
                   + f_0 * kk_1231[k];

        t_944[k] = -5.0 * hk_692[k]
                   + f_0 * kk_1232[k];

        t_945[k] = -5.0 * hk_693[k]
                   + f_0 * kk_1233[k];

        t_946[k] = -5.0 * hk_694[k]
                   + f_0 * kk_1234[k];

        t_947[k] = -5.0 * hk_695[k]
                   + f_0 * kk_1235[k];
    }

#pragma omp simd aligned(t_948, t_949, t_950, t_951, t_952, hk_696, hk_697, hk_698, hk_699, \
                         hk_700, kk_1236, kk_1237, kk_1238, kk_1239, \
                         kk_1240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = -5.0 * hk_696[k]
                   + f_0 * kk_1236[k];

        t_949[k] = -5.0 * hk_697[k]
                   + f_0 * kk_1237[k];

        t_950[k] = -5.0 * hk_698[k]
                   + f_0 * kk_1238[k];

        t_951[k] = -5.0 * hk_699[k]
                   + f_0 * kk_1239[k];

        t_952[k] = -5.0 * hk_700[k]
                   + f_0 * kk_1240[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, t_957, hk_701, hk_702, hk_703, hk_704, \
                         hk_705, kk_1241, kk_1242, kk_1243, kk_1244, \
                         kk_1245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = -5.0 * hk_701[k]
                   + f_0 * kk_1241[k];

        t_954[k] = -5.0 * hk_702[k]
                   + f_0 * kk_1242[k];

        t_955[k] = -5.0 * hk_703[k]
                   + f_0 * kk_1243[k];

        t_956[k] = -5.0 * hk_704[k]
                   + f_0 * kk_1244[k];

        t_957[k] = -5.0 * hk_705[k]
                   + f_0 * kk_1245[k];
    }

#pragma omp simd aligned(t_958, t_959, t_960, t_961, t_962, hk_706, hk_707, hk_708, hk_709, \
                         hk_710, kk_1246, kk_1247, kk_1248, kk_1249, \
                         kk_1250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_958[k] = -5.0 * hk_706[k]
                   + f_0 * kk_1246[k];

        t_959[k] = -5.0 * hk_707[k]
                   + f_0 * kk_1247[k];

        t_960[k] = -5.0 * hk_708[k]
                   + f_0 * kk_1248[k];

        t_961[k] = -5.0 * hk_709[k]
                   + f_0 * kk_1249[k];

        t_962[k] = -5.0 * hk_710[k]
                   + f_0 * kk_1250[k];
    }

#pragma omp simd aligned(t_963, t_964, t_965, t_966, t_967, hk_711, hk_712, hk_713, hk_714, \
                         hk_715, kk_1251, kk_1252, kk_1253, kk_1254, \
                         kk_1255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_963[k] = -5.0 * hk_711[k]
                   + f_0 * kk_1251[k];

        t_964[k] = -5.0 * hk_712[k]
                   + f_0 * kk_1252[k];

        t_965[k] = -5.0 * hk_713[k]
                   + f_0 * kk_1253[k];

        t_966[k] = -5.0 * hk_714[k]
                   + f_0 * kk_1254[k];

        t_967[k] = -5.0 * hk_715[k]
                   + f_0 * kk_1255[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, t_972, hk_716, hk_717, hk_718, hk_719, \
                         hk_720, kk_1256, kk_1257, kk_1258, kk_1259, \
                         kk_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = -5.0 * hk_716[k]
                   + f_0 * kk_1256[k];

        t_969[k] = -5.0 * hk_717[k]
                   + f_0 * kk_1257[k];

        t_970[k] = -5.0 * hk_718[k]
                   + f_0 * kk_1258[k];

        t_971[k] = -5.0 * hk_719[k]
                   + f_0 * kk_1259[k];

        t_972[k] = -6.0 * hk_720[k]
                   + f_0 * kk_1260[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, t_976, t_977, hk_721, hk_722, hk_723, hk_724, \
                         hk_725, kk_1261, kk_1262, kk_1263, kk_1264, \
                         kk_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = -6.0 * hk_721[k]
                   + f_0 * kk_1261[k];

        t_974[k] = -6.0 * hk_722[k]
                   + f_0 * kk_1262[k];

        t_975[k] = -6.0 * hk_723[k]
                   + f_0 * kk_1263[k];

        t_976[k] = -6.0 * hk_724[k]
                   + f_0 * kk_1264[k];

        t_977[k] = -6.0 * hk_725[k]
                   + f_0 * kk_1265[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, t_982, hk_726, hk_727, hk_728, hk_729, \
                         hk_730, kk_1266, kk_1267, kk_1268, kk_1269, \
                         kk_1270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = -6.0 * hk_726[k]
                   + f_0 * kk_1266[k];

        t_979[k] = -6.0 * hk_727[k]
                   + f_0 * kk_1267[k];

        t_980[k] = -6.0 * hk_728[k]
                   + f_0 * kk_1268[k];

        t_981[k] = -6.0 * hk_729[k]
                   + f_0 * kk_1269[k];

        t_982[k] = -6.0 * hk_730[k]
                   + f_0 * kk_1270[k];
    }
}

static auto
compute_prim_geom_10_ik_electron_repulsion_2_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t kk,
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

    const auto *hk_731 = buffer.data(hk + 731);
    const auto *hk_732 = buffer.data(hk + 732);
    const auto *hk_733 = buffer.data(hk + 733);
    const auto *hk_734 = buffer.data(hk + 734);
    const auto *hk_735 = buffer.data(hk + 735);
    const auto *hk_736 = buffer.data(hk + 736);
    const auto *hk_737 = buffer.data(hk + 737);
    const auto *hk_738 = buffer.data(hk + 738);
    const auto *hk_739 = buffer.data(hk + 739);
    const auto *hk_740 = buffer.data(hk + 740);
    const auto *hk_741 = buffer.data(hk + 741);
    const auto *hk_742 = buffer.data(hk + 742);
    const auto *hk_743 = buffer.data(hk + 743);
    const auto *hk_744 = buffer.data(hk + 744);
    const auto *hk_745 = buffer.data(hk + 745);
    const auto *hk_746 = buffer.data(hk + 746);
    const auto *hk_747 = buffer.data(hk + 747);
    const auto *hk_748 = buffer.data(hk + 748);
    const auto *hk_749 = buffer.data(hk + 749);
    const auto *hk_750 = buffer.data(hk + 750);
    const auto *hk_751 = buffer.data(hk + 751);
    const auto *hk_752 = buffer.data(hk + 752);
    const auto *hk_753 = buffer.data(hk + 753);
    const auto *hk_754 = buffer.data(hk + 754);
    const auto *hk_755 = buffer.data(hk + 755);

    const auto *kk_1271 = buffer.data(kk + 1271);
    const auto *kk_1272 = buffer.data(kk + 1272);
    const auto *kk_1273 = buffer.data(kk + 1273);
    const auto *kk_1274 = buffer.data(kk + 1274);
    const auto *kk_1275 = buffer.data(kk + 1275);
    const auto *kk_1276 = buffer.data(kk + 1276);
    const auto *kk_1277 = buffer.data(kk + 1277);
    const auto *kk_1278 = buffer.data(kk + 1278);
    const auto *kk_1279 = buffer.data(kk + 1279);
    const auto *kk_1280 = buffer.data(kk + 1280);
    const auto *kk_1281 = buffer.data(kk + 1281);
    const auto *kk_1282 = buffer.data(kk + 1282);
    const auto *kk_1283 = buffer.data(kk + 1283);
    const auto *kk_1284 = buffer.data(kk + 1284);
    const auto *kk_1285 = buffer.data(kk + 1285);
    const auto *kk_1286 = buffer.data(kk + 1286);
    const auto *kk_1287 = buffer.data(kk + 1287);
    const auto *kk_1288 = buffer.data(kk + 1288);
    const auto *kk_1289 = buffer.data(kk + 1289);
    const auto *kk_1290 = buffer.data(kk + 1290);
    const auto *kk_1291 = buffer.data(kk + 1291);
    const auto *kk_1292 = buffer.data(kk + 1292);
    const auto *kk_1293 = buffer.data(kk + 1293);
    const auto *kk_1294 = buffer.data(kk + 1294);
    const auto *kk_1295 = buffer.data(kk + 1295);

#pragma omp simd aligned(t_983, t_984, t_985, t_986, t_987, hk_731, hk_732, hk_733, hk_734, \
                         hk_735, kk_1271, kk_1272, kk_1273, kk_1274, \
                         kk_1275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = -6.0 * hk_731[k]
                   + f_0 * kk_1271[k];

        t_984[k] = -6.0 * hk_732[k]
                   + f_0 * kk_1272[k];

        t_985[k] = -6.0 * hk_733[k]
                   + f_0 * kk_1273[k];

        t_986[k] = -6.0 * hk_734[k]
                   + f_0 * kk_1274[k];

        t_987[k] = -6.0 * hk_735[k]
                   + f_0 * kk_1275[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, t_991, t_992, hk_736, hk_737, hk_738, hk_739, \
                         hk_740, kk_1276, kk_1277, kk_1278, kk_1279, \
                         kk_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = -6.0 * hk_736[k]
                   + f_0 * kk_1276[k];

        t_989[k] = -6.0 * hk_737[k]
                   + f_0 * kk_1277[k];

        t_990[k] = -6.0 * hk_738[k]
                   + f_0 * kk_1278[k];

        t_991[k] = -6.0 * hk_739[k]
                   + f_0 * kk_1279[k];

        t_992[k] = -6.0 * hk_740[k]
                   + f_0 * kk_1280[k];
    }

#pragma omp simd aligned(t_993, t_994, t_995, t_996, t_997, hk_741, hk_742, hk_743, hk_744, \
                         hk_745, kk_1281, kk_1282, kk_1283, kk_1284, \
                         kk_1285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_993[k] = -6.0 * hk_741[k]
                   + f_0 * kk_1281[k];

        t_994[k] = -6.0 * hk_742[k]
                   + f_0 * kk_1282[k];

        t_995[k] = -6.0 * hk_743[k]
                   + f_0 * kk_1283[k];

        t_996[k] = -6.0 * hk_744[k]
                   + f_0 * kk_1284[k];

        t_997[k] = -6.0 * hk_745[k]
                   + f_0 * kk_1285[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, t_1001, t_1002, hk_746, hk_747, hk_748, hk_749, \
                         hk_750, kk_1286, kk_1287, kk_1288, kk_1289, \
                         kk_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = -6.0 * hk_746[k]
                   + f_0 * kk_1286[k];

        t_999[k] = -6.0 * hk_747[k]
                   + f_0 * kk_1287[k];

        t_1000[k] = -6.0 * hk_748[k]
                    + f_0 * kk_1288[k];

        t_1001[k] = -6.0 * hk_749[k]
                    + f_0 * kk_1289[k];

        t_1002[k] = -6.0 * hk_750[k]
                    + f_0 * kk_1290[k];
    }

#pragma omp simd aligned(t_1003, t_1004, t_1005, t_1006, t_1007, hk_751, hk_752, hk_753, \
                         hk_754, hk_755, kk_1291, kk_1292, kk_1293, kk_1294, \
                         kk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1003[k] = -6.0 * hk_751[k]
                    + f_0 * kk_1291[k];

        t_1004[k] = -6.0 * hk_752[k]
                    + f_0 * kk_1292[k];

        t_1005[k] = -6.0 * hk_753[k]
                    + f_0 * kk_1293[k];

        t_1006[k] = -6.0 * hk_754[k]
                    + f_0 * kk_1294[k];

        t_1007[k] = -6.0 * hk_755[k]
                    + f_0 * kk_1295[k];
    }
}

auto
compute_prim_geom_10_ik_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t hk, const size_t kk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ik_electron_repulsion_2_piece0(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_2_piece1(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_2_piece2(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_2_piece3(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_2_piece4(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_2_piece5(buffer, target, hk, kk, ncols, alpha);

    compute_prim_geom_10_ik_electron_repulsion_2_piece6(buffer, target, hk, kk, ncols, alpha);
}

}  // namespace simdt2ceri
