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


#include "SimdElectronRepulsionGeom10VrrRecGK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_gk_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t hk,
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

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);
    const auto *fk_36 = buffer.data(fk + 36);
    const auto *fk_37 = buffer.data(fk + 37);
    const auto *fk_38 = buffer.data(fk + 38);
    const auto *fk_39 = buffer.data(fk + 39);
    const auto *fk_40 = buffer.data(fk + 40);
    const auto *fk_41 = buffer.data(fk + 41);
    const auto *fk_42 = buffer.data(fk + 42);
    const auto *fk_43 = buffer.data(fk + 43);
    const auto *fk_44 = buffer.data(fk + 44);
    const auto *fk_45 = buffer.data(fk + 45);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_47 = buffer.data(fk + 47);
    const auto *fk_48 = buffer.data(fk + 48);
    const auto *fk_49 = buffer.data(fk + 49);
    const auto *fk_50 = buffer.data(fk + 50);
    const auto *fk_51 = buffer.data(fk + 51);
    const auto *fk_52 = buffer.data(fk + 52);
    const auto *fk_53 = buffer.data(fk + 53);
    const auto *fk_54 = buffer.data(fk + 54);
    const auto *fk_55 = buffer.data(fk + 55);
    const auto *fk_56 = buffer.data(fk + 56);
    const auto *fk_57 = buffer.data(fk + 57);
    const auto *fk_58 = buffer.data(fk + 58);
    const auto *fk_59 = buffer.data(fk + 59);
    const auto *fk_60 = buffer.data(fk + 60);
    const auto *fk_61 = buffer.data(fk + 61);
    const auto *fk_62 = buffer.data(fk + 62);
    const auto *fk_63 = buffer.data(fk + 63);
    const auto *fk_64 = buffer.data(fk + 64);
    const auto *fk_65 = buffer.data(fk + 65);
    const auto *fk_66 = buffer.data(fk + 66);
    const auto *fk_67 = buffer.data(fk + 67);
    const auto *fk_68 = buffer.data(fk + 68);
    const auto *fk_69 = buffer.data(fk + 69);
    const auto *fk_70 = buffer.data(fk + 70);
    const auto *fk_71 = buffer.data(fk + 71);
    const auto *fk_72 = buffer.data(fk + 72);
    const auto *fk_73 = buffer.data(fk + 73);
    const auto *fk_74 = buffer.data(fk + 74);
    const auto *fk_75 = buffer.data(fk + 75);
    const auto *fk_76 = buffer.data(fk + 76);
    const auto *fk_77 = buffer.data(fk + 77);
    const auto *fk_78 = buffer.data(fk + 78);
    const auto *fk_79 = buffer.data(fk + 79);
    const auto *fk_80 = buffer.data(fk + 80);
    const auto *fk_81 = buffer.data(fk + 81);
    const auto *fk_82 = buffer.data(fk + 82);
    const auto *fk_83 = buffer.data(fk + 83);
    const auto *fk_84 = buffer.data(fk + 84);
    const auto *fk_85 = buffer.data(fk + 85);
    const auto *fk_86 = buffer.data(fk + 86);
    const auto *fk_87 = buffer.data(fk + 87);
    const auto *fk_88 = buffer.data(fk + 88);
    const auto *fk_89 = buffer.data(fk + 89);
    const auto *fk_90 = buffer.data(fk + 90);
    const auto *fk_91 = buffer.data(fk + 91);
    const auto *fk_92 = buffer.data(fk + 92);
    const auto *fk_93 = buffer.data(fk + 93);
    const auto *fk_94 = buffer.data(fk + 94);
    const auto *fk_95 = buffer.data(fk + 95);
    const auto *fk_96 = buffer.data(fk + 96);
    const auto *fk_97 = buffer.data(fk + 97);
    const auto *fk_98 = buffer.data(fk + 98);
    const auto *fk_99 = buffer.data(fk + 99);
    const auto *fk_100 = buffer.data(fk + 100);
    const auto *fk_101 = buffer.data(fk + 101);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_106 = buffer.data(fk + 106);
    const auto *fk_107 = buffer.data(fk + 107);
    const auto *fk_108 = buffer.data(fk + 108);
    const auto *fk_109 = buffer.data(fk + 109);
    const auto *fk_110 = buffer.data(fk + 110);
    const auto *fk_111 = buffer.data(fk + 111);
    const auto *fk_112 = buffer.data(fk + 112);
    const auto *fk_113 = buffer.data(fk + 113);
    const auto *fk_114 = buffer.data(fk + 114);
    const auto *fk_115 = buffer.data(fk + 115);
    const auto *fk_116 = buffer.data(fk + 116);
    const auto *fk_117 = buffer.data(fk + 117);
    const auto *fk_118 = buffer.data(fk + 118);
    const auto *fk_119 = buffer.data(fk + 119);
    const auto *fk_120 = buffer.data(fk + 120);
    const auto *fk_121 = buffer.data(fk + 121);
    const auto *fk_122 = buffer.data(fk + 122);
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_124 = buffer.data(fk + 124);
    const auto *fk_125 = buffer.data(fk + 125);
    const auto *fk_126 = buffer.data(fk + 126);
    const auto *fk_127 = buffer.data(fk + 127);
    const auto *fk_128 = buffer.data(fk + 128);
    const auto *fk_129 = buffer.data(fk + 129);
    const auto *fk_130 = buffer.data(fk + 130);
    const auto *fk_131 = buffer.data(fk + 131);
    const auto *fk_132 = buffer.data(fk + 132);
    const auto *fk_133 = buffer.data(fk + 133);
    const auto *fk_134 = buffer.data(fk + 134);
    const auto *fk_135 = buffer.data(fk + 135);
    const auto *fk_136 = buffer.data(fk + 136);
    const auto *fk_137 = buffer.data(fk + 137);
    const auto *fk_138 = buffer.data(fk + 138);
    const auto *fk_139 = buffer.data(fk + 139);
    const auto *fk_140 = buffer.data(fk + 140);
    const auto *fk_141 = buffer.data(fk + 141);
    const auto *fk_142 = buffer.data(fk + 142);
    const auto *fk_143 = buffer.data(fk + 143);
    const auto *fk_144 = buffer.data(fk + 144);
    const auto *fk_145 = buffer.data(fk + 145);
    const auto *fk_146 = buffer.data(fk + 146);
    const auto *fk_147 = buffer.data(fk + 147);
    const auto *fk_148 = buffer.data(fk + 148);
    const auto *fk_149 = buffer.data(fk + 149);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, fk_0, fk_1, fk_2, fk_3, fk_4, hk_0, hk_1, \
                         hk_2, hk_3, hk_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -4.0 * fk_0[k]
                 + f_0 * hk_0[k];

        t_1[k] = -4.0 * fk_1[k]
                 + f_0 * hk_1[k];

        t_2[k] = -4.0 * fk_2[k]
                 + f_0 * hk_2[k];

        t_3[k] = -4.0 * fk_3[k]
                 + f_0 * hk_3[k];

        t_4[k] = -4.0 * fk_4[k]
                 + f_0 * hk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, fk_5, fk_6, fk_7, fk_8, fk_9, hk_5, hk_6, \
                         hk_7, hk_8, hk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -4.0 * fk_5[k]
                 + f_0 * hk_5[k];

        t_6[k] = -4.0 * fk_6[k]
                 + f_0 * hk_6[k];

        t_7[k] = -4.0 * fk_7[k]
                 + f_0 * hk_7[k];

        t_8[k] = -4.0 * fk_8[k]
                 + f_0 * hk_8[k];

        t_9[k] = -4.0 * fk_9[k]
                 + f_0 * hk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, fk_10, fk_11, fk_12, fk_13, fk_14, \
                         hk_10, hk_11, hk_12, hk_13, hk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -4.0 * fk_10[k]
                  + f_0 * hk_10[k];

        t_11[k] = -4.0 * fk_11[k]
                  + f_0 * hk_11[k];

        t_12[k] = -4.0 * fk_12[k]
                  + f_0 * hk_12[k];

        t_13[k] = -4.0 * fk_13[k]
                  + f_0 * hk_13[k];

        t_14[k] = -4.0 * fk_14[k]
                  + f_0 * hk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, fk_15, fk_16, fk_17, fk_18, fk_19, \
                         hk_15, hk_16, hk_17, hk_18, hk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -4.0 * fk_15[k]
                  + f_0 * hk_15[k];

        t_16[k] = -4.0 * fk_16[k]
                  + f_0 * hk_16[k];

        t_17[k] = -4.0 * fk_17[k]
                  + f_0 * hk_17[k];

        t_18[k] = -4.0 * fk_18[k]
                  + f_0 * hk_18[k];

        t_19[k] = -4.0 * fk_19[k]
                  + f_0 * hk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, fk_20, fk_21, fk_22, fk_23, fk_24, \
                         hk_20, hk_21, hk_22, hk_23, hk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -4.0 * fk_20[k]
                  + f_0 * hk_20[k];

        t_21[k] = -4.0 * fk_21[k]
                  + f_0 * hk_21[k];

        t_22[k] = -4.0 * fk_22[k]
                  + f_0 * hk_22[k];

        t_23[k] = -4.0 * fk_23[k]
                  + f_0 * hk_23[k];

        t_24[k] = -4.0 * fk_24[k]
                  + f_0 * hk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, fk_25, fk_26, fk_27, fk_28, fk_29, \
                         hk_25, hk_26, hk_27, hk_28, hk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -4.0 * fk_25[k]
                  + f_0 * hk_25[k];

        t_26[k] = -4.0 * fk_26[k]
                  + f_0 * hk_26[k];

        t_27[k] = -4.0 * fk_27[k]
                  + f_0 * hk_27[k];

        t_28[k] = -4.0 * fk_28[k]
                  + f_0 * hk_28[k];

        t_29[k] = -4.0 * fk_29[k]
                  + f_0 * hk_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, fk_30, fk_31, fk_32, fk_33, fk_34, \
                         hk_30, hk_31, hk_32, hk_33, hk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -4.0 * fk_30[k]
                  + f_0 * hk_30[k];

        t_31[k] = -4.0 * fk_31[k]
                  + f_0 * hk_31[k];

        t_32[k] = -4.0 * fk_32[k]
                  + f_0 * hk_32[k];

        t_33[k] = -4.0 * fk_33[k]
                  + f_0 * hk_33[k];

        t_34[k] = -4.0 * fk_34[k]
                  + f_0 * hk_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, fk_35, fk_36, fk_37, fk_38, fk_39, \
                         hk_35, hk_36, hk_37, hk_38, hk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -4.0 * fk_35[k]
                  + f_0 * hk_35[k];

        t_36[k] = -3.0 * fk_36[k]
                  + f_0 * hk_36[k];

        t_37[k] = -3.0 * fk_37[k]
                  + f_0 * hk_37[k];

        t_38[k] = -3.0 * fk_38[k]
                  + f_0 * hk_38[k];

        t_39[k] = -3.0 * fk_39[k]
                  + f_0 * hk_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, fk_40, fk_41, fk_42, fk_43, fk_44, \
                         hk_40, hk_41, hk_42, hk_43, hk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * fk_40[k]
                  + f_0 * hk_40[k];

        t_41[k] = -3.0 * fk_41[k]
                  + f_0 * hk_41[k];

        t_42[k] = -3.0 * fk_42[k]
                  + f_0 * hk_42[k];

        t_43[k] = -3.0 * fk_43[k]
                  + f_0 * hk_43[k];

        t_44[k] = -3.0 * fk_44[k]
                  + f_0 * hk_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, fk_45, fk_46, fk_47, fk_48, fk_49, \
                         hk_45, hk_46, hk_47, hk_48, hk_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -3.0 * fk_45[k]
                  + f_0 * hk_45[k];

        t_46[k] = -3.0 * fk_46[k]
                  + f_0 * hk_46[k];

        t_47[k] = -3.0 * fk_47[k]
                  + f_0 * hk_47[k];

        t_48[k] = -3.0 * fk_48[k]
                  + f_0 * hk_48[k];

        t_49[k] = -3.0 * fk_49[k]
                  + f_0 * hk_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, fk_50, fk_51, fk_52, fk_53, fk_54, \
                         hk_50, hk_51, hk_52, hk_53, hk_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -3.0 * fk_50[k]
                  + f_0 * hk_50[k];

        t_51[k] = -3.0 * fk_51[k]
                  + f_0 * hk_51[k];

        t_52[k] = -3.0 * fk_52[k]
                  + f_0 * hk_52[k];

        t_53[k] = -3.0 * fk_53[k]
                  + f_0 * hk_53[k];

        t_54[k] = -3.0 * fk_54[k]
                  + f_0 * hk_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, fk_55, fk_56, fk_57, fk_58, fk_59, \
                         hk_55, hk_56, hk_57, hk_58, hk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -3.0 * fk_55[k]
                  + f_0 * hk_55[k];

        t_56[k] = -3.0 * fk_56[k]
                  + f_0 * hk_56[k];

        t_57[k] = -3.0 * fk_57[k]
                  + f_0 * hk_57[k];

        t_58[k] = -3.0 * fk_58[k]
                  + f_0 * hk_58[k];

        t_59[k] = -3.0 * fk_59[k]
                  + f_0 * hk_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, fk_60, fk_61, fk_62, fk_63, fk_64, \
                         hk_60, hk_61, hk_62, hk_63, hk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * fk_60[k]
                  + f_0 * hk_60[k];

        t_61[k] = -3.0 * fk_61[k]
                  + f_0 * hk_61[k];

        t_62[k] = -3.0 * fk_62[k]
                  + f_0 * hk_62[k];

        t_63[k] = -3.0 * fk_63[k]
                  + f_0 * hk_63[k];

        t_64[k] = -3.0 * fk_64[k]
                  + f_0 * hk_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, fk_65, fk_66, fk_67, fk_68, fk_69, \
                         hk_65, hk_66, hk_67, hk_68, hk_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * fk_65[k]
                  + f_0 * hk_65[k];

        t_66[k] = -3.0 * fk_66[k]
                  + f_0 * hk_66[k];

        t_67[k] = -3.0 * fk_67[k]
                  + f_0 * hk_67[k];

        t_68[k] = -3.0 * fk_68[k]
                  + f_0 * hk_68[k];

        t_69[k] = -3.0 * fk_69[k]
                  + f_0 * hk_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, fk_70, fk_71, fk_72, fk_73, fk_74, \
                         hk_70, hk_71, hk_72, hk_73, hk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -3.0 * fk_70[k]
                  + f_0 * hk_70[k];

        t_71[k] = -3.0 * fk_71[k]
                  + f_0 * hk_71[k];

        t_72[k] = -3.0 * fk_72[k]
                  + f_0 * hk_72[k];

        t_73[k] = -3.0 * fk_73[k]
                  + f_0 * hk_73[k];

        t_74[k] = -3.0 * fk_74[k]
                  + f_0 * hk_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, fk_75, fk_76, fk_77, fk_78, fk_79, \
                         hk_75, hk_76, hk_77, hk_78, hk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -3.0 * fk_75[k]
                  + f_0 * hk_75[k];

        t_76[k] = -3.0 * fk_76[k]
                  + f_0 * hk_76[k];

        t_77[k] = -3.0 * fk_77[k]
                  + f_0 * hk_77[k];

        t_78[k] = -3.0 * fk_78[k]
                  + f_0 * hk_78[k];

        t_79[k] = -3.0 * fk_79[k]
                  + f_0 * hk_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, fk_80, fk_81, fk_82, fk_83, fk_84, \
                         hk_80, hk_81, hk_82, hk_83, hk_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * fk_80[k]
                  + f_0 * hk_80[k];

        t_81[k] = -3.0 * fk_81[k]
                  + f_0 * hk_81[k];

        t_82[k] = -3.0 * fk_82[k]
                  + f_0 * hk_82[k];

        t_83[k] = -3.0 * fk_83[k]
                  + f_0 * hk_83[k];

        t_84[k] = -3.0 * fk_84[k]
                  + f_0 * hk_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, fk_85, fk_86, fk_87, fk_88, fk_89, \
                         hk_85, hk_86, hk_87, hk_88, hk_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -3.0 * fk_85[k]
                  + f_0 * hk_85[k];

        t_86[k] = -3.0 * fk_86[k]
                  + f_0 * hk_86[k];

        t_87[k] = -3.0 * fk_87[k]
                  + f_0 * hk_87[k];

        t_88[k] = -3.0 * fk_88[k]
                  + f_0 * hk_88[k];

        t_89[k] = -3.0 * fk_89[k]
                  + f_0 * hk_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, fk_90, fk_91, fk_92, fk_93, fk_94, \
                         hk_90, hk_91, hk_92, hk_93, hk_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * fk_90[k]
                  + f_0 * hk_90[k];

        t_91[k] = -3.0 * fk_91[k]
                  + f_0 * hk_91[k];

        t_92[k] = -3.0 * fk_92[k]
                  + f_0 * hk_92[k];

        t_93[k] = -3.0 * fk_93[k]
                  + f_0 * hk_93[k];

        t_94[k] = -3.0 * fk_94[k]
                  + f_0 * hk_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, fk_95, fk_96, fk_97, fk_98, fk_99, \
                         hk_95, hk_96, hk_97, hk_98, hk_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * fk_95[k]
                  + f_0 * hk_95[k];

        t_96[k] = -3.0 * fk_96[k]
                  + f_0 * hk_96[k];

        t_97[k] = -3.0 * fk_97[k]
                  + f_0 * hk_97[k];

        t_98[k] = -3.0 * fk_98[k]
                  + f_0 * hk_98[k];

        t_99[k] = -3.0 * fk_99[k]
                  + f_0 * hk_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, fk_100, fk_101, fk_102, fk_103, \
                         fk_104, hk_100, hk_101, hk_102, hk_103, \
                         hk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * fk_100[k]
                   + f_0 * hk_100[k];

        t_101[k] = -3.0 * fk_101[k]
                   + f_0 * hk_101[k];

        t_102[k] = -3.0 * fk_102[k]
                   + f_0 * hk_102[k];

        t_103[k] = -3.0 * fk_103[k]
                   + f_0 * hk_103[k];

        t_104[k] = -3.0 * fk_104[k]
                   + f_0 * hk_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, fk_105, fk_106, fk_107, fk_108, \
                         fk_109, hk_105, hk_106, hk_107, hk_108, \
                         hk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -3.0 * fk_105[k]
                   + f_0 * hk_105[k];

        t_106[k] = -3.0 * fk_106[k]
                   + f_0 * hk_106[k];

        t_107[k] = -3.0 * fk_107[k]
                   + f_0 * hk_107[k];

        t_108[k] = -2.0 * fk_108[k]
                   + f_0 * hk_108[k];

        t_109[k] = -2.0 * fk_109[k]
                   + f_0 * hk_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, fk_110, fk_111, fk_112, fk_113, \
                         fk_114, hk_110, hk_111, hk_112, hk_113, \
                         hk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * fk_110[k]
                   + f_0 * hk_110[k];

        t_111[k] = -2.0 * fk_111[k]
                   + f_0 * hk_111[k];

        t_112[k] = -2.0 * fk_112[k]
                   + f_0 * hk_112[k];

        t_113[k] = -2.0 * fk_113[k]
                   + f_0 * hk_113[k];

        t_114[k] = -2.0 * fk_114[k]
                   + f_0 * hk_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, fk_115, fk_116, fk_117, fk_118, \
                         fk_119, hk_115, hk_116, hk_117, hk_118, \
                         hk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * fk_115[k]
                   + f_0 * hk_115[k];

        t_116[k] = -2.0 * fk_116[k]
                   + f_0 * hk_116[k];

        t_117[k] = -2.0 * fk_117[k]
                   + f_0 * hk_117[k];

        t_118[k] = -2.0 * fk_118[k]
                   + f_0 * hk_118[k];

        t_119[k] = -2.0 * fk_119[k]
                   + f_0 * hk_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, fk_120, fk_121, fk_122, fk_123, \
                         fk_124, hk_120, hk_121, hk_122, hk_123, \
                         hk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -2.0 * fk_120[k]
                   + f_0 * hk_120[k];

        t_121[k] = -2.0 * fk_121[k]
                   + f_0 * hk_121[k];

        t_122[k] = -2.0 * fk_122[k]
                   + f_0 * hk_122[k];

        t_123[k] = -2.0 * fk_123[k]
                   + f_0 * hk_123[k];

        t_124[k] = -2.0 * fk_124[k]
                   + f_0 * hk_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, fk_125, fk_126, fk_127, fk_128, \
                         fk_129, hk_125, hk_126, hk_127, hk_128, \
                         hk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -2.0 * fk_125[k]
                   + f_0 * hk_125[k];

        t_126[k] = -2.0 * fk_126[k]
                   + f_0 * hk_126[k];

        t_127[k] = -2.0 * fk_127[k]
                   + f_0 * hk_127[k];

        t_128[k] = -2.0 * fk_128[k]
                   + f_0 * hk_128[k];

        t_129[k] = -2.0 * fk_129[k]
                   + f_0 * hk_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, fk_130, fk_131, fk_132, fk_133, \
                         fk_134, hk_130, hk_131, hk_132, hk_133, \
                         hk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -2.0 * fk_130[k]
                   + f_0 * hk_130[k];

        t_131[k] = -2.0 * fk_131[k]
                   + f_0 * hk_131[k];

        t_132[k] = -2.0 * fk_132[k]
                   + f_0 * hk_132[k];

        t_133[k] = -2.0 * fk_133[k]
                   + f_0 * hk_133[k];

        t_134[k] = -2.0 * fk_134[k]
                   + f_0 * hk_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, fk_135, fk_136, fk_137, fk_138, \
                         fk_139, hk_135, hk_136, hk_137, hk_138, \
                         hk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -2.0 * fk_135[k]
                   + f_0 * hk_135[k];

        t_136[k] = -2.0 * fk_136[k]
                   + f_0 * hk_136[k];

        t_137[k] = -2.0 * fk_137[k]
                   + f_0 * hk_137[k];

        t_138[k] = -2.0 * fk_138[k]
                   + f_0 * hk_138[k];

        t_139[k] = -2.0 * fk_139[k]
                   + f_0 * hk_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, fk_140, fk_141, fk_142, fk_143, \
                         fk_144, hk_140, hk_141, hk_142, hk_143, \
                         hk_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -2.0 * fk_140[k]
                   + f_0 * hk_140[k];

        t_141[k] = -2.0 * fk_141[k]
                   + f_0 * hk_141[k];

        t_142[k] = -2.0 * fk_142[k]
                   + f_0 * hk_142[k];

        t_143[k] = -2.0 * fk_143[k]
                   + f_0 * hk_143[k];

        t_144[k] = -2.0 * fk_144[k]
                   + f_0 * hk_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, fk_145, fk_146, fk_147, fk_148, \
                         fk_149, hk_145, hk_146, hk_147, hk_148, \
                         hk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -2.0 * fk_145[k]
                   + f_0 * hk_145[k];

        t_146[k] = -2.0 * fk_146[k]
                   + f_0 * hk_146[k];

        t_147[k] = -2.0 * fk_147[k]
                   + f_0 * hk_147[k];

        t_148[k] = -2.0 * fk_148[k]
                   + f_0 * hk_148[k];

        t_149[k] = -2.0 * fk_149[k]
                   + f_0 * hk_149[k];
    }
}

static auto
compute_prim_geom_10_gk_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t hk,
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

    const auto *fk_150 = buffer.data(fk + 150);
    const auto *fk_151 = buffer.data(fk + 151);
    const auto *fk_152 = buffer.data(fk + 152);
    const auto *fk_153 = buffer.data(fk + 153);
    const auto *fk_154 = buffer.data(fk + 154);
    const auto *fk_155 = buffer.data(fk + 155);
    const auto *fk_156 = buffer.data(fk + 156);
    const auto *fk_157 = buffer.data(fk + 157);
    const auto *fk_158 = buffer.data(fk + 158);
    const auto *fk_159 = buffer.data(fk + 159);
    const auto *fk_160 = buffer.data(fk + 160);
    const auto *fk_161 = buffer.data(fk + 161);
    const auto *fk_162 = buffer.data(fk + 162);
    const auto *fk_163 = buffer.data(fk + 163);
    const auto *fk_164 = buffer.data(fk + 164);
    const auto *fk_165 = buffer.data(fk + 165);
    const auto *fk_166 = buffer.data(fk + 166);
    const auto *fk_167 = buffer.data(fk + 167);
    const auto *fk_168 = buffer.data(fk + 168);
    const auto *fk_169 = buffer.data(fk + 169);
    const auto *fk_170 = buffer.data(fk + 170);
    const auto *fk_171 = buffer.data(fk + 171);
    const auto *fk_172 = buffer.data(fk + 172);
    const auto *fk_173 = buffer.data(fk + 173);
    const auto *fk_174 = buffer.data(fk + 174);
    const auto *fk_175 = buffer.data(fk + 175);
    const auto *fk_176 = buffer.data(fk + 176);
    const auto *fk_177 = buffer.data(fk + 177);
    const auto *fk_178 = buffer.data(fk + 178);
    const auto *fk_179 = buffer.data(fk + 179);
    const auto *fk_180 = buffer.data(fk + 180);
    const auto *fk_181 = buffer.data(fk + 181);
    const auto *fk_182 = buffer.data(fk + 182);
    const auto *fk_183 = buffer.data(fk + 183);
    const auto *fk_184 = buffer.data(fk + 184);
    const auto *fk_185 = buffer.data(fk + 185);
    const auto *fk_186 = buffer.data(fk + 186);
    const auto *fk_187 = buffer.data(fk + 187);
    const auto *fk_188 = buffer.data(fk + 188);
    const auto *fk_189 = buffer.data(fk + 189);
    const auto *fk_190 = buffer.data(fk + 190);
    const auto *fk_191 = buffer.data(fk + 191);
    const auto *fk_192 = buffer.data(fk + 192);
    const auto *fk_193 = buffer.data(fk + 193);
    const auto *fk_194 = buffer.data(fk + 194);
    const auto *fk_195 = buffer.data(fk + 195);
    const auto *fk_196 = buffer.data(fk + 196);
    const auto *fk_197 = buffer.data(fk + 197);
    const auto *fk_198 = buffer.data(fk + 198);
    const auto *fk_199 = buffer.data(fk + 199);
    const auto *fk_200 = buffer.data(fk + 200);
    const auto *fk_201 = buffer.data(fk + 201);
    const auto *fk_202 = buffer.data(fk + 202);
    const auto *fk_203 = buffer.data(fk + 203);
    const auto *fk_204 = buffer.data(fk + 204);
    const auto *fk_205 = buffer.data(fk + 205);
    const auto *fk_206 = buffer.data(fk + 206);
    const auto *fk_207 = buffer.data(fk + 207);
    const auto *fk_208 = buffer.data(fk + 208);
    const auto *fk_209 = buffer.data(fk + 209);
    const auto *fk_210 = buffer.data(fk + 210);
    const auto *fk_211 = buffer.data(fk + 211);
    const auto *fk_212 = buffer.data(fk + 212);
    const auto *fk_213 = buffer.data(fk + 213);
    const auto *fk_214 = buffer.data(fk + 214);
    const auto *fk_215 = buffer.data(fk + 215);
    const auto *fk_216 = buffer.data(fk + 216);
    const auto *fk_217 = buffer.data(fk + 217);
    const auto *fk_218 = buffer.data(fk + 218);
    const auto *fk_219 = buffer.data(fk + 219);
    const auto *fk_220 = buffer.data(fk + 220);
    const auto *fk_221 = buffer.data(fk + 221);
    const auto *fk_222 = buffer.data(fk + 222);
    const auto *fk_223 = buffer.data(fk + 223);
    const auto *fk_224 = buffer.data(fk + 224);
    const auto *fk_225 = buffer.data(fk + 225);
    const auto *fk_226 = buffer.data(fk + 226);
    const auto *fk_227 = buffer.data(fk + 227);
    const auto *fk_228 = buffer.data(fk + 228);
    const auto *fk_229 = buffer.data(fk + 229);
    const auto *fk_230 = buffer.data(fk + 230);
    const auto *fk_231 = buffer.data(fk + 231);
    const auto *fk_232 = buffer.data(fk + 232);
    const auto *fk_233 = buffer.data(fk + 233);
    const auto *fk_234 = buffer.data(fk + 234);
    const auto *fk_235 = buffer.data(fk + 235);
    const auto *fk_236 = buffer.data(fk + 236);
    const auto *fk_237 = buffer.data(fk + 237);
    const auto *fk_238 = buffer.data(fk + 238);
    const auto *fk_239 = buffer.data(fk + 239);
    const auto *fk_240 = buffer.data(fk + 240);
    const auto *fk_241 = buffer.data(fk + 241);
    const auto *fk_242 = buffer.data(fk + 242);
    const auto *fk_243 = buffer.data(fk + 243);
    const auto *fk_244 = buffer.data(fk + 244);
    const auto *fk_245 = buffer.data(fk + 245);
    const auto *fk_246 = buffer.data(fk + 246);
    const auto *fk_247 = buffer.data(fk + 247);
    const auto *fk_248 = buffer.data(fk + 248);
    const auto *fk_249 = buffer.data(fk + 249);
    const auto *fk_250 = buffer.data(fk + 250);
    const auto *fk_251 = buffer.data(fk + 251);
    const auto *fk_252 = buffer.data(fk + 252);
    const auto *fk_253 = buffer.data(fk + 253);
    const auto *fk_254 = buffer.data(fk + 254);
    const auto *fk_255 = buffer.data(fk + 255);
    const auto *fk_256 = buffer.data(fk + 256);
    const auto *fk_257 = buffer.data(fk + 257);
    const auto *fk_258 = buffer.data(fk + 258);
    const auto *fk_259 = buffer.data(fk + 259);
    const auto *fk_260 = buffer.data(fk + 260);
    const auto *fk_261 = buffer.data(fk + 261);
    const auto *fk_262 = buffer.data(fk + 262);
    const auto *fk_263 = buffer.data(fk + 263);
    const auto *fk_264 = buffer.data(fk + 264);
    const auto *fk_265 = buffer.data(fk + 265);
    const auto *fk_266 = buffer.data(fk + 266);
    const auto *fk_267 = buffer.data(fk + 267);
    const auto *fk_268 = buffer.data(fk + 268);
    const auto *fk_269 = buffer.data(fk + 269);
    const auto *fk_270 = buffer.data(fk + 270);
    const auto *fk_271 = buffer.data(fk + 271);
    const auto *fk_272 = buffer.data(fk + 272);
    const auto *fk_273 = buffer.data(fk + 273);
    const auto *fk_274 = buffer.data(fk + 274);
    const auto *fk_275 = buffer.data(fk + 275);
    const auto *fk_276 = buffer.data(fk + 276);
    const auto *fk_277 = buffer.data(fk + 277);
    const auto *fk_278 = buffer.data(fk + 278);
    const auto *fk_279 = buffer.data(fk + 279);
    const auto *fk_280 = buffer.data(fk + 280);
    const auto *fk_281 = buffer.data(fk + 281);
    const auto *fk_282 = buffer.data(fk + 282);
    const auto *fk_283 = buffer.data(fk + 283);
    const auto *fk_284 = buffer.data(fk + 284);
    const auto *fk_285 = buffer.data(fk + 285);
    const auto *fk_286 = buffer.data(fk + 286);
    const auto *fk_287 = buffer.data(fk + 287);
    const auto *fk_288 = buffer.data(fk + 288);
    const auto *fk_289 = buffer.data(fk + 289);
    const auto *fk_290 = buffer.data(fk + 290);
    const auto *fk_291 = buffer.data(fk + 291);
    const auto *fk_292 = buffer.data(fk + 292);
    const auto *fk_293 = buffer.data(fk + 293);
    const auto *fk_294 = buffer.data(fk + 294);
    const auto *fk_295 = buffer.data(fk + 295);
    const auto *fk_296 = buffer.data(fk + 296);
    const auto *fk_297 = buffer.data(fk + 297);
    const auto *fk_298 = buffer.data(fk + 298);
    const auto *fk_299 = buffer.data(fk + 299);

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

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, fk_150, fk_151, fk_152, fk_153, \
                         fk_154, hk_150, hk_151, hk_152, hk_153, \
                         hk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -2.0 * fk_150[k]
                   + f_0 * hk_150[k];

        t_151[k] = -2.0 * fk_151[k]
                   + f_0 * hk_151[k];

        t_152[k] = -2.0 * fk_152[k]
                   + f_0 * hk_152[k];

        t_153[k] = -2.0 * fk_153[k]
                   + f_0 * hk_153[k];

        t_154[k] = -2.0 * fk_154[k]
                   + f_0 * hk_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, fk_155, fk_156, fk_157, fk_158, \
                         fk_159, hk_155, hk_156, hk_157, hk_158, \
                         hk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -2.0 * fk_155[k]
                   + f_0 * hk_155[k];

        t_156[k] = -2.0 * fk_156[k]
                   + f_0 * hk_156[k];

        t_157[k] = -2.0 * fk_157[k]
                   + f_0 * hk_157[k];

        t_158[k] = -2.0 * fk_158[k]
                   + f_0 * hk_158[k];

        t_159[k] = -2.0 * fk_159[k]
                   + f_0 * hk_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, fk_160, fk_161, fk_162, fk_163, \
                         fk_164, hk_160, hk_161, hk_162, hk_163, \
                         hk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -2.0 * fk_160[k]
                   + f_0 * hk_160[k];

        t_161[k] = -2.0 * fk_161[k]
                   + f_0 * hk_161[k];

        t_162[k] = -2.0 * fk_162[k]
                   + f_0 * hk_162[k];

        t_163[k] = -2.0 * fk_163[k]
                   + f_0 * hk_163[k];

        t_164[k] = -2.0 * fk_164[k]
                   + f_0 * hk_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, fk_165, fk_166, fk_167, fk_168, \
                         fk_169, hk_165, hk_166, hk_167, hk_168, \
                         hk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -2.0 * fk_165[k]
                   + f_0 * hk_165[k];

        t_166[k] = -2.0 * fk_166[k]
                   + f_0 * hk_166[k];

        t_167[k] = -2.0 * fk_167[k]
                   + f_0 * hk_167[k];

        t_168[k] = -2.0 * fk_168[k]
                   + f_0 * hk_168[k];

        t_169[k] = -2.0 * fk_169[k]
                   + f_0 * hk_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, fk_170, fk_171, fk_172, fk_173, \
                         fk_174, hk_170, hk_171, hk_172, hk_173, \
                         hk_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -2.0 * fk_170[k]
                   + f_0 * hk_170[k];

        t_171[k] = -2.0 * fk_171[k]
                   + f_0 * hk_171[k];

        t_172[k] = -2.0 * fk_172[k]
                   + f_0 * hk_172[k];

        t_173[k] = -2.0 * fk_173[k]
                   + f_0 * hk_173[k];

        t_174[k] = -2.0 * fk_174[k]
                   + f_0 * hk_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, fk_175, fk_176, fk_177, fk_178, \
                         fk_179, hk_175, hk_176, hk_177, hk_178, \
                         hk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -2.0 * fk_175[k]
                   + f_0 * hk_175[k];

        t_176[k] = -2.0 * fk_176[k]
                   + f_0 * hk_176[k];

        t_177[k] = -2.0 * fk_177[k]
                   + f_0 * hk_177[k];

        t_178[k] = -2.0 * fk_178[k]
                   + f_0 * hk_178[k];

        t_179[k] = -2.0 * fk_179[k]
                   + f_0 * hk_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, fk_180, fk_181, fk_182, fk_183, \
                         fk_184, hk_180, hk_181, hk_182, hk_183, \
                         hk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * fk_180[k]
                   + f_0 * hk_180[k];

        t_181[k] = -2.0 * fk_181[k]
                   + f_0 * hk_181[k];

        t_182[k] = -2.0 * fk_182[k]
                   + f_0 * hk_182[k];

        t_183[k] = -2.0 * fk_183[k]
                   + f_0 * hk_183[k];

        t_184[k] = -2.0 * fk_184[k]
                   + f_0 * hk_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, fk_185, fk_186, fk_187, fk_188, \
                         fk_189, hk_185, hk_186, hk_187, hk_188, \
                         hk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * fk_185[k]
                   + f_0 * hk_185[k];

        t_186[k] = -2.0 * fk_186[k]
                   + f_0 * hk_186[k];

        t_187[k] = -2.0 * fk_187[k]
                   + f_0 * hk_187[k];

        t_188[k] = -2.0 * fk_188[k]
                   + f_0 * hk_188[k];

        t_189[k] = -2.0 * fk_189[k]
                   + f_0 * hk_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, fk_190, fk_191, fk_192, fk_193, \
                         fk_194, hk_190, hk_191, hk_192, hk_193, \
                         hk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * fk_190[k]
                   + f_0 * hk_190[k];

        t_191[k] = -2.0 * fk_191[k]
                   + f_0 * hk_191[k];

        t_192[k] = -2.0 * fk_192[k]
                   + f_0 * hk_192[k];

        t_193[k] = -2.0 * fk_193[k]
                   + f_0 * hk_193[k];

        t_194[k] = -2.0 * fk_194[k]
                   + f_0 * hk_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, fk_195, fk_196, fk_197, fk_198, \
                         fk_199, hk_195, hk_196, hk_197, hk_198, \
                         hk_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * fk_195[k]
                   + f_0 * hk_195[k];

        t_196[k] = -2.0 * fk_196[k]
                   + f_0 * hk_196[k];

        t_197[k] = -2.0 * fk_197[k]
                   + f_0 * hk_197[k];

        t_198[k] = -2.0 * fk_198[k]
                   + f_0 * hk_198[k];

        t_199[k] = -2.0 * fk_199[k]
                   + f_0 * hk_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, fk_200, fk_201, fk_202, fk_203, \
                         fk_204, hk_200, hk_201, hk_202, hk_203, \
                         hk_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * fk_200[k]
                   + f_0 * hk_200[k];

        t_201[k] = -2.0 * fk_201[k]
                   + f_0 * hk_201[k];

        t_202[k] = -2.0 * fk_202[k]
                   + f_0 * hk_202[k];

        t_203[k] = -2.0 * fk_203[k]
                   + f_0 * hk_203[k];

        t_204[k] = -2.0 * fk_204[k]
                   + f_0 * hk_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, fk_205, fk_206, fk_207, fk_208, \
                         fk_209, hk_205, hk_206, hk_207, hk_208, \
                         hk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * fk_205[k]
                   + f_0 * hk_205[k];

        t_206[k] = -2.0 * fk_206[k]
                   + f_0 * hk_206[k];

        t_207[k] = -2.0 * fk_207[k]
                   + f_0 * hk_207[k];

        t_208[k] = -2.0 * fk_208[k]
                   + f_0 * hk_208[k];

        t_209[k] = -2.0 * fk_209[k]
                   + f_0 * hk_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, fk_210, fk_211, fk_212, fk_213, \
                         fk_214, hk_210, hk_211, hk_212, hk_213, \
                         hk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * fk_210[k]
                   + f_0 * hk_210[k];

        t_211[k] = -2.0 * fk_211[k]
                   + f_0 * hk_211[k];

        t_212[k] = -2.0 * fk_212[k]
                   + f_0 * hk_212[k];

        t_213[k] = -2.0 * fk_213[k]
                   + f_0 * hk_213[k];

        t_214[k] = -2.0 * fk_214[k]
                   + f_0 * hk_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, fk_215, fk_216, fk_217, fk_218, \
                         fk_219, hk_215, hk_216, hk_217, hk_218, \
                         hk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * fk_215[k]
                   + f_0 * hk_215[k];

        t_216[k] = -fk_216[k]
                   + f_0 * hk_216[k];

        t_217[k] = -fk_217[k]
                   + f_0 * hk_217[k];

        t_218[k] = -fk_218[k]
                   + f_0 * hk_218[k];

        t_219[k] = -fk_219[k]
                   + f_0 * hk_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, fk_220, fk_221, fk_222, fk_223, \
                         fk_224, hk_220, hk_221, hk_222, hk_223, \
                         hk_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -fk_220[k]
                   + f_0 * hk_220[k];

        t_221[k] = -fk_221[k]
                   + f_0 * hk_221[k];

        t_222[k] = -fk_222[k]
                   + f_0 * hk_222[k];

        t_223[k] = -fk_223[k]
                   + f_0 * hk_223[k];

        t_224[k] = -fk_224[k]
                   + f_0 * hk_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, fk_225, fk_226, fk_227, fk_228, \
                         fk_229, hk_225, hk_226, hk_227, hk_228, \
                         hk_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -fk_225[k]
                   + f_0 * hk_225[k];

        t_226[k] = -fk_226[k]
                   + f_0 * hk_226[k];

        t_227[k] = -fk_227[k]
                   + f_0 * hk_227[k];

        t_228[k] = -fk_228[k]
                   + f_0 * hk_228[k];

        t_229[k] = -fk_229[k]
                   + f_0 * hk_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, fk_230, fk_231, fk_232, fk_233, \
                         fk_234, hk_230, hk_231, hk_232, hk_233, \
                         hk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -fk_230[k]
                   + f_0 * hk_230[k];

        t_231[k] = -fk_231[k]
                   + f_0 * hk_231[k];

        t_232[k] = -fk_232[k]
                   + f_0 * hk_232[k];

        t_233[k] = -fk_233[k]
                   + f_0 * hk_233[k];

        t_234[k] = -fk_234[k]
                   + f_0 * hk_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, fk_235, fk_236, fk_237, fk_238, \
                         fk_239, hk_235, hk_236, hk_237, hk_238, \
                         hk_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -fk_235[k]
                   + f_0 * hk_235[k];

        t_236[k] = -fk_236[k]
                   + f_0 * hk_236[k];

        t_237[k] = -fk_237[k]
                   + f_0 * hk_237[k];

        t_238[k] = -fk_238[k]
                   + f_0 * hk_238[k];

        t_239[k] = -fk_239[k]
                   + f_0 * hk_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, fk_240, fk_241, fk_242, fk_243, \
                         fk_244, hk_240, hk_241, hk_242, hk_243, \
                         hk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -fk_240[k]
                   + f_0 * hk_240[k];

        t_241[k] = -fk_241[k]
                   + f_0 * hk_241[k];

        t_242[k] = -fk_242[k]
                   + f_0 * hk_242[k];

        t_243[k] = -fk_243[k]
                   + f_0 * hk_243[k];

        t_244[k] = -fk_244[k]
                   + f_0 * hk_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, fk_245, fk_246, fk_247, fk_248, \
                         fk_249, hk_245, hk_246, hk_247, hk_248, \
                         hk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -fk_245[k]
                   + f_0 * hk_245[k];

        t_246[k] = -fk_246[k]
                   + f_0 * hk_246[k];

        t_247[k] = -fk_247[k]
                   + f_0 * hk_247[k];

        t_248[k] = -fk_248[k]
                   + f_0 * hk_248[k];

        t_249[k] = -fk_249[k]
                   + f_0 * hk_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, fk_250, fk_251, fk_252, fk_253, \
                         fk_254, hk_250, hk_251, hk_252, hk_253, \
                         hk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -fk_250[k]
                   + f_0 * hk_250[k];

        t_251[k] = -fk_251[k]
                   + f_0 * hk_251[k];

        t_252[k] = -fk_252[k]
                   + f_0 * hk_252[k];

        t_253[k] = -fk_253[k]
                   + f_0 * hk_253[k];

        t_254[k] = -fk_254[k]
                   + f_0 * hk_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, fk_255, fk_256, fk_257, fk_258, \
                         fk_259, hk_255, hk_256, hk_257, hk_258, \
                         hk_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -fk_255[k]
                   + f_0 * hk_255[k];

        t_256[k] = -fk_256[k]
                   + f_0 * hk_256[k];

        t_257[k] = -fk_257[k]
                   + f_0 * hk_257[k];

        t_258[k] = -fk_258[k]
                   + f_0 * hk_258[k];

        t_259[k] = -fk_259[k]
                   + f_0 * hk_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, fk_260, fk_261, fk_262, fk_263, \
                         fk_264, hk_260, hk_261, hk_262, hk_263, \
                         hk_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -fk_260[k]
                   + f_0 * hk_260[k];

        t_261[k] = -fk_261[k]
                   + f_0 * hk_261[k];

        t_262[k] = -fk_262[k]
                   + f_0 * hk_262[k];

        t_263[k] = -fk_263[k]
                   + f_0 * hk_263[k];

        t_264[k] = -fk_264[k]
                   + f_0 * hk_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, fk_265, fk_266, fk_267, fk_268, \
                         fk_269, hk_265, hk_266, hk_267, hk_268, \
                         hk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -fk_265[k]
                   + f_0 * hk_265[k];

        t_266[k] = -fk_266[k]
                   + f_0 * hk_266[k];

        t_267[k] = -fk_267[k]
                   + f_0 * hk_267[k];

        t_268[k] = -fk_268[k]
                   + f_0 * hk_268[k];

        t_269[k] = -fk_269[k]
                   + f_0 * hk_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, fk_270, fk_271, fk_272, fk_273, \
                         fk_274, hk_270, hk_271, hk_272, hk_273, \
                         hk_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -fk_270[k]
                   + f_0 * hk_270[k];

        t_271[k] = -fk_271[k]
                   + f_0 * hk_271[k];

        t_272[k] = -fk_272[k]
                   + f_0 * hk_272[k];

        t_273[k] = -fk_273[k]
                   + f_0 * hk_273[k];

        t_274[k] = -fk_274[k]
                   + f_0 * hk_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, fk_275, fk_276, fk_277, fk_278, \
                         fk_279, hk_275, hk_276, hk_277, hk_278, \
                         hk_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -fk_275[k]
                   + f_0 * hk_275[k];

        t_276[k] = -fk_276[k]
                   + f_0 * hk_276[k];

        t_277[k] = -fk_277[k]
                   + f_0 * hk_277[k];

        t_278[k] = -fk_278[k]
                   + f_0 * hk_278[k];

        t_279[k] = -fk_279[k]
                   + f_0 * hk_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, fk_280, fk_281, fk_282, fk_283, \
                         fk_284, hk_280, hk_281, hk_282, hk_283, \
                         hk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -fk_280[k]
                   + f_0 * hk_280[k];

        t_281[k] = -fk_281[k]
                   + f_0 * hk_281[k];

        t_282[k] = -fk_282[k]
                   + f_0 * hk_282[k];

        t_283[k] = -fk_283[k]
                   + f_0 * hk_283[k];

        t_284[k] = -fk_284[k]
                   + f_0 * hk_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, fk_285, fk_286, fk_287, fk_288, \
                         fk_289, hk_285, hk_286, hk_287, hk_288, \
                         hk_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -fk_285[k]
                   + f_0 * hk_285[k];

        t_286[k] = -fk_286[k]
                   + f_0 * hk_286[k];

        t_287[k] = -fk_287[k]
                   + f_0 * hk_287[k];

        t_288[k] = -fk_288[k]
                   + f_0 * hk_288[k];

        t_289[k] = -fk_289[k]
                   + f_0 * hk_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, fk_290, fk_291, fk_292, fk_293, \
                         fk_294, hk_290, hk_291, hk_292, hk_293, \
                         hk_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -fk_290[k]
                   + f_0 * hk_290[k];

        t_291[k] = -fk_291[k]
                   + f_0 * hk_291[k];

        t_292[k] = -fk_292[k]
                   + f_0 * hk_292[k];

        t_293[k] = -fk_293[k]
                   + f_0 * hk_293[k];

        t_294[k] = -fk_294[k]
                   + f_0 * hk_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, fk_295, fk_296, fk_297, fk_298, \
                         fk_299, hk_295, hk_296, hk_297, hk_298, \
                         hk_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -fk_295[k]
                   + f_0 * hk_295[k];

        t_296[k] = -fk_296[k]
                   + f_0 * hk_296[k];

        t_297[k] = -fk_297[k]
                   + f_0 * hk_297[k];

        t_298[k] = -fk_298[k]
                   + f_0 * hk_298[k];

        t_299[k] = -fk_299[k]
                   + f_0 * hk_299[k];
    }
}

static auto
compute_prim_geom_10_gk_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t hk,
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

    const auto *fk_300 = buffer.data(fk + 300);
    const auto *fk_301 = buffer.data(fk + 301);
    const auto *fk_302 = buffer.data(fk + 302);
    const auto *fk_303 = buffer.data(fk + 303);
    const auto *fk_304 = buffer.data(fk + 304);
    const auto *fk_305 = buffer.data(fk + 305);
    const auto *fk_306 = buffer.data(fk + 306);
    const auto *fk_307 = buffer.data(fk + 307);
    const auto *fk_308 = buffer.data(fk + 308);
    const auto *fk_309 = buffer.data(fk + 309);
    const auto *fk_310 = buffer.data(fk + 310);
    const auto *fk_311 = buffer.data(fk + 311);
    const auto *fk_312 = buffer.data(fk + 312);
    const auto *fk_313 = buffer.data(fk + 313);
    const auto *fk_314 = buffer.data(fk + 314);
    const auto *fk_315 = buffer.data(fk + 315);
    const auto *fk_316 = buffer.data(fk + 316);
    const auto *fk_317 = buffer.data(fk + 317);
    const auto *fk_318 = buffer.data(fk + 318);
    const auto *fk_319 = buffer.data(fk + 319);
    const auto *fk_320 = buffer.data(fk + 320);
    const auto *fk_321 = buffer.data(fk + 321);
    const auto *fk_322 = buffer.data(fk + 322);
    const auto *fk_323 = buffer.data(fk + 323);
    const auto *fk_324 = buffer.data(fk + 324);
    const auto *fk_325 = buffer.data(fk + 325);
    const auto *fk_326 = buffer.data(fk + 326);
    const auto *fk_327 = buffer.data(fk + 327);
    const auto *fk_328 = buffer.data(fk + 328);
    const auto *fk_329 = buffer.data(fk + 329);
    const auto *fk_330 = buffer.data(fk + 330);
    const auto *fk_331 = buffer.data(fk + 331);
    const auto *fk_332 = buffer.data(fk + 332);
    const auto *fk_333 = buffer.data(fk + 333);
    const auto *fk_334 = buffer.data(fk + 334);
    const auto *fk_335 = buffer.data(fk + 335);
    const auto *fk_336 = buffer.data(fk + 336);
    const auto *fk_337 = buffer.data(fk + 337);
    const auto *fk_338 = buffer.data(fk + 338);
    const auto *fk_339 = buffer.data(fk + 339);
    const auto *fk_340 = buffer.data(fk + 340);
    const auto *fk_341 = buffer.data(fk + 341);
    const auto *fk_342 = buffer.data(fk + 342);
    const auto *fk_343 = buffer.data(fk + 343);
    const auto *fk_344 = buffer.data(fk + 344);
    const auto *fk_345 = buffer.data(fk + 345);
    const auto *fk_346 = buffer.data(fk + 346);
    const auto *fk_347 = buffer.data(fk + 347);
    const auto *fk_348 = buffer.data(fk + 348);
    const auto *fk_349 = buffer.data(fk + 349);
    const auto *fk_350 = buffer.data(fk + 350);
    const auto *fk_351 = buffer.data(fk + 351);
    const auto *fk_352 = buffer.data(fk + 352);
    const auto *fk_353 = buffer.data(fk + 353);
    const auto *fk_354 = buffer.data(fk + 354);
    const auto *fk_355 = buffer.data(fk + 355);
    const auto *fk_356 = buffer.data(fk + 356);
    const auto *fk_357 = buffer.data(fk + 357);
    const auto *fk_358 = buffer.data(fk + 358);
    const auto *fk_359 = buffer.data(fk + 359);

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

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, fk_300, fk_301, fk_302, fk_303, \
                         fk_304, hk_300, hk_301, hk_302, hk_303, \
                         hk_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -fk_300[k]
                   + f_0 * hk_300[k];

        t_301[k] = -fk_301[k]
                   + f_0 * hk_301[k];

        t_302[k] = -fk_302[k]
                   + f_0 * hk_302[k];

        t_303[k] = -fk_303[k]
                   + f_0 * hk_303[k];

        t_304[k] = -fk_304[k]
                   + f_0 * hk_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, fk_305, fk_306, fk_307, fk_308, \
                         fk_309, hk_305, hk_306, hk_307, hk_308, \
                         hk_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -fk_305[k]
                   + f_0 * hk_305[k];

        t_306[k] = -fk_306[k]
                   + f_0 * hk_306[k];

        t_307[k] = -fk_307[k]
                   + f_0 * hk_307[k];

        t_308[k] = -fk_308[k]
                   + f_0 * hk_308[k];

        t_309[k] = -fk_309[k]
                   + f_0 * hk_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, fk_310, fk_311, fk_312, fk_313, \
                         fk_314, hk_310, hk_311, hk_312, hk_313, \
                         hk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -fk_310[k]
                   + f_0 * hk_310[k];

        t_311[k] = -fk_311[k]
                   + f_0 * hk_311[k];

        t_312[k] = -fk_312[k]
                   + f_0 * hk_312[k];

        t_313[k] = -fk_313[k]
                   + f_0 * hk_313[k];

        t_314[k] = -fk_314[k]
                   + f_0 * hk_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, fk_315, fk_316, fk_317, fk_318, \
                         fk_319, hk_315, hk_316, hk_317, hk_318, \
                         hk_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -fk_315[k]
                   + f_0 * hk_315[k];

        t_316[k] = -fk_316[k]
                   + f_0 * hk_316[k];

        t_317[k] = -fk_317[k]
                   + f_0 * hk_317[k];

        t_318[k] = -fk_318[k]
                   + f_0 * hk_318[k];

        t_319[k] = -fk_319[k]
                   + f_0 * hk_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, fk_320, fk_321, fk_322, fk_323, \
                         fk_324, hk_320, hk_321, hk_322, hk_323, \
                         hk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -fk_320[k]
                   + f_0 * hk_320[k];

        t_321[k] = -fk_321[k]
                   + f_0 * hk_321[k];

        t_322[k] = -fk_322[k]
                   + f_0 * hk_322[k];

        t_323[k] = -fk_323[k]
                   + f_0 * hk_323[k];

        t_324[k] = -fk_324[k]
                   + f_0 * hk_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, fk_325, fk_326, fk_327, fk_328, \
                         fk_329, hk_325, hk_326, hk_327, hk_328, \
                         hk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -fk_325[k]
                   + f_0 * hk_325[k];

        t_326[k] = -fk_326[k]
                   + f_0 * hk_326[k];

        t_327[k] = -fk_327[k]
                   + f_0 * hk_327[k];

        t_328[k] = -fk_328[k]
                   + f_0 * hk_328[k];

        t_329[k] = -fk_329[k]
                   + f_0 * hk_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, fk_330, fk_331, fk_332, fk_333, \
                         fk_334, hk_330, hk_331, hk_332, hk_333, \
                         hk_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -fk_330[k]
                   + f_0 * hk_330[k];

        t_331[k] = -fk_331[k]
                   + f_0 * hk_331[k];

        t_332[k] = -fk_332[k]
                   + f_0 * hk_332[k];

        t_333[k] = -fk_333[k]
                   + f_0 * hk_333[k];

        t_334[k] = -fk_334[k]
                   + f_0 * hk_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, fk_335, fk_336, fk_337, fk_338, \
                         fk_339, hk_335, hk_336, hk_337, hk_338, \
                         hk_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -fk_335[k]
                   + f_0 * hk_335[k];

        t_336[k] = -fk_336[k]
                   + f_0 * hk_336[k];

        t_337[k] = -fk_337[k]
                   + f_0 * hk_337[k];

        t_338[k] = -fk_338[k]
                   + f_0 * hk_338[k];

        t_339[k] = -fk_339[k]
                   + f_0 * hk_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, fk_340, fk_341, fk_342, fk_343, \
                         fk_344, hk_340, hk_341, hk_342, hk_343, \
                         hk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -fk_340[k]
                   + f_0 * hk_340[k];

        t_341[k] = -fk_341[k]
                   + f_0 * hk_341[k];

        t_342[k] = -fk_342[k]
                   + f_0 * hk_342[k];

        t_343[k] = -fk_343[k]
                   + f_0 * hk_343[k];

        t_344[k] = -fk_344[k]
                   + f_0 * hk_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, fk_345, fk_346, fk_347, fk_348, \
                         fk_349, hk_345, hk_346, hk_347, hk_348, \
                         hk_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -fk_345[k]
                   + f_0 * hk_345[k];

        t_346[k] = -fk_346[k]
                   + f_0 * hk_346[k];

        t_347[k] = -fk_347[k]
                   + f_0 * hk_347[k];

        t_348[k] = -fk_348[k]
                   + f_0 * hk_348[k];

        t_349[k] = -fk_349[k]
                   + f_0 * hk_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, fk_350, fk_351, fk_352, fk_353, \
                         fk_354, hk_350, hk_351, hk_352, hk_353, \
                         hk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -fk_350[k]
                   + f_0 * hk_350[k];

        t_351[k] = -fk_351[k]
                   + f_0 * hk_351[k];

        t_352[k] = -fk_352[k]
                   + f_0 * hk_352[k];

        t_353[k] = -fk_353[k]
                   + f_0 * hk_353[k];

        t_354[k] = -fk_354[k]
                   + f_0 * hk_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, fk_355, fk_356, fk_357, fk_358, \
                         fk_359, hk_355, hk_356, hk_357, hk_358, \
                         hk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -fk_355[k]
                   + f_0 * hk_355[k];

        t_356[k] = -fk_356[k]
                   + f_0 * hk_356[k];

        t_357[k] = -fk_357[k]
                   + f_0 * hk_357[k];

        t_358[k] = -fk_358[k]
                   + f_0 * hk_358[k];

        t_359[k] = -fk_359[k]
                   + f_0 * hk_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, t_365, t_366, t_367, hk_360, \
                         hk_361, hk_362, hk_363, hk_364, hk_365, hk_366, \
                         hk_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_0 * hk_360[k];

        t_361[k] = f_0 * hk_361[k];

        t_362[k] = f_0 * hk_362[k];

        t_363[k] = f_0 * hk_363[k];

        t_364[k] = f_0 * hk_364[k];

        t_365[k] = f_0 * hk_365[k];

        t_366[k] = f_0 * hk_366[k];

        t_367[k] = f_0 * hk_367[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, t_373, t_374, t_375, hk_368, \
                         hk_369, hk_370, hk_371, hk_372, hk_373, hk_374, \
                         hk_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_0 * hk_368[k];

        t_369[k] = f_0 * hk_369[k];

        t_370[k] = f_0 * hk_370[k];

        t_371[k] = f_0 * hk_371[k];

        t_372[k] = f_0 * hk_372[k];

        t_373[k] = f_0 * hk_373[k];

        t_374[k] = f_0 * hk_374[k];

        t_375[k] = f_0 * hk_375[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, t_380, t_381, t_382, t_383, hk_376, \
                         hk_377, hk_378, hk_379, hk_380, hk_381, hk_382, \
                         hk_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_0 * hk_376[k];

        t_377[k] = f_0 * hk_377[k];

        t_378[k] = f_0 * hk_378[k];

        t_379[k] = f_0 * hk_379[k];

        t_380[k] = f_0 * hk_380[k];

        t_381[k] = f_0 * hk_381[k];

        t_382[k] = f_0 * hk_382[k];

        t_383[k] = f_0 * hk_383[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, t_388, t_389, t_390, t_391, hk_384, \
                         hk_385, hk_386, hk_387, hk_388, hk_389, hk_390, \
                         hk_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_0 * hk_384[k];

        t_385[k] = f_0 * hk_385[k];

        t_386[k] = f_0 * hk_386[k];

        t_387[k] = f_0 * hk_387[k];

        t_388[k] = f_0 * hk_388[k];

        t_389[k] = f_0 * hk_389[k];

        t_390[k] = f_0 * hk_390[k];

        t_391[k] = f_0 * hk_391[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, t_397, t_398, t_399, hk_392, \
                         hk_393, hk_394, hk_395, hk_396, hk_397, hk_398, \
                         hk_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_0 * hk_392[k];

        t_393[k] = f_0 * hk_393[k];

        t_394[k] = f_0 * hk_394[k];

        t_395[k] = f_0 * hk_395[k];

        t_396[k] = f_0 * hk_396[k];

        t_397[k] = f_0 * hk_397[k];

        t_398[k] = f_0 * hk_398[k];

        t_399[k] = f_0 * hk_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, t_405, t_406, t_407, hk_400, \
                         hk_401, hk_402, hk_403, hk_404, hk_405, hk_406, \
                         hk_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_0 * hk_400[k];

        t_401[k] = f_0 * hk_401[k];

        t_402[k] = f_0 * hk_402[k];

        t_403[k] = f_0 * hk_403[k];

        t_404[k] = f_0 * hk_404[k];

        t_405[k] = f_0 * hk_405[k];

        t_406[k] = f_0 * hk_406[k];

        t_407[k] = f_0 * hk_407[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, t_413, t_414, t_415, hk_408, \
                         hk_409, hk_410, hk_411, hk_412, hk_413, hk_414, \
                         hk_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_0 * hk_408[k];

        t_409[k] = f_0 * hk_409[k];

        t_410[k] = f_0 * hk_410[k];

        t_411[k] = f_0 * hk_411[k];

        t_412[k] = f_0 * hk_412[k];

        t_413[k] = f_0 * hk_413[k];

        t_414[k] = f_0 * hk_414[k];

        t_415[k] = f_0 * hk_415[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, t_420, t_421, t_422, t_423, hk_416, \
                         hk_417, hk_418, hk_419, hk_420, hk_421, hk_422, \
                         hk_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_0 * hk_416[k];

        t_417[k] = f_0 * hk_417[k];

        t_418[k] = f_0 * hk_418[k];

        t_419[k] = f_0 * hk_419[k];

        t_420[k] = f_0 * hk_420[k];

        t_421[k] = f_0 * hk_421[k];

        t_422[k] = f_0 * hk_422[k];

        t_423[k] = f_0 * hk_423[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, t_428, t_429, t_430, t_431, hk_424, \
                         hk_425, hk_426, hk_427, hk_428, hk_429, hk_430, \
                         hk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_0 * hk_424[k];

        t_425[k] = f_0 * hk_425[k];

        t_426[k] = f_0 * hk_426[k];

        t_427[k] = f_0 * hk_427[k];

        t_428[k] = f_0 * hk_428[k];

        t_429[k] = f_0 * hk_429[k];

        t_430[k] = f_0 * hk_430[k];

        t_431[k] = f_0 * hk_431[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, t_436, t_437, t_438, t_439, hk_432, \
                         hk_433, hk_434, hk_435, hk_436, hk_437, hk_438, \
                         hk_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_0 * hk_432[k];

        t_433[k] = f_0 * hk_433[k];

        t_434[k] = f_0 * hk_434[k];

        t_435[k] = f_0 * hk_435[k];

        t_436[k] = f_0 * hk_436[k];

        t_437[k] = f_0 * hk_437[k];

        t_438[k] = f_0 * hk_438[k];

        t_439[k] = f_0 * hk_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, t_445, t_446, t_447, hk_440, \
                         hk_441, hk_442, hk_443, hk_444, hk_445, hk_446, \
                         hk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_0 * hk_440[k];

        t_441[k] = f_0 * hk_441[k];

        t_442[k] = f_0 * hk_442[k];

        t_443[k] = f_0 * hk_443[k];

        t_444[k] = f_0 * hk_444[k];

        t_445[k] = f_0 * hk_445[k];

        t_446[k] = f_0 * hk_446[k];

        t_447[k] = f_0 * hk_447[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, t_452, t_453, t_454, t_455, hk_448, \
                         hk_449, hk_450, hk_451, hk_452, hk_453, hk_454, \
                         hk_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = f_0 * hk_448[k];

        t_449[k] = f_0 * hk_449[k];

        t_450[k] = f_0 * hk_450[k];

        t_451[k] = f_0 * hk_451[k];

        t_452[k] = f_0 * hk_452[k];

        t_453[k] = f_0 * hk_453[k];

        t_454[k] = f_0 * hk_454[k];

        t_455[k] = f_0 * hk_455[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, t_460, t_461, t_462, t_463, hk_456, \
                         hk_457, hk_458, hk_459, hk_460, hk_461, hk_462, \
                         hk_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_0 * hk_456[k];

        t_457[k] = f_0 * hk_457[k];

        t_458[k] = f_0 * hk_458[k];

        t_459[k] = f_0 * hk_459[k];

        t_460[k] = f_0 * hk_460[k];

        t_461[k] = f_0 * hk_461[k];

        t_462[k] = f_0 * hk_462[k];

        t_463[k] = f_0 * hk_463[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, t_468, t_469, t_470, t_471, hk_464, \
                         hk_465, hk_466, hk_467, hk_468, hk_469, hk_470, \
                         hk_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_0 * hk_464[k];

        t_465[k] = f_0 * hk_465[k];

        t_466[k] = f_0 * hk_466[k];

        t_467[k] = f_0 * hk_467[k];

        t_468[k] = f_0 * hk_468[k];

        t_469[k] = f_0 * hk_469[k];

        t_470[k] = f_0 * hk_470[k];

        t_471[k] = f_0 * hk_471[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, t_477, t_478, t_479, hk_472, \
                         hk_473, hk_474, hk_475, hk_476, hk_477, hk_478, \
                         hk_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = f_0 * hk_472[k];

        t_473[k] = f_0 * hk_473[k];

        t_474[k] = f_0 * hk_474[k];

        t_475[k] = f_0 * hk_475[k];

        t_476[k] = f_0 * hk_476[k];

        t_477[k] = f_0 * hk_477[k];

        t_478[k] = f_0 * hk_478[k];

        t_479[k] = f_0 * hk_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, t_485, t_486, t_487, hk_480, \
                         hk_481, hk_482, hk_483, hk_484, hk_485, hk_486, \
                         hk_487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_0 * hk_480[k];

        t_481[k] = f_0 * hk_481[k];

        t_482[k] = f_0 * hk_482[k];

        t_483[k] = f_0 * hk_483[k];

        t_484[k] = f_0 * hk_484[k];

        t_485[k] = f_0 * hk_485[k];

        t_486[k] = f_0 * hk_486[k];

        t_487[k] = f_0 * hk_487[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, t_492, t_493, t_494, t_495, hk_488, \
                         hk_489, hk_490, hk_491, hk_492, hk_493, hk_494, \
                         hk_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_0 * hk_488[k];

        t_489[k] = f_0 * hk_489[k];

        t_490[k] = f_0 * hk_490[k];

        t_491[k] = f_0 * hk_491[k];

        t_492[k] = f_0 * hk_492[k];

        t_493[k] = f_0 * hk_493[k];

        t_494[k] = f_0 * hk_494[k];

        t_495[k] = f_0 * hk_495[k];
    }
}

static auto
compute_prim_geom_10_gk_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hk, const size_t ncols,
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

#pragma omp simd aligned(t_496, t_497, t_498, t_499, t_500, t_501, t_502, t_503, hk_496, \
                         hk_497, hk_498, hk_499, hk_500, hk_501, hk_502, \
                         hk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_0 * hk_496[k];

        t_497[k] = f_0 * hk_497[k];

        t_498[k] = f_0 * hk_498[k];

        t_499[k] = f_0 * hk_499[k];

        t_500[k] = f_0 * hk_500[k];

        t_501[k] = f_0 * hk_501[k];

        t_502[k] = f_0 * hk_502[k];

        t_503[k] = f_0 * hk_503[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, t_509, t_510, t_511, hk_504, \
                         hk_505, hk_506, hk_507, hk_508, hk_509, hk_510, \
                         hk_511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_0 * hk_504[k];

        t_505[k] = f_0 * hk_505[k];

        t_506[k] = f_0 * hk_506[k];

        t_507[k] = f_0 * hk_507[k];

        t_508[k] = f_0 * hk_508[k];

        t_509[k] = f_0 * hk_509[k];

        t_510[k] = f_0 * hk_510[k];

        t_511[k] = f_0 * hk_511[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, t_517, t_518, t_519, hk_512, \
                         hk_513, hk_514, hk_515, hk_516, hk_517, hk_518, \
                         hk_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_0 * hk_512[k];

        t_513[k] = f_0 * hk_513[k];

        t_514[k] = f_0 * hk_514[k];

        t_515[k] = f_0 * hk_515[k];

        t_516[k] = f_0 * hk_516[k];

        t_517[k] = f_0 * hk_517[k];

        t_518[k] = f_0 * hk_518[k];

        t_519[k] = f_0 * hk_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, t_525, t_526, t_527, hk_520, \
                         hk_521, hk_522, hk_523, hk_524, hk_525, hk_526, \
                         hk_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_0 * hk_520[k];

        t_521[k] = f_0 * hk_521[k];

        t_522[k] = f_0 * hk_522[k];

        t_523[k] = f_0 * hk_523[k];

        t_524[k] = f_0 * hk_524[k];

        t_525[k] = f_0 * hk_525[k];

        t_526[k] = f_0 * hk_526[k];

        t_527[k] = f_0 * hk_527[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, t_533, t_534, t_535, hk_528, \
                         hk_529, hk_530, hk_531, hk_532, hk_533, hk_534, \
                         hk_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_0 * hk_528[k];

        t_529[k] = f_0 * hk_529[k];

        t_530[k] = f_0 * hk_530[k];

        t_531[k] = f_0 * hk_531[k];

        t_532[k] = f_0 * hk_532[k];

        t_533[k] = f_0 * hk_533[k];

        t_534[k] = f_0 * hk_534[k];

        t_535[k] = f_0 * hk_535[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, hk_536, hk_537, hk_538, \
                         hk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_0 * hk_536[k];

        t_537[k] = f_0 * hk_537[k];

        t_538[k] = f_0 * hk_538[k];

        t_539[k] = f_0 * hk_539[k];
    }
}

auto
compute_prim_geom_10_gk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t fk, const size_t hk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gk_electron_repulsion_0_piece0(buffer, target, fk, hk, ncols, alpha);

    compute_prim_geom_10_gk_electron_repulsion_0_piece1(buffer, target, fk, hk, ncols, alpha);

    compute_prim_geom_10_gk_electron_repulsion_0_piece2(buffer, target, fk, hk, ncols, alpha);

    compute_prim_geom_10_gk_electron_repulsion_0_piece3(buffer, target, hk, ncols, alpha);
}

static auto
compute_prim_geom_10_gk_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t hk,
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

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);
    const auto *fk_36 = buffer.data(fk + 36);
    const auto *fk_37 = buffer.data(fk + 37);
    const auto *fk_38 = buffer.data(fk + 38);
    const auto *fk_39 = buffer.data(fk + 39);
    const auto *fk_40 = buffer.data(fk + 40);
    const auto *fk_41 = buffer.data(fk + 41);
    const auto *fk_42 = buffer.data(fk + 42);
    const auto *fk_43 = buffer.data(fk + 43);
    const auto *fk_44 = buffer.data(fk + 44);
    const auto *fk_45 = buffer.data(fk + 45);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_47 = buffer.data(fk + 47);
    const auto *fk_48 = buffer.data(fk + 48);
    const auto *fk_49 = buffer.data(fk + 49);
    const auto *fk_50 = buffer.data(fk + 50);
    const auto *fk_51 = buffer.data(fk + 51);
    const auto *fk_52 = buffer.data(fk + 52);
    const auto *fk_53 = buffer.data(fk + 53);
    const auto *fk_54 = buffer.data(fk + 54);
    const auto *fk_55 = buffer.data(fk + 55);
    const auto *fk_56 = buffer.data(fk + 56);
    const auto *fk_57 = buffer.data(fk + 57);
    const auto *fk_58 = buffer.data(fk + 58);
    const auto *fk_59 = buffer.data(fk + 59);
    const auto *fk_60 = buffer.data(fk + 60);
    const auto *fk_61 = buffer.data(fk + 61);
    const auto *fk_62 = buffer.data(fk + 62);
    const auto *fk_63 = buffer.data(fk + 63);
    const auto *fk_64 = buffer.data(fk + 64);
    const auto *fk_65 = buffer.data(fk + 65);
    const auto *fk_66 = buffer.data(fk + 66);
    const auto *fk_67 = buffer.data(fk + 67);
    const auto *fk_68 = buffer.data(fk + 68);
    const auto *fk_69 = buffer.data(fk + 69);
    const auto *fk_70 = buffer.data(fk + 70);
    const auto *fk_71 = buffer.data(fk + 71);
    const auto *fk_72 = buffer.data(fk + 72);
    const auto *fk_73 = buffer.data(fk + 73);
    const auto *fk_74 = buffer.data(fk + 74);
    const auto *fk_75 = buffer.data(fk + 75);
    const auto *fk_76 = buffer.data(fk + 76);
    const auto *fk_77 = buffer.data(fk + 77);
    const auto *fk_78 = buffer.data(fk + 78);
    const auto *fk_79 = buffer.data(fk + 79);
    const auto *fk_80 = buffer.data(fk + 80);
    const auto *fk_81 = buffer.data(fk + 81);
    const auto *fk_82 = buffer.data(fk + 82);
    const auto *fk_83 = buffer.data(fk + 83);
    const auto *fk_84 = buffer.data(fk + 84);
    const auto *fk_85 = buffer.data(fk + 85);
    const auto *fk_86 = buffer.data(fk + 86);
    const auto *fk_87 = buffer.data(fk + 87);
    const auto *fk_88 = buffer.data(fk + 88);
    const auto *fk_89 = buffer.data(fk + 89);
    const auto *fk_90 = buffer.data(fk + 90);
    const auto *fk_91 = buffer.data(fk + 91);
    const auto *fk_92 = buffer.data(fk + 92);
    const auto *fk_93 = buffer.data(fk + 93);
    const auto *fk_94 = buffer.data(fk + 94);
    const auto *fk_95 = buffer.data(fk + 95);
    const auto *fk_96 = buffer.data(fk + 96);
    const auto *fk_97 = buffer.data(fk + 97);
    const auto *fk_98 = buffer.data(fk + 98);
    const auto *fk_99 = buffer.data(fk + 99);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hk_36, hk_37, hk_38, hk_39, \
                         hk_40, hk_41, hk_42, hk_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hk_36[k];

        t_1[k] = f_0 * hk_37[k];

        t_2[k] = f_0 * hk_38[k];

        t_3[k] = f_0 * hk_39[k];

        t_4[k] = f_0 * hk_40[k];

        t_5[k] = f_0 * hk_41[k];

        t_6[k] = f_0 * hk_42[k];

        t_7[k] = f_0 * hk_43[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, hk_44, hk_45, hk_46, \
                         hk_47, hk_48, hk_49, hk_50, hk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hk_44[k];

        t_9[k] = f_0 * hk_45[k];

        t_10[k] = f_0 * hk_46[k];

        t_11[k] = f_0 * hk_47[k];

        t_12[k] = f_0 * hk_48[k];

        t_13[k] = f_0 * hk_49[k];

        t_14[k] = f_0 * hk_50[k];

        t_15[k] = f_0 * hk_51[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, hk_52, hk_53, hk_54, \
                         hk_55, hk_56, hk_57, hk_58, hk_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * hk_52[k];

        t_17[k] = f_0 * hk_53[k];

        t_18[k] = f_0 * hk_54[k];

        t_19[k] = f_0 * hk_55[k];

        t_20[k] = f_0 * hk_56[k];

        t_21[k] = f_0 * hk_57[k];

        t_22[k] = f_0 * hk_58[k];

        t_23[k] = f_0 * hk_59[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, hk_60, hk_61, hk_62, \
                         hk_63, hk_64, hk_65, hk_66, hk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * hk_60[k];

        t_25[k] = f_0 * hk_61[k];

        t_26[k] = f_0 * hk_62[k];

        t_27[k] = f_0 * hk_63[k];

        t_28[k] = f_0 * hk_64[k];

        t_29[k] = f_0 * hk_65[k];

        t_30[k] = f_0 * hk_66[k];

        t_31[k] = f_0 * hk_67[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, fk_0, fk_1, hk_68, hk_69, hk_70, \
                         hk_71, hk_108, hk_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * hk_68[k];

        t_33[k] = f_0 * hk_69[k];

        t_34[k] = f_0 * hk_70[k];

        t_35[k] = f_0 * hk_71[k];

        t_36[k] = -fk_0[k]
                  + f_0 * hk_108[k];

        t_37[k] = -fk_1[k]
                  + f_0 * hk_109[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, fk_2, fk_3, fk_4, fk_5, fk_6, hk_110, \
                         hk_111, hk_112, hk_113, hk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -fk_2[k]
                  + f_0 * hk_110[k];

        t_39[k] = -fk_3[k]
                  + f_0 * hk_111[k];

        t_40[k] = -fk_4[k]
                  + f_0 * hk_112[k];

        t_41[k] = -fk_5[k]
                  + f_0 * hk_113[k];

        t_42[k] = -fk_6[k]
                  + f_0 * hk_114[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, fk_7, fk_8, fk_9, fk_10, fk_11, hk_115, \
                         hk_116, hk_117, hk_118, hk_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -fk_7[k]
                  + f_0 * hk_115[k];

        t_44[k] = -fk_8[k]
                  + f_0 * hk_116[k];

        t_45[k] = -fk_9[k]
                  + f_0 * hk_117[k];

        t_46[k] = -fk_10[k]
                  + f_0 * hk_118[k];

        t_47[k] = -fk_11[k]
                  + f_0 * hk_119[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, fk_12, fk_13, fk_14, fk_15, fk_16, \
                         hk_120, hk_121, hk_122, hk_123, hk_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -fk_12[k]
                  + f_0 * hk_120[k];

        t_49[k] = -fk_13[k]
                  + f_0 * hk_121[k];

        t_50[k] = -fk_14[k]
                  + f_0 * hk_122[k];

        t_51[k] = -fk_15[k]
                  + f_0 * hk_123[k];

        t_52[k] = -fk_16[k]
                  + f_0 * hk_124[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, fk_17, fk_18, fk_19, fk_20, fk_21, \
                         hk_125, hk_126, hk_127, hk_128, hk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -fk_17[k]
                  + f_0 * hk_125[k];

        t_54[k] = -fk_18[k]
                  + f_0 * hk_126[k];

        t_55[k] = -fk_19[k]
                  + f_0 * hk_127[k];

        t_56[k] = -fk_20[k]
                  + f_0 * hk_128[k];

        t_57[k] = -fk_21[k]
                  + f_0 * hk_129[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, fk_22, fk_23, fk_24, fk_25, fk_26, \
                         hk_130, hk_131, hk_132, hk_133, hk_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -fk_22[k]
                  + f_0 * hk_130[k];

        t_59[k] = -fk_23[k]
                  + f_0 * hk_131[k];

        t_60[k] = -fk_24[k]
                  + f_0 * hk_132[k];

        t_61[k] = -fk_25[k]
                  + f_0 * hk_133[k];

        t_62[k] = -fk_26[k]
                  + f_0 * hk_134[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, fk_27, fk_28, fk_29, fk_30, fk_31, \
                         hk_135, hk_136, hk_137, hk_138, hk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -fk_27[k]
                  + f_0 * hk_135[k];

        t_64[k] = -fk_28[k]
                  + f_0 * hk_136[k];

        t_65[k] = -fk_29[k]
                  + f_0 * hk_137[k];

        t_66[k] = -fk_30[k]
                  + f_0 * hk_138[k];

        t_67[k] = -fk_31[k]
                  + f_0 * hk_139[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, fk_32, fk_33, fk_34, fk_35, \
                         hk_140, hk_141, hk_142, hk_143, hk_144, \
                         hk_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -fk_32[k]
                  + f_0 * hk_140[k];

        t_69[k] = -fk_33[k]
                  + f_0 * hk_141[k];

        t_70[k] = -fk_34[k]
                  + f_0 * hk_142[k];

        t_71[k] = -fk_35[k]
                  + f_0 * hk_143[k];

        t_72[k] = f_0 * hk_144[k];

        t_73[k] = f_0 * hk_145[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, t_80, t_81, hk_146, hk_147, \
                         hk_148, hk_149, hk_150, hk_151, hk_152, \
                         hk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * hk_146[k];

        t_75[k] = f_0 * hk_147[k];

        t_76[k] = f_0 * hk_148[k];

        t_77[k] = f_0 * hk_149[k];

        t_78[k] = f_0 * hk_150[k];

        t_79[k] = f_0 * hk_151[k];

        t_80[k] = f_0 * hk_152[k];

        t_81[k] = f_0 * hk_153[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, t_87, t_88, t_89, hk_154, hk_155, \
                         hk_156, hk_157, hk_158, hk_159, hk_160, \
                         hk_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * hk_154[k];

        t_83[k] = f_0 * hk_155[k];

        t_84[k] = f_0 * hk_156[k];

        t_85[k] = f_0 * hk_157[k];

        t_86[k] = f_0 * hk_158[k];

        t_87[k] = f_0 * hk_159[k];

        t_88[k] = f_0 * hk_160[k];

        t_89[k] = f_0 * hk_161[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, hk_162, hk_163, \
                         hk_164, hk_165, hk_166, hk_167, hk_168, \
                         hk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * hk_162[k];

        t_91[k] = f_0 * hk_163[k];

        t_92[k] = f_0 * hk_164[k];

        t_93[k] = f_0 * hk_165[k];

        t_94[k] = f_0 * hk_166[k];

        t_95[k] = f_0 * hk_167[k];

        t_96[k] = f_0 * hk_168[k];

        t_97[k] = f_0 * hk_169[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, t_104, t_105, hk_170, hk_171, \
                         hk_172, hk_173, hk_174, hk_175, hk_176, \
                         hk_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * hk_170[k];

        t_99[k] = f_0 * hk_171[k];

        t_100[k] = f_0 * hk_172[k];

        t_101[k] = f_0 * hk_173[k];

        t_102[k] = f_0 * hk_174[k];

        t_103[k] = f_0 * hk_175[k];

        t_104[k] = f_0 * hk_176[k];

        t_105[k] = f_0 * hk_177[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, fk_36, fk_37, fk_38, fk_39, \
                         hk_178, hk_179, hk_216, hk_217, hk_218, \
                         hk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * hk_178[k];

        t_107[k] = f_0 * hk_179[k];

        t_108[k] = -2.0 * fk_36[k]
                   + f_0 * hk_216[k];

        t_109[k] = -2.0 * fk_37[k]
                   + f_0 * hk_217[k];

        t_110[k] = -2.0 * fk_38[k]
                   + f_0 * hk_218[k];

        t_111[k] = -2.0 * fk_39[k]
                   + f_0 * hk_219[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, fk_40, fk_41, fk_42, fk_43, fk_44, \
                         hk_220, hk_221, hk_222, hk_223, hk_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -2.0 * fk_40[k]
                   + f_0 * hk_220[k];

        t_113[k] = -2.0 * fk_41[k]
                   + f_0 * hk_221[k];

        t_114[k] = -2.0 * fk_42[k]
                   + f_0 * hk_222[k];

        t_115[k] = -2.0 * fk_43[k]
                   + f_0 * hk_223[k];

        t_116[k] = -2.0 * fk_44[k]
                   + f_0 * hk_224[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, fk_45, fk_46, fk_47, fk_48, fk_49, \
                         hk_225, hk_226, hk_227, hk_228, hk_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -2.0 * fk_45[k]
                   + f_0 * hk_225[k];

        t_118[k] = -2.0 * fk_46[k]
                   + f_0 * hk_226[k];

        t_119[k] = -2.0 * fk_47[k]
                   + f_0 * hk_227[k];

        t_120[k] = -2.0 * fk_48[k]
                   + f_0 * hk_228[k];

        t_121[k] = -2.0 * fk_49[k]
                   + f_0 * hk_229[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, fk_50, fk_51, fk_52, fk_53, fk_54, \
                         hk_230, hk_231, hk_232, hk_233, hk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * fk_50[k]
                   + f_0 * hk_230[k];

        t_123[k] = -2.0 * fk_51[k]
                   + f_0 * hk_231[k];

        t_124[k] = -2.0 * fk_52[k]
                   + f_0 * hk_232[k];

        t_125[k] = -2.0 * fk_53[k]
                   + f_0 * hk_233[k];

        t_126[k] = -2.0 * fk_54[k]
                   + f_0 * hk_234[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, fk_55, fk_56, fk_57, fk_58, fk_59, \
                         hk_235, hk_236, hk_237, hk_238, hk_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * fk_55[k]
                   + f_0 * hk_235[k];

        t_128[k] = -2.0 * fk_56[k]
                   + f_0 * hk_236[k];

        t_129[k] = -2.0 * fk_57[k]
                   + f_0 * hk_237[k];

        t_130[k] = -2.0 * fk_58[k]
                   + f_0 * hk_238[k];

        t_131[k] = -2.0 * fk_59[k]
                   + f_0 * hk_239[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, fk_60, fk_61, fk_62, fk_63, fk_64, \
                         hk_240, hk_241, hk_242, hk_243, hk_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -2.0 * fk_60[k]
                   + f_0 * hk_240[k];

        t_133[k] = -2.0 * fk_61[k]
                   + f_0 * hk_241[k];

        t_134[k] = -2.0 * fk_62[k]
                   + f_0 * hk_242[k];

        t_135[k] = -2.0 * fk_63[k]
                   + f_0 * hk_243[k];

        t_136[k] = -2.0 * fk_64[k]
                   + f_0 * hk_244[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, fk_65, fk_66, fk_67, fk_68, fk_69, \
                         hk_245, hk_246, hk_247, hk_248, hk_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -2.0 * fk_65[k]
                   + f_0 * hk_245[k];

        t_138[k] = -2.0 * fk_66[k]
                   + f_0 * hk_246[k];

        t_139[k] = -2.0 * fk_67[k]
                   + f_0 * hk_247[k];

        t_140[k] = -2.0 * fk_68[k]
                   + f_0 * hk_248[k];

        t_141[k] = -2.0 * fk_69[k]
                   + f_0 * hk_249[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, fk_70, fk_71, fk_72, fk_73, fk_74, \
                         hk_250, hk_251, hk_252, hk_253, hk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -2.0 * fk_70[k]
                   + f_0 * hk_250[k];

        t_143[k] = -2.0 * fk_71[k]
                   + f_0 * hk_251[k];

        t_144[k] = -fk_72[k]
                   + f_0 * hk_252[k];

        t_145[k] = -fk_73[k]
                   + f_0 * hk_253[k];

        t_146[k] = -fk_74[k]
                   + f_0 * hk_254[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, fk_75, fk_76, fk_77, fk_78, fk_79, \
                         hk_255, hk_256, hk_257, hk_258, hk_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -fk_75[k]
                   + f_0 * hk_255[k];

        t_148[k] = -fk_76[k]
                   + f_0 * hk_256[k];

        t_149[k] = -fk_77[k]
                   + f_0 * hk_257[k];

        t_150[k] = -fk_78[k]
                   + f_0 * hk_258[k];

        t_151[k] = -fk_79[k]
                   + f_0 * hk_259[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, fk_80, fk_81, fk_82, fk_83, fk_84, \
                         hk_260, hk_261, hk_262, hk_263, hk_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -fk_80[k]
                   + f_0 * hk_260[k];

        t_153[k] = -fk_81[k]
                   + f_0 * hk_261[k];

        t_154[k] = -fk_82[k]
                   + f_0 * hk_262[k];

        t_155[k] = -fk_83[k]
                   + f_0 * hk_263[k];

        t_156[k] = -fk_84[k]
                   + f_0 * hk_264[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, fk_85, fk_86, fk_87, fk_88, fk_89, \
                         hk_265, hk_266, hk_267, hk_268, hk_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -fk_85[k]
                   + f_0 * hk_265[k];

        t_158[k] = -fk_86[k]
                   + f_0 * hk_266[k];

        t_159[k] = -fk_87[k]
                   + f_0 * hk_267[k];

        t_160[k] = -fk_88[k]
                   + f_0 * hk_268[k];

        t_161[k] = -fk_89[k]
                   + f_0 * hk_269[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, fk_90, fk_91, fk_92, fk_93, fk_94, \
                         hk_270, hk_271, hk_272, hk_273, hk_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -fk_90[k]
                   + f_0 * hk_270[k];

        t_163[k] = -fk_91[k]
                   + f_0 * hk_271[k];

        t_164[k] = -fk_92[k]
                   + f_0 * hk_272[k];

        t_165[k] = -fk_93[k]
                   + f_0 * hk_273[k];

        t_166[k] = -fk_94[k]
                   + f_0 * hk_274[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, fk_95, fk_96, fk_97, fk_98, fk_99, \
                         hk_275, hk_276, hk_277, hk_278, hk_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -fk_95[k]
                   + f_0 * hk_275[k];

        t_168[k] = -fk_96[k]
                   + f_0 * hk_276[k];

        t_169[k] = -fk_97[k]
                   + f_0 * hk_277[k];

        t_170[k] = -fk_98[k]
                   + f_0 * hk_278[k];

        t_171[k] = -fk_99[k]
                   + f_0 * hk_279[k];
    }
}

static auto
compute_prim_geom_10_gk_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t hk,
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

    const auto *fk_100 = buffer.data(fk + 100);
    const auto *fk_101 = buffer.data(fk + 101);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_106 = buffer.data(fk + 106);
    const auto *fk_107 = buffer.data(fk + 107);
    const auto *fk_108 = buffer.data(fk + 108);
    const auto *fk_109 = buffer.data(fk + 109);
    const auto *fk_110 = buffer.data(fk + 110);
    const auto *fk_111 = buffer.data(fk + 111);
    const auto *fk_112 = buffer.data(fk + 112);
    const auto *fk_113 = buffer.data(fk + 113);
    const auto *fk_114 = buffer.data(fk + 114);
    const auto *fk_115 = buffer.data(fk + 115);
    const auto *fk_116 = buffer.data(fk + 116);
    const auto *fk_117 = buffer.data(fk + 117);
    const auto *fk_118 = buffer.data(fk + 118);
    const auto *fk_119 = buffer.data(fk + 119);
    const auto *fk_120 = buffer.data(fk + 120);
    const auto *fk_121 = buffer.data(fk + 121);
    const auto *fk_122 = buffer.data(fk + 122);
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_124 = buffer.data(fk + 124);
    const auto *fk_125 = buffer.data(fk + 125);
    const auto *fk_126 = buffer.data(fk + 126);
    const auto *fk_127 = buffer.data(fk + 127);
    const auto *fk_128 = buffer.data(fk + 128);
    const auto *fk_129 = buffer.data(fk + 129);
    const auto *fk_130 = buffer.data(fk + 130);
    const auto *fk_131 = buffer.data(fk + 131);
    const auto *fk_132 = buffer.data(fk + 132);
    const auto *fk_133 = buffer.data(fk + 133);
    const auto *fk_134 = buffer.data(fk + 134);
    const auto *fk_135 = buffer.data(fk + 135);
    const auto *fk_136 = buffer.data(fk + 136);
    const auto *fk_137 = buffer.data(fk + 137);
    const auto *fk_138 = buffer.data(fk + 138);
    const auto *fk_139 = buffer.data(fk + 139);
    const auto *fk_140 = buffer.data(fk + 140);
    const auto *fk_141 = buffer.data(fk + 141);
    const auto *fk_142 = buffer.data(fk + 142);
    const auto *fk_143 = buffer.data(fk + 143);
    const auto *fk_144 = buffer.data(fk + 144);
    const auto *fk_145 = buffer.data(fk + 145);
    const auto *fk_146 = buffer.data(fk + 146);
    const auto *fk_147 = buffer.data(fk + 147);
    const auto *fk_148 = buffer.data(fk + 148);
    const auto *fk_149 = buffer.data(fk + 149);
    const auto *fk_150 = buffer.data(fk + 150);
    const auto *fk_151 = buffer.data(fk + 151);
    const auto *fk_152 = buffer.data(fk + 152);
    const auto *fk_153 = buffer.data(fk + 153);
    const auto *fk_154 = buffer.data(fk + 154);
    const auto *fk_155 = buffer.data(fk + 155);
    const auto *fk_156 = buffer.data(fk + 156);
    const auto *fk_157 = buffer.data(fk + 157);
    const auto *fk_158 = buffer.data(fk + 158);
    const auto *fk_159 = buffer.data(fk + 159);
    const auto *fk_160 = buffer.data(fk + 160);
    const auto *fk_161 = buffer.data(fk + 161);
    const auto *fk_162 = buffer.data(fk + 162);
    const auto *fk_163 = buffer.data(fk + 163);
    const auto *fk_164 = buffer.data(fk + 164);
    const auto *fk_165 = buffer.data(fk + 165);
    const auto *fk_166 = buffer.data(fk + 166);
    const auto *fk_167 = buffer.data(fk + 167);
    const auto *fk_168 = buffer.data(fk + 168);
    const auto *fk_169 = buffer.data(fk + 169);
    const auto *fk_170 = buffer.data(fk + 170);
    const auto *fk_171 = buffer.data(fk + 171);
    const auto *fk_172 = buffer.data(fk + 172);
    const auto *fk_173 = buffer.data(fk + 173);
    const auto *fk_174 = buffer.data(fk + 174);
    const auto *fk_175 = buffer.data(fk + 175);
    const auto *fk_176 = buffer.data(fk + 176);
    const auto *fk_177 = buffer.data(fk + 177);
    const auto *fk_178 = buffer.data(fk + 178);
    const auto *fk_179 = buffer.data(fk + 179);
    const auto *fk_180 = buffer.data(fk + 180);
    const auto *fk_181 = buffer.data(fk + 181);
    const auto *fk_182 = buffer.data(fk + 182);
    const auto *fk_183 = buffer.data(fk + 183);
    const auto *fk_184 = buffer.data(fk + 184);
    const auto *fk_185 = buffer.data(fk + 185);
    const auto *fk_186 = buffer.data(fk + 186);
    const auto *fk_187 = buffer.data(fk + 187);
    const auto *fk_188 = buffer.data(fk + 188);
    const auto *fk_189 = buffer.data(fk + 189);
    const auto *fk_190 = buffer.data(fk + 190);
    const auto *fk_191 = buffer.data(fk + 191);
    const auto *fk_192 = buffer.data(fk + 192);
    const auto *fk_193 = buffer.data(fk + 193);
    const auto *fk_194 = buffer.data(fk + 194);
    const auto *fk_195 = buffer.data(fk + 195);
    const auto *fk_196 = buffer.data(fk + 196);
    const auto *fk_197 = buffer.data(fk + 197);
    const auto *fk_198 = buffer.data(fk + 198);
    const auto *fk_199 = buffer.data(fk + 199);
    const auto *fk_200 = buffer.data(fk + 200);
    const auto *fk_201 = buffer.data(fk + 201);
    const auto *fk_202 = buffer.data(fk + 202);
    const auto *fk_203 = buffer.data(fk + 203);
    const auto *fk_204 = buffer.data(fk + 204);
    const auto *fk_205 = buffer.data(fk + 205);
    const auto *fk_206 = buffer.data(fk + 206);
    const auto *fk_207 = buffer.data(fk + 207);
    const auto *fk_208 = buffer.data(fk + 208);
    const auto *fk_209 = buffer.data(fk + 209);
    const auto *fk_210 = buffer.data(fk + 210);
    const auto *fk_211 = buffer.data(fk + 211);
    const auto *fk_212 = buffer.data(fk + 212);
    const auto *fk_213 = buffer.data(fk + 213);
    const auto *fk_214 = buffer.data(fk + 214);
    const auto *fk_215 = buffer.data(fk + 215);

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
    const auto *hk_477 = buffer.data(hk + 477);

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, fk_100, fk_101, fk_102, fk_103, \
                         fk_104, hk_280, hk_281, hk_282, hk_283, \
                         hk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -fk_100[k]
                   + f_0 * hk_280[k];

        t_173[k] = -fk_101[k]
                   + f_0 * hk_281[k];

        t_174[k] = -fk_102[k]
                   + f_0 * hk_282[k];

        t_175[k] = -fk_103[k]
                   + f_0 * hk_283[k];

        t_176[k] = -fk_104[k]
                   + f_0 * hk_284[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, t_182, fk_105, fk_106, fk_107, \
                         hk_285, hk_286, hk_287, hk_288, hk_289, \
                         hk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -fk_105[k]
                   + f_0 * hk_285[k];

        t_178[k] = -fk_106[k]
                   + f_0 * hk_286[k];

        t_179[k] = -fk_107[k]
                   + f_0 * hk_287[k];

        t_180[k] = f_0 * hk_288[k];

        t_181[k] = f_0 * hk_289[k];

        t_182[k] = f_0 * hk_290[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, t_188, t_189, t_190, hk_291, \
                         hk_292, hk_293, hk_294, hk_295, hk_296, hk_297, \
                         hk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_0 * hk_291[k];

        t_184[k] = f_0 * hk_292[k];

        t_185[k] = f_0 * hk_293[k];

        t_186[k] = f_0 * hk_294[k];

        t_187[k] = f_0 * hk_295[k];

        t_188[k] = f_0 * hk_296[k];

        t_189[k] = f_0 * hk_297[k];

        t_190[k] = f_0 * hk_298[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, t_196, t_197, t_198, hk_299, \
                         hk_300, hk_301, hk_302, hk_303, hk_304, hk_305, \
                         hk_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_0 * hk_299[k];

        t_192[k] = f_0 * hk_300[k];

        t_193[k] = f_0 * hk_301[k];

        t_194[k] = f_0 * hk_302[k];

        t_195[k] = f_0 * hk_303[k];

        t_196[k] = f_0 * hk_304[k];

        t_197[k] = f_0 * hk_305[k];

        t_198[k] = f_0 * hk_306[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, t_204, t_205, t_206, hk_307, \
                         hk_308, hk_309, hk_310, hk_311, hk_312, hk_313, \
                         hk_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_0 * hk_307[k];

        t_200[k] = f_0 * hk_308[k];

        t_201[k] = f_0 * hk_309[k];

        t_202[k] = f_0 * hk_310[k];

        t_203[k] = f_0 * hk_311[k];

        t_204[k] = f_0 * hk_312[k];

        t_205[k] = f_0 * hk_313[k];

        t_206[k] = f_0 * hk_314[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, t_213, t_214, hk_315, \
                         hk_316, hk_317, hk_318, hk_319, hk_320, hk_321, \
                         hk_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_0 * hk_315[k];

        t_208[k] = f_0 * hk_316[k];

        t_209[k] = f_0 * hk_317[k];

        t_210[k] = f_0 * hk_318[k];

        t_211[k] = f_0 * hk_319[k];

        t_212[k] = f_0 * hk_320[k];

        t_213[k] = f_0 * hk_321[k];

        t_214[k] = f_0 * hk_322[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, fk_108, fk_109, fk_110, fk_111, \
                         hk_323, hk_360, hk_361, hk_362, hk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_0 * hk_323[k];

        t_216[k] = -3.0 * fk_108[k]
                   + f_0 * hk_360[k];

        t_217[k] = -3.0 * fk_109[k]
                   + f_0 * hk_361[k];

        t_218[k] = -3.0 * fk_110[k]
                   + f_0 * hk_362[k];

        t_219[k] = -3.0 * fk_111[k]
                   + f_0 * hk_363[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, fk_112, fk_113, fk_114, fk_115, \
                         fk_116, hk_364, hk_365, hk_366, hk_367, \
                         hk_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -3.0 * fk_112[k]
                   + f_0 * hk_364[k];

        t_221[k] = -3.0 * fk_113[k]
                   + f_0 * hk_365[k];

        t_222[k] = -3.0 * fk_114[k]
                   + f_0 * hk_366[k];

        t_223[k] = -3.0 * fk_115[k]
                   + f_0 * hk_367[k];

        t_224[k] = -3.0 * fk_116[k]
                   + f_0 * hk_368[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, fk_117, fk_118, fk_119, fk_120, \
                         fk_121, hk_369, hk_370, hk_371, hk_372, \
                         hk_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -3.0 * fk_117[k]
                   + f_0 * hk_369[k];

        t_226[k] = -3.0 * fk_118[k]
                   + f_0 * hk_370[k];

        t_227[k] = -3.0 * fk_119[k]
                   + f_0 * hk_371[k];

        t_228[k] = -3.0 * fk_120[k]
                   + f_0 * hk_372[k];

        t_229[k] = -3.0 * fk_121[k]
                   + f_0 * hk_373[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, fk_122, fk_123, fk_124, fk_125, \
                         fk_126, hk_374, hk_375, hk_376, hk_377, \
                         hk_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -3.0 * fk_122[k]
                   + f_0 * hk_374[k];

        t_231[k] = -3.0 * fk_123[k]
                   + f_0 * hk_375[k];

        t_232[k] = -3.0 * fk_124[k]
                   + f_0 * hk_376[k];

        t_233[k] = -3.0 * fk_125[k]
                   + f_0 * hk_377[k];

        t_234[k] = -3.0 * fk_126[k]
                   + f_0 * hk_378[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, fk_127, fk_128, fk_129, fk_130, \
                         fk_131, hk_379, hk_380, hk_381, hk_382, \
                         hk_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -3.0 * fk_127[k]
                   + f_0 * hk_379[k];

        t_236[k] = -3.0 * fk_128[k]
                   + f_0 * hk_380[k];

        t_237[k] = -3.0 * fk_129[k]
                   + f_0 * hk_381[k];

        t_238[k] = -3.0 * fk_130[k]
                   + f_0 * hk_382[k];

        t_239[k] = -3.0 * fk_131[k]
                   + f_0 * hk_383[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, fk_132, fk_133, fk_134, fk_135, \
                         fk_136, hk_384, hk_385, hk_386, hk_387, \
                         hk_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * fk_132[k]
                   + f_0 * hk_384[k];

        t_241[k] = -3.0 * fk_133[k]
                   + f_0 * hk_385[k];

        t_242[k] = -3.0 * fk_134[k]
                   + f_0 * hk_386[k];

        t_243[k] = -3.0 * fk_135[k]
                   + f_0 * hk_387[k];

        t_244[k] = -3.0 * fk_136[k]
                   + f_0 * hk_388[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, fk_137, fk_138, fk_139, fk_140, \
                         fk_141, hk_389, hk_390, hk_391, hk_392, \
                         hk_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * fk_137[k]
                   + f_0 * hk_389[k];

        t_246[k] = -3.0 * fk_138[k]
                   + f_0 * hk_390[k];

        t_247[k] = -3.0 * fk_139[k]
                   + f_0 * hk_391[k];

        t_248[k] = -3.0 * fk_140[k]
                   + f_0 * hk_392[k];

        t_249[k] = -3.0 * fk_141[k]
                   + f_0 * hk_393[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, fk_142, fk_143, fk_144, fk_145, \
                         fk_146, hk_394, hk_395, hk_396, hk_397, \
                         hk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -3.0 * fk_142[k]
                   + f_0 * hk_394[k];

        t_251[k] = -3.0 * fk_143[k]
                   + f_0 * hk_395[k];

        t_252[k] = -2.0 * fk_144[k]
                   + f_0 * hk_396[k];

        t_253[k] = -2.0 * fk_145[k]
                   + f_0 * hk_397[k];

        t_254[k] = -2.0 * fk_146[k]
                   + f_0 * hk_398[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, fk_147, fk_148, fk_149, fk_150, \
                         fk_151, hk_399, hk_400, hk_401, hk_402, \
                         hk_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * fk_147[k]
                   + f_0 * hk_399[k];

        t_256[k] = -2.0 * fk_148[k]
                   + f_0 * hk_400[k];

        t_257[k] = -2.0 * fk_149[k]
                   + f_0 * hk_401[k];

        t_258[k] = -2.0 * fk_150[k]
                   + f_0 * hk_402[k];

        t_259[k] = -2.0 * fk_151[k]
                   + f_0 * hk_403[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, fk_152, fk_153, fk_154, fk_155, \
                         fk_156, hk_404, hk_405, hk_406, hk_407, \
                         hk_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * fk_152[k]
                   + f_0 * hk_404[k];

        t_261[k] = -2.0 * fk_153[k]
                   + f_0 * hk_405[k];

        t_262[k] = -2.0 * fk_154[k]
                   + f_0 * hk_406[k];

        t_263[k] = -2.0 * fk_155[k]
                   + f_0 * hk_407[k];

        t_264[k] = -2.0 * fk_156[k]
                   + f_0 * hk_408[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, fk_157, fk_158, fk_159, fk_160, \
                         fk_161, hk_409, hk_410, hk_411, hk_412, \
                         hk_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * fk_157[k]
                   + f_0 * hk_409[k];

        t_266[k] = -2.0 * fk_158[k]
                   + f_0 * hk_410[k];

        t_267[k] = -2.0 * fk_159[k]
                   + f_0 * hk_411[k];

        t_268[k] = -2.0 * fk_160[k]
                   + f_0 * hk_412[k];

        t_269[k] = -2.0 * fk_161[k]
                   + f_0 * hk_413[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, fk_162, fk_163, fk_164, fk_165, \
                         fk_166, hk_414, hk_415, hk_416, hk_417, \
                         hk_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * fk_162[k]
                   + f_0 * hk_414[k];

        t_271[k] = -2.0 * fk_163[k]
                   + f_0 * hk_415[k];

        t_272[k] = -2.0 * fk_164[k]
                   + f_0 * hk_416[k];

        t_273[k] = -2.0 * fk_165[k]
                   + f_0 * hk_417[k];

        t_274[k] = -2.0 * fk_166[k]
                   + f_0 * hk_418[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, fk_167, fk_168, fk_169, fk_170, \
                         fk_171, hk_419, hk_420, hk_421, hk_422, \
                         hk_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * fk_167[k]
                   + f_0 * hk_419[k];

        t_276[k] = -2.0 * fk_168[k]
                   + f_0 * hk_420[k];

        t_277[k] = -2.0 * fk_169[k]
                   + f_0 * hk_421[k];

        t_278[k] = -2.0 * fk_170[k]
                   + f_0 * hk_422[k];

        t_279[k] = -2.0 * fk_171[k]
                   + f_0 * hk_423[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, fk_172, fk_173, fk_174, fk_175, \
                         fk_176, hk_424, hk_425, hk_426, hk_427, \
                         hk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * fk_172[k]
                   + f_0 * hk_424[k];

        t_281[k] = -2.0 * fk_173[k]
                   + f_0 * hk_425[k];

        t_282[k] = -2.0 * fk_174[k]
                   + f_0 * hk_426[k];

        t_283[k] = -2.0 * fk_175[k]
                   + f_0 * hk_427[k];

        t_284[k] = -2.0 * fk_176[k]
                   + f_0 * hk_428[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, fk_177, fk_178, fk_179, fk_180, \
                         fk_181, hk_429, hk_430, hk_431, hk_432, \
                         hk_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -2.0 * fk_177[k]
                   + f_0 * hk_429[k];

        t_286[k] = -2.0 * fk_178[k]
                   + f_0 * hk_430[k];

        t_287[k] = -2.0 * fk_179[k]
                   + f_0 * hk_431[k];

        t_288[k] = -fk_180[k]
                   + f_0 * hk_432[k];

        t_289[k] = -fk_181[k]
                   + f_0 * hk_433[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, fk_182, fk_183, fk_184, fk_185, \
                         fk_186, hk_434, hk_435, hk_436, hk_437, \
                         hk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -fk_182[k]
                   + f_0 * hk_434[k];

        t_291[k] = -fk_183[k]
                   + f_0 * hk_435[k];

        t_292[k] = -fk_184[k]
                   + f_0 * hk_436[k];

        t_293[k] = -fk_185[k]
                   + f_0 * hk_437[k];

        t_294[k] = -fk_186[k]
                   + f_0 * hk_438[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, fk_187, fk_188, fk_189, fk_190, \
                         fk_191, hk_439, hk_440, hk_441, hk_442, \
                         hk_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -fk_187[k]
                   + f_0 * hk_439[k];

        t_296[k] = -fk_188[k]
                   + f_0 * hk_440[k];

        t_297[k] = -fk_189[k]
                   + f_0 * hk_441[k];

        t_298[k] = -fk_190[k]
                   + f_0 * hk_442[k];

        t_299[k] = -fk_191[k]
                   + f_0 * hk_443[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, fk_192, fk_193, fk_194, fk_195, \
                         fk_196, hk_444, hk_445, hk_446, hk_447, \
                         hk_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -fk_192[k]
                   + f_0 * hk_444[k];

        t_301[k] = -fk_193[k]
                   + f_0 * hk_445[k];

        t_302[k] = -fk_194[k]
                   + f_0 * hk_446[k];

        t_303[k] = -fk_195[k]
                   + f_0 * hk_447[k];

        t_304[k] = -fk_196[k]
                   + f_0 * hk_448[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, fk_197, fk_198, fk_199, fk_200, \
                         fk_201, hk_449, hk_450, hk_451, hk_452, \
                         hk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -fk_197[k]
                   + f_0 * hk_449[k];

        t_306[k] = -fk_198[k]
                   + f_0 * hk_450[k];

        t_307[k] = -fk_199[k]
                   + f_0 * hk_451[k];

        t_308[k] = -fk_200[k]
                   + f_0 * hk_452[k];

        t_309[k] = -fk_201[k]
                   + f_0 * hk_453[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, fk_202, fk_203, fk_204, fk_205, \
                         fk_206, hk_454, hk_455, hk_456, hk_457, \
                         hk_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -fk_202[k]
                   + f_0 * hk_454[k];

        t_311[k] = -fk_203[k]
                   + f_0 * hk_455[k];

        t_312[k] = -fk_204[k]
                   + f_0 * hk_456[k];

        t_313[k] = -fk_205[k]
                   + f_0 * hk_457[k];

        t_314[k] = -fk_206[k]
                   + f_0 * hk_458[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, fk_207, fk_208, fk_209, fk_210, \
                         fk_211, hk_459, hk_460, hk_461, hk_462, \
                         hk_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -fk_207[k]
                   + f_0 * hk_459[k];

        t_316[k] = -fk_208[k]
                   + f_0 * hk_460[k];

        t_317[k] = -fk_209[k]
                   + f_0 * hk_461[k];

        t_318[k] = -fk_210[k]
                   + f_0 * hk_462[k];

        t_319[k] = -fk_211[k]
                   + f_0 * hk_463[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, t_325, fk_212, fk_213, fk_214, \
                         fk_215, hk_464, hk_465, hk_466, hk_467, hk_468, \
                         hk_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -fk_212[k]
                   + f_0 * hk_464[k];

        t_321[k] = -fk_213[k]
                   + f_0 * hk_465[k];

        t_322[k] = -fk_214[k]
                   + f_0 * hk_466[k];

        t_323[k] = -fk_215[k]
                   + f_0 * hk_467[k];

        t_324[k] = f_0 * hk_468[k];

        t_325[k] = f_0 * hk_469[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, t_332, t_333, hk_470, \
                         hk_471, hk_472, hk_473, hk_474, hk_475, hk_476, \
                         hk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_0 * hk_470[k];

        t_327[k] = f_0 * hk_471[k];

        t_328[k] = f_0 * hk_472[k];

        t_329[k] = f_0 * hk_473[k];

        t_330[k] = f_0 * hk_474[k];

        t_331[k] = f_0 * hk_475[k];

        t_332[k] = f_0 * hk_476[k];

        t_333[k] = f_0 * hk_477[k];
    }
}

static auto
compute_prim_geom_10_gk_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t hk,
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

    const auto *fk_216 = buffer.data(fk + 216);
    const auto *fk_217 = buffer.data(fk + 217);
    const auto *fk_218 = buffer.data(fk + 218);
    const auto *fk_219 = buffer.data(fk + 219);
    const auto *fk_220 = buffer.data(fk + 220);
    const auto *fk_221 = buffer.data(fk + 221);
    const auto *fk_222 = buffer.data(fk + 222);
    const auto *fk_223 = buffer.data(fk + 223);
    const auto *fk_224 = buffer.data(fk + 224);
    const auto *fk_225 = buffer.data(fk + 225);
    const auto *fk_226 = buffer.data(fk + 226);
    const auto *fk_227 = buffer.data(fk + 227);
    const auto *fk_228 = buffer.data(fk + 228);
    const auto *fk_229 = buffer.data(fk + 229);
    const auto *fk_230 = buffer.data(fk + 230);
    const auto *fk_231 = buffer.data(fk + 231);
    const auto *fk_232 = buffer.data(fk + 232);
    const auto *fk_233 = buffer.data(fk + 233);
    const auto *fk_234 = buffer.data(fk + 234);
    const auto *fk_235 = buffer.data(fk + 235);
    const auto *fk_236 = buffer.data(fk + 236);
    const auto *fk_237 = buffer.data(fk + 237);
    const auto *fk_238 = buffer.data(fk + 238);
    const auto *fk_239 = buffer.data(fk + 239);
    const auto *fk_240 = buffer.data(fk + 240);
    const auto *fk_241 = buffer.data(fk + 241);
    const auto *fk_242 = buffer.data(fk + 242);
    const auto *fk_243 = buffer.data(fk + 243);
    const auto *fk_244 = buffer.data(fk + 244);
    const auto *fk_245 = buffer.data(fk + 245);
    const auto *fk_246 = buffer.data(fk + 246);
    const auto *fk_247 = buffer.data(fk + 247);
    const auto *fk_248 = buffer.data(fk + 248);
    const auto *fk_249 = buffer.data(fk + 249);
    const auto *fk_250 = buffer.data(fk + 250);
    const auto *fk_251 = buffer.data(fk + 251);
    const auto *fk_252 = buffer.data(fk + 252);
    const auto *fk_253 = buffer.data(fk + 253);
    const auto *fk_254 = buffer.data(fk + 254);
    const auto *fk_255 = buffer.data(fk + 255);
    const auto *fk_256 = buffer.data(fk + 256);
    const auto *fk_257 = buffer.data(fk + 257);
    const auto *fk_258 = buffer.data(fk + 258);
    const auto *fk_259 = buffer.data(fk + 259);
    const auto *fk_260 = buffer.data(fk + 260);
    const auto *fk_261 = buffer.data(fk + 261);
    const auto *fk_262 = buffer.data(fk + 262);
    const auto *fk_263 = buffer.data(fk + 263);
    const auto *fk_264 = buffer.data(fk + 264);
    const auto *fk_265 = buffer.data(fk + 265);
    const auto *fk_266 = buffer.data(fk + 266);
    const auto *fk_267 = buffer.data(fk + 267);
    const auto *fk_268 = buffer.data(fk + 268);
    const auto *fk_269 = buffer.data(fk + 269);
    const auto *fk_270 = buffer.data(fk + 270);
    const auto *fk_271 = buffer.data(fk + 271);
    const auto *fk_272 = buffer.data(fk + 272);
    const auto *fk_273 = buffer.data(fk + 273);
    const auto *fk_274 = buffer.data(fk + 274);
    const auto *fk_275 = buffer.data(fk + 275);
    const auto *fk_276 = buffer.data(fk + 276);
    const auto *fk_277 = buffer.data(fk + 277);
    const auto *fk_278 = buffer.data(fk + 278);
    const auto *fk_279 = buffer.data(fk + 279);
    const auto *fk_280 = buffer.data(fk + 280);
    const auto *fk_281 = buffer.data(fk + 281);
    const auto *fk_282 = buffer.data(fk + 282);
    const auto *fk_283 = buffer.data(fk + 283);
    const auto *fk_284 = buffer.data(fk + 284);
    const auto *fk_285 = buffer.data(fk + 285);
    const auto *fk_286 = buffer.data(fk + 286);
    const auto *fk_287 = buffer.data(fk + 287);
    const auto *fk_288 = buffer.data(fk + 288);
    const auto *fk_289 = buffer.data(fk + 289);
    const auto *fk_290 = buffer.data(fk + 290);
    const auto *fk_291 = buffer.data(fk + 291);
    const auto *fk_292 = buffer.data(fk + 292);
    const auto *fk_293 = buffer.data(fk + 293);
    const auto *fk_294 = buffer.data(fk + 294);
    const auto *fk_295 = buffer.data(fk + 295);
    const auto *fk_296 = buffer.data(fk + 296);
    const auto *fk_297 = buffer.data(fk + 297);
    const auto *fk_298 = buffer.data(fk + 298);
    const auto *fk_299 = buffer.data(fk + 299);
    const auto *fk_300 = buffer.data(fk + 300);
    const auto *fk_301 = buffer.data(fk + 301);
    const auto *fk_302 = buffer.data(fk + 302);
    const auto *fk_303 = buffer.data(fk + 303);
    const auto *fk_304 = buffer.data(fk + 304);
    const auto *fk_305 = buffer.data(fk + 305);
    const auto *fk_306 = buffer.data(fk + 306);
    const auto *fk_307 = buffer.data(fk + 307);
    const auto *fk_308 = buffer.data(fk + 308);
    const auto *fk_309 = buffer.data(fk + 309);
    const auto *fk_310 = buffer.data(fk + 310);
    const auto *fk_311 = buffer.data(fk + 311);
    const auto *fk_312 = buffer.data(fk + 312);
    const auto *fk_313 = buffer.data(fk + 313);
    const auto *fk_314 = buffer.data(fk + 314);
    const auto *fk_315 = buffer.data(fk + 315);
    const auto *fk_316 = buffer.data(fk + 316);
    const auto *fk_317 = buffer.data(fk + 317);
    const auto *fk_318 = buffer.data(fk + 318);
    const auto *fk_319 = buffer.data(fk + 319);
    const auto *fk_320 = buffer.data(fk + 320);
    const auto *fk_321 = buffer.data(fk + 321);
    const auto *fk_322 = buffer.data(fk + 322);
    const auto *fk_323 = buffer.data(fk + 323);
    const auto *fk_324 = buffer.data(fk + 324);
    const auto *fk_325 = buffer.data(fk + 325);
    const auto *fk_326 = buffer.data(fk + 326);
    const auto *fk_327 = buffer.data(fk + 327);
    const auto *fk_328 = buffer.data(fk + 328);
    const auto *fk_329 = buffer.data(fk + 329);
    const auto *fk_330 = buffer.data(fk + 330);
    const auto *fk_331 = buffer.data(fk + 331);
    const auto *fk_332 = buffer.data(fk + 332);
    const auto *fk_333 = buffer.data(fk + 333);
    const auto *fk_334 = buffer.data(fk + 334);
    const auto *fk_335 = buffer.data(fk + 335);
    const auto *fk_336 = buffer.data(fk + 336);
    const auto *fk_337 = buffer.data(fk + 337);
    const auto *fk_338 = buffer.data(fk + 338);
    const auto *fk_339 = buffer.data(fk + 339);
    const auto *fk_340 = buffer.data(fk + 340);
    const auto *fk_341 = buffer.data(fk + 341);
    const auto *fk_342 = buffer.data(fk + 342);
    const auto *fk_343 = buffer.data(fk + 343);
    const auto *fk_344 = buffer.data(fk + 344);
    const auto *fk_345 = buffer.data(fk + 345);
    const auto *fk_346 = buffer.data(fk + 346);
    const auto *fk_347 = buffer.data(fk + 347);
    const auto *fk_348 = buffer.data(fk + 348);
    const auto *fk_349 = buffer.data(fk + 349);

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

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, t_339, t_340, t_341, hk_478, \
                         hk_479, hk_480, hk_481, hk_482, hk_483, hk_484, \
                         hk_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_0 * hk_478[k];

        t_335[k] = f_0 * hk_479[k];

        t_336[k] = f_0 * hk_480[k];

        t_337[k] = f_0 * hk_481[k];

        t_338[k] = f_0 * hk_482[k];

        t_339[k] = f_0 * hk_483[k];

        t_340[k] = f_0 * hk_484[k];

        t_341[k] = f_0 * hk_485[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, t_347, t_348, t_349, hk_486, \
                         hk_487, hk_488, hk_489, hk_490, hk_491, hk_492, \
                         hk_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_0 * hk_486[k];

        t_343[k] = f_0 * hk_487[k];

        t_344[k] = f_0 * hk_488[k];

        t_345[k] = f_0 * hk_489[k];

        t_346[k] = f_0 * hk_490[k];

        t_347[k] = f_0 * hk_491[k];

        t_348[k] = f_0 * hk_492[k];

        t_349[k] = f_0 * hk_493[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, t_355, t_356, t_357, hk_494, \
                         hk_495, hk_496, hk_497, hk_498, hk_499, hk_500, \
                         hk_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_0 * hk_494[k];

        t_351[k] = f_0 * hk_495[k];

        t_352[k] = f_0 * hk_496[k];

        t_353[k] = f_0 * hk_497[k];

        t_354[k] = f_0 * hk_498[k];

        t_355[k] = f_0 * hk_499[k];

        t_356[k] = f_0 * hk_500[k];

        t_357[k] = f_0 * hk_501[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, t_363, fk_216, fk_217, fk_218, \
                         fk_219, hk_502, hk_503, hk_540, hk_541, hk_542, \
                         hk_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_0 * hk_502[k];

        t_359[k] = f_0 * hk_503[k];

        t_360[k] = -4.0 * fk_216[k]
                   + f_0 * hk_540[k];

        t_361[k] = -4.0 * fk_217[k]
                   + f_0 * hk_541[k];

        t_362[k] = -4.0 * fk_218[k]
                   + f_0 * hk_542[k];

        t_363[k] = -4.0 * fk_219[k]
                   + f_0 * hk_543[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, t_368, fk_220, fk_221, fk_222, fk_223, \
                         fk_224, hk_544, hk_545, hk_546, hk_547, \
                         hk_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = -4.0 * fk_220[k]
                   + f_0 * hk_544[k];

        t_365[k] = -4.0 * fk_221[k]
                   + f_0 * hk_545[k];

        t_366[k] = -4.0 * fk_222[k]
                   + f_0 * hk_546[k];

        t_367[k] = -4.0 * fk_223[k]
                   + f_0 * hk_547[k];

        t_368[k] = -4.0 * fk_224[k]
                   + f_0 * hk_548[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, t_373, fk_225, fk_226, fk_227, fk_228, \
                         fk_229, hk_549, hk_550, hk_551, hk_552, \
                         hk_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = -4.0 * fk_225[k]
                   + f_0 * hk_549[k];

        t_370[k] = -4.0 * fk_226[k]
                   + f_0 * hk_550[k];

        t_371[k] = -4.0 * fk_227[k]
                   + f_0 * hk_551[k];

        t_372[k] = -4.0 * fk_228[k]
                   + f_0 * hk_552[k];

        t_373[k] = -4.0 * fk_229[k]
                   + f_0 * hk_553[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, t_378, fk_230, fk_231, fk_232, fk_233, \
                         fk_234, hk_554, hk_555, hk_556, hk_557, \
                         hk_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = -4.0 * fk_230[k]
                   + f_0 * hk_554[k];

        t_375[k] = -4.0 * fk_231[k]
                   + f_0 * hk_555[k];

        t_376[k] = -4.0 * fk_232[k]
                   + f_0 * hk_556[k];

        t_377[k] = -4.0 * fk_233[k]
                   + f_0 * hk_557[k];

        t_378[k] = -4.0 * fk_234[k]
                   + f_0 * hk_558[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, fk_235, fk_236, fk_237, fk_238, \
                         fk_239, hk_559, hk_560, hk_561, hk_562, \
                         hk_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = -4.0 * fk_235[k]
                   + f_0 * hk_559[k];

        t_380[k] = -4.0 * fk_236[k]
                   + f_0 * hk_560[k];

        t_381[k] = -4.0 * fk_237[k]
                   + f_0 * hk_561[k];

        t_382[k] = -4.0 * fk_238[k]
                   + f_0 * hk_562[k];

        t_383[k] = -4.0 * fk_239[k]
                   + f_0 * hk_563[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, t_388, fk_240, fk_241, fk_242, fk_243, \
                         fk_244, hk_564, hk_565, hk_566, hk_567, \
                         hk_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = -4.0 * fk_240[k]
                   + f_0 * hk_564[k];

        t_385[k] = -4.0 * fk_241[k]
                   + f_0 * hk_565[k];

        t_386[k] = -4.0 * fk_242[k]
                   + f_0 * hk_566[k];

        t_387[k] = -4.0 * fk_243[k]
                   + f_0 * hk_567[k];

        t_388[k] = -4.0 * fk_244[k]
                   + f_0 * hk_568[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, fk_245, fk_246, fk_247, fk_248, \
                         fk_249, hk_569, hk_570, hk_571, hk_572, \
                         hk_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = -4.0 * fk_245[k]
                   + f_0 * hk_569[k];

        t_390[k] = -4.0 * fk_246[k]
                   + f_0 * hk_570[k];

        t_391[k] = -4.0 * fk_247[k]
                   + f_0 * hk_571[k];

        t_392[k] = -4.0 * fk_248[k]
                   + f_0 * hk_572[k];

        t_393[k] = -4.0 * fk_249[k]
                   + f_0 * hk_573[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, fk_250, fk_251, fk_252, fk_253, \
                         fk_254, hk_574, hk_575, hk_576, hk_577, \
                         hk_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = -4.0 * fk_250[k]
                   + f_0 * hk_574[k];

        t_395[k] = -4.0 * fk_251[k]
                   + f_0 * hk_575[k];

        t_396[k] = -3.0 * fk_252[k]
                   + f_0 * hk_576[k];

        t_397[k] = -3.0 * fk_253[k]
                   + f_0 * hk_577[k];

        t_398[k] = -3.0 * fk_254[k]
                   + f_0 * hk_578[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, t_403, fk_255, fk_256, fk_257, fk_258, \
                         fk_259, hk_579, hk_580, hk_581, hk_582, \
                         hk_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = -3.0 * fk_255[k]
                   + f_0 * hk_579[k];

        t_400[k] = -3.0 * fk_256[k]
                   + f_0 * hk_580[k];

        t_401[k] = -3.0 * fk_257[k]
                   + f_0 * hk_581[k];

        t_402[k] = -3.0 * fk_258[k]
                   + f_0 * hk_582[k];

        t_403[k] = -3.0 * fk_259[k]
                   + f_0 * hk_583[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, t_408, fk_260, fk_261, fk_262, fk_263, \
                         fk_264, hk_584, hk_585, hk_586, hk_587, \
                         hk_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = -3.0 * fk_260[k]
                   + f_0 * hk_584[k];

        t_405[k] = -3.0 * fk_261[k]
                   + f_0 * hk_585[k];

        t_406[k] = -3.0 * fk_262[k]
                   + f_0 * hk_586[k];

        t_407[k] = -3.0 * fk_263[k]
                   + f_0 * hk_587[k];

        t_408[k] = -3.0 * fk_264[k]
                   + f_0 * hk_588[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, fk_265, fk_266, fk_267, fk_268, \
                         fk_269, hk_589, hk_590, hk_591, hk_592, \
                         hk_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = -3.0 * fk_265[k]
                   + f_0 * hk_589[k];

        t_410[k] = -3.0 * fk_266[k]
                   + f_0 * hk_590[k];

        t_411[k] = -3.0 * fk_267[k]
                   + f_0 * hk_591[k];

        t_412[k] = -3.0 * fk_268[k]
                   + f_0 * hk_592[k];

        t_413[k] = -3.0 * fk_269[k]
                   + f_0 * hk_593[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, fk_270, fk_271, fk_272, fk_273, \
                         fk_274, hk_594, hk_595, hk_596, hk_597, \
                         hk_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = -3.0 * fk_270[k]
                   + f_0 * hk_594[k];

        t_415[k] = -3.0 * fk_271[k]
                   + f_0 * hk_595[k];

        t_416[k] = -3.0 * fk_272[k]
                   + f_0 * hk_596[k];

        t_417[k] = -3.0 * fk_273[k]
                   + f_0 * hk_597[k];

        t_418[k] = -3.0 * fk_274[k]
                   + f_0 * hk_598[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, t_423, fk_275, fk_276, fk_277, fk_278, \
                         fk_279, hk_599, hk_600, hk_601, hk_602, \
                         hk_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = -3.0 * fk_275[k]
                   + f_0 * hk_599[k];

        t_420[k] = -3.0 * fk_276[k]
                   + f_0 * hk_600[k];

        t_421[k] = -3.0 * fk_277[k]
                   + f_0 * hk_601[k];

        t_422[k] = -3.0 * fk_278[k]
                   + f_0 * hk_602[k];

        t_423[k] = -3.0 * fk_279[k]
                   + f_0 * hk_603[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, t_428, fk_280, fk_281, fk_282, fk_283, \
                         fk_284, hk_604, hk_605, hk_606, hk_607, \
                         hk_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = -3.0 * fk_280[k]
                   + f_0 * hk_604[k];

        t_425[k] = -3.0 * fk_281[k]
                   + f_0 * hk_605[k];

        t_426[k] = -3.0 * fk_282[k]
                   + f_0 * hk_606[k];

        t_427[k] = -3.0 * fk_283[k]
                   + f_0 * hk_607[k];

        t_428[k] = -3.0 * fk_284[k]
                   + f_0 * hk_608[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, fk_285, fk_286, fk_287, fk_288, \
                         fk_289, hk_609, hk_610, hk_611, hk_612, \
                         hk_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = -3.0 * fk_285[k]
                   + f_0 * hk_609[k];

        t_430[k] = -3.0 * fk_286[k]
                   + f_0 * hk_610[k];

        t_431[k] = -3.0 * fk_287[k]
                   + f_0 * hk_611[k];

        t_432[k] = -2.0 * fk_288[k]
                   + f_0 * hk_612[k];

        t_433[k] = -2.0 * fk_289[k]
                   + f_0 * hk_613[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, t_438, fk_290, fk_291, fk_292, fk_293, \
                         fk_294, hk_614, hk_615, hk_616, hk_617, \
                         hk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = -2.0 * fk_290[k]
                   + f_0 * hk_614[k];

        t_435[k] = -2.0 * fk_291[k]
                   + f_0 * hk_615[k];

        t_436[k] = -2.0 * fk_292[k]
                   + f_0 * hk_616[k];

        t_437[k] = -2.0 * fk_293[k]
                   + f_0 * hk_617[k];

        t_438[k] = -2.0 * fk_294[k]
                   + f_0 * hk_618[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, t_443, fk_295, fk_296, fk_297, fk_298, \
                         fk_299, hk_619, hk_620, hk_621, hk_622, \
                         hk_623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = -2.0 * fk_295[k]
                   + f_0 * hk_619[k];

        t_440[k] = -2.0 * fk_296[k]
                   + f_0 * hk_620[k];

        t_441[k] = -2.0 * fk_297[k]
                   + f_0 * hk_621[k];

        t_442[k] = -2.0 * fk_298[k]
                   + f_0 * hk_622[k];

        t_443[k] = -2.0 * fk_299[k]
                   + f_0 * hk_623[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, fk_300, fk_301, fk_302, fk_303, \
                         fk_304, hk_624, hk_625, hk_626, hk_627, \
                         hk_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = -2.0 * fk_300[k]
                   + f_0 * hk_624[k];

        t_445[k] = -2.0 * fk_301[k]
                   + f_0 * hk_625[k];

        t_446[k] = -2.0 * fk_302[k]
                   + f_0 * hk_626[k];

        t_447[k] = -2.0 * fk_303[k]
                   + f_0 * hk_627[k];

        t_448[k] = -2.0 * fk_304[k]
                   + f_0 * hk_628[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, fk_305, fk_306, fk_307, fk_308, \
                         fk_309, hk_629, hk_630, hk_631, hk_632, \
                         hk_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = -2.0 * fk_305[k]
                   + f_0 * hk_629[k];

        t_450[k] = -2.0 * fk_306[k]
                   + f_0 * hk_630[k];

        t_451[k] = -2.0 * fk_307[k]
                   + f_0 * hk_631[k];

        t_452[k] = -2.0 * fk_308[k]
                   + f_0 * hk_632[k];

        t_453[k] = -2.0 * fk_309[k]
                   + f_0 * hk_633[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, fk_310, fk_311, fk_312, fk_313, \
                         fk_314, hk_634, hk_635, hk_636, hk_637, \
                         hk_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = -2.0 * fk_310[k]
                   + f_0 * hk_634[k];

        t_455[k] = -2.0 * fk_311[k]
                   + f_0 * hk_635[k];

        t_456[k] = -2.0 * fk_312[k]
                   + f_0 * hk_636[k];

        t_457[k] = -2.0 * fk_313[k]
                   + f_0 * hk_637[k];

        t_458[k] = -2.0 * fk_314[k]
                   + f_0 * hk_638[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, t_463, fk_315, fk_316, fk_317, fk_318, \
                         fk_319, hk_639, hk_640, hk_641, hk_642, \
                         hk_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = -2.0 * fk_315[k]
                   + f_0 * hk_639[k];

        t_460[k] = -2.0 * fk_316[k]
                   + f_0 * hk_640[k];

        t_461[k] = -2.0 * fk_317[k]
                   + f_0 * hk_641[k];

        t_462[k] = -2.0 * fk_318[k]
                   + f_0 * hk_642[k];

        t_463[k] = -2.0 * fk_319[k]
                   + f_0 * hk_643[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, t_468, fk_320, fk_321, fk_322, fk_323, \
                         fk_324, hk_644, hk_645, hk_646, hk_647, \
                         hk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = -2.0 * fk_320[k]
                   + f_0 * hk_644[k];

        t_465[k] = -2.0 * fk_321[k]
                   + f_0 * hk_645[k];

        t_466[k] = -2.0 * fk_322[k]
                   + f_0 * hk_646[k];

        t_467[k] = -2.0 * fk_323[k]
                   + f_0 * hk_647[k];

        t_468[k] = -fk_324[k]
                   + f_0 * hk_648[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, t_473, fk_325, fk_326, fk_327, fk_328, \
                         fk_329, hk_649, hk_650, hk_651, hk_652, \
                         hk_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = -fk_325[k]
                   + f_0 * hk_649[k];

        t_470[k] = -fk_326[k]
                   + f_0 * hk_650[k];

        t_471[k] = -fk_327[k]
                   + f_0 * hk_651[k];

        t_472[k] = -fk_328[k]
                   + f_0 * hk_652[k];

        t_473[k] = -fk_329[k]
                   + f_0 * hk_653[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, t_478, fk_330, fk_331, fk_332, fk_333, \
                         fk_334, hk_654, hk_655, hk_656, hk_657, \
                         hk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = -fk_330[k]
                   + f_0 * hk_654[k];

        t_475[k] = -fk_331[k]
                   + f_0 * hk_655[k];

        t_476[k] = -fk_332[k]
                   + f_0 * hk_656[k];

        t_477[k] = -fk_333[k]
                   + f_0 * hk_657[k];

        t_478[k] = -fk_334[k]
                   + f_0 * hk_658[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, t_483, fk_335, fk_336, fk_337, fk_338, \
                         fk_339, hk_659, hk_660, hk_661, hk_662, \
                         hk_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = -fk_335[k]
                   + f_0 * hk_659[k];

        t_480[k] = -fk_336[k]
                   + f_0 * hk_660[k];

        t_481[k] = -fk_337[k]
                   + f_0 * hk_661[k];

        t_482[k] = -fk_338[k]
                   + f_0 * hk_662[k];

        t_483[k] = -fk_339[k]
                   + f_0 * hk_663[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, t_488, fk_340, fk_341, fk_342, fk_343, \
                         fk_344, hk_664, hk_665, hk_666, hk_667, \
                         hk_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = -fk_340[k]
                   + f_0 * hk_664[k];

        t_485[k] = -fk_341[k]
                   + f_0 * hk_665[k];

        t_486[k] = -fk_342[k]
                   + f_0 * hk_666[k];

        t_487[k] = -fk_343[k]
                   + f_0 * hk_667[k];

        t_488[k] = -fk_344[k]
                   + f_0 * hk_668[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, fk_345, fk_346, fk_347, fk_348, \
                         fk_349, hk_669, hk_670, hk_671, hk_672, \
                         hk_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = -fk_345[k]
                   + f_0 * hk_669[k];

        t_490[k] = -fk_346[k]
                   + f_0 * hk_670[k];

        t_491[k] = -fk_347[k]
                   + f_0 * hk_671[k];

        t_492[k] = -fk_348[k]
                   + f_0 * hk_672[k];

        t_493[k] = -fk_349[k]
                   + f_0 * hk_673[k];
    }
}

static auto
compute_prim_geom_10_gk_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t hk,
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

    const auto *fk_350 = buffer.data(fk + 350);
    const auto *fk_351 = buffer.data(fk + 351);
    const auto *fk_352 = buffer.data(fk + 352);
    const auto *fk_353 = buffer.data(fk + 353);
    const auto *fk_354 = buffer.data(fk + 354);
    const auto *fk_355 = buffer.data(fk + 355);
    const auto *fk_356 = buffer.data(fk + 356);
    const auto *fk_357 = buffer.data(fk + 357);
    const auto *fk_358 = buffer.data(fk + 358);
    const auto *fk_359 = buffer.data(fk + 359);

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

#pragma omp simd aligned(t_494, t_495, t_496, t_497, t_498, fk_350, fk_351, fk_352, fk_353, \
                         fk_354, hk_674, hk_675, hk_676, hk_677, \
                         hk_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = -fk_350[k]
                   + f_0 * hk_674[k];

        t_495[k] = -fk_351[k]
                   + f_0 * hk_675[k];

        t_496[k] = -fk_352[k]
                   + f_0 * hk_676[k];

        t_497[k] = -fk_353[k]
                   + f_0 * hk_677[k];

        t_498[k] = -fk_354[k]
                   + f_0 * hk_678[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, t_503, fk_355, fk_356, fk_357, fk_358, \
                         fk_359, hk_679, hk_680, hk_681, hk_682, \
                         hk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = -fk_355[k]
                   + f_0 * hk_679[k];

        t_500[k] = -fk_356[k]
                   + f_0 * hk_680[k];

        t_501[k] = -fk_357[k]
                   + f_0 * hk_681[k];

        t_502[k] = -fk_358[k]
                   + f_0 * hk_682[k];

        t_503[k] = -fk_359[k]
                   + f_0 * hk_683[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, t_509, t_510, t_511, hk_684, \
                         hk_685, hk_686, hk_687, hk_688, hk_689, hk_690, \
                         hk_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_0 * hk_684[k];

        t_505[k] = f_0 * hk_685[k];

        t_506[k] = f_0 * hk_686[k];

        t_507[k] = f_0 * hk_687[k];

        t_508[k] = f_0 * hk_688[k];

        t_509[k] = f_0 * hk_689[k];

        t_510[k] = f_0 * hk_690[k];

        t_511[k] = f_0 * hk_691[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, t_517, t_518, t_519, hk_692, \
                         hk_693, hk_694, hk_695, hk_696, hk_697, hk_698, \
                         hk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_0 * hk_692[k];

        t_513[k] = f_0 * hk_693[k];

        t_514[k] = f_0 * hk_694[k];

        t_515[k] = f_0 * hk_695[k];

        t_516[k] = f_0 * hk_696[k];

        t_517[k] = f_0 * hk_697[k];

        t_518[k] = f_0 * hk_698[k];

        t_519[k] = f_0 * hk_699[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, t_525, t_526, t_527, hk_700, \
                         hk_701, hk_702, hk_703, hk_704, hk_705, hk_706, \
                         hk_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_0 * hk_700[k];

        t_521[k] = f_0 * hk_701[k];

        t_522[k] = f_0 * hk_702[k];

        t_523[k] = f_0 * hk_703[k];

        t_524[k] = f_0 * hk_704[k];

        t_525[k] = f_0 * hk_705[k];

        t_526[k] = f_0 * hk_706[k];

        t_527[k] = f_0 * hk_707[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, t_533, t_534, t_535, hk_708, \
                         hk_709, hk_710, hk_711, hk_712, hk_713, hk_714, \
                         hk_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_0 * hk_708[k];

        t_529[k] = f_0 * hk_709[k];

        t_530[k] = f_0 * hk_710[k];

        t_531[k] = f_0 * hk_711[k];

        t_532[k] = f_0 * hk_712[k];

        t_533[k] = f_0 * hk_713[k];

        t_534[k] = f_0 * hk_714[k];

        t_535[k] = f_0 * hk_715[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, hk_716, hk_717, hk_718, \
                         hk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_0 * hk_716[k];

        t_537[k] = f_0 * hk_717[k];

        t_538[k] = f_0 * hk_718[k];

        t_539[k] = f_0 * hk_719[k];
    }
}

auto
compute_prim_geom_10_gk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t fk, const size_t hk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gk_electron_repulsion_1_piece0(buffer, target, fk, hk, ncols, alpha);

    compute_prim_geom_10_gk_electron_repulsion_1_piece1(buffer, target, fk, hk, ncols, alpha);

    compute_prim_geom_10_gk_electron_repulsion_1_piece2(buffer, target, fk, hk, ncols, alpha);

    compute_prim_geom_10_gk_electron_repulsion_1_piece3(buffer, target, fk, hk, ncols, alpha);
}

static auto
compute_prim_geom_10_gk_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t hk,
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

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_29 = buffer.data(fk + 29);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_34 = buffer.data(fk + 34);
    const auto *fk_35 = buffer.data(fk + 35);
    const auto *fk_36 = buffer.data(fk + 36);
    const auto *fk_37 = buffer.data(fk + 37);
    const auto *fk_38 = buffer.data(fk + 38);
    const auto *fk_39 = buffer.data(fk + 39);
    const auto *fk_40 = buffer.data(fk + 40);
    const auto *fk_41 = buffer.data(fk + 41);
    const auto *fk_42 = buffer.data(fk + 42);
    const auto *fk_43 = buffer.data(fk + 43);
    const auto *fk_44 = buffer.data(fk + 44);
    const auto *fk_45 = buffer.data(fk + 45);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_47 = buffer.data(fk + 47);
    const auto *fk_48 = buffer.data(fk + 48);
    const auto *fk_49 = buffer.data(fk + 49);
    const auto *fk_50 = buffer.data(fk + 50);
    const auto *fk_51 = buffer.data(fk + 51);
    const auto *fk_52 = buffer.data(fk + 52);
    const auto *fk_53 = buffer.data(fk + 53);
    const auto *fk_54 = buffer.data(fk + 54);
    const auto *fk_55 = buffer.data(fk + 55);
    const auto *fk_56 = buffer.data(fk + 56);
    const auto *fk_57 = buffer.data(fk + 57);
    const auto *fk_58 = buffer.data(fk + 58);
    const auto *fk_59 = buffer.data(fk + 59);
    const auto *fk_60 = buffer.data(fk + 60);
    const auto *fk_61 = buffer.data(fk + 61);
    const auto *fk_62 = buffer.data(fk + 62);
    const auto *fk_63 = buffer.data(fk + 63);
    const auto *fk_64 = buffer.data(fk + 64);
    const auto *fk_65 = buffer.data(fk + 65);
    const auto *fk_66 = buffer.data(fk + 66);
    const auto *fk_67 = buffer.data(fk + 67);
    const auto *fk_68 = buffer.data(fk + 68);
    const auto *fk_69 = buffer.data(fk + 69);
    const auto *fk_70 = buffer.data(fk + 70);
    const auto *fk_71 = buffer.data(fk + 71);
    const auto *fk_72 = buffer.data(fk + 72);
    const auto *fk_73 = buffer.data(fk + 73);
    const auto *fk_74 = buffer.data(fk + 74);
    const auto *fk_75 = buffer.data(fk + 75);
    const auto *fk_76 = buffer.data(fk + 76);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hk_72, hk_73, hk_74, hk_75, \
                         hk_76, hk_77, hk_78, hk_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hk_72[k];

        t_1[k] = f_0 * hk_73[k];

        t_2[k] = f_0 * hk_74[k];

        t_3[k] = f_0 * hk_75[k];

        t_4[k] = f_0 * hk_76[k];

        t_5[k] = f_0 * hk_77[k];

        t_6[k] = f_0 * hk_78[k];

        t_7[k] = f_0 * hk_79[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, hk_80, hk_81, hk_82, \
                         hk_83, hk_84, hk_85, hk_86, hk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hk_80[k];

        t_9[k] = f_0 * hk_81[k];

        t_10[k] = f_0 * hk_82[k];

        t_11[k] = f_0 * hk_83[k];

        t_12[k] = f_0 * hk_84[k];

        t_13[k] = f_0 * hk_85[k];

        t_14[k] = f_0 * hk_86[k];

        t_15[k] = f_0 * hk_87[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, hk_88, hk_89, hk_90, \
                         hk_91, hk_92, hk_93, hk_94, hk_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * hk_88[k];

        t_17[k] = f_0 * hk_89[k];

        t_18[k] = f_0 * hk_90[k];

        t_19[k] = f_0 * hk_91[k];

        t_20[k] = f_0 * hk_92[k];

        t_21[k] = f_0 * hk_93[k];

        t_22[k] = f_0 * hk_94[k];

        t_23[k] = f_0 * hk_95[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, hk_96, hk_97, hk_98, \
                         hk_99, hk_100, hk_101, hk_102, hk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * hk_96[k];

        t_25[k] = f_0 * hk_97[k];

        t_26[k] = f_0 * hk_98[k];

        t_27[k] = f_0 * hk_99[k];

        t_28[k] = f_0 * hk_100[k];

        t_29[k] = f_0 * hk_101[k];

        t_30[k] = f_0 * hk_102[k];

        t_31[k] = f_0 * hk_103[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, hk_104, hk_105, \
                         hk_106, hk_107, hk_144, hk_145, hk_146, \
                         hk_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * hk_104[k];

        t_33[k] = f_0 * hk_105[k];

        t_34[k] = f_0 * hk_106[k];

        t_35[k] = f_0 * hk_107[k];

        t_36[k] = f_0 * hk_144[k];

        t_37[k] = f_0 * hk_145[k];

        t_38[k] = f_0 * hk_146[k];

        t_39[k] = f_0 * hk_147[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, hk_148, hk_149, \
                         hk_150, hk_151, hk_152, hk_153, hk_154, \
                         hk_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * hk_148[k];

        t_41[k] = f_0 * hk_149[k];

        t_42[k] = f_0 * hk_150[k];

        t_43[k] = f_0 * hk_151[k];

        t_44[k] = f_0 * hk_152[k];

        t_45[k] = f_0 * hk_153[k];

        t_46[k] = f_0 * hk_154[k];

        t_47[k] = f_0 * hk_155[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, hk_156, hk_157, \
                         hk_158, hk_159, hk_160, hk_161, hk_162, \
                         hk_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * hk_156[k];

        t_49[k] = f_0 * hk_157[k];

        t_50[k] = f_0 * hk_158[k];

        t_51[k] = f_0 * hk_159[k];

        t_52[k] = f_0 * hk_160[k];

        t_53[k] = f_0 * hk_161[k];

        t_54[k] = f_0 * hk_162[k];

        t_55[k] = f_0 * hk_163[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, hk_164, hk_165, \
                         hk_166, hk_167, hk_168, hk_169, hk_170, \
                         hk_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * hk_164[k];

        t_57[k] = f_0 * hk_165[k];

        t_58[k] = f_0 * hk_166[k];

        t_59[k] = f_0 * hk_167[k];

        t_60[k] = f_0 * hk_168[k];

        t_61[k] = f_0 * hk_169[k];

        t_62[k] = f_0 * hk_170[k];

        t_63[k] = f_0 * hk_171[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, hk_172, hk_173, \
                         hk_174, hk_175, hk_176, hk_177, hk_178, \
                         hk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * hk_172[k];

        t_65[k] = f_0 * hk_173[k];

        t_66[k] = f_0 * hk_174[k];

        t_67[k] = f_0 * hk_175[k];

        t_68[k] = f_0 * hk_176[k];

        t_69[k] = f_0 * hk_177[k];

        t_70[k] = f_0 * hk_178[k];

        t_71[k] = f_0 * hk_179[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, fk_0, fk_1, fk_2, fk_3, fk_4, hk_180, \
                         hk_181, hk_182, hk_183, hk_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -fk_0[k]
                  + f_0 * hk_180[k];

        t_73[k] = -fk_1[k]
                  + f_0 * hk_181[k];

        t_74[k] = -fk_2[k]
                  + f_0 * hk_182[k];

        t_75[k] = -fk_3[k]
                  + f_0 * hk_183[k];

        t_76[k] = -fk_4[k]
                  + f_0 * hk_184[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, fk_5, fk_6, fk_7, fk_8, fk_9, hk_185, \
                         hk_186, hk_187, hk_188, hk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -fk_5[k]
                  + f_0 * hk_185[k];

        t_78[k] = -fk_6[k]
                  + f_0 * hk_186[k];

        t_79[k] = -fk_7[k]
                  + f_0 * hk_187[k];

        t_80[k] = -fk_8[k]
                  + f_0 * hk_188[k];

        t_81[k] = -fk_9[k]
                  + f_0 * hk_189[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, fk_10, fk_11, fk_12, fk_13, fk_14, \
                         hk_190, hk_191, hk_192, hk_193, hk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -fk_10[k]
                  + f_0 * hk_190[k];

        t_83[k] = -fk_11[k]
                  + f_0 * hk_191[k];

        t_84[k] = -fk_12[k]
                  + f_0 * hk_192[k];

        t_85[k] = -fk_13[k]
                  + f_0 * hk_193[k];

        t_86[k] = -fk_14[k]
                  + f_0 * hk_194[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, fk_15, fk_16, fk_17, fk_18, fk_19, \
                         hk_195, hk_196, hk_197, hk_198, hk_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -fk_15[k]
                  + f_0 * hk_195[k];

        t_88[k] = -fk_16[k]
                  + f_0 * hk_196[k];

        t_89[k] = -fk_17[k]
                  + f_0 * hk_197[k];

        t_90[k] = -fk_18[k]
                  + f_0 * hk_198[k];

        t_91[k] = -fk_19[k]
                  + f_0 * hk_199[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, fk_20, fk_21, fk_22, fk_23, fk_24, \
                         hk_200, hk_201, hk_202, hk_203, hk_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = -fk_20[k]
                  + f_0 * hk_200[k];

        t_93[k] = -fk_21[k]
                  + f_0 * hk_201[k];

        t_94[k] = -fk_22[k]
                  + f_0 * hk_202[k];

        t_95[k] = -fk_23[k]
                  + f_0 * hk_203[k];

        t_96[k] = -fk_24[k]
                  + f_0 * hk_204[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, fk_25, fk_26, fk_27, fk_28, fk_29, \
                         hk_205, hk_206, hk_207, hk_208, hk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -fk_25[k]
                  + f_0 * hk_205[k];

        t_98[k] = -fk_26[k]
                  + f_0 * hk_206[k];

        t_99[k] = -fk_27[k]
                  + f_0 * hk_207[k];

        t_100[k] = -fk_28[k]
                   + f_0 * hk_208[k];

        t_101[k] = -fk_29[k]
                   + f_0 * hk_209[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, fk_30, fk_31, fk_32, fk_33, fk_34, \
                         hk_210, hk_211, hk_212, hk_213, hk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -fk_30[k]
                   + f_0 * hk_210[k];

        t_103[k] = -fk_31[k]
                   + f_0 * hk_211[k];

        t_104[k] = -fk_32[k]
                   + f_0 * hk_212[k];

        t_105[k] = -fk_33[k]
                   + f_0 * hk_213[k];

        t_106[k] = -fk_34[k]
                   + f_0 * hk_214[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, fk_35, hk_215, \
                         hk_252, hk_253, hk_254, hk_255, hk_256, \
                         hk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -fk_35[k]
                   + f_0 * hk_215[k];

        t_108[k] = f_0 * hk_252[k];

        t_109[k] = f_0 * hk_253[k];

        t_110[k] = f_0 * hk_254[k];

        t_111[k] = f_0 * hk_255[k];

        t_112[k] = f_0 * hk_256[k];

        t_113[k] = f_0 * hk_257[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, t_120, t_121, hk_258, \
                         hk_259, hk_260, hk_261, hk_262, hk_263, hk_264, \
                         hk_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * hk_258[k];

        t_115[k] = f_0 * hk_259[k];

        t_116[k] = f_0 * hk_260[k];

        t_117[k] = f_0 * hk_261[k];

        t_118[k] = f_0 * hk_262[k];

        t_119[k] = f_0 * hk_263[k];

        t_120[k] = f_0 * hk_264[k];

        t_121[k] = f_0 * hk_265[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, t_128, t_129, hk_266, \
                         hk_267, hk_268, hk_269, hk_270, hk_271, hk_272, \
                         hk_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * hk_266[k];

        t_123[k] = f_0 * hk_267[k];

        t_124[k] = f_0 * hk_268[k];

        t_125[k] = f_0 * hk_269[k];

        t_126[k] = f_0 * hk_270[k];

        t_127[k] = f_0 * hk_271[k];

        t_128[k] = f_0 * hk_272[k];

        t_129[k] = f_0 * hk_273[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, t_135, t_136, t_137, hk_274, \
                         hk_275, hk_276, hk_277, hk_278, hk_279, hk_280, \
                         hk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * hk_274[k];

        t_131[k] = f_0 * hk_275[k];

        t_132[k] = f_0 * hk_276[k];

        t_133[k] = f_0 * hk_277[k];

        t_134[k] = f_0 * hk_278[k];

        t_135[k] = f_0 * hk_279[k];

        t_136[k] = f_0 * hk_280[k];

        t_137[k] = f_0 * hk_281[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, t_144, fk_36, hk_282, \
                         hk_283, hk_284, hk_285, hk_286, hk_287, \
                         hk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * hk_282[k];

        t_139[k] = f_0 * hk_283[k];

        t_140[k] = f_0 * hk_284[k];

        t_141[k] = f_0 * hk_285[k];

        t_142[k] = f_0 * hk_286[k];

        t_143[k] = f_0 * hk_287[k];

        t_144[k] = -fk_36[k]
                   + f_0 * hk_288[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, fk_37, fk_38, fk_39, fk_40, fk_41, \
                         hk_289, hk_290, hk_291, hk_292, hk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -fk_37[k]
                   + f_0 * hk_289[k];

        t_146[k] = -fk_38[k]
                   + f_0 * hk_290[k];

        t_147[k] = -fk_39[k]
                   + f_0 * hk_291[k];

        t_148[k] = -fk_40[k]
                   + f_0 * hk_292[k];

        t_149[k] = -fk_41[k]
                   + f_0 * hk_293[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, fk_42, fk_43, fk_44, fk_45, fk_46, \
                         hk_294, hk_295, hk_296, hk_297, hk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -fk_42[k]
                   + f_0 * hk_294[k];

        t_151[k] = -fk_43[k]
                   + f_0 * hk_295[k];

        t_152[k] = -fk_44[k]
                   + f_0 * hk_296[k];

        t_153[k] = -fk_45[k]
                   + f_0 * hk_297[k];

        t_154[k] = -fk_46[k]
                   + f_0 * hk_298[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, fk_47, fk_48, fk_49, fk_50, fk_51, \
                         hk_299, hk_300, hk_301, hk_302, hk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -fk_47[k]
                   + f_0 * hk_299[k];

        t_156[k] = -fk_48[k]
                   + f_0 * hk_300[k];

        t_157[k] = -fk_49[k]
                   + f_0 * hk_301[k];

        t_158[k] = -fk_50[k]
                   + f_0 * hk_302[k];

        t_159[k] = -fk_51[k]
                   + f_0 * hk_303[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, fk_52, fk_53, fk_54, fk_55, fk_56, \
                         hk_304, hk_305, hk_306, hk_307, hk_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -fk_52[k]
                   + f_0 * hk_304[k];

        t_161[k] = -fk_53[k]
                   + f_0 * hk_305[k];

        t_162[k] = -fk_54[k]
                   + f_0 * hk_306[k];

        t_163[k] = -fk_55[k]
                   + f_0 * hk_307[k];

        t_164[k] = -fk_56[k]
                   + f_0 * hk_308[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, fk_57, fk_58, fk_59, fk_60, fk_61, \
                         hk_309, hk_310, hk_311, hk_312, hk_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -fk_57[k]
                   + f_0 * hk_309[k];

        t_166[k] = -fk_58[k]
                   + f_0 * hk_310[k];

        t_167[k] = -fk_59[k]
                   + f_0 * hk_311[k];

        t_168[k] = -fk_60[k]
                   + f_0 * hk_312[k];

        t_169[k] = -fk_61[k]
                   + f_0 * hk_313[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, fk_62, fk_63, fk_64, fk_65, fk_66, \
                         hk_314, hk_315, hk_316, hk_317, hk_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -fk_62[k]
                   + f_0 * hk_314[k];

        t_171[k] = -fk_63[k]
                   + f_0 * hk_315[k];

        t_172[k] = -fk_64[k]
                   + f_0 * hk_316[k];

        t_173[k] = -fk_65[k]
                   + f_0 * hk_317[k];

        t_174[k] = -fk_66[k]
                   + f_0 * hk_318[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, fk_67, fk_68, fk_69, fk_70, fk_71, \
                         hk_319, hk_320, hk_321, hk_322, hk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -fk_67[k]
                   + f_0 * hk_319[k];

        t_176[k] = -fk_68[k]
                   + f_0 * hk_320[k];

        t_177[k] = -fk_69[k]
                   + f_0 * hk_321[k];

        t_178[k] = -fk_70[k]
                   + f_0 * hk_322[k];

        t_179[k] = -fk_71[k]
                   + f_0 * hk_323[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, fk_72, fk_73, fk_74, fk_75, fk_76, \
                         hk_324, hk_325, hk_326, hk_327, hk_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * fk_72[k]
                   + f_0 * hk_324[k];

        t_181[k] = -2.0 * fk_73[k]
                   + f_0 * hk_325[k];

        t_182[k] = -2.0 * fk_74[k]
                   + f_0 * hk_326[k];

        t_183[k] = -2.0 * fk_75[k]
                   + f_0 * hk_327[k];

        t_184[k] = -2.0 * fk_76[k]
                   + f_0 * hk_328[k];
    }
}

static auto
compute_prim_geom_10_gk_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t hk,
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

    const auto *fk_77 = buffer.data(fk + 77);
    const auto *fk_78 = buffer.data(fk + 78);
    const auto *fk_79 = buffer.data(fk + 79);
    const auto *fk_80 = buffer.data(fk + 80);
    const auto *fk_81 = buffer.data(fk + 81);
    const auto *fk_82 = buffer.data(fk + 82);
    const auto *fk_83 = buffer.data(fk + 83);
    const auto *fk_84 = buffer.data(fk + 84);
    const auto *fk_85 = buffer.data(fk + 85);
    const auto *fk_86 = buffer.data(fk + 86);
    const auto *fk_87 = buffer.data(fk + 87);
    const auto *fk_88 = buffer.data(fk + 88);
    const auto *fk_89 = buffer.data(fk + 89);
    const auto *fk_90 = buffer.data(fk + 90);
    const auto *fk_91 = buffer.data(fk + 91);
    const auto *fk_92 = buffer.data(fk + 92);
    const auto *fk_93 = buffer.data(fk + 93);
    const auto *fk_94 = buffer.data(fk + 94);
    const auto *fk_95 = buffer.data(fk + 95);
    const auto *fk_96 = buffer.data(fk + 96);
    const auto *fk_97 = buffer.data(fk + 97);
    const auto *fk_98 = buffer.data(fk + 98);
    const auto *fk_99 = buffer.data(fk + 99);
    const auto *fk_100 = buffer.data(fk + 100);
    const auto *fk_101 = buffer.data(fk + 101);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_106 = buffer.data(fk + 106);
    const auto *fk_107 = buffer.data(fk + 107);
    const auto *fk_108 = buffer.data(fk + 108);
    const auto *fk_109 = buffer.data(fk + 109);
    const auto *fk_110 = buffer.data(fk + 110);
    const auto *fk_111 = buffer.data(fk + 111);
    const auto *fk_112 = buffer.data(fk + 112);
    const auto *fk_113 = buffer.data(fk + 113);
    const auto *fk_114 = buffer.data(fk + 114);
    const auto *fk_115 = buffer.data(fk + 115);
    const auto *fk_116 = buffer.data(fk + 116);
    const auto *fk_117 = buffer.data(fk + 117);
    const auto *fk_118 = buffer.data(fk + 118);
    const auto *fk_119 = buffer.data(fk + 119);
    const auto *fk_120 = buffer.data(fk + 120);
    const auto *fk_121 = buffer.data(fk + 121);
    const auto *fk_122 = buffer.data(fk + 122);
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_124 = buffer.data(fk + 124);
    const auto *fk_125 = buffer.data(fk + 125);
    const auto *fk_126 = buffer.data(fk + 126);
    const auto *fk_127 = buffer.data(fk + 127);
    const auto *fk_128 = buffer.data(fk + 128);
    const auto *fk_129 = buffer.data(fk + 129);
    const auto *fk_130 = buffer.data(fk + 130);
    const auto *fk_131 = buffer.data(fk + 131);
    const auto *fk_132 = buffer.data(fk + 132);
    const auto *fk_133 = buffer.data(fk + 133);
    const auto *fk_134 = buffer.data(fk + 134);
    const auto *fk_135 = buffer.data(fk + 135);
    const auto *fk_136 = buffer.data(fk + 136);
    const auto *fk_137 = buffer.data(fk + 137);
    const auto *fk_138 = buffer.data(fk + 138);
    const auto *fk_139 = buffer.data(fk + 139);
    const auto *fk_140 = buffer.data(fk + 140);
    const auto *fk_141 = buffer.data(fk + 141);
    const auto *fk_142 = buffer.data(fk + 142);
    const auto *fk_143 = buffer.data(fk + 143);
    const auto *fk_144 = buffer.data(fk + 144);
    const auto *fk_145 = buffer.data(fk + 145);
    const auto *fk_146 = buffer.data(fk + 146);
    const auto *fk_147 = buffer.data(fk + 147);
    const auto *fk_148 = buffer.data(fk + 148);
    const auto *fk_149 = buffer.data(fk + 149);
    const auto *fk_150 = buffer.data(fk + 150);
    const auto *fk_151 = buffer.data(fk + 151);
    const auto *fk_152 = buffer.data(fk + 152);
    const auto *fk_153 = buffer.data(fk + 153);
    const auto *fk_154 = buffer.data(fk + 154);
    const auto *fk_155 = buffer.data(fk + 155);
    const auto *fk_156 = buffer.data(fk + 156);
    const auto *fk_157 = buffer.data(fk + 157);
    const auto *fk_158 = buffer.data(fk + 158);
    const auto *fk_159 = buffer.data(fk + 159);
    const auto *fk_160 = buffer.data(fk + 160);
    const auto *fk_161 = buffer.data(fk + 161);
    const auto *fk_162 = buffer.data(fk + 162);
    const auto *fk_163 = buffer.data(fk + 163);
    const auto *fk_164 = buffer.data(fk + 164);
    const auto *fk_165 = buffer.data(fk + 165);
    const auto *fk_166 = buffer.data(fk + 166);
    const auto *fk_167 = buffer.data(fk + 167);
    const auto *fk_168 = buffer.data(fk + 168);
    const auto *fk_169 = buffer.data(fk + 169);
    const auto *fk_170 = buffer.data(fk + 170);
    const auto *fk_171 = buffer.data(fk + 171);
    const auto *fk_172 = buffer.data(fk + 172);
    const auto *fk_173 = buffer.data(fk + 173);
    const auto *fk_174 = buffer.data(fk + 174);
    const auto *fk_175 = buffer.data(fk + 175);
    const auto *fk_176 = buffer.data(fk + 176);
    const auto *fk_177 = buffer.data(fk + 177);
    const auto *fk_178 = buffer.data(fk + 178);
    const auto *fk_179 = buffer.data(fk + 179);
    const auto *fk_180 = buffer.data(fk + 180);
    const auto *fk_181 = buffer.data(fk + 181);
    const auto *fk_182 = buffer.data(fk + 182);
    const auto *fk_183 = buffer.data(fk + 183);
    const auto *fk_184 = buffer.data(fk + 184);
    const auto *fk_185 = buffer.data(fk + 185);
    const auto *fk_186 = buffer.data(fk + 186);
    const auto *fk_187 = buffer.data(fk + 187);
    const auto *fk_188 = buffer.data(fk + 188);
    const auto *fk_189 = buffer.data(fk + 189);
    const auto *fk_190 = buffer.data(fk + 190);
    const auto *fk_191 = buffer.data(fk + 191);
    const auto *fk_192 = buffer.data(fk + 192);
    const auto *fk_193 = buffer.data(fk + 193);
    const auto *fk_194 = buffer.data(fk + 194);
    const auto *fk_195 = buffer.data(fk + 195);
    const auto *fk_196 = buffer.data(fk + 196);
    const auto *fk_197 = buffer.data(fk + 197);
    const auto *fk_198 = buffer.data(fk + 198);
    const auto *fk_199 = buffer.data(fk + 199);
    const auto *fk_200 = buffer.data(fk + 200);
    const auto *fk_201 = buffer.data(fk + 201);
    const auto *fk_202 = buffer.data(fk + 202);
    const auto *fk_203 = buffer.data(fk + 203);

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

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, fk_77, fk_78, fk_79, fk_80, fk_81, \
                         hk_329, hk_330, hk_331, hk_332, hk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * fk_77[k]
                   + f_0 * hk_329[k];

        t_186[k] = -2.0 * fk_78[k]
                   + f_0 * hk_330[k];

        t_187[k] = -2.0 * fk_79[k]
                   + f_0 * hk_331[k];

        t_188[k] = -2.0 * fk_80[k]
                   + f_0 * hk_332[k];

        t_189[k] = -2.0 * fk_81[k]
                   + f_0 * hk_333[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, fk_82, fk_83, fk_84, fk_85, fk_86, \
                         hk_334, hk_335, hk_336, hk_337, hk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * fk_82[k]
                   + f_0 * hk_334[k];

        t_191[k] = -2.0 * fk_83[k]
                   + f_0 * hk_335[k];

        t_192[k] = -2.0 * fk_84[k]
                   + f_0 * hk_336[k];

        t_193[k] = -2.0 * fk_85[k]
                   + f_0 * hk_337[k];

        t_194[k] = -2.0 * fk_86[k]
                   + f_0 * hk_338[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, fk_87, fk_88, fk_89, fk_90, fk_91, \
                         hk_339, hk_340, hk_341, hk_342, hk_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * fk_87[k]
                   + f_0 * hk_339[k];

        t_196[k] = -2.0 * fk_88[k]
                   + f_0 * hk_340[k];

        t_197[k] = -2.0 * fk_89[k]
                   + f_0 * hk_341[k];

        t_198[k] = -2.0 * fk_90[k]
                   + f_0 * hk_342[k];

        t_199[k] = -2.0 * fk_91[k]
                   + f_0 * hk_343[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, fk_92, fk_93, fk_94, fk_95, fk_96, \
                         hk_344, hk_345, hk_346, hk_347, hk_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * fk_92[k]
                   + f_0 * hk_344[k];

        t_201[k] = -2.0 * fk_93[k]
                   + f_0 * hk_345[k];

        t_202[k] = -2.0 * fk_94[k]
                   + f_0 * hk_346[k];

        t_203[k] = -2.0 * fk_95[k]
                   + f_0 * hk_347[k];

        t_204[k] = -2.0 * fk_96[k]
                   + f_0 * hk_348[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, fk_97, fk_98, fk_99, fk_100, \
                         fk_101, hk_349, hk_350, hk_351, hk_352, \
                         hk_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * fk_97[k]
                   + f_0 * hk_349[k];

        t_206[k] = -2.0 * fk_98[k]
                   + f_0 * hk_350[k];

        t_207[k] = -2.0 * fk_99[k]
                   + f_0 * hk_351[k];

        t_208[k] = -2.0 * fk_100[k]
                   + f_0 * hk_352[k];

        t_209[k] = -2.0 * fk_101[k]
                   + f_0 * hk_353[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, fk_102, fk_103, fk_104, fk_105, \
                         fk_106, hk_354, hk_355, hk_356, hk_357, \
                         hk_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * fk_102[k]
                   + f_0 * hk_354[k];

        t_211[k] = -2.0 * fk_103[k]
                   + f_0 * hk_355[k];

        t_212[k] = -2.0 * fk_104[k]
                   + f_0 * hk_356[k];

        t_213[k] = -2.0 * fk_105[k]
                   + f_0 * hk_357[k];

        t_214[k] = -2.0 * fk_106[k]
                   + f_0 * hk_358[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, t_220, t_221, fk_107, hk_359, \
                         hk_396, hk_397, hk_398, hk_399, hk_400, \
                         hk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * fk_107[k]
                   + f_0 * hk_359[k];

        t_216[k] = f_0 * hk_396[k];

        t_217[k] = f_0 * hk_397[k];

        t_218[k] = f_0 * hk_398[k];

        t_219[k] = f_0 * hk_399[k];

        t_220[k] = f_0 * hk_400[k];

        t_221[k] = f_0 * hk_401[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, t_228, t_229, hk_402, \
                         hk_403, hk_404, hk_405, hk_406, hk_407, hk_408, \
                         hk_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_0 * hk_402[k];

        t_223[k] = f_0 * hk_403[k];

        t_224[k] = f_0 * hk_404[k];

        t_225[k] = f_0 * hk_405[k];

        t_226[k] = f_0 * hk_406[k];

        t_227[k] = f_0 * hk_407[k];

        t_228[k] = f_0 * hk_408[k];

        t_229[k] = f_0 * hk_409[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, t_235, t_236, t_237, hk_410, \
                         hk_411, hk_412, hk_413, hk_414, hk_415, hk_416, \
                         hk_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_0 * hk_410[k];

        t_231[k] = f_0 * hk_411[k];

        t_232[k] = f_0 * hk_412[k];

        t_233[k] = f_0 * hk_413[k];

        t_234[k] = f_0 * hk_414[k];

        t_235[k] = f_0 * hk_415[k];

        t_236[k] = f_0 * hk_416[k];

        t_237[k] = f_0 * hk_417[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, t_243, t_244, t_245, hk_418, \
                         hk_419, hk_420, hk_421, hk_422, hk_423, hk_424, \
                         hk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_0 * hk_418[k];

        t_239[k] = f_0 * hk_419[k];

        t_240[k] = f_0 * hk_420[k];

        t_241[k] = f_0 * hk_421[k];

        t_242[k] = f_0 * hk_422[k];

        t_243[k] = f_0 * hk_423[k];

        t_244[k] = f_0 * hk_424[k];

        t_245[k] = f_0 * hk_425[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, t_251, t_252, fk_108, hk_426, \
                         hk_427, hk_428, hk_429, hk_430, hk_431, \
                         hk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_0 * hk_426[k];

        t_247[k] = f_0 * hk_427[k];

        t_248[k] = f_0 * hk_428[k];

        t_249[k] = f_0 * hk_429[k];

        t_250[k] = f_0 * hk_430[k];

        t_251[k] = f_0 * hk_431[k];

        t_252[k] = -fk_108[k]
                   + f_0 * hk_432[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, fk_109, fk_110, fk_111, fk_112, \
                         fk_113, hk_433, hk_434, hk_435, hk_436, \
                         hk_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = -fk_109[k]
                   + f_0 * hk_433[k];

        t_254[k] = -fk_110[k]
                   + f_0 * hk_434[k];

        t_255[k] = -fk_111[k]
                   + f_0 * hk_435[k];

        t_256[k] = -fk_112[k]
                   + f_0 * hk_436[k];

        t_257[k] = -fk_113[k]
                   + f_0 * hk_437[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, fk_114, fk_115, fk_116, fk_117, \
                         fk_118, hk_438, hk_439, hk_440, hk_441, \
                         hk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = -fk_114[k]
                   + f_0 * hk_438[k];

        t_259[k] = -fk_115[k]
                   + f_0 * hk_439[k];

        t_260[k] = -fk_116[k]
                   + f_0 * hk_440[k];

        t_261[k] = -fk_117[k]
                   + f_0 * hk_441[k];

        t_262[k] = -fk_118[k]
                   + f_0 * hk_442[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, fk_119, fk_120, fk_121, fk_122, \
                         fk_123, hk_443, hk_444, hk_445, hk_446, \
                         hk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = -fk_119[k]
                   + f_0 * hk_443[k];

        t_264[k] = -fk_120[k]
                   + f_0 * hk_444[k];

        t_265[k] = -fk_121[k]
                   + f_0 * hk_445[k];

        t_266[k] = -fk_122[k]
                   + f_0 * hk_446[k];

        t_267[k] = -fk_123[k]
                   + f_0 * hk_447[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, fk_124, fk_125, fk_126, fk_127, \
                         fk_128, hk_448, hk_449, hk_450, hk_451, \
                         hk_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = -fk_124[k]
                   + f_0 * hk_448[k];

        t_269[k] = -fk_125[k]
                   + f_0 * hk_449[k];

        t_270[k] = -fk_126[k]
                   + f_0 * hk_450[k];

        t_271[k] = -fk_127[k]
                   + f_0 * hk_451[k];

        t_272[k] = -fk_128[k]
                   + f_0 * hk_452[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, fk_129, fk_130, fk_131, fk_132, \
                         fk_133, hk_453, hk_454, hk_455, hk_456, \
                         hk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = -fk_129[k]
                   + f_0 * hk_453[k];

        t_274[k] = -fk_130[k]
                   + f_0 * hk_454[k];

        t_275[k] = -fk_131[k]
                   + f_0 * hk_455[k];

        t_276[k] = -fk_132[k]
                   + f_0 * hk_456[k];

        t_277[k] = -fk_133[k]
                   + f_0 * hk_457[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, fk_134, fk_135, fk_136, fk_137, \
                         fk_138, hk_458, hk_459, hk_460, hk_461, \
                         hk_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = -fk_134[k]
                   + f_0 * hk_458[k];

        t_279[k] = -fk_135[k]
                   + f_0 * hk_459[k];

        t_280[k] = -fk_136[k]
                   + f_0 * hk_460[k];

        t_281[k] = -fk_137[k]
                   + f_0 * hk_461[k];

        t_282[k] = -fk_138[k]
                   + f_0 * hk_462[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, fk_139, fk_140, fk_141, fk_142, \
                         fk_143, hk_463, hk_464, hk_465, hk_466, \
                         hk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = -fk_139[k]
                   + f_0 * hk_463[k];

        t_284[k] = -fk_140[k]
                   + f_0 * hk_464[k];

        t_285[k] = -fk_141[k]
                   + f_0 * hk_465[k];

        t_286[k] = -fk_142[k]
                   + f_0 * hk_466[k];

        t_287[k] = -fk_143[k]
                   + f_0 * hk_467[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, fk_144, fk_145, fk_146, fk_147, \
                         fk_148, hk_468, hk_469, hk_470, hk_471, \
                         hk_472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = -2.0 * fk_144[k]
                   + f_0 * hk_468[k];

        t_289[k] = -2.0 * fk_145[k]
                   + f_0 * hk_469[k];

        t_290[k] = -2.0 * fk_146[k]
                   + f_0 * hk_470[k];

        t_291[k] = -2.0 * fk_147[k]
                   + f_0 * hk_471[k];

        t_292[k] = -2.0 * fk_148[k]
                   + f_0 * hk_472[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, t_297, fk_149, fk_150, fk_151, fk_152, \
                         fk_153, hk_473, hk_474, hk_475, hk_476, \
                         hk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = -2.0 * fk_149[k]
                   + f_0 * hk_473[k];

        t_294[k] = -2.0 * fk_150[k]
                   + f_0 * hk_474[k];

        t_295[k] = -2.0 * fk_151[k]
                   + f_0 * hk_475[k];

        t_296[k] = -2.0 * fk_152[k]
                   + f_0 * hk_476[k];

        t_297[k] = -2.0 * fk_153[k]
                   + f_0 * hk_477[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, fk_154, fk_155, fk_156, fk_157, \
                         fk_158, hk_478, hk_479, hk_480, hk_481, \
                         hk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = -2.0 * fk_154[k]
                   + f_0 * hk_478[k];

        t_299[k] = -2.0 * fk_155[k]
                   + f_0 * hk_479[k];

        t_300[k] = -2.0 * fk_156[k]
                   + f_0 * hk_480[k];

        t_301[k] = -2.0 * fk_157[k]
                   + f_0 * hk_481[k];

        t_302[k] = -2.0 * fk_158[k]
                   + f_0 * hk_482[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, fk_159, fk_160, fk_161, fk_162, \
                         fk_163, hk_483, hk_484, hk_485, hk_486, \
                         hk_487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = -2.0 * fk_159[k]
                   + f_0 * hk_483[k];

        t_304[k] = -2.0 * fk_160[k]
                   + f_0 * hk_484[k];

        t_305[k] = -2.0 * fk_161[k]
                   + f_0 * hk_485[k];

        t_306[k] = -2.0 * fk_162[k]
                   + f_0 * hk_486[k];

        t_307[k] = -2.0 * fk_163[k]
                   + f_0 * hk_487[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, t_312, fk_164, fk_165, fk_166, fk_167, \
                         fk_168, hk_488, hk_489, hk_490, hk_491, \
                         hk_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = -2.0 * fk_164[k]
                   + f_0 * hk_488[k];

        t_309[k] = -2.0 * fk_165[k]
                   + f_0 * hk_489[k];

        t_310[k] = -2.0 * fk_166[k]
                   + f_0 * hk_490[k];

        t_311[k] = -2.0 * fk_167[k]
                   + f_0 * hk_491[k];

        t_312[k] = -2.0 * fk_168[k]
                   + f_0 * hk_492[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, fk_169, fk_170, fk_171, fk_172, \
                         fk_173, hk_493, hk_494, hk_495, hk_496, \
                         hk_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = -2.0 * fk_169[k]
                   + f_0 * hk_493[k];

        t_314[k] = -2.0 * fk_170[k]
                   + f_0 * hk_494[k];

        t_315[k] = -2.0 * fk_171[k]
                   + f_0 * hk_495[k];

        t_316[k] = -2.0 * fk_172[k]
                   + f_0 * hk_496[k];

        t_317[k] = -2.0 * fk_173[k]
                   + f_0 * hk_497[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, fk_174, fk_175, fk_176, fk_177, \
                         fk_178, hk_498, hk_499, hk_500, hk_501, \
                         hk_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = -2.0 * fk_174[k]
                   + f_0 * hk_498[k];

        t_319[k] = -2.0 * fk_175[k]
                   + f_0 * hk_499[k];

        t_320[k] = -2.0 * fk_176[k]
                   + f_0 * hk_500[k];

        t_321[k] = -2.0 * fk_177[k]
                   + f_0 * hk_501[k];

        t_322[k] = -2.0 * fk_178[k]
                   + f_0 * hk_502[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, fk_179, fk_180, fk_181, fk_182, \
                         fk_183, hk_503, hk_504, hk_505, hk_506, \
                         hk_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = -2.0 * fk_179[k]
                   + f_0 * hk_503[k];

        t_324[k] = -3.0 * fk_180[k]
                   + f_0 * hk_504[k];

        t_325[k] = -3.0 * fk_181[k]
                   + f_0 * hk_505[k];

        t_326[k] = -3.0 * fk_182[k]
                   + f_0 * hk_506[k];

        t_327[k] = -3.0 * fk_183[k]
                   + f_0 * hk_507[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, fk_184, fk_185, fk_186, fk_187, \
                         fk_188, hk_508, hk_509, hk_510, hk_511, \
                         hk_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = -3.0 * fk_184[k]
                   + f_0 * hk_508[k];

        t_329[k] = -3.0 * fk_185[k]
                   + f_0 * hk_509[k];

        t_330[k] = -3.0 * fk_186[k]
                   + f_0 * hk_510[k];

        t_331[k] = -3.0 * fk_187[k]
                   + f_0 * hk_511[k];

        t_332[k] = -3.0 * fk_188[k]
                   + f_0 * hk_512[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, t_337, fk_189, fk_190, fk_191, fk_192, \
                         fk_193, hk_513, hk_514, hk_515, hk_516, \
                         hk_517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = -3.0 * fk_189[k]
                   + f_0 * hk_513[k];

        t_334[k] = -3.0 * fk_190[k]
                   + f_0 * hk_514[k];

        t_335[k] = -3.0 * fk_191[k]
                   + f_0 * hk_515[k];

        t_336[k] = -3.0 * fk_192[k]
                   + f_0 * hk_516[k];

        t_337[k] = -3.0 * fk_193[k]
                   + f_0 * hk_517[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, fk_194, fk_195, fk_196, fk_197, \
                         fk_198, hk_518, hk_519, hk_520, hk_521, \
                         hk_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = -3.0 * fk_194[k]
                   + f_0 * hk_518[k];

        t_339[k] = -3.0 * fk_195[k]
                   + f_0 * hk_519[k];

        t_340[k] = -3.0 * fk_196[k]
                   + f_0 * hk_520[k];

        t_341[k] = -3.0 * fk_197[k]
                   + f_0 * hk_521[k];

        t_342[k] = -3.0 * fk_198[k]
                   + f_0 * hk_522[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, fk_199, fk_200, fk_201, fk_202, \
                         fk_203, hk_523, hk_524, hk_525, hk_526, \
                         hk_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = -3.0 * fk_199[k]
                   + f_0 * hk_523[k];

        t_344[k] = -3.0 * fk_200[k]
                   + f_0 * hk_524[k];

        t_345[k] = -3.0 * fk_201[k]
                   + f_0 * hk_525[k];

        t_346[k] = -3.0 * fk_202[k]
                   + f_0 * hk_526[k];

        t_347[k] = -3.0 * fk_203[k]
                   + f_0 * hk_527[k];
    }
}

static auto
compute_prim_geom_10_gk_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t hk,
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

    const auto *fk_204 = buffer.data(fk + 204);
    const auto *fk_205 = buffer.data(fk + 205);
    const auto *fk_206 = buffer.data(fk + 206);
    const auto *fk_207 = buffer.data(fk + 207);
    const auto *fk_208 = buffer.data(fk + 208);
    const auto *fk_209 = buffer.data(fk + 209);
    const auto *fk_210 = buffer.data(fk + 210);
    const auto *fk_211 = buffer.data(fk + 211);
    const auto *fk_212 = buffer.data(fk + 212);
    const auto *fk_213 = buffer.data(fk + 213);
    const auto *fk_214 = buffer.data(fk + 214);
    const auto *fk_215 = buffer.data(fk + 215);
    const auto *fk_216 = buffer.data(fk + 216);
    const auto *fk_217 = buffer.data(fk + 217);
    const auto *fk_218 = buffer.data(fk + 218);
    const auto *fk_219 = buffer.data(fk + 219);
    const auto *fk_220 = buffer.data(fk + 220);
    const auto *fk_221 = buffer.data(fk + 221);
    const auto *fk_222 = buffer.data(fk + 222);
    const auto *fk_223 = buffer.data(fk + 223);
    const auto *fk_224 = buffer.data(fk + 224);
    const auto *fk_225 = buffer.data(fk + 225);
    const auto *fk_226 = buffer.data(fk + 226);
    const auto *fk_227 = buffer.data(fk + 227);
    const auto *fk_228 = buffer.data(fk + 228);
    const auto *fk_229 = buffer.data(fk + 229);
    const auto *fk_230 = buffer.data(fk + 230);
    const auto *fk_231 = buffer.data(fk + 231);
    const auto *fk_232 = buffer.data(fk + 232);
    const auto *fk_233 = buffer.data(fk + 233);
    const auto *fk_234 = buffer.data(fk + 234);
    const auto *fk_235 = buffer.data(fk + 235);
    const auto *fk_236 = buffer.data(fk + 236);
    const auto *fk_237 = buffer.data(fk + 237);
    const auto *fk_238 = buffer.data(fk + 238);
    const auto *fk_239 = buffer.data(fk + 239);
    const auto *fk_240 = buffer.data(fk + 240);
    const auto *fk_241 = buffer.data(fk + 241);
    const auto *fk_242 = buffer.data(fk + 242);
    const auto *fk_243 = buffer.data(fk + 243);
    const auto *fk_244 = buffer.data(fk + 244);
    const auto *fk_245 = buffer.data(fk + 245);
    const auto *fk_246 = buffer.data(fk + 246);
    const auto *fk_247 = buffer.data(fk + 247);
    const auto *fk_248 = buffer.data(fk + 248);
    const auto *fk_249 = buffer.data(fk + 249);
    const auto *fk_250 = buffer.data(fk + 250);
    const auto *fk_251 = buffer.data(fk + 251);
    const auto *fk_252 = buffer.data(fk + 252);
    const auto *fk_253 = buffer.data(fk + 253);
    const auto *fk_254 = buffer.data(fk + 254);
    const auto *fk_255 = buffer.data(fk + 255);
    const auto *fk_256 = buffer.data(fk + 256);
    const auto *fk_257 = buffer.data(fk + 257);
    const auto *fk_258 = buffer.data(fk + 258);
    const auto *fk_259 = buffer.data(fk + 259);
    const auto *fk_260 = buffer.data(fk + 260);
    const auto *fk_261 = buffer.data(fk + 261);
    const auto *fk_262 = buffer.data(fk + 262);
    const auto *fk_263 = buffer.data(fk + 263);
    const auto *fk_264 = buffer.data(fk + 264);
    const auto *fk_265 = buffer.data(fk + 265);
    const auto *fk_266 = buffer.data(fk + 266);
    const auto *fk_267 = buffer.data(fk + 267);
    const auto *fk_268 = buffer.data(fk + 268);
    const auto *fk_269 = buffer.data(fk + 269);
    const auto *fk_270 = buffer.data(fk + 270);
    const auto *fk_271 = buffer.data(fk + 271);
    const auto *fk_272 = buffer.data(fk + 272);
    const auto *fk_273 = buffer.data(fk + 273);
    const auto *fk_274 = buffer.data(fk + 274);
    const auto *fk_275 = buffer.data(fk + 275);
    const auto *fk_276 = buffer.data(fk + 276);
    const auto *fk_277 = buffer.data(fk + 277);
    const auto *fk_278 = buffer.data(fk + 278);
    const auto *fk_279 = buffer.data(fk + 279);
    const auto *fk_280 = buffer.data(fk + 280);
    const auto *fk_281 = buffer.data(fk + 281);
    const auto *fk_282 = buffer.data(fk + 282);
    const auto *fk_283 = buffer.data(fk + 283);
    const auto *fk_284 = buffer.data(fk + 284);
    const auto *fk_285 = buffer.data(fk + 285);
    const auto *fk_286 = buffer.data(fk + 286);
    const auto *fk_287 = buffer.data(fk + 287);
    const auto *fk_288 = buffer.data(fk + 288);
    const auto *fk_289 = buffer.data(fk + 289);
    const auto *fk_290 = buffer.data(fk + 290);
    const auto *fk_291 = buffer.data(fk + 291);
    const auto *fk_292 = buffer.data(fk + 292);
    const auto *fk_293 = buffer.data(fk + 293);
    const auto *fk_294 = buffer.data(fk + 294);
    const auto *fk_295 = buffer.data(fk + 295);
    const auto *fk_296 = buffer.data(fk + 296);
    const auto *fk_297 = buffer.data(fk + 297);
    const auto *fk_298 = buffer.data(fk + 298);
    const auto *fk_299 = buffer.data(fk + 299);
    const auto *fk_300 = buffer.data(fk + 300);
    const auto *fk_301 = buffer.data(fk + 301);
    const auto *fk_302 = buffer.data(fk + 302);
    const auto *fk_303 = buffer.data(fk + 303);
    const auto *fk_304 = buffer.data(fk + 304);
    const auto *fk_305 = buffer.data(fk + 305);
    const auto *fk_306 = buffer.data(fk + 306);
    const auto *fk_307 = buffer.data(fk + 307);
    const auto *fk_308 = buffer.data(fk + 308);
    const auto *fk_309 = buffer.data(fk + 309);
    const auto *fk_310 = buffer.data(fk + 310);
    const auto *fk_311 = buffer.data(fk + 311);
    const auto *fk_312 = buffer.data(fk + 312);
    const auto *fk_313 = buffer.data(fk + 313);
    const auto *fk_314 = buffer.data(fk + 314);
    const auto *fk_315 = buffer.data(fk + 315);
    const auto *fk_316 = buffer.data(fk + 316);
    const auto *fk_317 = buffer.data(fk + 317);
    const auto *fk_318 = buffer.data(fk + 318);
    const auto *fk_319 = buffer.data(fk + 319);
    const auto *fk_320 = buffer.data(fk + 320);
    const auto *fk_321 = buffer.data(fk + 321);
    const auto *fk_322 = buffer.data(fk + 322);
    const auto *fk_323 = buffer.data(fk + 323);
    const auto *fk_324 = buffer.data(fk + 324);
    const auto *fk_325 = buffer.data(fk + 325);
    const auto *fk_326 = buffer.data(fk + 326);
    const auto *fk_327 = buffer.data(fk + 327);
    const auto *fk_328 = buffer.data(fk + 328);
    const auto *fk_329 = buffer.data(fk + 329);
    const auto *fk_330 = buffer.data(fk + 330);

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

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, fk_204, fk_205, fk_206, fk_207, \
                         fk_208, hk_528, hk_529, hk_530, hk_531, \
                         hk_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = -3.0 * fk_204[k]
                   + f_0 * hk_528[k];

        t_349[k] = -3.0 * fk_205[k]
                   + f_0 * hk_529[k];

        t_350[k] = -3.0 * fk_206[k]
                   + f_0 * hk_530[k];

        t_351[k] = -3.0 * fk_207[k]
                   + f_0 * hk_531[k];

        t_352[k] = -3.0 * fk_208[k]
                   + f_0 * hk_532[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, fk_209, fk_210, fk_211, fk_212, \
                         fk_213, hk_533, hk_534, hk_535, hk_536, \
                         hk_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = -3.0 * fk_209[k]
                   + f_0 * hk_533[k];

        t_354[k] = -3.0 * fk_210[k]
                   + f_0 * hk_534[k];

        t_355[k] = -3.0 * fk_211[k]
                   + f_0 * hk_535[k];

        t_356[k] = -3.0 * fk_212[k]
                   + f_0 * hk_536[k];

        t_357[k] = -3.0 * fk_213[k]
                   + f_0 * hk_537[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, t_363, t_364, fk_214, fk_215, \
                         hk_538, hk_539, hk_576, hk_577, hk_578, hk_579, \
                         hk_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = -3.0 * fk_214[k]
                   + f_0 * hk_538[k];

        t_359[k] = -3.0 * fk_215[k]
                   + f_0 * hk_539[k];

        t_360[k] = f_0 * hk_576[k];

        t_361[k] = f_0 * hk_577[k];

        t_362[k] = f_0 * hk_578[k];

        t_363[k] = f_0 * hk_579[k];

        t_364[k] = f_0 * hk_580[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, t_370, t_371, t_372, hk_581, \
                         hk_582, hk_583, hk_584, hk_585, hk_586, hk_587, \
                         hk_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_0 * hk_581[k];

        t_366[k] = f_0 * hk_582[k];

        t_367[k] = f_0 * hk_583[k];

        t_368[k] = f_0 * hk_584[k];

        t_369[k] = f_0 * hk_585[k];

        t_370[k] = f_0 * hk_586[k];

        t_371[k] = f_0 * hk_587[k];

        t_372[k] = f_0 * hk_588[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, t_378, t_379, t_380, hk_589, \
                         hk_590, hk_591, hk_592, hk_593, hk_594, hk_595, \
                         hk_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_0 * hk_589[k];

        t_374[k] = f_0 * hk_590[k];

        t_375[k] = f_0 * hk_591[k];

        t_376[k] = f_0 * hk_592[k];

        t_377[k] = f_0 * hk_593[k];

        t_378[k] = f_0 * hk_594[k];

        t_379[k] = f_0 * hk_595[k];

        t_380[k] = f_0 * hk_596[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, t_385, t_386, t_387, t_388, hk_597, \
                         hk_598, hk_599, hk_600, hk_601, hk_602, hk_603, \
                         hk_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_0 * hk_597[k];

        t_382[k] = f_0 * hk_598[k];

        t_383[k] = f_0 * hk_599[k];

        t_384[k] = f_0 * hk_600[k];

        t_385[k] = f_0 * hk_601[k];

        t_386[k] = f_0 * hk_602[k];

        t_387[k] = f_0 * hk_603[k];

        t_388[k] = f_0 * hk_604[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, t_394, t_395, hk_605, hk_606, \
                         hk_607, hk_608, hk_609, hk_610, hk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_0 * hk_605[k];

        t_390[k] = f_0 * hk_606[k];

        t_391[k] = f_0 * hk_607[k];

        t_392[k] = f_0 * hk_608[k];

        t_393[k] = f_0 * hk_609[k];

        t_394[k] = f_0 * hk_610[k];

        t_395[k] = f_0 * hk_611[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, t_400, fk_216, fk_217, fk_218, fk_219, \
                         fk_220, hk_612, hk_613, hk_614, hk_615, \
                         hk_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = -fk_216[k]
                   + f_0 * hk_612[k];

        t_397[k] = -fk_217[k]
                   + f_0 * hk_613[k];

        t_398[k] = -fk_218[k]
                   + f_0 * hk_614[k];

        t_399[k] = -fk_219[k]
                   + f_0 * hk_615[k];

        t_400[k] = -fk_220[k]
                   + f_0 * hk_616[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, t_405, fk_221, fk_222, fk_223, fk_224, \
                         fk_225, hk_617, hk_618, hk_619, hk_620, \
                         hk_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = -fk_221[k]
                   + f_0 * hk_617[k];

        t_402[k] = -fk_222[k]
                   + f_0 * hk_618[k];

        t_403[k] = -fk_223[k]
                   + f_0 * hk_619[k];

        t_404[k] = -fk_224[k]
                   + f_0 * hk_620[k];

        t_405[k] = -fk_225[k]
                   + f_0 * hk_621[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, t_410, fk_226, fk_227, fk_228, fk_229, \
                         fk_230, hk_622, hk_623, hk_624, hk_625, \
                         hk_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = -fk_226[k]
                   + f_0 * hk_622[k];

        t_407[k] = -fk_227[k]
                   + f_0 * hk_623[k];

        t_408[k] = -fk_228[k]
                   + f_0 * hk_624[k];

        t_409[k] = -fk_229[k]
                   + f_0 * hk_625[k];

        t_410[k] = -fk_230[k]
                   + f_0 * hk_626[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, t_415, fk_231, fk_232, fk_233, fk_234, \
                         fk_235, hk_627, hk_628, hk_629, hk_630, \
                         hk_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = -fk_231[k]
                   + f_0 * hk_627[k];

        t_412[k] = -fk_232[k]
                   + f_0 * hk_628[k];

        t_413[k] = -fk_233[k]
                   + f_0 * hk_629[k];

        t_414[k] = -fk_234[k]
                   + f_0 * hk_630[k];

        t_415[k] = -fk_235[k]
                   + f_0 * hk_631[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, t_420, fk_236, fk_237, fk_238, fk_239, \
                         fk_240, hk_632, hk_633, hk_634, hk_635, \
                         hk_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = -fk_236[k]
                   + f_0 * hk_632[k];

        t_417[k] = -fk_237[k]
                   + f_0 * hk_633[k];

        t_418[k] = -fk_238[k]
                   + f_0 * hk_634[k];

        t_419[k] = -fk_239[k]
                   + f_0 * hk_635[k];

        t_420[k] = -fk_240[k]
                   + f_0 * hk_636[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, t_425, fk_241, fk_242, fk_243, fk_244, \
                         fk_245, hk_637, hk_638, hk_639, hk_640, \
                         hk_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = -fk_241[k]
                   + f_0 * hk_637[k];

        t_422[k] = -fk_242[k]
                   + f_0 * hk_638[k];

        t_423[k] = -fk_243[k]
                   + f_0 * hk_639[k];

        t_424[k] = -fk_244[k]
                   + f_0 * hk_640[k];

        t_425[k] = -fk_245[k]
                   + f_0 * hk_641[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, t_430, fk_246, fk_247, fk_248, fk_249, \
                         fk_250, hk_642, hk_643, hk_644, hk_645, \
                         hk_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = -fk_246[k]
                   + f_0 * hk_642[k];

        t_427[k] = -fk_247[k]
                   + f_0 * hk_643[k];

        t_428[k] = -fk_248[k]
                   + f_0 * hk_644[k];

        t_429[k] = -fk_249[k]
                   + f_0 * hk_645[k];

        t_430[k] = -fk_250[k]
                   + f_0 * hk_646[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, fk_251, fk_252, fk_253, fk_254, \
                         fk_255, hk_647, hk_648, hk_649, hk_650, \
                         hk_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = -fk_251[k]
                   + f_0 * hk_647[k];

        t_432[k] = -2.0 * fk_252[k]
                   + f_0 * hk_648[k];

        t_433[k] = -2.0 * fk_253[k]
                   + f_0 * hk_649[k];

        t_434[k] = -2.0 * fk_254[k]
                   + f_0 * hk_650[k];

        t_435[k] = -2.0 * fk_255[k]
                   + f_0 * hk_651[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, fk_256, fk_257, fk_258, fk_259, \
                         fk_260, hk_652, hk_653, hk_654, hk_655, \
                         hk_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = -2.0 * fk_256[k]
                   + f_0 * hk_652[k];

        t_437[k] = -2.0 * fk_257[k]
                   + f_0 * hk_653[k];

        t_438[k] = -2.0 * fk_258[k]
                   + f_0 * hk_654[k];

        t_439[k] = -2.0 * fk_259[k]
                   + f_0 * hk_655[k];

        t_440[k] = -2.0 * fk_260[k]
                   + f_0 * hk_656[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, t_445, fk_261, fk_262, fk_263, fk_264, \
                         fk_265, hk_657, hk_658, hk_659, hk_660, \
                         hk_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = -2.0 * fk_261[k]
                   + f_0 * hk_657[k];

        t_442[k] = -2.0 * fk_262[k]
                   + f_0 * hk_658[k];

        t_443[k] = -2.0 * fk_263[k]
                   + f_0 * hk_659[k];

        t_444[k] = -2.0 * fk_264[k]
                   + f_0 * hk_660[k];

        t_445[k] = -2.0 * fk_265[k]
                   + f_0 * hk_661[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, t_450, fk_266, fk_267, fk_268, fk_269, \
                         fk_270, hk_662, hk_663, hk_664, hk_665, \
                         hk_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = -2.0 * fk_266[k]
                   + f_0 * hk_662[k];

        t_447[k] = -2.0 * fk_267[k]
                   + f_0 * hk_663[k];

        t_448[k] = -2.0 * fk_268[k]
                   + f_0 * hk_664[k];

        t_449[k] = -2.0 * fk_269[k]
                   + f_0 * hk_665[k];

        t_450[k] = -2.0 * fk_270[k]
                   + f_0 * hk_666[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, t_455, fk_271, fk_272, fk_273, fk_274, \
                         fk_275, hk_667, hk_668, hk_669, hk_670, \
                         hk_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = -2.0 * fk_271[k]
                   + f_0 * hk_667[k];

        t_452[k] = -2.0 * fk_272[k]
                   + f_0 * hk_668[k];

        t_453[k] = -2.0 * fk_273[k]
                   + f_0 * hk_669[k];

        t_454[k] = -2.0 * fk_274[k]
                   + f_0 * hk_670[k];

        t_455[k] = -2.0 * fk_275[k]
                   + f_0 * hk_671[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, t_460, fk_276, fk_277, fk_278, fk_279, \
                         fk_280, hk_672, hk_673, hk_674, hk_675, \
                         hk_676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = -2.0 * fk_276[k]
                   + f_0 * hk_672[k];

        t_457[k] = -2.0 * fk_277[k]
                   + f_0 * hk_673[k];

        t_458[k] = -2.0 * fk_278[k]
                   + f_0 * hk_674[k];

        t_459[k] = -2.0 * fk_279[k]
                   + f_0 * hk_675[k];

        t_460[k] = -2.0 * fk_280[k]
                   + f_0 * hk_676[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, t_465, fk_281, fk_282, fk_283, fk_284, \
                         fk_285, hk_677, hk_678, hk_679, hk_680, \
                         hk_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = -2.0 * fk_281[k]
                   + f_0 * hk_677[k];

        t_462[k] = -2.0 * fk_282[k]
                   + f_0 * hk_678[k];

        t_463[k] = -2.0 * fk_283[k]
                   + f_0 * hk_679[k];

        t_464[k] = -2.0 * fk_284[k]
                   + f_0 * hk_680[k];

        t_465[k] = -2.0 * fk_285[k]
                   + f_0 * hk_681[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, fk_286, fk_287, fk_288, fk_289, \
                         fk_290, hk_682, hk_683, hk_684, hk_685, \
                         hk_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = -2.0 * fk_286[k]
                   + f_0 * hk_682[k];

        t_467[k] = -2.0 * fk_287[k]
                   + f_0 * hk_683[k];

        t_468[k] = -3.0 * fk_288[k]
                   + f_0 * hk_684[k];

        t_469[k] = -3.0 * fk_289[k]
                   + f_0 * hk_685[k];

        t_470[k] = -3.0 * fk_290[k]
                   + f_0 * hk_686[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, t_475, fk_291, fk_292, fk_293, fk_294, \
                         fk_295, hk_687, hk_688, hk_689, hk_690, \
                         hk_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = -3.0 * fk_291[k]
                   + f_0 * hk_687[k];

        t_472[k] = -3.0 * fk_292[k]
                   + f_0 * hk_688[k];

        t_473[k] = -3.0 * fk_293[k]
                   + f_0 * hk_689[k];

        t_474[k] = -3.0 * fk_294[k]
                   + f_0 * hk_690[k];

        t_475[k] = -3.0 * fk_295[k]
                   + f_0 * hk_691[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, fk_296, fk_297, fk_298, fk_299, \
                         fk_300, hk_692, hk_693, hk_694, hk_695, \
                         hk_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = -3.0 * fk_296[k]
                   + f_0 * hk_692[k];

        t_477[k] = -3.0 * fk_297[k]
                   + f_0 * hk_693[k];

        t_478[k] = -3.0 * fk_298[k]
                   + f_0 * hk_694[k];

        t_479[k] = -3.0 * fk_299[k]
                   + f_0 * hk_695[k];

        t_480[k] = -3.0 * fk_300[k]
                   + f_0 * hk_696[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, fk_301, fk_302, fk_303, fk_304, \
                         fk_305, hk_697, hk_698, hk_699, hk_700, \
                         hk_701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = -3.0 * fk_301[k]
                   + f_0 * hk_697[k];

        t_482[k] = -3.0 * fk_302[k]
                   + f_0 * hk_698[k];

        t_483[k] = -3.0 * fk_303[k]
                   + f_0 * hk_699[k];

        t_484[k] = -3.0 * fk_304[k]
                   + f_0 * hk_700[k];

        t_485[k] = -3.0 * fk_305[k]
                   + f_0 * hk_701[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, fk_306, fk_307, fk_308, fk_309, \
                         fk_310, hk_702, hk_703, hk_704, hk_705, \
                         hk_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = -3.0 * fk_306[k]
                   + f_0 * hk_702[k];

        t_487[k] = -3.0 * fk_307[k]
                   + f_0 * hk_703[k];

        t_488[k] = -3.0 * fk_308[k]
                   + f_0 * hk_704[k];

        t_489[k] = -3.0 * fk_309[k]
                   + f_0 * hk_705[k];

        t_490[k] = -3.0 * fk_310[k]
                   + f_0 * hk_706[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, fk_311, fk_312, fk_313, fk_314, \
                         fk_315, hk_707, hk_708, hk_709, hk_710, \
                         hk_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = -3.0 * fk_311[k]
                   + f_0 * hk_707[k];

        t_492[k] = -3.0 * fk_312[k]
                   + f_0 * hk_708[k];

        t_493[k] = -3.0 * fk_313[k]
                   + f_0 * hk_709[k];

        t_494[k] = -3.0 * fk_314[k]
                   + f_0 * hk_710[k];

        t_495[k] = -3.0 * fk_315[k]
                   + f_0 * hk_711[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, t_500, fk_316, fk_317, fk_318, fk_319, \
                         fk_320, hk_712, hk_713, hk_714, hk_715, \
                         hk_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = -3.0 * fk_316[k]
                   + f_0 * hk_712[k];

        t_497[k] = -3.0 * fk_317[k]
                   + f_0 * hk_713[k];

        t_498[k] = -3.0 * fk_318[k]
                   + f_0 * hk_714[k];

        t_499[k] = -3.0 * fk_319[k]
                   + f_0 * hk_715[k];

        t_500[k] = -3.0 * fk_320[k]
                   + f_0 * hk_716[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, t_505, fk_321, fk_322, fk_323, fk_324, \
                         fk_325, hk_717, hk_718, hk_719, hk_720, \
                         hk_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = -3.0 * fk_321[k]
                   + f_0 * hk_717[k];

        t_502[k] = -3.0 * fk_322[k]
                   + f_0 * hk_718[k];

        t_503[k] = -3.0 * fk_323[k]
                   + f_0 * hk_719[k];

        t_504[k] = -4.0 * fk_324[k]
                   + f_0 * hk_720[k];

        t_505[k] = -4.0 * fk_325[k]
                   + f_0 * hk_721[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, t_509, t_510, fk_326, fk_327, fk_328, fk_329, \
                         fk_330, hk_722, hk_723, hk_724, hk_725, \
                         hk_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = -4.0 * fk_326[k]
                   + f_0 * hk_722[k];

        t_507[k] = -4.0 * fk_327[k]
                   + f_0 * hk_723[k];

        t_508[k] = -4.0 * fk_328[k]
                   + f_0 * hk_724[k];

        t_509[k] = -4.0 * fk_329[k]
                   + f_0 * hk_725[k];

        t_510[k] = -4.0 * fk_330[k]
                   + f_0 * hk_726[k];
    }
}

static auto
compute_prim_geom_10_gk_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fk, const size_t hk,
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

    const auto *fk_331 = buffer.data(fk + 331);
    const auto *fk_332 = buffer.data(fk + 332);
    const auto *fk_333 = buffer.data(fk + 333);
    const auto *fk_334 = buffer.data(fk + 334);
    const auto *fk_335 = buffer.data(fk + 335);
    const auto *fk_336 = buffer.data(fk + 336);
    const auto *fk_337 = buffer.data(fk + 337);
    const auto *fk_338 = buffer.data(fk + 338);
    const auto *fk_339 = buffer.data(fk + 339);
    const auto *fk_340 = buffer.data(fk + 340);
    const auto *fk_341 = buffer.data(fk + 341);
    const auto *fk_342 = buffer.data(fk + 342);
    const auto *fk_343 = buffer.data(fk + 343);
    const auto *fk_344 = buffer.data(fk + 344);
    const auto *fk_345 = buffer.data(fk + 345);
    const auto *fk_346 = buffer.data(fk + 346);
    const auto *fk_347 = buffer.data(fk + 347);
    const auto *fk_348 = buffer.data(fk + 348);
    const auto *fk_349 = buffer.data(fk + 349);
    const auto *fk_350 = buffer.data(fk + 350);
    const auto *fk_351 = buffer.data(fk + 351);
    const auto *fk_352 = buffer.data(fk + 352);
    const auto *fk_353 = buffer.data(fk + 353);
    const auto *fk_354 = buffer.data(fk + 354);
    const auto *fk_355 = buffer.data(fk + 355);
    const auto *fk_356 = buffer.data(fk + 356);
    const auto *fk_357 = buffer.data(fk + 357);
    const auto *fk_358 = buffer.data(fk + 358);
    const auto *fk_359 = buffer.data(fk + 359);

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
    const auto *hk_754 = buffer.data(hk + 754);
    const auto *hk_755 = buffer.data(hk + 755);

#pragma omp simd aligned(t_511, t_512, t_513, t_514, t_515, fk_331, fk_332, fk_333, fk_334, \
                         fk_335, hk_727, hk_728, hk_729, hk_730, \
                         hk_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = -4.0 * fk_331[k]
                   + f_0 * hk_727[k];

        t_512[k] = -4.0 * fk_332[k]
                   + f_0 * hk_728[k];

        t_513[k] = -4.0 * fk_333[k]
                   + f_0 * hk_729[k];

        t_514[k] = -4.0 * fk_334[k]
                   + f_0 * hk_730[k];

        t_515[k] = -4.0 * fk_335[k]
                   + f_0 * hk_731[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, t_520, fk_336, fk_337, fk_338, fk_339, \
                         fk_340, hk_732, hk_733, hk_734, hk_735, \
                         hk_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = -4.0 * fk_336[k]
                   + f_0 * hk_732[k];

        t_517[k] = -4.0 * fk_337[k]
                   + f_0 * hk_733[k];

        t_518[k] = -4.0 * fk_338[k]
                   + f_0 * hk_734[k];

        t_519[k] = -4.0 * fk_339[k]
                   + f_0 * hk_735[k];

        t_520[k] = -4.0 * fk_340[k]
                   + f_0 * hk_736[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, t_525, fk_341, fk_342, fk_343, fk_344, \
                         fk_345, hk_737, hk_738, hk_739, hk_740, \
                         hk_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = -4.0 * fk_341[k]
                   + f_0 * hk_737[k];

        t_522[k] = -4.0 * fk_342[k]
                   + f_0 * hk_738[k];

        t_523[k] = -4.0 * fk_343[k]
                   + f_0 * hk_739[k];

        t_524[k] = -4.0 * fk_344[k]
                   + f_0 * hk_740[k];

        t_525[k] = -4.0 * fk_345[k]
                   + f_0 * hk_741[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, t_530, fk_346, fk_347, fk_348, fk_349, \
                         fk_350, hk_742, hk_743, hk_744, hk_745, \
                         hk_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = -4.0 * fk_346[k]
                   + f_0 * hk_742[k];

        t_527[k] = -4.0 * fk_347[k]
                   + f_0 * hk_743[k];

        t_528[k] = -4.0 * fk_348[k]
                   + f_0 * hk_744[k];

        t_529[k] = -4.0 * fk_349[k]
                   + f_0 * hk_745[k];

        t_530[k] = -4.0 * fk_350[k]
                   + f_0 * hk_746[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, t_535, fk_351, fk_352, fk_353, fk_354, \
                         fk_355, hk_747, hk_748, hk_749, hk_750, \
                         hk_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = -4.0 * fk_351[k]
                   + f_0 * hk_747[k];

        t_532[k] = -4.0 * fk_352[k]
                   + f_0 * hk_748[k];

        t_533[k] = -4.0 * fk_353[k]
                   + f_0 * hk_749[k];

        t_534[k] = -4.0 * fk_354[k]
                   + f_0 * hk_750[k];

        t_535[k] = -4.0 * fk_355[k]
                   + f_0 * hk_751[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, fk_356, fk_357, fk_358, fk_359, hk_752, \
                         hk_753, hk_754, hk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = -4.0 * fk_356[k]
                   + f_0 * hk_752[k];

        t_537[k] = -4.0 * fk_357[k]
                   + f_0 * hk_753[k];

        t_538[k] = -4.0 * fk_358[k]
                   + f_0 * hk_754[k];

        t_539[k] = -4.0 * fk_359[k]
                   + f_0 * hk_755[k];
    }
}

auto
compute_prim_geom_10_gk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t fk, const size_t hk,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gk_electron_repulsion_2_piece0(buffer, target, fk, hk, ncols, alpha);

    compute_prim_geom_10_gk_electron_repulsion_2_piece1(buffer, target, fk, hk, ncols, alpha);

    compute_prim_geom_10_gk_electron_repulsion_2_piece2(buffer, target, fk, hk, ncols, alpha);

    compute_prim_geom_10_gk_electron_repulsion_2_piece3(buffer, target, fk, hk, ncols, alpha);
}

}  // namespace simdt2ceri
