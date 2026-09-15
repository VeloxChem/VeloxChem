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


#include "SimdElectronRepulsionGeom10VrrRecII.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_ii_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
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

    const auto *hi_0 = buffer.data(hi + 0);
    const auto *hi_1 = buffer.data(hi + 1);
    const auto *hi_2 = buffer.data(hi + 2);
    const auto *hi_3 = buffer.data(hi + 3);
    const auto *hi_4 = buffer.data(hi + 4);
    const auto *hi_5 = buffer.data(hi + 5);
    const auto *hi_6 = buffer.data(hi + 6);
    const auto *hi_7 = buffer.data(hi + 7);
    const auto *hi_8 = buffer.data(hi + 8);
    const auto *hi_9 = buffer.data(hi + 9);
    const auto *hi_10 = buffer.data(hi + 10);
    const auto *hi_11 = buffer.data(hi + 11);
    const auto *hi_12 = buffer.data(hi + 12);
    const auto *hi_13 = buffer.data(hi + 13);
    const auto *hi_14 = buffer.data(hi + 14);
    const auto *hi_15 = buffer.data(hi + 15);
    const auto *hi_16 = buffer.data(hi + 16);
    const auto *hi_17 = buffer.data(hi + 17);
    const auto *hi_18 = buffer.data(hi + 18);
    const auto *hi_19 = buffer.data(hi + 19);
    const auto *hi_20 = buffer.data(hi + 20);
    const auto *hi_21 = buffer.data(hi + 21);
    const auto *hi_22 = buffer.data(hi + 22);
    const auto *hi_23 = buffer.data(hi + 23);
    const auto *hi_24 = buffer.data(hi + 24);
    const auto *hi_25 = buffer.data(hi + 25);
    const auto *hi_26 = buffer.data(hi + 26);
    const auto *hi_27 = buffer.data(hi + 27);
    const auto *hi_28 = buffer.data(hi + 28);
    const auto *hi_29 = buffer.data(hi + 29);
    const auto *hi_30 = buffer.data(hi + 30);
    const auto *hi_31 = buffer.data(hi + 31);
    const auto *hi_32 = buffer.data(hi + 32);
    const auto *hi_33 = buffer.data(hi + 33);
    const auto *hi_34 = buffer.data(hi + 34);
    const auto *hi_35 = buffer.data(hi + 35);
    const auto *hi_36 = buffer.data(hi + 36);
    const auto *hi_37 = buffer.data(hi + 37);
    const auto *hi_38 = buffer.data(hi + 38);
    const auto *hi_39 = buffer.data(hi + 39);
    const auto *hi_40 = buffer.data(hi + 40);
    const auto *hi_41 = buffer.data(hi + 41);
    const auto *hi_42 = buffer.data(hi + 42);
    const auto *hi_43 = buffer.data(hi + 43);
    const auto *hi_44 = buffer.data(hi + 44);
    const auto *hi_45 = buffer.data(hi + 45);
    const auto *hi_46 = buffer.data(hi + 46);
    const auto *hi_47 = buffer.data(hi + 47);
    const auto *hi_48 = buffer.data(hi + 48);
    const auto *hi_49 = buffer.data(hi + 49);
    const auto *hi_50 = buffer.data(hi + 50);
    const auto *hi_51 = buffer.data(hi + 51);
    const auto *hi_52 = buffer.data(hi + 52);
    const auto *hi_53 = buffer.data(hi + 53);
    const auto *hi_54 = buffer.data(hi + 54);
    const auto *hi_55 = buffer.data(hi + 55);
    const auto *hi_56 = buffer.data(hi + 56);
    const auto *hi_57 = buffer.data(hi + 57);
    const auto *hi_58 = buffer.data(hi + 58);
    const auto *hi_59 = buffer.data(hi + 59);
    const auto *hi_60 = buffer.data(hi + 60);
    const auto *hi_61 = buffer.data(hi + 61);
    const auto *hi_62 = buffer.data(hi + 62);
    const auto *hi_63 = buffer.data(hi + 63);
    const auto *hi_64 = buffer.data(hi + 64);
    const auto *hi_65 = buffer.data(hi + 65);
    const auto *hi_66 = buffer.data(hi + 66);
    const auto *hi_67 = buffer.data(hi + 67);
    const auto *hi_68 = buffer.data(hi + 68);
    const auto *hi_69 = buffer.data(hi + 69);
    const auto *hi_70 = buffer.data(hi + 70);
    const auto *hi_71 = buffer.data(hi + 71);
    const auto *hi_72 = buffer.data(hi + 72);
    const auto *hi_73 = buffer.data(hi + 73);
    const auto *hi_74 = buffer.data(hi + 74);
    const auto *hi_75 = buffer.data(hi + 75);
    const auto *hi_76 = buffer.data(hi + 76);
    const auto *hi_77 = buffer.data(hi + 77);
    const auto *hi_78 = buffer.data(hi + 78);
    const auto *hi_79 = buffer.data(hi + 79);
    const auto *hi_80 = buffer.data(hi + 80);
    const auto *hi_81 = buffer.data(hi + 81);
    const auto *hi_82 = buffer.data(hi + 82);
    const auto *hi_83 = buffer.data(hi + 83);
    const auto *hi_84 = buffer.data(hi + 84);
    const auto *hi_85 = buffer.data(hi + 85);
    const auto *hi_86 = buffer.data(hi + 86);
    const auto *hi_87 = buffer.data(hi + 87);
    const auto *hi_88 = buffer.data(hi + 88);
    const auto *hi_89 = buffer.data(hi + 89);
    const auto *hi_90 = buffer.data(hi + 90);
    const auto *hi_91 = buffer.data(hi + 91);
    const auto *hi_92 = buffer.data(hi + 92);
    const auto *hi_93 = buffer.data(hi + 93);
    const auto *hi_94 = buffer.data(hi + 94);
    const auto *hi_95 = buffer.data(hi + 95);
    const auto *hi_96 = buffer.data(hi + 96);
    const auto *hi_97 = buffer.data(hi + 97);
    const auto *hi_98 = buffer.data(hi + 98);
    const auto *hi_99 = buffer.data(hi + 99);
    const auto *hi_100 = buffer.data(hi + 100);
    const auto *hi_101 = buffer.data(hi + 101);
    const auto *hi_102 = buffer.data(hi + 102);
    const auto *hi_103 = buffer.data(hi + 103);
    const auto *hi_104 = buffer.data(hi + 104);
    const auto *hi_105 = buffer.data(hi + 105);
    const auto *hi_106 = buffer.data(hi + 106);
    const auto *hi_107 = buffer.data(hi + 107);
    const auto *hi_108 = buffer.data(hi + 108);
    const auto *hi_109 = buffer.data(hi + 109);
    const auto *hi_110 = buffer.data(hi + 110);
    const auto *hi_111 = buffer.data(hi + 111);
    const auto *hi_112 = buffer.data(hi + 112);
    const auto *hi_113 = buffer.data(hi + 113);
    const auto *hi_114 = buffer.data(hi + 114);
    const auto *hi_115 = buffer.data(hi + 115);
    const auto *hi_116 = buffer.data(hi + 116);
    const auto *hi_117 = buffer.data(hi + 117);
    const auto *hi_118 = buffer.data(hi + 118);
    const auto *hi_119 = buffer.data(hi + 119);
    const auto *hi_120 = buffer.data(hi + 120);
    const auto *hi_121 = buffer.data(hi + 121);
    const auto *hi_122 = buffer.data(hi + 122);
    const auto *hi_123 = buffer.data(hi + 123);
    const auto *hi_124 = buffer.data(hi + 124);
    const auto *hi_125 = buffer.data(hi + 125);
    const auto *hi_126 = buffer.data(hi + 126);
    const auto *hi_127 = buffer.data(hi + 127);
    const auto *hi_128 = buffer.data(hi + 128);
    const auto *hi_129 = buffer.data(hi + 129);
    const auto *hi_130 = buffer.data(hi + 130);
    const auto *hi_131 = buffer.data(hi + 131);
    const auto *hi_132 = buffer.data(hi + 132);
    const auto *hi_133 = buffer.data(hi + 133);
    const auto *hi_134 = buffer.data(hi + 134);
    const auto *hi_135 = buffer.data(hi + 135);
    const auto *hi_136 = buffer.data(hi + 136);
    const auto *hi_137 = buffer.data(hi + 137);
    const auto *hi_138 = buffer.data(hi + 138);
    const auto *hi_139 = buffer.data(hi + 139);
    const auto *hi_140 = buffer.data(hi + 140);
    const auto *hi_141 = buffer.data(hi + 141);
    const auto *hi_142 = buffer.data(hi + 142);
    const auto *hi_143 = buffer.data(hi + 143);
    const auto *hi_144 = buffer.data(hi + 144);
    const auto *hi_145 = buffer.data(hi + 145);
    const auto *hi_146 = buffer.data(hi + 146);
    const auto *hi_147 = buffer.data(hi + 147);
    const auto *hi_148 = buffer.data(hi + 148);
    const auto *hi_149 = buffer.data(hi + 149);

    const auto *ki_0 = buffer.data(ki + 0);
    const auto *ki_1 = buffer.data(ki + 1);
    const auto *ki_2 = buffer.data(ki + 2);
    const auto *ki_3 = buffer.data(ki + 3);
    const auto *ki_4 = buffer.data(ki + 4);
    const auto *ki_5 = buffer.data(ki + 5);
    const auto *ki_6 = buffer.data(ki + 6);
    const auto *ki_7 = buffer.data(ki + 7);
    const auto *ki_8 = buffer.data(ki + 8);
    const auto *ki_9 = buffer.data(ki + 9);
    const auto *ki_10 = buffer.data(ki + 10);
    const auto *ki_11 = buffer.data(ki + 11);
    const auto *ki_12 = buffer.data(ki + 12);
    const auto *ki_13 = buffer.data(ki + 13);
    const auto *ki_14 = buffer.data(ki + 14);
    const auto *ki_15 = buffer.data(ki + 15);
    const auto *ki_16 = buffer.data(ki + 16);
    const auto *ki_17 = buffer.data(ki + 17);
    const auto *ki_18 = buffer.data(ki + 18);
    const auto *ki_19 = buffer.data(ki + 19);
    const auto *ki_20 = buffer.data(ki + 20);
    const auto *ki_21 = buffer.data(ki + 21);
    const auto *ki_22 = buffer.data(ki + 22);
    const auto *ki_23 = buffer.data(ki + 23);
    const auto *ki_24 = buffer.data(ki + 24);
    const auto *ki_25 = buffer.data(ki + 25);
    const auto *ki_26 = buffer.data(ki + 26);
    const auto *ki_27 = buffer.data(ki + 27);
    const auto *ki_28 = buffer.data(ki + 28);
    const auto *ki_29 = buffer.data(ki + 29);
    const auto *ki_30 = buffer.data(ki + 30);
    const auto *ki_31 = buffer.data(ki + 31);
    const auto *ki_32 = buffer.data(ki + 32);
    const auto *ki_33 = buffer.data(ki + 33);
    const auto *ki_34 = buffer.data(ki + 34);
    const auto *ki_35 = buffer.data(ki + 35);
    const auto *ki_36 = buffer.data(ki + 36);
    const auto *ki_37 = buffer.data(ki + 37);
    const auto *ki_38 = buffer.data(ki + 38);
    const auto *ki_39 = buffer.data(ki + 39);
    const auto *ki_40 = buffer.data(ki + 40);
    const auto *ki_41 = buffer.data(ki + 41);
    const auto *ki_42 = buffer.data(ki + 42);
    const auto *ki_43 = buffer.data(ki + 43);
    const auto *ki_44 = buffer.data(ki + 44);
    const auto *ki_45 = buffer.data(ki + 45);
    const auto *ki_46 = buffer.data(ki + 46);
    const auto *ki_47 = buffer.data(ki + 47);
    const auto *ki_48 = buffer.data(ki + 48);
    const auto *ki_49 = buffer.data(ki + 49);
    const auto *ki_50 = buffer.data(ki + 50);
    const auto *ki_51 = buffer.data(ki + 51);
    const auto *ki_52 = buffer.data(ki + 52);
    const auto *ki_53 = buffer.data(ki + 53);
    const auto *ki_54 = buffer.data(ki + 54);
    const auto *ki_55 = buffer.data(ki + 55);
    const auto *ki_56 = buffer.data(ki + 56);
    const auto *ki_57 = buffer.data(ki + 57);
    const auto *ki_58 = buffer.data(ki + 58);
    const auto *ki_59 = buffer.data(ki + 59);
    const auto *ki_60 = buffer.data(ki + 60);
    const auto *ki_61 = buffer.data(ki + 61);
    const auto *ki_62 = buffer.data(ki + 62);
    const auto *ki_63 = buffer.data(ki + 63);
    const auto *ki_64 = buffer.data(ki + 64);
    const auto *ki_65 = buffer.data(ki + 65);
    const auto *ki_66 = buffer.data(ki + 66);
    const auto *ki_67 = buffer.data(ki + 67);
    const auto *ki_68 = buffer.data(ki + 68);
    const auto *ki_69 = buffer.data(ki + 69);
    const auto *ki_70 = buffer.data(ki + 70);
    const auto *ki_71 = buffer.data(ki + 71);
    const auto *ki_72 = buffer.data(ki + 72);
    const auto *ki_73 = buffer.data(ki + 73);
    const auto *ki_74 = buffer.data(ki + 74);
    const auto *ki_75 = buffer.data(ki + 75);
    const auto *ki_76 = buffer.data(ki + 76);
    const auto *ki_77 = buffer.data(ki + 77);
    const auto *ki_78 = buffer.data(ki + 78);
    const auto *ki_79 = buffer.data(ki + 79);
    const auto *ki_80 = buffer.data(ki + 80);
    const auto *ki_81 = buffer.data(ki + 81);
    const auto *ki_82 = buffer.data(ki + 82);
    const auto *ki_83 = buffer.data(ki + 83);
    const auto *ki_84 = buffer.data(ki + 84);
    const auto *ki_85 = buffer.data(ki + 85);
    const auto *ki_86 = buffer.data(ki + 86);
    const auto *ki_87 = buffer.data(ki + 87);
    const auto *ki_88 = buffer.data(ki + 88);
    const auto *ki_89 = buffer.data(ki + 89);
    const auto *ki_90 = buffer.data(ki + 90);
    const auto *ki_91 = buffer.data(ki + 91);
    const auto *ki_92 = buffer.data(ki + 92);
    const auto *ki_93 = buffer.data(ki + 93);
    const auto *ki_94 = buffer.data(ki + 94);
    const auto *ki_95 = buffer.data(ki + 95);
    const auto *ki_96 = buffer.data(ki + 96);
    const auto *ki_97 = buffer.data(ki + 97);
    const auto *ki_98 = buffer.data(ki + 98);
    const auto *ki_99 = buffer.data(ki + 99);
    const auto *ki_100 = buffer.data(ki + 100);
    const auto *ki_101 = buffer.data(ki + 101);
    const auto *ki_102 = buffer.data(ki + 102);
    const auto *ki_103 = buffer.data(ki + 103);
    const auto *ki_104 = buffer.data(ki + 104);
    const auto *ki_105 = buffer.data(ki + 105);
    const auto *ki_106 = buffer.data(ki + 106);
    const auto *ki_107 = buffer.data(ki + 107);
    const auto *ki_108 = buffer.data(ki + 108);
    const auto *ki_109 = buffer.data(ki + 109);
    const auto *ki_110 = buffer.data(ki + 110);
    const auto *ki_111 = buffer.data(ki + 111);
    const auto *ki_112 = buffer.data(ki + 112);
    const auto *ki_113 = buffer.data(ki + 113);
    const auto *ki_114 = buffer.data(ki + 114);
    const auto *ki_115 = buffer.data(ki + 115);
    const auto *ki_116 = buffer.data(ki + 116);
    const auto *ki_117 = buffer.data(ki + 117);
    const auto *ki_118 = buffer.data(ki + 118);
    const auto *ki_119 = buffer.data(ki + 119);
    const auto *ki_120 = buffer.data(ki + 120);
    const auto *ki_121 = buffer.data(ki + 121);
    const auto *ki_122 = buffer.data(ki + 122);
    const auto *ki_123 = buffer.data(ki + 123);
    const auto *ki_124 = buffer.data(ki + 124);
    const auto *ki_125 = buffer.data(ki + 125);
    const auto *ki_126 = buffer.data(ki + 126);
    const auto *ki_127 = buffer.data(ki + 127);
    const auto *ki_128 = buffer.data(ki + 128);
    const auto *ki_129 = buffer.data(ki + 129);
    const auto *ki_130 = buffer.data(ki + 130);
    const auto *ki_131 = buffer.data(ki + 131);
    const auto *ki_132 = buffer.data(ki + 132);
    const auto *ki_133 = buffer.data(ki + 133);
    const auto *ki_134 = buffer.data(ki + 134);
    const auto *ki_135 = buffer.data(ki + 135);
    const auto *ki_136 = buffer.data(ki + 136);
    const auto *ki_137 = buffer.data(ki + 137);
    const auto *ki_138 = buffer.data(ki + 138);
    const auto *ki_139 = buffer.data(ki + 139);
    const auto *ki_140 = buffer.data(ki + 140);
    const auto *ki_141 = buffer.data(ki + 141);
    const auto *ki_142 = buffer.data(ki + 142);
    const auto *ki_143 = buffer.data(ki + 143);
    const auto *ki_144 = buffer.data(ki + 144);
    const auto *ki_145 = buffer.data(ki + 145);
    const auto *ki_146 = buffer.data(ki + 146);
    const auto *ki_147 = buffer.data(ki + 147);
    const auto *ki_148 = buffer.data(ki + 148);
    const auto *ki_149 = buffer.data(ki + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, hi_0, hi_1, hi_2, hi_3, hi_4, ki_0, ki_1, \
                         ki_2, ki_3, ki_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -6.0 * hi_0[k]
                 + f_0 * ki_0[k];

        t_1[k] = -6.0 * hi_1[k]
                 + f_0 * ki_1[k];

        t_2[k] = -6.0 * hi_2[k]
                 + f_0 * ki_2[k];

        t_3[k] = -6.0 * hi_3[k]
                 + f_0 * ki_3[k];

        t_4[k] = -6.0 * hi_4[k]
                 + f_0 * ki_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, hi_5, hi_6, hi_7, hi_8, hi_9, ki_5, ki_6, \
                         ki_7, ki_8, ki_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -6.0 * hi_5[k]
                 + f_0 * ki_5[k];

        t_6[k] = -6.0 * hi_6[k]
                 + f_0 * ki_6[k];

        t_7[k] = -6.0 * hi_7[k]
                 + f_0 * ki_7[k];

        t_8[k] = -6.0 * hi_8[k]
                 + f_0 * ki_8[k];

        t_9[k] = -6.0 * hi_9[k]
                 + f_0 * ki_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, hi_10, hi_11, hi_12, hi_13, hi_14, \
                         ki_10, ki_11, ki_12, ki_13, ki_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -6.0 * hi_10[k]
                  + f_0 * ki_10[k];

        t_11[k] = -6.0 * hi_11[k]
                  + f_0 * ki_11[k];

        t_12[k] = -6.0 * hi_12[k]
                  + f_0 * ki_12[k];

        t_13[k] = -6.0 * hi_13[k]
                  + f_0 * ki_13[k];

        t_14[k] = -6.0 * hi_14[k]
                  + f_0 * ki_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, hi_15, hi_16, hi_17, hi_18, hi_19, \
                         ki_15, ki_16, ki_17, ki_18, ki_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -6.0 * hi_15[k]
                  + f_0 * ki_15[k];

        t_16[k] = -6.0 * hi_16[k]
                  + f_0 * ki_16[k];

        t_17[k] = -6.0 * hi_17[k]
                  + f_0 * ki_17[k];

        t_18[k] = -6.0 * hi_18[k]
                  + f_0 * ki_18[k];

        t_19[k] = -6.0 * hi_19[k]
                  + f_0 * ki_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, hi_20, hi_21, hi_22, hi_23, hi_24, \
                         ki_20, ki_21, ki_22, ki_23, ki_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -6.0 * hi_20[k]
                  + f_0 * ki_20[k];

        t_21[k] = -6.0 * hi_21[k]
                  + f_0 * ki_21[k];

        t_22[k] = -6.0 * hi_22[k]
                  + f_0 * ki_22[k];

        t_23[k] = -6.0 * hi_23[k]
                  + f_0 * ki_23[k];

        t_24[k] = -6.0 * hi_24[k]
                  + f_0 * ki_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, hi_25, hi_26, hi_27, hi_28, hi_29, \
                         ki_25, ki_26, ki_27, ki_28, ki_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -6.0 * hi_25[k]
                  + f_0 * ki_25[k];

        t_26[k] = -6.0 * hi_26[k]
                  + f_0 * ki_26[k];

        t_27[k] = -6.0 * hi_27[k]
                  + f_0 * ki_27[k];

        t_28[k] = -5.0 * hi_28[k]
                  + f_0 * ki_28[k];

        t_29[k] = -5.0 * hi_29[k]
                  + f_0 * ki_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, hi_30, hi_31, hi_32, hi_33, hi_34, \
                         ki_30, ki_31, ki_32, ki_33, ki_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -5.0 * hi_30[k]
                  + f_0 * ki_30[k];

        t_31[k] = -5.0 * hi_31[k]
                  + f_0 * ki_31[k];

        t_32[k] = -5.0 * hi_32[k]
                  + f_0 * ki_32[k];

        t_33[k] = -5.0 * hi_33[k]
                  + f_0 * ki_33[k];

        t_34[k] = -5.0 * hi_34[k]
                  + f_0 * ki_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, hi_35, hi_36, hi_37, hi_38, hi_39, \
                         ki_35, ki_36, ki_37, ki_38, ki_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -5.0 * hi_35[k]
                  + f_0 * ki_35[k];

        t_36[k] = -5.0 * hi_36[k]
                  + f_0 * ki_36[k];

        t_37[k] = -5.0 * hi_37[k]
                  + f_0 * ki_37[k];

        t_38[k] = -5.0 * hi_38[k]
                  + f_0 * ki_38[k];

        t_39[k] = -5.0 * hi_39[k]
                  + f_0 * ki_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, hi_40, hi_41, hi_42, hi_43, hi_44, \
                         ki_40, ki_41, ki_42, ki_43, ki_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -5.0 * hi_40[k]
                  + f_0 * ki_40[k];

        t_41[k] = -5.0 * hi_41[k]
                  + f_0 * ki_41[k];

        t_42[k] = -5.0 * hi_42[k]
                  + f_0 * ki_42[k];

        t_43[k] = -5.0 * hi_43[k]
                  + f_0 * ki_43[k];

        t_44[k] = -5.0 * hi_44[k]
                  + f_0 * ki_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hi_45, hi_46, hi_47, hi_48, hi_49, \
                         ki_45, ki_46, ki_47, ki_48, ki_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -5.0 * hi_45[k]
                  + f_0 * ki_45[k];

        t_46[k] = -5.0 * hi_46[k]
                  + f_0 * ki_46[k];

        t_47[k] = -5.0 * hi_47[k]
                  + f_0 * ki_47[k];

        t_48[k] = -5.0 * hi_48[k]
                  + f_0 * ki_48[k];

        t_49[k] = -5.0 * hi_49[k]
                  + f_0 * ki_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hi_50, hi_51, hi_52, hi_53, hi_54, \
                         ki_50, ki_51, ki_52, ki_53, ki_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -5.0 * hi_50[k]
                  + f_0 * ki_50[k];

        t_51[k] = -5.0 * hi_51[k]
                  + f_0 * ki_51[k];

        t_52[k] = -5.0 * hi_52[k]
                  + f_0 * ki_52[k];

        t_53[k] = -5.0 * hi_53[k]
                  + f_0 * ki_53[k];

        t_54[k] = -5.0 * hi_54[k]
                  + f_0 * ki_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, hi_55, hi_56, hi_57, hi_58, hi_59, \
                         ki_55, ki_56, ki_57, ki_58, ki_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -5.0 * hi_55[k]
                  + f_0 * ki_55[k];

        t_56[k] = -5.0 * hi_56[k]
                  + f_0 * ki_56[k];

        t_57[k] = -5.0 * hi_57[k]
                  + f_0 * ki_57[k];

        t_58[k] = -5.0 * hi_58[k]
                  + f_0 * ki_58[k];

        t_59[k] = -5.0 * hi_59[k]
                  + f_0 * ki_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, hi_60, hi_61, hi_62, hi_63, hi_64, \
                         ki_60, ki_61, ki_62, ki_63, ki_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -5.0 * hi_60[k]
                  + f_0 * ki_60[k];

        t_61[k] = -5.0 * hi_61[k]
                  + f_0 * ki_61[k];

        t_62[k] = -5.0 * hi_62[k]
                  + f_0 * ki_62[k];

        t_63[k] = -5.0 * hi_63[k]
                  + f_0 * ki_63[k];

        t_64[k] = -5.0 * hi_64[k]
                  + f_0 * ki_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, hi_65, hi_66, hi_67, hi_68, hi_69, \
                         ki_65, ki_66, ki_67, ki_68, ki_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -5.0 * hi_65[k]
                  + f_0 * ki_65[k];

        t_66[k] = -5.0 * hi_66[k]
                  + f_0 * ki_66[k];

        t_67[k] = -5.0 * hi_67[k]
                  + f_0 * ki_67[k];

        t_68[k] = -5.0 * hi_68[k]
                  + f_0 * ki_68[k];

        t_69[k] = -5.0 * hi_69[k]
                  + f_0 * ki_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, hi_70, hi_71, hi_72, hi_73, hi_74, \
                         ki_70, ki_71, ki_72, ki_73, ki_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -5.0 * hi_70[k]
                  + f_0 * ki_70[k];

        t_71[k] = -5.0 * hi_71[k]
                  + f_0 * ki_71[k];

        t_72[k] = -5.0 * hi_72[k]
                  + f_0 * ki_72[k];

        t_73[k] = -5.0 * hi_73[k]
                  + f_0 * ki_73[k];

        t_74[k] = -5.0 * hi_74[k]
                  + f_0 * ki_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, hi_75, hi_76, hi_77, hi_78, hi_79, \
                         ki_75, ki_76, ki_77, ki_78, ki_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -5.0 * hi_75[k]
                  + f_0 * ki_75[k];

        t_76[k] = -5.0 * hi_76[k]
                  + f_0 * ki_76[k];

        t_77[k] = -5.0 * hi_77[k]
                  + f_0 * ki_77[k];

        t_78[k] = -5.0 * hi_78[k]
                  + f_0 * ki_78[k];

        t_79[k] = -5.0 * hi_79[k]
                  + f_0 * ki_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, hi_80, hi_81, hi_82, hi_83, hi_84, \
                         ki_80, ki_81, ki_82, ki_83, ki_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -5.0 * hi_80[k]
                  + f_0 * ki_80[k];

        t_81[k] = -5.0 * hi_81[k]
                  + f_0 * ki_81[k];

        t_82[k] = -5.0 * hi_82[k]
                  + f_0 * ki_82[k];

        t_83[k] = -5.0 * hi_83[k]
                  + f_0 * ki_83[k];

        t_84[k] = -4.0 * hi_84[k]
                  + f_0 * ki_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, hi_85, hi_86, hi_87, hi_88, hi_89, \
                         ki_85, ki_86, ki_87, ki_88, ki_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * hi_85[k]
                  + f_0 * ki_85[k];

        t_86[k] = -4.0 * hi_86[k]
                  + f_0 * ki_86[k];

        t_87[k] = -4.0 * hi_87[k]
                  + f_0 * ki_87[k];

        t_88[k] = -4.0 * hi_88[k]
                  + f_0 * ki_88[k];

        t_89[k] = -4.0 * hi_89[k]
                  + f_0 * ki_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, hi_90, hi_91, hi_92, hi_93, hi_94, \
                         ki_90, ki_91, ki_92, ki_93, ki_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -4.0 * hi_90[k]
                  + f_0 * ki_90[k];

        t_91[k] = -4.0 * hi_91[k]
                  + f_0 * ki_91[k];

        t_92[k] = -4.0 * hi_92[k]
                  + f_0 * ki_92[k];

        t_93[k] = -4.0 * hi_93[k]
                  + f_0 * ki_93[k];

        t_94[k] = -4.0 * hi_94[k]
                  + f_0 * ki_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, hi_95, hi_96, hi_97, hi_98, hi_99, \
                         ki_95, ki_96, ki_97, ki_98, ki_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -4.0 * hi_95[k]
                  + f_0 * ki_95[k];

        t_96[k] = -4.0 * hi_96[k]
                  + f_0 * ki_96[k];

        t_97[k] = -4.0 * hi_97[k]
                  + f_0 * ki_97[k];

        t_98[k] = -4.0 * hi_98[k]
                  + f_0 * ki_98[k];

        t_99[k] = -4.0 * hi_99[k]
                  + f_0 * ki_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, hi_100, hi_101, hi_102, hi_103, \
                         hi_104, ki_100, ki_101, ki_102, ki_103, \
                         ki_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -4.0 * hi_100[k]
                   + f_0 * ki_100[k];

        t_101[k] = -4.0 * hi_101[k]
                   + f_0 * ki_101[k];

        t_102[k] = -4.0 * hi_102[k]
                   + f_0 * ki_102[k];

        t_103[k] = -4.0 * hi_103[k]
                   + f_0 * ki_103[k];

        t_104[k] = -4.0 * hi_104[k]
                   + f_0 * ki_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, hi_105, hi_106, hi_107, hi_108, \
                         hi_109, ki_105, ki_106, ki_107, ki_108, \
                         ki_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -4.0 * hi_105[k]
                   + f_0 * ki_105[k];

        t_106[k] = -4.0 * hi_106[k]
                   + f_0 * ki_106[k];

        t_107[k] = -4.0 * hi_107[k]
                   + f_0 * ki_107[k];

        t_108[k] = -4.0 * hi_108[k]
                   + f_0 * ki_108[k];

        t_109[k] = -4.0 * hi_109[k]
                   + f_0 * ki_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, hi_110, hi_111, hi_112, hi_113, \
                         hi_114, ki_110, ki_111, ki_112, ki_113, \
                         ki_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -4.0 * hi_110[k]
                   + f_0 * ki_110[k];

        t_111[k] = -4.0 * hi_111[k]
                   + f_0 * ki_111[k];

        t_112[k] = -4.0 * hi_112[k]
                   + f_0 * ki_112[k];

        t_113[k] = -4.0 * hi_113[k]
                   + f_0 * ki_113[k];

        t_114[k] = -4.0 * hi_114[k]
                   + f_0 * ki_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, hi_115, hi_116, hi_117, hi_118, \
                         hi_119, ki_115, ki_116, ki_117, ki_118, \
                         ki_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -4.0 * hi_115[k]
                   + f_0 * ki_115[k];

        t_116[k] = -4.0 * hi_116[k]
                   + f_0 * ki_116[k];

        t_117[k] = -4.0 * hi_117[k]
                   + f_0 * ki_117[k];

        t_118[k] = -4.0 * hi_118[k]
                   + f_0 * ki_118[k];

        t_119[k] = -4.0 * hi_119[k]
                   + f_0 * ki_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, hi_120, hi_121, hi_122, hi_123, \
                         hi_124, ki_120, ki_121, ki_122, ki_123, \
                         ki_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -4.0 * hi_120[k]
                   + f_0 * ki_120[k];

        t_121[k] = -4.0 * hi_121[k]
                   + f_0 * ki_121[k];

        t_122[k] = -4.0 * hi_122[k]
                   + f_0 * ki_122[k];

        t_123[k] = -4.0 * hi_123[k]
                   + f_0 * ki_123[k];

        t_124[k] = -4.0 * hi_124[k]
                   + f_0 * ki_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, hi_125, hi_126, hi_127, hi_128, \
                         hi_129, ki_125, ki_126, ki_127, ki_128, \
                         ki_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -4.0 * hi_125[k]
                   + f_0 * ki_125[k];

        t_126[k] = -4.0 * hi_126[k]
                   + f_0 * ki_126[k];

        t_127[k] = -4.0 * hi_127[k]
                   + f_0 * ki_127[k];

        t_128[k] = -4.0 * hi_128[k]
                   + f_0 * ki_128[k];

        t_129[k] = -4.0 * hi_129[k]
                   + f_0 * ki_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, hi_130, hi_131, hi_132, hi_133, \
                         hi_134, ki_130, ki_131, ki_132, ki_133, \
                         ki_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -4.0 * hi_130[k]
                   + f_0 * ki_130[k];

        t_131[k] = -4.0 * hi_131[k]
                   + f_0 * ki_131[k];

        t_132[k] = -4.0 * hi_132[k]
                   + f_0 * ki_132[k];

        t_133[k] = -4.0 * hi_133[k]
                   + f_0 * ki_133[k];

        t_134[k] = -4.0 * hi_134[k]
                   + f_0 * ki_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, hi_135, hi_136, hi_137, hi_138, \
                         hi_139, ki_135, ki_136, ki_137, ki_138, \
                         ki_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -4.0 * hi_135[k]
                   + f_0 * ki_135[k];

        t_136[k] = -4.0 * hi_136[k]
                   + f_0 * ki_136[k];

        t_137[k] = -4.0 * hi_137[k]
                   + f_0 * ki_137[k];

        t_138[k] = -4.0 * hi_138[k]
                   + f_0 * ki_138[k];

        t_139[k] = -4.0 * hi_139[k]
                   + f_0 * ki_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, hi_140, hi_141, hi_142, hi_143, \
                         hi_144, ki_140, ki_141, ki_142, ki_143, \
                         ki_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -4.0 * hi_140[k]
                   + f_0 * ki_140[k];

        t_141[k] = -4.0 * hi_141[k]
                   + f_0 * ki_141[k];

        t_142[k] = -4.0 * hi_142[k]
                   + f_0 * ki_142[k];

        t_143[k] = -4.0 * hi_143[k]
                   + f_0 * ki_143[k];

        t_144[k] = -4.0 * hi_144[k]
                   + f_0 * ki_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, hi_145, hi_146, hi_147, hi_148, \
                         hi_149, ki_145, ki_146, ki_147, ki_148, \
                         ki_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -4.0 * hi_145[k]
                   + f_0 * ki_145[k];

        t_146[k] = -4.0 * hi_146[k]
                   + f_0 * ki_146[k];

        t_147[k] = -4.0 * hi_147[k]
                   + f_0 * ki_147[k];

        t_148[k] = -4.0 * hi_148[k]
                   + f_0 * ki_148[k];

        t_149[k] = -4.0 * hi_149[k]
                   + f_0 * ki_149[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
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

    const auto *hi_150 = buffer.data(hi + 150);
    const auto *hi_151 = buffer.data(hi + 151);
    const auto *hi_152 = buffer.data(hi + 152);
    const auto *hi_153 = buffer.data(hi + 153);
    const auto *hi_154 = buffer.data(hi + 154);
    const auto *hi_155 = buffer.data(hi + 155);
    const auto *hi_156 = buffer.data(hi + 156);
    const auto *hi_157 = buffer.data(hi + 157);
    const auto *hi_158 = buffer.data(hi + 158);
    const auto *hi_159 = buffer.data(hi + 159);
    const auto *hi_160 = buffer.data(hi + 160);
    const auto *hi_161 = buffer.data(hi + 161);
    const auto *hi_162 = buffer.data(hi + 162);
    const auto *hi_163 = buffer.data(hi + 163);
    const auto *hi_164 = buffer.data(hi + 164);
    const auto *hi_165 = buffer.data(hi + 165);
    const auto *hi_166 = buffer.data(hi + 166);
    const auto *hi_167 = buffer.data(hi + 167);
    const auto *hi_168 = buffer.data(hi + 168);
    const auto *hi_169 = buffer.data(hi + 169);
    const auto *hi_170 = buffer.data(hi + 170);
    const auto *hi_171 = buffer.data(hi + 171);
    const auto *hi_172 = buffer.data(hi + 172);
    const auto *hi_173 = buffer.data(hi + 173);
    const auto *hi_174 = buffer.data(hi + 174);
    const auto *hi_175 = buffer.data(hi + 175);
    const auto *hi_176 = buffer.data(hi + 176);
    const auto *hi_177 = buffer.data(hi + 177);
    const auto *hi_178 = buffer.data(hi + 178);
    const auto *hi_179 = buffer.data(hi + 179);
    const auto *hi_180 = buffer.data(hi + 180);
    const auto *hi_181 = buffer.data(hi + 181);
    const auto *hi_182 = buffer.data(hi + 182);
    const auto *hi_183 = buffer.data(hi + 183);
    const auto *hi_184 = buffer.data(hi + 184);
    const auto *hi_185 = buffer.data(hi + 185);
    const auto *hi_186 = buffer.data(hi + 186);
    const auto *hi_187 = buffer.data(hi + 187);
    const auto *hi_188 = buffer.data(hi + 188);
    const auto *hi_189 = buffer.data(hi + 189);
    const auto *hi_190 = buffer.data(hi + 190);
    const auto *hi_191 = buffer.data(hi + 191);
    const auto *hi_192 = buffer.data(hi + 192);
    const auto *hi_193 = buffer.data(hi + 193);
    const auto *hi_194 = buffer.data(hi + 194);
    const auto *hi_195 = buffer.data(hi + 195);
    const auto *hi_196 = buffer.data(hi + 196);
    const auto *hi_197 = buffer.data(hi + 197);
    const auto *hi_198 = buffer.data(hi + 198);
    const auto *hi_199 = buffer.data(hi + 199);
    const auto *hi_200 = buffer.data(hi + 200);
    const auto *hi_201 = buffer.data(hi + 201);
    const auto *hi_202 = buffer.data(hi + 202);
    const auto *hi_203 = buffer.data(hi + 203);
    const auto *hi_204 = buffer.data(hi + 204);
    const auto *hi_205 = buffer.data(hi + 205);
    const auto *hi_206 = buffer.data(hi + 206);
    const auto *hi_207 = buffer.data(hi + 207);
    const auto *hi_208 = buffer.data(hi + 208);
    const auto *hi_209 = buffer.data(hi + 209);
    const auto *hi_210 = buffer.data(hi + 210);
    const auto *hi_211 = buffer.data(hi + 211);
    const auto *hi_212 = buffer.data(hi + 212);
    const auto *hi_213 = buffer.data(hi + 213);
    const auto *hi_214 = buffer.data(hi + 214);
    const auto *hi_215 = buffer.data(hi + 215);
    const auto *hi_216 = buffer.data(hi + 216);
    const auto *hi_217 = buffer.data(hi + 217);
    const auto *hi_218 = buffer.data(hi + 218);
    const auto *hi_219 = buffer.data(hi + 219);
    const auto *hi_220 = buffer.data(hi + 220);
    const auto *hi_221 = buffer.data(hi + 221);
    const auto *hi_222 = buffer.data(hi + 222);
    const auto *hi_223 = buffer.data(hi + 223);
    const auto *hi_224 = buffer.data(hi + 224);
    const auto *hi_225 = buffer.data(hi + 225);
    const auto *hi_226 = buffer.data(hi + 226);
    const auto *hi_227 = buffer.data(hi + 227);
    const auto *hi_228 = buffer.data(hi + 228);
    const auto *hi_229 = buffer.data(hi + 229);
    const auto *hi_230 = buffer.data(hi + 230);
    const auto *hi_231 = buffer.data(hi + 231);
    const auto *hi_232 = buffer.data(hi + 232);
    const auto *hi_233 = buffer.data(hi + 233);
    const auto *hi_234 = buffer.data(hi + 234);
    const auto *hi_235 = buffer.data(hi + 235);
    const auto *hi_236 = buffer.data(hi + 236);
    const auto *hi_237 = buffer.data(hi + 237);
    const auto *hi_238 = buffer.data(hi + 238);
    const auto *hi_239 = buffer.data(hi + 239);
    const auto *hi_240 = buffer.data(hi + 240);
    const auto *hi_241 = buffer.data(hi + 241);
    const auto *hi_242 = buffer.data(hi + 242);
    const auto *hi_243 = buffer.data(hi + 243);
    const auto *hi_244 = buffer.data(hi + 244);
    const auto *hi_245 = buffer.data(hi + 245);
    const auto *hi_246 = buffer.data(hi + 246);
    const auto *hi_247 = buffer.data(hi + 247);
    const auto *hi_248 = buffer.data(hi + 248);
    const auto *hi_249 = buffer.data(hi + 249);
    const auto *hi_250 = buffer.data(hi + 250);
    const auto *hi_251 = buffer.data(hi + 251);
    const auto *hi_252 = buffer.data(hi + 252);
    const auto *hi_253 = buffer.data(hi + 253);
    const auto *hi_254 = buffer.data(hi + 254);
    const auto *hi_255 = buffer.data(hi + 255);
    const auto *hi_256 = buffer.data(hi + 256);
    const auto *hi_257 = buffer.data(hi + 257);
    const auto *hi_258 = buffer.data(hi + 258);
    const auto *hi_259 = buffer.data(hi + 259);
    const auto *hi_260 = buffer.data(hi + 260);
    const auto *hi_261 = buffer.data(hi + 261);
    const auto *hi_262 = buffer.data(hi + 262);
    const auto *hi_263 = buffer.data(hi + 263);
    const auto *hi_264 = buffer.data(hi + 264);
    const auto *hi_265 = buffer.data(hi + 265);
    const auto *hi_266 = buffer.data(hi + 266);
    const auto *hi_267 = buffer.data(hi + 267);
    const auto *hi_268 = buffer.data(hi + 268);
    const auto *hi_269 = buffer.data(hi + 269);
    const auto *hi_270 = buffer.data(hi + 270);
    const auto *hi_271 = buffer.data(hi + 271);
    const auto *hi_272 = buffer.data(hi + 272);
    const auto *hi_273 = buffer.data(hi + 273);
    const auto *hi_274 = buffer.data(hi + 274);
    const auto *hi_275 = buffer.data(hi + 275);
    const auto *hi_276 = buffer.data(hi + 276);
    const auto *hi_277 = buffer.data(hi + 277);
    const auto *hi_278 = buffer.data(hi + 278);
    const auto *hi_279 = buffer.data(hi + 279);
    const auto *hi_280 = buffer.data(hi + 280);
    const auto *hi_281 = buffer.data(hi + 281);
    const auto *hi_282 = buffer.data(hi + 282);
    const auto *hi_283 = buffer.data(hi + 283);
    const auto *hi_284 = buffer.data(hi + 284);
    const auto *hi_285 = buffer.data(hi + 285);
    const auto *hi_286 = buffer.data(hi + 286);
    const auto *hi_287 = buffer.data(hi + 287);
    const auto *hi_288 = buffer.data(hi + 288);
    const auto *hi_289 = buffer.data(hi + 289);
    const auto *hi_290 = buffer.data(hi + 290);
    const auto *hi_291 = buffer.data(hi + 291);
    const auto *hi_292 = buffer.data(hi + 292);
    const auto *hi_293 = buffer.data(hi + 293);
    const auto *hi_294 = buffer.data(hi + 294);
    const auto *hi_295 = buffer.data(hi + 295);
    const auto *hi_296 = buffer.data(hi + 296);
    const auto *hi_297 = buffer.data(hi + 297);
    const auto *hi_298 = buffer.data(hi + 298);
    const auto *hi_299 = buffer.data(hi + 299);

    const auto *ki_150 = buffer.data(ki + 150);
    const auto *ki_151 = buffer.data(ki + 151);
    const auto *ki_152 = buffer.data(ki + 152);
    const auto *ki_153 = buffer.data(ki + 153);
    const auto *ki_154 = buffer.data(ki + 154);
    const auto *ki_155 = buffer.data(ki + 155);
    const auto *ki_156 = buffer.data(ki + 156);
    const auto *ki_157 = buffer.data(ki + 157);
    const auto *ki_158 = buffer.data(ki + 158);
    const auto *ki_159 = buffer.data(ki + 159);
    const auto *ki_160 = buffer.data(ki + 160);
    const auto *ki_161 = buffer.data(ki + 161);
    const auto *ki_162 = buffer.data(ki + 162);
    const auto *ki_163 = buffer.data(ki + 163);
    const auto *ki_164 = buffer.data(ki + 164);
    const auto *ki_165 = buffer.data(ki + 165);
    const auto *ki_166 = buffer.data(ki + 166);
    const auto *ki_167 = buffer.data(ki + 167);
    const auto *ki_168 = buffer.data(ki + 168);
    const auto *ki_169 = buffer.data(ki + 169);
    const auto *ki_170 = buffer.data(ki + 170);
    const auto *ki_171 = buffer.data(ki + 171);
    const auto *ki_172 = buffer.data(ki + 172);
    const auto *ki_173 = buffer.data(ki + 173);
    const auto *ki_174 = buffer.data(ki + 174);
    const auto *ki_175 = buffer.data(ki + 175);
    const auto *ki_176 = buffer.data(ki + 176);
    const auto *ki_177 = buffer.data(ki + 177);
    const auto *ki_178 = buffer.data(ki + 178);
    const auto *ki_179 = buffer.data(ki + 179);
    const auto *ki_180 = buffer.data(ki + 180);
    const auto *ki_181 = buffer.data(ki + 181);
    const auto *ki_182 = buffer.data(ki + 182);
    const auto *ki_183 = buffer.data(ki + 183);
    const auto *ki_184 = buffer.data(ki + 184);
    const auto *ki_185 = buffer.data(ki + 185);
    const auto *ki_186 = buffer.data(ki + 186);
    const auto *ki_187 = buffer.data(ki + 187);
    const auto *ki_188 = buffer.data(ki + 188);
    const auto *ki_189 = buffer.data(ki + 189);
    const auto *ki_190 = buffer.data(ki + 190);
    const auto *ki_191 = buffer.data(ki + 191);
    const auto *ki_192 = buffer.data(ki + 192);
    const auto *ki_193 = buffer.data(ki + 193);
    const auto *ki_194 = buffer.data(ki + 194);
    const auto *ki_195 = buffer.data(ki + 195);
    const auto *ki_196 = buffer.data(ki + 196);
    const auto *ki_197 = buffer.data(ki + 197);
    const auto *ki_198 = buffer.data(ki + 198);
    const auto *ki_199 = buffer.data(ki + 199);
    const auto *ki_200 = buffer.data(ki + 200);
    const auto *ki_201 = buffer.data(ki + 201);
    const auto *ki_202 = buffer.data(ki + 202);
    const auto *ki_203 = buffer.data(ki + 203);
    const auto *ki_204 = buffer.data(ki + 204);
    const auto *ki_205 = buffer.data(ki + 205);
    const auto *ki_206 = buffer.data(ki + 206);
    const auto *ki_207 = buffer.data(ki + 207);
    const auto *ki_208 = buffer.data(ki + 208);
    const auto *ki_209 = buffer.data(ki + 209);
    const auto *ki_210 = buffer.data(ki + 210);
    const auto *ki_211 = buffer.data(ki + 211);
    const auto *ki_212 = buffer.data(ki + 212);
    const auto *ki_213 = buffer.data(ki + 213);
    const auto *ki_214 = buffer.data(ki + 214);
    const auto *ki_215 = buffer.data(ki + 215);
    const auto *ki_216 = buffer.data(ki + 216);
    const auto *ki_217 = buffer.data(ki + 217);
    const auto *ki_218 = buffer.data(ki + 218);
    const auto *ki_219 = buffer.data(ki + 219);
    const auto *ki_220 = buffer.data(ki + 220);
    const auto *ki_221 = buffer.data(ki + 221);
    const auto *ki_222 = buffer.data(ki + 222);
    const auto *ki_223 = buffer.data(ki + 223);
    const auto *ki_224 = buffer.data(ki + 224);
    const auto *ki_225 = buffer.data(ki + 225);
    const auto *ki_226 = buffer.data(ki + 226);
    const auto *ki_227 = buffer.data(ki + 227);
    const auto *ki_228 = buffer.data(ki + 228);
    const auto *ki_229 = buffer.data(ki + 229);
    const auto *ki_230 = buffer.data(ki + 230);
    const auto *ki_231 = buffer.data(ki + 231);
    const auto *ki_232 = buffer.data(ki + 232);
    const auto *ki_233 = buffer.data(ki + 233);
    const auto *ki_234 = buffer.data(ki + 234);
    const auto *ki_235 = buffer.data(ki + 235);
    const auto *ki_236 = buffer.data(ki + 236);
    const auto *ki_237 = buffer.data(ki + 237);
    const auto *ki_238 = buffer.data(ki + 238);
    const auto *ki_239 = buffer.data(ki + 239);
    const auto *ki_240 = buffer.data(ki + 240);
    const auto *ki_241 = buffer.data(ki + 241);
    const auto *ki_242 = buffer.data(ki + 242);
    const auto *ki_243 = buffer.data(ki + 243);
    const auto *ki_244 = buffer.data(ki + 244);
    const auto *ki_245 = buffer.data(ki + 245);
    const auto *ki_246 = buffer.data(ki + 246);
    const auto *ki_247 = buffer.data(ki + 247);
    const auto *ki_248 = buffer.data(ki + 248);
    const auto *ki_249 = buffer.data(ki + 249);
    const auto *ki_250 = buffer.data(ki + 250);
    const auto *ki_251 = buffer.data(ki + 251);
    const auto *ki_252 = buffer.data(ki + 252);
    const auto *ki_253 = buffer.data(ki + 253);
    const auto *ki_254 = buffer.data(ki + 254);
    const auto *ki_255 = buffer.data(ki + 255);
    const auto *ki_256 = buffer.data(ki + 256);
    const auto *ki_257 = buffer.data(ki + 257);
    const auto *ki_258 = buffer.data(ki + 258);
    const auto *ki_259 = buffer.data(ki + 259);
    const auto *ki_260 = buffer.data(ki + 260);
    const auto *ki_261 = buffer.data(ki + 261);
    const auto *ki_262 = buffer.data(ki + 262);
    const auto *ki_263 = buffer.data(ki + 263);
    const auto *ki_264 = buffer.data(ki + 264);
    const auto *ki_265 = buffer.data(ki + 265);
    const auto *ki_266 = buffer.data(ki + 266);
    const auto *ki_267 = buffer.data(ki + 267);
    const auto *ki_268 = buffer.data(ki + 268);
    const auto *ki_269 = buffer.data(ki + 269);
    const auto *ki_270 = buffer.data(ki + 270);
    const auto *ki_271 = buffer.data(ki + 271);
    const auto *ki_272 = buffer.data(ki + 272);
    const auto *ki_273 = buffer.data(ki + 273);
    const auto *ki_274 = buffer.data(ki + 274);
    const auto *ki_275 = buffer.data(ki + 275);
    const auto *ki_276 = buffer.data(ki + 276);
    const auto *ki_277 = buffer.data(ki + 277);
    const auto *ki_278 = buffer.data(ki + 278);
    const auto *ki_279 = buffer.data(ki + 279);
    const auto *ki_280 = buffer.data(ki + 280);
    const auto *ki_281 = buffer.data(ki + 281);
    const auto *ki_282 = buffer.data(ki + 282);
    const auto *ki_283 = buffer.data(ki + 283);
    const auto *ki_284 = buffer.data(ki + 284);
    const auto *ki_285 = buffer.data(ki + 285);
    const auto *ki_286 = buffer.data(ki + 286);
    const auto *ki_287 = buffer.data(ki + 287);
    const auto *ki_288 = buffer.data(ki + 288);
    const auto *ki_289 = buffer.data(ki + 289);
    const auto *ki_290 = buffer.data(ki + 290);
    const auto *ki_291 = buffer.data(ki + 291);
    const auto *ki_292 = buffer.data(ki + 292);
    const auto *ki_293 = buffer.data(ki + 293);
    const auto *ki_294 = buffer.data(ki + 294);
    const auto *ki_295 = buffer.data(ki + 295);
    const auto *ki_296 = buffer.data(ki + 296);
    const auto *ki_297 = buffer.data(ki + 297);
    const auto *ki_298 = buffer.data(ki + 298);
    const auto *ki_299 = buffer.data(ki + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, hi_150, hi_151, hi_152, hi_153, \
                         hi_154, ki_150, ki_151, ki_152, ki_153, \
                         ki_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -4.0 * hi_150[k]
                   + f_0 * ki_150[k];

        t_151[k] = -4.0 * hi_151[k]
                   + f_0 * ki_151[k];

        t_152[k] = -4.0 * hi_152[k]
                   + f_0 * ki_152[k];

        t_153[k] = -4.0 * hi_153[k]
                   + f_0 * ki_153[k];

        t_154[k] = -4.0 * hi_154[k]
                   + f_0 * ki_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, hi_155, hi_156, hi_157, hi_158, \
                         hi_159, ki_155, ki_156, ki_157, ki_158, \
                         ki_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -4.0 * hi_155[k]
                   + f_0 * ki_155[k];

        t_156[k] = -4.0 * hi_156[k]
                   + f_0 * ki_156[k];

        t_157[k] = -4.0 * hi_157[k]
                   + f_0 * ki_157[k];

        t_158[k] = -4.0 * hi_158[k]
                   + f_0 * ki_158[k];

        t_159[k] = -4.0 * hi_159[k]
                   + f_0 * ki_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, hi_160, hi_161, hi_162, hi_163, \
                         hi_164, ki_160, ki_161, ki_162, ki_163, \
                         ki_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -4.0 * hi_160[k]
                   + f_0 * ki_160[k];

        t_161[k] = -4.0 * hi_161[k]
                   + f_0 * ki_161[k];

        t_162[k] = -4.0 * hi_162[k]
                   + f_0 * ki_162[k];

        t_163[k] = -4.0 * hi_163[k]
                   + f_0 * ki_163[k];

        t_164[k] = -4.0 * hi_164[k]
                   + f_0 * ki_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, hi_165, hi_166, hi_167, hi_168, \
                         hi_169, ki_165, ki_166, ki_167, ki_168, \
                         ki_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -4.0 * hi_165[k]
                   + f_0 * ki_165[k];

        t_166[k] = -4.0 * hi_166[k]
                   + f_0 * ki_166[k];

        t_167[k] = -4.0 * hi_167[k]
                   + f_0 * ki_167[k];

        t_168[k] = -3.0 * hi_168[k]
                   + f_0 * ki_168[k];

        t_169[k] = -3.0 * hi_169[k]
                   + f_0 * ki_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, hi_170, hi_171, hi_172, hi_173, \
                         hi_174, ki_170, ki_171, ki_172, ki_173, \
                         ki_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * hi_170[k]
                   + f_0 * ki_170[k];

        t_171[k] = -3.0 * hi_171[k]
                   + f_0 * ki_171[k];

        t_172[k] = -3.0 * hi_172[k]
                   + f_0 * ki_172[k];

        t_173[k] = -3.0 * hi_173[k]
                   + f_0 * ki_173[k];

        t_174[k] = -3.0 * hi_174[k]
                   + f_0 * ki_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, hi_175, hi_176, hi_177, hi_178, \
                         hi_179, ki_175, ki_176, ki_177, ki_178, \
                         ki_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * hi_175[k]
                   + f_0 * ki_175[k];

        t_176[k] = -3.0 * hi_176[k]
                   + f_0 * ki_176[k];

        t_177[k] = -3.0 * hi_177[k]
                   + f_0 * ki_177[k];

        t_178[k] = -3.0 * hi_178[k]
                   + f_0 * ki_178[k];

        t_179[k] = -3.0 * hi_179[k]
                   + f_0 * ki_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, hi_180, hi_181, hi_182, hi_183, \
                         hi_184, ki_180, ki_181, ki_182, ki_183, \
                         ki_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * hi_180[k]
                   + f_0 * ki_180[k];

        t_181[k] = -3.0 * hi_181[k]
                   + f_0 * ki_181[k];

        t_182[k] = -3.0 * hi_182[k]
                   + f_0 * ki_182[k];

        t_183[k] = -3.0 * hi_183[k]
                   + f_0 * ki_183[k];

        t_184[k] = -3.0 * hi_184[k]
                   + f_0 * ki_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, hi_185, hi_186, hi_187, hi_188, \
                         hi_189, ki_185, ki_186, ki_187, ki_188, \
                         ki_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * hi_185[k]
                   + f_0 * ki_185[k];

        t_186[k] = -3.0 * hi_186[k]
                   + f_0 * ki_186[k];

        t_187[k] = -3.0 * hi_187[k]
                   + f_0 * ki_187[k];

        t_188[k] = -3.0 * hi_188[k]
                   + f_0 * ki_188[k];

        t_189[k] = -3.0 * hi_189[k]
                   + f_0 * ki_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, hi_190, hi_191, hi_192, hi_193, \
                         hi_194, ki_190, ki_191, ki_192, ki_193, \
                         ki_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * hi_190[k]
                   + f_0 * ki_190[k];

        t_191[k] = -3.0 * hi_191[k]
                   + f_0 * ki_191[k];

        t_192[k] = -3.0 * hi_192[k]
                   + f_0 * ki_192[k];

        t_193[k] = -3.0 * hi_193[k]
                   + f_0 * ki_193[k];

        t_194[k] = -3.0 * hi_194[k]
                   + f_0 * ki_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, hi_195, hi_196, hi_197, hi_198, \
                         hi_199, ki_195, ki_196, ki_197, ki_198, \
                         ki_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * hi_195[k]
                   + f_0 * ki_195[k];

        t_196[k] = -3.0 * hi_196[k]
                   + f_0 * ki_196[k];

        t_197[k] = -3.0 * hi_197[k]
                   + f_0 * ki_197[k];

        t_198[k] = -3.0 * hi_198[k]
                   + f_0 * ki_198[k];

        t_199[k] = -3.0 * hi_199[k]
                   + f_0 * ki_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, hi_200, hi_201, hi_202, hi_203, \
                         hi_204, ki_200, ki_201, ki_202, ki_203, \
                         ki_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -3.0 * hi_200[k]
                   + f_0 * ki_200[k];

        t_201[k] = -3.0 * hi_201[k]
                   + f_0 * ki_201[k];

        t_202[k] = -3.0 * hi_202[k]
                   + f_0 * ki_202[k];

        t_203[k] = -3.0 * hi_203[k]
                   + f_0 * ki_203[k];

        t_204[k] = -3.0 * hi_204[k]
                   + f_0 * ki_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, hi_205, hi_206, hi_207, hi_208, \
                         hi_209, ki_205, ki_206, ki_207, ki_208, \
                         ki_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -3.0 * hi_205[k]
                   + f_0 * ki_205[k];

        t_206[k] = -3.0 * hi_206[k]
                   + f_0 * ki_206[k];

        t_207[k] = -3.0 * hi_207[k]
                   + f_0 * ki_207[k];

        t_208[k] = -3.0 * hi_208[k]
                   + f_0 * ki_208[k];

        t_209[k] = -3.0 * hi_209[k]
                   + f_0 * ki_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, hi_210, hi_211, hi_212, hi_213, \
                         hi_214, ki_210, ki_211, ki_212, ki_213, \
                         ki_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -3.0 * hi_210[k]
                   + f_0 * ki_210[k];

        t_211[k] = -3.0 * hi_211[k]
                   + f_0 * ki_211[k];

        t_212[k] = -3.0 * hi_212[k]
                   + f_0 * ki_212[k];

        t_213[k] = -3.0 * hi_213[k]
                   + f_0 * ki_213[k];

        t_214[k] = -3.0 * hi_214[k]
                   + f_0 * ki_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, hi_215, hi_216, hi_217, hi_218, \
                         hi_219, ki_215, ki_216, ki_217, ki_218, \
                         ki_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -3.0 * hi_215[k]
                   + f_0 * ki_215[k];

        t_216[k] = -3.0 * hi_216[k]
                   + f_0 * ki_216[k];

        t_217[k] = -3.0 * hi_217[k]
                   + f_0 * ki_217[k];

        t_218[k] = -3.0 * hi_218[k]
                   + f_0 * ki_218[k];

        t_219[k] = -3.0 * hi_219[k]
                   + f_0 * ki_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, hi_220, hi_221, hi_222, hi_223, \
                         hi_224, ki_220, ki_221, ki_222, ki_223, \
                         ki_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -3.0 * hi_220[k]
                   + f_0 * ki_220[k];

        t_221[k] = -3.0 * hi_221[k]
                   + f_0 * ki_221[k];

        t_222[k] = -3.0 * hi_222[k]
                   + f_0 * ki_222[k];

        t_223[k] = -3.0 * hi_223[k]
                   + f_0 * ki_223[k];

        t_224[k] = -3.0 * hi_224[k]
                   + f_0 * ki_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, hi_225, hi_226, hi_227, hi_228, \
                         hi_229, ki_225, ki_226, ki_227, ki_228, \
                         ki_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -3.0 * hi_225[k]
                   + f_0 * ki_225[k];

        t_226[k] = -3.0 * hi_226[k]
                   + f_0 * ki_226[k];

        t_227[k] = -3.0 * hi_227[k]
                   + f_0 * ki_227[k];

        t_228[k] = -3.0 * hi_228[k]
                   + f_0 * ki_228[k];

        t_229[k] = -3.0 * hi_229[k]
                   + f_0 * ki_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, hi_230, hi_231, hi_232, hi_233, \
                         hi_234, ki_230, ki_231, ki_232, ki_233, \
                         ki_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -3.0 * hi_230[k]
                   + f_0 * ki_230[k];

        t_231[k] = -3.0 * hi_231[k]
                   + f_0 * ki_231[k];

        t_232[k] = -3.0 * hi_232[k]
                   + f_0 * ki_232[k];

        t_233[k] = -3.0 * hi_233[k]
                   + f_0 * ki_233[k];

        t_234[k] = -3.0 * hi_234[k]
                   + f_0 * ki_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, hi_235, hi_236, hi_237, hi_238, \
                         hi_239, ki_235, ki_236, ki_237, ki_238, \
                         ki_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -3.0 * hi_235[k]
                   + f_0 * ki_235[k];

        t_236[k] = -3.0 * hi_236[k]
                   + f_0 * ki_236[k];

        t_237[k] = -3.0 * hi_237[k]
                   + f_0 * ki_237[k];

        t_238[k] = -3.0 * hi_238[k]
                   + f_0 * ki_238[k];

        t_239[k] = -3.0 * hi_239[k]
                   + f_0 * ki_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, hi_240, hi_241, hi_242, hi_243, \
                         hi_244, ki_240, ki_241, ki_242, ki_243, \
                         ki_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * hi_240[k]
                   + f_0 * ki_240[k];

        t_241[k] = -3.0 * hi_241[k]
                   + f_0 * ki_241[k];

        t_242[k] = -3.0 * hi_242[k]
                   + f_0 * ki_242[k];

        t_243[k] = -3.0 * hi_243[k]
                   + f_0 * ki_243[k];

        t_244[k] = -3.0 * hi_244[k]
                   + f_0 * ki_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, hi_245, hi_246, hi_247, hi_248, \
                         hi_249, ki_245, ki_246, ki_247, ki_248, \
                         ki_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * hi_245[k]
                   + f_0 * ki_245[k];

        t_246[k] = -3.0 * hi_246[k]
                   + f_0 * ki_246[k];

        t_247[k] = -3.0 * hi_247[k]
                   + f_0 * ki_247[k];

        t_248[k] = -3.0 * hi_248[k]
                   + f_0 * ki_248[k];

        t_249[k] = -3.0 * hi_249[k]
                   + f_0 * ki_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, hi_250, hi_251, hi_252, hi_253, \
                         hi_254, ki_250, ki_251, ki_252, ki_253, \
                         ki_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -3.0 * hi_250[k]
                   + f_0 * ki_250[k];

        t_251[k] = -3.0 * hi_251[k]
                   + f_0 * ki_251[k];

        t_252[k] = -3.0 * hi_252[k]
                   + f_0 * ki_252[k];

        t_253[k] = -3.0 * hi_253[k]
                   + f_0 * ki_253[k];

        t_254[k] = -3.0 * hi_254[k]
                   + f_0 * ki_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, hi_255, hi_256, hi_257, hi_258, \
                         hi_259, ki_255, ki_256, ki_257, ki_258, \
                         ki_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -3.0 * hi_255[k]
                   + f_0 * ki_255[k];

        t_256[k] = -3.0 * hi_256[k]
                   + f_0 * ki_256[k];

        t_257[k] = -3.0 * hi_257[k]
                   + f_0 * ki_257[k];

        t_258[k] = -3.0 * hi_258[k]
                   + f_0 * ki_258[k];

        t_259[k] = -3.0 * hi_259[k]
                   + f_0 * ki_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, hi_260, hi_261, hi_262, hi_263, \
                         hi_264, ki_260, ki_261, ki_262, ki_263, \
                         ki_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -3.0 * hi_260[k]
                   + f_0 * ki_260[k];

        t_261[k] = -3.0 * hi_261[k]
                   + f_0 * ki_261[k];

        t_262[k] = -3.0 * hi_262[k]
                   + f_0 * ki_262[k];

        t_263[k] = -3.0 * hi_263[k]
                   + f_0 * ki_263[k];

        t_264[k] = -3.0 * hi_264[k]
                   + f_0 * ki_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, hi_265, hi_266, hi_267, hi_268, \
                         hi_269, ki_265, ki_266, ki_267, ki_268, \
                         ki_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -3.0 * hi_265[k]
                   + f_0 * ki_265[k];

        t_266[k] = -3.0 * hi_266[k]
                   + f_0 * ki_266[k];

        t_267[k] = -3.0 * hi_267[k]
                   + f_0 * ki_267[k];

        t_268[k] = -3.0 * hi_268[k]
                   + f_0 * ki_268[k];

        t_269[k] = -3.0 * hi_269[k]
                   + f_0 * ki_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, hi_270, hi_271, hi_272, hi_273, \
                         hi_274, ki_270, ki_271, ki_272, ki_273, \
                         ki_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -3.0 * hi_270[k]
                   + f_0 * ki_270[k];

        t_271[k] = -3.0 * hi_271[k]
                   + f_0 * ki_271[k];

        t_272[k] = -3.0 * hi_272[k]
                   + f_0 * ki_272[k];

        t_273[k] = -3.0 * hi_273[k]
                   + f_0 * ki_273[k];

        t_274[k] = -3.0 * hi_274[k]
                   + f_0 * ki_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, hi_275, hi_276, hi_277, hi_278, \
                         hi_279, ki_275, ki_276, ki_277, ki_278, \
                         ki_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * hi_275[k]
                   + f_0 * ki_275[k];

        t_276[k] = -3.0 * hi_276[k]
                   + f_0 * ki_276[k];

        t_277[k] = -3.0 * hi_277[k]
                   + f_0 * ki_277[k];

        t_278[k] = -3.0 * hi_278[k]
                   + f_0 * ki_278[k];

        t_279[k] = -3.0 * hi_279[k]
                   + f_0 * ki_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, hi_280, hi_281, hi_282, hi_283, \
                         hi_284, ki_280, ki_281, ki_282, ki_283, \
                         ki_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -2.0 * hi_280[k]
                   + f_0 * ki_280[k];

        t_281[k] = -2.0 * hi_281[k]
                   + f_0 * ki_281[k];

        t_282[k] = -2.0 * hi_282[k]
                   + f_0 * ki_282[k];

        t_283[k] = -2.0 * hi_283[k]
                   + f_0 * ki_283[k];

        t_284[k] = -2.0 * hi_284[k]
                   + f_0 * ki_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, hi_285, hi_286, hi_287, hi_288, \
                         hi_289, ki_285, ki_286, ki_287, ki_288, \
                         ki_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -2.0 * hi_285[k]
                   + f_0 * ki_285[k];

        t_286[k] = -2.0 * hi_286[k]
                   + f_0 * ki_286[k];

        t_287[k] = -2.0 * hi_287[k]
                   + f_0 * ki_287[k];

        t_288[k] = -2.0 * hi_288[k]
                   + f_0 * ki_288[k];

        t_289[k] = -2.0 * hi_289[k]
                   + f_0 * ki_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, hi_290, hi_291, hi_292, hi_293, \
                         hi_294, ki_290, ki_291, ki_292, ki_293, \
                         ki_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -2.0 * hi_290[k]
                   + f_0 * ki_290[k];

        t_291[k] = -2.0 * hi_291[k]
                   + f_0 * ki_291[k];

        t_292[k] = -2.0 * hi_292[k]
                   + f_0 * ki_292[k];

        t_293[k] = -2.0 * hi_293[k]
                   + f_0 * ki_293[k];

        t_294[k] = -2.0 * hi_294[k]
                   + f_0 * ki_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, hi_295, hi_296, hi_297, hi_298, \
                         hi_299, ki_295, ki_296, ki_297, ki_298, \
                         ki_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -2.0 * hi_295[k]
                   + f_0 * ki_295[k];

        t_296[k] = -2.0 * hi_296[k]
                   + f_0 * ki_296[k];

        t_297[k] = -2.0 * hi_297[k]
                   + f_0 * ki_297[k];

        t_298[k] = -2.0 * hi_298[k]
                   + f_0 * ki_298[k];

        t_299[k] = -2.0 * hi_299[k]
                   + f_0 * ki_299[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
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

    const auto *hi_300 = buffer.data(hi + 300);
    const auto *hi_301 = buffer.data(hi + 301);
    const auto *hi_302 = buffer.data(hi + 302);
    const auto *hi_303 = buffer.data(hi + 303);
    const auto *hi_304 = buffer.data(hi + 304);
    const auto *hi_305 = buffer.data(hi + 305);
    const auto *hi_306 = buffer.data(hi + 306);
    const auto *hi_307 = buffer.data(hi + 307);
    const auto *hi_308 = buffer.data(hi + 308);
    const auto *hi_309 = buffer.data(hi + 309);
    const auto *hi_310 = buffer.data(hi + 310);
    const auto *hi_311 = buffer.data(hi + 311);
    const auto *hi_312 = buffer.data(hi + 312);
    const auto *hi_313 = buffer.data(hi + 313);
    const auto *hi_314 = buffer.data(hi + 314);
    const auto *hi_315 = buffer.data(hi + 315);
    const auto *hi_316 = buffer.data(hi + 316);
    const auto *hi_317 = buffer.data(hi + 317);
    const auto *hi_318 = buffer.data(hi + 318);
    const auto *hi_319 = buffer.data(hi + 319);
    const auto *hi_320 = buffer.data(hi + 320);
    const auto *hi_321 = buffer.data(hi + 321);
    const auto *hi_322 = buffer.data(hi + 322);
    const auto *hi_323 = buffer.data(hi + 323);
    const auto *hi_324 = buffer.data(hi + 324);
    const auto *hi_325 = buffer.data(hi + 325);
    const auto *hi_326 = buffer.data(hi + 326);
    const auto *hi_327 = buffer.data(hi + 327);
    const auto *hi_328 = buffer.data(hi + 328);
    const auto *hi_329 = buffer.data(hi + 329);
    const auto *hi_330 = buffer.data(hi + 330);
    const auto *hi_331 = buffer.data(hi + 331);
    const auto *hi_332 = buffer.data(hi + 332);
    const auto *hi_333 = buffer.data(hi + 333);
    const auto *hi_334 = buffer.data(hi + 334);
    const auto *hi_335 = buffer.data(hi + 335);
    const auto *hi_336 = buffer.data(hi + 336);
    const auto *hi_337 = buffer.data(hi + 337);
    const auto *hi_338 = buffer.data(hi + 338);
    const auto *hi_339 = buffer.data(hi + 339);
    const auto *hi_340 = buffer.data(hi + 340);
    const auto *hi_341 = buffer.data(hi + 341);
    const auto *hi_342 = buffer.data(hi + 342);
    const auto *hi_343 = buffer.data(hi + 343);
    const auto *hi_344 = buffer.data(hi + 344);
    const auto *hi_345 = buffer.data(hi + 345);
    const auto *hi_346 = buffer.data(hi + 346);
    const auto *hi_347 = buffer.data(hi + 347);
    const auto *hi_348 = buffer.data(hi + 348);
    const auto *hi_349 = buffer.data(hi + 349);
    const auto *hi_350 = buffer.data(hi + 350);
    const auto *hi_351 = buffer.data(hi + 351);
    const auto *hi_352 = buffer.data(hi + 352);
    const auto *hi_353 = buffer.data(hi + 353);
    const auto *hi_354 = buffer.data(hi + 354);
    const auto *hi_355 = buffer.data(hi + 355);
    const auto *hi_356 = buffer.data(hi + 356);
    const auto *hi_357 = buffer.data(hi + 357);
    const auto *hi_358 = buffer.data(hi + 358);
    const auto *hi_359 = buffer.data(hi + 359);
    const auto *hi_360 = buffer.data(hi + 360);
    const auto *hi_361 = buffer.data(hi + 361);
    const auto *hi_362 = buffer.data(hi + 362);
    const auto *hi_363 = buffer.data(hi + 363);
    const auto *hi_364 = buffer.data(hi + 364);
    const auto *hi_365 = buffer.data(hi + 365);
    const auto *hi_366 = buffer.data(hi + 366);
    const auto *hi_367 = buffer.data(hi + 367);
    const auto *hi_368 = buffer.data(hi + 368);
    const auto *hi_369 = buffer.data(hi + 369);
    const auto *hi_370 = buffer.data(hi + 370);
    const auto *hi_371 = buffer.data(hi + 371);
    const auto *hi_372 = buffer.data(hi + 372);
    const auto *hi_373 = buffer.data(hi + 373);
    const auto *hi_374 = buffer.data(hi + 374);
    const auto *hi_375 = buffer.data(hi + 375);
    const auto *hi_376 = buffer.data(hi + 376);
    const auto *hi_377 = buffer.data(hi + 377);
    const auto *hi_378 = buffer.data(hi + 378);
    const auto *hi_379 = buffer.data(hi + 379);
    const auto *hi_380 = buffer.data(hi + 380);
    const auto *hi_381 = buffer.data(hi + 381);
    const auto *hi_382 = buffer.data(hi + 382);
    const auto *hi_383 = buffer.data(hi + 383);
    const auto *hi_384 = buffer.data(hi + 384);
    const auto *hi_385 = buffer.data(hi + 385);
    const auto *hi_386 = buffer.data(hi + 386);
    const auto *hi_387 = buffer.data(hi + 387);
    const auto *hi_388 = buffer.data(hi + 388);
    const auto *hi_389 = buffer.data(hi + 389);
    const auto *hi_390 = buffer.data(hi + 390);
    const auto *hi_391 = buffer.data(hi + 391);
    const auto *hi_392 = buffer.data(hi + 392);
    const auto *hi_393 = buffer.data(hi + 393);
    const auto *hi_394 = buffer.data(hi + 394);
    const auto *hi_395 = buffer.data(hi + 395);
    const auto *hi_396 = buffer.data(hi + 396);
    const auto *hi_397 = buffer.data(hi + 397);
    const auto *hi_398 = buffer.data(hi + 398);
    const auto *hi_399 = buffer.data(hi + 399);
    const auto *hi_400 = buffer.data(hi + 400);
    const auto *hi_401 = buffer.data(hi + 401);
    const auto *hi_402 = buffer.data(hi + 402);
    const auto *hi_403 = buffer.data(hi + 403);
    const auto *hi_404 = buffer.data(hi + 404);
    const auto *hi_405 = buffer.data(hi + 405);
    const auto *hi_406 = buffer.data(hi + 406);
    const auto *hi_407 = buffer.data(hi + 407);
    const auto *hi_408 = buffer.data(hi + 408);
    const auto *hi_409 = buffer.data(hi + 409);
    const auto *hi_410 = buffer.data(hi + 410);
    const auto *hi_411 = buffer.data(hi + 411);
    const auto *hi_412 = buffer.data(hi + 412);
    const auto *hi_413 = buffer.data(hi + 413);
    const auto *hi_414 = buffer.data(hi + 414);
    const auto *hi_415 = buffer.data(hi + 415);
    const auto *hi_416 = buffer.data(hi + 416);
    const auto *hi_417 = buffer.data(hi + 417);
    const auto *hi_418 = buffer.data(hi + 418);
    const auto *hi_419 = buffer.data(hi + 419);
    const auto *hi_420 = buffer.data(hi + 420);
    const auto *hi_421 = buffer.data(hi + 421);
    const auto *hi_422 = buffer.data(hi + 422);
    const auto *hi_423 = buffer.data(hi + 423);
    const auto *hi_424 = buffer.data(hi + 424);
    const auto *hi_425 = buffer.data(hi + 425);
    const auto *hi_426 = buffer.data(hi + 426);
    const auto *hi_427 = buffer.data(hi + 427);
    const auto *hi_428 = buffer.data(hi + 428);
    const auto *hi_429 = buffer.data(hi + 429);
    const auto *hi_430 = buffer.data(hi + 430);
    const auto *hi_431 = buffer.data(hi + 431);
    const auto *hi_432 = buffer.data(hi + 432);
    const auto *hi_433 = buffer.data(hi + 433);
    const auto *hi_434 = buffer.data(hi + 434);
    const auto *hi_435 = buffer.data(hi + 435);
    const auto *hi_436 = buffer.data(hi + 436);
    const auto *hi_437 = buffer.data(hi + 437);
    const auto *hi_438 = buffer.data(hi + 438);
    const auto *hi_439 = buffer.data(hi + 439);
    const auto *hi_440 = buffer.data(hi + 440);
    const auto *hi_441 = buffer.data(hi + 441);
    const auto *hi_442 = buffer.data(hi + 442);
    const auto *hi_443 = buffer.data(hi + 443);
    const auto *hi_444 = buffer.data(hi + 444);
    const auto *hi_445 = buffer.data(hi + 445);
    const auto *hi_446 = buffer.data(hi + 446);
    const auto *hi_447 = buffer.data(hi + 447);
    const auto *hi_448 = buffer.data(hi + 448);
    const auto *hi_449 = buffer.data(hi + 449);

    const auto *ki_300 = buffer.data(ki + 300);
    const auto *ki_301 = buffer.data(ki + 301);
    const auto *ki_302 = buffer.data(ki + 302);
    const auto *ki_303 = buffer.data(ki + 303);
    const auto *ki_304 = buffer.data(ki + 304);
    const auto *ki_305 = buffer.data(ki + 305);
    const auto *ki_306 = buffer.data(ki + 306);
    const auto *ki_307 = buffer.data(ki + 307);
    const auto *ki_308 = buffer.data(ki + 308);
    const auto *ki_309 = buffer.data(ki + 309);
    const auto *ki_310 = buffer.data(ki + 310);
    const auto *ki_311 = buffer.data(ki + 311);
    const auto *ki_312 = buffer.data(ki + 312);
    const auto *ki_313 = buffer.data(ki + 313);
    const auto *ki_314 = buffer.data(ki + 314);
    const auto *ki_315 = buffer.data(ki + 315);
    const auto *ki_316 = buffer.data(ki + 316);
    const auto *ki_317 = buffer.data(ki + 317);
    const auto *ki_318 = buffer.data(ki + 318);
    const auto *ki_319 = buffer.data(ki + 319);
    const auto *ki_320 = buffer.data(ki + 320);
    const auto *ki_321 = buffer.data(ki + 321);
    const auto *ki_322 = buffer.data(ki + 322);
    const auto *ki_323 = buffer.data(ki + 323);
    const auto *ki_324 = buffer.data(ki + 324);
    const auto *ki_325 = buffer.data(ki + 325);
    const auto *ki_326 = buffer.data(ki + 326);
    const auto *ki_327 = buffer.data(ki + 327);
    const auto *ki_328 = buffer.data(ki + 328);
    const auto *ki_329 = buffer.data(ki + 329);
    const auto *ki_330 = buffer.data(ki + 330);
    const auto *ki_331 = buffer.data(ki + 331);
    const auto *ki_332 = buffer.data(ki + 332);
    const auto *ki_333 = buffer.data(ki + 333);
    const auto *ki_334 = buffer.data(ki + 334);
    const auto *ki_335 = buffer.data(ki + 335);
    const auto *ki_336 = buffer.data(ki + 336);
    const auto *ki_337 = buffer.data(ki + 337);
    const auto *ki_338 = buffer.data(ki + 338);
    const auto *ki_339 = buffer.data(ki + 339);
    const auto *ki_340 = buffer.data(ki + 340);
    const auto *ki_341 = buffer.data(ki + 341);
    const auto *ki_342 = buffer.data(ki + 342);
    const auto *ki_343 = buffer.data(ki + 343);
    const auto *ki_344 = buffer.data(ki + 344);
    const auto *ki_345 = buffer.data(ki + 345);
    const auto *ki_346 = buffer.data(ki + 346);
    const auto *ki_347 = buffer.data(ki + 347);
    const auto *ki_348 = buffer.data(ki + 348);
    const auto *ki_349 = buffer.data(ki + 349);
    const auto *ki_350 = buffer.data(ki + 350);
    const auto *ki_351 = buffer.data(ki + 351);
    const auto *ki_352 = buffer.data(ki + 352);
    const auto *ki_353 = buffer.data(ki + 353);
    const auto *ki_354 = buffer.data(ki + 354);
    const auto *ki_355 = buffer.data(ki + 355);
    const auto *ki_356 = buffer.data(ki + 356);
    const auto *ki_357 = buffer.data(ki + 357);
    const auto *ki_358 = buffer.data(ki + 358);
    const auto *ki_359 = buffer.data(ki + 359);
    const auto *ki_360 = buffer.data(ki + 360);
    const auto *ki_361 = buffer.data(ki + 361);
    const auto *ki_362 = buffer.data(ki + 362);
    const auto *ki_363 = buffer.data(ki + 363);
    const auto *ki_364 = buffer.data(ki + 364);
    const auto *ki_365 = buffer.data(ki + 365);
    const auto *ki_366 = buffer.data(ki + 366);
    const auto *ki_367 = buffer.data(ki + 367);
    const auto *ki_368 = buffer.data(ki + 368);
    const auto *ki_369 = buffer.data(ki + 369);
    const auto *ki_370 = buffer.data(ki + 370);
    const auto *ki_371 = buffer.data(ki + 371);
    const auto *ki_372 = buffer.data(ki + 372);
    const auto *ki_373 = buffer.data(ki + 373);
    const auto *ki_374 = buffer.data(ki + 374);
    const auto *ki_375 = buffer.data(ki + 375);
    const auto *ki_376 = buffer.data(ki + 376);
    const auto *ki_377 = buffer.data(ki + 377);
    const auto *ki_378 = buffer.data(ki + 378);
    const auto *ki_379 = buffer.data(ki + 379);
    const auto *ki_380 = buffer.data(ki + 380);
    const auto *ki_381 = buffer.data(ki + 381);
    const auto *ki_382 = buffer.data(ki + 382);
    const auto *ki_383 = buffer.data(ki + 383);
    const auto *ki_384 = buffer.data(ki + 384);
    const auto *ki_385 = buffer.data(ki + 385);
    const auto *ki_386 = buffer.data(ki + 386);
    const auto *ki_387 = buffer.data(ki + 387);
    const auto *ki_388 = buffer.data(ki + 388);
    const auto *ki_389 = buffer.data(ki + 389);
    const auto *ki_390 = buffer.data(ki + 390);
    const auto *ki_391 = buffer.data(ki + 391);
    const auto *ki_392 = buffer.data(ki + 392);
    const auto *ki_393 = buffer.data(ki + 393);
    const auto *ki_394 = buffer.data(ki + 394);
    const auto *ki_395 = buffer.data(ki + 395);
    const auto *ki_396 = buffer.data(ki + 396);
    const auto *ki_397 = buffer.data(ki + 397);
    const auto *ki_398 = buffer.data(ki + 398);
    const auto *ki_399 = buffer.data(ki + 399);
    const auto *ki_400 = buffer.data(ki + 400);
    const auto *ki_401 = buffer.data(ki + 401);
    const auto *ki_402 = buffer.data(ki + 402);
    const auto *ki_403 = buffer.data(ki + 403);
    const auto *ki_404 = buffer.data(ki + 404);
    const auto *ki_405 = buffer.data(ki + 405);
    const auto *ki_406 = buffer.data(ki + 406);
    const auto *ki_407 = buffer.data(ki + 407);
    const auto *ki_408 = buffer.data(ki + 408);
    const auto *ki_409 = buffer.data(ki + 409);
    const auto *ki_410 = buffer.data(ki + 410);
    const auto *ki_411 = buffer.data(ki + 411);
    const auto *ki_412 = buffer.data(ki + 412);
    const auto *ki_413 = buffer.data(ki + 413);
    const auto *ki_414 = buffer.data(ki + 414);
    const auto *ki_415 = buffer.data(ki + 415);
    const auto *ki_416 = buffer.data(ki + 416);
    const auto *ki_417 = buffer.data(ki + 417);
    const auto *ki_418 = buffer.data(ki + 418);
    const auto *ki_419 = buffer.data(ki + 419);
    const auto *ki_420 = buffer.data(ki + 420);
    const auto *ki_421 = buffer.data(ki + 421);
    const auto *ki_422 = buffer.data(ki + 422);
    const auto *ki_423 = buffer.data(ki + 423);
    const auto *ki_424 = buffer.data(ki + 424);
    const auto *ki_425 = buffer.data(ki + 425);
    const auto *ki_426 = buffer.data(ki + 426);
    const auto *ki_427 = buffer.data(ki + 427);
    const auto *ki_428 = buffer.data(ki + 428);
    const auto *ki_429 = buffer.data(ki + 429);
    const auto *ki_430 = buffer.data(ki + 430);
    const auto *ki_431 = buffer.data(ki + 431);
    const auto *ki_432 = buffer.data(ki + 432);
    const auto *ki_433 = buffer.data(ki + 433);
    const auto *ki_434 = buffer.data(ki + 434);
    const auto *ki_435 = buffer.data(ki + 435);
    const auto *ki_436 = buffer.data(ki + 436);
    const auto *ki_437 = buffer.data(ki + 437);
    const auto *ki_438 = buffer.data(ki + 438);
    const auto *ki_439 = buffer.data(ki + 439);
    const auto *ki_440 = buffer.data(ki + 440);
    const auto *ki_441 = buffer.data(ki + 441);
    const auto *ki_442 = buffer.data(ki + 442);
    const auto *ki_443 = buffer.data(ki + 443);
    const auto *ki_444 = buffer.data(ki + 444);
    const auto *ki_445 = buffer.data(ki + 445);
    const auto *ki_446 = buffer.data(ki + 446);
    const auto *ki_447 = buffer.data(ki + 447);
    const auto *ki_448 = buffer.data(ki + 448);
    const auto *ki_449 = buffer.data(ki + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, hi_300, hi_301, hi_302, hi_303, \
                         hi_304, ki_300, ki_301, ki_302, ki_303, \
                         ki_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -2.0 * hi_300[k]
                   + f_0 * ki_300[k];

        t_301[k] = -2.0 * hi_301[k]
                   + f_0 * ki_301[k];

        t_302[k] = -2.0 * hi_302[k]
                   + f_0 * ki_302[k];

        t_303[k] = -2.0 * hi_303[k]
                   + f_0 * ki_303[k];

        t_304[k] = -2.0 * hi_304[k]
                   + f_0 * ki_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, hi_305, hi_306, hi_307, hi_308, \
                         hi_309, ki_305, ki_306, ki_307, ki_308, \
                         ki_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -2.0 * hi_305[k]
                   + f_0 * ki_305[k];

        t_306[k] = -2.0 * hi_306[k]
                   + f_0 * ki_306[k];

        t_307[k] = -2.0 * hi_307[k]
                   + f_0 * ki_307[k];

        t_308[k] = -2.0 * hi_308[k]
                   + f_0 * ki_308[k];

        t_309[k] = -2.0 * hi_309[k]
                   + f_0 * ki_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, hi_310, hi_311, hi_312, hi_313, \
                         hi_314, ki_310, ki_311, ki_312, ki_313, \
                         ki_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -2.0 * hi_310[k]
                   + f_0 * ki_310[k];

        t_311[k] = -2.0 * hi_311[k]
                   + f_0 * ki_311[k];

        t_312[k] = -2.0 * hi_312[k]
                   + f_0 * ki_312[k];

        t_313[k] = -2.0 * hi_313[k]
                   + f_0 * ki_313[k];

        t_314[k] = -2.0 * hi_314[k]
                   + f_0 * ki_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, hi_315, hi_316, hi_317, hi_318, \
                         hi_319, ki_315, ki_316, ki_317, ki_318, \
                         ki_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -2.0 * hi_315[k]
                   + f_0 * ki_315[k];

        t_316[k] = -2.0 * hi_316[k]
                   + f_0 * ki_316[k];

        t_317[k] = -2.0 * hi_317[k]
                   + f_0 * ki_317[k];

        t_318[k] = -2.0 * hi_318[k]
                   + f_0 * ki_318[k];

        t_319[k] = -2.0 * hi_319[k]
                   + f_0 * ki_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, hi_320, hi_321, hi_322, hi_323, \
                         hi_324, ki_320, ki_321, ki_322, ki_323, \
                         ki_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -2.0 * hi_320[k]
                   + f_0 * ki_320[k];

        t_321[k] = -2.0 * hi_321[k]
                   + f_0 * ki_321[k];

        t_322[k] = -2.0 * hi_322[k]
                   + f_0 * ki_322[k];

        t_323[k] = -2.0 * hi_323[k]
                   + f_0 * ki_323[k];

        t_324[k] = -2.0 * hi_324[k]
                   + f_0 * ki_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, hi_325, hi_326, hi_327, hi_328, \
                         hi_329, ki_325, ki_326, ki_327, ki_328, \
                         ki_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -2.0 * hi_325[k]
                   + f_0 * ki_325[k];

        t_326[k] = -2.0 * hi_326[k]
                   + f_0 * ki_326[k];

        t_327[k] = -2.0 * hi_327[k]
                   + f_0 * ki_327[k];

        t_328[k] = -2.0 * hi_328[k]
                   + f_0 * ki_328[k];

        t_329[k] = -2.0 * hi_329[k]
                   + f_0 * ki_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, hi_330, hi_331, hi_332, hi_333, \
                         hi_334, ki_330, ki_331, ki_332, ki_333, \
                         ki_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -2.0 * hi_330[k]
                   + f_0 * ki_330[k];

        t_331[k] = -2.0 * hi_331[k]
                   + f_0 * ki_331[k];

        t_332[k] = -2.0 * hi_332[k]
                   + f_0 * ki_332[k];

        t_333[k] = -2.0 * hi_333[k]
                   + f_0 * ki_333[k];

        t_334[k] = -2.0 * hi_334[k]
                   + f_0 * ki_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, hi_335, hi_336, hi_337, hi_338, \
                         hi_339, ki_335, ki_336, ki_337, ki_338, \
                         ki_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -2.0 * hi_335[k]
                   + f_0 * ki_335[k];

        t_336[k] = -2.0 * hi_336[k]
                   + f_0 * ki_336[k];

        t_337[k] = -2.0 * hi_337[k]
                   + f_0 * ki_337[k];

        t_338[k] = -2.0 * hi_338[k]
                   + f_0 * ki_338[k];

        t_339[k] = -2.0 * hi_339[k]
                   + f_0 * ki_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, hi_340, hi_341, hi_342, hi_343, \
                         hi_344, ki_340, ki_341, ki_342, ki_343, \
                         ki_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * hi_340[k]
                   + f_0 * ki_340[k];

        t_341[k] = -2.0 * hi_341[k]
                   + f_0 * ki_341[k];

        t_342[k] = -2.0 * hi_342[k]
                   + f_0 * ki_342[k];

        t_343[k] = -2.0 * hi_343[k]
                   + f_0 * ki_343[k];

        t_344[k] = -2.0 * hi_344[k]
                   + f_0 * ki_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, hi_345, hi_346, hi_347, hi_348, \
                         hi_349, ki_345, ki_346, ki_347, ki_348, \
                         ki_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * hi_345[k]
                   + f_0 * ki_345[k];

        t_346[k] = -2.0 * hi_346[k]
                   + f_0 * ki_346[k];

        t_347[k] = -2.0 * hi_347[k]
                   + f_0 * ki_347[k];

        t_348[k] = -2.0 * hi_348[k]
                   + f_0 * ki_348[k];

        t_349[k] = -2.0 * hi_349[k]
                   + f_0 * ki_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, hi_350, hi_351, hi_352, hi_353, \
                         hi_354, ki_350, ki_351, ki_352, ki_353, \
                         ki_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * hi_350[k]
                   + f_0 * ki_350[k];

        t_351[k] = -2.0 * hi_351[k]
                   + f_0 * ki_351[k];

        t_352[k] = -2.0 * hi_352[k]
                   + f_0 * ki_352[k];

        t_353[k] = -2.0 * hi_353[k]
                   + f_0 * ki_353[k];

        t_354[k] = -2.0 * hi_354[k]
                   + f_0 * ki_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, hi_355, hi_356, hi_357, hi_358, \
                         hi_359, ki_355, ki_356, ki_357, ki_358, \
                         ki_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * hi_355[k]
                   + f_0 * ki_355[k];

        t_356[k] = -2.0 * hi_356[k]
                   + f_0 * ki_356[k];

        t_357[k] = -2.0 * hi_357[k]
                   + f_0 * ki_357[k];

        t_358[k] = -2.0 * hi_358[k]
                   + f_0 * ki_358[k];

        t_359[k] = -2.0 * hi_359[k]
                   + f_0 * ki_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, hi_360, hi_361, hi_362, hi_363, \
                         hi_364, ki_360, ki_361, ki_362, ki_363, \
                         ki_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -2.0 * hi_360[k]
                   + f_0 * ki_360[k];

        t_361[k] = -2.0 * hi_361[k]
                   + f_0 * ki_361[k];

        t_362[k] = -2.0 * hi_362[k]
                   + f_0 * ki_362[k];

        t_363[k] = -2.0 * hi_363[k]
                   + f_0 * ki_363[k];

        t_364[k] = -2.0 * hi_364[k]
                   + f_0 * ki_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, hi_365, hi_366, hi_367, hi_368, \
                         hi_369, ki_365, ki_366, ki_367, ki_368, \
                         ki_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -2.0 * hi_365[k]
                   + f_0 * ki_365[k];

        t_366[k] = -2.0 * hi_366[k]
                   + f_0 * ki_366[k];

        t_367[k] = -2.0 * hi_367[k]
                   + f_0 * ki_367[k];

        t_368[k] = -2.0 * hi_368[k]
                   + f_0 * ki_368[k];

        t_369[k] = -2.0 * hi_369[k]
                   + f_0 * ki_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, hi_370, hi_371, hi_372, hi_373, \
                         hi_374, ki_370, ki_371, ki_372, ki_373, \
                         ki_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -2.0 * hi_370[k]
                   + f_0 * ki_370[k];

        t_371[k] = -2.0 * hi_371[k]
                   + f_0 * ki_371[k];

        t_372[k] = -2.0 * hi_372[k]
                   + f_0 * ki_372[k];

        t_373[k] = -2.0 * hi_373[k]
                   + f_0 * ki_373[k];

        t_374[k] = -2.0 * hi_374[k]
                   + f_0 * ki_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, hi_375, hi_376, hi_377, hi_378, \
                         hi_379, ki_375, ki_376, ki_377, ki_378, \
                         ki_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -2.0 * hi_375[k]
                   + f_0 * ki_375[k];

        t_376[k] = -2.0 * hi_376[k]
                   + f_0 * ki_376[k];

        t_377[k] = -2.0 * hi_377[k]
                   + f_0 * ki_377[k];

        t_378[k] = -2.0 * hi_378[k]
                   + f_0 * ki_378[k];

        t_379[k] = -2.0 * hi_379[k]
                   + f_0 * ki_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, hi_380, hi_381, hi_382, hi_383, \
                         hi_384, ki_380, ki_381, ki_382, ki_383, \
                         ki_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * hi_380[k]
                   + f_0 * ki_380[k];

        t_381[k] = -2.0 * hi_381[k]
                   + f_0 * ki_381[k];

        t_382[k] = -2.0 * hi_382[k]
                   + f_0 * ki_382[k];

        t_383[k] = -2.0 * hi_383[k]
                   + f_0 * ki_383[k];

        t_384[k] = -2.0 * hi_384[k]
                   + f_0 * ki_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, hi_385, hi_386, hi_387, hi_388, \
                         hi_389, ki_385, ki_386, ki_387, ki_388, \
                         ki_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * hi_385[k]
                   + f_0 * ki_385[k];

        t_386[k] = -2.0 * hi_386[k]
                   + f_0 * ki_386[k];

        t_387[k] = -2.0 * hi_387[k]
                   + f_0 * ki_387[k];

        t_388[k] = -2.0 * hi_388[k]
                   + f_0 * ki_388[k];

        t_389[k] = -2.0 * hi_389[k]
                   + f_0 * ki_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, hi_390, hi_391, hi_392, hi_393, \
                         hi_394, ki_390, ki_391, ki_392, ki_393, \
                         ki_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -2.0 * hi_390[k]
                   + f_0 * ki_390[k];

        t_391[k] = -2.0 * hi_391[k]
                   + f_0 * ki_391[k];

        t_392[k] = -2.0 * hi_392[k]
                   + f_0 * ki_392[k];

        t_393[k] = -2.0 * hi_393[k]
                   + f_0 * ki_393[k];

        t_394[k] = -2.0 * hi_394[k]
                   + f_0 * ki_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, hi_395, hi_396, hi_397, hi_398, \
                         hi_399, ki_395, ki_396, ki_397, ki_398, \
                         ki_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -2.0 * hi_395[k]
                   + f_0 * ki_395[k];

        t_396[k] = -2.0 * hi_396[k]
                   + f_0 * ki_396[k];

        t_397[k] = -2.0 * hi_397[k]
                   + f_0 * ki_397[k];

        t_398[k] = -2.0 * hi_398[k]
                   + f_0 * ki_398[k];

        t_399[k] = -2.0 * hi_399[k]
                   + f_0 * ki_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, hi_400, hi_401, hi_402, hi_403, \
                         hi_404, ki_400, ki_401, ki_402, ki_403, \
                         ki_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -2.0 * hi_400[k]
                   + f_0 * ki_400[k];

        t_401[k] = -2.0 * hi_401[k]
                   + f_0 * ki_401[k];

        t_402[k] = -2.0 * hi_402[k]
                   + f_0 * ki_402[k];

        t_403[k] = -2.0 * hi_403[k]
                   + f_0 * ki_403[k];

        t_404[k] = -2.0 * hi_404[k]
                   + f_0 * ki_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, hi_405, hi_406, hi_407, hi_408, \
                         hi_409, ki_405, ki_406, ki_407, ki_408, \
                         ki_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -2.0 * hi_405[k]
                   + f_0 * ki_405[k];

        t_406[k] = -2.0 * hi_406[k]
                   + f_0 * ki_406[k];

        t_407[k] = -2.0 * hi_407[k]
                   + f_0 * ki_407[k];

        t_408[k] = -2.0 * hi_408[k]
                   + f_0 * ki_408[k];

        t_409[k] = -2.0 * hi_409[k]
                   + f_0 * ki_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, hi_410, hi_411, hi_412, hi_413, \
                         hi_414, ki_410, ki_411, ki_412, ki_413, \
                         ki_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -2.0 * hi_410[k]
                   + f_0 * ki_410[k];

        t_411[k] = -2.0 * hi_411[k]
                   + f_0 * ki_411[k];

        t_412[k] = -2.0 * hi_412[k]
                   + f_0 * ki_412[k];

        t_413[k] = -2.0 * hi_413[k]
                   + f_0 * ki_413[k];

        t_414[k] = -2.0 * hi_414[k]
                   + f_0 * ki_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, hi_415, hi_416, hi_417, hi_418, \
                         hi_419, ki_415, ki_416, ki_417, ki_418, \
                         ki_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -2.0 * hi_415[k]
                   + f_0 * ki_415[k];

        t_416[k] = -2.0 * hi_416[k]
                   + f_0 * ki_416[k];

        t_417[k] = -2.0 * hi_417[k]
                   + f_0 * ki_417[k];

        t_418[k] = -2.0 * hi_418[k]
                   + f_0 * ki_418[k];

        t_419[k] = -2.0 * hi_419[k]
                   + f_0 * ki_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, hi_420, hi_421, hi_422, hi_423, \
                         hi_424, ki_420, ki_421, ki_422, ki_423, \
                         ki_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -hi_420[k]
                   + f_0 * ki_420[k];

        t_421[k] = -hi_421[k]
                   + f_0 * ki_421[k];

        t_422[k] = -hi_422[k]
                   + f_0 * ki_422[k];

        t_423[k] = -hi_423[k]
                   + f_0 * ki_423[k];

        t_424[k] = -hi_424[k]
                   + f_0 * ki_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, hi_425, hi_426, hi_427, hi_428, \
                         hi_429, ki_425, ki_426, ki_427, ki_428, \
                         ki_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -hi_425[k]
                   + f_0 * ki_425[k];

        t_426[k] = -hi_426[k]
                   + f_0 * ki_426[k];

        t_427[k] = -hi_427[k]
                   + f_0 * ki_427[k];

        t_428[k] = -hi_428[k]
                   + f_0 * ki_428[k];

        t_429[k] = -hi_429[k]
                   + f_0 * ki_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, hi_430, hi_431, hi_432, hi_433, \
                         hi_434, ki_430, ki_431, ki_432, ki_433, \
                         ki_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -hi_430[k]
                   + f_0 * ki_430[k];

        t_431[k] = -hi_431[k]
                   + f_0 * ki_431[k];

        t_432[k] = -hi_432[k]
                   + f_0 * ki_432[k];

        t_433[k] = -hi_433[k]
                   + f_0 * ki_433[k];

        t_434[k] = -hi_434[k]
                   + f_0 * ki_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, hi_435, hi_436, hi_437, hi_438, \
                         hi_439, ki_435, ki_436, ki_437, ki_438, \
                         ki_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -hi_435[k]
                   + f_0 * ki_435[k];

        t_436[k] = -hi_436[k]
                   + f_0 * ki_436[k];

        t_437[k] = -hi_437[k]
                   + f_0 * ki_437[k];

        t_438[k] = -hi_438[k]
                   + f_0 * ki_438[k];

        t_439[k] = -hi_439[k]
                   + f_0 * ki_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, hi_440, hi_441, hi_442, hi_443, \
                         hi_444, ki_440, ki_441, ki_442, ki_443, \
                         ki_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -hi_440[k]
                   + f_0 * ki_440[k];

        t_441[k] = -hi_441[k]
                   + f_0 * ki_441[k];

        t_442[k] = -hi_442[k]
                   + f_0 * ki_442[k];

        t_443[k] = -hi_443[k]
                   + f_0 * ki_443[k];

        t_444[k] = -hi_444[k]
                   + f_0 * ki_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, hi_445, hi_446, hi_447, hi_448, \
                         hi_449, ki_445, ki_446, ki_447, ki_448, \
                         ki_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -hi_445[k]
                   + f_0 * ki_445[k];

        t_446[k] = -hi_446[k]
                   + f_0 * ki_446[k];

        t_447[k] = -hi_447[k]
                   + f_0 * ki_447[k];

        t_448[k] = -hi_448[k]
                   + f_0 * ki_448[k];

        t_449[k] = -hi_449[k]
                   + f_0 * ki_449[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
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

    const auto *hi_450 = buffer.data(hi + 450);
    const auto *hi_451 = buffer.data(hi + 451);
    const auto *hi_452 = buffer.data(hi + 452);
    const auto *hi_453 = buffer.data(hi + 453);
    const auto *hi_454 = buffer.data(hi + 454);
    const auto *hi_455 = buffer.data(hi + 455);
    const auto *hi_456 = buffer.data(hi + 456);
    const auto *hi_457 = buffer.data(hi + 457);
    const auto *hi_458 = buffer.data(hi + 458);
    const auto *hi_459 = buffer.data(hi + 459);
    const auto *hi_460 = buffer.data(hi + 460);
    const auto *hi_461 = buffer.data(hi + 461);
    const auto *hi_462 = buffer.data(hi + 462);
    const auto *hi_463 = buffer.data(hi + 463);
    const auto *hi_464 = buffer.data(hi + 464);
    const auto *hi_465 = buffer.data(hi + 465);
    const auto *hi_466 = buffer.data(hi + 466);
    const auto *hi_467 = buffer.data(hi + 467);
    const auto *hi_468 = buffer.data(hi + 468);
    const auto *hi_469 = buffer.data(hi + 469);
    const auto *hi_470 = buffer.data(hi + 470);
    const auto *hi_471 = buffer.data(hi + 471);
    const auto *hi_472 = buffer.data(hi + 472);
    const auto *hi_473 = buffer.data(hi + 473);
    const auto *hi_474 = buffer.data(hi + 474);
    const auto *hi_475 = buffer.data(hi + 475);
    const auto *hi_476 = buffer.data(hi + 476);
    const auto *hi_477 = buffer.data(hi + 477);
    const auto *hi_478 = buffer.data(hi + 478);
    const auto *hi_479 = buffer.data(hi + 479);
    const auto *hi_480 = buffer.data(hi + 480);
    const auto *hi_481 = buffer.data(hi + 481);
    const auto *hi_482 = buffer.data(hi + 482);
    const auto *hi_483 = buffer.data(hi + 483);
    const auto *hi_484 = buffer.data(hi + 484);
    const auto *hi_485 = buffer.data(hi + 485);
    const auto *hi_486 = buffer.data(hi + 486);
    const auto *hi_487 = buffer.data(hi + 487);
    const auto *hi_488 = buffer.data(hi + 488);
    const auto *hi_489 = buffer.data(hi + 489);
    const auto *hi_490 = buffer.data(hi + 490);
    const auto *hi_491 = buffer.data(hi + 491);
    const auto *hi_492 = buffer.data(hi + 492);
    const auto *hi_493 = buffer.data(hi + 493);
    const auto *hi_494 = buffer.data(hi + 494);
    const auto *hi_495 = buffer.data(hi + 495);
    const auto *hi_496 = buffer.data(hi + 496);
    const auto *hi_497 = buffer.data(hi + 497);
    const auto *hi_498 = buffer.data(hi + 498);
    const auto *hi_499 = buffer.data(hi + 499);
    const auto *hi_500 = buffer.data(hi + 500);
    const auto *hi_501 = buffer.data(hi + 501);
    const auto *hi_502 = buffer.data(hi + 502);
    const auto *hi_503 = buffer.data(hi + 503);
    const auto *hi_504 = buffer.data(hi + 504);
    const auto *hi_505 = buffer.data(hi + 505);
    const auto *hi_506 = buffer.data(hi + 506);
    const auto *hi_507 = buffer.data(hi + 507);
    const auto *hi_508 = buffer.data(hi + 508);
    const auto *hi_509 = buffer.data(hi + 509);
    const auto *hi_510 = buffer.data(hi + 510);
    const auto *hi_511 = buffer.data(hi + 511);
    const auto *hi_512 = buffer.data(hi + 512);
    const auto *hi_513 = buffer.data(hi + 513);
    const auto *hi_514 = buffer.data(hi + 514);
    const auto *hi_515 = buffer.data(hi + 515);
    const auto *hi_516 = buffer.data(hi + 516);
    const auto *hi_517 = buffer.data(hi + 517);
    const auto *hi_518 = buffer.data(hi + 518);
    const auto *hi_519 = buffer.data(hi + 519);
    const auto *hi_520 = buffer.data(hi + 520);
    const auto *hi_521 = buffer.data(hi + 521);
    const auto *hi_522 = buffer.data(hi + 522);
    const auto *hi_523 = buffer.data(hi + 523);
    const auto *hi_524 = buffer.data(hi + 524);
    const auto *hi_525 = buffer.data(hi + 525);
    const auto *hi_526 = buffer.data(hi + 526);
    const auto *hi_527 = buffer.data(hi + 527);
    const auto *hi_528 = buffer.data(hi + 528);
    const auto *hi_529 = buffer.data(hi + 529);
    const auto *hi_530 = buffer.data(hi + 530);
    const auto *hi_531 = buffer.data(hi + 531);
    const auto *hi_532 = buffer.data(hi + 532);
    const auto *hi_533 = buffer.data(hi + 533);
    const auto *hi_534 = buffer.data(hi + 534);
    const auto *hi_535 = buffer.data(hi + 535);
    const auto *hi_536 = buffer.data(hi + 536);
    const auto *hi_537 = buffer.data(hi + 537);
    const auto *hi_538 = buffer.data(hi + 538);
    const auto *hi_539 = buffer.data(hi + 539);
    const auto *hi_540 = buffer.data(hi + 540);
    const auto *hi_541 = buffer.data(hi + 541);
    const auto *hi_542 = buffer.data(hi + 542);
    const auto *hi_543 = buffer.data(hi + 543);
    const auto *hi_544 = buffer.data(hi + 544);
    const auto *hi_545 = buffer.data(hi + 545);
    const auto *hi_546 = buffer.data(hi + 546);
    const auto *hi_547 = buffer.data(hi + 547);
    const auto *hi_548 = buffer.data(hi + 548);
    const auto *hi_549 = buffer.data(hi + 549);
    const auto *hi_550 = buffer.data(hi + 550);
    const auto *hi_551 = buffer.data(hi + 551);
    const auto *hi_552 = buffer.data(hi + 552);
    const auto *hi_553 = buffer.data(hi + 553);
    const auto *hi_554 = buffer.data(hi + 554);
    const auto *hi_555 = buffer.data(hi + 555);
    const auto *hi_556 = buffer.data(hi + 556);
    const auto *hi_557 = buffer.data(hi + 557);
    const auto *hi_558 = buffer.data(hi + 558);
    const auto *hi_559 = buffer.data(hi + 559);
    const auto *hi_560 = buffer.data(hi + 560);
    const auto *hi_561 = buffer.data(hi + 561);
    const auto *hi_562 = buffer.data(hi + 562);
    const auto *hi_563 = buffer.data(hi + 563);
    const auto *hi_564 = buffer.data(hi + 564);
    const auto *hi_565 = buffer.data(hi + 565);
    const auto *hi_566 = buffer.data(hi + 566);
    const auto *hi_567 = buffer.data(hi + 567);
    const auto *hi_568 = buffer.data(hi + 568);
    const auto *hi_569 = buffer.data(hi + 569);
    const auto *hi_570 = buffer.data(hi + 570);
    const auto *hi_571 = buffer.data(hi + 571);
    const auto *hi_572 = buffer.data(hi + 572);
    const auto *hi_573 = buffer.data(hi + 573);
    const auto *hi_574 = buffer.data(hi + 574);
    const auto *hi_575 = buffer.data(hi + 575);
    const auto *hi_576 = buffer.data(hi + 576);
    const auto *hi_577 = buffer.data(hi + 577);
    const auto *hi_578 = buffer.data(hi + 578);
    const auto *hi_579 = buffer.data(hi + 579);
    const auto *hi_580 = buffer.data(hi + 580);
    const auto *hi_581 = buffer.data(hi + 581);
    const auto *hi_582 = buffer.data(hi + 582);
    const auto *hi_583 = buffer.data(hi + 583);
    const auto *hi_584 = buffer.data(hi + 584);
    const auto *hi_585 = buffer.data(hi + 585);
    const auto *hi_586 = buffer.data(hi + 586);
    const auto *hi_587 = buffer.data(hi + 587);

    const auto *ki_450 = buffer.data(ki + 450);
    const auto *ki_451 = buffer.data(ki + 451);
    const auto *ki_452 = buffer.data(ki + 452);
    const auto *ki_453 = buffer.data(ki + 453);
    const auto *ki_454 = buffer.data(ki + 454);
    const auto *ki_455 = buffer.data(ki + 455);
    const auto *ki_456 = buffer.data(ki + 456);
    const auto *ki_457 = buffer.data(ki + 457);
    const auto *ki_458 = buffer.data(ki + 458);
    const auto *ki_459 = buffer.data(ki + 459);
    const auto *ki_460 = buffer.data(ki + 460);
    const auto *ki_461 = buffer.data(ki + 461);
    const auto *ki_462 = buffer.data(ki + 462);
    const auto *ki_463 = buffer.data(ki + 463);
    const auto *ki_464 = buffer.data(ki + 464);
    const auto *ki_465 = buffer.data(ki + 465);
    const auto *ki_466 = buffer.data(ki + 466);
    const auto *ki_467 = buffer.data(ki + 467);
    const auto *ki_468 = buffer.data(ki + 468);
    const auto *ki_469 = buffer.data(ki + 469);
    const auto *ki_470 = buffer.data(ki + 470);
    const auto *ki_471 = buffer.data(ki + 471);
    const auto *ki_472 = buffer.data(ki + 472);
    const auto *ki_473 = buffer.data(ki + 473);
    const auto *ki_474 = buffer.data(ki + 474);
    const auto *ki_475 = buffer.data(ki + 475);
    const auto *ki_476 = buffer.data(ki + 476);
    const auto *ki_477 = buffer.data(ki + 477);
    const auto *ki_478 = buffer.data(ki + 478);
    const auto *ki_479 = buffer.data(ki + 479);
    const auto *ki_480 = buffer.data(ki + 480);
    const auto *ki_481 = buffer.data(ki + 481);
    const auto *ki_482 = buffer.data(ki + 482);
    const auto *ki_483 = buffer.data(ki + 483);
    const auto *ki_484 = buffer.data(ki + 484);
    const auto *ki_485 = buffer.data(ki + 485);
    const auto *ki_486 = buffer.data(ki + 486);
    const auto *ki_487 = buffer.data(ki + 487);
    const auto *ki_488 = buffer.data(ki + 488);
    const auto *ki_489 = buffer.data(ki + 489);
    const auto *ki_490 = buffer.data(ki + 490);
    const auto *ki_491 = buffer.data(ki + 491);
    const auto *ki_492 = buffer.data(ki + 492);
    const auto *ki_493 = buffer.data(ki + 493);
    const auto *ki_494 = buffer.data(ki + 494);
    const auto *ki_495 = buffer.data(ki + 495);
    const auto *ki_496 = buffer.data(ki + 496);
    const auto *ki_497 = buffer.data(ki + 497);
    const auto *ki_498 = buffer.data(ki + 498);
    const auto *ki_499 = buffer.data(ki + 499);
    const auto *ki_500 = buffer.data(ki + 500);
    const auto *ki_501 = buffer.data(ki + 501);
    const auto *ki_502 = buffer.data(ki + 502);
    const auto *ki_503 = buffer.data(ki + 503);
    const auto *ki_504 = buffer.data(ki + 504);
    const auto *ki_505 = buffer.data(ki + 505);
    const auto *ki_506 = buffer.data(ki + 506);
    const auto *ki_507 = buffer.data(ki + 507);
    const auto *ki_508 = buffer.data(ki + 508);
    const auto *ki_509 = buffer.data(ki + 509);
    const auto *ki_510 = buffer.data(ki + 510);
    const auto *ki_511 = buffer.data(ki + 511);
    const auto *ki_512 = buffer.data(ki + 512);
    const auto *ki_513 = buffer.data(ki + 513);
    const auto *ki_514 = buffer.data(ki + 514);
    const auto *ki_515 = buffer.data(ki + 515);
    const auto *ki_516 = buffer.data(ki + 516);
    const auto *ki_517 = buffer.data(ki + 517);
    const auto *ki_518 = buffer.data(ki + 518);
    const auto *ki_519 = buffer.data(ki + 519);
    const auto *ki_520 = buffer.data(ki + 520);
    const auto *ki_521 = buffer.data(ki + 521);
    const auto *ki_522 = buffer.data(ki + 522);
    const auto *ki_523 = buffer.data(ki + 523);
    const auto *ki_524 = buffer.data(ki + 524);
    const auto *ki_525 = buffer.data(ki + 525);
    const auto *ki_526 = buffer.data(ki + 526);
    const auto *ki_527 = buffer.data(ki + 527);
    const auto *ki_528 = buffer.data(ki + 528);
    const auto *ki_529 = buffer.data(ki + 529);
    const auto *ki_530 = buffer.data(ki + 530);
    const auto *ki_531 = buffer.data(ki + 531);
    const auto *ki_532 = buffer.data(ki + 532);
    const auto *ki_533 = buffer.data(ki + 533);
    const auto *ki_534 = buffer.data(ki + 534);
    const auto *ki_535 = buffer.data(ki + 535);
    const auto *ki_536 = buffer.data(ki + 536);
    const auto *ki_537 = buffer.data(ki + 537);
    const auto *ki_538 = buffer.data(ki + 538);
    const auto *ki_539 = buffer.data(ki + 539);
    const auto *ki_540 = buffer.data(ki + 540);
    const auto *ki_541 = buffer.data(ki + 541);
    const auto *ki_542 = buffer.data(ki + 542);
    const auto *ki_543 = buffer.data(ki + 543);
    const auto *ki_544 = buffer.data(ki + 544);
    const auto *ki_545 = buffer.data(ki + 545);
    const auto *ki_546 = buffer.data(ki + 546);
    const auto *ki_547 = buffer.data(ki + 547);
    const auto *ki_548 = buffer.data(ki + 548);
    const auto *ki_549 = buffer.data(ki + 549);
    const auto *ki_550 = buffer.data(ki + 550);
    const auto *ki_551 = buffer.data(ki + 551);
    const auto *ki_552 = buffer.data(ki + 552);
    const auto *ki_553 = buffer.data(ki + 553);
    const auto *ki_554 = buffer.data(ki + 554);
    const auto *ki_555 = buffer.data(ki + 555);
    const auto *ki_556 = buffer.data(ki + 556);
    const auto *ki_557 = buffer.data(ki + 557);
    const auto *ki_558 = buffer.data(ki + 558);
    const auto *ki_559 = buffer.data(ki + 559);
    const auto *ki_560 = buffer.data(ki + 560);
    const auto *ki_561 = buffer.data(ki + 561);
    const auto *ki_562 = buffer.data(ki + 562);
    const auto *ki_563 = buffer.data(ki + 563);
    const auto *ki_564 = buffer.data(ki + 564);
    const auto *ki_565 = buffer.data(ki + 565);
    const auto *ki_566 = buffer.data(ki + 566);
    const auto *ki_567 = buffer.data(ki + 567);
    const auto *ki_568 = buffer.data(ki + 568);
    const auto *ki_569 = buffer.data(ki + 569);
    const auto *ki_570 = buffer.data(ki + 570);
    const auto *ki_571 = buffer.data(ki + 571);
    const auto *ki_572 = buffer.data(ki + 572);
    const auto *ki_573 = buffer.data(ki + 573);
    const auto *ki_574 = buffer.data(ki + 574);
    const auto *ki_575 = buffer.data(ki + 575);
    const auto *ki_576 = buffer.data(ki + 576);
    const auto *ki_577 = buffer.data(ki + 577);
    const auto *ki_578 = buffer.data(ki + 578);
    const auto *ki_579 = buffer.data(ki + 579);
    const auto *ki_580 = buffer.data(ki + 580);
    const auto *ki_581 = buffer.data(ki + 581);
    const auto *ki_582 = buffer.data(ki + 582);
    const auto *ki_583 = buffer.data(ki + 583);
    const auto *ki_584 = buffer.data(ki + 584);
    const auto *ki_585 = buffer.data(ki + 585);
    const auto *ki_586 = buffer.data(ki + 586);
    const auto *ki_587 = buffer.data(ki + 587);
    const auto *ki_588 = buffer.data(ki + 588);
    const auto *ki_589 = buffer.data(ki + 589);
    const auto *ki_590 = buffer.data(ki + 590);
    const auto *ki_591 = buffer.data(ki + 591);
    const auto *ki_592 = buffer.data(ki + 592);
    const auto *ki_593 = buffer.data(ki + 593);
    const auto *ki_594 = buffer.data(ki + 594);
    const auto *ki_595 = buffer.data(ki + 595);
    const auto *ki_596 = buffer.data(ki + 596);
    const auto *ki_597 = buffer.data(ki + 597);
    const auto *ki_598 = buffer.data(ki + 598);
    const auto *ki_599 = buffer.data(ki + 599);
    const auto *ki_600 = buffer.data(ki + 600);
    const auto *ki_601 = buffer.data(ki + 601);
    const auto *ki_602 = buffer.data(ki + 602);
    const auto *ki_603 = buffer.data(ki + 603);
    const auto *ki_604 = buffer.data(ki + 604);
    const auto *ki_605 = buffer.data(ki + 605);
    const auto *ki_606 = buffer.data(ki + 606);

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, hi_450, hi_451, hi_452, hi_453, \
                         hi_454, ki_450, ki_451, ki_452, ki_453, \
                         ki_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -hi_450[k]
                   + f_0 * ki_450[k];

        t_451[k] = -hi_451[k]
                   + f_0 * ki_451[k];

        t_452[k] = -hi_452[k]
                   + f_0 * ki_452[k];

        t_453[k] = -hi_453[k]
                   + f_0 * ki_453[k];

        t_454[k] = -hi_454[k]
                   + f_0 * ki_454[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, hi_455, hi_456, hi_457, hi_458, \
                         hi_459, ki_455, ki_456, ki_457, ki_458, \
                         ki_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -hi_455[k]
                   + f_0 * ki_455[k];

        t_456[k] = -hi_456[k]
                   + f_0 * ki_456[k];

        t_457[k] = -hi_457[k]
                   + f_0 * ki_457[k];

        t_458[k] = -hi_458[k]
                   + f_0 * ki_458[k];

        t_459[k] = -hi_459[k]
                   + f_0 * ki_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, hi_460, hi_461, hi_462, hi_463, \
                         hi_464, ki_460, ki_461, ki_462, ki_463, \
                         ki_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -hi_460[k]
                   + f_0 * ki_460[k];

        t_461[k] = -hi_461[k]
                   + f_0 * ki_461[k];

        t_462[k] = -hi_462[k]
                   + f_0 * ki_462[k];

        t_463[k] = -hi_463[k]
                   + f_0 * ki_463[k];

        t_464[k] = -hi_464[k]
                   + f_0 * ki_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, hi_465, hi_466, hi_467, hi_468, \
                         hi_469, ki_465, ki_466, ki_467, ki_468, \
                         ki_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -hi_465[k]
                   + f_0 * ki_465[k];

        t_466[k] = -hi_466[k]
                   + f_0 * ki_466[k];

        t_467[k] = -hi_467[k]
                   + f_0 * ki_467[k];

        t_468[k] = -hi_468[k]
                   + f_0 * ki_468[k];

        t_469[k] = -hi_469[k]
                   + f_0 * ki_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, hi_470, hi_471, hi_472, hi_473, \
                         hi_474, ki_470, ki_471, ki_472, ki_473, \
                         ki_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -hi_470[k]
                   + f_0 * ki_470[k];

        t_471[k] = -hi_471[k]
                   + f_0 * ki_471[k];

        t_472[k] = -hi_472[k]
                   + f_0 * ki_472[k];

        t_473[k] = -hi_473[k]
                   + f_0 * ki_473[k];

        t_474[k] = -hi_474[k]
                   + f_0 * ki_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, hi_475, hi_476, hi_477, hi_478, \
                         hi_479, ki_475, ki_476, ki_477, ki_478, \
                         ki_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -hi_475[k]
                   + f_0 * ki_475[k];

        t_476[k] = -hi_476[k]
                   + f_0 * ki_476[k];

        t_477[k] = -hi_477[k]
                   + f_0 * ki_477[k];

        t_478[k] = -hi_478[k]
                   + f_0 * ki_478[k];

        t_479[k] = -hi_479[k]
                   + f_0 * ki_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, hi_480, hi_481, hi_482, hi_483, \
                         hi_484, ki_480, ki_481, ki_482, ki_483, \
                         ki_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -hi_480[k]
                   + f_0 * ki_480[k];

        t_481[k] = -hi_481[k]
                   + f_0 * ki_481[k];

        t_482[k] = -hi_482[k]
                   + f_0 * ki_482[k];

        t_483[k] = -hi_483[k]
                   + f_0 * ki_483[k];

        t_484[k] = -hi_484[k]
                   + f_0 * ki_484[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, hi_485, hi_486, hi_487, hi_488, \
                         hi_489, ki_485, ki_486, ki_487, ki_488, \
                         ki_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -hi_485[k]
                   + f_0 * ki_485[k];

        t_486[k] = -hi_486[k]
                   + f_0 * ki_486[k];

        t_487[k] = -hi_487[k]
                   + f_0 * ki_487[k];

        t_488[k] = -hi_488[k]
                   + f_0 * ki_488[k];

        t_489[k] = -hi_489[k]
                   + f_0 * ki_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, hi_490, hi_491, hi_492, hi_493, \
                         hi_494, ki_490, ki_491, ki_492, ki_493, \
                         ki_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -hi_490[k]
                   + f_0 * ki_490[k];

        t_491[k] = -hi_491[k]
                   + f_0 * ki_491[k];

        t_492[k] = -hi_492[k]
                   + f_0 * ki_492[k];

        t_493[k] = -hi_493[k]
                   + f_0 * ki_493[k];

        t_494[k] = -hi_494[k]
                   + f_0 * ki_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, hi_495, hi_496, hi_497, hi_498, \
                         hi_499, ki_495, ki_496, ki_497, ki_498, \
                         ki_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -hi_495[k]
                   + f_0 * ki_495[k];

        t_496[k] = -hi_496[k]
                   + f_0 * ki_496[k];

        t_497[k] = -hi_497[k]
                   + f_0 * ki_497[k];

        t_498[k] = -hi_498[k]
                   + f_0 * ki_498[k];

        t_499[k] = -hi_499[k]
                   + f_0 * ki_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, hi_500, hi_501, hi_502, hi_503, \
                         hi_504, ki_500, ki_501, ki_502, ki_503, \
                         ki_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -hi_500[k]
                   + f_0 * ki_500[k];

        t_501[k] = -hi_501[k]
                   + f_0 * ki_501[k];

        t_502[k] = -hi_502[k]
                   + f_0 * ki_502[k];

        t_503[k] = -hi_503[k]
                   + f_0 * ki_503[k];

        t_504[k] = -hi_504[k]
                   + f_0 * ki_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, hi_505, hi_506, hi_507, hi_508, \
                         hi_509, ki_505, ki_506, ki_507, ki_508, \
                         ki_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -hi_505[k]
                   + f_0 * ki_505[k];

        t_506[k] = -hi_506[k]
                   + f_0 * ki_506[k];

        t_507[k] = -hi_507[k]
                   + f_0 * ki_507[k];

        t_508[k] = -hi_508[k]
                   + f_0 * ki_508[k];

        t_509[k] = -hi_509[k]
                   + f_0 * ki_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, hi_510, hi_511, hi_512, hi_513, \
                         hi_514, ki_510, ki_511, ki_512, ki_513, \
                         ki_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -hi_510[k]
                   + f_0 * ki_510[k];

        t_511[k] = -hi_511[k]
                   + f_0 * ki_511[k];

        t_512[k] = -hi_512[k]
                   + f_0 * ki_512[k];

        t_513[k] = -hi_513[k]
                   + f_0 * ki_513[k];

        t_514[k] = -hi_514[k]
                   + f_0 * ki_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, hi_515, hi_516, hi_517, hi_518, \
                         hi_519, ki_515, ki_516, ki_517, ki_518, \
                         ki_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -hi_515[k]
                   + f_0 * ki_515[k];

        t_516[k] = -hi_516[k]
                   + f_0 * ki_516[k];

        t_517[k] = -hi_517[k]
                   + f_0 * ki_517[k];

        t_518[k] = -hi_518[k]
                   + f_0 * ki_518[k];

        t_519[k] = -hi_519[k]
                   + f_0 * ki_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, hi_520, hi_521, hi_522, hi_523, \
                         hi_524, ki_520, ki_521, ki_522, ki_523, \
                         ki_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -hi_520[k]
                   + f_0 * ki_520[k];

        t_521[k] = -hi_521[k]
                   + f_0 * ki_521[k];

        t_522[k] = -hi_522[k]
                   + f_0 * ki_522[k];

        t_523[k] = -hi_523[k]
                   + f_0 * ki_523[k];

        t_524[k] = -hi_524[k]
                   + f_0 * ki_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, hi_525, hi_526, hi_527, hi_528, \
                         hi_529, ki_525, ki_526, ki_527, ki_528, \
                         ki_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -hi_525[k]
                   + f_0 * ki_525[k];

        t_526[k] = -hi_526[k]
                   + f_0 * ki_526[k];

        t_527[k] = -hi_527[k]
                   + f_0 * ki_527[k];

        t_528[k] = -hi_528[k]
                   + f_0 * ki_528[k];

        t_529[k] = -hi_529[k]
                   + f_0 * ki_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, hi_530, hi_531, hi_532, hi_533, \
                         hi_534, ki_530, ki_531, ki_532, ki_533, \
                         ki_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -hi_530[k]
                   + f_0 * ki_530[k];

        t_531[k] = -hi_531[k]
                   + f_0 * ki_531[k];

        t_532[k] = -hi_532[k]
                   + f_0 * ki_532[k];

        t_533[k] = -hi_533[k]
                   + f_0 * ki_533[k];

        t_534[k] = -hi_534[k]
                   + f_0 * ki_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, hi_535, hi_536, hi_537, hi_538, \
                         hi_539, ki_535, ki_536, ki_537, ki_538, \
                         ki_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -hi_535[k]
                   + f_0 * ki_535[k];

        t_536[k] = -hi_536[k]
                   + f_0 * ki_536[k];

        t_537[k] = -hi_537[k]
                   + f_0 * ki_537[k];

        t_538[k] = -hi_538[k]
                   + f_0 * ki_538[k];

        t_539[k] = -hi_539[k]
                   + f_0 * ki_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, hi_540, hi_541, hi_542, hi_543, \
                         hi_544, ki_540, ki_541, ki_542, ki_543, \
                         ki_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -hi_540[k]
                   + f_0 * ki_540[k];

        t_541[k] = -hi_541[k]
                   + f_0 * ki_541[k];

        t_542[k] = -hi_542[k]
                   + f_0 * ki_542[k];

        t_543[k] = -hi_543[k]
                   + f_0 * ki_543[k];

        t_544[k] = -hi_544[k]
                   + f_0 * ki_544[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, hi_545, hi_546, hi_547, hi_548, \
                         hi_549, ki_545, ki_546, ki_547, ki_548, \
                         ki_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -hi_545[k]
                   + f_0 * ki_545[k];

        t_546[k] = -hi_546[k]
                   + f_0 * ki_546[k];

        t_547[k] = -hi_547[k]
                   + f_0 * ki_547[k];

        t_548[k] = -hi_548[k]
                   + f_0 * ki_548[k];

        t_549[k] = -hi_549[k]
                   + f_0 * ki_549[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, hi_550, hi_551, hi_552, hi_553, \
                         hi_554, ki_550, ki_551, ki_552, ki_553, \
                         ki_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -hi_550[k]
                   + f_0 * ki_550[k];

        t_551[k] = -hi_551[k]
                   + f_0 * ki_551[k];

        t_552[k] = -hi_552[k]
                   + f_0 * ki_552[k];

        t_553[k] = -hi_553[k]
                   + f_0 * ki_553[k];

        t_554[k] = -hi_554[k]
                   + f_0 * ki_554[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, hi_555, hi_556, hi_557, hi_558, \
                         hi_559, ki_555, ki_556, ki_557, ki_558, \
                         ki_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -hi_555[k]
                   + f_0 * ki_555[k];

        t_556[k] = -hi_556[k]
                   + f_0 * ki_556[k];

        t_557[k] = -hi_557[k]
                   + f_0 * ki_557[k];

        t_558[k] = -hi_558[k]
                   + f_0 * ki_558[k];

        t_559[k] = -hi_559[k]
                   + f_0 * ki_559[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, hi_560, hi_561, hi_562, hi_563, \
                         hi_564, ki_560, ki_561, ki_562, ki_563, \
                         ki_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = -hi_560[k]
                   + f_0 * ki_560[k];

        t_561[k] = -hi_561[k]
                   + f_0 * ki_561[k];

        t_562[k] = -hi_562[k]
                   + f_0 * ki_562[k];

        t_563[k] = -hi_563[k]
                   + f_0 * ki_563[k];

        t_564[k] = -hi_564[k]
                   + f_0 * ki_564[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, hi_565, hi_566, hi_567, hi_568, \
                         hi_569, ki_565, ki_566, ki_567, ki_568, \
                         ki_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = -hi_565[k]
                   + f_0 * ki_565[k];

        t_566[k] = -hi_566[k]
                   + f_0 * ki_566[k];

        t_567[k] = -hi_567[k]
                   + f_0 * ki_567[k];

        t_568[k] = -hi_568[k]
                   + f_0 * ki_568[k];

        t_569[k] = -hi_569[k]
                   + f_0 * ki_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, hi_570, hi_571, hi_572, hi_573, \
                         hi_574, ki_570, ki_571, ki_572, ki_573, \
                         ki_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -hi_570[k]
                   + f_0 * ki_570[k];

        t_571[k] = -hi_571[k]
                   + f_0 * ki_571[k];

        t_572[k] = -hi_572[k]
                   + f_0 * ki_572[k];

        t_573[k] = -hi_573[k]
                   + f_0 * ki_573[k];

        t_574[k] = -hi_574[k]
                   + f_0 * ki_574[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, hi_575, hi_576, hi_577, hi_578, \
                         hi_579, ki_575, ki_576, ki_577, ki_578, \
                         ki_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = -hi_575[k]
                   + f_0 * ki_575[k];

        t_576[k] = -hi_576[k]
                   + f_0 * ki_576[k];

        t_577[k] = -hi_577[k]
                   + f_0 * ki_577[k];

        t_578[k] = -hi_578[k]
                   + f_0 * ki_578[k];

        t_579[k] = -hi_579[k]
                   + f_0 * ki_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, hi_580, hi_581, hi_582, hi_583, \
                         hi_584, ki_580, ki_581, ki_582, ki_583, \
                         ki_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -hi_580[k]
                   + f_0 * ki_580[k];

        t_581[k] = -hi_581[k]
                   + f_0 * ki_581[k];

        t_582[k] = -hi_582[k]
                   + f_0 * ki_582[k];

        t_583[k] = -hi_583[k]
                   + f_0 * ki_583[k];

        t_584[k] = -hi_584[k]
                   + f_0 * ki_584[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, t_590, hi_585, hi_586, hi_587, \
                         ki_585, ki_586, ki_587, ki_588, ki_589, \
                         ki_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -hi_585[k]
                   + f_0 * ki_585[k];

        t_586[k] = -hi_586[k]
                   + f_0 * ki_586[k];

        t_587[k] = -hi_587[k]
                   + f_0 * ki_587[k];

        t_588[k] = f_0 * ki_588[k];

        t_589[k] = f_0 * ki_589[k];

        t_590[k] = f_0 * ki_590[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, t_595, t_596, t_597, t_598, ki_591, \
                         ki_592, ki_593, ki_594, ki_595, ki_596, ki_597, \
                         ki_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_0 * ki_591[k];

        t_592[k] = f_0 * ki_592[k];

        t_593[k] = f_0 * ki_593[k];

        t_594[k] = f_0 * ki_594[k];

        t_595[k] = f_0 * ki_595[k];

        t_596[k] = f_0 * ki_596[k];

        t_597[k] = f_0 * ki_597[k];

        t_598[k] = f_0 * ki_598[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, t_604, t_605, t_606, ki_599, \
                         ki_600, ki_601, ki_602, ki_603, ki_604, ki_605, \
                         ki_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = f_0 * ki_599[k];

        t_600[k] = f_0 * ki_600[k];

        t_601[k] = f_0 * ki_601[k];

        t_602[k] = f_0 * ki_602[k];

        t_603[k] = f_0 * ki_603[k];

        t_604[k] = f_0 * ki_604[k];

        t_605[k] = f_0 * ki_605[k];

        t_606[k] = f_0 * ki_606[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ki, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ki_607 = buffer.data(ki + 607);
    const auto *ki_608 = buffer.data(ki + 608);
    const auto *ki_609 = buffer.data(ki + 609);
    const auto *ki_610 = buffer.data(ki + 610);
    const auto *ki_611 = buffer.data(ki + 611);
    const auto *ki_612 = buffer.data(ki + 612);
    const auto *ki_613 = buffer.data(ki + 613);
    const auto *ki_614 = buffer.data(ki + 614);
    const auto *ki_615 = buffer.data(ki + 615);
    const auto *ki_616 = buffer.data(ki + 616);
    const auto *ki_617 = buffer.data(ki + 617);
    const auto *ki_618 = buffer.data(ki + 618);
    const auto *ki_619 = buffer.data(ki + 619);
    const auto *ki_620 = buffer.data(ki + 620);
    const auto *ki_621 = buffer.data(ki + 621);
    const auto *ki_622 = buffer.data(ki + 622);
    const auto *ki_623 = buffer.data(ki + 623);
    const auto *ki_624 = buffer.data(ki + 624);
    const auto *ki_625 = buffer.data(ki + 625);
    const auto *ki_626 = buffer.data(ki + 626);
    const auto *ki_627 = buffer.data(ki + 627);
    const auto *ki_628 = buffer.data(ki + 628);
    const auto *ki_629 = buffer.data(ki + 629);
    const auto *ki_630 = buffer.data(ki + 630);
    const auto *ki_631 = buffer.data(ki + 631);
    const auto *ki_632 = buffer.data(ki + 632);
    const auto *ki_633 = buffer.data(ki + 633);
    const auto *ki_634 = buffer.data(ki + 634);
    const auto *ki_635 = buffer.data(ki + 635);
    const auto *ki_636 = buffer.data(ki + 636);
    const auto *ki_637 = buffer.data(ki + 637);
    const auto *ki_638 = buffer.data(ki + 638);
    const auto *ki_639 = buffer.data(ki + 639);
    const auto *ki_640 = buffer.data(ki + 640);
    const auto *ki_641 = buffer.data(ki + 641);
    const auto *ki_642 = buffer.data(ki + 642);
    const auto *ki_643 = buffer.data(ki + 643);
    const auto *ki_644 = buffer.data(ki + 644);
    const auto *ki_645 = buffer.data(ki + 645);
    const auto *ki_646 = buffer.data(ki + 646);
    const auto *ki_647 = buffer.data(ki + 647);
    const auto *ki_648 = buffer.data(ki + 648);
    const auto *ki_649 = buffer.data(ki + 649);
    const auto *ki_650 = buffer.data(ki + 650);
    const auto *ki_651 = buffer.data(ki + 651);
    const auto *ki_652 = buffer.data(ki + 652);
    const auto *ki_653 = buffer.data(ki + 653);
    const auto *ki_654 = buffer.data(ki + 654);
    const auto *ki_655 = buffer.data(ki + 655);
    const auto *ki_656 = buffer.data(ki + 656);
    const auto *ki_657 = buffer.data(ki + 657);
    const auto *ki_658 = buffer.data(ki + 658);
    const auto *ki_659 = buffer.data(ki + 659);
    const auto *ki_660 = buffer.data(ki + 660);
    const auto *ki_661 = buffer.data(ki + 661);
    const auto *ki_662 = buffer.data(ki + 662);
    const auto *ki_663 = buffer.data(ki + 663);
    const auto *ki_664 = buffer.data(ki + 664);
    const auto *ki_665 = buffer.data(ki + 665);
    const auto *ki_666 = buffer.data(ki + 666);
    const auto *ki_667 = buffer.data(ki + 667);
    const auto *ki_668 = buffer.data(ki + 668);
    const auto *ki_669 = buffer.data(ki + 669);
    const auto *ki_670 = buffer.data(ki + 670);
    const auto *ki_671 = buffer.data(ki + 671);
    const auto *ki_672 = buffer.data(ki + 672);
    const auto *ki_673 = buffer.data(ki + 673);
    const auto *ki_674 = buffer.data(ki + 674);
    const auto *ki_675 = buffer.data(ki + 675);
    const auto *ki_676 = buffer.data(ki + 676);
    const auto *ki_677 = buffer.data(ki + 677);
    const auto *ki_678 = buffer.data(ki + 678);
    const auto *ki_679 = buffer.data(ki + 679);
    const auto *ki_680 = buffer.data(ki + 680);
    const auto *ki_681 = buffer.data(ki + 681);
    const auto *ki_682 = buffer.data(ki + 682);
    const auto *ki_683 = buffer.data(ki + 683);
    const auto *ki_684 = buffer.data(ki + 684);
    const auto *ki_685 = buffer.data(ki + 685);
    const auto *ki_686 = buffer.data(ki + 686);
    const auto *ki_687 = buffer.data(ki + 687);
    const auto *ki_688 = buffer.data(ki + 688);
    const auto *ki_689 = buffer.data(ki + 689);
    const auto *ki_690 = buffer.data(ki + 690);
    const auto *ki_691 = buffer.data(ki + 691);
    const auto *ki_692 = buffer.data(ki + 692);
    const auto *ki_693 = buffer.data(ki + 693);
    const auto *ki_694 = buffer.data(ki + 694);
    const auto *ki_695 = buffer.data(ki + 695);
    const auto *ki_696 = buffer.data(ki + 696);
    const auto *ki_697 = buffer.data(ki + 697);
    const auto *ki_698 = buffer.data(ki + 698);
    const auto *ki_699 = buffer.data(ki + 699);
    const auto *ki_700 = buffer.data(ki + 700);
    const auto *ki_701 = buffer.data(ki + 701);
    const auto *ki_702 = buffer.data(ki + 702);
    const auto *ki_703 = buffer.data(ki + 703);
    const auto *ki_704 = buffer.data(ki + 704);
    const auto *ki_705 = buffer.data(ki + 705);
    const auto *ki_706 = buffer.data(ki + 706);
    const auto *ki_707 = buffer.data(ki + 707);
    const auto *ki_708 = buffer.data(ki + 708);
    const auto *ki_709 = buffer.data(ki + 709);
    const auto *ki_710 = buffer.data(ki + 710);
    const auto *ki_711 = buffer.data(ki + 711);
    const auto *ki_712 = buffer.data(ki + 712);
    const auto *ki_713 = buffer.data(ki + 713);
    const auto *ki_714 = buffer.data(ki + 714);
    const auto *ki_715 = buffer.data(ki + 715);
    const auto *ki_716 = buffer.data(ki + 716);
    const auto *ki_717 = buffer.data(ki + 717);
    const auto *ki_718 = buffer.data(ki + 718);
    const auto *ki_719 = buffer.data(ki + 719);
    const auto *ki_720 = buffer.data(ki + 720);
    const auto *ki_721 = buffer.data(ki + 721);
    const auto *ki_722 = buffer.data(ki + 722);
    const auto *ki_723 = buffer.data(ki + 723);
    const auto *ki_724 = buffer.data(ki + 724);
    const auto *ki_725 = buffer.data(ki + 725);
    const auto *ki_726 = buffer.data(ki + 726);
    const auto *ki_727 = buffer.data(ki + 727);
    const auto *ki_728 = buffer.data(ki + 728);
    const auto *ki_729 = buffer.data(ki + 729);
    const auto *ki_730 = buffer.data(ki + 730);
    const auto *ki_731 = buffer.data(ki + 731);
    const auto *ki_732 = buffer.data(ki + 732);
    const auto *ki_733 = buffer.data(ki + 733);
    const auto *ki_734 = buffer.data(ki + 734);
    const auto *ki_735 = buffer.data(ki + 735);
    const auto *ki_736 = buffer.data(ki + 736);
    const auto *ki_737 = buffer.data(ki + 737);
    const auto *ki_738 = buffer.data(ki + 738);
    const auto *ki_739 = buffer.data(ki + 739);
    const auto *ki_740 = buffer.data(ki + 740);
    const auto *ki_741 = buffer.data(ki + 741);
    const auto *ki_742 = buffer.data(ki + 742);
    const auto *ki_743 = buffer.data(ki + 743);
    const auto *ki_744 = buffer.data(ki + 744);
    const auto *ki_745 = buffer.data(ki + 745);
    const auto *ki_746 = buffer.data(ki + 746);
    const auto *ki_747 = buffer.data(ki + 747);
    const auto *ki_748 = buffer.data(ki + 748);
    const auto *ki_749 = buffer.data(ki + 749);
    const auto *ki_750 = buffer.data(ki + 750);
    const auto *ki_751 = buffer.data(ki + 751);
    const auto *ki_752 = buffer.data(ki + 752);
    const auto *ki_753 = buffer.data(ki + 753);
    const auto *ki_754 = buffer.data(ki + 754);
    const auto *ki_755 = buffer.data(ki + 755);
    const auto *ki_756 = buffer.data(ki + 756);
    const auto *ki_757 = buffer.data(ki + 757);
    const auto *ki_758 = buffer.data(ki + 758);
    const auto *ki_759 = buffer.data(ki + 759);
    const auto *ki_760 = buffer.data(ki + 760);
    const auto *ki_761 = buffer.data(ki + 761);
    const auto *ki_762 = buffer.data(ki + 762);
    const auto *ki_763 = buffer.data(ki + 763);
    const auto *ki_764 = buffer.data(ki + 764);
    const auto *ki_765 = buffer.data(ki + 765);
    const auto *ki_766 = buffer.data(ki + 766);
    const auto *ki_767 = buffer.data(ki + 767);
    const auto *ki_768 = buffer.data(ki + 768);
    const auto *ki_769 = buffer.data(ki + 769);
    const auto *ki_770 = buffer.data(ki + 770);
    const auto *ki_771 = buffer.data(ki + 771);
    const auto *ki_772 = buffer.data(ki + 772);
    const auto *ki_773 = buffer.data(ki + 773);
    const auto *ki_774 = buffer.data(ki + 774);
    const auto *ki_775 = buffer.data(ki + 775);
    const auto *ki_776 = buffer.data(ki + 776);
    const auto *ki_777 = buffer.data(ki + 777);
    const auto *ki_778 = buffer.data(ki + 778);
    const auto *ki_779 = buffer.data(ki + 779);
    const auto *ki_780 = buffer.data(ki + 780);
    const auto *ki_781 = buffer.data(ki + 781);
    const auto *ki_782 = buffer.data(ki + 782);
    const auto *ki_783 = buffer.data(ki + 783);

#pragma omp simd aligned(t_607, t_608, t_609, t_610, t_611, t_612, t_613, t_614, ki_607, \
                         ki_608, ki_609, ki_610, ki_611, ki_612, ki_613, \
                         ki_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_0 * ki_607[k];

        t_608[k] = f_0 * ki_608[k];

        t_609[k] = f_0 * ki_609[k];

        t_610[k] = f_0 * ki_610[k];

        t_611[k] = f_0 * ki_611[k];

        t_612[k] = f_0 * ki_612[k];

        t_613[k] = f_0 * ki_613[k];

        t_614[k] = f_0 * ki_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, t_620, t_621, t_622, ki_615, \
                         ki_616, ki_617, ki_618, ki_619, ki_620, ki_621, \
                         ki_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_0 * ki_615[k];

        t_616[k] = f_0 * ki_616[k];

        t_617[k] = f_0 * ki_617[k];

        t_618[k] = f_0 * ki_618[k];

        t_619[k] = f_0 * ki_619[k];

        t_620[k] = f_0 * ki_620[k];

        t_621[k] = f_0 * ki_621[k];

        t_622[k] = f_0 * ki_622[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, t_626, t_627, t_628, t_629, t_630, ki_623, \
                         ki_624, ki_625, ki_626, ki_627, ki_628, ki_629, \
                         ki_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_0 * ki_623[k];

        t_624[k] = f_0 * ki_624[k];

        t_625[k] = f_0 * ki_625[k];

        t_626[k] = f_0 * ki_626[k];

        t_627[k] = f_0 * ki_627[k];

        t_628[k] = f_0 * ki_628[k];

        t_629[k] = f_0 * ki_629[k];

        t_630[k] = f_0 * ki_630[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, t_635, t_636, t_637, t_638, ki_631, \
                         ki_632, ki_633, ki_634, ki_635, ki_636, ki_637, \
                         ki_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = f_0 * ki_631[k];

        t_632[k] = f_0 * ki_632[k];

        t_633[k] = f_0 * ki_633[k];

        t_634[k] = f_0 * ki_634[k];

        t_635[k] = f_0 * ki_635[k];

        t_636[k] = f_0 * ki_636[k];

        t_637[k] = f_0 * ki_637[k];

        t_638[k] = f_0 * ki_638[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, t_643, t_644, t_645, t_646, ki_639, \
                         ki_640, ki_641, ki_642, ki_643, ki_644, ki_645, \
                         ki_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_0 * ki_639[k];

        t_640[k] = f_0 * ki_640[k];

        t_641[k] = f_0 * ki_641[k];

        t_642[k] = f_0 * ki_642[k];

        t_643[k] = f_0 * ki_643[k];

        t_644[k] = f_0 * ki_644[k];

        t_645[k] = f_0 * ki_645[k];

        t_646[k] = f_0 * ki_646[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, t_651, t_652, t_653, t_654, ki_647, \
                         ki_648, ki_649, ki_650, ki_651, ki_652, ki_653, \
                         ki_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_0 * ki_647[k];

        t_648[k] = f_0 * ki_648[k];

        t_649[k] = f_0 * ki_649[k];

        t_650[k] = f_0 * ki_650[k];

        t_651[k] = f_0 * ki_651[k];

        t_652[k] = f_0 * ki_652[k];

        t_653[k] = f_0 * ki_653[k];

        t_654[k] = f_0 * ki_654[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, t_660, t_661, t_662, ki_655, \
                         ki_656, ki_657, ki_658, ki_659, ki_660, ki_661, \
                         ki_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_0 * ki_655[k];

        t_656[k] = f_0 * ki_656[k];

        t_657[k] = f_0 * ki_657[k];

        t_658[k] = f_0 * ki_658[k];

        t_659[k] = f_0 * ki_659[k];

        t_660[k] = f_0 * ki_660[k];

        t_661[k] = f_0 * ki_661[k];

        t_662[k] = f_0 * ki_662[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, t_667, t_668, t_669, t_670, ki_663, \
                         ki_664, ki_665, ki_666, ki_667, ki_668, ki_669, \
                         ki_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_0 * ki_663[k];

        t_664[k] = f_0 * ki_664[k];

        t_665[k] = f_0 * ki_665[k];

        t_666[k] = f_0 * ki_666[k];

        t_667[k] = f_0 * ki_667[k];

        t_668[k] = f_0 * ki_668[k];

        t_669[k] = f_0 * ki_669[k];

        t_670[k] = f_0 * ki_670[k];
    }

#pragma omp simd aligned(t_671, t_672, t_673, t_674, t_675, t_676, t_677, t_678, ki_671, \
                         ki_672, ki_673, ki_674, ki_675, ki_676, ki_677, \
                         ki_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_671[k] = f_0 * ki_671[k];

        t_672[k] = f_0 * ki_672[k];

        t_673[k] = f_0 * ki_673[k];

        t_674[k] = f_0 * ki_674[k];

        t_675[k] = f_0 * ki_675[k];

        t_676[k] = f_0 * ki_676[k];

        t_677[k] = f_0 * ki_677[k];

        t_678[k] = f_0 * ki_678[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, t_682, t_683, t_684, t_685, t_686, ki_679, \
                         ki_680, ki_681, ki_682, ki_683, ki_684, ki_685, \
                         ki_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_0 * ki_679[k];

        t_680[k] = f_0 * ki_680[k];

        t_681[k] = f_0 * ki_681[k];

        t_682[k] = f_0 * ki_682[k];

        t_683[k] = f_0 * ki_683[k];

        t_684[k] = f_0 * ki_684[k];

        t_685[k] = f_0 * ki_685[k];

        t_686[k] = f_0 * ki_686[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, t_691, t_692, t_693, t_694, ki_687, \
                         ki_688, ki_689, ki_690, ki_691, ki_692, ki_693, \
                         ki_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_0 * ki_687[k];

        t_688[k] = f_0 * ki_688[k];

        t_689[k] = f_0 * ki_689[k];

        t_690[k] = f_0 * ki_690[k];

        t_691[k] = f_0 * ki_691[k];

        t_692[k] = f_0 * ki_692[k];

        t_693[k] = f_0 * ki_693[k];

        t_694[k] = f_0 * ki_694[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, t_700, t_701, t_702, ki_695, \
                         ki_696, ki_697, ki_698, ki_699, ki_700, ki_701, \
                         ki_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_0 * ki_695[k];

        t_696[k] = f_0 * ki_696[k];

        t_697[k] = f_0 * ki_697[k];

        t_698[k] = f_0 * ki_698[k];

        t_699[k] = f_0 * ki_699[k];

        t_700[k] = f_0 * ki_700[k];

        t_701[k] = f_0 * ki_701[k];

        t_702[k] = f_0 * ki_702[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, t_706, t_707, t_708, t_709, t_710, ki_703, \
                         ki_704, ki_705, ki_706, ki_707, ki_708, ki_709, \
                         ki_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_0 * ki_703[k];

        t_704[k] = f_0 * ki_704[k];

        t_705[k] = f_0 * ki_705[k];

        t_706[k] = f_0 * ki_706[k];

        t_707[k] = f_0 * ki_707[k];

        t_708[k] = f_0 * ki_708[k];

        t_709[k] = f_0 * ki_709[k];

        t_710[k] = f_0 * ki_710[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, t_714, t_715, t_716, t_717, t_718, ki_711, \
                         ki_712, ki_713, ki_714, ki_715, ki_716, ki_717, \
                         ki_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = f_0 * ki_711[k];

        t_712[k] = f_0 * ki_712[k];

        t_713[k] = f_0 * ki_713[k];

        t_714[k] = f_0 * ki_714[k];

        t_715[k] = f_0 * ki_715[k];

        t_716[k] = f_0 * ki_716[k];

        t_717[k] = f_0 * ki_717[k];

        t_718[k] = f_0 * ki_718[k];
    }

#pragma omp simd aligned(t_719, t_720, t_721, t_722, t_723, t_724, t_725, t_726, ki_719, \
                         ki_720, ki_721, ki_722, ki_723, ki_724, ki_725, \
                         ki_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_719[k] = f_0 * ki_719[k];

        t_720[k] = f_0 * ki_720[k];

        t_721[k] = f_0 * ki_721[k];

        t_722[k] = f_0 * ki_722[k];

        t_723[k] = f_0 * ki_723[k];

        t_724[k] = f_0 * ki_724[k];

        t_725[k] = f_0 * ki_725[k];

        t_726[k] = f_0 * ki_726[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, t_730, t_731, t_732, t_733, t_734, ki_727, \
                         ki_728, ki_729, ki_730, ki_731, ki_732, ki_733, \
                         ki_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_0 * ki_727[k];

        t_728[k] = f_0 * ki_728[k];

        t_729[k] = f_0 * ki_729[k];

        t_730[k] = f_0 * ki_730[k];

        t_731[k] = f_0 * ki_731[k];

        t_732[k] = f_0 * ki_732[k];

        t_733[k] = f_0 * ki_733[k];

        t_734[k] = f_0 * ki_734[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, t_740, t_741, t_742, ki_735, \
                         ki_736, ki_737, ki_738, ki_739, ki_740, ki_741, \
                         ki_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_0 * ki_735[k];

        t_736[k] = f_0 * ki_736[k];

        t_737[k] = f_0 * ki_737[k];

        t_738[k] = f_0 * ki_738[k];

        t_739[k] = f_0 * ki_739[k];

        t_740[k] = f_0 * ki_740[k];

        t_741[k] = f_0 * ki_741[k];

        t_742[k] = f_0 * ki_742[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, t_748, t_749, t_750, ki_743, \
                         ki_744, ki_745, ki_746, ki_747, ki_748, ki_749, \
                         ki_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_0 * ki_743[k];

        t_744[k] = f_0 * ki_744[k];

        t_745[k] = f_0 * ki_745[k];

        t_746[k] = f_0 * ki_746[k];

        t_747[k] = f_0 * ki_747[k];

        t_748[k] = f_0 * ki_748[k];

        t_749[k] = f_0 * ki_749[k];

        t_750[k] = f_0 * ki_750[k];
    }

#pragma omp simd aligned(t_751, t_752, t_753, t_754, t_755, t_756, t_757, t_758, ki_751, \
                         ki_752, ki_753, ki_754, ki_755, ki_756, ki_757, \
                         ki_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_751[k] = f_0 * ki_751[k];

        t_752[k] = f_0 * ki_752[k];

        t_753[k] = f_0 * ki_753[k];

        t_754[k] = f_0 * ki_754[k];

        t_755[k] = f_0 * ki_755[k];

        t_756[k] = f_0 * ki_756[k];

        t_757[k] = f_0 * ki_757[k];

        t_758[k] = f_0 * ki_758[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, t_762, t_763, t_764, t_765, t_766, ki_759, \
                         ki_760, ki_761, ki_762, ki_763, ki_764, ki_765, \
                         ki_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_0 * ki_759[k];

        t_760[k] = f_0 * ki_760[k];

        t_761[k] = f_0 * ki_761[k];

        t_762[k] = f_0 * ki_762[k];

        t_763[k] = f_0 * ki_763[k];

        t_764[k] = f_0 * ki_764[k];

        t_765[k] = f_0 * ki_765[k];

        t_766[k] = f_0 * ki_766[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, t_770, t_771, t_772, t_773, t_774, ki_767, \
                         ki_768, ki_769, ki_770, ki_771, ki_772, ki_773, \
                         ki_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_0 * ki_767[k];

        t_768[k] = f_0 * ki_768[k];

        t_769[k] = f_0 * ki_769[k];

        t_770[k] = f_0 * ki_770[k];

        t_771[k] = f_0 * ki_771[k];

        t_772[k] = f_0 * ki_772[k];

        t_773[k] = f_0 * ki_773[k];

        t_774[k] = f_0 * ki_774[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, t_780, t_781, t_782, ki_775, \
                         ki_776, ki_777, ki_778, ki_779, ki_780, ki_781, \
                         ki_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_0 * ki_775[k];

        t_776[k] = f_0 * ki_776[k];

        t_777[k] = f_0 * ki_777[k];

        t_778[k] = f_0 * ki_778[k];

        t_779[k] = f_0 * ki_779[k];

        t_780[k] = f_0 * ki_780[k];

        t_781[k] = f_0 * ki_781[k];

        t_782[k] = f_0 * ki_782[k];
    }

#pragma omp simd aligned(t_783, ki_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = f_0 * ki_783[k];
    }
}

auto
compute_prim_geom_10_ii_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t hi, const size_t ki,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ii_electron_repulsion_0_piece0(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_0_piece1(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_0_piece2(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_0_piece3(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_0_piece4(buffer, target, ki, ncols, alpha);
}

static auto
compute_prim_geom_10_ii_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
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

    const auto *hi_0 = buffer.data(hi + 0);
    const auto *hi_1 = buffer.data(hi + 1);
    const auto *hi_2 = buffer.data(hi + 2);
    const auto *hi_3 = buffer.data(hi + 3);
    const auto *hi_4 = buffer.data(hi + 4);
    const auto *hi_5 = buffer.data(hi + 5);
    const auto *hi_6 = buffer.data(hi + 6);
    const auto *hi_7 = buffer.data(hi + 7);
    const auto *hi_8 = buffer.data(hi + 8);
    const auto *hi_9 = buffer.data(hi + 9);
    const auto *hi_10 = buffer.data(hi + 10);
    const auto *hi_11 = buffer.data(hi + 11);
    const auto *hi_12 = buffer.data(hi + 12);
    const auto *hi_13 = buffer.data(hi + 13);
    const auto *hi_14 = buffer.data(hi + 14);
    const auto *hi_15 = buffer.data(hi + 15);
    const auto *hi_16 = buffer.data(hi + 16);
    const auto *hi_17 = buffer.data(hi + 17);
    const auto *hi_18 = buffer.data(hi + 18);
    const auto *hi_19 = buffer.data(hi + 19);
    const auto *hi_20 = buffer.data(hi + 20);
    const auto *hi_21 = buffer.data(hi + 21);
    const auto *hi_22 = buffer.data(hi + 22);
    const auto *hi_23 = buffer.data(hi + 23);
    const auto *hi_24 = buffer.data(hi + 24);
    const auto *hi_25 = buffer.data(hi + 25);
    const auto *hi_26 = buffer.data(hi + 26);
    const auto *hi_27 = buffer.data(hi + 27);
    const auto *hi_28 = buffer.data(hi + 28);
    const auto *hi_29 = buffer.data(hi + 29);
    const auto *hi_30 = buffer.data(hi + 30);
    const auto *hi_31 = buffer.data(hi + 31);
    const auto *hi_32 = buffer.data(hi + 32);
    const auto *hi_33 = buffer.data(hi + 33);
    const auto *hi_34 = buffer.data(hi + 34);
    const auto *hi_35 = buffer.data(hi + 35);
    const auto *hi_36 = buffer.data(hi + 36);
    const auto *hi_37 = buffer.data(hi + 37);
    const auto *hi_38 = buffer.data(hi + 38);
    const auto *hi_39 = buffer.data(hi + 39);
    const auto *hi_40 = buffer.data(hi + 40);
    const auto *hi_41 = buffer.data(hi + 41);
    const auto *hi_42 = buffer.data(hi + 42);
    const auto *hi_43 = buffer.data(hi + 43);
    const auto *hi_44 = buffer.data(hi + 44);
    const auto *hi_45 = buffer.data(hi + 45);
    const auto *hi_46 = buffer.data(hi + 46);
    const auto *hi_47 = buffer.data(hi + 47);
    const auto *hi_48 = buffer.data(hi + 48);
    const auto *hi_49 = buffer.data(hi + 49);
    const auto *hi_50 = buffer.data(hi + 50);
    const auto *hi_51 = buffer.data(hi + 51);
    const auto *hi_52 = buffer.data(hi + 52);
    const auto *hi_53 = buffer.data(hi + 53);
    const auto *hi_54 = buffer.data(hi + 54);
    const auto *hi_55 = buffer.data(hi + 55);
    const auto *hi_56 = buffer.data(hi + 56);
    const auto *hi_57 = buffer.data(hi + 57);
    const auto *hi_58 = buffer.data(hi + 58);
    const auto *hi_59 = buffer.data(hi + 59);
    const auto *hi_60 = buffer.data(hi + 60);
    const auto *hi_61 = buffer.data(hi + 61);
    const auto *hi_62 = buffer.data(hi + 62);
    const auto *hi_63 = buffer.data(hi + 63);
    const auto *hi_64 = buffer.data(hi + 64);
    const auto *hi_65 = buffer.data(hi + 65);
    const auto *hi_66 = buffer.data(hi + 66);
    const auto *hi_67 = buffer.data(hi + 67);
    const auto *hi_68 = buffer.data(hi + 68);
    const auto *hi_69 = buffer.data(hi + 69);
    const auto *hi_70 = buffer.data(hi + 70);
    const auto *hi_71 = buffer.data(hi + 71);
    const auto *hi_72 = buffer.data(hi + 72);
    const auto *hi_73 = buffer.data(hi + 73);
    const auto *hi_74 = buffer.data(hi + 74);
    const auto *hi_75 = buffer.data(hi + 75);
    const auto *hi_76 = buffer.data(hi + 76);
    const auto *hi_77 = buffer.data(hi + 77);
    const auto *hi_78 = buffer.data(hi + 78);
    const auto *hi_79 = buffer.data(hi + 79);
    const auto *hi_80 = buffer.data(hi + 80);
    const auto *hi_81 = buffer.data(hi + 81);
    const auto *hi_82 = buffer.data(hi + 82);
    const auto *hi_83 = buffer.data(hi + 83);
    const auto *hi_84 = buffer.data(hi + 84);
    const auto *hi_85 = buffer.data(hi + 85);
    const auto *hi_86 = buffer.data(hi + 86);
    const auto *hi_87 = buffer.data(hi + 87);
    const auto *hi_88 = buffer.data(hi + 88);
    const auto *hi_89 = buffer.data(hi + 89);
    const auto *hi_90 = buffer.data(hi + 90);

    const auto *ki_28 = buffer.data(ki + 28);
    const auto *ki_29 = buffer.data(ki + 29);
    const auto *ki_30 = buffer.data(ki + 30);
    const auto *ki_31 = buffer.data(ki + 31);
    const auto *ki_32 = buffer.data(ki + 32);
    const auto *ki_33 = buffer.data(ki + 33);
    const auto *ki_34 = buffer.data(ki + 34);
    const auto *ki_35 = buffer.data(ki + 35);
    const auto *ki_36 = buffer.data(ki + 36);
    const auto *ki_37 = buffer.data(ki + 37);
    const auto *ki_38 = buffer.data(ki + 38);
    const auto *ki_39 = buffer.data(ki + 39);
    const auto *ki_40 = buffer.data(ki + 40);
    const auto *ki_41 = buffer.data(ki + 41);
    const auto *ki_42 = buffer.data(ki + 42);
    const auto *ki_43 = buffer.data(ki + 43);
    const auto *ki_44 = buffer.data(ki + 44);
    const auto *ki_45 = buffer.data(ki + 45);
    const auto *ki_46 = buffer.data(ki + 46);
    const auto *ki_47 = buffer.data(ki + 47);
    const auto *ki_48 = buffer.data(ki + 48);
    const auto *ki_49 = buffer.data(ki + 49);
    const auto *ki_50 = buffer.data(ki + 50);
    const auto *ki_51 = buffer.data(ki + 51);
    const auto *ki_52 = buffer.data(ki + 52);
    const auto *ki_53 = buffer.data(ki + 53);
    const auto *ki_54 = buffer.data(ki + 54);
    const auto *ki_55 = buffer.data(ki + 55);
    const auto *ki_84 = buffer.data(ki + 84);
    const auto *ki_85 = buffer.data(ki + 85);
    const auto *ki_86 = buffer.data(ki + 86);
    const auto *ki_87 = buffer.data(ki + 87);
    const auto *ki_88 = buffer.data(ki + 88);
    const auto *ki_89 = buffer.data(ki + 89);
    const auto *ki_90 = buffer.data(ki + 90);
    const auto *ki_91 = buffer.data(ki + 91);
    const auto *ki_92 = buffer.data(ki + 92);
    const auto *ki_93 = buffer.data(ki + 93);
    const auto *ki_94 = buffer.data(ki + 94);
    const auto *ki_95 = buffer.data(ki + 95);
    const auto *ki_96 = buffer.data(ki + 96);
    const auto *ki_97 = buffer.data(ki + 97);
    const auto *ki_98 = buffer.data(ki + 98);
    const auto *ki_99 = buffer.data(ki + 99);
    const auto *ki_100 = buffer.data(ki + 100);
    const auto *ki_101 = buffer.data(ki + 101);
    const auto *ki_102 = buffer.data(ki + 102);
    const auto *ki_103 = buffer.data(ki + 103);
    const auto *ki_104 = buffer.data(ki + 104);
    const auto *ki_105 = buffer.data(ki + 105);
    const auto *ki_106 = buffer.data(ki + 106);
    const auto *ki_107 = buffer.data(ki + 107);
    const auto *ki_108 = buffer.data(ki + 108);
    const auto *ki_109 = buffer.data(ki + 109);
    const auto *ki_110 = buffer.data(ki + 110);
    const auto *ki_111 = buffer.data(ki + 111);
    const auto *ki_112 = buffer.data(ki + 112);
    const auto *ki_113 = buffer.data(ki + 113);
    const auto *ki_114 = buffer.data(ki + 114);
    const auto *ki_115 = buffer.data(ki + 115);
    const auto *ki_116 = buffer.data(ki + 116);
    const auto *ki_117 = buffer.data(ki + 117);
    const auto *ki_118 = buffer.data(ki + 118);
    const auto *ki_119 = buffer.data(ki + 119);
    const auto *ki_120 = buffer.data(ki + 120);
    const auto *ki_121 = buffer.data(ki + 121);
    const auto *ki_122 = buffer.data(ki + 122);
    const auto *ki_123 = buffer.data(ki + 123);
    const auto *ki_124 = buffer.data(ki + 124);
    const auto *ki_125 = buffer.data(ki + 125);
    const auto *ki_126 = buffer.data(ki + 126);
    const auto *ki_127 = buffer.data(ki + 127);
    const auto *ki_128 = buffer.data(ki + 128);
    const auto *ki_129 = buffer.data(ki + 129);
    const auto *ki_130 = buffer.data(ki + 130);
    const auto *ki_131 = buffer.data(ki + 131);
    const auto *ki_132 = buffer.data(ki + 132);
    const auto *ki_133 = buffer.data(ki + 133);
    const auto *ki_134 = buffer.data(ki + 134);
    const auto *ki_135 = buffer.data(ki + 135);
    const auto *ki_136 = buffer.data(ki + 136);
    const auto *ki_137 = buffer.data(ki + 137);
    const auto *ki_138 = buffer.data(ki + 138);
    const auto *ki_139 = buffer.data(ki + 139);
    const auto *ki_168 = buffer.data(ki + 168);
    const auto *ki_169 = buffer.data(ki + 169);
    const auto *ki_170 = buffer.data(ki + 170);
    const auto *ki_171 = buffer.data(ki + 171);
    const auto *ki_172 = buffer.data(ki + 172);
    const auto *ki_173 = buffer.data(ki + 173);
    const auto *ki_174 = buffer.data(ki + 174);
    const auto *ki_175 = buffer.data(ki + 175);
    const auto *ki_176 = buffer.data(ki + 176);
    const auto *ki_177 = buffer.data(ki + 177);
    const auto *ki_178 = buffer.data(ki + 178);
    const auto *ki_179 = buffer.data(ki + 179);
    const auto *ki_180 = buffer.data(ki + 180);
    const auto *ki_181 = buffer.data(ki + 181);
    const auto *ki_182 = buffer.data(ki + 182);
    const auto *ki_183 = buffer.data(ki + 183);
    const auto *ki_184 = buffer.data(ki + 184);
    const auto *ki_185 = buffer.data(ki + 185);
    const auto *ki_186 = buffer.data(ki + 186);
    const auto *ki_187 = buffer.data(ki + 187);
    const auto *ki_188 = buffer.data(ki + 188);
    const auto *ki_189 = buffer.data(ki + 189);
    const auto *ki_190 = buffer.data(ki + 190);
    const auto *ki_191 = buffer.data(ki + 191);
    const auto *ki_192 = buffer.data(ki + 192);
    const auto *ki_193 = buffer.data(ki + 193);
    const auto *ki_194 = buffer.data(ki + 194);
    const auto *ki_195 = buffer.data(ki + 195);
    const auto *ki_196 = buffer.data(ki + 196);
    const auto *ki_197 = buffer.data(ki + 197);
    const auto *ki_198 = buffer.data(ki + 198);
    const auto *ki_199 = buffer.data(ki + 199);
    const auto *ki_200 = buffer.data(ki + 200);
    const auto *ki_201 = buffer.data(ki + 201);
    const auto *ki_202 = buffer.data(ki + 202);
    const auto *ki_203 = buffer.data(ki + 203);
    const auto *ki_204 = buffer.data(ki + 204);
    const auto *ki_205 = buffer.data(ki + 205);
    const auto *ki_206 = buffer.data(ki + 206);
    const auto *ki_207 = buffer.data(ki + 207);
    const auto *ki_208 = buffer.data(ki + 208);
    const auto *ki_209 = buffer.data(ki + 209);
    const auto *ki_210 = buffer.data(ki + 210);
    const auto *ki_211 = buffer.data(ki + 211);
    const auto *ki_212 = buffer.data(ki + 212);
    const auto *ki_213 = buffer.data(ki + 213);
    const auto *ki_214 = buffer.data(ki + 214);
    const auto *ki_215 = buffer.data(ki + 215);
    const auto *ki_216 = buffer.data(ki + 216);
    const auto *ki_217 = buffer.data(ki + 217);
    const auto *ki_218 = buffer.data(ki + 218);
    const auto *ki_219 = buffer.data(ki + 219);
    const auto *ki_220 = buffer.data(ki + 220);
    const auto *ki_221 = buffer.data(ki + 221);
    const auto *ki_222 = buffer.data(ki + 222);
    const auto *ki_223 = buffer.data(ki + 223);
    const auto *ki_224 = buffer.data(ki + 224);
    const auto *ki_225 = buffer.data(ki + 225);
    const auto *ki_226 = buffer.data(ki + 226);
    const auto *ki_227 = buffer.data(ki + 227);
    const auto *ki_228 = buffer.data(ki + 228);
    const auto *ki_229 = buffer.data(ki + 229);
    const auto *ki_230 = buffer.data(ki + 230);
    const auto *ki_231 = buffer.data(ki + 231);
    const auto *ki_232 = buffer.data(ki + 232);
    const auto *ki_233 = buffer.data(ki + 233);
    const auto *ki_234 = buffer.data(ki + 234);
    const auto *ki_235 = buffer.data(ki + 235);
    const auto *ki_236 = buffer.data(ki + 236);
    const auto *ki_237 = buffer.data(ki + 237);
    const auto *ki_238 = buffer.data(ki + 238);
    const auto *ki_239 = buffer.data(ki + 239);
    const auto *ki_240 = buffer.data(ki + 240);
    const auto *ki_241 = buffer.data(ki + 241);
    const auto *ki_242 = buffer.data(ki + 242);
    const auto *ki_243 = buffer.data(ki + 243);
    const auto *ki_244 = buffer.data(ki + 244);
    const auto *ki_245 = buffer.data(ki + 245);
    const auto *ki_246 = buffer.data(ki + 246);
    const auto *ki_247 = buffer.data(ki + 247);
    const auto *ki_248 = buffer.data(ki + 248);
    const auto *ki_249 = buffer.data(ki + 249);
    const auto *ki_250 = buffer.data(ki + 250);
    const auto *ki_251 = buffer.data(ki + 251);
    const auto *ki_280 = buffer.data(ki + 280);
    const auto *ki_281 = buffer.data(ki + 281);
    const auto *ki_282 = buffer.data(ki + 282);
    const auto *ki_283 = buffer.data(ki + 283);
    const auto *ki_284 = buffer.data(ki + 284);
    const auto *ki_285 = buffer.data(ki + 285);
    const auto *ki_286 = buffer.data(ki + 286);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ki_28, ki_29, ki_30, ki_31, \
                         ki_32, ki_33, ki_34, ki_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ki_28[k];

        t_1[k] = f_0 * ki_29[k];

        t_2[k] = f_0 * ki_30[k];

        t_3[k] = f_0 * ki_31[k];

        t_4[k] = f_0 * ki_32[k];

        t_5[k] = f_0 * ki_33[k];

        t_6[k] = f_0 * ki_34[k];

        t_7[k] = f_0 * ki_35[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, ki_36, ki_37, ki_38, \
                         ki_39, ki_40, ki_41, ki_42, ki_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ki_36[k];

        t_9[k] = f_0 * ki_37[k];

        t_10[k] = f_0 * ki_38[k];

        t_11[k] = f_0 * ki_39[k];

        t_12[k] = f_0 * ki_40[k];

        t_13[k] = f_0 * ki_41[k];

        t_14[k] = f_0 * ki_42[k];

        t_15[k] = f_0 * ki_43[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, ki_44, ki_45, ki_46, \
                         ki_47, ki_48, ki_49, ki_50, ki_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * ki_44[k];

        t_17[k] = f_0 * ki_45[k];

        t_18[k] = f_0 * ki_46[k];

        t_19[k] = f_0 * ki_47[k];

        t_20[k] = f_0 * ki_48[k];

        t_21[k] = f_0 * ki_49[k];

        t_22[k] = f_0 * ki_50[k];

        t_23[k] = f_0 * ki_51[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, hi_0, hi_1, ki_52, ki_53, ki_54, \
                         ki_55, ki_84, ki_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * ki_52[k];

        t_25[k] = f_0 * ki_53[k];

        t_26[k] = f_0 * ki_54[k];

        t_27[k] = f_0 * ki_55[k];

        t_28[k] = -hi_0[k]
                  + f_0 * ki_84[k];

        t_29[k] = -hi_1[k]
                  + f_0 * ki_85[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, hi_2, hi_3, hi_4, hi_5, hi_6, ki_86, \
                         ki_87, ki_88, ki_89, ki_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -hi_2[k]
                  + f_0 * ki_86[k];

        t_31[k] = -hi_3[k]
                  + f_0 * ki_87[k];

        t_32[k] = -hi_4[k]
                  + f_0 * ki_88[k];

        t_33[k] = -hi_5[k]
                  + f_0 * ki_89[k];

        t_34[k] = -hi_6[k]
                  + f_0 * ki_90[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, hi_7, hi_8, hi_9, hi_10, hi_11, ki_91, \
                         ki_92, ki_93, ki_94, ki_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -hi_7[k]
                  + f_0 * ki_91[k];

        t_36[k] = -hi_8[k]
                  + f_0 * ki_92[k];

        t_37[k] = -hi_9[k]
                  + f_0 * ki_93[k];

        t_38[k] = -hi_10[k]
                  + f_0 * ki_94[k];

        t_39[k] = -hi_11[k]
                  + f_0 * ki_95[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, hi_12, hi_13, hi_14, hi_15, hi_16, \
                         ki_96, ki_97, ki_98, ki_99, ki_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -hi_12[k]
                  + f_0 * ki_96[k];

        t_41[k] = -hi_13[k]
                  + f_0 * ki_97[k];

        t_42[k] = -hi_14[k]
                  + f_0 * ki_98[k];

        t_43[k] = -hi_15[k]
                  + f_0 * ki_99[k];

        t_44[k] = -hi_16[k]
                  + f_0 * ki_100[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hi_17, hi_18, hi_19, hi_20, hi_21, \
                         ki_101, ki_102, ki_103, ki_104, ki_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -hi_17[k]
                  + f_0 * ki_101[k];

        t_46[k] = -hi_18[k]
                  + f_0 * ki_102[k];

        t_47[k] = -hi_19[k]
                  + f_0 * ki_103[k];

        t_48[k] = -hi_20[k]
                  + f_0 * ki_104[k];

        t_49[k] = -hi_21[k]
                  + f_0 * ki_105[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hi_22, hi_23, hi_24, hi_25, hi_26, \
                         ki_106, ki_107, ki_108, ki_109, ki_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -hi_22[k]
                  + f_0 * ki_106[k];

        t_51[k] = -hi_23[k]
                  + f_0 * ki_107[k];

        t_52[k] = -hi_24[k]
                  + f_0 * ki_108[k];

        t_53[k] = -hi_25[k]
                  + f_0 * ki_109[k];

        t_54[k] = -hi_26[k]
                  + f_0 * ki_110[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, hi_27, ki_111, ki_112, \
                         ki_113, ki_114, ki_115, ki_116, ki_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -hi_27[k]
                  + f_0 * ki_111[k];

        t_56[k] = f_0 * ki_112[k];

        t_57[k] = f_0 * ki_113[k];

        t_58[k] = f_0 * ki_114[k];

        t_59[k] = f_0 * ki_115[k];

        t_60[k] = f_0 * ki_116[k];

        t_61[k] = f_0 * ki_117[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, t_67, t_68, t_69, ki_118, ki_119, \
                         ki_120, ki_121, ki_122, ki_123, ki_124, \
                         ki_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * ki_118[k];

        t_63[k] = f_0 * ki_119[k];

        t_64[k] = f_0 * ki_120[k];

        t_65[k] = f_0 * ki_121[k];

        t_66[k] = f_0 * ki_122[k];

        t_67[k] = f_0 * ki_123[k];

        t_68[k] = f_0 * ki_124[k];

        t_69[k] = f_0 * ki_125[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, t_75, t_76, t_77, ki_126, ki_127, \
                         ki_128, ki_129, ki_130, ki_131, ki_132, \
                         ki_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_0 * ki_126[k];

        t_71[k] = f_0 * ki_127[k];

        t_72[k] = f_0 * ki_128[k];

        t_73[k] = f_0 * ki_129[k];

        t_74[k] = f_0 * ki_130[k];

        t_75[k] = f_0 * ki_131[k];

        t_76[k] = f_0 * ki_132[k];

        t_77[k] = f_0 * ki_133[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, t_84, hi_28, ki_134, ki_135, \
                         ki_136, ki_137, ki_138, ki_139, ki_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * ki_134[k];

        t_79[k] = f_0 * ki_135[k];

        t_80[k] = f_0 * ki_136[k];

        t_81[k] = f_0 * ki_137[k];

        t_82[k] = f_0 * ki_138[k];

        t_83[k] = f_0 * ki_139[k];

        t_84[k] = -2.0 * hi_28[k]
                  + f_0 * ki_168[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, hi_29, hi_30, hi_31, hi_32, hi_33, \
                         ki_169, ki_170, ki_171, ki_172, ki_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * hi_29[k]
                  + f_0 * ki_169[k];

        t_86[k] = -2.0 * hi_30[k]
                  + f_0 * ki_170[k];

        t_87[k] = -2.0 * hi_31[k]
                  + f_0 * ki_171[k];

        t_88[k] = -2.0 * hi_32[k]
                  + f_0 * ki_172[k];

        t_89[k] = -2.0 * hi_33[k]
                  + f_0 * ki_173[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, hi_34, hi_35, hi_36, hi_37, hi_38, \
                         ki_174, ki_175, ki_176, ki_177, ki_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * hi_34[k]
                  + f_0 * ki_174[k];

        t_91[k] = -2.0 * hi_35[k]
                  + f_0 * ki_175[k];

        t_92[k] = -2.0 * hi_36[k]
                  + f_0 * ki_176[k];

        t_93[k] = -2.0 * hi_37[k]
                  + f_0 * ki_177[k];

        t_94[k] = -2.0 * hi_38[k]
                  + f_0 * ki_178[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, hi_39, hi_40, hi_41, hi_42, hi_43, \
                         ki_179, ki_180, ki_181, ki_182, ki_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * hi_39[k]
                  + f_0 * ki_179[k];

        t_96[k] = -2.0 * hi_40[k]
                  + f_0 * ki_180[k];

        t_97[k] = -2.0 * hi_41[k]
                  + f_0 * ki_181[k];

        t_98[k] = -2.0 * hi_42[k]
                  + f_0 * ki_182[k];

        t_99[k] = -2.0 * hi_43[k]
                  + f_0 * ki_183[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, hi_44, hi_45, hi_46, hi_47, hi_48, \
                         ki_184, ki_185, ki_186, ki_187, ki_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * hi_44[k]
                   + f_0 * ki_184[k];

        t_101[k] = -2.0 * hi_45[k]
                   + f_0 * ki_185[k];

        t_102[k] = -2.0 * hi_46[k]
                   + f_0 * ki_186[k];

        t_103[k] = -2.0 * hi_47[k]
                   + f_0 * ki_187[k];

        t_104[k] = -2.0 * hi_48[k]
                   + f_0 * ki_188[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, hi_49, hi_50, hi_51, hi_52, hi_53, \
                         ki_189, ki_190, ki_191, ki_192, ki_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * hi_49[k]
                   + f_0 * ki_189[k];

        t_106[k] = -2.0 * hi_50[k]
                   + f_0 * ki_190[k];

        t_107[k] = -2.0 * hi_51[k]
                   + f_0 * ki_191[k];

        t_108[k] = -2.0 * hi_52[k]
                   + f_0 * ki_192[k];

        t_109[k] = -2.0 * hi_53[k]
                   + f_0 * ki_193[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, hi_54, hi_55, hi_56, hi_57, hi_58, \
                         ki_194, ki_195, ki_196, ki_197, ki_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * hi_54[k]
                   + f_0 * ki_194[k];

        t_111[k] = -2.0 * hi_55[k]
                   + f_0 * ki_195[k];

        t_112[k] = -hi_56[k]
                   + f_0 * ki_196[k];

        t_113[k] = -hi_57[k]
                   + f_0 * ki_197[k];

        t_114[k] = -hi_58[k]
                   + f_0 * ki_198[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, hi_59, hi_60, hi_61, hi_62, hi_63, \
                         ki_199, ki_200, ki_201, ki_202, ki_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -hi_59[k]
                   + f_0 * ki_199[k];

        t_116[k] = -hi_60[k]
                   + f_0 * ki_200[k];

        t_117[k] = -hi_61[k]
                   + f_0 * ki_201[k];

        t_118[k] = -hi_62[k]
                   + f_0 * ki_202[k];

        t_119[k] = -hi_63[k]
                   + f_0 * ki_203[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, hi_64, hi_65, hi_66, hi_67, hi_68, \
                         ki_204, ki_205, ki_206, ki_207, ki_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -hi_64[k]
                   + f_0 * ki_204[k];

        t_121[k] = -hi_65[k]
                   + f_0 * ki_205[k];

        t_122[k] = -hi_66[k]
                   + f_0 * ki_206[k];

        t_123[k] = -hi_67[k]
                   + f_0 * ki_207[k];

        t_124[k] = -hi_68[k]
                   + f_0 * ki_208[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, hi_69, hi_70, hi_71, hi_72, hi_73, \
                         ki_209, ki_210, ki_211, ki_212, ki_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -hi_69[k]
                   + f_0 * ki_209[k];

        t_126[k] = -hi_70[k]
                   + f_0 * ki_210[k];

        t_127[k] = -hi_71[k]
                   + f_0 * ki_211[k];

        t_128[k] = -hi_72[k]
                   + f_0 * ki_212[k];

        t_129[k] = -hi_73[k]
                   + f_0 * ki_213[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, hi_74, hi_75, hi_76, hi_77, hi_78, \
                         ki_214, ki_215, ki_216, ki_217, ki_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -hi_74[k]
                   + f_0 * ki_214[k];

        t_131[k] = -hi_75[k]
                   + f_0 * ki_215[k];

        t_132[k] = -hi_76[k]
                   + f_0 * ki_216[k];

        t_133[k] = -hi_77[k]
                   + f_0 * ki_217[k];

        t_134[k] = -hi_78[k]
                   + f_0 * ki_218[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, hi_79, hi_80, hi_81, hi_82, hi_83, \
                         ki_219, ki_220, ki_221, ki_222, ki_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -hi_79[k]
                   + f_0 * ki_219[k];

        t_136[k] = -hi_80[k]
                   + f_0 * ki_220[k];

        t_137[k] = -hi_81[k]
                   + f_0 * ki_221[k];

        t_138[k] = -hi_82[k]
                   + f_0 * ki_222[k];

        t_139[k] = -hi_83[k]
                   + f_0 * ki_223[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, t_146, t_147, ki_224, \
                         ki_225, ki_226, ki_227, ki_228, ki_229, ki_230, \
                         ki_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * ki_224[k];

        t_141[k] = f_0 * ki_225[k];

        t_142[k] = f_0 * ki_226[k];

        t_143[k] = f_0 * ki_227[k];

        t_144[k] = f_0 * ki_228[k];

        t_145[k] = f_0 * ki_229[k];

        t_146[k] = f_0 * ki_230[k];

        t_147[k] = f_0 * ki_231[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, t_153, t_154, t_155, ki_232, \
                         ki_233, ki_234, ki_235, ki_236, ki_237, ki_238, \
                         ki_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * ki_232[k];

        t_149[k] = f_0 * ki_233[k];

        t_150[k] = f_0 * ki_234[k];

        t_151[k] = f_0 * ki_235[k];

        t_152[k] = f_0 * ki_236[k];

        t_153[k] = f_0 * ki_237[k];

        t_154[k] = f_0 * ki_238[k];

        t_155[k] = f_0 * ki_239[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, t_162, t_163, ki_240, \
                         ki_241, ki_242, ki_243, ki_244, ki_245, ki_246, \
                         ki_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * ki_240[k];

        t_157[k] = f_0 * ki_241[k];

        t_158[k] = f_0 * ki_242[k];

        t_159[k] = f_0 * ki_243[k];

        t_160[k] = f_0 * ki_244[k];

        t_161[k] = f_0 * ki_245[k];

        t_162[k] = f_0 * ki_246[k];

        t_163[k] = f_0 * ki_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, t_169, hi_84, hi_85, ki_248, \
                         ki_249, ki_250, ki_251, ki_280, ki_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_0 * ki_248[k];

        t_165[k] = f_0 * ki_249[k];

        t_166[k] = f_0 * ki_250[k];

        t_167[k] = f_0 * ki_251[k];

        t_168[k] = -3.0 * hi_84[k]
                   + f_0 * ki_280[k];

        t_169[k] = -3.0 * hi_85[k]
                   + f_0 * ki_281[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, hi_86, hi_87, hi_88, hi_89, hi_90, \
                         ki_282, ki_283, ki_284, ki_285, ki_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * hi_86[k]
                   + f_0 * ki_282[k];

        t_171[k] = -3.0 * hi_87[k]
                   + f_0 * ki_283[k];

        t_172[k] = -3.0 * hi_88[k]
                   + f_0 * ki_284[k];

        t_173[k] = -3.0 * hi_89[k]
                   + f_0 * ki_285[k];

        t_174[k] = -3.0 * hi_90[k]
                   + f_0 * ki_286[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_334 = buffer.data(target + 334);

    const auto *hi_91 = buffer.data(hi + 91);
    const auto *hi_92 = buffer.data(hi + 92);
    const auto *hi_93 = buffer.data(hi + 93);
    const auto *hi_94 = buffer.data(hi + 94);
    const auto *hi_95 = buffer.data(hi + 95);
    const auto *hi_96 = buffer.data(hi + 96);
    const auto *hi_97 = buffer.data(hi + 97);
    const auto *hi_98 = buffer.data(hi + 98);
    const auto *hi_99 = buffer.data(hi + 99);
    const auto *hi_100 = buffer.data(hi + 100);
    const auto *hi_101 = buffer.data(hi + 101);
    const auto *hi_102 = buffer.data(hi + 102);
    const auto *hi_103 = buffer.data(hi + 103);
    const auto *hi_104 = buffer.data(hi + 104);
    const auto *hi_105 = buffer.data(hi + 105);
    const auto *hi_106 = buffer.data(hi + 106);
    const auto *hi_107 = buffer.data(hi + 107);
    const auto *hi_108 = buffer.data(hi + 108);
    const auto *hi_109 = buffer.data(hi + 109);
    const auto *hi_110 = buffer.data(hi + 110);
    const auto *hi_111 = buffer.data(hi + 111);
    const auto *hi_112 = buffer.data(hi + 112);
    const auto *hi_113 = buffer.data(hi + 113);
    const auto *hi_114 = buffer.data(hi + 114);
    const auto *hi_115 = buffer.data(hi + 115);
    const auto *hi_116 = buffer.data(hi + 116);
    const auto *hi_117 = buffer.data(hi + 117);
    const auto *hi_118 = buffer.data(hi + 118);
    const auto *hi_119 = buffer.data(hi + 119);
    const auto *hi_120 = buffer.data(hi + 120);
    const auto *hi_121 = buffer.data(hi + 121);
    const auto *hi_122 = buffer.data(hi + 122);
    const auto *hi_123 = buffer.data(hi + 123);
    const auto *hi_124 = buffer.data(hi + 124);
    const auto *hi_125 = buffer.data(hi + 125);
    const auto *hi_126 = buffer.data(hi + 126);
    const auto *hi_127 = buffer.data(hi + 127);
    const auto *hi_128 = buffer.data(hi + 128);
    const auto *hi_129 = buffer.data(hi + 129);
    const auto *hi_130 = buffer.data(hi + 130);
    const auto *hi_131 = buffer.data(hi + 131);
    const auto *hi_132 = buffer.data(hi + 132);
    const auto *hi_133 = buffer.data(hi + 133);
    const auto *hi_134 = buffer.data(hi + 134);
    const auto *hi_135 = buffer.data(hi + 135);
    const auto *hi_136 = buffer.data(hi + 136);
    const auto *hi_137 = buffer.data(hi + 137);
    const auto *hi_138 = buffer.data(hi + 138);
    const auto *hi_139 = buffer.data(hi + 139);
    const auto *hi_140 = buffer.data(hi + 140);
    const auto *hi_141 = buffer.data(hi + 141);
    const auto *hi_142 = buffer.data(hi + 142);
    const auto *hi_143 = buffer.data(hi + 143);
    const auto *hi_144 = buffer.data(hi + 144);
    const auto *hi_145 = buffer.data(hi + 145);
    const auto *hi_146 = buffer.data(hi + 146);
    const auto *hi_147 = buffer.data(hi + 147);
    const auto *hi_148 = buffer.data(hi + 148);
    const auto *hi_149 = buffer.data(hi + 149);
    const auto *hi_150 = buffer.data(hi + 150);
    const auto *hi_151 = buffer.data(hi + 151);
    const auto *hi_152 = buffer.data(hi + 152);
    const auto *hi_153 = buffer.data(hi + 153);
    const auto *hi_154 = buffer.data(hi + 154);
    const auto *hi_155 = buffer.data(hi + 155);
    const auto *hi_156 = buffer.data(hi + 156);
    const auto *hi_157 = buffer.data(hi + 157);
    const auto *hi_158 = buffer.data(hi + 158);
    const auto *hi_159 = buffer.data(hi + 159);
    const auto *hi_160 = buffer.data(hi + 160);
    const auto *hi_161 = buffer.data(hi + 161);
    const auto *hi_162 = buffer.data(hi + 162);
    const auto *hi_163 = buffer.data(hi + 163);
    const auto *hi_164 = buffer.data(hi + 164);
    const auto *hi_165 = buffer.data(hi + 165);
    const auto *hi_166 = buffer.data(hi + 166);
    const auto *hi_167 = buffer.data(hi + 167);
    const auto *hi_168 = buffer.data(hi + 168);
    const auto *hi_169 = buffer.data(hi + 169);
    const auto *hi_170 = buffer.data(hi + 170);
    const auto *hi_171 = buffer.data(hi + 171);
    const auto *hi_172 = buffer.data(hi + 172);
    const auto *hi_173 = buffer.data(hi + 173);
    const auto *hi_174 = buffer.data(hi + 174);
    const auto *hi_175 = buffer.data(hi + 175);
    const auto *hi_176 = buffer.data(hi + 176);
    const auto *hi_177 = buffer.data(hi + 177);
    const auto *hi_178 = buffer.data(hi + 178);
    const auto *hi_179 = buffer.data(hi + 179);
    const auto *hi_180 = buffer.data(hi + 180);
    const auto *hi_181 = buffer.data(hi + 181);
    const auto *hi_182 = buffer.data(hi + 182);
    const auto *hi_183 = buffer.data(hi + 183);
    const auto *hi_184 = buffer.data(hi + 184);
    const auto *hi_185 = buffer.data(hi + 185);
    const auto *hi_186 = buffer.data(hi + 186);
    const auto *hi_187 = buffer.data(hi + 187);
    const auto *hi_188 = buffer.data(hi + 188);
    const auto *hi_189 = buffer.data(hi + 189);
    const auto *hi_190 = buffer.data(hi + 190);
    const auto *hi_191 = buffer.data(hi + 191);
    const auto *hi_192 = buffer.data(hi + 192);
    const auto *hi_193 = buffer.data(hi + 193);
    const auto *hi_194 = buffer.data(hi + 194);
    const auto *hi_195 = buffer.data(hi + 195);
    const auto *hi_196 = buffer.data(hi + 196);
    const auto *hi_197 = buffer.data(hi + 197);
    const auto *hi_198 = buffer.data(hi + 198);
    const auto *hi_199 = buffer.data(hi + 199);
    const auto *hi_200 = buffer.data(hi + 200);
    const auto *hi_201 = buffer.data(hi + 201);
    const auto *hi_202 = buffer.data(hi + 202);
    const auto *hi_203 = buffer.data(hi + 203);
    const auto *hi_204 = buffer.data(hi + 204);
    const auto *hi_205 = buffer.data(hi + 205);
    const auto *hi_206 = buffer.data(hi + 206);
    const auto *hi_207 = buffer.data(hi + 207);
    const auto *hi_208 = buffer.data(hi + 208);
    const auto *hi_209 = buffer.data(hi + 209);
    const auto *hi_210 = buffer.data(hi + 210);
    const auto *hi_211 = buffer.data(hi + 211);
    const auto *hi_212 = buffer.data(hi + 212);
    const auto *hi_213 = buffer.data(hi + 213);
    const auto *hi_214 = buffer.data(hi + 214);
    const auto *hi_215 = buffer.data(hi + 215);
    const auto *hi_216 = buffer.data(hi + 216);
    const auto *hi_217 = buffer.data(hi + 217);
    const auto *hi_218 = buffer.data(hi + 218);
    const auto *hi_219 = buffer.data(hi + 219);
    const auto *hi_220 = buffer.data(hi + 220);
    const auto *hi_221 = buffer.data(hi + 221);
    const auto *hi_222 = buffer.data(hi + 222);

    const auto *ki_287 = buffer.data(ki + 287);
    const auto *ki_288 = buffer.data(ki + 288);
    const auto *ki_289 = buffer.data(ki + 289);
    const auto *ki_290 = buffer.data(ki + 290);
    const auto *ki_291 = buffer.data(ki + 291);
    const auto *ki_292 = buffer.data(ki + 292);
    const auto *ki_293 = buffer.data(ki + 293);
    const auto *ki_294 = buffer.data(ki + 294);
    const auto *ki_295 = buffer.data(ki + 295);
    const auto *ki_296 = buffer.data(ki + 296);
    const auto *ki_297 = buffer.data(ki + 297);
    const auto *ki_298 = buffer.data(ki + 298);
    const auto *ki_299 = buffer.data(ki + 299);
    const auto *ki_300 = buffer.data(ki + 300);
    const auto *ki_301 = buffer.data(ki + 301);
    const auto *ki_302 = buffer.data(ki + 302);
    const auto *ki_303 = buffer.data(ki + 303);
    const auto *ki_304 = buffer.data(ki + 304);
    const auto *ki_305 = buffer.data(ki + 305);
    const auto *ki_306 = buffer.data(ki + 306);
    const auto *ki_307 = buffer.data(ki + 307);
    const auto *ki_308 = buffer.data(ki + 308);
    const auto *ki_309 = buffer.data(ki + 309);
    const auto *ki_310 = buffer.data(ki + 310);
    const auto *ki_311 = buffer.data(ki + 311);
    const auto *ki_312 = buffer.data(ki + 312);
    const auto *ki_313 = buffer.data(ki + 313);
    const auto *ki_314 = buffer.data(ki + 314);
    const auto *ki_315 = buffer.data(ki + 315);
    const auto *ki_316 = buffer.data(ki + 316);
    const auto *ki_317 = buffer.data(ki + 317);
    const auto *ki_318 = buffer.data(ki + 318);
    const auto *ki_319 = buffer.data(ki + 319);
    const auto *ki_320 = buffer.data(ki + 320);
    const auto *ki_321 = buffer.data(ki + 321);
    const auto *ki_322 = buffer.data(ki + 322);
    const auto *ki_323 = buffer.data(ki + 323);
    const auto *ki_324 = buffer.data(ki + 324);
    const auto *ki_325 = buffer.data(ki + 325);
    const auto *ki_326 = buffer.data(ki + 326);
    const auto *ki_327 = buffer.data(ki + 327);
    const auto *ki_328 = buffer.data(ki + 328);
    const auto *ki_329 = buffer.data(ki + 329);
    const auto *ki_330 = buffer.data(ki + 330);
    const auto *ki_331 = buffer.data(ki + 331);
    const auto *ki_332 = buffer.data(ki + 332);
    const auto *ki_333 = buffer.data(ki + 333);
    const auto *ki_334 = buffer.data(ki + 334);
    const auto *ki_335 = buffer.data(ki + 335);
    const auto *ki_336 = buffer.data(ki + 336);
    const auto *ki_337 = buffer.data(ki + 337);
    const auto *ki_338 = buffer.data(ki + 338);
    const auto *ki_339 = buffer.data(ki + 339);
    const auto *ki_340 = buffer.data(ki + 340);
    const auto *ki_341 = buffer.data(ki + 341);
    const auto *ki_342 = buffer.data(ki + 342);
    const auto *ki_343 = buffer.data(ki + 343);
    const auto *ki_344 = buffer.data(ki + 344);
    const auto *ki_345 = buffer.data(ki + 345);
    const auto *ki_346 = buffer.data(ki + 346);
    const auto *ki_347 = buffer.data(ki + 347);
    const auto *ki_348 = buffer.data(ki + 348);
    const auto *ki_349 = buffer.data(ki + 349);
    const auto *ki_350 = buffer.data(ki + 350);
    const auto *ki_351 = buffer.data(ki + 351);
    const auto *ki_352 = buffer.data(ki + 352);
    const auto *ki_353 = buffer.data(ki + 353);
    const auto *ki_354 = buffer.data(ki + 354);
    const auto *ki_355 = buffer.data(ki + 355);
    const auto *ki_356 = buffer.data(ki + 356);
    const auto *ki_357 = buffer.data(ki + 357);
    const auto *ki_358 = buffer.data(ki + 358);
    const auto *ki_359 = buffer.data(ki + 359);
    const auto *ki_360 = buffer.data(ki + 360);
    const auto *ki_361 = buffer.data(ki + 361);
    const auto *ki_362 = buffer.data(ki + 362);
    const auto *ki_363 = buffer.data(ki + 363);
    const auto *ki_364 = buffer.data(ki + 364);
    const auto *ki_365 = buffer.data(ki + 365);
    const auto *ki_366 = buffer.data(ki + 366);
    const auto *ki_367 = buffer.data(ki + 367);
    const auto *ki_368 = buffer.data(ki + 368);
    const auto *ki_369 = buffer.data(ki + 369);
    const auto *ki_370 = buffer.data(ki + 370);
    const auto *ki_371 = buffer.data(ki + 371);
    const auto *ki_372 = buffer.data(ki + 372);
    const auto *ki_373 = buffer.data(ki + 373);
    const auto *ki_374 = buffer.data(ki + 374);
    const auto *ki_375 = buffer.data(ki + 375);
    const auto *ki_376 = buffer.data(ki + 376);
    const auto *ki_377 = buffer.data(ki + 377);
    const auto *ki_378 = buffer.data(ki + 378);
    const auto *ki_379 = buffer.data(ki + 379);
    const auto *ki_380 = buffer.data(ki + 380);
    const auto *ki_381 = buffer.data(ki + 381);
    const auto *ki_382 = buffer.data(ki + 382);
    const auto *ki_383 = buffer.data(ki + 383);
    const auto *ki_384 = buffer.data(ki + 384);
    const auto *ki_385 = buffer.data(ki + 385);
    const auto *ki_386 = buffer.data(ki + 386);
    const auto *ki_387 = buffer.data(ki + 387);
    const auto *ki_388 = buffer.data(ki + 388);
    const auto *ki_389 = buffer.data(ki + 389);
    const auto *ki_390 = buffer.data(ki + 390);
    const auto *ki_391 = buffer.data(ki + 391);
    const auto *ki_420 = buffer.data(ki + 420);
    const auto *ki_421 = buffer.data(ki + 421);
    const auto *ki_422 = buffer.data(ki + 422);
    const auto *ki_423 = buffer.data(ki + 423);
    const auto *ki_424 = buffer.data(ki + 424);
    const auto *ki_425 = buffer.data(ki + 425);
    const auto *ki_426 = buffer.data(ki + 426);
    const auto *ki_427 = buffer.data(ki + 427);
    const auto *ki_428 = buffer.data(ki + 428);
    const auto *ki_429 = buffer.data(ki + 429);
    const auto *ki_430 = buffer.data(ki + 430);
    const auto *ki_431 = buffer.data(ki + 431);
    const auto *ki_432 = buffer.data(ki + 432);
    const auto *ki_433 = buffer.data(ki + 433);
    const auto *ki_434 = buffer.data(ki + 434);
    const auto *ki_435 = buffer.data(ki + 435);
    const auto *ki_436 = buffer.data(ki + 436);
    const auto *ki_437 = buffer.data(ki + 437);
    const auto *ki_438 = buffer.data(ki + 438);
    const auto *ki_439 = buffer.data(ki + 439);
    const auto *ki_440 = buffer.data(ki + 440);
    const auto *ki_441 = buffer.data(ki + 441);
    const auto *ki_442 = buffer.data(ki + 442);
    const auto *ki_443 = buffer.data(ki + 443);
    const auto *ki_444 = buffer.data(ki + 444);
    const auto *ki_445 = buffer.data(ki + 445);
    const auto *ki_446 = buffer.data(ki + 446);
    const auto *ki_447 = buffer.data(ki + 447);
    const auto *ki_448 = buffer.data(ki + 448);
    const auto *ki_449 = buffer.data(ki + 449);
    const auto *ki_450 = buffer.data(ki + 450);
    const auto *ki_451 = buffer.data(ki + 451);
    const auto *ki_452 = buffer.data(ki + 452);
    const auto *ki_453 = buffer.data(ki + 453);
    const auto *ki_454 = buffer.data(ki + 454);
    const auto *ki_455 = buffer.data(ki + 455);
    const auto *ki_456 = buffer.data(ki + 456);
    const auto *ki_457 = buffer.data(ki + 457);
    const auto *ki_458 = buffer.data(ki + 458);
    const auto *ki_459 = buffer.data(ki + 459);
    const auto *ki_460 = buffer.data(ki + 460);
    const auto *ki_461 = buffer.data(ki + 461);
    const auto *ki_462 = buffer.data(ki + 462);
    const auto *ki_463 = buffer.data(ki + 463);
    const auto *ki_464 = buffer.data(ki + 464);
    const auto *ki_465 = buffer.data(ki + 465);
    const auto *ki_466 = buffer.data(ki + 466);
    const auto *ki_467 = buffer.data(ki + 467);
    const auto *ki_468 = buffer.data(ki + 468);
    const auto *ki_469 = buffer.data(ki + 469);
    const auto *ki_470 = buffer.data(ki + 470);
    const auto *ki_471 = buffer.data(ki + 471);
    const auto *ki_472 = buffer.data(ki + 472);
    const auto *ki_473 = buffer.data(ki + 473);
    const auto *ki_474 = buffer.data(ki + 474);

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, hi_91, hi_92, hi_93, hi_94, hi_95, \
                         ki_287, ki_288, ki_289, ki_290, ki_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * hi_91[k]
                   + f_0 * ki_287[k];

        t_176[k] = -3.0 * hi_92[k]
                   + f_0 * ki_288[k];

        t_177[k] = -3.0 * hi_93[k]
                   + f_0 * ki_289[k];

        t_178[k] = -3.0 * hi_94[k]
                   + f_0 * ki_290[k];

        t_179[k] = -3.0 * hi_95[k]
                   + f_0 * ki_291[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, hi_96, hi_97, hi_98, hi_99, \
                         hi_100, ki_292, ki_293, ki_294, ki_295, \
                         ki_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * hi_96[k]
                   + f_0 * ki_292[k];

        t_181[k] = -3.0 * hi_97[k]
                   + f_0 * ki_293[k];

        t_182[k] = -3.0 * hi_98[k]
                   + f_0 * ki_294[k];

        t_183[k] = -3.0 * hi_99[k]
                   + f_0 * ki_295[k];

        t_184[k] = -3.0 * hi_100[k]
                   + f_0 * ki_296[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, hi_101, hi_102, hi_103, hi_104, \
                         hi_105, ki_297, ki_298, ki_299, ki_300, \
                         ki_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * hi_101[k]
                   + f_0 * ki_297[k];

        t_186[k] = -3.0 * hi_102[k]
                   + f_0 * ki_298[k];

        t_187[k] = -3.0 * hi_103[k]
                   + f_0 * ki_299[k];

        t_188[k] = -3.0 * hi_104[k]
                   + f_0 * ki_300[k];

        t_189[k] = -3.0 * hi_105[k]
                   + f_0 * ki_301[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, hi_106, hi_107, hi_108, hi_109, \
                         hi_110, ki_302, ki_303, ki_304, ki_305, \
                         ki_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * hi_106[k]
                   + f_0 * ki_302[k];

        t_191[k] = -3.0 * hi_107[k]
                   + f_0 * ki_303[k];

        t_192[k] = -3.0 * hi_108[k]
                   + f_0 * ki_304[k];

        t_193[k] = -3.0 * hi_109[k]
                   + f_0 * ki_305[k];

        t_194[k] = -3.0 * hi_110[k]
                   + f_0 * ki_306[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, hi_111, hi_112, hi_113, hi_114, \
                         hi_115, ki_307, ki_308, ki_309, ki_310, \
                         ki_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * hi_111[k]
                   + f_0 * ki_307[k];

        t_196[k] = -2.0 * hi_112[k]
                   + f_0 * ki_308[k];

        t_197[k] = -2.0 * hi_113[k]
                   + f_0 * ki_309[k];

        t_198[k] = -2.0 * hi_114[k]
                   + f_0 * ki_310[k];

        t_199[k] = -2.0 * hi_115[k]
                   + f_0 * ki_311[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, hi_116, hi_117, hi_118, hi_119, \
                         hi_120, ki_312, ki_313, ki_314, ki_315, \
                         ki_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * hi_116[k]
                   + f_0 * ki_312[k];

        t_201[k] = -2.0 * hi_117[k]
                   + f_0 * ki_313[k];

        t_202[k] = -2.0 * hi_118[k]
                   + f_0 * ki_314[k];

        t_203[k] = -2.0 * hi_119[k]
                   + f_0 * ki_315[k];

        t_204[k] = -2.0 * hi_120[k]
                   + f_0 * ki_316[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, hi_121, hi_122, hi_123, hi_124, \
                         hi_125, ki_317, ki_318, ki_319, ki_320, \
                         ki_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * hi_121[k]
                   + f_0 * ki_317[k];

        t_206[k] = -2.0 * hi_122[k]
                   + f_0 * ki_318[k];

        t_207[k] = -2.0 * hi_123[k]
                   + f_0 * ki_319[k];

        t_208[k] = -2.0 * hi_124[k]
                   + f_0 * ki_320[k];

        t_209[k] = -2.0 * hi_125[k]
                   + f_0 * ki_321[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, hi_126, hi_127, hi_128, hi_129, \
                         hi_130, ki_322, ki_323, ki_324, ki_325, \
                         ki_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * hi_126[k]
                   + f_0 * ki_322[k];

        t_211[k] = -2.0 * hi_127[k]
                   + f_0 * ki_323[k];

        t_212[k] = -2.0 * hi_128[k]
                   + f_0 * ki_324[k];

        t_213[k] = -2.0 * hi_129[k]
                   + f_0 * ki_325[k];

        t_214[k] = -2.0 * hi_130[k]
                   + f_0 * ki_326[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, hi_131, hi_132, hi_133, hi_134, \
                         hi_135, ki_327, ki_328, ki_329, ki_330, \
                         ki_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * hi_131[k]
                   + f_0 * ki_327[k];

        t_216[k] = -2.0 * hi_132[k]
                   + f_0 * ki_328[k];

        t_217[k] = -2.0 * hi_133[k]
                   + f_0 * ki_329[k];

        t_218[k] = -2.0 * hi_134[k]
                   + f_0 * ki_330[k];

        t_219[k] = -2.0 * hi_135[k]
                   + f_0 * ki_331[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, hi_136, hi_137, hi_138, hi_139, \
                         hi_140, ki_332, ki_333, ki_334, ki_335, \
                         ki_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * hi_136[k]
                   + f_0 * ki_332[k];

        t_221[k] = -2.0 * hi_137[k]
                   + f_0 * ki_333[k];

        t_222[k] = -2.0 * hi_138[k]
                   + f_0 * ki_334[k];

        t_223[k] = -2.0 * hi_139[k]
                   + f_0 * ki_335[k];

        t_224[k] = -hi_140[k]
                   + f_0 * ki_336[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, hi_141, hi_142, hi_143, hi_144, \
                         hi_145, ki_337, ki_338, ki_339, ki_340, \
                         ki_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -hi_141[k]
                   + f_0 * ki_337[k];

        t_226[k] = -hi_142[k]
                   + f_0 * ki_338[k];

        t_227[k] = -hi_143[k]
                   + f_0 * ki_339[k];

        t_228[k] = -hi_144[k]
                   + f_0 * ki_340[k];

        t_229[k] = -hi_145[k]
                   + f_0 * ki_341[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, hi_146, hi_147, hi_148, hi_149, \
                         hi_150, ki_342, ki_343, ki_344, ki_345, \
                         ki_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -hi_146[k]
                   + f_0 * ki_342[k];

        t_231[k] = -hi_147[k]
                   + f_0 * ki_343[k];

        t_232[k] = -hi_148[k]
                   + f_0 * ki_344[k];

        t_233[k] = -hi_149[k]
                   + f_0 * ki_345[k];

        t_234[k] = -hi_150[k]
                   + f_0 * ki_346[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, hi_151, hi_152, hi_153, hi_154, \
                         hi_155, ki_347, ki_348, ki_349, ki_350, \
                         ki_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -hi_151[k]
                   + f_0 * ki_347[k];

        t_236[k] = -hi_152[k]
                   + f_0 * ki_348[k];

        t_237[k] = -hi_153[k]
                   + f_0 * ki_349[k];

        t_238[k] = -hi_154[k]
                   + f_0 * ki_350[k];

        t_239[k] = -hi_155[k]
                   + f_0 * ki_351[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, hi_156, hi_157, hi_158, hi_159, \
                         hi_160, ki_352, ki_353, ki_354, ki_355, \
                         ki_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -hi_156[k]
                   + f_0 * ki_352[k];

        t_241[k] = -hi_157[k]
                   + f_0 * ki_353[k];

        t_242[k] = -hi_158[k]
                   + f_0 * ki_354[k];

        t_243[k] = -hi_159[k]
                   + f_0 * ki_355[k];

        t_244[k] = -hi_160[k]
                   + f_0 * ki_356[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, hi_161, hi_162, hi_163, hi_164, \
                         hi_165, ki_357, ki_358, ki_359, ki_360, \
                         ki_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -hi_161[k]
                   + f_0 * ki_357[k];

        t_246[k] = -hi_162[k]
                   + f_0 * ki_358[k];

        t_247[k] = -hi_163[k]
                   + f_0 * ki_359[k];

        t_248[k] = -hi_164[k]
                   + f_0 * ki_360[k];

        t_249[k] = -hi_165[k]
                   + f_0 * ki_361[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, t_255, t_256, hi_166, hi_167, \
                         ki_362, ki_363, ki_364, ki_365, ki_366, ki_367, \
                         ki_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -hi_166[k]
                   + f_0 * ki_362[k];

        t_251[k] = -hi_167[k]
                   + f_0 * ki_363[k];

        t_252[k] = f_0 * ki_364[k];

        t_253[k] = f_0 * ki_365[k];

        t_254[k] = f_0 * ki_366[k];

        t_255[k] = f_0 * ki_367[k];

        t_256[k] = f_0 * ki_368[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, t_262, t_263, t_264, ki_369, \
                         ki_370, ki_371, ki_372, ki_373, ki_374, ki_375, \
                         ki_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_0 * ki_369[k];

        t_258[k] = f_0 * ki_370[k];

        t_259[k] = f_0 * ki_371[k];

        t_260[k] = f_0 * ki_372[k];

        t_261[k] = f_0 * ki_373[k];

        t_262[k] = f_0 * ki_374[k];

        t_263[k] = f_0 * ki_375[k];

        t_264[k] = f_0 * ki_376[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, t_270, t_271, t_272, ki_377, \
                         ki_378, ki_379, ki_380, ki_381, ki_382, ki_383, \
                         ki_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_0 * ki_377[k];

        t_266[k] = f_0 * ki_378[k];

        t_267[k] = f_0 * ki_379[k];

        t_268[k] = f_0 * ki_380[k];

        t_269[k] = f_0 * ki_381[k];

        t_270[k] = f_0 * ki_382[k];

        t_271[k] = f_0 * ki_383[k];

        t_272[k] = f_0 * ki_384[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, t_278, t_279, ki_385, ki_386, \
                         ki_387, ki_388, ki_389, ki_390, ki_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_0 * ki_385[k];

        t_274[k] = f_0 * ki_386[k];

        t_275[k] = f_0 * ki_387[k];

        t_276[k] = f_0 * ki_388[k];

        t_277[k] = f_0 * ki_389[k];

        t_278[k] = f_0 * ki_390[k];

        t_279[k] = f_0 * ki_391[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, hi_168, hi_169, hi_170, hi_171, \
                         hi_172, ki_420, ki_421, ki_422, ki_423, \
                         ki_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -4.0 * hi_168[k]
                   + f_0 * ki_420[k];

        t_281[k] = -4.0 * hi_169[k]
                   + f_0 * ki_421[k];

        t_282[k] = -4.0 * hi_170[k]
                   + f_0 * ki_422[k];

        t_283[k] = -4.0 * hi_171[k]
                   + f_0 * ki_423[k];

        t_284[k] = -4.0 * hi_172[k]
                   + f_0 * ki_424[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, hi_173, hi_174, hi_175, hi_176, \
                         hi_177, ki_425, ki_426, ki_427, ki_428, \
                         ki_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -4.0 * hi_173[k]
                   + f_0 * ki_425[k];

        t_286[k] = -4.0 * hi_174[k]
                   + f_0 * ki_426[k];

        t_287[k] = -4.0 * hi_175[k]
                   + f_0 * ki_427[k];

        t_288[k] = -4.0 * hi_176[k]
                   + f_0 * ki_428[k];

        t_289[k] = -4.0 * hi_177[k]
                   + f_0 * ki_429[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, hi_178, hi_179, hi_180, hi_181, \
                         hi_182, ki_430, ki_431, ki_432, ki_433, \
                         ki_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -4.0 * hi_178[k]
                   + f_0 * ki_430[k];

        t_291[k] = -4.0 * hi_179[k]
                   + f_0 * ki_431[k];

        t_292[k] = -4.0 * hi_180[k]
                   + f_0 * ki_432[k];

        t_293[k] = -4.0 * hi_181[k]
                   + f_0 * ki_433[k];

        t_294[k] = -4.0 * hi_182[k]
                   + f_0 * ki_434[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, hi_183, hi_184, hi_185, hi_186, \
                         hi_187, ki_435, ki_436, ki_437, ki_438, \
                         ki_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * hi_183[k]
                   + f_0 * ki_435[k];

        t_296[k] = -4.0 * hi_184[k]
                   + f_0 * ki_436[k];

        t_297[k] = -4.0 * hi_185[k]
                   + f_0 * ki_437[k];

        t_298[k] = -4.0 * hi_186[k]
                   + f_0 * ki_438[k];

        t_299[k] = -4.0 * hi_187[k]
                   + f_0 * ki_439[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, hi_188, hi_189, hi_190, hi_191, \
                         hi_192, ki_440, ki_441, ki_442, ki_443, \
                         ki_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * hi_188[k]
                   + f_0 * ki_440[k];

        t_301[k] = -4.0 * hi_189[k]
                   + f_0 * ki_441[k];

        t_302[k] = -4.0 * hi_190[k]
                   + f_0 * ki_442[k];

        t_303[k] = -4.0 * hi_191[k]
                   + f_0 * ki_443[k];

        t_304[k] = -4.0 * hi_192[k]
                   + f_0 * ki_444[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, hi_193, hi_194, hi_195, hi_196, \
                         hi_197, ki_445, ki_446, ki_447, ki_448, \
                         ki_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * hi_193[k]
                   + f_0 * ki_445[k];

        t_306[k] = -4.0 * hi_194[k]
                   + f_0 * ki_446[k];

        t_307[k] = -4.0 * hi_195[k]
                   + f_0 * ki_447[k];

        t_308[k] = -3.0 * hi_196[k]
                   + f_0 * ki_448[k];

        t_309[k] = -3.0 * hi_197[k]
                   + f_0 * ki_449[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, hi_198, hi_199, hi_200, hi_201, \
                         hi_202, ki_450, ki_451, ki_452, ki_453, \
                         ki_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * hi_198[k]
                   + f_0 * ki_450[k];

        t_311[k] = -3.0 * hi_199[k]
                   + f_0 * ki_451[k];

        t_312[k] = -3.0 * hi_200[k]
                   + f_0 * ki_452[k];

        t_313[k] = -3.0 * hi_201[k]
                   + f_0 * ki_453[k];

        t_314[k] = -3.0 * hi_202[k]
                   + f_0 * ki_454[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, hi_203, hi_204, hi_205, hi_206, \
                         hi_207, ki_455, ki_456, ki_457, ki_458, \
                         ki_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -3.0 * hi_203[k]
                   + f_0 * ki_455[k];

        t_316[k] = -3.0 * hi_204[k]
                   + f_0 * ki_456[k];

        t_317[k] = -3.0 * hi_205[k]
                   + f_0 * ki_457[k];

        t_318[k] = -3.0 * hi_206[k]
                   + f_0 * ki_458[k];

        t_319[k] = -3.0 * hi_207[k]
                   + f_0 * ki_459[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, hi_208, hi_209, hi_210, hi_211, \
                         hi_212, ki_460, ki_461, ki_462, ki_463, \
                         ki_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -3.0 * hi_208[k]
                   + f_0 * ki_460[k];

        t_321[k] = -3.0 * hi_209[k]
                   + f_0 * ki_461[k];

        t_322[k] = -3.0 * hi_210[k]
                   + f_0 * ki_462[k];

        t_323[k] = -3.0 * hi_211[k]
                   + f_0 * ki_463[k];

        t_324[k] = -3.0 * hi_212[k]
                   + f_0 * ki_464[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, hi_213, hi_214, hi_215, hi_216, \
                         hi_217, ki_465, ki_466, ki_467, ki_468, \
                         ki_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -3.0 * hi_213[k]
                   + f_0 * ki_465[k];

        t_326[k] = -3.0 * hi_214[k]
                   + f_0 * ki_466[k];

        t_327[k] = -3.0 * hi_215[k]
                   + f_0 * ki_467[k];

        t_328[k] = -3.0 * hi_216[k]
                   + f_0 * ki_468[k];

        t_329[k] = -3.0 * hi_217[k]
                   + f_0 * ki_469[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, hi_218, hi_219, hi_220, hi_221, \
                         hi_222, ki_470, ki_471, ki_472, ki_473, \
                         ki_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -3.0 * hi_218[k]
                   + f_0 * ki_470[k];

        t_331[k] = -3.0 * hi_219[k]
                   + f_0 * ki_471[k];

        t_332[k] = -3.0 * hi_220[k]
                   + f_0 * ki_472[k];

        t_333[k] = -3.0 * hi_221[k]
                   + f_0 * ki_473[k];

        t_334[k] = -3.0 * hi_222[k]
                   + f_0 * ki_474[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hi_223 = buffer.data(hi + 223);
    const auto *hi_224 = buffer.data(hi + 224);
    const auto *hi_225 = buffer.data(hi + 225);
    const auto *hi_226 = buffer.data(hi + 226);
    const auto *hi_227 = buffer.data(hi + 227);
    const auto *hi_228 = buffer.data(hi + 228);
    const auto *hi_229 = buffer.data(hi + 229);
    const auto *hi_230 = buffer.data(hi + 230);
    const auto *hi_231 = buffer.data(hi + 231);
    const auto *hi_232 = buffer.data(hi + 232);
    const auto *hi_233 = buffer.data(hi + 233);
    const auto *hi_234 = buffer.data(hi + 234);
    const auto *hi_235 = buffer.data(hi + 235);
    const auto *hi_236 = buffer.data(hi + 236);
    const auto *hi_237 = buffer.data(hi + 237);
    const auto *hi_238 = buffer.data(hi + 238);
    const auto *hi_239 = buffer.data(hi + 239);
    const auto *hi_240 = buffer.data(hi + 240);
    const auto *hi_241 = buffer.data(hi + 241);
    const auto *hi_242 = buffer.data(hi + 242);
    const auto *hi_243 = buffer.data(hi + 243);
    const auto *hi_244 = buffer.data(hi + 244);
    const auto *hi_245 = buffer.data(hi + 245);
    const auto *hi_246 = buffer.data(hi + 246);
    const auto *hi_247 = buffer.data(hi + 247);
    const auto *hi_248 = buffer.data(hi + 248);
    const auto *hi_249 = buffer.data(hi + 249);
    const auto *hi_250 = buffer.data(hi + 250);
    const auto *hi_251 = buffer.data(hi + 251);
    const auto *hi_252 = buffer.data(hi + 252);
    const auto *hi_253 = buffer.data(hi + 253);
    const auto *hi_254 = buffer.data(hi + 254);
    const auto *hi_255 = buffer.data(hi + 255);
    const auto *hi_256 = buffer.data(hi + 256);
    const auto *hi_257 = buffer.data(hi + 257);
    const auto *hi_258 = buffer.data(hi + 258);
    const auto *hi_259 = buffer.data(hi + 259);
    const auto *hi_260 = buffer.data(hi + 260);
    const auto *hi_261 = buffer.data(hi + 261);
    const auto *hi_262 = buffer.data(hi + 262);
    const auto *hi_263 = buffer.data(hi + 263);
    const auto *hi_264 = buffer.data(hi + 264);
    const auto *hi_265 = buffer.data(hi + 265);
    const auto *hi_266 = buffer.data(hi + 266);
    const auto *hi_267 = buffer.data(hi + 267);
    const auto *hi_268 = buffer.data(hi + 268);
    const auto *hi_269 = buffer.data(hi + 269);
    const auto *hi_270 = buffer.data(hi + 270);
    const auto *hi_271 = buffer.data(hi + 271);
    const auto *hi_272 = buffer.data(hi + 272);
    const auto *hi_273 = buffer.data(hi + 273);
    const auto *hi_274 = buffer.data(hi + 274);
    const auto *hi_275 = buffer.data(hi + 275);
    const auto *hi_276 = buffer.data(hi + 276);
    const auto *hi_277 = buffer.data(hi + 277);
    const auto *hi_278 = buffer.data(hi + 278);
    const auto *hi_279 = buffer.data(hi + 279);
    const auto *hi_280 = buffer.data(hi + 280);
    const auto *hi_281 = buffer.data(hi + 281);
    const auto *hi_282 = buffer.data(hi + 282);
    const auto *hi_283 = buffer.data(hi + 283);
    const auto *hi_284 = buffer.data(hi + 284);
    const auto *hi_285 = buffer.data(hi + 285);
    const auto *hi_286 = buffer.data(hi + 286);
    const auto *hi_287 = buffer.data(hi + 287);
    const auto *hi_288 = buffer.data(hi + 288);
    const auto *hi_289 = buffer.data(hi + 289);
    const auto *hi_290 = buffer.data(hi + 290);
    const auto *hi_291 = buffer.data(hi + 291);
    const auto *hi_292 = buffer.data(hi + 292);
    const auto *hi_293 = buffer.data(hi + 293);
    const auto *hi_294 = buffer.data(hi + 294);
    const auto *hi_295 = buffer.data(hi + 295);
    const auto *hi_296 = buffer.data(hi + 296);
    const auto *hi_297 = buffer.data(hi + 297);
    const auto *hi_298 = buffer.data(hi + 298);
    const auto *hi_299 = buffer.data(hi + 299);
    const auto *hi_300 = buffer.data(hi + 300);
    const auto *hi_301 = buffer.data(hi + 301);
    const auto *hi_302 = buffer.data(hi + 302);
    const auto *hi_303 = buffer.data(hi + 303);
    const auto *hi_304 = buffer.data(hi + 304);
    const auto *hi_305 = buffer.data(hi + 305);
    const auto *hi_306 = buffer.data(hi + 306);
    const auto *hi_307 = buffer.data(hi + 307);
    const auto *hi_308 = buffer.data(hi + 308);
    const auto *hi_309 = buffer.data(hi + 309);
    const auto *hi_310 = buffer.data(hi + 310);
    const auto *hi_311 = buffer.data(hi + 311);
    const auto *hi_312 = buffer.data(hi + 312);
    const auto *hi_313 = buffer.data(hi + 313);
    const auto *hi_314 = buffer.data(hi + 314);
    const auto *hi_315 = buffer.data(hi + 315);
    const auto *hi_316 = buffer.data(hi + 316);
    const auto *hi_317 = buffer.data(hi + 317);
    const auto *hi_318 = buffer.data(hi + 318);
    const auto *hi_319 = buffer.data(hi + 319);
    const auto *hi_320 = buffer.data(hi + 320);
    const auto *hi_321 = buffer.data(hi + 321);
    const auto *hi_322 = buffer.data(hi + 322);
    const auto *hi_323 = buffer.data(hi + 323);
    const auto *hi_324 = buffer.data(hi + 324);
    const auto *hi_325 = buffer.data(hi + 325);
    const auto *hi_326 = buffer.data(hi + 326);
    const auto *hi_327 = buffer.data(hi + 327);
    const auto *hi_328 = buffer.data(hi + 328);
    const auto *hi_329 = buffer.data(hi + 329);
    const auto *hi_330 = buffer.data(hi + 330);
    const auto *hi_331 = buffer.data(hi + 331);
    const auto *hi_332 = buffer.data(hi + 332);
    const auto *hi_333 = buffer.data(hi + 333);
    const auto *hi_334 = buffer.data(hi + 334);
    const auto *hi_335 = buffer.data(hi + 335);
    const auto *hi_336 = buffer.data(hi + 336);
    const auto *hi_337 = buffer.data(hi + 337);
    const auto *hi_338 = buffer.data(hi + 338);
    const auto *hi_339 = buffer.data(hi + 339);
    const auto *hi_340 = buffer.data(hi + 340);
    const auto *hi_341 = buffer.data(hi + 341);
    const auto *hi_342 = buffer.data(hi + 342);
    const auto *hi_343 = buffer.data(hi + 343);
    const auto *hi_344 = buffer.data(hi + 344);
    const auto *hi_345 = buffer.data(hi + 345);
    const auto *hi_346 = buffer.data(hi + 346);
    const auto *hi_347 = buffer.data(hi + 347);
    const auto *hi_348 = buffer.data(hi + 348);
    const auto *hi_349 = buffer.data(hi + 349);
    const auto *hi_350 = buffer.data(hi + 350);
    const auto *hi_351 = buffer.data(hi + 351);
    const auto *hi_352 = buffer.data(hi + 352);
    const auto *hi_353 = buffer.data(hi + 353);
    const auto *hi_354 = buffer.data(hi + 354);

    const auto *ki_475 = buffer.data(ki + 475);
    const auto *ki_476 = buffer.data(ki + 476);
    const auto *ki_477 = buffer.data(ki + 477);
    const auto *ki_478 = buffer.data(ki + 478);
    const auto *ki_479 = buffer.data(ki + 479);
    const auto *ki_480 = buffer.data(ki + 480);
    const auto *ki_481 = buffer.data(ki + 481);
    const auto *ki_482 = buffer.data(ki + 482);
    const auto *ki_483 = buffer.data(ki + 483);
    const auto *ki_484 = buffer.data(ki + 484);
    const auto *ki_485 = buffer.data(ki + 485);
    const auto *ki_486 = buffer.data(ki + 486);
    const auto *ki_487 = buffer.data(ki + 487);
    const auto *ki_488 = buffer.data(ki + 488);
    const auto *ki_489 = buffer.data(ki + 489);
    const auto *ki_490 = buffer.data(ki + 490);
    const auto *ki_491 = buffer.data(ki + 491);
    const auto *ki_492 = buffer.data(ki + 492);
    const auto *ki_493 = buffer.data(ki + 493);
    const auto *ki_494 = buffer.data(ki + 494);
    const auto *ki_495 = buffer.data(ki + 495);
    const auto *ki_496 = buffer.data(ki + 496);
    const auto *ki_497 = buffer.data(ki + 497);
    const auto *ki_498 = buffer.data(ki + 498);
    const auto *ki_499 = buffer.data(ki + 499);
    const auto *ki_500 = buffer.data(ki + 500);
    const auto *ki_501 = buffer.data(ki + 501);
    const auto *ki_502 = buffer.data(ki + 502);
    const auto *ki_503 = buffer.data(ki + 503);
    const auto *ki_504 = buffer.data(ki + 504);
    const auto *ki_505 = buffer.data(ki + 505);
    const auto *ki_506 = buffer.data(ki + 506);
    const auto *ki_507 = buffer.data(ki + 507);
    const auto *ki_508 = buffer.data(ki + 508);
    const auto *ki_509 = buffer.data(ki + 509);
    const auto *ki_510 = buffer.data(ki + 510);
    const auto *ki_511 = buffer.data(ki + 511);
    const auto *ki_512 = buffer.data(ki + 512);
    const auto *ki_513 = buffer.data(ki + 513);
    const auto *ki_514 = buffer.data(ki + 514);
    const auto *ki_515 = buffer.data(ki + 515);
    const auto *ki_516 = buffer.data(ki + 516);
    const auto *ki_517 = buffer.data(ki + 517);
    const auto *ki_518 = buffer.data(ki + 518);
    const auto *ki_519 = buffer.data(ki + 519);
    const auto *ki_520 = buffer.data(ki + 520);
    const auto *ki_521 = buffer.data(ki + 521);
    const auto *ki_522 = buffer.data(ki + 522);
    const auto *ki_523 = buffer.data(ki + 523);
    const auto *ki_524 = buffer.data(ki + 524);
    const auto *ki_525 = buffer.data(ki + 525);
    const auto *ki_526 = buffer.data(ki + 526);
    const auto *ki_527 = buffer.data(ki + 527);
    const auto *ki_528 = buffer.data(ki + 528);
    const auto *ki_529 = buffer.data(ki + 529);
    const auto *ki_530 = buffer.data(ki + 530);
    const auto *ki_531 = buffer.data(ki + 531);
    const auto *ki_532 = buffer.data(ki + 532);
    const auto *ki_533 = buffer.data(ki + 533);
    const auto *ki_534 = buffer.data(ki + 534);
    const auto *ki_535 = buffer.data(ki + 535);
    const auto *ki_536 = buffer.data(ki + 536);
    const auto *ki_537 = buffer.data(ki + 537);
    const auto *ki_538 = buffer.data(ki + 538);
    const auto *ki_539 = buffer.data(ki + 539);
    const auto *ki_540 = buffer.data(ki + 540);
    const auto *ki_541 = buffer.data(ki + 541);
    const auto *ki_542 = buffer.data(ki + 542);
    const auto *ki_543 = buffer.data(ki + 543);
    const auto *ki_544 = buffer.data(ki + 544);
    const auto *ki_545 = buffer.data(ki + 545);
    const auto *ki_546 = buffer.data(ki + 546);
    const auto *ki_547 = buffer.data(ki + 547);
    const auto *ki_548 = buffer.data(ki + 548);
    const auto *ki_549 = buffer.data(ki + 549);
    const auto *ki_550 = buffer.data(ki + 550);
    const auto *ki_551 = buffer.data(ki + 551);
    const auto *ki_552 = buffer.data(ki + 552);
    const auto *ki_553 = buffer.data(ki + 553);
    const auto *ki_554 = buffer.data(ki + 554);
    const auto *ki_555 = buffer.data(ki + 555);
    const auto *ki_556 = buffer.data(ki + 556);
    const auto *ki_557 = buffer.data(ki + 557);
    const auto *ki_558 = buffer.data(ki + 558);
    const auto *ki_559 = buffer.data(ki + 559);
    const auto *ki_588 = buffer.data(ki + 588);
    const auto *ki_589 = buffer.data(ki + 589);
    const auto *ki_590 = buffer.data(ki + 590);
    const auto *ki_591 = buffer.data(ki + 591);
    const auto *ki_592 = buffer.data(ki + 592);
    const auto *ki_593 = buffer.data(ki + 593);
    const auto *ki_594 = buffer.data(ki + 594);
    const auto *ki_595 = buffer.data(ki + 595);
    const auto *ki_596 = buffer.data(ki + 596);
    const auto *ki_597 = buffer.data(ki + 597);
    const auto *ki_598 = buffer.data(ki + 598);
    const auto *ki_599 = buffer.data(ki + 599);
    const auto *ki_600 = buffer.data(ki + 600);
    const auto *ki_601 = buffer.data(ki + 601);
    const auto *ki_602 = buffer.data(ki + 602);
    const auto *ki_603 = buffer.data(ki + 603);
    const auto *ki_604 = buffer.data(ki + 604);
    const auto *ki_605 = buffer.data(ki + 605);
    const auto *ki_606 = buffer.data(ki + 606);
    const auto *ki_607 = buffer.data(ki + 607);
    const auto *ki_608 = buffer.data(ki + 608);
    const auto *ki_609 = buffer.data(ki + 609);
    const auto *ki_610 = buffer.data(ki + 610);
    const auto *ki_611 = buffer.data(ki + 611);
    const auto *ki_612 = buffer.data(ki + 612);
    const auto *ki_613 = buffer.data(ki + 613);
    const auto *ki_614 = buffer.data(ki + 614);
    const auto *ki_615 = buffer.data(ki + 615);
    const auto *ki_616 = buffer.data(ki + 616);
    const auto *ki_617 = buffer.data(ki + 617);
    const auto *ki_618 = buffer.data(ki + 618);
    const auto *ki_619 = buffer.data(ki + 619);
    const auto *ki_620 = buffer.data(ki + 620);
    const auto *ki_621 = buffer.data(ki + 621);
    const auto *ki_622 = buffer.data(ki + 622);
    const auto *ki_623 = buffer.data(ki + 623);
    const auto *ki_624 = buffer.data(ki + 624);
    const auto *ki_625 = buffer.data(ki + 625);
    const auto *ki_626 = buffer.data(ki + 626);
    const auto *ki_627 = buffer.data(ki + 627);
    const auto *ki_628 = buffer.data(ki + 628);
    const auto *ki_629 = buffer.data(ki + 629);
    const auto *ki_630 = buffer.data(ki + 630);
    const auto *ki_631 = buffer.data(ki + 631);
    const auto *ki_632 = buffer.data(ki + 632);
    const auto *ki_633 = buffer.data(ki + 633);
    const auto *ki_634 = buffer.data(ki + 634);
    const auto *ki_635 = buffer.data(ki + 635);
    const auto *ki_636 = buffer.data(ki + 636);
    const auto *ki_637 = buffer.data(ki + 637);
    const auto *ki_638 = buffer.data(ki + 638);
    const auto *ki_639 = buffer.data(ki + 639);
    const auto *ki_640 = buffer.data(ki + 640);
    const auto *ki_641 = buffer.data(ki + 641);
    const auto *ki_642 = buffer.data(ki + 642);
    const auto *ki_643 = buffer.data(ki + 643);
    const auto *ki_644 = buffer.data(ki + 644);
    const auto *ki_645 = buffer.data(ki + 645);
    const auto *ki_646 = buffer.data(ki + 646);
    const auto *ki_647 = buffer.data(ki + 647);
    const auto *ki_648 = buffer.data(ki + 648);
    const auto *ki_649 = buffer.data(ki + 649);
    const auto *ki_650 = buffer.data(ki + 650);
    const auto *ki_651 = buffer.data(ki + 651);
    const auto *ki_652 = buffer.data(ki + 652);
    const auto *ki_653 = buffer.data(ki + 653);
    const auto *ki_654 = buffer.data(ki + 654);
    const auto *ki_655 = buffer.data(ki + 655);
    const auto *ki_656 = buffer.data(ki + 656);
    const auto *ki_657 = buffer.data(ki + 657);
    const auto *ki_658 = buffer.data(ki + 658);
    const auto *ki_659 = buffer.data(ki + 659);
    const auto *ki_660 = buffer.data(ki + 660);
    const auto *ki_661 = buffer.data(ki + 661);
    const auto *ki_662 = buffer.data(ki + 662);

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, hi_223, hi_224, hi_225, hi_226, \
                         hi_227, ki_475, ki_476, ki_477, ki_478, \
                         ki_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -3.0 * hi_223[k]
                   + f_0 * ki_475[k];

        t_336[k] = -2.0 * hi_224[k]
                   + f_0 * ki_476[k];

        t_337[k] = -2.0 * hi_225[k]
                   + f_0 * ki_477[k];

        t_338[k] = -2.0 * hi_226[k]
                   + f_0 * ki_478[k];

        t_339[k] = -2.0 * hi_227[k]
                   + f_0 * ki_479[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, hi_228, hi_229, hi_230, hi_231, \
                         hi_232, ki_480, ki_481, ki_482, ki_483, \
                         ki_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * hi_228[k]
                   + f_0 * ki_480[k];

        t_341[k] = -2.0 * hi_229[k]
                   + f_0 * ki_481[k];

        t_342[k] = -2.0 * hi_230[k]
                   + f_0 * ki_482[k];

        t_343[k] = -2.0 * hi_231[k]
                   + f_0 * ki_483[k];

        t_344[k] = -2.0 * hi_232[k]
                   + f_0 * ki_484[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, hi_233, hi_234, hi_235, hi_236, \
                         hi_237, ki_485, ki_486, ki_487, ki_488, \
                         ki_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * hi_233[k]
                   + f_0 * ki_485[k];

        t_346[k] = -2.0 * hi_234[k]
                   + f_0 * ki_486[k];

        t_347[k] = -2.0 * hi_235[k]
                   + f_0 * ki_487[k];

        t_348[k] = -2.0 * hi_236[k]
                   + f_0 * ki_488[k];

        t_349[k] = -2.0 * hi_237[k]
                   + f_0 * ki_489[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, hi_238, hi_239, hi_240, hi_241, \
                         hi_242, ki_490, ki_491, ki_492, ki_493, \
                         ki_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * hi_238[k]
                   + f_0 * ki_490[k];

        t_351[k] = -2.0 * hi_239[k]
                   + f_0 * ki_491[k];

        t_352[k] = -2.0 * hi_240[k]
                   + f_0 * ki_492[k];

        t_353[k] = -2.0 * hi_241[k]
                   + f_0 * ki_493[k];

        t_354[k] = -2.0 * hi_242[k]
                   + f_0 * ki_494[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, hi_243, hi_244, hi_245, hi_246, \
                         hi_247, ki_495, ki_496, ki_497, ki_498, \
                         ki_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * hi_243[k]
                   + f_0 * ki_495[k];

        t_356[k] = -2.0 * hi_244[k]
                   + f_0 * ki_496[k];

        t_357[k] = -2.0 * hi_245[k]
                   + f_0 * ki_497[k];

        t_358[k] = -2.0 * hi_246[k]
                   + f_0 * ki_498[k];

        t_359[k] = -2.0 * hi_247[k]
                   + f_0 * ki_499[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, hi_248, hi_249, hi_250, hi_251, \
                         hi_252, ki_500, ki_501, ki_502, ki_503, \
                         ki_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -2.0 * hi_248[k]
                   + f_0 * ki_500[k];

        t_361[k] = -2.0 * hi_249[k]
                   + f_0 * ki_501[k];

        t_362[k] = -2.0 * hi_250[k]
                   + f_0 * ki_502[k];

        t_363[k] = -2.0 * hi_251[k]
                   + f_0 * ki_503[k];

        t_364[k] = -hi_252[k]
                   + f_0 * ki_504[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, hi_253, hi_254, hi_255, hi_256, \
                         hi_257, ki_505, ki_506, ki_507, ki_508, \
                         ki_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -hi_253[k]
                   + f_0 * ki_505[k];

        t_366[k] = -hi_254[k]
                   + f_0 * ki_506[k];

        t_367[k] = -hi_255[k]
                   + f_0 * ki_507[k];

        t_368[k] = -hi_256[k]
                   + f_0 * ki_508[k];

        t_369[k] = -hi_257[k]
                   + f_0 * ki_509[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, hi_258, hi_259, hi_260, hi_261, \
                         hi_262, ki_510, ki_511, ki_512, ki_513, \
                         ki_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -hi_258[k]
                   + f_0 * ki_510[k];

        t_371[k] = -hi_259[k]
                   + f_0 * ki_511[k];

        t_372[k] = -hi_260[k]
                   + f_0 * ki_512[k];

        t_373[k] = -hi_261[k]
                   + f_0 * ki_513[k];

        t_374[k] = -hi_262[k]
                   + f_0 * ki_514[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, hi_263, hi_264, hi_265, hi_266, \
                         hi_267, ki_515, ki_516, ki_517, ki_518, \
                         ki_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -hi_263[k]
                   + f_0 * ki_515[k];

        t_376[k] = -hi_264[k]
                   + f_0 * ki_516[k];

        t_377[k] = -hi_265[k]
                   + f_0 * ki_517[k];

        t_378[k] = -hi_266[k]
                   + f_0 * ki_518[k];

        t_379[k] = -hi_267[k]
                   + f_0 * ki_519[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, hi_268, hi_269, hi_270, hi_271, \
                         hi_272, ki_520, ki_521, ki_522, ki_523, \
                         ki_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -hi_268[k]
                   + f_0 * ki_520[k];

        t_381[k] = -hi_269[k]
                   + f_0 * ki_521[k];

        t_382[k] = -hi_270[k]
                   + f_0 * ki_522[k];

        t_383[k] = -hi_271[k]
                   + f_0 * ki_523[k];

        t_384[k] = -hi_272[k]
                   + f_0 * ki_524[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, hi_273, hi_274, hi_275, hi_276, \
                         hi_277, ki_525, ki_526, ki_527, ki_528, \
                         ki_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -hi_273[k]
                   + f_0 * ki_525[k];

        t_386[k] = -hi_274[k]
                   + f_0 * ki_526[k];

        t_387[k] = -hi_275[k]
                   + f_0 * ki_527[k];

        t_388[k] = -hi_276[k]
                   + f_0 * ki_528[k];

        t_389[k] = -hi_277[k]
                   + f_0 * ki_529[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, t_395, t_396, hi_278, hi_279, \
                         ki_530, ki_531, ki_532, ki_533, ki_534, ki_535, \
                         ki_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -hi_278[k]
                   + f_0 * ki_530[k];

        t_391[k] = -hi_279[k]
                   + f_0 * ki_531[k];

        t_392[k] = f_0 * ki_532[k];

        t_393[k] = f_0 * ki_533[k];

        t_394[k] = f_0 * ki_534[k];

        t_395[k] = f_0 * ki_535[k];

        t_396[k] = f_0 * ki_536[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, t_402, t_403, t_404, ki_537, \
                         ki_538, ki_539, ki_540, ki_541, ki_542, ki_543, \
                         ki_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_0 * ki_537[k];

        t_398[k] = f_0 * ki_538[k];

        t_399[k] = f_0 * ki_539[k];

        t_400[k] = f_0 * ki_540[k];

        t_401[k] = f_0 * ki_541[k];

        t_402[k] = f_0 * ki_542[k];

        t_403[k] = f_0 * ki_543[k];

        t_404[k] = f_0 * ki_544[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, ki_545, \
                         ki_546, ki_547, ki_548, ki_549, ki_550, ki_551, \
                         ki_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * ki_545[k];

        t_406[k] = f_0 * ki_546[k];

        t_407[k] = f_0 * ki_547[k];

        t_408[k] = f_0 * ki_548[k];

        t_409[k] = f_0 * ki_549[k];

        t_410[k] = f_0 * ki_550[k];

        t_411[k] = f_0 * ki_551[k];

        t_412[k] = f_0 * ki_552[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, ki_553, ki_554, \
                         ki_555, ki_556, ki_557, ki_558, ki_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * ki_553[k];

        t_414[k] = f_0 * ki_554[k];

        t_415[k] = f_0 * ki_555[k];

        t_416[k] = f_0 * ki_556[k];

        t_417[k] = f_0 * ki_557[k];

        t_418[k] = f_0 * ki_558[k];

        t_419[k] = f_0 * ki_559[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, hi_280, hi_281, hi_282, hi_283, \
                         hi_284, ki_588, ki_589, ki_590, ki_591, \
                         ki_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -5.0 * hi_280[k]
                   + f_0 * ki_588[k];

        t_421[k] = -5.0 * hi_281[k]
                   + f_0 * ki_589[k];

        t_422[k] = -5.0 * hi_282[k]
                   + f_0 * ki_590[k];

        t_423[k] = -5.0 * hi_283[k]
                   + f_0 * ki_591[k];

        t_424[k] = -5.0 * hi_284[k]
                   + f_0 * ki_592[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, hi_285, hi_286, hi_287, hi_288, \
                         hi_289, ki_593, ki_594, ki_595, ki_596, \
                         ki_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -5.0 * hi_285[k]
                   + f_0 * ki_593[k];

        t_426[k] = -5.0 * hi_286[k]
                   + f_0 * ki_594[k];

        t_427[k] = -5.0 * hi_287[k]
                   + f_0 * ki_595[k];

        t_428[k] = -5.0 * hi_288[k]
                   + f_0 * ki_596[k];

        t_429[k] = -5.0 * hi_289[k]
                   + f_0 * ki_597[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, hi_290, hi_291, hi_292, hi_293, \
                         hi_294, ki_598, ki_599, ki_600, ki_601, \
                         ki_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -5.0 * hi_290[k]
                   + f_0 * ki_598[k];

        t_431[k] = -5.0 * hi_291[k]
                   + f_0 * ki_599[k];

        t_432[k] = -5.0 * hi_292[k]
                   + f_0 * ki_600[k];

        t_433[k] = -5.0 * hi_293[k]
                   + f_0 * ki_601[k];

        t_434[k] = -5.0 * hi_294[k]
                   + f_0 * ki_602[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, hi_295, hi_296, hi_297, hi_298, \
                         hi_299, ki_603, ki_604, ki_605, ki_606, \
                         ki_607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -5.0 * hi_295[k]
                   + f_0 * ki_603[k];

        t_436[k] = -5.0 * hi_296[k]
                   + f_0 * ki_604[k];

        t_437[k] = -5.0 * hi_297[k]
                   + f_0 * ki_605[k];

        t_438[k] = -5.0 * hi_298[k]
                   + f_0 * ki_606[k];

        t_439[k] = -5.0 * hi_299[k]
                   + f_0 * ki_607[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, hi_300, hi_301, hi_302, hi_303, \
                         hi_304, ki_608, ki_609, ki_610, ki_611, \
                         ki_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -5.0 * hi_300[k]
                   + f_0 * ki_608[k];

        t_441[k] = -5.0 * hi_301[k]
                   + f_0 * ki_609[k];

        t_442[k] = -5.0 * hi_302[k]
                   + f_0 * ki_610[k];

        t_443[k] = -5.0 * hi_303[k]
                   + f_0 * ki_611[k];

        t_444[k] = -5.0 * hi_304[k]
                   + f_0 * ki_612[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, hi_305, hi_306, hi_307, hi_308, \
                         hi_309, ki_613, ki_614, ki_615, ki_616, \
                         ki_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -5.0 * hi_305[k]
                   + f_0 * ki_613[k];

        t_446[k] = -5.0 * hi_306[k]
                   + f_0 * ki_614[k];

        t_447[k] = -5.0 * hi_307[k]
                   + f_0 * ki_615[k];

        t_448[k] = -4.0 * hi_308[k]
                   + f_0 * ki_616[k];

        t_449[k] = -4.0 * hi_309[k]
                   + f_0 * ki_617[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, hi_310, hi_311, hi_312, hi_313, \
                         hi_314, ki_618, ki_619, ki_620, ki_621, \
                         ki_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -4.0 * hi_310[k]
                   + f_0 * ki_618[k];

        t_451[k] = -4.0 * hi_311[k]
                   + f_0 * ki_619[k];

        t_452[k] = -4.0 * hi_312[k]
                   + f_0 * ki_620[k];

        t_453[k] = -4.0 * hi_313[k]
                   + f_0 * ki_621[k];

        t_454[k] = -4.0 * hi_314[k]
                   + f_0 * ki_622[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, hi_315, hi_316, hi_317, hi_318, \
                         hi_319, ki_623, ki_624, ki_625, ki_626, \
                         ki_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -4.0 * hi_315[k]
                   + f_0 * ki_623[k];

        t_456[k] = -4.0 * hi_316[k]
                   + f_0 * ki_624[k];

        t_457[k] = -4.0 * hi_317[k]
                   + f_0 * ki_625[k];

        t_458[k] = -4.0 * hi_318[k]
                   + f_0 * ki_626[k];

        t_459[k] = -4.0 * hi_319[k]
                   + f_0 * ki_627[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, hi_320, hi_321, hi_322, hi_323, \
                         hi_324, ki_628, ki_629, ki_630, ki_631, \
                         ki_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -4.0 * hi_320[k]
                   + f_0 * ki_628[k];

        t_461[k] = -4.0 * hi_321[k]
                   + f_0 * ki_629[k];

        t_462[k] = -4.0 * hi_322[k]
                   + f_0 * ki_630[k];

        t_463[k] = -4.0 * hi_323[k]
                   + f_0 * ki_631[k];

        t_464[k] = -4.0 * hi_324[k]
                   + f_0 * ki_632[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, hi_325, hi_326, hi_327, hi_328, \
                         hi_329, ki_633, ki_634, ki_635, ki_636, \
                         ki_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -4.0 * hi_325[k]
                   + f_0 * ki_633[k];

        t_466[k] = -4.0 * hi_326[k]
                   + f_0 * ki_634[k];

        t_467[k] = -4.0 * hi_327[k]
                   + f_0 * ki_635[k];

        t_468[k] = -4.0 * hi_328[k]
                   + f_0 * ki_636[k];

        t_469[k] = -4.0 * hi_329[k]
                   + f_0 * ki_637[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, hi_330, hi_331, hi_332, hi_333, \
                         hi_334, ki_638, ki_639, ki_640, ki_641, \
                         ki_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -4.0 * hi_330[k]
                   + f_0 * ki_638[k];

        t_471[k] = -4.0 * hi_331[k]
                   + f_0 * ki_639[k];

        t_472[k] = -4.0 * hi_332[k]
                   + f_0 * ki_640[k];

        t_473[k] = -4.0 * hi_333[k]
                   + f_0 * ki_641[k];

        t_474[k] = -4.0 * hi_334[k]
                   + f_0 * ki_642[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, hi_335, hi_336, hi_337, hi_338, \
                         hi_339, ki_643, ki_644, ki_645, ki_646, \
                         ki_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -4.0 * hi_335[k]
                   + f_0 * ki_643[k];

        t_476[k] = -3.0 * hi_336[k]
                   + f_0 * ki_644[k];

        t_477[k] = -3.0 * hi_337[k]
                   + f_0 * ki_645[k];

        t_478[k] = -3.0 * hi_338[k]
                   + f_0 * ki_646[k];

        t_479[k] = -3.0 * hi_339[k]
                   + f_0 * ki_647[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, hi_340, hi_341, hi_342, hi_343, \
                         hi_344, ki_648, ki_649, ki_650, ki_651, \
                         ki_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -3.0 * hi_340[k]
                   + f_0 * ki_648[k];

        t_481[k] = -3.0 * hi_341[k]
                   + f_0 * ki_649[k];

        t_482[k] = -3.0 * hi_342[k]
                   + f_0 * ki_650[k];

        t_483[k] = -3.0 * hi_343[k]
                   + f_0 * ki_651[k];

        t_484[k] = -3.0 * hi_344[k]
                   + f_0 * ki_652[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, hi_345, hi_346, hi_347, hi_348, \
                         hi_349, ki_653, ki_654, ki_655, ki_656, \
                         ki_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -3.0 * hi_345[k]
                   + f_0 * ki_653[k];

        t_486[k] = -3.0 * hi_346[k]
                   + f_0 * ki_654[k];

        t_487[k] = -3.0 * hi_347[k]
                   + f_0 * ki_655[k];

        t_488[k] = -3.0 * hi_348[k]
                   + f_0 * ki_656[k];

        t_489[k] = -3.0 * hi_349[k]
                   + f_0 * ki_657[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, hi_350, hi_351, hi_352, hi_353, \
                         hi_354, ki_658, ki_659, ki_660, ki_661, \
                         ki_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -3.0 * hi_350[k]
                   + f_0 * ki_658[k];

        t_491[k] = -3.0 * hi_351[k]
                   + f_0 * ki_659[k];

        t_492[k] = -3.0 * hi_352[k]
                   + f_0 * ki_660[k];

        t_493[k] = -3.0 * hi_353[k]
                   + f_0 * ki_661[k];

        t_494[k] = -3.0 * hi_354[k]
                   + f_0 * ki_662[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hi_355 = buffer.data(hi + 355);
    const auto *hi_356 = buffer.data(hi + 356);
    const auto *hi_357 = buffer.data(hi + 357);
    const auto *hi_358 = buffer.data(hi + 358);
    const auto *hi_359 = buffer.data(hi + 359);
    const auto *hi_360 = buffer.data(hi + 360);
    const auto *hi_361 = buffer.data(hi + 361);
    const auto *hi_362 = buffer.data(hi + 362);
    const auto *hi_363 = buffer.data(hi + 363);
    const auto *hi_364 = buffer.data(hi + 364);
    const auto *hi_365 = buffer.data(hi + 365);
    const auto *hi_366 = buffer.data(hi + 366);
    const auto *hi_367 = buffer.data(hi + 367);
    const auto *hi_368 = buffer.data(hi + 368);
    const auto *hi_369 = buffer.data(hi + 369);
    const auto *hi_370 = buffer.data(hi + 370);
    const auto *hi_371 = buffer.data(hi + 371);
    const auto *hi_372 = buffer.data(hi + 372);
    const auto *hi_373 = buffer.data(hi + 373);
    const auto *hi_374 = buffer.data(hi + 374);
    const auto *hi_375 = buffer.data(hi + 375);
    const auto *hi_376 = buffer.data(hi + 376);
    const auto *hi_377 = buffer.data(hi + 377);
    const auto *hi_378 = buffer.data(hi + 378);
    const auto *hi_379 = buffer.data(hi + 379);
    const auto *hi_380 = buffer.data(hi + 380);
    const auto *hi_381 = buffer.data(hi + 381);
    const auto *hi_382 = buffer.data(hi + 382);
    const auto *hi_383 = buffer.data(hi + 383);
    const auto *hi_384 = buffer.data(hi + 384);
    const auto *hi_385 = buffer.data(hi + 385);
    const auto *hi_386 = buffer.data(hi + 386);
    const auto *hi_387 = buffer.data(hi + 387);
    const auto *hi_388 = buffer.data(hi + 388);
    const auto *hi_389 = buffer.data(hi + 389);
    const auto *hi_390 = buffer.data(hi + 390);
    const auto *hi_391 = buffer.data(hi + 391);
    const auto *hi_392 = buffer.data(hi + 392);
    const auto *hi_393 = buffer.data(hi + 393);
    const auto *hi_394 = buffer.data(hi + 394);
    const auto *hi_395 = buffer.data(hi + 395);
    const auto *hi_396 = buffer.data(hi + 396);
    const auto *hi_397 = buffer.data(hi + 397);
    const auto *hi_398 = buffer.data(hi + 398);
    const auto *hi_399 = buffer.data(hi + 399);
    const auto *hi_400 = buffer.data(hi + 400);
    const auto *hi_401 = buffer.data(hi + 401);
    const auto *hi_402 = buffer.data(hi + 402);
    const auto *hi_403 = buffer.data(hi + 403);
    const auto *hi_404 = buffer.data(hi + 404);
    const auto *hi_405 = buffer.data(hi + 405);
    const auto *hi_406 = buffer.data(hi + 406);
    const auto *hi_407 = buffer.data(hi + 407);
    const auto *hi_408 = buffer.data(hi + 408);
    const auto *hi_409 = buffer.data(hi + 409);
    const auto *hi_410 = buffer.data(hi + 410);
    const auto *hi_411 = buffer.data(hi + 411);
    const auto *hi_412 = buffer.data(hi + 412);
    const auto *hi_413 = buffer.data(hi + 413);
    const auto *hi_414 = buffer.data(hi + 414);
    const auto *hi_415 = buffer.data(hi + 415);
    const auto *hi_416 = buffer.data(hi + 416);
    const auto *hi_417 = buffer.data(hi + 417);
    const auto *hi_418 = buffer.data(hi + 418);
    const auto *hi_419 = buffer.data(hi + 419);
    const auto *hi_420 = buffer.data(hi + 420);
    const auto *hi_421 = buffer.data(hi + 421);
    const auto *hi_422 = buffer.data(hi + 422);
    const auto *hi_423 = buffer.data(hi + 423);
    const auto *hi_424 = buffer.data(hi + 424);
    const auto *hi_425 = buffer.data(hi + 425);
    const auto *hi_426 = buffer.data(hi + 426);
    const auto *hi_427 = buffer.data(hi + 427);
    const auto *hi_428 = buffer.data(hi + 428);
    const auto *hi_429 = buffer.data(hi + 429);
    const auto *hi_430 = buffer.data(hi + 430);
    const auto *hi_431 = buffer.data(hi + 431);
    const auto *hi_432 = buffer.data(hi + 432);
    const auto *hi_433 = buffer.data(hi + 433);
    const auto *hi_434 = buffer.data(hi + 434);
    const auto *hi_435 = buffer.data(hi + 435);
    const auto *hi_436 = buffer.data(hi + 436);
    const auto *hi_437 = buffer.data(hi + 437);
    const auto *hi_438 = buffer.data(hi + 438);
    const auto *hi_439 = buffer.data(hi + 439);
    const auto *hi_440 = buffer.data(hi + 440);
    const auto *hi_441 = buffer.data(hi + 441);
    const auto *hi_442 = buffer.data(hi + 442);
    const auto *hi_443 = buffer.data(hi + 443);
    const auto *hi_444 = buffer.data(hi + 444);
    const auto *hi_445 = buffer.data(hi + 445);
    const auto *hi_446 = buffer.data(hi + 446);
    const auto *hi_447 = buffer.data(hi + 447);
    const auto *hi_448 = buffer.data(hi + 448);
    const auto *hi_449 = buffer.data(hi + 449);
    const auto *hi_450 = buffer.data(hi + 450);
    const auto *hi_451 = buffer.data(hi + 451);
    const auto *hi_452 = buffer.data(hi + 452);
    const auto *hi_453 = buffer.data(hi + 453);
    const auto *hi_454 = buffer.data(hi + 454);
    const auto *hi_455 = buffer.data(hi + 455);
    const auto *hi_456 = buffer.data(hi + 456);
    const auto *hi_457 = buffer.data(hi + 457);
    const auto *hi_458 = buffer.data(hi + 458);
    const auto *hi_459 = buffer.data(hi + 459);
    const auto *hi_460 = buffer.data(hi + 460);
    const auto *hi_461 = buffer.data(hi + 461);
    const auto *hi_462 = buffer.data(hi + 462);
    const auto *hi_463 = buffer.data(hi + 463);
    const auto *hi_464 = buffer.data(hi + 464);
    const auto *hi_465 = buffer.data(hi + 465);
    const auto *hi_466 = buffer.data(hi + 466);
    const auto *hi_467 = buffer.data(hi + 467);
    const auto *hi_468 = buffer.data(hi + 468);
    const auto *hi_469 = buffer.data(hi + 469);
    const auto *hi_470 = buffer.data(hi + 470);
    const auto *hi_471 = buffer.data(hi + 471);
    const auto *hi_472 = buffer.data(hi + 472);
    const auto *hi_473 = buffer.data(hi + 473);
    const auto *hi_474 = buffer.data(hi + 474);
    const auto *hi_475 = buffer.data(hi + 475);
    const auto *hi_476 = buffer.data(hi + 476);
    const auto *hi_477 = buffer.data(hi + 477);
    const auto *hi_478 = buffer.data(hi + 478);
    const auto *hi_479 = buffer.data(hi + 479);
    const auto *hi_480 = buffer.data(hi + 480);
    const auto *hi_481 = buffer.data(hi + 481);
    const auto *hi_482 = buffer.data(hi + 482);
    const auto *hi_483 = buffer.data(hi + 483);
    const auto *hi_484 = buffer.data(hi + 484);
    const auto *hi_485 = buffer.data(hi + 485);
    const auto *hi_486 = buffer.data(hi + 486);

    const auto *ki_663 = buffer.data(ki + 663);
    const auto *ki_664 = buffer.data(ki + 664);
    const auto *ki_665 = buffer.data(ki + 665);
    const auto *ki_666 = buffer.data(ki + 666);
    const auto *ki_667 = buffer.data(ki + 667);
    const auto *ki_668 = buffer.data(ki + 668);
    const auto *ki_669 = buffer.data(ki + 669);
    const auto *ki_670 = buffer.data(ki + 670);
    const auto *ki_671 = buffer.data(ki + 671);
    const auto *ki_672 = buffer.data(ki + 672);
    const auto *ki_673 = buffer.data(ki + 673);
    const auto *ki_674 = buffer.data(ki + 674);
    const auto *ki_675 = buffer.data(ki + 675);
    const auto *ki_676 = buffer.data(ki + 676);
    const auto *ki_677 = buffer.data(ki + 677);
    const auto *ki_678 = buffer.data(ki + 678);
    const auto *ki_679 = buffer.data(ki + 679);
    const auto *ki_680 = buffer.data(ki + 680);
    const auto *ki_681 = buffer.data(ki + 681);
    const auto *ki_682 = buffer.data(ki + 682);
    const auto *ki_683 = buffer.data(ki + 683);
    const auto *ki_684 = buffer.data(ki + 684);
    const auto *ki_685 = buffer.data(ki + 685);
    const auto *ki_686 = buffer.data(ki + 686);
    const auto *ki_687 = buffer.data(ki + 687);
    const auto *ki_688 = buffer.data(ki + 688);
    const auto *ki_689 = buffer.data(ki + 689);
    const auto *ki_690 = buffer.data(ki + 690);
    const auto *ki_691 = buffer.data(ki + 691);
    const auto *ki_692 = buffer.data(ki + 692);
    const auto *ki_693 = buffer.data(ki + 693);
    const auto *ki_694 = buffer.data(ki + 694);
    const auto *ki_695 = buffer.data(ki + 695);
    const auto *ki_696 = buffer.data(ki + 696);
    const auto *ki_697 = buffer.data(ki + 697);
    const auto *ki_698 = buffer.data(ki + 698);
    const auto *ki_699 = buffer.data(ki + 699);
    const auto *ki_700 = buffer.data(ki + 700);
    const auto *ki_701 = buffer.data(ki + 701);
    const auto *ki_702 = buffer.data(ki + 702);
    const auto *ki_703 = buffer.data(ki + 703);
    const auto *ki_704 = buffer.data(ki + 704);
    const auto *ki_705 = buffer.data(ki + 705);
    const auto *ki_706 = buffer.data(ki + 706);
    const auto *ki_707 = buffer.data(ki + 707);
    const auto *ki_708 = buffer.data(ki + 708);
    const auto *ki_709 = buffer.data(ki + 709);
    const auto *ki_710 = buffer.data(ki + 710);
    const auto *ki_711 = buffer.data(ki + 711);
    const auto *ki_712 = buffer.data(ki + 712);
    const auto *ki_713 = buffer.data(ki + 713);
    const auto *ki_714 = buffer.data(ki + 714);
    const auto *ki_715 = buffer.data(ki + 715);
    const auto *ki_716 = buffer.data(ki + 716);
    const auto *ki_717 = buffer.data(ki + 717);
    const auto *ki_718 = buffer.data(ki + 718);
    const auto *ki_719 = buffer.data(ki + 719);
    const auto *ki_720 = buffer.data(ki + 720);
    const auto *ki_721 = buffer.data(ki + 721);
    const auto *ki_722 = buffer.data(ki + 722);
    const auto *ki_723 = buffer.data(ki + 723);
    const auto *ki_724 = buffer.data(ki + 724);
    const auto *ki_725 = buffer.data(ki + 725);
    const auto *ki_726 = buffer.data(ki + 726);
    const auto *ki_727 = buffer.data(ki + 727);
    const auto *ki_728 = buffer.data(ki + 728);
    const auto *ki_729 = buffer.data(ki + 729);
    const auto *ki_730 = buffer.data(ki + 730);
    const auto *ki_731 = buffer.data(ki + 731);
    const auto *ki_732 = buffer.data(ki + 732);
    const auto *ki_733 = buffer.data(ki + 733);
    const auto *ki_734 = buffer.data(ki + 734);
    const auto *ki_735 = buffer.data(ki + 735);
    const auto *ki_736 = buffer.data(ki + 736);
    const auto *ki_737 = buffer.data(ki + 737);
    const auto *ki_738 = buffer.data(ki + 738);
    const auto *ki_739 = buffer.data(ki + 739);
    const auto *ki_740 = buffer.data(ki + 740);
    const auto *ki_741 = buffer.data(ki + 741);
    const auto *ki_742 = buffer.data(ki + 742);
    const auto *ki_743 = buffer.data(ki + 743);
    const auto *ki_744 = buffer.data(ki + 744);
    const auto *ki_745 = buffer.data(ki + 745);
    const auto *ki_746 = buffer.data(ki + 746);
    const auto *ki_747 = buffer.data(ki + 747);
    const auto *ki_748 = buffer.data(ki + 748);
    const auto *ki_749 = buffer.data(ki + 749);
    const auto *ki_750 = buffer.data(ki + 750);
    const auto *ki_751 = buffer.data(ki + 751);
    const auto *ki_752 = buffer.data(ki + 752);
    const auto *ki_753 = buffer.data(ki + 753);
    const auto *ki_754 = buffer.data(ki + 754);
    const auto *ki_755 = buffer.data(ki + 755);
    const auto *ki_784 = buffer.data(ki + 784);
    const auto *ki_785 = buffer.data(ki + 785);
    const auto *ki_786 = buffer.data(ki + 786);
    const auto *ki_787 = buffer.data(ki + 787);
    const auto *ki_788 = buffer.data(ki + 788);
    const auto *ki_789 = buffer.data(ki + 789);
    const auto *ki_790 = buffer.data(ki + 790);
    const auto *ki_791 = buffer.data(ki + 791);
    const auto *ki_792 = buffer.data(ki + 792);
    const auto *ki_793 = buffer.data(ki + 793);
    const auto *ki_794 = buffer.data(ki + 794);
    const auto *ki_795 = buffer.data(ki + 795);
    const auto *ki_796 = buffer.data(ki + 796);
    const auto *ki_797 = buffer.data(ki + 797);
    const auto *ki_798 = buffer.data(ki + 798);
    const auto *ki_799 = buffer.data(ki + 799);
    const auto *ki_800 = buffer.data(ki + 800);
    const auto *ki_801 = buffer.data(ki + 801);
    const auto *ki_802 = buffer.data(ki + 802);
    const auto *ki_803 = buffer.data(ki + 803);
    const auto *ki_804 = buffer.data(ki + 804);
    const auto *ki_805 = buffer.data(ki + 805);
    const auto *ki_806 = buffer.data(ki + 806);
    const auto *ki_807 = buffer.data(ki + 807);
    const auto *ki_808 = buffer.data(ki + 808);
    const auto *ki_809 = buffer.data(ki + 809);
    const auto *ki_810 = buffer.data(ki + 810);
    const auto *ki_811 = buffer.data(ki + 811);
    const auto *ki_812 = buffer.data(ki + 812);
    const auto *ki_813 = buffer.data(ki + 813);
    const auto *ki_814 = buffer.data(ki + 814);
    const auto *ki_815 = buffer.data(ki + 815);
    const auto *ki_816 = buffer.data(ki + 816);
    const auto *ki_817 = buffer.data(ki + 817);
    const auto *ki_818 = buffer.data(ki + 818);
    const auto *ki_819 = buffer.data(ki + 819);
    const auto *ki_820 = buffer.data(ki + 820);
    const auto *ki_821 = buffer.data(ki + 821);
    const auto *ki_822 = buffer.data(ki + 822);
    const auto *ki_823 = buffer.data(ki + 823);
    const auto *ki_824 = buffer.data(ki + 824);
    const auto *ki_825 = buffer.data(ki + 825);
    const auto *ki_826 = buffer.data(ki + 826);
    const auto *ki_827 = buffer.data(ki + 827);
    const auto *ki_828 = buffer.data(ki + 828);
    const auto *ki_829 = buffer.data(ki + 829);
    const auto *ki_830 = buffer.data(ki + 830);
    const auto *ki_831 = buffer.data(ki + 831);
    const auto *ki_832 = buffer.data(ki + 832);
    const auto *ki_833 = buffer.data(ki + 833);
    const auto *ki_834 = buffer.data(ki + 834);
    const auto *ki_835 = buffer.data(ki + 835);
    const auto *ki_836 = buffer.data(ki + 836);
    const auto *ki_837 = buffer.data(ki + 837);
    const auto *ki_838 = buffer.data(ki + 838);
    const auto *ki_839 = buffer.data(ki + 839);
    const auto *ki_840 = buffer.data(ki + 840);
    const auto *ki_841 = buffer.data(ki + 841);
    const auto *ki_842 = buffer.data(ki + 842);
    const auto *ki_843 = buffer.data(ki + 843);
    const auto *ki_844 = buffer.data(ki + 844);
    const auto *ki_845 = buffer.data(ki + 845);
    const auto *ki_846 = buffer.data(ki + 846);
    const auto *ki_847 = buffer.data(ki + 847);
    const auto *ki_848 = buffer.data(ki + 848);
    const auto *ki_849 = buffer.data(ki + 849);
    const auto *ki_850 = buffer.data(ki + 850);

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, hi_355, hi_356, hi_357, hi_358, \
                         hi_359, ki_663, ki_664, ki_665, ki_666, \
                         ki_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -3.0 * hi_355[k]
                   + f_0 * ki_663[k];

        t_496[k] = -3.0 * hi_356[k]
                   + f_0 * ki_664[k];

        t_497[k] = -3.0 * hi_357[k]
                   + f_0 * ki_665[k];

        t_498[k] = -3.0 * hi_358[k]
                   + f_0 * ki_666[k];

        t_499[k] = -3.0 * hi_359[k]
                   + f_0 * ki_667[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, hi_360, hi_361, hi_362, hi_363, \
                         hi_364, ki_668, ki_669, ki_670, ki_671, \
                         ki_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -3.0 * hi_360[k]
                   + f_0 * ki_668[k];

        t_501[k] = -3.0 * hi_361[k]
                   + f_0 * ki_669[k];

        t_502[k] = -3.0 * hi_362[k]
                   + f_0 * ki_670[k];

        t_503[k] = -3.0 * hi_363[k]
                   + f_0 * ki_671[k];

        t_504[k] = -2.0 * hi_364[k]
                   + f_0 * ki_672[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, hi_365, hi_366, hi_367, hi_368, \
                         hi_369, ki_673, ki_674, ki_675, ki_676, \
                         ki_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -2.0 * hi_365[k]
                   + f_0 * ki_673[k];

        t_506[k] = -2.0 * hi_366[k]
                   + f_0 * ki_674[k];

        t_507[k] = -2.0 * hi_367[k]
                   + f_0 * ki_675[k];

        t_508[k] = -2.0 * hi_368[k]
                   + f_0 * ki_676[k];

        t_509[k] = -2.0 * hi_369[k]
                   + f_0 * ki_677[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, hi_370, hi_371, hi_372, hi_373, \
                         hi_374, ki_678, ki_679, ki_680, ki_681, \
                         ki_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -2.0 * hi_370[k]
                   + f_0 * ki_678[k];

        t_511[k] = -2.0 * hi_371[k]
                   + f_0 * ki_679[k];

        t_512[k] = -2.0 * hi_372[k]
                   + f_0 * ki_680[k];

        t_513[k] = -2.0 * hi_373[k]
                   + f_0 * ki_681[k];

        t_514[k] = -2.0 * hi_374[k]
                   + f_0 * ki_682[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, hi_375, hi_376, hi_377, hi_378, \
                         hi_379, ki_683, ki_684, ki_685, ki_686, \
                         ki_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -2.0 * hi_375[k]
                   + f_0 * ki_683[k];

        t_516[k] = -2.0 * hi_376[k]
                   + f_0 * ki_684[k];

        t_517[k] = -2.0 * hi_377[k]
                   + f_0 * ki_685[k];

        t_518[k] = -2.0 * hi_378[k]
                   + f_0 * ki_686[k];

        t_519[k] = -2.0 * hi_379[k]
                   + f_0 * ki_687[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, hi_380, hi_381, hi_382, hi_383, \
                         hi_384, ki_688, ki_689, ki_690, ki_691, \
                         ki_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -2.0 * hi_380[k]
                   + f_0 * ki_688[k];

        t_521[k] = -2.0 * hi_381[k]
                   + f_0 * ki_689[k];

        t_522[k] = -2.0 * hi_382[k]
                   + f_0 * ki_690[k];

        t_523[k] = -2.0 * hi_383[k]
                   + f_0 * ki_691[k];

        t_524[k] = -2.0 * hi_384[k]
                   + f_0 * ki_692[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, hi_385, hi_386, hi_387, hi_388, \
                         hi_389, ki_693, ki_694, ki_695, ki_696, \
                         ki_697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -2.0 * hi_385[k]
                   + f_0 * ki_693[k];

        t_526[k] = -2.0 * hi_386[k]
                   + f_0 * ki_694[k];

        t_527[k] = -2.0 * hi_387[k]
                   + f_0 * ki_695[k];

        t_528[k] = -2.0 * hi_388[k]
                   + f_0 * ki_696[k];

        t_529[k] = -2.0 * hi_389[k]
                   + f_0 * ki_697[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, hi_390, hi_391, hi_392, hi_393, \
                         hi_394, ki_698, ki_699, ki_700, ki_701, \
                         ki_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -2.0 * hi_390[k]
                   + f_0 * ki_698[k];

        t_531[k] = -2.0 * hi_391[k]
                   + f_0 * ki_699[k];

        t_532[k] = -hi_392[k]
                   + f_0 * ki_700[k];

        t_533[k] = -hi_393[k]
                   + f_0 * ki_701[k];

        t_534[k] = -hi_394[k]
                   + f_0 * ki_702[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, hi_395, hi_396, hi_397, hi_398, \
                         hi_399, ki_703, ki_704, ki_705, ki_706, \
                         ki_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -hi_395[k]
                   + f_0 * ki_703[k];

        t_536[k] = -hi_396[k]
                   + f_0 * ki_704[k];

        t_537[k] = -hi_397[k]
                   + f_0 * ki_705[k];

        t_538[k] = -hi_398[k]
                   + f_0 * ki_706[k];

        t_539[k] = -hi_399[k]
                   + f_0 * ki_707[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, hi_400, hi_401, hi_402, hi_403, \
                         hi_404, ki_708, ki_709, ki_710, ki_711, \
                         ki_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -hi_400[k]
                   + f_0 * ki_708[k];

        t_541[k] = -hi_401[k]
                   + f_0 * ki_709[k];

        t_542[k] = -hi_402[k]
                   + f_0 * ki_710[k];

        t_543[k] = -hi_403[k]
                   + f_0 * ki_711[k];

        t_544[k] = -hi_404[k]
                   + f_0 * ki_712[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, hi_405, hi_406, hi_407, hi_408, \
                         hi_409, ki_713, ki_714, ki_715, ki_716, \
                         ki_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -hi_405[k]
                   + f_0 * ki_713[k];

        t_546[k] = -hi_406[k]
                   + f_0 * ki_714[k];

        t_547[k] = -hi_407[k]
                   + f_0 * ki_715[k];

        t_548[k] = -hi_408[k]
                   + f_0 * ki_716[k];

        t_549[k] = -hi_409[k]
                   + f_0 * ki_717[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, hi_410, hi_411, hi_412, hi_413, \
                         hi_414, ki_718, ki_719, ki_720, ki_721, \
                         ki_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -hi_410[k]
                   + f_0 * ki_718[k];

        t_551[k] = -hi_411[k]
                   + f_0 * ki_719[k];

        t_552[k] = -hi_412[k]
                   + f_0 * ki_720[k];

        t_553[k] = -hi_413[k]
                   + f_0 * ki_721[k];

        t_554[k] = -hi_414[k]
                   + f_0 * ki_722[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, hi_415, hi_416, hi_417, hi_418, \
                         hi_419, ki_723, ki_724, ki_725, ki_726, \
                         ki_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -hi_415[k]
                   + f_0 * ki_723[k];

        t_556[k] = -hi_416[k]
                   + f_0 * ki_724[k];

        t_557[k] = -hi_417[k]
                   + f_0 * ki_725[k];

        t_558[k] = -hi_418[k]
                   + f_0 * ki_726[k];

        t_559[k] = -hi_419[k]
                   + f_0 * ki_727[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, t_565, t_566, t_567, ki_728, \
                         ki_729, ki_730, ki_731, ki_732, ki_733, ki_734, \
                         ki_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_0 * ki_728[k];

        t_561[k] = f_0 * ki_729[k];

        t_562[k] = f_0 * ki_730[k];

        t_563[k] = f_0 * ki_731[k];

        t_564[k] = f_0 * ki_732[k];

        t_565[k] = f_0 * ki_733[k];

        t_566[k] = f_0 * ki_734[k];

        t_567[k] = f_0 * ki_735[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, t_572, t_573, t_574, t_575, ki_736, \
                         ki_737, ki_738, ki_739, ki_740, ki_741, ki_742, \
                         ki_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_0 * ki_736[k];

        t_569[k] = f_0 * ki_737[k];

        t_570[k] = f_0 * ki_738[k];

        t_571[k] = f_0 * ki_739[k];

        t_572[k] = f_0 * ki_740[k];

        t_573[k] = f_0 * ki_741[k];

        t_574[k] = f_0 * ki_742[k];

        t_575[k] = f_0 * ki_743[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, t_581, t_582, t_583, ki_744, \
                         ki_745, ki_746, ki_747, ki_748, ki_749, ki_750, \
                         ki_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_0 * ki_744[k];

        t_577[k] = f_0 * ki_745[k];

        t_578[k] = f_0 * ki_746[k];

        t_579[k] = f_0 * ki_747[k];

        t_580[k] = f_0 * ki_748[k];

        t_581[k] = f_0 * ki_749[k];

        t_582[k] = f_0 * ki_750[k];

        t_583[k] = f_0 * ki_751[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, t_588, t_589, hi_420, hi_421, ki_752, \
                         ki_753, ki_754, ki_755, ki_784, ki_785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_0 * ki_752[k];

        t_585[k] = f_0 * ki_753[k];

        t_586[k] = f_0 * ki_754[k];

        t_587[k] = f_0 * ki_755[k];

        t_588[k] = -6.0 * hi_420[k]
                   + f_0 * ki_784[k];

        t_589[k] = -6.0 * hi_421[k]
                   + f_0 * ki_785[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, hi_422, hi_423, hi_424, hi_425, \
                         hi_426, ki_786, ki_787, ki_788, ki_789, \
                         ki_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -6.0 * hi_422[k]
                   + f_0 * ki_786[k];

        t_591[k] = -6.0 * hi_423[k]
                   + f_0 * ki_787[k];

        t_592[k] = -6.0 * hi_424[k]
                   + f_0 * ki_788[k];

        t_593[k] = -6.0 * hi_425[k]
                   + f_0 * ki_789[k];

        t_594[k] = -6.0 * hi_426[k]
                   + f_0 * ki_790[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, hi_427, hi_428, hi_429, hi_430, \
                         hi_431, ki_791, ki_792, ki_793, ki_794, \
                         ki_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -6.0 * hi_427[k]
                   + f_0 * ki_791[k];

        t_596[k] = -6.0 * hi_428[k]
                   + f_0 * ki_792[k];

        t_597[k] = -6.0 * hi_429[k]
                   + f_0 * ki_793[k];

        t_598[k] = -6.0 * hi_430[k]
                   + f_0 * ki_794[k];

        t_599[k] = -6.0 * hi_431[k]
                   + f_0 * ki_795[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, hi_432, hi_433, hi_434, hi_435, \
                         hi_436, ki_796, ki_797, ki_798, ki_799, \
                         ki_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -6.0 * hi_432[k]
                   + f_0 * ki_796[k];

        t_601[k] = -6.0 * hi_433[k]
                   + f_0 * ki_797[k];

        t_602[k] = -6.0 * hi_434[k]
                   + f_0 * ki_798[k];

        t_603[k] = -6.0 * hi_435[k]
                   + f_0 * ki_799[k];

        t_604[k] = -6.0 * hi_436[k]
                   + f_0 * ki_800[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, hi_437, hi_438, hi_439, hi_440, \
                         hi_441, ki_801, ki_802, ki_803, ki_804, \
                         ki_805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -6.0 * hi_437[k]
                   + f_0 * ki_801[k];

        t_606[k] = -6.0 * hi_438[k]
                   + f_0 * ki_802[k];

        t_607[k] = -6.0 * hi_439[k]
                   + f_0 * ki_803[k];

        t_608[k] = -6.0 * hi_440[k]
                   + f_0 * ki_804[k];

        t_609[k] = -6.0 * hi_441[k]
                   + f_0 * ki_805[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, hi_442, hi_443, hi_444, hi_445, \
                         hi_446, ki_806, ki_807, ki_808, ki_809, \
                         ki_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -6.0 * hi_442[k]
                   + f_0 * ki_806[k];

        t_611[k] = -6.0 * hi_443[k]
                   + f_0 * ki_807[k];

        t_612[k] = -6.0 * hi_444[k]
                   + f_0 * ki_808[k];

        t_613[k] = -6.0 * hi_445[k]
                   + f_0 * ki_809[k];

        t_614[k] = -6.0 * hi_446[k]
                   + f_0 * ki_810[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, hi_447, hi_448, hi_449, hi_450, \
                         hi_451, ki_811, ki_812, ki_813, ki_814, \
                         ki_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -6.0 * hi_447[k]
                   + f_0 * ki_811[k];

        t_616[k] = -5.0 * hi_448[k]
                   + f_0 * ki_812[k];

        t_617[k] = -5.0 * hi_449[k]
                   + f_0 * ki_813[k];

        t_618[k] = -5.0 * hi_450[k]
                   + f_0 * ki_814[k];

        t_619[k] = -5.0 * hi_451[k]
                   + f_0 * ki_815[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, hi_452, hi_453, hi_454, hi_455, \
                         hi_456, ki_816, ki_817, ki_818, ki_819, \
                         ki_820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -5.0 * hi_452[k]
                   + f_0 * ki_816[k];

        t_621[k] = -5.0 * hi_453[k]
                   + f_0 * ki_817[k];

        t_622[k] = -5.0 * hi_454[k]
                   + f_0 * ki_818[k];

        t_623[k] = -5.0 * hi_455[k]
                   + f_0 * ki_819[k];

        t_624[k] = -5.0 * hi_456[k]
                   + f_0 * ki_820[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, hi_457, hi_458, hi_459, hi_460, \
                         hi_461, ki_821, ki_822, ki_823, ki_824, \
                         ki_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -5.0 * hi_457[k]
                   + f_0 * ki_821[k];

        t_626[k] = -5.0 * hi_458[k]
                   + f_0 * ki_822[k];

        t_627[k] = -5.0 * hi_459[k]
                   + f_0 * ki_823[k];

        t_628[k] = -5.0 * hi_460[k]
                   + f_0 * ki_824[k];

        t_629[k] = -5.0 * hi_461[k]
                   + f_0 * ki_825[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, hi_462, hi_463, hi_464, hi_465, \
                         hi_466, ki_826, ki_827, ki_828, ki_829, \
                         ki_830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -5.0 * hi_462[k]
                   + f_0 * ki_826[k];

        t_631[k] = -5.0 * hi_463[k]
                   + f_0 * ki_827[k];

        t_632[k] = -5.0 * hi_464[k]
                   + f_0 * ki_828[k];

        t_633[k] = -5.0 * hi_465[k]
                   + f_0 * ki_829[k];

        t_634[k] = -5.0 * hi_466[k]
                   + f_0 * ki_830[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, hi_467, hi_468, hi_469, hi_470, \
                         hi_471, ki_831, ki_832, ki_833, ki_834, \
                         ki_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -5.0 * hi_467[k]
                   + f_0 * ki_831[k];

        t_636[k] = -5.0 * hi_468[k]
                   + f_0 * ki_832[k];

        t_637[k] = -5.0 * hi_469[k]
                   + f_0 * ki_833[k];

        t_638[k] = -5.0 * hi_470[k]
                   + f_0 * ki_834[k];

        t_639[k] = -5.0 * hi_471[k]
                   + f_0 * ki_835[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, hi_472, hi_473, hi_474, hi_475, \
                         hi_476, ki_836, ki_837, ki_838, ki_839, \
                         ki_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -5.0 * hi_472[k]
                   + f_0 * ki_836[k];

        t_641[k] = -5.0 * hi_473[k]
                   + f_0 * ki_837[k];

        t_642[k] = -5.0 * hi_474[k]
                   + f_0 * ki_838[k];

        t_643[k] = -5.0 * hi_475[k]
                   + f_0 * ki_839[k];

        t_644[k] = -4.0 * hi_476[k]
                   + f_0 * ki_840[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, hi_477, hi_478, hi_479, hi_480, \
                         hi_481, ki_841, ki_842, ki_843, ki_844, \
                         ki_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -4.0 * hi_477[k]
                   + f_0 * ki_841[k];

        t_646[k] = -4.0 * hi_478[k]
                   + f_0 * ki_842[k];

        t_647[k] = -4.0 * hi_479[k]
                   + f_0 * ki_843[k];

        t_648[k] = -4.0 * hi_480[k]
                   + f_0 * ki_844[k];

        t_649[k] = -4.0 * hi_481[k]
                   + f_0 * ki_845[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, hi_482, hi_483, hi_484, hi_485, \
                         hi_486, ki_846, ki_847, ki_848, ki_849, \
                         ki_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -4.0 * hi_482[k]
                   + f_0 * ki_846[k];

        t_651[k] = -4.0 * hi_483[k]
                   + f_0 * ki_847[k];

        t_652[k] = -4.0 * hi_484[k]
                   + f_0 * ki_848[k];

        t_653[k] = -4.0 * hi_485[k]
                   + f_0 * ki_849[k];

        t_654[k] = -4.0 * hi_486[k]
                   + f_0 * ki_850[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hi_487 = buffer.data(hi + 487);
    const auto *hi_488 = buffer.data(hi + 488);
    const auto *hi_489 = buffer.data(hi + 489);
    const auto *hi_490 = buffer.data(hi + 490);
    const auto *hi_491 = buffer.data(hi + 491);
    const auto *hi_492 = buffer.data(hi + 492);
    const auto *hi_493 = buffer.data(hi + 493);
    const auto *hi_494 = buffer.data(hi + 494);
    const auto *hi_495 = buffer.data(hi + 495);
    const auto *hi_496 = buffer.data(hi + 496);
    const auto *hi_497 = buffer.data(hi + 497);
    const auto *hi_498 = buffer.data(hi + 498);
    const auto *hi_499 = buffer.data(hi + 499);
    const auto *hi_500 = buffer.data(hi + 500);
    const auto *hi_501 = buffer.data(hi + 501);
    const auto *hi_502 = buffer.data(hi + 502);
    const auto *hi_503 = buffer.data(hi + 503);
    const auto *hi_504 = buffer.data(hi + 504);
    const auto *hi_505 = buffer.data(hi + 505);
    const auto *hi_506 = buffer.data(hi + 506);
    const auto *hi_507 = buffer.data(hi + 507);
    const auto *hi_508 = buffer.data(hi + 508);
    const auto *hi_509 = buffer.data(hi + 509);
    const auto *hi_510 = buffer.data(hi + 510);
    const auto *hi_511 = buffer.data(hi + 511);
    const auto *hi_512 = buffer.data(hi + 512);
    const auto *hi_513 = buffer.data(hi + 513);
    const auto *hi_514 = buffer.data(hi + 514);
    const auto *hi_515 = buffer.data(hi + 515);
    const auto *hi_516 = buffer.data(hi + 516);
    const auto *hi_517 = buffer.data(hi + 517);
    const auto *hi_518 = buffer.data(hi + 518);
    const auto *hi_519 = buffer.data(hi + 519);
    const auto *hi_520 = buffer.data(hi + 520);
    const auto *hi_521 = buffer.data(hi + 521);
    const auto *hi_522 = buffer.data(hi + 522);
    const auto *hi_523 = buffer.data(hi + 523);
    const auto *hi_524 = buffer.data(hi + 524);
    const auto *hi_525 = buffer.data(hi + 525);
    const auto *hi_526 = buffer.data(hi + 526);
    const auto *hi_527 = buffer.data(hi + 527);
    const auto *hi_528 = buffer.data(hi + 528);
    const auto *hi_529 = buffer.data(hi + 529);
    const auto *hi_530 = buffer.data(hi + 530);
    const auto *hi_531 = buffer.data(hi + 531);
    const auto *hi_532 = buffer.data(hi + 532);
    const auto *hi_533 = buffer.data(hi + 533);
    const auto *hi_534 = buffer.data(hi + 534);
    const auto *hi_535 = buffer.data(hi + 535);
    const auto *hi_536 = buffer.data(hi + 536);
    const auto *hi_537 = buffer.data(hi + 537);
    const auto *hi_538 = buffer.data(hi + 538);
    const auto *hi_539 = buffer.data(hi + 539);
    const auto *hi_540 = buffer.data(hi + 540);
    const auto *hi_541 = buffer.data(hi + 541);
    const auto *hi_542 = buffer.data(hi + 542);
    const auto *hi_543 = buffer.data(hi + 543);
    const auto *hi_544 = buffer.data(hi + 544);
    const auto *hi_545 = buffer.data(hi + 545);
    const auto *hi_546 = buffer.data(hi + 546);
    const auto *hi_547 = buffer.data(hi + 547);
    const auto *hi_548 = buffer.data(hi + 548);
    const auto *hi_549 = buffer.data(hi + 549);
    const auto *hi_550 = buffer.data(hi + 550);
    const auto *hi_551 = buffer.data(hi + 551);
    const auto *hi_552 = buffer.data(hi + 552);
    const auto *hi_553 = buffer.data(hi + 553);
    const auto *hi_554 = buffer.data(hi + 554);
    const auto *hi_555 = buffer.data(hi + 555);
    const auto *hi_556 = buffer.data(hi + 556);
    const auto *hi_557 = buffer.data(hi + 557);
    const auto *hi_558 = buffer.data(hi + 558);
    const auto *hi_559 = buffer.data(hi + 559);
    const auto *hi_560 = buffer.data(hi + 560);
    const auto *hi_561 = buffer.data(hi + 561);
    const auto *hi_562 = buffer.data(hi + 562);
    const auto *hi_563 = buffer.data(hi + 563);
    const auto *hi_564 = buffer.data(hi + 564);
    const auto *hi_565 = buffer.data(hi + 565);
    const auto *hi_566 = buffer.data(hi + 566);
    const auto *hi_567 = buffer.data(hi + 567);
    const auto *hi_568 = buffer.data(hi + 568);
    const auto *hi_569 = buffer.data(hi + 569);
    const auto *hi_570 = buffer.data(hi + 570);
    const auto *hi_571 = buffer.data(hi + 571);
    const auto *hi_572 = buffer.data(hi + 572);
    const auto *hi_573 = buffer.data(hi + 573);
    const auto *hi_574 = buffer.data(hi + 574);
    const auto *hi_575 = buffer.data(hi + 575);
    const auto *hi_576 = buffer.data(hi + 576);
    const auto *hi_577 = buffer.data(hi + 577);
    const auto *hi_578 = buffer.data(hi + 578);
    const auto *hi_579 = buffer.data(hi + 579);
    const auto *hi_580 = buffer.data(hi + 580);
    const auto *hi_581 = buffer.data(hi + 581);
    const auto *hi_582 = buffer.data(hi + 582);
    const auto *hi_583 = buffer.data(hi + 583);
    const auto *hi_584 = buffer.data(hi + 584);
    const auto *hi_585 = buffer.data(hi + 585);
    const auto *hi_586 = buffer.data(hi + 586);
    const auto *hi_587 = buffer.data(hi + 587);

    const auto *ki_851 = buffer.data(ki + 851);
    const auto *ki_852 = buffer.data(ki + 852);
    const auto *ki_853 = buffer.data(ki + 853);
    const auto *ki_854 = buffer.data(ki + 854);
    const auto *ki_855 = buffer.data(ki + 855);
    const auto *ki_856 = buffer.data(ki + 856);
    const auto *ki_857 = buffer.data(ki + 857);
    const auto *ki_858 = buffer.data(ki + 858);
    const auto *ki_859 = buffer.data(ki + 859);
    const auto *ki_860 = buffer.data(ki + 860);
    const auto *ki_861 = buffer.data(ki + 861);
    const auto *ki_862 = buffer.data(ki + 862);
    const auto *ki_863 = buffer.data(ki + 863);
    const auto *ki_864 = buffer.data(ki + 864);
    const auto *ki_865 = buffer.data(ki + 865);
    const auto *ki_866 = buffer.data(ki + 866);
    const auto *ki_867 = buffer.data(ki + 867);
    const auto *ki_868 = buffer.data(ki + 868);
    const auto *ki_869 = buffer.data(ki + 869);
    const auto *ki_870 = buffer.data(ki + 870);
    const auto *ki_871 = buffer.data(ki + 871);
    const auto *ki_872 = buffer.data(ki + 872);
    const auto *ki_873 = buffer.data(ki + 873);
    const auto *ki_874 = buffer.data(ki + 874);
    const auto *ki_875 = buffer.data(ki + 875);
    const auto *ki_876 = buffer.data(ki + 876);
    const auto *ki_877 = buffer.data(ki + 877);
    const auto *ki_878 = buffer.data(ki + 878);
    const auto *ki_879 = buffer.data(ki + 879);
    const auto *ki_880 = buffer.data(ki + 880);
    const auto *ki_881 = buffer.data(ki + 881);
    const auto *ki_882 = buffer.data(ki + 882);
    const auto *ki_883 = buffer.data(ki + 883);
    const auto *ki_884 = buffer.data(ki + 884);
    const auto *ki_885 = buffer.data(ki + 885);
    const auto *ki_886 = buffer.data(ki + 886);
    const auto *ki_887 = buffer.data(ki + 887);
    const auto *ki_888 = buffer.data(ki + 888);
    const auto *ki_889 = buffer.data(ki + 889);
    const auto *ki_890 = buffer.data(ki + 890);
    const auto *ki_891 = buffer.data(ki + 891);
    const auto *ki_892 = buffer.data(ki + 892);
    const auto *ki_893 = buffer.data(ki + 893);
    const auto *ki_894 = buffer.data(ki + 894);
    const auto *ki_895 = buffer.data(ki + 895);
    const auto *ki_896 = buffer.data(ki + 896);
    const auto *ki_897 = buffer.data(ki + 897);
    const auto *ki_898 = buffer.data(ki + 898);
    const auto *ki_899 = buffer.data(ki + 899);
    const auto *ki_900 = buffer.data(ki + 900);
    const auto *ki_901 = buffer.data(ki + 901);
    const auto *ki_902 = buffer.data(ki + 902);
    const auto *ki_903 = buffer.data(ki + 903);
    const auto *ki_904 = buffer.data(ki + 904);
    const auto *ki_905 = buffer.data(ki + 905);
    const auto *ki_906 = buffer.data(ki + 906);
    const auto *ki_907 = buffer.data(ki + 907);
    const auto *ki_908 = buffer.data(ki + 908);
    const auto *ki_909 = buffer.data(ki + 909);
    const auto *ki_910 = buffer.data(ki + 910);
    const auto *ki_911 = buffer.data(ki + 911);
    const auto *ki_912 = buffer.data(ki + 912);
    const auto *ki_913 = buffer.data(ki + 913);
    const auto *ki_914 = buffer.data(ki + 914);
    const auto *ki_915 = buffer.data(ki + 915);
    const auto *ki_916 = buffer.data(ki + 916);
    const auto *ki_917 = buffer.data(ki + 917);
    const auto *ki_918 = buffer.data(ki + 918);
    const auto *ki_919 = buffer.data(ki + 919);
    const auto *ki_920 = buffer.data(ki + 920);
    const auto *ki_921 = buffer.data(ki + 921);
    const auto *ki_922 = buffer.data(ki + 922);
    const auto *ki_923 = buffer.data(ki + 923);
    const auto *ki_924 = buffer.data(ki + 924);
    const auto *ki_925 = buffer.data(ki + 925);
    const auto *ki_926 = buffer.data(ki + 926);
    const auto *ki_927 = buffer.data(ki + 927);
    const auto *ki_928 = buffer.data(ki + 928);
    const auto *ki_929 = buffer.data(ki + 929);
    const auto *ki_930 = buffer.data(ki + 930);
    const auto *ki_931 = buffer.data(ki + 931);
    const auto *ki_932 = buffer.data(ki + 932);
    const auto *ki_933 = buffer.data(ki + 933);
    const auto *ki_934 = buffer.data(ki + 934);
    const auto *ki_935 = buffer.data(ki + 935);
    const auto *ki_936 = buffer.data(ki + 936);
    const auto *ki_937 = buffer.data(ki + 937);
    const auto *ki_938 = buffer.data(ki + 938);
    const auto *ki_939 = buffer.data(ki + 939);
    const auto *ki_940 = buffer.data(ki + 940);
    const auto *ki_941 = buffer.data(ki + 941);
    const auto *ki_942 = buffer.data(ki + 942);
    const auto *ki_943 = buffer.data(ki + 943);
    const auto *ki_944 = buffer.data(ki + 944);
    const auto *ki_945 = buffer.data(ki + 945);
    const auto *ki_946 = buffer.data(ki + 946);
    const auto *ki_947 = buffer.data(ki + 947);
    const auto *ki_948 = buffer.data(ki + 948);
    const auto *ki_949 = buffer.data(ki + 949);
    const auto *ki_950 = buffer.data(ki + 950);
    const auto *ki_951 = buffer.data(ki + 951);
    const auto *ki_952 = buffer.data(ki + 952);
    const auto *ki_953 = buffer.data(ki + 953);
    const auto *ki_954 = buffer.data(ki + 954);
    const auto *ki_955 = buffer.data(ki + 955);
    const auto *ki_956 = buffer.data(ki + 956);
    const auto *ki_957 = buffer.data(ki + 957);
    const auto *ki_958 = buffer.data(ki + 958);
    const auto *ki_959 = buffer.data(ki + 959);
    const auto *ki_960 = buffer.data(ki + 960);
    const auto *ki_961 = buffer.data(ki + 961);
    const auto *ki_962 = buffer.data(ki + 962);
    const auto *ki_963 = buffer.data(ki + 963);
    const auto *ki_964 = buffer.data(ki + 964);
    const auto *ki_965 = buffer.data(ki + 965);
    const auto *ki_966 = buffer.data(ki + 966);
    const auto *ki_967 = buffer.data(ki + 967);
    const auto *ki_968 = buffer.data(ki + 968);
    const auto *ki_969 = buffer.data(ki + 969);
    const auto *ki_970 = buffer.data(ki + 970);
    const auto *ki_971 = buffer.data(ki + 971);
    const auto *ki_972 = buffer.data(ki + 972);
    const auto *ki_973 = buffer.data(ki + 973);
    const auto *ki_974 = buffer.data(ki + 974);
    const auto *ki_975 = buffer.data(ki + 975);
    const auto *ki_976 = buffer.data(ki + 976);
    const auto *ki_977 = buffer.data(ki + 977);
    const auto *ki_978 = buffer.data(ki + 978);
    const auto *ki_979 = buffer.data(ki + 979);

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, hi_487, hi_488, hi_489, hi_490, \
                         hi_491, ki_851, ki_852, ki_853, ki_854, \
                         ki_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -4.0 * hi_487[k]
                   + f_0 * ki_851[k];

        t_656[k] = -4.0 * hi_488[k]
                   + f_0 * ki_852[k];

        t_657[k] = -4.0 * hi_489[k]
                   + f_0 * ki_853[k];

        t_658[k] = -4.0 * hi_490[k]
                   + f_0 * ki_854[k];

        t_659[k] = -4.0 * hi_491[k]
                   + f_0 * ki_855[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, hi_492, hi_493, hi_494, hi_495, \
                         hi_496, ki_856, ki_857, ki_858, ki_859, \
                         ki_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -4.0 * hi_492[k]
                   + f_0 * ki_856[k];

        t_661[k] = -4.0 * hi_493[k]
                   + f_0 * ki_857[k];

        t_662[k] = -4.0 * hi_494[k]
                   + f_0 * ki_858[k];

        t_663[k] = -4.0 * hi_495[k]
                   + f_0 * ki_859[k];

        t_664[k] = -4.0 * hi_496[k]
                   + f_0 * ki_860[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, hi_497, hi_498, hi_499, hi_500, \
                         hi_501, ki_861, ki_862, ki_863, ki_864, \
                         ki_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -4.0 * hi_497[k]
                   + f_0 * ki_861[k];

        t_666[k] = -4.0 * hi_498[k]
                   + f_0 * ki_862[k];

        t_667[k] = -4.0 * hi_499[k]
                   + f_0 * ki_863[k];

        t_668[k] = -4.0 * hi_500[k]
                   + f_0 * ki_864[k];

        t_669[k] = -4.0 * hi_501[k]
                   + f_0 * ki_865[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, hi_502, hi_503, hi_504, hi_505, \
                         hi_506, ki_866, ki_867, ki_868, ki_869, \
                         ki_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -4.0 * hi_502[k]
                   + f_0 * ki_866[k];

        t_671[k] = -4.0 * hi_503[k]
                   + f_0 * ki_867[k];

        t_672[k] = -3.0 * hi_504[k]
                   + f_0 * ki_868[k];

        t_673[k] = -3.0 * hi_505[k]
                   + f_0 * ki_869[k];

        t_674[k] = -3.0 * hi_506[k]
                   + f_0 * ki_870[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, hi_507, hi_508, hi_509, hi_510, \
                         hi_511, ki_871, ki_872, ki_873, ki_874, \
                         ki_875 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -3.0 * hi_507[k]
                   + f_0 * ki_871[k];

        t_676[k] = -3.0 * hi_508[k]
                   + f_0 * ki_872[k];

        t_677[k] = -3.0 * hi_509[k]
                   + f_0 * ki_873[k];

        t_678[k] = -3.0 * hi_510[k]
                   + f_0 * ki_874[k];

        t_679[k] = -3.0 * hi_511[k]
                   + f_0 * ki_875[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, hi_512, hi_513, hi_514, hi_515, \
                         hi_516, ki_876, ki_877, ki_878, ki_879, \
                         ki_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -3.0 * hi_512[k]
                   + f_0 * ki_876[k];

        t_681[k] = -3.0 * hi_513[k]
                   + f_0 * ki_877[k];

        t_682[k] = -3.0 * hi_514[k]
                   + f_0 * ki_878[k];

        t_683[k] = -3.0 * hi_515[k]
                   + f_0 * ki_879[k];

        t_684[k] = -3.0 * hi_516[k]
                   + f_0 * ki_880[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, hi_517, hi_518, hi_519, hi_520, \
                         hi_521, ki_881, ki_882, ki_883, ki_884, \
                         ki_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -3.0 * hi_517[k]
                   + f_0 * ki_881[k];

        t_686[k] = -3.0 * hi_518[k]
                   + f_0 * ki_882[k];

        t_687[k] = -3.0 * hi_519[k]
                   + f_0 * ki_883[k];

        t_688[k] = -3.0 * hi_520[k]
                   + f_0 * ki_884[k];

        t_689[k] = -3.0 * hi_521[k]
                   + f_0 * ki_885[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, hi_522, hi_523, hi_524, hi_525, \
                         hi_526, ki_886, ki_887, ki_888, ki_889, \
                         ki_890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -3.0 * hi_522[k]
                   + f_0 * ki_886[k];

        t_691[k] = -3.0 * hi_523[k]
                   + f_0 * ki_887[k];

        t_692[k] = -3.0 * hi_524[k]
                   + f_0 * ki_888[k];

        t_693[k] = -3.0 * hi_525[k]
                   + f_0 * ki_889[k];

        t_694[k] = -3.0 * hi_526[k]
                   + f_0 * ki_890[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, hi_527, hi_528, hi_529, hi_530, \
                         hi_531, ki_891, ki_892, ki_893, ki_894, \
                         ki_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -3.0 * hi_527[k]
                   + f_0 * ki_891[k];

        t_696[k] = -3.0 * hi_528[k]
                   + f_0 * ki_892[k];

        t_697[k] = -3.0 * hi_529[k]
                   + f_0 * ki_893[k];

        t_698[k] = -3.0 * hi_530[k]
                   + f_0 * ki_894[k];

        t_699[k] = -3.0 * hi_531[k]
                   + f_0 * ki_895[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, hi_532, hi_533, hi_534, hi_535, \
                         hi_536, ki_896, ki_897, ki_898, ki_899, \
                         ki_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -2.0 * hi_532[k]
                   + f_0 * ki_896[k];

        t_701[k] = -2.0 * hi_533[k]
                   + f_0 * ki_897[k];

        t_702[k] = -2.0 * hi_534[k]
                   + f_0 * ki_898[k];

        t_703[k] = -2.0 * hi_535[k]
                   + f_0 * ki_899[k];

        t_704[k] = -2.0 * hi_536[k]
                   + f_0 * ki_900[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, hi_537, hi_538, hi_539, hi_540, \
                         hi_541, ki_901, ki_902, ki_903, ki_904, \
                         ki_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -2.0 * hi_537[k]
                   + f_0 * ki_901[k];

        t_706[k] = -2.0 * hi_538[k]
                   + f_0 * ki_902[k];

        t_707[k] = -2.0 * hi_539[k]
                   + f_0 * ki_903[k];

        t_708[k] = -2.0 * hi_540[k]
                   + f_0 * ki_904[k];

        t_709[k] = -2.0 * hi_541[k]
                   + f_0 * ki_905[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, hi_542, hi_543, hi_544, hi_545, \
                         hi_546, ki_906, ki_907, ki_908, ki_909, \
                         ki_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -2.0 * hi_542[k]
                   + f_0 * ki_906[k];

        t_711[k] = -2.0 * hi_543[k]
                   + f_0 * ki_907[k];

        t_712[k] = -2.0 * hi_544[k]
                   + f_0 * ki_908[k];

        t_713[k] = -2.0 * hi_545[k]
                   + f_0 * ki_909[k];

        t_714[k] = -2.0 * hi_546[k]
                   + f_0 * ki_910[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, hi_547, hi_548, hi_549, hi_550, \
                         hi_551, ki_911, ki_912, ki_913, ki_914, \
                         ki_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -2.0 * hi_547[k]
                   + f_0 * ki_911[k];

        t_716[k] = -2.0 * hi_548[k]
                   + f_0 * ki_912[k];

        t_717[k] = -2.0 * hi_549[k]
                   + f_0 * ki_913[k];

        t_718[k] = -2.0 * hi_550[k]
                   + f_0 * ki_914[k];

        t_719[k] = -2.0 * hi_551[k]
                   + f_0 * ki_915[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, hi_552, hi_553, hi_554, hi_555, \
                         hi_556, ki_916, ki_917, ki_918, ki_919, \
                         ki_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -2.0 * hi_552[k]
                   + f_0 * ki_916[k];

        t_721[k] = -2.0 * hi_553[k]
                   + f_0 * ki_917[k];

        t_722[k] = -2.0 * hi_554[k]
                   + f_0 * ki_918[k];

        t_723[k] = -2.0 * hi_555[k]
                   + f_0 * ki_919[k];

        t_724[k] = -2.0 * hi_556[k]
                   + f_0 * ki_920[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, hi_557, hi_558, hi_559, hi_560, \
                         hi_561, ki_921, ki_922, ki_923, ki_924, \
                         ki_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -2.0 * hi_557[k]
                   + f_0 * ki_921[k];

        t_726[k] = -2.0 * hi_558[k]
                   + f_0 * ki_922[k];

        t_727[k] = -2.0 * hi_559[k]
                   + f_0 * ki_923[k];

        t_728[k] = -hi_560[k]
                   + f_0 * ki_924[k];

        t_729[k] = -hi_561[k]
                   + f_0 * ki_925[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, hi_562, hi_563, hi_564, hi_565, \
                         hi_566, ki_926, ki_927, ki_928, ki_929, \
                         ki_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -hi_562[k]
                   + f_0 * ki_926[k];

        t_731[k] = -hi_563[k]
                   + f_0 * ki_927[k];

        t_732[k] = -hi_564[k]
                   + f_0 * ki_928[k];

        t_733[k] = -hi_565[k]
                   + f_0 * ki_929[k];

        t_734[k] = -hi_566[k]
                   + f_0 * ki_930[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, hi_567, hi_568, hi_569, hi_570, \
                         hi_571, ki_931, ki_932, ki_933, ki_934, \
                         ki_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -hi_567[k]
                   + f_0 * ki_931[k];

        t_736[k] = -hi_568[k]
                   + f_0 * ki_932[k];

        t_737[k] = -hi_569[k]
                   + f_0 * ki_933[k];

        t_738[k] = -hi_570[k]
                   + f_0 * ki_934[k];

        t_739[k] = -hi_571[k]
                   + f_0 * ki_935[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, hi_572, hi_573, hi_574, hi_575, \
                         hi_576, ki_936, ki_937, ki_938, ki_939, \
                         ki_940 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -hi_572[k]
                   + f_0 * ki_936[k];

        t_741[k] = -hi_573[k]
                   + f_0 * ki_937[k];

        t_742[k] = -hi_574[k]
                   + f_0 * ki_938[k];

        t_743[k] = -hi_575[k]
                   + f_0 * ki_939[k];

        t_744[k] = -hi_576[k]
                   + f_0 * ki_940[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, hi_577, hi_578, hi_579, hi_580, \
                         hi_581, ki_941, ki_942, ki_943, ki_944, \
                         ki_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -hi_577[k]
                   + f_0 * ki_941[k];

        t_746[k] = -hi_578[k]
                   + f_0 * ki_942[k];

        t_747[k] = -hi_579[k]
                   + f_0 * ki_943[k];

        t_748[k] = -hi_580[k]
                   + f_0 * ki_944[k];

        t_749[k] = -hi_581[k]
                   + f_0 * ki_945[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, hi_582, hi_583, hi_584, hi_585, \
                         hi_586, ki_946, ki_947, ki_948, ki_949, \
                         ki_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -hi_582[k]
                   + f_0 * ki_946[k];

        t_751[k] = -hi_583[k]
                   + f_0 * ki_947[k];

        t_752[k] = -hi_584[k]
                   + f_0 * ki_948[k];

        t_753[k] = -hi_585[k]
                   + f_0 * ki_949[k];

        t_754[k] = -hi_586[k]
                   + f_0 * ki_950[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, t_760, t_761, hi_587, ki_951, \
                         ki_952, ki_953, ki_954, ki_955, ki_956, \
                         ki_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -hi_587[k]
                   + f_0 * ki_951[k];

        t_756[k] = f_0 * ki_952[k];

        t_757[k] = f_0 * ki_953[k];

        t_758[k] = f_0 * ki_954[k];

        t_759[k] = f_0 * ki_955[k];

        t_760[k] = f_0 * ki_956[k];

        t_761[k] = f_0 * ki_957[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, t_766, t_767, t_768, t_769, ki_958, \
                         ki_959, ki_960, ki_961, ki_962, ki_963, ki_964, \
                         ki_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_0 * ki_958[k];

        t_763[k] = f_0 * ki_959[k];

        t_764[k] = f_0 * ki_960[k];

        t_765[k] = f_0 * ki_961[k];

        t_766[k] = f_0 * ki_962[k];

        t_767[k] = f_0 * ki_963[k];

        t_768[k] = f_0 * ki_964[k];

        t_769[k] = f_0 * ki_965[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, t_775, t_776, t_777, ki_966, \
                         ki_967, ki_968, ki_969, ki_970, ki_971, ki_972, \
                         ki_973 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_0 * ki_966[k];

        t_771[k] = f_0 * ki_967[k];

        t_772[k] = f_0 * ki_968[k];

        t_773[k] = f_0 * ki_969[k];

        t_774[k] = f_0 * ki_970[k];

        t_775[k] = f_0 * ki_971[k];

        t_776[k] = f_0 * ki_972[k];

        t_777[k] = f_0 * ki_973[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, t_783, ki_974, ki_975, ki_976, \
                         ki_977, ki_978, ki_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_0 * ki_974[k];

        t_779[k] = f_0 * ki_975[k];

        t_780[k] = f_0 * ki_976[k];

        t_781[k] = f_0 * ki_977[k];

        t_782[k] = f_0 * ki_978[k];

        t_783[k] = f_0 * ki_979[k];
    }
}

auto
compute_prim_geom_10_ii_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t hi, const size_t ki,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ii_electron_repulsion_1_piece0(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_1_piece1(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_1_piece2(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_1_piece3(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_1_piece4(buffer, target, hi, ki, ncols, alpha);
}

static auto
compute_prim_geom_10_ii_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
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

    const auto *hi_0 = buffer.data(hi + 0);
    const auto *hi_1 = buffer.data(hi + 1);
    const auto *hi_2 = buffer.data(hi + 2);
    const auto *hi_3 = buffer.data(hi + 3);
    const auto *hi_4 = buffer.data(hi + 4);
    const auto *hi_5 = buffer.data(hi + 5);
    const auto *hi_6 = buffer.data(hi + 6);
    const auto *hi_7 = buffer.data(hi + 7);
    const auto *hi_8 = buffer.data(hi + 8);
    const auto *hi_9 = buffer.data(hi + 9);
    const auto *hi_10 = buffer.data(hi + 10);
    const auto *hi_11 = buffer.data(hi + 11);
    const auto *hi_12 = buffer.data(hi + 12);
    const auto *hi_13 = buffer.data(hi + 13);
    const auto *hi_14 = buffer.data(hi + 14);
    const auto *hi_15 = buffer.data(hi + 15);
    const auto *hi_16 = buffer.data(hi + 16);
    const auto *hi_17 = buffer.data(hi + 17);
    const auto *hi_18 = buffer.data(hi + 18);
    const auto *hi_19 = buffer.data(hi + 19);
    const auto *hi_20 = buffer.data(hi + 20);
    const auto *hi_21 = buffer.data(hi + 21);
    const auto *hi_22 = buffer.data(hi + 22);
    const auto *hi_23 = buffer.data(hi + 23);
    const auto *hi_24 = buffer.data(hi + 24);
    const auto *hi_25 = buffer.data(hi + 25);
    const auto *hi_26 = buffer.data(hi + 26);
    const auto *hi_27 = buffer.data(hi + 27);
    const auto *hi_28 = buffer.data(hi + 28);
    const auto *hi_29 = buffer.data(hi + 29);
    const auto *hi_30 = buffer.data(hi + 30);
    const auto *hi_31 = buffer.data(hi + 31);
    const auto *hi_32 = buffer.data(hi + 32);
    const auto *hi_33 = buffer.data(hi + 33);
    const auto *hi_34 = buffer.data(hi + 34);
    const auto *hi_35 = buffer.data(hi + 35);
    const auto *hi_36 = buffer.data(hi + 36);
    const auto *hi_37 = buffer.data(hi + 37);
    const auto *hi_38 = buffer.data(hi + 38);
    const auto *hi_39 = buffer.data(hi + 39);
    const auto *hi_40 = buffer.data(hi + 40);
    const auto *hi_41 = buffer.data(hi + 41);
    const auto *hi_42 = buffer.data(hi + 42);
    const auto *hi_43 = buffer.data(hi + 43);
    const auto *hi_44 = buffer.data(hi + 44);
    const auto *hi_45 = buffer.data(hi + 45);
    const auto *hi_46 = buffer.data(hi + 46);
    const auto *hi_47 = buffer.data(hi + 47);
    const auto *hi_48 = buffer.data(hi + 48);
    const auto *hi_49 = buffer.data(hi + 49);
    const auto *hi_50 = buffer.data(hi + 50);
    const auto *hi_51 = buffer.data(hi + 51);
    const auto *hi_52 = buffer.data(hi + 52);
    const auto *hi_53 = buffer.data(hi + 53);
    const auto *hi_54 = buffer.data(hi + 54);
    const auto *hi_55 = buffer.data(hi + 55);
    const auto *hi_56 = buffer.data(hi + 56);
    const auto *hi_57 = buffer.data(hi + 57);
    const auto *hi_58 = buffer.data(hi + 58);
    const auto *hi_59 = buffer.data(hi + 59);
    const auto *hi_60 = buffer.data(hi + 60);
    const auto *hi_61 = buffer.data(hi + 61);
    const auto *hi_62 = buffer.data(hi + 62);
    const auto *hi_63 = buffer.data(hi + 63);
    const auto *hi_64 = buffer.data(hi + 64);
    const auto *hi_65 = buffer.data(hi + 65);
    const auto *hi_66 = buffer.data(hi + 66);
    const auto *hi_67 = buffer.data(hi + 67);
    const auto *hi_68 = buffer.data(hi + 68);
    const auto *hi_69 = buffer.data(hi + 69);
    const auto *hi_70 = buffer.data(hi + 70);
    const auto *hi_71 = buffer.data(hi + 71);
    const auto *hi_72 = buffer.data(hi + 72);
    const auto *hi_73 = buffer.data(hi + 73);
    const auto *hi_74 = buffer.data(hi + 74);
    const auto *hi_75 = buffer.data(hi + 75);
    const auto *hi_76 = buffer.data(hi + 76);
    const auto *hi_77 = buffer.data(hi + 77);
    const auto *hi_78 = buffer.data(hi + 78);
    const auto *hi_79 = buffer.data(hi + 79);
    const auto *hi_80 = buffer.data(hi + 80);
    const auto *hi_81 = buffer.data(hi + 81);
    const auto *hi_82 = buffer.data(hi + 82);
    const auto *hi_83 = buffer.data(hi + 83);

    const auto *ki_56 = buffer.data(ki + 56);
    const auto *ki_57 = buffer.data(ki + 57);
    const auto *ki_58 = buffer.data(ki + 58);
    const auto *ki_59 = buffer.data(ki + 59);
    const auto *ki_60 = buffer.data(ki + 60);
    const auto *ki_61 = buffer.data(ki + 61);
    const auto *ki_62 = buffer.data(ki + 62);
    const auto *ki_63 = buffer.data(ki + 63);
    const auto *ki_64 = buffer.data(ki + 64);
    const auto *ki_65 = buffer.data(ki + 65);
    const auto *ki_66 = buffer.data(ki + 66);
    const auto *ki_67 = buffer.data(ki + 67);
    const auto *ki_68 = buffer.data(ki + 68);
    const auto *ki_69 = buffer.data(ki + 69);
    const auto *ki_70 = buffer.data(ki + 70);
    const auto *ki_71 = buffer.data(ki + 71);
    const auto *ki_72 = buffer.data(ki + 72);
    const auto *ki_73 = buffer.data(ki + 73);
    const auto *ki_74 = buffer.data(ki + 74);
    const auto *ki_75 = buffer.data(ki + 75);
    const auto *ki_76 = buffer.data(ki + 76);
    const auto *ki_77 = buffer.data(ki + 77);
    const auto *ki_78 = buffer.data(ki + 78);
    const auto *ki_79 = buffer.data(ki + 79);
    const auto *ki_80 = buffer.data(ki + 80);
    const auto *ki_81 = buffer.data(ki + 81);
    const auto *ki_82 = buffer.data(ki + 82);
    const auto *ki_83 = buffer.data(ki + 83);
    const auto *ki_112 = buffer.data(ki + 112);
    const auto *ki_113 = buffer.data(ki + 113);
    const auto *ki_114 = buffer.data(ki + 114);
    const auto *ki_115 = buffer.data(ki + 115);
    const auto *ki_116 = buffer.data(ki + 116);
    const auto *ki_117 = buffer.data(ki + 117);
    const auto *ki_118 = buffer.data(ki + 118);
    const auto *ki_119 = buffer.data(ki + 119);
    const auto *ki_120 = buffer.data(ki + 120);
    const auto *ki_121 = buffer.data(ki + 121);
    const auto *ki_122 = buffer.data(ki + 122);
    const auto *ki_123 = buffer.data(ki + 123);
    const auto *ki_124 = buffer.data(ki + 124);
    const auto *ki_125 = buffer.data(ki + 125);
    const auto *ki_126 = buffer.data(ki + 126);
    const auto *ki_127 = buffer.data(ki + 127);
    const auto *ki_128 = buffer.data(ki + 128);
    const auto *ki_129 = buffer.data(ki + 129);
    const auto *ki_130 = buffer.data(ki + 130);
    const auto *ki_131 = buffer.data(ki + 131);
    const auto *ki_132 = buffer.data(ki + 132);
    const auto *ki_133 = buffer.data(ki + 133);
    const auto *ki_134 = buffer.data(ki + 134);
    const auto *ki_135 = buffer.data(ki + 135);
    const auto *ki_136 = buffer.data(ki + 136);
    const auto *ki_137 = buffer.data(ki + 137);
    const auto *ki_138 = buffer.data(ki + 138);
    const auto *ki_139 = buffer.data(ki + 139);
    const auto *ki_140 = buffer.data(ki + 140);
    const auto *ki_141 = buffer.data(ki + 141);
    const auto *ki_142 = buffer.data(ki + 142);
    const auto *ki_143 = buffer.data(ki + 143);
    const auto *ki_144 = buffer.data(ki + 144);
    const auto *ki_145 = buffer.data(ki + 145);
    const auto *ki_146 = buffer.data(ki + 146);
    const auto *ki_147 = buffer.data(ki + 147);
    const auto *ki_148 = buffer.data(ki + 148);
    const auto *ki_149 = buffer.data(ki + 149);
    const auto *ki_150 = buffer.data(ki + 150);
    const auto *ki_151 = buffer.data(ki + 151);
    const auto *ki_152 = buffer.data(ki + 152);
    const auto *ki_153 = buffer.data(ki + 153);
    const auto *ki_154 = buffer.data(ki + 154);
    const auto *ki_155 = buffer.data(ki + 155);
    const auto *ki_156 = buffer.data(ki + 156);
    const auto *ki_157 = buffer.data(ki + 157);
    const auto *ki_158 = buffer.data(ki + 158);
    const auto *ki_159 = buffer.data(ki + 159);
    const auto *ki_160 = buffer.data(ki + 160);
    const auto *ki_161 = buffer.data(ki + 161);
    const auto *ki_162 = buffer.data(ki + 162);
    const auto *ki_163 = buffer.data(ki + 163);
    const auto *ki_164 = buffer.data(ki + 164);
    const auto *ki_165 = buffer.data(ki + 165);
    const auto *ki_166 = buffer.data(ki + 166);
    const auto *ki_167 = buffer.data(ki + 167);
    const auto *ki_196 = buffer.data(ki + 196);
    const auto *ki_197 = buffer.data(ki + 197);
    const auto *ki_198 = buffer.data(ki + 198);
    const auto *ki_199 = buffer.data(ki + 199);
    const auto *ki_200 = buffer.data(ki + 200);
    const auto *ki_201 = buffer.data(ki + 201);
    const auto *ki_202 = buffer.data(ki + 202);
    const auto *ki_203 = buffer.data(ki + 203);
    const auto *ki_204 = buffer.data(ki + 204);
    const auto *ki_205 = buffer.data(ki + 205);
    const auto *ki_206 = buffer.data(ki + 206);
    const auto *ki_207 = buffer.data(ki + 207);
    const auto *ki_208 = buffer.data(ki + 208);
    const auto *ki_209 = buffer.data(ki + 209);
    const auto *ki_210 = buffer.data(ki + 210);
    const auto *ki_211 = buffer.data(ki + 211);
    const auto *ki_212 = buffer.data(ki + 212);
    const auto *ki_213 = buffer.data(ki + 213);
    const auto *ki_214 = buffer.data(ki + 214);
    const auto *ki_215 = buffer.data(ki + 215);
    const auto *ki_216 = buffer.data(ki + 216);
    const auto *ki_217 = buffer.data(ki + 217);
    const auto *ki_218 = buffer.data(ki + 218);
    const auto *ki_219 = buffer.data(ki + 219);
    const auto *ki_220 = buffer.data(ki + 220);
    const auto *ki_221 = buffer.data(ki + 221);
    const auto *ki_222 = buffer.data(ki + 222);
    const auto *ki_223 = buffer.data(ki + 223);
    const auto *ki_224 = buffer.data(ki + 224);
    const auto *ki_225 = buffer.data(ki + 225);
    const auto *ki_226 = buffer.data(ki + 226);
    const auto *ki_227 = buffer.data(ki + 227);
    const auto *ki_228 = buffer.data(ki + 228);
    const auto *ki_229 = buffer.data(ki + 229);
    const auto *ki_230 = buffer.data(ki + 230);
    const auto *ki_231 = buffer.data(ki + 231);
    const auto *ki_232 = buffer.data(ki + 232);
    const auto *ki_233 = buffer.data(ki + 233);
    const auto *ki_234 = buffer.data(ki + 234);
    const auto *ki_235 = buffer.data(ki + 235);
    const auto *ki_236 = buffer.data(ki + 236);
    const auto *ki_237 = buffer.data(ki + 237);
    const auto *ki_238 = buffer.data(ki + 238);
    const auto *ki_239 = buffer.data(ki + 239);
    const auto *ki_240 = buffer.data(ki + 240);
    const auto *ki_241 = buffer.data(ki + 241);
    const auto *ki_242 = buffer.data(ki + 242);
    const auto *ki_243 = buffer.data(ki + 243);
    const auto *ki_244 = buffer.data(ki + 244);
    const auto *ki_245 = buffer.data(ki + 245);
    const auto *ki_246 = buffer.data(ki + 246);
    const auto *ki_247 = buffer.data(ki + 247);
    const auto *ki_248 = buffer.data(ki + 248);
    const auto *ki_249 = buffer.data(ki + 249);
    const auto *ki_250 = buffer.data(ki + 250);
    const auto *ki_251 = buffer.data(ki + 251);
    const auto *ki_252 = buffer.data(ki + 252);
    const auto *ki_253 = buffer.data(ki + 253);
    const auto *ki_254 = buffer.data(ki + 254);
    const auto *ki_255 = buffer.data(ki + 255);
    const auto *ki_256 = buffer.data(ki + 256);
    const auto *ki_257 = buffer.data(ki + 257);
    const auto *ki_258 = buffer.data(ki + 258);
    const auto *ki_259 = buffer.data(ki + 259);
    const auto *ki_260 = buffer.data(ki + 260);
    const auto *ki_261 = buffer.data(ki + 261);
    const auto *ki_262 = buffer.data(ki + 262);
    const auto *ki_263 = buffer.data(ki + 263);
    const auto *ki_264 = buffer.data(ki + 264);
    const auto *ki_265 = buffer.data(ki + 265);
    const auto *ki_266 = buffer.data(ki + 266);
    const auto *ki_267 = buffer.data(ki + 267);
    const auto *ki_268 = buffer.data(ki + 268);
    const auto *ki_269 = buffer.data(ki + 269);
    const auto *ki_270 = buffer.data(ki + 270);
    const auto *ki_271 = buffer.data(ki + 271);
    const auto *ki_272 = buffer.data(ki + 272);
    const auto *ki_273 = buffer.data(ki + 273);
    const auto *ki_274 = buffer.data(ki + 274);
    const auto *ki_275 = buffer.data(ki + 275);
    const auto *ki_276 = buffer.data(ki + 276);
    const auto *ki_277 = buffer.data(ki + 277);
    const auto *ki_278 = buffer.data(ki + 278);
    const auto *ki_279 = buffer.data(ki + 279);
    const auto *ki_308 = buffer.data(ki + 308);
    const auto *ki_309 = buffer.data(ki + 309);
    const auto *ki_310 = buffer.data(ki + 310);
    const auto *ki_311 = buffer.data(ki + 311);
    const auto *ki_312 = buffer.data(ki + 312);
    const auto *ki_313 = buffer.data(ki + 313);
    const auto *ki_314 = buffer.data(ki + 314);
    const auto *ki_315 = buffer.data(ki + 315);
    const auto *ki_316 = buffer.data(ki + 316);
    const auto *ki_317 = buffer.data(ki + 317);
    const auto *ki_318 = buffer.data(ki + 318);
    const auto *ki_319 = buffer.data(ki + 319);
    const auto *ki_320 = buffer.data(ki + 320);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ki_56, ki_57, ki_58, ki_59, \
                         ki_60, ki_61, ki_62, ki_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ki_56[k];

        t_1[k] = f_0 * ki_57[k];

        t_2[k] = f_0 * ki_58[k];

        t_3[k] = f_0 * ki_59[k];

        t_4[k] = f_0 * ki_60[k];

        t_5[k] = f_0 * ki_61[k];

        t_6[k] = f_0 * ki_62[k];

        t_7[k] = f_0 * ki_63[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, ki_64, ki_65, ki_66, \
                         ki_67, ki_68, ki_69, ki_70, ki_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ki_64[k];

        t_9[k] = f_0 * ki_65[k];

        t_10[k] = f_0 * ki_66[k];

        t_11[k] = f_0 * ki_67[k];

        t_12[k] = f_0 * ki_68[k];

        t_13[k] = f_0 * ki_69[k];

        t_14[k] = f_0 * ki_70[k];

        t_15[k] = f_0 * ki_71[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, ki_72, ki_73, ki_74, \
                         ki_75, ki_76, ki_77, ki_78, ki_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * ki_72[k];

        t_17[k] = f_0 * ki_73[k];

        t_18[k] = f_0 * ki_74[k];

        t_19[k] = f_0 * ki_75[k];

        t_20[k] = f_0 * ki_76[k];

        t_21[k] = f_0 * ki_77[k];

        t_22[k] = f_0 * ki_78[k];

        t_23[k] = f_0 * ki_79[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, ki_80, ki_81, ki_82, \
                         ki_83, ki_112, ki_113, ki_114, ki_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * ki_80[k];

        t_25[k] = f_0 * ki_81[k];

        t_26[k] = f_0 * ki_82[k];

        t_27[k] = f_0 * ki_83[k];

        t_28[k] = f_0 * ki_112[k];

        t_29[k] = f_0 * ki_113[k];

        t_30[k] = f_0 * ki_114[k];

        t_31[k] = f_0 * ki_115[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, ki_116, ki_117, \
                         ki_118, ki_119, ki_120, ki_121, ki_122, \
                         ki_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ki_116[k];

        t_33[k] = f_0 * ki_117[k];

        t_34[k] = f_0 * ki_118[k];

        t_35[k] = f_0 * ki_119[k];

        t_36[k] = f_0 * ki_120[k];

        t_37[k] = f_0 * ki_121[k];

        t_38[k] = f_0 * ki_122[k];

        t_39[k] = f_0 * ki_123[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, ki_124, ki_125, \
                         ki_126, ki_127, ki_128, ki_129, ki_130, \
                         ki_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * ki_124[k];

        t_41[k] = f_0 * ki_125[k];

        t_42[k] = f_0 * ki_126[k];

        t_43[k] = f_0 * ki_127[k];

        t_44[k] = f_0 * ki_128[k];

        t_45[k] = f_0 * ki_129[k];

        t_46[k] = f_0 * ki_130[k];

        t_47[k] = f_0 * ki_131[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, ki_132, ki_133, \
                         ki_134, ki_135, ki_136, ki_137, ki_138, \
                         ki_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * ki_132[k];

        t_49[k] = f_0 * ki_133[k];

        t_50[k] = f_0 * ki_134[k];

        t_51[k] = f_0 * ki_135[k];

        t_52[k] = f_0 * ki_136[k];

        t_53[k] = f_0 * ki_137[k];

        t_54[k] = f_0 * ki_138[k];

        t_55[k] = f_0 * ki_139[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, hi_0, hi_1, hi_2, hi_3, hi_4, ki_140, \
                         ki_141, ki_142, ki_143, ki_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -hi_0[k]
                  + f_0 * ki_140[k];

        t_57[k] = -hi_1[k]
                  + f_0 * ki_141[k];

        t_58[k] = -hi_2[k]
                  + f_0 * ki_142[k];

        t_59[k] = -hi_3[k]
                  + f_0 * ki_143[k];

        t_60[k] = -hi_4[k]
                  + f_0 * ki_144[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, hi_5, hi_6, hi_7, hi_8, hi_9, ki_145, \
                         ki_146, ki_147, ki_148, ki_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -hi_5[k]
                  + f_0 * ki_145[k];

        t_62[k] = -hi_6[k]
                  + f_0 * ki_146[k];

        t_63[k] = -hi_7[k]
                  + f_0 * ki_147[k];

        t_64[k] = -hi_8[k]
                  + f_0 * ki_148[k];

        t_65[k] = -hi_9[k]
                  + f_0 * ki_149[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, hi_10, hi_11, hi_12, hi_13, hi_14, \
                         ki_150, ki_151, ki_152, ki_153, ki_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -hi_10[k]
                  + f_0 * ki_150[k];

        t_67[k] = -hi_11[k]
                  + f_0 * ki_151[k];

        t_68[k] = -hi_12[k]
                  + f_0 * ki_152[k];

        t_69[k] = -hi_13[k]
                  + f_0 * ki_153[k];

        t_70[k] = -hi_14[k]
                  + f_0 * ki_154[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, hi_15, hi_16, hi_17, hi_18, hi_19, \
                         ki_155, ki_156, ki_157, ki_158, ki_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -hi_15[k]
                  + f_0 * ki_155[k];

        t_72[k] = -hi_16[k]
                  + f_0 * ki_156[k];

        t_73[k] = -hi_17[k]
                  + f_0 * ki_157[k];

        t_74[k] = -hi_18[k]
                  + f_0 * ki_158[k];

        t_75[k] = -hi_19[k]
                  + f_0 * ki_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, hi_20, hi_21, hi_22, hi_23, hi_24, \
                         ki_160, ki_161, ki_162, ki_163, ki_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -hi_20[k]
                  + f_0 * ki_160[k];

        t_77[k] = -hi_21[k]
                  + f_0 * ki_161[k];

        t_78[k] = -hi_22[k]
                  + f_0 * ki_162[k];

        t_79[k] = -hi_23[k]
                  + f_0 * ki_163[k];

        t_80[k] = -hi_24[k]
                  + f_0 * ki_164[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, hi_25, hi_26, hi_27, ki_165, \
                         ki_166, ki_167, ki_196, ki_197, ki_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -hi_25[k]
                  + f_0 * ki_165[k];

        t_82[k] = -hi_26[k]
                  + f_0 * ki_166[k];

        t_83[k] = -hi_27[k]
                  + f_0 * ki_167[k];

        t_84[k] = f_0 * ki_196[k];

        t_85[k] = f_0 * ki_197[k];

        t_86[k] = f_0 * ki_198[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, t_93, t_94, ki_199, ki_200, \
                         ki_201, ki_202, ki_203, ki_204, ki_205, \
                         ki_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_0 * ki_199[k];

        t_88[k] = f_0 * ki_200[k];

        t_89[k] = f_0 * ki_201[k];

        t_90[k] = f_0 * ki_202[k];

        t_91[k] = f_0 * ki_203[k];

        t_92[k] = f_0 * ki_204[k];

        t_93[k] = f_0 * ki_205[k];

        t_94[k] = f_0 * ki_206[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, t_100, t_101, t_102, ki_207, ki_208, \
                         ki_209, ki_210, ki_211, ki_212, ki_213, \
                         ki_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_0 * ki_207[k];

        t_96[k] = f_0 * ki_208[k];

        t_97[k] = f_0 * ki_209[k];

        t_98[k] = f_0 * ki_210[k];

        t_99[k] = f_0 * ki_211[k];

        t_100[k] = f_0 * ki_212[k];

        t_101[k] = f_0 * ki_213[k];

        t_102[k] = f_0 * ki_214[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, t_109, t_110, ki_215, \
                         ki_216, ki_217, ki_218, ki_219, ki_220, ki_221, \
                         ki_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_0 * ki_215[k];

        t_104[k] = f_0 * ki_216[k];

        t_105[k] = f_0 * ki_217[k];

        t_106[k] = f_0 * ki_218[k];

        t_107[k] = f_0 * ki_219[k];

        t_108[k] = f_0 * ki_220[k];

        t_109[k] = f_0 * ki_221[k];

        t_110[k] = f_0 * ki_222[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, hi_28, hi_29, hi_30, hi_31, \
                         ki_223, ki_224, ki_225, ki_226, ki_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_0 * ki_223[k];

        t_112[k] = -hi_28[k]
                   + f_0 * ki_224[k];

        t_113[k] = -hi_29[k]
                   + f_0 * ki_225[k];

        t_114[k] = -hi_30[k]
                   + f_0 * ki_226[k];

        t_115[k] = -hi_31[k]
                   + f_0 * ki_227[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, hi_32, hi_33, hi_34, hi_35, hi_36, \
                         ki_228, ki_229, ki_230, ki_231, ki_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = -hi_32[k]
                   + f_0 * ki_228[k];

        t_117[k] = -hi_33[k]
                   + f_0 * ki_229[k];

        t_118[k] = -hi_34[k]
                   + f_0 * ki_230[k];

        t_119[k] = -hi_35[k]
                   + f_0 * ki_231[k];

        t_120[k] = -hi_36[k]
                   + f_0 * ki_232[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, hi_37, hi_38, hi_39, hi_40, hi_41, \
                         ki_233, ki_234, ki_235, ki_236, ki_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -hi_37[k]
                   + f_0 * ki_233[k];

        t_122[k] = -hi_38[k]
                   + f_0 * ki_234[k];

        t_123[k] = -hi_39[k]
                   + f_0 * ki_235[k];

        t_124[k] = -hi_40[k]
                   + f_0 * ki_236[k];

        t_125[k] = -hi_41[k]
                   + f_0 * ki_237[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, hi_42, hi_43, hi_44, hi_45, hi_46, \
                         ki_238, ki_239, ki_240, ki_241, ki_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = -hi_42[k]
                   + f_0 * ki_238[k];

        t_127[k] = -hi_43[k]
                   + f_0 * ki_239[k];

        t_128[k] = -hi_44[k]
                   + f_0 * ki_240[k];

        t_129[k] = -hi_45[k]
                   + f_0 * ki_241[k];

        t_130[k] = -hi_46[k]
                   + f_0 * ki_242[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, hi_47, hi_48, hi_49, hi_50, hi_51, \
                         ki_243, ki_244, ki_245, ki_246, ki_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = -hi_47[k]
                   + f_0 * ki_243[k];

        t_132[k] = -hi_48[k]
                   + f_0 * ki_244[k];

        t_133[k] = -hi_49[k]
                   + f_0 * ki_245[k];

        t_134[k] = -hi_50[k]
                   + f_0 * ki_246[k];

        t_135[k] = -hi_51[k]
                   + f_0 * ki_247[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, hi_52, hi_53, hi_54, hi_55, hi_56, \
                         ki_248, ki_249, ki_250, ki_251, ki_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = -hi_52[k]
                   + f_0 * ki_248[k];

        t_137[k] = -hi_53[k]
                   + f_0 * ki_249[k];

        t_138[k] = -hi_54[k]
                   + f_0 * ki_250[k];

        t_139[k] = -hi_55[k]
                   + f_0 * ki_251[k];

        t_140[k] = -2.0 * hi_56[k]
                   + f_0 * ki_252[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, hi_57, hi_58, hi_59, hi_60, hi_61, \
                         ki_253, ki_254, ki_255, ki_256, ki_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = -2.0 * hi_57[k]
                   + f_0 * ki_253[k];

        t_142[k] = -2.0 * hi_58[k]
                   + f_0 * ki_254[k];

        t_143[k] = -2.0 * hi_59[k]
                   + f_0 * ki_255[k];

        t_144[k] = -2.0 * hi_60[k]
                   + f_0 * ki_256[k];

        t_145[k] = -2.0 * hi_61[k]
                   + f_0 * ki_257[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, hi_62, hi_63, hi_64, hi_65, hi_66, \
                         ki_258, ki_259, ki_260, ki_261, ki_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = -2.0 * hi_62[k]
                   + f_0 * ki_258[k];

        t_147[k] = -2.0 * hi_63[k]
                   + f_0 * ki_259[k];

        t_148[k] = -2.0 * hi_64[k]
                   + f_0 * ki_260[k];

        t_149[k] = -2.0 * hi_65[k]
                   + f_0 * ki_261[k];

        t_150[k] = -2.0 * hi_66[k]
                   + f_0 * ki_262[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, hi_67, hi_68, hi_69, hi_70, hi_71, \
                         ki_263, ki_264, ki_265, ki_266, ki_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = -2.0 * hi_67[k]
                   + f_0 * ki_263[k];

        t_152[k] = -2.0 * hi_68[k]
                   + f_0 * ki_264[k];

        t_153[k] = -2.0 * hi_69[k]
                   + f_0 * ki_265[k];

        t_154[k] = -2.0 * hi_70[k]
                   + f_0 * ki_266[k];

        t_155[k] = -2.0 * hi_71[k]
                   + f_0 * ki_267[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, hi_72, hi_73, hi_74, hi_75, hi_76, \
                         ki_268, ki_269, ki_270, ki_271, ki_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = -2.0 * hi_72[k]
                   + f_0 * ki_268[k];

        t_157[k] = -2.0 * hi_73[k]
                   + f_0 * ki_269[k];

        t_158[k] = -2.0 * hi_74[k]
                   + f_0 * ki_270[k];

        t_159[k] = -2.0 * hi_75[k]
                   + f_0 * ki_271[k];

        t_160[k] = -2.0 * hi_76[k]
                   + f_0 * ki_272[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, hi_77, hi_78, hi_79, hi_80, hi_81, \
                         ki_273, ki_274, ki_275, ki_276, ki_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = -2.0 * hi_77[k]
                   + f_0 * ki_273[k];

        t_162[k] = -2.0 * hi_78[k]
                   + f_0 * ki_274[k];

        t_163[k] = -2.0 * hi_79[k]
                   + f_0 * ki_275[k];

        t_164[k] = -2.0 * hi_80[k]
                   + f_0 * ki_276[k];

        t_165[k] = -2.0 * hi_81[k]
                   + f_0 * ki_277[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, t_171, t_172, hi_82, hi_83, \
                         ki_278, ki_279, ki_308, ki_309, ki_310, ki_311, \
                         ki_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = -2.0 * hi_82[k]
                   + f_0 * ki_278[k];

        t_167[k] = -2.0 * hi_83[k]
                   + f_0 * ki_279[k];

        t_168[k] = f_0 * ki_308[k];

        t_169[k] = f_0 * ki_309[k];

        t_170[k] = f_0 * ki_310[k];

        t_171[k] = f_0 * ki_311[k];

        t_172[k] = f_0 * ki_312[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, t_179, t_180, ki_313, \
                         ki_314, ki_315, ki_316, ki_317, ki_318, ki_319, \
                         ki_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_0 * ki_313[k];

        t_174[k] = f_0 * ki_314[k];

        t_175[k] = f_0 * ki_315[k];

        t_176[k] = f_0 * ki_316[k];

        t_177[k] = f_0 * ki_317[k];

        t_178[k] = f_0 * ki_318[k];

        t_179[k] = f_0 * ki_319[k];

        t_180[k] = f_0 * ki_320[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);

    const auto *hi_84 = buffer.data(hi + 84);
    const auto *hi_85 = buffer.data(hi + 85);
    const auto *hi_86 = buffer.data(hi + 86);
    const auto *hi_87 = buffer.data(hi + 87);
    const auto *hi_88 = buffer.data(hi + 88);
    const auto *hi_89 = buffer.data(hi + 89);
    const auto *hi_90 = buffer.data(hi + 90);
    const auto *hi_91 = buffer.data(hi + 91);
    const auto *hi_92 = buffer.data(hi + 92);
    const auto *hi_93 = buffer.data(hi + 93);
    const auto *hi_94 = buffer.data(hi + 94);
    const auto *hi_95 = buffer.data(hi + 95);
    const auto *hi_96 = buffer.data(hi + 96);
    const auto *hi_97 = buffer.data(hi + 97);
    const auto *hi_98 = buffer.data(hi + 98);
    const auto *hi_99 = buffer.data(hi + 99);
    const auto *hi_100 = buffer.data(hi + 100);
    const auto *hi_101 = buffer.data(hi + 101);
    const auto *hi_102 = buffer.data(hi + 102);
    const auto *hi_103 = buffer.data(hi + 103);
    const auto *hi_104 = buffer.data(hi + 104);
    const auto *hi_105 = buffer.data(hi + 105);
    const auto *hi_106 = buffer.data(hi + 106);
    const auto *hi_107 = buffer.data(hi + 107);
    const auto *hi_108 = buffer.data(hi + 108);
    const auto *hi_109 = buffer.data(hi + 109);
    const auto *hi_110 = buffer.data(hi + 110);
    const auto *hi_111 = buffer.data(hi + 111);
    const auto *hi_112 = buffer.data(hi + 112);
    const auto *hi_113 = buffer.data(hi + 113);
    const auto *hi_114 = buffer.data(hi + 114);
    const auto *hi_115 = buffer.data(hi + 115);
    const auto *hi_116 = buffer.data(hi + 116);
    const auto *hi_117 = buffer.data(hi + 117);
    const auto *hi_118 = buffer.data(hi + 118);
    const auto *hi_119 = buffer.data(hi + 119);
    const auto *hi_120 = buffer.data(hi + 120);
    const auto *hi_121 = buffer.data(hi + 121);
    const auto *hi_122 = buffer.data(hi + 122);
    const auto *hi_123 = buffer.data(hi + 123);
    const auto *hi_124 = buffer.data(hi + 124);
    const auto *hi_125 = buffer.data(hi + 125);
    const auto *hi_126 = buffer.data(hi + 126);
    const auto *hi_127 = buffer.data(hi + 127);
    const auto *hi_128 = buffer.data(hi + 128);
    const auto *hi_129 = buffer.data(hi + 129);
    const auto *hi_130 = buffer.data(hi + 130);
    const auto *hi_131 = buffer.data(hi + 131);
    const auto *hi_132 = buffer.data(hi + 132);
    const auto *hi_133 = buffer.data(hi + 133);
    const auto *hi_134 = buffer.data(hi + 134);
    const auto *hi_135 = buffer.data(hi + 135);
    const auto *hi_136 = buffer.data(hi + 136);
    const auto *hi_137 = buffer.data(hi + 137);
    const auto *hi_138 = buffer.data(hi + 138);
    const auto *hi_139 = buffer.data(hi + 139);
    const auto *hi_140 = buffer.data(hi + 140);
    const auto *hi_141 = buffer.data(hi + 141);
    const auto *hi_142 = buffer.data(hi + 142);
    const auto *hi_143 = buffer.data(hi + 143);
    const auto *hi_144 = buffer.data(hi + 144);
    const auto *hi_145 = buffer.data(hi + 145);
    const auto *hi_146 = buffer.data(hi + 146);
    const auto *hi_147 = buffer.data(hi + 147);
    const auto *hi_148 = buffer.data(hi + 148);
    const auto *hi_149 = buffer.data(hi + 149);
    const auto *hi_150 = buffer.data(hi + 150);
    const auto *hi_151 = buffer.data(hi + 151);
    const auto *hi_152 = buffer.data(hi + 152);
    const auto *hi_153 = buffer.data(hi + 153);
    const auto *hi_154 = buffer.data(hi + 154);
    const auto *hi_155 = buffer.data(hi + 155);
    const auto *hi_156 = buffer.data(hi + 156);
    const auto *hi_157 = buffer.data(hi + 157);
    const auto *hi_158 = buffer.data(hi + 158);
    const auto *hi_159 = buffer.data(hi + 159);
    const auto *hi_160 = buffer.data(hi + 160);
    const auto *hi_161 = buffer.data(hi + 161);
    const auto *hi_162 = buffer.data(hi + 162);
    const auto *hi_163 = buffer.data(hi + 163);
    const auto *hi_164 = buffer.data(hi + 164);
    const auto *hi_165 = buffer.data(hi + 165);
    const auto *hi_166 = buffer.data(hi + 166);
    const auto *hi_167 = buffer.data(hi + 167);
    const auto *hi_168 = buffer.data(hi + 168);
    const auto *hi_169 = buffer.data(hi + 169);
    const auto *hi_170 = buffer.data(hi + 170);
    const auto *hi_171 = buffer.data(hi + 171);
    const auto *hi_172 = buffer.data(hi + 172);
    const auto *hi_173 = buffer.data(hi + 173);
    const auto *hi_174 = buffer.data(hi + 174);
    const auto *hi_175 = buffer.data(hi + 175);
    const auto *hi_176 = buffer.data(hi + 176);
    const auto *hi_177 = buffer.data(hi + 177);
    const auto *hi_178 = buffer.data(hi + 178);
    const auto *hi_179 = buffer.data(hi + 179);
    const auto *hi_180 = buffer.data(hi + 180);
    const auto *hi_181 = buffer.data(hi + 181);
    const auto *hi_182 = buffer.data(hi + 182);
    const auto *hi_183 = buffer.data(hi + 183);
    const auto *hi_184 = buffer.data(hi + 184);
    const auto *hi_185 = buffer.data(hi + 185);
    const auto *hi_186 = buffer.data(hi + 186);
    const auto *hi_187 = buffer.data(hi + 187);
    const auto *hi_188 = buffer.data(hi + 188);
    const auto *hi_189 = buffer.data(hi + 189);
    const auto *hi_190 = buffer.data(hi + 190);
    const auto *hi_191 = buffer.data(hi + 191);
    const auto *hi_192 = buffer.data(hi + 192);
    const auto *hi_193 = buffer.data(hi + 193);
    const auto *hi_194 = buffer.data(hi + 194);
    const auto *hi_195 = buffer.data(hi + 195);
    const auto *hi_196 = buffer.data(hi + 196);
    const auto *hi_197 = buffer.data(hi + 197);
    const auto *hi_198 = buffer.data(hi + 198);
    const auto *hi_199 = buffer.data(hi + 199);
    const auto *hi_200 = buffer.data(hi + 200);
    const auto *hi_201 = buffer.data(hi + 201);
    const auto *hi_202 = buffer.data(hi + 202);
    const auto *hi_203 = buffer.data(hi + 203);
    const auto *hi_204 = buffer.data(hi + 204);
    const auto *hi_205 = buffer.data(hi + 205);
    const auto *hi_206 = buffer.data(hi + 206);

    const auto *ki_321 = buffer.data(ki + 321);
    const auto *ki_322 = buffer.data(ki + 322);
    const auto *ki_323 = buffer.data(ki + 323);
    const auto *ki_324 = buffer.data(ki + 324);
    const auto *ki_325 = buffer.data(ki + 325);
    const auto *ki_326 = buffer.data(ki + 326);
    const auto *ki_327 = buffer.data(ki + 327);
    const auto *ki_328 = buffer.data(ki + 328);
    const auto *ki_329 = buffer.data(ki + 329);
    const auto *ki_330 = buffer.data(ki + 330);
    const auto *ki_331 = buffer.data(ki + 331);
    const auto *ki_332 = buffer.data(ki + 332);
    const auto *ki_333 = buffer.data(ki + 333);
    const auto *ki_334 = buffer.data(ki + 334);
    const auto *ki_335 = buffer.data(ki + 335);
    const auto *ki_336 = buffer.data(ki + 336);
    const auto *ki_337 = buffer.data(ki + 337);
    const auto *ki_338 = buffer.data(ki + 338);
    const auto *ki_339 = buffer.data(ki + 339);
    const auto *ki_340 = buffer.data(ki + 340);
    const auto *ki_341 = buffer.data(ki + 341);
    const auto *ki_342 = buffer.data(ki + 342);
    const auto *ki_343 = buffer.data(ki + 343);
    const auto *ki_344 = buffer.data(ki + 344);
    const auto *ki_345 = buffer.data(ki + 345);
    const auto *ki_346 = buffer.data(ki + 346);
    const auto *ki_347 = buffer.data(ki + 347);
    const auto *ki_348 = buffer.data(ki + 348);
    const auto *ki_349 = buffer.data(ki + 349);
    const auto *ki_350 = buffer.data(ki + 350);
    const auto *ki_351 = buffer.data(ki + 351);
    const auto *ki_352 = buffer.data(ki + 352);
    const auto *ki_353 = buffer.data(ki + 353);
    const auto *ki_354 = buffer.data(ki + 354);
    const auto *ki_355 = buffer.data(ki + 355);
    const auto *ki_356 = buffer.data(ki + 356);
    const auto *ki_357 = buffer.data(ki + 357);
    const auto *ki_358 = buffer.data(ki + 358);
    const auto *ki_359 = buffer.data(ki + 359);
    const auto *ki_360 = buffer.data(ki + 360);
    const auto *ki_361 = buffer.data(ki + 361);
    const auto *ki_362 = buffer.data(ki + 362);
    const auto *ki_363 = buffer.data(ki + 363);
    const auto *ki_364 = buffer.data(ki + 364);
    const auto *ki_365 = buffer.data(ki + 365);
    const auto *ki_366 = buffer.data(ki + 366);
    const auto *ki_367 = buffer.data(ki + 367);
    const auto *ki_368 = buffer.data(ki + 368);
    const auto *ki_369 = buffer.data(ki + 369);
    const auto *ki_370 = buffer.data(ki + 370);
    const auto *ki_371 = buffer.data(ki + 371);
    const auto *ki_372 = buffer.data(ki + 372);
    const auto *ki_373 = buffer.data(ki + 373);
    const auto *ki_374 = buffer.data(ki + 374);
    const auto *ki_375 = buffer.data(ki + 375);
    const auto *ki_376 = buffer.data(ki + 376);
    const auto *ki_377 = buffer.data(ki + 377);
    const auto *ki_378 = buffer.data(ki + 378);
    const auto *ki_379 = buffer.data(ki + 379);
    const auto *ki_380 = buffer.data(ki + 380);
    const auto *ki_381 = buffer.data(ki + 381);
    const auto *ki_382 = buffer.data(ki + 382);
    const auto *ki_383 = buffer.data(ki + 383);
    const auto *ki_384 = buffer.data(ki + 384);
    const auto *ki_385 = buffer.data(ki + 385);
    const auto *ki_386 = buffer.data(ki + 386);
    const auto *ki_387 = buffer.data(ki + 387);
    const auto *ki_388 = buffer.data(ki + 388);
    const auto *ki_389 = buffer.data(ki + 389);
    const auto *ki_390 = buffer.data(ki + 390);
    const auto *ki_391 = buffer.data(ki + 391);
    const auto *ki_392 = buffer.data(ki + 392);
    const auto *ki_393 = buffer.data(ki + 393);
    const auto *ki_394 = buffer.data(ki + 394);
    const auto *ki_395 = buffer.data(ki + 395);
    const auto *ki_396 = buffer.data(ki + 396);
    const auto *ki_397 = buffer.data(ki + 397);
    const auto *ki_398 = buffer.data(ki + 398);
    const auto *ki_399 = buffer.data(ki + 399);
    const auto *ki_400 = buffer.data(ki + 400);
    const auto *ki_401 = buffer.data(ki + 401);
    const auto *ki_402 = buffer.data(ki + 402);
    const auto *ki_403 = buffer.data(ki + 403);
    const auto *ki_404 = buffer.data(ki + 404);
    const auto *ki_405 = buffer.data(ki + 405);
    const auto *ki_406 = buffer.data(ki + 406);
    const auto *ki_407 = buffer.data(ki + 407);
    const auto *ki_408 = buffer.data(ki + 408);
    const auto *ki_409 = buffer.data(ki + 409);
    const auto *ki_410 = buffer.data(ki + 410);
    const auto *ki_411 = buffer.data(ki + 411);
    const auto *ki_412 = buffer.data(ki + 412);
    const auto *ki_413 = buffer.data(ki + 413);
    const auto *ki_414 = buffer.data(ki + 414);
    const auto *ki_415 = buffer.data(ki + 415);
    const auto *ki_416 = buffer.data(ki + 416);
    const auto *ki_417 = buffer.data(ki + 417);
    const auto *ki_418 = buffer.data(ki + 418);
    const auto *ki_419 = buffer.data(ki + 419);
    const auto *ki_448 = buffer.data(ki + 448);
    const auto *ki_449 = buffer.data(ki + 449);
    const auto *ki_450 = buffer.data(ki + 450);
    const auto *ki_451 = buffer.data(ki + 451);
    const auto *ki_452 = buffer.data(ki + 452);
    const auto *ki_453 = buffer.data(ki + 453);
    const auto *ki_454 = buffer.data(ki + 454);
    const auto *ki_455 = buffer.data(ki + 455);
    const auto *ki_456 = buffer.data(ki + 456);
    const auto *ki_457 = buffer.data(ki + 457);
    const auto *ki_458 = buffer.data(ki + 458);
    const auto *ki_459 = buffer.data(ki + 459);
    const auto *ki_460 = buffer.data(ki + 460);
    const auto *ki_461 = buffer.data(ki + 461);
    const auto *ki_462 = buffer.data(ki + 462);
    const auto *ki_463 = buffer.data(ki + 463);
    const auto *ki_464 = buffer.data(ki + 464);
    const auto *ki_465 = buffer.data(ki + 465);
    const auto *ki_466 = buffer.data(ki + 466);
    const auto *ki_467 = buffer.data(ki + 467);
    const auto *ki_468 = buffer.data(ki + 468);
    const auto *ki_469 = buffer.data(ki + 469);
    const auto *ki_470 = buffer.data(ki + 470);
    const auto *ki_471 = buffer.data(ki + 471);
    const auto *ki_472 = buffer.data(ki + 472);
    const auto *ki_473 = buffer.data(ki + 473);
    const auto *ki_474 = buffer.data(ki + 474);
    const auto *ki_475 = buffer.data(ki + 475);
    const auto *ki_476 = buffer.data(ki + 476);
    const auto *ki_477 = buffer.data(ki + 477);
    const auto *ki_478 = buffer.data(ki + 478);
    const auto *ki_479 = buffer.data(ki + 479);
    const auto *ki_480 = buffer.data(ki + 480);
    const auto *ki_481 = buffer.data(ki + 481);
    const auto *ki_482 = buffer.data(ki + 482);
    const auto *ki_483 = buffer.data(ki + 483);
    const auto *ki_484 = buffer.data(ki + 484);
    const auto *ki_485 = buffer.data(ki + 485);
    const auto *ki_486 = buffer.data(ki + 486);
    const auto *ki_487 = buffer.data(ki + 487);
    const auto *ki_488 = buffer.data(ki + 488);
    const auto *ki_489 = buffer.data(ki + 489);
    const auto *ki_490 = buffer.data(ki + 490);
    const auto *ki_491 = buffer.data(ki + 491);
    const auto *ki_492 = buffer.data(ki + 492);
    const auto *ki_493 = buffer.data(ki + 493);
    const auto *ki_494 = buffer.data(ki + 494);
    const auto *ki_495 = buffer.data(ki + 495);
    const auto *ki_496 = buffer.data(ki + 496);
    const auto *ki_497 = buffer.data(ki + 497);
    const auto *ki_498 = buffer.data(ki + 498);
    const auto *ki_499 = buffer.data(ki + 499);
    const auto *ki_500 = buffer.data(ki + 500);
    const auto *ki_501 = buffer.data(ki + 501);
    const auto *ki_502 = buffer.data(ki + 502);
    const auto *ki_503 = buffer.data(ki + 503);
    const auto *ki_504 = buffer.data(ki + 504);
    const auto *ki_505 = buffer.data(ki + 505);
    const auto *ki_506 = buffer.data(ki + 506);
    const auto *ki_507 = buffer.data(ki + 507);
    const auto *ki_508 = buffer.data(ki + 508);
    const auto *ki_509 = buffer.data(ki + 509);
    const auto *ki_510 = buffer.data(ki + 510);
    const auto *ki_511 = buffer.data(ki + 511);
    const auto *ki_512 = buffer.data(ki + 512);
    const auto *ki_513 = buffer.data(ki + 513);
    const auto *ki_514 = buffer.data(ki + 514);

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, t_186, t_187, t_188, ki_321, \
                         ki_322, ki_323, ki_324, ki_325, ki_326, ki_327, \
                         ki_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_0 * ki_321[k];

        t_182[k] = f_0 * ki_322[k];

        t_183[k] = f_0 * ki_323[k];

        t_184[k] = f_0 * ki_324[k];

        t_185[k] = f_0 * ki_325[k];

        t_186[k] = f_0 * ki_326[k];

        t_187[k] = f_0 * ki_327[k];

        t_188[k] = f_0 * ki_328[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, ki_329, ki_330, \
                         ki_331, ki_332, ki_333, ki_334, ki_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * ki_329[k];

        t_190[k] = f_0 * ki_330[k];

        t_191[k] = f_0 * ki_331[k];

        t_192[k] = f_0 * ki_332[k];

        t_193[k] = f_0 * ki_333[k];

        t_194[k] = f_0 * ki_334[k];

        t_195[k] = f_0 * ki_335[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, hi_84, hi_85, hi_86, hi_87, hi_88, \
                         ki_336, ki_337, ki_338, ki_339, ki_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -hi_84[k]
                   + f_0 * ki_336[k];

        t_197[k] = -hi_85[k]
                   + f_0 * ki_337[k];

        t_198[k] = -hi_86[k]
                   + f_0 * ki_338[k];

        t_199[k] = -hi_87[k]
                   + f_0 * ki_339[k];

        t_200[k] = -hi_88[k]
                   + f_0 * ki_340[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, hi_89, hi_90, hi_91, hi_92, hi_93, \
                         ki_341, ki_342, ki_343, ki_344, ki_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = -hi_89[k]
                   + f_0 * ki_341[k];

        t_202[k] = -hi_90[k]
                   + f_0 * ki_342[k];

        t_203[k] = -hi_91[k]
                   + f_0 * ki_343[k];

        t_204[k] = -hi_92[k]
                   + f_0 * ki_344[k];

        t_205[k] = -hi_93[k]
                   + f_0 * ki_345[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, hi_94, hi_95, hi_96, hi_97, hi_98, \
                         ki_346, ki_347, ki_348, ki_349, ki_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = -hi_94[k]
                   + f_0 * ki_346[k];

        t_207[k] = -hi_95[k]
                   + f_0 * ki_347[k];

        t_208[k] = -hi_96[k]
                   + f_0 * ki_348[k];

        t_209[k] = -hi_97[k]
                   + f_0 * ki_349[k];

        t_210[k] = -hi_98[k]
                   + f_0 * ki_350[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, hi_99, hi_100, hi_101, hi_102, \
                         hi_103, ki_351, ki_352, ki_353, ki_354, \
                         ki_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = -hi_99[k]
                   + f_0 * ki_351[k];

        t_212[k] = -hi_100[k]
                   + f_0 * ki_352[k];

        t_213[k] = -hi_101[k]
                   + f_0 * ki_353[k];

        t_214[k] = -hi_102[k]
                   + f_0 * ki_354[k];

        t_215[k] = -hi_103[k]
                   + f_0 * ki_355[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, hi_104, hi_105, hi_106, hi_107, \
                         hi_108, ki_356, ki_357, ki_358, ki_359, \
                         ki_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = -hi_104[k]
                   + f_0 * ki_356[k];

        t_217[k] = -hi_105[k]
                   + f_0 * ki_357[k];

        t_218[k] = -hi_106[k]
                   + f_0 * ki_358[k];

        t_219[k] = -hi_107[k]
                   + f_0 * ki_359[k];

        t_220[k] = -hi_108[k]
                   + f_0 * ki_360[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, hi_109, hi_110, hi_111, hi_112, \
                         hi_113, ki_361, ki_362, ki_363, ki_364, \
                         ki_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = -hi_109[k]
                   + f_0 * ki_361[k];

        t_222[k] = -hi_110[k]
                   + f_0 * ki_362[k];

        t_223[k] = -hi_111[k]
                   + f_0 * ki_363[k];

        t_224[k] = -2.0 * hi_112[k]
                   + f_0 * ki_364[k];

        t_225[k] = -2.0 * hi_113[k]
                   + f_0 * ki_365[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, t_230, hi_114, hi_115, hi_116, hi_117, \
                         hi_118, ki_366, ki_367, ki_368, ki_369, \
                         ki_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = -2.0 * hi_114[k]
                   + f_0 * ki_366[k];

        t_227[k] = -2.0 * hi_115[k]
                   + f_0 * ki_367[k];

        t_228[k] = -2.0 * hi_116[k]
                   + f_0 * ki_368[k];

        t_229[k] = -2.0 * hi_117[k]
                   + f_0 * ki_369[k];

        t_230[k] = -2.0 * hi_118[k]
                   + f_0 * ki_370[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, hi_119, hi_120, hi_121, hi_122, \
                         hi_123, ki_371, ki_372, ki_373, ki_374, \
                         ki_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = -2.0 * hi_119[k]
                   + f_0 * ki_371[k];

        t_232[k] = -2.0 * hi_120[k]
                   + f_0 * ki_372[k];

        t_233[k] = -2.0 * hi_121[k]
                   + f_0 * ki_373[k];

        t_234[k] = -2.0 * hi_122[k]
                   + f_0 * ki_374[k];

        t_235[k] = -2.0 * hi_123[k]
                   + f_0 * ki_375[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, hi_124, hi_125, hi_126, hi_127, \
                         hi_128, ki_376, ki_377, ki_378, ki_379, \
                         ki_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = -2.0 * hi_124[k]
                   + f_0 * ki_376[k];

        t_237[k] = -2.0 * hi_125[k]
                   + f_0 * ki_377[k];

        t_238[k] = -2.0 * hi_126[k]
                   + f_0 * ki_378[k];

        t_239[k] = -2.0 * hi_127[k]
                   + f_0 * ki_379[k];

        t_240[k] = -2.0 * hi_128[k]
                   + f_0 * ki_380[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, hi_129, hi_130, hi_131, hi_132, \
                         hi_133, ki_381, ki_382, ki_383, ki_384, \
                         ki_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = -2.0 * hi_129[k]
                   + f_0 * ki_381[k];

        t_242[k] = -2.0 * hi_130[k]
                   + f_0 * ki_382[k];

        t_243[k] = -2.0 * hi_131[k]
                   + f_0 * ki_383[k];

        t_244[k] = -2.0 * hi_132[k]
                   + f_0 * ki_384[k];

        t_245[k] = -2.0 * hi_133[k]
                   + f_0 * ki_385[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, hi_134, hi_135, hi_136, hi_137, \
                         hi_138, ki_386, ki_387, ki_388, ki_389, \
                         ki_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = -2.0 * hi_134[k]
                   + f_0 * ki_386[k];

        t_247[k] = -2.0 * hi_135[k]
                   + f_0 * ki_387[k];

        t_248[k] = -2.0 * hi_136[k]
                   + f_0 * ki_388[k];

        t_249[k] = -2.0 * hi_137[k]
                   + f_0 * ki_389[k];

        t_250[k] = -2.0 * hi_138[k]
                   + f_0 * ki_390[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, hi_139, hi_140, hi_141, hi_142, \
                         hi_143, ki_391, ki_392, ki_393, ki_394, \
                         ki_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = -2.0 * hi_139[k]
                   + f_0 * ki_391[k];

        t_252[k] = -3.0 * hi_140[k]
                   + f_0 * ki_392[k];

        t_253[k] = -3.0 * hi_141[k]
                   + f_0 * ki_393[k];

        t_254[k] = -3.0 * hi_142[k]
                   + f_0 * ki_394[k];

        t_255[k] = -3.0 * hi_143[k]
                   + f_0 * ki_395[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, hi_144, hi_145, hi_146, hi_147, \
                         hi_148, ki_396, ki_397, ki_398, ki_399, \
                         ki_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = -3.0 * hi_144[k]
                   + f_0 * ki_396[k];

        t_257[k] = -3.0 * hi_145[k]
                   + f_0 * ki_397[k];

        t_258[k] = -3.0 * hi_146[k]
                   + f_0 * ki_398[k];

        t_259[k] = -3.0 * hi_147[k]
                   + f_0 * ki_399[k];

        t_260[k] = -3.0 * hi_148[k]
                   + f_0 * ki_400[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, hi_149, hi_150, hi_151, hi_152, \
                         hi_153, ki_401, ki_402, ki_403, ki_404, \
                         ki_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = -3.0 * hi_149[k]
                   + f_0 * ki_401[k];

        t_262[k] = -3.0 * hi_150[k]
                   + f_0 * ki_402[k];

        t_263[k] = -3.0 * hi_151[k]
                   + f_0 * ki_403[k];

        t_264[k] = -3.0 * hi_152[k]
                   + f_0 * ki_404[k];

        t_265[k] = -3.0 * hi_153[k]
                   + f_0 * ki_405[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, hi_154, hi_155, hi_156, hi_157, \
                         hi_158, ki_406, ki_407, ki_408, ki_409, \
                         ki_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = -3.0 * hi_154[k]
                   + f_0 * ki_406[k];

        t_267[k] = -3.0 * hi_155[k]
                   + f_0 * ki_407[k];

        t_268[k] = -3.0 * hi_156[k]
                   + f_0 * ki_408[k];

        t_269[k] = -3.0 * hi_157[k]
                   + f_0 * ki_409[k];

        t_270[k] = -3.0 * hi_158[k]
                   + f_0 * ki_410[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, hi_159, hi_160, hi_161, hi_162, \
                         hi_163, ki_411, ki_412, ki_413, ki_414, \
                         ki_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = -3.0 * hi_159[k]
                   + f_0 * ki_411[k];

        t_272[k] = -3.0 * hi_160[k]
                   + f_0 * ki_412[k];

        t_273[k] = -3.0 * hi_161[k]
                   + f_0 * ki_413[k];

        t_274[k] = -3.0 * hi_162[k]
                   + f_0 * ki_414[k];

        t_275[k] = -3.0 * hi_163[k]
                   + f_0 * ki_415[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, t_281, hi_164, hi_165, hi_166, \
                         hi_167, ki_416, ki_417, ki_418, ki_419, ki_448, \
                         ki_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = -3.0 * hi_164[k]
                   + f_0 * ki_416[k];

        t_277[k] = -3.0 * hi_165[k]
                   + f_0 * ki_417[k];

        t_278[k] = -3.0 * hi_166[k]
                   + f_0 * ki_418[k];

        t_279[k] = -3.0 * hi_167[k]
                   + f_0 * ki_419[k];

        t_280[k] = f_0 * ki_448[k];

        t_281[k] = f_0 * ki_449[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, t_287, t_288, t_289, ki_450, \
                         ki_451, ki_452, ki_453, ki_454, ki_455, ki_456, \
                         ki_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_0 * ki_450[k];

        t_283[k] = f_0 * ki_451[k];

        t_284[k] = f_0 * ki_452[k];

        t_285[k] = f_0 * ki_453[k];

        t_286[k] = f_0 * ki_454[k];

        t_287[k] = f_0 * ki_455[k];

        t_288[k] = f_0 * ki_456[k];

        t_289[k] = f_0 * ki_457[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, t_295, t_296, t_297, ki_458, \
                         ki_459, ki_460, ki_461, ki_462, ki_463, ki_464, \
                         ki_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_0 * ki_458[k];

        t_291[k] = f_0 * ki_459[k];

        t_292[k] = f_0 * ki_460[k];

        t_293[k] = f_0 * ki_461[k];

        t_294[k] = f_0 * ki_462[k];

        t_295[k] = f_0 * ki_463[k];

        t_296[k] = f_0 * ki_464[k];

        t_297[k] = f_0 * ki_465[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, t_303, t_304, t_305, ki_466, \
                         ki_467, ki_468, ki_469, ki_470, ki_471, ki_472, \
                         ki_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_0 * ki_466[k];

        t_299[k] = f_0 * ki_467[k];

        t_300[k] = f_0 * ki_468[k];

        t_301[k] = f_0 * ki_469[k];

        t_302[k] = f_0 * ki_470[k];

        t_303[k] = f_0 * ki_471[k];

        t_304[k] = f_0 * ki_472[k];

        t_305[k] = f_0 * ki_473[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, t_311, hi_168, hi_169, hi_170, \
                         hi_171, ki_474, ki_475, ki_476, ki_477, ki_478, \
                         ki_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_0 * ki_474[k];

        t_307[k] = f_0 * ki_475[k];

        t_308[k] = -hi_168[k]
                   + f_0 * ki_476[k];

        t_309[k] = -hi_169[k]
                   + f_0 * ki_477[k];

        t_310[k] = -hi_170[k]
                   + f_0 * ki_478[k];

        t_311[k] = -hi_171[k]
                   + f_0 * ki_479[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, hi_172, hi_173, hi_174, hi_175, \
                         hi_176, ki_480, ki_481, ki_482, ki_483, \
                         ki_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -hi_172[k]
                   + f_0 * ki_480[k];

        t_313[k] = -hi_173[k]
                   + f_0 * ki_481[k];

        t_314[k] = -hi_174[k]
                   + f_0 * ki_482[k];

        t_315[k] = -hi_175[k]
                   + f_0 * ki_483[k];

        t_316[k] = -hi_176[k]
                   + f_0 * ki_484[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, hi_177, hi_178, hi_179, hi_180, \
                         hi_181, ki_485, ki_486, ki_487, ki_488, \
                         ki_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = -hi_177[k]
                   + f_0 * ki_485[k];

        t_318[k] = -hi_178[k]
                   + f_0 * ki_486[k];

        t_319[k] = -hi_179[k]
                   + f_0 * ki_487[k];

        t_320[k] = -hi_180[k]
                   + f_0 * ki_488[k];

        t_321[k] = -hi_181[k]
                   + f_0 * ki_489[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, hi_182, hi_183, hi_184, hi_185, \
                         hi_186, ki_490, ki_491, ki_492, ki_493, \
                         ki_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = -hi_182[k]
                   + f_0 * ki_490[k];

        t_323[k] = -hi_183[k]
                   + f_0 * ki_491[k];

        t_324[k] = -hi_184[k]
                   + f_0 * ki_492[k];

        t_325[k] = -hi_185[k]
                   + f_0 * ki_493[k];

        t_326[k] = -hi_186[k]
                   + f_0 * ki_494[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, hi_187, hi_188, hi_189, hi_190, \
                         hi_191, ki_495, ki_496, ki_497, ki_498, \
                         ki_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = -hi_187[k]
                   + f_0 * ki_495[k];

        t_328[k] = -hi_188[k]
                   + f_0 * ki_496[k];

        t_329[k] = -hi_189[k]
                   + f_0 * ki_497[k];

        t_330[k] = -hi_190[k]
                   + f_0 * ki_498[k];

        t_331[k] = -hi_191[k]
                   + f_0 * ki_499[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, hi_192, hi_193, hi_194, hi_195, \
                         hi_196, ki_500, ki_501, ki_502, ki_503, \
                         ki_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -hi_192[k]
                   + f_0 * ki_500[k];

        t_333[k] = -hi_193[k]
                   + f_0 * ki_501[k];

        t_334[k] = -hi_194[k]
                   + f_0 * ki_502[k];

        t_335[k] = -hi_195[k]
                   + f_0 * ki_503[k];

        t_336[k] = -2.0 * hi_196[k]
                   + f_0 * ki_504[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, hi_197, hi_198, hi_199, hi_200, \
                         hi_201, ki_505, ki_506, ki_507, ki_508, \
                         ki_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -2.0 * hi_197[k]
                   + f_0 * ki_505[k];

        t_338[k] = -2.0 * hi_198[k]
                   + f_0 * ki_506[k];

        t_339[k] = -2.0 * hi_199[k]
                   + f_0 * ki_507[k];

        t_340[k] = -2.0 * hi_200[k]
                   + f_0 * ki_508[k];

        t_341[k] = -2.0 * hi_201[k]
                   + f_0 * ki_509[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, hi_202, hi_203, hi_204, hi_205, \
                         hi_206, ki_510, ki_511, ki_512, ki_513, \
                         ki_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -2.0 * hi_202[k]
                   + f_0 * ki_510[k];

        t_343[k] = -2.0 * hi_203[k]
                   + f_0 * ki_511[k];

        t_344[k] = -2.0 * hi_204[k]
                   + f_0 * ki_512[k];

        t_345[k] = -2.0 * hi_205[k]
                   + f_0 * ki_513[k];

        t_346[k] = -2.0 * hi_206[k]
                   + f_0 * ki_514[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hi_207 = buffer.data(hi + 207);
    const auto *hi_208 = buffer.data(hi + 208);
    const auto *hi_209 = buffer.data(hi + 209);
    const auto *hi_210 = buffer.data(hi + 210);
    const auto *hi_211 = buffer.data(hi + 211);
    const auto *hi_212 = buffer.data(hi + 212);
    const auto *hi_213 = buffer.data(hi + 213);
    const auto *hi_214 = buffer.data(hi + 214);
    const auto *hi_215 = buffer.data(hi + 215);
    const auto *hi_216 = buffer.data(hi + 216);
    const auto *hi_217 = buffer.data(hi + 217);
    const auto *hi_218 = buffer.data(hi + 218);
    const auto *hi_219 = buffer.data(hi + 219);
    const auto *hi_220 = buffer.data(hi + 220);
    const auto *hi_221 = buffer.data(hi + 221);
    const auto *hi_222 = buffer.data(hi + 222);
    const auto *hi_223 = buffer.data(hi + 223);
    const auto *hi_224 = buffer.data(hi + 224);
    const auto *hi_225 = buffer.data(hi + 225);
    const auto *hi_226 = buffer.data(hi + 226);
    const auto *hi_227 = buffer.data(hi + 227);
    const auto *hi_228 = buffer.data(hi + 228);
    const auto *hi_229 = buffer.data(hi + 229);
    const auto *hi_230 = buffer.data(hi + 230);
    const auto *hi_231 = buffer.data(hi + 231);
    const auto *hi_232 = buffer.data(hi + 232);
    const auto *hi_233 = buffer.data(hi + 233);
    const auto *hi_234 = buffer.data(hi + 234);
    const auto *hi_235 = buffer.data(hi + 235);
    const auto *hi_236 = buffer.data(hi + 236);
    const auto *hi_237 = buffer.data(hi + 237);
    const auto *hi_238 = buffer.data(hi + 238);
    const auto *hi_239 = buffer.data(hi + 239);
    const auto *hi_240 = buffer.data(hi + 240);
    const auto *hi_241 = buffer.data(hi + 241);
    const auto *hi_242 = buffer.data(hi + 242);
    const auto *hi_243 = buffer.data(hi + 243);
    const auto *hi_244 = buffer.data(hi + 244);
    const auto *hi_245 = buffer.data(hi + 245);
    const auto *hi_246 = buffer.data(hi + 246);
    const auto *hi_247 = buffer.data(hi + 247);
    const auto *hi_248 = buffer.data(hi + 248);
    const auto *hi_249 = buffer.data(hi + 249);
    const auto *hi_250 = buffer.data(hi + 250);
    const auto *hi_251 = buffer.data(hi + 251);
    const auto *hi_252 = buffer.data(hi + 252);
    const auto *hi_253 = buffer.data(hi + 253);
    const auto *hi_254 = buffer.data(hi + 254);
    const auto *hi_255 = buffer.data(hi + 255);
    const auto *hi_256 = buffer.data(hi + 256);
    const auto *hi_257 = buffer.data(hi + 257);
    const auto *hi_258 = buffer.data(hi + 258);
    const auto *hi_259 = buffer.data(hi + 259);
    const auto *hi_260 = buffer.data(hi + 260);
    const auto *hi_261 = buffer.data(hi + 261);
    const auto *hi_262 = buffer.data(hi + 262);
    const auto *hi_263 = buffer.data(hi + 263);
    const auto *hi_264 = buffer.data(hi + 264);
    const auto *hi_265 = buffer.data(hi + 265);
    const auto *hi_266 = buffer.data(hi + 266);
    const auto *hi_267 = buffer.data(hi + 267);
    const auto *hi_268 = buffer.data(hi + 268);
    const auto *hi_269 = buffer.data(hi + 269);
    const auto *hi_270 = buffer.data(hi + 270);
    const auto *hi_271 = buffer.data(hi + 271);
    const auto *hi_272 = buffer.data(hi + 272);
    const auto *hi_273 = buffer.data(hi + 273);
    const auto *hi_274 = buffer.data(hi + 274);
    const auto *hi_275 = buffer.data(hi + 275);
    const auto *hi_276 = buffer.data(hi + 276);
    const auto *hi_277 = buffer.data(hi + 277);
    const auto *hi_278 = buffer.data(hi + 278);
    const auto *hi_279 = buffer.data(hi + 279);
    const auto *hi_280 = buffer.data(hi + 280);
    const auto *hi_281 = buffer.data(hi + 281);
    const auto *hi_282 = buffer.data(hi + 282);
    const auto *hi_283 = buffer.data(hi + 283);
    const auto *hi_284 = buffer.data(hi + 284);
    const auto *hi_285 = buffer.data(hi + 285);
    const auto *hi_286 = buffer.data(hi + 286);
    const auto *hi_287 = buffer.data(hi + 287);
    const auto *hi_288 = buffer.data(hi + 288);
    const auto *hi_289 = buffer.data(hi + 289);
    const auto *hi_290 = buffer.data(hi + 290);
    const auto *hi_291 = buffer.data(hi + 291);
    const auto *hi_292 = buffer.data(hi + 292);
    const auto *hi_293 = buffer.data(hi + 293);
    const auto *hi_294 = buffer.data(hi + 294);
    const auto *hi_295 = buffer.data(hi + 295);
    const auto *hi_296 = buffer.data(hi + 296);
    const auto *hi_297 = buffer.data(hi + 297);
    const auto *hi_298 = buffer.data(hi + 298);
    const auto *hi_299 = buffer.data(hi + 299);
    const auto *hi_300 = buffer.data(hi + 300);
    const auto *hi_301 = buffer.data(hi + 301);
    const auto *hi_302 = buffer.data(hi + 302);
    const auto *hi_303 = buffer.data(hi + 303);
    const auto *hi_304 = buffer.data(hi + 304);
    const auto *hi_305 = buffer.data(hi + 305);
    const auto *hi_306 = buffer.data(hi + 306);
    const auto *hi_307 = buffer.data(hi + 307);
    const auto *hi_308 = buffer.data(hi + 308);
    const auto *hi_309 = buffer.data(hi + 309);
    const auto *hi_310 = buffer.data(hi + 310);
    const auto *hi_311 = buffer.data(hi + 311);
    const auto *hi_312 = buffer.data(hi + 312);
    const auto *hi_313 = buffer.data(hi + 313);
    const auto *hi_314 = buffer.data(hi + 314);
    const auto *hi_315 = buffer.data(hi + 315);
    const auto *hi_316 = buffer.data(hi + 316);
    const auto *hi_317 = buffer.data(hi + 317);
    const auto *hi_318 = buffer.data(hi + 318);
    const auto *hi_319 = buffer.data(hi + 319);
    const auto *hi_320 = buffer.data(hi + 320);
    const auto *hi_321 = buffer.data(hi + 321);
    const auto *hi_322 = buffer.data(hi + 322);
    const auto *hi_323 = buffer.data(hi + 323);
    const auto *hi_324 = buffer.data(hi + 324);
    const auto *hi_325 = buffer.data(hi + 325);
    const auto *hi_326 = buffer.data(hi + 326);
    const auto *hi_327 = buffer.data(hi + 327);
    const auto *hi_328 = buffer.data(hi + 328);
    const auto *hi_329 = buffer.data(hi + 329);
    const auto *hi_330 = buffer.data(hi + 330);
    const auto *hi_331 = buffer.data(hi + 331);
    const auto *hi_332 = buffer.data(hi + 332);
    const auto *hi_333 = buffer.data(hi + 333);
    const auto *hi_334 = buffer.data(hi + 334);
    const auto *hi_335 = buffer.data(hi + 335);
    const auto *hi_336 = buffer.data(hi + 336);
    const auto *hi_337 = buffer.data(hi + 337);
    const auto *hi_338 = buffer.data(hi + 338);

    const auto *ki_515 = buffer.data(ki + 515);
    const auto *ki_516 = buffer.data(ki + 516);
    const auto *ki_517 = buffer.data(ki + 517);
    const auto *ki_518 = buffer.data(ki + 518);
    const auto *ki_519 = buffer.data(ki + 519);
    const auto *ki_520 = buffer.data(ki + 520);
    const auto *ki_521 = buffer.data(ki + 521);
    const auto *ki_522 = buffer.data(ki + 522);
    const auto *ki_523 = buffer.data(ki + 523);
    const auto *ki_524 = buffer.data(ki + 524);
    const auto *ki_525 = buffer.data(ki + 525);
    const auto *ki_526 = buffer.data(ki + 526);
    const auto *ki_527 = buffer.data(ki + 527);
    const auto *ki_528 = buffer.data(ki + 528);
    const auto *ki_529 = buffer.data(ki + 529);
    const auto *ki_530 = buffer.data(ki + 530);
    const auto *ki_531 = buffer.data(ki + 531);
    const auto *ki_532 = buffer.data(ki + 532);
    const auto *ki_533 = buffer.data(ki + 533);
    const auto *ki_534 = buffer.data(ki + 534);
    const auto *ki_535 = buffer.data(ki + 535);
    const auto *ki_536 = buffer.data(ki + 536);
    const auto *ki_537 = buffer.data(ki + 537);
    const auto *ki_538 = buffer.data(ki + 538);
    const auto *ki_539 = buffer.data(ki + 539);
    const auto *ki_540 = buffer.data(ki + 540);
    const auto *ki_541 = buffer.data(ki + 541);
    const auto *ki_542 = buffer.data(ki + 542);
    const auto *ki_543 = buffer.data(ki + 543);
    const auto *ki_544 = buffer.data(ki + 544);
    const auto *ki_545 = buffer.data(ki + 545);
    const auto *ki_546 = buffer.data(ki + 546);
    const auto *ki_547 = buffer.data(ki + 547);
    const auto *ki_548 = buffer.data(ki + 548);
    const auto *ki_549 = buffer.data(ki + 549);
    const auto *ki_550 = buffer.data(ki + 550);
    const auto *ki_551 = buffer.data(ki + 551);
    const auto *ki_552 = buffer.data(ki + 552);
    const auto *ki_553 = buffer.data(ki + 553);
    const auto *ki_554 = buffer.data(ki + 554);
    const auto *ki_555 = buffer.data(ki + 555);
    const auto *ki_556 = buffer.data(ki + 556);
    const auto *ki_557 = buffer.data(ki + 557);
    const auto *ki_558 = buffer.data(ki + 558);
    const auto *ki_559 = buffer.data(ki + 559);
    const auto *ki_560 = buffer.data(ki + 560);
    const auto *ki_561 = buffer.data(ki + 561);
    const auto *ki_562 = buffer.data(ki + 562);
    const auto *ki_563 = buffer.data(ki + 563);
    const auto *ki_564 = buffer.data(ki + 564);
    const auto *ki_565 = buffer.data(ki + 565);
    const auto *ki_566 = buffer.data(ki + 566);
    const auto *ki_567 = buffer.data(ki + 567);
    const auto *ki_568 = buffer.data(ki + 568);
    const auto *ki_569 = buffer.data(ki + 569);
    const auto *ki_570 = buffer.data(ki + 570);
    const auto *ki_571 = buffer.data(ki + 571);
    const auto *ki_572 = buffer.data(ki + 572);
    const auto *ki_573 = buffer.data(ki + 573);
    const auto *ki_574 = buffer.data(ki + 574);
    const auto *ki_575 = buffer.data(ki + 575);
    const auto *ki_576 = buffer.data(ki + 576);
    const auto *ki_577 = buffer.data(ki + 577);
    const auto *ki_578 = buffer.data(ki + 578);
    const auto *ki_579 = buffer.data(ki + 579);
    const auto *ki_580 = buffer.data(ki + 580);
    const auto *ki_581 = buffer.data(ki + 581);
    const auto *ki_582 = buffer.data(ki + 582);
    const auto *ki_583 = buffer.data(ki + 583);
    const auto *ki_584 = buffer.data(ki + 584);
    const auto *ki_585 = buffer.data(ki + 585);
    const auto *ki_586 = buffer.data(ki + 586);
    const auto *ki_587 = buffer.data(ki + 587);
    const auto *ki_616 = buffer.data(ki + 616);
    const auto *ki_617 = buffer.data(ki + 617);
    const auto *ki_618 = buffer.data(ki + 618);
    const auto *ki_619 = buffer.data(ki + 619);
    const auto *ki_620 = buffer.data(ki + 620);
    const auto *ki_621 = buffer.data(ki + 621);
    const auto *ki_622 = buffer.data(ki + 622);
    const auto *ki_623 = buffer.data(ki + 623);
    const auto *ki_624 = buffer.data(ki + 624);
    const auto *ki_625 = buffer.data(ki + 625);
    const auto *ki_626 = buffer.data(ki + 626);
    const auto *ki_627 = buffer.data(ki + 627);
    const auto *ki_628 = buffer.data(ki + 628);
    const auto *ki_629 = buffer.data(ki + 629);
    const auto *ki_630 = buffer.data(ki + 630);
    const auto *ki_631 = buffer.data(ki + 631);
    const auto *ki_632 = buffer.data(ki + 632);
    const auto *ki_633 = buffer.data(ki + 633);
    const auto *ki_634 = buffer.data(ki + 634);
    const auto *ki_635 = buffer.data(ki + 635);
    const auto *ki_636 = buffer.data(ki + 636);
    const auto *ki_637 = buffer.data(ki + 637);
    const auto *ki_638 = buffer.data(ki + 638);
    const auto *ki_639 = buffer.data(ki + 639);
    const auto *ki_640 = buffer.data(ki + 640);
    const auto *ki_641 = buffer.data(ki + 641);
    const auto *ki_642 = buffer.data(ki + 642);
    const auto *ki_643 = buffer.data(ki + 643);
    const auto *ki_644 = buffer.data(ki + 644);
    const auto *ki_645 = buffer.data(ki + 645);
    const auto *ki_646 = buffer.data(ki + 646);
    const auto *ki_647 = buffer.data(ki + 647);
    const auto *ki_648 = buffer.data(ki + 648);
    const auto *ki_649 = buffer.data(ki + 649);
    const auto *ki_650 = buffer.data(ki + 650);
    const auto *ki_651 = buffer.data(ki + 651);
    const auto *ki_652 = buffer.data(ki + 652);
    const auto *ki_653 = buffer.data(ki + 653);
    const auto *ki_654 = buffer.data(ki + 654);
    const auto *ki_655 = buffer.data(ki + 655);
    const auto *ki_656 = buffer.data(ki + 656);
    const auto *ki_657 = buffer.data(ki + 657);
    const auto *ki_658 = buffer.data(ki + 658);
    const auto *ki_659 = buffer.data(ki + 659);
    const auto *ki_660 = buffer.data(ki + 660);
    const auto *ki_661 = buffer.data(ki + 661);
    const auto *ki_662 = buffer.data(ki + 662);
    const auto *ki_663 = buffer.data(ki + 663);
    const auto *ki_664 = buffer.data(ki + 664);
    const auto *ki_665 = buffer.data(ki + 665);
    const auto *ki_666 = buffer.data(ki + 666);
    const auto *ki_667 = buffer.data(ki + 667);
    const auto *ki_668 = buffer.data(ki + 668);
    const auto *ki_669 = buffer.data(ki + 669);
    const auto *ki_670 = buffer.data(ki + 670);
    const auto *ki_671 = buffer.data(ki + 671);
    const auto *ki_672 = buffer.data(ki + 672);
    const auto *ki_673 = buffer.data(ki + 673);
    const auto *ki_674 = buffer.data(ki + 674);
    const auto *ki_675 = buffer.data(ki + 675);
    const auto *ki_676 = buffer.data(ki + 676);
    const auto *ki_677 = buffer.data(ki + 677);
    const auto *ki_678 = buffer.data(ki + 678);
    const auto *ki_679 = buffer.data(ki + 679);
    const auto *ki_680 = buffer.data(ki + 680);
    const auto *ki_681 = buffer.data(ki + 681);
    const auto *ki_682 = buffer.data(ki + 682);
    const auto *ki_683 = buffer.data(ki + 683);
    const auto *ki_684 = buffer.data(ki + 684);
    const auto *ki_685 = buffer.data(ki + 685);
    const auto *ki_686 = buffer.data(ki + 686);
    const auto *ki_687 = buffer.data(ki + 687);
    const auto *ki_688 = buffer.data(ki + 688);
    const auto *ki_689 = buffer.data(ki + 689);
    const auto *ki_690 = buffer.data(ki + 690);
    const auto *ki_691 = buffer.data(ki + 691);
    const auto *ki_692 = buffer.data(ki + 692);
    const auto *ki_693 = buffer.data(ki + 693);
    const auto *ki_694 = buffer.data(ki + 694);
    const auto *ki_695 = buffer.data(ki + 695);
    const auto *ki_696 = buffer.data(ki + 696);
    const auto *ki_697 = buffer.data(ki + 697);
    const auto *ki_698 = buffer.data(ki + 698);
    const auto *ki_699 = buffer.data(ki + 699);
    const auto *ki_700 = buffer.data(ki + 700);
    const auto *ki_701 = buffer.data(ki + 701);
    const auto *ki_702 = buffer.data(ki + 702);

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, hi_207, hi_208, hi_209, hi_210, \
                         hi_211, ki_515, ki_516, ki_517, ki_518, \
                         ki_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -2.0 * hi_207[k]
                   + f_0 * ki_515[k];

        t_348[k] = -2.0 * hi_208[k]
                   + f_0 * ki_516[k];

        t_349[k] = -2.0 * hi_209[k]
                   + f_0 * ki_517[k];

        t_350[k] = -2.0 * hi_210[k]
                   + f_0 * ki_518[k];

        t_351[k] = -2.0 * hi_211[k]
                   + f_0 * ki_519[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, hi_212, hi_213, hi_214, hi_215, \
                         hi_216, ki_520, ki_521, ki_522, ki_523, \
                         ki_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -2.0 * hi_212[k]
                   + f_0 * ki_520[k];

        t_353[k] = -2.0 * hi_213[k]
                   + f_0 * ki_521[k];

        t_354[k] = -2.0 * hi_214[k]
                   + f_0 * ki_522[k];

        t_355[k] = -2.0 * hi_215[k]
                   + f_0 * ki_523[k];

        t_356[k] = -2.0 * hi_216[k]
                   + f_0 * ki_524[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, hi_217, hi_218, hi_219, hi_220, \
                         hi_221, ki_525, ki_526, ki_527, ki_528, \
                         ki_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -2.0 * hi_217[k]
                   + f_0 * ki_525[k];

        t_358[k] = -2.0 * hi_218[k]
                   + f_0 * ki_526[k];

        t_359[k] = -2.0 * hi_219[k]
                   + f_0 * ki_527[k];

        t_360[k] = -2.0 * hi_220[k]
                   + f_0 * ki_528[k];

        t_361[k] = -2.0 * hi_221[k]
                   + f_0 * ki_529[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, hi_222, hi_223, hi_224, hi_225, \
                         hi_226, ki_530, ki_531, ki_532, ki_533, \
                         ki_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -2.0 * hi_222[k]
                   + f_0 * ki_530[k];

        t_363[k] = -2.0 * hi_223[k]
                   + f_0 * ki_531[k];

        t_364[k] = -3.0 * hi_224[k]
                   + f_0 * ki_532[k];

        t_365[k] = -3.0 * hi_225[k]
                   + f_0 * ki_533[k];

        t_366[k] = -3.0 * hi_226[k]
                   + f_0 * ki_534[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, hi_227, hi_228, hi_229, hi_230, \
                         hi_231, ki_535, ki_536, ki_537, ki_538, \
                         ki_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -3.0 * hi_227[k]
                   + f_0 * ki_535[k];

        t_368[k] = -3.0 * hi_228[k]
                   + f_0 * ki_536[k];

        t_369[k] = -3.0 * hi_229[k]
                   + f_0 * ki_537[k];

        t_370[k] = -3.0 * hi_230[k]
                   + f_0 * ki_538[k];

        t_371[k] = -3.0 * hi_231[k]
                   + f_0 * ki_539[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, hi_232, hi_233, hi_234, hi_235, \
                         hi_236, ki_540, ki_541, ki_542, ki_543, \
                         ki_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -3.0 * hi_232[k]
                   + f_0 * ki_540[k];

        t_373[k] = -3.0 * hi_233[k]
                   + f_0 * ki_541[k];

        t_374[k] = -3.0 * hi_234[k]
                   + f_0 * ki_542[k];

        t_375[k] = -3.0 * hi_235[k]
                   + f_0 * ki_543[k];

        t_376[k] = -3.0 * hi_236[k]
                   + f_0 * ki_544[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, hi_237, hi_238, hi_239, hi_240, \
                         hi_241, ki_545, ki_546, ki_547, ki_548, \
                         ki_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -3.0 * hi_237[k]
                   + f_0 * ki_545[k];

        t_378[k] = -3.0 * hi_238[k]
                   + f_0 * ki_546[k];

        t_379[k] = -3.0 * hi_239[k]
                   + f_0 * ki_547[k];

        t_380[k] = -3.0 * hi_240[k]
                   + f_0 * ki_548[k];

        t_381[k] = -3.0 * hi_241[k]
                   + f_0 * ki_549[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, hi_242, hi_243, hi_244, hi_245, \
                         hi_246, ki_550, ki_551, ki_552, ki_553, \
                         ki_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -3.0 * hi_242[k]
                   + f_0 * ki_550[k];

        t_383[k] = -3.0 * hi_243[k]
                   + f_0 * ki_551[k];

        t_384[k] = -3.0 * hi_244[k]
                   + f_0 * ki_552[k];

        t_385[k] = -3.0 * hi_245[k]
                   + f_0 * ki_553[k];

        t_386[k] = -3.0 * hi_246[k]
                   + f_0 * ki_554[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, hi_247, hi_248, hi_249, hi_250, \
                         hi_251, ki_555, ki_556, ki_557, ki_558, \
                         ki_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -3.0 * hi_247[k]
                   + f_0 * ki_555[k];

        t_388[k] = -3.0 * hi_248[k]
                   + f_0 * ki_556[k];

        t_389[k] = -3.0 * hi_249[k]
                   + f_0 * ki_557[k];

        t_390[k] = -3.0 * hi_250[k]
                   + f_0 * ki_558[k];

        t_391[k] = -3.0 * hi_251[k]
                   + f_0 * ki_559[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, hi_252, hi_253, hi_254, hi_255, \
                         hi_256, ki_560, ki_561, ki_562, ki_563, \
                         ki_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -4.0 * hi_252[k]
                   + f_0 * ki_560[k];

        t_393[k] = -4.0 * hi_253[k]
                   + f_0 * ki_561[k];

        t_394[k] = -4.0 * hi_254[k]
                   + f_0 * ki_562[k];

        t_395[k] = -4.0 * hi_255[k]
                   + f_0 * ki_563[k];

        t_396[k] = -4.0 * hi_256[k]
                   + f_0 * ki_564[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, hi_257, hi_258, hi_259, hi_260, \
                         hi_261, ki_565, ki_566, ki_567, ki_568, \
                         ki_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -4.0 * hi_257[k]
                   + f_0 * ki_565[k];

        t_398[k] = -4.0 * hi_258[k]
                   + f_0 * ki_566[k];

        t_399[k] = -4.0 * hi_259[k]
                   + f_0 * ki_567[k];

        t_400[k] = -4.0 * hi_260[k]
                   + f_0 * ki_568[k];

        t_401[k] = -4.0 * hi_261[k]
                   + f_0 * ki_569[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, hi_262, hi_263, hi_264, hi_265, \
                         hi_266, ki_570, ki_571, ki_572, ki_573, \
                         ki_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -4.0 * hi_262[k]
                   + f_0 * ki_570[k];

        t_403[k] = -4.0 * hi_263[k]
                   + f_0 * ki_571[k];

        t_404[k] = -4.0 * hi_264[k]
                   + f_0 * ki_572[k];

        t_405[k] = -4.0 * hi_265[k]
                   + f_0 * ki_573[k];

        t_406[k] = -4.0 * hi_266[k]
                   + f_0 * ki_574[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, hi_267, hi_268, hi_269, hi_270, \
                         hi_271, ki_575, ki_576, ki_577, ki_578, \
                         ki_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -4.0 * hi_267[k]
                   + f_0 * ki_575[k];

        t_408[k] = -4.0 * hi_268[k]
                   + f_0 * ki_576[k];

        t_409[k] = -4.0 * hi_269[k]
                   + f_0 * ki_577[k];

        t_410[k] = -4.0 * hi_270[k]
                   + f_0 * ki_578[k];

        t_411[k] = -4.0 * hi_271[k]
                   + f_0 * ki_579[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, hi_272, hi_273, hi_274, hi_275, \
                         hi_276, ki_580, ki_581, ki_582, ki_583, \
                         ki_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -4.0 * hi_272[k]
                   + f_0 * ki_580[k];

        t_413[k] = -4.0 * hi_273[k]
                   + f_0 * ki_581[k];

        t_414[k] = -4.0 * hi_274[k]
                   + f_0 * ki_582[k];

        t_415[k] = -4.0 * hi_275[k]
                   + f_0 * ki_583[k];

        t_416[k] = -4.0 * hi_276[k]
                   + f_0 * ki_584[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, t_422, hi_277, hi_278, hi_279, \
                         ki_585, ki_586, ki_587, ki_616, ki_617, \
                         ki_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -4.0 * hi_277[k]
                   + f_0 * ki_585[k];

        t_418[k] = -4.0 * hi_278[k]
                   + f_0 * ki_586[k];

        t_419[k] = -4.0 * hi_279[k]
                   + f_0 * ki_587[k];

        t_420[k] = f_0 * ki_616[k];

        t_421[k] = f_0 * ki_617[k];

        t_422[k] = f_0 * ki_618[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, t_428, t_429, t_430, ki_619, \
                         ki_620, ki_621, ki_622, ki_623, ki_624, ki_625, \
                         ki_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_0 * ki_619[k];

        t_424[k] = f_0 * ki_620[k];

        t_425[k] = f_0 * ki_621[k];

        t_426[k] = f_0 * ki_622[k];

        t_427[k] = f_0 * ki_623[k];

        t_428[k] = f_0 * ki_624[k];

        t_429[k] = f_0 * ki_625[k];

        t_430[k] = f_0 * ki_626[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, t_436, t_437, t_438, ki_627, \
                         ki_628, ki_629, ki_630, ki_631, ki_632, ki_633, \
                         ki_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_0 * ki_627[k];

        t_432[k] = f_0 * ki_628[k];

        t_433[k] = f_0 * ki_629[k];

        t_434[k] = f_0 * ki_630[k];

        t_435[k] = f_0 * ki_631[k];

        t_436[k] = f_0 * ki_632[k];

        t_437[k] = f_0 * ki_633[k];

        t_438[k] = f_0 * ki_634[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, t_443, t_444, t_445, t_446, ki_635, \
                         ki_636, ki_637, ki_638, ki_639, ki_640, ki_641, \
                         ki_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_0 * ki_635[k];

        t_440[k] = f_0 * ki_636[k];

        t_441[k] = f_0 * ki_637[k];

        t_442[k] = f_0 * ki_638[k];

        t_443[k] = f_0 * ki_639[k];

        t_444[k] = f_0 * ki_640[k];

        t_445[k] = f_0 * ki_641[k];

        t_446[k] = f_0 * ki_642[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, t_451, hi_280, hi_281, hi_282, hi_283, \
                         ki_643, ki_644, ki_645, ki_646, ki_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_0 * ki_643[k];

        t_448[k] = -hi_280[k]
                   + f_0 * ki_644[k];

        t_449[k] = -hi_281[k]
                   + f_0 * ki_645[k];

        t_450[k] = -hi_282[k]
                   + f_0 * ki_646[k];

        t_451[k] = -hi_283[k]
                   + f_0 * ki_647[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, hi_284, hi_285, hi_286, hi_287, \
                         hi_288, ki_648, ki_649, ki_650, ki_651, \
                         ki_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = -hi_284[k]
                   + f_0 * ki_648[k];

        t_453[k] = -hi_285[k]
                   + f_0 * ki_649[k];

        t_454[k] = -hi_286[k]
                   + f_0 * ki_650[k];

        t_455[k] = -hi_287[k]
                   + f_0 * ki_651[k];

        t_456[k] = -hi_288[k]
                   + f_0 * ki_652[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, t_461, hi_289, hi_290, hi_291, hi_292, \
                         hi_293, ki_653, ki_654, ki_655, ki_656, \
                         ki_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = -hi_289[k]
                   + f_0 * ki_653[k];

        t_458[k] = -hi_290[k]
                   + f_0 * ki_654[k];

        t_459[k] = -hi_291[k]
                   + f_0 * ki_655[k];

        t_460[k] = -hi_292[k]
                   + f_0 * ki_656[k];

        t_461[k] = -hi_293[k]
                   + f_0 * ki_657[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, hi_294, hi_295, hi_296, hi_297, \
                         hi_298, ki_658, ki_659, ki_660, ki_661, \
                         ki_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = -hi_294[k]
                   + f_0 * ki_658[k];

        t_463[k] = -hi_295[k]
                   + f_0 * ki_659[k];

        t_464[k] = -hi_296[k]
                   + f_0 * ki_660[k];

        t_465[k] = -hi_297[k]
                   + f_0 * ki_661[k];

        t_466[k] = -hi_298[k]
                   + f_0 * ki_662[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, hi_299, hi_300, hi_301, hi_302, \
                         hi_303, ki_663, ki_664, ki_665, ki_666, \
                         ki_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = -hi_299[k]
                   + f_0 * ki_663[k];

        t_468[k] = -hi_300[k]
                   + f_0 * ki_664[k];

        t_469[k] = -hi_301[k]
                   + f_0 * ki_665[k];

        t_470[k] = -hi_302[k]
                   + f_0 * ki_666[k];

        t_471[k] = -hi_303[k]
                   + f_0 * ki_667[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, hi_304, hi_305, hi_306, hi_307, \
                         hi_308, ki_668, ki_669, ki_670, ki_671, \
                         ki_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = -hi_304[k]
                   + f_0 * ki_668[k];

        t_473[k] = -hi_305[k]
                   + f_0 * ki_669[k];

        t_474[k] = -hi_306[k]
                   + f_0 * ki_670[k];

        t_475[k] = -hi_307[k]
                   + f_0 * ki_671[k];

        t_476[k] = -2.0 * hi_308[k]
                   + f_0 * ki_672[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, hi_309, hi_310, hi_311, hi_312, \
                         hi_313, ki_673, ki_674, ki_675, ki_676, \
                         ki_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = -2.0 * hi_309[k]
                   + f_0 * ki_673[k];

        t_478[k] = -2.0 * hi_310[k]
                   + f_0 * ki_674[k];

        t_479[k] = -2.0 * hi_311[k]
                   + f_0 * ki_675[k];

        t_480[k] = -2.0 * hi_312[k]
                   + f_0 * ki_676[k];

        t_481[k] = -2.0 * hi_313[k]
                   + f_0 * ki_677[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, t_486, hi_314, hi_315, hi_316, hi_317, \
                         hi_318, ki_678, ki_679, ki_680, ki_681, \
                         ki_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = -2.0 * hi_314[k]
                   + f_0 * ki_678[k];

        t_483[k] = -2.0 * hi_315[k]
                   + f_0 * ki_679[k];

        t_484[k] = -2.0 * hi_316[k]
                   + f_0 * ki_680[k];

        t_485[k] = -2.0 * hi_317[k]
                   + f_0 * ki_681[k];

        t_486[k] = -2.0 * hi_318[k]
                   + f_0 * ki_682[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, hi_319, hi_320, hi_321, hi_322, \
                         hi_323, ki_683, ki_684, ki_685, ki_686, \
                         ki_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = -2.0 * hi_319[k]
                   + f_0 * ki_683[k];

        t_488[k] = -2.0 * hi_320[k]
                   + f_0 * ki_684[k];

        t_489[k] = -2.0 * hi_321[k]
                   + f_0 * ki_685[k];

        t_490[k] = -2.0 * hi_322[k]
                   + f_0 * ki_686[k];

        t_491[k] = -2.0 * hi_323[k]
                   + f_0 * ki_687[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, hi_324, hi_325, hi_326, hi_327, \
                         hi_328, ki_688, ki_689, ki_690, ki_691, \
                         ki_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = -2.0 * hi_324[k]
                   + f_0 * ki_688[k];

        t_493[k] = -2.0 * hi_325[k]
                   + f_0 * ki_689[k];

        t_494[k] = -2.0 * hi_326[k]
                   + f_0 * ki_690[k];

        t_495[k] = -2.0 * hi_327[k]
                   + f_0 * ki_691[k];

        t_496[k] = -2.0 * hi_328[k]
                   + f_0 * ki_692[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, t_500, t_501, hi_329, hi_330, hi_331, hi_332, \
                         hi_333, ki_693, ki_694, ki_695, ki_696, \
                         ki_697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = -2.0 * hi_329[k]
                   + f_0 * ki_693[k];

        t_498[k] = -2.0 * hi_330[k]
                   + f_0 * ki_694[k];

        t_499[k] = -2.0 * hi_331[k]
                   + f_0 * ki_695[k];

        t_500[k] = -2.0 * hi_332[k]
                   + f_0 * ki_696[k];

        t_501[k] = -2.0 * hi_333[k]
                   + f_0 * ki_697[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, hi_334, hi_335, hi_336, hi_337, \
                         hi_338, ki_698, ki_699, ki_700, ki_701, \
                         ki_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = -2.0 * hi_334[k]
                   + f_0 * ki_698[k];

        t_503[k] = -2.0 * hi_335[k]
                   + f_0 * ki_699[k];

        t_504[k] = -3.0 * hi_336[k]
                   + f_0 * ki_700[k];

        t_505[k] = -3.0 * hi_337[k]
                   + f_0 * ki_701[k];

        t_506[k] = -3.0 * hi_338[k]
                   + f_0 * ki_702[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
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

    const auto *hi_339 = buffer.data(hi + 339);
    const auto *hi_340 = buffer.data(hi + 340);
    const auto *hi_341 = buffer.data(hi + 341);
    const auto *hi_342 = buffer.data(hi + 342);
    const auto *hi_343 = buffer.data(hi + 343);
    const auto *hi_344 = buffer.data(hi + 344);
    const auto *hi_345 = buffer.data(hi + 345);
    const auto *hi_346 = buffer.data(hi + 346);
    const auto *hi_347 = buffer.data(hi + 347);
    const auto *hi_348 = buffer.data(hi + 348);
    const auto *hi_349 = buffer.data(hi + 349);
    const auto *hi_350 = buffer.data(hi + 350);
    const auto *hi_351 = buffer.data(hi + 351);
    const auto *hi_352 = buffer.data(hi + 352);
    const auto *hi_353 = buffer.data(hi + 353);
    const auto *hi_354 = buffer.data(hi + 354);
    const auto *hi_355 = buffer.data(hi + 355);
    const auto *hi_356 = buffer.data(hi + 356);
    const auto *hi_357 = buffer.data(hi + 357);
    const auto *hi_358 = buffer.data(hi + 358);
    const auto *hi_359 = buffer.data(hi + 359);
    const auto *hi_360 = buffer.data(hi + 360);
    const auto *hi_361 = buffer.data(hi + 361);
    const auto *hi_362 = buffer.data(hi + 362);
    const auto *hi_363 = buffer.data(hi + 363);
    const auto *hi_364 = buffer.data(hi + 364);
    const auto *hi_365 = buffer.data(hi + 365);
    const auto *hi_366 = buffer.data(hi + 366);
    const auto *hi_367 = buffer.data(hi + 367);
    const auto *hi_368 = buffer.data(hi + 368);
    const auto *hi_369 = buffer.data(hi + 369);
    const auto *hi_370 = buffer.data(hi + 370);
    const auto *hi_371 = buffer.data(hi + 371);
    const auto *hi_372 = buffer.data(hi + 372);
    const auto *hi_373 = buffer.data(hi + 373);
    const auto *hi_374 = buffer.data(hi + 374);
    const auto *hi_375 = buffer.data(hi + 375);
    const auto *hi_376 = buffer.data(hi + 376);
    const auto *hi_377 = buffer.data(hi + 377);
    const auto *hi_378 = buffer.data(hi + 378);
    const auto *hi_379 = buffer.data(hi + 379);
    const auto *hi_380 = buffer.data(hi + 380);
    const auto *hi_381 = buffer.data(hi + 381);
    const auto *hi_382 = buffer.data(hi + 382);
    const auto *hi_383 = buffer.data(hi + 383);
    const auto *hi_384 = buffer.data(hi + 384);
    const auto *hi_385 = buffer.data(hi + 385);
    const auto *hi_386 = buffer.data(hi + 386);
    const auto *hi_387 = buffer.data(hi + 387);
    const auto *hi_388 = buffer.data(hi + 388);
    const auto *hi_389 = buffer.data(hi + 389);
    const auto *hi_390 = buffer.data(hi + 390);
    const auto *hi_391 = buffer.data(hi + 391);
    const auto *hi_392 = buffer.data(hi + 392);
    const auto *hi_393 = buffer.data(hi + 393);
    const auto *hi_394 = buffer.data(hi + 394);
    const auto *hi_395 = buffer.data(hi + 395);
    const auto *hi_396 = buffer.data(hi + 396);
    const auto *hi_397 = buffer.data(hi + 397);
    const auto *hi_398 = buffer.data(hi + 398);
    const auto *hi_399 = buffer.data(hi + 399);
    const auto *hi_400 = buffer.data(hi + 400);
    const auto *hi_401 = buffer.data(hi + 401);
    const auto *hi_402 = buffer.data(hi + 402);
    const auto *hi_403 = buffer.data(hi + 403);
    const auto *hi_404 = buffer.data(hi + 404);
    const auto *hi_405 = buffer.data(hi + 405);
    const auto *hi_406 = buffer.data(hi + 406);
    const auto *hi_407 = buffer.data(hi + 407);
    const auto *hi_408 = buffer.data(hi + 408);
    const auto *hi_409 = buffer.data(hi + 409);
    const auto *hi_410 = buffer.data(hi + 410);
    const auto *hi_411 = buffer.data(hi + 411);
    const auto *hi_412 = buffer.data(hi + 412);
    const auto *hi_413 = buffer.data(hi + 413);
    const auto *hi_414 = buffer.data(hi + 414);
    const auto *hi_415 = buffer.data(hi + 415);
    const auto *hi_416 = buffer.data(hi + 416);
    const auto *hi_417 = buffer.data(hi + 417);
    const auto *hi_418 = buffer.data(hi + 418);
    const auto *hi_419 = buffer.data(hi + 419);
    const auto *hi_420 = buffer.data(hi + 420);
    const auto *hi_421 = buffer.data(hi + 421);
    const auto *hi_422 = buffer.data(hi + 422);
    const auto *hi_423 = buffer.data(hi + 423);
    const auto *hi_424 = buffer.data(hi + 424);
    const auto *hi_425 = buffer.data(hi + 425);
    const auto *hi_426 = buffer.data(hi + 426);
    const auto *hi_427 = buffer.data(hi + 427);
    const auto *hi_428 = buffer.data(hi + 428);
    const auto *hi_429 = buffer.data(hi + 429);
    const auto *hi_430 = buffer.data(hi + 430);
    const auto *hi_431 = buffer.data(hi + 431);
    const auto *hi_432 = buffer.data(hi + 432);
    const auto *hi_433 = buffer.data(hi + 433);
    const auto *hi_434 = buffer.data(hi + 434);
    const auto *hi_435 = buffer.data(hi + 435);
    const auto *hi_436 = buffer.data(hi + 436);
    const auto *hi_437 = buffer.data(hi + 437);
    const auto *hi_438 = buffer.data(hi + 438);
    const auto *hi_439 = buffer.data(hi + 439);
    const auto *hi_440 = buffer.data(hi + 440);
    const auto *hi_441 = buffer.data(hi + 441);
    const auto *hi_442 = buffer.data(hi + 442);
    const auto *hi_443 = buffer.data(hi + 443);
    const auto *hi_444 = buffer.data(hi + 444);
    const auto *hi_445 = buffer.data(hi + 445);
    const auto *hi_446 = buffer.data(hi + 446);
    const auto *hi_447 = buffer.data(hi + 447);
    const auto *hi_448 = buffer.data(hi + 448);
    const auto *hi_449 = buffer.data(hi + 449);
    const auto *hi_450 = buffer.data(hi + 450);
    const auto *hi_451 = buffer.data(hi + 451);
    const auto *hi_452 = buffer.data(hi + 452);
    const auto *hi_453 = buffer.data(hi + 453);
    const auto *hi_454 = buffer.data(hi + 454);
    const auto *hi_455 = buffer.data(hi + 455);
    const auto *hi_456 = buffer.data(hi + 456);
    const auto *hi_457 = buffer.data(hi + 457);
    const auto *hi_458 = buffer.data(hi + 458);
    const auto *hi_459 = buffer.data(hi + 459);
    const auto *hi_460 = buffer.data(hi + 460);
    const auto *hi_461 = buffer.data(hi + 461);
    const auto *hi_462 = buffer.data(hi + 462);
    const auto *hi_463 = buffer.data(hi + 463);
    const auto *hi_464 = buffer.data(hi + 464);
    const auto *hi_465 = buffer.data(hi + 465);
    const auto *hi_466 = buffer.data(hi + 466);
    const auto *hi_467 = buffer.data(hi + 467);
    const auto *hi_468 = buffer.data(hi + 468);
    const auto *hi_469 = buffer.data(hi + 469);
    const auto *hi_470 = buffer.data(hi + 470);

    const auto *ki_703 = buffer.data(ki + 703);
    const auto *ki_704 = buffer.data(ki + 704);
    const auto *ki_705 = buffer.data(ki + 705);
    const auto *ki_706 = buffer.data(ki + 706);
    const auto *ki_707 = buffer.data(ki + 707);
    const auto *ki_708 = buffer.data(ki + 708);
    const auto *ki_709 = buffer.data(ki + 709);
    const auto *ki_710 = buffer.data(ki + 710);
    const auto *ki_711 = buffer.data(ki + 711);
    const auto *ki_712 = buffer.data(ki + 712);
    const auto *ki_713 = buffer.data(ki + 713);
    const auto *ki_714 = buffer.data(ki + 714);
    const auto *ki_715 = buffer.data(ki + 715);
    const auto *ki_716 = buffer.data(ki + 716);
    const auto *ki_717 = buffer.data(ki + 717);
    const auto *ki_718 = buffer.data(ki + 718);
    const auto *ki_719 = buffer.data(ki + 719);
    const auto *ki_720 = buffer.data(ki + 720);
    const auto *ki_721 = buffer.data(ki + 721);
    const auto *ki_722 = buffer.data(ki + 722);
    const auto *ki_723 = buffer.data(ki + 723);
    const auto *ki_724 = buffer.data(ki + 724);
    const auto *ki_725 = buffer.data(ki + 725);
    const auto *ki_726 = buffer.data(ki + 726);
    const auto *ki_727 = buffer.data(ki + 727);
    const auto *ki_728 = buffer.data(ki + 728);
    const auto *ki_729 = buffer.data(ki + 729);
    const auto *ki_730 = buffer.data(ki + 730);
    const auto *ki_731 = buffer.data(ki + 731);
    const auto *ki_732 = buffer.data(ki + 732);
    const auto *ki_733 = buffer.data(ki + 733);
    const auto *ki_734 = buffer.data(ki + 734);
    const auto *ki_735 = buffer.data(ki + 735);
    const auto *ki_736 = buffer.data(ki + 736);
    const auto *ki_737 = buffer.data(ki + 737);
    const auto *ki_738 = buffer.data(ki + 738);
    const auto *ki_739 = buffer.data(ki + 739);
    const auto *ki_740 = buffer.data(ki + 740);
    const auto *ki_741 = buffer.data(ki + 741);
    const auto *ki_742 = buffer.data(ki + 742);
    const auto *ki_743 = buffer.data(ki + 743);
    const auto *ki_744 = buffer.data(ki + 744);
    const auto *ki_745 = buffer.data(ki + 745);
    const auto *ki_746 = buffer.data(ki + 746);
    const auto *ki_747 = buffer.data(ki + 747);
    const auto *ki_748 = buffer.data(ki + 748);
    const auto *ki_749 = buffer.data(ki + 749);
    const auto *ki_750 = buffer.data(ki + 750);
    const auto *ki_751 = buffer.data(ki + 751);
    const auto *ki_752 = buffer.data(ki + 752);
    const auto *ki_753 = buffer.data(ki + 753);
    const auto *ki_754 = buffer.data(ki + 754);
    const auto *ki_755 = buffer.data(ki + 755);
    const auto *ki_756 = buffer.data(ki + 756);
    const auto *ki_757 = buffer.data(ki + 757);
    const auto *ki_758 = buffer.data(ki + 758);
    const auto *ki_759 = buffer.data(ki + 759);
    const auto *ki_760 = buffer.data(ki + 760);
    const auto *ki_761 = buffer.data(ki + 761);
    const auto *ki_762 = buffer.data(ki + 762);
    const auto *ki_763 = buffer.data(ki + 763);
    const auto *ki_764 = buffer.data(ki + 764);
    const auto *ki_765 = buffer.data(ki + 765);
    const auto *ki_766 = buffer.data(ki + 766);
    const auto *ki_767 = buffer.data(ki + 767);
    const auto *ki_768 = buffer.data(ki + 768);
    const auto *ki_769 = buffer.data(ki + 769);
    const auto *ki_770 = buffer.data(ki + 770);
    const auto *ki_771 = buffer.data(ki + 771);
    const auto *ki_772 = buffer.data(ki + 772);
    const auto *ki_773 = buffer.data(ki + 773);
    const auto *ki_774 = buffer.data(ki + 774);
    const auto *ki_775 = buffer.data(ki + 775);
    const auto *ki_776 = buffer.data(ki + 776);
    const auto *ki_777 = buffer.data(ki + 777);
    const auto *ki_778 = buffer.data(ki + 778);
    const auto *ki_779 = buffer.data(ki + 779);
    const auto *ki_780 = buffer.data(ki + 780);
    const auto *ki_781 = buffer.data(ki + 781);
    const auto *ki_782 = buffer.data(ki + 782);
    const auto *ki_783 = buffer.data(ki + 783);
    const auto *ki_812 = buffer.data(ki + 812);
    const auto *ki_813 = buffer.data(ki + 813);
    const auto *ki_814 = buffer.data(ki + 814);
    const auto *ki_815 = buffer.data(ki + 815);
    const auto *ki_816 = buffer.data(ki + 816);
    const auto *ki_817 = buffer.data(ki + 817);
    const auto *ki_818 = buffer.data(ki + 818);
    const auto *ki_819 = buffer.data(ki + 819);
    const auto *ki_820 = buffer.data(ki + 820);
    const auto *ki_821 = buffer.data(ki + 821);
    const auto *ki_822 = buffer.data(ki + 822);
    const auto *ki_823 = buffer.data(ki + 823);
    const auto *ki_824 = buffer.data(ki + 824);
    const auto *ki_825 = buffer.data(ki + 825);
    const auto *ki_826 = buffer.data(ki + 826);
    const auto *ki_827 = buffer.data(ki + 827);
    const auto *ki_828 = buffer.data(ki + 828);
    const auto *ki_829 = buffer.data(ki + 829);
    const auto *ki_830 = buffer.data(ki + 830);
    const auto *ki_831 = buffer.data(ki + 831);
    const auto *ki_832 = buffer.data(ki + 832);
    const auto *ki_833 = buffer.data(ki + 833);
    const auto *ki_834 = buffer.data(ki + 834);
    const auto *ki_835 = buffer.data(ki + 835);
    const auto *ki_836 = buffer.data(ki + 836);
    const auto *ki_837 = buffer.data(ki + 837);
    const auto *ki_838 = buffer.data(ki + 838);
    const auto *ki_839 = buffer.data(ki + 839);
    const auto *ki_840 = buffer.data(ki + 840);
    const auto *ki_841 = buffer.data(ki + 841);
    const auto *ki_842 = buffer.data(ki + 842);
    const auto *ki_843 = buffer.data(ki + 843);
    const auto *ki_844 = buffer.data(ki + 844);
    const auto *ki_845 = buffer.data(ki + 845);
    const auto *ki_846 = buffer.data(ki + 846);
    const auto *ki_847 = buffer.data(ki + 847);
    const auto *ki_848 = buffer.data(ki + 848);
    const auto *ki_849 = buffer.data(ki + 849);
    const auto *ki_850 = buffer.data(ki + 850);
    const auto *ki_851 = buffer.data(ki + 851);
    const auto *ki_852 = buffer.data(ki + 852);
    const auto *ki_853 = buffer.data(ki + 853);
    const auto *ki_854 = buffer.data(ki + 854);
    const auto *ki_855 = buffer.data(ki + 855);
    const auto *ki_856 = buffer.data(ki + 856);
    const auto *ki_857 = buffer.data(ki + 857);
    const auto *ki_858 = buffer.data(ki + 858);
    const auto *ki_859 = buffer.data(ki + 859);
    const auto *ki_860 = buffer.data(ki + 860);
    const auto *ki_861 = buffer.data(ki + 861);
    const auto *ki_862 = buffer.data(ki + 862);
    const auto *ki_863 = buffer.data(ki + 863);
    const auto *ki_864 = buffer.data(ki + 864);
    const auto *ki_865 = buffer.data(ki + 865);
    const auto *ki_866 = buffer.data(ki + 866);
    const auto *ki_867 = buffer.data(ki + 867);
    const auto *ki_868 = buffer.data(ki + 868);
    const auto *ki_869 = buffer.data(ki + 869);
    const auto *ki_870 = buffer.data(ki + 870);
    const auto *ki_871 = buffer.data(ki + 871);
    const auto *ki_872 = buffer.data(ki + 872);
    const auto *ki_873 = buffer.data(ki + 873);
    const auto *ki_874 = buffer.data(ki + 874);
    const auto *ki_875 = buffer.data(ki + 875);
    const auto *ki_876 = buffer.data(ki + 876);
    const auto *ki_877 = buffer.data(ki + 877);
    const auto *ki_878 = buffer.data(ki + 878);
    const auto *ki_879 = buffer.data(ki + 879);
    const auto *ki_880 = buffer.data(ki + 880);
    const auto *ki_881 = buffer.data(ki + 881);
    const auto *ki_882 = buffer.data(ki + 882);
    const auto *ki_883 = buffer.data(ki + 883);
    const auto *ki_884 = buffer.data(ki + 884);
    const auto *ki_885 = buffer.data(ki + 885);
    const auto *ki_886 = buffer.data(ki + 886);
    const auto *ki_887 = buffer.data(ki + 887);
    const auto *ki_888 = buffer.data(ki + 888);
    const auto *ki_889 = buffer.data(ki + 889);
    const auto *ki_890 = buffer.data(ki + 890);

#pragma omp simd aligned(t_507, t_508, t_509, t_510, t_511, hi_339, hi_340, hi_341, hi_342, \
                         hi_343, ki_703, ki_704, ki_705, ki_706, \
                         ki_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = -3.0 * hi_339[k]
                   + f_0 * ki_703[k];

        t_508[k] = -3.0 * hi_340[k]
                   + f_0 * ki_704[k];

        t_509[k] = -3.0 * hi_341[k]
                   + f_0 * ki_705[k];

        t_510[k] = -3.0 * hi_342[k]
                   + f_0 * ki_706[k];

        t_511[k] = -3.0 * hi_343[k]
                   + f_0 * ki_707[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, hi_344, hi_345, hi_346, hi_347, \
                         hi_348, ki_708, ki_709, ki_710, ki_711, \
                         ki_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = -3.0 * hi_344[k]
                   + f_0 * ki_708[k];

        t_513[k] = -3.0 * hi_345[k]
                   + f_0 * ki_709[k];

        t_514[k] = -3.0 * hi_346[k]
                   + f_0 * ki_710[k];

        t_515[k] = -3.0 * hi_347[k]
                   + f_0 * ki_711[k];

        t_516[k] = -3.0 * hi_348[k]
                   + f_0 * ki_712[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, hi_349, hi_350, hi_351, hi_352, \
                         hi_353, ki_713, ki_714, ki_715, ki_716, \
                         ki_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = -3.0 * hi_349[k]
                   + f_0 * ki_713[k];

        t_518[k] = -3.0 * hi_350[k]
                   + f_0 * ki_714[k];

        t_519[k] = -3.0 * hi_351[k]
                   + f_0 * ki_715[k];

        t_520[k] = -3.0 * hi_352[k]
                   + f_0 * ki_716[k];

        t_521[k] = -3.0 * hi_353[k]
                   + f_0 * ki_717[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, hi_354, hi_355, hi_356, hi_357, \
                         hi_358, ki_718, ki_719, ki_720, ki_721, \
                         ki_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = -3.0 * hi_354[k]
                   + f_0 * ki_718[k];

        t_523[k] = -3.0 * hi_355[k]
                   + f_0 * ki_719[k];

        t_524[k] = -3.0 * hi_356[k]
                   + f_0 * ki_720[k];

        t_525[k] = -3.0 * hi_357[k]
                   + f_0 * ki_721[k];

        t_526[k] = -3.0 * hi_358[k]
                   + f_0 * ki_722[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, hi_359, hi_360, hi_361, hi_362, \
                         hi_363, ki_723, ki_724, ki_725, ki_726, \
                         ki_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = -3.0 * hi_359[k]
                   + f_0 * ki_723[k];

        t_528[k] = -3.0 * hi_360[k]
                   + f_0 * ki_724[k];

        t_529[k] = -3.0 * hi_361[k]
                   + f_0 * ki_725[k];

        t_530[k] = -3.0 * hi_362[k]
                   + f_0 * ki_726[k];

        t_531[k] = -3.0 * hi_363[k]
                   + f_0 * ki_727[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, hi_364, hi_365, hi_366, hi_367, \
                         hi_368, ki_728, ki_729, ki_730, ki_731, \
                         ki_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = -4.0 * hi_364[k]
                   + f_0 * ki_728[k];

        t_533[k] = -4.0 * hi_365[k]
                   + f_0 * ki_729[k];

        t_534[k] = -4.0 * hi_366[k]
                   + f_0 * ki_730[k];

        t_535[k] = -4.0 * hi_367[k]
                   + f_0 * ki_731[k];

        t_536[k] = -4.0 * hi_368[k]
                   + f_0 * ki_732[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, t_541, hi_369, hi_370, hi_371, hi_372, \
                         hi_373, ki_733, ki_734, ki_735, ki_736, \
                         ki_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = -4.0 * hi_369[k]
                   + f_0 * ki_733[k];

        t_538[k] = -4.0 * hi_370[k]
                   + f_0 * ki_734[k];

        t_539[k] = -4.0 * hi_371[k]
                   + f_0 * ki_735[k];

        t_540[k] = -4.0 * hi_372[k]
                   + f_0 * ki_736[k];

        t_541[k] = -4.0 * hi_373[k]
                   + f_0 * ki_737[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, hi_374, hi_375, hi_376, hi_377, \
                         hi_378, ki_738, ki_739, ki_740, ki_741, \
                         ki_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = -4.0 * hi_374[k]
                   + f_0 * ki_738[k];

        t_543[k] = -4.0 * hi_375[k]
                   + f_0 * ki_739[k];

        t_544[k] = -4.0 * hi_376[k]
                   + f_0 * ki_740[k];

        t_545[k] = -4.0 * hi_377[k]
                   + f_0 * ki_741[k];

        t_546[k] = -4.0 * hi_378[k]
                   + f_0 * ki_742[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, hi_379, hi_380, hi_381, hi_382, \
                         hi_383, ki_743, ki_744, ki_745, ki_746, \
                         ki_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = -4.0 * hi_379[k]
                   + f_0 * ki_743[k];

        t_548[k] = -4.0 * hi_380[k]
                   + f_0 * ki_744[k];

        t_549[k] = -4.0 * hi_381[k]
                   + f_0 * ki_745[k];

        t_550[k] = -4.0 * hi_382[k]
                   + f_0 * ki_746[k];

        t_551[k] = -4.0 * hi_383[k]
                   + f_0 * ki_747[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, t_556, hi_384, hi_385, hi_386, hi_387, \
                         hi_388, ki_748, ki_749, ki_750, ki_751, \
                         ki_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = -4.0 * hi_384[k]
                   + f_0 * ki_748[k];

        t_553[k] = -4.0 * hi_385[k]
                   + f_0 * ki_749[k];

        t_554[k] = -4.0 * hi_386[k]
                   + f_0 * ki_750[k];

        t_555[k] = -4.0 * hi_387[k]
                   + f_0 * ki_751[k];

        t_556[k] = -4.0 * hi_388[k]
                   + f_0 * ki_752[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, hi_389, hi_390, hi_391, hi_392, \
                         hi_393, ki_753, ki_754, ki_755, ki_756, \
                         ki_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = -4.0 * hi_389[k]
                   + f_0 * ki_753[k];

        t_558[k] = -4.0 * hi_390[k]
                   + f_0 * ki_754[k];

        t_559[k] = -4.0 * hi_391[k]
                   + f_0 * ki_755[k];

        t_560[k] = -5.0 * hi_392[k]
                   + f_0 * ki_756[k];

        t_561[k] = -5.0 * hi_393[k]
                   + f_0 * ki_757[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, hi_394, hi_395, hi_396, hi_397, \
                         hi_398, ki_758, ki_759, ki_760, ki_761, \
                         ki_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = -5.0 * hi_394[k]
                   + f_0 * ki_758[k];

        t_563[k] = -5.0 * hi_395[k]
                   + f_0 * ki_759[k];

        t_564[k] = -5.0 * hi_396[k]
                   + f_0 * ki_760[k];

        t_565[k] = -5.0 * hi_397[k]
                   + f_0 * ki_761[k];

        t_566[k] = -5.0 * hi_398[k]
                   + f_0 * ki_762[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, hi_399, hi_400, hi_401, hi_402, \
                         hi_403, ki_763, ki_764, ki_765, ki_766, \
                         ki_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -5.0 * hi_399[k]
                   + f_0 * ki_763[k];

        t_568[k] = -5.0 * hi_400[k]
                   + f_0 * ki_764[k];

        t_569[k] = -5.0 * hi_401[k]
                   + f_0 * ki_765[k];

        t_570[k] = -5.0 * hi_402[k]
                   + f_0 * ki_766[k];

        t_571[k] = -5.0 * hi_403[k]
                   + f_0 * ki_767[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, hi_404, hi_405, hi_406, hi_407, \
                         hi_408, ki_768, ki_769, ki_770, ki_771, \
                         ki_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = -5.0 * hi_404[k]
                   + f_0 * ki_768[k];

        t_573[k] = -5.0 * hi_405[k]
                   + f_0 * ki_769[k];

        t_574[k] = -5.0 * hi_406[k]
                   + f_0 * ki_770[k];

        t_575[k] = -5.0 * hi_407[k]
                   + f_0 * ki_771[k];

        t_576[k] = -5.0 * hi_408[k]
                   + f_0 * ki_772[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, hi_409, hi_410, hi_411, hi_412, \
                         hi_413, ki_773, ki_774, ki_775, ki_776, \
                         ki_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = -5.0 * hi_409[k]
                   + f_0 * ki_773[k];

        t_578[k] = -5.0 * hi_410[k]
                   + f_0 * ki_774[k];

        t_579[k] = -5.0 * hi_411[k]
                   + f_0 * ki_775[k];

        t_580[k] = -5.0 * hi_412[k]
                   + f_0 * ki_776[k];

        t_581[k] = -5.0 * hi_413[k]
                   + f_0 * ki_777[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, hi_414, hi_415, hi_416, hi_417, \
                         hi_418, ki_778, ki_779, ki_780, ki_781, \
                         ki_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = -5.0 * hi_414[k]
                   + f_0 * ki_778[k];

        t_583[k] = -5.0 * hi_415[k]
                   + f_0 * ki_779[k];

        t_584[k] = -5.0 * hi_416[k]
                   + f_0 * ki_780[k];

        t_585[k] = -5.0 * hi_417[k]
                   + f_0 * ki_781[k];

        t_586[k] = -5.0 * hi_418[k]
                   + f_0 * ki_782[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, t_591, t_592, t_593, hi_419, ki_783, \
                         ki_812, ki_813, ki_814, ki_815, ki_816, \
                         ki_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = -5.0 * hi_419[k]
                   + f_0 * ki_783[k];

        t_588[k] = f_0 * ki_812[k];

        t_589[k] = f_0 * ki_813[k];

        t_590[k] = f_0 * ki_814[k];

        t_591[k] = f_0 * ki_815[k];

        t_592[k] = f_0 * ki_816[k];

        t_593[k] = f_0 * ki_817[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, t_598, t_599, t_600, t_601, ki_818, \
                         ki_819, ki_820, ki_821, ki_822, ki_823, ki_824, \
                         ki_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_0 * ki_818[k];

        t_595[k] = f_0 * ki_819[k];

        t_596[k] = f_0 * ki_820[k];

        t_597[k] = f_0 * ki_821[k];

        t_598[k] = f_0 * ki_822[k];

        t_599[k] = f_0 * ki_823[k];

        t_600[k] = f_0 * ki_824[k];

        t_601[k] = f_0 * ki_825[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, t_607, t_608, t_609, ki_826, \
                         ki_827, ki_828, ki_829, ki_830, ki_831, ki_832, \
                         ki_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_0 * ki_826[k];

        t_603[k] = f_0 * ki_827[k];

        t_604[k] = f_0 * ki_828[k];

        t_605[k] = f_0 * ki_829[k];

        t_606[k] = f_0 * ki_830[k];

        t_607[k] = f_0 * ki_831[k];

        t_608[k] = f_0 * ki_832[k];

        t_609[k] = f_0 * ki_833[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, t_615, t_616, hi_420, ki_834, \
                         ki_835, ki_836, ki_837, ki_838, ki_839, \
                         ki_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_0 * ki_834[k];

        t_611[k] = f_0 * ki_835[k];

        t_612[k] = f_0 * ki_836[k];

        t_613[k] = f_0 * ki_837[k];

        t_614[k] = f_0 * ki_838[k];

        t_615[k] = f_0 * ki_839[k];

        t_616[k] = -hi_420[k]
                   + f_0 * ki_840[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, hi_421, hi_422, hi_423, hi_424, \
                         hi_425, ki_841, ki_842, ki_843, ki_844, \
                         ki_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = -hi_421[k]
                   + f_0 * ki_841[k];

        t_618[k] = -hi_422[k]
                   + f_0 * ki_842[k];

        t_619[k] = -hi_423[k]
                   + f_0 * ki_843[k];

        t_620[k] = -hi_424[k]
                   + f_0 * ki_844[k];

        t_621[k] = -hi_425[k]
                   + f_0 * ki_845[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, hi_426, hi_427, hi_428, hi_429, \
                         hi_430, ki_846, ki_847, ki_848, ki_849, \
                         ki_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = -hi_426[k]
                   + f_0 * ki_846[k];

        t_623[k] = -hi_427[k]
                   + f_0 * ki_847[k];

        t_624[k] = -hi_428[k]
                   + f_0 * ki_848[k];

        t_625[k] = -hi_429[k]
                   + f_0 * ki_849[k];

        t_626[k] = -hi_430[k]
                   + f_0 * ki_850[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, t_631, hi_431, hi_432, hi_433, hi_434, \
                         hi_435, ki_851, ki_852, ki_853, ki_854, \
                         ki_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = -hi_431[k]
                   + f_0 * ki_851[k];

        t_628[k] = -hi_432[k]
                   + f_0 * ki_852[k];

        t_629[k] = -hi_433[k]
                   + f_0 * ki_853[k];

        t_630[k] = -hi_434[k]
                   + f_0 * ki_854[k];

        t_631[k] = -hi_435[k]
                   + f_0 * ki_855[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, t_636, hi_436, hi_437, hi_438, hi_439, \
                         hi_440, ki_856, ki_857, ki_858, ki_859, \
                         ki_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = -hi_436[k]
                   + f_0 * ki_856[k];

        t_633[k] = -hi_437[k]
                   + f_0 * ki_857[k];

        t_634[k] = -hi_438[k]
                   + f_0 * ki_858[k];

        t_635[k] = -hi_439[k]
                   + f_0 * ki_859[k];

        t_636[k] = -hi_440[k]
                   + f_0 * ki_860[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, t_641, hi_441, hi_442, hi_443, hi_444, \
                         hi_445, ki_861, ki_862, ki_863, ki_864, \
                         ki_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = -hi_441[k]
                   + f_0 * ki_861[k];

        t_638[k] = -hi_442[k]
                   + f_0 * ki_862[k];

        t_639[k] = -hi_443[k]
                   + f_0 * ki_863[k];

        t_640[k] = -hi_444[k]
                   + f_0 * ki_864[k];

        t_641[k] = -hi_445[k]
                   + f_0 * ki_865[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, t_646, hi_446, hi_447, hi_448, hi_449, \
                         hi_450, ki_866, ki_867, ki_868, ki_869, \
                         ki_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = -hi_446[k]
                   + f_0 * ki_866[k];

        t_643[k] = -hi_447[k]
                   + f_0 * ki_867[k];

        t_644[k] = -2.0 * hi_448[k]
                   + f_0 * ki_868[k];

        t_645[k] = -2.0 * hi_449[k]
                   + f_0 * ki_869[k];

        t_646[k] = -2.0 * hi_450[k]
                   + f_0 * ki_870[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, t_651, hi_451, hi_452, hi_453, hi_454, \
                         hi_455, ki_871, ki_872, ki_873, ki_874, \
                         ki_875 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = -2.0 * hi_451[k]
                   + f_0 * ki_871[k];

        t_648[k] = -2.0 * hi_452[k]
                   + f_0 * ki_872[k];

        t_649[k] = -2.0 * hi_453[k]
                   + f_0 * ki_873[k];

        t_650[k] = -2.0 * hi_454[k]
                   + f_0 * ki_874[k];

        t_651[k] = -2.0 * hi_455[k]
                   + f_0 * ki_875[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, t_656, hi_456, hi_457, hi_458, hi_459, \
                         hi_460, ki_876, ki_877, ki_878, ki_879, \
                         ki_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = -2.0 * hi_456[k]
                   + f_0 * ki_876[k];

        t_653[k] = -2.0 * hi_457[k]
                   + f_0 * ki_877[k];

        t_654[k] = -2.0 * hi_458[k]
                   + f_0 * ki_878[k];

        t_655[k] = -2.0 * hi_459[k]
                   + f_0 * ki_879[k];

        t_656[k] = -2.0 * hi_460[k]
                   + f_0 * ki_880[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, hi_461, hi_462, hi_463, hi_464, \
                         hi_465, ki_881, ki_882, ki_883, ki_884, \
                         ki_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = -2.0 * hi_461[k]
                   + f_0 * ki_881[k];

        t_658[k] = -2.0 * hi_462[k]
                   + f_0 * ki_882[k];

        t_659[k] = -2.0 * hi_463[k]
                   + f_0 * ki_883[k];

        t_660[k] = -2.0 * hi_464[k]
                   + f_0 * ki_884[k];

        t_661[k] = -2.0 * hi_465[k]
                   + f_0 * ki_885[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, t_666, hi_466, hi_467, hi_468, hi_469, \
                         hi_470, ki_886, ki_887, ki_888, ki_889, \
                         ki_890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = -2.0 * hi_466[k]
                   + f_0 * ki_886[k];

        t_663[k] = -2.0 * hi_467[k]
                   + f_0 * ki_887[k];

        t_664[k] = -2.0 * hi_468[k]
                   + f_0 * ki_888[k];

        t_665[k] = -2.0 * hi_469[k]
                   + f_0 * ki_889[k];

        t_666[k] = -2.0 * hi_470[k]
                   + f_0 * ki_890[k];
    }
}

static auto
compute_prim_geom_10_ii_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ki,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hi_471 = buffer.data(hi + 471);
    const auto *hi_472 = buffer.data(hi + 472);
    const auto *hi_473 = buffer.data(hi + 473);
    const auto *hi_474 = buffer.data(hi + 474);
    const auto *hi_475 = buffer.data(hi + 475);
    const auto *hi_476 = buffer.data(hi + 476);
    const auto *hi_477 = buffer.data(hi + 477);
    const auto *hi_478 = buffer.data(hi + 478);
    const auto *hi_479 = buffer.data(hi + 479);
    const auto *hi_480 = buffer.data(hi + 480);
    const auto *hi_481 = buffer.data(hi + 481);
    const auto *hi_482 = buffer.data(hi + 482);
    const auto *hi_483 = buffer.data(hi + 483);
    const auto *hi_484 = buffer.data(hi + 484);
    const auto *hi_485 = buffer.data(hi + 485);
    const auto *hi_486 = buffer.data(hi + 486);
    const auto *hi_487 = buffer.data(hi + 487);
    const auto *hi_488 = buffer.data(hi + 488);
    const auto *hi_489 = buffer.data(hi + 489);
    const auto *hi_490 = buffer.data(hi + 490);
    const auto *hi_491 = buffer.data(hi + 491);
    const auto *hi_492 = buffer.data(hi + 492);
    const auto *hi_493 = buffer.data(hi + 493);
    const auto *hi_494 = buffer.data(hi + 494);
    const auto *hi_495 = buffer.data(hi + 495);
    const auto *hi_496 = buffer.data(hi + 496);
    const auto *hi_497 = buffer.data(hi + 497);
    const auto *hi_498 = buffer.data(hi + 498);
    const auto *hi_499 = buffer.data(hi + 499);
    const auto *hi_500 = buffer.data(hi + 500);
    const auto *hi_501 = buffer.data(hi + 501);
    const auto *hi_502 = buffer.data(hi + 502);
    const auto *hi_503 = buffer.data(hi + 503);
    const auto *hi_504 = buffer.data(hi + 504);
    const auto *hi_505 = buffer.data(hi + 505);
    const auto *hi_506 = buffer.data(hi + 506);
    const auto *hi_507 = buffer.data(hi + 507);
    const auto *hi_508 = buffer.data(hi + 508);
    const auto *hi_509 = buffer.data(hi + 509);
    const auto *hi_510 = buffer.data(hi + 510);
    const auto *hi_511 = buffer.data(hi + 511);
    const auto *hi_512 = buffer.data(hi + 512);
    const auto *hi_513 = buffer.data(hi + 513);
    const auto *hi_514 = buffer.data(hi + 514);
    const auto *hi_515 = buffer.data(hi + 515);
    const auto *hi_516 = buffer.data(hi + 516);
    const auto *hi_517 = buffer.data(hi + 517);
    const auto *hi_518 = buffer.data(hi + 518);
    const auto *hi_519 = buffer.data(hi + 519);
    const auto *hi_520 = buffer.data(hi + 520);
    const auto *hi_521 = buffer.data(hi + 521);
    const auto *hi_522 = buffer.data(hi + 522);
    const auto *hi_523 = buffer.data(hi + 523);
    const auto *hi_524 = buffer.data(hi + 524);
    const auto *hi_525 = buffer.data(hi + 525);
    const auto *hi_526 = buffer.data(hi + 526);
    const auto *hi_527 = buffer.data(hi + 527);
    const auto *hi_528 = buffer.data(hi + 528);
    const auto *hi_529 = buffer.data(hi + 529);
    const auto *hi_530 = buffer.data(hi + 530);
    const auto *hi_531 = buffer.data(hi + 531);
    const auto *hi_532 = buffer.data(hi + 532);
    const auto *hi_533 = buffer.data(hi + 533);
    const auto *hi_534 = buffer.data(hi + 534);
    const auto *hi_535 = buffer.data(hi + 535);
    const auto *hi_536 = buffer.data(hi + 536);
    const auto *hi_537 = buffer.data(hi + 537);
    const auto *hi_538 = buffer.data(hi + 538);
    const auto *hi_539 = buffer.data(hi + 539);
    const auto *hi_540 = buffer.data(hi + 540);
    const auto *hi_541 = buffer.data(hi + 541);
    const auto *hi_542 = buffer.data(hi + 542);
    const auto *hi_543 = buffer.data(hi + 543);
    const auto *hi_544 = buffer.data(hi + 544);
    const auto *hi_545 = buffer.data(hi + 545);
    const auto *hi_546 = buffer.data(hi + 546);
    const auto *hi_547 = buffer.data(hi + 547);
    const auto *hi_548 = buffer.data(hi + 548);
    const auto *hi_549 = buffer.data(hi + 549);
    const auto *hi_550 = buffer.data(hi + 550);
    const auto *hi_551 = buffer.data(hi + 551);
    const auto *hi_552 = buffer.data(hi + 552);
    const auto *hi_553 = buffer.data(hi + 553);
    const auto *hi_554 = buffer.data(hi + 554);
    const auto *hi_555 = buffer.data(hi + 555);
    const auto *hi_556 = buffer.data(hi + 556);
    const auto *hi_557 = buffer.data(hi + 557);
    const auto *hi_558 = buffer.data(hi + 558);
    const auto *hi_559 = buffer.data(hi + 559);
    const auto *hi_560 = buffer.data(hi + 560);
    const auto *hi_561 = buffer.data(hi + 561);
    const auto *hi_562 = buffer.data(hi + 562);
    const auto *hi_563 = buffer.data(hi + 563);
    const auto *hi_564 = buffer.data(hi + 564);
    const auto *hi_565 = buffer.data(hi + 565);
    const auto *hi_566 = buffer.data(hi + 566);
    const auto *hi_567 = buffer.data(hi + 567);
    const auto *hi_568 = buffer.data(hi + 568);
    const auto *hi_569 = buffer.data(hi + 569);
    const auto *hi_570 = buffer.data(hi + 570);
    const auto *hi_571 = buffer.data(hi + 571);
    const auto *hi_572 = buffer.data(hi + 572);
    const auto *hi_573 = buffer.data(hi + 573);
    const auto *hi_574 = buffer.data(hi + 574);
    const auto *hi_575 = buffer.data(hi + 575);
    const auto *hi_576 = buffer.data(hi + 576);
    const auto *hi_577 = buffer.data(hi + 577);
    const auto *hi_578 = buffer.data(hi + 578);
    const auto *hi_579 = buffer.data(hi + 579);
    const auto *hi_580 = buffer.data(hi + 580);
    const auto *hi_581 = buffer.data(hi + 581);
    const auto *hi_582 = buffer.data(hi + 582);
    const auto *hi_583 = buffer.data(hi + 583);
    const auto *hi_584 = buffer.data(hi + 584);
    const auto *hi_585 = buffer.data(hi + 585);
    const auto *hi_586 = buffer.data(hi + 586);
    const auto *hi_587 = buffer.data(hi + 587);

    const auto *ki_891 = buffer.data(ki + 891);
    const auto *ki_892 = buffer.data(ki + 892);
    const auto *ki_893 = buffer.data(ki + 893);
    const auto *ki_894 = buffer.data(ki + 894);
    const auto *ki_895 = buffer.data(ki + 895);
    const auto *ki_896 = buffer.data(ki + 896);
    const auto *ki_897 = buffer.data(ki + 897);
    const auto *ki_898 = buffer.data(ki + 898);
    const auto *ki_899 = buffer.data(ki + 899);
    const auto *ki_900 = buffer.data(ki + 900);
    const auto *ki_901 = buffer.data(ki + 901);
    const auto *ki_902 = buffer.data(ki + 902);
    const auto *ki_903 = buffer.data(ki + 903);
    const auto *ki_904 = buffer.data(ki + 904);
    const auto *ki_905 = buffer.data(ki + 905);
    const auto *ki_906 = buffer.data(ki + 906);
    const auto *ki_907 = buffer.data(ki + 907);
    const auto *ki_908 = buffer.data(ki + 908);
    const auto *ki_909 = buffer.data(ki + 909);
    const auto *ki_910 = buffer.data(ki + 910);
    const auto *ki_911 = buffer.data(ki + 911);
    const auto *ki_912 = buffer.data(ki + 912);
    const auto *ki_913 = buffer.data(ki + 913);
    const auto *ki_914 = buffer.data(ki + 914);
    const auto *ki_915 = buffer.data(ki + 915);
    const auto *ki_916 = buffer.data(ki + 916);
    const auto *ki_917 = buffer.data(ki + 917);
    const auto *ki_918 = buffer.data(ki + 918);
    const auto *ki_919 = buffer.data(ki + 919);
    const auto *ki_920 = buffer.data(ki + 920);
    const auto *ki_921 = buffer.data(ki + 921);
    const auto *ki_922 = buffer.data(ki + 922);
    const auto *ki_923 = buffer.data(ki + 923);
    const auto *ki_924 = buffer.data(ki + 924);
    const auto *ki_925 = buffer.data(ki + 925);
    const auto *ki_926 = buffer.data(ki + 926);
    const auto *ki_927 = buffer.data(ki + 927);
    const auto *ki_928 = buffer.data(ki + 928);
    const auto *ki_929 = buffer.data(ki + 929);
    const auto *ki_930 = buffer.data(ki + 930);
    const auto *ki_931 = buffer.data(ki + 931);
    const auto *ki_932 = buffer.data(ki + 932);
    const auto *ki_933 = buffer.data(ki + 933);
    const auto *ki_934 = buffer.data(ki + 934);
    const auto *ki_935 = buffer.data(ki + 935);
    const auto *ki_936 = buffer.data(ki + 936);
    const auto *ki_937 = buffer.data(ki + 937);
    const auto *ki_938 = buffer.data(ki + 938);
    const auto *ki_939 = buffer.data(ki + 939);
    const auto *ki_940 = buffer.data(ki + 940);
    const auto *ki_941 = buffer.data(ki + 941);
    const auto *ki_942 = buffer.data(ki + 942);
    const auto *ki_943 = buffer.data(ki + 943);
    const auto *ki_944 = buffer.data(ki + 944);
    const auto *ki_945 = buffer.data(ki + 945);
    const auto *ki_946 = buffer.data(ki + 946);
    const auto *ki_947 = buffer.data(ki + 947);
    const auto *ki_948 = buffer.data(ki + 948);
    const auto *ki_949 = buffer.data(ki + 949);
    const auto *ki_950 = buffer.data(ki + 950);
    const auto *ki_951 = buffer.data(ki + 951);
    const auto *ki_952 = buffer.data(ki + 952);
    const auto *ki_953 = buffer.data(ki + 953);
    const auto *ki_954 = buffer.data(ki + 954);
    const auto *ki_955 = buffer.data(ki + 955);
    const auto *ki_956 = buffer.data(ki + 956);
    const auto *ki_957 = buffer.data(ki + 957);
    const auto *ki_958 = buffer.data(ki + 958);
    const auto *ki_959 = buffer.data(ki + 959);
    const auto *ki_960 = buffer.data(ki + 960);
    const auto *ki_961 = buffer.data(ki + 961);
    const auto *ki_962 = buffer.data(ki + 962);
    const auto *ki_963 = buffer.data(ki + 963);
    const auto *ki_964 = buffer.data(ki + 964);
    const auto *ki_965 = buffer.data(ki + 965);
    const auto *ki_966 = buffer.data(ki + 966);
    const auto *ki_967 = buffer.data(ki + 967);
    const auto *ki_968 = buffer.data(ki + 968);
    const auto *ki_969 = buffer.data(ki + 969);
    const auto *ki_970 = buffer.data(ki + 970);
    const auto *ki_971 = buffer.data(ki + 971);
    const auto *ki_972 = buffer.data(ki + 972);
    const auto *ki_973 = buffer.data(ki + 973);
    const auto *ki_974 = buffer.data(ki + 974);
    const auto *ki_975 = buffer.data(ki + 975);
    const auto *ki_976 = buffer.data(ki + 976);
    const auto *ki_977 = buffer.data(ki + 977);
    const auto *ki_978 = buffer.data(ki + 978);
    const auto *ki_979 = buffer.data(ki + 979);
    const auto *ki_980 = buffer.data(ki + 980);
    const auto *ki_981 = buffer.data(ki + 981);
    const auto *ki_982 = buffer.data(ki + 982);
    const auto *ki_983 = buffer.data(ki + 983);
    const auto *ki_984 = buffer.data(ki + 984);
    const auto *ki_985 = buffer.data(ki + 985);
    const auto *ki_986 = buffer.data(ki + 986);
    const auto *ki_987 = buffer.data(ki + 987);
    const auto *ki_988 = buffer.data(ki + 988);
    const auto *ki_989 = buffer.data(ki + 989);
    const auto *ki_990 = buffer.data(ki + 990);
    const auto *ki_991 = buffer.data(ki + 991);
    const auto *ki_992 = buffer.data(ki + 992);
    const auto *ki_993 = buffer.data(ki + 993);
    const auto *ki_994 = buffer.data(ki + 994);
    const auto *ki_995 = buffer.data(ki + 995);
    const auto *ki_996 = buffer.data(ki + 996);
    const auto *ki_997 = buffer.data(ki + 997);
    const auto *ki_998 = buffer.data(ki + 998);
    const auto *ki_999 = buffer.data(ki + 999);
    const auto *ki_1000 = buffer.data(ki + 1000);
    const auto *ki_1001 = buffer.data(ki + 1001);
    const auto *ki_1002 = buffer.data(ki + 1002);
    const auto *ki_1003 = buffer.data(ki + 1003);
    const auto *ki_1004 = buffer.data(ki + 1004);
    const auto *ki_1005 = buffer.data(ki + 1005);
    const auto *ki_1006 = buffer.data(ki + 1006);
    const auto *ki_1007 = buffer.data(ki + 1007);

#pragma omp simd aligned(t_667, t_668, t_669, t_670, t_671, hi_471, hi_472, hi_473, hi_474, \
                         hi_475, ki_891, ki_892, ki_893, ki_894, \
                         ki_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = -2.0 * hi_471[k]
                   + f_0 * ki_891[k];

        t_668[k] = -2.0 * hi_472[k]
                   + f_0 * ki_892[k];

        t_669[k] = -2.0 * hi_473[k]
                   + f_0 * ki_893[k];

        t_670[k] = -2.0 * hi_474[k]
                   + f_0 * ki_894[k];

        t_671[k] = -2.0 * hi_475[k]
                   + f_0 * ki_895[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, t_676, hi_476, hi_477, hi_478, hi_479, \
                         hi_480, ki_896, ki_897, ki_898, ki_899, \
                         ki_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = -3.0 * hi_476[k]
                   + f_0 * ki_896[k];

        t_673[k] = -3.0 * hi_477[k]
                   + f_0 * ki_897[k];

        t_674[k] = -3.0 * hi_478[k]
                   + f_0 * ki_898[k];

        t_675[k] = -3.0 * hi_479[k]
                   + f_0 * ki_899[k];

        t_676[k] = -3.0 * hi_480[k]
                   + f_0 * ki_900[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, t_680, t_681, hi_481, hi_482, hi_483, hi_484, \
                         hi_485, ki_901, ki_902, ki_903, ki_904, \
                         ki_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = -3.0 * hi_481[k]
                   + f_0 * ki_901[k];

        t_678[k] = -3.0 * hi_482[k]
                   + f_0 * ki_902[k];

        t_679[k] = -3.0 * hi_483[k]
                   + f_0 * ki_903[k];

        t_680[k] = -3.0 * hi_484[k]
                   + f_0 * ki_904[k];

        t_681[k] = -3.0 * hi_485[k]
                   + f_0 * ki_905[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, t_686, hi_486, hi_487, hi_488, hi_489, \
                         hi_490, ki_906, ki_907, ki_908, ki_909, \
                         ki_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = -3.0 * hi_486[k]
                   + f_0 * ki_906[k];

        t_683[k] = -3.0 * hi_487[k]
                   + f_0 * ki_907[k];

        t_684[k] = -3.0 * hi_488[k]
                   + f_0 * ki_908[k];

        t_685[k] = -3.0 * hi_489[k]
                   + f_0 * ki_909[k];

        t_686[k] = -3.0 * hi_490[k]
                   + f_0 * ki_910[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, t_691, hi_491, hi_492, hi_493, hi_494, \
                         hi_495, ki_911, ki_912, ki_913, ki_914, \
                         ki_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = -3.0 * hi_491[k]
                   + f_0 * ki_911[k];

        t_688[k] = -3.0 * hi_492[k]
                   + f_0 * ki_912[k];

        t_689[k] = -3.0 * hi_493[k]
                   + f_0 * ki_913[k];

        t_690[k] = -3.0 * hi_494[k]
                   + f_0 * ki_914[k];

        t_691[k] = -3.0 * hi_495[k]
                   + f_0 * ki_915[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, t_696, hi_496, hi_497, hi_498, hi_499, \
                         hi_500, ki_916, ki_917, ki_918, ki_919, \
                         ki_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = -3.0 * hi_496[k]
                   + f_0 * ki_916[k];

        t_693[k] = -3.0 * hi_497[k]
                   + f_0 * ki_917[k];

        t_694[k] = -3.0 * hi_498[k]
                   + f_0 * ki_918[k];

        t_695[k] = -3.0 * hi_499[k]
                   + f_0 * ki_919[k];

        t_696[k] = -3.0 * hi_500[k]
                   + f_0 * ki_920[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, t_701, hi_501, hi_502, hi_503, hi_504, \
                         hi_505, ki_921, ki_922, ki_923, ki_924, \
                         ki_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = -3.0 * hi_501[k]
                   + f_0 * ki_921[k];

        t_698[k] = -3.0 * hi_502[k]
                   + f_0 * ki_922[k];

        t_699[k] = -3.0 * hi_503[k]
                   + f_0 * ki_923[k];

        t_700[k] = -4.0 * hi_504[k]
                   + f_0 * ki_924[k];

        t_701[k] = -4.0 * hi_505[k]
                   + f_0 * ki_925[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, t_706, hi_506, hi_507, hi_508, hi_509, \
                         hi_510, ki_926, ki_927, ki_928, ki_929, \
                         ki_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = -4.0 * hi_506[k]
                   + f_0 * ki_926[k];

        t_703[k] = -4.0 * hi_507[k]
                   + f_0 * ki_927[k];

        t_704[k] = -4.0 * hi_508[k]
                   + f_0 * ki_928[k];

        t_705[k] = -4.0 * hi_509[k]
                   + f_0 * ki_929[k];

        t_706[k] = -4.0 * hi_510[k]
                   + f_0 * ki_930[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, t_710, t_711, hi_511, hi_512, hi_513, hi_514, \
                         hi_515, ki_931, ki_932, ki_933, ki_934, \
                         ki_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = -4.0 * hi_511[k]
                   + f_0 * ki_931[k];

        t_708[k] = -4.0 * hi_512[k]
                   + f_0 * ki_932[k];

        t_709[k] = -4.0 * hi_513[k]
                   + f_0 * ki_933[k];

        t_710[k] = -4.0 * hi_514[k]
                   + f_0 * ki_934[k];

        t_711[k] = -4.0 * hi_515[k]
                   + f_0 * ki_935[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, t_716, hi_516, hi_517, hi_518, hi_519, \
                         hi_520, ki_936, ki_937, ki_938, ki_939, \
                         ki_940 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = -4.0 * hi_516[k]
                   + f_0 * ki_936[k];

        t_713[k] = -4.0 * hi_517[k]
                   + f_0 * ki_937[k];

        t_714[k] = -4.0 * hi_518[k]
                   + f_0 * ki_938[k];

        t_715[k] = -4.0 * hi_519[k]
                   + f_0 * ki_939[k];

        t_716[k] = -4.0 * hi_520[k]
                   + f_0 * ki_940[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, t_720, t_721, hi_521, hi_522, hi_523, hi_524, \
                         hi_525, ki_941, ki_942, ki_943, ki_944, \
                         ki_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = -4.0 * hi_521[k]
                   + f_0 * ki_941[k];

        t_718[k] = -4.0 * hi_522[k]
                   + f_0 * ki_942[k];

        t_719[k] = -4.0 * hi_523[k]
                   + f_0 * ki_943[k];

        t_720[k] = -4.0 * hi_524[k]
                   + f_0 * ki_944[k];

        t_721[k] = -4.0 * hi_525[k]
                   + f_0 * ki_945[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, t_725, t_726, hi_526, hi_527, hi_528, hi_529, \
                         hi_530, ki_946, ki_947, ki_948, ki_949, \
                         ki_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = -4.0 * hi_526[k]
                   + f_0 * ki_946[k];

        t_723[k] = -4.0 * hi_527[k]
                   + f_0 * ki_947[k];

        t_724[k] = -4.0 * hi_528[k]
                   + f_0 * ki_948[k];

        t_725[k] = -4.0 * hi_529[k]
                   + f_0 * ki_949[k];

        t_726[k] = -4.0 * hi_530[k]
                   + f_0 * ki_950[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, t_730, t_731, hi_531, hi_532, hi_533, hi_534, \
                         hi_535, ki_951, ki_952, ki_953, ki_954, \
                         ki_955 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = -4.0 * hi_531[k]
                   + f_0 * ki_951[k];

        t_728[k] = -5.0 * hi_532[k]
                   + f_0 * ki_952[k];

        t_729[k] = -5.0 * hi_533[k]
                   + f_0 * ki_953[k];

        t_730[k] = -5.0 * hi_534[k]
                   + f_0 * ki_954[k];

        t_731[k] = -5.0 * hi_535[k]
                   + f_0 * ki_955[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, t_735, t_736, hi_536, hi_537, hi_538, hi_539, \
                         hi_540, ki_956, ki_957, ki_958, ki_959, \
                         ki_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = -5.0 * hi_536[k]
                   + f_0 * ki_956[k];

        t_733[k] = -5.0 * hi_537[k]
                   + f_0 * ki_957[k];

        t_734[k] = -5.0 * hi_538[k]
                   + f_0 * ki_958[k];

        t_735[k] = -5.0 * hi_539[k]
                   + f_0 * ki_959[k];

        t_736[k] = -5.0 * hi_540[k]
                   + f_0 * ki_960[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, t_740, t_741, hi_541, hi_542, hi_543, hi_544, \
                         hi_545, ki_961, ki_962, ki_963, ki_964, \
                         ki_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = -5.0 * hi_541[k]
                   + f_0 * ki_961[k];

        t_738[k] = -5.0 * hi_542[k]
                   + f_0 * ki_962[k];

        t_739[k] = -5.0 * hi_543[k]
                   + f_0 * ki_963[k];

        t_740[k] = -5.0 * hi_544[k]
                   + f_0 * ki_964[k];

        t_741[k] = -5.0 * hi_545[k]
                   + f_0 * ki_965[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, t_746, hi_546, hi_547, hi_548, hi_549, \
                         hi_550, ki_966, ki_967, ki_968, ki_969, \
                         ki_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = -5.0 * hi_546[k]
                   + f_0 * ki_966[k];

        t_743[k] = -5.0 * hi_547[k]
                   + f_0 * ki_967[k];

        t_744[k] = -5.0 * hi_548[k]
                   + f_0 * ki_968[k];

        t_745[k] = -5.0 * hi_549[k]
                   + f_0 * ki_969[k];

        t_746[k] = -5.0 * hi_550[k]
                   + f_0 * ki_970[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, t_750, t_751, hi_551, hi_552, hi_553, hi_554, \
                         hi_555, ki_971, ki_972, ki_973, ki_974, \
                         ki_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = -5.0 * hi_551[k]
                   + f_0 * ki_971[k];

        t_748[k] = -5.0 * hi_552[k]
                   + f_0 * ki_972[k];

        t_749[k] = -5.0 * hi_553[k]
                   + f_0 * ki_973[k];

        t_750[k] = -5.0 * hi_554[k]
                   + f_0 * ki_974[k];

        t_751[k] = -5.0 * hi_555[k]
                   + f_0 * ki_975[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, t_756, hi_556, hi_557, hi_558, hi_559, \
                         hi_560, ki_976, ki_977, ki_978, ki_979, \
                         ki_980 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = -5.0 * hi_556[k]
                   + f_0 * ki_976[k];

        t_753[k] = -5.0 * hi_557[k]
                   + f_0 * ki_977[k];

        t_754[k] = -5.0 * hi_558[k]
                   + f_0 * ki_978[k];

        t_755[k] = -5.0 * hi_559[k]
                   + f_0 * ki_979[k];

        t_756[k] = -6.0 * hi_560[k]
                   + f_0 * ki_980[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, t_761, hi_561, hi_562, hi_563, hi_564, \
                         hi_565, ki_981, ki_982, ki_983, ki_984, \
                         ki_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = -6.0 * hi_561[k]
                   + f_0 * ki_981[k];

        t_758[k] = -6.0 * hi_562[k]
                   + f_0 * ki_982[k];

        t_759[k] = -6.0 * hi_563[k]
                   + f_0 * ki_983[k];

        t_760[k] = -6.0 * hi_564[k]
                   + f_0 * ki_984[k];

        t_761[k] = -6.0 * hi_565[k]
                   + f_0 * ki_985[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, t_766, hi_566, hi_567, hi_568, hi_569, \
                         hi_570, ki_986, ki_987, ki_988, ki_989, \
                         ki_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = -6.0 * hi_566[k]
                   + f_0 * ki_986[k];

        t_763[k] = -6.0 * hi_567[k]
                   + f_0 * ki_987[k];

        t_764[k] = -6.0 * hi_568[k]
                   + f_0 * ki_988[k];

        t_765[k] = -6.0 * hi_569[k]
                   + f_0 * ki_989[k];

        t_766[k] = -6.0 * hi_570[k]
                   + f_0 * ki_990[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, t_770, t_771, hi_571, hi_572, hi_573, hi_574, \
                         hi_575, ki_991, ki_992, ki_993, ki_994, \
                         ki_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = -6.0 * hi_571[k]
                   + f_0 * ki_991[k];

        t_768[k] = -6.0 * hi_572[k]
                   + f_0 * ki_992[k];

        t_769[k] = -6.0 * hi_573[k]
                   + f_0 * ki_993[k];

        t_770[k] = -6.0 * hi_574[k]
                   + f_0 * ki_994[k];

        t_771[k] = -6.0 * hi_575[k]
                   + f_0 * ki_995[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, hi_576, hi_577, hi_578, hi_579, \
                         hi_580, ki_996, ki_997, ki_998, ki_999, \
                         ki_1000 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = -6.0 * hi_576[k]
                   + f_0 * ki_996[k];

        t_773[k] = -6.0 * hi_577[k]
                   + f_0 * ki_997[k];

        t_774[k] = -6.0 * hi_578[k]
                   + f_0 * ki_998[k];

        t_775[k] = -6.0 * hi_579[k]
                   + f_0 * ki_999[k];

        t_776[k] = -6.0 * hi_580[k]
                   + f_0 * ki_1000[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, t_780, t_781, hi_581, hi_582, hi_583, hi_584, \
                         hi_585, ki_1001, ki_1002, ki_1003, ki_1004, \
                         ki_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = -6.0 * hi_581[k]
                   + f_0 * ki_1001[k];

        t_778[k] = -6.0 * hi_582[k]
                   + f_0 * ki_1002[k];

        t_779[k] = -6.0 * hi_583[k]
                   + f_0 * ki_1003[k];

        t_780[k] = -6.0 * hi_584[k]
                   + f_0 * ki_1004[k];

        t_781[k] = -6.0 * hi_585[k]
                   + f_0 * ki_1005[k];
    }

#pragma omp simd aligned(t_782, t_783, hi_586, hi_587, ki_1006, \
                         ki_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = -6.0 * hi_586[k]
                   + f_0 * ki_1006[k];

        t_783[k] = -6.0 * hi_587[k]
                   + f_0 * ki_1007[k];
    }
}

auto
compute_prim_geom_10_ii_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t hi, const size_t ki,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ii_electron_repulsion_2_piece0(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_2_piece1(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_2_piece2(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_2_piece3(buffer, target, hi, ki, ncols, alpha);

    compute_prim_geom_10_ii_electron_repulsion_2_piece4(buffer, target, hi, ki, ncols, alpha);
}

}  // namespace simdt2ceri
