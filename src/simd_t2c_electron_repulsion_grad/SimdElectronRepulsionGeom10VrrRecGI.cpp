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


#include "SimdElectronRepulsionGeom10VrrRecGI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_gi_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fi, const size_t hi,
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

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_41 = buffer.data(fi + 41);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_60 = buffer.data(fi + 60);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_63 = buffer.data(fi + 63);
    const auto *fi_64 = buffer.data(fi + 64);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_66 = buffer.data(fi + 66);
    const auto *fi_67 = buffer.data(fi + 67);
    const auto *fi_68 = buffer.data(fi + 68);
    const auto *fi_69 = buffer.data(fi + 69);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_71 = buffer.data(fi + 71);
    const auto *fi_72 = buffer.data(fi + 72);
    const auto *fi_73 = buffer.data(fi + 73);
    const auto *fi_74 = buffer.data(fi + 74);
    const auto *fi_75 = buffer.data(fi + 75);
    const auto *fi_76 = buffer.data(fi + 76);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_78 = buffer.data(fi + 78);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);
    const auto *fi_81 = buffer.data(fi + 81);
    const auto *fi_82 = buffer.data(fi + 82);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_84 = buffer.data(fi + 84);
    const auto *fi_85 = buffer.data(fi + 85);
    const auto *fi_86 = buffer.data(fi + 86);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_88 = buffer.data(fi + 88);
    const auto *fi_89 = buffer.data(fi + 89);
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_91 = buffer.data(fi + 91);
    const auto *fi_92 = buffer.data(fi + 92);
    const auto *fi_93 = buffer.data(fi + 93);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_95 = buffer.data(fi + 95);
    const auto *fi_96 = buffer.data(fi + 96);
    const auto *fi_97 = buffer.data(fi + 97);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_100 = buffer.data(fi + 100);
    const auto *fi_101 = buffer.data(fi + 101);
    const auto *fi_102 = buffer.data(fi + 102);
    const auto *fi_103 = buffer.data(fi + 103);
    const auto *fi_104 = buffer.data(fi + 104);
    const auto *fi_105 = buffer.data(fi + 105);
    const auto *fi_106 = buffer.data(fi + 106);
    const auto *fi_107 = buffer.data(fi + 107);
    const auto *fi_108 = buffer.data(fi + 108);
    const auto *fi_109 = buffer.data(fi + 109);
    const auto *fi_110 = buffer.data(fi + 110);
    const auto *fi_111 = buffer.data(fi + 111);
    const auto *fi_112 = buffer.data(fi + 112);
    const auto *fi_113 = buffer.data(fi + 113);
    const auto *fi_114 = buffer.data(fi + 114);
    const auto *fi_115 = buffer.data(fi + 115);
    const auto *fi_116 = buffer.data(fi + 116);
    const auto *fi_117 = buffer.data(fi + 117);
    const auto *fi_118 = buffer.data(fi + 118);
    const auto *fi_119 = buffer.data(fi + 119);
    const auto *fi_120 = buffer.data(fi + 120);
    const auto *fi_121 = buffer.data(fi + 121);
    const auto *fi_122 = buffer.data(fi + 122);
    const auto *fi_123 = buffer.data(fi + 123);
    const auto *fi_124 = buffer.data(fi + 124);
    const auto *fi_125 = buffer.data(fi + 125);
    const auto *fi_126 = buffer.data(fi + 126);
    const auto *fi_127 = buffer.data(fi + 127);
    const auto *fi_128 = buffer.data(fi + 128);
    const auto *fi_129 = buffer.data(fi + 129);
    const auto *fi_130 = buffer.data(fi + 130);
    const auto *fi_131 = buffer.data(fi + 131);
    const auto *fi_132 = buffer.data(fi + 132);
    const auto *fi_133 = buffer.data(fi + 133);
    const auto *fi_134 = buffer.data(fi + 134);
    const auto *fi_135 = buffer.data(fi + 135);
    const auto *fi_136 = buffer.data(fi + 136);
    const auto *fi_137 = buffer.data(fi + 137);
    const auto *fi_138 = buffer.data(fi + 138);
    const auto *fi_139 = buffer.data(fi + 139);
    const auto *fi_140 = buffer.data(fi + 140);
    const auto *fi_141 = buffer.data(fi + 141);
    const auto *fi_142 = buffer.data(fi + 142);
    const auto *fi_143 = buffer.data(fi + 143);
    const auto *fi_144 = buffer.data(fi + 144);
    const auto *fi_145 = buffer.data(fi + 145);
    const auto *fi_146 = buffer.data(fi + 146);
    const auto *fi_147 = buffer.data(fi + 147);
    const auto *fi_148 = buffer.data(fi + 148);
    const auto *fi_149 = buffer.data(fi + 149);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, fi_0, fi_1, fi_2, fi_3, fi_4, hi_0, hi_1, \
                         hi_2, hi_3, hi_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -4.0 * fi_0[k]
                 + f_0 * hi_0[k];

        t_1[k] = -4.0 * fi_1[k]
                 + f_0 * hi_1[k];

        t_2[k] = -4.0 * fi_2[k]
                 + f_0 * hi_2[k];

        t_3[k] = -4.0 * fi_3[k]
                 + f_0 * hi_3[k];

        t_4[k] = -4.0 * fi_4[k]
                 + f_0 * hi_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, fi_5, fi_6, fi_7, fi_8, fi_9, hi_5, hi_6, \
                         hi_7, hi_8, hi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -4.0 * fi_5[k]
                 + f_0 * hi_5[k];

        t_6[k] = -4.0 * fi_6[k]
                 + f_0 * hi_6[k];

        t_7[k] = -4.0 * fi_7[k]
                 + f_0 * hi_7[k];

        t_8[k] = -4.0 * fi_8[k]
                 + f_0 * hi_8[k];

        t_9[k] = -4.0 * fi_9[k]
                 + f_0 * hi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, fi_10, fi_11, fi_12, fi_13, fi_14, \
                         hi_10, hi_11, hi_12, hi_13, hi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -4.0 * fi_10[k]
                  + f_0 * hi_10[k];

        t_11[k] = -4.0 * fi_11[k]
                  + f_0 * hi_11[k];

        t_12[k] = -4.0 * fi_12[k]
                  + f_0 * hi_12[k];

        t_13[k] = -4.0 * fi_13[k]
                  + f_0 * hi_13[k];

        t_14[k] = -4.0 * fi_14[k]
                  + f_0 * hi_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, fi_15, fi_16, fi_17, fi_18, fi_19, \
                         hi_15, hi_16, hi_17, hi_18, hi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -4.0 * fi_15[k]
                  + f_0 * hi_15[k];

        t_16[k] = -4.0 * fi_16[k]
                  + f_0 * hi_16[k];

        t_17[k] = -4.0 * fi_17[k]
                  + f_0 * hi_17[k];

        t_18[k] = -4.0 * fi_18[k]
                  + f_0 * hi_18[k];

        t_19[k] = -4.0 * fi_19[k]
                  + f_0 * hi_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, fi_20, fi_21, fi_22, fi_23, fi_24, \
                         hi_20, hi_21, hi_22, hi_23, hi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -4.0 * fi_20[k]
                  + f_0 * hi_20[k];

        t_21[k] = -4.0 * fi_21[k]
                  + f_0 * hi_21[k];

        t_22[k] = -4.0 * fi_22[k]
                  + f_0 * hi_22[k];

        t_23[k] = -4.0 * fi_23[k]
                  + f_0 * hi_23[k];

        t_24[k] = -4.0 * fi_24[k]
                  + f_0 * hi_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, fi_25, fi_26, fi_27, fi_28, fi_29, \
                         hi_25, hi_26, hi_27, hi_28, hi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -4.0 * fi_25[k]
                  + f_0 * hi_25[k];

        t_26[k] = -4.0 * fi_26[k]
                  + f_0 * hi_26[k];

        t_27[k] = -4.0 * fi_27[k]
                  + f_0 * hi_27[k];

        t_28[k] = -3.0 * fi_28[k]
                  + f_0 * hi_28[k];

        t_29[k] = -3.0 * fi_29[k]
                  + f_0 * hi_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, fi_30, fi_31, fi_32, fi_33, fi_34, \
                         hi_30, hi_31, hi_32, hi_33, hi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -3.0 * fi_30[k]
                  + f_0 * hi_30[k];

        t_31[k] = -3.0 * fi_31[k]
                  + f_0 * hi_31[k];

        t_32[k] = -3.0 * fi_32[k]
                  + f_0 * hi_32[k];

        t_33[k] = -3.0 * fi_33[k]
                  + f_0 * hi_33[k];

        t_34[k] = -3.0 * fi_34[k]
                  + f_0 * hi_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, fi_35, fi_36, fi_37, fi_38, fi_39, \
                         hi_35, hi_36, hi_37, hi_38, hi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -3.0 * fi_35[k]
                  + f_0 * hi_35[k];

        t_36[k] = -3.0 * fi_36[k]
                  + f_0 * hi_36[k];

        t_37[k] = -3.0 * fi_37[k]
                  + f_0 * hi_37[k];

        t_38[k] = -3.0 * fi_38[k]
                  + f_0 * hi_38[k];

        t_39[k] = -3.0 * fi_39[k]
                  + f_0 * hi_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, fi_40, fi_41, fi_42, fi_43, fi_44, \
                         hi_40, hi_41, hi_42, hi_43, hi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * fi_40[k]
                  + f_0 * hi_40[k];

        t_41[k] = -3.0 * fi_41[k]
                  + f_0 * hi_41[k];

        t_42[k] = -3.0 * fi_42[k]
                  + f_0 * hi_42[k];

        t_43[k] = -3.0 * fi_43[k]
                  + f_0 * hi_43[k];

        t_44[k] = -3.0 * fi_44[k]
                  + f_0 * hi_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, fi_45, fi_46, fi_47, fi_48, fi_49, \
                         hi_45, hi_46, hi_47, hi_48, hi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -3.0 * fi_45[k]
                  + f_0 * hi_45[k];

        t_46[k] = -3.0 * fi_46[k]
                  + f_0 * hi_46[k];

        t_47[k] = -3.0 * fi_47[k]
                  + f_0 * hi_47[k];

        t_48[k] = -3.0 * fi_48[k]
                  + f_0 * hi_48[k];

        t_49[k] = -3.0 * fi_49[k]
                  + f_0 * hi_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, fi_50, fi_51, fi_52, fi_53, fi_54, \
                         hi_50, hi_51, hi_52, hi_53, hi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -3.0 * fi_50[k]
                  + f_0 * hi_50[k];

        t_51[k] = -3.0 * fi_51[k]
                  + f_0 * hi_51[k];

        t_52[k] = -3.0 * fi_52[k]
                  + f_0 * hi_52[k];

        t_53[k] = -3.0 * fi_53[k]
                  + f_0 * hi_53[k];

        t_54[k] = -3.0 * fi_54[k]
                  + f_0 * hi_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, fi_55, fi_56, fi_57, fi_58, fi_59, \
                         hi_55, hi_56, hi_57, hi_58, hi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -3.0 * fi_55[k]
                  + f_0 * hi_55[k];

        t_56[k] = -3.0 * fi_56[k]
                  + f_0 * hi_56[k];

        t_57[k] = -3.0 * fi_57[k]
                  + f_0 * hi_57[k];

        t_58[k] = -3.0 * fi_58[k]
                  + f_0 * hi_58[k];

        t_59[k] = -3.0 * fi_59[k]
                  + f_0 * hi_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, fi_60, fi_61, fi_62, fi_63, fi_64, \
                         hi_60, hi_61, hi_62, hi_63, hi_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * fi_60[k]
                  + f_0 * hi_60[k];

        t_61[k] = -3.0 * fi_61[k]
                  + f_0 * hi_61[k];

        t_62[k] = -3.0 * fi_62[k]
                  + f_0 * hi_62[k];

        t_63[k] = -3.0 * fi_63[k]
                  + f_0 * hi_63[k];

        t_64[k] = -3.0 * fi_64[k]
                  + f_0 * hi_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, fi_65, fi_66, fi_67, fi_68, fi_69, \
                         hi_65, hi_66, hi_67, hi_68, hi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * fi_65[k]
                  + f_0 * hi_65[k];

        t_66[k] = -3.0 * fi_66[k]
                  + f_0 * hi_66[k];

        t_67[k] = -3.0 * fi_67[k]
                  + f_0 * hi_67[k];

        t_68[k] = -3.0 * fi_68[k]
                  + f_0 * hi_68[k];

        t_69[k] = -3.0 * fi_69[k]
                  + f_0 * hi_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, fi_70, fi_71, fi_72, fi_73, fi_74, \
                         hi_70, hi_71, hi_72, hi_73, hi_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -3.0 * fi_70[k]
                  + f_0 * hi_70[k];

        t_71[k] = -3.0 * fi_71[k]
                  + f_0 * hi_71[k];

        t_72[k] = -3.0 * fi_72[k]
                  + f_0 * hi_72[k];

        t_73[k] = -3.0 * fi_73[k]
                  + f_0 * hi_73[k];

        t_74[k] = -3.0 * fi_74[k]
                  + f_0 * hi_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, fi_75, fi_76, fi_77, fi_78, fi_79, \
                         hi_75, hi_76, hi_77, hi_78, hi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -3.0 * fi_75[k]
                  + f_0 * hi_75[k];

        t_76[k] = -3.0 * fi_76[k]
                  + f_0 * hi_76[k];

        t_77[k] = -3.0 * fi_77[k]
                  + f_0 * hi_77[k];

        t_78[k] = -3.0 * fi_78[k]
                  + f_0 * hi_78[k];

        t_79[k] = -3.0 * fi_79[k]
                  + f_0 * hi_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, fi_80, fi_81, fi_82, fi_83, fi_84, \
                         hi_80, hi_81, hi_82, hi_83, hi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * fi_80[k]
                  + f_0 * hi_80[k];

        t_81[k] = -3.0 * fi_81[k]
                  + f_0 * hi_81[k];

        t_82[k] = -3.0 * fi_82[k]
                  + f_0 * hi_82[k];

        t_83[k] = -3.0 * fi_83[k]
                  + f_0 * hi_83[k];

        t_84[k] = -2.0 * fi_84[k]
                  + f_0 * hi_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, fi_85, fi_86, fi_87, fi_88, fi_89, \
                         hi_85, hi_86, hi_87, hi_88, hi_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * fi_85[k]
                  + f_0 * hi_85[k];

        t_86[k] = -2.0 * fi_86[k]
                  + f_0 * hi_86[k];

        t_87[k] = -2.0 * fi_87[k]
                  + f_0 * hi_87[k];

        t_88[k] = -2.0 * fi_88[k]
                  + f_0 * hi_88[k];

        t_89[k] = -2.0 * fi_89[k]
                  + f_0 * hi_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, fi_90, fi_91, fi_92, fi_93, fi_94, \
                         hi_90, hi_91, hi_92, hi_93, hi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * fi_90[k]
                  + f_0 * hi_90[k];

        t_91[k] = -2.0 * fi_91[k]
                  + f_0 * hi_91[k];

        t_92[k] = -2.0 * fi_92[k]
                  + f_0 * hi_92[k];

        t_93[k] = -2.0 * fi_93[k]
                  + f_0 * hi_93[k];

        t_94[k] = -2.0 * fi_94[k]
                  + f_0 * hi_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, fi_95, fi_96, fi_97, fi_98, fi_99, \
                         hi_95, hi_96, hi_97, hi_98, hi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * fi_95[k]
                  + f_0 * hi_95[k];

        t_96[k] = -2.0 * fi_96[k]
                  + f_0 * hi_96[k];

        t_97[k] = -2.0 * fi_97[k]
                  + f_0 * hi_97[k];

        t_98[k] = -2.0 * fi_98[k]
                  + f_0 * hi_98[k];

        t_99[k] = -2.0 * fi_99[k]
                  + f_0 * hi_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, fi_100, fi_101, fi_102, fi_103, \
                         fi_104, hi_100, hi_101, hi_102, hi_103, \
                         hi_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * fi_100[k]
                   + f_0 * hi_100[k];

        t_101[k] = -2.0 * fi_101[k]
                   + f_0 * hi_101[k];

        t_102[k] = -2.0 * fi_102[k]
                   + f_0 * hi_102[k];

        t_103[k] = -2.0 * fi_103[k]
                   + f_0 * hi_103[k];

        t_104[k] = -2.0 * fi_104[k]
                   + f_0 * hi_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, fi_105, fi_106, fi_107, fi_108, \
                         fi_109, hi_105, hi_106, hi_107, hi_108, \
                         hi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * fi_105[k]
                   + f_0 * hi_105[k];

        t_106[k] = -2.0 * fi_106[k]
                   + f_0 * hi_106[k];

        t_107[k] = -2.0 * fi_107[k]
                   + f_0 * hi_107[k];

        t_108[k] = -2.0 * fi_108[k]
                   + f_0 * hi_108[k];

        t_109[k] = -2.0 * fi_109[k]
                   + f_0 * hi_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, fi_110, fi_111, fi_112, fi_113, \
                         fi_114, hi_110, hi_111, hi_112, hi_113, \
                         hi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * fi_110[k]
                   + f_0 * hi_110[k];

        t_111[k] = -2.0 * fi_111[k]
                   + f_0 * hi_111[k];

        t_112[k] = -2.0 * fi_112[k]
                   + f_0 * hi_112[k];

        t_113[k] = -2.0 * fi_113[k]
                   + f_0 * hi_113[k];

        t_114[k] = -2.0 * fi_114[k]
                   + f_0 * hi_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, fi_115, fi_116, fi_117, fi_118, \
                         fi_119, hi_115, hi_116, hi_117, hi_118, \
                         hi_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * fi_115[k]
                   + f_0 * hi_115[k];

        t_116[k] = -2.0 * fi_116[k]
                   + f_0 * hi_116[k];

        t_117[k] = -2.0 * fi_117[k]
                   + f_0 * hi_117[k];

        t_118[k] = -2.0 * fi_118[k]
                   + f_0 * hi_118[k];

        t_119[k] = -2.0 * fi_119[k]
                   + f_0 * hi_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, fi_120, fi_121, fi_122, fi_123, \
                         fi_124, hi_120, hi_121, hi_122, hi_123, \
                         hi_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -2.0 * fi_120[k]
                   + f_0 * hi_120[k];

        t_121[k] = -2.0 * fi_121[k]
                   + f_0 * hi_121[k];

        t_122[k] = -2.0 * fi_122[k]
                   + f_0 * hi_122[k];

        t_123[k] = -2.0 * fi_123[k]
                   + f_0 * hi_123[k];

        t_124[k] = -2.0 * fi_124[k]
                   + f_0 * hi_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, fi_125, fi_126, fi_127, fi_128, \
                         fi_129, hi_125, hi_126, hi_127, hi_128, \
                         hi_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -2.0 * fi_125[k]
                   + f_0 * hi_125[k];

        t_126[k] = -2.0 * fi_126[k]
                   + f_0 * hi_126[k];

        t_127[k] = -2.0 * fi_127[k]
                   + f_0 * hi_127[k];

        t_128[k] = -2.0 * fi_128[k]
                   + f_0 * hi_128[k];

        t_129[k] = -2.0 * fi_129[k]
                   + f_0 * hi_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, fi_130, fi_131, fi_132, fi_133, \
                         fi_134, hi_130, hi_131, hi_132, hi_133, \
                         hi_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -2.0 * fi_130[k]
                   + f_0 * hi_130[k];

        t_131[k] = -2.0 * fi_131[k]
                   + f_0 * hi_131[k];

        t_132[k] = -2.0 * fi_132[k]
                   + f_0 * hi_132[k];

        t_133[k] = -2.0 * fi_133[k]
                   + f_0 * hi_133[k];

        t_134[k] = -2.0 * fi_134[k]
                   + f_0 * hi_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, fi_135, fi_136, fi_137, fi_138, \
                         fi_139, hi_135, hi_136, hi_137, hi_138, \
                         hi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -2.0 * fi_135[k]
                   + f_0 * hi_135[k];

        t_136[k] = -2.0 * fi_136[k]
                   + f_0 * hi_136[k];

        t_137[k] = -2.0 * fi_137[k]
                   + f_0 * hi_137[k];

        t_138[k] = -2.0 * fi_138[k]
                   + f_0 * hi_138[k];

        t_139[k] = -2.0 * fi_139[k]
                   + f_0 * hi_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, fi_140, fi_141, fi_142, fi_143, \
                         fi_144, hi_140, hi_141, hi_142, hi_143, \
                         hi_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -2.0 * fi_140[k]
                   + f_0 * hi_140[k];

        t_141[k] = -2.0 * fi_141[k]
                   + f_0 * hi_141[k];

        t_142[k] = -2.0 * fi_142[k]
                   + f_0 * hi_142[k];

        t_143[k] = -2.0 * fi_143[k]
                   + f_0 * hi_143[k];

        t_144[k] = -2.0 * fi_144[k]
                   + f_0 * hi_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, fi_145, fi_146, fi_147, fi_148, \
                         fi_149, hi_145, hi_146, hi_147, hi_148, \
                         hi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -2.0 * fi_145[k]
                   + f_0 * hi_145[k];

        t_146[k] = -2.0 * fi_146[k]
                   + f_0 * hi_146[k];

        t_147[k] = -2.0 * fi_147[k]
                   + f_0 * hi_147[k];

        t_148[k] = -2.0 * fi_148[k]
                   + f_0 * hi_148[k];

        t_149[k] = -2.0 * fi_149[k]
                   + f_0 * hi_149[k];
    }
}

static auto
compute_prim_geom_10_gi_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fi, const size_t hi,
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

    const auto *fi_150 = buffer.data(fi + 150);
    const auto *fi_151 = buffer.data(fi + 151);
    const auto *fi_152 = buffer.data(fi + 152);
    const auto *fi_153 = buffer.data(fi + 153);
    const auto *fi_154 = buffer.data(fi + 154);
    const auto *fi_155 = buffer.data(fi + 155);
    const auto *fi_156 = buffer.data(fi + 156);
    const auto *fi_157 = buffer.data(fi + 157);
    const auto *fi_158 = buffer.data(fi + 158);
    const auto *fi_159 = buffer.data(fi + 159);
    const auto *fi_160 = buffer.data(fi + 160);
    const auto *fi_161 = buffer.data(fi + 161);
    const auto *fi_162 = buffer.data(fi + 162);
    const auto *fi_163 = buffer.data(fi + 163);
    const auto *fi_164 = buffer.data(fi + 164);
    const auto *fi_165 = buffer.data(fi + 165);
    const auto *fi_166 = buffer.data(fi + 166);
    const auto *fi_167 = buffer.data(fi + 167);
    const auto *fi_168 = buffer.data(fi + 168);
    const auto *fi_169 = buffer.data(fi + 169);
    const auto *fi_170 = buffer.data(fi + 170);
    const auto *fi_171 = buffer.data(fi + 171);
    const auto *fi_172 = buffer.data(fi + 172);
    const auto *fi_173 = buffer.data(fi + 173);
    const auto *fi_174 = buffer.data(fi + 174);
    const auto *fi_175 = buffer.data(fi + 175);
    const auto *fi_176 = buffer.data(fi + 176);
    const auto *fi_177 = buffer.data(fi + 177);
    const auto *fi_178 = buffer.data(fi + 178);
    const auto *fi_179 = buffer.data(fi + 179);
    const auto *fi_180 = buffer.data(fi + 180);
    const auto *fi_181 = buffer.data(fi + 181);
    const auto *fi_182 = buffer.data(fi + 182);
    const auto *fi_183 = buffer.data(fi + 183);
    const auto *fi_184 = buffer.data(fi + 184);
    const auto *fi_185 = buffer.data(fi + 185);
    const auto *fi_186 = buffer.data(fi + 186);
    const auto *fi_187 = buffer.data(fi + 187);
    const auto *fi_188 = buffer.data(fi + 188);
    const auto *fi_189 = buffer.data(fi + 189);
    const auto *fi_190 = buffer.data(fi + 190);
    const auto *fi_191 = buffer.data(fi + 191);
    const auto *fi_192 = buffer.data(fi + 192);
    const auto *fi_193 = buffer.data(fi + 193);
    const auto *fi_194 = buffer.data(fi + 194);
    const auto *fi_195 = buffer.data(fi + 195);
    const auto *fi_196 = buffer.data(fi + 196);
    const auto *fi_197 = buffer.data(fi + 197);
    const auto *fi_198 = buffer.data(fi + 198);
    const auto *fi_199 = buffer.data(fi + 199);
    const auto *fi_200 = buffer.data(fi + 200);
    const auto *fi_201 = buffer.data(fi + 201);
    const auto *fi_202 = buffer.data(fi + 202);
    const auto *fi_203 = buffer.data(fi + 203);
    const auto *fi_204 = buffer.data(fi + 204);
    const auto *fi_205 = buffer.data(fi + 205);
    const auto *fi_206 = buffer.data(fi + 206);
    const auto *fi_207 = buffer.data(fi + 207);
    const auto *fi_208 = buffer.data(fi + 208);
    const auto *fi_209 = buffer.data(fi + 209);
    const auto *fi_210 = buffer.data(fi + 210);
    const auto *fi_211 = buffer.data(fi + 211);
    const auto *fi_212 = buffer.data(fi + 212);
    const auto *fi_213 = buffer.data(fi + 213);
    const auto *fi_214 = buffer.data(fi + 214);
    const auto *fi_215 = buffer.data(fi + 215);
    const auto *fi_216 = buffer.data(fi + 216);
    const auto *fi_217 = buffer.data(fi + 217);
    const auto *fi_218 = buffer.data(fi + 218);
    const auto *fi_219 = buffer.data(fi + 219);
    const auto *fi_220 = buffer.data(fi + 220);
    const auto *fi_221 = buffer.data(fi + 221);
    const auto *fi_222 = buffer.data(fi + 222);
    const auto *fi_223 = buffer.data(fi + 223);
    const auto *fi_224 = buffer.data(fi + 224);
    const auto *fi_225 = buffer.data(fi + 225);
    const auto *fi_226 = buffer.data(fi + 226);
    const auto *fi_227 = buffer.data(fi + 227);
    const auto *fi_228 = buffer.data(fi + 228);
    const auto *fi_229 = buffer.data(fi + 229);
    const auto *fi_230 = buffer.data(fi + 230);
    const auto *fi_231 = buffer.data(fi + 231);
    const auto *fi_232 = buffer.data(fi + 232);
    const auto *fi_233 = buffer.data(fi + 233);
    const auto *fi_234 = buffer.data(fi + 234);
    const auto *fi_235 = buffer.data(fi + 235);
    const auto *fi_236 = buffer.data(fi + 236);
    const auto *fi_237 = buffer.data(fi + 237);
    const auto *fi_238 = buffer.data(fi + 238);
    const auto *fi_239 = buffer.data(fi + 239);
    const auto *fi_240 = buffer.data(fi + 240);
    const auto *fi_241 = buffer.data(fi + 241);
    const auto *fi_242 = buffer.data(fi + 242);
    const auto *fi_243 = buffer.data(fi + 243);
    const auto *fi_244 = buffer.data(fi + 244);
    const auto *fi_245 = buffer.data(fi + 245);
    const auto *fi_246 = buffer.data(fi + 246);
    const auto *fi_247 = buffer.data(fi + 247);
    const auto *fi_248 = buffer.data(fi + 248);
    const auto *fi_249 = buffer.data(fi + 249);
    const auto *fi_250 = buffer.data(fi + 250);
    const auto *fi_251 = buffer.data(fi + 251);
    const auto *fi_252 = buffer.data(fi + 252);
    const auto *fi_253 = buffer.data(fi + 253);
    const auto *fi_254 = buffer.data(fi + 254);
    const auto *fi_255 = buffer.data(fi + 255);
    const auto *fi_256 = buffer.data(fi + 256);
    const auto *fi_257 = buffer.data(fi + 257);
    const auto *fi_258 = buffer.data(fi + 258);
    const auto *fi_259 = buffer.data(fi + 259);
    const auto *fi_260 = buffer.data(fi + 260);
    const auto *fi_261 = buffer.data(fi + 261);
    const auto *fi_262 = buffer.data(fi + 262);
    const auto *fi_263 = buffer.data(fi + 263);
    const auto *fi_264 = buffer.data(fi + 264);
    const auto *fi_265 = buffer.data(fi + 265);
    const auto *fi_266 = buffer.data(fi + 266);
    const auto *fi_267 = buffer.data(fi + 267);
    const auto *fi_268 = buffer.data(fi + 268);
    const auto *fi_269 = buffer.data(fi + 269);
    const auto *fi_270 = buffer.data(fi + 270);
    const auto *fi_271 = buffer.data(fi + 271);
    const auto *fi_272 = buffer.data(fi + 272);
    const auto *fi_273 = buffer.data(fi + 273);
    const auto *fi_274 = buffer.data(fi + 274);
    const auto *fi_275 = buffer.data(fi + 275);
    const auto *fi_276 = buffer.data(fi + 276);
    const auto *fi_277 = buffer.data(fi + 277);
    const auto *fi_278 = buffer.data(fi + 278);
    const auto *fi_279 = buffer.data(fi + 279);

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

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, fi_150, fi_151, fi_152, fi_153, \
                         fi_154, hi_150, hi_151, hi_152, hi_153, \
                         hi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -2.0 * fi_150[k]
                   + f_0 * hi_150[k];

        t_151[k] = -2.0 * fi_151[k]
                   + f_0 * hi_151[k];

        t_152[k] = -2.0 * fi_152[k]
                   + f_0 * hi_152[k];

        t_153[k] = -2.0 * fi_153[k]
                   + f_0 * hi_153[k];

        t_154[k] = -2.0 * fi_154[k]
                   + f_0 * hi_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, fi_155, fi_156, fi_157, fi_158, \
                         fi_159, hi_155, hi_156, hi_157, hi_158, \
                         hi_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -2.0 * fi_155[k]
                   + f_0 * hi_155[k];

        t_156[k] = -2.0 * fi_156[k]
                   + f_0 * hi_156[k];

        t_157[k] = -2.0 * fi_157[k]
                   + f_0 * hi_157[k];

        t_158[k] = -2.0 * fi_158[k]
                   + f_0 * hi_158[k];

        t_159[k] = -2.0 * fi_159[k]
                   + f_0 * hi_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, fi_160, fi_161, fi_162, fi_163, \
                         fi_164, hi_160, hi_161, hi_162, hi_163, \
                         hi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -2.0 * fi_160[k]
                   + f_0 * hi_160[k];

        t_161[k] = -2.0 * fi_161[k]
                   + f_0 * hi_161[k];

        t_162[k] = -2.0 * fi_162[k]
                   + f_0 * hi_162[k];

        t_163[k] = -2.0 * fi_163[k]
                   + f_0 * hi_163[k];

        t_164[k] = -2.0 * fi_164[k]
                   + f_0 * hi_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, fi_165, fi_166, fi_167, fi_168, \
                         fi_169, hi_165, hi_166, hi_167, hi_168, \
                         hi_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -2.0 * fi_165[k]
                   + f_0 * hi_165[k];

        t_166[k] = -2.0 * fi_166[k]
                   + f_0 * hi_166[k];

        t_167[k] = -2.0 * fi_167[k]
                   + f_0 * hi_167[k];

        t_168[k] = -fi_168[k]
                   + f_0 * hi_168[k];

        t_169[k] = -fi_169[k]
                   + f_0 * hi_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, fi_170, fi_171, fi_172, fi_173, \
                         fi_174, hi_170, hi_171, hi_172, hi_173, \
                         hi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -fi_170[k]
                   + f_0 * hi_170[k];

        t_171[k] = -fi_171[k]
                   + f_0 * hi_171[k];

        t_172[k] = -fi_172[k]
                   + f_0 * hi_172[k];

        t_173[k] = -fi_173[k]
                   + f_0 * hi_173[k];

        t_174[k] = -fi_174[k]
                   + f_0 * hi_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, fi_175, fi_176, fi_177, fi_178, \
                         fi_179, hi_175, hi_176, hi_177, hi_178, \
                         hi_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -fi_175[k]
                   + f_0 * hi_175[k];

        t_176[k] = -fi_176[k]
                   + f_0 * hi_176[k];

        t_177[k] = -fi_177[k]
                   + f_0 * hi_177[k];

        t_178[k] = -fi_178[k]
                   + f_0 * hi_178[k];

        t_179[k] = -fi_179[k]
                   + f_0 * hi_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, fi_180, fi_181, fi_182, fi_183, \
                         fi_184, hi_180, hi_181, hi_182, hi_183, \
                         hi_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -fi_180[k]
                   + f_0 * hi_180[k];

        t_181[k] = -fi_181[k]
                   + f_0 * hi_181[k];

        t_182[k] = -fi_182[k]
                   + f_0 * hi_182[k];

        t_183[k] = -fi_183[k]
                   + f_0 * hi_183[k];

        t_184[k] = -fi_184[k]
                   + f_0 * hi_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, fi_185, fi_186, fi_187, fi_188, \
                         fi_189, hi_185, hi_186, hi_187, hi_188, \
                         hi_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -fi_185[k]
                   + f_0 * hi_185[k];

        t_186[k] = -fi_186[k]
                   + f_0 * hi_186[k];

        t_187[k] = -fi_187[k]
                   + f_0 * hi_187[k];

        t_188[k] = -fi_188[k]
                   + f_0 * hi_188[k];

        t_189[k] = -fi_189[k]
                   + f_0 * hi_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, fi_190, fi_191, fi_192, fi_193, \
                         fi_194, hi_190, hi_191, hi_192, hi_193, \
                         hi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -fi_190[k]
                   + f_0 * hi_190[k];

        t_191[k] = -fi_191[k]
                   + f_0 * hi_191[k];

        t_192[k] = -fi_192[k]
                   + f_0 * hi_192[k];

        t_193[k] = -fi_193[k]
                   + f_0 * hi_193[k];

        t_194[k] = -fi_194[k]
                   + f_0 * hi_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, fi_195, fi_196, fi_197, fi_198, \
                         fi_199, hi_195, hi_196, hi_197, hi_198, \
                         hi_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -fi_195[k]
                   + f_0 * hi_195[k];

        t_196[k] = -fi_196[k]
                   + f_0 * hi_196[k];

        t_197[k] = -fi_197[k]
                   + f_0 * hi_197[k];

        t_198[k] = -fi_198[k]
                   + f_0 * hi_198[k];

        t_199[k] = -fi_199[k]
                   + f_0 * hi_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, fi_200, fi_201, fi_202, fi_203, \
                         fi_204, hi_200, hi_201, hi_202, hi_203, \
                         hi_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -fi_200[k]
                   + f_0 * hi_200[k];

        t_201[k] = -fi_201[k]
                   + f_0 * hi_201[k];

        t_202[k] = -fi_202[k]
                   + f_0 * hi_202[k];

        t_203[k] = -fi_203[k]
                   + f_0 * hi_203[k];

        t_204[k] = -fi_204[k]
                   + f_0 * hi_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, fi_205, fi_206, fi_207, fi_208, \
                         fi_209, hi_205, hi_206, hi_207, hi_208, \
                         hi_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -fi_205[k]
                   + f_0 * hi_205[k];

        t_206[k] = -fi_206[k]
                   + f_0 * hi_206[k];

        t_207[k] = -fi_207[k]
                   + f_0 * hi_207[k];

        t_208[k] = -fi_208[k]
                   + f_0 * hi_208[k];

        t_209[k] = -fi_209[k]
                   + f_0 * hi_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, fi_210, fi_211, fi_212, fi_213, \
                         fi_214, hi_210, hi_211, hi_212, hi_213, \
                         hi_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -fi_210[k]
                   + f_0 * hi_210[k];

        t_211[k] = -fi_211[k]
                   + f_0 * hi_211[k];

        t_212[k] = -fi_212[k]
                   + f_0 * hi_212[k];

        t_213[k] = -fi_213[k]
                   + f_0 * hi_213[k];

        t_214[k] = -fi_214[k]
                   + f_0 * hi_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, fi_215, fi_216, fi_217, fi_218, \
                         fi_219, hi_215, hi_216, hi_217, hi_218, \
                         hi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -fi_215[k]
                   + f_0 * hi_215[k];

        t_216[k] = -fi_216[k]
                   + f_0 * hi_216[k];

        t_217[k] = -fi_217[k]
                   + f_0 * hi_217[k];

        t_218[k] = -fi_218[k]
                   + f_0 * hi_218[k];

        t_219[k] = -fi_219[k]
                   + f_0 * hi_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, fi_220, fi_221, fi_222, fi_223, \
                         fi_224, hi_220, hi_221, hi_222, hi_223, \
                         hi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -fi_220[k]
                   + f_0 * hi_220[k];

        t_221[k] = -fi_221[k]
                   + f_0 * hi_221[k];

        t_222[k] = -fi_222[k]
                   + f_0 * hi_222[k];

        t_223[k] = -fi_223[k]
                   + f_0 * hi_223[k];

        t_224[k] = -fi_224[k]
                   + f_0 * hi_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, fi_225, fi_226, fi_227, fi_228, \
                         fi_229, hi_225, hi_226, hi_227, hi_228, \
                         hi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -fi_225[k]
                   + f_0 * hi_225[k];

        t_226[k] = -fi_226[k]
                   + f_0 * hi_226[k];

        t_227[k] = -fi_227[k]
                   + f_0 * hi_227[k];

        t_228[k] = -fi_228[k]
                   + f_0 * hi_228[k];

        t_229[k] = -fi_229[k]
                   + f_0 * hi_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, fi_230, fi_231, fi_232, fi_233, \
                         fi_234, hi_230, hi_231, hi_232, hi_233, \
                         hi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -fi_230[k]
                   + f_0 * hi_230[k];

        t_231[k] = -fi_231[k]
                   + f_0 * hi_231[k];

        t_232[k] = -fi_232[k]
                   + f_0 * hi_232[k];

        t_233[k] = -fi_233[k]
                   + f_0 * hi_233[k];

        t_234[k] = -fi_234[k]
                   + f_0 * hi_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, fi_235, fi_236, fi_237, fi_238, \
                         fi_239, hi_235, hi_236, hi_237, hi_238, \
                         hi_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -fi_235[k]
                   + f_0 * hi_235[k];

        t_236[k] = -fi_236[k]
                   + f_0 * hi_236[k];

        t_237[k] = -fi_237[k]
                   + f_0 * hi_237[k];

        t_238[k] = -fi_238[k]
                   + f_0 * hi_238[k];

        t_239[k] = -fi_239[k]
                   + f_0 * hi_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, fi_240, fi_241, fi_242, fi_243, \
                         fi_244, hi_240, hi_241, hi_242, hi_243, \
                         hi_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -fi_240[k]
                   + f_0 * hi_240[k];

        t_241[k] = -fi_241[k]
                   + f_0 * hi_241[k];

        t_242[k] = -fi_242[k]
                   + f_0 * hi_242[k];

        t_243[k] = -fi_243[k]
                   + f_0 * hi_243[k];

        t_244[k] = -fi_244[k]
                   + f_0 * hi_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, fi_245, fi_246, fi_247, fi_248, \
                         fi_249, hi_245, hi_246, hi_247, hi_248, \
                         hi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -fi_245[k]
                   + f_0 * hi_245[k];

        t_246[k] = -fi_246[k]
                   + f_0 * hi_246[k];

        t_247[k] = -fi_247[k]
                   + f_0 * hi_247[k];

        t_248[k] = -fi_248[k]
                   + f_0 * hi_248[k];

        t_249[k] = -fi_249[k]
                   + f_0 * hi_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, fi_250, fi_251, fi_252, fi_253, \
                         fi_254, hi_250, hi_251, hi_252, hi_253, \
                         hi_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -fi_250[k]
                   + f_0 * hi_250[k];

        t_251[k] = -fi_251[k]
                   + f_0 * hi_251[k];

        t_252[k] = -fi_252[k]
                   + f_0 * hi_252[k];

        t_253[k] = -fi_253[k]
                   + f_0 * hi_253[k];

        t_254[k] = -fi_254[k]
                   + f_0 * hi_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, fi_255, fi_256, fi_257, fi_258, \
                         fi_259, hi_255, hi_256, hi_257, hi_258, \
                         hi_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -fi_255[k]
                   + f_0 * hi_255[k];

        t_256[k] = -fi_256[k]
                   + f_0 * hi_256[k];

        t_257[k] = -fi_257[k]
                   + f_0 * hi_257[k];

        t_258[k] = -fi_258[k]
                   + f_0 * hi_258[k];

        t_259[k] = -fi_259[k]
                   + f_0 * hi_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, fi_260, fi_261, fi_262, fi_263, \
                         fi_264, hi_260, hi_261, hi_262, hi_263, \
                         hi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -fi_260[k]
                   + f_0 * hi_260[k];

        t_261[k] = -fi_261[k]
                   + f_0 * hi_261[k];

        t_262[k] = -fi_262[k]
                   + f_0 * hi_262[k];

        t_263[k] = -fi_263[k]
                   + f_0 * hi_263[k];

        t_264[k] = -fi_264[k]
                   + f_0 * hi_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, fi_265, fi_266, fi_267, fi_268, \
                         fi_269, hi_265, hi_266, hi_267, hi_268, \
                         hi_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -fi_265[k]
                   + f_0 * hi_265[k];

        t_266[k] = -fi_266[k]
                   + f_0 * hi_266[k];

        t_267[k] = -fi_267[k]
                   + f_0 * hi_267[k];

        t_268[k] = -fi_268[k]
                   + f_0 * hi_268[k];

        t_269[k] = -fi_269[k]
                   + f_0 * hi_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, fi_270, fi_271, fi_272, fi_273, \
                         fi_274, hi_270, hi_271, hi_272, hi_273, \
                         hi_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -fi_270[k]
                   + f_0 * hi_270[k];

        t_271[k] = -fi_271[k]
                   + f_0 * hi_271[k];

        t_272[k] = -fi_272[k]
                   + f_0 * hi_272[k];

        t_273[k] = -fi_273[k]
                   + f_0 * hi_273[k];

        t_274[k] = -fi_274[k]
                   + f_0 * hi_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, fi_275, fi_276, fi_277, fi_278, \
                         fi_279, hi_275, hi_276, hi_277, hi_278, \
                         hi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -fi_275[k]
                   + f_0 * hi_275[k];

        t_276[k] = -fi_276[k]
                   + f_0 * hi_276[k];

        t_277[k] = -fi_277[k]
                   + f_0 * hi_277[k];

        t_278[k] = -fi_278[k]
                   + f_0 * hi_278[k];

        t_279[k] = -fi_279[k]
                   + f_0 * hi_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, t_285, t_286, t_287, hi_280, \
                         hi_281, hi_282, hi_283, hi_284, hi_285, hi_286, \
                         hi_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_0 * hi_280[k];

        t_281[k] = f_0 * hi_281[k];

        t_282[k] = f_0 * hi_282[k];

        t_283[k] = f_0 * hi_283[k];

        t_284[k] = f_0 * hi_284[k];

        t_285[k] = f_0 * hi_285[k];

        t_286[k] = f_0 * hi_286[k];

        t_287[k] = f_0 * hi_287[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, t_293, t_294, t_295, hi_288, \
                         hi_289, hi_290, hi_291, hi_292, hi_293, hi_294, \
                         hi_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_0 * hi_288[k];

        t_289[k] = f_0 * hi_289[k];

        t_290[k] = f_0 * hi_290[k];

        t_291[k] = f_0 * hi_291[k];

        t_292[k] = f_0 * hi_292[k];

        t_293[k] = f_0 * hi_293[k];

        t_294[k] = f_0 * hi_294[k];

        t_295[k] = f_0 * hi_295[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, t_301, t_302, t_303, hi_296, \
                         hi_297, hi_298, hi_299, hi_300, hi_301, hi_302, \
                         hi_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_0 * hi_296[k];

        t_297[k] = f_0 * hi_297[k];

        t_298[k] = f_0 * hi_298[k];

        t_299[k] = f_0 * hi_299[k];

        t_300[k] = f_0 * hi_300[k];

        t_301[k] = f_0 * hi_301[k];

        t_302[k] = f_0 * hi_302[k];

        t_303[k] = f_0 * hi_303[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, t_309, t_310, t_311, hi_304, \
                         hi_305, hi_306, hi_307, hi_308, hi_309, hi_310, \
                         hi_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_0 * hi_304[k];

        t_305[k] = f_0 * hi_305[k];

        t_306[k] = f_0 * hi_306[k];

        t_307[k] = f_0 * hi_307[k];

        t_308[k] = f_0 * hi_308[k];

        t_309[k] = f_0 * hi_309[k];

        t_310[k] = f_0 * hi_310[k];

        t_311[k] = f_0 * hi_311[k];
    }
}

static auto
compute_prim_geom_10_gi_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hi, const size_t ncols,
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

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, t_317, t_318, t_319, hi_312, \
                         hi_313, hi_314, hi_315, hi_316, hi_317, hi_318, \
                         hi_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_0 * hi_312[k];

        t_313[k] = f_0 * hi_313[k];

        t_314[k] = f_0 * hi_314[k];

        t_315[k] = f_0 * hi_315[k];

        t_316[k] = f_0 * hi_316[k];

        t_317[k] = f_0 * hi_317[k];

        t_318[k] = f_0 * hi_318[k];

        t_319[k] = f_0 * hi_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, t_325, t_326, t_327, hi_320, \
                         hi_321, hi_322, hi_323, hi_324, hi_325, hi_326, \
                         hi_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_0 * hi_320[k];

        t_321[k] = f_0 * hi_321[k];

        t_322[k] = f_0 * hi_322[k];

        t_323[k] = f_0 * hi_323[k];

        t_324[k] = f_0 * hi_324[k];

        t_325[k] = f_0 * hi_325[k];

        t_326[k] = f_0 * hi_326[k];

        t_327[k] = f_0 * hi_327[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, t_333, t_334, t_335, hi_328, \
                         hi_329, hi_330, hi_331, hi_332, hi_333, hi_334, \
                         hi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_0 * hi_328[k];

        t_329[k] = f_0 * hi_329[k];

        t_330[k] = f_0 * hi_330[k];

        t_331[k] = f_0 * hi_331[k];

        t_332[k] = f_0 * hi_332[k];

        t_333[k] = f_0 * hi_333[k];

        t_334[k] = f_0 * hi_334[k];

        t_335[k] = f_0 * hi_335[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, t_340, t_341, t_342, t_343, hi_336, \
                         hi_337, hi_338, hi_339, hi_340, hi_341, hi_342, \
                         hi_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_0 * hi_336[k];

        t_337[k] = f_0 * hi_337[k];

        t_338[k] = f_0 * hi_338[k];

        t_339[k] = f_0 * hi_339[k];

        t_340[k] = f_0 * hi_340[k];

        t_341[k] = f_0 * hi_341[k];

        t_342[k] = f_0 * hi_342[k];

        t_343[k] = f_0 * hi_343[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, t_349, t_350, t_351, hi_344, \
                         hi_345, hi_346, hi_347, hi_348, hi_349, hi_350, \
                         hi_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_0 * hi_344[k];

        t_345[k] = f_0 * hi_345[k];

        t_346[k] = f_0 * hi_346[k];

        t_347[k] = f_0 * hi_347[k];

        t_348[k] = f_0 * hi_348[k];

        t_349[k] = f_0 * hi_349[k];

        t_350[k] = f_0 * hi_350[k];

        t_351[k] = f_0 * hi_351[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, t_357, t_358, t_359, hi_352, \
                         hi_353, hi_354, hi_355, hi_356, hi_357, hi_358, \
                         hi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_0 * hi_352[k];

        t_353[k] = f_0 * hi_353[k];

        t_354[k] = f_0 * hi_354[k];

        t_355[k] = f_0 * hi_355[k];

        t_356[k] = f_0 * hi_356[k];

        t_357[k] = f_0 * hi_357[k];

        t_358[k] = f_0 * hi_358[k];

        t_359[k] = f_0 * hi_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, t_365, t_366, t_367, hi_360, \
                         hi_361, hi_362, hi_363, hi_364, hi_365, hi_366, \
                         hi_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_0 * hi_360[k];

        t_361[k] = f_0 * hi_361[k];

        t_362[k] = f_0 * hi_362[k];

        t_363[k] = f_0 * hi_363[k];

        t_364[k] = f_0 * hi_364[k];

        t_365[k] = f_0 * hi_365[k];

        t_366[k] = f_0 * hi_366[k];

        t_367[k] = f_0 * hi_367[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, t_373, t_374, t_375, hi_368, \
                         hi_369, hi_370, hi_371, hi_372, hi_373, hi_374, \
                         hi_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_0 * hi_368[k];

        t_369[k] = f_0 * hi_369[k];

        t_370[k] = f_0 * hi_370[k];

        t_371[k] = f_0 * hi_371[k];

        t_372[k] = f_0 * hi_372[k];

        t_373[k] = f_0 * hi_373[k];

        t_374[k] = f_0 * hi_374[k];

        t_375[k] = f_0 * hi_375[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, t_380, t_381, t_382, t_383, hi_376, \
                         hi_377, hi_378, hi_379, hi_380, hi_381, hi_382, \
                         hi_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_0 * hi_376[k];

        t_377[k] = f_0 * hi_377[k];

        t_378[k] = f_0 * hi_378[k];

        t_379[k] = f_0 * hi_379[k];

        t_380[k] = f_0 * hi_380[k];

        t_381[k] = f_0 * hi_381[k];

        t_382[k] = f_0 * hi_382[k];

        t_383[k] = f_0 * hi_383[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, t_388, t_389, t_390, t_391, hi_384, \
                         hi_385, hi_386, hi_387, hi_388, hi_389, hi_390, \
                         hi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_0 * hi_384[k];

        t_385[k] = f_0 * hi_385[k];

        t_386[k] = f_0 * hi_386[k];

        t_387[k] = f_0 * hi_387[k];

        t_388[k] = f_0 * hi_388[k];

        t_389[k] = f_0 * hi_389[k];

        t_390[k] = f_0 * hi_390[k];

        t_391[k] = f_0 * hi_391[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, t_397, t_398, t_399, hi_392, \
                         hi_393, hi_394, hi_395, hi_396, hi_397, hi_398, \
                         hi_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_0 * hi_392[k];

        t_393[k] = f_0 * hi_393[k];

        t_394[k] = f_0 * hi_394[k];

        t_395[k] = f_0 * hi_395[k];

        t_396[k] = f_0 * hi_396[k];

        t_397[k] = f_0 * hi_397[k];

        t_398[k] = f_0 * hi_398[k];

        t_399[k] = f_0 * hi_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, t_405, t_406, t_407, hi_400, \
                         hi_401, hi_402, hi_403, hi_404, hi_405, hi_406, \
                         hi_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_0 * hi_400[k];

        t_401[k] = f_0 * hi_401[k];

        t_402[k] = f_0 * hi_402[k];

        t_403[k] = f_0 * hi_403[k];

        t_404[k] = f_0 * hi_404[k];

        t_405[k] = f_0 * hi_405[k];

        t_406[k] = f_0 * hi_406[k];

        t_407[k] = f_0 * hi_407[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, t_413, t_414, t_415, hi_408, \
                         hi_409, hi_410, hi_411, hi_412, hi_413, hi_414, \
                         hi_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_0 * hi_408[k];

        t_409[k] = f_0 * hi_409[k];

        t_410[k] = f_0 * hi_410[k];

        t_411[k] = f_0 * hi_411[k];

        t_412[k] = f_0 * hi_412[k];

        t_413[k] = f_0 * hi_413[k];

        t_414[k] = f_0 * hi_414[k];

        t_415[k] = f_0 * hi_415[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, hi_416, hi_417, hi_418, \
                         hi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_0 * hi_416[k];

        t_417[k] = f_0 * hi_417[k];

        t_418[k] = f_0 * hi_418[k];

        t_419[k] = f_0 * hi_419[k];
    }
}

auto
compute_prim_geom_10_gi_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t fi, const size_t hi,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gi_electron_repulsion_0_piece0(buffer, target, fi, hi, ncols, alpha);

    compute_prim_geom_10_gi_electron_repulsion_0_piece1(buffer, target, fi, hi, ncols, alpha);

    compute_prim_geom_10_gi_electron_repulsion_0_piece2(buffer, target, hi, ncols, alpha);
}

static auto
compute_prim_geom_10_gi_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fi, const size_t hi,
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

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_41 = buffer.data(fi + 41);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_60 = buffer.data(fi + 60);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_63 = buffer.data(fi + 63);
    const auto *fi_64 = buffer.data(fi + 64);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_66 = buffer.data(fi + 66);
    const auto *fi_67 = buffer.data(fi + 67);
    const auto *fi_68 = buffer.data(fi + 68);
    const auto *fi_69 = buffer.data(fi + 69);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_71 = buffer.data(fi + 71);
    const auto *fi_72 = buffer.data(fi + 72);
    const auto *fi_73 = buffer.data(fi + 73);
    const auto *fi_74 = buffer.data(fi + 74);
    const auto *fi_75 = buffer.data(fi + 75);
    const auto *fi_76 = buffer.data(fi + 76);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_78 = buffer.data(fi + 78);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);
    const auto *fi_81 = buffer.data(fi + 81);
    const auto *fi_82 = buffer.data(fi + 82);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_84 = buffer.data(fi + 84);
    const auto *fi_85 = buffer.data(fi + 85);
    const auto *fi_86 = buffer.data(fi + 86);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_88 = buffer.data(fi + 88);
    const auto *fi_89 = buffer.data(fi + 89);
    const auto *fi_90 = buffer.data(fi + 90);

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
    const auto *hi_280 = buffer.data(hi + 280);
    const auto *hi_281 = buffer.data(hi + 281);
    const auto *hi_282 = buffer.data(hi + 282);
    const auto *hi_283 = buffer.data(hi + 283);
    const auto *hi_284 = buffer.data(hi + 284);
    const auto *hi_285 = buffer.data(hi + 285);
    const auto *hi_286 = buffer.data(hi + 286);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hi_28, hi_29, hi_30, hi_31, \
                         hi_32, hi_33, hi_34, hi_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hi_28[k];

        t_1[k] = f_0 * hi_29[k];

        t_2[k] = f_0 * hi_30[k];

        t_3[k] = f_0 * hi_31[k];

        t_4[k] = f_0 * hi_32[k];

        t_5[k] = f_0 * hi_33[k];

        t_6[k] = f_0 * hi_34[k];

        t_7[k] = f_0 * hi_35[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, hi_36, hi_37, hi_38, \
                         hi_39, hi_40, hi_41, hi_42, hi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hi_36[k];

        t_9[k] = f_0 * hi_37[k];

        t_10[k] = f_0 * hi_38[k];

        t_11[k] = f_0 * hi_39[k];

        t_12[k] = f_0 * hi_40[k];

        t_13[k] = f_0 * hi_41[k];

        t_14[k] = f_0 * hi_42[k];

        t_15[k] = f_0 * hi_43[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, hi_44, hi_45, hi_46, \
                         hi_47, hi_48, hi_49, hi_50, hi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * hi_44[k];

        t_17[k] = f_0 * hi_45[k];

        t_18[k] = f_0 * hi_46[k];

        t_19[k] = f_0 * hi_47[k];

        t_20[k] = f_0 * hi_48[k];

        t_21[k] = f_0 * hi_49[k];

        t_22[k] = f_0 * hi_50[k];

        t_23[k] = f_0 * hi_51[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, fi_0, fi_1, hi_52, hi_53, hi_54, \
                         hi_55, hi_84, hi_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * hi_52[k];

        t_25[k] = f_0 * hi_53[k];

        t_26[k] = f_0 * hi_54[k];

        t_27[k] = f_0 * hi_55[k];

        t_28[k] = -fi_0[k]
                  + f_0 * hi_84[k];

        t_29[k] = -fi_1[k]
                  + f_0 * hi_85[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, fi_2, fi_3, fi_4, fi_5, fi_6, hi_86, \
                         hi_87, hi_88, hi_89, hi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -fi_2[k]
                  + f_0 * hi_86[k];

        t_31[k] = -fi_3[k]
                  + f_0 * hi_87[k];

        t_32[k] = -fi_4[k]
                  + f_0 * hi_88[k];

        t_33[k] = -fi_5[k]
                  + f_0 * hi_89[k];

        t_34[k] = -fi_6[k]
                  + f_0 * hi_90[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, fi_7, fi_8, fi_9, fi_10, fi_11, hi_91, \
                         hi_92, hi_93, hi_94, hi_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -fi_7[k]
                  + f_0 * hi_91[k];

        t_36[k] = -fi_8[k]
                  + f_0 * hi_92[k];

        t_37[k] = -fi_9[k]
                  + f_0 * hi_93[k];

        t_38[k] = -fi_10[k]
                  + f_0 * hi_94[k];

        t_39[k] = -fi_11[k]
                  + f_0 * hi_95[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, fi_12, fi_13, fi_14, fi_15, fi_16, \
                         hi_96, hi_97, hi_98, hi_99, hi_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -fi_12[k]
                  + f_0 * hi_96[k];

        t_41[k] = -fi_13[k]
                  + f_0 * hi_97[k];

        t_42[k] = -fi_14[k]
                  + f_0 * hi_98[k];

        t_43[k] = -fi_15[k]
                  + f_0 * hi_99[k];

        t_44[k] = -fi_16[k]
                  + f_0 * hi_100[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, fi_17, fi_18, fi_19, fi_20, fi_21, \
                         hi_101, hi_102, hi_103, hi_104, hi_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -fi_17[k]
                  + f_0 * hi_101[k];

        t_46[k] = -fi_18[k]
                  + f_0 * hi_102[k];

        t_47[k] = -fi_19[k]
                  + f_0 * hi_103[k];

        t_48[k] = -fi_20[k]
                  + f_0 * hi_104[k];

        t_49[k] = -fi_21[k]
                  + f_0 * hi_105[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, fi_22, fi_23, fi_24, fi_25, fi_26, \
                         hi_106, hi_107, hi_108, hi_109, hi_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -fi_22[k]
                  + f_0 * hi_106[k];

        t_51[k] = -fi_23[k]
                  + f_0 * hi_107[k];

        t_52[k] = -fi_24[k]
                  + f_0 * hi_108[k];

        t_53[k] = -fi_25[k]
                  + f_0 * hi_109[k];

        t_54[k] = -fi_26[k]
                  + f_0 * hi_110[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, fi_27, hi_111, hi_112, \
                         hi_113, hi_114, hi_115, hi_116, hi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -fi_27[k]
                  + f_0 * hi_111[k];

        t_56[k] = f_0 * hi_112[k];

        t_57[k] = f_0 * hi_113[k];

        t_58[k] = f_0 * hi_114[k];

        t_59[k] = f_0 * hi_115[k];

        t_60[k] = f_0 * hi_116[k];

        t_61[k] = f_0 * hi_117[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, t_67, t_68, t_69, hi_118, hi_119, \
                         hi_120, hi_121, hi_122, hi_123, hi_124, \
                         hi_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * hi_118[k];

        t_63[k] = f_0 * hi_119[k];

        t_64[k] = f_0 * hi_120[k];

        t_65[k] = f_0 * hi_121[k];

        t_66[k] = f_0 * hi_122[k];

        t_67[k] = f_0 * hi_123[k];

        t_68[k] = f_0 * hi_124[k];

        t_69[k] = f_0 * hi_125[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, t_75, t_76, t_77, hi_126, hi_127, \
                         hi_128, hi_129, hi_130, hi_131, hi_132, \
                         hi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_0 * hi_126[k];

        t_71[k] = f_0 * hi_127[k];

        t_72[k] = f_0 * hi_128[k];

        t_73[k] = f_0 * hi_129[k];

        t_74[k] = f_0 * hi_130[k];

        t_75[k] = f_0 * hi_131[k];

        t_76[k] = f_0 * hi_132[k];

        t_77[k] = f_0 * hi_133[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, t_84, fi_28, hi_134, hi_135, \
                         hi_136, hi_137, hi_138, hi_139, hi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * hi_134[k];

        t_79[k] = f_0 * hi_135[k];

        t_80[k] = f_0 * hi_136[k];

        t_81[k] = f_0 * hi_137[k];

        t_82[k] = f_0 * hi_138[k];

        t_83[k] = f_0 * hi_139[k];

        t_84[k] = -2.0 * fi_28[k]
                  + f_0 * hi_168[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, fi_29, fi_30, fi_31, fi_32, fi_33, \
                         hi_169, hi_170, hi_171, hi_172, hi_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * fi_29[k]
                  + f_0 * hi_169[k];

        t_86[k] = -2.0 * fi_30[k]
                  + f_0 * hi_170[k];

        t_87[k] = -2.0 * fi_31[k]
                  + f_0 * hi_171[k];

        t_88[k] = -2.0 * fi_32[k]
                  + f_0 * hi_172[k];

        t_89[k] = -2.0 * fi_33[k]
                  + f_0 * hi_173[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, fi_34, fi_35, fi_36, fi_37, fi_38, \
                         hi_174, hi_175, hi_176, hi_177, hi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * fi_34[k]
                  + f_0 * hi_174[k];

        t_91[k] = -2.0 * fi_35[k]
                  + f_0 * hi_175[k];

        t_92[k] = -2.0 * fi_36[k]
                  + f_0 * hi_176[k];

        t_93[k] = -2.0 * fi_37[k]
                  + f_0 * hi_177[k];

        t_94[k] = -2.0 * fi_38[k]
                  + f_0 * hi_178[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, fi_39, fi_40, fi_41, fi_42, fi_43, \
                         hi_179, hi_180, hi_181, hi_182, hi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * fi_39[k]
                  + f_0 * hi_179[k];

        t_96[k] = -2.0 * fi_40[k]
                  + f_0 * hi_180[k];

        t_97[k] = -2.0 * fi_41[k]
                  + f_0 * hi_181[k];

        t_98[k] = -2.0 * fi_42[k]
                  + f_0 * hi_182[k];

        t_99[k] = -2.0 * fi_43[k]
                  + f_0 * hi_183[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, fi_44, fi_45, fi_46, fi_47, fi_48, \
                         hi_184, hi_185, hi_186, hi_187, hi_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * fi_44[k]
                   + f_0 * hi_184[k];

        t_101[k] = -2.0 * fi_45[k]
                   + f_0 * hi_185[k];

        t_102[k] = -2.0 * fi_46[k]
                   + f_0 * hi_186[k];

        t_103[k] = -2.0 * fi_47[k]
                   + f_0 * hi_187[k];

        t_104[k] = -2.0 * fi_48[k]
                   + f_0 * hi_188[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, fi_49, fi_50, fi_51, fi_52, fi_53, \
                         hi_189, hi_190, hi_191, hi_192, hi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * fi_49[k]
                   + f_0 * hi_189[k];

        t_106[k] = -2.0 * fi_50[k]
                   + f_0 * hi_190[k];

        t_107[k] = -2.0 * fi_51[k]
                   + f_0 * hi_191[k];

        t_108[k] = -2.0 * fi_52[k]
                   + f_0 * hi_192[k];

        t_109[k] = -2.0 * fi_53[k]
                   + f_0 * hi_193[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, fi_54, fi_55, fi_56, fi_57, fi_58, \
                         hi_194, hi_195, hi_196, hi_197, hi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * fi_54[k]
                   + f_0 * hi_194[k];

        t_111[k] = -2.0 * fi_55[k]
                   + f_0 * hi_195[k];

        t_112[k] = -fi_56[k]
                   + f_0 * hi_196[k];

        t_113[k] = -fi_57[k]
                   + f_0 * hi_197[k];

        t_114[k] = -fi_58[k]
                   + f_0 * hi_198[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, fi_59, fi_60, fi_61, fi_62, fi_63, \
                         hi_199, hi_200, hi_201, hi_202, hi_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -fi_59[k]
                   + f_0 * hi_199[k];

        t_116[k] = -fi_60[k]
                   + f_0 * hi_200[k];

        t_117[k] = -fi_61[k]
                   + f_0 * hi_201[k];

        t_118[k] = -fi_62[k]
                   + f_0 * hi_202[k];

        t_119[k] = -fi_63[k]
                   + f_0 * hi_203[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, fi_64, fi_65, fi_66, fi_67, fi_68, \
                         hi_204, hi_205, hi_206, hi_207, hi_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -fi_64[k]
                   + f_0 * hi_204[k];

        t_121[k] = -fi_65[k]
                   + f_0 * hi_205[k];

        t_122[k] = -fi_66[k]
                   + f_0 * hi_206[k];

        t_123[k] = -fi_67[k]
                   + f_0 * hi_207[k];

        t_124[k] = -fi_68[k]
                   + f_0 * hi_208[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, fi_69, fi_70, fi_71, fi_72, fi_73, \
                         hi_209, hi_210, hi_211, hi_212, hi_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -fi_69[k]
                   + f_0 * hi_209[k];

        t_126[k] = -fi_70[k]
                   + f_0 * hi_210[k];

        t_127[k] = -fi_71[k]
                   + f_0 * hi_211[k];

        t_128[k] = -fi_72[k]
                   + f_0 * hi_212[k];

        t_129[k] = -fi_73[k]
                   + f_0 * hi_213[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, fi_74, fi_75, fi_76, fi_77, fi_78, \
                         hi_214, hi_215, hi_216, hi_217, hi_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -fi_74[k]
                   + f_0 * hi_214[k];

        t_131[k] = -fi_75[k]
                   + f_0 * hi_215[k];

        t_132[k] = -fi_76[k]
                   + f_0 * hi_216[k];

        t_133[k] = -fi_77[k]
                   + f_0 * hi_217[k];

        t_134[k] = -fi_78[k]
                   + f_0 * hi_218[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, fi_79, fi_80, fi_81, fi_82, fi_83, \
                         hi_219, hi_220, hi_221, hi_222, hi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -fi_79[k]
                   + f_0 * hi_219[k];

        t_136[k] = -fi_80[k]
                   + f_0 * hi_220[k];

        t_137[k] = -fi_81[k]
                   + f_0 * hi_221[k];

        t_138[k] = -fi_82[k]
                   + f_0 * hi_222[k];

        t_139[k] = -fi_83[k]
                   + f_0 * hi_223[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, t_146, t_147, hi_224, \
                         hi_225, hi_226, hi_227, hi_228, hi_229, hi_230, \
                         hi_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * hi_224[k];

        t_141[k] = f_0 * hi_225[k];

        t_142[k] = f_0 * hi_226[k];

        t_143[k] = f_0 * hi_227[k];

        t_144[k] = f_0 * hi_228[k];

        t_145[k] = f_0 * hi_229[k];

        t_146[k] = f_0 * hi_230[k];

        t_147[k] = f_0 * hi_231[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, t_153, t_154, t_155, hi_232, \
                         hi_233, hi_234, hi_235, hi_236, hi_237, hi_238, \
                         hi_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * hi_232[k];

        t_149[k] = f_0 * hi_233[k];

        t_150[k] = f_0 * hi_234[k];

        t_151[k] = f_0 * hi_235[k];

        t_152[k] = f_0 * hi_236[k];

        t_153[k] = f_0 * hi_237[k];

        t_154[k] = f_0 * hi_238[k];

        t_155[k] = f_0 * hi_239[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, t_162, t_163, hi_240, \
                         hi_241, hi_242, hi_243, hi_244, hi_245, hi_246, \
                         hi_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * hi_240[k];

        t_157[k] = f_0 * hi_241[k];

        t_158[k] = f_0 * hi_242[k];

        t_159[k] = f_0 * hi_243[k];

        t_160[k] = f_0 * hi_244[k];

        t_161[k] = f_0 * hi_245[k];

        t_162[k] = f_0 * hi_246[k];

        t_163[k] = f_0 * hi_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, t_169, fi_84, fi_85, hi_248, \
                         hi_249, hi_250, hi_251, hi_280, hi_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_0 * hi_248[k];

        t_165[k] = f_0 * hi_249[k];

        t_166[k] = f_0 * hi_250[k];

        t_167[k] = f_0 * hi_251[k];

        t_168[k] = -3.0 * fi_84[k]
                   + f_0 * hi_280[k];

        t_169[k] = -3.0 * fi_85[k]
                   + f_0 * hi_281[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, fi_86, fi_87, fi_88, fi_89, fi_90, \
                         hi_282, hi_283, hi_284, hi_285, hi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * fi_86[k]
                   + f_0 * hi_282[k];

        t_171[k] = -3.0 * fi_87[k]
                   + f_0 * hi_283[k];

        t_172[k] = -3.0 * fi_88[k]
                   + f_0 * hi_284[k];

        t_173[k] = -3.0 * fi_89[k]
                   + f_0 * hi_285[k];

        t_174[k] = -3.0 * fi_90[k]
                   + f_0 * hi_286[k];
    }
}

static auto
compute_prim_geom_10_gi_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fi, const size_t hi,
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

    const auto *fi_91 = buffer.data(fi + 91);
    const auto *fi_92 = buffer.data(fi + 92);
    const auto *fi_93 = buffer.data(fi + 93);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_95 = buffer.data(fi + 95);
    const auto *fi_96 = buffer.data(fi + 96);
    const auto *fi_97 = buffer.data(fi + 97);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_100 = buffer.data(fi + 100);
    const auto *fi_101 = buffer.data(fi + 101);
    const auto *fi_102 = buffer.data(fi + 102);
    const auto *fi_103 = buffer.data(fi + 103);
    const auto *fi_104 = buffer.data(fi + 104);
    const auto *fi_105 = buffer.data(fi + 105);
    const auto *fi_106 = buffer.data(fi + 106);
    const auto *fi_107 = buffer.data(fi + 107);
    const auto *fi_108 = buffer.data(fi + 108);
    const auto *fi_109 = buffer.data(fi + 109);
    const auto *fi_110 = buffer.data(fi + 110);
    const auto *fi_111 = buffer.data(fi + 111);
    const auto *fi_112 = buffer.data(fi + 112);
    const auto *fi_113 = buffer.data(fi + 113);
    const auto *fi_114 = buffer.data(fi + 114);
    const auto *fi_115 = buffer.data(fi + 115);
    const auto *fi_116 = buffer.data(fi + 116);
    const auto *fi_117 = buffer.data(fi + 117);
    const auto *fi_118 = buffer.data(fi + 118);
    const auto *fi_119 = buffer.data(fi + 119);
    const auto *fi_120 = buffer.data(fi + 120);
    const auto *fi_121 = buffer.data(fi + 121);
    const auto *fi_122 = buffer.data(fi + 122);
    const auto *fi_123 = buffer.data(fi + 123);
    const auto *fi_124 = buffer.data(fi + 124);
    const auto *fi_125 = buffer.data(fi + 125);
    const auto *fi_126 = buffer.data(fi + 126);
    const auto *fi_127 = buffer.data(fi + 127);
    const auto *fi_128 = buffer.data(fi + 128);
    const auto *fi_129 = buffer.data(fi + 129);
    const auto *fi_130 = buffer.data(fi + 130);
    const auto *fi_131 = buffer.data(fi + 131);
    const auto *fi_132 = buffer.data(fi + 132);
    const auto *fi_133 = buffer.data(fi + 133);
    const auto *fi_134 = buffer.data(fi + 134);
    const auto *fi_135 = buffer.data(fi + 135);
    const auto *fi_136 = buffer.data(fi + 136);
    const auto *fi_137 = buffer.data(fi + 137);
    const auto *fi_138 = buffer.data(fi + 138);
    const auto *fi_139 = buffer.data(fi + 139);
    const auto *fi_140 = buffer.data(fi + 140);
    const auto *fi_141 = buffer.data(fi + 141);
    const auto *fi_142 = buffer.data(fi + 142);
    const auto *fi_143 = buffer.data(fi + 143);
    const auto *fi_144 = buffer.data(fi + 144);
    const auto *fi_145 = buffer.data(fi + 145);
    const auto *fi_146 = buffer.data(fi + 146);
    const auto *fi_147 = buffer.data(fi + 147);
    const auto *fi_148 = buffer.data(fi + 148);
    const auto *fi_149 = buffer.data(fi + 149);
    const auto *fi_150 = buffer.data(fi + 150);
    const auto *fi_151 = buffer.data(fi + 151);
    const auto *fi_152 = buffer.data(fi + 152);
    const auto *fi_153 = buffer.data(fi + 153);
    const auto *fi_154 = buffer.data(fi + 154);
    const auto *fi_155 = buffer.data(fi + 155);
    const auto *fi_156 = buffer.data(fi + 156);
    const auto *fi_157 = buffer.data(fi + 157);
    const auto *fi_158 = buffer.data(fi + 158);
    const auto *fi_159 = buffer.data(fi + 159);
    const auto *fi_160 = buffer.data(fi + 160);
    const auto *fi_161 = buffer.data(fi + 161);
    const auto *fi_162 = buffer.data(fi + 162);
    const auto *fi_163 = buffer.data(fi + 163);
    const auto *fi_164 = buffer.data(fi + 164);
    const auto *fi_165 = buffer.data(fi + 165);
    const auto *fi_166 = buffer.data(fi + 166);
    const auto *fi_167 = buffer.data(fi + 167);
    const auto *fi_168 = buffer.data(fi + 168);
    const auto *fi_169 = buffer.data(fi + 169);
    const auto *fi_170 = buffer.data(fi + 170);
    const auto *fi_171 = buffer.data(fi + 171);
    const auto *fi_172 = buffer.data(fi + 172);
    const auto *fi_173 = buffer.data(fi + 173);
    const auto *fi_174 = buffer.data(fi + 174);
    const auto *fi_175 = buffer.data(fi + 175);
    const auto *fi_176 = buffer.data(fi + 176);
    const auto *fi_177 = buffer.data(fi + 177);
    const auto *fi_178 = buffer.data(fi + 178);
    const auto *fi_179 = buffer.data(fi + 179);
    const auto *fi_180 = buffer.data(fi + 180);
    const auto *fi_181 = buffer.data(fi + 181);
    const auto *fi_182 = buffer.data(fi + 182);
    const auto *fi_183 = buffer.data(fi + 183);
    const auto *fi_184 = buffer.data(fi + 184);
    const auto *fi_185 = buffer.data(fi + 185);
    const auto *fi_186 = buffer.data(fi + 186);
    const auto *fi_187 = buffer.data(fi + 187);
    const auto *fi_188 = buffer.data(fi + 188);
    const auto *fi_189 = buffer.data(fi + 189);
    const auto *fi_190 = buffer.data(fi + 190);
    const auto *fi_191 = buffer.data(fi + 191);
    const auto *fi_192 = buffer.data(fi + 192);
    const auto *fi_193 = buffer.data(fi + 193);
    const auto *fi_194 = buffer.data(fi + 194);
    const auto *fi_195 = buffer.data(fi + 195);
    const auto *fi_196 = buffer.data(fi + 196);
    const auto *fi_197 = buffer.data(fi + 197);
    const auto *fi_198 = buffer.data(fi + 198);
    const auto *fi_199 = buffer.data(fi + 199);
    const auto *fi_200 = buffer.data(fi + 200);
    const auto *fi_201 = buffer.data(fi + 201);
    const auto *fi_202 = buffer.data(fi + 202);
    const auto *fi_203 = buffer.data(fi + 203);
    const auto *fi_204 = buffer.data(fi + 204);
    const auto *fi_205 = buffer.data(fi + 205);
    const auto *fi_206 = buffer.data(fi + 206);
    const auto *fi_207 = buffer.data(fi + 207);
    const auto *fi_208 = buffer.data(fi + 208);
    const auto *fi_209 = buffer.data(fi + 209);
    const auto *fi_210 = buffer.data(fi + 210);
    const auto *fi_211 = buffer.data(fi + 211);
    const auto *fi_212 = buffer.data(fi + 212);
    const auto *fi_213 = buffer.data(fi + 213);
    const auto *fi_214 = buffer.data(fi + 214);
    const auto *fi_215 = buffer.data(fi + 215);
    const auto *fi_216 = buffer.data(fi + 216);
    const auto *fi_217 = buffer.data(fi + 217);
    const auto *fi_218 = buffer.data(fi + 218);
    const auto *fi_219 = buffer.data(fi + 219);
    const auto *fi_220 = buffer.data(fi + 220);
    const auto *fi_221 = buffer.data(fi + 221);
    const auto *fi_222 = buffer.data(fi + 222);

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

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, fi_91, fi_92, fi_93, fi_94, fi_95, \
                         hi_287, hi_288, hi_289, hi_290, hi_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * fi_91[k]
                   + f_0 * hi_287[k];

        t_176[k] = -3.0 * fi_92[k]
                   + f_0 * hi_288[k];

        t_177[k] = -3.0 * fi_93[k]
                   + f_0 * hi_289[k];

        t_178[k] = -3.0 * fi_94[k]
                   + f_0 * hi_290[k];

        t_179[k] = -3.0 * fi_95[k]
                   + f_0 * hi_291[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, fi_96, fi_97, fi_98, fi_99, \
                         fi_100, hi_292, hi_293, hi_294, hi_295, \
                         hi_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * fi_96[k]
                   + f_0 * hi_292[k];

        t_181[k] = -3.0 * fi_97[k]
                   + f_0 * hi_293[k];

        t_182[k] = -3.0 * fi_98[k]
                   + f_0 * hi_294[k];

        t_183[k] = -3.0 * fi_99[k]
                   + f_0 * hi_295[k];

        t_184[k] = -3.0 * fi_100[k]
                   + f_0 * hi_296[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, fi_101, fi_102, fi_103, fi_104, \
                         fi_105, hi_297, hi_298, hi_299, hi_300, \
                         hi_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * fi_101[k]
                   + f_0 * hi_297[k];

        t_186[k] = -3.0 * fi_102[k]
                   + f_0 * hi_298[k];

        t_187[k] = -3.0 * fi_103[k]
                   + f_0 * hi_299[k];

        t_188[k] = -3.0 * fi_104[k]
                   + f_0 * hi_300[k];

        t_189[k] = -3.0 * fi_105[k]
                   + f_0 * hi_301[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, fi_106, fi_107, fi_108, fi_109, \
                         fi_110, hi_302, hi_303, hi_304, hi_305, \
                         hi_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * fi_106[k]
                   + f_0 * hi_302[k];

        t_191[k] = -3.0 * fi_107[k]
                   + f_0 * hi_303[k];

        t_192[k] = -3.0 * fi_108[k]
                   + f_0 * hi_304[k];

        t_193[k] = -3.0 * fi_109[k]
                   + f_0 * hi_305[k];

        t_194[k] = -3.0 * fi_110[k]
                   + f_0 * hi_306[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, fi_111, fi_112, fi_113, fi_114, \
                         fi_115, hi_307, hi_308, hi_309, hi_310, \
                         hi_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * fi_111[k]
                   + f_0 * hi_307[k];

        t_196[k] = -2.0 * fi_112[k]
                   + f_0 * hi_308[k];

        t_197[k] = -2.0 * fi_113[k]
                   + f_0 * hi_309[k];

        t_198[k] = -2.0 * fi_114[k]
                   + f_0 * hi_310[k];

        t_199[k] = -2.0 * fi_115[k]
                   + f_0 * hi_311[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, fi_116, fi_117, fi_118, fi_119, \
                         fi_120, hi_312, hi_313, hi_314, hi_315, \
                         hi_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * fi_116[k]
                   + f_0 * hi_312[k];

        t_201[k] = -2.0 * fi_117[k]
                   + f_0 * hi_313[k];

        t_202[k] = -2.0 * fi_118[k]
                   + f_0 * hi_314[k];

        t_203[k] = -2.0 * fi_119[k]
                   + f_0 * hi_315[k];

        t_204[k] = -2.0 * fi_120[k]
                   + f_0 * hi_316[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, fi_121, fi_122, fi_123, fi_124, \
                         fi_125, hi_317, hi_318, hi_319, hi_320, \
                         hi_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * fi_121[k]
                   + f_0 * hi_317[k];

        t_206[k] = -2.0 * fi_122[k]
                   + f_0 * hi_318[k];

        t_207[k] = -2.0 * fi_123[k]
                   + f_0 * hi_319[k];

        t_208[k] = -2.0 * fi_124[k]
                   + f_0 * hi_320[k];

        t_209[k] = -2.0 * fi_125[k]
                   + f_0 * hi_321[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, fi_126, fi_127, fi_128, fi_129, \
                         fi_130, hi_322, hi_323, hi_324, hi_325, \
                         hi_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * fi_126[k]
                   + f_0 * hi_322[k];

        t_211[k] = -2.0 * fi_127[k]
                   + f_0 * hi_323[k];

        t_212[k] = -2.0 * fi_128[k]
                   + f_0 * hi_324[k];

        t_213[k] = -2.0 * fi_129[k]
                   + f_0 * hi_325[k];

        t_214[k] = -2.0 * fi_130[k]
                   + f_0 * hi_326[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, fi_131, fi_132, fi_133, fi_134, \
                         fi_135, hi_327, hi_328, hi_329, hi_330, \
                         hi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * fi_131[k]
                   + f_0 * hi_327[k];

        t_216[k] = -2.0 * fi_132[k]
                   + f_0 * hi_328[k];

        t_217[k] = -2.0 * fi_133[k]
                   + f_0 * hi_329[k];

        t_218[k] = -2.0 * fi_134[k]
                   + f_0 * hi_330[k];

        t_219[k] = -2.0 * fi_135[k]
                   + f_0 * hi_331[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, fi_136, fi_137, fi_138, fi_139, \
                         fi_140, hi_332, hi_333, hi_334, hi_335, \
                         hi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * fi_136[k]
                   + f_0 * hi_332[k];

        t_221[k] = -2.0 * fi_137[k]
                   + f_0 * hi_333[k];

        t_222[k] = -2.0 * fi_138[k]
                   + f_0 * hi_334[k];

        t_223[k] = -2.0 * fi_139[k]
                   + f_0 * hi_335[k];

        t_224[k] = -fi_140[k]
                   + f_0 * hi_336[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, fi_141, fi_142, fi_143, fi_144, \
                         fi_145, hi_337, hi_338, hi_339, hi_340, \
                         hi_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -fi_141[k]
                   + f_0 * hi_337[k];

        t_226[k] = -fi_142[k]
                   + f_0 * hi_338[k];

        t_227[k] = -fi_143[k]
                   + f_0 * hi_339[k];

        t_228[k] = -fi_144[k]
                   + f_0 * hi_340[k];

        t_229[k] = -fi_145[k]
                   + f_0 * hi_341[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, fi_146, fi_147, fi_148, fi_149, \
                         fi_150, hi_342, hi_343, hi_344, hi_345, \
                         hi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -fi_146[k]
                   + f_0 * hi_342[k];

        t_231[k] = -fi_147[k]
                   + f_0 * hi_343[k];

        t_232[k] = -fi_148[k]
                   + f_0 * hi_344[k];

        t_233[k] = -fi_149[k]
                   + f_0 * hi_345[k];

        t_234[k] = -fi_150[k]
                   + f_0 * hi_346[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, fi_151, fi_152, fi_153, fi_154, \
                         fi_155, hi_347, hi_348, hi_349, hi_350, \
                         hi_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -fi_151[k]
                   + f_0 * hi_347[k];

        t_236[k] = -fi_152[k]
                   + f_0 * hi_348[k];

        t_237[k] = -fi_153[k]
                   + f_0 * hi_349[k];

        t_238[k] = -fi_154[k]
                   + f_0 * hi_350[k];

        t_239[k] = -fi_155[k]
                   + f_0 * hi_351[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, fi_156, fi_157, fi_158, fi_159, \
                         fi_160, hi_352, hi_353, hi_354, hi_355, \
                         hi_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -fi_156[k]
                   + f_0 * hi_352[k];

        t_241[k] = -fi_157[k]
                   + f_0 * hi_353[k];

        t_242[k] = -fi_158[k]
                   + f_0 * hi_354[k];

        t_243[k] = -fi_159[k]
                   + f_0 * hi_355[k];

        t_244[k] = -fi_160[k]
                   + f_0 * hi_356[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, fi_161, fi_162, fi_163, fi_164, \
                         fi_165, hi_357, hi_358, hi_359, hi_360, \
                         hi_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -fi_161[k]
                   + f_0 * hi_357[k];

        t_246[k] = -fi_162[k]
                   + f_0 * hi_358[k];

        t_247[k] = -fi_163[k]
                   + f_0 * hi_359[k];

        t_248[k] = -fi_164[k]
                   + f_0 * hi_360[k];

        t_249[k] = -fi_165[k]
                   + f_0 * hi_361[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, t_255, t_256, fi_166, fi_167, \
                         hi_362, hi_363, hi_364, hi_365, hi_366, hi_367, \
                         hi_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -fi_166[k]
                   + f_0 * hi_362[k];

        t_251[k] = -fi_167[k]
                   + f_0 * hi_363[k];

        t_252[k] = f_0 * hi_364[k];

        t_253[k] = f_0 * hi_365[k];

        t_254[k] = f_0 * hi_366[k];

        t_255[k] = f_0 * hi_367[k];

        t_256[k] = f_0 * hi_368[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, t_262, t_263, t_264, hi_369, \
                         hi_370, hi_371, hi_372, hi_373, hi_374, hi_375, \
                         hi_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_0 * hi_369[k];

        t_258[k] = f_0 * hi_370[k];

        t_259[k] = f_0 * hi_371[k];

        t_260[k] = f_0 * hi_372[k];

        t_261[k] = f_0 * hi_373[k];

        t_262[k] = f_0 * hi_374[k];

        t_263[k] = f_0 * hi_375[k];

        t_264[k] = f_0 * hi_376[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, t_270, t_271, t_272, hi_377, \
                         hi_378, hi_379, hi_380, hi_381, hi_382, hi_383, \
                         hi_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_0 * hi_377[k];

        t_266[k] = f_0 * hi_378[k];

        t_267[k] = f_0 * hi_379[k];

        t_268[k] = f_0 * hi_380[k];

        t_269[k] = f_0 * hi_381[k];

        t_270[k] = f_0 * hi_382[k];

        t_271[k] = f_0 * hi_383[k];

        t_272[k] = f_0 * hi_384[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, t_278, t_279, hi_385, hi_386, \
                         hi_387, hi_388, hi_389, hi_390, hi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_0 * hi_385[k];

        t_274[k] = f_0 * hi_386[k];

        t_275[k] = f_0 * hi_387[k];

        t_276[k] = f_0 * hi_388[k];

        t_277[k] = f_0 * hi_389[k];

        t_278[k] = f_0 * hi_390[k];

        t_279[k] = f_0 * hi_391[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, fi_168, fi_169, fi_170, fi_171, \
                         fi_172, hi_420, hi_421, hi_422, hi_423, \
                         hi_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -4.0 * fi_168[k]
                   + f_0 * hi_420[k];

        t_281[k] = -4.0 * fi_169[k]
                   + f_0 * hi_421[k];

        t_282[k] = -4.0 * fi_170[k]
                   + f_0 * hi_422[k];

        t_283[k] = -4.0 * fi_171[k]
                   + f_0 * hi_423[k];

        t_284[k] = -4.0 * fi_172[k]
                   + f_0 * hi_424[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, fi_173, fi_174, fi_175, fi_176, \
                         fi_177, hi_425, hi_426, hi_427, hi_428, \
                         hi_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -4.0 * fi_173[k]
                   + f_0 * hi_425[k];

        t_286[k] = -4.0 * fi_174[k]
                   + f_0 * hi_426[k];

        t_287[k] = -4.0 * fi_175[k]
                   + f_0 * hi_427[k];

        t_288[k] = -4.0 * fi_176[k]
                   + f_0 * hi_428[k];

        t_289[k] = -4.0 * fi_177[k]
                   + f_0 * hi_429[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, fi_178, fi_179, fi_180, fi_181, \
                         fi_182, hi_430, hi_431, hi_432, hi_433, \
                         hi_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -4.0 * fi_178[k]
                   + f_0 * hi_430[k];

        t_291[k] = -4.0 * fi_179[k]
                   + f_0 * hi_431[k];

        t_292[k] = -4.0 * fi_180[k]
                   + f_0 * hi_432[k];

        t_293[k] = -4.0 * fi_181[k]
                   + f_0 * hi_433[k];

        t_294[k] = -4.0 * fi_182[k]
                   + f_0 * hi_434[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, fi_183, fi_184, fi_185, fi_186, \
                         fi_187, hi_435, hi_436, hi_437, hi_438, \
                         hi_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * fi_183[k]
                   + f_0 * hi_435[k];

        t_296[k] = -4.0 * fi_184[k]
                   + f_0 * hi_436[k];

        t_297[k] = -4.0 * fi_185[k]
                   + f_0 * hi_437[k];

        t_298[k] = -4.0 * fi_186[k]
                   + f_0 * hi_438[k];

        t_299[k] = -4.0 * fi_187[k]
                   + f_0 * hi_439[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, fi_188, fi_189, fi_190, fi_191, \
                         fi_192, hi_440, hi_441, hi_442, hi_443, \
                         hi_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * fi_188[k]
                   + f_0 * hi_440[k];

        t_301[k] = -4.0 * fi_189[k]
                   + f_0 * hi_441[k];

        t_302[k] = -4.0 * fi_190[k]
                   + f_0 * hi_442[k];

        t_303[k] = -4.0 * fi_191[k]
                   + f_0 * hi_443[k];

        t_304[k] = -4.0 * fi_192[k]
                   + f_0 * hi_444[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, fi_193, fi_194, fi_195, fi_196, \
                         fi_197, hi_445, hi_446, hi_447, hi_448, \
                         hi_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * fi_193[k]
                   + f_0 * hi_445[k];

        t_306[k] = -4.0 * fi_194[k]
                   + f_0 * hi_446[k];

        t_307[k] = -4.0 * fi_195[k]
                   + f_0 * hi_447[k];

        t_308[k] = -3.0 * fi_196[k]
                   + f_0 * hi_448[k];

        t_309[k] = -3.0 * fi_197[k]
                   + f_0 * hi_449[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, fi_198, fi_199, fi_200, fi_201, \
                         fi_202, hi_450, hi_451, hi_452, hi_453, \
                         hi_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * fi_198[k]
                   + f_0 * hi_450[k];

        t_311[k] = -3.0 * fi_199[k]
                   + f_0 * hi_451[k];

        t_312[k] = -3.0 * fi_200[k]
                   + f_0 * hi_452[k];

        t_313[k] = -3.0 * fi_201[k]
                   + f_0 * hi_453[k];

        t_314[k] = -3.0 * fi_202[k]
                   + f_0 * hi_454[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, fi_203, fi_204, fi_205, fi_206, \
                         fi_207, hi_455, hi_456, hi_457, hi_458, \
                         hi_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -3.0 * fi_203[k]
                   + f_0 * hi_455[k];

        t_316[k] = -3.0 * fi_204[k]
                   + f_0 * hi_456[k];

        t_317[k] = -3.0 * fi_205[k]
                   + f_0 * hi_457[k];

        t_318[k] = -3.0 * fi_206[k]
                   + f_0 * hi_458[k];

        t_319[k] = -3.0 * fi_207[k]
                   + f_0 * hi_459[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, fi_208, fi_209, fi_210, fi_211, \
                         fi_212, hi_460, hi_461, hi_462, hi_463, \
                         hi_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -3.0 * fi_208[k]
                   + f_0 * hi_460[k];

        t_321[k] = -3.0 * fi_209[k]
                   + f_0 * hi_461[k];

        t_322[k] = -3.0 * fi_210[k]
                   + f_0 * hi_462[k];

        t_323[k] = -3.0 * fi_211[k]
                   + f_0 * hi_463[k];

        t_324[k] = -3.0 * fi_212[k]
                   + f_0 * hi_464[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, fi_213, fi_214, fi_215, fi_216, \
                         fi_217, hi_465, hi_466, hi_467, hi_468, \
                         hi_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -3.0 * fi_213[k]
                   + f_0 * hi_465[k];

        t_326[k] = -3.0 * fi_214[k]
                   + f_0 * hi_466[k];

        t_327[k] = -3.0 * fi_215[k]
                   + f_0 * hi_467[k];

        t_328[k] = -3.0 * fi_216[k]
                   + f_0 * hi_468[k];

        t_329[k] = -3.0 * fi_217[k]
                   + f_0 * hi_469[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, fi_218, fi_219, fi_220, fi_221, \
                         fi_222, hi_470, hi_471, hi_472, hi_473, \
                         hi_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -3.0 * fi_218[k]
                   + f_0 * hi_470[k];

        t_331[k] = -3.0 * fi_219[k]
                   + f_0 * hi_471[k];

        t_332[k] = -3.0 * fi_220[k]
                   + f_0 * hi_472[k];

        t_333[k] = -3.0 * fi_221[k]
                   + f_0 * hi_473[k];

        t_334[k] = -3.0 * fi_222[k]
                   + f_0 * hi_474[k];
    }
}

static auto
compute_prim_geom_10_gi_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fi, const size_t hi,
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

    const auto *fi_223 = buffer.data(fi + 223);
    const auto *fi_224 = buffer.data(fi + 224);
    const auto *fi_225 = buffer.data(fi + 225);
    const auto *fi_226 = buffer.data(fi + 226);
    const auto *fi_227 = buffer.data(fi + 227);
    const auto *fi_228 = buffer.data(fi + 228);
    const auto *fi_229 = buffer.data(fi + 229);
    const auto *fi_230 = buffer.data(fi + 230);
    const auto *fi_231 = buffer.data(fi + 231);
    const auto *fi_232 = buffer.data(fi + 232);
    const auto *fi_233 = buffer.data(fi + 233);
    const auto *fi_234 = buffer.data(fi + 234);
    const auto *fi_235 = buffer.data(fi + 235);
    const auto *fi_236 = buffer.data(fi + 236);
    const auto *fi_237 = buffer.data(fi + 237);
    const auto *fi_238 = buffer.data(fi + 238);
    const auto *fi_239 = buffer.data(fi + 239);
    const auto *fi_240 = buffer.data(fi + 240);
    const auto *fi_241 = buffer.data(fi + 241);
    const auto *fi_242 = buffer.data(fi + 242);
    const auto *fi_243 = buffer.data(fi + 243);
    const auto *fi_244 = buffer.data(fi + 244);
    const auto *fi_245 = buffer.data(fi + 245);
    const auto *fi_246 = buffer.data(fi + 246);
    const auto *fi_247 = buffer.data(fi + 247);
    const auto *fi_248 = buffer.data(fi + 248);
    const auto *fi_249 = buffer.data(fi + 249);
    const auto *fi_250 = buffer.data(fi + 250);
    const auto *fi_251 = buffer.data(fi + 251);
    const auto *fi_252 = buffer.data(fi + 252);
    const auto *fi_253 = buffer.data(fi + 253);
    const auto *fi_254 = buffer.data(fi + 254);
    const auto *fi_255 = buffer.data(fi + 255);
    const auto *fi_256 = buffer.data(fi + 256);
    const auto *fi_257 = buffer.data(fi + 257);
    const auto *fi_258 = buffer.data(fi + 258);
    const auto *fi_259 = buffer.data(fi + 259);
    const auto *fi_260 = buffer.data(fi + 260);
    const auto *fi_261 = buffer.data(fi + 261);
    const auto *fi_262 = buffer.data(fi + 262);
    const auto *fi_263 = buffer.data(fi + 263);
    const auto *fi_264 = buffer.data(fi + 264);
    const auto *fi_265 = buffer.data(fi + 265);
    const auto *fi_266 = buffer.data(fi + 266);
    const auto *fi_267 = buffer.data(fi + 267);
    const auto *fi_268 = buffer.data(fi + 268);
    const auto *fi_269 = buffer.data(fi + 269);
    const auto *fi_270 = buffer.data(fi + 270);
    const auto *fi_271 = buffer.data(fi + 271);
    const auto *fi_272 = buffer.data(fi + 272);
    const auto *fi_273 = buffer.data(fi + 273);
    const auto *fi_274 = buffer.data(fi + 274);
    const auto *fi_275 = buffer.data(fi + 275);
    const auto *fi_276 = buffer.data(fi + 276);
    const auto *fi_277 = buffer.data(fi + 277);
    const auto *fi_278 = buffer.data(fi + 278);
    const auto *fi_279 = buffer.data(fi + 279);

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

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, fi_223, fi_224, fi_225, fi_226, \
                         fi_227, hi_475, hi_476, hi_477, hi_478, \
                         hi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -3.0 * fi_223[k]
                   + f_0 * hi_475[k];

        t_336[k] = -2.0 * fi_224[k]
                   + f_0 * hi_476[k];

        t_337[k] = -2.0 * fi_225[k]
                   + f_0 * hi_477[k];

        t_338[k] = -2.0 * fi_226[k]
                   + f_0 * hi_478[k];

        t_339[k] = -2.0 * fi_227[k]
                   + f_0 * hi_479[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, fi_228, fi_229, fi_230, fi_231, \
                         fi_232, hi_480, hi_481, hi_482, hi_483, \
                         hi_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * fi_228[k]
                   + f_0 * hi_480[k];

        t_341[k] = -2.0 * fi_229[k]
                   + f_0 * hi_481[k];

        t_342[k] = -2.0 * fi_230[k]
                   + f_0 * hi_482[k];

        t_343[k] = -2.0 * fi_231[k]
                   + f_0 * hi_483[k];

        t_344[k] = -2.0 * fi_232[k]
                   + f_0 * hi_484[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, fi_233, fi_234, fi_235, fi_236, \
                         fi_237, hi_485, hi_486, hi_487, hi_488, \
                         hi_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * fi_233[k]
                   + f_0 * hi_485[k];

        t_346[k] = -2.0 * fi_234[k]
                   + f_0 * hi_486[k];

        t_347[k] = -2.0 * fi_235[k]
                   + f_0 * hi_487[k];

        t_348[k] = -2.0 * fi_236[k]
                   + f_0 * hi_488[k];

        t_349[k] = -2.0 * fi_237[k]
                   + f_0 * hi_489[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, fi_238, fi_239, fi_240, fi_241, \
                         fi_242, hi_490, hi_491, hi_492, hi_493, \
                         hi_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * fi_238[k]
                   + f_0 * hi_490[k];

        t_351[k] = -2.0 * fi_239[k]
                   + f_0 * hi_491[k];

        t_352[k] = -2.0 * fi_240[k]
                   + f_0 * hi_492[k];

        t_353[k] = -2.0 * fi_241[k]
                   + f_0 * hi_493[k];

        t_354[k] = -2.0 * fi_242[k]
                   + f_0 * hi_494[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, fi_243, fi_244, fi_245, fi_246, \
                         fi_247, hi_495, hi_496, hi_497, hi_498, \
                         hi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * fi_243[k]
                   + f_0 * hi_495[k];

        t_356[k] = -2.0 * fi_244[k]
                   + f_0 * hi_496[k];

        t_357[k] = -2.0 * fi_245[k]
                   + f_0 * hi_497[k];

        t_358[k] = -2.0 * fi_246[k]
                   + f_0 * hi_498[k];

        t_359[k] = -2.0 * fi_247[k]
                   + f_0 * hi_499[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, fi_248, fi_249, fi_250, fi_251, \
                         fi_252, hi_500, hi_501, hi_502, hi_503, \
                         hi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -2.0 * fi_248[k]
                   + f_0 * hi_500[k];

        t_361[k] = -2.0 * fi_249[k]
                   + f_0 * hi_501[k];

        t_362[k] = -2.0 * fi_250[k]
                   + f_0 * hi_502[k];

        t_363[k] = -2.0 * fi_251[k]
                   + f_0 * hi_503[k];

        t_364[k] = -fi_252[k]
                   + f_0 * hi_504[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, fi_253, fi_254, fi_255, fi_256, \
                         fi_257, hi_505, hi_506, hi_507, hi_508, \
                         hi_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -fi_253[k]
                   + f_0 * hi_505[k];

        t_366[k] = -fi_254[k]
                   + f_0 * hi_506[k];

        t_367[k] = -fi_255[k]
                   + f_0 * hi_507[k];

        t_368[k] = -fi_256[k]
                   + f_0 * hi_508[k];

        t_369[k] = -fi_257[k]
                   + f_0 * hi_509[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, fi_258, fi_259, fi_260, fi_261, \
                         fi_262, hi_510, hi_511, hi_512, hi_513, \
                         hi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -fi_258[k]
                   + f_0 * hi_510[k];

        t_371[k] = -fi_259[k]
                   + f_0 * hi_511[k];

        t_372[k] = -fi_260[k]
                   + f_0 * hi_512[k];

        t_373[k] = -fi_261[k]
                   + f_0 * hi_513[k];

        t_374[k] = -fi_262[k]
                   + f_0 * hi_514[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, fi_263, fi_264, fi_265, fi_266, \
                         fi_267, hi_515, hi_516, hi_517, hi_518, \
                         hi_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -fi_263[k]
                   + f_0 * hi_515[k];

        t_376[k] = -fi_264[k]
                   + f_0 * hi_516[k];

        t_377[k] = -fi_265[k]
                   + f_0 * hi_517[k];

        t_378[k] = -fi_266[k]
                   + f_0 * hi_518[k];

        t_379[k] = -fi_267[k]
                   + f_0 * hi_519[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, fi_268, fi_269, fi_270, fi_271, \
                         fi_272, hi_520, hi_521, hi_522, hi_523, \
                         hi_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -fi_268[k]
                   + f_0 * hi_520[k];

        t_381[k] = -fi_269[k]
                   + f_0 * hi_521[k];

        t_382[k] = -fi_270[k]
                   + f_0 * hi_522[k];

        t_383[k] = -fi_271[k]
                   + f_0 * hi_523[k];

        t_384[k] = -fi_272[k]
                   + f_0 * hi_524[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, fi_273, fi_274, fi_275, fi_276, \
                         fi_277, hi_525, hi_526, hi_527, hi_528, \
                         hi_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -fi_273[k]
                   + f_0 * hi_525[k];

        t_386[k] = -fi_274[k]
                   + f_0 * hi_526[k];

        t_387[k] = -fi_275[k]
                   + f_0 * hi_527[k];

        t_388[k] = -fi_276[k]
                   + f_0 * hi_528[k];

        t_389[k] = -fi_277[k]
                   + f_0 * hi_529[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, t_395, t_396, fi_278, fi_279, \
                         hi_530, hi_531, hi_532, hi_533, hi_534, hi_535, \
                         hi_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -fi_278[k]
                   + f_0 * hi_530[k];

        t_391[k] = -fi_279[k]
                   + f_0 * hi_531[k];

        t_392[k] = f_0 * hi_532[k];

        t_393[k] = f_0 * hi_533[k];

        t_394[k] = f_0 * hi_534[k];

        t_395[k] = f_0 * hi_535[k];

        t_396[k] = f_0 * hi_536[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, t_402, t_403, t_404, hi_537, \
                         hi_538, hi_539, hi_540, hi_541, hi_542, hi_543, \
                         hi_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_0 * hi_537[k];

        t_398[k] = f_0 * hi_538[k];

        t_399[k] = f_0 * hi_539[k];

        t_400[k] = f_0 * hi_540[k];

        t_401[k] = f_0 * hi_541[k];

        t_402[k] = f_0 * hi_542[k];

        t_403[k] = f_0 * hi_543[k];

        t_404[k] = f_0 * hi_544[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, hi_545, \
                         hi_546, hi_547, hi_548, hi_549, hi_550, hi_551, \
                         hi_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * hi_545[k];

        t_406[k] = f_0 * hi_546[k];

        t_407[k] = f_0 * hi_547[k];

        t_408[k] = f_0 * hi_548[k];

        t_409[k] = f_0 * hi_549[k];

        t_410[k] = f_0 * hi_550[k];

        t_411[k] = f_0 * hi_551[k];

        t_412[k] = f_0 * hi_552[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, hi_553, hi_554, \
                         hi_555, hi_556, hi_557, hi_558, hi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * hi_553[k];

        t_414[k] = f_0 * hi_554[k];

        t_415[k] = f_0 * hi_555[k];

        t_416[k] = f_0 * hi_556[k];

        t_417[k] = f_0 * hi_557[k];

        t_418[k] = f_0 * hi_558[k];

        t_419[k] = f_0 * hi_559[k];
    }
}

auto
compute_prim_geom_10_gi_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t fi, const size_t hi,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gi_electron_repulsion_1_piece0(buffer, target, fi, hi, ncols, alpha);

    compute_prim_geom_10_gi_electron_repulsion_1_piece1(buffer, target, fi, hi, ncols, alpha);

    compute_prim_geom_10_gi_electron_repulsion_1_piece2(buffer, target, fi, hi, ncols, alpha);
}

static auto
compute_prim_geom_10_gi_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fi, const size_t hi,
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

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_41 = buffer.data(fi + 41);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_60 = buffer.data(fi + 60);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_63 = buffer.data(fi + 63);
    const auto *fi_64 = buffer.data(fi + 64);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_66 = buffer.data(fi + 66);
    const auto *fi_67 = buffer.data(fi + 67);
    const auto *fi_68 = buffer.data(fi + 68);
    const auto *fi_69 = buffer.data(fi + 69);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_71 = buffer.data(fi + 71);
    const auto *fi_72 = buffer.data(fi + 72);
    const auto *fi_73 = buffer.data(fi + 73);
    const auto *fi_74 = buffer.data(fi + 74);
    const auto *fi_75 = buffer.data(fi + 75);
    const auto *fi_76 = buffer.data(fi + 76);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_78 = buffer.data(fi + 78);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);
    const auto *fi_81 = buffer.data(fi + 81);
    const auto *fi_82 = buffer.data(fi + 82);
    const auto *fi_83 = buffer.data(fi + 83);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hi_56, hi_57, hi_58, hi_59, \
                         hi_60, hi_61, hi_62, hi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hi_56[k];

        t_1[k] = f_0 * hi_57[k];

        t_2[k] = f_0 * hi_58[k];

        t_3[k] = f_0 * hi_59[k];

        t_4[k] = f_0 * hi_60[k];

        t_5[k] = f_0 * hi_61[k];

        t_6[k] = f_0 * hi_62[k];

        t_7[k] = f_0 * hi_63[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, hi_64, hi_65, hi_66, \
                         hi_67, hi_68, hi_69, hi_70, hi_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hi_64[k];

        t_9[k] = f_0 * hi_65[k];

        t_10[k] = f_0 * hi_66[k];

        t_11[k] = f_0 * hi_67[k];

        t_12[k] = f_0 * hi_68[k];

        t_13[k] = f_0 * hi_69[k];

        t_14[k] = f_0 * hi_70[k];

        t_15[k] = f_0 * hi_71[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, hi_72, hi_73, hi_74, \
                         hi_75, hi_76, hi_77, hi_78, hi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * hi_72[k];

        t_17[k] = f_0 * hi_73[k];

        t_18[k] = f_0 * hi_74[k];

        t_19[k] = f_0 * hi_75[k];

        t_20[k] = f_0 * hi_76[k];

        t_21[k] = f_0 * hi_77[k];

        t_22[k] = f_0 * hi_78[k];

        t_23[k] = f_0 * hi_79[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, hi_80, hi_81, hi_82, \
                         hi_83, hi_112, hi_113, hi_114, hi_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * hi_80[k];

        t_25[k] = f_0 * hi_81[k];

        t_26[k] = f_0 * hi_82[k];

        t_27[k] = f_0 * hi_83[k];

        t_28[k] = f_0 * hi_112[k];

        t_29[k] = f_0 * hi_113[k];

        t_30[k] = f_0 * hi_114[k];

        t_31[k] = f_0 * hi_115[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, hi_116, hi_117, \
                         hi_118, hi_119, hi_120, hi_121, hi_122, \
                         hi_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * hi_116[k];

        t_33[k] = f_0 * hi_117[k];

        t_34[k] = f_0 * hi_118[k];

        t_35[k] = f_0 * hi_119[k];

        t_36[k] = f_0 * hi_120[k];

        t_37[k] = f_0 * hi_121[k];

        t_38[k] = f_0 * hi_122[k];

        t_39[k] = f_0 * hi_123[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, hi_124, hi_125, \
                         hi_126, hi_127, hi_128, hi_129, hi_130, \
                         hi_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * hi_124[k];

        t_41[k] = f_0 * hi_125[k];

        t_42[k] = f_0 * hi_126[k];

        t_43[k] = f_0 * hi_127[k];

        t_44[k] = f_0 * hi_128[k];

        t_45[k] = f_0 * hi_129[k];

        t_46[k] = f_0 * hi_130[k];

        t_47[k] = f_0 * hi_131[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, hi_132, hi_133, \
                         hi_134, hi_135, hi_136, hi_137, hi_138, \
                         hi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * hi_132[k];

        t_49[k] = f_0 * hi_133[k];

        t_50[k] = f_0 * hi_134[k];

        t_51[k] = f_0 * hi_135[k];

        t_52[k] = f_0 * hi_136[k];

        t_53[k] = f_0 * hi_137[k];

        t_54[k] = f_0 * hi_138[k];

        t_55[k] = f_0 * hi_139[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, fi_0, fi_1, fi_2, fi_3, fi_4, hi_140, \
                         hi_141, hi_142, hi_143, hi_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -fi_0[k]
                  + f_0 * hi_140[k];

        t_57[k] = -fi_1[k]
                  + f_0 * hi_141[k];

        t_58[k] = -fi_2[k]
                  + f_0 * hi_142[k];

        t_59[k] = -fi_3[k]
                  + f_0 * hi_143[k];

        t_60[k] = -fi_4[k]
                  + f_0 * hi_144[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, fi_5, fi_6, fi_7, fi_8, fi_9, hi_145, \
                         hi_146, hi_147, hi_148, hi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -fi_5[k]
                  + f_0 * hi_145[k];

        t_62[k] = -fi_6[k]
                  + f_0 * hi_146[k];

        t_63[k] = -fi_7[k]
                  + f_0 * hi_147[k];

        t_64[k] = -fi_8[k]
                  + f_0 * hi_148[k];

        t_65[k] = -fi_9[k]
                  + f_0 * hi_149[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, fi_10, fi_11, fi_12, fi_13, fi_14, \
                         hi_150, hi_151, hi_152, hi_153, hi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -fi_10[k]
                  + f_0 * hi_150[k];

        t_67[k] = -fi_11[k]
                  + f_0 * hi_151[k];

        t_68[k] = -fi_12[k]
                  + f_0 * hi_152[k];

        t_69[k] = -fi_13[k]
                  + f_0 * hi_153[k];

        t_70[k] = -fi_14[k]
                  + f_0 * hi_154[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, fi_15, fi_16, fi_17, fi_18, fi_19, \
                         hi_155, hi_156, hi_157, hi_158, hi_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -fi_15[k]
                  + f_0 * hi_155[k];

        t_72[k] = -fi_16[k]
                  + f_0 * hi_156[k];

        t_73[k] = -fi_17[k]
                  + f_0 * hi_157[k];

        t_74[k] = -fi_18[k]
                  + f_0 * hi_158[k];

        t_75[k] = -fi_19[k]
                  + f_0 * hi_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, fi_20, fi_21, fi_22, fi_23, fi_24, \
                         hi_160, hi_161, hi_162, hi_163, hi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -fi_20[k]
                  + f_0 * hi_160[k];

        t_77[k] = -fi_21[k]
                  + f_0 * hi_161[k];

        t_78[k] = -fi_22[k]
                  + f_0 * hi_162[k];

        t_79[k] = -fi_23[k]
                  + f_0 * hi_163[k];

        t_80[k] = -fi_24[k]
                  + f_0 * hi_164[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, fi_25, fi_26, fi_27, hi_165, \
                         hi_166, hi_167, hi_196, hi_197, hi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -fi_25[k]
                  + f_0 * hi_165[k];

        t_82[k] = -fi_26[k]
                  + f_0 * hi_166[k];

        t_83[k] = -fi_27[k]
                  + f_0 * hi_167[k];

        t_84[k] = f_0 * hi_196[k];

        t_85[k] = f_0 * hi_197[k];

        t_86[k] = f_0 * hi_198[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, t_93, t_94, hi_199, hi_200, \
                         hi_201, hi_202, hi_203, hi_204, hi_205, \
                         hi_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_0 * hi_199[k];

        t_88[k] = f_0 * hi_200[k];

        t_89[k] = f_0 * hi_201[k];

        t_90[k] = f_0 * hi_202[k];

        t_91[k] = f_0 * hi_203[k];

        t_92[k] = f_0 * hi_204[k];

        t_93[k] = f_0 * hi_205[k];

        t_94[k] = f_0 * hi_206[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, t_100, t_101, t_102, hi_207, hi_208, \
                         hi_209, hi_210, hi_211, hi_212, hi_213, \
                         hi_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_0 * hi_207[k];

        t_96[k] = f_0 * hi_208[k];

        t_97[k] = f_0 * hi_209[k];

        t_98[k] = f_0 * hi_210[k];

        t_99[k] = f_0 * hi_211[k];

        t_100[k] = f_0 * hi_212[k];

        t_101[k] = f_0 * hi_213[k];

        t_102[k] = f_0 * hi_214[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, t_109, t_110, hi_215, \
                         hi_216, hi_217, hi_218, hi_219, hi_220, hi_221, \
                         hi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_0 * hi_215[k];

        t_104[k] = f_0 * hi_216[k];

        t_105[k] = f_0 * hi_217[k];

        t_106[k] = f_0 * hi_218[k];

        t_107[k] = f_0 * hi_219[k];

        t_108[k] = f_0 * hi_220[k];

        t_109[k] = f_0 * hi_221[k];

        t_110[k] = f_0 * hi_222[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, fi_28, fi_29, fi_30, fi_31, \
                         hi_223, hi_224, hi_225, hi_226, hi_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_0 * hi_223[k];

        t_112[k] = -fi_28[k]
                   + f_0 * hi_224[k];

        t_113[k] = -fi_29[k]
                   + f_0 * hi_225[k];

        t_114[k] = -fi_30[k]
                   + f_0 * hi_226[k];

        t_115[k] = -fi_31[k]
                   + f_0 * hi_227[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, fi_32, fi_33, fi_34, fi_35, fi_36, \
                         hi_228, hi_229, hi_230, hi_231, hi_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = -fi_32[k]
                   + f_0 * hi_228[k];

        t_117[k] = -fi_33[k]
                   + f_0 * hi_229[k];

        t_118[k] = -fi_34[k]
                   + f_0 * hi_230[k];

        t_119[k] = -fi_35[k]
                   + f_0 * hi_231[k];

        t_120[k] = -fi_36[k]
                   + f_0 * hi_232[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, fi_37, fi_38, fi_39, fi_40, fi_41, \
                         hi_233, hi_234, hi_235, hi_236, hi_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -fi_37[k]
                   + f_0 * hi_233[k];

        t_122[k] = -fi_38[k]
                   + f_0 * hi_234[k];

        t_123[k] = -fi_39[k]
                   + f_0 * hi_235[k];

        t_124[k] = -fi_40[k]
                   + f_0 * hi_236[k];

        t_125[k] = -fi_41[k]
                   + f_0 * hi_237[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, fi_42, fi_43, fi_44, fi_45, fi_46, \
                         hi_238, hi_239, hi_240, hi_241, hi_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = -fi_42[k]
                   + f_0 * hi_238[k];

        t_127[k] = -fi_43[k]
                   + f_0 * hi_239[k];

        t_128[k] = -fi_44[k]
                   + f_0 * hi_240[k];

        t_129[k] = -fi_45[k]
                   + f_0 * hi_241[k];

        t_130[k] = -fi_46[k]
                   + f_0 * hi_242[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, fi_47, fi_48, fi_49, fi_50, fi_51, \
                         hi_243, hi_244, hi_245, hi_246, hi_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = -fi_47[k]
                   + f_0 * hi_243[k];

        t_132[k] = -fi_48[k]
                   + f_0 * hi_244[k];

        t_133[k] = -fi_49[k]
                   + f_0 * hi_245[k];

        t_134[k] = -fi_50[k]
                   + f_0 * hi_246[k];

        t_135[k] = -fi_51[k]
                   + f_0 * hi_247[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, fi_52, fi_53, fi_54, fi_55, fi_56, \
                         hi_248, hi_249, hi_250, hi_251, hi_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = -fi_52[k]
                   + f_0 * hi_248[k];

        t_137[k] = -fi_53[k]
                   + f_0 * hi_249[k];

        t_138[k] = -fi_54[k]
                   + f_0 * hi_250[k];

        t_139[k] = -fi_55[k]
                   + f_0 * hi_251[k];

        t_140[k] = -2.0 * fi_56[k]
                   + f_0 * hi_252[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, fi_57, fi_58, fi_59, fi_60, fi_61, \
                         hi_253, hi_254, hi_255, hi_256, hi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = -2.0 * fi_57[k]
                   + f_0 * hi_253[k];

        t_142[k] = -2.0 * fi_58[k]
                   + f_0 * hi_254[k];

        t_143[k] = -2.0 * fi_59[k]
                   + f_0 * hi_255[k];

        t_144[k] = -2.0 * fi_60[k]
                   + f_0 * hi_256[k];

        t_145[k] = -2.0 * fi_61[k]
                   + f_0 * hi_257[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, fi_62, fi_63, fi_64, fi_65, fi_66, \
                         hi_258, hi_259, hi_260, hi_261, hi_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = -2.0 * fi_62[k]
                   + f_0 * hi_258[k];

        t_147[k] = -2.0 * fi_63[k]
                   + f_0 * hi_259[k];

        t_148[k] = -2.0 * fi_64[k]
                   + f_0 * hi_260[k];

        t_149[k] = -2.0 * fi_65[k]
                   + f_0 * hi_261[k];

        t_150[k] = -2.0 * fi_66[k]
                   + f_0 * hi_262[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, fi_67, fi_68, fi_69, fi_70, fi_71, \
                         hi_263, hi_264, hi_265, hi_266, hi_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = -2.0 * fi_67[k]
                   + f_0 * hi_263[k];

        t_152[k] = -2.0 * fi_68[k]
                   + f_0 * hi_264[k];

        t_153[k] = -2.0 * fi_69[k]
                   + f_0 * hi_265[k];

        t_154[k] = -2.0 * fi_70[k]
                   + f_0 * hi_266[k];

        t_155[k] = -2.0 * fi_71[k]
                   + f_0 * hi_267[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, fi_72, fi_73, fi_74, fi_75, fi_76, \
                         hi_268, hi_269, hi_270, hi_271, hi_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = -2.0 * fi_72[k]
                   + f_0 * hi_268[k];

        t_157[k] = -2.0 * fi_73[k]
                   + f_0 * hi_269[k];

        t_158[k] = -2.0 * fi_74[k]
                   + f_0 * hi_270[k];

        t_159[k] = -2.0 * fi_75[k]
                   + f_0 * hi_271[k];

        t_160[k] = -2.0 * fi_76[k]
                   + f_0 * hi_272[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, fi_77, fi_78, fi_79, fi_80, fi_81, \
                         hi_273, hi_274, hi_275, hi_276, hi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = -2.0 * fi_77[k]
                   + f_0 * hi_273[k];

        t_162[k] = -2.0 * fi_78[k]
                   + f_0 * hi_274[k];

        t_163[k] = -2.0 * fi_79[k]
                   + f_0 * hi_275[k];

        t_164[k] = -2.0 * fi_80[k]
                   + f_0 * hi_276[k];

        t_165[k] = -2.0 * fi_81[k]
                   + f_0 * hi_277[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, t_171, t_172, fi_82, fi_83, \
                         hi_278, hi_279, hi_308, hi_309, hi_310, hi_311, \
                         hi_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = -2.0 * fi_82[k]
                   + f_0 * hi_278[k];

        t_167[k] = -2.0 * fi_83[k]
                   + f_0 * hi_279[k];

        t_168[k] = f_0 * hi_308[k];

        t_169[k] = f_0 * hi_309[k];

        t_170[k] = f_0 * hi_310[k];

        t_171[k] = f_0 * hi_311[k];

        t_172[k] = f_0 * hi_312[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, t_179, t_180, hi_313, \
                         hi_314, hi_315, hi_316, hi_317, hi_318, hi_319, \
                         hi_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_0 * hi_313[k];

        t_174[k] = f_0 * hi_314[k];

        t_175[k] = f_0 * hi_315[k];

        t_176[k] = f_0 * hi_316[k];

        t_177[k] = f_0 * hi_317[k];

        t_178[k] = f_0 * hi_318[k];

        t_179[k] = f_0 * hi_319[k];

        t_180[k] = f_0 * hi_320[k];
    }
}

static auto
compute_prim_geom_10_gi_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fi, const size_t hi,
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

    const auto *fi_84 = buffer.data(fi + 84);
    const auto *fi_85 = buffer.data(fi + 85);
    const auto *fi_86 = buffer.data(fi + 86);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_88 = buffer.data(fi + 88);
    const auto *fi_89 = buffer.data(fi + 89);
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_91 = buffer.data(fi + 91);
    const auto *fi_92 = buffer.data(fi + 92);
    const auto *fi_93 = buffer.data(fi + 93);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_95 = buffer.data(fi + 95);
    const auto *fi_96 = buffer.data(fi + 96);
    const auto *fi_97 = buffer.data(fi + 97);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_100 = buffer.data(fi + 100);
    const auto *fi_101 = buffer.data(fi + 101);
    const auto *fi_102 = buffer.data(fi + 102);
    const auto *fi_103 = buffer.data(fi + 103);
    const auto *fi_104 = buffer.data(fi + 104);
    const auto *fi_105 = buffer.data(fi + 105);
    const auto *fi_106 = buffer.data(fi + 106);
    const auto *fi_107 = buffer.data(fi + 107);
    const auto *fi_108 = buffer.data(fi + 108);
    const auto *fi_109 = buffer.data(fi + 109);
    const auto *fi_110 = buffer.data(fi + 110);
    const auto *fi_111 = buffer.data(fi + 111);
    const auto *fi_112 = buffer.data(fi + 112);
    const auto *fi_113 = buffer.data(fi + 113);
    const auto *fi_114 = buffer.data(fi + 114);
    const auto *fi_115 = buffer.data(fi + 115);
    const auto *fi_116 = buffer.data(fi + 116);
    const auto *fi_117 = buffer.data(fi + 117);
    const auto *fi_118 = buffer.data(fi + 118);
    const auto *fi_119 = buffer.data(fi + 119);
    const auto *fi_120 = buffer.data(fi + 120);
    const auto *fi_121 = buffer.data(fi + 121);
    const auto *fi_122 = buffer.data(fi + 122);
    const auto *fi_123 = buffer.data(fi + 123);
    const auto *fi_124 = buffer.data(fi + 124);
    const auto *fi_125 = buffer.data(fi + 125);
    const auto *fi_126 = buffer.data(fi + 126);
    const auto *fi_127 = buffer.data(fi + 127);
    const auto *fi_128 = buffer.data(fi + 128);
    const auto *fi_129 = buffer.data(fi + 129);
    const auto *fi_130 = buffer.data(fi + 130);
    const auto *fi_131 = buffer.data(fi + 131);
    const auto *fi_132 = buffer.data(fi + 132);
    const auto *fi_133 = buffer.data(fi + 133);
    const auto *fi_134 = buffer.data(fi + 134);
    const auto *fi_135 = buffer.data(fi + 135);
    const auto *fi_136 = buffer.data(fi + 136);
    const auto *fi_137 = buffer.data(fi + 137);
    const auto *fi_138 = buffer.data(fi + 138);
    const auto *fi_139 = buffer.data(fi + 139);
    const auto *fi_140 = buffer.data(fi + 140);
    const auto *fi_141 = buffer.data(fi + 141);
    const auto *fi_142 = buffer.data(fi + 142);
    const auto *fi_143 = buffer.data(fi + 143);
    const auto *fi_144 = buffer.data(fi + 144);
    const auto *fi_145 = buffer.data(fi + 145);
    const auto *fi_146 = buffer.data(fi + 146);
    const auto *fi_147 = buffer.data(fi + 147);
    const auto *fi_148 = buffer.data(fi + 148);
    const auto *fi_149 = buffer.data(fi + 149);
    const auto *fi_150 = buffer.data(fi + 150);
    const auto *fi_151 = buffer.data(fi + 151);
    const auto *fi_152 = buffer.data(fi + 152);
    const auto *fi_153 = buffer.data(fi + 153);
    const auto *fi_154 = buffer.data(fi + 154);
    const auto *fi_155 = buffer.data(fi + 155);
    const auto *fi_156 = buffer.data(fi + 156);
    const auto *fi_157 = buffer.data(fi + 157);
    const auto *fi_158 = buffer.data(fi + 158);
    const auto *fi_159 = buffer.data(fi + 159);
    const auto *fi_160 = buffer.data(fi + 160);
    const auto *fi_161 = buffer.data(fi + 161);
    const auto *fi_162 = buffer.data(fi + 162);
    const auto *fi_163 = buffer.data(fi + 163);
    const auto *fi_164 = buffer.data(fi + 164);
    const auto *fi_165 = buffer.data(fi + 165);
    const auto *fi_166 = buffer.data(fi + 166);
    const auto *fi_167 = buffer.data(fi + 167);
    const auto *fi_168 = buffer.data(fi + 168);
    const auto *fi_169 = buffer.data(fi + 169);
    const auto *fi_170 = buffer.data(fi + 170);
    const auto *fi_171 = buffer.data(fi + 171);
    const auto *fi_172 = buffer.data(fi + 172);
    const auto *fi_173 = buffer.data(fi + 173);
    const auto *fi_174 = buffer.data(fi + 174);
    const auto *fi_175 = buffer.data(fi + 175);
    const auto *fi_176 = buffer.data(fi + 176);
    const auto *fi_177 = buffer.data(fi + 177);
    const auto *fi_178 = buffer.data(fi + 178);
    const auto *fi_179 = buffer.data(fi + 179);
    const auto *fi_180 = buffer.data(fi + 180);
    const auto *fi_181 = buffer.data(fi + 181);
    const auto *fi_182 = buffer.data(fi + 182);
    const auto *fi_183 = buffer.data(fi + 183);
    const auto *fi_184 = buffer.data(fi + 184);
    const auto *fi_185 = buffer.data(fi + 185);
    const auto *fi_186 = buffer.data(fi + 186);
    const auto *fi_187 = buffer.data(fi + 187);
    const auto *fi_188 = buffer.data(fi + 188);
    const auto *fi_189 = buffer.data(fi + 189);
    const auto *fi_190 = buffer.data(fi + 190);
    const auto *fi_191 = buffer.data(fi + 191);
    const auto *fi_192 = buffer.data(fi + 192);
    const auto *fi_193 = buffer.data(fi + 193);
    const auto *fi_194 = buffer.data(fi + 194);
    const auto *fi_195 = buffer.data(fi + 195);
    const auto *fi_196 = buffer.data(fi + 196);
    const auto *fi_197 = buffer.data(fi + 197);
    const auto *fi_198 = buffer.data(fi + 198);
    const auto *fi_199 = buffer.data(fi + 199);
    const auto *fi_200 = buffer.data(fi + 200);
    const auto *fi_201 = buffer.data(fi + 201);
    const auto *fi_202 = buffer.data(fi + 202);
    const auto *fi_203 = buffer.data(fi + 203);
    const auto *fi_204 = buffer.data(fi + 204);
    const auto *fi_205 = buffer.data(fi + 205);
    const auto *fi_206 = buffer.data(fi + 206);

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

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, t_186, t_187, t_188, hi_321, \
                         hi_322, hi_323, hi_324, hi_325, hi_326, hi_327, \
                         hi_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_0 * hi_321[k];

        t_182[k] = f_0 * hi_322[k];

        t_183[k] = f_0 * hi_323[k];

        t_184[k] = f_0 * hi_324[k];

        t_185[k] = f_0 * hi_325[k];

        t_186[k] = f_0 * hi_326[k];

        t_187[k] = f_0 * hi_327[k];

        t_188[k] = f_0 * hi_328[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, hi_329, hi_330, \
                         hi_331, hi_332, hi_333, hi_334, hi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * hi_329[k];

        t_190[k] = f_0 * hi_330[k];

        t_191[k] = f_0 * hi_331[k];

        t_192[k] = f_0 * hi_332[k];

        t_193[k] = f_0 * hi_333[k];

        t_194[k] = f_0 * hi_334[k];

        t_195[k] = f_0 * hi_335[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, fi_84, fi_85, fi_86, fi_87, fi_88, \
                         hi_336, hi_337, hi_338, hi_339, hi_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -fi_84[k]
                   + f_0 * hi_336[k];

        t_197[k] = -fi_85[k]
                   + f_0 * hi_337[k];

        t_198[k] = -fi_86[k]
                   + f_0 * hi_338[k];

        t_199[k] = -fi_87[k]
                   + f_0 * hi_339[k];

        t_200[k] = -fi_88[k]
                   + f_0 * hi_340[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, fi_89, fi_90, fi_91, fi_92, fi_93, \
                         hi_341, hi_342, hi_343, hi_344, hi_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = -fi_89[k]
                   + f_0 * hi_341[k];

        t_202[k] = -fi_90[k]
                   + f_0 * hi_342[k];

        t_203[k] = -fi_91[k]
                   + f_0 * hi_343[k];

        t_204[k] = -fi_92[k]
                   + f_0 * hi_344[k];

        t_205[k] = -fi_93[k]
                   + f_0 * hi_345[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, fi_94, fi_95, fi_96, fi_97, fi_98, \
                         hi_346, hi_347, hi_348, hi_349, hi_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = -fi_94[k]
                   + f_0 * hi_346[k];

        t_207[k] = -fi_95[k]
                   + f_0 * hi_347[k];

        t_208[k] = -fi_96[k]
                   + f_0 * hi_348[k];

        t_209[k] = -fi_97[k]
                   + f_0 * hi_349[k];

        t_210[k] = -fi_98[k]
                   + f_0 * hi_350[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, fi_99, fi_100, fi_101, fi_102, \
                         fi_103, hi_351, hi_352, hi_353, hi_354, \
                         hi_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = -fi_99[k]
                   + f_0 * hi_351[k];

        t_212[k] = -fi_100[k]
                   + f_0 * hi_352[k];

        t_213[k] = -fi_101[k]
                   + f_0 * hi_353[k];

        t_214[k] = -fi_102[k]
                   + f_0 * hi_354[k];

        t_215[k] = -fi_103[k]
                   + f_0 * hi_355[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, fi_104, fi_105, fi_106, fi_107, \
                         fi_108, hi_356, hi_357, hi_358, hi_359, \
                         hi_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = -fi_104[k]
                   + f_0 * hi_356[k];

        t_217[k] = -fi_105[k]
                   + f_0 * hi_357[k];

        t_218[k] = -fi_106[k]
                   + f_0 * hi_358[k];

        t_219[k] = -fi_107[k]
                   + f_0 * hi_359[k];

        t_220[k] = -fi_108[k]
                   + f_0 * hi_360[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, fi_109, fi_110, fi_111, fi_112, \
                         fi_113, hi_361, hi_362, hi_363, hi_364, \
                         hi_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = -fi_109[k]
                   + f_0 * hi_361[k];

        t_222[k] = -fi_110[k]
                   + f_0 * hi_362[k];

        t_223[k] = -fi_111[k]
                   + f_0 * hi_363[k];

        t_224[k] = -2.0 * fi_112[k]
                   + f_0 * hi_364[k];

        t_225[k] = -2.0 * fi_113[k]
                   + f_0 * hi_365[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, t_230, fi_114, fi_115, fi_116, fi_117, \
                         fi_118, hi_366, hi_367, hi_368, hi_369, \
                         hi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = -2.0 * fi_114[k]
                   + f_0 * hi_366[k];

        t_227[k] = -2.0 * fi_115[k]
                   + f_0 * hi_367[k];

        t_228[k] = -2.0 * fi_116[k]
                   + f_0 * hi_368[k];

        t_229[k] = -2.0 * fi_117[k]
                   + f_0 * hi_369[k];

        t_230[k] = -2.0 * fi_118[k]
                   + f_0 * hi_370[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, fi_119, fi_120, fi_121, fi_122, \
                         fi_123, hi_371, hi_372, hi_373, hi_374, \
                         hi_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = -2.0 * fi_119[k]
                   + f_0 * hi_371[k];

        t_232[k] = -2.0 * fi_120[k]
                   + f_0 * hi_372[k];

        t_233[k] = -2.0 * fi_121[k]
                   + f_0 * hi_373[k];

        t_234[k] = -2.0 * fi_122[k]
                   + f_0 * hi_374[k];

        t_235[k] = -2.0 * fi_123[k]
                   + f_0 * hi_375[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, fi_124, fi_125, fi_126, fi_127, \
                         fi_128, hi_376, hi_377, hi_378, hi_379, \
                         hi_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = -2.0 * fi_124[k]
                   + f_0 * hi_376[k];

        t_237[k] = -2.0 * fi_125[k]
                   + f_0 * hi_377[k];

        t_238[k] = -2.0 * fi_126[k]
                   + f_0 * hi_378[k];

        t_239[k] = -2.0 * fi_127[k]
                   + f_0 * hi_379[k];

        t_240[k] = -2.0 * fi_128[k]
                   + f_0 * hi_380[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, fi_129, fi_130, fi_131, fi_132, \
                         fi_133, hi_381, hi_382, hi_383, hi_384, \
                         hi_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = -2.0 * fi_129[k]
                   + f_0 * hi_381[k];

        t_242[k] = -2.0 * fi_130[k]
                   + f_0 * hi_382[k];

        t_243[k] = -2.0 * fi_131[k]
                   + f_0 * hi_383[k];

        t_244[k] = -2.0 * fi_132[k]
                   + f_0 * hi_384[k];

        t_245[k] = -2.0 * fi_133[k]
                   + f_0 * hi_385[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, fi_134, fi_135, fi_136, fi_137, \
                         fi_138, hi_386, hi_387, hi_388, hi_389, \
                         hi_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = -2.0 * fi_134[k]
                   + f_0 * hi_386[k];

        t_247[k] = -2.0 * fi_135[k]
                   + f_0 * hi_387[k];

        t_248[k] = -2.0 * fi_136[k]
                   + f_0 * hi_388[k];

        t_249[k] = -2.0 * fi_137[k]
                   + f_0 * hi_389[k];

        t_250[k] = -2.0 * fi_138[k]
                   + f_0 * hi_390[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, fi_139, fi_140, fi_141, fi_142, \
                         fi_143, hi_391, hi_392, hi_393, hi_394, \
                         hi_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = -2.0 * fi_139[k]
                   + f_0 * hi_391[k];

        t_252[k] = -3.0 * fi_140[k]
                   + f_0 * hi_392[k];

        t_253[k] = -3.0 * fi_141[k]
                   + f_0 * hi_393[k];

        t_254[k] = -3.0 * fi_142[k]
                   + f_0 * hi_394[k];

        t_255[k] = -3.0 * fi_143[k]
                   + f_0 * hi_395[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, fi_144, fi_145, fi_146, fi_147, \
                         fi_148, hi_396, hi_397, hi_398, hi_399, \
                         hi_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = -3.0 * fi_144[k]
                   + f_0 * hi_396[k];

        t_257[k] = -3.0 * fi_145[k]
                   + f_0 * hi_397[k];

        t_258[k] = -3.0 * fi_146[k]
                   + f_0 * hi_398[k];

        t_259[k] = -3.0 * fi_147[k]
                   + f_0 * hi_399[k];

        t_260[k] = -3.0 * fi_148[k]
                   + f_0 * hi_400[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, fi_149, fi_150, fi_151, fi_152, \
                         fi_153, hi_401, hi_402, hi_403, hi_404, \
                         hi_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = -3.0 * fi_149[k]
                   + f_0 * hi_401[k];

        t_262[k] = -3.0 * fi_150[k]
                   + f_0 * hi_402[k];

        t_263[k] = -3.0 * fi_151[k]
                   + f_0 * hi_403[k];

        t_264[k] = -3.0 * fi_152[k]
                   + f_0 * hi_404[k];

        t_265[k] = -3.0 * fi_153[k]
                   + f_0 * hi_405[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, fi_154, fi_155, fi_156, fi_157, \
                         fi_158, hi_406, hi_407, hi_408, hi_409, \
                         hi_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = -3.0 * fi_154[k]
                   + f_0 * hi_406[k];

        t_267[k] = -3.0 * fi_155[k]
                   + f_0 * hi_407[k];

        t_268[k] = -3.0 * fi_156[k]
                   + f_0 * hi_408[k];

        t_269[k] = -3.0 * fi_157[k]
                   + f_0 * hi_409[k];

        t_270[k] = -3.0 * fi_158[k]
                   + f_0 * hi_410[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, fi_159, fi_160, fi_161, fi_162, \
                         fi_163, hi_411, hi_412, hi_413, hi_414, \
                         hi_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = -3.0 * fi_159[k]
                   + f_0 * hi_411[k];

        t_272[k] = -3.0 * fi_160[k]
                   + f_0 * hi_412[k];

        t_273[k] = -3.0 * fi_161[k]
                   + f_0 * hi_413[k];

        t_274[k] = -3.0 * fi_162[k]
                   + f_0 * hi_414[k];

        t_275[k] = -3.0 * fi_163[k]
                   + f_0 * hi_415[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, t_281, fi_164, fi_165, fi_166, \
                         fi_167, hi_416, hi_417, hi_418, hi_419, hi_448, \
                         hi_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = -3.0 * fi_164[k]
                   + f_0 * hi_416[k];

        t_277[k] = -3.0 * fi_165[k]
                   + f_0 * hi_417[k];

        t_278[k] = -3.0 * fi_166[k]
                   + f_0 * hi_418[k];

        t_279[k] = -3.0 * fi_167[k]
                   + f_0 * hi_419[k];

        t_280[k] = f_0 * hi_448[k];

        t_281[k] = f_0 * hi_449[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, t_287, t_288, t_289, hi_450, \
                         hi_451, hi_452, hi_453, hi_454, hi_455, hi_456, \
                         hi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_0 * hi_450[k];

        t_283[k] = f_0 * hi_451[k];

        t_284[k] = f_0 * hi_452[k];

        t_285[k] = f_0 * hi_453[k];

        t_286[k] = f_0 * hi_454[k];

        t_287[k] = f_0 * hi_455[k];

        t_288[k] = f_0 * hi_456[k];

        t_289[k] = f_0 * hi_457[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, t_295, t_296, t_297, hi_458, \
                         hi_459, hi_460, hi_461, hi_462, hi_463, hi_464, \
                         hi_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_0 * hi_458[k];

        t_291[k] = f_0 * hi_459[k];

        t_292[k] = f_0 * hi_460[k];

        t_293[k] = f_0 * hi_461[k];

        t_294[k] = f_0 * hi_462[k];

        t_295[k] = f_0 * hi_463[k];

        t_296[k] = f_0 * hi_464[k];

        t_297[k] = f_0 * hi_465[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, t_303, t_304, t_305, hi_466, \
                         hi_467, hi_468, hi_469, hi_470, hi_471, hi_472, \
                         hi_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_0 * hi_466[k];

        t_299[k] = f_0 * hi_467[k];

        t_300[k] = f_0 * hi_468[k];

        t_301[k] = f_0 * hi_469[k];

        t_302[k] = f_0 * hi_470[k];

        t_303[k] = f_0 * hi_471[k];

        t_304[k] = f_0 * hi_472[k];

        t_305[k] = f_0 * hi_473[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, t_311, fi_168, fi_169, fi_170, \
                         fi_171, hi_474, hi_475, hi_476, hi_477, hi_478, \
                         hi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_0 * hi_474[k];

        t_307[k] = f_0 * hi_475[k];

        t_308[k] = -fi_168[k]
                   + f_0 * hi_476[k];

        t_309[k] = -fi_169[k]
                   + f_0 * hi_477[k];

        t_310[k] = -fi_170[k]
                   + f_0 * hi_478[k];

        t_311[k] = -fi_171[k]
                   + f_0 * hi_479[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, fi_172, fi_173, fi_174, fi_175, \
                         fi_176, hi_480, hi_481, hi_482, hi_483, \
                         hi_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -fi_172[k]
                   + f_0 * hi_480[k];

        t_313[k] = -fi_173[k]
                   + f_0 * hi_481[k];

        t_314[k] = -fi_174[k]
                   + f_0 * hi_482[k];

        t_315[k] = -fi_175[k]
                   + f_0 * hi_483[k];

        t_316[k] = -fi_176[k]
                   + f_0 * hi_484[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, fi_177, fi_178, fi_179, fi_180, \
                         fi_181, hi_485, hi_486, hi_487, hi_488, \
                         hi_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = -fi_177[k]
                   + f_0 * hi_485[k];

        t_318[k] = -fi_178[k]
                   + f_0 * hi_486[k];

        t_319[k] = -fi_179[k]
                   + f_0 * hi_487[k];

        t_320[k] = -fi_180[k]
                   + f_0 * hi_488[k];

        t_321[k] = -fi_181[k]
                   + f_0 * hi_489[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, fi_182, fi_183, fi_184, fi_185, \
                         fi_186, hi_490, hi_491, hi_492, hi_493, \
                         hi_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = -fi_182[k]
                   + f_0 * hi_490[k];

        t_323[k] = -fi_183[k]
                   + f_0 * hi_491[k];

        t_324[k] = -fi_184[k]
                   + f_0 * hi_492[k];

        t_325[k] = -fi_185[k]
                   + f_0 * hi_493[k];

        t_326[k] = -fi_186[k]
                   + f_0 * hi_494[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, fi_187, fi_188, fi_189, fi_190, \
                         fi_191, hi_495, hi_496, hi_497, hi_498, \
                         hi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = -fi_187[k]
                   + f_0 * hi_495[k];

        t_328[k] = -fi_188[k]
                   + f_0 * hi_496[k];

        t_329[k] = -fi_189[k]
                   + f_0 * hi_497[k];

        t_330[k] = -fi_190[k]
                   + f_0 * hi_498[k];

        t_331[k] = -fi_191[k]
                   + f_0 * hi_499[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, fi_192, fi_193, fi_194, fi_195, \
                         fi_196, hi_500, hi_501, hi_502, hi_503, \
                         hi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -fi_192[k]
                   + f_0 * hi_500[k];

        t_333[k] = -fi_193[k]
                   + f_0 * hi_501[k];

        t_334[k] = -fi_194[k]
                   + f_0 * hi_502[k];

        t_335[k] = -fi_195[k]
                   + f_0 * hi_503[k];

        t_336[k] = -2.0 * fi_196[k]
                   + f_0 * hi_504[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, fi_197, fi_198, fi_199, fi_200, \
                         fi_201, hi_505, hi_506, hi_507, hi_508, \
                         hi_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -2.0 * fi_197[k]
                   + f_0 * hi_505[k];

        t_338[k] = -2.0 * fi_198[k]
                   + f_0 * hi_506[k];

        t_339[k] = -2.0 * fi_199[k]
                   + f_0 * hi_507[k];

        t_340[k] = -2.0 * fi_200[k]
                   + f_0 * hi_508[k];

        t_341[k] = -2.0 * fi_201[k]
                   + f_0 * hi_509[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, fi_202, fi_203, fi_204, fi_205, \
                         fi_206, hi_510, hi_511, hi_512, hi_513, \
                         hi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -2.0 * fi_202[k]
                   + f_0 * hi_510[k];

        t_343[k] = -2.0 * fi_203[k]
                   + f_0 * hi_511[k];

        t_344[k] = -2.0 * fi_204[k]
                   + f_0 * hi_512[k];

        t_345[k] = -2.0 * fi_205[k]
                   + f_0 * hi_513[k];

        t_346[k] = -2.0 * fi_206[k]
                   + f_0 * hi_514[k];
    }
}

static auto
compute_prim_geom_10_gi_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fi, const size_t hi,
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

    const auto *fi_207 = buffer.data(fi + 207);
    const auto *fi_208 = buffer.data(fi + 208);
    const auto *fi_209 = buffer.data(fi + 209);
    const auto *fi_210 = buffer.data(fi + 210);
    const auto *fi_211 = buffer.data(fi + 211);
    const auto *fi_212 = buffer.data(fi + 212);
    const auto *fi_213 = buffer.data(fi + 213);
    const auto *fi_214 = buffer.data(fi + 214);
    const auto *fi_215 = buffer.data(fi + 215);
    const auto *fi_216 = buffer.data(fi + 216);
    const auto *fi_217 = buffer.data(fi + 217);
    const auto *fi_218 = buffer.data(fi + 218);
    const auto *fi_219 = buffer.data(fi + 219);
    const auto *fi_220 = buffer.data(fi + 220);
    const auto *fi_221 = buffer.data(fi + 221);
    const auto *fi_222 = buffer.data(fi + 222);
    const auto *fi_223 = buffer.data(fi + 223);
    const auto *fi_224 = buffer.data(fi + 224);
    const auto *fi_225 = buffer.data(fi + 225);
    const auto *fi_226 = buffer.data(fi + 226);
    const auto *fi_227 = buffer.data(fi + 227);
    const auto *fi_228 = buffer.data(fi + 228);
    const auto *fi_229 = buffer.data(fi + 229);
    const auto *fi_230 = buffer.data(fi + 230);
    const auto *fi_231 = buffer.data(fi + 231);
    const auto *fi_232 = buffer.data(fi + 232);
    const auto *fi_233 = buffer.data(fi + 233);
    const auto *fi_234 = buffer.data(fi + 234);
    const auto *fi_235 = buffer.data(fi + 235);
    const auto *fi_236 = buffer.data(fi + 236);
    const auto *fi_237 = buffer.data(fi + 237);
    const auto *fi_238 = buffer.data(fi + 238);
    const auto *fi_239 = buffer.data(fi + 239);
    const auto *fi_240 = buffer.data(fi + 240);
    const auto *fi_241 = buffer.data(fi + 241);
    const auto *fi_242 = buffer.data(fi + 242);
    const auto *fi_243 = buffer.data(fi + 243);
    const auto *fi_244 = buffer.data(fi + 244);
    const auto *fi_245 = buffer.data(fi + 245);
    const auto *fi_246 = buffer.data(fi + 246);
    const auto *fi_247 = buffer.data(fi + 247);
    const auto *fi_248 = buffer.data(fi + 248);
    const auto *fi_249 = buffer.data(fi + 249);
    const auto *fi_250 = buffer.data(fi + 250);
    const auto *fi_251 = buffer.data(fi + 251);
    const auto *fi_252 = buffer.data(fi + 252);
    const auto *fi_253 = buffer.data(fi + 253);
    const auto *fi_254 = buffer.data(fi + 254);
    const auto *fi_255 = buffer.data(fi + 255);
    const auto *fi_256 = buffer.data(fi + 256);
    const auto *fi_257 = buffer.data(fi + 257);
    const auto *fi_258 = buffer.data(fi + 258);
    const auto *fi_259 = buffer.data(fi + 259);
    const auto *fi_260 = buffer.data(fi + 260);
    const auto *fi_261 = buffer.data(fi + 261);
    const auto *fi_262 = buffer.data(fi + 262);
    const auto *fi_263 = buffer.data(fi + 263);
    const auto *fi_264 = buffer.data(fi + 264);
    const auto *fi_265 = buffer.data(fi + 265);
    const auto *fi_266 = buffer.data(fi + 266);
    const auto *fi_267 = buffer.data(fi + 267);
    const auto *fi_268 = buffer.data(fi + 268);
    const auto *fi_269 = buffer.data(fi + 269);
    const auto *fi_270 = buffer.data(fi + 270);
    const auto *fi_271 = buffer.data(fi + 271);
    const auto *fi_272 = buffer.data(fi + 272);
    const auto *fi_273 = buffer.data(fi + 273);
    const auto *fi_274 = buffer.data(fi + 274);
    const auto *fi_275 = buffer.data(fi + 275);
    const auto *fi_276 = buffer.data(fi + 276);
    const auto *fi_277 = buffer.data(fi + 277);
    const auto *fi_278 = buffer.data(fi + 278);
    const auto *fi_279 = buffer.data(fi + 279);

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

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, fi_207, fi_208, fi_209, fi_210, \
                         fi_211, hi_515, hi_516, hi_517, hi_518, \
                         hi_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -2.0 * fi_207[k]
                   + f_0 * hi_515[k];

        t_348[k] = -2.0 * fi_208[k]
                   + f_0 * hi_516[k];

        t_349[k] = -2.0 * fi_209[k]
                   + f_0 * hi_517[k];

        t_350[k] = -2.0 * fi_210[k]
                   + f_0 * hi_518[k];

        t_351[k] = -2.0 * fi_211[k]
                   + f_0 * hi_519[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, fi_212, fi_213, fi_214, fi_215, \
                         fi_216, hi_520, hi_521, hi_522, hi_523, \
                         hi_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -2.0 * fi_212[k]
                   + f_0 * hi_520[k];

        t_353[k] = -2.0 * fi_213[k]
                   + f_0 * hi_521[k];

        t_354[k] = -2.0 * fi_214[k]
                   + f_0 * hi_522[k];

        t_355[k] = -2.0 * fi_215[k]
                   + f_0 * hi_523[k];

        t_356[k] = -2.0 * fi_216[k]
                   + f_0 * hi_524[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, fi_217, fi_218, fi_219, fi_220, \
                         fi_221, hi_525, hi_526, hi_527, hi_528, \
                         hi_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -2.0 * fi_217[k]
                   + f_0 * hi_525[k];

        t_358[k] = -2.0 * fi_218[k]
                   + f_0 * hi_526[k];

        t_359[k] = -2.0 * fi_219[k]
                   + f_0 * hi_527[k];

        t_360[k] = -2.0 * fi_220[k]
                   + f_0 * hi_528[k];

        t_361[k] = -2.0 * fi_221[k]
                   + f_0 * hi_529[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, fi_222, fi_223, fi_224, fi_225, \
                         fi_226, hi_530, hi_531, hi_532, hi_533, \
                         hi_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -2.0 * fi_222[k]
                   + f_0 * hi_530[k];

        t_363[k] = -2.0 * fi_223[k]
                   + f_0 * hi_531[k];

        t_364[k] = -3.0 * fi_224[k]
                   + f_0 * hi_532[k];

        t_365[k] = -3.0 * fi_225[k]
                   + f_0 * hi_533[k];

        t_366[k] = -3.0 * fi_226[k]
                   + f_0 * hi_534[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, fi_227, fi_228, fi_229, fi_230, \
                         fi_231, hi_535, hi_536, hi_537, hi_538, \
                         hi_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -3.0 * fi_227[k]
                   + f_0 * hi_535[k];

        t_368[k] = -3.0 * fi_228[k]
                   + f_0 * hi_536[k];

        t_369[k] = -3.0 * fi_229[k]
                   + f_0 * hi_537[k];

        t_370[k] = -3.0 * fi_230[k]
                   + f_0 * hi_538[k];

        t_371[k] = -3.0 * fi_231[k]
                   + f_0 * hi_539[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, fi_232, fi_233, fi_234, fi_235, \
                         fi_236, hi_540, hi_541, hi_542, hi_543, \
                         hi_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -3.0 * fi_232[k]
                   + f_0 * hi_540[k];

        t_373[k] = -3.0 * fi_233[k]
                   + f_0 * hi_541[k];

        t_374[k] = -3.0 * fi_234[k]
                   + f_0 * hi_542[k];

        t_375[k] = -3.0 * fi_235[k]
                   + f_0 * hi_543[k];

        t_376[k] = -3.0 * fi_236[k]
                   + f_0 * hi_544[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, fi_237, fi_238, fi_239, fi_240, \
                         fi_241, hi_545, hi_546, hi_547, hi_548, \
                         hi_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -3.0 * fi_237[k]
                   + f_0 * hi_545[k];

        t_378[k] = -3.0 * fi_238[k]
                   + f_0 * hi_546[k];

        t_379[k] = -3.0 * fi_239[k]
                   + f_0 * hi_547[k];

        t_380[k] = -3.0 * fi_240[k]
                   + f_0 * hi_548[k];

        t_381[k] = -3.0 * fi_241[k]
                   + f_0 * hi_549[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, fi_242, fi_243, fi_244, fi_245, \
                         fi_246, hi_550, hi_551, hi_552, hi_553, \
                         hi_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -3.0 * fi_242[k]
                   + f_0 * hi_550[k];

        t_383[k] = -3.0 * fi_243[k]
                   + f_0 * hi_551[k];

        t_384[k] = -3.0 * fi_244[k]
                   + f_0 * hi_552[k];

        t_385[k] = -3.0 * fi_245[k]
                   + f_0 * hi_553[k];

        t_386[k] = -3.0 * fi_246[k]
                   + f_0 * hi_554[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, fi_247, fi_248, fi_249, fi_250, \
                         fi_251, hi_555, hi_556, hi_557, hi_558, \
                         hi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -3.0 * fi_247[k]
                   + f_0 * hi_555[k];

        t_388[k] = -3.0 * fi_248[k]
                   + f_0 * hi_556[k];

        t_389[k] = -3.0 * fi_249[k]
                   + f_0 * hi_557[k];

        t_390[k] = -3.0 * fi_250[k]
                   + f_0 * hi_558[k];

        t_391[k] = -3.0 * fi_251[k]
                   + f_0 * hi_559[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, fi_252, fi_253, fi_254, fi_255, \
                         fi_256, hi_560, hi_561, hi_562, hi_563, \
                         hi_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -4.0 * fi_252[k]
                   + f_0 * hi_560[k];

        t_393[k] = -4.0 * fi_253[k]
                   + f_0 * hi_561[k];

        t_394[k] = -4.0 * fi_254[k]
                   + f_0 * hi_562[k];

        t_395[k] = -4.0 * fi_255[k]
                   + f_0 * hi_563[k];

        t_396[k] = -4.0 * fi_256[k]
                   + f_0 * hi_564[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, fi_257, fi_258, fi_259, fi_260, \
                         fi_261, hi_565, hi_566, hi_567, hi_568, \
                         hi_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -4.0 * fi_257[k]
                   + f_0 * hi_565[k];

        t_398[k] = -4.0 * fi_258[k]
                   + f_0 * hi_566[k];

        t_399[k] = -4.0 * fi_259[k]
                   + f_0 * hi_567[k];

        t_400[k] = -4.0 * fi_260[k]
                   + f_0 * hi_568[k];

        t_401[k] = -4.0 * fi_261[k]
                   + f_0 * hi_569[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, fi_262, fi_263, fi_264, fi_265, \
                         fi_266, hi_570, hi_571, hi_572, hi_573, \
                         hi_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -4.0 * fi_262[k]
                   + f_0 * hi_570[k];

        t_403[k] = -4.0 * fi_263[k]
                   + f_0 * hi_571[k];

        t_404[k] = -4.0 * fi_264[k]
                   + f_0 * hi_572[k];

        t_405[k] = -4.0 * fi_265[k]
                   + f_0 * hi_573[k];

        t_406[k] = -4.0 * fi_266[k]
                   + f_0 * hi_574[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, fi_267, fi_268, fi_269, fi_270, \
                         fi_271, hi_575, hi_576, hi_577, hi_578, \
                         hi_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -4.0 * fi_267[k]
                   + f_0 * hi_575[k];

        t_408[k] = -4.0 * fi_268[k]
                   + f_0 * hi_576[k];

        t_409[k] = -4.0 * fi_269[k]
                   + f_0 * hi_577[k];

        t_410[k] = -4.0 * fi_270[k]
                   + f_0 * hi_578[k];

        t_411[k] = -4.0 * fi_271[k]
                   + f_0 * hi_579[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, fi_272, fi_273, fi_274, fi_275, \
                         fi_276, hi_580, hi_581, hi_582, hi_583, \
                         hi_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -4.0 * fi_272[k]
                   + f_0 * hi_580[k];

        t_413[k] = -4.0 * fi_273[k]
                   + f_0 * hi_581[k];

        t_414[k] = -4.0 * fi_274[k]
                   + f_0 * hi_582[k];

        t_415[k] = -4.0 * fi_275[k]
                   + f_0 * hi_583[k];

        t_416[k] = -4.0 * fi_276[k]
                   + f_0 * hi_584[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, fi_277, fi_278, fi_279, hi_585, hi_586, \
                         hi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -4.0 * fi_277[k]
                   + f_0 * hi_585[k];

        t_418[k] = -4.0 * fi_278[k]
                   + f_0 * hi_586[k];

        t_419[k] = -4.0 * fi_279[k]
                   + f_0 * hi_587[k];
    }
}

auto
compute_prim_geom_10_gi_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t fi, const size_t hi,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gi_electron_repulsion_2_piece0(buffer, target, fi, hi, ncols, alpha);

    compute_prim_geom_10_gi_electron_repulsion_2_piece1(buffer, target, fi, hi, ncols, alpha);

    compute_prim_geom_10_gi_electron_repulsion_2_piece2(buffer, target, fi, hi, ncols, alpha);
}

}  // namespace simdt2ceri
