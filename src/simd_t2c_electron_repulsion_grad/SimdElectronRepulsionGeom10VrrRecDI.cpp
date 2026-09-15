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


#include "SimdElectronRepulsionGeom10VrrRecDI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_di_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t pi, const size_t fi,
                                             const size_t ncols, const double alpha) -> void
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

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_21 = buffer.data(pi + 21);
    const auto *pi_22 = buffer.data(pi + 22);
    const auto *pi_23 = buffer.data(pi + 23);
    const auto *pi_24 = buffer.data(pi + 24);
    const auto *pi_25 = buffer.data(pi + 25);
    const auto *pi_26 = buffer.data(pi + 26);
    const auto *pi_27 = buffer.data(pi + 27);
    const auto *pi_28 = buffer.data(pi + 28);
    const auto *pi_29 = buffer.data(pi + 29);
    const auto *pi_30 = buffer.data(pi + 30);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_32 = buffer.data(pi + 32);
    const auto *pi_33 = buffer.data(pi + 33);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_35 = buffer.data(pi + 35);
    const auto *pi_36 = buffer.data(pi + 36);
    const auto *pi_37 = buffer.data(pi + 37);
    const auto *pi_38 = buffer.data(pi + 38);
    const auto *pi_39 = buffer.data(pi + 39);
    const auto *pi_40 = buffer.data(pi + 40);
    const auto *pi_41 = buffer.data(pi + 41);
    const auto *pi_42 = buffer.data(pi + 42);
    const auto *pi_43 = buffer.data(pi + 43);
    const auto *pi_44 = buffer.data(pi + 44);
    const auto *pi_45 = buffer.data(pi + 45);
    const auto *pi_46 = buffer.data(pi + 46);
    const auto *pi_47 = buffer.data(pi + 47);
    const auto *pi_48 = buffer.data(pi + 48);
    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_50 = buffer.data(pi + 50);
    const auto *pi_51 = buffer.data(pi + 51);
    const auto *pi_52 = buffer.data(pi + 52);
    const auto *pi_53 = buffer.data(pi + 53);
    const auto *pi_54 = buffer.data(pi + 54);
    const auto *pi_55 = buffer.data(pi + 55);
    const auto *pi_56 = buffer.data(pi + 56);
    const auto *pi_57 = buffer.data(pi + 57);
    const auto *pi_58 = buffer.data(pi + 58);
    const auto *pi_59 = buffer.data(pi + 59);
    const auto *pi_60 = buffer.data(pi + 60);
    const auto *pi_61 = buffer.data(pi + 61);
    const auto *pi_62 = buffer.data(pi + 62);
    const auto *pi_63 = buffer.data(pi + 63);
    const auto *pi_64 = buffer.data(pi + 64);
    const auto *pi_65 = buffer.data(pi + 65);
    const auto *pi_66 = buffer.data(pi + 66);
    const auto *pi_67 = buffer.data(pi + 67);
    const auto *pi_68 = buffer.data(pi + 68);
    const auto *pi_69 = buffer.data(pi + 69);
    const auto *pi_70 = buffer.data(pi + 70);
    const auto *pi_71 = buffer.data(pi + 71);
    const auto *pi_72 = buffer.data(pi + 72);
    const auto *pi_73 = buffer.data(pi + 73);
    const auto *pi_74 = buffer.data(pi + 74);
    const auto *pi_75 = buffer.data(pi + 75);
    const auto *pi_76 = buffer.data(pi + 76);
    const auto *pi_77 = buffer.data(pi + 77);
    const auto *pi_78 = buffer.data(pi + 78);
    const auto *pi_79 = buffer.data(pi + 79);
    const auto *pi_80 = buffer.data(pi + 80);
    const auto *pi_81 = buffer.data(pi + 81);
    const auto *pi_82 = buffer.data(pi + 82);
    const auto *pi_83 = buffer.data(pi + 83);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pi_0, pi_1, pi_2, pi_3, pi_4, fi_0, fi_1, \
                         fi_2, fi_3, fi_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -2.0 * pi_0[k]
                 + f_0 * fi_0[k];

        t_1[k] = -2.0 * pi_1[k]
                 + f_0 * fi_1[k];

        t_2[k] = -2.0 * pi_2[k]
                 + f_0 * fi_2[k];

        t_3[k] = -2.0 * pi_3[k]
                 + f_0 * fi_3[k];

        t_4[k] = -2.0 * pi_4[k]
                 + f_0 * fi_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pi_5, pi_6, pi_7, pi_8, pi_9, fi_5, fi_6, \
                         fi_7, fi_8, fi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -2.0 * pi_5[k]
                 + f_0 * fi_5[k];

        t_6[k] = -2.0 * pi_6[k]
                 + f_0 * fi_6[k];

        t_7[k] = -2.0 * pi_7[k]
                 + f_0 * fi_7[k];

        t_8[k] = -2.0 * pi_8[k]
                 + f_0 * fi_8[k];

        t_9[k] = -2.0 * pi_9[k]
                 + f_0 * fi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pi_10, pi_11, pi_12, pi_13, pi_14, \
                         fi_10, fi_11, fi_12, fi_13, fi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -2.0 * pi_10[k]
                  + f_0 * fi_10[k];

        t_11[k] = -2.0 * pi_11[k]
                  + f_0 * fi_11[k];

        t_12[k] = -2.0 * pi_12[k]
                  + f_0 * fi_12[k];

        t_13[k] = -2.0 * pi_13[k]
                  + f_0 * fi_13[k];

        t_14[k] = -2.0 * pi_14[k]
                  + f_0 * fi_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pi_15, pi_16, pi_17, pi_18, pi_19, \
                         fi_15, fi_16, fi_17, fi_18, fi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -2.0 * pi_15[k]
                  + f_0 * fi_15[k];

        t_16[k] = -2.0 * pi_16[k]
                  + f_0 * fi_16[k];

        t_17[k] = -2.0 * pi_17[k]
                  + f_0 * fi_17[k];

        t_18[k] = -2.0 * pi_18[k]
                  + f_0 * fi_18[k];

        t_19[k] = -2.0 * pi_19[k]
                  + f_0 * fi_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pi_20, pi_21, pi_22, pi_23, pi_24, \
                         fi_20, fi_21, fi_22, fi_23, fi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -2.0 * pi_20[k]
                  + f_0 * fi_20[k];

        t_21[k] = -2.0 * pi_21[k]
                  + f_0 * fi_21[k];

        t_22[k] = -2.0 * pi_22[k]
                  + f_0 * fi_22[k];

        t_23[k] = -2.0 * pi_23[k]
                  + f_0 * fi_23[k];

        t_24[k] = -2.0 * pi_24[k]
                  + f_0 * fi_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pi_25, pi_26, pi_27, pi_28, pi_29, \
                         fi_25, fi_26, fi_27, fi_28, fi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -2.0 * pi_25[k]
                  + f_0 * fi_25[k];

        t_26[k] = -2.0 * pi_26[k]
                  + f_0 * fi_26[k];

        t_27[k] = -2.0 * pi_27[k]
                  + f_0 * fi_27[k];

        t_28[k] = -pi_28[k]
                  + f_0 * fi_28[k];

        t_29[k] = -pi_29[k]
                  + f_0 * fi_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pi_30, pi_31, pi_32, pi_33, pi_34, \
                         fi_30, fi_31, fi_32, fi_33, fi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -pi_30[k]
                  + f_0 * fi_30[k];

        t_31[k] = -pi_31[k]
                  + f_0 * fi_31[k];

        t_32[k] = -pi_32[k]
                  + f_0 * fi_32[k];

        t_33[k] = -pi_33[k]
                  + f_0 * fi_33[k];

        t_34[k] = -pi_34[k]
                  + f_0 * fi_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pi_35, pi_36, pi_37, pi_38, pi_39, \
                         fi_35, fi_36, fi_37, fi_38, fi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -pi_35[k]
                  + f_0 * fi_35[k];

        t_36[k] = -pi_36[k]
                  + f_0 * fi_36[k];

        t_37[k] = -pi_37[k]
                  + f_0 * fi_37[k];

        t_38[k] = -pi_38[k]
                  + f_0 * fi_38[k];

        t_39[k] = -pi_39[k]
                  + f_0 * fi_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pi_40, pi_41, pi_42, pi_43, pi_44, \
                         fi_40, fi_41, fi_42, fi_43, fi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -pi_40[k]
                  + f_0 * fi_40[k];

        t_41[k] = -pi_41[k]
                  + f_0 * fi_41[k];

        t_42[k] = -pi_42[k]
                  + f_0 * fi_42[k];

        t_43[k] = -pi_43[k]
                  + f_0 * fi_43[k];

        t_44[k] = -pi_44[k]
                  + f_0 * fi_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pi_45, pi_46, pi_47, pi_48, pi_49, \
                         fi_45, fi_46, fi_47, fi_48, fi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -pi_45[k]
                  + f_0 * fi_45[k];

        t_46[k] = -pi_46[k]
                  + f_0 * fi_46[k];

        t_47[k] = -pi_47[k]
                  + f_0 * fi_47[k];

        t_48[k] = -pi_48[k]
                  + f_0 * fi_48[k];

        t_49[k] = -pi_49[k]
                  + f_0 * fi_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pi_50, pi_51, pi_52, pi_53, pi_54, \
                         fi_50, fi_51, fi_52, fi_53, fi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -pi_50[k]
                  + f_0 * fi_50[k];

        t_51[k] = -pi_51[k]
                  + f_0 * fi_51[k];

        t_52[k] = -pi_52[k]
                  + f_0 * fi_52[k];

        t_53[k] = -pi_53[k]
                  + f_0 * fi_53[k];

        t_54[k] = -pi_54[k]
                  + f_0 * fi_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pi_55, pi_56, pi_57, pi_58, pi_59, \
                         fi_55, fi_56, fi_57, fi_58, fi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -pi_55[k]
                  + f_0 * fi_55[k];

        t_56[k] = -pi_56[k]
                  + f_0 * fi_56[k];

        t_57[k] = -pi_57[k]
                  + f_0 * fi_57[k];

        t_58[k] = -pi_58[k]
                  + f_0 * fi_58[k];

        t_59[k] = -pi_59[k]
                  + f_0 * fi_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pi_60, pi_61, pi_62, pi_63, pi_64, \
                         fi_60, fi_61, fi_62, fi_63, fi_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -pi_60[k]
                  + f_0 * fi_60[k];

        t_61[k] = -pi_61[k]
                  + f_0 * fi_61[k];

        t_62[k] = -pi_62[k]
                  + f_0 * fi_62[k];

        t_63[k] = -pi_63[k]
                  + f_0 * fi_63[k];

        t_64[k] = -pi_64[k]
                  + f_0 * fi_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pi_65, pi_66, pi_67, pi_68, pi_69, \
                         fi_65, fi_66, fi_67, fi_68, fi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -pi_65[k]
                  + f_0 * fi_65[k];

        t_66[k] = -pi_66[k]
                  + f_0 * fi_66[k];

        t_67[k] = -pi_67[k]
                  + f_0 * fi_67[k];

        t_68[k] = -pi_68[k]
                  + f_0 * fi_68[k];

        t_69[k] = -pi_69[k]
                  + f_0 * fi_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pi_70, pi_71, pi_72, pi_73, pi_74, \
                         fi_70, fi_71, fi_72, fi_73, fi_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -pi_70[k]
                  + f_0 * fi_70[k];

        t_71[k] = -pi_71[k]
                  + f_0 * fi_71[k];

        t_72[k] = -pi_72[k]
                  + f_0 * fi_72[k];

        t_73[k] = -pi_73[k]
                  + f_0 * fi_73[k];

        t_74[k] = -pi_74[k]
                  + f_0 * fi_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pi_75, pi_76, pi_77, pi_78, pi_79, \
                         fi_75, fi_76, fi_77, fi_78, fi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -pi_75[k]
                  + f_0 * fi_75[k];

        t_76[k] = -pi_76[k]
                  + f_0 * fi_76[k];

        t_77[k] = -pi_77[k]
                  + f_0 * fi_77[k];

        t_78[k] = -pi_78[k]
                  + f_0 * fi_78[k];

        t_79[k] = -pi_79[k]
                  + f_0 * fi_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, pi_80, pi_81, pi_82, pi_83, \
                         fi_80, fi_81, fi_82, fi_83, fi_84, fi_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -pi_80[k]
                  + f_0 * fi_80[k];

        t_81[k] = -pi_81[k]
                  + f_0 * fi_81[k];

        t_82[k] = -pi_82[k]
                  + f_0 * fi_82[k];

        t_83[k] = -pi_83[k]
                  + f_0 * fi_83[k];

        t_84[k] = f_0 * fi_84[k];

        t_85[k] = f_0 * fi_85[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, t_91, t_92, t_93, fi_86, fi_87, fi_88, \
                         fi_89, fi_90, fi_91, fi_92, fi_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_0 * fi_86[k];

        t_87[k] = f_0 * fi_87[k];

        t_88[k] = f_0 * fi_88[k];

        t_89[k] = f_0 * fi_89[k];

        t_90[k] = f_0 * fi_90[k];

        t_91[k] = f_0 * fi_91[k];

        t_92[k] = f_0 * fi_92[k];

        t_93[k] = f_0 * fi_93[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, t_99, t_100, t_101, fi_94, fi_95, \
                         fi_96, fi_97, fi_98, fi_99, fi_100, fi_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_0 * fi_94[k];

        t_95[k] = f_0 * fi_95[k];

        t_96[k] = f_0 * fi_96[k];

        t_97[k] = f_0 * fi_97[k];

        t_98[k] = f_0 * fi_98[k];

        t_99[k] = f_0 * fi_99[k];

        t_100[k] = f_0 * fi_100[k];

        t_101[k] = f_0 * fi_101[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, t_108, t_109, fi_102, \
                         fi_103, fi_104, fi_105, fi_106, fi_107, fi_108, \
                         fi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_0 * fi_102[k];

        t_103[k] = f_0 * fi_103[k];

        t_104[k] = f_0 * fi_104[k];

        t_105[k] = f_0 * fi_105[k];

        t_106[k] = f_0 * fi_106[k];

        t_107[k] = f_0 * fi_107[k];

        t_108[k] = f_0 * fi_108[k];

        t_109[k] = f_0 * fi_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, t_115, t_116, t_117, fi_110, \
                         fi_111, fi_112, fi_113, fi_114, fi_115, fi_116, \
                         fi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_0 * fi_110[k];

        t_111[k] = f_0 * fi_111[k];

        t_112[k] = f_0 * fi_112[k];

        t_113[k] = f_0 * fi_113[k];

        t_114[k] = f_0 * fi_114[k];

        t_115[k] = f_0 * fi_115[k];

        t_116[k] = f_0 * fi_116[k];

        t_117[k] = f_0 * fi_117[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, t_123, t_124, t_125, fi_118, \
                         fi_119, fi_120, fi_121, fi_122, fi_123, fi_124, \
                         fi_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_0 * fi_118[k];

        t_119[k] = f_0 * fi_119[k];

        t_120[k] = f_0 * fi_120[k];

        t_121[k] = f_0 * fi_121[k];

        t_122[k] = f_0 * fi_122[k];

        t_123[k] = f_0 * fi_123[k];

        t_124[k] = f_0 * fi_124[k];

        t_125[k] = f_0 * fi_125[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, t_131, t_132, t_133, fi_126, \
                         fi_127, fi_128, fi_129, fi_130, fi_131, fi_132, \
                         fi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_0 * fi_126[k];

        t_127[k] = f_0 * fi_127[k];

        t_128[k] = f_0 * fi_128[k];

        t_129[k] = f_0 * fi_129[k];

        t_130[k] = f_0 * fi_130[k];

        t_131[k] = f_0 * fi_131[k];

        t_132[k] = f_0 * fi_132[k];

        t_133[k] = f_0 * fi_133[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, t_139, t_140, t_141, fi_134, \
                         fi_135, fi_136, fi_137, fi_138, fi_139, fi_140, \
                         fi_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_0 * fi_134[k];

        t_135[k] = f_0 * fi_135[k];

        t_136[k] = f_0 * fi_136[k];

        t_137[k] = f_0 * fi_137[k];

        t_138[k] = f_0 * fi_138[k];

        t_139[k] = f_0 * fi_139[k];

        t_140[k] = f_0 * fi_140[k];

        t_141[k] = f_0 * fi_141[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, t_147, t_148, t_149, fi_142, \
                         fi_143, fi_144, fi_145, fi_146, fi_147, fi_148, \
                         fi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_0 * fi_142[k];

        t_143[k] = f_0 * fi_143[k];

        t_144[k] = f_0 * fi_144[k];

        t_145[k] = f_0 * fi_145[k];

        t_146[k] = f_0 * fi_146[k];

        t_147[k] = f_0 * fi_147[k];

        t_148[k] = f_0 * fi_148[k];

        t_149[k] = f_0 * fi_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, t_156, t_157, fi_150, \
                         fi_151, fi_152, fi_153, fi_154, fi_155, fi_156, \
                         fi_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_0 * fi_150[k];

        t_151[k] = f_0 * fi_151[k];

        t_152[k] = f_0 * fi_152[k];

        t_153[k] = f_0 * fi_153[k];

        t_154[k] = f_0 * fi_154[k];

        t_155[k] = f_0 * fi_155[k];

        t_156[k] = f_0 * fi_156[k];

        t_157[k] = f_0 * fi_157[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, t_163, t_164, t_165, fi_158, \
                         fi_159, fi_160, fi_161, fi_162, fi_163, fi_164, \
                         fi_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_0 * fi_158[k];

        t_159[k] = f_0 * fi_159[k];

        t_160[k] = f_0 * fi_160[k];

        t_161[k] = f_0 * fi_161[k];

        t_162[k] = f_0 * fi_162[k];

        t_163[k] = f_0 * fi_163[k];

        t_164[k] = f_0 * fi_164[k];

        t_165[k] = f_0 * fi_165[k];
    }

#pragma omp simd aligned(t_166, t_167, fi_166, fi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_0 * fi_166[k];

        t_167[k] = f_0 * fi_167[k];
    }
}

auto
compute_prim_geom_10_di_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t pi, const size_t fi,
                                             const size_t ncols, const double alpha) -> void
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

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_21 = buffer.data(pi + 21);
    const auto *pi_22 = buffer.data(pi + 22);
    const auto *pi_23 = buffer.data(pi + 23);
    const auto *pi_24 = buffer.data(pi + 24);
    const auto *pi_25 = buffer.data(pi + 25);
    const auto *pi_26 = buffer.data(pi + 26);
    const auto *pi_27 = buffer.data(pi + 27);
    const auto *pi_28 = buffer.data(pi + 28);
    const auto *pi_29 = buffer.data(pi + 29);
    const auto *pi_30 = buffer.data(pi + 30);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_32 = buffer.data(pi + 32);
    const auto *pi_33 = buffer.data(pi + 33);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_35 = buffer.data(pi + 35);
    const auto *pi_36 = buffer.data(pi + 36);
    const auto *pi_37 = buffer.data(pi + 37);
    const auto *pi_38 = buffer.data(pi + 38);
    const auto *pi_39 = buffer.data(pi + 39);
    const auto *pi_40 = buffer.data(pi + 40);
    const auto *pi_41 = buffer.data(pi + 41);
    const auto *pi_42 = buffer.data(pi + 42);
    const auto *pi_43 = buffer.data(pi + 43);
    const auto *pi_44 = buffer.data(pi + 44);
    const auto *pi_45 = buffer.data(pi + 45);
    const auto *pi_46 = buffer.data(pi + 46);
    const auto *pi_47 = buffer.data(pi + 47);
    const auto *pi_48 = buffer.data(pi + 48);
    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_50 = buffer.data(pi + 50);
    const auto *pi_51 = buffer.data(pi + 51);
    const auto *pi_52 = buffer.data(pi + 52);
    const auto *pi_53 = buffer.data(pi + 53);
    const auto *pi_54 = buffer.data(pi + 54);
    const auto *pi_55 = buffer.data(pi + 55);
    const auto *pi_56 = buffer.data(pi + 56);
    const auto *pi_57 = buffer.data(pi + 57);
    const auto *pi_58 = buffer.data(pi + 58);
    const auto *pi_59 = buffer.data(pi + 59);
    const auto *pi_60 = buffer.data(pi + 60);
    const auto *pi_61 = buffer.data(pi + 61);
    const auto *pi_62 = buffer.data(pi + 62);
    const auto *pi_63 = buffer.data(pi + 63);
    const auto *pi_64 = buffer.data(pi + 64);
    const auto *pi_65 = buffer.data(pi + 65);
    const auto *pi_66 = buffer.data(pi + 66);
    const auto *pi_67 = buffer.data(pi + 67);
    const auto *pi_68 = buffer.data(pi + 68);
    const auto *pi_69 = buffer.data(pi + 69);
    const auto *pi_70 = buffer.data(pi + 70);
    const auto *pi_71 = buffer.data(pi + 71);
    const auto *pi_72 = buffer.data(pi + 72);
    const auto *pi_73 = buffer.data(pi + 73);
    const auto *pi_74 = buffer.data(pi + 74);
    const auto *pi_75 = buffer.data(pi + 75);
    const auto *pi_76 = buffer.data(pi + 76);
    const auto *pi_77 = buffer.data(pi + 77);
    const auto *pi_78 = buffer.data(pi + 78);
    const auto *pi_79 = buffer.data(pi + 79);
    const auto *pi_80 = buffer.data(pi + 80);
    const auto *pi_81 = buffer.data(pi + 81);
    const auto *pi_82 = buffer.data(pi + 82);
    const auto *pi_83 = buffer.data(pi + 83);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, fi_28, fi_29, fi_30, fi_31, \
                         fi_32, fi_33, fi_34, fi_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_28[k];

        t_1[k] = f_0 * fi_29[k];

        t_2[k] = f_0 * fi_30[k];

        t_3[k] = f_0 * fi_31[k];

        t_4[k] = f_0 * fi_32[k];

        t_5[k] = f_0 * fi_33[k];

        t_6[k] = f_0 * fi_34[k];

        t_7[k] = f_0 * fi_35[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, fi_36, fi_37, fi_38, \
                         fi_39, fi_40, fi_41, fi_42, fi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * fi_36[k];

        t_9[k] = f_0 * fi_37[k];

        t_10[k] = f_0 * fi_38[k];

        t_11[k] = f_0 * fi_39[k];

        t_12[k] = f_0 * fi_40[k];

        t_13[k] = f_0 * fi_41[k];

        t_14[k] = f_0 * fi_42[k];

        t_15[k] = f_0 * fi_43[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, fi_44, fi_45, fi_46, \
                         fi_47, fi_48, fi_49, fi_50, fi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * fi_44[k];

        t_17[k] = f_0 * fi_45[k];

        t_18[k] = f_0 * fi_46[k];

        t_19[k] = f_0 * fi_47[k];

        t_20[k] = f_0 * fi_48[k];

        t_21[k] = f_0 * fi_49[k];

        t_22[k] = f_0 * fi_50[k];

        t_23[k] = f_0 * fi_51[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pi_0, pi_1, fi_52, fi_53, fi_54, \
                         fi_55, fi_84, fi_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * fi_52[k];

        t_25[k] = f_0 * fi_53[k];

        t_26[k] = f_0 * fi_54[k];

        t_27[k] = f_0 * fi_55[k];

        t_28[k] = -pi_0[k]
                  + f_0 * fi_84[k];

        t_29[k] = -pi_1[k]
                  + f_0 * fi_85[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pi_2, pi_3, pi_4, pi_5, pi_6, fi_86, \
                         fi_87, fi_88, fi_89, fi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -pi_2[k]
                  + f_0 * fi_86[k];

        t_31[k] = -pi_3[k]
                  + f_0 * fi_87[k];

        t_32[k] = -pi_4[k]
                  + f_0 * fi_88[k];

        t_33[k] = -pi_5[k]
                  + f_0 * fi_89[k];

        t_34[k] = -pi_6[k]
                  + f_0 * fi_90[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pi_7, pi_8, pi_9, pi_10, pi_11, fi_91, \
                         fi_92, fi_93, fi_94, fi_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -pi_7[k]
                  + f_0 * fi_91[k];

        t_36[k] = -pi_8[k]
                  + f_0 * fi_92[k];

        t_37[k] = -pi_9[k]
                  + f_0 * fi_93[k];

        t_38[k] = -pi_10[k]
                  + f_0 * fi_94[k];

        t_39[k] = -pi_11[k]
                  + f_0 * fi_95[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pi_12, pi_13, pi_14, pi_15, pi_16, \
                         fi_96, fi_97, fi_98, fi_99, fi_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -pi_12[k]
                  + f_0 * fi_96[k];

        t_41[k] = -pi_13[k]
                  + f_0 * fi_97[k];

        t_42[k] = -pi_14[k]
                  + f_0 * fi_98[k];

        t_43[k] = -pi_15[k]
                  + f_0 * fi_99[k];

        t_44[k] = -pi_16[k]
                  + f_0 * fi_100[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pi_17, pi_18, pi_19, pi_20, pi_21, \
                         fi_101, fi_102, fi_103, fi_104, fi_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -pi_17[k]
                  + f_0 * fi_101[k];

        t_46[k] = -pi_18[k]
                  + f_0 * fi_102[k];

        t_47[k] = -pi_19[k]
                  + f_0 * fi_103[k];

        t_48[k] = -pi_20[k]
                  + f_0 * fi_104[k];

        t_49[k] = -pi_21[k]
                  + f_0 * fi_105[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pi_22, pi_23, pi_24, pi_25, pi_26, \
                         fi_106, fi_107, fi_108, fi_109, fi_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -pi_22[k]
                  + f_0 * fi_106[k];

        t_51[k] = -pi_23[k]
                  + f_0 * fi_107[k];

        t_52[k] = -pi_24[k]
                  + f_0 * fi_108[k];

        t_53[k] = -pi_25[k]
                  + f_0 * fi_109[k];

        t_54[k] = -pi_26[k]
                  + f_0 * fi_110[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, pi_27, fi_111, fi_112, \
                         fi_113, fi_114, fi_115, fi_116, fi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -pi_27[k]
                  + f_0 * fi_111[k];

        t_56[k] = f_0 * fi_112[k];

        t_57[k] = f_0 * fi_113[k];

        t_58[k] = f_0 * fi_114[k];

        t_59[k] = f_0 * fi_115[k];

        t_60[k] = f_0 * fi_116[k];

        t_61[k] = f_0 * fi_117[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, t_67, t_68, t_69, fi_118, fi_119, \
                         fi_120, fi_121, fi_122, fi_123, fi_124, \
                         fi_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * fi_118[k];

        t_63[k] = f_0 * fi_119[k];

        t_64[k] = f_0 * fi_120[k];

        t_65[k] = f_0 * fi_121[k];

        t_66[k] = f_0 * fi_122[k];

        t_67[k] = f_0 * fi_123[k];

        t_68[k] = f_0 * fi_124[k];

        t_69[k] = f_0 * fi_125[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, t_75, t_76, t_77, fi_126, fi_127, \
                         fi_128, fi_129, fi_130, fi_131, fi_132, \
                         fi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_0 * fi_126[k];

        t_71[k] = f_0 * fi_127[k];

        t_72[k] = f_0 * fi_128[k];

        t_73[k] = f_0 * fi_129[k];

        t_74[k] = f_0 * fi_130[k];

        t_75[k] = f_0 * fi_131[k];

        t_76[k] = f_0 * fi_132[k];

        t_77[k] = f_0 * fi_133[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, t_84, pi_28, fi_134, fi_135, \
                         fi_136, fi_137, fi_138, fi_139, fi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * fi_134[k];

        t_79[k] = f_0 * fi_135[k];

        t_80[k] = f_0 * fi_136[k];

        t_81[k] = f_0 * fi_137[k];

        t_82[k] = f_0 * fi_138[k];

        t_83[k] = f_0 * fi_139[k];

        t_84[k] = -2.0 * pi_28[k]
                  + f_0 * fi_168[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pi_29, pi_30, pi_31, pi_32, pi_33, \
                         fi_169, fi_170, fi_171, fi_172, fi_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * pi_29[k]
                  + f_0 * fi_169[k];

        t_86[k] = -2.0 * pi_30[k]
                  + f_0 * fi_170[k];

        t_87[k] = -2.0 * pi_31[k]
                  + f_0 * fi_171[k];

        t_88[k] = -2.0 * pi_32[k]
                  + f_0 * fi_172[k];

        t_89[k] = -2.0 * pi_33[k]
                  + f_0 * fi_173[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pi_34, pi_35, pi_36, pi_37, pi_38, \
                         fi_174, fi_175, fi_176, fi_177, fi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * pi_34[k]
                  + f_0 * fi_174[k];

        t_91[k] = -2.0 * pi_35[k]
                  + f_0 * fi_175[k];

        t_92[k] = -2.0 * pi_36[k]
                  + f_0 * fi_176[k];

        t_93[k] = -2.0 * pi_37[k]
                  + f_0 * fi_177[k];

        t_94[k] = -2.0 * pi_38[k]
                  + f_0 * fi_178[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pi_39, pi_40, pi_41, pi_42, pi_43, \
                         fi_179, fi_180, fi_181, fi_182, fi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * pi_39[k]
                  + f_0 * fi_179[k];

        t_96[k] = -2.0 * pi_40[k]
                  + f_0 * fi_180[k];

        t_97[k] = -2.0 * pi_41[k]
                  + f_0 * fi_181[k];

        t_98[k] = -2.0 * pi_42[k]
                  + f_0 * fi_182[k];

        t_99[k] = -2.0 * pi_43[k]
                  + f_0 * fi_183[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pi_44, pi_45, pi_46, pi_47, pi_48, \
                         fi_184, fi_185, fi_186, fi_187, fi_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * pi_44[k]
                   + f_0 * fi_184[k];

        t_101[k] = -2.0 * pi_45[k]
                   + f_0 * fi_185[k];

        t_102[k] = -2.0 * pi_46[k]
                   + f_0 * fi_186[k];

        t_103[k] = -2.0 * pi_47[k]
                   + f_0 * fi_187[k];

        t_104[k] = -2.0 * pi_48[k]
                   + f_0 * fi_188[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pi_49, pi_50, pi_51, pi_52, pi_53, \
                         fi_189, fi_190, fi_191, fi_192, fi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * pi_49[k]
                   + f_0 * fi_189[k];

        t_106[k] = -2.0 * pi_50[k]
                   + f_0 * fi_190[k];

        t_107[k] = -2.0 * pi_51[k]
                   + f_0 * fi_191[k];

        t_108[k] = -2.0 * pi_52[k]
                   + f_0 * fi_192[k];

        t_109[k] = -2.0 * pi_53[k]
                   + f_0 * fi_193[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pi_54, pi_55, pi_56, pi_57, pi_58, \
                         fi_194, fi_195, fi_196, fi_197, fi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * pi_54[k]
                   + f_0 * fi_194[k];

        t_111[k] = -2.0 * pi_55[k]
                   + f_0 * fi_195[k];

        t_112[k] = -pi_56[k]
                   + f_0 * fi_196[k];

        t_113[k] = -pi_57[k]
                   + f_0 * fi_197[k];

        t_114[k] = -pi_58[k]
                   + f_0 * fi_198[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pi_59, pi_60, pi_61, pi_62, pi_63, \
                         fi_199, fi_200, fi_201, fi_202, fi_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -pi_59[k]
                   + f_0 * fi_199[k];

        t_116[k] = -pi_60[k]
                   + f_0 * fi_200[k];

        t_117[k] = -pi_61[k]
                   + f_0 * fi_201[k];

        t_118[k] = -pi_62[k]
                   + f_0 * fi_202[k];

        t_119[k] = -pi_63[k]
                   + f_0 * fi_203[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pi_64, pi_65, pi_66, pi_67, pi_68, \
                         fi_204, fi_205, fi_206, fi_207, fi_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -pi_64[k]
                   + f_0 * fi_204[k];

        t_121[k] = -pi_65[k]
                   + f_0 * fi_205[k];

        t_122[k] = -pi_66[k]
                   + f_0 * fi_206[k];

        t_123[k] = -pi_67[k]
                   + f_0 * fi_207[k];

        t_124[k] = -pi_68[k]
                   + f_0 * fi_208[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pi_69, pi_70, pi_71, pi_72, pi_73, \
                         fi_209, fi_210, fi_211, fi_212, fi_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -pi_69[k]
                   + f_0 * fi_209[k];

        t_126[k] = -pi_70[k]
                   + f_0 * fi_210[k];

        t_127[k] = -pi_71[k]
                   + f_0 * fi_211[k];

        t_128[k] = -pi_72[k]
                   + f_0 * fi_212[k];

        t_129[k] = -pi_73[k]
                   + f_0 * fi_213[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, pi_74, pi_75, pi_76, pi_77, pi_78, \
                         fi_214, fi_215, fi_216, fi_217, fi_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -pi_74[k]
                   + f_0 * fi_214[k];

        t_131[k] = -pi_75[k]
                   + f_0 * fi_215[k];

        t_132[k] = -pi_76[k]
                   + f_0 * fi_216[k];

        t_133[k] = -pi_77[k]
                   + f_0 * fi_217[k];

        t_134[k] = -pi_78[k]
                   + f_0 * fi_218[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pi_79, pi_80, pi_81, pi_82, pi_83, \
                         fi_219, fi_220, fi_221, fi_222, fi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -pi_79[k]
                   + f_0 * fi_219[k];

        t_136[k] = -pi_80[k]
                   + f_0 * fi_220[k];

        t_137[k] = -pi_81[k]
                   + f_0 * fi_221[k];

        t_138[k] = -pi_82[k]
                   + f_0 * fi_222[k];

        t_139[k] = -pi_83[k]
                   + f_0 * fi_223[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, t_146, t_147, fi_224, \
                         fi_225, fi_226, fi_227, fi_228, fi_229, fi_230, \
                         fi_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * fi_224[k];

        t_141[k] = f_0 * fi_225[k];

        t_142[k] = f_0 * fi_226[k];

        t_143[k] = f_0 * fi_227[k];

        t_144[k] = f_0 * fi_228[k];

        t_145[k] = f_0 * fi_229[k];

        t_146[k] = f_0 * fi_230[k];

        t_147[k] = f_0 * fi_231[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, t_153, t_154, t_155, fi_232, \
                         fi_233, fi_234, fi_235, fi_236, fi_237, fi_238, \
                         fi_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * fi_232[k];

        t_149[k] = f_0 * fi_233[k];

        t_150[k] = f_0 * fi_234[k];

        t_151[k] = f_0 * fi_235[k];

        t_152[k] = f_0 * fi_236[k];

        t_153[k] = f_0 * fi_237[k];

        t_154[k] = f_0 * fi_238[k];

        t_155[k] = f_0 * fi_239[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, t_162, t_163, fi_240, \
                         fi_241, fi_242, fi_243, fi_244, fi_245, fi_246, \
                         fi_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * fi_240[k];

        t_157[k] = f_0 * fi_241[k];

        t_158[k] = f_0 * fi_242[k];

        t_159[k] = f_0 * fi_243[k];

        t_160[k] = f_0 * fi_244[k];

        t_161[k] = f_0 * fi_245[k];

        t_162[k] = f_0 * fi_246[k];

        t_163[k] = f_0 * fi_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, fi_248, fi_249, fi_250, \
                         fi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_0 * fi_248[k];

        t_165[k] = f_0 * fi_249[k];

        t_166[k] = f_0 * fi_250[k];

        t_167[k] = f_0 * fi_251[k];
    }
}

auto
compute_prim_geom_10_di_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t pi, const size_t fi,
                                             const size_t ncols, const double alpha) -> void
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

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_21 = buffer.data(pi + 21);
    const auto *pi_22 = buffer.data(pi + 22);
    const auto *pi_23 = buffer.data(pi + 23);
    const auto *pi_24 = buffer.data(pi + 24);
    const auto *pi_25 = buffer.data(pi + 25);
    const auto *pi_26 = buffer.data(pi + 26);
    const auto *pi_27 = buffer.data(pi + 27);
    const auto *pi_28 = buffer.data(pi + 28);
    const auto *pi_29 = buffer.data(pi + 29);
    const auto *pi_30 = buffer.data(pi + 30);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_32 = buffer.data(pi + 32);
    const auto *pi_33 = buffer.data(pi + 33);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_35 = buffer.data(pi + 35);
    const auto *pi_36 = buffer.data(pi + 36);
    const auto *pi_37 = buffer.data(pi + 37);
    const auto *pi_38 = buffer.data(pi + 38);
    const auto *pi_39 = buffer.data(pi + 39);
    const auto *pi_40 = buffer.data(pi + 40);
    const auto *pi_41 = buffer.data(pi + 41);
    const auto *pi_42 = buffer.data(pi + 42);
    const auto *pi_43 = buffer.data(pi + 43);
    const auto *pi_44 = buffer.data(pi + 44);
    const auto *pi_45 = buffer.data(pi + 45);
    const auto *pi_46 = buffer.data(pi + 46);
    const auto *pi_47 = buffer.data(pi + 47);
    const auto *pi_48 = buffer.data(pi + 48);
    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_50 = buffer.data(pi + 50);
    const auto *pi_51 = buffer.data(pi + 51);
    const auto *pi_52 = buffer.data(pi + 52);
    const auto *pi_53 = buffer.data(pi + 53);
    const auto *pi_54 = buffer.data(pi + 54);
    const auto *pi_55 = buffer.data(pi + 55);
    const auto *pi_56 = buffer.data(pi + 56);
    const auto *pi_57 = buffer.data(pi + 57);
    const auto *pi_58 = buffer.data(pi + 58);
    const auto *pi_59 = buffer.data(pi + 59);
    const auto *pi_60 = buffer.data(pi + 60);
    const auto *pi_61 = buffer.data(pi + 61);
    const auto *pi_62 = buffer.data(pi + 62);
    const auto *pi_63 = buffer.data(pi + 63);
    const auto *pi_64 = buffer.data(pi + 64);
    const auto *pi_65 = buffer.data(pi + 65);
    const auto *pi_66 = buffer.data(pi + 66);
    const auto *pi_67 = buffer.data(pi + 67);
    const auto *pi_68 = buffer.data(pi + 68);
    const auto *pi_69 = buffer.data(pi + 69);
    const auto *pi_70 = buffer.data(pi + 70);
    const auto *pi_71 = buffer.data(pi + 71);
    const auto *pi_72 = buffer.data(pi + 72);
    const auto *pi_73 = buffer.data(pi + 73);
    const auto *pi_74 = buffer.data(pi + 74);
    const auto *pi_75 = buffer.data(pi + 75);
    const auto *pi_76 = buffer.data(pi + 76);
    const auto *pi_77 = buffer.data(pi + 77);
    const auto *pi_78 = buffer.data(pi + 78);
    const auto *pi_79 = buffer.data(pi + 79);
    const auto *pi_80 = buffer.data(pi + 80);
    const auto *pi_81 = buffer.data(pi + 81);
    const auto *pi_82 = buffer.data(pi + 82);
    const auto *pi_83 = buffer.data(pi + 83);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, fi_56, fi_57, fi_58, fi_59, \
                         fi_60, fi_61, fi_62, fi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_56[k];

        t_1[k] = f_0 * fi_57[k];

        t_2[k] = f_0 * fi_58[k];

        t_3[k] = f_0 * fi_59[k];

        t_4[k] = f_0 * fi_60[k];

        t_5[k] = f_0 * fi_61[k];

        t_6[k] = f_0 * fi_62[k];

        t_7[k] = f_0 * fi_63[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, fi_64, fi_65, fi_66, \
                         fi_67, fi_68, fi_69, fi_70, fi_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * fi_64[k];

        t_9[k] = f_0 * fi_65[k];

        t_10[k] = f_0 * fi_66[k];

        t_11[k] = f_0 * fi_67[k];

        t_12[k] = f_0 * fi_68[k];

        t_13[k] = f_0 * fi_69[k];

        t_14[k] = f_0 * fi_70[k];

        t_15[k] = f_0 * fi_71[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, fi_72, fi_73, fi_74, \
                         fi_75, fi_76, fi_77, fi_78, fi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * fi_72[k];

        t_17[k] = f_0 * fi_73[k];

        t_18[k] = f_0 * fi_74[k];

        t_19[k] = f_0 * fi_75[k];

        t_20[k] = f_0 * fi_76[k];

        t_21[k] = f_0 * fi_77[k];

        t_22[k] = f_0 * fi_78[k];

        t_23[k] = f_0 * fi_79[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, fi_80, fi_81, fi_82, \
                         fi_83, fi_112, fi_113, fi_114, fi_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * fi_80[k];

        t_25[k] = f_0 * fi_81[k];

        t_26[k] = f_0 * fi_82[k];

        t_27[k] = f_0 * fi_83[k];

        t_28[k] = f_0 * fi_112[k];

        t_29[k] = f_0 * fi_113[k];

        t_30[k] = f_0 * fi_114[k];

        t_31[k] = f_0 * fi_115[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, fi_116, fi_117, \
                         fi_118, fi_119, fi_120, fi_121, fi_122, \
                         fi_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * fi_116[k];

        t_33[k] = f_0 * fi_117[k];

        t_34[k] = f_0 * fi_118[k];

        t_35[k] = f_0 * fi_119[k];

        t_36[k] = f_0 * fi_120[k];

        t_37[k] = f_0 * fi_121[k];

        t_38[k] = f_0 * fi_122[k];

        t_39[k] = f_0 * fi_123[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, fi_124, fi_125, \
                         fi_126, fi_127, fi_128, fi_129, fi_130, \
                         fi_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * fi_124[k];

        t_41[k] = f_0 * fi_125[k];

        t_42[k] = f_0 * fi_126[k];

        t_43[k] = f_0 * fi_127[k];

        t_44[k] = f_0 * fi_128[k];

        t_45[k] = f_0 * fi_129[k];

        t_46[k] = f_0 * fi_130[k];

        t_47[k] = f_0 * fi_131[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, fi_132, fi_133, \
                         fi_134, fi_135, fi_136, fi_137, fi_138, \
                         fi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * fi_132[k];

        t_49[k] = f_0 * fi_133[k];

        t_50[k] = f_0 * fi_134[k];

        t_51[k] = f_0 * fi_135[k];

        t_52[k] = f_0 * fi_136[k];

        t_53[k] = f_0 * fi_137[k];

        t_54[k] = f_0 * fi_138[k];

        t_55[k] = f_0 * fi_139[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pi_0, pi_1, pi_2, pi_3, pi_4, fi_140, \
                         fi_141, fi_142, fi_143, fi_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -pi_0[k]
                  + f_0 * fi_140[k];

        t_57[k] = -pi_1[k]
                  + f_0 * fi_141[k];

        t_58[k] = -pi_2[k]
                  + f_0 * fi_142[k];

        t_59[k] = -pi_3[k]
                  + f_0 * fi_143[k];

        t_60[k] = -pi_4[k]
                  + f_0 * fi_144[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pi_5, pi_6, pi_7, pi_8, pi_9, fi_145, \
                         fi_146, fi_147, fi_148, fi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -pi_5[k]
                  + f_0 * fi_145[k];

        t_62[k] = -pi_6[k]
                  + f_0 * fi_146[k];

        t_63[k] = -pi_7[k]
                  + f_0 * fi_147[k];

        t_64[k] = -pi_8[k]
                  + f_0 * fi_148[k];

        t_65[k] = -pi_9[k]
                  + f_0 * fi_149[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pi_10, pi_11, pi_12, pi_13, pi_14, \
                         fi_150, fi_151, fi_152, fi_153, fi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -pi_10[k]
                  + f_0 * fi_150[k];

        t_67[k] = -pi_11[k]
                  + f_0 * fi_151[k];

        t_68[k] = -pi_12[k]
                  + f_0 * fi_152[k];

        t_69[k] = -pi_13[k]
                  + f_0 * fi_153[k];

        t_70[k] = -pi_14[k]
                  + f_0 * fi_154[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, pi_15, pi_16, pi_17, pi_18, pi_19, \
                         fi_155, fi_156, fi_157, fi_158, fi_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -pi_15[k]
                  + f_0 * fi_155[k];

        t_72[k] = -pi_16[k]
                  + f_0 * fi_156[k];

        t_73[k] = -pi_17[k]
                  + f_0 * fi_157[k];

        t_74[k] = -pi_18[k]
                  + f_0 * fi_158[k];

        t_75[k] = -pi_19[k]
                  + f_0 * fi_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pi_20, pi_21, pi_22, pi_23, pi_24, \
                         fi_160, fi_161, fi_162, fi_163, fi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -pi_20[k]
                  + f_0 * fi_160[k];

        t_77[k] = -pi_21[k]
                  + f_0 * fi_161[k];

        t_78[k] = -pi_22[k]
                  + f_0 * fi_162[k];

        t_79[k] = -pi_23[k]
                  + f_0 * fi_163[k];

        t_80[k] = -pi_24[k]
                  + f_0 * fi_164[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, pi_25, pi_26, pi_27, fi_165, \
                         fi_166, fi_167, fi_196, fi_197, fi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -pi_25[k]
                  + f_0 * fi_165[k];

        t_82[k] = -pi_26[k]
                  + f_0 * fi_166[k];

        t_83[k] = -pi_27[k]
                  + f_0 * fi_167[k];

        t_84[k] = f_0 * fi_196[k];

        t_85[k] = f_0 * fi_197[k];

        t_86[k] = f_0 * fi_198[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, t_93, t_94, fi_199, fi_200, \
                         fi_201, fi_202, fi_203, fi_204, fi_205, \
                         fi_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_0 * fi_199[k];

        t_88[k] = f_0 * fi_200[k];

        t_89[k] = f_0 * fi_201[k];

        t_90[k] = f_0 * fi_202[k];

        t_91[k] = f_0 * fi_203[k];

        t_92[k] = f_0 * fi_204[k];

        t_93[k] = f_0 * fi_205[k];

        t_94[k] = f_0 * fi_206[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, t_100, t_101, t_102, fi_207, fi_208, \
                         fi_209, fi_210, fi_211, fi_212, fi_213, \
                         fi_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_0 * fi_207[k];

        t_96[k] = f_0 * fi_208[k];

        t_97[k] = f_0 * fi_209[k];

        t_98[k] = f_0 * fi_210[k];

        t_99[k] = f_0 * fi_211[k];

        t_100[k] = f_0 * fi_212[k];

        t_101[k] = f_0 * fi_213[k];

        t_102[k] = f_0 * fi_214[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, t_109, t_110, fi_215, \
                         fi_216, fi_217, fi_218, fi_219, fi_220, fi_221, \
                         fi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_0 * fi_215[k];

        t_104[k] = f_0 * fi_216[k];

        t_105[k] = f_0 * fi_217[k];

        t_106[k] = f_0 * fi_218[k];

        t_107[k] = f_0 * fi_219[k];

        t_108[k] = f_0 * fi_220[k];

        t_109[k] = f_0 * fi_221[k];

        t_110[k] = f_0 * fi_222[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, pi_28, pi_29, pi_30, pi_31, \
                         fi_223, fi_224, fi_225, fi_226, fi_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_0 * fi_223[k];

        t_112[k] = -pi_28[k]
                   + f_0 * fi_224[k];

        t_113[k] = -pi_29[k]
                   + f_0 * fi_225[k];

        t_114[k] = -pi_30[k]
                   + f_0 * fi_226[k];

        t_115[k] = -pi_31[k]
                   + f_0 * fi_227[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, pi_32, pi_33, pi_34, pi_35, pi_36, \
                         fi_228, fi_229, fi_230, fi_231, fi_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = -pi_32[k]
                   + f_0 * fi_228[k];

        t_117[k] = -pi_33[k]
                   + f_0 * fi_229[k];

        t_118[k] = -pi_34[k]
                   + f_0 * fi_230[k];

        t_119[k] = -pi_35[k]
                   + f_0 * fi_231[k];

        t_120[k] = -pi_36[k]
                   + f_0 * fi_232[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pi_37, pi_38, pi_39, pi_40, pi_41, \
                         fi_233, fi_234, fi_235, fi_236, fi_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -pi_37[k]
                   + f_0 * fi_233[k];

        t_122[k] = -pi_38[k]
                   + f_0 * fi_234[k];

        t_123[k] = -pi_39[k]
                   + f_0 * fi_235[k];

        t_124[k] = -pi_40[k]
                   + f_0 * fi_236[k];

        t_125[k] = -pi_41[k]
                   + f_0 * fi_237[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, pi_42, pi_43, pi_44, pi_45, pi_46, \
                         fi_238, fi_239, fi_240, fi_241, fi_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = -pi_42[k]
                   + f_0 * fi_238[k];

        t_127[k] = -pi_43[k]
                   + f_0 * fi_239[k];

        t_128[k] = -pi_44[k]
                   + f_0 * fi_240[k];

        t_129[k] = -pi_45[k]
                   + f_0 * fi_241[k];

        t_130[k] = -pi_46[k]
                   + f_0 * fi_242[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, pi_47, pi_48, pi_49, pi_50, pi_51, \
                         fi_243, fi_244, fi_245, fi_246, fi_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = -pi_47[k]
                   + f_0 * fi_243[k];

        t_132[k] = -pi_48[k]
                   + f_0 * fi_244[k];

        t_133[k] = -pi_49[k]
                   + f_0 * fi_245[k];

        t_134[k] = -pi_50[k]
                   + f_0 * fi_246[k];

        t_135[k] = -pi_51[k]
                   + f_0 * fi_247[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pi_52, pi_53, pi_54, pi_55, pi_56, \
                         fi_248, fi_249, fi_250, fi_251, fi_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = -pi_52[k]
                   + f_0 * fi_248[k];

        t_137[k] = -pi_53[k]
                   + f_0 * fi_249[k];

        t_138[k] = -pi_54[k]
                   + f_0 * fi_250[k];

        t_139[k] = -pi_55[k]
                   + f_0 * fi_251[k];

        t_140[k] = -2.0 * pi_56[k]
                   + f_0 * fi_252[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, pi_57, pi_58, pi_59, pi_60, pi_61, \
                         fi_253, fi_254, fi_255, fi_256, fi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = -2.0 * pi_57[k]
                   + f_0 * fi_253[k];

        t_142[k] = -2.0 * pi_58[k]
                   + f_0 * fi_254[k];

        t_143[k] = -2.0 * pi_59[k]
                   + f_0 * fi_255[k];

        t_144[k] = -2.0 * pi_60[k]
                   + f_0 * fi_256[k];

        t_145[k] = -2.0 * pi_61[k]
                   + f_0 * fi_257[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, pi_62, pi_63, pi_64, pi_65, pi_66, \
                         fi_258, fi_259, fi_260, fi_261, fi_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = -2.0 * pi_62[k]
                   + f_0 * fi_258[k];

        t_147[k] = -2.0 * pi_63[k]
                   + f_0 * fi_259[k];

        t_148[k] = -2.0 * pi_64[k]
                   + f_0 * fi_260[k];

        t_149[k] = -2.0 * pi_65[k]
                   + f_0 * fi_261[k];

        t_150[k] = -2.0 * pi_66[k]
                   + f_0 * fi_262[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, pi_67, pi_68, pi_69, pi_70, pi_71, \
                         fi_263, fi_264, fi_265, fi_266, fi_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = -2.0 * pi_67[k]
                   + f_0 * fi_263[k];

        t_152[k] = -2.0 * pi_68[k]
                   + f_0 * fi_264[k];

        t_153[k] = -2.0 * pi_69[k]
                   + f_0 * fi_265[k];

        t_154[k] = -2.0 * pi_70[k]
                   + f_0 * fi_266[k];

        t_155[k] = -2.0 * pi_71[k]
                   + f_0 * fi_267[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pi_72, pi_73, pi_74, pi_75, pi_76, \
                         fi_268, fi_269, fi_270, fi_271, fi_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = -2.0 * pi_72[k]
                   + f_0 * fi_268[k];

        t_157[k] = -2.0 * pi_73[k]
                   + f_0 * fi_269[k];

        t_158[k] = -2.0 * pi_74[k]
                   + f_0 * fi_270[k];

        t_159[k] = -2.0 * pi_75[k]
                   + f_0 * fi_271[k];

        t_160[k] = -2.0 * pi_76[k]
                   + f_0 * fi_272[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, pi_77, pi_78, pi_79, pi_80, pi_81, \
                         fi_273, fi_274, fi_275, fi_276, fi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = -2.0 * pi_77[k]
                   + f_0 * fi_273[k];

        t_162[k] = -2.0 * pi_78[k]
                   + f_0 * fi_274[k];

        t_163[k] = -2.0 * pi_79[k]
                   + f_0 * fi_275[k];

        t_164[k] = -2.0 * pi_80[k]
                   + f_0 * fi_276[k];

        t_165[k] = -2.0 * pi_81[k]
                   + f_0 * fi_277[k];
    }

#pragma omp simd aligned(t_166, t_167, pi_82, pi_83, fi_278, fi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = -2.0 * pi_82[k]
                   + f_0 * fi_278[k];

        t_167[k] = -2.0 * pi_83[k]
                   + f_0 * fi_279[k];
    }
}

}  // namespace simdt2ceri
